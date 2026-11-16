"""Measure and cache selected-allele cross-locus sequence similarity.

This module compares one deterministically selected nucleotide and translated
amino-acid allele between different loci in an extracted MLST or cgMLST
schema. It is an intrinsic schema analysis: no genome assemblies are used, and
similarity must not be interpreted as proof of copy number or evolutionary
paralogy.

The input schemes do not provide curated representative sequences. Selection
therefore uses allele length as a reproducible approximation rather than
claiming that the chosen sequence is a biological or sequence-distance medoid.
The modal length is calculated from every non-empty allele in the locus; ties
between equally frequent lengths select the shortest length, matching the
repository's existing mode convention. The valid coding allele closest to
that modal length is selected, with original FASTA order breaking equal-
distance ties. A locus without any valid coding allele can still receive a
nucleotide-only selection: the eligible IUPAC DNA allele closest to the same
mode, again with FASTA order breaking ties.

NCBI BLAST+ performs independent nucleotide and protein searches. Identity is
the identical aligned character count divided by alignment length. Query and
subject coverage use their respective aligned coordinate spans. BLAST Score
Ratio (BSR) is directional: the best cross-locus raw score for a query is
divided by the selected query allele's self-alignment raw score. The highest-
raw-score HSP is retained for each directed locus pair. Separate local HSPs
are not merged.

Only BLAST-reported local similarities are returned. A missing row therefore
means that no hit passed the configured BLAST search floor; it is not reported
as zero similarity. The default E-value floor of 0.001 follows the search
setting used by Schema Refinery. It is a search sensitivity setting, not a
literature-derived pass/fail or candidate-paralogy threshold.

Results are cached atomically inside the extracted schema using the FASTA
source signature, BLAST versions, allele-selection method, and search
settings. Temporary FASTA files, BLAST databases, and raw output are removed
after calculation. The public API is independent of Marimo and includes a
typed progress callback for notebook and terminal monitoring.
"""

import argparse
import gzip
import json
import os
import subprocess
import sys
import tempfile
from collections import Counter
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import TextIO

from .alignment import IUPAC_DNA_BASES
from .allele_info import _evaluate_allele_sequence, _translate_valid_allele
from .info import FASTA_FILENAME_SUFFIXES

SCHEMA_SIMILARITY_CACHE_FILENAME = ".schema-similarity-cache.json"
SCHEMA_SIMILARITY_CACHE_VERSION = 2
SIMILARITY_METHOD_IDENTIFIER = "blast-modal-length-selected-alleles-bsr-v2"
DEFAULT_BLAST_EVALUE = 0.001
SEQUENCE_LEVEL_NUCLEOTIDE = "nucleotide"
SEQUENCE_LEVEL_AMINO_ACID = "amino_acid"
SIMILARITY_PROGRESS_STAGE_COUNT = 5

__all__ = [
    "DEFAULT_BLAST_EVALUE",
    "SCHEMA_SIMILARITY_CACHE_FILENAME",
    "SCHEMA_SIMILARITY_CACHE_VERSION",
    "SEQUENCE_LEVEL_AMINO_ACID",
    "SEQUENCE_LEVEL_NUCLEOTIDE",
    "SIMILARITY_METHOD_IDENTIFIER",
    "CrossLocusSimilarityMatch",
    "SchemaSimilarityAnalysis",
    "SimilarityProgress",
    "get_schema_similarity_analysis",
]


@dataclass(frozen=True)
class CrossLocusSimilarityMatch:
    """Best local similarity for one directed locus pair and sequence level.

    Attributes:
        sequence_level: ``nucleotide`` or ``amino_acid``.
        query_locus_name: Locus supplying the selected query allele.
        subject_locus_name: Different locus supplying the database sequence.
        query_allele_identifier: Original query FASTA identifier.
        subject_allele_identifier: Original subject FASTA identifier.
        identity_percentage: Identical aligned characters divided by HSP
            alignment length, as a percentage.
        alignment_length: Length of the selected local HSP.
        query_coverage_percentage: Covered query coordinate span divided by
            the complete query length.
        subject_coverage_percentage: Covered subject coordinate span divided
            by the complete subject length.
        minimum_coverage_percentage: Smaller of query and subject coverage.
        raw_score: BLAST raw score for the selected HSP.
        blast_score_ratio: Raw score divided by the query self-alignment raw
            score. Values are preserved rather than clamped to one.
        nucleotide_subject_strand: Subject strand for nucleotide matches;
            ``None`` for amino-acid matches.
    """

    sequence_level: str
    query_locus_name: str
    subject_locus_name: str
    query_allele_identifier: str
    subject_allele_identifier: str
    identity_percentage: float
    alignment_length: int
    query_coverage_percentage: float
    subject_coverage_percentage: float
    minimum_coverage_percentage: float
    raw_score: float
    blast_score_ratio: float
    nucleotide_subject_strand: str | None


@dataclass(frozen=True)
class SchemaSimilarityAnalysis:
    """Complete selected-allele cross-locus similarity analysis."""

    locus_count: int
    nucleotide_selected_locus_count: int
    amino_acid_selected_locus_count: int
    nucleotide_possible_directed_pair_count: int
    amino_acid_possible_directed_pair_count: int
    matches: tuple[CrossLocusSimilarityMatch, ...]
    blastn_version: str
    blastp_version: str
    method_identifier: str
    blast_evalue: float
    reused_from_cache: bool

    @property
    def nucleotide_match_count(self) -> int:
        """Return the number of directed nucleotide locus-pair hits."""
        return sum(
            match.sequence_level == SEQUENCE_LEVEL_NUCLEOTIDE
            for match in self.matches
        )

    @property
    def amino_acid_match_count(self) -> int:
        """Return the number of directed amino-acid locus-pair hits."""
        return sum(
            match.sequence_level == SEQUENCE_LEVEL_AMINO_ACID
            for match in self.matches
        )


@dataclass(frozen=True)
class SimilarityProgress:
    """Progress update emitted by :func:`get_schema_similarity_analysis`."""

    completed_stage_count: int
    total_stage_count: int
    stage_name: str
    message: str
    reused_from_cache: bool = False


@dataclass(frozen=True)
class _FastaRecord:
    """One normalized FASTA record."""

    identifier: str
    sequence: str


@dataclass(frozen=True)
class _SelectedAlleleSequence:
    """One selected allele encoded with a BLAST-safe internal identifier."""

    internal_identifier: str
    locus_name: str
    original_identifier: str
    sequence: str


@dataclass(frozen=True)
class _BlastHsp:
    """Parsed BLAST HSP fields used to derive similarity measurements."""

    query_identifier: str
    subject_identifier: str
    query_length: int
    subject_length: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    alignment_length: int
    identical_character_count: int
    raw_score: float
    subject_strand: str | None


@dataclass(frozen=True)
class _SimilarityLevelResult:
    """Matches returned by one nucleotide or protein BLAST search."""

    sequence_level: str
    matches: tuple[CrossLocusSimilarityMatch, ...]


type _SourceSignature = tuple[tuple[str, int, int], ...]
type SimilarityProgressCallback = Callable[[SimilarityProgress], None]


def _resolve_schema_directory_path(schema_directory_path: str | Path) -> Path:
    """Resolve and validate an extracted schema directory."""
    resolved_schema_directory_path = Path(schema_directory_path).resolve()
    if not resolved_schema_directory_path.exists():
        raise FileNotFoundError(
            f"Schema directory does not exist: {resolved_schema_directory_path}"
        )
    if not resolved_schema_directory_path.is_dir():
        raise NotADirectoryError(
            f"Schema path is not a directory: {resolved_schema_directory_path}"
        )
    return resolved_schema_directory_path


def _is_fasta_file(candidate_file_path: Path) -> bool:
    """Return whether a path is a recognized FASTA file."""
    return candidate_file_path.is_file() and candidate_file_path.name.lower().endswith(
        FASTA_FILENAME_SUFFIXES
    )


def _find_all_fasta_file_paths(schema_directory_path: Path) -> tuple[Path, ...]:
    """Return every recognized locus FASTA in deterministic path order."""
    return tuple(
        sorted(
            candidate_file_path
            for candidate_file_path in schema_directory_path.rglob("*")
            if _is_fasta_file(candidate_file_path)
        )
    )

def _remove_fasta_filename_suffix(relative_file_path: str) -> str:
    """Remove a recognized FASTA suffix from a relative path."""
    lowercase_relative_file_path = relative_file_path.lower()
    for filename_suffix in sorted(FASTA_FILENAME_SUFFIXES, key=len, reverse=True):
        if lowercase_relative_file_path.endswith(filename_suffix):
            return relative_file_path[: -len(filename_suffix)]
    return relative_file_path


def _locus_name_from_file_path(
    schema_directory_path: Path,
    locus_file_path: Path,
) -> str:
    """Return the schema-relative locus name without its FASTA suffix."""
    relative_file_path = locus_file_path.relative_to(schema_directory_path).as_posix()
    return _remove_fasta_filename_suffix(relative_file_path)

def _iterate_fasta_records(fasta_lines: Iterable[str]) -> Iterable[_FastaRecord]:
    """Yield normalized identifiers and sequences from FASTA text."""
    current_identifier: str | None = None
    current_sequence_parts: list[str] = []
    for fasta_line in fasta_lines:
        stripped_fasta_line = fasta_line.strip()
        if stripped_fasta_line.startswith(">"):
            if current_identifier is not None:
                yield _FastaRecord(
                    current_identifier,
                    "".join(current_sequence_parts).upper(),
                )
            current_identifier = stripped_fasta_line[1:].split(maxsplit=1)[0]
            current_sequence_parts = []
        elif stripped_fasta_line and current_identifier is not None:
            current_sequence_parts.append(stripped_fasta_line)

    if current_identifier is not None:
        yield _FastaRecord(
            current_identifier,
            "".join(current_sequence_parts).upper(),
        )


def _read_fasta_records(fasta_file_path: Path) -> tuple[_FastaRecord, ...]:
    """Read normalized FASTA records from plain text or gzip input."""
    if fasta_file_path.name.lower().endswith(".gz"):
        with gzip.open(fasta_file_path, mode="rt", encoding="utf-8") as fasta_file:
            return tuple(_iterate_fasta_records(fasta_file))
    with fasta_file_path.open(encoding="utf-8") as fasta_file:
        return tuple(_iterate_fasta_records(fasta_file))


def _is_valid_coding_sequence(normalized_sequence: str) -> bool:
    """Reuse the repository's four coding-sequence validation rules."""
    validation = _evaluate_allele_sequence(normalized_sequence)
    return not any(
        (
            validation.has_incomplete_orf,
            validation.has_ambiguous_bases,
            validation.is_missing_start_or_stop_codon,
            validation.has_in_frame_stop_codon,
        )
    )


def _is_eligible_nucleotide_sequence(normalized_sequence: str) -> bool:
    """Return whether a non-empty sequence contains standard IUPAC DNA."""
    return bool(normalized_sequence) and set(normalized_sequence).issubset(
        IUPAC_DNA_BASES
    )


def _create_selected_allele_sequence(
    *,
    sequence_index: int,
    sequence_level: str,
    locus_name: str,
    fasta_record: _FastaRecord,
) -> _SelectedAlleleSequence:
    """Create a BLAST-safe selected allele while preserving its original ID."""
    if sequence_level == SEQUENCE_LEVEL_AMINO_ACID:
        sequence = _translate_valid_allele(fasta_record.sequence)
    else:
        sequence = fasta_record.sequence
    return _SelectedAlleleSequence(
        internal_identifier=f"sequence_{sequence_index:09d}",
        locus_name=locus_name,
        original_identifier=fasta_record.identifier,
        sequence=sequence,
    )


def _select_allele_closest_to_modal_length(
    locus_records: tuple[_FastaRecord, ...],
    *,
    require_valid_coding_sequence: bool,
) -> _FastaRecord | None:
    """Select an eligible allele closest to the shortest tied length mode.

    The mode uses every non-empty allele in the locus, including records that
    are not valid coding sequences. The selected record must either pass all
    coding-sequence rules or, for the nucleotide-only fallback, contain only
    standard IUPAC DNA symbols. Equal distances retain original FASTA order.
    """
    non_empty_allele_lengths = [
        len(record.sequence) for record in locus_records if record.sequence
    ]
    if not non_empty_allele_lengths:
        return None
    length_counts = Counter(non_empty_allele_lengths)
    highest_length_frequency = max(length_counts.values())
    modal_allele_length = min(
        allele_length
        for allele_length, frequency in length_counts.items()
        if frequency == highest_length_frequency
    )
    for _, fasta_record in sorted(
        enumerate(locus_records),
        key=lambda indexed_record: (
            abs(len(indexed_record[1].sequence) - modal_allele_length),
            indexed_record[0],
        ),
    ):
        is_eligible = (
            _is_valid_coding_sequence(fasta_record.sequence)
            if require_valid_coding_sequence
            else _is_eligible_nucleotide_sequence(fasta_record.sequence)
        )
        if is_eligible:
            return fasta_record
    return None


def _collect_selected_allele_sequences(
    schema_directory_path: Path,
    all_fasta_file_paths: tuple[Path, ...],
) -> tuple[
    tuple[_SelectedAlleleSequence, ...],
    tuple[_SelectedAlleleSequence, ...],
    int,
]:
    """Select one modal-length-proximal allele per locus and sequence level."""
    selected_nucleotide_alleles: list[_SelectedAlleleSequence] = []
    selected_amino_acid_alleles: list[_SelectedAlleleSequence] = []
    sequence_index = 0

    for locus_file_path in all_fasta_file_paths:
        locus_name = _locus_name_from_file_path(
            schema_directory_path,
            locus_file_path,
        )
        locus_records = _read_fasta_records(locus_file_path)
        selected_valid_record = _select_allele_closest_to_modal_length(
            locus_records,
            require_valid_coding_sequence=True,
        )
        selected_nucleotide_record = selected_valid_record
        if selected_nucleotide_record is None:
            selected_nucleotide_record = _select_allele_closest_to_modal_length(
                locus_records,
                require_valid_coding_sequence=False,
            )

        if selected_nucleotide_record is not None:
            sequence_index += 1
            selected_nucleotide_alleles.append(
                _create_selected_allele_sequence(
                    sequence_index=sequence_index,
                    sequence_level=SEQUENCE_LEVEL_NUCLEOTIDE,
                    locus_name=locus_name,
                    fasta_record=selected_nucleotide_record,
                )
            )
        if selected_valid_record is not None:
            sequence_index += 1
            selected_amino_acid_alleles.append(
                _create_selected_allele_sequence(
                    sequence_index=sequence_index,
                    sequence_level=SEQUENCE_LEVEL_AMINO_ACID,
                    locus_name=locus_name,
                    fasta_record=selected_valid_record,
                )
            )

    return (
        tuple(selected_nucleotide_alleles),
        tuple(selected_amino_acid_alleles),
        len(all_fasta_file_paths),
    )


def _create_source_signature(
    schema_directory_path: Path,
    fasta_file_paths: Iterable[Path],
) -> _SourceSignature:
    """Create a path, size, and modification-time FASTA signature."""
    return tuple(
        (
            fasta_file_path.relative_to(schema_directory_path).as_posix(),
            fasta_file_path.stat().st_size,
            fasta_file_path.stat().st_mtime_ns,
        )
        for fasta_file_path in fasta_file_paths
    )


def _run_command(command_arguments: list[str], operation_name: str) -> None:
    """Run one BLAST+ command and raise a concise error on failure."""
    try:
        completed_process = subprocess.run(
            command_arguments,
            check=False,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as error:
        raise RuntimeError(
            f"{operation_name} executable was not found: {command_arguments[0]}"
        ) from error
    if completed_process.returncode != 0:
        error_message = completed_process.stderr.strip() or completed_process.stdout.strip()
        raise RuntimeError(
            f"{operation_name} failed with exit code "
            f"{completed_process.returncode}: {error_message}"
        )


def _get_blast_version(blast_executable: str) -> str:
    """Return the first version line emitted by a BLAST executable."""
    try:
        completed_process = subprocess.run(
            [blast_executable, "-version"],
            check=False,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as error:
        raise RuntimeError(f"BLAST executable was not found: {blast_executable}") from error
    if completed_process.returncode != 0:
        error_message = completed_process.stderr.strip() or completed_process.stdout.strip()
        raise RuntimeError(
            f"Could not determine the version of {blast_executable}: {error_message}"
        )
    version_lines = completed_process.stdout.strip().splitlines()
    return version_lines[0] if version_lines else f"{blast_executable} version unknown"


def _write_selected_allele_fasta(
    output_file_path: Path,
    selected_alleles: tuple[_SelectedAlleleSequence, ...],
) -> None:
    """Write internal BLAST-safe IDs and selected allele sequences."""
    with output_file_path.open("w", encoding="utf-8") as output_file:
        for selected_allele in selected_alleles:
            output_file.write(
                f">{selected_allele.internal_identifier}\n"
                f"{selected_allele.sequence}\n"
            )


def _parse_blast_output(
    blast_output_file_path: Path,
    *,
    sequence_level: str,
) -> tuple[_BlastHsp, ...]:
    """Parse the configured twelve-column BLAST tabular output."""
    parsed_hsps: list[_BlastHsp] = []
    with blast_output_file_path.open(encoding="utf-8") as blast_output_file:
        for line_number, output_line in enumerate(blast_output_file, start=1):
            if not output_line.strip():
                continue
            output_fields = output_line.rstrip("\n").split("\t")
            if len(output_fields) != 12:
                raise RuntimeError(
                    "Unexpected BLAST output at line "
                    f"{line_number}: expected 12 fields, found {len(output_fields)}"
                )
            (
                query_identifier,
                subject_identifier,
                query_length,
                subject_length,
                query_start,
                query_end,
                subject_start,
                subject_end,
                alignment_length,
                identical_character_count,
                raw_score,
                subject_strand,
            ) = output_fields
            parsed_hsps.append(
                _BlastHsp(
                    query_identifier=query_identifier,
                    subject_identifier=subject_identifier,
                    query_length=int(query_length),
                    subject_length=int(subject_length),
                    query_start=int(query_start),
                    query_end=int(query_end),
                    subject_start=int(subject_start),
                    subject_end=int(subject_end),
                    alignment_length=int(alignment_length),
                    identical_character_count=int(identical_character_count),
                    raw_score=float(raw_score),
                    subject_strand=(
                        subject_strand
                        if sequence_level == SEQUENCE_LEVEL_NUCLEOTIDE
                        else None
                    ),
                )
            )
    return tuple(parsed_hsps)


def _create_similarity_matches(
    *,
    sequence_level: str,
    selected_alleles: tuple[_SelectedAlleleSequence, ...],
    blast_hsps: tuple[_BlastHsp, ...],
) -> tuple[CrossLocusSimilarityMatch, ...]:
    """Calculate directional BSR and retain the best HSP per locus pair."""
    selected_alleles_by_identifier = {
        selected_allele.internal_identifier: selected_allele
        for selected_allele in selected_alleles
    }
    self_scores: dict[str, float] = {}
    for blast_hsp in blast_hsps:
        if blast_hsp.query_identifier == blast_hsp.subject_identifier:
            self_scores[blast_hsp.query_identifier] = max(
                blast_hsp.raw_score,
                self_scores.get(blast_hsp.query_identifier, 0.0),
            )

    best_matches_by_locus_pair: dict[
        tuple[str, str], CrossLocusSimilarityMatch
    ] = {}
    for blast_hsp in blast_hsps:
        query_selected_allele = selected_alleles_by_identifier.get(
            blast_hsp.query_identifier
        )
        subject_selected_allele = selected_alleles_by_identifier.get(
            blast_hsp.subject_identifier
        )
        if query_selected_allele is None or subject_selected_allele is None:
            continue
        if query_selected_allele.locus_name == subject_selected_allele.locus_name:
            continue
        query_self_score = self_scores.get(blast_hsp.query_identifier)
        if query_self_score is None or query_self_score <= 0:
            continue

        identity_percentage = (
            blast_hsp.identical_character_count
            / blast_hsp.alignment_length
            * 100
            if blast_hsp.alignment_length
            else 0.0
        )
        query_coverage_percentage = (
            abs(blast_hsp.query_end - blast_hsp.query_start) + 1
        ) / blast_hsp.query_length * 100
        subject_coverage_percentage = (
            abs(blast_hsp.subject_end - blast_hsp.subject_start) + 1
        ) / blast_hsp.subject_length * 100
        similarity_match = CrossLocusSimilarityMatch(
            sequence_level=sequence_level,
            query_locus_name=query_selected_allele.locus_name,
            subject_locus_name=subject_selected_allele.locus_name,
            query_allele_identifier=query_selected_allele.original_identifier,
            subject_allele_identifier=subject_selected_allele.original_identifier,
            identity_percentage=identity_percentage,
            alignment_length=blast_hsp.alignment_length,
            query_coverage_percentage=query_coverage_percentage,
            subject_coverage_percentage=subject_coverage_percentage,
            minimum_coverage_percentage=min(
                query_coverage_percentage,
                subject_coverage_percentage,
            ),
            raw_score=blast_hsp.raw_score,
            blast_score_ratio=blast_hsp.raw_score / query_self_score,
            nucleotide_subject_strand=blast_hsp.subject_strand,
        )
        locus_pair = (
            query_selected_allele.locus_name,
            subject_selected_allele.locus_name,
        )
        current_best_match = best_matches_by_locus_pair.get(locus_pair)
        candidate_sort_key = (
            similarity_match.raw_score,
            similarity_match.blast_score_ratio,
            similarity_match.minimum_coverage_percentage,
            similarity_match.identity_percentage,
            similarity_match.query_allele_identifier,
            similarity_match.subject_allele_identifier,
        )
        current_sort_key = (
            current_best_match.raw_score,
            current_best_match.blast_score_ratio,
            current_best_match.minimum_coverage_percentage,
            current_best_match.identity_percentage,
            current_best_match.query_allele_identifier,
            current_best_match.subject_allele_identifier,
        ) if current_best_match is not None else None
        if current_sort_key is None or candidate_sort_key > current_sort_key:
            best_matches_by_locus_pair[locus_pair] = similarity_match

    return tuple(
        best_matches_by_locus_pair[locus_pair]
        for locus_pair in sorted(best_matches_by_locus_pair)
    )


def _run_similarity_level(
    *,
    sequence_level: str,
    selected_alleles: tuple[_SelectedAlleleSequence, ...],
    working_directory_path: Path,
    blast_executable: str,
    makeblastdb_executable: str,
    blast_evalue: float,
    thread_count: int,
) -> _SimilarityLevelResult:
    """Build one temporary database and run one batched BLAST search."""
    if not selected_alleles:
        return _SimilarityLevelResult(sequence_level, ())

    database_type = (
        "nucl" if sequence_level == SEQUENCE_LEVEL_NUCLEOTIDE else "prot"
    )
    input_fasta_file_path = working_directory_path / f"{sequence_level}.fasta"
    database_path = working_directory_path / f"{sequence_level}_database"
    blast_output_file_path = working_directory_path / f"{sequence_level}.tsv"
    _write_selected_allele_fasta(input_fasta_file_path, selected_alleles)
    _run_command(
        [
            makeblastdb_executable,
            "-in",
            str(input_fasta_file_path),
            "-dbtype",
            database_type,
            "-out",
            str(database_path),
            "-parse_seqids",
        ],
        f"makeblastdb ({sequence_level})",
    )
    _run_command(
        [
            blast_executable,
            "-query",
            str(input_fasta_file_path),
            "-db",
            str(database_path),
            "-out",
            str(blast_output_file_path),
            "-outfmt",
            (
                "6 qseqid sseqid qlen slen qstart qend sstart send length "
                "nident score sstrand"
            ),
            "-task",
            (
                "blastn"
                if sequence_level == SEQUENCE_LEVEL_NUCLEOTIDE
                else "blastp"
            ),
            "-evalue",
            str(blast_evalue),
            "-num_threads",
            str(thread_count),
            "-max_target_seqs",
            str(len(selected_alleles)),
            "-max_hsps",
            "1",
        ],
        f"BLAST ({sequence_level})",
    )
    blast_hsps = _parse_blast_output(
        blast_output_file_path,
        sequence_level=sequence_level,
    )
    return _SimilarityLevelResult(
        sequence_level=sequence_level,
        matches=_create_similarity_matches(
            sequence_level=sequence_level,
            selected_alleles=selected_alleles,
            blast_hsps=blast_hsps,
        ),
    )


def _deserialize_similarity_match(
    match_data: object,
) -> CrossLocusSimilarityMatch:
    """Convert one validated cache dictionary into a typed match."""
    if not isinstance(match_data, dict):
        raise TypeError("Cached similarity match is not an object")
    return CrossLocusSimilarityMatch(
        sequence_level=str(match_data["sequence_level"]),
        query_locus_name=str(match_data["query_locus_name"]),
        subject_locus_name=str(match_data["subject_locus_name"]),
        query_allele_identifier=str(match_data["query_allele_identifier"]),
        subject_allele_identifier=str(match_data["subject_allele_identifier"]),
        identity_percentage=float(match_data["identity_percentage"]),
        alignment_length=int(match_data["alignment_length"]),
        query_coverage_percentage=float(match_data["query_coverage_percentage"]),
        subject_coverage_percentage=float(
            match_data["subject_coverage_percentage"]
        ),
        minimum_coverage_percentage=float(
            match_data["minimum_coverage_percentage"]
        ),
        raw_score=float(match_data["raw_score"]),
        blast_score_ratio=float(match_data["blast_score_ratio"]),
        nucleotide_subject_strand=(
            str(match_data["nucleotide_subject_strand"])
            if match_data.get("nucleotide_subject_strand") is not None
            else None
        ),
    )


def _read_cached_analysis(
    *,
    schema_directory_path: Path,
    source_signature: _SourceSignature,
    blastn_version: str | None,
    blastp_version: str | None,
    blast_evalue: float,
) -> SchemaSimilarityAnalysis | None:
    """Read a compatible similarity cache or return ``None``."""
    cache_file_path = schema_directory_path / SCHEMA_SIMILARITY_CACHE_FILENAME
    try:
        cache_data = json.loads(cache_file_path.read_text(encoding="utf-8"))
        if int(cache_data["cache_version"]) != SCHEMA_SIMILARITY_CACHE_VERSION:
            return None
        if str(cache_data["method_identifier"]) != SIMILARITY_METHOD_IDENTIFIER:
            return None
        cached_blastn_version = str(cache_data["blastn_version"])
        cached_blastp_version = str(cache_data["blastp_version"])
        if (
            blastn_version is not None
            and cached_blastn_version != blastn_version
        ):
            return None
        if (
            blastp_version is not None
            and cached_blastp_version != blastp_version
        ):
            return None
        if float(cache_data["blast_evalue"]) != blast_evalue:
            return None
        cached_source_signature = tuple(
            (str(path), int(file_size), int(modification_time))
            for path, file_size, modification_time in cache_data["source_signature"]
        )
        if cached_source_signature != source_signature:
            return None
        return SchemaSimilarityAnalysis(
            locus_count=int(cache_data["locus_count"]),
            nucleotide_selected_locus_count=int(
                cache_data["nucleotide_selected_locus_count"]
            ),
            amino_acid_selected_locus_count=int(
                cache_data["amino_acid_selected_locus_count"]
            ),
            nucleotide_possible_directed_pair_count=int(
                cache_data["nucleotide_possible_directed_pair_count"]
            ),
            amino_acid_possible_directed_pair_count=int(
                cache_data["amino_acid_possible_directed_pair_count"]
            ),
            matches=tuple(
                _deserialize_similarity_match(match_data)
                for match_data in cache_data["matches"]
            ),
            blastn_version=cached_blastn_version,
            blastp_version=cached_blastp_version,
            method_identifier=SIMILARITY_METHOD_IDENTIFIER,
            blast_evalue=blast_evalue,
            reused_from_cache=True,
        )
    except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError):
        return None


def _write_cached_analysis(
    *,
    schema_directory_path: Path,
    source_signature: _SourceSignature,
    analysis: SchemaSimilarityAnalysis,
) -> None:
    """Atomically persist one similarity analysis when the schema is writable."""
    cache_data = {
        "cache_version": SCHEMA_SIMILARITY_CACHE_VERSION,
        "method_identifier": analysis.method_identifier,
        "blastn_version": analysis.blastn_version,
        "blastp_version": analysis.blastp_version,
        "blast_evalue": analysis.blast_evalue,
        "source_signature": source_signature,
        "locus_count": analysis.locus_count,
        "nucleotide_selected_locus_count": (
            analysis.nucleotide_selected_locus_count
        ),
        "amino_acid_selected_locus_count": (
            analysis.amino_acid_selected_locus_count
        ),
        "nucleotide_possible_directed_pair_count": (
            analysis.nucleotide_possible_directed_pair_count
        ),
        "amino_acid_possible_directed_pair_count": (
            analysis.amino_acid_possible_directed_pair_count
        ),
        "matches": [asdict(match) for match in analysis.matches],
    }
    temporary_cache_file_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=schema_directory_path,
            prefix=f"{SCHEMA_SIMILARITY_CACHE_FILENAME}.",
            suffix=".tmp",
            delete=False,
        ) as temporary_cache_file:
            json.dump(cache_data, temporary_cache_file, separators=(",", ":"))
            temporary_cache_file_path = Path(temporary_cache_file.name)
        temporary_cache_file_path.replace(
            schema_directory_path / SCHEMA_SIMILARITY_CACHE_FILENAME
        )
    except OSError:
        if temporary_cache_file_path is not None:
            temporary_cache_file_path.unlink(missing_ok=True)


def _resolve_maximum_threads(maximum_threads: int | None) -> int:
    """Resolve a positive thread budget from a caller value or the machine."""
    if maximum_threads is not None:
        if maximum_threads < 1:
            raise ValueError("maximum_threads must be at least 1")
        return maximum_threads
    return max(1, os.process_cpu_count() or 1)


def _emit_progress(
    progress: SimilarityProgress,
    *,
    progress_callback: SimilarityProgressCallback | None,
    progress_output_stream: TextIO | None,
) -> None:
    """Send a progress update to the notebook callback and terminal stream."""
    if progress_callback is not None:
        progress_callback(progress)
    if progress_output_stream is not None:
        print(
            f"[BLAST similarity {progress.completed_stage_count}/"
            f"{progress.total_stage_count}] {progress.message}",
            file=progress_output_stream,
            flush=True,
        )


def get_schema_similarity_analysis(
    schema_directory_path: str | Path,
    *,
    maximum_threads: int | None = None,
    blast_evalue: float = DEFAULT_BLAST_EVALUE,
    blastn_executable: str = "blastn",
    blastp_executable: str = "blastp",
    makeblastdb_executable: str = "makeblastdb",
    progress_callback: SimilarityProgressCallback | None = None,
    progress_output_stream: TextIO | None = None,
) -> SchemaSimilarityAnalysis:
    """Measure selected-allele nucleotide and protein cross-locus hits.

    Args:
        schema_directory_path: Extracted scheme directory containing one locus
            FASTA per target.
        maximum_threads: Total BLAST thread budget. Defaults to all logical
            CPUs available to the process and is split between simultaneous
            nucleotide and protein searches.
        blast_evalue: BLAST reporting floor. This controls search sensitivity;
            it is not a candidate-paralogy classification threshold.
        blastn_executable: ``blastn`` executable name or path.
        blastp_executable: ``blastp`` executable name or path.
        makeblastdb_executable: ``makeblastdb`` executable name or path.
        progress_callback: Optional typed callback for notebook progress.
        progress_output_stream: Optional terminal stream for progress lines.

    Returns:
        A cache-aware :class:`SchemaSimilarityAnalysis` containing only
        directed cross-locus local hits reported by BLAST.
    """
    if blast_evalue <= 0:
        raise ValueError("blast_evalue must be greater than zero")
    resolved_schema_directory_path = _resolve_schema_directory_path(
        schema_directory_path
    )
    all_fasta_file_paths = _find_all_fasta_file_paths(
        resolved_schema_directory_path
    )
    source_signature = _create_source_signature(
        resolved_schema_directory_path,
        all_fasta_file_paths,
    )
    try:
        blastn_version = _get_blast_version(blastn_executable)
        blastp_version = _get_blast_version(blastp_executable)
    except RuntimeError:
        cached_analysis_without_installed_blast = _read_cached_analysis(
            schema_directory_path=resolved_schema_directory_path,
            source_signature=source_signature,
            blastn_version=None,
            blastp_version=None,
            blast_evalue=blast_evalue,
        )
        if cached_analysis_without_installed_blast is None:
            raise
        _emit_progress(
            SimilarityProgress(
                completed_stage_count=SIMILARITY_PROGRESS_STAGE_COUNT,
                total_stage_count=SIMILARITY_PROGRESS_STAGE_COUNT,
                stage_name="cache",
                message=(
                    "Reused the complete cross-locus similarity cache; "
                    "BLAST executables were not needed."
                ),
                reused_from_cache=True,
            ),
            progress_callback=progress_callback,
            progress_output_stream=progress_output_stream,
        )
        return cached_analysis_without_installed_blast
    cached_analysis = _read_cached_analysis(
        schema_directory_path=resolved_schema_directory_path,
        source_signature=source_signature,
        blastn_version=blastn_version,
        blastp_version=blastp_version,
        blast_evalue=blast_evalue,
    )
    if cached_analysis is not None:
        _emit_progress(
            SimilarityProgress(
                completed_stage_count=SIMILARITY_PROGRESS_STAGE_COUNT,
                total_stage_count=SIMILARITY_PROGRESS_STAGE_COUNT,
                stage_name="cache",
                message="Reused the complete cross-locus similarity cache.",
                reused_from_cache=True,
            ),
            progress_callback=progress_callback,
            progress_output_stream=progress_output_stream,
        )
        return cached_analysis

    (
        selected_nucleotide_alleles,
        selected_amino_acid_alleles,
        locus_count,
    ) = _collect_selected_allele_sequences(
        resolved_schema_directory_path,
        all_fasta_file_paths,
    )
    _emit_progress(
        SimilarityProgress(
            completed_stage_count=1,
            total_stage_count=SIMILARITY_PROGRESS_STAGE_COUNT,
            stage_name="selected_alleles",
            message=(
                f"Selected modal-length-proximal alleles for {locus_count:,} loci "
                f"({len(selected_nucleotide_alleles):,} nucleotide; "
                f"{len(selected_amino_acid_alleles):,} amino-acid sequences)."
            ),
        ),
        progress_callback=progress_callback,
        progress_output_stream=progress_output_stream,
    )

    maximum_thread_count = _resolve_maximum_threads(maximum_threads)
    available_levels = sum(
        bool(selected_alleles)
        for selected_alleles in (
            selected_nucleotide_alleles,
            selected_amino_acid_alleles,
        )
    )
    if available_levels > 1:
        nucleotide_thread_count = max(1, maximum_thread_count // 2)
        amino_acid_thread_count = max(
            1,
            maximum_thread_count - nucleotide_thread_count,
        )
    else:
        nucleotide_thread_count = maximum_thread_count
        amino_acid_thread_count = maximum_thread_count

    level_results: dict[str, _SimilarityLevelResult] = {}
    with tempfile.TemporaryDirectory(prefix="schema-similarity-") as working_directory:
        working_directory_path = Path(working_directory)
        level_jobs = (
            (
                SEQUENCE_LEVEL_NUCLEOTIDE,
                selected_nucleotide_alleles,
                blastn_executable,
                nucleotide_thread_count,
            ),
            (
                SEQUENCE_LEVEL_AMINO_ACID,
                selected_amino_acid_alleles,
                blastp_executable,
                amino_acid_thread_count,
            ),
        )
        with ThreadPoolExecutor(max_workers=max(1, available_levels)) as executor:
            future_to_sequence_level = {
                executor.submit(
                    _run_similarity_level,
                    sequence_level=sequence_level,
                    selected_alleles=selected_alleles,
                    working_directory_path=working_directory_path,
                    blast_executable=blast_executable,
                    makeblastdb_executable=makeblastdb_executable,
                    blast_evalue=blast_evalue,
                    thread_count=thread_count,
                ): sequence_level
                for (
                    sequence_level,
                    selected_alleles,
                    blast_executable,
                    thread_count,
                ) in level_jobs
            }
            for completed_level_count, completed_future in enumerate(
                as_completed(future_to_sequence_level),
                start=1,
            ):
                sequence_level = future_to_sequence_level[completed_future]
                level_results[sequence_level] = completed_future.result()
                _emit_progress(
                    SimilarityProgress(
                        completed_stage_count=1 + completed_level_count,
                        total_stage_count=SIMILARITY_PROGRESS_STAGE_COUNT,
                        stage_name=sequence_level,
                        message=(
                            f"Completed {sequence_level.replace('_', '-')} BLAST "
                            f"with {len(level_results[sequence_level].matches):,} "
                            "directed cross-locus hits."
                        ),
                    ),
                    progress_callback=progress_callback,
                    progress_output_stream=progress_output_stream,
                )

    matches = tuple(
        sorted(
            (
                *level_results[SEQUENCE_LEVEL_NUCLEOTIDE].matches,
                *level_results[SEQUENCE_LEVEL_AMINO_ACID].matches,
            ),
            key=lambda match: (
                match.sequence_level,
                match.query_locus_name,
                match.subject_locus_name,
            ),
        )
    )
    _emit_progress(
        SimilarityProgress(
            completed_stage_count=4,
            total_stage_count=SIMILARITY_PROGRESS_STAGE_COUNT,
            stage_name="results",
            message=f"Assembled {len(matches):,} cross-locus similarity rows.",
        ),
        progress_callback=progress_callback,
        progress_output_stream=progress_output_stream,
    )
    nucleotide_selected_locus_count = len(
        {
            selected_allele.locus_name
            for selected_allele in selected_nucleotide_alleles
        }
    )
    amino_acid_selected_locus_count = len(
        {
            selected_allele.locus_name
            for selected_allele in selected_amino_acid_alleles
        }
    )
    analysis = SchemaSimilarityAnalysis(
        locus_count=locus_count,
        nucleotide_selected_locus_count=nucleotide_selected_locus_count,
        amino_acid_selected_locus_count=amino_acid_selected_locus_count,
        nucleotide_possible_directed_pair_count=(
            nucleotide_selected_locus_count
            * max(0, nucleotide_selected_locus_count - 1)
        ),
        amino_acid_possible_directed_pair_count=(
            amino_acid_selected_locus_count
            * max(0, amino_acid_selected_locus_count - 1)
        ),
        matches=matches,
        blastn_version=blastn_version,
        blastp_version=blastp_version,
        method_identifier=SIMILARITY_METHOD_IDENTIFIER,
        blast_evalue=blast_evalue,
        reused_from_cache=False,
    )
    _write_cached_analysis(
        schema_directory_path=resolved_schema_directory_path,
        source_signature=source_signature,
        analysis=analysis,
    )
    _emit_progress(
        SimilarityProgress(
            completed_stage_count=SIMILARITY_PROGRESS_STAGE_COUNT,
            total_stage_count=SIMILARITY_PROGRESS_STAGE_COUNT,
            stage_name="cache",
            message="Cross-locus similarity analysis complete and cached.",
        ),
        progress_callback=progress_callback,
        progress_output_stream=progress_output_stream,
    )
    return replace(analysis, reused_from_cache=False)


def _parse_positive_integer(argument_value: str) -> int:
    """Parse a positive command-line integer."""
    parsed_value = int(argument_value)
    if parsed_value < 1:
        raise argparse.ArgumentTypeError("value must be at least 1")
    return parsed_value


def main() -> None:
    """Run cross-locus similarity analysis from the command line."""
    argument_parser = argparse.ArgumentParser(description=__doc__)
    argument_parser.add_argument(
        "schema_directory",
        type=Path,
        help="Extracted schema directory to analyze.",
    )
    argument_parser.add_argument(
        "--threads",
        type=_parse_positive_integer,
        default=None,
        help="Total BLAST thread budget (default: all available CPUs).",
    )
    argument_parser.add_argument(
        "--json",
        action="store_true",
        help="Emit the complete analysis as JSON.",
    )
    parsed_arguments = argument_parser.parse_args()
    analysis = get_schema_similarity_analysis(
        parsed_arguments.schema_directory,
        maximum_threads=parsed_arguments.threads,
        progress_output_stream=sys.stderr,
    )
    if parsed_arguments.json:
        print(json.dumps(asdict(analysis), indent=2))
        return
    print(f"Loci: {analysis.locus_count}")
    print(
        "Loci with selected nucleotide alleles: "
        f"{analysis.nucleotide_selected_locus_count}"
    )
    print(
        "Loci with selected amino-acid alleles: "
        f"{analysis.amino_acid_selected_locus_count}"
    )
    print(f"Nucleotide directed hits: {analysis.nucleotide_match_count}")
    print(f"Amino-acid directed hits: {analysis.amino_acid_match_count}")
    print(f"Cache reused: {analysis.reused_from_cache}")


if __name__ == "__main__":
    main()
