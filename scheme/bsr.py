"""Measure representative-to-allele protein BSR coverage for a scheme.

This module implements the protein-level representative validation used by
Schema Refinery without requiring a curated ``short/`` directory.  The current
scheme collection does not contain curated representatives, so one valid coding
allele is selected per locus with the repository's existing deterministic rule:
the allele closest to the shortest tied modal sequence length, with FASTA order
breaking equal-distance ties.

Every valid coding allele is translated and placed in one BLASTp database.  The
selected representatives are searched against that database in a single
multi-query run.  Directional BLAST Score Ratio (BSR) is the best local raw
score for a representative-to-allele match divided by that representative's
self-alignment raw score.  Missing BLAST hits remain distinguishable from
reported low-BSR hits.

The raw, threshold-independent measurements support two views.  Within-locus
coverage reports how many valid translated alleles are represented at a chosen
BSR threshold.  Cross-locus matches are reduced to the strongest directional
match for each locus pair, collapsed into unordered pairs, and connected into
high-BSR similarity groups.  These groups are sequence-similarity observations,
not confirmation of copy number or evolutionary paralogy.

Results are cached atomically inside the extracted scheme using a FASTA source
signature, BLASTp version, selection method, and E-value.  Changing the display
threshold reclassifies cached measurements without rerunning BLAST or rereading
the FASTA files.  The public API is independent of Marimo and includes progress
callbacks plus a command-line interface.
"""

import argparse
import json
import sys
import tempfile
from collections.abc import Callable, Iterable
from dataclasses import asdict, dataclass, replace
from math import isfinite
from pathlib import Path
from typing import TextIO

from .blast import create_blast_database, get_blast_version, run_blast_command
from .cache_io import write_json_atomically
from .fasta import (
    FastaRecord,
    create_fasta_source_signature,
    find_locus_fasta_paths,
    locus_name_from_file_path,
    read_fasta_records,
    resolve_scheme_directory_path,
)
from .resources import resolve_worker_count
from .selection import select_allele_closest_to_modal_length
from .sequence import evaluate_allele_sequence, translate_valid_allele
from .similarity import (
    DEFAULT_BLAST_EVALUE,
    DEFAULT_SAME_GENE_BLAST_SCORE_RATIO_THRESHOLD,
)

# Private aliases preserve existing test and downstream patch points while the
# implementation now lives in the shared operational module.
_run_command = run_blast_command
_get_blast_version = get_blast_version

SCHEME_REPRESENTATIVE_BSR_CACHE_FILENAME = ".scheme-representative-bsr-cache.json"
LEGACY_REPRESENTATIVE_BSR_CACHE_FILENAME = (
    ".schema-representative-bsr-cache.json"
)
SCHEME_REPRESENTATIVE_BSR_CACHE_FILENAMES = (
    SCHEME_REPRESENTATIVE_BSR_CACHE_FILENAME,
    LEGACY_REPRESENTATIVE_BSR_CACHE_FILENAME,
)
SCHEME_REPRESENTATIVE_BSR_CACHE_VERSION = 1
REPRESENTATIVE_BSR_METHOD_IDENTIFIER = (
    "blastp-modal-length-representative-vs-valid-alleles-bsr-v1"
)
REPRESENTATIVE_BSR_PROGRESS_STAGE_COUNT = 5

__all__ = [
    "DEFAULT_BLAST_EVALUE",
    "DEFAULT_SAME_GENE_BLAST_SCORE_RATIO_THRESHOLD",
    "REPRESENTATIVE_BSR_METHOD_IDENTIFIER",
    "REPRESENTATIVE_BSR_PROGRESS_STAGE_COUNT",
    "SCHEME_REPRESENTATIVE_BSR_CACHE_FILENAME",
    "SCHEME_REPRESENTATIVE_BSR_CACHE_FILENAMES",
    "SCHEME_REPRESENTATIVE_BSR_CACHE_VERSION",
    "AlleleBlastScoreRatio",
    "HighBsrLocusGroup",
    "HighBsrLocusPair",
    "InterLocusBlastScoreRatio",
    "LocusRepresentativeBsrCoverage",
    "LocusRepresentativeBsrMeasurements",
    "RepresentativeBsrClassification",
    "RepresentativeBsrProgress",
    "SchemeRepresentativeBsrAnalysis",
    "classify_representative_bsr_analysis",
    "get_scheme_representative_bsr_analysis",
]

type _SourceSignature = tuple[tuple[str, int, int], ...]

_create_source_signature = create_fasta_source_signature
_find_all_fasta_file_paths = find_locus_fasta_paths
_locus_name_from_file_path = locus_name_from_file_path
_resolve_scheme_directory_path = resolve_scheme_directory_path
_translate_valid_allele = translate_valid_allele


def _read_fasta_records(fasta_file_path: Path) -> tuple[FastaRecord, ...]:
    """Read FASTA records with the established uppercase BLAST normalization."""
    return tuple(
        FastaRecord(record.identifier, record.sequence.upper())
        for record in read_fasta_records(fasta_file_path)
    )


def _is_valid_coding_sequence(normalized_sequence: str) -> bool:
    """Return whether a normalized allele passes every CDS rule."""
    return evaluate_allele_sequence(normalized_sequence).is_valid


def _select_allele_closest_to_modal_length(
    locus_records: tuple[FastaRecord, ...],
    *,
    require_valid_coding_sequence: bool,
) -> FastaRecord | None:
    """Select a record using the shared deterministic rule."""
    return select_allele_closest_to_modal_length(
        locus_records,
        require_valid_coding_sequence=require_valid_coding_sequence,
    )


@dataclass(frozen=True)
class AlleleBlastScoreRatio:
    """Directional BSR from a locus representative to one valid allele.

    ``blast_score_ratio`` is ``None`` when BLASTp reported no hit or when the
    representative self score was unavailable.  Those cases are not silently
    converted into a measured score of zero.
    """

    allele_identifier: str
    blast_score_ratio: float | None


@dataclass(frozen=True)
class LocusRepresentativeBsrMeasurements:
    """Threshold-independent representative-to-allele measurements per locus."""

    locus_name: str
    representative_allele_identifier: str | None
    representative_self_score: float | None
    total_allele_count: int
    valid_translated_allele_count: int
    allele_measurements: tuple[AlleleBlastScoreRatio, ...]

    @property
    def unevaluable_allele_count(self) -> int:
        """Return records excluded because they are not valid coding alleles."""
        return self.total_allele_count - self.valid_translated_allele_count


@dataclass(frozen=True)
class InterLocusBlastScoreRatio:
    """Strongest representative-to-allele hit for one directed locus pair."""

    query_locus_name: str
    subject_locus_name: str
    representative_allele_identifier: str
    subject_allele_identifier: str
    raw_score: float
    blast_score_ratio: float


@dataclass(frozen=True)
class SchemeRepresentativeBsrAnalysis:
    """Complete threshold-independent representative BSR analysis."""

    locus_measurements: tuple[LocusRepresentativeBsrMeasurements, ...]
    inter_locus_measurements: tuple[InterLocusBlastScoreRatio, ...]
    blastp_version: str
    method_identifier: str
    blast_evalue: float
    reused_from_cache: bool

    @property
    def locus_count(self) -> int:
        """Return the number of locus FASTA files in the analysis."""
        return len(self.locus_measurements)

    @property
    def representative_locus_count(self) -> int:
        """Return loci for which a valid representative could be selected."""
        return sum(
            measurement.representative_allele_identifier is not None
            for measurement in self.locus_measurements
        )

    @property
    def self_score_locus_count(self) -> int:
        """Return loci with a usable representative BLASTp self score."""
        return sum(
            measurement.representative_self_score is not None
            for measurement in self.locus_measurements
        )

    @property
    def total_allele_count(self) -> int:
        """Return all FASTA records, including invalid coding sequences."""
        return sum(
            measurement.total_allele_count
            for measurement in self.locus_measurements
        )

    @property
    def valid_translated_allele_count(self) -> int:
        """Return valid alleles included in protein-level BSR evaluation."""
        return sum(
            measurement.valid_translated_allele_count
            for measurement in self.locus_measurements
        )


@dataclass(frozen=True)
class LocusRepresentativeBsrCoverage:
    """Threshold-classified representative coverage for one locus."""

    locus_name: str
    representative_allele_identifier: str | None
    total_allele_count: int
    valid_translated_allele_count: int
    unevaluable_allele_count: int
    covered_allele_count: int
    below_threshold_hit_count: int
    no_reported_hit_count: int
    coverage_percentage: float | None


@dataclass(frozen=True)
class HighBsrLocusPair:
    """One unordered cross-locus pair reaching the selected BSR threshold."""

    first_locus_name: str
    second_locus_name: str
    evidence: InterLocusBlastScoreRatio


@dataclass(frozen=True)
class HighBsrLocusGroup:
    """Connected component of loci linked by qualifying high-BSR pairs."""

    group_number: int
    locus_names: tuple[str, ...]
    qualifying_pair_count: int
    maximum_blast_score_ratio: float


@dataclass(frozen=True)
class RepresentativeBsrClassification:
    """Within- and inter-locus results at one configurable BSR threshold."""

    blast_score_ratio_threshold: float
    locus_coverages: tuple[LocusRepresentativeBsrCoverage, ...]
    high_bsr_pairs: tuple[HighBsrLocusPair, ...]
    high_bsr_groups: tuple[HighBsrLocusGroup, ...]

    @property
    def measurable_allele_count(self) -> int:
        """Return valid alleles in loci with a usable representative score."""
        return sum(
            coverage.valid_translated_allele_count
            for coverage in self.locus_coverages
            if coverage.coverage_percentage is not None
        )

    @property
    def covered_allele_count(self) -> int:
        """Return measurable alleles reaching the selected BSR threshold."""
        return sum(
            coverage.covered_allele_count for coverage in self.locus_coverages
        )

    @property
    def covered_allele_percentage(self) -> float | None:
        """Return covered alleles as a percentage of measurable alleles."""
        if self.measurable_allele_count == 0:
            return None
        return self.covered_allele_count / self.measurable_allele_count * 100

    @property
    def unavailable_valid_allele_count(self) -> int:
        """Return valid alleles lacking a usable representative self score."""
        return sum(
            coverage.valid_translated_allele_count
            for coverage in self.locus_coverages
            if coverage.coverage_percentage is None
        )

    @property
    def incomplete_coverage_locus_count(self) -> int:
        """Return measurable loci with at least one uncovered valid allele."""
        return sum(
            coverage.coverage_percentage is not None
            and coverage.covered_allele_count
            < coverage.valid_translated_allele_count
            for coverage in self.locus_coverages
        )

    @property
    def high_bsr_locus_count(self) -> int:
        """Return unique loci occurring in at least one qualifying pair."""
        return len(
            {
                locus_name
                for pair in self.high_bsr_pairs
                for locus_name in (pair.first_locus_name, pair.second_locus_name)
            }
        )


@dataclass(frozen=True)
class RepresentativeBsrProgress:
    """One progress update emitted by the representative BSR analysis."""

    completed_stage_count: int
    total_stage_count: int
    stage_name: str
    message: str
    reused_from_cache: bool = False


type RepresentativeBsrProgressCallback = Callable[[RepresentativeBsrProgress], None]


@dataclass(frozen=True)
class _ProteinAllele:
    """One valid translated allele with a BLAST-safe internal identifier."""

    internal_identifier: str
    locus_name: str
    original_identifier: str
    protein_sequence: str


@dataclass(frozen=True)
class _CollectedLocus:
    """Private collection metadata retained before BLAST calculation."""

    locus_name: str
    representative_internal_identifier: str | None
    representative_allele_identifier: str | None
    total_allele_count: int
    protein_alleles: tuple[_ProteinAllele, ...]


def _collect_protein_alleles(
    scheme_directory_path: Path,
    fasta_file_paths: tuple[Path, ...],
) -> tuple[tuple[_CollectedLocus, ...], tuple[_ProteinAllele, ...]]:
    """Translate valid alleles and select one approximate representative."""
    collected_loci: list[_CollectedLocus] = []
    all_protein_alleles: list[_ProteinAllele] = []
    allele_index = 0
    for fasta_file_path in fasta_file_paths:
        locus_name = _locus_name_from_file_path(
            scheme_directory_path,
            fasta_file_path,
        )
        fasta_records = _read_fasta_records(fasta_file_path)
        representative_record = _select_allele_closest_to_modal_length(
            fasta_records,
            require_valid_coding_sequence=True,
        )
        representative_internal_identifier: str | None = None
        locus_protein_alleles: list[_ProteinAllele] = []
        for fasta_record in fasta_records:
            if not _is_valid_coding_sequence(fasta_record.sequence):
                continue
            allele_index += 1
            protein_allele = _ProteinAllele(
                internal_identifier=f"allele_{allele_index:012d}",
                locus_name=locus_name,
                original_identifier=fasta_record.identifier,
                protein_sequence=_translate_valid_allele(fasta_record.sequence),
            )
            locus_protein_alleles.append(protein_allele)
            all_protein_alleles.append(protein_allele)
            if fasta_record is representative_record:
                representative_internal_identifier = (
                    protein_allele.internal_identifier
                )
        collected_loci.append(
            _CollectedLocus(
                locus_name=locus_name,
                representative_internal_identifier=(
                    representative_internal_identifier
                ),
                representative_allele_identifier=(
                    representative_record.identifier
                    if representative_record is not None
                    else None
                ),
                total_allele_count=len(fasta_records),
                protein_alleles=tuple(locus_protein_alleles),
            )
        )
    return tuple(collected_loci), tuple(all_protein_alleles)


def _write_protein_fasta(
    output_file_path: Path,
    protein_alleles: Iterable[_ProteinAllele],
) -> None:
    """Write internal identifiers and protein sequences to a temporary FASTA."""
    with output_file_path.open("w", encoding="utf-8") as output_file:
        for protein_allele in protein_alleles:
            output_file.write(
                f">{protein_allele.internal_identifier}\n"
                f"{protein_allele.protein_sequence}\n"
            )


def _parse_score_fields(output_line: str, line_number: int) -> tuple[str, str, float]:
    """Parse one three-column BLAST score output row."""
    output_fields = output_line.rstrip("\n").split("\t")
    if len(output_fields) != 3:
        raise RuntimeError(
            "Unexpected BLASTp output at line "
            f"{line_number}: expected 3 fields, found {len(output_fields)}"
        )
    query_identifier, subject_identifier, raw_score = output_fields
    return query_identifier, subject_identifier, float(raw_score)


def _create_bsr_measurements(
    *,
    collected_loci: tuple[_CollectedLocus, ...],
    all_protein_alleles: tuple[_ProteinAllele, ...],
    blast_output_file_path: Path,
) -> tuple[
    tuple[LocusRepresentativeBsrMeasurements, ...],
    tuple[InterLocusBlastScoreRatio, ...],
]:
    """Stream BLAST output and reduce it to reusable BSR measurements."""
    protein_alleles_by_identifier = {
        protein_allele.internal_identifier: protein_allele
        for protein_allele in all_protein_alleles
    }
    collected_loci_by_representative_identifier = {
        locus.representative_internal_identifier: locus
        for locus in collected_loci
        if locus.representative_internal_identifier is not None
    }
    self_scores: dict[str, float] = {}
    with blast_output_file_path.open(encoding="utf-8") as blast_output_file:
        for line_number, output_line in enumerate(blast_output_file, start=1):
            if not output_line.strip():
                continue
            query_identifier, subject_identifier, raw_score = _parse_score_fields(
                output_line,
                line_number,
            )
            if query_identifier == subject_identifier:
                self_scores[query_identifier] = max(
                    raw_score,
                    self_scores.get(query_identifier, 0.0),
                )

    best_intra_locus_scores: dict[tuple[str, str], float] = {}
    best_inter_locus_measurements: dict[
        tuple[str, str], InterLocusBlastScoreRatio
    ] = {}
    with blast_output_file_path.open(encoding="utf-8") as blast_output_file:
        for line_number, output_line in enumerate(blast_output_file, start=1):
            if not output_line.strip():
                continue
            query_identifier, subject_identifier, raw_score = _parse_score_fields(
                output_line,
                line_number,
            )
            query_locus = collected_loci_by_representative_identifier.get(
                query_identifier
            )
            subject_allele = protein_alleles_by_identifier.get(subject_identifier)
            query_self_score = self_scores.get(query_identifier)
            if (
                query_locus is None
                or subject_allele is None
                or query_self_score is None
                or query_self_score <= 0
            ):
                continue
            blast_score_ratio = raw_score / query_self_score
            if query_locus.locus_name == subject_allele.locus_name:
                intra_locus_key = (
                    query_locus.locus_name,
                    subject_allele.internal_identifier,
                )
                best_intra_locus_scores[intra_locus_key] = max(
                    blast_score_ratio,
                    best_intra_locus_scores.get(intra_locus_key, 0.0),
                )
                continue
            inter_locus_measurement = InterLocusBlastScoreRatio(
                query_locus_name=query_locus.locus_name,
                subject_locus_name=subject_allele.locus_name,
                representative_allele_identifier=(
                    query_locus.representative_allele_identifier or ""
                ),
                subject_allele_identifier=subject_allele.original_identifier,
                raw_score=raw_score,
                blast_score_ratio=blast_score_ratio,
            )
            directed_locus_pair = (
                query_locus.locus_name,
                subject_allele.locus_name,
            )
            current_measurement = best_inter_locus_measurements.get(
                directed_locus_pair
            )
            candidate_sort_key = (
                inter_locus_measurement.blast_score_ratio,
                inter_locus_measurement.raw_score,
                inter_locus_measurement.representative_allele_identifier,
                inter_locus_measurement.subject_allele_identifier,
            )
            current_sort_key = (
                (
                    current_measurement.blast_score_ratio,
                    current_measurement.raw_score,
                    current_measurement.representative_allele_identifier,
                    current_measurement.subject_allele_identifier,
                )
                if current_measurement is not None
                else None
            )
            if current_sort_key is None or candidate_sort_key > current_sort_key:
                best_inter_locus_measurements[directed_locus_pair] = (
                    inter_locus_measurement
                )

    locus_measurements: list[LocusRepresentativeBsrMeasurements] = []
    for collected_locus in collected_loci:
        representative_identifier = (
            collected_locus.representative_internal_identifier
        )
        representative_self_score = (
            self_scores.get(representative_identifier)
            if representative_identifier is not None
            else None
        )
        allele_measurements = tuple(
            AlleleBlastScoreRatio(
                allele_identifier=protein_allele.original_identifier,
                blast_score_ratio=(
                    best_intra_locus_scores.get(
                        (
                            collected_locus.locus_name,
                            protein_allele.internal_identifier,
                        )
                    )
                    if representative_self_score is not None
                    else None
                ),
            )
            for protein_allele in collected_locus.protein_alleles
        )
        locus_measurements.append(
            LocusRepresentativeBsrMeasurements(
                locus_name=collected_locus.locus_name,
                representative_allele_identifier=(
                    collected_locus.representative_allele_identifier
                ),
                representative_self_score=representative_self_score,
                total_allele_count=collected_locus.total_allele_count,
                valid_translated_allele_count=len(
                    collected_locus.protein_alleles
                ),
                allele_measurements=allele_measurements,
            )
        )
    return (
        tuple(locus_measurements),
        tuple(
            best_inter_locus_measurements[directed_locus_pair]
            for directed_locus_pair in sorted(best_inter_locus_measurements)
        ),
    )


def _run_blastp_analysis(
    *,
    collected_loci: tuple[_CollectedLocus, ...],
    all_protein_alleles: tuple[_ProteinAllele, ...],
    working_directory_path: Path,
    blastp_executable: str,
    makeblastdb_executable: str,
    blast_evalue: float,
    thread_count: int,
) -> tuple[
    tuple[LocusRepresentativeBsrMeasurements, ...],
    tuple[InterLocusBlastScoreRatio, ...],
]:
    """Build one protein database and run representatives as batched queries."""
    representative_identifiers = {
        locus.representative_internal_identifier
        for locus in collected_loci
        if locus.representative_internal_identifier is not None
    }
    representative_alleles = tuple(
        protein_allele
        for protein_allele in all_protein_alleles
        if protein_allele.internal_identifier in representative_identifiers
    )
    if not representative_alleles or not all_protein_alleles:
        return (
            tuple(
                LocusRepresentativeBsrMeasurements(
                    locus_name=locus.locus_name,
                    representative_allele_identifier=(
                        locus.representative_allele_identifier
                    ),
                    representative_self_score=None,
                    total_allele_count=locus.total_allele_count,
                    valid_translated_allele_count=len(locus.protein_alleles),
                    allele_measurements=tuple(
                        AlleleBlastScoreRatio(
                            allele_identifier=protein_allele.original_identifier,
                            blast_score_ratio=None,
                        )
                        for protein_allele in locus.protein_alleles
                    ),
                )
                for locus in collected_loci
            ),
            (),
        )

    all_proteins_file_path = working_directory_path / "all_valid_proteins.fasta"
    representatives_file_path = working_directory_path / "representatives.fasta"
    database_path = working_directory_path / "all_valid_proteins_database"
    blast_output_file_path = working_directory_path / "representative_hits.tsv"
    _write_protein_fasta(all_proteins_file_path, all_protein_alleles)
    _write_protein_fasta(representatives_file_path, representative_alleles)
    create_blast_database(
        input_fasta_file_path=all_proteins_file_path,
        database_path=database_path,
        database_type="prot",
        makeblastdb_executable=makeblastdb_executable,
        operation_name="makeblastdb (representative BSR)",
    )
    _run_command(
        [
            blastp_executable,
            "-query",
            str(representatives_file_path),
            "-db",
            str(database_path),
            "-out",
            str(blast_output_file_path),
            "-outfmt",
            "6 qseqid sseqid score",
            "-task",
            "blastp",
            "-evalue",
            str(blast_evalue),
            "-num_threads",
            str(thread_count),
            "-max_target_seqs",
            str(len(all_protein_alleles)),
            "-max_hsps",
            "1",
        ],
        "BLASTp (representative BSR)",
    )
    return _create_bsr_measurements(
        collected_loci=collected_loci,
        all_protein_alleles=all_protein_alleles,
        blast_output_file_path=blast_output_file_path,
    )


def _inter_locus_measurement_sort_key(
    measurement: InterLocusBlastScoreRatio,
) -> tuple[float, float, str, str, str, str]:
    """Return a deterministic strongest-evidence ordering key."""
    return (
        measurement.blast_score_ratio,
        measurement.raw_score,
        measurement.query_locus_name,
        measurement.subject_locus_name,
        measurement.representative_allele_identifier,
        measurement.subject_allele_identifier,
    )


def _create_high_bsr_groups(
    high_bsr_pairs: tuple[HighBsrLocusPair, ...],
) -> tuple[HighBsrLocusGroup, ...]:
    """Create stable connected components from qualifying unordered pairs."""
    adjacency: dict[str, set[str]] = {}
    for pair in high_bsr_pairs:
        adjacency.setdefault(pair.first_locus_name, set()).add(
            pair.second_locus_name
        )
        adjacency.setdefault(pair.second_locus_name, set()).add(
            pair.first_locus_name
        )
    components: list[tuple[str, ...]] = []
    unvisited_loci = set(adjacency)
    while unvisited_loci:
        starting_locus = min(unvisited_loci)
        component_loci: set[str] = set()
        pending_loci = [starting_locus]
        while pending_loci:
            current_locus = pending_loci.pop()
            if current_locus in component_loci:
                continue
            component_loci.add(current_locus)
            unvisited_loci.discard(current_locus)
            pending_loci.extend(
                sorted(adjacency.get(current_locus, ()), reverse=True)
            )
        components.append(tuple(sorted(component_loci)))
    components.sort()

    high_bsr_groups: list[HighBsrLocusGroup] = []
    for group_number, locus_names in enumerate(components, start=1):
        locus_name_set = set(locus_names)
        component_pairs = tuple(
            pair
            for pair in high_bsr_pairs
            if pair.first_locus_name in locus_name_set
            and pair.second_locus_name in locus_name_set
        )
        high_bsr_groups.append(
            HighBsrLocusGroup(
                group_number=group_number,
                locus_names=locus_names,
                qualifying_pair_count=len(component_pairs),
                maximum_blast_score_ratio=max(
                    pair.evidence.blast_score_ratio for pair in component_pairs
                ),
            )
        )
    return tuple(high_bsr_groups)


def classify_representative_bsr_analysis(
    analysis: SchemeRepresentativeBsrAnalysis,
    *,
    blast_score_ratio_threshold: float = (
        DEFAULT_SAME_GENE_BLAST_SCORE_RATIO_THRESHOLD
    ),
) -> RepresentativeBsrClassification:
    """Classify cached BSR measurements at one advisory threshold.

    The comparison is inclusive (``BSR >= threshold``).  Threshold changes do
    not invoke BLAST or read the source FASTA files.
    """
    if (
        not isfinite(blast_score_ratio_threshold)
        or blast_score_ratio_threshold < 0
    ):
        raise ValueError("blast_score_ratio_threshold must be at least 0")

    locus_coverages: list[LocusRepresentativeBsrCoverage] = []
    for measurements in analysis.locus_measurements:
        has_usable_representative = (
            measurements.representative_allele_identifier is not None
            and measurements.representative_self_score is not None
            and measurements.representative_self_score > 0
        )
        covered_allele_count = sum(
            allele_measurement.blast_score_ratio is not None
            and allele_measurement.blast_score_ratio
            >= blast_score_ratio_threshold
            for allele_measurement in measurements.allele_measurements
        )
        below_threshold_hit_count = sum(
            allele_measurement.blast_score_ratio is not None
            and allele_measurement.blast_score_ratio
            < blast_score_ratio_threshold
            for allele_measurement in measurements.allele_measurements
        )
        no_reported_hit_count = sum(
            allele_measurement.blast_score_ratio is None
            for allele_measurement in measurements.allele_measurements
        )
        coverage_percentage = (
            covered_allele_count
            / measurements.valid_translated_allele_count
            * 100
            if has_usable_representative
            and measurements.valid_translated_allele_count > 0
            else None
        )
        locus_coverages.append(
            LocusRepresentativeBsrCoverage(
                locus_name=measurements.locus_name,
                representative_allele_identifier=(
                    measurements.representative_allele_identifier
                ),
                total_allele_count=measurements.total_allele_count,
                valid_translated_allele_count=(
                    measurements.valid_translated_allele_count
                ),
                unevaluable_allele_count=measurements.unevaluable_allele_count,
                covered_allele_count=covered_allele_count,
                below_threshold_hit_count=below_threshold_hit_count,
                no_reported_hit_count=no_reported_hit_count,
                coverage_percentage=coverage_percentage,
            )
        )

    best_measurements_by_unordered_pair: dict[
        tuple[str, str], InterLocusBlastScoreRatio
    ] = {}
    for measurement in analysis.inter_locus_measurements:
        if measurement.blast_score_ratio < blast_score_ratio_threshold:
            continue
        unordered_pair = tuple(
            sorted((measurement.query_locus_name, measurement.subject_locus_name))
        )
        current_measurement = best_measurements_by_unordered_pair.get(
            unordered_pair
        )
        if current_measurement is None or _inter_locus_measurement_sort_key(
            measurement
        ) > _inter_locus_measurement_sort_key(current_measurement):
            best_measurements_by_unordered_pair[unordered_pair] = measurement
    high_bsr_pairs = tuple(
        HighBsrLocusPair(
            first_locus_name=unordered_pair[0],
            second_locus_name=unordered_pair[1],
            evidence=best_measurements_by_unordered_pair[unordered_pair],
        )
        for unordered_pair in sorted(best_measurements_by_unordered_pair)
    )
    return RepresentativeBsrClassification(
        blast_score_ratio_threshold=blast_score_ratio_threshold,
        locus_coverages=tuple(locus_coverages),
        high_bsr_pairs=high_bsr_pairs,
        high_bsr_groups=_create_high_bsr_groups(high_bsr_pairs),
    )


def _deserialize_locus_measurements(
    data: dict[str, object],
) -> LocusRepresentativeBsrMeasurements:
    """Deserialize one cached per-locus measurement object."""
    allele_measurement_data = data["allele_measurements"]
    if not isinstance(allele_measurement_data, list):
        raise TypeError("allele_measurements must be a list")
    return LocusRepresentativeBsrMeasurements(
        locus_name=str(data["locus_name"]),
        representative_allele_identifier=(
            str(data["representative_allele_identifier"])
            if data.get("representative_allele_identifier") is not None
            else None
        ),
        representative_self_score=(
            float(data["representative_self_score"])
            if data.get("representative_self_score") is not None
            else None
        ),
        total_allele_count=int(data["total_allele_count"]),
        valid_translated_allele_count=int(data["valid_translated_allele_count"]),
        allele_measurements=tuple(
            AlleleBlastScoreRatio(
                allele_identifier=str(allele_data["allele_identifier"]),
                blast_score_ratio=(
                    float(allele_data["blast_score_ratio"])
                    if allele_data.get("blast_score_ratio") is not None
                    else None
                ),
            )
            for allele_data in allele_measurement_data
            if isinstance(allele_data, dict)
        ),
    )


def _deserialize_inter_locus_measurement(
    data: dict[str, object],
) -> InterLocusBlastScoreRatio:
    """Deserialize one cached directed inter-locus measurement."""
    return InterLocusBlastScoreRatio(
        query_locus_name=str(data["query_locus_name"]),
        subject_locus_name=str(data["subject_locus_name"]),
        representative_allele_identifier=str(
            data["representative_allele_identifier"]
        ),
        subject_allele_identifier=str(data["subject_allele_identifier"]),
        raw_score=float(data["raw_score"]),
        blast_score_ratio=float(data["blast_score_ratio"]),
    )


def _read_cached_analysis(
    *,
    scheme_directory_path: Path,
    source_signature: _SourceSignature,
    blastp_version: str | None,
    blast_evalue: float,
) -> SchemeRepresentativeBsrAnalysis | None:
    """Return a compatible cached analysis, or ``None`` when unavailable."""
    cache_file_path = next(
        (
            scheme_directory_path / cache_filename
            for cache_filename in SCHEME_REPRESENTATIVE_BSR_CACHE_FILENAMES
            if (scheme_directory_path / cache_filename).is_file()
        ),
        scheme_directory_path / SCHEME_REPRESENTATIVE_BSR_CACHE_FILENAME,
    )
    try:
        cache_data = json.loads(cache_file_path.read_text(encoding="utf-8"))
        if (
            int(cache_data["cache_version"])
            != SCHEME_REPRESENTATIVE_BSR_CACHE_VERSION
        ):
            return None
        if (
            str(cache_data["method_identifier"])
            != REPRESENTATIVE_BSR_METHOD_IDENTIFIER
        ):
            return None
        cached_blastp_version = str(cache_data["blastp_version"])
        if blastp_version is not None and cached_blastp_version != blastp_version:
            return None
        if float(cache_data["blast_evalue"]) != blast_evalue:
            return None
        cached_source_signature = tuple(
            (str(path), int(file_size), int(modification_time))
            for path, file_size, modification_time in cache_data[
                "source_signature"
            ]
        )
        if cached_source_signature != source_signature:
            return None
        return SchemeRepresentativeBsrAnalysis(
            locus_measurements=tuple(
                _deserialize_locus_measurements(locus_data)
                for locus_data in cache_data["locus_measurements"]
            ),
            inter_locus_measurements=tuple(
                _deserialize_inter_locus_measurement(measurement_data)
                for measurement_data in cache_data["inter_locus_measurements"]
            ),
            blastp_version=cached_blastp_version,
            method_identifier=REPRESENTATIVE_BSR_METHOD_IDENTIFIER,
            blast_evalue=blast_evalue,
            reused_from_cache=True,
        )
    except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError):
        return None


def _write_cached_analysis(
    *,
    scheme_directory_path: Path,
    source_signature: _SourceSignature,
    analysis: SchemeRepresentativeBsrAnalysis,
) -> None:
    """Atomically write the raw BSR analysis when the scheme is writable."""
    cache_data = {
        "cache_version": SCHEME_REPRESENTATIVE_BSR_CACHE_VERSION,
        "method_identifier": analysis.method_identifier,
        "blastp_version": analysis.blastp_version,
        "blast_evalue": analysis.blast_evalue,
        "source_signature": source_signature,
        "locus_measurements": [
            asdict(measurement) for measurement in analysis.locus_measurements
        ],
        "inter_locus_measurements": [
            asdict(measurement)
            for measurement in analysis.inter_locus_measurements
        ],
    }
    try:
        write_json_atomically(
            scheme_directory_path / SCHEME_REPRESENTATIVE_BSR_CACHE_FILENAME,
            cache_data,
        )
    except OSError:
        pass


def _emit_progress(
    progress: RepresentativeBsrProgress,
    *,
    progress_callback: RepresentativeBsrProgressCallback | None,
    progress_output_stream: TextIO | None,
) -> None:
    """Send progress to an optional callback and terminal output stream."""
    if progress_callback is not None:
        progress_callback(progress)
    if progress_output_stream is not None:
        print(
            f"[Representative BSR {progress.completed_stage_count}/"
            f"{progress.total_stage_count}] {progress.message}",
            file=progress_output_stream,
            flush=True,
        )


def get_scheme_representative_bsr_analysis(
    scheme_directory_path: str | Path,
    *,
    maximum_threads: int | None = None,
    blast_evalue: float = DEFAULT_BLAST_EVALUE,
    blastp_executable: str = "blastp",
    makeblastdb_executable: str = "makeblastdb",
    progress_callback: RepresentativeBsrProgressCallback | None = None,
    progress_output_stream: TextIO | None = None,
) -> SchemeRepresentativeBsrAnalysis:
    """Measure representative-to-all-valid-allele protein BSR values.

    Args:
        scheme_directory_path: Extracted scheme directory containing locus
            FASTA files.
        maximum_threads: BLASTp thread budget. The default uses 80% of CPUs
            available to the process, rounded down to an even count when
            possible.
        blast_evalue: BLASTp reporting floor.  It controls search sensitivity,
            not BSR classification.
        blastp_executable: ``blastp`` executable name or path.
        makeblastdb_executable: ``makeblastdb`` executable name or path.
        progress_callback: Optional callback for notebook progress updates.
        progress_output_stream: Optional stream for terminal progress lines.

    Returns:
        Threshold-independent cached measurements.  Apply
        :func:`classify_representative_bsr_analysis` to select a BSR threshold.
    """
    if not isfinite(blast_evalue) or blast_evalue <= 0:
        raise ValueError("blast_evalue must be greater than zero")
    resolved_scheme_directory_path = _resolve_scheme_directory_path(
        scheme_directory_path
    )
    fasta_file_paths = _find_all_fasta_file_paths(
        resolved_scheme_directory_path
    )
    source_signature = _create_source_signature(
        resolved_scheme_directory_path,
        fasta_file_paths,
    )
    try:
        blastp_version = _get_blast_version(blastp_executable)
    except RuntimeError:
        cached_analysis_without_blast = _read_cached_analysis(
            scheme_directory_path=resolved_scheme_directory_path,
            source_signature=source_signature,
            blastp_version=None,
            blast_evalue=blast_evalue,
        )
        if cached_analysis_without_blast is None:
            raise
        _emit_progress(
            RepresentativeBsrProgress(
                completed_stage_count=REPRESENTATIVE_BSR_PROGRESS_STAGE_COUNT,
                total_stage_count=REPRESENTATIVE_BSR_PROGRESS_STAGE_COUNT,
                stage_name="cache",
                message=(
                    "Reused the complete representative BSR cache; BLASTp "
                    "was not needed."
                ),
                reused_from_cache=True,
            ),
            progress_callback=progress_callback,
            progress_output_stream=progress_output_stream,
        )
        return cached_analysis_without_blast

    cached_analysis = _read_cached_analysis(
        scheme_directory_path=resolved_scheme_directory_path,
        source_signature=source_signature,
        blastp_version=blastp_version,
        blast_evalue=blast_evalue,
    )
    if cached_analysis is not None:
        _emit_progress(
            RepresentativeBsrProgress(
                completed_stage_count=REPRESENTATIVE_BSR_PROGRESS_STAGE_COUNT,
                total_stage_count=REPRESENTATIVE_BSR_PROGRESS_STAGE_COUNT,
                stage_name="cache",
                message="Reused the complete representative BSR cache.",
                reused_from_cache=True,
            ),
            progress_callback=progress_callback,
            progress_output_stream=progress_output_stream,
        )
        return cached_analysis

    collected_loci, all_protein_alleles = _collect_protein_alleles(
        resolved_scheme_directory_path,
        fasta_file_paths,
    )
    representative_count = sum(
        locus.representative_internal_identifier is not None
        for locus in collected_loci
    )
    _emit_progress(
        RepresentativeBsrProgress(
            completed_stage_count=1,
            total_stage_count=REPRESENTATIVE_BSR_PROGRESS_STAGE_COUNT,
            stage_name="sequences",
            message=(
                f"Selected {representative_count:,} approximate representatives "
                f"and translated {len(all_protein_alleles):,} valid alleles "
                f"across {len(collected_loci):,} loci."
            ),
        ),
        progress_callback=progress_callback,
        progress_output_stream=progress_output_stream,
    )
    _emit_progress(
        RepresentativeBsrProgress(
            completed_stage_count=2,
            total_stage_count=REPRESENTATIVE_BSR_PROGRESS_STAGE_COUNT,
            stage_name="database",
            message="Preparing the combined valid-protein BLAST database.",
        ),
        progress_callback=progress_callback,
        progress_output_stream=progress_output_stream,
    )
    with tempfile.TemporaryDirectory(prefix="scheme-representative-bsr-") as directory:
        locus_measurements, inter_locus_measurements = _run_blastp_analysis(
            collected_loci=collected_loci,
            all_protein_alleles=all_protein_alleles,
            working_directory_path=Path(directory),
            blastp_executable=blastp_executable,
            makeblastdb_executable=makeblastdb_executable,
            blast_evalue=blast_evalue,
            thread_count=resolve_worker_count(maximum_threads),
        )
    _emit_progress(
        RepresentativeBsrProgress(
            completed_stage_count=3,
            total_stage_count=REPRESENTATIVE_BSR_PROGRESS_STAGE_COUNT,
            stage_name="blastp",
            message="Completed representative-versus-allele BLASTp.",
        ),
        progress_callback=progress_callback,
        progress_output_stream=progress_output_stream,
    )
    _emit_progress(
        RepresentativeBsrProgress(
            completed_stage_count=4,
            total_stage_count=REPRESENTATIVE_BSR_PROGRESS_STAGE_COUNT,
            stage_name="measurements",
            message=(
                f"Reduced BLASTp output to {len(inter_locus_measurements):,} "
                "directed inter-locus maximum-BSR measurements."
            ),
        ),
        progress_callback=progress_callback,
        progress_output_stream=progress_output_stream,
    )
    analysis = SchemeRepresentativeBsrAnalysis(
        locus_measurements=locus_measurements,
        inter_locus_measurements=inter_locus_measurements,
        blastp_version=blastp_version,
        method_identifier=REPRESENTATIVE_BSR_METHOD_IDENTIFIER,
        blast_evalue=blast_evalue,
        reused_from_cache=False,
    )
    _write_cached_analysis(
        scheme_directory_path=resolved_scheme_directory_path,
        source_signature=source_signature,
        analysis=analysis,
    )
    _emit_progress(
        RepresentativeBsrProgress(
            completed_stage_count=REPRESENTATIVE_BSR_PROGRESS_STAGE_COUNT,
            total_stage_count=REPRESENTATIVE_BSR_PROGRESS_STAGE_COUNT,
            stage_name="cache",
            message="Representative BSR analysis complete and cached.",
        ),
        progress_callback=progress_callback,
        progress_output_stream=progress_output_stream,
    )
    return replace(analysis, reused_from_cache=False)


def _parse_non_negative_float(argument_value: str) -> float:
    """Parse a finite, non-negative command-line floating-point value."""
    parsed_value = float(argument_value)
    if not isfinite(parsed_value) or parsed_value < 0:
        raise argparse.ArgumentTypeError("value must be a finite number at least 0")
    return parsed_value


def _parse_positive_integer(argument_value: str) -> int:
    """Parse a positive command-line integer."""
    parsed_value = int(argument_value)
    if parsed_value < 1:
        raise argparse.ArgumentTypeError("value must be at least 1")
    return parsed_value


def main() -> None:
    """Run representative BSR coverage analysis from the command line."""
    argument_parser = argparse.ArgumentParser(description=__doc__)
    argument_parser.add_argument(
        "scheme_directory",
        type=Path,
        help="Extracted scheme directory to analyze.",
    )
    argument_parser.add_argument(
        "--bsr",
        type=_parse_non_negative_float,
        default=DEFAULT_SAME_GENE_BLAST_SCORE_RATIO_THRESHOLD,
        help="Inclusive BSR reference threshold (default: 0.6).",
    )
    argument_parser.add_argument(
        "--threads",
        type=_parse_positive_integer,
        default=None,
        help=(
            "BLASTp thread budget (default: 80%% of available CPUs, rounded "
            "down to an even count when possible)."
        ),
    )
    argument_parser.add_argument(
        "--json",
        action="store_true",
        help="Emit the complete raw analysis and classified result as JSON.",
    )
    parsed_arguments = argument_parser.parse_args()
    analysis = get_scheme_representative_bsr_analysis(
        parsed_arguments.scheme_directory,
        maximum_threads=parsed_arguments.threads,
        progress_output_stream=sys.stderr,
    )
    classification = classify_representative_bsr_analysis(
        analysis,
        blast_score_ratio_threshold=parsed_arguments.bsr,
    )
    if parsed_arguments.json:
        print(
            json.dumps(
                {
                    "analysis": asdict(analysis),
                    "classification": asdict(classification),
                },
                indent=2,
            )
        )
        return
    covered_percentage = classification.covered_allele_percentage
    print(f"Loci: {analysis.locus_count}")
    print(f"Loci with representatives: {analysis.representative_locus_count}")
    print(f"Valid translated alleles: {analysis.valid_translated_allele_count}")
    print(
        "Representative-covered alleles: "
        f"{classification.covered_allele_count}/"
        f"{classification.measurable_allele_count}"
        + (
            f" ({covered_percentage:.2f}%)"
            if covered_percentage is not None
            else " (unavailable)"
        )
    )
    print(f"High-BSR inter-locus pairs: {len(classification.high_bsr_pairs)}")
    print(f"High-BSR groups: {len(classification.high_bsr_groups)}")
    print(f"Cache reused: {analysis.reused_from_cache}")


if __name__ == "__main__":
    main()
