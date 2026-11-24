"""Inspect an extracted MLST scheme and report basic sequence statistics.

This module counts loci and evaluates allele coding sequences in a directory
produced by ``scheme.cache.cache_archive``. A locus is represented by one FASTA
file, and an allele is represented by one FASTA record within those files. The
evaluation reports valid alleles, individual invalidity categories, and alleles
shorter than a configurable minimum sequence length. It also reports
scheme-wide counts and percentages of alleles below 80% or above 120% of their
respective locus's modal length, together with pooled scheme-wide GC content.
GC content uses only canonical A/C/G/T bases as its denominator; ambiguous and
other noncanonical characters are counted separately. Common FASTA filename
extensions are recognized in both uncompressed and gzip-compressed form.

Notebook code can call :func:`count_loci`, :func:`count_alleles`, or
:func:`get_scheme_information`. The module can also be executed as a
command-line program and can print either a readable summary or JSON.
"""

import argparse
import gzip
import json
from collections import Counter
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path

from .cache_io import write_json_atomically
from .fasta import (
    FASTA_FILENAME_SUFFIXES,
    find_locus_fasta_paths,
    is_fasta_file,
    iterate_fasta_records,
    resolve_scheme_directory_path,
)
from .sequence import (
    BACTERIAL_START_CODONS,
    STOP_CODONS,
    VALID_DNA_BASES,
    NucleotideBaseComposition,
    contains_in_frame_stop_codon,
    count_nucleotide_base_composition,
    evaluate_allele_sequence,
)

SCHEME_INFORMATION_CACHE_FILENAME = ".scheme-information-cache.json"
LEGACY_INFORMATION_CACHE_FILENAME = ".schema-information-cache.json"
SCHEME_INFORMATION_CACHE_FILENAMES = (
    SCHEME_INFORMATION_CACHE_FILENAME,
    LEGACY_INFORMATION_CACHE_FILENAME,
)
SCHEME_INFORMATION_CACHE_VERSION = 4

__all__ = [
    "BACTERIAL_START_CODONS",
    "FASTA_FILENAME_SUFFIXES",
    "SCHEME_INFORMATION_CACHE_FILENAME",
    "SCHEME_INFORMATION_CACHE_FILENAMES",
    "SCHEME_INFORMATION_CACHE_VERSION",
    "STOP_CODONS",
    "VALID_DNA_BASES",
    "NucleotideBaseComposition",
    "SchemeInformation",
    "count_alleles",
    "count_alleles_below_minimum_length",
    "count_invalid_alleles",
    "count_loci",
    "count_nucleotide_base_composition",
    "count_valid_alleles",
    "get_scheme_information",
    "is_valid_allele_sequence",
]


@dataclass(frozen=True)
class SchemeInformation:
    """Summary counts for an extracted MLST scheme.

    Attributes:
        locus_count: Number of recognized FASTA files in the scheme.
        allele_count: Number of FASTA records across all recognized files.
        valid_allele_count: Number of alleles that satisfy every coding-sequence
            validation rule.
        invalid_allele_count: Sum of the incomplete ORF, ambiguous bases,
            missing start/stop codon, and in-frame stop codon counts. An allele
            can contribute to more than one category.
        incomplete_orf_allele_count: Number of alleles whose length is not a
            multiple of three.
        ambiguous_bases_allele_count: Number of alleles containing characters
            other than ``A``, ``C``, ``G``, and ``T``.
        missing_start_or_stop_codon_allele_count: Number of alleles without a
            recognized start codon, terminal stop codon, or both.
        in_frame_stop_codon_allele_count: Number of alleles containing a stop
            codon before the end of the coding sequence.
        alleles_below_minimum_length_count: Number of alleles shorter than
            ``minimum_sequence_length``.
        minimum_sequence_length: Minimum allele length used for the short-allele
            count, in nucleotides.
        alleles_below_locus_mode_threshold_count: Total alleles across the
            scheme strictly shorter than 80% of their locus's modal length.
        alleles_above_locus_mode_threshold_count: Total alleles across the
            scheme strictly longer than 120% of their locus's modal length.
        alleles_below_locus_mode_threshold_percentage: Below-threshold count as
            a percentage of all alleles in the scheme.
        alleles_above_locus_mode_threshold_percentage: Above-threshold count as
            a percentage of all alleles in the scheme.
        guanine_cytosine_content_percentage: Pooled scheme-wide GC percentage
            using canonical bases only, or ``None`` when the scheme contains no
            canonical bases.
        canonical_nucleotide_base_count: Total A/C/G/T bases used as the
            scheme-wide GC-content denominator.
        excluded_nucleotide_base_count: Total noncanonical, non-whitespace
            characters excluded from the GC-content denominator.
    """

    locus_count: int
    allele_count: int
    valid_allele_count: int
    invalid_allele_count: int
    incomplete_orf_allele_count: int
    ambiguous_bases_allele_count: int
    missing_start_or_stop_codon_allele_count: int
    in_frame_stop_codon_allele_count: int
    alleles_below_minimum_length_count: int
    minimum_sequence_length: int
    alleles_below_locus_mode_threshold_count: int
    alleles_above_locus_mode_threshold_count: int
    alleles_below_locus_mode_threshold_percentage: float
    alleles_above_locus_mode_threshold_percentage: float
    guanine_cytosine_content_percentage: float | None
    canonical_nucleotide_base_count: int
    excluded_nucleotide_base_count: int


@dataclass(frozen=True)
class _AlleleValidationResult:
    """Validation flags for one normalized allele sequence."""

    has_incomplete_orf: bool
    has_ambiguous_bases: bool
    is_missing_start_or_stop_codon: bool
    has_in_frame_stop_codon: bool
    is_below_minimum_length: bool

    @property
    def is_valid(self) -> bool:
        """Return whether no coding-sequence invalidity flag is set."""
        return not any(
            (
                self.has_incomplete_orf,
                self.has_ambiguous_bases,
                self.is_missing_start_or_stop_codon,
                self.has_in_frame_stop_codon,
            )
        )


@dataclass(frozen=True)
class _AlleleCounts:
    """Aggregated allele validation counts from parsed FASTA records."""

    allele_count: int
    valid_allele_count: int
    incomplete_orf_allele_count: int
    ambiguous_bases_allele_count: int
    missing_start_or_stop_codon_allele_count: int
    in_frame_stop_codon_allele_count: int
    alleles_below_minimum_length_count: int
    allele_lengths: tuple[int, ...]
    guanine_cytosine_base_count: int
    canonical_nucleotide_base_count: int
    excluded_nucleotide_base_count: int

    @property
    def invalid_allele_count(self) -> int:
        """Return the requested sum of all invalidity category counts."""
        return (
            self.incomplete_orf_allele_count
            + self.ambiguous_bases_allele_count
            + self.missing_start_or_stop_codon_allele_count
            + self.in_frame_stop_codon_allele_count
        )


@dataclass(frozen=True)
class _CachedSchemeAnalysis:
    """Scheme metrics and allele lengths persisted in the extracted scheme."""

    locus_count: int
    allele_counts: _AlleleCounts
    alleles_below_locus_mode_threshold_count: int
    alleles_above_locus_mode_threshold_count: int


_resolve_scheme_directory_path = resolve_scheme_directory_path
_is_fasta_file = is_fasta_file
_find_allele_file_paths = find_locus_fasta_paths


def _iterate_fasta_sequences(fasta_lines: Iterable[str]) -> Iterable[str]:
    """Yield one assembled nucleotide sequence per FASTA header."""
    return (record.sequence for record in iterate_fasta_records(fasta_lines))


_contains_in_frame_stop_codon = contains_in_frame_stop_codon


def _evaluate_allele_sequence(
    allele_sequence: str,
    minimum_sequence_length: int,
) -> _AlleleValidationResult:
    """Evaluate all coding-sequence and minimum-length rules for one allele."""
    normalized_allele_sequence = "".join(allele_sequence.split()).upper()
    shared_validation = evaluate_allele_sequence(normalized_allele_sequence)

    return _AlleleValidationResult(
        has_incomplete_orf=shared_validation.has_incomplete_orf,
        has_ambiguous_bases=shared_validation.has_ambiguous_bases,
        is_missing_start_or_stop_codon=(
            shared_validation.is_missing_start_or_stop_codon
        ),
        has_in_frame_stop_codon=shared_validation.has_in_frame_stop_codon,
        is_below_minimum_length=(
            len(normalized_allele_sequence) < minimum_sequence_length
        ),
    )


def is_valid_allele_sequence(allele_sequence: str) -> bool:
    """Return whether an allele is a valid complete bacterial coding sequence.

    A valid allele contains only ``A``, ``C``, ``G``, and ``T``; has a length
    divisible by three; begins with ``ATG``, ``GTG``, or ``TTG``; ends with
    ``TAA``, ``TAG``, or ``TGA``; and contains no in-frame stop codon before
    its terminal stop. Whitespace and nucleotide letter case are ignored.

    Args:
        allele_sequence: Nucleotide sequence to validate.

    Returns:
        ``True`` when every validation rule passes, otherwise ``False``.
    """
    allele_validation_result = _evaluate_allele_sequence(
        allele_sequence,
        minimum_sequence_length=0,
    )
    return allele_validation_result.is_valid


def _count_alleles_in_lines(
    fasta_lines: Iterable[str],
    minimum_sequence_length: int,
) -> _AlleleCounts:
    """Count allele validation outcomes from FASTA text lines."""
    allele_count = 0
    valid_allele_count = 0
    incomplete_orf_allele_count = 0
    ambiguous_bases_allele_count = 0
    missing_start_or_stop_codon_allele_count = 0
    in_frame_stop_codon_allele_count = 0
    alleles_below_minimum_length_count = 0
    allele_lengths: list[int] = []
    guanine_cytosine_base_count = 0
    canonical_nucleotide_base_count = 0
    excluded_nucleotide_base_count = 0

    for allele_sequence in _iterate_fasta_sequences(fasta_lines):
        normalized_allele_sequence = "".join(allele_sequence.split())
        allele_lengths.append(len(normalized_allele_sequence))
        allele_count += 1
        nucleotide_base_composition = count_nucleotide_base_composition(
            normalized_allele_sequence
        )
        guanine_cytosine_base_count += (
            nucleotide_base_composition.guanine_cytosine_base_count
        )
        canonical_nucleotide_base_count += (
            nucleotide_base_composition.canonical_nucleotide_base_count
        )
        excluded_nucleotide_base_count += (
            nucleotide_base_composition.excluded_nucleotide_base_count
        )
        allele_validation_result = _evaluate_allele_sequence(
            allele_sequence,
            minimum_sequence_length,
        )
        if allele_validation_result.is_valid:
            valid_allele_count += 1
        if allele_validation_result.has_incomplete_orf:
            incomplete_orf_allele_count += 1
        if allele_validation_result.has_ambiguous_bases:
            ambiguous_bases_allele_count += 1
        if allele_validation_result.is_missing_start_or_stop_codon:
            missing_start_or_stop_codon_allele_count += 1
        if allele_validation_result.has_in_frame_stop_codon:
            in_frame_stop_codon_allele_count += 1
        if allele_validation_result.is_below_minimum_length:
            alleles_below_minimum_length_count += 1

    return _AlleleCounts(
        allele_count=allele_count,
        valid_allele_count=valid_allele_count,
        incomplete_orf_allele_count=incomplete_orf_allele_count,
        ambiguous_bases_allele_count=ambiguous_bases_allele_count,
        missing_start_or_stop_codon_allele_count=(
            missing_start_or_stop_codon_allele_count
        ),
        in_frame_stop_codon_allele_count=in_frame_stop_codon_allele_count,
        alleles_below_minimum_length_count=alleles_below_minimum_length_count,
        allele_lengths=tuple(allele_lengths),
        guanine_cytosine_base_count=guanine_cytosine_base_count,
        canonical_nucleotide_base_count=canonical_nucleotide_base_count,
        excluded_nucleotide_base_count=excluded_nucleotide_base_count,
    )


def _count_alleles_in_file(
    allele_file_path: Path,
    minimum_sequence_length: int,
) -> _AlleleCounts:
    """Count total and valid alleles in one FASTA file."""
    if allele_file_path.name.lower().endswith(".gz"):
        with gzip.open(
            allele_file_path,
            mode="rt",
            encoding="utf-8",
        ) as compressed_allele_file:
            return _count_alleles_in_lines(
                compressed_allele_file,
                minimum_sequence_length,
            )

    with allele_file_path.open(encoding="utf-8") as allele_file:
        return _count_alleles_in_lines(
            allele_file,
            minimum_sequence_length,
        )


def _count_alleles_in_files(
    allele_file_paths: Iterable[Path],
    minimum_sequence_length: int,
) -> tuple[_AlleleCounts, int, int]:
    """Count alleles and aggregate modal-length outliers across loci.

    Each FASTA file is treated as one locus. Its modal length is used only
    while accumulating the two scheme-level outlier counts; no per-locus
    result is retained or exposed.

    Args:
        allele_file_paths: FASTA files to parse as loci.
        minimum_sequence_length: Minimum accepted allele length in
            nucleotides.

    Returns:
        The aggregate allele validation counts, followed by the scheme-wide
        counts below 80% and above 120% of the respective locus modes.
    """
    allele_count = 0
    valid_allele_count = 0
    incomplete_orf_allele_count = 0
    ambiguous_bases_allele_count = 0
    missing_start_or_stop_codon_allele_count = 0
    in_frame_stop_codon_allele_count = 0
    alleles_below_minimum_length_count = 0
    allele_lengths: list[int] = []
    alleles_below_locus_mode_threshold_count = 0
    alleles_above_locus_mode_threshold_count = 0
    guanine_cytosine_base_count = 0
    canonical_nucleotide_base_count = 0
    excluded_nucleotide_base_count = 0

    for allele_file_path in allele_file_paths:
        allele_file_counts = _count_alleles_in_file(
            allele_file_path,
            minimum_sequence_length,
        )
        allele_count += allele_file_counts.allele_count
        valid_allele_count += allele_file_counts.valid_allele_count
        incomplete_orf_allele_count += (
            allele_file_counts.incomplete_orf_allele_count
        )
        ambiguous_bases_allele_count += (
            allele_file_counts.ambiguous_bases_allele_count
        )
        missing_start_or_stop_codon_allele_count += (
            allele_file_counts.missing_start_or_stop_codon_allele_count
        )
        in_frame_stop_codon_allele_count += (
            allele_file_counts.in_frame_stop_codon_allele_count
        )
        alleles_below_minimum_length_count += (
            allele_file_counts.alleles_below_minimum_length_count
        )
        allele_lengths.extend(allele_file_counts.allele_lengths)
        guanine_cytosine_base_count += (
            allele_file_counts.guanine_cytosine_base_count
        )
        canonical_nucleotide_base_count += (
            allele_file_counts.canonical_nucleotide_base_count
        )
        excluded_nucleotide_base_count += (
            allele_file_counts.excluded_nucleotide_base_count
        )
        (
            locus_alleles_below_threshold_count,
            locus_alleles_above_threshold_count,
        ) = _count_locus_length_outliers(allele_file_counts.allele_lengths)
        alleles_below_locus_mode_threshold_count += (
            locus_alleles_below_threshold_count
        )
        alleles_above_locus_mode_threshold_count += (
            locus_alleles_above_threshold_count
        )

    allele_counts = _AlleleCounts(
        allele_count=allele_count,
        valid_allele_count=valid_allele_count,
        incomplete_orf_allele_count=incomplete_orf_allele_count,
        ambiguous_bases_allele_count=ambiguous_bases_allele_count,
        missing_start_or_stop_codon_allele_count=(
            missing_start_or_stop_codon_allele_count
        ),
        in_frame_stop_codon_allele_count=in_frame_stop_codon_allele_count,
        alleles_below_minimum_length_count=alleles_below_minimum_length_count,
        allele_lengths=tuple(allele_lengths),
        guanine_cytosine_base_count=guanine_cytosine_base_count,
        canonical_nucleotide_base_count=canonical_nucleotide_base_count,
        excluded_nucleotide_base_count=excluded_nucleotide_base_count,
    )
    return (
        allele_counts,
        alleles_below_locus_mode_threshold_count,
        alleles_above_locus_mode_threshold_count,
    )


def _count_locus_length_outliers(
    allele_lengths: tuple[int, ...],
) -> tuple[int, int]:
    """Count length outliers for one locus for use in scheme-wide totals.

    When multiple lengths have the same highest frequency, the shortest tied
    length is selected as the mode. Empty loci contribute zero to both totals.

    Args:
        allele_lengths: Sequence lengths for all alleles in one locus.

    Returns:
        A pair containing the counts strictly below 80% and strictly above
        120% of the selected modal length.
    """
    if not allele_lengths:
        return 0, 0

    length_frequencies = Counter(allele_lengths)
    highest_length_frequency = max(length_frequencies.values())
    modal_allele_length_candidates = sorted(
        allele_length
        for allele_length, allele_length_frequency in length_frequencies.items()
        if allele_length_frequency == highest_length_frequency
    )
    modal_allele_length = modal_allele_length_candidates[0]
    lower_length_threshold = modal_allele_length * 0.8
    upper_length_threshold = modal_allele_length * 1.2
    alleles_below_threshold_count = sum(
        allele_length < lower_length_threshold for allele_length in allele_lengths
    )
    alleles_above_threshold_count = sum(
        allele_length > upper_length_threshold for allele_length in allele_lengths
    )

    return alleles_below_threshold_count, alleles_above_threshold_count


def _read_cached_scheme_analysis(
    scheme_directory_path: Path,
) -> _CachedSchemeAnalysis | None:
    """Read a compatible cached analysis, returning ``None`` if unavailable."""
    scheme_information_cache_file_path = next(
        (
            scheme_directory_path / cache_filename
            for cache_filename in SCHEME_INFORMATION_CACHE_FILENAMES
            if (scheme_directory_path / cache_filename).is_file()
        ),
        scheme_directory_path / SCHEME_INFORMATION_CACHE_FILENAME,
    )

    try:
        cached_analysis_data = json.loads(
            scheme_information_cache_file_path.read_text(encoding="utf-8")
        )
        if (
            cached_analysis_data["cache_version"]
            != SCHEME_INFORMATION_CACHE_VERSION
        ):
            return None

        allele_lengths = tuple(
            int(allele_length)
            for allele_length in cached_analysis_data["allele_lengths"]
        )
        alleles_below_locus_mode_threshold_count = int(
            cached_analysis_data[
                "alleles_below_locus_mode_threshold_count"
            ]
        )
        alleles_above_locus_mode_threshold_count = int(
            cached_analysis_data[
                "alleles_above_locus_mode_threshold_count"
            ]
        )
        allele_counts = _AlleleCounts(
            allele_count=int(cached_analysis_data["allele_count"]),
            valid_allele_count=int(cached_analysis_data["valid_allele_count"]),
            incomplete_orf_allele_count=int(
                cached_analysis_data["incomplete_orf_allele_count"]
            ),
            ambiguous_bases_allele_count=int(
                cached_analysis_data["ambiguous_bases_allele_count"]
            ),
            missing_start_or_stop_codon_allele_count=int(
                cached_analysis_data[
                    "missing_start_or_stop_codon_allele_count"
                ]
            ),
            in_frame_stop_codon_allele_count=int(
                cached_analysis_data["in_frame_stop_codon_allele_count"]
            ),
            alleles_below_minimum_length_count=0,
            allele_lengths=allele_lengths,
            guanine_cytosine_base_count=int(
                cached_analysis_data["guanine_cytosine_base_count"]
            ),
            canonical_nucleotide_base_count=int(
                cached_analysis_data["canonical_nucleotide_base_count"]
            ),
            excluded_nucleotide_base_count=int(
                cached_analysis_data["excluded_nucleotide_base_count"]
            ),
        )
        locus_count = int(cached_analysis_data["locus_count"])
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return None

    if allele_counts.allele_count != len(allele_lengths):
        return None
    if alleles_below_locus_mode_threshold_count < 0:
        return None
    if alleles_above_locus_mode_threshold_count < 0:
        return None
    if (
        alleles_below_locus_mode_threshold_count
        + alleles_above_locus_mode_threshold_count
        > allele_counts.allele_count
    ):
        return None
    if min(
        allele_counts.guanine_cytosine_base_count,
        allele_counts.canonical_nucleotide_base_count,
        allele_counts.excluded_nucleotide_base_count,
    ) < 0:
        return None
    if (
        allele_counts.guanine_cytosine_base_count
        > allele_counts.canonical_nucleotide_base_count
    ):
        return None

    return _CachedSchemeAnalysis(
        locus_count=locus_count,
        allele_counts=allele_counts,
        alleles_below_locus_mode_threshold_count=(
            alleles_below_locus_mode_threshold_count
        ),
        alleles_above_locus_mode_threshold_count=(
            alleles_above_locus_mode_threshold_count
        ),
    )


def _write_cached_scheme_analysis(
    scheme_directory_path: Path,
    cached_scheme_analysis: _CachedSchemeAnalysis,
) -> None:
    """Atomically persist scheme metrics and allele lengths as JSON."""
    allele_counts = cached_scheme_analysis.allele_counts
    cached_analysis_data = {
        "cache_version": SCHEME_INFORMATION_CACHE_VERSION,
        "locus_count": cached_scheme_analysis.locus_count,
        "allele_count": allele_counts.allele_count,
        "valid_allele_count": allele_counts.valid_allele_count,
        "incomplete_orf_allele_count": (
            allele_counts.incomplete_orf_allele_count
        ),
        "ambiguous_bases_allele_count": (
            allele_counts.ambiguous_bases_allele_count
        ),
        "missing_start_or_stop_codon_allele_count": (
            allele_counts.missing_start_or_stop_codon_allele_count
        ),
        "in_frame_stop_codon_allele_count": (
            allele_counts.in_frame_stop_codon_allele_count
        ),
        "allele_lengths": allele_counts.allele_lengths,
        "guanine_cytosine_base_count": (
            allele_counts.guanine_cytosine_base_count
        ),
        "canonical_nucleotide_base_count": (
            allele_counts.canonical_nucleotide_base_count
        ),
        "excluded_nucleotide_base_count": (
            allele_counts.excluded_nucleotide_base_count
        ),
        "alleles_below_locus_mode_threshold_count": (
            cached_scheme_analysis.alleles_below_locus_mode_threshold_count
        ),
        "alleles_above_locus_mode_threshold_count": (
            cached_scheme_analysis.alleles_above_locus_mode_threshold_count
        ),
    }
    scheme_information_cache_file_path = (
        scheme_directory_path / SCHEME_INFORMATION_CACHE_FILENAME
    )

    write_json_atomically(
        scheme_information_cache_file_path,
        cached_analysis_data,
    )


def _get_or_create_cached_scheme_analysis(
    scheme_directory_path: Path,
) -> _CachedSchemeAnalysis:
    """Return cached scheme metrics, calculating and persisting them once."""
    cached_scheme_analysis = _read_cached_scheme_analysis(scheme_directory_path)
    if cached_scheme_analysis is not None:
        return cached_scheme_analysis

    from .static_analysis import get_shared_scheme_static_analysis

    shared_static_analysis = get_shared_scheme_static_analysis(
        scheme_directory_path
    )
    allele_lengths = tuple(
        allele_length
        for locus_analysis in shared_static_analysis.loci
        for allele_length in locus_analysis.allele_lengths
    )
    allele_counts = _AlleleCounts(
        allele_count=len(allele_lengths),
        valid_allele_count=sum(
            locus_analysis.valid_allele_count
            for locus_analysis in shared_static_analysis.loci
        ),
        incomplete_orf_allele_count=sum(
            locus_analysis.incomplete_orf_allele_count
            for locus_analysis in shared_static_analysis.loci
        ),
        ambiguous_bases_allele_count=sum(
            locus_analysis.ambiguous_bases_allele_count
            for locus_analysis in shared_static_analysis.loci
        ),
        missing_start_or_stop_codon_allele_count=sum(
            locus_analysis.missing_start_or_stop_codon_allele_count
            for locus_analysis in shared_static_analysis.loci
        ),
        in_frame_stop_codon_allele_count=sum(
            locus_analysis.in_frame_stop_codon_allele_count
            for locus_analysis in shared_static_analysis.loci
        ),
        alleles_below_minimum_length_count=0,
        allele_lengths=allele_lengths,
        guanine_cytosine_base_count=sum(
            locus_analysis.guanine_cytosine_base_count
            for locus_analysis in shared_static_analysis.loci
        ),
        canonical_nucleotide_base_count=sum(
            locus_analysis.canonical_nucleotide_base_count
            for locus_analysis in shared_static_analysis.loci
        ),
        excluded_nucleotide_base_count=sum(
            locus_analysis.excluded_nucleotide_base_count
            for locus_analysis in shared_static_analysis.loci
        ),
    )

    cached_scheme_analysis = _CachedSchemeAnalysis(
        locus_count=len(shared_static_analysis.loci),
        allele_counts=allele_counts,
        alleles_below_locus_mode_threshold_count=(
            sum(
                locus_analysis.alleles_below_mode_threshold_count
                for locus_analysis in shared_static_analysis.loci
            )
        ),
        alleles_above_locus_mode_threshold_count=(
            sum(
                locus_analysis.alleles_above_mode_threshold_count
                for locus_analysis in shared_static_analysis.loci
            )
        ),
    )

    try:
        _write_cached_scheme_analysis(
            scheme_directory_path,
            cached_scheme_analysis,
        )
    except OSError:
        # Read-only scheme directories can still be analyzed without caching.
        pass

    return cached_scheme_analysis


def count_loci(scheme_directory_path: str | Path) -> int:
    """Count locus FASTA files in an extracted scheme directory.

    Args:
        scheme_directory_path: Directory containing an extracted MLST scheme.

    Returns:
        The number of recognized FASTA files found recursively.

    Raises:
        FileNotFoundError: If ``scheme_directory_path`` does not exist.
        NotADirectoryError: If ``scheme_directory_path`` is not a directory.
    """
    resolved_scheme_directory_path = _resolve_scheme_directory_path(
        scheme_directory_path
    )
    cached_scheme_analysis = _get_or_create_cached_scheme_analysis(
        resolved_scheme_directory_path
    )
    return cached_scheme_analysis.locus_count


def count_alleles(scheme_directory_path: str | Path) -> int:
    """Count FASTA records across all loci in an extracted scheme directory.

    Args:
        scheme_directory_path: Directory containing an extracted MLST scheme.

    Returns:
        The total number of FASTA records across all recognized locus files.

    Raises:
        FileNotFoundError: If ``scheme_directory_path`` does not exist.
        NotADirectoryError: If ``scheme_directory_path`` is not a directory.
        OSError: If an allele file cannot be read or decompressed.
        UnicodeError: If an allele file is not valid UTF-8 text.
    """
    resolved_scheme_directory_path = _resolve_scheme_directory_path(
        scheme_directory_path
    )
    cached_scheme_analysis = _get_or_create_cached_scheme_analysis(
        resolved_scheme_directory_path
    )
    return cached_scheme_analysis.allele_counts.allele_count


def count_valid_alleles(scheme_directory_path: str | Path) -> int:
    """Count valid coding alleles across all loci in an extracted scheme.

    Args:
        scheme_directory_path: Directory containing an extracted MLST scheme.

    Returns:
        The number of alleles accepted by :func:`is_valid_allele_sequence`.

    Raises:
        FileNotFoundError: If ``scheme_directory_path`` does not exist.
        NotADirectoryError: If ``scheme_directory_path`` is not a directory.
        OSError: If an allele file cannot be read or decompressed.
        UnicodeError: If an allele file is not valid UTF-8 text.
    """
    resolved_scheme_directory_path = _resolve_scheme_directory_path(
        scheme_directory_path
    )
    cached_scheme_analysis = _get_or_create_cached_scheme_analysis(
        resolved_scheme_directory_path
    )
    return cached_scheme_analysis.allele_counts.valid_allele_count


def count_invalid_alleles(scheme_directory_path: str | Path) -> int:
    """Return the sum of all coding-sequence invalidity category counts.

    An allele that fails multiple rules contributes once to each failed
    category, so this value can exceed the number of unique invalid alleles.

    Args:
        scheme_directory_path: Directory containing an extracted MLST scheme.

    Returns:
        The sum of the incomplete ORF, ambiguous bases, missing start/stop
        codon, and in-frame stop codon counts.

    Raises:
        FileNotFoundError: If ``scheme_directory_path`` does not exist.
        NotADirectoryError: If ``scheme_directory_path`` is not a directory.
        OSError: If an allele file cannot be read or decompressed.
        UnicodeError: If an allele file is not valid UTF-8 text.
    """
    resolved_scheme_directory_path = _resolve_scheme_directory_path(
        scheme_directory_path
    )
    cached_scheme_analysis = _get_or_create_cached_scheme_analysis(
        resolved_scheme_directory_path
    )
    return cached_scheme_analysis.allele_counts.invalid_allele_count


def count_alleles_below_minimum_length(
    scheme_directory_path: str | Path,
    minimum_sequence_length: int,
) -> int:
    """Count alleles shorter than a specified nucleotide length.

    This focused function only evaluates sequence lengths. It is useful for
    reactive interfaces where the threshold changes independently of the
    scheme and the other validation metrics should not be recalculated.

    Args:
        scheme_directory_path: Directory containing an extracted MLST scheme.
        minimum_sequence_length: Minimum accepted allele length in nucleotides.
            The value must be zero or greater.

    Returns:
        The number of alleles shorter than ``minimum_sequence_length``.

    Raises:
        FileNotFoundError: If ``scheme_directory_path`` does not exist.
        NotADirectoryError: If ``scheme_directory_path`` is not a directory.
        OSError: If an allele file cannot be read or decompressed.
        UnicodeError: If an allele file is not valid UTF-8 text.
        ValueError: If ``minimum_sequence_length`` is negative.
    """
    if minimum_sequence_length < 0:
        raise ValueError("minimum_sequence_length cannot be negative")

    resolved_scheme_directory_path = _resolve_scheme_directory_path(
        scheme_directory_path
    )
    cached_scheme_analysis = _get_or_create_cached_scheme_analysis(
        resolved_scheme_directory_path
    )
    return sum(
        allele_length < minimum_sequence_length
        for allele_length in cached_scheme_analysis.allele_counts.allele_lengths
    )


def get_scheme_information(
    scheme_directory_path: str | Path,
    minimum_sequence_length: int = 0,
) -> SchemeInformation:
    """Collect locus and allele counts for an extracted scheme directory.

    The first call scans the directory and persists the static metrics and
    allele lengths. Later calls reuse that analysis, including when only the
    minimum sequence length changes.

    Args:
        scheme_directory_path: Directory containing an extracted MLST scheme.
        minimum_sequence_length: Minimum accepted allele length in nucleotides.
            Alleles shorter than this value are counted separately. The value
            must be zero or greater.

    Returns:
        An immutable :class:`SchemeInformation` containing all scheme-level
        counts and percentages.

    Raises:
        FileNotFoundError: If ``scheme_directory_path`` does not exist.
        NotADirectoryError: If ``scheme_directory_path`` is not a directory.
        OSError: If an allele file cannot be read or decompressed.
        UnicodeError: If an allele file is not valid UTF-8 text.
        ValueError: If ``minimum_sequence_length`` is negative.
    """
    if minimum_sequence_length < 0:
        raise ValueError("minimum_sequence_length cannot be negative")

    resolved_scheme_directory_path = _resolve_scheme_directory_path(
        scheme_directory_path
    )
    cached_scheme_analysis = _get_or_create_cached_scheme_analysis(
        resolved_scheme_directory_path
    )
    allele_counts = cached_scheme_analysis.allele_counts
    alleles_below_minimum_length_count = sum(
        allele_length < minimum_sequence_length
        for allele_length in allele_counts.allele_lengths
    )
    if allele_counts.allele_count == 0:
        alleles_below_locus_mode_threshold_percentage = 0.0
        alleles_above_locus_mode_threshold_percentage = 0.0
    else:
        alleles_below_locus_mode_threshold_percentage = (
            cached_scheme_analysis.alleles_below_locus_mode_threshold_count
            / allele_counts.allele_count
            * 100
        )
        alleles_above_locus_mode_threshold_percentage = (
            cached_scheme_analysis.alleles_above_locus_mode_threshold_count
            / allele_counts.allele_count
            * 100
        )
    scheme_nucleotide_base_composition = NucleotideBaseComposition(
        guanine_cytosine_base_count=(
            allele_counts.guanine_cytosine_base_count
        ),
        canonical_nucleotide_base_count=(
            allele_counts.canonical_nucleotide_base_count
        ),
        excluded_nucleotide_base_count=(
            allele_counts.excluded_nucleotide_base_count
        ),
    )

    return SchemeInformation(
        locus_count=cached_scheme_analysis.locus_count,
        allele_count=allele_counts.allele_count,
        valid_allele_count=allele_counts.valid_allele_count,
        invalid_allele_count=allele_counts.invalid_allele_count,
        incomplete_orf_allele_count=allele_counts.incomplete_orf_allele_count,
        ambiguous_bases_allele_count=(
            allele_counts.ambiguous_bases_allele_count
        ),
        missing_start_or_stop_codon_allele_count=(
            allele_counts.missing_start_or_stop_codon_allele_count
        ),
        in_frame_stop_codon_allele_count=(
            allele_counts.in_frame_stop_codon_allele_count
        ),
        alleles_below_minimum_length_count=alleles_below_minimum_length_count,
        minimum_sequence_length=minimum_sequence_length,
        alleles_below_locus_mode_threshold_count=(
            cached_scheme_analysis.alleles_below_locus_mode_threshold_count
        ),
        alleles_above_locus_mode_threshold_count=(
            cached_scheme_analysis.alleles_above_locus_mode_threshold_count
        ),
        alleles_below_locus_mode_threshold_percentage=(
            alleles_below_locus_mode_threshold_percentage
        ),
        alleles_above_locus_mode_threshold_percentage=(
            alleles_above_locus_mode_threshold_percentage
        ),
        guanine_cytosine_content_percentage=(
            scheme_nucleotide_base_composition.guanine_cytosine_content_percentage
        ),
        canonical_nucleotide_base_count=(
            scheme_nucleotide_base_composition.canonical_nucleotide_base_count
        ),
        excluded_nucleotide_base_count=(
            scheme_nucleotide_base_composition.excluded_nucleotide_base_count
        ),
    )


def _parse_non_negative_integer(argument_value: str) -> int:
    """Parse a non-negative integer for a command-line argument."""
    parsed_integer = int(argument_value)
    if parsed_integer < 0:
        raise argparse.ArgumentTypeError("value must be zero or greater")
    return parsed_integer


def main() -> None:
    """Parse command-line arguments and print scheme information."""
    argument_parser = argparse.ArgumentParser(
        description="Count loci and alleles in an extracted MLST scheme."
    )
    argument_parser.add_argument(
        "scheme_directory",
        type=Path,
        help="directory containing an extracted MLST scheme",
    )
    argument_parser.add_argument(
        "--json",
        action="store_true",
        dest="print_json",
        help="print the result as JSON",
    )
    argument_parser.add_argument(
        "--ml",
        type=_parse_non_negative_integer,
        default=0,
        dest="minimum_sequence_length",
        metavar="BASE_PAIRS",
        help="minimum allele length in nucleotides (default: 0)",
    )
    command_line_arguments = argument_parser.parse_args()

    resolved_scheme_directory_path = _resolve_scheme_directory_path(
        command_line_arguments.scheme_directory
    )
    scheme_information = get_scheme_information(
        resolved_scheme_directory_path,
        minimum_sequence_length=command_line_arguments.minimum_sequence_length,
    )

    if command_line_arguments.print_json:
        json_output = {
            "scheme_directory": str(resolved_scheme_directory_path),
            **asdict(scheme_information),
        }
        print(json.dumps(json_output, indent=2, sort_keys=True))
        return

    print(f"Scheme directory: {resolved_scheme_directory_path}")
    print(f"Loci: {scheme_information.locus_count}")
    print(f"Alleles: {scheme_information.allele_count}")
    print(f"Valid alleles: {scheme_information.valid_allele_count}")
    print(f"Invalid alleles: {scheme_information.invalid_allele_count}")
    print(f"Incomplete ORF: {scheme_information.incomplete_orf_allele_count}")
    print(f"Ambiguous bases: {scheme_information.ambiguous_bases_allele_count}")
    print(
        "Missing start/stop codon: "
        f"{scheme_information.missing_start_or_stop_codon_allele_count}"
    )
    print(
        "In-frame stop codon: "
        f"{scheme_information.in_frame_stop_codon_allele_count}"
    )
    guanine_cytosine_content_label = (
        "unavailable"
        if scheme_information.guanine_cytosine_content_percentage is None
        else f"{scheme_information.guanine_cytosine_content_percentage:.2f}%"
    )
    print(f"GC content: {guanine_cytosine_content_label}")
    print(
        "Canonical bases in GC denominator: "
        f"{scheme_information.canonical_nucleotide_base_count}"
    )
    print(
        "Noncanonical bases excluded from GC denominator: "
        f"{scheme_information.excluded_nucleotide_base_count}"
    )
    print(
        f"Alleles <{scheme_information.minimum_sequence_length} bp: "
        f"{scheme_information.alleles_below_minimum_length_count}"
    )
    print(
        "Alleles <80% of locus mode: "
        f"{scheme_information.alleles_below_locus_mode_threshold_count} "
        "("
        f"{scheme_information.alleles_below_locus_mode_threshold_percentage:.2f}%"
        ")"
    )
    print(
        "Alleles >120% of locus mode: "
        f"{scheme_information.alleles_above_locus_mode_threshold_count} "
        "("
        f"{scheme_information.alleles_above_locus_mode_threshold_percentage:.2f}%"
        ")"
    )


if __name__ == "__main__":
    main()
