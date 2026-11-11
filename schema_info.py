"""Inspect an extracted MLST scheme and report basic sequence statistics.

This module counts loci and evaluates allele coding sequences in a directory
produced by ``schema_cache.cache_archive``. A locus is represented by one FASTA
file, and an allele is represented by one FASTA record within those files. The
evaluation reports valid alleles, individual invalidity categories, and alleles
shorter than a configurable minimum sequence length. It also reports
scheme-wide counts and percentages of alleles below 80% or above 120% of their
respective locus's modal length. Common FASTA filename extensions are
recognized in both uncompressed and gzip-compressed form.

Notebook code can call :func:`count_loci`, :func:`count_alleles`, or
:func:`get_schema_information`. The module can also be executed as a
command-line program and can print either a readable summary or JSON.
"""

import argparse
import gzip
import json
import tempfile
from collections import Counter
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path

FASTA_FILENAME_SUFFIXES = (
    ".fasta",
    ".fas",
    ".fa",
    ".fna",
    ".fasta.gz",
    ".fas.gz",
    ".fa.gz",
    ".fna.gz",
)
VALID_DNA_BASES = frozenset("ACGT")
BACTERIAL_START_CODONS = frozenset({"ATG", "GTG", "TTG"})
STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
SCHEMA_INFORMATION_CACHE_FILENAME = ".schema-information-cache.json"
SCHEMA_INFORMATION_CACHE_VERSION = 3

__all__ = [
    "BACTERIAL_START_CODONS",
    "FASTA_FILENAME_SUFFIXES",
    "SCHEMA_INFORMATION_CACHE_FILENAME",
    "SCHEMA_INFORMATION_CACHE_VERSION",
    "STOP_CODONS",
    "VALID_DNA_BASES",
    "SchemaInformation",
    "count_alleles",
    "count_alleles_below_minimum_length",
    "count_invalid_alleles",
    "count_loci",
    "count_valid_alleles",
    "get_schema_information",
    "is_valid_allele_sequence",
]


@dataclass(frozen=True)
class SchemaInformation:
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
class _CachedSchemaAnalysis:
    """Schema metrics and allele lengths persisted in the extracted schema."""

    locus_count: int
    allele_counts: _AlleleCounts
    alleles_below_locus_mode_threshold_count: int
    alleles_above_locus_mode_threshold_count: int


def _resolve_schema_directory_path(schema_directory_path: str | Path) -> Path:
    """Resolve and validate an extracted schema directory path.

    Args:
        schema_directory_path: Directory containing the extracted schema.

    Returns:
        The absolute path to the schema directory.

    Raises:
        FileNotFoundError: If the supplied path does not exist.
        NotADirectoryError: If the supplied path is not a directory.
    """
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


def _is_fasta_file(allele_file_path: Path) -> bool:
    """Return whether a path has a recognized FASTA filename suffix."""
    lowercase_filename = allele_file_path.name.lower()
    return allele_file_path.is_file() and lowercase_filename.endswith(
        FASTA_FILENAME_SUFFIXES
    )


def _find_allele_file_paths(schema_directory_path: Path) -> tuple[Path, ...]:
    """Return all recognized FASTA files below a schema directory."""
    return tuple(
        sorted(
            candidate_file_path
            for candidate_file_path in schema_directory_path.rglob("*")
            if _is_fasta_file(candidate_file_path)
        )
    )


def _iterate_fasta_sequences(fasta_lines: Iterable[str]) -> Iterable[str]:
    """Yield one assembled nucleotide sequence per FASTA header."""
    current_sequence_parts: list[str] | None = None

    for fasta_line in fasta_lines:
        stripped_fasta_line = fasta_line.strip()

        if stripped_fasta_line.startswith(">"):
            if current_sequence_parts is not None:
                yield "".join(current_sequence_parts)
            current_sequence_parts = []
            continue

        if current_sequence_parts is not None:
            current_sequence_parts.append(stripped_fasta_line)

    if current_sequence_parts is not None:
        yield "".join(current_sequence_parts)


def _contains_in_frame_stop_codon(normalized_allele_sequence: str) -> bool:
    """Return whether an allele has a stop codon before its valid endpoint."""
    final_complete_codon_end_position = len(normalized_allele_sequence) - (
        len(normalized_allele_sequence) % 3
    )
    in_frame_codon_end_position = final_complete_codon_end_position

    has_terminal_stop_codon = (
        len(normalized_allele_sequence) % 3 == 0
        and normalized_allele_sequence[-3:] in STOP_CODONS
    )
    if has_terminal_stop_codon:
        in_frame_codon_end_position -= 3

    return any(
        normalized_allele_sequence[
            codon_start_position : codon_start_position + 3
        ]
        in STOP_CODONS
        for codon_start_position in range(
            3,
            in_frame_codon_end_position,
            3,
        )
    )


def _evaluate_allele_sequence(
    allele_sequence: str,
    minimum_sequence_length: int,
) -> _AlleleValidationResult:
    """Evaluate all coding-sequence and minimum-length rules for one allele."""
    normalized_allele_sequence = "".join(allele_sequence.split()).upper()
    has_incomplete_orf = len(normalized_allele_sequence) % 3 != 0
    has_ambiguous_bases = not set(normalized_allele_sequence).issubset(
        VALID_DNA_BASES
    )
    is_missing_start_or_stop_codon = (
        normalized_allele_sequence[:3] not in BACTERIAL_START_CODONS
        or normalized_allele_sequence[-3:] not in STOP_CODONS
    )

    return _AlleleValidationResult(
        has_incomplete_orf=has_incomplete_orf,
        has_ambiguous_bases=has_ambiguous_bases,
        is_missing_start_or_stop_codon=is_missing_start_or_stop_codon,
        has_in_frame_stop_codon=_contains_in_frame_stop_codon(
            normalized_allele_sequence
        ),
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

    for allele_sequence in _iterate_fasta_sequences(fasta_lines):
        normalized_allele_sequence = "".join(allele_sequence.split())
        allele_lengths.append(len(normalized_allele_sequence))
        allele_count += 1
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


def _read_cached_schema_analysis(
    schema_directory_path: Path,
) -> _CachedSchemaAnalysis | None:
    """Read a compatible cached analysis, returning ``None`` if unavailable."""
    schema_information_cache_file_path = (
        schema_directory_path / SCHEMA_INFORMATION_CACHE_FILENAME
    )

    try:
        cached_analysis_data = json.loads(
            schema_information_cache_file_path.read_text(encoding="utf-8")
        )
        if (
            cached_analysis_data["cache_version"]
            != SCHEMA_INFORMATION_CACHE_VERSION
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

    return _CachedSchemaAnalysis(
        locus_count=locus_count,
        allele_counts=allele_counts,
        alleles_below_locus_mode_threshold_count=(
            alleles_below_locus_mode_threshold_count
        ),
        alleles_above_locus_mode_threshold_count=(
            alleles_above_locus_mode_threshold_count
        ),
    )


def _write_cached_schema_analysis(
    schema_directory_path: Path,
    cached_schema_analysis: _CachedSchemaAnalysis,
) -> None:
    """Atomically persist schema metrics and allele lengths as JSON."""
    allele_counts = cached_schema_analysis.allele_counts
    cached_analysis_data = {
        "cache_version": SCHEMA_INFORMATION_CACHE_VERSION,
        "locus_count": cached_schema_analysis.locus_count,
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
        "alleles_below_locus_mode_threshold_count": (
            cached_schema_analysis.alleles_below_locus_mode_threshold_count
        ),
        "alleles_above_locus_mode_threshold_count": (
            cached_schema_analysis.alleles_above_locus_mode_threshold_count
        ),
    }
    schema_information_cache_file_path = (
        schema_directory_path / SCHEMA_INFORMATION_CACHE_FILENAME
    )

    temporary_cache_file_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=schema_directory_path,
            prefix=f"{SCHEMA_INFORMATION_CACHE_FILENAME}.",
            delete=False,
        ) as temporary_cache_file:
            json.dump(cached_analysis_data, temporary_cache_file, separators=(",", ":"))
            temporary_cache_file_path = Path(temporary_cache_file.name)

        temporary_cache_file_path.replace(schema_information_cache_file_path)
    finally:
        if temporary_cache_file_path is not None:
            temporary_cache_file_path.unlink(missing_ok=True)


def _get_or_create_cached_schema_analysis(
    schema_directory_path: Path,
) -> _CachedSchemaAnalysis:
    """Return cached schema metrics, calculating and persisting them once."""
    cached_schema_analysis = _read_cached_schema_analysis(schema_directory_path)
    if cached_schema_analysis is not None:
        return cached_schema_analysis

    allele_file_paths = _find_allele_file_paths(schema_directory_path)
    (
        allele_counts,
        alleles_below_locus_mode_threshold_count,
        alleles_above_locus_mode_threshold_count,
    ) = _count_alleles_in_files(
        allele_file_paths,
        minimum_sequence_length=0,
    )

    cached_schema_analysis = _CachedSchemaAnalysis(
        locus_count=len(allele_file_paths),
        allele_counts=allele_counts,
        alleles_below_locus_mode_threshold_count=(
            alleles_below_locus_mode_threshold_count
        ),
        alleles_above_locus_mode_threshold_count=(
            alleles_above_locus_mode_threshold_count
        ),
    )

    try:
        _write_cached_schema_analysis(
            schema_directory_path,
            cached_schema_analysis,
        )
    except OSError:
        # Read-only schema directories can still be analyzed without caching.
        pass

    return cached_schema_analysis


def count_loci(schema_directory_path: str | Path) -> int:
    """Count locus FASTA files in an extracted schema directory.

    Args:
        schema_directory_path: Directory containing an extracted MLST scheme.

    Returns:
        The number of recognized FASTA files found recursively.

    Raises:
        FileNotFoundError: If ``schema_directory_path`` does not exist.
        NotADirectoryError: If ``schema_directory_path`` is not a directory.
    """
    resolved_schema_directory_path = _resolve_schema_directory_path(
        schema_directory_path
    )
    cached_schema_analysis = _get_or_create_cached_schema_analysis(
        resolved_schema_directory_path
    )
    return cached_schema_analysis.locus_count


def count_alleles(schema_directory_path: str | Path) -> int:
    """Count FASTA records across all loci in an extracted schema directory.

    Args:
        schema_directory_path: Directory containing an extracted MLST scheme.

    Returns:
        The total number of FASTA records across all recognized locus files.

    Raises:
        FileNotFoundError: If ``schema_directory_path`` does not exist.
        NotADirectoryError: If ``schema_directory_path`` is not a directory.
        OSError: If an allele file cannot be read or decompressed.
        UnicodeError: If an allele file is not valid UTF-8 text.
    """
    resolved_schema_directory_path = _resolve_schema_directory_path(
        schema_directory_path
    )
    cached_schema_analysis = _get_or_create_cached_schema_analysis(
        resolved_schema_directory_path
    )
    return cached_schema_analysis.allele_counts.allele_count


def count_valid_alleles(schema_directory_path: str | Path) -> int:
    """Count valid coding alleles across all loci in an extracted schema.

    Args:
        schema_directory_path: Directory containing an extracted MLST scheme.

    Returns:
        The number of alleles accepted by :func:`is_valid_allele_sequence`.

    Raises:
        FileNotFoundError: If ``schema_directory_path`` does not exist.
        NotADirectoryError: If ``schema_directory_path`` is not a directory.
        OSError: If an allele file cannot be read or decompressed.
        UnicodeError: If an allele file is not valid UTF-8 text.
    """
    resolved_schema_directory_path = _resolve_schema_directory_path(
        schema_directory_path
    )
    cached_schema_analysis = _get_or_create_cached_schema_analysis(
        resolved_schema_directory_path
    )
    return cached_schema_analysis.allele_counts.valid_allele_count


def count_invalid_alleles(schema_directory_path: str | Path) -> int:
    """Return the sum of all coding-sequence invalidity category counts.

    An allele that fails multiple rules contributes once to each failed
    category, so this value can exceed the number of unique invalid alleles.

    Args:
        schema_directory_path: Directory containing an extracted MLST scheme.

    Returns:
        The sum of the incomplete ORF, ambiguous bases, missing start/stop
        codon, and in-frame stop codon counts.

    Raises:
        FileNotFoundError: If ``schema_directory_path`` does not exist.
        NotADirectoryError: If ``schema_directory_path`` is not a directory.
        OSError: If an allele file cannot be read or decompressed.
        UnicodeError: If an allele file is not valid UTF-8 text.
    """
    resolved_schema_directory_path = _resolve_schema_directory_path(
        schema_directory_path
    )
    cached_schema_analysis = _get_or_create_cached_schema_analysis(
        resolved_schema_directory_path
    )
    return cached_schema_analysis.allele_counts.invalid_allele_count


def count_alleles_below_minimum_length(
    schema_directory_path: str | Path,
    minimum_sequence_length: int,
) -> int:
    """Count alleles shorter than a specified nucleotide length.

    This focused function only evaluates sequence lengths. It is useful for
    reactive interfaces where the threshold changes independently of the
    schema and the other validation metrics should not be recalculated.

    Args:
        schema_directory_path: Directory containing an extracted MLST scheme.
        minimum_sequence_length: Minimum accepted allele length in nucleotides.
            The value must be zero or greater.

    Returns:
        The number of alleles shorter than ``minimum_sequence_length``.

    Raises:
        FileNotFoundError: If ``schema_directory_path`` does not exist.
        NotADirectoryError: If ``schema_directory_path`` is not a directory.
        OSError: If an allele file cannot be read or decompressed.
        UnicodeError: If an allele file is not valid UTF-8 text.
        ValueError: If ``minimum_sequence_length`` is negative.
    """
    if minimum_sequence_length < 0:
        raise ValueError("minimum_sequence_length cannot be negative")

    resolved_schema_directory_path = _resolve_schema_directory_path(
        schema_directory_path
    )
    cached_schema_analysis = _get_or_create_cached_schema_analysis(
        resolved_schema_directory_path
    )
    return sum(
        allele_length < minimum_sequence_length
        for allele_length in cached_schema_analysis.allele_counts.allele_lengths
    )


def get_schema_information(
    schema_directory_path: str | Path,
    minimum_sequence_length: int = 0,
) -> SchemaInformation:
    """Collect locus and allele counts for an extracted schema directory.

    The first call scans the directory and persists the static metrics and
    allele lengths. Later calls reuse that analysis, including when only the
    minimum sequence length changes.

    Args:
        schema_directory_path: Directory containing an extracted MLST scheme.
        minimum_sequence_length: Minimum accepted allele length in nucleotides.
            Alleles shorter than this value are counted separately. The value
            must be zero or greater.

    Returns:
        An immutable :class:`SchemaInformation` containing all scheme-level
        counts and percentages.

    Raises:
        FileNotFoundError: If ``schema_directory_path`` does not exist.
        NotADirectoryError: If ``schema_directory_path`` is not a directory.
        OSError: If an allele file cannot be read or decompressed.
        UnicodeError: If an allele file is not valid UTF-8 text.
        ValueError: If ``minimum_sequence_length`` is negative.
    """
    if minimum_sequence_length < 0:
        raise ValueError("minimum_sequence_length cannot be negative")

    resolved_schema_directory_path = _resolve_schema_directory_path(
        schema_directory_path
    )
    cached_schema_analysis = _get_or_create_cached_schema_analysis(
        resolved_schema_directory_path
    )
    allele_counts = cached_schema_analysis.allele_counts
    alleles_below_minimum_length_count = sum(
        allele_length < minimum_sequence_length
        for allele_length in allele_counts.allele_lengths
    )
    if allele_counts.allele_count == 0:
        alleles_below_locus_mode_threshold_percentage = 0.0
        alleles_above_locus_mode_threshold_percentage = 0.0
    else:
        alleles_below_locus_mode_threshold_percentage = (
            cached_schema_analysis.alleles_below_locus_mode_threshold_count
            / allele_counts.allele_count
            * 100
        )
        alleles_above_locus_mode_threshold_percentage = (
            cached_schema_analysis.alleles_above_locus_mode_threshold_count
            / allele_counts.allele_count
            * 100
        )

    return SchemaInformation(
        locus_count=cached_schema_analysis.locus_count,
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
            cached_schema_analysis.alleles_below_locus_mode_threshold_count
        ),
        alleles_above_locus_mode_threshold_count=(
            cached_schema_analysis.alleles_above_locus_mode_threshold_count
        ),
        alleles_below_locus_mode_threshold_percentage=(
            alleles_below_locus_mode_threshold_percentage
        ),
        alleles_above_locus_mode_threshold_percentage=(
            alleles_above_locus_mode_threshold_percentage
        ),
    )


def _parse_non_negative_integer(argument_value: str) -> int:
    """Parse a non-negative integer for a command-line argument."""
    parsed_integer = int(argument_value)
    if parsed_integer < 0:
        raise argparse.ArgumentTypeError("value must be zero or greater")
    return parsed_integer


def main() -> None:
    """Parse command-line arguments and print schema information."""
    argument_parser = argparse.ArgumentParser(
        description="Count loci and alleles in an extracted MLST scheme."
    )
    argument_parser.add_argument(
        "schema_directory",
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

    resolved_schema_directory_path = _resolve_schema_directory_path(
        command_line_arguments.schema_directory
    )
    schema_information = get_schema_information(
        resolved_schema_directory_path,
        minimum_sequence_length=command_line_arguments.minimum_sequence_length,
    )

    if command_line_arguments.print_json:
        json_output = {
            "schema_directory": str(resolved_schema_directory_path),
            **asdict(schema_information),
        }
        print(json.dumps(json_output, indent=2, sort_keys=True))
        return

    print(f"Schema directory: {resolved_schema_directory_path}")
    print(f"Loci: {schema_information.locus_count}")
    print(f"Alleles: {schema_information.allele_count}")
    print(f"Valid alleles: {schema_information.valid_allele_count}")
    print(f"Invalid alleles: {schema_information.invalid_allele_count}")
    print(f"Incomplete ORF: {schema_information.incomplete_orf_allele_count}")
    print(f"Ambiguous bases: {schema_information.ambiguous_bases_allele_count}")
    print(
        "Missing start/stop codon: "
        f"{schema_information.missing_start_or_stop_codon_allele_count}"
    )
    print(
        "In-frame stop codon: "
        f"{schema_information.in_frame_stop_codon_allele_count}"
    )
    print(
        f"Alleles <{schema_information.minimum_sequence_length} bp: "
        f"{schema_information.alleles_below_minimum_length_count}"
    )
    print(
        "Alleles <80% of locus mode: "
        f"{schema_information.alleles_below_locus_mode_threshold_count} "
        "("
        f"{schema_information.alleles_below_locus_mode_threshold_percentage:.2f}%"
        ")"
    )
    print(
        "Alleles >120% of locus mode: "
        f"{schema_information.alleles_above_locus_mode_threshold_count} "
        "("
        f"{schema_information.alleles_above_locus_mode_threshold_percentage:.2f}%"
        ")"
    )


if __name__ == "__main__":
    main()
