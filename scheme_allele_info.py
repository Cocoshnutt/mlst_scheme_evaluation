"""Calculate and cache chewBBACA-style allele information per locus.

This library analyzes every FASTA locus in an extracted MLST scheme and
returns the integrity and diversity fields used by the chewBBACA Schema
Evaluator allele-analysis table. Static sequence metrics and allele lengths
are cached inside the extracted scheme. Changing the minimum allele length
therefore recalculates only the per-locus short-allele counts and does not
reread FASTA sequences.

The public :func:`get_scheme_allele_information` function is intended for the
Marimo table in ``notebook.py`` and can also be imported by other Python code.
"""

import argparse
import gzip
import json
import re
import tempfile
from collections import Counter
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path

from schema_info import (
    BACTERIAL_START_CODONS,
    FASTA_FILENAME_SUFFIXES,
    STOP_CODONS,
    VALID_DNA_BASES,
)

SCHEME_ALLELE_INFORMATION_CACHE_FILENAME = ".scheme-allele-information-cache.json"
SCHEME_ALLELE_INFORMATION_CACHE_VERSION = 1

CODON_TO_AMINO_ACID = {
    "AAA": "K",
    "AAC": "N",
    "AAG": "K",
    "AAT": "N",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AGA": "R",
    "AGC": "S",
    "AGG": "R",
    "AGT": "S",
    "ATA": "I",
    "ATC": "I",
    "ATG": "M",
    "ATT": "I",
    "CAA": "Q",
    "CAC": "H",
    "CAG": "Q",
    "CAT": "H",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "GAA": "E",
    "GAC": "D",
    "GAG": "E",
    "GAT": "D",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "TAA": "*",
    "TAC": "Y",
    "TAG": "*",
    "TAT": "Y",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TGA": "*",
    "TGC": "C",
    "TGG": "W",
    "TGT": "C",
    "TTA": "L",
    "TTC": "F",
    "TTG": "L",
    "TTT": "F",
}

__all__ = [
    "CODON_TO_AMINO_ACID",
    "SCHEME_ALLELE_INFORMATION_CACHE_FILENAME",
    "SCHEME_ALLELE_INFORMATION_CACHE_VERSION",
    "LocusAlleleInformation",
    "get_scheme_allele_information",
]


@dataclass(frozen=True)
class LocusAlleleInformation:
    """Integrity and diversity metrics for one locus.

    Attributes:
        locus_name: Schema-relative locus name without a FASTA suffix.
        total_allele_count: Total number of FASTA records in the locus.
        valid_allele_count: Alleles satisfying every coding-sequence rule.
        invalid_allele_count: Sum of the four invalidity category counts.
        validated_allele_proportion: Valid alleles divided by total alleles.
        distinct_protein_allele_count: Unique proteins translated from valid
            allele sequences.
        incomplete_orf_allele_count: Alleles whose size is not divisible by
            three.
        ambiguous_bases_allele_count: Alleles containing non-ACGT bases.
        missing_start_or_stop_codon_allele_count: Alleles missing a recognized
            start codon, terminal stop codon, or both.
        in_frame_stop_codon_allele_count: Alleles with a stop codon before the
            terminal codon.
        alleles_below_minimum_length_count: Alleles shorter than
            ``minimum_sequence_length``.
        minimum_sequence_length: Minimum length used for the short-allele
            count, in nucleotides.
        alleles_below_mode_threshold_count: Alleles strictly below 80% of the
            locus's modal allele length.
        alleles_above_mode_threshold_count: Alleles strictly above 120% of the
            locus's modal allele length.
        missing_allele_identifier_count: Missing positive integer identifiers
            in the range from one through the largest observed identifier.
            ``None`` indicates that at least one FASTA identifier could not be
            interpreted using chewBBACA's identifier convention.
    """

    locus_name: str
    total_allele_count: int
    valid_allele_count: int
    invalid_allele_count: int
    validated_allele_proportion: float
    distinct_protein_allele_count: int
    incomplete_orf_allele_count: int
    ambiguous_bases_allele_count: int
    missing_start_or_stop_codon_allele_count: int
    in_frame_stop_codon_allele_count: int
    alleles_below_minimum_length_count: int
    minimum_sequence_length: int
    alleles_below_mode_threshold_count: int
    alleles_above_mode_threshold_count: int
    missing_allele_identifier_count: int | None


@dataclass(frozen=True)
class _AlleleSequenceValidation:
    """Coding-sequence validation flags for one normalized allele."""

    has_incomplete_orf: bool
    has_ambiguous_bases: bool
    is_missing_start_or_stop_codon: bool
    has_in_frame_stop_codon: bool

    @property
    def is_valid(self) -> bool:
        """Return whether the allele has no coding-sequence exception."""
        return not any(
            (
                self.has_incomplete_orf,
                self.has_ambiguous_bases,
                self.is_missing_start_or_stop_codon,
                self.has_in_frame_stop_codon,
            )
        )


@dataclass(frozen=True)
class _CachedLocusAlleleInformation:
    """Static per-locus metrics and lengths stored in the persistent cache."""

    locus_name: str
    total_allele_count: int
    valid_allele_count: int
    distinct_protein_allele_count: int
    incomplete_orf_allele_count: int
    ambiguous_bases_allele_count: int
    missing_start_or_stop_codon_allele_count: int
    in_frame_stop_codon_allele_count: int
    alleles_below_mode_threshold_count: int
    alleles_above_mode_threshold_count: int
    missing_allele_identifier_count: int | None
    allele_lengths: tuple[int, ...]

    @property
    def invalid_allele_count(self) -> int:
        """Return the sum of all four invalidity category counts."""
        return (
            self.incomplete_orf_allele_count
            + self.ambiguous_bases_allele_count
            + self.missing_start_or_stop_codon_allele_count
            + self.in_frame_stop_codon_allele_count
        )


def _resolve_schema_directory_path(schema_directory_path: str | Path) -> Path:
    """Resolve and validate an extracted scheme directory."""
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


def _find_allele_file_paths(
    schema_directory_path: Path,
) -> tuple[Path, ...]:
    """Return recognized locus FASTA files in deterministic order."""
    return tuple(
        sorted(
            candidate_file_path
            for candidate_file_path in schema_directory_path.rglob("*")
            if _is_fasta_file(candidate_file_path)
        )
    )


def _remove_fasta_filename_suffix(relative_file_path: str) -> str:
    """Remove a recognized FASTA suffix from a schema-relative path."""
    lowercase_relative_file_path = relative_file_path.lower()
    for filename_suffix in sorted(FASTA_FILENAME_SUFFIXES, key=len, reverse=True):
        if lowercase_relative_file_path.endswith(filename_suffix):
            return relative_file_path[: -len(filename_suffix)]
    return relative_file_path


def _iterate_fasta_records(
    fasta_lines: Iterable[str],
) -> Iterable[tuple[str, str]]:
    """Yield identifier and assembled sequence pairs from FASTA text."""
    current_identifier: str | None = None
    current_sequence_parts: list[str] = []

    for fasta_line in fasta_lines:
        stripped_fasta_line = fasta_line.strip()
        if stripped_fasta_line.startswith(">"):
            if current_identifier is not None:
                yield current_identifier, "".join(current_sequence_parts)
            current_identifier = stripped_fasta_line[1:].split(maxsplit=1)[0]
            current_sequence_parts = []
        elif stripped_fasta_line and current_identifier is not None:
            current_sequence_parts.append(stripped_fasta_line)

    if current_identifier is not None:
        yield current_identifier, "".join(current_sequence_parts)


def _read_fasta_records(
    allele_file_path: Path,
) -> tuple[tuple[str, str], ...]:
    """Read every identifier and sequence from one locus FASTA file."""
    if allele_file_path.name.lower().endswith(".gz"):
        with gzip.open(
            allele_file_path,
            mode="rt",
            encoding="utf-8",
        ) as allele_file:
            return tuple(_iterate_fasta_records(allele_file))

    with allele_file_path.open(encoding="utf-8") as allele_file:
        return tuple(_iterate_fasta_records(allele_file))


def _contains_in_frame_stop_codon(normalized_allele_sequence: str) -> bool:
    """Return whether a stop codon occurs before the terminal codon."""
    final_complete_codon_end_position = len(normalized_allele_sequence) - (
        len(normalized_allele_sequence) % 3
    )
    in_frame_codon_end_position = final_complete_codon_end_position
    if (
        len(normalized_allele_sequence) % 3 == 0
        and normalized_allele_sequence[-3:] in STOP_CODONS
    ):
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
    normalized_allele_sequence: str,
) -> _AlleleSequenceValidation:
    """Evaluate the four chewBBACA-style coding-sequence rules."""
    return _AlleleSequenceValidation(
        has_incomplete_orf=len(normalized_allele_sequence) % 3 != 0,
        has_ambiguous_bases=not set(normalized_allele_sequence).issubset(
            VALID_DNA_BASES
        ),
        is_missing_start_or_stop_codon=(
            normalized_allele_sequence[:3] not in BACTERIAL_START_CODONS
            or normalized_allele_sequence[-3:] not in STOP_CODONS
        ),
        has_in_frame_stop_codon=_contains_in_frame_stop_codon(
            normalized_allele_sequence
        ),
    )


def _translate_valid_allele(normalized_allele_sequence: str) -> str:
    """Translate a valid bacterial CDS, normalizing its start codon to M."""
    translated_amino_acids = ["M"]
    translated_amino_acids.extend(
        CODON_TO_AMINO_ACID[
            normalized_allele_sequence[
                codon_start_position : codon_start_position + 3
            ]
        ]
        for codon_start_position in range(
            3,
            len(normalized_allele_sequence) - 3,
            3,
        )
    )
    return "".join(translated_amino_acids)


def _extract_allele_identifier(fasta_identifier: str) -> int | None:
    """Extract a positive integer allele ID using chewBBACA conventions."""
    identifier_match = re.search(r"(?:^|[_-])\*?(\d+)$", fasta_identifier)
    if identifier_match is None:
        return None
    allele_identifier = int(identifier_match.group(1))
    return allele_identifier if allele_identifier > 0 else None


def _count_missing_allele_identifiers(
    fasta_identifiers: Iterable[str],
) -> int | None:
    """Count absent IDs from one through the largest observed allele ID."""
    allele_identifiers: list[int] = []
    for fasta_identifier in fasta_identifiers:
        allele_identifier = _extract_allele_identifier(fasta_identifier)
        if allele_identifier is None:
            return None
        allele_identifiers.append(allele_identifier)

    if not allele_identifiers:
        return 0
    unique_allele_identifiers = set(allele_identifiers)
    return max(unique_allele_identifiers) - len(unique_allele_identifiers)


def _count_length_outliers(allele_lengths: tuple[int, ...]) -> tuple[int, int]:
    """Count lengths strictly outside 80–120% of the shortest tied mode."""
    if not allele_lengths:
        return 0, 0
    length_frequencies = Counter(allele_lengths)
    highest_length_frequency = max(length_frequencies.values())
    modal_allele_length = min(
        allele_length
        for allele_length, allele_length_frequency in length_frequencies.items()
        if allele_length_frequency == highest_length_frequency
    )
    lower_length_threshold = modal_allele_length * 0.8
    upper_length_threshold = modal_allele_length * 1.2
    return (
        sum(
            allele_length < lower_length_threshold
            for allele_length in allele_lengths
        ),
        sum(
            allele_length > upper_length_threshold
            for allele_length in allele_lengths
        ),
    )


def _analyze_locus(
    schema_directory_path: Path,
    allele_file_path: Path,
) -> _CachedLocusAlleleInformation:
    """Calculate every static allele-analysis metric for one locus."""
    fasta_records = _read_fasta_records(allele_file_path)
    allele_lengths: list[int] = []
    valid_allele_count = 0
    incomplete_orf_allele_count = 0
    ambiguous_bases_allele_count = 0
    missing_start_or_stop_codon_allele_count = 0
    in_frame_stop_codon_allele_count = 0
    distinct_protein_alleles: set[str] = set()

    for _fasta_identifier, allele_sequence in fasta_records:
        normalized_allele_sequence = "".join(allele_sequence.split()).upper()
        allele_lengths.append(len(normalized_allele_sequence))
        allele_sequence_validation = _evaluate_allele_sequence(
            normalized_allele_sequence
        )
        if allele_sequence_validation.is_valid:
            valid_allele_count += 1
            distinct_protein_alleles.add(
                _translate_valid_allele(normalized_allele_sequence)
            )
        if allele_sequence_validation.has_incomplete_orf:
            incomplete_orf_allele_count += 1
        if allele_sequence_validation.has_ambiguous_bases:
            ambiguous_bases_allele_count += 1
        if allele_sequence_validation.is_missing_start_or_stop_codon:
            missing_start_or_stop_codon_allele_count += 1
        if allele_sequence_validation.has_in_frame_stop_codon:
            in_frame_stop_codon_allele_count += 1

    immutable_allele_lengths = tuple(allele_lengths)
    (
        alleles_below_mode_threshold_count,
        alleles_above_mode_threshold_count,
    ) = _count_length_outliers(immutable_allele_lengths)
    relative_allele_file_path = allele_file_path.relative_to(
        schema_directory_path
    ).as_posix()
    return _CachedLocusAlleleInformation(
        locus_name=_remove_fasta_filename_suffix(relative_allele_file_path),
        total_allele_count=len(fasta_records),
        valid_allele_count=valid_allele_count,
        distinct_protein_allele_count=len(distinct_protein_alleles),
        incomplete_orf_allele_count=incomplete_orf_allele_count,
        ambiguous_bases_allele_count=ambiguous_bases_allele_count,
        missing_start_or_stop_codon_allele_count=(
            missing_start_or_stop_codon_allele_count
        ),
        in_frame_stop_codon_allele_count=in_frame_stop_codon_allele_count,
        alleles_below_mode_threshold_count=alleles_below_mode_threshold_count,
        alleles_above_mode_threshold_count=alleles_above_mode_threshold_count,
        missing_allele_identifier_count=_count_missing_allele_identifiers(
            fasta_identifier for fasta_identifier, _sequence in fasta_records
        ),
        allele_lengths=immutable_allele_lengths,
    )


def _create_source_signature(
    schema_directory_path: Path,
    allele_file_paths: Iterable[Path],
) -> tuple[tuple[str, int, int], ...]:
    """Create a file metadata signature for cache invalidation."""
    source_signature: list[tuple[str, int, int]] = []
    for allele_file_path in allele_file_paths:
        allele_file_status = allele_file_path.stat()
        source_signature.append(
            (
                allele_file_path.relative_to(
                    schema_directory_path
                ).as_posix(),
                allele_file_status.st_size,
                allele_file_status.st_mtime_ns,
            )
        )
    return tuple(source_signature)


def _optional_integer(value: object) -> int | None:
    """Convert a cached optional integer value."""
    return None if value is None else int(value)


def _read_cached_analysis(
    schema_directory_path: Path,
    source_signature: tuple[tuple[str, int, int], ...],
) -> tuple[_CachedLocusAlleleInformation, ...] | None:
    """Read a compatible cached allele analysis, or return ``None``."""
    cache_file_path = (
        schema_directory_path / SCHEME_ALLELE_INFORMATION_CACHE_FILENAME
    )
    try:
        cached_data = json.loads(cache_file_path.read_text(encoding="utf-8"))
        if (
            cached_data["cache_version"]
            != SCHEME_ALLELE_INFORMATION_CACHE_VERSION
        ):
            return None
        cached_source_signature = tuple(
            (str(relative_path), int(file_size), int(modification_time))
            for relative_path, file_size, modification_time in cached_data[
                "source_signature"
            ]
        )
        if cached_source_signature != source_signature:
            return None

        return tuple(
            _CachedLocusAlleleInformation(
                locus_name=str(locus_data["locus_name"]),
                total_allele_count=int(locus_data["total_allele_count"]),
                valid_allele_count=int(locus_data["valid_allele_count"]),
                distinct_protein_allele_count=int(
                    locus_data["distinct_protein_allele_count"]
                ),
                incomplete_orf_allele_count=int(
                    locus_data["incomplete_orf_allele_count"]
                ),
                ambiguous_bases_allele_count=int(
                    locus_data["ambiguous_bases_allele_count"]
                ),
                missing_start_or_stop_codon_allele_count=int(
                    locus_data["missing_start_or_stop_codon_allele_count"]
                ),
                in_frame_stop_codon_allele_count=int(
                    locus_data["in_frame_stop_codon_allele_count"]
                ),
                alleles_below_mode_threshold_count=int(
                    locus_data["alleles_below_mode_threshold_count"]
                ),
                alleles_above_mode_threshold_count=int(
                    locus_data["alleles_above_mode_threshold_count"]
                ),
                missing_allele_identifier_count=_optional_integer(
                    locus_data["missing_allele_identifier_count"]
                ),
                allele_lengths=tuple(
                    int(allele_length)
                    for allele_length in locus_data["allele_lengths"]
                ),
            )
            for locus_data in cached_data["locus_information"]
        )
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return None


def _write_cached_analysis(
    schema_directory_path: Path,
    source_signature: tuple[tuple[str, int, int], ...],
    locus_information: tuple[_CachedLocusAlleleInformation, ...],
) -> None:
    """Atomically persist static allele-analysis results as JSON."""
    cached_data = {
        "cache_version": SCHEME_ALLELE_INFORMATION_CACHE_VERSION,
        "source_signature": source_signature,
        "locus_information": [
            asdict(locus_information_item)
            for locus_information_item in locus_information
        ],
    }
    temporary_cache_file_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=schema_directory_path,
            prefix=f"{SCHEME_ALLELE_INFORMATION_CACHE_FILENAME}.",
            delete=False,
        ) as temporary_cache_file:
            json.dump(cached_data, temporary_cache_file, separators=(",", ":"))
            temporary_cache_file_path = Path(temporary_cache_file.name)

        temporary_cache_file_path.replace(
            schema_directory_path / SCHEME_ALLELE_INFORMATION_CACHE_FILENAME
        )
    finally:
        if temporary_cache_file_path is not None:
            temporary_cache_file_path.unlink(missing_ok=True)


def _get_or_create_cached_analysis(
    schema_directory_path: Path,
) -> tuple[_CachedLocusAlleleInformation, ...]:
    """Return cached static metrics, scanning FASTA files only when needed."""
    allele_file_paths = _find_allele_file_paths(schema_directory_path)
    source_signature = _create_source_signature(
        schema_directory_path,
        allele_file_paths,
    )
    cached_analysis = _read_cached_analysis(
        schema_directory_path,
        source_signature,
    )
    if cached_analysis is not None:
        return cached_analysis

    cached_analysis = tuple(
        _analyze_locus(schema_directory_path, allele_file_path)
        for allele_file_path in allele_file_paths
    )
    try:
        _write_cached_analysis(
            schema_directory_path,
            source_signature,
            cached_analysis,
        )
    except OSError:
        # Read-only schema directories can still be analyzed without caching.
        pass
    return cached_analysis


def get_scheme_allele_information(
    schema_directory_path: str | Path,
    minimum_sequence_length: int = 0,
) -> tuple[LocusAlleleInformation, ...]:
    """Return chewBBACA-style allele-analysis fields for every locus.

    Args:
        schema_directory_path: Directory containing an extracted MLST scheme.
        minimum_sequence_length: Minimum accepted allele length in nucleotides.
            Only the short-allele counts depend on this value.

    Returns:
        Per-locus allele information ordered by schema-relative FASTA path.

    Raises:
        FileNotFoundError: If ``schema_directory_path`` does not exist.
        NotADirectoryError: If ``schema_directory_path`` is not a directory.
        ValueError: If ``minimum_sequence_length`` is negative.
        OSError: If a FASTA file cannot be read or decompressed.
        UnicodeError: If a FASTA file is not valid UTF-8 text.
    """
    if minimum_sequence_length < 0:
        raise ValueError("minimum_sequence_length cannot be negative")

    resolved_schema_directory_path = _resolve_schema_directory_path(
        schema_directory_path
    )
    cached_analysis = _get_or_create_cached_analysis(
        resolved_schema_directory_path
    )
    return tuple(
        LocusAlleleInformation(
            locus_name=cached_locus_information.locus_name,
            total_allele_count=cached_locus_information.total_allele_count,
            valid_allele_count=cached_locus_information.valid_allele_count,
            invalid_allele_count=(
                cached_locus_information.invalid_allele_count
            ),
            validated_allele_proportion=(
                cached_locus_information.valid_allele_count
                / cached_locus_information.total_allele_count
                if cached_locus_information.total_allele_count
                else 0.0
            ),
            distinct_protein_allele_count=(
                cached_locus_information.distinct_protein_allele_count
            ),
            incomplete_orf_allele_count=(
                cached_locus_information.incomplete_orf_allele_count
            ),
            ambiguous_bases_allele_count=(
                cached_locus_information.ambiguous_bases_allele_count
            ),
            missing_start_or_stop_codon_allele_count=(
                cached_locus_information.missing_start_or_stop_codon_allele_count
            ),
            in_frame_stop_codon_allele_count=(
                cached_locus_information.in_frame_stop_codon_allele_count
            ),
            alleles_below_minimum_length_count=sum(
                allele_length < minimum_sequence_length
                for allele_length in cached_locus_information.allele_lengths
            ),
            minimum_sequence_length=minimum_sequence_length,
            alleles_below_mode_threshold_count=(
                cached_locus_information.alleles_below_mode_threshold_count
            ),
            alleles_above_mode_threshold_count=(
                cached_locus_information.alleles_above_mode_threshold_count
            ),
            missing_allele_identifier_count=(
                cached_locus_information.missing_allele_identifier_count
            ),
        )
        for cached_locus_information in cached_analysis
    )


def _parse_non_negative_integer(argument_value: str) -> int:
    """Parse a non-negative integer command-line argument."""
    parsed_integer = int(argument_value)
    if parsed_integer < 0:
        raise argparse.ArgumentTypeError("value must be zero or greater")
    return parsed_integer


def main() -> None:
    """Print per-locus allele information as JSON."""
    argument_parser = argparse.ArgumentParser(
        description="Calculate per-locus allele information for an MLST scheme."
    )
    argument_parser.add_argument(
        "schema_directory",
        type=Path,
        help="directory containing an extracted MLST scheme",
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
    locus_information = get_scheme_allele_information(
        command_line_arguments.schema_directory,
        minimum_sequence_length=command_line_arguments.minimum_sequence_length,
    )
    print(
        json.dumps(
            [asdict(locus_information_item) for locus_information_item in locus_information],
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
