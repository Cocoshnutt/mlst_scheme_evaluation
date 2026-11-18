"""Calculate and cache the shared static FASTA analysis in one traversal.

Scheme information, locus statistics, and allele information require the same
FASTA records, sequence lengths, base composition, and coding-sequence
validation outcomes.  This module calculates that threshold-independent union
once and stores it in a versioned, source-signature-aware cache.  The three
established public modules project their existing result and compatibility
cache types from this analysis, preserving their APIs and scientific rules.

The cache contains measurements only.  Minimum-length and frequent-peak
thresholds remain outside it and are applied by the consuming modules.
"""

import json
import tempfile
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path

from .allele_info import (
    _evaluate_allele_sequence,
    _find_allele_file_paths,
    _read_fasta_records,
    _remove_fasta_filename_suffix,
    _translate_valid_allele,
)
from .info import NucleotideBaseComposition, count_nucleotide_base_composition
from .stats import _calculate_locus_statistics

SCHEME_STATIC_ANALYSIS_CACHE_FILENAME = ".scheme-static-analysis-cache.json"
SCHEME_STATIC_ANALYSIS_CACHE_FILENAMES = (
    SCHEME_STATIC_ANALYSIS_CACHE_FILENAME,
)
SCHEME_STATIC_ANALYSIS_CACHE_VERSION = 1

__all__ = [
    "SCHEME_STATIC_ANALYSIS_CACHE_FILENAME",
    "SCHEME_STATIC_ANALYSIS_CACHE_FILENAMES",
    "SCHEME_STATIC_ANALYSIS_CACHE_VERSION",
    "SharedLocusStaticAnalysis",
    "SharedSchemeStaticAnalysis",
    "get_shared_scheme_static_analysis",
]

type SourceSignature = tuple[tuple[str, int, int], ...]


@dataclass(frozen=True)
class SharedLocusStaticAnalysis:
    """Threshold-independent measurements for one locus FASTA file."""

    relative_file_path: str
    locus_name: str
    allele_lengths: tuple[int, ...]
    valid_allele_count: int
    distinct_protein_allele_count: int
    incomplete_orf_allele_count: int
    ambiguous_bases_allele_count: int
    missing_start_or_stop_codon_allele_count: int
    in_frame_stop_codon_allele_count: int
    modal_allele_length: int | None
    modal_allele_lengths: tuple[int, ...]
    mean_allele_length: float | None
    allele_length_frequencies: tuple[tuple[int, int], ...]
    alleles_below_mode_threshold_count: int
    alleles_above_mode_threshold_count: int
    minimum_allele_length: int | None
    first_quartile_allele_length: float | None
    median_allele_length: float | None
    third_quartile_allele_length: float | None
    maximum_allele_length: int | None
    lower_whisker_allele_length: int | None
    upper_whisker_allele_length: int | None
    guanine_cytosine_base_count: int
    canonical_nucleotide_base_count: int
    excluded_nucleotide_base_count: int

    @property
    def allele_count(self) -> int:
        """Return the number of FASTA records represented by this locus."""
        return len(self.allele_lengths)


@dataclass(frozen=True)
class SharedSchemeStaticAnalysis:
    """Source signature and ordered shared measurements for one scheme."""

    source_signature: SourceSignature
    loci: tuple[SharedLocusStaticAnalysis, ...]


def _create_source_signature(
    scheme_directory_path: Path,
    allele_file_paths: Iterable[Path],
) -> SourceSignature:
    """Return the path, size, and modification time of every locus file."""
    return tuple(
        (
            allele_file_path.relative_to(scheme_directory_path).as_posix(),
            allele_file_path.stat().st_size,
            allele_file_path.stat().st_mtime_ns,
        )
        for allele_file_path in allele_file_paths
    )


def _resolve_scheme_directory_path(
    scheme_directory_path: str | Path,
) -> Path:
    """Resolve and validate a directory containing scheme FASTA files."""
    resolved_scheme_directory_path = Path(scheme_directory_path).resolve()
    if not resolved_scheme_directory_path.exists():
        raise FileNotFoundError(
            f"Scheme directory does not exist: {resolved_scheme_directory_path}"
        )
    if not resolved_scheme_directory_path.is_dir():
        raise NotADirectoryError(
            f"Scheme path is not a directory: {resolved_scheme_directory_path}"
        )
    return resolved_scheme_directory_path


def _analyze_locus(
    scheme_directory_path: Path,
    allele_file_path: Path,
) -> SharedLocusStaticAnalysis:
    """Calculate the union of all current static metrics for one locus."""
    allele_lengths: list[int] = []
    valid_allele_count = 0
    incomplete_orf_allele_count = 0
    ambiguous_bases_allele_count = 0
    missing_start_or_stop_codon_allele_count = 0
    in_frame_stop_codon_allele_count = 0
    distinct_protein_alleles: set[str] = set()
    guanine_cytosine_base_count = 0
    canonical_nucleotide_base_count = 0
    excluded_nucleotide_base_count = 0

    for _fasta_identifier, allele_sequence in _read_fasta_records(
        allele_file_path
    ):
        normalized_allele_sequence = "".join(allele_sequence.split()).upper()
        allele_lengths.append(len(normalized_allele_sequence))
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
    nucleotide_base_composition = NucleotideBaseComposition(
        guanine_cytosine_base_count=guanine_cytosine_base_count,
        canonical_nucleotide_base_count=canonical_nucleotide_base_count,
        excluded_nucleotide_base_count=excluded_nucleotide_base_count,
    )
    relative_file_path = allele_file_path.relative_to(
        scheme_directory_path
    ).as_posix()
    locus_statistics = _calculate_locus_statistics(
        relative_file_path,
        immutable_allele_lengths,
        nucleotide_base_composition,
    )
    return SharedLocusStaticAnalysis(
        relative_file_path=relative_file_path,
        locus_name=_remove_fasta_filename_suffix(relative_file_path),
        allele_lengths=immutable_allele_lengths,
        valid_allele_count=valid_allele_count,
        distinct_protein_allele_count=len(distinct_protein_alleles),
        incomplete_orf_allele_count=incomplete_orf_allele_count,
        ambiguous_bases_allele_count=ambiguous_bases_allele_count,
        missing_start_or_stop_codon_allele_count=(
            missing_start_or_stop_codon_allele_count
        ),
        in_frame_stop_codon_allele_count=in_frame_stop_codon_allele_count,
        modal_allele_length=locus_statistics.modal_allele_length,
        modal_allele_lengths=locus_statistics.modal_allele_lengths,
        mean_allele_length=locus_statistics.mean_allele_length,
        allele_length_frequencies=locus_statistics.allele_length_frequencies,
        alleles_below_mode_threshold_count=(
            locus_statistics.alleles_below_mode_threshold_count
        ),
        alleles_above_mode_threshold_count=(
            locus_statistics.alleles_above_mode_threshold_count
        ),
        minimum_allele_length=locus_statistics.minimum_allele_length,
        first_quartile_allele_length=(
            locus_statistics.first_quartile_allele_length
        ),
        median_allele_length=locus_statistics.median_allele_length,
        third_quartile_allele_length=(
            locus_statistics.third_quartile_allele_length
        ),
        maximum_allele_length=locus_statistics.maximum_allele_length,
        lower_whisker_allele_length=(
            locus_statistics.lower_whisker_allele_length
        ),
        upper_whisker_allele_length=(
            locus_statistics.upper_whisker_allele_length
        ),
        guanine_cytosine_base_count=guanine_cytosine_base_count,
        canonical_nucleotide_base_count=canonical_nucleotide_base_count,
        excluded_nucleotide_base_count=excluded_nucleotide_base_count,
    )


def _optional_integer(value: object) -> int | None:
    """Convert a cached optional value to an integer."""
    return None if value is None else int(value)


def _optional_float(value: object) -> float | None:
    """Convert a cached optional value to a float."""
    return None if value is None else float(value)


def _read_cached_analysis(
    scheme_directory_path: Path,
    source_signature: SourceSignature,
) -> SharedSchemeStaticAnalysis | None:
    """Read a compatible shared cache, returning ``None`` when unusable."""
    cache_file_path = (
        scheme_directory_path / SCHEME_STATIC_ANALYSIS_CACHE_FILENAME
    )
    try:
        cached_data = json.loads(cache_file_path.read_text(encoding="utf-8"))
        if cached_data["cache_version"] != SCHEME_STATIC_ANALYSIS_CACHE_VERSION:
            return None
        cached_source_signature = tuple(
            (str(relative_path), int(file_size), int(modification_time_ns))
            for relative_path, file_size, modification_time_ns in cached_data[
                "source_signature"
            ]
        )
        if cached_source_signature != source_signature:
            return None
        loci = tuple(
            SharedLocusStaticAnalysis(
                relative_file_path=str(locus_data["relative_file_path"]),
                locus_name=str(locus_data["locus_name"]),
                allele_lengths=tuple(
                    int(allele_length)
                    for allele_length in locus_data["allele_lengths"]
                ),
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
                modal_allele_length=_optional_integer(
                    locus_data["modal_allele_length"]
                ),
                modal_allele_lengths=tuple(
                    int(allele_length)
                    for allele_length in locus_data["modal_allele_lengths"]
                ),
                mean_allele_length=_optional_float(
                    locus_data["mean_allele_length"]
                ),
                allele_length_frequencies=tuple(
                    (int(allele_length), int(allele_count))
                    for allele_length, allele_count in locus_data[
                        "allele_length_frequencies"
                    ]
                ),
                alleles_below_mode_threshold_count=int(
                    locus_data["alleles_below_mode_threshold_count"]
                ),
                alleles_above_mode_threshold_count=int(
                    locus_data["alleles_above_mode_threshold_count"]
                ),
                minimum_allele_length=_optional_integer(
                    locus_data["minimum_allele_length"]
                ),
                first_quartile_allele_length=_optional_float(
                    locus_data["first_quartile_allele_length"]
                ),
                median_allele_length=_optional_float(
                    locus_data["median_allele_length"]
                ),
                third_quartile_allele_length=_optional_float(
                    locus_data["third_quartile_allele_length"]
                ),
                maximum_allele_length=_optional_integer(
                    locus_data["maximum_allele_length"]
                ),
                lower_whisker_allele_length=_optional_integer(
                    locus_data["lower_whisker_allele_length"]
                ),
                upper_whisker_allele_length=_optional_integer(
                    locus_data["upper_whisker_allele_length"]
                ),
                guanine_cytosine_base_count=int(
                    locus_data["guanine_cytosine_base_count"]
                ),
                canonical_nucleotide_base_count=int(
                    locus_data["canonical_nucleotide_base_count"]
                ),
                excluded_nucleotide_base_count=int(
                    locus_data["excluded_nucleotide_base_count"]
                ),
            )
            for locus_data in cached_data["loci"]
        )
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return None

    if tuple(locus.relative_file_path for locus in loci) != tuple(
        relative_path for relative_path, _file_size, _modification_time in source_signature
    ):
        return None
    if any(
        min(
            locus.valid_allele_count,
            locus.distinct_protein_allele_count,
            locus.incomplete_orf_allele_count,
            locus.ambiguous_bases_allele_count,
            locus.missing_start_or_stop_codon_allele_count,
            locus.in_frame_stop_codon_allele_count,
            locus.alleles_below_mode_threshold_count,
            locus.alleles_above_mode_threshold_count,
            locus.guanine_cytosine_base_count,
            locus.canonical_nucleotide_base_count,
            locus.excluded_nucleotide_base_count,
        )
        < 0
        or locus.valid_allele_count > locus.allele_count
        or locus.guanine_cytosine_base_count
        > locus.canonical_nucleotide_base_count
        for locus in loci
    ):
        return None
    return SharedSchemeStaticAnalysis(
        source_signature=cached_source_signature,
        loci=loci,
    )


def _write_cached_analysis(
    scheme_directory_path: Path,
    analysis: SharedSchemeStaticAnalysis,
) -> None:
    """Atomically persist a shared static analysis."""
    cached_data = {
        "cache_version": SCHEME_STATIC_ANALYSIS_CACHE_VERSION,
        "source_signature": analysis.source_signature,
        "loci": [asdict(locus) for locus in analysis.loci],
    }
    temporary_cache_file_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=scheme_directory_path,
            prefix=f"{SCHEME_STATIC_ANALYSIS_CACHE_FILENAME}.",
            suffix=".tmp",
            delete=False,
        ) as temporary_cache_file:
            json.dump(cached_data, temporary_cache_file, separators=(",", ":"))
            temporary_cache_file_path = Path(temporary_cache_file.name)
        temporary_cache_file_path.replace(
            scheme_directory_path / SCHEME_STATIC_ANALYSIS_CACHE_FILENAME
        )
    finally:
        if temporary_cache_file_path is not None:
            temporary_cache_file_path.unlink(missing_ok=True)


def get_shared_scheme_static_analysis(
    scheme_directory_path: str | Path,
) -> SharedSchemeStaticAnalysis:
    """Return the cached or newly calculated one-pass static analysis.

    Args:
        scheme_directory_path: Directory containing scheme locus FASTA files.

    Returns:
        Source provenance and ordered threshold-independent locus measurements.

    Raises:
        FileNotFoundError: If the supplied path does not exist.
        NotADirectoryError: If the supplied path is not a directory.
        OSError: If a FASTA file cannot be read or decompressed.
        UnicodeError: If a FASTA file is not valid UTF-8 text.

    Cache-write failures are non-fatal so readable, read-only schemes remain
    analyzable.
    """
    resolved_scheme_directory_path = _resolve_scheme_directory_path(
        scheme_directory_path
    )
    allele_file_paths = _find_allele_file_paths(resolved_scheme_directory_path)
    source_signature = _create_source_signature(
        resolved_scheme_directory_path,
        allele_file_paths,
    )
    cached_analysis = _read_cached_analysis(
        resolved_scheme_directory_path,
        source_signature,
    )
    if cached_analysis is not None:
        return cached_analysis

    analysis = SharedSchemeStaticAnalysis(
        source_signature=source_signature,
        loci=tuple(
            _analyze_locus(resolved_scheme_directory_path, allele_file_path)
            for allele_file_path in allele_file_paths
        ),
    )
    try:
        _write_cached_analysis(resolved_scheme_directory_path, analysis)
    except OSError:
        pass
    return analysis
