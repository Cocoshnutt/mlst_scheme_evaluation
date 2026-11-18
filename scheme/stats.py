"""Calculate, cache, and visualize per-locus MLST scheme statistics.

This module treats each recognized FASTA file in an extracted scheme as one
locus and summarizes the lengths and nucleotide composition of its allele
sequences. The public library API returns explicit per-locus statistics and
creates the four established allele-length Altair panels plus a GC-content
distribution:
the allele-count distribution, modal-length distribution, minimum/maximum/
mean/median lengths per locus, and an interactive box-plot overview.

The cache stores threshold-independent exact-length frequency distributions.
The notebook can therefore change the descriptive frequent-peak threshold and
reclassify single-mode, multimodal, and outlier loci without rereading allele
sequences. A frequent peak is an exact allele length whose count is at least a
configurable percentage of the dominant exact-length count. This is a project
heuristic, not a literature-derived quality threshold.
"""

import gzip
import json
import math
import tempfile
from collections import Counter
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path

import altair as alt

from .info import (
    FASTA_FILENAME_SUFFIXES,
    NucleotideBaseComposition,
    count_nucleotide_base_composition,
)

SCHEME_STATISTICS_CACHE_FILENAME = ".scheme-statistics-cache.json"
LEGACY_STATISTICS_CACHE_FILENAME = ".schema-statistics-cache.json"
SCHEME_STATISTICS_CACHE_FILENAMES = (
    SCHEME_STATISTICS_CACHE_FILENAME,
    LEGACY_STATISTICS_CACHE_FILENAME,
)
SCHEME_STATISTICS_CACHE_VERSION = 3
DEFAULT_FREQUENT_PEAK_RELATIVE_THRESHOLD_PERCENTAGE = 50.0

type AltairChart = alt.Chart | alt.LayerChart | alt.VConcatChart

__all__ = [
    "DEFAULT_FREQUENT_PEAK_RELATIVE_THRESHOLD_PERCENTAGE",
    "SCHEME_STATISTICS_CACHE_FILENAME",
    "SCHEME_STATISTICS_CACHE_FILENAMES",
    "SCHEME_STATISTICS_CACHE_VERSION",
    "AlleleLengthPeak",
    "LocusStatistics",
    "SchemeStatisticsCharts",
    "create_scheme_statistics_charts",
    "get_scheme_statistics",
]


@dataclass(frozen=True)
class AlleleLengthPeak:
    """One exact allele length meeting the frequent-peak heuristic.

    Attributes:
        allele_length: Exact sequence length in nucleotides.
        allele_count: Number of allele records with this length.
        allele_percentage: Percentage of the locus's allele records with this
            length.
    """

    allele_length: int
    allele_count: int
    allele_percentage: float


@dataclass(frozen=True)
class LocusStatistics:
    """Allele count and sequence-length statistics for one locus.

    Attributes:
        locus_name: Scheme-relative path of the locus FASTA file.
        allele_count: Number of allele sequences in the locus.
        modal_allele_length: Most frequently occurring allele length. The
            shortest length is selected when modes are tied.
        modal_allele_lengths: Every length tied for the highest frequency, in
            ascending order.
        mean_allele_length: Arithmetic mean allele length.
        frequent_allele_length_peaks: Exact lengths whose counts meet the
            selected proportion of the dominant exact-length count.
        alleles_below_mode_threshold_count: Alleles strictly below 80% of the
            shortest tied mode.
        alleles_above_mode_threshold_count: Alleles strictly above 120% of the
            shortest tied mode.
        length_distribution_category: Descriptive category based on the
            number of frequent peaks and whether modal-length outliers exist.
        minimum_allele_length: Length of the shortest allele.
        first_quartile_allele_length: 25th percentile allele length.
        median_allele_length: Median allele length.
        third_quartile_allele_length: 75th percentile allele length.
        maximum_allele_length: Length of the longest allele.
        lower_whisker_allele_length: Smallest length inside the lower Tukey
            fence.
        upper_whisker_allele_length: Largest length inside the upper Tukey
            fence.
        guanine_cytosine_content_percentage: Pooled locus GC percentage using
            canonical bases only, or ``None`` when no canonical bases exist.
        canonical_nucleotide_base_count: Total A/C/G/T bases used as the
            locus GC-content denominator.
        excluded_nucleotide_base_count: Total noncanonical, non-whitespace
            characters excluded from the GC-content denominator.
    """

    locus_name: str
    allele_count: int
    modal_allele_length: int | None
    modal_allele_lengths: tuple[int, ...]
    mean_allele_length: float | None
    frequent_allele_length_peaks: tuple[AlleleLengthPeak, ...]
    alleles_below_mode_threshold_count: int
    alleles_above_mode_threshold_count: int
    length_distribution_category: str
    minimum_allele_length: int | None
    first_quartile_allele_length: float | None
    median_allele_length: float | None
    third_quartile_allele_length: float | None
    maximum_allele_length: int | None
    lower_whisker_allele_length: int | None
    upper_whisker_allele_length: int | None
    guanine_cytosine_content_percentage: float | None
    canonical_nucleotide_base_count: int
    excluded_nucleotide_base_count: int


@dataclass(frozen=True)
class _CachedLocusStatistics:
    """Threshold-independent per-locus values stored in the cache."""

    locus_name: str
    allele_count: int
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


@dataclass(frozen=True)
class SchemeStatisticsCharts:
    """The locus-statistics charts displayed by the Marimo report."""

    total_alleles_chart: AltairChart
    allele_mode_size_chart: AltairChart
    locus_statistics_chart: AltairChart
    allele_size_variation_chart: AltairChart
    guanine_cytosine_content_chart: AltairChart


def _resolve_scheme_directory_path(scheme_directory_path: str | Path) -> Path:
    """Resolve and validate an extracted scheme directory.

    Args:
        scheme_directory_path: Directory containing an extracted MLST scheme.

    Returns:
        Absolute path to the scheme directory.

    Raises:
        FileNotFoundError: If the supplied path does not exist.
        NotADirectoryError: If the supplied path is not a directory.
    """
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


def _is_fasta_file(candidate_file_path: Path) -> bool:
    """Return whether a path has a recognized FASTA filename suffix."""
    lowercase_filename = candidate_file_path.name.lower()
    return candidate_file_path.is_file() and lowercase_filename.endswith(
        FASTA_FILENAME_SUFFIXES
    )


def _find_allele_file_paths(
    scheme_directory_path: Path,
) -> tuple[Path, ...]:
    """Return recognized locus FASTA files in deterministic path order."""
    return tuple(
        sorted(
            candidate_file_path
            for candidate_file_path in scheme_directory_path.rglob("*")
            if _is_fasta_file(candidate_file_path)
        )
    )


def _summarize_fasta_lines(
    fasta_lines: Iterable[str],
) -> tuple[tuple[int, ...], NucleotideBaseComposition]:
    """Summarize lengths and pooled base composition from FASTA text lines."""
    current_sequence_length: int | None = None
    allele_lengths: list[int] = []
    guanine_cytosine_base_count = 0
    canonical_nucleotide_base_count = 0
    excluded_nucleotide_base_count = 0

    for fasta_line in fasta_lines:
        stripped_fasta_line = fasta_line.strip()
        if stripped_fasta_line.startswith(">"):
            if current_sequence_length is not None:
                allele_lengths.append(current_sequence_length)
            current_sequence_length = 0
        elif stripped_fasta_line and current_sequence_length is not None:
            current_sequence_length += len(stripped_fasta_line)
            nucleotide_base_composition = count_nucleotide_base_composition(
                stripped_fasta_line
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

    if current_sequence_length is not None:
        allele_lengths.append(current_sequence_length)

    return (
        tuple(allele_lengths),
        NucleotideBaseComposition(
            guanine_cytosine_base_count=guanine_cytosine_base_count,
            canonical_nucleotide_base_count=canonical_nucleotide_base_count,
            excluded_nucleotide_base_count=excluded_nucleotide_base_count,
        ),
    )


def _read_locus_sequence_summary(
    allele_file_path: Path,
) -> tuple[tuple[int, ...], NucleotideBaseComposition]:
    """Read allele lengths and pooled base composition from one FASTA file."""
    if allele_file_path.name.lower().endswith(".gz"):
        with gzip.open(
            allele_file_path,
            mode="rt",
            encoding="utf-8",
        ) as allele_file:
            return _summarize_fasta_lines(allele_file)

    with allele_file_path.open(encoding="utf-8") as allele_file:
        return _summarize_fasta_lines(allele_file)


def _calculate_percentile(
    sorted_allele_lengths: tuple[int, ...],
    percentile: float,
) -> float:
    """Calculate a linearly interpolated percentile from sorted lengths."""
    percentile_position = (len(sorted_allele_lengths) - 1) * percentile
    lower_position = math.floor(percentile_position)
    upper_position = math.ceil(percentile_position)
    if lower_position == upper_position:
        return float(sorted_allele_lengths[lower_position])

    interpolation_fraction = percentile_position - lower_position
    lower_value = sorted_allele_lengths[lower_position]
    upper_value = sorted_allele_lengths[upper_position]
    return lower_value + (upper_value - lower_value) * interpolation_fraction


def _calculate_locus_statistics(
    locus_name: str,
    allele_lengths: tuple[int, ...],
    nucleotide_base_composition: NucleotideBaseComposition,
) -> _CachedLocusStatistics:
    """Calculate threshold-independent length and GC statistics for one locus."""
    if not allele_lengths:
        return _CachedLocusStatistics(
            locus_name=locus_name,
            allele_count=0,
            modal_allele_length=None,
            modal_allele_lengths=(),
            mean_allele_length=None,
            allele_length_frequencies=(),
            alleles_below_mode_threshold_count=0,
            alleles_above_mode_threshold_count=0,
            minimum_allele_length=None,
            first_quartile_allele_length=None,
            median_allele_length=None,
            third_quartile_allele_length=None,
            maximum_allele_length=None,
            lower_whisker_allele_length=None,
            upper_whisker_allele_length=None,
            guanine_cytosine_base_count=(
                nucleotide_base_composition.guanine_cytosine_base_count
            ),
            canonical_nucleotide_base_count=(
                nucleotide_base_composition.canonical_nucleotide_base_count
            ),
            excluded_nucleotide_base_count=(
                nucleotide_base_composition.excluded_nucleotide_base_count
            ),
        )

    sorted_allele_lengths = tuple(sorted(allele_lengths))
    length_frequencies = Counter(sorted_allele_lengths)
    highest_length_frequency = max(length_frequencies.values())
    modal_allele_lengths = tuple(
        allele_length
        for allele_length, allele_length_frequency in sorted(
            length_frequencies.items()
        )
        if allele_length_frequency == highest_length_frequency
    )
    modal_allele_length = modal_allele_lengths[0]
    lower_mode_threshold = modal_allele_length * 0.8
    upper_mode_threshold = modal_allele_length * 1.2
    alleles_below_mode_threshold_count = sum(
        allele_length < lower_mode_threshold
        for allele_length in sorted_allele_lengths
    )
    alleles_above_mode_threshold_count = sum(
        allele_length > upper_mode_threshold
        for allele_length in sorted_allele_lengths
    )
    first_quartile_allele_length = _calculate_percentile(
        sorted_allele_lengths,
        0.25,
    )
    median_allele_length = _calculate_percentile(sorted_allele_lengths, 0.5)
    third_quartile_allele_length = _calculate_percentile(
        sorted_allele_lengths,
        0.75,
    )
    interquartile_range = (
        third_quartile_allele_length - first_quartile_allele_length
    )
    lower_fence = first_quartile_allele_length - 1.5 * interquartile_range
    upper_fence = third_quartile_allele_length + 1.5 * interquartile_range
    lower_whisker_allele_length = next(
        allele_length
        for allele_length in sorted_allele_lengths
        if allele_length >= lower_fence
    )
    upper_whisker_allele_length = next(
        allele_length
        for allele_length in reversed(sorted_allele_lengths)
        if allele_length <= upper_fence
    )

    return _CachedLocusStatistics(
        locus_name=locus_name,
        allele_count=len(sorted_allele_lengths),
        modal_allele_length=modal_allele_length,
        modal_allele_lengths=modal_allele_lengths,
        mean_allele_length=sum(sorted_allele_lengths) / len(sorted_allele_lengths),
        allele_length_frequencies=tuple(sorted(length_frequencies.items())),
        alleles_below_mode_threshold_count=(
            alleles_below_mode_threshold_count
        ),
        alleles_above_mode_threshold_count=(
            alleles_above_mode_threshold_count
        ),
        minimum_allele_length=sorted_allele_lengths[0],
        first_quartile_allele_length=first_quartile_allele_length,
        median_allele_length=median_allele_length,
        third_quartile_allele_length=third_quartile_allele_length,
        maximum_allele_length=sorted_allele_lengths[-1],
        lower_whisker_allele_length=lower_whisker_allele_length,
        upper_whisker_allele_length=upper_whisker_allele_length,
        guanine_cytosine_base_count=(
            nucleotide_base_composition.guanine_cytosine_base_count
        ),
        canonical_nucleotide_base_count=(
            nucleotide_base_composition.canonical_nucleotide_base_count
        ),
        excluded_nucleotide_base_count=(
            nucleotide_base_composition.excluded_nucleotide_base_count
        ),
    )


def _validate_frequent_peak_relative_threshold_percentage(
    frequent_peak_relative_threshold_percentage: float,
) -> None:
    """Validate the descriptive frequent-peak threshold."""
    if (
        not math.isfinite(frequent_peak_relative_threshold_percentage)
        or not 0 <= frequent_peak_relative_threshold_percentage <= 100
    ):
        raise ValueError(
            "frequent_peak_relative_threshold_percentage must be between "
            "0 and 100"
        )


def _materialize_locus_statistics(
    cached_locus_statistics: _CachedLocusStatistics,
    frequent_peak_relative_threshold_percentage: float,
) -> LocusStatistics:
    """Apply a peak threshold to cached exact-length frequencies."""
    if cached_locus_statistics.allele_count == 0:
        frequent_allele_length_peaks: tuple[AlleleLengthPeak, ...] = ()
        length_distribution_category = "Empty"
    else:
        highest_length_frequency = max(
            allele_count
            for _, allele_count in (
                cached_locus_statistics.allele_length_frequencies
            )
        )
        minimum_peak_count = (
            highest_length_frequency
            * frequent_peak_relative_threshold_percentage
            / 100
        )
        frequent_allele_length_peaks = tuple(
            AlleleLengthPeak(
                allele_length=allele_length,
                allele_count=allele_count,
                allele_percentage=(
                    allele_count / cached_locus_statistics.allele_count * 100
                ),
            )
            for allele_length, allele_count in (
                cached_locus_statistics.allele_length_frequencies
            )
            if allele_count >= minimum_peak_count
        )
        has_length_outliers = (
            cached_locus_statistics.alleles_below_mode_threshold_count
            + cached_locus_statistics.alleles_above_mode_threshold_count
            > 0
        )
        modality_label = (
            "Multimodal"
            if len(frequent_allele_length_peaks) > 1
            else "Single-mode"
        )
        length_distribution_category = (
            f"{modality_label} with outliers"
            if has_length_outliers
            else modality_label
        )

    nucleotide_base_composition = NucleotideBaseComposition(
        guanine_cytosine_base_count=(
            cached_locus_statistics.guanine_cytosine_base_count
        ),
        canonical_nucleotide_base_count=(
            cached_locus_statistics.canonical_nucleotide_base_count
        ),
        excluded_nucleotide_base_count=(
            cached_locus_statistics.excluded_nucleotide_base_count
        ),
    )
    return LocusStatistics(
        locus_name=cached_locus_statistics.locus_name,
        allele_count=cached_locus_statistics.allele_count,
        modal_allele_length=cached_locus_statistics.modal_allele_length,
        modal_allele_lengths=cached_locus_statistics.modal_allele_lengths,
        mean_allele_length=cached_locus_statistics.mean_allele_length,
        frequent_allele_length_peaks=frequent_allele_length_peaks,
        alleles_below_mode_threshold_count=(
            cached_locus_statistics.alleles_below_mode_threshold_count
        ),
        alleles_above_mode_threshold_count=(
            cached_locus_statistics.alleles_above_mode_threshold_count
        ),
        length_distribution_category=length_distribution_category,
        minimum_allele_length=cached_locus_statistics.minimum_allele_length,
        first_quartile_allele_length=(
            cached_locus_statistics.first_quartile_allele_length
        ),
        median_allele_length=cached_locus_statistics.median_allele_length,
        third_quartile_allele_length=(
            cached_locus_statistics.third_quartile_allele_length
        ),
        maximum_allele_length=cached_locus_statistics.maximum_allele_length,
        lower_whisker_allele_length=(
            cached_locus_statistics.lower_whisker_allele_length
        ),
        upper_whisker_allele_length=(
            cached_locus_statistics.upper_whisker_allele_length
        ),
        guanine_cytosine_content_percentage=(
            nucleotide_base_composition.guanine_cytosine_content_percentage
        ),
        canonical_nucleotide_base_count=(
            nucleotide_base_composition.canonical_nucleotide_base_count
        ),
        excluded_nucleotide_base_count=(
            nucleotide_base_composition.excluded_nucleotide_base_count
        ),
    )


def _create_source_signature(
    scheme_directory_path: Path,
    allele_file_paths: Iterable[Path],
) -> tuple[tuple[str, int, int], ...]:
    """Create a metadata signature used to invalidate stale statistics."""
    source_signature: list[tuple[str, int, int]] = []
    for allele_file_path in allele_file_paths:
        allele_file_status = allele_file_path.stat()
        source_signature.append(
            (
                allele_file_path.relative_to(
                    scheme_directory_path
                ).as_posix(),
                allele_file_status.st_size,
                allele_file_status.st_mtime_ns,
            )
        )
    return tuple(source_signature)


def _read_cached_scheme_statistics(
    scheme_directory_path: Path,
    source_signature: tuple[tuple[str, int, int], ...],
) -> tuple[_CachedLocusStatistics, ...] | None:
    """Read a compatible statistics cache, or return ``None``."""
    cache_file_path = next(
        (
            scheme_directory_path / cache_filename
            for cache_filename in SCHEME_STATISTICS_CACHE_FILENAMES
            if (scheme_directory_path / cache_filename).is_file()
        ),
        scheme_directory_path / SCHEME_STATISTICS_CACHE_FILENAME,
    )
    try:
        cached_data = json.loads(cache_file_path.read_text(encoding="utf-8"))
        if cached_data["cache_version"] != SCHEME_STATISTICS_CACHE_VERSION:
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
            _CachedLocusStatistics(
                locus_name=str(locus_data["locus_name"]),
                allele_count=int(locus_data["allele_count"]),
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
            for locus_data in cached_data["locus_statistics"]
        )
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return None


def _optional_integer(value: object) -> int | None:
    """Convert a cached optional value to an integer."""
    return None if value is None else int(value)


def _optional_float(value: object) -> float | None:
    """Convert a cached optional value to a float."""
    return None if value is None else float(value)


def _write_cached_scheme_statistics(
    scheme_directory_path: Path,
    source_signature: tuple[tuple[str, int, int], ...],
    locus_statistics: tuple[_CachedLocusStatistics, ...],
) -> None:
    """Atomically write calculated locus statistics to the scheme cache."""
    cached_data = {
        "cache_version": SCHEME_STATISTICS_CACHE_VERSION,
        "source_signature": source_signature,
        "locus_statistics": [
            asdict(locus_information) for locus_information in locus_statistics
        ],
    }
    temporary_cache_file_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=scheme_directory_path,
            prefix=f"{SCHEME_STATISTICS_CACHE_FILENAME}.",
            delete=False,
        ) as temporary_cache_file:
            json.dump(cached_data, temporary_cache_file, separators=(",", ":"))
            temporary_cache_file_path = Path(temporary_cache_file.name)

        temporary_cache_file_path.replace(
            scheme_directory_path / SCHEME_STATISTICS_CACHE_FILENAME
        )
    finally:
        if temporary_cache_file_path is not None:
            temporary_cache_file_path.unlink(missing_ok=True)


def get_scheme_statistics(
    scheme_directory_path: str | Path,
    *,
    frequent_peak_relative_threshold_percentage: float = (
        DEFAULT_FREQUENT_PEAK_RELATIVE_THRESHOLD_PERCENTAGE
    ),
) -> tuple[LocusStatistics, ...]:
    """Return cached or newly calculated statistics for every locus.

    Args:
        scheme_directory_path: Directory containing an extracted MLST scheme.
        frequent_peak_relative_threshold_percentage: Minimum exact-length
            count as a percentage of the dominant exact-length count. This
            descriptive threshold controls frequent peaks and the modality
            category without invalidating the underlying cache.

    Returns:
        Locus statistics ordered by scheme-relative FASTA path.

    Raises:
        FileNotFoundError: If ``scheme_directory_path`` does not exist.
        NotADirectoryError: If ``scheme_directory_path`` is not a directory.
        OSError: If a locus file cannot be read or decompressed.
        UnicodeError: If a locus file is not valid UTF-8 text.
        ValueError: If the frequent-peak threshold is outside 0–100.
    """
    _validate_frequent_peak_relative_threshold_percentage(
        frequent_peak_relative_threshold_percentage
    )
    resolved_scheme_directory_path = _resolve_scheme_directory_path(
        scheme_directory_path
    )
    allele_file_paths = _find_allele_file_paths(resolved_scheme_directory_path)
    source_signature = _create_source_signature(
        resolved_scheme_directory_path,
        allele_file_paths,
    )
    cached_scheme_statistics = _read_cached_scheme_statistics(
        resolved_scheme_directory_path,
        source_signature,
    )
    if cached_scheme_statistics is not None:
        return tuple(
            _materialize_locus_statistics(
                cached_locus_statistics,
                frequent_peak_relative_threshold_percentage,
            )
            for cached_locus_statistics in cached_scheme_statistics
        )

    from .static_analysis import get_shared_scheme_static_analysis

    shared_static_analysis = get_shared_scheme_static_analysis(
        resolved_scheme_directory_path
    )
    cached_locus_statistics = tuple(
        _CachedLocusStatistics(
            locus_name=locus_analysis.relative_file_path,
            allele_count=locus_analysis.allele_count,
            modal_allele_length=locus_analysis.modal_allele_length,
            modal_allele_lengths=locus_analysis.modal_allele_lengths,
            mean_allele_length=locus_analysis.mean_allele_length,
            allele_length_frequencies=(
                locus_analysis.allele_length_frequencies
            ),
            alleles_below_mode_threshold_count=(
                locus_analysis.alleles_below_mode_threshold_count
            ),
            alleles_above_mode_threshold_count=(
                locus_analysis.alleles_above_mode_threshold_count
            ),
            minimum_allele_length=locus_analysis.minimum_allele_length,
            first_quartile_allele_length=(
                locus_analysis.first_quartile_allele_length
            ),
            median_allele_length=locus_analysis.median_allele_length,
            third_quartile_allele_length=(
                locus_analysis.third_quartile_allele_length
            ),
            maximum_allele_length=locus_analysis.maximum_allele_length,
            lower_whisker_allele_length=(
                locus_analysis.lower_whisker_allele_length
            ),
            upper_whisker_allele_length=(
                locus_analysis.upper_whisker_allele_length
            ),
            guanine_cytosine_base_count=(
                locus_analysis.guanine_cytosine_base_count
            ),
            canonical_nucleotide_base_count=(
                locus_analysis.canonical_nucleotide_base_count
            ),
            excluded_nucleotide_base_count=(
                locus_analysis.excluded_nucleotide_base_count
            ),
        )
        for locus_analysis in shared_static_analysis.loci
    )
    try:
        _write_cached_scheme_statistics(
            resolved_scheme_directory_path,
            source_signature,
            cached_locus_statistics,
        )
    except OSError:
        # Read-only scheme directories can still be analyzed without caching.
        pass
    return tuple(
        _materialize_locus_statistics(
            locus_statistics,
            frequent_peak_relative_threshold_percentage,
        )
        for locus_statistics in cached_locus_statistics
    )


def _apply_chart_style(chart: AltairChart) -> AltairChart:
    """Apply shared sizing and typography to a statistics chart."""
    return (
        chart.configure_view(stroke=None)
        .configure_axis(labelFont="Inter", titleFont="Inter")
        .configure_title(font="Inter", anchor="start", fontSize=18)
        .configure_legend(labelFont="Inter", titleFont="Inter")
    )


def create_scheme_statistics_charts(
    locus_statistics: Iterable[LocusStatistics],
) -> SchemeStatisticsCharts:
    """Create the four Altair locus-statistics panels.

    Args:
        locus_statistics: Per-locus statistics returned by
            :func:`get_scheme_statistics`.

    Returns:
        Styled Altair charts corresponding to panels A through D and the
        per-locus GC-content distribution.
    """
    ordered_locus_statistics = tuple(locus_statistics)
    # Keep chart datasets deliberately narrow. Embedding the complete public
    # dataclass would copy tied modes and peak records into every chart layer,
    # producing multi-megabyte Marimo outputs for large schemes.
    locus_rows = [
        {
            "locus_index": locus_index,
            "locus_name": locus_information.locus_name,
            "allele_count": locus_information.allele_count,
            "modal_allele_length": locus_information.modal_allele_length,
            "mean_allele_length": locus_information.mean_allele_length,
            "minimum_allele_length": locus_information.minimum_allele_length,
            "first_quartile_allele_length": (
                locus_information.first_quartile_allele_length
            ),
            "median_allele_length": locus_information.median_allele_length,
            "third_quartile_allele_length": (
                locus_information.third_quartile_allele_length
            ),
            "maximum_allele_length": locus_information.maximum_allele_length,
            "lower_whisker_allele_length": (
                locus_information.lower_whisker_allele_length
            ),
            "upper_whisker_allele_length": (
                locus_information.upper_whisker_allele_length
            ),
            "guanine_cytosine_content_percentage": (
                locus_information.guanine_cytosine_content_percentage
            ),
        }
        for locus_index, locus_information in enumerate(
            ordered_locus_statistics,
            start=1,
        )
    ]
    nonempty_locus_rows = [
        locus_row
        for locus_row in locus_rows
        if locus_row["allele_count"] > 0
    ]
    total_allele_rows = [
        {"allele_count": locus_row["allele_count"]}
        for locus_row in locus_rows
    ]
    modal_allele_length_rows = [
        {"modal_allele_length": locus_row["modal_allele_length"]}
        for locus_row in nonempty_locus_rows
    ]
    guanine_cytosine_content_rows = [
        {
            "guanine_cytosine_content_percentage": locus_row[
                "guanine_cytosine_content_percentage"
            ]
        }
        for locus_row in locus_rows
        if locus_row["guanine_cytosine_content_percentage"] is not None
    ]

    total_alleles_chart = _apply_chart_style(
        alt.Chart(alt.Data(values=total_allele_rows))
        .mark_bar(color="#4C78A8")
        .encode(
            x=alt.X(
                "allele_count:Q",
                bin=alt.Bin(maxbins=30),
                title="Alleles per locus",
            ),
            y=alt.Y("count():Q", title="Number of loci"),
            tooltip=[
                alt.Tooltip("allele_count:Q", bin=True, title="Allele range"),
                alt.Tooltip("count():Q", title="Loci"),
            ],
        )
        .properties(title="Panel A · Total Alleles", width=1100, height=280)
    )

    allele_mode_size_chart = _apply_chart_style(
        alt.Chart(alt.Data(values=modal_allele_length_rows))
        .mark_bar(color="#72B7B2")
        .encode(
            x=alt.X(
                "modal_allele_length:Q",
                bin=alt.Bin(maxbins=30),
                title="Modal allele length (bp)",
            ),
            y=alt.Y("count():Q", title="Number of loci"),
            tooltip=[
                alt.Tooltip(
                    "modal_allele_length:Q",
                    bin=True,
                    title="Mode-size range",
                ),
                alt.Tooltip("count():Q", title="Loci"),
            ],
        )
        .properties(title="Panel B · Allele Mode Size", width=1100, height=280)
    )

    locus_length_rows = [
        {
            "locus_index": locus_row["locus_index"],
            "locus_name": locus_row["locus_name"],
            "allele_count": locus_row["allele_count"],
            "length_statistic": length_statistic,
            "allele_length": locus_row[field_name],
        }
        for locus_row in nonempty_locus_rows
        for length_statistic, field_name in (
            ("Minimum", "minimum_allele_length"),
            ("Maximum", "maximum_allele_length"),
            ("Mean", "mean_allele_length"),
            ("Median", "median_allele_length"),
        )
    ]
    locus_statistics_chart = _apply_chart_style(
        alt.Chart(alt.Data(values=locus_length_rows))
        .mark_circle(size=35, opacity=0.8)
        .encode(
            x=alt.X("allele_length:Q", title="Allele size (bp)"),
            y=alt.Y("allele_count:Q", title="Number of alleles"),
            color=alt.Color(
                "length_statistic:N",
                title="Statistic",
                scale=alt.Scale(
                    domain=["Minimum", "Maximum", "Mean", "Median"],
                    range=["#4C78A8", "#F58518", "#B279A2", "#54A24B"],
                ),
            ),
            tooltip=[
                alt.Tooltip("locus_name:N", title="Locus"),
                alt.Tooltip("length_statistic:N", title="Statistic"),
                alt.Tooltip("allele_length:Q", title="Allele size (bp)"),
                alt.Tooltip("allele_count:Q", title="Number of alleles"),
            ],
        )
        .properties(title="Panel C · Locus Statistics", width=1100, height=320)
        .interactive(bind_y=False)
    )

    visible_locus_interval = alt.selection_interval(
        name="visible_locus_interval",
        encodings=["x"],
    )
    box_plot_base = alt.Chart(alt.Data(values=nonempty_locus_rows)).encode(
        x=alt.X(
            "locus_index:Q",
            title="Locus index",
            scale=alt.Scale(domain=visible_locus_interval),
        ),
        tooltip=[
            alt.Tooltip("locus_name:N", title="Locus"),
            alt.Tooltip("allele_count:Q", title="Alleles"),
            alt.Tooltip(
                "minimum_allele_length:Q",
                title="Minimum (bp)",
            ),
            alt.Tooltip(
                "first_quartile_allele_length:Q",
                title="Q1 (bp)",
                format=".1f",
            ),
            alt.Tooltip(
                "median_allele_length:Q",
                title="Median (bp)",
                format=".1f",
            ),
            alt.Tooltip(
                "third_quartile_allele_length:Q",
                title="Q3 (bp)",
                format=".1f",
            ),
            alt.Tooltip(
                "maximum_allele_length:Q",
                title="Maximum (bp)",
            ),
        ],
    )
    whisker_chart = box_plot_base.mark_rule(
        color="#00BFFF",
        strokeWidth=2,
    ).encode(
        y=alt.Y(
            "lower_whisker_allele_length:Q",
            title="Allele length (bp)",
        ),
        y2="upper_whisker_allele_length:Q",
    )
    quartile_box_chart = box_plot_base.mark_bar(
        color="#72B7B2",
        size=10,
    ).encode(
        y="first_quartile_allele_length:Q",
        y2="third_quartile_allele_length:Q",
    )
    median_tick_chart = box_plot_base.mark_tick(
        color="#FFB000",
        size=12,
        thickness=2,
    ).encode(y="median_allele_length:Q")
    box_plot_chart = (
        (whisker_chart + quartile_box_chart + median_tick_chart)
        .transform_filter(visible_locus_interval)
        .properties(width=1100, height=320)
    )
    range_selector_chart = (
        alt.Chart(alt.Data(values=nonempty_locus_rows))
        .mark_area(color="#BFD7EA")
        .encode(
            x=alt.X("locus_index:Q", title="Drag to select visible loci"),
            y=alt.Y("allele_count:Q", title=None, axis=None),
        )
        .add_params(visible_locus_interval)
        .properties(width=900, height=55)
    )
    allele_size_variation_chart = _apply_chart_style(
        alt.vconcat(
            box_plot_chart,
            range_selector_chart,
            title="Panel D · Allele Size Variation",
        ).resolve_scale(x="independent")
    )

    guanine_cytosine_content_chart = _apply_chart_style(
        alt.Chart(alt.Data(values=guanine_cytosine_content_rows))
        .mark_bar(color="#59A14F")
        .encode(
            x=alt.X(
                "guanine_cytosine_content_percentage:Q",
                bin=alt.Bin(maxbins=30),
                scale=alt.Scale(domain=[0, 100]),
                title="GC content (%)",
            ),
            y=alt.Y("count():Q", title="Number of loci"),
            tooltip=[
                alt.Tooltip(
                    "guanine_cytosine_content_percentage:Q",
                    bin=True,
                    title="GC-content range (%)",
                    format=".2f",
                ),
                alt.Tooltip("count():Q", title="Loci"),
            ],
        )
        .properties(
            title="Per-locus GC Content",
            width=1100,
            height=280,
        )
    )

    return SchemeStatisticsCharts(
        total_alleles_chart=total_alleles_chart,
        allele_mode_size_chart=allele_mode_size_chart,
        locus_statistics_chart=locus_statistics_chart,
        allele_size_variation_chart=allele_size_variation_chart,
        guanine_cytosine_content_chart=guanine_cytosine_content_chart,
    )
