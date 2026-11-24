"""Calculate and cache per-locus MLST scheme statistics.

This module treats each recognized FASTA file in an extracted scheme as one
locus and summarizes the lengths and nucleotide composition of its allele
sequences. The public library API returns explicit per-locus statistics;
Altair specifications live separately in :mod:`scheme.charts`.

The cache stores threshold-independent exact-length frequency distributions.
The notebook can therefore change the descriptive frequent-peak threshold and
reclassify single-mode, multimodal, and outlier loci without rereading allele
sequences. A frequent peak is an exact allele length whose count is at least a
configurable percentage of the dominant exact-length count. This is a project
heuristic, not a literature-derived quality threshold.
"""

import json
import math
from collections import Counter
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path

from .cache_io import write_json_atomically
from .fasta import (
    create_fasta_source_signature,
    find_locus_fasta_paths,
    is_fasta_file,
    iterate_fasta_records,
    read_fasta_records,
    resolve_scheme_directory_path,
)
from .sequence import (
    NucleotideBaseComposition,
    count_nucleotide_base_composition,
    normalize_nucleotide_sequence,
)

SCHEME_STATISTICS_CACHE_FILENAME = ".scheme-statistics-cache.json"
LEGACY_STATISTICS_CACHE_FILENAME = ".schema-statistics-cache.json"
SCHEME_STATISTICS_CACHE_FILENAMES = (
    SCHEME_STATISTICS_CACHE_FILENAME,
    LEGACY_STATISTICS_CACHE_FILENAME,
)
SCHEME_STATISTICS_CACHE_VERSION = 3
DEFAULT_FREQUENT_PEAK_RELATIVE_THRESHOLD_PERCENTAGE = 50.0

__all__ = [
    "DEFAULT_FREQUENT_PEAK_RELATIVE_THRESHOLD_PERCENTAGE",
    "SCHEME_STATISTICS_CACHE_FILENAME",
    "SCHEME_STATISTICS_CACHE_FILENAMES",
    "SCHEME_STATISTICS_CACHE_VERSION",
    "AlleleLengthPeak",
    "LocusStatistics",
    "LocusStatisticsMeasurements",
    "calculate_locus_statistics_measurements",
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
class LocusStatisticsMeasurements:
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


_resolve_scheme_directory_path = resolve_scheme_directory_path
_is_fasta_file = is_fasta_file
_find_allele_file_paths = find_locus_fasta_paths


def _summarize_fasta_lines(
    fasta_lines: Iterable[str],
) -> tuple[tuple[int, ...], NucleotideBaseComposition]:
    """Summarize lengths and pooled base composition from FASTA text lines."""
    allele_lengths: list[int] = []
    guanine_cytosine_base_count = 0
    canonical_nucleotide_base_count = 0
    excluded_nucleotide_base_count = 0

    for fasta_record in iterate_fasta_records(fasta_lines):
        normalized_sequence = normalize_nucleotide_sequence(
            fasta_record.sequence
        )
        allele_lengths.append(len(normalized_sequence))
        nucleotide_base_composition = count_nucleotide_base_composition(
            normalized_sequence
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
    records = read_fasta_records(allele_file_path)
    return _summarize_fasta_lines(
        line
        for record in records
        for line in (f">{record.identifier}", record.sequence)
    )


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


def calculate_locus_statistics_measurements(
    locus_name: str,
    allele_lengths: tuple[int, ...],
    nucleotide_base_composition: NucleotideBaseComposition,
) -> LocusStatisticsMeasurements:
    """Calculate threshold-independent length and GC statistics for one locus."""
    if not allele_lengths:
        return LocusStatisticsMeasurements(
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

    return LocusStatisticsMeasurements(
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
    cached_locus_statistics: LocusStatisticsMeasurements,
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


_create_source_signature = create_fasta_source_signature


def _read_cached_scheme_statistics(
    scheme_directory_path: Path,
    source_signature: tuple[tuple[str, int, int], ...],
) -> tuple[LocusStatisticsMeasurements, ...] | None:
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
            LocusStatisticsMeasurements(
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
    locus_statistics: tuple[LocusStatisticsMeasurements, ...],
) -> None:
    """Atomically write calculated locus statistics to the scheme cache."""
    cached_data = {
        "cache_version": SCHEME_STATISTICS_CACHE_VERSION,
        "source_signature": source_signature,
        "locus_statistics": [
            asdict(locus_information) for locus_information in locus_statistics
        ],
    }
    write_json_atomically(
        scheme_directory_path / SCHEME_STATISTICS_CACHE_FILENAME,
        cached_data,
    )


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
        LocusStatisticsMeasurements(
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
