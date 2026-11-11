"""Calculate, cache, and visualize per-locus MLST scheme statistics.

This module treats each recognized FASTA file in an extracted scheme as one
locus and summarizes the lengths of its allele sequences. The public library
API returns explicit per-locus statistics and creates four Altair panels:
the allele-count distribution, modal-length distribution, minimum/maximum/
median lengths per locus, and an interactive box-plot overview.

Calculated statistics are cached in the extracted scheme so rerunning the
statistics section in ``notebook.py`` does not reread allele sequences.
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

from schema_info import FASTA_FILENAME_SUFFIXES

SCHEMA_STATISTICS_CACHE_FILENAME = ".schema-statistics-cache.json"
SCHEMA_STATISTICS_CACHE_VERSION = 1

type AltairChart = alt.Chart | alt.LayerChart | alt.VConcatChart

__all__ = [
    "SCHEMA_STATISTICS_CACHE_FILENAME",
    "SCHEMA_STATISTICS_CACHE_VERSION",
    "LocusStatistics",
    "SchemaStatisticsCharts",
    "create_schema_statistics_charts",
    "get_schema_statistics",
]


@dataclass(frozen=True)
class LocusStatistics:
    """Allele count and sequence-length statistics for one locus.

    Attributes:
        locus_name: Schema-relative path of the locus FASTA file.
        allele_count: Number of allele sequences in the locus.
        modal_allele_length: Most frequently occurring allele length. The
            shortest length is selected when modes are tied.
        minimum_allele_length: Length of the shortest allele.
        first_quartile_allele_length: 25th percentile allele length.
        median_allele_length: Median allele length.
        third_quartile_allele_length: 75th percentile allele length.
        maximum_allele_length: Length of the longest allele.
        lower_whisker_allele_length: Smallest length inside the lower Tukey
            fence.
        upper_whisker_allele_length: Largest length inside the upper Tukey
            fence.
    """

    locus_name: str
    allele_count: int
    modal_allele_length: int | None
    minimum_allele_length: int | None
    first_quartile_allele_length: float | None
    median_allele_length: float | None
    third_quartile_allele_length: float | None
    maximum_allele_length: int | None
    lower_whisker_allele_length: int | None
    upper_whisker_allele_length: int | None


@dataclass(frozen=True)
class SchemaStatisticsCharts:
    """The four locus-statistics charts displayed by the Marimo report."""

    total_alleles_chart: AltairChart
    allele_mode_size_chart: AltairChart
    locus_statistics_chart: AltairChart
    allele_size_variation_chart: AltairChart


def _resolve_schema_directory_path(schema_directory_path: str | Path) -> Path:
    """Resolve and validate an extracted schema directory.

    Args:
        schema_directory_path: Directory containing an extracted MLST scheme.

    Returns:
        Absolute path to the schema directory.

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


def _is_fasta_file(candidate_file_path: Path) -> bool:
    """Return whether a path has a recognized FASTA filename suffix."""
    lowercase_filename = candidate_file_path.name.lower()
    return candidate_file_path.is_file() and lowercase_filename.endswith(
        FASTA_FILENAME_SUFFIXES
    )


def _find_allele_file_paths(
    schema_directory_path: Path,
) -> tuple[Path, ...]:
    """Return recognized locus FASTA files in deterministic path order."""
    return tuple(
        sorted(
            candidate_file_path
            for candidate_file_path in schema_directory_path.rglob("*")
            if _is_fasta_file(candidate_file_path)
        )
    )


def _iterate_fasta_sequence_lengths(fasta_lines: Iterable[str]) -> Iterable[int]:
    """Yield assembled sequence lengths from FASTA-formatted text lines."""
    current_sequence_length: int | None = None

    for fasta_line in fasta_lines:
        stripped_fasta_line = fasta_line.strip()
        if stripped_fasta_line.startswith(">"):
            if current_sequence_length is not None:
                yield current_sequence_length
            current_sequence_length = 0
        elif stripped_fasta_line and current_sequence_length is not None:
            current_sequence_length += len(stripped_fasta_line)

    if current_sequence_length is not None:
        yield current_sequence_length


def _read_allele_lengths(allele_file_path: Path) -> tuple[int, ...]:
    """Read every assembled allele length from one FASTA file."""
    if allele_file_path.name.lower().endswith(".gz"):
        with gzip.open(
            allele_file_path,
            mode="rt",
            encoding="utf-8",
        ) as allele_file:
            return tuple(_iterate_fasta_sequence_lengths(allele_file))

    with allele_file_path.open(encoding="utf-8") as allele_file:
        return tuple(_iterate_fasta_sequence_lengths(allele_file))


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
) -> LocusStatistics:
    """Calculate allele-count and length statistics for one locus."""
    if not allele_lengths:
        return LocusStatistics(
            locus_name=locus_name,
            allele_count=0,
            modal_allele_length=None,
            minimum_allele_length=None,
            first_quartile_allele_length=None,
            median_allele_length=None,
            third_quartile_allele_length=None,
            maximum_allele_length=None,
            lower_whisker_allele_length=None,
            upper_whisker_allele_length=None,
        )

    sorted_allele_lengths = tuple(sorted(allele_lengths))
    length_frequencies = Counter(sorted_allele_lengths)
    highest_length_frequency = max(length_frequencies.values())
    modal_allele_length = min(
        allele_length
        for allele_length, allele_length_frequency in length_frequencies.items()
        if allele_length_frequency == highest_length_frequency
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

    return LocusStatistics(
        locus_name=locus_name,
        allele_count=len(sorted_allele_lengths),
        modal_allele_length=modal_allele_length,
        minimum_allele_length=sorted_allele_lengths[0],
        first_quartile_allele_length=first_quartile_allele_length,
        median_allele_length=median_allele_length,
        third_quartile_allele_length=third_quartile_allele_length,
        maximum_allele_length=sorted_allele_lengths[-1],
        lower_whisker_allele_length=lower_whisker_allele_length,
        upper_whisker_allele_length=upper_whisker_allele_length,
    )


def _create_source_signature(
    schema_directory_path: Path,
    allele_file_paths: Iterable[Path],
) -> tuple[tuple[str, int, int], ...]:
    """Create a metadata signature used to invalidate stale statistics."""
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


def _read_cached_schema_statistics(
    schema_directory_path: Path,
    source_signature: tuple[tuple[str, int, int], ...],
) -> tuple[LocusStatistics, ...] | None:
    """Read a compatible statistics cache, or return ``None``."""
    cache_file_path = schema_directory_path / SCHEMA_STATISTICS_CACHE_FILENAME
    try:
        cached_data = json.loads(cache_file_path.read_text(encoding="utf-8"))
        if cached_data["cache_version"] != SCHEMA_STATISTICS_CACHE_VERSION:
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
            LocusStatistics(
                locus_name=str(locus_data["locus_name"]),
                allele_count=int(locus_data["allele_count"]),
                modal_allele_length=_optional_integer(
                    locus_data["modal_allele_length"]
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


def _write_cached_schema_statistics(
    schema_directory_path: Path,
    source_signature: tuple[tuple[str, int, int], ...],
    locus_statistics: tuple[LocusStatistics, ...],
) -> None:
    """Atomically write calculated locus statistics to the scheme cache."""
    cached_data = {
        "cache_version": SCHEMA_STATISTICS_CACHE_VERSION,
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
            dir=schema_directory_path,
            prefix=f"{SCHEMA_STATISTICS_CACHE_FILENAME}.",
            delete=False,
        ) as temporary_cache_file:
            json.dump(cached_data, temporary_cache_file, separators=(",", ":"))
            temporary_cache_file_path = Path(temporary_cache_file.name)

        temporary_cache_file_path.replace(
            schema_directory_path / SCHEMA_STATISTICS_CACHE_FILENAME
        )
    finally:
        if temporary_cache_file_path is not None:
            temporary_cache_file_path.unlink(missing_ok=True)


def get_schema_statistics(
    schema_directory_path: str | Path,
) -> tuple[LocusStatistics, ...]:
    """Return cached or newly calculated statistics for every locus.

    Args:
        schema_directory_path: Directory containing an extracted MLST scheme.

    Returns:
        Locus statistics ordered by schema-relative FASTA path.

    Raises:
        FileNotFoundError: If ``schema_directory_path`` does not exist.
        NotADirectoryError: If ``schema_directory_path`` is not a directory.
        OSError: If a locus file cannot be read or decompressed.
        UnicodeError: If a locus file is not valid UTF-8 text.
    """
    resolved_schema_directory_path = _resolve_schema_directory_path(
        schema_directory_path
    )
    allele_file_paths = _find_allele_file_paths(resolved_schema_directory_path)
    source_signature = _create_source_signature(
        resolved_schema_directory_path,
        allele_file_paths,
    )
    cached_schema_statistics = _read_cached_schema_statistics(
        resolved_schema_directory_path,
        source_signature,
    )
    if cached_schema_statistics is not None:
        return cached_schema_statistics

    locus_statistics = tuple(
        _calculate_locus_statistics(
            allele_file_path.relative_to(
                resolved_schema_directory_path
            ).as_posix(),
            _read_allele_lengths(allele_file_path),
        )
        for allele_file_path in allele_file_paths
    )
    try:
        _write_cached_schema_statistics(
            resolved_schema_directory_path,
            source_signature,
            locus_statistics,
        )
    except OSError:
        # Read-only schema directories can still be analyzed without caching.
        pass
    return locus_statistics


def _apply_chart_style(chart: AltairChart) -> AltairChart:
    """Apply shared sizing and typography to a statistics chart."""
    return (
        chart.configure_view(stroke=None)
        .configure_axis(labelFont="Inter", titleFont="Inter")
        .configure_title(font="Inter", anchor="start", fontSize=18)
        .configure_legend(labelFont="Inter", titleFont="Inter")
    )


def create_schema_statistics_charts(
    locus_statistics: Iterable[LocusStatistics],
) -> SchemaStatisticsCharts:
    """Create the four Altair locus-statistics panels.

    Args:
        locus_statistics: Per-locus statistics returned by
            :func:`get_schema_statistics`.

    Returns:
        Four styled Altair charts corresponding to panels A through D.
    """
    ordered_locus_statistics = tuple(locus_statistics)
    locus_rows = [
        {
            "locus_index": locus_index,
            **asdict(locus_information),
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

    total_alleles_chart = _apply_chart_style(
        alt.Chart(alt.Data(values=locus_rows))
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
        alt.Chart(alt.Data(values=nonempty_locus_rows))
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
                    domain=["Minimum", "Maximum", "Median"],
                    range=["#4C78A8", "#F58518", "#54A24B"],
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

    return SchemaStatisticsCharts(
        total_alleles_chart=total_alleles_chart,
        allele_mode_size_chart=allele_mode_size_chart,
        locus_statistics_chart=locus_statistics_chart,
        allele_size_variation_chart=allele_size_variation_chart,
    )
