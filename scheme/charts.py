"""Altair chart construction for intrinsic scheme-evaluation results.

Scientific measurements remain in their owning analysis modules. This module
only projects typed results into bounded chart datasets and specifications for
the Marimo report.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import TYPE_CHECKING

import altair as alt

if TYPE_CHECKING:
    from .similarity import SchemeSimilarityAnalysis
    from .stats import LocusStatistics

type AltairChart = alt.Chart | alt.LayerChart | alt.VConcatChart

__all__ = [
    "SchemeStatisticsCharts",
    "SimilarityScatterChartResult",
    "create_scheme_statistics_charts",
    "create_similarity_scatter_chart",
]


@dataclass(frozen=True)
class SchemeStatisticsCharts:
    """The locus-statistics charts displayed by the Marimo report."""

    total_alleles_chart: AltairChart
    allele_mode_size_chart: AltairChart
    locus_statistics_chart: AltairChart
    allele_size_variation_chart: AltairChart
    guanine_cytosine_content_chart: AltairChart


@dataclass(frozen=True)
class SimilarityScatterChartResult:
    """A bounded chart plus complete and displayed directed-hit counts."""

    chart: AltairChart | None
    displayed_point_count: int
    total_point_count: int


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
            :func:`scheme.stats.get_scheme_statistics`.

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
        locus_row for locus_row in locus_rows if locus_row["allele_count"] > 0
    ]
    total_allele_rows = [
        {"allele_count": locus_row["allele_count"]} for locus_row in locus_rows
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


def create_similarity_scatter_chart(
    scheme_similarity_analysis: SchemeSimilarityAnalysis,
    *,
    maximum_plot_point_count: int = 20_000,
) -> SimilarityScatterChartResult:
    """Create the bounded identity-versus-coverage similarity chart.

    The complete directed-hit count remains available on the analysis result.
    Only the embedded chart dataset is deterministically sampled when it
    exceeds ``maximum_plot_point_count``.

    Args:
        scheme_similarity_analysis: Cached cross-locus BLAST measurements.
        maximum_plot_point_count: Maximum rows embedded in the chart.

    Returns:
        A chart result that retains complete and displayed point counts. Its
        chart is ``None`` when BLAST reported no cross-locus matches.

    Raises:
        ValueError: If ``maximum_plot_point_count`` is less than one.
    """
    if maximum_plot_point_count < 1:
        raise ValueError("maximum_plot_point_count must be at least 1")
    all_plot_rows = [
        {
            "Sequence Level": (
                "Nucleotide"
                if similarity_match.sequence_level == "nucleotide"
                else "Amino acid"
            ),
            "Query Locus": similarity_match.query_locus_name,
            "Subject Locus": similarity_match.subject_locus_name,
            "Identity (%)": round(similarity_match.identity_percentage, 3),
            "Minimum Coverage (%)": round(
                similarity_match.minimum_coverage_percentage,
                3,
            ),
            "Directional BSR": round(
                similarity_match.blast_score_ratio,
                6,
            ),
        }
        for similarity_match in scheme_similarity_analysis.matches
    ]
    if not all_plot_rows:
        return SimilarityScatterChartResult(
            chart=None,
            displayed_point_count=0,
            total_point_count=0,
        )
    if len(all_plot_rows) > maximum_plot_point_count:
        sampling_interval = len(all_plot_rows) / maximum_plot_point_count
        plot_rows = [
            all_plot_rows[int(plot_index * sampling_interval)]
            for plot_index in range(maximum_plot_point_count)
        ]
    else:
        plot_rows = all_plot_rows
    chart = (
        alt.Chart(alt.Data(values=plot_rows))
        .mark_circle(opacity=0.65, size=65)
        .encode(
            x=alt.X(
                "Identity (%):Q",
                title="Identity (%)",
                scale=alt.Scale(domain=[0, 100]),
            ),
            y=alt.Y(
                "Minimum Coverage (%):Q",
                title="Minimum query/subject coverage (%)",
                scale=alt.Scale(domain=[0, 100]),
            ),
            color=alt.Color(
                "Sequence Level:N",
                title="Sequence level",
                scale=alt.Scale(
                    domain=["Nucleotide", "Amino acid"],
                    range=["#3b82f6", "#f97316"],
                ),
            ),
            tooltip=[
                alt.Tooltip("Sequence Level:N"),
                alt.Tooltip("Query Locus:N"),
                alt.Tooltip("Subject Locus:N"),
                alt.Tooltip("Identity (%):Q", format=".3f"),
                alt.Tooltip("Minimum Coverage (%):Q", format=".3f"),
                alt.Tooltip("Directional BSR:Q", format=".6f"),
            ],
        )
        .properties(
            title="Identity and coverage of BLAST-reported cross-locus hits",
            width=1100,
            height=420,
        )
    )
    return SimilarityScatterChartResult(
        chart=chart,
        displayed_point_count=len(plot_rows),
        total_point_count=len(all_plot_rows),
    )
