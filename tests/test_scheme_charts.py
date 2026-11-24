"""Tests for intrinsic-report Altair chart projections."""

import unittest

from scheme.charts import create_similarity_scatter_chart
from scheme.similarity import CrossLocusSimilarityMatch, SchemeSimilarityAnalysis


def _similarity_match(match_index: int) -> CrossLocusSimilarityMatch:
    """Create one directed match with distinguishable plotted values."""
    return CrossLocusSimilarityMatch(
        sequence_level="nucleotide" if match_index % 2 == 0 else "amino_acid",
        query_locus_name=f"query_{match_index}",
        subject_locus_name=f"subject_{match_index}",
        query_allele_identifier=f"query_allele_{match_index}",
        subject_allele_identifier=f"subject_allele_{match_index}",
        identity_percentage=80.0 + match_index,
        alignment_length=100 + match_index,
        query_coverage_percentage=70.0 + match_index,
        subject_coverage_percentage=60.0 + match_index,
        minimum_coverage_percentage=60.0 + match_index,
        raw_score=50.0 + match_index,
        blast_score_ratio=0.5 + match_index / 100,
        nucleotide_subject_strand="plus" if match_index % 2 == 0 else None,
    )


def _similarity_analysis(
    matches: tuple[CrossLocusSimilarityMatch, ...],
) -> SchemeSimilarityAnalysis:
    """Create a minimal immutable similarity result for chart tests."""
    return SchemeSimilarityAnalysis(
        locus_count=10,
        nucleotide_selected_locus_count=10,
        amino_acid_selected_locus_count=10,
        nucleotide_possible_directed_pair_count=90,
        amino_acid_possible_directed_pair_count=90,
        matches=matches,
        blastn_version="blastn: test",
        blastp_version="blastp: test",
        method_identifier="test",
        blast_evalue=0.001,
        reused_from_cache=False,
    )


class SimilarityScatterChartTest(unittest.TestCase):
    """Protect empty, bounded, and deterministic chart projection behavior."""

    def test_empty_analysis_returns_explicit_empty_chart_result(self) -> None:
        """Let the notebook render its established no-hit callout."""
        chart_result = create_similarity_scatter_chart(_similarity_analysis(()))

        self.assertIsNone(chart_result.chart)
        self.assertEqual(chart_result.displayed_point_count, 0)
        self.assertEqual(chart_result.total_point_count, 0)

    def test_large_result_is_deterministically_bounded(self) -> None:
        """Embed an evenly spaced preview while reporting the complete count."""
        analysis = _similarity_analysis(
            tuple(_similarity_match(match_index) for match_index in range(5))
        )

        chart_result = create_similarity_scatter_chart(
            analysis,
            maximum_plot_point_count=2,
        )

        self.assertIsNotNone(chart_result.chart)
        self.assertEqual(chart_result.displayed_point_count, 2)
        self.assertEqual(chart_result.total_point_count, 5)
        chart_specification = chart_result.chart.to_dict()
        plotted_rows = chart_specification["data"]["values"]
        self.assertEqual(
            [row["Query Locus"] for row in plotted_rows],
            ["query_0", "query_2"],
        )
        self.assertEqual(
            set(plotted_rows[0]),
            {
                "Sequence Level",
                "Query Locus",
                "Subject Locus",
                "Identity (%)",
                "Minimum Coverage (%)",
                "Directional BSR",
            },
        )

    def test_nonpositive_plot_limit_is_rejected(self) -> None:
        """Avoid a division-by-zero or silently empty bounded chart."""
        with self.assertRaisesRegex(
            ValueError,
            "maximum_plot_point_count must be at least 1",
        ):
            create_similarity_scatter_chart(
                _similarity_analysis((_similarity_match(0),)),
                maximum_plot_point_count=0,
            )


if __name__ == "__main__":
    unittest.main()
