"""Tests for cached per-locus allele-length statistics."""

import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scheme import stats
from scheme.charts import create_scheme_statistics_charts
from scheme.stats import (
    SCHEME_STATISTICS_CACHE_FILENAME,
    SCHEME_STATISTICS_CACHE_FILENAMES,
    SCHEME_STATISTICS_CACHE_VERSION,
    get_scheme_statistics,
)


def _write_length_fixture(
    fasta_file_path: Path,
    allele_lengths: tuple[int, ...],
) -> None:
    """Write deterministic synthetic sequences with the requested lengths."""
    fasta_file_path.write_text(
        "".join(
            f">allele_{allele_index}\n{'A' * allele_length}\n"
            for allele_index, allele_length in enumerate(allele_lengths, start=1)
        ),
        encoding="utf-8",
    )


class LocusLengthStatisticsTest(unittest.TestCase):
    """Test means, modes, peaks, outliers, and categories."""

    def test_empty_locus_has_explicit_empty_category(self) -> None:
        """Represent a zero-record FASTA without fabricated length values."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "empty.fasta").write_text("", encoding="utf-8")

            locus_statistics = get_scheme_statistics(scheme_directory_path)[0]

        self.assertEqual(locus_statistics.allele_count, 0)
        self.assertIsNone(locus_statistics.mean_allele_length)
        self.assertEqual(locus_statistics.modal_allele_lengths, ())
        self.assertEqual(locus_statistics.frequent_allele_length_peaks, ())
        self.assertEqual(locus_statistics.length_distribution_category, "Empty")

    def test_mean_all_tied_modes_and_peak_measurements_are_reported(self) -> None:
        """Keep all global modes while retaining the shortest primary mode."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            _write_length_fixture(
                scheme_directory_path / "tied.fasta",
                (90, 90, 100, 100, 120),
            )

            locus_statistics = get_scheme_statistics(scheme_directory_path)[0]

        self.assertEqual(locus_statistics.modal_allele_length, 90)
        self.assertEqual(locus_statistics.modal_allele_lengths, (90, 100))
        self.assertEqual(locus_statistics.mean_allele_length, 100.0)
        self.assertEqual(
            tuple(
                (
                    peak.allele_length,
                    peak.allele_count,
                    peak.allele_percentage,
                )
                for peak in locus_statistics.frequent_allele_length_peaks
            ),
            ((90, 2, 40.0), (100, 2, 40.0), (120, 1, 20.0)),
        )
        self.assertEqual(locus_statistics.alleles_above_mode_threshold_count, 1)
        self.assertEqual(
            locus_statistics.length_distribution_category,
            "Multimodal with outliers",
        )

    def test_peak_threshold_boundary_is_inclusive_and_reactive(self) -> None:
        """Include exactly 50%-height peaks and remove them above that value."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            _write_length_fixture(
                scheme_directory_path / "peaks.fasta",
                (100, 100, 100, 100, 110, 110),
            )

            default_statistics = get_scheme_statistics(scheme_directory_path)[0]
            stricter_statistics = get_scheme_statistics(
                scheme_directory_path,
                frequent_peak_relative_threshold_percentage=50.1,
            )[0]

        self.assertEqual(
            tuple(
                peak.allele_length
                for peak in default_statistics.frequent_allele_length_peaks
            ),
            (100, 110),
        )
        self.assertEqual(
            default_statistics.length_distribution_category,
            "Multimodal",
        )
        self.assertEqual(
            tuple(
                peak.allele_length
                for peak in stricter_statistics.frequent_allele_length_peaks
            ),
            (100,),
        )
        self.assertEqual(
            stricter_statistics.length_distribution_category,
            "Single-mode",
        )

    def test_outlier_boundaries_are_strict(self) -> None:
        """Do not flag lengths exactly at 80% or 120% of the primary mode."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            _write_length_fixture(
                scheme_directory_path / "boundaries.fasta",
                (80, 100, 100, 120),
            )

            locus_statistics = get_scheme_statistics(
                scheme_directory_path,
                frequent_peak_relative_threshold_percentage=100,
            )[0]

        self.assertEqual(locus_statistics.alleles_below_mode_threshold_count, 0)
        self.assertEqual(locus_statistics.alleles_above_mode_threshold_count, 0)
        self.assertEqual(
            locus_statistics.length_distribution_category,
            "Single-mode",
        )

    def test_single_mode_with_outlier_category(self) -> None:
        """Separate a dominant length with rare outliers from multimodality."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            _write_length_fixture(
                scheme_directory_path / "outlier.fasta",
                (70, 100, 100, 100, 100),
            )

            locus_statistics = get_scheme_statistics(scheme_directory_path)[0]

        self.assertEqual(locus_statistics.alleles_below_mode_threshold_count, 1)
        self.assertEqual(
            locus_statistics.length_distribution_category,
            "Single-mode with outliers",
        )

    def test_panel_c_includes_mean_length_series(self) -> None:
        """Expose the new mean statistic without removing existing series."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            _write_length_fixture(
                scheme_directory_path / "chart.fasta",
                (90, 100, 110),
            )
            locus_statistics = get_scheme_statistics(scheme_directory_path)

            chart_specification = create_scheme_statistics_charts(
                locus_statistics
            ).locus_statistics_chart.to_dict()

        color_domain = chart_specification["encoding"]["color"]["scale"]["domain"]
        self.assertEqual(
            color_domain,
            ["Minimum", "Maximum", "Mean", "Median"],
        )

    def test_locus_gc_content_uses_pooled_canonical_bases(self) -> None:
        """Exclude ambiguous bases and weight the locus value by canonical bases."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "composition.fasta").write_text(
                ">allele_1\nACGTNN\n>allele_2\nGGAA\n",
                encoding="utf-8",
            )

            locus_statistics = get_scheme_statistics(scheme_directory_path)[0]

        self.assertEqual(locus_statistics.canonical_nucleotide_base_count, 8)
        self.assertEqual(locus_statistics.excluded_nucleotide_base_count, 2)
        self.assertEqual(
            locus_statistics.guanine_cytosine_content_percentage,
            50.0,
        )

    def test_locus_without_canonical_bases_has_no_gc_percentage(self) -> None:
        """Keep ambiguous-only loci visible without inventing a zero value."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "ambiguous.fasta").write_text(
                ">allele_1\nNN--\n",
                encoding="utf-8",
            )

            locus_statistics = get_scheme_statistics(scheme_directory_path)[0]

        self.assertIsNone(locus_statistics.guanine_cytosine_content_percentage)
        self.assertEqual(locus_statistics.canonical_nucleotide_base_count, 0)
        self.assertEqual(locus_statistics.excluded_nucleotide_base_count, 4)

    def test_gc_chart_has_canonical_zero_to_one_hundred_domain(self) -> None:
        """Render the per-locus GC distribution on its full percentage domain."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "chart-gc.fasta").write_text(
                ">allele_1\nACGT\n",
                encoding="utf-8",
            )
            locus_statistics = get_scheme_statistics(scheme_directory_path)

            chart_specification = create_scheme_statistics_charts(
                locus_statistics
            ).guanine_cytosine_content_chart.to_dict()

        self.assertEqual(
            chart_specification["encoding"]["x"]["scale"]["domain"],
            [0, 100],
        )

    def test_charts_do_not_embed_unplotted_peak_or_mode_collections(self) -> None:
        """Keep large-scheme chart payloads limited to plotted fields."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            _write_length_fixture(
                scheme_directory_path / "payload.fasta",
                tuple(range(1, 101)),
            )
            locus_statistics = get_scheme_statistics(scheme_directory_path)

            scheme_charts = create_scheme_statistics_charts(locus_statistics)

        for chart in (
            scheme_charts.total_alleles_chart,
            scheme_charts.allele_mode_size_chart,
            scheme_charts.locus_statistics_chart,
            scheme_charts.allele_size_variation_chart,
            scheme_charts.guanine_cytosine_content_chart,
        ):
            serialized_chart = chart.to_json()
            self.assertNotIn("frequent_allele_length_peaks", serialized_chart)
            self.assertNotIn("modal_allele_lengths", serialized_chart)

    def test_invalid_peak_threshold_is_rejected(self) -> None:
        """Reject percentages outside their meaningful domain."""
        with (
            tempfile.TemporaryDirectory() as temporary_directory,
            self.assertRaises(ValueError),
        ):
            get_scheme_statistics(
                temporary_directory,
                frequent_peak_relative_threshold_percentage=100.1,
            )


class SchemeStatisticsCacheTest(unittest.TestCase):
    """Test gzip support and threshold-independent cache reuse."""

    def test_gzip_fasta_lengths_are_analyzed(self) -> None:
        """Read compressed loci with the same length semantics as plain FASTA."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            with gzip.open(
                scheme_directory_path / "compressed.fasta.gz",
                mode="wt",
                encoding="utf-8",
            ) as compressed_fasta_file:
                compressed_fasta_file.write(
                    ">allele_1\nAAA\n>allele_2\nAAAAAA\n"
                )

            locus_statistics = get_scheme_statistics(scheme_directory_path)[0]

        self.assertEqual(locus_statistics.mean_allele_length, 4.5)
        self.assertEqual(locus_statistics.modal_allele_lengths, (3, 6))
        self.assertEqual(locus_statistics.canonical_nucleotide_base_count, 9)
        self.assertEqual(locus_statistics.guanine_cytosine_content_percentage, 0.0)

    def test_threshold_change_reuses_cached_length_frequencies(self) -> None:
        """Reclassify cached histograms without reading FASTA sequences again."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            _write_length_fixture(
                scheme_directory_path / "cached.fasta",
                (100, 100, 100, 100, 110, 110),
            )
            get_scheme_statistics(scheme_directory_path)
            legacy_cache_file_path = (
                scheme_directory_path / SCHEME_STATISTICS_CACHE_FILENAMES[1]
            )
            (scheme_directory_path / SCHEME_STATISTICS_CACHE_FILENAME).replace(
                legacy_cache_file_path
            )

            with patch.object(
                stats,
                "_read_locus_sequence_summary",
                side_effect=AssertionError("FASTA should not be rescanned"),
            ):
                reclassified_statistics = get_scheme_statistics(
                    scheme_directory_path,
                    frequent_peak_relative_threshold_percentage=75,
                )[0]

            cached_data = json.loads(
                legacy_cache_file_path.read_text(encoding="utf-8")
            )

        self.assertEqual(
            reclassified_statistics.length_distribution_category,
            "Single-mode",
        )
        self.assertEqual(cached_data["cache_version"], SCHEME_STATISTICS_CACHE_VERSION)
        self.assertEqual(
            cached_data["locus_statistics"][0]["allele_length_frequencies"],
            [[100, 4], [110, 2]],
        )


if __name__ == "__main__":
    unittest.main()
