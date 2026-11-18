"""Tests for the one-pass shared static FASTA analysis."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scheme import allele_info, info, static_analysis, stats


class SharedStaticAnalysisTest(unittest.TestCase):
    """Protect single traversal, projections, and cache resilience."""

    def test_three_public_analyses_share_one_fasta_traversal(self) -> None:
        """Populate all established projections from one shared scan."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "alpha.fasta").write_text(
                ">alpha_1\nATGAAATAA\n>alpha_2\nATGNNNTAA\n",
                encoding="utf-8",
            )
            (scheme_directory_path / "beta.fna").write_text(
                ">beta_1\nGTGCCCTAG\n",
                encoding="utf-8",
            )
            original_reader = static_analysis._read_fasta_records
            with patch.object(
                static_analysis,
                "_read_fasta_records",
                wraps=original_reader,
            ) as read_fasta_records:
                scheme_information = info.get_scheme_information(
                    scheme_directory_path
                )
                locus_statistics = stats.get_scheme_statistics(
                    scheme_directory_path
                )
                locus_information = allele_info.get_scheme_allele_information(
                    scheme_directory_path
                )

            cache_filenames = (
                static_analysis.SCHEME_STATIC_ANALYSIS_CACHE_FILENAME,
                info.SCHEME_INFORMATION_CACHE_FILENAME,
                stats.SCHEME_STATISTICS_CACHE_FILENAME,
                allele_info.SCHEME_ALLELE_INFORMATION_CACHE_FILENAME,
            )
            self.assertTrue(
                all(
                    (scheme_directory_path / cache_filename).is_file()
                    for cache_filename in cache_filenames
                )
            )

        self.assertEqual(read_fasta_records.call_count, 2)
        self.assertEqual(scheme_information.locus_count, 2)
        self.assertEqual(scheme_information.allele_count, 3)
        self.assertEqual(
            scheme_information.allele_count,
            sum(item.allele_count for item in locus_statistics),
        )
        self.assertEqual(
            scheme_information.valid_allele_count,
            sum(item.valid_allele_count for item in locus_information),
        )

    def test_shared_cache_reuse_does_not_rescan_fasta(self) -> None:
        """Reuse compatible shared measurements directly."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "locus.fasta").write_text(
                ">allele_1\nATGTAA\n",
                encoding="utf-8",
            )
            first_analysis = static_analysis.get_shared_scheme_static_analysis(
                scheme_directory_path
            )
            with patch.object(
                static_analysis,
                "_analyze_locus",
                side_effect=AssertionError("FASTA was rescanned"),
            ):
                cached_analysis = (
                    static_analysis.get_shared_scheme_static_analysis(
                        scheme_directory_path
                    )
                )

        self.assertEqual(cached_analysis, first_analysis)

    def test_changed_source_signature_invalidates_shared_cache(self) -> None:
        """Recalculate after a locus file changes."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            locus_file_path = scheme_directory_path / "locus.fasta"
            locus_file_path.write_text(
                ">allele_1\nATGTAA\n",
                encoding="utf-8",
            )
            static_analysis.get_shared_scheme_static_analysis(
                scheme_directory_path
            )
            locus_file_path.write_text(
                ">allele_1\nATGTAA\n>allele_2\nGTGTAG\n",
                encoding="utf-8",
            )
            changed_analysis = static_analysis.get_shared_scheme_static_analysis(
                scheme_directory_path
            )

        self.assertEqual(changed_analysis.loci[0].allele_count, 2)

    def test_corrupt_cache_is_replaced_with_valid_measurements(self) -> None:
        """Treat malformed JSON as a cache miss rather than an analysis error."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "locus.fasta").write_text(
                ">allele_1\nATGTAA\n",
                encoding="utf-8",
            )
            cache_file_path = (
                scheme_directory_path
                / static_analysis.SCHEME_STATIC_ANALYSIS_CACHE_FILENAME
            )
            cache_file_path.write_text("{not-json", encoding="utf-8")
            analysis = static_analysis.get_shared_scheme_static_analysis(
                scheme_directory_path
            )
            cached_data = json.loads(cache_file_path.read_text(encoding="utf-8"))

        self.assertEqual(analysis.loci[0].valid_allele_count, 1)
        self.assertEqual(
            cached_data["cache_version"],
            static_analysis.SCHEME_STATIC_ANALYSIS_CACHE_VERSION,
        )

    def test_cache_write_failure_is_non_fatal(self) -> None:
        """Analyze readable FASTA files even when persistence is unavailable."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "locus.fasta").write_text(
                ">allele_1\nATGTAA\n",
                encoding="utf-8",
            )
            with patch.object(
                static_analysis,
                "_write_cached_analysis",
                side_effect=PermissionError("read-only scheme"),
            ):
                analysis = static_analysis.get_shared_scheme_static_analysis(
                    scheme_directory_path
                )

        self.assertEqual(analysis.loci[0].allele_count, 1)

    def test_public_entry_point_validates_the_scheme_path(self) -> None:
        """Report missing paths and files with the established error types."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_directory_path = Path(temporary_directory)
            missing_path = temporary_directory_path / "missing"
            file_path = temporary_directory_path / "not-a-directory"
            file_path.write_text("data", encoding="utf-8")

            with self.assertRaises(FileNotFoundError):
                static_analysis.get_shared_scheme_static_analysis(missing_path)
            with self.assertRaises(NotADirectoryError):
                static_analysis.get_shared_scheme_static_analysis(file_path)


if __name__ == "__main__":
    unittest.main()
