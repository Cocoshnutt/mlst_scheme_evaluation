"""Tests that analysis remains usable when cache persistence is unavailable."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scheme import alignment, allele_info, bsr, info, similarity, stats
from scheme.bsr import SchemeRepresentativeBsrAnalysis
from scheme.similarity import SchemeSimilarityAnalysis


class StaticAnalysisCacheFailureTest(unittest.TestCase):
    """Protect non-fatal cache-write behavior for FASTA-only analyses."""

    def test_scheme_information_survives_cache_write_failure(self) -> None:
        """Return calculated aggregate information from a read-only input."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "locus.fasta").write_text(
                ">allele_1\nATGTAA\n",
                encoding="utf-8",
            )
            with patch.object(
                info,
                "_write_cached_scheme_analysis",
                side_effect=PermissionError("read-only scheme"),
            ):
                scheme_information = info.get_scheme_information(
                    scheme_directory_path
                )

        self.assertEqual(scheme_information.locus_count, 1)
        self.assertEqual(scheme_information.allele_count, 1)

    def test_locus_statistics_survive_cache_write_failure(self) -> None:
        """Return calculated length statistics from a read-only input."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "locus.fasta").write_text(
                ">allele_1\nATGTAA\n",
                encoding="utf-8",
            )
            with patch.object(
                stats,
                "_write_cached_scheme_statistics",
                side_effect=PermissionError("read-only scheme"),
            ):
                locus_statistics = stats.get_scheme_statistics(
                    scheme_directory_path
                )

        self.assertEqual(len(locus_statistics), 1)
        self.assertEqual(locus_statistics[0].allele_count, 1)

    def test_allele_information_survives_cache_write_failure(self) -> None:
        """Return calculated per-locus information from a read-only input."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "locus.fasta").write_text(
                ">allele_1\nATGTAA\n",
                encoding="utf-8",
            )
            with patch.object(
                allele_info,
                "_write_cached_analysis",
                side_effect=PermissionError("read-only scheme"),
            ):
                locus_information = (
                    allele_info.get_scheme_allele_information(
                        scheme_directory_path
                    )
                )

        self.assertEqual(len(locus_information), 1)
        self.assertEqual(locus_information[0].valid_allele_count, 1)


class ExpensiveAnalysisCacheFailureTest(unittest.TestCase):
    """Protect non-fatal cache checkpoint and final-write behavior."""

    def test_alignment_checkpoint_write_failure_is_ignored(self) -> None:
        """Do not fail a completed alignment because persistence is unavailable."""
        with tempfile.TemporaryDirectory() as temporary_directory, patch.object(
            alignment,
            "_write_cached_alignment_cache",
            side_effect=PermissionError("read-only scheme"),
        ):
            alignment._try_write_cached_alignment_cache(
                Path(temporary_directory),
                {},
            )

    def test_similarity_cache_write_failure_is_ignored(self) -> None:
        """Do not fail completed BLAST measurements when JSON cannot be written."""
        analysis = SchemeSimilarityAnalysis(
            locus_count=0,
            nucleotide_selected_locus_count=0,
            amino_acid_selected_locus_count=0,
            nucleotide_possible_directed_pair_count=0,
            amino_acid_possible_directed_pair_count=0,
            matches=(),
            blastn_version="blastn: test",
            blastp_version="blastp: test",
            method_identifier="test",
            blast_evalue=0.001,
            reused_from_cache=False,
        )
        with tempfile.TemporaryDirectory() as temporary_directory, patch.object(
            similarity.tempfile,
            "NamedTemporaryFile",
            side_effect=PermissionError("read-only scheme"),
        ):
            similarity._write_cached_analysis(
                scheme_directory_path=Path(temporary_directory),
                source_signature=(),
                analysis=analysis,
            )

    def test_representative_bsr_cache_write_failure_is_ignored(self) -> None:
        """Do not fail completed BSR measurements when JSON cannot be written."""
        analysis = SchemeRepresentativeBsrAnalysis(
            locus_measurements=(),
            inter_locus_measurements=(),
            blastp_version="blastp: test",
            method_identifier="test",
            blast_evalue=0.001,
            reused_from_cache=False,
        )
        with tempfile.TemporaryDirectory() as temporary_directory, patch.object(
            bsr.tempfile,
            "NamedTemporaryFile",
            side_effect=PermissionError("read-only scheme"),
        ):
            bsr._write_cached_analysis(
                scheme_directory_path=Path(temporary_directory),
                source_signature=(),
                analysis=analysis,
            )


if __name__ == "__main__":
    unittest.main()
