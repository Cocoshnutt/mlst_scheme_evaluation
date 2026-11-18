"""Tests for the retained per-locus allele-information behavior."""

import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scheme import allele_info
from scheme.allele_info import (
    SCHEME_ALLELE_INFORMATION_CACHE_FILENAME,
    SCHEME_ALLELE_INFORMATION_CACHE_VERSION,
    get_scheme_allele_information,
)


class SchemeAlleleInformationTest(unittest.TestCase):
    """Exercise cache reuse and supported allele-information fields."""

    def test_counts_nucleotide_alleles_and_distinct_valid_translations(
        self,
    ) -> None:
        """Collapse synonymous valid alleles only at the protein level."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "locus.fasta").write_text(
                (
                    ">locus_1\nATGGCTTAA\n"
                    ">locus_2\nATGGCCTAA\n"
                    ">locus_3\nATGAAATAA\n"
                    ">locus_4\nATGNNNTAA\n"
                ),
                encoding="utf-8",
            )

            locus_information = get_scheme_allele_information(
                scheme_directory_path
            )[0]

        self.assertEqual(locus_information.total_allele_count, 4)
        self.assertEqual(locus_information.valid_allele_count, 3)
        self.assertEqual(locus_information.distinct_protein_allele_count, 2)

    def test_minimum_length_reuses_cached_static_analysis(self) -> None:
        """Recalculate only the minimum-length-dependent count."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "locus.fasta").write_text(
                ">locus_1\nATGTAA\n>locus_2\nGTGTAG\n",
                encoding="utf-8",
            )
            no_minimum = get_scheme_allele_information(
                scheme_directory_path,
                minimum_sequence_length=0,
            )[0]
            with patch.object(
                allele_info,
                "_analyze_locus",
                side_effect=AssertionError("FASTA was rescanned"),
            ):
                longer_minimum = get_scheme_allele_information(
                    scheme_directory_path,
                    minimum_sequence_length=7,
                )[0]

        self.assertEqual(no_minimum.alleles_below_minimum_length_count, 0)
        self.assertEqual(longer_minimum.alleles_below_minimum_length_count, 2)
        self.assertEqual(no_minimum.total_allele_count, 2)
        self.assertEqual(longer_minimum.total_allele_count, 2)

    def test_cache_and_public_result_omit_removed_sequence_criteria(self) -> None:
        """Do not expose empty, duplicate, or shared-sequence measurements."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "locus.fasta").write_text(
                ">locus_1\nATGTAA\n",
                encoding="utf-8",
            )
            locus_information = get_scheme_allele_information(
                scheme_directory_path
            )[0]
            cached_data = json.loads(
                (
                    scheme_directory_path
                    / SCHEME_ALLELE_INFORMATION_CACHE_FILENAME
                ).read_text(encoding="utf-8")
            )

        self.assertEqual(
            cached_data["cache_version"],
            SCHEME_ALLELE_INFORMATION_CACHE_VERSION,
        )
        self.assertNotIn(
            "empty_sequence_count",
            cached_data["locus_information"][0],
        )
        self.assertNotIn(
            "within_locus_duplicate_sequence_count",
            cached_data["locus_information"][0],
        )
        self.assertNotIn(
            "missing_allele_identifier_count",
            cached_data["locus_information"][0],
        )
        self.assertNotIn("shared_sequence_groups", cached_data)
        self.assertFalse(hasattr(locus_information, "empty_sequence_count"))

    def test_reads_gzip_compressed_fasta(self) -> None:
        """Retain support for gzip-compressed locus files."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            with gzip.open(
                scheme_directory_path / "locus.fasta.gz",
                mode="wt",
                encoding="utf-8",
            ) as compressed_fasta_file:
                compressed_fasta_file.write(">locus_1\nGTGTAG\n")

            locus_information = get_scheme_allele_information(
                scheme_directory_path
            )[0]

        self.assertEqual(locus_information.total_allele_count, 1)
        self.assertEqual(locus_information.valid_allele_count, 1)


if __name__ == "__main__":
    unittest.main()
