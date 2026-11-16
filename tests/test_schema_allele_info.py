"""Tests for per-locus empty and duplicate nucleotide sequence metrics."""

import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from schema import allele_info
from schema.allele_info import (
    SCHEME_ALLELE_INFORMATION_CACHE_FILENAME,
    SCHEME_ALLELE_INFORMATION_CACHE_VERSION,
    get_scheme_allele_information,
    get_scheme_sequence_integrity_information,
)


class SchemeSequenceIntegrityTest(unittest.TestCase):
    """Exercise threshold-independent sequence-integrity analysis."""

    def test_counts_empty_and_within_locus_duplicate_sequences(self) -> None:
        """Normalize records and keep cross-locus matches independent."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            (schema_directory_path / "empty_locus.fasta").write_text(
                "",
                encoding="utf-8",
            )
            (schema_directory_path / "integrity.fasta").write_text(
                (
                    ">integrity_1\n"
                    "ATG\n"
                    "TAA\n"
                    ">integrity_2\n"
                    "atgtaa\n"
                    ">integrity_3\n"
                    "\n"
                    ">integrity_4\n"
                    ">integrity_5\n"
                    "GTGTAG\n"
                ),
                encoding="utf-8",
            )
            (schema_directory_path / "unique.fasta").write_text(
                ">unique_1\nATGTAA\n>unique_2\nTTGTGA\n",
                encoding="utf-8",
            )

            locus_information = {
                information.locus_name: information
                for information in get_scheme_allele_information(
                    schema_directory_path
                )
            }

            self.assertEqual(locus_information["empty_locus"].total_allele_count, 0)
            self.assertEqual(locus_information["empty_locus"].empty_sequence_count, 0)
            self.assertEqual(locus_information["integrity"].total_allele_count, 5)
            self.assertEqual(locus_information["integrity"].empty_sequence_count, 2)
            self.assertEqual(
                locus_information["integrity"].within_locus_duplicate_sequence_count,
                1,
            )
            self.assertEqual(
                locus_information["unique"].within_locus_duplicate_sequence_count,
                0,
            )

            scheme_information = get_scheme_sequence_integrity_information(
                schema_directory_path
            )
            self.assertEqual(scheme_information.locus_count, 3)
            self.assertEqual(scheme_information.empty_sequence_count, 2)
            self.assertEqual(scheme_information.loci_with_empty_sequences_count, 2)
            self.assertAlmostEqual(
                scheme_information.loci_with_empty_sequences_percentage,
                200 / 3,
            )
            self.assertEqual(
                scheme_information.loci_with_empty_sequences,
                ("empty_locus", "integrity"),
            )
            self.assertEqual(
                scheme_information.within_locus_duplicate_sequence_count,
                1,
            )
            self.assertEqual(
                scheme_information.loci_with_within_locus_duplicate_sequences_count,
                1,
            )
            self.assertAlmostEqual(
                scheme_information.loci_with_within_locus_duplicate_sequences_percentage,
                100 / 3,
            )
            self.assertEqual(
                scheme_information.loci_with_within_locus_duplicate_sequences,
                ("integrity",),
            )

    def test_three_identical_records_count_two_duplicates(self) -> None:
        """Count redundant records rather than sequence groups or pairs."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            (schema_directory_path / "triplicate.fasta").write_text(
                (
                    ">triplicate_1\nATGTAA\n"
                    ">triplicate_2\natgtaa\n"
                    ">triplicate_3\nATG\nTAA\n"
                ),
                encoding="utf-8",
            )

            locus_information = get_scheme_allele_information(
                schema_directory_path
            )[0]

            self.assertEqual(
                locus_information.within_locus_duplicate_sequence_count,
                2,
            )

    def test_minimum_length_does_not_change_integrity_metrics(self) -> None:
        """Keep the new metrics independent of the reactive length setting."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            (schema_directory_path / "locus.fasta").write_text(
                ">locus_1\nATGTAA\n>locus_2\natgtaa\n>locus_3\n",
                encoding="utf-8",
            )

            no_minimum = get_scheme_allele_information(
                schema_directory_path,
                minimum_sequence_length=0,
            )[0]
            long_minimum = get_scheme_allele_information(
                schema_directory_path,
                minimum_sequence_length=100,
            )[0]

            self.assertEqual(no_minimum.alleles_below_minimum_length_count, 0)
            self.assertEqual(long_minimum.alleles_below_minimum_length_count, 3)
            self.assertEqual(
                no_minimum.empty_sequence_count,
                long_minimum.empty_sequence_count,
            )
            self.assertEqual(
                no_minimum.within_locus_duplicate_sequence_count,
                long_minimum.within_locus_duplicate_sequence_count,
            )

    def test_reads_gzip_compressed_fasta(self) -> None:
        """Apply the same normalized duplicate rule to gzip input."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            compressed_fasta_file_path = schema_directory_path / "locus.fasta.gz"
            with gzip.open(
                compressed_fasta_file_path,
                mode="wt",
                encoding="utf-8",
            ) as compressed_fasta_file:
                compressed_fasta_file.write(
                    ">locus_1\nGTGTAG\n>locus_2\ngtgtag\n"
                )

            locus_information = get_scheme_allele_information(
                schema_directory_path
            )[0]

            self.assertEqual(
                locus_information.within_locus_duplicate_sequence_count,
                1,
            )

    def test_empty_scheme_has_zero_percentages(self) -> None:
        """Avoid division by zero when no locus FASTA files exist."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_information = get_scheme_sequence_integrity_information(
                temporary_directory
            )

            self.assertEqual(scheme_information.locus_count, 0)
            self.assertEqual(
                scheme_information.loci_with_empty_sequences_percentage,
                0.0,
            )
            self.assertEqual(
                scheme_information.loci_with_within_locus_duplicate_sequences_percentage,
                0.0,
            )

    def test_integrity_metrics_are_loaded_from_cache(self) -> None:
        """Persist the new static metrics and reuse them without rescanning."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            (schema_directory_path / "locus.fasta").write_text(
                ">locus_1\nATGTAA\n>locus_2\nATGTAA\n",
                encoding="utf-8",
            )
            get_scheme_allele_information(schema_directory_path)

            cache_file_path = (
                schema_directory_path
                / SCHEME_ALLELE_INFORMATION_CACHE_FILENAME
            )
            cached_data = json.loads(cache_file_path.read_text(encoding="utf-8"))
            self.assertEqual(
                cached_data["cache_version"],
                SCHEME_ALLELE_INFORMATION_CACHE_VERSION,
            )
            self.assertEqual(
                cached_data["locus_information"][0][
                    "within_locus_duplicate_sequence_count"
                ],
                1,
            )

            with patch.object(
                allele_info,
                "_analyze_locus",
                side_effect=AssertionError("FASTA was rescanned"),
            ):
                cached_information = get_scheme_allele_information(
                    schema_directory_path
                )[0]

            self.assertEqual(
                cached_information.within_locus_duplicate_sequence_count,
                1,
            )


if __name__ == "__main__":
    unittest.main()
