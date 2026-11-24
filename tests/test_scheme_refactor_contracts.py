"""Characterization contracts for the shared-primitives refactor.

These tests deliberately describe behavior that is currently spread across
multiple intrinsic-analysis modules.  They provide a stable comparison point
while FASTA parsing, sequence validation, translation, and representative
selection move to shared modules.
"""

import gzip
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path

from scheme.allele_info import get_scheme_allele_information
from scheme.fasta import iterate_fasta_records, read_fasta_records
from scheme.info import get_scheme_information
from scheme.selection import select_allele_closest_to_modal_length
from scheme.sequence import evaluate_allele_sequence, translate_valid_allele
from scheme.stats import get_scheme_statistics


class FastaParsingContractTest(unittest.TestCase):
    """Protect record boundaries, identifiers, and empty-record handling."""

    def test_header_defined_records_have_consistent_boundaries(self) -> None:
        """Ignore preamble text and retain empty records between headers."""
        fasta_lines = (
            "text before the first header is ignored\n",
            ">first description with spaces\n",
            "atg\n",
            "aaa\n",
            ">empty description\n",
            "\n",
            ">third\n",
            "ttgcccTGA\n",
        )

        fasta_records = tuple(iterate_fasta_records(fasta_lines))

        self.assertEqual(
            tuple((record.identifier, record.sequence) for record in fasta_records),
            (
                ("first", "atgaaa"),
                ("empty", ""),
                ("third", "ttgcccTGA"),
            ),
        )
    def test_plain_and_gzip_readers_return_the_same_records(self) -> None:
        """Treat compressed and uncompressed FASTA content identically."""
        fasta_text = ">one label\nATGAAA\nTAA\n>empty\n"
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory_path = Path(temporary_directory)
            plain_path = directory_path / "locus.fasta"
            gzip_path = directory_path / "locus.fasta.gz"
            plain_path.write_text(fasta_text, encoding="utf-8")
            with gzip.open(gzip_path, mode="wt", encoding="utf-8") as output_file:
                output_file.write(fasta_text)

            plain_records = read_fasta_records(plain_path)
            gzip_records = read_fasta_records(gzip_path)

        self.assertEqual(
            tuple((record.identifier, record.sequence) for record in plain_records),
            (("one", "ATGAAATAA"), ("empty", "")),
        )
        self.assertEqual(gzip_records, plain_records)


class SequenceContractTest(unittest.TestCase):
    """Protect exact validation flags and bacterial-CDS translation."""

    def test_invalidity_categories_can_overlap(self) -> None:
        """Retain each category independently rather than deduplicating records."""
        validation_cases = {
            "": (False, False, True, False),
            "ATGNNNTAA": (False, True, False, False),
            "ATGTAACCCTAA": (False, False, False, True),
            "CTGNNNTA": (True, True, True, False),
        }

        for sequence, expected_flags in validation_cases.items():
            with self.subTest(sequence=sequence):
                validation = evaluate_allele_sequence(sequence)
                self.assertEqual(
                    (
                        validation.has_incomplete_orf,
                        validation.has_ambiguous_bases,
                        validation.is_missing_start_or_stop_codon,
                        validation.has_in_frame_stop_codon,
                    ),
                    expected_flags,
                )

    def test_translation_normalizes_bacterial_starts_and_removes_stop(self) -> None:
        """Translate all accepted starts as methionine and omit terminal stops."""
        self.assertEqual(translate_valid_allele("ATGGCCTAA"), "MA")
        self.assertEqual(translate_valid_allele("GTGGCCTAG"), "MA")
        self.assertEqual(translate_valid_allele("TTGGGGTGA"), "MG")
        self.assertEqual(translate_valid_allele("ATGTAA"), "M")

    def test_shared_selection_uses_shortest_mode_and_fasta_order(self) -> None:
        """Keep the representative rule independent of BLAST modules."""
        records = tuple(
            iterate_fasta_records(
                (
                    ">long\n",
                    "ATGAAAAAATAA\n",
                    ">first_short\n",
                    "ATGCCCTAA\n",
                    ">second_short\n",
                    "ATGGGGTAA\n",
                )
            )
        )

        selected_record = select_allele_closest_to_modal_length(
            records,
            require_valid_coding_sequence=True,
        )

        self.assertIsNotNone(selected_record)
        self.assertEqual(selected_record.identifier, "first_short")


class PublicResultGoldenContractTest(unittest.TestCase):
    """Lock representative static public results before moving primitives."""

    def test_static_public_results_match_golden_fixture(self) -> None:
        """Preserve counts, denominators, names, and threshold boundaries."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "alpha.fasta").write_text(
                (
                    ">valid_atg description\nATGGCTTAA\n"
                    ">valid_gtg\nGTGGCCTAG\n"
                    ">ambiguous\nATGNNNTAA\n"
                    ">internal\nATGTAACCCTAA\n"
                    ">short\nATGTAA\n"
                ),
                encoding="utf-8",
            )
            with gzip.open(
                scheme_directory_path / "nested.fna.gz",
                mode="wt",
                encoding="utf-8",
            ) as compressed_fasta_file:
                compressed_fasta_file.write(
                    ">valid_ttg\nTTGGGGTGA\n"
                    ">missing\nCTGAAATAA\n"
                    ">incomplete\nATGAAATA\n"
                    ">empty\n"
                )

            scheme_information = get_scheme_information(
                scheme_directory_path,
                minimum_sequence_length=9,
            )
            locus_statistics = get_scheme_statistics(scheme_directory_path)
            locus_information = get_scheme_allele_information(
                scheme_directory_path,
                minimum_sequence_length=9,
            )

        self.assertEqual(
            asdict(scheme_information),
            {
                "locus_count": 2,
                "allele_count": 9,
                "valid_allele_count": 4,
                "invalid_allele_count": 6,
                "incomplete_orf_allele_count": 1,
                "ambiguous_bases_allele_count": 1,
                "missing_start_or_stop_codon_allele_count": 3,
                "in_frame_stop_codon_allele_count": 1,
                "alleles_below_minimum_length_count": 3,
                "minimum_sequence_length": 9,
                "alleles_below_locus_mode_threshold_count": 2,
                "alleles_above_locus_mode_threshold_count": 1,
                "alleles_below_locus_mode_threshold_percentage": 200 / 9,
                "alleles_above_locus_mode_threshold_percentage": 100 / 9,
                "guanine_cytosine_content_percentage": 33.82352941176471,
                "canonical_nucleotide_base_count": 68,
                "excluded_nucleotide_base_count": 3,
            },
        )
        self.assertEqual(
            [
                (
                    statistics.locus_name,
                    statistics.allele_count,
                    statistics.modal_allele_lengths,
                    statistics.mean_allele_length,
                    statistics.alleles_below_mode_threshold_count,
                    statistics.alleles_above_mode_threshold_count,
                    statistics.length_distribution_category,
                    statistics.canonical_nucleotide_base_count,
                    statistics.excluded_nucleotide_base_count,
                )
                for statistics in locus_statistics
            ],
            [
                (
                    "alpha.fasta",
                    5,
                    (9,),
                    9.0,
                    1,
                    1,
                    "Single-mode with outliers",
                    42,
                    3,
                ),
                (
                    "nested.fna.gz",
                    4,
                    (9,),
                    6.5,
                    1,
                    0,
                    "Multimodal with outliers",
                    26,
                    0,
                ),
            ],
        )
        self.assertEqual(
            [asdict(information) for information in locus_information],
            [
                {
                    "locus_name": "alpha",
                    "total_allele_count": 5,
                    "valid_allele_count": 3,
                    "invalid_allele_count": 2,
                    "validated_allele_proportion": 0.6,
                    "distinct_protein_allele_count": 2,
                    "incomplete_orf_allele_count": 0,
                    "ambiguous_bases_allele_count": 1,
                    "missing_start_or_stop_codon_allele_count": 0,
                    "in_frame_stop_codon_allele_count": 1,
                    "alleles_below_minimum_length_count": 1,
                    "minimum_sequence_length": 9,
                    "alleles_below_mode_threshold_count": 1,
                    "alleles_above_mode_threshold_count": 1,
                },
                {
                    "locus_name": "nested",
                    "total_allele_count": 4,
                    "valid_allele_count": 1,
                    "invalid_allele_count": 4,
                    "validated_allele_proportion": 0.25,
                    "distinct_protein_allele_count": 1,
                    "incomplete_orf_allele_count": 1,
                    "ambiguous_bases_allele_count": 0,
                    "missing_start_or_stop_codon_allele_count": 3,
                    "in_frame_stop_codon_allele_count": 0,
                    "alleles_below_minimum_length_count": 2,
                    "minimum_sequence_length": 9,
                    "alleles_below_mode_threshold_count": 1,
                    "alleles_above_mode_threshold_count": 0,
                },
            ],
        )


if __name__ == "__main__":
    unittest.main()
