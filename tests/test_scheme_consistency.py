"""Cross-module tests for scientific rules shared by intrinsic analyses."""

import gzip
import tempfile
import unittest
from pathlib import Path

from scheme import alignment, allele_info, info, similarity, stats
from scheme.bsr import _collect_protein_alleles


class FastaDiscoveryConsistencyTest(unittest.TestCase):
    """Ensure every analysis discovers the same recognized locus files."""

    def test_recognized_fasta_paths_are_identical_across_analyses(self) -> None:
        """Protect recursive, compressed, case-insensitive FASTA discovery."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            nested_directory_path = scheme_directory_path / "nested"
            nested_directory_path.mkdir()
            (scheme_directory_path / "first.fasta").write_text(
                ">first_1\nATGTAA\n",
                encoding="utf-8",
            )
            with gzip.open(
                nested_directory_path / "second.FNA.GZ",
                mode="wt",
                encoding="utf-8",
            ) as compressed_fasta_file:
                compressed_fasta_file.write(">second_1\nGTGTAG\n")
            (nested_directory_path / "ignored.txt").write_text(
                ">not_a_locus\nATGTAA\n",
                encoding="utf-8",
            )

            discovered_path_sets = (
                info._find_allele_file_paths(scheme_directory_path),
                stats._find_allele_file_paths(scheme_directory_path),
                allele_info._find_allele_file_paths(scheme_directory_path),
                alignment._find_allele_file_paths(scheme_directory_path),
                similarity._find_all_fasta_file_paths(scheme_directory_path),
            )

        expected_relative_paths = ("first.fasta", "nested/second.FNA.GZ")
        for discovered_paths in discovered_path_sets:
            with self.subTest(discovery_function=discovered_paths):
                self.assertEqual(
                    tuple(
                        path.relative_to(scheme_directory_path).as_posix()
                        for path in discovered_paths
                    ),
                    expected_relative_paths,
                )


class CodingSequenceConsistencyTest(unittest.TestCase):
    """Ensure modules agree on the coding-sequence validity definition."""

    def test_validation_corpus_has_identical_outcomes(self) -> None:
        """Protect all starts, invalidity categories, and normalization."""
        validation_cases = {
            "atg_valid": ("ATGAAATAA", True),
            "gtg_valid": ("GTGCCCTAG", True),
            "ttg_valid": ("TTGGGGTGA", True),
            "lowercase_whitespace": (" atg aaa \n taa ", True),
            "incomplete_orf": ("ATGAAATA", False),
            "ambiguous_base": ("ATGNNNTAA", False),
            "missing_start": ("CTGAAATAA", False),
            "missing_stop": ("ATGAAACCC", False),
            "internal_stop": ("ATGTAACCCTAA", False),
            "empty": ("", False),
        }

        for case_name, (allele_sequence, expected_validity) in (
            validation_cases.items()
        ):
            normalized_sequence = "".join(allele_sequence.split()).upper()
            with self.subTest(case=case_name):
                self.assertEqual(
                    info.is_valid_allele_sequence(allele_sequence),
                    expected_validity,
                )
                self.assertEqual(
                    allele_info._evaluate_allele_sequence(
                        normalized_sequence
                    ).is_valid,
                    expected_validity,
                )
                self.assertEqual(
                    similarity._is_valid_coding_sequence(normalized_sequence),
                    expected_validity,
                )


class RepresentativeSelectionConsistencyTest(unittest.TestCase):
    """Ensure similarity and representative BSR select the same valid CDS."""

    def test_selected_protein_alleles_match_bsr_representatives(self) -> None:
        """Protect shared modal-length, validity, and FASTA-order tie rules."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "gene_a.fasta").write_text(
                (
                    ">invalid_first\nATGNNNTAA\n"
                    ">valid_first\nATGCCCTAA\n"
                    ">valid_second\nATGGGGTAA\n"
                ),
                encoding="utf-8",
            )
            (scheme_directory_path / "gene_b.fasta").write_text(
                ">nucleotide_only\nATGNNNTAA\n",
                encoding="utf-8",
            )
            fasta_file_paths = similarity._find_all_fasta_file_paths(
                scheme_directory_path
            )

            (
                selected_nucleotide_alleles,
                selected_amino_acid_alleles,
                locus_count,
            ) = similarity._collect_selected_allele_sequences(
                scheme_directory_path,
                fasta_file_paths,
            )
            collected_loci, _ = _collect_protein_alleles(
                scheme_directory_path,
                fasta_file_paths,
            )

        similarity_protein_selection = {
            selected_allele.locus_name: selected_allele.original_identifier
            for selected_allele in selected_amino_acid_alleles
        }
        bsr_representative_selection = {
            collected_locus.locus_name: (
                collected_locus.representative_allele_identifier
            )
            for collected_locus in collected_loci
        }
        self.assertEqual(locus_count, 2)
        self.assertEqual(
            similarity_protein_selection,
            {"gene_a": "valid_first"},
        )
        self.assertEqual(
            bsr_representative_selection,
            {"gene_a": "valid_first", "gene_b": None},
        )
        self.assertEqual(
            {
                selected_allele.locus_name: selected_allele.original_identifier
                for selected_allele in selected_nucleotide_alleles
            },
            {
                "gene_a": "valid_first",
                "gene_b": "nucleotide_only",
            },
        )


if __name__ == "__main__":
    unittest.main()
