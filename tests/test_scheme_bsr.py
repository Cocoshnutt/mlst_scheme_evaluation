"""Tests for representative-versus-allele protein BSR analysis."""

import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scheme import bsr
from scheme.bsr import (
    SCHEME_REPRESENTATIVE_BSR_CACHE_FILENAME,
    SCHEME_REPRESENTATIVE_BSR_CACHE_FILENAMES,
    AlleleBlastScoreRatio,
    InterLocusBlastScoreRatio,
    LocusRepresentativeBsrMeasurements,
    SchemeRepresentativeBsrAnalysis,
    _collect_protein_alleles,
    _CollectedLocus,
    _create_bsr_measurements,
    _create_high_bsr_groups,
    _ProteinAllele,
    classify_representative_bsr_analysis,
    get_scheme_representative_bsr_analysis,
)
from scheme.similarity import _find_all_fasta_file_paths


class ProteinAlleleCollectionTest(unittest.TestCase):
    """Test representative selection and the protein evaluation denominator."""

    def test_collects_only_valid_translated_alleles(self) -> None:
        """Exclude invalid records and select the modal-length valid CDS."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "gene_a.fasta").write_text(
                ">outlier\nATGAAAAAATAA\n"
                ">modal_first\nATGCCCTAA\n"
                ">modal_second\nATGGGGTAA\n"
                ">invalid\nATGNNNTAA\n",
                encoding="utf-8",
            )

            collected_loci, protein_alleles = _collect_protein_alleles(
                scheme_directory_path,
                _find_all_fasta_file_paths(scheme_directory_path),
            )

        self.assertEqual(len(collected_loci), 1)
        self.assertEqual(collected_loci[0].total_allele_count, 4)
        self.assertEqual(len(protein_alleles), 3)
        self.assertEqual(
            collected_loci[0].representative_allele_identifier,
            "modal_first",
        )
        self.assertEqual(
            [protein_allele.protein_sequence for protein_allele in protein_alleles],
            ["MKK", "MP", "MG"],
        )


class BsrMeasurementTest(unittest.TestCase):
    """Test directional normalization, no-hit handling, and pair reduction."""

    def setUp(self) -> None:
        """Create two small loci with stable BLAST-safe identifiers."""
        self.gene_a_first = _ProteinAllele("a1", "gene_a", "gene_a_1", "MP")
        self.gene_a_second = _ProteinAllele("a2", "gene_a", "gene_a_2", "MA")
        self.gene_b_first = _ProteinAllele("b1", "gene_b", "gene_b_1", "MP")
        self.collected_loci = (
            _CollectedLocus(
                locus_name="gene_a",
                representative_internal_identifier="a1",
                representative_allele_identifier="gene_a_1",
                total_allele_count=3,
                protein_alleles=(self.gene_a_first, self.gene_a_second),
            ),
            _CollectedLocus(
                locus_name="gene_b",
                representative_internal_identifier="b1",
                representative_allele_identifier="gene_b_1",
                total_allele_count=1,
                protein_alleles=(self.gene_b_first,),
            ),
        )
        self.protein_alleles = (
            self.gene_a_first,
            self.gene_a_second,
            self.gene_b_first,
        )

    def test_uses_representative_self_score_and_retains_missing_hits(self) -> None:
        """Normalize by each query self score and keep absent hits as None."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            blast_output_file_path = Path(temporary_directory) / "hits.tsv"
            blast_output_file_path.write_text(
                "a1\ta1\t100\n"
                "a1\ta2\t60\n"
                "a1\tb1\t70\n"
                "b1\tb1\t80\n"
                "b1\ta1\t40\n"
                "b1\ta2\t56\n",
                encoding="utf-8",
            )

            locus_measurements, inter_locus_measurements = (
                _create_bsr_measurements(
                    collected_loci=self.collected_loci,
                    all_protein_alleles=self.protein_alleles,
                    blast_output_file_path=blast_output_file_path,
                )
            )

        gene_a_measurements = locus_measurements[0]
        self.assertEqual(gene_a_measurements.representative_self_score, 100.0)
        self.assertEqual(
            [
                measurement.blast_score_ratio
                for measurement in gene_a_measurements.allele_measurements
            ],
            [1.0, 0.6],
        )
        self.assertEqual(gene_a_measurements.unevaluable_allele_count, 1)
        self.assertEqual(len(inter_locus_measurements), 2)
        gene_b_query = next(
            measurement
            for measurement in inter_locus_measurements
            if measurement.query_locus_name == "gene_b"
        )
        self.assertEqual(gene_b_query.subject_allele_identifier, "gene_a_2")
        self.assertEqual(gene_b_query.blast_score_ratio, 0.7)

    def test_no_reported_hit_is_not_converted_to_zero(self) -> None:
        """Preserve unavailable local similarity separately from low BSR."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            blast_output_file_path = Path(temporary_directory) / "hits.tsv"
            blast_output_file_path.write_text(
                "a1\ta1\t100\n"
                "b1\tb1\t80\n",
                encoding="utf-8",
            )

            locus_measurements, _ = _create_bsr_measurements(
                collected_loci=self.collected_loci,
                all_protein_alleles=self.protein_alleles,
                blast_output_file_path=blast_output_file_path,
            )

        self.assertIsNone(
            locus_measurements[0].allele_measurements[1].blast_score_ratio
        )


class BsrClassificationTest(unittest.TestCase):
    """Test inclusive thresholding, pair collapse, and connected groups."""

    @staticmethod
    def _analysis() -> SchemeRepresentativeBsrAnalysis:
        """Return raw measurements with one unavailable locus and three edges."""
        return SchemeRepresentativeBsrAnalysis(
            locus_measurements=(
                LocusRepresentativeBsrMeasurements(
                    locus_name="gene_a",
                    representative_allele_identifier="gene_a_1",
                    representative_self_score=100.0,
                    total_allele_count=3,
                    valid_translated_allele_count=2,
                    allele_measurements=(
                        AlleleBlastScoreRatio("gene_a_1", 1.0),
                        AlleleBlastScoreRatio("gene_a_2", 0.6),
                    ),
                ),
                LocusRepresentativeBsrMeasurements(
                    locus_name="gene_b",
                    representative_allele_identifier="gene_b_1",
                    representative_self_score=80.0,
                    total_allele_count=2,
                    valid_translated_allele_count=2,
                    allele_measurements=(
                        AlleleBlastScoreRatio("gene_b_1", 1.0),
                        AlleleBlastScoreRatio("gene_b_2", None),
                    ),
                ),
                LocusRepresentativeBsrMeasurements(
                    locus_name="gene_d",
                    representative_allele_identifier="gene_d_1",
                    representative_self_score=None,
                    total_allele_count=1,
                    valid_translated_allele_count=1,
                    allele_measurements=(
                        AlleleBlastScoreRatio("gene_d_1", None),
                    ),
                ),
            ),
            inter_locus_measurements=(
                InterLocusBlastScoreRatio(
                    "gene_a", "gene_b", "gene_a_1", "gene_b_2", 70.0, 0.7
                ),
                InterLocusBlastScoreRatio(
                    "gene_b", "gene_a", "gene_b_1", "gene_a_1", 56.0, 0.7
                ),
                InterLocusBlastScoreRatio(
                    "gene_b", "gene_c", "gene_b_1", "gene_c_1", 64.0, 0.8
                ),
                InterLocusBlastScoreRatio(
                    "gene_d", "gene_e", "gene_d_1", "gene_e_1", 50.0, 0.5
                ),
            ),
            blastp_version="blastp: test",
            method_identifier="test",
            blast_evalue=0.001,
            reused_from_cache=False,
        )

    def test_threshold_is_inclusive_and_invalid_alleles_are_not_denominator(self) -> None:
        """Count BSR exactly 0.6 as covered and omit invalid records."""
        classification = classify_representative_bsr_analysis(
            self._analysis(),
            blast_score_ratio_threshold=0.6,
        )

        gene_a_coverage = classification.locus_coverages[0]
        self.assertEqual(gene_a_coverage.covered_allele_count, 2)
        self.assertEqual(gene_a_coverage.coverage_percentage, 100.0)
        self.assertEqual(gene_a_coverage.unevaluable_allele_count, 1)
        self.assertEqual(classification.measurable_allele_count, 4)
        self.assertEqual(classification.covered_allele_count, 3)
        self.assertEqual(classification.unavailable_valid_allele_count, 1)

    def test_reciprocal_pairs_collapse_and_transitive_edges_form_group(self) -> None:
        """Create one A-B-C component from two qualifying unordered pairs."""
        classification = classify_representative_bsr_analysis(
            self._analysis(),
            blast_score_ratio_threshold=0.6,
        )

        self.assertEqual(len(classification.high_bsr_pairs), 2)
        gene_a_gene_b_pair = classification.high_bsr_pairs[0]
        self.assertEqual(
            (gene_a_gene_b_pair.first_locus_name, gene_a_gene_b_pair.second_locus_name),
            ("gene_a", "gene_b"),
        )
        self.assertEqual(
            gene_a_gene_b_pair.evidence.representative_allele_identifier,
            "gene_a_1",
        )
        self.assertEqual(len(classification.high_bsr_groups), 1)
        self.assertEqual(
            classification.high_bsr_groups[0].locus_names,
            ("gene_a", "gene_b", "gene_c"),
        )
        self.assertEqual(classification.high_bsr_groups[0].qualifying_pair_count, 2)
        self.assertEqual(classification.high_bsr_groups[0].maximum_blast_score_ratio, 0.8)

    def test_threshold_change_reclassifies_same_raw_analysis(self) -> None:
        """Change classification without changing the raw measurement object."""
        analysis = self._analysis()
        default_result = classify_representative_bsr_analysis(
            analysis,
            blast_score_ratio_threshold=0.6,
        )
        strict_result = classify_representative_bsr_analysis(
            analysis,
            blast_score_ratio_threshold=0.75,
        )

        self.assertEqual(len(default_result.high_bsr_pairs), 2)
        self.assertEqual(len(strict_result.high_bsr_pairs), 1)
        self.assertIs(default_result.locus_coverages[0].__class__, strict_result.locus_coverages[0].__class__)

    def test_empty_pairs_create_no_groups(self) -> None:
        """Return an empty group tuple without a maximum-on-empty error."""
        self.assertEqual(_create_high_bsr_groups(()), ())


class RepresentativeBsrCacheTest(unittest.TestCase):
    """Test cache reuse and FASTA-signature invalidation."""

    @staticmethod
    def _empty_analysis_result(*args, **kwargs):
        """Return a deterministic no-hit analysis for a collected scheme."""
        collected_loci = kwargs["collected_loci"]
        return (
            tuple(
                LocusRepresentativeBsrMeasurements(
                    locus_name=locus.locus_name,
                    representative_allele_identifier=(
                        locus.representative_allele_identifier
                    ),
                    representative_self_score=10.0,
                    total_allele_count=locus.total_allele_count,
                    valid_translated_allele_count=len(locus.protein_alleles),
                    allele_measurements=tuple(
                        AlleleBlastScoreRatio(
                            protein_allele.original_identifier,
                            1.0,
                        )
                        for protein_allele in locus.protein_alleles
                    ),
                )
                for locus in collected_loci
            ),
            (),
        )

    def test_reuses_cache_and_invalidates_when_source_changes(self) -> None:
        """Avoid rerunning BLAST until the FASTA source signature changes."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            locus_file_path = scheme_directory_path / "gene_a.fasta"
            locus_file_path.write_text(
                ">gene_a_1\nATGCCCTAA\n",
                encoding="utf-8",
            )
            with (
                patch.object(
                    bsr,
                    "_get_blast_version",
                    return_value="blastp: test",
                ),
                patch.object(
                    bsr,
                    "_run_blastp_analysis",
                    side_effect=self._empty_analysis_result,
                ) as run_blastp_analysis,
            ):
                first_analysis = get_scheme_representative_bsr_analysis(
                    scheme_directory_path
                )
                legacy_cache_file_path = (
                    scheme_directory_path
                    / SCHEME_REPRESENTATIVE_BSR_CACHE_FILENAMES[1]
                )
                (
                    scheme_directory_path
                    / SCHEME_REPRESENTATIVE_BSR_CACHE_FILENAME
                ).replace(legacy_cache_file_path)
                cached_analysis = get_scheme_representative_bsr_analysis(
                    scheme_directory_path
                )
                locus_file_path.write_text(
                    ">gene_a_1\nATGCCCTAA\n"
                    ">gene_a_2\nATGGGGTAA\n",
                    encoding="utf-8",
                )
                recalculated_analysis = get_scheme_representative_bsr_analysis(
                    scheme_directory_path
                )

        self.assertFalse(first_analysis.reused_from_cache)
        self.assertTrue(cached_analysis.reused_from_cache)
        self.assertFalse(recalculated_analysis.reused_from_cache)
        self.assertEqual(run_blastp_analysis.call_count, 2)


@unittest.skipUnless(
    shutil.which("blastp") and shutil.which("makeblastdb"),
    "BLAST+ executables are required for the integration test",
)
class RepresentativeBsrIntegrationTest(unittest.TestCase):
    """Exercise the public API with a real batched protein BLAST search."""

    def test_tiny_scheme_analysis_and_cache_reuse(self) -> None:
        """Measure identical proteins across loci and then reuse the cache."""
        coding_codons = (
            "GCT",
            "TGT",
            "GAT",
            "GAA",
            "TTT",
            "GGT",
            "CAT",
            "ATT",
            "AAA",
            "CTT",
            "AAT",
            "CCT",
            "CAA",
            "CGT",
            "TCT",
            "ACT",
            "GTT",
            "TGG",
            "TAT",
        ) * 4
        coding_sequence = "ATG" + "".join(coding_codons) + "TAA"
        synonymous_sequence = coding_sequence.replace("GCT", "GCC", 1)
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "gene_a.fasta").write_text(
                f">gene_a_1\n{coding_sequence}\n"
                f">gene_a_2\n{synonymous_sequence}\n",
                encoding="utf-8",
            )
            (scheme_directory_path / "gene_b.fasta").write_text(
                f">gene_b_1\n{coding_sequence}\n",
                encoding="utf-8",
            )

            analysis = get_scheme_representative_bsr_analysis(
                scheme_directory_path,
                maximum_threads=1,
            )
            cached_analysis = get_scheme_representative_bsr_analysis(
                scheme_directory_path,
                maximum_threads=1,
            )
            classification = classify_representative_bsr_analysis(analysis)

        self.assertFalse(analysis.reused_from_cache)
        self.assertTrue(cached_analysis.reused_from_cache)
        self.assertEqual(classification.covered_allele_count, 3)
        self.assertEqual(classification.measurable_allele_count, 3)
        self.assertEqual(len(classification.high_bsr_pairs), 1)
        self.assertEqual(classification.high_bsr_groups[0].locus_names, ("gene_a", "gene_b"))


if __name__ == "__main__":
    unittest.main()
