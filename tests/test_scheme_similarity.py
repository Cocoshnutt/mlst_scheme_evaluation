"""Tests for selected-allele cross-locus BLAST similarity analysis."""

import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scheme import similarity
from scheme.similarity import (
    SCHEME_SIMILARITY_CACHE_FILENAME,
    SCHEME_SIMILARITY_CACHE_FILENAMES,
    SEQUENCE_LEVEL_AMINO_ACID,
    SEQUENCE_LEVEL_NUCLEOTIDE,
    CandidateParalogThresholds,
    CrossLocusSimilarityMatch,
    SchemeSimilarityAnalysis,
    _BlastHsp,
    _collect_selected_allele_sequences,
    _create_similarity_matches,
    _find_all_fasta_file_paths,
    _SelectedAlleleSequence,
    _SimilarityLevelResult,
    classify_candidate_paralogous_locus_pairs,
    get_scheme_similarity_analysis,
)


class SelectedAlleleTest(unittest.TestCase):
    """Test modal-length selection and deterministic fallbacks."""

    def test_valid_allele_closest_to_modal_length_is_selected(self) -> None:
        """Prefer a valid modal-length allele over an earlier length outlier."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "gene_a.fasta").write_text(
                ">gene_a_1\nATGAAAAAATAA\n"
                ">gene_a_2\nATGCCCTAA\n"
                ">gene_a_3\nATGGGGTAA\n",
                encoding="utf-8",
            )

            (
                selected_nucleotide_alleles,
                selected_amino_acid_alleles,
                locus_count,
            ) = _collect_selected_allele_sequences(
                scheme_directory_path,
                _find_all_fasta_file_paths(scheme_directory_path),
            )

        self.assertEqual(locus_count, 1)
        self.assertEqual(selected_nucleotide_alleles[0].original_identifier, "gene_a_2")
        self.assertEqual(
            selected_amino_acid_alleles[0].original_identifier,
            "gene_a_2",
        )

    def test_shortest_tied_modal_length_and_fasta_order_break_ties(self) -> None:
        """Use the shortest tied mode, then retain FASTA order at equal distance."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "gene_b.fasta").write_text(
                ">gene_b_1\nATGAAAAAATAA\n"
                ">gene_b_2\nATGCCCTAA\n"
                ">gene_b_3\nATGGGGTAA\n"
                ">gene_b_4\nATGCCCC CCTAA\n".replace(" ", ""),
                encoding="utf-8",
            )

            (
                selected_nucleotide_alleles,
                selected_amino_acid_alleles,
                locus_count,
            ) = _collect_selected_allele_sequences(
                scheme_directory_path,
                _find_all_fasta_file_paths(scheme_directory_path),
            )

        self.assertEqual(locus_count, 1)
        self.assertEqual(
            selected_nucleotide_alleles[0].original_identifier,
            "gene_b_2",
        )
        self.assertEqual(
            selected_amino_acid_alleles[0].original_identifier,
            "gene_b_2",
        )
        self.assertEqual(selected_amino_acid_alleles[0].sequence, "MP")

    def test_invalid_cds_can_supply_only_a_nucleotide_fallback(self) -> None:
        """Keep IUPAC DNA measurable when no allele can be translated."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "fragment.fasta").write_text(
                ">fragment_1\nACGTR\n",
                encoding="utf-8",
            )

            (
                selected_nucleotide_alleles,
                selected_amino_acid_alleles,
                _,
            ) = _collect_selected_allele_sequences(
                scheme_directory_path,
                _find_all_fasta_file_paths(scheme_directory_path),
            )

        self.assertEqual(len(selected_nucleotide_alleles), 1)
        self.assertEqual(selected_amino_acid_alleles, ())


class SimilarityMetricTest(unittest.TestCase):
    """Test identity, asymmetric coverage, best-HSP, and directional BSR."""

    def setUp(self) -> None:
        """Create two selected nucleotide alleles with stable internal IDs."""
        self.selected_alleles = (
            _SelectedAlleleSequence(
                internal_identifier="sequence_000000001",
                locus_name="gene_a",
                original_identifier="gene_a_1",
                sequence="ACGTACGTAA",
            ),
            _SelectedAlleleSequence(
                internal_identifier="sequence_000000002",
                locus_name="gene_b",
                original_identifier="gene_b_1",
                sequence="TTACGTACGTAACCGGTTAA",
            ),
        )

    def test_metrics_are_directional_and_use_coordinate_coverage(self) -> None:
        """Normalize reciprocal raw scores to each query's own self score."""
        blast_hsps = (
            _BlastHsp(
                "sequence_000000001",
                "sequence_000000001",
                10,
                10,
                1,
                10,
                1,
                10,
                10,
                10,
                100.0,
                "plus",
            ),
            _BlastHsp(
                "sequence_000000002",
                "sequence_000000002",
                20,
                20,
                1,
                20,
                1,
                20,
                20,
                20,
                80.0,
                "plus",
            ),
            _BlastHsp(
                "sequence_000000001",
                "sequence_000000002",
                10,
                20,
                1,
                8,
                10,
                3,
                8,
                6,
                60.0,
                "minus",
            ),
            _BlastHsp(
                "sequence_000000002",
                "sequence_000000001",
                20,
                10,
                3,
                10,
                1,
                8,
                8,
                6,
                60.0,
                "plus",
            ),
        )

        matches = _create_similarity_matches(
            sequence_level=SEQUENCE_LEVEL_NUCLEOTIDE,
            selected_alleles=self.selected_alleles,
            blast_hsps=blast_hsps,
        )

        self.assertEqual(len(matches), 2)
        gene_a_query_match = next(
            match for match in matches if match.query_locus_name == "gene_a"
        )
        gene_b_query_match = next(
            match for match in matches if match.query_locus_name == "gene_b"
        )
        self.assertEqual(gene_a_query_match.identity_percentage, 75.0)
        self.assertEqual(gene_a_query_match.query_coverage_percentage, 80.0)
        self.assertEqual(gene_a_query_match.subject_coverage_percentage, 40.0)
        self.assertEqual(gene_a_query_match.minimum_coverage_percentage, 40.0)
        self.assertEqual(gene_a_query_match.blast_score_ratio, 0.6)
        self.assertEqual(gene_b_query_match.blast_score_ratio, 0.75)
        self.assertEqual(gene_a_query_match.nucleotide_subject_strand, "minus")

    def test_highest_raw_score_hsp_is_retained(self) -> None:
        """Choose raw score before identity, coverage, or BSR tie breakers."""
        blast_hsps = (
            _BlastHsp(
                "sequence_000000001",
                "sequence_000000001",
                10,
                10,
                1,
                10,
                1,
                10,
                10,
                10,
                100.0,
                "plus",
            ),
            _BlastHsp(
                "sequence_000000001",
                "sequence_000000002",
                10,
                20,
                1,
                10,
                1,
                10,
                10,
                10,
                40.0,
                "plus",
            ),
            _BlastHsp(
                "sequence_000000001",
                "sequence_000000002",
                10,
                20,
                1,
                6,
                1,
                6,
                6,
                4,
                50.0,
                "minus",
            ),
        )

        matches = _create_similarity_matches(
            sequence_level=SEQUENCE_LEVEL_NUCLEOTIDE,
            selected_alleles=self.selected_alleles,
            blast_hsps=blast_hsps,
        )

        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0].raw_score, 50.0)
        self.assertEqual(matches[0].alignment_length, 6)

    def test_amino_acid_matches_do_not_report_a_strand(self) -> None:
        """Discard BLAST's non-applicable strand placeholder for proteins."""
        selected_protein_alleles = tuple(
            _SelectedAlleleSequence(
                internal_identifier=selected_allele.internal_identifier,
                locus_name=selected_allele.locus_name,
                original_identifier=selected_allele.original_identifier,
                sequence="MK",
            )
            for selected_allele in self.selected_alleles
        )
        blast_hsps = (
            _BlastHsp(
                "sequence_000000001",
                "sequence_000000001",
                2,
                2,
                1,
                2,
                1,
                2,
                2,
                2,
                10.0,
                None,
            ),
            _BlastHsp(
                "sequence_000000001",
                "sequence_000000002",
                2,
                2,
                1,
                2,
                1,
                2,
                2,
                2,
                8.0,
                None,
            ),
        )

        matches = _create_similarity_matches(
            sequence_level=SEQUENCE_LEVEL_AMINO_ACID,
            selected_alleles=selected_protein_alleles,
            blast_hsps=blast_hsps,
        )

        self.assertIsNone(matches[0].nucleotide_subject_strand)


class CandidateParalogClassificationTest(unittest.TestCase):
    """Test advisory classification without repeating BLAST searches."""

    @staticmethod
    def _match(
        *,
        query_locus_name: str,
        subject_locus_name: str,
        sequence_level: str,
        identity_percentage: float,
        alignment_length: int,
        blast_score_ratio: float = 0.1,
    ) -> CrossLocusSimilarityMatch:
        """Create one directed similarity row with irrelevant fields fixed."""
        return CrossLocusSimilarityMatch(
            sequence_level=sequence_level,
            query_locus_name=query_locus_name,
            subject_locus_name=subject_locus_name,
            query_allele_identifier=f"{query_locus_name}_1",
            subject_allele_identifier=f"{subject_locus_name}_1",
            identity_percentage=identity_percentage,
            alignment_length=alignment_length,
            query_coverage_percentage=80.0,
            subject_coverage_percentage=75.0,
            minimum_coverage_percentage=75.0,
            raw_score=50.0,
            blast_score_ratio=blast_score_ratio,
            nucleotide_subject_strand=(
                "plus" if sequence_level == SEQUENCE_LEVEL_NUCLEOTIDE else None
            ),
        )

    @staticmethod
    def _analysis(
        *matches: CrossLocusSimilarityMatch,
    ) -> SchemeSimilarityAnalysis:
        """Wrap synthetic matches in the public analysis result type."""
        return SchemeSimilarityAnalysis(
            locus_count=4,
            nucleotide_selected_locus_count=4,
            amino_acid_selected_locus_count=4,
            nucleotide_possible_directed_pair_count=12,
            amino_acid_possible_directed_pair_count=12,
            matches=matches,
            blastn_version="blastn: test",
            blastp_version="blastp: test",
            method_identifier="test",
            blast_evalue=0.001,
            reused_from_cache=True,
        )

    def test_empty_analysis_has_no_candidate_or_evidence_pairs(self) -> None:
        """Return zero counts when BLAST reported no cross-locus matches."""
        result = classify_candidate_paralogous_locus_pairs(self._analysis())

        self.assertEqual(result.evaluated_locus_pair_count, 0)
        self.assertEqual(result.blast_score_ratio_evidence_pair_count, 0)
        self.assertEqual(result.candidate_pairs, ())

    def test_strict_amino_acid_identity_boundary_is_respected(self) -> None:
        """Require protein identity to be greater than, rather than equal to, 70%."""
        analysis = self._analysis(
            self._match(
                query_locus_name="gene_a",
                subject_locus_name="gene_b",
                sequence_level=SEQUENCE_LEVEL_AMINO_ACID,
                identity_percentage=70.0,
                alignment_length=80,
            ),
            self._match(
                query_locus_name="gene_c",
                subject_locus_name="gene_d",
                sequence_level=SEQUENCE_LEVEL_AMINO_ACID,
                identity_percentage=70.001,
                alignment_length=80,
            ),
        )

        result = classify_candidate_paralogous_locus_pairs(analysis)

        self.assertEqual(result.candidate_pair_count, 1)
        self.assertEqual(
            (
                result.candidate_pairs[0].first_locus_name,
                result.candidate_pairs[0].second_locus_name,
            ),
            ("gene_c", "gene_d"),
        )

    def test_nucleotide_identity_is_inclusive_but_length_is_strict(self) -> None:
        """Accept 90% identity only when the HSP is longer than 100 bp."""
        analysis = self._analysis(
            self._match(
                query_locus_name="gene_a",
                subject_locus_name="gene_b",
                sequence_level=SEQUENCE_LEVEL_NUCLEOTIDE,
                identity_percentage=90.0,
                alignment_length=100,
            ),
            self._match(
                query_locus_name="gene_c",
                subject_locus_name="gene_d",
                sequence_level=SEQUENCE_LEVEL_NUCLEOTIDE,
                identity_percentage=90.0,
                alignment_length=101,
            ),
        )

        result = classify_candidate_paralogous_locus_pairs(analysis)

        self.assertEqual(result.candidate_pair_count, 1)
        self.assertTrue(
            result.candidate_pairs[0].meets_nucleotide_repeated_homology_criterion
        )

    def test_reciprocal_and_sequence_level_hits_collapse_to_one_pair(self) -> None:
        """Retain strongest directed evidence but count one unordered pair."""
        analysis = self._analysis(
            self._match(
                query_locus_name="gene_a",
                subject_locus_name="gene_b",
                sequence_level=SEQUENCE_LEVEL_AMINO_ACID,
                identity_percentage=71.0,
                alignment_length=50,
                blast_score_ratio=0.7,
            ),
            self._match(
                query_locus_name="gene_b",
                subject_locus_name="gene_a",
                sequence_level=SEQUENCE_LEVEL_AMINO_ACID,
                identity_percentage=75.0,
                alignment_length=45,
                blast_score_ratio=0.8,
            ),
            self._match(
                query_locus_name="gene_a",
                subject_locus_name="gene_b",
                sequence_level=SEQUENCE_LEVEL_NUCLEOTIDE,
                identity_percentage=92.0,
                alignment_length=110,
                blast_score_ratio=0.65,
            ),
            self._match(
                query_locus_name="gene_b",
                subject_locus_name="gene_a",
                sequence_level=SEQUENCE_LEVEL_NUCLEOTIDE,
                identity_percentage=95.0,
                alignment_length=105,
                blast_score_ratio=0.6,
            ),
        )

        result = classify_candidate_paralogous_locus_pairs(analysis)
        candidate_pair = result.candidate_pairs[0]

        self.assertEqual(result.evaluated_locus_pair_count, 1)
        self.assertEqual(result.candidate_pair_count, 1)
        self.assertEqual(result.candidate_locus_count, 2)
        self.assertEqual(result.both_criteria_candidate_pair_count, 1)
        self.assertEqual(result.blast_score_ratio_evidence_pair_count, 1)
        self.assertEqual(
            candidate_pair.amino_acid_evidence_match.query_locus_name,
            "gene_b",
        )
        self.assertEqual(
            candidate_pair.nucleotide_evidence_match.query_locus_name,
            "gene_b",
        )
        self.assertEqual(
            candidate_pair.maximum_blast_score_ratio_match.sequence_level,
            SEQUENCE_LEVEL_AMINO_ACID,
        )

    def test_bsr_only_evidence_does_not_create_a_candidate(self) -> None:
        """Count the BSR reference while keeping it non-classifying."""
        analysis = self._analysis(
            self._match(
                query_locus_name="gene_a",
                subject_locus_name="gene_b",
                sequence_level=SEQUENCE_LEVEL_AMINO_ACID,
                identity_percentage=40.0,
                alignment_length=80,
                blast_score_ratio=0.6,
            )
        )

        result = classify_candidate_paralogous_locus_pairs(analysis)

        self.assertEqual(result.blast_score_ratio_evidence_pair_count, 1)
        self.assertEqual(result.candidate_pair_count, 0)

    def test_nonqualifying_sequence_level_match_remains_visible(self) -> None:
        """Distinguish a below-threshold hit from an absent BLAST hit."""
        analysis = self._analysis(
            self._match(
                query_locus_name="gene_a",
                subject_locus_name="gene_b",
                sequence_level=SEQUENCE_LEVEL_AMINO_ACID,
                identity_percentage=75.0,
                alignment_length=80,
            ),
            self._match(
                query_locus_name="gene_a",
                subject_locus_name="gene_b",
                sequence_level=SEQUENCE_LEVEL_NUCLEOTIDE,
                identity_percentage=70.61,
                alignment_length=2_215,
            ),
        )

        result = classify_candidate_paralogous_locus_pairs(analysis)
        candidate_pair = result.candidate_pairs[0]

        self.assertTrue(candidate_pair.meets_amino_acid_identity_criterion)
        self.assertFalse(
            candidate_pair.meets_nucleotide_repeated_homology_criterion
        )
        self.assertEqual(
            candidate_pair.nucleotide_evidence_match.identity_percentage,
            70.61,
        )
        self.assertEqual(
            candidate_pair.nucleotide_evidence_match.alignment_length,
            2_215,
        )

    def test_configurable_thresholds_reclassify_existing_measurements(self) -> None:
        """Apply new thresholds to the same immutable BLAST measurements."""
        analysis = self._analysis(
            self._match(
                query_locus_name="gene_a",
                subject_locus_name="gene_b",
                sequence_level=SEQUENCE_LEVEL_AMINO_ACID,
                identity_percentage=75.0,
                alignment_length=80,
            )
        )

        default_result = classify_candidate_paralogous_locus_pairs(analysis)
        stricter_result = classify_candidate_paralogous_locus_pairs(
            analysis,
            thresholds=CandidateParalogThresholds(
                amino_acid_identity_threshold_percentage=80.0
            ),
        )

        self.assertEqual(default_result.candidate_pair_count, 1)
        self.assertEqual(stricter_result.candidate_pair_count, 0)
        self.assertEqual(
            stricter_result.thresholds.amino_acid_identity_threshold_percentage,
            80.0,
        )


class SimilarityCacheTest(unittest.TestCase):
    """Test whole-analysis cache reuse and source invalidation."""

    @staticmethod
    def _empty_level_result(
        *,
        sequence_level: str,
        **_: object,
    ) -> _SimilarityLevelResult:
        """Return a deterministic no-hit result without launching BLAST."""
        return _SimilarityLevelResult(sequence_level=sequence_level, matches=())

    def test_cache_is_reused_and_invalidated_by_fasta_change(self) -> None:
        """Avoid repeated searches until a source FASTA signature changes."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            first_locus_path = scheme_directory_path / "gene_a.fasta"
            first_locus_path.write_text(
                ">gene_a_1\nATGAAATAA\n",
                encoding="utf-8",
            )
            (scheme_directory_path / "gene_b.fasta").write_text(
                ">gene_b_1\nATGCCCTAA\n",
                encoding="utf-8",
            )
            progress_updates = []

            with (
                patch.object(
                    similarity,
                    "_get_blast_version",
                    side_effect=lambda executable: f"{executable}: 1.0",
                ),
                patch.object(
                    similarity,
                    "_run_similarity_level",
                    side_effect=self._empty_level_result,
                ) as run_similarity_level,
            ):
                first_analysis = get_scheme_similarity_analysis(
                    scheme_directory_path,
                    progress_callback=progress_updates.append,
                )
                cached_analysis = get_scheme_similarity_analysis(
                    scheme_directory_path,
                    progress_callback=progress_updates.append,
                )
                first_locus_path.write_text(
                    ">gene_a_1\nATGAAATAA\n>gene_a_2\nATGGGGTAA\n",
                    encoding="utf-8",
                )
                recalculated_analysis = get_scheme_similarity_analysis(
                    scheme_directory_path,
                )

            self.assertEqual(run_similarity_level.call_count, 4)
            self.assertFalse(first_analysis.reused_from_cache)
            self.assertTrue(cached_analysis.reused_from_cache)
            self.assertFalse(recalculated_analysis.reused_from_cache)
            self.assertTrue(
                (scheme_directory_path / SCHEME_SIMILARITY_CACHE_FILENAME).is_file()
            )
            self.assertTrue(
                any(progress.reused_from_cache for progress in progress_updates)
            )

    def test_empty_scheme_is_reported_without_locus_pairs(self) -> None:
        """Handle a scheme without any recognized locus FASTA files."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            with (
                patch.object(
                    similarity,
                    "_get_blast_version",
                    side_effect=lambda executable: f"{executable}: 1.0",
                ),
                patch.object(
                    similarity,
                    "_run_similarity_level",
                    side_effect=self._empty_level_result,
                ),
            ):
                analysis = get_scheme_similarity_analysis(scheme_directory_path)

        self.assertEqual(analysis.locus_count, 0)
        self.assertEqual(analysis.matches, ())
        self.assertEqual(analysis.nucleotide_possible_directed_pair_count, 0)
        self.assertEqual(analysis.amino_acid_possible_directed_pair_count, 0)

    def test_complete_cache_does_not_require_installed_blast(self) -> None:
        """Keep cached measurements available outside a BLAST-enabled shell."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "gene_a.fasta").write_text(
                ">gene_a_1\nATGAAATAA\n",
                encoding="utf-8",
            )
            with (
                patch.object(
                    similarity,
                    "_get_blast_version",
                    side_effect=lambda executable: f"{executable}: 1.0",
                ),
                patch.object(
                    similarity,
                    "_run_similarity_level",
                    side_effect=self._empty_level_result,
                ),
            ):
                get_scheme_similarity_analysis(scheme_directory_path)
            legacy_cache_file_path = (
                scheme_directory_path / SCHEME_SIMILARITY_CACHE_FILENAMES[1]
            )
            (scheme_directory_path / SCHEME_SIMILARITY_CACHE_FILENAME).replace(
                legacy_cache_file_path
            )

            with (
                patch.object(
                    similarity,
                    "_get_blast_version",
                    side_effect=RuntimeError("BLAST unavailable"),
                ),
                patch.object(
                    similarity,
                    "_run_similarity_level",
                ) as run_similarity_level,
            ):
                cached_analysis = get_scheme_similarity_analysis(
                    scheme_directory_path
                )

        self.assertTrue(cached_analysis.reused_from_cache)
        run_similarity_level.assert_not_called()


@unittest.skipUnless(
    shutil.which("blastn")
    and shutil.which("blastp")
    and shutil.which("makeblastdb"),
    "NCBI BLAST+ is not available",
)
class BlastIntegrationTest(unittest.TestCase):
    """Exercise the public API with real nucleotide and protein BLAST runs."""

    def test_tiny_scheme_similarity_and_cache_reuse(self) -> None:
        """Find reciprocal cross-locus hits and then load them from cache."""
        diverse_coding_codons = [
            "GCT",
            "TTT",
            "GAA",
            "AAA",
            "CCG",
            "GTT",
            "TGG",
            "CAT",
            "AAT",
            "GAC",
        ] * 4
        variant_coding_codons = diverse_coding_codons.copy()
        variant_coding_codons[21] = "TAT"
        first_sequence = f"ATG{''.join(diverse_coding_codons)}TAA"
        second_sequence = f"ATG{''.join(variant_coding_codons)}TAA"
        with tempfile.TemporaryDirectory() as temporary_directory:
            scheme_directory_path = Path(temporary_directory)
            (scheme_directory_path / "gene_a.fasta").write_text(
                f">gene_a_1\n{first_sequence}\n",
                encoding="utf-8",
            )
            (scheme_directory_path / "gene_b.fasta").write_text(
                f">gene_b_1\n{second_sequence}\n",
                encoding="utf-8",
            )

            analysis = get_scheme_similarity_analysis(
                scheme_directory_path,
                maximum_threads=2,
            )
            cached_analysis = get_scheme_similarity_analysis(
                scheme_directory_path,
                maximum_threads=1,
            )

        self.assertEqual(analysis.nucleotide_match_count, 2)
        self.assertEqual(analysis.amino_acid_match_count, 2)
        self.assertTrue(
            all(match.minimum_coverage_percentage > 95 for match in analysis.matches)
        )
        self.assertTrue(cached_analysis.reused_from_cache)


if __name__ == "__main__":
    unittest.main()
