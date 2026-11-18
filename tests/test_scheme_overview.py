"""Tests for progressive intrinsic scheme-overview aggregation."""

import unittest

from scheme.alignment import (
    ALIGNMENT_STATUS_ALIGNED,
    ALIGNMENT_STATUS_INSUFFICIENT_SEQUENCES,
    LocusAlignmentStatistics,
    SchemeAlignmentAnalysis,
    SchemeAlignmentSummary,
)
from scheme.allele_info import LocusAlleleInformation
from scheme.bsr import (
    RepresentativeBsrClassification,
    SchemeRepresentativeBsrAnalysis,
)
from scheme.info import SchemeInformation
from scheme.overview import create_intrinsic_scheme_overview
from scheme.similarity import (
    CandidateParalogAnalysis,
    CandidateParalogThresholds,
    SchemeSimilarityAnalysis,
)
from scheme.stats import AlleleLengthPeak, LocusStatistics


def _scheme_information() -> SchemeInformation:
    """Return a compact consistent scheme-wide fixture."""
    return SchemeInformation(
        locus_count=3,
        allele_count=6,
        valid_allele_count=4,
        invalid_allele_count=2,
        incomplete_orf_allele_count=1,
        ambiguous_bases_allele_count=1,
        missing_start_or_stop_codon_allele_count=0,
        in_frame_stop_codon_allele_count=0,
        alleles_below_minimum_length_count=0,
        minimum_sequence_length=0,
        alleles_below_locus_mode_threshold_count=1,
        alleles_above_locus_mode_threshold_count=1,
        alleles_below_locus_mode_threshold_percentage=100 / 6,
        alleles_above_locus_mode_threshold_percentage=100 / 6,
        guanine_cytosine_content_percentage=50.0,
        canonical_nucleotide_base_count=60,
        excluded_nucleotide_base_count=2,
    )


def _locus_statistics(
    locus_name: str,
    *,
    allele_count: int,
    frequent_peak_count: int,
    has_outlier: bool = False,
) -> LocusStatistics:
    """Return one locus-statistics fixture with selected overview traits."""
    peaks = tuple(
        AlleleLengthPeak(allele_length=90 + index, allele_count=1, allele_percentage=50)
        for index in range(frequent_peak_count)
    )
    return LocusStatistics(
        locus_name=locus_name,
        allele_count=allele_count,
        modal_allele_length=90 if allele_count else None,
        modal_allele_lengths=(90,) if allele_count else (),
        mean_allele_length=90.0 if allele_count else None,
        frequent_allele_length_peaks=peaks,
        alleles_below_mode_threshold_count=int(has_outlier),
        alleles_above_mode_threshold_count=0,
        length_distribution_category="fixture",
        minimum_allele_length=90 if allele_count else None,
        first_quartile_allele_length=90.0 if allele_count else None,
        median_allele_length=90.0 if allele_count else None,
        third_quartile_allele_length=90.0 if allele_count else None,
        maximum_allele_length=90 if allele_count else None,
        lower_whisker_allele_length=90 if allele_count else None,
        upper_whisker_allele_length=90 if allele_count else None,
        guanine_cytosine_content_percentage=50.0 if allele_count else None,
        canonical_nucleotide_base_count=allele_count * 90,
        excluded_nucleotide_base_count=0,
    )


def _allele_information(
    locus_name: str,
    *,
    total_count: int,
    valid_count: int,
    invalid_count: int,
    short_count: int,
) -> LocusAlleleInformation:
    """Return one per-locus allele-information fixture."""
    return LocusAlleleInformation(
        locus_name=locus_name,
        total_allele_count=total_count,
        valid_allele_count=valid_count,
        invalid_allele_count=invalid_count,
        validated_allele_proportion=(
            valid_count / total_count if total_count else 0.0
        ),
        distinct_protein_allele_count=valid_count,
        incomplete_orf_allele_count=invalid_count,
        ambiguous_bases_allele_count=0,
        missing_start_or_stop_codon_allele_count=0,
        in_frame_stop_codon_allele_count=0,
        alleles_below_minimum_length_count=short_count,
        minimum_sequence_length=80,
        alleles_below_mode_threshold_count=0,
        alleles_above_mode_threshold_count=0,
    )


def _alignment_locus(
    locus_name: str,
    *,
    callable_sites: int | None,
    polymorphic_sites: int | None,
) -> LocusAlignmentStatistics:
    """Return one aligned or unavailable locus fixture."""
    is_aligned = callable_sites is not None
    return LocusAlignmentStatistics(
        locus_name=locus_name,
        total_allele_count=2,
        nucleotide_aligned_sequence_count=2 if is_aligned else 1,
        nucleotide_excluded_sequence_count=0,
        nucleotide_alignment_status=(
            ALIGNMENT_STATUS_ALIGNED
            if is_aligned
            else ALIGNMENT_STATUS_INSUFFICIENT_SEQUENCES
        ),
        nucleotide_alignment_length=10 if is_aligned else None,
        nucleotide_gap_character_count=0,
        nucleotide_gap_proportion_percentage=0.0 if is_aligned else None,
        nucleotide_callable_site_count=callable_sites,
        nucleotide_polymorphic_site_count=polymorphic_sites,
        nucleotide_polymorphic_site_percentage=(
            polymorphic_sites / callable_sites * 100
            if callable_sites and polymorphic_sites is not None
            else None
        ),
        nucleotide_alignment_error=None,
        amino_acid_aligned_sequence_count=2 if is_aligned else 1,
        amino_acid_excluded_sequence_count=0,
        amino_acid_alignment_status=(
            ALIGNMENT_STATUS_ALIGNED
            if is_aligned
            else ALIGNMENT_STATUS_INSUFFICIENT_SEQUENCES
        ),
        amino_acid_alignment_length=3 if is_aligned else None,
        amino_acid_gap_character_count=0,
        amino_acid_gap_proportion_percentage=0.0 if is_aligned else None,
        amino_acid_alignment_error=None,
    )


class IntrinsicSchemeOverviewTest(unittest.TestCase):
    """Protect static aggregation and optional expensive-stage summaries."""

    def setUp(self) -> None:
        """Create shared three-locus static fixtures."""
        self.locus_statistics = (
            _locus_statistics(
                "single",
                allele_count=2,
                frequent_peak_count=1,
            ),
            _locus_statistics(
                "multi",
                allele_count=4,
                frequent_peak_count=2,
                has_outlier=True,
            ),
            _locus_statistics(
                "empty",
                allele_count=0,
                frequent_peak_count=0,
            ),
        )
        self.allele_information = (
            _allele_information(
                "single",
                total_count=2,
                valid_count=2,
                invalid_count=0,
                short_count=0,
            ),
            _allele_information(
                "multi",
                total_count=4,
                valid_count=2,
                invalid_count=2,
                short_count=1,
            ),
            _allele_information(
                "empty",
                total_count=0,
                valid_count=0,
                invalid_count=0,
                short_count=0,
            ),
        )

    def test_static_overview_is_available_before_expensive_stages(self) -> None:
        """Report automatic metrics while leaving requested stages pending."""
        overview = create_intrinsic_scheme_overview(
            scheme_information=_scheme_information(),
            locus_statistics=self.locus_statistics,
            scheme_allele_information=self.allele_information,
            minimum_sequence_length=80,
            alleles_below_minimum_length_count=1,
            frequent_peak_relative_threshold_percentage=50.0,
        )

        self.assertEqual(overview.static.empty_locus_count, 1)
        self.assertEqual(overview.static.single_mode_locus_count, 1)
        self.assertEqual(overview.static.multimodal_locus_count, 1)
        self.assertEqual(overview.static.length_outlier_locus_count, 1)
        self.assertEqual(overview.static.invalid_allele_locus_count, 1)
        self.assertEqual(overview.static.no_valid_allele_locus_count, 1)
        self.assertEqual(overview.static.short_allele_locus_count, 1)
        self.assertIsNone(overview.alignment)
        self.assertIsNone(overview.similarity)
        self.assertIsNone(overview.representative_bsr)

    def test_completed_stages_add_pooled_and_classified_metrics(self) -> None:
        """Attach optional stages without recalculating their measurements."""
        alignment_analysis = SchemeAlignmentAnalysis(
            mafft_version="MAFFT fixture",
            alignment_method_identifier="fixture-method",
            locus_statistics=(
                _alignment_locus(
                    "first",
                    callable_sites=10,
                    polymorphic_sites=2,
                ),
                _alignment_locus(
                    "second",
                    callable_sites=30,
                    polymorphic_sites=1,
                ),
                _alignment_locus(
                    "third",
                    callable_sites=None,
                    polymorphic_sites=None,
                ),
            ),
        )
        alignment_summary = SchemeAlignmentSummary(
            locus_count=3,
            nucleotide_aligned_locus_count=2,
            nucleotide_unavailable_locus_count=1,
            nucleotide_failed_locus_count=0,
            amino_acid_aligned_locus_count=2,
            amino_acid_unavailable_locus_count=1,
            amino_acid_failed_locus_count=0,
            amino_acid_gap_warning_threshold_percentage=5.0,
            amino_acid_gap_warning_locus_count=0,
            amino_acid_gap_warning_locus_percentage=0.0,
            amino_acid_gap_warning_loci=(),
        )
        similarity_analysis = SchemeSimilarityAnalysis(
            locus_count=3,
            nucleotide_selected_locus_count=3,
            amino_acid_selected_locus_count=2,
            nucleotide_possible_directed_pair_count=6,
            amino_acid_possible_directed_pair_count=2,
            matches=(),
            blastn_version="blastn fixture",
            blastp_version="blastp fixture",
            method_identifier="fixture-method",
            blast_evalue=0.001,
            reused_from_cache=True,
        )
        candidate_analysis = CandidateParalogAnalysis(
            thresholds=CandidateParalogThresholds(),
            evaluated_locus_pair_count=0,
            blast_score_ratio_evidence_pair_count=0,
            candidate_pairs=(),
        )
        representative_analysis = SchemeRepresentativeBsrAnalysis(
            locus_measurements=(),
            inter_locus_measurements=(),
            blastp_version="blastp fixture",
            method_identifier="fixture-method",
            blast_evalue=0.001,
            reused_from_cache=True,
        )
        representative_classification = RepresentativeBsrClassification(
            blast_score_ratio_threshold=0.6,
            locus_coverages=(),
            high_bsr_pairs=(),
            high_bsr_groups=(),
        )

        overview = create_intrinsic_scheme_overview(
            scheme_information=_scheme_information(),
            locus_statistics=self.locus_statistics,
            scheme_allele_information=self.allele_information,
            minimum_sequence_length=80,
            alleles_below_minimum_length_count=1,
            frequent_peak_relative_threshold_percentage=50.0,
            scheme_alignment_analysis=alignment_analysis,
            scheme_alignment_summary=alignment_summary,
            scheme_similarity_analysis=similarity_analysis,
            candidate_paralog_analysis=candidate_analysis,
            scheme_representative_bsr_analysis=representative_analysis,
            representative_bsr_classification=(
                representative_classification
            ),
        )

        self.assertIsNotNone(overview.alignment)
        assert overview.alignment is not None
        self.assertEqual(overview.alignment.nucleotide_callable_site_count, 40)
        self.assertEqual(overview.alignment.nucleotide_polymorphic_site_count, 3)
        self.assertEqual(
            overview.alignment.nucleotide_polymorphic_site_percentage,
            7.5,
        )
        self.assertEqual(overview.alignment.polymorphic_locus_count, 2)
        self.assertIs(
            overview.similarity.candidate_paralog_analysis,
            candidate_analysis,
        )
        self.assertIs(
            overview.representative_bsr.classification,
            representative_classification,
        )


if __name__ == "__main__":
    unittest.main()
