"""Aggregate existing intrinsic-analysis results for the notebook overview.

This module performs no sequence parsing and invokes no external tools. It
combines already calculated typed results into compact scheme-wide summaries.
Expensive analysis summaries are optional so the Marimo overview can render
static measurements immediately and update progressively after requested
MAFFT or BLAST stages finish.
"""

from dataclasses import dataclass

from .alignment import SchemeAlignmentAnalysis, SchemeAlignmentSummary
from .allele_info import LocusAlleleInformation
from .bsr import (
    RepresentativeBsrClassification,
    SchemeRepresentativeBsrAnalysis,
)
from .info import SchemeInformation
from .similarity import CandidateParalogAnalysis, SchemeSimilarityAnalysis
from .stats import LocusStatistics

__all__ = [
    "AlignmentOverview",
    "IntrinsicSchemeOverview",
    "RepresentativeBsrOverview",
    "SimilarityOverview",
    "StaticSchemeOverview",
    "create_intrinsic_scheme_overview",
]


@dataclass(frozen=True)
class StaticSchemeOverview:
    """Automatically calculated scheme-wide and per-locus measurements."""

    scheme_information: SchemeInformation
    minimum_sequence_length: int
    alleles_below_minimum_length_count: int
    frequent_peak_relative_threshold_percentage: float
    empty_locus_count: int
    single_mode_locus_count: int
    multimodal_locus_count: int
    length_outlier_locus_count: int
    invalid_allele_locus_count: int
    no_valid_allele_locus_count: int
    short_allele_locus_count: int


@dataclass(frozen=True)
class AlignmentOverview:
    """Scheme-wide alignment availability and polymorphism measurements."""

    analysis: SchemeAlignmentAnalysis
    summary: SchemeAlignmentSummary
    nucleotide_callable_site_count: int
    nucleotide_polymorphic_site_count: int
    nucleotide_polymorphic_site_percentage: float | None
    polymorphic_locus_count: int


@dataclass(frozen=True)
class SimilarityOverview:
    """Selected-allele similarity and threshold-classified candidate results."""

    analysis: SchemeSimilarityAnalysis
    candidate_paralog_analysis: CandidateParalogAnalysis | None


@dataclass(frozen=True)
class RepresentativeBsrOverview:
    """Representative coverage and high-BSR inter-locus evidence."""

    analysis: SchemeRepresentativeBsrAnalysis
    classification: RepresentativeBsrClassification | None


@dataclass(frozen=True)
class IntrinsicSchemeOverview:
    """Complete progressive overview for one intrinsic scheme evaluation."""

    static: StaticSchemeOverview
    alignment: AlignmentOverview | None
    similarity: SimilarityOverview | None
    representative_bsr: RepresentativeBsrOverview | None


def _create_static_overview(
    *,
    scheme_information: SchemeInformation,
    locus_statistics: tuple[LocusStatistics, ...],
    scheme_allele_information: tuple[LocusAlleleInformation, ...],
    minimum_sequence_length: int,
    alleles_below_minimum_length_count: int,
    frequent_peak_relative_threshold_percentage: float,
) -> StaticSchemeOverview:
    """Aggregate automatic static measurements without rereading FASTA files."""
    return StaticSchemeOverview(
        scheme_information=scheme_information,
        minimum_sequence_length=minimum_sequence_length,
        alleles_below_minimum_length_count=(
            alleles_below_minimum_length_count
        ),
        frequent_peak_relative_threshold_percentage=(
            frequent_peak_relative_threshold_percentage
        ),
        empty_locus_count=sum(
            locus.allele_count == 0 for locus in locus_statistics
        ),
        single_mode_locus_count=sum(
            locus.allele_count > 0
            and len(locus.frequent_allele_length_peaks) == 1
            for locus in locus_statistics
        ),
        multimodal_locus_count=sum(
            len(locus.frequent_allele_length_peaks) > 1
            for locus in locus_statistics
        ),
        length_outlier_locus_count=sum(
            locus.alleles_below_mode_threshold_count
            + locus.alleles_above_mode_threshold_count
            > 0
            for locus in locus_statistics
        ),
        invalid_allele_locus_count=sum(
            locus.invalid_allele_count > 0
            for locus in scheme_allele_information
        ),
        no_valid_allele_locus_count=sum(
            locus.valid_allele_count == 0
            for locus in scheme_allele_information
        ),
        short_allele_locus_count=sum(
            locus.alleles_below_minimum_length_count > 0
            for locus in scheme_allele_information
        ),
    )


def _create_alignment_overview(
    analysis: SchemeAlignmentAnalysis,
    summary: SchemeAlignmentSummary,
) -> AlignmentOverview:
    """Pool callable and polymorphic nucleotide sites across aligned loci."""
    nucleotide_callable_site_count = sum(
        locus.nucleotide_callable_site_count or 0
        for locus in analysis.locus_statistics
    )
    nucleotide_polymorphic_site_count = sum(
        locus.nucleotide_polymorphic_site_count or 0
        for locus in analysis.locus_statistics
    )
    nucleotide_polymorphic_site_percentage = (
        nucleotide_polymorphic_site_count
        / nucleotide_callable_site_count
        * 100
        if nucleotide_callable_site_count
        else None
    )
    return AlignmentOverview(
        analysis=analysis,
        summary=summary,
        nucleotide_callable_site_count=nucleotide_callable_site_count,
        nucleotide_polymorphic_site_count=(
            nucleotide_polymorphic_site_count
        ),
        nucleotide_polymorphic_site_percentage=(
            nucleotide_polymorphic_site_percentage
        ),
        polymorphic_locus_count=sum(
            (locus.nucleotide_polymorphic_site_count or 0) > 0
            for locus in analysis.locus_statistics
        ),
    )


def create_intrinsic_scheme_overview(
    *,
    scheme_information: SchemeInformation,
    locus_statistics: tuple[LocusStatistics, ...],
    scheme_allele_information: tuple[LocusAlleleInformation, ...],
    minimum_sequence_length: int,
    alleles_below_minimum_length_count: int,
    frequent_peak_relative_threshold_percentage: float,
    scheme_alignment_analysis: SchemeAlignmentAnalysis | None = None,
    scheme_alignment_summary: SchemeAlignmentSummary | None = None,
    scheme_similarity_analysis: SchemeSimilarityAnalysis | None = None,
    candidate_paralog_analysis: CandidateParalogAnalysis | None = None,
    scheme_representative_bsr_analysis: (
        SchemeRepresentativeBsrAnalysis | None
    ) = None,
    representative_bsr_classification: (
        RepresentativeBsrClassification | None
    ) = None,
) -> IntrinsicSchemeOverview:
    """Create a progressive overview from existing analysis results.

    An expensive stage is unavailable until its raw analysis exists. Its
    threshold-dependent classification may briefly remain unavailable while a
    reactive threshold cell updates; raw measurements are still retained.
    """
    alignment_overview = (
        _create_alignment_overview(
            scheme_alignment_analysis,
            scheme_alignment_summary,
        )
        if scheme_alignment_analysis is not None
        and scheme_alignment_summary is not None
        else None
    )
    similarity_overview = (
        SimilarityOverview(
            analysis=scheme_similarity_analysis,
            candidate_paralog_analysis=candidate_paralog_analysis,
        )
        if scheme_similarity_analysis is not None
        else None
    )
    representative_bsr_overview = (
        RepresentativeBsrOverview(
            analysis=scheme_representative_bsr_analysis,
            classification=representative_bsr_classification,
        )
        if scheme_representative_bsr_analysis is not None
        else None
    )
    return IntrinsicSchemeOverview(
        static=_create_static_overview(
            scheme_information=scheme_information,
            locus_statistics=locus_statistics,
            scheme_allele_information=scheme_allele_information,
            minimum_sequence_length=minimum_sequence_length,
            alleles_below_minimum_length_count=(
                alleles_below_minimum_length_count
            ),
            frequent_peak_relative_threshold_percentage=(
                frequent_peak_relative_threshold_percentage
            ),
        ),
        alignment=alignment_overview,
        similarity=similarity_overview,
        representative_bsr=representative_bsr_overview,
    )
