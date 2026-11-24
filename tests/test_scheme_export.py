"""Tests for FAIR-oriented intrinsic evaluation bundles."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import tempfile
import unittest
import zipfile
from datetime import UTC, datetime
from pathlib import Path

from scheme.alignment import SchemeAlignmentAnalysis, SchemeAlignmentSummary
from scheme.allele_info import LocusAlleleInformation
from scheme.bsr import (
    RepresentativeBsrClassification,
    SchemeRepresentativeBsrAnalysis,
)
from scheme.export import (
    IntrinsicOperationalSettings,
    create_intrinsic_evaluation_bundle,
    intrinsic_evaluation_bundle_filename,
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
    return SchemeInformation(
        locus_count=1,
        allele_count=2,
        valid_allele_count=2,
        invalid_allele_count=0,
        incomplete_orf_allele_count=0,
        ambiguous_bases_allele_count=0,
        missing_start_or_stop_codon_allele_count=0,
        in_frame_stop_codon_allele_count=0,
        alleles_below_minimum_length_count=0,
        minimum_sequence_length=50,
        alleles_below_locus_mode_threshold_count=0,
        alleles_above_locus_mode_threshold_count=0,
        alleles_below_locus_mode_threshold_percentage=0.0,
        alleles_above_locus_mode_threshold_percentage=0.0,
        guanine_cytosine_content_percentage=50.0,
        canonical_nucleotide_base_count=18,
        excluded_nucleotide_base_count=0,
    )


def _locus_statistics() -> tuple[LocusStatistics, ...]:
    return (
        LocusStatistics(
            locus_name="locus_a",
            allele_count=2,
            modal_allele_length=9,
            modal_allele_lengths=(9,),
            mean_allele_length=9.0,
            frequent_allele_length_peaks=(
                AlleleLengthPeak(9, 2, 100.0),
            ),
            alleles_below_mode_threshold_count=0,
            alleles_above_mode_threshold_count=0,
            length_distribution_category="Single-mode",
            minimum_allele_length=9,
            first_quartile_allele_length=9.0,
            median_allele_length=9.0,
            third_quartile_allele_length=9.0,
            maximum_allele_length=9,
            lower_whisker_allele_length=9,
            upper_whisker_allele_length=9,
            guanine_cytosine_content_percentage=None,
            canonical_nucleotide_base_count=18,
            excluded_nucleotide_base_count=0,
        ),
    )


def _allele_information() -> tuple[LocusAlleleInformation, ...]:
    return (
        LocusAlleleInformation(
            locus_name="locus_a",
            total_allele_count=2,
            valid_allele_count=2,
            invalid_allele_count=0,
            validated_allele_proportion=1.0,
            distinct_protein_allele_count=1,
            incomplete_orf_allele_count=0,
            ambiguous_bases_allele_count=0,
            missing_start_or_stop_codon_allele_count=0,
            in_frame_stop_codon_allele_count=0,
            alleles_below_minimum_length_count=0,
            minimum_sequence_length=50,
            alleles_below_mode_threshold_count=0,
            alleles_above_mode_threshold_count=0,
        ),
    )


class IntrinsicEvaluationBundleTest(unittest.TestCase):
    """Protect bundle contents, stage status, provenance, and checksums."""

    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.repository_root = Path(self.temporary_directory.name)
        self.archive_file_path = (
            self.repository_root / "archive" / "mlst" / "fixture.tar.zst"
        )
        self.archive_file_path.parent.mkdir(parents=True)
        self.archive_file_path.write_bytes(b"fixture archive bytes")
        self.inventory_file_path = (
            self.repository_root / "metadata" / "scheme_sources.toml"
        )
        self.inventory_file_path.parent.mkdir()
        self.inventory_file_path.write_text(
            """
schema_version = 1
provenance_status = "incomplete"

[[schemes]]
scheme_id = "fixture-scheme"
archive_path = "archive/mlst/fixture.tar.zst"
source_provider = "TODO"
archive_sha256 = "TODO"
notes = "Historical source details are not yet known."
""".strip()
            + "\n",
            encoding="utf-8",
        )
        self.locus_statistics = _locus_statistics()
        self.allele_information = _allele_information()
        self.scheme_information = _scheme_information()
        self.operational_settings = IntrinsicOperationalSettings(
            detected_available_cpu_count=12,
            submitted_worker_budget=8,
            selected_mafft_strategy="fft-ns-2",
            amino_acid_gap_warning_threshold_percentage=5.0,
            candidate_paralog_thresholds=CandidateParalogThresholds(),
            representative_bsr_threshold=0.6,
        )

    def _create_overview(self, **optional_results: object):
        return create_intrinsic_scheme_overview(
            scheme_information=self.scheme_information,
            locus_statistics=self.locus_statistics,
            scheme_allele_information=self.allele_information,
            minimum_sequence_length=50,
            alleles_below_minimum_length_count=0,
            frequent_peak_relative_threshold_percentage=50.0,
            **optional_results,
        )

    def _create_bundle(self, overview, **optional_results: object) -> bytes:
        return create_intrinsic_evaluation_bundle(
            repository_root_directory_path=self.repository_root,
            archive_file_path=self.archive_file_path,
            scheme_source_inventory_file_path=self.inventory_file_path,
            operational_settings=self.operational_settings,
            locus_statistics=self.locus_statistics,
            scheme_allele_information=self.allele_information,
            intrinsic_scheme_overview=overview,
            scheme_alignment_analysis=optional_results.get("alignment"),
            scheme_alignment_summary=optional_results.get("alignment_summary"),
            scheme_similarity_analysis=optional_results.get("similarity"),
            candidate_paralog_analysis=optional_results.get("candidates"),
            scheme_representative_bsr_analysis=optional_results.get(
                "representative_bsr"
            ),
            representative_bsr_classification=optional_results.get(
                "representative_classification"
            ),
            created_at=datetime(2026, 11, 19, 12, 0, tzinfo=UTC),
        )

    def test_static_bundle_is_explicitly_partial_and_self_verifying(self) -> None:
        bundle = self._create_bundle(self._create_overview())

        with zipfile.ZipFile(io.BytesIO(bundle)) as bundle_zip:
            file_names = set(bundle_zip.namelist())
            self.assertEqual(
                file_names,
                {
                    "checksums.sha256",
                    "input_provenance.json",
                    "locus_allele_information.tsv",
                    "locus_statistics.tsv",
                    "manifest.json",
                    "software_versions.json",
                    "summary.json",
                },
            )
            manifest = json.loads(bundle_zip.read("manifest.json"))
            self.assertEqual(manifest["stage_status"]["static"], "complete")
            self.assertEqual(manifest["stage_status"]["alignment"], "not_run")
            self.assertEqual(manifest["software_license"], "EUPL-1.2")
            self.assertEqual(manifest["missing_value_representation"]["json"], "null")
            operational_provenance = manifest["operational_provenance"]
            self.assertEqual(
                operational_provenance["submitted_settings"]["worker_budget"],
                8,
            )
            self.assertEqual(
                operational_provenance["submitted_settings"]["mafft_strategy"],
                "fft-ns-2",
            )
            self.assertEqual(
                operational_provenance["input_identities"]["scheme_id"],
                "fixture-scheme",
            )
            self.assertEqual(
                operational_provenance["input_identities"][
                    "scheme_source_inventory_path"
                ],
                "metadata/scheme_sources.toml",
            )
            self.assertFalse(
                operational_provenance["cache_key_policy"][
                    "worker_budget_is_scientific_cache_key"
                ]
            )

            provenance = json.loads(bundle_zip.read("input_provenance.json"))
            self.assertFalse(provenance["provenance_complete"])
            self.assertEqual(
                provenance["scheme_source_placeholder_fields"],
                ["archive_sha256", "source_provider"],
            )
            self.assertEqual(
                provenance["archive_sha256"],
                hashlib.sha256(b"fixture archive bytes").hexdigest(),
            )
            self.assertIsNone(provenance["recorded_archive_sha256_matches"])

            statistics_rows = tuple(
                csv.DictReader(
                    io.StringIO(
                        bundle_zip.read("locus_statistics.tsv").decode("utf-8")
                    ),
                    dialect="excel-tab",
                )
            )
            self.assertEqual(
                statistics_rows[0]["guanine_cytosine_content_percentage"],
                "unavailable",
            )
            self.assertEqual(
                json.loads(statistics_rows[0]["frequent_allele_length_peaks"]),
                [{"allele_length": 9, "allele_count": 2, "allele_percentage": 100.0}],
            )

            checksum_lines = bundle_zip.read("checksums.sha256").decode().splitlines()
            recorded_checksums = {
                line.split("  ", maxsplit=1)[1]: line.split("  ", maxsplit=1)[0]
                for line in checksum_lines
            }
            self.assertEqual(
                set(recorded_checksums),
                file_names - {"checksums.sha256"},
            )
            for file_name, recorded_checksum in recorded_checksums.items():
                self.assertEqual(
                    recorded_checksum,
                    hashlib.sha256(bundle_zip.read(file_name)).hexdigest(),
                )

    def test_complete_results_emit_every_optional_table_and_tool_provenance(self) -> None:
        alignment = SchemeAlignmentAnalysis(
            mafft_version="MAFFT fixture",
            alignment_method_identifier="mafft-fixture",
            locus_statistics=(),
        )
        alignment_summary = SchemeAlignmentSummary(
            locus_count=0,
            nucleotide_aligned_locus_count=0,
            nucleotide_unavailable_locus_count=0,
            nucleotide_failed_locus_count=0,
            amino_acid_aligned_locus_count=0,
            amino_acid_unavailable_locus_count=0,
            amino_acid_failed_locus_count=0,
            amino_acid_gap_warning_threshold_percentage=5.0,
            amino_acid_gap_warning_locus_count=0,
            amino_acid_gap_warning_locus_percentage=0.0,
            amino_acid_gap_warning_loci=(),
        )
        similarity = SchemeSimilarityAnalysis(
            locus_count=0,
            nucleotide_selected_locus_count=0,
            amino_acid_selected_locus_count=0,
            nucleotide_possible_directed_pair_count=0,
            amino_acid_possible_directed_pair_count=0,
            matches=(),
            blastn_version="blastn fixture",
            blastp_version="blastp fixture",
            method_identifier="similarity-fixture",
            blast_evalue=1e-3,
            reused_from_cache=True,
        )
        candidates = CandidateParalogAnalysis(
            thresholds=CandidateParalogThresholds(),
            evaluated_locus_pair_count=0,
            blast_score_ratio_evidence_pair_count=0,
            candidate_pairs=(),
        )
        representative_bsr = SchemeRepresentativeBsrAnalysis(
            locus_measurements=(),
            inter_locus_measurements=(),
            blastp_version="blastp fixture",
            method_identifier="representative-fixture",
            blast_evalue=1e-3,
            reused_from_cache=False,
        )
        representative_classification = RepresentativeBsrClassification(
            blast_score_ratio_threshold=0.6,
            locus_coverages=(),
            high_bsr_pairs=(),
            high_bsr_groups=(),
        )
        overview = self._create_overview(
            scheme_alignment_analysis=alignment,
            scheme_alignment_summary=alignment_summary,
            scheme_similarity_analysis=similarity,
            candidate_paralog_analysis=candidates,
            scheme_representative_bsr_analysis=representative_bsr,
            representative_bsr_classification=representative_classification,
        )
        bundle = self._create_bundle(
            overview,
            alignment=alignment,
            alignment_summary=alignment_summary,
            similarity=similarity,
            candidates=candidates,
            representative_bsr=representative_bsr,
            representative_classification=representative_classification,
        )

        with zipfile.ZipFile(io.BytesIO(bundle)) as bundle_zip:
            file_names = set(bundle_zip.namelist())
            self.assertTrue(
                {
                    "alignment_loci.tsv",
                    "cross_locus_similarity.tsv",
                    "candidate_paralog_pairs.tsv",
                    "representative_bsr_loci.tsv",
                    "representative_bsr_inter_locus.tsv",
                    "representative_bsr_coverage.tsv",
                    "representative_bsr_pairs.tsv",
                    "representative_bsr_groups.tsv",
                }.issubset(file_names)
            )
            manifest = json.loads(bundle_zip.read("manifest.json"))
            self.assertTrue(
                all(status == "complete" for status in manifest["stage_status"].values())
            )
            analysis_provenance = manifest["operational_provenance"]["analyses"]
            self.assertEqual(
                analysis_provenance["alignment"]["method_identifier"],
                "mafft-fixture",
            )
            self.assertEqual(
                analysis_provenance["alignment"]["cache_reuse"]["status"],
                "not_applicable",
            )
            self.assertEqual(
                analysis_provenance["cross_locus_similarity"]["cache_reuse"],
                "reused",
            )
            self.assertEqual(
                analysis_provenance["representative_bsr"]["cache_reuse"],
                "not_reused",
            )
            summary = json.loads(bundle_zip.read("summary.json"))
            self.assertEqual(
                summary["alignment"]["alignment_method_identifier"],
                "mafft-fixture",
            )
            self.assertEqual(
                summary["similarity"]["method_identifier"],
                "similarity-fixture",
            )
            self.assertEqual(
                summary["representative_bsr"]["method_identifier"],
                "representative-fixture",
            )

    def test_bundle_filename_removes_known_archive_suffix(self) -> None:
        self.assertEqual(
            intrinsic_evaluation_bundle_filename(
                Path("Yersinia example.tar.zst")
            ),
            "Yersinia-example-intrinsic-evaluation.zip",
        )


if __name__ == "__main__":
    unittest.main()
