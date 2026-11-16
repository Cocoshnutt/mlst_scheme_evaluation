"""Tests for cached within-locus nucleotide and protein alignments."""

import json
import shutil
import subprocess
import tempfile
import unittest
from dataclasses import asdict
from io import StringIO
from pathlib import Path
from threading import Event
from unittest.mock import call, patch

from schema import alignment
from schema.alignment import (
    ALIGNMENT_METHOD_IDENTIFIER,
    ALIGNMENT_METHOD_IDENTIFIERS,
    ALIGNMENT_STATUS_ALIGNED,
    ALIGNMENT_STATUS_FAILED,
    ALIGNMENT_STATUS_INSUFFICIENT_SEQUENCES,
    ALIGNMENT_STRATEGY_AUTO,
    ALIGNMENT_STRATEGY_FFT_NS_1,
    ALIGNMENT_STRATEGY_FFT_NS_2,
    SCHEMA_ALIGNMENT_CACHE_FILENAME,
    SCHEMA_ALIGNMENT_CACHE_VERSION,
    AlignmentProgress,
    LocusAlignmentStatistics,
    SchemaAlignmentAnalysis,
    _AlignmentMetrics,
    _calculate_gap_metrics_from_aligned_sequences,
    _calculate_nucleotide_polymorphism_metrics,
    _get_mafft_strategy_arguments,
    _run_mafft_alignment,
    calculate_schema_alignment_summary,
    get_schema_alignment_analysis,
)


def _make_locus_statistics(
    locus_name: str,
    *,
    nucleotide_status: str = ALIGNMENT_STATUS_ALIGNED,
    nucleotide_gap_percentage: float | None = 0.0,
    amino_acid_status: str = ALIGNMENT_STATUS_ALIGNED,
    amino_acid_gap_percentage: float | None = 0.0,
) -> LocusAlignmentStatistics:
    """Create concise locus statistics for summary and cache tests."""
    return LocusAlignmentStatistics(
        locus_name=locus_name,
        total_allele_count=2,
        nucleotide_aligned_sequence_count=2,
        nucleotide_excluded_sequence_count=0,
        nucleotide_alignment_status=nucleotide_status,
        nucleotide_alignment_length=(
            6 if nucleotide_status == ALIGNMENT_STATUS_ALIGNED else None
        ),
        nucleotide_gap_character_count=0,
        nucleotide_gap_proportion_percentage=nucleotide_gap_percentage,
        nucleotide_callable_site_count=(
            6 if nucleotide_status == ALIGNMENT_STATUS_ALIGNED else None
        ),
        nucleotide_polymorphic_site_count=(
            0 if nucleotide_status == ALIGNMENT_STATUS_ALIGNED else None
        ),
        nucleotide_polymorphic_site_percentage=(
            0.0 if nucleotide_status == ALIGNMENT_STATUS_ALIGNED else None
        ),
        nucleotide_alignment_error=(
            "nucleotide failure"
            if nucleotide_status == ALIGNMENT_STATUS_FAILED
            else None
        ),
        amino_acid_aligned_sequence_count=2,
        amino_acid_excluded_sequence_count=0,
        amino_acid_alignment_status=amino_acid_status,
        amino_acid_alignment_length=(
            2 if amino_acid_status == ALIGNMENT_STATUS_ALIGNED else None
        ),
        amino_acid_gap_character_count=0,
        amino_acid_gap_proportion_percentage=amino_acid_gap_percentage,
        amino_acid_alignment_error=(
            "amino-acid failure"
            if amino_acid_status == ALIGNMENT_STATUS_FAILED
            else None
        ),
    )


class AlignmentMetricTest(unittest.TestCase):
    """Test alignment eligibility, gap calculation, and warning summaries."""

    def test_gap_proportion_counts_terminal_and_internal_gaps(self) -> None:
        """Use all gap characters and all alignment cells as the denominator."""
        alignment_metrics = _calculate_gap_metrics_from_aligned_sequences(
            ("A-C-", "ATCG")
        )

        self.assertEqual(alignment_metrics.alignment_length, 4)
        self.assertEqual(alignment_metrics.gap_character_count, 2)
        self.assertEqual(alignment_metrics.gap_proportion_percentage, 25.0)

    def test_polymorphism_uses_only_callable_canonical_observations(self) -> None:
        """Ignore gaps and ambiguity while retaining comparable substitutions."""
        (
            callable_site_count,
            polymorphic_site_count,
            polymorphic_site_percentage,
        ) = _calculate_nucleotide_polymorphism_metrics(
            ("ACGT", "atgt", "ANGT", "A-GT")
        )

        self.assertEqual(callable_site_count, 4)
        self.assertEqual(polymorphic_site_count, 1)
        self.assertEqual(polymorphic_site_percentage, 25.0)

    def test_indel_only_column_is_not_a_polymorphic_nucleotide_site(self) -> None:
        """Keep insertion/deletion variation in the separate gap statistic."""
        (
            callable_site_count,
            polymorphic_site_count,
            polymorphic_site_percentage,
        ) = _calculate_nucleotide_polymorphism_metrics(("AC-G", "ACAG"))

        self.assertEqual(callable_site_count, 3)
        self.assertEqual(polymorphic_site_count, 0)
        self.assertEqual(polymorphic_site_percentage, 0.0)

    def test_no_callable_columns_has_unavailable_percentage(self) -> None:
        """Avoid presenting zero percent when no site can be compared."""
        (
            callable_site_count,
            polymorphic_site_count,
            polymorphic_site_percentage,
        ) = _calculate_nucleotide_polymorphism_metrics(("NN", "--"))

        self.assertEqual(callable_site_count, 0)
        self.assertEqual(polymorphic_site_count, 0)
        self.assertIsNone(polymorphic_site_percentage)

    def test_strategy_arguments_match_documented_mafft_modes(self) -> None:
        """Use explicit fast strategies while retaining MAFFT Auto."""
        self.assertEqual(
            _get_mafft_strategy_arguments(ALIGNMENT_STRATEGY_FFT_NS_2),
            ("--retree", "2", "--maxiterate", "0"),
        )
        self.assertEqual(
            _get_mafft_strategy_arguments(ALIGNMENT_STRATEGY_FFT_NS_1),
            ("--retree", "1", "--maxiterate", "0"),
        )
        self.assertEqual(
            _get_mafft_strategy_arguments(ALIGNMENT_STRATEGY_AUTO),
            ("--auto",),
        )

    def test_selected_strategy_is_passed_to_mafft(self) -> None:
        """Construct the expected command for each selectable strategy."""
        successful_process = subprocess.CompletedProcess(
            args=["mafft"],
            returncode=0,
            stdout=">sequence_1\nATG-\n>sequence_2\nATGG\n",
            stderr="",
        )
        for alignment_strategy, expected_arguments in (
            (
                ALIGNMENT_STRATEGY_FFT_NS_2,
                ["--retree", "2", "--maxiterate", "0"],
            ),
            (
                ALIGNMENT_STRATEGY_FFT_NS_1,
                ["--retree", "1", "--maxiterate", "0"],
            ),
            (ALIGNMENT_STRATEGY_AUTO, ["--auto"]),
        ):
            with (
                self.subTest(alignment_strategy=alignment_strategy),
                patch(
                    "schema.alignment.subprocess.run",
                    return_value=successful_process,
                ) as subprocess_run,
            ):
                _run_mafft_alignment(
                    ("ATG", "ATGG"),
                    sequence_type="nucleotide",
                    mafft_executable="mafft",
                    alignment_strategy=alignment_strategy,
                )

            mafft_command = subprocess_run.call_args.args[0]
            self.assertEqual(mafft_command[:2], ["mafft", "--quiet"])
            self.assertEqual(
                mafft_command[2 : 2 + len(expected_arguments)], expected_arguments
            )
            self.assertIn("--thread", mafft_command)
            self.assertIn("--nuc", mafft_command)

    def test_unsupported_strategy_is_rejected(self) -> None:
        """Reject misspelled strategies before invoking MAFFT."""
        with (
            patch("schema.alignment.subprocess.run") as subprocess_run,
            self.assertRaisesRegex(ValueError, "alignment_strategy"),
        ):
            _run_mafft_alignment(
                ("ATG", "ATGG"),
                sequence_type="nucleotide",
                mafft_executable="mafft",
                alignment_strategy="unsupported",
            )
        subprocess_run.assert_not_called()

    def test_less_than_two_sequences_is_unavailable(self) -> None:
        """Do not report a misleading zero for a non-comparative locus."""
        with patch("schema.alignment.subprocess.run") as subprocess_run:
            alignment_metrics = _run_mafft_alignment(
                ("ATGTAA",),
                sequence_type="nucleotide",
                mafft_executable="mafft",
            )

        self.assertEqual(
            alignment_metrics.status,
            ALIGNMENT_STATUS_INSUFFICIENT_SEQUENCES,
        )
        self.assertIsNone(alignment_metrics.gap_proportion_percentage)
        subprocess_run.assert_not_called()

    def test_identical_sequences_skip_mafft_exactly(self) -> None:
        """Return exact zero-gap metrics without launching a subprocess."""
        for sequence_type, sequence in (
            ("nucleotide", "ATGAAATAA"),
            ("amino_acid", "MK"),
        ):
            with (
                self.subTest(sequence_type=sequence_type),
                patch("schema.alignment.subprocess.run") as subprocess_run,
            ):
                alignment_metrics = _run_mafft_alignment(
                    (sequence, sequence, sequence),
                    sequence_type=sequence_type,
                    mafft_executable="mafft",
                )

            self.assertEqual(alignment_metrics.status, ALIGNMENT_STATUS_ALIGNED)
            self.assertEqual(alignment_metrics.alignment_length, len(sequence))
            self.assertEqual(alignment_metrics.gap_character_count, 0)
            self.assertEqual(alignment_metrics.gap_proportion_percentage, 0.0)
            if sequence_type == "nucleotide":
                self.assertEqual(
                    alignment_metrics.nucleotide_callable_site_count,
                    len(sequence),
                )
                self.assertEqual(
                    alignment_metrics.nucleotide_polymorphic_site_count,
                    0,
                )
                self.assertEqual(
                    alignment_metrics.nucleotide_polymorphic_site_percentage,
                    0.0,
                )
            else:
                self.assertIsNone(
                    alignment_metrics.nucleotide_callable_site_count
                )
            subprocess_run.assert_not_called()

    def test_empty_locus_is_unavailable_for_both_alignments(self) -> None:
        """Report an empty FASTA locus without invoking MAFFT."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            allele_file_path = schema_directory_path / "empty.fasta"
            allele_file_path.write_text("", encoding="utf-8")
            with patch("schema.alignment.subprocess.run") as subprocess_run:
                locus_statistics = alignment._analyze_locus(
                    schema_directory_path,
                    allele_file_path,
                    "mafft",
                )

        self.assertEqual(
            locus_statistics.nucleotide_alignment_status,
            ALIGNMENT_STATUS_INSUFFICIENT_SEQUENCES,
        )
        self.assertEqual(
            locus_statistics.amino_acid_alignment_status,
            ALIGNMENT_STATUS_INSUFFICIENT_SEQUENCES,
        )
        subprocess_run.assert_not_called()

    def test_locus_eligibility_and_duplicate_record_inclusion(self) -> None:
        """Apply the approved nucleotide and valid-CDS inclusion rules."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            allele_file_path = schema_directory_path / "locus.fasta"
            allele_file_path.write_text(
                (
                    ">locus_1\nATGAAATAA\n"
                    ">locus_2\nATGAAATAA\n"
                    ">locus_3\nATGAAAAAATAA\n"
                    ">locus_4\nATGNNNTAA\n"
                    ">locus_5\nATG?TAA\n"
                    ">locus_6\n"
                ),
                encoding="utf-8",
            )
            aligned_metrics = _AlignmentMetrics(
                status=ALIGNMENT_STATUS_ALIGNED,
                alignment_length=9,
                gap_character_count=0,
                gap_proportion_percentage=0.0,
                error=None,
            )

            with patch(
                "schema.alignment._run_mafft_alignment",
                return_value=aligned_metrics,
            ) as run_mafft_alignment:
                locus_statistics = alignment._analyze_locus(
                    schema_directory_path,
                    allele_file_path,
                    "mafft",
                )

        self.assertEqual(locus_statistics.total_allele_count, 6)
        self.assertEqual(locus_statistics.nucleotide_aligned_sequence_count, 4)
        self.assertEqual(locus_statistics.nucleotide_excluded_sequence_count, 2)
        self.assertEqual(locus_statistics.amino_acid_aligned_sequence_count, 3)
        self.assertEqual(locus_statistics.amino_acid_excluded_sequence_count, 3)
        self.assertEqual(
            run_mafft_alignment.call_args_list,
            [
                call(
                    (
                        "ATGAAATAA",
                        "ATGAAATAA",
                        "ATGAAAAAATAA",
                        "ATGNNNTAA",
                    ),
                    sequence_type="nucleotide",
                    mafft_executable="mafft",
                    alignment_strategy=ALIGNMENT_STRATEGY_FFT_NS_2,
                ),
                call(
                    ("MK", "MK", "MKK"),
                    sequence_type="amino_acid",
                    mafft_executable="mafft",
                    alignment_strategy=ALIGNMENT_STRATEGY_FFT_NS_2,
                ),
            ],
        )

    def test_mafft_failure_is_reported_per_alignment(self) -> None:
        """Keep a locus-level tool failure visible without fabricating a value."""
        failed_process = subprocess.CompletedProcess(
            args=["mafft"],
            returncode=1,
            stdout="",
            stderr="alignment failed",
        )
        with patch("schema.alignment.subprocess.run", return_value=failed_process):
            alignment_metrics = _run_mafft_alignment(
                ("ATGTAA", "ATGAAATAA"),
                sequence_type="nucleotide",
                mafft_executable="mafft",
            )

        self.assertEqual(alignment_metrics.status, ALIGNMENT_STATUS_FAILED)
        self.assertIsNone(alignment_metrics.gap_proportion_percentage)
        self.assertEqual(alignment_metrics.error, "alignment failed")

    def test_warning_threshold_is_strict_and_uses_aligned_loci(self) -> None:
        """Flag only values greater than the advisory threshold."""
        alignment_analysis = SchemaAlignmentAnalysis(
            mafft_version="test",
            alignment_method_identifier=ALIGNMENT_METHOD_IDENTIFIER,
            locus_statistics=(
                _make_locus_statistics("boundary", amino_acid_gap_percentage=5.0),
                _make_locus_statistics("above", amino_acid_gap_percentage=5.1),
                _make_locus_statistics(
                    "unavailable",
                    nucleotide_status=ALIGNMENT_STATUS_INSUFFICIENT_SEQUENCES,
                    nucleotide_gap_percentage=None,
                    amino_acid_status=ALIGNMENT_STATUS_INSUFFICIENT_SEQUENCES,
                    amino_acid_gap_percentage=None,
                ),
                _make_locus_statistics(
                    "failed",
                    nucleotide_status=ALIGNMENT_STATUS_FAILED,
                    nucleotide_gap_percentage=None,
                    amino_acid_status=ALIGNMENT_STATUS_FAILED,
                    amino_acid_gap_percentage=None,
                ),
            ),
        )

        alignment_summary = calculate_schema_alignment_summary(
            alignment_analysis,
            amino_acid_gap_warning_threshold_percentage=5.0,
        )

        self.assertEqual(alignment_summary.amino_acid_aligned_locus_count, 2)
        self.assertEqual(alignment_summary.amino_acid_unavailable_locus_count, 1)
        self.assertEqual(alignment_summary.amino_acid_failed_locus_count, 1)
        self.assertEqual(alignment_summary.amino_acid_gap_warning_locus_count, 1)
        self.assertEqual(alignment_summary.amino_acid_gap_warning_loci, ("above",))
        self.assertEqual(
            alignment_summary.amino_acid_gap_warning_locus_percentage,
            50.0,
        )


class AlignmentCacheTest(unittest.TestCase):
    """Test persistent alignment-summary cache behavior."""

    def test_strategies_have_independent_reusable_cache_entries(self) -> None:
        """Keep completed work when switching strategies and switching back."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            (schema_directory_path / "locus.fasta").write_text(
                ">locus_1\nATGTAA\n>locus_2\nATGAAATAA\n",
                encoding="utf-8",
            )

            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    side_effect=lambda schema_path, allele_path, executable, strategy, cached: (
                        _make_locus_statistics(
                            allele_path.stem,
                            nucleotide_gap_percentage=(
                                2.0 if strategy == ALIGNMENT_STRATEGY_FFT_NS_2 else 1.0
                            ),
                        )
                    ),
                ) as analyze_locus,
            ):
                fft_ns_2_analysis = get_schema_alignment_analysis(
                    schema_directory_path,
                    alignment_strategy=ALIGNMENT_STRATEGY_FFT_NS_2,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )
                fft_ns_1_analysis = get_schema_alignment_analysis(
                    schema_directory_path,
                    alignment_strategy=ALIGNMENT_STRATEGY_FFT_NS_1,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    side_effect=AssertionError("MAFFT version was requested"),
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    side_effect=AssertionError("locus was realigned"),
                ),
            ):
                reused_fft_ns_2_analysis = get_schema_alignment_analysis(
                    schema_directory_path,
                    alignment_strategy=ALIGNMENT_STRATEGY_FFT_NS_2,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

            cached_data = json.loads(
                (schema_directory_path / SCHEMA_ALIGNMENT_CACHE_FILENAME).read_text(
                    encoding="utf-8"
                )
            )

        self.assertEqual(analyze_locus.call_count, 2)
        self.assertEqual(
            fft_ns_2_analysis.alignment_strategy, ALIGNMENT_STRATEGY_FFT_NS_2
        )
        self.assertEqual(
            fft_ns_1_analysis.alignment_strategy, ALIGNMENT_STRATEGY_FFT_NS_1
        )
        self.assertEqual(reused_fft_ns_2_analysis, fft_ns_2_analysis)
        self.assertEqual(
            set(cached_data["methods"]),
            {
                ALIGNMENT_METHOD_IDENTIFIERS[ALIGNMENT_STRATEGY_FFT_NS_1],
                ALIGNMENT_METHOD_IDENTIFIERS[ALIGNMENT_STRATEGY_FFT_NS_2],
            },
        )

    def test_cache_reuse_does_not_require_mafft(self) -> None:
        """Load derived metrics without rerunning or locating the aligner."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            (schema_directory_path / "locus.fasta").write_text(
                ">locus_1\nATGTAA\n>locus_2\nATGAAATAA\n",
                encoding="utf-8",
            )
            expected_locus_statistics = _make_locus_statistics("locus")
            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    return_value=expected_locus_statistics,
                ),
            ):
                initial_analysis = get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    side_effect=AssertionError("MAFFT version was requested"),
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    side_effect=AssertionError("locus was realigned"),
                ),
            ):
                cached_analysis = get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

        self.assertEqual(cached_analysis, initial_analysis)

    def test_invalid_worker_count_is_rejected_on_cache_miss(self) -> None:
        """Reject nonsensical concurrency without invoking MAFFT."""
        with (
            tempfile.TemporaryDirectory() as temporary_directory,
            self.assertRaisesRegex(ValueError, "at least one"),
        ):
            get_schema_alignment_analysis(
                temporary_directory,
                maximum_parallel_loci=0,
                progress_output_stream=None,
            )

    def test_source_change_invalidates_cached_alignment_metrics(self) -> None:
        """Realign when a locus FASTA path, size, or modification time changes."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            allele_file_path = schema_directory_path / "locus.fasta"
            allele_file_path.write_text(
                ">locus_1\nATGTAA\n>locus_2\nATGAAATAA\n",
                encoding="utf-8",
            )
            initial_locus_statistics = _make_locus_statistics("locus")
            changed_locus_statistics = _make_locus_statistics(
                "locus",
                nucleotide_gap_percentage=10.0,
            )
            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    return_value=initial_locus_statistics,
                ),
            ):
                get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

            allele_file_path.write_text(
                ">locus_1\nATGTAA\n>locus_2\nATGAAATAA\n>locus_3\nATGAAA\n",
                encoding="utf-8",
            )
            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    return_value=changed_locus_statistics,
                ) as analyze_locus,
            ):
                changed_analysis = get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

        analyze_locus.assert_called_once()
        self.assertEqual(
            changed_analysis.locus_statistics[0].nucleotide_gap_proportion_percentage,
            10.0,
        )

    def test_only_changed_and_added_loci_are_analyzed(self) -> None:
        """Reuse unchanged entries while processing changed and new loci."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            alpha_path = schema_directory_path / "alpha.fasta"
            beta_path = schema_directory_path / "beta.fasta"
            alpha_path.write_text(">a1\nATGTAA\n>a2\nATGAAATAA\n", encoding="utf-8")
            beta_path.write_text(">b1\nATGTAA\n>b2\nATGAAATAA\n", encoding="utf-8")
            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    side_effect=lambda schema_path, allele_path, executable, strategy, cached: (
                        _make_locus_statistics(allele_path.stem)
                    ),
                ),
            ):
                get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

            beta_path.write_text(
                ">b1\nATGTAA\n>b2\nATGAAATAA\n>b3\nATGAAA\n",
                encoding="utf-8",
            )
            gamma_path = schema_directory_path / "gamma.fasta"
            gamma_path.write_text(
                ">g1\nATGTAA\n>g2\nATGAAATAA\n",
                encoding="utf-8",
            )
            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    side_effect=lambda schema_path, allele_path, executable, strategy, cached: (
                        _make_locus_statistics(allele_path.stem)
                    ),
                ) as analyze_locus,
            ):
                alignment_analysis = get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

        self.assertEqual(
            {
                call_arguments.args[1].name
                for call_arguments in analyze_locus.call_args_list
            },
            {"beta.fasta", "gamma.fasta"},
        )
        self.assertEqual(
            tuple(
                locus_statistics.locus_name
                for locus_statistics in alignment_analysis.locus_statistics
            ),
            ("alpha", "beta", "gamma"),
        )

    def test_deleted_locus_is_removed_without_realigning(self) -> None:
        """Clean obsolete entries while reusing all remaining current loci."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            for locus_name in ("alpha", "beta"):
                (schema_directory_path / f"{locus_name}.fasta").write_text(
                    ">allele_1\nATGTAA\n>allele_2\nATGAAATAA\n",
                    encoding="utf-8",
                )
            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    side_effect=lambda schema_path, allele_path, executable, strategy, cached: (
                        _make_locus_statistics(allele_path.stem)
                    ),
                ),
            ):
                get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

            (schema_directory_path / "beta.fasta").unlink()
            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    side_effect=AssertionError("MAFFT version was requested"),
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    side_effect=AssertionError("a locus was realigned"),
                ),
            ):
                alignment_analysis = get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

            cached_data = json.loads(
                (schema_directory_path / SCHEMA_ALIGNMENT_CACHE_FILENAME).read_text(
                    encoding="utf-8"
                )
            )

        self.assertEqual(
            tuple(
                locus_statistics.locus_name
                for locus_statistics in alignment_analysis.locus_statistics
            ),
            ("alpha",),
        )
        self.assertEqual(
            tuple(cached_data["methods"][ALIGNMENT_METHOD_IDENTIFIER]["loci"]),
            ("alpha.fasta",),
        )

    def test_failed_locus_is_retried_on_next_call(self) -> None:
        """Do not treat a tool failure as a reusable biological result."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            (schema_directory_path / "locus.fasta").write_text(
                ">allele_1\nATGTAA\n>allele_2\nATGAAATAA\n",
                encoding="utf-8",
            )
            failed_statistics = _make_locus_statistics(
                "locus",
                nucleotide_status=ALIGNMENT_STATUS_FAILED,
                nucleotide_gap_percentage=None,
            )
            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    return_value=failed_statistics,
                ),
            ):
                get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    return_value=_make_locus_statistics("locus"),
                ) as analyze_locus,
            ):
                retried_analysis = get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

        analyze_locus.assert_called_once()
        self.assertEqual(
            retried_analysis.locus_statistics[0].nucleotide_alignment_status,
            ALIGNMENT_STATUS_ALIGNED,
        )

    def test_interrupted_run_checkpoints_completed_locus(self) -> None:
        """Resume after an exception without repeating completed locus work."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            for locus_name in ("alpha", "beta"):
                (schema_directory_path / f"{locus_name}.fasta").write_text(
                    ">allele_1\nATGTAA\n>allele_2\nATGAAATAA\n",
                    encoding="utf-8",
                )
            alpha_progress_reported = Event()

            def analyze_then_fail(
                schema_path, allele_path, executable, strategy, cached
            ):
                if allele_path.stem == "alpha":
                    return _make_locus_statistics("alpha")
                alpha_progress_reported.wait(timeout=2)
                raise RuntimeError("interrupted test run")

            def record_progress(progress: AlignmentProgress) -> None:
                if progress.current_locus_name == "alpha":
                    alpha_progress_reported.set()

            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    side_effect=analyze_then_fail,
                ),
                self.assertRaisesRegex(RuntimeError, "interrupted test run"),
            ):
                get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                    progress_callback=record_progress,
                )

            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    return_value=_make_locus_statistics("beta"),
                ) as analyze_locus,
            ):
                resumed_analysis = get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

        analyze_locus.assert_called_once()
        self.assertEqual(analyze_locus.call_args.args[1].name, "beta.fasta")
        self.assertEqual(
            tuple(
                locus_statistics.locus_name
                for locus_statistics in resumed_analysis.locus_statistics
            ),
            ("alpha", "beta"),
        )

    def test_version_one_cache_is_migrated_and_upgraded(self) -> None:
        """Pass legacy results into the one-time polymorphism upgrade."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            allele_file_path = schema_directory_path / "locus.fasta"
            allele_file_path.write_text(
                ">allele_1\nATGTAA\n>allele_2\nATGAAATAA\n",
                encoding="utf-8",
            )
            allele_file_status = allele_file_path.stat()
            legacy_locus_statistics = asdict(_make_locus_statistics("locus"))
            for field_name in (
                "nucleotide_callable_site_count",
                "nucleotide_polymorphic_site_count",
                "nucleotide_polymorphic_site_percentage",
            ):
                legacy_locus_statistics.pop(field_name)
            version_one_cache = {
                "cache_version": 1,
                "alignment_method_identifier": ALIGNMENT_METHOD_IDENTIFIERS[
                    ALIGNMENT_STRATEGY_AUTO
                ],
                "mafft_version": "MAFFT test version",
                "source_signature": [
                    [
                        "locus.fasta",
                        allele_file_status.st_size,
                        allele_file_status.st_mtime_ns,
                    ]
                ],
                "locus_statistics": [legacy_locus_statistics],
            }
            (schema_directory_path / SCHEMA_ALIGNMENT_CACHE_FILENAME).write_text(
                json.dumps(version_one_cache),
                encoding="utf-8",
            )

            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    return_value=_make_locus_statistics("locus"),
                ) as analyze_locus,
            ):
                alignment_analysis = get_schema_alignment_analysis(
                    schema_directory_path,
                    alignment_strategy=ALIGNMENT_STRATEGY_AUTO,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

            migrated_cache = json.loads(
                (schema_directory_path / SCHEMA_ALIGNMENT_CACHE_FILENAME).read_text(
                    encoding="utf-8"
                )
            )

        self.assertEqual(len(alignment_analysis.locus_statistics), 1)
        analyze_locus.assert_called_once()
        self.assertIsNotNone(analyze_locus.call_args.args[4])
        self.assertEqual(
            migrated_cache["cache_version"],
            SCHEMA_ALIGNMENT_CACHE_VERSION,
        )
        self.assertIn(
            "locus.fasta",
            migrated_cache["methods"][
                ALIGNMENT_METHOD_IDENTIFIERS[ALIGNMENT_STRATEGY_AUTO]
            ]["loci"],
        )

    def test_version_two_cache_is_migrated_and_upgraded(self) -> None:
        """Preserve legacy Auto metrics during the polymorphism upgrade."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            allele_file_path = schema_directory_path / "locus.fasta"
            allele_file_path.write_text(
                ">allele_1\nATGTAA\n>allele_2\nATGAAATAA\n",
                encoding="utf-8",
            )
            allele_file_status = allele_file_path.stat()
            legacy_locus_statistics = asdict(_make_locus_statistics("locus"))
            for field_name in (
                "nucleotide_callable_site_count",
                "nucleotide_polymorphic_site_count",
                "nucleotide_polymorphic_site_percentage",
            ):
                legacy_locus_statistics.pop(field_name)
            version_two_cache = {
                "cache_version": 2,
                "alignment_method_identifier": ALIGNMENT_METHOD_IDENTIFIERS[
                    ALIGNMENT_STRATEGY_AUTO
                ],
                "mafft_version": "MAFFT test version",
                "loci": {
                    "locus.fasta": {
                        "file_size": allele_file_status.st_size,
                        "modification_time_ns": allele_file_status.st_mtime_ns,
                        "statistics": legacy_locus_statistics,
                    }
                },
            }
            (schema_directory_path / SCHEMA_ALIGNMENT_CACHE_FILENAME).write_text(
                json.dumps(version_two_cache),
                encoding="utf-8",
            )

            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    return_value=_make_locus_statistics("locus"),
                ) as analyze_locus,
            ):
                alignment_analysis = get_schema_alignment_analysis(
                    schema_directory_path,
                    alignment_strategy=ALIGNMENT_STRATEGY_AUTO,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

            migrated_cache = json.loads(
                (schema_directory_path / SCHEMA_ALIGNMENT_CACHE_FILENAME).read_text(
                    encoding="utf-8"
                )
            )

        self.assertEqual(alignment_analysis.alignment_strategy, ALIGNMENT_STRATEGY_AUTO)
        analyze_locus.assert_called_once()
        self.assertIsNotNone(analyze_locus.call_args.args[4])
        self.assertEqual(
            migrated_cache["cache_version"],
            SCHEMA_ALIGNMENT_CACHE_VERSION,
        )
        self.assertIn(
            ALIGNMENT_METHOD_IDENTIFIERS[ALIGNMENT_STRATEGY_AUTO],
            migrated_cache["methods"],
        )

    def test_version_three_cache_reuses_protein_while_upgrading_nucleotide(self) -> None:
        """Run only nucleotide MAFFT when adding polymorphism to a v3 locus."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            allele_file_path = schema_directory_path / "locus.fasta"
            allele_file_path.write_text(
                ">allele_1\nATGTAA\n>allele_2\nATGAAATAA\n",
                encoding="utf-8",
            )
            allele_file_status = allele_file_path.stat()
            legacy_locus_statistics = asdict(_make_locus_statistics("locus"))
            for field_name in (
                "nucleotide_callable_site_count",
                "nucleotide_polymorphic_site_count",
                "nucleotide_polymorphic_site_percentage",
            ):
                legacy_locus_statistics.pop(field_name)
            version_three_cache = {
                "cache_version": 3,
                "methods": {
                    ALIGNMENT_METHOD_IDENTIFIER: {
                        "mafft_version": "MAFFT test version",
                        "loci": {
                            "locus.fasta": {
                                "file_size": allele_file_status.st_size,
                                "modification_time_ns": (
                                    allele_file_status.st_mtime_ns
                                ),
                                "statistics": legacy_locus_statistics,
                            }
                        },
                    }
                },
            }
            (schema_directory_path / SCHEMA_ALIGNMENT_CACHE_FILENAME).write_text(
                json.dumps(version_three_cache),
                encoding="utf-8",
            )
            nucleotide_metrics = _AlignmentMetrics(
                status=ALIGNMENT_STATUS_ALIGNED,
                alignment_length=9,
                gap_character_count=3,
                gap_proportion_percentage=100 / 6,
                error=None,
                nucleotide_callable_site_count=9,
                nucleotide_polymorphic_site_count=1,
                nucleotide_polymorphic_site_percentage=100 / 9,
            )

            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._run_mafft_alignment",
                    return_value=nucleotide_metrics,
                ) as run_mafft_alignment,
            ):
                alignment_analysis = get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

            upgraded_cache = json.loads(
                (schema_directory_path / SCHEMA_ALIGNMENT_CACHE_FILENAME).read_text(
                    encoding="utf-8"
                )
            )

        run_mafft_alignment.assert_called_once()
        self.assertEqual(
            run_mafft_alignment.call_args.kwargs["sequence_type"],
            "nucleotide",
        )
        locus_statistics = alignment_analysis.locus_statistics[0]
        self.assertEqual(locus_statistics.nucleotide_polymorphic_site_count, 1)
        self.assertEqual(locus_statistics.amino_acid_alignment_length, 2)
        upgraded_locus_cache = upgraded_cache["methods"][
            ALIGNMENT_METHOD_IDENTIFIER
        ]["loci"]["locus.fasta"]
        self.assertTrue(
            upgraded_locus_cache[
                "nucleotide_polymorphism_metrics_calculated"
            ]
        )

    def test_empty_schema_returns_empty_analysis_and_summary(self) -> None:
        """Handle a schema with no recognized locus FASTA files."""
        with (
            tempfile.TemporaryDirectory() as temporary_directory,
            patch(
                "schema.alignment._get_mafft_version",
                return_value="MAFFT test version",
            ),
            patch("schema.alignment._analyze_locus") as analyze_locus,
        ):
            alignment_analysis = get_schema_alignment_analysis(
                temporary_directory,
                maximum_parallel_loci=1,
                progress_output_stream=None,
            )
            alignment_summary = calculate_schema_alignment_summary(alignment_analysis)

        analyze_locus.assert_not_called()
        self.assertEqual(alignment_analysis.locus_statistics, ())
        self.assertEqual(alignment_summary.locus_count, 0)
        self.assertEqual(
            alignment_summary.amino_acid_gap_warning_locus_percentage,
            0.0,
        )

    def test_default_worker_count_uses_every_process_available_cpu(self) -> None:
        """Respect CPU affinity and cgroup limits when resolving concurrency."""
        with patch("schema.alignment.os.process_cpu_count", return_value=12):
            worker_count = alignment._resolve_maximum_parallel_loci(None)

        self.assertEqual(worker_count, 12)

    def test_progress_reports_completed_loci_and_worker_count(self) -> None:
        """Expose cache-miss progress in the caller-provided terminal stream."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            for locus_name in ("alpha", "beta"):
                (schema_directory_path / f"{locus_name}.fasta").write_text(
                    ">allele_1\nATGTAA\n>allele_2\nATGAAATAA\n",
                    encoding="utf-8",
                )
            progress_output_stream = StringIO()

            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    side_effect=lambda schema_path, allele_path, executable, strategy, cached: (
                        _make_locus_statistics(allele_path.stem)
                    ),
                ),
            ):
                alignment_analysis = get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=2,
                    progress_output_stream=progress_output_stream,
                )

        progress_output = progress_output_stream.getvalue()
        self.assertIn("0/2 loci", progress_output)
        self.assertIn("2/2 loci", progress_output)
        self.assertIn("FFT-NS-2", progress_output)
        self.assertIn("2 workers", progress_output)
        self.assertIn("loci/s", progress_output)
        self.assertIn("complete", progress_output)
        self.assertEqual(
            tuple(
                locus_statistics.locus_name
                for locus_statistics in alignment_analysis.locus_statistics
            ),
            ("alpha", "beta"),
        )

    def test_progress_callback_reports_cache_reuse_and_new_work(self) -> None:
        """Provide notebook progress for partial and complete cache reuse."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            alpha_path = schema_directory_path / "alpha.fasta"
            alpha_path.write_text(
                ">allele_1\nATGTAA\n>allele_2\nATGAAATAA\n",
                encoding="utf-8",
            )
            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    return_value=_make_locus_statistics("alpha"),
                ),
            ):
                get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )

            beta_path = schema_directory_path / "beta.fasta"
            beta_path.write_text(
                ">allele_1\nATGTAA\n>allele_2\nATGAAATAA\n",
                encoding="utf-8",
            )
            partial_progress: list[AlignmentProgress] = []
            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    return_value=_make_locus_statistics("beta"),
                ),
            ):
                get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                    progress_callback=partial_progress.append,
                )

            cached_progress: list[AlignmentProgress] = []
            with patch(
                "schema.alignment._get_mafft_version",
                side_effect=AssertionError("MAFFT version was requested"),
            ):
                get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                    progress_callback=cached_progress.append,
                )

        self.assertEqual(partial_progress[0].completed_locus_count, 1)
        self.assertEqual(partial_progress[0].cached_locus_count, 1)
        self.assertIsNone(partial_progress[0].current_locus_name)
        self.assertEqual(partial_progress[-1].completed_locus_count, 2)
        self.assertEqual(partial_progress[-1].current_locus_name, "beta")
        self.assertEqual(len(cached_progress), 1)
        self.assertEqual(cached_progress[0].completed_locus_count, 2)
        self.assertEqual(cached_progress[0].cached_locus_count, 2)

    def test_scheduler_handles_five_hundred_loci(self) -> None:
        """Schedule a cgMLST-sized batch while preserving deterministic order."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            for locus_index in range(500):
                (schema_directory_path / f"locus_{locus_index:03d}.fasta").write_text(
                    ">allele_1\nATGTAA\n>allele_2\nATGAAATAA\n",
                    encoding="utf-8",
                )

            with (
                patch(
                    "schema.alignment._get_mafft_version",
                    return_value="MAFFT test version",
                ),
                patch(
                    "schema.alignment._analyze_locus",
                    side_effect=lambda schema_path, allele_path, executable, strategy, cached: (
                        _make_locus_statistics(allele_path.stem)
                    ),
                ),
            ):
                alignment_analysis = get_schema_alignment_analysis(
                    schema_directory_path,
                    maximum_parallel_loci=8,
                    progress_output_stream=None,
                )

        self.assertEqual(len(alignment_analysis.locus_statistics), 500)
        self.assertEqual(
            alignment_analysis.locus_statistics[0].locus_name,
            "locus_000",
        )
        self.assertEqual(
            alignment_analysis.locus_statistics[-1].locus_name,
            "locus_499",
        )


@unittest.skipUnless(shutil.which("mafft"), "MAFFT is not available")
class MafftIntegrationTest(unittest.TestCase):
    """Run one tiny real alignment when the development dependency exists."""

    def test_tiny_schema_alignment(self) -> None:
        """Verify all strategies with real nucleotide and protein alignments."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            (schema_directory_path / "locus.fasta").write_text(
                (
                    ">locus_1\nATGAAATAA\n"
                    ">locus_2\nATGAAAAAATAA\n"
                    ">locus_3\nATGGAATAA\n"
                ),
                encoding="utf-8",
            )

            alignment_analyses = tuple(
                get_schema_alignment_analysis(
                    schema_directory_path,
                    alignment_strategy=alignment_strategy,
                    maximum_parallel_loci=1,
                    progress_output_stream=None,
                )
                for alignment_strategy in (
                    ALIGNMENT_STRATEGY_FFT_NS_2,
                    ALIGNMENT_STRATEGY_FFT_NS_1,
                    ALIGNMENT_STRATEGY_AUTO,
                )
            )

        for alignment_analysis in alignment_analyses:
            with self.subTest(
                alignment_strategy=alignment_analysis.alignment_strategy
            ):
                locus_statistics = alignment_analysis.locus_statistics[0]
                self.assertEqual(
                    locus_statistics.nucleotide_alignment_status,
                    ALIGNMENT_STATUS_ALIGNED,
                )
                self.assertEqual(
                    locus_statistics.amino_acid_alignment_status,
                    ALIGNMENT_STATUS_ALIGNED,
                )
                self.assertGreater(
                    locus_statistics.nucleotide_gap_proportion_percentage or 0,
                    0,
                )
                self.assertGreater(
                    locus_statistics.amino_acid_gap_proportion_percentage or 0,
                    0,
                )
                self.assertGreater(
                    locus_statistics.nucleotide_callable_site_count or 0,
                    0,
                )
                self.assertGreater(
                    locus_statistics.nucleotide_polymorphic_site_count or 0,
                    0,
                )
                self.assertGreater(
                    locus_statistics.nucleotide_polymorphic_site_percentage or 0,
                    0,
                )


if __name__ == "__main__":
    unittest.main()
