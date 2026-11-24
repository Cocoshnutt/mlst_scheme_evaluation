"""Tests for shared BLAST operations and analysis command contracts."""

import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import call, patch

from scheme import bsr, similarity
from scheme.blast import (
    create_blast_database,
    get_blast_version,
    run_blast_command,
)
from scheme.bsr import (
    _CollectedLocus,
    _parse_score_fields,
    _ProteinAllele,
    _run_blastp_analysis,
)
from scheme.similarity import (
    SEQUENCE_LEVEL_AMINO_ACID,
    SEQUENCE_LEVEL_NUCLEOTIDE,
    _parse_blast_output,
    _run_similarity_level,
    _SelectedAlleleSequence,
)


class BlastProcessTest(unittest.TestCase):
    """Protect shared subprocess and version-discovery behavior."""

    def test_checked_command_uses_the_supplied_vector(self) -> None:
        """Run commands without a shell and preserve every argument."""
        completed_process = subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout="",
            stderr="",
        )
        with patch(
            "scheme.blast.subprocess.run",
            return_value=completed_process,
        ) as subprocess_run:
            run_blast_command(["blastp", "-query", "input.fasta"], "BLASTp")

        subprocess_run.assert_called_once_with(
            ["blastp", "-query", "input.fasta"],
            check=False,
            capture_output=True,
            text=True,
        )

    def test_missing_executable_preserves_operational_error(self) -> None:
        """Name both the failed operation and missing executable."""
        with (
            patch(
                "scheme.blast.subprocess.run",
                side_effect=FileNotFoundError,
            ),
            self.assertRaisesRegex(
                RuntimeError,
                "makeblastdb .* executable was not found: missing-db-tool",
            ),
        ):
            run_blast_command(["missing-db-tool"], "makeblastdb (protein)")

    def test_nonzero_command_prefers_stderr_and_reports_exit_code(self) -> None:
        """Retain the established concise BLAST failure message."""
        completed_process = subprocess.CompletedProcess(
            args=[],
            returncode=2,
            stdout="stdout detail",
            stderr="stderr detail",
        )
        with (
            patch(
                "scheme.blast.subprocess.run",
                return_value=completed_process,
            ),
            self.assertRaisesRegex(
                RuntimeError,
                "BLASTn failed with exit code 2: stderr detail",
            ),
        ):
            run_blast_command(["blastn"], "BLASTn")

    def test_version_uses_first_output_line(self) -> None:
        """Keep cache provenance independent of extra version text."""
        completed_process = subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout="blastp: 2.17.0+\nPackage details\n",
            stderr="",
        )
        with patch(
            "scheme.blast.subprocess.run",
            return_value=completed_process,
        ) as subprocess_run:
            version = get_blast_version("custom-blastp")

        self.assertEqual(version, "blastp: 2.17.0+")
        subprocess_run.assert_called_once_with(
            ["custom-blastp", "-version"],
            check=False,
            capture_output=True,
            text=True,
        )

    def test_version_failures_preserve_the_executable_name(self) -> None:
        """Differentiate absent tools from unsuccessful version commands."""
        with (
            patch(
                "scheme.blast.subprocess.run",
                side_effect=FileNotFoundError,
            ),
            self.assertRaisesRegex(
                RuntimeError,
                "BLAST executable was not found: unavailable-blastp",
            ),
        ):
            get_blast_version("unavailable-blastp")

        failed_process = subprocess.CompletedProcess(
            args=[],
            returncode=1,
            stdout="version stdout",
            stderr="version stderr",
        )
        with (
            patch(
                "scheme.blast.subprocess.run",
                return_value=failed_process,
            ),
            self.assertRaisesRegex(
                RuntimeError,
                "Could not determine the version of broken-blastn: version stderr",
            ),
        ):
            get_blast_version("broken-blastn")

    def test_database_creation_preserves_makeblastdb_arguments(self) -> None:
        """Use parsed sequence IDs for both supported database types."""
        input_file_path = Path("selected alleles.fasta")
        database_path = Path("temporary database")
        with patch("scheme.blast.run_blast_command") as run_command:
            create_blast_database(
                input_fasta_file_path=input_file_path,
                database_path=database_path,
                database_type="nucl",
                makeblastdb_executable="custom-makeblastdb",
                operation_name="makeblastdb (nucleotide)",
            )

        run_command.assert_called_once_with(
            [
                "custom-makeblastdb",
                "-in",
                str(input_file_path),
                "-dbtype",
                "nucl",
                "-out",
                str(database_path),
                "-parse_seqids",
            ],
            "makeblastdb (nucleotide)",
        )


class AnalysisBlastCommandContractTest(unittest.TestCase):
    """Protect analysis-owned BLAST search formats and arguments."""

    def test_analysis_parsers_reject_the_wrong_field_count(self) -> None:
        """Keep the twelve- and three-column output contracts independent."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_file_path = Path(temporary_directory) / "malformed.tsv"
            output_file_path.write_text("query\tsubject\t10\n", encoding="utf-8")
            with self.assertRaisesRegex(
                RuntimeError,
                "expected 12 fields, found 3",
            ):
                _parse_blast_output(
                    output_file_path,
                    sequence_level=SEQUENCE_LEVEL_NUCLEOTIDE,
                )

        with self.assertRaisesRegex(
            RuntimeError,
            "expected 3 fields, found 2",
        ):
            _parse_score_fields("query\tsubject\n", 7)

    def test_similarity_search_commands_remain_sequence_level_specific(self) -> None:
        """Keep blastn/blastp tasks, output fields, and thread counts exact."""
        selected_allele = _SelectedAlleleSequence(
            internal_identifier="sequence_000000001",
            locus_name="gene_a",
            original_identifier="gene_a_1",
            sequence="ATGAAATAA",
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            working_directory_path = Path(temporary_directory)
            (working_directory_path / "nucleotide.tsv").write_text(
                "",
                encoding="utf-8",
            )
            (working_directory_path / "amino_acid.tsv").write_text(
                "",
                encoding="utf-8",
            )
            with (
                patch.object(similarity, "create_blast_database") as create_database,
                patch.object(similarity, "_run_command") as run_command,
            ):
                _run_similarity_level(
                    sequence_level=SEQUENCE_LEVEL_NUCLEOTIDE,
                    selected_alleles=(selected_allele,),
                    working_directory_path=working_directory_path,
                    blast_executable="custom-blastn",
                    makeblastdb_executable="custom-makeblastdb",
                    blast_evalue=0.001,
                    thread_count=3,
                )
                _run_similarity_level(
                    sequence_level=SEQUENCE_LEVEL_AMINO_ACID,
                    selected_alleles=(selected_allele,),
                    working_directory_path=working_directory_path,
                    blast_executable="custom-blastp",
                    makeblastdb_executable="custom-makeblastdb",
                    blast_evalue=0.001,
                    thread_count=4,
                )

        self.assertEqual(
            [
                database_call.kwargs["database_type"]
                for database_call in create_database.call_args_list
            ],
            ["nucl", "prot"],
        )
        common_output_format = (
            "6 qseqid sseqid qlen slen qstart qend sstart send length "
            "nident score sstrand"
        )
        self.assertEqual(
            run_command.call_args_list,
            [
                call(
                    [
                        "custom-blastn",
                        "-query",
                        str(working_directory_path / "nucleotide.fasta"),
                        "-db",
                        str(working_directory_path / "nucleotide_database"),
                        "-out",
                        str(working_directory_path / "nucleotide.tsv"),
                        "-outfmt",
                        common_output_format,
                        "-task",
                        "blastn",
                        "-evalue",
                        "0.001",
                        "-num_threads",
                        "3",
                        "-max_target_seqs",
                        "1",
                        "-max_hsps",
                        "1",
                    ],
                    "BLAST (nucleotide)",
                ),
                call(
                    [
                        "custom-blastp",
                        "-query",
                        str(working_directory_path / "amino_acid.fasta"),
                        "-db",
                        str(working_directory_path / "amino_acid_database"),
                        "-out",
                        str(working_directory_path / "amino_acid.tsv"),
                        "-outfmt",
                        common_output_format,
                        "-task",
                        "blastp",
                        "-evalue",
                        "0.001",
                        "-num_threads",
                        "4",
                        "-max_target_seqs",
                        "1",
                        "-max_hsps",
                        "1",
                    ],
                    "BLAST (amino_acid)",
                ),
            ],
        )

    def test_representative_bsr_search_command_remains_unchanged(self) -> None:
        """Keep the three-column protein-score contract owned by BSR."""
        protein_allele = _ProteinAllele("a1", "gene_a", "gene_a_1", "MK")
        collected_locus = _CollectedLocus(
            locus_name="gene_a",
            representative_internal_identifier="a1",
            representative_allele_identifier="gene_a_1",
            total_allele_count=1,
            protein_alleles=(protein_allele,),
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            working_directory_path = Path(temporary_directory)
            (working_directory_path / "representative_hits.tsv").write_text(
                "",
                encoding="utf-8",
            )
            with (
                patch.object(bsr, "create_blast_database") as create_database,
                patch.object(bsr, "_run_command") as run_command,
            ):
                _run_blastp_analysis(
                    collected_loci=(collected_locus,),
                    all_protein_alleles=(protein_allele,),
                    working_directory_path=working_directory_path,
                    blastp_executable="custom-blastp",
                    makeblastdb_executable="custom-makeblastdb",
                    blast_evalue=0.005,
                    thread_count=5,
                )

        self.assertEqual(create_database.call_args.kwargs["database_type"], "prot")
        run_command.assert_called_once_with(
            [
                "custom-blastp",
                "-query",
                str(working_directory_path / "representatives.fasta"),
                "-db",
                str(working_directory_path / "all_valid_proteins_database"),
                "-out",
                str(working_directory_path / "representative_hits.tsv"),
                "-outfmt",
                "6 qseqid sseqid score",
                "-task",
                "blastp",
                "-evalue",
                "0.005",
                "-num_threads",
                "5",
                "-max_target_seqs",
                "1",
                "-max_hsps",
                "1",
            ],
            "BLASTp (representative BSR)",
        )


if __name__ == "__main__":
    unittest.main()
