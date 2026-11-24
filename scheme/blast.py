"""Shared operational helpers for invoking NCBI BLAST+ executables.

This module owns process execution, executable-version discovery, and the
standard ``makeblastdb`` command used by scheme analyses.  Scientific search
arguments, tabular output formats, parsing, and interpretation remain in the
analysis modules that define those contracts.
"""

import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import Literal

type BlastDatabaseType = Literal["nucl", "prot"]

__all__ = [
    "BlastDatabaseType",
    "create_blast_database",
    "get_blast_version",
    "run_blast_command",
]


def run_blast_command(
    command_arguments: Sequence[str],
    operation_name: str,
) -> None:
    """Run one BLAST+ command and raise a concise operational error.

    Args:
        command_arguments: Executable followed by its command-line arguments.
        operation_name: Human-readable operation included in error messages.

    Raises:
        ValueError: If no executable is supplied.
        RuntimeError: If the executable is absent or exits unsuccessfully.
    """
    if not command_arguments:
        raise ValueError("command_arguments must include an executable")
    command_vector = list(command_arguments)
    try:
        completed_process = subprocess.run(
            command_vector,
            check=False,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as error:
        raise RuntimeError(
            f"{operation_name} executable was not found: {command_vector[0]}"
        ) from error
    if completed_process.returncode != 0:
        error_message = (
            completed_process.stderr.strip() or completed_process.stdout.strip()
        )
        raise RuntimeError(
            f"{operation_name} failed with exit code "
            f"{completed_process.returncode}: {error_message}"
        )


def get_blast_version(blast_executable: str) -> str:
    """Return the first version line emitted by a BLAST executable.

    Raises:
        RuntimeError: If the executable is absent or version discovery fails.
    """
    try:
        completed_process = subprocess.run(
            [blast_executable, "-version"],
            check=False,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as error:
        raise RuntimeError(
            f"BLAST executable was not found: {blast_executable}"
        ) from error
    if completed_process.returncode != 0:
        error_message = (
            completed_process.stderr.strip() or completed_process.stdout.strip()
        )
        raise RuntimeError(
            f"Could not determine the version of {blast_executable}: {error_message}"
        )
    version_lines = completed_process.stdout.strip().splitlines()
    return version_lines[0] if version_lines else f"{blast_executable} version unknown"


def create_blast_database(
    *,
    input_fasta_file_path: Path,
    database_path: Path,
    database_type: BlastDatabaseType,
    makeblastdb_executable: str,
    operation_name: str,
) -> None:
    """Create one nucleotide or protein BLAST database with parsed IDs."""
    run_blast_command(
        [
            makeblastdb_executable,
            "-in",
            str(input_fasta_file_path),
            "-dbtype",
            database_type,
            "-out",
            str(database_path),
            "-parse_seqids",
        ],
        operation_name,
    )
