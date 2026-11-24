"""Shared FASTA discovery, parsing, naming, and source-signature helpers.

The intrinsic analyses treat every recognized FASTA file below an extracted
scheme directory as one locus.  This module provides the single implementation
of that file contract while deliberately leaving biological normalization and
validation to :mod:`scheme.sequence`.
"""

import gzip
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

FASTA_FILENAME_SUFFIXES = (
    ".fasta",
    ".fas",
    ".fa",
    ".fna",
    ".fasta.gz",
    ".fas.gz",
    ".fa.gz",
    ".fna.gz",
)

type FastaSourceSignature = tuple[tuple[str, int, int], ...]

__all__ = [
    "FASTA_FILENAME_SUFFIXES",
    "FastaRecord",
    "FastaSourceSignature",
    "create_fasta_source_signature",
    "find_locus_fasta_paths",
    "is_fasta_file",
    "iterate_fasta_records",
    "locus_name_from_file_path",
    "read_fasta_records",
    "remove_fasta_filename_suffix",
    "resolve_scheme_directory_path",
]


@dataclass(frozen=True)
class FastaRecord:
    """One header-defined FASTA record in original sequence letter case."""

    identifier: str
    sequence: str


def resolve_scheme_directory_path(scheme_directory_path: str | Path) -> Path:
    """Resolve and validate a directory containing scheme locus FASTA files."""
    resolved_scheme_directory_path = Path(scheme_directory_path).resolve()
    if not resolved_scheme_directory_path.exists():
        raise FileNotFoundError(
            f"Scheme directory does not exist: {resolved_scheme_directory_path}"
        )
    if not resolved_scheme_directory_path.is_dir():
        raise NotADirectoryError(
            f"Scheme path is not a directory: {resolved_scheme_directory_path}"
        )
    return resolved_scheme_directory_path


def is_fasta_file(candidate_file_path: Path) -> bool:
    """Return whether a file has a recognized FASTA suffix."""
    return candidate_file_path.is_file() and candidate_file_path.name.lower().endswith(
        FASTA_FILENAME_SUFFIXES
    )


def find_locus_fasta_paths(scheme_directory_path: Path) -> tuple[Path, ...]:
    """Return every locus FASTA recursively in deterministic path order."""
    return tuple(
        sorted(
            candidate_file_path
            for candidate_file_path in scheme_directory_path.rglob("*")
            if is_fasta_file(candidate_file_path)
        )
    )


def remove_fasta_filename_suffix(relative_file_path: str) -> str:
    """Remove the longest recognized FASTA suffix from a relative path."""
    lowercase_relative_file_path = relative_file_path.lower()
    for filename_suffix in sorted(FASTA_FILENAME_SUFFIXES, key=len, reverse=True):
        if lowercase_relative_file_path.endswith(filename_suffix):
            return relative_file_path[: -len(filename_suffix)]
    return relative_file_path


def locus_name_from_file_path(
    scheme_directory_path: Path,
    locus_file_path: Path,
) -> str:
    """Return a scheme-relative locus name without its FASTA suffix."""
    relative_file_path = locus_file_path.relative_to(
        scheme_directory_path
    ).as_posix()
    return remove_fasta_filename_suffix(relative_file_path)


def iterate_fasta_records(fasta_lines: Iterable[str]) -> Iterator[FastaRecord]:
    """Yield header-defined records while ignoring text before the first header.

    The identifier is the first whitespace-delimited header token.  Blank
    lines are ignored, header-defined empty records are retained, and sequence
    letter case is preserved for callers to normalize according to their
    scientific contract.
    """
    current_identifier: str | None = None
    current_sequence_parts: list[str] = []
    for fasta_line in fasta_lines:
        stripped_fasta_line = fasta_line.strip()
        if stripped_fasta_line.startswith(">"):
            if current_identifier is not None:
                yield FastaRecord(
                    identifier=current_identifier,
                    sequence="".join(current_sequence_parts),
                )
            current_identifier = stripped_fasta_line[1:].split(maxsplit=1)[0]
            current_sequence_parts = []
        elif stripped_fasta_line and current_identifier is not None:
            current_sequence_parts.append(stripped_fasta_line)

    if current_identifier is not None:
        yield FastaRecord(
            identifier=current_identifier,
            sequence="".join(current_sequence_parts),
        )


def read_fasta_records(fasta_file_path: Path) -> tuple[FastaRecord, ...]:
    """Read header-defined records from plain-text or gzip-compressed FASTA."""
    if fasta_file_path.name.lower().endswith(".gz"):
        with gzip.open(
            fasta_file_path,
            mode="rt",
            encoding="utf-8",
        ) as fasta_file:
            return tuple(iterate_fasta_records(fasta_file))
    with fasta_file_path.open(encoding="utf-8") as fasta_file:
        return tuple(iterate_fasta_records(fasta_file))


def create_fasta_source_signature(
    scheme_directory_path: Path,
    fasta_file_paths: Iterable[Path],
) -> FastaSourceSignature:
    """Return relative path, size, and mtime for every supplied FASTA file."""
    source_signature: list[tuple[str, int, int]] = []
    for fasta_file_path in fasta_file_paths:
        file_status = fasta_file_path.stat()
        source_signature.append(
            (
                fasta_file_path.relative_to(scheme_directory_path).as_posix(),
                file_status.st_size,
                file_status.st_mtime_ns,
            )
        )
    return tuple(source_signature)
