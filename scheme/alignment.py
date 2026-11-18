"""Align alleles within loci and report nucleotide and amino-acid gaps.

This module runs independent MAFFT multiple-sequence alignments for the
nucleotide alleles and translated protein alleles of each locus. The default
FFT-NS-2 strategy prioritizes throughput while retaining two guide-tree
rebuilds. FFT-NS-1 provides a faster, rougher guide tree, and MAFFT ``--auto``
remains available when automatic selection of a potentially more expensive
method is preferred. Nucleotide input includes every non-empty record
containing only standard IUPAC DNA
symbols. Protein input includes every allele that satisfies the repository's
existing coding-sequence validation rules. Exact duplicate records remain in
both analyses so the proportions describe the scheme as stored.

Gap proportion is the number of ``-`` characters divided by the total number
of cells in the alignment. Leading and trailing gaps count because they
capture allele-length variation. A locus with fewer than two eligible
sequences is reported as unavailable rather than as zero. An input containing
only identical eligible sequences has an exact zero-gap result and bypasses
the external aligner.

The nucleotide MSA also provides substitution-polymorphism metrics. Gaps and
ambiguous IUPAC symbols are missing observations. A column is callable when at
least two sequences provide canonical A/C/G/T bases and polymorphic when those
observations contain at least two distinct bases. The percentage uses callable
columns as its denominator, keeping indel variation in the separate gap metric.

Full alignments are not persisted. Source signatures, MAFFT provenance,
eligibility counts, and derived metrics are cached independently per locus and
alignment strategy. Unchanged successful entries are reused, failed entries
are retried, and periodic atomic checkpoints allow interrupted analyses to
resume. Version-three caches reuse their protein results while nucleotide
alignments are upgraded once with polymorphism metrics. Changing strategies
retains the other strategies' cached results.

The public functions are suitable for the Marimo notebook and command-line
use. A typed progress callback supports notebook-native status components
without importing Marimo here. The amino-acid warning threshold is applied
only to cached proportions, so changing it never reruns MAFFT.
"""

import argparse
import gzip
import json
import subprocess
import sys
import tempfile
import time
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TextIO

from .allele_info import (
    _evaluate_allele_sequence,
    _translate_valid_allele,
)
from .info import FASTA_FILENAME_SUFFIXES
from .resources import resolve_worker_count

SCHEME_ALIGNMENT_CACHE_FILENAME = ".scheme-alignment-cache.json"
LEGACY_ALIGNMENT_CACHE_FILENAME = ".schema-alignment-cache.json"
SCHEME_ALIGNMENT_CACHE_FILENAMES = (
    SCHEME_ALIGNMENT_CACHE_FILENAME,
    LEGACY_ALIGNMENT_CACHE_FILENAME,
)
SCHEME_ALIGNMENT_CACHE_VERSION = 4
SCHEME_ALIGNMENT_CACHE_CHECKPOINT_LOCUS_INTERVAL = 10
SCHEME_ALIGNMENT_CACHE_CHECKPOINT_SECONDS = 5.0
DEFAULT_AMINO_ACID_GAP_WARNING_THRESHOLD_PERCENTAGE = 5.0
IUPAC_DNA_BASES = frozenset("ACGTRYSWKMBDHVN")
ALIGNMENT_STRATEGY_FFT_NS_2 = "fft-ns-2"
ALIGNMENT_STRATEGY_FFT_NS_1 = "fft-ns-1"
ALIGNMENT_STRATEGY_AUTO = "auto"
DEFAULT_ALIGNMENT_STRATEGY = ALIGNMENT_STRATEGY_FFT_NS_2
ALIGNMENT_STRATEGY_LABELS = {
    ALIGNMENT_STRATEGY_FFT_NS_2: "FFT-NS-2 (fast, recommended)",
    ALIGNMENT_STRATEGY_FFT_NS_1: "FFT-NS-1 (fastest, rougher guide tree)",
    ALIGNMENT_STRATEGY_AUTO: "Auto (MAFFT-selected, potentially slower)",
}
ALIGNMENT_METHOD_IDENTIFIERS = {
    alignment_strategy: (
        f"mafft-{alignment_strategy}-independent-nucleotide-protein-v1"
    )
    for alignment_strategy in ALIGNMENT_STRATEGY_LABELS
}
ALIGNMENT_METHOD_IDENTIFIER = ALIGNMENT_METHOD_IDENTIFIERS[DEFAULT_ALIGNMENT_STRATEGY]
LEGACY_AUTO_ALIGNMENT_METHOD_IDENTIFIER = "mafft-auto-independent-nucleotide-protein-v1"
ALIGNMENT_STATUS_ALIGNED = "aligned"
ALIGNMENT_STATUS_INSUFFICIENT_SEQUENCES = "insufficient_sequences"
ALIGNMENT_STATUS_FAILED = "failed"

__all__ = [
    "ALIGNMENT_METHOD_IDENTIFIER",
    "ALIGNMENT_METHOD_IDENTIFIERS",
    "ALIGNMENT_STATUS_ALIGNED",
    "ALIGNMENT_STATUS_FAILED",
    "ALIGNMENT_STATUS_INSUFFICIENT_SEQUENCES",
    "ALIGNMENT_STRATEGY_AUTO",
    "ALIGNMENT_STRATEGY_FFT_NS_1",
    "ALIGNMENT_STRATEGY_FFT_NS_2",
    "ALIGNMENT_STRATEGY_LABELS",
    "DEFAULT_ALIGNMENT_STRATEGY",
    "DEFAULT_AMINO_ACID_GAP_WARNING_THRESHOLD_PERCENTAGE",
    "IUPAC_DNA_BASES",
    "SCHEME_ALIGNMENT_CACHE_FILENAME",
    "SCHEME_ALIGNMENT_CACHE_FILENAMES",
    "SCHEME_ALIGNMENT_CACHE_VERSION",
    "AlignmentProgress",
    "LocusAlignmentStatistics",
    "SchemeAlignmentAnalysis",
    "SchemeAlignmentSummary",
    "calculate_scheme_alignment_summary",
    "get_scheme_alignment_analysis",
]


@dataclass(frozen=True)
class LocusAlignmentStatistics:
    """Nucleotide and amino-acid alignment statistics for one locus.

    Attributes:
        locus_name: Scheme-relative locus name without a FASTA suffix.
        total_allele_count: Total FASTA records in the locus.
        nucleotide_aligned_sequence_count: Non-empty IUPAC DNA records passed
            to MAFFT.
        nucleotide_excluded_sequence_count: Records excluded from nucleotide
            alignment because they were empty or contained non-IUPAC symbols.
        nucleotide_alignment_status: ``aligned``, ``insufficient_sequences``,
            or ``failed``.
        nucleotide_alignment_length: Number of columns in the nucleotide MSA.
        nucleotide_gap_character_count: Total ``-`` characters in that MSA.
        nucleotide_gap_proportion_percentage: Gap characters divided by all
            nucleotide alignment cells, expressed as a percentage.
        nucleotide_callable_site_count: Alignment columns containing at least
            two canonical A/C/G/T observations. Gaps and ambiguous IUPAC
            symbols are treated as missing observations.
        nucleotide_polymorphic_site_count: Callable columns containing at
            least two distinct canonical nucleotide states.
        nucleotide_polymorphic_site_percentage: Polymorphic sites divided by
            callable sites, expressed as a percentage.
        nucleotide_alignment_error: MAFFT or alignment-parsing failure, when
            the nucleotide alignment status is ``failed``.
        amino_acid_aligned_sequence_count: Valid CDS translations passed to
            MAFFT. Alternative starts are normalized to methionine and the
            terminal stop codon is omitted.
        amino_acid_excluded_sequence_count: Records excluded because they did
            not satisfy the existing CDS validation rules.
        amino_acid_alignment_status: ``aligned``,
            ``insufficient_sequences``, or ``failed``.
        amino_acid_alignment_length: Number of columns in the protein MSA.
        amino_acid_gap_character_count: Total ``-`` characters in that MSA.
        amino_acid_gap_proportion_percentage: Gap characters divided by all
            protein alignment cells, expressed as a percentage.
        amino_acid_alignment_error: MAFFT or alignment-parsing failure, when
            the amino-acid alignment status is ``failed``.
    """

    locus_name: str
    total_allele_count: int
    nucleotide_aligned_sequence_count: int
    nucleotide_excluded_sequence_count: int
    nucleotide_alignment_status: str
    nucleotide_alignment_length: int | None
    nucleotide_gap_character_count: int
    nucleotide_gap_proportion_percentage: float | None
    nucleotide_callable_site_count: int | None
    nucleotide_polymorphic_site_count: int | None
    nucleotide_polymorphic_site_percentage: float | None
    nucleotide_alignment_error: str | None
    amino_acid_aligned_sequence_count: int
    amino_acid_excluded_sequence_count: int
    amino_acid_alignment_status: str
    amino_acid_alignment_length: int | None
    amino_acid_gap_character_count: int
    amino_acid_gap_proportion_percentage: float | None
    amino_acid_alignment_error: str | None


@dataclass(frozen=True)
class SchemeAlignmentAnalysis:
    """Cached alignment provenance and ordered per-locus statistics."""

    mafft_version: str
    alignment_method_identifier: str
    locus_statistics: tuple[LocusAlignmentStatistics, ...]
    alignment_strategy: str = DEFAULT_ALIGNMENT_STRATEGY


@dataclass(frozen=True)
class SchemeAlignmentSummary:
    """Scheme-level availability and amino-acid warning-threshold summary."""

    locus_count: int
    nucleotide_aligned_locus_count: int
    nucleotide_unavailable_locus_count: int
    nucleotide_failed_locus_count: int
    amino_acid_aligned_locus_count: int
    amino_acid_unavailable_locus_count: int
    amino_acid_failed_locus_count: int
    amino_acid_gap_warning_threshold_percentage: float
    amino_acid_gap_warning_locus_count: int
    amino_acid_gap_warning_locus_percentage: float
    amino_acid_gap_warning_loci: tuple[str, ...]


@dataclass(frozen=True)
class AlignmentProgress:
    """Progress reported while cached and uncached loci are assembled.

    Attributes:
        total_locus_count: Total loci in the current scheme.
        completed_locus_count: Cached or newly analyzed loci now available.
        cached_locus_count: Loci reused from a compatible per-locus cache.
        current_locus_name: Most recently analyzed locus, or ``None`` for the
            initial cache-reuse update.
        alignment_strategy: MAFFT strategy used for uncached work.
        elapsed_seconds: Wall-clock seconds spent processing uncached loci.
        newly_completed_loci_per_second: Completion rate for uncached loci,
            or ``None`` until a locus has completed.
    """

    total_locus_count: int
    completed_locus_count: int
    cached_locus_count: int
    current_locus_name: str | None
    alignment_strategy: str
    elapsed_seconds: float
    newly_completed_loci_per_second: float | None


@dataclass(frozen=True)
class _AlignmentMetrics:
    """Internal status and gap measurements for one MAFFT invocation."""

    status: str
    alignment_length: int | None
    gap_character_count: int
    gap_proportion_percentage: float | None
    error: str | None
    nucleotide_callable_site_count: int | None = None
    nucleotide_polymorphic_site_count: int | None = None
    nucleotide_polymorphic_site_percentage: float | None = None


@dataclass(frozen=True)
class _CachedLocusAlignment:
    """Source signature and derived alignment metrics for one locus."""

    file_size: int
    modification_time_ns: int
    locus_statistics: LocusAlignmentStatistics
    nucleotide_polymorphism_metrics_calculated: bool


@dataclass(frozen=True)
class _CachedAlignmentMethod:
    """MAFFT provenance and per-locus entries for one strategy."""

    mafft_version: str
    locus_entries: dict[str, _CachedLocusAlignment]


@dataclass(frozen=True)
class _AlignmentCache:
    """Parsed versioned cache containing independent strategy namespaces."""

    cache_version: int
    method_caches: dict[str, _CachedAlignmentMethod]


class _MafftProgressReporter:
    """Write throttled scheme-level MAFFT progress to a text stream."""

    _BAR_WIDTH = 30

    def __init__(
        self,
        total_locus_count: int,
        initial_completed_locus_count: int,
        worker_count: int,
        alignment_strategy: str,
        output_stream: TextIO | None,
    ) -> None:
        self.total_locus_count = total_locus_count
        self.initial_completed_locus_count = initial_completed_locus_count
        self.worker_count = worker_count
        self.alignment_strategy = alignment_strategy
        self.output_stream = output_stream
        self.completed_locus_count = initial_completed_locus_count
        self.started_at = time.monotonic()
        self.output_is_terminal = bool(
            output_stream is not None
            and getattr(output_stream, "isatty", lambda: False)()
        )
        self.non_terminal_update_interval = max(1, total_locus_count // 20)

    def start(self) -> None:
        """Show the initial progress state for a non-empty analysis."""
        if self.output_stream is None or self.total_locus_count == 0:
            return
        self._write_progress(final_update=False)

    def advance(self) -> None:
        """Record one locus whose nucleotide and protein work has finished."""
        self.completed_locus_count += 1
        if self.output_stream is None:
            return
        should_write_non_terminal_update = (
            self.completed_locus_count == self.total_locus_count
            or self.completed_locus_count % self.non_terminal_update_interval == 0
        )
        if self.output_is_terminal or should_write_non_terminal_update:
            self._write_progress(
                final_update=(self.completed_locus_count == self.total_locus_count)
            )

    def finish(self) -> None:
        """Terminate an in-place terminal progress line after interruption."""
        if (
            self.output_stream is not None
            and self.output_is_terminal
            and self.completed_locus_count < self.total_locus_count
        ):
            self.output_stream.write("\n")
            self.output_stream.flush()

    def _write_progress(self, *, final_update: bool) -> None:
        """Render the current count, elapsed time, rate, and estimated time."""
        if self.output_stream is None:
            return
        completed_fraction = self.completed_locus_count / self.total_locus_count
        completed_bar_width = round(completed_fraction * self._BAR_WIDTH)
        progress_bar = "#" * completed_bar_width + "-" * (
            self._BAR_WIDTH - completed_bar_width
        )
        elapsed_seconds = max(0.0, time.monotonic() - self.started_at)
        newly_completed_locus_count = (
            self.completed_locus_count - self.initial_completed_locus_count
        )
        loci_per_second = (
            newly_completed_locus_count / elapsed_seconds
            if elapsed_seconds > 0 and newly_completed_locus_count
            else 0.0
        )
        remaining_locus_count = self.total_locus_count - self.completed_locus_count
        estimated_remaining_seconds = (
            remaining_locus_count / loci_per_second if loci_per_second else None
        )
        progress_text = (
            f"MAFFT alignments [{progress_bar}] "
            f"{self.completed_locus_count}/{self.total_locus_count} loci "
            f"({completed_fraction * 100:5.1f}%) | "
            f"{ALIGNMENT_STRATEGY_LABELS[self.alignment_strategy]} | "
            f"{self.worker_count} workers | "
            f"{self.initial_completed_locus_count} cached | "
            f"elapsed {_format_duration(elapsed_seconds)}"
        )
        if loci_per_second:
            progress_text += f" | {loci_per_second:.2f} loci/s"
        if estimated_remaining_seconds is not None and not final_update:
            progress_text += f" | ETA {_format_duration(estimated_remaining_seconds)}"
        if final_update:
            progress_text += " | complete"

        if self.output_is_terminal:
            line_ending = "\n" if final_update else ""
            self.output_stream.write(f"\r{progress_text}\033[K{line_ending}")
        else:
            self.output_stream.write(f"{progress_text}\n")
        self.output_stream.flush()


def _format_duration(duration_seconds: float) -> str:
    """Format a progress duration as seconds or hours, minutes, and seconds."""
    rounded_duration_seconds = max(0, round(duration_seconds))
    hours, remaining_seconds = divmod(rounded_duration_seconds, 3600)
    minutes, seconds = divmod(remaining_seconds, 60)
    if hours:
        return f"{hours:d}h{minutes:02d}m{seconds:02d}s"
    if minutes:
        return f"{minutes:d}m{seconds:02d}s"
    return f"{seconds:d}s"


def _resolve_scheme_directory_path(scheme_directory_path: str | Path) -> Path:
    """Resolve and validate an extracted scheme directory."""
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


def _is_fasta_file(candidate_file_path: Path) -> bool:
    """Return whether a path is a recognized FASTA file."""
    return candidate_file_path.is_file() and candidate_file_path.name.lower().endswith(
        FASTA_FILENAME_SUFFIXES
    )


def _find_allele_file_paths(
    scheme_directory_path: Path,
) -> tuple[Path, ...]:
    """Return recognized locus FASTA files in deterministic order."""
    return tuple(
        sorted(
            candidate_file_path
            for candidate_file_path in scheme_directory_path.rglob("*")
            if _is_fasta_file(candidate_file_path)
        )
    )


def _remove_fasta_filename_suffix(relative_file_path: str) -> str:
    """Remove a recognized FASTA suffix from a scheme-relative path."""
    lowercase_relative_file_path = relative_file_path.lower()
    for filename_suffix in sorted(FASTA_FILENAME_SUFFIXES, key=len, reverse=True):
        if lowercase_relative_file_path.endswith(filename_suffix):
            return relative_file_path[: -len(filename_suffix)]
    return relative_file_path


def _iterate_fasta_sequences(fasta_lines: Iterable[str]) -> Iterable[str]:
    """Yield assembled sequences, including empty header-defined records."""
    current_sequence_parts: list[str] | None = None
    for fasta_line in fasta_lines:
        stripped_fasta_line = fasta_line.strip()
        if stripped_fasta_line.startswith(">"):
            if current_sequence_parts is not None:
                yield "".join(current_sequence_parts)
            current_sequence_parts = []
        elif stripped_fasta_line and current_sequence_parts is not None:
            current_sequence_parts.append(stripped_fasta_line)

    if current_sequence_parts is not None:
        yield "".join(current_sequence_parts)


def _read_fasta_sequences(allele_file_path: Path) -> tuple[str, ...]:
    """Read all assembled allele sequences from one locus FASTA file."""
    if allele_file_path.name.lower().endswith(".gz"):
        with gzip.open(
            allele_file_path,
            mode="rt",
            encoding="utf-8",
        ) as allele_file:
            return tuple(_iterate_fasta_sequences(allele_file))

    with allele_file_path.open(encoding="utf-8") as allele_file:
        return tuple(_iterate_fasta_sequences(allele_file))


def _create_source_signature(
    scheme_directory_path: Path,
    allele_file_paths: Iterable[Path],
) -> tuple[tuple[str, int, int], ...]:
    """Create the FASTA metadata signature used for cache invalidation."""
    source_signature: list[tuple[str, int, int]] = []
    for allele_file_path in allele_file_paths:
        allele_file_status = allele_file_path.stat()
        source_signature.append(
            (
                allele_file_path.relative_to(scheme_directory_path).as_posix(),
                allele_file_status.st_size,
                allele_file_status.st_mtime_ns,
            )
        )
    return tuple(source_signature)


def _get_mafft_version(mafft_executable: str) -> str:
    """Return the MAFFT version or raise a descriptive runtime error."""
    try:
        completed_process = subprocess.run(
            [mafft_executable, "--version"],
            check=False,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as error:
        raise RuntimeError(
            f"MAFFT executable was not found: {mafft_executable}. "
            "Enter the devenv shell after installing the updated dependencies."
        ) from error

    version_output = (
        completed_process.stdout.strip() or completed_process.stderr.strip()
    )
    if completed_process.returncode != 0:
        raise RuntimeError(
            "Unable to determine the MAFFT version: "
            f"{version_output or 'no diagnostic output'}"
        )
    return version_output or "unknown"


def _serialize_fasta_sequences(sequences: tuple[str, ...]) -> str:
    """Serialize sequences with safe synthetic identifiers for MAFFT."""
    return "".join(
        f">sequence_{sequence_index}\n{sequence}\n"
        for sequence_index, sequence in enumerate(sequences, start=1)
    )


def _parse_aligned_fasta(aligned_fasta_text: str) -> tuple[str, ...]:
    """Parse assembled aligned sequences from MAFFT FASTA output."""
    aligned_sequences = tuple(_iterate_fasta_sequences(aligned_fasta_text.splitlines()))
    if not aligned_sequences:
        raise ValueError("MAFFT returned no aligned sequences")
    alignment_lengths = {len(sequence) for sequence in aligned_sequences}
    if len(alignment_lengths) != 1:
        raise ValueError("MAFFT returned sequences with inconsistent alignment lengths")
    return aligned_sequences


def _calculate_gap_metrics_from_aligned_sequences(
    aligned_sequences: tuple[str, ...],
    *,
    calculate_nucleotide_polymorphism: bool = False,
) -> _AlignmentMetrics:
    """Calculate gap and optional substitution-polymorphism MSA metrics."""
    alignment_length = len(aligned_sequences[0])
    gap_character_count = sum(
        aligned_sequence.count("-") for aligned_sequence in aligned_sequences
    )
    alignment_cell_count = len(aligned_sequences) * alignment_length
    gap_proportion_percentage = (
        gap_character_count / alignment_cell_count * 100
        if alignment_cell_count
        else 0.0
    )
    nucleotide_callable_site_count: int | None = None
    nucleotide_polymorphic_site_count: int | None = None
    nucleotide_polymorphic_site_percentage: float | None = None
    if calculate_nucleotide_polymorphism:
        (
            nucleotide_callable_site_count,
            nucleotide_polymorphic_site_count,
            nucleotide_polymorphic_site_percentage,
        ) = _calculate_nucleotide_polymorphism_metrics(aligned_sequences)
    return _AlignmentMetrics(
        status=ALIGNMENT_STATUS_ALIGNED,
        alignment_length=alignment_length,
        gap_character_count=gap_character_count,
        gap_proportion_percentage=gap_proportion_percentage,
        error=None,
        nucleotide_callable_site_count=nucleotide_callable_site_count,
        nucleotide_polymorphic_site_count=nucleotide_polymorphic_site_count,
        nucleotide_polymorphic_site_percentage=(
            nucleotide_polymorphic_site_percentage
        ),
    )


def _calculate_nucleotide_polymorphism_metrics(
    aligned_sequences: tuple[str, ...],
) -> tuple[int, int, float | None]:
    """Count callable and substitution-polymorphic nucleotide columns.

    Gaps and ambiguous IUPAC symbols are missing observations. A column is
    callable when it contains at least two canonical A/C/G/T observations and
    polymorphic when those observations include at least two distinct states.
    """
    canonical_nucleotides = frozenset("ACGT")
    callable_site_count = 0
    polymorphic_site_count = 0
    for alignment_column in zip(*aligned_sequences, strict=True):
        canonical_observations = tuple(
            nucleotide.upper()
            for nucleotide in alignment_column
            if nucleotide.upper() in canonical_nucleotides
        )
        if len(canonical_observations) < 2:
            continue
        callable_site_count += 1
        if len(set(canonical_observations)) >= 2:
            polymorphic_site_count += 1

    polymorphic_site_percentage = (
        polymorphic_site_count / callable_site_count * 100
        if callable_site_count
        else None
    )
    return (
        callable_site_count,
        polymorphic_site_count,
        polymorphic_site_percentage,
    )


def _resolve_alignment_strategy(alignment_strategy: str) -> str:
    """Validate and return a supported MAFFT alignment strategy."""
    if alignment_strategy not in ALIGNMENT_STRATEGY_LABELS:
        supported_strategies = ", ".join(ALIGNMENT_STRATEGY_LABELS)
        raise ValueError(f"alignment_strategy must be one of: {supported_strategies}")
    return alignment_strategy


def _get_mafft_strategy_arguments(alignment_strategy: str) -> tuple[str, ...]:
    """Return the MAFFT options implementing one documented strategy."""
    resolved_alignment_strategy = _resolve_alignment_strategy(alignment_strategy)
    if resolved_alignment_strategy == ALIGNMENT_STRATEGY_AUTO:
        return ("--auto",)
    if resolved_alignment_strategy == ALIGNMENT_STRATEGY_FFT_NS_1:
        return ("--retree", "1", "--maxiterate", "0")
    return ("--retree", "2", "--maxiterate", "0")


def _run_mafft_alignment(
    sequences: tuple[str, ...],
    sequence_type: str,
    mafft_executable: str,
    alignment_strategy: str = DEFAULT_ALIGNMENT_STRATEGY,
) -> _AlignmentMetrics:
    """Run one MAFFT alignment and return gap metrics or a failure result."""
    strategy_arguments = _get_mafft_strategy_arguments(alignment_strategy)
    if len(sequences) < 2:
        return _AlignmentMetrics(
            status=ALIGNMENT_STATUS_INSUFFICIENT_SEQUENCES,
            alignment_length=None,
            gap_character_count=0,
            gap_proportion_percentage=None,
            error=None,
        )

    if all(sequence == sequences[0] for sequence in sequences[1:]):
        identical_alignment_metrics = _calculate_gap_metrics_from_aligned_sequences(
            sequences,
            calculate_nucleotide_polymorphism=(
                sequence_type == "nucleotide"
            ),
        )
        return _AlignmentMetrics(
            status=ALIGNMENT_STATUS_ALIGNED,
            alignment_length=identical_alignment_metrics.alignment_length,
            gap_character_count=identical_alignment_metrics.gap_character_count,
            gap_proportion_percentage=(
                identical_alignment_metrics.gap_proportion_percentage
            ),
            error=None,
            nucleotide_callable_site_count=(
                identical_alignment_metrics.nucleotide_callable_site_count
            ),
            nucleotide_polymorphic_site_count=(
                identical_alignment_metrics.nucleotide_polymorphic_site_count
            ),
            nucleotide_polymorphic_site_percentage=(
                identical_alignment_metrics.nucleotide_polymorphic_site_percentage
            ),
        )

    sequence_type_option = "--nuc" if sequence_type == "nucleotide" else "--amino"
    try:
        completed_process = subprocess.run(
            [
                mafft_executable,
                "--quiet",
                *strategy_arguments,
                "--thread",
                "1",
                sequence_type_option,
                "-",
            ],
            input=_serialize_fasta_sequences(sequences),
            check=False,
            capture_output=True,
            text=True,
        )
        if completed_process.returncode != 0:
            diagnostic_output = completed_process.stderr.strip()
            raise RuntimeError(
                diagnostic_output or "MAFFT exited without diagnostic output"
            )
        aligned_sequences = _parse_aligned_fasta(completed_process.stdout)
        if len(aligned_sequences) != len(sequences):
            raise ValueError(
                "MAFFT returned a different number of sequences than it received"
            )
        return _calculate_gap_metrics_from_aligned_sequences(
            aligned_sequences,
            calculate_nucleotide_polymorphism=(sequence_type == "nucleotide"),
        )
    except (OSError, RuntimeError, ValueError) as error:
        return _AlignmentMetrics(
            status=ALIGNMENT_STATUS_FAILED,
            alignment_length=None,
            gap_character_count=0,
            gap_proportion_percentage=None,
            error=str(error),
        )


def _restore_cached_alignment_metrics(
    locus_statistics: LocusAlignmentStatistics,
    sequence_type: str,
) -> _AlignmentMetrics:
    """Restore one sequence type's metrics for partial cache reuse."""
    if sequence_type == "nucleotide":
        return _AlignmentMetrics(
            status=locus_statistics.nucleotide_alignment_status,
            alignment_length=locus_statistics.nucleotide_alignment_length,
            gap_character_count=locus_statistics.nucleotide_gap_character_count,
            gap_proportion_percentage=(
                locus_statistics.nucleotide_gap_proportion_percentage
            ),
            error=locus_statistics.nucleotide_alignment_error,
            nucleotide_callable_site_count=(
                locus_statistics.nucleotide_callable_site_count
            ),
            nucleotide_polymorphic_site_count=(
                locus_statistics.nucleotide_polymorphic_site_count
            ),
            nucleotide_polymorphic_site_percentage=(
                locus_statistics.nucleotide_polymorphic_site_percentage
            ),
        )
    return _AlignmentMetrics(
        status=locus_statistics.amino_acid_alignment_status,
        alignment_length=locus_statistics.amino_acid_alignment_length,
        gap_character_count=locus_statistics.amino_acid_gap_character_count,
        gap_proportion_percentage=(
            locus_statistics.amino_acid_gap_proportion_percentage
        ),
        error=locus_statistics.amino_acid_alignment_error,
    )


def _analyze_locus(
    scheme_directory_path: Path,
    allele_file_path: Path,
    mafft_executable: str,
    alignment_strategy: str = DEFAULT_ALIGNMENT_STRATEGY,
    cached_locus_alignment: _CachedLocusAlignment | None = None,
) -> LocusAlignmentStatistics:
    """Align required sequence types, reusing compatible partial results."""
    allele_sequences = _read_fasta_sequences(allele_file_path)
    normalized_allele_sequences = tuple(
        "".join(allele_sequence.split()).upper() for allele_sequence in allele_sequences
    )
    nucleotide_sequences = tuple(
        normalized_allele_sequence
        for normalized_allele_sequence in normalized_allele_sequences
        if normalized_allele_sequence
        and set(normalized_allele_sequence).issubset(IUPAC_DNA_BASES)
    )
    amino_acid_sequences = tuple(
        _translate_valid_allele(normalized_allele_sequence)
        for normalized_allele_sequence in normalized_allele_sequences
        if _evaluate_allele_sequence(normalized_allele_sequence).is_valid
    )
    cached_locus_statistics = (
        cached_locus_alignment.locus_statistics
        if cached_locus_alignment is not None
        else None
    )
    can_reuse_nucleotide_alignment = bool(
        cached_locus_alignment is not None
        and cached_locus_statistics is not None
        and cached_locus_statistics.nucleotide_alignment_status
        != ALIGNMENT_STATUS_FAILED
        and (
            cached_locus_alignment.nucleotide_polymorphism_metrics_calculated
            or cached_locus_statistics.nucleotide_alignment_status
            == ALIGNMENT_STATUS_INSUFFICIENT_SEQUENCES
        )
    )
    if can_reuse_nucleotide_alignment and cached_locus_statistics is not None:
        nucleotide_alignment_metrics = _restore_cached_alignment_metrics(
            cached_locus_statistics,
            sequence_type="nucleotide",
        )
    else:
        nucleotide_alignment_metrics = _run_mafft_alignment(
            nucleotide_sequences,
            sequence_type="nucleotide",
            mafft_executable=mafft_executable,
            alignment_strategy=alignment_strategy,
        )

    can_reuse_amino_acid_alignment = bool(
        cached_locus_statistics is not None
        and cached_locus_statistics.amino_acid_alignment_status
        != ALIGNMENT_STATUS_FAILED
    )
    if can_reuse_amino_acid_alignment and cached_locus_statistics is not None:
        amino_acid_alignment_metrics = _restore_cached_alignment_metrics(
            cached_locus_statistics,
            sequence_type="amino_acid",
        )
    else:
        amino_acid_alignment_metrics = _run_mafft_alignment(
            amino_acid_sequences,
            sequence_type="amino_acid",
            mafft_executable=mafft_executable,
            alignment_strategy=alignment_strategy,
        )
    locus_name = _remove_fasta_filename_suffix(
        allele_file_path.relative_to(scheme_directory_path).as_posix()
    )
    return LocusAlignmentStatistics(
        locus_name=locus_name,
        total_allele_count=len(normalized_allele_sequences),
        nucleotide_aligned_sequence_count=len(nucleotide_sequences),
        nucleotide_excluded_sequence_count=(
            len(normalized_allele_sequences) - len(nucleotide_sequences)
        ),
        nucleotide_alignment_status=nucleotide_alignment_metrics.status,
        nucleotide_alignment_length=nucleotide_alignment_metrics.alignment_length,
        nucleotide_gap_character_count=(
            nucleotide_alignment_metrics.gap_character_count
        ),
        nucleotide_gap_proportion_percentage=(
            nucleotide_alignment_metrics.gap_proportion_percentage
        ),
        nucleotide_callable_site_count=(
            nucleotide_alignment_metrics.nucleotide_callable_site_count
        ),
        nucleotide_polymorphic_site_count=(
            nucleotide_alignment_metrics.nucleotide_polymorphic_site_count
        ),
        nucleotide_polymorphic_site_percentage=(
            nucleotide_alignment_metrics.nucleotide_polymorphic_site_percentage
        ),
        nucleotide_alignment_error=nucleotide_alignment_metrics.error,
        amino_acid_aligned_sequence_count=len(amino_acid_sequences),
        amino_acid_excluded_sequence_count=(
            len(normalized_allele_sequences) - len(amino_acid_sequences)
        ),
        amino_acid_alignment_status=amino_acid_alignment_metrics.status,
        amino_acid_alignment_length=amino_acid_alignment_metrics.alignment_length,
        amino_acid_gap_character_count=(
            amino_acid_alignment_metrics.gap_character_count
        ),
        amino_acid_gap_proportion_percentage=(
            amino_acid_alignment_metrics.gap_proportion_percentage
        ),
        amino_acid_alignment_error=amino_acid_alignment_metrics.error,
    )


def _optional_integer(value: object) -> int | None:
    """Convert a cached optional value to an integer."""
    return None if value is None else int(value)


def _optional_float(value: object) -> float | None:
    """Convert a cached optional value to a float."""
    return None if value is None else float(value)


def _optional_string(value: object) -> str | None:
    """Convert a cached optional value to a string."""
    return None if value is None else str(value)


def _deserialize_locus_statistics(
    locus_data: dict[str, object],
) -> LocusAlignmentStatistics:
    """Create typed locus statistics from validated cache-shaped data."""
    return LocusAlignmentStatistics(
        locus_name=str(locus_data["locus_name"]),
        total_allele_count=int(locus_data["total_allele_count"]),
        nucleotide_aligned_sequence_count=int(
            locus_data["nucleotide_aligned_sequence_count"]
        ),
        nucleotide_excluded_sequence_count=int(
            locus_data["nucleotide_excluded_sequence_count"]
        ),
        nucleotide_alignment_status=str(locus_data["nucleotide_alignment_status"]),
        nucleotide_alignment_length=_optional_integer(
            locus_data["nucleotide_alignment_length"]
        ),
        nucleotide_gap_character_count=int(
            locus_data["nucleotide_gap_character_count"]
        ),
        nucleotide_gap_proportion_percentage=_optional_float(
            locus_data["nucleotide_gap_proportion_percentage"]
        ),
        nucleotide_callable_site_count=_optional_integer(
            locus_data.get("nucleotide_callable_site_count")
        ),
        nucleotide_polymorphic_site_count=_optional_integer(
            locus_data.get("nucleotide_polymorphic_site_count")
        ),
        nucleotide_polymorphic_site_percentage=_optional_float(
            locus_data.get("nucleotide_polymorphic_site_percentage")
        ),
        nucleotide_alignment_error=_optional_string(
            locus_data["nucleotide_alignment_error"]
        ),
        amino_acid_aligned_sequence_count=int(
            locus_data["amino_acid_aligned_sequence_count"]
        ),
        amino_acid_excluded_sequence_count=int(
            locus_data["amino_acid_excluded_sequence_count"]
        ),
        amino_acid_alignment_status=str(locus_data["amino_acid_alignment_status"]),
        amino_acid_alignment_length=_optional_integer(
            locus_data["amino_acid_alignment_length"]
        ),
        amino_acid_gap_character_count=int(
            locus_data["amino_acid_gap_character_count"]
        ),
        amino_acid_gap_proportion_percentage=_optional_float(
            locus_data["amino_acid_gap_proportion_percentage"]
        ),
        amino_acid_alignment_error=_optional_string(
            locus_data["amino_acid_alignment_error"]
        ),
    )


def _read_cached_alignment_cache(
    scheme_directory_path: Path,
) -> _AlignmentCache | None:
    """Read the strategy cache or migrate compatible version-one/two data."""
    cache_file_path = next(
        (
            scheme_directory_path / cache_filename
            for cache_filename in SCHEME_ALIGNMENT_CACHE_FILENAMES
            if (scheme_directory_path / cache_filename).is_file()
        ),
        scheme_directory_path / SCHEME_ALIGNMENT_CACHE_FILENAME,
    )
    try:
        cached_data = json.loads(cache_file_path.read_text(encoding="utf-8"))
        cache_version = int(cached_data["cache_version"])
        if cache_version in (3, SCHEME_ALIGNMENT_CACHE_VERSION):
            method_caches = {
                str(alignment_method_identifier): _CachedAlignmentMethod(
                    mafft_version=str(dict(method_data)["mafft_version"]),
                    locus_entries=_deserialize_cached_locus_entries(
                        dict(dict(method_data)["loci"]),
                        cache_version=cache_version,
                    ),
                )
                for alignment_method_identifier, method_data in dict(
                    cached_data["methods"]
                ).items()
                if alignment_method_identifier in ALIGNMENT_METHOD_IDENTIFIERS.values()
            }
            return _AlignmentCache(
                cache_version=cache_version,
                method_caches=method_caches,
            )

        alignment_method_identifier = str(cached_data["alignment_method_identifier"])
        if alignment_method_identifier != LEGACY_AUTO_ALIGNMENT_METHOD_IDENTIFIER:
            return None
        mafft_version = str(cached_data["mafft_version"])
        if cache_version == 2:
            locus_entries = _deserialize_cached_locus_entries(
                dict(cached_data["loci"]),
                cache_version=cache_version,
            )
        elif cache_version == 1:
            locus_entries: dict[str, _CachedLocusAlignment] = {}
            source_signature = tuple(cached_data["source_signature"])
            locus_statistics_data = tuple(cached_data["locus_statistics"])
            if len(source_signature) != len(locus_statistics_data):
                return None
            for signature_data, locus_data in zip(
                source_signature,
                locus_statistics_data,
                strict=True,
            ):
                relative_path, file_size, modification_time_ns = signature_data
                locus_entries[str(relative_path)] = _CachedLocusAlignment(
                    file_size=int(file_size),
                    modification_time_ns=int(modification_time_ns),
                    locus_statistics=_deserialize_locus_statistics(dict(locus_data)),
                    nucleotide_polymorphism_metrics_calculated=False,
                )
        else:
            return None
        return _AlignmentCache(
            cache_version=cache_version,
            method_caches={
                LEGACY_AUTO_ALIGNMENT_METHOD_IDENTIFIER: _CachedAlignmentMethod(
                    mafft_version=mafft_version,
                    locus_entries=locus_entries,
                )
            },
        )
    except KeyError, OSError, TypeError, ValueError, json.JSONDecodeError:
        return None


def _deserialize_cached_locus_entries(
    locus_entries_data: dict[object, object],
    *,
    cache_version: int,
) -> dict[str, _CachedLocusAlignment]:
    """Deserialize path-keyed locus entries from a strategy cache."""
    locus_entries: dict[str, _CachedLocusAlignment] = {}
    for relative_path, locus_entry_data in locus_entries_data.items():
        parsed_locus_entry_data = dict(locus_entry_data)
        locus_entries[str(relative_path)] = _CachedLocusAlignment(
            file_size=int(parsed_locus_entry_data["file_size"]),
            modification_time_ns=int(parsed_locus_entry_data["modification_time_ns"]),
            locus_statistics=_deserialize_locus_statistics(
                dict(parsed_locus_entry_data["statistics"])
            ),
            nucleotide_polymorphism_metrics_calculated=bool(
                parsed_locus_entry_data.get(
                    "nucleotide_polymorphism_metrics_calculated",
                    cache_version >= SCHEME_ALIGNMENT_CACHE_VERSION,
                )
            ),
        )
    return locus_entries


def _write_cached_alignment_cache(
    scheme_directory_path: Path,
    method_caches: dict[str, _CachedAlignmentMethod],
) -> None:
    """Atomically persist all strategy-specific per-locus metrics."""
    cached_data = {
        "cache_version": SCHEME_ALIGNMENT_CACHE_VERSION,
        "methods": {
            alignment_method_identifier: {
                "mafft_version": method_caches[
                    alignment_method_identifier
                ].mafft_version,
                "loci": {
                    relative_path: {
                        "file_size": cached_locus_alignment.file_size,
                        "modification_time_ns": (
                            cached_locus_alignment.modification_time_ns
                        ),
                        "statistics": asdict(cached_locus_alignment.locus_statistics),
                        "nucleotide_polymorphism_metrics_calculated": (
                            cached_locus_alignment.nucleotide_polymorphism_metrics_calculated
                        ),
                    }
                    for relative_path, cached_locus_alignment in sorted(
                        method_caches[alignment_method_identifier].locus_entries.items()
                    )
                },
            }
            for alignment_method_identifier in sorted(method_caches)
        },
    }
    temporary_cache_file_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=scheme_directory_path,
            prefix=f"{SCHEME_ALIGNMENT_CACHE_FILENAME}.",
            delete=False,
        ) as temporary_cache_file:
            json.dump(cached_data, temporary_cache_file, separators=(",", ":"))
            temporary_cache_file_path = Path(temporary_cache_file.name)

        temporary_cache_file_path.replace(
            scheme_directory_path / SCHEME_ALIGNMENT_CACHE_FILENAME
        )
    finally:
        if temporary_cache_file_path is not None:
            temporary_cache_file_path.unlink(missing_ok=True)


def _try_write_cached_alignment_cache(
    scheme_directory_path: Path,
    method_caches: dict[str, _CachedAlignmentMethod],
) -> None:
    """Persist a checkpoint without making a read-only directory fatal."""
    try:
        _write_cached_alignment_cache(
            scheme_directory_path,
            method_caches,
        )
    except OSError:
        pass


def _is_reusable_cached_locus(
    cached_locus_alignment: _CachedLocusAlignment,
    file_size: int,
    modification_time_ns: int,
) -> bool:
    """Return whether an unchanged locus has non-failed cached measurements."""
    locus_statistics = cached_locus_alignment.locus_statistics
    return (
        _is_source_compatible_cached_locus(
            cached_locus_alignment,
            file_size,
            modification_time_ns,
        )
        and locus_statistics.nucleotide_alignment_status != ALIGNMENT_STATUS_FAILED
        and locus_statistics.amino_acid_alignment_status != ALIGNMENT_STATUS_FAILED
        and (
            cached_locus_alignment.nucleotide_polymorphism_metrics_calculated
            or locus_statistics.nucleotide_alignment_status
            == ALIGNMENT_STATUS_INSUFFICIENT_SEQUENCES
        )
    )


def _is_source_compatible_cached_locus(
    cached_locus_alignment: _CachedLocusAlignment,
    file_size: int,
    modification_time_ns: int,
) -> bool:
    """Return whether cached metrics describe the unchanged locus file."""
    return (
        cached_locus_alignment.file_size == file_size
        and cached_locus_alignment.modification_time_ns == modification_time_ns
    )


def _mark_nucleotide_polymorphism_metrics_calculated(
    cached_locus_alignment: _CachedLocusAlignment,
) -> _CachedLocusAlignment:
    """Normalize a reusable unavailable result during cache migration."""
    if cached_locus_alignment.nucleotide_polymorphism_metrics_calculated:
        return cached_locus_alignment
    return _CachedLocusAlignment(
        file_size=cached_locus_alignment.file_size,
        modification_time_ns=cached_locus_alignment.modification_time_ns,
        locus_statistics=cached_locus_alignment.locus_statistics,
        nucleotide_polymorphism_metrics_calculated=True,
    )


def _resolve_maximum_parallel_loci(
    maximum_parallel_loci: int | None,
) -> int:
    """Resolve and validate the number of concurrent single-threaded jobs."""
    if maximum_parallel_loci is not None:
        if maximum_parallel_loci < 1:
            raise ValueError("maximum_parallel_loci must be at least one")
        return maximum_parallel_loci
    return resolve_worker_count(None)


def get_scheme_alignment_analysis(
    scheme_directory_path: str | Path,
    *,
    mafft_executable: str = "mafft",
    alignment_strategy: str = DEFAULT_ALIGNMENT_STRATEGY,
    maximum_parallel_loci: int | None = None,
    progress_output_stream: TextIO | None = sys.stdout,
    progress_callback: Callable[[AlignmentProgress], None] | None = None,
) -> SchemeAlignmentAnalysis:
    """Return cached or newly calculated within-locus alignment statistics.

    Unchanged, non-failed loci are reused from the selected strategy's
    per-locus cache. Missing, changed, and previously failed loci run
    independent nucleotide and protein alignments. Each MAFFT process uses one
    thread; ``maximum_parallel_loci`` bounds the number of loci analyzed
    concurrently.

    Args:
        scheme_directory_path: Directory containing an extracted MLST scheme.
        mafft_executable: MAFFT executable name or path.
        alignment_strategy: ``fft-ns-2`` (the fast default), ``fft-ns-1``
            (the fastest and roughest option), or ``auto`` (MAFFT-selected).
        maximum_parallel_loci: Maximum concurrent locus analyses. The default
            uses 80% of CPUs available to the current process, rounded down to
            an even count when possible and with a one-worker minimum.
        progress_output_stream: Stream receiving completion-based MAFFT
            progress on a cache miss. Pass ``None`` to disable progress.
        progress_callback: Optional function called first with cache-reuse
            progress and then after each newly analyzed locus. This keeps the
            library independent from notebook-specific UI components.

    Returns:
        Alignment provenance and statistics ordered by locus path.

    Raises:
        FileNotFoundError: If the scheme directory does not exist.
        NotADirectoryError: If the scheme path is not a directory.
        RuntimeError: If MAFFT is unavailable or its version cannot be read.
        ValueError: If the strategy is unsupported or
            ``maximum_parallel_loci`` is less than one.
    """
    resolved_alignment_strategy = _resolve_alignment_strategy(alignment_strategy)
    alignment_method_identifier = ALIGNMENT_METHOD_IDENTIFIERS[
        resolved_alignment_strategy
    ]
    resolved_scheme_directory_path = _resolve_scheme_directory_path(
        scheme_directory_path
    )
    allele_file_paths = _find_allele_file_paths(resolved_scheme_directory_path)
    source_signature = _create_source_signature(
        resolved_scheme_directory_path,
        allele_file_paths,
    )
    source_signature_by_relative_path = {
        relative_path: (file_size, modification_time_ns)
        for relative_path, file_size, modification_time_ns in source_signature
    }
    allele_file_path_by_relative_path = {
        allele_file_path.relative_to(
            resolved_scheme_directory_path
        ).as_posix(): allele_file_path
        for allele_file_path in allele_file_paths
    }
    cached_alignment_cache = _read_cached_alignment_cache(
        resolved_scheme_directory_path
    )
    method_caches = (
        dict(cached_alignment_cache.method_caches)
        if cached_alignment_cache is not None
        else {}
    )
    selected_method_cache = method_caches.get(alignment_method_identifier)
    source_compatible_locus_entries: dict[str, _CachedLocusAlignment] = {}
    reusable_locus_entries: dict[str, _CachedLocusAlignment] = {}
    if selected_method_cache is not None:
        for relative_path, (
            file_size,
            modification_time_ns,
        ) in source_signature_by_relative_path.items():
            cached_locus_alignment = selected_method_cache.locus_entries.get(
                relative_path
            )
            if cached_locus_alignment is None:
                continue
            if _is_source_compatible_cached_locus(
                cached_locus_alignment,
                file_size,
                modification_time_ns,
            ):
                source_compatible_locus_entries[relative_path] = (
                    cached_locus_alignment
                )
            if _is_reusable_cached_locus(
                cached_locus_alignment,
                file_size,
                modification_time_ns,
            ):
                reusable_locus_entries[relative_path] = (
                    _mark_nucleotide_polymorphism_metrics_calculated(
                        cached_locus_alignment
                    )
                )

    total_locus_count = len(allele_file_paths)
    if (
        selected_method_cache is not None
        and len(reusable_locus_entries) == total_locus_count
    ):
        if progress_callback is not None:
            progress_callback(
                AlignmentProgress(
                    total_locus_count=total_locus_count,
                    completed_locus_count=total_locus_count,
                    cached_locus_count=total_locus_count,
                    current_locus_name=None,
                    alignment_strategy=resolved_alignment_strategy,
                    elapsed_seconds=0.0,
                    newly_completed_loci_per_second=None,
                )
            )
        if (
            cached_alignment_cache is None
            or cached_alignment_cache.cache_version != SCHEME_ALIGNMENT_CACHE_VERSION
            or set(selected_method_cache.locus_entries)
            != set(source_signature_by_relative_path)
        ):
            method_caches[alignment_method_identifier] = _CachedAlignmentMethod(
                mafft_version=selected_method_cache.mafft_version,
                locus_entries=reusable_locus_entries,
            )
            _try_write_cached_alignment_cache(
                resolved_scheme_directory_path,
                method_caches,
            )
        return SchemeAlignmentAnalysis(
            mafft_version=selected_method_cache.mafft_version,
            alignment_strategy=resolved_alignment_strategy,
            alignment_method_identifier=alignment_method_identifier,
            locus_statistics=tuple(
                reusable_locus_entries[relative_path].locus_statistics
                for relative_path in source_signature_by_relative_path
            ),
        )

    resolved_maximum_parallel_loci = _resolve_maximum_parallel_loci(
        maximum_parallel_loci
    )
    mafft_version = _get_mafft_version(mafft_executable)
    if (
        selected_method_cache is None
        or selected_method_cache.mafft_version != mafft_version
    ):
        source_compatible_locus_entries = {}
        reusable_locus_entries = {}
    relative_paths_requiring_analysis = tuple(
        relative_path
        for relative_path in source_signature_by_relative_path
        if relative_path not in reusable_locus_entries
    )
    cached_locus_count = len(reusable_locus_entries)
    effective_worker_count = min(
        resolved_maximum_parallel_loci,
        max(1, len(relative_paths_requiring_analysis)),
    )
    if progress_callback is not None:
        progress_callback(
            AlignmentProgress(
                total_locus_count=total_locus_count,
                completed_locus_count=cached_locus_count,
                cached_locus_count=cached_locus_count,
                current_locus_name=None,
                alignment_strategy=resolved_alignment_strategy,
                elapsed_seconds=0.0,
                newly_completed_loci_per_second=None,
            )
        )
    progress_reporter = _MafftProgressReporter(
        total_locus_count=total_locus_count,
        initial_completed_locus_count=cached_locus_count,
        worker_count=effective_worker_count,
        alignment_strategy=resolved_alignment_strategy,
        output_stream=progress_output_stream,
    )
    completed_locus_count = cached_locus_count
    completed_since_checkpoint = 0
    last_checkpoint_time = time.monotonic()
    current_locus_entries = dict(reusable_locus_entries)
    progress_reporter.start()
    executor = ThreadPoolExecutor(max_workers=effective_worker_count)
    try:
        future_to_relative_path = {
            executor.submit(
                _analyze_locus,
                resolved_scheme_directory_path,
                allele_file_path_by_relative_path[relative_path],
                mafft_executable,
                resolved_alignment_strategy,
                source_compatible_locus_entries.get(relative_path),
            ): relative_path
            for relative_path in relative_paths_requiring_analysis
        }
        for completed_future in as_completed(future_to_relative_path):
            relative_path = future_to_relative_path[completed_future]
            locus_statistics = completed_future.result()
            file_size, modification_time_ns = source_signature_by_relative_path[
                relative_path
            ]
            current_locus_entries[relative_path] = _CachedLocusAlignment(
                file_size=file_size,
                modification_time_ns=modification_time_ns,
                locus_statistics=locus_statistics,
                nucleotide_polymorphism_metrics_calculated=True,
            )
            completed_locus_count += 1
            completed_since_checkpoint += 1
            progress_reporter.advance()
            if progress_callback is not None:
                elapsed_seconds = max(
                    0.0,
                    time.monotonic() - progress_reporter.started_at,
                )
                newly_completed_locus_count = completed_locus_count - cached_locus_count
                progress_callback(
                    AlignmentProgress(
                        total_locus_count=total_locus_count,
                        completed_locus_count=completed_locus_count,
                        cached_locus_count=cached_locus_count,
                        current_locus_name=locus_statistics.locus_name,
                        alignment_strategy=resolved_alignment_strategy,
                        elapsed_seconds=elapsed_seconds,
                        newly_completed_loci_per_second=(
                            newly_completed_locus_count / elapsed_seconds
                            if elapsed_seconds > 0
                            else None
                        ),
                    )
                )
            current_time = time.monotonic()
            if (
                completed_since_checkpoint
                >= SCHEME_ALIGNMENT_CACHE_CHECKPOINT_LOCUS_INTERVAL
                or current_time - last_checkpoint_time
                >= SCHEME_ALIGNMENT_CACHE_CHECKPOINT_SECONDS
            ):
                method_caches[alignment_method_identifier] = _CachedAlignmentMethod(
                    mafft_version=mafft_version,
                    locus_entries=dict(current_locus_entries),
                )
                _try_write_cached_alignment_cache(
                    resolved_scheme_directory_path,
                    method_caches,
                )
                completed_since_checkpoint = 0
                last_checkpoint_time = current_time
    finally:
        executor.shutdown(wait=True, cancel_futures=True)
        progress_reporter.finish()
        method_caches[alignment_method_identifier] = _CachedAlignmentMethod(
            mafft_version=mafft_version,
            locus_entries=dict(current_locus_entries),
        )
        _try_write_cached_alignment_cache(
            resolved_scheme_directory_path,
            method_caches,
        )

    alignment_analysis = SchemeAlignmentAnalysis(
        mafft_version=mafft_version,
        alignment_strategy=resolved_alignment_strategy,
        alignment_method_identifier=alignment_method_identifier,
        locus_statistics=tuple(
            current_locus_entries[relative_path].locus_statistics
            for relative_path in source_signature_by_relative_path
        ),
    )
    return alignment_analysis


def calculate_scheme_alignment_summary(
    alignment_analysis: SchemeAlignmentAnalysis,
    amino_acid_gap_warning_threshold_percentage: float = (
        DEFAULT_AMINO_ACID_GAP_WARNING_THRESHOLD_PERCENTAGE
    ),
) -> SchemeAlignmentSummary:
    """Summarize availability and apply an advisory amino-acid gap threshold.

    The comparison is strictly greater than the supplied threshold. The
    warning percentage uses successfully aligned amino-acid loci as its
    denominator. This threshold is an inspection aid, not a scheme pass/fail
    criterion.
    """
    if not 0 <= amino_acid_gap_warning_threshold_percentage <= 100:
        raise ValueError(
            "amino_acid_gap_warning_threshold_percentage must be between 0 and 100"
        )
    locus_statistics = alignment_analysis.locus_statistics
    nucleotide_aligned_locus_count = sum(
        locus.nucleotide_alignment_status == ALIGNMENT_STATUS_ALIGNED
        for locus in locus_statistics
    )
    nucleotide_unavailable_locus_count = sum(
        locus.nucleotide_alignment_status == ALIGNMENT_STATUS_INSUFFICIENT_SEQUENCES
        for locus in locus_statistics
    )
    nucleotide_failed_locus_count = sum(
        locus.nucleotide_alignment_status == ALIGNMENT_STATUS_FAILED
        for locus in locus_statistics
    )
    amino_acid_aligned_locus_count = sum(
        locus.amino_acid_alignment_status == ALIGNMENT_STATUS_ALIGNED
        for locus in locus_statistics
    )
    amino_acid_unavailable_locus_count = sum(
        locus.amino_acid_alignment_status == ALIGNMENT_STATUS_INSUFFICIENT_SEQUENCES
        for locus in locus_statistics
    )
    amino_acid_failed_locus_count = sum(
        locus.amino_acid_alignment_status == ALIGNMENT_STATUS_FAILED
        for locus in locus_statistics
    )
    amino_acid_gap_warning_loci = tuple(
        locus.locus_name
        for locus in locus_statistics
        if locus.amino_acid_gap_proportion_percentage is not None
        and locus.amino_acid_gap_proportion_percentage
        > amino_acid_gap_warning_threshold_percentage
    )
    amino_acid_gap_warning_locus_percentage = (
        len(amino_acid_gap_warning_loci) / amino_acid_aligned_locus_count * 100
        if amino_acid_aligned_locus_count
        else 0.0
    )
    return SchemeAlignmentSummary(
        locus_count=len(locus_statistics),
        nucleotide_aligned_locus_count=nucleotide_aligned_locus_count,
        nucleotide_unavailable_locus_count=nucleotide_unavailable_locus_count,
        nucleotide_failed_locus_count=nucleotide_failed_locus_count,
        amino_acid_aligned_locus_count=amino_acid_aligned_locus_count,
        amino_acid_unavailable_locus_count=(amino_acid_unavailable_locus_count),
        amino_acid_failed_locus_count=amino_acid_failed_locus_count,
        amino_acid_gap_warning_threshold_percentage=(
            amino_acid_gap_warning_threshold_percentage
        ),
        amino_acid_gap_warning_locus_count=len(amino_acid_gap_warning_loci),
        amino_acid_gap_warning_locus_percentage=(
            amino_acid_gap_warning_locus_percentage
        ),
        amino_acid_gap_warning_loci=amino_acid_gap_warning_loci,
    )


def _parse_positive_integer(argument_value: str) -> int:
    """Parse a positive integer command-line argument."""
    parsed_integer = int(argument_value)
    if parsed_integer < 1:
        raise argparse.ArgumentTypeError("value must be at least one")
    return parsed_integer


def _parse_percentage(argument_value: str) -> float:
    """Parse a percentage between zero and one hundred."""
    parsed_percentage = float(argument_value)
    if not 0 <= parsed_percentage <= 100:
        raise argparse.ArgumentTypeError("percentage must be between 0 and 100")
    return parsed_percentage


def main() -> None:
    """Print cached or calculated alignment statistics as JSON."""
    argument_parser = argparse.ArgumentParser(
        description=(
            "Align alleles within each locus and report nucleotide and "
            "amino-acid gap proportions."
        )
    )
    argument_parser.add_argument(
        "scheme_directory",
        type=Path,
        help="directory containing an extracted MLST scheme",
    )
    argument_parser.add_argument(
        "--mafft",
        default="mafft",
        dest="mafft_executable",
        help="MAFFT executable name or path (default: mafft)",
    )
    argument_parser.add_argument(
        "--strategy",
        choices=tuple(ALIGNMENT_STRATEGY_LABELS),
        default=DEFAULT_ALIGNMENT_STRATEGY,
        dest="alignment_strategy",
        help=(
            "MAFFT strategy: fft-ns-2 is the fast default; fft-ns-1 is "
            "fastest with a rougher guide tree; auto lets MAFFT choose"
        ),
    )
    argument_parser.add_argument(
        "--workers",
        type=_parse_positive_integer,
        default=None,
        dest="maximum_parallel_loci",
        metavar="COUNT",
        help=(
            "maximum concurrent single-threaded locus alignments (default: "
            "80%% of available CPUs, rounded down to an even count when "
            "possible)"
        ),
    )
    argument_parser.add_argument(
        "--amino-acid-gap-threshold",
        type=_parse_percentage,
        default=DEFAULT_AMINO_ACID_GAP_WARNING_THRESHOLD_PERCENTAGE,
        dest="amino_acid_gap_warning_threshold_percentage",
        metavar="PERCENT",
        help=("advisory amino-acid gap threshold (default: 5; comparison is >)"),
    )
    command_line_arguments = argument_parser.parse_args()
    alignment_analysis = get_scheme_alignment_analysis(
        command_line_arguments.scheme_directory,
        mafft_executable=command_line_arguments.mafft_executable,
        alignment_strategy=command_line_arguments.alignment_strategy,
        maximum_parallel_loci=command_line_arguments.maximum_parallel_loci,
        progress_output_stream=sys.stderr,
    )
    alignment_summary = calculate_scheme_alignment_summary(
        alignment_analysis,
        amino_acid_gap_warning_threshold_percentage=(
            command_line_arguments.amino_acid_gap_warning_threshold_percentage
        ),
    )
    print(
        json.dumps(
            {
                "mafft_version": alignment_analysis.mafft_version,
                "alignment_strategy": alignment_analysis.alignment_strategy,
                "alignment_method_identifier": (
                    alignment_analysis.alignment_method_identifier
                ),
                "summary": asdict(alignment_summary),
                "locus_statistics": [
                    asdict(locus_statistics)
                    for locus_statistics in alignment_analysis.locus_statistics
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
