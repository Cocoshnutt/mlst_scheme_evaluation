"""Create FAIR-oriented intrinsic evaluation bundles from typed results.

Internal analysis caches are deliberately not exported. The bundle contains
open JSON and TSV representations of the currently available typed results,
input and software provenance, explicit stage status, and SHA-256 checksums.
Creating a bundle never invokes MAFFT, BLAST, or a scientific calculation.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import platform
import tomllib
import zipfile
from dataclasses import asdict, dataclass, fields
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import TYPE_CHECKING

from .alignment import LocusAlignmentStatistics
from .allele_info import LocusAlleleInformation
from .bsr import (
    HighBsrLocusGroup,
    HighBsrLocusPair,
    InterLocusBlastScoreRatio,
    LocusRepresentativeBsrCoverage,
    LocusRepresentativeBsrMeasurements,
)
from .similarity import (
    CandidateParalogousLocusPair,
    CandidateParalogThresholds,
    CrossLocusSimilarityMatch,
)
from .static_analysis import SCHEME_STATIC_ANALYSIS_METHOD_IDENTIFIER
from .stats import LocusStatistics

if TYPE_CHECKING:
    from collections.abc import Iterable

    from .alignment import SchemeAlignmentAnalysis, SchemeAlignmentSummary
    from .bsr import (
        RepresentativeBsrClassification,
        SchemeRepresentativeBsrAnalysis,
    )
    from .overview import IntrinsicSchemeOverview
    from .similarity import CandidateParalogAnalysis, SchemeSimilarityAnalysis

INTRINSIC_EVALUATION_BUNDLE_VERSION = 1
UNAVAILABLE_TSV_VALUE = "unavailable"

__all__ = [
    "INTRINSIC_EVALUATION_BUNDLE_VERSION",
    "IntrinsicOperationalSettings",
    "create_intrinsic_evaluation_bundle",
    "intrinsic_evaluation_bundle_filename",
    "load_scheme_source_record",
]


@dataclass(frozen=True)
class IntrinsicOperationalSettings:
    """Submitted notebook resource, method, and classification settings.

    These values describe how the evaluation was requested. They are exported
    for reproducibility but remain separate from scientific cache identities.
    """

    detected_available_cpu_count: int
    submitted_worker_budget: int
    selected_mafft_strategy: str
    amino_acid_gap_warning_threshold_percentage: float
    candidate_paralog_thresholds: CandidateParalogThresholds
    representative_bsr_threshold: float


def _json_bytes(value: object) -> bytes:
    """Serialize one human-readable JSON artifact with a terminal newline."""
    return (
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    ).encode("utf-8")


def _tsv_value(value: object) -> object:
    """Convert nested and unavailable values to explicit TSV cells."""
    if value is None:
        return UNAVAILABLE_TSV_VALUE
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, tuple | list | dict):
        return json.dumps(value, separators=(",", ":"), ensure_ascii=False)
    return value


def _dataclass_rows_to_tsv(
    row_type: type[object],
    rows: Iterable[object],
) -> bytes:
    """Serialize immutable result rows as a rectangular UTF-8 TSV table."""
    field_names = tuple(field.name for field in fields(row_type))
    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output,
        fieldnames=field_names,
        dialect="excel-tab",
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    for row in rows:
        row_data = asdict(row)
        writer.writerow(
            {
                field_name: _tsv_value(row_data[field_name])
                for field_name in field_names
            }
        )
    return output.getvalue().encode("utf-8")


def _sha256_file(file_path: Path) -> str:
    """Return a streaming SHA-256 digest for one input file."""
    digest = hashlib.sha256()
    with file_path.open("rb") as input_file:
        for file_chunk in iter(lambda: input_file.read(1024 * 1024), b""):
            digest.update(file_chunk)
    return digest.hexdigest()


def _contains_placeholder(value: object) -> bool:
    """Return whether metadata explicitly contains an unfinished TODO value."""
    return isinstance(value, str) and "TODO" in value.upper()


def _load_scheme_source_inventory(
    inventory_file_path: Path,
    archive_file_path: Path,
    repository_root_directory_path: Path,
) -> tuple[dict[str, object], dict[str, object] | None, tuple[str, ...]]:
    """Load inventory metadata, one source record, and its metadata gaps."""
    inventory_data = tomllib.loads(
        inventory_file_path.read_text(encoding="utf-8")
    )
    inventory_metadata = {
        str(key): value
        for key, value in inventory_data.items()
        if key != "schemes"
    }
    relative_archive_path = archive_file_path.resolve().relative_to(
        repository_root_directory_path.resolve()
    ).as_posix()
    for raw_record in inventory_data.get("schemes", []):
        if not isinstance(raw_record, dict):
            continue
        if raw_record.get("archive_path") != relative_archive_path:
            continue
        source_record = {str(key): value for key, value in raw_record.items()}
        placeholder_fields = tuple(
            sorted(
                field_name
                for field_name, field_value in source_record.items()
                if _contains_placeholder(field_value)
            )
        )
        return inventory_metadata, source_record, placeholder_fields
    return inventory_metadata, None, ("record",)


def load_scheme_source_record(
    inventory_file_path: Path,
    archive_file_path: Path,
    repository_root_directory_path: Path,
) -> tuple[dict[str, object] | None, tuple[str, ...]]:
    """Load the matching scheme-source record and report placeholder fields.

    The function never guesses metadata. A missing record returns ``None`` and
    one explicit ``record`` gap; unfinished values are returned unchanged and
    identified by field name.
    """
    _, source_record, placeholder_fields = _load_scheme_source_inventory(
        inventory_file_path,
        archive_file_path,
        repository_root_directory_path,
    )
    return source_record, placeholder_fields


def intrinsic_evaluation_bundle_filename(archive_file_path: Path) -> str:
    """Return a stable download filename for one scheme archive."""
    archive_name = archive_file_path.name
    for suffix in (".tar.zst", ".tar.gz", ".tar.xz", ".tgz", ".zip"):
        if archive_name.lower().endswith(suffix):
            archive_name = archive_name[: -len(suffix)]
            break
    safe_name = "".join(
        character if character.isalnum() or character in "-_" else "-"
        for character in archive_name
    ).strip("-")
    return f"{safe_name or 'scheme'}-intrinsic-evaluation.zip"


def _create_summary(
    overview: IntrinsicSchemeOverview,
) -> dict[str, object]:
    """Create a compact scheme-wide summary without duplicating detail rows."""
    static_summary = asdict(overview.static)
    alignment_summary = None
    if overview.alignment is not None:
        alignment_summary = {
            **asdict(overview.alignment.summary),
            "alignment_method_identifier": (
                overview.alignment.analysis.alignment_method_identifier
            ),
            "alignment_strategy": (
                overview.alignment.analysis.alignment_strategy
            ),
            "mafft_version": overview.alignment.analysis.mafft_version,
            "nucleotide_callable_site_count": (
                overview.alignment.nucleotide_callable_site_count
            ),
            "nucleotide_polymorphic_site_count": (
                overview.alignment.nucleotide_polymorphic_site_count
            ),
            "nucleotide_polymorphic_site_percentage": (
                overview.alignment.nucleotide_polymorphic_site_percentage
            ),
            "polymorphic_locus_count": overview.alignment.polymorphic_locus_count,
        }
    similarity_summary = None
    if overview.similarity is not None:
        analysis = overview.similarity.analysis
        classification = overview.similarity.candidate_paralog_analysis
        similarity_summary = {
            "method_identifier": analysis.method_identifier,
            "blastn_version": analysis.blastn_version,
            "blastp_version": analysis.blastp_version,
            "blast_evalue": analysis.blast_evalue,
            "reused_from_cache": analysis.reused_from_cache,
            "locus_count": analysis.locus_count,
            "nucleotide_selected_locus_count": analysis.nucleotide_selected_locus_count,
            "amino_acid_selected_locus_count": analysis.amino_acid_selected_locus_count,
            "nucleotide_possible_directed_pair_count": (
                analysis.nucleotide_possible_directed_pair_count
            ),
            "amino_acid_possible_directed_pair_count": (
                analysis.amino_acid_possible_directed_pair_count
            ),
            "nucleotide_match_count": analysis.nucleotide_match_count,
            "amino_acid_match_count": analysis.amino_acid_match_count,
            "candidate_classification": (
                None
                if classification is None
                else {
                    "thresholds": asdict(classification.thresholds),
                    "evaluated_locus_pair_count": (
                        classification.evaluated_locus_pair_count
                    ),
                    "candidate_pair_count": classification.candidate_pair_count,
                    "candidate_locus_count": classification.candidate_locus_count,
                    "amino_acid_candidate_pair_count": (
                        classification.amino_acid_candidate_pair_count
                    ),
                    "nucleotide_candidate_pair_count": (
                        classification.nucleotide_candidate_pair_count
                    ),
                    "both_criteria_candidate_pair_count": (
                        classification.both_criteria_candidate_pair_count
                    ),
                    "blast_score_ratio_evidence_pair_count": (
                        classification.blast_score_ratio_evidence_pair_count
                    ),
                }
            ),
        }
    representative_bsr_summary = None
    if overview.representative_bsr is not None:
        analysis = overview.representative_bsr.analysis
        classification = overview.representative_bsr.classification
        representative_bsr_summary = {
            "method_identifier": analysis.method_identifier,
            "blastp_version": analysis.blastp_version,
            "blast_evalue": analysis.blast_evalue,
            "reused_from_cache": analysis.reused_from_cache,
            "locus_count": analysis.locus_count,
            "representative_locus_count": analysis.representative_locus_count,
            "self_score_locus_count": analysis.self_score_locus_count,
            "total_allele_count": analysis.total_allele_count,
            "valid_translated_allele_count": analysis.valid_translated_allele_count,
            "classification": (
                None
                if classification is None
                else {
                    "blast_score_ratio_threshold": (
                        classification.blast_score_ratio_threshold
                    ),
                    "measurable_allele_count": classification.measurable_allele_count,
                    "covered_allele_count": classification.covered_allele_count,
                    "covered_allele_percentage": (
                        classification.covered_allele_percentage
                    ),
                    "unavailable_valid_allele_count": (
                        classification.unavailable_valid_allele_count
                    ),
                    "incomplete_coverage_locus_count": (
                        classification.incomplete_coverage_locus_count
                    ),
                    "high_bsr_pair_count": len(classification.high_bsr_pairs),
                    "high_bsr_locus_count": classification.high_bsr_locus_count,
                    "high_bsr_group_count": len(classification.high_bsr_groups),
                }
            ),
        }
    return {
        "static": static_summary,
        "alignment": alignment_summary,
        "similarity": similarity_summary,
        "representative_bsr": representative_bsr_summary,
    }


def _project_version() -> str:
    """Return the installed project version or an explicit unavailable value."""
    try:
        return version("mlst-scheme-evaluation")
    except PackageNotFoundError:
        return "unavailable"


def _alignment_cache_reuse(
    analysis: SchemeAlignmentAnalysis | None,
) -> dict[str, object] | None:
    """Describe per-locus MAFFT reuse without collapsing mixed requests."""
    if analysis is None:
        return None
    locus_count = len(analysis.locus_statistics)
    if locus_count == 0:
        reuse_status = "not_applicable"
    elif analysis.reused_locus_count == locus_count:
        reuse_status = "all_reused"
    elif analysis.reused_locus_count:
        reuse_status = "partially_reused"
    else:
        reuse_status = "not_reused"
    return {
        "status": reuse_status,
        "reused_locus_count": analysis.reused_locus_count,
        "calculated_locus_count": analysis.calculated_locus_count,
        "total_locus_count": locus_count,
    }


def _whole_analysis_cache_reuse(
    reused_from_cache: bool | None,
) -> str | None:
    """Describe whole-analysis BLAST cache reuse when a stage is available."""
    if reused_from_cache is None:
        return None
    return "reused" if reused_from_cache else "not_reused"


def _create_operational_provenance(
    *,
    settings: IntrinsicOperationalSettings,
    overview: IntrinsicSchemeOverview,
    scheme_id: object,
    archive_path: str,
    archive_sha256: str,
    inventory_path: str,
    inventory_metadata: dict[str, object],
    alignment_analysis: SchemeAlignmentAnalysis | None,
    similarity_analysis: SchemeSimilarityAnalysis | None,
    representative_bsr_analysis: SchemeRepresentativeBsrAnalysis | None,
) -> dict[str, object]:
    """Assemble submitted settings, input identities, and run provenance."""
    return {
        "submitted_settings": {
            "detected_available_cpu_count": (
                settings.detected_available_cpu_count
            ),
            "worker_budget": settings.submitted_worker_budget,
            "mafft_strategy": settings.selected_mafft_strategy,
            "thresholds": {
                "minimum_sequence_length_bp": (
                    overview.static.minimum_sequence_length
                ),
                "frequent_peak_relative_threshold_percentage": (
                    overview.static.frequent_peak_relative_threshold_percentage
                ),
                "amino_acid_gap_warning_threshold_percentage": (
                    settings.amino_acid_gap_warning_threshold_percentage
                ),
                "candidate_paralog": asdict(
                    settings.candidate_paralog_thresholds
                ),
                "representative_bsr": settings.representative_bsr_threshold,
            },
        },
        "input_identities": {
            "scheme_id": scheme_id,
            "scheme_archive_path": archive_path,
            "scheme_archive_sha256": archive_sha256,
            "scheme_source_inventory_path": inventory_path,
            "scheme_source_inventory_schema_version": inventory_metadata.get(
                "schema_version"
            ),
            "scheme_source_inventory_provenance_status": (
                inventory_metadata.get("provenance_status")
            ),
            "genome_set_manifest": None,
        },
        "analyses": {
            "static": {
                "method_identifier": SCHEME_STATIC_ANALYSIS_METHOD_IDENTIFIER,
                "tool_versions": {},
                "cache_reuse": "not_recorded_by_current_result_types",
            },
            "alignment": {
                "method_identifier": (
                    alignment_analysis.alignment_method_identifier
                    if alignment_analysis is not None
                    else None
                ),
                "mafft_strategy_used": (
                    alignment_analysis.alignment_strategy
                    if alignment_analysis is not None
                    else None
                ),
                "tool_versions": {
                    "mafft": (
                        alignment_analysis.mafft_version
                        if alignment_analysis is not None
                        else None
                    )
                },
                "parameters": {
                    "mafft_strategy": (
                        alignment_analysis.alignment_strategy
                        if alignment_analysis is not None
                        else settings.selected_mafft_strategy
                    )
                },
                "cache_reuse": _alignment_cache_reuse(alignment_analysis),
            },
            "cross_locus_similarity": {
                "method_identifier": (
                    similarity_analysis.method_identifier
                    if similarity_analysis is not None
                    else None
                ),
                "tool_versions": {
                    "blastn": (
                        similarity_analysis.blastn_version
                        if similarity_analysis is not None
                        else None
                    ),
                    "blastp": (
                        similarity_analysis.blastp_version
                        if similarity_analysis is not None
                        else None
                    ),
                },
                "parameters": {
                    "blast_evalue": (
                        similarity_analysis.blast_evalue
                        if similarity_analysis is not None
                        else None
                    )
                },
                "cache_reuse": _whole_analysis_cache_reuse(
                    similarity_analysis.reused_from_cache
                    if similarity_analysis is not None
                    else None
                ),
            },
            "representative_bsr": {
                "method_identifier": (
                    representative_bsr_analysis.method_identifier
                    if representative_bsr_analysis is not None
                    else None
                ),
                "tool_versions": {
                    "blastp": (
                        representative_bsr_analysis.blastp_version
                        if representative_bsr_analysis is not None
                        else None
                    )
                },
                "parameters": {
                    "blast_evalue": (
                        representative_bsr_analysis.blast_evalue
                        if representative_bsr_analysis is not None
                        else None
                    )
                },
                "cache_reuse": _whole_analysis_cache_reuse(
                    representative_bsr_analysis.reused_from_cache
                    if representative_bsr_analysis is not None
                    else None
                ),
            },
        },
        "cache_key_policy": {
            "worker_budget_is_scientific_cache_key": False,
            "note": (
                "Worker count affects scheduling only and is excluded from "
                "scientific cache identities."
            ),
        },
    }


def create_intrinsic_evaluation_bundle(
    *,
    repository_root_directory_path: Path,
    archive_file_path: Path,
    scheme_source_inventory_file_path: Path,
    operational_settings: IntrinsicOperationalSettings,
    locus_statistics: tuple[LocusStatistics, ...],
    scheme_allele_information: tuple[LocusAlleleInformation, ...],
    intrinsic_scheme_overview: IntrinsicSchemeOverview,
    scheme_alignment_analysis: SchemeAlignmentAnalysis | None,
    scheme_alignment_summary: SchemeAlignmentSummary | None,
    scheme_similarity_analysis: SchemeSimilarityAnalysis | None,
    candidate_paralog_analysis: CandidateParalogAnalysis | None,
    scheme_representative_bsr_analysis: SchemeRepresentativeBsrAnalysis | None,
    representative_bsr_classification: RepresentativeBsrClassification | None,
    created_at: datetime | None = None,
) -> bytes:
    """Return a ZIP bundle for all currently available intrinsic results.

    Optional analysis arguments may be ``None``. Their manifest status is then
    ``not_run`` and no corresponding detail table is emitted. The function
    serializes supplied results only and never performs an analysis.
    """
    archive_path = archive_file_path.resolve()
    repository_root = repository_root_directory_path.resolve()
    inventory_metadata, source_record, placeholder_fields = (
        _load_scheme_source_inventory(
            scheme_source_inventory_file_path,
            archive_path,
            repository_root,
        )
    )
    created_timestamp = created_at or datetime.now(UTC)
    if created_timestamp.tzinfo is None:
        raise ValueError("created_at must include a timezone")
    stage_status = {
        "static": "complete",
        "alignment": (
            "complete"
            if scheme_alignment_analysis is not None
            and scheme_alignment_summary is not None
            else "not_run"
        ),
        "similarity": (
            "complete" if scheme_similarity_analysis is not None else "not_run"
        ),
        "candidate_paralog_classification": (
            "complete" if candidate_paralog_analysis is not None else "not_run"
        ),
        "representative_bsr": (
            "complete"
            if scheme_representative_bsr_analysis is not None
            else "not_run"
        ),
        "representative_bsr_classification": (
            "complete"
            if representative_bsr_classification is not None
            else "not_run"
        ),
    }
    relative_archive_path = archive_path.relative_to(repository_root).as_posix()
    relative_inventory_path = (
        scheme_source_inventory_file_path.resolve()
        .relative_to(repository_root)
        .as_posix()
    )
    archive_sha256 = _sha256_file(archive_path)
    recorded_archive_sha256 = (
        source_record.get("archive_sha256")
        if source_record is not None
        else None
    )
    recorded_archive_sha256_matches = (
        archive_sha256 == recorded_archive_sha256
        if isinstance(recorded_archive_sha256, str)
        and not _contains_placeholder(recorded_archive_sha256)
        else None
    )
    artifacts: dict[str, bytes] = {
        "input_provenance.json": _json_bytes(
            {
                "archive_path": relative_archive_path,
                "archive_size_bytes": archive_path.stat().st_size,
                "archive_sha256": archive_sha256,
                "recorded_archive_sha256_matches": (
                    recorded_archive_sha256_matches
                ),
                "scheme_source_inventory": inventory_metadata,
                "scheme_source_record": source_record,
                "scheme_source_placeholder_fields": placeholder_fields,
                "provenance_complete": (
                    not placeholder_fields
                    and inventory_metadata.get("provenance_status")
                    == "complete"
                ),
            }
        ),
        "summary.json": _json_bytes(_create_summary(intrinsic_scheme_overview)),
        "locus_statistics.tsv": _dataclass_rows_to_tsv(
            LocusStatistics,
            locus_statistics,
        ),
        "locus_allele_information.tsv": _dataclass_rows_to_tsv(
            LocusAlleleInformation,
            scheme_allele_information,
        ),
    }
    if scheme_alignment_analysis is not None:
        artifacts["alignment_loci.tsv"] = _dataclass_rows_to_tsv(
            LocusAlignmentStatistics,
            scheme_alignment_analysis.locus_statistics,
        )
    if scheme_similarity_analysis is not None:
        artifacts["cross_locus_similarity.tsv"] = _dataclass_rows_to_tsv(
            CrossLocusSimilarityMatch,
            scheme_similarity_analysis.matches,
        )
    if candidate_paralog_analysis is not None:
        artifacts["candidate_paralog_pairs.tsv"] = _dataclass_rows_to_tsv(
            CandidateParalogousLocusPair,
            candidate_paralog_analysis.candidate_pairs,
        )
    if scheme_representative_bsr_analysis is not None:
        artifacts["representative_bsr_loci.tsv"] = _dataclass_rows_to_tsv(
            LocusRepresentativeBsrMeasurements,
            scheme_representative_bsr_analysis.locus_measurements,
        )
        artifacts["representative_bsr_inter_locus.tsv"] = _dataclass_rows_to_tsv(
            InterLocusBlastScoreRatio,
            scheme_representative_bsr_analysis.inter_locus_measurements,
        )
    if representative_bsr_classification is not None:
        artifacts["representative_bsr_coverage.tsv"] = _dataclass_rows_to_tsv(
            LocusRepresentativeBsrCoverage,
            representative_bsr_classification.locus_coverages,
        )
        artifacts["representative_bsr_pairs.tsv"] = _dataclass_rows_to_tsv(
            HighBsrLocusPair,
            representative_bsr_classification.high_bsr_pairs,
        )
        artifacts["representative_bsr_groups.tsv"] = _dataclass_rows_to_tsv(
            HighBsrLocusGroup,
            representative_bsr_classification.high_bsr_groups,
        )
    software_versions = {
        "mlst_scheme_evaluation": _project_version(),
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "mafft": (
            scheme_alignment_analysis.mafft_version
            if scheme_alignment_analysis is not None
            else None
        ),
        "blastn": (
            scheme_similarity_analysis.blastn_version
            if scheme_similarity_analysis is not None
            else None
        ),
        "blastp_similarity": (
            scheme_similarity_analysis.blastp_version
            if scheme_similarity_analysis is not None
            else None
        ),
        "blastp_representative_bsr": (
            scheme_representative_bsr_analysis.blastp_version
            if scheme_representative_bsr_analysis is not None
            else None
        ),
    }
    artifacts["software_versions.json"] = _json_bytes(software_versions)

    artifact_checksums = {
        file_name: hashlib.sha256(file_content).hexdigest()
        for file_name, file_content in sorted(artifacts.items())
    }
    manifest = {
        "bundle_format": "mlst-scheme-intrinsic-evaluation",
        "bundle_format_version": INTRINSIC_EVALUATION_BUNDLE_VERSION,
        "created_at": created_timestamp.astimezone(UTC).isoformat(),
        "scheme_id": (
            source_record.get("scheme_id") if source_record is not None else None
        ),
        "stage_status": stage_status,
        "missing_value_representation": {
            "json": "null",
            "tsv": UNAVAILABLE_TSV_VALUE,
        },
        "software_license": "EUPL-1.2",
        "operational_provenance": _create_operational_provenance(
            settings=operational_settings,
            overview=intrinsic_scheme_overview,
            scheme_id=(
                source_record.get("scheme_id")
                if source_record is not None
                else None
            ),
            archive_path=relative_archive_path,
            archive_sha256=archive_sha256,
            inventory_path=relative_inventory_path,
            inventory_metadata=inventory_metadata,
            alignment_analysis=scheme_alignment_analysis,
            similarity_analysis=scheme_similarity_analysis,
            representative_bsr_analysis=scheme_representative_bsr_analysis,
        ),
        "files": [
            {
                "path": file_name,
                "sha256": artifact_checksums[file_name],
            }
            for file_name in sorted(artifacts)
        ],
        "notes": (
            "Screening classifications are descriptive and must not be "
            "interpreted as an overall scheme pass/fail decision."
        ),
    }
    artifacts["manifest.json"] = _json_bytes(manifest)
    checksums = {
        file_name: hashlib.sha256(file_content).hexdigest()
        for file_name, file_content in sorted(artifacts.items())
    }
    artifacts["checksums.sha256"] = "".join(
        f"{file_checksum}  {file_name}\n"
        for file_name, file_checksum in checksums.items()
    ).encode("utf-8")

    bundle_buffer = io.BytesIO()
    with zipfile.ZipFile(
        bundle_buffer,
        mode="w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=6,
    ) as bundle_zip:
        for file_name, file_content in sorted(artifacts.items()):
            bundle_zip.writestr(file_name, file_content)
    return bundle_buffer.getvalue()
