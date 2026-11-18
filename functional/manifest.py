"""Load and validate versioned representative-genome set manifests.

Phase 2 analyses require more provenance than a collection of assembly files.
This module defines the small, typed contract used by the functional Marimo
notebook: dataset identity, the intrinsic scheme being evaluated, versioned
assembly identifiers, local input paths and checksums, quality metadata,
lineage context, and inclusion rationale.

Validation distinguishes blocking errors from non-blocking warnings. Unknown
optional provenance may be recorded explicitly, but missing analysis inputs,
placeholder values, duplicate identifiers, and an empty genome set prevent the
manifest from being considered ready for computation.
"""

import hashlib
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

GENOME_MANIFEST_SCHEMA_VERSION = 1
PLACEHOLDER_VALUE_PREFIX = "TODO"

__all__ = [
    "GENOME_MANIFEST_SCHEMA_VERSION",
    "GenomeManifestIssue",
    "GenomeManifestValidation",
    "GenomeRecord",
    "GenomeSetManifest",
    "load_genome_set_manifest",
    "validate_genome_set_manifest",
]

type IssueSeverity = Literal["error", "warning"]


@dataclass(frozen=True)
class GenomeRecord:
    """Provenance, local input, quality, and context for one assembly."""

    genome_id: str
    scientific_name: str
    taxonomy_identifier: str
    assembly_accession: str
    source_url: str
    retrieved_date: str
    local_path: str
    sha256: str
    assembly_level: str
    contig_count: int | None
    total_length_bp: int | None
    n50_bp: int | None
    completeness_percentage: float | None
    contamination_percentage: float | None
    lineage: str
    inclusion_rationale: str
    source_license: str
    source_citation: str
    notes: str


@dataclass(frozen=True)
class GenomeSetManifest:
    """Versioned genome-set metadata loaded from one TOML document."""

    schema_version: int
    dataset_id: str
    title: str
    description: str
    scheme_id: str
    created_date: str
    dataset_version: str
    provenance_status: str
    genomes: tuple[GenomeRecord, ...]
    manifest_path: Path
    manifest_sha256: str


@dataclass(frozen=True)
class GenomeManifestIssue:
    """One actionable manifest validation finding."""

    severity: IssueSeverity
    location: str
    message: str


@dataclass(frozen=True)
class GenomeManifestValidation:
    """A loaded manifest together with ordered readiness findings."""

    manifest: GenomeSetManifest
    issues: tuple[GenomeManifestIssue, ...]

    @property
    def error_count(self) -> int:
        """Return the number of computation-blocking findings."""
        return sum(issue.severity == "error" for issue in self.issues)

    @property
    def warning_count(self) -> int:
        """Return the number of non-blocking provenance findings."""
        return sum(issue.severity == "warning" for issue in self.issues)

    @property
    def is_ready(self) -> bool:
        """Return whether the manifest can identify real Phase 2 inputs."""
        return self.error_count == 0


def _text_value(raw_value: object) -> str:
    """Return a TOML scalar as stripped text, or an empty string."""
    return raw_value.strip() if isinstance(raw_value, str) else ""


def _optional_integer(raw_value: object) -> int | None:
    """Return an integer value while rejecting booleans and other types."""
    if isinstance(raw_value, bool) or not isinstance(raw_value, int):
        return None
    return raw_value


def _optional_float(raw_value: object) -> float | None:
    """Return a numeric value while rejecting booleans and other types."""
    if isinstance(raw_value, bool) or not isinstance(raw_value, int | float):
        return None
    return float(raw_value)


def _is_missing_or_placeholder(value: str) -> bool:
    """Return whether text is empty or intentionally marked for completion."""
    return not value or value.upper().startswith(PLACEHOLDER_VALUE_PREFIX)


def _load_genome_record(raw_record: object) -> GenomeRecord:
    """Convert one parsed TOML table into a typed genome record."""
    record = raw_record if isinstance(raw_record, dict) else {}
    return GenomeRecord(
        genome_id=_text_value(record.get("genome_id")),
        scientific_name=_text_value(record.get("scientific_name")),
        taxonomy_identifier=_text_value(record.get("taxonomy_identifier")),
        assembly_accession=_text_value(record.get("assembly_accession")),
        source_url=_text_value(record.get("source_url")),
        retrieved_date=_text_value(record.get("retrieved_date")),
        local_path=_text_value(record.get("local_path")),
        sha256=_text_value(record.get("sha256")),
        assembly_level=_text_value(record.get("assembly_level")),
        contig_count=_optional_integer(record.get("contig_count")),
        total_length_bp=_optional_integer(record.get("total_length_bp")),
        n50_bp=_optional_integer(record.get("n50_bp")),
        completeness_percentage=_optional_float(
            record.get("completeness_percentage")
        ),
        contamination_percentage=_optional_float(
            record.get("contamination_percentage")
        ),
        lineage=_text_value(record.get("lineage")),
        inclusion_rationale=_text_value(record.get("inclusion_rationale")),
        source_license=_text_value(record.get("source_license")),
        source_citation=_text_value(record.get("source_citation")),
        notes=_text_value(record.get("notes")),
    )


def load_genome_set_manifest(
    manifest_file_path: str | Path,
) -> GenomeSetManifest:
    """Load a Phase 2 genome-set manifest from TOML.

    Args:
        manifest_file_path: Path to a UTF-8 TOML manifest.

    Returns:
        The typed manifest, including its SHA-256 content identity.

    Raises:
        FileNotFoundError: If the manifest does not exist.
        IsADirectoryError: If the path does not identify a regular file.
        tomllib.TOMLDecodeError: If the document is not valid TOML.
    """
    resolved_manifest_file_path = Path(manifest_file_path).resolve()
    manifest_bytes = resolved_manifest_file_path.read_bytes()
    parsed_manifest = tomllib.loads(manifest_bytes.decode("utf-8"))
    raw_genomes = parsed_manifest.get("genomes", ())
    if not isinstance(raw_genomes, list):
        raw_genomes = ()
    raw_schema_version = parsed_manifest.get("schema_version")
    schema_version = (
        raw_schema_version
        if isinstance(raw_schema_version, int)
        and not isinstance(raw_schema_version, bool)
        else 0
    )
    return GenomeSetManifest(
        schema_version=schema_version,
        dataset_id=_text_value(parsed_manifest.get("dataset_id")),
        title=_text_value(parsed_manifest.get("title")),
        description=_text_value(parsed_manifest.get("description")),
        scheme_id=_text_value(parsed_manifest.get("scheme_id")),
        created_date=_text_value(parsed_manifest.get("created_date")),
        dataset_version=_text_value(parsed_manifest.get("dataset_version")),
        provenance_status=_text_value(parsed_manifest.get("provenance_status")),
        genomes=tuple(_load_genome_record(record) for record in raw_genomes),
        manifest_path=resolved_manifest_file_path,
        manifest_sha256=hashlib.sha256(manifest_bytes).hexdigest(),
    )


def _require_text(
    issues: list[GenomeManifestIssue],
    *,
    location: str,
    value: str,
) -> None:
    """Add a blocking issue when required text is absent or a placeholder."""
    if _is_missing_or_placeholder(value):
        issues.append(
            GenomeManifestIssue(
                severity="error",
                location=location,
                message="Provide a verified value; TODO and empty values are not ready.",
            )
        )


def _warn_for_missing_text(
    issues: list[GenomeManifestIssue],
    *,
    location: str,
    value: str,
    message: str,
) -> None:
    """Add a warning for optional but reusable provenance."""
    if _is_missing_or_placeholder(value):
        issues.append(
            GenomeManifestIssue(
                severity="warning",
                location=location,
                message=message,
            )
        )


def validate_genome_set_manifest(
    manifest: GenomeSetManifest,
) -> GenomeManifestValidation:
    """Return deterministic readiness errors and provenance warnings."""
    issues: list[GenomeManifestIssue] = []
    if manifest.schema_version != GENOME_MANIFEST_SCHEMA_VERSION:
        issues.append(
            GenomeManifestIssue(
                severity="error",
                location="schema_version",
                message=(
                    f"Expected schema version {GENOME_MANIFEST_SCHEMA_VERSION}; "
                    f"found {manifest.schema_version}."
                ),
            )
        )
    for field_name, field_value in (
        ("dataset_id", manifest.dataset_id),
        ("title", manifest.title),
        ("description", manifest.description),
        ("scheme_id", manifest.scheme_id),
        ("created_date", manifest.created_date),
        ("dataset_version", manifest.dataset_version),
        ("provenance_status", manifest.provenance_status),
    ):
        _require_text(issues, location=field_name, value=field_value)
    if (
        not _is_missing_or_placeholder(manifest.provenance_status)
        and manifest.provenance_status.lower() != "complete"
    ):
        issues.append(
            GenomeManifestIssue(
                severity="warning",
                location="provenance_status",
                message="Mark provenance complete only after verifying every record.",
            )
        )
    if not manifest.genomes:
        issues.append(
            GenomeManifestIssue(
                severity="error",
                location="genomes",
                message="Add at least one representative genome record.",
            )
        )

    seen_genome_ids: set[str] = set()
    for genome_index, genome in enumerate(manifest.genomes):
        location_prefix = f"genomes[{genome_index}]"
        for field_name, field_value in (
            ("genome_id", genome.genome_id),
            ("scientific_name", genome.scientific_name),
            ("assembly_accession", genome.assembly_accession),
            ("source_url", genome.source_url),
            ("retrieved_date", genome.retrieved_date),
            ("local_path", genome.local_path),
            ("sha256", genome.sha256),
            ("inclusion_rationale", genome.inclusion_rationale),
        ):
            _require_text(
                issues,
                location=f"{location_prefix}.{field_name}",
                value=field_value,
            )
        if genome.genome_id in seen_genome_ids and genome.genome_id:
            issues.append(
                GenomeManifestIssue(
                    severity="error",
                    location=f"{location_prefix}.genome_id",
                    message=f"Duplicate genome identifier: {genome.genome_id}",
                )
            )
        seen_genome_ids.add(genome.genome_id)
        if (
            genome.sha256
            and not _is_missing_or_placeholder(genome.sha256)
            and (
                len(genome.sha256) != 64
                or any(
                    character not in "0123456789abcdefABCDEF"
                    for character in genome.sha256
                )
            )
        ):
            issues.append(
                GenomeManifestIssue(
                    severity="error",
                    location=f"{location_prefix}.sha256",
                    message="Use a 64-character hexadecimal SHA-256 digest.",
                )
            )
        if genome.local_path and not _is_missing_or_placeholder(
            genome.local_path
        ):
            local_path = Path(genome.local_path).expanduser()
            if not local_path.is_absolute():
                local_path = manifest.manifest_path.parent / local_path
            if not local_path.is_file():
                issues.append(
                    GenomeManifestIssue(
                        severity="error",
                        location=f"{location_prefix}.local_path",
                        message=(
                            "The manifest-relative assembly path is not an "
                            "accessible file."
                        ),
                    )
                )
        for field_name, field_value, message in (
            (
                "taxonomy_identifier",
                genome.taxonomy_identifier,
                "Record the NCBI Taxonomy identifier when available.",
            ),
            (
                "assembly_level",
                genome.assembly_level,
                "Record the assembly level or an explicit unknown status.",
            ),
            (
                "lineage",
                genome.lineage,
                "Record lineage context or explain why it is unavailable.",
            ),
            (
                "source_license",
                genome.source_license,
                "Record source terms or an explicit unknown status.",
            ),
            (
                "source_citation",
                genome.source_citation,
                "Record the recommended citation or an explicit unknown status.",
            ),
        ):
            _warn_for_missing_text(
                issues,
                location=f"{location_prefix}.{field_name}",
                value=field_value,
                message=message,
            )
        for field_name, field_value in (
            ("contig_count", genome.contig_count),
            ("total_length_bp", genome.total_length_bp),
            ("n50_bp", genome.n50_bp),
            ("completeness_percentage", genome.completeness_percentage),
            ("contamination_percentage", genome.contamination_percentage),
        ):
            if field_value is None:
                issues.append(
                    GenomeManifestIssue(
                        severity="warning",
                        location=f"{location_prefix}.{field_name}",
                        message=(
                            "Record the assembly-quality value or explain its "
                            "absence in notes."
                        ),
                    )
                )
        for field_name, field_value in (
            ("contig_count", genome.contig_count),
            ("total_length_bp", genome.total_length_bp),
            ("n50_bp", genome.n50_bp),
        ):
            if field_value is not None and field_value < 0:
                issues.append(
                    GenomeManifestIssue(
                        severity="error",
                        location=f"{location_prefix}.{field_name}",
                        message="Use a non-negative integer.",
                    )
                )
        for field_name, field_value in (
            ("completeness_percentage", genome.completeness_percentage),
            ("contamination_percentage", genome.contamination_percentage),
        ):
            if field_value is not None and not 0 <= field_value <= 100:
                issues.append(
                    GenomeManifestIssue(
                        severity="error",
                        location=f"{location_prefix}.{field_name}",
                        message="Use a percentage between 0 and 100.",
                    )
                )
    return GenomeManifestValidation(manifest=manifest, issues=tuple(issues))
