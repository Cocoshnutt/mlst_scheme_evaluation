"""Tests for Phase 2 representative-genome manifest readiness."""

import json
import tempfile
import unittest
from pathlib import Path

from functional.manifest import (
    GENOME_MANIFEST_SCHEMA_VERSION,
    load_genome_set_manifest,
    validate_genome_set_manifest,
)

REPOSITORY_ROOT_DIRECTORY = Path(__file__).resolve().parents[1]


def _valid_manifest_text(*, genome_id: str = "genome-1") -> str:
    """Return a complete one-genome manifest for a temporary fixture."""
    return f'''schema_version = 1
dataset_id = "yersinia-validation-v1"
title = "Yersinia representative genomes"
description = "A temporary validation fixture."
scheme_id = "cgmlst-yersinia-genus"
created_date = "2026-11-18"
dataset_version = "1"
provenance_status = "complete"

[[genomes]]
genome_id = "{genome_id}"
scientific_name = "Yersinia enterocolitica"
taxonomy_identifier = "NCBITaxon:630"
assembly_accession = "GCF_000000001.1"
source_url = "https://example.org/assembly"
retrieved_date = "2026-11-18"
local_path = "genome.fna"
sha256 = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
assembly_level = "Complete Genome"
contig_count = 1
total_length_bp = 4500000
n50_bp = 4500000
completeness_percentage = 99.5
contamination_percentage = 0.2
lineage = "lineage-1"
inclusion_rationale = "Representative temporary test genome"
source_license = "test license"
source_citation = "test citation"
notes = ""
'''


class GenomeManifestTest(unittest.TestCase):
    """Protect typed loading, readiness rules, and the FAIR schema artifact."""

    def test_complete_manifest_is_ready(self) -> None:
        """Accept verified fields and an accessible manifest-relative input."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory_path = Path(temporary_directory)
            (directory_path / "genome.fna").write_text(
                ">contig\nACGT\n",
                encoding="utf-8",
            )
            manifest_path = directory_path / "genome-set.toml"
            manifest_path.write_text(_valid_manifest_text(), encoding="utf-8")

            manifest = load_genome_set_manifest(manifest_path)
            validation = validate_genome_set_manifest(manifest)

        self.assertEqual(manifest.schema_version, GENOME_MANIFEST_SCHEMA_VERSION)
        self.assertEqual(manifest.genomes[0].genome_id, "genome-1")
        self.assertEqual(len(manifest.manifest_sha256), 64)
        self.assertTrue(validation.is_ready)
        self.assertEqual(validation.error_count, 0)
        self.assertEqual(validation.warning_count, 0)

    def test_repository_template_is_deliberately_not_ready(self) -> None:
        """Keep missing real genomes visible instead of fabricating inputs."""
        template_path = (
            REPOSITORY_ROOT_DIRECTORY
            / "metadata"
            / "genome_sets"
            / "template.toml"
        )

        validation = validate_genome_set_manifest(
            load_genome_set_manifest(template_path)
        )

        self.assertFalse(validation.is_ready)
        issue_locations = {issue.location for issue in validation.issues}
        self.assertIn("dataset_id", issue_locations)
        self.assertIn("scheme_id", issue_locations)
        self.assertIn("genomes", issue_locations)

    def test_invalid_genome_fields_are_actionable(self) -> None:
        """Reject duplicate IDs, bad checksums, bounds, and missing inputs."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory_path = Path(temporary_directory)
            second_genome_table = "\n[[genomes]]" + _valid_manifest_text(
                genome_id="genome-1"
            ).split("[[genomes]]", maxsplit=1)[1]
            manifest_text = _valid_manifest_text() + second_genome_table.replace(
                'sha256 = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"',
                'sha256 = "not-a-checksum"',
            ).replace(
                "completeness_percentage = 99.5",
                "completeness_percentage = 101.0",
            )
            manifest_path = directory_path / "invalid.toml"
            manifest_path.write_text(manifest_text, encoding="utf-8")

            validation = validate_genome_set_manifest(
                load_genome_set_manifest(manifest_path)
            )

        issue_locations = {issue.location for issue in validation.issues}
        self.assertIn("genomes[1].genome_id", issue_locations)
        self.assertIn("genomes[1].sha256", issue_locations)
        self.assertIn("genomes[1].completeness_percentage", issue_locations)
        self.assertIn("genomes[0].local_path", issue_locations)
        self.assertFalse(validation.is_ready)

    def test_json_schema_is_machine_readable_and_versioned(self) -> None:
        """Keep the interoperability contract parseable without dependencies."""
        schema_path = (
            REPOSITORY_ROOT_DIRECTORY
            / "metadata"
            / "genome_manifest.schema.json"
        )
        schema = json.loads(schema_path.read_text(encoding="utf-8"))

        self.assertEqual(schema["properties"]["schema_version"]["const"], 1)
        self.assertIn("genomes", schema["required"])
        self.assertFalse(schema["additionalProperties"])


if __name__ == "__main__":
    unittest.main()
