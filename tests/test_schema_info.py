"""Tests for scheme-wide nucleotide-composition statistics."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from schema import info
from schema.info import (
    SCHEMA_INFORMATION_CACHE_FILENAME,
    SCHEMA_INFORMATION_CACHE_VERSION,
    count_nucleotide_base_composition,
    get_schema_information,
)


class NucleotideBaseCompositionTest(unittest.TestCase):
    """Test the canonical-base GC-content definition."""

    def test_whitespace_case_and_noncanonical_characters_are_handled(self) -> None:
        """Ignore whitespace and exclude ambiguous characters from the denominator."""
        nucleotide_base_composition = count_nucleotide_base_composition(
            "aC gT\nN-"
        )

        self.assertEqual(
            nucleotide_base_composition.guanine_cytosine_base_count,
            2,
        )
        self.assertEqual(
            nucleotide_base_composition.canonical_nucleotide_base_count,
            4,
        )
        self.assertEqual(
            nucleotide_base_composition.excluded_nucleotide_base_count,
            2,
        )
        self.assertEqual(
            nucleotide_base_composition.guanine_cytosine_content_percentage,
            50.0,
        )

    def test_no_canonical_bases_has_unavailable_percentage(self) -> None:
        """Do not represent an undefined GC denominator as zero percent."""
        nucleotide_base_composition = count_nucleotide_base_composition("NN--")

        self.assertIsNone(
            nucleotide_base_composition.guanine_cytosine_content_percentage
        )
        self.assertEqual(
            nucleotide_base_composition.excluded_nucleotide_base_count,
            4,
        )


class SchemaGuanineCytosineContentTest(unittest.TestCase):
    """Test pooled scheme-wide GC content and cache reuse."""

    def test_scheme_content_pools_canonical_bases_across_loci(self) -> None:
        """Weight the scheme value by bases instead of averaging locus values."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            (schema_directory_path / "first.fasta").write_text(
                ">allele_1\nACGTNN\n>allele_2\nGG\n",
                encoding="utf-8",
            )
            (schema_directory_path / "second.fasta").write_text(
                ">allele_1\nAAAA\n",
                encoding="utf-8",
            )

            schema_information = get_schema_information(schema_directory_path)

        self.assertEqual(schema_information.canonical_nucleotide_base_count, 10)
        self.assertEqual(schema_information.excluded_nucleotide_base_count, 2)
        self.assertEqual(schema_information.guanine_cytosine_content_percentage, 40.0)

    def test_cached_composition_is_reused_without_rescanning_fasta(self) -> None:
        """Persist raw base counts so later summary requests stay inexpensive."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            schema_directory_path = Path(temporary_directory)
            (schema_directory_path / "cached.fasta").write_text(
                ">allele_1\nACGTN\n",
                encoding="utf-8",
            )
            get_schema_information(schema_directory_path)

            with patch.object(
                info,
                "_count_alleles_in_files",
                side_effect=AssertionError("FASTA should not be rescanned"),
            ):
                cached_schema_information = get_schema_information(
                    schema_directory_path,
                    minimum_sequence_length=100,
                )

            cached_data = json.loads(
                (
                    schema_directory_path / SCHEMA_INFORMATION_CACHE_FILENAME
                ).read_text(encoding="utf-8")
            )

        self.assertEqual(cached_schema_information.canonical_nucleotide_base_count, 4)
        self.assertEqual(cached_schema_information.excluded_nucleotide_base_count, 1)
        self.assertEqual(cached_data["cache_version"], SCHEMA_INFORMATION_CACHE_VERSION)
        self.assertEqual(cached_data["guanine_cytosine_base_count"], 2)


if __name__ == "__main__":
    unittest.main()
