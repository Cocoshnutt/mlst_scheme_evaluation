"""Tests for shared atomic JSON cache persistence."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scheme.cache_io import write_json_atomically


class AtomicJsonWriteTest(unittest.TestCase):
    """Protect replacement, formatting, propagation, and cleanup behavior."""

    def test_compact_json_atomically_replaces_an_existing_cache(self) -> None:
        """Write the established compact representation to the target path."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_file_path = Path(temporary_directory) / ".analysis-cache.json"
            output_file_path.write_text('{"old":true}', encoding="utf-8")

            write_json_atomically(
                output_file_path,
                {"cache_version": 3, "values": [1, 2]},
            )

            self.assertEqual(
                output_file_path.read_text(encoding="utf-8"),
                '{"cache_version":3,"values":[1,2]}',
            )
            self.assertEqual(
                json.loads(output_file_path.read_text(encoding="utf-8")),
                {"cache_version": 3, "values": [1, 2]},
            )
            self.assertEqual(
                tuple(Path(temporary_directory).iterdir()), (output_file_path,)
            )

    def test_replacement_failure_preserves_target_and_removes_temporary_file(
        self,
    ) -> None:
        """Do not leave partial JSON beside a previously valid cache."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_file_path = Path(temporary_directory) / ".analysis-cache.json"
            output_file_path.write_text('{"old":true}', encoding="utf-8")

            with (
                patch.object(
                    Path,
                    "replace",
                    side_effect=PermissionError("read-only cache"),
                ),
                self.assertRaises(PermissionError),
            ):
                write_json_atomically(output_file_path, {"new": True})

            self.assertEqual(
                output_file_path.read_text(encoding="utf-8"),
                '{"old":true}',
            )
            self.assertEqual(
                tuple(Path(temporary_directory).iterdir()), (output_file_path,)
            )

    def test_serialization_failure_removes_temporary_file(self) -> None:
        """Clean up when a future cache payload is not JSON serializable."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_file_path = Path(temporary_directory) / ".analysis-cache.json"

            with self.assertRaises(TypeError):
                write_json_atomically(output_file_path, {"invalid": object()})

            self.assertFalse(output_file_path.exists())
            self.assertEqual(tuple(Path(temporary_directory).iterdir()), ())

    def test_temporary_file_creation_errors_propagate_to_the_owner(self) -> None:
        """Let each analysis retain its established cache-failure policy."""
        with (
            patch(
                "scheme.cache_io.tempfile.NamedTemporaryFile",
                side_effect=PermissionError("read-only cache"),
            ),
            self.assertRaises(PermissionError),
        ):
            write_json_atomically(Path("cache.json"), {"cache_version": 1})


if __name__ == "__main__":
    unittest.main()
