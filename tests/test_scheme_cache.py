"""Tests for atomic scheme-archive extraction and cache clearing."""

import io
import tarfile
import tempfile
import unittest
from pathlib import Path

from scheme.cache import (
    COMPLETE_CACHE_MARKER_FILENAME,
    cache_archive,
    clear_cache,
)


def _write_tar_archive(
    archive_file_path: Path,
    members: dict[str, bytes],
) -> None:
    """Create a small gzip-compressed tar archive for one test."""
    with tarfile.open(archive_file_path, mode="w:gz") as archive_file:
        for member_name, member_content in members.items():
            member_information = tarfile.TarInfo(member_name)
            member_information.size = len(member_content)
            archive_file.addfile(member_information, io.BytesIO(member_content))


class SchemeArchiveCacheTest(unittest.TestCase):
    """Protect extraction identity, atomicity, and reuse semantics."""

    def test_complete_cached_extraction_is_reused(self) -> None:
        """Do not overwrite a completed extraction on a repeated request."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_directory_path = Path(temporary_directory)
            archive_file_path = temporary_directory_path / "scheme.tar.gz"
            cache_root_directory_path = temporary_directory_path / "cache"
            _write_tar_archive(
                archive_file_path,
                {"locus.fasta": b">allele_1\nATGTAA\n"},
            )

            first_cached_directory_path = cache_archive(
                archive_file_path,
                cache_root_directory_path,
            )
            extracted_locus_file_path = (
                first_cached_directory_path / "locus.fasta"
            )
            extracted_locus_file_path.write_text(
                ">local_marker\nGTGTAG\n",
                encoding="utf-8",
            )

            second_cached_directory_path = cache_archive(
                archive_file_path,
                cache_root_directory_path,
            )

            self.assertEqual(
                second_cached_directory_path,
                first_cached_directory_path,
            )
            self.assertEqual(
                extracted_locus_file_path.read_text(encoding="utf-8"),
                ">local_marker\nGTGTAG\n",
            )
            self.assertTrue(
                (
                    first_cached_directory_path
                    / COMPLETE_CACHE_MARKER_FILENAME
                ).is_file()
            )

    def test_changed_archive_metadata_creates_a_new_cache_identity(self) -> None:
        """Use a distinct extraction after the source archive changes."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_directory_path = Path(temporary_directory)
            archive_file_path = temporary_directory_path / "scheme.tar.gz"
            cache_root_directory_path = temporary_directory_path / "cache"
            _write_tar_archive(
                archive_file_path,
                {"locus.fasta": b">allele_1\nATGTAA\n"},
            )
            first_cached_directory_path = cache_archive(
                archive_file_path,
                cache_root_directory_path,
            )

            _write_tar_archive(
                archive_file_path,
                {"locus.fasta": b">allele_1\nATGCCCTAA\n"},
            )
            second_cached_directory_path = cache_archive(
                archive_file_path,
                cache_root_directory_path,
            )

            self.assertNotEqual(
                second_cached_directory_path,
                first_cached_directory_path,
            )
            self.assertTrue(first_cached_directory_path.is_dir())
            self.assertEqual(
                (second_cached_directory_path / "locus.fasta").read_text(
                    encoding="utf-8"
                ),
                ">allele_1\nATGCCCTAA\n",
            )

    def test_incomplete_cached_extraction_is_replaced(self) -> None:
        """Never reuse a destination whose completion marker is missing."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_directory_path = Path(temporary_directory)
            archive_file_path = temporary_directory_path / "scheme.tar.gz"
            cache_root_directory_path = temporary_directory_path / "cache"
            _write_tar_archive(
                archive_file_path,
                {"locus.fasta": b">allele_1\nATGTAA\n"},
            )
            cached_directory_path = cache_archive(
                archive_file_path,
                cache_root_directory_path,
            )
            (
                cached_directory_path / COMPLETE_CACHE_MARKER_FILENAME
            ).unlink()
            (cached_directory_path / "stale.txt").write_text(
                "partial extraction",
                encoding="utf-8",
            )

            refreshed_directory_path = cache_archive(
                archive_file_path,
                cache_root_directory_path,
            )

            self.assertEqual(refreshed_directory_path, cached_directory_path)
            self.assertFalse((refreshed_directory_path / "stale.txt").exists())
            self.assertTrue(
                (
                    refreshed_directory_path
                    / COMPLETE_CACHE_MARKER_FILENAME
                ).is_file()
            )

    def test_failed_extraction_removes_temporary_directory(self) -> None:
        """Do not leave a partial extraction after an invalid archive."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_directory_path = Path(temporary_directory)
            archive_file_path = temporary_directory_path / "invalid.tar.gz"
            cache_root_directory_path = temporary_directory_path / "cache"
            archive_file_path.write_bytes(b"not a tar archive")

            with self.assertRaises(tarfile.ReadError):
                cache_archive(archive_file_path, cache_root_directory_path)

            self.assertEqual(tuple(cache_root_directory_path.iterdir()), ())

    def test_archive_member_cannot_escape_the_cache_directory(self) -> None:
        """Retain tarfile's data-filter path traversal protection."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_directory_path = Path(temporary_directory)
            archive_file_path = temporary_directory_path / "unsafe.tar.gz"
            cache_root_directory_path = temporary_directory_path / "cache"
            escaped_file_path = temporary_directory_path / "escaped.txt"
            _write_tar_archive(
                archive_file_path,
                {"../../escaped.txt": b"must not be extracted"},
            )

            with self.assertRaises(tarfile.FilterError):
                cache_archive(archive_file_path, cache_root_directory_path)

            self.assertFalse(escaped_file_path.exists())
            self.assertEqual(tuple(cache_root_directory_path.iterdir()), ())


class ClearSchemeCacheTest(unittest.TestCase):
    """Protect clear-cache scope and safety checks."""

    def test_clear_cache_removes_only_the_requested_cache_root(self) -> None:
        """Remove a custom cache tree while retaining its parent directory."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_directory_path = Path(temporary_directory)
            cache_root_directory_path = temporary_directory_path / "cache"
            cache_root_directory_path.mkdir()
            (cache_root_directory_path / "result.json").write_text(
                "{}",
                encoding="utf-8",
            )

            clear_cache(cache_root_directory_path)

            self.assertFalse(cache_root_directory_path.exists())
            self.assertTrue(temporary_directory_path.is_dir())

    def test_clear_cache_refuses_the_current_working_directory(self) -> None:
        """Do not permit the broadest accidental project deletion target."""
        with self.assertRaises(ValueError):
            clear_cache(Path.cwd())


if __name__ == "__main__":
    unittest.main()
