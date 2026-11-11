"""Extract and cache MLST scheme archives for repeated use.

This module avoids extracting an unchanged scheme archive more than once. Each
archive is extracted into a cache directory identified by its source path, file
size, and modification time. Extraction first happens in a temporary directory
so an interrupted operation cannot be mistaken for a complete cache entry.

Notebook code can import :func:`cache_archive` and :func:`clear_cache`. The
module can also be executed as a command-line program to populate or clear the
same cache.
"""

import argparse
import hashlib
import shutil
import tarfile
import tempfile
from pathlib import Path

DEFAULT_CACHE_ROOT_DIRECTORY = Path(".cache/schemes")
COMPLETE_CACHE_MARKER_FILENAME = ".complete"
ARCHIVE_FILENAME_SUFFIXES = (
    ".tar.gz",
    ".tar.bz2",
    ".tar.xz",
    ".tar.zst",
    ".tgz",
    ".tbz2",
    ".txz",
    ".tzst",
    ".tar",
)

__all__ = [
    "DEFAULT_CACHE_ROOT_DIRECTORY",
    "cache_archive",
    "clear_cache",
]


def _create_archive_fingerprint(archive_file_path: Path) -> str:
    """Create a short cache identifier from an archive's path and metadata.

    Args:
        archive_file_path: Absolute path to the archive being identified.

    Returns:
        The first twelve hexadecimal characters of the metadata hash.
    """
    archive_file_status = archive_file_path.stat()
    archive_identity = (
        f"{archive_file_path}:"
        f"{archive_file_status.st_size}:"
        f"{archive_file_status.st_mtime_ns}"
    )
    return hashlib.sha256(archive_identity.encode()).hexdigest()[:12]


def _remove_archive_filename_suffix(archive_filename: str) -> str:
    """Remove a recognized tar archive suffix from a filename.

    Args:
        archive_filename: Filename whose archive suffix should be removed.

    Returns:
        The filename without its recognized archive suffix. An unrecognized
        filename is returned unchanged.
    """
    for archive_filename_suffix in ARCHIVE_FILENAME_SUFFIXES:
        if archive_filename.endswith(archive_filename_suffix):
            return archive_filename.removesuffix(archive_filename_suffix)
    return archive_filename


def cache_archive(
    archive_file_path: str | Path,
    cache_root_directory_path: str | Path = DEFAULT_CACHE_ROOT_DIRECTORY,
) -> Path:
    """Extract an archive once and return its reusable cache directory.

    An unchanged archive returns the existing directory without extracting it
    again. Supported compression formats are detected by :mod:`tarfile`.

    Args:
        archive_file_path: Path to the compressed tar archive to extract.
        cache_root_directory_path: Directory in which extracted archives are
            cached.

    Returns:
        The absolute path to the archive's extracted cache directory.

    Raises:
        FileNotFoundError: If ``archive_file_path`` is not an existing file.
        tarfile.TarError: If the file is not a readable tar archive.
    """
    resolved_archive_file_path = Path(archive_file_path).resolve()
    resolved_cache_root_directory_path = Path(cache_root_directory_path).resolve()

    if not resolved_archive_file_path.is_file():
        raise FileNotFoundError(
            f"Archive does not exist: {resolved_archive_file_path}"
        )

    resolved_cache_root_directory_path.mkdir(parents=True, exist_ok=True)

    archive_fingerprint = _create_archive_fingerprint(resolved_archive_file_path)
    archive_base_name = _remove_archive_filename_suffix(
        resolved_archive_file_path.name
    )
    cached_archive_directory_path = resolved_cache_root_directory_path / (
        f"{archive_base_name}-{archive_fingerprint}"
    )
    complete_cache_marker_path = (
        cached_archive_directory_path / COMPLETE_CACHE_MARKER_FILENAME
    )

    if complete_cache_marker_path.exists():
        return cached_archive_directory_path

    temporary_extraction_directory_path = Path(
        tempfile.mkdtemp(
            prefix=".extract-",
            dir=resolved_cache_root_directory_path,
        )
    )

    try:
        with tarfile.open(resolved_archive_file_path, "r:*") as archive_file:
            archive_file.extractall(
                temporary_extraction_directory_path,
                filter="data",
            )

        temporary_complete_cache_marker_path = (
            temporary_extraction_directory_path / COMPLETE_CACHE_MARKER_FILENAME
        )
        temporary_complete_cache_marker_path.touch()

        if cached_archive_directory_path.exists():
            shutil.rmtree(cached_archive_directory_path)
        temporary_extraction_directory_path.rename(cached_archive_directory_path)

    except Exception:
        shutil.rmtree(temporary_extraction_directory_path, ignore_errors=True)
        raise

    return cached_archive_directory_path


def clear_cache(
    cache_root_directory_path: str | Path = DEFAULT_CACHE_ROOT_DIRECTORY,
) -> None:
    """Remove every extracted scheme from the cache directory.

    Args:
        cache_root_directory_path: Cache directory to remove. The default is
            ``.cache/schemes`` relative to the current working directory.

    Raises:
        ValueError: If the resolved cache path is a filesystem root or the
            current working directory.
    """
    resolved_cache_root_directory_path = Path(cache_root_directory_path).resolve()
    filesystem_root_directory_path = Path(
        resolved_cache_root_directory_path.anchor
    )
    current_working_directory_path = Path.cwd().resolve()

    if resolved_cache_root_directory_path in (
        filesystem_root_directory_path,
        current_working_directory_path,
    ):
        raise ValueError(
            "Refusing to clear unsafe cache path: "
            f"{resolved_cache_root_directory_path}"
        )

    shutil.rmtree(resolved_cache_root_directory_path, ignore_errors=True)


def main() -> None:
    """Parse command-line arguments and populate or clear the scheme cache."""
    argument_parser = argparse.ArgumentParser(
        description="Extract an archive into a reusable local cache."
    )
    argument_parser.add_argument(
        "archive_file",
        nargs="?",
        type=Path,
        help="tar archive to cache (compression is detected automatically)",
    )
    argument_parser.add_argument(
        "--cache-root",
        dest="cache_root_directory",
        type=Path,
        default=DEFAULT_CACHE_ROOT_DIRECTORY,
        help="cache directory (default: .cache/schemes)",
    )
    argument_parser.add_argument(
        "--clear-cache",
        action="store_true",
        help="remove the cache instead of extracting an archive",
    )
    command_line_arguments = argument_parser.parse_args()

    if command_line_arguments.clear_cache:
        if command_line_arguments.archive_file is not None:
            argument_parser.error("archive cannot be used with --clear-cache")
        clear_cache(command_line_arguments.cache_root_directory)
        print(
            "Cleared cache: "
            f"{command_line_arguments.cache_root_directory.resolve()}"
        )
        return

    if command_line_arguments.archive_file is None:
        argument_parser.error("archive is required unless --clear-cache is used")

    cached_archive_directory_path = cache_archive(
        command_line_arguments.archive_file,
        command_line_arguments.cache_root_directory,
    )
    print(cached_archive_directory_path)


if __name__ == "__main__":
    main()
