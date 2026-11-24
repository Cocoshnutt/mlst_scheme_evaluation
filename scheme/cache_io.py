"""Shared low-level persistence for versioned JSON analysis caches.

Cache ownership, schemas, compatibility rules, validation, and deserialization
remain in their analysis modules. This module only provides the common
same-directory temporary-file and atomic-replacement write operation.
"""

import json
import tempfile
from pathlib import Path

__all__ = ["write_json_atomically"]


def write_json_atomically(
    output_file_path: Path,
    cache_data: object,
) -> None:
    """Write compact JSON and atomically replace ``output_file_path``.

    The temporary file is created beside the destination so replacement stays
    on the same filesystem. Any temporary file is removed when serialization
    or replacement fails. Operational errors propagate so each owning
    analysis can retain its established fatal or non-fatal cache policy.

    Args:
        output_file_path: Final JSON cache path.
        cache_data: JSON-serializable cache payload.

    Raises:
        OSError: If the temporary file cannot be written or replaced.
        TypeError: If ``cache_data`` is not JSON serializable.
        ValueError: If JSON serialization rejects a value.
    """
    temporary_file_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=output_file_path.parent,
            prefix=f"{output_file_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary_file:
            temporary_file_path = Path(temporary_file.name)
            json.dump(cache_data, temporary_file, separators=(",", ":"))
        temporary_file_path.replace(output_file_path)
    finally:
        if temporary_file_path is not None:
            temporary_file_path.unlink(missing_ok=True)
