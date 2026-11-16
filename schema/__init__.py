"""Reusable schema-evaluation analyses and command-line modules.

The package keeps scientific calculations separate from the Marimo notebook.
Individual modules expose typed library APIs and, where applicable, can be
invoked with ``python -m schema.<module>``.
"""

__all__ = [
    "alignment",
    "allele_info",
    "cache",
    "info",
    "similarity",
    "stats",
]
