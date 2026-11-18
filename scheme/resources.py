"""Resolve shared compute-resource budgets for expensive scheme analyses.

The default preserves capacity for the operating system and other processes by
using 80% of the logical CPUs available to the current process.  The result is
rounded down to an even number when possible so paired or concurrent workloads
divide cleanly, with a one-worker floor for small or constrained machines.

Callers may provide an explicit positive worker count.  Explicit values are
not capped because command-line users may intentionally coordinate resources
outside this process, while the Marimo control limits its choices to the CPUs
currently visible to the process.
"""

import os

DEFAULT_CPU_UTILIZATION_PERCENTAGE = 80

__all__ = [
    "DEFAULT_CPU_UTILIZATION_PERCENTAGE",
    "calculate_default_worker_count",
    "get_available_cpu_count",
    "resolve_worker_count",
]


def get_available_cpu_count() -> int:
    """Return the logical CPU count available to the current process."""
    return max(1, os.process_cpu_count() or 1)


def calculate_default_worker_count(
    available_cpu_count: int,
    *,
    utilization_percentage: int = DEFAULT_CPU_UTILIZATION_PERCENTAGE,
) -> int:
    """Return the conservative default worker count for an available CPU count.

    The percentage is applied with integer floor division.  Counts of two or
    more are then rounded down to an even value; a one-worker floor prevents a
    zero budget on small machines.
    """
    if available_cpu_count < 1:
        raise ValueError("available_cpu_count must be at least 1")
    if not 1 <= utilization_percentage <= 100:
        raise ValueError("utilization_percentage must be between 1 and 100")

    worker_count = max(
        1,
        available_cpu_count * utilization_percentage // 100,
    )
    if worker_count >= 2:
        worker_count -= worker_count % 2
    return max(1, worker_count)


def resolve_worker_count(maximum_workers: int | None) -> int:
    """Return an explicit positive budget or the conservative default."""
    if maximum_workers is not None:
        if maximum_workers < 1:
            raise ValueError("maximum_workers must be at least 1")
        return maximum_workers
    return calculate_default_worker_count(get_available_cpu_count())
