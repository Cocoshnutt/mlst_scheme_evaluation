"""Tests for shared compute-resource budgeting."""

import unittest
from unittest.mock import patch

from scheme.resources import (
    calculate_default_worker_count,
    get_available_cpu_count,
    resolve_worker_count,
)


class ResourceBudgetTest(unittest.TestCase):
    """Protect conservative defaults and explicit caller overrides."""

    def test_default_worker_count_is_low_core_safe_and_even_when_possible(
        self,
    ) -> None:
        """Apply the 80% policy without ever returning zero workers."""
        expected_worker_counts = {
            1: 1,
            2: 1,
            3: 2,
            4: 2,
            5: 4,
            6: 4,
            7: 4,
            8: 6,
            10: 8,
            12: 8,
            16: 12,
        }

        for available_cpu_count, expected_worker_count in (
            expected_worker_counts.items()
        ):
            with self.subTest(available_cpu_count=available_cpu_count):
                self.assertEqual(
                    calculate_default_worker_count(available_cpu_count),
                    expected_worker_count,
                )

    def test_process_available_cpu_count_has_a_one_cpu_fallback(self) -> None:
        """Handle runtimes that cannot report a process-visible CPU count."""
        with patch("scheme.resources.os.process_cpu_count", return_value=None):
            self.assertEqual(get_available_cpu_count(), 1)

    def test_implicit_budget_uses_process_available_cpu_count(self) -> None:
        """Respect affinity and cgroup limits exposed by Python."""
        with patch("scheme.resources.os.process_cpu_count", return_value=12):
            self.assertEqual(resolve_worker_count(None), 8)

    def test_explicit_positive_budget_is_preserved(self) -> None:
        """Leave deliberate library and CLI overrides authoritative."""
        self.assertEqual(resolve_worker_count(7), 7)

    def test_invalid_counts_and_percentages_are_rejected(self) -> None:
        """Reject configurations that could create an empty executor."""
        with self.assertRaisesRegex(ValueError, "available_cpu_count"):
            calculate_default_worker_count(0)
        with self.assertRaisesRegex(ValueError, "utilization_percentage"):
            calculate_default_worker_count(4, utilization_percentage=0)
        with self.assertRaisesRegex(ValueError, "maximum_workers"):
            resolve_worker_count(0)


if __name__ == "__main__":
    unittest.main()
