"""Tests for the generated Conda environment specification."""

import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts import sync_conda_environment
from scripts.sync_conda_environment import (
    CondaConfiguration,
    CondaEnvironmentGenerationError,
    export_locked_python_requirements,
    load_conda_configuration,
    render_conda_environment,
    write_if_changed,
)


class CondaEnvironmentRenderingTest(unittest.TestCase):
    """Test project metadata parsing and deterministic YAML rendering."""

    def test_loads_project_and_conda_metadata(self) -> None:
        """Use the project name, Python constraint, channels, and native tools."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            pyproject_path = Path(temporary_directory) / "pyproject.toml"
            pyproject_path.write_text(
                """
[project]
name = "example-project"
requires-python = ">=3.14,<3.15"
dependencies = ["example>=1"]

[tool.mlst_scheme_evaluation.conda]
channels = ["conda-forge", "bioconda", "nodefaults"]
dependencies = ["conda-forge::mafft", "bioconda::blast"]
""".lstrip(),
                encoding="utf-8",
            )

            configuration = load_conda_configuration(pyproject_path)

        self.assertEqual(configuration.environment_name, "example-project")
        self.assertEqual(configuration.python_requirement, ">=3.14,<3.15")
        self.assertEqual(
            configuration.conda_dependencies,
            ("conda-forge::mafft", "bioconda::blast"),
        )

    def test_renders_native_tools_and_locked_pip_requirements(self) -> None:
        """Keep Conda tools separate from exact, marker-aware Python pins."""
        content = render_conda_environment(
            CondaConfiguration(
                environment_name="example-project",
                python_requirement=">=3.14,<3.15",
                channels=("conda-forge", "bioconda", "nodefaults"),
                conda_dependencies=("conda-forge::mafft", "bioconda::blast"),
            ),
            (
                "altair==6.3.0",
                "click==8.5.1 ; sys_platform != 'emscripten'",
            ),
        )

        self.assertIn('name: "example-project"', content)
        self.assertIn('  - "python >=3.14,<3.15"', content)
        self.assertIn('  - "conda-forge::mafft"', content)
        self.assertIn('  - "bioconda::blast"', content)
        self.assertIn('      - "altair==6.3.0"', content)
        self.assertIn(
            '      - "click==8.5.1 ; sys_platform != \'emscripten\'"',
            content,
        )
        self.assertTrue(content.endswith("\n"))


class CondaEnvironmentSynchronizationTest(unittest.TestCase):
    """Test no-op, replacement, and uv-export failure behavior."""

    def test_identical_content_is_not_rewritten(self) -> None:
        """Preserve inode and timestamps when the generated file is current."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_path = Path(temporary_directory) / "environment.yml"
            output_path.write_text("current\n", encoding="utf-8")
            original_stat = output_path.stat()

            changed = write_if_changed(output_path, "current\n")
            current_stat = output_path.stat()

        self.assertFalse(changed)
        self.assertEqual(current_stat.st_ino, original_stat.st_ino)
        self.assertEqual(current_stat.st_mtime_ns, original_stat.st_mtime_ns)

    def test_changed_content_is_atomically_replaced(self) -> None:
        """Replace stale output and leave no temporary sibling behind."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_path = Path(temporary_directory) / "environment.yml"
            output_path.write_text("stale\n", encoding="utf-8")

            changed = write_if_changed(output_path, "current\n")

            self.assertTrue(changed)
            self.assertEqual(output_path.read_text(encoding="utf-8"), "current\n")
            self.assertEqual(tuple(Path(temporary_directory).glob("*.tmp")), ())

    def test_uv_export_failure_is_reported(self) -> None:
        """Reject a stale or otherwise unexportable uv lock file."""
        failed_process = subprocess.CompletedProcess(
            args=["uv", "export"],
            returncode=1,
            stdout="",
            stderr="The lock file needs to be updated",
        )
        with (
            patch(
                "scripts.sync_conda_environment.subprocess.run",
                return_value=failed_process,
            ),
            self.assertRaisesRegex(
                CondaEnvironmentGenerationError,
                "lock file needs to be updated",
            ),
        ):
            export_locked_python_requirements(Path.cwd())

    def test_check_mode_does_not_modify_stale_output(self) -> None:
        """Return a failure without rewriting in freshness-check mode."""
        with tempfile.TemporaryDirectory() as temporary_directory:
            project_root = Path(temporary_directory)
            output_path = project_root / "environment.yml"
            output_path.write_text("stale\n", encoding="utf-8")
            with (
                patch.object(sync_conda_environment, "PROJECT_ROOT", project_root),
                patch.object(
                    sync_conda_environment,
                    "generate_conda_environment",
                    return_value="current\n",
                ),
            ):
                exit_code = sync_conda_environment.main(["--check"])

            self.assertEqual(exit_code, 1)
            self.assertEqual(output_path.read_text(encoding="utf-8"), "stale\n")


if __name__ == "__main__":
    unittest.main()
