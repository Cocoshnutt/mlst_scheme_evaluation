"""Generate a Conda environment specification from project and uv metadata."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import tomllib
from dataclasses import dataclass
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PYPROJECT_FILENAME = "pyproject.toml"
CONDA_ENVIRONMENT_FILENAME = "environment.yml"


@dataclass(frozen=True)
class CondaConfiguration:
    """Project fields required to render the Conda environment specification."""

    environment_name: str
    python_requirement: str
    channels: tuple[str, ...]
    conda_dependencies: tuple[str, ...]


class CondaEnvironmentGenerationError(RuntimeError):
    """Raised when the Conda environment specification cannot be generated."""


def _require_non_empty_string(value: object, field_name: str) -> str:
    """Return a configured string or raise a concise metadata error."""
    if not isinstance(value, str) or not value.strip():
        raise CondaEnvironmentGenerationError(
            f"{field_name} must be a non-empty string"
        )
    return value.strip()


def _require_string_tuple(value: object, field_name: str) -> tuple[str, ...]:
    """Return a non-empty tuple of configured strings."""
    if not isinstance(value, list) or not value:
        raise CondaEnvironmentGenerationError(
            f"{field_name} must be a non-empty list of strings"
        )
    values = tuple(
        _require_non_empty_string(item, f"{field_name} item") for item in value
    )
    if len(set(values)) != len(values):
        raise CondaEnvironmentGenerationError(
            f"{field_name} must not contain duplicate values"
        )
    return values


def load_conda_configuration(pyproject_path: Path) -> CondaConfiguration:
    """Load project and Conda-only dependencies from ``pyproject.toml``."""
    try:
        pyproject_data = tomllib.loads(pyproject_path.read_text(encoding="utf-8"))
        project_data = pyproject_data["project"]
        conda_data = pyproject_data["tool"]["mlst_scheme_evaluation"]["conda"]
    except FileNotFoundError as error:
        raise CondaEnvironmentGenerationError(
            f"Project metadata was not found: {pyproject_path}"
        ) from error
    except tomllib.TOMLDecodeError as error:
        raise CondaEnvironmentGenerationError(
            f"Project metadata is not valid TOML: {error}"
        ) from error
    except (KeyError, TypeError) as error:
        raise CondaEnvironmentGenerationError(
            "pyproject.toml must define [project] and "
            "[tool.mlst_scheme_evaluation.conda] metadata"
        ) from error

    return CondaConfiguration(
        environment_name=_require_non_empty_string(
            project_data.get("name"), "project.name"
        ),
        python_requirement=_require_non_empty_string(
            project_data.get("requires-python"), "project.requires-python"
        ),
        channels=_require_string_tuple(conda_data.get("channels"), "conda.channels"),
        conda_dependencies=_require_string_tuple(
            conda_data.get("dependencies"), "conda.dependencies"
        ),
    )


def export_locked_python_requirements(
    project_root: Path,
    *,
    uv_executable: str = "uv",
) -> tuple[str, ...]:
    """Export exact Python requirements from the current uv lock file."""
    command = (
        uv_executable,
        "export",
        "--locked",
        "--format",
        "requirements.txt",
        "--no-dev",
        "--no-emit-project",
        "--no-hashes",
        "--no-annotate",
        "--no-header",
    )
    try:
        completed_process = subprocess.run(
            command,
            cwd=project_root,
            check=False,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as error:
        raise CondaEnvironmentGenerationError(
            f"uv executable was not found: {uv_executable}"
        ) from error

    if completed_process.returncode != 0:
        error_message = completed_process.stderr.strip() or completed_process.stdout.strip()
        raise CondaEnvironmentGenerationError(
            f"uv could not export the locked Python dependencies: {error_message}"
        )

    requirements = tuple(
        line.strip()
        for line in completed_process.stdout.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    )
    if not requirements:
        raise CondaEnvironmentGenerationError(
            "uv exported no Python dependencies from the project lock file"
        )
    return requirements


def _yaml_string(value: str) -> str:
    """Return a JSON-quoted string, which is also a valid YAML scalar."""
    return json.dumps(value, ensure_ascii=False)


def render_conda_environment(
    configuration: CondaConfiguration,
    python_requirements: tuple[str, ...],
) -> str:
    """Render a deterministic Conda environment YAML document."""
    lines = [
        "# Generated by scripts/sync_conda_environment.py.",
        "# Sources: pyproject.toml and uv.lock. Do not edit this file manually.",
        f"name: {_yaml_string(configuration.environment_name)}",
        "channels:",
    ]
    lines.extend(f"  - {_yaml_string(channel)}" for channel in configuration.channels)
    lines.append("dependencies:")
    lines.append(
        f"  - {_yaml_string(f'python {configuration.python_requirement}')}"
    )
    lines.append('  - "pip"')
    lines.extend(
        f"  - {_yaml_string(dependency)}"
        for dependency in configuration.conda_dependencies
    )
    lines.append("  - pip:")
    lines.extend(
        f"      - {_yaml_string(requirement)}"
        for requirement in python_requirements
    )
    return "\n".join(lines) + "\n"


def write_if_changed(output_path: Path, content: str) -> bool:
    """Atomically replace ``output_path`` only when its content differs."""
    try:
        if output_path.read_text(encoding="utf-8") == content:
            return False
    except FileNotFoundError:
        pass

    temporary_file_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=output_path.parent,
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary_file:
            temporary_file.write(content)
            temporary_file.flush()
            os.fsync(temporary_file.fileno())
            temporary_file_path = Path(temporary_file.name)
        temporary_file_path.replace(output_path)
    except OSError:
        if temporary_file_path is not None:
            temporary_file_path.unlink(missing_ok=True)
        raise
    return True


def generate_conda_environment(project_root: Path = PROJECT_ROOT) -> str:
    """Generate the expected Conda environment document for a project root."""
    configuration = load_conda_configuration(project_root / PYPROJECT_FILENAME)
    python_requirements = export_locked_python_requirements(project_root)
    return render_conda_environment(configuration, python_requirements)


def _build_argument_parser() -> argparse.ArgumentParser:
    """Create the command-line parser."""
    argument_parser = argparse.ArgumentParser(
        description=(
            "Synchronize environment.yml with pyproject.toml and the uv lock file."
        )
    )
    argument_parser.add_argument(
        "--check",
        action="store_true",
        help="check that environment.yml is current without changing it",
    )
    return argument_parser


def main(arguments: list[str] | None = None) -> int:
    """Synchronize or check the generated Conda environment specification."""
    command_line_arguments = _build_argument_parser().parse_args(arguments)
    output_path = PROJECT_ROOT / CONDA_ENVIRONMENT_FILENAME
    try:
        expected_content = generate_conda_environment(PROJECT_ROOT)
        if command_line_arguments.check:
            try:
                is_current = output_path.read_text(encoding="utf-8") == expected_content
            except FileNotFoundError:
                is_current = False
            if not is_current:
                print(
                    f"{CONDA_ENVIRONMENT_FILENAME} is missing or out of date; "
                    "run scripts/sync_conda_environment.py",
                    file=sys.stderr,
                )
                return 1
            print(f"{CONDA_ENVIRONMENT_FILENAME} is current")
            return 0

        changed = write_if_changed(output_path, expected_content)
    except (CondaEnvironmentGenerationError, OSError) as error:
        print(f"Could not synchronize {CONDA_ENVIRONMENT_FILENAME}: {error}", file=sys.stderr)
        return 2

    if changed:
        print(f"Updated {CONDA_ENVIRONMENT_FILENAME}")
    else:
        print(f"{CONDA_ENVIRONMENT_FILENAME} is already current")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
