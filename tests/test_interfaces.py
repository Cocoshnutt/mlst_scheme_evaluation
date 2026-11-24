"""Tests for current command-line and Marimo notebook interface contracts."""

import ast
import contextlib
import io
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

from scheme import alignment, allele_info, bsr, cache, info, similarity

REPOSITORY_ROOT_DIRECTORY = Path(__file__).resolve().parents[1]
NOTEBOOK_FILE_PATH = (
    REPOSITORY_ROOT_DIRECTORY / "notebooks" / "intrinsic.py"
)
FUNCTIONAL_NOTEBOOK_FILE_PATH = (
    REPOSITORY_ROOT_DIRECTORY / "notebooks" / "functional.py"
)


class CommandLineHelpTest(unittest.TestCase):
    """Ensure every documented module CLI exposes help without doing analysis."""

    def test_documented_module_help_contracts(self) -> None:
        """Protect command availability and the principal option names."""
        module_help_contracts = (
            (cache, ("--cache-root", "--clear-cache")),
            (info, ("--json", "--ml")),
            (allele_info, ("--ml",)),
            (alignment, ("--strategy", "--workers")),
            (similarity, ("--threads", "--json")),
            (bsr, ("--bsr", "--threads", "--json")),
        )

        for command_line_module, expected_options in module_help_contracts:
            with self.subTest(module=command_line_module.__name__):
                output_stream = io.StringIO()
                with (
                    patch.object(
                        sys,
                        "argv",
                        [command_line_module.__name__, "--help"],
                    ),
                    contextlib.redirect_stdout(output_stream),
                    self.assertRaises(SystemExit) as exit_context,
                ):
                    command_line_module.main()

                self.assertEqual(exit_context.exception.code, 0)
                help_output = output_stream.getvalue()
                self.assertIn("usage:", help_output)
                for expected_option in expected_options:
                    self.assertIn(expected_option, help_output)


def _is_marimo_cell(function_node: ast.FunctionDef) -> bool:
    """Return whether a function has an ``@app.cell`` decorator."""
    for decorator in function_node.decorator_list:
        decorator_target = (
            decorator.func if isinstance(decorator, ast.Call) else decorator
        )
        if (
            isinstance(decorator_target, ast.Attribute)
            and isinstance(decorator_target.value, ast.Name)
            and decorator_target.value.id == "app"
            and decorator_target.attr == "cell"
        ):
            return True
    return False


def _called_function_names(function_node: ast.FunctionDef) -> set[str]:
    """Collect directly named calls made anywhere in one cell body."""
    return {
        call_node.func.id
        for call_node in ast.walk(function_node)
        if isinstance(call_node, ast.Call)
        and isinstance(call_node.func, ast.Name)
    }


def _called_attribute_names(function_node: ast.FunctionDef) -> set[str]:
    """Collect attribute-call names made anywhere in one cell body."""
    return {
        call_node.func.attr
        for call_node in ast.walk(function_node)
        if isinstance(call_node, ast.Call)
        and isinstance(call_node.func, ast.Attribute)
    }


def _notebook_cell_nodes(
    notebook_file_path: Path = NOTEBOOK_FILE_PATH,
) -> tuple[ast.FunctionDef, ...]:
    """Parse and return the Marimo cell functions in notebook order."""
    notebook_syntax_tree = ast.parse(
        notebook_file_path.read_text(encoding="utf-8"),
        filename=str(notebook_file_path),
    )
    return tuple(
        function_node
        for function_node in notebook_syntax_tree.body
        if isinstance(function_node, ast.FunctionDef)
        and _is_marimo_cell(function_node)
    )


class NotebookDependencyContractTest(unittest.TestCase):
    """Protect threshold-independent work from reactive classification controls."""

    def test_expensive_analysis_and_classification_use_separate_cells(self) -> None:
        """Keep threshold changes from directly rerunning MAFFT or BLAST."""
        cell_call_sets = [
            _called_function_names(function_node)
            for function_node in _notebook_cell_nodes()
        ]
        analysis_classification_pairs = (
            (
                "get_scheme_alignment_analysis",
                "calculate_scheme_alignment_summary",
            ),
            (
                "get_scheme_similarity_analysis",
                "classify_candidate_paralogous_locus_pairs",
            ),
            (
                "get_scheme_representative_bsr_analysis",
                "classify_representative_bsr_analysis",
            ),
        )

        for analysis_function, classification_function in (
            analysis_classification_pairs
        ):
            with self.subTest(analysis=analysis_function):
                analysis_cells = [
                    call_set
                    for call_set in cell_call_sets
                    if analysis_function in call_set
                ]
                classification_cells = [
                    call_set
                    for call_set in cell_call_sets
                    if classification_function in call_set
                ]
                self.assertEqual(len(analysis_cells), 1)
                self.assertEqual(len(classification_cells), 1)
                self.assertNotIn(
                    classification_function,
                    analysis_cells[0],
                )
                self.assertNotIn(
                    analysis_function,
                    classification_cells[0],
                )

    def test_expensive_analyses_require_individual_or_complete_run(self) -> None:
        """Keep MAFFT and BLAST idle while exposing pending overview state."""
        analysis_run_contracts = (
            (
                "get_scheme_alignment_analysis",
                "run_alignment_analysis_button",
                "scheme_alignment_analysis",
            ),
            (
                "get_scheme_similarity_analysis",
                "run_similarity_analysis_button",
                "scheme_similarity_analysis",
            ),
            (
                "get_scheme_representative_bsr_analysis",
                "run_representative_bsr_analysis_button",
                "scheme_representative_bsr_analysis",
            ),
        )

        for analysis_function, individual_run_button, result_name in (
            analysis_run_contracts
        ):
            analysis_cells = [
                function_node
                for function_node in _notebook_cell_nodes()
                if analysis_function in _called_function_names(function_node)
            ]
            with self.subTest(analysis=analysis_function):
                self.assertEqual(len(analysis_cells), 1)
                analysis_cell = analysis_cells[0]
                parameter_names = {
                    argument.arg for argument in analysis_cell.args.args
                }
                self.assertIn(individual_run_button, parameter_names)
                self.assertIn(
                    "run_complete_evaluation_button",
                    parameter_names,
                )
                self.assertNotIn("stop", _called_attribute_names(analysis_cell))
                pending_assignments = [
                    syntax_node
                    for syntax_node in ast.walk(analysis_cell)
                    if isinstance(syntax_node, ast.Assign)
                    and any(
                        isinstance(target, ast.Name)
                        and target.id == result_name
                        for target in syntax_node.targets
                    )
                    and isinstance(syntax_node.value, ast.Constant)
                    and syntax_node.value.value is None
                ]
                self.assertEqual(len(pending_assignments), 1)

    def test_expensive_analyses_share_the_notebook_worker_budget(self) -> None:
        """Pass one visible resource limit to MAFFT and both BLAST stages."""
        analysis_keyword_contracts = (
            ("get_scheme_alignment_analysis", "maximum_parallel_loci"),
            ("get_scheme_similarity_analysis", "maximum_threads"),
            ("get_scheme_representative_bsr_analysis", "maximum_threads"),
        )

        for analysis_function, expected_keyword in analysis_keyword_contracts:
            analysis_cells = [
                function_node
                for function_node in _notebook_cell_nodes()
                if analysis_function in _called_function_names(function_node)
            ]
            with self.subTest(analysis=analysis_function):
                self.assertEqual(len(analysis_cells), 1)
                analysis_cell = analysis_cells[0]
                self.assertIn(
                    "maximum_worker_count",
                    {argument.arg for argument in analysis_cell.args.args},
                )
                analysis_calls = [
                    call_node
                    for call_node in ast.walk(analysis_cell)
                    if isinstance(call_node, ast.Call)
                    and isinstance(call_node.func, ast.Name)
                    and call_node.func.id == analysis_function
                ]
                self.assertEqual(len(analysis_calls), 1)
                worker_keywords = [
                    keyword
                    for keyword in analysis_calls[0].keywords
                    if keyword.arg == expected_keyword
                ]
                self.assertEqual(len(worker_keywords), 1)
                self.assertIsInstance(worker_keywords[0].value, ast.Name)
                self.assertEqual(
                    worker_keywords[0].value.id,
                    "maximum_worker_count",
                )

    def test_static_analyses_do_not_depend_on_run_controls(self) -> None:
        """Load the inexpensive first-stage results automatically."""
        static_analysis_functions = (
            "get_scheme_information",
            "get_scheme_statistics",
            "get_scheme_allele_information",
        )

        for analysis_function in static_analysis_functions:
            analysis_cells = [
                function_node
                for function_node in _notebook_cell_nodes()
                if analysis_function in _called_function_names(function_node)
            ]
            with self.subTest(analysis=analysis_function):
                self.assertEqual(len(analysis_cells), 1)
                parameter_names = {
                    argument.arg for argument in analysis_cells[0].args.args
                }
                self.assertFalse(
                    any("run_" in parameter_name for parameter_name in parameter_names)
                )
                self.assertNotIn(
                    "stop",
                    _called_attribute_names(analysis_cells[0]),
                )


class FunctionalNotebookContractTest(unittest.TestCase):
    """Keep Phase 2 honest while only its input boundary is implemented."""

    def test_notebook_validates_manifests_without_analysis_controls(self) -> None:
        """Do not imply that genome-dependent calculations already exist."""
        functional_cells = _notebook_cell_nodes(FUNCTIONAL_NOTEBOOK_FILE_PATH)
        called_functions = set().union(
            *(_called_function_names(cell) for cell in functional_cells)
        )
        called_attributes = set().union(
            *(_called_attribute_names(cell) for cell in functional_cells)
        )

        self.assertIn("load_genome_set_manifest", called_functions)
        self.assertIn("validate_genome_set_manifest", called_functions)
        self.assertNotIn("run_button", called_attributes)
        self.assertFalse(
            any(
                function_name.startswith("get_scheme_")
                for function_name in called_functions
            )
        )

    def test_expensive_controls_ignore_reactive_static_thresholds(self) -> None:
        """Keep static display thresholds from invalidating expensive results."""
        readiness_cells = [
            function_node
            for function_node in _notebook_cell_nodes()
            if any(
                isinstance(syntax_node, ast.Name)
                and isinstance(syntax_node.ctx, ast.Store)
                and syntax_node.id == "static_analysis_ready"
                for syntax_node in ast.walk(function_node)
            )
        ]

        self.assertEqual(len(readiness_cells), 1)
        self.assertEqual(
            tuple(
                argument.arg
                for argument in readiness_cells[0].args.args
            ),
            ("scheme_information",),
        )

    def test_presentation_proxy_dependencies_are_not_reintroduced(self) -> None:
        """Use scientific results, not rendered-section flags, as dependencies."""
        notebook_names = {
            syntax_node.id
            for function_node in _notebook_cell_nodes()
            for syntax_node in ast.walk(function_node)
            if isinstance(syntax_node, ast.Name)
        }
        parameter_names = {
            argument.arg
            for function_node in _notebook_cell_nodes()
            for argument in function_node.args.args
        }

        for name in notebook_names | parameter_names:
            self.assertFalse(name.endswith("_displayed"), name)
            self.assertFalse(name.endswith("_section_ready"), name)

    def test_progressive_overview_uses_one_library_aggregation_cell(self) -> None:
        """Keep overview aggregation separate from rendering and analyses."""
        aggregation_cells = [
            function_node
            for function_node in _notebook_cell_nodes()
            if "create_intrinsic_scheme_overview"
            in _called_function_names(function_node)
        ]
        self.assertEqual(len(aggregation_cells), 1)
        aggregation_cell = aggregation_cells[0]
        parameter_names = {
            argument.arg for argument in aggregation_cell.args.args
        }
        for optional_result_name in (
            "scheme_alignment_analysis",
            "scheme_alignment_summary",
            "scheme_similarity_analysis",
            "candidate_paralog_analysis",
            "scheme_representative_bsr_analysis",
            "representative_bsr_classification",
        ):
            self.assertIn(optional_result_name, parameter_names)

        overview_render_cells = [
            function_node
            for function_node in _notebook_cell_nodes()
            if "intrinsic_scheme_overview"
            in {argument.arg for argument in function_node.args.args}
            and "create_intrinsic_evaluation_bundle"
            not in _called_function_names(function_node)
        ]
        self.assertEqual(len(overview_render_cells), 1)
        first_expensive_control_line = min(
            function_node.lineno
            for function_node in _notebook_cell_nodes()
            if any(
                isinstance(syntax_node, ast.Name)
                and isinstance(syntax_node.ctx, ast.Store)
                and syntax_node.id == "run_complete_evaluation_button"
                for syntax_node in ast.walk(function_node)
            )
        )
        self.assertLess(
            overview_render_cells[0].lineno,
            first_expensive_control_line,
        )


if __name__ == "__main__":
    unittest.main()
