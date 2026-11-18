import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full")


@app.cell
def _():
    import importlib as _importlib
    import sys as _sys
    from pathlib import Path

    import marimo as mo

    repository_root_directory_path = Path(__file__).resolve().parents[1]
    repository_root_text = str(repository_root_directory_path)
    if repository_root_text not in _sys.path:
        _sys.path.insert(0, repository_root_text)

    from functional import manifest as _manifest

    manifest_module_source = mo.watch.file(
        repository_root_directory_path / "functional" / "manifest.py"
    )
    manifest_module_source.read_bytes()
    _manifest = _importlib.reload(_manifest)

    load_genome_set_manifest = _manifest.load_genome_set_manifest
    validate_genome_set_manifest = _manifest.validate_genome_set_manifest
    genome_manifest_paths = tuple(
        sorted(
            (repository_root_directory_path / "metadata" / "genome_sets").glob(
                "*.toml"
            )
        )
    )
    return (
        Path,
        genome_manifest_paths,
        load_genome_set_manifest,
        mo,
        repository_root_directory_path,
        validate_genome_set_manifest,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        # Functional Scheme Performance

        Phase 2 evaluates how a scheme performs on representative genomes. It
        is intentionally separate from the intrinsic FASTA-only evaluation:
        recovery, call completeness, reproducibility, multicopy evidence, and
        tool agreement depend on a versioned genome set and explicit typing-run
        provenance.

        This notebook currently provides the **input-readiness boundary**. It
        validates representative-genome manifests and identifies the data still
        required before genome-dependent calculations can be implemented. It
        does not run typing tools or report placeholder performance metrics.
        """
    )


@app.cell
def _(genome_manifest_paths, mo):
    manifest_options = {
        manifest_path.name: str(manifest_path)
        for manifest_path in genome_manifest_paths
    }
    selected_manifest_input = mo.ui.dropdown(
        options=manifest_options,
        value=(next(iter(manifest_options)) if manifest_options else None),
        label="Representative genome-set manifest",
        searchable=True,
    )
    mo.vstack(
        [
            mo.md(
                """
                ## Genome-set readiness

                Select a TOML manifest from `metadata/genome_sets/`. The
                included template is deliberately incomplete and shows the
                required provenance fields. Restart the application after
                adding a new manifest so it appears in this list.
                """
            ),
            selected_manifest_input,
        ],
        gap=1,
    )
    return (selected_manifest_input,)


@app.cell
def _(
    Path,
    load_genome_set_manifest,
    mo,
    selected_manifest_input,
    validate_genome_set_manifest,
):
    mo.stop(
        selected_manifest_input.value is None,
        mo.md("No genome-set manifest is available."),
    )
    selected_manifest_path = Path(selected_manifest_input.value)
    selected_manifest_source = mo.watch.file(selected_manifest_path)
    selected_manifest_source.read_bytes()
    genome_set_manifest = load_genome_set_manifest(selected_manifest_path)
    genome_manifest_validation = validate_genome_set_manifest(
        genome_set_manifest
    )
    return genome_manifest_validation, genome_set_manifest


@app.cell(hide_code=True)
def _(genome_manifest_validation, genome_set_manifest, mo):
    readiness_label = (
        "Ready for a defined analysis importer"
        if genome_manifest_validation.is_ready
        else "Preparation required"
    )
    issue_rows = [
        {
            "Severity": issue.severity.title(),
            "Location": issue.location,
            "Action": issue.message,
        }
        for issue in genome_manifest_validation.issues
    ]
    issue_table = (
        mo.ui.table(
            issue_rows,
            selection=None,
            pagination=True,
            page_size=15,
            show_search=True,
            show_download=False,
            max_columns=None,
            max_height=500,
        )
        if issue_rows
        else mo.md("No manifest-readiness issues were found.")
    )
    mo.vstack(
        [
            mo.md(
                f"""
                ### {readiness_label}

                - **Dataset ID:** `{genome_set_manifest.dataset_id or "Not supplied"}`
                - **Scheme ID:** `{genome_set_manifest.scheme_id or "Not supplied"}`
                - **Dataset version:** `{genome_set_manifest.dataset_version or "Not supplied"}`
                - **Genome records:** {len(genome_set_manifest.genomes):,}
                - **Blocking errors:** {genome_manifest_validation.error_count:,}
                - **Provenance warnings:** {genome_manifest_validation.warning_count:,}
                - **Manifest SHA-256:** `{genome_set_manifest.manifest_sha256}`

                A structurally ready genome manifest is necessary but not
                sufficient to run Phase 2. The allele-call input format and
                typing-run provenance must also be defined.
                """
            ),
            issue_table,
        ],
        gap=1,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        ---

        ## Planned analytical sequence

        1. **Input and provenance checks** — verify the genome-set version,
           assembly checksums, scheme identity, typing tools, tool versions,
           parameters, and allele-call state vocabulary.
        2. **Recovery and completeness** — distinguish per-locus prevalence
           across genomes from the percentage of scheme targets recovered in
           one genome. Missing, failed, partial, and intentionally excluded
           calls must remain distinct.
        3. **Multicopy evidence** — summarize genome-level multiple-hit or
           multiple-call states per locus. These results can support or refute
           intrinsic candidate-paralogy flags; sequence similarity alone is
           not copy-number evidence.
        4. **Reproducibility** — compare explicitly paired repeat runs or
           depth conditions. The repeat unit and denominator must be declared
           before a percentage is calculated.
        5. **Typing-tool comparison** — compare Juno, chewBBACA, SeqSphere,
           and source-implementation calls only after their call states and
           allele identifiers have a documented normalization contract.
        6. **Combined interpretation** — join Phase 2 findings to intrinsic
           scheme measurements by stable scheme and locus identifiers without
           merging their independent caches or notebooks.

        Detailed input, cache, missingness, and FAIR export decisions are in
        `docs/functional-workflow.md`.
        """
    )


@app.cell(hide_code=True)
def _(mo):
    mo.callout(
        mo.md(
            """
            **Information needed before implementing calculations**

            - Which assemblies form each representative set and why.
            - Where the versioned assembly files will live at runtime.
            - Which typing-tool outputs already exist and their exact formats.
            - The call-state vocabulary used by each tool.
            - Which runs are genuine repeats or sequencing-depth comparisons.
            - Which lineage and assembly-quality fields are available.
            - Whether SeqSphere or other restricted outputs can be redistributed.
            """
        ),
        kind="warn",
    )


if __name__ == "__main__":
    app.run()
