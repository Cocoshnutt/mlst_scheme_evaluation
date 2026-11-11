import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full")


@app.cell
def _():
    import importlib as _importlib
    from pathlib import Path

    import marimo as mo

    import schema_cache as _schema_cache
    import schema_info as _schema_info
    import schema_stats as _schema_stats
    import scheme_allele_info as _scheme_allele_info

    _scheme_allele_info_source_file = mo.watch.file("scheme_allele_info.py")
    _schema_cache_source_file = mo.watch.file("schema_cache.py")
    _schema_info_source_file = mo.watch.file("schema_info.py")
    _schema_stats_source_file = mo.watch.file("schema_stats.py")
    _scheme_allele_info_source_file.read_bytes()
    _schema_cache_source_file.read_bytes()
    _schema_info_source_file.read_bytes()
    _schema_stats_source_file.read_bytes()

    _scheme_allele_info = _importlib.reload(_scheme_allele_info)
    _schema_cache = _importlib.reload(_schema_cache)
    _schema_info = _importlib.reload(_schema_info)
    _schema_stats = _importlib.reload(_schema_stats)

    cache_archive = _schema_cache.cache_archive
    clear_cache = _schema_cache.clear_cache
    count_alleles_below_minimum_length = (
        _schema_info.count_alleles_below_minimum_length
    )
    create_schema_statistics_charts = (
        _schema_stats.create_schema_statistics_charts
    )
    get_scheme_allele_information = (
        _scheme_allele_info.get_scheme_allele_information
    )
    get_schema_information = _schema_info.get_schema_information
    get_schema_statistics = _schema_stats.get_schema_statistics
    analysis_cache_filenames = (
        _schema_info.SCHEMA_INFORMATION_CACHE_FILENAME,
        _schema_stats.SCHEMA_STATISTICS_CACHE_FILENAME,
        _scheme_allele_info.SCHEME_ALLELE_INFORMATION_CACHE_FILENAME,
    )

    return (
        Path,
        analysis_cache_filenames,
        cache_archive,
        clear_cache,
        count_alleles_below_minimum_length,
        create_schema_statistics_charts,
        get_scheme_allele_information,
        get_schema_information,
        get_schema_statistics,
        mo,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.Html(r"""
    <style>
    @import url("https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap");

    :root,
    .marimo {
        --marimo-text-font: "Inter", sans-serif;
        --marimo-heading-font: "Inter", sans-serif;
    }
    </style>
    """)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Schema Evaluation

    ## Schema Summary Data
    - Loci: Total number of loci that were evaluated.
    - Alleles: Total number of alleles.
    - Valid Alleles: Total number of valid alleles. An allele is considered valid if its sequence size is a multiple of 3, if it has a single start and stop codon, if has no in-frame stop codons, and if it contains no ambiguous bases.
    - Invalid Alleles: Sum of the Incomplete ORF, Ambiguous Bases, Missing Start/Stop Codon, and In-frame Stop Codon counts.
    - Incomplete ORF: Alleles whose sequence size is not a multiple of 3.
    - Ambiguous Bases: Alleles containing characters other than A, C, T, or G.
    - Missing Start/Stop Codon: Alleles missing a recognized start codon, stop codon, or both.
    - In-frame Stop Codon: Alleles containing an in-frame stop codon before the sequence endpoint.
    - Alleles below minimum length: Alleles shorter than the selected minimum sequence length in nucleotides.
    - Length outliers: Scheme-wide totals of alleles strictly below 80% or above 120% of their respective locus's modal allele length.
    """)


@app.cell
def _(mo):
    get_schema_reset, set_schema_reset = mo.state(0)
    return get_schema_reset, set_schema_reset


@app.cell
def _(Path, get_schema_reset, mo):
    get_schema_reset()
    schema_paths = sorted(
        archive_file_path
        for archive_file_path in Path("schemes").rglob("*.tar.*")
        if "BUSCO_bac" not in archive_file_path.parts
    )
    schema_options = {
        str(archive_file_path): archive_file_path
        for archive_file_path in schema_paths
    }

    schema = mo.ui.dropdown(
        options=schema_options,
        label="Scheme archive",
        searchable=True,
    )
    mo.vstack([schema])
    return (schema,)


@app.cell
def _(analysis_cache_filenames, cache_archive, mo, schema):
    mo.stop(schema.value is None, mo.md("Select a scheme archive to continue."))

    scheme_path = cache_archive(schema.value)
    unavailable_analysis_cache_count = sum(
        not (scheme_path / cache_filename).is_file()
        for cache_filename in analysis_cache_filenames
    )
    mo.md(f"Using cached scheme at `{scheme_path}`")
    return scheme_path, unavailable_analysis_cache_count


@app.cell(hide_code=True)
def _(mo, unavailable_analysis_cache_count):
    if unavailable_analysis_cache_count:
        mo.callout(
            "Analysis results are not fully cached for this scheme. The "
            "initial computation may take longer for schemes with many loci "
            "or alleles. Results will be cached for subsequent sessions.",
            kind="warn",
        )
    else:
        mo.callout(
            "Cached analysis results are available for this scheme, so the "
            "statistics and allele table can be loaded without rescanning "
            "the FASTA files.",
            kind="info",
        )


@app.cell
def _(mo):
    minimum_sequence_length_input = mo.ui.number(
        start=0,
        step=1,
        value=0,
        debounce=True,
        label="Minimum allele length (bp)",
    )
    mo.vstack([minimum_sequence_length_input])
    return (minimum_sequence_length_input,)


@app.cell
def _(get_schema_information, scheme_path):
    schema_information = get_schema_information(scheme_path)
    return (schema_information,)


@app.cell
def _(
    count_alleles_below_minimum_length,
    minimum_sequence_length_input,
    scheme_path,
):
    minimum_sequence_length = int(minimum_sequence_length_input.value or 0)
    alleles_below_minimum_length_count = count_alleles_below_minimum_length(
        scheme_path,
        minimum_sequence_length=minimum_sequence_length,
    )
    return alleles_below_minimum_length_count, minimum_sequence_length


@app.cell
def _(
    alleles_below_minimum_length_count,
    minimum_sequence_length,
    mo,
    schema_information,
):
    mo.md(
        f"""
        ## Selected schema summary data

        - **Loci:** {schema_information.locus_count:,}
        - **Alleles:** {schema_information.allele_count:,}
        - **Valid alleles:** {schema_information.valid_allele_count:,}
        - **Invalid alleles:** {schema_information.invalid_allele_count:,}
        - **Incomplete ORF:** {schema_information.incomplete_orf_allele_count:,}
        - **Ambiguous bases:** {schema_information.ambiguous_bases_allele_count:,}
        - **Missing start/stop codon:** {schema_information.missing_start_or_stop_codon_allele_count:,}
        - **In-frame stop codon:** {schema_information.in_frame_stop_codon_allele_count:,}
        - **Alleles <{minimum_sequence_length:,} bp:** {alleles_below_minimum_length_count:,}
        - **Alleles <80% of locus mode:** {schema_information.alleles_below_locus_mode_threshold_count:,} ({schema_information.alleles_below_locus_mode_threshold_percentage:.2f}%)
        - **Alleles >120% of locus mode:** {schema_information.alleles_above_locus_mode_threshold_count:,} ({schema_information.alleles_above_locus_mode_threshold_percentage:.2f}%)
        """
    )


@app.cell
def _(get_schema_statistics, scheme_path):
    locus_statistics = get_schema_statistics(scheme_path)
    return (locus_statistics,)


@app.cell
def _(create_schema_statistics_charts, locus_statistics):
    schema_statistics_charts = create_schema_statistics_charts(locus_statistics)
    return (schema_statistics_charts,)


@app.cell(hide_code=True)
def _(locus_statistics, mo, schema_statistics_charts):
    total_alleles_chart = mo.ui.altair_chart(
        schema_statistics_charts.total_alleles_chart,
        chart_selection=False,
        legend_selection=False,
    )
    allele_mode_size_chart = mo.ui.altair_chart(
        schema_statistics_charts.allele_mode_size_chart,
        chart_selection=False,
        legend_selection=False,
    )
    locus_statistics_chart = mo.ui.altair_chart(
        schema_statistics_charts.locus_statistics_chart,
        chart_selection=False,
        legend_selection=False,
    )
    allele_size_variation_chart = mo.ui.altair_chart(
        schema_statistics_charts.allele_size_variation_chart,
        chart_selection=False,
        legend_selection=False,
    )
    mo.vstack(
        [
            mo.md(
                f"""
                ---

                # Locus Statistics
                The following panels summarize allele counts and sequence-length
                variation across the **{len(locus_statistics):,} evaluated loci**.
                """
            ),
            total_alleles_chart,
            allele_mode_size_chart,
            locus_statistics_chart,
            allele_size_variation_chart,
        ],
        gap=2,
    )


@app.cell
def _(
    get_scheme_allele_information,
    minimum_sequence_length,
    scheme_path,
):
    scheme_allele_information = get_scheme_allele_information(
        scheme_path,
        minimum_sequence_length=minimum_sequence_length,
    )
    return (scheme_allele_information,)


@app.cell(hide_code=True)
def _(minimum_sequence_length, mo, scheme_allele_information):
    minimum_length_column_name = f"Alleles <{minimum_sequence_length}bp"
    allele_information_rows = [
        {
            "Locus": locus_information.locus_name,
            "Total Alleles": locus_information.total_allele_count,
            "Valid Alleles": locus_information.valid_allele_count,
            "Invalid Alleles": locus_information.invalid_allele_count,
            "Proportion of Validated Alleles": round(
                locus_information.validated_allele_proportion,
                3,
            ),
            "Distinct Protein Alleles": (
                locus_information.distinct_protein_allele_count
            ),
            "Incomplete ORF": locus_information.incomplete_orf_allele_count,
            "Ambiguous Bases": (
                locus_information.ambiguous_bases_allele_count
            ),
            "Missing Start/Stop Codon": (
                locus_information.missing_start_or_stop_codon_allele_count
            ),
            "In-frame Stop Codon": (
                locus_information.in_frame_stop_codon_allele_count
            ),
            minimum_length_column_name: (
                locus_information.alleles_below_minimum_length_count
            ),
            "Alleles below threshold": (
                locus_information.alleles_below_mode_threshold_count
            ),
            "Alleles above threshold": (
                locus_information.alleles_above_mode_threshold_count
            ),
            "Missing Allele IDs": (
                locus_information.missing_allele_identifier_count
            ),
        }
        for locus_information in scheme_allele_information
    ]
    allele_information_table = mo.ui.table(
        allele_information_rows,
        pagination=True,
        page_size=10,
        selection=None,
        show_column_summaries=False,
        show_data_types=False,
        format_mapping={
            "Invalid Alleles": "{:d}",
            "Proportion of Validated Alleles": "{:.3f}",
        },
        freeze_columns_left=["Locus"],
        text_justify_columns={
            column_name: "right"
            for column_name in allele_information_rows[0]
            if column_name != "Locus"
        }
        if allele_information_rows
        else None,
        column_widths={
            "Locus": 240,
            "Total Alleles": 115,
            "Valid Alleles": 115,
            "Invalid Alleles": 120,
            "Proportion of Validated Alleles": 215,
            "Distinct Protein Alleles": 165,
            "Incomplete ORF": 125,
            "Ambiguous Bases": 130,
            "Missing Start/Stop Codon": 190,
            "In-frame Stop Codon": 145,
            minimum_length_column_name: 135,
            "Alleles below threshold": 165,
            "Alleles above threshold": 165,
            "Missing Allele IDs": 145,
        },
        header_tooltip={
            "Invalid Alleles": (
                "Sum of Incomplete ORF, Ambiguous Bases, Missing Start/Stop "
                "Codon, and In-frame Stop Codon counts."
            ),
            "Distinct Protein Alleles": (
                "Distinct proteins translated from valid allele sequences."
            ),
            "Alleles below threshold": (
                "Alleles strictly below 80% of the locus modal length."
            ),
            "Alleles above threshold": (
                "Alleles strictly above 120% of the locus modal length."
            ),
            "Missing Allele IDs": (
                "Missing positive integer IDs from 1 through the largest "
                "observed ID; blank means an ID could not be interpreted."
            ),
        },
        show_search=True,
        show_download=True,
        max_columns=None,
        max_height=600,
    )
    mo.vstack(
        [
            mo.md(
                """
                ---

                # Allele Analysis

                Per-locus allele integrity and diversity metrics.
                """
            ),
            allele_information_table,
        ],
        gap=2,
    )


@app.cell
def _(mo):
    clear_cache_button = mo.ui.run_button(
        kind="danger",
        label="Clear scheme cache",
    )
    mo.vstack(
        [
            mo.callout(
                "Clearing the cache removes extracted schemes and all computed "
                "analysis results. The next session or scheme selection may "
                "therefore take longer while those results are rebuilt.",
                kind="warn",
            ),
            clear_cache_button,
        ]
    )
    return (clear_cache_button,)


@app.cell
def _(clear_cache, clear_cache_button, mo, set_schema_reset):
    if clear_cache_button.value:
        clear_cache()
        set_schema_reset(lambda value: value + 1)
        mo.status.toast(
            "Scheme cache cleared",
            "The archive selection was reset.",
        )


if __name__ == "__main__":
    app.run()
