import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full")


@app.cell
def _():
    import importlib as _importlib
    import sys as _sys
    from pathlib import Path

    import altair as alt
    import marimo as mo

    from schema import alignment as _alignment
    from schema import allele_info as _allele_info
    from schema import cache as _cache
    from schema import info as _info
    from schema import similarity as _similarity
    from schema import stats as _stats

    _alignment_source_file = mo.watch.file("schema/alignment.py")
    _allele_info_source_file = mo.watch.file("schema/allele_info.py")
    _cache_source_file = mo.watch.file("schema/cache.py")
    _info_source_file = mo.watch.file("schema/info.py")
    _similarity_source_file = mo.watch.file("schema/similarity.py")
    _stats_source_file = mo.watch.file("schema/stats.py")
    _alignment_source_file.read_bytes()
    _allele_info_source_file.read_bytes()
    _cache_source_file.read_bytes()
    _info_source_file.read_bytes()
    _similarity_source_file.read_bytes()
    _stats_source_file.read_bytes()

    _cache = _importlib.reload(_cache)
    _info = _importlib.reload(_info)
    _allele_info = _importlib.reload(_allele_info)
    _stats = _importlib.reload(_stats)
    _alignment = _importlib.reload(_alignment)
    _similarity = _importlib.reload(_similarity)

    calculate_schema_alignment_summary = (
        _alignment.calculate_schema_alignment_summary
    )
    cache_archive = _cache.cache_archive
    clear_cache = _cache.clear_cache
    count_alleles_below_minimum_length = _info.count_alleles_below_minimum_length
    create_schema_statistics_charts = _stats.create_schema_statistics_charts
    default_amino_acid_gap_warning_threshold_percentage = (
        _alignment.DEFAULT_AMINO_ACID_GAP_WARNING_THRESHOLD_PERCENTAGE
    )
    default_alignment_strategy = _alignment.DEFAULT_ALIGNMENT_STRATEGY
    alignment_strategy_labels = _alignment.ALIGNMENT_STRATEGY_LABELS
    get_scheme_allele_information = _allele_info.get_scheme_allele_information
    get_scheme_sequence_integrity_information = (
        _allele_info.get_scheme_sequence_integrity_information
    )
    get_schema_information = _info.get_schema_information
    get_schema_alignment_analysis = _alignment.get_schema_alignment_analysis
    get_schema_similarity_analysis = (
        _similarity.get_schema_similarity_analysis
    )
    get_schema_statistics = _stats.get_schema_statistics
    mafft_progress_output_stream = _sys.__stdout__ or _sys.stdout
    similarity_progress_output_stream = _sys.__stdout__ or _sys.stdout
    analysis_cache_filenames = (
        _alignment.SCHEMA_ALIGNMENT_CACHE_FILENAME,
        _info.SCHEMA_INFORMATION_CACHE_FILENAME,
        _similarity.SCHEMA_SIMILARITY_CACHE_FILENAME,
        _stats.SCHEMA_STATISTICS_CACHE_FILENAME,
        _allele_info.SCHEME_ALLELE_INFORMATION_CACHE_FILENAME,
    )

    return (
        Path,
        alt,
        alignment_strategy_labels,
        analysis_cache_filenames,
        calculate_schema_alignment_summary,
        cache_archive,
        clear_cache,
        count_alleles_below_minimum_length,
        create_schema_statistics_charts,
        default_amino_acid_gap_warning_threshold_percentage,
        default_alignment_strategy,
        get_scheme_allele_information,
        get_scheme_sequence_integrity_information,
        get_schema_alignment_analysis,
        get_schema_information,
        get_schema_similarity_analysis,
        get_schema_statistics,
        mafft_progress_output_stream,
        mo,
        similarity_progress_output_stream,
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
    - Empty Sequences: FASTA records containing no nucleotides. A locus with no FASTA records is also counted as affected.
    - Within-locus Duplicate Nucleotide Sequences: Redundant non-empty allele records with an exact normalized sequence match inside the same locus. This does not compare different loci.
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
        for archive_file_path in Path("archive").rglob("*.tar.*")
        if "BUSCO_bac" not in archive_file_path.parts
    )
    schema_options = {
        str(archive_file_path): archive_file_path for archive_file_path in schema_paths
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

    with mo.status.spinner(
        title="Loading selected scheme archive",
        subtitle="Reusing the extraction cache when available.",
    ):
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
        _analysis_cache_notice = mo.callout(
            "Analysis results are not fully cached for this scheme. The "
            "initial computation may take longer for schemes with many loci "
            "or alleles. Results will be cached for subsequent sessions.",
            kind="warn",
        )
    else:
        _analysis_cache_notice = mo.callout(
            "Analysis cache files are available for this scheme. Statistics "
            "and the allele table can be loaded without recalculating their "
            "static metrics; the alignment stage will still verify each "
            "per-locus cache entry and process changed, incomplete, or "
            "previously failed loci.",
            kind="info",
        )
    mo.vstack([_analysis_cache_notice])


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
def _(get_schema_information, mo, scheme_path):
    with mo.status.spinner(
        title="Calculating scheme summary",
        subtitle="Scanning allele records when no summary cache is available.",
    ):
        schema_information = get_schema_information(scheme_path)
    return (schema_information,)


@app.cell
def _(
    count_alleles_below_minimum_length,
    minimum_sequence_length_input,
    schema_information,
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
    schema_summary_displayed = True
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
    return (schema_summary_displayed,)


@app.cell(hide_code=True)
def _(mo, schema_summary_displayed):
    locus_statistics_section_ready = schema_summary_displayed
    mo.md(
        """
        ---

        # Locus Statistics

        These panels summarize allele counts and sequence-length variation
        across the evaluated loci.
        """
    )
    return (locus_statistics_section_ready,)


@app.cell
def _(
    create_schema_statistics_charts,
    get_schema_statistics,
    locus_statistics_section_ready,
    mo,
    scheme_path,
):
    with mo.status.spinner(
        title="Calculating locus statistics",
        subtitle="Loading per-locus allele-count and length distributions.",
    ):
        locus_statistics = get_schema_statistics(scheme_path)
        schema_statistics_charts = create_schema_statistics_charts(locus_statistics)
    return locus_statistics, schema_statistics_charts


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
    locus_statistics_displayed = True
    mo.vstack(
        [
            mo.md(f"**{len(locus_statistics):,} loci evaluated.**"),
            total_alleles_chart,
            allele_mode_size_chart,
            locus_statistics_chart,
            allele_size_variation_chart,
        ],
        gap=2,
    )
    return (locus_statistics_displayed,)


@app.cell(hide_code=True)
def _(locus_statistics_displayed, mo):
    allele_analysis_section_ready = locus_statistics_displayed
    mo.md(
        """
        ---

        # Allele Analysis

        Per-locus allele integrity and diversity metrics. Scheme-wide empty and
        duplicate sequence measurements are shown here because they are derived
        by the same cached allele-record analysis.
        """
    )
    return (allele_analysis_section_ready,)


@app.cell
def _(
    allele_analysis_section_ready,
    get_scheme_allele_information,
    get_scheme_sequence_integrity_information,
    minimum_sequence_length,
    mo,
    scheme_path,
):
    with mo.status.spinner(
        title="Calculating allele analysis",
        subtitle="Loading per-locus integrity and diversity measurements.",
    ):
        scheme_allele_information = get_scheme_allele_information(
            scheme_path,
            minimum_sequence_length=minimum_sequence_length,
        )
        scheme_sequence_integrity_information = (
            get_scheme_sequence_integrity_information(scheme_path)
        )
    return scheme_allele_information, scheme_sequence_integrity_information


@app.cell(hide_code=True)
def _(
    minimum_sequence_length,
    mo,
    scheme_allele_information,
    scheme_sequence_integrity_information,
):
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
            "Empty Sequences": locus_information.empty_sequence_count,
            "Within-locus Duplicate Sequences": (
                locus_information.within_locus_duplicate_sequence_count
            ),
            "Incomplete ORF": locus_information.incomplete_orf_allele_count,
            "Ambiguous Bases": (locus_information.ambiguous_bases_allele_count),
            "Missing Start/Stop Codon": (
                locus_information.missing_start_or_stop_codon_allele_count
            ),
            "In-frame Stop Codon": (locus_information.in_frame_stop_codon_allele_count),
            minimum_length_column_name: (
                locus_information.alleles_below_minimum_length_count
            ),
            "Alleles below threshold": (
                locus_information.alleles_below_mode_threshold_count
            ),
            "Alleles above threshold": (
                locus_information.alleles_above_mode_threshold_count
            ),
            "Missing Allele IDs": (locus_information.missing_allele_identifier_count),
        }
        for locus_information in scheme_allele_information
    ]
    allele_information_table = mo.ui.table(
        allele_information_rows,
        pagination=True,
        page_size=50,
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
            "Locus": 280,
            "Total Alleles": 130,
            "Valid Alleles": 130,
            "Invalid Alleles": 145,
            "Proportion of Validated Alleles": 245,
            "Distinct Protein Alleles": 200,
            "Empty Sequences": 165,
            "Within-locus Duplicate Sequences": 270,
            "Incomplete ORF": 150,
            "Ambiguous Bases": 165,
            "Missing Start/Stop Codon": 215,
            "In-frame Stop Codon": 185,
            minimum_length_column_name: 155,
            "Alleles below threshold": 200,
            "Alleles above threshold": 200,
            "Missing Allele IDs": 165,
        },
        header_tooltip={
            "Invalid Alleles": (
                "Sum of Incomplete ORF, Ambiguous Bases, Missing Start/Stop "
                "Codon, and In-frame Stop Codon counts."
            ),
            "Distinct Protein Alleles": (
                "Distinct proteins translated from valid allele sequences."
            ),
            "Empty Sequences": (
                "FASTA records containing no nucleotides. A locus with zero "
                "total alleles is also included in the scheme-level affected "
                "locus count."
            ),
            "Within-locus Duplicate Sequences": (
                "Redundant non-empty allele records beyond the first exact "
                "normalized nucleotide sequence within this locus. This is "
                "not a cross-locus comparison."
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
    allele_analysis_displayed = True
    mo.vstack(
        [
            mo.md(
                f"""
                ## Sequence Integrity

                - **Empty sequences:** {scheme_sequence_integrity_information.empty_sequence_count:,} records; {scheme_sequence_integrity_information.loci_with_empty_sequences_count:,} of {scheme_sequence_integrity_information.locus_count:,} loci affected ({scheme_sequence_integrity_information.loci_with_empty_sequences_percentage:.2f}%).
                - **Within-locus duplicate nucleotide sequences:** {scheme_sequence_integrity_information.within_locus_duplicate_sequence_count:,} redundant records; {scheme_sequence_integrity_information.loci_with_within_locus_duplicate_sequences_count:,} of {scheme_sequence_integrity_information.locus_count:,} loci affected ({scheme_sequence_integrity_information.loci_with_within_locus_duplicate_sequences_percentage:.2f}%).

                These duplicate measurements compare exact normalized sequences
                only within each locus; they are not cross-locus similarity or
                paralogy detection.
                """
            ),
            allele_information_table,
        ],
        gap=2,
    )
    return (allele_analysis_displayed,)


@app.cell
def _(
    allele_analysis_displayed,
    alignment_strategy_labels,
    default_amino_acid_gap_warning_threshold_percentage,
    default_alignment_strategy,
    mo,
):
    amino_acid_gap_warning_threshold_input = mo.ui.number(
        start=0,
        stop=100,
        step=0.1,
        value=default_amino_acid_gap_warning_threshold_percentage,
        debounce=True,
        label="Advisory amino-acid gap threshold (%)",
    )
    alignment_strategy_input = mo.ui.dropdown(
        options={
            alignment_strategy_label: alignment_strategy
            for alignment_strategy, alignment_strategy_label in (
                alignment_strategy_labels.items()
            )
        },
        value=alignment_strategy_labels[default_alignment_strategy],
        label="MAFFT alignment strategy",
    )
    alignment_section_ready = allele_analysis_displayed
    mo.vstack(
        [
            mo.md(
                f"""
                ---

                # Within-locus Alignment Gap Statistics

                These statistics align nucleotide alleles and translated
                protein alleles independently inside each locus. Independent
                alignments are used because nucleotide and protein-level
                variation answer different quality questions.

                **Strategy choice:** FFT-NS-2 is the fast default and rebuilds
                its guide tree twice. FFT-NS-1 is faster but uses a rougher,
                single guide tree. Auto lets MAFFT choose and may select a much
                slower iterative method. Because strategy choice can change an
                alignment and its measured gap proportion, cached results are
                kept separately for every strategy and the chosen strategy is
                reported with the results.

                **Method and scientific choices**

                - Loci are processed concurrently across every CPU available
                  to the Marimo process. Each MAFFT subprocess uses one thread
                  to avoid nested oversubscription.
                - Unchanged loci are reused from a versioned per-locus cache.
                  New, modified, interrupted, or previously failed loci are
                  processed and checkpointed incrementally.
                - If every eligible sequence is exactly identical, the exact
                  zero-gap result is returned without launching MAFFT.
                - Nucleotide alignments include every non-empty record with
                  standard IUPAC DNA symbols. Protein alignments include only
                  alleles that pass the existing CDS-integrity rules.
                - Exact duplicate records remain represented; no sequence is
                  deduplicated before a non-trivial alignment.
                - Gap proportion is `gap characters / (sequence count ×
                  alignment length) × 100`. Leading and trailing gaps count.
                - A nucleotide alignment column is **callable** when at least
                  two sequences contain canonical A/C/G/T observations. Gaps
                  and ambiguous IUPAC symbols are treated as missing data.
                - A callable column is **polymorphic** when its canonical
                  observations include at least two distinct bases. The
                  percentage is `polymorphic sites / callable sites × 100`.
                  This substitution-only statistic deliberately keeps indel
                  variation in the separate gap metrics above.
                - At least two eligible sequences are required; otherwise the
                  measurement is unavailable rather than zero.
                - These metrics do not depend on the **Minimum allele length**
                  setting. The advisory threshold below is applied to cached
                  proportions and never triggers realignment.
                - The default **greater than
                  {default_amino_acid_gap_warning_threshold_percentage:.1f}%**
                  amino-acid gap flag is the study-specific Guglielmini et al.
                  (2021) screening criterion, not a universal pass/fail rule.
                  The Savin et al. (2019) contiguous-gap criterion is not
                  calculated by this proportion metric.
                """
            ),
            alignment_strategy_input,
            amino_acid_gap_warning_threshold_input,
        ],
        gap=1,
    )
    return (
        alignment_section_ready,
        alignment_strategy_input,
        amino_acid_gap_warning_threshold_input,
    )


@app.cell
def _(
    alignment_section_ready,
    alignment_strategy_input,
    alignment_strategy_labels,
    default_alignment_strategy,
    get_schema_alignment_analysis,
    locus_statistics,
    mafft_progress_output_stream,
    mo,
    scheme_path,
):
    selected_alignment_strategy = (
        alignment_strategy_input.value or default_alignment_strategy
    )
    selected_alignment_strategy_label = alignment_strategy_labels[
        selected_alignment_strategy
    ]
    _completed_locus_count = [0]
    with mo.status.progress_bar(
        total=len(locus_statistics),
        title=(
            "Preparing nucleotide and amino-acid alignments with "
            f"{selected_alignment_strategy_label}"
        ),
        subtitle="Checking the per-locus cache.",
        completion_title="Alignment analysis complete",
        completion_subtitle="All current loci have alignment measurements.",
        show_rate=True,
        show_eta=True,
        remove_on_exit=False,
    ) as _alignment_progress_bar:

        def _update_alignment_progress(alignment_progress):
            progress_increment = max(
                0,
                alignment_progress.completed_locus_count - _completed_locus_count[0],
            )
            if alignment_progress.current_locus_name is None:
                progress_subtitle = (
                    f"Reused {alignment_progress.cached_locus_count:,} cached "
                    f"loci; "
                    f"{alignment_progress.total_locus_count - alignment_progress.cached_locus_count:,} "
                    "require processing."
                )
            else:
                progress_rate = (
                    f" at {alignment_progress.newly_completed_loci_per_second:.2f} loci/s"
                    if alignment_progress.newly_completed_loci_per_second is not None
                    else ""
                )
                progress_subtitle = (
                    f"Completed {alignment_progress.current_locus_name}; "
                    f"{alignment_progress.cached_locus_count:,} loci came "
                    f"from cache{progress_rate}; elapsed "
                    f"{alignment_progress.elapsed_seconds:.1f}s."
                )
            _alignment_progress_bar.update(
                increment=progress_increment,
                subtitle=progress_subtitle,
            )
            _completed_locus_count[0] = alignment_progress.completed_locus_count

        schema_alignment_analysis = get_schema_alignment_analysis(
            scheme_path,
            alignment_strategy=selected_alignment_strategy,
            progress_output_stream=mafft_progress_output_stream,
            progress_callback=_update_alignment_progress,
        )
    return (schema_alignment_analysis,)


@app.cell
def _(
    amino_acid_gap_warning_threshold_input,
    calculate_schema_alignment_summary,
    schema_alignment_analysis,
):
    amino_acid_gap_warning_threshold_percentage = float(
        amino_acid_gap_warning_threshold_input.value or 0
    )
    schema_alignment_summary = calculate_schema_alignment_summary(
        schema_alignment_analysis,
        amino_acid_gap_warning_threshold_percentage=(
            amino_acid_gap_warning_threshold_percentage
        ),
    )
    return (
        amino_acid_gap_warning_threshold_percentage,
        schema_alignment_summary,
    )


@app.cell(hide_code=True)
def _(
    amino_acid_gap_warning_threshold_percentage,
    alignment_strategy_labels,
    mo,
    schema_alignment_analysis,
    schema_alignment_summary,
):
    alignment_information_rows = [
        {
            "Locus": locus_statistics.locus_name,
            "Nucleotide Sequences Aligned": (
                locus_statistics.nucleotide_aligned_sequence_count
            ),
            "Nucleotide Sequences Excluded": (
                locus_statistics.nucleotide_excluded_sequence_count
            ),
            "Nucleotide Status": (locus_statistics.nucleotide_alignment_status),
            "Nucleotide Alignment Length": (
                locus_statistics.nucleotide_alignment_length
            ),
            "Nucleotide Gap Characters": (
                locus_statistics.nucleotide_gap_character_count
            ),
            "Nucleotide Gap Proportion (%)": (
                round(
                    locus_statistics.nucleotide_gap_proportion_percentage,
                    3,
                )
                if locus_statistics.nucleotide_gap_proportion_percentage is not None
                else None
            ),
            "Nucleotide Callable Sites": (
                locus_statistics.nucleotide_callable_site_count
            ),
            "Nucleotide Polymorphic Sites": (
                locus_statistics.nucleotide_polymorphic_site_count
            ),
            "Nucleotide Polymorphic Sites (%)": (
                round(
                    locus_statistics.nucleotide_polymorphic_site_percentage,
                    3,
                )
                if locus_statistics.nucleotide_polymorphic_site_percentage
                is not None
                else None
            ),
            "Amino-acid Sequences Aligned": (
                locus_statistics.amino_acid_aligned_sequence_count
            ),
            "Amino-acid Sequences Excluded": (
                locus_statistics.amino_acid_excluded_sequence_count
            ),
            "Amino-acid Status": (locus_statistics.amino_acid_alignment_status),
            "Amino-acid Alignment Length": (
                locus_statistics.amino_acid_alignment_length
            ),
            "Amino-acid Gap Characters": (
                locus_statistics.amino_acid_gap_character_count
            ),
            "Amino-acid Gap Proportion (%)": (
                round(
                    locus_statistics.amino_acid_gap_proportion_percentage,
                    3,
                )
                if locus_statistics.amino_acid_gap_proportion_percentage is not None
                else None
            ),
            "Above Advisory Amino-acid Threshold": (
                locus_statistics.amino_acid_gap_proportion_percentage
                > amino_acid_gap_warning_threshold_percentage
                if locus_statistics.amino_acid_gap_proportion_percentage is not None
                else None
            ),
            "Alignment Error": " | ".join(
                alignment_error
                for alignment_error in (
                    locus_statistics.nucleotide_alignment_error,
                    locus_statistics.amino_acid_alignment_error,
                )
                if alignment_error
            )
            or None,
        }
        for locus_statistics in schema_alignment_analysis.locus_statistics
    ]
    alignment_information_table = mo.ui.table(
        alignment_information_rows,
        pagination=True,
        page_size=50,
        selection=None,
        show_column_summaries=False,
        show_data_types=False,
        freeze_columns_left=["Locus"],
        text_justify_columns={
            column_name: "right"
            for column_name in alignment_information_rows[0]
            if column_name not in {"Locus", "Alignment Error"}
        }
        if alignment_information_rows
        else None,
        column_widths={
            "Locus": 280,
            "Nucleotide Sequences Aligned": 240,
            "Nucleotide Sequences Excluded": 245,
            "Nucleotide Status": 190,
            "Nucleotide Alignment Length": 240,
            "Nucleotide Gap Characters": 220,
            "Nucleotide Gap Proportion (%)": 245,
            "Nucleotide Callable Sites": 215,
            "Nucleotide Polymorphic Sites": 235,
            "Nucleotide Polymorphic Sites (%)": 255,
            "Amino-acid Sequences Aligned": 240,
            "Amino-acid Sequences Excluded": 245,
            "Amino-acid Status": 190,
            "Amino-acid Alignment Length": 235,
            "Amino-acid Gap Characters": 220,
            "Amino-acid Gap Proportion (%)": 240,
            "Above Advisory Amino-acid Threshold": 285,
            "Alignment Error": 360,
        },
        header_tooltip={
            "Nucleotide Gap Proportion (%)": (
                "All '-' characters divided by all cells in the nucleotide "
                "alignment; leading and trailing gaps are included."
            ),
            "Nucleotide Callable Sites": (
                "Columns with at least two canonical A/C/G/T observations; "
                "gaps and ambiguous symbols are treated as missing."
            ),
            "Nucleotide Polymorphic Sites": (
                "Callable columns containing at least two distinct canonical "
                "nucleotide states."
            ),
            "Nucleotide Polymorphic Sites (%)": (
                "Polymorphic nucleotide sites divided by callable nucleotide "
                "sites. Indel-only differences are excluded."
            ),
            "Amino-acid Gap Proportion (%)": (
                "All '-' characters divided by all cells in the amino-acid "
                "alignment; leading and trailing gaps are included."
            ),
            "Above Advisory Amino-acid Threshold": (
                "True only when the measured amino-acid gap proportion is "
                "strictly greater than the selected study-specific threshold."
            ),
        },
        show_search=True,
        show_download=True,
        max_columns=None,
        max_height=600,
    )
    alignment_analysis_displayed = True
    mo.vstack(
        [
            mo.md(
                f"""
                ## Alignment Results

                The analysis used **{alignment_strategy_labels[schema_alignment_analysis.alignment_strategy]}**
                (`{schema_alignment_analysis.alignment_method_identifier}`) with
                `{schema_alignment_analysis.mafft_version}`.

                **Availability:** nucleotide alignments were calculated for
                {schema_alignment_summary.nucleotide_aligned_locus_count:,} of
                {schema_alignment_summary.locus_count:,} loci;
                {schema_alignment_summary.nucleotide_unavailable_locus_count:,}
                had fewer than two eligible sequences and
                {schema_alignment_summary.nucleotide_failed_locus_count:,} failed.
                Amino-acid alignments were calculated for
                {schema_alignment_summary.amino_acid_aligned_locus_count:,} loci;
                {schema_alignment_summary.amino_acid_unavailable_locus_count:,}
                were unavailable and
                {schema_alignment_summary.amino_acid_failed_locus_count:,} failed.

                **Advisory amino-acid flag:**
                {schema_alignment_summary.amino_acid_gap_warning_locus_count:,}
                of {schema_alignment_summary.amino_acid_aligned_locus_count:,}
                aligned loci
                ({schema_alignment_summary.amino_acid_gap_warning_locus_percentage:.2f}%)
                are strictly above the selected threshold. Affected loci can be
                identified and downloaded from the table. The current advisory
                threshold is **>{amino_acid_gap_warning_threshold_percentage:.1f}%**.
                """
            ),
            alignment_information_table,
        ],
        gap=2,
    )
    return (alignment_analysis_displayed,)


@app.cell(hide_code=True)
def _(alignment_analysis_displayed, mo):
    similarity_section_ready = alignment_analysis_displayed
    mo.md(
        """
        ---

        # Cross-locus Sequence Similarity

        This analysis compares one deterministically selected allele
        **between different loci** at nucleotide and amino-acid level. It is
        a cross-locus similarity screen, but not the candidate-paralog
        classification from the next roadmap item. Similarity alone does not
        demonstrate copy number or evolutionary paralogy.

        **Why and how one allele is selected**

        - These schemes do not contain curated `short/` representative
          directories. The analysis therefore selects directly from each
          locus FASTA and does not label the chosen allele as a curated or
          biological representative.
        - The modal length is calculated from all non-empty alleles. If
          multiple lengths are equally frequent, the shortest tied length is
          used. The valid coding allele closest to that length is selected;
          original FASTA order breaks equal-distance ties.
        - If a locus has no valid coding allele, its eligible IUPAC DNA allele
          closest to the same modal length can still be used for the
          nucleotide comparison. Such a locus has no amino-acid comparison.
        - This one-selected-allele-per-locus screen is a fast approximation,
          not a sequence medoid calculation. The later all-alleles cross-locus
          roadmap item will provide broader, substantially more expensive
          coverage.
        - BLAST nucleotide and protein searches run concurrently and divide
          the CPUs available to the notebook. Results and provenance are
          cached with the source FASTA signature. The standard `blastn` and
          `blastp` tasks and their standard low-complexity masking behavior
          are used with an E-value floor of 0.001.
        - Identity is calculated within the selected local alignment. Query
          and subject coverage use their respective aligned coordinate spans;
          minimum coverage is the smaller of those values.
        - BSR is directional: `cross-locus raw score / query self score`. Each
          row therefore has a distinct query and subject, and reciprocal rows
          may differ. The best raw-score local alignment is retained without
          merging separate HSPs or clamping BSR values above 1.
        - Only BLAST-reported local hits are listed. No row means no hit passed
          the search E-value floor; it does not mean measured zero similarity.
        - These measurements do not depend on the **Minimum allele length**
          setting. No literature screening threshold is applied in this
          section.
        """
    )
    return (similarity_section_ready,)


@app.cell
def _(
    get_schema_similarity_analysis,
    mo,
    scheme_path,
    similarity_progress_output_stream,
    similarity_section_ready,
):
    _last_similarity_progress_count = [0]
    with mo.status.progress_bar(
        total=5,
        title="Measuring cross-locus nucleotide and amino-acid similarity",
        subtitle="Checking the versioned BLAST result cache.",
        completion_title="Cross-locus similarity analysis complete",
        completion_subtitle="Similarity measurements are ready.",
        show_rate=False,
        show_eta=False,
        remove_on_exit=False,
    ) as _similarity_progress_bar:

        def _update_similarity_progress(similarity_progress):
            _progress_increment = max(
                0,
                similarity_progress.completed_stage_count
                - _last_similarity_progress_count[0],
            )
            _similarity_progress_bar.update(
                increment=_progress_increment,
                subtitle=similarity_progress.message,
            )
            _last_similarity_progress_count[0] = (
                similarity_progress.completed_stage_count
            )

        schema_similarity_analysis = get_schema_similarity_analysis(
            scheme_path,
            progress_callback=_update_similarity_progress,
            progress_output_stream=similarity_progress_output_stream,
        )
    return (schema_similarity_analysis,)


@app.cell(hide_code=True)
def _(alt, mo, schema_similarity_analysis):
    similarity_information_rows = [
        {
            "Locus Pair": (
                f"{similarity_match.query_locus_name} → "
                f"{similarity_match.subject_locus_name}"
            ),
            "Sequence Level": (
                "Nucleotide"
                if similarity_match.sequence_level == "nucleotide"
                else "Amino acid"
            ),
            "Identity (%)": round(similarity_match.identity_percentage, 3),
            "Alignment Length": similarity_match.alignment_length,
            "Query Coverage (%)": round(
                similarity_match.query_coverage_percentage,
                3,
            ),
            "Subject Coverage (%)": round(
                similarity_match.subject_coverage_percentage,
                3,
            ),
            "Minimum Coverage (%)": round(
                similarity_match.minimum_coverage_percentage,
                3,
            ),
            "Raw Score": similarity_match.raw_score,
            "Directional BSR": round(similarity_match.blast_score_ratio, 6),
            "Query Allele ID": similarity_match.query_allele_identifier,
            "Subject Allele ID": similarity_match.subject_allele_identifier,
            "Subject Strand": similarity_match.nucleotide_subject_strand,
        }
        for similarity_match in schema_similarity_analysis.matches
    ]
    _maximum_similarity_plot_point_count = 20_000
    _all_similarity_plot_rows = [
        {
            "Sequence Level": (
                "Nucleotide"
                if similarity_match.sequence_level == "nucleotide"
                else "Amino acid"
            ),
            "Query Locus": similarity_match.query_locus_name,
            "Subject Locus": similarity_match.subject_locus_name,
            "Identity (%)": round(similarity_match.identity_percentage, 3),
            "Minimum Coverage (%)": round(
                similarity_match.minimum_coverage_percentage,
                3,
            ),
            "Directional BSR": round(
                similarity_match.blast_score_ratio,
                6,
            ),
        }
        for similarity_match in schema_similarity_analysis.matches
    ]
    if len(_all_similarity_plot_rows) > _maximum_similarity_plot_point_count:
        _plot_sampling_interval = (
            len(_all_similarity_plot_rows) / _maximum_similarity_plot_point_count
        )
        similarity_plot_rows = [
            _all_similarity_plot_rows[
                int(_plot_index * _plot_sampling_interval)
            ]
            for _plot_index in range(_maximum_similarity_plot_point_count)
        ]
    else:
        similarity_plot_rows = _all_similarity_plot_rows
    similarity_scatter_chart = alt.Chart(
        alt.Data(values=similarity_plot_rows)
    ).mark_circle(opacity=0.65, size=65).encode(
        x=alt.X(
            "Identity (%):Q",
            title="Identity (%)",
            scale=alt.Scale(domain=[0, 100]),
        ),
        y=alt.Y(
            "Minimum Coverage (%):Q",
            title="Minimum query/subject coverage (%)",
            scale=alt.Scale(domain=[0, 100]),
        ),
        color=alt.Color(
            "Sequence Level:N",
            title="Sequence level",
            scale=alt.Scale(
                domain=["Nucleotide", "Amino acid"],
                range=["#3b82f6", "#f97316"],
            ),
        ),
        tooltip=[
            alt.Tooltip("Sequence Level:N"),
            alt.Tooltip("Query Locus:N"),
            alt.Tooltip("Subject Locus:N"),
            alt.Tooltip("Identity (%):Q", format=".3f"),
            alt.Tooltip("Minimum Coverage (%):Q", format=".3f"),
            alt.Tooltip("Directional BSR:Q", format=".6f"),
        ],
    ).properties(
        title="Identity and coverage of BLAST-reported cross-locus hits",
        width=1100,
        height=420,
    )
    if similarity_plot_rows:
        similarity_chart = mo.ui.altair_chart(
            similarity_scatter_chart,
            chart_selection=False,
            legend_selection=False,
        )
    else:
        similarity_chart = mo.callout(
            "BLAST reported no cross-locus local similarities at the "
            "configured E-value floor.",
            kind="info",
        )
    if similarity_information_rows:
        similarity_information_table = mo.ui.table(
            similarity_information_rows,
            pagination=True,
            page_size=50,
            selection=None,
            show_column_summaries=False,
            show_data_types=False,
            freeze_columns_left=["Locus Pair"],
            text_justify_columns={
                column_name: "right"
                for column_name in (
                    "Identity (%)",
                    "Alignment Length",
                    "Query Coverage (%)",
                    "Subject Coverage (%)",
                    "Minimum Coverage (%)",
                    "Raw Score",
                    "Directional BSR",
                )
            },
            column_widths={
                "Locus Pair": 430,
                "Sequence Level": 180,
                "Identity (%)": 165,
                "Alignment Length": 190,
                "Query Coverage (%)": 205,
                "Subject Coverage (%)": 215,
                "Minimum Coverage (%)": 225,
                "Raw Score": 150,
                "Directional BSR": 190,
                "Query Allele ID": 240,
                "Subject Allele ID": 240,
                "Subject Strand": 180,
            },
            header_tooltip={
                "Locus Pair": (
                    "Directed query → subject locus pair. Direction matters "
                    "because BSR is normalized by the query self score."
                ),
                "Identity (%)": (
                    "Identical characters divided by the selected local HSP "
                    "alignment length."
                ),
                "Minimum Coverage (%)": (
                    "The smaller of the query and subject aligned-coordinate "
                    "coverage percentages."
                ),
                "Directional BSR": (
                    "Selected cross-locus raw score divided by the query "
                    "allele's self-alignment raw score."
                ),
                "Query Allele ID": (
                    "Original FASTA identifier of the selected query allele."
                ),
                "Subject Allele ID": (
                    "Original FASTA identifier of the selected subject allele."
                ),
            },
            show_search=True,
            show_download=True,
            max_columns=None,
            max_height=600,
        )
    else:
        similarity_information_table = mo.md(
            "No result rows are available for the searchable table."
        )
    mo.vstack(
        [
            mo.md(
                f"""
                ## Cross-locus Similarity Results

                **Selected alleles:** a nucleotide allele was available for
                {schema_similarity_analysis.nucleotide_selected_locus_count:,}
                of {schema_similarity_analysis.locus_count:,} loci; a valid
                coding allele and amino-acid translation were available for
                {schema_similarity_analysis.amino_acid_selected_locus_count:,}.

                **BLAST-reported directed hits:**
                {schema_similarity_analysis.nucleotide_match_count:,} nucleotide
                hits from
                {schema_similarity_analysis.nucleotide_possible_directed_pair_count:,}
                possible selected-allele directed pairs, and
                {schema_similarity_analysis.amino_acid_match_count:,} amino-acid
                hits from
                {schema_similarity_analysis.amino_acid_possible_directed_pair_count:,}
                possible selected-allele directed pairs.

                **Provenance:** `{schema_similarity_analysis.blastn_version}`;
                `{schema_similarity_analysis.blastp_version}`; E-value floor
                `{schema_similarity_analysis.blast_evalue:g}`; method
                `{schema_similarity_analysis.method_identifier}`. Cache reused:
                **{schema_similarity_analysis.reused_from_cache}**.

                The scatter plot displays
                **{len(similarity_plot_rows):,}** of
                **{len(similarity_information_rows):,}** result rows. When more
                than 20,000 hits are present, the plot uses a deterministic
                evenly spaced sample for browser performance; the searchable
                and downloadable table remains complete.
                """
            ),
            similarity_chart,
            similarity_information_table,
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
