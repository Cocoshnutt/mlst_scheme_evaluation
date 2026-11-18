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
    from schema import bsr as _bsr
    from schema import cache as _cache
    from schema import info as _info
    from schema import similarity as _similarity
    from schema import stats as _stats

    _alignment_source_file = mo.watch.file("schema/alignment.py")
    _allele_info_source_file = mo.watch.file("schema/allele_info.py")
    _bsr_source_file = mo.watch.file("schema/bsr.py")
    _cache_source_file = mo.watch.file("schema/cache.py")
    _info_source_file = mo.watch.file("schema/info.py")
    _similarity_source_file = mo.watch.file("schema/similarity.py")
    _stats_source_file = mo.watch.file("schema/stats.py")
    _alignment_source_file.read_bytes()
    _allele_info_source_file.read_bytes()
    _bsr_source_file.read_bytes()
    _cache_source_file.read_bytes()
    _info_source_file.read_bytes()
    _similarity_source_file.read_bytes()
    _stats_source_file.read_bytes()

    _cache = _importlib.reload(_cache)
    _info = _importlib.reload(_info)
    _allele_info = _importlib.reload(_allele_info)
    _bsr = _importlib.reload(_bsr)
    _stats = _importlib.reload(_stats)
    _alignment = _importlib.reload(_alignment)
    _similarity = _importlib.reload(_similarity)

    calculate_schema_alignment_summary = (
        _alignment.calculate_schema_alignment_summary
    )
    classify_representative_bsr_analysis = (
        _bsr.classify_representative_bsr_analysis
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
    get_schema_information = _info.get_schema_information
    get_schema_alignment_analysis = _alignment.get_schema_alignment_analysis
    get_schema_representative_bsr_analysis = (
        _bsr.get_schema_representative_bsr_analysis
    )
    get_schema_similarity_analysis = (
        _similarity.get_schema_similarity_analysis
    )
    classify_candidate_paralogous_locus_pairs = (
        _similarity.classify_candidate_paralogous_locus_pairs
    )
    CandidateParalogThresholds = _similarity.CandidateParalogThresholds
    default_candidate_paralog_amino_acid_identity_threshold_percentage = (
        _similarity.DEFAULT_CANDIDATE_PARALOG_AMINO_ACID_IDENTITY_THRESHOLD_PERCENTAGE
    )
    default_repeated_homology_nucleotide_identity_threshold_percentage = (
        _similarity.DEFAULT_REPEATED_HOMOLOGY_NUCLEOTIDE_IDENTITY_THRESHOLD_PERCENTAGE
    )
    default_repeated_homology_nucleotide_alignment_length_threshold = (
        _similarity.DEFAULT_REPEATED_HOMOLOGY_NUCLEOTIDE_ALIGNMENT_LENGTH_THRESHOLD
    )
    default_same_gene_blast_score_ratio_threshold = (
        _similarity.DEFAULT_SAME_GENE_BLAST_SCORE_RATIO_THRESHOLD
    )
    default_representative_bsr_threshold = (
        _bsr.DEFAULT_SAME_GENE_BLAST_SCORE_RATIO_THRESHOLD
    )
    default_frequent_peak_relative_threshold_percentage = (
        _stats.DEFAULT_FREQUENT_PEAK_RELATIVE_THRESHOLD_PERCENTAGE
    )
    get_schema_statistics = _stats.get_schema_statistics
    mafft_progress_output_stream = _sys.__stdout__ or _sys.stdout
    representative_bsr_progress_output_stream = _sys.__stdout__ or _sys.stdout
    similarity_progress_output_stream = _sys.__stdout__ or _sys.stdout
    analysis_cache_filenames = (
        _alignment.SCHEMA_ALIGNMENT_CACHE_FILENAME,
        _bsr.SCHEMA_REPRESENTATIVE_BSR_CACHE_FILENAME,
        _info.SCHEMA_INFORMATION_CACHE_FILENAME,
        _similarity.SCHEMA_SIMILARITY_CACHE_FILENAME,
        _stats.SCHEMA_STATISTICS_CACHE_FILENAME,
        _allele_info.SCHEME_ALLELE_INFORMATION_CACHE_FILENAME,
    )

    return (
        Path,
        CandidateParalogThresholds,
        alt,
        alignment_strategy_labels,
        analysis_cache_filenames,
        calculate_schema_alignment_summary,
        cache_archive,
        classify_candidate_paralogous_locus_pairs,
        classify_representative_bsr_analysis,
        clear_cache,
        count_alleles_below_minimum_length,
        create_schema_statistics_charts,
        default_candidate_paralog_amino_acid_identity_threshold_percentage,
        default_amino_acid_gap_warning_threshold_percentage,
        default_alignment_strategy,
        default_frequent_peak_relative_threshold_percentage,
        default_repeated_homology_nucleotide_alignment_length_threshold,
        default_repeated_homology_nucleotide_identity_threshold_percentage,
        default_representative_bsr_threshold,
        default_same_gene_blast_score_ratio_threshold,
        get_scheme_allele_information,
        get_schema_alignment_analysis,
        get_schema_information,
        get_schema_representative_bsr_analysis,
        get_schema_similarity_analysis,
        get_schema_statistics,
        mafft_progress_output_stream,
        mo,
        representative_bsr_progress_output_stream,
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
    - GC Content: Percentage of pooled canonical bases that are G or C. A, C, G, and T form the denominator; ambiguous and other noncanonical characters are excluded and reported separately.
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
    _guanine_cytosine_content_label = (
        "Unavailable (no canonical bases)"
        if schema_information.guanine_cytosine_content_percentage is None
        else f"{schema_information.guanine_cytosine_content_percentage:.2f}%"
    )
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
        - **GC content:** {_guanine_cytosine_content_label}
        - **Canonical bases in GC denominator:** {schema_information.canonical_nucleotide_base_count:,}
        - **Noncanonical bases excluded from GC denominator:** {schema_information.excluded_nucleotide_base_count:,}
        - **Alleles <{minimum_sequence_length:,} bp:** {alleles_below_minimum_length_count:,}
        - **Alleles <80% of locus mode:** {schema_information.alleles_below_locus_mode_threshold_count:,} ({schema_information.alleles_below_locus_mode_threshold_percentage:.2f}%)
        - **Alleles >120% of locus mode:** {schema_information.alleles_above_locus_mode_threshold_count:,} ({schema_information.alleles_above_locus_mode_threshold_percentage:.2f}%)
        """
    )
    return (schema_summary_displayed,)


@app.cell(hide_code=True)
def _(
    default_frequent_peak_relative_threshold_percentage,
    mo,
    schema_summary_displayed,
):
    frequent_peak_relative_threshold_input = mo.ui.number(
        start=0,
        stop=100,
        step=0.1,
        value=default_frequent_peak_relative_threshold_percentage,
        debounce=True,
        label="Frequent peak minimum relative height (%)",
    )
    locus_statistics_section_ready = schema_summary_displayed
    mo.vstack(
        [
            mo.md(
                """
                ---

                # Locus Statistics

                These panels summarize allele counts and sequence-length
                variation across the evaluated loci.

                **Length-distribution definitions**

                - Mean length is the arithmetic mean across every FASTA record
                  in the locus. Existing zero-length record handling is
                  preserved.
                - All exact lengths tied for the highest record count are
                  reported as modes. The shortest tied mode remains the
                  primary mode used by the existing charts and outlier rule.
                - A frequent length peak is an exact length whose record count
                  is at least the selected percentage of the dominant exact-
                  length count. The default 50% value is a configurable
                  descriptive project heuristic, not a literature-derived
                  quality threshold.
                - A locus is single-mode when one frequent peak remains and
                  multimodal when more than one remains. An outlier suffix is
                  added when any allele is strictly below 80% or strictly
                  above 120% of the primary mode. Values exactly on either
                  boundary are not outliers. A zero-record locus is Empty.
                - Changing the frequent-peak percentage reclassifies cached
                  exact-length counts without rereading FASTA files. It is
                  independent of the **Minimum allele length** control.

                **GC-content definition**

                - Per-locus GC content pools every allele record in that locus,
                  including duplicate records, and calculates
                  `(G + C) / (A + C + G + T) × 100`.
                - Ambiguous and other noncanonical characters are excluded from
                  the denominator and counted separately. A locus without any
                  canonical bases is reported as unavailable rather than 0%.
                - The scheme-wide value in the top summary pools canonical
                  bases across the whole scheme; it is not an unweighted mean
                  of locus percentages. GC content is independent of the
                  **Minimum allele length** control.
                """
            ),
            frequent_peak_relative_threshold_input,
        ]
    )
    return frequent_peak_relative_threshold_input, locus_statistics_section_ready


@app.cell
def _(
    create_schema_statistics_charts,
    frequent_peak_relative_threshold_input,
    get_schema_statistics,
    locus_statistics_section_ready,
    mo,
    scheme_path,
):
    with mo.status.spinner(
        title="Calculating locus statistics",
        subtitle="Loading per-locus allele-count and length distributions.",
    ):
        frequent_peak_relative_threshold_percentage = float(
            frequent_peak_relative_threshold_input.value
        )
        locus_statistics = get_schema_statistics(
            scheme_path,
            frequent_peak_relative_threshold_percentage=(
                frequent_peak_relative_threshold_percentage
            ),
        )
        schema_statistics_charts = create_schema_statistics_charts(locus_statistics)
    return (
        frequent_peak_relative_threshold_percentage,
        locus_statistics,
        schema_statistics_charts,
    )


@app.cell(hide_code=True)
def _(
    frequent_peak_relative_threshold_percentage,
    locus_statistics,
    mo,
    schema_statistics_charts,
):
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
    guanine_cytosine_content_chart = mo.ui.altair_chart(
        schema_statistics_charts.guanine_cytosine_content_chart,
        chart_selection=False,
        legend_selection=False,
    )
    _length_distribution_categories = (
        "Empty",
        "Single-mode",
        "Single-mode with outliers",
        "Multimodal",
        "Multimodal with outliers",
    )
    _category_counts = {
        category: sum(
            locus_information.length_distribution_category == category
            for locus_information in locus_statistics
        )
        for category in _length_distribution_categories
    }
    _maximum_length_value_preview_count = 10

    def _format_length_value_preview(formatted_values):
        _visible_values = formatted_values[
            :_maximum_length_value_preview_count
        ]
        _hidden_value_count = len(formatted_values) - len(_visible_values)
        _preview = "; ".join(_visible_values)
        if _hidden_value_count > 0:
            _preview = f"{_preview}; … (+{_hidden_value_count:,} more)"
        return _preview or None

    _locus_length_distribution_rows = [
        {
            "Locus": locus_information.locus_name,
            "Alleles": locus_information.allele_count,
            "Mean Length (bp)": (
                round(locus_information.mean_allele_length, 3)
                if locus_information.mean_allele_length is not None
                else None
            ),
            "Primary Mode (bp)": locus_information.modal_allele_length,
            "Tied Mode Count": len(locus_information.modal_allele_lengths),
            "Tied Modes Preview (bp)": _format_length_value_preview(
                [
                    f"{allele_length:,}"
                    for allele_length in locus_information.modal_allele_lengths
                ]
            ),
            "Frequent Peak Count": len(
                locus_information.frequent_allele_length_peaks
            ),
            "Frequent Peaks Preview": _format_length_value_preview(
                [
                    f"{peak.allele_length:,} bp: {peak.allele_count:,} "
                    f"({peak.allele_percentage:.2f}%)"
                    for peak in sorted(
                        locus_information.frequent_allele_length_peaks,
                        key=lambda peak: (
                            -peak.allele_count,
                            peak.allele_length,
                        ),
                    )
                ]
            ),
            "Alleles <80% Mode": (
                locus_information.alleles_below_mode_threshold_count
            ),
            "Alleles >120% Mode": (
                locus_information.alleles_above_mode_threshold_count
            ),
            "Length Distribution Category": (
                locus_information.length_distribution_category
            ),
            "GC Content (%)": (
                round(
                    locus_information.guanine_cytosine_content_percentage,
                    3,
                )
                if locus_information.guanine_cytosine_content_percentage
                is not None
                else None
            ),
            "Canonical Bases": (
                locus_information.canonical_nucleotide_base_count
            ),
            "Excluded Noncanonical Bases": (
                locus_information.excluded_nucleotide_base_count
            ),
        }
        for locus_information in locus_statistics
    ]
    if _locus_length_distribution_rows:
        locus_length_distribution_table = mo.ui.table(
            _locus_length_distribution_rows,
            pagination=True,
            page_size=50,
            selection=None,
            show_column_summaries=False,
            show_data_types=False,
            freeze_columns_left=["Locus"],
            text_justify_columns={
                column_name: "right"
                for column_name in (
                    "Alleles",
                    "Mean Length (bp)",
                    "Primary Mode (bp)",
                    "Tied Mode Count",
                    "Frequent Peak Count",
                    "Alleles <80% Mode",
                    "Alleles >120% Mode",
                    "GC Content (%)",
                    "Canonical Bases",
                    "Excluded Noncanonical Bases",
                )
            },
            column_widths={
                "Locus": 430,
                "Alleles": 140,
                "Mean Length (bp)": 190,
                "Primary Mode (bp)": 190,
                "Tied Mode Count": 180,
                "Tied Modes Preview (bp)": 320,
                "Frequent Peak Count": 190,
                "Frequent Peaks Preview": 480,
                "Alleles <80% Mode": 210,
                "Alleles >120% Mode": 220,
                "Length Distribution Category": 285,
                "GC Content (%)": 180,
                "Canonical Bases": 190,
                "Excluded Noncanonical Bases": 290,
            },
            header_tooltip={
                "Mean Length (bp)": (
                    "Arithmetic mean across every allele record in the locus."
                ),
                "Primary Mode (bp)": (
                    "Shortest length among all lengths tied for the highest "
                    "frequency."
                ),
                "Tied Mode Count": (
                    "Number of exact allele lengths tied for the highest "
                    "frequency."
                ),
                "Tied Modes Preview (bp)": (
                    "Up to ten tied modes in ascending order. The full count "
                    "is reported separately."
                ),
                "Frequent Peak Count": (
                    "Number of exact lengths meeting the selected percentage "
                    "of the dominant length count."
                ),
                "Frequent Peaks Preview": (
                    "Exact lengths meeting the selected percentage of the "
                    "dominant length count, ordered by count then length. Up "
                    "to ten values are shown."
                ),
                "Length Distribution Category": (
                    "Single-mode or multimodal under the selected peak rule, "
                    "plus whether strict 80/120% modal-length outliers exist."
                ),
                "GC Content (%)": (
                    "Pooled (G + C) / (A + C + G + T) across every allele "
                    "record in the locus."
                ),
                "Canonical Bases": (
                    "A/C/G/T bases used as the locus GC-content denominator."
                ),
                "Excluded Noncanonical Bases": (
                    "Non-whitespace characters excluded from the GC-content "
                    "denominator because they are not A, C, G, or T."
                ),
            },
            show_search=True,
            show_download=True,
            max_columns=None,
            max_height=600,
        )
    else:
        locus_length_distribution_table = mo.callout(
            "No locus FASTA files are available for the length-distribution "
            "table.",
            kind="info",
        )
    locus_statistics_displayed = True
    mo.vstack(
        [
            mo.md(
                f"""
                **{len(locus_statistics):,} loci evaluated.** With frequent
                peaks defined at
                **≥{frequent_peak_relative_threshold_percentage:g}% of the
                dominant exact-length count**:

                - **Single-mode:** {_category_counts["Single-mode"]:,}
                - **Single-mode with outliers:** {_category_counts["Single-mode with outliers"]:,}
                - **Multimodal:** {_category_counts["Multimodal"]:,}
                - **Multimodal with outliers:** {_category_counts["Multimodal with outliers"]:,}
                - **Empty:** {_category_counts["Empty"]:,}
                """
            ),
            total_alleles_chart,
            allele_mode_size_chart,
            locus_statistics_chart,
            allele_size_variation_chart,
            mo.md(
                """
                ## Per-locus GC-content Distribution

                The histogram uses the pooled canonical-base GC percentage for
                each locus. Loci without canonical bases have no percentage and
                are omitted from the histogram, but remain visible in the table
                below with their canonical and excluded base counts.
                """
            ),
            guanine_cytosine_content_chart,
            mo.md(
                """
                ## Per-locus Length-distribution and GC Statistics

                To keep large schemes responsive, tied-mode and frequent-peak
                cells show at most ten values plus the number omitted. Their
                complete counts are reported in separate columns; the typed
                library results retain every value.
                """
            ),
            locus_length_distribution_table,
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

        Per-locus allele integrity and diversity metrics derived from the
        cached allele-record analysis. Each FASTA record is treated as a
        distinct nucleotide allele because the current scheme collection is
        expected to contain distinct allele records. Valid CDS alleles are
        translated and collapsed by protein sequence, so synonymous nucleotide
        alleles can contribute one distinct translated protein. Invalid or
        untranslatable alleles do not contribute a protein. These counts are
        independent of the **Minimum allele length** control.
        """
    )
    return (allele_analysis_section_ready,)


@app.cell
def _(
    allele_analysis_section_ready,
    get_scheme_allele_information,
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
    return (scheme_allele_information,)


@app.cell(hide_code=True)
def _(
    minimum_sequence_length,
    mo,
    scheme_allele_information,
):
    minimum_length_column_name = f"Alleles <{minimum_sequence_length}bp"
    allele_information_rows = [
        {
            "Locus": locus_information.locus_name,
            "Nucleotide Alleles": locus_information.total_allele_count,
            "Valid Alleles": locus_information.valid_allele_count,
            "Invalid Alleles": locus_information.invalid_allele_count,
            "Proportion of Validated Alleles": round(
                locus_information.validated_allele_proportion,
                3,
            ),
            "Distinct Translated Proteins": (
                locus_information.distinct_protein_allele_count
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
            "Nucleotide Alleles": 190,
            "Valid Alleles": 130,
            "Invalid Alleles": 145,
            "Proportion of Validated Alleles": 245,
            "Distinct Translated Proteins": 255,
            "Incomplete ORF": 150,
            "Ambiguous Bases": 165,
            "Missing Start/Stop Codon": 215,
            "In-frame Stop Codon": 185,
            minimum_length_column_name: 155,
            "Alleles below threshold": 200,
            "Alleles above threshold": 200,
        },
        header_tooltip={
            "Invalid Alleles": (
                "Sum of Incomplete ORF, Ambiguous Bases, Missing Start/Stop "
                "Codon, and In-frame Stop Codon counts."
            ),
            "Nucleotide Alleles": (
                "FASTA records in this locus. Records are treated as distinct "
                "nucleotide alleles; no sequence-deduplication criterion is "
                "applied."
            ),
            "Distinct Translated Proteins": (
                "Unique protein sequences translated from alleles that pass "
                "the existing CDS-validity rules. Synonymous nucleotide "
                "alleles can encode the same protein."
            ),
            "Alleles below threshold": (
                "Alleles strictly below 80% of the locus modal length."
            ),
            "Alleles above threshold": (
                "Alleles strictly above 120% of the locus modal length."
            ),
        },
        show_search=True,
        show_download=True,
        max_columns=None,
        max_height=600,
    )
    allele_analysis_displayed = True
    mo.vstack([allele_information_table])
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
        the measurement stage used by the candidate-paralog screen below.
        Similarity alone does not demonstrate copy number or evolutionary
        paralogy.

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
          setting. This raw-result section remains unfiltered; the next
          section applies configurable literature-derived screening values.
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
    similarity_analysis_displayed = True
    return (similarity_analysis_displayed,)


@app.cell(hide_code=True)
def _(
    default_candidate_paralog_amino_acid_identity_threshold_percentage,
    default_repeated_homology_nucleotide_alignment_length_threshold,
    default_repeated_homology_nucleotide_identity_threshold_percentage,
    default_same_gene_blast_score_ratio_threshold,
    mo,
    similarity_analysis_displayed,
):
    candidate_amino_acid_identity_threshold_input = mo.ui.number(
        start=0,
        stop=100,
        step=0.1,
        value=(
            default_candidate_paralog_amino_acid_identity_threshold_percentage
        ),
        debounce=True,
        label="Amino-acid identity threshold (%)",
    )
    candidate_nucleotide_identity_threshold_input = mo.ui.number(
        start=0,
        stop=100,
        step=0.1,
        value=(
            default_repeated_homology_nucleotide_identity_threshold_percentage
        ),
        debounce=True,
        label="Nucleotide identity threshold (%)",
    )
    candidate_nucleotide_alignment_length_threshold_input = mo.ui.number(
        start=0,
        step=1,
        value=(
            default_repeated_homology_nucleotide_alignment_length_threshold
        ),
        debounce=True,
        label="Nucleotide alignment-length threshold (bp)",
    )
    candidate_blast_score_ratio_threshold_input = mo.ui.number(
        start=0,
        step=0.01,
        value=default_same_gene_blast_score_ratio_threshold,
        debounce=True,
        label="Supplementary directional BSR threshold",
    )
    candidate_paralog_section_ready = similarity_analysis_displayed
    mo.vstack(
        [
            mo.md(
                """
                ---

                # Candidate Paralogous Loci

                Reciprocal BLAST rows are collapsed into one unordered locus
                pair. A pair is flagged when it has **amino-acid identity
                strictly greater than the selected threshold**, or
                **nucleotide identity at least the selected threshold across
                an alignment strictly longer than the selected base-pair
                threshold**. No additional coverage gate is applied; coverage
                remains visible for interpretation. For every flagged pair,
                the table shows the strongest available measurement at each
                sequence level. A measurement can therefore be shown even
                when that level did not classify the pair; a blank means
                BLAST reported no hit at that sequence level.

                These defaults reproduce study-specific screening references:
                greater than 70% amino-acid identity from Guglielmini et al.
                (2021), and at least 90% nucleotide identity across more than
                100 bp from Pires et al. (2024). Directional BSR at least 0.6,
                from Silva et al. (2018), is shown only as supplementary
                same-gene evidence and **cannot flag a candidate by itself**.

                This is an intrinsic selected-allele similarity screen. It
                reports **candidate** paralogous loci, not confirmed paralogy
                or copy number. Confirmation requires genome-level evidence,
                and the future all-alleles cross-locus analysis will provide
                broader sequence coverage. Threshold changes reclassify the
                cached BLAST measurements without rerunning BLAST and are
                independent of the **Minimum allele length** setting.
                """
            ),
            candidate_amino_acid_identity_threshold_input,
            candidate_nucleotide_identity_threshold_input,
            candidate_nucleotide_alignment_length_threshold_input,
            candidate_blast_score_ratio_threshold_input,
        ],
        gap=1,
    )
    return (
        candidate_amino_acid_identity_threshold_input,
        candidate_blast_score_ratio_threshold_input,
        candidate_nucleotide_alignment_length_threshold_input,
        candidate_nucleotide_identity_threshold_input,
        candidate_paralog_section_ready,
    )


@app.cell
def _(
    CandidateParalogThresholds,
    candidate_amino_acid_identity_threshold_input,
    candidate_blast_score_ratio_threshold_input,
    candidate_nucleotide_alignment_length_threshold_input,
    candidate_nucleotide_identity_threshold_input,
    candidate_paralog_section_ready,
    classify_candidate_paralogous_locus_pairs,
    schema_similarity_analysis,
):
    candidate_paralog_thresholds = CandidateParalogThresholds(
        amino_acid_identity_threshold_percentage=float(
            candidate_amino_acid_identity_threshold_input.value
        ),
        nucleotide_identity_threshold_percentage=float(
            candidate_nucleotide_identity_threshold_input.value
        ),
        nucleotide_alignment_length_threshold=int(
            candidate_nucleotide_alignment_length_threshold_input.value
        ),
        blast_score_ratio_threshold=float(
            candidate_blast_score_ratio_threshold_input.value
        ),
    )
    candidate_paralog_analysis = classify_candidate_paralogous_locus_pairs(
        schema_similarity_analysis,
        thresholds=candidate_paralog_thresholds,
    )
    return candidate_paralog_analysis, candidate_paralog_thresholds


@app.cell(hide_code=True)
def _(candidate_paralog_analysis, candidate_paralog_thresholds, mo):
    candidate_paralog_information_rows = []
    for candidate_pair in candidate_paralog_analysis.candidate_pairs:
        _amino_acid_match = candidate_pair.amino_acid_evidence_match
        _nucleotide_match = candidate_pair.nucleotide_evidence_match
        _maximum_bsr_match = candidate_pair.maximum_blast_score_ratio_match
        candidate_paralog_information_rows.append(
            {
                "Locus Pair": (
                    f"{candidate_pair.first_locus_name} ↔ "
                    f"{candidate_pair.second_locus_name}"
                ),
                "Classifying Evidence": candidate_pair.evidence_label,
                "Amino-acid Identity (%)": (
                    round(_amino_acid_match.identity_percentage, 3)
                    if _amino_acid_match is not None
                    else None
                ),
                "Amino-acid Alignment Length": (
                    _amino_acid_match.alignment_length
                    if _amino_acid_match is not None
                    else None
                ),
                "Amino-acid Minimum Coverage (%)": (
                    round(_amino_acid_match.minimum_coverage_percentage, 3)
                    if _amino_acid_match is not None
                    else None
                ),
                "Amino-acid Evidence Direction": (
                    f"{_amino_acid_match.query_locus_name} → "
                    f"{_amino_acid_match.subject_locus_name}"
                    if _amino_acid_match is not None
                    else None
                ),
                "Nucleotide Identity (%)": (
                    round(_nucleotide_match.identity_percentage, 3)
                    if _nucleotide_match is not None
                    else None
                ),
                "Nucleotide Alignment Length (bp)": (
                    _nucleotide_match.alignment_length
                    if _nucleotide_match is not None
                    else None
                ),
                "Nucleotide Minimum Coverage (%)": (
                    round(_nucleotide_match.minimum_coverage_percentage, 3)
                    if _nucleotide_match is not None
                    else None
                ),
                "Nucleotide Evidence Direction": (
                    f"{_nucleotide_match.query_locus_name} → "
                    f"{_nucleotide_match.subject_locus_name}"
                    if _nucleotide_match is not None
                    else None
                ),
                "Maximum Directional BSR": round(
                    _maximum_bsr_match.blast_score_ratio,
                    6,
                ),
                "Maximum BSR Direction": (
                    f"{_maximum_bsr_match.query_locus_name} → "
                    f"{_maximum_bsr_match.subject_locus_name}"
                ),
                "BSR Reference Reached": (
                    candidate_pair.meets_blast_score_ratio_reference
                ),
            }
        )

    if candidate_paralog_information_rows:
        candidate_paralog_information_table = mo.ui.table(
            candidate_paralog_information_rows,
            pagination=True,
            page_size=50,
            selection=None,
            show_column_summaries=False,
            show_data_types=False,
            freeze_columns_left=["Locus Pair"],
            text_justify_columns={
                column_name: "right"
                for column_name in (
                    "Amino-acid Identity (%)",
                    "Amino-acid Alignment Length",
                    "Amino-acid Minimum Coverage (%)",
                    "Nucleotide Identity (%)",
                    "Nucleotide Alignment Length (bp)",
                    "Nucleotide Minimum Coverage (%)",
                    "Maximum Directional BSR",
                )
            },
            column_widths={
                "Locus Pair": 430,
                "Classifying Evidence": 230,
                "Amino-acid Identity (%)": 220,
                "Amino-acid Alignment Length": 245,
                "Amino-acid Minimum Coverage (%)": 285,
                "Amino-acid Evidence Direction": 340,
                "Nucleotide Identity (%)": 220,
                "Nucleotide Alignment Length (bp)": 285,
                "Nucleotide Minimum Coverage (%)": 285,
                "Nucleotide Evidence Direction": 340,
                "Maximum Directional BSR": 235,
                "Maximum BSR Direction": 340,
                "BSR Reference Reached": 220,
            },
            header_tooltip={
                "Locus Pair": (
                    "One unordered pair after collapsing reciprocal directed "
                    "BLAST rows."
                ),
                "Classifying Evidence": (
                    "The amino-acid identity and/or nucleotide repeated-"
                    "homology criterion that flags this pair."
                ),
                "Amino-acid Identity (%)": (
                    "Strongest qualifying amino-acid identity when that "
                    "criterion is met; otherwise the strongest reported "
                    "amino-acid match."
                ),
                "Amino-acid Alignment Length": (
                    "Length of the displayed protein HSP in amino acids."
                ),
                "Nucleotide Identity (%)": (
                    "Strongest qualifying nucleotide identity when that "
                    "criterion is met; otherwise the strongest reported "
                    "nucleotide match."
                ),
                "Nucleotide Alignment Length (bp)": (
                    "Length of the displayed nucleotide HSP in base pairs."
                ),
                "Maximum Directional BSR": (
                    "Highest BSR across both directions and sequence levels; "
                    "supplementary evidence only."
                ),
            },
            show_search=True,
            show_download=True,
            max_columns=None,
            max_height=600,
        )
    else:
        candidate_paralog_information_table = mo.callout(
            "No BLAST-reported locus pair meets either selected candidate "
            "screening criterion.",
            kind="info",
        )

    candidate_paralog_analysis_displayed = True
    mo.vstack(
        [
            mo.md(
                f"""
                ## Candidate-paralog Summary Statistics

                **{candidate_paralog_analysis.candidate_pair_count:,} unique
                candidate pairs** involve
                **{candidate_paralog_analysis.candidate_locus_count:,} unique
                loci**. Of these pairs,
                **{candidate_paralog_analysis.amino_acid_candidate_pair_count:,}**
                meet the amino-acid criterion,
                **{candidate_paralog_analysis.nucleotide_candidate_pair_count:,}**
                meet the nucleotide criterion, and
                **{candidate_paralog_analysis.both_criteria_candidate_pair_count:,}**
                meet both.

                **Supplementary BSR evidence:**
                **{candidate_paralog_analysis.blast_score_ratio_evidence_pair_count:,}**
                of
                **{candidate_paralog_analysis.evaluated_locus_pair_count:,}**
                unique BLAST-reported locus pairs reach directional BSR
                ≥ {candidate_paralog_thresholds.blast_score_ratio_threshold:g}.
                A BSR-only pair is not included in the candidate table.

                **Applied screening values:** amino-acid identity
                > {candidate_paralog_thresholds.amino_acid_identity_threshold_percentage:g}%;
                or nucleotide identity
                ≥ {candidate_paralog_thresholds.nucleotide_identity_threshold_percentage:g}%
                across >
                {candidate_paralog_thresholds.nucleotide_alignment_length_threshold:,}
                bp. These are configurable literature references, not a
                scheme-level pass/fail classification.
                """
            ),
            candidate_paralog_information_table,
        ],
        gap=2,
    )
    return (candidate_paralog_analysis_displayed,)


@app.cell(hide_code=True)
def _(
    candidate_paralog_analysis_displayed,
    default_representative_bsr_threshold,
    mo,
):
    representative_bsr_threshold_input = mo.ui.number(
        start=0,
        step=0.01,
        value=default_representative_bsr_threshold,
        debounce=True,
        label="Representative-to-allele BSR threshold",
    )
    representative_bsr_section_ready = candidate_paralog_analysis_displayed
    mo.vstack(
        [
            mo.md(
                """
                ---

                # Representative-versus-allele BSR Coverage

                This protein-level analysis asks whether the valid translated
                alleles in each locus are represented by a sufficiently
                similar selected allele, and whether that selected allele also
                has high-BSR matches to alleles assigned to other loci.

                **How and why it is calculated**

                - The selected schemes do not contain curated `short/`
                  representative directories. One **approximate
                  representative** is therefore selected from each locus using
                  the same reproducible rule as the earlier similarity
                  analysis: the valid coding allele closest to the shortest
                  tied modal nucleotide length, with FASTA order breaking
                  equal-distance ties.
                - Every allele that passes the repository's CDS-validity rules
                  is translated. Invalid or untranslatable records are shown
                  as unevaluable and are excluded from the protein coverage
                  denominator.
                - The selected representatives are searched against all valid
                  translated alleles with BLASTp. Directional BSR is
                  `representative-to-allele raw score / representative self
                  raw score`. A missing BLAST hit is kept separate from a
                  reported hit below the threshold.
                - Coverage uses the inclusive comparison **BSR ≥ threshold**.
                  The default 0.6 value is a study-specific same-gene
                  similarity reference described by Silva et al. (2018), not
                  a universal scheme pass/fail threshold.
                - Cross-locus hits meeting the threshold are collapsed into
                  unordered locus pairs and then into transitive connected
                  groups. These are **high-BSR similarity groups**, not
                  confirmed paralogs. Copy number and evolutionary paralogy
                  require genome-level evidence.
                - BLASTp uses the CPUs available to the notebook. Raw BSR
                  measurements are cached, so changing this threshold only
                  reclassifies existing values. The analysis is independent
                  of the **Minimum allele length** setting.
                """
            ),
            representative_bsr_threshold_input,
        ],
        gap=1,
    )
    return representative_bsr_section_ready, representative_bsr_threshold_input


@app.cell
def _(
    get_schema_representative_bsr_analysis,
    mo,
    representative_bsr_progress_output_stream,
    representative_bsr_section_ready,
    scheme_path,
):
    _last_representative_bsr_progress_count = [0]
    with mo.status.progress_bar(
        total=5,
        title="Measuring representative-versus-allele protein BSR",
        subtitle="Checking the versioned representative BSR cache.",
        completion_title="Representative BSR analysis complete",
        completion_subtitle="Coverage and high-BSR groups are ready.",
        show_rate=False,
        show_eta=False,
        remove_on_exit=False,
    ) as _representative_bsr_progress_bar:

        def _update_representative_bsr_progress(representative_bsr_progress):
            _progress_increment = max(
                0,
                representative_bsr_progress.completed_stage_count
                - _last_representative_bsr_progress_count[0],
            )
            _representative_bsr_progress_bar.update(
                increment=_progress_increment,
                subtitle=representative_bsr_progress.message,
            )
            _last_representative_bsr_progress_count[0] = (
                representative_bsr_progress.completed_stage_count
            )

        schema_representative_bsr_analysis = (
            get_schema_representative_bsr_analysis(
                scheme_path,
                progress_callback=_update_representative_bsr_progress,
                progress_output_stream=representative_bsr_progress_output_stream,
            )
        )
    return (schema_representative_bsr_analysis,)


@app.cell
def _(
    classify_representative_bsr_analysis,
    representative_bsr_threshold_input,
    schema_representative_bsr_analysis,
):
    representative_bsr_classification = classify_representative_bsr_analysis(
        schema_representative_bsr_analysis,
        blast_score_ratio_threshold=float(
            representative_bsr_threshold_input.value
        ),
    )
    return (representative_bsr_classification,)


@app.cell(hide_code=True)
def _(
    mo,
    representative_bsr_classification,
    schema_representative_bsr_analysis,
):
    _maximum_locus_preview_count = 2_000
    _maximum_pair_preview_count = 1_000
    _maximum_group_member_preview_count = 25

    _all_locus_coverage_rows = [
        {
            "Locus": locus_coverage.locus_name,
            "Total Alleles": locus_coverage.total_allele_count,
            "Valid Translated Alleles": (
                locus_coverage.valid_translated_allele_count
            ),
            "Unevaluable Alleles": locus_coverage.unevaluable_allele_count,
            "Covered Alleles": locus_coverage.covered_allele_count,
            "Reported Hits Below Threshold": (
                locus_coverage.below_threshold_hit_count
            ),
            "Representative Coverage (%)": (
                round(locus_coverage.coverage_percentage, 3)
                if locus_coverage.coverage_percentage is not None
                else None
            ),
        }
        for locus_coverage in representative_bsr_classification.locus_coverages
    ]
    representative_bsr_locus_rows = _all_locus_coverage_rows[
        :_maximum_locus_preview_count
    ]
    if representative_bsr_locus_rows:
        representative_bsr_locus_table = mo.ui.table(
            representative_bsr_locus_rows,
            pagination=True,
            page_size=50,
            selection=None,
            show_column_summaries=False,
            show_data_types=False,
            freeze_columns_left=["Locus"],
            text_justify_columns={
                _column_name: "right"
                for _column_name in representative_bsr_locus_rows[0]
                if _column_name != "Locus"
            },
            column_widths={
                "Locus": 300,
                "Total Alleles": 165,
                "Valid Translated Alleles": 235,
                "Unevaluable Alleles": 205,
                "Covered Alleles": 180,
                "Reported Hits Below Threshold": 285,
                "Representative Coverage (%)": 285,
            },
            header_tooltip={
                "Valid Translated Alleles": (
                    "Alleles passing the CDS-validity rules and included in "
                    "the protein-level coverage denominator."
                ),
                "Unevaluable Alleles": (
                    "Allele records excluded from protein BSR because they "
                    "are not valid translatable CDS sequences."
                ),
                "Representative Coverage (%)": (
                    "Valid translated alleles with representative-to-allele "
                    "BSR at least the selected threshold."
                ),
            },
            show_search=True,
            show_download=False,
            max_columns=None,
            max_height=600,
        )
    else:
        representative_bsr_locus_table = mo.callout(
            "No locus FASTA files are available for BSR coverage.",
            kind="info",
        )

    _sorted_high_bsr_pairs = sorted(
        representative_bsr_classification.high_bsr_pairs,
        key=lambda _pair: (
            -_pair.evidence.blast_score_ratio,
            _pair.first_locus_name,
            _pair.second_locus_name,
        ),
    )
    representative_bsr_pair_rows = [
        {
            "Locus Pair": (
                f"{high_bsr_pair.first_locus_name} ↔ "
                f"{high_bsr_pair.second_locus_name}"
            ),
            "Maximum Directional BSR": round(
                high_bsr_pair.evidence.blast_score_ratio,
                6,
            ),
            "Evidence Direction": (
                f"{high_bsr_pair.evidence.query_locus_name} → "
                f"{high_bsr_pair.evidence.subject_locus_name}"
            ),
            "Supporting Subject Allele": (
                high_bsr_pair.evidence.subject_allele_identifier
            ),
            "Raw Score": high_bsr_pair.evidence.raw_score,
        }
        for high_bsr_pair in _sorted_high_bsr_pairs[
            :_maximum_pair_preview_count
        ]
    ]
    if representative_bsr_pair_rows:
        representative_bsr_pair_table = mo.ui.table(
            representative_bsr_pair_rows,
            pagination=True,
            page_size=50,
            selection=None,
            show_column_summaries=False,
            show_data_types=False,
            freeze_columns_left=["Locus Pair"],
            text_justify_columns={
                "Maximum Directional BSR": "right",
                "Raw Score": "right",
            },
            column_widths={
                "Locus Pair": 430,
                "Maximum Directional BSR": 245,
                "Evidence Direction": 360,
                "Supporting Subject Allele": 285,
                "Raw Score": 150,
            },
            header_tooltip={
                "Locus Pair": (
                    "One unordered pair after collapsing reciprocal "
                    "representative-to-allele directions."
                ),
                "Maximum Directional BSR": (
                    "Strongest representative-normalized protein raw score "
                    "observed for this locus pair."
                ),
                "Supporting Subject Allele": (
                    "Allele from the other locus supporting the displayed "
                    "maximum directional BSR."
                ),
            },
            show_search=True,
            show_download=False,
            max_columns=None,
            max_height=600,
        )
    else:
        representative_bsr_pair_table = mo.callout(
            "No inter-locus representative-to-allele match reaches the "
            "selected BSR threshold.",
            kind="info",
        )

    representative_bsr_group_rows = []
    for high_bsr_group in representative_bsr_classification.high_bsr_groups:
        _visible_group_loci = high_bsr_group.locus_names[
            :_maximum_group_member_preview_count
        ]
        _remaining_group_locus_count = (
            len(high_bsr_group.locus_names) - len(_visible_group_loci)
        )
        _group_locus_label = ", ".join(_visible_group_loci)
        if _remaining_group_locus_count:
            _group_locus_label += (
                f" … (+{_remaining_group_locus_count:,} more)"
            )
        representative_bsr_group_rows.append(
            {
                "Group": high_bsr_group.group_number,
                "Locus Count": len(high_bsr_group.locus_names),
                "Qualifying Pair Count": high_bsr_group.qualifying_pair_count,
                "Maximum Directional BSR": round(
                    high_bsr_group.maximum_blast_score_ratio,
                    6,
                ),
                "Member Loci (preview)": _group_locus_label,
            }
        )
    if representative_bsr_group_rows:
        representative_bsr_group_table = mo.ui.table(
            representative_bsr_group_rows,
            pagination=True,
            page_size=25,
            selection=None,
            show_column_summaries=False,
            show_data_types=False,
            text_justify_columns={
                "Group": "right",
                "Locus Count": "right",
                "Qualifying Pair Count": "right",
                "Maximum Directional BSR": "right",
            },
            column_widths={
                "Group": 110,
                "Locus Count": 150,
                "Qualifying Pair Count": 210,
                "Maximum Directional BSR": 245,
                "Member Loci (preview)": 700,
            },
            show_search=True,
            show_download=False,
            max_columns=None,
            max_height=500,
        )
    else:
        representative_bsr_group_table = mo.callout(
            "No high-BSR inter-locus groups are present at this threshold.",
            kind="info",
        )

    _covered_allele_percentage = (
        f"{representative_bsr_classification.covered_allele_percentage:.2f}%"
        if representative_bsr_classification.covered_allele_percentage
        is not None
        else "unavailable"
    )
    representative_bsr_analysis_displayed = True
    mo.vstack(
        [
            mo.md(
                f"""
                ## Representative BSR Summary Statistics

                **Availability:** an approximate representative was selected
                for
                **{schema_representative_bsr_analysis.representative_locus_count:,}**
                of **{schema_representative_bsr_analysis.locus_count:,}** loci;
                **{schema_representative_bsr_analysis.self_score_locus_count:,}**
                had a usable BLASTp self score. The protein analysis evaluated
                **{schema_representative_bsr_analysis.valid_translated_allele_count:,}**
                valid alleles from
                **{schema_representative_bsr_analysis.total_allele_count:,}**
                total records.

                **Within-locus coverage:**
                **{representative_bsr_classification.covered_allele_count:,}**
                of
                **{representative_bsr_classification.measurable_allele_count:,}**
                measurable valid alleles are covered at BSR ≥
                **{representative_bsr_classification.blast_score_ratio_threshold:g}**
                (**{_covered_allele_percentage}**).
                **{representative_bsr_classification.incomplete_coverage_locus_count:,}**
                measurable loci contain at least one uncovered valid allele;
                **{representative_bsr_classification.unavailable_valid_allele_count:,}**
                otherwise-valid alleles lack a usable representative self
                score.

                **Inter-locus high-BSR evidence:**
                **{len(representative_bsr_classification.high_bsr_pairs):,}**
                unique locus pairs involving
                **{representative_bsr_classification.high_bsr_locus_count:,}**
                loci form
                **{len(representative_bsr_classification.high_bsr_groups):,}**
                connected similarity groups. These are descriptive
                similarity results, not confirmed paralogy.

                **Provenance:**
                `{schema_representative_bsr_analysis.blastp_version}`;
                E-value floor
                `{schema_representative_bsr_analysis.blast_evalue:g}`; method
                `{schema_representative_bsr_analysis.method_identifier}`.
                Cache reused:
                **{schema_representative_bsr_analysis.reused_from_cache}**.

                The locus table displays
                **{len(representative_bsr_locus_rows):,}** of
                **{len(_all_locus_coverage_rows):,}** rows, and the evidence
                table displays
                **{len(representative_bsr_pair_rows):,}** of
                **{len(representative_bsr_classification.high_bsr_pairs):,}**
                qualifying pairs. Previews are bounded to protect Marimo from
                oversized serialized outputs; no CSV is generated.
                """
            ),
            mo.md("### Per-locus representative coverage"),
            representative_bsr_locus_table,
            mo.md("### Inter-locus high-BSR similarity groups"),
            representative_bsr_group_table,
            mo.md("### Strongest supporting high-BSR locus pairs"),
            representative_bsr_pair_table,
        ],
        gap=2,
    )
    return (representative_bsr_analysis_displayed,)


@app.cell
def _(mo, representative_bsr_analysis_displayed):
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
        _ = mo.status.toast(
            "Scheme cache cleared",
            "The archive selection was reset.",
        )


if __name__ == "__main__":
    app.run()
