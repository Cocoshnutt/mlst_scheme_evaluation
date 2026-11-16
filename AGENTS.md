# Repository agent guide for `mlst_scheme_evaluation`

Use this file to acquire the repository context quickly without reading the
scheme archives, generated caches, or lock file. This is a Marimo-based MLST
and cgMLST scheme evaluator. It extracts a selected archive, calculates
scheme-wide and per-locus FASTA metrics, caches expensive static analysis, and
renders the results in `notebook.py`.

## Fast context bootstrap

Start every task with these read-only checks:

```console
git status --short
rg --files -g '!uv.lock' -g '!.cache/**' -g '!__pycache__/**' -g '!archive/**'
sed -n '1,220p' README.md
sed -n '1,220p' pyproject.toml
```

The worktree may already contain user changes or untracked files. Do not
discard, reset, format, or delete them. Avoid opening `uv.lock`, `.cache/`, or
the binary files under `archive/` unless the task specifically concerns
them.

Next, inspect only the module relevant to the task. A compact API outline is
usually sufficient before reading implementations:

```console
rg -n '^(class |def |async def |[A-Z][A-Z0-9_]+\s*=|__all__|@app\.cell)' \
  schema/cache.py schema/info.py schema/stats.py schema/allele_info.py \
  schema/alignment.py schema/similarity.py notebook.py
```

| Task | Read first | Main public API |
| --- | --- | --- |
| Archive extraction or clearing | `schema/cache.py` | `cache_archive`, `clear_cache` |
| Scheme-wide sequence metrics | `schema/info.py` | `get_schema_information` and count helpers |
| Locus distributions and charts | `schema/stats.py` | `get_schema_statistics`, `create_schema_statistics_charts` |
| Per-locus allele table | `schema/allele_info.py` | `get_scheme_allele_information` |
| Within-locus alignments | `schema/alignment.py` | `get_schema_alignment_analysis`, `calculate_schema_alignment_summary` |
| Cross-locus similarity | `schema/similarity.py` | `get_schema_similarity_analysis` |
| UI, controls, or layout | Relevant library plus `notebook.py` | Marimo cells |
| Dependencies or development environment | `README.md`, `pyproject.toml`, `devenv.nix` | `uv sync --locked`; generated `environment.yml`; optional `devenv shell` |
| Planned metrics and research criteria | `README.md` TODO | Treat thresholds as advisory flags |

## Architecture and data flow

```text
archive/**/*.tar.*
        |
        v
schema.cache.cache_archive()
        |
        v
.cache/schemes/<archive-name>-<12-character fingerprint>/
        |
        +--> schema/info.py --------> scheme summary
        +--> schema/stats.py -------> four Altair panels
        +--> schema/allele_info.py -> searchable per-locus table
        +--> schema/alignment.py ---> nucleotide/protein gap statistics
        +--> schema/similarity.py --> nucleotide/protein cross-locus similarity
                                      |
                                      v
                                  notebook.py
```

Keep analysis in importable library modules and presentation in
`notebook.py`. Do not duplicate FASTA parsing or metric calculations in the
notebook. When adding a substantial analysis, expose a typed public function
and let the notebook consume its return value. Existing command-line entry
points are available in `schema/cache.py`, `schema/info.py`, and
`schema/allele_info.py`.

### Module responsibilities

- `schema/cache.py` validates and extracts tar archives atomically. Cache
  identity is based on the archive's resolved path, size, and modification
  time. A `.complete` marker distinguishes reusable extractions from partial
  work.
- `schema/info.py` computes aggregate scheme metrics: loci, alleles, CDS
  validity categories, configurable short-allele counts, and scheme-wide
  80%/120% modal-length outlier totals.
- `schema/stats.py` computes per-locus allele-count and length summaries and
  builds four vertically stacked Altair charts.
- `schema/allele_info.py` computes chewBBACA-style per-locus integrity and
  diversity fields for the Marimo table.
- `schema/alignment.py` runs cached, parallel MAFFT nucleotide and protein
  alignments and derives per-locus gap proportions without retaining full
  alignments.
- `schema/similarity.py` selects one modal-length-proximal allele per locus,
  runs cached nucleotide and protein BLAST searches, and derives directional
  inter-locus identity, coverage, and BSR measurements.
- `notebook.py` owns selection, status notices, the minimum-length control,
  display formatting, charts, the allele table, and the clear-cache action.
- `README.md` documents use and contains the phased analysis roadmap. Preserve
  the distinction between intrinsic FASTA-only quality and genome-dependent
  functional performance.

## Current behavioral invariants

Preserve these unless the user explicitly requests a semantic change:

- Python is `>=3.14,<3.15`; modern syntax such as `type` aliases is expected.
- Recognized FASTA suffixes are defined once in `schema/info.py`, including
  gzip-compressed variants.
- Library scanners currently treat every recognized FASTA file below an
  extracted schema recursively as a locus. The current scheme collection has
  no `short/` directories. Supporting curated representatives in future
  inputs remains a separate TODO, not current special behavior.
- BUSCO archives under `archive/BUSCO_bac/` remain on disk but are deliberately
  excluded from the notebook dropdown.
- A valid coding allele has only A/C/G/T bases, a length divisible by three,
  a bacterial start codon (`ATG`, `GTG`, or `TTG`), a terminal stop codon, and
  no earlier in-frame stop codon.
- `invalid_allele_count` is the requested sum of invalidity-category counts.
  One allele may contribute to multiple categories, so this value is not
  necessarily the number of distinct invalid records.
- Modal-length ties select the shortest tied length. Outliers are strictly
  below 80% or strictly above 120% of that mode; values on the boundaries are
  not outliers. Empty loci contribute no outliers.
- Changing the minimum allele length must not rescan every FASTA. Static
  metrics and lengths are cached; only the threshold-dependent count is
  recalculated.
- Empty-sequence and within-locus duplicate-sequence metrics are static and do
  not depend on the minimum allele-length setting. An empty sequence is a
  FASTA record with no normalized nucleotides; a zero-record FASTA is an
  affected locus but contributes no empty sequence record.
- Within-locus duplicate counts include redundant non-empty records beyond the
  first exact match after whitespace removal and uppercasing. They do not
  compare different loci or reverse complements. Cross-locus shared-sequence
  detection is a separate future analysis.
- The notebook's four locus charts are displayed underneath one another.
  Panel C plots allele size on x against allele count on y. Panel D uses a
  lower interval selector to update the main plot's visible x-domain.
- Clearing the cache removes extracted schemes and embedded analysis caches,
  resets the selected archive, shows a toast, and makes the next analysis
  slower. Keep both the initial-computation and clear-cache warnings.
- The notebook uses Inter through a Google Fonts import with a sans-serif
  fallback.
- Literature values in the README are study-specific screening references.
  Display measured values and flagged loci; do not silently classify an
  entire scheme as suitable or unsuitable.
- Similarity-only results must be labelled candidate paralogous loci. Confirmed
  copy number or paralogy requires genome-level evidence.
- Cross-locus BLAST currently selects one allele directly from each locus
  FASTA rather than using a curated representative. It chooses the valid CDS
  closest to the shortest tied modal allele length, with FASTA order breaking
  equal-distance ties. If no valid CDS exists, an eligible IUPAC DNA allele is
  selected by the same rule for nucleotide-only analysis. This is a
  deterministic length-based approximation, not a sequence medoid.
- Within-locus alignments default to MAFFT FFT-NS-2. FFT-NS-1 is the optional
  fastest strategy with a rougher guide tree; Auto may select slower iterative
  methods. Strategy choice can change the alignment and measured gap
  proportions, so always report it with results.
- Each MAFFT subprocess uses one thread while loci are scheduled concurrently
  across the CPUs available to the process. Inputs with at least two eligible,
  exactly identical sequences use an exact zero-gap fast path.
- Nucleotide polymorphism is substitution-only. A nucleotide MSA column is
  callable when at least two sequences contain canonical A/C/G/T observations
  after treating gaps and ambiguous IUPAC symbols as missing. It is
  polymorphic when those observations contain at least two distinct bases;
  percentage is polymorphic sites divided by callable sites. Indel variation
  remains in the separate gap metrics.

## Cache model

The extraction cache root is `.cache/schemes`. Each extracted archive contains
independent JSON analysis caches:

| File | Owner | Invalidation |
| --- | --- | --- |
| `.schema-information-cache.json` | `schema/info.py` | Cache version; extracted archives are treated as immutable |
| `.schema-statistics-cache.json` | `schema/stats.py` | Cache version plus FASTA path/size/mtime signature |
| `.scheme-allele-information-cache.json` | `schema/allele_info.py` | Cache version plus FASTA path/size/mtime signature |
| `.schema-alignment-cache.json` | `schema/alignment.py` | Cache version, MAFFT version, strategy, and per-locus FASTA path/size/mtime signature |
| `.schema-similarity-cache.json` | `schema/similarity.py` | Cache version, BLAST versions, method, E-value, and FASTA path/size/mtime signature |

The alignment cache v4 keeps FFT-NS-2, FFT-NS-1, and Auto results in separate
method namespaces. Version-one and version-two alignment caches migrate into
the Auto namespace. Version-three caches retain their protein results and run
only the nucleotide alignment needed to add polymorphism metrics. Preserve
other strategy namespaces when checkpointing one strategy, and retry failed
entries instead of treating them as reusable biological results.

When cached fields or their meaning change, increment the corresponding cache
version. Keep writes atomic with a temporary file followed by replacement.
Failure to write a cache should not prevent analysis of a readable scheme.
Do not commit `.cache/` contents.

## Marimo conventions

Marimo derives dependencies from cell parameters and returned names. When
editing `notebook.py`:

- Return every value consumed by another cell and include it as a parameter in
  the consumer.
- Avoid defining the same global name in multiple cells.
- Keep expensive, threshold-independent work in a separate cell from reactive
  controls.
- Add a library module to the first cell's `mo.watch.file(...)` and reload flow
  when the notebook must hot-reload it during editing.
- Update `analysis_cache_filenames` when adding a persistent analysis cache so
  the initial-computation notice remains accurate.
- Build table rows explicitly and keep column widths/tooltips in sync with
  added fields; wide automatic headers have previously broken the table.
- Keep user-facing sections clearly separated and keep all primary plots at a
  consistent width unless a compact selector intentionally differs.

## Implementation style

- Use explicit signatures, type hints, variable names, module docstrings, and
  public-function docstrings. Prefer `output_file_path` over shorthand such as
  `outfile`.
- Preserve immutable result dataclasses and tuples used by the current public
  APIs unless a change requires otherwise.
- Reuse the shared FASTA constants and validation semantics instead of creating
  nearly equivalent definitions.
- Handle empty loci, malformed FASTA records, gzip input, tied modes, and
  read-only extracted directories.
- Use atomic caches and source signatures for expensive new static analyses.
- Keep threshold-independent intermediates cached so UI threshold changes are
  cheap.
- Do not use the real scheme archives for unit fixtures. Create small temporary
  FASTA files that isolate empty, malformed, tied-mode, boundary, and duplicate
  cases.
- Update README usage or TODO checkboxes only when behavior is actually
  implemented.

## Dependency environment exports

`environment.yml` is generated from `pyproject.toml` and `uv.lock` by
`scripts/sync_conda_environment.py`. Do not edit it directly. Devenv runs the
generator after its uv synchronization task and before entering the shell. For
dependency changes made outside devenv, regenerate it with:

```console
uv run --locked python scripts/sync_conda_environment.py
```

Use the same command with `--check` in validation. The generated Conda
specification installs MAFFT and BLAST+ from Conda channels and exact locked
Python packages through pip; it does not create or mutate a user's Conda
environment.

## Validation

Run the smallest relevant checks first. The baseline static validation is:

```console
uv run --locked python scripts/sync_conda_environment.py --check
uv run --locked python -m py_compile \
  schema/cache.py schema/info.py schema/stats.py schema/allele_info.py \
  schema/alignment.py schema/similarity.py notebook.py
uv run --locked ruff check \
  schema/cache.py schema/info.py schema/stats.py schema/allele_info.py \
  schema/alignment.py schema/similarity.py notebook.py
```

For notebook changes, also run:

```console
uv run --locked marimo check notebook.py
```

`marimo check` has occasionally stalled in this environment. Do not report it
as passed unless it exits successfully; stop it and report the limitation if
it remains silent for an unreasonable time.

For library changes, add focused temporary-directory tests or exercise the
relevant CLI without modifying a real archive:

```console
uv run --locked python -m schema.cache --help
uv run --locked python -m schema.info --help
uv run --locked python -m schema.allele_info --help
uv run --locked python -m schema.alignment --help
uv run --locked python -m schema.similarity --help
```

After validation, inspect `git status --short` again and report only files and
artifacts changed for the task. Never clear the user's cache merely to test a
change unless they explicitly requested it.
