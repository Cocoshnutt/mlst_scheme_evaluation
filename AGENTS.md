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
  schema/alignment.py schema/similarity.py schema/bsr.py notebook.py
```

| Task | Read first | Main public API |
| --- | --- | --- |
| Archive extraction or clearing | `schema/cache.py` | `cache_archive`, `clear_cache` |
| Scheme-wide sequence metrics | `schema/info.py` | `get_schema_information` and count helpers |
| Locus distributions and charts | `schema/stats.py` | `get_schema_statistics`, `create_schema_statistics_charts` |
| Per-locus allele table | `schema/allele_info.py` | `get_scheme_allele_information` |
| Within-locus alignments | `schema/alignment.py` | `get_schema_alignment_analysis`, `calculate_schema_alignment_summary` |
| Cross-locus similarity | `schema/similarity.py` | `get_schema_similarity_analysis` |
| Representative-to-allele protein BSR | `schema/bsr.py` | `get_schema_representative_bsr_analysis`, `classify_representative_bsr_analysis` |
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
        +--> schema/bsr.py ---------> representative coverage/high-BSR groups
                                      |
                                      v
                                  notebook.py
```

Keep analysis in importable library modules and presentation in
`notebook.py`. Do not duplicate FASTA parsing or metric calculations in the
notebook. When adding a substantial analysis, expose a typed public function
and let the notebook consume its return value. Existing command-line entry
points are available in `schema/cache.py`, `schema/info.py`,
`schema/allele_info.py`, `schema/alignment.py`, `schema/similarity.py`, and
`schema/bsr.py`.

### Module responsibilities

- `schema/cache.py` validates and extracts tar archives atomically. Cache
  identity is based on the archive's resolved path, size, and modification
  time. A `.complete` marker distinguishes reusable extractions from partial
  work.
- `schema/info.py` computes aggregate scheme metrics: loci, alleles, CDS
  validity categories, configurable short-allele counts, scheme-wide 80%/120%
  modal-length outlier totals, and pooled scheme-wide GC content.
- `schema/stats.py` computes per-locus allele-count, length, and GC-content
  summaries, classifies configurable frequent length peaks from cached
  exact-length histograms, and builds the vertically stacked Altair charts.
- `schema/allele_info.py` computes chewBBACA-style per-locus integrity and
  diversity fields for the Marimo table.
- `schema/alignment.py` runs cached, parallel MAFFT nucleotide and protein
  alignments and derives per-locus gap proportions without retaining full
  alignments.
- `schema/similarity.py` selects one modal-length-proximal allele per locus,
  runs cached nucleotide and protein BLAST searches, and derives directional
  inter-locus identity, coverage, and BSR measurements. It also collapses
  reciprocal hits into unique candidate-paralog locus pairs using separately
  configurable literature screening values.
- `schema/bsr.py` selects the same modal-length-proximal valid CDS as an
  approximate representative, searches all selected representatives against
  all valid translated alleles in one cached BLASTp analysis, and reports
  within-locus coverage plus connected inter-locus high-BSR similarity groups.
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
- Locus statistics retain all tied modes and use the shortest as the primary
  mode. A frequent length peak is an exact length whose count is at least the
  configurable percentage of the dominant exact-length count, defaulting to
  50%. One frequent peak is single-mode and more than one is multimodal;
  categories separately retain whether strict 80/120% primary-mode outliers
  exist. This peak rule is a descriptive project heuristic, not a literature
  threshold, and changing it reclassifies cached histograms without rereading
  FASTA files.
- Changing the minimum allele length must not rescan every FASTA. Static
  metrics and lengths are cached; only the threshold-dependent count is
  recalculated.
- Empty-sequence and exact duplicate-sequence detection within or between loci
  are deliberately not evaluation criteria for the current scheme collection.
- Each FASTA record is reported as one nucleotide allele because the current
  scheme collection is expected to contain distinct allele records; do not
  reintroduce sequence hashing for this count. Distinct translated proteins
  are unique translations of alleles that pass the existing CDS-validity
  rules. Synonymous valid nucleotide alleles collapse to one protein, invalid
  or untranslatable alleles contribute no protein, and neither count depends
  on the minimum allele-length control.
- GC content is `(G + C) / (A + C + G + T) * 100`. Whitespace and case are
  normalized; ambiguous and other noncanonical characters are excluded from
  the denominator and counted separately. Duplicate allele records remain in
  the pooled counts. Scheme-wide GC pools bases across all loci rather than
  averaging locus percentages. A scheme or locus without canonical bases has
  an unavailable percentage, not 0%, and GC metrics do not depend on the
  minimum allele-length setting.
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
- Candidate-paralog screening collapses reciprocal hits into one unordered
  locus pair. A pair is flagged by amino-acid identity strictly greater than
  70%, or nucleotide identity of at least 90% across an alignment strictly
  longer than 100 bp. No coverage gate is applied. Directional BSR of at least
  0.6 is supplementary evidence only and cannot flag a candidate by itself.
  These values remain configurable study-specific references.
- Candidate-pair output retains the strongest available nucleotide and
  amino-acid measurement even when that sequence level does not meet its
  criterion. A missing sequence-level measurement means BLAST reported no hit,
  not merely that a hit was below threshold. Reclassification uses cached
  BLAST measurements, does not rerun BLAST, and is independent of the minimum
  allele-length control.
- Cross-locus BLAST currently selects one allele directly from each locus
  FASTA rather than using a curated representative. It chooses the valid CDS
  closest to the shortest tied modal allele length, with FASTA order breaking
  equal-distance ties. If no valid CDS exists, an eligible IUPAC DNA allele is
  selected by the same rule for nucleotide-only analysis. This is a
  deterministic length-based approximation, not a sequence medoid.
- Representative-to-allele BSR is protein-only. It uses the same selected
  modal-length-proximal valid CDS per locus and defines directional BSR as the
  representative-to-allele raw score divided by the representative self raw
  score. Only valid translated alleles enter the coverage denominator;
  invalid or untranslatable records remain unevaluable. Missing BLAST hits are
  distinct from reported below-threshold hits.
- Representative BSR coverage uses an inclusive configurable threshold,
  defaulting to the study-specific 0.6 same-gene reference. Cross-locus hits
  meeting it are collapsed into unordered pairs and transitive connected
  high-BSR similarity groups. These groups are not confirmed paralogs, and
  neither their classification nor the raw analysis depends on the minimum
  allele-length control.
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
- Keep Altair datasets limited to fields actually plotted. Large-schema tables
  must report complete counts but use bounded previews for potentially long
  collections; the locus-length table currently previews at most ten tied
  modes or frequent peaks per locus. Do not raise Marimo's output-size limit to
  mask an unnecessarily large serialized component.

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
| `.schema-representative-bsr-cache.json` | `schema/bsr.py` | Cache version, BLASTp version, method, E-value, and FASTA path/size/mtime signature |

The schema-information cache is currently version 4 and the schema-statistics
cache is currently version 3; these versions add the raw nucleotide counts
needed to reuse scheme-wide and per-locus GC-content calculations.

The scheme-allele-information cache is currently version 6. It stores
per-locus allele metrics, but no empty-sequence counts, within-locus duplicate
counts, cross-locus shared-sequence hashes, or allele-identifier integrity
metrics.

The representative BSR cache is currently version 1. It stores
threshold-independent per-allele BSR values and strongest directed
inter-locus measurements so threshold changes reclassify results without
rerunning BLASTp. Notebook tables use bounded previews and do not generate CSV
downloads for this analysis.

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

The schema-statistics cache v2 stores threshold-independent exact-length
frequency histograms. Frequent-peak threshold changes must reuse these
histograms instead of invalidating the cache or rescanning FASTA files.

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
  schema/alignment.py schema/similarity.py schema/bsr.py notebook.py
uv run --locked ruff check \
  schema/cache.py schema/info.py schema/stats.py schema/allele_info.py \
  schema/alignment.py schema/similarity.py schema/bsr.py notebook.py
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
uv run --locked python -m schema.bsr --help
uv run --locked python -m schema.similarity --help
```

After validation, inspect `git status --short` again and report only files and
artifacts changed for the task. Never clear the user's cache merely to test a
change unless they explicitly requested it.
