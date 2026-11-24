# Architecture

## Scope and design boundary

This repository contains two related but independently executable evaluation
workflows:

- `notebooks/intrinsic.py` evaluates properties measurable directly from an
  MLST or cgMLST allele scheme.
- `notebooks/functional.py` prepares the genome-dependent Phase 2 workflow.
  It currently validates representative-genome manifests; functional
  calculations are not implemented yet.

The notebooks share small library modules and stable identifiers, but they do
not import or execute one another. Their caches remain independent. A future
combined report should join their versioned result artifacts by scheme and
locus identifiers rather than reaching into either notebook or cache.

## Repository data flow

```text
Intrinsic workflow                             Functional workflow
==================                             ===================

archive/**/*.tar.*                             metadata/genome_sets/*.toml
        |                                                  |
        v                                                  v
scheme.cache                                      functional.manifest
atomic extraction                                load + readiness validation
        |                                                  |
        v                                                  v
.cache/schemes/<fingerprint>/                  typed GenomeSetManifest
        |                                                  |
        +--> shared FASTA/sequence rules                    +--> [planned]
        |    scheme.fasta                                   |    call importers
        |    scheme.sequence                                |    run provenance
        |    scheme.selection                               |    normalization
        |                                                   |
        +--> static intrinsic measurements                  +--> [planned]
        |    info / stats / allele_info                     |    recovery
        |                                                   |    completeness
        +--> requested expensive measurements               |    multicopy
        |    alignment / similarity / bsr                   |    reproducibility
        |                                                   |    tool agreement
        v                                                   v
scheme.overview                                  typed functional results
scheme.charts                                    [not implemented]
        |                                                   |
        v                                                   v
notebooks/intrinsic.py                           notebooks/functional.py
        |                                                   |
        v                                                   v
versioned intrinsic ZIP                         versioned functional bundle
JSON + TSV + provenance                         [planned]
        \                                                   /
         +------ future artifact-based combined report ----+
```

`.cache/` accelerates computation but is never a research-output interface.
The intrinsic ZIP is the current portable result boundary. The corresponding
functional artifact contract is documented before implementation so Phase 2
does not become coupled to Phase 1 internals.

## Intrinsic workflow

### Input and extraction

The selected repository archive is extracted by `scheme.cache` into
`.cache/schemes/` using an archive fingerprint and an atomic completion marker.
The extracted directory is treated as immutable for the evaluation session.
Input source metadata is matched by repository-relative archive path in
`metadata/scheme_sources.toml`.

### Shared scientific primitives

- `scheme.fasta` owns recognized FASTA suffixes, recursive locus discovery,
  plain/gzip reading, and record boundaries.
- `scheme.sequence` owns normalization, nucleotide composition, CDS validity,
  and translation.
- `scheme.selection` owns deterministic modal-length-proximal representative
  selection.

Scientific analyses consume these modules rather than importing private
helpers from one another.

### Measurements and classification

- `scheme.static_analysis` performs the shared one-pass static FASTA scan.
  `scheme.info`, `scheme.stats`, and `scheme.allele_info` project its typed
  measurements into their established public APIs.
- `scheme.alignment` performs requested per-locus MAFFT nucleotide and protein
  alignments. Its cache is per strategy and per locus, so one request may mix
  reused and newly calculated loci.
- `scheme.similarity` performs requested selected-allele BLASTn/BLASTp
  measurements. Candidate-paralog thresholds classify the cached raw results
  without rerunning BLAST.
- `scheme.bsr` performs requested representative-to-allele BLASTp
  measurements. Its threshold-dependent coverage and high-BSR groups are also
  derived without rerunning BLAST.
- `scheme.blast` owns checked BLAST process and temporary database mechanics;
  each analysis retains its own command arguments, parser, cache, and
  scientific interpretation.

### Notebook dependency model

Static analyses run automatically after scheme selection. MAFFT and BLAST
cells are gated by individual controls and one complete-evaluation control.
`scheme.overview` aggregates already available typed results, allowing the
top-of-notebook overview to show static values immediately and requested
results progressively. Rendering cells do not act as execution barriers.

Reactive classification thresholds are separate from expensive raw analyses.
Changing a threshold updates overview and detail presentation while reusing
the same threshold-independent cached measurements.

### Result artifact and operational provenance

`scheme.export` creates the intrinsic ZIP lazily when the download is
requested. `manifest.json` records:

- stage completion;
- the submitted worker budget and selected MAFFT strategy;
- active static, alignment, similarity, and BSR thresholds;
- scheme/archive and scheme-source-inventory identities;
- method identifiers and external-tool versions; and
- whole-analysis BLAST cache reuse or per-locus MAFFT reuse counts.

Worker budget affects scheduling only. It is intentionally present in
operational provenance and absent from scientific cache keys. Tool method,
version, input identity, and scientifically meaningful parameters continue to
control cache compatibility where applicable.

## Functional workflow

### Current implementation

`functional.manifest` loads a versioned genome-set TOML document, calculates
its SHA-256 identity, and returns typed genome and dataset metadata. Validation
separates blocking input errors from non-blocking provenance warnings.
`notebooks/functional.py` presents this readiness state only; it does not yet
run typing tools or display fabricated performance results.

The manifest supplies the future join to the intrinsic workflow through its
stable `scheme_id`. Versioned assembly accessions, checksums, and genome IDs
will identify functional inputs.

### Planned calculation boundary

Before calculations are added, Phase 2 still needs real example outputs and a
typed normalization contract for exact, novel, partial, absent, multicopy,
failed, excluded, and unevaluated allele calls. Once defined, importers should
produce normalized typed tables before metrics calculate:

1. per-locus recovery/prevalence;
2. per-genome completeness and missing-call patterns;
3. genome-level multicopy evidence;
4. explicitly paired reproducibility comparisons; and
5. typing-tool concordance and discordance.

Phase 2 will use `.cache/functional/` with genome-manifest, scheme, normalized
input, method, and tool-version identities. Worker budgets and display settings
will remain outside cache keys. Its publishable output will be an independent
open bundle with manifest, provenance, typed tables, and checksums.

## Shared versus phase-specific code

| Concern | Location | Boundary |
| --- | --- | --- |
| FASTA parsing, CDS rules, representative selection | `scheme/fasta.py`, `sequence.py`, `selection.py` | Shared foundational sequence behavior. |
| Intrinsic measurements and caches | `scheme/info.py`, `stats.py`, `allele_info.py`, `alignment.py`, `similarity.py`, `bsr.py` | Scheme-archive inputs only. |
| Intrinsic aggregation and presentation projections | `scheme/overview.py`, `charts.py`, `export.py` | No external-tool execution. |
| Genome-set input contract | `functional/manifest.py` | Phase 2 manifest loading and validation. |
| Functional calculations | `functional/` | Planned; must not import notebook code. |
| User interfaces | `notebooks/intrinsic.py`, `notebooks/functional.py` | Orchestration and interpretation, not scientific parsing. |

## Cache and artifact rules

1. Cache raw, threshold-independent measurements when classification settings
   can change.
2. Include scientifically meaningful method and input identities in cache
   compatibility; exclude worker counts and presentation settings.
3. Keep writes atomic and tolerate persistence failure when readable inputs can
   still be analyzed.
4. Treat failed external-tool work as retryable, not as a reusable biological
   result.
5. Export typed open results and explicit missing states; never publish cache
   JSON as the result contract.
6. Join phases through stable identifiers and versioned artifacts, not shared
   mutable state.

## Related documentation

- User setup, implemented analyses, and roadmap: `README.md`
- Intrinsic output bundle and FAIR publication checklist:
  `docs/fair-outputs.md`
- Phase 2 inputs, missingness, caches, and planned artifacts:
  `docs/functional-workflow.md`
- Adding a typed analysis: `docs/adding-an-analysis.md`
- Shared static-cache decision and benchmark: `docs/static-analysis-cache.md`
