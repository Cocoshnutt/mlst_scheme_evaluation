# Phase 2 functional-performance workflow

## Boundary

Phase 1 measures intrinsic properties of scheme allele FASTA files. Phase 2
measures functional performance against representative genomes and typing-tool
outputs. They use separate Marimo notebooks and independent caches because
their inputs, failure modes, runtimes, and provenance differ. Stable scheme and
locus identifiers provide the join for a later combined report.

`notebooks/functional.py` currently validates the genome-set boundary only. It
does not run genome typing or calculate performance metrics until the actual
input formats and scientific denominators are specified.

## Required inputs

Every evaluation dataset needs:

- a completed scheme record in `metadata/scheme_sources.toml`;
- a versioned genome-set TOML manifest conforming to
  `metadata/genome_manifest.schema.json`;
- accessible assembly files whose SHA-256 digests match the manifest;
- a documented selection rationale and relevant lineage representation;
- assembly-quality metadata, with unavailable fields identified explicitly;
- one run manifest per typing tool and configuration;
- normalized allele-call tables with a documented missing/call-state
  vocabulary; and
- explicit pairings for repeated runs or sequencing-depth comparisons.

The template at `metadata/genome_sets/template.toml` is a place to record the
genome provenance requested during the refactor. It must be copied to a stable
dataset-specific filename and completed with verified values. Unknown optional
information should be written as `unknown` and explained rather than inferred.

## Analysis contracts that still require real input examples

The following choices must be made from actual tool outputs rather than
guessed during architectural work:

- the canonical representation of exact, novel, partial, missing, failed, and
  multiple allele calls;
- how source-tool and cross-tool allele identifiers map to scheme sequences;
- whether recovery counts a sequence hit, a valid allele call, or both as
  separate measurements;
- the denominators for genomes excluded from a locus or loci excluded from a
  tool run;
- what constitutes a reproducibility pair;
- how multicopy evidence is represented by every selected tool; and
- whether proprietary outputs may be stored or only derived open tables may be
  published.

These are scientific data-contract decisions, not notebook presentation
details.

## Cache design

Phase 2 caches will live separately under `.cache/functional/` and be keyed by
at least:

- genome-manifest SHA-256;
- scheme/archive identity;
- normalized call-input checksum;
- analysis cache version;
- typing tool and version where measurements depend on it; and
- scientifically meaningful parameters.

Threshold-independent measurements must be cached separately from screening or
display classifications. Worker counts and presentation settings must not be
cache keys. No Phase 2 cache is created by the current readiness notebook.

## Missingness and denominators

Phase 2 result types must distinguish at least:

- measured present or called;
- measured absent;
- partial or below the tool's callable criteria;
- multiple hits or calls;
- tool failure;
- input unavailable;
- excluded by an explicit rule; and
- not evaluated.

Percentages must name their numerator, denominator, and excluded states. Locus
prevalence across genomes must not be conflated with scheme completeness within
one genome.

## FAIR result artifacts

Caches are implementation details. A publishable functional evaluation should
export a versioned open bundle such as:

```text
evaluation/
├── manifest.json
├── genome_set.toml
├── run_provenance.json
├── summary.json
├── locus_prevalence.tsv.gz
├── genome_completeness.tsv.gz
├── multicopy_evidence.tsv.gz
├── reproducibility.tsv.gz
├── tool_concordance.tsv.gz
├── software_versions.json
└── checksums.sha256
```

Each table requires stable metric identifiers, units, denominator fields, and
explicit missing-value states. TSV or compressed TSV should remain available
even if Parquet is added for scale. A combined report may consume these
artifacts, but should not require access to internal caches.

## Phase 3

Do not introduce Phase 3 abstractions yet. Cross-organism comparison should
reuse proven Phase 2 result artifacts after at least one real functional
evaluation establishes the data model.
