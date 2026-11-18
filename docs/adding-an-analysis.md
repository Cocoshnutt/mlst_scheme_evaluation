# Adding an evaluation analysis

Keep new analyses small, typed, and scientifically explicit. Add ordinary
library code first; the Marimo notebook should orchestrate and present it.

## Developer checklist

- Define immutable raw measurement and threshold-dependent classification
  result types where both concepts exist.
- State missing, unavailable, failed, excluded, and unevaluated cases
  explicitly.
- Name every numerator and denominator used by a count or percentage.
- Keep screening thresholds separate from expensive raw measurements.
- Define the cache version, scientifically meaningful key fields, provenance,
  atomic-write behavior, and invalidation rules.
- Add typed progress reporting when work is expensive enough to require user
  feedback.
- Add scheme-level summary fields to `scheme/overview.py` when the measurement
  belongs in the intrinsic overview.
- Keep detailed rows and plots bounded without changing complete reported
  counts.
- Add persistent cache filenames to the relevant notebook cache-status group.
- Add focused temporary-input tests, boundary cases, cache resilience tests,
  and cross-module consistency tests.
- Update scientific methodology and interpretation documentation when the
  implemented behavior changes.

Do not mark an analysis implemented in the roadmap until its scientific result
and notebook presentation both exist. A parser, manifest, or placeholder UI is
infrastructure rather than a completed measurement.

## Documentation placement

Use the notebook for concise purpose, current settings, key formulas, threshold
sources, limitations, and result interpretation. It should contain enough
context to interpret the visible output without becoming the canonical home of
every technical detail.

Use module docstrings for module responsibility and the scientific contract.
Use public-function docstrings for exact inputs, outputs, errors, missingness,
and cache behavior. Document private functions when an algorithmic choice is
non-obvious or scientifically meaningful.

Keep CLI `--help` focused on operational usage. Put complete methodology,
provenance requirements, cache design, and developer guidance under `docs/`.

Do not make a notebook import large Markdown documents merely to eliminate all
repetition. Small audience-appropriate repetition improves usability, while
normative definitions should still have one canonical source.

## Scientific boundaries

- Do not change a scientific definition while moving or reorganizing code.
- Label descriptive project heuristics separately from literature-derived
  screening values.
- Do not convert a screening flag into a scheme-level pass/fail decision.
- Keep intrinsic sequence similarity separate from genome-level copy-number or
  functional evidence.
- Keep threshold-independent intermediates reusable when only classification
  settings change.
- Record method and input provenance sufficient to explain a result without
  relying on internal cache files.

## Phase 2 additions

Before adding a functional-performance metric, define the genome manifest,
typing-run provenance, normalized call states, excluded states, and join keys.
Follow `docs/functional-workflow.md`. Do not make Phase 2 calculations depend
on executing or importing the intrinsic notebook; combine exported typed
results by stable identifiers instead.
