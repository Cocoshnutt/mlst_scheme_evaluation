# FAIR intrinsic evaluation outputs

## Purpose and scope

The intrinsic Marimo application can download the results currently present in
its reactive graph as a versioned ZIP bundle. The bundle is a research-output
artifact, unlike `.cache/`, which is an implementation detail used to avoid
recomputation. Creating the download control does not run MAFFT or BLAST. The
archive checksum and ZIP are calculated only after the user requests the
download.

A bundle may contain only the automatically calculated static results. If some
or all requested expensive analyses have finished, their result tables are
included too. `manifest.json` records each stage as `complete` or `not_run`, so
a partial bundle cannot be mistaken for a complete intrinsic evaluation.

## Bundle contract

Bundle format version 1 always contains:

| File | Meaning |
| --- | --- |
| `manifest.json` | Bundle format version, creation time, scheme identifier, stage status, submitted operational settings, input identities, tool/method and cache-reuse provenance, file hashes, missing-value conventions, and interpretation warning. |
| `input_provenance.json` | Repository-relative archive path, archive size and SHA-256, the matching source-inventory record, explicit placeholder fields, and recorded-checksum comparison when available. |
| `summary.json` | Compact scheme-wide static and available expensive-stage summaries, including active screening settings and tool/method provenance. |
| `locus_statistics.tsv` | Complete typed per-locus length-distribution and nucleotide-composition results. |
| `locus_allele_information.tsv` | Complete typed per-locus integrity and protein-diversity results. |
| `software_versions.json` | Project, Python, platform, and available MAFFT/BLAST versions. |
| `checksums.sha256` | SHA-256 digest for every other file in the ZIP. |

Available requested stages additionally contribute:

- `alignment_loci.tsv`;
- `cross_locus_similarity.tsv` and, when classified,
  `candidate_paralog_pairs.tsv`;
- `representative_bsr_loci.tsv` and
  `representative_bsr_inter_locus.tsv`; and
- when representative BSR has been classified,
  `representative_bsr_coverage.tsv`, `representative_bsr_pairs.tsv`, and
  `representative_bsr_groups.tsv`.

JSON uses `null` for unavailable values. TSV uses the literal `unavailable`,
not an empty cell. Tuple, list, and nested dataclass fields are represented as
compact JSON inside a TSV cell. This retains typed evidence without creating a
large collection of sparsely useful secondary tables. Field names are stable,
machine-readable identifiers taken from the public result dataclasses.

The bundle records the screening thresholds currently selected in the
notebook. Threshold-dependent tables are therefore a reproducible view of the
cached raw measurements. Candidate-paralog and high-BSR results remain
screening evidence; neither is an overall scheme pass/fail decision.

The manifest also records the submitted shared worker budget. This is
operational provenance only: worker count affects scheduling and remains
outside every scientific cache key. MAFFT provenance reports reused and newly
calculated locus counts because one request may reuse only part of its
per-locus cache. BLAST stages report whether the complete compatible analysis
was reused.

## FAIR application

### Findable

- Stable `scheme_id` values come from `metadata/scheme_sources.toml`.
- Bundle and source versions, file names, and SHA-256 digests support indexing
  and exact identification.
- Published releases should be archived in a persistent repository and cite a
  release DOI when one exists.

### Accessible

- JSON, TSV, and SHA-256 text files can be read without this software.
- Data that cannot be redistributed should retain a source URL, access terms,
  and retrieval instructions in the source inventory rather than being copied
  into the bundle unlawfully.

### Interoperable

- FASTA remains the sequence-input convention; JSON and TSV are the result
  conventions.
- Organisms should use verified NCBI Taxonomy identifiers, and later genome
  inputs should use versioned assembly accessions.
- Adding Parquet for scale is acceptable later, but it must not replace the
  open rectangular TSV representation.

### Reusable

- The bundle includes input identity, method identifiers, tool versions,
  settings, explicit missingness, and checksums.
- The software is licensed under EUPL-1.2. External scheme licenses remain
  separate and must be recorded per source.
- Scientific definitions and denominators remain documented in the notebook,
  module contracts, and tests; cache data is not treated as the normative
  result format.

## Publication checklist and known metadata gaps

Before publishing an evaluation bundle:

1. Replace every `TODO` in `metadata/scheme_sources.toml` with verified
   information, or use `unknown` with an explanation when a historical fact
   cannot be recovered.
2. Set the inventory's `provenance_status` to `complete` only after reviewing
   every record. The exporter reports provenance as incomplete until both this
   status and the selected record are complete.
3. Verify that the recorded archive SHA-256 matches the digest exported in
   `input_provenance.json`.
4. Confirm whether the selected archive may be redistributed and preserve the
   original provider's citation and license terms.
5. Run the complete intrinsic evaluation if the publication claims to include
   MAFFT and BLAST results; check `manifest.json` rather than inferring
   completeness from filenames.
6. Record the release version and persistent identifier used to publish the
   bundle.

`CITATION.cff` records the confirmed project author and contact email. Add an
ORCID or affiliation only after the project owner confirms it, and add a
release date or persistent identifier when a corresponding release is
published rather than guessing those values in advance.

## Extending the bundle

When a new intrinsic analysis is added, update `scheme.export` in the same
change. Export its typed rows, add an explicit stage status, include its compact
scheme-level metrics in `summary.json`, and test both absent and available
states. Increment the bundle format version if a field is removed, renamed, or
changes meaning. Adding a new optional file or field may remain backward
compatible, but it still requires documentation and tests.
