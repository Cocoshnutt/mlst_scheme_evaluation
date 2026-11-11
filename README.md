# mlst_scheme_evaluation

Enter the reproducible development shell and open the notebook with:

```console
devenv shell
marimo edit notebook.py
```

Python, uv, marimo, and the notebook dependencies are declared in
`pyproject.toml` and locked in `uv.lock`. `devenv shell` runs `uv sync`
automatically and activates the devenv-managed virtual environment. Use
`uv add <package>` to add or update a dependency.

On a machine without a usable NVIDIA or AMD GPU, marimo 0.25.x may emit a
warning that it cannot collect GPU statistics. It is harmless. To suppress
that upstream warning when desired, start marimo with:

```console
marimo --log-level error edit notebook.py
```

Cache and extract a scheme archive from a notebook cell with:

```python
from schema_cache import cache_archive, clear_cache
from schema_info import get_schema_information

scheme_directory = cache_archive("schemes/mlst/corynebacterium.tar.zst")
schema_information = get_schema_information(
    scheme_directory,
    minimum_sequence_length=300,
)
```

Calling `cache_archive` again reuses the extracted directory. Call
`clear_cache()` to remove all cached schemes. Schema statistics and allele
lengths are cached inside each extracted schema, so repeated evaluations and
new minimum-length thresholds do not reread the FASTA files. Clearing the
scheme cache also clears these calculated results.

The schema information module can also be used from the command line:

```console
python schema_info.py PATH_TO_EXTRACTED_SCHEMA
python schema_info.py PATH_TO_EXTRACTED_SCHEMA --ml 300 --json
```

The standard output includes scheme-wide counts and percentages for alleles
below 80% and above 120% of their respective locus's modal length.

Selecting a scheme in `notebook.py` displays four locus-statistics panels in a
clearly separated section beneath the scheme summary. `schema_stats.py`
provides the reusable calculations and Altair chart construction. Its results
are cached inside the extracted scheme.

The notebook also displays a searchable, paginated allele-analysis table built
by `scheme_allele_info.py`. The table follows the chewBBACA per-locus fields and
reuses cached static metrics when the minimum allele length changes. The module
can also emit the table data as JSON:

```console
python scheme_allele_info.py PATH_TO_EXTRACTED_SCHEMA --ml 201
```

## TODO: scheme evaluation roadmap

Keep intrinsic scheme quality separate from functional scheme performance.
The former can be calculated directly from allele FASTA files; the latter
requires representative genome assemblies and, for some analyses, reference
typing or epidemiological data. Report measured values, the threshold and
literature source used to flag them, and the affected loci. Do not turn
study-specific screening thresholds into a universal scheme pass/fail result.

### Phase 1 — intrinsic scheme quality (no genomes required)

- [x] Report scheme size as the total number of loci.
- [x] Report allelic diversity, including allele counts and their distribution
  across loci.
- [x] Check ambiguous bases, coding-sequence length, start and stop codons, and
  in-frame stop codons.
- [x] Report minimum, maximum, median, modal, and outlier allele lengths.
- [ ] Detect empty and duplicate nucleotide sequences, with the number and
  percentage of affected loci.
- [ ] Align alleles within each locus and report nucleotide and amino-acid gap
  proportions.
- [ ] Calculate the number and percentage of polymorphic nucleotide sites per
  locus.
- [ ] Measure inter-locus nucleotide and amino-acid similarity, including
  identity, coverage, and BLAST Score Ratio (BSR).
- [ ] Report highly similar inter-locus pairs as **candidate paralogous loci**.
  Sequence similarity alone must not be presented as proof of copy number or
  evolutionary paralogy.

#### Schema Refinery statistics integration

- [ ] Add mean allele length, all tied modes, frequent length peaks, and clear
  single-mode/multimodal/outlier locus categories.
- [ ] Add per-locus and scheme-wide GC content.
- [ ] Detect exact duplicate alleles within a locus and exact sequences shared
  between loci using cached sequence hashes.
- [ ] Report distinct nucleotide alleles and distinct translated proteins.
- [ ] Check missing, duplicate, and malformed allele identifiers.
- [ ] Validate representative sequences in `short/` without counting them as
  ordinary loci.
- [ ] Add representative-versus-allele BSR coverage and inter-locus high-BSR
  groups.
- [ ] Support comparison against a second schema (Schema Refinery
  `MatchSchemas`-style analysis).
- [ ] Support optional annotation imports while keeping FASTA-only analysis
  usable without annotations.
- [ ] Cache threshold-independent intermediates by archive fingerprint and
  analysis version so changing display thresholds does not rescan FASTA files.

#### Literature-derived screening criteria

Implement these as configurable inspection flags with the cited value shown in
the interface. They are study-specific criteria, not universal acceptance
thresholds.

| Criterion | Literature screening value | Planned use | Reference |
| --- | --- | --- | --- |
| Potential paralogy | Greater than 70% amino-acid identity between loci | Flag candidate locus pairs | Guglielmini et al. (2021) |
| Repeated homologous regions | At least 90% nucleotide identity across more than 100 bp | Flag repeated inter-locus regions | Pires et al. (2024) |
| Alignment gaps | More than 5% gaps in amino-acid alignments | Flag affected loci | Guglielmini et al. (2021) |
| Alignment gaps | Gaps longer than 6 nucleotides | Flag affected alleles and loci | Savin et al. (2019) |
| Minimum gene length | 50 bp | Flag short loci or alleles | Pires et al. (2024) |
| Sequence similarity | BSR of at least 0.6 | Identify sequences that may belong to the same gene; do not treat this as a definitive paralogy threshold | Silva et al. (2018) |

### Phase 2 — functional scheme performance (genomes required)

- [ ] Define representative genome manifests and record assembly quality and
  lineage metadata.
- [ ] Calculate core-locus prevalence per locus. Offer the literature-used 95%
  prevalence value as an optional reference, clearly distinguished from
  recovering 95% of loci in one genome.
- [ ] Calculate target recovery/call rate per genome.
- [ ] Report missing calls per genome, locus, and lineage.
- [ ] Measure multicopy-hit frequency per locus to provide genome-level
  evidence for potential paralogy.
- [ ] Measure allele-call reproducibility across repeated analyses and, where
  available, different sequencing depths.
- [ ] Compare calls from Juno, chewBBACA, SeqSphere, and the scheme's source
  implementation, including discordant calls and allele distances.
- [ ] Evaluate phylogenetic congruence against an appropriate reference tree.
- [ ] Assess discriminatory power using allelic-distance distributions.
- [ ] Assess epidemiological concordance using supported outbreak or
  transmission groups, including false grouping.
- [ ] Add Schema Refinery spurious-gene classifications when the required
  genome allele-call data are available.

### Phase 3 — typing performance and epidemiological validation

- [ ] Compare Juno with each source implementation on a versioned validation
  genome set.
- [ ] Investigate allele-call and allelic-distance discrepancies.
- [ ] Evaluate recovery of epidemiologically supported clusters and document
  the limits of the available ground truth.

### Scheme documentation and compatibility

- [ ] Record taxonomic scope and intended application.
- [ ] Record the development dataset and its representation of genetic
  diversity.
- [ ] Record development-assembly quality criteria such as N50, contig count,
  and completeness.
- [ ] Record scheme source, version, retrieval date, and allele-database
  snapshot.
- [ ] Validate Juno-typing inputs: FASTA files, sequence-type profile table,
  matching locus identifiers, and KMA indices.
- [ ] Validate Juno-cgMLST inputs: allele FASTA files, Prodigal training file,
  and organism configuration.
