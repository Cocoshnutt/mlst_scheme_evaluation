# mlst_scheme_evaluation

## Project overview

This project develops a standardized, reproducible framework and checklist for
evaluating publicly available MLST and cgMLST schemes. It initially focuses on
schemes for *Yersinia* and *Corynebacterium*. Evaluation criteria derived from
the literature review are implemented incrementally to measure intrinsic scheme
quality, including allelic diversity, sequence integrity, length variation,
polymorphic sites, alignment characteristics, and inter-locus sequence
similarity.

The next stage will use representative genome datasets to assess functional
performance, including locus recovery, allele-call completeness,
reproducibility, and discrepancies between typing tools and pipelines. The
methods will then be validated with schemes from additional organisms, such as
*Salmonella* and *Escherichia coli*, to test their broader applicability and
identify organism- or scheme-specific limitations. The intended outputs are a
final report and a draft validation report that support the assessment of
existing schemes, their integration into the Juno typing pipelines, and the
validation of newly developed schemes.

### Implemented so far

- Archive extraction and versioned caching for the included MLST and cgMLST
  scheme collection.
- Scheme-wide and per-locus allele counts, length distributions, modal-length
  outliers and coding-sequence integrity checks.
- A progressive scheme-wide overview, interactive summary statistics, charts,
  and a searchable per-locus table in the Marimo application.
- Cached nucleotide and protein MAFFT analyses for alignment-gap statistics and
  nucleotide polymorphic-site measurements.
- Cached inter-locus BLAST analyses reporting nucleotide and protein identity,
  coverage, and BLAST Score Ratio measurements.
- Cached representative-to-allele protein BSR coverage and inter-locus
  high-BSR similarity groups.
- A separate Phase 2 functional-performance notebook foundation with a typed,
  machine-readable representative-genome manifest and explicit readiness
  validation.

Genome-based functional evaluation, cross-organism validation, standardized
screening flags, and the final reporting deliverables remain planned work; see
the roadmap below for the detailed status.

## Setup without Nix

The evaluation runs on Linux and macOS, and on Windows through
[WSL](https://learn.microsoft.com/windows/wsl/install). Native Windows is not
currently supported because the analysis calls MAFFT and NCBI BLAST+ command-line
tools.

Install these prerequisites:

- [uv](https://docs.astral.sh/uv/getting-started/installation/) for Python and
  Python package management
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) for within-locus alignments
- [NCBI BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) for cross-locus
  nucleotide and protein similarity

On Ubuntu, Debian, or an Ubuntu-based WSL distribution, install the system tools
with:

```console
sudo apt update
sudo apt install mafft ncbi-blast+
```

On macOS with Homebrew, install them with:

```console
brew install mafft blast
```

For another Linux distribution, install packages providing the `mafft`,
`blastn`, `blastp`, and `makeblastdb` executables. Confirm that all prerequisites
are available on `PATH`:

```console
uv --version
mafft --version
blastn -version
blastp -version
makeblastdb -version
```

From the repository root, create the project environment from the committed
lock file and run the evaluation application:

```console
uv sync --locked
uv run --locked marimo run notebooks/intrinsic.py
```

`marimo run` serves the evaluation as an application without exposing the
notebook editor.

The genome-dependent Phase 2 readiness notebook is a separate application:

```console
uv run --locked marimo run notebooks/functional.py
```

It currently validates manifests from `metadata/genome_sets/` and identifies
missing inputs. It does not yet run typing tools or calculate functional
performance metrics. Complete a dataset-specific copy of
`metadata/genome_sets/template.toml` before implementing those calculations;
see [`docs/functional-workflow.md`](docs/functional-workflow.md). Developer
steps for extending either phase are in
[`docs/adding-an-analysis.md`](docs/adding-an-analysis.md).

After selecting a scheme, static summary, locus-statistics, and allele-table
results load automatically. MAFFT and BLAST analyses run only when requested
with their section controls. A single **Run complete intrinsic evaluation**
control requests all currently implemented expensive intrinsic analyses, and
compatible cached results are reused. The overview at the beginning updates
progressively as requested stages finish; pending overview sections do not
trigger computation.

The expensive-analysis panel also provides one shared worker budget for MAFFT,
BLASTn, and BLASTp. Its default is 80% of the logical CPUs available to the
Marimo process, rounded down to an even number when possible and never below
one. This leaves some capacity for the operating system and behaves safely on
small or CPU-constrained machines. The control can be reduced for shared or
memory-limited systems, or increased to all process-visible CPUs when maximum
throughput is appropriate.

The command-line interfaces use the same default. Pass `--workers COUNT` to
`scheme.alignment` or `--threads COUNT` to `scheme.similarity` and `scheme.bsr`
for an explicit positive override. Explicit CLI and library values are not
silently capped. Worker counts affect scheduling only, so they are deliberately
excluded from scientific cache keys and changing them reuses compatible
results.

uv automatically obtains a compatible Python 3.14 interpreter when one is not
already installed. `uv sync --locked` creates `.venv` and installs the exact
Python dependencies recorded in `uv.lock`; it fails instead of changing the
lock file when `pyproject.toml` and `uv.lock` disagree.

If an analysis reports that MAFFT or a BLAST executable was not found, rerun the
version commands above in the same shell used to start Marimo. The system-tool
installation directory must be on `PATH` before Marimo starts.

## Conda setup

Conda users can create a complete environment from the committed generated
specification. Python and pip come from conda-forge, MAFFT comes from
conda-forge, BLAST+ comes from Bioconda, and the remaining Python packages use
the exact versions exported from `uv.lock`:

```console
conda env create --file environment.yml
conda activate mlst-scheme-evaluation
marimo run notebooks/intrinsic.py
```

After pulling dependency changes, update an existing environment with:

```console
conda env update \
  --name mlst-scheme-evaluation \
  --file environment.yml \
  --prune
```

`environment.yml` is generated; do not edit it directly. Entering
`devenv shell` synchronizes it after uv validates the lock file, but generation
does not require Nix. Run the shared command manually after changing
`pyproject.toml` or `uv.lock`:

```console
uv run --locked python scripts/sync_conda_environment.py
```

To check freshness without changing the file, add `--check`. Commit an updated
`environment.yml` together with its dependency-source changes. The generator
does not create, update, activate, or remove anyone's installed Conda
environment.

## Optional Nix setup

The existing devenv configuration remains available as an alternative. It
provides Python, uv, MAFFT, and BLAST+, synchronizes the Python environment, and
activates it:

```console
devenv shell
marimo run notebooks/intrinsic.py
```

Nix is not required when the setup in the previous section is used.

## Dependency maintenance

Python dependencies are declared in `pyproject.toml` and locked in `uv.lock`.
Use `uv add <package>` to add a dependency, or
`uv lock --upgrade-package <package>` to update one. MAFFT and BLAST+ are system
dependencies; keep the non-Nix installation instructions and `devenv.nix` in
sync if those requirements change. Regenerate `environment.yml` after changing
`pyproject.toml` or `uv.lock`; devenv performs this step automatically.

On a machine without a usable NVIDIA or AMD GPU, marimo 0.25.x may emit a
warning that it cannot collect GPU statistics. It is harmless. To suppress
that upstream warning when desired, start marimo with:

```console
uv run --locked marimo --log-level error run notebooks/intrinsic.py
```

Cache and extract a scheme archive from a notebook cell with:

```python
from scheme.cache import cache_archive, clear_cache
from scheme.info import get_scheme_information

scheme_directory = cache_archive("archive/mlst/corynebacterium.tar.zst")
scheme_information = get_scheme_information(
    scheme_directory,
    minimum_sequence_length=300,
)
```

Calling `cache_archive` again reuses the extracted directory. Call
`clear_cache()` to remove all cached schemes. Scheme statistics and allele
lengths are cached inside each extracted scheme, so repeated evaluations and
new minimum-length thresholds do not reread the FASTA files. Clearing the
scheme cache also clears these calculated results.

On a cold evaluation, scheme information, locus statistics, and per-locus
allele information now come from one source-signature-aware static analysis.
The existing APIs and compatibility cache formats remain available, while the
shared traversal avoids validating and summarizing every FASTA record three
times. The measurement and migration rationale are documented in
[`docs/static-analysis-cache.md`](docs/static-analysis-cache.md).

Analysis cache filenames now use the official `scheme` terminology. Existing
`.schema-*-cache.json` files remain readable, while newly calculated results
are written as `.scheme-*-cache.json` files.

The scheme information module can also be used from the command line:

```console
uv run --locked python -m scheme.info PATH_TO_EXTRACTED_SCHEME
uv run --locked python -m scheme.info PATH_TO_EXTRACTED_SCHEME --ml 300 --json
```

The standard output includes scheme-wide counts and percentages for alleles
below 80% and above 120% of their respective locus's modal length.

Selecting a scheme in `notebooks/intrinsic.py` displays four locus-statistics
panels in a clearly separated section beneath the scheme summary.
`scheme/stats.py` provides the reusable calculations and Altair chart
construction. Its results are cached inside the extracted scheme.

The notebook also displays a searchable, paginated allele-analysis table built
by `scheme/allele_info.py`. The table follows the chewBBACA per-locus fields and
reuses cached static metrics when the minimum allele length changes. The module
can also emit the table data as JSON:

```console
uv run --locked python -m scheme.allele_info PATH_TO_EXTRACTED_SCHEME --ml 201
```

The remaining command-line modules retain their interfaces through the
package entry points:

```console
uv run --locked python -m scheme.cache --help
uv run --locked python -m scheme.alignment --help
uv run --locked python -m scheme.bsr --help
uv run --locked python -m scheme.similarity --help
```

`scheme.stats` remains a library-only chart and statistics module.

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
- [x] Align alleles within each locus and report nucleotide and amino-acid gap
  proportions.
- [x] Calculate the number and percentage of polymorphic nucleotide sites per
  locus.
- [x] Measure inter-locus nucleotide and amino-acid similarity, including
  identity, coverage, and BLAST Score Ratio (BSR).
- [x] Report highly similar inter-locus pairs as **candidate paralogous loci**.
  Sequence similarity alone must not be presented as proof of copy number or
  evolutionary paralogy.

#### Schema Refinery statistics integration

- [x] Add mean allele length, all tied modes, frequent length peaks, and clear
  single-mode/multimodal/outlier locus categories.
- [x] Add per-locus and scheme-wide GC content.
- [x] Report distinct nucleotide alleles and distinct translated proteins.
- [x] Add representative-versus-allele BSR coverage and inter-locus high-BSR
  groups.
- [x] Cache threshold-independent intermediates for implemented analyses by
  archive fingerprint and analysis version so changing display thresholds does
  not rescan FASTA files. Future analyses must preserve this behavior.

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

Implemented foundation:

- [x] Keep Phase 2 in `notebooks/functional.py`, separate from the intrinsic
  notebook while sharing stable repository conventions.
- [x] Define a typed representative-genome manifest, machine-readable JSON
  Schema, template, SHA-256 identity, and readiness validation.
- [x] Document independent cache boundaries, missing-value requirements, and
  FAIR result-artifact expectations.

- [ ] Define representative genome manifests and record assembly quality and
  lineage metadata for the real evaluation collections.
- [ ] Calculate core-locus prevalence per locus. Offer the literature-used 95%
  prevalence value as an optional reference, clearly distinguished from
  recovering 95% of loci in one genome.
- [ ] Calculate allele-call completeness, including target recovery/call rate
  per genome and missing-call patterns by genome, locus, and lineage.
- [ ] Measure multicopy-hit frequency per locus to provide genome-level
  evidence for potential paralogy.
- [ ] Measure allele-call reproducibility across repeated analyses and, where
  available, different sequencing depths.
- [ ] Compare calls from Juno, chewBBACA, SeqSphere, and the scheme's source
  implementation, including discordant calls and allele distances.
- [ ] Add Schema Refinery spurious-gene classifications when the required
  genome allele-call data are available.

#### Extended validation analyses

- [ ] Evaluate phylogenetic congruence against an appropriate reference tree.
- [ ] Assess discriminatory power using allelic-distance distributions.
- [ ] Assess epidemiological concordance using supported outbreak or
  transmission groups, including false grouping.

### Phase 3 — cross-organism and typing-pipeline validation

- [ ] Compare Juno with each source implementation on a versioned validation
  genome set.
- [ ] Investigate allele-call and allelic-distance discrepancies.
- [ ] Evaluate recovery of epidemiologically supported clusters and document
  the limits of the available ground truth.
- [ ] Apply the framework to additional MLST and cgMLST schemes, including
  schemes for *Salmonella* and *Escherichia coli*.
- [ ] Compare findings across organisms and document organism- or
  scheme-specific limitations in the evaluation criteria.

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

### Documentation and deliverables

- [ ] Produce a versioned, standardized evaluation checklist that records each
  criterion, its required inputs, calculation, literature basis, and intended
  interpretation.
- [ ] Complete every pending field in `metadata/scheme_sources.toml` after the
  refactor and before publishing a release or research output.
- [ ] Document the evaluation methodologies and findings in the final report.
- [ ] Prepare a draft validation report covering cross-organism and
  typing-pipeline validation.

## License

The software and project documentation are licensed under the European Union
Public Licence, version 1.2 (`EUPL-1.2`). See [LICENSE](LICENSE).

The scheme archives under `archive/` originate from external data providers and
are not relicensed by this project. Their source, licensing, and transformation
details must be recorded separately in `metadata/scheme_sources.toml` before a
release or research output is published.
