# Shared static analysis cache decision

## Decision

Scheme information, locus statistics, and per-locus allele information share
one threshold-independent FASTA traversal and versioned cache. Their existing
public APIs, result dataclasses, CLI output, and compatibility cache files are
retained as projections.

The shared cache is `.scheme-static-analysis-cache.json`, owned by
`scheme/static_analysis.py`. It is invalidated by its schema version and the
ordered relative path, size, and modification time of every recognized FASTA
file. Writes are atomic. An unreadable, corrupt, incompatible, or unwritable
cache is treated as a non-fatal cache miss.

Minimum allele length and frequent-length-peak settings are deliberately not
part of the cache key. They continue to reclassify stored lengths and exact
length frequencies without rereading FASTA files.

## Measurement

The change was justified by measuring the actual static notebook sequence:

1. Extract an archive into an isolated temporary cache.
2. Call `get_scheme_information()`.
3. Call `get_scheme_statistics()`.
4. Call `get_scheme_allele_information()`.
5. Repeat the three calls to measure warm-cache behavior.

No repository cache was read, cleared, or modified, and MAFFT and BLAST were
not invoked. Each measurement below is one indicative run on 2026-11-18 under
WSL2 Linux with Python 3.14.8 and an Intel Core Ultra 7 265H exposed as 16
logical CPUs. Operating-system file caching and other machine load were not
controlled, so these values are engineering evidence rather than a formal
performance benchmark.

| Scheme archive | Loci | Alleles | Before, cold | Shared cache, cold | Before, warm | Shared cache, warm |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `mlst/corynebacterium.tar.zst` | 7 | 905 | 0.029 s | 0.019 s | 0.0007 s | 0.0007 s |
| `cgmlst/yersinia_cgmlst_alleles.tar.zst` | 500 | 180,652 | 17.05 s | 12.19 s | 0.046 s | 0.045 s |
| `cgmlst/yersinia_enterobase_cgmlst_v1_hiercc.tar.zst` | 1,553 | 684,276 | 72.34 s | 47.86 s | 0.179 s | 0.207 s |

The largest measured cold startup improved by approximately 24.5 seconds, or
34%. This was considered material enough to accept the cache migration. Warm
behavior remains effectively unchanged.

## Compatibility and scientific boundaries

- Existing `.scheme-information-cache.json`,
  `.scheme-statistics-cache.json`, and
  `.scheme-allele-information-cache.json` files are still read first.
- A new shared analysis is projected back into those established cache formats
  without changing their versions or meanings.
- The shared calculation reuses the existing FASTA parser, coding-sequence
  validator, translation function, base-composition definition, and locus
  statistics calculator.
- Missingness, invalidity-category overlap, modal-length ties, strict 80/120%
  boundaries, GC denominators, and distinct-protein rules are unchanged.
- The shared cache is an implementation detail, not a FAIR result artifact.
  Published evaluation outputs still require explicit, documented exports.
