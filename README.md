# sporos

Fast DNA seeding in [Mojo](https://mojolang.org): (w, k)-minimizers, open and closed syncmers,
maximal exact matches (MEMs) and maximal unique matches (MUMs).

Every kernel is tested against a brute-force oracle, and `sporos mem` reproduces MUMmer's
`mummer` output exactly (same match sets; `-mum` and `-mumreference` byte-identical).

## Performance (preliminary)

Single thread unless noted. 100 Mbp random DNA for seeds; ns per base, lower is better.

| Seeds, k=21 w=11 | sporos (xaros) | simd-minimizers 3.0 (xaros) | sporos (M2) | simd-minimizers (M2) |
|---|---|---|---|---|
| minimizers, forward | **1.14** | 1.79 | **1.13** | 1.49 |
| minimizers, canonical | **1.21** | 2.12 | **1.26** | 1.97 |
| closed syncmers, forward | **1.12** | 1.99 | **1.17** | 1.58 |
| closed syncmers, canonical | **1.20–1.32** | 2.55 | **1.29** | 2.05 |

minimap2's `mm_sketch` takes 7.8 ns/bp on the M2. xaros is an Intel i9-7900X (AVX-512); the M2 numbers
were taken on a loaded machine (thread CPU time, best of 7).

| MEMs, hg38 chr21 vs panTro chr21 (M2) | L=50 | L=100 |
|---|---|---|
| **sporos, 8 threads** | **0.18 s** | **0.17 s** |
| sporos, 1 thread | 0.46 s | 0.35 s |
| copMEM2 `-mf` | 0.85 s | 0.42 s |
| MUMmer4 `mummer -maxmatch -n` | 15.7 s | 12.1 s |

Wall time including reading FASTA and writing output. MUMmer's `-mum` and `-mumreference` take
11.6–12.1 s on the same pair; sporos's take the same time as `-maxmatch`.

## Design

- **Minimizers and syncmers** share one engine (`sporos/window.mojo`). The sequence is split into 8 SIMD
  lanes (16 with AVX-512), each rolling a 32-bit ntHash (rotation 7; canonical = forward + reverse
  complement). Bases reach the lanes through a byte transpose that also precomputes seed-table indices,
  so every seed lookup is one byte shuffle (NEON `tbl` / x86 `pshufb`). The sliding-window minimum uses
  prefix/suffix blocks with single u32 keys (hash bits above a block offset), which keeps leftmost
  tie-breaking exact. Dedup is fused into the loop: one vertical compare per step and a
  count-trailing-zeros walk over survivor bitmasks.
- **MEMs** (`sporos/mems.mojo`) use sampled k-mer anchors: reference K-mers every `k1` positions,
  query K-mers every `k2`, with `k1`, `k2` coprime and `k1 * k2 <= L - K + 1`. Anchors of one MEM sit
  `k1 * k2` apart on its diagonal, so an anchor reports its MEM iff the match does not reach `k1 * k2`
  bases to the left: each MEM is found and extended exactly once, with no deduplication pass. The
  index is built by a parallel two-level radix sort; probes are batched with prefetching.
- **MUMs** are the MEMs not contained, in reference or query coordinates, by another MEM of the same
  query sequence: exactly the matches whose string occurs once in each (MUMmer's `-mum`).
  Query-coordinate containment alone gives `-mumreference`.
- Threads come from a dependency-free pthreads `parallel_for` (`sporos/parallel.mojo`).

## Command line

```sh
pixi run build-cli
build/sporos mem -mum -l 50 -t 8 reference.fa query.fa     # mummer's flags and output format
build/sporos mem -maxmatch -b -l 100 reference.fa query.fa
build/sporos minimizers -k 21 -w 11 -t 8 seqs.fa           # name<TAB>position
build/sporos syncmers -k 21 -s 11 [--open 5] seqs.fa
```

## Library

```mojo
from sporos.minimizers import minimizer_positions
from sporos.syncmers import syncmer_positions

var seq = "ACGTTGCA...".as_bytes()
var mins = minimizer_positions(seq, k=21, w=11)                  # canonical by default
var syncs = syncmer_positions(seq, k=21, s=11)                   # closed syncmers
var opens = syncmer_positions(seq, k=21, s=11, open_offset=5)    # open syncmers
```

## Develop

Requires [pixi](https://pixi.sh). `pixi run test` builds and runs the test suites, `pixi run bench-seeds`
times the seed kernels, `pixi run format` formats the sources. Coding agents: see `CLAUDE.md` (the
Modular Mojo skills are vendored in `.claude/skills` and injected by hooks).
