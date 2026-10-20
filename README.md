# sporos

Fast DNA seeding in [Mojo](https://mojolang.org): (w, k)-minimizers, open and closed syncmers,
maximal exact matches (MEMs) and maximal unique matches (MUMs).

Every kernel is tested against a brute-force oracle, and `sporos mem` reproduces MUMmer's
`mummer` output exactly (same match sets; `-mum` and `-mumreference` byte-identical).

## Performance

On xaros (Intel i9-7900X, AVX-512), against each tool's best settings; the full method, every table and
the correctness checks (60/60 match sets equal to MUMmer and copMEM2) are in
[`benchmarks/competitors`](benchmarks/competitors/README.md).

| Seeds, CHM13 chr1, 1 thread, ns/bp | sporos | simd-minimizers 3.0 | minimap2 `mm_sketch` |
|---|---|---|---|
| minimizers k=21 w=11, forward / canonical | **1.18 / 1.11** | 1.81 / 2.17 | — / 17.1 |
| closed syncmers k=31 s=21, forward / canonical | **1.09 / 1.25** | 2.00 / 2.58 | |
| minimizers k=31 w=19, forward / canonical | **1.04 / 0.98** | 1.64 / 1.93 | |

| MEMs, best wall s (1 thread / 10 threads) | sporos | copMEM2 `-mf` | MUMmer4 |
|---|---|---|---|
| hg38 chr21 vs chimp chr21, L=50 | **1.08 / 0.47** | 1.96 / 1.95 | 17.4 |
| hg38 chr1 vs CHM13 chr1, L=100 | **8.95 / 2.72** | 11.4 / 11.1 | 132 |
| hg38 chr1 vs CHM13 chr1, L=50 | **36.8 / 7.6** | 46.9 / 46.5 | 159 |

Seed-index lookups (`SeedIndex.lookup_batch`, hg38 chr21 index, k=19 w=12, occurrence cap 20, simulated
reads): 117 Mbp/s per thread and 560 Mbp/s on 8 threads for read minimizers on 150 bp reads, and 37 Mbp/s per
thread for every read k-mer (GraphAligner's mode). `sporos seeds` on 1 M 150 bp reads takes 1.2 s end to end
with 10 threads, including index build and 27 M output lines.

Known gaps: with 10 threads simd-minimizers' chunked harness is 1.1-1.8x faster on seeds, and `-mum` on
repeat-dense input (chr1, L=50) spends most of its time filtering MEMs to MUMs.

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
build/sporos syncmers -k 21 -s 11 [--open 5[,8...]] seqs.fa   # closed, or open at one or more offsets
build/sporos seeds -k 19 -w 12 ref.fa reads.fq.gz           # minimizer seed index + lookups
build/sporos seeds -k 20 --syncmers -s 16 [--open 2] ref.fa reads.fq
```

## Library

```mojo
from sporos.minimizers import minimizer_positions
from sporos.syncmers import syncmer_positions

var seq = "ACGTTGCA...".as_bytes()
var mins = minimizer_positions(seq, k=21, w=11)                  # canonical by default
var syncs = syncmer_positions(seq, k=21, s=11)                   # closed syncmers
var opens = syncmer_positions(seq, k=21, s=11, open_offset=5)    # open syncmers
var multi = syncmer_positions(seq, k=21, s=11, open_offsets=[2, 5, 8])  # open at any of several offsets
```

Syncmers come in both kinds everywhere (library, batches, seed index, CLI, C, Python): a k-mer is a
**closed** syncmer when its smallest s-mer starts at offset 0 or k - s, and an **open** syncmer when it
starts at offset `open_offset` (0-based; strobealign's k=20, s=16 is `open_offset=2`), or at any of a set of
offsets (`open_offsets`; `{0, k - s}` gives the closed syncmers). One offset is a single vector compare; a set
is a bitmask test.

## Develop

Requires [pixi](https://pixi.sh). `pixi run test` builds and runs the test suites, `pixi run bench-seeds`
times the seed kernels, `pixi run format` formats the sources. Coding agents: see `CLAUDE.md` (the
Modular Mojo skills are vendored in `.claude/skills` and injected by hooks).

## C and Python

```sh
pixi run build-c        # build/c: libsporos, the Mojo runtime libraries it loads, and c/sporos.h
pixi run build-python   # build/python/sporos: the package and its Mojo extension (Python >= 3.10)
pixi run test-c && pixi run test-python
```

Both folders run without pixi or Mojo. From C (`-Lbuild/c -lsporos`, build/c on the rpath):

```c
#include "sporos.h"
uint32_t pos[1000];
int64_t n = sporos_minimizers(seq, len, 21, 11, /*canonical*/ 1, /*threads*/ 1, pos, 1000);  // len >= n
sporos_index *index = sporos_index_new(refs, ref_lens, n_refs, /*min_len*/ 50, /*threads*/ 0);
sporos_matches *mums = sporos_index_matches(index, query, qlen, SPOROS_MUM, /*reverse*/ 0, 0);
const sporos_match *m = sporos_matches_data(mums);  // ref_seq, ref_pos, query_pos, length
sporos_matches_free(mums);
sporos_index_free(index);
```

From Python (`PYTHONPATH=build/python`):

```python
import sporos
sporos.minimizers("ACGT...", k=21, w=11)               # array('I') of positions
sporos.syncmers("ACGT...", k=21, s=11, open_offset=5)          # or open_offset=[2, 5, 8]
index = sporos.Index(["ACGT...", "GGCA..."], min_len=50, threads=8)
index.matches(query, mode="mum", reverse=False)        # [(ref_index, ref_pos, query_pos, length), ...]
```
