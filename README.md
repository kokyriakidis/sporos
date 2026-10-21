# sporos

Fast DNA seeding in [Mojo](https://mojolang.org): (w, k)-minimizers, open and closed syncmers, k-mer seed
indexes, maximal exact matches (MEMs) and maximal unique matches (MUMs). sporos ships as a Mojo library, a
command-line tool, a C library and a Python package.

- **Fast.** On one thread, minimizers and syncmers run 1.3–2.4x faster than simd-minimizers. On hg38 chr1 vs
  CHM13 chr1, MEMs are found 2x faster than copMEM2 (8–12x with 10 threads), and MUMs 4–18x faster than MUMmer
  on one thread. See [Performance](#performance).
- **Exact.** Every kernel is tested against a brute-force oracle. `sporos mem` gives MUMmer's `mummer` output:
  the same match sets, and `-mum` and `-mumreference` byte for byte. Results do not depend on the thread count.
- **Practical.** Reads FASTA and FASTQ (plain, gzip or block-parallel BGZF) from files or pipes. Multithreaded
  throughout, with no dependencies beyond the Mojo runtime and libdeflate.

## What it computes

| | Definition |
|---|---|
| **Minimizers** | For each window of `w` consecutive k-mers, the start of its smallest k-mer under a 32-bit ntHash, ties broken leftmost. Canonical by default (forward and reverse-complement hashes combined). Each position is reported once. |
| **Closed syncmers** | k-mers whose smallest s-mer starts at offset 0 or `k - s`. |
| **Open syncmers** | k-mers whose smallest s-mer starts at a chosen offset `t`, or at any offset of a set. strobealign's default is `k=20, s=16, t=2`. |
| **Seed index** | Reference minimizers, syncmers or all k-mers (k ≤ 32), keyed by canonical k-mer so lookups find both strands. Hits carry occurrence counts, with an optional abundance cutoff (GraphAligner-style `--keep`). |
| **MEMs** | Every maximal exact match of length ≥ L between a query sequence and the reference (MUMmer's `-maxmatch`). |
| **MUMs** | MEMs whose string occurs exactly once in the reference and once in the query sequence (MUMmer's `-mum`). `-mumreference` requires uniqueness in the reference only. |

MEMs and MUMs are **pairwise**, as in MUMmer. The reference can hold many sequences, indexed together, and
each query record is matched against all of them on the forward strand, the reverse strand or both.
Multi-genome MUMs, unique in each of N genomes at once as in mumemto or Mumento2, are not supported.

Non-ACGT bases (`N` and other IUPAC codes) never match and break k-mers. Lowercase is accepted, and `U` reads
as `T`. Positions are 0-based in the Mojo, C and Python APIs and in the `minimizers`, `syncmers` and `seeds`
output. `sporos mem` prints MUMmer's 1-based coordinates.

## Install

sporos builds with [pixi](https://pixi.sh), which fetches the Mojo toolchain (nightly channel) and libdeflate:

```sh
git clone https://github.com/kokyriakidis/sporos && cd sporos
pixi run build-cli        # builds build/sporos
pixi run test             # optional: checks every kernel against its brute-force oracle
```

pixi targets Linux x86-64, Linux arm64 and macOS on Apple silicon. sporos is tested on Linux x86-64 with
AVX-512, including 8-lane (AVX2-width) builds, and on macOS arm64. Linux arm64 is untested. The C library and
the Python package build into self-contained folders that run without pixi or Mojo (see
[C and Python](#c-and-python)).

## Command line

```sh
sporos minimizers -k 21 -w 11 -t 8 genome.fa            # name<TAB>position, one line per minimizer
sporos syncmers -k 21 -s 11 genome.fa                   # closed syncmers
sporos syncmers -k 20 -s 16 --open 2 genome.fa          # open syncmers (--open 1,3 for a set)
sporos seeds -k 19 -w 12 --keep 0.9998 ref.fa reads.fq.gz
sporos mem -mum -l 50 -t 8 reference.fa query.fa        # MUMmer's flags and output format
sporos mem -maxmatch -b -l 100 reference.fa query.fa    # all MEMs, both strands
```

Every command takes `-t` for the thread count (default: every physical core). `seeds` and `mem` also take
`-v`, which prints phase timings and peak memory to stderr. Any input can be gzipped, BGZF-compressed, or a
pipe (`/dev/stdin`).

### `minimizers` and `syncmers`

```
sporos minimizers -k K -w W [--forward] [--count] [-t threads] <seqs.fa>
sporos syncmers   -k K -s S [--open T[,T...]] [--forward] [--count] [-t threads] <seqs.fa>
```

| Option | Meaning |
|---|---|
| `-k`, `-w` | k-mer length (1 to 64, default 15) and window length in k-mers (default 10) |
| `-s` | s-mer length for syncmers (`1 ≤ s ≤ k`) |
| `--open T[,T...]` | open syncmers at offset `T`, or at any offset in the list (each in `[0, k - s]` and below 64) |
| `--forward` | forward-strand hashes instead of canonical ones |
| `--count` | print only the totals |

Each output line is `name<TAB>position`, 0-based and increasing within each record. Files with many records,
such as reads, are scanned with one read per SIMD lane.

### `seeds`

```
sporos seeds -k K [-w W | --syncmers -s S [--open T[,T...]] | --all-kmers] [--query-kmers]
             [--max-occ N | --keep F] [-t threads] [-v] <reference.fa> <reads.fa>
```

`seeds` builds a seed index of the reference (minimizers by default; k ≤ 32, default k=19, w=12) and looks up
every read's seeds, or every read k-mer with `--query-kmers`. `--max-occ N` skips seeds that occur more than
`N` times in the reference, and `--keep F` sets that cutoff at the `F` quantile of the indexed k-mers. Each
output line is one hit:

```
read  query_pos  ref_name  ref_pos  strand  occurrences
```

`strand` is `-` when the read k-mer matches the reverse complement of the reference.

### `mem`

```
sporos mem [-mum | -mumreference | -maxmatch] [-l L] [-b | -r] [-F] [-n] [-t threads] [-v]
           [-K anchor_k] [-k1 ref_step] [-k2 query_step] <reference.fa> <query.fa>
```

| Option | Meaning |
|---|---|
| `-maxmatch`, `-mum`, `-mumreference` | all MEMs, MUMs, or MEMs unique in the reference (the default, as in MUMmer) |
| `-l L` | minimum match length (default 20) |
| `-b`, `-r` | both strands, or the reverse complement only (default: forward) |
| `-F` | four-column output, reference name first (always on with several reference sequences) |
| `-n` | accepted for MUMmer compatibility; non-ACGT bases never match anyway |
| `-K`, `-k1`, `-k2` | anchor length and the coprime reference and query sampling steps (chosen automatically) |

The output follows `mummer`: a `> query` header (`> query Reverse` for the reverse strand), then one
`ref_pos query_pos length` line per match, 1-based.

## Library

Put the repository on Mojo's include path (`mojo build -I path/to/sporos ...`):

```mojo
from sporos.batch import minimizers_batch
from sporos.fasta import read_fasta
from sporos.index import MINIMIZER_SEEDS, SeedIndex
from sporos.mems import MemIndex, MemParams, filter_mums, find_mems
from sporos.minimizers import minimizer_positions
from sporos.syncmers import syncmer_positions
from sporos.text import SeqText


def main() raises:
    var genome = read_fasta("ref.fa", threads=4)  # FASTA/FASTQ, plain, gzip or BGZF
    var seq = genome.seq(0)

    # Seeds: k-mer start positions, increasing.
    var mins = minimizer_positions(seq, k=21, w=11, threads=4)  # canonical by default
    var closed = syncmer_positions(seq, k=21, s=11)  # closed syncmers
    var opened = syncmer_positions(seq, k=20, s=16, open_offset=2)  # open; or open_offsets=[1, 3]

    # Many short reads: one read per SIMD lane.
    var reads = read_fasta("reads.fq")
    var spans = List[Span[UInt8, origin_of(reads.data)]]()
    for r in range(len(reads)):
        spans.append(reads.seq(r))
    var batch = minimizers_batch(spans, k=15, w=10, threads=4)  # batch.seeds(i): read i's positions

    # A seed index of the reference, probed with every read's minimizers.
    var refs = List[Span[UInt8, origin_of(genome.data)]]()
    refs.append(seq)
    var index = SeedIndex(refs, k=19, w=12, scheme=MINIMIZER_SEEDS, threads=4)
    var hits = index.lookup_batch(spans, max_occurrences=index.occurrence_cutoff(0.9998), threads=4)
    for h in hits.of(0):  # read 0's hits
        print(h.query_pos, h.ref_pos, h.reverse, h.occurrences)

    # MEMs and MUMs of a query against the reference.
    var rtext = SeqText(refs)
    var qseqs = List[Span[UInt8, origin_of(reads.data)]]()
    qseqs.append(reads.seq(0))
    var qtext = SeqText(qseqs)
    var mem_index = MemIndex(rtext, MemParams.choose(30, len(rtext), len(qtext)), threads=4)
    var mems = find_mems(mem_index, rtext, qtext, threads=4)  # every MEM of length >= 30
    var mums = filter_mums(mems, qtext, threads=4)  # MUMmer's -mum
    for m in mums:
        print(m.ref_pos, m.query_pos, m.length)
```

Related calls:
- `syncmers_batch` is the batch counterpart for syncmers.
- `SeedIndex.lookup` probes a single read.
- `filter_mums(..., unique_in_query=False)` gives `-mumreference`.

The docstrings in `sporos/` describe every parameter and its limits.

## C and Python

```sh
pixi run build-c        # build/c: libsporos, the Mojo runtime libraries it loads, and sporos.h
pixi run build-python   # build/python/sporos: the package and its Mojo extension (Python >= 3.10)
pixi run test-c && pixi run test-python
```

C (`cc -I build/c app.c -L build/c -lsporos -Wl,-rpath,build/c`; the API is in [`c/sporos.h`](c/sporos.h)):

```c
#include "sporos.h"

uint32_t positions[1024];
int64_t n = sporos_minimizers(seq, len, 21, 11, /*canonical*/ 1, /*threads*/ 1, positions, 1024);
// n is the full count even when it exceeds the room given; n <= len always.

const uint8_t *refs[] = {ref};
int64_t lens[] = {ref_len};
sporos_index *index = sporos_index_new(refs, lens, 1, /*min_len*/ 50, /*threads, 0 = all cores*/ 0);
sporos_matches *mums = sporos_index_matches(index, query, query_len, SPOROS_MUM, /*reverse*/ 0, 0);
const sporos_match *m = sporos_matches_data(mums);  // ref_seq, ref_pos, query_pos, length
for (int64_t i = 0; i < sporos_matches_count(mums); i++) { /* ... */ }
sporos_matches_free(mums);
sporos_index_free(index);
```

The C API also has `sporos_syncmers` (closed, or open at one offset) and `sporos_syncmers_open` (open at a set
of offsets). Invalid arguments return `SPOROS_INVALID_ARGUMENT` or `NULL`.

Python (`PYTHONPATH=build/python`):

```python
import sporos

sporos.minimizers(seq, k=21, w=11)                       # array('I') of positions
sporos.syncmers(seq, k=20, s=16, open_offset=2)          # open_offset: None (closed), an int, or a list
index = sporos.Index([ref1, ref2], min_len=50, threads=8)
index.matches(query, mode="mum", reverse=False)          # [(ref_index, ref_pos, query_pos, length), ...]
```

The modes are `"maxmatch"`, `"mum"` and `"mumreference"`. Sequences can be `str` or `bytes`. Invalid
arguments raise `ValueError` or `TypeError`.

## Performance

Measured on an Intel i9-7900X (Skylake-X, AVX-512, 10 cores, turbo off) against each tool's best settings, as
best wall time including input and output. The full method, competitor versions, every table and the
correctness checks are in [`benchmarks/competitors`](benchmarks/competitors/README.md); all 60 match sets equal
MUMmer's and copMEM2's. Competitor rows come from that suite; sporos rows are from the current code.

| MEMs (1 / 10 threads) | sporos | copMEM2 `-mf` | MUMmer4 (1 thread) |
|---|---|---|---|
| hg38 vs CHM13 chr1, L=100 | **6.0 s / 1.4 s** | 11.4 s / 11.1 s | 132 s |
| hg38 vs CHM13 chr1, L=50 (65.6 M MEMs) | **22.5 s / 3.9 s** | 46.9 s / 46.5 s | 159 s |
| hg38 vs chimp chr21, L=50 | **0.81 s / 0.24 s** | 1.96 s / 1.95 s | 17.4 s |

| MUMs, `-mum` | sporos (1 / 10 threads) | MUMmer4 (1 thread) |
|---|---|---|
| hg38 vs CHM13 chr1, L=100 | **6.0 s / 1.5 s** | 93.8 s |
| hg38 vs CHM13 chr1, L=50 | **23.2 s / 4.8 s** | 87.5 s |
| hg38 vs chimp chr21, L=50 | **0.82 s / 0.26 s** | 14.6 s |

| Seeds, CHM13 chr1, ns/bp (forward / canonical) | sporos | simd-minimizers 3.0 | minimap2 `mm_sketch` |
|---|---|---|---|
| minimizers k=21 w=11, 1 thread | **1.18 / 1.11** | 1.81 / 2.17 | — / 17.1 |
| minimizers k=31 w=19, 1 thread | **1.04 / 0.98** | 1.64 / 1.93 | |
| closed syncmers k=31 s=21, 1 thread | **1.09 / 1.25** | 2.00 / 2.58 | |
| minimizers k=21 w=11, 10 threads | **0.18** / 0.21 | 0.19 / 0.20 | — / 1.96 |
| minimizers k=31 w=19, 10 threads | **0.14 / 0.16** | 0.16 / 0.19 | |

Notes on the comparison:
- With 10 threads, sporos returns one merged list of positions; simd-minimizers' harness only counts each
  thread's positions.
- MUMmer 4.0.1 runs single-threaded, since its threads do not help with one query record.
- copMEM2 gains nothing from threads on one query record.

**Seed-index lookups** (hg38 chr21 index, k=19, w=12, occurrence cap 20, simulated 150 bp reads):
- read minimizers: 133 Mbp/s per thread, 560 Mbp/s on 8 threads;
- every read k-mer (GraphAligner's mode): 36 Mbp/s per thread.

`sporos seeds` on 1 M reads takes 1.2 s end to end with 10 threads, including the index build and 27 M output
lines.

**Peak memory** for MEMs on chr1 is 1.2 GB at L=100 and about 2.5 GB at L=50; copMEM2 `-mf` uses 1.7 GB.

## How it works

- **Seed kernel** (`sporos/window.mojo`).
  - One engine serves minimizers and syncmers. A sequence splits into 8 SIMD lanes (16 with AVX-512), each
    rolling a 32-bit ntHash.
  - A byte transpose feeds the lanes and precomputes seed-table indices, so each hash update is one byte
    shuffle (NEON `tbl` or x86 `pshufb`).
  - The sliding-window minimum uses prefix/suffix blocks with single u32 keys, which keeps leftmost
    tie-breaking exact. Deduplication is fused into the loop.
  - Batches give each read its own lane.
  - Threads take small segments in order and copy their positions once into the final list, after the
    previous segment publishes where it ends.
- **MEMs** (`sporos/mems.mojo`, `sporos/table.mojo`).
  - The reference is indexed by K-mer anchors every `k1` bases, and the query probes with anchors every `k2`
    bases, with `k1` and `k2` coprime and `k1 · k2 ≤ L - K + 1`.
  - Anchors of one MEM sit `k1 · k2` apart on its diagonal, so the anchor whose match does not reach that far
    to the left reports it: each MEM is found exactly once, with no deduplication pass.
  - The hash table is built by a parallel two-level radix sort, with about four anchors per bucket.
  - Probes stream through a software pipeline: hash, bucket bounds, a vector fingerprint compare, then a
    prefetch of the candidate's reference text before extension, 32 bases per compare.
  - Matches are radix-sorted per query chunk.
- **MUMs.** A MEM is a MUM unless another MEM of the same query contains it in reference or query
  coordinates. Both checks are linear sweeps, parallelized by blocks, after a parallel radix sort by
  reference start.
- **Seed index** (`sporos/index.mojo`).
  - Canonical k-mers are hashed by a bijection on 2k bits. Bucket and fingerprint therefore identify a k-mer
    exactly whenever their bits suffice; otherwise candidates are checked against the text.
  - Batched lookups stream every read through one prefetch pipeline.
- **Input** (`sporos/fasta.mojo`).
  - FASTA is parsed in place.
  - BGZF inflates block by block in parallel with libdeflate.
  - Plain gzip goes through libdeflate, with zlib as a fallback.
- **Threads** come from a small pthreads `parallel_for` (`sporos/parallel.mojo`).

## Limits

- Each sequence, and the concatenated reference, must be shorter than 2^32 bases.
- Minimizers and syncmers take k ≤ 64.
- The seed index and MEM anchors take k ≤ 32.
- Open-syncmer offsets must be below 64.
- The performance numbers above are from x86-64 Linux. Apple-silicon performance has not been re-measured
  since the latest optimizations.

## Develop

| Task | What it does |
|---|---|
| `pixi run test` | builds and runs the Mojo test suites (oracles in `tests/reference.mojo`) |
| `pixi run test-c`, `pixi run test-python` | builds the bindings and checks them against the library's outputs |
| `pixi run bench-seeds`, `bench-reads`, `bench-index` | measures seed kernel, read batch and seed-index throughput |
| `pixi run format` | formats the sources (120 columns) |
| `bash benchmarks/competitors/run_*.sh` | runs the competitor suite (see its README) |

Pass `-D SPOROS_LANES=8` to `mojo build` to force AVX2-width lanes. Coding agents should read
[`CLAUDE.md`](CLAUDE.md): the Modular Mojo skills are vendored in `.claude/skills` and take precedence over
pretrained Mojo habits.

## License

[MIT](LICENSE).
