# Competitor benchmarks

sporos against MUMmer4, copMEM2, mumemto, BWA (`fastmap` SMEMs), ropebwt3 (`mem` SMEMs), simd-minimizers,
minimap2's `mm_sketch` and seed lookups, and strobealign's seed lookups, on one Linux box. Raw results are TSVs in [`results/`](results/); the tables
below come from `python3 -I tables.py`.

## Reproduce

The benchmark box needs `pixi`, `cargo`, `gcc`, `git` and `curl`. From the dev machine, mirror the checkout
(the repo is private, the box has no GitHub credentials), then run everything there:

```sh
rsync -az --exclude .pixi --exclude build ./ xaros.tail2f508.ts.net:sporos/
ssh xaros.tail2f508.ts.net
cd ~/sporos
bash benchmarks/competitors/setup_xaros.sh   # idempotent; tools, harnesses, sporos builds, data, versions.tsv
bash benchmarks/competitors/run_seeds.sh     # ~4 min
bash benchmarks/competitors/run_mems.sh      # ~1 h 45 min; PAIRS="ecoli chr21" for the ~25 min subset
python3 -I benchmarks/competitors/tables.py  # markdown tables from results/*.tsv
```

`setup_xaros.sh` puts everything under `~/sporos-bench`: a pixi env with bioconda `mummer4`, `mumemto`, `bwa`
and `ropebwt3` (bioconda has no linux-64 `e-mem`, so E-MEM is skipped); copMEM2, minimap2 and simd-minimizers
at pinned commits, built with `-O3 -march=native` / `target-cpu=native`; the two seed harnesses in
[`harness/`](harness/); and snapshots of `build/sporos`, `bench_seeds` and
[`bench_seeds_fasta`](bench_seeds_fasta.mojo) for `SPOROS_LANES=8` and `16` in `~/sporos-bench/bin`, so a
rebuild of the checkout cannot change binaries mid-run.

Data (`~/sporos-bench/data`):

| pair | reference | query |
|---|---|---|
| ecoli | E. coli K-12 MG1655 `NC_000913.3` (4.6 Mbp) | E. coli O157:H7 Sakai `NC_002695.2` (5.5 Mbp) |
| chr21 | hg38 chr21 (46.7 Mbp, 6.6 Mbp N) | chimpanzee Clint_PTRv2 chr21 `NC_036900.1` (33.2 Mbp) |
| chr1 | hg38 chr1 (248.9 Mbp, 18.5 Mbp N) | T2T CHM13v2.0 chr1 (248.4 Mbp) |

CHM13 chr1 is also the input for the seed benchmarks.

## Method

- **Timing** (`lib.sh: measure`): `/usr/bin/time -v`; best wall time of 3 runs, with the user/system time of
  that run and the largest max RSS. A run longer than 60 s is not repeated (`reps` column), which only affects
  the chr1 pair. Before every measurement the script samples `/proc/stat` for 0.5 s and waits while other jobs
  keep more than 1.5 CPUs busy; the busy-CPU sample and the 1-minute load average are stored with each row
  (`busy_cpus`, `load1`). The box idles at about 0.5-1 busy CPU (desktop session), and `load1` includes
  the previous benchmark's own threads.
- **MEMs** (`run_mems.sh`), for L = 50 and 100, forward strand:
  - `sporos mem -maxmatch|-mum|-mumreference -l L -n -t T ref query`, T = 1 and 10.
  - `mummer -maxmatch|-mum|-mumreference -n -l L ref query`, single-threaded: in MUMmer 4.0.1 `-qthreads`
    only splits work across query records (no gain with one query; checked) and `-k >1 -threads` aborts
    (`terminate called after throwing an instance of 'char const*'`, empty output).
  - `copmem2 -q [-mf] -t T -l L -o out ref query` (all MEMs; copMEM2 has no MUM modes).
  - `mumemto -l L -o prefix ref query`: strict multi-MUMs, which mumemto computes over **both strands**
    (a match must be unique among both strands of both inputs), single-threaded, building its PFP/SA index
    every run. Not the same set as MUMmer's `-mum`; compared for information only. Capped at
    `MUMEMTO_TIMEOUT` (300 s) and reported as DNF beyond that. (In the recorded run the chr1 cap was
    applied by an external watchdog: the L=50 run was killed at 300 s and the L=100 run skipped.)
  - `bwa fastmap -l L idx query` (SMEMs, both strands, positions for up to 20 hits) and
    `ropebwt3 mem -t T -l L idx query` (SMEMs, both strands, no positions) from prebuilt indexes; the
    index builds (`bwa index`, `ropebwt3 build -t T -do`) are timed separately. SMEMs are a different match
    class, so these rows are context, not like-for-like.
  - Every tool writes its matches to a file (that I/O is in the timing). The sorted match sets stay in
    `~/sporos-bench/out/<pair>/` for re-checking (about 15 GB for chr1; the recorded run deleted them).
- **Correctness**: each output is reduced to the set of (query, strand, ref pos, query pos, length) and
  sorted (`LC_ALL=C sort -u`). Checks: sporos `-t 1` vs MUMmer for each mode; sporos `-t 10` vs `-t 1`;
  sporos `-maxmatch` vs copMEM2 (default and `-mf`, `-t 1` and `-t 10`); mumemto's forward-strand MUMs vs
  MUMmer `-mum` (info).
- **Seeds** (`run_seeds.sh`): kernel throughput, no FASTA parsing or output. sporos runs on 1-byte 2-bit
  codes (`CodeBuffer`), simd-minimizers on its packed 2-bit `PackedSeqVec`, `mm_sketch` on ASCII
  (and supports k ≤ 28, so its k=31 row is k=28). Single thread: best of 7 runs (sporos on chr1 and
  `mm_sketch`: 5), thread CPU time for sporos and wall time for the others; 10 threads: best wall time, with
  simd-minimizers fed 40 overlapping chunks by 10 scoped threads and `mm_sketch` 10 chunks (one per
  thread). The random inputs differ per tool (each harness draws its own uniform ACGT). Closed syncmers
  use simd-minimizers' parameterisation (s = k, window w). Output counts differ slightly between tools
  because hash functions differ.
- **CLI**: `sporos minimizers -k 21 -w 11 -t T chm13_chr1.fa > /dev/null`, end to end.

## Results (xaros, 2026-10-20; seeds and lookups re-run 2026-10-21)

The MEM and MUM rows below predate later sporos speedups (single-threaded chr1 MEMs at L=50 now take
19.5 s, MUMs 21.6 s); the competitors' rows stand. Current sporos numbers are in the main
[README](../../README.md#benchmarks).

**Machine**: Intel Core i9-7900X (10 cores, SMT off so 10 hardware threads, AVX-512), 91 GB RAM, Ubuntu, kernel 7.0.0-38,
gcc 15.2.0. Shared desktop box: it idles at about 0.5-1 busy CPU, and at 11:46 and 12:25 a separate
`mojo build` (~1.7 CPUs) ran alongside a single-threaded measurement (chr1 L=50 `sporos -maxmatch -t 1`,
chr1 L=50 `bwa fastmap`). 1-minute load average when the scripts started: 0.96 (`run_seeds.sh`, 11:10) and
1.73 (`run_mems.sh`, 11:13, right after the seed runs); per-row values are in the TSVs. Total runtime:
`run_seeds.sh` 3.5 min, `run_mems.sh` 1 h 42 min.

**Versions** ([`results/versions.tsv`](results/versions.tsv)):

| tool | version |
|---|---|
| sporos | `4db0ae6` plus uncommitted `cli/sporos_cli.mojo` edits (`-dirty`), Mojo 1.2.0.dev2026102005 |
| MUMmer | 4.0.1 (bioconda `mummer4`) |
| copMEM2 | git `8040165`, g++ -O3 -march=native -funroll-loops |
| mumemto | 1.4.2 (bioconda) |
| bwa | 0.7.19-r1273 (bioconda) |
| ropebwt3 | 3.10-r281 (bioconda) |
| simd-minimizers | git `a173cbe`, rustc 1.93.1, `-C target-cpu=native` |
| minimap2 (`mm_sketch`) | git `3c28777`, gcc -O3 -march=native |

### Summary

- **Correctness**: sporos equals MUMmer as a set in all 18 mode/L/pair combinations (up to 65.6 M MEMs on
  chr1 L=50), `-t 10` equals `-t 1` everywhere, and sporos `-maxmatch` equals copMEM2 (default and `-mf`,
  1 and 10 threads) in all 24 comparisons: 60/60 PASS. mumemto's forward-strand MUMs are a near-subset of
  MUMmer `-mum` (it also requires uniqueness against the reverse strands): 21887 of 21942 on ecoli L=50.
- **MEMs vs MUMmer**: 12.5-25x faster single-threaded on ecoli/chr21, 1.6x (`-mum` L=50) to 14.7x
  (`-maxmatch` L=100) on chr1; sporos' 10 threads add another 1.6-4.8x (MUMmer has no working threads).
- **MEMs vs copMEM2** (the strongest competitor): single-threaded 1.8-3.2x on chr21, 6.8-20x on ecoli
  (copMEM2 spends 0.5-1 s of system time on its fixed 1-2 GB tables), but only **1.3x on chr1**. copMEM2's
  `-t 10` gives no speedup on any pair (presumably because each query is a single record), so at 10 threads
  sporos leads 4.1-6.2x on chr1.
- **Minimizers/syncmers** (re-run 2026-10-21): one thread, sporos (16 lanes) is 1.7-2.3x faster than
  simd-minimizers on random and chr1 input, and 8 lanes (AVX2 width) 1.3-1.7x faster on chr1; with 10
  threads sporos leads 1.04-1.47x. `mm_sketch` is 15-18x slower than sporos.

### Where competitors win or sporos is weak

- **Multi-threaded seeds: the closest race.** On chr1 with 10 threads sporos leads by only 4% on forward
  (w=11) minimizers (0.182 vs 0.188 ns/bp), 1.2-1.5x elsewhere. sporos scales 5.8-7.5x from 1 to 10
  threads, simd-minimizers 9.5-10.4x. (Its harness keeps each thread's positions in a private vector and
  only counts them; sporos returns one merged, ordered position list.) Until 2026-10-21 simd-minimizers
  won here; the ordered-output chain changed that.
- **`-mum` / `-mumreference` on repeat-dense input.** On chr1 L=50 (65.6 M MEMs, vs 7.5 M at L=100)
  sporos `-mum` takes 56 s single-threaded vs MUMmer's 87 s (only 1.6x) and 32 s with 10
  threads (1.7x scaling); `-mumreference` 13 s at 10 threads vs 7.6 s for `-maxmatch`, which also writes
  280x more lines. Filtering all MEMs down to MUMs, not MEM finding, looks like the cost there.
- **Memory on chr1**: sporos 5.6 GB at L=50 (MUMmer 6.2 GB, copMEM2 2.8 GB, `-mf` 1.8 GB); 2.3 GB at
  L=100 (copMEM2 `-mf` 1.7 GB). On ecoli/chr21 sporos uses less than MUMmer at L=100 and far less than
  copMEM2.
- **Threads on small inputs**: sporos `-t 10` gains only ~2x on chr21 (0.47 vs 1.08 s at L=50); input
  reading and other fixed costs dominate at that size.
- Other notes: MUMmer 4.0.1's `-threads` needs `-k > 1`, which aborts; mumemto takes 13 s on E. coli,
  165-169 s on chr21 and did not finish chr1 within 300 s (DNF; killed); `ropebwt3 mem -t 10` and copMEM2
  `-t 10` do not speed up a single query record; `bwa fastmap` (SMEMs) is ~1.5x slower than sporos
  `-maxmatch -t 1` on chr1 L=50 but needs a 215 s index build first.

The `sporos -t 10` vs `-t 1` ratios and all ratios in the tables compare at equal thread counts.

### Seed-index lookups (2026-10-21)

hg38 chr21 index, 1 M simulated 150 bp reads, one thread, reads in memory, best of 3; harnesses in
[`harness/lookup`](harness/lookup/). sporos and minimap2 index the same (w, k)-minimizers (different
hashes, near-identical density) and skip seeds with more than 20 hits; strobealign indexes randstrobes
with its defaults for 150 bp reads and its own frequency cutoff.

| | Mbp/s | hits per read | index build (1 / 8 threads) |
|---|---|---|---|
| sporos `lookup_batch`, k=19 w=12 | **136** | 27.7 | 0.34 / 0.07 s |
| sporos `lookup` (one read per call), k=19 w=12 | 47 | 27.7 | |
| minimap2 `mm_sketch` + `mm_idx_get`, k=19 w=12 | 34.5 | 27.6 | 1.63 / 1.22 s |
| sporos `lookup_batch`, k=21 w=11 | **130** | 29.0 | 0.39 / 0.07 s |
| minimap2, k=21 w=11 | 32.1 | 29.2 | 1.68 / 1.25 s |
| strobealign `find_hits` (one hit per seed) | 14.5 | 51.0 | 1.78 / 1.34 s |
| strobealign `find_hits` + every position | 13.1 | 97.2 | |

- sporos `lookup_batch` on 8 threads: 570 Mbp/s.
- strobealign does more per seed (both orientations, partial-strobe lookups, rescue), so its row is a
  reference point rather than the same work.
- minimap2 and strobealign build times include reading the FASTA; sporos's do not (about 0.03 s here).
- The harness is a lighter path than minimap2's own seeding (`mm_collect_matches` also sorts and filters).
- One-read `lookup` calls were 37 Mbp/s before short reads skipped the SIMD lanes (whose warm-up dominates
  on 150 bp): about the same as minimap2's per-read loop.

### MEMs and MUMs

**ecoli, L=50** (best wall s; user s; max RSS MB)

| tool | T | -maxmatch | -mum | -mumreference | other | max RSS MB |
|---|---:|---:|---:|---:|---|---:|
| sporos | 1 | 0.13 (0.08u) | 0.13 (0.08u) | 0.13 (0.07u) |  | 88 |
| sporos | 10 | 0.06 (0.08u) | 0.07 (0.09u) | 0.07 (0.08u) |  | 88 |
| mummer | 1 | 1.97 (1.93u) = 15.2x sporos | 1.70 (1.66u) = 13.1x sporos | 1.69 (1.66u) = 13.0x sporos |  | 72 |
| copmem2 | 1 | 1.67 (0.74u) = 12.8x sporos |  |  |  | 2095 |
| copmem2 -mf | 1 | 0.88 (0.43u) = 6.8x sporos |  |  |  | 1071 |
| copmem2 | 10 | 1.67 (0.74u) = 27.8x sporos |  |  |  | 2095 |
| copmem2 -mf | 10 | 0.89 (0.43u) = 14.8x sporos |  |  |  | 1071 |
| mumemto | 1 |  |  |  | mum (both strands): 13.19 (12.94u) | 348 |
| bwa fastmap | 1 |  |  |  | SMEM (fwd+rc): 1.29 (1.27u) | 15 |
| ropebwt3 mem | 1 |  |  |  | SMEM (fwd+rc, no pos): 2.12 (2.11u) | 22 |
| ropebwt3 mem | 10 |  |  |  | SMEM (fwd+rc, no pos): 2.12 (2.10u) | 22 |

**ecoli, L=100** (best wall s; user s; max RSS MB)

| tool | T | -maxmatch | -mum | -mumreference | other | max RSS MB |
|---|---:|---:|---:|---:|---|---:|
| sporos | 1 | 0.08 (0.05u) | 0.08 (0.05u) | 0.08 (0.04u) |  | 61 |
| sporos | 10 | 0.05 (0.05u) | 0.05 (0.05u) | 0.05 (0.05u) |  | 61 |
| mummer | 1 | 1.90 (1.86u) = 23.8x sporos | 1.70 (1.66u) = 21.2x sporos | 1.67 (1.63u) = 20.9x sporos |  | 72 |
| copmem2 | 1 | 1.58 (0.66u) = 19.8x sporos |  |  |  | 2092 |
| copmem2 -mf | 1 | 0.81 (0.34u) = 10.1x sporos |  |  |  | 1068 |
| copmem2 | 10 | 1.58 (0.66u) = 31.6x sporos |  |  |  | 2092 |
| copmem2 -mf | 10 | 0.82 (0.35u) = 16.4x sporos |  |  |  | 1068 |
| mumemto | 1 |  |  |  | mum (both strands): 13.02 (12.77u) | 348 |
| bwa fastmap | 1 |  |  |  | SMEM (fwd+rc): 1.27 (1.24u) | 15 |
| ropebwt3 mem | 1 |  |  |  | SMEM (fwd+rc, no pos): 1.65 (1.64u) | 21 |
| ropebwt3 mem | 10 |  |  |  | SMEM (fwd+rc, no pos): 1.64 (1.63u) | 21 |

Index builds: bwa index -t1 1.76 s (49 MB); ropebwt3 build -t1 0.78 s (51 MB); ropebwt3 build -t10 0.56 s (51 MB)

**chr21, L=50** (best wall s; user s; max RSS MB)

| tool | T | -maxmatch | -mum | -mumreference | other | max RSS MB |
|---|---:|---:|---:|---:|---|---:|
| sporos | 1 | 1.08 (0.73u) | 1.17 (0.80u) | 1.10 (0.73u) |  | 522 |
| sporos | 10 | 0.47 (0.83u) | 0.57 (0.90u) | 0.50 (0.83u) |  | 522 |
| mummer | 1 | 17.42 (17.03u) = 16.1x sporos | 14.60 (14.19u) = 12.5x sporos | 14.56 (14.15u) = 13.2x sporos |  | 610 |
| copmem2 | 1 | 2.75 (1.74u) = 2.5x sporos |  |  |  | 2202 |
| copmem2 -mf | 1 | 1.96 (1.42u) = 1.8x sporos |  |  |  | 1169 |
| copmem2 | 10 | 2.75 (1.77u) = 5.9x sporos |  |  |  | 2203 |
| copmem2 -mf | 10 | 1.95 (1.41u) = 4.1x sporos |  |  |  | 1169 |
| mumemto | 1 |  |  |  | mum (both strands): 169.45 (166.86u) | 3130 |
| bwa fastmap | 1 |  |  |  | SMEM (fwd+rc): 11.55 (11.45u) | 101 |
| ropebwt3 mem | 1 |  |  |  | SMEM (fwd+rc, no pos): 25.04 (24.96u) | 130 |
| ropebwt3 mem | 10 |  |  |  | SMEM (fwd+rc, no pos): 25.07 (24.97u) | 130 |

**chr21, L=100** (best wall s; user s; max RSS MB)

| tool | T | -maxmatch | -mum | -mumreference | other | max RSS MB |
|---|---:|---:|---:|---:|---|---:|
| sporos | 1 | 0.66 (0.40u) | 0.68 (0.40u) | 0.67 (0.40u) |  | 328 |
| sporos | 10 | 0.40 (0.45u) | 0.42 (0.49u) | 0.41 (0.47u) |  | 329 |
| mummer | 1 | 16.52 (16.14u) = 25.0x sporos | 14.55 (14.18u) = 21.4x sporos | 14.56 (14.16u) = 21.7x sporos |  | 610 |
| copmem2 | 1 | 2.08 (1.13u) = 3.2x sporos |  |  |  | 2159 |
| copmem2 -mf | 1 | 1.27 (0.78u) = 1.9x sporos |  |  |  | 1132 |
| copmem2 | 10 | 2.06 (1.11u) = 5.1x sporos |  |  |  | 2158 |
| copmem2 -mf | 10 | 1.26 (0.75u) = 3.1x sporos |  |  |  | 1132 |
| mumemto | 1 |  |  |  | mum (both strands): 165.38 (163.14u) | 3130 |
| bwa fastmap | 1 |  |  |  | SMEM (fwd+rc): 11.00 (10.91u) | 101 |
| ropebwt3 mem | 1 |  |  |  | SMEM (fwd+rc, no pos): 17.87 (17.80u) | 118 |
| ropebwt3 mem | 10 |  |  |  | SMEM (fwd+rc, no pos): 17.87 (17.80u) | 118 |

Index builds: bwa index -t1 33.35 s (74 MB); ropebwt3 build -t1 8.29 s (493 MB); ropebwt3 build -t10 6.59 s (493 MB)

**chr1, L=50** (best wall s; user s; max RSS MB)

| tool | T | -maxmatch | -mum | -mumreference | other | max RSS MB |
|---|---:|---:|---:|---:|---|---:|
| sporos | 1 | 36.77 (32.25u) | 56.10 (52.87u) | 37.69 (34.41u) |  | 5557 |
| sporos | 10 | 7.64 (34.92u) | 32.17 (56.18u) | 13.02 (37.14u) |  | 5218 |
| mummer | 1 | 159.21 (153.81u) = 4.3x sporos | 87.49 (85.61u) = 1.6x sporos | 87.86 (86.05u) = 2.3x sporos |  | 6246 |
| copmem2 | 1 | 47.07 (45.06u) = 1.3x sporos |  |  |  | 2794 |
| copmem2 -mf | 1 | 46.91 (45.36u) = 1.3x sporos |  |  |  | 1758 |
| copmem2 | 10 | 47.12 (45.08u) = 6.2x sporos |  |  |  | 2794 |
| copmem2 -mf | 10 | 46.51 (44.99u) = 6.1x sporos |  |  |  | 1758 |
| mumemto | 1 |  |  |  | mum (both strands): DNF (killed at 300 s) |  |
| bwa fastmap | 1 |  |  |  | SMEM (fwd+rc): 53.47 (53.02u) | 596 |
| ropebwt3 mem | 1 |  |  |  | SMEM (fwd+rc, no pos): 184.05 (183.57u) | 771 |
| ropebwt3 mem | 10 |  |  |  | SMEM (fwd+rc, no pos): 188.89 (188.39u) | 771 |

**chr1, L=100** (best wall s; user s; max RSS MB)

| tool | T | -maxmatch | -mum | -mumreference | other | max RSS MB |
|---|---:|---:|---:|---:|---|---:|
| sporos | 1 | 8.95 (7.25u) | 10.48 (8.91u) | 8.96 (7.39u) |  | 2291 |
| sporos | 10 | 2.72 (9.37u) | 4.75 (10.87u) | 3.25 (9.74u) |  | 2291 |
| mummer | 1 | 131.76 (129.41u) = 14.7x sporos | 93.84 (91.92u) = 9.0x sporos | 87.56 (85.64u) = 9.8x sporos |  | 3123 |
| copmem2 | 1 | 11.95 (10.50u) = 1.3x sporos |  |  |  | 2723 |
| copmem2 -mf | 1 | 11.43 (10.46u) = 1.3x sporos |  |  |  | 1674 |
| copmem2 | 10 | 11.30 (9.83u) = 4.2x sporos |  |  |  | 2723 |
| copmem2 -mf | 10 | 11.05 (10.02u) = 4.1x sporos |  |  |  | 1674 |
| mumemto | 1 |  |  |  | mum (both strands): DNF (skipped: L=50 run passed 300 s) |  |
| bwa fastmap | 1 |  |  |  | SMEM (fwd+rc): 52.92 (52.51u) | 596 |
| ropebwt3 mem | 1 |  |  |  | SMEM (fwd+rc, no pos): 185.78 (185.32u) | 749 |
| ropebwt3 mem | 10 |  |  |  | SMEM (fwd+rc, no pos): 192.65 (192.23u) | 749 |

Index builds: bwa index -t1 214.53 s (359 MB); ropebwt3 build -t1 50.82 s (2614 MB); ropebwt3 build -t10 41.27 s (2614 MB)

### Correctness checks

| pair | L | check | n (a) | n (b) | only a | only b | result |
|---|---:|---|---:|---:|---:|---:|---|
| ecoli | 50 | sporos=mummer -maxmatch | 23008 | 23008 | 0 | 0 | PASS |
| ecoli | 50 | sporos -t10=-t1 -maxmatch | 23008 | 23008 | 0 | 0 | PASS |
| ecoli | 50 | sporos=mummer -mum | 21942 | 21942 | 0 | 0 | PASS |
| ecoli | 50 | sporos -t10=-t1 -mum | 21942 | 21942 | 0 | 0 | PASS |
| ecoli | 50 | sporos=mummer -mumreference | 22158 | 22158 | 0 | 0 | PASS |
| ecoli | 50 | sporos -t10=-t1 -mumreference | 22158 | 22158 | 0 | 0 | PASS |
| ecoli | 50 | sporos=copmem2 -t1 | 23008 | 23008 | 0 | 0 | PASS |
| ecoli | 50 | sporos=copmem2 -mf -t1 | 23008 | 23008 | 0 | 0 | PASS |
| ecoli | 50 | sporos=copmem2 -t10 | 23008 | 23008 | 0 | 0 | PASS |
| ecoli | 50 | sporos=copmem2 -mf -t10 | 23008 | 23008 | 0 | 0 | PASS |
| ecoli | 50 | mumemto(+,+) vs mummer -mum | 21887 | 21942 | 1 | 56 | info |
| ecoli | 100 | sporos=mummer -maxmatch | 12048 | 12048 | 0 | 0 | PASS |
| ecoli | 100 | sporos -t10=-t1 -maxmatch | 12048 | 12048 | 0 | 0 | PASS |
| ecoli | 100 | sporos=mummer -mum | 11497 | 11497 | 0 | 0 | PASS |
| ecoli | 100 | sporos -t10=-t1 -mum | 11497 | 11497 | 0 | 0 | PASS |
| ecoli | 100 | sporos=mummer -mumreference | 11609 | 11609 | 0 | 0 | PASS |
| ecoli | 100 | sporos -t10=-t1 -mumreference | 11609 | 11609 | 0 | 0 | PASS |
| ecoli | 100 | sporos=copmem2 -t1 | 12048 | 12048 | 0 | 0 | PASS |
| ecoli | 100 | sporos=copmem2 -mf -t1 | 12048 | 12048 | 0 | 0 | PASS |
| ecoli | 100 | sporos=copmem2 -t10 | 12048 | 12048 | 0 | 0 | PASS |
| ecoli | 100 | sporos=copmem2 -mf -t10 | 12048 | 12048 | 0 | 0 | PASS |
| ecoli | 100 | mumemto(+,+) vs mummer -mum | 11475 | 11497 | 1 | 23 | info |
| chr21 | 50 | sporos=mummer -maxmatch | 508649 | 508649 | 0 | 0 | PASS |
| chr21 | 50 | sporos -t10=-t1 -maxmatch | 508649 | 508649 | 0 | 0 | PASS |
| chr21 | 50 | sporos=mummer -mum | 220360 | 220360 | 0 | 0 | PASS |
| chr21 | 50 | sporos -t10=-t1 -mum | 220360 | 220360 | 0 | 0 | PASS |
| chr21 | 50 | sporos=mummer -mumreference | 230872 | 230872 | 0 | 0 | PASS |
| chr21 | 50 | sporos -t10=-t1 -mumreference | 230872 | 230872 | 0 | 0 | PASS |
| chr21 | 50 | sporos=copmem2 -t1 | 508649 | 508649 | 0 | 0 | PASS |
| chr21 | 50 | sporos=copmem2 -mf -t1 | 508649 | 508649 | 0 | 0 | PASS |
| chr21 | 50 | sporos=copmem2 -t10 | 508649 | 508649 | 0 | 0 | PASS |
| chr21 | 50 | sporos=copmem2 -mf -t10 | 508649 | 508649 | 0 | 0 | PASS |
| chr21 | 50 | mumemto(+,+) vs mummer -mum | 214881 | 220360 | 2 | 5481 | info |
| chr21 | 100 | sporos=mummer -maxmatch | 113188 | 113188 | 0 | 0 | PASS |
| chr21 | 100 | sporos -t10=-t1 -maxmatch | 113188 | 113188 | 0 | 0 | PASS |
| chr21 | 100 | sporos=mummer -mum | 103127 | 103127 | 0 | 0 | PASS |
| chr21 | 100 | sporos -t10=-t1 -mum | 103127 | 103127 | 0 | 0 | PASS |
| chr21 | 100 | sporos=mummer -mumreference | 103893 | 103893 | 0 | 0 | PASS |
| chr21 | 100 | sporos -t10=-t1 -mumreference | 103893 | 103893 | 0 | 0 | PASS |
| chr21 | 100 | sporos=copmem2 -t1 | 113188 | 113188 | 0 | 0 | PASS |
| chr21 | 100 | sporos=copmem2 -mf -t1 | 113188 | 113188 | 0 | 0 | PASS |
| chr21 | 100 | sporos=copmem2 -t10 | 113188 | 113188 | 0 | 0 | PASS |
| chr21 | 100 | sporos=copmem2 -mf -t10 | 113188 | 113188 | 0 | 0 | PASS |
| chr21 | 100 | mumemto(+,+) vs mummer -mum | 101765 | 103127 | 2 | 1364 | info |
| chr1 | 50 | sporos=mummer -maxmatch | 65640548 | 65640548 | 0 | 0 | PASS |
| chr1 | 50 | sporos -t10=-t1 -maxmatch | 65640548 | 65640548 | 0 | 0 | PASS |
| chr1 | 50 | sporos=mummer -mum | 232880 | 232880 | 0 | 0 | PASS |
| chr1 | 50 | sporos -t10=-t1 -mum | 232880 | 232880 | 0 | 0 | PASS |
| chr1 | 50 | sporos=mummer -mumreference | 431409 | 431409 | 0 | 0 | PASS |
| chr1 | 50 | sporos -t10=-t1 -mumreference | 431409 | 431409 | 0 | 0 | PASS |
| chr1 | 50 | sporos=copmem2 -t1 | 65640548 | 65640548 | 0 | 0 | PASS |
| chr1 | 50 | sporos=copmem2 -mf -t1 | 65640548 | 65640548 | 0 | 0 | PASS |
| chr1 | 50 | sporos=copmem2 -t10 | 65640548 | 65640548 | 0 | 0 | PASS |
| chr1 | 50 | sporos=copmem2 -mf -t10 | 65640548 | 65640548 | 0 | 0 | PASS |
| chr1 | 50 | mumemto(+,+) vs mummer -mum | NA | 232880 | NA | 232880 | DNF |
| chr1 | 100 | sporos=mummer -maxmatch | 7469747 | 7469747 | 0 | 0 | PASS |
| chr1 | 100 | sporos -t10=-t1 -maxmatch | 7469747 | 7469747 | 0 | 0 | PASS |
| chr1 | 100 | sporos=mummer -mum | 202839 | 202839 | 0 | 0 | PASS |
| chr1 | 100 | sporos -t10=-t1 -mum | 202839 | 202839 | 0 | 0 | PASS |
| chr1 | 100 | sporos=mummer -mumreference | 334935 | 334935 | 0 | 0 | PASS |
| chr1 | 100 | sporos -t10=-t1 -mumreference | 334935 | 334935 | 0 | 0 | PASS |
| chr1 | 100 | sporos=copmem2 -t1 | 7469747 | 7469747 | 0 | 0 | PASS |
| chr1 | 100 | sporos=copmem2 -mf -t1 | 7469747 | 7469747 | 0 | 0 | PASS |
| chr1 | 100 | sporos=copmem2 -t10 | 7469747 | 7469747 | 0 | 0 | PASS |
| chr1 | 100 | sporos=copmem2 -mf -t10 | 7469747 | 7469747 | 0 | 0 | PASS |
| chr1 | 100 | mumemto(+,+) vs mummer -mum | NA | 202839 | NA | 202839 | DNF |

### Minimizers and syncmers

**random100M, 1 thread** (ns/bp, lower is better)

| mode | k | w | sporos (lanes=8) | sporos (lanes=16) | simd-minimizers | minimap2 mm_sketch |
|---|---:|---:|---:|---:|---:|---:|
| fwd | 21 | 11 | 1.3480 | 0.9505 | 1.7910 |  |
| can | 21 | 11 | 1.4355 | 1.0248 | 2.2000 | 18.2650 |
| closed-fwd | 21 | 11 | 1.5353 | 1.1569 | 1.9990 |  |
| closed-can | 21 | 11 | 1.5047 | 1.1582 | 2.5770 |  |
| fwd | 31 | 19 | 1.2507 | 0.8448 | 1.6510 |  |
| can | 31 | 19 | 1.3285 | 0.9225 | 1.9670 |  |
| closed-fwd | 31 | 19 | 1.4186 | 1.0040 | 1.8740 |  |
| closed-can | 31 | 19 | 1.3914 | 1.0130 | 2.3370 |  |
| open-fwd | 16 | 5 | 1.3829 | 1.0656 | 2.0360 |  |
| open-can | 15 | 5 |  |  | 2.5320 |  |
| can | 28 | 19 |  |  |  | 15.9120 |

**chm13_chr1, 1 thread** (ns/bp, lower is better)

| mode | k | w | sporos (lanes=8) | sporos (lanes=16) | simd-minimizers | minimap2 mm_sketch |
|---|---:|---:|---:|---:|---:|---:|
| fwd | 21 | 11 | 1.3286 | 1.0513 | 1.7900 |  |
| can | 21 | 11 | 1.4303 | 1.0463 | 2.1260 | 17.0560 |
| closed-fwd | 21 | 11 | 1.3928 | 1.1521 | 1.9870 |  |
| closed-can | 21 | 11 | 1.5185 | 1.1644 | 2.5490 |  |
| fwd | 31 | 19 | 1.2240 | 0.9353 | 1.6250 |  |
| can | 31 | 19 | 1.3218 | 0.9328 | 1.9380 |  |
| closed-fwd | 31 | 19 | 1.2671 | 0.9984 | 1.8380 |  |
| closed-can | 31 | 19 | 1.4134 | 1.0243 | 2.3520 |  |
| open-fwd | 16 | 5 |  |  | 1.9750 |  |
| open-can | 15 | 5 |  |  | 2.5340 |  |
| can | 28 | 19 |  |  |  | 15.2630 |

**chm13_chr1, 10 threads** (ns/bp, lower is better)

| mode | k | w | sporos (lanes=16) | simd-minimizers | minimap2 mm_sketch |
|---|---:|---:|---:|---:|---:|
| fwd | 21 | 11 | 0.1816 | 0.1880 |  |
| can | 21 | 11 | 0.1730 | 0.2130 | 1.9850 |
| closed-fwd | 21 | 11 | 0.1745 | 0.2060 |  |
| closed-can | 21 | 11 | 0.1773 | 0.2550 |  |
| fwd | 31 | 19 | 0.1254 | 0.1590 |  |
| can | 31 | 19 | 0.1383 | 0.1920 |  |
| closed-fwd | 31 | 19 | 0.1350 | 0.1780 |  |
| closed-can | 31 | 19 | 0.1547 | 0.2270 |  |
| open-fwd | 16 | 5 |  | 0.2000 |  |
| open-can | 15 | 5 |  | 0.2480 |  |
| can | 28 | 19 |  |  | 1.9350 |

**sporos minimizers CLI, CHM13 chr1** (`-k 21 -w 11`, output to /dev/null)

| threads | wall s | user s | max RSS MB |
|---:|---:|---:|---:|
| 1 | 0.90 | 0.64 | 569 |
| 10 | 0.27 | 1.00 | 488 |
