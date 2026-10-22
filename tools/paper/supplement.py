"""The supplement, written once; render.py turns it into paper/supplement.tex and build/paper/supplement.html.
Same markup as main_text.py; its tables come from tables.STABLES."""
from tables import STABLES, V

# The title is the paper's (render.py: SUPPLEMENT_TITLE).
TABLES = STABLES

BODY = [
('sec', 'Setup'),
('p',
 '`benchmarks/competitors/setup_xaros.sh` in the sporos repository installs every tool at the version of Table S5, '
 'downloads the genomes from NCBI and UCSC, simulates the reads (`sim_reads.py`, fixed seeds) and builds the '
 'harnesses; `run_seeds.sh`, `run_lookups.sh`, `run_mems.sh`, `run_capped.sh` and `run_smems.sh` run the '
 'benchmarks and write the raw measurements that the main tables are made from. Each timed command ran under '
 '`/usr/bin/time -v` after the machine had kept at most 1.5 processors busy for half a second; wall time is the best '
 'of three runs (one for runs past 60 s), memory the largest peak resident set size of those runs. Kernel '
 'benchmarks (Table 1, Table 2) time only the computation, inside each tool\'s own harness.'),
('table', 'versions'),
('p', 'The commands, with `R` and `Q` the reference and query FASTA files, `L` the minimum length and `T` the threads:'),
('bullets', [
 'sporos: `sporos mem -maxmatch|-mum|-mumreference -l L -n -t T R Q`; capped: `sporos mem -maxmatch|-mum -b '
 '--max-occ N -l 50 -t T R Q`; SMEMs: `sporos index -t T R`, then `sporos smem -l 19 -s f -t T R reads.fa`.',
 'MUMmer4: `mummer -maxmatch|-mum|-mumreference -n -l L R Q`. Its `-qthreads` divides only between query records, '
 'and `-threads` requires `-k` above 1, which aborted, so it ran on one thread.',
 'copMEM2: `copmem2 -q [-mf] -t T -l L -o out R Q`; `-t 10` did not change its time on a single query record.',
 'Mumemto: `mumemto -k 2 -f N -l 50` (per-genome cap) and `mumemto -k 2 -f 0 -F N -l 50` (total cap); MUMs: '
 '`mumemto -l L`.',
 'minibwa: `minibwa index -t T R prefix`, then `minibwa fastmap -l 19 -s f -b 64 prefix reads.fa` (batches of 64, '
 'its fastest setting).',
 'BWA: `bwa index`, then `bwa fastmap -l 19 index reads.fa`.',
 'ropebwt3: `ropebwt3 build -t T -do index.fmd R`, then `ropebwt3 mem -t T -l 19 -c f index.fmd reads.fa`.',
 'Minimizers and syncmers: `benchmarks/bench_seeds.mojo` and `bench_seeds_fasta.mojo` (sporos), the '
 '`harness/smbench` program on simd-minimizers\' packed sequences, and `harness/mm2` on minimap2\'s `mm_sketch`.',
 'Seed lookups: `benchmarks/bench_index.mojo` (sporos), `harness/lookup/mm2_lookup.c` (minimap2\'s `mm_sketch` and '
 '`mm_idx_get`) and `harness/lookup/sa_bench` (strobealign\'s `randstrobes_query` and `find_hits`).']),

('sec', 'Additional results'),
('sub', 'Minimizers and syncmers of random DNA'),
('table', 'seeds_random'),
('sub', 'MEMs and MUMs'),
('p',
 'Tables S2 to S4 complete main Table 3: MEMs and MUMs of at least 100 bp, MEMs unique in the reference, and peak '
 'memory. copMEM2 holds fixed tables of 1 to 2 GB whatever the input size, which dominates its time on the small '
 'pair; sporos\'s anchor index grows with the reference and the sampling, and its peak on chromosome 1 at L = 50 '
 'comes from the 65 million MEMs it holds before filtering.'),
('table', 'mems100'),
('table', 'mumreference'),
('table', 'memory'),
('sub', 'FM-index builds'),
('table', 'builds'),

('sec', 'SMEM outputs'),
('p',
 'Table S7 compares the SMEMs of sporos with those of `bwa fastmap` and minibwa, by read, interval, number of '
 f'occurrences and listed positions. Of the {V["smem_only_other"]} SMEMs that only the other tool reports, '
 f'{V["smem_through_n"]} list a position whose reference span holds an N, which BWA and minibwa replace by random '
 f'bases when they build their index and sporos never matches, and {V["smem_unlisted"]} list no positions, so '
 f'that this check cannot tell; {V["smem_unexplained"]} remain. '
 +('Ten threads gave the same output as one for every read set and f.' if V["smem_threads_same"] else '')),
('table', 'smem_checks'),
]
