"""The supplement, written once; render.py turns it into paper/supplement.tex and build/paper/supplement.html.
Same markup as main_text.py; its tables come from tables.STABLES."""
from tables import STABLES, V

# The title is the paper's (render.py: SUPPLEMENT_TITLE).
TABLES = STABLES

BODY = [
('sec', 'Setup'),
('p',
 'In the sporos repository, `benchmarks/competitors/setup_xaros.sh` installs every tool at the version listed in '
 'Table S1, downloads the genomes from NCBI and UCSC, simulates the reads (`sim_reads.py`, fixed seeds) and builds '
 'the benchmark harnesses; `run_all.sh` runs every benchmark and writes the raw measurements from which all tables '
 'are generated. Each timed command ran under `/usr/bin/time -v` once at most 1.5 processors had been busy over the '
 'preceding half second. Wall time is the best of three runs (one for runs exceeding 60 s), and memory the largest '
 'peak resident set size. The kernel benchmarks of main Tables 1 and 2 time only the computation, within each '
 'tool\'s harness; the minimizer and syncmer harnesses perform one untimed run first, so that threads, memory pages '
 'and clock frequency are warm for every tool.'),
('table', 'versions'),
('p', 'Commands, with `R` and `Q` the reference and query FASTA files, `L` the minimum length and `T` the threads:'),
('bullets', [
 'sporos: `sporos mem -maxmatch|-mum|-mumreference -l L -n -t T R Q`; capped: `sporos mem -maxmatch|-mum -b '
 '--max-occ N -l 50 -t T R Q`; SMEMs: `sporos index -t T R`, then `sporos smem -l 19 -s f -t T R reads.fa`.',
 'MUMmer4: `mummer -maxmatch|-mum|-mumreference -n -l L R Q`. `-qthreads` divides work only between query records, '
 'and `-threads` requires `-k` > 1, a combination that aborted in our runs; MUMmer4 therefore ran on one thread.',
 'copMEM2: `copmem2 -q [-mf] -t T -l L -o out R Q`; `-t 10` did not reduce its time on a single query record.',
 'Mumemto: `mumemto -k 2 -f N -l 50` (per-genome cap) and `mumemto -k 2 -f 0 -F N -l 50` (total cap); MUMs: '
 '`mumemto -l L`.',
 'minibwa: `minibwa index -t T R prefix`, then `minibwa fastmap -l 19 -s f -b 64 prefix reads.fa`.',
 'BWA: `bwa index`, then `bwa fastmap -l 19 index reads.fa`.',
 'ropebwt3: `ropebwt3 build -t T -do index.fmd R`, then `ropebwt3 mem -t T -l 19 -c f index.fmd reads.fa`.',
 'Minimizers and syncmers: `benchmarks/bench_seeds.mojo` and `bench_seeds_fasta.mojo` (sporos), '
 '`harness/smbench` on packed sequences (simd-minimizers) and `harness/mm2` (minimap2\'s `mm_sketch`).',
 'Seed lookups: `benchmarks/bench_index.mojo` (sporos), `harness/lookup/mm2_lookup.c` (minimap2\'s `mm_sketch` and '
 '`mm_idx_get`) and `harness/lookup/sa_bench` (strobealign\'s `randstrobes_query` and `find_hits`).']),

('sec', 'Additional results'),
('sub', 'Minimizers and syncmers of random DNA'),
('table', 'seeds_random'),
('sub', 'MEMs and MUMs'),
('p',
 'Tables S3 to S5 complement main Table 3 with MEMs and MUMs of at least 100 bp, reference-unique MEMs and peak '
 f'memory. copMEM2 used at least {V["copmem2_rss_min"]} GB on every input. The memory of sporos grows with the '
 f'reference and with the number of MEMs it holds, {V["checks_max_mems"]} million on chromosome 1 at L = 50.'),
('table', 'mems100'),
('table', 'mumreference'),
('table', 'memory'),
('sub', 'FM-index construction'),
('table', 'builds'),

('sec', 'SMEM outputs'),
('p',
 'Table S7 compares the SMEMs of sporos with those of `bwa fastmap` and minibwa by read, interval, number of '
 f'occurrences and listed positions. Of the {V["smem_only_other"]} SMEMs reported only by the other tool, '
 f'{V["smem_through_n"]} list a position whose reference span contains an N, which BWA and minibwa replace with '
 f'random bases during indexing and sporos never matches; {V["smem_unlisted"]} list no positions, as they have more '
 f'than 20 occurrences, and cannot be checked this way; {V["smem_unexplained"]} remain. '
 + ('Ten threads gave the same output as one for every read set and f.' if V["smem_threads_same"] else '')),
('table', 'smem_checks'),
]
