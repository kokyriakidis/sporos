"""The main paper, written once; render.py turns it into paper/ (LaTeX, an Overleaf project) and build/paper/.

Inline markup: **bold**, *italic*, `code`, [@Key; @Key] parenthetical citation, [@!Key] narrative citation,
\\( \\) inline maths, \\[ \\] display maths, [#tab:x] / [#fig:x] cross-references.
Blocks: ("sec", title), ("sub", title), ("p", text), ("list", [items]) numbered, ("bullets", [items]), ("table", key),
("figure", key). Every number comes from tables.V or the tables, which read the benchmark results.
"""
from tables import TABLES, V, span

TITLE = "sporos: one fast, hardware-aware library to seed them all"
# Authors in order, each with an email; both are corresponding authors.
AUTHORS = [("Konstantinos Kyriakidis", "kkyriaki@ucsc.edu"), ("Benedict Paten", "bpaten@ucsc.edu")]
AFFILIATION = "UC Santa Cruz Genomics Institute, Santa Cruz, CA, USA"
KEYWORDS = ["seeding", "minimizers", "syncmers", "maximal exact matches", "MUMs", "SMEMs", "FM-index", "SIMD"]
WIDE_TABLES = ("seeds", "mems", "capped", "smems")

ABSTRACT = [
 ("Motivation",
  "Nearly every comparison of DNA sequences, from mapping reads to aligning genomes and building pangenomes, starts "
  "by finding seeds: exact matches between a query and a reference. Seeding is a substantial part of the work, "
  "40% to 60% of single-threaded BWA-MEM's time on short reads, and the data it runs on keep growing, from complete "
  "telomere-to-telomere "
  "genomes to pangenomes of many individuals. Yet seeds are computed by a patchwork of tools, one for each kind of "
  "seed, each with its own index and definitions; several use a single processor core per query sequence, and the "
  "newest kinds of seed, frequency-capped maximal exact matches and f-SMEMs, are offered by one or two tools each."),
 ("Results",
  "We present sporos, a library that computes all of these seeds exactly, with one set of definitions and one "
  "interface. Its main algorithmic contribution counts the copies of a MEM's string from the MEMs themselves: a "
  "string occurs in the reference once for every MEM whose query interval contains its own, so frequency-capped "
  "MEMs need no index of the genomes beyond the one that finds the MEMs. On one Intel Skylake-X processor, sporos "
  f"computed minimizers and syncmers {V['seed_t1_min']}–{V['seed_t1_max']} times faster than simd-minimizers on "
  f"one core, looked seeds up {V['lookup_vs_mm2']} times faster than minimap2, found MEMs between two human "
  f"chromosomes 1 {V['copmem2_chr1_min']}–{V['copmem2_chr1_max']} times faster than copMEM2 and MUMs up to "
  f"{V['mummer_max']} times faster than MUMmer4 with identical output, frequency-capped MEMs "
  f"{V['capped_min']}–{V['capped_max']} times faster than Mumemto, and SMEMs of short reads "
  f"{V['smem_150bp_mb']} times faster than minibwa."),
 ("Availability and implementation",
  "sporos is free software under the MIT license, written in Mojo with C, Python and command-line interfaces, "
  "at https://github.com/kokyriakidis/sporos."),
]

BODY = [
('sec', 'Introduction'),
('p',
 'Almost every comparison of DNA sequences starts from seeds: exact matches between a query and a reference that '
 'locate where the two may align. A read mapper chains the seeds of a read before aligning it [@Li2018; @Li2013], '
 'a whole-genome aligner anchors its alignment on matches unique to both genomes [@Delcher1999; @Marcais2018], '
 'and pangenome tools describe the shared structure of many genomes by the matches they have in common '
 '[@Shivakumar2025]. Two families of seeds serve these tasks. *Sampled k-mers* keep a subset of the k-mers of '
 'each sequence, chosen by their content so that two copies of a region keep the same ones: minimizers, the '
 'smallest k-mer of every window of w consecutive k-mers under a hash order [@Schleimer2003; @Roberts2004], and '
 'syncmers, k-mers whose smallest s-mer lies at a chosen offset [@Edgar2021]. *Maximal exact matches* (MEMs) '
 'extend a match in both directions until the sequences differ; their unique subset, the maximal unique matches '
 '(MUMs), anchors whole-genome alignment [@Delcher1999], and the super-maximal exact matches (SMEMs) of a read, '
 'the MEMs no longer one contains, seed BWA-MEM [@Li2012; @Li2013]. Seeding is not a negligible step: on short '
 'reads, single-threaded BWA-MEM spends 40% to 60% of its time finding SMEMs and locating their occurrences '
 '[@Vasimuddin2019], and accelerating that step has been the subject of dedicated hardware [@Subramaniyan2021]. '
 'The data that seeding runs on keep growing: complete telomere-to-telomere assemblies [@Nurk2022] and pangenomes '
 'of dozens of individuals [@Liao2023] make whole-genome comparisons and repeat-aware matching routine. The '
 'processors that run it have changed too, with vector units of up to 512 bits and tens of cores, and a seeding '
 'method keeps pace only if it uses both.'),
('p',
 'Each kind of seed is computed by its own tool, on its own data structure. Minimizers are computed inside read '
 'mappers such as minimap2 [@Li2018], with syncmers and strobemers in strobealign [@Sahlin2021; @Sahlin2022]; '
 'simd-minimizers computes them with SIMD instructions several times faster than earlier code '
 '[@GrootKoerkamp2025simd]. MEMs and MUMs between genomes come from suffix trees and suffix arrays in MUMmer '
 '[@Kurtz2004; @Marcais2018], from sparse suffix arrays [@Khan2009; @Vyverman2013] and from hash tables of '
 'sampled k-mers in E-MEM and copMEM [@Khiste2015; @Grabowski2019; @Grabowski2023]. SMEMs come from FM-indexes '
 '[@Ferragina2000] in BWA [@Li2013], bwa-mem2 [@Vasimuddin2019], ropebwt3 [@Li2024] and minibwa [@Li2026]. '
 'Recent variants each live in a single tool. Mumemto finds multi-MEMs across a collection of genomes and caps how '
 'often a match may occur, relaxing uniqueness to rare repeats, from a suffix array, LCP array and document array '
 'built over the whole collection [@Shivakumar2025]. f-SMEMs, the maximal query intervals that occur at least f '
 'times, are k-MEMs in MONI [@Tatarnikov2022], f-MEMs in Gagie\'s terms [@Gagie2026] and (L, f)-SMEMs in ropebwt3 '
 'and minibwa, the two tools that compute them for whole read sets [@Li2024; @Li2026].'),
('p',
 'We present sporos, a library that computes all of these seeds exactly: minimizers and syncmers, an index of '
 'them for looking seeds up, MEMs, MUMs and frequency-capped MEMs between genomes, and SMEMs and f-SMEMs of reads. '
 'It follows each established tool\'s definitions and output, so that it can replace them, and it was faster '
 'than each of them on our benchmarks. Our contributions are:'),
('list', [
 '**Frequency-capped MEMs from the MEMs themselves** (Section 3.4): a MEM\'s string occurs in the reference once '
 'for every MEM whose query interval contains its own, and in the query once for every MEM whose reference '
 'interval does (Lemma 1), so a sweep over the MEMs that keeps the N largest interval ends counts the copies up '
 'to a cap N. This needs no suffix array, LCP or document array of the genomes.',
 '**MEMs found once from a sparse anchor index** (Section 3.3): sampled anchors as in copMEM [@Grabowski2019], '
 'each MEM reported by exactly one anchor, so that no pass removes duplicates, with probes streamed through a '
 'prefetching pipeline.',
 '**An SMEM search that skips unique stretches** (Section 3.5): once a match occurs once, its occurrence is '
 'located and the match extended by comparing the reference 32 bases at a time, instead of one FM-index step per '
 'base, with 16 reads searched together so that their memory accesses overlap.',
 '**SIMD sampled k-mers in order** (Section 3.2): minimizers and syncmers on 8 or 16 lanes with exact leftmost '
 'ties, written in sequence order by many threads, and a seed index probed by whole batches of reads.',
 '**One tested library** ([#fig:overview]): one definition per seed for the Mojo, C, Python and command-line '
 'interfaces, output identical to MUMmer4, copMEM2 and BWA, and every kernel checked against brute-force '
 'oracles.']),
('p',
 f'On one Intel Skylake-X processor sporos was the fastest tool for every seed we measured ([#fig:speedups]), '
 f'from {V["seed_t1_min"]} times faster than simd-minimizers on minimizers to {V["capped_max"]} times faster than '
 'Mumemto on frequency-capped MEMs.'),

('sec', 'Definitions'),
('p',
 'Sequences are strings over \\(\\{A, C, G, T\\}\\); any other symbol, such as N, matches nothing and ends every '
 'k-mer and match that reaches it. The reference \\(R\\) may hold several sequences; a match never crosses from one '
 'into the next. For a string \\(S\\), \\(S[i..j)\\) is its substring from position \\(i\\) up to but excluding '
 '\\(j\\), and \\(\\overline{S}\\) its reverse complement.'),
('p',
 '**Sampled k-mers.** Given a hash \\(h\\) of k-mers, a (w, k)-minimizer is the k-mer of smallest hash in a window '
 'of w consecutive k-mers, the leftmost one on ties; a sequence\'s minimizers are those of all its windows, each '
 'position reported once [@Roberts2004]. A *canonical* minimizer orders k-mers by a hash that gives a k-mer and its '
 'reverse complement the same value, so that both strands select the same positions. A k-mer is a *closed '
 'syncmer* if its smallest s-mer starts at offset 0 or \\(k - s\\), and an *open syncmer* at offset t if it starts '
 'at t [@Edgar2021]. Where we compare with simd-minimizers we follow its parameters: k is the length of the '
 'smaller k-mers ranked and w their number, so that a syncmer spans \\(k + w - 1\\) bases [@GrootKoerkamp2025simd].'),
('p',
 '**Exact matches between genomes.** A *MEM* of a query \\(Q\\) and the reference is a triple \\((r, q, \\ell)\\) with '
 '\\(R[r..r+\\ell) = Q[q..q+\\ell)\\) that extends neither left nor right, of length \\(\\ell\\) at least L. It is a '
 '*MUM* if its string occurs exactly once in \\(R\\) and once in \\(Q\\), and unique in the reference '
 '(MUMmer\'s `-mumreference`) if it occurs once in \\(R\\) [@Kurtz2004]. A MEM is *frequency-capped* at N if its string '
 'occurs at most N times in \\(R\\), and, in the stricter form, at most N times in \\(Q\\) as well; a floor M asks for at '
 'least M copies in \\(R\\). With N = 1 they are the MUMs. Matches of the reverse strand are those of \\(\\overline{Q}\\), '
 'whose positions count along \\(\\overline{Q}\\), as MUMmer reports them.'),
('p',
 '**Exact matches of reads.** An interval \\(Q[i..j)\\) of a read *occurs* if \\(Q[i..j)\\) or its reverse complement '
 'is a substring of \\(R\\); an *SMEM* is an interval of length at least L that occurs and is contained in no longer '
 'interval that occurs [@Li2012]. An *f-SMEM* replaces "occurs" by "occurs at least f times" [@Gagie2026]; f = 1 '
 'gives the SMEMs. A cap N drops the SMEMs that occur more than N times.'),

('sec', 'Methods'),
('sub', 'Overview'),
('p',
 'sporos is written in Mojo [@Mojo], which compiles through MLIR [@Lattner2021] to native code and exposes '
 'SIMD vectors of any width as ordinary types. Input is FASTA or FASTQ, plain, gzip or BGZF, inflated in parallel '
 'block by block; bases are coded in two bits, one per byte, and any other symbol as a fifth code. Three engines '
 'compute the seeds ([#fig:overview]), all multi-threaded by a small pool of POSIX threads, behind a Mojo '
 'interface, a command-line tool whose MEM output is MUMmer\'s and SMEM output `bwa fastmap`\'s, and C and Python '
 'bindings.'),
('figure', 'overview'),
('sub', 'Minimizers and syncmers'),
('p',
 'One kernel serves minimizers and syncmers. A sequence is split into 16 contiguous chunks on processors with '
 'AVX-512 and 8 otherwise, each scanned by one lane of a SIMD register, as simd-minimizers does [@GrootKoerkamp2025simd]. '
 'Each lane rolls a 32-bit ntHash [@Mohamadi2016] one base at a time; bases reach the lanes through a byte '
 'transpose, and each hash update is a single byte shuffle that looks up the base\'s seed. Each k-mer is keyed in '
 'one 32-bit word: the high bits of its hash, then its offset within a block of w consecutive k-mers. The minimum '
 'of the last w keys follows from prefix and suffix minima of these blocks, three `min` instructions per k-mer '
 'amortised, and within a block the offsets break ties to the leftmost k-mer, so that the result is exact. A '
 'k-mer covering any other symbol than A, C, G and T is keyed above every valid one. Positions equal to the '
 'previous window\'s are dropped in the loop with one vector compare, and every eight steps the surviving lanes '
 'are appended by walking their bit mask. Syncmers use the same windows over s-mers.'),
('p',
 'For many short reads, each read takes its own lane. For one long sequence on many threads, each thread takes '
 'small segments in order and copies its positions into the final list once the previous segment has published '
 'where its positions end, so that the result is one list in sequence order without a merge.'),
('p',
 '**Seed index.** The seed index holds the minimizers, syncmers or all k-mers (k at most 32) of a reference, keyed '
 'by canonical k-mer through an invertible mixing of its 2k bits. The high bits of the mix choose a bucket and '
 'the rest are kept as a fingerprint, so that the two identify the k-mer exactly whenever they hold 2k bits '
 'together; otherwise a candidate is checked against the reference. The table is built by a parallel radix sort. '
 'A batch of reads is seeded and probed through one software pipeline that hashes each seed, loads its bucket '
 'bounds, compares fingerprints and prefetches the next probes, crossing read boundaries without draining. '
 'Seeds may be skipped above an occurrence count, chosen as a quantile of the counts as in GraphAligner '
 '[@Rautiainen2020].'),
('sub', 'MEMs and MUMs'),
('p',
 'MEMs are found from short anchors, as in copMEM [@Grabowski2019]: the reference K-mers starting every \\(k_1\\) '
 'bases (K at most 32) are indexed, and the query K-mers starting every \\(k_2\\) bases probe them, with \\(k_1\\) '
 'and \\(k_2\\) coprime and \\(k_1 k_2 \\le L - K + 1\\). A MEM of length at least L then contains a pair of sampled '
 'positions on its diagonal, since among any \\(k_1 k_2\\) consecutive query positions one is sampled in both '
 'sequences (the Chinese remainder theorem). sporos chooses K, \\(k_1\\) and \\(k_2\\) for each L to minimise the '
 'estimated work, one insertion per sampled reference position and three for each probe. Each anchor hit is '
 'extended in both directions, 32 bases per comparison, as far as the match goes. Anchors of one MEM lie '
 '\\(k_1 k_2\\) apart on its diagonal, and only the anchor whose leftward extension does not reach that far '
 'reports the MEM: each MEM is found exactly once, with no pass to remove duplicates. The anchor table is built '
 'by a two-level radix sort with about four anchors per bucket, and probes stream through a pipeline that hashes, '
 'loads bucket bounds, compares fingerprints in a vector and prefetches the reference text the extension will '
 'read. MEMs are radix-sorted by query position.'),
('p',
 'A MEM is unique in the reference unless another MEM of the same query contains it in query coordinates, and '
 'unique in the query unless another contains it in reference coordinates (Lemma 1 below with N = 1). sporos '
 'therefore finds MUMs by two linear sweeps over the MEMs, one in query order and one after a parallel radix sort '
 'by reference position, each parallel over blocks of MEMs.'),
('sub', 'Frequency-capped MEMs'),
('p',
 '**Lemma 1.** Let \\(\\mathcal{M}\\) be the MEMs of length at least L of a query sequence \\(Q\\) against \\(R\\), and '
 '\\(m = (r, q, \\ell) \\in \\mathcal{M}\\). The occurrences of \\(Q[q..q+\\ell)\\) in \\(R\\) correspond one to one to '
 'the MEMs \\((r\', q\', \\ell\') \\in \\mathcal{M}\\) with \\(q\' \\le q\\) and \\(q + \\ell \\le q\' + \\ell\'\\), and '
 'the occurrences of \\(R[r..r+\\ell)\\) in \\(Q\\) to those with \\(r\' \\le r\\) and \\(r + \\ell \\le r\' + \\ell\'\\).'),
('p',
 '*Proof.* An occurrence of \\(Q[q..q+\\ell)\\) at position \\(p\\) of \\(R\\) is a match on diagonal \\(p - q\\), '
 'which extends left and right to a MEM whose query interval contains \\([q, q+\\ell)\\); it is at least \\(\\ell \\ge '
 'L\\) long, so it belongs to \\(\\mathcal{M}\\). Distinct occurrences lie on distinct diagonals and give distinct '
 'MEMs. Conversely, a MEM \\((r\', q\', \\ell\')\\) whose query interval contains \\(m\\)\'s holds an occurrence at '
 '\\(p = r\' + q - q\'\\). The second statement is symmetric. \\(\\square\\)'),
('p',
 'The copies of every MEM\'s string are therefore counts of containing intervals among the MEMs themselves. sporos '
 'counts them, up to the cap N, in one sweep over the MEMs of a query sequence ordered by start, longest first '
 'among equal starts. The sweep keeps the N largest interval ends seen so far in a sorted buffer, dropping those '
 'that end before the current start; an interval\'s count of containing intervals, capped at N, is the number of '
 'kept ends that reach its own end, found by binary search. With N = 1 the sweep is the MUM test. The query-side '
 'count is the same sweep over reference coordinates. The sweep is parallel over blocks: each block starts from '
 'the N largest ends of the blocks before it, which a first pass computes. Finding the MEMs and sorting them '
 'dominates; the sweep takes \\(O(\\log N)\\) time per MEM, plus the insertions into the buffer of N ends. A floor M '
 'counts up to \\(\\max(N, M)\\).'),
('p',
 'Mumemto instead builds the suffix array, LCP array and document array of the whole collection by prefix-free '
 'parsing and scans them [@Shivakumar2025]; it defines matches across all genomes at once and counts copies on '
 'both strands of each, so its matches are not exactly sporos\'s (Section 4.4).'),
('sub', 'SMEMs and f-SMEMs'),
('p',
 '**Index.** sporos builds an FMD-index [@Li2012]: the Burrows–Wheeler transform [@BurrowsWheeler1994] of the reference\'s '
 'runs of A, C, G and T and their reverse complements, each run ended by a separator, so that the interval of a '
 'pattern comes with that of its reverse complement and no match crosses a sequence end or a non-ACGT base. The '
 'suffix array is built in parallel by libsais [@libsais]. The transform is stored as in minibwa, in 64-byte '
 'blocks of 128 rows holding four cumulative counts and 128 two-bit symbols [@Li2026]; separators are stored as '
 'A and listed apart, and a flag marks the few blocks that hold one. Ranks of all four bases take three popcounts '
 'per word without branches. The index also holds the two-bit reference, suffix-array samples of 40 bits and the '
 'interval of every k-mer up to a length of 12, which starts every search with one lookup.'),
('p',
 '**Search.** SMEMs are found by the forward–backward algorithm of ropebwt3 and minibwa [@Gagie2024; @Li2024; @Li2026]: '
 'from each start, a backward search over L bases, a forward extension, and a backward pass from the match end '
 'that finds the next start. A match that occurs once is not followed further in the index: sporos locates its '
 'occurrence, which the output needs in any case, and extends it in both directions by comparing 32 bases at a '
 'time against the two-bit reference. Most of a read lies in such unique stretches, which other FM-index tools '
 'walk one base at a time. Sixteen reads are searched together, their steps interleaved round-robin, each step '
 'prefetching the block the next will read, and the occurrences of repeated SMEMs are located in batches in the '
 'same way. f-SMEMs replace every test that an interval is empty by a test that it holds fewer than f rows, and a '
 'cap N drops SMEMs with more than N. Threads take slices of reads.'),
('p',
 'The same SMEMs also follow from the MEM engine without an index: computed on both strands, the SMEMs are the '
 'MEM query intervals that no other one contains. This suits a single long query, such as a genome.'),
('sub', 'Validation'),
('p',
 'Every kernel is checked in continuous integration against brute-force oracles that share no code with the '
 'library: minimizers and syncmers against a direct scan of every window, MEMs, MUMs and capped MEMs against '
 'enumeration of all pairs of starts on repetitive inputs, and SMEMs and f-SMEMs against string search for each '
 'interval. The C and Python bindings are checked against the Mojo library\'s output. On the benchmarks below, '
 'every output of sporos equalled MUMmer4\'s, its MEMs equalled copMEM2\'s, and ten threads gave the same output as '
 'one.'),

('sec', 'Results'),
('sub', 'Experimental setup'),
('p',
 f'All measurements were taken on an {V["machine_cpu"]} (Skylake-X, AVX-512, ten cores, hyper-threading off) with '
 '91 GB of memory, under Ubuntu Linux. '
 'sporos, simd-minimizers, minimap2, copMEM2 and strobealign were compiled for the native instruction set, minibwa '
 'with its own build, and MUMmer4, Mumemto, BWA and ropebwt3 installed from Bioconda. The genome pairs are '
 '*Escherichia coli* K-12 MG1655 '
 'against O157:H7 Sakai [@Blattner1997; @Hayashi2001], human chromosome 21 of GRCh38 [@Schneider2017] against '
 'chimpanzee chromosome 21 of the complete Clint assembly [@Yoo2025], and chromosome 1 of GRCh38 against that of '
 'T2T-CHM13 [@Nurk2022]. Reads were simulated from GRCh38 chromosome 21: one million of 150 bp and 4,000 of 10 kbp, '
 'each base substituted with probability 0.002 and half the reads reverse-complemented. Each command was run '
 'three times and its best wall time kept (once when it exceeded 60 s), after waiting for the machine to be idle; '
 'every tool wrote its output to a file. Supplementary Section S1 gives the versions and commands.'),
('sub', 'Minimizers and syncmers'),
('p',
 f'sporos computed minimizers and closed syncmers {V["seed_t1_min"]}–{V["seed_t1_max"]} times faster than '
 f'simd-minimizers on one thread, over every scheme on CHM13 chromosome 1 and on random DNA ([#tab:seeds]; '
 f'Supplementary Table S2), at {V["seed_chr1_ns"]} ns per base for canonical (21, 11)-minimizers of chromosome 1, '
 f'or {V["seed_chr1_gbps"]} Gbp per second. Restricted to 8 lanes, the vector width simd-minimizers uses, it was '
 f'{V["seed8_min"]}–{V["seed8_max"]} times faster on chromosome 1, and the 16-lane build was '
 f'{V["seed16_vs_8_min"]}–{V["seed16_vs_8_max"]} times faster than the 8-lane one, where simd-minimizers found its '
 f'own kernel slower with AVX-512 [@GrootKoerkamp2025simd]. minimap2\'s `mm_sketch` was '
 f'{span(V["seed_mm2_min"], V["seed_mm2_max"])} times slower than sporos. On ten threads the race is closer: sporos was '
 f'{V["seed_t10_phrase"]} simd-minimizers, whose harness keeps each thread\'s '
 'positions apart and only counts them, while sporos returns one list in sequence order.'),
('table', 'seeds'),
('sub', 'Seed-index lookups'),
('p',
 f'Seeding 150 bp reads and looking their minimizers up, sporos processed {V["lookup_sporos"]} Mbp of reads per '
 f'second on one thread, {V["lookup_vs_mm2"]} times minimap2\'s {V["lookup_mm2"]} Mbp/s for the same minimizers '
 f'and occurrence limit, and {V["lookup_vs_sa"]} times strobealign\'s rate for its randstrobes ([#tab:lookups]). '
 f'Probed one read per call, sporos was still {V["lookup_one_vs_mm2"]} times faster than minimap2.'),
('table', 'lookups'),
('sub', 'MEMs and MUMs'),
('p',
 f'sporos found every MEM of the three genome pairs {V["copmem2_t1_min"]}–{V["copmem2_t1_max"]} times faster than '
 f'copMEM2 on one thread, and {V["copmem2_chr1_min"]}–{V["copmem2_chr1_max"]} times faster on the pair of human '
 f'chromosomes 1 ([#tab:mems]; Supplementary Table S3). copMEM2 and MUMmer4 do not divide a single query sequence '
 f'among threads, so on ten threads sporos was {V["copmem2_chr1_t10_min"]}–{V["copmem2_chr1_t10_max"]} times '
 f'faster than copMEM2 on chromosome 1, finding its {V["checks_max_mems"]} million MEMs of at least 50 bp in {V["chr1_mems_t10"]} s. '
 f'Against MUMmer4, sporos was {V["mummer_min"]}–{V["mummer_max"]} times faster on one thread for MEMs, MUMs and '
 'MEMs unique in the reference (Supplementary Table S4); the smallest margin is for MUMs of the chromosome 1 pair, '
 f'where MUMmer4 never enumerates the {V["checks_max_mems"]} million MEMs that sporos filters. In all '
 f'{V["checks_total"]} comparisons, every output of sporos equalled MUMmer4\'s or copMEM2\'s match for match. '
 f'sporos used {V["chr1_rss_gb"]} GB at most on chromosome 1 (Supplementary Table S5).'),
('table', 'mems'),
('sub', 'Frequency-capped MEMs'),
('p',
 f'sporos computed frequency-capped MEMs {V["capped_min"]}–{V["capped_max"]} times faster than Mumemto on one '
 'thread, and in a small fraction of its memory ([#tab:capped]). The comparison is between the closest settings of '
 'tools whose definitions differ: Mumemto treats the two genomes as one collection, counts copies on both strands '
 'of each and calls a match maximal if its copies cannot all be extended together, whereas sporos reports the '
 'MEMs of the query against the reference, MUMmer\'s pairs of positions, and counts copies in the reference as '
 'given. On the *E. coli* pair the outputs largely agree: sporos reported '
 f'{span(V["share_mm_ecoli_min"], V["share_mm_ecoli_max"])}% of Mumemto\'s forward-strand position pairs and Mumemto '
 f'{span(V["share_sp_ecoli_min"], V["share_sp_ecoli_max"])}% of sporos\'s. On chromosome 21 the shares were '
 f'{span(V["share_mm_chr21_min"], V["share_mm_chr21_max"])}% and {span(V["share_sp_chr21_min"], V["share_sp_chr21_max"])}%: '
 'in repeats, two copies of a string that is maximal across the collection often extend further as a pair, and '
 'sporos then reports that longer match instead ([#tab:capped]).'),
('table', 'capped'),
('sub', 'SMEMs and f-SMEMs'),
('p',
 f'On one thread, sporos found the SMEMs of 150 bp reads {V["smem_150bp_mb"]} times faster than minibwa, '
 f'{V["smem_150bp_bwa"]} times faster than `bwa fastmap` and {V["smem_150bp_rb"]} times faster than ropebwt3 '
 f'without positions, and those of 10 kbp reads {V["smem_10kb_mb"]}, {V["smem_10kb_bwa"]} and '
 f'{V["smem_10kb_rb"]} times faster ([#tab:smems]). For f-SMEMs it was {V["fsmem2_mb"]} times faster than minibwa '
 f'at f = 2 and {V["fsmem5_mb"]} times at f = 5, and {V["fsmem2_rb"]} and {V["fsmem5_rb"]} times faster than '
 f'ropebwt3. Asking for two copies made the search {V["fsmem_cost"]} times slower than for SMEMs, since a match '
 'of one copy no longer ends the search. The SMEMs of sporos equalled those of `bwa fastmap` except where BWA, '
 'which replaces each N of the reference by a random base, matched through one (Supplementary Section S3).'),
('table', 'smems'),
('figure', 'speedups'),

('sec', 'Discussion'),
('p',
 'sporos computes the seeds that read mappers, genome aligners and pangenome tools start from, exactly as the '
 'established tools define them and faster than each of them on our benchmarks. Its frequency-capped MEMs show '
 'that a property of a set of matches, how often each match\'s string occurs, can follow from the matches alone: '
 'counting containing intervals replaces an index of every genome, which makes capped matches between two '
 'chromosomes a matter of seconds.'),
('p',
 'One library also removes the conversions between tools. A pipeline can index a reference once and ask it for '
 'MEMs, MUMs or SMEMs, from Mojo, C or Python, with one coordinate system and one set of occurrence limits, and '
 'the command-line tool writes MUMmer\'s and `bwa fastmap`\'s formats for existing scripts.'),
('p',
 '**Limitations.** All measurements come from one x86-64 processor; the 8-lane rows of [#tab:seeds] show the '
 'kernel at the width of AVX2, but other processors, including ARM, were not measured. The comparison with '
 'Mumemto covers pairs of genomes, not the collections of many genomes Mumemto was designed for, and the two '
 'tools\' definitions differ. Positions are 32-bit, so each sequence and the concatenated reference must be shorter '
 'than \\(2^{32}\\) bases, and building an FM-index takes about ten bytes of memory per reference base. '
 f'{V["seed_t10_limit"]}.'),
]

BACK = [
    ("Author contributions", "Konstantinos Kyriakidis (Conceptualization [lead], Methodology [lead], Software [lead], "
     "Validation [lead], Formal analysis [lead], Investigation [lead], Data curation [lead], Visualization [lead], "
     "Writing—original draft [lead]), Benedict Paten (Supervision [lead], Writing—review & editing [equal])."),
    ("Supplementary data", "Supplementary data are available at *Bioinformatics* online."),
    ("Conflict of interest", "None declared."),
    ("Funding", "This work received no specific funding."),
    ("Data availability", "The software and the benchmark scripts are available at "
     "https://github.com/kokyriakidis/sporos; `benchmarks/competitors/` sets up every tool at a pinned version, "
     "downloads the genomes, simulates the reads and runs each benchmark, and its `results/` directory holds the "
     "raw measurements behind every table."),
]

FIGURES = {
"overview": "The sporos library. Input is parsed and coded once; three engines compute sampled k-mers and the seed "
            "index that looks them up (top), MEMs, MUMs and frequency-capped MEMs between genomes (middle), and SMEMs "
            "and f-SMEMs of reads (bottom), behind Mojo, command-line, C and Python interfaces.",
"speedups": "Speedup of sporos over the fastest other tool for each task, on one thread (filled) and ten "
            "(hollow), on a log scale; the dashed line marks parity, and the values are at the right (Tables 1 to 5). "
            "Each tool runs at its best on that many threads: copMEM2, MUMmer4, minibwa, BWA and Mumemto do not gain "
            "from more threads on these inputs, and the seed-lookup harnesses run on one.",
}
