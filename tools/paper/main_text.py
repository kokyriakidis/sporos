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
  "As sequencing data and reference collections grow, from complete genomes to pangenomes, finding seeds, the "
  "exact matches between a query and a reference from which reads are mapped and genomes are aligned, has become a "
  "central step of sequence analysis. Many kinds of seeds are in use, from sampled k-mers such as minimizers to "
  "maximal exact matches, and each is computed by its own specialized tool, with its own index and output format. "
  "And as processors gain wider vector units and deeper memory hierarchies, the speed of seeding depends "
  "increasingly on how well an implementation uses each core."),
 ("Results",
  "We present sporos, one library that computes all of these seeds exactly, built around each core's vector "
  "units and memory parallelism. Its key contribution counts the copies of a MEM's string from the MEMs "
  "themselves: a string occurs in the reference once for every MEM whose query interval contains its own, so "
  "frequency-capped MEMs need no additional index. On a single core of an Intel Skylake-X processor, sporos "
  f"computed minimizers and syncmers {V['seed_t1_min']}–{V['seed_t1_max']} times faster than simd-minimizers, looked "
  f"up seeds {V['lookup_vs_mm2']} times faster than minimap2, found MEMs between two human chromosomes 1 "
  f"{V['copmem2_chr1_min']}–{V['copmem2_chr1_max']} times faster than copMEM2, MEMs and MUMs "
  f"{V['mummer_min']}–{V['mummer_max']} times faster than MUMmer4 with identical output, frequency-capped MEMs "
  f"{V['capped_min']}–{V['capped_max']} times faster than Mumemto, and SMEMs of short reads "
  f"{V['smem_150bp_mb']} times faster than minibwa."),
 ("Availability and implementation",
  "sporos is free software under the MIT license, written in Mojo with C, Python and command-line interfaces, "
  "at https://github.com/kokyriakidis/sporos."),
]

BODY = [
('sec', 'Introduction'),
('p',
 'Most analyses that compare DNA sequences begin by finding seeds, exact matches that indicate where a query and a '
 'reference may align. Read mappers chain seeds before alignment [@Li2013; @Li2018], whole-genome aligners anchor '
 'on matches unique to both genomes [@Delcher1999; @Marcais2018], and pangenome tools summarize shared sequence '
 'through matches common to many genomes [@Shivakumar2025]. Two families of seeds are in wide use. *Sampled '
 'k-mers* retain a content-dependent subset of k-mers, so that copies of a region retain the same ones: minimizers, '
 'the k-mer of smallest hash in each window of w consecutive k-mers [@Schleimer2003; @Roberts2004], and syncmers, '
 'k-mers whose smallest s-mer lies at a fixed offset [@Edgar2021]. *Maximal exact matches* (MEMs) extend until the '
 'sequences differ; maximal unique matches (MUMs) anchor whole-genome alignment [@Delcher1999], and super-maximal '
 'exact matches (SMEMs), the MEMs of a read not contained in a longer one, seed BWA-MEM [@Li2012; @Li2013].'),
('p',
 'Seeding is a substantial cost. On short reads, single-threaded BWA-MEM spends 40–60% of its time finding SMEMs '
 'and locating their occurrences [@Vasimuddin2019], and dedicated hardware has been proposed to accelerate this '
 'step [@Subramaniyan2021]. The inputs keep growing, from complete telomere-to-telomere assemblies [@Nurk2022] to '
 'pangenomes of dozens of individuals [@Liao2023], while the throughput of a processor core grows mainly through '
 'wider vector units and greater memory-level parallelism.'),
('p',
 'Each kind of seed is computed by a different tool on a different data structure. Minimizers are computed within '
 'read mappers such as minimap2 [@Li2018], syncmers and strobemers within strobealign [@Sahlin2021; @Sahlin2022], '
 'and simd-minimizers computes random minimizers with SIMD instructions [@GrootKoerkamp2025simd]. MEMs and MUMs '
 'between genomes are found with suffix trees and suffix arrays in MUMmer [@Kurtz2004; @Marcais2018], with sparse '
 'suffix arrays [@Khan2009; @Vyverman2013] and with hash tables of sampled k-mers in E-MEM and copMEM '
 '[@Khiste2015; @Grabowski2019; @Grabowski2023]. SMEMs are found with FM-indexes [@Ferragina2000] in BWA '
 '[@Li2013], bwa-mem2 [@Vasimuddin2019], ropebwt3 [@Li2024] and minibwa [@Li2026]. Newer seed types are supported '
 'by few tools. Mumemto relaxes multi-MUMs to matches that occur at most a given number of times, using the suffix '
 'array, LCP array and document array of the whole collection [@Shivakumar2025]. f-SMEMs, the maximal query '
 'intervals that occur at least f times, were introduced as k-MEMs [@Tatarnikov2022] and f-MEMs [@Gagie2026], and '
 'ropebwt3 and minibwa compute them for read sets [@Li2024; @Li2026]. Most of these tools parallelize across query '
 'sequences and process each sequence on one core, as a program built on a library would. We therefore compare the '
 'tools primarily on a single core, and also report ten threads, which matter most when a single long query, such '
 'as a chromosome, cannot be divided across sequences.'),
('p',
 'Here we present sporos, a library that computes all of these seeds exactly: minimizers and syncmers, a seed index '
 'for lookups, MEMs, MUMs and frequency-capped MEMs between genomes, and SMEMs and f-SMEMs of reads. It reproduces '
 'the definitions and output formats of the established tools, and on a single core it was faster than each of '
 'them. Our contributions are:'),
('list', [
 '**Frequency-capped MEMs without a collection index** (Section 3.4). The occurrences of a MEM\'s string in the '
 'reference correspond one to one to the MEMs whose query interval contains it (Lemma 1), so a sweep that retains the '
 'N largest interval ends counts occurrences up to a cap N.',
 '**Duplicate-free MEM enumeration** (Section 3.3). With the sparse sampling of copMEM [@Grabowski2019], each MEM is '
 'reported by exactly one anchor, and anchor probes stream through a prefetching pipeline.',
 '**SMEM search with direct extension of unique matches** (Section 3.5). A match that occurs once is located and '
 'extended by comparing 32 bases at a time against the reference, and 16 reads are searched in interleaved fashion '
 'so that their memory accesses overlap.',
 '**Vectorized k-mer sampling** (Section 3.2). Minimizers and syncmers are computed on 8 or 16 SIMD lanes with exact '
 'leftmost tie-breaking and returned in sequence order from any number of threads; a seed index answers batched '
 'lookups through one prefetching pipeline.',
 '**A tested library** ([#fig:overview]). One implementation per seed type serves the Mojo, C, Python and '
 'command-line interfaces; outputs equal those of MUMmer4 and copMEM2, and those of BWA except at ambiguous '
 'reference bases, and every kernel is tested against a brute-force oracle.']),
('p',
 f'On a single core of an Intel Skylake-X processor, sporos was faster than every other tool for every seed type we '
 f'measured ([#fig:speedups]), from {V["seed_t1_min"]} times faster than simd-minimizers for minimizers to '
 f'{V["capped_max"]} times faster than Mumemto for frequency-capped MEMs.'),

('sec', 'Definitions'),
('p',
 'Sequences are strings over \\(\\{A, C, G, T\\}\\); any other symbol, such as N, matches nothing, so no k-mer or match '
 'spans it. The reference \\(R\\) may contain several sequences, and matches do not cross sequence boundaries. '
 '\\(S[i..j)\\) denotes the substring of \\(S\\) from position \\(i\\) to \\(j - 1\\), and \\(\\overline{S}\\) its '
 'reverse complement.'),
('p',
 '**Sampled k-mers.** For a hash \\(h\\) on k-mers, the (w, k)-minimizer of a window of w consecutive k-mers is its '
 'k-mer of minimum hash, with ties broken towards the leftmost; the minimizers of a sequence are those of all its '
 'windows, each reported once [@Roberts2004]. Canonical minimizers use a hash that is invariant under reverse '
 'complementation. A k-mer is a *closed syncmer* if its minimal s-mer starts at offset 0 or \\(k - s\\), and an '
 '*open syncmer* with offset t if it starts at t [@Edgar2021]. For comparison with simd-minimizers we use its '
 'parameterization, in which a syncmer is defined by the length k of the ranked substrings and their number w, and '
 'spans \\(k + w - 1\\) bases [@GrootKoerkamp2025simd].'),
('p',
 '**Exact matches between genomes.** A *MEM* of a query sequence \\(Q\\) and \\(R\\) is a triple \\((r, q, \\ell)\\) '
 'with \\(R[r..r+\\ell) = Q[q..q+\\ell)\\) and \\(\\ell \\ge L\\) that cannot be extended to the left or right. It is a '
 '*MUM* if its string occurs exactly once in \\(R\\) and once in \\(Q\\), and *reference-unique* (MUMmer\'s '
 '`-mumreference`) if it occurs once in \\(R\\) [@Kurtz2004]. It is *frequency-capped* at N if its string occurs at most '
 'N times in \\(R\\) or, in the strict form, at most N times in both \\(R\\) and \\(Q\\); a floor M additionally requires '
 'at least M occurrences in \\(R\\). N = 1 yields the MUMs. Reverse-strand matches are those of \\(\\overline{Q}\\), with '
 'positions along \\(\\overline{Q}\\), as in MUMmer.'),
('p',
 '**Exact matches of reads.** An interval \\(Q[i..j)\\) of a read *occurs* if it or its reverse complement is a '
 'substring of \\(R\\). An *SMEM* is an interval of length at least L that occurs and is not contained in a longer '
 'occurring interval [@Li2012]. An *f-SMEM* replaces "occurs" by "occurs at least f times" [@Tatarnikov2022; '
 '@Gagie2026]; f = 1 gives the SMEMs, and a cap N removes SMEMs that occur more than N times.'),

('sec', 'Methods'),
('sub', 'Overview'),
('p',
 'sporos is implemented in Mojo [@Mojo], which compiles through MLIR [@Lattner2021] to native code and exposes SIMD '
 'vectors as first-class types. Input may be FASTA or FASTQ, plain, gzip- or BGZF-compressed, with BGZF blocks '
 'decompressed in parallel. Bases are stored as 2-bit codes, one per byte, with a fifth code for other symbols. '
 'Three engines compute the seeds ([#fig:overview]); all are multithreaded with POSIX threads and exposed through a '
 'Mojo API, C and Python bindings and a command-line tool that writes MUMmer\'s MEM format and `bwa fastmap`\'s SMEM '
 'format.'),
('figure', 'overview'),
('sub', 'Minimizers and syncmers'),
('p',
 'A single kernel computes minimizers and syncmers. As in simd-minimizers [@GrootKoerkamp2025simd], the sequence is '
 'divided into contiguous chunks, one per SIMD lane: 16 with AVX-512 and 8 otherwise. Each lane rolls a 32-bit '
 'ntHash [@Mohamadi2016] with a rotation step of 7 instead of 1, which decorrelates the high bits of consecutive '
 'hashes; the canonical hash is the sum of the forward and reverse-complement hashes. Bases reach the lanes through '
 'a byte transpose, and each hash update is one byte-shuffle table lookup. Each k-mer is represented by a 32-bit key '
 'holding the high bits of its hash followed by its offset within a block of w consecutive k-mers, so that the '
 'sliding-window minimum, computed from prefix and suffix minima of these blocks with three `min` operations per '
 'k-mer (amortized), resolves ties to the leftmost k-mer exactly. k-mers containing other symbols receive a key above '
 'every valid one. Repeated minimizers are removed in the loop with one vector comparison, and surviving positions '
 'are collected every eight steps from a bit mask. Syncmers apply the same window minimum to s-mers.'),
('p',
 'For batches of short reads, each read occupies its own lane. For a long sequence on several threads, threads '
 'process short segments in order, and each copies its positions into the output once its predecessor has published '
 'where its positions end, yielding a single position list in sequence order without a merge.'),
('p',
 '**Seed index.** The seed index stores the minimizers, syncmers or all k-mers (k ≤ 32) of a reference, keyed by '
 'canonical k-mer through an invertible mixing of its 2k bits. The high bits of the mixed key select a bucket and the '
 'remaining bits form a fingerprint; when bucket and fingerprint together hold 2k bits they identify the k-mer '
 'exactly, and otherwise candidates are verified against the reference. The table is built by a parallel radix sort. '
 'A batch of reads passes through one software pipeline that computes seeds, loads bucket bounds, compares '
 'fingerprints and prefetches subsequent probes across read boundaries. Seeds above an occurrence threshold, chosen '
 'as a quantile of the count distribution as in GraphAligner [@Rautiainen2020], can be skipped.'),
('sub', 'MEMs and MUMs'),
('p',
 'Following copMEM [@Grabowski2019], reference K-mers (K ≤ 32) sampled every \\(k_1\\) positions are indexed and '
 'query K-mers sampled every \\(k_2\\) positions are looked up, with \\(k_1\\) and \\(k_2\\) coprime and \\(k_1 k_2 \\le '
 'L - K + 1\\). By the Chinese remainder theorem, any \\(k_1 k_2\\) consecutive positions of a diagonal include one '
 'sampled in both sequences, so every MEM of length at least L contains an anchor. K, \\(k_1\\) and \\(k_2\\) are '
 'chosen for each L to minimize the estimated work, counting one unit per sampled reference position and three per '
 'probe. Each anchor hit is extended in both directions, 32 bases per comparison. Because the anchors of a MEM lie '
 '\\(k_1 k_2\\) apart along its diagonal, only the anchor whose leftward extension is shorter than \\(k_1 k_2\\) reports '
 'it, so each MEM is found exactly once without deduplication. The anchor table is built by a two-level radix sort '
 'with about four anchors per bucket, and probes are pipelined through hashing, bucket lookup, vectorized fingerprint '
 'comparison and prefetching of the reference text to be extended. MEMs are radix-sorted by query position.'),
('p',
 'By Lemma 1 (Section 3.4) with N = 1, a MEM is unique in the reference unless another MEM of the same query sequence '
 'contains it in query coordinates, and unique in the query unless another contains it in reference coordinates. '
 'MUMs therefore follow from two linear sweeps, one in query order and one after a parallel radix sort by reference '
 'position, each parallelized over blocks of MEMs.'),
('sub', 'Frequency-capped MEMs'),
('p',
 '**Lemma 1.** Let \\(\\mathcal{M}\\) be the MEMs of length at least L of a query sequence \\(Q\\) against \\(R\\), and '
 '\\(m = (r, q, \\ell) \\in \\mathcal{M}\\). The occurrences of \\(Q[q..q+\\ell)\\) in \\(R\\) correspond one to one to '
 'the MEMs \\((r\', q\', \\ell\') \\in \\mathcal{M}\\) with \\(q\' \\le q\\) and \\(q + \\ell \\le q\' + \\ell\'\\), and '
 'the occurrences of \\(R[r..r+\\ell)\\) in \\(Q\\) to those with \\(r\' \\le r\\) and \\(r + \\ell \\le r\' + \\ell\'\\).'),
('p',
 '*Proof.* An occurrence of \\(Q[q..q+\\ell)\\) at position \\(p\\) of \\(R\\) is a match on diagonal \\(p - q\\), which '
 'extends in both directions to a MEM whose query interval contains \\([q, q+\\ell)\\); its length is at least '
 '\\(\\ell \\ge L\\), so it belongs to \\(\\mathcal{M}\\). Distinct occurrences lie on distinct diagonals and give '
 'distinct MEMs. Conversely, a MEM \\((r\', q\', \\ell\')\\) whose query interval contains that of \\(m\\) holds an '
 'occurrence at \\(p = r\' + q - q\'\\). The second statement is symmetric. \\(\\square\\)'),
('p',
 'Occurrence counts thus reduce to counting containing intervals among the MEMs. sporos processes the MEMs of each '
 'query sequence in order of start position, handling equal starts together, and maintains the N largest interval '
 'ends seen so far in a sorted buffer, discarding ends at or before the current start. The number of intervals '
 'containing an interval, capped at N, is the number of retained ends at or beyond its own end, found by binary '
 'search; with N = 1 this is the MUM test. Query-side counts use the same sweep in reference coordinates. The sweep is '
 'parallelized over blocks of MEMs, each block starting from the N largest ends of the preceding blocks, computed in '
 'a first pass. Each MEM costs \\(O(\\log N)\\) for the search and at most \\(O(N)\\) for an insertion into the buffer. '
 'A floor M is handled by counting up to \\(\\max(N, M)\\).'),
('p',
 'Mumemto instead traverses the suffix array, LCP array, BWT and document array of the whole collection, built by '
 'prefix-free parsing [@Shivakumar2025]. It defines matches over all genomes jointly and counts occurrences on both '
 'strands, so its output differs from that of sporos by definition (Section 4.5).'),
('sub', 'SMEMs and f-SMEMs'),
('p',
 '**Index.** sporos builds an FMD-index [@Li2012]: the Burrows–Wheeler transform [@BurrowsWheeler1994] of the maximal '
 'ACGT runs of the reference and their reverse complements, each terminated by a separator, so that the interval of '
 'a pattern is obtained together with that of its reverse complement and no match crosses a sequence end or a '
 'non-ACGT base. The suffix array is constructed in parallel with libsais [@libsais]. As in minibwa [@Li2026], the '
 'transform is stored in 64-byte blocks of 128 rows, each holding four cumulative counts and 128 two-bit symbols; '
 'separators are stored as A and listed separately, with a flag on the few blocks that contain one. The ranks of all '
 'four bases are computed without branches using three population counts per word. The index also stores the 2-bit '
 'reference, 40-bit suffix-array samples and the interval of every k-mer for a k of up to 12, chosen by reference '
 'size, so that each search starts with one table lookup.'),
('p',
 '**Search.** SMEMs are found with the forward–backward algorithm [@Gagie2024] used by ropebwt3 and minibwa '
 '[@Li2024; @Li2026]: from each start position, a backward search over L bases, a forward extension and a backward '
 'pass from the end of the match that determines the next start. When a match occurs only once, sporos stops '
 'traversing the index: it locates the occurrence, which the output requires anyway, and extends the match in both '
 'directions by comparing 32 bases at a time against the 2-bit reference. Sixteen reads are searched concurrently, '
 'their steps interleaved round-robin and each step prefetching the block its successor will access; occurrences of '
 'repeated SMEMs are located in batches in the same way. For f-SMEMs, each test for an empty interval becomes a test '
 'for fewer than f rows, and a cap N removes SMEMs with more than N occurrences. Threads process slices of reads.'),
('p',
 'The MEM engine yields the same SMEMs without an index: computed on both strands, the SMEMs are the MEM query '
 'intervals not contained in another, which suits single long queries such as whole chromosomes.'),
('sub', 'Validation'),
('p',
 'All kernels are tested in continuous integration against brute-force oracles that share no code with the library: '
 'minimizers and syncmers against exhaustive window scans, MEMs, MUMs and frequency-capped MEMs against enumeration '
 'of all pairs of start positions on repetitive sequences, and SMEMs and f-SMEMs against direct string search. The C '
 'and Python bindings are tested against the Mojo library.'),

('sec', 'Results'),
('sub', 'Experimental setup'),
('p',
 f'Benchmarks ran on an {V["machine_cpu"]} (Skylake-X, AVX-512, ten cores) with turbo boost and hyper-threading '
 'disabled and 91 GB of memory, under Ubuntu Linux. sporos, simd-minimizers, minimap2, copMEM2 and strobealign were '
 'compiled for the native instruction set, minibwa with its default build, and MUMmer4, Mumemto, BWA and ropebwt3 '
 'were installed from Bioconda. The genome pairs were *Escherichia coli* K-12 MG1655 and O157:H7 Sakai '
 '[@Blattner1997; @Hayashi2001], human (GRCh38 [@Schneider2017]) and chimpanzee (Clint [@Yoo2025]) chromosome 21, '
 'and GRCh38 and T2T-CHM13 [@Nurk2022] chromosome 1. Reads were simulated from GRCh38 chromosome 21: 1,000,000 of '
 '150 bp and 4,000 of 10 kbp, each base replaced by a random base with probability 0.002 and half of the reads '
 'reverse-complemented. Unless stated otherwise, each command was run three times (once if it exceeded 60 s) after '
 'the machine became idle, and the best wall time is reported; MEM and SMEM tools wrote their output to files, '
 'except ropebwt3, whose output was discarded. Supplementary Section S1 lists versions and commands.'),
('sub', 'Minimizers and syncmers'),
('p',
 f'On one thread, sporos computed minimizers and closed syncmers {V["seed_t1_min"]}–{V["seed_t1_max"]} times faster '
 'than simd-minimizers across all schemes on CHM13 chromosome 1 and on random DNA ([#tab:seeds]; Supplementary '
 f'Table S2); canonical (21, 11)-minimizers of chromosome 1 took {V["seed_chr1_ns"]} ns per base '
 f'({V["seed_chr1_gbps"]} Gbp/s). With 8 lanes, the vector width of simd-minimizers, sporos was '
 f'{V["seed8_min"]}–{V["seed8_max"]} times faster on chromosome 1, and the 16-lane AVX-512 build was another '
 f'{V["seed16_vs_8_min"]}–{V["seed16_vs_8_max"]} times faster, whereas simd-minimizers reported a slowdown with '
 f'AVX-512 [@GrootKoerkamp2025simd]. minimap2\'s `mm_sketch` was {span(V["seed_mm2_min"], V["seed_mm2_max"])} times '
 f'slower than sporos. On ten threads, sporos was {V["seed_t10_phrase"]} simd-minimizers; our harness for '
 'simd-minimizers only counts each thread\'s positions, whereas sporos returns one ordered list.'),
('table', 'seeds'),
('sub', 'Seed-index lookups'),
('p',
 f'Seeding 150 bp reads and looking up their minimizers, sporos processed {V["lookup_sporos"]} Mbp per second on one '
 f'thread, {V["lookup_vs_mm2"]} times the throughput of minimap2 ({V["lookup_mm2"]} Mbp/s) for the same minimizers '
 f'and occurrence cap, and {V["lookup_vs_sa"]} times that of strobealign for its randstrobes ([#tab:lookups]). With '
 f'one read per call, sporos remained {V["lookup_one_vs_mm2"]} times faster than minimap2.'),
('table', 'lookups'),
('sub', 'MEMs and MUMs'),
('p',
 f'On one thread, sporos found all MEMs of the three genome pairs {V["copmem2_t1_min"]}–{V["copmem2_t1_max"]} times '
 f'faster than copMEM2, and {V["copmem2_chr1_min"]}–{V["copmem2_chr1_max"]} times faster on the chromosome 1 pair '
 '([#tab:mems]; Supplementary Table S3). As neither copMEM2 nor MUMmer4 divides a single query sequence among '
 f'threads, sporos was {V["copmem2_chr1_t10_min"]}–{V["copmem2_chr1_t10_max"]} times faster than copMEM2 on ten '
 f'threads, finding the {V["checks_max_mems"]} million MEMs of at least 50 bp between the chromosomes 1 in '
 f'{V["chr1_mems_t10"]} s. sporos was {V["mummer_min"]}–{V["mummer_max"]} times faster than MUMmer4 for MEMs, MUMs '
 'and reference-unique MEMs (Supplementary Table S4); the smallest margin was for MUMs of the chromosome 1 pair, for '
 f'which sporos first enumerates all {V["checks_max_mems"]} million MEMs. All {V["checks_total"]} output checks '
 'passed: the outputs of sporos equalled those of MUMmer4 and copMEM2, and its ten-thread outputs its one-thread '
 f'outputs. Peak memory on chromosome 1 was at most {V["chr1_rss_max_gb"]} GB (Supplementary Table S5).'),
('table', 'mems'),
('sub', 'Frequency-capped MEMs'),
('p',
 f'sporos computed frequency-capped MEMs {V["capped_min"]}–{V["capped_max"]} times faster than Mumemto on one '
 f'thread, using {V["capped_mem_pct_min"]}–{V["capped_mem_pct_max"]}% of its memory ([#tab:capped]). The definitions '
 'differ, so the comparison is between the closest settings: Mumemto treats both genomes as one collection, counts '
 'occurrences on both strands and considers a match maximal if its occurrences cannot all be extended together, '
 'whereas sporos reports MUMmer-style position pairs of the query against the reference and counts occurrences in '
 'the reference as given. On the *E. coli* pair, sporos reported '
 f'{span(V["share_mm_ecoli_min"], V["share_mm_ecoli_max"])}% of the forward-strand position pairs of Mumemto, and '
 f'Mumemto {span(V["share_sp_ecoli_min"], V["share_sp_ecoli_max"])}% of those of sporos. On chromosome 21 these '
 f'shares were {span(V["share_mm_chr21_min"], V["share_mm_chr21_max"])}% and '
 f'{span(V["share_sp_chr21_min"], V["share_sp_chr21_max"])}%, which we attribute to the two notions of maximality: in '
 'repeats, two copies of a string that is maximal across the collection can extend further as a pair, and sporos '
 'then reports the longer match.'),
('table', 'capped'),
('sub', 'SMEMs and f-SMEMs'),
('p',
 f'On one thread, sporos found the SMEMs of 150 bp reads {V["smem_150bp_mb"]} times faster than minibwa, '
 f'{V["smem_150bp_bwa"]} times faster than `bwa fastmap` and {V["smem_150bp_rb"]} times faster than ropebwt3 without '
 f'position output, and those of 10 kbp reads {V["smem_10kb_mb"]}, {V["smem_10kb_bwa"]} and {V["smem_10kb_rb"]} '
 f'times faster ([#tab:smems]). For f-SMEMs, sporos was {V["fsmem2_mb"]} times faster than minibwa at f = 2 and '
 f'{V["fsmem5_mb"]} times at f = 5, and {V["fsmem2_rb"]} and {V["fsmem5_rb"]} times faster than ropebwt3. Requiring '
 f'two occurrences slowed the search {V["fsmem_cost"]} times relative to SMEMs, since a match that occurs once no '
 f'longer ends it. sporos and `bwa fastmap` differed in {V["bwa_diff_total"]} of {V["bwa_smems_m"]} million SMEMs of '
 f'the short reads; all {V["bwa_only_bwa"]} reported only by BWA pass through an N of the reference, which BWA '
 'replaces with random bases (Supplementary Section S3). Ten threads gave the same SMEMs as one.'),
('table', 'smems'),
('figure', 'speedups'),

('sec', 'Discussion'),
('p',
 'sporos computes the seeds used by read mappers, genome aligners and pangenome tools, following the established '
 'definitions, and on a single core it was faster than each tool it can replace. Its frequency-capped MEMs show that '
 'how often each match\'s string occurs can be derived from the set of matches itself: counting containing intervals '
 'replaces an index of the genomes and reduces capped matching between two chromosomes to about a second.'),
('p',
 'A single library also avoids converting between tools. A pipeline loads a reference once and requests MEMs, MUMs '
 'or SMEMs from Mojo, C or Python, with one coordinate system and one set of occurrence limits, and the command-line '
 'tool writes the MUMmer and `bwa fastmap` formats that existing scripts expect.'),
('p',
 '**Limitations.** All measurements were taken on one x86-64 processor; the 8-lane rows of [#tab:seeds] indicate '
 'performance at AVX2 width, but other processors, including ARM, were not benchmarked. Mumemto was compared on '
 'pairs of genomes rather than the large collections it targets, and its definitions differ from those of sporos. '
 'Positions are 32-bit, limiting each sequence and the concatenated reference to fewer than \\(2^{32}\\) bases, and '
 f'FM-index construction needs about ten bytes of memory per reference base. {V["seed_t10_limit"]}.'),
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
"overview": "The sporos library. Input is parsed and encoded once; three engines compute sampled k-mers and the seed "
            "index that looks them up (top), MEMs, MUMs and frequency-capped MEMs between genomes (middle), and SMEMs "
            "and f-SMEMs of reads (bottom), behind Mojo, command-line, C and Python interfaces.",
"speedups": "Speedup of sporos over the fastest other tool for each task, on one thread (filled) and ten threads "
            "(hollow), on a log scale; the dashed line marks parity and the values are listed on the right (Tables 1 "
            "to 5). Each tool uses its best setting for the thread count: copMEM2 and MUMmer4 gain nothing from threads "
            "on a single query, the SMEM searches of minibwa and BWA have no threads option, Mumemto ran on one "
            "thread, and the seed-lookup harnesses run on one.",
}
