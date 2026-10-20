"""Seed-index build and lookup throughput (GraphAligner-style seeding).

    pixi run bench-index <reference.fa> <reads.fa> [k] [w] [max_occurrences]

Builds a (w, k)-minimizer index of the reference, then looks up every read's
k-mers (GraphAligner's semantics) and, separately, only its minimizers, on
one thread. Prints build time and lookup throughput in Mbp/s.
"""

from std.sys import argv

from std.time import perf_counter_ns

from timing import thread_cpu_ns
from sporos.fasta import read_fasta
from sporos.index import MINIMIZER_SEEDS, SeedHit, SeedIndex


def main() raises:
    var args = argv()
    if len(args) < 3:
        print("usage: bench_index <reference.fa> <reads.fa> [k] [w] [max_occurrences]")
        return
    var k = Int(args[3]) if len(args) > 3 else 19
    var w = Int(args[4]) if len(args) > 4 else 12
    var max_occ = Int(args[5]) if len(args) > 5 else 0
    var reference = read_fasta(String(args[1]))
    var reads = read_fasta(String(args[2]))
    var seqs = List[Span[UInt8, origin_of(reference.data)]]()
    for r in range(len(reference)):
        seqs.append(reference.seq(r))
    for threads in [1, 8]:
        var start = thread_cpu_ns()
        var wall = perf_counter_ns()
        var index = SeedIndex(seqs, k, w, MINIMIZER_SEEDS, threads)
        print(t"build threads={threads}: {Float64(perf_counter_ns() - wall) / 1e9} s wall, {len(index.table)} seeds")
        _ = start
    var index = SeedIndex(seqs, k, w, MINIMIZER_SEEDS, 1)
    var bases = reads.total_length()
    var spans = List[Span[UInt8, origin_of(reads.data)]](capacity=len(reads))
    for r in range(len(reads)):
        spans.append(reads.seq(r))
    for all_kmers in [True, False]:
        var mode = "all k-mers" if all_kmers else "minimizers"
        var best = Int.MAX
        var count = 0
        for _ in range(3):
            var hits = List[SeedHit]()
            count = 0
            var start = thread_cpu_ns()
            for r in range(len(reads)):
                hits.clear()
                index.lookup(reads.seq(r), hits, all_kmers, max_occ)
                count += len(hits)
            best = min(best, thread_cpu_ns() - start)
        var mbps = Float64(bases) / (Float64(best) / 1e9) / 1e6
        print(t"lookup {mode}, per read: {mbps} Mbp/s, {Float64(count) / Float64(len(reads))} hits/read")
        best = Int.MAX
        var batch_count = 0
        for _ in range(3):
            var start = thread_cpu_ns()
            var batch = index.lookup_batch(spans, all_kmers, max_occ)
            best = min(best, thread_cpu_ns() - start)
            batch_count = len(batch.hits)
        mbps = Float64(bases) / (Float64(best) / 1e9) / 1e6
        print(t"lookup {mode}, batch: {mbps} Mbp/s, {Float64(batch_count) / Float64(len(reads))} hits/read")
        if batch_count != count:
            raise Error("batch and per-read lookups disagree")
        for threads in [8]:
            var wall = perf_counter_ns()
            var batch = index.lookup_batch(spans, all_kmers, max_occ, threads)
            var seconds = Float64(perf_counter_ns() - wall) / 1e9
            print(t"lookup {mode}, batch, {threads} threads: {Float64(bases) / seconds / 1e6} Mbp/s wall")
            _ = len(batch)
