"""Minimizer/syncmer throughput on a real FASTA record, single- or multi-threaded.

`bench_seeds_fasta <seqs.fa> [threads]` runs the configurations of
`benchmarks/bench_seeds.mojo` on the first record. One thread reports the
minimum thread CPU time; more threads report the minimum wall time, as the
simd-minimizers and minimap2 harnesses do.
"""

from std.sys import argv
from std.time import perf_counter_ns

from timing import thread_cpu_ns
from sporos.alphabet import CodeBuffer
from sporos.fasta import read_fasta
from sporos.minimizers import minimizer_positions_from_codes
from sporos.syncmers import syncmer_positions_from_codes
from sporos.window import LANES, PADDING

comptime REPS = 5


def _now(threads: Int) -> Int:
    return thread_cpu_ns() if threads == 1 else perf_counter_ns()


def _report(name: String, k: Int, w: Int, n: Int, threads: Int, count: Int, best_ns: Int):
    var ns_per_bp = Float64(best_ns) / Float64(n)
    print(t"sporos {name} k={k} w={w} T={threads} lanes={LANES} n={n} out={count} ns/bp={ns_per_bp}")


def main() raises:
    var args = argv()
    if len(args) < 2:
        raise Error("usage: bench_seeds_fasta <seqs.fa> [threads]")
    var threads = Int(args[2]) if len(args) > 2 else 1
    var records = read_fasta(String(args[1]))
    var seq = records.seq(0)
    var n = len(seq)
    var start = perf_counter_ns()
    var codes = CodeBuffer.encode(seq, PADDING)
    print(t"sporos encode n={n} ns/bp={Float64(perf_counter_ns() - start) / Float64(n)}")
    for cfg in [(21, 11), (31, 19)]:
        var k = cfg[0]
        var w = cfg[1]
        for canonical in [False, True]:
            var best = Int.MAX
            var count = 0
            for _ in range(REPS):
                var t0 = _now(threads)
                var m = minimizer_positions_from_codes(codes.unsafe_ptr(), n, k, w, canonical, threads)
                best = min(best, _now(threads) - t0)
                count = len(m)
            _report("can" if canonical else "fwd", k, w, n, threads, count, best)
        for canonical in [False, True]:
            # simd-minimizers' closed syncmers take (k=s, w=k-s+1).
            var best = Int.MAX
            var count = 0
            for _ in range(REPS):
                var t0 = _now(threads)
                var m = syncmer_positions_from_codes(codes.unsafe_ptr(), n, k + w - 1, k, canonical, None, threads)
                best = min(best, _now(threads) - t0)
                count = len(m)
            _report("closed-can" if canonical else "closed-fwd", k, w, n, threads, count, best)
