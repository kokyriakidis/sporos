"""Single-thread throughput of minimizers and syncmers on random DNA.

`pixi run bench-seeds [length]` prints the minimum per-base CPU time of
several runs for the configurations simd-minimizers reports.
"""

from std.sys import argv

from timing import thread_cpu_ns
from sporos.alphabet import CodeBuffer
from sporos.minimizers import minimizer_positions_from_codes
from sporos.syncmers import syncmer_positions_from_codes
from sporos.window import PADDING
from reference import random_dna

comptime REPS = 7


def _report(name: String, k: Int, w: Int, n: Int, count: Int, best_ns: Int):
    var density = Float64(count) / Float64(n)
    var ns_per_bp = Float64(best_ns) / Float64(n)
    print(name, t"k={k} w={w} n={n} out={count} density={density} ns/bp={ns_per_bp}")


def main() raises:
    var n = 100_000_000
    var args = argv()
    if len(args) > 1:
        n = Int(args[1])
    var seq = random_dna(n, 5)
    var codes = CodeBuffer.encode(Span(seq), PADDING)
    for cfg in [(21, 11), (31, 19)]:
        var k = cfg[0]
        var w = cfg[1]
        for canonical in [False, True]:
            var best = Int.MAX
            var count = 0
            for _ in range(REPS):
                var start = thread_cpu_ns()
                var m = minimizer_positions_from_codes(codes.unsafe_ptr(), n, k, w, canonical, 1)
                best = min(best, thread_cpu_ns() - start)
                count = len(m)
            _report("can" if canonical else "fwd", k, w, n, count, best)
        for canonical in [False, True]:
            # simd-minimizers' closed syncmers take (k=s, w=k-s+1).
            var best = Int.MAX
            var count = 0
            for _ in range(REPS):
                var start = thread_cpu_ns()
                var m = syncmer_positions_from_codes(codes.unsafe_ptr(), n, k + w - 1, k, canonical, None, 1)
                best = min(best, thread_cpu_ns() - start)
                count = len(m)
            _report("closed-can" if canonical else "closed-fwd", k, w, n, count, best)
