"""Per-call overhead of the seed kernels on read-sized sequences.

    pixi run bench-reads

Times minimizers (k21/w11, canonical) over 20 Mbp of random DNA cut into
reads of several lengths, one call per read, against one call over all of it.
"""

from timing import thread_cpu_ns
from sporos.batch import minimizers_batch
from sporos.minimizers import minimizer_positions
from reference import random_dna

comptime TOTAL = 20_000_000
comptime REPS = 5


def main() raises:
    var seq = random_dna(TOTAL, 11)
    var whole = Int.MAX
    for _ in range(REPS):
        var start = thread_cpu_ns()
        var m = minimizer_positions(Span(seq), 21, 11, True, 1)
        whole = min(whole, thread_cpu_ns() - start)
        _ = len(m)
    print(t"whole sequence: {Float64(whole) / Float64(TOTAL)} ns/bp")
    for read_len in [150, 1000, 10_000, 100_000]:
        var best = Int.MAX
        var reads = TOTAL // read_len
        for _ in range(REPS):
            var start = thread_cpu_ns()
            var count = 0
            for r in range(reads):
                var m = minimizer_positions(Span(seq)[r * read_len : (r + 1) * read_len], 21, 11, True, 1)
                count += len(m)
            best = min(best, thread_cpu_ns() - start)
        print(t"reads of {read_len}: {Float64(best) / Float64(reads * read_len)} ns/bp (one call per read)")
        var spans = List[Span[UInt8, origin_of(seq)]]()
        for r in range(reads):
            spans.append(Span(seq)[r * read_len : (r + 1) * read_len])
        best = Int.MAX
        for _ in range(REPS):
            var start = thread_cpu_ns()
            var batch = minimizers_batch(spans, 21, 11, True, 1)
            best = min(best, thread_cpu_ns() - start)
            _ = len(batch)
        print(t"reads of {read_len}: {Float64(best) / Float64(reads * read_len)} ns/bp (batch)")
