"""Single-thread throughput of `sporos.yak`'s compressed k-mer hashing and hifiasm sketch on random DNA.

`pixi run bench-yak [length]` prints the minimum per-base CPU time of several runs and a checksum of each output,
so that builds for different targets (`--target-cpu`, `-D SPOROS_HASH_LANES=n`) can be compared for speed and
for identical answers.
"""

from std.sys import argv

from timing import thread_cpu_ns
from sporos.yak import FilterTable, Sketch, hpc_kmer_hashes, sketch_read
from sporos.yak.hash import HASH_LANES
from reference import random_dna

comptime REPS = 5


def main() raises:
    var n = 20_000_000
    var args = argv()
    if len(args) > 1:
        n = Int(args[1])
    var seq = random_dna(n, 7, 0.0001)
    print(t"HASH_LANES={HASH_LANES} n={n}")
    var hashes = List[UInt64]()
    var best = Int.MAX
    for _ in range(REPS):
        hashes.clear()
        var t0 = thread_cpu_ns()
        hpc_kmer_hashes(Span(seq), 51, hashes)
        best = min(best, thread_cpu_ns() - t0)
    var check: UInt64 = 0
    for h in hashes:
        check = check * 31 + h
    print(t"hpc_kmer_hashes k=51: {len(hashes)} k-mers, ns/bp={Float64(best) / Float64(n)}, checksum {check}")
    var none = FilterTable()
    var sketch = Sketch()
    best = Int.MAX
    for _ in range(REPS):
        var t0 = thread_cpu_ns()
        sketch_read(Span(seq), 51, 51, True, none, 0, 1000, sketch)
        best = min(best, thread_cpu_ns() - t0)
    check = 0
    for m in sketch.mz:
        check = check * 31 + m.x + UInt64(m.pos)
    print(t"sketch_read k=w=51: {len(sketch.mz)} minimizers, ns/bp={Float64(best) / Float64(n)}, checksum {check}")
