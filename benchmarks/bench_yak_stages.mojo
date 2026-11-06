"""Every stage of `sporos.yak` on a read file, as hifiasm's candidate detection runs them: counting (exact and
through Bloom filters), peaks, filter, sketches, position table. Prints each stage's wall and CPU time.

`pixi run bench-yak-stages <reads.fq[.gz]> [threads] [bloom 0/1]`.
"""

from std.sys import argv
from std.time import perf_counter_ns

from timing import thread_cpu_ns
from sporos import read_fasta
from sporos.yak import FilterTable, KmerCounts, PositionTable, analyze_count, sketch_all


def _report(name: String, wall_ns: Int, extra: String = ""):
    print(t"{name}\t{Float64(wall_ns) / 1e9} s wall\t{extra}")


def main() raises:
    var args = argv()
    var path = String(args[1])
    var threads = Int(args[2]) if len(args) > 2 else 1
    var bloom = Int(args[3]) if len(args) > 3 else 0
    var t0 = perf_counter_ns()
    var records = read_fasta(path, threads)
    var t1 = perf_counter_ns()
    _report("read", t1 - t0, String(t"{len(records)} reads {records.total_length()} bp"))
    var counts = KmerCounts(records, 51, threads, 1 << 24, 37 if bloom == 1 else 0)
    var t2 = perf_counter_ns()
    _report("count", t2 - t1, String(t"{counts.distinct()} distinct"))
    var peaks = analyze_count(counts.histogram(threads), 5)
    var t3 = perf_counter_ns()
    _report("peaks", t3 - t2, String(t"hom {peaks.hom} het {peaks.het}"))
    var filter = FilterTable(counts, min(Int(Float64(peaks.hom) * 5.0), 4094), 2000, threads)
    var t4 = perf_counter_ns()
    _report("filter", t4 - t3, String(t"{filter.size} k-mers"))
    var sketches = sketch_all(records, filter, 51, 51, threads, 500, 1000)
    var t5 = perf_counter_ns()
    var minimizers = 0
    for s in range(len(sketches)):
        minimizers += len(sketches[s].mz)
    _report("sketch", t5 - t4, String(t"{minimizers} minimizers"))
    var table = PositionTable(sketches, threads)
    var t6 = perf_counter_ns()
    _report("table", t6 - t5, String(t"{table.indexed} positions of {table.distinct} hashes"))
    _report("total", t6 - t0)
