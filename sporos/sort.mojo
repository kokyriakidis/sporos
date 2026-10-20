"""Parallel sorting of 64-bit keys.

Keys are partitioned by their top bits into buckets (a counting pass and a
scatter, both split across threads), then every bucket is sorted on its own.
"""

from std.bit import count_leading_zeros
from std.math import ceildiv
from std.memory import Layout, alloc, dealloc

from sporos.parallel import ParallelTask, parallel_for

comptime _BUCKET_BITS = 12
comptime _SERIAL = 1 << 16
"""Below this many keys, one `sort` call."""


def parallel_sort(mut keys: List[UInt64], threads: Int):
    """Sorts `keys` in increasing order using up to `threads` threads."""
    var n = len(keys)
    if threads <= 1 or n < _SERIAL:
        sort(Span(keys))
        return
    var top = UInt64(0)
    for k in keys:
        top = max(top, k)
    var bits = 64 - Int(count_leading_zeros(top | 1))
    var shift = UInt64(max(0, bits - _BUCKET_BITS))
    var buckets = 1 << min(bits, _BUCKET_BITS)
    var chunks = threads * 4
    var histogram_alloc = alloc(Layout[Int](length=chunks * buckets))
    var histogram = histogram_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
    var src = keys.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
    var count_task = _CountTask(src, n, chunks, histogram, buckets, shift)
    parallel_for(count_task, chunks, threads)
    # Bucket-major exclusive prefix: each chunk's slot in every bucket.
    var starts_alloc = alloc(Layout[Int](length=buckets + 1))
    var starts = starts_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
    var at = 0
    for b in range(buckets):
        starts[unsafe_offset=b] = at
        for c in range(chunks):
            var cell = c * buckets + b
            var count = histogram[unsafe_offset=cell]
            histogram[unsafe_offset=cell] = at
            at += count
    starts[unsafe_offset=buckets] = at
    var sorted_alloc = alloc(Layout[UInt64](length=n))
    var dst = sorted_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
    var scatter_task = _ScatterTask(src, dst, n, chunks, histogram, buckets, shift)
    parallel_for(scatter_task, chunks, threads)
    var sort_task = _BucketSortTask(dst, starts)
    parallel_for(sort_task, buckets, threads)
    for i in range(n):
        src[unsafe_offset=i] = dst[unsafe_offset=i]
    dealloc(histogram_alloc^)
    dealloc(starts_alloc^)
    dealloc(sorted_alloc^)


struct _CountTask(ParallelTask):
    var keys: MutPointer[UInt64, MutUntrackedOrigin]
    var n: Int
    var chunks: Int
    var histogram: MutPointer[Int, MutUntrackedOrigin]
    var buckets: Int
    var shift: UInt64

    def __init__(
        out self,
        keys: MutPointer[UInt64, MutUntrackedOrigin],
        n: Int,
        chunks: Int,
        histogram: MutPointer[Int, MutUntrackedOrigin],
        buckets: Int,
        shift: UInt64,
    ):
        self.keys = keys
        self.n = n
        self.chunks = chunks
        self.histogram = histogram
        self.buckets = buckets
        self.shift = shift

    def run(self, index: Int):
        var per = ceildiv(self.n, self.chunks)
        var histogram = self.histogram.unsafe_offset(index * self.buckets)
        for b in range(self.buckets):
            histogram[unsafe_offset=b] = 0
        for i in range(index * per, min(self.n, (index + 1) * per)):
            histogram[unsafe_offset=Int(self.keys[unsafe_offset=i] >> self.shift)] += 1


struct _ScatterTask(ParallelTask):
    var keys: MutPointer[UInt64, MutUntrackedOrigin]
    var out: MutPointer[UInt64, MutUntrackedOrigin]
    var n: Int
    var chunks: Int
    var histogram: MutPointer[Int, MutUntrackedOrigin]
    var buckets: Int
    var shift: UInt64

    def __init__(
        out self,
        keys: MutPointer[UInt64, MutUntrackedOrigin],
        out_keys: MutPointer[UInt64, MutUntrackedOrigin],
        n: Int,
        chunks: Int,
        histogram: MutPointer[Int, MutUntrackedOrigin],
        buckets: Int,
        shift: UInt64,
    ):
        self.keys = keys
        self.out = out_keys
        self.n = n
        self.chunks = chunks
        self.histogram = histogram
        self.buckets = buckets
        self.shift = shift

    def run(self, index: Int):
        var per = ceildiv(self.n, self.chunks)
        var slots = self.histogram.unsafe_offset(index * self.buckets)
        for i in range(index * per, min(self.n, (index + 1) * per)):
            var k = self.keys[unsafe_offset=i]
            var b = Int(k >> self.shift)
            self.out[unsafe_offset=slots[unsafe_offset=b]] = k
            slots[unsafe_offset=b] += 1


struct _BucketSortTask(ParallelTask):
    var keys: MutPointer[UInt64, MutUntrackedOrigin]
    var starts: MutPointer[Int, MutUntrackedOrigin]

    def __init__(out self, keys: MutPointer[UInt64, MutUntrackedOrigin], starts: MutPointer[Int, MutUntrackedOrigin]):
        self.keys = keys
        self.starts = starts

    def run(self, index: Int):
        var lo = self.starts[unsafe_offset=index]
        var hi = self.starts[unsafe_offset=index + 1]
        if hi - lo > 1:
            sort(Span(unsafe_ptr=self.keys.unsafe_offset(lo), length=hi - lo))
