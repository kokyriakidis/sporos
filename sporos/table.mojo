"""A parallel-built hash table of 64-bit entries in CSR layout.

Anchors come from an `AnchorSource` as (hash, entry) pairs, chunk by chunk.
Bucket `b` holds the entries whose hash has `b` as its top bits, in
`entries[offsets[b] : offsets[b + 1]]`. Construction is a deterministic
two-level radix sort: chunks count and fill their anchors in parallel, scatter
them to 256 parts by the bucket's top bits, then each part is counting-sorted
by bucket. Entries keep their source order within a bucket unless
`sort_buckets` asks for them sorted by value.
"""

from std.math import ceildiv
from std.memory import Layout, alloc, dealloc

from sporos.parallel import ParallelTask, parallel_for


trait AnchorSource(Deinitable):
    """Produces a table's anchors in independent chunks."""

    def chunks(self) -> Int:
        ...

    def count(self, chunk: Int) -> Int:
        """Anchors chunk `chunk` will write."""
        ...

    def fill(
        self,
        chunk: Int,
        hashes: MutPointer[UInt64, MutUntrackedOrigin],
        entries: MutPointer[UInt64, MutUntrackedOrigin],
    ):
        """Writes the chunk's `count(chunk)` anchors."""
        ...


struct HashTable(Movable, Sized):
    var offsets: List[UInt32]
    var entries: List[UInt64]
    var shift: UInt64
    """Bucket of a hash: `hash >> shift`."""

    @inline(.always)
    def bucket(self, hash: UInt64) -> Int:
        return Int(hash >> self.shift)

    def __len__(self) -> Int:
        return len(self.entries)

    def __init__[S: AnchorSource](out self, source: S, threads: Int, sort_buckets: Bool = False):
        var chunks = source.chunks()
        var firsts_alloc = alloc(Layout[Int](length=chunks + 1))
        var firsts = firsts_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        var src = Pointer(to=source).unsafe_origin_cast[ImmUntrackedOrigin]()
        var count_task = _CountTask[S](src, firsts)
        parallel_for(count_task, chunks, threads)
        firsts[unsafe_offset=0] = 0
        for c in range(chunks):
            firsts[unsafe_offset=c + 1] += firsts[unsafe_offset=c]
        var total = firsts[unsafe_offset=chunks]

        var bits = 1
        while (1 << bits) < max(1, total):
            bits += 1
        self.shift = UInt64(64 - bits)
        var num_buckets = 1 << bits
        var part_bits = min(8, bits)
        var parts = 1 << part_bits
        var part_shift = UInt64(64 - part_bits)

        # Fill each chunk's anchors and count them per part.
        var hashes_alloc = alloc(Layout[UInt64](length=max(1, total)))
        var values_alloc = alloc(Layout[UInt64](length=max(1, total)))
        var histogram_alloc = alloc(Layout[Int](length=chunks * parts))
        var hashes = hashes_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        var values = values_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        var histogram = histogram_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        var fill_task = _FillTask[S](src, firsts, hashes, values, histogram, parts, part_shift)
        parallel_for(fill_task, chunks, threads)

        # Part-major start of every (chunk, part) cell, in place.
        var starts_alloc = alloc(Layout[Int](length=parts + 1))
        var part_starts = starts_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        var at = 0
        for p in range(parts):
            part_starts[unsafe_offset=p] = at
            for c in range(chunks):
                var cell = c * parts + p
                var count = histogram[unsafe_offset=cell]
                histogram[unsafe_offset=cell] = at
                at += count
        part_starts[unsafe_offset=parts] = at
        var part_hashes_alloc = alloc(Layout[UInt64](length=max(1, total)))
        var part_values_alloc = alloc(Layout[UInt64](length=max(1, total)))
        var part_hashes = part_hashes_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        var part_values = part_values_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        var scatter_task = _ScatterTask(firsts, hashes, values, histogram, part_hashes, part_values, parts, part_shift)
        parallel_for(scatter_task, chunks, threads)
        dealloc(hashes_alloc^)
        dealloc(values_alloc^)
        dealloc(histogram_alloc^)

        self.offsets = List[UInt32](unsafe_uninit_length=num_buckets + 1)
        self.entries = List[UInt64](unsafe_uninit_length=total)
        var sort_task = _PartSortTask(
            part_starts,
            part_hashes,
            part_values,
            self.offsets.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
            self.entries.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
            num_buckets // parts,
            self.shift,
            sort_buckets,
        )
        parallel_for(sort_task, parts, threads)
        self.offsets[num_buckets] = UInt32(total)
        dealloc(firsts_alloc^)
        dealloc(starts_alloc^)
        dealloc(part_hashes_alloc^)
        dealloc(part_values_alloc^)


struct _CountTask[S: AnchorSource](ParallelTask):
    var source: ImmPointer[Self.S, ImmUntrackedOrigin]
    var firsts: MutPointer[Int, MutUntrackedOrigin]

    def __init__(out self, source: ImmPointer[Self.S, ImmUntrackedOrigin], firsts: MutPointer[Int, MutUntrackedOrigin]):
        self.source = source
        self.firsts = firsts

    def run(self, index: Int):
        self.firsts[unsafe_offset=index + 1] = self.source[].count(index)


struct _FillTask[S: AnchorSource](ParallelTask):
    var source: ImmPointer[Self.S, ImmUntrackedOrigin]
    var firsts: ImmPointer[Int, ImmUntrackedOrigin]
    var hashes: MutPointer[UInt64, MutUntrackedOrigin]
    var values: MutPointer[UInt64, MutUntrackedOrigin]
    var histogram: MutPointer[Int, MutUntrackedOrigin]
    var parts: Int
    var part_shift: UInt64

    def __init__(
        out self,
        source: ImmPointer[Self.S, ImmUntrackedOrigin],
        firsts: ImmPointer[Int, ImmUntrackedOrigin],
        hashes: MutPointer[UInt64, MutUntrackedOrigin],
        values: MutPointer[UInt64, MutUntrackedOrigin],
        histogram: MutPointer[Int, MutUntrackedOrigin],
        parts: Int,
        part_shift: UInt64,
    ):
        self.source = source
        self.firsts = firsts
        self.hashes = hashes
        self.values = values
        self.histogram = histogram
        self.parts = parts
        self.part_shift = part_shift

    def run(self, index: Int):
        var lo = self.firsts[unsafe_offset=index]
        var hi = self.firsts[unsafe_offset=index + 1]
        self.source[].fill(index, self.hashes.unsafe_offset(lo), self.values.unsafe_offset(lo))
        var histogram = self.histogram.unsafe_offset(index * self.parts)
        for p in range(self.parts):
            histogram[unsafe_offset=p] = 0
        for i in range(lo, hi):
            histogram[unsafe_offset=Int(self.hashes[unsafe_offset=i] >> self.part_shift)] += 1


struct _ScatterTask(ParallelTask):
    """Moves one chunk's anchors to their parts, keeping order within parts;
    `cursors` holds the chunk's start in each part."""

    var firsts: ImmPointer[Int, ImmUntrackedOrigin]
    var hashes: ImmPointer[UInt64, ImmUntrackedOrigin]
    var values: ImmPointer[UInt64, ImmUntrackedOrigin]
    var cursors: MutPointer[Int, MutUntrackedOrigin]
    var part_hashes: MutPointer[UInt64, MutUntrackedOrigin]
    var part_values: MutPointer[UInt64, MutUntrackedOrigin]
    var parts: Int
    var part_shift: UInt64

    def __init__(
        out self,
        firsts: ImmPointer[Int, ImmUntrackedOrigin],
        hashes: ImmPointer[UInt64, ImmUntrackedOrigin],
        values: ImmPointer[UInt64, ImmUntrackedOrigin],
        cursors: MutPointer[Int, MutUntrackedOrigin],
        part_hashes: MutPointer[UInt64, MutUntrackedOrigin],
        part_values: MutPointer[UInt64, MutUntrackedOrigin],
        parts: Int,
        part_shift: UInt64,
    ):
        self.firsts = firsts
        self.hashes = hashes
        self.values = values
        self.cursors = cursors
        self.part_hashes = part_hashes
        self.part_values = part_values
        self.parts = parts
        self.part_shift = part_shift

    def run(self, index: Int):
        var cursors = self.cursors.unsafe_offset(index * self.parts)
        for i in range(self.firsts[unsafe_offset=index], self.firsts[unsafe_offset=index + 1]):
            var h = self.hashes[unsafe_offset=i]
            var p = Int(h >> self.part_shift)
            var slot = cursors[unsafe_offset=p]
            self.part_hashes[unsafe_offset=slot] = h
            self.part_values[unsafe_offset=slot] = self.values[unsafe_offset=i]
            cursors[unsafe_offset=p] = slot + 1


struct _PartSortTask(ParallelTask):
    """Counting sort of one part's anchors by bucket into the final table."""

    var part_starts: ImmPointer[Int, ImmUntrackedOrigin]
    var part_hashes: ImmPointer[UInt64, ImmUntrackedOrigin]
    var part_values: ImmPointer[UInt64, ImmUntrackedOrigin]
    var offsets: MutPointer[UInt32, MutUntrackedOrigin]
    var entries: MutPointer[UInt64, MutUntrackedOrigin]
    var buckets_per_part: Int
    var shift: UInt64
    var sort_buckets: Bool

    def __init__(
        out self,
        part_starts: ImmPointer[Int, ImmUntrackedOrigin],
        part_hashes: ImmPointer[UInt64, ImmUntrackedOrigin],
        part_values: ImmPointer[UInt64, ImmUntrackedOrigin],
        offsets: MutPointer[UInt32, MutUntrackedOrigin],
        entries: MutPointer[UInt64, MutUntrackedOrigin],
        buckets_per_part: Int,
        shift: UInt64,
        sort_buckets: Bool,
    ):
        self.part_starts = part_starts
        self.part_hashes = part_hashes
        self.part_values = part_values
        self.offsets = offsets
        self.entries = entries
        self.buckets_per_part = buckets_per_part
        self.shift = shift
        self.sort_buckets = sort_buckets

    def run(self, index: Int):
        var lo = self.part_starts[unsafe_offset=index]
        var hi = self.part_starts[unsafe_offset=index + 1]
        var first_bucket = index * self.buckets_per_part
        var offsets = self.offsets.unsafe_offset(first_bucket)
        for b in range(self.buckets_per_part):
            offsets[unsafe_offset=b] = 0
        for i in range(lo, hi):
            offsets[unsafe_offset=Int(self.part_hashes[unsafe_offset=i] >> self.shift) - first_bucket] += 1
        var at = UInt32(lo)
        for b in range(self.buckets_per_part):
            var count = offsets[unsafe_offset=b]
            offsets[unsafe_offset=b] = at
            at += count
        # offsets[b] advances to the end of bucket b while scattering, which
        # is the start of bucket b + 1; shift back afterwards.
        for i in range(lo, hi):
            var b = Int(self.part_hashes[unsafe_offset=i] >> self.shift) - first_bucket
            var slot = offsets[unsafe_offset=b]
            self.entries[unsafe_offset=Int(slot)] = self.part_values[unsafe_offset=i]
            offsets[unsafe_offset=b] = slot + 1
        for b in reversed(range(1, self.buckets_per_part)):
            offsets[unsafe_offset=b] = offsets[unsafe_offset=b - 1]
        offsets[unsafe_offset=0] = UInt32(lo)
        if self.sort_buckets:
            # Buckets are small: insertion sort each by entry value.
            for b in range(self.buckets_per_part):
                var start = Int(offsets[unsafe_offset=b])
                var end = Int(offsets[unsafe_offset=b + 1]) if b + 1 < self.buckets_per_part else hi
                for i in range(start + 1, end):
                    var v = self.entries[unsafe_offset=i]
                    var j = i
                    while j > start and self.entries[unsafe_offset=j - 1] > v:
                        self.entries[unsafe_offset=j] = self.entries[unsafe_offset=j - 1]
                        j -= 1
                    self.entries[unsafe_offset=j] = v
