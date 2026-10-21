"""A parallel-built hash table of 64-bit entries in CSR layout.

Anchors come from an `AnchorSource` as (hash, entry) pairs, chunk by chunk.
Bucket `b` holds the entries whose hash has `b` as its top bits, in
`entries[offsets[b] : offsets[b + 1]]`. Construction is a deterministic
two-level radix sort: chunks count and fill their anchors in parallel, scatter
them to 1024 parts by the bucket's top bits, then each part is counting-sorted
by bucket in place. Entries keep their source order within a bucket unless
`sort_buckets` asks for them sorted by value.
"""

from std.math import ceildiv
from std.memory import Layout, alloc, dealloc

from sporos.pages import advise_huge_pages
from sporos.parallel import ParallelTask, parallel_for


trait AnchorSource(Deinitable):
    """Produces a table's anchors in independent chunks."""

    def chunks(self) -> Int:
        ...

    def capacity(self, chunk: Int) -> Int:
        """At least as many as the anchors chunk `chunk` writes."""
        ...

    def fill(
        self,
        chunk: Int,
        hashes: MutPointer[UInt64, MutUntrackedOrigin],
        entries: MutPointer[UInt64, MutUntrackedOrigin],
    ) -> Int:
        """Writes the chunk's anchors (at most `capacity(chunk)`) and returns
        how many."""
        ...


comptime _PART_BITS = 10
"""First-level radix: anchors are scattered to 2^_PART_BITS parts, small
enough that each part's counting sort stays in L2 (best of 8-12 bits on
hg38 chr1)."""

comptime ENTRY_PADDING = 8
"""Readable (garbage) entries past the end, for whole-vector bucket scans."""


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

    def __init__[S: AnchorSource](out self, source: S, threads: Int, sort_buckets: Bool = False, bucket_load: Int = 1):
        """Hashes `source`'s anchors into buckets of about `bucket_load`
        entries. `entries` has `ENTRY_PADDING` readable slots past its end.

        Chunks are filled twice (to count, then to scatter) rather than kept,
        so the build holds 12 bytes per anchor beyond the table itself."""
        var chunks = source.chunks()
        var src = Pointer(to=source).unsafe_origin_cast[ImmUntrackedOrigin]()
        # Count every chunk's anchors by the top _PART_BITS bits of their hash.
        comptime max_parts = 1 << _PART_BITS
        var histogram_alloc = alloc(Layout[Int](length=chunks * max_parts))
        var histogram = histogram_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        var count_task = _PartCountTask[S](src, histogram, max_parts, UInt64(64 - _PART_BITS))
        parallel_for(count_task, chunks, threads)
        var total = 0
        for cell in range(chunks * max_parts):
            total += histogram[unsafe_offset=cell]

        var bits = 1
        while (1 << bits) * bucket_load < max(1, total):
            bits += 1
        self.shift = UInt64(64 - bits)
        var num_buckets = 1 << bits
        var part_bits = min(_PART_BITS, bits)
        var parts = 1 << part_bits
        var part_shift = UInt64(64 - part_bits)
        if parts < max_parts:
            # A tiny table: fold the counts into its fewer parts, in place.
            var fold = max_parts // parts
            for c in range(chunks):
                for p in range(parts):
                    var sum = 0
                    for q in range(p * fold, (p + 1) * fold):
                        sum += histogram[unsafe_offset=c * max_parts + q]
                    histogram[unsafe_offset=c * parts + p] = sum
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

        # Scatter straight into the table's entries, with each anchor's bucket;
        # sorting each part by bucket in place then finishes the table.
        self.entries = List[UInt64](capacity=total + ENTRY_PADDING)
        self.entries.resize(unsafe_uninit_length=total)
        self.offsets = List[UInt32](unsafe_uninit_length=num_buckets + 1)
        # Probes hit both arrays at random: huge pages spare them a page walk each.
        advise_huge_pages(Int(self.offsets.unsafe_ptr()), (num_buckets + 1) * 4)
        advise_huge_pages(Int(self.entries.unsafe_ptr()), total * 8)
        var entries = self.entries.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        var buckets_alloc = alloc(Layout[UInt32](length=max(1, total)))
        var buckets = buckets_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        var scatter_task = _ScatterTask[S](src, histogram, entries, buckets, parts, part_shift, self.shift)
        parallel_for(scatter_task, chunks, threads)
        dealloc(histogram_alloc^)

        var sort_task = _PartSortTask(
            part_starts,
            buckets,
            self.offsets.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
            entries,
            num_buckets // parts,
            sort_buckets,
        )
        parallel_for(sort_task, parts, threads)
        self.offsets[num_buckets] = UInt32(total)
        dealloc(starts_alloc^)
        dealloc(buckets_alloc^)


struct _ChunkAnchors(Movable):
    """One chunk's anchors, filled into its own buffers."""

    var hashes: List[UInt64]
    var values: List[UInt64]

    def __init__[S: AnchorSource](out self, source: S, chunk: Int):
        var room = source.capacity(chunk)
        self.hashes = List[UInt64](unsafe_uninit_length=room)
        self.values = List[UInt64](unsafe_uninit_length=room)
        var count = source.fill(
            chunk,
            self.hashes.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
            self.values.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        )
        self.hashes.resize(unsafe_uninit_length=count)
        self.values.resize(unsafe_uninit_length=count)


struct _PartCountTask[S: AnchorSource](ParallelTask):
    var source: ImmPointer[Self.S, ImmUntrackedOrigin]
    var histogram: MutPointer[Int, MutUntrackedOrigin]
    var parts: Int
    var part_shift: UInt64

    def __init__(
        out self,
        source: ImmPointer[Self.S, ImmUntrackedOrigin],
        histogram: MutPointer[Int, MutUntrackedOrigin],
        parts: Int,
        part_shift: UInt64,
    ):
        self.source = source
        self.histogram = histogram
        self.parts = parts
        self.part_shift = part_shift

    def run(self, index: Int):
        var anchors = _ChunkAnchors(self.source[], index)
        var histogram = self.histogram.unsafe_offset(index * self.parts)
        for p in range(self.parts):
            histogram[unsafe_offset=p] = 0
        for h in anchors.hashes:
            histogram[unsafe_offset=Int(h >> self.part_shift)] += 1


struct _ScatterTask[S: AnchorSource](ParallelTask):
    """Refills one chunk and moves its anchors to their parts, keeping order
    within parts; `cursors` holds the chunk's start in each part."""

    var source: ImmPointer[Self.S, ImmUntrackedOrigin]
    var cursors: MutPointer[Int, MutUntrackedOrigin]
    var entries: MutPointer[UInt64, MutUntrackedOrigin]
    var buckets: MutPointer[UInt32, MutUntrackedOrigin]
    var parts: Int
    var part_shift: UInt64
    var shift: UInt64

    def __init__(
        out self,
        source: ImmPointer[Self.S, ImmUntrackedOrigin],
        cursors: MutPointer[Int, MutUntrackedOrigin],
        entries: MutPointer[UInt64, MutUntrackedOrigin],
        buckets: MutPointer[UInt32, MutUntrackedOrigin],
        parts: Int,
        part_shift: UInt64,
        shift: UInt64,
    ):
        self.source = source
        self.cursors = cursors
        self.entries = entries
        self.buckets = buckets
        self.parts = parts
        self.part_shift = part_shift
        self.shift = shift

    def run(self, index: Int):
        var anchors = _ChunkAnchors(self.source[], index)
        var cursors = self.cursors.unsafe_offset(index * self.parts)
        for i in range(len(anchors.hashes)):
            var h = anchors.hashes[i]
            var p = Int(h >> self.part_shift)
            var slot = cursors[unsafe_offset=p]
            self.entries[unsafe_offset=slot] = anchors.values[i]
            self.buckets[unsafe_offset=slot] = UInt32(h >> self.shift)
            cursors[unsafe_offset=p] = slot + 1


struct _PartSortTask(ParallelTask):
    """Counting sort of one part's entries by bucket, in place."""

    var part_starts: ImmPointer[Int, ImmUntrackedOrigin]
    var buckets: ImmPointer[UInt32, ImmUntrackedOrigin]
    var offsets: MutPointer[UInt32, MutUntrackedOrigin]
    var entries: MutPointer[UInt64, MutUntrackedOrigin]
    var buckets_per_part: Int
    var sort_buckets: Bool

    def __init__(
        out self,
        part_starts: ImmPointer[Int, ImmUntrackedOrigin],
        buckets: ImmPointer[UInt32, ImmUntrackedOrigin],
        offsets: MutPointer[UInt32, MutUntrackedOrigin],
        entries: MutPointer[UInt64, MutUntrackedOrigin],
        buckets_per_part: Int,
        sort_buckets: Bool,
    ):
        self.part_starts = part_starts
        self.buckets = buckets
        self.offsets = offsets
        self.entries = entries
        self.buckets_per_part = buckets_per_part
        self.sort_buckets = sort_buckets

    def run(self, index: Int):
        var lo = self.part_starts[unsafe_offset=index]
        var hi = self.part_starts[unsafe_offset=index + 1]
        var first_bucket = index * self.buckets_per_part
        var offsets = self.offsets.unsafe_offset(first_bucket)
        for b in range(self.buckets_per_part):
            offsets[unsafe_offset=b] = 0
        for i in range(lo, hi):
            offsets[unsafe_offset=Int(self.buckets[unsafe_offset=i]) - first_bucket] += 1
        var at = UInt32(lo)
        for b in range(self.buckets_per_part):
            var count = offsets[unsafe_offset=b]
            offsets[unsafe_offset=b] = at
            at += count
        # The part's entries, by bucket, through a copy; offsets[b] advances to
        # the end of bucket b while scattering, which is the start of bucket
        # b + 1, so shift back afterwards.
        var part = List[UInt64](capacity=hi - lo)
        for i in range(lo, hi):
            part.append(self.entries[unsafe_offset=i])
        for i in range(lo, hi):
            var b = Int(self.buckets[unsafe_offset=i]) - first_bucket
            var slot = offsets[unsafe_offset=b]
            self.entries[unsafe_offset=Int(slot)] = part[i - lo]
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
