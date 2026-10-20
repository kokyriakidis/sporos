"""Maximal exact matches (MEMs) and maximal unique matches (MUMs).

Sampled k-mer anchoring (as in copMEM): reference K-mers are indexed every
`k1` positions and query K-mers probed every `k2` positions, with `k1`, `k2`
coprime and `k1 * k2 <= L - K + 1`. By the Chinese remainder theorem every
match of length >= L contains an anchor pair, and the anchor offsets inside a
MEM form a single residue class modulo `k1 * k2`. A MEM is therefore reported
exactly once, from the anchor whose left extension is shorter than `k1 * k2`,
with no hashing or sorting of results.

A MEM is a MUM when its string occurs exactly once in the reference and once
in its query sequence. Any other occurrence would be covered by a different
MEM (of length >= L) whose reference or query interval contains this one, so
MUMs are exactly the MEMs not contained, in either coordinate, by another MEM
of the same query sequence.
"""

from std.bit import count_leading_zeros, count_trailing_zeros
from std.math import ceildiv, gcd
from std.memory import Layout, alloc, dealloc
from std.sys.intrinsics import PrefetchOptions, prefetch

from sporos.parallel import ParallelTask, parallel_for
from sporos.text import SeqText


@fieldwise_init
struct Mem(Comparable, ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """A maximal exact match between reference and query text offsets."""

    var ref_pos: UInt32
    var query_pos: UInt32
    var length: UInt32

    def __lt__(self, other: Self) -> Bool:
        if self.query_pos != other.query_pos:
            return self.query_pos < other.query_pos
        if self.ref_pos != other.ref_pos:
            return self.ref_pos < other.ref_pos
        return self.length < other.length

    def __eq__(self, other: Self) -> Bool:
        return (
            self.ref_pos == other.ref_pos
            and self.query_pos == other.query_pos
            and self.length == other.length
        )


@fieldwise_init
struct MemParams(ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """Anchor sampling: K-mer length and reference/query sampling steps."""

    var min_len: Int
    var k: Int
    var ref_step: Int
    var query_step: Int

    @staticmethod
    def choose(min_len: Int, ref_len: Int, query_len: Int) -> Self:
        """Picks K and coprime steps for `min_len`, minimizing estimated work.

        Index construction costs one sequential insert per sampled reference
        position; each sampled query position costs a random-access probe,
        weighted 3x.
        """
        var k = min(min_len, max(16, min(32, (min_len + 1) // 2 + 4)))
        var span = min_len - k + 1
        var best = Self(min_len, k, 1, 1)
        var best_cost = Float64.MAX
        for k1 in range(1, span + 1):
            for k2 in range(1, span // k1 + 1):
                if gcd(k1, k2) != 1:
                    continue
                var cost = Float64(ref_len) / Float64(k1) + 3.0 * Float64(
                    query_len
                ) / Float64(k2)
                if cost < best_cost:
                    best_cost = cost
                    best = Self(min_len, k, k1, k2)
        return best

    def validate(self) raises:
        if self.min_len < 1:
            raise Error(t"minimum match length must be >= 1, got {self.min_len}")
        if self.k < 1 or self.k > 32 or self.k > self.min_len:
            raise Error(t"anchor k must be in [1, min(32, L)], got {self.k}")
        if self.ref_step < 1 or self.query_step < 1:
            raise Error("sampling steps must be >= 1")
        if gcd(self.ref_step, self.query_step) != 1:
            raise Error("sampling steps must be coprime")
        if self.ref_step * self.query_step > self.min_len - self.k + 1:
            raise Error("ref_step * query_step must be <= L - k + 1")


@inline(.always)
def mix64(x: UInt64) -> UInt64:
    """MurmurHash3 finalizer: an invertible 64-bit mixer."""
    var h = x
    h ^= h >> 33
    h *= 0xFF51AFD7ED558CCD
    h ^= h >> 33
    h *= 0xC4CEB9FE1A85EC53
    h ^= h >> 33
    return h


@inline(.always)
def match_forward(a: SeqText, i: Int, b: SeqText, j: Int) -> Int:
    """Length of the common prefix of `a[i:]` and `b[j:]`, 32 bases at a time."""
    var n = 0
    while True:
        var d = (a.window(i + n) ^ b.window(j + n)) | a.mask(i + n) | b.mask(j + n)
        if d != 0:
            return n + Int(count_trailing_zeros(d) >> 1)
        n += 32


@inline(.always)
def match_backward(a: SeqText, i: Int, b: SeqText, j: Int, limit: Int) -> Int:
    """Length of the common suffix of `a[:i]` and `b[:j]`, or some value
    >= `limit` once it reaches `limit`."""
    var n = 0
    while n < limit:
        var at = i - n - 32
        var bt = j - n - 32
        var d = (a.window(at) ^ b.window(bt)) | a.mask(at) | b.mask(bt)
        if d != 0:
            return n + Int(count_leading_zeros(d) >> 1)
        n += 32
    return n


comptime _AHEAD = PrefetchOptions().for_read().high_locality()
comptime _BUILD_CHUNK = 1 << 20
"""Reference positions per anchor-collection task."""
comptime _PROBE_BATCH = 32
"""Query anchors hashed (and their buckets prefetched) before any is probed."""


@inline(.always)
def _next_multiple(position: Int, step: Int) -> Int:
    """The smallest multiple of `step` that is >= `position` (>= 0)."""
    return (position + step - 1) // step * step


@inline(.always)
def _kmer_mask(k: Int) -> UInt64:
    return UInt64.MAX if k == 32 else (UInt64(1) << UInt64(2 * k)) - 1


struct _AnchorCounts(ParallelTask):
    """Counts the indexable reference anchors of each chunk."""

    var text: ImmPointer[SeqText, ImmUntrackedOrigin]
    var counts: MutPointer[Int, MutUntrackedOrigin]
    var n: Int
    var k: Int
    var step: Int

    def __init__(
        out self,
        text: ImmPointer[SeqText, ImmUntrackedOrigin],
        counts: MutPointer[Int, MutUntrackedOrigin],
        n: Int,
        k: Int,
        step: Int,
    ):
        self.text = text
        self.counts = counts
        self.n = n
        self.k = k
        self.step = step

    def run(self, index: Int):
        var first = index * _BUILD_CHUNK
        var last = min(self.n - self.k + 1, first + _BUILD_CHUNK)
        var kmask = _kmer_mask(self.k)
        var count = 0
        for start in range(_next_multiple(first, self.step), last, self.step):
            count += Int(self.text[].mask(start) & kmask == 0)
        self.counts[unsafe_offset=index] = count


struct _AnchorFill(ParallelTask):
    """Writes each chunk's anchors (bucket, fingerprint << 32 | position) and
    counts them per bucket part (the bucket's top `part_bits` bits)."""

    var text: ImmPointer[SeqText, ImmUntrackedOrigin]
    var firsts: ImmPointer[Int, ImmUntrackedOrigin]
    var buckets: MutPointer[UInt32, MutUntrackedOrigin]
    var anchors: MutPointer[UInt64, MutUntrackedOrigin]
    var histogram: MutPointer[Int, MutUntrackedOrigin]
    var n: Int
    var k: Int
    var step: Int
    var shift: UInt64
    var parts: Int
    var part_shift: UInt32

    def __init__(
        out self,
        text: ImmPointer[SeqText, ImmUntrackedOrigin],
        firsts: ImmPointer[Int, ImmUntrackedOrigin],
        buckets: MutPointer[UInt32, MutUntrackedOrigin],
        anchors: MutPointer[UInt64, MutUntrackedOrigin],
        histogram: MutPointer[Int, MutUntrackedOrigin],
        n: Int,
        k: Int,
        step: Int,
        shift: UInt64,
        parts: Int,
        part_shift: UInt32,
    ):
        self.text = text
        self.firsts = firsts
        self.buckets = buckets
        self.anchors = anchors
        self.histogram = histogram
        self.n = n
        self.k = k
        self.step = step
        self.shift = shift
        self.parts = parts
        self.part_shift = part_shift

    def run(self, index: Int):
        var first = index * _BUILD_CHUNK
        var last = min(self.n - self.k + 1, first + _BUILD_CHUNK)
        var mask = _kmer_mask(self.k)
        var slot = self.firsts[unsafe_offset=index]
        var histogram = self.histogram.unsafe_offset(index * self.parts)
        for p in range(self.parts):
            histogram[unsafe_offset=p] = 0
        ref text = self.text[]
        for start in range(_next_multiple(first, self.step), last, self.step):
            if text.mask(start) & mask != 0:
                continue
            var h = mix64(text.window(start) & mask)
            var bucket = UInt32(h >> self.shift)
            self.buckets[unsafe_offset=slot] = bucket
            self.anchors[unsafe_offset=slot] = (h << 32) | UInt64(start)
            histogram[unsafe_offset=Int(bucket >> self.part_shift)] += 1
            slot += 1


struct _PartScatter(ParallelTask):
    """Moves one chunk's anchors to their parts (part-major, position order
    kept within each part); `cursors` holds the chunk's start in each part."""

    var firsts: ImmPointer[Int, ImmUntrackedOrigin]
    var buckets: ImmPointer[UInt32, ImmUntrackedOrigin]
    var anchors: ImmPointer[UInt64, ImmUntrackedOrigin]
    var cursors: MutPointer[Int, MutUntrackedOrigin]
    var part_buckets: MutPointer[UInt32, MutUntrackedOrigin]
    var part_anchors: MutPointer[UInt64, MutUntrackedOrigin]
    var parts: Int
    var part_shift: UInt32

    def __init__(
        out self,
        firsts: ImmPointer[Int, ImmUntrackedOrigin],
        buckets: ImmPointer[UInt32, ImmUntrackedOrigin],
        anchors: ImmPointer[UInt64, ImmUntrackedOrigin],
        cursors: MutPointer[Int, MutUntrackedOrigin],
        part_buckets: MutPointer[UInt32, MutUntrackedOrigin],
        part_anchors: MutPointer[UInt64, MutUntrackedOrigin],
        parts: Int,
        part_shift: UInt32,
    ):
        self.firsts = firsts
        self.buckets = buckets
        self.anchors = anchors
        self.cursors = cursors
        self.part_buckets = part_buckets
        self.part_anchors = part_anchors
        self.parts = parts
        self.part_shift = part_shift

    def run(self, index: Int):
        var cursors = self.cursors.unsafe_offset(index * self.parts)
        for i in range(self.firsts[unsafe_offset=index], self.firsts[unsafe_offset=index + 1]):
            var bucket = self.buckets[unsafe_offset=i]
            var p = Int(bucket >> self.part_shift)
            var slot = cursors[unsafe_offset=p]
            self.part_buckets[unsafe_offset=slot] = bucket
            self.part_anchors[unsafe_offset=slot] = self.anchors[unsafe_offset=i]
            cursors[unsafe_offset=p] = slot + 1


struct _PartSort(ParallelTask):
    """Counting sort of one part's anchors by bucket into the final table."""

    var part_starts: ImmPointer[Int, ImmUntrackedOrigin]
    var part_buckets: ImmPointer[UInt32, ImmUntrackedOrigin]
    var part_anchors: ImmPointer[UInt64, ImmUntrackedOrigin]
    var offsets: MutPointer[UInt32, MutUntrackedOrigin]
    var entries: MutPointer[UInt64, MutUntrackedOrigin]
    var buckets_per_part: Int

    def __init__(
        out self,
        part_starts: ImmPointer[Int, ImmUntrackedOrigin],
        part_buckets: ImmPointer[UInt32, ImmUntrackedOrigin],
        part_anchors: ImmPointer[UInt64, ImmUntrackedOrigin],
        offsets: MutPointer[UInt32, MutUntrackedOrigin],
        entries: MutPointer[UInt64, MutUntrackedOrigin],
        buckets_per_part: Int,
    ):
        self.part_starts = part_starts
        self.part_buckets = part_buckets
        self.part_anchors = part_anchors
        self.offsets = offsets
        self.entries = entries
        self.buckets_per_part = buckets_per_part

    def run(self, index: Int):
        var lo = self.part_starts[unsafe_offset=index]
        var hi = self.part_starts[unsafe_offset=index + 1]
        var first_bucket = index * self.buckets_per_part
        var offsets = self.offsets.unsafe_offset(first_bucket)
        for b in range(self.buckets_per_part):
            offsets[unsafe_offset=b] = 0
        for i in range(lo, hi):
            offsets[unsafe_offset=Int(self.part_buckets[unsafe_offset=i]) - first_bucket] += 1
        var at = UInt32(lo)
        for b in range(self.buckets_per_part):
            var count = offsets[unsafe_offset=b]
            offsets[unsafe_offset=b] = at
            at += count
        # offsets[b] advances to the end of bucket b while scattering, which
        # is the start of bucket b + 1; shift back afterwards.
        for i in range(lo, hi):
            var b = Int(self.part_buckets[unsafe_offset=i]) - first_bucket
            var slot = offsets[unsafe_offset=b]
            self.entries[unsafe_offset=Int(slot)] = self.part_anchors[unsafe_offset=i]
            offsets[unsafe_offset=b] = slot + 1
        for b in reversed(range(1, self.buckets_per_part)):
            offsets[unsafe_offset=b] = offsets[unsafe_offset=b - 1]
        offsets[unsafe_offset=0] = UInt32(lo)


struct MemIndex(Movable):
    """Sampled K-mer index of a reference text.

    Buckets are the top bits of each anchor's `mix64` hash; an entry holds
    the hash's low 32 bits (a fingerprint) above the anchor's position.
    """

    var params: MemParams
    var offsets: List[UInt32]
    var entries: List[UInt64]
    var shift: UInt64

    def __init__(
        out self, reference: SeqText, params: MemParams, threads: Int = 1
    ) raises:
        params.validate()
        if len(reference) >= Int(UInt32.MAX):
            raise Error("reference text must be shorter than 2^32")
        self.params = params
        var k = params.k
        var step = params.ref_step
        var n = len(reference)
        var text = Pointer(to=reference).unsafe_origin_cast[ImmUntrackedOrigin]()

        # Anchors per chunk, then each chunk's anchors at its offset. The
        # tasks reach these buffers through untracked pointers, so they are
        # explicit allocations freed at the end.
        var chunks = max(1, ceildiv(max(0, n - k + 1), _BUILD_CHUNK))
        var firsts_alloc = alloc(Layout[Int](length=chunks + 1))
        var firsts = firsts_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        for c in range(chunks + 1):
            firsts[unsafe_offset=c] = 0
        var count_task = _AnchorCounts(text, firsts.unsafe_offset(1), n, k, step)
        parallel_for(count_task, chunks, threads)
        for c in range(chunks):
            firsts[unsafe_offset=c + 1] += firsts[unsafe_offset=c]
        var total = firsts[unsafe_offset=chunks]

        var bits = 1
        while (1 << bits) < max(1, total):
            bits += 1
        self.shift = UInt64(64 - bits)
        var num_buckets = 1 << bits

        # Radix by the top bits of the bucket: fill chunks while counting
        # their parts, scatter to part-major order, then sort each part.
        var part_bits = min(8, bits)
        var parts = 1 << part_bits
        var part_shift = UInt32(bits - part_bits)
        var buckets_alloc = alloc(Layout[UInt32](length=max(1, total)))
        var anchors_alloc = alloc(Layout[UInt64](length=max(1, total)))
        var histogram_alloc = alloc(Layout[Int](length=chunks * parts))
        var buckets = buckets_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        var anchors = anchors_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        var histogram = histogram_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        var fill_task = _AnchorFill(
            text, firsts, buckets, anchors, histogram, n, k, step, self.shift, parts, part_shift
        )
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
        var part_buckets_alloc = alloc(Layout[UInt32](length=max(1, total)))
        var part_anchors_alloc = alloc(Layout[UInt64](length=max(1, total)))
        var part_buckets = part_buckets_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        var part_anchors = part_anchors_alloc.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        var scatter_task = _PartScatter(
            firsts, buckets, anchors, histogram, part_buckets, part_anchors, parts, part_shift
        )
        parallel_for(scatter_task, chunks, threads)
        dealloc(buckets_alloc^)
        dealloc(anchors_alloc^)
        dealloc(histogram_alloc^)

        self.offsets = List[UInt32](unsafe_uninit_length=num_buckets + 1)
        self.entries = List[UInt64](unsafe_uninit_length=total)
        var sort_task = _PartSort(
            part_starts,
            part_buckets,
            part_anchors,
            self.offsets.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
            self.entries.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
            num_buckets // parts,
        )
        parallel_for(sort_task, parts, threads)
        self.offsets[num_buckets] = UInt32(total)
        dealloc(firsts_alloc^)
        dealloc(starts_alloc^)
        dealloc(part_buckets_alloc^)
        dealloc(part_anchors_alloc^)

    def find(
        self,
        reference: SeqText,
        query: SeqText,
        begin: Int,
        end: Int,
        mut out: List[Mem],
    ):
        """Appends MEMs whose sampled query anchor lies in `[begin, end)`."""
        var k = self.params.k
        var qstep = self.params.query_step
        var mask = _kmer_mask(k)
        var batch_pos = Array[Int, _PROBE_BATCH](fill=0)
        var batch_hash = Array[UInt64, _PROBE_BATCH](fill=0)
        var batched = 0
        for j in range(_next_multiple(begin, qstep), min(end, len(query) - k + 1), qstep):
            if query.mask(j) & mask != 0:
                continue
            var h = mix64(query.window(j) & mask)
            prefetch[_AHEAD](self.offsets.unsafe_ptr().unsafe_offset(Int(h >> self.shift)))
            batch_pos[batched] = j
            batch_hash[batched] = h
            batched += 1
            if batched == _PROBE_BATCH:
                self._probe_batch(reference, query, batch_pos, batch_hash, batched, out)
                batched = 0
        self._probe_batch(reference, query, batch_pos, batch_hash, batched, out)

    @inline(.always)
    def _probe_batch(
        self,
        reference: SeqText,
        query: SeqText,
        positions: Array[Int, _PROBE_BATCH],
        hashes: Array[UInt64, _PROBE_BATCH],
        count: Int,
        mut out: List[Mem],
    ):
        """Probes a batch: bucket bounds (prefetching entries), then matches."""
        var off = self.offsets.unsafe_ptr()
        var ent = self.entries.unsafe_ptr()
        var lows = Array[UInt32, _PROBE_BATCH](fill=0)
        var highs = Array[UInt32, _PROBE_BATCH](fill=0)
        for a in range(count):
            var b = Int(hashes[a] >> self.shift)
            lows[a] = off[unsafe_offset=b]
            highs[a] = off[unsafe_offset=b + 1]
            prefetch[_AHEAD](ent.unsafe_offset(Int(lows[a])))
        var k = self.params.k
        var period = self.params.ref_step * self.params.query_step
        var min_len = self.params.min_len
        for a in range(count):
            var j = positions[a]
            var fp = hashes[a] << 32
            for e in range(Int(lows[a]), Int(highs[a])):
                var entry = ent[unsafe_offset=e]
                if (entry ^ fp) >> 32 != 0:
                    continue
                var r = Int(entry & 0xFFFFFFFF)
                # Anchors of one MEM sit `period` apart on its diagonal, so
                # this anchor is the MEM's first iff the match does not reach
                # `period` bases to the left: stop looking there.
                var left = match_backward(reference, r, query, j, period)
                if left >= period:
                    continue
                var right = match_forward(reference, r, query, j)
                if right < k or left + right < min_len:
                    continue
                out.append(
                    Mem(UInt32(r - left), UInt32(j - left), UInt32(left + right))
                )


struct _FindTask(ParallelTask):
    var index: ImmPointer[MemIndex, ImmUntrackedOrigin]
    var reference: ImmPointer[SeqText, ImmUntrackedOrigin]
    var query: ImmPointer[SeqText, ImmUntrackedOrigin]
    var results: MutPointer[List[Mem], MutUntrackedOrigin]
    var chunk: Int

    def __init__(
        out self,
        index: ImmPointer[MemIndex, ImmUntrackedOrigin],
        reference: ImmPointer[SeqText, ImmUntrackedOrigin],
        query: ImmPointer[SeqText, ImmUntrackedOrigin],
        results: MutPointer[List[Mem], MutUntrackedOrigin],
        chunk: Int,
    ):
        self.index = index
        self.reference = reference
        self.query = query
        self.results = results
        self.chunk = chunk

    def run(self, index: Int):
        var begin = index * self.chunk
        var end = min(len(self.query[]), begin + self.chunk)
        ref results = self.results[unsafe_offset=index]
        self.index[].find(self.reference[], self.query[], begin, end, results)
        sort(Span(results))


def find_mems(
    index: MemIndex, reference: SeqText, query: SeqText, threads: Int = 1
) -> List[Mem]:
    """All MEMs of length >= `index.params.min_len`, sorted by query offset."""
    comptime chunk_target = 1 << 18
    var n = len(query)
    var chunks = max(1, ceildiv(n, chunk_target)) if threads > 1 else 1
    var chunk = ceildiv(max(n, 1), chunks)
    var results = List[List[Mem]](length=chunks, fill=List[Mem]())
    var task = _FindTask(
        Pointer(to=index).unsafe_origin_cast[ImmUntrackedOrigin](),
        Pointer(to=reference).unsafe_origin_cast[ImmUntrackedOrigin](),
        Pointer(to=query).unsafe_origin_cast[ImmUntrackedOrigin](),
        results.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        chunk,
    )
    parallel_for(task, chunks, threads)
    var out = List[Mem]()
    for r in results:
        out.extend(Span(r))
    # Chunks are sorted and cover consecutive query ranges; only MEMs that
    # start left of their anchor's chunk can be out of place, by little.
    for i in range(1, len(out)):
        var m = out[i]
        var j = i
        while j > 0 and m < out[j - 1]:
            out[j] = out[j - 1]
            j -= 1
        out[j] = m
    return out^


def filter_mums(
    mems: List[Mem], query: SeqText, unique_in_query: Bool = True
) -> List[Mem]:
    """MEMs unique in the reference and (by default) in their query sequence.

    These are MUMmer's `-mum` matches; with `unique_in_query=False` they are
    its `-mumreference` matches (unique in the reference only). `mems` must
    come from one `find_mems` call (sorted by query offset).
    """
    var out = List[Mem]()
    var lo = 0
    while lo < len(mems):
        var seq = query.sequence_of(Int(mems[lo].query_pos))
        var seq_end = query.starts[seq] + query.lengths[seq]
        var hi = lo
        while hi < len(mems) and Int(mems[hi].query_pos) < seq_end:
            hi += 1
        var group = List[Mem](capacity=hi - lo)
        for i in range(lo, hi):
            group.append(mems[i])
        var keep = List[Bool](length=len(group), fill=True)
        _mark_contained[by_ref=False](group, keep)
        if unique_in_query:
            _mark_contained[by_ref=True](group, keep)
        for i in range(len(group)):
            if keep[i]:
                out.append(group[i])
        lo = hi
    sort(Span(out))
    return out^


@inline(.always)
def _interval_start[by_ref: Bool](m: Mem) -> Int:
    comptime if by_ref:
        return Int(m.ref_pos)
    else:
        return Int(m.query_pos)


def _mark_contained[by_ref: Bool](group: List[Mem], mut keep: List[Bool]):
    """Clears `keep` for intervals contained in (or equal to) another one."""
    var order = List[Int](capacity=len(group))
    for i in range(len(group)):
        order.append(i)

    def before(a: Int, b: Int) {imm group} -> Bool:
        var sa = _interval_start[by_ref](group[a])
        var sb = _interval_start[by_ref](group[b])
        if sa != sb:
            return sa < sb
        return group[a].length > group[b].length

    sort(Span(order), before)
    var max_end = -1
    for idx in range(len(order)):
        var m = group[order[idx]]
        var s = _interval_start[by_ref](m)
        var e = s + Int(m.length)
        if e <= max_end:
            keep[order[idx]] = False
        # An identical interval right after this one makes both non-unique.
        if idx + 1 < len(order):
            var n = group[order[idx + 1]]
            if _interval_start[by_ref](n) == s and n.length == m.length:
                keep[order[idx]] = False
        max_end = max(max_end, e)
