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
from std.math import ceildiv, gcd, iota
from std.memory import pack_bits
from std.sys.intrinsics import PrefetchOptions, prefetch

from sporos.parallel import ParallelTask, parallel_for
from sporos.sort import parallel_sort
from sporos.table import AnchorSource, HashTable
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
        return self.ref_pos == other.ref_pos and self.query_pos == other.query_pos and self.length == other.length


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
            # The cost falls as k2 grows: only the largest coprime k2 counts.
            var k2 = span // k1
            while gcd(k1, k2) != 1:
                k2 -= 1
            var cost = Float64(ref_len) / Float64(k1) + 3.0 * Float64(query_len) / Float64(k2)
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
comptime _SCAN = 8
"""Entries a bucket's vector prefilter covers."""
comptime _BUCKET_LOAD = 4
"""Anchors per bucket: a quarter of the offsets, with no loss of speed
(measured on hg38 chr1 and chr21)."""
comptime _LAG = 16
"""Anchors between pipeline stages (hash -> bucket bounds -> match)."""
comptime _RING = 64
"""Pipeline ring size (a power of two above 2 * _LAG)."""


@inline(.always)
def _next_multiple(position: Int, step: Int) -> Int:
    """The smallest multiple of `step` that is >= `position` (>= 0)."""
    return (position + step - 1) // step * step


@inline(.always)
def _kmer_mask(k: Int) -> UInt64:
    return UInt64.MAX if k == 32 else (UInt64(1) << UInt64(2 * k)) - 1


struct _SampledAnchors(AnchorSource):
    """Every `step`-th reference K-mer without masked bases, keyed by its
    `mix64` hash; entries hold the hash's low 32 bits above the position."""

    var text: ImmPointer[SeqText, ImmUntrackedOrigin]
    var n: Int
    var k: Int
    var step: Int

    def __init__(out self, text: ImmPointer[SeqText, ImmUntrackedOrigin], k: Int, step: Int):
        self.text = text
        self.n = len(text[])
        self.k = k
        self.step = step

    def chunks(self) -> Int:
        return max(1, ceildiv(max(0, self.n - self.k + 1), _BUILD_CHUNK))

    def capacity(self, chunk: Int) -> Int:
        var first = chunk * _BUILD_CHUNK
        var last = min(self.n - self.k + 1, first + _BUILD_CHUNK)
        return max(0, ceildiv(last - _next_multiple(first, self.step), self.step))

    def fill(
        self,
        chunk: Int,
        hashes: MutPointer[UInt64, MutUntrackedOrigin],
        entries: MutPointer[UInt64, MutUntrackedOrigin],
    ) -> Int:
        var first = chunk * _BUILD_CHUNK
        var last = min(self.n - self.k + 1, first + _BUILD_CHUNK)
        var mask = _kmer_mask(self.k)
        ref text = self.text[]
        var slot = 0
        for start in range(_next_multiple(first, self.step), last, self.step):
            if text.mask(start) & mask != 0:
                continue
            var h = mix64(text.window(start) & mask)
            hashes[unsafe_offset=slot] = h
            entries[unsafe_offset=slot] = (h << 32) | UInt64(start)
            slot += 1
        return slot


struct MemIndex(Movable):
    """Sampled K-mer index of a reference text.

    Buckets are the top bits of each anchor's `mix64` hash; an entry holds
    the hash's low 32 bits (a fingerprint) above the anchor's position.
    """

    var params: MemParams
    var table: HashTable

    def __init__(out self, reference: SeqText, params: MemParams, threads: Int = 1) raises:
        params.validate()
        if len(reference) >= Int(UInt32.MAX):
            raise Error("reference text must be shorter than 2^32")
        self.params = params
        var source = _SampledAnchors(
            Pointer(to=reference).unsafe_origin_cast[ImmUntrackedOrigin](), params.k, params.ref_step
        )
        self.table = HashTable(source, threads, bucket_load=_BUCKET_LOAD)

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
        var off = self.table.offsets.unsafe_ptr()
        # Software pipeline over anchors: anchor n is hashed (its bucket bounds
        # prefetched), anchor n - _LAG has its bounds read (entries prefetched),
        # and anchor n - 2 * _LAG is matched.
        var positions = Array[Int, _RING](fill=0)
        var hashes = Array[UInt64, _RING](fill=0)
        var lows = Array[UInt32, _RING](fill=0)
        var highs = Array[UInt32, _RING](fill=0)
        # Candidates (fingerprint matches) wait _LAG more for their prefetched
        # reference text before being extended.
        var cands = _Candidates()
        var n = 0
        for j in range(_next_multiple(begin, qstep), min(end, len(query) - k + 1), qstep):
            if query.mask(j) & mask != 0:
                continue
            var h = mix64(query.window(j) & mask)
            var slot = n & (_RING - 1)
            positions[slot] = j
            hashes[slot] = h
            prefetch[_AHEAD](off.unsafe_offset(Int(h >> self.table.shift)))
            if n >= _LAG:
                self._bounds(n - _LAG, hashes, lows, highs)
            if n >= 2 * _LAG:
                self._match((n - 2 * _LAG) & (_RING - 1), positions, hashes, lows, highs, reference, query, cands, out)
            n += 1
        for m in range(max(0, n - _LAG), n):
            self._bounds(m, hashes, lows, highs)
        for m in range(max(0, n - 2 * _LAG), n):
            self._match(m & (_RING - 1), positions, hashes, lows, highs, reference, query, cands, out)
        while cands.tail < cands.head:
            self._extend(cands, reference, query, out)

    @inline(.always)
    def _bounds(
        self,
        n: Int,
        hashes: Array[UInt64, _RING],
        mut lows: Array[UInt32, _RING],
        mut highs: Array[UInt32, _RING],
    ):
        var slot = n & (_RING - 1)
        var b = Int(hashes[slot] >> self.table.shift)
        var off = self.table.offsets.unsafe_ptr()
        lows[slot] = off[unsafe_offset=b]
        highs[slot] = off[unsafe_offset=b + 1]
        var first = self.table.entries.unsafe_ptr().unsafe_offset(Int(lows[slot]))
        prefetch[_AHEAD](first)
        # The bucket's vector prefilter usually spans two lines.
        prefetch[_AHEAD](first.unsafe_offset(_SCAN - 1))

    @inline(.always)
    def _match(
        self,
        slot: Int,
        positions: Array[Int, _RING],
        hashes: Array[UInt64, _RING],
        lows: Array[UInt32, _RING],
        highs: Array[UInt32, _RING],
        reference: SeqText,
        query: SeqText,
        mut cands: _Candidates,
        mut out: List[Mem],
    ):
        var ent = self.table.entries.unsafe_ptr()
        var j = positions[slot]
        var fp = hashes[slot] << 32
        var lo = Int(lows[slot])
        var hi = Int(highs[slot])
        # Small buckets: one vector compare of the fingerprints picks the
        # candidates (`table.ENTRY_PADDING` keeps the load in bounds).
        var candidates = UInt64.MAX
        if hi - lo <= _SCAN:
            var lanes = ent.unsafe_offset(lo).unsafe_load[width=_SCAN]() >> 32
            var hit = lanes.eq(fp >> 32) & iota[DType.uint64, _SCAN]().lt(UInt64(hi - lo))
            candidates = UInt64(pack_bits[DType.uint8](hit))
            if candidates == 0:
                return
        for e in range(lo, hi):
            if e - lo < _SCAN and (candidates >> UInt64(e - lo)) & 1 == 0:
                continue
            var entry = ent[unsafe_offset=e]
            if (entry ^ fp) >> 32 != 0:
                continue
            var r = Int(entry & 0xFFFFFFFF)
            reference.prefetch(r)
            if cands.head - cands.tail == _CANDS:
                self._extend(cands, reference, query, out)
            cands.push(j, r)
            if cands.head - cands.tail > _LAG:
                self._extend(cands, reference, query, out)

    @inline(.always)
    def _extend(self, mut cands: _Candidates, reference: SeqText, query: SeqText, mut out: List[Mem]):
        """Extends the oldest candidate, reporting it if it starts a MEM."""
        var k = self.params.k
        var period = self.params.ref_step * self.params.query_step
        var j = cands.queries[cands.tail & (_CANDS - 1)]
        var r = cands.refs[cands.tail & (_CANDS - 1)]
        cands.tail += 1
        # Anchors of one MEM sit `period` apart on its diagonal, so this
        # anchor is the MEM's first iff the match does not reach `period`
        # bases to the left: stop looking there.
        var left = match_backward(reference, r, query, j, period)
        if left >= period:
            return
        var right = match_forward(reference, r, query, j)
        if right < k or left + right < self.params.min_len:
            return
        out.append(Mem(UInt32(r - left), UInt32(j - left), UInt32(left + right)))


comptime _CANDS = 64
"""Candidate ring size (a power of two above _LAG)."""


struct _Candidates(Movable):
    """A ring of (query, reference) anchor pairs awaiting extension."""

    var queries: Array[Int, _CANDS]
    var refs: Array[Int, _CANDS]
    var head: Int
    var tail: Int

    def __init__(out self):
        self.queries = Array[Int, _CANDS](fill=0)
        self.refs = Array[Int, _CANDS](fill=0)
        self.head = 0
        self.tail = 0

    @inline(.always)
    def push(mut self, j: Int, r: Int):
        self.queries[self.head & (_CANDS - 1)] = j
        self.refs[self.head & (_CANDS - 1)] = r
        self.head += 1


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
        sort_mems(results)


comptime _RADIX_BITS = 11


def sort_mems(mut mems: List[Mem]):
    """Sorts MEMs by (query, reference) position, which identifies each MEM:
    an LSD radix sort on `(query - lowest) << 32 | reference`, skipping the
    digits every key shares."""
    var n = len(mems)
    if n < 64:
        sort(Span(mems))
        return
    var low = Int(mems[0].query_pos)
    var high = low
    for m in mems:
        low = min(low, Int(m.query_pos))
        high = max(high, Int(m.query_pos))
    var key_bits = 32 + (64 - Int(count_leading_zeros(UInt64(high - low) | 1)))
    var base = UInt64(low)
    comptime digits = 1 << _RADIX_BITS
    var counts = List[Int](length=digits, fill=0)
    var other = List[Mem](unsafe_uninit_length=n)
    for shift in range(0, key_bits, _RADIX_BITS):
        for d in range(digits):
            counts[d] = 0
        var src = mems.unsafe_ptr()
        for i in range(n):
            counts[_digit(src[unsafe_offset=i], base, shift)] += 1
        if counts[_digit(src[unsafe_offset=0], base, shift)] == n:
            continue  # every key has this digit
        var at = 0
        for d in range(digits):
            var c = counts[d]
            counts[d] = at
            at += c
        var dst = other.unsafe_ptr()
        for i in range(n):
            var m = src[unsafe_offset=i]
            var d = _digit(m, base, shift)
            dst[unsafe_offset=counts[d]] = m
            counts[d] += 1
        var sorted = other^
        other = mems^
        mems = sorted^


@inline(.always)
def _digit(m: Mem, base: UInt64, shift: Int) -> Int:
    var key = ((UInt64(m.query_pos) - base) << 32) | UInt64(m.ref_pos)
    return Int((key >> UInt64(shift)) & ((1 << _RADIX_BITS) - 1))


def find_mems(index: MemIndex, reference: SeqText, query: SeqText, threads: Int = 1) -> List[Mem]:
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
    var total = 0
    for r in results:
        total += len(r)
    var out = List[Mem](capacity=total)
    for i in range(len(results)):
        out.extend(Span(results[i]))
        results[i] = List[Mem]()  # freed as it is copied: one copy at a time
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


def filter_mums(mems: List[Mem], query: SeqText, unique_in_query: Bool = True, threads: Int = 1) -> List[Mem]:
    """MEMs unique in the reference and (by default) in their query sequence.

    These are MUMmer's `-mum` matches; with `unique_in_query=False` they are
    its `-mumreference` matches (unique in the reference only). `mems` must
    come from one `find_mems` call (sorted by query offset); so is the result.
    """
    var keep = List[Bool](length=len(mems), fill=True)
    var lo = 0
    while lo < len(mems):
        var seq = query.sequence_of(Int(mems[lo].query_pos))
        var seq_end = query.starts[seq] + query.lengths[seq]
        var hi = lo
        while hi < len(mems) and Int(mems[hi].query_pos) < seq_end:
            hi += 1
        # Query intervals: `mems` is already in query order.
        var count = hi - lo
        var starts = List[UInt64](unsafe_uninit_length=count)
        var lengths = List[UInt32](unsafe_uninit_length=count)
        var fill = _GroupFillTask(
            mems.unsafe_ptr().unsafe_offset(lo), starts.unsafe_ptr(), lengths.unsafe_ptr(), count, by_ref=False
        )
        parallel_for(fill, _REF_ORDER_TASKS, threads)
        _mark_contained(starts, 0, lengths, keep.unsafe_ptr().unsafe_offset(lo), threads)
        if unique_in_query:
            # Reference intervals, in reference order: start << 32 | member.
            var keys = _GroupFillTask(
                mems.unsafe_ptr().unsafe_offset(lo), starts.unsafe_ptr(), lengths.unsafe_ptr(), count, by_ref=True
            )
            parallel_for(keys, _REF_ORDER_TASKS, threads)
            # Order by start only: the sweep groups equal starts in any order.
            parallel_sort(starts, threads, low_bits=32)
            # Gather lengths, sweep, scatter the verdicts: the random accesses
            # run in parallel, the sweep over contiguous arrays.
            var flags = List[Bool](length=count, fill=True)
            var gather = _RefOrderTask(
                mems.unsafe_ptr().unsafe_offset(lo), starts.unsafe_ptr().as_imm(), lengths.unsafe_ptr(), count
            )
            parallel_for(gather, _REF_ORDER_TASKS, threads)
            _mark_contained(starts, 32, lengths, flags.unsafe_ptr(), threads)
            var scatter = _ClearTask(
                starts.unsafe_ptr().as_imm(), flags.unsafe_ptr().as_imm(), keep.unsafe_ptr().unsafe_offset(lo), count
            )
            parallel_for(scatter, _REF_ORDER_TASKS, threads)
        lo = hi
    # Collect the survivors: count per block, then copy into place.
    var kept = List[Int](length=_REF_ORDER_TASKS + 1, fill=0)
    var counting = _CollectTask(
        mems.unsafe_ptr(),
        keep.unsafe_ptr().as_imm(),
        kept.unsafe_ptr(),
        Optional[MutPointer[Mem, MutUntrackedOrigin]](None),
        len(mems),
    )
    parallel_for(counting, _REF_ORDER_TASKS, threads)
    for b in range(_REF_ORDER_TASKS):
        kept[b + 1] += kept[b]
    var out = List[Mem](unsafe_uninit_length=kept[_REF_ORDER_TASKS])
    var copying = _CollectTask(
        mems.unsafe_ptr(),
        keep.unsafe_ptr().as_imm(),
        kept.unsafe_ptr(),
        Optional(out.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()),
        len(mems),
    )
    parallel_for(copying, _REF_ORDER_TASKS, threads)
    _ = len(out)
    return out^


struct _GroupFillTask[mems_origin: ImmOrigin, starts_origin: MutOrigin, lengths_origin: MutOrigin](ParallelTask):
    """Fills a group's interval arrays: query starts and lengths, or (`by_ref`)
    reference keys `ref_pos << 32 | member`."""

    var mems: ImmPointer[Mem, Self.mems_origin]
    var starts: MutPointer[UInt64, Self.starts_origin]
    var lengths: MutPointer[UInt32, Self.lengths_origin]
    var count: Int
    var by_ref: Bool

    def __init__(
        out self,
        mems: ImmPointer[Mem, Self.mems_origin],
        starts: MutPointer[UInt64, Self.starts_origin],
        lengths: MutPointer[UInt32, Self.lengths_origin],
        count: Int,
        *,
        by_ref: Bool,
    ):
        self.mems = mems
        self.starts = starts
        self.lengths = lengths
        self.count = count
        self.by_ref = by_ref

    def run(self, index: Int):
        var per = ceildiv(self.count, _REF_ORDER_TASKS)
        for i in range(index * per, min(self.count, (index + 1) * per)):
            var m = self.mems[unsafe_offset=i]
            if self.by_ref:
                self.starts[unsafe_offset=i] = (UInt64(m.ref_pos) << 32) | UInt64(i)
            else:
                self.starts[unsafe_offset=i] = UInt64(m.query_pos)
                self.lengths[unsafe_offset=i] = m.length


struct _CollectTask[mems_origin: ImmOrigin, keep_origin: ImmOrigin, kept_origin: MutOrigin](ParallelTask):
    """Counts a block's survivors (no `out`) or copies them to `out` from
    the block's offset in `kept`."""

    var mems: ImmPointer[Mem, Self.mems_origin]
    var keep: ImmPointer[Bool, Self.keep_origin]
    var kept: MutPointer[Int, Self.kept_origin]
    var out: Optional[MutPointer[Mem, MutUntrackedOrigin]]
    var count: Int

    def __init__(
        out self,
        mems: ImmPointer[Mem, Self.mems_origin],
        keep: ImmPointer[Bool, Self.keep_origin],
        kept: MutPointer[Int, Self.kept_origin],
        out_ptr: Optional[MutPointer[Mem, MutUntrackedOrigin]],
        count: Int,
    ):
        self.mems = mems
        self.keep = keep
        self.kept = kept
        self.out = out_ptr
        self.count = count

    def run(self, index: Int):
        var per = ceildiv(self.count, _REF_ORDER_TASKS)
        var lo = index * per
        var hi = min(self.count, lo + per)
        if not self.out:
            var n = 0
            for i in range(lo, hi):
                n += Int(self.keep[unsafe_offset=i])
            self.kept[unsafe_offset=index + 1] = n
            return
        var dst = self.out.value()
        var at = self.kept[unsafe_offset=index]
        for i in range(lo, hi):
            if self.keep[unsafe_offset=i]:
                dst[unsafe_offset=at] = self.mems[unsafe_offset=i]
                at += 1


comptime _REF_ORDER_TASKS = 64


struct _RefOrderTask[mems_origin: ImmOrigin, order_origin: ImmOrigin, out_origin: MutOrigin](ParallelTask):
    """Copies the lengths of a group's MEMs into reference order. Tracked
    origins keep the buffers alive while the tasks run."""

    var mems: ImmPointer[Mem, Self.mems_origin]
    var order: ImmPointer[UInt64, Self.order_origin]
    """Reference start << 32 | member, in reference order."""
    var lengths: MutPointer[UInt32, Self.out_origin]
    var count: Int

    def __init__(
        out self,
        mems: ImmPointer[Mem, Self.mems_origin],
        order: ImmPointer[UInt64, Self.order_origin],
        lengths: MutPointer[UInt32, Self.out_origin],
        count: Int,
    ):
        self.mems = mems
        self.order = order
        self.lengths = lengths
        self.count = count

    def run(self, index: Int):
        var per = ceildiv(self.count, _REF_ORDER_TASKS)
        for i in range(index * per, min(self.count, (index + 1) * per)):
            var member = Int(self.order[unsafe_offset=i] & 0xFFFFFFFF)
            self.lengths[unsafe_offset=i] = self.mems[unsafe_offset=member].length


struct _ClearTask[order_origin: ImmOrigin, flags_origin: ImmOrigin, keep_origin: MutOrigin](ParallelTask):
    """Clears `keep` for the members whose `flags`, in reference order, are
    cleared."""

    var order: ImmPointer[UInt64, Self.order_origin]
    var flags: ImmPointer[Bool, Self.flags_origin]
    var keep: MutPointer[Bool, Self.keep_origin]
    var count: Int

    def __init__(
        out self,
        order: ImmPointer[UInt64, Self.order_origin],
        flags: ImmPointer[Bool, Self.flags_origin],
        keep: MutPointer[Bool, Self.keep_origin],
        count: Int,
    ):
        self.order = order
        self.flags = flags
        self.keep = keep
        self.count = count

    def run(self, index: Int):
        var per = ceildiv(self.count, _REF_ORDER_TASKS)
        for i in range(index * per, min(self.count, (index + 1) * per)):
            if not self.flags[unsafe_offset=i]:
                self.keep[unsafe_offset=Int(self.order[unsafe_offset=i] & 0xFFFFFFFF)] = False


def _mark_contained(
    starts: List[UInt64], shift: UInt64, lengths: List[UInt32], keep: MutPointer[Bool, _], threads: Int = 1
):
    """Clears `keep[i]` for the intervals `[starts[i] >> shift, ... + lengths[i])`,
    in increasing start, contained in (or equal to) another one.

    Of the intervals sharing a start, all but the longest are contained in
    it, and the longest is when it has a copy or an earlier start reaches as
    far. In parallel: blocks starting at group boundaries find their largest
    end, a prefix maximum gives each block what precedes it, and the blocks
    then sweep independently."""
    var count = len(starts)
    if threads <= 1 or count < 1 << 16:
        _ = _sweep(Span(starts), shift, Span(lengths), keep, 0, count, -1, True)
        return
    var blocks = threads * 4
    var bounds = List[Int](length=blocks + 1, fill=count)
    bounds[0] = 0
    for b in range(1, blocks):
        var at = max(bounds[b - 1], b * count // blocks)
        while at < count and at > 0 and starts[at] >> shift == starts[at - 1] >> shift:
            at += 1
        bounds[b] = at
    var ends = List[Int](length=blocks, fill=-1)
    var scan = _SweepTask(
        starts.unsafe_ptr().as_imm(),
        shift,
        lengths.unsafe_ptr().as_imm(),
        keep.unsafe_origin_cast[MutUntrackedOrigin](),
        bounds.unsafe_ptr().as_imm(),
        ends.unsafe_ptr(),
        count,
        mark=False,
    )
    parallel_for(scan, blocks, threads)
    # Each block's incoming max end: the largest end of all earlier blocks.
    var before = -1
    for b in range(blocks):
        var mine = ends[b]
        ends[b] = before
        before = max(before, mine)
    var mark = _SweepTask(
        starts.unsafe_ptr().as_imm(),
        shift,
        lengths.unsafe_ptr().as_imm(),
        keep.unsafe_origin_cast[MutUntrackedOrigin](),
        bounds.unsafe_ptr().as_imm(),
        ends.unsafe_ptr(),
        count,
        mark=True,
    )
    parallel_for(mark, blocks, threads)
    _ = len(ends)


def _sweep(
    starts: Span[UInt64, _],
    shift: UInt64,
    lengths: Span[UInt32, _],
    keep: MutPointer[Bool, _],
    lo: Int,
    hi: Int,
    incoming: Int,
    mark: Bool,
) -> Int:
    """Sweeps groups `[lo, hi)` (whole groups) after intervals ending at
    `incoming`; with `mark`, clears `keep` for the contained ones. Returns
    the largest end seen."""
    var max_end = incoming
    var i = lo
    while i < hi:
        var start = Int(starts[i] >> shift)
        var j = i
        var longest = 0
        var copies = 0
        while j < hi and Int(starts[j] >> shift) == start:
            var length = Int(lengths[j])
            if length > longest:
                longest = length
                copies = 0
            copies += Int(length == longest)
            j += 1
        var end = start + longest
        if mark:
            var dropped = copies > 1 or end <= max_end
            for t in range(i, j):
                if dropped or Int(lengths[t]) < longest:
                    keep[unsafe_offset=t] = False
        max_end = max(max_end, end)
        i = j
    return max_end


struct _SweepTask[
    starts_origin: ImmOrigin, lengths_origin: ImmOrigin, bounds_origin: ImmOrigin, ends_origin: MutOrigin
](ParallelTask):
    """One block of a parallel `_mark_contained`: its largest end, or (`mark`)
    its sweep after the block's incoming end in `ends`."""

    var starts: ImmPointer[UInt64, Self.starts_origin]
    var shift: UInt64
    var lengths: ImmPointer[UInt32, Self.lengths_origin]
    var keep: MutPointer[Bool, MutUntrackedOrigin]
    var bounds: ImmPointer[Int, Self.bounds_origin]
    var ends: MutPointer[Int, Self.ends_origin]
    var count: Int
    var mark: Bool

    def __init__(
        out self,
        starts: ImmPointer[UInt64, Self.starts_origin],
        shift: UInt64,
        lengths: ImmPointer[UInt32, Self.lengths_origin],
        keep: MutPointer[Bool, MutUntrackedOrigin],
        bounds: ImmPointer[Int, Self.bounds_origin],
        ends: MutPointer[Int, Self.ends_origin],
        count: Int,
        *,
        mark: Bool,
    ):
        self.starts = starts
        self.shift = shift
        self.lengths = lengths
        self.keep = keep
        self.bounds = bounds
        self.ends = ends
        self.count = count
        self.mark = mark

    def run(self, index: Int):
        var lo = self.bounds[unsafe_offset=index]
        var hi = self.bounds[unsafe_offset=index + 1]
        var starts = Span(unsafe_ptr=self.starts, length=self.count)
        var lengths = Span(unsafe_ptr=self.lengths, length=self.count)
        var end = _sweep(starts, self.shift, lengths, self.keep, lo, hi, self.ends[unsafe_offset=index], self.mark)
        if not self.mark:
            self.ends[unsafe_offset=index] = end
