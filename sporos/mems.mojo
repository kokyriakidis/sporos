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


@inline(.always)
def _ref_forward(reference: SeqText, i: Int, query: SeqText, j: Int, stop: Int) -> Int:
    """`match_forward` from reference `i`, whose unmasked stretch ends at
    `stop`: reference masks are replaced by that bound."""
    var n = 0
    while True:
        var d = (reference.window(i + n) ^ query.window(j + n)) | query.mask(j + n)
        if d != 0:
            return min(n + Int(count_trailing_zeros(d) >> 1), stop - i)
        n += 32
        if i + n >= stop:
            return stop - i


@inline(.always)
def _ref_backward(reference: SeqText, i: Int, query: SeqText, j: Int, limit: Int, start: Int) -> Int:
    """`match_backward` to reference `i`, whose unmasked stretch starts at
    `start`."""
    var n = 0
    while n < limit and i - n > start:
        var at = i - n - 32
        var bt = j - n - 32
        var d = (reference.window(at) ^ query.window(bt)) | query.mask(bt)
        if d != 0:
            return min(n + Int(count_leading_zeros(d) >> 1), i - start)
        n += 32
    return min(n, i - start)


comptime _BLOCK_BITS = 12


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
    var runs: List[Int]
    """The reference's masked runs (`SeqText.masked_runs`): extensions clamp
    to the unmasked stretch around an anchor instead of reading masks."""
    var length: Int
    var block_lo: List[Int]
    """Per 2^_BLOCK_BITS bases: the unmasked stretch holding the block
    (`block_lo`, `block_hi`), or `block_lo = -1` if masked bases split it."""
    var block_hi: List[Int]

    def __init__(out self, reference: SeqText, params: MemParams, threads: Int = 1) raises:
        params.validate()
        if len(reference) >= Int(UInt32.MAX):
            raise Error("reference text must be shorter than 2^32")
        self.params = params
        var source = _SampledAnchors(
            Pointer(to=reference).unsafe_origin_cast[ImmUntrackedOrigin](), params.k, params.ref_step
        )
        self.table = HashTable(source, threads, bucket_load=_BUCKET_LOAD)
        self.runs = reference.masked_runs()
        self.length = len(reference)
        var blocks = (self.length >> _BLOCK_BITS) + 1
        self.block_lo = List[Int](length=blocks, fill=-1)
        self.block_hi = List[Int](length=blocks, fill=-1)
        for b in range(blocks):
            var first = b << _BLOCK_BITS
            var last = min(self.length, first + (1 << _BLOCK_BITS)) - 1
            if first > last:
                continue
            var stretch = self._search(first)
            if stretch[0] <= first and stretch[1] > last:
                self.block_lo[b] = stretch[0]
                self.block_hi[b] = stretch[1]

    @inline(.always)
    def _unmasked(self, r: Int) -> Tuple[Int, Int]:
        """The unmasked stretch `[lo, hi)` of the reference around `r`
        (unmasked)."""
        var b = r >> _BLOCK_BITS
        var lo = self.block_lo[b]
        if lo >= 0:
            return (lo, self.block_hi[b])
        return self._search(r)

    def _search(self, r: Int) -> Tuple[Int, Int]:
        """`_unmasked` by binary search over the masked runs. For a masked
        `r`, a stretch that does not hold it."""
        var count = len(self.runs) // 2
        var a = 0
        var b = count
        while a < b:
            var mid = (a + b) // 2
            if self.runs[2 * mid] > r:
                b = mid
            else:
                a = mid + 1
        var hi = self.runs[2 * a] if a < count else self.length
        var lo = self.runs[2 * a - 1] if a > 0 else 0
        return (lo, hi)

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
            reference.prefetch_words(r)
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
        var stretch = self._unmasked(r)
        var left = _ref_backward(reference, r, query, j, period, stretch[0])
        if left >= period:
            return
        var right = _ref_forward(reference, r, query, j, stretch[1])
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
comptime _RADIX_MIN = 1 << 8
"""Below this many MEMs an insertion sort beats the radix passes."""


def sort_mems(mut mems: List[Mem]):
    """Sorts MEMs by (query, reference) position, which identifies each MEM:
    an LSD radix sort on `(query - lowest) << 32 | reference`, skipping the
    digits every key shares."""
    var n = len(mems)
    if n < _RADIX_MIN:
        # Few MEMs, found nearly in order: an insertion sort.
        for i in range(1, n):
            var m = mems[i]
            var j = i
            while j > 0 and m < mems[j - 1]:
                mems[j] = mems[j - 1]
                j -= 1
            mems[j] = m
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
    # Threads pay off only for many MEMs; per-read calls stay serial.
    var workers = threads if len(mems) >= _PARALLEL_MEMS else 1
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
        parallel_for(fill, _REF_ORDER_TASKS, workers)
        _mark_contained(starts, 0, lengths, keep.unsafe_ptr().unsafe_offset(lo), workers)
        if unique_in_query:
            # Reference intervals, in reference order: start << 32 | member.
            var keys = _GroupFillTask(
                mems.unsafe_ptr().unsafe_offset(lo), starts.unsafe_ptr(), lengths.unsafe_ptr(), count, by_ref=True
            )
            parallel_for(keys, _REF_ORDER_TASKS, workers)
            # Order by start only: the sweep groups equal starts in any order.
            parallel_sort(starts, workers, low_bits=32)
            # Gather lengths, sweep, scatter the verdicts: the random accesses
            # run in parallel, the sweep over contiguous arrays.
            var flags = List[Bool](length=count, fill=True)
            var gather = _RefOrderTask(
                mems.unsafe_ptr().unsafe_offset(lo), starts.unsafe_ptr().as_imm(), lengths.unsafe_ptr(), count
            )
            parallel_for(gather, _REF_ORDER_TASKS, workers)
            _mark_contained(starts, 32, lengths, flags.unsafe_ptr(), workers)
            var scatter = _ClearTask(
                starts.unsafe_ptr().as_imm(), flags.unsafe_ptr().as_imm(), keep.unsafe_ptr().unsafe_offset(lo), count
            )
            parallel_for(scatter, _REF_ORDER_TASKS, workers)
        lo = hi
    return _collect(mems, keep, workers)


def _collect(mems: List[Mem], keep: List[Bool], workers: Int) -> List[Mem]:
    """The MEMs whose `keep` is set, in order: counted per block, then
    copied into place."""
    var kept = List[Int](length=_REF_ORDER_TASKS + 1, fill=0)
    var counting = _CollectTask(
        mems.unsafe_ptr(),
        keep.unsafe_ptr().as_imm(),
        kept.unsafe_ptr(),
        Optional[MutPointer[Mem, MutUntrackedOrigin]](None),
        len(mems),
    )
    parallel_for(counting, _REF_ORDER_TASKS, workers)
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
    parallel_for(copying, _REF_ORDER_TASKS, workers)
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


@fieldwise_init
struct Smem(Comparable, ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """A supermaximal exact match: query bases `[start, end)` (forward-strand
    text offsets) occur in the reference at `Smems.hits[first : first + count]`."""

    var start: UInt32
    var end: UInt32
    var first: Int
    var count: Int

    def __lt__(self, other: Self) -> Bool:
        return self.start < other.start or (self.start == other.start and self.end < other.end)

    def __eq__(self, other: Self) -> Bool:
        return self.start == other.start and self.end == other.end and self.count == other.count


@fieldwise_init
struct SmemHit(Comparable, ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """Where an SMEM occurs: the reference text offset of its leftmost base,
    and whether the reference holds the query's reverse complement there."""

    var ref_pos: UInt32
    var reverse: Bool

    def __lt__(self, other: Self) -> Bool:
        return self.ref_pos < other.ref_pos or (self.ref_pos == other.ref_pos and not self.reverse and other.reverse)

    def __eq__(self, other: Self) -> Bool:
        return self.ref_pos == other.ref_pos and self.reverse == other.reverse


struct Smems(Movable, Sized):
    """SMEMs in increasing query start, each with its hits sorted by reference
    position (forward before reverse at the same position)."""

    var smems: List[Smem]
    var hits: List[SmemHit]

    def __init__(out self, var smems: List[Smem], var hits: List[SmemHit]):
        self.smems = smems^
        self.hits = hits^

    def __len__(self) -> Int:
        return len(self.smems)

    def hits_of(ref self, i: Int) -> Span[SmemHit, origin_of(self.hits)]:
        var smem = self.smems[i]
        return Span(self.hits)[smem.first : smem.first + smem.count]


def find_smems(index: MemIndex, reference: SeqText, query: SeqText, query_rc: SeqText, threads: Int = 1) -> Smems:
    """The SMEMs (supermaximal exact matches, as in BWA) of length >=
    `index.params.min_len` of every query sequence against the reference on
    both strands. `query_rc` packs the same sequences reverse-complemented
    (`SeqText(seqs, reverse_complement=True)`).

    An SMEM is the query interval of a MEM (on either strand) that no other
    MEM's query interval strictly contains, with every reference position
    where that whole interval matches. Any interval containing a MEM of
    length >= L is longer, so filtering at L gives exactly the SMEMs of
    length >= L."""
    var forward = find_mems(index, reference, query, threads)
    var backward = find_mems(index, reference, query_rc, threads)
    # Reverse-strand MEMs in forward query coordinates, sorted by start.
    var keys = List[UInt64](capacity=len(backward))
    var seq = 0
    for i in range(len(backward)):
        var m = backward[i]
        var q = Int(m.query_pos)
        while seq + 1 < len(query_rc.starts) and query_rc.starts[seq + 1] <= q:
            seq += 1
        var base = query_rc.starts[seq]
        var start = base + query_rc.lengths[seq] - (q - base) - Int(m.length)
        keys.append((UInt64(start) << 32) | UInt64(i))
    parallel_sort(keys, threads, low_bits=32)
    # Sweep both lists by query start: of the MEMs starting at one position,
    # only the longest can be an SMEM, and it is one unless an earlier start
    # reaches as far.
    var smems = List[Smem]()
    var hits = List[SmemHit]()
    var reach = -1
    var f = 0
    var b = 0
    while f < len(forward) or b < len(keys):
        var start: Int
        if b == len(keys) or (f < len(forward) and Int(forward[f].query_pos) <= Int(keys[b] >> 32)):
            start = Int(forward[f].query_pos)
        else:
            start = Int(keys[b] >> 32)
        var longest = 0
        var f_end = f
        while f_end < len(forward) and Int(forward[f_end].query_pos) == start:
            longest = max(longest, Int(forward[f_end].length))
            f_end += 1
        var b_end = b
        while b_end < len(keys) and Int(keys[b_end] >> 32) == start:
            longest = max(longest, Int(backward[Int(keys[b_end] & 0xFFFFFFFF)].length))
            b_end += 1
        if start + longest > reach:
            var first = len(hits)
            for i in range(f, f_end):
                if Int(forward[i].length) == longest:
                    hits.append(SmemHit(forward[i].ref_pos, False))
            for i in range(b, b_end):
                var m = backward[Int(keys[i] & 0xFFFFFFFF)]
                if Int(m.length) == longest:
                    hits.append(SmemHit(m.ref_pos, True))
            var count = len(hits) - first
            sort(Span(hits)[first:])
            smems.append(Smem(UInt32(start), UInt32(start + longest), first, count))
            reach = start + longest
        f = f_end
        b = b_end
    return Smems(smems^, hits^)


comptime _REF_ORDER_TASKS = 64
comptime _PARALLEL_MEMS = 1 << 16
"""MEMs below which MUM filtering runs on one thread."""


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


def filter_by_occurrences(
    mems: List[Mem],
    query: SeqText,
    max_in_reference: Int,
    max_in_query: Int = 0,
    min_in_reference: Int = 1,
    threads: Int = 1,
) -> List[Mem]:
    """MEMs whose string occurs at most `max_in_reference` times in the
    reference, at most `max_in_query` times in its query sequence (0: no
    limit, for either) and at least `min_in_reference` times in the
    reference: frequency-capped (and floored) MEMs, as Mumemto's `-F` caps
    multi-MEMs. With caps of 1 they are `filter_mums`'s MUMs.

    The string of a MEM occurs in the reference once per MEM whose query
    interval contains its own (itself included), and in the query once per
    MEM whose reference interval does: sweeps that count containing
    intervals, up to the cap, give both. `mems` must come from one
    `find_mems` call; the result keeps its order."""
    if min_in_reference <= 1 and max_in_reference <= 1 and max_in_query <= 1:
        if max_in_reference == 0 and max_in_query == 0:
            return mems.copy()
        if max_in_reference == 1:
            return filter_mums(mems, query, unique_in_query=max_in_query == 1, threads=threads)
    var keep = List[Bool](length=len(mems), fill=True)
    var workers = threads if len(mems) >= _PARALLEL_MEMS else 1
    var ref_cap = max(max_in_reference, min_in_reference)
    var lo = 0
    while lo < len(mems):
        var seq = query.sequence_of(Int(mems[lo].query_pos))
        var seq_end = query.starts[seq] + query.lengths[seq]
        var hi = lo
        while hi < len(mems) and Int(mems[hi].query_pos) < seq_end:
            hi += 1
        var count = hi - lo
        var starts = List[UInt64](unsafe_uninit_length=count)
        var lengths = List[UInt32](unsafe_uninit_length=count)
        var others = List[UInt32](unsafe_uninit_length=count)
        if max_in_reference > 0 or min_in_reference > 1:
            # Copies in the reference: MEMs containing it in query coordinates.
            var fill = _GroupFillTask(
                mems.unsafe_ptr().unsafe_offset(lo), starts.unsafe_ptr(), lengths.unsafe_ptr(), count, by_ref=False
            )
            parallel_for(fill, _REF_ORDER_TASKS, workers)
            _count_containing(starts, 0, lengths, max(ref_cap, 1), others, workers)
            for i in range(count):
                var copies = Int(others[i]) + 1
                if (max_in_reference > 0 and copies > max_in_reference) or copies < min_in_reference:
                    keep[lo + i] = False
        if max_in_query > 0:
            # Copies in the query: MEMs containing it in reference coordinates.
            var keys = _GroupFillTask(
                mems.unsafe_ptr().unsafe_offset(lo), starts.unsafe_ptr(), lengths.unsafe_ptr(), count, by_ref=True
            )
            parallel_for(keys, _REF_ORDER_TASKS, workers)
            parallel_sort(starts, workers, low_bits=32)
            var gather = _RefOrderTask(
                mems.unsafe_ptr().unsafe_offset(lo), starts.unsafe_ptr().as_imm(), lengths.unsafe_ptr(), count
            )
            parallel_for(gather, _REF_ORDER_TASKS, workers)
            _count_containing(starts, 32, lengths, max_in_query, others, workers)
            for i in range(count):
                if Int(others[i]) + 1 > max_in_query:
                    keep[lo + Int(starts[i] & 0xFFFFFFFF)] = False
        lo = hi
    return _collect(mems, keep, workers)


struct _TopEnds(Copyable, Movable):
    """The `cap` largest interval ends seen so far, ascending in
    `ends[head:]`."""

    var ends: List[Int]
    var head: Int
    var cap: Int

    def __init__(out self, cap: Int):
        self.ends = List[Int](capacity=min(cap, 1 << 16) + 16)
        self.head = 0
        self.cap = cap

    @inline(.always)
    def at_least(self, end: Int) -> Int:
        """How many of the kept ends reach `end`: the number of seen ends
        that do, capped at `cap`."""
        var lo = self.head
        var hi = len(self.ends)
        var p = self.ends.unsafe_ptr()
        # Most intervals are contained in none: no end reaches theirs.
        if hi == lo or p[unsafe_offset=hi - 1] < end:
            return 0
        while lo < hi:
            var mid = (lo + hi) >> 1
            if p[unsafe_offset=mid] < end:
                lo = mid + 1
            else:
                hi = mid
        return len(self.ends) - lo

    @inline(.always)
    def add(mut self, end: Int):
        if len(self.ends) - self.head == self.cap:
            if end <= self.ends[self.head]:
                return
            self.head += 1
            if self.head >= max(self.cap, 64):
                # Drop the evicted prefix now and then.
                var kept = List[Int](capacity=self.cap + 16)
                for i in range(self.head, len(self.ends)):
                    kept.append(self.ends[i])
                self.ends = kept^
                self.head = 0
        # Sweeps meet ends mostly in increasing order: insert from the back.
        self.ends.append(end)
        var i = len(self.ends) - 1
        var p = self.ends.unsafe_ptr()
        while i > self.head and p[unsafe_offset=i - 1] > end:
            p[unsafe_offset=i] = p[unsafe_offset=i - 1]
            i -= 1
        p[unsafe_offset=i] = end


def _count_containing(
    starts: List[UInt64], shift: UInt64, lengths: List[UInt32], cap: Int, mut others: List[UInt32], threads: Int
):
    """`others[i]`: how many other intervals `[starts[j] >> shift, ... +
    lengths[j])`, in increasing start, contain (or equal) interval `i`,
    capped at `cap`. Parallel by blocks as `_mark_contained`: each block's
    sweep starts from the largest ends of the blocks before it."""
    var count = len(starts)
    if threads <= 1 or count < 1 << 16:
        var top = _TopEnds(cap)
        _count_sweep(Span(starts), shift, Span(lengths), others.unsafe_ptr(), 0, count, top, cap, True)
        return
    var blocks = threads * 4
    var bounds = List[Int](length=blocks + 1, fill=count)
    bounds[0] = 0
    for b in range(1, blocks):
        var at = max(bounds[b - 1], b * count // blocks)
        while at < count and at > 0 and starts[at] >> shift == starts[at - 1] >> shift:
            at += 1
        bounds[b] = at
    # Each block's own largest ends, then what precedes each block.
    var tops = List[_TopEnds](length=blocks, fill=_TopEnds(cap))
    var own = _CountTask(
        starts.unsafe_ptr().as_imm(),
        shift,
        lengths.unsafe_ptr().as_imm(),
        others.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        bounds.unsafe_ptr().as_imm(),
        tops.unsafe_ptr(),
        count,
        cap,
        write=False,
    )
    parallel_for(own, blocks, threads)
    var running = _TopEnds(cap)
    for b in range(blocks):
        var mine = tops[b].copy()
        tops[b] = running.copy()
        for i in range(mine.head, len(mine.ends)):
            running.add(mine.ends[i])
    var sweep = _CountTask(
        starts.unsafe_ptr().as_imm(),
        shift,
        lengths.unsafe_ptr().as_imm(),
        others.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        bounds.unsafe_ptr().as_imm(),
        tops.unsafe_ptr(),
        count,
        cap,
        write=True,
    )
    parallel_for(sweep, blocks, threads)
    _ = len(tops)


def _count_sweep(
    starts: Span[UInt64, _],
    shift: UInt64,
    lengths: Span[UInt32, _],
    others: MutPointer[UInt32, _],
    lo: Int,
    hi: Int,
    mut top: _TopEnds,
    cap: Int,
    write: Bool,
):
    """Sweeps groups `[lo, hi)` (whole groups) after the ends in `top`,
    adding their ends to it; with `write`, sets `others`. Of a group sharing
    a start, every member at least as long as an interval contains it."""
    var group = List[UInt64]()
    var i = lo
    while i < hi:
        var start = Int(starts[i] >> shift)
        var j = i + 1
        while j < hi and Int(starts[j] >> shift) == start:
            j += 1
        if write:
            if j - i == 1:
                others[unsafe_offset=i] = UInt32(top.at_least(start + Int(lengths[i])))
            else:
                # Longest first: a member's copies within the group are the
                # members at least as long, ranked before it or tied with it.
                group.clear()
                for t in range(i, j):
                    group.append((UInt64(lengths[t]) << 32) | UInt64(t - i))
                sort(Span(group))
                var at = len(group)
                while at > 0:
                    var length = group[at - 1] >> 32
                    var run = at - 1
                    while run > 0 and group[run - 1] >> 32 == length:
                        run -= 1
                    # Members [run, len(group)) are at least this long.
                    var copies = min(cap, top.at_least(start + Int(length)) + len(group) - run - 1)
                    for u in range(run, at):
                        others[unsafe_offset=i + Int(group[u] & 0xFFFFFFFF)] = UInt32(copies)
                    at = run
        for t in range(i, j):
            top.add(start + Int(lengths[t]))
        i = j


struct _CountTask[
    starts_origin: ImmOrigin, lengths_origin: ImmOrigin, bounds_origin: ImmOrigin, tops_origin: MutOrigin
](ParallelTask):
    """One block of a parallel `_count_containing`: its own largest ends, or
    (`write`) its counts after the ends before it, in `tops`."""

    var starts: ImmPointer[UInt64, Self.starts_origin]
    var shift: UInt64
    var lengths: ImmPointer[UInt32, Self.lengths_origin]
    var others: MutPointer[UInt32, MutUntrackedOrigin]
    var bounds: ImmPointer[Int, Self.bounds_origin]
    var tops: MutPointer[_TopEnds, Self.tops_origin]
    var count: Int
    var cap: Int
    var write: Bool

    def __init__(
        out self,
        starts: ImmPointer[UInt64, Self.starts_origin],
        shift: UInt64,
        lengths: ImmPointer[UInt32, Self.lengths_origin],
        others: MutPointer[UInt32, MutUntrackedOrigin],
        bounds: ImmPointer[Int, Self.bounds_origin],
        tops: MutPointer[_TopEnds, Self.tops_origin],
        count: Int,
        cap: Int,
        *,
        write: Bool,
    ):
        self.starts = starts
        self.shift = shift
        self.lengths = lengths
        self.others = others
        self.bounds = bounds
        self.tops = tops
        self.count = count
        self.cap = cap
        self.write = write

    def run(self, index: Int):
        var lo = self.bounds[unsafe_offset=index]
        var hi = self.bounds[unsafe_offset=index + 1]
        var starts = Span(unsafe_ptr=self.starts, length=self.count)
        var lengths = Span(unsafe_ptr=self.lengths, length=self.count)
        var top = self.tops[unsafe_offset=index].copy()
        _count_sweep(starts, self.shift, lengths, self.others, lo, hi, top, self.cap, self.write)
        if not self.write:
            self.tops[unsafe_offset=index] = top^


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
