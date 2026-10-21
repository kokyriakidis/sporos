"""Maximal exact matches (MEMs) between a reference and a query text.

Sampled k-mer anchoring (as in copMEM): reference K-mers are indexed every
`k1` positions and query K-mers probed every `k2` positions, with `k1`, `k2`
coprime and `k1 * k2 <= L - K + 1`. By the Chinese remainder theorem every
match of length >= L contains an anchor pair, and the anchor offsets inside a
MEM form a single residue class modulo `k1 * k2`. A MEM is therefore reported
exactly once, from the anchor whose left extension is shorter than `k1 * k2`,
with no hashing or sorting of results.
"""

from std.bit import count_leading_zeros, count_trailing_zeros
from std.math import ceildiv, gcd, iota
from std.memory import pack_bits
from std.sys.intrinsics import PrefetchOptions, prefetch

from sporos.parallel import ParallelTask, parallel_for
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
