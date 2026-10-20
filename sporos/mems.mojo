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

    def count(self, chunk: Int) -> Int:
        var first = chunk * _BUILD_CHUNK
        var last = min(self.n - self.k + 1, first + _BUILD_CHUNK)
        var mask = _kmer_mask(self.k)
        var count = 0
        for start in range(_next_multiple(first, self.step), last, self.step):
            count += Int(self.text[].mask(start) & mask == 0)
        return count

    def fill(
        self,
        chunk: Int,
        hashes: MutPointer[UInt64, MutUntrackedOrigin],
        entries: MutPointer[UInt64, MutUntrackedOrigin],
    ):
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


struct MemIndex(Movable):
    """Sampled K-mer index of a reference text.

    Buckets are the top bits of each anchor's `mix64` hash; an entry holds
    the hash's low 32 bits (a fingerprint) above the anchor's position.
    """

    var params: MemParams
    var table: HashTable

    def __init__(
        out self, reference: SeqText, params: MemParams, threads: Int = 1
    ) raises:
        params.validate()
        if len(reference) >= Int(UInt32.MAX):
            raise Error("reference text must be shorter than 2^32")
        self.params = params
        var source = _SampledAnchors(
            Pointer(to=reference).unsafe_origin_cast[ImmUntrackedOrigin](), params.k, params.ref_step
        )
        self.table = HashTable(source, threads)

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
        var ent = self.table.entries.unsafe_ptr()
        # Software pipeline over anchors: anchor n is hashed (its bucket bounds
        # prefetched), anchor n - _LAG has its bounds read (entries prefetched),
        # and anchor n - 2 * _LAG is matched.
        var positions = Array[Int, _RING](fill=0)
        var hashes = Array[UInt64, _RING](fill=0)
        var lows = Array[UInt32, _RING](fill=0)
        var highs = Array[UInt32, _RING](fill=0)
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
                self._match((n - 2 * _LAG) & (_RING - 1), positions, hashes, lows, highs, reference, query, out)
            n += 1
        for m in range(max(0, n - _LAG), n):
            self._bounds(m, hashes, lows, highs)
        for m in range(max(0, n - 2 * _LAG), n):
            self._match(m & (_RING - 1), positions, hashes, lows, highs, reference, query, out)

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
        prefetch[_AHEAD](self.table.entries.unsafe_ptr().unsafe_offset(Int(lows[slot])))

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
        mut out: List[Mem],
    ):
        var k = self.params.k
        var period = self.params.ref_step * self.params.query_step
        var min_len = self.params.min_len
        var ent = self.table.entries.unsafe_ptr()
        var j = positions[slot]
        var fp = hashes[slot] << 32
        for e in range(Int(lows[slot]), Int(highs[slot])):
            var entry = ent[unsafe_offset=e]
            if (entry ^ fp) >> 32 != 0:
                continue
            var r = Int(entry & 0xFFFFFFFF)
            # Anchors of one MEM sit `period` apart on its diagonal, so this
            # anchor is the MEM's first iff the match does not reach `period`
            # bases to the left: stop looking there.
            var left = match_backward(reference, r, query, j, period)
            if left >= period:
                continue
            var right = match_forward(reference, r, query, j)
            if right < k or left + right < min_len:
                continue
            out.append(Mem(UInt32(r - left), UInt32(j - left), UInt32(left + right)))


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
