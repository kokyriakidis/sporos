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
def match_forward(a: ImmPointer[UInt8, _], b: ImmPointer[UInt8, _]) -> Int:
    """Length of the common prefix of `a` and `b` (stops at a mismatch)."""
    var n = 0
    while True:
        var d = a.unsafe_offset(n).unsafe_bitcast[UInt64]()[] ^ b.unsafe_offset(
            n
        ).unsafe_bitcast[UInt64]()[]
        if d != 0:
            return n + Int(count_trailing_zeros(d) >> 3)
        n += 8


@inline(.always)
def match_backward(
    a: ImmPointer[UInt8, _], b: ImmPointer[UInt8, _], limit: Int
) -> Int:
    """Length of the common suffix of the text ending just before `a`, `b`,
    or some value >= `limit` once it reaches `limit`."""
    var n = 0
    while n < limit:
        var d = a.unsafe_offset(-n - 8).unsafe_bitcast[UInt64]()[] ^ b.unsafe_offset(
            -n - 8
        ).unsafe_bitcast[UInt64]()[]
        if d != 0:
            return n + Int(count_leading_zeros(d) >> 3)
        n += 8
    return n


struct MemIndex(Movable):
    """Sampled K-mer index of a reference text."""

    var params: MemParams
    var offsets: List[UInt32]
    var entries: List[UInt64]
    var shift: UInt64

    def __init__(out self, reference: SeqText, params: MemParams) raises:
        params.validate()
        if len(reference) >= Int(UInt32.MAX):
            raise Error("reference text must be shorter than 2^32")
        self.params = params
        var k = params.k
        var step = params.ref_step
        var bits = 1
        while (1 << bits) < max(1, len(reference) // step):
            bits += 1
        self.shift = UInt64(64 - bits)
        var buckets = 1 << bits
        self.offsets = List[UInt32](length=buckets + 1, fill=0)
        self.entries = List[UInt64]()

        var text = reference.ptr()
        var mask = UInt64.MAX if k == 32 else (UInt64(1) << UInt64(2 * k)) - 1
        var off = self.offsets.unsafe_ptr()

        # Pass 1: count entries per bucket.
        var kmer: UInt64 = 0
        var run = 0
        for i in range(len(reference)):
            var c = text[unsafe_offset=i]
            kmer = ((kmer << 2) | UInt64(c & 3)) & mask
            run = run + 1 if c < 4 else 0
            var start = i - k + 1
            if run >= k and start % step == 0:
                var b = Int(mix64(kmer) >> self.shift)
                off[unsafe_offset=b + 1] += 1
        for b in range(buckets):
            off[unsafe_offset=b + 1] += off[unsafe_offset=b]
        self.entries.resize(Int(off[unsafe_offset=buckets]), 0)

        # Pass 2: fill; offsets[b] advances to the end of bucket b.
        var ent = self.entries.unsafe_ptr()
        kmer = 0
        run = 0
        for i in range(len(reference)):
            var c = text[unsafe_offset=i]
            kmer = ((kmer << 2) | UInt64(c & 3)) & mask
            run = run + 1 if c < 4 else 0
            var start = i - k + 1
            if run >= k and start % step == 0:
                var h = mix64(kmer)
                var b = Int(h >> self.shift)
                var slot = Int(off[unsafe_offset=b])
                off[unsafe_offset=b] += 1
                ent.unsafe_offset(slot).unsafe_write(
                    (h << 32) | UInt64(start)
                )
        for b in reversed(range(buckets)):
            off[unsafe_offset=b + 1] = off[unsafe_offset=b]
        off[unsafe_offset=0] = 0

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
        var period = self.params.ref_step * qstep
        var min_len = self.params.min_len
        var mask = UInt64.MAX if k == 32 else (UInt64(1) << UInt64(2 * k)) - 1
        var rtext = reference.ptr()
        var qtext = query.ptr()
        var off = self.offsets.unsafe_ptr()
        var ent = self.entries.unsafe_ptr()
        var first = begin + (qstep - begin % qstep) % qstep
        var kmer: UInt64 = 0
        var run = 0
        # Warm up the rolling K-mer so the first anchor ends at first + k - 1.
        var i = max(0, first - k + 1) if first < end else end
        while i < end + k - 1 and i < len(query):
            var c = qtext[unsafe_offset=i]
            kmer = ((kmer << 2) | UInt64(c & 3)) & mask
            run = run + 1 if c < 4 else 0
            var j = i - k + 1
            i += 1
            if j < first or run < k or j % qstep != 0:
                continue
            var h = mix64(kmer)
            var b = Int(h >> self.shift)
            var fp = h << 32
            for e in range(Int(off[unsafe_offset=b]), Int(off[unsafe_offset=b + 1])):
                var entry = ent[unsafe_offset=e]
                if (entry ^ fp) >> 32 != 0:
                    continue
                var r = Int(entry & 0xFFFFFFFF)
                # Anchors of one MEM sit `period` apart on its diagonal, so
                # this anchor is the MEM's first iff the match does not reach
                # `period` bases to the left: stop looking there.
                var left = match_backward(
                    rtext.unsafe_offset(r), qtext.unsafe_offset(j), period
                )
                if left >= period:
                    continue
                var right = match_forward(
                    rtext.unsafe_offset(r), qtext.unsafe_offset(j)
                )
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
        self.index[].find(
            self.reference[], self.query[], begin, end, self.results[unsafe_offset=index]
        )


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
    sort(Span(out))
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
