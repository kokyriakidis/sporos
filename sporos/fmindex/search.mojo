"""The SMEM search: Gagie's algorithm, many queries stepped round-robin as
state machines, each step prefetching what the next one reads; unique
matches are located and extended by comparing bases (see `index.mojo`)."""

from std.sys.intrinsics import prefetch

from sporos.mems import Smem, SmemHit, Smems
from sporos.parallel import ParallelTask, parallel_for

from .index import FMIndex, Interval, _AHEAD, _BATCH, _SEPARATOR, encode_bases


@fieldwise_init
struct _Search(ImplicitlyCopyable, TrivialRegisterPassable):
    """One query's SMEM search, stepped one BWT access at a time."""

    var query: Int
    var begin: Int
    var length: Int
    var stage: Int
    var x: Int
    var i: Int
    var key: Int
    var p: Interval
    var row: Int
    """While locating: the current row; then the text offset of the match."""
    var steps: Int

    @staticmethod
    def start(query: Int, bounds: List[Int]) -> Self:
        return Self(query, bounds[query], bounds[query + 1] - bounds[query], _SETUP, 0, 0, 0, Interval(0, 0, 0), 0, 0)


comptime _SETUP = 0
comptime _CACHE_FIRST = 1
comptime _BACK_FIRST = 2
comptime _FORWARD = 3
comptime _CACHE_SECOND = 4
comptime _BACK_SECOND = 5
comptime _LOCATE = 6
comptime _DIRECT = 7
comptime _DIRECT_BACK = 8


@fieldwise_init
struct _Found(ImplicitlyCopyable, TrivialRegisterPassable):
    var query: Int
    var start: Int
    var end: Int
    var lo: Int
    var size: Int
    var pos: Int
    """Text offset of the only hit, when the search located it; else -1."""


def search_smems[
    mut: Bool, //, origin: Origin[mut=mut]
](
    fm: FMIndex,
    queries: List[Span[UInt8, origin]],
    min_len: Int,
    max_hits: Int = 20,
    min_occurrences: Int = 1,
    max_occurrences: Int = 0,
    threads: Int = 1,
) -> Smems:
    """The SMEMs of length >= `min_len` of each query against the
    reference, both strands: query intervals of maximal exact matches
    not contained in another's. Hits are listed for SMEMs with at most
    `max_hits` of them.

    With `min_occurrences` f > 1 they are f-SMEMs (Gagie's f-MEMs,
    minibwa's `-s`): the maximal query intervals that occur at least f
    times. `max_occurrences` (0: no limit) drops those occurring more
    often, as Mumemto's `-F` caps matches."""
    var bounds = List[Int](capacity=len(queries) + 1)
    bounds.append(0)
    for q in queries:
        bounds.append(bounds[len(bounds) - 1] + len(q))
    var codes = List[UInt8](unsafe_uninit_length=bounds[len(bounds) - 1] + 32)
    for i in range(len(queries)):
        encode_bases(queries[i], codes.unsafe_ptr().unsafe_offset(bounds[i]))
    return search_encoded(fm, codes, bounds, min_len, max_hits, min_occurrences, max_occurrences, threads)


def search_encoded(
    fm: FMIndex,
    codes: List[UInt8],
    bounds: List[Int],
    min_len: Int,
    max_hits: Int,
    min_occurrences: Int = 1,
    max_occurrences: Int = 0,
    threads: Int = 1,
) -> Smems:
    """`smems` on queries already encoded with `encode_bases`: query `q`
    is `codes[bounds[q] : bounds[q + 1]]`. The search reads up to 32 bytes
    past a query, so `codes` is copied when it does not hold them. With
    `threads`, slices of queries are searched in parallel."""
    if len(codes) < bounds[len(bounds) - 1] + 32:
        var padded = codes.copy()
        padded.resize(bounds[len(bounds) - 1] + 32, UInt8(4))
        return search_encoded(fm, padded, bounds, min_len, max_hits, min_occurrences, max_occurrences, threads)
    if threads > 1 and len(bounds) - 1 >= 2 * _SLICE_QUERIES:
        return _search_parallel(fm, codes, bounds, min_len, max_hits, min_occurrences, max_occurrences, threads)
    var found = List[_Found]()
    _search(fm, codes.unsafe_ptr(), bounds, max(min_len, 1), max(min_occurrences, 1), max_hits >= 1, found)
    if max_occurrences > 0:
        var kept = 0
        for i in range(len(found)):
            if found[i].size <= max_occurrences:
                found[kept] = found[i]
                kept += 1
        found.resize(unsafe_uninit_length=kept)
    var nq = len(bounds) - 1
    var result = Smems()
    # Each query's SMEMs come out in increasing start; group them by query.
    var per = List[Int](length=nq + 1, fill=0)
    for f in found:
        per[f.query + 1] += 1
    for q in range(nq):
        per[q + 1] += per[q]
    result.offsets = per.copy()
    var order = List[_Found](unsafe_uninit_length=len(found))
    for f in found:
        order[per[f.query]] = f
        per[f.query] += 1
    var rows = List[Int]()
    for f in order:
        if f.size <= max_hits and f.pos < 0:
            for r in range(f.lo, f.lo + f.size):
                rows.append(r)
    fm.locate(rows)
    result.smems = List[Smem](capacity=len(order))
    result.hits = List[SmemHit](capacity=len(rows))
    var next = 0
    for f in order:
        var first = len(result.hits)
        if f.pos >= 0:
            result.hits.append(fm.reference_hit(f.pos, f.end - f.start))
        elif f.size <= max_hits:
            for _ in range(f.size):
                result.hits.append(fm.reference_hit(rows[next], f.end - f.start))
                next += 1
            _sort_hits(result.hits, first)
        result.smems.append(Smem(UInt32(f.start), UInt32(f.end), first, f.size))
    return result^


comptime _SLICE_QUERIES = 1 << 11
"""Queries per task when many are searched with several threads."""


def _search_parallel(
    fm: FMIndex,
    codes: List[UInt8],
    bounds: List[Int],
    min_len: Int,
    max_hits: Int,
    min_occurrences: Int,
    max_occurrences: Int,
    threads: Int,
) -> Smems:
    """`search_encoded` over slices of queries in parallel, results joined
    in query order."""
    var nq = len(bounds) - 1
    var slices = (nq + _SLICE_QUERIES - 1) // _SLICE_QUERIES
    var parts = List[Smems](capacity=slices)
    for _ in range(slices):
        parts.append(Smems())
    var task = _SliceTask(
        Pointer(to=fm).as_imm(),
        Pointer(to=codes).as_imm(),
        Pointer(to=bounds).as_imm(),
        parts.unsafe_ptr(),
        nq,
        min_len,
        max_hits,
        min_occurrences,
        max_occurrences,
    )
    parallel_for(task, slices, threads)
    var out = Smems()
    for p in range(slices):
        ref part = parts[p]
        var base_smems = len(out.smems)
        var base_hits = len(out.hits)
        for smem in part.smems:
            out.smems.append(Smem(smem.start, smem.end, smem.first + base_hits, smem.count))
        out.hits.extend(Span(part.hits))
        for q in range(1, len(part.offsets)):
            out.offsets.append(base_smems + part.offsets[q])
    return out^


struct _SliceTask[fm_origin: ImmOrigin, codes_origin: ImmOrigin, bounds_origin: ImmOrigin, out_origin: MutOrigin](
    ParallelTask
):
    """Searches one slice of `_SLICE_QUERIES` queries into its own `Smems`."""

    var fm: ImmPointer[FMIndex, Self.fm_origin]
    var codes: ImmPointer[List[UInt8], Self.codes_origin]
    var bounds: ImmPointer[List[Int], Self.bounds_origin]
    var parts: MutPointer[Smems, Self.out_origin]
    var queries: Int
    var min_len: Int
    var max_hits: Int
    var min_occurrences: Int
    var max_occurrences: Int

    def __init__(
        out self,
        fm: ImmPointer[FMIndex, Self.fm_origin],
        codes: ImmPointer[List[UInt8], Self.codes_origin],
        bounds: ImmPointer[List[Int], Self.bounds_origin],
        parts: MutPointer[Smems, Self.out_origin],
        queries: Int,
        min_len: Int,
        max_hits: Int,
        min_occurrences: Int,
        max_occurrences: Int,
    ):
        self.fm = fm
        self.codes = codes
        self.bounds = bounds
        self.parts = parts
        self.queries = queries
        self.min_len = min_len
        self.max_hits = max_hits
        self.min_occurrences = min_occurrences
        self.max_occurrences = max_occurrences

    def run(self, index: Int):
        var lo = index * _SLICE_QUERIES
        var hi = min(self.queries, lo + _SLICE_QUERIES)
        # Bounds stay absolute offsets into the shared codes.
        var sub = List[Int](capacity=hi - lo + 1)
        for q in range(lo, hi + 1):
            sub.append(self.bounds[][q])
        self.parts[unsafe_offset=index] = search_encoded(
            self.fm[], self.codes[], sub, self.min_len, self.max_hits, self.min_occurrences, self.max_occurrences
        )


def _search(
    fm: FMIndex,
    codes: ImmPointer[UInt8, _],
    bounds: List[Int],
    min_len: Int,
    min_occ: Int,
    locate: Bool,
    mut found: List[_Found],
):
    """Runs the searches of every query, `_BATCH` at a time round-robin.
    With `locate`, unique matches are located as soon as they are found
    and extended by comparing bases."""
    var nq = len(bounds) - 1
    var active = List[_Search](capacity=_BATCH)
    var next = 0
    while next < nq and len(active) < _BATCH:
        active.append(_Search.start(next, bounds))
        next += 1
    var j = 0
    while len(active) > 0:
        if j >= len(active):
            j = 0
        if _step(fm, active.unsafe_ptr()[unsafe_offset=j], codes, min_len, min_occ, locate, found):
            if next < nq:
                active[j] = _Search.start(next, bounds)
                next += 1
            else:
                active[j] = active[len(active) - 1]
                _ = active.pop()
                continue
        j += 1


@inline(.always)
def _step(
    fm: FMIndex,
    mut s: _Search,
    codes: ImmPointer[UInt8, _],
    min_len: Int,
    min_occ: Int,
    locate: Bool,
    mut found: List[_Found],
) -> Bool:
    """One step of a search; True once the query is done."""
    var q = codes.unsafe_offset(s.begin)
    if s.stage == _FORWARD:
        if s.p.size == 1 and locate:
            # A unique match: find where it is, then extend it there.
            s.row = s.p.lo
            s.steps = 0
            s.key = _DIRECT
            s.stage = _LOCATE
            fm._prefetch_row(s.row)
            return False
        if s.i == s.length:
            found.append(_Found(s.query, s.x, s.i, s.p.lo, s.p.size, -1))
            return True
        var c = Int(q[unsafe_offset=s.i])
        if c < 4:
            var ok = fm.forward(s.p, c)
            if ok.size >= min_occ:
                s.p = ok
                s.i += 1
                fm._prefetch_interval(ok.rc, ok.size)
                return False
        found.append(_Found(s.query, s.x, s.i, s.p.lo, s.p.size, -1))
        _next_start(fm, s, q)
        return False
    if s.stage == _LOCATE:
        var row = s.row
        var pos = -1
        if fm._sampled(row):
            pos = fm._sample(row)
        else:
            var next = fm._lf(row)
            if next < 0:
                pos = fm.separator_pos[fm._separator_index(row)]
            else:
                s.row = next
                s.steps += 1
                fm._prefetch_row(next)
                return False
        s.row = pos + s.steps
        s.stage = s.key
        # Prefetch the first reference base the comparison reads.
        var f: Int
        if s.stage == _DIRECT:
            var m = s.i - s.x
            f = s.row + m if s.row < fm.half else 2 * fm.half - 2 - s.row - m
        else:
            f = s.row - 1 if s.row < fm.half else 2 * fm.half - 1 - s.row
        prefetch[_AHEAD](fm.packed.unsafe_ptr().unsafe_offset(max(f, 0) >> 5))
        return False
    if s.stage == _DIRECT_BACK:
        # q[i + 1, x + min_len) occurs once, at text offset `row`: compare
        # backwards down to x; then the forward extension by comparison.
        var t = s.row
        var m = s.x + min_len - 1 - s.i
        var want = s.i - s.x + 1
        var n: Int
        if t < fm.half:
            var first = fm.run_text.unsafe_ptr()[unsafe_offset=fm._run_of(t)]
            n = fm._match_left(q, s.i, want, s.begin, t - 1, first, False)
            t -= n
        else:
            var f = 2 * fm.half - 1 - t - m
            var end = fm.run_text.unsafe_ptr()[unsafe_offset=fm._run_of(f) + 1] - 1
            n = fm._match_left(q, s.i, want, s.begin, f + m, end, True)
            t -= n
        s.i -= n
        if s.i >= s.x:
            s.x = s.i + 1
            s.stage = _SETUP
        else:
            s.row = t
            s.i = s.x + min_len
            s.stage = _DIRECT
        return False
    if s.stage == _DIRECT:
        # q[x, i) occurs once, at text offset `row`: compare on from there.
        var t = s.row
        var m = s.i - s.x
        if t < fm.half:
            var end = fm.run_text.unsafe_ptr()[unsafe_offset=fm._run_of(t) + 1] - 1
            s.i += fm._match_right(q, s.i, s.length - s.i, t + m, end, False)
        else:
            # On the reverse strand: the match is the reverse complement
            # of forward bases [f, f + m), which extends leftwards.
            var f = 2 * fm.half - 1 - t - m
            var first = fm.run_text.unsafe_ptr()[unsafe_offset=fm._run_of(f)]
            s.i += fm._match_right(q, s.i, s.length - s.i, f - 1, first, True)
        found.append(_Found(s.query, s.x, s.i, 0, 1, t))
        if s.i == s.length:
            return True
        _next_start(fm, s, q)
        return False
    if s.stage == _BACK_FIRST or s.stage == _BACK_SECOND:
        if s.p.size < min_occ:
            # The base at i + 1 occurs too rarely.
            s.x = s.i + 2
            s.stage = _SETUP
            return False
        if s.stage == _BACK_FIRST and s.p.size == 1 and locate:
            # Unique already: locate it now (the forward extension needs
            # it) and finish this pass by comparing bases.
            s.row = s.p.lo
            s.steps = 0
            s.key = _DIRECT_BACK
            s.stage = _LOCATE
            fm._prefetch_row(s.row)
            return False
        var floor = s.x if s.stage == _BACK_FIRST else s.x + 1
        if s.i < floor:
            if s.stage == _BACK_FIRST:
                fm._prefetch_interval(s.p.rc, s.p.size)
                s.i = s.x + min_len
                s.stage = _FORWARD
            else:
                s.x = s.i + 1
                s.stage = _SETUP
            return False
        var ok = fm.backward(s.p, Int(q[unsafe_offset=s.i]))
        if ok.size < min_occ:
            s.x = s.i + 1
            s.stage = _SETUP
        else:
            s.p = ok
            s.i -= 1
            fm._prefetch_interval(ok.lo, ok.size)
        return False
    if s.stage == _CACHE_FIRST or s.stage == _CACHE_SECOND:
        s.p = fm._cached(s.key)
        if s.p.size < min_occ:
            # The k-mer occurs too rarely: back off to base-by-base search,
            # which finds where it fails.
            s.i += fm.cache_k
            s.p = fm.base_interval(Int(q[unsafe_offset=s.i]))
            s.i -= 1
        fm._prefetch_interval(s.p.lo, s.p.size)
        s.stage += 1
        return False
    # _SETUP: a backward search over [x, x + min_len), unless an N is in it.
    if s.length - s.x < min_len:
        return True
    var last = -1
    for t in range(s.x, s.x + min_len):
        if q[unsafe_offset=t] > 3:
            last = t
    if last >= 0:
        s.x = last + 1
        return False
    s.i = s.x + min_len - 1
    if fm.cache_k > 0 and min_len >= fm.cache_k:
        s.key = fm._key(q, s.i)
        s.i -= fm.cache_k
        fm._prefetch_cache(s.key)
        s.stage = _CACHE_FIRST
    else:
        s.p = fm.base_interval(Int(q[unsafe_offset=s.i]))
        s.i -= 1
        fm._prefetch_interval(s.p.lo, s.p.size)
        s.stage = _BACK_FIRST
    return False


@inline(.always)
def _next_start(fm: FMIndex, mut s: _Search, q: ImmPointer[UInt8, _]):
    """After an SMEM ending at `s.i` (inside the query): the backward pass
    from `s.i` that finds the next start."""
    var c = Int(q[unsafe_offset=s.i])
    if c >= 4:
        s.x = s.i + 1
        s.stage = _SETUP
    elif fm.cache_k > 0 and s.i - s.x - 1 >= fm.cache_k:
        s.key = fm._key(q, s.i)
        s.i -= fm.cache_k
        fm._prefetch_cache(s.key)
        s.stage = _CACHE_SECOND
    else:
        s.p = fm.base_interval(c)
        s.i -= 1
        fm._prefetch_interval(s.p.lo, s.p.size)
        s.stage = _BACK_SECOND


def _sort_hits(mut hits: List[SmemHit], first: Int):
    for i in range(first + 1, len(hits)):
        var h = hits[i]
        var j = i
        while j > first and h < hits[j - 1]:
            hits[j] = hits[j - 1]
            j -= 1
        hits[j] = h
