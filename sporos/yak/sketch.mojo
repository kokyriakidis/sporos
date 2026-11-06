"""hifiasm's minimizer sketch (`mz1_ha_sketch`, sketch.cpp) and its high-occurrence subsampling (`select_mz_h`).

A (w, k)-minimizer sketch of a homopolymer-compressed read where a minimizer is the smallest k-mer of its window
under `(count, hash)`: the count is what the filter table knows of the k-mer (0 for an ordinary k-mer, its count
for one over the cut-off), so the rarest k-mer of a window wins and a filtered k-mer (`INT32_MAX`) never does.
Positions are the k-mer's END position on the raw read (the last letter of its last homopolymer run) and `span`
how many raw letters the compressed k-mer covers. Then `select_mz_h`: in every stretch of consecutive
high-occurrence minimizers long enough (`sample_dist`), those that are minima of a wider window (`rewin`) are
kept, else the few smallest, and the rest of the high-occurrence ones are dropped.
"""

from sporos.yak.filter import FilterTable
from sporos.yak.hash import HASH_LANES, Kmer, nt4, yak_hash64_64, yak_hash64_64_lanes

comptime _DUMMY_X = UInt64.MAX
comptime _DUMMY_RID = (UInt32(1) << 28) - 1
comptime _FILTERED = 1 << 28
"""`ha_ft_cnt` answers at least this for a k-mer over `max_kmer_cnt`: the sketch drops it."""
comptime _MAX_HIGH_OCC = 8
comptime _MAX_MAX_HIGH_OCC = 16
comptime _MARK: UInt64 = 1 << 32
"""`mt` entries carry the minimizer's `l` in the low 32 bits and this mark when `select_mz_h` keeps it."""


@fieldwise_init
struct Minimizer(ImplicitlyCopyable, TrivialRegisterPassable):
    """`ha_mz1_t`: the hash `x`, `rid` (the k-mer's count while sketching, the read's id after), the END position,
    the strand (`rev` 1 when the reverse complement is canonical) and the span in raw letters."""

    var x: UInt64
    var rid: UInt32
    var pos: UInt32
    var rev: UInt8
    var span: UInt8

    @staticmethod
    @inline(.always)
    def dummy() -> Minimizer:
        return Minimizer(_DUMMY_X, _DUMMY_RID, 0, 0, 0)

    @inline(.always)
    def is_dummy(self) -> Bool:
        return self.x == _DUMMY_X

    @inline(.always)
    def start(self) -> Int:
        return Int(self.pos) + 1 - Int(self.span)


@inline(.always)
def mzcmp(a: Minimizer, b: Minimizer) -> Int:
    """`mz1_mzcmp`: by count, then by hash."""
    if a.rid < b.rid:
        return -1
    if a.rid > b.rid:
        return 1
    if a.x > b.x:
        return 1
    if a.x < b.x:
        return -1
    return 0


@inline(.always)
def mz_lt(a: Minimizer, b: Minimizer) -> Bool:
    return mzcmp(a, b) < 0


struct _TinyQueue:
    """`tiny_queue_t`: the raw lengths of the last k homopolymer runs."""

    var a: Array[Int32, 64]
    var front: Int
    var count: Int

    def __init__(out self):
        self.a = Array[Int32, 64](fill=0)
        self.front = 0
        self.count = 0

    @inline(.always)
    def push(mut self, x: Int):
        self.a[(self.count + self.front) & 0x3F] = Int32(x)
        self.count += 1

    @inline(.always)
    def shift(mut self) -> Int:
        if self.count == 0:
            return -1
        var x = Int(self.a[self.front])
        self.front = (self.front + 1) & 0x3F
        self.count -= 1
        return x

    @inline(.always)
    def reset(mut self):
        self.front = 0
        self.count = 0


struct Sketch(Movable):
    """One read's minimizers (`p`) and, for each, its `l` (letters since the last ambiguous one when it was
    buffered) with `select_mz_h`'s mark (`mt`)."""

    var mz: List[Minimizer]
    var mt: List[UInt64]
    var compressed_length: Int
    """`tl`: how many letters the compressed read has (symmetric k-mers' letters counted too)."""

    def __init__(out self):
        self.mz = List[Minimizer]()
        self.mt = List[UInt64]()
        self.compressed_length = 0

    def clear(mut self):
        self.mz.clear()
        self.mt.clear()
        self.compressed_length = 0

    @inline(.always)
    def push(mut self, m: Minimizer, l: Int):
        self.mz.append(m)
        self.mt.append(UInt64(l))


struct _Events:
    """What the serial scan of a read leaves for the hashing and the window: one event per position the window
    sees (a homopolymer run of A/C/G/T, or an ambiguous letter; symmetric k-mers are not events), with its `l`
    and, for a whole k-mer, its index among the k-mers (-1 otherwise); per k-mer, its END position, strand,
    span and canonical planes, then its hash. Sized for the read once; `events` and `kmers` say how much is
    used."""

    var l: List[Int32]
    var kmer: List[Int32]
    var pos: List[UInt32]
    var rev: List[UInt8]
    var span: List[UInt8]
    var plane0: List[UInt64]
    var plane1: List[UInt64]
    var hash: List[UInt64]
    var events: Int
    var kmers: Int

    def __init__(out self, capacity: Int):
        self.l = List[Int32](length=capacity, fill=0)
        self.kmer = List[Int32](length=capacity, fill=0)
        self.pos = List[UInt32](length=capacity, fill=0)
        self.rev = List[UInt8](length=capacity, fill=0)
        self.span = List[UInt8](length=capacity, fill=0)
        self.plane0 = List[UInt64](length=capacity, fill=0)
        self.plane1 = List[UInt64](length=capacity, fill=0)
        self.hash = List[UInt64](length=capacity, fill=0)
        self.events = 0
        self.kmers = 0


def _scan(seq: Span[UInt8, _], k: Int, is_hpc: Bool, mut ev: _Events) -> Int:
    """The serial part of `mz1_ha_sketch`: homopolymer runs, the rolling planes, `l`, the span; returns `tl`."""
    var length = len(seq)
    var src = seq.unsafe_ptr()
    var el = ev.l.unsafe_ptr()
    var ek = ev.kmer.unsafe_ptr()
    var epos = ev.pos.unsafe_ptr()
    var erev = ev.rev.unsafe_ptr()
    var espan = ev.span.unsafe_ptr()
    var e0 = ev.plane0.unsafe_ptr()
    var e1 = ev.plane1.unsafe_ptr()
    var kmer = Kmer.empty(k)
    var tq = _TinyQueue()
    var l = 0
    var tl = 0
    var kmer_span = 0
    var events = 0
    var kmers = 0
    var i = 0
    while i < length:
        var c = nt4(src[unsafe_offset=i])
        if c < 4:
            if is_hpc:
                var skip_len = 1
                if i + 1 < length and nt4(src[unsafe_offset=i + 1]) == c:
                    skip_len = 2
                    while i + skip_len < length and nt4(src[unsafe_offset=i + skip_len]) == c:
                        skip_len += 1
                    i += skip_len - 1
                tq.push(skip_len)
                kmer_span += skip_len
                if tq.count > k:
                    kmer_span -= tq.shift()
            else:
                kmer_span = l + 1 if l + 1 < k else k
            kmer.push(c)
            if kmer.symmetric():
                # `continue` in the C: a symmetric k-mer neither enters the buffer nor advances it.
                i += 1
                continue
            l += 1
            tl += 1
            if l >= k and kmer_span < 256:
                var planes = kmer.canonical()
                ek[unsafe_offset=events] = Int32(kmers)
                epos[unsafe_offset=kmers] = UInt32(i)
                erev[unsafe_offset=kmers] = UInt8(1 if kmer.reverse_is_canonical() else 0)
                espan[unsafe_offset=kmers] = UInt8(kmer_span)
                e0[unsafe_offset=kmers] = planes[0]
                e1[unsafe_offset=kmers] = planes[1]
                kmers += 1
            else:
                ek[unsafe_offset=events] = -1
        else:
            l = 0
            tq.reset()
            kmer_span = 0
            ek[unsafe_offset=events] = -1
        el[unsafe_offset=events] = Int32(l)
        events += 1
        i += 1
    ev.events = events
    ev.kmers = kmers
    return tl


def _hash_all(mut ev: _Events):
    """Every k-mer's `yak_hash_long`, eight at a time."""
    comptime LANES = HASH_LANES
    var n = ev.kmers
    var p0 = ev.plane0.unsafe_ptr()
    var p1 = ev.plane1.unsafe_ptr()
    var out = ev.hash.unsafe_ptr()
    var i = 0
    while i + LANES <= n:
        var h = yak_hash64_64_lanes(p0.unsafe_load[width=LANES](i)) + yak_hash64_64_lanes(
            p1.unsafe_load[width=LANES](i)
        )
        out.unsafe_store(i, h)
        i += LANES
    while i < n:
        out[unsafe_offset=i] = yak_hash64_64(p0[unsafe_offset=i]) + yak_hash64_64(p1[unsafe_offset=i])
        i += 1


def sketch_read(
    seq: Span[UInt8, _],
    w: Int,
    k: Int,
    is_hpc: Bool,
    filter: FilterTable,
    sample_dist: Int,
    rewin: Int,
    mut out: Sketch,
):
    """`mz1_ha_sketch` with `hf` the filter (an empty `FilterTable()` for none), `sample_dist` and `ws = rewin`, `pt` none, `rid` left
    as the k-mer counts' residue of 0 (the caller sets it). In three passes with the C's result: the serial scan,
    the hashes eight k-mers at a time, then the window over the events."""
    out.clear()
    var length = len(seq)
    if length <= 0:
        return
    var ev = _Events(length + 16)
    var tl = _scan(seq, k, is_hpc, ev)
    _hash_all(ev)
    var buf = List[Minimizer](length=w, fill=Minimizer.dummy())
    var buf_p = List[Int](length=w, fill=0)
    var buf_pos = 0
    var min_pos = 0
    var min_m = Minimizer.dummy()
    var min_s = Int(UInt32.MAX)
    for e in range(ev.events):
        var l = Int(ev.l[e])
        var info = Minimizer.dummy()
        var km = Int(ev.kmer[e])
        if km >= 0:
            var y = ev.hash[km]
            var cnt = filter.count(y)
            if cnt < _FILTERED:
                info = Minimizer(y, UInt32(cnt), ev.pos[km], ev.rev[km], ev.span[km])
        buf[buf_pos] = info
        buf_p[buf_pos] = l
        if l == w + k - 1 and not min_m.is_dummy():
            # The first window: identical k-mers were not stored yet.
            for j in range(buf_pos + 1, w):
                if mzcmp(min_m, buf[j]) == 0 and buf[j].pos != min_m.pos:
                    out.push(buf[j], buf_p[j])
            for j in range(0, buf_pos):
                if mzcmp(min_m, buf[j]) == 0 and buf[j].pos != min_m.pos:
                    out.push(buf[j], buf_p[j])
        if mzcmp(min_m, info) >= 0:
            # A new minimum: write the old one.
            if l >= w + k and not min_m.is_dummy():
                out.push(min_m, min_s)
            min_m = info
            min_pos = buf_pos
            min_s = buf_p[buf_pos]
        elif buf_pos == min_pos:
            # The old minimum left the window: write it, find the new one over the whole buffer.
            if l >= w + k - 1 and not min_m.is_dummy():
                out.push(min_m, min_s)
            min_m = Minimizer.dummy()
            for j in range(buf_pos + 1, w):
                if mzcmp(min_m, buf[j]) >= 0:
                    min_m = buf[j]
                    min_pos = j
                    min_s = buf_p[j]
            for j in range(0, buf_pos + 1):
                if mzcmp(min_m, buf[j]) >= 0:
                    min_m = buf[j]
                    min_pos = j
                    min_s = buf_p[j]
            if l >= w + k - 1 and not min_m.is_dummy():
                for j in range(buf_pos + 1, w):
                    if mzcmp(min_m, buf[j]) == 0 and min_m.pos != buf[j].pos:
                        out.push(buf[j], buf_p[j])
                for j in range(0, buf_pos + 1):
                    if mzcmp(min_m, buf[j]) == 0 and min_m.pos != buf[j].pos:
                        out.push(buf[j], buf_p[j])
        buf_pos += 1
        if buf_pos == w:
            buf_pos = 0
    if not min_m.is_dummy():
        out.push(min_m, min_s)
    out.compressed_length = tl
    if sample_dist > w:
        select_mz_h(out, length, sample_dist, rewin, k, tl)


@inline(.always)
def _gl(mt: List[UInt64], i: Int) -> Int:
    return Int(UInt32(mt[i] & 0xFFFFFFFF))


@inline(.always)
def _is_high(mz: List[Minimizer], m: Int) -> Bool:
    """`A_M`: a high-occurrence entry (its count is kept in `rid`)."""
    return m >= 0 and mz[m].rid > 0


def _mzcmp_l(mz: List[Minimizer], ai: Int, bi: Int) -> Int:
    """`mz1_mzcmp_l`: two entries by index, high-occurrence ones ordered by `mzcmp`, an ordinary one after
    them, an absent one (-1) first."""
    if ai >= 0 and bi >= 0:
        var a = mz[ai]
        var b = mz[bi]
        if a.rid > 0 and b.rid > 0:
            return mzcmp(a, b)
        return Int(a.rid == 0) - Int(b.rid == 0)
    return Int(ai < 0) - Int(bi < 0)


def _qfw(mz: List[Minimizer], mut mt: List[UInt64], n: Int, tot_l: Int, ws: Int, start: Int) -> Tuple[Int, Int]:
    """`mz1_qfw`: the first wider window's minimum among the high-occurrence entries, marked; returns the index
    it stopped at and the minimum's index (`mi`)."""
    var si = start
    var i = start
    var mi = -1
    while i < n:
        if (
            _gl(mt, i) >= ws
            or (i + 1 < n and _gl(mt, i) < ws and _gl(mt, i + 1) > ws)
            or (i + 1 == n and tot_l >= ws and _gl(mt, i) < ws)
        ):
            for m in range(si, i + 1):
                if not _is_high(mz, m):
                    continue
                if _mzcmp_l(mz, mi, m) >= 0:
                    mi = m
            if mi >= 0 and _is_high(mz, mi):
                for m in range(si, i + 1):
                    if not _is_high(mz, m):
                        continue
                    if _mzcmp_l(mz, mi, m) == 0:
                        mt[m] |= _MARK
            break
        i += 1
    return (i, mi)


def _hf_select(mut mz: List[Minimizer], si: Int, ei: Int, n: Int, length: Int, sample_dist: Int, force: Bool):
    """`mz1_hf_select`: in the high-occurrence stretch `si + 1 .. ei`, keep up to `(pe - ps) / sample_dist` of
    the smallest entries (by count, then hash) whose count is below the stretch's length."""
    if ei - si <= 1:
        return
    var ps = 0 if si < 0 else Int(mz[si].pos)
    var pe = length if ei == n else Int(mz[ei].pos)
    var st = si + 1
    var en = ei
    var max_high_occ = Int(Float64(pe - ps) / Float64(sample_dist) + 0.499)
    if max_high_occ > _MAX_MAX_HIGH_OCC:
        max_high_occ = _MAX_MAX_HIGH_OCC
    # A max-heap (by `mz_lt`) of the smallest `max_high_occ` entries, each remembering its index in `pos`.
    var heap = List[Minimizer](capacity=_MAX_MAX_HIGH_OCC)
    var j = st
    while j < en and len(heap) < max_high_occ:
        var b = mz[j]
        b.pos = UInt32(j)
        heap.append(b)
        j += 1
    _heap_make(heap)
    while j < en:
        if mz_lt(mz[j], heap[0]):
            var b = mz[j]
            b.pos = UInt32(j)
            heap[0] = b
            _heap_down(heap, 0)
        j += 1
    for b in heap:
        if Int(b.rid) < pe - ps or force:
            mz[Int(b.pos)].rid = 0


def _heap_down(mut heap: List[Minimizer], start: Int):
    """`ks_heapdown` for `mz_lt`: a max-heap, the largest at the root."""
    var n = len(heap)
    var i = start
    var tmp = heap[i]
    var k = (i << 1) + 1
    while k < n:
        if k + 1 < n and mz_lt(heap[k], heap[k + 1]):
            k += 1
        if mz_lt(heap[k], tmp):
            break
        heap[i] = heap[k]
        i = k
        k = (i << 1) + 1
    heap[i] = tmp


def _heap_make(mut heap: List[Minimizer]):
    var n = len(heap)
    var i = (n >> 1) - 1
    while i >= 0:
        _heap_down(heap, i)
        i -= 1


def select_mz_h(mut out: Sketch, length: Int, sample_dist: Int, rewin: Int, k: Int, tot_l: Int):
    """`mz1_select_mz_h`: thins the high-occurrence minimizers (those with a count, `rid > 0`) and drops the rest
    of them, keeping every ordinary one; `rewin` is the wider window (`w` inside the C function)."""
    var n = len(out.mz)
    if n == 0:
        return
    var w = rewin
    var ws = w + k - 1
    # Any stretch of high-occurrence entries long enough to subsample?
    var any_streak = False
    var last0 = -1
    for i in range(n + 1):
        if i == n or out.mz[i].rid == 0:
            if i - last0 > 1:
                var ps = 0 if last0 < 0 else Int(out.mz[last0].pos)
                var pe = length if i == n else Int(out.mz[i].pos)
                if Int(Float64(pe - ps) / Float64(sample_dist) + 0.499) > 0:
                    any_streak = True
                    break
            last0 = i
    if not any_streak:
        return
    var first = _qfw(out.mz, out.mt, n, tot_l, ws, 0)
    var i = first[0]
    var mi = first[1]
    var si = 0
    if i < n:
        i += 1
        while i < n:
            while si < i:
                if _gl(out.mt, si) + w > _gl(out.mt, i):
                    break
                si += 1
            if _mzcmp_l(out.mz, i, mi) <= 0:
                if _is_high(out.mz, mi):
                    out.mt[mi] |= _MARK
                mi = i
            elif si > mi:
                if _is_high(out.mz, mi):
                    out.mt[mi] |= _MARK
                mi = -1
                for m in range(si, i + 1):
                    if _mzcmp_l(out.mz, mi, m) >= 0:
                        mi = m
                if _is_high(out.mz, mi):
                    for m in range(si, i + 1):
                        if not _is_high(out.mz, m):
                            continue
                        if _mzcmp_l(out.mz, mi, m) == 0:
                            out.mt[m] |= _MARK
            i += 1
        if _is_high(out.mz, mi):
            out.mt[mi] |= _MARK
        i = n - 1
        while si < n and _gl(out.mt, si) + w <= tot_l + 1:
            if si > mi:
                if _is_high(out.mz, mi):
                    out.mt[mi] |= _MARK
                mi = -1
                for m in range(si, i + 1):
                    if _mzcmp_l(out.mz, mi, m) >= 0:
                        mi = m
                if _is_high(out.mz, mi):
                    for m in range(si, i + 1):
                        if not _is_high(out.mz, m):
                            continue
                        if _mzcmp_l(out.mz, mi, m) == 0:
                            out.mt[m] |= _MARK
            si += 1
        # In each long stretch, keep the marked entries; a stretch with none keeps the smallest few.
        last0 = -1
        for j in range(n + 1):
            if j == n or out.mz[j].rid == 0:
                if j - last0 > 1:
                    var ps = 0 if last0 < 0 else Int(out.mz[last0].pos)
                    var pe = length if j == n else Int(out.mz[j].pos)
                    if Int(Float64(pe - ps) / Float64(sample_dist) + 0.499) > 0:
                        var kept = 0
                        for m in range(last0 + 1, j):
                            if out.mt[m] & _MARK != 0:
                                out.mz[m].rid = 0
                                kept += 1
                        if kept == 0:
                            _hf_select(out.mz, last0, j, n, length, sample_dist, False)
                last0 = j
    # Squeeze out what is still high-occurrence.
    var mz = List[Minimizer](capacity=n)
    var mt = List[UInt64](capacity=n)
    for j in range(n):
        if out.mz[j].rid == 0:
            mz.append(out.mz[j])
            mt.append(out.mt[j])
    out.mz = mz^
    out.mt = mt^
