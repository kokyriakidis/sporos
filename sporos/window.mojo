"""SIMD sliding-window minimum over rolling ntHash values.

This is the shared engine behind minimizers and syncmers. A segment of window
start positions is split into `LANES` contiguous chunks and every lane runs
the same scalar algorithm in a SIMD register:

1. roll the 32-bit forward (and reverse-complement) ntHash one base,
2. key each k-mer by hash bits 30..b followed by its b-bit offset
   inside a w-aligned block (`b = offset_bits(w)`); k-mers with a non-ACGT
   base get `HASH_NONE` and lose to every valid k-mer,
3. maintain the minimum of the last `w` keys with the branch-free
   prefix/suffix-block method: 3 `min` operations per element, amortized.

Single u32 keys make every `min` one instruction while keeping exact
leftmost tie-breaking: offsets order k-mers inside a block, and the
cross-block comparison looks at hash bits only. Bases reach the
lanes through a 16-step byte transpose, and seed lookups are single byte
shuffles (NEON `tbl` / x86 `pshufb`).

Window results are deduplicated in the loop with one vertical compare per
step; every 8 steps the survivors' bitmask is walked with count-trailing-zeros
and each survivor is appended to its lane, then lanes are concatenated.
"""

from std.bit import count_trailing_zeros, rotate_bits_left, rotate_bits_right
from std.math import ceildiv
from std.memory import Layout, alloc, bitcast, dealloc, pack_bits
from std.math import iota
from std.sys.defines import get_defined_int
from std.sys.info import CompilationTarget

from sporos.alphabet import N_CODE
from sporos.nthash import (
    ROTATION,
    SEED_A,
    SEED_C,
    SEED_G,
    SEED_T,
    canonical_hash,
    forward_hash,
    rotl32,
    rotr32,
)
from sporos.parallel import ParallelTask, parallel_for

comptime LANES = get_defined_int[
    "SPOROS_LANES", 16 if CompilationTarget.has_avx512f() else 8
]()
"""Sequence chunks scanned side by side, one u32 SIMD lane each."""
comptime _Mask = SIMD[DType.bool, LANES]
comptime _U32V = SIMD[DType.uint32, LANES]
comptime _U64V = SIMD[DType.uint64, LANES]

comptime SENTINEL: UInt32 = 0xFFFFFFFF
"""Marks a window without output."""
comptime _Bits = DType.uint8 if LANES == 8 else DType.uint16
comptime _BATCH = 64 // LANES
"""Steps per output bitmask (one bit per lane and step)."""
comptime PADDING = 96
"""Bytes of `N_CODE` padding kernels need after the last base."""
comptime _MIN_SEGMENT = 1 << 20

comptime MINIMIZERS = 0
comptime CLOSED_SYNCMERS = 1
comptime OPEN_SYNCMERS = 2


comptime _U8x16 = SIMD[DType.uint8, 16]
comptime _TRANSPOSE = 16
"""Steps per transposed input block (bytes loaded per lane at once)."""
comptime _ROUND = 256
"""Steps transposed ahead of the main loop at a time."""


@inline(.always)
def _seed_table(s0: UInt32, s1: UInt32, s2: UInt32, s3: UInt32) -> _U8x16:
    return bitcast[DType.uint8, 16](SIMD[DType.uint32, 4](s0, s1, s2, s3))


@inline(.always)
def _lookup(table: _U8x16, index: SIMD[DType.uint8, 4 * LANES]) -> _U32V:
    """Per-lane 32-bit table lookup with one byte shuffle (NEON tbl / pshufb).
    """
    return bitcast[DType.uint32, LANES](table._dynamic_shuffle(index))


@inline(.always)
def _transpose(
    src: ImmPointer[UInt8, _], chunk: Int, dst: MutPointer[UInt8, _]
):
    """Writes the seed-table indices of `_TRANSPOSE` steps of `LANES` codes,
    step-major, into `dst` (see `_store_indices`).

    Lane `l` reads `src[l * chunk : l * chunk + 16]`; rounds of interleaves
    at doubling element widths turn the lane vectors into step vectors.
    """
    comptime assert LANES == 8 or LANES == 16, "LANES must be 8 or 16"
    var v = List[_U8x16](capacity=LANES)
    comptime for l in range(LANES):
        v.append(src.unsafe_load[width=16](l * chunk))
    var a0 = bitcast[DType.uint16, 16](v[0].interleave(v[1]))
    var a1 = bitcast[DType.uint16, 16](v[2].interleave(v[3]))
    var a2 = bitcast[DType.uint16, 16](v[4].interleave(v[5]))
    var a3 = bitcast[DType.uint16, 16](v[6].interleave(v[7]))
    var b0 = bitcast[DType.uint32, 16](a0.interleave(a1))
    var b1 = bitcast[DType.uint32, 16](a2.interleave(a3))
    comptime if LANES == 8:
        _store_indices(dst, bitcast[DType.uint8, _TRANSPOSE * LANES](b0.interleave(b1)))
    else:
        var a4 = bitcast[DType.uint16, 16](v[8].interleave(v[9]))
        var a5 = bitcast[DType.uint16, 16](v[10].interleave(v[11]))
        var a6 = bitcast[DType.uint16, 16](v[12].interleave(v[13]))
        var a7 = bitcast[DType.uint16, 16](v[14].interleave(v[15]))
        var b2 = bitcast[DType.uint32, 16](a4.interleave(a5))
        var b3 = bitcast[DType.uint32, 16](a6.interleave(a7))
        var c0 = bitcast[DType.uint64, 16](b0.interleave(b1))
        var c1 = bitcast[DType.uint64, 16](b2.interleave(b3))
        _store_indices(dst, bitcast[DType.uint8, _TRANSPOSE * LANES](c0.interleave(c1)))


@fieldwise_init
struct ScanParams(ImplicitlyCopyable, TrivialRegisterPassable):
    """Parameters of one sliding-window scan.

    `k` is the length of the hashed k-mer (the s-mer for syncmers), `w` the
    number of consecutive k-mers per window, `span` the number of bases one
    output must cover with valid bases (`k` for minimizers, the syncmer
    length for syncmers), and `offset` the open-syncmer position.
    """

    var k: Int
    var w: Int
    var span: Int
    var offset: Int


comptime HASH_NONE: UInt32 = 0xFFFFFFFF
"""Key of a k-mer containing a non-ACGT base; valid keys are below 2^31."""


def offset_bits(w: Int) -> Int:
    """Bits needed for an offset in `[0, w)`."""
    var bits = 0
    while (1 << bits) < w:
        bits += 1
    return bits


@inline(.always)
def order_key(hash: UInt32, w: Int) -> UInt32:
    """The part of a k-mer hash that orders minimizers: bits 30 down to
    `offset_bits(w)` (ties are broken toward the leftmost k-mer)."""
    return (hash & 0x7FFFFFFF) >> UInt32(offset_bits(w))


@fieldwise_init
struct _Consts(ImplicitlyCopyable, TrivialRegisterPassable):
    """Loop-invariant vectors of one scan."""

    var fin: _U8x16
    var fout: _U8x16
    var rin: _U8x16
    var rout: _U8x16
    var k: _U32V
    var span: _U32V
    var hash_mask: _U32V
    var pos_mask: _U32V
    var w: _U32V
    var back: _U32V
    var offset: _U32V


@fieldwise_init
struct _Lanes(ImplicitlyCopyable, TrivialRegisterPassable):
    """Per-lane state carried from step to step (kept in registers)."""

    var fh: _U32V
    var rh: _U32V
    var run: _U32V
    var pos: _U32V
    var block: _U32V
    var prefix: _U32V
    var prev: _U32V


comptime _TIMES4 = _U8x16(0, 4, 8, 12, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80)
"""Maps a code c to 4c, and N to 0x80, which both tbl and pshufb read as zero."""
comptime _BYTE_IN_WORD = iota[DType.uint8, 4 * _TRANSPOSE * LANES]() & 3


@inline(.always)
def _store_indices(
    dst: MutPointer[UInt8, _], codes: SIMD[DType.uint8, _TRANSPOSE * LANES]
):
    """Stores, for each code c, the byte indices `4c .. 4c+3` that `_lookup`
    needs (N maps to 0x80.., which looks up zero and sets the lane's sign)."""
    var scaled = _TIMES4._dynamic_shuffle(codes)
    var twice = scaled.interleave(scaled)
    comptime n = _TRANSPOSE * LANES
    var four = bitcast[DType.uint8, 4 * n](
        bitcast[DType.uint16, n](twice).interleave(bitcast[DType.uint16, n](twice))
    )
    dst.unsafe_store(0, four + _BYTE_IN_WORD)


@inline(.always)
def _step[
    mode: Int, canonical: Bool, emit: Bool
](
    mut lanes: _Lanes,
    consts: _Consts,
    index: SIMD[DType.uint8, 4 * LANES],
    ri: Int,
    wi: Int,
    fring: MutPointer[_U32V, _],
    rring: MutPointer[_U32V, _],
    ring: MutPointer[_U32V, _],
    obuf: MutPointer[UInt32, _],
) -> UInt32:
    """Advances every lane one base; returns the lanes' keep bits if `emit`.

    `wi` is the k-mer's offset in its w-block: the block's prefix minimum
    is in `lanes.prefix` and the previous block's suffix minima are in
    `ring[wi + 1 ..]` (`ring[w]` is a permanent `HASH_NONE`).
    """
    # An N's index bytes are >= 0x80, so its lane is negative as a u32.
    var valid = bitcast[DType.uint32, LANES](index).lt(0x80000000)
    lanes.run = valid.select(lanes.run + 1, _U32V(0))

    lanes.fh = (
        rotate_bits_left[ROTATION](lanes.fh)
        ^ fring[unsafe_offset=ri]
        ^ _lookup(consts.fin, index)
    )
    fring.unsafe_offset(ri).unsafe_write(_lookup(consts.fout, index))
    var h: _U32V
    comptime if canonical:
        lanes.rh = (
            rotate_bits_right[ROTATION](lanes.rh)
            ^ rring[unsafe_offset=ri]
            ^ _lookup(consts.rin, index)
        )
        rring.unsafe_offset(ri).unsafe_write(_lookup(consts.rout, index))
        h = lanes.fh + lanes.rh
    else:
        h = lanes.fh

    # Hash bits 30..b above the b-bit block offset: `min` inside a
    # block is leftmost on ties; across blocks only hash bits are compared,
    # so the older (left) k-mer wins ties there too.
    var key = lanes.run.ge(consts.k).select(
        (h & consts.hash_mask) | UInt32(wi), _U32V(HASH_NONE)
    )
    lanes.prefix = min(lanes.prefix, key)
    var suffix = ring[unsafe_offset=wi + 1]
    var from_prefix = (lanes.prefix | consts.pos_mask).lt(
        suffix & ~consts.pos_mask
    )
    var win = from_prefix.select(lanes.prefix, suffix)
    ring.unsafe_offset(wi).unsafe_write(key)

    var bits: UInt32 = 0
    comptime if emit:
        # `lanes.block` is the position of the current block's first k-mer.
        var minpos = (win & consts.pos_mask) + from_prefix.select(
            lanes.block, lanes.block - consts.w
        )
        var value: _U32V
        var keep: _Mask
        comptime if mode == MINIMIZERS:
            var found = win.ne(HASH_NONE)
            value = found.select(minpos, _U32V(SENTINEL))
            keep = value.ne(lanes.prev) & found
            lanes.prev = value
        else:
            value = lanes.pos - consts.back
            var hit: _Mask
            comptime if mode == CLOSED_SYNCMERS:
                hit = minpos.eq(value) | minpos.eq(lanes.pos)
            else:
                hit = minpos.eq(value + consts.offset)
            keep = hit & lanes.run.ge(consts.span)
        obuf.unsafe_store(0, value)
        bits = UInt32(pack_bits[_Bits](keep))
    comptime if mode != MINIMIZERS:
        lanes.pos += 1
    return bits


def scan_segment[
    mode: Int, canonical: Bool
](
    codes: ImmPointer[UInt8, _],
    seg_start: Int,
    seg_windows: Int,
    params: ScanParams,
    mut out: List[UInt32],
):
    """Appends the outputs of windows `[seg_start, seg_start + seg_windows)`.

    Window `j` covers k-mers `j .. j + w - 1`. For minimizers the output is
    the minimizer position (deduplicated across consecutive windows); for
    syncmers it is `j` itself when the k-mer at `j` is a syncmer. Lanes get
    `seg_windows // LANES` windows each; the few left over run scalar.
    """
    var k = params.k
    var w = params.w
    var chunk = seg_windows // LANES
    if chunk > 0:
        _scan_lanes[mode, canonical](codes, seg_start, chunk, params, out)
    _scan_scalar[mode, canonical](
        codes, seg_start + chunk * LANES, seg_windows - chunk * LANES, params, out
    )


def _scan_lanes[
    mode: Int, canonical: Bool
](
    codes: ImmPointer[UInt8, _],
    seg_start: Int,
    chunk: Int,
    params: ScanParams,
    mut out: List[UInt32],
):
    var k = params.k
    var w = params.w
    var steps = chunk + k + w - 2
    var first_out = k + w - 2
    var bits_w = offset_bits(w)
    var pos_mask = _U32V(UInt32((1 << bits_w) - 1))
    var consts = _Consts(
        fin=_seed_table(SEED_A, SEED_C, SEED_G, SEED_T),
        fout=_seed_table(
            rotl32(SEED_A, ROTATION * k),
            rotl32(SEED_C, ROTATION * k),
            rotl32(SEED_G, ROTATION * k),
            rotl32(SEED_T, ROTATION * k),
        ),
        rin=_seed_table(
            rotl32(SEED_T, ROTATION * (k - 1)),
            rotl32(SEED_G, ROTATION * (k - 1)),
            rotl32(SEED_C, ROTATION * (k - 1)),
            rotl32(SEED_A, ROTATION * (k - 1)),
        ),
        rout=_seed_table(
            rotr32(SEED_T, ROTATION),
            rotr32(SEED_G, ROTATION),
            rotr32(SEED_C, ROTATION),
            rotr32(SEED_A, ROTATION),
        ),
        k=_U32V(UInt32(k)),
        span=_U32V(UInt32(params.span)),
        hash_mask=_U32V(0x7FFFFFFF) & ~pos_mask,
        pos_mask=pos_mask,
        w=_U32V(UInt32(w)),
        back=_U32V(UInt32(w - 1)),
        offset=_U32V(UInt32(params.offset)),
    )
    var lane_ids = iota[DType.uint32, LANES]()
    var first_kmer = _U32V(UInt32(seg_start)) + lane_ids * UInt32(chunk) - UInt32(k - 1)
    var lanes = _Lanes(
        fh=_U32V(0),
        rh=_U32V(0),
        run=_U32V(0),
        pos=first_kmer,
        block=first_kmer,
        prefix=_U32V(HASH_NONE),
        prev=_U32V(SENTINEL),
    )

    # Rings: outgoing hash contributions (k) and window keys (w + sentinel).
    var fring = List[_U32V](length=k, fill=_U32V(0))
    var rring = List[_U32V](length=k, fill=_U32V(0))
    var ring = List[_U32V](length=w + 1, fill=_U32V(HASH_NONE))
    var fring_p = fring.unsafe_ptr()
    var rring_p = rring.unsafe_ptr()
    var ring_p = ring.unsafe_ptr()

    var tbuf_alloc = alloc(Layout[UInt8](length=_ROUND * 4 * LANES))
    var tbuf = tbuf_alloc.unsafe_ptr()
    var obuf_alloc = alloc(Layout[UInt32](length=_BATCH * LANES))
    var obuf = obuf_alloc.unsafe_ptr()
    # Lane l compacts into lanes_out[l * chunk ...] (capacity chunk).
    var lanes_alloc = alloc(Layout[UInt32](length=LANES * chunk))
    var lanes_out = lanes_alloc.unsafe_ptr()
    var fill_alloc = alloc(Layout[Int](length=LANES))
    var lane_fill = fill_alloc.unsafe_ptr()
    for l in range(LANES):
        lane_fill.unsafe_offset(l).unsafe_write(0)
    var batch_bits: UInt64 = 0
    var batch = 0

    var src = codes.unsafe_offset(seg_start)
    var ri = 0
    var wi = 0
    var t = 0
    while t < steps:
        var round_len = min(_ROUND, steps - t)
        for r in range(0, round_len, _TRANSPOSE):
            _transpose(src.unsafe_offset(t + r), chunk, tbuf.unsafe_offset(r * 4 * LANES))
        var i = 0
        while i < round_len:
            var n = min(w - wi, round_len - i)
            if t + i >= first_out:
                for _ in range(n):
                    var index = tbuf.unsafe_load[width=4 * LANES](i * 4 * LANES)
                    var bits = _step[mode, canonical, True](
                        lanes,
                        consts,
                        index,
                        ri,
                        wi,
                        fring_p,
                        rring_p,
                        ring_p,
                        obuf.unsafe_offset(batch * LANES),
                    )
                    batch_bits |= UInt64(bits) << UInt64(batch * LANES)
                    batch += 1
                    if batch == _BATCH:
                        _emit(batch_bits, obuf, chunk, lane_fill, lanes_out)
                        batch_bits = 0
                        batch = 0
                    ri = 0 if ri + 1 == k else ri + 1
                    wi += 1
                    i += 1
            else:
                n = min(n, first_out - (t + i))
                for _ in range(n):
                    var index = tbuf.unsafe_load[width=4 * LANES](i * 4 * LANES)
                    _ = _step[mode, canonical, False](
                        lanes, consts, index, ri, wi, fring_p, rring_p, ring_p, obuf
                    )
                    ri = 0 if ri + 1 == k else ri + 1
                    wi += 1
                    i += 1
            if wi == w:
                # Close the block: its keys become suffix minima in place.
                for j in reversed(range(w - 1)):
                    ring_p.unsafe_offset(j).unsafe_write(
                        min(ring_p[unsafe_offset=j], ring_p[unsafe_offset=j + 1])
                    )
                lanes.prefix = _U32V(HASH_NONE)
                lanes.block += UInt32(w)
                wi = 0
        t += round_len

    _emit(batch_bits, obuf, chunk, lane_fill, lanes_out)
    for l in range(LANES):
        append_dedup(
            out,
            Span(
                unsafe_ptr=lanes_out.unsafe_offset(l * chunk),
                length=lane_fill[unsafe_offset=l],
            ),
        )
    dealloc(tbuf_alloc^)
    dealloc(obuf_alloc^)
    dealloc(lanes_alloc^)
    dealloc(fill_alloc^)


def _scan_scalar[
    mode: Int, canonical: Bool
](
    codes: ImmPointer[UInt8, _],
    first: Int,
    windows: Int,
    params: ScanParams,
    mut out: List[UInt32],
):
    """Windows `[first, first + windows)` one by one (the lanes' remainder)."""
    var k = params.k
    var w = params.w
    var tail = List[UInt32]()
    for j in range(first, first + windows):
        var best = UInt64.MAX
        for i in range(j, j + w):
            var ok = True
            for b in range(i, i + k):
                if codes[unsafe_offset=b] >= N_CODE:
                    ok = False
            if not ok:
                continue
            var p = codes.unsafe_offset(i)
            var h = canonical_hash(p, k) if canonical else forward_hash(p, k)
            best = min(best, (UInt64(order_key(h, w)) << 32) | UInt64(i))
        var minpos = Int(best & 0xFFFFFFFF)
        comptime if mode == MINIMIZERS:
            if best != UInt64.MAX:
                tail.append(UInt32(minpos))
        else:
            var span_ok = True
            for b in range(j, j + params.span):
                if codes[unsafe_offset=b] >= N_CODE:
                    span_ok = False
            var hit: Bool
            comptime if mode == CLOSED_SYNCMERS:
                hit = minpos == j or minpos == j + w - 1
            else:
                hit = minpos == j + params.offset
            if span_ok and hit:
                tail.append(UInt32(j))
    var deduped = List[UInt32]()
    for v in tail:
        if len(deduped) == 0 or deduped[len(deduped) - 1] != v:
            deduped.append(v)
    append_dedup(out, Span(deduped))


@inline(.always)
def _emit(
    bits: UInt64,
    values: ImmPointer[UInt32, _],
    chunk: Int,
    lane_fill: MutPointer[Int, _],
    lanes_out: MutPointer[UInt32, _],
):
    """Appends the values whose bit is set to their lanes' regions.

    Bit `s * LANES + l` stands for step `s` of lane `l`; walking bits in
    increasing order keeps each lane's values in step order."""
    var remaining = bits
    while remaining != 0:
        var i = Int(count_trailing_zeros(remaining))
        var lane = i % LANES
        var n = lane_fill[unsafe_offset=lane]
        lanes_out.unsafe_offset(lane * chunk + n).unsafe_write(
            values[unsafe_offset=i]
        )
        lane_fill.unsafe_offset(lane).unsafe_write(n + 1)
        remaining &= remaining - 1


def append_dedup(mut out: List[UInt32], src: Span[UInt32, _]):
    """Appends `src` to `out`, dropping a leading duplicate of `out`'s tail."""
    var begin = 0
    if len(out) > 0 and len(src) > 0 and src[0] == out[len(out) - 1]:
        begin = 1
    out.extend(
        Span(unsafe_ptr=src.unsafe_ptr().unsafe_offset(begin), length=len(src) - begin)
    )


struct _ScanTask[mode: Int, canonical: Bool](ParallelTask):
    var codes: ImmPointer[UInt8, ImmUntrackedOrigin]
    var results: MutPointer[List[UInt32], MutUntrackedOrigin]
    var segment: Int
    var windows: Int
    var params: ScanParams

    def __init__(
        out self,
        codes: ImmPointer[UInt8, ImmUntrackedOrigin],
        results: MutPointer[List[UInt32], MutUntrackedOrigin],
        segment: Int,
        windows: Int,
        params: ScanParams,
    ):
        self.codes = codes
        self.results = results
        self.segment = segment
        self.windows = windows
        self.params = params

    def run(self, index: Int):
        var start = index * self.segment
        var count = min(self.segment, self.windows - start)
        scan_segment[Self.mode, Self.canonical](
            self.codes, start, count, self.params, self.results[unsafe_offset=index]
        )


def scan[
    mode: Int, canonical: Bool
](
    codes: ImmPointer[UInt8, _],
    windows: Int,
    params: ScanParams,
    threads: Int,
) -> List[UInt32]:
    """Scans windows `[0, windows)`, split into segments across `threads`."""
    var out = List[UInt32]()
    if windows <= 0:
        return out^
    var num_segments = 1
    if threads > 1 and windows >= 2 * _MIN_SEGMENT:
        num_segments = min(threads * 4, windows // _MIN_SEGMENT)
    var segment = ceildiv(windows, num_segments)
    num_segments = ceildiv(windows, segment)
    if num_segments == 1:
        scan_segment[mode, canonical](codes, 0, windows, params, out)
        return out^
    var results = List[List[UInt32]](length=num_segments, fill=List[UInt32]())
    var task = _ScanTask[mode, canonical](
        codes.unsafe_origin_cast[ImmUntrackedOrigin](),
        results.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        segment,
        windows,
        params,
    )
    parallel_for(task, num_segments, threads)
    var total = 0
    for r in results:
        total += len(r)
    out.reserve(total)
    for r in results:
        append_dedup(out, Span(r))
    return out^
