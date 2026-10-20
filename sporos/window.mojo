"""SIMD sliding-window minimum over rolling ntHash values.

This is the shared engine behind minimizers and syncmers. A segment of window
start positions is split into `LANES` contiguous chunks and every lane runs
the same scalar algorithm in a SIMD register:

1. roll the 32-bit forward (and reverse-complement) ntHash one base,
2. key each k-mer by its top `31 - b` hash bits followed by its b-bit offset
   inside a w-aligned block (`b = offset_bits(w)`); k-mers with a non-ACGT
   base get `HASH_NONE` and lose to every valid k-mer,
3. maintain the minimum of the last `w` keys with the branch-free
   prefix/suffix-block method: 3 `min` operations per element, amortized.

Single u32 keys make every `min` one instruction while keeping exact
leftmost tie-breaking: offsets order k-mers inside a block, and the
cross-block comparison looks at hash bits only. Bases reach the
lanes through a 16-step byte transpose, and seed lookups are single byte
shuffles (NEON `tbl` / x86 `pshufb`).

Window results go to a small transposed block buffer that is compacted per
lane, then lanes are concatenated in order.
"""

from std.bit import rotate_bits_left, rotate_bits_right
from std.math import ceildiv
from std.memory import Layout, alloc, bitcast, dealloc

from sporos.alphabet import N_CODE
from sporos.nthash import SEED_A, SEED_C, SEED_G, SEED_T, rotl32, rotr32
from sporos.parallel import ParallelTask, parallel_for

comptime LANES = 8
comptime _Mask = SIMD[DType.bool, LANES]
comptime _U32V = SIMD[DType.uint32, LANES]
comptime _U64V = SIMD[DType.uint64, LANES]

comptime SENTINEL: UInt32 = 0xFFFFFFFF
"""Marks a window without output."""
comptime _BLOCK = 512
comptime PADDING = 96
"""Bytes of `N_CODE` padding kernels need after the last base."""
comptime _MIN_SEGMENT = 1 << 20

comptime MINIMIZERS = 0
comptime CLOSED_SYNCMERS = 1
comptime OPEN_SYNCMERS = 2


comptime _U8x16 = SIMD[DType.uint8, 16]
comptime _TRANSPOSE = 16
"""Steps per transposed input block (bytes loaded per lane at once)."""


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
    """Writes `_TRANSPOSE` steps of `LANES` codes, step-major, into `dst`.

    Lane `l` reads `src[l * chunk : l * chunk + 16]`; three rounds of
    interleaves turn the 8 lane vectors into 16 vectors of 8 lanes.
    """
    var v0 = src.unsafe_load[width=16](0)
    var v1 = src.unsafe_load[width=16](chunk)
    var v2 = src.unsafe_load[width=16](2 * chunk)
    var v3 = src.unsafe_load[width=16](3 * chunk)
    var v4 = src.unsafe_load[width=16](4 * chunk)
    var v5 = src.unsafe_load[width=16](5 * chunk)
    var v6 = src.unsafe_load[width=16](6 * chunk)
    var v7 = src.unsafe_load[width=16](7 * chunk)
    var a0 = bitcast[DType.uint16, 16](v0.interleave(v1))
    var a1 = bitcast[DType.uint16, 16](v2.interleave(v3))
    var a2 = bitcast[DType.uint16, 16](v4.interleave(v5))
    var a3 = bitcast[DType.uint16, 16](v6.interleave(v7))
    var b0 = bitcast[DType.uint32, 16](a0.interleave(a1))
    var b1 = bitcast[DType.uint32, 16](a2.interleave(a3))
    dst.unsafe_store(0, bitcast[DType.uint8, 128](b0.interleave(b1)))


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
    """The part of a k-mer hash that orders minimizers: its top
    `31 - offset_bits(w)` bits (ties are broken toward the leftmost k-mer)."""
    return (hash >> 1) >> UInt32(offset_bits(w))


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
    syncmers it is `j` itself when the k-mer at `j` is a syncmer.
    """
    var k = params.k
    var w = params.w
    var chunk = ceildiv(seg_windows, LANES)
    var steps = chunk + k + w - 2
    var first_out = k + w - 2

    var fin = _seed_table(SEED_A, SEED_C, SEED_G, SEED_T)
    var fout = _seed_table(
        rotl32(SEED_A, k), rotl32(SEED_C, k), rotl32(SEED_G, k), rotl32(SEED_T, k)
    )
    var rin = _seed_table(
        rotl32(SEED_T, k - 1),
        rotl32(SEED_G, k - 1),
        rotl32(SEED_C, k - 1),
        rotl32(SEED_A, k - 1),
    )
    var rout = _seed_table(
        rotr32(SEED_T, 1), rotr32(SEED_G, 1), rotr32(SEED_C, 1), rotr32(SEED_A, 1)
    )

    # Rings: outgoing hash contributions (size k) and window keys (size w).
    var fring = List[_U32V](length=k, fill=_U32V(0))
    var rring = List[_U32V](length=k, fill=_U32V(0))
    var ring = List[_U32V](length=w, fill=_U32V(HASH_NONE))
    var fring_p = fring.unsafe_ptr()
    var rring_p = rring.unsafe_ptr()
    var ring_p = ring.unsafe_ptr()

    var buf_alloc = alloc(Layout[UInt32](length=_BLOCK * LANES))
    var buf = buf_alloc.unsafe_ptr()
    var tbuf_alloc = alloc(Layout[UInt8](length=_TRANSPOSE * LANES))
    var tbuf = tbuf_alloc.unsafe_ptr()

    # Each lane compacts into its own region of `lanes_out` (capacity chunk).
    var lanes_alloc = alloc(Layout[UInt32](length=LANES * chunk + 1))
    var lanes_out = lanes_alloc.unsafe_ptr()
    var lane_counts = List[Int](length=LANES, fill=0)
    var lane_fill = List[Int](length=LANES, fill=0)
    var lane_last = List[UInt32](length=LANES, fill=SENTINEL)
    for l in range(LANES):
        lane_counts[l] = max(0, min(chunk, seg_windows - l * chunk))

    var lane_ids = _U32V(0, 1, 2, 3, 4, 5, 6, 7)
    var pos = (
        _U32V(UInt32(seg_start)) + lane_ids * UInt32(chunk) - UInt32(k - 1)
    )
    var fh = _U32V(0)
    var rh = _U32V(0)
    var run = _U32V(0)
    var prefix = _U32V(HASH_NONE)
    var bits = offset_bits(w)
    var pos_mask_v = _U32V(UInt32((1 << bits) - 1))
    var hash_mask_v = _U32V(0x7FFFFFFF) & ~pos_mask_v
    var w_v = _U32V(UInt32(w))
    var k_v = _U32V(UInt32(k))
    var span_v = _U32V(UInt32(params.span))
    var back_v = _U32V(UInt32(w - 1))
    var offset_v = _U32V(UInt32(params.offset))
    var none_v = _U32V(HASH_NONE)
    var ri = 0
    var wi = 0
    var filled = 0
    var block_base = 0

    var src = codes.unsafe_offset(seg_start)
    var t = 0
    while t < steps:
        _transpose(src.unsafe_offset(t), chunk, tbuf)
        var block_end = min(steps, t + _TRANSPOSE)
        for tt in range(block_end - t):
            var c8 = tbuf.unsafe_load[width=LANES](tt * LANES)
            var valid = c8.lt(UInt8(N_CODE))
            var c = c8.cast[DType.uint32]()
            # Byte indices 4c..4c+3 per lane; N (c=4) gets bit 7 set -> 0.
            var index = bitcast[DType.uint8, 4 * LANES](
                c * 0x04040404 + 0x03020100 + (c >> 2) * 0x80808080
            )
            run = valid.select(run + 1, _U32V(0))

            fh = rotate_bits_left[1](fh) ^ fring_p[unsafe_offset=ri] ^ _lookup(
                fin, index
            )
            fring_p.unsafe_offset(ri).unsafe_write(_lookup(fout, index))
            var h: _U32V
            comptime if canonical:
                rh = (
                    rotate_bits_right[1](rh)
                    ^ rring_p[unsafe_offset=ri]
                    ^ _lookup(rin, index)
                )
                rring_p.unsafe_offset(ri).unsafe_write(_lookup(rout, index))
                h = fh + rh
            else:
                h = fh
            ri = 0 if ri + 1 == k else ri + 1
            # Key: top (31 - b) hash bits above the k-mer's offset within its
            # w-aligned block, so `min` inside a block is leftmost on ties.
            var key = run.ge(k_v).select(
                ((h >> 1) & hash_mask_v) | UInt32(wi), none_v
            )

            # Sliding minimum over the last w keys: prefix of the current
            # block vs suffix of the previous one. Across blocks only the
            # hash bits are compared, so the older k-mer wins ties.
            if wi == 0:
                prefix = key
            else:
                prefix = min(prefix, key)
            var cur_start = pos - UInt32(wi)
            var win: _U32V
            var from_prefix: _Mask
            if wi == w - 1:
                win = prefix
                from_prefix = _Mask(fill=True)
            else:
                var suffix = ring_p[unsafe_offset=wi + 1]
                from_prefix = (prefix | pos_mask_v).lt(suffix & ~pos_mask_v)
                win = from_prefix.select(prefix, suffix)
            ring_p.unsafe_offset(wi).unsafe_write(key)
            if wi == w - 1:
                for j in reversed(range(w - 1)):
                    ring_p.unsafe_offset(j).unsafe_write(
                        min(ring_p[unsafe_offset=j], ring_p[unsafe_offset=j + 1])
                    )
                wi = 0
            else:
                wi += 1

            if t + tt >= first_out:
                var minpos = (win & pos_mask_v) + from_prefix.select(
                    cur_start, cur_start - w_v
                )
                var value: _U32V
                comptime if mode == MINIMIZERS:
                    value = win.eq(none_v).select(_U32V(SENTINEL), minpos)
                elif mode == CLOSED_SYNCMERS:
                    var start = pos - back_v
                    value = (
                        (minpos.eq(start) | minpos.eq(pos)) & run.ge(span_v)
                    ).select(start, _U32V(SENTINEL))
                else:
                    var start = pos - back_v
                    value = (
                        minpos.eq(start + offset_v) & run.ge(span_v)
                    ).select(start, _U32V(SENTINEL))
                buf.unsafe_store(filled * LANES, value)
                filled += 1
                if filled == _BLOCK:
                    _flush(
                        buf,
                        filled,
                        block_base,
                        chunk,
                        lane_counts,
                        lane_last,
                        lane_fill,
                        lanes_out,
                    )
                    block_base += filled
                    filled = 0
            pos += 1
        t = block_end

    _flush(
        buf,
        filled,
        block_base,
        chunk,
        lane_counts,
        lane_last,
        lane_fill,
        lanes_out,
    )
    dealloc(buf_alloc^)
    dealloc(tbuf_alloc^)

    for l in range(LANES):
        append_dedup(
            out,
            Span(unsafe_ptr=lanes_out.unsafe_offset(l * chunk), length=lane_fill[l]),
        )
    dealloc(lanes_alloc^)


@inline(.always)
def _flush(
    buf: ImmPointer[UInt32, _],
    filled: Int,
    block_base: Int,
    chunk: Int,
    lane_counts: List[Int],
    mut lane_last: List[UInt32],
    mut lane_fill: List[Int],
    lanes_out: MutPointer[UInt32, _],
):
    """Branch-free per-lane compaction: every value is written, and the
    cursor only advances past values that are neither `SENTINEL` nor a
    repeat of the previous window's value."""
    for l in range(LANES):
        var limit = min(filled, lane_counts[l] - block_base)
        var last = lane_last[l]
        var dst = lanes_out.unsafe_offset(l * chunk)
        var n = lane_fill[l]
        for s in range(limit):
            var v = buf[unsafe_offset=s * LANES + l]
            dst.unsafe_offset(n).unsafe_write(v)
            n += Int((v != SENTINEL) & (v != last))
            last = v
        lane_last[l] = last
        lane_fill[l] = n


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
