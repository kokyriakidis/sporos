"""Building an FM-index: the text of ACGT runs and their reverse
complements, its suffix array by libsais, the BWT blocks, samples and
separators, the two-bit copy and the k-mer cache."""

from std.bit import count_trailing_zeros
from std.ffi import OwnedDLHandle
from std.memory import Layout, pack_bits
from std.os import getenv, setenv as os_setenv
from std.sys.info import CompilationTarget, size_of
from std.sys.intrinsics import prefetch

from sporos.pages import advise_huge_pages
from sporos.parallel import ParallelTask, parallel_for
from sporos.text import _pack32

from .index import (
    FMIndex,
    Interval,
    _AHEAD,
    _BUILD_BLOCKS,
    _COUNT_MASK,
    _ENCODE_AHEAD,
    _MAX_CACHE_K,
    _SAMPLE_BYTES,
    _SEPARATOR,
    _SEPARATOR_FLAG,
    _code,
)


def build_index[
    mut: Bool, //, origin: Origin[mut=mut]
](
    seqs: List[Span[UInt8, origin]],
    names: List[String],
    threads: Int = 1,
    sample_bits: Int = 0,
    wide_suffix_array: Bool = False,
) raises -> FMIndex:
    """Indexes `seqs` (named `names`) with up to `threads` threads, with a
    suffix array sample every 2^`sample_bits` rows (0: every 4 rows up to
    2^29 bases, 8 up to 2^31, then 16). `wide_suffix_array` builds with
    64-bit suffix array values even when 32 bits suffice (for tests)."""
    var index = FMIndex()
    index.names = names.copy()
    var offset = 0
    for s in seqs:
        index.lengths.append(len(s))
        index.offsets.append(offset)
        offset += len(s) + 1
    # ACGT runs of the forward strand, then the text: runs and separators,
    # then the reverse complement of all that but the last separator, then one.
    if offset > Int(UInt32.MAX):
        raise Error("FM-index references must be shorter than 2^32 bases in all")
    var half = 0
    for s in range(len(seqs)):
        var runs = _acgt_runs(seqs[s])
        for r in range(0, len(runs), 2):
            index.run_text.append(half)
            index.run_seq.append(s)
            index.run_ref.append(runs[r])
            half += runs[r + 1] - runs[r] + 1
    index.run_text.append(half)
    index.half = half
    var n = 2 * half
    index.rows = n
    var nblocks = (n >> 7) + 1
    index._blocks = _huge_zeros[DType.uint64](nblocks * 8 + 8)
    index._first = _aligned_first(index._blocks)
    index.sample_bits = sample_bits if sample_bits > 0 else (2 if n <= 1 << 30 else (3 if n <= 1 << 32 else 4))
    index.samples = _huge_zeros[DType.uint8](index._sample_count() * _SAMPLE_BYTES + 8)
    index.packed = _huge_zeros[DType.uint64]((half >> 5) + 2)
    if n == 0:
        _build_cache(index, threads)
        return index^
    var text = List[UInt8](unsafe_uninit_length=n + 32)
    var sources = List[Int](capacity=len(index.run_seq))
    for r in range(len(index.run_seq)):
        sources.append(Int(seqs[index.run_seq[r]].unsafe_ptr()) + index.run_ref[r])
    var chunks = (half + _FILL_CHUNK - 1) // _FILL_CHUNK
    for phase in range(2):
        var fill = _FillTask(
            text.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
            index.packed.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
            index.run_text.unsafe_ptr().as_imm().unsafe_origin_cast[ImmUntrackedOrigin](),
            sources.unsafe_ptr().as_imm().unsafe_origin_cast[ImmUntrackedOrigin](),
            len(index.run_seq),
            half,
            phase,
        )
        parallel_for(fill, chunks, threads)
    _ = len(sources)
    _ = len(index.packed)

    var lib = _load_libsais()
    var separators = List[UInt64](length=(n + 63) >> 6, fill=0)
    if n < Int(Int32.MAX) and not wide_suffix_array:
        var sa = List[Int32](unsafe_uninit_length=n)
        var status: Int32
        if threads > 1:
            status = lib.get_function[Int32]("libsais_omp")(
                text.unsafe_ptr(), sa.unsafe_ptr(), Int32(n), Int32(0), Int(0), Int32(threads)
            )
        else:
            status = lib.get_function[Int32]("libsais")(text.unsafe_ptr(), sa.unsafe_ptr(), Int32(n), Int32(0), Int(0))
        if status != 0:
            raise Error("libsais failed")
        _encode(index, text, sa, separators, threads)
    else:
        var sa = List[Int64](unsafe_uninit_length=n)
        var status: Int64
        if threads > 1:
            status = lib.get_function[Int64]("libsais64_omp")(
                text.unsafe_ptr(), sa.unsafe_ptr(), Int64(n), Int64(0), Int(0), Int64(threads)
            )
        else:
            status = lib.get_function[Int64]("libsais64")(
                text.unsafe_ptr(), sa.unsafe_ptr(), Int64(n), Int64(0), Int(0)
            )
        if status != 0:
            raise Error("libsais64 failed")
        _encode(index, text, sa, separators, threads)
    _build_cache(index, threads)
    return index^


def _encode[
    dtype: DType
](mut fm: FMIndex, text: List[UInt8], sa: List[Scalar[dtype]], mut separators: List[UInt64], threads: Int):
    """Fills the blocks and samples from the suffix array, then lists the
    separator rows and turns block counts into cumulative ones."""
    var n = fm.rows
    var nblocks = (n >> 7) + 1
    var task = _EncodeTask[dtype](
        text.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
        sa.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
        fm._blocks.unsafe_ptr().unsafe_offset(fm._first).unsafe_origin_cast[MutUntrackedOrigin](),
        fm.samples.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        separators.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        n,
        nblocks,
        fm.sample_bits,
    )
    parallel_for(task, (nblocks + _BUILD_BLOCKS - 1) // _BUILD_BLOCKS, threads)
    for w in range(len(separators)):
        var bits = separators[w]
        while bits != 0:
            var b = Int(count_trailing_zeros(bits))
            bits &= bits - 1
            var row = (w << 6) + b
            fm.separators.append(row)
            fm.separator_pos.append(Int(sa[row]))
    _ = len(text)
    var blocks = fm._blocks.unsafe_ptr().unsafe_offset(fm._first)
    var running = SIMD[DType.uint64, 4](0)
    for b in range(nblocks):
        var p = blocks.unsafe_offset(b << 3)
        var word = p.unsafe_load[width=4]()
        var local = word & _COUNT_MASK
        p.unsafe_store(running | (word & ~_COUNT_MASK))
        running += local
    fm.counts = running.cast[DType.int64]()
    var first = Int64(len(fm.separators))
    fm.starts = SIMD[DType.int64, 4](
        first,
        first + fm.counts[0],
        first + fm.counts[0] + fm.counts[1],
        first + fm.counts[0] + fm.counts[1] + fm.counts[2],
    )


def _build_cache(mut fm: FMIndex, threads: Int):
    """Every k-mer's interval, k growing with the reference (none for tiny ones)."""
    var k = 0
    while k < _MAX_CACHE_K and (1 << (2 * (k + 3))) <= fm.rows:
        k += 1
    fm.cache_k = k if k >= 4 else 0
    if fm.cache_k == 0:
        fm.cache = List[SIMD[DType.uint32, 4]]()
        return
    fm.cache = _huge_cache(fm.cache_k)
    var top = min(3, fm.cache_k)
    var cache = fm.cache.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
    var task = _CacheTask(Pointer(to=fm).as_imm(), cache, top)
    parallel_for(task, 1 << (2 * top), threads)


def _acgt_runs(seq: Span[UInt8, _]) -> List[Int]:
    """`[start0, end0, start1, end1, ...]`: the maximal ACGT runs of `seq`
    (either case, U as T), 32 bases per step."""
    comptime W = 32
    var out = List[Int]()
    var n = len(seq)
    var src = seq.unsafe_ptr()
    var inside = False
    var i = 0
    while i < n:
        var take = min(W, n - i)
        var bits: UInt32
        if take == W:
            var lower = src.unsafe_load[width=W](i) | 0x20
            var valid = lower.eq(0x61) | lower.eq(0x63) | lower.eq(0x67) | lower.eq(0x74) | lower.eq(0x75)
            bits = pack_bits[DType.uint32](valid)
            if bits == (UInt32.MAX if inside else 0):
                i += W
                continue
        else:
            bits = 0
            for j in range(take):
                if _code(src[unsafe_offset=i + j]) != _SEPARATOR:
                    bits |= UInt32(1) << UInt32(j)
        for j in range(take):
            var valid = (bits >> UInt32(j)) & 1 != 0
            if valid != inside:
                out.append(i + j)
                inside = valid
        i += take
    if inside:
        out.append(n)
    return out^


comptime _FILL_CHUNK = 1 << 20
"""Forward-half text bases per task when the text is written."""


struct _FillTask(ParallelTask):
    """Writes one chunk of the text: in phase 0 the forward half (run bases
    as codes 1-4, a 0 after each run), in phase 1 the reverse-complement
    half and the two-bit copy of the forward half."""

    var text: MutPointer[UInt8, MutUntrackedOrigin]
    var packed: MutPointer[UInt64, MutUntrackedOrigin]
    var run_text: ImmPointer[Int, ImmUntrackedOrigin]
    var sources: ImmPointer[Int, ImmUntrackedOrigin]
    var runs: Int
    var half: Int
    var phase: Int

    def __init__(
        out self,
        text: MutPointer[UInt8, MutUntrackedOrigin],
        packed: MutPointer[UInt64, MutUntrackedOrigin],
        run_text: ImmPointer[Int, ImmUntrackedOrigin],
        sources: ImmPointer[Int, ImmUntrackedOrigin],
        runs: Int,
        half: Int,
        phase: Int,
    ):
        self.text = text
        self.packed = packed
        self.run_text = run_text
        self.sources = sources
        self.runs = runs
        self.half = half
        self.phase = phase

    def run(self, index: Int):
        comptime W = 32
        var lo = index * _FILL_CHUNK
        var hi = min(self.half, lo + _FILL_CHUNK)
        if self.phase == 0:
            # The run holding lo, then run by run.
            var r = 0
            var count = self.runs
            while count > 1:
                var step = count >> 1
                r += step * Int(self.run_text[unsafe_offset=r + step] <= lo)
                count -= step
            var pos = lo
            while pos < hi:
                var start = self.run_text[unsafe_offset=r]
                var end = self.run_text[unsafe_offset=r + 1] - 1  # its separator
                var stop = min(hi, end)
                var src = ImmPointer[UInt8, ImmUntrackedOrigin](
                    unsafe_from_address=self.sources[unsafe_offset=r] + pos - start
                )
                var dst = self.text.unsafe_offset(pos)
                var i = 0
                while i + W <= stop - pos:
                    var raw = src.unsafe_load[width=W](i)
                    dst.unsafe_offset(i).unsafe_store((((raw >> 1) ^ (raw >> 2)) & 3) + 1)
                    i += W
                while i < stop - pos:
                    dst[unsafe_offset=i] = _code(src[unsafe_offset=i]) + 1
                    i += 1
                pos = stop
                if pos == end and pos < hi:
                    self.text[unsafe_offset=pos] = 0
                    pos += 1
                    r += 1
            return
        # Phase 1: G[g] = complement of F[half - 2 - g] (separators stay 0),
        # G[half - 1] = 0; and the two-bit words of F (separators as 0).
        var half = self.half
        var g = lo
        while g < hi:
            if g + W <= half - 1 and half - 33 - g >= 0:
                var f = self.text.unsafe_offset(half - 33 - g).unsafe_load[width=W]().reversed()
                var comp = f.eq(0).select(SIMD[DType.uint8, W](0), SIMD[DType.uint8, W](5) - f)
                self.text.unsafe_offset(half + g).unsafe_store(comp)
                g += W
            else:
                var f = self.text[unsafe_offset=half - 2 - g] if g < half - 1 else UInt8(0)
                self.text[unsafe_offset=half + g] = 0 if f == 0 else 5 - f
                g += 1
        # Words fully inside [lo, hi); chunks are multiples of 32 bases.
        for word in range(lo >> 5, (hi + 31) >> 5):
            var base = word << 5
            var f = self.text.unsafe_offset(base).unsafe_load[width=W]()
            var codes = f.eq(0).select(SIMD[DType.uint8, W](0), f - 1)
            var value = _pack32(codes)
            if base + W > half:
                value &= (UInt64(1) << UInt64(2 * (half - base))) - 1
            self.packed[unsafe_offset=word] = value


struct _EncodeTask[dtype: DType](ParallelTask):
    """Encodes one run of BWT blocks from the suffix array: symbols, the
    blocks' own counts (made cumulative afterwards), samples and separator
    bits."""

    var text: ImmPointer[UInt8, ImmUntrackedOrigin]
    var sa: ImmPointer[Scalar[Self.dtype], ImmUntrackedOrigin]
    var blocks: MutPointer[UInt64, MutUntrackedOrigin]
    var samples: MutPointer[UInt8, MutUntrackedOrigin]
    var separators: MutPointer[UInt64, MutUntrackedOrigin]
    var n: Int
    var nblocks: Int
    var sample_bits: Int

    def __init__(
        out self,
        text: ImmPointer[UInt8, ImmUntrackedOrigin],
        sa: ImmPointer[Scalar[Self.dtype], ImmUntrackedOrigin],
        blocks: MutPointer[UInt64, MutUntrackedOrigin],
        samples: MutPointer[UInt8, MutUntrackedOrigin],
        separators: MutPointer[UInt64, MutUntrackedOrigin],
        n: Int,
        nblocks: Int,
        sample_bits: Int,
    ):
        self.text = text
        self.sa = sa
        self.blocks = blocks
        self.samples = samples
        self.separators = separators
        self.n = n
        self.nblocks = nblocks
        self.sample_bits = sample_bits

    def run(self, index: Int):
        for b in range(index * _BUILD_BLOCKS, min(self.nblocks, (index + 1) * _BUILD_BLOCKS)):
            var counts = SIMD[DType.uint64, 4](0)
            var first_half = SIMD[DType.uint64, 4](0)
            var words = SIMD[DType.uint64, 4](0)
            var flag = False
            var lo = b << 7
            var hi = min(self.n, lo + 128)
            for row in range(lo, hi):
                var s = Int(self.sa[unsafe_offset=row])
                # The text is read at random suffix array positions: fetch ahead.
                if row + _ENCODE_AHEAD < self.n:
                    var ahead = Int(self.sa[unsafe_offset=row + _ENCODE_AHEAD])
                    prefetch[_AHEAD](self.text.unsafe_offset(max(ahead - 1, 0)))
                if row & ((1 << self.sample_bits) - 1) == 0:
                    # Byte stores: a wider one would reach the next task's samples.
                    var at = (row >> self.sample_bits) * _SAMPLE_BYTES
                    for byte in range(_SAMPLE_BYTES):
                        self.samples[unsafe_offset=at + byte] = UInt8((s >> (8 * byte)) & 0xFF)
                var sym = 0 if s == 0 else Int(self.text[unsafe_offset=s - 1])
                if sym == 0:
                    flag = True
                    self.separators[unsafe_offset=row >> 6] |= UInt64(1) << UInt64(row & 63)
                else:
                    var c = sym - 1
                    words[(row & 127) >> 5] |= UInt64(c) << UInt64(2 * (row & 31))
                    counts[c] += 1
                if row & 127 == 63:
                    first_half = counts
            var head = counts | (first_half << 56)
            if flag:
                head[0] |= _SEPARATOR_FLAG
            var p = self.blocks.unsafe_offset(b << 3)
            p.unsafe_store(head)
            p.unsafe_offset(4).unsafe_store(words)


struct _CacheTask[index_origin: ImmOrigin](ParallelTask):
    """Fills the cache for the k-mers ending with one `top`-base suffix,
    depth first from that suffix's interval."""

    var index: ImmPointer[FMIndex, Self.index_origin]
    var cache: MutPointer[SIMD[DType.uint32, 4], MutUntrackedOrigin]
    var top: Int

    def __init__(
        out self,
        index: ImmPointer[FMIndex, Self.index_origin],
        cache: MutPointer[SIMD[DType.uint32, 4], MutUntrackedOrigin],
        top: Int,
    ):
        self.index = index
        self.cache = cache
        self.top = top

    def run(self, task: Int):
        ref index = self.index[]
        var k = index.cache_k
        # The suffix: base `top - 1 - t` of `task` (two bits each, last base highest).
        var key = 0
        var p = Interval(0, 0, 0)
        for t in range(self.top):
            var c = (task >> (2 * (self.top - 1 - t))) & 3
            key |= c << (2 * (k - 1 - t))
            p = index.base_interval(c) if t == 0 else index.backward(p, c)
            if p.size == 0:
                return
        var stack = List[Interval]()
        var keys = List[Int]()
        var depths = List[Int]()
        stack.append(p)
        keys.append(key)
        depths.append(self.top)
        while len(stack) > 0:
            var ik = stack.pop()
            var kk = keys.pop()
            var d = depths.pop()
            if d == k:
                var e = SIMD[DType.uint32, 4](
                    UInt32(ik.lo & 0xFFFFFFFF),
                    UInt32(ik.rc & 0xFFFFFFFF),
                    UInt32(ik.size & 0xFFFFFFFF),
                    UInt32((ik.lo >> 32) | ((ik.rc >> 32) << 8) | ((ik.size >> 32) << 16)),
                )
                self.cache[unsafe_offset=kk] = e
                continue
            for c in range(4):
                var ok = index.backward(ik, c)
                if ok.size > 0:
                    stack.append(ok)
                    keys.append(kk | (c << (2 * (k - 1 - d))))
                    depths.append(d + 1)


def _load_libsais() raises -> OwnedDLHandle:
    # OpenMP threads that spin after each libsais call compete with sporos's
    # own threads that run next: let them sleep (unless the user chose).
    if not getenv("OMP_WAIT_POLICY"):
        _ = os_setenv("OMP_WAIT_POLICY", "passive")
    comptime suffix = ".dylib" if CompilationTarget.is_macos() else ".so"
    var names: List[String] = [
        "@rpath/libsais.2.dylib" if CompilationTarget.is_macos() else "libsais.so.2",
        String("libsais") + suffix,
        "/opt/homebrew/lib/libsais.dylib",
    ]
    var prefix = getenv("CONDA_PREFIX")
    if prefix:
        names.append(prefix + "/lib/libsais" + suffix)
    for name in names:
        try:
            var lib = OwnedDLHandle(name)
            # A second handle, never closed, keeps the library loaded: its
            # OpenMP worker threads outlive each call, and unmapping their
            # code under them crashes.
            var second = OwnedDLHandle(name)
            alloc(Layout[OwnedDLHandle](length=1)).unsafe_leak().unsafe_write(second^)
            return lib^
        except:
            pass
    raise Error("building an FM-index needs libsais (conda-forge: libsais)")


def _huge_zeros[dtype: DType](length: Int) -> List[Scalar[dtype]]:
    """A zeroed list backed by transparent huge pages where available: the
    index is read at random, and huge pages spare most TLB misses."""
    var list = List[Scalar[dtype]](length=length, fill=0)
    advise_huge_pages(Int(list.unsafe_ptr()), length * size_of[Scalar[dtype]]())
    return list^


def _huge_uninit[dtype: DType](length: Int) -> List[Scalar[dtype]]:
    """`_huge_zeros` without the zeros, for a list about to be overwritten."""
    var list = List[Scalar[dtype]](unsafe_uninit_length=length)
    advise_huge_pages(Int(list.unsafe_ptr()), length * size_of[Scalar[dtype]]())
    return list^


def _huge_cache(k: Int) -> List[SIMD[DType.uint32, 4]]:
    var cache = List[SIMD[DType.uint32, 4]](length=1 << (2 * k), fill=SIMD[DType.uint32, 4](0))
    advise_huge_pages(Int(cache.unsafe_ptr()), len(cache) * 16)
    return cache^


def _aligned_first(blocks: List[UInt64]) -> Int:
    """Word offset of the first 64-byte boundary in `blocks`."""
    return ((-Int(blocks.unsafe_ptr())) & 63) >> 3
