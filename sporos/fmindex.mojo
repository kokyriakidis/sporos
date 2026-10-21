"""FM-index SMEM engine for many short queries (reads) against one reference.

The index is the BWT of the reference's ACGT runs followed by their reverse
complements, each run ended by a separator, so that one index serves both
strands and a pattern's interval comes with its reverse complement's
(a bidirectional FMD-index, as in BWA, ropebwt3 and minibwa). Separators
split the text at sequence ends and at every non-ACGT base, so matches never
cross them: the SMEMs are exact with respect to the reference as written.

Layout (minibwa's): the BWT is cut into 128-row blocks of 64 bytes, four
cumulative base counts (low 56 bits, plus the count in the block's first 64
rows above them) and 128 two-bit symbols. Separators are stored as A and
listed apart; bit 63 of a block's A count flags the rare blocks that hold
one, where ranks subtract them. Ranks count all four bases at once with
three popcounts per word.

SMEMs follow Gagie's algorithm (minibwa's `mb_bwt_smem`): from a start `x`,
a backward search over `[x, x + L)`, a forward extension to the longest
match, then a backward pass from its end for the next start. Many queries
are searched together as state machines stepped round-robin, each step
prefetching the blocks the next one reads, so the cache misses of different
queries overlap. The first `k` bases of each backward search come from a
cache of every k-mer's interval. Once a match is unique, its one occurrence
is located (suffix array samples, five bytes each, every 8 rows, or 16 for
references over 2^32 rows) and the match is extended by comparing bases
against a two-bit copy of the reference instead of base by base in the
index. Hits of repeated SMEMs are located in batches.

The suffix array is built with libsais (Apache-2.0, conda-forge), loaded at
run time.
"""

from std.bit import count_leading_zeros, count_trailing_zeros, pop_count
from std.ffi import OwnedDLHandle
from std.memory import Layout, bitcast, pack_bits, unsafe_memcpy
from std.os import getenv
from std.sys.info import CompilationTarget, size_of
from std.sys.intrinsics import PrefetchOptions, prefetch

from sporos.mems import Smem, SmemHit
from sporos.pages import advise_huge_pages
from sporos.parallel import ParallelTask, parallel_for
from sporos.text import _pack32

comptime _COUNT_MASK: UInt64 = (UInt64(1) << 56) - 1
comptime _SEPARATOR_FLAG: UInt64 = UInt64(1) << 63
comptime _EVEN: UInt64 = 0x5555555555555555
comptime _ODD: UInt64 = 0xAAAAAAAAAAAAAAAA
comptime _SAMPLE_BYTES = 5
"""Suffix array samples are 40-bit values."""
comptime _MAX_CACHE_K = 12
comptime _SEPARATOR = 4
"""The BWT symbol of a separator, and the query code of a base that never matches."""
comptime _AHEAD = PrefetchOptions().for_read().high_locality()
comptime _BATCH = 64
"""Queries searched together (state machines stepped round-robin)."""
comptime _LOCATE_BATCH = 64
comptime _BUILD_BLOCKS = 1 << 13
"""BWT blocks per task when the index is built."""
comptime _MAGIC = "SPOROSFM"
comptime _VERSION = 3


@fieldwise_init
struct Interval(ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """A bidirectional interval: the pattern's suffixes are BWT rows
    `[lo, lo + size)`, its reverse complement's are `[rc, rc + size)`."""

    var lo: Int
    var rc: Int
    var size: Int


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


struct ReadSmems(Movable, Sized):
    """SMEMs of a batch of queries: query `q`'s are `smems[offsets[q] :
    offsets[q + 1]]`, in increasing start. Each SMEM's hits, sorted by
    reference position (forward before reverse), are listed when it has at
    most `max_hits` of them and left out otherwise; `Smem.count` is always
    the number of hits."""

    var offsets: List[Int]
    var smems: List[Smem]
    var hits: List[SmemHit]

    def __init__(out self):
        self.offsets = [0]
        self.smems = List[Smem]()
        self.hits = List[SmemHit]()

    def __len__(self) -> Int:
        return len(self.offsets) - 1

    def hits_of(ref self, i: Int) -> Span[SmemHit, origin_of(self.hits)]:
        """The listed hits of SMEM `i` (empty when it has more than `max_hits`)."""
        var end = self.smems[i + 1].first if i + 1 < len(self.smems) else len(self.hits)
        return Span(self.hits)[self.smems[i].first : end]


struct FMIndex(Movable):
    """The FMD-index of a reference (see the module notes). Build it with
    `FMIndex.build` or read one written by `save` with `FMIndex.load`."""

    var _blocks: List[UInt64]
    var _first: Int
    """Word offset of block 0 in `_blocks` (blocks are cache-line aligned)."""
    var rows: Int
    """BWT length: twice the ACGT bases plus one separator per run and strand."""
    var starts: SIMD[DType.int64, 4]
    """First row of the suffixes starting with A, C, G, T."""
    var counts: SIMD[DType.int64, 4]
    var separators: List[Int]
    """Rows whose BWT symbol is a separator, increasing."""
    var separator_pos: List[Int]
    """Their suffix array values: the text offsets of run starts."""
    var samples: List[UInt8]
    """Suffix array values of the rows divisible by 2^sample_bits, 40 bits each."""
    var sample_bits: Int
    var half: Int
    """Text offset of the reverse-complement half."""
    var run_text: List[Int]
    """Text offset of each run (forward half), with a final `half` sentinel."""
    var run_seq: List[Int]
    var run_ref: List[Int]
    """Sequence and offset in it of each run's first base."""
    var names: List[String]
    var lengths: List[Int]
    var offsets: List[Int]
    """Reference text offset of each sequence, as `SeqText.starts`."""
    var packed: List[UInt64]
    """The forward half of the text, two bits per base (separators as A),
    for extending unique matches by direct comparison."""
    var cache: List[SIMD[DType.uint32, 4]]
    var cache_k: Int

    def __init__(out self):
        self._blocks = List[UInt64]()
        self._first = 0
        self.rows = 0
        self.starts = SIMD[DType.int64, 4](0)
        self.counts = SIMD[DType.int64, 4](0)
        self.separators = List[Int]()
        self.separator_pos = List[Int]()
        self.samples = List[UInt8]()
        self.sample_bits = 3
        self.half = 0
        self.run_text = List[Int]()
        self.run_seq = List[Int]()
        self.run_ref = List[Int]()
        self.names = List[String]()
        self.lengths = List[Int]()
        self.offsets = List[Int]()
        self.packed = List[UInt64]()
        self.cache = List[SIMD[DType.uint32, 4]]()
        self.cache_k = 0

    # ---------------------------------------------------------------- build

    @staticmethod
    def build[
        mut: Bool, //, origin: Origin[mut=mut]
    ](seqs: List[Span[UInt8, origin]], names: List[String], threads: Int = 1, sample_bits: Int = 0) raises -> FMIndex:
        """Indexes `seqs` (named `names`) with up to `threads` threads, with a
        suffix array sample every 2^`sample_bits` rows (0: 8 rows, or 16 for
        references of over 2^31 bases)."""
        var index = FMIndex()
        index.names = names.copy()
        var offset = 0
        for s in seqs:
            index.lengths.append(len(s))
            index.offsets.append(offset)
            offset += len(s) + 1
        # ACGT runs of the forward strand, then the text: runs and separators,
        # then the reverse complement of all that but the last separator, then one.
        var half = 0
        for s in range(len(seqs)):
            var seq = seqs[s]
            var i = 0
            while i < len(seq):
                if _code(seq[i]) == _SEPARATOR:
                    i += 1
                    continue
                var j = i
                while j < len(seq) and _code(seq[j]) != _SEPARATOR:
                    j += 1
                index.run_text.append(half)
                index.run_seq.append(s)
                index.run_ref.append(i)
                half += j - i + 1
                i = j
        index.run_text.append(half)
        index.half = half
        var n = 2 * half
        index.rows = n
        var nblocks = (n >> 7) + 1
        index._blocks = _huge_zeros[DType.uint64](nblocks * 8 + 8)
        index._first = _aligned_first(index._blocks)
        index.sample_bits = sample_bits if sample_bits > 0 else (3 if n <= 1 << 32 else 4)
        index.samples = _huge_zeros[DType.uint8](index._sample_count() * _SAMPLE_BYTES + 8)
        if n == 0:
            index._build_cache(threads)
            return index^
        var text = List[UInt8](unsafe_uninit_length=n)
        var at = 0
        for r in range(len(index.run_seq)):
            var seq = seqs[index.run_seq[r]]
            var length = index.run_text[r + 1] - index.run_text[r] - 1
            for b in range(length):
                text[at + b] = _code(seq[index.run_ref[r] + b]) + 1
            at += length
            text[at] = 0
            at += 1
        for g in range(half - 1):
            var f = text[half - 2 - g]
            text[half + g] = 0 if f == 0 else 5 - f
        text[n - 1] = 0
        index.packed = _huge_zeros[DType.uint64]((half >> 5) + 2)
        for i in range(half):
            var code = text[i]
            if code != 0:
                index.packed[i >> 5] |= UInt64(code - 1) << UInt64(2 * (i & 31))

        var lib = _load_libsais()
        var separators = List[UInt64](length=(n + 63) >> 6, fill=0)
        if n < Int(Int32.MAX):
            var sa = List[Int32](unsafe_uninit_length=n)
            var status: Int32
            if threads > 1:
                status = lib.get_function[Int32]("libsais_omp")(
                    text.unsafe_ptr(), sa.unsafe_ptr(), Int32(n), Int32(0), Int(0), Int32(threads)
                )
            else:
                status = lib.get_function[Int32]("libsais")(
                    text.unsafe_ptr(), sa.unsafe_ptr(), Int32(n), Int32(0), Int(0)
                )
            if status != 0:
                raise Error("libsais failed")
            index._encode(text, sa, separators, threads)
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
            index._encode(text, sa, separators, threads)
        index._build_cache(threads)
        return index^

    def _encode[
        dtype: DType
    ](mut self, text: List[UInt8], sa: List[Scalar[dtype]], mut separators: List[UInt64], threads: Int):
        """Fills the blocks and samples from the suffix array, then lists the
        separator rows and turns block counts into cumulative ones."""
        var n = self.rows
        var nblocks = (n >> 7) + 1
        var task = _EncodeTask[dtype](
            text.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
            sa.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
            self._blocks.unsafe_ptr().unsafe_offset(self._first).unsafe_origin_cast[MutUntrackedOrigin](),
            self.samples.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
            separators.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
            n,
            nblocks,
            self.sample_bits,
        )
        parallel_for(task, (nblocks + _BUILD_BLOCKS - 1) // _BUILD_BLOCKS, threads)
        for w in range(len(separators)):
            var bits = separators[w]
            while bits != 0:
                var b = Int(count_trailing_zeros(bits))
                bits &= bits - 1
                var row = (w << 6) + b
                self.separators.append(row)
                self.separator_pos.append(Int(sa[row]))
        _ = len(text)
        var blocks = self._blocks.unsafe_ptr().unsafe_offset(self._first)
        var running = SIMD[DType.uint64, 4](0)
        for b in range(nblocks):
            var p = blocks.unsafe_offset(b << 3)
            var word = p.unsafe_load[width=4]()
            var local = word & _COUNT_MASK
            p.unsafe_store(running | (word & ~_COUNT_MASK))
            running += local
        self.counts = running.cast[DType.int64]()
        var first = Int64(len(self.separators))
        self.starts = SIMD[DType.int64, 4](
            first,
            first + self.counts[0],
            first + self.counts[0] + self.counts[1],
            first + self.counts[0] + self.counts[1] + self.counts[2],
        )

    def _build_cache(mut self, threads: Int):
        """Every k-mer's interval, k growing with the reference (none for tiny ones)."""
        var k = 0
        while k < _MAX_CACHE_K and (1 << (2 * (k + 3))) <= self.rows:
            k += 1
        self.cache_k = k if k >= 4 else 0
        if self.cache_k == 0:
            self.cache = List[SIMD[DType.uint32, 4]]()
            return
        self.cache = _huge_cache(self.cache_k)
        var top = min(3, self.cache_k)
        var cache = self.cache.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        var task = _CacheTask(Pointer(to=self).as_imm(), cache, top)
        parallel_for(task, 1 << (2 * top), threads)

    # ---------------------------------------------------------------- ranks

    @inline(.always)
    def _block_ptr(self) -> ImmPointer[UInt64, origin_of(self._blocks)]:
        return self._blocks.unsafe_ptr().as_imm().unsafe_offset(self._first)

    @inline(.always)
    def _prefetch(self, row: Int):
        prefetch[_AHEAD](self._block_ptr().unsafe_offset((row >> 7) << 3))

    @inline(.always)
    def _prefetch_cache(self, key: Int):
        prefetch[_AHEAD](self.cache.unsafe_ptr().unsafe_offset(key).unsafe_bitcast[UInt32]())

    @inline(.always)
    def _prefetch_interval(self, lo: Int, size: Int):
        self._prefetch(lo)
        if size > 1:
            self._prefetch(lo + size)

    @inline(.always)
    def rank4(self, k: Int) -> SIMD[DType.int64, 4]:
        """Occurrences of A, C, G, T in BWT rows `[0, k)`."""
        var p = self._block_ptr().unsafe_offset((k >> 7) << 3)
        var word = p.unsafe_load[width=4]()
        var off = k & 127
        var high = off >= 64
        var r = off & 63
        var base = word & _COUNT_MASK
        if high:
            base += (word >> 56) & 0x7F
        var w = 6 if high else 4
        var full = r >= 32
        var a = p[unsafe_offset=w] & (UInt64.MAX if full else (UInt64(1) << UInt64(2 * r)) - 1)
        var b = p[unsafe_offset=w + 1] & (((UInt64(1) << UInt64(2 * (r & 31))) - 1) if full else UInt64(0))
        var hi = pop_count(a & _ODD) + pop_count(b & _ODD)
        var lo = pop_count(a & _EVEN) + pop_count(b & _EVEN)
        var both = pop_count(a & (a >> 1) & _EVEN) + pop_count(b & (b >> 1) & _EVEN)
        var result = (base + SIMD[DType.uint64, 4](UInt64(r) - hi - lo + both, lo - both, hi - both, both)).cast[
            DType.int64
        ]()
        if word[0] & _SEPARATOR_FLAG != 0:
            result[0] -= Int64(self._separators_in((k >> 6) << 6, k))
        return result

    @inline(.always)
    def _rank1(self, p: ImmPointer[UInt64, _], k: Int, c: Int) -> Int:
        """Occurrences of base `c` in rows `[0, k)`; `p` is k's block."""
        var off = k & 127
        var high = off >= 64
        var r = off & 63
        var count = p[unsafe_offset=c]
        var base = count & _COUNT_MASK
        if high:
            base += (count >> 56) & 0x7F
        var w = 6 if high else 4
        var full = r >= 32
        var pattern = UInt64(c) * _EVEN
        var x = p[unsafe_offset=w] ^ pattern
        var y = p[unsafe_offset=w + 1] ^ pattern
        var a = ~(x | (x >> 1)) & _EVEN & (UInt64.MAX if full else (UInt64(1) << UInt64(2 * r)) - 1)
        var b = ~(y | (y >> 1)) & _EVEN & (((UInt64(1) << UInt64(2 * (r & 31))) - 1) if full else UInt64(0))
        var n = Int(base + pop_count(a) + pop_count(b))
        if c == 0 and p[] & _SEPARATOR_FLAG != 0:
            n -= self._separators_in((k >> 6) << 6, k)
        return n

    @inline(.always)
    def _lf(self, k: Int) -> Int:
        """The row of the suffix one base longer than row `k`'s (-1 when that
        base is a separator)."""
        var p = self._block_ptr().unsafe_offset((k >> 7) << 3)
        var c = Int((p[unsafe_offset=4 + ((k & 127) >> 5)] >> UInt64(2 * (k & 31))) & 3)
        if c == 0 and p[] & _SEPARATOR_FLAG != 0 and self._separator_index(k) >= 0:
            return -1
        return Int(self.starts[c]) + self._rank1(p, k, c)

    @inline(.always)
    def symbol(self, k: Int) -> Int:
        """BWT symbol of row `k`: 0-3 for ACGT, 4 for a separator."""
        var p = self._block_ptr().unsafe_offset((k >> 7) << 3)
        var c = Int((p[unsafe_offset=4 + ((k & 127) >> 5)] >> UInt64(2 * (k & 31))) & 3)
        if c == 0 and p[] & _SEPARATOR_FLAG != 0 and self._separator_index(k) >= 0:
            return _SEPARATOR
        return c

    def _separators_in(self, lo: Int, hi: Int) -> Int:
        var i = self._lower_bound(lo)
        var count = 0
        while i + count < len(self.separators) and self.separators[i + count] < hi:
            count += 1
        return count

    def _separator_index(self, row: Int) -> Int:
        var i = self._lower_bound(row)
        return i if i < len(self.separators) and self.separators[i] == row else -1

    def _lower_bound(self, row: Int) -> Int:
        var lo = 0
        var hi = len(self.separators)
        while lo < hi:
            var mid = (lo + hi) >> 1
            if self.separators[mid] < row:
                lo = mid + 1
            else:
                hi = mid
        return lo

    # ----------------------------------------------------------- intervals

    @inline(.always)
    def base_interval(self, c: Int) -> Interval:
        """The interval of the one-base pattern `c` (0-3)."""
        return Interval(Int(self.starts[c]), Int(self.starts[3 - c]), Int(self.counts[c]))

    @inline(.always)
    def backward(self, ik: Interval, c: Int) -> Interval:
        """The interval of `c` followed by the pattern of `ik`."""
        if ik.size == 1:
            # One suffix: it extends only by its own BWT symbol, and the
            # reverse complement's interval keeps its start.
            var p = self._block_ptr().unsafe_offset((ik.lo >> 7) << 3)
            var k = ik.lo
            var sym = Int((p[unsafe_offset=4 + ((k & 127) >> 5)] >> UInt64(2 * (k & 31))) & 3)
            if sym != c or (c == 0 and p[] & _SEPARATOR_FLAG != 0 and self._separator_index(k) >= 0):
                return Interval(0, 0, 0)
            return Interval(Int(self.starts[c]) + self._rank1(p, k, c), ik.rc, 1)
        var tk = self.rank4(ik.lo)
        var tl = self.rank4(ik.lo + ik.size)
        var sizes = tl - tk
        # Suffixes of the reverse complement are ordered by their next base:
        # separator first, then T, G, C, A (complements of the bases before).
        var ends = Int64(ik.size) - sizes.reduce_add()
        var lanes = SIMD[DType.int64, 4](0, 1, 2, 3)
        var after = lanes.gt(Int64(c)).select(sizes, SIMD[DType.int64, 4](0)).reduce_add()
        return Interval(Int(self.starts[c] + tk[c]), ik.rc + Int(ends + after), Int(sizes[c]))

    @inline(.always)
    def forward(self, ik: Interval, c: Int) -> Interval:
        """The interval of the pattern of `ik` followed by `c`."""
        var r = self.backward(Interval(ik.rc, ik.lo, ik.size), 3 - c)
        return Interval(r.rc, r.lo, r.size)

    @inline(.always)
    def _cached(self, key: Int) -> Interval:
        var e = self.cache.unsafe_ptr()[unsafe_offset=key]
        var high = e[3]
        return Interval(
            Int(e[0]) | Int(high & 0xFF) << 32,
            Int(e[1]) | Int((high >> 8) & 0xFF) << 32,
            Int(e[2]) | Int((high >> 16) & 0xFF) << 32,
        )

    @inline(.always)
    def _key(self, q: ImmPointer[UInt8, _], end: Int) -> Int:
        """Cache key of the k-mer ending at `end` (its last base highest)."""
        var key = 0
        for m in range(self.cache_k):
            key = (key << 2) | Int(q[unsafe_offset=end - m])
        return key

    @inline(.always)
    def _base(self, f: Int) -> Int:
        """Base `f` of the forward half of the text (0 for a separator)."""
        return Int((self.packed.unsafe_ptr()[unsafe_offset=f >> 5] >> UInt64(2 * (f & 31))) & 3)

    @inline(.always)
    def _window(self, e: Int) -> UInt64:
        """Forward-half bases `e .. e + 31` (base `e` lowest), for `e >= 0`."""
        var w = self.packed.unsafe_ptr()
        var at = e >> 5
        var shift = UInt64(2 * (e & 31))
        # (hi << 1) << (63 - shift) avoids an undefined shift by 64 when shift = 0.
        return (w[unsafe_offset=at] >> shift) | ((w[unsafe_offset=at + 1] << 1) << (63 - shift))

    def _match_right(self, q: ImmPointer[UInt8, _], i: Int, limit: Int, e: Int, bound: Int, reverse: Bool) -> Int:
        """How many query bases from `q[i]` rightwards (at most `limit`; the
        query must be readable 32 bytes past them) equal the forward half's
        from base `e` rightwards, before `bound`; with `reverse`, the
        complements of its bases from `e` leftwards, down to `bound`."""
        var room = min(limit, (e - bound + 1) if reverse else (bound - e))
        var n = 0
        while n < room:
            var raw = q.unsafe_offset(i + n).unsafe_load[width=32]()
            var stop = Int(count_trailing_zeros(pack_bits[DType.uint32](raw.gt(3))))
            var same: Int
            if not reverse:
                same = Int(count_trailing_zeros(_pack32(raw & 3) ^ self._window(e + n))) >> 1
            elif e - n >= 31:
                var rc = (SIMD[DType.uint8, 32](3) - raw).reversed() & 3
                same = Int(count_leading_zeros(_pack32(rc) ^ self._window(e - n - 31))) >> 1
            else:
                while n < room and Int(q[unsafe_offset=i + n]) == 3 - self._base(e - n):
                    n += 1
                return n
            var good = min(same, stop)
            n += good
            if good < 32:
                break
        return min(n, room)

    def _match_left(
        self, q: ImmPointer[UInt8, _], i: Int, limit: Int, begin: Int, e: Int, bound: Int, reverse: Bool
    ) -> Int:
        """How many query bases from `q[i]` leftwards (at most `limit`, all
        ACGT; `q[-begin]` is the first readable byte) equal the forward
        half's from base `e` leftwards, down to `bound`; with `reverse`, the
        complements of its bases from `e` rightwards, before `bound`."""
        var room = min(limit, (bound - e) if reverse else (e - bound + 1))
        var n = 0
        while n < room:
            if i - n - 31 < -begin or (not reverse and e - n < 31):
                if reverse:
                    while n < room and Int(q[unsafe_offset=i - n]) == 3 - self._base(e + n):
                        n += 1
                else:
                    while n < room and Int(q[unsafe_offset=i - n]) == self._base(e - n):
                        n += 1
                return n
            var raw = q.unsafe_offset(i - n - 31).unsafe_load[width=32]()
            var same: Int
            if not reverse:
                same = Int(count_leading_zeros(_pack32(raw & 3) ^ self._window(e - n - 31))) >> 1
            else:
                var rc = (SIMD[DType.uint8, 32](3) - raw).reversed() & 3
                same = Int(count_trailing_zeros(_pack32(rc) ^ self._window(e + n))) >> 1
            n += same
            if same < 32:
                break
        return min(n, room)

    @inline(.always)
    def _run_of(self, f: Int) -> Int:
        """The run holding forward-half text offset `f`."""
        var runs = self.run_text.unsafe_ptr()
        var lo = 0
        var hi = len(self.run_seq) - 1
        while lo < hi:
            var mid = (lo + hi + 1) >> 1
            if runs[unsafe_offset=mid] <= f:
                lo = mid
            else:
                hi = mid - 1
        return lo

    # ------------------------------------------------------------- locate

    def locate(self, mut rows: List[Int]):
        """Replaces each BWT row by its suffix array value (a text offset),
        walking LF from many rows at once."""
        var cur = List[Int](length=_LOCATE_BATCH, fill=0)
        var steps = List[Int](length=_LOCATE_BATCH, fill=0)
        var which = List[Int](length=_LOCATE_BATCH, fill=0)
        for first in range(0, len(rows), _LOCATE_BATCH):
            var live = min(_LOCATE_BATCH, len(rows) - first)
            for t in range(live):
                cur[t] = rows[first + t]
                steps[t] = 0
                which[t] = first + t
                self._prefetch(cur[t])
            while live > 0:
                var t = 0
                while t < live:
                    var row = cur[t]
                    var value = -1
                    if self._sampled(row):
                        value = self._sample(row)
                    else:
                        var next = self._lf(row)
                        if next < 0:
                            value = self.separator_pos[self._separator_index(row)]
                        else:
                            row = next
                    if value >= 0:
                        rows[which[t]] = value + steps[t]
                        live -= 1
                        cur[t] = cur[live]
                        steps[t] = steps[live]
                        which[t] = which[live]
                        continue
                    cur[t] = row
                    steps[t] += 1
                    self._prefetch_row(row)
                    t += 1

    @inline(.always)
    def reference_hit(self, text_pos: Int, length: Int) -> SmemHit:
        """The reference position (`SeqText` offset of the leftmost base) and
        strand of a match of `length` bases at text offset `text_pos`."""
        var reverse = text_pos >= self.half
        var f = text_pos if not reverse else 2 * self.half - 1 - text_pos - length
        var lo = self._run_of(f)
        var seq = self.run_seq[lo]
        return SmemHit(UInt32(self.offsets[seq] + self.run_ref[lo] + f - self.run_text[lo]), reverse)

    def sequence_of(self, offset: Int) -> Int:
        """Index of the sequence containing reference text offset `offset`."""
        var lo = 0
        var hi = len(self.offsets) - 1
        while lo < hi:
            var mid = (lo + hi + 1) >> 1
            if self.offsets[mid] <= offset:
                lo = mid
            else:
                hi = mid - 1
        return lo

    # --------------------------------------------------------------- SMEMs

    def smems[
        mut: Bool, //, origin: Origin[mut=mut]
    ](self, queries: List[Span[UInt8, origin]], min_len: Int, max_hits: Int = 20) -> ReadSmems:
        """The SMEMs of length >= `min_len` of each query against the
        reference, both strands: query intervals of maximal exact matches
        not contained in another's. Hits are listed for SMEMs with at most
        `max_hits` of them."""
        var bounds = List[Int](capacity=len(queries) + 1)
        bounds.append(0)
        for q in queries:
            bounds.append(bounds[len(bounds) - 1] + len(q))
        var codes = List[UInt8](unsafe_uninit_length=bounds[len(bounds) - 1] + 32)
        for i in range(len(queries)):
            encode_bases(queries[i], codes.unsafe_ptr().unsafe_offset(bounds[i]))
        return self.find_smems(codes, bounds, min_len, max_hits)

    def find_smems(self, codes: List[UInt8], bounds: List[Int], min_len: Int, max_hits: Int) -> ReadSmems:
        """`smems` on queries already encoded with `encode_bases`: query `q`
        is `codes[bounds[q] : bounds[q + 1]]`, and `codes` holds 32 more
        bytes (of any value) after the last query."""
        var found = List[_Found]()
        self._search(codes.unsafe_ptr(), bounds, max(min_len, 1), max_hits >= 1, found)
        var nq = len(bounds) - 1
        var result = ReadSmems()
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
        self.locate(rows)
        result.smems = List[Smem](capacity=len(order))
        result.hits = List[SmemHit](capacity=len(rows))
        var next = 0
        for f in order:
            var first = len(result.hits)
            if f.pos >= 0:
                result.hits.append(self.reference_hit(f.pos, f.end - f.start))
            elif f.size <= max_hits:
                for _ in range(f.size):
                    result.hits.append(self.reference_hit(rows[next], f.end - f.start))
                    next += 1
                _sort_hits(result.hits, first)
            result.smems.append(Smem(UInt32(f.start), UInt32(f.end), first, f.size))
        return result^

    def _search(
        self, codes: ImmPointer[UInt8, _], bounds: List[Int], min_len: Int, locate: Bool, mut found: List[_Found]
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
            if self._step(active[j], codes, min_len, locate, found):
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
        self, mut s: _Search, codes: ImmPointer[UInt8, _], min_len: Int, locate: Bool, mut found: List[_Found]
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
                self._prefetch_row(s.row)
                return False
            if s.i == s.length:
                found.append(_Found(s.query, s.x, s.i, s.p.lo, s.p.size, -1))
                return True
            var c = Int(q[unsafe_offset=s.i])
            if c < 4:
                var ok = self.forward(s.p, c)
                if ok.size > 0:
                    s.p = ok
                    s.i += 1
                    self._prefetch_interval(ok.rc, ok.size)
                    return False
            found.append(_Found(s.query, s.x, s.i, s.p.lo, s.p.size, -1))
            self._next_start(s, q)
            return False
        if s.stage == _LOCATE:
            var row = s.row
            var pos = -1
            if self._sampled(row):
                pos = self._sample(row)
            else:
                var next = self._lf(row)
                if next < 0:
                    pos = self.separator_pos[self._separator_index(row)]
                else:
                    s.row = next
                    s.steps += 1
                    self._prefetch_row(next)
                    return False
            s.row = pos + s.steps
            s.stage = s.key
            # Prefetch the first reference base the comparison reads.
            var f: Int
            if s.stage == _DIRECT:
                var m = s.i - s.x
                f = s.row + m if s.row < self.half else 2 * self.half - 2 - s.row - m
            else:
                f = s.row - 1 if s.row < self.half else 2 * self.half - 1 - s.row
            prefetch[_AHEAD](self.packed.unsafe_ptr().unsafe_offset(max(f, 0) >> 5))
            return False
        if s.stage == _DIRECT_BACK:
            # q[i + 1, x + min_len) occurs once, at text offset `row`: compare
            # backwards down to x; then the forward extension by comparison.
            var t = s.row
            var m = s.x + min_len - 1 - s.i
            var want = s.i - s.x + 1
            var n: Int
            if t < self.half:
                var first = self.run_text[self._run_of(t)]
                n = self._match_left(q, s.i, want, s.begin, t - 1, first, False)
                t -= n
            else:
                var f = 2 * self.half - 1 - t - m
                var end = self.run_text[self._run_of(f) + 1] - 1
                n = self._match_left(q, s.i, want, s.begin, f + m, end, True)
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
            if t < self.half:
                var end = self.run_text[self._run_of(t) + 1] - 1
                s.i += self._match_right(q, s.i, s.length - s.i, t + m, end, False)
            else:
                # On the reverse strand: the match is the reverse complement
                # of forward bases [f, f + m), which extends leftwards.
                var f = 2 * self.half - 1 - t - m
                var first = self.run_text[self._run_of(f)]
                s.i += self._match_right(q, s.i, s.length - s.i, f - 1, first, True)
            found.append(_Found(s.query, s.x, s.i, 0, 1, t))
            if s.i == s.length:
                return True
            self._next_start(s, q)
            return False
        if s.stage == _BACK_FIRST or s.stage == _BACK_SECOND:
            if s.p.size == 0:
                # The base at i + 1 does not occur at all.
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
                self._prefetch_row(s.row)
                return False
            var floor = s.x if s.stage == _BACK_FIRST else s.x + 1
            if s.i < floor:
                if s.stage == _BACK_FIRST:
                    self._prefetch_interval(s.p.rc, s.p.size)
                    s.i = s.x + min_len
                    s.stage = _FORWARD
                else:
                    s.x = s.i + 1
                    s.stage = _SETUP
                return False
            var ok = self.backward(s.p, Int(q[unsafe_offset=s.i]))
            if ok.size == 0:
                s.x = s.i + 1
                s.stage = _SETUP
            else:
                s.p = ok
                s.i -= 1
                self._prefetch_interval(ok.lo, ok.size)
            return False
        if s.stage == _CACHE_FIRST or s.stage == _CACHE_SECOND:
            s.p = self._cached(s.key)
            if s.p.size == 0:
                # The k-mer does not occur: back off to base-by-base search,
                # which finds where it fails.
                s.i += self.cache_k
                s.p = self.base_interval(Int(q[unsafe_offset=s.i]))
                s.i -= 1
            self._prefetch_interval(s.p.lo, s.p.size)
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
        if self.cache_k > 0 and min_len >= self.cache_k:
            s.key = self._key(q, s.i)
            s.i -= self.cache_k
            self._prefetch_cache(s.key)
            s.stage = _CACHE_FIRST
        else:
            s.p = self.base_interval(Int(q[unsafe_offset=s.i]))
            s.i -= 1
            self._prefetch_interval(s.p.lo, s.p.size)
            s.stage = _BACK_FIRST
        return False

    @inline(.always)
    def _next_start(self, mut s: _Search, q: ImmPointer[UInt8, _]):
        """After an SMEM ending at `s.i` (inside the query): the backward pass
        from `s.i` that finds the next start."""
        var c = Int(q[unsafe_offset=s.i])
        if c >= 4:
            s.x = s.i + 1
            s.stage = _SETUP
        elif self.cache_k > 0 and s.i - s.x - 1 >= self.cache_k:
            s.key = self._key(q, s.i)
            s.i -= self.cache_k
            self._prefetch_cache(s.key)
            s.stage = _CACHE_SECOND
        else:
            s.p = self.base_interval(c)
            s.i -= 1
            self._prefetch_interval(s.p.lo, s.p.size)
            s.stage = _BACK_SECOND

    @inline(.always)
    def _prefetch_row(self, row: Int):
        """Prefetches what locating from `row` reads next."""
        if self._sampled(row):
            prefetch[_AHEAD](self.samples.unsafe_ptr().unsafe_offset((row >> self.sample_bits) * _SAMPLE_BYTES))
        else:
            self._prefetch(row)

    @inline(.always)
    def _sampled(self, row: Int) -> Bool:
        return row & ((1 << self.sample_bits) - 1) == 0

    @inline(.always)
    def _sample(self, row: Int) -> Int:
        var p = self.samples.unsafe_ptr().unsafe_offset((row >> self.sample_bits) * _SAMPLE_BYTES)
        return Int(bitcast[DType.uint64, 1](p.unsafe_load[width=8]())[0] & ((UInt64(1) << 40) - 1))

    def _sample_count(self) -> Int:
        return (self.rows + (1 << self.sample_bits) - 1) >> self.sample_bits

    # ---------------------------------------------------------- save, load

    def save(self, path: String) raises:
        """Writes the index to `path` (read it back with `FMIndex.load`)."""
        with open(path, "w") as f:
            f.write_bytes(_MAGIC.as_bytes())
            var header: List[UInt64] = [
                UInt64(_VERSION),
                UInt64(self.rows),
                UInt64(self.half),
                UInt64(len(self.separators)),
                UInt64(len(self.run_seq)),
                UInt64(len(self.names)),
                UInt64(self.cache_k),
                UInt64(self.sample_bits),
            ]
            _write_words(f, Span(header))
            _write_words(f, Span(self._blocks)[self._first : self._first + ((self.rows >> 7) + 1) * 8])
            _write_bytes(f, self.samples.unsafe_ptr(), len(self.samples))
            _write_words(f, Span(self.packed))
            _write_bytes(f, self.cache.unsafe_ptr().unsafe_bitcast[UInt8](), len(self.cache) * 16)
            _write_ints(f, self.separators)
            _write_ints(f, self.separator_pos)
            _write_ints(f, self.run_text)
            _write_ints(f, self.run_seq)
            _write_ints(f, self.run_ref)
            _write_ints(f, self.lengths)
            for name in self.names:
                var length: List[UInt64] = [UInt64(name.byte_length())]
                _write_words(f, Span(length))
                f.write_bytes(name.as_bytes())

    @staticmethod
    def load(path: String, threads: Int = 1) raises -> FMIndex:
        """Reads an index written by `save`, with up to `threads` threads."""
        var index = FMIndex()
        var ranges = List[Tuple[Int, Int, Int]]()
        with open(path, "r") as f:
            var magic = List[UInt8](length=8, fill=0)
            _read_exact(f, magic.unsafe_ptr(), 8)
            var expected = _MAGIC.as_bytes()
            for i in range(8):
                if magic[i] != expected[i]:
                    raise Error(t"{path} is not a sporos FM-index")
            var header = List[UInt64](length=8, fill=0)
            _read_words(f, header)
            if header[0] != UInt64(_VERSION):
                raise Error(t"{path}: unsupported FM-index version {header[0]}")
            index.rows = Int(header[1])
            index.half = Int(header[2])
            var nsep = Int(header[3])
            var nruns = Int(header[4])
            var nseqs = Int(header[5])
            # The large arrays are read in parallel below; skip them here.
            var nblocks = (index.rows >> 7) + 1
            index._blocks = _huge_uninit[DType.uint64](nblocks * 8 + 8)
            index._first = _aligned_first(index._blocks)
            index.sample_bits = Int(header[7])
            index.samples = _huge_uninit[DType.uint8](index._sample_count() * _SAMPLE_BYTES + 8)
            index.packed = _huge_uninit[DType.uint64]((index.half >> 5) + 2)
            index.cache_k = Int(header[6])
            if index.cache_k > 0:
                index.cache = _huge_cache(index.cache_k)
            var at = 8 + 8 * len(header)
            ranges.append((at, Int(index._blocks.unsafe_ptr().unsafe_offset(index._first)), nblocks * 64))
            at += nblocks * 64
            ranges.append((at, Int(index.samples.unsafe_ptr()), len(index.samples)))
            at += len(index.samples)
            ranges.append((at, Int(index.packed.unsafe_ptr()), len(index.packed) * 8))
            at += len(index.packed) * 8
            ranges.append((at, Int(index.cache.unsafe_ptr()), len(index.cache) * 16))
            at += len(index.cache) * 16
            _ = f.seek(at)
            index.separators = _read_ints(f, nsep)
            index.separator_pos = _read_ints(f, nsep)
            index.run_text = _read_ints(f, nruns + 1)
            index.run_seq = _read_ints(f, nruns)
            index.run_ref = _read_ints(f, nruns)
            index.lengths = _read_ints(f, nseqs)
            var offset = 0
            for s in range(nseqs):
                index.offsets.append(offset)
                offset += index.lengths[s] + 1
                var length = List[UInt64](length=1, fill=0)
                _read_words(f, length)
                var bytes = List[UInt8](length=Int(length[0]), fill=0)
                _read_exact(f, bytes.unsafe_ptr(), len(bytes))
                index.names.append(String(unsafe_from_utf8=Span(bytes)))
        var failed = List[Int](length=1, fill=0)
        for r in ranges:
            var slices = max(1, min(threads, r[2] >> 24))
            var task = _RangeReadTask(path, r[0], r[1], r[2], slices, failed.unsafe_ptr())
            parallel_for(task, slices, threads)
        if failed[0] != 0:
            raise Error(t"{path}: truncated FM-index file")
        _ = len(index.cache)  # keep every array alive while tasks write through addresses
        index.counts = index.rank4(index.rows)
        var first = Int64(len(index.separators))
        index.starts = SIMD[DType.int64, 4](
            first,
            first + index.counts[0],
            first + index.counts[0] + index.counts[1],
            first + index.counts[0] + index.counts[1] + index.counts[2],
        )
        return index^


def encode_bases(seq: Span[UInt8, _], dest: MutPointer[UInt8, _]):
    """Writes `seq` as query codes: 0-3 for ACGT (either case, U as T), 4
    for any other byte."""
    comptime W = 32
    var n = len(seq)
    var src = seq.unsafe_ptr()
    var i = 0
    while i + W <= n:
        var raw = src.unsafe_load[width=W](i)
        var lower = raw | 0x20
        var valid = lower.eq(0x61) | lower.eq(0x63) | lower.eq(0x67) | lower.eq(0x74) | lower.eq(0x75)
        var code = ((raw >> 1) ^ (raw >> 2)) & 3
        dest.unsafe_offset(i).unsafe_store(valid.select(code, SIMD[DType.uint8, W](_SEPARATOR)))
        i += W
    while i < n:
        dest[unsafe_offset=i] = _code(src[unsafe_offset=i])
        i += 1


@inline(.always)
def _code(b: UInt8) -> UInt8:
    var lower = b | 0x20
    if lower == 0x61 or lower == 0x63 or lower == 0x67 or lower == 0x74 or lower == 0x75:
        return ((b >> 1) ^ (b >> 2)) & 3
    return _SEPARATOR


def _sort_hits(mut hits: List[SmemHit], first: Int):
    for i in range(first + 1, len(hits)):
        var h = hits[i]
        var j = i
        while j > first and h < hits[j - 1]:
            hits[j] = hits[j - 1]
            j -= 1
        hits[j] = h


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
                if row & ((1 << self.sample_bits) - 1) == 0:
                    # Byte stores: a wider one would reach the next task's samples.
                    var at = (row >> self.sample_bits) * _SAMPLE_BYTES
                    for b in range(_SAMPLE_BYTES):
                        self.samples[unsafe_offset=at + b] = UInt8((s >> (8 * b)) & 0xFF)
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


comptime _IO_CHUNK = 1 << 30


def _write_words(mut f: FileHandle, words: Span[UInt64, _]) raises:
    var bytes = Span[UInt8, ImmUntrackedOrigin](
        unsafe_ptr=words.unsafe_ptr().unsafe_bitcast[UInt8]().as_imm().unsafe_origin_cast[ImmUntrackedOrigin](),
        length=len(words) * 8,
    )
    for at in range(0, len(bytes), _IO_CHUNK):
        f.write_bytes(bytes[at : min(len(bytes), at + _IO_CHUNK)])


def _write_bytes(mut f: FileHandle, data: ImmPointer[UInt8, _], length: Int) raises:
    for at in range(0, length, _IO_CHUNK):
        f.write_bytes(
            Span[UInt8, ImmUntrackedOrigin](
                unsafe_ptr=data.unsafe_offset(at).unsafe_origin_cast[ImmUntrackedOrigin](),
                length=min(_IO_CHUNK, length - at),
            )
        )


def _write_ints(mut f: FileHandle, values: List[Int]) raises:
    var words = List[UInt64](capacity=len(values))
    for v in values:
        words.append(UInt64(v))
    _write_words(f, words)


def _read_exact(mut f: FileHandle, dest: MutPointer[UInt8, _], length: Int) raises:
    var done = 0
    while done < length:
        var room = min(length - done, _IO_CHUNK)
        var got = Int(
            f.read(
                Span[UInt8, MutUntrackedOrigin](
                    unsafe_ptr=dest.unsafe_offset(done).unsafe_origin_cast[MutUntrackedOrigin](), length=room
                )
            )
        )
        if got <= 0:
            raise Error("truncated FM-index file")
        done += got


def _read_words(mut f: FileHandle, mut words: List[UInt64]) raises:
    _read_exact(f, words.unsafe_ptr().unsafe_bitcast[UInt8](), len(words) * 8)


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


struct _RangeReadTask(ParallelTask):
    """Reads slice `index` of `length` bytes at file offset `offset` to
    address `dest`."""

    var path: String
    var offset: Int
    var dest: Int
    var length: Int
    var slices: Int
    var failed: MutPointer[Int, MutUntrackedOrigin]

    def __init__(out self, path: String, offset: Int, dest: Int, length: Int, slices: Int, failed: MutPointer[Int, _]):
        self.path = path.copy()
        self.offset = offset
        self.dest = dest
        self.length = length
        self.slices = slices
        self.failed = failed.unsafe_origin_cast[MutUntrackedOrigin]()

    def run(self, index: Int):
        var per = (self.length + self.slices - 1) // self.slices
        var lo = index * per
        var hi = min(self.length, lo + per)
        if lo >= hi:
            return
        try:
            with open(self.path, "r") as f:
                _ = f.seek(self.offset + lo)
                var dest = MutPointer[UInt8, MutUntrackedOrigin](unsafe_from_address=self.dest)
                _read_exact(f, dest.unsafe_offset(lo), hi - lo)
        except:
            self.failed[] = 1


def _huge_cache(k: Int) -> List[SIMD[DType.uint32, 4]]:
    var cache = List[SIMD[DType.uint32, 4]](length=1 << (2 * k), fill=SIMD[DType.uint32, 4](0))
    advise_huge_pages(Int(cache.unsafe_ptr()), len(cache) * 16)
    return cache^


def _aligned_first(blocks: List[UInt64]) -> Int:
    """Word offset of the first 64-byte boundary in `blocks`."""
    return ((-Int(blocks.unsafe_ptr())) & 63) >> 3


def _read_ints(mut f: FileHandle, count: Int) raises -> List[Int]:
    var words = List[UInt64](length=count, fill=0)
    _read_words(f, words)
    var values = List[Int](capacity=count)
    for w in words:
        values.append(Int(w))
    return values^
