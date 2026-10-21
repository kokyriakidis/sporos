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
is located (suffix array samples, five bytes each, every 4, 8 or 16 rows
as the reference grows) and the match is extended by comparing bases
against a two-bit copy of the reference instead of base by base in the
index. Hits of repeated SMEMs are located in batches.

The suffix array is built with libsais (Apache-2.0, conda-forge), loaded at
run time.
"""

from std.bit import count_leading_zeros, count_trailing_zeros, pop_count
from std.memory import bitcast, pack_bits
from std.sys.intrinsics import PrefetchOptions, prefetch

from sporos.mems import SmemHit, Smems
from sporos.text import _pack32

from .build import build_index
from .io import load_index, save_index
from .search import search_encoded, search_smems

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
comptime _BATCH = 16
"""Queries searched together (state machines stepped round-robin)."""
comptime _LOCATE_BATCH = 64
comptime _ENCODE_AHEAD = 32
"""Suffix array rows between a text prefetch and its use when encoding."""
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
        return build_index(seqs, names, threads, sample_bits, wide_suffix_array)

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
        var r = off & 63
        var high = UInt64(0) - UInt64(off >> 6)
        var base = (word & _COUNT_MASK) + ((word >> 56) & (0x7F & high))
        var w = 4 + ((off >> 5) & 2)
        # Branch-free masks of the first r symbols across the two words.
        var low = (UInt64(1) << UInt64((2 * r) & 63)) - 1
        var full = UInt64(0) - UInt64(r >> 5)
        var a = p[unsafe_offset=w] & (low | full)
        var b = p[unsafe_offset=w + 1] & (low & full)
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
        var r = off & 63
        var high = UInt64(0) - UInt64(off >> 6)
        var count = p[unsafe_offset=c]
        var base = (count & _COUNT_MASK) + ((count >> 56) & (0x7F & high))
        var w = 4 + ((off >> 5) & 2)
        var low = (UInt64(1) << UInt64((2 * r) & 63)) - 1
        var full = UInt64(0) - UInt64(r >> 5)
        var pattern = UInt64(c) * _EVEN
        var x = p[unsafe_offset=w] ^ pattern
        var y = p[unsafe_offset=w + 1] ^ pattern
        var a = ~(x | (x >> 1)) & _EVEN & (low | full)
        var b = ~(y | (y >> 1)) & _EVEN & (low & full)
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
        var k = self._rank_le(ik.lo, c)
        var l = self._rank_le(ik.lo + ik.size, c)
        # Suffixes of the reverse complement are ordered by their next base:
        # separator first, then T, G, C, A (complements of the bases before).
        # Those before c's complement are the ones preceded by a separator or
        # by a base above c: all but those preceded by bases up to c.
        return Interval(Int(self.starts[c]) + k[0], ik.rc + ik.size - (l[1] - k[1]), l[0] - k[0])

    @inline(.always)
    def _rank_le(self, k: Int, c: Int) -> Tuple[Int, Int]:
        """Occurrences in rows `[0, k)` of base `c`, and of bases up to `c`."""
        var p = self._block_ptr().unsafe_offset((k >> 7) << 3)
        var word = p.unsafe_load[width=4]()
        var off = k & 127
        var r = off & 63
        var high = UInt64(0) - UInt64(off >> 6)
        var counts = (word & _COUNT_MASK) + ((word >> 56) & (0x7F & high))
        var lanes = SIMD[DType.int64, 4](0, 1, 2, 3)
        var up_to = Int(lanes.le(Int64(c)).select(counts, SIMD[DType.uint64, 4](0)).reduce_add())
        var w = 4 + ((off >> 5) & 2)
        var low = (UInt64(1) << UInt64((2 * r) & 63)) - 1
        var full = UInt64(0) - UInt64(r >> 5)
        # Bit 2i of `ch`/`cl` is c's high/low bit, matched per symbol.
        var ch = UInt64(0) - UInt64(c >> 1)
        var cl = UInt64(0) - UInt64(c & 1)
        var eq = 0
        var above = 0
        comptime for half in range(2):
            var x = p[unsafe_offset=w + half] & ((low | full) if half == 0 else (low & full))
            var h = (x >> 1) & _EVEN
            var l = x & _EVEN
            var mask = (low | full) if half == 0 else (low & full)
            eq += Int(pop_count(~((h ^ ch) | (l ^ cl)) & _EVEN & mask))
            above += Int(pop_count(((h & ~ch) | (~(h ^ ch) & l & ~cl)) & _EVEN))
        var rank = Int(counts[c]) + eq
        var le = up_to + r - above
        if word[0] & _SEPARATOR_FLAG != 0:
            # Separators are stored as A: counted by both, wrongly.
            var seps = self._separators_in((k >> 6) << 6, k)
            le -= seps
            if c == 0:
                rank -= seps
        return (rank, le)

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
        # Branch-free: halve the candidate range around `lo` each round.
        var lo = 0
        var n = len(self.run_seq)
        while n > 1:
            var half = n >> 1
            lo += half * Int(runs[unsafe_offset=lo + half] <= f)
            n -= half
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
        var run = self._run_of(f)
        return SmemHit(UInt32(self.run_seq[run]), UInt32(self.run_ref[run] + f - self.run_text[run]), reverse)

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
    ](
        self,
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
        often, as Mumemto's `-F` caps matches. `threads` search slices of
        queries in parallel."""
        return search_smems(self, queries, min_len, max_hits, min_occurrences, max_occurrences, threads)

    def smems_encoded(
        self,
        codes: List[UInt8],
        bounds: List[Int],
        min_len: Int,
        max_hits: Int = 20,
        min_occurrences: Int = 1,
        max_occurrences: Int = 0,
        threads: Int = 1,
    ) -> Smems:
        """`smems` on queries already encoded with `encode_bases`: query `q`
        is `codes[bounds[q] : bounds[q + 1]]` (no copy when 32 bytes follow
        the last query)."""
        return search_encoded(self, codes, bounds, min_len, max_hits, min_occurrences, max_occurrences, threads)

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
        save_index(self, path)

    @staticmethod
    def load(path: String, threads: Int = 1) raises -> FMIndex:
        """Reads an index written by `save`, with up to `threads` threads."""
        return load_index(path, threads)


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
    if i < n and n >= W:
        # The tail as the last full window, overlapping bases already done.
        i = n - W
        var raw = src.unsafe_load[width=W](i)
        var lower = raw | 0x20
        var valid = lower.eq(0x61) | lower.eq(0x63) | lower.eq(0x67) | lower.eq(0x74) | lower.eq(0x75)
        var code = ((raw >> 1) ^ (raw >> 2)) & 3
        dest.unsafe_offset(i).unsafe_store(valid.select(code, SIMD[DType.uint8, W](_SEPARATOR)))
        return
    while i < n:
        dest[unsafe_offset=i] = _code(src[unsafe_offset=i])
        i += 1


@inline(.always)
def _code(b: UInt8) -> UInt8:
    var lower = b | 0x20
    if lower == 0x61 or lower == 0x63 or lower == 0x67 or lower == 0x74 or lower == 0x75:
        return ((b >> 1) ^ (b >> 2)) & 3
    return _SEPARATOR
