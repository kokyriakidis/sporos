"""Padded, concatenated sequence text used by the exact-match engines.

Sequences are stored as 2-bit codes (0..3) in one buffer, separated by a
sentinel byte and surrounded by `TEXT_PADDING` sentinel bytes, so that 8-byte
comparisons can run off either end without bounds checks. The reference and
the query use different codes for N and for sentinels, so a match can never
extend through them.
"""

from std.memory import Allocation, Layout, alloc, dealloc

from sporos.alphabet import N_CODE, encode_into

comptime TEXT_PADDING = 32

comptime REF_N: UInt8 = N_CODE
comptime REF_SEP: UInt8 = 6
comptime QUERY_N: UInt8 = 5
comptime QUERY_SEP: UInt8 = 7


struct SeqText(Movable, Sized):
    """One or more sequences, concatenated with sentinels between them.

    `start(i)` is the text offset of sequence `i`; text offsets of bases are
    what the match engines report.
    """

    var _alloc: Allocation[UInt8]
    var _len: Int
    var starts: List[Int]
    var lengths: List[Int]
    var is_query: Bool

    def __init__[
        mut: Bool, //, origin: Origin[mut=mut]
    ](out self, seqs: List[Span[UInt8, origin]], is_query: Bool):
        var total = 0
        for s in seqs:
            total += len(s) + 1
        self._len = total
        self._alloc = alloc(Layout[UInt8](length=total + 2 * TEXT_PADDING))
        self.starts = List[Int](capacity=len(seqs))
        self.lengths = List[Int](capacity=len(seqs))
        self.is_query = is_query
        var n_code = QUERY_N if is_query else REF_N
        var sep = QUERY_SEP if is_query else REF_SEP
        var base = self._alloc.unsafe_ptr()
        for i in range(TEXT_PADDING):
            base.unsafe_offset(i).unsafe_write(sep)
            base.unsafe_offset(TEXT_PADDING + total + i).unsafe_write(sep)
        var p = base.unsafe_offset(TEXT_PADDING)
        var off = 0
        for s in seqs:
            self.starts.append(off)
            self.lengths.append(len(s))
            encode_into(s, p.unsafe_offset(off))
            if is_query:
                for j in range(off, off + len(s)):
                    if p[unsafe_offset=j] == N_CODE:
                        p.unsafe_offset(j).unsafe_write(n_code)
            off += len(s)
            p.unsafe_offset(off).unsafe_write(sep)
            off += 1

    def __init__(out self, *, reverse_complement_of: SeqText):
        """The reverse complement of a single-sequence query text."""
        ref src = reverse_complement_of
        self._len = src._len
        self._alloc = alloc(Layout[UInt8](length=src._len + 2 * TEXT_PADDING))
        self.starts = src.starts.copy()
        self.lengths = src.lengths.copy()
        self.is_query = src.is_query
        var sep = QUERY_SEP if src.is_query else REF_SEP
        var base = self._alloc.unsafe_ptr()
        for i in range(TEXT_PADDING):
            base.unsafe_offset(i).unsafe_write(sep)
            base.unsafe_offset(TEXT_PADDING + self._len + i).unsafe_write(sep)
        var dst = base.unsafe_offset(TEXT_PADDING)
        var s = src.ptr()
        for i in range(len(src.starts)):
            var start = src.starts[i]
            var n = src.lengths[i]
            for j in range(n):
                var c = s[unsafe_offset=start + n - 1 - j]
                dst.unsafe_offset(start + j).unsafe_write(3 - c if c < 4 else c)
            dst.unsafe_offset(start + n).unsafe_write(sep)

    def __deinit__(deinit self):
        dealloc(self._alloc^)

    @inline(.always)
    def ptr(ref self) -> Pointer[UInt8, origin_of(self._alloc._alloc)]:
        """Pointer to text offset 0 (padding precedes it)."""
        return self._alloc.unsafe_ptr().unsafe_offset(TEXT_PADDING)

    def __len__(self) -> Int:
        return self._len

    def sequence_of(self, offset: Int) -> Int:
        """Index of the sequence containing text offset `offset`."""
        var lo = 0
        var hi = len(self.starts) - 1
        while lo < hi:
            var mid = (lo + hi + 1) // 2
            if self.starts[mid] <= offset:
                lo = mid
            else:
                hi = mid - 1
        return lo
