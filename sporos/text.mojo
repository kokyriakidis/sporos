"""Packed, concatenated sequence text used by the exact-match engines.

Sequences are stored 2 bits per base, 32 bases per `UInt64` word (base `i`
in bits `2 * (i % 32)`), with a parallel mask holding `0b11` for every base
that must never match: non-ACGT bases, the separator after each sequence,
and the padding words around the text. Any 32 consecutive bases, and their
mask, are one shifted 64-bit window, so k-mers are extracted in O(1) and
matches extend 32 bases per comparison.
"""

from std.memory import bitcast

comptime _PAD_WORDS = 4
"""All-masked words before and after the text (128 bases each side)."""
comptime _BASES = 32


struct SeqText(Movable, Sized):
    """One or more sequences, concatenated with a masked separator after each.

    `starts[i]` is the text offset of sequence `i`; text offsets of bases are
    what the match engines report.
    """

    var words: List[UInt64]
    var masks: List[UInt64]
    var _len: Int
    var starts: List[Int]
    var lengths: List[Int]

    def __init__[
        mut: Bool, //, origin: Origin[mut=mut]
    ](out self, seqs: List[Span[UInt8, origin]], reverse_complement: Bool = False):
        """Packs `seqs`; with `reverse_complement`, each sequence is stored as
        its reverse complement instead of its forward strand."""
        var total = 0
        for s in seqs:
            total += len(s) + 1
        self._len = total
        var nwords = (total + _BASES - 1) // _BASES + 2 * _PAD_WORDS
        self.words = List[UInt64](length=nwords, fill=0)
        self.masks = List[UInt64](length=nwords, fill=0)
        for w in range(_PAD_WORDS):
            self.masks[w] = UInt64.MAX
            self.masks[nwords - 1 - w] = UInt64.MAX
        self.starts = List[Int](capacity=len(seqs))
        self.lengths = List[Int](capacity=len(seqs))
        var off = 0
        for s in seqs:
            self.starts.append(off)
            self.lengths.append(len(s))
            _pack(s, off, reverse_complement, self.words, self.masks)
            off += len(s) + 1

    def __len__(self) -> Int:
        return self._len

    @inline(.always)
    def window(self, i: Int) -> UInt64:
        """Bases `i .. i + 31` (base `i` in the low bits), for
        `-96 <= i <= len(self) + 64`."""
        return _window(self.words.unsafe_ptr(), i)

    @inline(.always)
    def mask(self, i: Int) -> UInt64:
        """The never-match mask of the same 32 bases as `window(i)`."""
        return _window(self.masks.unsafe_ptr(), i)

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


@inline(.always)
def _window(words: ImmPointer[UInt64, _], i: Int) -> UInt64:
    var at = i + _PAD_WORDS * _BASES
    var w = at >> 5
    var shift = UInt64(2 * (at & 31))
    var lo = words[unsafe_offset=w] >> shift
    # (hi << 1) << (63 - shift) avoids an undefined shift by 64 when shift = 0.
    var hi = (words[unsafe_offset=w + 1] << 1) << (63 - shift)
    return lo | hi


@inline(.always)
def _pack32(codes: SIMD[DType.uint8, 32]) -> UInt64:
    """Packs 32 two-bit values (one per byte) into a word, value i at bit 2i."""
    var x = bitcast[DType.uint16, 16](codes)
    var y = ((x | (x >> 6)) & 0x000F).cast[DType.uint8]()
    var z = bitcast[DType.uint16, 8](y)
    var b = ((z | (z >> 4)) & 0x00FF).cast[DType.uint8]()
    return bitcast[DType.uint64, 1](b)[0]


def _pack(
    seq: Span[UInt8, _],
    offset: Int,
    reverse_complement: Bool,
    mut words: List[UInt64],
    mut masks: List[UInt64],
):
    """ORs `seq` (or its reverse complement) into the text at offset `offset`,
    32 bases at a time; `masks` must be zero there and gets 0b11 per non-ACGT
    base. The separator after the sequence is masked too."""
    comptime W = 32
    var n = len(seq)
    var src = seq.unsafe_ptr()
    var base = offset + _PAD_WORDS * _BASES
    var wp = words.unsafe_ptr()
    var mp = masks.unsafe_ptr()
    var i = 0
    while i < n:
        var take = min(W, n - i)
        var raw = SIMD[DType.uint8, W](ord("N"))
        if reverse_complement:
            # Bases n-1-i, n-2-i, ...: load the 32 bytes ending at n - i.
            if take == W:
                raw = src.unsafe_load[width=W](n - i - W).reversed()
            else:
                for j in range(take):
                    raw[j] = src[unsafe_offset=n - 1 - i - j]
        else:
            if take == W:
                raw = src.unsafe_load[width=W](i)
            else:
                for j in range(take):
                    raw[j] = src[unsafe_offset=i + j]
        var lower = raw | 0x20
        var valid = (
            lower.eq(0x61) | lower.eq(0x63) | lower.eq(0x67) | lower.eq(0x74) | lower.eq(0x75)
        )
        var code = ((raw >> 1) ^ (raw >> 2)) & 3
        if reverse_complement:
            code = 3 - code
        var keep: UInt64 = UInt64.MAX if take == W else (UInt64(1) << UInt64(2 * take)) - 1
        var value = _pack32(valid.select(code, SIMD[DType.uint8, W](0))) & keep
        var bad = _pack32(valid.select(SIMD[DType.uint8, W](0), SIMD[DType.uint8, W](3))) & keep
        var at = base + i
        var w = at >> 5
        var shift = UInt64(2 * (at & 31))
        wp[unsafe_offset=w] |= value << shift
        mp[unsafe_offset=w] |= bad << shift
        if shift != 0:
            wp[unsafe_offset=w + 1] |= value >> (64 - shift)
            mp[unsafe_offset=w + 1] |= bad >> (64 - shift)
        i += W
    var sep = base + n
    mp[unsafe_offset=sep >> 5] |= UInt64(3) << UInt64(2 * (sep & 31))
