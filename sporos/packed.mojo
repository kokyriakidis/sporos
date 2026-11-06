"""Reads kept two bits per base, as hifiasm keeps them (`All_reads`, `ha_compress_base`): A/C/G/T (and U, any case)
as the codes 0 to 3, every other letter remembered only as a position, which decodes to the code 4. That is all
yak's and hifiasm's scans read of a sequence, at a quarter of its bytes.

A read starts on a 64-bit word, 32 bases to a word, base `L` of a word at bit `8 * (L % 8) + 2 * (L / 8)`. That
order makes both directions a few vector instructions: a word decodes by being repeated over four lanes, viewed
as 32 bytes and shifted right lane by lane (byte `L` then holds base `L`'s code in its low bits), and 32 codes
encode by the reverse shift and an OR of the four lanes.
"""

from std.bit import count_trailing_zeros
from std.memory import bitcast, pack_bits

from sporos.alphabet import encode_into
from sporos.fasta import FastaRecords

comptime BASES_PER_WORD = 32
comptime AMBIGUOUS: UInt8 = 4
"""The code of a letter other than A/C/G/T/U."""

comptime _Bytes = SIMD[DType.uint8, BASES_PER_WORD]


@inline(.always)
def _field_shifts() -> _Bytes:
    """Lane `L`'s shift: `2 * (L / 8)`."""
    var s = _Bytes(0)
    comptime for lane in range(BASES_PER_WORD):
        s[lane] = UInt8(2 * (lane // 8))
    return s


@inline(.always)
def pack_word(codes: _Bytes) -> UInt64:
    """32 codes (0 to 3) as one word."""
    return bitcast[DType.uint64, 4](codes << _field_shifts()).reduce_or()


@inline(.always)
def unpack_word(word: UInt64) -> _Bytes:
    """The 32 codes of a word."""
    return (bitcast[DType.uint8, BASES_PER_WORD](SIMD[DType.uint64, 4](word)) >> _field_shifts()) & 3


struct PackedReads(Movable, Sized):
    """Named reads, two bits per base: `words[word_starts[r]:]` hold read `r`'s codes and
    `ambiguous[ambiguous_starts[r]:ambiguous_starts[r + 1]]` the positions of its other letters, in order."""

    var names: List[String]
    var lengths: List[Int]
    var word_starts: List[Int]
    var words: List[UInt64]
    var ambiguous: List[UInt32]
    var ambiguous_starts: List[Int]
    var _codes: List[UInt8]

    def __init__(out self):
        self.names = List[String]()
        self.lengths = List[Int]()
        self.word_starts = List[Int]()
        self.words = List[UInt64]()
        self.ambiguous = List[UInt32]()
        self.ambiguous_starts = [0]
        self._codes = List[UInt8]()

    def __init__(out self, records: FastaRecords):
        """`records`, packed."""
        self = Self()
        self.reserve(len(records), records.total_length())
        for r in range(len(records)):
            self.append(records.names[r], records.seq(r))

    def __len__(self) -> Int:
        return len(self.lengths)

    def reserve(mut self, reads: Int, bases: Int):
        """Room for `reads` more reads of `bases` more letters in all."""
        self.names.reserve(len(self.names) + reads)
        self.lengths.reserve(len(self.lengths) + reads)
        self.word_starts.reserve(len(self.word_starts) + reads)
        self.ambiguous_starts.reserve(len(self.ambiguous_starts) + reads)
        self.words.reserve(len(self.words) + bases // BASES_PER_WORD + reads)

    def append(mut self, name: String, seq: Span[UInt8, _]):
        """Adds a read: its letters coded as `encode_into` codes them, 32 at a time."""
        var n = len(seq)
        var padded = (n + BASES_PER_WORD - 1) // BASES_PER_WORD * BASES_PER_WORD
        if len(self._codes) < padded:
            self._codes.resize(unsafe_uninit_length=padded)
        var cp = self._codes.unsafe_ptr()
        encode_into(seq, cp, AMBIGUOUS)
        for i in range(n, padded):
            cp[unsafe_offset=i] = 0
        self.names.append(name)
        self.lengths.append(n)
        self.word_starts.append(len(self.words))
        var at = 0
        while at < padded:
            var codes = cp.unsafe_load[width=BASES_PER_WORD](at)
            var other = codes.eq(AMBIGUOUS)
            if other.reduce_or():
                var bits = pack_bits[DType.uint32](other)
                while bits != 0:
                    self.ambiguous.append(UInt32(at + Int(count_trailing_zeros(bits))))
                    bits &= bits - 1
                codes = other.select(_Bytes(0), codes)
            self.words.append(pack_word(codes))
            at += BASES_PER_WORD
        self.ambiguous_starts.append(len(self.ambiguous))

    def total_length(self) -> Int:
        var total = 0
        for length in self.lengths:
            total += length
        return total

    @inline(.always)
    def padded_length(self, r: Int) -> Int:
        """Read `r`'s length rounded up to whole words: what `decode` writes."""
        return (self.lengths[r] + BASES_PER_WORD - 1) // BASES_PER_WORD * BASES_PER_WORD

    def decode(self, r: Int, dst: MutPointer[UInt8, _]):
        """Writes read `r`'s codes to `dst`, which must hold `padded_length(r)` bytes (past the read's end, 0s)."""
        var words = self.words.unsafe_ptr().unsafe_offset(self.word_starts[r])
        var n_words = self.padded_length(r) // BASES_PER_WORD
        for i in range(n_words):
            dst.unsafe_store(i * BASES_PER_WORD, unpack_word(words[unsafe_offset=i]))
        for j in range(self.ambiguous_starts[r], self.ambiguous_starts[r + 1]):
            dst[unsafe_offset=Int(self.ambiguous[j])] = AMBIGUOUS

    def bytes_used(self) -> Int:
        """The sequence storage: the words and the ambiguous positions."""
        return 8 * len(self.words) + 4 * len(self.ambiguous)
