"""hifiasm's k-mer representation and hash (htab.h, htab.cpp).

A k-mer of up to 63 letters is four 64-bit bit planes: `x0`, `x1` the forward k-mer's low and high bit of
every letter (A=00, C=01, G=10, T=11, newest letter in bit 0), `x2`, `x3` the reverse complement's (newest
letter in bit k-1). The canonical hash is `yak_hash64_64` of the two planes of the smaller strand, summed.
"""


from std.bit import count_trailing_zeros
from std.memory import pack_bits
from std.sys import get_defined_int, simd_width_of

from sporos.alphabet import encode_into


comptime _LETTERS: UInt32 = (1 << 0) | (1 << 2) | (1 << 6) | (1 << 19) | (1 << 20)
"""The lowercase letters a, c, g, t and u as bits counted from `a`."""


@inline(.always)
def nt4(c: UInt8) -> UInt8:
    """`seq_nt4_table`: A/a 0, C/c 1, G/g 2, T/t/U/u 3, anything else 4. Without branches: `(c >> 1) ^ (c >> 2)`
    maps the letters' ASCII codes to their 2-bit codes, and a bit mask tells the letters from anything else."""
    var offset = UInt32((c | 0x20) - UInt8(ord("a")))
    var valid = offset < 32 and ((_LETTERS >> (offset & 31)) & 1) == 1
    var code = ((c >> 1) ^ (c >> 2)) & 3
    return code if valid else 4


@inline(.always)
def yak_hash64_64(key: UInt64) -> UInt64:
    """htab.h:150, Thomas Wang's 64-bit mix, wrapping."""
    var k = key
    k = ~k + (k << 21)
    k = k ^ (k >> 24)
    k = (k + (k << 3)) + (k << 8)
    k = k ^ (k >> 14)
    k = (k + (k << 2)) + (k << 4)
    k = k ^ (k >> 28)
    k = k + (k << 31)
    return k


comptime HASH_LANES = get_defined_int["SPOROS_HASH_LANES", simd_width_of[DType.uint64]()]()
"""Keys hashed at once: the target's native 64-bit lanes (8 with AVX-512, 4 with AVX2, 2 with NEON) unless the
build sets `-D SPOROS_HASH_LANES=n`. The hash is wrapping integer arithmetic, the same in every lane width."""


@inline(.always)
def yak_hash64_64_lanes(key: SIMD[DType.uint64, HASH_LANES]) -> SIMD[DType.uint64, HASH_LANES]:
    """`yak_hash64_64` of `HASH_LANES` keys at once, lane by lane the same wrapping arithmetic."""
    var k = key
    k = ~k + (k << 21)
    k = k ^ (k >> 24)
    k = (k + (k << 3)) + (k << 8)
    k = k ^ (k >> 14)
    k = (k + (k << 2)) + (k << 4)
    k = k ^ (k >> 28)
    k = k + (k << 31)
    return k


@fieldwise_init
struct Kmer(ImplicitlyCopyable, TrivialRegisterPassable):
    """The four bit planes of a rolling k-mer (`uint64_t x[4]`), plus the mask and shift of its length."""

    var x0: UInt64
    var x1: UInt64
    var x2: UInt64
    var x3: UInt64
    var mask: UInt64
    var shift: UInt64

    @staticmethod
    def empty(k: Int) -> Kmer:
        return Kmer(0, 0, 0, 0, (UInt64(1) << UInt64(k)) - 1, UInt64(k - 1))

    @inline(.always)
    def push(mut self, c: UInt8):
        """Appends letter code `c` (0..3): the forward planes shift up, the reverse ones shift down."""
        var low = UInt64(c & 1)
        var high = UInt64(c >> 1)
        self.x0 = ((self.x0 << 1) | low) & self.mask
        self.x1 = ((self.x1 << 1) | high) & self.mask
        self.x2 = (self.x2 >> 1) | ((1 - low) << self.shift)
        self.x3 = (self.x3 >> 1) | ((1 - high) << self.shift)

    @inline(.always)
    def clear(mut self):
        self.x0 = 0
        self.x1 = 0
        self.x2 = 0
        self.x3 = 0

    @inline(.always)
    def symmetric(self) -> Bool:
        """`kmer[1] == kmer[3]`: the sketch skips these, their strand being unknown."""
        return self.x1 == self.x3

    @inline(.always)
    def reverse_is_canonical(self) -> Bool:
        """`z = kmer[1] < kmer[3] ? 0 : 1`."""
        return not (self.x1 < self.x3)

    @inline(.always)
    def canonical(self) -> Tuple[UInt64, UInt64]:
        """The canonical strand's two planes, the ones `hash` hashes."""
        if self.x1 < self.x3:
            return (self.x0, self.x1)
        return (self.x2, self.x3)

    @inline(.always)
    def hash(self) -> UInt64:
        """`yak_hash_long`: the hash of the canonical strand's two planes."""
        if self.x1 < self.x3:
            return yak_hash64_64(self.x0) + yak_hash64_64(self.x1)
        return yak_hash64_64(self.x2) + yak_hash64_64(self.x3)


comptime _RUN_WIDTH = 32


def hpc_run_starts(seq: Span[UInt8, _], mut codes: List[UInt8], mut starts: List[UInt32], is_hpc: Bool = True) -> Int:
    """The runs a homopolymer-compressed scan of `seq` steps through: `codes[1 + i]` is base `i`'s yak code (4
    for anything but A/C/G/T/U), and `starts[r]` where run `r` begins, `starts[runs] = len(seq)`; returns the
    number of runs. A run is a stretch of one letter (with `is_hpc`; else every base is its own run), and every
    ambiguous base is a run of its own, as hifiasm steps over them one by one."""
    var n = len(seq)
    fit_codes(codes, n)
    encode_into(seq, codes.unsafe_ptr().unsafe_offset(1), 4)
    return run_starts(codes, n, starts, is_hpc)


def fit_codes(mut codes: List[UInt8], n: Int):
    """Sizes a codes buffer for `run_starts` on `n` bases: a sentinel, the codes, and a padded tail (packed reads
    decode whole words)."""
    var need = n + 1 + 2 * _RUN_WIDTH
    if len(codes) < need:
        codes.resize(unsafe_uninit_length=need)
    codes.unsafe_ptr()[unsafe_offset=0] = 255


def run_starts(codes: List[UInt8], n: Int, mut starts: List[UInt32], is_hpc: Bool = True) -> Int:
    """`hpc_run_starts` on codes already in `codes[1 : n + 1]` (`codes[0]` a sentinel no code equals, as
    `fit_codes` leaves it). A run starts wherever a code differs from the one before it (or is 4), found
    `_RUN_WIDTH` bases at a time as the bits of a vector comparison."""
    if len(starts) < n + 1:
        starts.resize(unsafe_uninit_length=n + 1)
    var cp = codes.unsafe_ptr()
    var sp = starts.unsafe_ptr()
    var runs = 0
    if not is_hpc:
        for i in range(n):
            sp[unsafe_offset=i] = UInt32(i)
        sp[unsafe_offset=n] = UInt32(n)
        return n
    var i = 0
    while i + _RUN_WIDTH <= n:
        var cur = cp.unsafe_load[width=_RUN_WIDTH](i + 1)
        var prev = cp.unsafe_load[width=_RUN_WIDTH](i)
        var bits = pack_bits[DType.uint32](cur.ne(prev) | cur.eq(4))
        while bits != 0:
            sp[unsafe_offset=runs] = UInt32(i) + UInt32(count_trailing_zeros(bits))
            runs += 1
            bits &= bits - 1
        i += _RUN_WIDTH
    while i < n:
        var c = cp[unsafe_offset=i + 1]
        if c != cp[unsafe_offset=i] or c == 4:
            sp[unsafe_offset=runs] = UInt32(i)
            runs += 1
        i += 1
    sp[unsafe_offset=runs] = UInt32(n)
    return runs
