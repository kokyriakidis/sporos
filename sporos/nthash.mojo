"""32-bit ntHash: a rolling, strand-aware hash of nucleotide k-mers.

The forward hash of a k-mer `b_0..b_{k-1}` is `XOR_j rotl(S[b_j], R(k-1-j))`
and the reverse-complement hash is `XOR_j rotl(S[3-b_j], R j)`, with rotation
step `R = ROTATION`. A step of 7 rather than ntHash's 1 decorrelates the high
bits of consecutive k-mers' hashes, which the minimizer order compares. The canonical hash is
their wrapping sum, which is identical for a k-mer and its reverse complement.
Bases coded `N_CODE` contribute nothing; kernels mask such k-mers out.
"""

from std.sys.defines import get_defined_int

from sporos.alphabet import N_CODE

comptime ROTATION = 7
"""Bits each base's seed rotates by per position."""

comptime SEED_A: UInt32 = 0x3C8BFBB3
comptime SEED_C: UInt32 = 0x3193C185
comptime SEED_G: UInt32 = 0x20323ED0
comptime SEED_T: UInt32 = 0x295549F5


@inline(.always)
def seed(code: UInt8) -> UInt32:
    """Returns the ntHash seed of a 2-bit code (0 for `N_CODE`)."""
    if code == 0:
        return SEED_A
    if code == 1:
        return SEED_C
    if code == 2:
        return SEED_G
    if code == 3:
        return SEED_T
    return 0


@inline(.always)
def rotl32(x: UInt32, r: Int) -> UInt32:
    var s = r & 31
    if s == 0:
        return x
    return (x << UInt32(s)) | (x >> UInt32(32 - s))


@inline(.always)
def rotr32(x: UInt32, r: Int) -> UInt32:
    return rotl32(x, 32 - (r & 31))


@inline(.always)
def complement(code: UInt8) -> UInt8:
    return 3 - code if code < N_CODE else code


def forward_hash(codes: ImmPointer[UInt8, _], k: Int) -> UInt32:
    """Non-rolling forward ntHash of the k-mer at `codes`."""
    var h: UInt32 = 0
    for j in range(k):
        h ^= rotl32(seed(codes[unsafe_offset=j]), ROTATION * (k - 1 - j))
    return h


def reverse_hash(codes: ImmPointer[UInt8, _], k: Int) -> UInt32:
    """Non-rolling reverse-complement ntHash of the k-mer at `codes`."""
    var h: UInt32 = 0
    for j in range(k):
        h ^= rotl32(seed(complement(codes[unsafe_offset=j])), ROTATION * j)
    return h


def canonical_hash(codes: ImmPointer[UInt8, _], k: Int) -> UInt32:
    return forward_hash(codes, k) + reverse_hash(codes, k)


comptime SMALL_K = get_defined_int["SPOROS_SMALL_K", 0]()
"""K-mers up to this length are hashed by `small_hash` instead of ntHash.

Off by default: on AVX-512 it was no faster than ntHash, and its minimizer
order is slightly biased (k=15, w=10 density 0.1828 vs 0.1820 for ntHash;
ideal 0.1818). `-D SPOROS_SMALL_K=16` enables it, e.g. to try NEON, where
32-bit vector multiplies are cheap."""
comptime MUL_SEED: UInt32 = 0x9E3779B1


def small_hash(codes: ImmPointer[UInt8, _], k: Int, canonical: Bool) -> UInt32:
    """Multiplicative hash of a 2-bit k-mer (k <= 16, first base in the high
    bits): the k-mer, or with `canonical` the smaller of it and its reverse
    complement, times an odd constant (wrapping)."""
    var fwd: UInt32 = 0
    var rc: UInt32 = 0
    for j in range(k):
        var c = UInt32(codes[unsafe_offset=j] & 3)
        fwd = (fwd << 2) | c
        rc = (rc >> 2) | ((c ^ 3) << UInt32(2 * (k - 1)))
    if canonical:
        return min(fwd, rc) * MUL_SEED
    return fwd * MUL_SEED


def seed_hash(codes: ImmPointer[UInt8, _], k: Int, canonical: Bool) -> UInt32:
    """The hash that orders seeds: `small_hash` for k <= `SMALL_K`, else
    ntHash (canonical: forward + reverse complement)."""
    if k <= SMALL_K:
        return small_hash(codes, k, canonical)
    return canonical_hash(codes, k) if canonical else forward_hash(codes, k)
