"""32-bit ntHash: a rolling, strand-aware hash of nucleotide k-mers.

The forward hash of a k-mer `b_0..b_{k-1}` is `XOR_j rotl(S[b_j], k-1-j)` and
the reverse-complement hash is `XOR_j rotl(S[3-b_j], j)`. The canonical hash is
their wrapping sum, which is identical for a k-mer and its reverse complement.
Bases coded `N_CODE` contribute nothing; kernels mask such k-mers out.
"""

from sporos.alphabet import N_CODE

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
        h ^= rotl32(seed(codes[unsafe_offset=j]), k - 1 - j)
    return h


def reverse_hash(codes: ImmPointer[UInt8, _], k: Int) -> UInt32:
    """Non-rolling reverse-complement ntHash of the k-mer at `codes`."""
    var h: UInt32 = 0
    for j in range(k):
        h ^= rotl32(seed(complement(codes[unsafe_offset=j])), j)
    return h


def canonical_hash(codes: ImmPointer[UInt8, _], k: Int) -> UInt32:
    return forward_hash(codes, k) + reverse_hash(codes, k)
