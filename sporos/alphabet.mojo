"""Nucleotide alphabet: ASCII to 2-bit codes (A=0, C=1, G=2, T/U=3, other=4)."""

from std.memory import Allocation, Layout, alloc, dealloc

comptime N_CODE: UInt8 = 4
"""Code for any byte that is not A/C/G/T/U in either case."""

comptime _ENCODE_WIDTH = 32


@inline(.always)
def encode_base(c: UInt8) -> UInt8:
    """Returns the 2-bit code of an ASCII nucleotide, or `N_CODE`."""
    var lower = c | 0x20
    var valid = lower == 0x61 or lower == 0x63 or lower == 0x67 or lower == 0x74 or lower == 0x75
    # (c >> 1) ^ (c >> 2) maps a/A->0, c/C->1, g/G->2, t/T/u/U->3.
    return ((c >> 1) ^ (c >> 2)) & 3 if valid else N_CODE


def encode_into(seq: Span[UInt8, _], dst: MutPointer[UInt8, _], n_code: UInt8 = N_CODE):
    """Encodes `seq` into `dst`, which must hold `len(seq)` bytes; non-ACGT
    bytes become `n_code`."""
    comptime V = SIMD[DType.uint8, _ENCODE_WIDTH]
    var src = seq.unsafe_ptr()
    var n = len(seq)
    var i = 0
    while i + _ENCODE_WIDTH <= n:
        var c = src.unsafe_load[width=_ENCODE_WIDTH](i)
        var lower = c | 0x20
        var valid = lower.eq(0x61) | lower.eq(0x63) | lower.eq(0x67) | lower.eq(0x74) | lower.eq(0x75)
        var code = ((c >> 1) ^ (c >> 2)) & 3
        dst.unsafe_store(i, valid.select(code, V(n_code)))
        i += _ENCODE_WIDTH
    while i < n:
        var c = encode_base(src[unsafe_offset=i])
        dst.unsafe_offset(i).unsafe_write(n_code if c == N_CODE else c)
        i += 1


struct CodeBuffer(Movable, Sized):
    """An owned, padded array of 2-bit nucleotide codes.

    `padding` bytes of `N_CODE` follow the sequence so that SIMD kernels can
    read past the end without bounds checks.
    """

    var _alloc: Allocation[UInt8]
    var length: Int
    var padding: Int

    def __init__(out self, length: Int, padding: Int):
        self._alloc = alloc(Layout[UInt8](length=length + padding))
        self.length = length
        self.padding = padding
        var p = self._alloc.unsafe_ptr()
        for i in range(length, length + padding):
            p.unsafe_offset(i).unsafe_write(N_CODE)

    @staticmethod
    def encode(seq: Span[UInt8, _], padding: Int) -> CodeBuffer:
        var buf = CodeBuffer(len(seq), padding)
        encode_into(seq, buf._alloc.unsafe_ptr())
        return buf^

    def __deinit__(deinit self):
        dealloc(self._alloc^)

    @inline(.always)
    def unsafe_ptr(ref self) -> Pointer[UInt8, origin_of(self._alloc._alloc)]:
        return self._alloc.unsafe_ptr()

    @inline(.always)
    def __getitem__(self, i: Int) -> UInt8:
        return self.unsafe_ptr()[unsafe_offset=i]

    def __len__(self) -> Int:
        return self.length
