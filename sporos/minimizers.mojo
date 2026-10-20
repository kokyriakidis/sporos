"""(w, k)-minimizers.

The minimizer of a window of `w` consecutive k-mers is the k-mer with the
smallest `window.order_key` (bits 30 down to `ceil(log2(w))` of its 32-bit
ntHash, canonical by default), ties broken toward the leftmost k-mer. K-mers containing a non-ACGT base are never selected. Positions are
reported once per run of consecutive windows sharing a minimizer, in
increasing order. A sequence with fewer than `w` k-mers forms one window.
"""

from sporos.alphabet import CodeBuffer
from sporos.nthash import canonical_hash, forward_hash
from sporos.window import (
    MINIMIZERS,
    PADDING,
    ScanParams,
    order_key,
    scan,
)


def _check(k: Int, w: Int, length: Int) raises:
    if length >= Int(UInt32.MAX):
        raise Error("sequences must be shorter than 2^32 bases")
    if k < 1 or k > 64:
        raise Error(t"k must be in [1, 64], got {k}")
    if w < 1:
        raise Error(t"w must be >= 1, got {w}")


def minimizer_positions_from_codes(
    codes: ImmPointer[UInt8, _],
    length: Int,
    k: Int,
    w: Int,
    canonical: Bool = True,
    threads: Int = 1,
) raises -> List[UInt32]:
    """Minimizer positions of a pre-encoded sequence.

    `codes` must be followed by at least `window.PADDING` bytes of `N_CODE`.
    """
    _check(k, w, length)
    var kmers = length - k + 1
    if kmers <= 0:
        return List[UInt32]()
    if kmers < w:
        return _short_minimizer(codes, kmers, k, w, canonical)
    var params = ScanParams(k=k, w=w, span=k, offset=0)
    if canonical:
        return scan[MINIMIZERS, True](codes, kmers - w + 1, params, threads)
    return scan[MINIMIZERS, False](codes, kmers - w + 1, params, threads)


def minimizer_positions(
    seq: Span[UInt8, _],
    k: Int,
    w: Int,
    canonical: Bool = True,
    threads: Int = 1,
) raises -> List[UInt32]:
    """Start positions of the (w, k)-minimizers of an ASCII DNA sequence."""
    _check(k, w, len(seq))
    var kmers = len(seq) - k + 1
    if kmers <= 0:
        return List[UInt32]()
    if kmers < w:
        var codes = CodeBuffer.encode(seq, PADDING)
        return _short_minimizer(codes.unsafe_ptr(), kmers, k, w, canonical)
    var params = ScanParams(k=k, w=w, span=k, offset=0)
    if canonical:
        return scan[MINIMIZERS, True, ascii=True](seq.unsafe_ptr(), kmers - w + 1, params, threads, len(seq))
    return scan[MINIMIZERS, False, ascii=True](seq.unsafe_ptr(), kmers - w + 1, params, threads, len(seq))


def _short_minimizer(
    codes: ImmPointer[UInt8, _], kmers: Int, k: Int, w: Int, canonical: Bool
) -> List[UInt32]:
    var best = UInt64.MAX
    var run = 0
    for i in range(kmers + k - 1):
        run = run + 1 if codes[unsafe_offset=i] < 4 else 0
        var start = i - k + 1
        if run < k:
            continue
        var p = codes.unsafe_offset(start)
        var h = canonical_hash(p, k) if canonical else forward_hash(p, k)
        best = min(best, (UInt64(order_key(h, w)) << 32) | UInt64(start))
    var out = List[UInt32]()
    if best != UInt64.MAX:
        out.append(UInt32(best & 0xFFFFFFFF))
    return out^
