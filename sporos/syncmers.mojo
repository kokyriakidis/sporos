"""Closed and open syncmers.

A k-mer is a closed syncmer when the smallest of its `k - s + 1` s-mers
(by `window.order_key` of the 32-bit ntHash with `w = k - s + 1`, canonical
by default, ties broken leftmost) starts at offset 0
or `k - s`, and an open syncmer when it starts at offset `t`. K-mers containing
a non-ACGT base are never syncmers. Positions are reported in increasing order.
"""

from sporos.alphabet import CodeBuffer
from sporos.window import (
    CLOSED_SYNCMERS,
    OPEN_SYNCMERS,
    PADDING,
    ScanParams,
    scan,
)


def _check(k: Int, s: Int) raises:
    if s < 1 or s > k:
        raise Error(t"need 1 <= s <= k, got k={k} s={s}")


def syncmer_positions_from_codes(
    codes: ImmPointer[UInt8, _],
    length: Int,
    k: Int,
    s: Int,
    canonical: Bool = True,
    open_offset: Optional[Int] = None,
    threads: Int = 1,
) raises -> List[UInt32]:
    """Syncmer positions of a pre-encoded sequence.

    `codes` must be followed by at least `window.PADDING` bytes of `N_CODE`.
    Closed syncmers unless `open_offset` is given.
    """
    _check(k, s)
    var kmers = length - k + 1
    if kmers <= 0:
        return List[UInt32]()
    var w = k - s + 1
    if case .Some(t) = open_offset:
        if t < 0 or t > k - s:
            raise Error(t"open syncmer offset must be in [0, {k - s}], got {t}")
        var params = ScanParams(k=s, w=w, span=k, offset=t)
        if canonical:
            return scan[OPEN_SYNCMERS, True](codes, kmers, params, threads)
        return scan[OPEN_SYNCMERS, False](codes, kmers, params, threads)
    var params = ScanParams(k=s, w=w, span=k, offset=0)
    if canonical:
        return scan[CLOSED_SYNCMERS, True](codes, kmers, params, threads)
    return scan[CLOSED_SYNCMERS, False](codes, kmers, params, threads)


def syncmer_positions(
    seq: Span[UInt8, _],
    k: Int,
    s: Int,
    canonical: Bool = True,
    open_offset: Optional[Int] = None,
    threads: Int = 1,
) raises -> List[UInt32]:
    """Start positions of the syncmers of an ASCII DNA sequence."""
    var codes = CodeBuffer.encode(seq, PADDING)
    return syncmer_positions_from_codes(
        codes.unsafe_ptr(), len(seq), k, s, canonical, open_offset, threads
    )
