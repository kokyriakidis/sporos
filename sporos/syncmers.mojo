"""Closed and open syncmers.

A k-mer is a closed syncmer when the smallest of its `k - s + 1` s-mers
(by `window.order_key` of the 32-bit seed hash with `w = k - s + 1`,
canonical by default, ties broken leftmost) starts at offset 0 or `k - s`,
and an open syncmer when it starts at one of the chosen offsets `t` (0-based;
one offset, or a set of them). K-mers containing a non-ACGT base are never
syncmers. Positions are reported in increasing order.
"""

from sporos.window import (
    CLOSED_SYNCMERS,
    OPEN_SYNCMER_SET,
    OPEN_SYNCMERS,
    ScanParams,
    scan,
    scan_ascii,
)


def _check(k: Int, s: Int, length: Int) raises:
    if length >= Int(UInt32.MAX):
        raise Error("sequences must be shorter than 2^32 bases")
    if s < 1 or s > k:
        raise Error(t"need 1 <= s <= k, got k={k} s={s}")


def open_offset_mask(k: Int, s: Int, open_offset: Optional[Int], open_offsets: List[Int]) raises -> UInt64:
    """The open-syncmer offsets as a bitmask (0 for closed syncmers)."""
    if open_offset and len(open_offsets) > 0:
        raise Error("give open_offset or open_offsets, not both")
    var mask: UInt64 = 0
    var offsets = open_offsets.copy()
    if case .Some(t) = open_offset:
        offsets.append(t)
    for t in offsets:
        if t < 0 or t > k - s:
            raise Error(t"open syncmer offset must be in [0, {k - s}], got {t}")
        if t >= 64:
            raise Error(t"open syncmer offsets must be below 64, got {t}")
        mask |= UInt64(1) << UInt64(t)
    return mask


def syncmer_params(k: Int, s: Int, mask: UInt64) -> ScanParams:
    return ScanParams(k=s, w=k - s + 1, span=k, offsets=mask)


def syncmer_mode(mask: UInt64) -> Int:
    """The scan mode for an offset mask: one offset is a cheaper test than a set."""
    if mask == 0:
        return CLOSED_SYNCMERS
    return OPEN_SYNCMERS if mask & (mask - 1) == 0 else OPEN_SYNCMER_SET


def syncmer_positions_from_codes(
    codes: ImmPointer[UInt8, _],
    length: Int,
    k: Int,
    s: Int,
    canonical: Bool = True,
    open_offset: Optional[Int] = None,
    threads: Int = 1,
    *,
    open_offsets: List[Int] = List[Int](),
) raises -> List[UInt32]:
    """Syncmer positions of a pre-encoded sequence.

    `codes` must be followed by at least `window.PADDING` bytes of `N_CODE`.
    Closed syncmers unless `open_offset` or `open_offsets` is given.
    """
    _check(k, s, length)
    var kmers = length - k + 1
    if kmers <= 0:
        return List[UInt32]()
    var mask = open_offset_mask(k, s, open_offset, open_offsets)
    var params = syncmer_params(k, s, mask)
    if syncmer_mode(mask) == OPEN_SYNCMER_SET:
        if canonical:
            return scan[OPEN_SYNCMER_SET, True](codes, kmers, params, threads)
        return scan[OPEN_SYNCMER_SET, False](codes, kmers, params, threads)
    if mask != 0:
        if canonical:
            return scan[OPEN_SYNCMERS, True](codes, kmers, params, threads)
        return scan[OPEN_SYNCMERS, False](codes, kmers, params, threads)
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
    *,
    open_offsets: List[Int] = List[Int](),
) raises -> List[UInt32]:
    """Start positions of the syncmers of an ASCII DNA sequence: closed, or
    open at `open_offset` or any of `open_offsets`."""
    _check(k, s, len(seq))
    var kmers = len(seq) - k + 1
    if kmers <= 0:
        return List[UInt32]()
    var mask = open_offset_mask(k, s, open_offset, open_offsets)
    var params = syncmer_params(k, s, mask)
    if syncmer_mode(mask) == OPEN_SYNCMER_SET:
        if canonical:
            return scan_ascii[OPEN_SYNCMER_SET, True](seq, kmers, params, threads)
        return scan_ascii[OPEN_SYNCMER_SET, False](seq, kmers, params, threads)
    if mask != 0:
        if canonical:
            return scan_ascii[OPEN_SYNCMERS, True](seq, kmers, params, threads)
        return scan_ascii[OPEN_SYNCMERS, False](seq, kmers, params, threads)
    if canonical:
        return scan_ascii[CLOSED_SYNCMERS, True](seq, kmers, params, threads)
    return scan_ascii[CLOSED_SYNCMERS, False](seq, kmers, params, threads)
