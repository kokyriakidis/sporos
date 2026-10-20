"""Brute-force reference implementations used as test oracles."""

from sporos.alphabet import encode_base
from sporos.nthash import seed_hash
from sporos.window import order_key


def random_dna(n: Int, rng_seed: Int, n_rate: Float64 = 0.0) -> List[UInt8]:
    """Random ACGT bytes (mixed case) with occasional N runs (xorshift64*)."""
    comptime bases = "ACGTacgt"
    var state = UInt64(rng_seed) * 0x9E3779B97F4A7C15 + 0x632BE59BD9B4E019
    var threshold = UInt64(n_rate * 18446744073709551615.0)
    var out = List[UInt8](capacity=n)
    var i = 0
    while i < n:
        state ^= state >> 12
        state ^= state << 25
        state ^= state >> 27
        var r = state * 0x2545F4914F6CDD1D
        if n_rate > 0 and r < threshold:
            var run = Int((r >> 20) % 40) + 1
            for _ in range(min(run, n - i)):
                out.append(UInt8(ord("N")))
                i += 1
            continue
        out.append(bases.as_bytes()[Int(r >> 61)])
        i += 1
    return out^


def encode_all(seq: List[UInt8]) -> List[UInt8]:
    var codes = List[UInt8](capacity=len(seq))
    for c in seq:
        codes.append(encode_base(c))
    return codes^


def _kmer_keys(
    codes: List[UInt8], k: Int, w: Int, canonical: Bool
) -> List[UInt64]:
    var keys = List[UInt64]()
    for i in range(len(codes) - k + 1):
        var ok = True
        for j in range(k):
            if codes[i + j] > 3:
                ok = False
        if not ok:
            keys.append(UInt64.MAX)
            continue
        var p = codes.unsafe_ptr().unsafe_offset(i)
        var h = seed_hash(p, k, canonical)
        keys.append((UInt64(order_key(h, w)) << 32) | UInt64(i))
    return keys^


def naive_minimizers(
    seq: List[UInt8], k: Int, w: Int, canonical: Bool
) -> List[UInt32]:
    var codes = encode_all(seq)
    var out = List[UInt32]()
    if len(codes) < k:
        return out^
    var keys = _kmer_keys(codes, k, w, canonical)
    var windows = max(1, len(keys) - w + 1)
    for j in range(windows):
        var best = UInt64.MAX
        for i in range(j, min(j + w, len(keys))):
            best = min(best, keys[i])
        if best == UInt64.MAX:
            continue
        var p = UInt32(best & 0xFFFFFFFF)
        if len(out) == 0 or out[len(out) - 1] != p:
            out.append(p)
    return out^


def naive_syncmers(
    seq: List[UInt8], k: Int, s: Int, canonical: Bool, open_offset: Int = -1
) -> List[UInt32]:
    var codes = encode_all(seq)
    var out = List[UInt32]()
    if len(codes) < k:
        return out^
    var w = k - s + 1
    var keys = _kmer_keys(codes, s, w, canonical)
    for j in range(len(codes) - k + 1):
        var ok = True
        for i in range(j, j + k):
            if codes[i] > 3:
                ok = False
        if not ok:
            continue
        var best = UInt64.MAX
        for i in range(j, j + w):
            best = min(best, keys[i])
        var p = Int(best & 0xFFFFFFFF) - j
        var hit = (p == 0 or p == w - 1) if open_offset < 0 else p == open_offset
        if hit:
            out.append(UInt32(j))
    return out^


def naive_mems(
    reference: List[UInt8], query: List[UInt8], min_len: Int
) -> List[Tuple[Int, Int, Int]]:
    """All (ref_pos, query_pos, length) maximal matches of ACGT bases, sorted
    by query then reference position. Non-ACGT bases never match."""
    var r = encode_all(reference)
    var q = encode_all(query)
    var out = List[Tuple[Int, Int, Int]]()
    for j in range(len(q)):
        for i in range(len(r)):
            if r[i] > 3 or r[i] != q[j]:
                continue
            if i > 0 and j > 0 and r[i - 1] < 4 and r[i - 1] == q[j - 1]:
                continue
            var n = 0
            while i + n < len(r) and j + n < len(q) and r[i + n] < 4 and r[i + n] == q[j + n]:
                n += 1
            if n >= min_len:
                out.append((i, j, n))
    return out^


def naive_count(text: List[UInt8], pattern_start: Int, length: Int, source: List[UInt8]) -> Int:
    """Occurrences of source[pattern_start : pattern_start + length] in text."""
    var t = encode_all(text)
    var s = encode_all(source)
    var count = 0
    for i in range(len(t) - length + 1):
        var ok = True
        for n in range(length):
            if t[i + n] != s[pattern_start + n]:
                ok = False
                break
        if ok:
            count += 1
    return count


def canonical_kmer_of(codes: List[UInt8], start: Int, k: Int) -> Tuple[UInt64, Bool]:
    """(canonical 2-bit k-mer, reverse) with base i in bits 2i, as SeedIndex."""
    var fwd = UInt64(0)
    var rc = UInt64(0)
    for i in range(k):
        var c = UInt64(codes[start + i])
        fwd |= c << UInt64(2 * i)
        rc |= (3 - c) << UInt64(2 * (k - 1 - i))
    if rc < fwd:
        return (rc, True)
    return (fwd, False)
