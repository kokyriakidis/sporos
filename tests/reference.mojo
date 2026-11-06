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


def _kmer_keys(codes: List[UInt8], k: Int, w: Int, canonical: Bool) -> List[UInt64]:
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


def naive_minimizers(seq: List[UInt8], k: Int, w: Int, canonical: Bool) -> List[UInt32]:
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
    seq: List[UInt8], k: Int, s: Int, canonical: Bool, open_offset: Int = -1, *, open_offsets: List[Int] = []
) -> List[UInt32]:
    var opened = open_offsets.copy()
    if open_offset >= 0:
        opened.append(open_offset)
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
        var hit = (p == 0 or p == w - 1) if len(opened) == 0 else p in opened
        if hit:
            out.append(UInt32(j))
    return out^


def naive_mems(reference: List[UInt8], query: List[UInt8], min_len: Int) -> List[Tuple[Int, Int, Int]]:
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


def canonical_kmer_of(codes: List[UInt8], start: Int, k: Int) -> Tuple[UInt128, Bool]:
    """(canonical 2-bit k-mer, reverse) with base i in bits 2i, as SeedIndex (k <= 64)."""
    var fwd = UInt128(0)
    var rc = UInt128(0)
    for i in range(k):
        var c = UInt128(codes[start + i])
        fwd |= c << UInt128(2 * i)
        rc |= (3 - c) << UInt128(2 * (k - 1 - i))
    if rc < fwd:
        return (rc, True)
    return (fwd, False)


def yak_test_sequence(n: Int, seed: Int) -> List[UInt8]:
    """The sequences the hifiasm sketch digests in `test_yak.mojo` were computed on: an LCG's ACGT letters with
    homopolymer runs and the odd `N`, the same generator as the C program that ran hifiasm's own sketch."""
    comptime letters = "ACGT"
    var out = List[UInt8](capacity=n)
    var state = UInt64(seed)
    while len(out) < n:
        state = state * 6364136223846793005 + 1442695040888963407
        var c = letters.as_bytes()[Int((state >> 33) & 3)]
        if (state >> 40) & 511 == 0:
            c = UInt8(ord("N"))
        var run = 2 + Int((state >> 24) % 5) if (state >> 20) & 15 == 0 else 1
        for _ in range(run):
            if len(out) < n:
                out.append(c)
    return out^


def _yak_code(c: UInt8) -> Int:
    var u = c | 0x20
    if u == UInt8(ord("a")):
        return 0
    if u == UInt8(ord("c")):
        return 1
    if u == UInt8(ord("g")):
        return 2
    if u == UInt8(ord("t")) or u == UInt8(ord("u")):
        return 3
    return 4


def naive_hpc_kmer_planes(seq: List[UInt8], k: Int) -> List[Tuple[UInt64, UInt64, UInt64, UInt64]]:
    """Every homopolymer-compressed k-mer of `seq` (runs of non-ACGT letters cut it) as yak's four planes, built
    from scratch for each window: `x0`, `x1` the forward low and high bits with the oldest letter in bit k-1,
    `x2`, `x3` the reverse complement's with the oldest letter in bit 0."""
    var out = List[Tuple[UInt64, UInt64, UInt64, UInt64]]()
    var segment = List[Int]()

    def flush(segment: List[Int]) {mut out, k}:
        for start in range(len(segment) - k + 1):
            var x0: UInt64 = 0
            var x1: UInt64 = 0
            var x2: UInt64 = 0
            var x3: UInt64 = 0
            for j in range(k):
                var c = segment[start + j]
                x0 |= UInt64(c & 1) << UInt64(k - 1 - j)
                x1 |= UInt64(c >> 1) << UInt64(k - 1 - j)
                x2 |= UInt64(1 - (c & 1)) << UInt64(j)
                x3 |= UInt64(1 - (c >> 1)) << UInt64(j)
            out.append((x0, x1, x2, x3))

    for byte in seq:
        var c = _yak_code(byte)
        if c == 4:
            flush(segment)
            segment.clear()
        elif len(segment) == 0 or segment[len(segment) - 1] != c:
            segment.append(c)
    flush(segment)
    return out^
