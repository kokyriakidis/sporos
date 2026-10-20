from std.testing import assert_equal, assert_true, TestSuite

from sporos.index import ALL_KMERS, MINIMIZER_SEEDS, SYNCMER_SEEDS, SeedHit, SeedIndex, reverse_complement
from reference import canonical_kmer_of, encode_all, naive_minimizers, naive_syncmers, random_dna


def _expected_hits(
    reference: List[UInt8], ref_positions: List[UInt32], query: List[UInt8], query_positions: List[UInt32], k: Int, max_occ: Int
) -> List[SeedHit]:
    var rcodes = encode_all(reference)
    var qcodes = encode_all(query)
    var out = List[SeedHit]()
    for qp in query_positions:
        var q = canonical_kmer_of(qcodes, Int(qp), k)
        var matches = List[Tuple[Int, Bool]]()
        for rp in ref_positions:
            var r = canonical_kmer_of(rcodes, Int(rp), k)
            if r[0] == q[0]:
                matches.append((Int(rp), r[1]))
        if len(matches) == 0 or (max_occ > 0 and len(matches) > max_occ):
            continue
        for m in matches:
            out.append(SeedHit(qp, UInt32(m[0]), UInt32(len(matches)), m[1] != q[1]))
    return out^


def _sorted(var hits: List[SeedHit]) -> List[SeedHit]:
    sort(Span(hits))
    return hits^


def _check(scheme: Int, k: Int, w: Int, n: Int, case_id: Int, all_kmers: Bool, max_occ: Int) raises:
    var reference = random_dna(n, case_id, 0.003)
    # Plant repeats, reverse-complement copies and the reference itself in the query.
    var query = random_dna(n // 2, case_id + 1, 0.003)
    for i in range(min(300, n // 3)):
        query[i] = reference[i + 50]
    var seqs = List[Span[UInt8, origin_of(reference)]]()
    seqs.append(Span(reference))
    var index = SeedIndex(seqs, k, w, scheme, threads=2)
    var ref_positions: List[UInt32]
    if scheme == MINIMIZER_SEEDS:
        ref_positions = naive_minimizers(reference, k, w, True)
    elif scheme == SYNCMER_SEEDS:
        ref_positions = naive_syncmers(reference, k, w, True)
    else:
        ref_positions = List[UInt32]()
        var codes = encode_all(reference)
        for i in range(len(codes) - k + 1):
            var ok = True
            for j in range(k):
                if codes[i + j] > 3:
                    ok = False
            if ok:
                ref_positions.append(UInt32(i))
    var query_positions = index.query_positions(Span(query), all_kmers)
    var expected = _sorted(_expected_hits(reference, ref_positions, query, query_positions, k, max_occ))
    var got = List[SeedHit]()
    index.lookup(Span(query), got, all_kmers, max_occ)
    got = _sorted(got^)
    if len(got) != len(expected):
        raise Error(String(t"scheme={scheme} k={k} w={w} n={n}: got {len(got)} hits, expected {len(expected)}"))
    for i in range(len(got)):
        if not (got[i] == expected[i]):
            raise Error(String(t"scheme={scheme} k={k} w={w}: hit {i} differs"))


def test_reverse_complement_kmer() raises:
    var codes: List[UInt8] = [0, 1, 2, 3, 3, 0, 2]
    var fwd = UInt64(0)
    for i in range(7):
        fwd |= UInt64(codes[i]) << UInt64(2 * i)
    var rc = UInt64(0)
    for i in range(7):
        rc |= UInt64(3 - codes[6 - i]) << UInt64(2 * i)
    assert_equal(reverse_complement(fwd, 7), rc)
    assert_equal(reverse_complement(reverse_complement(fwd, 7), 7), fwd)


def test_seed_index_matches_reference() raises:
    var case_id = 500
    for cfg in [(15, 10), (19, 12), (21, 11), (31, 19), (32, 5), (12, 6)]:
        for n in [200, 3000]:
            case_id += 3
            _check(MINIMIZER_SEEDS, cfg[0], cfg[1], n, case_id, False, 0)
            _check(MINIMIZER_SEEDS, cfg[0], cfg[1], n, case_id, True, 0)
            _check(MINIMIZER_SEEDS, cfg[0], cfg[1], n, case_id, True, 2)
            if cfg[1] < cfg[0]:
                _check(SYNCMER_SEEDS, cfg[0], cfg[1], n, case_id, False, 0)
            _check(ALL_KMERS, cfg[0], cfg[1], n, case_id, False, 3)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
