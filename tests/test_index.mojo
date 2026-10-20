from std.testing import assert_equal, assert_true, TestSuite

from sporos.index import (
    ALL_KMERS,
    MINIMIZER_SEEDS,
    SYNCMER_SEEDS,
    SeedHit,
    SeedIndex,
    fingerprint,
    kmer_hash,
    kmer_mask,
    reverse_complement,
)
from reference import canonical_kmer_of, encode_all, naive_minimizers, naive_syncmers, random_dna


def _expected_hits(
    reference: List[UInt8],
    ref_positions: List[UInt32],
    query: List[UInt8],
    query_positions: List[UInt32],
    k: Int,
    max_occ: Int,
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


def _check(
    scheme: Int,
    k: Int,
    w: Int,
    n: Int,
    case_id: Int,
    all_kmers: Bool,
    max_occ: Int,
    open_offset: Int = -1,
    *,
    open_offsets: List[Int] = [],
) raises:
    var reference = random_dna(n, case_id, 0.003)
    # Plant repeats, reverse-complement copies and the reference itself in the query.
    var query = random_dna(n // 2, case_id + 1, 0.003)
    for i in range(min(300, n // 3)):
        query[i] = reference[i + 50]
    var seqs = List[Span[UInt8, origin_of(reference)]]()
    seqs.append(Span(reference))
    var opened: Optional[Int] = None
    if open_offset >= 0:
        opened = open_offset
    var index = SeedIndex(seqs, k, w, scheme, threads=2, open_offset=opened, open_offsets=open_offsets)
    var ref_positions: List[UInt32]
    if scheme == MINIMIZER_SEEDS:
        ref_positions = naive_minimizers(reference, k, w, True)
    elif scheme == SYNCMER_SEEDS:
        ref_positions = naive_syncmers(reference, k, w, True, open_offset, open_offsets=open_offsets)
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
    # The query's seeds, from the oracles (or every valid k-mer).
    var query_positions: List[UInt32]
    if scheme == MINIMIZER_SEEDS and not all_kmers:
        query_positions = naive_minimizers(query, k, w, True)
    elif scheme == SYNCMER_SEEDS and not all_kmers:
        query_positions = naive_syncmers(query, k, w, True, open_offset, open_offsets=open_offsets)
    else:
        query_positions = List[UInt32]()
        var codes = encode_all(query)
        for i in range(len(codes) - k + 1):
            var ok = True
            for j in range(k):
                if codes[i + j] > 3:
                    ok = False
            if ok:
                query_positions.append(UInt32(i))
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
                _check(SYNCMER_SEEDS, cfg[0], cfg[1], n, case_id, False, 0, (cfg[0] - cfg[1]) // 2)
                _check(SYNCMER_SEEDS, cfg[0], cfg[1], n, case_id, False, 0, open_offsets=[0, cfg[0] - cfg[1] - 1])
            _check(ALL_KMERS, cfg[0], cfg[1], n, case_id, False, 3)


def test_lookup_batch_matches_lookup() raises:
    var reference = random_dna(60_000, 77, 0.002)
    # Repeats, so some k-mers pass and some fail an occurrence cap.
    for copy in range(6):
        for i in range(500):
            reference[20_000 + copy * 2000 + i] = reference[1000 + i]
    var seqs = List[Span[UInt8, origin_of(reference)]]()
    seqs.append(Span(reference))
    # Reads of varied lengths (some shorter than k, some long), from the
    # reference and its noisy copy, enough to span several lookup tasks.
    var noisy = random_dna(60_000, 78, 0.01)
    for i in range(0, 60_000, 3):
        noisy[i] = reference[i]
    var reads = List[Span[UInt8, origin_of(noisy)]]()
    var at = 0
    for r in range(5000):
        var length = [5, 40, 150, 151, 999][r % 5] if r % 97 != 0 else 0
        if r % 500 == 1:
            length = 20_000
        if at + length > len(noisy):
            at = (r * 7919) % 1000
        reads.append(Span(noisy)[at : at + length])
        at += 37 + length // 3
    for scheme in [MINIMIZER_SEEDS, SYNCMER_SEEDS, ALL_KMERS]:
        var index = SeedIndex(seqs, 19, 12 if scheme != SYNCMER_SEEDS else 15, scheme, threads=2)
        for all_kmers in [False, True]:
            for max_occ in [0, 3]:
                for threads in [1, 3]:
                    var batch = index.lookup_batch(reads, all_kmers, max_occ, threads)
                    assert_equal(len(batch), len(reads))
                    var hits = List[SeedHit]()
                    for r in range(len(reads)):
                        hits.clear()
                        index.lookup(reads[r], hits, all_kmers, max_occ)
                        var got = batch.of(r)
                        if len(got) != len(hits):
                            raise Error(String(t"scheme={scheme} read {r}: {len(got)} hits, expected {len(hits)}"))
                        for i in range(len(hits)):
                            if not (got[i] == hits[i]):
                                raise Error(String(t"scheme={scheme} read {r}: hit {i} differs"))


def test_bucket_and_fingerprint_identify_kmers() raises:
    # With bucket bits + 31 >= 2k, (bucket, fingerprint) must be unique per k-mer.
    var state = UInt64(99)
    for k in [5, 12, 16, 19, 25, 31, 32]:
        var bucket_bits = max(1, 2 * k - 31)
        var keys = List[UInt64]()
        var seen = List[UInt64]()
        for _ in range(100_000):
            state ^= state << 13
            state ^= state >> 7
            state ^= state << 17
            var kmer = state & kmer_mask(k)
            var h = kmer_hash(kmer, k)
            seen.append(kmer)
            keys.append(((h >> UInt64(64 - bucket_bits)) << 31) | fingerprint(h, k))
        sort(Span(seen))
        sort(Span(keys))
        var distinct_kmers = 0
        var distinct_keys = 0
        for i in range(len(seen)):
            distinct_kmers += Int(i == 0 or seen[i] != seen[i - 1])
            distinct_keys += Int(i == 0 or keys[i] != keys[i - 1])
        assert_equal(distinct_keys, distinct_kmers)


def test_no_hits_for_absent_kmers() raises:
    # These 16-mers shared a bucket and fingerprint under the previous hash.
    var reference = List[UInt8]("CCGCGCGGACAGCGAG".as_bytes())
    var query = List[UInt8]("CAGGCGTTAATGAATA".as_bytes())
    var seqs = List[Span[UInt8, origin_of(reference)]]()
    seqs.append(Span(reference))
    var index = SeedIndex(seqs, 16, 1, ALL_KMERS)
    var hits = List[SeedHit]()
    index.lookup(Span(query), hits, True)
    assert_equal(len(hits), 0)
    index.lookup(Span(reference), hits, True)
    assert_equal(len(hits), 1)


def test_occurrence_cutoff() raises:
    # Ten copies of one 40-mer among random sequence: its k-mers occur 10x.
    var reference = random_dna(5000, 901)
    for copy in range(10):
        for i in range(40):
            reference[1000 + copy * 300 + i] = reference[100 + i]
    var seqs = List[Span[UInt8, origin_of(reference)]]()
    seqs.append(Span(reference))
    var index = SeedIndex(seqs, 15, 1, ALL_KMERS)
    assert_equal(index.occurrence_cutoff(0.5), 1)
    assert_true(index.occurrence_cutoff(1.0) >= 10)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
