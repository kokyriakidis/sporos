"""The `sporos.yak` kernels against oracles: yak's hash against values its C code computes, compressed k-mers and their counts
against brute force, the minimizer sketch against digests of hifiasm's own sketch (the hifiasmCandidates fork's
`hifiasm_sketch_minimizers`, unfiltered) on generated sequences, the position table against a brute-force map.
hifiasm's filtered sketch with its re-seeding is held to the C code on real reads by hifiasmCandidatesMojo."""

from std.testing import assert_equal, assert_true, TestSuite

from sporos.fasta import FastaRecords
from sporos.yak import (
    BloomFilters,
    FilterTable,
    KmerCounts,
    PositionTable,
    Sketch,
    analyze_count,
    hpc_kmer_hashes,
    nt4,
    position_pos,
    position_rid,
    sketch_all,
    sketch_read,
    yak_hash64_64,
    yak_hash64_64_lanes,
)
from reference import naive_hpc_kmer_planes, yak_test_sequence


def _records(seqs: List[List[UInt8]]) -> FastaRecords:
    var records = FastaRecords()
    for i in range(len(seqs)):
        records.names.append(String(t"r{i}"))
        records.starts.append(len(records.data))
        records.lengths.append(len(seqs[i]))
        records.data.extend(Span(seqs[i]))
    return records^


def _naive_hashes(seq: List[UInt8], k: Int) -> List[UInt64]:
    var out = List[UInt64]()
    for p in naive_hpc_kmer_planes(seq, k):
        if p[1] < p[3]:
            out.append(yak_hash64_64(p[0]) + yak_hash64_64(p[1]))
        else:
            out.append(yak_hash64_64(p[2]) + yak_hash64_64(p[3]))
    return out^


def test_yak_hash_matches_c() raises:
    assert_equal(yak_hash64_64(0), UInt64(8633297058295171728))
    assert_equal(yak_hash64_64(1), UInt64(6614235796240398542))
    assert_equal(yak_hash64_64(123456789), UInt64(16581954974024456952))
    assert_equal(yak_hash64_64(UInt64.MAX), UInt64(2272383144869939092))
    var keys = SIMD[DType.uint64, 8](0, 1, 123456789, UInt64.MAX, 5, 6, 7, 8)
    var lanes = yak_hash64_64_lanes(keys)
    for i in range(8):
        assert_equal(lanes[i], yak_hash64_64(keys[i]))


def test_nt4() raises:
    var expected: List[UInt8] = [0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4]
    var letters = "AaCcGgTtUuNn-"
    for i in range(len(expected)):
        assert_equal(nt4(letters.as_bytes()[i]), expected[i])
    for c in range(256):
        var code = nt4(UInt8(c))
        assert_true(code <= 4)


def test_hpc_kmer_hashes_against_brute_force() raises:
    var ks: List[Int] = [15, 21, 31, 51]
    for index in range(4):
        var seq = yak_test_sequence(3000 + 700 * index, 11 + index)
        var k = ks[index]
        var got = List[UInt64]()
        hpc_kmer_hashes(Span(seq), k, got)
        var expected = _naive_hashes(seq, k)
        assert_equal(len(got), len(expected))
        for i in range(len(got)):
            assert_equal(got[i], expected[i])


def test_counts_against_brute_force() raises:
    var base = yak_test_sequence(3000, 21)
    var seqs = List[List[UInt8]]()
    for _ in range(3):
        seqs.append(base.copy())
    seqs.append(yak_test_sequence(2000, 22))
    var records = _records(seqs)
    var expected = Dict[UInt64, Int]()
    for s in seqs:
        for h in _naive_hashes(s, 21):
            expected[h] = expected.get(h, 0) + 1
    # Exact (two blocks), then through Bloom filters: these count from the second sighting, so a k-mer seen once
    # is absent and the rest count the same (the filters are large enough for no false positive here).
    var exact = KmerCounts(records, 21, 3, 6000, 0)
    var bloomed = KmerCounts(records, 21, 3, 6000, 26)
    assert_equal(exact.distinct(), len(expected))
    var repeated = 0
    for entry in expected.items():
        assert_equal(exact.count(entry.key), entry.value)
        assert_equal(bloomed.count(entry.key), entry.value if entry.value >= 2 else 0)
        if entry.value >= 2:
            repeated += 1
    assert_equal(bloomed.distinct(), repeated)
    var hist = exact.histogram()
    var threes = 0
    for entry in expected.items():
        if entry.value == 3:
            threes += 1
    assert_equal(hist[3], threes)


def test_bloom_insert() raises:
    var bf = BloomFilters(2, 16)
    assert_equal(bf.insert(0, 12345), 0)
    assert_equal(bf.insert(0, 12345), 4)
    assert_equal(bf.insert(1, 12345), 0)


def test_analyze_count_peaks() raises:
    var cnt = List[Int](length=4096, fill=0)
    cnt[1] = 100000
    cnt[2] = 20000
    cnt[3] = 3000
    cnt[4] = 600
    for i in range(5, 200):
        cnt[i] = 100
    cnt[14] = 2500
    cnt[15] = 3000
    cnt[16] = 2600
    cnt[29] = 8000
    cnt[30] = 9000
    cnt[31] = 8500
    var peaks = analyze_count(cnt, 5)
    assert_equal(peaks.hom, 30)
    assert_equal(peaks.het, 15)
    var flat = List[Int](length=4096, fill=0)
    flat[1] = 10
    assert_equal(analyze_count(flat, 5).hom, -1)


def test_filter_table() raises:
    var base = yak_test_sequence(500, 31)
    var seqs = List[List[UInt8]]()
    for _ in range(7):
        seqs.append(base.copy())
    seqs.append(yak_test_sequence(400, 32))
    var counts = KmerCounts(_records(seqs), 15, 2, 100000, 0)
    var hashes = _naive_hashes(base, 15)
    assert_equal(FilterTable(counts, 5).count(hashes[0]), 7)
    assert_equal(FilterTable(counts, 5).count(_naive_hashes(seqs[7], 15)[0]), 0)
    assert_equal(FilterTable(counts, 5, 6).count(hashes[0]), 2147483647)
    assert_equal(FilterTable().count(hashes[0]), 0)


def _digest(sketch: Sketch) -> Tuple[Int, UInt64]:
    """As the C program digests `hifiasm_sketch_minimizers`: START positions, sorted and deduplicated, each
    `(start, span, rev, hash)` folded into FNV-1a."""
    var keyed = List[Tuple[Int, Int]]()
    for i in range(len(sketch.mz)):
        keyed.append((sketch.mz[i].start(), i))

    def less(a: Tuple[Int, Int], b: Tuple[Int, Int]) -> Bool:
        return a[0] < b[0] or (a[0] == b[0] and a[1] < b[1])

    sort(Span(keyed), less)
    var h: UInt64 = 1469598103934665603
    var count = 0
    var last = -1
    for entry in keyed:
        if entry[0] == last:
            continue
        last = entry[0]
        var m = sketch.mz[entry[1]]
        var values = [UInt64(entry[0]), UInt64(m.span), UInt64(m.rev), m.x]
        for v in values:
            h ^= v
            h *= 1099511628211
        count += 1
    return (count, h)


def test_sketch_matches_hifiasm() raises:
    # (length, seed, k, w, hpc, count, digest), from hifiasm's `hifiasm_sketch_minimizers` on the same sequences.
    var cases: List[Tuple[Int, Int, Int, Int, Bool, Int, UInt64]] = [
        (5000, 1, 15, 10, True, 543, UInt64(8472492873556547778)),
        (20000, 2, 51, 51, True, 359, UInt64(10156247427611996938)),
        (20000, 3, 31, 19, False, 1843, UInt64(6277190111758445410)),
        (3000, 4, 63, 5, True, 524, UInt64(9299979190830843069)),
        (12000, 5, 21, 11, True, 1211, UInt64(9410452428482244271)),
    ]
    var none = FilterTable()
    for c in cases:
        var seq = yak_test_sequence(c[0], c[1])
        var out = Sketch()
        sketch_read(Span(seq), c[3], c[2], c[4], none, 0, 1000, out)
        var d = _digest(out)
        assert_equal(d[0], c[5])
        assert_equal(d[1], c[6])


def test_position_table_against_brute_force() raises:
    var seqs = List[List[UInt8]]()
    var shared = yak_test_sequence(4000, 41)
    for r in range(6):
        var s = shared.copy()
        s.extend(Span(yak_test_sequence(1000, 50 + r)))
        seqs.append(s^)
    var records = _records(seqs)
    var sketches = sketch_all(records, FilterTable(), 11, 15, 2, 0, 1000)
    var table = PositionTable(sketches, 2)
    var occurrences = Dict[UInt64, List[Tuple[Int, Int]]]()
    for r in range(len(sketches)):
        for m in sketches[r].mz:
            if m.x not in occurrences:
                occurrences[m.x] = List[Tuple[Int, Int]]()
            occurrences[m.x].append((Int(m.rid), Int(m.pos)))
    var kept = 0
    var positions = 0
    for entry in occurrences.items():
        var found = table.get(entry.key)
        var n = len(entry.value)
        if n < 2 or n > 4094:
            assert_equal(found[1], 0)
            continue
        kept += 1
        positions += n
        assert_equal(found[1], n)
        for i in range(n):
            var p = table.position(entry.key, found[0] + i)
            assert_equal(position_rid(p), entry.value[i][0])
            assert_equal(position_pos(p), entry.value[i][1])
    assert_true(kept > 100)
    assert_equal(table.distinct, kept)
    assert_equal(table.indexed, positions)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
