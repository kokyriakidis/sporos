from std.testing import assert_equal, assert_raises, assert_true, TestSuite

from sporos.alphabet import encode_base, N_CODE
from sporos.minimizers import minimizer_positions
from sporos.syncmers import syncmer_positions
from sporos.nthash import canonical_hash, forward_hash, reverse_hash
from reference import naive_minimizers, naive_syncmers, random_dna, encode_all


def _assert_same(got: List[UInt32], expected: List[UInt32], label: String) raises:
    if len(got) != len(expected):
        var n = min(len(got), len(expected))
        for i in range(n):
            if got[i] != expected[i]:
                raise Error(t"{label}: first mismatch at {i}: got {got[i]} expected {expected[i]}")
        raise Error(t"{label}: length {len(got)} != expected {len(expected)}")
    for i in range(len(got)):
        if got[i] != expected[i]:
            raise Error(t"{label}: mismatch at {i}: got {got[i]} expected {expected[i]}")


def test_encoding() raises:
    assert_equal(encode_base(UInt8(ord("A"))), 0)
    assert_equal(encode_base(UInt8(ord("c"))), 1)
    assert_equal(encode_base(UInt8(ord("G"))), 2)
    assert_equal(encode_base(UInt8(ord("t"))), 3)
    assert_equal(encode_base(UInt8(ord("U"))), 3)
    assert_equal(encode_base(UInt8(ord("N"))), N_CODE)
    assert_equal(encode_base(UInt8(ord("-"))), N_CODE)


def test_canonical_hash_is_strand_symmetric() raises:
    var seq = random_dna(200, 7)
    var codes = encode_all(seq)
    var rc = List[UInt8]()
    for i in reversed(range(len(codes))):
        rc.append(3 - codes[i])
    var k = 21
    for i in range(len(codes) - k + 1):
        var j = len(codes) - k - i
        var a = codes.unsafe_ptr().unsafe_offset(i)
        var b = rc.unsafe_ptr().unsafe_offset(j)
        assert_equal(canonical_hash(a, k), canonical_hash(b, k))
        assert_equal(forward_hash(a, k), reverse_hash(b, k))


def test_minimizers_match_reference() raises:
    var configs = [(5, 1), (11, 5), (15, 10), (21, 11), (31, 19), (40, 7), (64, 3)]
    var lengths = [1, 4, 20, 63, 64, 65, 1000, 4099]
    var case_id = 0
    for cfg in configs:
        for n in lengths:
            for canonical in [True, False]:
                for n_rate in [0.0, 0.01]:
                    case_id += 1
                    var seq = random_dna(n, case_id, n_rate)
                    var got = minimizer_positions(Span(seq), cfg[0], cfg[1], canonical)
                    var expected = naive_minimizers(seq, cfg[0], cfg[1], canonical)
                    _assert_same(
                        got,
                        expected,
                        String(
                            t"minimizers k={cfg[0]} w={cfg[1]} n={n} canonical={canonical} nrate={n_rate}",
                        ),
                    )


def test_minimizers_multithreaded_match_single() raises:
    var seq = random_dna(5_000_000, 99, 0.0001)
    var single = minimizer_positions(Span(seq), 21, 11, True, 1)
    var multi = minimizer_positions(Span(seq), 21, 11, True, 8)
    _assert_same(multi, single, "threads")
    assert_true(len(single) > 400_000)


def test_syncmers_match_reference() raises:
    var configs = [(15, 5), (21, 11), (31, 15), (12, 12), (8, 1)]
    var lengths = [10, 63, 64, 65, 2000]
    var case_id = 1000
    for cfg in configs:
        for n in lengths:
            for canonical in [True, False]:
                case_id += 1
                var seq = random_dna(n, case_id, 0.01)
                var got = syncmer_positions(Span(seq), cfg[0], cfg[1], canonical)
                var expected = naive_syncmers(seq, cfg[0], cfg[1], canonical)
                _assert_same(got, expected, String(t"closed syncmers k={cfg[0]} s={cfg[1]} n={n}"))
                var t = (cfg[0] - cfg[1]) // 2
                var got_open = syncmer_positions(Span(seq), cfg[0], cfg[1], canonical, t)
                var expected_open = naive_syncmers(seq, cfg[0], cfg[1], canonical, t)
                _assert_same(
                    got_open,
                    expected_open,
                    String(
                        t"open syncmers k={cfg[0]} s={cfg[1]} t={t} n={n}",
                    ),
                )


def test_open_syncmers_with_several_offsets() raises:
    var configs = [(15, 5), (21, 11), (20, 16), (40, 5)]
    var lengths = [10, 65, 3000]
    var case_id = 3000
    for cfg in configs:
        var k = cfg[0]
        var s = cfg[1]
        var sets: List[List[Int]] = [[0, k - s], [1, (k - s) // 2, k - s - 1], [k - s]]
        if k - s >= 32:
            sets.append([3, 31, 32, 33])
        for n in lengths:
            for canonical in [True, False]:
                for ref offsets in sets:
                    case_id += 1
                    var seq = random_dna(n, case_id, 0.01)
                    var got = syncmer_positions(Span(seq), k, s, canonical, open_offsets=offsets)
                    var expected = naive_syncmers(seq, k, s, canonical, open_offsets=offsets)
                    _assert_same(got, expected, String(t"open syncmers k={k} s={s} n={n} offsets={len(offsets)}"))
    # Offsets {0, k - s} pick exactly the closed syncmers.
    var seq = random_dna(10_000, 77, 0.01)
    _assert_same(
        syncmer_positions(Span(seq), 21, 11, open_offsets=[0, 10]),
        syncmer_positions(Span(seq), 21, 11),
        "offsets {0, k - s} vs closed",
    )
    with assert_raises():
        _ = syncmer_positions(Span(seq), 21, 11, open_offsets=[2, 11])
    with assert_raises():
        _ = syncmer_positions(Span(seq), 21, 11, True, 2, open_offsets=[3])


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
