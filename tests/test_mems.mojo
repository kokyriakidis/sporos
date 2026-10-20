from std.testing import assert_equal, assert_true, TestSuite

from sporos.mems import Mem, MemIndex, MemParams, filter_mums, find_mems
from sporos.text import SeqText
from reference import naive_count, naive_mems, random_dna


def _mutate(src: List[UInt8], rate_per_million: Int, rng_seed: Int) -> List[UInt8]:
    """Copy of `src` with substitutions, plus a duplicated block so repeats exist."""
    var out = src.copy()
    var state = UInt64(rng_seed) * 0x9E3779B97F4A7C15 + 1
    for i in range(len(out)):
        state ^= state >> 12
        state ^= state << 25
        state ^= state >> 27
        var r = state * 0x2545F4914F6CDD1D
        if Int(r % 1_000_000) < rate_per_million:
            out[i] = UInt8(ord("ACGT"[byte=Int((r >> 40) % 4)]))
    if len(out) > 400:
        for i in range(150):
            out[len(out) - 200 + i] = out[100 + i]
    return out^


def _check_mems(reference: List[UInt8], query: List[UInt8], params: MemParams, threads: Int) raises:
    var seqs = List[Span[UInt8, origin_of(reference)]]()
    seqs.append(Span(reference))
    var rtext = SeqText(seqs, is_query=False)
    var qseqs = List[Span[UInt8, origin_of(query)]]()
    qseqs.append(Span(query))
    var qtext = SeqText(qseqs, is_query=True)
    var index = MemIndex(rtext, params)
    var got = find_mems(index, rtext, qtext, threads)
    var expected = naive_mems(reference, query, params.min_len)
    if len(got) != len(expected):
        raise Error(String(t"{params}: got {len(got)} MEMs, expected {len(expected)}"))
    for i in range(len(got)):
        var e = expected[i]
        assert_equal(Int(got[i].ref_pos), e[0])
        assert_equal(Int(got[i].query_pos), e[1])
        assert_equal(Int(got[i].length), e[2])

    var mums = filter_mums(got, qtext)
    var expected_mums = 0
    for m in got:
        var a = naive_count(reference, Int(m.ref_pos), Int(m.length), reference)
        var b = naive_count(query, Int(m.query_pos), Int(m.length), query)
        if a == 1 and b == 1:
            expected_mums += 1
    assert_equal(len(mums), expected_mums)
    for m in mums:
        assert_equal(naive_count(reference, Int(m.ref_pos), Int(m.length), reference), 1)
        assert_equal(naive_count(query, Int(m.query_pos), Int(m.length), query), 1)


def test_mems_match_reference() raises:
    var case_id = 0
    for n in [50, 300, 1000]:
        for min_len in [12, 20, 33, 60]:
            for n_rate in [0.0, 0.002]:
                case_id += 1
                var reference = random_dna(n, case_id, n_rate)
                var query = _mutate(reference, 30_000, case_id)
                var params = MemParams.choose(min_len, n, n)
                _check_mems(reference, query, params, 1)
                _check_mems(reference, query, params, 4)
                # Every coprime sampling must give the same answer.
                var k = min(min_len, 12)
                _check_mems(reference, query, MemParams(min_len, k, 1, 1), 1)
                if min_len - k + 1 >= 6:
                    _check_mems(reference, query, MemParams(min_len, k, 3, 2), 1)
                    _check_mems(reference, query, MemParams(min_len, k, 2, 3), 1)


def test_low_complexity_repeats() raises:
    var reference = List[UInt8]()
    for i in range(600):
        reference.append(UInt8(ord("ACACAGT"[byte=i % 7])))
    var query = List[UInt8]()
    for i in range(500):
        query.append(UInt8(ord("ACACAGT"[byte=(i + 3) % 7])))
    _check_mems(reference, query, MemParams.choose(20, 600, 500), 1)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
