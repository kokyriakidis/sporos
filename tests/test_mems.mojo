from std.testing import assert_equal, assert_true, TestSuite

from sporos.mems import Mem, MemIndex, MemParams, filter_mums, find_mems
from sporos.sort import parallel_sort
from sporos.text import SeqText
from sporos.alphabet import encode_base
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
    var rtext = SeqText(seqs)
    var qseqs = List[Span[UInt8, origin_of(query)]]()
    qseqs.append(Span(query))
    var qtext = SeqText(qseqs)
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


def test_reverse_complement_query() raises:
    var reference = random_dna(1500, 77, 0.002)
    var query = _mutate(reference, 20_000, 78)
    var rc = List[UInt8](capacity=len(query))
    for i in reversed(range(len(query))):
        var code = encode_base(query[i])
        rc.append(UInt8(ord("TGCAN"[byte=Int(min(code, 4))])))
    var seqs = List[Span[UInt8, origin_of(reference)]]()
    seqs.append(Span(reference))
    var rtext = SeqText(seqs)
    var qseqs = List[Span[UInt8, origin_of(query)]]()
    qseqs.append(Span(query))
    var qtext = SeqText(qseqs, reverse_complement=True)
    var params = MemParams.choose(20, len(reference), len(query))
    var got = find_mems(MemIndex(rtext, params), rtext, qtext, 2)
    var expected = naive_mems(reference, rc, 20)
    assert_equal(len(got), len(expected))
    for i in range(len(got)):
        assert_equal(Int(got[i].ref_pos), expected[i][0])
        assert_equal(Int(got[i].query_pos), expected[i][1])
        assert_equal(Int(got[i].length), expected[i][2])


def test_low_complexity_repeats() raises:
    var reference = List[UInt8]()
    for i in range(600):
        reference.append(UInt8(ord("ACACAGT"[byte=i % 7])))
    var query = List[UInt8]()
    for i in range(500):
        query.append(UInt8(ord("ACACAGT"[byte=(i + 3) % 7])))
    _check_mems(reference, query, MemParams.choose(20, 600, 500), 1)


def test_parallel_sort() raises:
    var state = UInt64(12345)
    for n in [0, 1, 1000, 70_000, 1_000_000]:
        for bits in [8, 40, 64]:
            var keys = List[UInt64](capacity=n)
            for _ in range(n):
                state ^= state << 13
                state ^= state >> 7
                state ^= state << 17
                keys.append(state >> UInt64(64 - bits))
            var expected = keys.copy()
            sort(Span(expected))
            for threads in [1, 4]:
                var got = keys.copy()
                parallel_sort(got, threads)
                assert_equal(len(got), n)
                for i in range(n):
                    if got[i] != expected[i]:
                        raise Error(String(t"n={n} bits={bits} threads={threads}: key {i} differs"))
                # By the high bits only: those in order, and the same keys.
                var partial = keys.copy()
                parallel_sort(partial, threads, low_bits=bits // 2)
                for i in range(1, n):
                    if partial[i - 1] >> UInt64(bits // 2) > partial[i] >> UInt64(bits // 2):
                        raise Error(String(t"n={n} bits={bits} threads={threads}: high bits out of order at {i}"))
                sort(Span(partial))
                for i in range(n):
                    if partial[i] != expected[i]:
                        raise Error(String(t"n={n} bits={bits} threads={threads}: partial sort lost key {i}"))


def test_mum_filter_threads_agree() raises:
    # Sixty noisy copies of one unit: enough MEMs for the parallel sorts.
    var unit = random_dna(2500, 31)
    var reference = List[UInt8]()
    for copy in range(60):
        reference.extend(Span(_mutate(unit, 20_000, 100 + copy)))
    var query = _mutate(reference, 20_000, 7)
    var rseqs = List[Span[UInt8, origin_of(reference)]]()
    rseqs.append(Span(reference))
    var rtext = SeqText(rseqs)
    var qseqs = List[Span[UInt8, origin_of(query)]]()
    qseqs.append(Span(query))
    var qtext = SeqText(qseqs)
    var params = MemParams.choose(20, len(reference), len(query))
    var index = MemIndex(rtext, params)
    var mems = find_mems(index, rtext, qtext, 4)
    if len(mems) <= 70_000:
        raise Error(String(t"only {len(mems)} MEMs"))
    for unique_in_query in [True, False]:
        var serial = filter_mums(mems, qtext, unique_in_query, threads=1)
        var parallel = filter_mums(mems, qtext, unique_in_query, threads=4)
        assert_equal(len(parallel), len(serial))
        for i in range(len(serial)):
            assert_true(parallel[i] == serial[i])


def test_params_for_long_min_len() raises:
    # The step search is linear in L - k: long minimum lengths stay instant.
    var params = MemParams.choose(10_000_000, 1000, 1000)
    params.validate()
    var short = MemParams.choose(50, 1_000_000, 1_000_000)
    short.validate()
    assert_equal(short.ref_step * short.query_step <= 50 - short.k + 1, True)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
