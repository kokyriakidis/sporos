from std.collections import Set
from std.testing import assert_equal, assert_true, TestSuite

from sporos.mems import Mem, MemIndex, MemParams, SmemHit, filter_by_occurrences, filter_mums, find_mems, find_smems
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


def _reverse_complement(seq: List[UInt8]) -> List[UInt8]:
    var out = List[UInt8](capacity=len(seq))
    for i in reversed(range(len(seq))):
        var c = seq[i] | 0x20
        var rc = UInt8(ord("N"))
        if c == UInt8(ord("a")):
            rc = UInt8(ord("T"))
        elif c == UInt8(ord("c")):
            rc = UInt8(ord("G"))
        elif c == UInt8(ord("g")):
            rc = UInt8(ord("C"))
        elif c == UInt8(ord("t")):
            rc = UInt8(ord("A"))
        out.append(rc)
    return out^


def _naive_smems(reference: List[UInt8], query: List[UInt8], min_len: Int) -> List[Tuple[Int, Int, Int, Bool]]:
    """(start, end, ref_pos, reverse) for every hit of every SMEM, by definition:
    MEM intervals on both strands that no other MEM interval strictly contains."""
    var all = List[Tuple[Int, Int, Int, Bool]]()
    for m in naive_mems(reference, query, min_len):
        all.append((m[1], m[1] + m[2], m[0], False))
    var n = len(query)
    for m in naive_mems(reference, _reverse_complement(query), min_len):
        all.append((n - m[1] - m[2], n - m[1], m[0], True))
    var out = List[Tuple[Int, Int, Int, Bool]]()
    for a in all:
        var contained = False
        for b in all:
            if b[0] <= a[0] and a[1] <= b[1] and (b[0] < a[0] or a[1] < b[1]):
                contained = True
                break
        if not contained:
            out.append(a)
    return out^


def test_smems_match_definition() raises:
    var case_id = 4000
    for n in [300, 1500]:
        for min_len in [12, 20, 31]:
            for threads in [1, 3]:
                case_id += 1
                # Two reference sequences with repeats, a reverse-complement copy and N runs;
                # two query sequences of their pieces.
                var ref1 = _mutate(random_dna(n, case_id, 0.003), 30_000, case_id)
                var ref2 = random_dna(n // 2, case_id + 1, 0.003)
                var piece = List[UInt8](Span(ref1)[50 : 50 + n // 4])
                ref2.extend(Span(_reverse_complement(piece)))
                var reference = ref1.copy()
                reference.append(UInt8(ord("N")))  # the separator SeqText puts between sequences
                reference.extend(Span(ref2))
                var q1 = _mutate(List[UInt8](Span(ref1)[n // 3 :]), 20_000, case_id + 2)
                var q2 = _reverse_complement(List[UInt8](Span(ref2)[10 : n // 3]))
                # Sequences as slices of one buffer each (one origin per list).
                var ref_buffer = ref1.copy()
                ref_buffer.extend(Span(ref2))
                var rspans = List[Span[UInt8, origin_of(ref_buffer)]]()
                rspans.append(Span(ref_buffer)[: len(ref1)])
                rspans.append(Span(ref_buffer)[len(ref1) :])
                var rtext = SeqText(rspans)
                var queries = List[List[UInt8]]()
                queries.append(q1.copy())
                queries.append(q2.copy())
                var query_buffer = q1.copy()
                query_buffer.extend(Span(q2))
                var qspans = List[Span[UInt8, origin_of(query_buffer)]]()
                qspans.append(Span(query_buffer)[: len(q1)])
                qspans.append(Span(query_buffer)[len(q1) :])
                var qtext = SeqText(qspans)
                var qrc = SeqText(qspans, reverse_complement=True)
                var index = MemIndex(rtext, MemParams.choose(min_len, len(rtext), len(qtext)), threads)
                var got = find_smems(index, rtext, qtext, qrc, threads)
                var flat = List[Tuple[Int, Int, Int, Bool]]()
                for i in range(len(got)):
                    var smem = got.smems[i]
                    for h in got.hits_of(i):
                        flat.append((Int(smem.start), Int(smem.end), Int(h.ref_pos), h.reverse))
                var expected = List[Tuple[Int, Int, Int, Bool]]()
                for qi in range(len(queries)):
                    for e in _naive_smems(reference, queries[qi], min_len):
                        var offset = qtext.starts[qi]
                        expected.append((e[0] + offset, e[1] + offset, e[2], e[3]))
                if len(flat) != len(expected):
                    raise Error(String(t"n={n} L={min_len}: {len(flat)} SMEM hits, expected {len(expected)}"))
                # Same multiset: compare as sorted keys.
                var a = List[Int]()
                var b = List[Int]()
                for x in flat:
                    a.append(((x[0] * 4096 + (x[1] - x[0])) * 65536 + x[2]) * 2 + Int(x[3]))
                for x in expected:
                    b.append(((x[0] * 4096 + (x[1] - x[0])) * 65536 + x[2]) * 2 + Int(x[3]))
                sort(Span(a))
                sort(Span(b))
                for i in range(len(a)):
                    if a[i] != b[i]:
                        raise Error(String(t"n={n} L={min_len} threads={threads}: SMEM hit {i} differs"))
                assert_true(len(flat) > 0)


def test_occurrence_filter_parallel_matches_serial() raises:
    # Enough MEMs in one query sequence for the block-parallel sweeps, with
    # large equal-start groups that straddle block bounds; a sample is also
    # checked against a direct count of containing intervals.
    var n = 150_000
    var seq = random_dna(2 * n + 1000, 77)  # longer than every synthetic MEM reaches
    var spans = List[Span[UInt8, origin_of(seq)]]()
    spans.append(Span(seq))
    var qtext = SeqText(spans)
    var state = UInt64(12345)
    var mems = List[Mem](capacity=n)
    var q = 0
    while len(mems) < n:
        state ^= state >> 12
        state ^= state << 25
        state ^= state >> 27
        var r = state * 0x2545F4914F6CDD1D
        var group = 1 + Int((r >> 8) % 40) if (r >> 50) % 16 == 0 else 1
        for g in range(min(group, n - len(mems))):
            var length = 20 + Int((r >> UInt64(13 + g % 30)) % 180)
            mems.append(Mem(UInt32(Int((r >> 20) % 1_000_000) + 7 * g), UInt32(q), UInt32(length)))
        q += Int((r >> 40) % 4)
    for limits in [(2, 0, 1), (5, 0, 1), (1_000_000, 0, 1), (0, 0, 2), (3, 3, 1), (0, 2, 1), (10, 0, 3)]:
        var serial = filter_by_occurrences(mems.copy(), qtext, limits[0], limits[1], limits[2], threads=1)
        var parallel = filter_by_occurrences(mems.copy(), qtext, limits[0], limits[1], limits[2], threads=4)
        assert_equal(len(serial), len(parallel))
        for i in range(len(serial)):
            assert_true(serial[i] == parallel[i])
        if limits[1] == 0:
            # Reference copies: MEMs whose query interval contains (or equals) this one's.
            var kept = Set[Int]()
            for m in serial:
                kept.add((Int(m.ref_pos) << 32) | (Int(m.query_pos) << 8) ^ Int(m.length))
            for s in range(0, n, 997):
                var a = mems[s]
                var copies = 0
                for b in mems:
                    if b.query_pos <= a.query_pos and b.query_pos + b.length >= a.query_pos + a.length:
                        copies += 1
                var keep = (limits[0] == 0 or copies <= limits[0]) and copies >= limits[2]
                var key = (Int(a.ref_pos) << 32) | (Int(a.query_pos) << 8) ^ Int(a.length)
                assert_equal(keep, key in kept)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
