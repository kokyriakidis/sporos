from std.testing import assert_equal, assert_true, TestSuite

from sporos.fmindex import FMIndex, ReadSmems
from sporos.mems import MemIndex, MemParams, find_smems
from sporos.text import SeqText
from reference import random_dna


def _next(mut state: UInt64) -> UInt64:
    state ^= state >> 12
    state ^= state << 25
    state ^= state >> 27
    return state * 0x2545F4914F6CDD1D


def _reverse_complement(seq: Span[UInt8, _]) -> List[UInt8]:
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


def _reads(reference: List[UInt8], count: Int, seed: Int) -> List[List[UInt8]]:
    """Reads cut from `reference` (either strand, with substitutions and the
    odd N), random reads, and reads shorter than any seed."""
    var state = UInt64(seed) * 0x9E3779B97F4A7C15 + 7
    var reads = List[List[UInt8]]()
    for r in range(count):
        var length = 20 + Int(_next(state) % 160)
        if r % 7 == 6:
            reads.append(random_dna(length, seed + r))
            continue
        if r % 11 == 10:
            reads.append(random_dna(Int(_next(state) % 6), seed + r))
            continue
        var at = Int(_next(state) % UInt64(max(1, len(reference) - length)))
        var read = List[UInt8](Span(reference)[at : min(len(reference), at + length)])
        for i in range(len(read)):
            var x = _next(state) % 1000
            if x < 8:
                read[i] = UInt8(ord("ACGT"[byte=Int(x % 4)]))
            elif x == 8:
                read[i] = UInt8(ord("N"))
        if r % 2 == 1:
            read = _reverse_complement(Span(read))
        reads.append(read^)
    return reads^


def _check(
    reference: List[UInt8], cuts: List[Int], queries: List[List[UInt8]], threads: Int, save: Bool, sample_bits: Int = 0
) raises:
    """FM-index SMEMs equal the MEM engine's (`find_smems`, itself held to a
    brute-force oracle) for every query and several lengths."""
    var rspans = List[Span[UInt8, origin_of(reference)]]()
    var names = List[String]()
    for i in range(len(cuts) - 1):
        rspans.append(Span(reference)[cuts[i] : cuts[i + 1]])
        names.append(String(t"r{i}"))
    var rtext = SeqText(rspans)
    var fm = FMIndex.build(rspans, names, threads, sample_bits)
    if save:
        fm.save("build/test_fmindex.sfm")
        fm = FMIndex.load("build/test_fmindex.sfm", threads)
        assert_equal(fm.names[len(fm.names) - 1], names[len(names) - 1])
    var buffer = List[UInt8]()
    var bounds: List[Int] = [0]
    for q in queries:
        buffer.extend(Span(q))
        bounds.append(len(buffer))
    var qspans = List[Span[UInt8, origin_of(buffer)]]()
    for i in range(len(queries)):
        qspans.append(Span(buffer)[bounds[i] : bounds[i + 1]])
    for min_len in [3, 12, 19, 40]:
        var got = fm.smems(qspans, min_len, max_hits=1 << 30)
        var capped = fm.smems(qspans, min_len, max_hits=2)
        # With no hits wanted, nothing is located: unique matches extend in the index.
        var unlocated = fm.smems(qspans, min_len, max_hits=0)
        assert_equal(len(unlocated.hits), 0)
        assert_equal(len(unlocated.smems), len(got.smems))
        for i in range(len(got.smems)):
            assert_true(unlocated.smems[i] == got.smems[i])
        assert_equal(len(got), len(queries))
        var total = 0
        for qi in range(len(queries)):
            var one = List[Span[UInt8, origin_of(buffer)]]()
            one.append(qspans[qi])
            var expected_count = 0
            if len(qspans[qi]) > 0:
                var qtext = SeqText(one)
                var qrc = SeqText(one, reverse_complement=True)
                var mem = MemIndex(rtext, MemParams.choose(min_len, len(rtext), len(qtext)), 1)
                var want = find_smems(mem, rtext, qtext, qrc, 1)
                expected_count = len(want)
                if got.offsets[qi + 1] - got.offsets[qi] != len(want):
                    raise Error(
                        String(
                            t"L={min_len} query {qi}: {got.offsets[qi + 1] - got.offsets[qi]} SMEMs, expected"
                            t" {len(want)}"
                        )
                    )
                for s in range(len(want)):
                    var a = got.smems[got.offsets[qi] + s]
                    var b = want.smems[s]
                    if a.start != b.start or a.end != b.end or a.count != b.count:
                        raise Error(String(t"L={min_len} query {qi} SMEM {s}: {a} expected {b}"))
                    var ha = got.hits_of(got.offsets[qi] + s)
                    var hb = want.hits_of(s)
                    assert_equal(len(ha), len(hb))
                    for h in range(len(ha)):
                        if ha[h] != hb[h]:
                            raise Error(String(t"L={min_len} query {qi} SMEM {s} hit {h}: {ha[h]} expected {hb[h]}"))
                    var hc = capped.hits_of(capped.offsets[qi] + s)
                    assert_equal(len(hc), len(hb) if len(hb) <= 2 else 0)
            assert_equal(got.offsets[qi + 1] - got.offsets[qi], expected_count)
            total += expected_count
        assert_true(total > 0)


def test_fm_smems_small() raises:
    # Tiny: no k-mer cache. Two sequences, one holding a reverse-complement copy of the other's piece.
    var ref1 = random_dna(300, 11, 0.004)
    var ref2 = random_dna(200, 12)
    ref2.extend(Span(_reverse_complement(Span(ref1)[40:120])))
    var reference = ref1.copy()
    reference.extend(Span(ref2))
    var cuts: List[Int] = [0, len(ref1), len(reference)]
    _check(reference, cuts, _reads(reference, 60, 5), 1, False)


def test_fm_smems_repeats() raises:
    # Repeats (a duplicated block) and N runs, the k-mer cache on, several threads, save and load.
    var reference = random_dna(30_000, 21, 0.0005)
    for i in range(2000):
        reference[20_000 + i] = reference[3_000 + i]
    var extra = random_dna(5_000, 22)
    reference.extend(Span(extra))
    var cuts: List[Int] = [0, 12_000, 12_001, 30_000, len(reference)]
    _check(reference, cuts, _reads(reference, 300, 6), 3, True)


def test_fm_sample_spacing() raises:
    # Every row sampled, and one in 256 (long locate walks through separators).
    var reference = random_dna(20_000, 41, 0.002)
    var cuts: List[Int] = [0, 7_000, len(reference)]
    for bits in [1, 8]:
        _check(reference, cuts, _reads(reference, 120, 7 + bits), 2, False, bits)


def test_fm_wide_suffix_array() raises:
    # The 64-bit suffix array path (references past 2^30 bases) writes the same index.
    var reference = random_dna(40_000, 31, 0.001)
    var spans = List[Span[UInt8, origin_of(reference)]]()
    spans.append(Span(reference)[:25_000])
    spans.append(Span(reference)[25_000:])
    var names: List[String] = ["a", "b"]
    for bits in [2, 4]:
        FMIndex.build(spans, names, 3, bits).save("build/test_fm_narrow.sfm")
        FMIndex.build(spans, names, 3, bits, wide_suffix_array=True).save("build/test_fm_wide.sfm")
        with open("build/test_fm_narrow.sfm", "r") as f:
            var narrow = f.read_bytes()
            with open("build/test_fm_wide.sfm", "r") as g:
                var wide = g.read_bytes()
                assert_equal(len(narrow), len(wide))
                for i in range(len(narrow)):
                    if narrow[i] != wide[i]:
                        raise Error(String(t"sample_bits={bits}: indexes differ at byte {i}"))


def test_fm_smems_no_bases() raises:
    # A reference with no ACGT at all: nothing matches.
    var reference = List[UInt8](length=50, fill=UInt8(ord("N")))
    var spans = List[Span[UInt8, origin_of(reference)]]()
    spans.append(Span(reference))
    var names: List[String] = ["n"]
    var fm = FMIndex.build(spans, names, 1)
    var query = random_dna(100, 3)
    var qspans = List[Span[UInt8, origin_of(query)]]()
    qspans.append(Span(query))
    var got = fm.smems(qspans, 1)
    assert_equal(len(got.smems), 0)
    # Saved and loaded, it is still empty.
    fm.save("build/test_fm_empty.sfm")
    var loaded = FMIndex.load("build/test_fm_empty.sfm")
    assert_equal(len(loaded.smems(qspans, 1).smems), 0)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
