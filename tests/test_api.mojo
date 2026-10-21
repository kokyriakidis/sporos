"""The public API (`from sporos import ...`): `Reference` against the
building blocks it wraps, and threaded FM-index SMEMs against one thread."""

from std.testing import assert_equal, assert_raises, assert_true, TestSuite

from sporos import (
    FMIndex,
    Match,
    MemIndex,
    MemParams,
    Reference,
    SeqText,
    filter_by_occurrences,
    find_mems,
    minimizer_positions,
)
from reference import random_dna


def _reverse_complement(seq: Span[UInt8, _]) -> List[UInt8]:
    var out = List[UInt8](capacity=len(seq))
    for i in reversed(range(len(seq))):
        var c = seq[i] | 0x20
        if c == UInt8(ord("a")):
            out.append(UInt8(ord("T")))
        elif c == UInt8(ord("c")):
            out.append(UInt8(ord("G")))
        elif c == UInt8(ord("g")):
            out.append(UInt8(ord("C")))
        elif c == UInt8(ord("t")):
            out.append(UInt8(ord("A")))
        else:
            out.append(UInt8(ord("N")))
    return out^


def test_reference_mems_match_building_blocks() raises:
    # Two reference sequences; a query of their pieces, one reverse-complemented.
    var ref1 = random_dna(20_000, 1)
    var ref2 = random_dna(15_000, 2)
    for i in range(500):
        ref2[9_000 + i] = ref1[3_000 + i]  # a repeat across sequences
    var query = List[UInt8](Span(ref1)[1_000:6_000])
    query.extend(Span(_reverse_complement(Span(ref2)[2_000:4_000])))
    query.extend(Span(random_dna(300, 3)))
    var spans = List[Span[UInt8, origin_of(ref1)]]()
    spans.append(Span(ref1))
    var both = ref1.copy()
    both.extend(Span(ref2))
    var seqs = List[Span[UInt8, origin_of(both)]]()
    seqs.append(Span(both)[: len(ref1)])
    seqs.append(Span(both)[len(ref1) :])
    var names: List[String] = ["chrA", "chrB"]
    var reference = Reference(seqs, names)
    assert_equal(len(reference), 2)
    assert_equal(reference.name(1), "chrB")
    var rtext = SeqText(seqs)
    for mode in ["maxmatch", "mum", "mumreference"]:
        for cap in [0, 3]:
            var got = reference.mems(Span(query), 30, mode, forward=True, reverse=True, max_occurrences=cap, threads=2)
            # The same through the building blocks, mapped to (sequence, position).
            var want = List[Match]()
            var index = MemIndex(rtext, MemParams.choose(30, len(rtext), len(query)), 1)
            var limit = cap if cap > 0 else (0 if mode == "maxmatch" else 1)
            for strand in range(2):
                var q = List[Span[UInt8, origin_of(query)]]()
                q.append(Span(query))
                var qtext = SeqText(q, reverse_complement=strand == 1)
                var kept = filter_by_occurrences(
                    find_mems(index, rtext, qtext), qtext, limit, limit if mode == "mum" else 0, 1
                )
                for m in kept:
                    var seq = rtext.sequence_of(Int(m.ref_pos))
                    want.append(
                        Match(seq, Int(m.ref_pos) - rtext.starts[seq], Int(m.query_pos), Int(m.length), strand == 1)
                    )
            assert_equal(len(got), len(want))
            for i in range(len(got)):
                assert_true(got[i] == want[i])
            assert_true(len(got) > 0)
    # The planted pieces: query[0:5000] = chrA[1000:6000]; query[5000:7000] =
    # reverse complement of chrB[2000:4000].
    var found_forward = False
    var found_reverse = False
    for m in reference.mems(Span(query), 100, "maxmatch", reverse=True):
        # (A match may run a few bases past a planted piece by chance.)
        if not m.reverse and m.sequence == 0 and m.ref_pos <= 1_000 and m.ref_pos + m.length >= 6_000:
            found_forward = True
        if m.reverse and m.sequence == 1 and m.ref_pos <= 2_000 and m.ref_pos + m.length >= 4_000:
            found_reverse = True
    assert_true(found_forward)
    assert_true(found_reverse)
    with assert_raises():
        _ = reference.mems(Span(query), 30, "mums")
    with assert_raises():
        _ = reference.mems(Span(query), 30, "mum", min_occurrences=2)


def test_reference_smems_and_threads() raises:
    var genome = random_dna(60_000, 5, 0.0005)
    var seqs = List[Span[UInt8, origin_of(genome)]]()
    seqs.append(Span(genome)[:30_000])
    seqs.append(Span(genome)[30_000:])
    var reference = Reference(seqs)
    assert_equal(reference.name(0), "seq0")
    # Many reads: enough for several slices when threaded.
    var reads = List[List[UInt8]]()
    var state = UInt64(99)
    for r in range(6_000):
        state = state * 6364136223846793005 + 1442695040888963407
        var at = Int((state >> 20) % 59_800)
        var read = List[UInt8](Span(genome)[at : at + 150])
        if r % 3 == 0:
            read[75] = UInt8(ord("A")) if read[75] != UInt8(ord("A")) else UInt8(ord("C"))
        if r % 2 == 1:
            read = _reverse_complement(Span(read))
        reads.append(read^)
    var buffer = List[UInt8]()
    var bounds: List[Int] = [0]
    for r in reads:
        buffer.extend(Span(r))
        bounds.append(len(buffer))
    var spans = List[Span[UInt8, origin_of(buffer)]]()
    for i in range(len(reads)):
        spans.append(Span(buffer)[bounds[i] : bounds[i + 1]])
    var one = reference.smems(spans, 19, threads=1)
    var many = reference.smems(spans, 19, threads=4)
    assert_equal(len(one), len(reads))
    assert_equal(len(many), len(reads))
    assert_equal(len(one.smems), len(many.smems))
    assert_equal(len(one.hits), len(many.hits))
    for q in range(len(one)):
        assert_equal(len(one.of(q)), len(many.of(q)))
    for i in range(len(one.smems)):
        assert_true(one.smems[i] == many.smems[i])
        assert_equal(len(one.hits_of(i)), len(many.hits_of(i)))
        for h in range(len(one.hits_of(i))):
            assert_true(one.hits_of(i)[h] == many.hits_of(i)[h])
    # A read cut from sequence 1 forward hits it at its own position.
    var read = List[UInt8](Span(genome)[40_000:40_150])
    var single = List[Span[UInt8, origin_of(read)]]()
    single.append(Span(read))
    var hit = reference.smems(single, 19)
    assert_equal(len(hit.smems), 1)
    var h = hit.hits_of(0)[0]
    assert_equal(Int(h.sequence), 1)
    assert_equal(Int(h.position), 10_000)
    assert_true(not h.reverse)
    # The FM-index saves and loads with the same answers.
    reference.save_fm_index("build/test_api.sfm")
    var loaded = FMIndex.load("build/test_api.sfm")
    var again = loaded.smems(spans, 19, threads=2)
    assert_equal(len(again.smems), len(one.smems))


def test_seed_functions_are_exported() raises:
    var seq = random_dna(500, 9)
    assert_true(len(minimizer_positions(Span(seq), 15, 10)) > 0)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
