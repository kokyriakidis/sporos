"""Writes the Mojo library's answers on fixed inputs, which the C and Python tests hold their bindings to.

`tests/test_c_api.cpp` and `tests/test_python.py` read the file:

    api_expected <out.txt>

One record per line, fields separated by spaces:

    seq <name> <bases>
    minimizers <seq> <k> <w> <canonical> <threads> <count> <positions...>
    syncmers <seq> <k> <s> <open offsets, comma-separated; -1 for closed> <canonical> <threads> <count> <positions...>
    matches <query> <min_len> <mode> <reverse> <max_occ> <min_occ> <threads> <count> (<ref_seq> <ref_pos> <query_pos> <length>)...
    smems <read> <reads> <min_len> <max_hits> <min_occ> <max_occ> <threads> (<count> (<start> <end> <occurrences> <hits> (<ref_seq> <ref_pos> <reverse>)...)...)...
    seeds <read> <reads> <k> <w> <scheme> <open offset; -1 for none> <max_occ> <all_kmers> <threads> (<count> (<query_pos> <ref_seq> <ref_pos> <occurrences> <reverse>)...)...
    cutoff <k> <w> <scheme> <open offset> <keep_fraction> <cutoff>

`matches`, `smems`, `seeds` and `cutoff` index `ref0` and `ref1`, in that order; mode 0 is maxmatch,
1 mum, 2 mumreference; scheme 0 is minimizers, 1 syncmers, 2 every k-mer. `smems` and `seeds` take the
reads `<read>0` to `<read><reads - 1>`, each one's answers in turn.
"""

from std.sys import argv

from sporos.alphabet import encode_base
from sporos.fmindex import FMIndex
from sporos.index import SeedIndex
from sporos.mems import MemIndex, MemParams, filter_by_occurrences, find_mems
from sporos.minimizers import minimizer_positions
from sporos.syncmers import syncmer_positions
from sporos.text import SeqText
from reference import random_dna


def reverse_complement(seq: List[UInt8]) -> List[UInt8]:
    """`seq`'s reverse complement, upper case, anything but a base as N."""
    var out = List[UInt8](capacity=len(seq))
    for i in reversed(range(len(seq))):
        var code = Int(encode_base(seq[i]))
        out.append(UInt8(ord("TGCAN"[byte=code])))
    return out^


def slice_of(seq: List[UInt8], start: Int, end: Int) -> List[UInt8]:
    var out = List[UInt8](capacity=end - start)
    for i in range(start, end):
        out.append(seq[i])
    return out^


def write_seq(mut out: String, name: String, seq: List[UInt8]):
    out.write("seq ", name, " ", StringSlice(unsafe_from_utf8=Span(seq)), "\n")


def write_positions(mut out: String, positions: List[UInt32]):
    out.write(len(positions))
    for p in positions:
        out.write(" ", p)
    out.write("\n")


def write_minimizers(
    mut out: String, name: String, seq: List[UInt8], k: Int, w: Int, canonical: Bool, threads: Int
) raises:
    out.write("minimizers ", name, " ", k, " ", w, " ", Int(canonical), " ", threads, " ")
    write_positions(out, minimizer_positions(Span(seq), k, w, canonical, threads))


def write_syncmers(
    mut out: String,
    name: String,
    seq: List[UInt8],
    k: Int,
    s: Int,
    open_offsets: List[Int],
    canonical: Bool,
    threads: Int,
) raises:
    out.write("syncmers ", name, " ", k, " ", s, " ")
    if len(open_offsets) == 0:
        out.write("-1")
    for i in range(len(open_offsets)):
        out.write("," if i > 0 else "", open_offsets[i])
    out.write(" ", Int(canonical), " ", threads, " ")
    write_positions(out, syncmer_positions(Span(seq), k, s, canonical, None, threads, open_offsets=open_offsets))


def write_matches(
    mut out: String,
    refs: List[List[UInt8]],
    name: String,
    query: List[UInt8],
    min_len: Int,
    mode: Int,
    reverse: Bool,
    threads: Int,
    max_occurrences: Int = 0,
    min_occurrences: Int = 1,
) raises:
    var spans = List[Span[UInt8, origin_of(refs)]]()
    for r in refs:
        spans.append(Span(r))
    var rtext = SeqText(spans)
    var index = MemIndex(rtext, MemParams.choose(min_len, len(rtext), len(rtext)), threads)
    var qspans = List[Span[UInt8, origin_of(query)]]()
    qspans.append(Span(query))
    var qtext = SeqText(qspans, reverse_complement=reverse)
    var cap = max_occurrences if max_occurrences > 0 else (0 if mode == 0 else 1)
    var kept = filter_by_occurrences(
        find_mems(index, rtext, qtext, threads), qtext, cap, cap if mode == 1 else 0, min_occurrences, threads
    )
    out.write("matches ", name, " ", min_len, " ", mode, " ", Int(reverse), " ", max_occurrences, " ")
    out.write(min_occurrences, " ", threads, " ", len(kept))
    for m in kept:
        var seq = rtext.sequence_of(Int(m.ref_pos))
        out.write(" ", seq, " ", Int(m.ref_pos) - rtext.starts[seq], " ", m.query_pos, " ", m.length)
    out.write("\n")


def write_smems(
    mut out: String,
    fm: FMIndex,
    reads: List[List[UInt8]],
    min_len: Int,
    max_hits: Int,
    min_occurrences: Int,
    max_occurrences: Int,
    threads: Int,
):
    var spans = List[Span[UInt8, origin_of(reads)]]()
    for r in reads:
        spans.append(Span(r))
    var found = fm.smems(spans, min_len, max_hits, min_occurrences, max_occurrences, threads)
    out.write("smems read ", len(reads), " ", min_len, " ", max_hits, " ", min_occurrences, " ")
    out.write(max_occurrences, " ", threads)
    for q in range(len(reads)):
        out.write(" ", found.offsets[q + 1] - found.offsets[q])
        for i in range(found.offsets[q], found.offsets[q + 1]):
            var smem = found.smems[i]
            out.write(" ", smem.start, " ", smem.end, " ", smem.count, " ", len(found.hits_of(i)))
            for h in found.hits_of(i):
                out.write(" ", h.sequence, " ", h.position, " ", Int(h.reverse))
    out.write("\n")


def fm_index(refs: List[List[UInt8]]) raises -> FMIndex:
    var spans = List[Span[UInt8, origin_of(refs)]]()
    for r in refs:
        spans.append(Span(r))
    var names: List[String] = ["ref0", "ref1"]
    return FMIndex.build(spans, names)


def seed_index(refs: List[List[UInt8]], k: Int, w: Int, scheme: Int, open_offset: Int) raises -> SeedIndex:
    var spans = List[Span[UInt8, origin_of(refs)]]()
    for r in refs:
        spans.append(Span(r))
    return SeedIndex(spans, k, w, scheme, 1, Optional[Int](open_offset) if open_offset >= 0 else None)


def write_seeds(
    mut out: String,
    refs: List[List[UInt8]],
    reads: List[List[UInt8]],
    k: Int,
    w: Int,
    scheme: Int,
    open_offset: Int,
    max_occurrences: Int,
    all_kmers: Bool,
    threads: Int,
) raises:
    var index = seed_index(refs, k, w, scheme, open_offset)
    var spans = List[Span[UInt8, origin_of(reads)]]()
    for r in reads:
        spans.append(Span(r))
    var found = index.lookup_batch(spans, all_kmers, max_occurrences, threads)
    out.write("seeds read ", len(reads), " ", k, " ", w, " ", scheme, " ", open_offset, " ", max_occurrences, " ")
    out.write(Int(all_kmers), " ", threads)
    for r in range(len(reads)):
        out.write(" ", len(found.of(r)))
        for h in found.of(r):
            var seq = index.text.sequence_of(Int(h.ref_pos))
            out.write(" ", h.query_pos, " ", seq, " ", Int(h.ref_pos) - index.text.starts[seq], " ")
            out.write(h.occurrences, " ", Int(h.reverse))
    out.write("\n")


def write_cutoff(
    mut out: String, refs: List[List[UInt8]], k: Int, w: Int, scheme: Int, open_offset: Int, keep: Float64
) raises:
    var index = seed_index(refs, k, w, scheme, open_offset)
    out.write("cutoff ", k, " ", w, " ", scheme, " ", open_offset, " ", keep, " ", index.occurrence_cutoff(keep), "\n")


def main() raises:
    var args = argv()
    if len(args) != 2:
        raise Error("usage: api_expected <out.txt>")
    var out = String()

    var mixed = random_dna(5_000, 11, n_rate=0.002)
    var short = random_dna(37, 3)
    var long = random_dna(200_000, 5, n_rate=0.0005)
    write_seq(out, "mixed", mixed)
    write_seq(out, "short", short)
    write_seq(out, "long", long)
    write_minimizers(out, "mixed", mixed, 21, 11, True, 1)
    write_minimizers(out, "mixed", mixed, 15, 5, False, 1)
    write_minimizers(out, "mixed", mixed, 64, 3, True, 1)
    write_minimizers(out, "short", short, 5, 40, True, 1)
    write_minimizers(out, "long", long, 31, 19, True, 4)
    write_minimizers(out, "long", long, 19, 10, False, 3)
    write_syncmers(out, "mixed", mixed, 21, 11, [], True, 1)
    write_syncmers(out, "mixed", mixed, 21, 11, [5], True, 1)
    write_syncmers(out, "mixed", mixed, 21, 11, [1, 5, 9], True, 1)
    write_syncmers(out, "mixed", mixed, 15, 7, [], False, 1)
    write_syncmers(out, "short", short, 9, 4, [2], False, 1)
    write_syncmers(out, "long", long, 31, 15, [], True, 4)
    write_syncmers(out, "long", long, 25, 12, [6], True, 2)
    write_syncmers(out, "long", long, 25, 12, [0, 6, 13], False, 2)

    # Two references and a query of their pieces: forward, reverse-complemented, repeated, with noise.
    var ref0 = random_dna(20_000, 21)
    var ref1 = random_dna(15_000, 22, n_rate=0.001)
    var query = slice_of(ref0, 1_000, 3_000)
    query.extend(Span(random_dna(500, 23)))
    query.extend(Span(reverse_complement(slice_of(ref1, 2_000, 4_000))))
    query.extend(Span(slice_of(ref0, 1_500, 2_200)))
    query.extend(Span(slice_of(ref1, 9_000, 9_800)))
    query.extend(Span(slice_of(ref0, 15_000, 15_600)))
    # A piece repeated across the references, so occurrence limits have something to count.
    for i in range(600):
        ref1[12_000 + i] = ref0[5_000 + i]
    query.extend(Span(slice_of(ref0, 5_000, 5_600)))
    for i in range(len(query)):
        if i % 97 == 0:
            query[i] = UInt8(ord("A")) if query[i] != UInt8(ord("A")) else UInt8(ord("C"))
    write_seq(out, "ref0", ref0)
    write_seq(out, "ref1", ref1)
    write_seq(out, "query", query)
    var refs = List[List[UInt8]]()
    refs.append(ref0^)
    refs.append(ref1^)
    for mode in range(3):
        for reverse in range(2):
            write_matches(out, refs, "query", query, 20, mode, reverse == 1, 1)
    write_matches(out, refs, "query", query, 50, 0, False, 4)
    write_matches(out, refs, "query", query, 50, 1, True, 4)
    write_matches(out, refs, "query", query, 20, 0, False, 1, max_occurrences=1)
    write_matches(out, refs, "query", query, 20, 0, False, 2, min_occurrences=2)
    write_matches(out, refs, "query", query, 20, 1, False, 1, max_occurrences=2)
    write_matches(out, refs, "query", query, 20, 2, False, 1, max_occurrences=3, min_occurrences=2)

    # Reads cut from the references, some reverse-complemented, with a substitution every 37 bases.
    var reads = List[List[UInt8]]()
    for r in range(40):
        var source = Span(refs[r % 2])
        var at = (r * 7_919) % (len(source) - 200)
        if r == 39:
            at = 12_000  # the repeat
        var read = slice_of(refs[r % 2], at, at + 150)
        for i in range(r % 5, len(read), 37):
            read[i] = UInt8(ord("A")) if read[i] != UInt8(ord("A")) else UInt8(ord("G"))
        if r % 3 == 0:
            read = reverse_complement(read)
        write_seq(out, String(t"read{r}"), read)
        reads.append(read^)
    var fm = fm_index(refs)
    write_smems(out, fm, reads, 19, 20, 1, 0, 1)
    write_smems(out, fm, reads, 19, 1, 1, 0, 2)
    write_smems(out, fm, reads, 15, 20, 2, 0, 2)
    write_smems(out, fm, reads, 19, 20, 1, 1, 1)
    write_seeds(out, refs, reads, 15, 10, 0, -1, 0, False, 1)
    write_seeds(out, refs, reads, 15, 10, 0, -1, 1, False, 2)
    write_seeds(out, refs, reads, 15, 10, 0, -1, 0, True, 1)
    write_seeds(out, refs, reads, 17, 8, 1, -1, 0, False, 1)
    write_seeds(out, refs, reads, 17, 9, 1, 2, 0, False, 1)
    write_seeds(out, refs, reads, 13, 1, 2, -1, 0, False, 3)
    write_cutoff(out, refs, 15, 10, 0, -1, 0.9)
    write_cutoff(out, refs, 13, 1, 2, -1, 0.999)

    with open(args[1], "w") as f:
        f.write(out)
