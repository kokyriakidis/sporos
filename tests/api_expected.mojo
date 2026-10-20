"""Writes the Mojo library's answers on fixed inputs, which the C and Python tests hold their bindings to.

`tests/test_c_api.cpp` and `tests/test_python.py` read the file:

    api_expected <out.txt>

One record per line, fields separated by spaces:

    seq <name> <bases>
    minimizers <seq> <k> <w> <canonical> <threads> <count> <positions...>
    syncmers <seq> <k> <s> <open offsets, comma-separated; -1 for closed> <canonical> <threads> <count> <positions...>
    matches <query> <min_len> <mode> <reverse> <threads> <count> (<ref_seq> <ref_pos> <query_pos> <length>)...

`matches` index every `ref*` sequence, in their order; mode 0 is maxmatch, 1 mum, 2 mumreference.
"""

from std.sys import argv

from sporos.alphabet import encode_base
from sporos.mems import Mem, MemIndex, MemParams, filter_mums, find_mems
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
) raises:
    var spans = List[Span[UInt8, origin_of(refs)]]()
    for r in refs:
        spans.append(Span(r))
    var rtext = SeqText(spans)
    var index = MemIndex(rtext, MemParams.choose(min_len, len(rtext), len(rtext)), threads)
    var qspans = List[Span[UInt8, origin_of(query)]]()
    qspans.append(Span(query))
    var qtext = SeqText(qspans, reverse_complement=reverse)
    var mems = find_mems(index, rtext, qtext, threads)
    var kept: List[Mem]
    match mode:
        case 1:
            kept = filter_mums(mems, qtext)
        case 2:
            kept = filter_mums(mems, qtext, unique_in_query=False)
        case _:
            kept = mems^
    out.write("matches ", name, " ", min_len, " ", mode, " ", Int(reverse), " ", threads, " ", len(kept))
    for m in kept:
        var seq = rtext.sequence_of(Int(m.ref_pos))
        out.write(" ", seq, " ", Int(m.ref_pos) - rtext.starts[seq], " ", m.query_pos, " ", m.length)
    out.write("\n")


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

    with open(args[1], "w") as f:
        f.write(out)
