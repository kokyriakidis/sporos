"""The `sporos` command line.

    sporos mem [-mum | -mumreference | -maxmatch] [-l L] [-b | -r] [-F] [-n]
               [-t threads] [-v] <reference.fa> <query.fa>
    sporos minimizers -k K -w W [--forward] [-t threads] <seqs.fa>
    sporos syncmers -k K -s S [--open T] [--forward] [-t threads] <seqs.fa>

`sporos mem` takes MUMmer's `mummer` flags and prints its output format;
`-mumreference` is the default, as in MUMmer. Non-ACGT bases never match
(MUMmer's `-n`, accepted and implied).
"""

from std.sys import argv, stderr
from std.time import perf_counter_ns

from sporos.alphabet import CodeBuffer
from sporos.fasta import FastaRecords, read_fasta
from sporos.mems import Mem, MemIndex, MemParams, filter_mums, find_mems
from sporos.minimizers import minimizer_positions_from_codes
from sporos.output import OutputBuffer
from sporos.parallel import default_threads
from sporos.syncmers import syncmer_positions_from_codes
from sporos.text import SeqText
from sporos.window import PADDING

comptime USAGE = """usage:
  sporos mem [-mum | -mumreference | -maxmatch] [-l L] [-b | -r] [-F] [-n] [-t threads] [-v] <reference.fa> <query.fa>
  sporos minimizers -k K -w W [--forward] [-t threads] <seqs.fa>
  sporos syncmers -k K -s S [--open T] [--forward] [-t threads] <seqs.fa>"""


def _log(verbose: Bool, message: String):
    if verbose:
        print(message, file=stderr)


def _seconds(start: Int) -> Float64:
    return Float64(perf_counter_ns() - start) / 1e9


def run_mem(args: List[String]) raises:
    var mode = "mumreference"
    var min_len = 20
    var forward = True
    var reverse = False
    var four_column = False
    var threads = default_threads()
    var verbose = False
    var paths = List[String]()
    var i = 0
    while i < len(args):
        var a = args[i]
        if a == "-mum" or a == "-mumreference" or a == "-maxmatch":
            mode = String(a.removeprefix("-"))
        elif a == "-l":
            i += 1
            min_len = Int(args[i])
        elif a == "-b":
            reverse = True
        elif a == "-r":
            forward = False
            reverse = True
        elif a == "-F":
            four_column = True
        elif a == "-n":
            pass
        elif a == "-t" or a == "-threads":
            i += 1
            threads = Int(args[i])
        elif a == "-v":
            verbose = True
        elif a.startswith("-"):
            raise Error(t"unknown option {a}\n{USAGE}")
        else:
            paths.append(a)
        i += 1
    if len(paths) != 2:
        raise Error(USAGE)

    var start = perf_counter_ns()
    var reference = read_fasta(paths[0])
    var query = read_fasta(paths[1])
    _log(verbose, String(t"read inputs: {_seconds(start)} s"))

    start = perf_counter_ns()
    var ref_seqs = List[Span[UInt8, origin_of(reference.data)]]()
    for r in range(len(reference)):
        ref_seqs.append(reference.seq(r))
    var rtext = SeqText(ref_seqs, is_query=False)
    _log(verbose, String(t"  text: {_seconds(start)} s"))
    var params = MemParams.choose(min_len, len(rtext), query.total_length())
    var index = MemIndex(rtext, params, threads)
    _log(
        verbose,
        String(
            t"index: {_seconds(start)} s ({len(index.entries)} anchors, K={params.k}"
            t" k1={params.ref_step} k2={params.query_step})"
        ),
    )

    four_column = four_column or len(reference) > 1
    var name_width = 0
    for name in reference.names:
        name_width = max(name_width, name.byte_length())
    var out = OutputBuffer()
    var query_start = perf_counter_ns()
    var total = 0
    for q in range(len(query)):
        var qseqs = List[Span[UInt8, origin_of(query.data)]]()
        qseqs.append(query.seq(q))
        var qtext = SeqText(qseqs, is_query=True)
        if forward:
            total += _report(
                out, index, rtext, qtext, query.names[q], False, mode,
                four_column, reference, name_width, threads,
            )
        if reverse:
            var rc = SeqText(reverse_complement_of=qtext)
            total += _report(
                out, index, rtext, rc, query.names[q], True, mode,
                four_column, reference, name_width, threads,
            )
    out.flush()
    _log(verbose, String(t"matches: {_seconds(query_start)} s ({total} reported)"))


def _report(
    mut out: OutputBuffer,
    index: MemIndex,
    rtext: SeqText,
    qtext: SeqText,
    name: String,
    reverse: Bool,
    mode: String,
    four_column: Bool,
    reference: FastaRecords,
    name_width: Int,
    threads: Int,
) -> Int:
    var mems = find_mems(index, rtext, qtext, threads)
    var hits: List[Mem]
    if mode == "maxmatch":
        hits = mems^
    elif mode == "mum":
        hits = filter_mums(mems, qtext)
        # MUMmer prints -mum matches in reference order.
        def by_ref(a: Mem, b: Mem) -> Bool:
            return a.ref_pos < b.ref_pos

        sort(Span(hits), by_ref)
    else:
        hits = filter_mums(mems, qtext, unique_in_query=False)
    out.write("> ")
    out.write(name)
    if reverse:
        out.write(" Reverse")
    out.write_byte(ord("\n"))
    for m in hits:
        var ref_pos = Int(m.ref_pos)
        if four_column:
            var seq = rtext.sequence_of(ref_pos)
            ref_pos -= rtext.starts[seq]
            out.write("  ")
            out.write(reference.names[seq])
            for _ in range(name_width + 1 - reference.names[seq].byte_length()):
                out.write_byte(ord(" "))
            out.write_byte(ord(" "))
        out.write_int(ref_pos + 1, 8)
        out.write("  ")
        out.write_int(Int(m.query_pos) + 1, 8)
        out.write("  ")
        out.write_int(Int(m.length), 8)
        out.write_byte(ord("\n"))
    return len(hits)


def run_seeds(command: String, args: List[String]) raises:
    var k = 15
    var w = 10
    var s = 0
    var open_offset: Optional[Int] = None
    var canonical = True
    var threads = default_threads()
    var paths = List[String]()
    var i = 0
    while i < len(args):
        var a = args[i]
        if a == "-k":
            i += 1
            k = Int(args[i])
        elif a == "-w":
            i += 1
            w = Int(args[i])
        elif a == "-s":
            i += 1
            s = Int(args[i])
        elif a == "--open":
            i += 1
            open_offset = Int(args[i])
        elif a == "--forward":
            canonical = False
        elif a == "-t" or a == "-threads":
            i += 1
            threads = Int(args[i])
        elif a.startswith("-"):
            raise Error(t"unknown option {a}\n{USAGE}")
        else:
            paths.append(a)
        i += 1
    if len(paths) != 1:
        raise Error(USAGE)
    var records = read_fasta(paths[0])
    var out = OutputBuffer()
    for r in range(len(records)):
        var seq = records.seq(r)
        var codes = CodeBuffer.encode(seq, PADDING)
        var positions: List[UInt32]
        if command == "minimizers":
            positions = minimizer_positions_from_codes(
                codes.unsafe_ptr(), len(seq), k, w, canonical, threads
            )
        else:
            if s == 0:
                raise Error("syncmers need -s")
            positions = syncmer_positions_from_codes(
                codes.unsafe_ptr(), len(seq), k, s, canonical, open_offset, threads
            )
        for p in positions:
            out.write(records.names[r])
            out.write_byte(ord("\t"))
            out.write_int(Int(p))
            out.write_byte(ord("\n"))
    out.flush()


def main() raises:
    var raw = argv()
    if len(raw) < 2:
        print(USAGE)
        return
    var command = String(raw[1])
    var args = List[String]()
    for i in range(2, len(raw)):
        args.append(String(raw[i]))
    if command == "mem":
        run_mem(args)
    elif command == "minimizers" or command == "syncmers":
        run_seeds(command, args)
    else:
        print(USAGE)
