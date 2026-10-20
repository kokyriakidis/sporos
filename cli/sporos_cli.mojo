"""The `sporos` command line.

    sporos mem [-mum | -mumreference | -maxmatch] [-l L] [-b | -r] [-F] [-n]
               [-t threads] [-v] <reference.fa> <query.fa>
    sporos minimizers -k K -w W [--forward] [-t threads] <seqs.fa>
    sporos syncmers -k K -s S [--open T] [--forward] [-t threads] <seqs.fa>

`sporos mem` takes MUMmer's `mummer` flags and prints its output format;
`-mumreference` is the default, as in MUMmer. Non-ACGT bases never match
(MUMmer's `-n`, accepted and implied).
"""

from std.math import ceildiv
from std.sys import argv, stderr
from std.time import perf_counter_ns

from sporos.alphabet import CodeBuffer
from sporos.fasta import FastaRecords, read_fasta
from sporos.mems import Mem, MemIndex, MemParams, filter_mums, find_mems
from sporos.minimizers import minimizer_positions_from_codes
from sporos.output import OutputBuffer, put_padded
from sporos.parallel import ParallelTask, default_threads, parallel_for
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
    var rtext = SeqText(ref_seqs)
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
        var qtext = SeqText(qseqs)
        if forward:
            total += _report(
                out, index, rtext, qtext, query.names[q], False, mode,
                four_column, reference, name_width, threads,
            )
        if reverse:
            var rc = SeqText(qseqs, reverse_complement=True)
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
    # Format slices of the matches in parallel, then write them in order.
    var slices = max(1, min(threads * 4, len(hits) // 50_000))
    var texts = List[List[UInt8]](length=slices, fill=List[UInt8]())
    var task = _FormatTask(
        Pointer(to=hits).unsafe_origin_cast[MutUntrackedOrigin](),
        Pointer(to=rtext).unsafe_origin_cast[ImmUntrackedOrigin](),
        Pointer(to=reference).unsafe_origin_cast[ImmUntrackedOrigin](),
        texts.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        slices,
        four_column,
        name_width,
    )
    parallel_for(task, slices, threads)
    for text in texts:
        out.write_bytes(Span(text))
    return len(hits)


struct _FormatTask(ParallelTask):
    """Formats one slice of matches as MUMmer lines."""

    var hits: ImmPointer[List[Mem], ImmUntrackedOrigin]
    var rtext: ImmPointer[SeqText, ImmUntrackedOrigin]
    var reference: ImmPointer[FastaRecords, ImmUntrackedOrigin]
    var texts: MutPointer[List[UInt8], MutUntrackedOrigin]
    var slices: Int
    var four_column: Bool
    var name_width: Int

    def __init__(
        out self,
        hits: ImmPointer[List[Mem], ImmUntrackedOrigin],
        rtext: ImmPointer[SeqText, ImmUntrackedOrigin],
        reference: ImmPointer[FastaRecords, ImmUntrackedOrigin],
        texts: MutPointer[List[UInt8], MutUntrackedOrigin],
        slices: Int,
        four_column: Bool,
        name_width: Int,
    ):
        self.hits = hits
        self.rtext = rtext
        self.reference = reference
        self.texts = texts
        self.slices = slices
        self.four_column = four_column
        self.name_width = name_width

    def run(self, index: Int):
        ref hits = self.hits[]
        var per = ceildiv(len(hits), self.slices)
        var lo = index * per
        var hi = min(len(hits), lo + per)
        ref text = self.texts[unsafe_offset=index]
        var line_max = 3 * 22 + 8 + (self.name_width + 4 if self.four_column else 0)
        text.resize(unsafe_uninit_length=max(0, hi - lo) * line_max)
        var p = text.unsafe_ptr()
        var at = 0
        for i in range(lo, hi):
            var m = hits[i]
            var ref_pos = Int(m.ref_pos)
            if self.four_column:
                var seq = self.rtext[].sequence_of(ref_pos)
                ref_pos -= self.rtext[].starts[seq]
                ref name = self.reference[].names[seq]
                p[unsafe_offset=at] = UInt8(ord(" "))
                p[unsafe_offset=at + 1] = UInt8(ord(" "))
                at += 2
                var bytes = name.as_bytes()
                for b in range(len(bytes)):
                    p[unsafe_offset=at + b] = bytes[b]
                at += len(bytes)
                for _ in range(self.name_width + 2 - len(bytes)):
                    p[unsafe_offset=at] = UInt8(ord(" "))
                    at += 1
            at += put_padded(p.unsafe_offset(at), ref_pos + 1, 8)
            p[unsafe_offset=at] = UInt8(ord(" "))
            p[unsafe_offset=at + 1] = UInt8(ord(" "))
            at += 2
            at += put_padded(p.unsafe_offset(at), Int(m.query_pos) + 1, 8)
            p[unsafe_offset=at] = UInt8(ord(" "))
            p[unsafe_offset=at + 1] = UInt8(ord(" "))
            at += 2
            at += put_padded(p.unsafe_offset(at), Int(m.length), 8)
            p[unsafe_offset=at] = UInt8(ord("\n"))
            at += 1
        text.resize(unsafe_uninit_length=at)


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
