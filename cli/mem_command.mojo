"""`sporos mem`: MEMs, MUMs and frequency-capped MEMs, in MUMmer's format."""

from std.ffi import external_call
from std.math import ceildiv
from std.memory import unsafe_memcpy
from std.os.path import exists
from std.sys import argv, stderr
from std.sys.info import CompilationTarget
from std.time import perf_counter_ns

from sporos.batch import SeedBatch, minimizers_batch, syncmers_batch
from sporos.fasta import FastaRecords, read_fasta
from sporos.fmindex import FMIndex, encode_bases
from sporos.index import ALL_KMERS, MINIMIZER_SEEDS, SYNCMER_SEEDS, SeedIndex
from sporos.mems import Mem, MemIndex, MemParams, Smems, filter_by_occurrences, find_mems, find_smems
from sporos.minimizers import minimizer_positions
from sporos.output import OutputBuffer, _PAIRS, put_padded
from sporos.parallel import ParallelTask, default_threads, parallel_for
from sporos.syncmers import syncmer_positions
from sporos.text import SeqText
from sporos.window import LANES

from cli.common import (
    USAGE,
    _FORMAT_SLICE,
    _SEED_SLICE,
    _SEED_SLICE_BASES,
    _count_option,
    _int_list,
    _log,
    _padded_field,
    _peak_memory,
    _put_field,
    _seconds,
    _threads,
    _value,
)


def run_mem(args: List[String]) raises:
    var kind = "mumreference"
    var max_occ = 0
    var min_occ = 1
    var min_len = 20
    var forward = True
    var reverse = False
    var four_column = False
    var threads = default_threads()
    var verbose = False
    var anchor_k = 0
    var ref_step = 0
    var query_step = 0
    var paths = List[String]()
    var i = 0
    while i < len(args):
        var a = args[i]
        if a == "-K":
            i += 1
            anchor_k = Int(_value(args, i))
        elif a == "-k1":
            i += 1
            ref_step = Int(_value(args, i))
        elif a == "-k2":
            i += 1
            query_step = Int(_value(args, i))
        elif a == "-mum" or a == "-mumreference" or a == "-maxmatch":
            kind = String(a.removeprefix("-"))
        elif a == "--max-occ":
            i += 1
            max_occ = _count_option(args, i, "--max-occ")
        elif a == "--min-occ":
            i += 1
            min_occ = _count_option(args, i, "--min-occ")
        elif a == "-l":
            i += 1
            min_len = Int(_value(args, i))
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
            threads = _threads(args, i)
        elif a == "-v":
            verbose = True
        elif a.startswith("-"):
            raise Error(t"unknown option {a}\n{USAGE}")
        else:
            paths.append(a)
        i += 1
    if len(paths) != 2:
        raise Error(USAGE)
    var cap = max_occ if max_occ > 0 else (0 if kind == "maxmatch" else 1)
    if cap > 0 and min_occ > cap:
        raise Error(t"--min-occ {min_occ} is above the most occurrences allowed ({cap}); raise it with --max-occ")
    var mode = _MemFilter(kind, max_occ, min_occ)

    var start = perf_counter_ns()
    var reference = read_fasta(paths[0], threads)
    var query = read_fasta(paths[1], threads)
    _log(verbose, String(t"read inputs: {_seconds(start)} s"))

    start = perf_counter_ns()
    var ref_seqs = List[Span[UInt8, origin_of(reference.data)]]()
    for r in range(len(reference)):
        ref_seqs.append(reference.seq(r))
    var rtext = SeqText(ref_seqs)
    # Matching reads the packed text; output needs only the names.
    reference.data = List[UInt8]()
    _log(verbose, String(t"  text: {_seconds(start)} s"))
    for q in range(len(query)):
        if len(query.seq(q)) >= Int(UInt32.MAX):
            raise Error("query sequences must be shorter than 2^32 bases")
    # One query record (a genome): pack it now and drop its raw bytes before
    # the index is built.
    var packed = List[SeqText]()
    if len(query) == 1:
        var qseqs = List[Span[UInt8, origin_of(query.data)]]()
        qseqs.append(query.seq(0))
        if forward:
            packed.append(SeqText(qseqs))
        if reverse:
            packed.append(SeqText(qseqs, reverse_complement=True))
        query.data = List[UInt8]()
    var params = MemParams.choose(min_len, len(rtext), query.total_length())
    if anchor_k > 0 or ref_step > 0 or query_step > 0:
        # Tuning overrides: anchor length and coprime sampling steps.
        params = MemParams(
            min_len,
            anchor_k if anchor_k > 0 else params.k,
            ref_step if ref_step > 0 else params.ref_step,
            query_step if query_step > 0 else params.query_step,
        )
    var index = MemIndex(rtext, params, threads)
    _log(
        verbose,
        String(
            t"index: {_seconds(start)} s ({len(index.table)} anchors, K={params.k}"
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
    if len(query) > 1:
        # Several query records (reads): slices of records (up to 4096 reads or
        # 256 Kbp) in parallel rounds, each round written before the next.
        var bounds: List[Int] = [0]
        var bases = 0
        for q in range(len(query)):
            bases += query.lengths[q]
            if bases >= _SEED_SLICE_BASES or q + 1 - bounds[len(bounds) - 1] >= _SEED_SLICE:
                bounds.append(q + 1)
                bases = 0
        if bounds[len(bounds) - 1] != len(query):
            bounds.append(len(query))
        var slices = len(bounds) - 1
        var round = threads * 4
        var texts = List[List[UInt8]](length=round, fill=List[UInt8]())
        var totals = List[Int](length=round, fill=0)
        for first in range(0, slices, round):
            var count = min(round, slices - first)
            var task = _RecordsTask(
                Pointer(to=index).as_imm(),
                Pointer(to=rtext).as_imm(),
                Pointer(to=query).as_imm(),
                Pointer(to=reference).as_imm(),
                Pointer(to=bounds).as_imm(),
                texts.unsafe_ptr(),
                totals.unsafe_ptr(),
                first,
                mode,
                forward,
                reverse,
                four_column,
                name_width,
            )
            parallel_for(task, count, threads)
            for s in range(count):
                out.write_bytes(Span(texts[s]))
                total += totals[s]
        out.flush()
        _log(verbose, String(t"matches: {_seconds(query_start)} s ({total} reported)"))
        return
    for q in range(len(query)):
        # The record's strands: packed above for a single record.
        var texts = packed^
        packed = List[SeqText]()
        if len(texts) == 0:
            var qseqs = List[Span[UInt8, origin_of(query.data)]]()
            qseqs.append(query.seq(q))
            if forward:
                texts.append(SeqText(qseqs))
            if reverse:
                texts.append(SeqText(qseqs, reverse_complement=True))
        if forward:
            total += _report(
                out,
                index,
                rtext,
                texts[0],
                query.names[q],
                False,
                mode,
                four_column,
                reference,
                name_width,
                threads,
            )
        if reverse:
            total += _report(
                out,
                index,
                rtext,
                texts[len(texts) - 1],
                query.names[q],
                True,
                mode,
                four_column,
                reference,
                name_width,
                threads,
            )
    out.flush()
    _log(verbose, String(t"matches: {_seconds(query_start)} s ({total} reported)"))


struct _RecordsTask[
    index_origin: ImmOrigin,
    text_origin: ImmOrigin,
    query_origin: ImmOrigin,
    ref_origin: ImmOrigin,
    bounds_origin: ImmOrigin,
    out_origin: MutOrigin,
    total_origin: MutOrigin,
](ParallelTask):
    """Reports the matches of one slice of query records, single-threaded."""

    var index: ImmPointer[MemIndex, Self.index_origin]
    var rtext: ImmPointer[SeqText, Self.text_origin]
    var query: ImmPointer[FastaRecords, Self.query_origin]
    var reference: ImmPointer[FastaRecords, Self.ref_origin]
    var bounds: ImmPointer[List[Int], Self.bounds_origin]
    var texts: MutPointer[List[UInt8], Self.out_origin]
    var totals: MutPointer[Int, Self.total_origin]
    var first: Int
    var mode: _MemFilter
    var forward: Bool
    var reverse: Bool
    var four_column: Bool
    var name_width: Int

    def __init__(
        out self,
        index: ImmPointer[MemIndex, Self.index_origin],
        rtext: ImmPointer[SeqText, Self.text_origin],
        query: ImmPointer[FastaRecords, Self.query_origin],
        reference: ImmPointer[FastaRecords, Self.ref_origin],
        bounds: ImmPointer[List[Int], Self.bounds_origin],
        texts: MutPointer[List[UInt8], Self.out_origin],
        totals: MutPointer[Int, Self.total_origin],
        first: Int,
        mode: _MemFilter,
        forward: Bool,
        reverse: Bool,
        four_column: Bool,
        name_width: Int,
    ):
        self.index = index
        self.rtext = rtext
        self.query = query
        self.reference = reference
        self.bounds = bounds
        self.texts = texts
        self.totals = totals
        self.first = first
        self.mode = mode.copy()
        self.forward = forward
        self.reverse = reverse
        self.four_column = four_column
        self.name_width = name_width

    def run(self, index: Int):
        ref query = self.query[]
        var lo = self.bounds[][self.first + index]
        var hi = self.bounds[][self.first + index + 1]
        var out = OutputBuffer(fd=-1)
        var total = 0
        for q in range(lo, hi):
            var qseqs = List[Span[UInt8, origin_of(query.data)]]()
            qseqs.append(query.seq(q))
            if self.forward:
                var qtext = SeqText(qseqs)
                total += _report(
                    out,
                    self.index[],
                    self.rtext[],
                    qtext,
                    query.names[q],
                    False,
                    self.mode,
                    self.four_column,
                    self.reference[],
                    self.name_width,
                    1,
                )
            if self.reverse:
                var rc = SeqText(qseqs, reverse_complement=True)
                total += _report(
                    out,
                    self.index[],
                    self.rtext[],
                    rc,
                    query.names[q],
                    True,
                    self.mode,
                    self.four_column,
                    self.reference[],
                    self.name_width,
                    1,
                )
        self.texts[unsafe_offset=index] = out.take()
        self.totals[unsafe_offset=index] = total


@fieldwise_init
struct _MemFilter(Copyable, Movable):
    """Which MEMs `mem` prints: the mode, and the most (0: no limit) and
    fewest times their string may occur (in the reference; for `-mum`, the
    most also in the query)."""

    var kind: String
    var max_occ: Int
    var min_occ: Int


def _report(
    mut out: OutputBuffer,
    index: MemIndex,
    rtext: SeqText,
    qtext: SeqText,
    name: String,
    reverse: Bool,
    mode: _MemFilter,
    four_column: Bool,
    reference: FastaRecords,
    name_width: Int,
    threads: Int,
) -> Int:
    var mems = find_mems(index, rtext, qtext, threads)
    var hits: List[Mem]
    if mode.kind == "mum":
        hits = filter_by_occurrences(mems^, qtext, max(mode.max_occ, 1), max(mode.max_occ, 1), mode.min_occ, threads)

        # MUMmer prints -mum matches in reference order.
        def by_ref(a: Mem, b: Mem) -> Bool:
            return a.ref_pos < b.ref_pos

        sort(Span(hits), by_ref)
    elif mode.kind == "mumreference":
        hits = filter_by_occurrences(mems^, qtext, max(mode.max_occ, 1), 0, mode.min_occ, threads)
    else:
        # -maxmatch: every MEM, unless limits were asked for.
        hits = filter_by_occurrences(mems^, qtext, mode.max_occ, 0, mode.min_occ, threads)
    out.write("> ")
    out.write(name)
    if reverse:
        out.write(" Reverse")
    out.write_byte(ord("\n"))
    # Format rounds of slices in parallel and write each round in order, so
    # the formatted text in memory stays small.
    comptime per = 1 << 16
    var slices = ceildiv(len(hits), per)
    var round = threads * 4
    var texts = List[List[UInt8]](length=min(round, max(slices, 1)), fill=List[UInt8]())
    for first in range(0, slices, round):
        var count = min(round, slices - first)
        var task = _FormatTask(
            Pointer(to=hits).as_imm(),
            Pointer(to=rtext),
            Pointer(to=reference),
            texts.unsafe_ptr(),
            first,
            per,
            four_column,
            name_width,
        )
        parallel_for(task, count, threads)
        for t in range(count):
            out.write_bytes(Span(texts[t]))
    return len(hits)


struct _FormatTask[hit_origin: ImmOrigin, text_origin: ImmOrigin, fasta_origin: ImmOrigin, out_origin: MutOrigin](
    ParallelTask
):
    """Formats slice `first + index` (`per` matches each) as MUMmer lines.

    Tracked origins keep the inputs alive while the task runs."""

    var hits: ImmPointer[List[Mem], Self.hit_origin]
    var rtext: ImmPointer[SeqText, Self.text_origin]
    var reference: ImmPointer[FastaRecords, Self.fasta_origin]
    var texts: MutPointer[List[UInt8], Self.out_origin]
    var first: Int
    var per: Int
    var four_column: Bool
    var name_width: Int

    def __init__(
        out self,
        hits: ImmPointer[List[Mem], Self.hit_origin],
        rtext: ImmPointer[SeqText, Self.text_origin],
        reference: ImmPointer[FastaRecords, Self.fasta_origin],
        texts: MutPointer[List[UInt8], Self.out_origin],
        first: Int,
        per: Int,
        four_column: Bool,
        name_width: Int,
    ):
        self.hits = hits
        self.rtext = rtext
        self.reference = reference
        self.texts = texts
        self.first = first
        self.per = per
        self.four_column = four_column
        self.name_width = name_width

    def run(self, index: Int):
        ref hits = self.hits[]
        var lo = (self.first + index) * self.per
        var hi = min(len(hits), lo + self.per)
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
