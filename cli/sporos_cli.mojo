"""The `sporos` command line.

    sporos mem [-mum | -mumreference | -maxmatch] [-l L] [-b | -r] [-F] [-n]
               [-t threads] [-v] <reference.fa> <query.fa>
    sporos minimizers -k K -w W [--forward] [-t threads] <seqs.fa>
    sporos syncmers -k K -s S [--open T[,T...]] [--forward] [-t threads] <seqs.fa>

`sporos mem` takes MUMmer's `mummer` flags and prints its output format;
`-mumreference` is the default, as in MUMmer. Non-ACGT bases never match
(MUMmer's `-n`, accepted and implied).
"""

from std.ffi import external_call
from std.math import ceildiv
from std.memory import unsafe_memcpy
from std.sys import argv, stderr
from std.sys.info import CompilationTarget
from std.time import perf_counter_ns

from sporos.batch import SeedBatch, minimizers_batch, syncmers_batch
from sporos.fasta import FastaRecords, read_fasta
from sporos.index import ALL_KMERS, MINIMIZER_SEEDS, SYNCMER_SEEDS, SeedIndex
from sporos.mems import Mem, MemIndex, MemParams, Smems, filter_mums, find_mems, find_smems
from sporos.minimizers import minimizer_positions
from sporos.output import OutputBuffer, _PAIRS, put_padded
from sporos.parallel import ParallelTask, default_threads, parallel_for
from sporos.syncmers import syncmer_positions
from sporos.text import SeqText
from sporos.window import LANES

comptime USAGE = """usage:
  sporos seeds -k K [-w W | --syncmers -s S [--open T[,T...]] | --all-kmers] [--query-kmers]
               [--max-occ N | --keep F] [-t threads] [-v] <reference.fa> <reads.fa>
  sporos mem [-mum | -mumreference | -maxmatch] [-l L] [-b | -r] [-F] [-n] [-t threads] [-v]
             [-K anchor_k] [-k1 ref_step] [-k2 query_step] <reference.fa> <query.fa>
  sporos smem [-l L] [-w W] [-t threads] [-v] <reference.fa> <query.fa>
  sporos minimizers -k K -w W [--forward] [--count] [-t threads] <seqs.fa>
  sporos syncmers -k K -s S [--open T[,T...]] [--forward] [--count] [-t threads] <seqs.fa>"""


def _log(verbose: Bool, message: String):
    if verbose:
        print(message, _peak_memory(), file=stderr)


def _peak_memory() -> String:
    """`[peak N MB]`: the process's peak resident memory (getrusage)."""
    # struct rusage: two timevals, then ru_maxrss (KB on Linux, bytes on macOS).
    var usage = List[Int64](length=18, fill=0)
    if external_call["getrusage", Int32](Int32(0), usage.unsafe_ptr()) != 0:
        return String()
    var peak = Int(usage[4])
    comptime if CompilationTarget.is_macos():
        peak //= 1024
    return String(t"[peak {peak // 1024} MB]")


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
            mode = String(a.removeprefix("-"))
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
            threads = Int(_value(args, i))
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
    var mode: String
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
        mode: String,
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
        self.mode = mode
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
        hits = filter_mums(mems, qtext, threads=threads)

        # MUMmer prints -mum matches in reference order.
        def by_ref(a: Mem, b: Mem) -> Bool:
            return a.ref_pos < b.ref_pos

        sort(Span(hits), by_ref)
    else:
        hits = filter_mums(mems, qtext, unique_in_query=False, threads=threads)
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


struct _SeedFormatTask[seeds: ImmOrigin, names: ImmOrigin, out_origin: MutOrigin](ParallelTask):
    """Formats one slice of seed positions as `name<TAB>position` lines.

    Tracked origins keep the positions alive while the task runs."""

    var positions: ImmPointer[List[UInt32], Self.seeds]
    var name: ImmPointer[String, Self.names]
    var texts: MutPointer[List[UInt8], Self.out_origin]
    var first: Int
    var per: Int

    def __init__(
        out self,
        positions: ImmPointer[List[UInt32], Self.seeds],
        name: ImmPointer[String, Self.names],
        texts: MutPointer[List[UInt8], Self.out_origin],
        first: Int,
        per: Int,
    ):
        self.positions = positions
        self.name = name
        self.texts = texts
        self.first = first
        self.per = per

    def run(self, index: Int):
        ref positions = self.positions[]
        var lo = (self.first + index) * self.per
        var hi = min(len(positions), lo + self.per)
        ref text = self.texts[unsafe_offset=index]
        if lo >= hi:
            text.resize(unsafe_uninit_length=0)
            return
        # The line prefix ("name<TAB>") padded for 16-byte copies; the buffer
        # keeps 16 bytes of slack for the last copy's overshoot.
        var name = self.name[].as_bytes()
        var prefix_len = len(name) + 1
        var prefix = List[UInt8](length=prefix_len + 16, fill=0)
        unsafe_memcpy(dest=prefix.unsafe_ptr(), src=name.unsafe_ptr(), count=len(name))
        prefix[len(name)] = UInt8(ord("\t"))
        text.resize(unsafe_uninit_length=(hi - lo) * (prefix_len + 21) + 16)
        var p = text.unsafe_ptr()
        # A position prints as its high digits (value // 10000, kept as text and
        # rebuilt only when it changes) and four low digits from the pair table.
        var high_text = List[UInt8](length=32, fill=0)
        var h = high_text.unsafe_ptr()
        var table = _PAIRS.as_bytes().unsafe_ptr()
        var high = -1
        var high_width = 0
        var at = 0
        for i in range(lo, hi):
            var value = Int(positions[i])
            var q = prefix.unsafe_ptr()
            var copied = 0
            while copied < prefix_len:
                p.unsafe_offset(at + copied).unsafe_store(q.unsafe_offset(copied).unsafe_load[width=16]())
                copied += 16
            at += prefix_len
            if value < 10000:
                at += put_padded(p.unsafe_offset(at), value, 0)
            else:
                var top = value // 10000
                var low = value - 10000 * top
                if top != high:
                    high = top
                    high_width = put_padded(h, top, 0)
                p.unsafe_offset(at).unsafe_store(h.unsafe_load[width=16]())
                at += high_width
                var upper = low // 100
                var lower = low - 100 * upper
                p[unsafe_offset=at] = table[unsafe_offset=2 * upper]
                p[unsafe_offset=at + 1] = table[unsafe_offset=2 * upper + 1]
                p[unsafe_offset=at + 2] = table[unsafe_offset=2 * lower]
                p[unsafe_offset=at + 3] = table[unsafe_offset=2 * lower + 1]
                at += 4
            p[unsafe_offset=at] = UInt8(ord("\n"))
            at += 1
        text.resize(unsafe_uninit_length=at)


comptime _BATCH_FORMAT_READS = 1 << 12
"""Reads per batch-formatting task."""


struct _BatchFormatTask[
    batch_origin: ImmOrigin, names_origin: ImmOrigin, bounds_origin: ImmOrigin, out_origin: MutOrigin
](ParallelTask):
    """Formats reads `[bounds[first + index], bounds[first + index + 1])` of
    a batch as `name<TAB>position` lines."""

    var batch: ImmPointer[SeedBatch, Self.batch_origin]
    var names: ImmPointer[List[String], Self.names_origin]
    var bounds: ImmPointer[List[Int], Self.bounds_origin]
    var texts: MutPointer[List[UInt8], Self.out_origin]
    var first: Int

    def __init__(
        out self,
        batch: ImmPointer[SeedBatch, Self.batch_origin],
        names: ImmPointer[List[String], Self.names_origin],
        bounds: ImmPointer[List[Int], Self.bounds_origin],
        texts: MutPointer[List[UInt8], Self.out_origin],
        first: Int,
    ):
        self.batch = batch
        self.names = names
        self.bounds = bounds
        self.texts = texts
        self.first = first

    def run(self, index: Int):
        ref batch = self.batch[]
        ref names = self.names[]
        var lo = self.bounds[][self.first + index]
        var hi = self.bounds[][self.first + index + 1]
        ref text = self.texts[unsafe_offset=index]
        text.resize(unsafe_uninit_length=0)
        var prefix = List[UInt8]()
        var at = 0
        for r in range(lo, hi):
            var seeds = batch.seeds(r)
            if len(seeds) == 0:
                continue
            # "name<TAB>", padded for 16-byte copies.
            var name = names[r].as_bytes()
            var prefix_len = len(name) + 1
            prefix.resize(unsafe_uninit_length=ceildiv(prefix_len, 16) * 16)
            unsafe_memcpy(dest=prefix.unsafe_ptr(), src=name.unsafe_ptr(), count=len(name))
            prefix[len(name)] = UInt8(ord("\t"))
            var room = at + len(seeds) * (prefix_len + 11) + 16
            if room > text.capacity():
                text.reserve(max(room, 2 * text.capacity()))
            text.resize(unsafe_uninit_length=room)
            var p = text.unsafe_ptr()
            var q = prefix.unsafe_ptr()
            for position in seeds:
                var copied = 0
                while copied < prefix_len:
                    p.unsafe_offset(at + copied).unsafe_store(q.unsafe_offset(copied).unsafe_load[width=16]())
                    copied += 16
                at += prefix_len
                at += put_padded(p.unsafe_offset(at), Int(position), 0)
                p[unsafe_offset=at] = UInt8(ord("\n"))
                at += 1
        text.resize(unsafe_uninit_length=at)


comptime _FORMAT_SLICE = 1 << 16
"""Positions per formatting task."""


def _value(args: List[String], i: Int) raises -> String:
    """The value of the option before `args[i]`."""
    if i >= len(args):
        raise Error(t"{args[i - 1]} needs a value\n{USAGE}")
    return args[i]


def _int_list(arg: String) raises -> List[Int]:
    """The integers of a comma-separated list."""
    var out = List[Int]()
    for field in arg.split(","):
        out.append(Int(field))
    return out^


def run_seeds(command: String, args: List[String]) raises:
    var k = 15
    var w = 10
    var s = 0
    var open_offsets = List[Int]()
    var canonical = True
    var count_only = False
    var threads = default_threads()
    var paths = List[String]()
    var i = 0
    while i < len(args):
        var a = args[i]
        if a == "--count":
            count_only = True
        elif a == "-k":
            i += 1
            k = Int(_value(args, i))
        elif a == "-w":
            i += 1
            w = Int(_value(args, i))
        elif a == "-s":
            i += 1
            s = Int(_value(args, i))
        elif a == "--open":
            i += 1
            open_offsets = _int_list(_value(args, i))
        elif a == "--forward":
            canonical = False
        elif a == "-t" or a == "-threads":
            i += 1
            threads = Int(_value(args, i))
        elif a.startswith("-"):
            raise Error(t"unknown option {a}\n{USAGE}")
        else:
            paths.append(a)
        i += 1
    if len(paths) != 1:
        raise Error(USAGE)
    var start = perf_counter_ns()
    var records = read_fasta(paths[0], threads)
    var read_seconds = _seconds(start)
    start = perf_counter_ns()
    if command == "syncmers" and s == 0:
        raise Error("syncmers need -s")
    # Many records (reads): scan them in batches, one read per SIMD lane.
    var batched = len(records) >= 4 * LANES
    var batch = SeedBatch(List[UInt32](), [0])
    if batched:
        var spans = List[Span[UInt8, origin_of(records.data)]](capacity=len(records))
        for r in range(len(records)):
            spans.append(records.seq(r))
        if command == "minimizers":
            batch = minimizers_batch(spans, k, w, canonical, threads)
        else:
            batch = syncmers_batch(spans, k, s, canonical, None, threads, open_offsets=open_offsets)
    var seeds = 0
    var out = OutputBuffer()
    var round = max(1, threads * 4)
    var texts = List[List[UInt8]](length=round, fill=List[UInt8]())
    if batched:
        # Many reads: format runs of reads straight from the batch, in
        # parallel rounds.
        seeds = len(batch.positions)
        if not count_only:
            # Tasks of up to _BATCH_FORMAT_READS reads or _FORMAT_SLICE positions.
            var bounds: List[Int] = [0]
            for r in range(len(records)):
                var start = bounds[len(bounds) - 1]
                if r + 1 - start >= _BATCH_FORMAT_READS or batch.offsets[r + 1] - batch.offsets[start] >= _FORMAT_SLICE:
                    bounds.append(r + 1)
            if bounds[len(bounds) - 1] != len(records):
                bounds.append(len(records))
            var tasks = len(bounds) - 1
            for first in range(0, tasks, round):
                var count = min(round, tasks - first)
                var task = _BatchFormatTask(
                    Pointer(to=batch).as_imm(),
                    Pointer(to=records.names).as_imm(),
                    Pointer(to=bounds).as_imm(),
                    texts.unsafe_ptr(),
                    first,
                )
                parallel_for(task, count, threads)
                for t in range(count):
                    out.write_bytes(Span(texts[t]))
    for r in range(len(records) if not batched else 0):
        var seq = records.seq(r)
        var positions: List[UInt32]
        if command == "minimizers":
            positions = minimizer_positions(seq, k, w, canonical, threads)
        else:
            positions = syncmer_positions(seq, k, s, canonical, None, threads, open_offsets=open_offsets)
        seeds += len(positions)
        if count_only:
            continue
        # Rounds of slices formatted in parallel into reused buffers, each
        # round written before the next: little formatted text is held.
        var slices = ceildiv(len(positions), _FORMAT_SLICE)
        for first in range(0, slices, round):
            var count = min(round, slices - first)
            var task = _SeedFormatTask(
                Pointer(to=positions).as_imm(),
                Pointer(to=records.names[r]).as_imm(),
                texts.unsafe_ptr(),
                first,
                _FORMAT_SLICE,
            )
            parallel_for(task, count, threads)
            for t in range(count):
                out.write_bytes(Span(texts[t]))
    if count_only:
        out.write(String(t"{records.total_length()} bases\t{seeds} {command}\n"))
        print(
            t"read {read_seconds} s, seeds {_seconds(start)} s ({len(records)} sequences)",
            file=stderr,
        )
    out.flush()


struct _SeedLookupTask[index_origin: ImmOrigin, read_origin: ImmOrigin, name_origin: ImmOrigin, out_origin: MutOrigin](
    ParallelTask
):
    """Looks up one slice of reads (one batch) and formats their hits as TSV lines."""

    var index: ImmPointer[SeedIndex, Self.index_origin]
    var reads: ImmPointer[FastaRecords, Self.read_origin]
    var ref_names: ImmPointer[List[String], Self.name_origin]
    var texts: MutPointer[List[UInt8], Self.out_origin]
    var bounds: ImmPointer[Int, ImmUntrackedOrigin]
    var first: Int
    var all_kmers: Bool
    var max_occ: Int

    def __init__(
        out self,
        index: ImmPointer[SeedIndex, Self.index_origin],
        reads: ImmPointer[FastaRecords, Self.read_origin],
        ref_names: ImmPointer[List[String], Self.name_origin],
        texts: MutPointer[List[UInt8], Self.out_origin],
        bounds: ImmPointer[Int, ImmUntrackedOrigin],
        first: Int,
        all_kmers: Bool,
        max_occ: Int,
    ):
        self.index = index
        self.reads = reads
        self.ref_names = ref_names
        self.texts = texts
        self.bounds = bounds
        self.first = first
        self.all_kmers = all_kmers
        self.max_occ = max_occ

    def run(self, index: Int):
        ref reads = self.reads[]
        ref seed_index = self.index[]
        var lo = self.bounds[unsafe_offset=self.first + index]
        var hi = self.bounds[unsafe_offset=self.first + index + 1]
        ref text = self.texts[unsafe_offset=index]
        text.resize(unsafe_uninit_length=0)
        if lo >= hi:
            return
        var spans = List[Span[UInt8, origin_of(reads.data)]](capacity=hi - lo)
        for r in range(lo, hi):
            spans.append(reads.seq(r))
        var hits = seed_index.lookup_batch(spans, self.all_kmers, self.max_occ)
        # Name fields ("name<TAB>") padded for 16-byte copies.
        var ref_fields = List[List[UInt8]]()
        var ref_room = 0
        for name in self.ref_names[]:
            ref_fields.append(_padded_field(name))
            ref_room = max(ref_room, name.byte_length() + 1)
        var at = 0
        for r in range(lo, hi):
            var read_hits = hits.of(r - lo)
            if len(read_hits) == 0:
                continue
            var read_field = _padded_field(reads.names[r])
            var read_room = reads.names[r].byte_length() + 1
            # Room for every line (two fields, three numbers, four more bytes)
            # and the 16-byte copies' overshoot.
            var room = at + len(read_hits) * (read_room + ref_room + 3 * 10 + 4) + 16
            if room > text.capacity():
                text.reserve(max(room, 2 * text.capacity()))
            text.resize(unsafe_uninit_length=room)
            var p = text.unsafe_ptr()
            for h in read_hits:
                var ref_pos = Int(h.ref_pos)
                var seq = seed_index.text.sequence_of(ref_pos) if len(ref_fields) > 1 else 0
                at += _put_field(p.unsafe_offset(at), read_field, read_room)
                at += put_padded(p.unsafe_offset(at), Int(h.query_pos), 0)
                p[unsafe_offset=at] = UInt8(ord("\t"))
                at += 1
                at += _put_field(p.unsafe_offset(at), ref_fields[seq], self.ref_names[][seq].byte_length() + 1)
                at += put_padded(p.unsafe_offset(at), ref_pos - seed_index.text.starts[seq], 0)
                p[unsafe_offset=at] = UInt8(ord("\t"))
                p[unsafe_offset=at + 1] = UInt8(ord("-")) if h.reverse else UInt8(ord("+"))
                p[unsafe_offset=at + 2] = UInt8(ord("\t"))
                at += 3
                at += put_padded(p.unsafe_offset(at), Int(h.occurrences), 0)
                p[unsafe_offset=at] = UInt8(ord("\n"))
                at += 1
        text.resize(unsafe_uninit_length=at)


comptime _SEED_SLICE = 1 << 12
comptime _SEED_SLICE_BASES = 1 << 18
"""A seed-lookup task's reads: up to this many reads or bases."""


def _padded_field(name: String) -> List[UInt8]:
    """`name<TAB>`, zero-padded to a multiple of 16 bytes for `_put_field`."""
    var bytes = name.as_bytes()
    var field = List[UInt8](length=ceildiv(len(bytes) + 1, 16) * 16, fill=0)
    unsafe_memcpy(dest=field.unsafe_ptr(), src=bytes.unsafe_ptr(), count=len(bytes))
    field[len(bytes)] = UInt8(ord("\t"))
    return field^


@inline(.always)
def _put_field(dst: MutPointer[UInt8, _], field: List[UInt8], length: Int) -> Int:
    """Copies a `_padded_field` in 16-byte stores (up to 15 bytes past
    `length` are scratch); returns `length`."""
    var src = field.unsafe_ptr()
    var copied = 0
    while copied < length:
        dst.unsafe_offset(copied).unsafe_store(src.unsafe_offset(copied).unsafe_load[width=16]())
        copied += 16
    return length


def run_seed_index(args: List[String]) raises:
    var k = 19
    var w = 12
    var s = 0
    var open_offsets = List[Int]()
    var scheme = MINIMIZER_SEEDS
    var all_kmers = False
    var max_occ = 0
    var keep = 0.0
    var threads = default_threads()
    var verbose = False
    var paths = List[String]()
    var i = 0
    while i < len(args):
        var a = args[i]
        if a == "-k":
            i += 1
            k = Int(_value(args, i))
        elif a == "-w":
            i += 1
            w = Int(_value(args, i))
        elif a == "--keep":
            i += 1
            keep = Float64(_value(args, i))
        elif a == "-s":
            i += 1
            s = Int(_value(args, i))
        elif a == "--open":
            i += 1
            open_offsets = _int_list(_value(args, i))
        elif a == "--syncmers":
            scheme = SYNCMER_SEEDS
        elif a == "--all-kmers":
            scheme = ALL_KMERS
        elif a == "--query-kmers":
            all_kmers = True
        elif a == "--max-occ":
            i += 1
            max_occ = Int(_value(args, i))
        elif a == "-t" or a == "-threads":
            i += 1
            threads = Int(_value(args, i))
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
    var reference = read_fasta(paths[0], threads)
    var reads = read_fasta(paths[1], threads)
    _log(verbose, String(t"read inputs: {_seconds(start)} s"))
    start = perf_counter_ns()
    var seqs = List[Span[UInt8, origin_of(reference.data)]]()
    for r in range(len(reference)):
        seqs.append(reference.seq(r))
    if scheme == SYNCMER_SEEDS:
        if s == 0:
            raise Error("--syncmers needs -s (the s-mer length)")
        w = s  # SeedIndex takes the syncmers' s-mer length as `w`
    elif len(open_offsets) > 0:
        raise Error("--open applies to --syncmers")
    var index = SeedIndex(seqs, k, w, scheme, threads, open_offsets=open_offsets)
    if keep > 0.0:
        max_occ = index.occurrence_cutoff(keep)
    _log(verbose, String(t"index: {_seconds(start)} s ({len(index.table)} seeds, max occurrences {max_occ})"))
    start = perf_counter_ns()
    # Rounds of read slices looked up and formatted in parallel into reused
    # buffers, each round written before the next.
    var bounds: List[Int] = [0]
    var bases = 0
    for r in range(len(reads)):
        bases += reads.lengths[r]
        if bases >= _SEED_SLICE_BASES or r + 1 - bounds[len(bounds) - 1] >= _SEED_SLICE:
            bounds.append(r + 1)
            bases = 0
    if bounds[len(bounds) - 1] != len(reads):
        bounds.append(len(reads))
    var slices = len(bounds) - 1
    var round = max(1, threads * 4)
    var texts = List[List[UInt8]](length=round, fill=List[UInt8]())
    var out = OutputBuffer()
    for first in range(0, slices, round):
        var count = min(round, slices - first)
        var task = _SeedLookupTask(
            Pointer(to=index).as_imm(),
            Pointer(to=reads).as_imm(),
            Pointer(to=reference.names).as_imm(),
            texts.unsafe_ptr(),
            bounds.unsafe_ptr().as_imm().unsafe_origin_cast[ImmUntrackedOrigin](),
            first,
            all_kmers,
            max_occ,
        )
        parallel_for(task, count, threads)
        for t in range(count):
            out.write_bytes(Span(texts[t]))
    out.flush()
    _ = len(bounds)  # read by the tasks through an untracked pointer
    _log(verbose, String(t"lookups and output: {_seconds(start)} s"))


def run_smem(args: List[String]) raises:
    """SMEMs of every query record against the reference on both strands, in
    `bwa fastmap`'s output format."""
    var min_len = 17
    var max_hits = 20
    var threads = default_threads()
    var verbose = False
    var paths = List[String]()
    var i = 0
    while i < len(args):
        var a = args[i]
        if a == "-l":
            i += 1
            min_len = Int(_value(args, i))
        elif a == "-w":
            i += 1
            max_hits = Int(_value(args, i))
        elif a == "-t" or a == "-threads":
            i += 1
            threads = Int(_value(args, i))
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
    var reference = read_fasta(paths[0], threads)
    var query = read_fasta(paths[1], threads)
    _log(verbose, String(t"read inputs: {_seconds(start)} s"))
    start = perf_counter_ns()
    var ref_seqs = List[Span[UInt8, origin_of(reference.data)]]()
    for r in range(len(reference)):
        ref_seqs.append(reference.seq(r))
    var rtext = SeqText(ref_seqs)
    for q in range(len(query)):
        if len(query.seq(q)) >= Int(UInt32.MAX):
            raise Error("query sequences must be shorter than 2^32 bases")
    var index = MemIndex(rtext, MemParams.choose(min_len, len(rtext), query.total_length()), threads)
    _log(verbose, String(t"index: {_seconds(start)} s ({len(index.table)} anchors)"))
    start = perf_counter_ns()
    var out = OutputBuffer()
    # One record per slice task (single-threaded each) for many records; a
    # single record uses every thread inside the matcher.
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
    var inner = threads if len(query) == 1 else 1
    var total = 0
    for first in range(0, slices, round):
        var count = min(round, slices - first)
        var task = _SmemTask(
            Pointer(to=index).as_imm(),
            Pointer(to=rtext).as_imm(),
            Pointer(to=query).as_imm(),
            Pointer(to=reference.names).as_imm(),
            Pointer(to=bounds).as_imm(),
            texts.unsafe_ptr(),
            first,
            max_hits,
            inner,
        )
        parallel_for(task, count, threads if inner == 1 else 1)
        for s in range(count):
            out.write_bytes(Span(texts[s]))
    out.flush()
    _log(verbose, String(t"smems: {_seconds(start)} s"))


struct _SmemTask[
    index_origin: ImmOrigin,
    text_origin: ImmOrigin,
    query_origin: ImmOrigin,
    names_origin: ImmOrigin,
    bounds_origin: ImmOrigin,
    out_origin: MutOrigin,
](ParallelTask):
    """Finds and formats (as `bwa fastmap`) the SMEMs of one slice of query
    records."""

    var index: ImmPointer[MemIndex, Self.index_origin]
    var rtext: ImmPointer[SeqText, Self.text_origin]
    var query: ImmPointer[FastaRecords, Self.query_origin]
    var names: ImmPointer[List[String], Self.names_origin]
    var bounds: ImmPointer[List[Int], Self.bounds_origin]
    var texts: MutPointer[List[UInt8], Self.out_origin]
    var first: Int
    var max_hits: Int
    var threads: Int

    def __init__(
        out self,
        index: ImmPointer[MemIndex, Self.index_origin],
        rtext: ImmPointer[SeqText, Self.text_origin],
        query: ImmPointer[FastaRecords, Self.query_origin],
        names: ImmPointer[List[String], Self.names_origin],
        bounds: ImmPointer[List[Int], Self.bounds_origin],
        texts: MutPointer[List[UInt8], Self.out_origin],
        first: Int,
        max_hits: Int,
        threads: Int,
    ):
        self.index = index
        self.rtext = rtext
        self.query = query
        self.names = names
        self.bounds = bounds
        self.texts = texts
        self.first = first
        self.max_hits = max_hits
        self.threads = threads

    def run(self, index: Int):
        ref query = self.query[]
        ref rtext = self.rtext[]
        var lo = self.bounds[][self.first + index]
        var hi = self.bounds[][self.first + index + 1]
        var out = OutputBuffer(fd=-1)
        for q in range(lo, hi):
            var qseqs = List[Span[UInt8, origin_of(query.data)]]()
            qseqs.append(query.seq(q))
            var qtext = SeqText(qseqs)
            var qrc = SeqText(qseqs, reverse_complement=True)
            var found = find_smems(self.index[], rtext, qtext, qrc, self.threads)
            out.write("SQ\t")
            out.write(query.names[q])
            out.write_byte(ord("\t"))
            out.write_int(len(query.seq(q)))
            out.write_byte(ord("\n"))
            for s in range(len(found)):
                var smem = found.smems[s]
                out.write("EM\t")
                out.write_int(Int(smem.start))
                out.write_byte(ord("\t"))
                out.write_int(Int(smem.end))
                out.write_byte(ord("\t"))
                out.write_int(smem.count)
                if smem.count > self.max_hits:
                    out.write("\t*\n")  # bwa fastmap follows these lines with an empty one
                else:
                    for h in found.hits_of(s):
                        var pos = Int(h.ref_pos)
                        var seq = rtext.sequence_of(pos)
                        out.write_byte(ord("\t"))
                        out.write(self.names[][seq])
                        out.write(":-" if h.reverse else ":+")
                        out.write_int(pos - rtext.starts[seq] + 1)
                out.write_byte(ord("\n"))
            out.write("//\n")
        self.texts[unsafe_offset=index] = out.take()


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
    elif command == "seeds":
        run_seed_index(args)
    elif command == "smem":
        run_smem(args)
    elif command == "minimizers" or command == "syncmers":
        run_seeds(command, args)
    else:
        print(USAGE)
