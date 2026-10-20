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

from sporos.fasta import FastaRecords, read_fasta
from sporos.index import ALL_KMERS, MINIMIZER_SEEDS, SYNCMER_SEEDS, SeedHit, SeedIndex
from sporos.mems import Mem, MemIndex, MemParams, filter_mums, find_mems
from sporos.minimizers import minimizer_positions
from sporos.output import OutputBuffer, put_padded
from sporos.parallel import ParallelTask, default_threads, parallel_for
from sporos.syncmers import syncmer_positions
from sporos.text import SeqText

comptime USAGE = """usage:
  sporos seeds -k K -w W [--syncmers | --all-kmers] [--query-kmers] [--max-occ N | --keep F] [-t threads]
               <reference.fa> <reads.fa>
  sporos mem [-mum | -mumreference | -maxmatch] [-l L] [-b | -r] [-F] [-n] [-t threads] [-v]
             [-K anchor_k] [-k1 ref_step] [-k2 query_step] <reference.fa> <query.fa>
  sporos minimizers -k K -w W [--forward] [--count] [-t threads] <seqs.fa>
  sporos syncmers -k K -s S [--open T] [--forward] [--count] [-t threads] <seqs.fa>"""


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
    var anchor_k = 0
    var ref_step = 0
    var query_step = 0
    var paths = List[String]()
    var i = 0
    while i < len(args):
        var a = args[i]
        if a == "-K":
            i += 1
            anchor_k = Int(args[i])
        elif a == "-k1":
            i += 1
            ref_step = Int(args[i])
        elif a == "-k2":
            i += 1
            query_step = Int(args[i])
        elif a == "-mum" or a == "-mumreference" or a == "-maxmatch":
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
    if threads > 1 and len(query) >= 2 * threads:
        # Many query records (reads): one thread per slice of records.
        var slices = min(len(query), threads * 8)
        var texts = List[List[UInt8]](length=slices, fill=List[UInt8]())
        var totals = List[Int](length=slices, fill=0)
        var task = _RecordsTask(
            Pointer(to=index).as_imm(),
            Pointer(to=rtext).as_imm(),
            Pointer(to=query).as_imm(),
            Pointer(to=reference).as_imm(),
            texts.unsafe_ptr(),
            totals.unsafe_ptr(),
            slices,
            mode,
            forward,
            reverse,
            four_column,
            name_width,
        )
        parallel_for(task, slices, threads)
        for s in range(slices):
            out.write_bytes(Span(texts[s]))
            total += totals[s]
        out.flush()
        _log(verbose, String(t"matches: {_seconds(query_start)} s ({total} reported)"))
        return
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


struct _RecordsTask[
    index_origin: ImmOrigin,
    text_origin: ImmOrigin,
    query_origin: ImmOrigin,
    ref_origin: ImmOrigin,
    out_origin: MutOrigin,
    total_origin: MutOrigin,
](ParallelTask):
    """Reports the matches of one slice of query records, single-threaded."""

    var index: ImmPointer[MemIndex, Self.index_origin]
    var rtext: ImmPointer[SeqText, Self.text_origin]
    var query: ImmPointer[FastaRecords, Self.query_origin]
    var reference: ImmPointer[FastaRecords, Self.ref_origin]
    var texts: MutPointer[List[UInt8], Self.out_origin]
    var totals: MutPointer[Int, Self.total_origin]
    var slices: Int
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
        texts: MutPointer[List[UInt8], Self.out_origin],
        totals: MutPointer[Int, Self.total_origin],
        slices: Int,
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
        self.texts = texts
        self.totals = totals
        self.slices = slices
        self.mode = mode
        self.forward = forward
        self.reverse = reverse
        self.four_column = four_column
        self.name_width = name_width

    def run(self, index: Int):
        ref query = self.query[]
        var per = ceildiv(len(query), self.slices)
        var lo = index * per
        var hi = min(len(query), lo + per)
        var out = OutputBuffer(fd=-1)
        var total = 0
        for q in range(lo, hi):
            var qseqs = List[Span[UInt8, origin_of(query.data)]]()
            qseqs.append(query.seq(q))
            if self.forward:
                var qtext = SeqText(qseqs)
                total += _report(
                    out, self.index[], self.rtext[], qtext, query.names[q], False, self.mode,
                    self.four_column, self.reference[], self.name_width, 1,
                )
            if self.reverse:
                var rc = SeqText(qseqs, reverse_complement=True)
                total += _report(
                    out, self.index[], self.rtext[], rc, query.names[q], True, self.mode,
                    self.four_column, self.reference[], self.name_width, 1,
                )
        self.texts[unsafe_offset=index] = out.buf.copy()
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
        Pointer(to=hits).as_imm(),
        Pointer(to=rtext),
        Pointer(to=reference),
        texts.unsafe_ptr(),
        slices,
        four_column,
        name_width,
    )
    parallel_for(task, slices, threads)
    for text in texts:
        out.write_bytes(Span(text))
    return len(hits)


struct _FormatTask[
    hit_origin: ImmOrigin, text_origin: ImmOrigin, fasta_origin: ImmOrigin, out_origin: MutOrigin
](ParallelTask):
    """Formats one slice of matches as MUMmer lines.

    Tracked origins keep the inputs alive while the task runs."""

    var hits: ImmPointer[List[Mem], Self.hit_origin]
    var rtext: ImmPointer[SeqText, Self.text_origin]
    var reference: ImmPointer[FastaRecords, Self.fasta_origin]
    var texts: MutPointer[List[UInt8], Self.out_origin]
    var slices: Int
    var four_column: Bool
    var name_width: Int

    def __init__(
        out self,
        hits: ImmPointer[List[Mem], Self.hit_origin],
        rtext: ImmPointer[SeqText, Self.text_origin],
        reference: ImmPointer[FastaRecords, Self.fasta_origin],
        texts: MutPointer[List[UInt8], Self.out_origin],
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


struct _SeedFormatTask[
    seeds: ImmOrigin, names: ImmOrigin, out_origin: MutOrigin
](ParallelTask):
    """Formats one slice of seed positions as `name<TAB>position` lines.

    Tracked origins keep the positions alive while the task runs."""

    var positions: ImmPointer[List[UInt32], Self.seeds]
    var name: ImmPointer[String, Self.names]
    var texts: MutPointer[List[UInt8], Self.out_origin]
    var slices: Int

    def __init__(
        out self,
        positions: ImmPointer[List[UInt32], Self.seeds],
        name: ImmPointer[String, Self.names],
        texts: MutPointer[List[UInt8], Self.out_origin],
        slices: Int,
    ):
        self.positions = positions
        self.name = name
        self.texts = texts
        self.slices = slices

    def run(self, index: Int):
        ref positions = self.positions[]
        var per = ceildiv(len(positions), self.slices)
        var lo = index * per
        var hi = min(len(positions), lo + per)
        var name = self.name[].as_bytes()
        ref text = self.texts[unsafe_offset=index]
        text.resize(unsafe_uninit_length=max(0, hi - lo) * (len(name) + 12))
        var p = text.unsafe_ptr()
        var at = 0
        for i in range(lo, hi):
            for b in range(len(name)):
                p[unsafe_offset=at + b] = name[b]
            at += len(name)
            p[unsafe_offset=at] = UInt8(ord("\t"))
            at += 1
            at += put_padded(p.unsafe_offset(at), Int(positions[i]), 0)
            p[unsafe_offset=at] = UInt8(ord("\n"))
            at += 1
        text.resize(unsafe_uninit_length=at)


def run_seeds(command: String, args: List[String]) raises:
    var k = 15
    var w = 10
    var s = 0
    var open_offset: Optional[Int] = None
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
    var start = perf_counter_ns()
    var records = read_fasta(paths[0])
    var read_seconds = _seconds(start)
    start = perf_counter_ns()
    var seeds = 0
    var out = OutputBuffer()
    for r in range(len(records)):
        var seq = records.seq(r)
        var positions: List[UInt32]
        if command == "minimizers":
            positions = minimizer_positions(seq, k, w, canonical, threads)
        else:
            if s == 0:
                raise Error("syncmers need -s")
            positions = syncmer_positions(seq, k, s, canonical, open_offset, threads)
        seeds += len(positions)
        if count_only:
            continue
        var slices = max(1, min(threads * 4, len(positions) // 100_000))
        var texts = List[List[UInt8]](length=slices, fill=List[UInt8]())
        var task = _SeedFormatTask(
            Pointer(to=positions).as_imm(),
            Pointer(to=records.names[r]).as_imm(),
            texts.unsafe_ptr(),
            slices,
        )
        parallel_for(task, slices, threads)
        for text in texts:
            out.write_bytes(Span(text))
    if count_only:
        out.write(String(t"{records.total_length()} bases\t{seeds} {command}\n"))
        print(
            t"read {read_seconds} s, seeds {_seconds(start)} s ({len(records)} sequences)",
            file=stderr,
        )
    out.flush()


struct _SeedLookupTask[
    index_origin: ImmOrigin, read_origin: ImmOrigin, name_origin: ImmOrigin, out_origin: MutOrigin
](ParallelTask):
    """Looks up one slice of reads and formats their hits as TSV lines."""

    var index: ImmPointer[SeedIndex, Self.index_origin]
    var reads: ImmPointer[FastaRecords, Self.read_origin]
    var ref_names: ImmPointer[List[String], Self.name_origin]
    var texts: MutPointer[List[UInt8], Self.out_origin]
    var slices: Int
    var all_kmers: Bool
    var max_occ: Int

    def __init__(
        out self,
        index: ImmPointer[SeedIndex, Self.index_origin],
        reads: ImmPointer[FastaRecords, Self.read_origin],
        ref_names: ImmPointer[List[String], Self.name_origin],
        texts: MutPointer[List[UInt8], Self.out_origin],
        slices: Int,
        all_kmers: Bool,
        max_occ: Int,
    ):
        self.index = index
        self.reads = reads
        self.ref_names = ref_names
        self.texts = texts
        self.slices = slices
        self.all_kmers = all_kmers
        self.max_occ = max_occ

    def run(self, index: Int):
        ref reads = self.reads[]
        ref seed_index = self.index[]
        var per = ceildiv(len(reads), self.slices)
        var lo = index * per
        var hi = min(len(reads), lo + per)
        ref text = self.texts[unsafe_offset=index]
        var hits = List[SeedHit]()
        var digits = List[UInt8](length=32, fill=0)
        for r in range(lo, hi):
            hits.clear()
            try:
                seed_index.lookup(reads.seq(r), hits, self.all_kmers, self.max_occ)
            except:
                continue
            var read_name = reads.names[r].as_bytes()
            for h in hits:
                var ref_pos = Int(h.ref_pos)
                var seq = seed_index.text.sequence_of(ref_pos)
                text.extend(read_name)
                text.append(UInt8(ord("\t")))
                _append_int(text, Int(h.query_pos), digits)
                text.append(UInt8(ord("\t")))
                text.extend(self.ref_names[][seq].as_bytes())
                text.append(UInt8(ord("\t")))
                _append_int(text, ref_pos - seed_index.text.starts[seq], digits)
                text.append(UInt8(ord("\t")))
                text.append(UInt8(ord("-")) if h.reverse else UInt8(ord("+")))
                text.append(UInt8(ord("\t")))
                _append_int(text, Int(h.occurrences), digits)
                text.append(UInt8(ord("\n")))


@inline(.always)
def _append_int(mut text: List[UInt8], value: Int, mut digits: List[UInt8]):
    var n = put_padded(digits.unsafe_ptr(), value, 0)
    text.extend(Span(digits)[:n])


def run_seed_index(args: List[String]) raises:
    var k = 19
    var w = 12
    var scheme = MINIMIZER_SEEDS
    var all_kmers = False
    var max_occ = 0
    var keep = 0.0
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
        elif a == "--keep":
            i += 1
            keep = Float64(args[i])
        elif a == "--syncmers":
            scheme = SYNCMER_SEEDS
        elif a == "--all-kmers":
            scheme = ALL_KMERS
        elif a == "--query-kmers":
            all_kmers = True
        elif a == "--max-occ":
            i += 1
            max_occ = Int(args[i])
        elif a == "-t" or a == "-threads":
            i += 1
            threads = Int(args[i])
        elif a.startswith("-"):
            raise Error(t"unknown option {a}\n{USAGE}")
        else:
            paths.append(a)
        i += 1
    if len(paths) != 2:
        raise Error(USAGE)
    var reference = read_fasta(paths[0])
    var reads = read_fasta(paths[1])
    var seqs = List[Span[UInt8, origin_of(reference.data)]]()
    for r in range(len(reference)):
        seqs.append(reference.seq(r))
    var index = SeedIndex(seqs, k, w, scheme, threads)
    if keep > 0.0:
        max_occ = index.occurrence_cutoff(keep)
    var slices = max(1, min(threads * 8, len(reads)))
    var texts = List[List[UInt8]](length=slices, fill=List[UInt8]())
    var task = _SeedLookupTask(
        Pointer(to=index).as_imm(),
        Pointer(to=reads).as_imm(),
        Pointer(to=reference.names).as_imm(),
        texts.unsafe_ptr(),
        slices,
        all_kmers,
        max_occ,
    )
    parallel_for(task, slices, threads)
    var out = OutputBuffer()
    for text in texts:
        out.write_bytes(Span(text))
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
    elif command == "seeds":
        run_seed_index(args)
    elif command == "minimizers" or command == "syncmers":
        run_seeds(command, args)
    else:
        print(USAGE)
