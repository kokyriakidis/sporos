"""`sporos smem` and `sporos index`: SMEMs in `bwa fastmap`'s format, from
the FM-index (reads) or the MEM engine (one long query)."""

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


def run_index(args: List[String]) raises:
    """Builds the FM-index `sporos smem` uses for reads, by default beside the
    reference as `<reference>.sfm`. `-s` samples the suffix array every
    2^s rows (default 2 up to 2^29 bases, 3 up to 2^31, then 4): smaller is
    faster, larger is smaller."""
    var threads = default_threads()
    var verbose = False
    var sample_bits = 0
    var paths = List[String]()
    var i = 0
    while i < len(args):
        var a = args[i]
        if a == "-s":
            i += 1
            sample_bits = Int(_value(args, i))
            if sample_bits < 1 or sample_bits > 8:
                raise Error("-s takes 1 to 8")
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
    if len(paths) < 1 or len(paths) > 2:
        raise Error(USAGE)
    var out = paths[1] if len(paths) == 2 else paths[0] + ".sfm"
    var start = perf_counter_ns()
    var reference = read_fasta(paths[0], threads)
    _log(verbose, String(t"read reference: {_seconds(start)} s"))
    start = perf_counter_ns()
    var fm = _build_fm(reference, threads, sample_bits)
    _log(verbose, String(t"index: {_seconds(start)} s ({fm.rows} rows)"))
    start = perf_counter_ns()
    fm.save(out)
    _log(verbose, String(t"wrote {out}: {_seconds(start)} s"))


def _build_fm(reference: FastaRecords, threads: Int, sample_bits: Int = 0) raises -> FMIndex:
    var seqs = List[Span[UInt8, origin_of(reference.data)]]()
    for r in range(len(reference)):
        seqs.append(reference.seq(r))
    return FMIndex.build(seqs, reference.names, threads, sample_bits)


def _fm_index_path(reference: String) -> String:
    """The FM-index to use for `reference`: the path itself when it is one,
    else `<reference>.sfm` if that exists, else none."""
    if reference.endswith(".sfm"):
        return reference
    var path = reference + ".sfm"
    return path if exists(path) else String()


def run_smem(args: List[String]) raises:
    """SMEMs of every query record against the reference on both strands, in
    `bwa fastmap`'s output format. Many queries (reads) are searched in an
    FM-index, loaded from `<reference>.sfm` (`sporos index`) or built on the
    spot; a single query (a genome) goes through the MEM engine. `--engine`
    picks one."""
    var min_len = 17
    var max_hits = 20
    var min_occ = 1
    var max_occ = 0
    var engine = String("auto")
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
        elif a == "-s":
            i += 1
            min_occ = _count_option(args, i, "-s")
        elif a == "--max-occ":
            i += 1
            max_occ = _count_option(args, i, "--max-occ")
        elif a == "--engine":
            i += 1
            engine = _value(args, i)
            if engine != "fm" and engine != "mem":
                raise Error(t"--engine takes fm or mem, not {engine}")
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
    var start = perf_counter_ns()
    var index_path = _fm_index_path(paths[0]) if engine != "mem" else String()
    var query = read_fasta(paths[1], threads)
    for q in range(len(query)):
        if query.lengths[q] >= Int(UInt32.MAX):
            raise Error("query sequences must be shorter than 2^32 bases")
    var counted = min_occ > 1 or max_occ > 0
    if counted and engine == "mem":
        raise Error("-s and --max-occ need the FM engine (--engine fm)")
    if engine == "fm" or (engine == "auto" and (index_path or len(query) > 1 or counted)):
        var fm: FMIndex
        if index_path:
            fm = FMIndex.load(index_path, threads)
            _log(verbose, String(t"read inputs and index: {_seconds(start)} s"))
        else:
            var reference = read_fasta(paths[0], threads)
            _log(verbose, String(t"read inputs: {_seconds(start)} s"))
            start = perf_counter_ns()
            fm = _build_fm(reference, threads)
            _log(verbose, String(t"index: {_seconds(start)} s ({fm.rows} rows, k-mer cache k={fm.cache_k})"))
        start = perf_counter_ns()
        _fm_smems(fm, query, min_len, max_hits, min_occ, max_occ, threads)
        _log(verbose, String(t"smems: {_seconds(start)} s"))
        return
    var reference = read_fasta(paths[0], threads)
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


def _fm_smems(
    fm: FMIndex, query: FastaRecords, min_len: Int, max_hits: Int, min_occ: Int, max_occ: Int, threads: Int
) raises:
    var out = OutputBuffer()
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
    for first in range(0, slices, round):
        var count = min(round, slices - first)
        var task = _FmSmemTask(
            Pointer(to=fm).as_imm(),
            Pointer(to=query).as_imm(),
            Pointer(to=bounds).as_imm(),
            texts.unsafe_ptr(),
            first,
            min_len,
            max_hits,
            min_occ,
            max_occ,
        )
        parallel_for(task, count, threads)
        for s in range(count):
            out.write_bytes(Span(texts[s]))
    out.flush()


struct _FmSmemTask[
    index_origin: ImmOrigin,
    query_origin: ImmOrigin,
    bounds_origin: ImmOrigin,
    out_origin: MutOrigin,
](ParallelTask):
    """Finds and formats (as `bwa fastmap`) the SMEMs of one slice of query
    records in the FM-index."""

    var fm: ImmPointer[FMIndex, Self.index_origin]
    var query: ImmPointer[FastaRecords, Self.query_origin]
    var bounds: ImmPointer[List[Int], Self.bounds_origin]
    var texts: MutPointer[List[UInt8], Self.out_origin]
    var first: Int
    var min_len: Int
    var max_hits: Int
    var min_occ: Int
    var max_occ: Int

    def __init__(
        out self,
        fm: ImmPointer[FMIndex, Self.index_origin],
        query: ImmPointer[FastaRecords, Self.query_origin],
        bounds: ImmPointer[List[Int], Self.bounds_origin],
        texts: MutPointer[List[UInt8], Self.out_origin],
        first: Int,
        min_len: Int,
        max_hits: Int,
        min_occ: Int,
        max_occ: Int,
    ):
        self.fm = fm
        self.query = query
        self.bounds = bounds
        self.texts = texts
        self.first = first
        self.min_len = min_len
        self.max_hits = max_hits
        self.min_occ = min_occ
        self.max_occ = max_occ

    def run(self, index: Int):
        ref query = self.query[]
        ref fm = self.fm[]
        var lo = self.bounds[][self.first + index]
        var hi = self.bounds[][self.first + index + 1]
        var starts = List[Int](capacity=hi - lo + 1)
        starts.append(0)
        for q in range(lo, hi):
            starts.append(starts[len(starts) - 1] + query.lengths[q])
        var codes = List[UInt8](unsafe_uninit_length=starts[len(starts) - 1] + 32)
        for q in range(lo, hi):
            encode_bases(query.seq(q), codes.unsafe_ptr().unsafe_offset(starts[q - lo]))
        var found = fm.smems_encoded(codes, starts, self.min_len, self.max_hits, self.min_occ, self.max_occ)
        var out = OutputBuffer(fd=-1)
        for q in range(lo, hi):
            out.write("SQ\t")
            out.write(query.names[q])
            out.write_byte(ord("\t"))
            out.write_int(query.lengths[q])
            out.write_byte(ord("\n"))
            for s in range(found.offsets[q - lo], found.offsets[q - lo + 1]):
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
                        out.write_byte(ord("\t"))
                        out.write(fm.names[Int(h.sequence)])
                        out.write(":-" if h.reverse else ":+")
                        out.write_int(Int(h.position) + 1)
                out.write_byte(ord("\n"))
            out.write("//\n")
        self.texts[unsafe_offset=index] = out.take()


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
            for s in range(found.offsets[0], found.offsets[1]):
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
                        out.write_byte(ord("\t"))
                        out.write(self.names[][Int(h.sequence)])
                        out.write(":-" if h.reverse else ":+")
                        out.write_int(Int(h.position) + 1)
                out.write_byte(ord("\n"))
            out.write("//\n")
        self.texts[unsafe_offset=index] = out.take()
