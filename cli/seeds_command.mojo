"""`sporos seeds`: a seed index of the reference, looked up with every read."""

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
            # Per line: the two names, three numbers, a tab, a strand, a tab and a newline.
            var room = at + len(read_hits) * (read_room + ref_room + 3 * 10 + 5) + 16
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
