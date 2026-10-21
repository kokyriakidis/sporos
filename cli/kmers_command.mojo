"""`sporos minimizers` and `sporos syncmers`: seed positions per record."""

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
            threads = _threads(args, i)
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
