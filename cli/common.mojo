"""Shared pieces of the command line: usage, options, logging, slicing."""

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


comptime USAGE = """usage:
  sporos seeds -k K [-w W | --syncmers -s S [--open T[,T...]] | --all-kmers] [--query-kmers]
               [--max-occ N | --keep F] [-t threads] [-v] <reference.fa> <reads.fa>
  sporos mem [-mum | -mumreference | -maxmatch] [-l L] [--max-occ N] [--min-occ M] [-b | -r] [-F] [-n]
             [-t threads] [-v] [-K anchor_k] [-k1 ref_step] [-k2 query_step] <reference.fa> <query.fa>
  sporos smem [-l L] [-w W] [-s min_occ] [--max-occ N] [--engine fm|mem] [-t threads] [-v]
              <reference.fa | index.sfm> <query.fa>
  sporos index [-s sample_bits] [-t threads] [-v] <reference.fa> [index.sfm]
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


def _count_option(args: List[String], i: Int, name: String) raises -> Int:
    """An occurrence count after `name`: at least 1."""
    var value = Int(_value(args, i))
    if value < 1:
        raise Error(t"{name} takes a count of at least 1, not {value}")
    return value


def _threads(args: List[String], i: Int) raises -> Int:
    """The thread count after `-t`: at least 1."""
    var threads = Int(_value(args, i))
    if threads < 1:
        raise Error(t"-t takes a thread count of at least 1, not {threads}")
    return threads


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


comptime _FORMAT_SLICE = 1 << 16
"""Positions per formatting task."""


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
