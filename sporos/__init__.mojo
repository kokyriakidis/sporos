"""sporos: fast, exact DNA seeding in Mojo.

The public API is what this package exports; everything else in its modules
is implementation and may change.

- Sequences: `read_fasta`, `FastaRecords` (FASTA/FASTQ, plain, gzip or BGZF); `PackedReads` (two bits per
  base, other letters as positions); `ByteStream` (a file's inflated bytes a chunk at a time, for read sets too
  large to hold whole).
- Seeds: `minimizer_positions`, `syncmer_positions` (one sequence),
  `minimizers_batch`, `syncmers_batch` (many reads, `SeedBatch`).
- Seed index: `SeedIndex` (`MINIMIZER_SEEDS`, `SYNCMER_SEEDS`, `ALL_KMERS`),
  with `SeedHit`, `SeedHits`.
- Matches against a reference: `Reference` (`mems` -> `Match`, `smems` ->
  `Smems` of `Smem` and `SmemHit`).
- Building blocks: `SeqText`, `MemParams`, `MemIndex`, `find_mems` -> `Mem`,
  `sort_mems`, `filter_mums`, `filter_by_occurrences`, `occurrence_counts`, `find_smems` (MEM
  engine); `FMIndex` (SMEMs of reads, saved and loaded indexes).
- yak's k-mer counting and hifiasm's minimizer seeds, bit for bit as their C code (`sporos.yak`):
  `KmerCounts`, `FilterTable`, `analyze_count`, `sketch_read` -> `Sketch`, `sketch_all`, `PositionTable`.
- `default_threads()`: one thread per physical core.
"""

from .batch import SeedBatch, minimizers_batch, syncmers_batch
from .fasta import FastaRecords, read_fasta
from .packed import PackedReads
from .stream import ByteStream
from .fmindex import FMIndex
from .index import ALL_KMERS, MINIMIZER_SEEDS, SYNCMER_SEEDS, SeedHit, SeedHits, SeedIndex
from .mems import (
    Mem,
    MemIndex,
    MemParams,
    Smem,
    SmemHit,
    Smems,
    filter_by_occurrences,
    occurrence_counts,
    filter_mums,
    find_mems,
    find_smems,
    sort_mems,
)
from .minimizers import minimizer_positions
from .parallel import default_threads
from .reference import Match, Reference
from .syncmers import syncmer_positions
from .text import SeqText
from .yak import FilterTable, KmerCounts, PositionTable, Sketch, analyze_count, sketch_all, sketch_read
