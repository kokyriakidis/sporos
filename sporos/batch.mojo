"""Seeds of many short sequences (reads), one read per SIMD lane.

Splitting a short read across lanes wastes most of each lane on warm-up
(k + w - 2 bases), so batches give every lane its own read instead: groups of
`LANES` reads are copied into one `N`-padded scratch buffer (reused across
groups) and scanned together. Results are flat, with per-read offsets.
"""

from std.math import ceildiv

from std.memory import unsafe_memcpy

from sporos.minimizers import minimizer_positions
from sporos.syncmers import syncmer_positions
from sporos.parallel import ParallelTask, parallel_for
from sporos.window import (
    CLOSED_SYNCMERS,
    LANES,
    MINIMIZERS,
    OPEN_SYNCMERS,
    PADDING,
    ScanParams,
    scan_reads,
)

comptime _GROUPS_PER_TASK = 64
comptime _LONG = 1 << 13
"""Reads at least this long are scanned alone, split across all lanes."""


struct SeedBatch(Movable, Sized):
    """Seed positions of a batch of reads: read `i`'s are
    `positions[offsets[i] : offsets[i + 1]]`, relative to the read."""

    var positions: List[UInt32]
    var offsets: List[Int]

    def __init__(out self, var positions: List[UInt32], var offsets: List[Int]):
        self.positions = positions^
        self.offsets = offsets^

    def __len__(self) -> Int:
        return len(self.offsets) - 1

    def seeds(ref self, i: Int) -> Span[UInt32, origin_of(self.positions)]:
        return Span(self.positions)[self.offsets[i] : self.offsets[i + 1]]


@fieldwise_init
struct _Config(ImplicitlyCopyable, TrivialRegisterPassable):
    var mode: Int
    var canonical: Bool
    var params: ScanParams


def _windows(cfg: _Config, length: Int) -> Int:
    """Lane windows of a read (minimizer reads with fewer than w k-mers have
    one partial window, handled outside the lanes: -1)."""
    var kmers = length - cfg.params.k + 1 if cfg.mode == MINIMIZERS else length - (cfg.params.span) + 1
    if kmers <= 0:
        return 0
    if cfg.mode == MINIMIZERS:
        return kmers - cfg.params.w + 1 if kmers >= cfg.params.w else -1
    return kmers


struct _BatchTask[mut: Bool, //, origin: Origin[mut=mut]](ParallelTask):
    var reads: ImmPointer[Span[UInt8, Self.origin], ImmUntrackedOrigin]
    var count: Int
    var cfg: _Config
    var results: MutPointer[List[UInt32], MutUntrackedOrigin]
    var counts: MutPointer[Int, MutUntrackedOrigin]

    def __init__(
        out self,
        reads: ImmPointer[Span[UInt8, Self.origin], ImmUntrackedOrigin],
        count: Int,
        cfg: _Config,
        results: MutPointer[List[UInt32], MutUntrackedOrigin],
        counts: MutPointer[Int, MutUntrackedOrigin],
    ):
        self.reads = reads
        self.count = count
        self.cfg = cfg
        self.results = results
        self.counts = counts

    def run(self, index: Int):
        var first = index * _GROUPS_PER_TASK * LANES
        var last = min(self.count, first + _GROUPS_PER_TASK * LANES)
        ref out = self.results[unsafe_offset=index]
        var scratch = List[UInt8]()
        var lane_outs = List[List[UInt32]](length=LANES, fill=List[UInt32]())
        var group = List[Int](capacity=LANES)
        for r in range(first, last):
            if len(self.reads[unsafe_offset=r]) >= _LONG:
                # Long reads scan faster split across all lanes.
                self._group(group, scratch, lane_outs, out)
                group.clear()
                var positions = self._single(r)
                out.extend(Span(positions))
                self.counts[unsafe_offset=r] = len(positions)
                continue
            group.append(r)
            if len(group) == LANES:
                self._group(group, scratch, lane_outs, out)
                group.clear()
        self._group(group, scratch, lane_outs, out)

    def _single(self, r: Int) -> List[UInt32]:
        var read = self.reads[unsafe_offset=r]
        var params = self.cfg.params
        try:
            if self.cfg.mode == MINIMIZERS:
                return minimizer_positions(read, params.k, params.w, self.cfg.canonical)
            var opened: Optional[Int] = None
            if self.cfg.mode == OPEN_SYNCMERS:
                opened = params.offset
            return syncmer_positions(read, params.span, params.k, self.cfg.canonical, opened)
        except:
            return List[UInt32]()

    def _group(
        self,
        group: List[Int],
        mut scratch: List[UInt8],
        mut lane_outs: List[List[UInt32]],
        mut out: List[UInt32],
    ):
        """Scans up to `LANES` reads, one per lane, through a padded copy."""
        var n = len(group)
        if n == 0:
            return
        var max_len = 0
        for l in range(n):
            max_len = max(max_len, len(self.reads[unsafe_offset=group[l]]))
        var stride = max_len + PADDING
        if len(scratch) < LANES * stride:
            scratch.resize(LANES * stride, UInt8(ord("N")))
        var windows = SIMD[DType.uint32, LANES](0)
        var chunk = 0
        for l in range(LANES):
            var dst = scratch.unsafe_ptr().unsafe_offset(l * stride)
            var length = 0
            if l < n:
                var read = self.reads[unsafe_offset=group[l]]
                length = len(read)
                unsafe_memcpy(dest=dst, src=read.unsafe_ptr(), count=length)
                var lw = _windows(self.cfg, length)
                windows[l] = UInt32(max(0, lw))
                chunk = max(chunk, lw)
            # N padding after the read: only up to what the lanes can read.
            for i in range(length, min(stride, length + PADDING + chunk + 16)):
                dst[unsafe_offset=i] = UInt8(ord("N"))
            lane_outs[l].clear()
        if chunk > 0:
            scan_reads(self.cfg.mode, self.cfg.canonical, scratch.unsafe_ptr(), stride, windows, chunk, self.cfg.params, lane_outs)
        for l in range(n):
            var r = group[l]
            if _windows(self.cfg, len(self.reads[unsafe_offset=r])) < 0:
                # Fewer than w k-mers: one partial window, scalar.
                lane_outs[l] = self._single(r)
            out.extend(Span(lane_outs[l]))
            self.counts[unsafe_offset=r] = len(lane_outs[l])


def _batch[
    mut: Bool, //, origin: Origin[mut=mut]
](reads: List[Span[UInt8, origin]], cfg: _Config, threads: Int) -> SeedBatch:
    var groups = ceildiv(len(reads), LANES)
    var tasks = max(1, ceildiv(groups, _GROUPS_PER_TASK))
    var results = List[List[UInt32]](length=tasks, fill=List[UInt32]())
    var counts = List[Int](length=len(reads), fill=0)
    var task = _BatchTask[origin](
        reads.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
        len(reads),
        cfg,
        results.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        counts.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
    )
    parallel_for(task, tasks, threads)
    var offsets = List[Int](capacity=len(reads) + 1)
    offsets.append(0)
    for c in counts:
        offsets.append(offsets[len(offsets) - 1] + c)
    var positions = List[UInt32](capacity=offsets[len(offsets) - 1])
    for r in results:
        positions.extend(Span(r))
    return SeedBatch(positions^, offsets^)


def minimizers_batch[
    mut: Bool, //, origin: Origin[mut=mut]
](reads: List[Span[UInt8, origin]], k: Int, w: Int, canonical: Bool = True, threads: Int = 1) raises -> SeedBatch:
    """The (w, k)-minimizer positions of every read (as `minimizer_positions`)."""
    if k < 1 or k > 64 or w < 1:
        raise Error(t"need 1 <= k <= 64 and w >= 1, got k={k} w={w}")
    return _batch(reads, _Config(MINIMIZERS, canonical, ScanParams(k=k, w=w, span=k, offset=0)), threads)


def syncmers_batch[
    mut: Bool, //, origin: Origin[mut=mut]
](
    reads: List[Span[UInt8, origin]],
    k: Int,
    s: Int,
    canonical: Bool = True,
    open_offset: Optional[Int] = None,
    threads: Int = 1,
) raises -> SeedBatch:
    """The syncmer positions of every read (as `syncmer_positions`)."""
    if s < 1 or s > k:
        raise Error(t"need 1 <= s <= k, got k={k} s={s}")
    var w = k - s + 1
    if case .Some(t) = open_offset:
        if t < 0 or t > k - s:
            raise Error(t"open syncmer offset must be in [0, {k - s}], got {t}")
        return _batch(reads, _Config(OPEN_SYNCMERS, canonical, ScanParams(k=s, w=w, span=k, offset=t)), threads)
    return _batch(reads, _Config(CLOSED_SYNCMERS, canonical, ScanParams(k=s, w=w, span=k, offset=0)), threads)
