"""hifiasm's position table (`ha_pt_gen`, htab.cpp): every minimizer occurrence of every read, by k-mer hash.

Built as hifiasm builds it, without ever holding the reads' sketches: two passes over the reads, a batch of
reads at a time, each sketching the batch afresh. The first counts each hash's occurrences (a yak count table),
from whose histogram the coverage peaks are taken again (`ha_pt_gen`'s `peak_hom` is what the anchoring uses as
`hom_cov`), keeping the hashes occurring 2 to 4094 times, whose occurrences are then given their places; the
second sketches the reads again and puts each occurrence in its place, in read order. An occurrence
(`ha_idxpos_t`) packs the read, the END position, the strand and the span.
"""

from std.atomic import Atomic

from sporos.packed import PackedReads
from sporos.parallel import ParallelTask, parallel_for
from sporos.yak.count import (
    EMPTY,
    PARTITIONS,
    YAK_COUNTER_BITS,
    YAK_MAX_COUNT,
    CountTable,
    log2_floor,
    partition_histogram,
    share_by_bases,
)
from sporos.yak.filter import FilterTable
from sporos.yak.histogram import Peaks, analyze_count
from sporos.yak.sketch import Minimizer, Sketch, SketchScratch, sketch_read


@inline(.always)
def pack_position(rid: Int, pos: Int, rev: Int, span: Int) -> UInt64:
    """`ha_idxpos_t`: rid in 28 bits, pos in 27, rev in 1, span in 8."""
    return (UInt64(rid) << 36) | (UInt64(rev) << 35) | (UInt64(pos) << 8) | UInt64(span)


@inline(.always)
def position_rid(p: UInt64) -> Int:
    return Int(p >> 36)


@inline(.always)
def position_rev(p: UInt64) -> Int:
    return Int((p >> 35) & 1)


@inline(.always)
def position_pos(p: UInt64) -> Int:
    return Int((p >> 8) & 0x7FFFFFF)


@inline(.always)
def position_span(p: UInt64) -> Int:
    return Int(p & 0xFF)


struct _SketchTask[
    records_origin: ImmOrigin,
    filter_origin: ImmOrigin,
    order_origin: ImmOrigin,
    next_origin: MutOrigin,
    out_origin: MutOrigin,
](ParallelTask):
    """One worker's share of the sketching: it takes reads off a shared counter, longest first so that no long
    read is left for the end, and sketches them in one working memory. Tracked origins keep the inputs alive
    while the task runs."""

    var records: ImmPointer[PackedReads, Self.records_origin]
    var filter: ImmPointer[FilterTable, Self.filter_origin]
    var order: ImmPointer[Int, Self.order_origin]
    var next: MutPointer[Atomic[Int64], Self.next_origin]
    var w: Int
    var k: Int
    var is_hpc: Bool
    var sample_dist: Int
    var rewin: Int
    var out: MutPointer[Sketch, Self.out_origin]

    def __init__(
        out self,
        records: ImmPointer[PackedReads, Self.records_origin],
        filter: ImmPointer[FilterTable, Self.filter_origin],
        order: ImmPointer[Int, Self.order_origin],
        next: MutPointer[Atomic[Int64], Self.next_origin],
        w: Int,
        k: Int,
        is_hpc: Bool,
        sample_dist: Int,
        rewin: Int,
        out_ptr: MutPointer[Sketch, Self.out_origin],
    ):
        self.records = records
        self.filter = filter
        self.order = order
        self.next = next
        self.w = w
        self.k = k
        self.is_hpc = is_hpc
        self.sample_dist = sample_dist
        self.rewin = rewin
        self.out = out_ptr

    def run(self, worker: Int):
        var n = len(self.records[])
        var scratch = SketchScratch()
        while True:
            var i = Int(self.next[].fetch_add(1))
            if i >= n:
                break
            var r = self.order[unsafe_offset=i]
            ref sketch = self.out[unsafe_offset=r]
            sketch_read(
                self.records[],
                r,
                self.w,
                self.k,
                self.is_hpc,
                self.filter[],
                self.sample_dist,
                self.rewin,
                sketch,
                scratch,
            )
            for j in range(len(sketch.mz)):
                sketch.mz[j].rid = UInt32(r)


def _longest_first(records: PackedReads) -> List[Int]:
    """The reads' indices by length, longest first."""
    var n = len(records)
    var keys = List[UInt64](capacity=n)
    for r in range(n):
        keys.append((UInt64(records.lengths[r]) << 32) | UInt64(n - 1 - r))
    sort(keys)
    var order = List[Int](capacity=n)
    for i in range(n - 1, -1, -1):
        order.append(n - 1 - Int(keys[i] & 0xFFFFFFFF))
    return order^


def sketch_all(
    records: PackedReads,
    filter: FilterTable,
    w: Int,
    k: Int,
    threads: Int,
    sample_dist: Int,
    rewin: Int,
    is_hpc: Bool = True,
) -> List[Sketch]:
    """Every read's filtered, subsampled sketch (`sketch_read`, `rid` set to the read's index), in parallel.
    hifiasm's overlap path uses `sample_dist` 500 and `rewin` 1000."""
    var n = len(records)
    var out = List[Sketch](capacity=n)
    for _ in range(n):
        out.append(Sketch())
    var order = _longest_first(records)
    var next = Atomic[Int64](0)
    var workers = max(1, min(threads, n))
    var task = _SketchTask(
        Pointer(to=records),
        Pointer(to=filter),
        order.unsafe_ptr().as_imm(),
        Pointer(to=next),
        w,
        k,
        is_hpc,
        sample_dist,
        rewin,
        out.unsafe_ptr(),
    )
    parallel_for(task, workers, threads)
    return out^


struct _Partition(Movable):
    """One of 4096 maps from a hash's high bits to its occurrences: `keys` hold `high << 12 | filled` while the
    occurrences are inserted, `starts` and `counts` where in `positions` a hash's occurrences lie. A slot is in
    use when its count is (a kept hash occurs at least twice), so no key value is reserved."""

    var keys: List[UInt64]
    var starts: List[UInt32]
    var counts: List[UInt16]
    var positions: List[UInt64]
    var shift: UInt64
    var kept: Int

    def __init__(out self):
        self.keys = List[UInt64]()
        self.starts = List[UInt32]()
        self.counts = List[UInt16]()
        self.positions = List[UInt64]()
        self.shift = 63
        self.kept = 0

    @inline(.always)
    def _slot(self, high: UInt64) -> Int:
        return Int((high * 0x9E3779B97F4A7C15) >> self.shift)

    def build(mut self, counted: CountTable, low: Int, high_count: Int):
        """Lays out the hashes of `counted` whose count lies in `low ..= high_count`."""
        var n = 0
        var total = 0
        for key in counted.keys:
            if key != EMPTY:
                var c = Int(key & YAK_MAX_COUNT)
                if c >= low and c <= high_count:
                    n += 1
                    total += c
        var capacity = 16
        while capacity < 2 * n + 1:
            capacity *= 2
        self.keys = List[UInt64](length=capacity, fill=0)
        self.shift = UInt64(64 - log2_floor(capacity))
        self.kept = n
        self.starts = List[UInt32](length=capacity, fill=0)
        self.counts = List[UInt16](length=capacity, fill=0)
        self.positions = List[UInt64](length=total, fill=0)
        var offset = 0
        for key in counted.keys:
            if key == EMPTY:
                continue
            var c = Int(key & YAK_MAX_COUNT)
            if c < low or c > high_count:
                continue
            var high = key >> YAK_COUNTER_BITS
            var i = self._slot(high)
            while self.counts[i] != 0:
                i = (i + 1) & (len(self.keys) - 1)
            self.keys[i] = high << YAK_COUNTER_BITS
            self.starts[i] = UInt32(offset)
            self.counts[i] = UInt16(c)
            offset += c

    @inline(.always)
    def _find(self, hash: UInt64) -> Int:
        """The slot of `hash`, or -1."""
        if len(self.keys) == 0:
            return -1
        var high = hash >> YAK_COUNTER_BITS
        var i = self._slot(high)
        while self.counts[i] != 0:
            if self.keys[i] >> YAK_COUNTER_BITS == high:
                return i
            i = (i + 1) & (len(self.keys) - 1)
        return -1

    def insert(mut self, hash: UInt64, position: UInt64) -> Bool:
        """`ha_pt_insert_list`: one occurrence, if the hash is kept and has room; whether it was placed."""
        var i = self._find(hash)
        if i < 0:
            return False
        var filled = Int(self.keys[i] & YAK_MAX_COUNT)
        if filled >= Int(self.counts[i]):
            return False
        self.positions[Int(self.starts[i]) + filled] = position
        self.keys[i] += 1
        return True

    def get(self, hash: UInt64) -> Tuple[Int, Int]:
        """Where the hash's occurrences lie in `positions` and how many (0 for an unknown hash)."""
        var i = self._find(hash)
        if i < 0:
            return (0, 0)
        return (Int(self.starts[i]), Int(self.counts[i]))


comptime DEFAULT_BATCH = 1 << 26
"""Bases of reads the position table sketches at a time (hifiasm's `chunk_size` is 20 Mbp)."""


struct _MinimizerBlock(Movable):
    """One run's minimizers of a batch, grouped by partition, each group in read order: hashes (and packed
    positions, in the second pass) `offsets[p] ..< offsets[p + 1]` for partition `p`."""

    var hashes: List[UInt64]
    var positions: List[UInt64]
    var offsets: List[Int]

    def __init__(out self):
        self.hashes = List[UInt64]()
        self.positions = List[UInt64]()
        self.offsets = List[Int]()


struct _MinimizerScanTask[reads_origin: ImmOrigin, filter_origin: ImmOrigin, blocks_origin: MutOrigin](ParallelTask):
    """Sketches one run of a batch's reads (its equal share of the bases, in order) and groups the minimizers by
    partition with a counting sort, which keeps each partition's in read order."""

    var reads: ImmPointer[PackedReads, Self.reads_origin]
    var filter: ImmPointer[FilterTable, Self.filter_origin]
    var first: Int
    var last: Int
    var runs: Int
    var w: Int
    var k: Int
    var is_hpc: Bool
    var sample_dist: Int
    var rewin: Int
    var with_positions: Bool
    var blocks: MutPointer[_MinimizerBlock, Self.blocks_origin]

    def __init__(
        out self,
        reads: ImmPointer[PackedReads, Self.reads_origin],
        filter: ImmPointer[FilterTable, Self.filter_origin],
        first: Int,
        last: Int,
        runs: Int,
        w: Int,
        k: Int,
        is_hpc: Bool,
        sample_dist: Int,
        rewin: Int,
        with_positions: Bool,
        blocks: MutPointer[_MinimizerBlock, Self.blocks_origin],
    ):
        self.reads = reads
        self.filter = filter
        self.first = first
        self.last = last
        self.runs = runs
        self.w = w
        self.k = k
        self.is_hpc = is_hpc
        self.sample_dist = sample_dist
        self.rewin = rewin
        self.with_positions = with_positions
        self.blocks = blocks

    def run(self, index: Int):
        ref reads = self.reads[]
        var share = share_by_bases(reads.lengths, self.first, self.last, self.runs, index)
        var bases = 0
        for r in range(share[0], share[1]):
            bases += reads.lengths[r]
        # About one minimizer per 27 bases on hifiasm's settings: room for twice that, growth covers the rest.
        var hashes = List[UInt64](capacity=bases // 13 + 16)
        var positions = List[UInt64](capacity=bases // 13 + 16 if self.with_positions else 0)
        var offsets = List[Int](length=PARTITIONS + 1, fill=0)
        var op = offsets.unsafe_ptr()
        var scratch = SketchScratch()
        var sketch = Sketch()
        for r in range(share[0], share[1]):
            sketch_read(
                reads,
                r,
                self.w,
                self.k,
                self.is_hpc,
                self.filter[],
                self.sample_dist,
                self.rewin,
                sketch,
                scratch,
            )
            for m in sketch.mz:
                hashes.append(m.x)
                op[unsafe_offset=Int(m.x & (PARTITIONS - 1)) + 1] += 1
                if self.with_positions:
                    positions.append(pack_position(r, Int(m.pos), Int(m.rev), Int(m.span)))
        for p in range(PARTITIONS):
            op[unsafe_offset=p + 1] += op[unsafe_offset=p]
        var total = len(hashes)
        var at = offsets.copy()
        var ap = at.unsafe_ptr()
        var src = hashes.unsafe_ptr()
        var grouped = List[UInt64](unsafe_uninit_length=total)
        var gp = grouped.unsafe_ptr()
        ref block = self.blocks[unsafe_offset=index]
        if self.with_positions:
            var pos = positions.unsafe_ptr()
            var grouped_positions = List[UInt64](unsafe_uninit_length=total)
            var gq = grouped_positions.unsafe_ptr()
            for i in range(total):
                var h = src[unsafe_offset=i]
                var p = Int(h & (PARTITIONS - 1))
                var to = ap[unsafe_offset=p]
                gp[unsafe_offset=to] = h
                gq[unsafe_offset=to] = pos[unsafe_offset=i]
                ap[unsafe_offset=p] = to + 1
            block.positions = grouped_positions^
        else:
            for i in range(total):
                var h = src[unsafe_offset=i]
                var p = Int(h & (PARTITIONS - 1))
                gp[unsafe_offset=ap[unsafe_offset=p]] = h
                ap[unsafe_offset=p] += 1
        block.hashes = grouped^
        block.offsets = offsets^


struct _CountTask[blocks_origin: ImmOrigin, counts_origin: MutOrigin](ParallelTask):
    """Counts one partition's minimizer hashes of a batch, from every run's block."""

    var blocks: ImmPointer[_MinimizerBlock, Self.blocks_origin]
    var runs: Int
    var counts: MutPointer[CountTable, Self.counts_origin]

    def __init__(
        out self,
        blocks: ImmPointer[_MinimizerBlock, Self.blocks_origin],
        runs: Int,
        counts: MutPointer[CountTable, Self.counts_origin],
    ):
        self.blocks = blocks
        self.runs = runs
        self.counts = counts

    def run(self, p: Int):
        ref table = self.counts[unsafe_offset=p]
        for t in range(self.runs):
            ref block = self.blocks[unsafe_offset=t]
            var hashes = block.hashes.unsafe_ptr()
            for i in range(block.offsets[p], block.offsets[p + 1]):
                table.add(hashes[unsafe_offset=i])


struct _LayoutTask[counts_origin: ImmOrigin, partitions_origin: MutOrigin](ParallelTask):
    """Lays out one partition's kept hashes: each one's place for its occurrences."""

    var counts: ImmPointer[CountTable, Self.counts_origin]
    var partitions: MutPointer[_Partition, Self.partitions_origin]

    def __init__(
        out self,
        counts: ImmPointer[CountTable, Self.counts_origin],
        partitions: MutPointer[_Partition, Self.partitions_origin],
    ):
        self.counts = counts
        self.partitions = partitions

    def run(self, p: Int):
        self.partitions[unsafe_offset=p].build(self.counts[unsafe_offset=p], 2, YAK_MAX_COUNT - 1)


struct _FillTask[blocks_origin: ImmOrigin, partitions_origin: MutOrigin, placed_origin: MutOrigin](ParallelTask):
    """Puts one partition's occurrences of a batch in their places, from every run's block in order."""

    var blocks: ImmPointer[_MinimizerBlock, Self.blocks_origin]
    var runs: Int
    var partitions: MutPointer[_Partition, Self.partitions_origin]
    var placed: MutPointer[Int, Self.placed_origin]

    def __init__(
        out self,
        blocks: ImmPointer[_MinimizerBlock, Self.blocks_origin],
        runs: Int,
        partitions: MutPointer[_Partition, Self.partitions_origin],
        placed: MutPointer[Int, Self.placed_origin],
    ):
        self.blocks = blocks
        self.runs = runs
        self.partitions = partitions
        self.placed = placed

    def run(self, p: Int):
        ref part = self.partitions[unsafe_offset=p]
        if part.kept == 0:
            return
        var placed = 0
        for t in range(self.runs):
            ref block = self.blocks[unsafe_offset=t]
            var hashes = block.hashes.unsafe_ptr()
            var positions = block.positions.unsafe_ptr()
            for i in range(block.offsets[p], block.offsets[p + 1]):
                placed += Int(part.insert(hashes[unsafe_offset=i], positions[unsafe_offset=i]))
        self.placed[unsafe_offset=p] += placed


struct PositionTable(Movable):
    """`ha_pt_t`: the occurrences of every kept minimizer hash, and the coverage peaks of the minimizer counts."""

    var partitions: List[_Partition]
    var peaks: Peaks
    var collected: Int
    """Minimizers of all the reads."""
    var counted: Int
    """Distinct minimizer hashes among them."""
    var distinct: Int
    """Distinct hashes kept (2 to 4094 occurrences)."""
    var indexed: Int
    """Occurrences laid out."""

    def __init__(
        out self,
        reads: PackedReads,
        filter: FilterTable,
        w: Int,
        k: Int,
        threads: Int,
        sample_dist: Int,
        rewin: Int,
        min_hist_count: Int = 5,
        is_hpc: Bool = True,
        batch: Int = DEFAULT_BATCH,
    ) raises:
        """The table of the reads' sketches (`sketch_read` with these settings, `rid` the read's index): the
        sketches' minimizers counted, the peaks taken (the histogram's low point searched from `min_hist_count`,
        hifiasm's `min_hist_kmer_cnt`), the hashes occurring 2 to 4094 times kept and their occurrences laid out.
        The reads are sketched twice, `batch` bases of them at a time, and no sketch outlives its batch. Each
        batch is sketched on `threads` threads, each its equal share of the bases, and each partition counted or
        filled on its own."""
        var runs = max(1, threads)
        var counts = List[CountTable](capacity=PARTITIONS)
        for _ in range(PARTITIONS):
            counts.append(CountTable())
        self.collected = 0
        var n = len(reads)
        var first = 0
        while first < n:
            var last = _batch_end(reads, first, batch)
            var blocks = _blocks(runs)
            var scan = _MinimizerScanTask(
                Pointer(to=reads),
                Pointer(to=filter),
                first,
                last,
                runs,
                w,
                k,
                is_hpc,
                sample_dist,
                rewin,
                False,
                blocks.unsafe_ptr(),
            )
            parallel_for(scan, runs, threads)
            for t in range(runs):
                self.collected += len(blocks[t].hashes)
            var count = _CountTask(blocks.unsafe_ptr().as_imm(), runs, counts.unsafe_ptr())
            parallel_for(count, PARTITIONS, threads)
            first = last
        self.peaks = analyze_count(partition_histogram(counts, threads), min_hist_count)
        self.counted = 0
        for p in range(PARTITIONS):
            self.counted += counts[p].size
        self.partitions = List[_Partition](capacity=PARTITIONS)
        for _ in range(PARTITIONS):
            self.partitions.append(_Partition())
        var layout = _LayoutTask(counts.unsafe_ptr().as_imm(), self.partitions.unsafe_ptr())
        parallel_for(layout, PARTITIONS, threads)
        # The counts end here (their last use), before the occurrences are placed.
        self.distinct = 0
        self.indexed = 0
        for p in range(PARTITIONS):
            self.distinct += self.partitions[p].kept
            self.indexed += len(self.partitions[p].positions)
        var placed = List[Int](length=PARTITIONS, fill=0)
        first = 0
        while first < n:
            var last = _batch_end(reads, first, batch)
            var blocks = _blocks(runs)
            var scan = _MinimizerScanTask(
                Pointer(to=reads),
                Pointer(to=filter),
                first,
                last,
                runs,
                w,
                k,
                is_hpc,
                sample_dist,
                rewin,
                True,
                blocks.unsafe_ptr(),
            )
            parallel_for(scan, runs, threads)
            var fill = _FillTask(blocks.unsafe_ptr().as_imm(), runs, self.partitions.unsafe_ptr(), placed.unsafe_ptr())
            parallel_for(fill, PARTITIONS, threads)
            first = last
        # hifiasm's `assert(tot_cnt == pt->tot_pos)`: the second pass places exactly the occurrences counted.
        var total_placed = 0
        for p in range(PARTITIONS):
            total_placed += placed[p]
        if total_placed != self.indexed:
            raise Error(t"position table: {total_placed} occurrences placed, {self.indexed} counted")

    @inline(.always)
    def get(self, hash: UInt64) -> Tuple[Int, Int]:
        """`ha_pt_get`: the start of the hash's occurrences in its partition's `positions` and their number."""
        return self.partitions[Int(hash & (PARTITIONS - 1))].get(hash)

    @inline(.always)
    def position(self, hash: UInt64, i: Int) -> UInt64:
        return self.partitions[Int(hash & (PARTITIONS - 1))].positions[i]


def _batch_end(reads: PackedReads, first: Int, batch: Int) -> Int:
    """The end of the batch of reads starting at `first`: whole reads until `batch` bases are reached."""
    var n = len(reads)
    var bases = 0
    var last = first
    while last < n and bases < batch:
        bases += reads.lengths[last]
        last += 1
    return last


def _blocks(runs: Int) -> List[_MinimizerBlock]:
    var blocks = List[_MinimizerBlock](capacity=runs)
    for _ in range(runs):
        blocks.append(_MinimizerBlock())
    return blocks^
