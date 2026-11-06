"""hifiasm's position table (`ha_pt_gen`, htab.cpp): every minimizer occurrence of every read, by k-mer hash.

Two passes over the sketches: the first counts each hash's occurrences (a yak count table), from whose histogram
the coverage peaks are taken again (`ha_pt_gen`'s `peak_hom` is what the anchoring uses as `hom_cov`), keeping
the hashes occurring 2 to 4094 times; the second lays the occurrences out contiguously per hash. An occurrence
(`ha_idxpos_t`) packs the read, the END position, the strand and the span.
"""

from sporos.fasta import FastaRecords
from sporos.parallel import ParallelTask, parallel_for
from sporos.yak.count import EMPTY, PARTITIONS, YAK_COUNTER_BITS, YAK_MAX_COUNT, YAK_N_COUNTS, CountTable, log2_floor
from sporos.yak.filter import FilterTable
from sporos.yak.histogram import Peaks, analyze_count
from sporos.yak.sketch import Minimizer, Sketch, sketch_read


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


struct _SketchTask[records_origin: ImmOrigin, filter_origin: ImmOrigin, out_origin: MutOrigin](ParallelTask):
    """Sketches one read. Tracked origins keep the inputs alive while the task runs."""

    var records: ImmPointer[FastaRecords, Self.records_origin]
    var filter: ImmPointer[FilterTable, Self.filter_origin]
    var w: Int
    var k: Int
    var is_hpc: Bool
    var sample_dist: Int
    var rewin: Int
    var out: MutPointer[Sketch, Self.out_origin]

    def __init__(
        out self,
        records: ImmPointer[FastaRecords, Self.records_origin],
        filter: ImmPointer[FilterTable, Self.filter_origin],
        w: Int,
        k: Int,
        is_hpc: Bool,
        sample_dist: Int,
        rewin: Int,
        out_ptr: MutPointer[Sketch, Self.out_origin],
    ):
        self.records = records
        self.filter = filter
        self.w = w
        self.k = k
        self.is_hpc = is_hpc
        self.sample_dist = sample_dist
        self.rewin = rewin
        self.out = out_ptr

    def run(self, r: Int):
        var sketch = Sketch()
        sketch_read(
            self.records[].seq(r), self.w, self.k, self.is_hpc, self.filter[], self.sample_dist, self.rewin, sketch
        )
        for i in range(len(sketch.mz)):
            sketch.mz[i].rid = UInt32(r)
        self.out[unsafe_offset=r] = sketch^


def sketch_all(
    records: FastaRecords,
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
    var task = _SketchTask(Pointer(to=records), Pointer(to=filter), w, k, is_hpc, sample_dist, rewin, out.unsafe_ptr())
    parallel_for(task, n, threads)
    return out^


struct _Partition(Movable):
    """One of 4096 maps from a hash's high bits to its occurrences: `keys` hold `high << 12 | filled` while the
    occurrences are inserted, `starts` and `counts` where in `positions` a hash's occurrences lie."""

    var keys: List[UInt64]
    var starts: List[UInt32]
    var counts: List[UInt16]
    var positions: List[UInt64]
    var shift: UInt64

    def __init__(out self):
        self.keys = List[UInt64]()
        self.starts = List[UInt32]()
        self.counts = List[UInt16]()
        self.positions = List[UInt64]()
        self.shift = 63

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
        self.keys = List[UInt64](length=capacity, fill=EMPTY)
        self.shift = UInt64(64 - log2_floor(capacity))
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
            while self.keys[i] != EMPTY:
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
        while True:
            var key = self.keys[i]
            if key == EMPTY:
                return -1
            if key >> YAK_COUNTER_BITS == high:
                return i
            i = (i + 1) & (len(self.keys) - 1)

    def insert(mut self, hash: UInt64, position: UInt64):
        """`ha_pt_insert_list`: one occurrence, if the hash is kept and has room."""
        var i = self._find(hash)
        if i < 0:
            return
        var filled = Int(self.keys[i] & YAK_MAX_COUNT)
        if filled >= Int(self.counts[i]):
            return
        self.positions[Int(self.starts[i]) + filled] = position
        self.keys[i] += 1

    def get(self, hash: UInt64) -> Tuple[Int, Int]:
        """Where the hash's occurrences lie in `positions` and how many (0 for an unknown hash)."""
        var i = self._find(hash)
        if i < 0:
            return (0, 0)
        return (Int(self.starts[i]), Int(self.counts[i]))


struct _BucketTask[sketches_origin: ImmOrigin, buckets_origin: MutOrigin](ParallelTask):
    """Sorts one run of reads' minimizers into this thread's per-partition buckets of packed
    (hash, position) pairs."""

    var sketches: ImmPointer[Sketch, Self.sketches_origin]
    var n_reads: Int
    var runs: Int
    var buckets: MutPointer[List[List[UInt64]], Self.buckets_origin]

    def __init__(
        out self,
        sketches: ImmPointer[Sketch, Self.sketches_origin],
        n_reads: Int,
        runs: Int,
        buckets: MutPointer[List[List[UInt64]], Self.buckets_origin],
    ):
        self.sketches = sketches
        self.n_reads = n_reads
        self.runs = runs
        self.buckets = buckets

    def run(self, index: Int):
        var lo = self.n_reads * index // self.runs
        var hi = self.n_reads * (index + 1) // self.runs
        ref mine = self.buckets[unsafe_offset=index]
        for r in range(lo, hi):
            ref sketch = self.sketches[unsafe_offset=r]
            for m in sketch.mz:
                ref bucket = mine[Int(m.x & (PARTITIONS - 1))]
                bucket.append(m.x)
                bucket.append(pack_position(Int(m.rid), Int(m.pos), Int(m.rev), Int(m.span)))


struct _InsertTask[buckets_origin: ImmOrigin, partitions_origin: MutOrigin](ParallelTask):
    """Inserts every thread's bucket of one partition, in read order, into that partition."""

    var buckets: ImmPointer[List[List[UInt64]], Self.buckets_origin]
    var runs: Int
    var partitions: MutPointer[_Partition, Self.partitions_origin]

    def __init__(
        out self,
        buckets: ImmPointer[List[List[UInt64]], Self.buckets_origin],
        runs: Int,
        partitions: MutPointer[_Partition, Self.partitions_origin],
    ):
        self.buckets = buckets
        self.runs = runs
        self.partitions = partitions

    def run(self, p: Int):
        ref part = self.partitions[unsafe_offset=p]
        if len(part.keys) == 0:
            return
        for t in range(self.runs):
            ref bucket = self.buckets[unsafe_offset=t][p]
            var i = 0
            while i < len(bucket):
                part.insert(bucket[i], bucket[i + 1])
                i += 2


struct PositionTable(Movable):
    """`ha_pt_t`: the occurrences of every kept minimizer hash, and the coverage peaks of the minimizer counts."""

    var partitions: List[_Partition]
    var peaks: Peaks
    var collected: Int
    """Minimizers in the sketches."""
    var counted: Int
    """Distinct minimizer hashes among them."""
    var distinct: Int
    """Distinct hashes kept (2 to 4094 occurrences)."""
    var indexed: Int
    """Occurrences laid out."""

    def __init__(out self, sketches: List[Sketch], threads: Int = 1, min_hist_count: Int = 5):
        """Counts the sketches' minimizers, takes the peaks (the histogram's low point searched from
        `min_hist_count`, hifiasm's `min_hist_kmer_cnt`), keeps the hashes occurring 2 to 4094 times and lays
        their occurrences out."""
        var counts = List[CountTable](capacity=PARTITIONS)
        for _ in range(PARTITIONS):
            counts.append(CountTable())
        var collected = 0
        for s in range(len(sketches)):
            for m in sketches[s].mz:
                counts[Int(m.x & (PARTITIONS - 1))].add(m.x)
                collected += 1
        var cnt = List[Int](length=YAK_N_COUNTS, fill=0)
        var total_distinct = 0
        for p in range(PARTITIONS):
            total_distinct += counts[p].size
            for key in counts[p].keys:
                if key != EMPTY:
                    cnt[Int(key & YAK_MAX_COUNT)] += 1
        self.peaks = analyze_count(cnt, min_hist_count)
        self.collected = collected
        self.counted = total_distinct
        self.partitions = List[_Partition](capacity=PARTITIONS)
        self.distinct = 0
        self.indexed = 0
        for p in range(PARTITIONS):
            var part = _Partition()
            part.build(counts[p], 2, YAK_MAX_COUNT - 1)
            self.indexed += len(part.positions)
            self.partitions.append(part^)
        for p in range(PARTITIONS):
            for key in self.partitions[p].keys:
                if key != EMPTY:
                    self.distinct += 1
        var runs = max(1, threads)
        var buckets = List[List[List[UInt64]]](capacity=runs)
        for _ in range(runs):
            var mine = List[List[UInt64]](capacity=PARTITIONS)
            for _ in range(PARTITIONS):
                mine.append(List[UInt64]())
            buckets.append(mine^)
        var bucket_task = _BucketTask(sketches.unsafe_ptr(), len(sketches), runs, buckets.unsafe_ptr())
        parallel_for(bucket_task, runs, threads)
        var task = _InsertTask(buckets.unsafe_ptr().as_imm(), runs, self.partitions.unsafe_ptr())
        parallel_for(task, PARTITIONS, threads)

    @inline(.always)
    def get(self, hash: UInt64) -> Tuple[Int, Int]:
        """`ha_pt_get`: the start of the hash's occurrences in its partition's `positions` and their number."""
        return self.partitions[Int(hash & (PARTITIONS - 1))].get(hash)

    @inline(.always)
    def position(self, hash: UInt64, i: Int) -> UInt64:
        return self.partitions[Int(hash & (PARTITIONS - 1))].positions[i]
