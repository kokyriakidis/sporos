"""Counting of homopolymer-compressed k-mers (hifiasm's `ha_ft_gen` through `ha_count`, htab.cpp).

Each table key is `hash >> 12 << 12 | count` as in yak, and the hash's low 12 bits choose one of 4096 tables,
which are filled in parallel: a block of reads is scanned on every thread, each thread sorting its k-mers into
per-partition buffers, then every partition inserts its buffers, through its Bloom filter when there is one
(`-f`, see `bloom.mojo`); `-f0` counts exactly."""

from std.bit import count_leading_zeros

from sporos.fasta import FastaRecords
from sporos.parallel import ParallelTask, parallel_for

from sporos.yak.bloom import BF_N_HASHES, DEFAULT_BF_SHIFT, BloomFilters
from sporos.yak.hash import HASH_LANES, Kmer, hpc_run_starts, nt4, yak_hash64_64, yak_hash64_64_lanes


comptime YAK_COUNTER_BITS = 12
"""yak keeps a k-mer's count in the low 12 bits of its table key."""
comptime YAK_MAX_COUNT = (1 << YAK_COUNTER_BITS) - 1
"""Counts saturate at 4095."""
comptime YAK_N_COUNTS = 1 << YAK_COUNTER_BITS
comptime PARTITIONS = 1 << YAK_COUNTER_BITS
"""The count tables are split 4096 ways by a hash's low 12 bits (`pre = YAK_COUNTER_BITS`, htab.cpp:62)."""
comptime EMPTY: UInt64 = 0
"""An unused table slot."""


struct CountTable(Movable):
    """An open-addressing set of yak keys (`hash >> 12 << 12 | count`), linear probing, load at most a half."""

    var keys: List[UInt64]
    var size: Int
    var shift: UInt64
    """`64 - log2(capacity)`: a slot is the top bits of the mixed key."""

    def __init__(out self, capacity: Int = 16):
        self.keys = List[UInt64](length=capacity, fill=EMPTY)
        self.size = 0
        self.shift = UInt64(64 - log2_floor(capacity))

    @inline(.always)
    def _slot(self, high: UInt64) -> Int:
        # The high bits are a hash already (yak_ct_hash is the key >> 12); fold them into the capacity.
        return Int((high * 0x9E3779B97F4A7C15) >> self.shift)

    def _grow(mut self):
        var old = self.keys^
        self.keys = List[UInt64](length=len(old) * 2, fill=EMPTY)
        self.shift = UInt64(64 - log2_floor(len(self.keys)))
        for key in old:
            if key != EMPTY:
                var i = self._slot(key >> YAK_COUNTER_BITS)
                while self.keys[i] != EMPTY:
                    i = (i + 1) & (len(self.keys) - 1)
                self.keys[i] = key

    def reserve(mut self, entries: Int):
        """Grows the table, if needed, to hold `entries` at most half full."""
        var capacity = len(self.keys)
        while 2 * entries > capacity:
            capacity *= 2
        if capacity == len(self.keys):
            return
        var old = self.keys^
        self.keys = List[UInt64](length=capacity, fill=EMPTY)
        self.shift = UInt64(64 - log2_floor(capacity))
        for key in old:
            if key != EMPTY:
                var i = self._slot(key >> YAK_COUNTER_BITS)
                while self.keys[i] != EMPTY:
                    i = (i + 1) & (capacity - 1)
                self.keys[i] = key

    def add(mut self, hash: UInt64):
        """Counts one occurrence of the k-mer with this hash (saturating at 4095)."""
        if 2 * (self.size + 1) > len(self.keys):
            self._grow()
        var high = hash >> YAK_COUNTER_BITS
        var i = self._slot(high)
        while True:
            var key = self.keys[i]
            if key == EMPTY:
                self.keys[i] = (high << YAK_COUNTER_BITS) | 1
                self.size += 1
                return
            if key >> YAK_COUNTER_BITS == high:
                if key & YAK_MAX_COUNT < YAK_MAX_COUNT:
                    self.keys[i] = key + 1
                return
            i = (i + 1) & (len(self.keys) - 1)

    def add_seen_twice(mut self, hash: UInt64, bloom: BloomFilters, filter: Int):
        """`ha_ct_insert_list` with a Bloom filter: the k-mer is counted only once the filter has seen it, and
        enters the table at 2 then, its first sighting made up for."""
        var high = hash >> YAK_COUNTER_BITS
        if bloom.insert(filter, high) != BF_N_HASHES:
            return
        if 2 * (self.size + 1) > len(self.keys):
            self._grow()
        var i = self._slot(high)
        while True:
            var key = self.keys[i]
            if key == EMPTY:
                self.keys[i] = (high << YAK_COUNTER_BITS) | 2
                self.size += 1
                return
            if key >> YAK_COUNTER_BITS == high:
                if key & YAK_MAX_COUNT < YAK_MAX_COUNT:
                    self.keys[i] = key + 1
                return
            i = (i + 1) & (len(self.keys) - 1)

    def count(self, hash: UInt64) -> Int:
        """The count of this k-mer, 0 if absent."""
        var high = hash >> YAK_COUNTER_BITS
        var i = self._slot(high)
        while True:
            var key = self.keys[i]
            if key == EMPTY:
                return 0
            if key >> YAK_COUNTER_BITS == high:
                return Int(key & YAK_MAX_COUNT)
            i = (i + 1) & (len(self.keys) - 1)


@inline(.always)
def log2_floor(n: Int) -> Int:
    """`floor(log2(n))` for `n >= 1`."""
    return 63 - Int(count_leading_zeros(UInt64(n)))


def hpc_kmer_hashes(seq: Span[UInt8, _], k: Int, mut out: List[UInt64]):
    """`count_seq_buf_HPC`: the canonical hash of every homopolymer-compressed k-mer of `seq` (symmetric k-mers
    included), appended to `out`."""
    var scratch = HpcScratch()
    hpc_kmer_hashes(seq, k, out, scratch)


struct HpcScratch(Movable):
    """What hashing a read's compressed k-mers needs besides its output, kept by a caller hashing many reads."""

    var codes: List[UInt8]
    var starts: List[UInt32]
    var planes: List[UInt64]

    def __init__(out self):
        self.codes = List[UInt8]()
        self.starts = List[UInt32]()
        self.planes = List[UInt64]()


@no_inline
def hpc_kmer_hashes(seq: Span[UInt8, _], k: Int, mut out: List[UInt64], mut scratch: HpcScratch):
    """`hpc_kmer_hashes` through `scratch`. The read's homopolymer runs are found a vector at a time
    (`hpc_run_starts`), the rolling k-mer steps once per run, and the k-mers' canonical planes are hashed
    `HASH_LANES` at a time."""
    var runs = hpc_run_starts(seq, scratch.codes, scratch.starts)
    var start = len(out)
    out.resize(unsafe_uninit_length=start + runs)
    var need = 2 * runs + 2 * HASH_LANES
    if len(scratch.planes) < need:
        scratch.planes.resize(unsafe_uninit_length=need)
    var pp = scratch.planes.unsafe_ptr()
    var codes = scratch.codes.unsafe_ptr().unsafe_offset(1)
    var starts = scratch.starts.unsafe_ptr()
    var half = runs + HASH_LANES
    var kmer = Kmer.empty(k)
    var l = 0
    var n = 0
    for r in range(runs):
        var c = codes[unsafe_offset=Int(starts[unsafe_offset=r])]
        if c < 4:
            kmer.push(c)
            l += 1
            if l >= k:
                var canonical = kmer.canonical()
                pp[unsafe_offset=n] = canonical[0]
                pp[unsafe_offset=half + n] = canonical[1]
                n += 1
        else:
            l = 0
            kmer.clear()
    var p0 = pp
    var p1 = pp.unsafe_offset(half)
    var dst = out.unsafe_ptr().unsafe_offset(start)
    var i = 0
    while i + HASH_LANES <= n:
        dst.unsafe_store(
            i,
            yak_hash64_64_lanes(p0.unsafe_load[width=HASH_LANES](i))
            + yak_hash64_64_lanes(p1.unsafe_load[width=HASH_LANES](i)),
        )
        i += HASH_LANES
    while i < n:
        dst[unsafe_offset=i] = yak_hash64_64(p0[unsafe_offset=i]) + yak_hash64_64(p1[unsafe_offset=i])
        i += 1
    out.resize(unsafe_uninit_length=start + n)


struct _Block(Movable):
    """One thread's share of a block of reads: the hashes of its k-mers grouped by partition, each group in read
    order (`hashes[offsets[p] : offsets[p + 1]]` for partition `p`)."""

    var hashes: List[UInt64]
    var offsets: List[Int]

    def __init__(out self):
        self.hashes = List[UInt64]()
        self.offsets = List[Int]()


struct _ScanTask[records_origin: ImmOrigin, blocks_origin: MutOrigin](ParallelTask):
    """Hashes the k-mers of one run of reads, then groups them by partition with a counting sort, which keeps
    each partition's hashes in read order. Tracked origins keep the inputs alive while the task runs."""

    var records: ImmPointer[FastaRecords, Self.records_origin]
    var first: Int
    var last: Int
    var runs: Int
    var k: Int
    var blocks: MutPointer[_Block, Self.blocks_origin]

    def __init__(
        out self,
        records: ImmPointer[FastaRecords, Self.records_origin],
        first: Int,
        last: Int,
        runs: Int,
        k: Int,
        blocks: MutPointer[_Block, Self.blocks_origin],
    ):
        self.records = records
        self.first = first
        self.last = last
        self.runs = runs
        self.k = k
        self.blocks = blocks

    def _share(self, index: Int) -> Tuple[Int, Int]:
        """Run `index`'s reads: the contiguous range whose bases are its equal share of the block's (a range
        starts at the first read starting at or past its share), so the runs take the same time however the
        read lengths vary, and the runs in order are the block's reads in order."""
        ref lengths = self.records[].lengths
        var total = 0
        for r in range(self.first, self.last):
            total += lengths[r]
        var lo = self.last
        var hi = self.last
        var before = 0
        var r = self.first
        while r < self.last:
            if before * self.runs >= total * index:
                lo = r
                break
            before += lengths[r]
            r += 1
        while r < self.last:
            if before * self.runs >= total * (index + 1):
                hi = r
                break
            before += lengths[r]
            r += 1
        if index == self.runs - 1:
            hi = self.last
        return (lo, hi)

    def run(self, index: Int):
        var range_ = self._share(index)
        var lo = range_[0]
        var hi = range_[1]
        # A read has at most as many k-mers as letters: room for all of them up front, as `List.resize` grows
        # to exactly the length asked for and would copy the list once per read.
        var bases = 0
        for r in range(lo, hi):
            bases += self.records[].lengths[r]
        var hashes = List[UInt64](capacity=bases)
        var scratch = HpcScratch()
        # The partition sizes are counted read by read, while each read's hashes are still in cache.
        var offsets = List[Int](length=PARTITIONS + 1, fill=0)
        var op = offsets.unsafe_ptr()
        for r in range(lo, hi):
            var before = len(hashes)
            hpc_kmer_hashes(self.records[].seq(r), self.k, hashes, scratch)
            var fresh = hashes.unsafe_ptr()
            for i in range(before, len(hashes)):
                op[unsafe_offset=Int(fresh[unsafe_offset=i] & (PARTITIONS - 1)) + 1] += 1
        var total = len(hashes)
        var src = hashes.unsafe_ptr()
        for p in range(PARTITIONS):
            op[unsafe_offset=p + 1] += op[unsafe_offset=p]
        var grouped = List[UInt64](unsafe_uninit_length=total)
        var dst = grouped.unsafe_ptr()
        var at = offsets.copy()
        var ap = at.unsafe_ptr()
        for i in range(total):
            var h = src[unsafe_offset=i]
            var p = Int(h & (PARTITIONS - 1))
            dst[unsafe_offset=ap[unsafe_offset=p]] = h
            ap[unsafe_offset=p] += 1
        ref block = self.blocks[unsafe_offset=index]
        block.hashes = grouped^
        block.offsets = offsets^


comptime _AHEAD = 16
"""How many inserts ahead the insert loop prefetches."""


struct _InsertTask[tables_origin: MutOrigin, bloom_origin: ImmOrigin, blocks_origin: ImmOrigin](ParallelTask):
    """Inserts every thread's hashes of one partition into that partition's table, in read order, through the
    partition's Bloom filter when there is one."""

    var tables: MutPointer[CountTable, Self.tables_origin]
    var blooms: ImmPointer[BloomFilters, Self.bloom_origin]
    var with_bloom: Bool
    var blocks: ImmPointer[_Block, Self.blocks_origin]
    var runs: Int

    def __init__(
        out self,
        tables: MutPointer[CountTable, Self.tables_origin],
        blooms: ImmPointer[BloomFilters, Self.bloom_origin],
        with_bloom: Bool,
        blocks: ImmPointer[_Block, Self.blocks_origin],
        runs: Int,
    ):
        self.tables = tables
        self.blooms = blooms
        self.with_bloom = with_bloom
        self.blocks = blocks
        self.runs = runs

    def run(self, partition: Int):
        ref table = self.tables[unsafe_offset=partition]
        var incoming = 0
        for t in range(self.runs):
            ref block = self.blocks[unsafe_offset=t]
            incoming += block.offsets[partition + 1] - block.offsets[partition]
        if not self.with_bloom:
            # Room for a quarter of the incoming k-mers being new, at most half full; growth covers the rest.
            table.reserve(table.size + incoming // 4)
        for t in range(self.runs):
            ref block = self.blocks[unsafe_offset=t]
            var hashes = block.hashes.unsafe_ptr()
            var end = block.offsets[partition + 1]
            # Each Bloom insert is a cache miss in the partition's 4 MB filter: the block `_AHEAD` inserts on
            # is started now. (The table itself stays in L2, where prefetching it only costs.)
            if self.with_bloom:
                ref bloom = self.blooms[]
                for i in range(block.offsets[partition], end):
                    if i + _AHEAD < end:
                        bloom.prefetch(partition, hashes[unsafe_offset=i + _AHEAD] >> YAK_COUNTER_BITS)
                    table.add_seen_twice(hashes[unsafe_offset=i], bloom, partition)
            else:
                for i in range(block.offsets[partition], end):
                    table.add(hashes[unsafe_offset=i])


struct KmerCounts(Movable):
    """The counts of every compressed k-mer of a read set, in 4096 partitioned tables."""

    var tables: List[CountTable]
    var k: Int
    var total_bases: Int

    def __init__(
        out self,
        records: FastaRecords,
        k: Int,
        threads: Int = 1,
        block: Int = 1 << 24,
        bf_shift: Int = DEFAULT_BF_SHIFT,
    ) raises:
        """Counts the k-mers of `records`, `block` bases of reads at a time, through Bloom filters of
        `2^bf_shift` bits in all (0: none, exact counts)."""
        self.k = k
        self.total_bases = records.total_length()
        self.tables = List[CountTable](capacity=PARTITIONS)
        for _ in range(PARTITIONS):
            self.tables.append(CountTable())
        # `ha_ct_init`: a filter per partition only when it is bigger than the partition index itself.
        var with_bloom = bf_shift > YAK_COUNTER_BITS
        var blooms = BloomFilters(PARTITIONS if with_bloom else 1, bf_shift - YAK_COUNTER_BITS if with_bloom else 12)
        var runs = max(1, threads)
        var n = len(records)
        var first = 0
        while first < n:
            var bases = 0
            var last = first
            while last < n and bases < block:
                bases += records.lengths[last]
                last += 1
            var blocks = List[_Block](capacity=runs)
            for _ in range(runs):
                blocks.append(_Block())
            var scan = _ScanTask(Pointer(to=records), first, last, runs, k, blocks.unsafe_ptr())
            parallel_for(scan, runs, threads)
            var insert = _InsertTask(
                self.tables.unsafe_ptr(), Pointer(to=blooms).as_imm(), with_bloom, blocks.unsafe_ptr().as_imm(), runs
            )
            parallel_for(insert, PARTITIONS, threads)
            first = last

    def count(self, hash: UInt64) -> Int:
        return self.tables[Int(hash & (PARTITIONS - 1))].count(hash)

    def distinct(self) -> Int:
        var total = 0
        for p in range(len(self.tables)):
            total += self.tables[p].size
        return total

    def histogram(self, threads: Int = 1) -> List[Int]:
        """`ha_ct_hist`: how many distinct k-mers have each count, 0 .. 4095; the partitions in groups, one
        histogram per group, in parallel."""
        return partition_histogram(self.tables, threads)


struct _HistogramTask[tables_origin: ImmOrigin, out_origin: MutOrigin](ParallelTask):
    """One group of partitions' histogram."""

    var tables: ImmPointer[CountTable, Self.tables_origin]
    var groups: Int
    var out: MutPointer[List[Int], Self.out_origin]

    def __init__(
        out self,
        tables: ImmPointer[CountTable, Self.tables_origin],
        groups: Int,
        out_ptr: MutPointer[List[Int], Self.out_origin],
    ):
        self.tables = tables
        self.groups = groups
        self.out = out_ptr

    def run(self, g: Int):
        var cnt = List[Int](length=YAK_N_COUNTS, fill=0)
        for p in range(PARTITIONS * g // self.groups, PARTITIONS * (g + 1) // self.groups):
            ref table = self.tables[unsafe_offset=p]
            var keys = table.keys.unsafe_ptr()
            for i in range(len(table.keys)):
                var key = keys[unsafe_offset=i]
                if key != EMPTY:
                    cnt[Int(key & YAK_MAX_COUNT)] += 1
        self.out[unsafe_offset=g] = cnt^


def partition_histogram(tables: List[CountTable], threads: Int) -> List[Int]:
    """The count histogram of `PARTITIONS` count tables, summed from per-group histograms made in parallel."""
    var groups = max(1, min(PARTITIONS, 4 * threads))
    var parts = List[List[Int]](capacity=groups)
    for _ in range(groups):
        parts.append(List[Int]())
    var task = _HistogramTask(tables.unsafe_ptr(), groups, parts.unsafe_ptr())
    parallel_for(task, groups, threads)
    var cnt = List[Int](length=YAK_N_COUNTS, fill=0)
    for g in range(groups):
        for c in range(YAK_N_COUNTS):
            cnt[c] += parts[g][c]
    return cnt^
