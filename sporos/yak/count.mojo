"""Counting of homopolymer-compressed k-mers (hifiasm's `ha_ft_gen` through `ha_count`, htab.cpp).

Each table key is `hash >> 12 << 12 | count` as in yak, and the hash's low 12 bits choose one of 4096 tables,
which are filled in parallel: a block of reads is scanned on every thread, each thread sorting its k-mers into
per-partition buffers, then every partition inserts its buffers, through its Bloom filter when there is one
(`-f`, see `bloom.mojo`); `-f0` counts exactly."""

from std.bit import count_leading_zeros

from sporos.fasta import FastaRecords
from sporos.parallel import ParallelTask, parallel_for

from sporos.yak.bloom import BF_N_HASHES, DEFAULT_BF_SHIFT, BloomFilters
from sporos.yak.hash import HASH_LANES, Kmer, nt4, yak_hash64_64, yak_hash64_64_lanes


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
    included), appended to `out`. The k-mers' canonical planes are collected first, then hashed
    `HASH_LANES` at a time."""
    var length = len(seq)
    var start = len(out)
    out.resize(start + length, 0)
    var planes = List[UInt64](length=2 * length + 2 * HASH_LANES, fill=0)
    var pp = planes.unsafe_ptr()
    var kmer = Kmer.empty(k)
    var l = 0
    var n = 0
    var last: UInt8 = 255
    var src = seq.unsafe_ptr()
    for i in range(length):
        var c = nt4(src[unsafe_offset=i])
        if c < 4:
            if c != last:
                kmer.push(c)
                l += 1
                if l >= k:
                    var canonical = kmer.canonical()
                    pp[unsafe_offset=n] = canonical[0]
                    pp[unsafe_offset=length + HASH_LANES + n] = canonical[1]
                    n += 1
                last = c
        else:
            l = 0
            last = 255
            kmer.clear()
    var p0 = pp
    var p1 = pp.unsafe_offset(length + HASH_LANES)
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
    out.resize(start + n, 0)


struct _ScanTask[records_origin: ImmOrigin, buffers_origin: MutOrigin](ParallelTask):
    """Hashes the k-mers of one run of reads into this thread's partition buffers. Tracked origins keep the
    inputs alive while the task runs."""

    var records: ImmPointer[FastaRecords, Self.records_origin]
    var first: Int
    var last: Int
    var runs: Int
    var k: Int
    var buffers: MutPointer[List[List[UInt64]], Self.buffers_origin]

    def __init__(
        out self,
        records: ImmPointer[FastaRecords, Self.records_origin],
        first: Int,
        last: Int,
        runs: Int,
        k: Int,
        buffers: MutPointer[List[List[UInt64]], Self.buffers_origin],
    ):
        self.records = records
        self.first = first
        self.last = last
        self.runs = runs
        self.k = k
        self.buffers = buffers

    def run(self, index: Int):
        var n = self.last - self.first
        var lo = self.first + n * index // self.runs
        var hi = self.first + n * (index + 1) // self.runs
        var hashes = List[UInt64]()
        ref mine = self.buffers[unsafe_offset=index]
        for r in range(lo, hi):
            hashes.clear()
            hpc_kmer_hashes(self.records[].seq(r), self.k, hashes)
            for h in hashes:
                mine[Int(h & (PARTITIONS - 1))].append(h)


struct _InsertTask[tables_origin: MutOrigin, bloom_origin: ImmOrigin, buffers_origin: ImmOrigin](ParallelTask):
    """Inserts every thread's buffer of one partition into that partition's table, through the partition's
    Bloom filter when there is one."""

    var tables: MutPointer[CountTable, Self.tables_origin]
    var blooms: ImmPointer[BloomFilters, Self.bloom_origin]
    var with_bloom: Bool
    var buffers: ImmPointer[List[List[UInt64]], Self.buffers_origin]
    var runs: Int

    def __init__(
        out self,
        tables: MutPointer[CountTable, Self.tables_origin],
        blooms: ImmPointer[BloomFilters, Self.bloom_origin],
        with_bloom: Bool,
        buffers: ImmPointer[List[List[UInt64]], Self.buffers_origin],
        runs: Int,
    ):
        self.tables = tables
        self.blooms = blooms
        self.with_bloom = with_bloom
        self.buffers = buffers
        self.runs = runs

    def run(self, partition: Int):
        ref table = self.tables[unsafe_offset=partition]
        if self.with_bloom:
            for t in range(self.runs):
                for h in self.buffers[unsafe_offset=t][partition]:
                    table.add_seen_twice(h, self.blooms[], partition)
            return
        for t in range(self.runs):
            for h in self.buffers[unsafe_offset=t][partition]:
                table.add(h)


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
            var buffers = List[List[List[UInt64]]](capacity=runs)
            for _ in range(runs):
                var mine = List[List[UInt64]](capacity=PARTITIONS)
                for _ in range(PARTITIONS):
                    mine.append(List[UInt64]())
                buffers.append(mine^)
            var scan = _ScanTask(Pointer(to=records), first, last, runs, k, buffers.unsafe_ptr())
            parallel_for(scan, runs, threads)
            var insert = _InsertTask(
                self.tables.unsafe_ptr(), Pointer(to=blooms).as_imm(), with_bloom, buffers.unsafe_ptr().as_imm(), runs
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

    def histogram(self) -> List[Int]:
        """`ha_ct_hist`: how many distinct k-mers have each count, 0 .. 4095."""
        var cnt = List[Int](length=YAK_N_COUNTS, fill=0)
        for p in range(len(self.tables)):
            for key in self.tables[p].keys:
                if key != EMPTY:
                    cnt[Int(key & YAK_MAX_COUNT)] += 1
        return cnt^
