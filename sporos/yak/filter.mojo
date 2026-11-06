"""hifiasm's high-occurrence filter (`gen_hh`, `ha_ft_cnt`, htab.cpp:1038): the k-mers counted `cutoff` times or
more, keyed by their full hash, each with its count, or `INT16_MAX` when over `max_kmer_cnt` (2000)."""

from sporos.parallel import ParallelTask, parallel_for
from sporos.yak.count import EMPTY, PARTITIONS, YAK_COUNTER_BITS, YAK_MAX_COUNT, CountTable, KmerCounts, log2_floor


comptime MAX_KMER_CNT = 2000
"""hifiasm's `max_kmer_cnt` (CommandLines.cpp:169): a filtered k-mer counted more often is never a seed."""
comptime _INT16_MAX = 32767
comptime _INT32_MAX = 2147483647


struct _CollectTask[tables_origin: ImmOrigin, out_origin: MutOrigin](ParallelTask):
    """One group of partitions' k-mers counted at least `cutoff` times, as `(hash, count)`."""

    var tables: ImmPointer[CountTable, Self.tables_origin]
    var groups: Int
    var cutoff: Int
    var out: MutPointer[List[Tuple[UInt64, Int]], Self.out_origin]

    def __init__(
        out self,
        tables: ImmPointer[CountTable, Self.tables_origin],
        groups: Int,
        cutoff: Int,
        out_ptr: MutPointer[List[Tuple[UInt64, Int]], Self.out_origin],
    ):
        self.tables = tables
        self.groups = groups
        self.cutoff = cutoff
        self.out = out_ptr

    def run(self, g: Int):
        var found = List[Tuple[UInt64, Int]]()
        for p in range(PARTITIONS * g // self.groups, PARTITIONS * (g + 1) // self.groups):
            ref table = self.tables[unsafe_offset=p]
            var keys = table.keys.unsafe_ptr()
            for i in range(len(table.keys)):
                var key = keys[unsafe_offset=i]
                if key == EMPTY:
                    continue
                var c = Int(key & YAK_MAX_COUNT)
                if c >= self.cutoff:
                    # The original hash: the table's high bits with the partition as the low 12 bits.
                    found.append(((key >> YAK_COUNTER_BITS << YAK_COUNTER_BITS) | UInt64(p), c))
        self.out[unsafe_offset=g] = found^


struct FilterTable(Movable):
    """hifiasm's high-occurrence filter (`gen_hh`): the k-mers counted `cutoff` times or more, keyed by their full
    hash, each with its count, or `INT16_MAX` when over `max_cnt` (2000). `count` answers as `ha_ft_cnt`: 0 for
    a k-mer not in the table, the count, or `INT32_MAX` for one over the maximum.

    The sketch asks it once per k-mer of every read and nearly always for one it does not hold, so the table is
    kept at most an eighth full, and an occupancy bitmap (one bit a slot, a sixty-fourth of the keys' memory)
    answers most of those questions alone: a k-mer whose home slot is empty is not in the table."""

    var keys: List[UInt64]
    var values: List[Int16]
    var used: List[UInt64]
    """One bit per slot: set when the slot holds a k-mer."""
    var size: Int
    var shift: UInt64
    var mask: Int

    def __init__(out self):
        """An empty filter: every k-mer counts 0, none is filtered."""
        self.keys = List[UInt64](length=64, fill=0)
        self.values = List[Int16](length=64, fill=0)
        self.used = List[UInt64](length=1, fill=0)
        self.size = 0
        self.shift = UInt64(64 - log2_floor(64))
        self.mask = 63

    def __init__(out self, counts: KmerCounts, cutoff: Int, max_cnt: Int = MAX_KMER_CNT, threads: Int = 1):
        """The k-mers of `counts` counted at least `cutoff` times; the count tables are scanned in parallel."""
        var limit = min(max_cnt, YAK_MAX_COUNT - 1, _INT16_MAX - 1)
        var groups = max(1, min(PARTITIONS, 4 * threads))
        var found = List[List[Tuple[UInt64, Int]]](capacity=groups)
        for _ in range(groups):
            found.append(List[Tuple[UInt64, Int]]())
        var task = _CollectTask(counts.tables.unsafe_ptr(), groups, cutoff, found.unsafe_ptr())
        parallel_for(task, groups, threads)
        var n = 0
        for g in range(groups):
            n += len(found[g])
        var capacity = 64
        while capacity < 8 * n:
            capacity *= 2
        self.keys = List[UInt64](length=capacity, fill=0)
        self.values = List[Int16](length=capacity, fill=0)
        self.used = List[UInt64](length=capacity // 64, fill=0)
        self.size = 0
        self.shift = UInt64(64 - log2_floor(capacity))
        self.mask = capacity - 1
        for g in range(groups):
            for entry in found[g]:
                var c = entry[1]
                self._put(entry[0], Int16(_INT16_MAX) if c > limit else Int16(c))

    @inline(.always)
    def _slot(self, y: UInt64) -> Int:
        return Int((y * 0x9E3779B97F4A7C15) >> self.shift)

    @inline(.always)
    def _occupied(self, i: Int) -> Bool:
        return (self.used.unsafe_ptr()[unsafe_offset=i >> 6] >> UInt64(i & 63)) & 1 == 1

    def _put(mut self, y: UInt64, value: Int16):
        var i = self._slot(y)
        while self._occupied(i) and self.keys[i] != y:
            i = (i + 1) & self.mask
        if not self._occupied(i):
            self.size += 1
            self.used[i >> 6] |= UInt64(1) << UInt64(i & 63)
        self.keys[i] = y
        self.values[i] = value

    @inline(.always)
    def count(self, y: UInt64) -> Int:
        """`ha_ft_cnt`."""
        var i = self._slot(y)
        if not self._occupied(i):
            return 0
        var keys = self.keys.unsafe_ptr()
        while self._occupied(i):
            if keys[unsafe_offset=i] == y:
                var v = Int(self.values.unsafe_ptr()[unsafe_offset=i])
                return _INT32_MAX if v == _INT16_MAX else v
            i = (i + 1) & self.mask
        return 0
