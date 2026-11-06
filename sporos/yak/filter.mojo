"""hifiasm's high-occurrence filter (`gen_hh`, `ha_ft_cnt`, htab.cpp:1038): the k-mers counted `cutoff` times or
more, keyed by their full hash, each with its count, or `INT16_MAX` when over `max_kmer_cnt` (2000)."""

from sporos.yak.count import EMPTY, YAK_COUNTER_BITS, YAK_MAX_COUNT, KmerCounts, log2_floor


comptime MAX_KMER_CNT = 2000
"""hifiasm's `max_kmer_cnt` (CommandLines.cpp:169): a filtered k-mer counted more often is never a seed."""
comptime _INT16_MAX = 32767
comptime _INT32_MAX = 2147483647


struct FilterTable(Movable):
    """hifiasm's high-occurrence filter (`gen_hh`): the k-mers counted `cutoff` times or more, keyed by their full
    hash, each with its count, or `INT16_MAX` when over `max_cnt` (2000). `count` answers as `ha_ft_cnt`: 0 for
    a k-mer not in the table, the count, or `INT32_MAX` for one over the maximum."""

    var keys: List[UInt64]
    var values: List[Int16]
    var size: Int
    var shift: UInt64

    def __init__(out self):
        """An empty filter: every k-mer counts 0, none is filtered."""
        self.keys = List[UInt64](length=16, fill=EMPTY)
        self.values = List[Int16](length=16, fill=0)
        self.size = 0
        self.shift = UInt64(64 - log2_floor(16))

    def __init__(out self, counts: KmerCounts, cutoff: Int, max_cnt: Int = MAX_KMER_CNT):
        var limit = min(max_cnt, YAK_MAX_COUNT - 1, _INT16_MAX - 1)
        var n = 0
        for p in range(len(counts.tables)):
            for key in counts.tables[p].keys:
                if key != EMPTY and Int(key & YAK_MAX_COUNT) >= cutoff:
                    n += 1
        var capacity = 16
        while capacity < 2 * n + 1:
            capacity *= 2
        self.keys = List[UInt64](length=capacity, fill=EMPTY)
        self.values = List[Int16](length=capacity, fill=0)
        self.size = 0
        self.shift = UInt64(64 - log2_floor(capacity))
        for p in range(len(counts.tables)):
            for key in counts.tables[p].keys:
                if key == EMPTY:
                    continue
                var c = Int(key & YAK_MAX_COUNT)
                if c < cutoff:
                    continue
                # The original hash: the table's high bits with the partition as the low 12 bits.
                var y = (key >> YAK_COUNTER_BITS << YAK_COUNTER_BITS) | UInt64(p)
                self._put(y, Int16(_INT16_MAX) if c > limit else Int16(c))

    @inline(.always)
    def _slot(self, y: UInt64) -> Int:
        return Int((y * 0x9E3779B97F4A7C15) >> self.shift)

    def _put(mut self, y: UInt64, value: Int16):
        var i = self._slot(y)
        while self.keys[i] != EMPTY and self.keys[i] != y:
            i = (i + 1) & (len(self.keys) - 1)
        if self.keys[i] == EMPTY:
            self.size += 1
        self.keys[i] = y
        self.values[i] = value

    def count(self, y: UInt64) -> Int:
        """`ha_ft_cnt`."""
        if self.size == 0:
            return 0
        var i = self._slot(y)
        while True:
            var key = self.keys[i]
            if key == EMPTY:
                return 0
            if key == y:
                var v = Int(self.values[i])
                return _INT32_MAX if v == _INT16_MAX else v
            i = (i + 1) & (len(self.keys) - 1)
