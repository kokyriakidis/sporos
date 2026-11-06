"""yak's blocked Bloom filter (`yak_bf_t`, htab.cpp:79), through which hifiasm's first k-mer counting pass goes:
`2^bf_shift` bits in all (`-f`, 37 by default, 16 GB), split over the 4096 count partitions, so a k-mer enters a
count table on its second sighting only and the error k-mers seen once never take table memory."""


from std.ffi import external_call
from std.sys.info import CompilationTarget
from std.sys.intrinsics import PrefetchOptions, prefetch


comptime BLK_SHIFT = 9
"""`YAK_BLK_SHIFT`: a Bloom block is 512 bits, one cache line."""
comptime BLK_MASK = (1 << BLK_SHIFT) - 1
comptime BF_N_HASHES = 4
"""`opt->bf_n_hash = 4` (htab.cpp:59)."""
comptime DEFAULT_BF_SHIFT = 37
"""`asm_opt->bf_shift = 37` (CommandLines.cpp:168): 2^37 bits in all, 2^25 per partition."""


struct BloomFilters(Movable):
    """`count` of yak's blocked Bloom filters (`yak_bf_t`, htab.cpp:79), `2^n_shift` bits each in 512-bit blocks:
    a k-mer picks a block by its low bits and `BF_N_HASHES` bit positions in it by two more hashes.

    All of them lie in one anonymous mapping: its pages are zero until first written, so nothing is cleared up
    front and the threads that insert fault them in, in parallel; on Linux the mapping is 2 MB aligned and asks
    for transparent huge pages, as the inserts touch it at random. Filter `f` is the bytes
    `f << (n_shift - 3) ..` of it."""

    var base: Int
    var raw: Int
    var raw_bytes: Int
    var n_shift: Int
    var count: Int

    def __init__(out self, count: Int, n_shift: Int) raises:
        self.count = count
        self.n_shift = n_shift
        var bytes = count << (n_shift - 3)
        comptime HUGE = 1 << 21
        self.raw_bytes = bytes + HUGE
        comptime PROT_READ_WRITE: Int32 = 3
        var flags: Int32
        comptime if CompilationTarget.is_linux():
            flags = 0x22
        else:
            flags = 0x1002
        self.raw = external_call["mmap", Int](0, self.raw_bytes, PROT_READ_WRITE, flags, Int32(-1), 0)
        if self.raw == -1:
            raise Error(t"cannot map {bytes} bytes for the Bloom filters")
        self.base = (self.raw + HUGE - 1) & ~(HUGE - 1)
        # Huge pages: the 16 GB of a full-size filter fault in 2 MB at a time. With 4 KB pages, xaros took about
        # five times longer over the eight million faults.
        comptime if CompilationTarget.is_linux():
            comptime MADV_HUGEPAGE: Int32 = 14
            _ = external_call["madvise", Int32](self.base, bytes, MADV_HUGEPAGE)

    def __deinit__(deinit self):
        _ = external_call["munmap", Int32](self.raw, self.raw_bytes)

    @inline(.always)
    @inline(.always)
    def prefetch(self, filter: Int, hash: UInt64):
        """Starts loading the block `insert(filter, hash)` will touch."""
        var x = UInt64(self.n_shift - BLK_SHIFT)
        var y = Int(hash & ((UInt64(1) << x) - 1))
        var p = MutPointer[UInt8, MutUntrackedOrigin](
            unsafe_from_address=self.base + (filter << (self.n_shift - 3)) + (y << (BLK_SHIFT - 3))
        )
        prefetch[PrefetchOptions().for_write().high_locality()](p)

    def insert(self, filter: Int, hash: UInt64) -> Int:
        """`yak_bf_insert` into filter `filter`: sets the k-mer's bits and returns how many were already set (all
        `BF_N_HASHES` when the k-mer was seen before)."""
        var x = UInt64(self.n_shift - BLK_SHIFT)
        var y = Int(hash & ((UInt64(1) << x) - 1))
        var h1 = Int((hash >> x) & UInt64(BLK_MASK))
        var h2 = Int((hash >> UInt64(self.n_shift)) & UInt64(BLK_MASK))
        if (h2 & 31) == 0:
            h2 = (h2 + 1) & BLK_MASK
        var p = MutPointer[UInt8, MutUntrackedOrigin](
            unsafe_from_address=self.base + (filter << (self.n_shift - 3)) + (y << (BLK_SHIFT - 3))
        )
        var z = h1
        var cnt = 0
        for _ in range(BF_N_HASHES):
            var at = z >> 3
            var u = UInt8(1) << UInt8(z & 7)
            var q = p[unsafe_offset=at]
            if q & u != 0:
                cnt += 1
            p[unsafe_offset=at] = q | u
            z = (z + h2) & BLK_MASK
        return cnt
