"""Saving and loading an FM-index (one file; the large arrays are read in
parallel into huge pages)."""

from sporos.parallel import ParallelTask, parallel_for

from .build import _aligned_first, _huge_cache, _huge_uninit
from .index import FMIndex, _MAGIC, _SAMPLE_BYTES, _VERSION


def save_index(fm: FMIndex, path: String) raises:
    """Writes the index to `path` (read it back with `FMIndex.load`)."""
    with open(path, "w") as f:
        f.write_bytes(_MAGIC.as_bytes())
        var header: List[UInt64] = [
            UInt64(_VERSION),
            UInt64(fm.rows),
            UInt64(fm.half),
            UInt64(len(fm.separators)),
            UInt64(len(fm.run_seq)),
            UInt64(len(fm.names)),
            UInt64(fm.cache_k),
            UInt64(fm.sample_bits),
        ]
        _write_words(f, Span(header))
        _write_words(f, Span(fm._blocks)[fm._first : fm._first + ((fm.rows >> 7) + 1) * 8])
        _write_bytes(f, fm.samples.unsafe_ptr(), len(fm.samples))
        _write_words(f, Span(fm.packed))
        _write_bytes(f, fm.cache.unsafe_ptr().unsafe_bitcast[UInt8](), len(fm.cache) * 16)
        _write_ints(f, fm.separators)
        _write_ints(f, fm.separator_pos)
        _write_ints(f, fm.run_text)
        _write_ints(f, fm.run_seq)
        _write_ints(f, fm.run_ref)
        _write_ints(f, fm.lengths)
        for name in fm.names:
            var length: List[UInt64] = [UInt64(name.byte_length())]
            _write_words(f, Span(length))
            f.write_bytes(name.as_bytes())


def load_index(path: String, threads: Int = 1) raises -> FMIndex:
    """Reads an index written by `save`, with up to `threads` threads."""
    var index = FMIndex()
    var ranges = List[Tuple[Int, Int, Int]]()
    with open(path, "r") as f:
        var magic = List[UInt8](length=8, fill=0)
        _read_exact(f, magic.unsafe_ptr(), 8)
        var expected = _MAGIC.as_bytes()
        for i in range(8):
            if magic[i] != expected[i]:
                raise Error(t"{path} is not a sporos FM-index")
        var header = List[UInt64](length=8, fill=0)
        _read_words(f, header)
        if header[0] != UInt64(_VERSION):
            raise Error(t"{path}: unsupported FM-index version {header[0]}")
        index.rows = Int(header[1])
        index.half = Int(header[2])
        var nsep = Int(header[3])
        var nruns = Int(header[4])
        var nseqs = Int(header[5])
        # The large arrays are read in parallel below; skip them here.
        var nblocks = (index.rows >> 7) + 1
        index._blocks = _huge_uninit[DType.uint64](nblocks * 8 + 8)
        index._first = _aligned_first(index._blocks)
        index.sample_bits = Int(header[7])
        index.samples = _huge_uninit[DType.uint8](index._sample_count() * _SAMPLE_BYTES + 8)
        index.packed = _huge_uninit[DType.uint64]((index.half >> 5) + 2)
        index.cache_k = Int(header[6])
        if index.cache_k > 0:
            index.cache = _huge_cache(index.cache_k)
        var at = 8 + 8 * len(header)
        ranges.append((at, Int(index._blocks.unsafe_ptr().unsafe_offset(index._first)), nblocks * 64))
        at += nblocks * 64
        ranges.append((at, Int(index.samples.unsafe_ptr()), len(index.samples)))
        at += len(index.samples)
        ranges.append((at, Int(index.packed.unsafe_ptr()), len(index.packed) * 8))
        at += len(index.packed) * 8
        ranges.append((at, Int(index.cache.unsafe_ptr()), len(index.cache) * 16))
        at += len(index.cache) * 16
        _ = f.seek(at)
        index.separators = _read_ints(f, nsep)
        index.separator_pos = _read_ints(f, nsep)
        index.run_text = _read_ints(f, nruns + 1)
        index.run_seq = _read_ints(f, nruns)
        index.run_ref = _read_ints(f, nruns)
        index.lengths = _read_ints(f, nseqs)
        var offset = 0
        for s in range(nseqs):
            index.offsets.append(offset)
            offset += index.lengths[s] + 1
            var length = List[UInt64](length=1, fill=0)
            _read_words(f, length)
            var bytes = List[UInt8](length=Int(length[0]), fill=0)
            _read_exact(f, bytes.unsafe_ptr(), len(bytes))
            index.names.append(String(unsafe_from_utf8=Span(bytes)))
    var failed = List[Int](length=1, fill=0)
    for r in ranges:
        var slices = max(1, min(threads, r[2] >> 24))
        var task = _RangeReadTask(path, r[0], r[1], r[2], slices, failed.unsafe_ptr())
        parallel_for(task, slices, threads)
    if failed[0] != 0:
        raise Error(t"{path}: truncated FM-index file")
    _ = len(index.cache)  # keep every array alive while tasks write through addresses
    index.counts = index.rank4(index.rows)
    var first = Int64(len(index.separators))
    index.starts = SIMD[DType.int64, 4](
        first,
        first + index.counts[0],
        first + index.counts[0] + index.counts[1],
        first + index.counts[0] + index.counts[1] + index.counts[2],
    )
    return index^


comptime _IO_CHUNK = 1 << 30


def _write_words(mut f: FileHandle, words: Span[UInt64, _]) raises:
    var bytes = Span[UInt8, ImmUntrackedOrigin](
        unsafe_ptr=words.unsafe_ptr().unsafe_bitcast[UInt8]().as_imm().unsafe_origin_cast[ImmUntrackedOrigin](),
        length=len(words) * 8,
    )
    for at in range(0, len(bytes), _IO_CHUNK):
        f.write_bytes(bytes[at : min(len(bytes), at + _IO_CHUNK)])


def _write_bytes(mut f: FileHandle, data: ImmPointer[UInt8, _], length: Int) raises:
    for at in range(0, length, _IO_CHUNK):
        f.write_bytes(
            Span[UInt8, ImmUntrackedOrigin](
                unsafe_ptr=data.unsafe_offset(at).unsafe_origin_cast[ImmUntrackedOrigin](),
                length=min(_IO_CHUNK, length - at),
            )
        )


def _write_ints(mut f: FileHandle, values: List[Int]) raises:
    var words = List[UInt64](capacity=len(values))
    for v in values:
        words.append(UInt64(v))
    _write_words(f, words)


def _read_exact(mut f: FileHandle, dest: MutPointer[UInt8, _], length: Int) raises:
    var done = 0
    while done < length:
        var room = min(length - done, _IO_CHUNK)
        var got = Int(
            f.read(
                Span[UInt8, MutUntrackedOrigin](
                    unsafe_ptr=dest.unsafe_offset(done).unsafe_origin_cast[MutUntrackedOrigin](), length=room
                )
            )
        )
        if got <= 0:
            raise Error("truncated FM-index file")
        done += got


def _read_words(mut f: FileHandle, mut words: List[UInt64]) raises:
    _read_exact(f, words.unsafe_ptr().unsafe_bitcast[UInt8](), len(words) * 8)


struct _RangeReadTask(ParallelTask):
    """Reads slice `index` of `length` bytes at file offset `offset` to
    address `dest`."""

    var path: String
    var offset: Int
    var dest: Int
    var length: Int
    var slices: Int
    var failed: MutPointer[Int, MutUntrackedOrigin]

    def __init__(out self, path: String, offset: Int, dest: Int, length: Int, slices: Int, failed: MutPointer[Int, _]):
        self.path = path.copy()
        self.offset = offset
        self.dest = dest
        self.length = length
        self.slices = slices
        self.failed = failed.unsafe_origin_cast[MutUntrackedOrigin]()

    def run(self, index: Int):
        var per = (self.length + self.slices - 1) // self.slices
        var lo = index * per
        var hi = min(self.length, lo + per)
        if lo >= hi:
            return
        try:
            with open(self.path, "r") as f:
                _ = f.seek(self.offset + lo)
                var dest = MutPointer[UInt8, MutUntrackedOrigin](unsafe_from_address=self.dest)
                _read_exact(f, dest.unsafe_offset(lo), hi - lo)
        except:
            self.failed[] = 1


def _read_ints(mut f: FileHandle, count: Int) raises -> List[Int]:
    var words = List[UInt64](length=count, fill=0)
    _read_words(f, words)
    var values = List[Int](capacity=count)
    for w in words:
        values.append(Int(w))
    return values^
