"""Fast FASTA/FASTQ reading: whole file in (gzip transparently), records out.

Gzipped input is detected by its magic bytes. BGZF files (bgzip, the common
case for genomes) are inflated block-parallel with libdeflate; other gzip
files with libdeflate in one shot; without libdeflate, with the system zlib.
Both libraries are loaded at run time, so neither is a build dependency.
"""

from std.bit import count_trailing_zeros
from std.ffi import OwnedDLHandle
from std.math import ceildiv
from std.os import getenv, stat

from sporos.parallel import ParallelTask, parallel_for
from std.memory import pack_bits, unsafe_memcpy, unsafe_memmove
from std.sys.info import CompilationTarget

comptime _HEADER = UInt8(ord(">"))
comptime _AT = UInt8(ord("@"))
comptime _NEWLINE = UInt8(ord("\n"))
comptime _RETURN = UInt8(ord("\r"))
comptime _SPACE = UInt8(ord(" "))
comptime _TAB = UInt8(ord("\t"))


struct FastaRecords(Movable, Sized):
    """FASTA records: names (first header token) and raw sequence bytes."""

    var names: List[String]
    var data: List[UInt8]
    var starts: List[Int]
    var lengths: List[Int]

    def __init__(out self):
        self.names = List[String]()
        self.data = List[UInt8]()
        self.starts = List[Int]()
        self.lengths = List[Int]()

    def __len__(self) -> Int:
        return len(self.names)

    def seq(ref self, i: Int) -> Span[UInt8, origin_of(self.data)]:
        return Span(self.data)[self.starts[i] : self.starts[i] + self.lengths[i]]

    def total_length(self) -> Int:
        return len(self.data)


@inline(.always)
def _line_end(src: ImmPointer[UInt8, _], start: Int, n: Int) -> Int:
    """Index of the first newline at or after `start`, or `n`."""
    comptime W = 32
    var i = start
    while i + W <= n:
        var hits = src.unsafe_load[width=W](i).eq(_NEWLINE)
        if hits.reduce_or():
            return i + Int(count_trailing_zeros(pack_bits[DType.uint32](hits)))
        i += W
    while i < n and src[unsafe_offset=i] != _NEWLINE:
        i += 1
    return i


def parse_fastq(bytes: Span[UInt8, _]) raises -> FastaRecords:
    """Parses four-line FASTQ records (qualities are skipped)."""
    var records = FastaRecords()
    var n = len(bytes)
    var src = bytes.unsafe_ptr()
    var i = 0
    while i < n:
        var c = src[unsafe_offset=i]
        if c == _NEWLINE or c == _RETURN:
            i += 1
            continue
        if c != _AT:
            raise Error("FASTQ record must start with '@'")
        var j = i + 1
        while (
            j < n
            and src[unsafe_offset=j] != _SPACE
            and src[unsafe_offset=j] != _TAB
            and src[unsafe_offset=j] != _NEWLINE
            and src[unsafe_offset=j] != _RETURN
        ):
            j += 1
        records.names.append(String(StringSlice(unsafe_from_utf8=bytes[i + 1 : j])))
        # A header at the end of the input has an empty sequence.
        var seq_start = min(_line_end(src, j, n) + 1, n)
        var seq_end = min(_line_end(src, seq_start, n), n)
        var trimmed = seq_end
        if trimmed > seq_start and src[unsafe_offset=trimmed - 1] == _RETURN:
            trimmed -= 1
        records.starts.append(len(records.data))
        records.lengths.append(trimmed - seq_start)
        records.data.extend(bytes[seq_start:trimmed])
        var plus_end = _line_end(src, seq_end + 1, n)
        i = _line_end(src, plus_end + 1, n) + 1
    return records^


def parse_sequences(bytes: Span[UInt8, _]) raises -> FastaRecords:
    """Parses FASTA or FASTQ, telling them apart by the first record."""
    for i in range(len(bytes)):
        var c = bytes[i]
        if c == _AT:
            return parse_fastq(bytes)
        if c != _NEWLINE and c != _RETURN and c != _SPACE:
            break
    return parse_fasta(bytes)


def parse_fasta(bytes: Span[UInt8, _]) raises -> FastaRecords:
    """Parses FASTA text; line breaks (LF or CRLF) inside sequences are dropped."""
    var records = FastaRecords()
    records.data.reserve(len(bytes))
    var n = len(bytes)
    var src = bytes.unsafe_ptr()
    var i = 0
    while i < n:
        var c = src[unsafe_offset=i]
        if c == _HEADER:
            var name_start = i + 1
            var j = name_start
            while (
                j < n
                and src[unsafe_offset=j] != _SPACE
                and src[unsafe_offset=j] != _TAB
                and src[unsafe_offset=j] != _NEWLINE
                and src[unsafe_offset=j] != _RETURN
            ):
                j += 1
            records.names.append(String(StringSlice(unsafe_from_utf8=bytes[name_start:j])))
            records.starts.append(len(records.data))
            records.lengths.append(0)
            i = _line_end(src, j, n) + 1
            continue
        if len(records.names) == 0:
            if c == _NEWLINE or c == _RETURN or c == _SPACE:
                i += 1
                continue
            raise Error("FASTA input must start with a '>' header line")
        # A sequence line: copy up to the line break in one block.
        var j = _line_end(src, i, n)
        var end = j
        if end > i and src[unsafe_offset=end - 1] == _RETURN:
            end -= 1
        var length = end - i
        if length > 0:
            var old = len(records.data)
            records.data.resize(unsafe_uninit_length=old + length)
            unsafe_memcpy(
                dest=records.data.unsafe_ptr().unsafe_offset(old),
                src=src.unsafe_offset(i),
                count=length,
            )
            records.lengths[len(records.lengths) - 1] += length
        i = j + 1
    return records^


def read_fasta(path: String, threads: Int = 1) raises -> FastaRecords:
    """Reads a FASTA or FASTQ file, gzipped or not (BGZF inflates on
    `threads` threads)."""
    var bytes = read_file(path, threads)
    for i in range(len(bytes)):
        var c = bytes[i]
        if c == _AT:
            return parse_fastq(Span(bytes))
        if c != _NEWLINE and c != _RETURN and c != _SPACE:
            break
    return parse_fasta_in_place(bytes^)


def parse_fasta_in_place(var bytes: List[UInt8]) raises -> FastaRecords:
    """Parses FASTA, compacting the sequence bytes to the front of `bytes`
    (which becomes the records' data) instead of copying them."""
    var records = FastaRecords()
    var n = len(bytes)
    var src = bytes.unsafe_ptr()
    var at = 0
    var i = 0
    while i < n:
        var c = src[unsafe_offset=i]
        if c == _HEADER:
            var name_start = i + 1
            var j = name_start
            while (
                j < n
                and src[unsafe_offset=j] != _SPACE
                and src[unsafe_offset=j] != _TAB
                and src[unsafe_offset=j] != _NEWLINE
                and src[unsafe_offset=j] != _RETURN
            ):
                j += 1
            records.names.append(String(StringSlice(unsafe_from_utf8=Span(bytes)[name_start:j])))
            records.starts.append(at)
            records.lengths.append(0)
            i = _line_end(src, j, n) + 1
            continue
        if len(records.names) == 0:
            if c == _NEWLINE or c == _RETURN or c == _SPACE:
                i += 1
                continue
            raise Error("FASTA input must start with a '>' header line")
        var j = _line_end(src, i, n)
        var end = j
        if end > i and src[unsafe_offset=end - 1] == _RETURN:
            end -= 1
        var length = end - i
        if length > 0:
            if at != i:
                var raw = src.unsafe_origin_cast[MutUntrackedOrigin]()
                unsafe_memmove(dest=raw.unsafe_offset(at), src=raw.unsafe_offset(i), count=length)
            at += length
            records.lengths[len(records.lengths) - 1] += length
        i = j + 1
    bytes.resize(unsafe_uninit_length=at)
    records.data = bytes^
    return records^


def read_file(path: String, threads: Int = 1) raises -> List[UInt8]:
    """The file's (or pipe's) bytes, inflated if they are gzipped."""
    var raw = _read_raw(path)
    if len(raw) >= 2 and raw[0] == 0x1F and raw[1] == 0x8B:
        var libdeflate = _load_libdeflate()
        if libdeflate:
            if _is_bgzf(Span(raw)):
                return _inflate_bgzf(Span(raw), libdeflate.value(), threads)
            return _inflate_gzip(Span(raw), libdeflate.value())
        return _gunzip(path)
    return raw^


def _read_raw(path: String) raises -> List[UInt8]:
    """Every byte of `path` up to end of file. The size `stat` reports is
    only a first guess (pipes report 0). A single read() call cannot return
    2 GiB or more on macOS, so reads go in blocks."""
    var guess = max(0, Int(stat(path).st_size))
    var out = List[UInt8](unsafe_uninit_length=max(guess, 1 << 16))
    var probe = List[UInt8](unsafe_uninit_length=1 << 16)
    var done = 0
    with open(path, "r") as f:
        while True:
            if done == len(out):
                # Full: probe for more before growing (a regular file ends here).
                var more = f.read(Span(probe))
                if more < 0:
                    raise Error(t"cannot read {path}")
                if more == 0:
                    break
                out.resize(unsafe_uninit_length=2 * len(out))
                unsafe_memcpy(dest=out.unsafe_ptr().unsafe_offset(done), src=probe.unsafe_ptr(), count=Int(more))
                done += Int(more)
                continue
            var room = min(1 << 30, len(out) - done)
            var got = f.read(Span(out)[done : done + room])
            if got < 0:
                raise Error(t"cannot read {path}")
            if got == 0:
                break
            done += Int(got)
    out.resize(unsafe_uninit_length=done)
    return out^


def _load_libdeflate() -> Optional[OwnedDLHandle]:
    comptime suffix = ".dylib" if CompilationTarget.is_macos() else ".so"
    var names: List[String] = [
        "@rpath/libdeflate.0.dylib" if CompilationTarget.is_macos() else "libdeflate.so.0",
        String("libdeflate") + suffix,
        "/opt/homebrew/lib/libdeflate.dylib",
    ]
    var prefix = getenv("CONDA_PREFIX")
    if prefix:
        names.append(prefix + "/lib/libdeflate" + suffix)
    for name in names:
        try:
            return OwnedDLHandle(name)
        except:
            pass
    return None


def _is_bgzf(data: Span[UInt8, _]) -> Bool:
    """Whether a gzip stream starts with a BGZF block (an FEXTRA `BC` field)."""
    return len(data) >= 18 and data[3] & 4 != 0 and data[12] == 66 and data[13] == 67


struct _InflateTask[data_origin: ImmOrigin, blocks_origin: ImmOrigin, lib_origin: ImmOrigin](ParallelTask):
    """Inflates one run of BGZF blocks into their slots of the output."""

    var data: ImmPointer[UInt8, Self.data_origin]
    var blocks: ImmPointer[Int, Self.blocks_origin]
    var dst: MutPointer[UInt8, MutUntrackedOrigin]
    var lib: ImmPointer[OwnedDLHandle, Self.lib_origin]
    var per_task: Int
    var count: Int
    var failed: MutPointer[Int, MutUntrackedOrigin]

    def __init__(
        out self,
        data: ImmPointer[UInt8, Self.data_origin],
        blocks: ImmPointer[Int, Self.blocks_origin],
        dst: MutPointer[UInt8, MutUntrackedOrigin],
        lib: ImmPointer[OwnedDLHandle, Self.lib_origin],
        per_task: Int,
        count: Int,
        failed: MutPointer[Int, MutUntrackedOrigin],
    ):
        self.data = data
        self.blocks = blocks
        self.dst = dst
        self.lib = lib
        self.per_task = per_task
        self.count = count
        self.failed = failed

    def run(self, index: Int):
        try:
            ref lib = self.lib[]
            var decompressor = lib.get_function[Int]("libdeflate_alloc_decompressor")()
            # One symbol lookup per task, not per block.
            var inflate = lib.get_function[Int32]("libdeflate_deflate_decompress")
            for b in range(index * self.per_task, min(self.count, (index + 1) * self.per_task)):
                # Block b: compressed deflate data, its size, output offset and size.
                var block = self.blocks.unsafe_offset(4 * b)
                var status = inflate(
                    decompressor,
                    self.data.unsafe_offset(block[unsafe_offset=0]),
                    block[unsafe_offset=1],
                    self.dst.unsafe_offset(block[unsafe_offset=2]),
                    block[unsafe_offset=3],
                    Int(0),
                )
                if status != 0:
                    self.failed[unsafe_offset=0] = 1
            lib.get_function[NoneType]("libdeflate_free_decompressor")(decompressor)
        except:
            self.failed[unsafe_offset=0] = 1


def _inflate_bgzf(data: Span[UInt8, _], lib: OwnedDLHandle, threads: Int) raises -> List[UInt8]:
    """Inflates a BGZF file block-parallel: block headers give every block's
    compressed range and output size up front."""
    var blocks = List[Int]()
    var total = 0
    var at = 0
    var n = len(data)
    while at < n:
        if at + 18 > n or data[at] != 0x1F or data[at + 1] != 0x8B:
            raise Error("corrupt BGZF block header")
        if not _is_bgzf(data[at:]):
            break  # a plain gzip member: the rest inflates serially
        var xlen = Int(data[at + 10]) | (Int(data[at + 11]) << 8)
        var bsize = -1
        var x = at + 12
        while x + 4 <= at + 12 + xlen and x + 4 <= n:
            var sub_len = Int(data[x + 2]) | (Int(data[x + 3]) << 8)
            if data[x] == 66 and data[x + 1] == 67 and sub_len == 2 and x + 6 <= n:
                bsize = Int(data[x + 4]) | (Int(data[x + 5]) << 8)
            x += 4 + sub_len
        # A block holds its header, the extra field, deflate data and an 8-byte trailer.
        if bsize < xlen + 19 or at + bsize + 1 > n:
            raise Error("corrupt BGZF block")
        var end = at + bsize + 1
        var isize = (
            Int(data[end - 4]) | (Int(data[end - 3]) << 8) | (Int(data[end - 2]) << 16) | (Int(data[end - 1]) << 24)
        )
        blocks.append(at + 12 + xlen)
        blocks.append(end - 8 - (at + 12 + xlen))
        blocks.append(total)
        blocks.append(isize)
        total += isize
        at = end
    var out = List[UInt8](unsafe_uninit_length=total)
    var count = len(blocks) // 4
    var per_task = max(1, ceildiv(count, max(1, threads) * 16))
    var failed = List[Int](length=1, fill=0)
    var task = _InflateTask(
        data.unsafe_ptr().as_imm(),
        blocks.unsafe_ptr().as_imm(),
        out.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        Pointer(to=lib),
        per_task,
        count,
        failed.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
    )
    if count > 0:
        parallel_for(task, ceildiv(count, per_task), threads)
    if failed[0] != 0 or (count == 0 and at >= n):
        raise Error("BGZF inflation failed")
    if at < n:
        out.extend(Span(_inflate_gzip(data[at:], lib)))
    return out^


def _inflate_gzip(data: Span[UInt8, _], lib: OwnedDLHandle) raises -> List[UInt8]:
    """Inflates every gzip member of `data` with libdeflate."""
    var out = List[UInt8]()
    var decompressor = lib.get_function[Int]("libdeflate_alloc_decompressor")()
    var at = 0
    var room = max(1 << 20, 4 * len(data))
    while at < len(data):
        var used = 0
        var produced = 0
        while True:
            var old = len(out)
            out.resize(unsafe_uninit_length=old + room)
            var status = lib.get_function[Int32]("libdeflate_gzip_decompress_ex")(
                decompressor,
                data.unsafe_ptr().unsafe_offset(at),
                len(data) - at,
                out.unsafe_ptr().unsafe_offset(old),
                room,
                Pointer(to=used),
                Pointer(to=produced),
            )
            out.resize(unsafe_uninit_length=old + (produced if status == 0 else 0))
            if status == 0:
                break
            if status != 3:  # LIBDEFLATE_INSUFFICIENT_SPACE
                lib.get_function[NoneType]("libdeflate_free_decompressor")(decompressor)
                raise Error("gzip inflation failed")
            room *= 2
        at += used
    lib.get_function[NoneType]("libdeflate_free_decompressor")(decompressor)
    return out^


def _gunzip(path: String) raises -> List[UInt8]:
    comptime library = "libz.dylib" if CompilationTarget.is_macos() else "libz.so.1"
    var zlib = OwnedDLHandle(library)
    var c_path = List[UInt8](path.as_bytes())
    c_path.append(0)
    var mode: List[UInt8] = [UInt8(ord("r")), UInt8(ord("b")), 0]
    var gz = zlib.get_function[Int]("gzopen")(c_path.unsafe_ptr(), mode.unsafe_ptr())
    if gz == 0:
        raise Error(t"cannot open {path}")
    _ = zlib.get_function[Int32]("gzbuffer")(gz, UInt32(1 << 20))
    comptime CHUNK = 1 << 24
    var out = List[UInt8]()
    while True:
        var old = len(out)
        out.resize(unsafe_uninit_length=old + CHUNK)
        var got = Int(zlib.get_function[Int32]("gzread")(gz, out.unsafe_ptr().unsafe_offset(old), UInt32(CHUNK)))
        if got < 0:
            _ = zlib.get_function[Int32]("gzclose")(gz)
            raise Error(t"error inflating {path}")
        out.resize(unsafe_uninit_length=old + got)
        if got < CHUNK:
            break
    _ = zlib.get_function[Int32]("gzclose")(gz)
    return out^
