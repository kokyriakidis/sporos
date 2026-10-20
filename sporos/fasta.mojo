"""Fast FASTA/FASTQ reading: whole file in (gzip transparently), records out.

Gzipped input is detected by its magic bytes and inflated with the system
zlib, loaded at run time, so sporos has no build-time dependency on it.
"""

from std.bit import count_trailing_zeros
from std.ffi import OwnedDLHandle
from std.memory import pack_bits, unsafe_memcpy
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
        while j < n and src[unsafe_offset=j] != _SPACE and src[unsafe_offset=j] != _TAB and src[
            unsafe_offset=j
        ] != _NEWLINE and src[unsafe_offset=j] != _RETURN:
            j += 1
        records.names.append(String(StringSlice(unsafe_from_utf8=bytes[i + 1 : j])))
        var seq_start = _line_end(src, j, n) + 1
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
            while j < n and src[unsafe_offset=j] != _SPACE and src[
                unsafe_offset=j
            ] != _TAB and src[unsafe_offset=j] != _NEWLINE and src[
                unsafe_offset=j
            ] != _RETURN:
                j += 1
            records.names.append(
                String(StringSlice(unsafe_from_utf8=bytes[name_start:j]))
            )
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


def read_fasta(path: String) raises -> FastaRecords:
    """Reads a FASTA or FASTQ file, gzipped or not."""
    var bytes = read_file(path)
    return parse_sequences(Span(bytes))


def read_file(path: String) raises -> List[UInt8]:
    """The file's bytes, inflated if it is gzipped."""
    var magic: List[UInt8]
    with open(path, "r") as f:
        magic = f.read_bytes(2)
    if len(magic) == 2 and magic[0] == 0x1F and magic[1] == 0x8B:
        return _gunzip(path)
    with open(path, "r") as f:
        return f.read_bytes()


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
        var got = Int(
            zlib.get_function[Int32]("gzread")(gz, out.unsafe_ptr().unsafe_offset(old), UInt32(CHUNK))
        )
        if got < 0:
            _ = zlib.get_function[Int32]("gzclose")(gz)
            raise Error(t"error inflating {path}")
        out.resize(unsafe_uninit_length=old + got)
        if got < CHUNK:
            break
    _ = zlib.get_function[Int32]("gzclose")(gz)
    return out^
