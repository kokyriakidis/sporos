"""Minimal, fast FASTA reading: whole file in, records out."""

from std.bit import count_trailing_zeros
from std.memory import pack_bits, unsafe_memcpy

comptime _HEADER = UInt8(ord(">"))
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
    var bytes: List[UInt8]
    with open(path, "r") as f:
        bytes = f.read_bytes()
    return parse_fasta(Span(bytes))
