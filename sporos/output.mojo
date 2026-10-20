"""Buffered byte output to a file descriptor with fast integer formatting."""

from std.ffi import external_call
from std.sys import exit, stderr
from std.sys.info import CompilationTarget

comptime _FLUSH_AT = 1 << 20


struct OutputBuffer(Movable):
    """Accumulates bytes and writes them to `fd` in large blocks; with
    `fd = -1` it only accumulates (read them from `buf`)."""

    var buf: List[UInt8]
    var fd: Int32

    def __init__(out self, fd: Int32 = 1):
        self.buf = List[UInt8](capacity=_FLUSH_AT + 4096)
        self.fd = fd

    def write(mut self, text: StringSlice):
        self.write_bytes(text.as_bytes())

    def write_bytes(mut self, bytes: Span[UInt8, _]):
        if self.fd >= 0 and len(bytes) >= _FLUSH_AT:
            # Large blocks go out directly, without a copy into the buffer.
            self.flush()
            _write_all(self.fd, bytes.unsafe_ptr(), len(bytes))
            return
        self.buf.extend(bytes)
        if len(self.buf) >= _FLUSH_AT and self.fd >= 0:
            self.flush()

    @inline(.always)
    def write_byte(mut self, b: Int):
        self.buf.append(UInt8(b))

    def write_int(mut self, value: Int, width: Int = 0, left: Bool = False):
        """Writes `value` in decimal, padded with spaces to `width`."""
        var digits = InlineDigits(value)
        var pad = max(0, width - digits.count)
        if not left:
            for _ in range(pad):
                self.buf.append(UInt8(ord(" ")))
        for i in range(digits.count):
            self.buf.append(digits.bytes[digits.count - 1 - i])
        if left:
            for _ in range(pad):
                self.buf.append(UInt8(ord(" ")))
        if len(self.buf) >= _FLUSH_AT and self.fd >= 0:
            self.flush()

    def flush(mut self):
        if self.fd < 0:
            return
        _write_all(self.fd, self.buf.unsafe_ptr(), len(self.buf))
        self.buf.clear()


def _write_all(fd: Int32, data: ImmPointer[UInt8, _], length: Int):
    """Writes every byte, retrying interrupted writes; any other failure (a
    full disk, a closed pipe) ends the program with an error, never silently
    truncated output."""
    var done = 0
    while done < length:
        var n = external_call["write", Int](fd, data.unsafe_offset(done), length - done)
        if n > 0:
            done += n
            continue
        if n < 0 and _errno() == _EINTR:
            continue
        print("sporos: write failed (errno ", _errno(), ")", sep="", file=stderr)
        exit(1)


comptime _EINTR: Int32 = 4


def _errno() -> Int32:
    comptime if CompilationTarget.is_macos():
        return external_call["__error", MutPointer[Int32, MutAnyOrigin]]()[]
    else:
        return external_call["__errno_location", MutPointer[Int32, MutAnyOrigin]]()[]


@fieldwise_init
struct InlineDigits(TrivialRegisterPassable):
    """Decimal digits of an integer, least significant first."""

    var bytes: SIMD[DType.uint8, 32]
    var count: Int

    def __init__(out self, value: Int):
        self.bytes = SIMD[DType.uint8, 32](0)
        self.count = 0
        var v = value
        var negative = v < 0
        if negative:
            v = -v
        while True:
            self.bytes[self.count] = UInt8(ord("0") + v % 10)
            self.count += 1
            v //= 10
            if v == 0:
                break
        if negative:
            self.bytes[self.count] = UInt8(ord("-"))
            self.count += 1


@inline(.always)
def put_padded(dst: MutPointer[UInt8, _], value: Int, width: Int) -> Int:
    """Writes non-negative `value` right-aligned in `width` columns (wider
    when it needs more digits); returns the bytes written."""
    var digits = 1
    var v = value
    while v >= 10:
        v //= 10
        digits += 1
    var n = max(width, digits)
    for i in range(n - digits):
        dst[unsafe_offset=i] = UInt8(ord(" "))
    v = value
    for i in reversed(range(n - digits, n)):
        dst[unsafe_offset=i] = UInt8(ord("0") + v % 10)
        v //= 10
    return n
