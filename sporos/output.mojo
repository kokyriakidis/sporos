"""Buffered byte output to a file descriptor with fast integer formatting."""

from std.ffi import external_call

comptime _FLUSH_AT = 1 << 20


struct OutputBuffer(Movable):
    """Accumulates bytes and writes them to `fd` in large blocks."""

    var buf: List[UInt8]
    var fd: Int32

    def __init__(out self, fd: Int32 = 1):
        self.buf = List[UInt8](capacity=_FLUSH_AT + 4096)
        self.fd = fd

    def write(mut self, text: StringSlice):
        self.write_bytes(text.as_bytes())

    def write_bytes(mut self, bytes: Span[UInt8, _]):
        self.buf.extend(bytes)
        if len(self.buf) >= _FLUSH_AT:
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
        if len(self.buf) >= _FLUSH_AT:
            self.flush()

    def flush(mut self):
        var done = 0
        while done < len(self.buf):
            var n = external_call["write", Int](
                self.fd, self.buf.unsafe_ptr().unsafe_offset(done), len(self.buf) - done
            )
            if n <= 0:
                break
            done += n
        self.buf.clear()


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
