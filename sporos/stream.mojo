"""A file's bytes a chunk at a time, inflated if gzipped, so that a read set never has to fit in memory whole.

How each kind of file streams:
- plain: read a chunk at a time;
- BGZF (bgzip): compressed chunks, whose complete blocks are inflated in parallel with libdeflate, the partial
  last block carried into the next chunk; a plain gzip member after them continues through zlib;
- other gzip at most `whole_limit` bytes compressed: inflated in one call by libdeflate;
- larger gzip: decoded on `threads` threads by rapidgzip (rapidgzip-core, Rust, through the C interface in
  rapidgzip-ffi: `libsporos_rapidgzip`), which finds DEFLATE blocks mid-stream and decodes them in parallel;
  without that library, through zlib's `gzread`, one thread, as hifiasm's kseq reads it;
- pipes and anything `stat` gives no size for: read whole (`read_file`).
"""

from std.ffi import OwnedDLHandle, external_call
from std.os import getenv, stat
from std.os.path import dirname
from std.sys import argv
from std.sys.info import CompilationTarget

from sporos.fasta import _inflate_bgzf_into, _is_bgzf, _load_libdeflate, read_file

comptime DEFAULT_CHUNK = 1 << 26
"""Bytes of inflated data a `next` call aims to deliver."""
comptime WHOLE_GZIP_LIMIT = 1 << 30
"""A plain gzip file at most this large (compressed) is inflated whole, by libdeflate."""

comptime _PLAIN = 0
comptime _BGZF = 1
comptime _ZLIB = 2
comptime _WHOLE = 3
comptime _DONE = 4
comptime _RAPIDGZIP = 5

comptime RAPIDGZIP_LIBRARY = "libsporos_rapidgzip" + (".dylib" if CompilationTarget.is_macos() else ".so")
comptime _RAPIDGZIP_ABI = 1


def load_rapidgzip() -> Optional[OwnedDLHandle]:
    """`libsporos_rapidgzip` (`pixi run build-rapidgzip`), if it can be found: at $SPOROS_RAPIDGZIP_LIB, beside the
    running program, in the environment pixi sets up ($CONDA_PREFIX/lib), or on the loader's path; and only if its
    interface is the one this code speaks."""
    var names = List[String]()
    var env = getenv("SPOROS_RAPIDGZIP_LIB")
    if env:
        names.append(env)
    var args = argv()
    if len(args) > 0:
        var dir = dirname(String(args[0]))
        names.append((dir if dir else ".") + "/" + RAPIDGZIP_LIBRARY)
    var prefix = getenv("CONDA_PREFIX")
    if prefix:
        names.append(prefix + "/lib/" + RAPIDGZIP_LIBRARY)
    names.append(RAPIDGZIP_LIBRARY)
    for name in names:
        try:
            var lib = OwnedDLHandle(name)
            if Int(lib.get_function[UInt32]("srgz_abi")()) == _RAPIDGZIP_ABI:
                return lib^
        except:
            pass
    return None


struct ByteStream(Movable):
    """The bytes of `path`, delivered by `next` a chunk at a time."""

    var path: String
    var threads: Int
    var chunk: Int
    var kind: Int
    var file: Optional[FileHandle]
    var libdeflate: Optional[OwnedDLHandle]
    var zlib: Optional[OwnedDLHandle]
    var gz: Int
    """zlib's `gzFile`, 0 when not open."""
    var rapidgzip: Optional[OwnedDLHandle]
    var rgz: Int
    """rapidgzip's stream, 0 when not open."""
    var carry: List[UInt8]
    """BGZF bytes read but not inflated yet: the start of a block."""
    var carry_offset: Int
    """Where in the file `carry` starts."""
    var at_eof: Bool

    def __init__(
        out self,
        path: String,
        threads: Int = 1,
        chunk: Int = DEFAULT_CHUNK,
        whole_limit: Int = WHOLE_GZIP_LIMIT,
        use_rapidgzip: Bool = True,
    ) raises:
        """`path`'s bytes, `chunk` of them per `next`, on `threads` threads where the kind of file allows;
        `use_rapidgzip=False` keeps large gzip on zlib even where rapidgzip is found. $SPOROS_GZIP_WHOLE_LIMIT, when
        set, replaces `whole_limit` (0: every gzip file streams)."""
        self.path = path
        self.threads = max(1, threads)
        self.chunk = max(1 << 12, chunk)
        self.file = None
        self.libdeflate = None
        self.zlib = None
        self.gz = 0
        self.rapidgzip = None
        self.rgz = 0
        self.carry = List[UInt8]()
        self.carry_offset = 0
        self.at_eof = False
        var limit = whole_limit
        var env = getenv("SPOROS_GZIP_WHOLE_LIMIT")
        if env:
            limit = Int(env)
        var size = Int(stat(path).st_size)
        if size <= 0:
            self.kind = _WHOLE
            return
        var head = List[UInt8](length=18, fill=0)
        var got: Int
        with open(path, "r") as f:
            got = Int(f.read(Span(head)))
        var gzipped = got >= 2 and head[0] == 0x1F and head[1] == 0x8B
        if not gzipped:
            self.kind = _PLAIN
            self.file = open(path, "r")
            return
        self.libdeflate = _load_libdeflate()
        if self.libdeflate and got >= 18 and _is_bgzf(Span(head)[:got]):
            self.kind = _BGZF
            self.file = open(path, "r")
        elif self.libdeflate and size <= limit:
            self.kind = _WHOLE
        else:
            if use_rapidgzip:
                self.rapidgzip = load_rapidgzip()
            if self.rapidgzip:
                self.kind = _RAPIDGZIP
                self._open_rapidgzip()
            else:
                self.kind = _ZLIB
                self._open_zlib(0)

    def __deinit__(deinit self):
        if self.gz != 0 and self.zlib:
            try:
                _ = self.zlib.value().get_function[Int32]("gzclose")(self.gz)
            except:
                pass
        if self.rgz != 0 and self.rapidgzip:
            try:
                self.rapidgzip.value().get_function[NoneType]("srgz_close")(self.rgz)
            except:
                pass

    def _rapidgzip_error(self) raises -> String:
        var length = 0
        var text = self.rapidgzip.value().get_function[MutPointer[UInt8, MutUntrackedOrigin]]("srgz_error")(
            Pointer(to=length)
        )
        return String(unsafe_from_utf8=Span(unsafe_ptr=text, length=length))

    def _open_rapidgzip(mut self) raises:
        var path = self.path.as_bytes()
        var handle = 0
        var status = self.rapidgzip.value().get_function[Int64]("srgz_open")(
            path.unsafe_ptr(), len(path), UInt32(self.threads), Pointer(to=handle)
        )
        if status < 0:
            var message = self._rapidgzip_error()
            raise Error(message)
        self.rgz = handle
        self.kind = _RAPIDGZIP

    def _next_rapidgzip(mut self, mut out: List[UInt8]) raises -> Bool:
        if self.rgz == 0:
            return False
        ref lib = self.rapidgzip.value()
        var old = len(out)
        out.resize(unsafe_uninit_length=old + self.chunk)
        var got = Int(lib.get_function[Int64]("srgz_read")(self.rgz, out.unsafe_ptr().unsafe_offset(old), self.chunk))
        if got < 0:
            out.resize(unsafe_uninit_length=old)
            var message = self._rapidgzip_error()
            raise Error(message)
        out.resize(unsafe_uninit_length=old + got)
        if got < self.chunk:
            # The end, verified by the library: every member's CRC32 and size.
            lib.get_function[NoneType]("srgz_close")(self.rgz)
            self.rgz = 0
        return got > 0

    def _open_zlib(mut self, offset: Int) raises:
        """Continues the file through zlib from byte `offset`."""
        comptime library = "libz.dylib" if CompilationTarget.is_macos() else "libz.so.1"
        self.zlib = OwnedDLHandle(library)
        # zlib gets its own descriptor, positioned at `offset`, so that it can carry on mid-file.
        var f = open(self.path, "r")
        _ = f.seek(offset)
        var fd = external_call["dup", Int32](Int32(f._get_raw_fd()))
        f.close()
        if fd < 0:
            raise Error(t"cannot open {self.path}")
        var mode: List[UInt8] = [UInt8(ord("r")), UInt8(ord("b")), 0]
        self.gz = self.zlib.value().get_function[Int]("gzdopen")(fd, mode.unsafe_ptr())
        if self.gz == 0:
            raise Error(t"cannot open {self.path} through zlib")
        _ = self.zlib.value().get_function[Int32]("gzbuffer")(self.gz, UInt32(1 << 20))
        self.kind = _ZLIB

    def decoder(self) -> String:
        """What reads the file: "plain", "bgzf" (libdeflate, block-parallel), "libdeflate" (whole), "zlib" or
        "rapidgzip" ("done" once used up)."""
        if self.kind == _PLAIN:
            return "plain"
        if self.kind == _BGZF:
            return "bgzf"
        if self.kind == _ZLIB:
            return "zlib"
        if self.kind == _WHOLE:
            return "libdeflate"
        if self.kind == _RAPIDGZIP:
            return "rapidgzip"
        return "done"

    def next(mut self, mut out: List[UInt8]) raises -> Bool:
        """Appends the next chunk of bytes to `out`; False once the file is used up (nothing appended)."""
        if self.kind == _PLAIN:
            return self._next_plain(out)
        if self.kind == _BGZF:
            return self._next_bgzf(out)
        if self.kind == _ZLIB:
            return self._next_zlib(out)
        if self.kind == _RAPIDGZIP:
            return self._next_rapidgzip(out)
        if self.kind == _WHOLE:
            self.kind = _DONE
            if len(out) == 0:
                out = read_file(self.path, self.threads)
            else:
                out.extend(Span(read_file(self.path, self.threads)))
            return True
        return False

    @staticmethod
    def _read(mut file: FileHandle, path: String, mut out: List[UInt8], want: Int, mut at_eof: Bool) raises -> Int:
        """Appends up to `want` bytes of the file to `out`; how many (fewer only at its end)."""
        var old = len(out)
        out.resize(unsafe_uninit_length=old + want)
        var done = 0
        while done < want:
            var got = Int(file.read(Span(out)[old + done : old + want]))
            if got < 0:
                raise Error(t"cannot read {path}")
            if got == 0:
                at_eof = True
                break
            done += got
        out.resize(unsafe_uninit_length=old + done)
        return done

    def _next_plain(mut self, mut out: List[UInt8]) raises -> Bool:
        if self.at_eof:
            return False
        return Self._read(self.file.value(), self.path, out, self.chunk, self.at_eof) > 0

    def _next_zlib(mut self, mut out: List[UInt8]) raises -> Bool:
        if self.gz == 0:
            return False
        ref zlib = self.zlib.value()
        var old = len(out)
        out.resize(unsafe_uninit_length=old + self.chunk)
        var done = 0
        while done < self.chunk:
            var room = min(self.chunk - done, 1 << 30)
            var got = Int(
                zlib.get_function[Int32]("gzread")(self.gz, out.unsafe_ptr().unsafe_offset(old + done), UInt32(room))
            )
            if got < 0:
                raise Error(t"error inflating {self.path}")
            if got == 0:
                _ = zlib.get_function[Int32]("gzclose")(self.gz)
                self.gz = 0
                break
            done += got
        out.resize(unsafe_uninit_length=old + done)
        return done > 0

    def _next_bgzf(mut self, mut out: List[UInt8]) raises -> Bool:
        """Inflates whole blocks, in parallel, until they make a chunk: the blocks' trailers give their inflated
        sizes, so a chunk stays near `chunk` however well the data compresses. More of the file is read only
        when the blocks at hand fall short."""
        while True:
            var whole = _complete_bgzf_blocks(Span(self.carry), self.chunk)
            if not whole[2] and not whole[1] and not self.at_eof:
                # Compressed reads of a third of the chunk: BGZF reads inflate about three- to fourfold.
                _ = Self._read(self.file.value(), self.path, self.carry, max(1 << 16, self.chunk // 3), self.at_eof)
                continue
            var appended = False
            if whole[0] > 0:
                var before = len(out)
                _inflate_bgzf_into(Span(self.carry)[: whole[0]], self.libdeflate.value(), self.threads, out)
                appended = len(out) > before
                self.carry_offset += whole[0]
                var rest = len(self.carry) - whole[0]
                var cp = self.carry.unsafe_ptr()
                for i in range(rest):
                    cp[unsafe_offset=i] = cp[unsafe_offset=whole[0] + i]
                self.carry.resize(unsafe_uninit_length=rest)
            if whole[1]:
                # A plain gzip member: zlib takes the file from there.
                self.file = None
                self.carry = List[UInt8]()
                self._open_zlib(self.carry_offset)
                return self._next_zlib(out) or appended
            if appended:
                return True
            if whole[0] > 0:
                continue  # only empty blocks (the end-of-file marker): look further
            if len(self.carry) > 0:
                raise Error(t"{self.path}: truncated or corrupt BGZF data")
            self.kind = _DONE
            return False


def _complete_bgzf_blocks(data: Span[UInt8, _], limit: Int) -> Tuple[Int, Bool, Bool]:
    """How many leading bytes of `data` are whole BGZF blocks, taking blocks until they inflate to `limit` bytes;
    whether a gzip member that is not BGZF starts right after them; and whether `limit` was reached."""
    var n = len(data)
    var at = 0
    var inflated = 0
    while at + 18 <= n:
        if inflated >= limit:
            return (at, False, True)
        if data[at] != 0x1F or data[at + 1] != 0x8B:
            return (at, False, False)
        if not _is_bgzf(data[at:]):
            return (at, True, False)
        var xlen = Int(data[at + 10]) | (Int(data[at + 11]) << 8)
        var bsize = -1
        var x = at + 12
        while x + 4 <= at + 12 + xlen and x + 4 <= n:
            var sub_len = Int(data[x + 2]) | (Int(data[x + 3]) << 8)
            if data[x] == 66 and data[x + 1] == 67 and sub_len == 2 and x + 6 <= n:
                bsize = Int(data[x + 4]) | (Int(data[x + 5]) << 8)
            x += 4 + sub_len
        if bsize < 0 or at + bsize + 1 > n:
            return (at, False, False)
        var end = at + bsize + 1
        inflated += (
            Int(data[end - 4]) | (Int(data[end - 3]) << 8) | (Int(data[end - 2]) << 16) | (Int(data[end - 1]) << 24)
        )
        at = end
    return (at, False, inflated >= limit)
