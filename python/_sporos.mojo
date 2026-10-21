"""The Python extension under the `sporos` package (see `sporos/__init__.py`).

Each function takes sequences as ASCII `str`, read in place, and returns positions and matches as
`bytes` of native `uint32`s (SMEMs and seed hits as `uint64`s, each query's offsets first), which the
package wraps as `array`s and tuples; this module is not meant to be called directly.

    pixi run build-python   # build/python: the package, this module and the runtime libraries it loads
"""

from std.os import abort
from std.python import Python, PythonObject
from std.python._cpython import PyObjectPtr, Py_ssize_t
from std.python.bindings import PythonModuleBuilder
from std.sys import size_of

from sporos.fmindex import FMIndex
from sporos.index import SeedIndex
from sporos.mems import MemIndex, MemParams, filter_by_occurrences, find_mems
from sporos.minimizers import minimizer_positions
from sporos.parallel import MAX_THREADS, default_threads
from sporos.syncmers import syncmer_positions
from sporos.text import SeqText

comptime MODE_MAXMATCH = 0
"""A matches call's mode for every MEM (`SPOROS_MAXMATCH`)."""
comptime MODE_MUM = 1
"""A matches call's mode for MUMs (`SPOROS_MUM`)."""
comptime MODE_MUMREFERENCE = 2
"""A matches call's mode for MEMs unique in the reference (`SPOROS_MUMREFERENCE`)."""


@export
def PyInit__sporos() abi("C") -> PythonObject:
    """The `_sporos` module, its functions registered, which Python calls on import; a failure aborts."""
    try:
        var module = PythonModuleBuilder("_sporos")
        module.def_function[py_minimizers]("minimizers")
        module.def_function[py_syncmers]("syncmers")
        _ = (
            module.add_type[PyIndex]("Index")
            .def_py_init[PyIndex.py_init]()
            .def_method[PyIndex.matches]("matches")
            .def_method[PyIndex.sequences]("sequences")
        )
        _ = (
            module.add_type[PyFMIndex]("FMIndex")
            .def_py_init[PyFMIndex.py_init]()
            .def_method[PyFMIndex.sequences]("sequences")
            .def_method[PyFMIndex.name]("name")
            .def_method[PyFMIndex.save]("save")
            .def_method[PyFMIndex.smems]("smems")
        )
        _ = (
            module.add_type[PySeedIndex]("SeedIndex")
            .def_py_init[PySeedIndex.py_init]()
            .def_method[PySeedIndex.occurrence_cutoff]("occurrence_cutoff")
            .def_method[PySeedIndex.lookup]("lookup")
        )
        return module.finalize()
    except error:
        abort(String("sporos: the extension failed to load: ", error))


def bases(seq: PythonObject) -> Span[UInt8, ImmutAnyOrigin]:
    """An ASCII `str`'s bytes, in place: alive while the caller holds `seq`."""
    return Python().as_string_slice(seq).as_bytes()


def spans_of(seqs: PythonObject) raises -> List[Span[UInt8, ImmutAnyOrigin]]:
    """A list of ASCII `str`s as spans, alive while the caller holds the list."""
    var spans = List[Span[UInt8, ImmutAnyOrigin]]()
    for seq in seqs:
        spans.append(bases(seq))
    return spans^


def threads_of(count: PythonObject) raises -> Int:
    """`count`, or every physical core for zero or less."""
    var threads = Int(py=count)
    return min(threads, MAX_THREADS) if threads > 0 else default_threads()


def bytes_of[T: TrivialRegisterPassable](items: List[T]) raises -> PythonObject:
    """`items`' memory as a Python `bytes`."""
    ref cpython = Python().cpython()
    var made = cpython.lib.call["PyBytes_FromStringAndSize", PyObjectPtr](
        items.unsafe_ptr().unsafe_bitcast[UInt8](), Py_ssize_t(len(items) * size_of[T]())
    )
    if not made:
        raise cpython.get_error()
    return PythonObject(from_owned=made)


def py_minimizers(
    seq: PythonObject, k: PythonObject, w: PythonObject, canonical: PythonObject, threads: PythonObject
) raises -> PythonObject:
    """The (w, k)-minimizer positions of `seq`, as `uint32` bytes."""
    return bytes_of(minimizer_positions(bases(seq), Int(py=k), Int(py=w), Bool(py=canonical), threads_of(threads)))


def py_syncmers(
    seq: PythonObject,
    k: PythonObject,
    s: PythonObject,
    open_offsets: PythonObject,
    canonical: PythonObject,
    threads: PythonObject,
) raises -> PythonObject:
    """The syncmer positions of `seq`, as `uint32` bytes: closed for an empty list of `open_offsets`,
    else open at any of them."""
    var opened = List[Int]()
    for t in open_offsets:
        opened.append(Int(py=t))
    return bytes_of(
        syncmer_positions(
            bases(seq),
            Int(py=k),
            Int(py=s),
            Bool(py=canonical),
            None,
            threads_of(threads),
            open_offsets=opened,
        )
    )


struct PyIndex(Movable, Writable):
    """A reference's packed text and its anchor index, which reads it (see `sporos.Index`)."""

    var text: SeqText
    var index: MemIndex

    def __init__(out self, var text: SeqText, min_len: Int, threads: Int) raises:
        var params = MemParams.choose(min_len, len(text), len(text))
        self.index = MemIndex(text, params, threads)
        self.text = text^

    @staticmethod
    def py_init(out self: PyIndex, args: PythonObject, kwargs: PythonObject) raises:
        """`Index(seqs, min_len, threads)`, the sequences a list of ASCII `str`s."""
        if len(args) != 3:
            raise Error("Index takes the reference sequences, the minimum match length and the threads")
        var spans = spans_of(args[0])
        if not spans:
            raise Error("an index needs at least one reference sequence")
        # `args` holds every `str` the spans read while the text copies them.
        self = PyIndex(SeqText(spans), Int(py=args[1]), threads_of(args[2]))

    def write_to(self, mut writer: Some[Writer]):
        writer.write("Index(", len(self.text.starts), " sequences, min_len=", self.index.params.min_len, ")")

    def write_repr_to(self, mut writer: Some[Writer]):
        self.write_to(writer)

    @staticmethod
    def sequences(self_ptr: MutPointer[Self, MutAnyOrigin]) -> PythonObject:
        """How many reference sequences the index holds."""
        return PythonObject(len(self_ptr[].text.starts))

    @staticmethod
    def matches(
        self_ptr: MutPointer[Self, MutAnyOrigin],
        query: PythonObject,
        mode: PythonObject,
        reverse: PythonObject,
        limits: PythonObject,
        threads: PythonObject,
    ) raises -> PythonObject:
        """The matches of `query` (or its reverse complement) under `mode`, with `limits` the
        `[max_occurrences, min_occurrences]` of `Reference.mems`, as `uint32` bytes, four per match:
        reference sequence, position in it, query position, length."""
        ref reference = self_ptr[].text
        var spans = List[Span[UInt8, ImmutAnyOrigin]]()
        spans.append(bases(query))
        if len(spans[0]) >= Int(UInt32.MAX):
            raise Error("a query holds fewer than 2^32 bases")
        var kind = Int(py=mode)
        var max_occurrences = Int(py=limits[0])
        var min_occurrences = Int(py=limits[1])
        if max_occurrences < 0 or min_occurrences < 1:
            raise Error("max_occurrences is at least 0, min_occurrences at least 1")
        var cap = max_occurrences if max_occurrences > 0 else (0 if kind == MODE_MAXMATCH else 1)
        if cap > 0 and min_occurrences > cap:
            raise Error(t"min_occurrences {min_occurrences} is above the most allowed ({cap})")
        var text = SeqText(spans, reverse_complement=Bool(py=reverse))
        var mems = find_mems(self_ptr[].index, reference, text, threads_of(threads))
        var kept = filter_by_occurrences(
            mems^, text, cap, cap if kind == MODE_MUM else 0, min_occurrences, threads_of(threads)
        )
        var fields = List[UInt32](capacity=4 * len(kept))
        for m in kept:
            var seq = reference.sequence_of(Int(m.ref_pos))
            fields.append(UInt32(seq))
            fields.append(m.ref_pos - UInt32(reference.starts[seq]))
            fields.append(m.query_pos)
            fields.append(m.length)
        return bytes_of(fields)


struct PyFMIndex(Movable, Writable):
    """An FM-index of reference sequences, for SMEMs (see `sporos.FMIndex`)."""

    var index: FMIndex

    def __init__(out self, var index: FMIndex):
        self.index = index^

    @staticmethod
    def py_init(out self: PyFMIndex, args: PythonObject, kwargs: PythonObject) raises:
        """`FMIndex(seqs, names, threads)` builds, `FMIndex(path, threads)` loads a saved index."""
        if len(args) == 2:
            self = PyFMIndex(FMIndex.load(String(py=args[0]), threads_of(args[1])))
            return
        if len(args) != 3:
            raise Error("FMIndex takes the sequences, their names and the threads, or a path and the threads")
        var spans = spans_of(args[0])
        if not spans:
            raise Error("an index needs at least one reference sequence")
        var names = List[String]()
        for name in args[1]:
            names.append(String(py=name))
        self = PyFMIndex(FMIndex.build(spans, names, threads_of(args[2])))

    def write_to(self, mut writer: Some[Writer]):
        writer.write("FMIndex(", len(self.index.names), " sequences)")

    def write_repr_to(self, mut writer: Some[Writer]):
        self.write_to(writer)

    @staticmethod
    def sequences(self_ptr: MutPointer[Self, MutAnyOrigin]) -> PythonObject:
        """How many reference sequences the index holds."""
        return PythonObject(len(self_ptr[].index.names))

    @staticmethod
    def name(self_ptr: MutPointer[Self, MutAnyOrigin], i: PythonObject) raises -> PythonObject:
        """Sequence `i`'s name."""
        return PythonObject(self_ptr[].index.names[Int(py=i)])

    @staticmethod
    def save(self_ptr: MutPointer[Self, MutAnyOrigin], path: PythonObject) raises -> PythonObject:
        """Writes the index to `path`."""
        self_ptr[].index.save(String(py=path))
        return PythonObject(None)

    @staticmethod
    def smems(
        self_ptr: MutPointer[Self, MutAnyOrigin], queries: PythonObject, settings: PythonObject
    ) raises -> PythonObject:
        """The SMEMs of every query, `settings` being `[min_len, max_hits, min_occurrences,
        max_occurrences, threads]`, as `uint64` bytes: the queries' SMEM offsets (one more than the
        queries), four fields per SMEM (start, end, count, first hit), three per hit (sequence,
        position, reverse)."""
        var spans = spans_of(queries)
        for span in spans:
            if len(span) >= Int(UInt32.MAX):
                raise Error("a query holds fewer than 2^32 bases")
        var found = self_ptr[].index.smems(
            spans,
            Int(py=settings[0]),
            Int(py=settings[1]),
            Int(py=settings[2]),
            Int(py=settings[3]),
            threads_of(settings[4]),
        )
        var fields = List[UInt64](capacity=len(found.offsets) + 4 * len(found.smems) + 3 * len(found.hits))
        for o in found.offsets:
            fields.append(UInt64(o))
        for smem in found.smems:
            fields.append(UInt64(smem.start))
            fields.append(UInt64(smem.end))
            fields.append(UInt64(smem.count))
            fields.append(UInt64(smem.first))
        for h in found.hits:
            fields.append(UInt64(h.sequence))
            fields.append(UInt64(h.position))
            fields.append(UInt64(h.reverse))
        return bytes_of(fields)


struct PySeedIndex(Movable, Writable):
    """A seed index of reference sequences (see `sporos.SeedIndex`)."""

    var index: SeedIndex

    def __init__(out self, var index: SeedIndex):
        self.index = index^

    @staticmethod
    def py_init(out self: PySeedIndex, args: PythonObject, kwargs: PythonObject) raises:
        """`SeedIndex(seqs, k, w, scheme, open_offsets, threads)`."""
        if len(args) != 6:
            raise Error("SeedIndex takes the sequences, k, w, the scheme, the open offsets and the threads")
        var spans = spans_of(args[0])
        if not spans:
            raise Error("an index needs at least one reference sequence")
        var opened = List[Int]()
        for t in args[4]:
            opened.append(Int(py=t))
        self = PySeedIndex(
            SeedIndex(
                spans,
                Int(py=args[1]),
                Int(py=args[2]),
                Int(py=args[3]),
                threads_of(args[5]),
                open_offsets=opened,
            )
        )

    def write_to(self, mut writer: Some[Writer]):
        writer.write("SeedIndex(", len(self.index.text.starts), " sequences, k=", self.index.k, ")")

    def write_repr_to(self, mut writer: Some[Writer]):
        self.write_to(writer)

    @staticmethod
    def occurrence_cutoff(self_ptr: MutPointer[Self, MutAnyOrigin], keep_fraction: PythonObject) raises -> PythonObject:
        """The smallest occurrence count that `keep_fraction` of the indexed k-mers stay within."""
        return PythonObject(self_ptr[].index.occurrence_cutoff(Float64(py=keep_fraction)))

    @staticmethod
    def lookup(
        self_ptr: MutPointer[Self, MutAnyOrigin],
        reads: PythonObject,
        all_kmers: PythonObject,
        max_occurrences: PythonObject,
        threads: PythonObject,
    ) raises -> PythonObject:
        """The hits of every read, as `uint64` bytes: the reads' hit offsets (one more than the
        reads), then five fields per hit (query position, sequence, position, occurrences,
        reverse)."""
        var spans = spans_of(reads)
        for span in spans:
            if len(span) >= Int(UInt32.MAX):
                raise Error("a read holds fewer than 2^32 bases")
        var found = self_ptr[].index.lookup_batch(
            spans, Bool(py=all_kmers), Int(py=max_occurrences), threads_of(threads)
        )
        ref text = self_ptr[].index.text
        var fields = List[UInt64](capacity=len(found.offsets) + 5 * len(found.hits))
        for o in found.offsets:
            fields.append(UInt64(o))
        for h in found.hits:
            var seq = text.sequence_of(Int(h.ref_pos))
            fields.append(UInt64(h.query_pos))
            fields.append(UInt64(seq))
            fields.append(UInt64(Int(h.ref_pos) - text.starts[seq]))
            fields.append(UInt64(h.occurrences))
            fields.append(UInt64(h.reverse))
        return bytes_of(fields)
