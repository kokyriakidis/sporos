"""The Python extension under the `sporos` package (see `sporos/__init__.py`).

Each function takes sequences as ASCII `str`, read in place, and returns positions and matches as
`bytes` of native `uint32`s, which the package wraps as `array('I')` and tuples; this module is not
meant to be called directly.

    pixi run build-python   # build/python: the package, this module and the runtime libraries it loads
"""

from std.os import abort
from std.python import Python, PythonObject
from std.python._cpython import PyObjectPtr, Py_ssize_t
from std.python.bindings import PythonModuleBuilder
from std.sys import size_of

from sporos.mems import Mem, MemIndex, MemParams, filter_mums, find_mems
from sporos.minimizers import minimizer_positions
from sporos.parallel import default_threads
from sporos.syncmers import syncmer_positions
from sporos.text import SeqText

comptime MODE_MUM = 1
"""A matches call's mode for MUMs (`SPOROS_MUM`); any but these two is every MEM."""
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
        return module.finalize()
    except error:
        abort(String("sporos: the extension failed to load: ", error))


def bases(seq: PythonObject) -> Span[UInt8, ImmutAnyOrigin]:
    """An ASCII `str`'s bytes, in place: alive while the caller holds `seq`."""
    return Python().as_string_slice(seq).as_bytes()


def threads_of(count: PythonObject) raises -> Int:
    """`count`, or every physical core for zero or less."""
    var threads = Int(py=count)
    return threads if threads > 0 else default_threads()


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
    open_offset: PythonObject,
    canonical: PythonObject,
    threads: PythonObject,
) raises -> PythonObject:
    """The syncmer positions of `seq`, closed for a negative `open_offset`, as `uint32` bytes."""
    var offset = Int(py=open_offset)
    return bytes_of(
        syncmer_positions(
            bases(seq),
            Int(py=k),
            Int(py=s),
            Bool(py=canonical),
            Optional[Int](offset) if offset >= 0 else None,
            threads_of(threads),
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
        var seqs = args[0]
        var spans = List[Span[UInt8, ImmutAnyOrigin]]()
        for seq in seqs:
            spans.append(bases(seq))
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
        threads: PythonObject,
    ) raises -> PythonObject:
        """The matches of `query` (or its reverse complement) under `mode`, as `uint32` bytes, four per
        match: reference sequence, position in it, query position, length."""
        ref reference = self_ptr[].text
        var spans = List[Span[UInt8, ImmutAnyOrigin]]()
        spans.append(bases(query))
        if len(spans[0]) >= Int(UInt32.MAX):
            raise Error("a query holds fewer than 2^32 bases")
        var text = SeqText(spans, reverse_complement=Bool(py=reverse))
        var mems = find_mems(self_ptr[].index, reference, text, threads_of(threads))
        var kept: List[Mem]
        match Int(py=mode):
            case == MODE_MUM:
                kept = filter_mums(mems, text)
            case == MODE_MUMREFERENCE:
                kept = filter_mums(mems, text, unique_in_query=False)
            case _:
                kept = mems^
        var fields = List[UInt32](capacity=4 * len(kept))
        for m in kept:
            var seq = reference.sequence_of(Int(m.ref_pos))
            fields.append(UInt32(seq))
            fields.append(m.ref_pos - UInt32(reference.starts[seq]))
            fields.append(m.query_pos)
            fields.append(m.length)
        return bytes_of(fields)
