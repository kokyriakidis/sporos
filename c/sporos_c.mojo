"""
The C API of sporos: (w, k)-minimizers, closed and open syncmers, and the MEMs and MUMs of a query
against an indexed reference, for C, C++ and any language with a C foreign-function interface.
`sporos.h` declares these.

    pixi run build-c [target-cpu]   # build/c: libsporos, its runtime libraries and the header

Sequences come as ASCII bytes and lengths. Seed positions go into the caller's buffer; an index and
a set of matches are opaque handles the caller frees. Failure is a negative result or a null handle
(`SPOROS_*` in the header). Threads are sporos's own pthreads, so no Mojo runtime is set up.
"""

from std.memory import Allocation, Layout, alloc, dealloc, unsafe_memcpy

from sporos.mems import Mem, MemIndex, MemParams, filter_mums, find_mems
from sporos.minimizers import minimizer_positions
from sporos.parallel import default_threads
from sporos.syncmers import syncmer_positions
from sporos.text import SeqText

comptime INVALID_ARGUMENT = -1
"""`SPOROS_INVALID_ARGUMENT`: a negative length, or parameters the kernel refuses."""

comptime MODE_MAXMATCH = 0
"""`SPOROS_MAXMATCH`: every MEM (MUMmer's `-maxmatch`)."""
comptime MODE_MUM = 1
"""`SPOROS_MUM`: MEMs unique in the reference and the query (MUMmer's `-mum`)."""
comptime MODE_MUMREFERENCE = 2
"""`SPOROS_MUMREFERENCE`: MEMs unique in the reference (MUMmer's `-mumreference`)."""

comptime CBytes = ImmPointer[UInt8, ImmutAnyOrigin]
"""A C `const uint8_t*`."""


@fieldwise_init
struct Match(ImplicitlyCopyable, TrivialRegisterPassable):
    """`sporos_match`: four `uint32_t`, laid out as the C struct."""

    var ref_seq: UInt32
    var ref_pos: UInt32
    var query_pos: UInt32
    var length: UInt32


struct Index(Movable):
    """`sporos_index`: the packed reference and its anchor index, which reads it."""

    var text: SeqText
    var index: MemIndex

    def __init__(out self, var text: SeqText, min_len: Int, threads: Int) raises:
        var params = MemParams.choose(min_len, len(text), len(text))
        self.index = MemIndex(text, params, threads)
        self.text = text^


struct Matches(Movable):
    """`sporos_matches`: one query's matches, sorted by query position."""

    var matches: List[Match]

    def __init__(out self, var matches: List[Match]):
        self.matches = matches^


@inline(.always)
def thread_count(threads: Int32) -> Int:
    """`threads`, or every physical core for zero or less."""
    return Int(threads) if threads > 0 else default_threads()


def bytes_of(seq: CBytes, length: Int) -> Span[UInt8, ImmutAnyOrigin]:
    """A C buffer as a span; never read when empty, so it may be null then."""
    return Span(unsafe_ptr=seq, length=length)


def copy_out(positions: List[UInt32], dest: MutPointer[UInt32, MutAnyOrigin], capacity: Int) -> Int:
    """Writes the first `capacity` positions to `dest`; returns how many there are."""
    unsafe_memcpy(dest=dest, src=positions.unsafe_ptr(), count=min(len(positions), capacity))
    return len(positions)


@export("sporos_minimizers")
def sporos_minimizers(
    seq: CBytes,
    length: Int,
    k: Int32,
    w: Int32,
    canonical: Int32,
    threads: Int32,
    positions_out: MutPointer[UInt32, MutAnyOrigin],
    out_capacity: Int,
) abi("C") -> Int:
    """The (w, k)-minimizer positions of `seq` into `positions_out`, at most `out_capacity` of them; returns how
    many there are, or `INVALID_ARGUMENT`."""
    if length < 0 or out_capacity < 0:
        return INVALID_ARGUMENT
    try:
        var positions = minimizer_positions(
            bytes_of(seq, length), Int(k), Int(w), canonical != 0, thread_count(threads)
        )
        return copy_out(positions, positions_out, out_capacity)
    except:
        return INVALID_ARGUMENT


@export("sporos_syncmers")
def sporos_syncmers(
    seq: CBytes,
    length: Int,
    k: Int32,
    s: Int32,
    open_offset: Int32,
    canonical: Int32,
    threads: Int32,
    positions_out: MutPointer[UInt32, MutAnyOrigin],
    out_capacity: Int,
) abi("C") -> Int:
    """The syncmer positions of `seq` into `positions_out`: closed for a negative `open_offset`, else open at
    that offset; returns how many there are, or `INVALID_ARGUMENT`."""
    if length < 0 or out_capacity < 0:
        return INVALID_ARGUMENT
    try:
        var offset = Optional[Int](Int(open_offset)) if open_offset >= 0 else None
        var positions = syncmer_positions(
            bytes_of(seq, length), Int(k), Int(s), canonical != 0, offset, thread_count(threads)
        )
        return copy_out(positions, positions_out, out_capacity)
    except:
        return INVALID_ARGUMENT


@export("sporos_syncmers_open")
def sporos_syncmers_open(
    seq: CBytes,
    length: Int,
    k: Int32,
    s: Int32,
    offsets: ImmPointer[Int32, ImmutAnyOrigin],
    offset_count: Int,
    canonical: Int32,
    threads: Int32,
    positions_out: MutPointer[UInt32, MutAnyOrigin],
    out_capacity: Int,
) abi("C") -> Int:
    """The open syncmer positions of `seq` for any of `offset_count` s-mer offsets into `positions_out`;
    returns how many there are, or `INVALID_ARGUMENT`."""
    if length < 0 or out_capacity < 0 or offset_count < 1:
        return INVALID_ARGUMENT
    try:
        var opened = List[Int](capacity=offset_count)
        for i in range(offset_count):
            opened.append(Int(offsets[unsafe_offset=i]))
        var positions = syncmer_positions(
            bytes_of(seq, length), Int(k), Int(s), canonical != 0, None, thread_count(threads), open_offsets=opened
        )
        return copy_out(positions, positions_out, out_capacity)
    except:
        return INVALID_ARGUMENT


@export("sporos_index_new")
def sporos_index_new(
    seqs: ImmPointer[CBytes, ImmutAnyOrigin],
    lengths: ImmPointer[Int, ImmutAnyOrigin],
    count: Int,
    min_len: Int32,
    threads: Int32,
) abi("C") -> OptionalPointer[Index, MutAnyOrigin]:
    """An index of `count` reference sequences for matches of at least `min_len`, for
    `sporos_index_free`; null for bad arguments or a reference of 2^32 bases or more."""
    if count < 1:
        return None
    var spans = List[Span[UInt8, ImmutAnyOrigin]](capacity=count)
    for i in range(count):
        var length = lengths[unsafe_offset=i]
        if length < 0:
            return None
        spans.append(bytes_of(seqs[unsafe_offset=i], length))
    try:
        var index = Index(SeqText(spans), Int(min_len), thread_count(threads))
        var place = alloc(Layout[Index].single())
        place.unsafe_ptr().unsafe_write(index^)
        return place^.unsafe_leak().unsafe_origin_cast[MutAnyOrigin]()
    except:
        return None


@export("sporos_index_free")
def sporos_index_free(index: OptionalPointer[Index, MutAnyOrigin]) abi("C"):
    """Frees an index `sporos_index_new` made; null is nothing to free."""
    if case .Some(at) = index:
        var owned = Allocation(
            unsafe_owned_ptr=at.unsafe_origin_cast[MutUntrackedOrigin](), layout=Layout[Index].single()
        )
        _ = owned.unsafe_ptr().unsafe_take_pointee()
        dealloc(owned^)


@export("sporos_index_sequences")
def sporos_index_sequences(index: ImmPointer[Index, ImmutAnyOrigin]) abi("C") -> Int:
    """How many reference sequences `index` holds."""
    return len(index[].text.starts)


@export("sporos_index_matches")
def sporos_index_matches(
    index: ImmPointer[Index, ImmutAnyOrigin],
    query: CBytes,
    length: Int,
    mode: Int32,
    reverse_complement: Int32,
    threads: Int32,
) abi("C") -> OptionalPointer[Matches, MutAnyOrigin]:
    """The matches of `query` (its reverse complement when `reverse_complement` is nonzero, positions
    then counted on it) against `index` under `mode`, for `sporos_matches_free`; null for bad
    arguments. The index is only read, so many threads may query it at once."""
    if length < 0 or mode < MODE_MAXMATCH or mode > MODE_MUMREFERENCE:
        return None
    # Query offsets are 32-bit.
    if length >= Int(UInt32.MAX):
        return None
    var spans = List[Span[UInt8, ImmutAnyOrigin]]()
    spans.append(bytes_of(query, length))
    var text = SeqText(spans, reverse_complement=reverse_complement != 0)
    ref reference = index[].text
    var mems = find_mems(index[].index, reference, text, thread_count(threads))
    var kept: List[Mem]
    match Int(mode):
        case == MODE_MUM:
            kept = filter_mums(mems, text)
        case == MODE_MUMREFERENCE:
            kept = filter_mums(mems, text, unique_in_query=False)
        case _:
            kept = mems^
    var matches = List[Match](capacity=len(kept))
    for m in kept:
        var seq = reference.sequence_of(Int(m.ref_pos))
        var start = UInt32(reference.starts[seq])
        matches.append(Match(UInt32(seq), m.ref_pos - start, m.query_pos, m.length))
    var place = alloc(Layout[Matches].single())
    place.unsafe_ptr().unsafe_write(Matches(matches^))
    return place^.unsafe_leak().unsafe_origin_cast[MutAnyOrigin]()


@export("sporos_matches_count")
def sporos_matches_count(matches: ImmPointer[Matches, ImmutAnyOrigin]) abi("C") -> Int:
    """How many matches there are."""
    return len(matches[].matches)


@export("sporos_matches_data")
def sporos_matches_data(matches: ImmPointer[Matches, ImmutAnyOrigin]) abi("C") -> ImmPointer[Match, ImmutAnyOrigin]:
    """The matches, `sporos_matches_count` of them, alive until `sporos_matches_free`."""
    return matches[].matches.unsafe_ptr()


@export("sporos_matches_free")
def sporos_matches_free(matches: OptionalPointer[Matches, MutAnyOrigin]) abi("C"):
    """Frees matches `sporos_index_matches` made; null is nothing to free."""
    if case .Some(at) = matches:
        var owned = Allocation(
            unsafe_owned_ptr=at.unsafe_origin_cast[MutUntrackedOrigin](), layout=Layout[Matches].single()
        )
        _ = owned.unsafe_ptr().unsafe_take_pointee()
        dealloc(owned^)
