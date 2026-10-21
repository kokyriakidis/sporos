"""
The C API of sporos: (w, k)-minimizers, closed and open syncmers, the MEMs, MUMs and frequency-capped
MEMs of a query against an indexed reference, the SMEMs of reads against an FM-index, and seed-index
lookups, for C, C++ and any language with a C foreign-function interface. `sporos.h` declares these.

    pixi run build-c [target-cpu]   # build/c: libsporos, its runtime libraries and the header

Sequences come as ASCII bytes and lengths. Seed positions go into the caller's buffer; indexes and
their answers are opaque handles the caller frees. Failure is a negative result or a null handle
(`SPOROS_*` in the header). Threads are sporos's own pthreads, so no Mojo runtime is set up.
"""

from std.memory import Allocation, Layout, alloc, dealloc, unsafe_memcpy

from sporos.fmindex import FMIndex
from sporos.index import ALL_KMERS, SeedIndex
from sporos.mems import MemIndex, MemParams, filter_by_occurrences, find_mems
from sporos.minimizers import minimizer_positions
from sporos.parallel import MAX_THREADS, default_threads
from sporos.syncmers import syncmer_positions
from sporos.text import SeqText

comptime INVALID_ARGUMENT = -1
"""`SPOROS_INVALID_ARGUMENT`: a negative length, or parameters the kernel refuses."""
comptime IO_ERROR = -2
"""`SPOROS_IO_ERROR`: a file that could not be written."""

comptime MODE_MAXMATCH = 0
"""`SPOROS_MAXMATCH`: every MEM (MUMmer's `-maxmatch`)."""
comptime MODE_MUM = 1
"""`SPOROS_MUM`: MEMs unique in the reference and the query (MUMmer's `-mum`)."""
comptime MODE_MUMREFERENCE = 2
"""`SPOROS_MUMREFERENCE`: MEMs unique in the reference (MUMmer's `-mumreference`)."""

comptime CBytes = ImmPointer[UInt8, ImmutAnyOrigin]
"""A C `const uint8_t*`."""
comptime CSeqs = ImmPointer[CBytes, ImmutAnyOrigin]
"""A C `const uint8_t *const *`: sequences, with their lengths beside."""
comptime CLengths = ImmPointer[Int, ImmutAnyOrigin]
"""A C `const int64_t *`: sequence lengths."""


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


@fieldwise_init
struct CSmem(ImplicitlyCopyable, TrivialRegisterPassable):
    """`sporos_smem`, laid out as the C struct."""

    var start: UInt32
    var end: UInt32
    var count: Int64
    var first_hit: Int64
    var hit_count: Int64


@fieldwise_init
struct CSmemHit(ImplicitlyCopyable, TrivialRegisterPassable):
    """`sporos_smem_hit`, laid out as the C struct."""

    var ref_seq: UInt32
    var ref_pos: UInt32
    var reverse: Int32


struct FM(Movable):
    """`sporos_fm_index`: the FM-index and its sequence names as C strings."""

    var index: FMIndex
    var names: List[List[UInt8]]

    def __init__(out self, var index: FMIndex):
        self.names = List[List[UInt8]](capacity=len(index.names))
        for name in index.names:
            var bytes = List[UInt8](name.as_bytes())
            bytes.append(0)
            self.names.append(bytes^)
        self.index = index^


struct CSmems(Movable):
    """`sporos_smems`: every query's SMEMs, `smems[offsets[q]:offsets[q + 1]]`, and their hits."""

    var offsets: List[Int]
    var smems: List[CSmem]
    var hits: List[CSmemHit]

    def __init__(out self, var offsets: List[Int], var smems: List[CSmem], var hits: List[CSmemHit]):
        self.offsets = offsets^
        self.smems = smems^
        self.hits = hits^


@fieldwise_init
struct CSeedHit(ImplicitlyCopyable, TrivialRegisterPassable):
    """`sporos_seed_hit`, laid out as the C struct."""

    var query_pos: UInt32
    var ref_seq: UInt32
    var ref_pos: UInt32
    var occurrences: UInt32
    var reverse: Int32


struct CSeedHits(Movable):
    """`sporos_seed_hits`: every read's hits, `hits[offsets[r]:offsets[r + 1]]`."""

    var offsets: List[Int]
    var hits: List[CSeedHit]

    def __init__(out self, var offsets: List[Int], var hits: List[CSeedHit]):
        self.offsets = offsets^
        self.hits = hits^


@inline(.always)
def thread_count(threads: Int32) -> Int:
    """`threads`, or every physical core for zero or less."""
    return min(Int(threads), MAX_THREADS) if threads > 0 else default_threads()


def bytes_of(seq: CBytes, length: Int) -> Span[UInt8, ImmutAnyOrigin]:
    """A C buffer as a span; never read when empty, so it may be null then."""
    return Span(unsafe_ptr=seq, length=length)


def spans_of(seqs: CSeqs, lengths: CLengths, count: Int) -> Optional[List[Span[UInt8, ImmutAnyOrigin]]]:
    """`count` C buffers as spans; none for a negative length."""
    var spans = List[Span[UInt8, ImmutAnyOrigin]](capacity=count)
    for i in range(count):
        var length = lengths[unsafe_offset=i]
        if length < 0:
            return None
        spans.append(bytes_of(seqs[unsafe_offset=i], length))
    return spans^


def string_of(text: CBytes) -> String:
    """A NUL-terminated C string."""
    var length = 0
    while text[unsafe_offset=length] != 0:
        length += 1
    return String(StringSlice(unsafe_from_utf8=bytes_of(text, length)))


def handle[T: Movable](var value: T) -> OptionalPointer[T, MutAnyOrigin]:
    """`value` on the heap, for the caller to free with `free_handle`."""
    var place = alloc(Layout[T].single())
    place.unsafe_ptr().unsafe_write(value^)
    return place^.unsafe_leak().unsafe_origin_cast[MutAnyOrigin]()


def free_handle[T: Movable & Deinitable](at: OptionalPointer[T, MutAnyOrigin]):
    """Frees a `handle`; null is nothing to free."""
    if case .Some(p) = at:
        var owned = Allocation(unsafe_owned_ptr=p.unsafe_origin_cast[MutUntrackedOrigin](), layout=Layout[T].single())
        _ = owned.unsafe_ptr().unsafe_take_pointee()
        dealloc(owned^)


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
    seqs: CSeqs, lengths: CLengths, count: Int, min_len: Int32, threads: Int32
) abi("C") -> OptionalPointer[Index, MutAnyOrigin]:
    """An index of `count` reference sequences for matches of at least `min_len`, for
    `sporos_index_free`; null for bad arguments or a reference of 2^32 bases or more."""
    if count < 1:
        return None
    var spans = spans_of(seqs, lengths, count)
    if not spans:
        return None
    try:
        return handle(Index(SeqText(spans.value()), Int(min_len), thread_count(threads)))
    except:
        return None


@export("sporos_index_free")
def sporos_index_free(index: OptionalPointer[Index, MutAnyOrigin]) abi("C"):
    """Frees an index `sporos_index_new` made; null is nothing to free."""
    free_handle(index)


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
    max_occurrences: Int32,
    min_occurrences: Int32,
    threads: Int32,
) abi("C") -> OptionalPointer[Matches, MutAnyOrigin]:
    """The matches of `query` (its reverse complement when `reverse_complement` is nonzero, positions
    then counted on it) against `index` under `mode`, for `sporos_matches_free`; null for bad
    arguments. `max_occurrences` N > 0 relaxes the modes' "once" to "at most N times" (for maxmatch it
    caps the reference count) and `min_occurrences` asks for at least that many reference copies, as
    `Reference.mems` does. The index is only read, so many threads may query it at once."""
    if length < 0 or mode < MODE_MAXMATCH or mode > MODE_MUMREFERENCE or max_occurrences < 0:
        return None
    var cap = Int(max_occurrences) if max_occurrences > 0 else (0 if mode == MODE_MAXMATCH else 1)
    if min_occurrences < 1 or (cap > 0 and Int(min_occurrences) > cap):
        return None
    # Query offsets are 32-bit.
    if length >= Int(UInt32.MAX):
        return None
    var spans = List[Span[UInt8, ImmutAnyOrigin]]()
    spans.append(bytes_of(query, length))
    var text = SeqText(spans, reverse_complement=reverse_complement != 0)
    ref reference = index[].text
    var mems = find_mems(index[].index, reference, text, thread_count(threads))
    var kept = filter_by_occurrences(
        mems^, text, cap, cap if mode == MODE_MUM else 0, Int(min_occurrences), thread_count(threads)
    )
    var matches = List[Match](capacity=len(kept))
    for m in kept:
        var seq = reference.sequence_of(Int(m.ref_pos))
        var start = UInt32(reference.starts[seq])
        matches.append(Match(UInt32(seq), m.ref_pos - start, m.query_pos, m.length))
    return handle(Matches(matches^))


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
    free_handle(matches)


# ------------------------------------------------------------------ FM-index


@export("sporos_fm_index_new")
def sporos_fm_index_new(
    seqs: CSeqs, lengths: CLengths, names: OptionalPointer[CBytes, ImmutAnyOrigin], count: Int, threads: Int32
) abi("C") -> OptionalPointer[FM, MutAnyOrigin]:
    """An FM-index of `count` reference sequences, named by `names` (or `seq0`, `seq1`, ... when it
    is null), for `sporos_fm_index_free`; null for bad arguments."""
    if count < 1:
        return None
    var spans = spans_of(seqs, lengths, count)
    if not spans:
        return None
    var named = List[String](capacity=count)
    for i in range(count):
        if case .Some(given) = names:
            named.append(string_of(given[unsafe_offset=i]))
        else:
            named.append(String(t"seq{i}"))
    try:
        return handle(FM(FMIndex.build(spans.value(), named, thread_count(threads))))
    except:
        return None


@export("sporos_fm_index_load")
def sporos_fm_index_load(path: CBytes, threads: Int32) abi("C") -> OptionalPointer[FM, MutAnyOrigin]:
    """The FM-index saved at `path`, for `sporos_fm_index_free`; null when it cannot be read."""
    try:
        return handle(FM(FMIndex.load(string_of(path), thread_count(threads))))
    except:
        return None


@export("sporos_fm_index_save")
def sporos_fm_index_save(index: ImmPointer[FM, ImmutAnyOrigin], path: CBytes) abi("C") -> Int32:
    """Writes `index` to `path`: 0, or `IO_ERROR`."""
    try:
        index[].index.save(string_of(path))
        return 0
    except:
        return IO_ERROR


@export("sporos_fm_index_free")
def sporos_fm_index_free(index: OptionalPointer[FM, MutAnyOrigin]) abi("C"):
    """Frees an FM-index; null is nothing to free."""
    free_handle(index)


@export("sporos_fm_index_sequences")
def sporos_fm_index_sequences(index: ImmPointer[FM, ImmutAnyOrigin]) abi("C") -> Int:
    """How many reference sequences `index` holds."""
    return len(index[].names)


@export("sporos_fm_index_name")
def sporos_fm_index_name(
    index: ImmPointer[FM, ImmutAnyOrigin], i: Int
) abi("C") -> OptionalPointer[UInt8, ImmutAnyOrigin]:
    """Sequence `i`'s name, NUL-terminated, alive with the index; null for no such sequence."""
    if i < 0 or i >= len(index[].names):
        return None
    return index[].names[i].unsafe_ptr().unsafe_origin_cast[ImmutAnyOrigin]()


@export("sporos_fm_index_smems")
def sporos_fm_index_smems(
    index: ImmPointer[FM, ImmutAnyOrigin],
    queries: CSeqs,
    lengths: CLengths,
    count: Int,
    min_len: Int32,
    max_hits: Int32,
    min_occurrences: Int32,
    max_occurrences: Int32,
    threads: Int32,
) abi("C") -> OptionalPointer[CSmems, MutAnyOrigin]:
    """The SMEMs of `count` queries, both strands, for `sporos_smems_free` (see `FMIndex.smems`);
    null for bad arguments."""
    if count < 0 or min_len < 1 or max_hits < 0 or min_occurrences < 1 or max_occurrences < 0:
        return None
    var spans = spans_of(queries, lengths, count)
    if not spans:
        return None
    for span in spans.value():
        if len(span) >= Int(UInt32.MAX):
            return None
    var found = index[].index.smems(
        spans.value(), Int(min_len), Int(max_hits), Int(min_occurrences), Int(max_occurrences), thread_count(threads)
    )
    var smems = List[CSmem](capacity=len(found.smems))
    for i in range(len(found.smems)):
        var smem = found.smems[i]
        var hits = len(found.hits_of(i))
        smems.append(CSmem(smem.start, smem.end, Int64(smem.count), Int64(smem.first), Int64(hits)))
    var hits = List[CSmemHit](capacity=len(found.hits))
    for h in found.hits:
        hits.append(CSmemHit(h.sequence, h.position, Int32(h.reverse)))
    return handle(CSmems(found.offsets.copy(), smems^, hits^))


@export("sporos_smems_count")
def sporos_smems_count(smems: ImmPointer[CSmems, ImmutAnyOrigin], query: Int) abi("C") -> Int:
    """How many SMEMs query `query` has."""
    return smems[].offsets[query + 1] - smems[].offsets[query]


@export("sporos_smems_data")
def sporos_smems_data(
    smems: ImmPointer[CSmems, ImmutAnyOrigin], query: Int
) abi("C") -> ImmPointer[CSmem, ImmutAnyOrigin]:
    """Query `query`'s SMEMs, `sporos_smems_count` of them, alive until `sporos_smems_free`."""
    return smems[].smems.unsafe_ptr().unsafe_offset(smems[].offsets[query])


@export("sporos_smems_hits")
def sporos_smems_hits(smems: ImmPointer[CSmems, ImmutAnyOrigin]) abi("C") -> ImmPointer[CSmemHit, ImmutAnyOrigin]:
    """Every SMEM's hits; an SMEM's are `hit_count` of them from `first_hit`."""
    return smems[].hits.unsafe_ptr()


@export("sporos_smems_free")
def sporos_smems_free(smems: OptionalPointer[CSmems, MutAnyOrigin]) abi("C"):
    """Frees SMEMs; null is nothing to free."""
    free_handle(smems)


# ---------------------------------------------------------------- seed index


@export("sporos_seed_index_new")
def sporos_seed_index_new(
    seqs: CSeqs,
    lengths: CLengths,
    count: Int,
    k: Int32,
    w: Int32,
    scheme: Int32,
    open_offset: Int32,
    threads: Int32,
) abi("C") -> OptionalPointer[SeedIndex, MutAnyOrigin]:
    """A seed index of `count` reference sequences (see `SeedIndex`), for `sporos_seed_index_free`;
    null for bad arguments."""
    if count < 1 or scheme < 0 or scheme > ALL_KMERS:
        return None
    var spans = spans_of(seqs, lengths, count)
    if not spans:
        return None
    try:
        var offset = Optional[Int](Int(open_offset)) if open_offset >= 0 else None
        return handle(SeedIndex(spans.value(), Int(k), Int(w), Int(scheme), thread_count(threads), offset))
    except:
        return None


@export("sporos_seed_index_free")
def sporos_seed_index_free(index: OptionalPointer[SeedIndex, MutAnyOrigin]) abi("C"):
    """Frees a seed index; null is nothing to free."""
    free_handle(index)


@export("sporos_seed_index_occurrence_cutoff")
def sporos_seed_index_occurrence_cutoff(
    index: ImmPointer[SeedIndex, ImmutAnyOrigin], keep_fraction: Float64
) abi("C") -> Int:
    """The smallest occurrence count that `keep_fraction` of the indexed k-mers stay within."""
    return index[].occurrence_cutoff(keep_fraction)


@export("sporos_seed_index_lookup")
def sporos_seed_index_lookup(
    index: ImmPointer[SeedIndex, ImmutAnyOrigin],
    reads: CSeqs,
    lengths: CLengths,
    count: Int,
    all_kmers: Int32,
    max_occurrences: Int,
    threads: Int32,
) abi("C") -> OptionalPointer[CSeedHits, MutAnyOrigin]:
    """The hits of every read's seeds (every k-mer's with `all_kmers`), skipping k-mers occurring
    more than `max_occurrences` times (0: no limit), for `sporos_seed_hits_free`; null for bad
    arguments."""
    if count < 0 or max_occurrences < 0:
        return None
    var spans = spans_of(reads, lengths, count)
    if not spans:
        return None
    for span in spans.value():
        if len(span) >= Int(UInt32.MAX):
            return None
    var found = index[].lookup_batch(spans.value(), all_kmers != 0, max_occurrences, thread_count(threads))
    ref text = index[].text
    var hits = List[CSeedHit](capacity=len(found.hits))
    for h in found.hits:
        var seq = text.sequence_of(Int(h.ref_pos))
        hits.append(
            CSeedHit(h.query_pos, UInt32(seq), h.ref_pos - UInt32(text.starts[seq]), h.occurrences, Int32(h.reverse))
        )
    return handle(CSeedHits(found.offsets.copy(), hits^))


@export("sporos_seed_hits_count")
def sporos_seed_hits_count(hits: ImmPointer[CSeedHits, ImmutAnyOrigin], r: Int) abi("C") -> Int:
    """How many hits read `r` has."""
    return hits[].offsets[r + 1] - hits[].offsets[r]


@export("sporos_seed_hits_data")
def sporos_seed_hits_data(
    hits: ImmPointer[CSeedHits, ImmutAnyOrigin], r: Int
) abi("C") -> ImmPointer[CSeedHit, ImmutAnyOrigin]:
    """Read `r`'s hits in query order, `sporos_seed_hits_count` of them, alive until
    `sporos_seed_hits_free`."""
    return hits[].hits.unsafe_ptr().unsafe_offset(hits[].offsets[r])


@export("sporos_seed_hits_free")
def sporos_seed_hits_free(hits: OptionalPointer[CSeedHits, MutAnyOrigin]) abi("C"):
    """Frees seed hits; null is nothing to free."""
    free_handle(hits)
