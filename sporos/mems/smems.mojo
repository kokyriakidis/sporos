"""SMEMs from MEMs: the query intervals of MEMs (either strand) that no
other MEM's query interval contains. Needs no index: suited to one long query
such as a genome (reads go to the FM-index, `sporos.fmindex`)."""

from sporos.sort import parallel_sort
from sporos.text import SeqText

from .engine import Mem, MemIndex, find_mems


@fieldwise_init
struct Smem(Comparable, ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """A supermaximal exact match: bases `[start, end)` of its query occur
    `count` times in the reference (both strands); its listed hits start at
    `Smems.hits[first]` (see `Smems.hits_of`)."""

    var start: UInt32
    var end: UInt32
    var first: Int
    var count: Int

    def __lt__(self, other: Self) -> Bool:
        return self.start < other.start or (self.start == other.start and self.end < other.end)

    def __eq__(self, other: Self) -> Bool:
        return self.start == other.start and self.end == other.end and self.count == other.count


@fieldwise_init
struct SmemHit(Comparable, ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """Where an SMEM occurs: reference sequence `sequence`, 0-based
    `position` of the match's leftmost base in it, and whether the reference
    holds the query's reverse complement there."""

    var sequence: UInt32
    var position: UInt32
    var reverse: Bool

    def __lt__(self, other: Self) -> Bool:
        if self.sequence != other.sequence:
            return self.sequence < other.sequence
        if self.position != other.position:
            return self.position < other.position
        return not self.reverse and other.reverse

    def __eq__(self, other: Self) -> Bool:
        return self.sequence == other.sequence and self.position == other.position and self.reverse == other.reverse


struct Smems(Movable, Sized):
    """SMEMs of a batch of queries: query `q`'s are `of(q)`, in increasing
    start. Each SMEM's hits, sorted by reference sequence and position
    (forward before reverse), are `hits_of(i)`; an engine may leave them out
    for SMEMs with more than its `max_hits`, while `Smem.count` always
    counts them all."""

    var offsets: List[Int]
    var smems: List[Smem]
    var hits: List[SmemHit]

    def __init__(out self):
        self.offsets = [0]
        self.smems = List[Smem]()
        self.hits = List[SmemHit]()

    def __len__(self) -> Int:
        """The number of queries."""
        return len(self.offsets) - 1

    def of(ref self, query: Int) -> Span[Smem, origin_of(self.smems)]:
        """The SMEMs of query `query`."""
        return Span(self.smems)[self.offsets[query] : self.offsets[query + 1]]

    def hits_of(ref self, i: Int) -> Span[SmemHit, origin_of(self.hits)]:
        """The listed hits of SMEM `i` (an index into `smems`)."""
        var end = self.smems[i + 1].first if i + 1 < len(self.smems) else len(self.hits)
        return Span(self.hits)[self.smems[i].first : end]


def find_smems(index: MemIndex, reference: SeqText, query: SeqText, query_rc: SeqText, threads: Int = 1) -> Smems:
    """The SMEMs (supermaximal exact matches, as in BWA) of length >=
    `index.params.min_len` of every query sequence against the reference on
    both strands. `query_rc` packs the same sequences reverse-complemented
    (`SeqText(seqs, reverse_complement=True)`).

    An SMEM is the query interval of a MEM (on either strand) that no other
    MEM's query interval strictly contains, with every reference position
    where that whole interval matches. Any interval containing a MEM of
    length >= L is longer, so filtering at L gives exactly the SMEMs of
    length >= L."""
    var forward = find_mems(index, reference, query, threads)
    var backward = find_mems(index, reference, query_rc, threads)
    # Reverse-strand MEMs in forward query coordinates, sorted by start.
    var keys = List[UInt64](capacity=len(backward))
    var seq = 0
    for i in range(len(backward)):
        var m = backward[i]
        var q = Int(m.query_pos)
        while seq + 1 < len(query_rc.starts) and query_rc.starts[seq + 1] <= q:
            seq += 1
        var base = query_rc.starts[seq]
        var start = base + query_rc.lengths[seq] - (q - base) - Int(m.length)
        keys.append((UInt64(start) << 32) | UInt64(i))
    parallel_sort(keys, threads, low_bits=32)
    # Sweep both lists by query start: of the MEMs starting at one position,
    # only the longest can be an SMEM, and it is one unless an earlier start
    # reaches as far.
    var result = Smems()
    var counts = List[Int](length=len(query.starts) + 1, fill=0)
    ref smems = result.smems
    ref hits = result.hits
    var reach = -1
    var f = 0
    var b = 0
    while f < len(forward) or b < len(keys):
        var start: Int
        if b == len(keys) or (f < len(forward) and Int(forward[f].query_pos) <= Int(keys[b] >> 32)):
            start = Int(forward[f].query_pos)
        else:
            start = Int(keys[b] >> 32)
        var longest = 0
        var f_end = f
        while f_end < len(forward) and Int(forward[f_end].query_pos) == start:
            longest = max(longest, Int(forward[f_end].length))
            f_end += 1
        var b_end = b
        while b_end < len(keys) and Int(keys[b_end] >> 32) == start:
            longest = max(longest, Int(backward[Int(keys[b_end] & 0xFFFFFFFF)].length))
            b_end += 1
        if start + longest > reach:
            var first = len(hits)
            for i in range(f, f_end):
                if Int(forward[i].length) == longest:
                    hits.append(_hit(reference, Int(forward[i].ref_pos), False))
            for i in range(b, b_end):
                var m = backward[Int(keys[i] & 0xFFFFFFFF)]
                if Int(m.length) == longest:
                    hits.append(_hit(reference, Int(m.ref_pos), True))
            var count = len(hits) - first
            sort(Span(hits)[first:])
            # In its own query's coordinates (MEMs never cross a sequence end).
            var q = query.sequence_of(start)
            var base = query.starts[q]
            smems.append(Smem(UInt32(start - base), UInt32(start - base + longest), first, count))
            counts[q + 1] += 1
            reach = start + longest
        f = f_end
        b = b_end
    for q in range(len(query.starts)):
        counts[q + 1] += counts[q]
    result.offsets = counts^
    return result^


@always_inline
def _hit(reference: SeqText, offset: Int, reverse: Bool) -> SmemHit:
    var seq = reference.sequence_of(offset)
    return SmemHit(UInt32(seq), UInt32(offset - reference.starts[seq]), reverse)
