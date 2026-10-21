"""SMEMs from MEMs: the query intervals of MEMs (either strand) that no
other MEM's query interval contains. Needs no index: suited to one long query
such as a genome (reads go to the FM-index, `sporos.fmindex`)."""

from sporos.sort import parallel_sort
from sporos.text import SeqText

from .engine import Mem, MemIndex, find_mems


@fieldwise_init
struct Smem(Comparable, ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """A supermaximal exact match: query bases `[start, end)` (forward-strand
    text offsets) occur in the reference at `Smems.hits[first : first + count]`."""

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
    """Where an SMEM occurs: the reference text offset of its leftmost base,
    and whether the reference holds the query's reverse complement there."""

    var ref_pos: UInt32
    var reverse: Bool

    def __lt__(self, other: Self) -> Bool:
        return self.ref_pos < other.ref_pos or (self.ref_pos == other.ref_pos and not self.reverse and other.reverse)

    def __eq__(self, other: Self) -> Bool:
        return self.ref_pos == other.ref_pos and self.reverse == other.reverse


struct Smems(Movable, Sized):
    """SMEMs in increasing query start, each with its hits sorted by reference
    position (forward before reverse at the same position)."""

    var smems: List[Smem]
    var hits: List[SmemHit]

    def __init__(out self, var smems: List[Smem], var hits: List[SmemHit]):
        self.smems = smems^
        self.hits = hits^

    def __len__(self) -> Int:
        return len(self.smems)

    def hits_of(ref self, i: Int) -> Span[SmemHit, origin_of(self.hits)]:
        var smem = self.smems[i]
        return Span(self.hits)[smem.first : smem.first + smem.count]


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
    var smems = List[Smem]()
    var hits = List[SmemHit]()
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
                    hits.append(SmemHit(forward[i].ref_pos, False))
            for i in range(b, b_end):
                var m = backward[Int(keys[i] & 0xFFFFFFFF)]
                if Int(m.length) == longest:
                    hits.append(SmemHit(m.ref_pos, True))
            var count = len(hits) - first
            sort(Span(hits)[first:])
            smems.append(Smem(UInt32(start), UInt32(start + longest), first, count))
            reach = start + longest
        f = f_end
        b = b_end
    return Smems(smems^, hits^)
