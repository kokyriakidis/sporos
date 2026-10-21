"""Frequency-capped (and floored) MEMs: MUMs relaxed to rare repeats.

A MEM's string occurs in the reference once per MEM whose query interval
contains its own, and in the query once per MEM whose reference interval
does; the MUM sweeps count those containers, up to a cap, instead of marking
them.
"""

from std.memory import unsafe_memmove

from sporos.parallel import ParallelTask, parallel_for
from sporos.sort import parallel_sort
from sporos.text import SeqText

from .engine import Mem
from .mums import _GroupFillTask, _PARALLEL_MEMS, _REF_ORDER_TASKS, _RefOrderTask, _collect, filter_mums


def filter_by_occurrences(
    var mems: List[Mem],
    query: SeqText,
    max_in_reference: Int,
    max_in_query: Int = 0,
    min_in_reference: Int = 1,
    threads: Int = 1,
) -> List[Mem]:
    """MEMs whose string occurs at most `max_in_reference` times in the
    reference, at most `max_in_query` times in its query sequence (0: no
    limit, for either) and at least `min_in_reference` times in the
    reference: frequency-capped (and floored) MEMs, as Mumemto's `-F` caps
    multi-MEMs. With caps of 1 they are `filter_mums`'s MUMs.

    The string of a MEM occurs in the reference once per MEM whose query
    interval contains its own (itself included), and in the query once per
    MEM whose reference interval does: sweeps that count containing
    intervals, up to the cap, give both. `mems` must come from one
    `find_mems` call (it is consumed: without limits it is returned as is);
    the result keeps its order."""
    if min_in_reference <= 1 and max_in_reference <= 1 and max_in_query <= 1:
        if max_in_reference == 0 and max_in_query == 0:
            return mems^
        if max_in_reference == 1:
            return filter_mums(mems, query, unique_in_query=max_in_query == 1, threads=threads)
    var keep = List[Bool](length=len(mems), fill=True)
    var workers = threads if len(mems) >= _PARALLEL_MEMS else 1
    var ref_cap = max(max_in_reference, min_in_reference)
    var lo = 0
    while lo < len(mems):
        var seq = query.sequence_of(Int(mems[lo].query_pos))
        var seq_end = query.starts[seq] + query.lengths[seq]
        var hi = lo + 1  # (a MEM past its sequence's end, from no find_mems, still moves on)
        while hi < len(mems) and Int(mems[hi].query_pos) < seq_end:
            hi += 1
        var count = hi - lo
        var starts = List[UInt64](unsafe_uninit_length=count)
        var lengths = List[UInt32](unsafe_uninit_length=count)
        var others = List[UInt32](unsafe_uninit_length=count)
        if max_in_reference > 0 or min_in_reference > 1:
            # Copies in the reference: MEMs containing it in query coordinates.
            var fill = _GroupFillTask(
                mems.unsafe_ptr().unsafe_offset(lo), starts.unsafe_ptr(), lengths.unsafe_ptr(), count, by_ref=False
            )
            parallel_for(fill, _REF_ORDER_TASKS, workers)
            _count_containing(starts, 0, lengths, max(ref_cap, 1), others, workers)
            for i in range(count):
                var copies = Int(others[i]) + 1
                if (max_in_reference > 0 and copies > max_in_reference) or copies < min_in_reference:
                    keep[lo + i] = False
        if max_in_query > 0:
            # Copies in the query: MEMs containing it in reference coordinates.
            var keys = _GroupFillTask(
                mems.unsafe_ptr().unsafe_offset(lo), starts.unsafe_ptr(), lengths.unsafe_ptr(), count, by_ref=True
            )
            parallel_for(keys, _REF_ORDER_TASKS, workers)
            parallel_sort(starts, workers, low_bits=32)
            var gather = _RefOrderTask(
                mems.unsafe_ptr().unsafe_offset(lo), starts.unsafe_ptr().as_imm(), lengths.unsafe_ptr(), count
            )
            parallel_for(gather, _REF_ORDER_TASKS, workers)
            _count_containing(starts, 32, lengths, max_in_query, others, workers)
            for i in range(count):
                if Int(others[i]) + 1 > max_in_query:
                    keep[lo + Int(starts[i] & 0xFFFFFFFF)] = False
        lo = hi
    return _collect(mems, keep, workers)


struct _TopEnds(Copyable, Movable):
    """The `cap` largest interval ends seen so far, ascending in
    `ends[head:]`."""

    var ends: List[Int]
    var head: Int
    var cap: Int

    def __init__(out self, cap: Int):
        self.ends = List[Int](capacity=min(cap, 64) + 16)
        self.head = 0
        self.cap = cap

    @inline(.always)
    def drop_through(mut self, start: Int):
        """Forgets the ends at or before `start`: they reach no interval
        starting there or later."""
        var p = self.ends.unsafe_ptr()
        while self.head < len(self.ends) and p[unsafe_offset=self.head] <= start:
            self.head += 1
        self._compact()

    @inline(.always)
    def _compact(mut self):
        """Drops the forgotten prefix once it outgrows the kept ends:
        amortized O(1) per end, and memory bounded by twice the kept ends."""
        var live = len(self.ends) - self.head
        if self.head >= 64 and self.head >= live:
            # (memmove allows the overlap; one untracked view of the buffer.)
            var p = self.ends.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
            unsafe_memmove(dest=p, src=p.unsafe_offset(self.head).as_imm(), count=live)
            self.ends.resize(unsafe_uninit_length=live)
            self.head = 0

    @inline(.always)
    def _first_at_least(self, end: Int) -> Int:
        """Index of the first kept end that reaches `end` (or the length)."""
        var lo = self.head
        var hi = len(self.ends)
        var p = self.ends.unsafe_ptr()
        while lo < hi:
            var mid = (lo + hi) >> 1
            if p[unsafe_offset=mid] < end:
                lo = mid + 1
            else:
                hi = mid
        return lo

    @inline(.always)
    def at_least(self, end: Int) -> Int:
        """How many of the kept ends reach `end`: the number of seen ends
        that do, capped at `cap`."""
        var hi = len(self.ends)
        # Most intervals are contained in none: no end reaches theirs.
        if hi == self.head or self.ends.unsafe_ptr()[unsafe_offset=hi - 1] < end:
            return 0
        return hi - self._first_at_least(end)

    @inline(.always)
    def add(mut self, end: Int):
        if len(self.ends) - self.head == self.cap:
            if end <= self.ends.unsafe_ptr()[unsafe_offset=self.head]:
                return
            self.head += 1
            self._compact()
        var n = len(self.ends)
        # Sweeps meet ends mostly in increasing order: usually an append;
        # otherwise one block move above the insertion point.
        if n == self.head or self.ends.unsafe_ptr()[unsafe_offset=n - 1] <= end:
            self.ends.append(end)
            return
        var at = self._first_at_least(end + 1)
        self.ends.append(end)
        var p = self.ends.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        unsafe_memmove(dest=p.unsafe_offset(at + 1), src=p.unsafe_offset(at).as_imm(), count=n - at)
        p[unsafe_offset=at] = end


def _count_containing(
    starts: List[UInt64], shift: UInt64, lengths: List[UInt32], var cap: Int, mut others: List[UInt32], threads: Int
):
    """`others[i]`: how many other intervals `[starts[j] >> shift, ... +
    lengths[j])`, in increasing start, contain (or equal) interval `i`,
    capped at `cap`. Parallel by blocks as `_mark_contained`: each block's
    sweep starts from the largest ends of the blocks before it."""
    var count = len(starts)
    # No interval has more than count - 1 others.
    cap = max(1, min(cap, count))
    if threads <= 1 or count < 1 << 16:
        var top = _TopEnds(cap)
        _count_sweep(Span(starts), shift, Span(lengths), others.unsafe_ptr(), 0, count, top, cap, True)
        return
    var blocks = threads * 4
    var bounds = List[Int](length=blocks + 1, fill=count)
    bounds[0] = 0
    for b in range(1, blocks):
        var at = max(bounds[b - 1], b * count // blocks)
        while at < count and at > 0 and starts[at] >> shift == starts[at - 1] >> shift:
            at += 1
        bounds[b] = at
    # Each block's own largest ends, then what precedes each block.
    var tops = List[_TopEnds](length=blocks, fill=_TopEnds(cap))
    var own = _CountTask(
        starts.unsafe_ptr().as_imm(),
        shift,
        lengths.unsafe_ptr().as_imm(),
        others.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        bounds.unsafe_ptr().as_imm(),
        tops.unsafe_ptr(),
        count,
        cap,
        write=False,
    )
    parallel_for(own, blocks, threads)
    var running = _TopEnds(cap)
    for b in range(blocks):
        var mine = tops[b].copy()
        if bounds[b] < count:
            running.drop_through(Int(starts[bounds[b]] >> shift) - 1)
        tops[b] = running.copy()
        for i in range(mine.head, len(mine.ends)):
            running.add(mine.ends[i])
    var sweep = _CountTask(
        starts.unsafe_ptr().as_imm(),
        shift,
        lengths.unsafe_ptr().as_imm(),
        others.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        bounds.unsafe_ptr().as_imm(),
        tops.unsafe_ptr(),
        count,
        cap,
        write=True,
    )
    parallel_for(sweep, blocks, threads)
    _ = len(tops)


def _count_sweep(
    starts: Span[UInt64, _],
    shift: UInt64,
    lengths: Span[UInt32, _],
    others: MutPointer[UInt32, _],
    lo: Int,
    hi: Int,
    mut top: _TopEnds,
    cap: Int,
    write: Bool,
):
    """Sweeps groups `[lo, hi)` (whole groups) after the ends in `top`,
    adding their ends to it; with `write`, sets `others`. Of a group sharing
    a start, every member at least as long as an interval contains it."""
    var group = List[UInt64]()
    var i = lo
    while i < hi:
        var start = Int(starts[i] >> shift)
        var j = i + 1
        while j < hi and Int(starts[j] >> shift) == start:
            j += 1
        # Ends at or before this start reach none of the intervals from here on.
        top.drop_through(start)
        if write:
            if j - i == 1:
                others[unsafe_offset=i] = UInt32(top.at_least(start + Int(lengths[i])))
            else:
                # Longest first: a member's copies within the group are the
                # members at least as long, ranked before it or tied with it.
                group.clear()
                for t in range(i, j):
                    group.append((UInt64(lengths[t]) << 32) | UInt64(t - i))
                sort(Span(group))
                var at = len(group)
                while at > 0:
                    var length = group[at - 1] >> 32
                    var run = at - 1
                    while run > 0 and group[run - 1] >> 32 == length:
                        run -= 1
                    # Members [run, len(group)) are at least this long.
                    var copies = min(cap, top.at_least(start + Int(length)) + len(group) - run - 1)
                    for u in range(run, at):
                        others[unsafe_offset=i + Int(group[u] & 0xFFFFFFFF)] = UInt32(copies)
                    at = run
        for t in range(i, j):
            top.add(start + Int(lengths[t]))
        i = j


struct _CountTask[
    starts_origin: ImmOrigin, lengths_origin: ImmOrigin, bounds_origin: ImmOrigin, tops_origin: MutOrigin
](ParallelTask):
    """One block of a parallel `_count_containing`: its own largest ends, or
    (`write`) its counts after the ends before it, in `tops`."""

    var starts: ImmPointer[UInt64, Self.starts_origin]
    var shift: UInt64
    var lengths: ImmPointer[UInt32, Self.lengths_origin]
    var others: MutPointer[UInt32, MutUntrackedOrigin]
    var bounds: ImmPointer[Int, Self.bounds_origin]
    var tops: MutPointer[_TopEnds, Self.tops_origin]
    var count: Int
    var cap: Int
    var write: Bool

    def __init__(
        out self,
        starts: ImmPointer[UInt64, Self.starts_origin],
        shift: UInt64,
        lengths: ImmPointer[UInt32, Self.lengths_origin],
        others: MutPointer[UInt32, MutUntrackedOrigin],
        bounds: ImmPointer[Int, Self.bounds_origin],
        tops: MutPointer[_TopEnds, Self.tops_origin],
        count: Int,
        cap: Int,
        *,
        write: Bool,
    ):
        self.starts = starts
        self.shift = shift
        self.lengths = lengths
        self.others = others
        self.bounds = bounds
        self.tops = tops
        self.count = count
        self.cap = cap
        self.write = write

    def run(self, index: Int):
        var lo = self.bounds[unsafe_offset=index]
        var hi = self.bounds[unsafe_offset=index + 1]
        var starts = Span(unsafe_ptr=self.starts, length=self.count)
        var lengths = Span(unsafe_ptr=self.lengths, length=self.count)
        var top = self.tops[unsafe_offset=index].copy()
        _count_sweep(starts, self.shift, lengths, self.others, lo, hi, top, self.cap, self.write)
        if not self.write:
            self.tops[unsafe_offset=index] = top^
