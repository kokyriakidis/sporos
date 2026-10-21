"""Maximal unique matches (MUMs) among MEMs.

A MEM is a MUM when its string occurs exactly once in the reference and once
in its query sequence. Any other occurrence would be covered by a different
MEM (of length >= L) whose reference or query interval contains this one, so
MUMs are exactly the MEMs not contained, in either coordinate, by another MEM
of the same query sequence: linear sweeps, parallel by blocks.
"""

from std.math import ceildiv

from sporos.parallel import ParallelTask, parallel_for
from sporos.sort import parallel_sort
from sporos.text import SeqText

from .engine import Mem


def filter_mums(mems: List[Mem], query: SeqText, unique_in_query: Bool = True, threads: Int = 1) -> List[Mem]:
    """MEMs unique in the reference and (by default) in their query sequence.

    These are MUMmer's `-mum` matches; with `unique_in_query=False` they are
    its `-mumreference` matches (unique in the reference only). `mems` must
    come from one `find_mems` call (sorted by query offset); so is the result.
    """
    var keep = List[Bool](length=len(mems), fill=True)
    # Threads pay off only for many MEMs; per-read calls stay serial.
    var workers = threads if len(mems) >= _PARALLEL_MEMS else 1
    var lo = 0
    while lo < len(mems):
        var seq = query.sequence_of(Int(mems[lo].query_pos))
        var seq_end = query.starts[seq] + query.lengths[seq]
        var hi = lo + 1  # (a MEM past its sequence's end, from no find_mems, still moves on)
        while hi < len(mems) and Int(mems[hi].query_pos) < seq_end:
            hi += 1
        # Query intervals: `mems` is already in query order.
        var count = hi - lo
        var starts = List[UInt64](unsafe_uninit_length=count)
        var lengths = List[UInt32](unsafe_uninit_length=count)
        var fill = _GroupFillTask(
            mems.unsafe_ptr().unsafe_offset(lo), starts.unsafe_ptr(), lengths.unsafe_ptr(), count, by_ref=False
        )
        parallel_for(fill, _REF_ORDER_TASKS, workers)
        _mark_contained(starts, 0, lengths, keep.unsafe_ptr().unsafe_offset(lo), workers)
        if unique_in_query:
            # Reference intervals, in reference order: start << 32 | member.
            var keys = _GroupFillTask(
                mems.unsafe_ptr().unsafe_offset(lo), starts.unsafe_ptr(), lengths.unsafe_ptr(), count, by_ref=True
            )
            parallel_for(keys, _REF_ORDER_TASKS, workers)
            # Order by start only: the sweep groups equal starts in any order.
            parallel_sort(starts, workers, low_bits=32)
            # Gather lengths, sweep, scatter the verdicts: the random accesses
            # run in parallel, the sweep over contiguous arrays.
            var flags = List[Bool](length=count, fill=True)
            var gather = _RefOrderTask(
                mems.unsafe_ptr().unsafe_offset(lo), starts.unsafe_ptr().as_imm(), lengths.unsafe_ptr(), count
            )
            parallel_for(gather, _REF_ORDER_TASKS, workers)
            _mark_contained(starts, 32, lengths, flags.unsafe_ptr(), workers)
            var scatter = _ClearTask(
                starts.unsafe_ptr().as_imm(), flags.unsafe_ptr().as_imm(), keep.unsafe_ptr().unsafe_offset(lo), count
            )
            parallel_for(scatter, _REF_ORDER_TASKS, workers)
        lo = hi
    return _collect(mems, keep, workers)


def _collect(mems: List[Mem], keep: List[Bool], workers: Int) -> List[Mem]:
    """The MEMs whose `keep` is set, in order: counted per block, then
    copied into place."""
    var kept = List[Int](length=_REF_ORDER_TASKS + 1, fill=0)
    var counting = _CollectTask(
        mems.unsafe_ptr(),
        keep.unsafe_ptr().as_imm(),
        kept.unsafe_ptr(),
        Optional[MutPointer[Mem, MutUntrackedOrigin]](None),
        len(mems),
    )
    parallel_for(counting, _REF_ORDER_TASKS, workers)
    for b in range(_REF_ORDER_TASKS):
        kept[b + 1] += kept[b]
    var out = List[Mem](unsafe_uninit_length=kept[_REF_ORDER_TASKS])
    var copying = _CollectTask(
        mems.unsafe_ptr(),
        keep.unsafe_ptr().as_imm(),
        kept.unsafe_ptr(),
        Optional(out.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()),
        len(mems),
    )
    parallel_for(copying, _REF_ORDER_TASKS, workers)
    _ = len(out)
    return out^


struct _GroupFillTask[mems_origin: ImmOrigin, starts_origin: MutOrigin, lengths_origin: MutOrigin](ParallelTask):
    """Fills a group's interval arrays: query starts and lengths, or (`by_ref`)
    reference keys `ref_pos << 32 | member`."""

    var mems: ImmPointer[Mem, Self.mems_origin]
    var starts: MutPointer[UInt64, Self.starts_origin]
    var lengths: MutPointer[UInt32, Self.lengths_origin]
    var count: Int
    var by_ref: Bool

    def __init__(
        out self,
        mems: ImmPointer[Mem, Self.mems_origin],
        starts: MutPointer[UInt64, Self.starts_origin],
        lengths: MutPointer[UInt32, Self.lengths_origin],
        count: Int,
        *,
        by_ref: Bool,
    ):
        self.mems = mems
        self.starts = starts
        self.lengths = lengths
        self.count = count
        self.by_ref = by_ref

    def run(self, index: Int):
        var per = ceildiv(self.count, _REF_ORDER_TASKS)
        for i in range(index * per, min(self.count, (index + 1) * per)):
            var m = self.mems[unsafe_offset=i]
            if self.by_ref:
                self.starts[unsafe_offset=i] = (UInt64(m.ref_pos) << 32) | UInt64(i)
            else:
                self.starts[unsafe_offset=i] = UInt64(m.query_pos)
                self.lengths[unsafe_offset=i] = m.length


struct _CollectTask[mems_origin: ImmOrigin, keep_origin: ImmOrigin, kept_origin: MutOrigin](ParallelTask):
    """Counts a block's survivors (no `out`) or copies them to `out` from
    the block's offset in `kept`."""

    var mems: ImmPointer[Mem, Self.mems_origin]
    var keep: ImmPointer[Bool, Self.keep_origin]
    var kept: MutPointer[Int, Self.kept_origin]
    var out: Optional[MutPointer[Mem, MutUntrackedOrigin]]
    var count: Int

    def __init__(
        out self,
        mems: ImmPointer[Mem, Self.mems_origin],
        keep: ImmPointer[Bool, Self.keep_origin],
        kept: MutPointer[Int, Self.kept_origin],
        out_ptr: Optional[MutPointer[Mem, MutUntrackedOrigin]],
        count: Int,
    ):
        self.mems = mems
        self.keep = keep
        self.kept = kept
        self.out = out_ptr
        self.count = count

    def run(self, index: Int):
        var per = ceildiv(self.count, _REF_ORDER_TASKS)
        var lo = index * per
        var hi = min(self.count, lo + per)
        if not self.out:
            var n = 0
            for i in range(lo, hi):
                n += Int(self.keep[unsafe_offset=i])
            self.kept[unsafe_offset=index + 1] = n
            return
        var dst = self.out.value()
        var at = self.kept[unsafe_offset=index]
        for i in range(lo, hi):
            if self.keep[unsafe_offset=i]:
                dst[unsafe_offset=at] = self.mems[unsafe_offset=i]
                at += 1


comptime _REF_ORDER_TASKS = 64
comptime _PARALLEL_MEMS = 1 << 16
"""MEMs below which MUM filtering runs on one thread."""


struct _RefOrderTask[mems_origin: ImmOrigin, order_origin: ImmOrigin, out_origin: MutOrigin](ParallelTask):
    """Copies the lengths of a group's MEMs into reference order. Tracked
    origins keep the buffers alive while the tasks run."""

    var mems: ImmPointer[Mem, Self.mems_origin]
    var order: ImmPointer[UInt64, Self.order_origin]
    """Reference start << 32 | member, in reference order."""
    var lengths: MutPointer[UInt32, Self.out_origin]
    var count: Int

    def __init__(
        out self,
        mems: ImmPointer[Mem, Self.mems_origin],
        order: ImmPointer[UInt64, Self.order_origin],
        lengths: MutPointer[UInt32, Self.out_origin],
        count: Int,
    ):
        self.mems = mems
        self.order = order
        self.lengths = lengths
        self.count = count

    def run(self, index: Int):
        var per = ceildiv(self.count, _REF_ORDER_TASKS)
        for i in range(index * per, min(self.count, (index + 1) * per)):
            var member = Int(self.order[unsafe_offset=i] & 0xFFFFFFFF)
            self.lengths[unsafe_offset=i] = self.mems[unsafe_offset=member].length


struct _ClearTask[order_origin: ImmOrigin, flags_origin: ImmOrigin, keep_origin: MutOrigin](ParallelTask):
    """Clears `keep` for the members whose `flags`, in reference order, are
    cleared."""

    var order: ImmPointer[UInt64, Self.order_origin]
    var flags: ImmPointer[Bool, Self.flags_origin]
    var keep: MutPointer[Bool, Self.keep_origin]
    var count: Int

    def __init__(
        out self,
        order: ImmPointer[UInt64, Self.order_origin],
        flags: ImmPointer[Bool, Self.flags_origin],
        keep: MutPointer[Bool, Self.keep_origin],
        count: Int,
    ):
        self.order = order
        self.flags = flags
        self.keep = keep
        self.count = count

    def run(self, index: Int):
        var per = ceildiv(self.count, _REF_ORDER_TASKS)
        for i in range(index * per, min(self.count, (index + 1) * per)):
            if not self.flags[unsafe_offset=i]:
                self.keep[unsafe_offset=Int(self.order[unsafe_offset=i] & 0xFFFFFFFF)] = False


def _mark_contained(
    starts: List[UInt64], shift: UInt64, lengths: List[UInt32], keep: MutPointer[Bool, _], threads: Int = 1
):
    """Clears `keep[i]` for the intervals `[starts[i] >> shift, ... + lengths[i])`,
    in increasing start, contained in (or equal to) another one.

    Of the intervals sharing a start, all but the longest are contained in
    it, and the longest is when it has a copy or an earlier start reaches as
    far. In parallel: blocks starting at group boundaries find their largest
    end, a prefix maximum gives each block what precedes it, and the blocks
    then sweep independently."""
    var count = len(starts)
    if threads <= 1 or count < 1 << 16:
        _ = _sweep(Span(starts), shift, Span(lengths), keep, 0, count, -1, True)
        return
    var blocks = threads * 4
    var bounds = List[Int](length=blocks + 1, fill=count)
    bounds[0] = 0
    for b in range(1, blocks):
        var at = max(bounds[b - 1], b * count // blocks)
        while at < count and at > 0 and starts[at] >> shift == starts[at - 1] >> shift:
            at += 1
        bounds[b] = at
    var ends = List[Int](length=blocks, fill=-1)
    var scan = _SweepTask(
        starts.unsafe_ptr().as_imm(),
        shift,
        lengths.unsafe_ptr().as_imm(),
        keep.unsafe_origin_cast[MutUntrackedOrigin](),
        bounds.unsafe_ptr().as_imm(),
        ends.unsafe_ptr(),
        count,
        mark=False,
    )
    parallel_for(scan, blocks, threads)
    # Each block's incoming max end: the largest end of all earlier blocks.
    var before = -1
    for b in range(blocks):
        var mine = ends[b]
        ends[b] = before
        before = max(before, mine)
    var mark = _SweepTask(
        starts.unsafe_ptr().as_imm(),
        shift,
        lengths.unsafe_ptr().as_imm(),
        keep.unsafe_origin_cast[MutUntrackedOrigin](),
        bounds.unsafe_ptr().as_imm(),
        ends.unsafe_ptr(),
        count,
        mark=True,
    )
    parallel_for(mark, blocks, threads)
    _ = len(ends)


def _sweep(
    starts: Span[UInt64, _],
    shift: UInt64,
    lengths: Span[UInt32, _],
    keep: MutPointer[Bool, _],
    lo: Int,
    hi: Int,
    incoming: Int,
    mark: Bool,
) -> Int:
    """Sweeps groups `[lo, hi)` (whole groups) after intervals ending at
    `incoming`; with `mark`, clears `keep` for the contained ones. Returns
    the largest end seen."""
    var max_end = incoming
    var i = lo
    while i < hi:
        var start = Int(starts[i] >> shift)
        var j = i
        var longest = 0
        var copies = 0
        while j < hi and Int(starts[j] >> shift) == start:
            var length = Int(lengths[j])
            if length > longest:
                longest = length
                copies = 0
            copies += Int(length == longest)
            j += 1
        var end = start + longest
        if mark:
            var dropped = copies > 1 or end <= max_end
            for t in range(i, j):
                if dropped or Int(lengths[t]) < longest:
                    keep[unsafe_offset=t] = False
        max_end = max(max_end, end)
        i = j
    return max_end


struct _SweepTask[
    starts_origin: ImmOrigin, lengths_origin: ImmOrigin, bounds_origin: ImmOrigin, ends_origin: MutOrigin
](ParallelTask):
    """One block of a parallel `_mark_contained`: its largest end, or (`mark`)
    its sweep after the block's incoming end in `ends`."""

    var starts: ImmPointer[UInt64, Self.starts_origin]
    var shift: UInt64
    var lengths: ImmPointer[UInt32, Self.lengths_origin]
    var keep: MutPointer[Bool, MutUntrackedOrigin]
    var bounds: ImmPointer[Int, Self.bounds_origin]
    var ends: MutPointer[Int, Self.ends_origin]
    var count: Int
    var mark: Bool

    def __init__(
        out self,
        starts: ImmPointer[UInt64, Self.starts_origin],
        shift: UInt64,
        lengths: ImmPointer[UInt32, Self.lengths_origin],
        keep: MutPointer[Bool, MutUntrackedOrigin],
        bounds: ImmPointer[Int, Self.bounds_origin],
        ends: MutPointer[Int, Self.ends_origin],
        count: Int,
        *,
        mark: Bool,
    ):
        self.starts = starts
        self.shift = shift
        self.lengths = lengths
        self.keep = keep
        self.bounds = bounds
        self.ends = ends
        self.count = count
        self.mark = mark

    def run(self, index: Int):
        var lo = self.bounds[unsafe_offset=index]
        var hi = self.bounds[unsafe_offset=index + 1]
        var starts = Span(unsafe_ptr=self.starts, length=self.count)
        var lengths = Span(unsafe_ptr=self.lengths, length=self.count)
        var end = _sweep(starts, self.shift, lengths, self.keep, lo, hi, self.ends[unsafe_offset=index], self.mark)
        if not self.mark:
            self.ends[unsafe_offset=index] = end
