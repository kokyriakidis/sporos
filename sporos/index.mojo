"""A k-mer seed index over reference sequences, for seed-and-extend aligners.

The index stores selected reference k-mers (minimizers, closed syncmers, or
every k-mer) keyed by canonical k-mer, so a query k-mer finds reference
occurrences on both strands. Lookups take the query's minimizers (or every
query k-mer), probe the table through a software prefetch pipeline, and
return hits with occurrence counts for abundance filtering.

Entries are `fingerprint << 33 | reverse << 32 | position`, where the
fingerprint is the low 31 bits of the canonical k-mer's `mix64` hash and
`reverse` says the reference k-mer is the reverse complement of its
canonical form. Buckets are sorted, so the occurrences of one k-mer are
contiguous. When bucket and fingerprint bits cannot tell all k-mers apart
(`2k > bucket bits + 31`), candidates are checked against the text.
"""

from std.bit import bit_reverse
from std.math import ceildiv

from sporos.alphabet import CodeBuffer, N_CODE, encode_into
from sporos.mems import mix64
from sporos.minimizers import minimizer_positions_from_codes
from sporos.syncmers import syncmer_positions_from_codes
from sporos.table import AnchorSource, HashTable
from sporos.text import SeqText
from sporos.window import PADDING
from std.sys.intrinsics import PrefetchOptions, prefetch

comptime _AHEAD = PrefetchOptions().for_read().high_locality()
comptime _CHUNK = 1 << 16
"""Selected positions per index-build chunk."""
comptime _LAG = 16
comptime _RING = 64
comptime _PAIRS: UInt64 = 0x5555555555555555


@fieldwise_init
struct SeedHit(Comparable, ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """A query k-mer that occurs in the reference.

    `query_pos` and `ref_pos` are k-mer start positions (`ref_pos` in the
    index text; see `SeedIndex.text`); `occurrences` counts the k-mer's
    reference occurrences on both strands; `reverse` means the query k-mer
    matches the reverse complement of the reference at `ref_pos`.
    """

    var query_pos: UInt32
    var ref_pos: UInt32
    var occurrences: UInt32
    var reverse: Bool

    def __lt__(self, other: Self) -> Bool:
        if self.query_pos != other.query_pos:
            return self.query_pos < other.query_pos
        return self.ref_pos < other.ref_pos

    def __eq__(self, other: Self) -> Bool:
        return (
            self.query_pos == other.query_pos
            and self.ref_pos == other.ref_pos
            and self.occurrences == other.occurrences
            and self.reverse == other.reverse
        )


comptime MINIMIZER_SEEDS = 0
comptime SYNCMER_SEEDS = 1
comptime ALL_KMERS = 2


@inline(.always)
def kmer_mask(k: Int) -> UInt64:
    return UInt64.MAX if k == 32 else (UInt64(1) << UInt64(2 * k)) - 1


@inline(.always)
def reverse_complement(kmer: UInt64, k: Int) -> UInt64:
    """Reverse complement of a 2-bit k-mer (base i in bits 2i, A=0..T=3)."""
    var r = bit_reverse(~kmer)
    # bit_reverse also swapped the two bits of every base; swap them back.
    r = ((r >> 1) & _PAIRS) | ((r & _PAIRS) << 1)
    return r >> UInt64(64 - 2 * k)


@inline(.always)
def canonical(kmer: UInt64, k: Int) -> Tuple[UInt64, Bool]:
    """The smaller of a k-mer and its reverse complement, and whether it was
    the reverse complement."""
    var rc = reverse_complement(kmer, k)
    if rc < kmer:
        return (rc, True)
    return (kmer, False)


def select_positions[
    mut: Bool, //, origin: Origin[mut=mut]
](
    seqs: List[Span[UInt8, origin]],
    k: Int,
    w: Int,
    scheme: Int,
    threads: Int,
    open_offset: Optional[Int] = None,
    *,
    open_offsets: List[Int] = List[Int](),
) raises -> List[UInt32]:
    """Selected k-mer start positions in the concatenated text (sequences
    separated by one base, as in `SeqText`)."""
    var total = 0
    for s in seqs:
        total += len(s) + 1
    var codes = CodeBuffer(total, PADDING)
    var p = codes.unsafe_ptr()
    var off = 0
    for s in seqs:
        encode_into(s, p.unsafe_offset(off))
        p[unsafe_offset=off + len(s)] = N_CODE
        off += len(s) + 1
    if scheme == MINIMIZER_SEEDS:
        return minimizer_positions_from_codes(p, total, k, w, True, threads)
    if scheme == SYNCMER_SEEDS:
        return syncmer_positions_from_codes(p, total, k, w, True, open_offset, threads, open_offsets=open_offsets)
    var out = List[UInt32]()
    var run = 0
    for i in range(total):
        run = run + 1 if p[unsafe_offset=i] < N_CODE else 0
        if run >= k:
            out.append(UInt32(i - k + 1))
    return out^


struct _SelectedKmers[origin: ImmOrigin](AnchorSource):
    """The selected reference positions as table anchors (tracked, so the
    positions stay alive while the table is built)."""

    var text: ImmPointer[SeqText, ImmUntrackedOrigin]
    var positions: ImmPointer[UInt32, Self.origin]
    var count_all: Int
    var k: Int

    def __init__(
        out self,
        text: ImmPointer[SeqText, ImmUntrackedOrigin],
        positions: ImmPointer[UInt32, Self.origin],
        count_all: Int,
        k: Int,
    ):
        self.text = text
        self.positions = positions
        self.count_all = count_all
        self.k = k

    def chunks(self) -> Int:
        return max(1, ceildiv(self.count_all, _CHUNK))

    def count(self, chunk: Int) -> Int:
        return max(0, min(_CHUNK, self.count_all - chunk * _CHUNK))

    def fill(
        self,
        chunk: Int,
        hashes: MutPointer[UInt64, MutUntrackedOrigin],
        entries: MutPointer[UInt64, MutUntrackedOrigin],
    ):
        var mask = kmer_mask(self.k)
        var lo = chunk * _CHUNK
        for i in range(lo, lo + self.count(chunk)):
            var pos = Int(self.positions[unsafe_offset=i])
            var form = canonical(self.text[].window(pos) & mask, self.k)
            var h = mix64(form[0])
            hashes[unsafe_offset=i - lo] = h
            entries[unsafe_offset=i - lo] = ((h & 0x7FFFFFFF) << 33) | (UInt64(Int(form[1])) << 32) | UInt64(pos)


struct SeedIndex(Movable):
    """Selected reference k-mers (k <= 32), keyed by canonical k-mer."""

    var text: SeqText
    """The packed reference; `ref_pos` values are offsets into it, and
    `text.sequence_of` / `text.starts` map them to sequences."""
    var table: HashTable
    var k: Int
    var w: Int
    var scheme: Int
    var open_offsets: List[Int]
    """For syncmer seeds: open syncmers with these s-mer offsets, else closed."""
    var verify: Bool

    def __init__[
        mut: Bool, //, origin: Origin[mut=mut]
    ](
        out self,
        seqs: List[Span[UInt8, origin]],
        k: Int,
        w: Int,
        scheme: Int = MINIMIZER_SEEDS,
        threads: Int = 1,
        open_offset: Optional[Int] = None,
        *,
        open_offsets: List[Int] = List[Int](),
    ) raises:
        """Indexes the (w, k)-minimizers, the syncmers with s-mer length `w`
        (closed, or open at `open_offset` or any of `open_offsets`), or (`ALL_KMERS`) every k-mer of
        `seqs`."""
        if k < 1 or k > 32:
            raise Error(t"seed index k must be in [1, 32], got {k}")
        self.k = k
        self.w = w
        self.scheme = scheme
        if open_offset and len(open_offsets) > 0:
            raise Error("give open_offset or open_offsets, not both")
        self.open_offsets = open_offsets.copy()
        if case .Some(t) = open_offset:
            self.open_offsets.append(t)
        self.text = SeqText(seqs)
        if len(self.text) >= Int(UInt32.MAX):
            raise Error("reference text must be shorter than 2^32")
        var positions = select_positions(seqs, k, w, scheme, threads, open_offsets=self.open_offsets)
        var source = _SelectedKmers(
            Pointer(to=self.text).unsafe_origin_cast[MutUntrackedOrigin](),
            positions.unsafe_ptr().as_imm(),
            len(positions),
            k,
        )
        self.table = HashTable(source, threads, sort_buckets=True)
        var bucket_bits = 64 - Int(self.table.shift)
        self.verify = 2 * k > bucket_bits + 31

    def occurrence_cutoff(self, keep_fraction: Float64) -> Int:
        """The smallest occurrence count such that `keep_fraction` of the
        distinct indexed k-mers occur at most that often (for `lookup`'s
        `max_occurrences`, like GraphAligner's quantile cutoff). Counts are
        exact when `verify` is false; otherwise k-mers sharing a fingerprint
        are counted together."""
        var counts = List[Int]()
        var ent = self.table.entries.unsafe_ptr()
        var off = self.table.offsets.unsafe_ptr()
        for b in range(len(self.table.offsets) - 1):
            var e = Int(off[unsafe_offset=b])
            var end = Int(off[unsafe_offset=b + 1])
            while e < end:
                var fp = ent[unsafe_offset=e] >> 33
                var run = 0
                while e < end and ent[unsafe_offset=e] >> 33 == fp:
                    run += 1
                    e += 1
                if len(counts) <= run:
                    counts.resize(run + 1, 0)
                counts[run] += 1
        var distinct = 0
        for c in counts:
            distinct += c
        var keep = Int(Float64(distinct) * keep_fraction)
        var seen = 0
        for occ in range(len(counts)):
            seen += counts[occ]
            if seen >= keep:
                return max(1, occ)
        return max(1, len(counts) - 1)

    def query_positions(self, query: Span[UInt8, _], all_kmers: Bool, threads: Int = 1) raises -> List[UInt32]:
        """The query k-mers to look up: its seeds under the index's scheme,
        or every k-mer."""
        var seqs = List[Span[UInt8, query.origin]]()
        seqs.append(query)
        return select_positions(
            seqs, self.k, self.w, ALL_KMERS if all_kmers else self.scheme, threads, open_offsets=self.open_offsets
        )

    def lookup(
        self,
        query: Span[UInt8, _],
        mut hits: List[SeedHit],
        all_kmers: Bool = False,
        max_occurrences: Int = 0,
    ) raises:
        """Appends the hits of `query`'s seeds (or of all its k-mers), skipping
        k-mers with more than `max_occurrences` reference occurrences
        (0 = no limit). Hits come in query order."""
        var seqs = List[Span[UInt8, query.origin]]()
        seqs.append(query)
        var qtext = SeqText(seqs)
        var positions = self.query_positions(query, all_kmers)
        var mask = kmer_mask(self.k)
        var off = self.table.offsets.unsafe_ptr()
        var qpos = Array[UInt32, _RING](fill=0)
        var hashes = Array[UInt64, _RING](fill=0)
        var strands = Array[Bool, _RING](fill=False)
        var lows = Array[UInt32, _RING](fill=0)
        var highs = Array[UInt32, _RING](fill=0)
        var n = 0
        for p in positions:
            var form = canonical(qtext.window(Int(p)) & mask, self.k)
            var h = mix64(form[0])
            var slot = n & (_RING - 1)
            qpos[slot] = p
            hashes[slot] = h
            strands[slot] = form[1]
            prefetch[_AHEAD](off.unsafe_offset(self.table.bucket(h)))
            if n >= _LAG:
                self._bounds(n - _LAG, hashes, lows, highs)
            if n >= 2 * _LAG:
                self._report(
                    (n - 2 * _LAG) & (_RING - 1), qpos, hashes, strands, lows, highs, qtext, max_occurrences, hits
                )
            n += 1
        for m in range(max(0, n - _LAG), n):
            self._bounds(m, hashes, lows, highs)
        for m in range(max(0, n - 2 * _LAG), n):
            self._report(m & (_RING - 1), qpos, hashes, strands, lows, highs, qtext, max_occurrences, hits)

    @inline(.always)
    def _bounds(
        self,
        n: Int,
        hashes: Array[UInt64, _RING],
        mut lows: Array[UInt32, _RING],
        mut highs: Array[UInt32, _RING],
    ):
        var slot = n & (_RING - 1)
        var b = self.table.bucket(hashes[slot])
        var off = self.table.offsets.unsafe_ptr()
        lows[slot] = off[unsafe_offset=b]
        highs[slot] = off[unsafe_offset=b + 1]
        prefetch[_AHEAD](self.table.entries.unsafe_ptr().unsafe_offset(Int(lows[slot])))

    @inline(.always)
    def _report(
        self,
        slot: Int,
        qpos: Array[UInt32, _RING],
        hashes: Array[UInt64, _RING],
        strands: Array[Bool, _RING],
        lows: Array[UInt32, _RING],
        highs: Array[UInt32, _RING],
        qtext: SeqText,
        max_occurrences: Int,
        mut hits: List[SeedHit],
    ):
        var ent = self.table.entries.unsafe_ptr()
        var fp = hashes[slot] & 0x7FFFFFFF
        # Entries are sorted, so this k-mer's occurrences are one run.
        var first = Int(lows[slot])
        var end = Int(highs[slot])
        while first < end and (ent[unsafe_offset=first] >> 33) < fp:
            first += 1
        var last = first
        while last < end and (ent[unsafe_offset=last] >> 33) == fp:
            last += 1
        if first == last:
            return
        var query_kmer = UInt64(0)
        var mask = kmer_mask(self.k)
        if self.verify:
            query_kmer = canonical(qtext.window(Int(qpos[slot])) & mask, self.k)[0]
            var kept = first
            for e in range(first, last):
                var pos = Int(ent[unsafe_offset=e] & 0xFFFFFFFF)
                if canonical(self.text.window(pos) & mask, self.k)[0] == query_kmer:
                    kept += 1
            if kept == first:
                return
        var occurrences = last - first
        if self.verify:
            occurrences = 0
            for e in range(first, last):
                var pos = Int(ent[unsafe_offset=e] & 0xFFFFFFFF)
                occurrences += Int(canonical(self.text.window(pos) & mask, self.k)[0] == query_kmer)
        if max_occurrences > 0 and occurrences > max_occurrences:
            return
        for e in range(first, last):
            var entry = ent[unsafe_offset=e]
            var pos = Int(entry & 0xFFFFFFFF)
            if self.verify and canonical(self.text.window(pos) & mask, self.k)[0] != query_kmer:
                continue
            var ref_reverse = Bool((entry >> 32) & 1)
            hits.append(SeedHit(qpos[slot], UInt32(pos), UInt32(occurrences), ref_reverse != strands[slot]))
