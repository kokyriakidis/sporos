"""A k-mer seed index over reference sequences, for seed-and-extend aligners.

The index stores selected reference k-mers (minimizers, closed syncmers, or
every k-mer) keyed by canonical k-mer, so a query k-mer finds reference
occurrences on both strands. Lookups take the query's minimizers (or every
query k-mer), probe the table through a software prefetch pipeline, and
return hits with occurrence counts for abundance filtering.

K-mers of up to 32 bases are held in a `UInt64` and hashed by `kmer_hash`,
a bijection on 2k bits: its top bits pick the bucket and its low 31 bits are
the fingerprint. Longer k-mers, up to 64 bases, are held in a `UInt128` and
hashed to 64 bits by `wide_kmer_hash`, which no bucket and fingerprint can
invert. Entries are `fingerprint << 33 | reverse << 32 | position`, where
`reverse` says the reference k-mer is the reverse complement of its canonical
form. Buckets are sorted, so the occurrences of one k-mer are contiguous.
Bucket and fingerprint identify a k-mer exactly when `k <= 32` and
`bucket bits + 31 >= 2k`; otherwise (`verify`) candidates are checked against
the text.
"""

from std.bit import bit_reverse, count_trailing_zeros, pop_count
from std.math import ceildiv, iota
from std.memory import pack_bits

from std.os import abort
from std.sys import size_of

from sporos.alphabet import CodeBuffer, N_CODE, encode_into
from sporos.batch import SeedBatch, minimizers_batch, syncmers_batch
from sporos.parallel import ParallelTask, parallel_for
from sporos.minimizers import minimizer_positions, minimizer_positions_from_codes
from sporos.syncmers import syncmer_positions, syncmer_positions_from_codes
from sporos.table import AnchorSource, HashTable
from sporos.text import SeqText
from sporos.window import PADDING
from std.sys.intrinsics import PrefetchOptions, prefetch

comptime _AHEAD = PrefetchOptions().for_read().high_locality()
comptime _CHUNK = 1 << 16
"""Selected positions per index-build chunk."""
comptime _LAG = 16
comptime _SCAN = 8
"""Entries one vector compare covers (`table.ENTRY_PADDING`)."""
comptime _BUCKET_LOAD = 4
"""Entries per bucket: a quarter of the offsets keeps more of them cached,
and a vector compare scans a bucket (measured best on hg38 chr21)."""
comptime _RING = 64
comptime _TASK_BASES = 1 << 20
comptime _TASK_READS = 1 << 12
"""A lookup task's reads: up to this many bases or reads."""


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


comptime MAX_K = 64
"""The longest seed-index k-mer: 32 bases fit a `UInt64`, 64 a `UInt128`."""


@inline(.always)
def kmer_mask[D: DType = DType.uint64](k: Int) -> Scalar[D]:
    """The low 2k bits of a `Scalar[D]` (all of them when they hold k bases)."""
    comptime bits = size_of[Scalar[D]]() * 8
    return ~Scalar[D](0) if 2 * k >= bits else (Scalar[D](1) << Scalar[D](2 * k)) - 1


@inline(.always)
def reverse_complement[D: DType, //](kmer: Scalar[D], k: Int) -> Scalar[D]:
    """Reverse complement of a 2-bit k-mer (base i in bits 2i, A=0..T=3)."""
    comptime bits = size_of[Scalar[D]]() * 8
    comptime pairs = ~Scalar[D](0) // 3  # 0x55...5: the low bit of every base
    var r = bit_reverse(~kmer)
    # bit_reverse also swapped the two bits of every base; swap them back.
    r = ((r >> 1) & pairs) | ((r & pairs) << 1)
    return r >> Scalar[D](bits - 2 * k)


@inline(.always)
def kmer_hash(kmer: UInt64, k: Int) -> UInt64:
    """An invertible mix of a 2k-bit k-mer, in the top 2k bits of the result
    (xorshifts and odd multiplies modulo 2^2k), so any of its top bits
    together with its low 31 (`fingerprint`) that cover all 2k identify the
    k-mer."""
    var mask = kmer_mask(k)
    var half = UInt64(k)
    var h = kmer
    h ^= h >> half
    h = (h * 0xFF51AFD7ED558CCD) & mask
    h ^= h >> half
    h = (h * 0xC4CEB9FE1A85EC53) & mask
    h ^= h >> half
    return h << UInt64(64 - 2 * k)


@inline(.always)
def wide_kmer_hash(kmer: UInt128) -> UInt64:
    """A 64-bit hash of a k-mer of more than 32 bases (MurmurHash3's finalizer
    over both halves). It is not invertible, so lookups check candidates
    against the text."""
    var h = UInt64(kmer >> 64) ^ 0x9E3779B97F4A7C15
    h = _murmur_mix(h) ^ UInt64(kmer & 0xFFFFFFFFFFFFFFFF)
    return _murmur_mix(h)


@inline(.always)
def _murmur_mix(x: UInt64) -> UInt64:
    var h = x
    h ^= h >> 33
    h *= 0xFF51AFD7ED558CCD
    h ^= h >> 33
    h *= 0xC4CEB9FE1A85EC53
    h ^= h >> 33
    return h


@inline(.always)
def fingerprint(hash: UInt64, k: Int) -> UInt64:
    """The low 31 bits of a `kmer_hash` value's 2k-bit mix, or of a
    `wide_kmer_hash` value."""
    return (hash >> UInt64(max(0, 64 - 2 * k))) & 0x7FFFFFFF


@inline(.always)
def canonical[D: DType, //](kmer: Scalar[D], k: Int) -> Tuple[Scalar[D], Bool]:
    """The smaller of a k-mer and its reverse complement, and whether it was
    the reverse complement."""
    var rc = reverse_complement(kmer, k)
    if rc < kmer:
        return (rc, True)
    return (kmer, False)


@inline(.always)
def _hash[D: DType, //](kmer: Scalar[D], k: Int) -> UInt64:
    """The table hash of a canonical k-mer: `kmer_hash` up to 32 bases, else
    `wide_kmer_hash`."""
    comptime if D == DType.uint64:
        return kmer_hash(rebind[UInt64](kmer), k)
    else:
        return wide_kmer_hash(rebind[UInt128](kmer))


@inline(.always)
def _kmer[D: DType](text: SeqText, pos: Int, k: Int) -> Scalar[D]:
    """The k-mer at `pos` of `text`: one 32-base window, or two for `UInt128`."""
    comptime if D == DType.uint64:
        return rebind[Scalar[D]](text.window(pos) & kmer_mask(k))
    else:
        var both = UInt128(text.window(pos)) | (UInt128(text.window(pos + 32)) << 64)
        return rebind[Scalar[D]](both & kmer_mask[DType.uint128](k))


@inline(.always)
def _covers_masked[D: DType](text: SeqText, pos: Int, k: Int) -> Bool:
    """Whether the k-mer at `pos` covers a non-ACGT base."""
    comptime if D == DType.uint64:
        return text.mask(pos) & kmer_mask(k) != 0
    else:
        var both = UInt128(text.mask(pos)) | (UInt128(text.mask(pos + 32)) << 64)
        return both & kmer_mask[DType.uint128](k) != 0


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


struct _SelectedKmers[D: DType, origin: ImmOrigin](AnchorSource):
    """The selected reference positions as table anchors (tracked, so the
    positions stay alive while the table is built), their k-mers held in
    `Scalar[D]`."""

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

    def capacity(self, chunk: Int) -> Int:
        return max(0, min(_CHUNK, self.count_all - chunk * _CHUNK))

    def fill(
        self,
        chunk: Int,
        hashes: MutPointer[UInt64, MutUntrackedOrigin],
        entries: MutPointer[UInt64, MutUntrackedOrigin],
    ) -> Int:
        var lo = chunk * _CHUNK
        for i in range(lo, lo + self.capacity(chunk)):
            var pos = Int(self.positions[unsafe_offset=i])
            var form = canonical(_kmer[Self.D](self.text[], pos, self.k), self.k)
            var h = _hash(form[0], self.k)
            hashes[unsafe_offset=i - lo] = h
            entries[unsafe_offset=i - lo] = (fingerprint(h, self.k) << 33) | (UInt64(Int(form[1])) << 32) | UInt64(pos)
        return self.capacity(chunk)


struct SeedIndex(Movable):
    """Selected reference k-mers (k <= 64), keyed by canonical k-mer."""

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
        if k < 1 or k > MAX_K:
            raise Error(t"seed index k must be in [1, {MAX_K}], got {k}")
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
        var text = Pointer(to=self.text).unsafe_origin_cast[MutUntrackedOrigin]()
        if k <= 32:
            self.table = _table[DType.uint64](text, positions, k, threads)
        else:
            self.table = _table[DType.uint128](text, positions, k, threads)
        var bucket_bits = 64 - Int(self.table.shift)
        self.verify = k > 32 or 2 * k > bucket_bits + 31

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
        if self.scheme == MINIMIZER_SEEDS and not all_kmers:
            return minimizer_positions(query, self.k, self.w, True, threads)
        if self.scheme == SYNCMER_SEEDS and not all_kmers:
            return syncmer_positions(query, self.k, self.w, True, None, threads, open_offsets=self.open_offsets)
        var seqs = List[Span[UInt8, query.origin]]()
        seqs.append(query)
        return select_positions(seqs, self.k, self.w, ALL_KMERS, threads)

    def lookup(
        self,
        query: Span[UInt8, _],
        mut hits: List[SeedHit],
        all_kmers: Bool = False,
        max_occurrences: Int = 0,
    ) raises:
        """Appends the hits of `query`'s seeds (or of all its k-mers), skipping
        k-mers with more than `max_occurrences` reference occurrences
        (0 = no limit). Hits come in query order. For many reads,
        `lookup_batch` is much faster."""
        var reads = List[Span[UInt8, query.origin]]()
        reads.append(query)
        var count = 0
        self._lookup_reads(reads.unsafe_ptr(), 1, all_kmers, max_occurrences, hits, Pointer(to=count))

    def lookup_batch[
        mut: Bool, //, origin: Origin[mut=mut]
    ](
        self, reads: List[Span[UInt8, origin]], all_kmers: Bool = False, max_occurrences: Int = 0, threads: Int = 1
    ) -> SeedHits:
        """The hits of every read, as `lookup` gives them. Reads stream through
        one probe pipeline, so its prefetches cover read boundaries."""
        var bounds: List[Int] = [0]
        var bases = 0
        for i in range(len(reads)):
            bases += len(reads[i])
            if bases >= _TASK_BASES or i + 1 - bounds[len(bounds) - 1] >= _TASK_READS:
                bounds.append(i + 1)
                bases = 0
        if bounds[len(bounds) - 1] != len(reads):
            bounds.append(len(reads))
        var tasks = len(bounds) - 1
        var counts = List[Int](length=len(reads), fill=0)
        var results = List[List[SeedHit]](length=tasks, fill=List[SeedHit]())
        if tasks > 0:
            var task = _LookupTask[origin](
                Pointer(to=self).unsafe_origin_cast[ImmUntrackedOrigin](),
                reads.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
                bounds^,
                all_kmers,
                max_occurrences,
                results.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
                counts.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
            )
            parallel_for(task, tasks, threads)
        var offsets = _prefix_sums(counts)
        var hits = List[SeedHit](capacity=offsets[len(offsets) - 1])
        for r in results:
            hits.extend(Span(r))
        return SeedHits(hits^, offsets^)

    def _lookup_reads[
        mut: Bool, //, origin: Origin[mut=mut]
    ](
        self,
        reads: ImmPointer[Span[UInt8, origin], _],
        count: Int,
        all_kmers: Bool,
        max_occurrences: Int,
        mut hits: List[SeedHit],
        counts: MutPointer[Int, _],
    ) raises:
        """Appends the hits of `count` reads, counting each read's hits."""
        if self.k <= 32:
            self._lookup_reads_of[DType.uint64](reads, count, all_kmers, max_occurrences, hits, counts)
        else:
            self._lookup_reads_of[DType.uint128](reads, count, all_kmers, max_occurrences, hits, counts)

    def _lookup_reads_of[
        mut: Bool, //, D: DType, origin: Origin[mut=mut]
    ](
        self,
        reads: ImmPointer[Span[UInt8, origin], _],
        count: Int,
        all_kmers: Bool,
        max_occurrences: Int,
        mut hits: List[SeedHit],
        counts: MutPointer[Int, _],
    ) raises:
        """`_lookup_reads` with k-mers held in `Scalar[D]`."""
        var spans = List[Span[UInt8, origin]](capacity=count)
        for r in range(count):
            spans.append(reads[unsafe_offset=r])
        var every = all_kmers or self.scheme == ALL_KMERS
        var seeds = SeedBatch(List[UInt32](), [0])
        if not every:
            if self.scheme == MINIMIZER_SEEDS:
                seeds = minimizers_batch(spans, self.k, self.w, True)
            else:
                seeds = syncmers_batch(spans, self.k, self.w, True, None, open_offsets=self.open_offsets)
        var selected = seeds.positions.unsafe_ptr()
        var text = SeqText(spans)
        var ring = Array[_Probe[D], _RING](fill=_Probe[D](0, 0, 0, 0, False, 0, 0))
        var off = self.table.offsets.unsafe_ptr()
        var n = 0
        for r in range(count):
            counts[unsafe_offset=r] = 0
            var base = text.starts[r]
            var probes = len(spans[r]) - self.k + 1 if every else seeds.offsets[r + 1] - seeds.offsets[r]
            for i in range(probes):
                var start: Int
                if every:
                    start = i
                    if _covers_masked[D](text, base + start, self.k):
                        continue
                else:
                    start = Int(selected[unsafe_offset=seeds.offsets[r] + i])
                var form = canonical(_kmer[D](text, base + start, self.k), self.k)
                var h = _hash(form[0], self.k)
                ring[n & (_RING - 1)] = _Probe[D](r, UInt32(start), form[0], h, form[1], 0, 0)
                prefetch[_AHEAD](off.unsafe_offset(self.table.bucket(h)))
                if n >= _LAG:
                    self._bounds(ring[(n - _LAG) & (_RING - 1)])
                if n >= 2 * _LAG:
                    self._report(ring[(n - 2 * _LAG) & (_RING - 1)], max_occurrences, hits, counts)
                n += 1
        for m in range(max(0, n - _LAG), n):
            self._bounds(ring[m & (_RING - 1)])
        for m in range(max(0, n - 2 * _LAG), n):
            self._report(ring[m & (_RING - 1)], max_occurrences, hits, counts)

    @inline(.always)
    def _bounds[D: DType](self, mut probe: _Probe[D]):
        var b = self.table.bucket(probe.hash)
        var off = self.table.offsets.unsafe_ptr()
        probe.low = off[unsafe_offset=b]
        probe.high = off[unsafe_offset=b + 1]
        var first = self.table.entries.unsafe_ptr().unsafe_offset(Int(probe.low))
        prefetch[_AHEAD](first)
        # The bucket scan's vector load usually spans two lines.
        prefetch[_AHEAD](first.unsafe_offset(_SCAN - 1))

    @inline(.always)
    def _report[
        D: DType
    ](self, probe: _Probe[D], max_occurrences: Int, mut hits: List[SeedHit], counts: MutPointer[Int, _]):
        var ent = self.table.entries.unsafe_ptr()
        var fp = fingerprint(probe.hash, self.k)
        # Entries are sorted, so this k-mer's occurrences are one run.
        var first = Int(probe.low)
        var end = Int(probe.high)
        var last: Int
        if end - first <= _SCAN:
            # Small bucket: one vector compare finds the run.
            var lanes = ent.unsafe_offset(first).unsafe_load[width=_SCAN]() >> 33
            var hit = lanes.eq(fp) & iota[DType.uint64, _SCAN]().lt(UInt64(end - first))
            var bits = pack_bits[DType.uint8](hit)
            if bits == 0:
                return
            first += Int(count_trailing_zeros(bits))
            last = first + Int(pop_count(bits))
        else:
            while first < end and (ent[unsafe_offset=first] >> 33) < fp:
                first += 1
            if (
                not self.verify
                and max_occurrences > 0
                and first + max_occurrences < end
                and (ent[unsafe_offset=first + max_occurrences] >> 33) == fp
            ):
                return  # more than `max_occurrences`, without counting them all
            last = first
            while last < end and (ent[unsafe_offset=last] >> 33) == fp:
                last += 1
            if first == last:
                return
        var occurrences = last - first
        if self.verify:
            occurrences = 0
            for e in range(first, last):
                var pos = Int(ent[unsafe_offset=e] & 0xFFFFFFFF)
                occurrences += Int(canonical(_kmer[D](self.text, pos, self.k), self.k)[0] == probe.kmer)
            if occurrences == 0:
                return
        if max_occurrences > 0 and occurrences > max_occurrences:
            return
        var before = len(hits)
        for e in range(first, last):
            var entry = ent[unsafe_offset=e]
            var pos = Int(entry & 0xFFFFFFFF)
            if self.verify and canonical(_kmer[D](self.text, pos, self.k), self.k)[0] != probe.kmer:
                continue
            var ref_reverse = Bool((entry >> 32) & 1)
            hits.append(SeedHit(probe.qpos, UInt32(pos), UInt32(occurrences), ref_reverse != probe.reverse))
        counts[unsafe_offset=probe.read] += len(hits) - before


def _table[
    D: DType
](text: ImmPointer[SeqText, ImmUntrackedOrigin], positions: List[UInt32], k: Int, threads: Int) raises -> HashTable:
    """The hash table of the k-mers at `positions`, held in `Scalar[D]`."""
    var source = _SelectedKmers[D](text, positions.unsafe_ptr().as_imm(), len(positions), k)
    return HashTable(source, threads, sort_buckets=True, bucket_load=_BUCKET_LOAD)


def _prefix_sums(counts: List[Int]) -> List[Int]:
    """`[0, c0, c0 + c1, ...]`: each read's first hit."""
    var offsets = List[Int](capacity=len(counts) + 1)
    offsets.append(0)
    for c in counts:
        offsets.append(offsets[len(offsets) - 1] + c)
    return offsets^


@fieldwise_init
struct _Probe[D: DType](ImplicitlyCopyable, TrivialRegisterPassable):
    """A query k-mer in the probe pipeline."""

    var read: Int
    var qpos: UInt32
    var kmer: Scalar[Self.D]
    """The canonical query k-mer."""
    var hash: UInt64
    var reverse: Bool
    var low: UInt32
    var high: UInt32


struct SeedHits(Movable, Sized):
    """The hits of a batch of reads: read `i`'s are
    `hits[offsets[i] : offsets[i + 1]]`, in query order."""

    var hits: List[SeedHit]
    var offsets: List[Int]

    def __init__(out self, var hits: List[SeedHit], var offsets: List[Int]):
        self.hits = hits^
        self.offsets = offsets^

    def __len__(self) -> Int:
        return len(self.offsets) - 1

    def of(ref self, i: Int) -> Span[SeedHit, origin_of(self.hits)]:
        return Span(self.hits)[self.offsets[i] : self.offsets[i + 1]]


struct _LookupTask[mut: Bool, //, origin: Origin[mut=mut]](ParallelTask):
    """Looks up one run of reads (`bounds[t] .. bounds[t + 1]`)."""

    var index: ImmPointer[SeedIndex, ImmUntrackedOrigin]
    var reads: ImmPointer[Span[UInt8, Self.origin], ImmUntrackedOrigin]
    var bounds: List[Int]
    var all_kmers: Bool
    var max_occurrences: Int
    var results: MutPointer[List[SeedHit], MutUntrackedOrigin]
    var counts: MutPointer[Int, MutUntrackedOrigin]

    def __init__(
        out self,
        index: ImmPointer[SeedIndex, ImmUntrackedOrigin],
        reads: ImmPointer[Span[UInt8, Self.origin], ImmUntrackedOrigin],
        var bounds: List[Int],
        all_kmers: Bool,
        max_occurrences: Int,
        results: MutPointer[List[SeedHit], MutUntrackedOrigin],
        counts: MutPointer[Int, MutUntrackedOrigin],
    ):
        self.index = index
        self.reads = reads
        self.bounds = bounds^
        self.all_kmers = all_kmers
        self.max_occurrences = max_occurrences
        self.results = results
        self.counts = counts

    def run(self, t: Int):
        var lo = self.bounds[t]
        var hi = self.bounds[t + 1]
        ref out = self.results[unsafe_offset=t]
        try:
            self.index[]._lookup_reads(
                self.reads.unsafe_offset(lo),
                hi - lo,
                self.all_kmers,
                self.max_occurrences,
                out,
                self.counts.unsafe_offset(lo),
            )
        except e:
            abort(String(t"seed lookup failed: {e}"))
