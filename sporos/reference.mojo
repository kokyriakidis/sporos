"""The high-level entry point: a reference of named sequences to match
queries against.

    var genome = Reference.read("genome.fa", threads=8)
    for m in genome.mems(query, min_len=50, mode="mum", reverse=True):
        print(genome.name(m.sequence), m.ref_pos, m.query_pos, m.length, m.reverse)
    var smems = genome.smems(reads, min_len=19)

The indexes behind it are built on first use and kept: the MEM index for each
minimum length in turn, and one FM-index for SMEMs. `sporos.mems` and
`sporos.fmindex` hold the building blocks for finer control.
"""

from sporos.fasta import FastaRecords, read_fasta
from sporos.fmindex import FMIndex
from sporos.mems import Mem, MemIndex, MemParams, Smems, filter_by_occurrences, find_mems
from sporos.text import SeqText


@fieldwise_init
struct Match(Comparable, ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """An exact match: `length` bases at 0-based `ref_pos` in reference
    sequence `sequence` equal the query's at `query_pos`. On the reverse
    strand (`reverse`), `query_pos` is in the query's reverse complement, as
    MUMmer reports it."""

    var sequence: Int
    var ref_pos: Int
    var query_pos: Int
    var length: Int
    var reverse: Bool

    def __lt__(self, other: Self) -> Bool:
        if self.reverse != other.reverse:
            return not self.reverse
        if self.query_pos != other.query_pos:
            return self.query_pos < other.query_pos
        if self.sequence != other.sequence:
            return self.sequence < other.sequence
        return self.ref_pos < other.ref_pos

    def __eq__(self, other: Self) -> Bool:
        return (
            self.sequence == other.sequence
            and self.ref_pos == other.ref_pos
            and self.query_pos == other.query_pos
            and self.length == other.length
            and self.reverse == other.reverse
        )


struct Reference(Movable, Sized):
    """One or more named DNA sequences (ASCII; non-ACGT bases never match)."""

    var records: FastaRecords
    """The sequences: `names`, `seq(i)`, `lengths`."""
    var _text: Optional[SeqText]
    var _index: Optional[MemIndex]
    var _index_min_len: Int
    var _fm: Optional[FMIndex]

    def __init__(out self, var records: FastaRecords):
        """A reference of the records of a parsed FASTA/FASTQ file."""
        self.records = records^
        self._text = None
        self._index = None
        self._index_min_len = 0
        self._fm = None

    def __init__[
        mut: Bool, //, origin: Origin[mut=mut]
    ](out self, seqs: List[Span[UInt8, origin]], names: List[String] = List[String]()):
        """A reference of `seqs` (copied), named `names` or `seq0`, `seq1`, ..."""
        var records = FastaRecords()
        for i in range(len(seqs)):
            records.names.append(names[i] if i < len(names) else String(t"seq{i}"))
            records.starts.append(len(records.data))
            records.lengths.append(len(seqs[i]))
            records.data.extend(seqs[i])
        self = Self(records^)

    @staticmethod
    def read(path: String, threads: Int = 1) raises -> Reference:
        """The sequences of a FASTA or FASTQ file (plain, gzip or BGZF)."""
        return Reference(read_fasta(path, threads))

    def __len__(self) -> Int:
        """The number of sequences."""
        return len(self.records)

    def name(self, i: Int) -> String:
        """The name of sequence `i` (all of them: `records.names`)."""
        return self.records.names[i]

    def seq(ref self, i: Int) -> Span[UInt8, origin_of(self.records.data)]:
        """Sequence `i`."""
        return self.records.seq(i)

    # ----------------------------------------------------------------- MEMs

    def mems[
        mut: Bool, //, origin: Origin[mut=mut]
    ](
        mut self,
        query: Span[UInt8, origin],
        min_len: Int,
        mode: String = "maxmatch",
        forward: Bool = True,
        reverse: Bool = False,
        max_occurrences: Int = 0,
        min_occurrences: Int = 1,
        threads: Int = 1,
    ) raises -> List[Match]:
        """The exact matches of length >= `min_len` between `query` and the
        reference, on the forward strand, the reverse strand or both.

        `mode` (MUMmer's): `"maxmatch"`, every MEM; `"mumreference"`, MEMs
        whose string occurs once in the reference; `"mum"`, once in the
        reference and once in the query. `max_occurrences` N relaxes "once"
        to "at most N times" (for `"maxmatch"` it caps the reference count),
        and `min_occurrences` asks for at least that many reference copies.
        Matches come strand by strand, in query order."""
        if mode != "maxmatch" and mode != "mum" and mode != "mumreference":
            raise Error(t"mode must be maxmatch, mum or mumreference, not {mode}")
        if len(query) >= Int(UInt32.MAX):
            raise Error("queries must be shorter than 2^32 bases")
        var cap = max_occurrences if max_occurrences > 0 else (0 if mode == "maxmatch" else 1)
        if cap > 0 and min_occurrences > cap:
            raise Error(t"min_occurrences {min_occurrences} is above the most allowed ({cap})")
        self._prepare_mems(min_len, len(query), threads)
        var out = List[Match]()
        for strand in range(2):
            if (strand == 0 and not forward) or (strand == 1 and not reverse):
                continue
            var qseqs = List[Span[UInt8, origin]]()
            qseqs.append(query)
            var qtext = SeqText(qseqs, reverse_complement=strand == 1)
            ref text = self._text.value()
            var found = find_mems(self._index.value(), text, qtext, threads)
            var kept = filter_by_occurrences(found^, qtext, cap, cap if mode == "mum" else 0, min_occurrences, threads)
            out.reserve(len(out) + len(kept))
            for m in kept:
                var seq = text.sequence_of(Int(m.ref_pos))
                out.append(Match(seq, Int(m.ref_pos) - text.starts[seq], Int(m.query_pos), Int(m.length), strand == 1))
        return out^

    def _prepare_mems(mut self, min_len: Int, query_len: Int, threads: Int) raises:
        """The packed text and a MEM index for `min_len` (any sampling chosen
        for one minimum length is exact for every query)."""
        if not self._text:
            var spans = List[Span[UInt8, origin_of(self.records.data)]]()
            for i in range(len(self.records)):
                spans.append(self.records.seq(i))
            self._text = SeqText(spans)
        if not self._index or self._index_min_len != min_len:
            var params = MemParams.choose(min_len, len(self._text.value()), query_len)
            params.validate()
            self._index = MemIndex(self._text.value(), params, threads)
            self._index_min_len = min_len

    # ---------------------------------------------------------------- SMEMs

    def smems[
        mut: Bool, //, origin: Origin[mut=mut]
    ](
        mut self,
        queries: List[Span[UInt8, origin]],
        min_len: Int,
        max_hits: Int = 20,
        min_occurrences: Int = 1,
        max_occurrences: Int = 0,
        threads: Int = 1,
    ) raises -> Smems:
        """The SMEMs of length >= `min_len` of every query, both strands, from
        an FM-index of the reference (built on first use). See
        `FMIndex.smems`; hits give reference sequence indices and 0-based
        positions."""
        self._prepare_fm(threads)
        return self._fm.value().smems(queries, min_len, max_hits, min_occurrences, max_occurrences, threads)

    def save_fm_index(mut self, path: String, threads: Int = 1) raises:
        """Writes the reference's FM-index (built if need be) for `FMIndex.load`
        and `sporos smem`."""
        self._prepare_fm(threads)
        self._fm.value().save(path)

    def _prepare_fm(mut self, threads: Int) raises:
        if not self._fm:
            var spans = List[Span[UInt8, origin_of(self.records.data)]]()
            for i in range(len(self.records)):
                spans.append(self.records.seq(i))
            self._fm = FMIndex.build(spans, self.records.names, threads)
