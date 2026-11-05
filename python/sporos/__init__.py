"""Fast DNA seeding from Python: sporos's Mojo library.

>>> import sporos
>>> list(sporos.minimizers("ACGTACGTAC", k=5, w=1))
[0, 1, 2, 3, 4, 5]
>>> list(sporos.syncmers("ACGTNACGTAC", k=4, s=4))
[0, 5, 6, 7]
>>> index = sporos.Index(["TTTTTTTTACGTACCGTAGGCATCAATTTTTTTT", "GGGGGGGCCGATTACAGGATCAGGGGGGGG"], min_len=10)
>>> index.matches("ACGTACCGTAGGCATCAA")
[(0, 8, 0, 18)]
>>> index.matches("GATCCTGTAATCGG", reverse=True)
[(1, 7, 0, 14)]
>>> fm = sporos.FMIndex(["TTTTTTTTACGTACCGTAGGCATCAATTTTTTTT"], names=["chrA"])
>>> fm.smems(["GGACGTACCGTAGGCATCAAGG"], min_len=10)
[[(2, 20, 1, [(0, 8, False)])]]

Sequences are `str` or `bytes` of ASCII letters: A, C, G, T and U in either case are bases, anything
else (N, say) a base no seed covers and no match crosses. Positions are 0-based. Seeds come back as an
`array('I')` of positions, matches as a list of `(ref_index, ref_pos, query_pos, length)` tuples,
SMEMs per query as `(start, end, count, hits)` with hits `(ref_index, ref_pos, reverse)`, and seed
hits per read as `(query_pos, ref_index, ref_pos, occurrences, reverse)`. `max_occurrences` and
`min_occurrences` count a string's occurrences in the reference throughout. `threads` of zero or less
uses every physical core.
"""

from __future__ import annotations

from array import array
from typing import Iterable, List, Optional, Tuple, Union

from . import _sporos

__all__ = ["FMIndex", "Index", "SeedIndex", "minimizers", "syncmers"]

Text = Union[str, bytes, bytearray, memoryview]
"""A sequence as `str`, or as ASCII bytes."""

Match = Tuple[int, int, int, int]
"""`(ref_index, ref_pos, query_pos, length)`: `length` bases of reference `ref_index` from `ref_pos`
equal the query's from `query_pos`."""

Smem = Tuple[int, int, int, List[Tuple[int, int, bool]]]
"""`(start, end, count, hits)`: query bases `[start, end)` occur `count` times in the reference, at
`hits` `(ref_index, ref_pos, reverse)` (none listed past `max_hits`)."""

SeedHit = Tuple[int, int, int, int, bool]
"""`(query_pos, ref_index, ref_pos, occurrences, reverse)`: the read's k-mer at `query_pos` occurs at
`ref_pos` of reference `ref_index` (on its reverse strand with `reverse`), one of `occurrences`."""

# The array type code of a native `uint32_t`, which the extension's bytes hold.
_UINT32 = "I" if array("I").itemsize == 4 else "L"

_MODES = {"maxmatch": 0, "mum": 1, "mumreference": 2}

_SCHEMES = {"minimizers": 0, "syncmers": 1, "kmers": 2}


def _text(seq: Text) -> str:
    """`seq` as the ASCII `str` the extension reads in place; `TypeError` unless text or bytes,
    `ValueError` unless ASCII, as positions count bytes."""
    if isinstance(seq, (bytes, bytearray, memoryview)):
        return bytes(seq).decode("ascii")
    if not isinstance(seq, str):
        raise TypeError(f"a sequence is a str or bytes, not {type(seq).__name__}")
    if not seq.isascii():
        raise ValueError("a sequence holds ASCII letters alone")
    return seq


def _int(value, name: str) -> int:
    """`value`, refused with `TypeError` unless an `int`."""
    if not isinstance(value, int) or isinstance(value, bool):
        raise TypeError(f"{name} is an int, not {type(value).__name__}")
    return value


def _texts(seqs: Iterable[Text], what: str) -> List[str]:
    """Many sequences; `TypeError` for a single one, which would iterate as letters."""
    if isinstance(seqs, (str, bytes, bytearray, memoryview)):
        raise TypeError(f"{what} is a list of sequences, not one sequence")
    return [_text(seq) for seq in seqs]


def _offsets_and_fields(found: bytes, count: int) -> Tuple[array, array]:
    """The extension's `uint64` answer: `count + 1` offsets, then the records' fields."""
    fields = array("Q")
    fields.frombytes(found)
    return fields[: count + 1], fields[count + 1 :]


def _positions(found: bytes) -> array:
    out = array(_UINT32)
    out.frombytes(found)
    return out


def _call(function, *args):
    """`function(*args)` with the library's refusals as `ValueError`."""
    try:
        return function(*args)
    except (TypeError, ValueError, MemoryError):
        raise
    except Exception as error:
        raise ValueError(str(error)) from None


def minimizers(seq: Text, k: int, w: int, canonical: bool = True, threads: int = 1) -> array:
    """The start positions of the (w, k)-minimizers of `seq`, in increasing order, each once per run
    of windows sharing it; canonical (strand-independent) unless `canonical` is false. `k` lies in
    [1, 64], `w` is at least 1; a sequence of fewer than `w` k-mers is one window."""
    seq = _text(seq)
    found = _call(
        _sporos.minimizers, seq, _int(k, "k"), _int(w, "w"), bool(canonical), _int(threads, "threads")
    )
    return _positions(found)


def syncmers(
    seq: Text,
    k: int,
    s: int,
    open_offset: Union[int, Iterable[int], None] = None,
    canonical: bool = True,
    threads: int = 1,
) -> array:
    """The start positions of the syncmers of `seq`, in increasing order: k-mers whose smallest s-mer
    starts at offset 0 or `k - s` (closed), or with `open_offset` at that offset, or any of a set of
    offsets (open). `1 <= s <= k` and every offset lies in `[0, k - s]`."""
    seq = _text(seq)
    if open_offset is None:
        offsets = []
    elif isinstance(open_offset, int) and not isinstance(open_offset, bool):
        offsets = [open_offset]
    else:
        offsets = [_int(t, "open_offset") for t in open_offset]
        if not offsets:
            raise ValueError("open_offset needs at least one offset")
    found = _call(
        _sporos.syncmers, seq, _int(k, "k"), _int(s, "s"), offsets, bool(canonical), _int(threads, "threads")
    )
    return _positions(found)


class Index:
    """An index of reference sequences for exact matches of at least `min_len` bases: every maximal
    exact match (MEM), or the maximal unique ones (MUMs), of a query against all of them at once.
    The sequences are copied in; fewer than 2^32 bases in all."""

    def __init__(self, seqs: Iterable[Text], min_len: int, threads: int = 1):
        texts = [_text(seq) for seq in seqs]
        if not texts:
            raise ValueError("an index needs at least one reference sequence")
        self.min_len = _int(min_len, "min_len")
        self._index = _call(_sporos.Index, texts, self.min_len, _int(threads, "threads"))

    def __len__(self) -> int:
        """How many reference sequences the index holds."""
        return self._index.sequences()

    def __repr__(self) -> str:
        return f"sporos.Index({len(self)} sequences, min_len={self.min_len})"

    def matches(
        self,
        query: Text,
        mode: str = "mum",
        reverse: bool = False,
        threads: int = 1,
        max_occurrences: int = 0,
        min_occurrences: int = 1,
    ) -> List[Match]:
        """The matches of `query` as `(ref_index, ref_pos, query_pos, length)`, sorted by query
        position, then reference: `mode` "maxmatch" every MEM, "mum" the MEMs unique in the references
        and the query, "mumreference" those unique in the references (MUMmer's `-maxmatch`, `-mum`,
        `-mumreference`). `max_occurrences` N relaxes "unique" to "at most N times" (Mumemto's `-F`;
        for "maxmatch" it caps the reference count) and `min_occurrences` keeps matches whose string
        occurs at least that often in the references. With `reverse` the query's reverse complement
        is matched instead, and query positions count along it, as MUMmer reports reverse matches."""
        if mode not in _MODES:
            raise ValueError(f"mode is one of {', '.join(_MODES)}, not {mode!r}")
        query = _text(query)
        limits = [_int(max_occurrences, "max_occurrences"), _int(min_occurrences, "min_occurrences")]
        found = _call(self._index.matches, query, _MODES[mode], bool(reverse), limits, _int(threads, "threads"))
        fields = iter(_positions(found))
        return list(zip(fields, fields, fields, fields))


class FMIndex:
    """An FM-index of reference sequences for the super-maximal exact matches (SMEMs) of many reads,
    as `bwa fastmap` and minibwa find them, and f-SMEMs. Build one from sequences (copied in), named
    `names` or seq0, seq1, ..., or `FMIndex.load` one saved by `save` or `sporos index`."""

    def __init__(self, seqs: Iterable[Text], names: Optional[Iterable[str]] = None, threads: int = 1):
        texts = _texts(seqs, "seqs")
        if not texts:
            raise ValueError("an index needs at least one reference sequence")
        named = [str(name) for name in names] if names is not None else []
        if len(named) not in (0, len(texts)):
            raise ValueError(f"{len(named)} names for {len(texts)} sequences")
        named = named or [f"seq{i}" for i in range(len(texts))]
        self._index = _call(_sporos.FMIndex, texts, named, _int(threads, "threads"))

    @classmethod
    def load(cls, path: str, threads: int = 1) -> "FMIndex":
        """The index saved at `path`."""
        index = cls.__new__(cls)
        index._index = _call(_sporos.FMIndex, str(path), _int(threads, "threads"))
        return index

    def save(self, path: str) -> None:
        """Writes the index to `path`, for `FMIndex.load` and `sporos smem`."""
        _call(self._index.save, str(path))

    def __len__(self) -> int:
        """How many reference sequences the index holds."""
        return self._index.sequences()

    @property
    def names(self) -> List[str]:
        """The reference sequences' names."""
        return [self._index.name(i) for i in range(len(self))]

    def __repr__(self) -> str:
        return f"sporos.FMIndex({len(self)} sequences)"

    def smems(
        self,
        queries: Iterable[Text],
        min_len: int,
        max_hits: int = 20,
        min_occurrences: int = 1,
        max_occurrences: int = 0,
        threads: int = 1,
    ) -> List[List[Smem]]:
        """Each query's SMEMs of at least `min_len` bases, both strands, in query order: maximal query
        intervals that occur in the reference and lie in no longer one that does, as
        `(start, end, count, hits)`. Hits are listed for SMEMs occurring at most `max_hits` times.
        With `min_occurrences` f > 1 they are f-SMEMs, the maximal intervals occurring at least f times
        (minibwa's `-s`); `max_occurrences` > 0 drops SMEMs occurring more often."""
        texts = _texts(queries, "queries")
        settings = [
            _int(min_len, "min_len"),
            _int(max_hits, "max_hits"),
            _int(min_occurrences, "min_occurrences"),
            _int(max_occurrences, "max_occurrences"),
            _int(threads, "threads"),
        ]
        if settings[0] < 1 or settings[1] < 0 or settings[2] < 1 or settings[3] < 0:
            raise ValueError("min_len and min_occurrences are at least 1, max_hits and max_occurrences at least 0")
        offsets, fields = _offsets_and_fields(_call(self._index.smems, texts, settings), len(texts))
        total = offsets[-1]
        records, hit_fields = fields[: 4 * total], fields[4 * total :]
        hits = [
            (hit_fields[i], hit_fields[i + 1], bool(hit_fields[i + 2])) for i in range(0, len(hit_fields), 3)
        ]
        smems = []
        for i in range(total):
            start, end, count, first = records[4 * i : 4 * i + 4]
            last = records[4 * i + 7] if i + 1 < total else len(hits)
            smems.append((start, end, count, hits[first:last]))
        return [smems[offsets[q] : offsets[q + 1]] for q in range(len(texts))]


class SeedIndex:
    """A hash index of reference k-mers (k <= 64) for seed lookups: the (w, k)-minimizers
    (`scheme` "minimizers"), the syncmers with s-mer length `w` ("syncmers": closed, or open at
    `open_offset`, an offset or a list of them), or every k-mer ("kmers"). The sequences are copied
    in; fewer than 2^32 bases in all."""

    def __init__(
        self,
        seqs: Iterable[Text],
        k: int,
        w: int,
        scheme: str = "minimizers",
        open_offset: Union[int, Iterable[int], None] = None,
        threads: int = 1,
    ):
        if scheme not in _SCHEMES:
            raise ValueError(f"scheme is one of {', '.join(_SCHEMES)}, not {scheme!r}")
        texts = _texts(seqs, "seqs")
        if not texts:
            raise ValueError("an index needs at least one reference sequence")
        if open_offset is None:
            offsets = []
        elif isinstance(open_offset, int) and not isinstance(open_offset, bool):
            offsets = [open_offset]
        else:
            offsets = [_int(t, "open_offset") for t in open_offset]
        self.k, self.w, self.scheme = _int(k, "k"), _int(w, "w"), scheme
        self._index = _call(
            _sporos.SeedIndex, texts, self.k, self.w, _SCHEMES[scheme], offsets, _int(threads, "threads")
        )

    def __repr__(self) -> str:
        return f"sporos.SeedIndex(k={self.k}, w={self.w}, scheme={self.scheme!r})"

    def occurrence_cutoff(self, keep_fraction: float) -> int:
        """The smallest occurrence count that `keep_fraction` of the indexed k-mers stay within: a
        `max_occurrences` for `lookup`, like GraphAligner's quantile cutoff."""
        return self._index.occurrence_cutoff(float(keep_fraction))

    def lookup(
        self, reads: Iterable[Text], max_occurrences: int = 0, all_kmers: bool = False, threads: int = 1
    ) -> List[List[SeedHit]]:
        """Each read's seed hits, in read order, as `(query_pos, ref_index, ref_pos, occurrences,
        reverse)`: its seeds under the index's scheme (every k-mer with `all_kmers`), skipping k-mers
        that occur more than `max_occurrences` times (0 for no limit)."""
        texts = _texts(reads, "reads")
        limit = _int(max_occurrences, "max_occurrences")
        if limit < 0:
            raise ValueError("max_occurrences is at least 0")
        found = _call(self._index.lookup, texts, bool(all_kmers), limit, _int(threads, "threads"))
        offsets, fields = _offsets_and_fields(found, len(texts))
        hits = [
            (fields[i], fields[i + 1], fields[i + 2], fields[i + 3], bool(fields[i + 4]))
            for i in range(0, len(fields), 5)
        ]
        return [hits[offsets[r] : offsets[r + 1]] for r in range(len(texts))]
