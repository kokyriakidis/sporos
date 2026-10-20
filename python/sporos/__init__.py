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

Sequences are `str` or `bytes` of ASCII letters: A, C, G, T and U in either case are bases, anything
else (N, say) a base no seed covers and no match crosses. Positions are 0-based. Seeds come back as an
`array('I')` of positions, matches as a list of `(ref_index, ref_pos, query_pos, length)` tuples.
`threads` of zero or less uses every physical core.
"""

from __future__ import annotations

from array import array
from typing import Iterable, List, Optional, Tuple, Union

from . import _sporos

__all__ = ["Index", "minimizers", "syncmers"]

Text = Union[str, bytes, bytearray, memoryview]
"""A sequence as `str`, or as ASCII bytes."""

Match = Tuple[int, int, int, int]
"""`(ref_index, ref_pos, query_pos, length)`: `length` bases of reference `ref_index` from `ref_pos`
equal the query's from `query_pos`."""

# The array type code of a native `uint32_t`, which the extension's bytes hold.
_UINT32 = "I" if array("I").itemsize == 4 else "L"

_MODES = {"maxmatch": 0, "mum": 1, "mumreference": 2}


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

    def matches(self, query: Text, mode: str = "mum", reverse: bool = False, threads: int = 1) -> List[Match]:
        """The matches of `query` as `(ref_index, ref_pos, query_pos, length)`, sorted by query
        position, then reference: `mode` "maxmatch" every MEM, "mum" the MEMs unique in the references
        and the query, "mumreference" those unique in the references (MUMmer's `-maxmatch`, `-mum`,
        `-mumreference`). With `reverse` the query's reverse complement is matched instead, and query
        positions count along it, as MUMmer reports reverse matches."""
        if mode not in _MODES:
            raise ValueError(f"mode is one of {', '.join(_MODES)}, not {mode!r}")
        query = _text(query)
        found = _call(self._index.matches, query, _MODES[mode], bool(reverse), _int(threads, "threads"))
        fields = iter(_positions(found))
        return list(zip(fields, fields, fields, fields))
