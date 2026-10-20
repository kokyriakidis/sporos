"""The Python package (python/sporos) against its own docstring's examples, what the Mojo library
answers on the same inputs (tests/api_expected.mojo writes those), and a brute-force match oracle
sharing no code with the library.

    pixi run test-python
"""

import doctest
import random
import sys

import sporos

failures = 0


def check(condition: bool, label: str) -> None:
    """Reports a failed check without stopping the run."""
    global failures
    if not condition:
        failures += 1
        print(f"failed: {label}", file=sys.stderr)


_CODES = {"a": 0, "c": 1, "g": 2, "t": 3, "u": 3}


def same(a: str, b: str) -> bool:
    code = _CODES.get(a.lower())
    return code is not None and code == _CODES.get(b.lower())


def reverse_complement(seq: str) -> str:
    return "".join("TGCA"[_CODES[c.lower()]] if c.lower() in _CODES else "N" for c in reversed(seq))


def occurrences(text: str, query: str, start: int, length: int) -> int:
    """How often `query[start:start + length]` occurs in `text`."""
    return sum(
        all(same(text[i + j], query[start + j]) for j in range(length)) for i in range(len(text) - length + 1)
    )


def oracle(refs, query, min_len, mode, reverse):
    """Every MEM of at least `min_len` by trying every pair of starts, kept by `mode` from how often its
    string occurs in the references and the query, sorted as the library sorts them."""
    if reverse:
        query = reverse_complement(query)
    out = []
    for r, ref in enumerate(refs):
        for i in range(len(ref)):
            for j in range(len(query)):
                if i > 0 and j > 0 and same(ref[i - 1], query[j - 1]):
                    continue
                length = 0
                while i + length < len(ref) and j + length < len(query) and same(ref[i + length], query[j + length]):
                    length += 1
                if length < min_len:
                    continue
                if mode != "maxmatch":
                    if sum(occurrences(other, query, j, length) for other in refs) != 1:
                        continue
                    if mode == "mum" and occurrences(query, query, j, length) != 1:
                        continue
                out.append((r, i, j, length))
    out.sort(key=lambda m: (m[2], m[0], m[1], m[3]))
    return out


def random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGTacgt") for _ in range(length))


def test_known() -> None:
    check(list(sporos.minimizers("ACGTNACGTAC", 4, 1, canonical=False)) == [0, 5, 6, 7], "N breaks k-mers")
    check(list(sporos.minimizers(b"ACGTACGTAC", 5, 1)) == [0, 1, 2, 3, 4, 5], "bytes input")
    check(list(sporos.minimizers("ACG", 5, 3)) == [], "shorter than k")
    check(list(sporos.minimizers("", 5, 3)) == [], "empty")
    check(list(sporos.syncmers("ACGTACGTAC", 5, 5)) == [0, 1, 2, 3, 4, 5], "s = k")
    check(len(sporos.syncmers("ACGTACGTACGTTGCA", 7, 3, open_offset=2)) >= 0, "open syncmers")
    for call, error in [
        (lambda: sporos.minimizers("ACGT", 0, 3), ValueError),
        (lambda: sporos.minimizers("ACGT", 65, 3), ValueError),
        (lambda: sporos.minimizers("ACGT", 5, 0), ValueError),
        (lambda: sporos.minimizers("ACGTé", 3, 3), ValueError),
        (lambda: sporos.minimizers(42, 3, 3), TypeError),
        (lambda: sporos.minimizers("ACGT", 3.0, 3), TypeError),
        (lambda: sporos.syncmers("ACGT", 3, 5), ValueError),
        (lambda: sporos.syncmers("ACGTACGT", 5, 3, open_offset=3), ValueError),
        (lambda: sporos.syncmers("ACGTACGT", 5, 3, open_offset=-1), ValueError),
        (lambda: sporos.Index([], 10), ValueError),
        (lambda: sporos.Index(["ACGT"], 0), ValueError),
        (lambda: sporos.Index(["ACGTACGTACGT"], 5).matches("ACGT", mode="mems"), ValueError),
    ]:
        try:
            call()
            check(False, f"no {error.__name__}")
        except error:
            pass

    refs = ["TTTTTTTTACGTACCGTAGGCATCAATTTTTTTT", b"GGGGGGGCCGATTACAGGATCAGGGGGGGG"]
    index = sporos.Index(refs, 10)
    check(len(index) == 2, "len")
    check(repr(index) == "sporos.Index(2 sequences, min_len=10)", "repr")
    query = "ACGTACCGTAGGCATCAA"
    check(index.matches(query + "CC" + query) == [], "repeat in the query is no MUM")
    check(len(index.matches(query + "CC" + query, mode="mumreference")) == 2, "mumreference")
    check(len(index.matches(query + "CC" + query, mode="maxmatch")) == 2, "maxmatch")
    check(index.matches("") == [], "empty query")


def test_oracle() -> None:
    rng = random.Random(7)
    for round in range(4):
        refs = [random_dna(160, rng), random_dna(120, rng)]
        refs[1] = refs[1][:60] + "N" + refs[1][61:]
        query = (
            refs[0][10:50]
            + random_dna(20, rng)
            + reverse_complement(refs[1][30:80])
            + refs[0][20:40]
            + refs[1][90:120]
        )
        min_len = 10 + 3 * round
        index = sporos.Index(refs, min_len, threads=1 + round % 2)
        for mode in ("maxmatch", "mum", "mumreference"):
            for reverse in (False, True):
                got = index.matches(query, mode=mode, reverse=reverse, threads=round % 3)
                expected = oracle(refs, query, min_len, mode, reverse)
                check(got == expected, f"oracle round {round} {mode} reverse={reverse}: {len(got)} vs {len(expected)}")


def test_against_mojo(path: str) -> None:
    seqs = {}
    indexes = {}
    records = 0
    with open(path) as file:
        for line in file:
            fields = line.split()
            kind, name = fields[0], fields[1]
            if kind == "seq":
                seqs[name] = fields[2]
                continue
            records += 1
            numbers = [int(field) for field in fields[2:]]
            seq = seqs[name]
            if kind == "minimizers":
                k, w, canonical, threads, count = numbers[:5]
                got = sporos.minimizers(seq, k, w, canonical=bool(canonical), threads=threads)
                check(list(got) == numbers[5:] and len(got) == count, line[:60])
                check(got == sporos.minimizers(seq.encode(), k, w, bool(canonical), threads), "bytes " + line[:60])
            elif kind == "syncmers":
                k, s, offset, canonical, threads, count = numbers[:6]
                got = sporos.syncmers(
                    seq, k, s, open_offset=None if offset < 0 else offset, canonical=bool(canonical), threads=threads
                )
                check(list(got) == numbers[6:] and len(got) == count, line[:60])
            else:
                min_len, mode, reverse, threads, count = numbers[:5]
                if min_len not in indexes:
                    indexes[min_len] = sporos.Index([seqs["ref0"], seqs["ref1"]], min_len, threads=threads)
                flat = numbers[5:]
                expected = [tuple(flat[i : i + 4]) for i in range(0, len(flat), 4)]
                got = indexes[min_len].matches(
                    seq, mode=("maxmatch", "mum", "mumreference")[mode], reverse=bool(reverse), threads=threads
                )
                check(got == expected and len(got) == count, line[:60])
    check(records > 0, "expected records")


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_python.py <expected.txt from tests/api_expected.mojo>", file=sys.stderr)
        return 2
    tested = doctest.testmod(sporos)
    check(tested.failed == 0 and tested.attempted > 0, "docstring examples")
    test_known()
    test_oracle()
    test_against_mojo(sys.argv[1])
    if failures:
        print(f"test_python: {failures} checks failed", file=sys.stderr)
        return 1
    print("test_python: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
