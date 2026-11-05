"""The Python package (python/sporos) against its own docstring's examples, what the Mojo library
answers on the same inputs (seeds, matches, SMEMs and seed-index hits; tests/api_expected.mojo writes
those), and a brute-force match oracle sharing no code with the library.

    pixi run test-python
"""

import doctest
import os
import random
import sys
import tempfile

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


def oracle(refs, query, min_len, mode, reverse, max_occurrences=0, min_occurrences=1):
    """Every MEM of at least `min_len` by trying every pair of starts, kept by `mode` and the occurrence
    limits from how often its string occurs in the references and the query, sorted as the library
    sorts them."""
    cap = max_occurrences if max_occurrences > 0 else (0 if mode == "maxmatch" else 1)
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
                if cap > 0 or min_occurrences > 1:
                    in_refs = sum(occurrences(other, query, j, length) for other in refs)
                    if (cap > 0 and in_refs > cap) or in_refs < min_occurrences:
                        continue
                    if mode == "mum" and occurrences(query, query, j, length) > cap:
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
    seq = "ACGTTGCAGGCATTACGATCGATCCGATAGCTAGCTTAGCAGTCGA" * 5
    check(sporos.syncmers(seq, 9, 4, open_offset=[0, 5]) == sporos.syncmers(seq, 9, 4), "offsets {0, k - s} are closed")
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
        (lambda: sporos.syncmers("ACGTACGT", 5, 3, open_offset=[0, 3]), ValueError),
        (lambda: sporos.syncmers("ACGTACGT", 5, 3, open_offset=[]), ValueError),
        (lambda: sporos.syncmers("ACGTACGT", 5, 3, open_offset=[1.0]), TypeError),
        (lambda: sporos.Index([], 10), ValueError),
        (lambda: sporos.Index(["ACGT"], 0), ValueError),
        (lambda: sporos.Index(["ACGTACGTACGT"], 5).matches("ACGT", mode="mems"), ValueError),
        (lambda: sporos.Index(["ACGTACGTACGT"], 5).matches("ACGT", min_occurrences=2), ValueError),
        (lambda: sporos.Index(["ACGTACGTACGT"], 5).matches("ACGT", mode="maxmatch", min_occurrences=0), ValueError),
        (lambda: sporos.FMIndex([]), ValueError),
        (lambda: sporos.FMIndex("ACGT"), TypeError),
        (lambda: sporos.FMIndex(["ACGT"], names=["a", "b"]), ValueError),
        (lambda: sporos.FMIndex(["ACGTACGT"]).smems("ACGT", 3), TypeError),
        (lambda: sporos.FMIndex(["ACGTACGT"]).smems(["ACGT"], 0), ValueError),
        (lambda: sporos.FMIndex.load("/nonexistent/dir/x.sfm"), ValueError),
        (lambda: sporos.SeedIndex(["ACGTACGT"], 65, 1, scheme="kmers"), ValueError),
        (lambda: sporos.SeedIndex(["ACGTACGT"], 5, 1, scheme="strobemers"), ValueError),
        (lambda: sporos.SeedIndex(["ACGTACGT"], 5, 1).lookup(["ACGT"], max_occurrences=-1), ValueError),
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
    check(len(index.matches(query + "CC" + query, max_occurrences=2)) == 2, "query cap of 2")
    check(index.matches(query, mode="maxmatch", min_occurrences=2) == [], "floor of 2")

    fm = sporos.FMIndex(refs, names=["chrA", "chrB"])
    check(len(fm) == 2 and fm.names == ["chrA", "chrB"], "FMIndex names")
    check(sporos.FMIndex(refs).names == ["seq0", "seq1"], "FMIndex default names")
    check(fm.smems(["GG" + query + "GG"], 10) == [[(2, 20, 1, [(0, 8, False)])]], "SMEM")
    check(fm.smems([reverse_complement(query)], 10) == [[(0, 18, 1, [(0, 8, True)])]], "reverse SMEM")
    check(fm.smems([query], 10, min_occurrences=2) == [[]], "nothing occurs twice")
    check(fm.smems([], 10) == [], "no queries")
    seeds = sporos.SeedIndex(refs[:1], 5, 1, scheme="kmers")
    # ACGTA also occurs as TACGT at 7.
    expected = [(0, 0, 8, 2, False), (0, 0, 7, 2, True), (1, 0, 9, 1, False), (2, 0, 10, 1, False)]
    check(seeds.lookup(["ACGTACC"]) == [expected], "seed hits")
    check(seeds.occurrence_cutoff(1.0) >= 1, "cutoff")


def test_oracle() -> None:
    rng = random.Random(7)
    for round in range(4):
        refs = [random_dna(160, rng), random_dna(120, rng)]
        refs[1] = refs[1][:60] + "N" + refs[1][61:90] + refs[0][15:45]
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
                for cap, floor in ((2, 1), (3, 2)):
                    got = index.matches(query, mode, reverse, 1, max_occurrences=cap, min_occurrences=floor)
                    expected = oracle(refs, query, min_len, mode, reverse, cap, floor)
                    check(got == expected, f"oracle round {round} {mode} reverse={reverse} limits {cap} {floor}")


def test_against_mojo(path: str) -> None:
    seqs = {}
    indexes = {}
    records = 0
    fm = None
    with open(path) as file:
        for line in file:
            fields = line.split()
            kind, name = fields[0], fields[1]
            if kind == "seq":
                seqs[name] = fields[2]
                continue
            records += 1
            refs = [seqs["ref0"], seqs["ref1"]] if "ref1" in seqs else []
            if kind in ("smems", "seeds"):
                numbers = iter(int(field) for field in fields[2:])
                reads = [seqs[f"{name}{r}"] for r in range(next(numbers))]
                if kind == "smems":
                    min_len, max_hits, min_occ, max_occ, threads = (next(numbers) for _ in range(5))
                    if fm is None:
                        fm = sporos.FMIndex(refs, names=["ref0", "ref1"], threads=2)
                    got = fm.smems(reads, min_len, max_hits, min_occ, max_occ, threads)
                    expected = []
                    for _ in reads:
                        smems = []
                        for _ in range(next(numbers)):
                            start, end, count, listed = (next(numbers) for _ in range(4))
                            hits = [(next(numbers), next(numbers), bool(next(numbers))) for _ in range(listed)]
                            smems.append((start, end, count, hits))
                        expected.append(smems)
                else:
                    k, w, scheme, open_offset, max_occ, all_kmers, threads = (next(numbers) for _ in range(7))
                    index = sporos.SeedIndex(
                        refs,
                        k,
                        w,
                        scheme=("minimizers", "syncmers", "kmers")[scheme],
                        open_offset=None if open_offset < 0 else open_offset,
                        threads=threads,
                    )
                    got = index.lookup(reads, max_occurrences=max_occ, all_kmers=bool(all_kmers), threads=threads)
                    expected = [
                        [
                            tuple(next(numbers) for _ in range(4)) + (bool(next(numbers)),)
                            for _ in range(next(numbers))
                        ]
                        for _ in reads
                    ]
                check(got == expected, line[:60])
                continue
            if kind == "cutoff":
                k, w, scheme, open_offset = int(name), int(fields[2]), int(fields[3]), int(fields[4])
                index = sporos.SeedIndex(
                    refs,
                    k,
                    w,
                    scheme=("minimizers", "syncmers", "kmers")[scheme],
                    open_offset=None if open_offset < 0 else open_offset,
                )
                check(index.occurrence_cutoff(float(fields[5])) == int(fields[6]), line[:60])
                continue
            # Syncmers' open offsets (field 4) are a comma-separated list.
            numbers = [int(field) for field in (fields[2:4] + fields[5:] if kind == "syncmers" else fields[2:])]
            seq = seqs[name]
            if kind == "minimizers":
                k, w, canonical, threads, count = numbers[:5]
                got = sporos.minimizers(seq, k, w, canonical=bool(canonical), threads=threads)
                check(list(got) == numbers[5:] and len(got) == count, line[:60])
                check(got == sporos.minimizers(seq.encode(), k, w, bool(canonical), threads), "bytes " + line[:60])
            elif kind == "syncmers":
                offsets = [int(t) for t in fields[4].split(",")]
                k, s, canonical, threads, count = numbers[:5]
                opened = None if offsets == [-1] else offsets[0] if len(offsets) == 1 else offsets
                got = sporos.syncmers(seq, k, s, open_offset=opened, canonical=bool(canonical), threads=threads)
                check(list(got) == numbers[5:] and len(got) == count, line[:60])
            else:
                min_len, mode, reverse, max_occ, min_occ, threads, count = numbers[:7]
                if min_len not in indexes:
                    indexes[min_len] = sporos.Index(refs, min_len, threads=threads)
                flat = numbers[7:]
                expected = [tuple(flat[i : i + 4]) for i in range(0, len(flat), 4)]
                got = indexes[min_len].matches(
                    seq,
                    mode=("maxmatch", "mum", "mumreference")[mode],
                    reverse=bool(reverse),
                    threads=threads,
                    max_occurrences=max_occ,
                    min_occurrences=min_occ,
                )
                check(got == expected and len(got) == count, line[:60])
    check(records > 0, "expected records")

    # A saved FM-index loads with its names and the same SMEMs.
    check(fm is not None, "smems records")
    if fm is None:
        return
    reads = [seqs[f"read{r}"] for r in range(sum(name.startswith("read") for name in seqs))]
    with tempfile.TemporaryDirectory() as directory:
        saved = os.path.join(directory, "index.sfm")
        fm.save(saved)
        loaded = sporos.FMIndex.load(saved, threads=2)
    check(loaded.names == ["ref0", "ref1"], "loaded names")
    check(loaded.smems(reads, 19, threads=3) == fm.smems(reads, 19), "loaded SMEMs")


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
