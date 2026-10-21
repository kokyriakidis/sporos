"""Checks `sporos smem` against SMEMs computed from their definition by plain string search.

    python3 tests/smem_oracle.py <sporos> <work-dir>

An SMEM of length >= L is a query interval [s, e) that occurs in the reference on either strand (never
across a non-ACGT base or a sequence end) while neither [s - 1, e) nor [s, e + 1) does; its hits are all
its occurrences, forward and reverse-complement. The references and reads are built to be hard: homopolymers,
tandem repeats, palindromes, N and IUPAC runs, lowercase, one-base and empty sequences, duplicate records,
reads across sequence ends and N runs, reads shorter than L. Both engines, several L and suffix-array
sampling rates are compared on every SMEM, its count and every hit.
"""

import random
import subprocess
import sys

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def runs_of(seq):
    """(offset, run) for each maximal ACGT run of `seq` (case-insensitive, U as T)."""
    seq = seq.upper().replace("U", "T")
    out, i = [], 0
    while i < len(seq):
        if seq[i] not in "ACGT":
            i += 1
            continue
        j = i
        while j < len(seq) and seq[j] in "ACGT":
            j += 1
        out.append((i, seq[i:j]))
        i = j
    return out


class Reference:
    def __init__(self, records):
        self.runs = []  # (name, offset, run)
        for name, seq in records:
            for off, run in runs_of(seq):
                self.runs.append((name, off, run))
        self.text = "#" + "#".join(r[2] for r in self.runs) + "#"

    def occurs(self, sub, at_least=1):
        if at_least == 1:
            return sub in self.text or revcomp(sub) in self.text
        return self.count(sub) >= at_least

    def count(self, sub):
        n = 0
        for pattern in (sub, revcomp(sub)):
            at = self.text.find(pattern)
            while at >= 0:
                n += 1
                at = self.text.find(pattern, at + 1)
        return n

    def hits(self, sub):
        out = []
        for strand, pattern in (("+", sub), ("-", revcomp(sub))):
            for name, off, run in self.runs:
                at = run.find(pattern)
                while at >= 0:
                    out.append(f"{name}:{strand}{off + at + 1}")
                    at = run.find(pattern, at + 1)
        return out


def oracle(ref, query, min_len, at_least=1):
    """SMEMs, or with `at_least` f > 1 f-SMEMs: maximal intervals occurring at least f times."""
    q = query.upper().replace("U", "T")
    n = len(q)
    ends = [0] * n
    e = 0
    for s in range(n):
        e = max(e, s)
        while e < n and ref.occurs(q[s : e + 1], at_least):
            e += 1
        ends[s] = e
    smems = []
    for s in range(n):
        if ends[s] - s >= min_len and (s == 0 or ends[s - 1] < ends[s]):
            smems.append((s, ends[s]))
    return smems


def parse(text):
    """{read: [(start, end, count, sorted hits)]} from fastmap-style output."""
    out, read = {}, None
    for line in text.splitlines():
        f = line.split("\t")
        if f[0] == "SQ":
            read = f[1]
            out[read] = []
        elif f[0] == "EM":
            out[read].append((int(f[1]), int(f[2]), int(f[3]), tuple(sorted(f[4:]))))
    return out


def dna(rnd, n):
    return "".join(rnd.choice("ACGT") for _ in range(n))


def mutate(rnd, s, rate):
    return "".join(rnd.choice("ACGT") if rnd.random() < rate else c for c in s)


def build(rnd, size):
    unit = dna(rnd, 7)
    pal = dna(rnd, 40)
    blocks = [
        dna(rnd, size),
        "A" * 300,
        "ACACACACAC" * 30,
        mutate(rnd, unit * 60, 0.02),
        pal + revcomp(pal),
        "NNNN",
        dna(rnd, 200).lower(),
        "N",
        dna(rnd, size // 2),
        "RYKMSW" + dna(rnd, 50) + "n" * 7,
    ]
    main = "".join(blocks)
    copy = mutate(rnd, main[1000 : 1000 + size // 3], 0.01)
    records = [
        ("chrA", main),
        ("chrB", revcomp(copy) + "N" * 30 + dna(rnd, 300)),
        ("tiny", "G"),
        ("empty", ""),
        ("allN", "N" * 40),
        ("dup1", main[200:700]),
        ("dup2", main[200:700]),
    ]
    reads = []
    flat = [(name, seq) for name, seq in records if seq]
    for i in range(400):
        name, seq = rnd.choice(flat)
        length = rnd.choice([1, 5, 18, 19, 31, 32, 33, 64, 100, 150, 600])
        at = rnd.randrange(-20, max(1, len(seq) - length + 20))
        piece = seq[max(0, at) : max(0, at) + length] or "A"
        piece = mutate(rnd, piece, rnd.choice([0, 0, 0.01, 0.05]))
        if rnd.random() < 0.5:
            piece = revcomp(piece.upper())
        if rnd.random() < 0.1:
            k = rnd.randrange(len(piece))
            piece = piece[:k] + "N" + piece[k + 1 :]
        reads.append(piece)
    reads += [pal + revcomp(pal), "A" * 120, "ACAC" * 40, unit * 20, "N" * 50, "acgtnACGT" * 9, records[5][1]]
    reads.append(revcomp(main[50:2050]))
    return records, reads


def write_fasta(path, records):
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n{seq}\n")


def main():
    sporos, work = sys.argv[1], sys.argv[2]
    failures = 0
    for case, size in ((1, 600), (2, 30_000)):
        rnd = random.Random(1000 + case)
        records, reads = build(rnd, size)
        ref_path, reads_path = f"{work}/oracle_ref{case}.fa", f"{work}/oracle_reads{case}.fa"
        write_fasta(ref_path, records)
        named = [(f"q{i}", r) for i, r in enumerate(reads)]
        write_fasta(reads_path, named)
        ref = Reference(records)
        indexes = []
        for bits in (1, 2, 4):
            path = f"{work}/oracle{case}_s{bits}.sfm"
            subprocess.run([sporos, "index", "-s", str(bits), "-t", "3", ref_path, path], check=True)
            indexes.append(path)
        for min_len in (1, 2, 5, 11, 19, 31, 32, 33, 64):
            expected = {}
            for name, r in named:
                expected[name] = [
                    (s, e, len(h := ref.hits(r.upper().replace("U", "T")[s:e])), tuple(sorted(h)))
                    for s, e in oracle(ref, r, min_len)
                ]
            runs = [("fm", ["--engine", "fm", ref_path]), ("mem", ["--engine", "mem", ref_path])]
            runs += [(f"index {p.rsplit('_', 1)[1]}", [p]) for p in indexes]
            for label, args in runs:
                out = subprocess.run(
                    [sporos, "smem", "-l", str(min_len), "-w", "100000000", "-t", "3"] + args[:-1] + [args[-1], reads_path],
                    check=True,
                    capture_output=True,
                    text=True,
                ).stdout
                got = parse(out)
                bad = [name for name, _ in named if got.get(name, []) != expected[name]]
                total = sum(len(v) for v in expected.values())
                if bad:
                    failures += 1
                    name = bad[0]
                    print(f"  FAIL  case {case} L={min_len} {label}: {len(bad)} reads differ, e.g. {name}")
                    print(f"        got      {got.get(name)}")
                    print(f"        expected {expected[name]}")
                else:
                    print(f"  ok    case {case} L={min_len} {label}: {total} SMEMs")
        # f-SMEMs (-s f: at least f occurrences) and capped SMEMs (--max-occ N: at most N).
        for min_len in (1, 5, 19, 33):
            for at_least, at_most in ((2, 0), (3, 0), (8, 0), (1, 5), (2, 10), (1, 1)):
                expected = {}
                for name, r in named:
                    rows = []
                    for s, e in oracle(ref, r, min_len, at_least):
                        h = ref.hits(r.upper().replace("U", "T")[s:e])
                        if at_most == 0 or len(h) <= at_most:
                            rows.append((s, e, len(h), tuple(sorted(h))))
                    expected[name] = rows
                for label, args in (("fm", [ref_path]), (f"index {indexes[0].rsplit('_', 1)[1]}", [indexes[0]])):
                    cmd = [sporos, "smem", "-l", str(min_len), "-w", "100000000", "-t", "3", "-s", str(at_least)]
                    if at_most:
                        cmd += ["--max-occ", str(at_most)]
                    out = subprocess.run(cmd + args + [reads_path], check=True, capture_output=True, text=True).stdout
                    got = parse(out)
                    bad = [name for name, _ in named if got.get(name, []) != expected[name]]
                    total = sum(len(v) for v in expected.values())
                    what = f"case {case} L={min_len} -s {at_least} --max-occ {at_most} {label}"
                    if bad:
                        failures += 1
                        name = bad[0]
                        print(f"  FAIL  {what}: {len(bad)} reads differ, e.g. {name}")
                        print(f"        got      {got.get(name)}")
                        print(f"        expected {expected[name]}")
                    else:
                        print(f"  ok    {what}: {total} SMEMs")
    if failures:
        print(f"smem_oracle: {failures} check(s) failed")
        sys.exit(1)
    print("smem_oracle: all checks passed")


if __name__ == "__main__":
    main()
