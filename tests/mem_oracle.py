"""Checks `sporos mem` with occurrence limits against brute force.

    python3 tests/mem_oracle.py <sporos> <work-dir>

Every MEM of length >= L between the query and the reference (forward strand; `-r` on the reverse
complement), found by seeding and extending; then each MEM's string counted in the reference and in the
query by plain string search. `--max-occ N` keeps MEMs whose string occurs at most N times in the
reference (and, with `-mum`, at most N times in the query); `--min-occ M` those occurring at least M times
in the reference. Inputs: repeats of several copy numbers, tandem repeats, homopolymers, N runs, two
reference sequences.
"""

import random
import subprocess
import sys

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def mems(ref, query, min_len):
    """(ref_pos, query_pos, length), 0-based; `ref` uses '#' between sequences, N never matches."""
    k = min(min_len, 8)
    index = {}
    for r in range(len(ref) - k + 1):
        kmer = ref[r : r + k]
        if "N" in kmer or "#" in kmer:
            continue
        index.setdefault(kmer, []).append(r)
    out = []
    for q in range(len(query) - k + 1):
        for r in index.get(query[q : q + k], ()):
            if q > 0 and r > 0 and ref[r - 1] == query[q - 1] and query[q - 1] in "ACGT":
                continue  # not left-maximal
            n = k
            while r + n < len(ref) and q + n < len(query) and ref[r + n] == query[q + n] and ref[r + n] in "ACGT":
                n += 1
            if n >= min_len:
                out.append((r, q, n))
    return out


def count(text, s):
    n, at = 0, text.find(s)
    while at >= 0:
        n += 1
        at = text.find(s, at + 1)
    return n


def expected(ref, query, min_len, mode, max_occ, min_occ):
    kept = []
    for r, q, n in mems(ref, query, min_len):
        s = query[q : q + n]
        in_ref = count(ref, s)
        cap = max_occ if max_occ else (1 if mode != "maxmatch" else 0)
        if cap and in_ref > cap:
            continue
        if in_ref < min_occ:
            continue
        if mode == "mum" and count(query, s) > max(cap, 1):
            continue
        kept.append((r, q, n))
    return sorted(kept)


def parse(text, offsets):
    """Matches as 0-based (text offset of the reference position, query position, length)."""
    out = []
    for line in text.splitlines():
        if line.startswith(">") or not line.strip():
            continue
        f = line.split()
        if len(f) == 4:  # name ref_pos query_pos length
            out.append((offsets[f[0]] + int(f[1]) - 1, int(f[2]) - 1, int(f[3])))
        else:
            out.append((int(f[0]) - 1, int(f[1]) - 1, int(f[2])))
    return sorted(out)


def main():
    sporos, work = sys.argv[1], sys.argv[2]
    rnd = random.Random(31)

    def dna(n):
        return "".join(rnd.choice("ACGT") for _ in range(n))

    def mutate(s, rate):
        return "".join(rnd.choice("ACGT") if rnd.random() < rate else c for c in s)

    unit = dna(60)
    copies = "".join(mutate(unit, 0.01) + dna(40) for _ in range(6))  # a repeat in ~6 copies
    ref1 = dna(1500) + copies + "A" * 80 + "ACAC" * 30 + "N" * 20 + dna(800) + unit + dna(300)
    ref2 = dna(700) + unit + unit + dna(200)
    ref = ref1 + "#" + ref2
    query = mutate(ref1[200:1400], 0.02) + unit + dna(100) + unit + copies[:300] + "A" * 60 + mutate(ref2[100:600], 0.01)
    # Reverse-complement pieces, for -r: a repeat's copies and unique sequence.
    query += revcomp(copies[100:500]) + dna(50) + revcomp(mutate(ref1[2400:3000], 0.02)) + revcomp(unit)
    with open(f"{work}/occ_ref.fa", "w") as f:
        f.write(f">r1\n{ref1}\n>r2\n{ref2}\n")
    with open(f"{work}/occ_query.fa", "w") as f:
        f.write(f">q\n{query}\n")
    offsets = {"r1": 0, "r2": len(ref1) + 1}
    failures = 0
    for strand, q in (("", query), ("-r", revcomp(query))):
        for min_len in (12, 20, 40):
            for mode in ("maxmatch", "mumreference", "mum"):
                for max_occ, min_occ in ((0, 1), (1, 1), (2, 1), (5, 1), (100, 1), (0, 2), (5, 3)):
                    if mode == "maxmatch" and max_occ == 1:
                        continue
                    args = [sporos, "mem", f"-{mode}", "-l", str(min_len), "-t", "3", "-F"]
                    if strand:
                        args.append(strand)
                    if max_occ:
                        args += ["--max-occ", str(max_occ)]
                    if min_occ > 1:
                        args += ["--min-occ", str(min_occ)]
                    out = subprocess.run(
                        args + [f"{work}/occ_ref.fa", f"{work}/occ_query.fa"], check=True, capture_output=True, text=True
                    ).stdout
                    got = parse(out, offsets)
                    want = expected(ref, q, min_len, mode, max_occ, min_occ)
                    label = f"{mode} {strand or 'fwd'} L={min_len} --max-occ {max_occ} --min-occ {min_occ}"
                    if got != want:
                        failures += 1
                        print(f"  FAIL  {label}: {len(got)} matches, expected {len(want)}")
                        print("        only got:", sorted(set(got) - set(want))[:5])
                        print("        only expected:", sorted(set(want) - set(got))[:5])
                    else:
                        print(f"  ok    {label}: {len(got)} matches")
    if failures:
        print(f"mem_oracle: {failures} check(s) failed")
        sys.exit(1)
    print("mem_oracle: all checks passed")


if __name__ == "__main__":
    main()
