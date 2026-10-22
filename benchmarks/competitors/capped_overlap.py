"""Forward-strand overlap of Mumemto's capped multi-MEMs with sporos's capped MEMs between two genomes.

    python3 -I capped_overlap.py <mumemto.mems> <sporos mem output>

Mumemto prints each matched string once with all its copies (0-based starts, the genome of each, strands);
sporos prints MUMmer-style position pairs (1-based). Every forward copy in the reference (genome 0) paired
with every forward copy in the query (genome 1) gives Mumemto's pairs. Prints, tab-separated: Mumemto's
lines, its forward pairs, sporos's forward pairs, and the pairs in both.
"""
import sys

mm, lines = set(), 0
for line in open(sys.argv[1]):
    f = line.split("\t")
    if len(f) < 4:
        continue
    lines += 1
    length, starts, genomes, strands = int(f[0]), f[1].split(","), f[2].split(","), f[3].strip().split(",")
    ref = [int(s) + 1 for s, g, d in zip(starts, genomes, strands) if g == "0" and d == "+"]
    qry = [int(s) + 1 for s, g, d in zip(starts, genomes, strands) if g == "1" and d == "+"]
    mm.update((r, q, length) for r in ref for q in qry)

sp, forward = set(), False
for line in open(sys.argv[2]):
    if line.startswith(">"):
        forward = "Reverse" not in line
        continue
    f = line.split()
    if forward and len(f) >= 3:
        sp.add((int(f[-3]), int(f[-2]), int(f[-1])))

print(lines, len(mm), len(sp), len(mm & sp), sep="\t")
