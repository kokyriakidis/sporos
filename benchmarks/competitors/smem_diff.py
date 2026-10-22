"""Compares two SMEM outputs in `bwa fastmap` format (SQ/EM/// lines) over a single-sequence reference.

    python3 -I smem_diff.py <ref.fa> <sporos.txt> <other.txt>

Prints, tab-separated: SMEMs in each output, SMEMs only in each (by read, interval, count and listed hits),
how many of the other tool's own SMEMs list a hit whose reference span holds a non-ACGT base (tools that fill
N with random bases, bwa and minibwa, match through it, which sporos never does), and how many of the rest list
no positions (more occurrences than the tools print), so that the check cannot tell.
"""
import sys


def smems(path):
    out, read = set(), None
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] == "SQ":
            read = f[1]
        elif f[0] == "EM":
            out.add((read, int(f[1]), int(f[2]), int(f[3]), tuple(sorted(f[4:]))))
    return out


seq = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">"))
a, b = smems(sys.argv[2]), smems(sys.argv[3])
only_a, only_b = a - b, b - a


def through_n(record):
    _, start, end, _, hits = record
    for h in hits:
        if h in ("*", "."):
            continue
        pos = int(h.rsplit(":", 1)[1][1:]) - 1
        if any(c not in "ACGTacgt" for c in seq[pos:pos + end - start]):
            return True
    return False


n = [r for r in only_b if through_n(r)]
unlisted = [r for r in only_b if not through_n(r) and all(h in ("*", ".") for h in r[4])]
print(len(a), len(b), len(only_a), len(only_b), len(n), len(unlisted), sep="\t")
