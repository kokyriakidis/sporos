"""Simulated reads for the SMEM benchmark: substitutions only, half reverse-complemented.

    python3 -I sim_reads.py <ref.fa> <reads> <length> <substitution rate> <seed> <out.fa>

Reads start uniformly at random; a window with more than a tenth Ns is drawn again. The benchmark's two sets
(setup_xaros.sh) are 1 M reads of 150 bp and 4,000 of 10 kbp from hg38 chr21 at a 0.2% substitution rate.
"""
import random, sys
ref_path, n, length, err, seed, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
seq = "".join(l.strip() for l in open(ref_path) if not l.startswith(">")).upper()
rnd = random.Random(seed)
comp = str.maketrans("ACGTN", "TGCAN")
with open(out, "w") as f:
    i = 0
    while i < n:
        p = rnd.randrange(0, len(seq) - length)
        r = seq[p:p + length]
        if r.count("N") > length // 10:
            continue
        r = list(r)
        for j in range(length):
            if rnd.random() < err:
                r[j] = rnd.choice("ACGT")
        r = "".join(r)
        if rnd.random() < 0.5:
            r = r.translate(comp)[::-1]
        f.write(f">r{i}\n{r}\n")
        i += 1
