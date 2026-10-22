"""The paper's tables, figure data and every number its text quotes, from the competitor benchmark's results
(benchmarks/competitors/results/: run_seeds.sh, run_lookups.sh, run_mems.sh, run_capped.sh, run_smems.sh).
main_text.py and supplement.py import TABLES, STABLES, V and SPEEDUPS from here, so no number in the paper
is typed by hand.

    python3 tools/paper/tables.py    # prints V, the numbers the text quotes
"""
import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
RES = os.environ.get("SPOROS_RESULTS", os.path.join(ROOT, "benchmarks", "competitors", "results"))


def tsv(name):
    path = os.path.join(RES, name)
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ------------------------------------------------------------------------------- formatting
def num(x, digits=3):
    """Three significant figures, thousands separated; a dash for a missing value."""
    if x is None:
        return "–"
    if x >= 10 ** digits:
        return f"{x:,.0f}"
    return f"{x:#.{digits}g}".rstrip(".")


def secs(x):
    """A wall time as /usr/bin/time measures it (to 0.01 s): two decimals below 10 s, three figures above."""
    if x is None:
        return "–"
    return f"{x:.2f}" if x < 10 else num(x)


def ratio(x):
    return f"{x:.1f}" if x < 10 else f"{x:.0f}"


def bold_min(cells, values):
    """Bold the lowest of `values` (None skipped)."""
    present = [v for v in values if v is not None]
    best = min(present) if present else None
    return [f"**{c}**" if v is not None and v == best else c for c, v in zip(cells, values)]


def bold_max(cells, values):
    present = [v for v in values if v is not None]
    best = max(present) if present else None
    return [f"**{c}**" if v is not None and v == best else c for c, v in zip(cells, values)]


# ------------------------------------------------------------------------------- seeds
def seed_rows():
    """Every kernel timing in seeds_raw.txt: tool, input, mode, k, w, threads, ns/bp. (Parsed from the raw lines,
    which name every variant; run_seeds.sh's own parser once left the ASCII and open-set variants unnamed.)"""
    out, tool, inp = [], None, None
    path = os.path.join(RES, "seeds_raw.txt")
    if not os.path.exists(path):
        return out
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("## "):
            head = line[3:].split(":", 1)[0]
            tool, inp = head.rsplit(" ", 1)
            continue
        if "ns/bp" not in line or "k=" not in line or "out=" not in line:
            continue
        flat = re.sub(r"=\s+", "=", line)
        fields = flat.split()
        mode = next((f for f in fields if re.fullmatch(r"(fwd|can|closed-fwd|closed-can|open-fwd|open-can|fwd-ascii|"
                                                       r"can-ascii|open-set-fwd|open-set-can)", f)), "can" if fields[0] == "mm_sketch" else "")
        kv = dict(f.split("=", 1) for f in fields if "=" in f)
        ns = float(kv["ns/bp"]) if "ns/bp" in kv else float(fields[fields.index("ns/bp") - 1])
        out.append(dict(tool=tool, input=inp, mode=mode, k=int(kv["k"]), w=int(kv["w"]), T=int(kv.get("T", 1)), ns=ns))
    return out


SEEDS = seed_rows()


def seed(tool, inp, mode, k, w, T):
    for r in SEEDS:
        if r["tool"] == tool and r["input"] == inp and r["mode"] == mode and r["k"] == k and r["w"] == w and r["T"] == T:
            return r["ns"]
    return None


SP16, SP8, SMIN, MM2 = "sporos (lanes=16)", "sporos (lanes=8)", "simd-minimizers", "minimap2 mm_sketch"
SEED_SCHEMES = [("fwd", 21, 11, "Minimizers, forward"), ("can", 21, 11, "Minimizers, canonical"),
                ("closed-fwd", 21, 11, "Closed syncmers, forward"), ("closed-can", 21, 11, "Closed syncmers, canonical"),
                ("fwd", 31, 19, "Minimizers, forward"), ("can", 31, 19, "Minimizers, canonical"),
                ("closed-fwd", 31, 19, "Closed syncmers, forward"), ("closed-can", 31, 19, "Closed syncmers, canonical")]


def mm2_seed(inp, mode, k, w, T):
    """mm_sketch computes canonical minimizers only, with k at most 28 (k = 31 runs as 28)."""
    if mode != "can":
        return None
    return seed(MM2, inp, "can", 28 if k == 31 else k, w, T)


def seeds_table(inp):
    rows = []
    for mode, k, w, label in SEED_SCHEMES:
        one = [seed(SP16, inp, mode, k, w, 1), seed(SP8, inp, mode, k, w, 1), seed(SMIN, inp, mode, k, w, 1),
               mm2_seed(inp, mode, k, w, 1)]
        ten = [seed(SP16, inp, mode, k, w, 10), seed(SMIN, inp, mode, k, w, 10), mm2_seed(inp, mode, k, w, 10)]
        if inp != "chm13_chr1":
            ten = []
        cells = [label, f"{k}, {w}"]
        cells += bold_min([num(v) if v else "–" for v in one], one)
        if ten:
            cells += bold_min([num(v) if v else "–" for v in ten], ten)
        rows.append(cells)
    return rows


def seed_speedup(mode, k, w, T, inp="chm13_chr1"):
    s, o = seed(SP16, inp, mode, k, w, T), seed(SMIN, inp, mode, k, w, T)
    return o / s if s and o else None


# ------------------------------------------------------------------------------- lookups
LOOKUPS = {(r["tool"], r["measure"]): (float(r["value"]), r["hits_per_read"]) for r in tsv("lookups.tsv")}


def lookup(tool, measure):
    v = LOOKUPS.get((tool, measure))
    return v[0] if v else None


def hits(tool, measure):
    v = LOOKUPS.get((tool, measure))
    return float(v[1]) if v and v[1] else None


def lookups_table():
    rows = []
    for k, w in ((19, 12), (21, 11)):
        sp, mm = f"sporos k={k} w={w}", f"minimap2 k={k} w={w}"
        rows.append([f"sporos, batch (k = {k}, w = {w})", f"**{num(lookup(sp, 'lookup minimizers, batch'))}**",
                     num(hits(sp, "lookup minimizers, batch")),
                     f"{num(lookup(sp, 'build, 1 threads'))} / {num(lookup(sp, 'build, 8 threads'))}"])
        rows.append([f"sporos, one read per call (k = {k}, w = {w})", num(lookup(sp, "lookup minimizers, per read")),
                     num(hits(sp, "lookup minimizers, per read")), ""])
        rows.append([f"minimap2 (k = {k}, w = {w})", num(lookup(mm, "lookup minimizers, per read")),
                     num(hits(mm, "lookup minimizers, per read")),
                     f"{num(lookup(mm, 'build, 1 threads'))} / {num(lookup(mm, 'build, 8 threads'))}"])
    sa = "strobealign defaults"
    rows.append(["strobealign, one hit per seed", num(lookup(sa, "lookup find_hits (one per seed)")),
                 num(hits(sa, "lookup find_hits (one per seed)")),
                 f"{num(lookup(sa, 'build, 1 threads'))} / {num(lookup(sa, 'build, 8 threads'))}"])
    rows.append(["strobealign, every position", num(lookup(sa, "lookup find_hits + every position")),
                 num(hits(sa, "lookup find_hits + every position")), ""])
    return rows


# ------------------------------------------------------------------------------- MEMs and MUMs
MEMS = tsv("mems.tsv")
PAIRS = [("ecoli", "E. coli (4.6 / 5.5 Mbp)"), ("chr21", "Human / chimpanzee chr21 (47 / 33 Mbp)"),
         ("chr1", "hg38 / CHM13 chr1 (249 / 248 Mbp)")]
PAIR_SHORT = {"ecoli": "E. coli", "chr21": "chr21", "chr1": "chr1"}


def timed(rows, **match):
    """The wall time and peak memory of the one row matching every field, or (None, None) (a timeout reads as
    no time)."""
    for r in rows:
        if all(r.get(k) == str(v) for k, v in match.items()):
            try:
                return float(r["wall_s"]), float(r["max_rss_mb"])
            except ValueError:
                return None, None
    return None, None


def wall(rows, **match):
    return timed(rows, **match)[0]


def rss(rows, **match):
    return timed(rows, **match)[1]


def copmem2(pair, L, T):
    """copMEM2's faster setting at this thread count (default or -mf)."""
    times = [wall(MEMS, pair=pair, L=L, tool=t, mode="maxmatch", threads=T) for t in ("copmem2", "copmem2 -mf")]
    times = [t for t in times if t is not None]
    return min(times) if times else None


def mems_table(L):
    rows = []
    for pair, label in PAIRS:
        s1, s10 = wall(MEMS, pair=pair, L=L, tool="sporos", mode="maxmatch", threads=1), wall(MEMS, pair=pair, L=L, tool="sporos", mode="maxmatch", threads=10)
        c1, c10, mu = copmem2(pair, L, 1), copmem2(pair, L, 10), wall(MEMS, pair=pair, L=L, tool="mummer", mode="maxmatch", threads=1)
        m1, m10 = wall(MEMS, pair=pair, L=L, tool="sporos", mode="mum", threads=1), wall(MEMS, pair=pair, L=L, tool="sporos", mode="mum", threads=10)
        mm = wall(MEMS, pair=pair, L=L, tool="mummer", mode="mum", threads=1)
        one = bold_min([secs(v) if v else "–" for v in (s1, c1, mu)], [s1, c1, mu])
        ten = bold_min([secs(v) if v else "–" for v in (s10, c10)], [s10, c10])
        mum1 = bold_min([secs(v) if v else "–" for v in (m1, mm)], [m1, mm])
        rows.append([label] + [one[0], ten[0], one[1], ten[1], one[2], mum1[0], secs(m10) if m10 else "–", mum1[1]])
    return rows


def mems_memory_table():
    rows = []
    for pair, label in PAIRS:
        for L in ("50", "100"):
            vals = [rss(MEMS, pair=pair, L=L, tool="sporos", mode="maxmatch", threads=1),
                    rss(MEMS, pair=pair, L=L, tool="copmem2", mode="maxmatch", threads=1),
                    rss(MEMS, pair=pair, L=L, tool="copmem2 -mf", mode="maxmatch", threads=1),
                    rss(MEMS, pair=pair, L=L, tool="mummer", mode="maxmatch", threads=1),
                    rss(MEMS, pair=pair, L=L, tool="sporos", mode="mum", threads=1),
                    rss(MEMS, pair=pair, L=L, tool="mummer", mode="mum", threads=1)]
            rows.append([f"{PAIR_SHORT[pair]}, L = {L}"] + [f"{v:,.0f}" if v else "–" for v in vals])
    return rows


def mumreference_table():
    rows = []
    for pair, label in PAIRS:
        for L in ("50", "100"):
            s1 = wall(MEMS, pair=pair, L=L, tool="sporos", mode="mumreference", threads=1)
            s10 = wall(MEMS, pair=pair, L=L, tool="sporos", mode="mumreference", threads=10)
            mu = wall(MEMS, pair=pair, L=L, tool="mummer", mode="mumreference", threads=1)
            rows.append([f"{PAIR_SHORT[pair]}, L = {L}", secs(s1) if s1 else "–", secs(s10) if s10 else "–",
                         secs(mu) if mu else "–", ratio(mu / s1) if mu and s1 else "–"])
    return rows


def mumemto_table():
    rows = []
    for pair, label in PAIRS:
        for L in ("50", "100"):
            r = next((r for r in MEMS if r["pair"] == pair and r["L"] == L and r["tool"] == "mumemto"), None)
            if not r:
                continue
            t = r["wall_s"]
            rows.append([f"{PAIR_SHORT[pair]}, L = {L}", t if not t[0].isdigit() else num(float(t)),
                         num(wall(MEMS, pair=pair, L=L, tool="sporos", mode="mum", threads=1) or 0)])
    return rows


CHECKS = tsv("mems_checks.tsv")


def check_counts():
    exact = [r for r in CHECKS if r["result"] in ("PASS", "FAIL")]
    return sum(r["result"] == "PASS" for r in exact), len(exact)


# ------------------------------------------------------------------------------- capped MEMs
CAPPED = tsv("capped.tsv")
CAPPED_CHECKS = tsv("capped_checks.tsv")


def capped_table():
    rows = []
    for pair, label in PAIRS:
        for kind, name in (("per-genome", "per genome"), ("total", "in total")):
            for N in sorted({r["mode"].split()[-1] for r in CAPPED if r["pair"] == pair}, key=int):
                mode = f"{kind} cap {N}"
                s1, s10 = wall(CAPPED, pair=pair, tool="sporos", mode=mode, threads=1), wall(CAPPED, pair=pair, tool="sporos", mode=mode, threads=10)
                r1, mem = timed(CAPPED, pair=pair, tool="mumemto", mode=mode, threads=1)
                if s1 is None:
                    continue
                chk = next((c for c in CAPPED_CHECKS if c["pair"] == pair and c["mode"] == mode), None)
                share = (f"{100 * int(chk['common']) / int(chk['mumemto_fwd_pairs']):.0f} / "
                         f"{100 * int(chk['common']) / int(chk['sporos_fwd_pairs']):.0f}"
                         if chk and chk["common"] not in ("NA", "") and int(chk["mumemto_fwd_pairs"]) else "–")
                rows.append([f"{PAIR_SHORT[pair]}, at most {N} {name}", f"**{secs(s1)}**", secs(s10) if s10 else "–",
                             secs(r1) if r1 else "DNF", ratio(r1 / s1) if r1 else "–",
                             f"{rss(CAPPED, pair=pair, tool='sporos', mode=mode, threads=1):,.0f}",
                             f"{mem:,.0f}" if mem else "–", share])
    return rows


# ------------------------------------------------------------------------------- SMEMs
SMEMS = tsv("smems.tsv")
SMEM_CHECKS = tsv("smems_checks.tsv")
READSETS = [("150bp", "1 M reads of 150 bp"), ("10kb", "4,000 reads of 10 kbp")]


def smem_mode(f):
    return "SMEM" if f == 1 else f"f-SMEM f={f}"


def smem_times(reads, f):
    m = smem_mode(f)
    return dict(
        s1=wall(SMEMS, reads=reads, tool="sporos", mode=m, threads=1),
        s10=wall(SMEMS, reads=reads, tool="sporos", mode=m, threads=10),
        mb=wall(SMEMS, reads=reads, tool="minibwa", mode=m, threads=1),
        bwa=wall(SMEMS, reads=reads, tool="bwa fastmap", mode=m, threads=1) if f == 1 else None,
        rb1=wall(SMEMS, reads=reads, tool="ropebwt3 mem", mode=f"{m} (no positions)", threads=1),
        rb10=wall(SMEMS, reads=reads, tool="ropebwt3 mem", mode=f"{m} (no positions)", threads=10))


def smems_table():
    rows = []
    for reads, label in READSETS:
        for f in (1, 2, 5):
            t = smem_times(reads, f)
            if t["s1"] is None:
                continue
            one = bold_min([secs(v) if v else "–" for v in (t["s1"], t["mb"], t["bwa"], t["rb1"])], [t["s1"], t["mb"], t["bwa"], t["rb1"]])
            ten = bold_min([secs(v) if v else "–" for v in (t["s10"], t["rb10"])], [t["s10"], t["rb10"]])
            rows.append([f"{label}, {'SMEMs' if f == 1 else f'f = {f}'}", one[0], ten[0], one[1], one[2], one[3], ten[1]])
    return rows


def smem_builds():
    out = {}
    for tool in ("sporos index", "minibwa index", "ropebwt3 build", "bwa index"):
        out[tool] = (timed(SMEMS, reads="chr21", tool=tool, mode="build", threads=1),
                     timed(SMEMS, reads="chr21", tool=tool, mode="build", threads=10))
    return out


def smem_check(reads, f, name):
    return next((c for c in SMEM_CHECKS if c["reads"] == reads and c["f"] == str(f) and c["check"] == name), None)


# ------------------------------------------------------------------------------- speedups (figure 2)
def _fastest(pairs):
    pairs = [(n, v) for n, v in pairs if v]
    return min(pairs, key=lambda p: p[1]) if pairs else (None, None)


def speedups():
    out = []

    def add(label, rival, t1, t10):
        if t1 or t10:
            out.append(dict(label=label, rival=rival, t1=t1, t10=t10))

    add("Minimizers, canonical (chr1)", "simd-minimizers", seed_speedup("can", 21, 11, 1), seed_speedup("can", 21, 11, 10))
    add("Closed syncmers, canonical (chr1)", "simd-minimizers", seed_speedup("closed-can", 21, 11, 1), seed_speedup("closed-can", 21, 11, 10))
    sp, mm = lookup("sporos k=19 w=12", "lookup minimizers, batch"), lookup("minimap2 k=19 w=12", "lookup minimizers, per read")
    add("Seed-index lookups (150 bp)", "minimap2", sp / mm if sp and mm else None, None)
    for mode, label in (("maxmatch", "MEMs"), ("mum", "MUMs")):
        s1 = wall(MEMS, pair="chr1", L="50", tool="sporos", mode=mode, threads=1)
        s10 = wall(MEMS, pair="chr1", L="50", tool="sporos", mode=mode, threads=10)
        if mode == "maxmatch":
            name1, o1 = _fastest([("copMEM2", copmem2("chr1", "50", 1)), ("MUMmer4", wall(MEMS, pair="chr1", L="50", tool="mummer", mode=mode, threads=1))])
            name10, o10 = _fastest([("copMEM2", copmem2("chr1", "50", 10)), ("MUMmer4", wall(MEMS, pair="chr1", L="50", tool="mummer", mode=mode, threads=1))])
        else:
            name1 = name10 = "MUMmer4"
            o1 = o10 = wall(MEMS, pair="chr1", L="50", tool="mummer", mode=mode, threads=1)
        add(f"{label}, L = 50 (chr1)", name1, o1 / s1 if o1 and s1 else None, o10 / s10 if o10 and s10 else None)
    s1 = wall(CAPPED, pair="chr21", tool="sporos", mode="total cap 10", threads=1)
    s10 = wall(CAPPED, pair="chr21", tool="sporos", mode="total cap 10", threads=10)
    mu = wall(CAPPED, pair="chr21", tool="mumemto", mode="total cap 10", threads=1)
    add("Capped MEMs, ≤ 10 copies (chr21)", "Mumemto", mu / s1 if mu and s1 else None, mu / s10 if mu and s10 else None)
    for reads, f, label in (("150bp", 1, "SMEMs, 150 bp reads"), ("10kb", 1, "SMEMs, 10 kbp reads"),
                            ("150bp", 2, "f-SMEMs, f = 2 (150 bp)")):
        t = smem_times(reads, f)
        n1, o1 = _fastest([("minibwa", t["mb"]), ("bwa", t["bwa"]), ("ropebwt3", t["rb1"])])
        n10, o10 = _fastest([("minibwa", t["mb"]), ("bwa", t["bwa"]), ("ropebwt3", t["rb10"])])
        add(label, n1 if n1 == n10 else f"{n1}; {n10}", o1 / t["s1"] if o1 and t["s1"] else None,
            o10 / t["s10"] if o10 and t["s10"] else None)
    return out


SPEEDUPS = speedups()

# ------------------------------------------------------------------------------- versions and machine
VERSIONS = {r["tool"]: r["version"] for r in tsv("versions.tsv")}


def machine():
    path = os.path.join(RES, "machine.txt")
    if not os.path.exists(path):
        return {}
    out = {}
    for line in open(path):
        if ": " in line:
            k, v = line.rstrip("\n").split(": ", 1)
            out.setdefault(k, v)
    return out


MACHINE = machine()


# ------------------------------------------------------------------------------- the numbers the text quotes
def _range(values, fmt=ratio):
    """The least and greatest of `values` as the text quotes them: "a–b", or one value when they agree."""
    values = [v for v in values if v]
    return (fmt(min(values)), fmt(max(values))) if values else ("?", "?")


def span(lo, hi):
    return lo if lo == hi else f"{lo}–{hi}"


def quoted():
    V = {}
    # Seeds: sporos (16 lanes) against simd-minimizers, every scheme, both inputs, one thread; chr1, ten threads.
    one = [seed_speedup(m, k, w, 1, inp) for m, k, w, _ in SEED_SCHEMES for inp in ("random100M", "chm13_chr1")]
    V["seed_t1_min"], V["seed_t1_max"] = _range(one)
    ten = [seed_speedup(m, k, w, 10) for m, k, w, _ in SEED_SCHEMES]
    V["seed_t10_min"], V["seed_t10_max"] = _range(ten)
    eight = [seed(SMIN, "chm13_chr1", m, k, w, 1) / seed(SP8, "chm13_chr1", m, k, w, 1) for m, k, w, _ in SEED_SCHEMES
             if seed(SMIN, "chm13_chr1", m, k, w, 1) and seed(SP8, "chm13_chr1", m, k, w, 1)]
    V["seed8_min"], V["seed8_max"] = _range(eight)
    wide = [seed(SP8, inp, m, k, w, 1) / seed(SP16, inp, m, k, w, 1) for m, k, w, _ in SEED_SCHEMES
            for inp in ("random100M", "chm13_chr1") if seed(SP8, inp, m, k, w, 1) and seed(SP16, inp, m, k, w, 1)]
    V["seed16_vs_8_min"], V["seed16_vs_8_max"] = _range(wide, lambda x: f"{x:.2f}")
    mm = [mm2_seed("chm13_chr1", m, k, w, 1) / seed(SP16, "chm13_chr1", m, k, w, 1) for m, k, w, _ in SEED_SCHEMES
          if mm2_seed("chm13_chr1", m, k, w, 1) and seed(SP16, "chm13_chr1", m, k, w, 1)]
    V["seed_mm2_min"], V["seed_mm2_max"] = _range(mm)
    v = seed(SP16, "chm13_chr1", "can", 21, 11, 1)
    V["seed_chr1_ns"] = num(v) if v else "?"
    V["seed_chr1_gbps"] = f"{1 / v:.2f}" if v else "?"
    # Lookups.
    sp, m2 = lookup("sporos k=19 w=12", "lookup minimizers, batch"), lookup("minimap2 k=19 w=12", "lookup minimizers, per read")
    sa = lookup("strobealign defaults", "lookup find_hits (one per seed)")
    one_read = lookup("sporos k=19 w=12", "lookup minimizers, per read")
    V["lookup_sporos"], V["lookup_mm2"], V["lookup_sa"] = (num(x) if x else "?" for x in (sp, m2, sa))
    V["lookup_vs_mm2"] = ratio(sp / m2) if sp and m2 else "?"
    V["lookup_vs_sa"] = ratio(sp / sa) if sp and sa else "?"
    V["lookup_one_vs_mm2"] = ratio(one_read / m2) if one_read and m2 else "?"
    # MEMs.
    cm = [copmem2(p, L, 1) / wall(MEMS, pair=p, L=L, tool="sporos", mode="maxmatch", threads=1)
          for p, _ in PAIRS for L in ("50", "100") if copmem2(p, L, 1) and wall(MEMS, pair=p, L=L, tool="sporos", mode="maxmatch", threads=1)]
    V["copmem2_t1_min"], V["copmem2_t1_max"] = _range(cm)
    cm1 = [copmem2("chr1", L, 1) / wall(MEMS, pair="chr1", L=L, tool="sporos", mode="maxmatch", threads=1) for L in ("50", "100")
           if copmem2("chr1", L, 1) and wall(MEMS, pair="chr1", L=L, tool="sporos", mode="maxmatch", threads=1)]
    V["copmem2_chr1_min"], V["copmem2_chr1_max"] = _range(cm1)
    cm10 = [copmem2("chr1", L, 10) / wall(MEMS, pair="chr1", L=L, tool="sporos", mode="maxmatch", threads=10) for L in ("50", "100")
            if copmem2("chr1", L, 10) and wall(MEMS, pair="chr1", L=L, tool="sporos", mode="maxmatch", threads=10)]
    V["copmem2_chr1_t10_min"], V["copmem2_chr1_t10_max"] = _range(cm10)
    mu = [wall(MEMS, pair=p, L=L, tool="mummer", mode=m, threads=1) / wall(MEMS, pair=p, L=L, tool="sporos", mode=m, threads=1)
          for p, _ in PAIRS for L in ("50", "100") for m in ("maxmatch", "mum", "mumreference")
          if wall(MEMS, pair=p, L=L, tool="mummer", mode=m, threads=1) and wall(MEMS, pair=p, L=L, tool="sporos", mode=m, threads=1)]
    V["mummer_min"], V["mummer_max"] = _range(mu)
    V["checks_pass"], V["checks_total"] = check_counts()
    big = max((int(r["n_a"]) for r in CHECKS if r["result"] == "PASS" and r["n_a"].isdigit()), default=0)
    V["checks_max_mems"] = f"{big / 1e6:.1f}"
    s = wall(MEMS, pair="chr1", L="50", tool="sporos", mode="maxmatch", threads=1)
    V["chr1_mems_t1"] = num(s) if s else "?"
    s = wall(MEMS, pair="chr1", L="50", tool="sporos", mode="maxmatch", threads=10)
    V["chr1_mems_t10"] = num(s) if s else "?"
    V["chr1_rss_gb"] = f"{rss(MEMS, pair='chr1', L='50', tool='sporos', mode='maxmatch', threads=1) / 1024:.1f}" if rss(MEMS, pair="chr1", L="50", tool="sporos", mode="maxmatch", threads=1) else "?"
    # Capped MEMs.
    cap = [wall(CAPPED, pair=p, tool="mumemto", mode=r["mode"], threads=1) / wall(CAPPED, pair=p, tool="sporos", mode=r["mode"], threads=1)
           for p, _ in PAIRS for r in CAPPED if r["pair"] == p and r["tool"] == "sporos" and r["threads"] == "1"
           and wall(CAPPED, pair=p, tool="mumemto", mode=r["mode"], threads=1)]
    V["capped_min"], V["capped_max"] = _range(cap, lambda x: f"{x:.0f}")
    for pair in ("ecoli", "chr21"):
        rows = [c for c in CAPPED_CHECKS if c["pair"] == pair and c["common"] not in ("NA", "")]
        V[f"share_mm_{pair}_min"], V[f"share_mm_{pair}_max"] = _range(
            [100 * int(c["common"]) / int(c["mumemto_fwd_pairs"]) for c in rows], lambda x: f"{x:.0f}")
        V[f"share_sp_{pair}_min"], V[f"share_sp_{pair}_max"] = _range(
            [100 * int(c["common"]) / int(c["sporos_fwd_pairs"]) for c in rows], lambda x: f"{x:.0f}")
    # SMEMs.
    for reads in ("150bp", "10kb"):
        t = smem_times(reads, 1)
        V[f"smem_{reads}_mb"] = ratio(t["mb"] / t["s1"]) if t["mb"] and t["s1"] else "?"
        V[f"smem_{reads}_bwa"] = ratio(t["bwa"] / t["s1"]) if t["bwa"] and t["s1"] else "?"
        V[f"smem_{reads}_rb"] = ratio(t["rb1"] / t["s1"]) if t["rb1"] and t["s1"] else "?"
        V[f"smem_{reads}_s1"] = num(t["s1"]) if t["s1"] else "?"
        V[f"smem_{reads}_s10"] = num(t["s10"]) if t["s10"] else "?"
    for f in (2, 5):
        t = smem_times("150bp", f)
        V[f"fsmem{f}_mb"] = ratio(t["mb"] / t["s1"]) if t["mb"] and t["s1"] else "?"
        V[f"fsmem{f}_rb"] = ratio(t["rb1"] / t["s1"]) if t["rb1"] and t["s1"] else "?"
    plain, two = smem_times("150bp", 1)["s1"], smem_times("150bp", 2)["s1"]
    V["fsmem_cost"] = ratio(two / plain) if two and plain else "?"
    c = smem_check("150bp", 1, "sporos vs bwa fastmap")
    V["bwa_only_bwa"], V["bwa_only_bwa_n"] = (c["only_other"], c["other_through_N"]) if c else ("?", "?")
    V["bwa_only_sporos"] = c["only_sporos"] if c else "?"
    V["bwa_smems"] = f"{int(c['n_other']):,}" if c else "?"
    cpu = MACHINE.get("cpu", "")
    model = re.search(r"Core\(TM\) (\S+) CPU @ ([\d.]+)GHz", cpu)
    V["machine_cpu"] = f"Intel Core {model.group(1)} at {float(model.group(2)):.1f} GHz" if model else cpu
    return V


V = quoted()

# ------------------------------------------------------------------------------- tables
TABLES = {
    "seeds": {
        "caption": "Minimizers and closed syncmers of CHM13 chr1 (248 Mbp): ns per base, the best of 5 to 7 runs on one "
                   "thread and of 3 on ten. Every tool reads its input already encoded (simd-minimizers packed at two bits per "
                   "base, sporos at one byte per base, minimap2's mm_sketch in ASCII) and writes positions to memory. sporos "
                   "runs 16 lanes of 32 bits (AVX-512) or 8 (the width simd-minimizers uses, AVX2). mm_sketch computes canonical "
                   "minimizers only, of k at most 28 (the k = 31 rows run k = 28). Bold: fastest at each thread count.",
        "header": ["Scheme", "k, w", "sporos-16", "sporos-8", "simd-min.", "mm_sketch", "sporos-16", "simd-min.", "mm_sketch"],
        "groups": [("One thread", 2, 5), ("Ten threads", 6, 8)],
        "align": "llrrrrrrr",
        "fit": True,
        "rows": seeds_table("chm13_chr1"),
    },
    "lookups": {
        "caption": "Seed-index lookups of 1 M simulated 150 bp reads against hg38 chr21 on one thread: throughput in Mbp of "
                   "reads per second, the best of three runs, reads in memory, timing only the loop that seeds each read and "
                   "copies out its hits; index build time on 1 and 8 threads. sporos and minimap2 index the same (w, k)-minimizers "
                   "under different hashes and skip seeds with more than 20 hits; strobealign indexes randstrobes with its defaults "
                   "for 150 bp reads and does more per seed (both orientations, partial-strobe lookups), so its rows are a reference "
                   "point. minimap2's and strobealign's build times include reading the FASTA file.",
        "header": ["Seeding", "Mbp/s", "Hits per read", "Build (s), 1 / 8 threads"],
        "align": "lrrr",
        "fit": True,
        "rows": lookups_table(),
    },
    "mems": {
        "caption": "MEMs and MUMs of at least 50 bp between whole chromosomes, forward strand: wall time in seconds, the best of "
                   "three runs (one run past 60 s), output written to a file. copMEM2 is the faster of its default and `-mf` "
                   "settings; it and MUMmer4 gain nothing from more threads on a single query sequence (Supplementary Section S1). "
                   "Every sporos output equals MUMmer4's, and its MEMs equal copMEM2's. Bold: fastest at each thread count.",
        "header": ["Pair", "sporos, 1 thr.", "sporos, 10 thr.", "copMEM2, 1 thr.", "copMEM2, 10 thr.", "MUMmer4",
                   "sporos, 1 thr.", "sporos, 10 thr.", "MUMmer4"],
        "groups": [("MEMs (`-maxmatch`)", 1, 5), ("MUMs (`-mum`)", 6, 8)],
        "align": "lrrrrrrrr",
        "fit": True,
        "rows": mems_table("50"),
    },
    "capped": {
        "caption": "Frequency-capped MEMs of at least 50 bp on both strands: wall time in seconds and peak memory in MB. "
                   "Mumemto (one thread) finds the multi-MEMs of the two genomes as one collection, capped per genome (`-k 2 -f N`) "
                   "or in total (`-k 2 -f 0 -F N`); sporos finds the MEMs of the query against the reference whose string occurs at "
                   "most N times in the reference (`-maxmatch --max-occ N`, for the total cap) or in the reference and the query "
                   "(`-mum --max-occ N`, per genome). The definitions differ (Section 4.4): the last column gives the share of "
                   "Mumemto's forward-strand position pairs that sporos also reports, and of sporos's that Mumemto reports.",
        "header": ["Pair and cap", "sporos, 1 thr.", "sporos, 10 thr.", "Mumemto", "Speedup", "sporos MB", "Mumemto MB", "Shared (%)"],
        "align": "lrrrrrrr",
        "fit": True,
        "rows": capped_table(),
    },
    "smems": {
        "caption": "SMEMs (f = 1) and f-SMEMs of at least 19 bp of simulated reads against hg38 chr21, both strands: wall time in "
                   "seconds, the best of three runs, the index loaded from disk. sporos, minibwa and bwa fastmap list the positions of "
                   "SMEMs with at most 20 hits; ropebwt3 is run without positions, its default, which favours it. minibwa and bwa "
                   "fastmap are single-threaded; bwa computes no f-SMEMs. Bold: fastest at each thread count.",
        "header": ["Reads", "sporos", "minibwa", "bwa", "ropebwt3", "sporos", "ropebwt3"],
        "groups": [("One thread", 1, 4), ("Ten threads", 5, 6)],
        "align": "lrrrrrr",
        "fit": True,
        "rows": smems_table(),
    },
}

STABLES = {
    "seeds_random": {
        "label": "S2",
        "caption": "Minimizers and closed syncmers of 100 Mbp of uniformly random DNA on one thread: ns per base, as in main Table 1. "
                   "Each harness draws its own random sequence.",
        "header": ["Scheme", "k, w", "sporos-16", "sporos-8", "simd-min.", "mm_sketch"],
        "align": "llrrrr",
        "rows": seeds_table("random100M"),
    },
    "mems100": {
        "label": "S3",
        "caption": "MEMs and MUMs of at least 100 bp, as main Table 3.",
        "header": ["Pair", "sporos, 1 thr.", "sporos, 10 thr.", "copMEM2, 1 thr.", "copMEM2, 10 thr.", "MUMmer4",
                   "sporos, 1 thr.", "sporos, 10 thr.", "MUMmer4"],
        "groups": [("MEMs (`-maxmatch`)", 1, 5), ("MUMs (`-mum`)", 6, 8)],
        "align": "lrrrrrrrr",
        "fit": True,
        "rows": mems_table("100"),
    },
    "mumreference": {
        "label": "S4",
        "caption": "MEMs unique in the reference (`-mumreference`): wall time in seconds and the speedup over MUMmer4 on one thread.",
        "header": ["Pair", "sporos, 1 thr.", "sporos, 10 thr.", "MUMmer4", "Speedup"],
        "align": "lrrrr",
        "rows": mumreference_table(),
    },
    "memory": {
        "label": "S5",
        "caption": "Peak resident memory in MB of the one-thread runs of main Table 3 and Supplementary Table S3.",
        "header": ["Pair", "sporos MEMs", "copMEM2", "copMEM2 -mf", "MUMmer4 MEMs", "sporos MUMs", "MUMmer4 MUMs"],
        "align": "lrrrrrr",
        "rows": mems_memory_table(),
    },
}


def versions_table():
    names = {"sporos": "sporos", "mummer": "MUMmer4", "mumemto": "Mumemto", "bwa": "BWA", "ropebwt3": "ropebwt3",
             "minibwa": "minibwa", "copMEM2": "copMEM2", "minimap2 mm_sketch": "minimap2", "simd-minimizers":
             "simd-minimizers", "strobealign": "strobealign", "gcc": "GCC"}
    return [[names.get(k, k), v] for k, v in VERSIONS.items() if k != "sporos"]


def builds_table():
    rows = []
    labels = {"sporos index": "sporos", "minibwa index": "minibwa", "ropebwt3 build": "ropebwt3", "bwa index": "BWA"}
    for tool, ((t1, m1), (t10, m10)) in smem_builds().items():
        rows.append([labels[tool], secs(t1) if t1 else "–", secs(t10) if t10 else "–",
                     f"{max(m for m in (m1, m10) if m):,.0f}" if (m1 or m10) else "–"])
    return rows


def smem_checks_table():
    rows = []
    for c in SMEM_CHECKS:
        if c["check"] == "sporos -t10=-t1":
            continue
        unlisted = c.get("other_unlisted") or ""
        rows.append([f"{dict(READSETS)[c['reads']]}, {'SMEMs' if c['f'] == '1' else 'f = ' + c['f']}",
                     c["check"].replace("sporos vs ", ""), f"{int(c['n_sporos']):,}", f"{int(c['n_other']):,}",
                     f"{int(c['only_sporos']):,}", f"{int(c['only_other']):,}", f"{int(c['other_through_N']):,}",
                     f"{int(unlisted):,}" if unlisted.isdigit() else "–"])
    return rows


def smem_differences():
    """Over every comparison: the other tools' own SMEMs, those through an N, those listing no positions."""
    rows = [c for c in SMEM_CHECKS if c["check"] != "sporos -t10=-t1"]
    only = sum(int(c["only_other"]) for c in rows)
    through = sum(int(c["other_through_N"]) for c in rows)
    unlisted = sum(int(c["other_unlisted"]) for c in rows if (c.get("other_unlisted") or "").isdigit())
    return only, through, unlisted


def threads_same():
    """Whether every ten-thread SMEM output equalled the one-thread one."""
    same = [c["other_through_N"] for c in SMEM_CHECKS if c["check"] == "sporos -t10=-t1"]
    return bool(same) and all(x == "PASS" for x in same)


STABLES["versions"] = {
    "label": "S1",
    "caption": "Versions of the tools compared. Tools built from source were compiled for the benchmark machine; MUMmer4, "
               "Mumemto, BWA and ropebwt3 are Bioconda packages.",
    "header": ["Tool", "Version"],
    "align": "ll",
    "rows": versions_table(),
}
STABLES["builds"] = {
    "label": "S6",
    "caption": "Building the FM-indexes of hg38 chr21 for Table 5: wall time in seconds on 1 and 10 threads, and the larger "
               "peak memory in MB. BWA's index builds on one thread.",
    "header": ["Tool", "1 thread", "10 threads", "Peak MB"],
    "align": "lrrr",
    "rows": builds_table(),
}
STABLES["smem_checks"] = {
    "label": "S7",
    "caption": "SMEM outputs compared, by read, interval, number of occurrences and listed positions. The last column counts "
               "the other tool's own SMEMs that list a position whose reference span holds a non-ACGT base under Through N "
               "(BWA and minibwa replace each N of the reference by a random base, and match through it), and those of the rest "
               "that list no positions, having more than 20 occurrences, under No positions.",
    "header": ["Reads", "Against", "sporos", "Other", "Only sporos", "Only other", "Through N", "No positions"],
    "align": "llrrrrrr",
    "fit": True,
    "rows": smem_checks_table(),
}
V["smem_threads_same"] = threads_same()
_cm = [rss(MEMS, pair=p, L=L, tool=t, mode="maxmatch", threads=1) for p, _ in PAIRS for L in ("50", "100")
       for t in ("copmem2", "copmem2 -mf")]
V["copmem2_rss_min"] = f"{min(m for m in _cm if m) / 1024:.1f}" if any(_cm) else "?"
_only, _through, _unlisted = smem_differences()
V["smem_only_other"], V["smem_through_n"], V["smem_unlisted"] = f"{_only:,}", f"{_through:,}", f"{_unlisted:,}"
V["smem_unexplained"] = f"{_only - _through - _unlisted:,}"


if __name__ == "__main__":
    for k, v in V.items():
        print(f"{k:24} {v}")
    print()
    for r in SPEEDUPS:
        print(r)
