"""The README's measured numbers, from the same results as the paper (tables.py), so the two cannot drift apart:
the speed bullets under the title, the benchmark chart's alt text and the Benchmarks section's tables. Each is
written between a pair of markers in README.md (`<!-- speed:begin -->` ... `<!-- speed:end -->`, and `chart`,
`results`); the rest of the README stays as written. render.py calls update(); scripts/benchmark_chart.py draws
the chart from panels().
"""
import os
import re

import tables as T

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
README = os.path.join(ROOT, "README.md")
V = T.V


def s(x):
    """Seconds as the README prints them."""
    return "–" if x is None else (f"{x:.2f} s" if x < 10 else f"{x:.1f} s" if x < 100 else f"{x:.0f} s")


def mb(x):
    return "–" if x is None else (f"{x / 1024:.1f} GB" if x >= 1024 else f"{x:.0f} MB")


def x(v):
    return f"{v:.1f}" if v < 10 else f"{v:.0f}"


# ------------------------------------------------------------------------------- the speed bullets
def ten_thread_seeds():
    ten = [T.seed_speedup(m, k, w, 10) for m, k, w, _ in T.SEED_SCHEMES]
    ten = [v for v in ten if v]
    slow = [v for v in ten if v < 1]
    if slow:
        return f"from {1 / min(slow):.1f}x slower to {x(max(ten))}x faster on ten ({len(slow)} of {len(ten)} schemes slower)"
    return f"{x(min(ten))}–{x(max(ten))}x faster on ten"


def speed():
    builds = T.smem_builds()
    sp1, sp10 = builds["sporos index"][0][0], builds["sporos index"][1][0]
    mb1, mb10 = builds["minibwa index"][0][0], builds["minibwa index"][1][0]
    faster_build = sp1 and mb1 and sp10 and mb10 and sp1 < mb1 and sp10 < mb10
    capped = T.wall(T.CAPPED, pair="chr21", tool="sporos", mode="total cap 10", threads=1)
    mumemto = T.wall(T.CAPPED, pair="chr21", tool="mumemto", mode="total cap 10", threads=1)
    lines = [
        f"  - Minimizers and syncmers run {T.span(V['seed_t1_min'], V['seed_t1_max'])}x faster than simd-minimizers on one "
        f"thread, and {ten_thread_seeds()}.",
        f"  - Seed-index lookups of reads run {V['lookup_vs_mm2']}x faster than minimap2's ({V['lookup_one_vs_mm2']}x one "
        f"read per call) and {V['lookup_vs_sa']}x faster than strobealign's.",
        f"  - MEMs between hg38 and CHM13 chr1 run {T.span(V['copmem2_chr1_min'], V['copmem2_chr1_max'])}x faster than "
        f"copMEM2 on one thread and {T.span(V['copmem2_chr1_t10_min'], V['copmem2_chr1_t10_max'])}x faster with 10 threads.",
        f"  - MEMs and MUMs run {T.span(V['mummer_min'], V['mummer_max'])}x faster than MUMmer4, with identical output.",
        f"  - Frequency-capped MEMs and MUMs (Mumemto's `-f`/`-F`) run {T.span(V['capped_min'], V['capped_max'])}x faster "
        f"than Mumemto ({s(capped)} against {s(mumemto)} between human and chimpanzee chr21); f-SMEMs (at least f "
        f"occurrences) {V['fsmem2_mb']}x faster than minibwa's `-s`.",
        f"  - SMEMs of reads run {V['smem_150bp_mb']}x faster than minibwa and {V['smem_150bp_bwa']}x faster than "
        f"`bwa fastmap` on one thread" + (", from an FM-index that builds faster than minibwa's." if faster_build else "."),
    ]
    return "\n".join(lines)


# ------------------------------------------------------------------------------- the chart
def panels():
    """The README chart's panels: (title, subtitle, [(tool, value)], unit), one benchmark each, lower is better."""
    mm1 = T.wall(T.MEMS, pair="chr1", L="50", tool="mummer", mode="maxmatch", threads=1)
    t150 = T.smem_times("150bp", 1)
    return [
        ("MEMs, hg38 vs CHM13 chr1, L = 50", "1 thread, wall time in seconds",
         [("sporos", T.wall(T.MEMS, pair="chr1", L="50", tool="sporos", mode="maxmatch", threads=1)),
          ("copMEM2", T.copmem2("chr1", "50", 1)), ("MUMmer4", mm1)], " s"),
        ("MEMs with 10 threads, same pair", "wall time in seconds (MUMmer4 runs single-threaded)",
         [("sporos", T.wall(T.MEMS, pair="chr1", L="50", tool="sporos", mode="maxmatch", threads=10)),
          ("copMEM2", T.copmem2("chr1", "50", 10))], " s"),
        ("MUMs (-mum), same pair", "1 thread, wall time in seconds",
         [("sporos", T.wall(T.MEMS, pair="chr1", L="50", tool="sporos", mode="mum", threads=1)),
          ("MUMmer4", T.wall(T.MEMS, pair="chr1", L="50", tool="mummer", mode="mum", threads=1))], " s"),
        ("SMEMs of 1 M 150 bp reads against hg38 chr21, L = 19", "1 thread, wall time in seconds (index built beforehand)",
         [("sporos", t150["s1"]), ("minibwa", t150["mb"]), ("bwa fastmap", t150["bwa"])], " s"),
        ("Canonical minimizers, CHM13 chr1, k = 21, w = 11", "1 thread, nanoseconds per base",
         [("sporos", T.seed(T.SP16, "chm13_chr1", "can", 21, 11, 1)), ("simd-minimizers", T.seed(T.SMIN, "chm13_chr1", "can", 21, 11, 1)),
          ("minimap2", T.mm2_seed("chm13_chr1", "can", 21, 11, 1))], " ns"),
    ]


def chart():
    def v(value, unit):
        return f"{value:.3g}" + (" s" if unit == " s" else " ns per base")

    alt = "Bar charts, lower is better. " + " ".join(
        f"{title}, {sub.split(',')[0]}: " + ", ".join(f"{tool} {v(value, unit)}" for tool, value in bars) + "."
        for title, sub, bars, unit in panels())
    return ('<p align="center">\n  <picture>\n'
            '    <source media="(prefers-color-scheme: dark)" srcset="docs/benchmarks-dark.svg">\n'
            f'    <img alt="{alt}" src="docs/benchmarks-light.svg" width="820">\n  </picture>\n</p>')


# ------------------------------------------------------------------------------- the results tables
def results():
    out = []
    pairs = [("chr1", "100", "hg38 vs CHM13 chr1, L=100"), ("chr1", "50", f"hg38 vs CHM13 chr1, L=50 ({V['checks_max_mems']} M MEMs)"),
             ("chr21", "50", "hg38 vs chimp chr21, L=50")]
    out += ["| MEMs (1 / 10 threads) | sporos | copMEM2, its faster setting | MUMmer4 (1 thread) |", "|---|---|---|---|"]
    for pair, L, label in pairs:
        w = lambda tool, mode, t: T.wall(T.MEMS, pair=pair, L=L, tool=tool, mode=mode, threads=t)
        out.append(f"| {label} | **{s(w('sporos', 'maxmatch', 1))} / {s(w('sporos', 'maxmatch', 10))}** | "
                   f"{s(T.copmem2(pair, L, 1))} / {s(T.copmem2(pair, L, 10))} | {s(w('mummer', 'maxmatch', 1))} |")
    out += ["", "| MUMs, `-mum` | sporos (1 / 10 threads) | MUMmer4 (1 thread) |", "|---|---|---|"]
    for pair, L, label in pairs:
        w = lambda tool, t: T.wall(T.MEMS, pair=pair, L=L, tool=tool, mode="mum", threads=t)
        out.append(f"| {label.split(' (')[0]} | **{s(w('sporos', 1))} / {s(w('sporos', 10))}** | {s(w('mummer', 1))} |")
    builds = T.smem_builds()
    out += ["", "| SMEMs vs hg38 chr21, L=19 | sporos (1 / 10 threads) | minibwa `fastmap -b 64` (1 thread) | "
                "`bwa fastmap` (1 thread) |", "|---|---|---|---|"]
    for reads, label in (("150bp", "1 M simulated 150 bp reads"), ("10kb", "4 K simulated 10 kb reads")):
        t = T.smem_times(reads, 1)
        out.append(f"| {label} | **{s(t['s1'])} / {s(t['s10'])}** | {s(t['mb'])} | {s(t['bwa'])} |")
    (sp1, spm1), (sp10, spm10) = builds["sporos index"]
    (mb1, _), (mb10, mbm10) = builds["minibwa index"]
    (bw1, bwm), _ = builds["bwa index"]
    out.append(f"| index build (1 / 10 threads) | {s(sp1)} / **{s(sp10)}** ({mb(max(spm1 or 0, spm10 or 0))}) | "
               f"{s(mb1)} / {s(mb10)} ({mb(mbm10)}) | {s(bw1)} (one thread) |")
    out += ["", "| f-SMEMs, 1 M reads vs chr21, L=19 | sporos `-s f` (1 / 10 threads) | minibwa `fastmap -s f` | "
                "ropebwt3 `mem -c f` (1 thread) |", "|---|---|---|---|"]
    for f in (2, 5):
        t = T.smem_times("150bp", f)
        out.append(f"| f = {f} | **{s(t['s1'])} / {s(t['s10'])}** | {s(t['mb'])} | {s(t['rb1'])} |")
    out += ["", "| Frequency-capped matches, L=50, both strands | sporos `mem --max-occ N -b` (1 / 10 threads) | "
                "Mumemto `-k 2 -f N` / `-f 0 -F N` |", "|---|---|---|"]
    for pair, label in (("ecoli", "E. coli K-12 vs O157"), ("chr21", "hg38 vs chimp chr21")):
        rows = [r for r in T.CAPPED if r["pair"] == pair]
        caps = sorted({r["mode"].split()[-1] for r in rows}, key=int)

        def rng(tool, t, field="wall_s"):
            vals = [float(r[field]) for r in rows if r["tool"] == tool and r["threads"] == str(t)
                    and r[field].replace(".", "", 1).isdigit()]
            return (min(vals), max(vals)) if vals else (None, None)

        def srng(lo_hi):
            lo, hi = lo_hi
            return s(lo) if s(lo) == s(hi) else f"{s(lo)[:-2]}–{s(hi)}"

        out.append(f"| {label}, N = {' or '.join(caps)} | **{srng(rng('sporos', 1))} / {srng(rng('sporos', 10))}** "
                   f"({mb(rng('sporos', 1, 'max_rss_mb')[1])}) | {srng(rng('mumemto', 1))} "
                   f"({mb(rng('mumemto', 1, 'max_rss_mb')[1])}) |")
    out += ["", "| Seeds, CHM13 chr1, ns/bp (forward / canonical) | sporos (16 lanes) | simd-minimizers | minimap2 `mm_sketch` |",
            "|---|---|---|---|"]
    for mode, k, w, label in (("", 21, 11, "minimizers k=21 w=11"), ("", 31, 19, "minimizers k=31 w=19"),
                              ("closed-", 21, 11, "closed syncmers k=21 w=11"), ("closed-", 31, 19, "closed syncmers k=31 w=19")):
        for T_ in (1, 10):
            sp = [T.seed(T.SP16, "chm13_chr1", mode + m, k, w, T_) for m in ("fwd", "can")]
            sm = [T.seed(T.SMIN, "chm13_chr1", mode + m, k, w, T_) for m in ("fwd", "can")]
            mm = T.mm2_seed("chm13_chr1", "can", k, w, T_) if not mode else None
            fmt = lambda a: " / ".join("–" if v is None else f"{v:.2f}" for v in a)
            out.append(f"| {label}, {T_} thread{'s' if T_ > 1 else ''} | {fmt(sp)} | {fmt(sm)} | "
                       f"{'— / ' + format(mm, '.2f') + (' (k=28)' if k == 31 else '') if mm else ''} |")
    sp = lambda m: T.lookup("sporos k=19 w=12", m)
    out += ["",
            "- **Seed-index lookups** (hg38 chr21 index, k=19, w=12, occurrence cap 20, 1 M simulated 150 bp reads, every",
            "  hit returned):",
            f"  - read minimizers: {sp('lookup minimizers, batch'):.0f} Mbp/s per thread batched "
            f"({sp('lookup minimizers, per read'):.0f} Mbp/s one read per call), "
            f"{sp('lookup minimizers, batch, 8 threads'):.0f} Mbp/s on 8 threads; minimap2's `mm_sketch` + `mm_idx_get` "
            f"on the same minimizers: {V['lookup_mm2']} Mbp/s; strobealign's randstrobe `find_hits`: {V['lookup_sa']} Mbp/s;",
            f"  - every read k-mer: {sp('lookup all k-mers, batch'):.0f} Mbp/s per thread;",
            f"  - index build: {sp('build, 1 threads'):.2f} s on one thread (minimap2 "
            f"{T.lookup('minimap2 k=19 w=12', 'build, 1 threads'):.2f} s, strobealign "
            f"{T.lookup('strobealign defaults', 'build, 1 threads'):.2f} s)."]
    cli = {r["threads"]: r for r in T.tsv("seeds_cli.tsv")}
    if cli:
        out.append(f"- **End to end**, `sporos minimizers -k 21 -w 11` on CHM13 chr1, reading the FASTA and writing every "
                   f"line: {s(float(cli['1']['wall_s']))} on one thread, {s(float(cli['10']['wall_s']))} on ten.")
    peak = lambda L: T.rss(T.MEMS, pair="chr1", L=L, tool="sporos", mode="maxmatch", threads=1)
    out.append(f"- **Peak memory** for MEMs on chr1 is {mb(peak('100'))} at L=100 and {mb(peak('50'))} at L=50; copMEM2 "
               f"`-mf` uses {mb(T.rss(T.MEMS, pair='chr1', L='50', tool='copmem2 -mf', mode='maxmatch', threads=1))}.")
    c = T.smem_check("150bp", 1, "sporos vs minibwa")
    c10 = T.smem_check("10kb", 1, "sporos vs minibwa")
    out += ["- **SMEMs:**",
            "  - minibwa and `bwa fastmap` run on one thread.",
            f"  - sporos peaks at {mb(T.rss(T.SMEMS, reads='150bp', tool='sporos', mode='SMEM', threads=1))} on 1 M reads "
            f"(minibwa: {mb(T.rss(T.SMEMS, reads='150bp', tool='minibwa', mode='SMEM', threads=1))}), as it reads the "
            "whole query file at once."]
    if c and c10:
        out.append(f"  - Its SMEMs equal minibwa's except where minibwa matches through reference `N`s: "
                   f"{int(c['only_other']):,} of minibwa's {int(c['n_other']) / 1e6:.2f} M SMEMs on the short reads "
                   f"(and {int(c['only_sporos'])} of sporos's are not in minibwa's), and {int(c10['only_other']):,} on the long reads.")
    return "\n".join(out)


def update():
    text = open(README).read()
    for name, block in (("speed", speed()), ("chart", chart()), ("results", results())):
        begin, end = f"<!-- {name}:begin (generated by pixi run paper; edit tools/paper/readme.py) -->", f"<!-- {name}:end -->"
        pattern = re.compile(re.escape(begin) + r".*?" + re.escape(end), re.S)
        if not pattern.search(text):
            raise SystemExit(f"README.md has no {name} markers")
        text = pattern.sub(lambda _: f"{begin}\n{block}\n{end}", text)
    open(README, "w").write(text)


if __name__ == "__main__":
    update()
    print("README.md: speed, chart and results blocks rewritten")
