"""Markdown tables from results/*.tsv: `python3 -I tables.py [results_dir]`."""

import csv
import sys
from collections import defaultdict
from pathlib import Path

RES = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / "results"


def num(x):
    try:
        return float(x)
    except ValueError:
        return None


def fmt(r):
    wall, user = num(r["wall_s"]), num(r["user_s"])
    if wall is None:
        return r["wall_s"]
    return f"{wall:.2f} ({user:.2f}u)"


def rows(name):
    path = RES / name
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def mems_tables():
    data = rows("mems.tsv")
    out = []
    builds = [r for r in data if r["mode"] == "build"]
    for pair in dict.fromkeys(r["pair"] for r in data):
        for L in dict.fromkeys(r["L"] for r in data if r["pair"] == pair and r["L"] != "-"):
            out.append(f"\n**{pair}, L={L}** (best wall s; user s; max RSS MB)\n")
            out.append("| tool | T | -maxmatch | -mum | -mumreference | other | max RSS MB |")
            out.append("|---|---:|---:|---:|---:|---|---:|")
            cells = defaultdict(dict)
            rss = defaultdict(int)
            sporos = {}
            for r in data:
                if r["pair"] != pair or r["L"] != L:
                    continue
                key = (r["tool"], r["threads"])
                cell = fmt(r)
                if r["tool"] == "sporos":
                    sporos[(r["mode"], r["threads"])] = num(r["wall_s"])
                cells[key][r["mode"]] = (cell, r["mode"], num(r["wall_s"]), r["threads"])
                rss[key] = max(rss[key], int(num(r["max_rss_mb"]) or 0))
            for key, modes in cells.items():
                cols = []
                other = []
                for m in ("maxmatch", "mum", "mumreference"):
                    if m in modes:
                        cell, _, wall, t = modes[m]
                        ref = sporos.get((m, t))
                        if key[0] != "sporos" and ref and wall:
                            cell += f" = {wall / ref:.1f}x sporos"
                        cols.append(cell)
                    else:
                        cols.append("")
                for m, (cell, _, _, _) in modes.items():
                    if m not in ("maxmatch", "mum", "mumreference"):
                        other.append(f"{m}: {cell}")
                out.append(f"| {key[0]} | {key[1]} | " + " | ".join(cols) + f" | {'; '.join(other)} | {rss[key] or ''} |")
        b = [r for r in builds if r["pair"] == pair]
        if b:
            out.append(
                "\nIndex builds: "
                + "; ".join(f"{r['tool']} -t{r['threads']} {float(r['wall_s']):.2f} s ({r['max_rss_mb']} MB)" for r in b)
            )
    return "\n".join(out)


def checks_table():
    data = rows("mems_checks.tsv")
    out = ["| pair | L | check | n (a) | n (b) | only a | only b | result |", "|---|---:|---|---:|---:|---:|---:|---|"]
    for r in data:
        out.append(
            f"| {r['pair']} | {r['L']} | {r['check']} | {r['n_a']} | {r['n_b']} | {r['only_a']} | {r['only_b']} | {r['result']} |"
        )
    return "\n".join(out)


def seeds_table():
    data = rows("seeds.tsv")
    cols = list(dict.fromkeys((r["tool"], r["input"], r["threads"]) for r in data))
    groups = defaultdict(dict)
    for r in data:
        groups[(r["mode"], r["k"], r["w"])][(r["tool"], r["input"], r["threads"])] = r["ns_per_bp"]
    out = []
    for inp, t in dict.fromkeys((c[1], c[2]) for c in cols):
        sub = [c for c in cols if c[1] == inp and c[2] == t]
        out.append(f"\n**{inp}, {t} thread{'s' if t != '1' else ''}** (ns/bp, lower is better)\n")
        out.append("| mode | k | w | " + " | ".join(c[0] for c in sub) + " |")
        out.append("|---|---:|---:|" + "---:|" * len(sub))
        for (mode, k, w), vals in groups.items():
            if not any(c in vals for c in sub):
                continue
            out.append(f"| {mode} | {k} | {w} | " + " | ".join(vals.get(c, "") for c in sub) + " |")
    cli = rows("seeds_cli.tsv")
    if cli:
        out.append("\n**sporos minimizers CLI, CHM13 chr1** (`-k 21 -w 11`, output to /dev/null)\n")
        out.append("| threads | wall s | user s | max RSS MB |")
        out.append("|---:|---:|---:|---:|")
        for r in cli:
            out.append(f"| {r['threads']} | {r['wall_s']} | {r['user_s']} | {r['max_rss_mb']} |")
    return "\n".join(out)


if __name__ == "__main__":
    print("## MEMs and MUMs")
    print(mems_tables())
    print("\n## Correctness checks\n")
    print(checks_table())
    print("\n## Minimizers and syncmers")
    print(seeds_table())
