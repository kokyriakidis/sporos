"""The paper's figures, each defined once and drawn twice: TikZ for paper/ and SVG for the HTML draft.

  overview   the library: input, the three engines and what each reports, the interfaces
  speedups   per task, sporos's speedup over the fastest other tool, one and ten threads (log scale)
"""
import math

from tables import SPEEDUPS

WIDE = {"overview": True, "speedups": True}

TIKZ_STYLES = "\\usetikzlibrary{plotmarks}"

# ------------------------------------------------------------------------------- overview
# Coordinates in cm, y upwards. A lane is (x0, y0, x1, y1, fill, label); a box is (x, y-centre, width, height,
# kind, lines), anchored at its left edge.
LANES = {
    "kmers": (3.0, 3.55, 17.0, 5.15, "E3EEF9", "sampled k-mers"),
    "mems": (3.0, 1.75, 17.0, 3.35, "FCEBDD", "exact matches between genomes"),
    "smems": (3.0, -0.05, 17.0, 1.55, "E5F3E8", "exact matches of reads"),
}
BOXES = {
    "input": (0.0, 2.55, 2.6, 1.3, "data", ["FASTA/FASTQ", "plain, gzip, BGZF", "→ 2-bit codes"]),
    "window": (3.3, 4.25, 4.2, 0.95, "step", ["Minimizers, syncmers:", "SIMD lanes, sliding minimum"]),
    "seedindex": (8.0, 4.25, 4.2, 0.95, "step", ["Seed index: canonical", "k-mer hash, fingerprints"]),
    "lookups": (12.7, 4.25, 4.0, 0.95, "step", ["Batched lookups:", "one prefetch pipeline"]),
    "anchors": (3.3, 2.45, 4.2, 0.95, "step", ["Sparse anchor index:", "k1 · k2 ≤ L − K + 1"]),
    "extend": (8.0, 2.45, 4.2, 0.95, "step", ["MEMs: each found once,", "32-base extension"]),
    "sweeps": (12.7, 2.45, 4.0, 0.95, "step", ["Containment sweeps:", "MUMs, capped MEMs"]),
    "fmd": (3.3, 0.65, 4.2, 0.95, "step", ["FMD-index: libsais,", "128-row blocks, k-mer cache"]),
    "search": (8.0, 0.65, 4.2, 0.95, "step", ["SMEMs, f-SMEMs: 16 reads", "interleaved, prefetched"]),
    "unique": (12.7, 0.65, 4.0, 0.95, "step", ["Unique stretches:", "located, then compared"]),
    "api": (17.5, 2.45, 2.6, 1.3, "out", ["Mojo API, CLI,", "C and Python"]),
}
ARROWS = [("input", "window"), ("input", "anchors"), ("input", "fmd"), ("window", "seedindex"),
          ("anchors", "extend"), ("extend", "sweeps"), ("fmd", "search"), ("search", "unique"),
          ("seedindex", "lookups"), ("lookups", "api"), ("sweeps", "api"), ("unique", "api")]
FILL = {"data": "EDEDED", "step": "FFFFFF", "out": "FBE3D6"}


def _turn(a, b):
    """Where an elbow arrow turns: midway, or just before the interfaces box so as to pass no other box."""
    (xa, _), (xb, _) = _anchor(a, "right"), _anchor(b, "left")
    return xb - 0.25 if b == "api" else (xa + xb) / 2


def _anchor(box, side):
    x, y, w, h, *_ = BOXES[box]
    return (x + w, y) if side == "right" else (x, y)


def _overview_tikz():
    out = []
    for x0, y0, x1, y1, fill, label in LANES.values():
        out.append(f"\\fill[fill={{rgb,255:red,{int(fill[0:2], 16)};green,{int(fill[2:4], 16)};blue,{int(fill[4:6], 16)}}}, "
                   f"rounded corners=3pt] ({x0},{y0}) rectangle ({x1},{y1});"
                   f"\\node[anchor=north west, font=\\scriptsize\\itshape, inner sep=2pt] at ({x0 + 0.05},{y1 - 0.02}) {{{label}}};")
    for key, (x, y, w, h, kind, lines) in BOXES.items():
        fill = FILL[kind]
        rgb = f"{{rgb,255:red,{int(fill[0:2], 16)};green,{int(fill[2:4], 16)};blue,{int(fill[4:6], 16)}}}"
        style = {"data": "draw=black!60", "step": "draw=black!70", "out": "draw=black!80, thick"}[kind]
        text = "\\\\".join(_tex_text(line) for line in lines)
        out.append(f"\\node[{style}, fill={rgb}, rounded corners=2pt, align=center, font=\\footnotesize, inner sep=2pt, "
                   f"minimum width={w}cm, minimum height={h}cm, anchor=west] ({key}) at ({x},{y}) {{{text}}};")
    for a, b in ARROWS:
        (xa, ya), (xb, yb) = _anchor(a, "right"), _anchor(b, "left")
        if abs(ya - yb) < 0.01:
            out.append(f"\\draw[-latex, thick, black!70] ({xa:.2f},{ya:.2f}) -- ({xb:.2f},{yb:.2f});")
        else:
            mid = _turn(a, b)
            out.append(f"\\draw[-latex, thick, black!70] ({xa:.2f},{ya:.2f}) -- ({mid:.2f},{ya:.2f}) -- "
                       f"({mid:.2f},{yb:.2f}) -- ({xb:.2f},{yb:.2f});")
    return "\n".join(out)


def _tex_text(line):
    return (line.replace("·", "$\\cdot$").replace("≤", "$\\le$").replace("−", "$-$").replace("→", "$\\to$")
            .replace("k1", "$k_1$").replace("k2", "$k_2$"))


def _overview_svg():
    s = 44
    x_max, y_min, y_max = 20.3, -0.2, 5.3
    W, H = (x_max + 0.3) * s, (y_max - y_min) * s

    def X(x):
        return (x + 0.15) * s

    def Y(y):
        return (y_max - y) * s

    out = [f'<svg viewBox="0 0 {W:.0f} {H:.0f}" class="chart" role="img" aria-label="The sporos library">',
           '<defs><marker id="ah" markerWidth="8" markerHeight="8" refX="7" refY="3" orient="auto">'
           '<path d="M0,0 L7,3 L0,6 z" fill="currentColor"/></marker></defs>']
    for x0, y0, x1, y1, fill, label in LANES.values():
        out.append(f'<rect x="{X(x0):.1f}" y="{Y(y1):.1f}" width="{(x1 - x0) * s:.1f}" height="{(y1 - y0) * s:.1f}" '
                   f'rx="5" fill="#{fill}"/><text x="{X(x0 + 0.12):.1f}" y="{Y(y1 - 0.33):.1f}" font-size="10" '
                   f'font-style="italic" fill="#333">{label}</text>')
    for key, (x, y, w, h, kind, lines) in BOXES.items():
        stroke = {"data": "#777", "step": "#555", "out": "#222"}[kind]
        width = 2 if kind == "out" else 1
        out.append(f'<rect x="{X(x):.1f}" y="{Y(y + h / 2):.1f}" width="{w * s:.1f}" height="{h * s:.1f}" rx="4" '
                   f'fill="#{FILL[kind]}" stroke="{stroke}" stroke-width="{width}"/>')
        for i, line in enumerate(lines):
            ty = Y(y) + (i - (len(lines) - 1) / 2) * 13 + 4
            text = line.replace("k1", "k₁").replace("k2", "k₂")
            out.append(f'<text x="{X(x + w / 2):.1f}" y="{ty:.1f}" text-anchor="middle" font-size="11" '
                       f'fill="#111">{text}</text>')
    for a, b in ARROWS:
        (xa, ya), (xb, yb) = _anchor(a, "right"), _anchor(b, "left")
        mid = _turn(a, b)
        pts = [(xa, ya), (xb, yb)] if abs(ya - yb) < 0.01 else [(xa, ya), (mid, ya), (mid, yb), (xb, yb)]
        path = " ".join(f"{X(px):.1f},{Y(py):.1f}" for px, py in pts[:-1]) + f" {X(pts[-1][0]) - 2:.1f},{Y(pts[-1][1]):.1f}"
        out.append(f'<polyline points="{path}" fill="none" stroke="currentColor" stroke-width="1.3" marker-end="url(#ah)"/>')
    out.append("</svg>")
    return "\n".join(out)


# ------------------------------------------------------------------------------- speedups
# A dot plot across the page: rows grouped by kind of seed, each labelled with its task and the fastest other tool;
# a filled dot at one thread and a hollow one at ten on a log axis from XMIN to XMAX, with parity (1x) dashed and
# the slower side shaded; the speedups themselves in two columns at the right, whose headings carry the legend.
XMIN, XMAX = 0.5, 1000
TICKS = [(0.5, "0.5"), (1, "1"), (2, "2"), (5, "5"), (10, "10"), (20, "20"), (50, "50"), (100, "100"), (200, "200"),
         (500, "500"), (1000, "1,000")]
# Geometry in cm: plot width, row and group-heading heights, value columns to the right of the plot.
PW, ROW, HEAD, COL1, COL2 = 8.6, 0.46, 0.55, 1.3, 3.2


def _xfrac(v):
    v = min(max(v, XMIN), XMAX)
    return (math.log10(v) - math.log10(XMIN)) / (math.log10(XMAX) - math.log10(XMIN))


def _times(v):
    if not v:
        return "–"
    return (f"{v:.2f}" if v < 1 else f"{v:.1f}" if v < 10 else f"{v:,.0f}") + "×"


def _times_tex(v):
    return _times(v).replace("×", "$\\times$")


def _speedup_rows():
    """(y of each row's centre, from the top, in cm; the row or group heading) top to bottom, and the total height."""
    out, y, group = [], 0.0, None
    for r in SPEEDUPS:
        if r["group"] != group:
            group = r["group"]
            y += HEAD
            out.append((y - HEAD / 2 + 0.06, "group", group))
        y += ROW
        out.append((y - ROW / 2, "row", r))
    return out, y


def _speedups_tikz():
    rows, height = _speedup_rows()
    top = height
    out = [f"\\fill[gray!12] (0,0) rectangle ({_xfrac(1) * PW:.3f},{top:.3f});",
           f"\\node[anchor=south, font=\\scriptsize, black!60] at ({_xfrac(0.7) * PW:.3f},{top + 0.02:.3f}) {{slower}};"]
    for v, lab in TICKS:
        x = _xfrac(v) * PW
        out.append(f"\\draw[gray!25] ({x:.3f},0) -- ({x:.3f},{top:.3f});"
                   f"\\node[below, font=\\scriptsize] at ({x:.3f},0) {{{lab}$\\times$}};")
    out.append(f"\\draw[black!55, dashed] ({_xfrac(1) * PW:.3f},0) -- ({_xfrac(1) * PW:.3f},{top:.3f});")
    out.append(f"\\draw[black!60] (0,0) -- ({PW},0);")
    out.append(f"\\node[font=\\footnotesize] at ({PW / 2:.2f},-0.75) {{speedup of sporos over the fastest other tool (log scale)}};")
    hy = top + 0.3
    # Each value column's heading, right-aligned with its values, with its marker to the left of the text.
    out.append(f"\\node[anchor=east, font=\\footnotesize\\bfseries] (h1) at ({PW + COL1 + 0.35:.3f},{hy:.3f}) {{1 thread}};"
               f"\\fill[black] ([xshift=-5pt]h1.west) circle (2.4pt);"
               f"\\node[anchor=east, font=\\footnotesize\\bfseries] (h10) at ({PW + COL2 + 0.35:.3f},{hy:.3f}) {{10 threads}};"
               f"\\draw[black, fill=white, thick] ([xshift=-5pt]h10.west) circle (2.4pt);")
    for yc, kind, item in rows:
        y = top - yc
        if kind == "group":
            out.append(f"\\node[anchor=west, font=\\footnotesize\\bfseries] at (-5.6,{y:.3f}) {{{item}}};")
            continue
        r = item
        out.append(f"\\node[anchor=east, font=\\footnotesize] at (-0.12,{y:.3f}) "
                   f"{{{_tex_text(r['label'])} \\textcolor{{black!60}}{{vs {r['rival']}}}}};")
        xs = [v for v in (r["t1"], r["t10"]) if v]
        if len(xs) == 2:
            out.append(f"\\draw[black!45, line width=0.8pt] ({_xfrac(min(xs)) * PW:.3f},{y:.3f}) -- ({_xfrac(max(xs)) * PW:.3f},{y:.3f});")
        if r["t10"]:
            out.append(f"\\draw[black, fill=white, thick] ({_xfrac(r['t10']) * PW:.3f},{y:.3f}) circle (2.4pt);")
        if r["t1"]:
            out.append(f"\\fill[black] ({_xfrac(r['t1']) * PW:.3f},{y:.3f}) circle (2.4pt);")
        out.append(f"\\node[anchor=east, font=\\footnotesize] at ({PW + COL1 + 0.35:.3f},{y:.3f}) {{{_times_tex(r['t1'])}}};"
                   f"\\node[anchor=east, font=\\footnotesize] at ({PW + COL2 + 0.35:.3f},{y:.3f}) {{{_times_tex(r['t10'])}}};")
    return "\n".join(out)


def _speedups_svg():
    rows, height = _speedup_rows()
    s = 40  # px per cm
    left, top = 5.8 * s, 1.0 * s
    W, H = left + (PW + COL2 + 0.6) * s, top + height * s + 1.4 * s

    def X(v):
        return left + _xfrac(v) * PW * s

    base = top + height * s
    out = [f'<svg viewBox="0 0 {W:.0f} {H:.0f}" class="chart" role="img" aria-label="Speedup of sporos over the fastest other tool">',
           f'<rect x="{left:.1f}" y="{top:.1f}" width="{X(1) - left:.1f}" height="{height * s:.1f}" fill="currentColor" fill-opacity=".07"/>',
           f'<text x="{X(0.7):.1f}" y="{top - 4:.1f}" text-anchor="middle" font-size="11" fill="currentColor" fill-opacity=".75">slower</text>']
    for v, lab in TICKS:
        out.append(f'<line x1="{X(v):.1f}" y1="{top:.1f}" x2="{X(v):.1f}" y2="{base:.1f}" stroke="currentColor" stroke-opacity=".12"/>'
                   f'<text x="{X(v):.1f}" y="{base + 15:.1f}" text-anchor="middle" font-size="11" fill="currentColor">{lab}×</text>')
    out.append(f'<line x1="{X(1):.1f}" y1="{top:.1f}" x2="{X(1):.1f}" y2="{base:.1f}" stroke="currentColor" stroke-opacity=".6" stroke-dasharray="4 3"/>'
               f'<line x1="{left:.1f}" y1="{base:.1f}" x2="{left + PW * s:.1f}" y2="{base:.1f}" stroke="currentColor" stroke-opacity=".6"/>'
               f'<text x="{left + PW * s / 2:.1f}" y="{base + 36:.1f}" text-anchor="middle" font-size="12" fill="currentColor">'
               'speedup of sporos over the fastest other tool (log scale)</text>')
    c1, c2 = left + (PW + COL1 + 0.35) * s, left + (PW + COL2 + 0.35) * s
    out.append(f'<circle cx="{c1 - 76:.1f}" cy="{top - 18:.1f}" r="5" fill="currentColor"/>'
               f'<text x="{c1:.1f}" y="{top - 14:.1f}" text-anchor="end" font-size="12" font-weight="700" fill="currentColor">1 thread</text>'
               f'<circle cx="{c2 - 88:.1f}" cy="{top - 18:.1f}" r="5" fill="none" stroke="currentColor" stroke-width="1.6"/>'
               f'<text x="{c2:.1f}" y="{top - 14:.1f}" text-anchor="end" font-size="12" font-weight="700" fill="currentColor">10 threads</text>')
    for yc, kind, item in rows:
        y = top + yc * s
        if kind == "group":
            out.append(f'<text x="8" y="{y + 4:.1f}" font-size="12" font-weight="700" fill="currentColor">{item}</text>')
            continue
        r = item
        out.append(f'<text x="{left - 6:.1f}" y="{y + 4:.1f}" text-anchor="end" font-size="12" fill="currentColor">{r["label"]} '
                   f'<tspan fill-opacity=".6">vs {r["rival"]}</tspan></text>')
        xs = [v for v in (r["t1"], r["t10"]) if v]
        if len(xs) == 2:
            out.append(f'<line x1="{X(min(xs)):.1f}" y1="{y:.1f}" x2="{X(max(xs)):.1f}" y2="{y:.1f}" stroke="currentColor" stroke-opacity=".45" stroke-width="1.5"/>')
        if r["t10"]:
            out.append(f'<circle cx="{X(r["t10"]):.1f}" cy="{y:.1f}" r="5" fill="var(--paper, #fff)" stroke="currentColor" stroke-width="1.6">'
                       f'<title>{r["label"]}, 10 threads: {_times(r["t10"])} over {r["rival"]}</title></circle>')
        if r["t1"]:
            out.append(f'<circle cx="{X(r["t1"]):.1f}" cy="{y:.1f}" r="5" fill="currentColor"><title>{r["label"]}, 1 thread: '
                       f'{_times(r["t1"])} over {r["rival"]}</title></circle>')
        out.append(f'<text x="{c1:.1f}" y="{y + 4:.1f}" text-anchor="end" font-size="12" fill="currentColor">{_times(r["t1"])}</text>'
                   f'<text x="{c2:.1f}" y="{y + 4:.1f}" text-anchor="end" font-size="12" fill="currentColor">{_times(r["t10"])}</text>')
    out.append("</svg>")
    return "\n".join(out)


def tikz(key):
    return {"overview": _overview_tikz, "speedups": _speedups_tikz}[key]()


def svg(key):
    return {"overview": _overview_svg, "speedups": _speedups_svg}[key]()
