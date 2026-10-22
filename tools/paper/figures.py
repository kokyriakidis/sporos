"""The paper's figures, each defined once and drawn twice: TikZ for paper/ and SVG for the HTML draft.

  overview   the library: input, the three engines and what each reports, the interfaces
  speedups   per task, sporos's speedup over the fastest other tool, one and ten threads (log scale)
"""
import math

from tables import SPEEDUPS

WIDE = {"overview": True, "speedups": False}

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
          ("seedindex", "api"), ("sweeps", "api"), ("unique", "api")]
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
# One row per task: filled marker at one thread, hollow at ten (each against the other tool's best at that
# thread count); the other tool's name at the right. Log-scale x from XMIN to XMAX.
XMIN, XMAX = 1, 1000
TICKS = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]


def _speedups_layout():
    rows = SPEEDUPS
    return rows, len(rows)


def _xfrac(v):
    return (math.log10(v) - math.log10(XMIN)) / (math.log10(XMAX) - math.log10(XMIN))


def _speedups_tikz():
    rows, n = _speedups_layout()
    width, step, left = 6.0, 0.42, 0.0
    height = n * step
    out = [f"\\draw[black!60] (0,0) -- ({width},0);"]
    for t in TICKS:
        x = _xfrac(t) * width
        out.append(f"\\draw[gray!25] ({x:.3f},0) -- ({x:.3f},{height:.3f});"
                   f"\\node[below, font=\\scriptsize] at ({x:.3f},0) {{{t}$\\times$}};")
    out.append(f"\\node[font=\\small] at ({width / 2:.2f},-0.65) {{speedup over the fastest other tool}};")
    for i, r in enumerate(rows):
        y = height - (i + 0.5) * step
        out.append(f"\\node[anchor=east, font=\\scriptsize] at ({left - 0.1:.2f},{y:.3f}) {{{_tex_text(r['label'])}}};")
        out.append(f"\\node[anchor=west, font=\\scriptsize\\itshape, black!70] at ({width + 0.1:.2f},{y:.3f}) {{{r['rival']}}};")
        xs = [v for v in (r["t1"], r["t10"]) if v]
        if len(xs) == 2:
            out.append(f"\\draw[black!40] ({_xfrac(min(xs)) * width:.3f},{y:.3f}) -- ({_xfrac(max(xs)) * width:.3f},{y:.3f});")
        if r["t1"]:
            out.append(f"\\fill[black] ({_xfrac(r['t1']) * width:.3f},{y:.3f}) circle (2.2pt);")
        if r["t10"]:
            out.append(f"\\draw[black, fill=white, thick] ({_xfrac(r['t10']) * width:.3f},{y:.3f}) circle (2.2pt);")
    out.append(f"\\fill[black] (0.2,{height + 0.35:.3f}) circle (2.2pt) node[right, font=\\scriptsize] {{1 thread}};"
               f"\\draw[black, fill=white, thick] (2.0,{height + 0.35:.3f}) circle (2.2pt) node[right, font=\\scriptsize] {{10 threads}};")
    return "\n".join(out)


def _speedups_svg():
    rows, n = _speedups_layout()
    left, width, step, top = 210, 360, 26, 34
    W, H = left + width + 130, top + n * step + 50

    def X(v):
        return left + _xfrac(v) * width

    out = [f'<svg viewBox="0 0 {W} {H}" class="chart" role="img" aria-label="Speedup of sporos over the fastest other tool">',
           f'<circle cx="{left + 6}" cy="14" r="5" fill="currentColor"/><text x="{left + 16}" y="18" font-size="11" fill="currentColor">1 thread</text>',
           f'<circle cx="{left + 96}" cy="14" r="5" fill="none" stroke="currentColor" stroke-width="1.6"/>'
           f'<text x="{left + 106}" y="18" font-size="11" fill="currentColor">10 threads</text>']
    base = top + n * step
    for t in TICKS:
        out.append(f'<line x1="{X(t):.1f}" y1="{top}" x2="{X(t):.1f}" y2="{base}" stroke="currentColor" stroke-opacity=".15"/>'
                   f'<text x="{X(t):.1f}" y="{base + 16}" text-anchor="middle" font-size="11" fill="currentColor">{t}×</text>')
    out.append(f'<line x1="{left}" y1="{base}" x2="{left + width}" y2="{base}" stroke="currentColor" stroke-opacity=".6"/>'
               f'<text x="{left + width / 2:.0f}" y="{base + 36}" text-anchor="middle" font-size="12" fill="currentColor">'
               'speedup over the fastest other tool</text>')
    for i, r in enumerate(rows):
        y = top + (i + 0.5) * step
        out.append(f'<text x="{left - 8}" y="{y + 4:.1f}" text-anchor="end" font-size="11" fill="currentColor">{r["label"]}</text>'
                   f'<text x="{left + width + 8}" y="{y + 4:.1f}" font-size="11" font-style="italic" fill="currentColor" '
                   f'fill-opacity=".75">{r["rival"]}</text>')
        xs = [v for v in (r["t1"], r["t10"]) if v]
        if len(xs) == 2:
            out.append(f'<line x1="{X(min(xs)):.1f}" y1="{y:.1f}" x2="{X(max(xs)):.1f}" y2="{y:.1f}" stroke="currentColor" stroke-opacity=".35"/>')
        if r["t1"]:
            out.append(f'<circle cx="{X(r["t1"]):.1f}" cy="{y:.1f}" r="5" fill="currentColor"><title>{r["label"]}, 1 thread: '
                       f'{r["t1"]:.1f}× over {r["rival"]}</title></circle>')
        if r["t10"]:
            out.append(f'<circle cx="{X(r["t10"]):.1f}" cy="{y:.1f}" r="5" fill="var(--paper, #fff)" stroke="currentColor" '
                       f'stroke-width="1.6"><title>{r["label"]}, 10 threads: {r["t10"]:.1f}× over {r["rival"]}</title></circle>')
    out.append("</svg>")
    return "\n".join(out)


def tikz(key):
    return {"overview": _overview_tikz, "speedups": _speedups_tikz}[key]()


def svg(key):
    return {"overview": _overview_svg, "speedups": _speedups_svg}[key]()
