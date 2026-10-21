"""Draws the README's benchmark chart as light and dark SVGs.

    python3 scripts/benchmark_chart.py      # writes docs/benchmarks-{light,dark}.svg

The numbers are the README's Performance tables (xaros, Intel i9-7900X); update both together.
Each panel is one benchmark on its own scale, lower is better; sporos is the accent, rivals gray.
"""

from pathlib import Path

PANELS = [
    (
        "MEMs, hg38 vs CHM13 chr1, L = 50",
        "1 thread, wall time in seconds",
        [("sporos", 19.5), ("copMEM2", 46.5), ("MUMmer4", 159.6)],
        " s",
    ),
    (
        "MEMs with 10 threads, same pair",
        "wall time in seconds (MUMmer4 runs single-threaded)",
        [("sporos", 3.3), ("copMEM2", 46.4)],
        " s",
    ),
    (
        "MUMs (-mum), same pair",
        "1 thread, wall time in seconds",
        [("sporos", 21.6), ("MUMmer4", 87.1)],
        " s",
    ),
    (
        "SMEMs of 1 M 150 bp reads against hg38 chr21, L = 19",
        "1 thread, wall time in seconds (index built beforehand)",
        [("sporos", 1.85), ("minibwa", 8.85), ("bwa fastmap", 27.3)],
        " s",
    ),
    (
        "Canonical minimizers, CHM13 chr1, k = 21, w = 11",
        "1 thread, nanoseconds per base",
        [("sporos", 1.05), ("simd-minimizers", 2.13), ("minimap2", 17.1)],
        " ns",
    ),
]

THEMES = {
    # Transparent surface: the page is GitHub's (#ffffff light, #0d1117 dark).
    "light": {"accent": "#2a78d6", "context": "#898781", "ink": "#0b0b0b", "ink2": "#52514e", "rule": "#e1e0d9"},
    "dark": {"accent": "#3987e5", "context": "#898781", "ink": "#f0efec", "ink2": "#c3c2b7", "rule": "#2c2c2a"},
}

WIDTH = 820
LABEL_W = 150  # tool names
VALUE_W = 70  # value at the bar tip
BAR_H = 18
ROW = 26
FONT = 'system-ui, -apple-system, "Segoe UI", Helvetica, Arial, sans-serif'


def _bar(x, y, length, height, radius):
    """A bar from the baseline at x: square there, rounded at its data end."""
    r = min(radius, length / 2, height / 2)
    return (
        f"M{x:.1f},{y:.1f} H{x + length - r:.1f} A{r},{r} 0 0 1 {x + length:.1f},{y + r:.1f} "
        f"V{y + height - r:.1f} A{r},{r} 0 0 1 {x + length - r:.1f},{y + height:.1f} H{x:.1f} Z"
    )


def _value(v, unit):
    if unit == " s":
        return (f"{v:.1f}" if v < 100 else f"{v:.0f}") + unit
    return (f"{v:.2f}" if v < 10 else f"{v:.1f}") + unit


def render(theme):
    t = THEMES[theme]
    parts = []
    y = 8
    area = WIDTH - LABEL_W - VALUE_W
    for title, subtitle, rows, unit in PANELS:
        parts.append(f'<text x="0" y="{y + 14}" class="title">{title}</text>')
        parts.append(f'<text x="0" y="{y + 31}" class="sub">{subtitle} · lower is better</text>')
        y += 44
        top = max(v for _, v in rows)
        best_rival = min(rows[1:], key=lambda row: row[1])
        for name, v in rows:
            length = max(2.0, area * v / top)
            fill = t["accent"] if name == "sporos" else t["context"]
            cls = "name strong" if name == "sporos" else "name"
            parts.append(f'<text x="{LABEL_W - 12}" y="{y + 13}" class="{cls}" text-anchor="end">{name}</text>')
            parts.append(f'<path d="{_bar(LABEL_W, y, length, BAR_H, 4)}" fill="{fill}"/>')
            label = _value(v, unit)
            parts.append(f'<text x="{LABEL_W + length + 8}" y="{y + 13}" class="val">{label}')
            if name == "sporos":
                factor = best_rival[1] / v
                parts.append(f'<tspan class="note" dx="10">{factor:.1f}× faster than {best_rival[0]}</tspan>')
            parts.append("</text>")
            y += ROW
        # The panel's baseline, a hairline.
        first = y - ROW * len(rows) - 4
        parts.append(f'<line x1="{LABEL_W}" y1="{first}" x2="{LABEL_W}" y2="{y - 4}" stroke="{t["rule"]}" stroke-width="1"/>')
        y += 22
    height = y - 14
    style = (
        f".title{{font:600 15px {FONT};fill:{t['ink']}}}"
        f".sub{{font:400 12.5px {FONT};fill:{t['ink2']}}}"
        f".name{{font:400 13px {FONT};fill:{t['ink2']}}}"
        f".strong{{font-weight:600;fill:{t['ink']}}}"
        f".val{{font:500 13px {FONT};fill:{t['ink']};font-variant-numeric:tabular-nums}}"
        f".note{{font:400 12.5px {FONT};fill:{t['ink2']}}}"
    )
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{height}" viewBox="0 0 {WIDTH} {height}" '
        f'role="img" aria-label="sporos benchmark results: sporos is fastest in every panel">'
        f"<title>sporos against copMEM2, MUMmer4, minibwa, BWA, simd-minimizers and minimap2</title>"
        f"<style>{style}</style>" + "".join(parts) + "</svg>\n"
    )


def main():
    out = Path(__file__).resolve().parent.parent / "docs"
    out.mkdir(exist_ok=True)
    for theme in THEMES:
        (out / f"benchmarks-{theme}.svg").write_text(render(theme))


if __name__ == "__main__":
    main()
