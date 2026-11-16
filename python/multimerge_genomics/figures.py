"""Diagnostic figures for MMGEN results, as standalone SVG.

Every function takes plain Python data (usually a result dictionary of the
engine, or part of one) and returns an SVG string: the web interface shows it,
and a script or notebook can save it for a paper (`save(svg, "fig.svg")`).
No plotting library is needed. Text is drawn in `currentColor`, so a figure
reads on a light or a dark page; marks use a colour-blind-checked palette, and
when more than three groups share a plot they also differ by marker shape.
Each mark carries a `<title>`, which browsers show as a tooltip.

The figures are tools to judge a result, not decoration: each one is the view
in which the usual artefact shows (a second peak in a k-mer spectrum, an
inflated QQ plot, a deletion that is really missing data...).
"""
from __future__ import annotations

import html as _html
import math
import statistics

# Categorical slots (fixed order; a colour follows its entity, never its rank).
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
NEUTRAL = "#898781"
LIGHT_BLUE = "#86b6ef"
STATUS = {"good": "#0ca30c", "warning": "#fab219", "serious": "#ec835a", "critical": "#d03b3b"}
FONT = "font: 12px system-ui, -apple-system, 'Segoe UI', sans-serif"
SHAPES = ["circle", "square", "triangle", "diamond", "cross", "triangle-down", "star", "plus"]

# GRCh38 primary assembly: lengths and approximate centromere midpoints (bp).
GRCH38 = {
    "1": (248956422, 123400000), "2": (242193529, 93900000), "3": (198295559, 90900000), "4": (190214555, 50000000),
    "5": (181538259, 48800000), "6": (170805979, 59800000), "7": (159345973, 60100000), "8": (145138636, 45200000),
    "9": (138394717, 43000000), "10": (133797422, 39800000), "11": (135086622, 53400000), "12": (133275309, 35500000),
    "13": (114364328, 17700000), "14": (107043718, 17200000), "15": (101991189, 19000000), "16": (90338345, 36800000),
    "17": (83257441, 25100000), "18": (80373285, 18500000), "19": (58617616, 26200000), "20": (64444167, 28100000),
    "21": (46709983, 12000000), "22": (50818468, 15000000), "X": (156040895, 60600000), "Y": (57227415, 10400000),
}


# ------------------------------------------------------------------ basics
def esc(text) -> str:
    return _html.escape(str(text), quote=True)


def save(svg: str, path: str) -> str:
    """Writes the SVG to a file (UTF-8) and returns the path."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(svg)
    return path


def fmt(v, digits=3) -> str:
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "—"
    if isinstance(v, int):
        return f"{v:,}"
    if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e6):
        return f"{v:.2e}"
    return f"{v:.{digits}g}"


def ticks(a: float, b: float, n: int = 6) -> list:
    if not (math.isfinite(a) and math.isfinite(b)) or a == b:
        return [a]
    raw = (b - a) / n
    mag = 10 ** math.floor(math.log10(raw))
    step = next(m * mag for m in (1, 2, 2.5, 5, 10) if m * mag >= raw)
    first = math.ceil(a / step - 1e-9) * step
    out, v = [], first
    while v <= b + 1e-9 * abs(step):
        out.append(round(v, 12))
        v += step
    return out


def tick_label(v) -> str:
    if v == 0:
        return "0"
    if abs(v) >= 1e6 and v == int(v):
        return f"{v / 1e6:g}M"
    if abs(v) >= 1e4 and v == int(v):
        return f"{v / 1e3:g}k"
    return f"{v:g}"


class Svg:
    """A small SVG builder: shapes in user units, text in currentColor."""

    def __init__(self, width: int, height: int, title: str):
        self.w, self.h, self.title, self.parts = width, height, title, []

    def add(self, s: str):
        self.parts.append(s)

    def text(self, x, y, s, anchor="start", size=12, weight=None, opacity=1.0, rotate=None, cls=None):
        tr = f' transform="rotate({rotate} {x:.1f} {y:.1f})"' if rotate else ""
        wt = f' font-weight="{weight}"' if weight else ""
        self.add(f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" font-size="{size}" fill="currentColor"'
                 f' fill-opacity="{opacity}"{wt}{tr}>{esc(s)}</text>')

    def line(self, x1, y1, x2, y2, color="currentColor", width=1.0, opacity=1.0, dash=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" stroke-width="{width}"'
                 f' stroke-opacity="{opacity}"{d}/>')

    def rect(self, x, y, w, h, fill, opacity=1.0, rx=0, tip=None, stroke=None):
        st = f' stroke="{stroke}" stroke-width="1"' if stroke else ""
        t = f"<title>{esc(tip)}</title>" if tip else ""
        self.add(f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(w, 0):.1f}" height="{max(h, 0):.1f}" rx="{rx}" fill="{fill}"'
                 f' fill-opacity="{opacity}"{st}>{t}</rect>')

    def path(self, d, fill="none", stroke="currentColor", width=1.5, opacity=1.0, fill_opacity=None, dash=None, tip=None):
        fo = f' fill-opacity="{fill_opacity}"' if fill_opacity is not None else ""
        ds = f' stroke-dasharray="{dash}"' if dash else ""
        t = f"<title>{esc(tip)}</title>" if tip else ""
        self.add(f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{width}" stroke-opacity="{opacity}"{fo}{ds}>{t}</path>')

    def marker(self, x, y, shape="circle", r=4.0, fill=SERIES[0], opacity=0.85, tip=None):
        t = f"<title>{esc(tip)}</title>" if tip else ""
        ring = ' stroke="white" stroke-opacity="0.55" stroke-width="0.8"'
        if shape == "circle":
            self.add(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}" fill-opacity="{opacity}"{ring}>{t}</circle>')
            return
        pts = {
            "square": [(-r * .85, -r * .85), (r * .85, -r * .85), (r * .85, r * .85), (-r * .85, r * .85)],
            "triangle": [(0, -r * 1.1), (r, r * .8), (-r, r * .8)],
            "triangle-down": [(0, r * 1.1), (r, -r * .8), (-r, -r * .8)],
            "diamond": [(0, -r * 1.2), (r, 0), (0, r * 1.2), (-r, 0)],
            "star": [(r * math.cos(a), r * math.sin(a)) if i % 2 == 0 else (r * .45 * math.cos(a), r * .45 * math.sin(a))
                     for i, a in enumerate(-math.pi / 2 + k * math.pi / 5 for k in range(10))],
        }
        if shape in ("cross", "plus"):
            a = r if shape == "plus" else r * .8
            seg = [(-a, 0, a, 0), (0, -a, 0, a)] if shape == "plus" else [(-a, -a, a, a), (-a, a, a, -a)]
            d = " ".join(f"M{x + p:.1f} {y + q:.1f}L{x + s:.1f} {y + u:.1f}" for p, q, s, u in seg)
            self.add(f'<path d="{d}" stroke="{fill}" stroke-width="2.2" stroke-opacity="{opacity}">{t}</path>')
            return
        d = "M" + "L".join(f"{x + p:.1f} {y + q:.1f}" for p, q in pts[shape]) + "Z"
        self.add(f'<path d="{d}" fill="{fill}" fill-opacity="{opacity}"{ring}>{t}</path>')

    def render(self) -> str:
        return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.w} {self.h}" width="100%" '
                f'style="max-width: {self.w}px; {FONT}; color: inherit;" role="img" aria-label="{esc(self.title)}">'
                + "".join(self.parts) + "</svg>")


def legend(svg: Svg, items, x, y, max_width=None):
    """items: (label, colour, shape or 'bar' or 'line'). Wraps at max_width."""
    cx, cy = x, y
    limit = x + (max_width or svg.w - x - 10)
    for label, colour, shape in items:
        w = 22 + 6.6 * len(str(label))
        if cx + w > limit and cx > x:
            cx, cy = x, cy + 18
        if shape == "bar":
            svg.rect(cx, cy - 8, 12, 10, colour, rx=2)
        elif shape == "line":
            svg.line(cx, cy - 3, cx + 14, cy - 3, colour, 2)
        elif shape == "dash":
            svg.line(cx, cy - 3, cx + 14, cy - 3, colour, 1.5, dash="4 3")
        else:
            svg.marker(cx + 6, cy - 3, shape, 4.5, colour, 0.95)
        svg.text(cx + 18, cy + 1, label, size=12, opacity=0.85)
        cx += w + 12
    return cy


class Axes:
    """Linear (or log10) x/y axes inside a frame, with recessive grid lines."""

    def __init__(self, svg: Svg, left, top, width, height, x0, x1, y0, y1, xlog=False, ylog=False):
        self.svg, self.L, self.T, self.W, self.H = svg, left, top, width, height
        self.xlog, self.ylog = xlog, ylog
        self.x0, self.x1 = (math.log10(x0), math.log10(x1)) if xlog else (x0, x1)
        self.y0, self.y1 = (math.log10(y0), math.log10(y1)) if ylog else (y0, y1)
        if self.x1 == self.x0:
            self.x1 += 1
        if self.y1 == self.y0:
            self.y1 += 1

    def X(self, v):
        v = math.log10(v) if self.xlog else v
        return self.L + (v - self.x0) / (self.x1 - self.x0) * self.W

    def Y(self, v):
        v = math.log10(v) if self.ylog else v
        return self.T + self.H - (v - self.y0) / (self.y1 - self.y0) * self.H

    def inside(self, x, y):
        return self.L - 0.5 <= x <= self.L + self.W + 0.5 and self.T - 0.5 <= y <= self.T + self.H + 0.5

    def draw(self, x_label="", y_label="", xt=None, yt=None, x_fmt=tick_label, y_fmt=tick_label, grid=True):
        s = self.svg
        if yt is None:
            yt = [10 ** e for e in range(math.floor(self.y0), math.ceil(self.y1) + 1)] if self.ylog else ticks(self.y0, self.y1)
        if xt is None:
            xt = [10 ** e for e in range(math.floor(self.x0), math.ceil(self.x1) + 1)] if self.xlog else ticks(self.x0, self.x1)
        for v in yt:
            y = self.Y(v)
            if self.T - 0.5 <= y <= self.T + self.H + 0.5:
                if grid:
                    s.line(self.L, y, self.L + self.W, y, opacity=0.1)
                s.text(self.L - 6, y + 4, y_fmt(v), "end", 11, opacity=0.7)
        for v in xt:
            x = self.X(v)
            if self.L - 0.5 <= x <= self.L + self.W + 0.5:
                s.text(x, self.T + self.H + 16, x_fmt(v), "middle", 11, opacity=0.7)
        s.line(self.L, self.T + self.H, self.L + self.W, self.T + self.H, opacity=0.35)
        s.line(self.L, self.T, self.L, self.T + self.H, opacity=0.35)
        if x_label:
            s.text(self.L + self.W / 2, self.T + self.H + 34, x_label, "middle", 12, opacity=0.85)
        if y_label:
            s.text(self.L - 44, self.T + self.H / 2, y_label, "middle", 12, opacity=0.85, rotate=-90)


def _span(values, pad=0.05, floor=None):
    vals = [v for v in values if v is not None and math.isfinite(v)]
    if not vals:
        return 0.0, 1.0
    a, b = min(vals), max(vals)
    if floor is not None:
        a = min(a, floor)
    if a == b:
        a, b = a - 1, b + 1
    p = (b - a) * pad
    return a - p, b + p


def natural_key(name: str):
    n = str(name)
    if n.lower().startswith("chr"):
        n = n[3:]
    order = {"X": 1000, "Y": 1001, "M": 1002, "MT": 1002}
    if n.isdigit():
        return (0, int(n), "")
    if n.upper() in order:
        return (0, order[n.upper()], "")
    return (1, 0, n)


# ------------------------------------------------------------------ association (screen 30)
def _norm_inv(p: float) -> float:
    return statistics.NormalDist().inv_cdf(p)


def genomic_inflation(pvalues) -> float | None:
    """lambda_GC: the median of the 1-df chi-square statistics of the p-values
    over its expected value under the null (0.4549)."""
    ps = sorted(p for p in pvalues if p is not None and 0 < p <= 1)
    if len(ps) < 5:
        return None
    med = statistics.median(ps)
    z = _norm_inv(max(med / 2, 1e-300))
    return (z * z) / 0.454936423119572


def qq_plot(pvalues, title="QQ plot of the p-values", width=460, height=400) -> str:
    """Observed against expected -log10 p under the null, with the 95% band of
    the order statistics and lambda_GC. Points lifted along the whole line mean
    inflation (population structure, relatives, a miscalibrated test); a straight
    line that rises only at the end is what true signals look like."""
    ps = sorted(p for p in pvalues if p is not None and 0 < p <= 1)
    n = len(ps)
    s = Svg(width, height, title)
    s.text(10, 18, title, size=13, weight=600)
    if n == 0:
        s.text(10, 44, "No p-values to plot.", opacity=0.7)
        return s.render()
    exp = [-math.log10((i + 0.5) / n) for i in range(n)]
    obs = [-math.log10(p) for p in ps]
    top = max(max(exp), max(obs)) * 1.05
    ax = Axes(s, 60, 34, width - 80, height - 84, 0, top, 0, top)
    ax.draw("expected −log10 p", "observed −log10 p")
    band = []
    try:
        from scipy.stats import beta
        for i in range(n):
            band.append((exp[i], -math.log10(beta.ppf(0.975, i + 1, n - i)), -math.log10(beta.ppf(0.025, i + 1, n - i))))
    except Exception:
        band = []
    if band:
        up = "".join(f"{'M' if i == 0 else 'L'}{ax.X(e):.1f} {ax.Y(min(hi, top)):.1f}" for i, (e, lo, hi) in enumerate(band))
        dn = "".join(f"L{ax.X(e):.1f} {ax.Y(lo):.1f}" for e, lo, hi in reversed(band))
        s.path(up + dn + "Z", fill=NEUTRAL, stroke="none", fill_opacity=0.18)
    s.line(ax.X(0), ax.Y(0), ax.X(top), ax.Y(top), NEUTRAL, 1.2, 0.9, "4 3")
    for e, o, p in zip(exp, obs, ps):
        s.marker(ax.X(e), ax.Y(o), "circle", 3.2, SERIES[0], 0.85, tip=f"p = {fmt(p)} (expected {fmt(10 ** -e)})")
    lam = genomic_inflation(ps)
    s.text(ax.L + 10, ax.T + 14, f"λGC = {lam:.3f}" if lam is not None else "λGC: too few tests", "start", 12, 600)
    s.text(ax.L + 10, ax.T + 30, f"{n:,} tests", "start", 11, opacity=0.7)
    if band:
        legend(s, [("95% band under the null", NEUTRAL, "bar")], ax.L + 10, ax.T + 50)
    return s.render()


def manhattan(rows, title="Association along the genome", fdr_key="fdr", width=900, height=330, labels=8) -> str:
    """-log10 p of each tested allele by chromosome and position; the alleles with
    FDR < 0.05 in orange, a dashed line at the Bonferroni threshold."""
    pts = [r for r in rows if r.get("p") is not None and r["p"] > 0]
    s = Svg(width, height, title)
    s.text(10, 18, title, size=13, weight=600)
    if not pts:
        s.text(10, 44, "No p-values to plot.", opacity=0.7)
        return s.render()
    chroms = sorted({str(r["chrom"]) for r in pts}, key=natural_key)

    def bare(c):
        c = c[3:] if c.lower().startswith("chr") else c
        return "X" if c.upper() == "X" else ("Y" if c.upper() == "Y" else c)

    if all(bare(c) in GRCH38 for c in chroms):
        # the whole chromosome, as a genome-wide plot draws it (GRCh38 lengths)
        spans = {c: (1, GRCH38[bare(c)][0]) for c in chroms}
    else:
        spans = {c: (min(r["pos"] for r in pts if str(r["chrom"]) == c), max(r["pos"] for r in pts if str(r["chrom"]) == c)) for c in chroms}
    total = sum(max(b - a, 1) for a, b in spans.values())
    share = [max((spans[c][1] - spans[c][0]) / total, 0.6 / len(chroms)) for c in chroms]
    norm = sum(share)
    L, R, T, B = 60, 16, 34, 52
    pw = width - L - R
    offs, x = {}, L
    for c, sh in zip(chroms, share):
        w = pw * sh / norm
        offs[c] = (x, w)
        x += w
    ymax = max(-math.log10(r["p"]) for r in pts) * 1.1
    bonf = -math.log10(0.05 / len(pts))
    ymax = max(ymax, bonf * 1.1)
    ax = Axes(s, L, T, pw, height - T - B, 0, 1, 0, ymax)
    ax.draw("", "−log10 p", xt=[])
    for i, c in enumerate(chroms):
        x0, w = offs[c]
        if i % 2 == 1:
            s.rect(x0, T, w, height - T - B, "currentColor", 0.04)
        s.text(x0 + w / 2, height - B + 16, c.replace("chr", ""), "middle", 11, opacity=0.75)
    s.text(L + pw / 2, height - 12, "chromosome", "middle", 12, opacity=0.85)
    s.line(L, ax.Y(bonf), L + pw, ax.Y(bonf), NEUTRAL, 1.2, 0.9, "4 3")
    s.text(L + pw - 4, ax.Y(bonf) - 4, f"Bonferroni 0.05/{len(pts)}", "end", 10, opacity=0.7)
    hits = []
    for r in pts:
        c = str(r["chrom"])
        x0, w = offs[c]
        a, b = spans[c]
        px = x0 + 6 + (w - 12) * ((r["pos"] - a) / (b - a) if b > a else 0.5)
        py = ax.Y(-math.log10(r["p"]))
        sig = r.get(fdr_key) is not None and r[fdr_key] < 0.05
        col = SERIES[1] if sig else (SERIES[0] if chroms.index(c) % 2 == 0 else LIGHT_BLUE)
        s.marker(px, py, "circle", 3.6 if sig else 3, col, 0.9,
                 tip=f"{r.get('id', '')} {c}:{r['pos']:,} {r.get('ref', '')}>{r.get('alt', '')}  p = {fmt(r['p'])}  FDR = {fmt(r.get(fdr_key))}")
        hits.append((r["p"], px, py, r.get("id") or f"{c}:{r['pos']}"))
    placed = []
    for p, px, py, lab in sorted(hits)[:labels]:
        if p >= 0.05:
            break
        w = 6 * len(str(lab))
        x0 = px + 6 if px + 6 + w < width - R else px - 6 - w
        box = (x0, py - 15, x0 + w, py - 3)
        if all(box[2] < b[0] or box[0] > b[2] or box[3] < b[1] or box[1] > b[3] for b in placed):
            placed.append(box)
            s.text(x0, py - 5, lab, size=10, opacity=0.85)
    legend(s, [("FDR < 0.05", SERIES[1], "circle"), ("FDR ≥ 0.05", SERIES[0], "circle")], width - 250, 20)
    return s.render()


def forest(rows, title="Effect of each allele (95% CI)", null=0.0, log=False, x_label="effect per copy", width=640, max_rows=20) -> str:
    """rows: (label, estimate, low, high[, note]). One line per allele, the
    interval as a bar and the estimate as a dot; the dashed line is no effect."""
    rows = [r for r in rows if r[1] is not None and r[2] is not None and r[3] is not None
            and all(math.isfinite(v) and (not log or v > 0) for v in r[1:4])][:max_rows]
    h = 60 + 22 * max(len(rows), 1) + 40
    s = Svg(width, h, title)
    s.text(10, 18, title, size=13, weight=600)
    if not rows:
        s.text(10, 44, "No finite interval to show.", opacity=0.7)
        return s.render()
    lab_w = min(220, 12 + 6.6 * max(len(str(r[0])) for r in rows))
    vals = [v for r in rows for v in r[1:4]] + [null]
    if log:
        lo, hi = min(vals), max(vals)
        lo, hi = lo / 1.15, hi * 1.15
    else:
        lo, hi = _span(vals, 0.06)
    ax = Axes(s, lab_w + 10, 34, width - lab_w - 30, 22 * len(rows), lo, hi, 0, len(rows), xlog=log)
    xt = None
    if log:
        cand = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
        xt = [v for v in cand if lo <= v <= hi] or [null]
    ax.draw(x_label, "", xt=xt, yt=[], grid=False)
    s.line(ax.X(null), ax.T, ax.X(null), ax.T + ax.H, NEUTRAL, 1.2, 0.9, "4 3")
    for i, r in enumerate(rows):
        y = ax.T + 22 * i + 11
        s.text(ax.L - 8, y + 4, r[0], "end", 11, opacity=0.85)
        x1, x2 = ax.X(max(r[2], lo)), ax.X(min(r[3], hi))
        s.line(x1, y, x2, y, SERIES[0], 2)
        s.marker(ax.X(r[1]), y, "square", 4, SERIES[0], 0.95,
                 tip=f"{r[0]}: {fmt(r[1])} [{fmt(r[2])}, {fmt(r[3])}]" + (f" — {r[4]}" if len(r) > 4 and r[4] else ""))
    return s.render()


# ------------------------------------------------------------------ relatedness (screen 16)
KING_CUTS = [(0.354, "duplicate / twin"), (0.177, "1st degree"), (0.0884, "2nd degree"), (0.0442, "3rd degree")]


def king_plot(pairs, title="Kinship against IBS0 (the KING plot)", width=560, height=400) -> str:
    """pairs: (kinship, ibs0, label, relationship). Parent and child share an
    allele at every site (IBS0 near 0) with kinship near 0.25; full siblings
    have the same kinship but some IBS0. Points between the bands are the ones
    to look at again (contamination, a sample swap, or a missing relative)."""
    pairs = [p for p in pairs if p[0] is not None]
    s = Svg(width, height, title)
    s.text(10, 18, title, size=13, weight=600)
    if not pairs:
        s.text(10, 44, "No pairs to plot.", opacity=0.7)
        return s.render()
    xs = [p[1] for p in pairs if p[1] is not None] or [0.0]
    x1 = max(max(xs) * 1.1, 0.02)
    y0 = min(min(p[0] for p in pairs) - 0.03, -0.05)
    y1 = max(max(p[0] for p in pairs) + 0.03, 0.4)
    ax = Axes(s, 64, 34, width - 200, height - 84, 0, x1, y0, y1)
    ax.draw("IBS0 (sites sharing no allele, per site where both are heterozygous)", "kinship")
    for cut, name in KING_CUTS:
        if y0 < cut < y1:
            s.line(ax.L, ax.Y(cut), ax.L + ax.W, ax.Y(cut), NEUTRAL, 1, 0.8, "4 3")
            s.text(ax.L + ax.W + 6, ax.Y(cut) + 4, f"{name} ≥ {cut:g}", size=10, opacity=0.75)
    pairs = [(k, ibs0, lab, str(rel).replace("_", " ")) for k, ibs0, lab, rel in pairs]
    classes = []
    for p in pairs:
        if p[3] not in classes:
            classes.append(p[3])
    style = {c: (SERIES[i % len(SERIES)], SHAPES[i % len(SHAPES)]) for i, c in enumerate(classes)}
    for k, ibs0, lab, rel in pairs:
        col, shp = style[rel]
        x = ax.X(ibs0 if ibs0 is not None else 0)
        s.marker(x, ax.Y(k), shp, 4.5, col, 0.85, tip=f"{lab}: kinship {k:.4f}, IBS0 {fmt(ibs0)} — {rel}")
    legend(s, [(c, style[c][0], style[c][1]) for c in classes], ax.L + ax.W + 6, ax.T + ax.H - 18 * len(classes) + 6, 150)
    return s.render()


# ------------------------------------------------------------------ homozygosity (screen 20)
def karyogram(runs, title="Runs of homozygosity along the chromosomes", sizes=None, width=820) -> str:
    """runs: (chrom, start0, end). The autosomes and X, each drawn to scale, with
    the runs as orange blocks and the total per chromosome on the right. A few
    long runs on several chromosomes: recent shared ancestry of the parents; many
    short runs everywhere: an old, small population; one chromosome covered end
    to end: uniparental disomy. Chromosomes absent from the file show no runs:
    no data is not the same as no homozygosity."""
    sizes = sizes or GRCH38
    chroms = [c for c in sizes if c != "Y"]
    per = {c: [] for c in chroms}
    for c, a, b in runs:
        key = str(c)[3:] if str(c).lower().startswith("chr") else str(c)
        key = "X" if key.upper() == "X" else key
        if key in per:
            per[key].append((a, b))
    row = 19
    h = 50 + row * len(chroms) + 30
    s = Svg(width, h, title)
    s.text(10, 18, title, size=13, weight=600)
    L, R = 44, 110
    longest = max(sizes[c][0] for c in chroms)
    scale = (width - L - R) / longest
    for i, c in enumerate(chroms):
        y = 38 + i * row
        length, cen = sizes[c]
        s.text(L - 8, y + 10, c, "end", 11, opacity=0.8)
        s.rect(L, y + 3, length * scale, 9, NEUTRAL, 0.28, rx=4, tip=f"chr{c}: {length / 1e6:.1f} Mb")
        if cen:
            s.line(L + cen * scale, y + 1, L + cen * scale, y + 14, "currentColor", 1, 0.35)
        tot = 0
        for a, b in sorted(per[c]):
            tot += b - a
            s.rect(L + a * scale, y + 1, max((b - a) * scale, 1.5), 13, SERIES[1], 0.9, rx=2,
                   tip=f"chr{c}:{a + 1:,}-{b:,} ({(b - a) / 1e6:.2f} Mb)")
        if tot:
            s.text(width - R + 8, y + 11, f"{tot / 1e6:.1f} Mb ({100 * tot / length:.1f}%)", size=10, opacity=0.8)
    ticks_mb = [0, 50, 100, 150, 200, 250]
    yb = 38 + len(chroms) * row + 6
    for t in ticks_mb:
        if t * 1e6 <= longest:
            x = L + t * 1e6 * scale
            s.line(x, yb, x, yb + 4, opacity=0.4)
            s.text(x, yb + 16, f"{t} Mb", "middle", 10, opacity=0.6)
    legend(s, [("run of homozygosity", SERIES[1], "bar"), ("chromosome (GRCh38 length)", NEUTRAL, "bar")], width - 420, 20)
    if not any(per.values()):
        s.text(width - R - 20, 60, "no run found", "end", 12, 600, opacity=0.8)
    return s.render()


# ------------------------------------------------------------------ ancestry (screen 17)
def structure_plot(people, groups, title="Ancestry proportions (STRUCTURE-style)", width=820, colors=None) -> str:
    """people: (name, [proportion per group]) in the order of `groups`. One
    horizontal bar per person, sorted by the main group and its share; the
    proportions sum to one."""
    colors = colors or {g: SERIES[i % len(SERIES)] for i, g in enumerate(groups)}

    def main(p):
        j = max(range(len(groups)), key=lambda q: p[1][q])
        return (j, -p[1][j])

    order = sorted(people, key=main)
    row = 22
    lab_w = min(200, 14 + 6.6 * max((len(str(p[0])) for p in people), default=4))
    h = 78 + row * max(len(order), 1) + 30
    s = Svg(width, h, title)
    s.text(10, 18, title, size=13, weight=600)
    L, R = lab_w + 10, 20
    pw = width - L - R
    for i, (name, props) in enumerate(order):
        y = 58 + i * row
        s.text(L - 8, y + 13, name, "end", 11, opacity=0.85)
        x = L
        for g, p in zip(groups, props):
            w = pw * max(p, 0)
            if w >= 0.5:
                s.rect(x, y + 2, max(w - 2, 0.5), row - 6, colors[g], 0.95, rx=2, tip=f"{name}: {g} {100 * p:.1f}%")
                if w > 44:
                    s.add(f'<text x="{x + w / 2 - 1:.1f}" y="{y + 14:.1f}" text-anchor="middle" font-size="10" fill="white">{100 * p:.0f}%</text>')
            x += w
    yb = 58 + len(order) * row + 4
    for t in (0, 0.25, 0.5, 0.75, 1):
        s.text(L + pw * t, yb + 10, f"{int(t * 100)}%", "middle", 10, opacity=0.6)
    legend(s, [(g, colors[g], "bar") for g in groups], L, 40, pw)
    return s.render()


def group_scatter(points, x_label, y_label, title, width=560, height=430, ellipses=True, label_points=False, colors=None,
                  extra_points=None) -> str:
    """points: (x, y, group, label). Groups get a colour and, past three, a
    marker shape too; with `ellipses`, the 95% normal ellipse of each group of
    three or more points. `extra_points`: (x, y, label) drawn on top, larger,
    as the samples being examined (e.g. people projected on a reference)."""
    pts = [p for p in points if p[0] is not None and p[1] is not None and math.isfinite(p[0]) and math.isfinite(p[1])]
    extra = [p for p in (extra_points or []) if p[0] is not None and p[1] is not None]
    s = Svg(width, height, title)
    s.text(10, 18, title, size=13, weight=600)
    if not pts and not extra:
        s.text(10, 44, "Nothing to plot.", opacity=0.7)
        return s.render()
    groups = []
    for p in pts:
        if p[2] not in groups:
            groups.append(p[2])
    colors = colors or {g: SERIES[i % len(SERIES)] for i, g in enumerate(groups)}
    shapes = {g: ("circle" if len(groups) <= 3 else SHAPES[i % len(SHAPES)]) for i, g in enumerate(groups)}
    allx = [p[0] for p in pts + extra]
    ally = [p[1] for p in pts + extra]
    ell = {}
    if ellipses:
        for g in groups:
            gp = [(p[0], p[1]) for p in pts if p[2] == g]
            e = _ellipse(gp)
            if e:
                ell[g] = e
                allx += [q[0] for q in e]
                ally += [q[1] for q in e]
    x0, x1 = _span(allx)
    y0, y1 = _span(ally)
    legend_h = 18 * (1 + (sum(22 + 6.6 * len(str(g)) + 12 for g in groups) // (width - 80)))
    ax = Axes(s, 64, 34 + legend_h, width - 84, height - 84 - legend_h, x0, x1, y0, y1)
    ax.draw(x_label, y_label)
    for g, e in ell.items():
        d = "M" + "L".join(f"{ax.X(a):.1f} {ax.Y(b):.1f}" for a, b in e) + "Z"
        s.path(d, fill=colors[g], stroke=colors[g], width=1.2, opacity=0.8, fill_opacity=0.07, tip=f"{g}: 95% ellipse")
    for x, y, g, lab in pts:
        s.marker(ax.X(x), ax.Y(y), shapes[g], 3.6, colors[g], 0.8, tip=f"{lab} ({g}): {fmt(x)}, {fmt(y)}")
        if label_points:
            s.text(ax.X(x) + 6, ax.Y(y) - 5, lab, size=10, opacity=0.8)
    for x, y, lab in extra:
        s.add(f'<circle cx="{ax.X(x):.1f}" cy="{ax.Y(y):.1f}" r="7" fill="none" stroke="currentColor" stroke-width="2"><title>{esc(lab)}</title></circle>')
        s.marker(ax.X(x), ax.Y(y), "circle", 3, "currentColor", 1, tip=lab)
        s.text(ax.X(x) + 10, ax.Y(y) - 8, lab, size=11, weight=600)
    legend(s, [(g, colors[g], shapes[g]) for g in groups] + ([("examined", "currentColor", "circle")] if extra else []),
           64, 34, width - 80)
    return s.render()


def _ellipse(points, conf_chi2=5.991, n=60):
    if len(points) < 3:
        return None
    mx = sum(p[0] for p in points) / len(points)
    my = sum(p[1] for p in points) / len(points)
    sxx = sum((p[0] - mx) ** 2 for p in points) / (len(points) - 1)
    syy = sum((p[1] - my) ** 2 for p in points) / (len(points) - 1)
    sxy = sum((p[0] - mx) * (p[1] - my) for p in points) / (len(points) - 1)
    tr, det = sxx + syy, sxx * syy - sxy * sxy
    disc = math.sqrt(max(tr * tr / 4 - det, 0))
    l1, l2 = tr / 2 + disc, tr / 2 - disc
    if l1 <= 0:
        return None
    ang = 0.5 * math.atan2(2 * sxy, sxx - syy)
    a, b = math.sqrt(conf_chi2 * max(l1, 0)), math.sqrt(conf_chi2 * max(l2, 0))
    ca, sa = math.cos(ang), math.sin(ang)
    return [(mx + a * math.cos(t) * ca - b * math.sin(t) * sa, my + a * math.cos(t) * sa + b * math.sin(t) * ca)
            for t in (2 * math.pi * i / n for i in range(n))]


# ------------------------------------------------------------------ k-mer spectrum (screens 1 and 13)
def spectrum_model(histogram, k=None) -> dict:
    """Reads a k-mer abundance spectrum as GenomeScope does by eye: the error
    k-mers (counts below the first valley), the coverage peak, a second peak at
    half of it (a heterozygous diploid genome), and from them the genome size
    (solid k-mers / k-mer coverage of one copy), the share of repeats and a
    rough heterozygosity. Rough: GenomeScope fits a model; this reads the peaks."""
    h = sorted((int(c), int(n)) for c, n in histogram if int(n) > 0 and int(c) > 0)
    out = {"valley": None, "peak": None, "het_peak": None, "haploid_coverage": None, "genome_size": None,
           "error_fraction": None, "repeat_fraction": None, "heterozygosity": None, "notes": []}
    if len(h) < 3:
        out["notes"].append("too few points in the spectrum")
        return out
    top = h[-1][0]
    dense = [0] * (min(top, 10001) + 2)
    for c, n in h:
        if c < len(dense):
            dense[c] += n
    sm = [0.0] * len(dense)
    for i in range(1, len(dense) - 1):
        sm[i] = (dense[i - 1] + 2 * dense[i] + dense[i + 1]) / 4
    valley = None
    for i in range(2, len(sm) - 1):
        if sm[i] <= sm[i - 1] and sm[i] < sm[i + 1]:
            valley = i
            break
    if valley is None:
        out["notes"].append("no valley after the error k-mers: the coverage may be too low to separate them")
        valley = 1
    peak = max(range(valley, len(sm) - 1), key=lambda i: sm[i] * (1 if i > 1 else 0)) if valley < len(sm) - 1 else None
    out["valley"], out["peak"] = valley, peak
    tot = sum(c * n for c, n in h)
    if tot == 0 or peak is None:
        return out
    err = sum(c * n for c, n in h if c < valley)
    out["error_fraction"] = err / tot
    # a heterozygous peak at about half the main peak?
    half = peak / 2
    cand = [i for i in range(max(valley, int(half * 0.7)), min(len(sm) - 1, int(half * 1.3) + 1))
            if sm[i] >= sm[i - 1] and sm[i] >= sm[i + 1] and i != peak]
    het = max(cand, key=lambda i: sm[i]) if cand else None
    if het is not None and sm[het] >= 0.1 * sm[peak]:
        out["het_peak"] = het
        cov = het                      # coverage of one copy (one haplotype)
        hom = peak
    else:
        # or the main peak is the heterozygous one and a smaller one sits at twice it
        dbl = [i for i in range(int(peak * 1.7), min(len(sm) - 1, int(peak * 2.3) + 1)) if sm[i] >= sm[i - 1] and sm[i] >= sm[i + 1]]
        dpk = max(dbl, key=lambda i: sm[i]) if dbl else None
        if dpk is not None and sm[dpk] >= 0.1 * sm[peak]:
            out["het_peak"] = peak
            cov, hom = peak, dpk
        else:
            cov, hom = peak, None
    solid = sum(c * n for c, n in h if c >= valley)
    if out["het_peak"] is not None:
        out["haploid_coverage"] = cov
        size = solid / (2 * cov)       # diploid: the genome of one haplotype
        out["genome_size"] = size
        het_kmers = sum(n for c, n in h if valley <= c <= 1.5 * cov)
        if k:
            out["heterozygosity"] = min(1.0, het_kmers / (2 * k * size)) if size else None
        out["notes"].append(f"two peaks ({cov}x and {hom}x): a diploid genome with heterozygous sites; the genome size is one haplotype's")
    else:
        out["haploid_coverage"] = cov
        out["genome_size"] = solid / cov
        out["notes"].append(f"one peak at {cov}x: a haploid genome, or a diploid one with few heterozygous sites")
    rep = sum(c * n for c, n in h if c > 2.5 * (hom or cov))
    out["repeat_fraction"] = rep / solid if solid else None
    return out


def kmer_spectrum(histogram, title="K-mer spectrum", log_y=True, model=None, width=640, height=360, x_max=None) -> str:
    """Distinct k-mers against the times each was seen. Linear y shows the
    coverage peaks; log y also shows the error k-mers at the left and the
    repeats on the right. `model` (from spectrum_model) marks the valley, the
    peaks and the estimates."""
    pts = sorted((int(c), int(n)) for c, n in histogram if int(n) > 0 and int(c) > 0)
    s = Svg(width, height, title)
    s.text(10, 18, title, size=13, weight=600)
    if not pts:
        s.text(10, 44, "Empty spectrum.", opacity=0.7)
        return s.render()
    model = model if model is not None else spectrum_model(pts)
    pk = model.get("peak") or pts[-1][0]
    xm = x_max or min(pts[-1][0], max(3 * pk, 20))
    shown = [(c, n) for c, n in pts if c <= xm]
    if not log_y and model.get("valley"):
        ymax = max((n for c, n in shown if c >= model["valley"]), default=max(n for _, n in shown)) * 1.15
    else:
        ymax = max(n for _, n in shown) * (2 if log_y else 1.1)
    ax = Axes(s, 70, 34, width - 90, height - 84, 0.5, xm + 0.5, 0.8 if log_y else 0, ymax, ylog=log_y)
    ax.draw("times a k-mer was seen", "distinct k-mers" + (" (log)" if log_y else ""))
    bw = max(1.0, ax.W / (xm + 1) - 1)
    for c, n in shown:
        y = ax.Y(max(n, 0.8) if log_y else n)
        y = max(y, ax.T)
        s.rect(ax.X(c) - bw / 2, y, bw, ax.T + ax.H - y, SERIES[0], 0.85, tip=f"seen {c}×: {n:,} distinct k-mers")
    for key, name, col in (("valley", "valley", NEUTRAL), ("het_peak", "½ peak", SERIES[1]), ("peak", "peak", SERIES[1])):
        v = model.get(key)
        if v and v <= xm and not (key == "het_peak" and v == model.get("peak")):
            s.line(ax.X(v), ax.T, ax.X(v), ax.T + ax.H, col, 1.2, 0.9, "4 3")
            s.text(ax.X(v) + 4, ax.T + 12, f"{name} {v}×", size=10, opacity=0.85)
    lines = []
    if model.get("genome_size"):
        lines.append(f"genome size ≈ {fmt_bp(model['genome_size'])}")
    if model.get("error_fraction") is not None:
        lines.append(f"error k-mers ≈ {100 * model['error_fraction']:.1f}% of all")
    if model.get("repeat_fraction") is not None:
        lines.append(f"repeats ≈ {100 * model['repeat_fraction']:.1f}%")
    if model.get("heterozygosity") is not None:
        lines.append(f"heterozygosity ≈ {100 * model['heterozygosity']:.2f}% (rough)")
    for i, t in enumerate(lines):
        s.text(width - 22, 52 + 15 * i, t, "end", 11, opacity=0.85)
    return s.render()


def fmt_bp(n: float) -> str:
    if n >= 1e9:
        return f"{n / 1e9:.2f} Gb"
    if n >= 1e6:
        return f"{n / 1e6:.2f} Mb"
    if n >= 1e3:
        return f"{n / 1e3:.1f} kb"
    return f"{n:.0f} bp"


# ------------------------------------------------------------------ viral genomes (screen 26)
def viral_genomes(result, title="Mutations along the genome", width=900, max_genomes=40) -> str:
    """One row per genome, aligned to the reference, as Nextclade draws it: the
    genes on top; substitutions (blue ticks), deletions (red), insertions
    (violet triangles) and missing or ambiguous bases (dark grey blocks).
    A long grey block is a stretch the sequencing did not read (an amplicon
    dropout), not a deletion."""
    refs = result.get("reference_sequences") or []
    length = refs[0]["length"] if refs else max((g.get("alignment_end") or 0 for g in result.get("genomes", [])), default=1)
    genomes = [g for g in result.get("genomes", []) if g.get("aligned")][:max_genomes]
    genes = result.get("genes") or []
    lab_w = min(260, 14 + 6.4 * max((len(g["name"]) for g in genomes), default=6))
    row = 15
    gene_h = 34 if genes else 0
    h = 54 + gene_h + row * max(len(genomes), 1) + 40
    s = Svg(width, h, title)
    s.text(10, 18, title, size=13, weight=600)
    L, R = lab_w + 8, 16
    pw = width - L - R
    X = lambda p: L + (p - 1) / max(length - 1, 1) * pw
    y = 52
    if genes:
        for i, g in enumerate(genes):
            for a, b in g.get("parts", []):
                s.rect(X(a), y + (0 if i % 2 == 0 else 12), max(X(b) - X(a), 1), 10, NEUTRAL, 0.45, rx=2,
                       tip=f"{g['name']} {a:,}-{b:,} ({g.get('product', '')})")
                if X(b) - X(a) > 7 * len(g["name"]) + 4:
                    s.text((X(a) + X(b)) / 2, y + (0 if i % 2 == 0 else 12) + 8.5, g["name"], "middle", 9, opacity=0.9)
        y += gene_h
    for gi, g in enumerate(genomes):
        yy = y + gi * row
        s.text(L - 6, yy + 9, g["name"], "end", 10, opacity=0.85)
        a0, a1 = g.get("alignment_start") or 1, g.get("alignment_end") or length
        s.line(X(a0), yy + 6, X(a1), yy + 6, NEUTRAL, 1, 0.6)
        for a, b in (g.get("missing") or []) + (g.get("not_covered") or []):
            s.rect(X(a), yy + 1, max(X(b) - X(a), 1.5), 11, "#3d3d3a", 0.75, tip=f"{g['name']}: not read {a:,}-{b:,}")
        for d in g.get("deletions") or []:
            s.rect(X(d["start"]), yy + 1, max(X(d["end"]) - X(d["start"]), 2), 11, SERIES[7], 0.95,
                   tip=f"{g['name']}: deletion {d['label']}" + (f" ({', '.join(d['aa'])})" if d.get("aa") else ""))
        for m in g.get("substitutions") or []:
            s.rect(X(m["pos"]) - 0.75, yy + 1, 1.5, 11, SERIES[0], 0.95,
                   tip=f"{g['name']}: {m['label']}" + (f" ({', '.join(m['aa'])})" if m.get("aa") else ""))
        for m in g.get("insertions") or []:
            x = X(m.get("pos", 0))
            s.path(f"M{x - 3.5:.1f} {yy:.1f}L{x + 3.5:.1f} {yy:.1f}L{x:.1f} {yy + 6:.1f}Z", fill=SERIES[6], stroke="none",
                   tip=f"{g['name']}: insertion {m.get('label', '')}")
    yb = y + row * len(genomes) + 6
    for t in ticks(1, length, 8):
        if 1 <= t <= length:
            s.text(X(t), yb + 12, f"{int(t):,}", "middle", 10, opacity=0.6)
    legend(s, [("substitution", SERIES[0], "bar"), ("deletion", SERIES[7], "bar"), ("insertion", SERIES[6], "bar"),
               ("not read (N, gap)", "#3d3d3a", "bar")], width - 470, 32)
    return s.render()


# ------------------------------------------------------------------ taxonomy (screens 24 and 27)
def sunburst(lineages, title="Taxonomic composition", width=560, ranks=None, min_share=0.004) -> str:
    """lineages: (path, count) with path from the top rank down (e.g. phylum,
    class, order, family, genus). Rings from the centre outwards; each arc's
    angle is its share of the reads, coloured by the top-level group."""
    total = sum(c for _, c in lineages if c > 0)
    s = Svg(width, width + 40, title)
    s.text(10, 18, title, size=13, weight=600)
    if not total:
        s.text(10, 44, "Nothing classified.", opacity=0.7)
        return s.render()
    tree = {}
    for path, c in lineages:
        if c <= 0:
            continue
        node = tree
        for name in path:
            name = name or "unclassified"
            node = node.setdefault(name, {"_n": 0, "_c": {}})
            node["_n"] += c
            node = node["_c"]
    depth = max((len(p) for p, _ in lineages), default=1)
    cx, cy = width / 2, width / 2 + 34
    r0, ring = width * 0.07, (width * 0.44 - width * 0.07) / max(depth, 1)
    tops = sorted(tree.items(), key=lambda kv: -kv[1]["_n"])
    colors = {}
    for i, (name, _) in enumerate(tops):
        colors[name] = SERIES[i] if i < 7 else NEUTRAL

    def arc(a0, a1, ra, rb):
        if a1 - a0 >= 2 * math.pi - 1e-9:
            a1 = a0 + 2 * math.pi - 1e-4
        large = 1 if a1 - a0 > math.pi else 0
        p = lambda r, a: (cx + r * math.sin(a), cy - r * math.cos(a))
        x0, y0 = p(rb, a0)
        x1, y1 = p(rb, a1)
        x2, y2 = p(ra, a1)
        x3, y3 = p(ra, a0)
        return (f"M{x0:.1f} {y0:.1f}A{rb:.1f} {rb:.1f} 0 {large} 1 {x1:.1f} {y1:.1f}L{x2:.1f} {y2:.1f}"
                f"A{ra:.1f} {ra:.1f} 0 {large} 0 {x3:.1f} {y3:.1f}Z")

    def walk(children, a0, a1, level, top, path):
        tot = sum(v["_n"] for v in children.values())
        a = a0
        for name, v in sorted(children.items(), key=lambda kv: -kv[1]["_n"]):
            span = (a1 - a0) * v["_n"] / tot if tot else 0
            if v["_n"] / total < min_share:
                a += span
                continue
            t = top or name
            ra, rb = r0 + level * ring, r0 + (level + 1) * ring - 1.5
            share = 100 * v["_n"] / total
            s.path(arc(a, a + span, ra, rb), fill=colors.get(t, NEUTRAL), stroke="white", width=1, opacity=0.6,
                   fill_opacity=max(0.95 - 0.13 * level, 0.35), tip=" › ".join(path + [name]) + f": {share:.1f}%")
            mid = a + span / 2
            rr = (ra + rb) / 2
            lab = name if len(name) <= 16 else name[:15] + "…"
            # labels run along the arc (tangential), turned so they never read upside down
            if span * rr > 6.0 * len(lab) + 8 and rb - ra >= 13:
                tx, ty = cx + rr * math.sin(mid), cy - rr * math.cos(mid)
                rot = math.degrees(mid)
                if math.pi / 2 < mid < 3 * math.pi / 2:
                    rot -= 180
                s.add(f'<text x="{tx:.1f}" y="{ty:.1f}" text-anchor="middle" dominant-baseline="central" font-size="10" '
                      f'fill="white" transform="rotate({rot:.1f} {tx:.1f} {ty:.1f})">{esc(lab)}</text>')
            walk(v["_c"], a, a + span, level + 1, t, path + [name])
            a += span

    walk(tree, 0.0, 2 * math.pi, 0, None, [])
    s.text(cx, cy + 4, f"{int(total):,}" if float(total).is_integer() else f"{total:,.1f}", "middle", 12, 600)
    if ranks:
        s.text(10, width + 34, "rings, inside to outside: " + ", ".join(ranks[:depth]), size=10, opacity=0.7)
    legend(s, [(n, colors[n], "bar") for n, _ in tops[:7]], 10, 36, width - 20)
    return s.render()


# ------------------------------------------------------------------ pharmacogenomics (screen 19)
def pgx_evidence(gene_result, gene_table, title=None, width=760) -> str:
    """Why a pharmacogene call has alternatives: rows are the positions that
    define the candidate star alleles (the call and its alternatives), columns
    the candidates, then the person's genotype. A filled cell is a base the
    allele requires; the person's column shows what was read: reference, one
    copy, two copies, or not read (grey). Without phase, two variants of one
    person can sit on the same chromosome or on opposite ones: that is where
    the alternatives come from."""
    g = gene_result
    title = title or f"{g.get('gene', '')}: the positions behind {g.get('diplotype') or 'the call'} and its alternatives"
    names = []
    for d in [g.get("diplotype") or ""] + [a.get("diplotype") or "" for a in g.get("x_alternatives") or []]:
        for part in d.replace(" OR ", "/").split("/"):
            part = part.strip()
            if part and part not in names:
                names.append(part)
    alleles = {a["name"]: a for a in gene_table.get("alleles", [])}
    names = [n for n in names if n in alleles][:10]
    pos_by_id = {p["id"]: p for p in gene_table.get("positions", [])}
    used = []
    for n in names:
        for v in alleles[n].get("variants", []):
            if v["pos_id"] not in used:
                used.append(v["pos_id"])
    used.sort(key=lambda i: pos_by_id.get(i, {}).get("pos", 0))
    read = {}
    for e in g.get("evidence") or []:
        read[(str(e.get("chrom")), e.get("pos"))] = e
    row = 20
    col = 70
    lab_w = 190
    h = 70 + row * max(len(used), 1) + 40
    w = max(width, int(lab_w + col * (len(names) + 1) + 20))
    s = Svg(w, h, title)
    s.text(10, 18, title, size=13, weight=600)
    if not names or not used:
        s.text(10, 44, "The candidates differ by no defined position (or are not in the tables).", opacity=0.7)
        return s.render()
    chosen = set((g.get("diplotype") or "").replace(" OR ", "/").split("/"))
    for j, n in enumerate(names + ["the person"]):
        x = lab_w + j * col + col / 2
        s.text(x, 50, n, "middle", 11, 700 if n in chosen or n == "the person" else None, opacity=0.95)
    s.line(lab_w + len(names) * col + 2, 38, lab_w + len(names) * col + 2, 58 + row * len(used), opacity=0.35)
    for i, pid in enumerate(used):
        p = pos_by_id.get(pid, {})
        y = 58 + i * row
        lab = f"{p.get('rsid') or ''} {p.get('chrom', '')}:{p.get('pos', 0):,}".strip()
        s.text(lab_w - 8, y + 14, lab, "end", 10, opacity=0.85)
        ref = p.get("ref_cpic") or p.get("grch38_ref") or ""
        for j, n in enumerate(names):
            need = next((v.get("accepted") for v in alleles[n].get("variants", []) if v["pos_id"] == pid), None)
            x = lab_w + j * col
            if need:
                txt = "/".join(need)
                s.rect(x + 4, y + 2, col - 8, row - 4, SERIES[0], 0.85, rx=3, tip=f"{n} requires {txt} at {lab} (reference {ref})")
                s.add(f'<text x="{x + col / 2:.1f}" y="{y + 14:.1f}" text-anchor="middle" font-size="10" fill="white">{esc(txt[:9])}</text>')
            else:
                s.rect(x + 4, y + 2, col - 8, row - 4, NEUTRAL, 0.12, rx=3, tip=f"{n}: reference ({ref}) at {lab}")
        e = read.get((str(p.get("chrom")), p.get("pos")))
        x = lab_w + len(names) * col
        if e is None or not e.get("genotype") or "." in str(e.get("genotype")):
            s.rect(x + 6, y + 2, col - 8, row - 4, NEUTRAL, 0.45, rx=3, tip=f"{lab}: not read")
            s.text(x + col / 2 + 2, y + 14, "not read", "middle", 9, opacity=0.9)
        else:
            gt = str(e["genotype"])
            al = gt.replace("|", "/").split("/")
            n_alt = sum(1 for a in al if a != ref)
            colr = [NEUTRAL, SERIES[1], SERIES[1]][min(n_alt, 2)]
            op = [0.15, 0.45, 0.9][min(n_alt, 2)]
            s.rect(x + 6, y + 2, col - 8, row - 4, colr, op, rx=3,
                   tip=f"{lab}: {gt} ({['reference', 'one copy', 'two copies'][min(n_alt, 2)]}; {e.get('origin', '')})")
            s.text(x + col / 2 + 2, y + 14, gt, "middle", 10, 600 if n_alt else None)
    legend(s, [("base the allele requires", SERIES[0], "bar"), ("one copy read", "#f3a37f", "bar"), ("two copies read", SERIES[1], "bar"),
               ("not read", NEUTRAL, "bar")], 10, h - 12, w - 20)
    return s.render()


# ------------------------------------------------------------------ assembly graph (screens 3 and 7)
GRAPH_KINDS = {"": ("unitig", SERIES[0]), "tip": ("tip (dead end: often an error)", SERIES[7]),
               "bubble": ("bubble branch (SNP or error)", SERIES[3]), "repeat": ("repeat (two ways in, two out)", SERIES[6]),
               "isolated": ("isolated (a whole sequence)", SERIES[2])}


def assembly_graph(graph, title=None, width=900, max_height=720, label_nodes=40) -> str:
    """Draws the result of `unitig_graph` or `overlap_graph` as Bandage does:
    each node a thick line through its laid-out points (longer for longer
    sequences), an edge from the end of one to the start of the next. Colours
    say what the shape means — tips (dead ends, usually sequencing errors),
    bubbles (two parallel paths: a SNP or an error), repeats (a sequence
    entered from two places and left to two: where an assembly breaks) —
    and a closed ring is a circular molecule (a plasmid, a phage, a small
    chromosome). Overlap graphs colour each contig by its scaffold, and draw
    the joins the scaffolder did not make in grey."""
    nodes = graph.get("nodes") or []
    lay = graph.get("layout") or {}
    u = lay.get("unit", 10) or 10
    W0, H0 = (lay.get("width") or 1) + 4 * u, (lay.get("height") or 1) + 4 * u
    scale = min((width - 20) / W0, (max_height - 70) / H0)
    width_px = width
    height_px = int(H0 * scale + 70)
    s = Svg(width_px, height_px, title or "Assembly graph")
    s.text(10, 18, title or "Assembly graph", size=13, weight=600)
    if not nodes:
        s.text(10, 44, "No node to draw.", opacity=0.7)
        return s.render()
    ox, oy = 10 + 2 * u * scale, 44 + 2 * u * scale
    P = lambda p: (ox + p[0] * scale, oy + p[1] * scale)
    byid = {n["id"]: n for n in nodes}
    overlap = "scaffold" in nodes[0]
    lw = max(1.6, min(7.0, 0.45 * u * scale))
    for e in graph.get("edges") or []:
        a, b = byid.get(e["from"]), byid.get(e["to"])
        if not a or not b:
            continue
        (x1, y1), (x2, y2) = P(a["path"][-1]), P(b["path"][0])
        grey = overlap and not e.get("used", True)
        s.line(x1, y1, x2, y2, NEUTRAL if grey else "currentColor", max(0.8, lw * 0.25), 0.55 if grey else 0.7, "3 2" if grey else None)
    for n in nodes:
        pts = [P(p) for p in n.get("path") or [[n["x"], n["y"]]]]
        if overlap:
            col = SERIES[n.get("scaffold", 0) % len(SERIES)]
            what = f"contig {n['id']}, scaffold {n.get('scaffold')}"
        else:
            name, col = GRAPH_KINDS.get(n.get("kind", ""), GRAPH_KINDS[""])
            what = f"unitig {n['id']} ({name})"
        tip = f"{what}: {n.get('length', 0):,} bases, GC {100 * (n.get('gc') or 0):.0f}%" + (", circular" if n.get("circular") else "")
        if len(pts) == 1:
            s.marker(pts[0][0], pts[0][1], "circle", lw / 2 + 1, col, 0.9, tip=tip)
        else:
            d = "M" + "L".join(f"{x:.1f} {y:.1f}" for x, y in pts)
            if n.get("circular") and len(pts) > 2:
                d += "Z"
            s.add(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="{lw:.1f}" stroke-linecap="round" stroke-linejoin="round" '
                  f'stroke-opacity="0.9"><title>{esc(tip)}</title></path>')
            ex, ey = pts[-1]
            s.add(f'<circle cx="{ex:.1f}" cy="{ey:.1f}" r="{max(1.0, lw * 0.32):.1f}" fill="currentColor" fill-opacity="0.75"/>')
    if len(nodes) <= label_nodes:
        for n in nodes:
            x, y = P([n["x"], n["y"]])
            s.text(x, y - lw - 2, str(n["id"]), "middle", 10, opacity=0.8)
    if overlap:
        items = [("contig (colour: its scaffold)", SERIES[0], "bar"), ("join made", "currentColor", "line"), ("overlap not joined", NEUTRAL, "dash")]
    else:
        present = {n.get("kind", "") for n in nodes}
        items = [(GRAPH_KINDS[k][0], GRAPH_KINDS[k][1], "bar") for k in ("", "isolated", "tip", "bubble", "repeat") if k in present]
    legend(s, items, 10, 34, width_px - 20)
    return s.render()


# ------------------------------------------------------------------ genome browser (screens 9 and 21)
# The nucleotide colours genome browsers use (IGV): A green, C blue, G orange, T red.
BASE_COLOURS = {"A": "#0ca30c", "C": "#2a78d6", "G": "#d48a00", "T": "#d03b3b", "N": NEUTRAL}
FORWARD, REVERSE = "#b8d3f3", "#f4c6b3"


def genome_browser(view, features=None, title=None, width=920, max_rows=60) -> str:
    """Draws the result of `region_reads` as a genome browser does: the
    reference, the annotation (`features`: (name, start, end, strand, kind),
    kind "exon"/"CDS" drawn thick and "gene" thin), the coverage with the
    share of each base where it differs from the reference, and the reads in
    rows — blue forward, orange reverse; mismatches in the colour of their
    base, faint when the base quality is low; soft clips faded; deletions as
    a line, introns as a thin line, insertions as a violet mark. A variant
    seen only on one strand, only near the ends of the reads, only in clipped
    or low-quality bases, or only in reads with mapping quality 0 is probably
    an artefact."""
    reg = view["region"]
    a, b = reg["start"], reg["end"]
    n = b - a + 1
    title = title or f"{reg['chrom']}:{a:,}-{b:,}"
    reads = view.get("reads") or []
    rows = min(view.get("rows") or 0, max_rows)
    L, R = 70, 14
    pw = width - L - R
    bp = pw / n
    X = lambda p: L + (p - a) * bp
    ref = view.get("reference") or ""
    feats = [f for f in (features or []) if f[2] >= a and f[1] <= b]
    y_ruler, y_ref = 34, 52
    y_ann = y_ref + (16 if ref else 0)
    ann_h = 26 if feats else 0
    y_cov = y_ann + ann_h + 6
    cov_h = 64
    y_reads = y_cov + cov_h + 14
    row_h = 9 if rows <= 30 else 6
    h = y_reads + row_h * max(rows, 1) + 40
    s = Svg(width, h, title)
    s.text(10, 18, title, size=13, weight=600)
    info = f"{view.get('reads_total', len(reads)):,} reads, mean depth {fmt(view.get('mean_depth'))}"
    if view.get("downsampled"):
        info += f" — {view.get('reads_shown', len(reads)):,} drawn (a fixed sample; the coverage counts all)"
    s.text(width - R, 18, info, "end", 11, opacity=0.75)
    # ruler
    for t in ticks(a, b, 8):
        if a <= t <= b:
            s.line(X(t), y_ruler, X(t), y_ruler + 5, opacity=0.5)
            s.text(X(t), y_ruler - 3, f"{int(t):,}", "middle", 10, opacity=0.7)
    s.line(L, y_ruler + 5, L + pw, y_ruler + 5, opacity=0.3)
    focus = reg.get("focus")
    if focus:
        s.line(X(focus) + bp / 2, y_ruler + 5, X(focus) + bp / 2, h - 30, "currentColor", 1, 0.45, "3 3")
    # reference
    if ref:
        s.text(L - 6, y_ref + 9, "reference", "end", 10, opacity=0.7)
        for i, c in enumerate(ref[:n]):
            col = BASE_COLOURS.get(c.upper(), NEUTRAL)
            if bp >= 7:
                s.add(f'<text x="{X(a + i) + bp / 2:.1f}" y="{y_ref + 10:.1f}" text-anchor="middle" font-size="{min(11, bp * 0.9):.0f}" '
                      f'font-family="ui-monospace, monospace" fill="{col}">{esc(c)}</text>')
            else:
                s.rect(X(a + i), y_ref + 2, max(bp, 0.6), 9, col, 0.8)
    # annotation
    if feats:
        s.text(L - 6, y_ann + 14, "genes", "end", 10, opacity=0.7)
        for name, fs, fe, strand, kind in feats:
            x0, x1 = X(max(fs, a)), X(min(fe, b) + 1)
            if kind in ("exon", "CDS"):
                s.rect(x0, y_ann + 4, max(x1 - x0, 1), 14, SERIES[6], 0.85 if kind == "CDS" else 0.55, rx=2,
                       tip=f"{name} {kind} {fs:,}-{fe:,} ({strand})")
            else:
                s.line(x0, y_ann + 11, x1, y_ann + 11, SERIES[6], 1.5, 0.9)
                step = 40
                arrow = "›" if strand == "+" else ("‹" if strand == "-" else "")
                if arrow:
                    xx = x0 + 12
                    while xx < x1 - 6:
                        s.text(xx, y_ann + 14.5, arrow, "middle", 11, opacity=0.8)
                        xx += step
                if x1 - x0 > 7 * len(name) + 10:
                    s.text((x0 + x1) / 2, y_ann + 26, name, "middle", 10, opacity=0.85)
    # coverage
    cov = view.get("coverage") or []
    covf = {c[0]: c for c in (view.get("coverage_forward") or [])}
    maxd = max((c[1] for c in cov), default=1) or 1
    s.text(L - 6, y_cov + 12, "coverage", "end", 10, opacity=0.7)
    s.text(L - 6, y_cov + 24, f"0–{maxd}", "end", 9, opacity=0.6)
    for c in cov:
        pos, depth, A, C, G, T, dele, ins = c[:8]
        if depth <= 0:
            continue
        hgt = cov_h * depth / maxd
        x = X(pos)
        refb = ref[pos - a].upper() if ref and 0 <= pos - a < len(ref) else None
        counts = {"A": A, "C": C, "G": G, "T": T}
        alt = {k: v for k, v in counts.items() if refb and k != refb and v > 0}
        alt_total = sum(alt.values()) + dele
        tip = f"{reg['chrom']}:{pos:,}  depth {depth}  A {A}  C {C}  G {G}  T {T}  del {dele}  ins {ins}"
        if pos in covf:
            fA, fC, fG, fT = covf[pos][1:5]
            tip += f"  (forward: A {fA} C {fC} G {fG} T {fT})"
        if refb and alt_total / depth >= 0.2:
            # IGV colours the bar where at least 20% of the bases differ from the reference
            y = y_cov + cov_h
            for base in ("A", "C", "G", "T"):
                v = counts[base]
                if v:
                    hh = cov_h * v / maxd
                    s.rect(x, y - hh, max(bp - 0.3, 0.6), hh, BASE_COLOURS[base], 0.9, tip=tip)
                    y -= hh
        else:
            s.rect(x, y_cov + cov_h - hgt, max(bp - 0.3, 0.6), hgt, NEUTRAL, 0.55, tip=tip)
    s.line(L, y_cov + cov_h, L + pw, y_cov + cov_h, opacity=0.3)
    # reads
    s.text(L - 6, y_reads + 8, "reads", "end", 10, opacity=0.7)
    for r in reads:
        row = r.get("row", 0)
        if row >= max_rows:
            continue
        y = y_reads + row * row_h
        col = FORWARD if r.get("strand") == "+" else REVERSE
        op = 0.45 if (r.get("mapq") or 0) == 0 else 1.0
        tip = (f"{r.get('name')}  {r.get('strand')}  MAPQ {r.get('mapq')}  {r.get('cigar', '')}  "
               f"{r.get('start'):,}-{r.get('end'):,}  mean quality {fmt(r.get('mean_qual'))}")
        blocks = r.get("blocks") or [[r["start"], r["end"]]]
        for i, (bs, be) in enumerate(blocks):
            if be < a or bs > b:
                continue
            x0, x1 = X(max(bs, a)), X(min(be, b) + 1)
            s.rect(x0, y, max(x1 - x0, 0.8), row_h - 2, col, op, tip=tip)
        for (ds, dl) in r.get("deletions") or []:
            s.line(X(ds), y + (row_h - 2) / 2, X(ds + dl), y + (row_h - 2) / 2, "currentColor", 1.2, 0.85)
        for (ks, kl) in r.get("skips") or []:
            s.line(X(max(ks, a)), y + (row_h - 2) / 2, X(min(ks + kl, b + 1)), y + (row_h - 2) / 2, NEUTRAL, 0.8, 0.8)
        cl, cr = r.get("clip_left") or 0, r.get("clip_right") or 0
        if cl and r["start"] - cl <= b:
            s.rect(X(max(r["start"] - cl, a)), y, max(X(r["start"]) - X(max(r["start"] - cl, a)), 0.5), row_h - 2, col, 0.25,
                   tip=f"{r.get('name')}: {cl} soft-clipped bases")
        if cr and r["end"] + 1 >= a:
            s.rect(X(r["end"] + 1), y, max(X(min(r["end"] + 1 + cr, b + 1)) - X(r["end"] + 1), 0.5), row_h - 2, col, 0.25,
                   tip=f"{r.get('name')}: {cr} soft-clipped bases")
        for pos, base, q in r.get("mismatches") or []:
            if a <= pos <= b:
                faint = q is not None and q < 20
                s.rect(X(pos), y, max(bp, 1.2), row_h - 2, BASE_COLOURS.get(str(base).upper(), NEUTRAL), 0.3 if faint else 1.0,
                       tip=f"{r.get('name')}: {base} at {pos:,}, base quality {q}")
        for pos, ln, seq in r.get("insertions") or []:
            if a <= pos <= b:
                s.rect(X(pos + 1) - 1, y - 1, 2, row_h, SERIES[6], 1.0, tip=f"{r.get('name')}: insertion of {ln} after {pos:,} ({seq})")
    if (view.get("rows") or 0) > max_rows:
        s.text(L, h - 24, f"{view['rows'] - max_rows} more row(s) of reads not drawn", size=10, opacity=0.7)
    legend(s, [("forward read", FORWARD, "bar"), ("reverse read", REVERSE, "bar"), ("A", BASE_COLOURS["A"], "bar"),
               ("C", BASE_COLOURS["C"], "bar"), ("G", BASE_COLOURS["G"], "bar"), ("T", BASE_COLOURS["T"], "bar"),
               ("insertion", SERIES[6], "bar")], L, h - 8, pw)
    return s.render()


def strand_and_position(view, position, alt=None) -> dict:
    """The checks an analyst makes on a candidate variant at `position`: how
    many forward and reverse reads carry the ALT base and the reference
    (Fisher's exact test and the strand odds ratio, as GATK's SOR), and where
    in the reads the ALT base sits (a variant only within the last bases of
    the reads is usually an artefact). From the coverage of every read and
    the reads drawn."""
    cov = {c[0]: c for c in view.get("coverage") or []}
    fwd = {c[0]: c for c in view.get("coverage_forward") or []}
    ref = view.get("reference") or ""
    a = view["region"]["start"]
    if position not in cov:
        return {}
    idx = {"A": 2, "C": 3, "G": 4, "T": 5}
    refb = ref[position - a].upper() if ref and 0 <= position - a < len(ref) else None
    c = cov[position]
    counts = {k: c[i] for k, i in idx.items()}
    if alt is None:
        alt = max((k for k in counts if k != refb), key=lambda k: counts[k], default=None)
    out = {"position": position, "ref": refb, "alt": alt, "depth": c[1]}
    if not refb or not alt or alt not in idx:
        return out
    f = fwd.get(position)
    fi = {"A": 1, "C": 2, "G": 3, "T": 4}
    ref_f = f[fi[refb]] if f else None
    alt_f = f[fi[alt]] if f else None
    ref_r = counts[refb] - ref_f if ref_f is not None else None
    alt_r = counts[alt] - alt_f if alt_f is not None else None
    out.update(ref_forward=ref_f, ref_reverse=ref_r, alt_forward=alt_f, alt_reverse=alt_r)
    if None not in (ref_f, ref_r, alt_f, alt_r):
        try:
            from scipy.stats import fisher_exact
            out["fisher_p"] = float(fisher_exact([[ref_f, ref_r], [alt_f, alt_r]])[1])
        except Exception:
            out["fisher_p"] = None
        # GATK's StrandOddsRatio (pseudocount 1)
        rf, rr, af, ar = ref_f + 1, ref_r + 1, alt_f + 1, alt_r + 1
        sym = (rf * ar) / (rr * af) + (rr * af) / (rf * ar)
        out["sor"] = math.log(sym) + math.log(min(rf, rr) / max(rf, rr)) - math.log(min(af, ar) / max(af, ar))
    # where the ALT base sits in the drawn reads
    dist = []
    for r in view.get("reads") or []:
        for pos, base, q in r.get("mismatches") or []:
            if pos == position and str(base).upper() == alt:
                dist.append(min(position - r["start"], r["end"] - position))
    if dist:
        out["alt_reads_drawn"] = len(dist)
        out["alt_median_distance_to_end"] = statistics.median(dist)
        out["alt_within_10_of_end"] = sum(1 for d in dist if d < 10) / len(dist)
    return out


def copy_number_profile(profile, controls_median, features=(), title="Copy number along the locus", width=900, height=330,
                        expected=2.0) -> str:
    """profile: (start, end, mean depth) bins; copies = expected × depth /
    the median depth of control regions (diploid). Dashed lines at 1, 2 and
    3 copies; `features`: (name, start, end) boxes under the plot (genes,
    exons). A whole gene at 1 copy: a heterozygous deletion; at 3: a
    duplication; only part of a gene changed: a hybrid allele."""
    pts = [(s0, e0, expected * d / controls_median) for s0, e0, d in profile if controls_median]
    s = Svg(width, height, title)
    s.text(10, 18, title, size=13, weight=600)
    if not pts:
        s.text(10, 44, "No depth to plot (no control depth).", opacity=0.7)
        return s.render()
    x0, x1 = min(p[0] for p in pts), max(p[1] for p in pts)
    ymax = max(4.2, max(p[2] for p in pts) * 1.1)
    feat_h = 22 * (1 if features else 0)
    ax = Axes(s, 60, 34, width - 80, height - 84 - feat_h, x0, x1, 0, ymax)
    ax.draw("", "copies", xt=[])
    for t in ticks(x0, x1, 7):
        if x0 <= t <= x1:
            s.text(ax.X(t), ax.T + ax.H + 14, f"{int(t):,}", "middle", 10, opacity=0.65)
    for c in (1, 2, 3, 4):
        if c < ymax:
            s.line(ax.L, ax.Y(c), ax.L + ax.W, ax.Y(c), NEUTRAL, 1 if c != 2 else 1.4, 0.7 if c != 2 else 0.95, "4 3")
            s.text(ax.L + ax.W + 4, ax.Y(c) + 4, str(c), size=10, opacity=0.7)
    for s0, e0, cn in pts:
        if cn <= 0:
            # no read at all: outside the data (a cut BAM) or unmappable, not a deletion
            s.rect(ax.X(s0), ax.T, max(ax.X(e0) - ax.X(s0) - 0.6, 0.8), ax.H, NEUTRAL, 0.12, tip=f"{s0:,}-{e0:,}: no read")
            continue
        col = SERIES[0] if 1.5 <= cn < 2.5 else (SERIES[7] if cn < 1.5 else SERIES[1])
        s.rect(ax.X(s0), ax.Y(min(cn, ymax)), max(ax.X(e0) - ax.X(s0) - 0.6, 0.8), ax.Y(0) - ax.Y(min(cn, ymax)), col, 0.85,
               tip=f"{s0:,}-{e0:,}: {cn:.2f} copies")
    y = ax.T + ax.H + 24
    for i, (name, fs, fe) in enumerate(features):
        xa, xb = ax.X(max(fs, x0)), ax.X(min(fe, x1))
        s.rect(xa, y + (i % 2) * 0, max(xb - xa, 1), 10, SERIES[6], 0.55, rx=2, tip=f"{name} {fs:,}-{fe:,}")
        s.text((xa + xb) / 2, y + 20, name, "middle", 10, opacity=0.85)
    legend(s, [("about 2 copies", SERIES[0], "bar"), ("fewer", SERIES[7], "bar"), ("more", SERIES[1], "bar"), ("no read", NEUTRAL, "bar")],
           width - 400, 20)
    return s.render()


# ------------------------------------------------------------------ contaminant screen (screen 31)
CALL_COLOURS = {"detected": STATUS["critical"], "background": NEUTRAL, "below_threshold": STATUS["warning"], "not_detected": NEUTRAL}


def genome_windows(organisms, title="Where the reads fall on each genome", width=860, max_rows=24) -> str:
    """organisms: dicts of a contaminant_screen sample (name, call, role,
    window_map, genome_windows, window_bp, unique_reads...). One strip per
    organism, one cell per genome window: dark where its own reads fell, light
    where only shared reads did. Reads of a real infection are spread over the
    genome; a handful of windows lit in one place is a vector element, an
    integrated fragment or a gene shared with other organisms."""
    rows = [o for o in organisms if o.get("window_map") and ("2" in o["window_map"] or "1" in o["window_map"])][:max_rows]
    lab_w = min(300, 14 + 6.3 * max((len(o["name"]) for o in rows), default=8))
    row = 26
    h = 56 + row * max(len(rows), 1) + 24
    s = Svg(width, h, title)
    s.text(10, 18, title, size=13, weight=600)
    if not rows:
        s.text(10, 44, "No read fell on any agent or control.", opacity=0.7)
        return s.render()
    L, R = lab_w + 8, 150
    pw = width - L - R
    for i, o in enumerate(rows):
        y = 44 + i * row
        m = o["window_map"]
        n = len(m)
        call = o.get("call", "")
        label = o["name"]
        s.text(L - 8, y + 12, label, "end", 11, opacity=0.9)
        s.rect(L, y + 3, pw, 12, NEUTRAL, 0.18, rx=2, tip=f"{o['name']}: {n} windows of {o.get('window_bp', '?')} bp")
        cw = pw / max(n, 1)
        j = 0
        while j < n:
            c = m[j]
            k = j
            while k < n and m[k] == c:
                k += 1
            if c != "0":
                col = SERIES[0] if c == "2" else LIGHT_BLUE
                s.rect(L + j * cw, y + 3, max((k - j) * cw, 1.2), 12, col, 0.95,
                       tip=f"{o['name']}: windows {j + 1}–{k} " + ("(its own reads)" if c == "2" else "(shared reads only)"))
            j = k
        hit = sum(1 for ch in m if ch != "0")
        tag = "control" if o.get("role") == "control" else call.replace("_", " ")
        s.rect(L + pw + 10, y + 4, 10, 10, CALL_COLOURS.get(call, NEUTRAL) if o.get("role") != "control" else STATUS["good"], 1, rx=2)
        s.text(L + pw + 26, y + 13, f"{tag}: {hit}/{n} windows", size=10, opacity=0.85)
    legend(s, [("own reads", SERIES[0], "bar"), ("shared reads only", LIGHT_BLUE, "bar"), ("no read", NEUTRAL, "bar")], L, 36 - 4)
    return s.render()


# ------------------------------------------------------------------ generic small charts
def histogram(values, bins=20, title="Histogram", x_label="", lo=None, hi=None, width=460, height=320, ref_line=None) -> str:
    vals = [v for v in values if v is not None and math.isfinite(v)]
    s = Svg(width, height, title)
    s.text(10, 18, title, size=13, weight=600)
    if not vals:
        s.text(10, 44, "No values.", opacity=0.7)
        return s.render()
    lo = min(vals) if lo is None else lo
    hi = max(vals) if hi is None else hi
    if hi == lo:
        hi = lo + 1
    counts = [0] * bins
    for v in vals:
        i = min(bins - 1, max(0, int((v - lo) / (hi - lo) * bins)))
        counts[i] += 1
    ax = Axes(s, 60, 34, width - 80, height - 84, lo, hi, 0, max(counts) * 1.12)
    ax.draw(x_label, "count")
    bw = ax.W / bins
    for i, c in enumerate(counts):
        y = ax.Y(c)
        a = lo + (hi - lo) * i / bins
        s.rect(ax.L + i * bw + 1, y, bw - 2, ax.T + ax.H - y, SERIES[0], 0.85, rx=0,
               tip=f"{fmt(a)} to {fmt(a + (hi - lo) / bins)}: {c:,}")
    if ref_line is not None:
        s.line(ax.L, ax.Y(ref_line), ax.L + ax.W, ax.Y(ref_line), NEUTRAL, 1.2, 0.9, "4 3")
    return s.render()


def heatmap(matrix, row_labels, col_labels, title="Heatmap", width=520, value_label="", reverse=False) -> str:
    """A sequential (one-hue) heatmap; each cell's value is in its tooltip."""
    n, m = len(row_labels), len(col_labels)
    lab_w = min(160, 12 + 6.4 * max((len(str(r)) for r in row_labels), default=4))
    cell = max(10, min(36, (width - lab_w - 40) / max(m, 1)))
    top_h = min(140, 12 + 6.4 * max((len(str(c)) for c in col_labels), default=4))
    h = 40 + top_h + cell * n + 30
    s = Svg(int(lab_w + cell * m + 40), int(h), title)
    s.text(10, 18, title, size=13, weight=600)
    vals = [v for row in matrix for v in row if v is not None and math.isfinite(v)]
    lo, hi = (min(vals), max(vals)) if vals else (0, 1)
    ramp = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
    for j, c in enumerate(col_labels):
        x = lab_w + j * cell + cell / 2
        s.text(x, 34 + top_h - 6, c, "start", 10, opacity=0.8, rotate=-60)
    for i, r in enumerate(row_labels):
        y = 34 + top_h + i * cell
        s.text(lab_w - 6, y + cell / 2 + 4, r, "end", 10, opacity=0.8)
        for j in range(m):
            v = matrix[i][j]
            if v is None:
                continue
            t = 0 if hi == lo else (v - lo) / (hi - lo)
            t = 1 - t if reverse else t
            col = ramp[min(len(ramp) - 1, int(t * (len(ramp) - 1) + 0.5))]
            s.rect(lab_w + j * cell + 1, y + 1, cell - 2, cell - 2, col, 1, rx=2, tip=f"{r} × {col_labels[j]}: {fmt(v)}")
    s.text(lab_w, h - 8, f"{value_label}: light {fmt(lo if not reverse else hi)} → dark {fmt(hi if not reverse else lo)}", size=10, opacity=0.7)
    return s.render()
