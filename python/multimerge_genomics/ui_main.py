import marimo

__generated_with = "0.23.11"
# app_title is the browser tab title.
app = marimo.App(width="full", app_title="MMGEN Suite")


@app.cell
def imports():
    import marimo as mo
    import polars as pl
    import numpy as np
    import os
    import asyncio
    import re
    import multimerge_genomics as _engine

    # Every call to the engine goes through the session record: the AI
    # explains "the last result of this screen" from it, and screen 29 saves
    # the session as a Python script. An engine from before the AI layer is
    # used as it is.
    try:
        from multimerge_genomics.ai import Recorder as _Recorder
        mmg = _Recorder(_engine)
    except ImportError:
        mmg = _engine
    # The diagnostic figures (SVG), shared by the screens and usable in scripts.
    from multimerge_genomics import figures as figs
    return asyncio, figs, mmg, mo, np, os, pl, re


@app.cell
def global_state(mmg, mo):
    # allow_self_loops: the cell that calls the setter is re-run too. Without it
    # the menu (which calls set_page) never left the page when a screen opened,
    # and the screen you left through "Back to Menu" stayed on the page.
    get_page, _set_page = mo.state("home", allow_self_loops=True)

    def set_page(page):
        # the session record notes the screen each call is run from
        if hasattr(mmg, "context"):
            mmg.context = page
        _set_page(page)

    # The AI settings (screen 29). A plain dict, read when a question is sent,
    # so that changing a setting never sends anything by itself; the key is
    # kept in memory only, never written to disk. get_ai_rev redraws what
    # shows the settings.
    ai_settings = {"key": "", "model": "claude-opus-5-5", "effort": "medium", "level": "student",
                   "consent": False, "human": False}
    get_ai_rev, set_ai_rev = mo.state(0)
    # The agent of screen 29, as (agent, round counter): a new tuple redraws it.
    get_agent, set_agent = mo.state(None)

    # Message of the "Generate simulated patient" button (screen 4). It lives
    # here and not in the screen's own cell: a state created in the very cell
    # that shows it is reborn empty on every run.
    get_patient_msg, set_patient_msg = mo.state("", allow_self_loops=True)
    return ai_settings, get_agent, get_ai_rev, get_page, get_patient_msg, set_agent, set_ai_rev, set_page, set_patient_msg


@app.cell
def file_helpers(mmg, mo, os, re):
    # Shared by the screens: the files of the folder each screen offers (and
    # those of the examples folder), the default output names, a look at the
    # first records of a result, and the Python call that does the same thing
    # (shown under every result of screens 9-14, so each screen is also a
    # recipe for a script).
    NO_FILE = "(no matching file in this folder)"
    EXAMPLES_DIR = "mmgen_examples"
    FILE_KINDS = {
        "vcf": (".vcf", ".vcf.gz", ".vcf.bgz", ".vcf.lz4", ".vcf.parquet"),
        "bed": (".bed", ".bed.gz", ".bed.bgz", ".bed.lz4", ".bed.parquet"),
        "gff": (".gff", ".gff.gz", ".gff.lz4", ".gff.parquet", ".gff3", ".gff3.gz", ".gff3.lz4", ".gff3.parquet",
                ".gtf", ".gtf.gz", ".gtf.lz4", ".gtf.parquet"),
        "sam": (".sam", ".sam.gz", ".sam.lz4", ".sam.parquet", ".bam"),
        "fasta": (".fa", ".fa.gz", ".fa.lz4", ".fasta", ".fasta.gz", ".fasta.lz4", ".fna", ".fna.gz", ".fna.lz4"),
        "fastq": (".fq", ".fq.gz", ".fq.lz4", ".fastq", ".fastq.gz", ".fastq.lz4"),
        "genome": (".fai", ".genome", ".sizes"),
        "mapping": (".txt", ".tsv", ".map"),
        "indexed": (".vcf.gz", ".vcf.bgz", ".bed.gz", ".bed.bgz", ".gff.gz", ".gff3.gz", ".gtf.gz", ".bam"),
        "parquet": (".parquet",),
        # Screens 1-8.
        "kmer_input": (".fa", ".fa.gz", ".fasta", ".fasta.gz", ".fna", ".fna.gz", ".fq", ".fq.gz", ".fastq", ".fastq.gz"),
        "plain_fasta": (".fa", ".fasta", ".fna"),
        "arrow": (".arrow",),
        "kmer_text": (".txt",),
        "vcf_merge": (".vcf", ".vcf.gz"),
        # Human-genome screens.
        "bed_any": (".bed", ".bed.gz"),
        "jar": (".jar",),
        "bam": (".bam",),
        "table": (".tsv", ".tsv.gz", ".txt", ".txt.gz", ".csv", ".csv.gz"),
        "sig": (".sig",),
        # Research screens.
        "sheet": (".tsv", ".txt", ".files", ".csv"),
        "model": (".model",),
        "geneset": (".gmt", ".gmt.gz"),
        "tx_index": (".idx",),
        "screen_db": (".mmgscreen",),
    }
    SKIP_DIRS = {"target", ".git", ".venv", "venv", "node_modules", "__pycache__", "__marimo__", "dist", "dist_colab", "dist_windows", "dist_windows_cuda"}

    def files_under(folder, *kinds, depth=4, limit=4000):
        """The files under `folder` (and its subfolders, a few levels deep) with
        an extension of `kinds`, then those of the examples folder."""
        ends = tuple(s for k in kinds for s in FILE_KINDS[k])
        found = []
        root = (folder or ".").strip().rstrip("/\\") or "."

        def walk(path, level):
            try:
                entries = sorted(os.scandir(path), key=lambda e: e.name.lower())
            except OSError:
                return
            for e in entries:
                if len(found) >= limit:
                    return
                if e.is_dir(follow_symlinks=False):
                    if level < depth and e.name not in SKIP_DIRS and not e.name.startswith(("_backup", ".")):
                        walk(e.path, level + 1)
                elif e.name.lower().endswith(ends):
                    rel = e.path if root != "." else os.path.relpath(e.path, ".")
                    found.append(rel.replace("\\", "/"))

        if os.path.isdir(root):
            walk(root, 0)
        if os.path.isdir(EXAMPLES_DIR) and not os.path.abspath(EXAMPLES_DIR).startswith(os.path.abspath(root) + os.sep):
            walk(EXAMPLES_DIR, 0)
        return found or [NO_FILE]
    COMPRESSIONS = (".gz", ".bgz", ".lz4", ".parquet")
    FORMAT_EXTENSIONS = ("vcf", "bed", "gff", "gff3", "gtf", "sam", "bam", "fa", "fasta", "fna", "fq", "fastq", "tsv", "txt")

    def files_with(*kinds):
        """The files of the current folder with an extension of `kinds`, then
        those of the examples folder (as 'mmgen_examples/name')."""
        ends = tuple(s for k in kinds for s in FILE_KINDS[k])
        found = sorted(n for n in os.listdir(".") if n.lower().endswith(ends) and os.path.isfile(n))
        if os.path.isdir(EXAMPLES_DIR):
            found += sorted(f"{EXAMPLES_DIR}/{n}" for n in os.listdir(EXAMPLES_DIR)
                            if n.lower().endswith(ends) and os.path.isfile(os.path.join(EXAMPLES_DIR, n)))
        return found or [NO_FILE]

    def real_files(*kinds):
        return [n for n in files_with(*kinds) if n != NO_FILE]

    def chosen(value):
        """A real file was picked, not the "no file" placeholder."""
        return bool(value) and value != NO_FILE

    def split_name(path):
        """'data/x.vcf.gz' -> ('x', '.vcf.gz'); 'x.bam' -> ('x', '.bam')."""
        base, ext = os.path.basename(path), ""
        for comp in COMPRESSIONS:
            if base.lower().endswith(comp):
                base, ext = base[: -len(comp)], base[-len(comp):]
                break
        stem, dot, fmt = base.rpartition(".")
        if dot and stem and fmt.lower() in FORMAT_EXTENSIONS:
            return stem, "." + fmt + ext
        return base, ext

    def beside(source, name):
        """`name` in the folder of `source` (a result goes next to its input)."""
        folder = os.path.dirname(source)
        return f"{folder}/{name}" if folder else name

    def output_name(typed, source, tag, ext=None):
        """The name typed in the form, or <source>.<tag><extension> when it was left blank."""
        typed = (typed or "").strip()
        if typed:
            return typed
        stem, source_ext = split_name(source)
        return beside(source, f"{stem}.{tag}{source_ext if ext is None else ext}")

    def clashes(output, inputs):
        """The output would overwrite one of the inputs."""
        out = os.path.abspath(output)
        return any(os.path.abspath(p) == out for p in inputs if p)

    def size_of(path):
        n = float(os.path.getsize(path))
        for unit in ("bytes", "KB", "MB", "GB"):
            if n < 1024 or unit == "GB":
                return f"{n:,.0f} {unit}" if unit == "bytes" else f"{n:,.1f} {unit}"
            n /= 1024

    def facts(pairs):
        """A Markdown list of `label: value` lines."""
        return "\n".join(f"* **{label}:** {value}" for label, value in pairs)

    def warn(text):
        return mo.callout(mo.md("⚠️ " + text), kind="warn")

    def error_panel(err):
        return mo.callout(mo.md(f"**🔥 ERROR:** {err}"), kind="danger")

    def skipped_panel(paths, skipped, kept):
        """Why the engine skipped lines: the validator's first errors. Without it
        a file in the wrong format just gave a green panel full of zeros."""
        if not skipped:
            return None
        found = []
        for path in [paths] if isinstance(paths, str) else paths:
            try:
                # the engine itself: this check is not a call the session should record
                report = getattr(mmg, "_module", mmg).validate_file(path, max_lines=100_000)
            except Exception:
                continue
            for p in report["problems"]:
                if p["severity"] == "error" and len(found) < 3:
                    found.append(f"* `{os.path.basename(path)}`, line {p['line']:,}: {p['message']}"
                                 + (f"  \n  ↳ fix: {p['fix']}" if p["fix"] else ""))
        if kept:
            head = f"**{skipped:,} line(s) were skipped**: they are not valid records."
        else:
            head = f"**Nothing was read: all {skipped:,} line(s) were skipped**, because none is a valid record."
        why = "\n\nWhat the validator says:\n\n" + "\n".join(found) if found else ""
        return mo.callout(mo.md(head + why + "\n\nScreen **10** (Validate a file) gives the full report."),
                          kind="warn" if kept else "danger")

    # Example gene panels, for teaching (the same lists as the demo data's
    # panels.json). Clinical panels are curated by each laboratory.
    GENE_PANELS = {
        "Hereditary breast and ovarian cancer": ["BRCA1", "BRCA2", "PALB2", "TP53", "PTEN", "CDH1", "STK11", "ATM", "CHEK2", "BARD1", "BRIP1", "RAD51C", "RAD51D"],
        "Lynch syndrome (hereditary colorectal cancer)": ["MLH1", "MSH2", "MSH6", "PMS2", "EPCAM"],
        "Familial hypercholesterolemia": ["LDLR", "APOB", "PCSK9", "LDLRAP1"],
        "Cardiomyopathy": ["MYH7", "MYBPC3", "TNNT2", "TNNI3", "TPM1", "ACTC1", "MYL2", "MYL3", "LMNA", "TTN", "PLN", "DSP", "PKP2"],
        "Inherited arrhythmia": ["KCNQ1", "KCNH2", "SCN5A", "RYR2", "CASQ2"],
        "Bleeding and clotting": ["F2", "F5", "F8", "F9", "F11", "VWF"],
        "Carrier screening (common recessive conditions)": ["CFTR", "HBB", "HBA1", "HBA2", "HEXA", "GBA1", "ASPA", "FANCC", "BLM", "PAH", "SERPINA1", "GJB2", "MEFV"],
        "Pharmacogenes": ["CYP2C19", "CYP2C9", "CYP2D6", "CYP3A5", "CYP4F2", "DPYD", "TPMT", "NUDT15", "SLCO1B1", "UGT1A1", "VKORC1", "IFNL3", "NAT2", "CYP2B6", "RYR1", "CACNA1S"],
    }

    def bed_next_to(vcf):
        """The callable-regions BED of a VCF: x.callable.bed next to x.vcf.gz,
        or the Genome in a Bottle *_benchmark*.bed of the same name."""
        folder, stem = os.path.dirname(vcf) or ".", split_name(vcf)[0]
        cands = [os.path.join(folder, stem + ".callable.bed")]
        if "_benchmark" in stem:
            prefix = stem.split("_benchmark")[0] + "_benchmark"
            cands += [os.path.join(folder, prefix + "_noinconsistent.bed"), os.path.join(folder, prefix + ".bed")]
        return next((c.replace("\\", "/") for c in cands if os.path.exists(c)), None)

    def pair_fastqs(paths):
        """FASTQ files grouped into samples: x_R1.fastq.gz goes with
        x_R2.fastq.gz (also _1/_2, .R1/.R2, _R1_001); any other file is
        single-end. Returns [(sample name, r1, r2 or None)], in the order given."""
        mate = re.compile(r"^(.*?)[._-]R?([12])(?:_001)?\.(?:fastq|fq)(?:\.gz|\.lz4)?$", re.I)
        groups, order = {}, []
        for path in paths:
            m = mate.match(os.path.basename(path))
            name, which = (m.group(1), m.group(2)) if m else (split_name(path)[0], "1")
            key = (os.path.dirname(path), name)
            if key not in groups or which in groups[key]:
                key = key if key not in groups else (key[0], f"{name}_{len(order) + 1}")
                groups[key] = {}
                order.append(key)
            groups[key][which] = path
        return [(key[1], g.get("1") or g.get("2"), g.get("2") if "1" in g else None) for key, g in ((k, groups[k]) for k in order)]

    def words(text):
        """'BRCA1, BRCA2; TP53' -> ['BRCA1', 'BRCA2', 'TP53']."""
        return [w for w in (text or "").replace(",", " ").replace(";", " ").split() if w]

    def human_notice():
        """The notice on every screen that reads a person's genome."""
        return mo.callout(mo.md(
            "**⚠️ Experimental software, for teaching and research only.** These results are not a medical diagnosis, "
            "not a paternity test with legal value and not advice on any treatment or medicine. The authors and "
            "contributors take no responsibility for decisions made on the basis of them. A finding must be confirmed "
            "by an accredited laboratory and interpreted by a qualified professional (a physician or a genetic "
            "counselor). Analyze only the genomes of people who have agreed to it: genetic data are sensitive "
            "personal data, about the person and about their relatives."), kind="danger")

    def svg_download(svg, name, title=None):
        """A button that saves a chart as a standalone SVG file (vector, for a paper)."""
        body = svg if "xmlns=" in svg else svg.replace("<svg ", '<svg xmlns="http://www.w3.org/2000/svg" ', 1)
        if title and "<title>" not in body[:300]:
            body = body.replace(">", f"><title>{title}</title>", 1)
        safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)[:60] or "figure"
        return mo.download(data=body.encode("utf-8"), filename=f"{safe}.svg", mimetype="image/svg+xml", label="⬇️ SVG")

    def figure_panel(svg, name, note=None):
        """A figure of `figs` with a button that saves it as SVG (vector, for a
        paper) and, below it, what to look for."""
        if not svg:
            return None
        safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)[:60] or "figure"
        parts = [mo.Html(svg), mo.download(data=svg.encode("utf-8"), filename=f"{safe}.svg", mimetype="image/svg+xml",
                                           label="⬇️ SVG")]
        if note:
            parts.append(mo.md(f"*{note}*"))
        return mo.vstack(parts, gap=0.25)

    def make_form(widgets, submit):
        """A form with one field per line and plain labels (no JSON tree)."""
        template = "\n\n".join("{" + key + "}" for key in widgets)
        return mo.md(template).batch(**widgets).form(submit_button_label=submit, bordered=True)

    def python_call(func, *args, **kwargs):
        parts = [repr(a) for a in args] + [f"{k}={v!r}" for k, v in kwargs.items()]
        return f"mmg.{func}({', '.join(parts)})"

    def code_panel(*lines):
        return mo.md("**The same from Python** (paste it into a script or a notebook):\n\n```python\n"
                     "import multimerge_genomics as mmg\n\n" + "\n".join(lines) + "\n```")

    def preview_panel(path, n=8):
        """The first records of a result written as text or gzip."""
        import gzip
        name = os.path.basename(path)
        low = path.lower()
        if low.endswith((".bam", ".parquet", ".lz4")):
            kind = "BAM" if low.endswith(".bam") else "Parquet" if low.endswith(".parquet") else "LZ4"
            return mo.md(f"*`{name}` is a {kind} file ({size_of(path)}): this preview only opens text and gzip.*")
        is_sam = ".sam" in low
        rows = []
        try:
            opener = gzip.open if low.endswith((".gz", ".bgz")) else open
            with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    line = line.rstrip("\r\n")
                    if line.startswith("#") or (is_sam and line.startswith("@")):
                        continue
                    rows.append(line if len(line) <= 160 else line[:157] + "...")
                    if len(rows) >= n:
                        break
        except Exception as err:
            return mo.md(f"*Could not preview `{name}`: {err}*")
        if not rows:
            return mo.md(f"*`{name}` has no records.*")
        return mo.vstack([mo.md(f"**The first {len(rows)} line(s) of `{name}`:**"), mo.plain_text("\n".join(rows))])

    def scatter_chart(points, x_label, y_label, title, labels=(), hlines=(), vlines=(), max_background=4000, names=None):
        """A scatter plot as SVG: `points` are (x, y, kind, tip) with kind "up",
        "down" or "ns" (not significant, drawn grey underneath; sampled when
        there are many). `labels`: (x, y, text) written next to their point.
        `names`: the legend of each kind, when not up / down / not significant.
        Colours of the FASTQ report (validated); text in currentColor."""
        import math
        pts = [p for p in points if p[0] is not None and p[1] is not None and math.isfinite(p[0]) and math.isfinite(p[1])]
        if not pts:
            return None
        ns = [p for p in pts if p[2] == "ns"]
        if len(ns) > max_background:
            step = len(ns) / max_background
            ns = [ns[int(i * step)] for i in range(max_background)]
        shown = ns + [p for p in pts if p[2] != "ns"]
        xs, ys = [p[0] for p in shown], [p[1] for p in shown]
        x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
        if x0 == x1:
            x0, x1 = x0 - 1, x1 + 1
        if y0 == y1:
            y0, y1 = y0 - 1, y1 + 1
        padx, pady = 0.04 * (x1 - x0), 0.06 * (y1 - y0)
        x0, x1, y0, y1 = x0 - padx, x1 + padx, y0 - pady, y1 + pady
        W, H, L, R, T, B = 640, 380, 58, 16, 16, 44
        pw, ph = W - L - R, H - T - B
        X = lambda v: L + (v - x0) / (x1 - x0) * pw
        Y = lambda v: T + ph - (v - y0) / (y1 - y0) * ph

        def ticks(a, b):
            raw = (b - a) / 6
            mag = 10 ** math.floor(math.log10(raw))
            step = next(m * mag for m in (1, 2, 5, 10) if m * mag >= raw)
            first = math.ceil(a / step) * step
            return [first + i * step for i in range(int((b - first) / step) + 1)]

        parts = []
        for v in ticks(y0, y1):
            parts.append(f'<line x1="{L}" x2="{W - R}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" stroke="currentColor" stroke-opacity="0.12"/>')
            parts.append(f'<text x="{L - 6}" y="{Y(v) + 4:.1f}" text-anchor="end" fill="currentColor" fill-opacity="0.7">{v:g}</text>')
        for v in ticks(x0, x1):
            parts.append(f'<text x="{X(v):.1f}" y="{T + ph + 16}" text-anchor="middle" fill="currentColor" fill-opacity="0.7">{v:g}</text>')
        for v in hlines:
            if y0 < v < y1:
                parts.append(f'<line x1="{L}" x2="{W - R}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" stroke="currentColor" stroke-opacity="0.45" stroke-dasharray="4 3"/>')
        for v in vlines:
            if x0 < v < x1:
                parts.append(f'<line y1="{T}" y2="{T + ph}" x1="{X(v):.1f}" x2="{X(v):.1f}" stroke="currentColor" stroke-opacity="0.45" stroke-dasharray="4 3"/>')
        colour = {"up": "#eb6834", "down": "#2a78d6", "ns": "#898781"}
        for x, y, kind, tip in ns:
            parts.append(f'<circle cx="{X(x):.1f}" cy="{Y(y):.1f}" r="2" fill="{colour["ns"]}" fill-opacity="0.35"/>')
        for i, (x, y, kind, tip) in enumerate(p for p in shown if p[2] != "ns"):
            t = f"<title>{tip}</title>" if tip and i < 1500 else ""
            parts.append(f'<circle cx="{X(x):.1f}" cy="{Y(y):.1f}" r="3" fill="{colour[kind]}" fill-opacity="0.85">{t}</circle>')
        placed = []  # label boxes already drawn: a label that would overlap one moves down, or is left out
        for x, y, text in labels:
            if not (x0 < x < x1 and y0 < y < y1):
                continue
            anchor = "end" if X(x) > L + 0.75 * pw else "start"
            dx = -6 if anchor == "end" else 6
            width = 7.2 * len(str(text))
            left = X(x) + dx - (width if anchor == "end" else 0)
            for shift in (0, 13, 26, -13):
                top = Y(y) - 4 - 11 + shift
                box = (left, top, left + width, top + 13)
                if all(box[2] < b[0] or box[0] > b[2] or box[3] < b[1] or box[1] > b[3] for b in placed):
                    placed.append(box)
                    parts.append(f'<text x="{X(x) + dx:.1f}" y="{Y(y) - 4 + shift:.1f}" text-anchor="{anchor}" fill="currentColor" font-weight="600">{text}</text>')
                    break
        parts.append(f'<text x="{L + pw / 2:.1f}" y="{H - 6}" text-anchor="middle" fill="currentColor" fill-opacity="0.8">{x_label}</text>')
        parts.append(f'<text transform="translate(14,{T + ph / 2:.1f}) rotate(-90)" text-anchor="middle" fill="currentColor" fill-opacity="0.8">{y_label}</text>')
        n_up = sum(1 for p in pts if p[2] == "up")
        n_down = sum(1 for p in pts if p[2] == "down")
        n_ns = sum(1 for p in pts if p[2] == "ns")
        name = dict({"up": "up", "down": "down", "ns": "not significant"}, **(names or {}))
        legend = (f'<div style="font: 12px system-ui, sans-serif; display: flex; flex-wrap: wrap; gap: 4px 16px; margin: 4px 0;">'
                  + (f'<span><span style="color:{colour["up"]}">●</span> {name["up"]} ({n_up:,})</span>' if n_up or not names else "")
                  + (f'<span><span style="color:{colour["down"]}">●</span> {name["down"]} ({n_down:,})</span>' if n_down or not names else "")
                  + (f'<span><span style="color:{colour["ns"]}">●</span> {name["ns"]} ({n_ns:,}'
                     + (f", {len(ns):,} drawn" if len(ns) < n_ns else "") + ")</span>" if n_ns or not names else "")
                  + "</div>")
        svg = (f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width: {W}px; font: 12px system-ui, sans-serif;" role="img" '
               f'aria-label="{title}">' + "".join(parts) + "</svg>")
        return mo.vstack([mo.md(f"**{title}**"), mo.Html(legend + svg), svg_download(svg, title.split(":")[0].lower(), title)], gap=0.25)

    def spectrum_chart(histogram):
        """Bars of the abundance spectrum: how many distinct k-mers were seen once,
        twice... The y axis is logarithmic: the k-mers seen once (in reads, mostly
        sequencing errors) outnumber the rest by orders of magnitude."""
        import math
        points = [(int(c), int(n)) for c, n in histogram if int(n) > 0]
        if not points:
            return None
        total = sum(n for _, n in points)
        x_max, running = points[-1][0], 0
        for c, n in points:
            running += n
            if running >= 0.999 * total:
                x_max = c
                break
        x_max = max(10, min(300, x_max))
        shown = {c: n for c, n in points if c <= x_max}
        beyond = total - sum(shown.values())
        y_top = 10 ** max(1, math.ceil(math.log10(max(shown.values()))))
        W, H, LEFT, RIGHT, TOP, BOTTOM = 720, 290, 64, 16, 28, 46
        pw, ph = W - LEFT - RIGHT, H - TOP - BOTTOM
        base = TOP + ph

        def y_of(v):
            return base - ph * math.log10(v) / math.log10(y_top)

        def short(v):
            for div, suffix in ((10**9, "G"), (10**6, "M"), (10**3, "k")):
                if v >= div:
                    return f"{v // div}{suffix}"
            return str(v)

        slot = pw / x_max
        bar = max(1.0, slot - 2)  # a 2 px gap between bars
        parts = []
        for i in range(0, round(math.log10(y_top)) + 1):
            yy = y_of(10**i)
            parts.append(f'<line x1="{LEFT}" x2="{W - RIGHT}" y1="{yy:.1f}" y2="{yy:.1f}" stroke="currentColor" '
                         f'stroke-opacity="{0.4 if i == 0 else 0.12}"/>')
            parts.append(f'<text x="{LEFT - 8}" y="{yy + 4:.1f}" text-anchor="end" fill="currentColor" fill-opacity="0.7">{short(10**i)}</text>')
        step = next(s for s in (1, 2, 5, 10, 20, 25, 50, 100) if x_max / s <= 12)
        for c in sorted({1, *range(step, x_max + 1, step)}):
            parts.append(f'<text x="{LEFT + (c - 0.5) * slot:.1f}" y="{base + 17}" text-anchor="middle" fill="currentColor" fill-opacity="0.7">{c}</text>')
        for c, n in shown.items():
            x0 = LEFT + (c - 1) * slot + (slot - bar) / 2
            top = min(y_of(n), base - 1.5)  # a single k-mer still shows a stub
            r = 4 if bar >= 8 and base - top >= 4 else 0
            if r:
                d = (f"M{x0:.1f},{base:.1f} V{top + r:.1f} Q{x0:.1f},{top:.1f} {x0 + r:.1f},{top:.1f} "
                     f"H{x0 + bar - r:.1f} Q{x0 + bar:.1f},{top:.1f} {x0 + bar:.1f},{top + r:.1f} V{base:.1f} Z")
            else:
                d = f"M{x0:.1f},{base:.1f} V{top:.1f} H{x0 + bar:.1f} V{base:.1f} Z"
            parts.append(f'<path d="{d}" fill="#2a78d6"/>')
            # The hover target is the whole column, wider and taller than the bar.
            parts.append(f'<rect x="{LEFT + (c - 1) * slot:.1f}" y="{TOP}" width="{slot:.1f}" height="{ph}" fill="transparent">'
                         f'<title>{n:,} distinct k-mer(s) seen {c:,} time(s)</title></rect>')
        # The coverage peak of sequencing reads: the highest bar after the first valley.
        valley = next((c for c in range(1, x_max) if shown.get(c, 0) < shown.get(c + 1, 0)), None)
        if valley is not None:
            peak = max((c for c in shown if c > valley), key=lambda c: shown[c], default=None)
            if peak is not None and shown[peak] > shown.get(valley, 0):
                parts.append(f'<text x="{LEFT + (peak - 0.5) * slot:.1f}" y="{y_of(shown[peak]) - 8:.1f}" text-anchor="middle" '
                             f'fill="currentColor" font-weight="600">peak: {peak}×</text>')
        parts.append(f'<text x="{LEFT + pw / 2:.1f}" y="{H - 6}" text-anchor="middle" fill="currentColor" fill-opacity="0.8">times a k-mer was seen</text>')
        parts.append(f'<text transform="translate(14,{TOP + ph / 2:.1f}) rotate(-90)" text-anchor="middle" fill="currentColor" '
                     f'fill-opacity="0.8">distinct k-mers (log scale)</text>')
        svg = (f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width: {W}px; font: 12px system-ui, sans-serif;" role="img" '
               f'aria-label="Abundance spectrum: the number of distinct k-mers seen each number of times">' + "".join(parts) + "</svg>")
        items = [mo.md("**Abundance spectrum**"), mo.Html(svg), svg_download(svg, "kmer_spectrum", "Abundance spectrum")]
        if beyond:
            items.append(mo.md(f"*{beyond:,} distinct k-mer(s) seen more than {x_max:,} times are left out of the chart; the table below has them.*"))
        return mo.vstack(items)

    def tree_chart(tree, title, unit="Mash distance"):
        """A rectangular tree as SVG, from the layout the engine gives: x is
        the distance from the root, the leaves one per row in `leaf_order`,
        an inner node halfway between its children. Branches in currentColor
        (both themes), leaves as blue dots, a scale bar below."""
        import html as _html
        import math
        nodes = tree.get("nodes") or []
        leaves = list(tree.get("leaf_order") or [])
        if not nodes or not leaves:
            return None
        row = {n: i for i, n in enumerate(leaves)}
        ys = {}

        def y_of(i):
            stack = [i]
            while stack:  # without recursion: a tree of thousands of leaves is deep
                u = stack[-1]
                kids = nodes[u]["children"]
                if u in ys:
                    stack.pop()
                elif not kids:
                    ys[u] = row.get(u, 0)
                    stack.pop()
                elif all(k in ys for k in kids):
                    ys[u] = sum(ys[k] for k in kids) / len(kids)
                    stack.pop()
                else:
                    stack.extend(k for k in kids if k not in ys)
            return ys[i]

        depth = max(n["x"] for n in nodes) or 1.0
        longest = max(len(nodes[i]["name"]) for i in leaves)
        W, L, T, rowh = 760, 16, 14, 22
        R = min(330, 7 * longest + 20)
        H = T + rowh * len(leaves) + 46
        pw = W - L - R
        X = lambda v: L + v / depth * pw
        Y = lambda r: T + rowh * r + rowh / 2
        parts = []
        for i, n in enumerate(nodes):
            if n["parent"] is not None:
                px = nodes[n["parent"]]["x"]
                parts.append(f'<line x1="{X(px):.1f}" x2="{X(n["x"]):.1f}" y1="{Y(y_of(i)):.1f}" y2="{Y(y_of(i)):.1f}" '
                             f'stroke="currentColor" stroke-opacity="0.8" stroke-width="2"><title>branch {n["branch_length"]:.5g}</title></line>')
            if n["children"]:
                cy = [y_of(c) for c in n["children"]]
                parts.append(f'<line x1="{X(n["x"]):.1f}" x2="{X(n["x"]):.1f}" y1="{Y(min(cy)):.1f}" y2="{Y(max(cy)):.1f}" '
                             f'stroke="currentColor" stroke-opacity="0.8" stroke-width="2"/>')
        for i in leaves:
            n = nodes[i]
            name = _html.escape(n["name"])
            parts.append(f'<circle cx="{X(n["x"]):.1f}" cy="{Y(y_of(i)):.1f}" r="4" fill="#2a78d6"><title>{name}: {n["x"]:.5g} from the root</title></circle>')
            parts.append(f'<text x="{X(n["x"]) + 8:.1f}" y="{Y(y_of(i)) + 4:.1f}" fill="currentColor">{name}</text>')
        raw = depth / 4
        mag = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1
        bar = next(m * mag for m in (1, 2, 5, 10) if m * mag >= raw)
        yb = H - 22
        parts.append(f'<line x1="{L}" x2="{X(bar):.1f}" y1="{yb}" y2="{yb}" stroke="currentColor" stroke-width="2"/>')
        parts.append(f'<text x="{X(bar) + 6:.1f}" y="{yb + 4}" fill="currentColor" fill-opacity="0.8">{bar:g} ({unit})</text>')
        svg = (f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width: {W}px; font: 12px system-ui, sans-serif;" role="img" '
               f'aria-label="{_html.escape(title)}">' + "".join(parts) + "</svg>")
        return mo.vstack([mo.md(f"**{title}**"), mo.Html(svg), svg_download(svg, "tree", title)], gap=0.25)

    return (
        EXAMPLES_DIR,
        GENE_PANELS,
        NO_FILE,
        bed_next_to,
        beside,
        words,
        chosen,
        clashes,
        code_panel,
        error_panel,
        facts,
        files_under,
        files_with,
        human_notice,
        make_form,
        figure_panel,
        output_name,
        svg_download,
        pair_fastqs,
        preview_panel,
        python_call,
        real_files,
        scatter_chart,
        size_of,
        skipped_panel,
        spectrum_chart,
        split_name,
        tree_chart,
        warn,
    )


@app.cell
def example_files(EXAMPLES_DIR, mmg, mo, os, pl):
    # The examples of every screen. The files are small, written on demand into
    # mmgen_examples/ by make_examples (the same bytes every time: nothing is
    # downloaded), and each example is the set of values that fills in the form
    # of a screen and task.
    NO_EXAMPLE = "(none: fill in the form yourself)"
    EXAMPLES_VERSION = "6"
    HUMAN_DIR = EXAMPLES_DIR + "/human"
    RESEARCH_DIR = EXAMPLES_DIR + "/research"
    # The demo packs: (folder, version, zip name, URL where it is published —
    # None until it is, e.g. as a GitHub Release asset).
    # The demo packs: (folder, version, zip, where to download it, what it holds).
    # Downloaded from the GitHub release the first time an example needs them,
    # unless the zip is already in this folder. The 16S data (the mothur MiSeq
    # SOP) carry no explicit licence and are not published.
    _PACKS_URL = "https://github.com/fbcouto/multimerge-genomics-research/releases/download/v2.3.0/"
    DEMO_PACKS = {
        "human": (HUMAN_DIR, "3", "mmgen_demo_human_v3.zip", _PACKS_URL + "mmgen_demo_human_v3.zip",
                  "49 MB of public data cut from Genome in a Bottle, ClinVar, RefSeq, 1000 Genomes and the PGS Catalog"),
        "research": (RESEARCH_DIR, "4", "mmgen_demo_research_v4.zip", _PACKS_URL + "mmgen_demo_research_v4.zip",
                     "46 MB of public data: the airway RNA-seq experiment, the MAGeCK demos, reference sketches, "
                     "mock-community reads, the AMRFinderPlus catalog, four isolates, 20 SARS-CoV-2 lineages, a PCR panel "
                     "and tumor calls"),
        "pgx": (EXAMPLES_DIR + "/pgx", "1", "mmgen_demo_pgx_v1.zip", _PACKS_URL + "mmgen_demo_pgx_v1.zip",
                "33 MB of public data: the reads of the HG002 genome (Genome in a Bottle) over CYP2D6, the regions its copy "
                "number is measured against, and the HLA class I genes"),
        "amplicon16s": (EXAMPLES_DIR + "/research_16s", "2", "mmgen_demo_16s_v2.zip", _PACKS_URL + "mmgen_demo_16s_v2.zip",
                        "17 MB of public data: the 16S reads of the mothur MiSeq SOP (Kozich et al. 2013; redistributed with the "
                        "permission of Patrick D. Schloss) and an RDP classifier"),
    }
    EXAMPLE_FILES = (
        "toy_genome.fa", "toy_genome.fa.fai", "cohort.vcf", "cohort.vcf.gz", "sample_a.vcf", "sample_b.vcf.gz",
        "sample_c.vcf", "indels.vcf", "with_problems.vcf", "five_columns.vcf", "ensembl_names.vcf",
        "patient_grch38.vcf", "cohort_grch37.vcf", "grch38_ensembl.sizes", "genes.bed", "peaks.bed",
        "annotation.gff3", "alignments.sam", "alignments.sorted.bam", "pileup_demo.bam", "roh_demo.vcf", "reads_R1.fastq.gz", "reads_R2.fastq.gz",
        "rnaseq_ctrl_1.bam", "rnaseq_treat_2.bam", "rnaseq_samples.tsv", "rnaseq_counts.tsv",
        "toy_bacterium.fasta", "toy_patient.fasta", "toy_bacterium_reads.fastq.gz", "toy_bacterium_K31.arrow",
        "toy_patient_K31.arrow", "toy_bacterium_reads_K31.arrow", "kmers_gc_rich.txt", "kmers_at_rich.txt",
        "plasmid_contigs.fasta", "chromosome_contigs.fasta", "three_species.fasta", "README.md",
    )

    def make_examples():
        """Writes the example files into mmgen_examples/, unless they are all
        there already, from this version of the examples. Returns True when it
        wrote them."""
        import gzip
        import random

        # the engine itself, not the session record: writing the examples is
        # not something the person ran
        engine = getattr(mmg, "_module", mmg)
        folder = EXAMPLES_DIR
        stamp = os.path.join(folder, "VERSION")
        try:
            with open(stamp, encoding="utf-8") as fh:
                if fh.read().strip() == EXAMPLES_VERSION and all(os.path.exists(os.path.join(folder, n)) for n in EXAMPLE_FILES):
                    return False
        except OSError:
            pass
        os.makedirs(folder, exist_ok=True)
        rng = random.Random(2026)

        def path(name):
            return os.path.join(folder, name)

        def write(name, text):
            if name.endswith(".gz"):  # mtime=0: the same bytes on every run
                with open(path(name), "wb") as fh:
                    fh.write(gzip.compress(text.encode("utf-8"), mtime=0))
            else:
                with open(path(name), "w", encoding="utf-8", newline="\n") as fh:
                    fh.write(text)

        def dna(n, gc=0.5):
            return "".join(rng.choices("ACGT", weights=[1 - gc, gc, gc, 1 - gc], k=n))

        def revcomp(s):
            return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

        def other(base):
            return rng.choice([b for b in "ACGT" if b != base])

        def fasta_text(records):
            return "".join(f">{name}\n" + "".join(seq[i:i + 60] + "\n" for i in range(0, len(seq), 60)) for name, seq in records)

        # A small genome with three chromosomes (chr10 shows lexical against
        # natural order) and three repeats where an indel can slide.
        genome = {"chr1": dna(50_000), "chr2": dna(30_000), "chr10": dna(20_000)}

        def plant(chrom, pos, motif):  # pos is 1-based
            s = genome[chrom]
            genome[chrom] = s[:pos - 1] + motif + s[pos - 1 + len(motif):]

        def ref_at(chrom, pos, n=1):
            return genome[chrom][pos - 1:pos - 1 + n]

        plant("chr1", 20_000, "GTTTTTTC")
        plant("chr1", 30_000, "G" + "CA" * 6 + "T")
        plant("chr2", 12_000, "C" + "AG" * 4 + "T")
        write("toy_genome.fa", fasta_text(genome.items()))
        fai, offset = [], 0
        for name, seq in genome.items():
            offset += len(name) + 2
            fai.append(f"{name}\t{len(seq)}\t{offset}\t60\t61\n")
            offset += len(seq) + (len(seq) + 59) // 60
        write("toy_genome.fa.fai", "".join(fai))

        # Variants on that genome: REF is always the genome's base.
        contigs = "".join(f"##contig=<ID={c},length={len(s)}>\n" for c, s in genome.items())
        head = ("##fileformat=VCFv4.2\n##source=MMGEN examples\n##reference=toy_genome.fa\n" + contigs
                + '##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
                + '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n'
                + '##FILTER=<ID=q10,Description="Quality below 10">\n'
                + '##FILTER=<ID=LowDP,Description="Depth below 8">\n'
                + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                + '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
        columns_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"

        def site(chrom, pos):
            ref = ref_at(chrom, pos)
            kind = rng.choices(["snp", "del", "ins", "multi"], weights=[78, 8, 8, 6])[0]
            if kind == "snp":
                return ref, other(ref)
            if kind == "multi":
                return ref, ",".join(rng.sample([b for b in "ACGT" if b != ref], 2))
            if kind == "del":
                return ref_at(chrom, pos, rng.randint(2, 6)), ref
            return ref, ref + dna(rng.randint(1, 5))

        def columns(chrom, pos):
            ref, alt = site(chrom, pos)
            n_alt = alt.count(",") + 1
            qual = "." if rng.random() < 0.03 else f"{rng.uniform(3, 99):.1f}"
            dp = rng.randint(3, 80)
            flt = "q10" if qual != "." and float(qual) < 10 else "LowDP" if dp < 8 else "PASS"
            af = ",".join(f"{rng.uniform(0.001, 0.05) if rng.random() < 0.25 else rng.uniform(0.05, 0.95):.3f}" for _ in range(n_alt))
            return [chrom, str(pos), ".", ref, alt, qual, flt, f"DP={dp};AF={af}"], n_alt, dp

        def positions(n):
            spots, names = set(), list(genome)
            while len(spots) < n:
                chrom = rng.choices(names, weights=[len(genome[c]) for c in names])[0]
                spots.add((chrom, rng.randint(100, len(genome[chrom]) - 100)))
            return sorted(spots, key=lambda cp: (names.index(cp[0]), cp[1]))

        def gt(alleles):
            return "/".join(str(a) for a in sorted(alleles))

        trio = []
        for chrom, pos in positions(500):
            cols, n_alt, dp = columns(chrom, pos)
            pool = [0, 0, 1] if n_alt == 1 else [0, 1, 2]
            mother = [rng.choice(pool), rng.choice(pool)]
            father = [rng.choice(pool), rng.choice(pool)]
            child = [rng.choice(mother), rng.choice(father)]
            trio.append("\t".join(cols + ["GT:DP"] + [f"{gt(g)}:{rng.randint(2, dp)}" for g in (mother, father, child)]))
        write("cohort.vcf", head + columns_line + "\tFORMAT\tmother\tfather\tchild\n" + "\n".join(trio) + "\n")
        engine.convert_file(path("cohort.vcf"), path("cohort.vcf.gz"))  # sorted, BGZF: ready for an index

        for name, n in (("sample_a.vcf", 250), ("sample_b.vcf.gz", 200), ("sample_c.vcf", 150)):
            lines = []
            for chrom, pos in positions(n):
                cols, n_alt, _ = columns(chrom, pos)
                lines.append("\t".join(cols + ["GT", rng.choice(["0/1", "0/1", "1/1"] if n_alt == 1 else ["1/2", "0/1", "0/2"])]))
            rng.shuffle(lines)
            write(name, head + columns_line + "\tFORMAT\t" + name.split(".")[0] + "\n" + "\n".join(lines) + "\n")

        # Indels that are not left-aligned, a variant to trim and one to split.
        b = ref_at("chr2", 5_000, 3)
        indels = [
            ("chr1", 20_005, "homopolymer_del", ref_at("chr1", 20_005, 2), ref_at("chr1", 20_005)),
            ("chr1", 25_000, "snp", ref_at("chr1", 25_000), other(ref_at("chr1", 25_000))),
            ("chr1", 30_010, "repeat_del", ref_at("chr1", 30_010, 3), ref_at("chr1", 30_010)),
            ("chr2", 5_000, "trim_me", b, b[0] + other(b[1]) + b[2]),
            ("chr2", 12_008, "repeat_ins", ref_at("chr2", 12_008), ref_at("chr2", 12_008) + "AG"),
            ("chr10", 8_000, "multiallelic", ref_at("chr10", 8_000), ",".join(rng.sample([x for x in "ACGT" if x != ref_at("chr10", 8_000)], 2))),
        ]
        write("indels.vcf", head + columns_line + "\n" + "".join(f"{c}\t{p}\t{i}\t{r}\t{a}\t50\tPASS\tDP=30\n" for c, p, i, r, a in indels))

        def good(chrom, pos):
            return f"{chrom}\t{pos}\t.\t{ref_at(chrom, pos)}\t{other(ref_at(chrom, pos))}\t50\tPASS\tDP=20"

        problems = [
            good("chr1", 1_200), good("chr1", 3_400),
            good("chr1", 900),                                # out of order
            "2" + good("chr2", 1_500)[4:],                    # `2` among `chr2`
            "chr2\t2x00\t.\tA\tG\t50\tPASS\tDP=20",           # a position that is not a number
            "chr2\t3000\t.\tA\tG\t50",                        # six columns
            good("chr2", 4_000),
            "chr10\t500\t.\tAXG\tA\t50\tPASS\tDP=20",         # X is not a base
            good("chr10", 900),
        ]
        write("with_problems.vcf", head + columns_line + "\n" + "\n".join(problems) + "\n")
        write("five_columns.vcf", "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tQUAL\tGT\n"
              + "".join(f"chr1\t{p}\t.\t{q}\t{g}\n" for p, q, g in ((100, 10.5, "0/0"), (200, 50.2, "0/1"), (300, 33.0, "1/1"))))

        ensembl = []
        for chrom, pos in positions(30):
            cols, _, _ = columns(chrom, pos)
            ensembl.append("\t".join([chrom[3:]] + cols[1:]))
        ensembl += [f"MT\t{p}\t.\t{r}\t{a}\t60.0\tPASS\tDP=40" for p, r, a in ((73, "A", "G"), (263, "A", "G"), (16519, "T", "C"))]
        write("ensembl_names.vcf", head.replace(contigs, "".join(f"##contig=<ID={c[3:]},length={len(s)}>\n" for c, s in genome.items())
                                                + "##contig=<ID=MT,length=16569>\n") + columns_line + "\n" + "\n".join(ensembl) + "\n")

        # The lengths of the human chromosomes in the two assemblies.
        human = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
        grch38 = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
                  133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                  58617616, 64444167, 46709983, 50818468, 156040895, 57227415, 16569]
        grch37 = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431,
                  135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                  59128983, 63025520, 48129895, 51304566, 155270560, 59373566, 16569]

        def human_vcf(name, build, names, lengths):
            recs = []
            for _ in range(12):
                i = rng.randrange(3)
                ref = rng.choice("ACGT")
                recs.append((i, rng.randrange(10_000, 1_000_000), ref, other(ref)))
            recs.sort()
            write(name, f"##fileformat=VCFv4.2\n##reference={build}\n"
                  + "".join(f"##contig=<ID={n},length={l}>\n" for n, l in zip(names, lengths)) + columns_line + "\n"
                  + "".join(f"{names[i]}\t{p}\t.\t{r}\t{a}\t60\tPASS\t.\n" for i, p, r, a in recs))

        human_vcf("patient_grch38.vcf", "GRCh38", ["chr" + n if n != "MT" else "chrM" for n in human], grch38)
        human_vcf("cohort_grch37.vcf", "GRCh37", human, grch37)
        write("grch38_ensembl.sizes", "".join(f"{n}\t{l}\n" for n, l in zip(human, grch38)))

        # Genes and peaks (BED), and the same genes as a GFF3 annotation.
        genes = []
        for chrom, s in genome.items():
            x = rng.randint(500, 2_000)
            while x < len(s) - 7_000:
                length = rng.randint(1_500, 6_000)
                genes.append((chrom, x, x + length, rng.choice("+-")))
                x += length + rng.randint(500, 3_000)
        write("genes.bed", "".join(f"{c}\t{a}\t{e}\tgene{i + 1:02d}\t0\t{st}\n" for i, (c, a, e, st) in enumerate(genes)))
        peaks = []
        for _ in range(120):
            if rng.random() < 0.6:
                chrom, a, e, _st = rng.choice(genes)
                start = max(0, rng.randint(a - 300, e))
            else:
                chrom = rng.choice(list(genome))
                start = rng.randint(0, len(genome[chrom]) - 1_000)
            end = start + rng.randint(100, 600)
            peaks.append((chrom, start, end))
            if rng.random() < 0.25:  # a neighbour a few bases away, for merge -d
                gap = end + rng.randint(10, 150)
                peaks.append((chrom, gap, gap + rng.randint(100, 400)))
        peaks = [p for p in peaks if p[2] <= len(genome[p[0]])]
        rng.shuffle(peaks)
        write("peaks.bed", "".join(f"{c}\t{a}\t{e}\tpeak{i + 1:03d}\t{rng.randint(10, 1000)}\n" for i, (c, a, e) in enumerate(peaks)))
        gff = []
        exons_of = {}
        for i, (chrom, a, e, st) in enumerate(genes):
            gid = f"gene{i + 1:02d}"
            m = a + rng.randint(1, 50)  # the transcript starts a little after the gene
            gff.append(f"{chrom}\tMMGEN\tgene\t{a + 1}\t{e}\t.\t{st}\t.\tID={gid};Name={gid}")
            gff.append(f"{chrom}\tMMGEN\tmRNA\t{m + 1}\t{e}\t.\t{st}\t.\tID={gid}.t1;Parent={gid}")
            n_ex = rng.randint(2, 4)
            bounds = [m] + sorted(rng.sample(range(m + 100, e - 100), 2 * n_ex - 2)) + [e]
            for j in range(n_ex):
                gff.append(f"{chrom}\tMMGEN\texon\t{bounds[2 * j] + 1}\t{bounds[2 * j + 1]}\t.\t{st}\t.\tID={gid}.t1.exon{j + 1};Parent={gid}.t1")
                exons_of.setdefault(gid, []).append((chrom, bounds[2 * j], bounds[2 * j + 1]))
        rng.shuffle(gff)
        write("annotation.gff3", "##gff-version 3\n" + "".join(f"##sequence-region {c} 1 {len(s)}\n" for c, s in genome.items()) + "\n".join(gff) + "\n")

        # Alignments (SAM), out of order, and the same sorted as BAM.
        sam = ["@HD\tVN:1.6\tSO:unsorted"] + [f"@SQ\tSN:{c}\tLN:{len(s)}" for c, s in genome.items()] + ["@PG\tID:mmgen\tPN:mmgen_examples"]
        names = list(genome)
        for i in range(600):
            chrom = rng.choices(names, weights=[len(genome[c]) for c in names])[0]
            pos = rng.randint(1, len(genome[chrom]) - 100)
            seq, cigar = ref_at(chrom, pos, 100), "100M"
            if rng.random() < 0.1:  # the first 5 bases do not align
                seq, cigar, pos = dna(5) + seq[5:], "5S95M", pos + 5
            qual = "".join(chr(33 + q) for q in rng.choices(range(25, 41), k=100))
            sam.append(f"read{i + 1:04d}\t{rng.choice([0, 16])}\t{chrom}\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}\tNM:i:0")
        write("alignments.sam", "\n".join(sam) + "\n")
        engine.sort_genomic([path("alignments.sam")], path("alignments.sorted.bam"))

        # Reads piled on chr1:25,000 (the SNP of indels.vcf): about 40% carry
        # its ALT allele, as a heterozygous variant looks; a few have a base of
        # low quality there, which the pileup's quality filter leaves out.
        snp_ref = ref_at("chr1", 25_000)
        snp_alt = indels[1][4]
        pile = ["@HD\tVN:1.6\tSO:unsorted"] + [f"@SQ\tSN:{c}\tLN:{len(s)}" for c, s in genome.items()]
        for i in range(80):
            start = rng.randint(25_000 - 95, 25_000)
            seq = list(ref_at("chr1", start, 100))
            qual = [chr(33 + rng.randint(30, 40)) for _ in range(100)]
            k = 25_000 - start
            if rng.random() < 0.4:
                seq[k] = snp_alt
            if i % 16 == 0:
                qual[k] = "#"  # quality 2
            pile.append(f"pile{i + 1:03d}\t{rng.choice([0, 16])}\tchr1\t{start}\t60\t100M\t*\t0\t0\t{''.join(seq)}\t{''.join(qual)}")
        write("pileup_demo.sam", "\n".join(pile) + "\n")
        engine.sort_genomic([path("pileup_demo.sam")], path("pileup_demo.bam"), index=True)
        os.remove(path("pileup_demo.sam"))
        assert snp_ref != snp_alt

        # A child of first cousins (simulated): SNV calls every ~2 kb on two
        # chromosomes of 60 and 40 Mb, with no heterozygous call in chr1:20-32 Mb
        # and chr2:5-9 Mb, the stretches inherited twice from one ancestor.
        roh_lines = []
        for chrom, length, runs in (("chr1", 60_000_000, [(20_000_000, 32_000_000)]), ("chr2", 40_000_000, [(5_000_000, 9_000_000)])):
            pos = 5_000
            while pos < length - 5_000:
                in_run = any(a <= pos < b for a, b in runs)
                ref = rng.choice("ACGT")
                gt = "1/1" if in_run else rng.choices(["0/1", "1/1"], weights=[0.6, 0.4])[0]
                roh_lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{other(ref)}\t50\tPASS\t.\tGT\t{gt}")
                pos += rng.randint(1_000, 3_000)
        write("roh_demo.vcf", "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=60000000>\n##contig=<ID=chr2,length=40000000>\n"
              + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              + columns_line + "\tFORMAT\tchild\n" + "\n".join(roh_lines) + "\n")

        # Paired-end reads with the Illumina adapter in the short fragments, a
        # quality drop at the 3' end in some reads, and a few N.
        def quality_strings(length, count, drop):
            out = []
            for _ in range(count):
                tail = rng.randint(10, 30) if drop else 0
                out.append("".join(chr(33 + max(2, min(41, int((12 if i >= length - tail else 38 - 6 * i / length) + rng.uniform(-4, 3)))))
                                   for i in range(length)))
            return out

        good_q, bad_q = quality_strings(100, 60, False), quality_strings(100, 30, True)
        adapter1, adapter2 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        r1, r2 = [], []
        for i in range(3_000):
            chrom = rng.choice(names)
            size = rng.randint(50, 99) if rng.random() < 0.35 else rng.randint(100, 300)
            start = rng.randint(0, len(genome[chrom]) - size)
            frag = genome[chrom][start:start + size]
            for mate, seq, adapter in ((r1, frag, adapter1), (r2, revcomp(frag), adapter2)):
                read = list((seq + adapter + dna(100))[:100])
                qual = list(rng.choice(bad_q if rng.random() < 0.3 else good_q))
                if rng.random() < 0.02:
                    j = rng.randrange(100)
                    read[j], qual[j] = "N", "#"
                mate.append(f"@pair{i + 1:05d} {1 if mate is r1 else 2}:N:0:1\n{''.join(read)}\n+\n{''.join(qual)}\n")
        write("reads_R1.fastq.gz", "".join(r1))
        write("reads_R2.fastq.gz", "".join(r2))

        # A 30 kb "bacterium", a patient with 3 SNPs, reads at 30x with 0.5% of
        # errors, and their 31-mers as screen 1 writes them.
        bacterium = dna(30_000)
        write("toy_bacterium.fasta", fasta_text([("toy_bacterium", bacterium)]))
        patient = list(bacterium)
        for p in (5_000, 15_000, 25_000):
            patient[p] = {"A": "C", "C": "G", "G": "T", "T": "A"}[patient[p]]
        write("toy_patient.fasta", fasta_text([("toy_patient", "".join(patient))]))
        quals150 = quality_strings(150, 60, False)
        reads = []
        for i in range(6_000):
            start = rng.randint(0, len(bacterium) - 150)
            read = bacterium[start:start + 150]
            read = list(revcomp(read) if rng.random() < 0.5 else read)
            for j in range(150):
                if rng.random() < 0.005:
                    read[j] = other(read[j])
            reads.append(f"@read{i + 1:05d}\n{''.join(read)}\n+\n{rng.choice(quals150)}\n")
        write("toy_bacterium_reads.fastq.gz", "".join(reads))
        for name in ("toy_bacterium", "toy_patient"):
            pl.DataFrame({"K-mers": [engine.read_fasta_batch([path(name + ".fasta")], 31)]}).write_ipc(path(name + "_K31.arrow"))
        pl.DataFrame({"K-mers": [engine.read_fastq_batch([path("toy_bacterium_reads.fastq.gz")], 20, 31)]}).write_ipc(path("toy_bacterium_reads_K31.arrow"))

        # Screen 2: 31-mers of two genomes with different GC content.
        for name, gc in (("kmers_gc_rich.txt", 0.65), ("kmers_at_rich.txt", 0.35)):
            s = dna(2_030, gc)
            write(name, "\n".join(s[i:i + 31] for i in range(2_000)) + "\n")

        # Screen 7: a circular plasmid in 4 pieces and a linear chromosome in 5.
        plasmid = dna(20_000)
        cuts = [0, 4_800, 9_900, 15_200, 20_000]
        pieces = [plasmid[cuts[i]:cuts[i + 1]] + (plasmid + plasmid)[cuts[i + 1]:cuts[i + 1] + 54] for i in range(4)]
        rng.shuffle(pieces)
        write("plasmid_contigs.fasta", fasta_text((f"contig_{i + 1}", s) for i, s in enumerate(pieces)))
        linear = dna(24_000)
        cuts = [0, 5_000, 10_400, 15_100, 19_800, 24_000]
        pieces = [linear[cuts[i]:cuts[i + 1] + 30] for i in range(5)]
        rng.shuffle(pieces)
        write("chromosome_contigs.fasta", fasta_text((f"contig_{i + 1}", s) for i, s in enumerate(pieces)))

        # RNA-seq: four samples (two control, two treated) of reads from the
        # exons of the toy genome's genes, with three genes planted as changed
        # (8x up, 3x up, 8x down); counts with a little extra-Poisson spread.
        # A random generator of its own, so the files above do not change.
        rng2 = random.Random(7)
        rna_genes = sorted(exons_of)
        base_expr = {g: rng2.randint(60, 400) for g in rna_genes}
        changed = {rna_genes[2]: 8.0, rna_genes[6]: 3.0, rna_genes[9]: 0.125}
        rna_rows = []
        for sample, cond in (("ctrl_1", "control"), ("ctrl_2", "control"), ("treat_1", "treated"), ("treat_2", "treated")):
            depth = rng2.uniform(0.8, 1.25)
            lines = ["@HD\tVN:1.6\tSO:unsorted"] + [f"@SQ\tSN:{c}\tLN:{len(s)}" for c, s in genome.items()]
            for g in rna_genes:
                mean = base_expr[g] * depth * (changed.get(g, 1.0) if cond == "treated" else 1.0)
                n = max(0, int(round(rng2.gauss(mean, (mean + 0.02 * mean * mean) ** 0.5))))
                usable = [x for x in exons_of[g] if x[2] - x[1] >= 120]
                for k in range(n if usable else 0):
                    c, a, e = rng2.choice(usable)
                    pos = rng2.randint(a + 1, e - 99)
                    lines.append(f"{sample}.{g}.{k}\t{rng2.choice([0, 16])}\t{c}\t{pos}\t60\t100M\t*\t0\t0\t{ref_at(c, pos, 100)}\t{'I' * 100}\tNH:i:1")
            write(f"rnaseq_{sample}.sam", "\n".join(lines) + "\n")
            engine.sort_genomic([path(f"rnaseq_{sample}.sam")], path(f"rnaseq_{sample}.bam"))
            os.remove(path(f"rnaseq_{sample}.sam"))
            rna_rows.append(f"{sample}\t{cond}")
        write("rnaseq_samples.tsv", "sample\tcondition\n" + "\n".join(rna_rows) + "\n")
        engine.count_reads([path(f"rnaseq_{s.split(chr(9))[0]}.bam") for s in rna_rows], path("annotation.gff3"),
                        path("rnaseq_counts.tsv"), sample_names=[s.split("\t")[0] for s in rna_rows])

        # Three "species" for screen 2: 100 kb each from its own third-order
        # Markov model (each 3-base context has its own next-base
        # probabilities), at about the same GC, so only the 4-base word
        # preferences tell them apart. A generator of its own, so the files
        # above do not change.
        rng3 = random.Random(13)
        species = []
        for sp in ("species_A", "species_B", "species_C"):
            table = {}
            for ctx in range(64):
                w = [rng3.gammavariate(0.6, 1.0) for _ in range(4)]
                tot = sum(w)
                table[ctx] = [x / tot for x in w]
            ctx, out = 0, []
            for _ in range(100_000):
                r, acc, b = rng3.random(), 0.0, 3
                for i, pr in enumerate(table[ctx]):
                    acc += pr
                    if r < acc:
                        b = i
                        break
                out.append("ACGT"[b])
                ctx = ((ctx << 2) | b) & 63
            species.append((sp, "".join(out)))
        write("three_species.fasta", fasta_text(species))

        write("README.md", "# MMGEN examples\n\nSmall simulated files, written by the web interface (`ui_main.py`) when an example "
              "is picked on a screen. They are made up: no real person or organism. Delete the folder at any time; "
              "it is written again when needed.\n\n"
              "* `toy_genome.fa` (+ `.fai`): a genome of 100,000 bases in three chromosomes (chr1, chr2, chr10), with three repeats.\n"
              "* `cohort.vcf` (+ `cohort.vcf.gz`, sorted BGZF): 500 variants of a trio (mother, father, child) on that genome.\n"
              "* `sample_a.vcf`, `sample_b.vcf.gz`, `sample_c.vcf`: three single-sample VCF files, out of order.\n"
              "* `indels.vcf`: indels that are not left-aligned, a variant to trim and one to split.\n"
              "* `with_problems.vcf`, `five_columns.vcf`: files with planted defects, for the validator.\n"
              "* `ensembl_names.vcf`: chromosomes named 1, 2, 10, MT.\n"
              "* `patient_grch38.vcf`, `cohort_grch37.vcf`, `grch38_ensembl.sizes`: headers of the two human assemblies.\n"
              "* `genes.bed`, `peaks.bed`, `annotation.gff3`: intervals and an annotation on the toy genome.\n"
              "* `alignments.sam`, `alignments.sorted.bam`: 600 aligned reads, out of order and sorted.\n"
              "* `pileup_demo.bam` (+ `.bai`): 80 reads piled on chr1:25,000, about 40% with the ALT allele of the SNP there.\n"
              "* `roh_demo.vcf`: a simulated child of first cousins, with two long runs of homozygosity.\n"
              "* `rnaseq_*.bam`, `rnaseq_samples.tsv`, `rnaseq_counts.tsv`: four simulated RNA-seq samples (two control, two "
              "treated) on the toy genome's genes, with three genes changed by treatment.\n"
              "* `human/`: the human-genome demo data (public data cut from Genome in a Bottle, ClinVar, RefSeq and the "
              "PGS Catalog), unpacked from `mmgen_demo_human_v1.zip` when a human-genome example is picked.\n"
              "* `reads_R1.fastq.gz`, `reads_R2.fastq.gz`: 3,000 read pairs with an adapter and low-quality ends.\n"
              "* `toy_bacterium.fasta`, `toy_patient.fasta` (3 SNPs), `toy_bacterium_reads.fastq.gz` (30x, 0.5% errors) "
              "and their 31-mers (`*_K31.arrow`).\n"
              "* `kmers_gc_rich.txt`, `kmers_at_rich.txt`: 31-mers of two genomes with different GC content.\n"
              "* `three_species.fasta`: three simulated 100 kb genomes with different 4-base word preferences (screen 2).\n"
              "* `plasmid_contigs.fasta`, `chromosome_contigs.fasta`: contigs that overlap at their ends.\n")
        with open(stamp, "w", encoding="utf-8") as fh:
            fh.write(EXAMPLES_VERSION)
        return True

    _E = EXAMPLES_DIR + "/"
    _TRIO = "`cohort.vcf` holds 500 variants of a trio (mother, father and child) on a small simulated genome: SNPs, insertions, deletions and multiallelic sites."
    EXAMPLES = {
        ("qc", None): [
            ("A 30 kb bacterial genome (FASTA), k = 31",
             "Every 31-base window of a simulated genome of 30,000 bases becomes a k-mer: expect 29,970 distinct k-mers in `toy_bacterium_K31.arrow`, the file that screens 3, 4 and 6 read.",
             {"files": [_E + "toy_bacterium.fasta"], "kmer_size": 31}),
            ("Sequencing reads at 30× coverage (FASTQ.gz), mean quality ≥ 20",
             "6,000 reads of 150 bases from the same genome, from both strands, with 0.5% of sequencing errors. Each error makes k-mers that are not in the genome, and the reverse strand doubles the rest: expect several times the 29,970 k-mers of the genome.",
             {"files": [_E + "toy_bacterium_reads.fastq.gz"], "kmer_size": 31, "min_quality": 20}),
            ("A patient with 3 mutations (FASTA), k = 31",
             "The same genome with three single-base changes. Its k-mers, `toy_patient_K31.arrow`, are compared with the healthy genome's on screen 4.",
             {"files": [_E + "toy_patient.fasta"], "kmer_size": 31}),
        ],
        ("cluster", None): [
            ("Three simulated species in one file",
             "Three 100 kb genomes, each drawn from its own composition model (the same GC content, different "
             "4-base word preferences), cut into 5 kb windows. Expect 3 groups, each holding the 20 windows of one "
             "species: the tetranucleotide signature tells them apart even at equal GC.",
             {"fastas": [_E + "three_species.fasta"], "window": 5000, "min_cluster": 5, "dims": 5}),
        ],
        ("assembly", None): [
            ("The 31-mers of the 30 kb genome",
             "A genome without repeats is rebuilt as a single contig of 30,000 bases, written to `MMGEN_Contigs_K31.fasta` next to the input.",
             {"kmer_file": _E + "toy_bacterium_K31.arrow"}),
            ("The 31-mers of the reads, sequencing errors included",
             "Every sequencing error opens a branch in the graph, and the reverse strand makes a second copy of it, so the assembly breaks into many short contigs. This module does not correct errors: assemble error-free sequences, or look at the errors first with the k-mer spectrum (screen 13).",
             {"kmer_file": _E + "toy_bacterium_reads_K31.arrow"}),
        ],
        ("variant", None): [
            ("A patient with 3 SNPs against the healthy genome",
             "Each SNP changes the 31 k-mers that cover it: expect 3 × 31 = 93 k-mers found only in the patient.",
             {"reference": _E + "toy_bacterium_K31.arrow", "patient": _E + "toy_patient_K31.arrow", "filter_mb": 16, "num_hashes": 5}),
        ],
        ("annotation", None): [
            ("A short gene: start codon, 14 codons, stop codon",
             "Translates to MSTHEGENEISREAD followed by the stop codon (the red asterisk).",
             {"dna": "ATGTCTACTCATGAAGGTGAAAATGAAATTTCTCGAGAAGCTGATTAA"}),
            ("The same gene with a nonsense mutation",
             "One base changed (C → T) turns the arginine codon CGA into the stop codon TGA: the protein ends after 11 amino acids.",
             {"dna": "ATGTCTACTCATGAAGGTGAAAATGAAATTTCTTGAGAAGCTGATTAA"}),
            ("The same gene with one base deleted (a frameshift)",
             "Without the C of the fourth codon, every codon after it is read out of frame: the protein changes from there on, and the stop codon is lost.",
             {"dna": "ATGTCTACTATGAAGGTGAAAATGAAATTTCTCGAGAAGCTGATTAA"}),
            ("A mutated stretch of E. coli (a result of screen 4)",
             "One of the k-mers that screen 4 reported in a simulated patient made from E. coli.",
             {"dna": "ATTTACCCTTGCTTTCACATAATCTGCCCGCCGATTTTGCTGCGTTGCGTAAATA"}),
        ],
        ("bloom", None): [
            ("The 31-mers of the 30 kb genome, 1 MB filter, 5 hashes",
             "29,970 k-mers × 5 hashes set about 150 thousand of the 8.4 million bits of the filter: a saturation under 2%, so false positives are rare.",
             {"kmer_file": _E + "toy_bacterium_K31.arrow", "filter_mb": 1, "num_hashes": 5}),
            ("The 31-mers of the reads (many more), same filter",
             "The sequencing errors multiply the k-mers: the same filter fills up much more, and the chance of a false positive grows with it.",
             {"kmer_file": _E + "toy_bacterium_reads_K31.arrow", "filter_mb": 1, "num_hashes": 5}),
        ],
        ("scaffold", None): [
            ("A circular plasmid in 4 contigs that overlap by 54 bases",
             "The four pieces of a plasmid of 20,000 bases, out of order: each overlaps the next by 54 bases, and the last overlaps the first. Expect one scaffold of 20,000 bases, recognized as circular.",
             {"contigs_fasta": _E + "plasmid_contigs.fasta", "overlap": 54}),
            ("A linear chromosome in 5 contigs that overlap by 30 bases",
             "Five pieces of a sequence of 24,000 bases: expect one scaffold of 24,000 bases, not circular.",
             {"contigs_fasta": _E + "chromosome_contigs.fasta", "overlap": 30}),
        ],
        ("vcf", None): [
            ("Three samples (one gzipped), merged by position",
             "Three VCF files out of order, 600 variants in all, become one file sorted by chromosome (in lexical order: chr1, chr10, chr2) and position. Variants at the same position keep the order of the input files: the sort is stable.",
             {"files": [_E + "sample_a.vcf", _E + "sample_b.vcf.gz", _E + "sample_c.vcf"], "sort_by": "POS", "output_name": _E + "merged_by_position.vcf"}),
            ("The same three, by quality (highest first)",
             "The same 600 variants ordered by QUAL, from the highest; a missing QUAL (`.`) counts as 0.",
             {"files": [_E + "sample_a.vcf", _E + "sample_b.vcf.gz", _E + "sample_c.vcf"], "sort_by": "QUAL", "output_name": _E + "merged_by_quality.vcf"}),
        ],
        ("sort", "Sort & merge"): [
            ("Two VCF samples out of order → one sorted .vcf.gz with its index",
             "The output keeps the chromosome order of the `##contig` lines (chr1, chr2, chr10), as bcftools does, and is written as BGZF with a `.tbi` index, ready for **Index & query a region**.",
             {"files": [_E + "sample_a.vcf", _E + "sample_b.vcf.gz"], "index": True}),
            ("Peaks (BED) in natural chromosome order",
             "BED files are sorted in lexical order by default (chr1, chr10, chr2), as bedtools does; the natural order puts chr2 before chr10.",
             {"files": [_E + "peaks.bed"], "order": "natural"}),
            ("A GFF3 annotation: each gene before its transcript and exons",
             "The lines of `annotation.gff3` come shuffled. Sorted by position, a feature that starts at the same base as another but is longer comes first, so a parent precedes its children.",
             {"files": [_E + "annotation.gff3"]}),
            ("SAM alignments → a sorted BAM with its .bai index",
             "600 reads out of order. With the index asked for, the output is a BAM; its `@HD` line now says `SO:coordinate`.",
             {"files": [_E + "alignments.sam"], "index": True}),
        ],
        ("sort", "Index & query a region"): [
            ("One region of a VCF (the .tbi is built first)",
             "`cohort.vcf.gz` is sorted and compressed as BGZF. The first run writes `cohort.vcf.gz.tbi`; after that, only the blocks that hold the region are read.",
             {"file": _E + "cohort.vcf.gz", "regions": "chr1:10000-20000"}),
            ("The reads of a BAM that touch a region",
             "The `.bai` index is built first; then the reads of chr2 between 5,000 and 6,000 are shown.",
             {"file": _E + "alignments.sorted.bam", "regions": "chr2:5000-6000"}),
            ("Several regions, saved to a file",
             "A whole chromosome (chr10) and the first 10,000 bases of chr2, written to `cohort.regions.vcf` with the header.",
             {"file": _E + "cohort.vcf.gz", "regions": "chr2:1-10000, chr10", "output": _E + "cohort.regions.vcf"}),
        ],
        ("sort", "merge"): [
            ("Join peaks that overlap or are up to 100 bases apart",
             "`peaks.bed` holds about 150 peaks, some of them pairs a few bases apart: with a distance of 100, those pairs become one interval.",
             {"a": [_E + "peaks.bed"], "distance": 100}),
        ],
        ("sort", "intersect"): [
            ("The peaks that touch a gene (-u)",
             "Each peak that overlaps at least one gene of `genes.bed` is written once, as it is.",
             {"a": [_E + "peaks.bed"], "b": [_E + "genes.bed"], "mode": "unique"}),
        ],
        ("sort", "subtract"): [
            ("The parts of the genes that no peak covers",
             "Each gene is clipped where a peak overlaps it; a gene with a peak in the middle becomes two pieces.",
             {"a": [_E + "genes.bed"], "b": [_E + "peaks.bed"]}),
        ],
        ("sort", "complement"): [
            ("The stretches of the genome outside every gene",
             "The intergenic regions of the toy genome, from the chromosome lengths in `toy_genome.fa.fai`.",
             {"a": [_E + "genes.bed"], "genome": _E + "toy_genome.fa.fai"}),
        ],
        ("sort", "closest"): [
            ("The closest gene to each peak, with the distance",
             "A peak inside a gene is at distance 0; a peak between two genes gets the nearer one.",
             {"a": [_E + "peaks.bed"], "b": [_E + "genes.bed"]}),
        ],
        ("sort", "coverage"): [
            ("How much of each gene the peaks cover",
             "For each gene: how many peaks overlap it, how many of its bases they cover, its length and the covered fraction.",
             {"a": [_E + "genes.bed"], "b": [_E + "peaks.bed"]}),
        ],
        ("sort", "slop"): [
            ("Each gene with 1,000 bases added on both sides",
             "A common way to include the promoter and the downstream region; intervals stop at the chromosome ends given by `toy_genome.fa.fai`.",
             {"a": [_E + "genes.bed"], "genome": _E + "toy_genome.fa.fai", "left": 1000, "right": 1000}),
        ],
        ("validate", "Validate a file"): [
            ("A VCF with planted problems",
             "`with_problems.vcf` has a position out of order, `2` among `chr2`, a position that is not a number, a line with too few columns and an invalid base. Each kind of problem is reported once, at the line where it first shows up, with how to fix it.",
             {"file": _E + "with_problems.vcf"}),
            ("A file that is not a VCF: only 5 columns",
             "`five_columns.vcf` has the columns CHROM, POS, ID, QUAL and GT. A VCF needs at least eight: CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO.",
             {"file": _E + "five_columns.vcf"}),
            ("A clean multi-sample VCF",
             _TRIO + " Expect no problem.",
             {"file": _E + "cohort.vcf"}),
            ("Sequencing reads (FASTQ.gz)",
             "3,000 reads of the first file of a pair: the validator checks the four lines of every read.",
             {"file": _E + "reads_R1.fastq.gz"}),
        ],
        ("validate", "Rename chromosomes (chr1 ↔ 1)"): [
            ("Ensembl names (1, 2, 10, MT) → UCSC (chr1, chr2, chr10, chrM)",
             "The records and the `##contig` lines are renamed; MT becomes chrM.",
             {"file": _E + "ensembl_names.vcf", "style": "chr"}),
            ("UCSC names → Ensembl, in a BED file",
             "chr1, chr2 and chr10 become 1, 2 and 10.",
             {"file": _E + "peaks.bed", "style": "no_chr"}),
        ],
        ("validate", "Compare reference genomes"): [
            ("GRCh38 against GRCh37: different genomes",
             "The two headers declare the human chromosomes with their lengths in GRCh38 and in GRCh37. The lengths differ, so coordinates from one file mean other bases in the other.",
             {"a": _E + "patient_grch38.vcf", "b": _E + "cohort_grch37.vcf"}),
            ("The same genome (GRCh38), named two ways",
             "A VCF with `chr1` against a chromosome-sizes file with `1`: the lengths match, so it is the same genome; only the names must be made equal before combining them.",
             {"a": _E + "patient_grch38.vcf", "b": _E + "grch38_ensembl.sizes"}),
        ],
        ("vcftools", "Statistics"): [
            ("A trio: mother, father and child", _TRIO, {"file": _E + "cohort.vcf"}),
            ("A file that is not a VCF (5 columns)",
             "Nothing can be counted: the screen says why, with the validator's report.",
             {"file": _E + "five_columns.vcf"}),
        ],
        ("vcftools", "Filter"): [
            ("High-quality SNPs that passed every filter",
             _TRIO + " Keeps the sites whose ALT alleles are all SNPs, with QUAL of at least 30 and FILTER PASS.",
             {"file": _E + "cohort.vcf", "expression": 'TYPE == "snp" && QUAL >= 30 && FILTER == "PASS"'}),
            ("Rare variants (allele frequency under 5%)",
             "Keeps the sites whose `INFO/AF` is below 0.05.",
             {"file": _E + "cohort.vcf", "expression": "INFO/AF < 0.05"}),
            ("One chromosome, and only the child",
             "The records of chr2, with the genotype column of the child alone.",
             {"file": _E + "cohort.vcf", "regions": "chr2", "samples": "child"}),
        ],
        ("vcftools", "Normalize"): [
            ("Split the multiallelic sites of the trio",
             "Each site with two ALT alleles becomes two records, one per allele, with the genotypes rewritten for each.",
             {"file": _E + "cohort.vcf", "split": True}),
            ("Left-align indels in repeats, with the genome",
             "`indels.vcf` has a deletion in a run of T, a deletion and an insertion in dinucleotide repeats, all written at the right end of the repeat. With `toy_genome.fa`, each moves to the leftmost position where it can be written, as `bcftools norm -f` does; a variant with extra bases is trimmed and a multiallelic site is split.",
             {"file": _E + "indels.vcf", "split": True, "reference": _E + "toy_genome.fa"}),
        ],
        ("fastq", "Quality report"): [
            ("Paired-end reads with an adapter and low-quality ends (R1)",
             "3,000 reads of 100 bases. About a third come from fragments shorter than the read, so the Illumina adapter shows at their 3' end; some reads lose quality at the end, and a few have an N.",
             {"file": _E + "reads_R1.fastq.gz"}),
            ("Reads at 30× from a bacterial genome",
             "6,000 reads of 150 bases with 0.5% of errors and no adapter.",
             {"file": _E + "toy_bacterium_reads.fastq.gz"}),
        ],
        ("fastq", "Trim reads"): [
            ("Single-end: the adapter found automatically, quality 20",
             "The adapter is detected in the first reads; the ends below quality 20 and the adapter are cut, and the reads that end up shorter than 20 bases are dropped.",
             {"file": _E + "reads_R1.fastq.gz", "adapter": "auto", "quality": 20, "min_length": 20}),
            ("Paired-end: R1 and R2 trimmed in step",
             "When one read of a pair is dropped, so is its mate: the two outputs keep the same reads in the same order.",
             {"file": _E + "reads_R1.fastq.gz", "file2": _E + "reads_R2.fastq.gz", "adapter": "auto", "quality": 20, "min_length": 20}),
            ("Strict: quality 30, N removed, reads under 50 bases dropped",
             "Harder settings cut more and drop more reads: compare the counts with the first example.",
             {"file": _E + "reads_R1.fastq.gz", "adapter": "AGATCGGAAGAGC", "quality": 30, "min_length": 50, "trim_n": True}),
        ],
        ("kmercount", None): [
            ("Reads at 30× coverage (k = 31): the error peak and the coverage peak",
             "Most k-mers seen once are sequencing errors; the k-mers of the genome form a second peak near 20×. A read of 150 bases holds 120 k-mers of 31, so a coverage of 30× in bases is 30 × 120 / 150 = 24× in k-mers, and the k-mers that hold an error bring it down a little more.",
             {"files": [_E + "toy_bacterium_reads.fastq.gz"], "k": 31}),
            ("The genome itself (k = 21): almost every k-mer seen once",
             "A random genome of 30,000 bases has almost no repeated 21-mer.",
             {"files": [_E + "toy_bacterium.fasta"], "k": 21}),
            ("The reads, with the table of the k-mers seen at least 5 times",
             "Leaving out the k-mers seen fewer than 5 times drops most of the errors from the table.",
             {"files": [_E + "toy_bacterium_reads.fastq.gz"], "k": 31, "min_count": 5, "output": _E + "reads_31mers.tsv.gz"}),
        ],
        ("convert", None): [
            ("VCF → Parquet, QUAL as numbers",
             _TRIO + " As Parquet it opens as a table in pandas, polars or DuckDB.",
             {"file": _E + "cohort.vcf", "target": ".parquet", "numbers": "number"}),
            ("VCF → Parquet, QUAL kept as text",
             "QUAL is stored exactly as written, so converting back gives the same bytes.",
             {"file": _E + "cohort.vcf", "target": ".parquet", "numbers": "text", "output": _E + "cohort.text.vcf.parquet"}),
            ("SAM → BAM", "The 600 alignments of `alignments.sam` in the binary format.",
             {"file": _E + "alignments.sam", "target": ".bam"}),
            ("BED → gzip (BGZF)", "The peaks compressed as BGZF, which any gzip reader opens and which can be indexed.",
             {"file": _E + "peaks.bed", "target": ".gz"}),
        ],
    }
    # The human-genome screens use the demo data in mmgen_examples/human: real,
    # public genomes of the Genome in a Bottle project, cut to chr22, the genes
    # of the example panels and the PharmCAT positions.
    _H = HUMAN_DIR + "/"
    _P = EXAMPLES_DIR + "/pgx/"
    _R = RESEARCH_DIR + "/"
    _R16 = EXAMPLES_DIR + "/research_16s/"
    _GIAB = ("Genome in a Bottle (GIAB) reference genomes are public data from consenting volunteers of the Personal Genome "
             "Project, used worldwide to test sequencing. The demo keeps chr22, the example panel genes and the "
             "pharmacogene positions.")
    EXAMPLES[("cluster", None)].append(
        ("Four bacterial isolates: two Klebsiella pneumoniae, two Staphylococcus aureus",
         "The assemblies of four isolates of screen 25, in 10 kb windows. Expect the two species apart (about 57% and "
         "33% GC) and the two strains of each species together: composition separates species, not strains. Plasmid "
         "windows may form small groups of their own or stay ungrouped.",
         {"fastas": [_R + f"microbio/isolates/{_g}.fna.gz" for _g in ("kpneumoniae_KPNIH1_KPC-3", "kpneumoniae_ATCC_BAA-2146_NDM-1",
                                                                        "saureus_USA300_MRSA", "saureus_NCTC8325_MSSA")],
          "window": 10000, "min_cluster": 10, "dims": 5}))
    EXAMPLES.update({
        ("findings", "ClinVar findings"): [
            ("HG002: pathogenic and likely pathogenic, at least 1 review star",
             _GIAB + " HG002 is the son of the Ashkenazi trio. Expect one finding in one copy: a variant ClinVar classifies as "
             "pathogenic, with 2 stars, in F11 (coagulation factor XI). One copy of a variant for a recessive condition is a "
             "carrier state; most people carry a few.",
             {"vcf": _H + "giab/HG002.vcf.gz", "callable": _H + "giab/HG002.callable.bed", "clinvar": _H + "clinvar_demo.vcf.gz",
              "classes": ["plp"], "min_stars": 1}),
            ("HG002, a hereditary-cancer panel: every classification",
             "The same person, only the genes of the hereditary breast and ovarian cancer panel, and every ClinVar class with any "
             "number of stars: most matches are benign; a few are uncertain or conflicting.",
             {"vcf": _H + "giab/HG002.vcf.gz", "callable": _H + "giab/HG002.callable.bed", "clinvar": _H + "clinvar_demo.vcf.gz",
              "classes": ["plp", "vus", "conflicting", "blb", "risk"], "min_stars": 0, "panel": "Hereditary breast and ovarian cancer"}),
            ("HG005: pathogenic and likely pathogenic, at least 1 review star",
             _GIAB + " HG005 is the son of the Han Chinese trio.",
             {"vcf": _H + "giab/HG005.vcf.gz", "callable": _H + "giab/HG005.callable.bed", "clinvar": _H + "clinvar_demo.vcf.gz",
              "classes": ["plp"], "min_stars": 1}),
        ],
        ("findings", "Gene panel from an annotation"): [
            ("Hereditary breast and ovarian cancer: the gene coordinates as BED",
             "The 13 genes of the panel, found in the RefSeq GRCh38 annotation (its accessions NC_0000.. are turned into chr names).",
             {"gff": _H + "panel_genes.gff3.gz", "panel": "Hereditary breast and ovarian cancer", "padding": 0}),
            ("Pharmacogenes, with 1,000 bases on each side",
             "The genes that drug-response tests look at, padded to include their promoters.",
             {"gff": _H + "panel_genes.gff3.gz", "panel": "Pharmacogenes", "padding": 1000}),
        ],
        ("pgx", "PharmCAT"): [
            ("HG002 (son of the Ashkenazi trio)",
             _GIAB + " PharmCAT reads the genotype at each of its 1,226 positions; the X chromosome is not in this release of "
             "the GIAB calls, so its genes stay without a call.",
             {"vcf": _H + "giab/HG002.vcf.gz", "callable": _H + "giab/HG002.callable.bed"}),
            ("HG005 (son of the Han Chinese trio)", _GIAB,
             {"vcf": _H + "giab/HG005.vcf.gz", "callable": _H + "giab/HG005.callable.bed"}),
            ("HG001 (NA12878)", _GIAB + " NA12878 is one of the most studied human genomes.",
             {"vcf": _H + "giab/HG001.vcf.gz", "callable": _H + "giab/HG001.callable.bed"}),
        ],
        ("pgx", "MMGEN's engine (mmgen-pgx)"): [
            ("HG002 (son of the Ashkenazi trio)",
             _GIAB + " MMGEN's engine reads its 19 genes from the VCF; chrX (G6PD) and chrM (MT-RNR1) are not in this release "
             "of the GIAB calls, so those two stay without a call, with the reason. Run the PharmCAT engine on it too and compare.",
             {"vcf": _H + "giab/HG002.vcf.gz", "callable": _H + "giab/HG002.callable.bed"}),
            ("HG002 with its reads: CYP2D6 and the HLA genes",
             _GIAB + " The reads of HG002 over CYP2D6 and the HLA class I genes (a cut of its public 2x250 Illumina BAM). "
             "Expect CYP2D6 \\*2/\\*4 with two copies, the HLA-A, -B and -C types, and an outside-call file that PharmCAT can "
             "use. It needs the HLA reference sequences (the download button above).",
             {"vcf": _H + "giab/HG002.vcf.gz", "callable": _H + "giab/HG002.callable.bed",
              "bam": _P + "HG002.GRCh38.2x250.pgx_v2.bam"}),
            ("HG001 (NA12878)",
             _GIAB + " NA12878 is one of the most studied human genomes. Its VCF is not phased, so in CYP2C19 the positions read "
             "fit both \\*1/\\*2 and \\*2/\\*35. The engine calls the highest-scoring diplotype, as PharmCAT does: \\*1/\\*2 "
             "(also the GeT-RM consensus for this sample), and lists \\*2/\\*35 (poor metabolizer) as *also possible*, with a "
             "warning. Compare the two engines here.",
             {"vcf": _H + "giab/HG001.vcf.gz", "callable": _H + "giab/HG001.callable.bed"}),
        ],
        ("pgx", "A group of people and a protocol's medicines (mmgen-pgx)"): [
            ("Seven people and the medicines given around a spinal cord surgery",
             _GIAB + " The seven Genome in a Bottle people, as if they were the participants of a clinical study, against an "
             "illustrative list of the medicines given around a spinal cord surgery: analgesics, anti-inflammatories, "
             "antiemetics, anticoagulants, stomach protection and anaesthetics. Expect: clopidogrel and the proton-pump "
             "inhibitors guided by CYP2C19, the anti-inflammatories by CYP2C9, warfarin by CYP2C9, VKORC1 and CYP4F2, the "
             "volatile anaesthetics and succinylcholine by RYR1 and CACNA1S; codeine, tramadol, ondansetron and amitriptyline "
             "need CYP2D6, which a VCF alone does not give (*gene not called*); gabapentin, pregabalin, enoxaparin and morphine "
             "have no CPIC guideline. Each person's ancestry is estimated with the 1000 Genomes panel and weights the population "
             "frequencies NAT2 is chosen by.",
             {"vcfs": [_H + f"giab/HG00{_i}.vcf.gz" for _i in range(1, 8)], "callable": "auto",
              "ancestry_panel": _H + "ancestry_panel_genomewide_1000g.tsv.gz",
              "drugs": "codeine, tramadol, morphine, ondansetron, celecoxib, ibuprofen, amitriptyline, gabapentin, pregabalin, "
                       "warfarin, clopidogrel, enoxaparin, omeprazole, pantoprazole, sevoflurane, succinylcholine"}),
        ],
        ("lookup", "Look up variants in relatives (cascade testing)"): [
            ("The F11 variant of HG002 in his parents",
             _GIAB + " HG002 carries one copy of a pathogenic F11 variant (screen 18). Which parent passed it on?",
             {"variants": "chr4:186274193:G:T", "vcfs": [_H + "giab/HG002.vcf.gz", _H + "giab/HG003.vcf.gz", _H + "giab/HG004.vcf.gz"],
              "use_beds": True}),
            ("Four pharmacogene variants in the Han Chinese trio",
             "CYP2C19*2, CYP2C19*17, CYP3A5*3 and the VKORC1 promoter variant (GRCh38 positions from dbSNP) in HG005, HG006 and HG007.",
             {"variants": "chr10:94781859:G:A\nchr10:94761900:C:T\nchr7:99672916:T:C\nchr16:31096368:C:T",
              "vcfs": [_H + "giab/HG005.vcf.gz", _H + "giab/HG006.vcf.gz", _H + "giab/HG007.vcf.gz"], "use_beds": True}),
            ("By rsID, through ClinVar",
             "rs121965063 is resolved to its coordinates with the ClinVar file, then looked up in the trio.",
             {"variants": "rs121965063", "vcfs": [_H + "giab/HG002.vcf.gz", _H + "giab/HG003.vcf.gz", _H + "giab/HG004.vcf.gz"],
              "use_beds": True, "clinvar": _H + "clinvar_demo.vcf.gz"}),
        ],
        ("kin", "Relatedness of several people"): [
            ("The seven GIAB genomes: two trios and NA12878",
             _GIAB + " Expect four parent-offspring pairs (HG002 with HG003 and HG004, HG005 with HG006 and HG007) and every "
             "other pair unrelated; the two couples are not related to each other. Between people of different ancestries "
             "the kinship comes out negative, a known effect of the KING estimator.",
             {"vcfs": [_H + f"giab/HG00{_i}.vcf.gz" for _i in range(1, 8)], "use_beds": True}),
        ],
        ("kin", "Trio and paternity test"): [
            ("HG002 with his parents (the Ashkenazi trio)",
             _GIAB + " The father is not excluded: the few opposite homozygotes between him and the son are at the level of "
             "genotyping errors, and the paternity index is very high. Expect a few de novo candidates.",
             {"child": _H + "giab/HG002.vcf.gz", "father": _H + "giab/HG003.vcf.gz", "mother": _H + "giab/HG004.vcf.gz",
              "use_beds": True, "af": _H + "af_1000g_demo.vcf.gz", "af_field": "EUR_AF"}),
            ("A wrong father: HG006 tested as the father of HG002",
             "The same son and mother, with a man of the other trio as the alleged father: opposite homozygotes at a rate "
             "far above the error level, so he is excluded.",
             {"child": _H + "giab/HG002.vcf.gz", "father": _H + "giab/HG006.vcf.gz", "mother": _H + "giab/HG004.vcf.gz",
              "use_beds": True, "af": _H + "af_1000g_demo.vcf.gz", "af_field": "AF"}),
            ("A duo: HG005 and his father, no mother",
             "Without the mother, fewer sites are informative and the index is smaller, but a true father is still not excluded.",
             {"child": _H + "giab/HG005.vcf.gz", "father": _H + "giab/HG006.vcf.gz", "mother": "none",
              "use_beds": True, "af": _H + "af_1000g_demo.vcf.gz", "af_field": "EAS_AF"}),
        ],
        ("anc", "Ancestry"): [
            ("The seven GIAB genomes against the 1000 Genomes groups (chr22)",
             _GIAB + " The panel holds the allele frequencies of the five 1000 Genomes super-populations on chr22. Expect "
             "HG001 European, the Han Chinese trio East Asian, and the Ashkenazi trio mostly European with a South Asian share: "
             "Ashkenazi and Middle Eastern ancestry fall between the panel's groups. One chromosome gives rough errors.",
             {"panel": _H + "ancestry_panel_chr22_1000g.tsv.gz", "vcfs": [_H + f"giab/HG00{_i}.vcf.gz" for _i in range(1, 8)], "use_beds": True}),
        ],
        ("anc", "Polygenic score"): [
            ("A breast-cancer score (PGS000004, 313 variants) for HG002",
             _GIAB + " The demo data keeps the score's variants. The result is where HG002's weighted sum falls among the panel's "
             "groups — a statistical summary that transfers poorly across ancestries, not a risk for a person.",
             {"scoring": _H + "pgs/PGS000004_hmPOS_GRCh38.txt.gz", "vcfs": [_H + "giab/HG002.vcf.gz"], "use_beds": True,
              "panel": _H + "pgs_panel_1000g.tsv.gz"}),
            ("A coronary-artery-disease score (PGS003438) for the seven GIAB genomes",
             "The same score for seven people: compare each one with the group of their own ancestry.",
             {"scoring": _H + "pgs/PGS003438_hmPOS_GRCh38.txt.gz", "vcfs": [_H + f"giab/HG00{_i}.vcf.gz" for _i in range(1, 8)],
              "use_beds": True, "panel": _H + "pgs_panel_1000g.tsv.gz"}),
        ],
        ("anc", "Build a reference panel"): [
            ("From the frequencies of a sites VCF (1000 Genomes)",
             "The 1000 Genomes sites file gives each super-population's frequency in INFO (AFR_AF, EUR_AF...): one panel column each.",
             {"source": _H + "af_1000g_demo.vcf.gz", "kind": "freq", "groups": "AFR, AMR, EAS, EUR, SAS", "min_maf": 0.05}),
        ],
        ("rnaseq", "Count reads per gene"): [
            ("Four simulated RNA-seq samples on the toy genome",
             "Two control and two treated samples of reads from the exons of the toy genome's genes (`annotation.gff3`). "
             "Each read that falls in the exons of one gene is counted for it: expect one column per sample, one row per gene.",
             {"bams": [_E + f"rnaseq_{_s}.bam" for _s in ("ctrl_1", "ctrl_2", "treat_1", "treat_2")], "annotation": _E + "annotation.gff3"}),
        ],
        ("rnaseq", "Differential expression"): [
            ("The four simulated samples: treated against control",
             "Three genes were planted as changed by the treatment: gene03 (8 times up), gene07 (3 times up) and gene10 "
             "(8 times down). With only two samples per group, expect those three, and maybe a false one.",
             {"counts": _E + "rnaseq_counts.tsv", "samples": _E + "rnaseq_samples.tsv", "condition": "condition", "reference": "control"}),
            ("Airway: dexamethasone against untreated (Himes et al. 2014)",
             "Airway smooth-muscle cells of four donors, treated with the corticosteroid dexamethasone or not (GEO GSE52778), "
             "compared within each donor (`~ cell_line + dex`). Expect thousands of genes changed, among them the classic "
             "glucocorticoid-response genes FKBP5, TSC22D3, DUSP1, KLF15 and ZBTB16, strongly up.",
             {"counts": _R + "rnaseq/airway_counts.tsv.gz", "samples": _R + "rnaseq/airway_samples.tsv", "condition": "dex",
              "reference": "untrt", "factors": "cell_line", "names": _R + "rnaseq/airway_gene_names.tsv"}),
        ],
        ("rnaseq", "Test gene sets"): [
            ("Airway: a positive and a negative control set",
             "The dexamethasone result of the airway experiment (computed here the first time) against a GMT file of two sets: "
             "classic glucocorticoid-receptor target genes (FKBP5, TSC22D3, PER1, KLF15, DUSP1, ZBTB16...) and housekeeping genes. "
             "Expect the first strongly **up** (z near 5.6) and the second not significant: what a gene-set test should do. Then "
             "untick *Rank only the genes with an adjusted p-value* and submit again: with the ~9,800 barely expressed genes back "
             "in the ranking, the housekeeping set turns 'significant' although it hardly changes (mean log2 fold change 0.15). "
             "That is why lowly expressed genes are left out first.",
             {"results": _R + "rnaseq/airway_dex.de.tsv", "sets": _R + "rnaseq/airway_controls.gmt"}),
            ("Airway against the built-in spinal cord injury sets",
             "The same result against the 11 built-in sets (glial scar, inflammation, axon regeneration, myelin, laminin "
             "receptors...). Airway smooth-muscle cells are not spinal cord: the neuronal and myelin sets have fewer genes "
             "expressed, and few sets move in one direction. It shows how the screen reads; on RNA-seq of an injured cord "
             "(treated against control) the same sets read the scar and the regeneration programs.",
             {"results": _R + "rnaseq/airway_dex.de.tsv", "sets": "spinal_cord_injury"}),
        ],
        ("crispr", "Count guides in reads"): [
            ("The MAGeCK demo: two samples of 2,500 reads",
             "Reads of 30 bases whose first 20 are a guide of the library (2,550 guides). MAGeCK itself maps 1,453 and 1,471 of "
             "the 2,500 reads of each sample; the counts here are the same, guide by guide.",
             {"library": _R + "crispr/demo2_library.txt", "fastqs": [_R + "crispr/demo2_sample1.fastq.gz", _R + "crispr/demo2_sample2.fastq.gz"],
              "labels": "sample1, sample2"}),
        ],
        ("crispr", "Rank genes (MAGeCK RRA)"): [
            ("The MAGeCK demo screen (100 genes)",
             "Guide counts before and after a screen in HL60 and KBM7 cells. The genes whose guides drop out are the ones the "
             "cells need; the result is the same, line by line, as MAGeCK's.",
             {"counts": _R + "crispr/demo1_counts.txt", "treatment": "HL60.final, KBM7.final", "control": "HL60.initial, KBM7.initial"}),
            ("A leukemia screen: 73,151 guides, 7,115 genes (about a minute)",
             "The genome-wide screen of MAGeCK's tutorial in two leukemia cell lines. Expect ribosomal-protein genes (RPL3, RPL8, "
             "RPS2...) at the top of the negative selection: knocking out an essential gene kills the cell, and its guides vanish.",
             {"counts": _R + "crispr/leukemia_counts.csv.gz", "treatment": "HL60.final, KBM7.final", "control": "HL60.initial, KBM7.initial"}),
        ],
        ("meta", "Which species, and how much"): [
            ("ZymoBIOMICS EVEN community: 8 bacteria and 2 yeasts",
             "30,000 read pairs of a public run of a mock community with known composition (8 bacteria at 12% and 2 yeasts at "
             "2% of the DNA), against the sketches of 43 reference genomes and a human decoy. Expect the ten species; the "
             "proportions come out near the truth, less so for strains far from the reference genome (Listeria, Bacillus).",
             {"samples": [_R + "microbio/ERR2984773_zymo_even_D6300_first30k_R1.fastq.gz", _R + "microbio/ERR2984773_zymo_even_D6300_first30k_R2.fastq.gz"],
              "databases": [_R + "microbio/references_k31_s1000.sig"]}),
            ("ZymoBIOMICS LOG community: from 89% down to 0.0001%",
             "The same species in amounts that fall tenfold each step. With 30,000 read pairs only the abundant ones are seen: "
             "detection depends on depth.",
             {"samples": [_R + "microbio/ERR2935805_zymo_log_D6310_first30k_R1.fastq.gz", _R + "microbio/ERR2935805_zymo_log_D6310_first30k_R2.fastq.gz"],
              "databases": [_R + "microbio/references_k31_s1000.sig"]}),
        ],
        ("meta", "Classify the reads"): [
            ("EVEN community, read by read, against the ten species",
             "Each read pair is assigned to the species whose k-mers it shares most, with a denser sketch (scaled 20); host "
             "(human) reads are counted apart.",
             {"input": _R + "microbio/ERR2984773_zymo_even_D6300_first30k_R1.fastq.gz", "input2": _R + "microbio/ERR2984773_zymo_even_D6300_first30k_R2.fastq.gz",
              "databases": [_R + "microbio/zymo_species_k31_s20.sig"]}),
        ],
        ("meta", "Sketch reference genomes"): [
            ("The four isolates of screen 25",
             "One sketch per genome, with the same hashes as sourmash: small files that stand for whole genomes.",
             {"genomes": [_R + f"microbio/isolates/{_g}.fna.gz" for _g in ("kpneumoniae_KPNIH1_KPC-3", "kpneumoniae_ATCC_BAA-2146_NDM-1",
                                                                            "saureus_USA300_MRSA", "saureus_NCTC8325_MSSA")],
              "output": _R + "microbio/isolates.sig"}),
        ],
        ("amr", None): [
            ("K. pneumoniae KPNIH1: a KPC carbapenemase outbreak strain",
             "The isolate of a hospital outbreak (ST258). Expect blaKPC-3, a carbapenemase, among its resistance genes.",
             {"inputs": [_R + "microbio/isolates/kpneumoniae_KPNIH1_KPC-3.fna.gz"]}),
            ("K. pneumoniae ATCC BAA-2146: NDM-1",
             "A reference strain carrying the metallo-beta-lactamase NDM-1 and several other resistance genes.",
             {"inputs": [_R + "microbio/isolates/kpneumoniae_ATCC_BAA-2146_NDM-1.fna.gz"]}),
            ("S. aureus USA300 (MRSA)",
             "A community-acquired MRSA: expect mecA, whose copy here differs slightly from the catalog's alleles (identity just "
             "under 100%).",
             {"inputs": [_R + "microbio/isolates/saureus_USA300_MRSA.fna.gz"]}),
            ("S. aureus NCTC 8325 (no mecA)",
             "A methicillin-susceptible strain: no mecA, but the genes S. aureus carries by nature.",
             {"inputs": [_R + "microbio/isolates/saureus_NCTC8325_MSSA.fna.gz"]}),
        ],
        ("outbreak", "Distances and a tree"): [
            ("SARS-CoV-2: 20 lineages, rooted on lineage A",
             "One genome per Pango lineage, from lineage A (January 2020) to KP.3.1.1 (2024), in one multi-FASTA: each "
             "record is a genome. Expect the two Delta genomes side by side and the Omicron lineages (BA.1 to KP.3.1.1) "
             "in one clade, far from the rest. Rooted on lineage A, the oldest; at the midpoint the root would fall inside "
             "Omicron, the longest branch.",
             {"genomes": [_R + "outbreak/sars_cov_2_20_lineages.fasta"], "records": True, "outgroup": "A_MN985325.1"}),
            ("Klebsiella pneumoniae: 13 isolates of known sequence types",
             "Sketches of 13 public genomes, named by strain and sequence type (ST). Isolates of the same ST should be "
             "neighbors in the tree, at about 99.7-99.9% ANI: ST258 (the KPC outbreak clone), ST11, ST15, ST307. ST258 and "
             "ST11 then group together: both belong to clonal group 258.",
             {"genomes": [_R + "outbreak/klebsiella_pneumoniae_13_k31_s100.sig"]}),
            ("Four bacterial isolates, from their genomes (FASTA)",
             "Two K. pneumoniae and two S. aureus, sketched here from their assemblies (a few seconds): two species, so two "
             "distant pairs; within each pair, the ANI says how close the strains are.",
             {"genomes": [_R + f"microbio/isolates/{_g}.fna.gz" for _g in ("kpneumoniae_KPNIH1_KPC-3", "kpneumoniae_ATCC_BAA-2146_NDM-1",
                                                                            "saureus_USA300_MRSA", "saureus_NCTC8325_MSSA")]}),
        ],
        ("outbreak", "Mutations of viral genomes"): [
            ("SARS-CoV-2: 20 lineages against the Wuhan-Hu-1 reference",
             "Each genome is aligned to NC_045512.2 and its changes are named as Nextclade names them (the same changes, "
             "genome by genome, as Nextclade 3.29). Look for S:D614G and ORF1b:P314L in every lineage after B.1, S:N501Y "
             "in Alpha, Beta, Gamma and Omicron, S:L452R and S:P681R in Delta, and the many spike changes of Omicron.",
             {"reference": _R + "outbreak/NC_045512.2.fna.gz", "annotation": _R + "outbreak/NC_045512.2.gff3.gz",
              "genomes": [_R + "outbreak/sars_cov_2_20_lineages.fasta"]}),
        ],
        ("amplicon", "Amplicon variants (ASVs)"): [
            ("mothur MiSeq SOP: 19 mouse stool samples and a Mock community",
             "16S V4 reads (2 x 250) of the mothur tutorial: mice sampled early (days 0-9) and late (days 141-150) after "
             "weaning, 2,000 read pairs each, and a Mock community of 20 known strains. Each pair is merged, reads with more "
             "than one expected error are dropped, the rest denoised into ASVs and chimeras removed. Expect about three "
             "quarters of the reads kept, and the Mock's 20 strains among the most abundant ASVs.",
             {"sheet": _R16 + "samples.tsv"}),
        ],
        ("amplicon", "Taxonomy"): [
            ("The ASVs of the SOP, with the RDP 19 classifier (V4)",
             "The naive Bayesian classifier of the RDP, trained on the V4 region of the RDP training set 19. Expect the "
             "mouse gut families: Muribaculaceae first, then Lachnospiraceae, Bacteroidaceae, Lactobacillaceae; few ASVs "
             "reach a genus at confidence 0.8, as is usual for 250-base V4 reads.",
             {"fasta": _R16 + "sop_asvs.fasta", "model": _R16 + "rdp19_v4_classifier.model", "table": _R16 + "sop_asv_table.tsv"}),
        ],
        ("amplicon", "Diversity"): [
            ("Early against late after weaning (Bray-Curtis)",
             "The 19 mouse samples: alpha diversity per sample, and how far apart the samples are (Bray-Curtis), drawn by "
             "principal coordinates. The gut community changes as the mice grow: expect early and late apart, and "
             "PERMANOVA p near 0.001.",
             {"table": _R16 + "sop_asv_table_mouse.tsv", "metadata": _R16 + "samples_design.tsv", "column": "time"}),
            ("The same with Jaccard (which ASVs are present, not how many)",
             "Presence and absence only: the first axis alone separates early from late.",
             {"table": _R16 + "sop_asv_table_mouse.tsv", "metadata": _R16 + "samples_design.tsv", "column": "time", "metric": "jaccard"}),
        ],
        ("amplicon", "Train a classifier"): [
            ("RDP training set 19, V4 region only (a few seconds)",
             "Only the stretch between the 515F and 806R primers of each reference sequence is used: a smaller model "
             "that suits V4 reads better than one trained on whole 16S genes.",
             {"reference": _R16 + "rdp_19_toGenus_trainset.fa.gz", "primer_f": "GTGYCAGCMGCCGCGGTAA", "primer_r": "GGACTACNVGGGTWTCTAAT",
              "output": _R16 + "my_v4_classifier.model"}),
        ],
        ("tmb", None): [
            ("HCC1395 breast-cancer cell line, chromosome 17",
             "The SEQC2 somatic reference calls of a triple-negative breast-cancer cell line (high-confidence SNVs and "
             "indels) counted over its high-confidence region of chromosome 17. These calls carry no depth or allele "
             "fraction, so those filters cannot apply (the screen says so). Expect about 15-20 mutations per megabase: a "
             "high burden.",
             {"vcf": _R + "tumor/HCC1395_somatic_chr17.vcf.gz", "bed": _R + "tumor/HCC1395_high_confidence_chr17.bed"}),
            ("The same, coding sequences only (as an exome or a panel counts)",
             "Only the variants inside the coding sequences (CDS) of the RefSeq annotation, over the coding bases of the "
             "region: what a whole-exome or panel test would count.",
             {"vcf": _R + "tumor/HCC1395_somatic_chr17.vcf.gz", "bed": _R + "tumor/HCC1395_high_confidence_chr17.bed",
              "annotation": _R + "tumor/refseq_chr17_cds.gff3.gz", "coding": True}),
        ],
        ("outbreak", "In-silico PCR"): [
            ("CDC N1 and N2 assays on 20 SARS-CoV-2 genomes",
             "The two primer pairs of the CDC test (with their probes) against every lineage. Both should give their "
             "product in all 20 genomes (72 and 67 bases); the table shows where a newer lineage has a mismatch under a "
             "primer or the probe, what assay designers watch for.",
             {"panel": _R + "outbreak/pcr_demo_panel.tsv", "inputs": [_R + "outbreak/sars_cov_2_20_lineages.fasta"], "records": True}),
            ("16S and mecA on four bacterial isolates",
             "Universal 16S primers give one product per rRNA operon (a bacterium has several copies); the mecA pair "
             "(methicillin resistance) gives its 533-base product only in the MRSA, not in the susceptible S. aureus.",
             {"panel": _R + "outbreak/pcr_demo_panel.tsv",
              "inputs": [_R + f"microbio/isolates/{_g}.fna.gz" for _g in ("kpneumoniae_KPNIH1_KPC-3", "kpneumoniae_ATCC_BAA-2146_NDM-1",
                                                                           "saureus_USA300_MRSA", "saureus_NCTC8325_MSSA")]}),
        ],
        ("chrom", "Genetic sex"): [
            ("HG001 (NA12878), chromosome X",
             "Two X chromosomes give heterozygous calls on chrX at about the rate of the autosomes: expect XX.",
             {"vcf": _H + "giab/HG001.chrX.vcf.gz"}),
            ("HG002, chromosomes X and Y",
             "One X gives no heterozygous call outside the pseudoautosomal regions, and there are calls on Y: expect XY.",
             {"vcf": _H + "giab/HG002.chrXY.vcf.gz"}),
            ("A file without chrX",
             "The demo cut of HG002 has no sex chromosomes: the screen says it cannot tell.",
             {"vcf": _H + "giab/HG002.vcf.gz"}),
        ],
        ("chrom", "Runs of homozygosity"): [
            ("A child of first cousins (simulated)",
             "A child of related parents inherits some stretches twice from the same ancestor: there, (almost) no heterozygous call. "
             "Expect two runs, chr1:20-32 Mb and chr2:5-9 Mb (F_ROH is computed against the whole human autosome length).",
             {"vcf": _E + "roh_demo.vcf", "min_length": 1_000_000}),
            ("HG002, with runs from 300 kb",
             _GIAB + " A person whose parents are not related has few and short runs.",
             {"vcf": _H + "giab/HG002.vcf.gz", "callable": _H + "giab/HG002.callable.bed", "min_length": 300_000}),
            ("The same child, strict method",
             "With the strict method any heterozygous call ends a run. On this simulated child (no genotyping errors) both "
             "methods find the same two runs; on real data, the strict one splits runs at every genotyping error.",
             {"vcf": _E + "roh_demo.vcf", "min_length": 1_000_000, "method": "strict"}),
        ],
        ("lookup", "Allele counts at a position (BAM)"): [
            ("Reads over a heterozygous SNP (simulated)",
             "80 simulated reads cover chr1:25,000 of the toy genome; about 40% carry the ALT allele, as a heterozygous variant "
             "looks. A few bases there have quality 2 and are left out by the quality filter. Below the counts, the reads as a "
             "genome browser draws them, against the toy genome, with the checks of strand bias and of the place of the ALT "
             "base in the reads.",
             {"bam": _E + "pileup_demo.bam", "position": "chr1:25000", "min_baseq": 13, "reference": _E + "toy_genome.fa",
              "annotation": _E + "annotation.gff3"}),
        ],
    })

    def make_demo_pack(kind):
        """Unpacks a demo pack ("human" or "research") into mmgen_examples/,
        from its zip in this folder (or downloaded, once it is published).
        Returns True when it unpacked it."""
        import shutil
        import urllib.request
        import zipfile
        folder, version, pack_name, url, what = DEMO_PACKS[kind]
        try:
            with open(os.path.join(folder, "VERSION"), encoding="utf-8") as fh:
                if fh.read().strip() == version:
                    return False
        except OSError:
            pass
        shutil.rmtree(folder, ignore_errors=True)  # an older version, or half unpacked
        os.makedirs(EXAMPLES_DIR, exist_ok=True)
        pack = next((c for c in (pack_name, os.path.join(EXAMPLES_DIR, pack_name)) if os.path.exists(c)), None)
        if pack is None and url:
            pack = os.path.join(EXAMPLES_DIR, pack_name)
            urllib.request.urlretrieve(url, pack)
        if pack is None:
            raise RuntimeError(
                f"this example needs the demo data `{pack_name}` ({what}). Put it in this folder: it is unpacked "
                f"into `{folder}/` the first time.")
        with zipfile.ZipFile(pack) as z:
            z.extractall(EXAMPLES_DIR)
        return True

    def make_human_demo():
        return make_demo_pack("human")

    # Simulated data sets kept in the package (multimerge_genomics.demos), each
    # in its own folder: written the first time an example needs them, and
    # again only when the package brings a newer version.
    SYNTHETIC_DEMOS = {
        "association": EXAMPLES_DIR + "/association",
        "contaminants": EXAMPLES_DIR + "/contaminants",
        "rnaseq_fastq": EXAMPLES_DIR + "/rnaseq_fastq",
    }

    def make_synthetic_demo(kind):
        """Writes a simulated demo into its folder. Returns True when it wrote it."""
        import importlib
        folder = SYNTHETIC_DEMOS[kind]
        demo = importlib.import_module("multimerge_genomics.demos." + kind)
        try:
            with open(os.path.join(folder, "VERSION"), encoding="utf-8") as fh:
                if fh.read().strip() == str(demo.DEMO_VERSION):
                    return False
        except OSError:
            pass
        import shutil
        shutil.rmtree(folder, ignore_errors=True)  # an older version: what was made from it goes too
        demo.write_demo(folder)
        with open(os.path.join(folder, "VERSION"), "w", encoding="utf-8") as fh:
            fh.write(f"{demo.DEMO_VERSION}\n")
        return True

    # Screen 29 (AI): problems and tasks on the simulated example files. They
    # need an API key; nothing is sent until the form is submitted.
    EXAMPLES.update({
        ("ai", "Describe your problem"): [
            ("Two VCFs that share no variant",
             "Two files of the same simulated genome, where one names chromosomes `chr1` and the other `1`, the most common "
             "first stumble. Expect the AI to spot the naming, send you to screen 10 and write a script that renames, sorts "
             "and compares them.",
             {"problem": "I have two VCF files of the same genome, cohort.vcf and ensembl_names.vcf, from two different "
                         "pipelines. When I look for the variants they have in common I find none, although they should "
                         "share many. What is going on, and how do I compare them?",
              "folder": EXAMPLES_DIR, "level": "student"}),
            ("A file a tool calls invalid",
             "`with_problems.vcf` has defects planted on purpose. Expect a path that starts with the validator (screen 10), "
             "what each problem means and how to fix it.",
             {"problem": "A program refuses my file with_problems.vcf and says it is invalid, and I do not understand its "
                         "error message. What is wrong with the file, and how do I fix it and get a sorted, indexed copy?",
              "folder": EXAMPLES_DIR, "level": "student"}),
            ("Reads from a sequencing facility",
             "Paired reads with an adapter and low-quality ends. Expect the quality report and the trimming of screen 12, "
             "with the reason behind each parameter.",
             {"problem": "A sequencing facility sent me paired reads, reads_R1.fastq.gz and reads_R2.fastq.gz. How do I "
                         "check their quality and remove the adapters before I use them, and which numbers tell me they "
                         "are good enough?",
              "folder": EXAMPLES_DIR, "level": "researcher"}),
        ],
        ("ai", "Let the AI run a task"): [
            ("Check a broken VCF and make a clean copy",
             "The agent validates `with_problems.vcf`, explains what it finds and tries to write a clean, sorted copy. You "
             "approve each round of calls before it runs.",
             {"task": "Validate with_problems.vcf, tell me what is wrong with it, then write a sorted copy next to it if "
                      "the engine can, and check the copy with the validator.",
              "folder": EXAMPLES_DIR, "rounds": 6}),
            ("Merge three samples and describe the result",
             "Three small VCFs out of order: the agent should sort and merge them into one indexed file and give its "
             "statistics.",
             {"task": "Merge sample_a.vcf, sample_b.vcf.gz and sample_c.vcf into one sorted file merged_samples.vcf.gz "
                      "with its index, then give me the statistics of the merged file.",
              "folder": EXAMPLES_DIR, "rounds": 6}),
        ],
    })

    _A = EXAMPLES_DIR + "/association/"
    _ASSOC_COHORT = ("A simulated study of 122 people after a treatment (120 genotyped at 40 real SNVs of five candidate "
                     "genes; the genotypes and outcomes are drawn by a program and say nothing about these variants in real people). ")
    _COVS = "age, sex, center, anc_afr, anc_eas"
    EXAMPLES.update({
        ("assoc", "Test variants against an outcome"): [
            ("Recovery score, adjusted for age, sex, hospital and ancestry",
             _ASSOC_COHORT + "A real effect was planted in OPRM1 rs12190259 (+6 points per copy). Expect 115 people analysed "
             "(the report says who was left out and why) and that variant as the only one with FDR < 0.05, at about +7.6 "
             "points [4.4, 10.9]. Two of the three ancestry proportions are covariates: the three sum to 1.",
             {"vcfs": [_A + "cohort.vcf.gz"], "phenotypes": _A + "phenotypes.tsv", "outcome": "recovery_score", "covariates": _COVS}),
            ("The same without ancestry: a false hit from population structure",
             _ASSOC_COHORT + "BDNF rs2030324 has no effect, but it is more frequent in people with more AFR-like ancestry, "
             "and ancestry changes the outcome here. Adjusted only for age, sex and hospital it comes second (p ≈ 0.003); "
             "compare with the first example, where it has p ≈ 0.87.",
             {"vcfs": [_A + "cohort.vcf.gz"], "phenotypes": _A + "phenotypes.tsv", "outcome": "recovery_score", "covariates": "age, sex, center"}),
            ("Recovered yes/no: less power, and Firth's method for a rare variant",
             _ASSOC_COHORT + "The same planted variant against a yes/no outcome: odds ratio about 2.0 [0.9, 4.4], p ≈ 0.07 — "
             "missed, and the power line says why (about 4.4 would be needed). IL6 rs2069830 has 3 carriers, all "
             "recovered: the plain logistic fit separates and Firth's method is used.",
             {"vcfs": [_A + "cohort.vcf.gz"], "phenotypes": _A + "phenotypes.tsv", "outcome": "recovered", "case_value": "yes", "covariates": _COVS}),
            ("With 1,000 permutations: empirical and family-wise p",
             _ASSOC_COHORT + "As the first example, with 1,000 permutations (Freedman-Lane): the planted variant keeps an "
             "empirical and a family-wise p of 0.001 — none of the shuffled data sets gave a stronger result across the 40 tests.",
             {"vcfs": [_A + "cohort.vcf.gz"], "phenotypes": _A + "phenotypes.tsv", "outcome": "recovery_score", "covariates": _COVS,
              "permutations": 1000}),
        ],
        ("assoc", "Study power: how many people?"): [
            ("Five patients: can anything be detected?",
             "A first trial with 5 people, a variant at frequency 0.25 and an odds ratio of 2. Expect a power under 10% "
             "(the 2,000 simulated studies agree), about 160 people needed for 80%, and no odds ratio up to 10 detectable with 5.",
             {"outcome": "binary", "n": 5, "allele_freq": 0.25, "effect": 2.0, "simulations": 2000}),
            ("120 people, 40 candidate variants, a measured outcome",
             "The demo cohort: a quantitative outcome, an allele at frequency 0.16, half a standard deviation per copy, "
             "40 tests (α = 0.05/40) and 5 covariates.",
             {"outcome": "quantitative", "n": 120, "allele_freq": 0.16, "effect": 0.5, "tests": 40, "covariates": 5}),
            ("A genome-wide study (α = 5×10⁻⁸)",
             "A common variant (frequency 0.3) with a modest odds ratio of 1.2, as most found by genome-wide studies, "
             "tested at the genome-wide threshold: 10,000 people give about 70% power, and about 11,000 are needed for 80%.",
             {"outcome": "binary", "n": 10000, "allele_freq": 0.3, "effect": 1.2, "alpha": 5e-8}),
        ],
    })

    _QF = EXAMPLES_DIR + "/rnaseq_fastq/"
    _QREADS = [_QF + f"{_s}_R{_m}.fastq.gz" for _s in ("control_1", "control_2", "control_3", "treated_1", "treated_2", "treated_3") for _m in (1, 2)]
    _QDEMO = ("A simulated transcriptome (269 transcripts of 140 genes, many genes with isoforms that share exons) and six "
              "paired-end samples, 3 control and 3 treated, about 20,000 pairs each. ")
    _CF = EXAMPLES_DIR + "/contaminants/"
    EXAMPLES.update({
        ("rnaseq", "Build a transcript index"): [
            ("The demo transcriptome: 269 transcripts of 140 genes",
             _QDEMO + "The genes come from the FASTA headers (`gene:` and `gene_symbol:`, as Ensembl writes them). Expect about "
             "120,000 k-mers in about 450 equivalence classes: a class is a set of transcripts that share a k-mer, so the "
             "isoforms of a gene make classes of their own.",
             {"fasta": _QF + "transcripts.fa"}),
        ],
        ("rnaseq", "Quantify from FASTQ (pseudoalignment)"): [
            ("The six samples, paired-end, with their conditions",
             _QDEMO + "Expect 100% of the pairs pseudoaligned (the reads come only from these transcripts), about 37% to a "
             "single transcript (the others fit several isoforms: the EM shares them out) and a fragment length measured "
             "near 200 ± 25. Then pick the task **Differential expression** and its example *The six samples quantified "
             "from FASTQ*: Gene7, Gene19 and Gene33 go up 4×, Gene12, Gene48 and Gene71 down.",
             {"index": _QF + "transcripts.idx", "fastqs": _QREADS, "conditions": _QF + "samples.tsv", "output": _QF + "quant"}),
            ("The same reads as single-end (R1 only, fragments of 200 ± 25)",
             _QDEMO + "Only the first read of each pair: the fragment length cannot be measured from one read, so it is given. "
             "Fewer reads land on a single transcript than with the pairs.",
             {"index": _QF + "transcripts.idx", "fastqs": _QREADS[0::2], "conditions": _QF + "samples.tsv",
              "fragment_length": 200, "fragment_sd": 25, "output": _QF + "quant_single"}),
        ],
        ("rnaseq", "See how one read is pseudoaligned"): [
            ("A read that fits two isoforms, and its mate, that fits one",
             "The read alone holds k-mers shared by both isoforms of Gene80 (DEMO0080-201 and -202): it is compatible with the "
             "two. Its mate covers an exon only -202 has, so the pair goes to one transcript. That is why pairs resolve "
             "isoforms better than single reads.",
             {"index": _QF + "transcripts.idx",
              "read": "GGTCTCCCAACCGATCACGTTGTGTTAAGTGTTTTAACGCGATACGGGTCATTGAGATAGATTTACACGAGACTA",
              "read2": "AGTTGCCCGCAGAATTAGCCTTGGCCTCTCTTATCCTAGAGTGTCTATGATTCGGTAGACGCGCTCGCTTGTGCG"}),
        ],
        ("contam", "Screen samples for adventitious agents"): [
            ("A raw material and its reagent blank (simulated)",
             "A simulated cell-culture harvest (20,121 read pairs, mostly host) and its reagent blank (3,000 pairs), against a "
             "database of synthetic genomes: a host, two viruses, two bacteria, a mycoplasma, a yeast and a process control. "
             "Expect, in the raw material: the 36 kb virus (15 pairs) and the mycoplasma (12 pairs) **detected** — few reads, "
             "but their own and spread over their genomes; the reagent bacterium **background** (it is in the blank at a much "
             "higher rate); the small virus **below threshold** — 40 pairs, all on one 500-base region (4 of its 20 genome "
             "windows), as a vector element would give; and the process control recovered at 100%.",
             {"database": _CF + "demo.mmgscreen", "fastqs": [_CF + "raw_material_R1.fq.gz", _CF + "raw_material_R2.fq.gz",
                                                             _CF + "reagent_blank_R1.fq.gz", _CF + "reagent_blank_R2.fq.gz"],
              "negative_controls": "reagent_blank", "spike_ins": _CF + "spike_ins.tsv"}),
            ("The raw material without its blank",
             "The same raw material, screened alone: with no negative control to compare with, the reagent bacterium (4 pairs, "
             "spread over its genome) is now **detected**. A blank processed with the samples is what tells a contaminant of "
             "the reagents from one of the material.",
             {"database": _CF + "demo.mmgscreen", "fastqs": [_CF + "raw_material_R1.fq.gz", _CF + "raw_material_R2.fq.gz"],
              "spike_ins": _CF + "spike_ins.tsv", "report": _CF + "screen_without_blank"}),
        ],
        ("contam", "Build a screening database"): [
            ("The demo references: a host, viruses, bacteria, a mycoplasma, a yeast and a control",
             "`references.tsv` lists each genome with its kind and file. The agents and the process control are indexed with "
             "their genome windows (for the breadth of coverage); the host goes into a compact set. Expect about 170,000 "
             "agent minimizers, 125,000 host minimizers and a SHA-256 for every file.",
             {"table": _CF + "references.tsv"}),
        ],
    })
    EXAMPLES[("rnaseq", "Differential expression")] = EXAMPLES.get(("rnaseq", "Differential expression"), []) + [
        ("The six samples quantified from FASTQ",
         _QDEMO + "The gene counts estimated from the reads (task *Quantify from FASTQ*; made now if they are not there yet). "
         "Expect exactly six genes with padj < 0.1: Gene7, Gene19 and Gene33 up (log2 fold change near +2), Gene12, Gene48 and "
         "Gene71 down (near −2). Gene90 switches isoforms without changing its total, so the gene-level test does not see it: "
         "compare its two isoforms in `transcript_tpm.tsv`.",
         {"counts": _QF + "quant/gene_counts.tsv", "samples": _QF + "quant/samples.tsv", "condition": "condition", "reference": "control",
          "names": _QF + "quant/gene_names.tsv"}),
    ]

    def example_picker(page, task=None):
        """The list of examples of a screen and task."""
        options = {NO_EXAMPLE: ""}
        for i, item in enumerate(EXAMPLES.get((page, task), [])):
            options[item[0]] = str(i)
        return mo.ui.dropdown(options=options, value=NO_EXAMPLE, label="🧪 **Try an example:**", full_width=True)

    # Example files made from others the first time they are needed (after the
    # demo packs are unpacked): a result an example starts from.
    def _airway_de():
        # the engine itself: preparing an example is not a call of the session
        getattr(mmg, "_module", mmg).differential_expression(
            _R + "rnaseq/airway_counts.tsv.gz", _R + "rnaseq/airway_samples.tsv", "dex", "untrt",
            output=_R + "rnaseq/airway_dex.de.tsv", factors=["cell_line"], gene_names=_R + "rnaseq/airway_gene_names.tsv")

    def _airway_controls():
        _sets = [
            ("GLUCOCORTICOID_RESPONSE", "Classic glucocorticoid-receptor target genes (a positive control for dexamethasone)",
             ["FKBP5", "TSC22D3", "PER1", "KLF15", "DUSP1", "ZBTB16", "KLF9", "ERRFI1", "SGK1", "PDK4", "CRISPLD2"]),
            ("HOUSEKEEPING", "Housekeeping genes (a negative control)",
             ["GAPDH", "ACTB", "RPLP0", "B2M", "PPIA", "TBP", "HPRT1", "PGK1", "GUSB", "YWHAZ"]),
        ]
        with open(_R + "rnaseq/airway_controls.gmt", "w", encoding="utf-8", newline="\n") as _fh:
            for _name, _desc, _genes in _sets:
                _fh.write("\t".join([_name, _desc, *_genes]) + "\n")

    _C = EXAMPLES_DIR + "/contaminants/"
    _Q = EXAMPLES_DIR + "/rnaseq_fastq/"

    def _screen_db():
        import csv
        with open(_C + "references.tsv", encoding="utf-8", newline="") as _fh:
            _refs = list(csv.DictReader(_fh, delimiter="\t"))
        getattr(mmg, "_module", mmg).build_screen_database([_C + _r["file"] for _r in _refs], _C + "demo.mmgscreen",
                                                          kinds=[_r["kind"] for _r in _refs], names=[_r["name"] for _r in _refs])

    def _transcript_index():
        getattr(mmg, "_module", mmg).build_transcript_index(_Q + "transcripts.fa", _Q + "transcripts.idx")

    def _quantified():
        if not os.path.exists(_Q + "transcripts.idx"):
            _transcript_index()
        _names = ("control_1", "control_2", "control_3", "treated_1", "treated_2", "treated_3")
        getattr(mmg, "_module", mmg).quantify_samples(_Q + "transcripts.idx", [
            {"name": _n, "r1": _Q + f"{_n}_R1.fastq.gz", "r2": _Q + f"{_n}_R2.fastq.gz", "condition": _n.split("_")[0]}
            for _n in _names], _Q + "quant")

    DERIVED_EXAMPLE_FILES = {_R + "rnaseq/airway_dex.de.tsv": _airway_de, _R + "rnaseq/airway_controls.gmt": _airway_controls,
                             _C + "demo.mmgscreen": _screen_db, _Q + "transcripts.idx": _transcript_index,
                             _Q + "quant/gene_counts.tsv": _quantified}

    def example_for(page, task, picked):
        """The form values of the picked example and a note about it; writes
        the example files the first time. ({}, None) when none is picked."""
        items = EXAMPLES.get((page, task), [])
        if not picked or not str(picked).isdigit() or int(picked) >= len(items):
            return {}, None
        label, about, values = items[int(picked)]
        flat = [x for v in values.values() for x in (v if isinstance(v, list) else [v])]
        try:
            with mo.status.spinner("Writing the example files..."):
                make_examples()
                for _kind, (_folder, *_rest) in DEMO_PACKS.items():
                    if any(isinstance(x, str) and x.startswith(_folder + "/") for x in flat):
                        make_demo_pack(_kind)
                for _kind, _folder in SYNTHETIC_DEMOS.items():
                    if any(isinstance(x, str) and x.startswith(_folder + "/") for x in flat):
                        make_synthetic_demo(_kind)
                for x in flat:
                    if x in DERIVED_EXAMPLE_FILES and not os.path.exists(x):
                        DERIVED_EXAMPLE_FILES[x]()
        except Exception as err:
            return {}, mo.callout(mo.md(f"**🔥 ERROR:** the example files could not be written in `{EXAMPLES_DIR}/`: {err}"), kind="danger")
        return dict(values), mo.callout(mo.md(
            f"**🧪 {label}**\n\n{about}\n\nThe form below is filled in with this example (its files are in `{EXAMPLES_DIR}/`): "
            "submit it as it is, or change a field first."), kind="info")

    def label_of(options, value, default):
        """The label a dropdown shows for one of its values (examples keep values)."""
        return next((k for k, v in options.items() if v == value), default)

    return EXAMPLES, HUMAN_DIR, RESEARCH_DIR, example_for, example_picker, label_of, make_demo_pack, make_examples, make_human_demo


@app.cell
def main_menu(get_page, mmg, mo, set_page):
    mo.stop(get_page() != "home")

    btn_menu_qc = mo.ui.button(label="🔬 1. QC & Ingestion (FASTQ / FASTA)", on_change=lambda _: set_page("qc"))
    btn_menu_cluster = mo.ui.button(label="🧬 2. Taxonomic Clustering", on_change=lambda _: set_page("cluster"))
    btn_menu_assembly = mo.ui.button(label="🧩 3. Genome Assembly (De Bruijn)", on_change=lambda _: set_page("assembly"))
    btn_menu_variant = mo.ui.button(label="⚠️ 4. Variant Calling (SNPs)", on_change=lambda _: set_page("variant"))
    btn_menu_annotation = mo.ui.button(label="📘 5. Functional Annotation", on_change=lambda _: set_page("annotation"))
    btn_menu_bloom = mo.ui.button(label="🔍 6. Bloom Filters", on_change=lambda _: set_page("bloom"))
    btn_menu_scaffold = mo.ui.button(label="🖇️ 7. Genome Scaffolding", on_change=lambda _: set_page("scaffold"))
    btn_menu_vcf = mo.ui.button(label="🗄️ 8. VCF Multimerge & Sorting", on_change=lambda _: set_page("vcf"))
    btn_menu_sort = mo.ui.button(label="📏 9. Sort, Index & Intervals", on_change=lambda _: set_page("sort"))
    btn_menu_validate = mo.ui.button(label="🩺 10. Validation & Chromosome Names", on_change=lambda _: set_page("validate"))
    btn_menu_vcftools = mo.ui.button(label="🧪 11. VCF Stats, Filter & Normalize", on_change=lambda _: set_page("vcftools"))
    btn_menu_fastq = mo.ui.button(label="📈 12. FASTQ Quality & Trimming", on_change=lambda _: set_page("fastq"))
    btn_menu_kmercount = mo.ui.button(label="🔢 13. K-mer Counting", on_change=lambda _: set_page("kmercount"))
    btn_menu_convert = mo.ui.button(label="📦 14. Formats & Conversion", on_change=lambda _: set_page("convert"))
    btn_menu_guide = mo.ui.button(label="📚 15. Which Test Sees What (a guide)", on_change=lambda _: set_page("guide"))
    btn_menu_kin = mo.ui.button(label="👪 16. Relatedness, Trios & Paternity", on_change=lambda _: set_page("kin"))
    btn_menu_anc = mo.ui.button(label="🌍 17. Ancestry & Polygenic Scores", on_change=lambda _: set_page("anc"))
    btn_menu_findings = mo.ui.button(label="🧾 18. Disease Findings & Gene Panels", on_change=lambda _: set_page("findings"))
    btn_menu_pgx = mo.ui.button(label="💊 19. Pharmacogenomics", on_change=lambda _: set_page("pgx"))
    btn_menu_chrom = mo.ui.button(label="🧫 20. Sex Chromosomes & Homozygosity", on_change=lambda _: set_page("chrom"))
    btn_menu_lookup = mo.ui.button(label="🔎 21. Variant Lookup & Allele Counts", on_change=lambda _: set_page("lookup"))
    btn_menu_assoc = mo.ui.button(label="📐 30. Genotype & Outcome: Association and Power", on_change=lambda _: set_page("assoc"))
    btn_menu_rnaseq = mo.ui.button(label="📊 22. RNA-seq: Gene Counts & Differential Expression", on_change=lambda _: set_page("rnaseq"))
    btn_menu_crispr = mo.ui.button(label="✂️ 23. CRISPR Screens (MAGeCK)", on_change=lambda _: set_page("crispr"))
    btn_menu_meta = mo.ui.button(label="🦠 24. Metagenomics: Species & Abundance", on_change=lambda _: set_page("meta"))
    btn_menu_amr = mo.ui.button(label="💉 25. Antimicrobial Resistance Genes", on_change=lambda _: set_page("amr"))
    btn_menu_outbreak = mo.ui.button(label="🌳 26. Outbreaks: Trees, Viral Mutations & PCR", on_change=lambda _: set_page("outbreak"))
    btn_menu_amplicon = mo.ui.button(label="🌿 27. 16S Amplicons: Taxonomy & Diversity", on_change=lambda _: set_page("amplicon"))
    btn_menu_tmb = mo.ui.button(label="🎯 28. Tumor Mutational Burden (TMB)", on_change=lambda _: set_page("tmb"))
    btn_menu_contam = mo.ui.button(label="🛡️ 31. Contaminant Screening: Adventitious Agents", on_change=lambda _: set_page("contam"))
    btn_menu_ai = mo.ui.button(label="🤖 29. AI Assistant: Describe Your Problem", on_change=lambda _: set_page("ai"))

    # What this installation really has: the GPU modules fall back to the CPU.
    _engine = mmg.engine_info()
    _where = "NVIDIA GPU (CUDA)" if _engine["gpu"] else "CPU"
    _grid = 'style="display: grid; grid-template-columns: 1fr 1fr; gap: 15px; max-width: 850px; margin: 0 auto;"'
    _section = 'style="color: gray; margin: 30px 0 12px 0; text-align: center;"'

    menu_screen = mo.md(f"""
    <div style="text-align: center; margin-top: 40px; font-family: sans-serif;">
        <h1 style="font-size: 3em; margin-bottom: 0px; text-align: center;">🧬 MMGEN</h1>
        <h3 style="color: gray; margin-top: 5px;">High-Performance Genomics Framework</h3>
        <p style="margin-bottom: 10px;">Rust engine with optional NVIDIA CUDA acceleration. Pick a module to start:</p>

        <h4 {_section}>K-mer modules</h4>
        <div {_grid}>
            {btn_menu_qc}
            {btn_menu_cluster}
            {btn_menu_assembly}
            {btn_menu_variant}
            {btn_menu_annotation}
            {btn_menu_bloom}
            {btn_menu_scaffold}
        </div>

        <h4 {_section}>File tools, for files of any size</h4>
        <div {_grid}>
            {btn_menu_vcf}
            {btn_menu_sort}
            {btn_menu_validate}
            {btn_menu_vcftools}
            {btn_menu_fastq}
            {btn_menu_kmercount}
            {btn_menu_convert}
        </div>

        <h4 {_section}>Human genetics — experimental, for teaching and research</h4>
        <div {_grid}>
            {btn_menu_guide}
            {btn_menu_kin}
            {btn_menu_anc}
            {btn_menu_findings}
            {btn_menu_pgx}
            {btn_menu_chrom}
            {btn_menu_lookup}
            {btn_menu_assoc}
        </div>

        <h4 {_section}>Research, microbiology and oncology — experimental</h4>
        <div {_grid}>
            {btn_menu_rnaseq}
            {btn_menu_crispr}
            {btn_menu_meta}
            {btn_menu_amr}
            {btn_menu_outbreak}
            {btn_menu_amplicon}
            {btn_menu_tmb}
            {btn_menu_contam}
        </div>

        <h4 {_section}>AI assistant — optional, uses the Claude API</h4>
        <div style="display: grid; grid-template-columns: 1fr; gap: 15px; max-width: 420px; margin: 0 auto;">
            {btn_menu_ai}
        </div>
        <p style="color: gray; font-size: 0.85em; margin-top: 40px;">Engine {_engine["version"]} · running on: {_where} · {_engine["threads"]} threads</p>
    </div>
    """)

    menu_screen
    return (
        btn_menu_ai,
        btn_menu_amr,
        btn_menu_outbreak,
        btn_menu_amplicon,
        btn_menu_tmb,
        btn_menu_anc,
        btn_menu_annotation,
        btn_menu_assoc,
        btn_menu_crispr,
        btn_menu_meta,
        btn_menu_assembly,
        btn_menu_bloom,
        btn_menu_chrom,
        btn_menu_contam,
        btn_menu_kin,
        btn_menu_cluster,
        btn_menu_convert,
        btn_menu_fastq,
        btn_menu_findings,
        btn_menu_guide,
        btn_menu_lookup,
        btn_menu_pgx,
        btn_menu_kmercount,
        btn_menu_qc,
        btn_menu_rnaseq,
        btn_menu_scaffold,
        btn_menu_sort,
        btn_menu_validate,
        btn_menu_variant,
        btn_menu_vcf,
        btn_menu_vcftools,
        menu_screen,
    )


@app.cell
def qc_page_ui(example_picker, get_page, mo, set_page):
    mo.stop(get_page() != "qc")

    btn_back_qc = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    qc_example = example_picker("qc")

    _ui = mo.vstack([
        btn_back_qc,
        mo.md("<br>\n\n# 🧬 1. Data Ingestion & QC (Batch / GZIP)"),
        mo.callout(mo.md("The engine reads **several files at once** (for example paired-end R1 and R2) and decompresses `.gz` files on the fly, straight into memory, so nothing extra is written to disk."), kind="info"),
        qc_example,
    ])
    _ui
    return btn_back_qc, qc_example


@app.cell
def qc_page_form(example_for, files_with, get_page, make_form, mo, qc_example):
    mo.stop(get_page() != "qc")
    _ex, _about = example_for("qc", None, qc_example.value)
    form_qc = make_form({
        "files": mo.ui.multiselect(options=files_with("kmer_input"), value=_ex.get("files", []), label="🧬 DNA files (FASTA or FASTQ; several allowed, .gz too):", full_width=True),
        "kmer_size": mo.ui.slider(start=21, stop=63, step=2, value=_ex.get("kmer_size", 55), show_value=True, label="K-mer size:"),
        "min_quality": mo.ui.slider(start=0, stop=40, step=1, value=_ex.get("min_quality", 20), show_value=True, label="Minimum mean quality (FASTQ only):"),
    }, "🚀 Start batch ingestion")
    mo.vstack([_x for _x in (_about, form_qc) if _x is not None])
    return (form_qc,)


@app.cell
async def qc_page_backend(asyncio, beside, chosen, facts, figs, figure_panel, form_qc, get_page, mmg, mo, pl, split_name, warn):
    mo.stop(get_page() != "qc")
    mo.stop(form_qc.value is None, mo.md("⏳ *Waiting for the settings...*"))

    _paths = [_p for _p in form_qc.value.get("files", []) if chosen(_p)]
    mo.stop(not _paths, warn("Select at least one file."))

    _k = form_qc.value["kmer_size"]
    _min_quality = form_qc.value["min_quality"]

    _is_fasta = any(_c.lower().endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")) for _c in _paths)

    try:
        with mo.status.spinner(f"Processing {len(_paths)} file(s) in parallel (K={_k})..."):
            await asyncio.sleep(0)

            if _is_fasta:
                _kmer_bytes = await asyncio.to_thread(mmg.read_fasta_batch, _paths, _k)
            else:
                _kmer_bytes = await asyncio.to_thread(mmg.read_fastq_batch, _paths, _min_quality, _k)

        with mo.status.spinner("Saving the k-mers (Arrow file)..."):
            await asyncio.sleep(0)
            _df = pl.DataFrame({"K-mers": [_kmer_bytes]})

            # Next to the (first) input: x.fastq.gz -> x_K55.arrow.
            _base_name = "Multi_Batch" if len(_paths) > 1 else split_name(_paths[0])[0]
            _out_name = beside(_paths[0], f"{_base_name}_K{_k}.arrow")
            _df.write_ipc(_out_name)

        _blocks = [
            mo.callout(mo.md(f"""
            ### ✅ Batch processed successfully!
            * **Files read:** `{len(_paths)}`
            * **K-mer size:** `{_k}`
            * **Distinct k-mers extracted:** `{len(_kmer_bytes) // 16:,}`
            * **Arrow file written:** `{_out_name}`
            """), kind="success")
        ]
        if not _is_fasta:
            # Reads: the k-mer spectrum says, before any assembly, the error rate, the coverage, the genome size and
            # whether the genome is diploid with heterozygous sites.
            with mo.status.spinner("Counting how often each k-mer was seen (the spectrum)..."):
                _kk = min(int(_k), 31)
                _kk -= 1 if _kk % 2 == 0 else 0  # canonical k-mers need an odd k
                _cnt = await asyncio.to_thread(mmg.count_kmers, _paths, None, k=_kk)
            _model = figs.spectrum_model(_cnt["histogram"], k=_cnt["k"])
            _rows = []
            if _model.get("genome_size"):
                _rows = [("Coverage peak", f"`{_model['peak']}×`" + (f" (heterozygous peak at `{_model['het_peak']}×`)" if _model.get("het_peak") and _model["het_peak"] != _model["peak"] else "")),
                         ("Genome size, from the spectrum", f"`{figs.fmt_bp(_model['genome_size'])}`"),
                         ("Error k-mers", f"`{100 * (_model['error_fraction'] or 0):.1f}%` of the k-mers counted")]
            _blocks.append(mo.hstack([
                figure_panel(figs.kmer_spectrum(_cnt["histogram"], title=f"K-mer spectrum of the reads (k = {_cnt['k']}, log scale)", model=_model),
                             "reads_kmer_spectrum",
                             "The tall bars at the left are k-mers seen once or twice: sequencing errors. The hump is the genome, "
                             "at its coverage. Too few reads put the hump into the errors; a second hump at half the coverage is a "
                             "diploid genome with heterozygous sites. Screen 13 counts larger files and shows the linear scale."),
                mo.md(facts(_rows) if _rows else "*No coverage peak apart from the errors: too few reads, or too much variation, to read the spectrum.*"),
            ], widths=[2, 1], wrap=True))
        _panel = mo.vstack(_blocks)

    except Exception as e:
        _panel = mo.callout(mo.md(f"**🔥 ERROR:** `{str(e)}`"), kind="danger")

    _panel
    return


@app.cell
def cluster_page_ui(example_picker, get_page, mo, set_page):
    mo.stop(get_page() != "cluster")

    # A global variable (no leading "_"), like on the other screens: marimo only
    # keeps a component alive when it is bound to a global variable.
    btn_back_cluster = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    cluster_example = example_picker("cluster")

    _ui = mo.vstack([
        btn_back_cluster,
        mo.md("<br>\n\n# 🧬 Taxonomic Clustering"),
        mo.callout(mo.md(
            "Groups DNA by its **composition**, as metagenome binning tools (MetaBAT, CONCOCT) do: every sequence is cut "
            "into windows, each window is described by the frequency of its 136 tetranucleotides (4-base words, both "
            "strands together), the windows are projected with PCA and grouped with HDBSCAN. Windows of the same genome "
            "share a signature, so the groups tend to be species. Strains of one species are not separated, and "
            "plasmids or phages may form groups of their own. **Experimental**, for teaching."), kind="info"),
        cluster_example,
    ])
    _ui
    return btn_back_cluster, cluster_example


@app.cell
def cluster_page_form(cluster_example, example_for, files_with, get_page, make_form, mo):
    mo.stop(get_page() != "cluster")
    _ex, _about = example_for("cluster", None, cluster_example.value)
    _files = [f for f in files_with("fasta") if not f.startswith("(")]
    _files += [f for f in _ex.get("fastas", []) if f not in _files]
    form_cluster = make_form({
        "fastas": mo.ui.multiselect(options=_files, value=_ex.get("fastas", []), label="🧬 Sequences (FASTA: genomes, contigs or long reads):", full_width=True),
        "window": mo.ui.number(start=1000, stop=100000, step=1000, value=_ex.get("window", 5000), label="Window length (bases):"),
        "min_cluster": mo.ui.number(start=2, stop=1000, step=1, value=_ex.get("min_cluster", 5), label="Smallest group (windows):"),
        "dims": mo.ui.number(start=2, stop=20, step=1, value=_ex.get("dims", 5), label="PCA dimensions kept:"),
    }, "🧬 Cluster the sequences")
    mo.vstack([_x for _x in (_about, form_cluster) if _x is not None])
    return (form_cluster,)


@app.cell
async def cluster_page_backend(asyncio, error_panel, facts, form_cluster, get_page, mo, np, os, warn):
    mo.stop(get_page() != "cluster")
    mo.stop(form_cluster.value is None, mo.md("🧬 *Pick the sequences and submit the form.*"))
    _v = form_cluster.value
    _paths = list(_v["fastas"] or [])
    mo.stop(not _paths, warn("Select at least one FASTA file."))
    _W, _min, _dims = int(_v["window"]), int(_v["min_cluster"]), int(_v["dims"])

    def _compute():
        import gzip
        from sklearn.cluster import HDBSCAN
        from sklearn.decomposition import PCA

        code = np.full(256, -1, dtype=np.int64)
        for _i, _b in enumerate(b"ACGT"):
            code[_b] = _i
            code[ord(chr(_b).lower())] = _i
        # canonical tetranucleotides: a word and its reverse complement count together
        words = np.arange(256)
        rc = np.zeros(256, dtype=np.int64)
        for _w in range(256):
            _d = [(_w >> (2 * (3 - _j))) & 3 for _j in range(4)]
            rc[_w] = sum((3 - _d[3 - _j]) << (2 * (3 - _j)) for _j in range(4))
        canon = np.minimum(words, rc)
        keep = np.unique(canon)
        col = {int(_c): _j for _j, _c in enumerate(keep)}
        to_col = np.array([col[int(_c)] for _c in canon])

        def records(path):
            opener = gzip.open if path.endswith(".gz") else open
            name, parts = None, []
            with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    line = line.strip()
                    if line.startswith(">"):
                        if name is not None:
                            yield name, "".join(parts)
                        name, parts = line[1:].split()[0] if line[1:].split() else "sequence", []
                    elif line:
                        parts.append(line)
            if name is not None:
                yield name, "".join(parts)

        rows, labels, gcs = [], [], []
        n_seq = 0
        for path in _paths:
            for name, seq in records(path):
                n_seq += 1
                label = os.path.basename(path) if len(_paths) > 1 else name
                arr = code[np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)]
                for a in range(0, len(arr) - _W // 2, _W):
                    w = arr[a:a + _W]
                    if len(w) < _W // 2 or (w < 0).mean() > 0.1:
                        continue
                    ok = (w[:-3] >= 0) & (w[1:-2] >= 0) & (w[2:-1] >= 0) & (w[3:] >= 0)
                    four = (w[:-3] * 64 + w[1:-2] * 16 + w[2:-1] * 4 + w[3:])[ok]
                    counts = np.bincount(to_col[four], minlength=len(keep)).astype(float)
                    rows.append(counts / max(counts.sum(), 1))
                    labels.append(label)
                    valid = w[w >= 0]
                    gcs.append(float(((valid == 1) | (valid == 2)).mean()) if len(valid) else 0.0)
        if len(rows) < max(4, _min):
            return {"error": f"only {len(rows)} window(s) of {_W:,} bases: give longer sequences, or a shorter window."}
        X = np.log(np.array(rows) + 1e-4)
        X = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-9)
        pca = PCA(n_components=min(_dims, X.shape[0] - 1, X.shape[1]))
        Z = pca.fit_transform(X)
        lab = HDBSCAN(min_cluster_size=_min, copy=True).fit_predict(Z)
        return {"Z": Z, "labels": lab, "sources": labels, "gc": gcs, "sequences": n_seq,
                "explained": pca.explained_variance_ratio_}

    try:
        with mo.status.spinner("Cutting the sequences into windows and measuring their composition..."):
            _r = await asyncio.to_thread(_compute)
    except Exception as _err:
        _r = {"error": f"{type(_err).__name__}: {_err}"}
    mo.stop("error" in _r, error_panel(_r.get("error")))

    _lab, _src, _gc, _Z = _r["labels"], _r["sources"], _r["gc"], _r["Z"]
    _clusters = sorted(set(int(_c) for _c in _lab) - {-1})
    _noise = int((_lab == -1).sum())
    _names = {_c: f"group {_c + 1}" for _c in _clusters}
    _names[-1] = "ungrouped"

    # one row per group: size, GC and where its windows came from
    _summary = []
    for _c in _clusters + ([-1] if _noise else []):
        _idx = [_i for _i, _x in enumerate(_lab) if _x == _c]
        _from = {}
        for _i in _idx:
            _from[_src[_i]] = _from.get(_src[_i], 0) + 1
        _summary.append({"group": _names[_c], "windows": len(_idx), "bases": len(_idx) * _W,
                         "GC %": round(100 * sum(_gc[_i] for _i in _idx) / len(_idx), 1),
                         "from": ", ".join(f"{_k} ({_n})" for _k, _n in sorted(_from.items(), key=lambda _t: -_t[1]))})

    # the scatter: PC1 x PC2, one colour per group (validated palette), grey = ungrouped
    _PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
    _GREY = "#898781"
    _col = {_c: (_PAL[_i] if _i < len(_PAL) else _GREY) for _i, _c in enumerate(_clusters)}
    _col[-1] = _GREY
    _xs, _ys = _Z[:, 0], _Z[:, 1]
    _x0, _x1, _y0, _y1 = float(_xs.min()), float(_xs.max()), float(_ys.min()), float(_ys.max())
    _dx, _dy = (_x1 - _x0) or 1, (_y1 - _y0) or 1
    _x0, _x1, _y0, _y1 = _x0 - 0.06 * _dx, _x1 + 0.06 * _dx, _y0 - 0.06 * _dy, _y1 + 0.06 * _dy
    W, H, L, R, T, B = 600, 420, 50, 16, 16, 40
    X = lambda v: L + (v - _x0) / (_x1 - _x0) * (W - L - R)
    Y = lambda v: T + (H - T - B) - (v - _y0) / (_y1 - _y0) * (H - T - B)
    _order = sorted(range(len(_lab)), key=lambda _i: _lab[_i] != -1)  # grey first, underneath
    _dots = "".join(f'<circle cx="{X(float(_xs[_i])):.1f}" cy="{Y(float(_ys[_i])):.1f}" r="4" fill="{_col[int(_lab[_i])]}" '
                    f'stroke="var(--background, #fff)" stroke-width="1.5"><title>{_src[_i]}: {_names[int(_lab[_i])]}</title></circle>'
                    for _i in _order)
    _ev = _r["explained"]
    _axes = (f'<text x="{L + (W - L - R) / 2:.1f}" y="{H - 8}" text-anchor="middle" fill="currentColor" fill-opacity="0.8">PC1 ({100 * _ev[0]:.0f}%)</text>'
             f'<text transform="translate(14,{T + (H - T - B) / 2:.1f}) rotate(-90)" text-anchor="middle" fill="currentColor" fill-opacity="0.8">PC2 ({100 * _ev[1]:.0f}%)</text>')
    _legend = ('<div style="display:flex; flex-wrap:wrap; gap: 6px 16px; font: 12px system-ui, sans-serif; margin: 6px 0;">'
               + "".join(f'<span><span style="display:inline-block; width:10px; height:10px; border-radius:50%; background:{_col[_c]}; '
                         f'margin-right:5px; vertical-align:-1px"></span>{_names[_c]}</span>' for _c in _clusters + ([-1] if _noise else []))
               + "</div>")
    _chart = mo.Html(f'<div style="max-width:{W}px">{_legend}<svg viewBox="0 0 {W} {H}" width="100%" role="img" '
                     f'aria-label="Windows on the first two principal components, coloured by group" '
                     f'style="font: 12px system-ui, sans-serif; color: inherit;">{_dots}{_axes}</svg></div>')

    _panel = mo.vstack([
        mo.callout(mo.md("### 🧬 Clustered\n" + facts([
            ("Groups found", f"`{len(_clusters)}`"),
            ("Windows", f"`{len(_lab):,}` of {_W:,} bases, from {_r['sequences']:,} sequence(s)"),
            ("Ungrouped windows", f"`{_noise:,}`"),
        ])), kind="success"),
        mo.ui.table(_summary, selection=None, label="The groups and where their windows came from"),
        _chart,
        mo.md("A group made of windows from a single source is a clean bin; a source split across groups, or windows "
              "left ungrouped, often means a plasmid, a phage or a region of unusual composition. Try a longer window "
              "or a larger smallest group for steadier groups."),
    ])
    _panel
    return


@app.cell
def assembly_page_ui(example_picker, get_page, mo, set_page):
    mo.stop(get_page() != "assembly")
    btn_back_assembly = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    assembly_example = example_picker("assembly")

    _ui = mo.vstack([
        btn_back_assembly,
        mo.md("<br>\n\n# 🧩 Genome Assembly (De Bruijn)"),
        mo.callout(mo.md("This module takes the k-mers from Module 1, splits each one into a prefix (K-1) and a suffix (K-1), and maps the overlaps in memory, in Rust. The unbranched paths of that graph are the contigs."), kind="info"),
        assembly_example,
    ])
    _ui
    return assembly_example, btn_back_assembly


@app.cell
def assembly_page_form(assembly_example, example_for, files_with, get_page, make_form, mo):
    mo.stop(get_page() != "assembly")
    _ex, _about = example_for("assembly", None, assembly_example.value)
    # The .arrow files written by Module 1.
    form_assembly = make_form({
        "kmer_file": mo.ui.dropdown(options=files_with("arrow"), value=_ex.get("kmer_file"), label="🧬 K-mer file (.arrow, written by Module 1):", full_width=True),
    }, "🧩 Build the De Bruijn graph and assemble")
    mo.vstack([_x for _x in (_about, form_assembly) if _x is not None])
    return (form_assembly,)


@app.cell
async def assembly_page_backend(asyncio, beside, chosen, facts, figs, figure_panel, form_assembly, get_page, mmg, mo, os, pl, warn):
    mo.stop(get_page() != "assembly")
    mo.stop(form_assembly.value is None, mo.md("🧩 *Waiting for the assembly button...*"))

    _path = form_assembly.value.get("kmer_file")
    mo.stop(not chosen(_path), warn("Select a valid .arrow file."))

    try:
        # K comes from the file name (file_K55.arrow -> 55).
        _k_text = _path.split("_K")[-1].split(".arrow")[0]
        _k_size = int(_k_text)

        with mo.status.spinner(f"1/3 Reading the k-mers (K={_k_size})..."):
            await asyncio.sleep(0)
            _df = pl.read_ipc(_path)
            _kmer_bytes = _df.get_column("K-mers").item(0)

        with mo.status.spinner("2/3 Rust: building the graph and walking its paths..."):
            await asyncio.sleep(0)
            _nodes, _edges = await asyncio.to_thread(mmg.build_debruijn_graph, _kmer_bytes, _k_size)
            contig_list = await asyncio.to_thread(mmg.assemble_debruijn, _kmer_bytes, _k_size)

        with mo.status.spinner("3/3 Saving to disk and preparing the table..."):
            await asyncio.sleep(0)

            _out_name = beside(_path, f"MMGEN_Contigs_K{_k_size}.fasta")
            with open(_out_name, "w", encoding="utf-8") as fasta_file:
                for i, c in enumerate(contig_list):
                    fasta_file.write(f">NODE_{i+1}_length_{len(c)}\n{c}\n")

            _full_path = os.path.abspath(_out_name)
            _top_contigs = contig_list[:100]

            df_contigs = pl.DataFrame({
                "Contig ID": [f"NODE_{i+1}" for i in range(len(_top_contigs))],
                "Length (bp)": [len(c) for c in _top_contigs],
                "Start of sequence": [c[:50] + "..." for c in _top_contigs]
            })
            _table_ui = mo.ui.table(df_contigs, pagination=True, page_size=10)

        _panel = mo.vstack([
            mo.callout(mo.md(f"""
            ### 🗺️ Graph resolved successfully (128-bit k-mers)
            * **Nodes:** `{_nodes:,}` | **Edges:** `{_edges:,}`
            * **Contigs produced:** `{len(contig_list):,}`
            * **Longest contig:** `{len(contig_list[0]):,} bp`
            """), kind="success"),
            mo.callout(mo.md(f"💾 **File saved:** `{_full_path}`"), kind="info"),
            mo.md(f"#### 📋 Top {len(_top_contigs)} longest contigs (unitigs):"),
            _table_ui
        ])
        if hasattr(mmg, "unitig_graph"):
            with mo.status.spinner("Laying out the assembly graph..."):
                _g = await asyncio.to_thread(mmg.unitig_graph, _kmer_bytes, _k_size, 2000)
            _sm = _g["summary"]
            _facts_g = [("Unitigs / edges / components", f"`{_sm['nodes']:,}` / `{_sm['edges']:,}` / `{_sm['components']:,}`"
                         + (f" ({_sm['circular_components']} circular)" if _sm.get("circular_components") else "")),
                        ("N50 / longest", f"`{_sm['n50']:,}` / `{_sm['longest']:,}` bases"),
                        ("Tips / bubbles / repeats", f"`{_sm['tips']:,}` / `{_sm['bubbles']:,}` / `{_sm['repeats']:,}`")]
            if _sm.get("dropped_nodes"):
                _facts_g.append(("Not drawn", f"`{_sm['dropped_nodes']:,}` unitigs of the smallest components"))
            _panel = mo.vstack([_panel, mo.md("#### 🕸️ The assembly graph\n" + facts(_facts_g)), figure_panel(
                figs.assembly_graph(_g, title=f"Assembly graph (k = {_k_size})"), "assembly_graph",
                "As Bandage draws it. One long line: the genome assembled in one piece; a ring: a circular molecule. "
                "Tips are dead ends, usually sequencing errors; bubbles are two parallel paths (a SNP between two "
                "copies, or an error); repeats are where the path can go two ways — where an assembly breaks. Reads "
                "from both strands give two mirror-image components: these k-mers are not canonical.")])

    except Exception as e:
        _panel = mo.callout(mo.md(f"**🔥 ERROR:** `{str(e)}`"), kind="danger")

    _panel
    return


@app.cell
def variant_page_ui(EXAMPLES_DIR, NO_FILE, beside, example_picker, files_with, get_page, get_patient_msg, make_examples, mo, os, set_page, set_patient_msg):
    mo.stop(get_page() != "variant")

    def _back(_):
        set_patient_msg("")
        set_page("home")

    btn_back_variant = mo.ui.button(label="⬅️ Back to Menu", on_change=_back)

    # The simulated patient: a copy of a "healthy" genome with a few SNPs
    # planted at positions it reports, so the search below can be checked.
    _example_genome = f"{EXAMPLES_DIR}/toy_bacterium.fasta"
    _genomes = {f: f for f in files_with("plain_fasta") if f != NO_FILE and not f.endswith("mutated_patient.fasta")}
    _genomes.setdefault(_example_genome, _example_genome)
    patient_source = mo.ui.dropdown(options=_genomes, value=next(iter(_genomes)), label="Healthy genome to copy (FASTA):", full_width=True)
    patient_snps = mo.ui.number(start=1, stop=50, step=1, value=3, label="SNPs to plant:")

    def generate_patient():
        import random

        source = patient_source.value
        try:
            if source == _example_genome and not os.path.exists(source):
                make_examples()
            records, name, seq = [], None, []
            with open(source, encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    line = line.strip()
                    if line.startswith(">"):
                        if name is not None:
                            records.append([name, list("".join(seq))])
                        name, seq = line[1:], []
                    elif line:
                        seq.append(line)
            if name is not None:
                records.append([name, list("".join(seq))])
        except OSError as err:
            set_patient_msg(f"❌ Could not read `{source}`: {err}")
            return
        # Positions at least 100 bases from each other and from the ends of
        # their sequence: each SNP then changes exactly k k-mers, for any k up
        # to 63 (screen 1 asks for k later).
        n = int(patient_snps.value)
        rng = random.Random(2026)  # the same genome gives the same patient
        slots = [(i, pos) for i, (_, s) in enumerate(records) for pos in range(100, len(s) - 100, 100)
                 if s[pos].upper() in "ACGT"]
        if len(slots) < n:
            set_patient_msg(f"❌ `{source}` is too short for {n} SNPs at least 100 bases apart.")
            return
        planted = []
        for i, pos in sorted(rng.sample(slots, n)):
            ref = records[i][1][pos].upper()
            alt = rng.choice([b for b in "ACGT" if b != ref])
            records[i][1][pos] = alt
            planted.append(f"{records[i][0].split()[0]}:{pos + 1} {ref}>{alt}")
        out = beside(source, "mutated_patient.fasta")
        with open(out, "w", newline="\n") as fh:
            for name, s in records:
                fh.write(f">{name} (simulated patient)\n")
                text = "".join(s)
                for j in range(0, len(text), 80):
                    fh.write(text[j:j + 80] + "\n")
        set_patient_msg(f"✅ Done! `{out}` was written with {n} SNP(s): " + ", ".join(planted)
                        + ". Run both genomes through Module 1, then compare them below.")

    btn_generate = mo.ui.button(label="💉 Generate simulated patient", on_click=lambda _: generate_patient())
    variant_example = example_picker("variant")

    _ui = mo.vstack([
        btn_back_variant,
        mo.md("<br>\n\n# ⚠️ Variant Calling (SNPs)"),
        mo.callout(mo.md("To find mutations we compare a patient against a reference. The button below makes a simulated patient: a copy of the genome you pick, with SNPs planted at positions it tells you; run both through Module 1 to get their `.arrow` files. Or pick the example: a small genome and a patient, ready as k-mers."), kind="info"),
        mo.hstack([patient_source, patient_snps, btn_generate], justify="start", align="end", widths=[4, 1, 1.4]),
        # The button's feedback (an empty one gave "****", which Markdown draws
        # as a horizontal rule).
        mo.md(f"**{get_patient_msg()}**" if get_patient_msg() else ""),
        variant_example,
    ])
    _ui
    return btn_back_variant, btn_generate, generate_patient, patient_snps, patient_source, variant_example


@app.cell
def variant_page_form(example_for, files_with, get_page, make_form, mo, variant_example):
    mo.stop(get_page() != "variant")
    _ex, _about = example_for("variant", None, variant_example.value)
    _files = files_with("arrow")
    form_variant = make_form({
        "reference": mo.ui.dropdown(options=_files, value=_ex.get("reference"), label="🧬 1. Healthy reference genome (.arrow):", full_width=True),
        "patient": mo.ui.dropdown(options=_files, value=_ex.get("patient"), label="🩸 2. Patient sample (.arrow):", full_width=True),
        "filter_mb": mo.ui.slider(start=1, stop=128, step=1, value=_ex.get("filter_mb", 16), show_value=True, label="Bloom filter size (MB):"),
        "num_hashes": mo.ui.slider(start=1, stop=10, step=1, value=_ex.get("num_hashes", 5), show_value=True, label="Hashes per k-mer:"),
    }, "🔍 3. Run the search")
    mo.vstack([_x for _x in (_about, form_variant) if _x is not None])
    return (form_variant,)


@app.cell
async def variant_page_backend(asyncio, chosen, form_variant, get_page, mmg, mo, pl, warn):
    mo.stop(get_page() != "variant")
    mo.stop(form_variant.value is None, mo.md("🧪 *Waiting for the settings...*"))

    _ref = form_variant.value.get("reference")
    _pat = form_variant.value.get("patient")

    mo.stop(not chosen(_ref) or not chosen(_pat), warn("Select both files."))
    mo.stop(_ref == _pat, mo.callout(mo.md("⚠️ The reference and the patient cannot be the same file!"), kind="warn"))

    _k_text = _ref.split("_K")[-1].split(".arrow")[0]
    _k_size = int(_k_text)
    _mb = form_variant.value["filter_mb"]
    _hashes = form_variant.value["num_hashes"]
    # Where the Bloom filter really runs on this machine.
    _where = "GPU" if mmg.engine_info()["gpu"] else "CPU"

    try:
        with mo.status.spinner("1/3 Loading the k-mers..."):
            await asyncio.sleep(0)
            _ref_bytes = pl.read_ipc(_ref).get_column("K-mers").item(0)
            _pat_bytes = pl.read_ipc(_pat).get_column("K-mers").item(0)

        with mo.status.spinner(f"2/3 Building the reference Bloom filter on the {_where}..."):
            await asyncio.sleep(0)
            _bloom_bytes = await asyncio.to_thread(mmg.build_bloom_filter, _ref_bytes, _mb, _hashes)

        with mo.status.spinner(f"3/3 Looking for k-mers found only in the patient ({_where})..."):
            await asyncio.sleep(0)
            _mutation_bytes = await asyncio.to_thread(mmg.extract_mutations, _pat_bytes, _bloom_bytes, _hashes)
            _mutated_seqs = await asyncio.to_thread(mmg.decode_kmers_128, _mutation_bytes, _k_size)

        _total_mutations = len(_mutated_seqs)
        _preview = _mutated_seqs[:15] if _total_mutations > 0 else ["No mutation found."]

        _panel = mo.vstack([
            mo.callout(mo.md(f"""### 🏆 Genome scan complete (on the {_where})
The patient's {len(_pat_bytes) // 16:,} k-mers were checked against the reference filter.

* **Mutations detected (k-mers found only in the patient):** `{_total_mutations:,}`

**Sample of mutated sequences (k={_k_size}):**
```text
{chr(10).join(_preview)}
```"""), kind="success")
        ])

    except Exception as e:
        _panel = mo.callout(mo.md(f"**🔥 ERROR:** `{str(e)}`"), kind="danger")

    _panel
    return


@app.cell
def annotation_page_ui(example_picker, get_page, mo, set_page):
    mo.stop(get_page() != "annotation")
    btn_back_annotation = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    annotation_example = example_picker("annotation")

    _ui = mo.vstack([
        btn_back_annotation,
        mo.md("<br>\n\n# 📘 Functional Annotation (Translation)"),
        mo.callout(mo.md("DNA alone is only an instruction manual. To see whether a mutation matters, we translate the code three letters at a time (codons) and check whether the **protein** changed. Paste a sequence below to try it."), kind="info"),
        annotation_example,
    ])
    _ui
    return annotation_example, btn_back_annotation


@app.cell
def annotation_page_form(annotation_example, example_for, get_page, make_form, mo):
    mo.stop(get_page() != "annotation")
    _ex, _about = example_for("annotation", None, annotation_example.value)
    form_annotation = make_form({
        "dna": mo.ui.text_area(
            label="🧬 Paste the mutated DNA from Module 4 (or a whole contig):",
            # One of the mutated sequences found in the E. coli example.
            value=_ex.get("dna", "ATTTACCCTTGCTTTCACATAATCTGCCCGCCGATTTTGCTGCGTTGCGTAAATA"),
            full_width=True,
        ),
    }, "⚙️ Translate to protein")
    mo.vstack([_x for _x in (_about, form_annotation) if _x is not None])
    return (form_annotation,)


@app.cell
async def annotation_page_backend(asyncio, form_annotation, get_page, mmg, mo):
    mo.stop(get_page() != "annotation")
    mo.stop(form_annotation.value is None, mo.md("🧬 *Waiting for the DNA to translate...*"))

    _dna = form_annotation.value.get("dna", "").strip().upper()
    mo.stop(not _dna, mo.callout(mo.md("⚠️ Enter a valid DNA sequence."), kind="warn"))

    try:
        with mo.status.spinner("Translating DNA into protein in Rust (Rayon, multi-core)..."):
            await asyncio.sleep(0)
            # Rust takes a list of sequences and returns a list of proteins.
            _results = await asyncio.to_thread(mmg.translate_proteins, [_dna])
            _protein = _results[0]

        # Make the stop codons (*) stand out. (&#42; and not "*": with asterisks
        # around it Markdown did not make it bold and the screen showed "***".)
        _protein_html = _protein.replace("*", "<b style='color:red'>&#42;</b>")

        _panel = mo.vstack([
            mo.callout(mo.md(f"""
            ### 🧬 Translation complete!
            * **DNA length:** `{len(_dna)} bases`
            * **Protein length:** `{len(_protein)} amino acids`

            **Peptide sequence (amino acids):**
            > {_protein_html}

            *(Note: the red asterisk marks a stop codon. If a mutation created an early stop (a nonsense mutation), the protein is truncated and has probably lost its function.)*
            """), kind="success")
        ])
    except Exception as e:
        _panel = mo.callout(mo.md(f"**🔥 ERROR:** `{str(e)}`"), kind="danger")

    _panel
    return


@app.cell
def bloom_page_ui(example_picker, get_page, mo, set_page):
    mo.stop(get_page() != "bloom")
    btn_back_bloom = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    bloom_example = example_picker("bloom")

    _ui = mo.vstack([
        btn_back_bloom,
        mo.md("<br>\n\n# 🔍 Probabilistic Bloom Filters"),
        mo.callout(mo.md("A Bloom filter answers \"have I seen this k-mer?\" in a fixed amount of memory: each k-mer sets a few bits, chosen by its hashes. It never misses a k-mer it holds, but the fuller it gets, the more often it says yes to one it does not hold (a false positive)."), kind="info"),
        bloom_example,
    ])
    _ui
    return bloom_example, btn_back_bloom


@app.cell
def bloom_page_form(bloom_example, example_for, files_with, get_page, make_form, mo):
    mo.stop(get_page() != "bloom")
    _ex, _about = example_for("bloom", None, bloom_example.value)
    form_bloom = make_form({
        "kmer_file": mo.ui.dropdown(options=files_with("arrow"), value=_ex.get("kmer_file"), label="🧬 Genome (.arrow):", full_width=True),
        "filter_mb": mo.ui.slider(start=1, stop=128, step=1, value=_ex.get("filter_mb", 16), show_value=True, label="Filter size (MB):"),
        "num_hashes": mo.ui.slider(start=1, stop=10, step=1, value=_ex.get("num_hashes", 5), show_value=True, label="Number of hashes:"),
    }, "🔍 Build the Bloom filter")
    mo.vstack([_x for _x in (_about, form_bloom) if _x is not None])
    return (form_bloom,)


@app.cell
async def bloom_page_backend(asyncio, chosen, form_bloom, get_page, mmg, mo, np, pl, warn):
    mo.stop(get_page() != "bloom")
    mo.stop(form_bloom.value is None, mo.md("🔍 *Set the parameters above and click 'Build'...*"))

    _path = form_bloom.value.get("kmer_file")
    mo.stop(not chosen(_path), warn("Select a valid .arrow file."))

    _filter_mb = form_bloom.value["filter_mb"]
    _num_hashes = form_bloom.value["num_hashes"]
    # Where the filter really runs on this machine.
    _where = "GPU" if mmg.engine_info()["gpu"] else "CPU"

    with mo.status.spinner(f"🔥 Building the filter on the {_where}..."):
        await asyncio.sleep(0)
        _df = pl.read_ipc(_path)
        _kmer_bytes = _df.get_column("K-mers").item(0)
        _total_kmers = len(_kmer_bytes) // 16

        _bloom_bytes = await asyncio.to_thread(mmg.build_bloom_filter, _kmer_bytes, _filter_mb, _num_hashes)

        # Statistics
        _arr_8 = np.frombuffer(_bloom_bytes, dtype=np.uint8)
        _bits_set = np.unpackbits(_arr_8).sum()
        _bits_total = _filter_mb * 1024 * 1024 * 8
        _saturation = (_bits_set / _bits_total) * 100

    _panel = mo.callout(mo.md(f"""
    ### 🚀 Bloom filter built successfully!
    The {_where} ran **{_total_kmers * _num_hashes:,}** hash operations.
    * **Saturation:** `{_saturation:.2f}%`
    * **Bits set:** `{_bits_set:,}` of `{_bits_total:,}`
    """), kind="success")

    _panel
    return


@app.cell
def scaffold_page_ui(example_picker, get_page, mo, set_page):
    mo.stop(get_page() != "scaffold")
    btn_back_scaffold = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    scaffold_example = example_picker("scaffold")

    _ui = mo.vstack([
        btn_back_scaffold,
        mo.md("<br>\n\n# 🧬 Genome Scaffolding"),
        mo.callout(mo.md("The De Bruijn graph produced contigs. Scaffolding takes those blocks and stitches them together wherever their ends overlap exactly, trying to rebuild the original chromosome, circular ones included."), kind="info"),
        scaffold_example,
    ])
    _ui
    return btn_back_scaffold, scaffold_example


@app.cell
def scaffold_page_form(example_for, files_with, get_page, make_form, mo, scaffold_example):
    mo.stop(get_page() != "scaffold")
    _ex, _about = example_for("scaffold", None, scaffold_example.value)
    form_scaffold = make_form({
        "contigs_fasta": mo.ui.dropdown(options=files_with("plain_fasta"), value=_ex.get("contigs_fasta"), label="🧩 Contigs file (FASTA):", full_width=True),
        # With K=55 in the assembly, the contigs overlap by exactly K-1 (54).
        "overlap": mo.ui.slider(start=21, stop=100, step=1, value=_ex.get("overlap", 54), show_value=True, label="Exact overlap at the ends (bases):"),
    }, "🏗️ Close the genome (scaffolding)")
    mo.vstack([_x for _x in (_about, form_scaffold) if _x is not None])
    return (form_scaffold,)


@app.cell
async def scaffold_page_backend(asyncio, chosen, figs, figure_panel, form_scaffold, get_page, mmg, mo, warn):
    mo.stop(get_page() != "scaffold")
    mo.stop(form_scaffold.value is None, mo.md("🏗️ *Waiting for a FASTA file...*"))

    _path = form_scaffold.value.get("contigs_fasta")
    _overlap = form_scaffold.value.get("overlap")

    mo.stop(not chosen(_path), warn("Select a valid FASTA file produced by Module 3."))

    try:
        with mo.status.spinner("Reading the contigs..."):
            await asyncio.sleep(0)
            with open(_path, "r") as f:
                lines = f.readlines()

            contigs = []
            current = []
            for line in lines:
                if line.startswith(">"):
                    if current:
                        contigs.append("".join(current))
                        current = []
                else:
                    current.append(line.strip())
            if current:
                contigs.append("".join(current))

        with mo.status.spinner(f"Stitching ends with an exact overlap of {_overlap} bases, in Rust..."):
            await asyncio.sleep(0)
            scaffolds, circular = await asyncio.to_thread(mmg.build_scaffolds, contigs, _overlap)

        _panel = mo.vstack([
            mo.callout(mo.md(f"""
            ### 🏆 Scaffolding complete!
            The algorithm compared the ends of **{len(contigs):,}** contigs.

            * **Resulting scaffolds:** `{len(scaffolds):,}`
            * **Longest continuous sequence:** `{len(scaffolds[0]):,} bases`
            * **Circular genomes detected:** `{circular}` 🔄

            *(A circular genome means the tail of a sequence joined its own head exactly.)*
            """), kind="success")
        ])
        if hasattr(mmg, "overlap_graph"):
            _og = await asyncio.to_thread(mmg.overlap_graph, contigs, int(_overlap), 2000)
            _panel = mo.vstack([_panel, figure_panel(
                figs.assembly_graph(_og, title=f"Contigs joined by their ends (overlap {_overlap})"), "overlap_graph",
                "Each contig a line, coloured by the scaffold it ends up in; an arrow from the end of one contig to the start "
                "of the next. A loop of contigs is a circular molecule. Grey dashed links overlap but were not joined (the "
                "scaffolder needs an overlap of exactly the length given)."), mo.md(_og.get("text", "").replace("\n", "  \n"))])

    except Exception as e:
        _panel = mo.callout(mo.md(f"**🔥 ERROR:** `{str(e)}`"), kind="danger")

    _panel
    return


@app.cell
def vcf_page_ui(example_picker, get_page, mo, set_page):
    mo.stop(get_page() != "vcf")
    btn_back_vcf = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    vcf_example = example_picker("vcf")

    _ui = mo.vstack([
        btn_back_vcf,
        mo.md("<br>\n\n# 🗄️ 8. VCF Multimerge & Sorting"),
        mo.callout(mo.md("Bulk sorting of genetic variants. Select one or more `.vcf` or `.vcf.gz` files: the engine merges them all into a single output, with a stable sort, for files of any size (what does not fit in memory is sorted in runs on disk)."), kind="info"),
        vcf_example,
    ])
    _ui
    return btn_back_vcf, vcf_example


@app.cell
def vcf_page_form(example_for, files_with, get_page, make_form, mo, vcf_example):
    mo.stop(get_page() != "vcf")
    _ex, _about = example_for("vcf", None, vcf_example.value)
    form_vcf = make_form({
        "files": mo.ui.multiselect(options=files_with("vcf_merge"), value=_ex.get("files", []), label="🗄️ Files (VCF, VCF.GZ):", full_width=True),
        "sort_by": mo.ui.dropdown(options=["POS", "QUAL"], value=_ex.get("sort_by", "POS"), label="📊 Sort by:", full_width=True),
        "output_name": mo.ui.text(value=_ex.get("output_name", "Merged_Sorted.vcf"), label="💾 Output file name:", full_width=True),
    }, "🚀 Start the Multimerge engine")
    mo.vstack([_x for _x in (_about, form_vcf) if _x is not None])
    return (form_vcf,)


@app.cell
async def vcf_page_backend(asyncio, chosen, form_vcf, get_page, mmg, mo, os, warn):
    mo.stop(get_page() != "vcf")
    mo.stop(form_vcf.value is None, mo.md("⏳ *Waiting for the settings...*"))

    _paths = [_p for _p in form_vcf.value.get("files", []) if chosen(_p)]
    _sort_by = form_vcf.value.get("sort_by")
    _output = (form_vcf.value.get("output_name") or "").strip() or "Merged_Sorted.vcf"

    mo.stop(not _paths, warn("Select at least one file."))
    # The result of an earlier run shows up in the list (it is a .vcf in the same
    # folder): reading and writing the same file would destroy the data.
    mo.stop(
        os.path.abspath(_output) in [os.path.abspath(_c) for _c in _paths],
        mo.callout(mo.md("⚠️ The output file cannot be one of the input files. Choose another name."), kind="warn")
    )

    try:
        with mo.status.spinner(f"Merging and sorting {len(_paths)} file(s) in the Multimerge engine..."):
            await asyncio.sleep(0)

            # ONE call with ALL the files: the engine merges and sorts them with a
            # stable sort. (It used to be one CLI call per file, all with the same
            # output: each erased the previous one and only the last file
            # survived.) It runs in a thread (Rust releases the GIL) so the
            # interface does not freeze.
            _variants = await asyncio.to_thread(mmg.multimerge_vcf, _paths, _sort_by, _output)

            _full_path = os.path.abspath(_output)

        _panel = mo.vstack([
            mo.callout(mo.md(f"""
            ### ✅ VCF Multimerge complete!
            The engine read, merged and sorted the files.

            * **Files merged:** `{len(_paths)}`
            * **Variants written:** `{_variants:,}`
            * **Sorted by:** `{_sort_by}`
            * **Saved to:** `{_full_path}`
            """), kind="success")
        ])

    except Exception as e:
        _panel = mo.callout(mo.md(f"**🔥 MULTIMERGE ENGINE ERROR:** `{str(e)}`"), kind="danger")

    _panel
    return


# ---------------------------------------------------------------------------
# Screens 9-14: the file tools. Each has three cells: the screen (title and the
# choice of task), the form (rebuilt when the task changes or the file list is
# refreshed) and the run. Every result ends with the equivalent Python call.
# ---------------------------------------------------------------------------


@app.cell
def sort_page_ui(get_page, mo, set_page):
    mo.stop(get_page() != "sort")
    btn_back_sort = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_sort = mo.ui.button(label="🔄 Refresh the file list")
    sort_task = mo.ui.radio(
        options=["Sort & merge", "Index & query a region", "Interval operations (BED)"],
        value="Sort & merge", inline=True, label="**What to do:**",
    )
    interval_op = mo.ui.dropdown(
        options={
            "merge: join overlapping or nearby intervals (bedtools merge)": "merge",
            "intersect: A against B (bedtools intersect)": "intersect",
            "subtract: remove from A what B covers (bedtools subtract)": "subtract",
            "complement: the stretches no interval covers (bedtools complement)": "complement",
            "closest: the nearest B to each A (bedtools closest -d)": "closest",
            "coverage: how much of each A is covered by B (bedtools coverage)": "coverage",
            "slop: grow every interval (bedtools slop)": "slop",
        },
        value="merge: join overlapping or nearby intervals (bedtools merge)", label="Operation:", full_width=True,
    )

    _ui = mo.vstack([
        btn_back_sort,
        mo.md("<br>\n\n# 📏 9. Sorting, Indexing & Intervals"),
        mo.callout(mo.md("Sort and merge VCF, BED, GFF/GTF and SAM/BAM files of any size (what does not fit in memory is sorted in runs on disk); build the `.tbi` / `.bai` index of a sorted file and read only the region you ask for; or run the interval operations of bedtools on BED files, which do not need to be sorted. Results are written next to their input."), kind="info"),
        mo.hstack([sort_task, btn_refresh_sort], justify="space-between", align="center"),
    ])
    _ui
    return btn_back_sort, btn_refresh_sort, interval_op, sort_task


@app.cell
def sort_page_example(example_picker, get_page, interval_op, mo, sort_task):
    mo.stop(get_page() != "sort")
    # The operation is picked here, above the examples of that operation.
    _intervals = sort_task.value == "Interval operations (BED)"
    sort_example = example_picker("sort", interval_op.value if _intervals else sort_task.value)
    mo.vstack([interval_op, sort_example] if _intervals else [sort_example])
    return (sort_example,)


@app.cell
def sort_page_form(btn_refresh_sort, example_for, files_with, get_page, interval_op, label_of, make_form, mo, real_files, sort_example, sort_task):
    mo.stop(get_page() != "sort")
    btn_refresh_sort  # a click re-runs this cell, which lists the folder again
    _task = sort_task.value
    _ex, _about = example_for("sort", interval_op.value if _task == "Interval operations (BED)" else _task, sort_example.value)

    if _task == "Sort & merge":
        _orders = {
            "default: the header for VCF and SAM, lexical for BED and GFF": "default",
            "lexical: chr1, chr10, chr2 (bedtools)": "lexical",
            "natural: chr1, chr2, chr10": "natural",
            "header: the ##contig / @SQ lines (bcftools, samtools)": "header",
        }
        for _fai in real_files("genome"):
            _orders[f"the order of {_fai}"] = _fai
        _by = {"position": "position", "quality (VCF, highest first)": "quality", "read name (SAM)": "name"}
        form_sort = make_form({
            "files": mo.ui.multiselect(options=files_with("vcf", "bed", "gff", "sam"), value=_ex.get("files", []), label="Files (several are merged into one output):", full_width=True),
            "by": mo.ui.dropdown(options=_by, value=label_of(_by, _ex.get("by"), "position"), label="Sort by:", full_width=True),
            "order": mo.ui.dropdown(options=_orders, value=label_of(_orders, _ex.get("order"), "default: the header for VCF and SAM, lexical for BED and GFF"), label="Chromosome order:", full_width=True),
            "memory": mo.ui.number(start=0, stop=1_048_576, step=64, value=_ex.get("memory", 0), label="Memory budget in MB (0 = half of the free RAM):"),
            "index": mo.ui.checkbox(value=_ex.get("index", False), label="Also write the index (.tbi for a .gz output, .bai for a .bam output)"),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: <first file>.sorted.<same extension>", label="Output file:", full_width=True),
        }, "🚀 Sort")
    elif _task == "Index & query a region":
        _indexed = files_with("indexed")
        form_sort = make_form({
            "file": mo.ui.dropdown(options=_indexed, value=_ex.get("file"), label="Sorted file (a .gz written as BGZF, or a .bam):", full_width=True),
            "regions": mo.ui.text(value=_ex.get("regions", ""), placeholder="chr1:1000000-2000000, chr2", label="Region(s), 1-based, end included (blank: only build the index):", full_width=True),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: show the records here", label="Save the records to:", full_width=True),
        }, "🔎 Index and query")
    else:
        _op = interval_op.value
        _bed = files_with("bed")
        _fields = {"a": mo.ui.multiselect(options=_bed, value=_ex.get("a", []), label="A: BED file(s):", full_width=True)}
        if _op in ("intersect", "subtract", "closest", "coverage"):
            _fields["b"] = mo.ui.multiselect(options=_bed, value=_ex.get("b", []), label="B: BED file(s):", full_width=True)
        if _op == "intersect":
            _modes = {
                "regions: the part of each A that overlaps B": "regions",
                "unique: the lines of A that touch B (-u)": "unique",
                "absent: the lines of A that do not touch B (-v)": "absent",
            }
            _fields["mode"] = mo.ui.dropdown(options=_modes, value=label_of(_modes, _ex.get("mode"), "regions: the part of each A that overlaps B"), label="Write:", full_width=True)
        if _op == "subtract":
            _fields["whole"] = mo.ui.checkbox(value=_ex.get("whole", False), label="Drop every A that B touches, instead of clipping it (-A)")
        if _op == "merge":
            _fields["distance"] = mo.ui.number(start=0, stop=10**9, step=1, value=_ex.get("distance", 0), label="Also join intervals up to this many bases apart (-d):")
        if _op in ("complement", "slop"):
            _fields["genome"] = mo.ui.dropdown(options=files_with("genome"), value=_ex.get("genome"), label="Genome file (.fai, or chromosome sizes: name and length per line):", full_width=True)
        if _op == "slop":
            _fields["left"] = mo.ui.number(start=0, stop=10**9, step=1, value=_ex.get("left", 100), label="Bases to add before the start:")
            _fields["right"] = mo.ui.number(start=0, stop=10**9, step=1, value=_ex.get("right", 100), label="Bases to add after the end:")
        _fields["output"] = mo.ui.text(value=_ex.get("output", ""), placeholder=f"blank: <first A file>.{_op}.<same extension>", label="Output file:", full_width=True)
        form_sort = make_form(_fields, f"🚀 Run {_op}")

    mo.vstack([_x for _x in (_about, form_sort) if _x is not None])
    return (form_sort,)


@app.cell
async def sort_page_backend(asyncio, chosen, clashes, code_panel, error_panel, facts, form_sort, get_page, interval_op, mmg, mo, os, output_name, preview_panel, python_call, size_of, skipped_panel, sort_task, split_name, warn):
    mo.stop(get_page() != "sort")
    mo.stop(form_sort.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_sort.value
    _task = sort_task.value

    if _task == "Sort & merge":
        _paths = [_p for _p in _v["files"] if chosen(_p)]
        mo.stop(not _paths, warn("Select at least one file."))
        _ext = None
        if _v["index"]:
            # An index needs a BGZF .gz or a .bam output.
            _fmt = (split_name(_paths[0])[1].lower().split(".") + ["", ""])[1]
            _ext = ".bam" if _fmt in ("sam", "bam") else f".{_fmt}.gz" if _fmt else None
        _out = output_name(_v["output"], _paths[0], "sorted", _ext)
        mo.stop(clashes(_out, _paths), warn("The output file cannot be one of the input files: choose another name."))
        _kw = {}
        if _v["by"] != "position":
            _kw["by"] = _v["by"]
        if _v["order"] != "default":
            _kw["chrom_order"] = _v["order"]
        if int(_v["memory"] or 0) > 0:
            _kw["memory_mb"] = int(_v["memory"])
        if _v["index"]:
            _kw["index"] = True
        try:
            with mo.status.spinner(f"Sorting {len(_paths)} file(s)..."):
                _records, _runs, _skipped = await asyncio.to_thread(mmg.sort_genomic, _paths, _out, **_kw)
            _rows = [
                ("Records written", f"`{_records:,}`"),
                ("Runs on disk", f"`{_runs}`" + (" (everything fit in memory)" if _runs == 0 else "")),
                ("Malformed lines skipped", f"`{_skipped:,}`"),
                ("Saved to", f"`{os.path.abspath(_out)}` ({size_of(_out)})"),
            ]
            if _v["index"]:
                _rows.append(("Index", f"`{_out}{'.bai' if _out.lower().endswith('.bam') else '.tbi'}`"))
            _panel = mo.vstack([_x for _x in (
                mo.callout(mo.md(("### ✅ Sorted\n" if _records else "### ⚠️ Nothing to sort\n") + facts(_rows)),
                           kind="success" if _records else "warn"),
                skipped_panel(_paths, _skipped, _records),
                preview_panel(_out),
                code_panel(f"records, runs, skipped = {python_call('sort_genomic', _paths, _out, **_kw)}"),
            ) if _x is not None])
        except Exception as _err:
            _panel = error_panel(_err)

    elif _task == "Index & query a region":
        _path = _v["file"]
        mo.stop(not chosen(_path), warn("Select a file. Only a sorted `.gz` written as BGZF (every `.gz` this suite writes is one; so is `bgzip`'s) or a `.bam` can be indexed: sort it first with **Sort & merge**."))
        _regions = [_r for _r in (_v["regions"] or "").replace(";", " ").replace(",", " ").split() if _r]
        _out = (_v["output"] or "").strip()
        mo.stop(bool(_out) and not _regions, warn("Type the region(s) whose records should be saved."))
        mo.stop(bool(_out) and clashes(_out, [_path]), warn("The output file cannot be the input file: choose another name."))
        _index = _path + (".bai" if _path.lower().endswith(".bam") else ".tbi")
        _arg = _regions[0] if len(_regions) == 1 else _regions
        _calls, _items = [], []
        try:
            if not os.path.exists(_index) or os.path.getmtime(_index) < os.path.getmtime(_path):
                with mo.status.spinner("Building the index..."):
                    _ix = await asyncio.to_thread(mmg.index_file, _path)
                _items.append(mo.callout(mo.md("### ✅ Index written\n" + facts([
                    ("Index", f"`{_ix['index']}`"),
                    ("Sequences", f"`{_ix['sequences']:,}`"),
                    ("Records", f"`{_ix['records']:,}`"),
                    ("Lines skipped", f"`{_ix['skipped']:,}`"),
                ])), kind="success"))
                if _ix["skipped"]:
                    _items.append(skipped_panel(_path, _ix["skipped"], _ix["records"]))
                _calls.append(python_call("index_file", _path))
            else:
                _items.append(mo.md(f"Using the index `{_index}`, which is newer than the file."))
            if _regions:
                with mo.status.spinner(f"Reading {', '.join(_regions)}..."):
                    if _out:
                        _n = await asyncio.to_thread(mmg.query_region, _path, _arg, _out)
                    else:
                        _records = await asyncio.to_thread(mmg.query_region, _path, _arg)
                        _n = len(_records)
                _items.append(mo.callout(mo.md(f"### ✅ {_n:,} record(s) in {', '.join(_regions)}"
                                               + (f"\n* **Saved to:** `{os.path.abspath(_out)}`" if _out else "")), kind="success"))
                if _out:
                    _items.append(preview_panel(_out))
                    _calls.append(python_call("query_region", _path, _arg, _out))
                else:
                    if _records:
                        _items.append(mo.md(f"**{'The first 50' if _n > 50 else 'The'} record(s):**"))
                        _items.append(mo.plain_text("\n".join(_l if len(_l) <= 160 else _l[:157] + "..." for _l in _records[:50])))
                    _calls.append(f"records = {python_call('query_region', _path, _arg)}")
            elif not _calls:
                _items.append(mo.md("The index is up to date: type a region to read it."))
            if _calls:
                _items.append(code_panel(*_calls))
            _panel = mo.vstack(_items)
        except Exception as _err:
            _panel = error_panel(_err)

    else:
        _op = interval_op.value
        _a = [_p for _p in _v["a"] if chosen(_p)]
        _b = [_p for _p in _v.get("b", []) if chosen(_p)]
        _genome = _v.get("genome")
        mo.stop(not _a, warn("Select at least one A file."))
        mo.stop("b" in _v and not _b, warn("Select at least one B file."))
        mo.stop("genome" in _v and not chosen(_genome), warn("Select the genome file: a `.fai` index, or a chromosome-sizes file (`.genome` or `.sizes`: name and length per line, tab-separated)."))
        _out = output_name(_v["output"], _a[0], _op)
        mo.stop(clashes(_out, _a + _b + ([_genome] if _genome else [])), warn("The output file cannot be one of the input files: choose another name."))
        _notes = {
            "closest": "Each line is the line of A, the line of its closest B, and the distance (0 = they overlap; two intervals that touch are at distance 1). Ties are all written.",
            "coverage": "Each line is the line of A plus four columns: how many B overlap it, how many of its bases they cover, its length, and the covered fraction.",
            "complement": "The output follows the order of the genome file.",
        }
        try:
            with mo.status.spinner(f"Running {_op}..."):
                _read = None
                if _op == "merge":
                    _kw = {"distance": int(_v["distance"])} if int(_v["distance"] or 0) > 0 else {}
                    _read, _n = await asyncio.to_thread(mmg.bed_merge, _a, _out, **_kw)
                    _call = f"read, written = {python_call('bed_merge', _a, _out, **_kw)}"
                elif _op == "intersect":
                    _kw = {"mode": _v["mode"]} if _v["mode"] != "regions" else {}
                    _n = await asyncio.to_thread(mmg.bed_intersect, _a, _b, _out, **_kw)
                    _call = python_call("bed_intersect", _a, _b, _out, **_kw)
                elif _op == "subtract":
                    _kw = {"remove_whole": True} if _v["whole"] else {}
                    _n = await asyncio.to_thread(mmg.bed_subtract, _a, _b, _out, **_kw)
                    _call = python_call("bed_subtract", _a, _b, _out, **_kw)
                elif _op == "complement":
                    _n = await asyncio.to_thread(mmg.bed_complement, _a, _genome, _out)
                    _call = python_call("bed_complement", _a, _genome, _out)
                elif _op == "closest":
                    _n = await asyncio.to_thread(mmg.bed_closest, _a, _b, _out)
                    _call = python_call("bed_closest", _a, _b, _out)
                elif _op == "coverage":
                    _n = await asyncio.to_thread(mmg.bed_coverage, _a, _b, _out)
                    _call = python_call("bed_coverage", _a, _b, _out)
                else:
                    _kw = {"left": int(_v["left"]), "right": int(_v["right"])}
                    _n = await asyncio.to_thread(mmg.bed_slop, _a, _genome, _out, **_kw)
                    _call = python_call("bed_slop", _a, _genome, _out, **_kw)
            _rows = ([("Intervals read", f"`{_read:,}`")] if _read is not None else []) + [
                ("Lines written", f"`{_n:,}`"),
                ("Saved to", f"`{os.path.abspath(_out)}` ({size_of(_out)})"),
            ]
            _panel = mo.vstack([
                mo.callout(mo.md(f"### ✅ {_op.capitalize()} done\n" + facts(_rows) + (f"\n\n{_notes[_op]}" if _op in _notes else "")), kind="success"),
                preview_panel(_out),
                code_panel(_call),
            ])
        except Exception as _err:
            _panel = error_panel(_err)

    _panel
    return


@app.cell
def validate_page_ui(get_page, mo, set_page):
    mo.stop(get_page() != "validate")
    btn_back_validate = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_validate = mo.ui.button(label="🔄 Refresh the file list")
    validate_task = mo.ui.radio(
        options=["Validate a file", "Rename chromosomes (chr1 ↔ 1)", "Compare reference genomes"],
        value="Validate a file", inline=True, label="**What to do:**",
    )
    _ui = mo.vstack([
        btn_back_validate,
        mo.md("<br>\n\n# 🩺 10. Validation & Chromosome Names"),
        mo.callout(mo.md("Find out what is wrong with a file before another tool fails on it: the validator reads VCF, BED, GFF/GTF, SAM/BAM, FASTA and FASTQ files of any size and gives, for each kind of problem, the line where it first shows up, how often it repeats and how to fix it. Rename the chromosomes of a whole file (`chr1` ↔ `1`), and check that two files are on the same reference genome (GRCh37 or GRCh38) before combining them."), kind="info"),
        mo.hstack([validate_task, btn_refresh_validate], justify="space-between", align="center"),
    ])
    _ui
    return btn_back_validate, btn_refresh_validate, validate_task


@app.cell
def validate_page_example(example_picker, get_page, mo, validate_task):
    mo.stop(get_page() != "validate")
    validate_example = example_picker("validate", validate_task.value)
    validate_example
    return (validate_example,)


@app.cell
def validate_page_form(btn_refresh_validate, example_for, files_with, get_page, label_of, make_form, mo, real_files, validate_example, validate_task):
    mo.stop(get_page() != "validate")
    btn_refresh_validate  # a click re-runs this cell, which lists the folder again
    _task = validate_task.value
    _ex, _about = example_for("validate", _task, validate_example.value)

    if _task == "Validate a file":
        form_validate = make_form({
            "file": mo.ui.dropdown(options=files_with("vcf", "bed", "gff", "sam", "fasta", "fastq"), value=_ex.get("file"), label="File:", full_width=True),
            "max_lines": mo.ui.number(start=0, stop=10**12, step=1000, value=_ex.get("max_lines", 0), label="Check only the first N lines (0 = the whole file):"),
        }, "🩺 Validate")
    elif _task.startswith("Rename"):
        _maps = {"(none)": "none"}
        for _m in real_files("mapping"):
            _maps[_m] = _m
        _styles = {"chr1, chrX, chrM (UCSC)": "chr", "1, X, MT (Ensembl, NCBI)": "no_chr"}
        form_validate = make_form({
            "file": mo.ui.dropdown(options=files_with("vcf", "bed", "gff", "sam"), value=_ex.get("file"), label="File:", full_width=True),
            "style": mo.ui.dropdown(options=_styles, value=label_of(_styles, _ex.get("style"), "chr1, chrX, chrM (UCSC)"), label="Rename to:", full_width=True),
            "mapping": mo.ui.dropdown(options=_maps, value=label_of(_maps, _ex.get("mapping"), "(none)"), label="Mapping file for other names (two columns: name, new name):", full_width=True),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: <file>.chr.<same extension> (or .nochr)", label="Output file:", full_width=True),
        }, "🔁 Rename")
    else:
        _refs = files_with("vcf", "sam", "gff", "genome")
        form_validate = make_form({
            "a": mo.ui.dropdown(options=_refs, value=_ex.get("a"), label="File A:", full_width=True),
            "b": mo.ui.dropdown(options=_refs, value=_ex.get("b"), label="File B:", full_width=True),
        }, "⚖️ Compare")

    mo.vstack([_x for _x in (_about, form_validate) if _x is not None])
    return (form_validate,)


@app.cell
async def validate_page_backend(asyncio, chosen, clashes, code_panel, error_panel, facts, form_validate, get_page, mmg, mo, os, output_name, preview_panel, python_call, size_of, validate_task, warn):
    mo.stop(get_page() != "validate")
    mo.stop(form_validate.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_validate.value
    _task = validate_task.value

    if _task == "Validate a file":
        _path = _v["file"]
        mo.stop(not chosen(_path), warn("Select a file."))
        _lines = int(_v["max_lines"] or 0)
        _kw = {"max_lines": _lines} if _lines > 0 else {}
        try:
            with mo.status.spinner(f"Reading {_path}..."):
                _r = await asyncio.to_thread(mmg.validate_file, _path, **_kw)
            if _r["valid"] and not _r["warnings"]:
                _head, _kind = "### ✅ No problem found", "success"
            elif _r["valid"]:
                _head, _kind = f"### ✅ Valid, with {_r['warnings']} warning(s)", "warn"
            else:
                _head, _kind = f"### ❌ {_r['errors']} error(s), {_r['warnings']} warning(s)", "danger"
            _naming = {"chr": "UCSC (`chr1`)", "no_chr": "Ensembl / NCBI (`1`)", "mixed": "**mixed**: `chr1` and `1` in the same file",
                       "other": "other names", "-": "none"}.get(_r["naming"], _r["naming"])
            _rows = [
                ("Format", f"`{_r['format'].upper()}`"),
                ("Lines read", f"`{_r['lines']:,}`" + ("" if _r["complete"] else " (only the beginning of the file)")),
                ("Records", f"`{_r['records']:,}`"),
                ("Sequences", f"`{_r['sequences']:,}`"),
                ("Chromosome names", _naming),
            ]
            if _r["format"] in ("vcf", "bed", "gff", "gtf", "sam"):
                _rows.append(("Sorted by position", "yes" if _r["sorted"] else "**no**"))
            _icons = {"error": "❌", "warning": "⚠️", "info": "ℹ️"}
            _problems = []
            for _p in _r["problems"]:
                _times = f" *({_p['occurrences']:,} times)*" if _p["occurrences"] > 1 else ""
                _problems.append(f"* {_icons.get(_p['severity'], '•')} **line {_p['line']:,}**: {_p['message']}{_times}"
                                 + (f"  \n  ↳ fix: {_p['fix']}" if _p["fix"] else ""))
            _next = []
            if any("normalize_chromosomes" in (_p["fix"] or "") for _p in _r["problems"]):
                _next.append("rename the chromosomes with **Rename chromosomes** above")
            if any("sort_genomic" in (_p["fix"] or "") for _p in _r["problems"]):
                _next.append("sort the file on screen **9. Sorting, Indexing & Intervals**")
            _items = [mo.callout(mo.md(_head + "\n" + facts(_rows)), kind=_kind)]
            if _problems:
                _items.append(mo.md("**What the validator found** (the line of the first occurrence of each kind of problem):\n\n" + "\n".join(_problems)))
            if _next:
                _items.append(mo.md("**Next step:** " + "; ".join(_next) + "."))
            _items.append(code_panel(f"report = {python_call('validate_file', _path, **_kw)}", 'print(report["text"])'))
            _panel = mo.vstack(_items)
        except Exception as _err:
            _panel = error_panel(_err)

    elif _task.startswith("Rename"):
        _path = _v["file"]
        mo.stop(not chosen(_path), warn("Select a file."))
        _style = _v["style"]
        _map = None if _v["mapping"] == "none" else _v["mapping"]
        _out = output_name(_v["output"], _path, "chr" if _style == "chr" else "nochr")
        mo.stop(clashes(_out, [_path, _map]), warn("The output file cannot be one of the input files: choose another name."))
        _kw = {"style": _style}
        if _map:
            _kw["mapping"] = _map
        try:
            with mo.status.spinner(f"Rewriting {_path}..."):
                _r = await asyncio.to_thread(mmg.normalize_chromosomes, _path, _out, **_kw)
            _pairs = ", ".join(f"`{_a}` → `{_b}`" for _a, _b in _r["converted"][:40]) + (" ..." if len(_r["converted"]) > 40 else "")
            _kept = ", ".join(f"`{_k}`" for _k in _r["kept"][:40]) + (" ..." if len(_r["kept"]) > 40 else "")
            _rows = [
                ("Records", f"`{_r['records']:,}`"),
                ("Records changed", f"`{_r['changed']:,}`"),
                ("Names converted", _pairs or "none"),
            ]
            if _r["kept"]:
                _rows.append(("Names kept as they were", _kept + " (scaffolds and accessions need a mapping file)"))
            _rows.append(("Saved to", f"`{os.path.abspath(_out)}` ({size_of(_out)})"))
            _panel = mo.vstack([
                mo.callout(mo.md("### ✅ Chromosomes renamed\n" + facts(_rows)), kind="success"),
                preview_panel(_out),
                code_panel(python_call("normalize_chromosomes", _path, _out, **_kw)),
            ])
        except Exception as _err:
            _panel = error_panel(_err)

    else:
        _a, _b = _v["a"], _v["b"]
        mo.stop(not chosen(_a) or not chosen(_b), warn("Select the two files."))
        mo.stop(_a == _b, warn("Select two different files."))
        try:
            with mo.status.spinner("Comparing the sequence lengths..."):
                _c = await asyncio.to_thread(mmg.compare_references, _a, _b)
            _kind = {True: "success", False: "danger"}.get(_c["compatible"], "warn")
            _icon = {True: "✅", False: "❌"}.get(_c["compatible"], "❔")
            _different = ", ".join(f"`{_d[0]}` ({_d[1]:,} × {_d[2]:,})" for _d in _c["different"][:10]) + (" ..." if len(_c["different"]) > 10 else "")
            _rows = [
                ("A", f"`{_a}`: " + (_c["assembly_a"] or "assembly not recognized") + (f" (declares {_c['declared_a']})" if _c["declared_a"] else "")),
                ("B", f"`{_b}`: " + (_c["assembly_b"] or "assembly not recognized") + (f" (declares {_c['declared_b']})" if _c["declared_b"] else "")),
                ("Sequences with the same length", f"`{_c['matching']:,}`"),
                ("With different lengths", f"`{len(_c['different']):,}`" + (f": {_different}" if _different else "")),
                ("Only in A / only in B", f"`{_c['only_in_a']:,}` / `{_c['only_in_b']:,}`"),
            ]
            _items = [mo.callout(mo.md(f"### {_icon} {_c['verdict']}\n" + facts(_rows)), kind=_kind)]
            if _c["different_naming"]:
                _items.append(mo.md("The two files name the chromosomes differently (`chr1` in one, `1` in the other); they were matched anyway. Before combining them, rename one with **Rename chromosomes** above."))
            _items.append(code_panel(f"result = {python_call('compare_references', _a, _b)}", 'print(result["verdict"])'))
            _panel = mo.vstack(_items)
        except Exception as _err:
            _panel = error_panel(_err)

    _panel
    return


@app.cell
def vcftools_page_ui(get_page, mo, set_page):
    mo.stop(get_page() != "vcftools")
    btn_back_vcftools = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_vcftools = mo.ui.button(label="🔄 Refresh the file list")
    vcftools_task = mo.ui.radio(options=["Statistics", "Filter", "Normalize"], value="Statistics", inline=True, label="**What to do:**")
    _ui = mo.vstack([
        btn_back_vcftools,
        mo.md("<br>\n\n# 🧪 11. VCF Statistics, Filtering & Normalization"),
        mo.callout(mo.md("Count what a VCF holds; keep the records that match an expression, as in bcftools; and normalize the variants: split multiallelic sites, trim them, and, given the genome, left-align the indels as `bcftools norm -f` does. The file is read as a stream, so its size does not matter."), kind="info"),
        mo.hstack([vcftools_task, btn_refresh_vcftools], justify="space-between", align="center"),
    ])
    _ui
    return btn_back_vcftools, btn_refresh_vcftools, vcftools_task


@app.cell
def vcftools_page_example(example_picker, get_page, mo, vcftools_task):
    mo.stop(get_page() != "vcftools")
    vcftools_example = example_picker("vcftools", vcftools_task.value)
    vcftools_example
    return (vcftools_example,)


@app.cell
def vcftools_page_form(btn_refresh_vcftools, example_for, files_with, get_page, label_of, make_form, mo, real_files, vcftools_example, vcftools_task):
    mo.stop(get_page() != "vcftools")
    btn_refresh_vcftools  # a click re-runs this cell, which lists the folder again
    _task = vcftools_task.value
    _ex, _about = example_for("vcftools", _task, vcftools_example.value)
    _vcfs = files_with("vcf")
    _above = [_about] if _about is not None else []

    if _task == "Statistics":
        form_vcftools = make_form({"file": mo.ui.dropdown(options=_vcfs, value=_ex.get("file"), label="VCF file:", full_width=True)}, "📊 Count")
    elif _task == "Filter":
        _above += [mo.accordion({"How to write an expression": mo.md(
            "Fields: `QUAL`, `POS`, `N_ALT` (numbers); `CHROM`, `ID`, `REF`, `ALT`, `FILTER`, `TYPE` (text); `INFO/KEY` "
            "(alone, it tests that the key is present). Operators: `==` `!=` `<` `<=` `>` `>=`, `&&` (and), `||` (or), `!` (not) "
            "and parentheses. A comparison with a missing value is false.\n\n"
            "`TYPE` is `\"snp\"`, `\"mnp\"`, `\"ins\"`, `\"del\"`, `\"indel\"`, `\"other\"` or `\"ref\"`, and works as in bcftools: "
            "`TYPE == \"snp\"` keeps a site when **every** ALT allele is a SNP, `TYPE ~ \"snp\"` when **at least one** is. "
            "`FILTER ~ \"q10\"` keeps a site whose filters include q10; `FILTER == \"q10\"` asks for exactly q10.\n\n"
            "Examples: `QUAL >= 30 && FILTER == \"PASS\"` · `TYPE == \"snp\" && INFO/DP > 10` · `!(INFO/AF < 0.01)`"
        )})]
        form_vcftools = make_form({
            "file": mo.ui.dropdown(options=_vcfs, value=_ex.get("file"), label="VCF file:", full_width=True),
            "expression": mo.ui.text(value=_ex.get("expression", ""), placeholder='QUAL >= 30 && FILTER == "PASS"', label="Expression (blank: no expression):", full_width=True),
            "regions": mo.ui.text(value=_ex.get("regions", ""), placeholder="chr1:1000000-2000000, chrX", label="Regions (blank: the whole genome):", full_width=True),
            "samples": mo.ui.text(value=_ex.get("samples", ""), placeholder="sample1, sample2", label="Samples to keep, in this order (blank: all of them):", full_width=True),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: <file>.filtered.<same extension>", label="Output file:", full_width=True),
        }, "🔍 Filter")
    else:
        _refs = {"(none: split and trim only, no left-alignment)": "none"}
        for _f in real_files("fasta"):
            _refs[_f] = _f
        form_vcftools = make_form({
            "file": mo.ui.dropdown(options=_vcfs, value=_ex.get("file"), label="VCF file:", full_width=True),
            "split": mo.ui.checkbox(value=_ex.get("split", True), label="Split multiallelic sites into one record per allele"),
            "reference": mo.ui.dropdown(options=_refs, value=label_of(_refs, _ex.get("reference"), "(none: split and trim only, no left-alignment)"),
                                        label="Reference genome (FASTA), to check REF and left-align the indels:", full_width=True),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: <file>.norm.<same extension>", label="Output file:", full_width=True),
        }, "🧹 Normalize")

    mo.vstack(_above + [form_vcftools])
    return (form_vcftools,)


@app.cell
async def vcftools_page_backend(asyncio, chosen, clashes, code_panel, error_panel, facts, form_vcftools, get_page, mmg, mo, os, output_name, preview_panel, python_call, size_of, skipped_panel, vcftools_task, warn):
    mo.stop(get_page() != "vcftools")
    mo.stop(form_vcftools.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_vcftools.value
    _task = vcftools_task.value
    _path = _v["file"]
    mo.stop(not chosen(_path), warn("Select a VCF file."))

    if _task == "Statistics":
        try:
            with mo.status.spinner(f"Reading {_path}..."):
                _r = await asyncio.to_thread(mmg.vcf_stats, _path)
            _qual = f"mean {_r['mean_qual']:.2f}" if _r["mean_qual"] is not None else "none"
            _rows = [
                ("Records", f"`{_r['records']:,}` in `{len(_r['sequences']):,}` sequence(s), with `{_r['samples']:,}` sample(s)"),
                ("Sites by type", f"SNP `{_r['snps']:,}` · MNP `{_r['mnps']:,}` · insertion `{_r['insertions']:,}` · deletion `{_r['deletions']:,}` · other `{_r['other']:,}` · no ALT `{_r['no_alt']:,}`"),
                ("Multiallelic sites", f"`{_r['multiallelic']:,}`"),
                ("SNP alleles", f"`{_r['transitions']:,}` transitions, `{_r['transversions']:,}` transversions"
                 + (f", Ts/Tv `{_r['ts_tv']:.3f}`" if _r["ts_tv"] is not None else "")),
                ("FILTER", f"`{_r['pass']:,}` PASS, `{_r['filter_missing']:,}` missing"),
                ("QUAL", f"{_qual}, `{_r['qual_missing']:,}` missing"),
                ("Lines skipped (not valid records)", f"`{_r['skipped']:,}`"),
            ]
            _items = [mo.callout(mo.md(f"### 📊 {_path}\n" + facts(_rows)), kind="success" if _r["records"] else "warn")]
            if _r["skipped"]:
                _items.append(skipped_panel(_path, _r["skipped"], _r["records"]))
            if _r["sequences"]:
                _items.append(mo.ui.table([{"sequence": _s, "records": _n} for _s, _n in _r["sequences"]],
                                          selection=None, pagination=True, page_size=10, label="Records per sequence"))
            _items.append(code_panel(f'print({python_call("vcf_stats", _path)}["text"])'))
            _panel = mo.vstack(_items)
        except Exception as _err:
            _panel = error_panel(_err)

    elif _task == "Filter":
        _expr = (_v["expression"] or "").strip()
        _regions = [_r for _r in (_v["regions"] or "").replace(";", " ").replace(",", " ").split() if _r]
        _samples = [_s for _s in (_v["samples"] or "").replace(";", " ").replace(",", " ").split() if _s]
        mo.stop(not (_expr or _regions or _samples), warn("Type an expression, the regions or the samples to keep."))
        _out = output_name(_v["output"], _path, "filtered")
        mo.stop(clashes(_out, [_path]), warn("The output file cannot be the input file: choose another name."))
        _kw = {}
        if _expr:
            _kw["expression"] = _expr
        if _regions:
            _kw["regions"] = _regions[0] if len(_regions) == 1 else _regions
        if _samples:
            _kw["samples"] = _samples
        try:
            with mo.status.spinner(f"Filtering {_path}..."):
                _r = await asyncio.to_thread(mmg.vcf_filter, _path, _out, **_kw)
            _pct = f" ({100 * _r['kept'] / _r['read']:.1f}%)" if _r["read"] else ""
            _rows = [
                ("Records read", f"`{_r['read']:,}`"),
                ("Records kept", f"`{_r['kept']:,}`{_pct}"),
                ("Samples", f"`{_r['samples_before']:,}` → `{_r['samples_after']:,}`"),
                ("Lines skipped (not valid records)", f"`{_r['skipped']:,}`"),
                ("Saved to", f"`{os.path.abspath(_out)}` ({size_of(_out)})"),
            ]
            _panel = mo.vstack([_x for _x in (
                mo.callout(mo.md("### ✅ Filtered\n" + facts(_rows)), kind="success" if _r["read"] else "warn"),
                skipped_panel(_path, _r["skipped"], _r["read"]),
                preview_panel(_out),
                code_panel(python_call("vcf_filter", _path, _out, **_kw)),
            ) if _x is not None])
        except Exception as _err:
            _panel = error_panel(_err)

    else:
        _ref = None if _v["reference"] == "none" else _v["reference"]
        _out = output_name(_v["output"], _path, "norm")
        mo.stop(clashes(_out, [_path, _ref]), warn("The output file cannot be one of the input files: choose another name."))
        _kw = {}
        if not _v["split"]:
            _kw["split_multiallelic"] = False
        if _ref:
            _kw["reference"] = _ref
        try:
            with mo.status.spinner(f"Normalizing {_path}..."):
                _r = await asyncio.to_thread(mmg.vcf_normalize, _path, _out, **_kw)
            _rows = [
                ("Records read / written", f"`{_r['read']:,}` / `{_r['written']:,}`"),
                ("Multiallelic sites split", f"`{_r['split']:,}`"),
                ("Variants trimmed", f"`{_r['trimmed']:,}`"),
                ("Position changed", f"`{_r['position_changed']:,}`"),
            ]
            if _ref:
                _rows += [
                    ("Indels left-aligned", f"`{_r['left_aligned']:,}`"),
                    ("REF different from the genome (only trimmed)", f"`{_r['ref_mismatch']:,}`"),
                    ("Sequence missing from the genome (only trimmed)", f"`{_r['not_in_reference']:,}`"),
                ]
            _rows += [
                ("Sorted again (positions moved)", "yes" if _r["resorted"] else "no"),
                ("Lines skipped (not valid records)", f"`{_r['skipped']:,}`"),
                ("Saved to", f"`{os.path.abspath(_out)}` ({size_of(_out)})"),
            ]
            _items = [mo.callout(mo.md("### ✅ Normalized\n" + facts(_rows)), kind="success" if _r["read"] else "warn")]
            if _r["skipped"]:
                _items.append(skipped_panel(_path, _r["skipped"], _r["read"]))
            if _ref and (_r["ref_mismatch"] or _r["not_in_reference"]):
                _items.append(warn("Some records do not match the genome. Are the VCF and the FASTA on the same reference? "
                                   "Screen **10** compares the two (Compare reference genomes)."))
            _items += [preview_panel(_out), code_panel(python_call("vcf_normalize", _path, _out, **_kw))]
            _panel = mo.vstack(_items)
        except Exception as _err:
            _panel = error_panel(_err)

    _panel
    return


@app.cell
def fastq_page_ui(get_page, mo, set_page):
    mo.stop(get_page() != "fastq")
    btn_back_fastq = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_fastq = mo.ui.button(label="🔄 Refresh the file list")
    fastq_task = mo.ui.radio(options=["Quality report", "Trim reads"], value="Quality report", inline=True, label="**What to do:**")
    _ui = mo.vstack([
        btn_back_fastq,
        mo.md("<br>\n\n# 📈 12. FASTQ Quality Report & Trimming"),
        mo.callout(mo.md("Measure the quality of the reads in one pass over a file of any size: quality along the read, base composition, GC, N bases, duplication, lengths and known adapters, in a report page that also opens offline. Then trim them as `cutadapt --no-indels` does: low-quality ends, the 3' adapter and N at the ends, dropping the reads that end up too short. The two files of a pair stay in step."), kind="info"),
        mo.hstack([fastq_task, btn_refresh_fastq], justify="space-between", align="center"),
    ])
    _ui
    return btn_back_fastq, btn_refresh_fastq, fastq_task


@app.cell
def fastq_page_example(example_picker, fastq_task, get_page, mo):
    mo.stop(get_page() != "fastq")
    fastq_example = example_picker("fastq", fastq_task.value)
    fastq_example
    return (fastq_example,)


@app.cell
def fastq_page_form(btn_refresh_fastq, example_for, fastq_example, fastq_task, files_with, get_page, label_of, make_form, mo, real_files):
    mo.stop(get_page() != "fastq")
    btn_refresh_fastq  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("fastq", fastq_task.value, fastq_example.value)
    _fastqs = files_with("fastq")

    if fastq_task.value == "Quality report":
        form_fastq = make_form({
            "file": mo.ui.dropdown(options=_fastqs, value=_ex.get("file"), label="FASTQ file:", full_width=True),
            "max_reads": mo.ui.number(start=0, stop=10**12, step=10000, value=_ex.get("max_reads", 0), label="Look only at the first N reads (0 = the whole file):"),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: <file>.report.html", label="Report page:", full_width=True),
        }, "📈 Measure")
    else:
        _second = {"(single-end: no second file)": "none"}
        for _f in real_files("fastq"):
            _second[_f] = _f
        _adapters = {
            "auto: the known adapter found in the first reads": "auto",
            "none: no adapter trimming": "none",
            "Illumina Universal: AGATCGGAAGAGC": "AGATCGGAAGAGC",
            "Nextera Transposase: CTGTCTCTTATACACATCT": "CTGTCTCTTATACACATCT",
            "Illumina Small RNA 3': TGGAATTCTCGGGTGCCAAGG": "TGGAATTCTCGGGTGCCAAGG",
            "BGI/MGI: AAGTCGGAGGCCAAGCGGTC": "AAGTCGGAGGCCAAGCGGTC",
            "another one: type it below": "custom",
        }
        form_fastq = make_form({
            "file": mo.ui.dropdown(options=_fastqs, value=_ex.get("file"), label="Reads (R1):", full_width=True),
            "file2": mo.ui.dropdown(options=_second, value=label_of(_second, _ex.get("file2"), "(single-end: no second file)"), label="Mate reads (R2), for paired-end data:", full_width=True),
            "adapter": mo.ui.dropdown(options=_adapters, value=label_of(_adapters, _ex.get("adapter"), "auto: the known adapter found in the first reads"), label="3' adapter:", full_width=True),
            "custom": mo.ui.text(value=_ex.get("custom", ""), placeholder="ACGT... (only with 'another one')", label="Adapter sequence:", full_width=True),
            "quality": mo.ui.slider(start=0, stop=40, step=1, value=_ex.get("quality", 20), show_value=True, label="Trim the 3' end below this Phred quality (0 = off):"),
            "min_length": mo.ui.number(start=0, stop=100000, step=1, value=_ex.get("min_length", 20), label="Drop reads shorter than (bases):"),
            "trim_n": mo.ui.checkbox(value=_ex.get("trim_n", False), label="Remove N bases from both ends"),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: <R1 file>.trimmed.<same extension>", label="Output for R1:", full_width=True),
            "output2": mo.ui.text(value=_ex.get("output2", ""), placeholder="blank: <R2 file>.trimmed.<same extension>", label="Output for R2:", full_width=True),
        }, "✂️ Trim")

    mo.vstack([_x for _x in (_about, form_fastq) if _x is not None])
    return (form_fastq,)


@app.cell
async def fastq_page_backend(asyncio, chosen, clashes, code_panel, error_panel, facts, fastq_task, form_fastq, get_page, mmg, mo, os, output_name, python_call, size_of, warn):
    mo.stop(get_page() != "fastq")
    mo.stop(form_fastq.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_fastq.value
    _path = _v["file"]
    mo.stop(not chosen(_path), warn("Select a FASTQ file."))

    if fastq_task.value == "Quality report":
        _out = output_name(_v["output"], _path, "report", ".html")
        mo.stop(clashes(_out, [_path]), warn("The report cannot overwrite the reads: choose another name."))
        _reads = int(_v["max_reads"] or 0)
        _kw = {"max_reads": _reads} if _reads > 0 else {}
        try:
            with mo.status.spinner(f"Measuring {_path}..."):
                _r = await asyncio.to_thread(mmg.fastq_report, _path, _out, **_kw)
            _adapter = _r["likely_adapter"]
            _rows = [
                ("Reads", f"`{_r['reads']:,}`" + ("" if _r["complete"] else " (only the beginning of the file)")),
                ("Bases", f"`{_r['bases']:,}`"),
                ("Read length", f"min `{_r['min_length']:,}`, mean `{_r['mean_length']:,.1f}`, max `{_r['max_length']:,}`"),
                ("Quality", f"mean `{_r['mean_quality']:.1f}` (Phred+{_r['phred_offset']})"),
                ("GC", f"`{_r['gc_percent']:.1f}%`"),
                ("N bases", f"`{_r['n_percent']:.3f}%`"),
                ("Duplicated reads", f"`{_r['duplication_percent']:.1f}%` (of the first {_r['duplication_sample']:,})"),
                ("Adapter", f"{_adapter['name']}, `{_adapter['sequence']}`" if _adapter else "none of the known ones"),
                ("Report page", f"`{os.path.abspath(_out)}` ({size_of(_out)})"),
            ]
            _items = [mo.callout(mo.md(f"### 📈 {_path}\n" + facts(_rows)), kind="success")]
            for _w in _r["warnings"]:
                _items.append(warn(_w.replace("fastq_trim(", "**Trim reads** above, or fastq_trim(")))
            with open(_out, encoding="utf-8") as _fh:
                _items.append(mo.iframe(_fh.read(), height="1100px"))
            _items.append(code_panel(python_call("fastq_report", _path, _out, **_kw)))
            _panel = mo.vstack(_items)
        except Exception as _err:
            _panel = error_panel(_err)

    else:
        _path2 = None if _v["file2"] == "none" else _v["file2"]
        mo.stop(_path2 == _path, warn("R1 and R2 must be different files."))
        _adapter = _v["adapter"]
        if _adapter == "custom":
            _adapter = (_v["custom"] or "").strip().upper()
            mo.stop(not _adapter or set(_adapter) - set("ACGTN"), warn("Type the adapter sequence, with A, C, G, T and N only."))
        _out = output_name(_v["output"], _path, "trimmed")
        _out2 = output_name(_v["output2"], _path2, "trimmed") if _path2 else None
        mo.stop(clashes(_out, [_path, _path2]) or (_out2 is not None and (clashes(_out2, [_path, _path2]) or _out2 == _out)),
                warn("The outputs must be two new files, different from the inputs."))
        _kw = {}
        if _path2:
            _kw.update(input2=_path2, output2=_out2)
        if _adapter != "none":
            _kw["adapter"] = _adapter
        if int(_v["quality"]) > 0:
            _kw["quality_cutoff"] = int(_v["quality"])
        if _v["trim_n"]:
            _kw["trim_n"] = True
        if int(_v["min_length"] or 0) > 0:
            _kw["min_length"] = int(_v["min_length"])
        try:
            with mo.status.spinner(f"Trimming {_path}" + (f" and {_path2}" if _path2 else "") + "..."):
                _r = await asyncio.to_thread(mmg.fastq_trim, _path, _out, **_kw)
            _unit = "pairs" if _r["paired"] else "reads"
            _kept = f" ({100 * _r['reads_out'] / _r['reads_in']:.1f}%)" if _r["reads_in"] else ""
            _used = f"`{_r['adapter']}`" if _r["adapter"] else "none"
            if _r["paired"] and _r["adapter2"] and _r["adapter2"] != _r["adapter"]:
                _used += f" (R2: `{_r['adapter2']}`)"
            _rows = [
                (f"{_unit.capitalize()} in / kept", f"`{_r['reads_in']:,}` / `{_r['reads_out']:,}`{_kept}"),
                ("Bases in / kept", f"`{_r['bases_in']:,}` / `{_r['bases_out']:,}`"),
                ("Adapter used", _used + f", found in `{_r['reads_with_adapter']:,}` read(s)"
                 + (f" of R1 and `{_r['reads_with_adapter2']:,}` of R2" if _r["paired"] else "")),
                ("Bases removed", f"quality `{_r['bases_removed_quality']:,}` · adapter `{_r['bases_removed_adapter']:,}` · N `{_r['bases_removed_n']:,}`"),
                (f"{_unit.capitalize()} dropped", f"too short `{_r['discarded_too_short']:,}` · too many N `{_r['discarded_too_many_n']:,}` · low quality `{_r['discarded_low_quality']:,}`"),
                ("Saved to", f"`{os.path.abspath(_out)}`" + (f" and `{os.path.abspath(_out2)}`" if _out2 else "")),
            ]
            _panel = mo.vstack([
                mo.callout(mo.md("### ✂️ Trimmed\n" + facts(_rows)), kind="success"),
                code_panel(python_call("fastq_trim", _path, _out, **_kw)),
            ])
        except Exception as _err:
            _panel = error_panel(_err)

    _panel
    return


@app.cell
def kmercount_page_ui(example_picker, get_page, mo, set_page):
    mo.stop(get_page() != "kmercount")
    btn_back_kmercount = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_kmercount = mo.ui.button(label="🔄 Refresh the file list")
    kmercount_example = example_picker("kmercount")
    _ui = mo.vstack([
        btn_back_kmercount,
        mo.md("<br>\n\n# 🔢 13. K-mer Counting"),
        mo.callout(mo.md("Count every k-mer (k from 1 to 63) of FASTA and FASTQ files of any size: what does not fit in memory is sorted in runs on disk and merged. Canonical by default, so a k-mer and its reverse complement count as one and reads from the two strands add up. The abundance spectrum shows how many distinct k-mers were seen once, twice, and so on: in sequencing reads, the k-mers seen once are mostly errors, and the peak sits at the coverage."), kind="info"),
        mo.hstack([kmercount_example, btn_refresh_kmercount], justify="space-between", align="center"),
    ])
    _ui
    return btn_back_kmercount, btn_refresh_kmercount, kmercount_example


@app.cell
def kmercount_page_form(btn_refresh_kmercount, example_for, files_with, get_page, kmercount_example, make_form, mo):
    mo.stop(get_page() != "kmercount")
    btn_refresh_kmercount  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("kmercount", None, kmercount_example.value)
    form_kmercount = make_form({
        "files": mo.ui.multiselect(options=files_with("fasta", "fastq"), value=_ex.get("files", []), label="FASTA / FASTQ files:", full_width=True),
        "k": mo.ui.number(start=1, stop=63, step=1, value=_ex.get("k", 31), label="k (1 to 63):"),
        "canonical": mo.ui.checkbox(value=_ex.get("canonical", True), label="Canonical: count a k-mer and its reverse complement as one"),
        "min_count": mo.ui.number(start=1, stop=10**9, step=1, value=_ex.get("min_count", 1), label="Leave out of the table the k-mers seen fewer times than:"),
        "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: the summary only, no table (for example kmers.tsv.gz)", label="Table of counts:", full_width=True),
    }, "🔢 Count")
    mo.vstack([_x for _x in (_about, form_kmercount) if _x is not None])
    return (form_kmercount,)


@app.cell
async def kmercount_page_backend(asyncio, chosen, clashes, code_panel, error_panel, facts, figs, figure_panel, form_kmercount, get_page, mmg, mo, os, python_call, size_of, spectrum_chart, warn):
    mo.stop(get_page() != "kmercount")
    mo.stop(form_kmercount.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_kmercount.value
    _paths = [_p for _p in _v["files"] if chosen(_p)]
    mo.stop(not _paths, warn("Select at least one file."))
    _k = int(_v["k"])
    mo.stop(not 1 <= _k <= 63, warn("k goes from 1 to 63."))
    _out = (_v["output"] or "").strip() or None
    mo.stop(_out is not None and clashes(_out, _paths), warn("The table cannot overwrite an input file: choose another name."))
    _kw = {"k": _k}
    if not _v["canonical"]:
        _kw["canonical"] = False
    if int(_v["min_count"] or 1) > 1:
        _kw["min_count"] = int(_v["min_count"])
    try:
        with mo.status.spinner(f"Counting the {_k}-mers of {len(_paths)} file(s)..."):
            _r = await asyncio.to_thread(mmg.count_kmers, _paths, _out, **_kw)
        _once = f" ({100 * _r['unique_kmers'] / _r['distinct_kmers']:.1f}% of the distinct)" if _r["distinct_kmers"] else ""
        _rows = [
            ("k", f"`{_r['k']}`, " + ("canonical" if _r["canonical"] else "as read (one strand)")),
            ("Sequences / bases read", f"`{_r['sequences']:,}` / `{_r['bases']:,}`"),
            ("K-mers counted", f"`{_r['total_kmers']:,}`"),
            ("Distinct k-mers", f"`{_r['distinct_kmers']:,}`"),
            ("Seen only once", f"`{_r['unique_kmers']:,}`{_once}"),
            ("Highest count", f"`{_r['max_count']:,}`"),
            ("Runs on disk", f"`{_r['runs']}`" + (" (everything fit in memory)" if _r["runs"] == 0 else "")),
        ]
        if _out:
            _rows.append(("Table", f"`{_r['written']:,}` line(s) in `{os.path.abspath(_out)}` ({size_of(_out)})"))
        _items = [mo.callout(mo.md("### 🔢 Counted\n" + facts(_rows)), kind="success")]
        _model = figs.spectrum_model(_r["histogram"], k=_r["k"])
        if _model.get("genome_size") and _r["canonical"]:
            _est = [("Coverage peak", f"`{_model['peak']}×`" + (f", with a heterozygous peak at `{_model['het_peak']}×`" if _model.get("het_peak") and _model["het_peak"] != _model["peak"] else "")),
                    ("Genome size (from the spectrum)", f"`{figs.fmt_bp(_model['genome_size'])}`" + (" (one haplotype)" if _model.get("het_peak") else "")),
                    ("Error k-mers", f"`{100 * (_model['error_fraction'] or 0):.1f}%` of all k-mers counted (below the valley at {_model['valley']}×)"),
                    ("Repeats", f"`{100 * (_model['repeat_fraction'] or 0):.1f}%` of the solid k-mers are seen more than 2.5 × the peak")]
            if _model.get("heterozygosity") is not None:
                _est.append(("Heterozygosity (rough)", f"`{100 * _model['heterozygosity']:.2f}%`"))
            _items.append(mo.callout(mo.md("### 🔭 What the spectrum says\n" + facts(_est) + "\n\n" + "; ".join(_model["notes"]) +
                                           ". Read from the peaks, as one reads GenomeScope's plot; GenomeScope fits a model and "
                                           "is the reference for a published estimate."), kind="info"))
        _chart = spectrum_chart(_r["histogram"])
        if _chart is not None:
            _items.append(mo.hstack([
                figure_panel(figs.kmer_spectrum(_r["histogram"], title="Spectrum (log scale): errors on the left, repeats on the right", model=_model),
                             "kmer_spectrum_log"),
                figure_panel(figs.kmer_spectrum(_r["histogram"], title="Spectrum (linear): the coverage peaks", log_y=False, model=_model),
                             "kmer_spectrum_linear"),
            ], widths="equal", wrap=True))
            _items.append(mo.ui.table([{"times seen": _c, "distinct k-mers": _n} for _c, _n in _r["histogram"]],
                                      selection=None, pagination=True, page_size=10, label="The spectrum as a table"))
        _items.append(code_panel(f"result = {python_call('count_kmers', _paths, _out, **_kw)}", 'print(result["text"])'))
        _panel = mo.vstack(_items)
    except Exception as _err:
        _panel = error_panel(_err)

    _panel
    return


@app.cell
def convert_page_ui(example_picker, get_page, mo, set_page):
    mo.stop(get_page() != "convert")
    btn_back_convert = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_convert = mo.ui.button(label="🔄 Refresh the file list")
    convert_example = example_picker("convert")
    _ui = mo.vstack([
        btn_back_convert,
        mo.md("<br>\n\n# 📦 14. Formats & Conversion"),
        mo.callout(mo.md("Change only the container of a VCF, BED, GFF/GTF or SAM file, keeping every line and its order: plain text, gzip (written as BGZF, which can also be indexed), LZ4, BAM (from and to SAM) and Parquet, a table of typed columns that pandas, polars, DuckDB or Spark open directly, with the header kept in the file, so it converts back to the same text."), kind="info"),
        mo.hstack([convert_example, btn_refresh_convert], justify="space-between", align="center"),
    ])
    _ui
    return btn_back_convert, btn_refresh_convert, convert_example


@app.cell
def convert_page_form(btn_refresh_convert, convert_example, example_for, files_with, get_page, label_of, make_form, mo):
    mo.stop(get_page() != "convert")
    btn_refresh_convert  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("convert", None, convert_example.value)
    _targets = {
        "plain text": "text",
        "gzip, written as BGZF (.gz)": ".gz",
        "LZ4 (.lz4)": ".lz4",
        "BAM (.bam), for SAM": ".bam",
        "Parquet (.parquet)": ".parquet",
    }
    _numbers = {
        "as numbers (filter and sum without conversion)": "number",
        "as text, exactly as written (a byte-for-byte round trip)": "text",
    }
    form_convert = make_form({
        "file": mo.ui.dropdown(options=files_with("vcf", "bed", "gff", "sam", "parquet"), value=_ex.get("file"), label="File:", full_width=True),
        "target": mo.ui.dropdown(options=_targets, value=label_of(_targets, _ex.get("target"), "gzip, written as BGZF (.gz)"), label="Convert to:", full_width=True),
        "numbers": mo.ui.radio(options=_numbers, value=label_of(_numbers, _ex.get("numbers"), "as numbers (filter and sum without conversion)"),
                               label="QUAL (VCF) and score (GFF) in a Parquet output:"),
        "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: the same name with the new extension", label="Output file:", full_width=True),
    }, "📦 Convert")
    mo.vstack([_x for _x in (_about, form_convert) if _x is not None])
    return (form_convert,)


@app.cell
async def convert_page_backend(asyncio, beside, chosen, clashes, code_panel, error_panel, facts, form_convert, get_page, mmg, mo, os, preview_panel, python_call, size_of, split_name, warn):
    mo.stop(get_page() != "convert")
    mo.stop(form_convert.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_convert.value
    _path = _v["file"]
    mo.stop(not chosen(_path), warn("Select a file."))
    _target = _v["target"]
    _out = (_v["output"] or "").strip()
    if not _out:
        # 'x.vcf.gz' -> vcf; 'x.bam' -> sam; a bare 'x.parquet' says nothing.
        _stem, _ext = split_name(_path)
        _fmt = (_ext.lower().split(".") + ["", ""])[1]
        _fmt = {"bam": "sam", "parquet": ""}.get(_fmt, _fmt)
        mo.stop(not _fmt, warn("The name of this file does not say its format (`x.vcf.parquet` would): type the output name, for example `x.vcf.gz`."))
        mo.stop(_target == ".bam" and _fmt != "sam", warn("BAM holds alignments (SAM) only."))
        _out = beside(_path, _stem + {"text": f".{_fmt}", ".bam": ".bam"}.get(_target, f".{_fmt}{_target}"))
    mo.stop(clashes(_out, [_path]), warn("The file is already in that container: choose another one, or type another output name."))
    _parquet = _out.lower().endswith(".parquet")
    _text_numbers = _parquet and _v["numbers"] == "text"
    try:
        with mo.status.spinner(f"Converting {_path}..."):
            # parquet_options is a setting of the whole process: set it for this
            # conversion only, and put it back.
            _before = mmg.parquet_options()["numbers"]
            mmg.parquet_options(numbers="text" if _text_numbers else "number")
            try:
                _r = await asyncio.to_thread(mmg.convert_file, _path, _out)
            finally:
                mmg.parquet_options(numbers=_before)
        _a, _b = os.path.getsize(_path), os.path.getsize(_out)
        _ratio = (f", {_a / _b:.1f}× smaller" if _b and _a > _b * 1.05 else f", {_b / _a:.1f}× larger" if _a and _b > _a * 1.05 else "")
        _rows = [
            ("Lines", f"`{_r['lines']:,}` (`{_r['records']:,}` records)"),
            ("Size", f"{size_of(_path)} → {size_of(_out)}{_ratio}"),
            ("Saved to", f"`{os.path.abspath(_out)}`"),
        ]
        _items = [mo.callout(mo.md("### ✅ Converted\n" + facts(_rows)), kind="success")]
        if _parquet:
            _items.append(mo.md(f"Open it as a table with `pandas.read_parquet({_out!r})` or `polars.read_parquet({_out!r})`."))
        else:
            _items.append(preview_panel(_out))
        _calls = ['mmg.parquet_options(numbers="text")'] if _text_numbers else []
        _items.append(code_panel(*_calls, python_call("convert_file", _path, _out)))
        _panel = mo.vstack(_items)
    except Exception as _err:
        _panel = error_panel(_err)

    _panel
    return


@app.cell
def guide_page_ui(get_page, human_notice, mo, set_page):
    mo.stop(get_page() != "guide")
    btn_back_guide = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))

    # (test, scale, what it sees well, what it does not see, files, what MMGEN does with them)
    _tests = [
        ("G-banded karyotype", "Whole chromosomes (> 5–10 Mb)",
         "Aneuploidies (Down, Turner, Klinefelter syndromes); balanced translocations and inversions",
         "Microdeletions; single-base variants", "Images; ISCN text", "Nothing: no sequence data"),
        ("Chromosomal microarray (CGH or SNP array)", "Kilobases to megabases",
         "Microdeletions and microduplications (CNVs); a SNP array also shows runs of homozygosity (consanguinity, uniparental disomy)",
         "Balanced rearrangements; single-base variants", "CEL / IDAT, then genotype tables or VCF",
         "Runs of homozygosity and genetic sex from a VCF (screen 20)"),
        ("MLPA", "Exons (tens of bases)", "Deletions and duplications of target genes (SMN1, DMD, 22q11.2)",
         "Point variants outside the probes", "Fragment-analysis tables", "Planned: copy number from the read depth of a BAM"),
        ("Gene panel or exome (short reads)", "Single base", "SNVs and small indels in the targeted genes or coding regions",
         "Non-coding regions; long repeat expansions; complex rearrangements; many CNVs", "FASTQ → BAM → VCF",
         "Quality and trimming (12), sorting and indexing (9), VCF tools (11), gene panels and ClinVar findings (18)"),
        ("Whole genome, short reads (WGS)", "Single base, plus structural", "Almost everything: SNVs, indels, CNVs, regulatory regions",
         "Extreme repeat expansions; regions with near-identical copies (pseudogenes, e.g. CYP2D6, SMN1); long-range phasing",
         "FASTQ, BAM / CRAM, VCF", "The file tools (9–14), for files of any size, and the human-genome screens (16–21)"),
        ("Trio whole genome (child + parents)", "As WGS, for three people", "De novo variants; inheritance checks; parentage",
         "As WGS", "Three VCFs (or one with three samples)", "Trio analysis, de novo candidates, parentage and relatedness (screen 16)"),
        ("Whole genome, long reads (PacBio HiFi, Nanopore)", "Single base over tens of kilobases",
         "Structural variants, repeat expansions, pseudogene regions, native methylation, phasing",
         "Very low-frequency variants; it costs much more", "FASTQ / BAM (methylation tags), VCF",
         "The file tools read its BAM and VCF; no long-read-specific analysis yet"),
        ("Ancestry and polygenic-score tests (genotyping array or genome)", "Hundreds of thousands of SNPs",
         "Continental ancestry proportions; polygenic scores", "Fine-scale or ethnic origins; scores carry over poorly to other ancestries",
         "Genotype tables, VCF", "Ancestry proportions, polygenic scores and reference panels (screen 17)"),
        ("Optical genome mapping", "Large structures (from ~500 bases up)", "Balanced translocations, inversions, large insertions",
         "Single-base variants, small indels", "BNX / CMAP", "Nothing"),
        ("Tumor profiling / liquid biopsy", "Single base at 500–1000× depth",
         "Low-frequency somatic variants, gene fusions, tumor mutational burden (TMB), microsatellite instability (MSI)",
         "Very low tumor fraction; inherited variants are filtered out", "BAM, somatic VCF", "Tumor mutational burden from the somatic VCF and the panel's BED (screen 28)"),
        ("Sanger sequencing, qPCR", "One amplicon or target", "Confirming a known variant; testing relatives for it (cascade testing)",
         "Anything outside the amplicon", "AB1 traces; Ct values", "A variant looked up in relatives, and the allele counts at a position of a BAM (screen 21)"),
        ("Pharmacogenomic panel (+ CYP2D6 copy number)", "Selected variants and copy number",
         "Star alleles and metabolizer phenotype (CYP2C19, CYP2C9, DPYD, TPMT...); HLA risk alleles",
         "Rare alleles not on the panel; CYP2D6 duplications and hybrids without a copy-number assay", "VCF or genotype table",
         "Star alleles, phenotypes and guideline texts by MMGEN's engine or PharmCAT (screen 19); from a BAM, CYP2D6 copy number and HLA-A, -B, -C"),
        ("Shotgun metagenomics (mNGS)", "All the DNA of a sample", "Bacteria, viruses, fungi and parasites at once; unknown pathogens",
         "Organisms at very low abundance under the host DNA; RNA viruses unless RNA is sequenced", "FASTQ",
         "Species and relative abundance, read classification, host reads apart (screen 24); adventitious agents in a biological material, with blanks, spike-ins and breadth of coverage (screen 31)"),
        ("16S rRNA / ITS amplicons", "One barcode gene", "The bacterial or fungal composition, usually down to genus",
         "Resistance genes, strains, viruses", "FASTQ → ASV / OTU tables", "Amplicon variants, taxonomy (RDP classifier), alpha and beta diversity, PERMANOVA (screen 27)"),
        ("Pathogen whole genome", "Single base", "Outbreak tracking (SNP distances), resistance genes and mutations, viral variants",
         "Mixed infections at a low level", "FASTQ, FASTA assemblies, VCF", "Resistance genes of the AMRFinderPlus catalog (screen 25); distance trees and viral mutations against the reference (screen 26)"),
        ("Syndromic multiplex PCR panel", "15–40 targets", "Fast detection of the listed pathogens and resistance genes (mecA, KPC, NDM)",
         "Anything not on the panel", "Ct values / detected-or-not", "In-silico PCR of a panel's primers against genomes or reads (screen 26)"),
        ("RNA-seq", "Gene expression", "Which genes are switched on or off (for example after a drug)",
         "DNA variants in genes that are not expressed", "FASTQ → BAM → count table",
         "Quantification straight from the FASTQ (as kallisto), or counts per gene from a BAM; differential expression, as DESeq2, and gene-set tests (screen 22)"),
        ("CRISPR screen", "One gene knocked out per cell", "Which genes cause resistance or sensitivity to a drug or a virus",
         "Genes whose loss changes nothing the screen measures", "FASTQ → guide counts", "Guide counts and gene ranking, the same as MAGeCK (screen 23)"),
    ]
    _areas = [
        ("Human genetics (inherited)", "The diploid human genome, as inherited", "Exome, panels, array, karyotype",
         "Pathogenic variants (ACMG criteria), inherited risk"),
        ("Somatic oncology", "The tumor genome (acquired mutations)", "Tumor panels, liquid biopsy",
         "Biomarkers (TMB, MSI), mutations with a targeted therapy"),
        ("Microbiology and metagenomics", "Haploid and viral genomes, plasmids", "Shotgun mNGS, 16S rRNA, pathogen WGS",
         "Taxonomy (relative abundance), phylogenetic trees, resistance genes (AMR)"),
        ("Pharmacogenomics", "Polymorphic drug-response genes (CYPs, HLA)", "SNP/NGS panels + copy number (CYP2D6)",
         "Star alleles (e.g. *1/*4), metabolizer phenotype, CPIC guidelines"),
    ]

    def _table(head, rows):
        _md = "| " + " | ".join(head) + " |\n|" + " :--- |" * len(head) + "\n"
        return _md + "\n".join("| " + " | ".join(r) + " |" for r in rows)

    _ui = mo.vstack([
        btn_back_guide,
        mo.md("<br>\n\n# 📚 15. Which Test Sees What"),
        mo.callout(mo.md("No genetic test sees everything. **A normal result rules out only what that test can see**: a normal exome does not rule out a balanced translocation, and a normal karyotype does not rule out a single-base variant. This guide lists the common tests, from the cheapest to the most expensive, with their blind spots, the files they produce, and what MMGEN does with those files. The screens marked *planned* are being built."), kind="info"),
        human_notice(),
        mo.md("### The tests and their blind spots\n\n" + _table(["Test", "Scale", "Sees well", "Does not see", "Files", "In MMGEN"], _tests)),
        mo.md("### Four areas, four kinds of result\n\n" + _table(["Area", "Genome studied", "Main tests", "Kind of result"], _areas)
              + "\n\nIn pharmacogenomics the names change: instead of a coordinate (chr10:94,781,859 G>A) a result is a **star allele** "
              "(CYP2C19 \\*2) and a pair of them, the diplotype (\\*1/\\*2), read against the CPIC guidelines. In metagenomics the result "
              "is a table of abundances (OTUs or ASVs) and a list of resistance genes."),
    ])
    _ui
    return (btn_back_guide,)


@app.cell
def kin_page_ui(get_page, human_notice, mo, set_page):
    mo.stop(get_page() != "kin")
    btn_back_kin = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_kin = mo.ui.button(label="🔄 Refresh the file list")
    kin_task = mo.ui.radio(options=["Relatedness of several people", "Trio and paternity test"], value="Relatedness of several people",
                           inline=True, label="**What to do:**")
    kin_folder = mo.ui.text(value=".", label="Data folder (its subfolders are searched too):", full_width=True)
    _ui = mo.vstack([
        btn_back_kin,
        mo.md("<br>\n\n# 👪 16. Relatedness, Trios & Paternity"),
        human_notice(),
        mo.callout(mo.md(
            "**Relatedness**: for every pair of people, the KING kinship coefficient — about 0.5 for the same person or identical "
            "twins, 0.25 for parent and child or full siblings, 0.125 for second-degree relatives, 0.0625 for first cousins, 0 for "
            "unrelated people — and the share of sites where the two are opposite homozygotes (IBS0), near zero only between a "
            "parent and a child. **Trio**: a child's genotypes checked against the father's and the mother's (Mendelian "
            "consistency), the variants found in the child only (de novo candidates), and a parentage test: how often the "
            "alleged father and the child are opposite homozygotes, which a true father explains only by genotyping errors, and "
            "the combined paternity index from population allele frequencies.\n\n"
            "**This is not a paternity test with legal value**: that requires an accredited laboratory, verified identities, a "
            "chain of custody and the consent of everyone tested."), kind="info"),
        kin_folder,
        mo.hstack([kin_task, btn_refresh_kin], justify="space-between", align="center"),
    ])
    _ui
    return btn_back_kin, btn_refresh_kin, kin_folder, kin_task


@app.cell
def kin_page_example(example_picker, get_page, kin_task, mo):
    mo.stop(get_page() != "kin")
    kin_example = example_picker("kin", kin_task.value)
    kin_example
    return (kin_example,)


@app.cell
def kin_page_form(NO_FILE, btn_refresh_kin, example_for, files_under, get_page, kin_example, kin_folder, kin_task, label_of, make_form, mo):
    mo.stop(get_page() != "kin")
    btn_refresh_kin  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("kin", kin_task.value, kin_example.value)
    _vcfs = files_under(kin_folder.value, "vcf")
    _beds_label = "Use the callable-regions BED next to each VCF (`x.callable.bed`, or the Genome in a Bottle `*_benchmark*.bed`)"
    if kin_task.value.startswith("Relatedness"):
        form_kin = make_form({
            "vcfs": mo.ui.multiselect(options=_vcfs, value=_ex.get("vcfs", []), label="People's VCFs:", full_width=True),
            "samples": mo.ui.text(value=_ex.get("samples", ""), placeholder="blank: every sample", label="Samples, for a multi-sample VCF:", full_width=True),
            "use_beds": mo.ui.checkbox(value=_ex.get("use_beds", True), label=_beds_label),
        }, "👪 Compare")
    else:
        _people = list(_vcfs)
        _mothers = {"(none: a duo, child and alleged father only)": "none", **{_f: _f for _f in _vcfs if _f != NO_FILE}}
        _afs = {"(none: no paternity index)": "none", **{_f: _f for _f in _vcfs if _f != NO_FILE}}
        _fields = {"AF (all populations)": "AF", "AFR_AF (African)": "AFR_AF", "AMR_AF (admixed American)": "AMR_AF",
                   "EAS_AF (East Asian)": "EAS_AF", "EUR_AF (European)": "EUR_AF", "SAS_AF (South Asian)": "SAS_AF"}
        form_kin = make_form({
            "child": mo.ui.dropdown(options=_people, value=_ex.get("child"), label="Child:", full_width=True),
            "father": mo.ui.dropdown(options=_people, value=_ex.get("father"), label="Alleged father:", full_width=True),
            "mother": mo.ui.dropdown(options=_mothers, value=label_of(_mothers, _ex.get("mother", "none"), "(none: a duo, child and alleged father only)"), label="Mother:", full_width=True),
            "use_beds": mo.ui.checkbox(value=_ex.get("use_beds", True), label=_beds_label),
            "af": mo.ui.dropdown(options=_afs, value=label_of(_afs, _ex.get("af", "none"), "(none: no paternity index)"),
                                 label="Population allele frequencies (e.g. the 1000 Genomes sites VCF), for the paternity index:", full_width=True),
            "af_field": mo.ui.dropdown(options=_fields, value=label_of(_fields, _ex.get("af_field", "AF"), "AF (all populations)"),
                                       label="Frequencies of which population (the closest to the family's ancestry):", full_width=True),
            "de_novo": mo.ui.text(value=_ex.get("de_novo", ""), placeholder="blank: only count them", label="Write the de novo candidates to (VCF):", full_width=True),
        }, "👪 Analyze the trio")
    mo.vstack([_x for _x in (_about, form_kin) if _x is not None])
    return (form_kin,)


@app.cell
async def kin_page_backend(asyncio, bed_next_to, chosen, code_panel, error_panel, facts, figs, figure_panel, form_kin, get_page, kin_task, mmg, mo, python_call, warn, words):
    mo.stop(get_page() != "kin")
    mo.stop(form_kin.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_kin.value
    _names = {"duplicate_or_mz_twin": "the same person or identical twins", "parent_offspring": "parent and child",
              "full_siblings": "full siblings", "second_degree": "second degree (grandparent, uncle or aunt, half-sibling)",
              "third_degree": "third degree (first cousins)", "unrelated": "unrelated", "insufficient_data": "too few sites to tell"}
    if kin_task.value.startswith("Relatedness"):
        _vcfs = [_p for _p in _v["vcfs"] if chosen(_p)]
        mo.stop(not _vcfs, warn("Select the VCFs (one per person), or one VCF with several samples."))
        _kw = {}
        _samples = words(_v["samples"])
        if _samples:
            _kw["samples"] = _samples
        if _v["use_beds"]:
            _beds = [bed_next_to(_p) for _p in _vcfs] if len(_vcfs) > 1 else [bed_next_to(_vcfs[0])]
            if any(_beds):
                _kw["callable_beds"] = _beds if len(_vcfs) > 1 else _beds[0]
        mo.stop(len(_vcfs) == 1 and len(_samples) == 1, warn("Two people at least: select more VCFs, or more samples."))
        try:
            with mo.status.spinner("Comparing every pair..."):
                _r = await asyncio.to_thread(mmg.relatedness, _vcfs if len(_vcfs) > 1 else _vcfs[0], **_kw)
            _lab = [_p["label"] for _p in _r["persons"]]
            _who = lambda x: _lab[x] if isinstance(x, int) else x
            _pairs = sorted(_r["pairs"], key=lambda _p: -_p["kinship"])
            _rows = [("People", f"`{len(_lab)}`: {', '.join(_lab)}"), ("SNV sites used", f"`{_r['snv_sites']:,}` on `{_r['autosomes']}` autosome(s)")]
            _items = [mo.callout(mo.md("### 👪 Relatedness\n" + facts(_rows)), kind="success")]
            _items.append(mo.ui.table([{
                "person A": _who(_p["a"]), "person B": _who(_p["b"]), "kinship": round(_p["kinship"], 4),
                "relationship": _names.get(_p["relationship"], _p["relationship"]), "sites": _p["sites"],
                "IBS0 / both heterozygous": round(_p["ibs0_ratio"], 4) if _p["ibs0_ratio"] is not None else "",
            } for _p in _pairs], selection=None, pagination=True, page_size=15, label="Every pair, most related first"))
            _items.append(figure_panel(
                figs.king_plot([(_p["kinship"], _p["ibs0_ratio"], f"{_who(_p['a'])} – {_who(_p['b'])}",
                                 _names.get(_p["relationship"], _p["relationship"])) for _p in _pairs]),
                "kinship_ibs0",
                "Parent and child share an allele at nearly every site (IBS0 near 0) with kinship near 0.25; full siblings have "
                "the same kinship but some IBS0. A pair between the bands deserves a second look (a contaminated or swapped "
                "sample, or ancestry far from the others: KING-robust kinship turns negative between populations)."))
            if _r["autosomes"] < 10:
                _items.append(warn(f"Only {_r['autosomes']} autosome(s): the relationship classes are tentative. A whole genome tells them apart much better."))
            _items.append(code_panel(f"r = {python_call('relatedness', _vcfs if len(_vcfs) > 1 else _vcfs[0], **_kw)}", 'print(r["text"])'))
            _panel = mo.vstack(_items)
        except Exception as _err:
            _panel = error_panel(_err)
    else:
        mo.stop(not chosen(_v["child"]) or not chosen(_v["father"]), warn("Select the child and the alleged father."))
        _mother = None if _v["mother"] == "none" else _v["mother"]
        mo.stop(len({_v["child"], _v["father"], _mother} - {None}) < (3 if _mother else 2), warn("The child, the father and the mother must be different files."))
        _kw = {}
        if _mother:
            _kw["mother"] = _mother
        if _v["use_beds"]:
            for _role, _path in (("child", _v["child"]), ("father", _v["father"]), ("mother", _mother)):
                if _path and bed_next_to(_path):
                    _kw[f"{_role}_bed"] = bed_next_to(_path)
        if _v["af"] != "none":
            _kw["allele_frequencies"] = _v["af"]
            _kw["af_field"] = _v["af_field"]
        if (_v["de_novo"] or "").strip():
            _kw["de_novo_output"] = _v["de_novo"].strip()
        try:
            with mo.status.spinner("Checking the trio (a minute for whole genomes)..."):
                _r = await asyncio.to_thread(mmg.trio_analysis, _v["child"], _v["father"], **_kw)
            _verdict = {"not_excluded": ("✅ The alleged father is not excluded", "success"),
                        "excluded": ("❌ The alleged father is excluded", "danger"),
                        "inconclusive": ("❔ Inconclusive", "warn")}.get(_r["verdict"], (_r["verdict"], "warn"))
            _p = _r["probability_of_paternity"]
            _prob = ("" if _p is None else "> 99.9999%" if _p > 0.999999 else "< 0.0001%" if _p < 0.000001 else f"{100 * _p:.4f}%")
            _rows = [
                ("Exclusions (father and child opposite homozygotes)", f"`{_r['father_exclusions']:,}` of `{_r['father_informative']:,}` informative sites "
                 f"(rate `{_r['exclusion_rate']:.2e}`; genotyping errors alone: `{_r['error_level']:.1e}`"
                 + (f"; an unrelated man: about `{_r['unrelated_expected_rate']:.2e}`" if _r.get("unrelated_expected_rate") else "") + ")"),
                ("Mendelian errors", f"`{_r['errors']:,}` of `{_r['sites_checked']:,}` sites ({100 * _r['error_rate']:.3f}%)"),
                ("De novo candidates", f"`{_r['de_novo']:,}`" + (f", written to `{_kw['de_novo_output']}`" if "de_novo_output" in _kw else "")),
            ]
            if _r["log10_cpi"] is not None:
                _rows.append(("Combined paternity index", f"10^{_r['log10_cpi']:,.0f} over `{_r['cpi_sites']:,}` sites; probability of paternity {_prob} (prior {_r['prior']})"))
            _items = [mo.callout(mo.md(f"### {_verdict[0]}\n{_r['reason']}\n\n" + facts(_rows)), kind=_verdict[1])]
            _items.append(mo.callout(mo.md(
                "**Why the numbers are so large**: the paternity index treats every site as independent, which they are not along a "
                "chromosome, so it is far larger than the 15-25 markers of a forensic test would give. Read it as a teaching tool. "
                "De novo candidates are mostly genotyping errors; real de novo variants are confirmed in the laboratory."), kind="info"))
            _items.append(code_panel(f"r = {python_call('trio_analysis', _v['child'], _v['father'], **_kw)}", 'print(r["text"])'))
            _panel = mo.vstack(_items)
        except Exception as _err:
            _panel = error_panel(_err)
    _panel
    return


@app.cell
def anc_page_ui(get_page, human_notice, mo, set_page):
    mo.stop(get_page() != "anc")
    btn_back_anc = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_anc = mo.ui.button(label="🔄 Refresh the file list")
    anc_task = mo.ui.radio(options=["Ancestry", "Polygenic score", "Build a reference panel"], value="Ancestry", inline=True, label="**What to do:**")
    anc_folder = mo.ui.text(value=".", label="Data folder (its subfolders are searched too):", full_width=True)
    _ui = mo.vstack([
        btn_back_anc,
        mo.md("<br>\n\n# 🌍 17. Ancestry & Polygenic Scores"),
        human_notice(),
        mo.callout(mo.md(
            "**Ancestry**: the proportions of a person's genome that best match each reference group of a **panel** of allele "
            "frequencies (for example the five continental groups of the 1000 Genomes Project), as ADMIXTURE's supervised mode "
            "computes them. The groups are broad samples of present-day people, not ethnic or national identities, and a "
            "person's answer depends on which groups the panel has. **Polygenic score**: the weighted sum of many small-effect "
            "variants from a published score (the PGS Catalog), and where it falls among the panel's groups. A score is a "
            "statistical summary for research and teaching: it is not a diagnosis, and it carries over poorly to ancestries "
            "unlike those of the study that made it."), kind="info"),
        anc_folder,
        mo.hstack([anc_task, btn_refresh_anc], justify="space-between", align="center"),
    ])
    _ui
    return anc_folder, anc_task, btn_back_anc, btn_refresh_anc


@app.cell
def anc_page_example(anc_task, example_picker, get_page, mo):
    mo.stop(get_page() != "anc")
    anc_example = example_picker("anc", anc_task.value)
    anc_example
    return (anc_example,)


@app.cell
def anc_page_form(NO_FILE, anc_example, anc_folder, anc_task, btn_refresh_anc, example_for, files_under, get_page, label_of, make_form, mo):
    mo.stop(get_page() != "anc")
    btn_refresh_anc  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("anc", anc_task.value, anc_example.value)
    _folder = anc_folder.value
    _vcfs = files_under(_folder, "vcf")
    _beds_label = "Use the callable-regions BED next to each VCF (`x.callable.bed`, or the Genome in a Bottle `*_benchmark*.bed`)"
    _tsvgz =[_f for _f in files_under(_folder, "vcf", "table") if _f != NO_FILE]
    if anc_task.value == "Ancestry":
        _pan = [_f for _f in _tsvgz if "panel" in _f.lower() and _f.endswith((".tsv", ".tsv.gz", ".txt"))] or [NO_FILE]
        form_anc = make_form({
            "panel": mo.ui.dropdown(options=_pan, value=_ex.get("panel"), label="Reference panel (a frequency panel TSV):", full_width=True),
            "vcfs": mo.ui.multiselect(options=_vcfs, value=_ex.get("vcfs", []), label="People's VCFs:", full_width=True),
            "samples": mo.ui.text(value=_ex.get("samples", ""), placeholder="blank: every sample", label="Samples, for a multi-sample VCF:", full_width=True),
            "use_beds": mo.ui.checkbox(value=_ex.get("use_beds", True), label=_beds_label),
        }, "🌍 Estimate")
    elif anc_task.value == "Polygenic score":
        _scores = [_f for _f in _tsvgz if "pgs" in _f.lower() and not _f.lower().endswith(".vcf.gz")] or [NO_FILE]
        _pan = {"(none: the score only, no comparison)": "none", **{_f: _f for _f in _tsvgz if "panel" in _f.lower()}}
        _imp = {"pooled frequency of the panel": "pooled", "AFR": "AFR", "AMR": "AMR", "EAS": "EAS", "EUR": "EUR", "SAS": "SAS"}
        form_anc = make_form({
            "scoring": mo.ui.dropdown(options=_scores, value=_ex.get("scoring"), label="Scoring file (PGS Catalog, harmonized GRCh38):", full_width=True),
            "vcfs": mo.ui.multiselect(options=_vcfs, value=_ex.get("vcfs", []), label="People's VCFs:", full_width=True),
            "use_beds": mo.ui.checkbox(value=_ex.get("use_beds", True), label=_beds_label),
            "panel": mo.ui.dropdown(options=_pan, value=label_of(_pan, _ex.get("panel", "none"), "(none: the score only, no comparison)"),
                                    label="Frequency panel, to place the score among its groups:", full_width=True),
            "impute": mo.ui.dropdown(options=_imp, value=label_of(_imp, _ex.get("impute", "pooled"), "pooled frequency of the panel"),
                                     label="Variants a file cannot call are filled with the frequency of:", full_width=True),
            "keep_ambiguous": mo.ui.checkbox(value=_ex.get("keep_ambiguous", False), label="Keep A/T and C/G variants (their strand cannot be checked)"),
        }, "📈 Score")
    else:
        _kinds = {"a sites VCF with each group's frequency in INFO (AFR_AF, EUR_AF...)": "freq",
                  "a multi-sample VCF of genotypes, with a table of each sample's group": "geno"}
        _tables = {"(not needed for a sites VCF)": "none", **{_f: _f for _f in files_under(_folder, "table") if _f != NO_FILE}}
        form_anc = make_form({
            "source": mo.ui.dropdown(options=_vcfs, value=_ex.get("source"), label="Source VCF:", full_width=True),
            "kind": mo.ui.dropdown(options=_kinds, value=label_of(_kinds, _ex.get("kind", "freq"), "a sites VCF with each group's frequency in INFO (AFR_AF, EUR_AF...)"), label="It is:", full_width=True),
            "groups": mo.ui.text(value=_ex.get("groups", "AFR, AMR, EAS, EUR, SAS"), label="Groups (for a sites VCF: read from <group>_AF):", full_width=True),
            "table": mo.ui.dropdown(options=_tables, value="(not needed for a sites VCF)", label="Table of samples and groups (for genotypes):", full_width=True),
            "columns": mo.ui.text(value=_ex.get("columns", "SampleID, Superpopulation"), label="Its sample and group columns:", full_width=True),
            "min_maf": mo.ui.number(start=0.0, stop=0.5, step=0.01, value=_ex.get("min_maf", 0.05), label="Keep sites whose minor-allele frequency reaches, in some group:"),
            "thin_bp": mo.ui.number(start=0, stop=10**7, step=1000, value=_ex.get("thin_bp", 0), label="At most one site per (bases; 0 = all):"),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: <source>.panel.tsv.gz next to it", label="Panel file:", full_width=True),
        }, "🧱 Build the panel")
    mo.vstack([_x for _x in (_about, form_anc) if _x is not None])
    return (form_anc,)


@app.cell
async def anc_page_backend(anc_task, asyncio, bed_next_to, beside, chosen, code_panel, error_panel, facts, figs, figure_panel, form_anc, get_page, mmg, mo, python_call, split_name, warn, words):
    mo.stop(get_page() != "anc")
    mo.stop(form_anc.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_anc.value
    try:
        if anc_task.value in ("Ancestry", "Polygenic score"):
            _vcfs = [_p for _p in _v["vcfs"] if chosen(_p)]
            mo.stop(not _vcfs, warn("Select the people's VCFs."))
            _target = _vcfs if len(_vcfs) > 1 else _vcfs[0]
            _kw = {}
            if anc_task.value == "Ancestry" and words(_v["samples"]):
                _kw["samples"] = words(_v["samples"])
            if _v["use_beds"]:
                _beds = [bed_next_to(_p) for _p in _vcfs]
                if any(_beds) and len(_vcfs) > 1:
                    _kw["callable_beds"] = _beds
                elif _beds[0] and "samples" not in _kw:
                    _kw["callable_beds"] = _beds[0]
        if anc_task.value == "Ancestry":
            mo.stop(not chosen(_v["panel"]), warn("Select a reference panel (build one with **Build a reference panel**)."))
            with mo.status.spinner("Estimating the proportions (with a bootstrap for the errors)..."):
                _r = await asyncio.to_thread(mmg.estimate_ancestry, _v["panel"], _target, **_kw)
            _groups = _r[0]["groups"] if _r else []
            _rows = []
            for _x in _r:
                _row = {"person": _x["sample"] or _x["vcf"]}
                for _g, _p, _s in zip(_x["groups"], _x["proportions"], _x["standard_errors"]):
                    _row[_g] = f"{100 * _p:.1f}%" + (f" ± {100 * _s:.1f}" if _s is not None else "")
                _row["best single group"] = _x["best_group"]
                _row["sites"] = _x["sites_used"]
                _rows.append(_row)
            _warns = sorted({_w for _x in _r for _w in _x["warnings"]})
            _items = [mo.callout(mo.md("### 🌍 Ancestry proportions\n" + facts([("People", f"`{len(_r)}`"), ("Panel groups", ", ".join(_groups))])), kind="success"),
                      figure_panel(figs.structure_plot([(_x["sample"] or _x["vcf"], _x["proportions"]) for _x in _r], _groups),
                                   "ancestry_structure",
                                   "One bar per person, sorted by the main group. Proportions are relative to the panel's groups: "
                                   "ancestry from outside them (e.g. Middle Eastern against the five 1000 Genomes groups) shows as a "
                                   "mixture of the nearest ones."),
                      mo.ui.table(_rows, selection=None, pagination=True, page_size=15, label="Proportions (± bootstrap standard error)")]
            for _w in _warns[:5]:
                _items.append(warn(_w))
            # The PCA picture: the people projected on clouds simulated from the panel's frequencies.
            if hasattr(mmg, "ancestry_pca"):
                try:
                    with mo.status.spinner("Projecting the people on the panel's groups (PCA)..."):
                        _pca = await asyncio.to_thread(mmg.ancestry_pca, _v["panel"], _target, **_kw)
                    _ve = _pca.get("variance_explained") or [0, 0]
                    _ref_pts = [(_x["pcs"][0], _x["pcs"][1], _x["group"], _x["group"]) for _x in _pca["reference"] if _x.get("pcs")]
                    _people = [(_x["pcs"][0], _x["pcs"][1], _x["name"]) for _x in _pca["people"] if _x.get("pcs")]
                    _items.append(figure_panel(
                        figs.group_scatter(_ref_pts, f"PC1 ({100 * _ve[0]:.1f}%)", f"PC2 ({100 * _ve[1]:.1f}%)" if len(_ve) > 1 else "PC2",
                                           f"The people on the panel's groups (PCA of {_pca['sites_used']:,} sites)", extra_points=_people),
                        "ancestry_pca",
                        "Each cloud is a group simulated from the panel's allele frequencies (real people of a group spread "
                        "2 to 4 times more: the clouds ignore linkage and the structure inside each group). A person between two "
                        "clouds has ancestry from both; one far from every cloud has ancestry the panel does not hold."))
                    _items.append(mo.ui.table([{"person": _x["name"], "closest group": _x.get("closest_group") or "",
                                                "distance (cloud spreads)": round(_x["distance"], 1) if _x.get("distance") is not None else "",
                                                "sites observed": _x.get("sites_observed")} for _x in _pca["people"]],
                                              selection=None, label="Each person's place in the PCA"))
                except Exception as _perr:
                    _items.append(mo.md(f"*The PCA could not be drawn: {_perr}*"))
            _items.append(code_panel(f"r = {python_call('estimate_ancestry', _v['panel'], _target, **_kw)}", 'for x in r:\n    print(x["text"])',
                                     f"pca = {python_call('ancestry_pca', _v['panel'], _target, **_kw)}"))
        elif anc_task.value == "Polygenic score":
            mo.stop(not chosen(_v["scoring"]), warn("Select a scoring file (download one from pgscatalog.org: the harmonized GRCh38 file)."))
            if _v["panel"] != "none":
                _kw["panel"] = _v["panel"]
            if _v["impute"] != "pooled":
                _kw["impute"] = _v["impute"]
            if _v["keep_ambiguous"]:
                _kw["keep_ambiguous"] = True
            with mo.status.spinner("Scoring..."):
                _r = await asyncio.to_thread(mmg.polygenic_score, _v["scoring"], _target, **_kw)
            _first = _r[0] if _r else {}
            _rows = []
            for _x in _r:
                _row = {"person": _x["sample"] or _x["vcf"], "score": round(_x["score"], 4), "coverage": f"{100 * _x['coverage']:.1f}%"}
                for _g in _x["groups"]:
                    _row[f"percentile in {_g['group']}"] = f"{_g['percentile']:.0f}" if _g.get("percentile") is not None else ""
                _rows.append(_row)
            _facts = [("Score", f"{_first.get('pgs_id') or ''} {_first.get('trait_reported') or ''}".strip() or "—"),
                      ("Variants", f"`{_first.get('variants', 0):,}` in the file, `{_first.get('usable', 0):,}` usable")]
            _items = [mo.callout(mo.md("### 📈 Polygenic score\n" + facts(_facts)), kind="success"),
                      mo.ui.table(_rows, selection=None, pagination=True, page_size=15, label="The score of each person")]
            _items.append(mo.callout(mo.md(
                "**Reading a percentile**: the score of a group's average person is at the 50th; the comparison assumes independent "
                "variants in Hardy-Weinberg equilibrium and is only meaningful against the group of the person's own ancestry. "
                "**A polygenic score is not a diagnosis** and does not say whether a person will develop a condition."), kind="warn"))
            for _w in sorted({_w for _x in _r for _w in _x["warnings"]})[:5]:
                _items.append(warn(_w))
            _items.append(code_panel(f"r = {python_call('polygenic_score', _v['scoring'], _target, **_kw)}", 'for x in r:\n    print(x["text"])'))
        else:
            mo.stop(not chosen(_v["source"]), warn("Select the source VCF."))
            _groups = words(_v["groups"])
            _out = (_v["output"] or "").strip() or beside(_v["source"], split_name(_v["source"])[0] + ".panel.tsv.gz")
            _kw = {}
            if float(_v["min_maf"]) != 0.05:
                _kw["min_maf"] = float(_v["min_maf"])
            if int(_v["thin_bp"] or 0):
                _kw["thin_bp"] = int(_v["thin_bp"])
            if _v["kind"] == "freq":
                mo.stop(not _groups, warn("Type the groups."))
                with mo.status.spinner("Building the panel..."):
                    _r = await asyncio.to_thread(mmg.panel_from_frequencies, _v["source"], _groups, [_g + "_AF" for _g in _groups], _out, **_kw)
                _call = python_call("panel_from_frequencies", _v["source"], _groups, [_g + "_AF" for _g in _groups], _out, **_kw)
            else:
                _cols = words(_v["columns"])
                mo.stop(_v["table"] == "none" or len(_cols) != 2, warn("Select the table of samples and type its sample and group columns."))
                with mo.status.spinner("Counting the alleles of every group..."):
                    _r = await asyncio.to_thread(mmg.build_panel, _v["source"], _v["table"], _out, _cols[0], _cols[1], **_kw)
                _call = python_call("build_panel", _v["source"], _v["table"], _out, _cols[0], _cols[1], **_kw)
            _items = [mo.callout(mo.md("### 🧱 Panel written\n" + facts([
                ("Groups", ", ".join(_r["groups"])),
                ("Sites written", f"`{_r['written']:,}` of `{_r['records']:,}` records"),
                ("Saved to", f"`{_out}`"),
            ])), kind="success"), code_panel(f"r = {_call}", 'print(r["text"])')]
        _panel = mo.vstack(_items)
    except Exception as _err:
        _panel = error_panel(_err)
    _panel
    return


@app.cell
def chrom_page_ui(get_page, human_notice, mo, set_page):
    mo.stop(get_page() != "chrom")
    btn_back_chrom = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_chrom = mo.ui.button(label="🔄 Refresh the file list")
    chrom_task = mo.ui.radio(options=["Genetic sex", "Runs of homozygosity"], value="Genetic sex", inline=True, label="**What to do:**")
    chrom_folder = mo.ui.text(value=".", label="Data folder (its subfolders are searched too):", full_width=True)
    _ui = mo.vstack([
        btn_back_chrom,
        mo.md("<br>\n\n# 🧫 20. Sex Chromosomes & Homozygosity"),
        human_notice(),
        mo.callout(mo.md(
            "Two things a SNP array or a karyotype shows that a VCF can show too. **Genetic sex**: with two X chromosomes, the "
            "calls on X outside the pseudoautosomal regions are heterozygous about as often as the autosomes' calls; with one X, "
            "almost never; and only a Y gives calls on Y. Patterns in between can suggest X0 (Turner) or XXY (Klinefelter), which "
            "only a karyotype confirms. **Runs of homozygosity**: long stretches (almost) without heterozygous calls, inherited "
            "twice from the same ancestor — many and long when the parents are related; one long run on a single chromosome can "
            "point to uniparental disomy. Runs are found as `plink --homozyg` finds them (windows of 50 calls that tolerate one "
            "heterozygous call, so a genotyping error does not split a run), or strictly (any heterozygous call ends a run). "
            "Genetic sex is not gender identity."), kind="info"),
        chrom_folder,
        mo.hstack([chrom_task, btn_refresh_chrom], justify="space-between", align="center"),
    ])
    _ui
    return btn_back_chrom, btn_refresh_chrom, chrom_folder, chrom_task


@app.cell
def chrom_page_example(chrom_task, example_picker, get_page, mo):
    mo.stop(get_page() != "chrom")
    chrom_example = example_picker("chrom", chrom_task.value)
    chrom_example
    return (chrom_example,)


@app.cell
def chrom_page_form(NO_FILE, btn_refresh_chrom, chrom_example, chrom_folder, chrom_task, example_for, files_under, get_page, label_of, make_form, mo):
    mo.stop(get_page() != "chrom")
    btn_refresh_chrom  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("chrom", chrom_task.value, chrom_example.value)
    _vcfs = files_under(chrom_folder.value, "vcf")
    _beds = {"(none)": "none", **{_b: _b for _b in files_under(chrom_folder.value, "bed_any") if _b != NO_FILE}}
    _fields = {
        "vcf": mo.ui.dropdown(options=_vcfs, value=_ex.get("vcf"), label="The person's VCF:", full_width=True),
        "sample": mo.ui.text(value=_ex.get("sample", ""), placeholder="blank: the first sample", label="Sample, in a multi-sample VCF:"),
        "callable": mo.ui.dropdown(options=_beds, value=label_of(_beds, _ex.get("callable", "none"), "(none)"), label="Callable regions (BED):", full_width=True),
    }
    if chrom_task.value == "Runs of homozygosity":
        _fields["min_length"] = mo.ui.number(start=10_000, stop=10**8, step=100_000, value=_ex.get("min_length", 1_000_000), label="Minimum length of a run (bases):")
        _fields["min_sites"] = mo.ui.number(start=5, stop=10**6, step=5, value=_ex.get("min_sites", 50), label="Minimum SNV calls in a run:")
        _methods = {"Sliding windows, as plink --homozyg (tolerates genotyping errors)": "window", "Strict: any heterozygous call ends a run": "strict"}
        _fields["method"] = mo.ui.dropdown(options=_methods, value=label_of(_methods, _ex.get("method", "window"), "Sliding windows, as plink --homozyg (tolerates genotyping errors)"),
                                           label="Method:", full_width=True)
        _fields["window_het"] = mo.ui.number(start=0, stop=10, step=1, value=_ex.get("window_het", 1),
                                             label="Windows: heterozygous calls a window of 50 may hold (plink's default: 1):")
        _fields["output"] = mo.ui.text(value=_ex.get("output", ""), placeholder="blank: only shown here", label="Save the runs to (BED):", full_width=True)
        form_chrom = make_form(_fields, "🧫 Find the runs")
    else:
        form_chrom = make_form(_fields, "🧫 Infer")
    mo.vstack([_x for _x in (_about, form_chrom) if _x is not None])
    return (form_chrom,)


@app.cell
async def chrom_page_backend(asyncio, chosen, chrom_task, code_panel, error_panel, facts, figs, figure_panel, form_chrom, get_page, mmg, mo, python_call, warn):
    mo.stop(get_page() != "chrom")
    mo.stop(form_chrom.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_chrom.value
    mo.stop(not chosen(_v["vcf"]), warn("Select the person's VCF."))
    _kw = {}
    if (_v["sample"] or "").strip():
        _kw["sample"] = _v["sample"].strip()
    if _v["callable"] != "none":
        _kw["callable_bed"] = _v["callable"]
    try:
        if chrom_task.value == "Genetic sex":
            with mo.status.spinner("Reading the sex chromosomes..."):
                _r = await asyncio.to_thread(mmg.genetic_sex, _v["vcf"], **_kw)
            _head = {"XX": "♀ XX: two X chromosomes", "XY": "♂ XY: one X and a Y", "undetermined": "❔ Undetermined"}.get(_r["verdict"], _r["verdict"])
            _rows = [
                ("X outside the pseudoautosomal regions", f"`{_r['x_het']:,}` heterozygous, `{_r['x_hom_alt']:,}` homozygous, `{_r['x_haploid']:,}` haploid calls"
                 + (f"; heterozygosity `{_r['x_het_ratio']:.2f}` of the autosomes'" if _r.get("x_het_ratio") is not None else "")),
                ("Y outside the pseudoautosomal regions", f"`{_r['y_calls']:,}` calls (`{_r['y_non_ref']:,}` non-reference)"),
                ("Assembly", f"{_r['build']} ({_r['build_source']})"),
            ]
            _flags = ", ".join(_r["flags"]) if _r["flags"] else "none"
            _items = [mo.callout(mo.md(f"### {_head}\n{_r['reason']}\n\n" + facts(_rows) + f"\n* **Flags:** {_flags}"),
                                 kind="success" if _r["verdict"] in ("XX", "XY") else "warn")]
            _items.append(code_panel(f"r = {python_call('genetic_sex', _v['vcf'], **_kw)}", 'print(r["text"])'))
        else:
            if int(_v["min_length"]) != 1_000_000:
                _kw["min_length"] = int(_v["min_length"])
            if int(_v["min_sites"]) != 50:
                _kw["min_sites"] = int(_v["min_sites"])
            if _v["method"] == "strict":
                _kw["window_snps"] = 0
            elif int(_v["window_het"]) != 1:
                _kw["window_het"] = int(_v["window_het"])
            if (_v["output"] or "").strip():
                _kw["output"] = _v["output"].strip()
            with mo.status.spinner("Looking for runs of homozygosity..."):
                _r = await asyncio.to_thread(mmg.runs_of_homozygosity, _v["vcf"], **_kw)
            _rows = [
                ("Runs", f"`{_r['count']}`, `{_r['total_length'] / 1e6:.2f}` Mb in all, the longest `{_r['longest'] / 1e6:.2f}` Mb"),
                ("F_ROH", f"`{_r['f_roh']:.4f}` of the autosomes" + (f" (`{_r['f_roh_covered']:.4f}` of the span the file covers)" if _r.get("f_roh_covered") is not None else "")),
                ("SNV calls", f"`{_r['snv_calls']:,}` (`{_r['het_calls']:,}` heterozygous)"),
                ("Method", f"sliding windows of {_r['window_snps']} calls with at most {_r['window_het']} heterozygous, as plink --homozyg"
                 if _r["method"] == "window" else "strict: any heterozygous call ends a run"),
            ]
            _items = [mo.callout(mo.md("### 🧫 Runs of homozygosity\n" + facts(_rows)), kind="success")]
            _runs_k = [(_x["chrom"], _x["start0"], _x["end"]) if isinstance(_x, dict) else (_x[0], _x[1], _x[2]) for _x in _r["runs"]]
            _items.append(figure_panel(
                figs.karyogram(_runs_k), "runs_of_homozygosity",
                "A few long runs on several chromosomes: the parents share a recent ancestor; many short runs everywhere: an old, "
                "small population; one chromosome covered end to end: possibly uniparental disomy. A chromosome the file does not "
                "cover shows no runs — no data is not the same as no homozygosity."))
            if _r["runs"]:
                _runs = []
                for _x in _r["runs"]:
                    _c, _s, _e, _n = (_x["chrom"], _x["start0"], _x["end"], _x["sites"]) if isinstance(_x, dict) else _x
                    _runs.append({"chromosome": _c, "start (1-based)": _s + 1, "end": _e, "length (Mb)": round((_e - _s) / 1e6, 2), "SNV calls": _n})
                _items.append(mo.ui.table(_runs, selection=None, pagination=True, page_size=15, label="The runs"))
            _items.append(mo.md("*For orientation: a child of first cousins has F_ROH around 1/16 (0.0625) on average; most people with "
                                "unrelated parents have little, in short runs. Consanguinity is common in many families and communities; "
                                "a genetic counselor can say what it means for a family.*"))
            _items.append(code_panel(f"r = {python_call('runs_of_homozygosity', _v['vcf'], **_kw)}", 'print(r["text"])'))
        _panel = mo.vstack(_items)
    except Exception as _err:
        _panel = error_panel(_err)
    _panel
    return


@app.cell
def findings_page_ui(get_page, human_notice, mo, set_page):
    mo.stop(get_page() != "findings")
    btn_back_findings = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_findings = mo.ui.button(label="🔄 Refresh the file list")
    findings_task = mo.ui.radio(options=["ClinVar findings", "Gene panel from an annotation"], value="ClinVar findings",
                                inline=True, label="**What to do:**")
    findings_folder = mo.ui.text(value=".", label="Data folder (its subfolders are searched too):", full_width=True)
    _ui = mo.vstack([
        btn_back_findings,
        mo.md("<br>\n\n# 🧾 18. Disease Findings & Gene Panels"),
        human_notice(),
        mo.callout(mo.md("Which of a person's variants does **ClinVar** (NCBI's public archive of variant classifications) list, and how is each one classified: pathogenic, likely pathogenic, uncertain, benign. Each finding comes with its review stars (how well supported the classification is), the conditions and the zygosity. Restrict it to the genes of a **panel**, as the gene panels of clinical laboratories do, or get the coordinates of a panel's genes from an annotation (GFF). The person's VCF and ClinVar must be on the same assembly; download ClinVar from ftp.ncbi.nlm.nih.gov/pub/clinvar/."), kind="info"),
        findings_folder,
        mo.hstack([findings_task, btn_refresh_findings], justify="space-between", align="center"),
    ])
    _ui
    return btn_back_findings, btn_refresh_findings, findings_folder, findings_task


@app.cell
def findings_page_example(example_picker, findings_task, get_page, mo):
    mo.stop(get_page() != "findings")
    findings_example = example_picker("findings", findings_task.value)
    findings_example
    return (findings_example,)


@app.cell
def findings_page_form(GENE_PANELS, NO_FILE, btn_refresh_findings, example_for, files_under, findings_example, findings_folder, findings_task, get_page, label_of, make_form, mo):
    mo.stop(get_page() != "findings")
    btn_refresh_findings  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("findings", findings_task.value, findings_example.value)
    _folder = findings_folder.value
    _panels = {"(no panel: the whole file)": "none", **{_n: _n for _n in GENE_PANELS}}
    if findings_task.value == "ClinVar findings":
        _vcfs = files_under(_folder, "vcf")
        _beds = {"(none: compare every call in the file)": "none"}
        for _b in files_under(_folder, "bed_any"):
            if _b != NO_FILE:
                _beds[_b] = _b
        _classes = {
            "pathogenic and likely pathogenic": "plp",
            "uncertain significance (VUS)": "vus",
            "conflicting classifications": "conflicting",
            "benign and likely benign": "blb",
            "risk factors, protective and other": "risk",
        }
        _clinvar = _ex.get("clinvar") or next((_f for _f in _vcfs if "clinvar" in _f.lower()), None)
        form_findings = make_form({
            "vcf": mo.ui.dropdown(options=_vcfs, value=_ex.get("vcf"), label="The person's VCF:", full_width=True),
            "sample": mo.ui.text(value=_ex.get("sample", ""), placeholder="blank: the first sample", label="Sample, in a multi-sample VCF:"),
            "callable": mo.ui.dropdown(options=_beds, value=label_of(_beds, _ex.get("callable", "none"), "(none: compare every call in the file)"),
                                       label="Callable regions (BED): calls outside them are not compared:", full_width=True),
            "clinvar": mo.ui.dropdown(options=_vcfs, value=_clinvar, label="ClinVar VCF (the same assembly as the person's):", full_width=True),
            "classes": mo.ui.multiselect(options=_classes, value=[label_of(_classes, _c, None) for _c in _ex.get("classes", ["plp"])],
                                         label="Classifications to show:", full_width=True),
            "min_stars": mo.ui.slider(start=0, stop=4, step=1, value=_ex.get("min_stars", 1), show_value=True,
                                      label="Minimum review stars (0 = no assertion, 4 = practice guideline):"),
            "panel": mo.ui.dropdown(options=_panels, value=label_of(_panels, _ex.get("panel", "none"), "(no panel: the whole file)"),
                                    label="Only the genes of an example panel:", full_width=True),
            "genes": mo.ui.text(value=_ex.get("genes", ""), placeholder="BRCA1, BRCA2", label="... and/or these genes:", full_width=True),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: shown here only", label="Also save the table (TSV):", full_width=True),
        }, "🧾 Find")
    else:
        _namings = {"as the annotation says": "keep", "chr1, chr2... (UCSC)": "chr", "1, 2... (Ensembl, NCBI)": "no_chr"}
        form_findings = make_form({
            "gff": mo.ui.dropdown(options=files_under(_folder, "gff"), value=_ex.get("gff"), label="Annotation (GFF3 / GTF, e.g. the RefSeq GRCh38 GFF):", full_width=True),
            "panel": mo.ui.dropdown(options=_panels, value=label_of(_panels, _ex.get("panel", "none"), "(no panel: the whole file)"),
                                    label="An example panel:", full_width=True),
            "genes": mo.ui.text(value=_ex.get("genes", ""), placeholder="BRCA1, BRCA2", label="... and/or these genes:", full_width=True),
            "padding": mo.ui.number(start=0, stop=10**7, step=100, value=_ex.get("padding", 0), label="Bases to add on each side:"),
            "naming": mo.ui.dropdown(options=_namings, value=label_of(_namings, _ex.get("naming", "chr"), "chr1, chr2... (UCSC)"), label="Chromosome names:", full_width=True),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: panel.bed next to the annotation", label="Save the regions to (BED):", full_width=True),
        }, "🧬 Get the regions")
    mo.vstack([_x for _x in (_about, form_findings) if _x is not None])
    return (form_findings,)


@app.cell
async def findings_page_backend(GENE_PANELS, asyncio, beside, chosen, code_panel, error_panel, facts, findings_task, form_findings, get_page, mmg, mo, os, python_call, warn, words):
    mo.stop(get_page() != "findings")
    mo.stop(form_findings.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_findings.value
    _genes = list(dict.fromkeys(GENE_PANELS.get(_v["panel"], []) + [_g.upper() for _g in words(_v["genes"])]))

    if findings_task.value == "ClinVar findings":
        mo.stop(not chosen(_v["vcf"]), warn("Select the person's VCF."))
        mo.stop(not chosen(_v["clinvar"]), warn("Select the ClinVar VCF (clinvar.vcf.gz, from ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/)."))
        mo.stop(_v["vcf"] == _v["clinvar"], warn("The person's VCF and ClinVar must be two different files."))
        mo.stop(not _v["classes"], warn("Select at least one classification."))
        _kw = {"classes": list(_v["classes"]), "min_stars": int(_v["min_stars"])}
        if (_v["sample"] or "").strip():
            _kw["sample"] = _v["sample"].strip()
        if _v["callable"] != "none":
            _kw["callable_bed"] = _v["callable"]
        if _genes:
            _kw["genes"] = _genes
        if (_v["output"] or "").strip():
            _kw["output"] = _v["output"].strip()
        try:
            with mo.status.spinner("Comparing the person's variants with ClinVar..."):
                _r = await asyncio.to_thread(mmg.clinvar_findings, _v["vcf"], _v["clinvar"], **_kw)
            _cats = ", ".join(f"{_n:,} {_c}" for _c, _n in _r["by_category"] if _n) or "none"
            _rows = [
                ("Sample", f"`{_r['sample']}`"),
                ("Alleles checked", f"`{_r['person_alleles']:,}`" + (f" (inside the callable regions; `{_r['person_outside_callable']:,}` outside were left out)" if _r["callable"] else "")),
                ("Found in ClinVar", f"`{_r['matched']:,}`: {_cats}"),
                ("Shown here", f"`{len(_r['findings']):,}` (filters: {_r['filters']})"),
            ]
            _plp = _r["plp_heterozygous"] + _r["plp_homozygous"] + _r["plp_hemizygous"] + _r["plp_unknown_zygosity"]
            if "plp" in _kw["classes"]:
                _rows.append(("Pathogenic / likely pathogenic", f"`{_r['plp_heterozygous']}` in one copy, `{_r['plp_homozygous']}` in two copies, `{_r['plp_hemizygous']}` hemizygous"))
            _items = [mo.callout(mo.md("### 🧾 ClinVar findings\n" + facts(_rows)), kind="success")]
            if _r["findings"]:
                _items.append(mo.ui.table([{
                    "gene": ", ".join(_f["genes"]),
                    "variant": f"{_f['chrom']}:{_f['pos']} {_f['ref']}>{_f['alt']}",
                    "zygosity": _f["zygosity"],
                    "classification": _f["significance"].replace("_", " "),
                    "stars": "★" * _f["stars"] + "☆" * (4 - _f["stars"]),
                    "conditions": "; ".join(_f["conditions"][:3]) + (" ..." if len(_f["conditions"]) > 3 else ""),
                    "rsID": _f["rsid"] or "",
                    "ClinVar": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{_f['variation_id']}/",
                } for _f in _r["findings"]], selection=None, pagination=True, page_size=10, label="The findings"))
            if _plp:
                _items.append(mo.callout(mo.md(
                    "**Reading a finding.** ClinVar gathers what laboratories submitted; the stars say how well a classification is "
                    "supported. One copy of a variant for a *recessive* condition is a carrier state, which most people have for a few "
                    "conditions; one copy for a *dominant* condition is a different matter. Which applies depends on the condition, "
                    "which ClinVar's file does not say. **Only an accredited laboratory and a genetic counselor can say what a finding "
                    "means for a person.**"), kind="warn"))
            _args = {_k: _x for _k, _x in _kw.items() if not (_k == "classes" and _x == ["plp"]) and not (_k == "min_stars" and _x == 1)}
            _items.append(code_panel(f"r = {python_call('clinvar_findings', _v['vcf'], _v['clinvar'], **_args)}", 'print(r["text"])'))
            _panel = mo.vstack(_items)
        except Exception as _err:
            _panel = error_panel(_err)
    else:
        mo.stop(not chosen(_v["gff"]), warn("Select an annotation (GFF3 or GTF)."))
        mo.stop(not _genes, warn("Choose a panel or type the genes."))
        _out = (_v["output"] or "").strip() or beside(_v["gff"], "panel.bed")
        _kw = {"output": _out}
        if int(_v["padding"] or 0):
            _kw["padding"] = int(_v["padding"])
        if _v["naming"] != "keep":
            _kw["naming"] = _v["naming"]
        try:
            with mo.status.spinner(f"Reading {_v['gff']}..."):
                _r = await asyncio.to_thread(mmg.gene_regions, _v["gff"], _genes, **_kw)
            _rows = [
                ("Genes found", f"`{len(_r['regions'])}` of `{len(_genes)}`"),
                ("Not found", ", ".join(f"`{_g}`" for _g in _r["not_found"]) or "none"),
                ("Saved to", f"`{os.path.abspath(_out)}`"),
            ]
            _items = [mo.callout(mo.md("### 🧬 The panel's regions\n" + facts(_rows)), kind="success" if _r["regions"] else "warn")]
            if _r["synonyms"]:
                _items.append(mo.md("Other names to try: " + ", ".join(f"`{_s}`" for _s in _r["synonyms"][:20])))
            _items.append(mo.ui.table([{
                "gene": _g["gene"], "chrom": _g["chrom"], "start (1-based)": _g["start"] + 1, "end": _g["end"],
                "length": _g["end"] - _g["start"], "strand": _g["strand"],
            } for _g in _r["regions"]], selection=None, pagination=True, page_size=15, label="Regions"))
            _items.append(code_panel(f"r = {python_call('gene_regions', _v['gff'], _genes, **_kw)}",
                                     "# then, for example: mmg.vcf_filter(vcf, 'panel_variants.vcf', regions=[...])"))
            _panel = mo.vstack(_items)
        except Exception as _err:
            _panel = error_panel(_err)
    _panel
    return


@app.cell
def pgx_page_ui(get_page, human_notice, mo, set_page):
    mo.stop(get_page() != "pgx")
    import shutil as _shutil
    _java = _shutil.which("java")
    btn_back_pgx = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_pgx = mo.ui.button(label="🔄 Refresh the file list")
    btn_get_pharmcat = mo.ui.button(label="⬇️ Download PharmCAT 3.4.0 (32 MB, from GitHub)", value=0, on_click=lambda v: v + 1)
    btn_get_pgx_tables = mo.ui.button(label="⬇️ Get the mmgen-pgx tables (CPIC v1.61.0)", value=0, on_click=lambda v: v + 1)
    btn_get_hla = mo.ui.button(label="⬇️ Download the HLA reference sequences (68 MB, IPD-IMGT/HLA)", value=0, on_click=lambda v: v + 1)
    # MMGEN's own engine is the default (decision of the author, 2026-11-16, after it
    # met the acceptance criterion against PharmCAT); PharmCAT is optional, for comparison.
    pgx_engine = mo.ui.radio(options=["MMGEN's engine (mmgen-pgx)", "PharmCAT", "A group of people and a protocol's medicines (mmgen-pgx)"],
                             value="MMGEN's engine (mmgen-pgx)", inline=True, label="**Engine:**")
    pgx_folder = mo.ui.text(value=".", label="Data folder (its subfolders are searched too):", full_width=True)
    _ui = mo.vstack([
        btn_back_pgx,
        mo.md("<br>\n\n# 💊 19. Pharmacogenomics"),
        human_notice(),
        mo.callout(mo.md(
            "How a person's genes may change the response to medicines: **star alleles** (e.g. CYP2C19 \\*1/\\*2), the **metabolizer "
            "phenotype** and the matching guideline texts. MMGEN's own engine is the default; PharmCAT gives the same kind of "
            "result, for comparison:\n\n"
            "- **MMGEN's engine (mmgen-pgx)**: no Java needed. It reads 19 genes from the VCF with the CPIC tables (v1.61.0, public "
            "domain) and, given the person's reads (a BAM), also calls **CYP2D6 with its copy number** and the **HLA-A, -B and -C** "
            "types. When the positions it can read fit more than one diplotype, it calls the one PharmCAT would (the highest "
            "Named Allele Matcher score; for NAT2, the most common in CPIC's population frequencies), lists the others and "
            "warns when one of them would change the phenotype or a recommendation.\n"
            "- **A group of people and a protocol's medicines**: MMGEN's engine on several people at once (a clinical study, a "
            "class), with the medicines of a protocol, e.g. those given around a surgery. For each person and medicine: the "
            "recommendation, or why there is none (a gene not called, or no CPIC guideline for that medicine).\n"
            "- **PharmCAT** (optional, for comparison), the open-source tool of PharmGKB and CPIC (pharmcat.org, MPL-2.0; Java 17 "
            "or newer): MMGEN prepares its input from the VCF. It does not call CYP2D6 or HLA from a VCF: give them as an outside-call file — the one MMGEN's "
            "engine writes from a BAM, or a laboratory's."), kind="info"),
        mo.md(f"Java (for PharmCAT): {'found at `' + _java + '`' if _java else '**not found**: install Java 17 or newer (adoptium.net; on Colab: `!apt-get install -y openjdk-21-jre-headless`)'}"),
        pgx_folder,
        mo.hstack([pgx_engine, btn_refresh_pgx], justify="space-between", align="center"),
        mo.hstack([btn_get_pgx_tables, btn_get_hla, btn_get_pharmcat], justify="start", wrap=True),
    ])
    _ui
    return btn_back_pgx, btn_get_hla, btn_get_pgx_tables, btn_get_pharmcat, btn_refresh_pgx, pgx_engine, pgx_folder


@app.cell
def pgx_page_download(btn_get_hla, btn_get_pgx_tables, btn_get_pharmcat, error_panel, get_page, mo, os):
    mo.stop(get_page() != "pgx")
    # What the two engines need, each in its folder under mmgen_tools/:
    # PharmCAT with the positions file of its release; the mmgen-pgx tables
    # (their zip, from this folder or from the release where it is published);
    # and the IPD-IMGT/HLA sequences, used unmodified and never redistributed
    # with MMGEN (their licence, CC BY-ND), fetched from the official mirror of
    # IPD-IMGT/HLA on GitHub (the EBI server takes minutes for each file).
    PHARMCAT_DIR = "mmgen_tools/pharmcat"
    PGX_TABLES_VERSION = "cpic-v1.61.0"
    PGX_TABLES_DIR = f"mmgen_tools/pgx/{PGX_TABLES_VERSION}"
    PGX_TABLES_ZIP = f"mmgen_pgx_tables_{PGX_TABLES_VERSION}.zip"
    PGX_TABLES_URL = "https://github.com/fbcouto/multimerge-genomics-research/releases/download/v2.3.0/" + PGX_TABLES_ZIP
    PGX_TABLES_SHA256 = "e329b8ec6e8567dc2570cd3823a96401e54135b3a5545305fb05c7265caa1f85"
    PGX_HLA_FILES = ("A_gen.fasta", "B_gen.fasta", "C_gen.fasta")
    # The IPD-IMGT/HLA release the tables pin in bam_genes.json (the engine uses no other):
    # its branch of the official mirror, and the SHA-256 of each file. Keep equal to the tables.
    PGX_HLA_VERSION = "3.66.0"
    PGX_HLA_URL = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/3660/fasta/"
    PGX_HLA_SHA256 = {
        "A_gen.fasta": "f265427f646d20da7f26f63f30620f1787c15417dfd731a39e8818026455d931",
        "B_gen.fasta": "3748309e4aa560f8a0c31a8dfd6f24ddce2f3c9a52bd66db193de5b2f1bcd8f8",
        "C_gen.fasta": "edec46685befd8d0f47871a6e4f916f4195a8b3d5930ec55ae3fa207e117893a",
    }

    def pgx_tables_dir():
        """The folder of the mmgen-pgx tables (MMGEN_PGX_TABLES, or the one the button fills), or None."""
        for _d in (os.environ.get("MMGEN_PGX_TABLES"), PGX_TABLES_DIR):
            if _d and os.path.exists(os.path.join(_d, "manifest.json")):
                return os.path.abspath(_d)
        return None

    def pgx_hla_dir():
        """The HLA folder (MMGEN_HLA_DIR, or mmgen_tools/hla) when the three files are in it or in
        its <release>/ subfolder, or None. The engine checks their SHA-256 when it reads them."""
        _d = os.environ.get("MMGEN_HLA_DIR") or "mmgen_tools/hla"
        for _sub in (os.path.join(_d, PGX_HLA_VERSION), _d):
            if all(os.path.exists(os.path.join(_sub, _f)) for _f in PGX_HLA_FILES):
                return os.path.abspath(_d)
        return None

    import hashlib as _hashlib
    import urllib.request as _url
    _out = []
    if btn_get_pharmcat.value:
        _base = "https://github.com/PharmGKB/PharmCAT/releases/download/v3.4.0/"
        try:
            os.makedirs(PHARMCAT_DIR, exist_ok=True)
            with mo.status.spinner("Downloading PharmCAT 3.4.0..."):
                for _f in ("pharmcat-3.4.0-all.jar", "pharmcat_positions_3.4.0.vcf"):
                    if not os.path.exists(os.path.join(PHARMCAT_DIR, _f)):
                        _url.urlretrieve(_base + _f, os.path.join(PHARMCAT_DIR, _f) + ".part")
                        os.replace(os.path.join(PHARMCAT_DIR, _f) + ".part", os.path.join(PHARMCAT_DIR, _f))
            _out.append(mo.callout(mo.md(f"✅ PharmCAT 3.4.0 is in `{os.path.abspath(PHARMCAT_DIR)}`."), kind="success"))
        except Exception as _err:
            _out.append(error_panel(f"the download failed ({_err}). Download pharmcat-3.4.0-all.jar and pharmcat_positions_3.4.0.vcf "
                                    f"from https://github.com/PharmGKB/PharmCAT/releases into `{PHARMCAT_DIR}/` by hand."))
    if btn_get_pgx_tables.value:
        import zipfile as _zipfile
        try:
            _zip = next((_c for _c in (PGX_TABLES_ZIP, os.path.join("mmgen_tools", "pgx", PGX_TABLES_ZIP)) if os.path.exists(_c)), None)
            if _zip is None and PGX_TABLES_URL:
                os.makedirs(os.path.join("mmgen_tools", "pgx"), exist_ok=True)
                _zip = os.path.join("mmgen_tools", "pgx", PGX_TABLES_ZIP)
                with mo.status.spinner("Downloading the mmgen-pgx tables..."):
                    _url.urlretrieve(PGX_TABLES_URL, _zip + ".part")
                with open(_zip + ".part", "rb") as _fh:
                    _got = _hashlib.sha256(_fh.read()).hexdigest()
                if _got != PGX_TABLES_SHA256:
                    os.remove(_zip + ".part")
                    raise RuntimeError(f"the downloaded `{PGX_TABLES_ZIP}` is not the published one (SHA-256 {_got[:12]}..., "
                                       f"expected {PGX_TABLES_SHA256[:12]}...): try again, or download it by hand from {PGX_TABLES_URL}")
                os.replace(_zip + ".part", _zip)
            if _zip is None:
                raise RuntimeError(f"`{PGX_TABLES_ZIP}` is not in this folder: put it here, or set MMGEN_PGX_TABLES to the "
                                   "folder of the tables")
            with _zipfile.ZipFile(_zip) as _z:
                if "manifest.json" not in _z.namelist():
                    raise RuntimeError(f"`{_zip}` has no manifest.json: it is not an mmgen-pgx table bundle")
                _z.extractall(PGX_TABLES_DIR)
            # The engine checks the size and the SHA-256 of every table against the manifest before each run.
            _out.append(mo.callout(mo.md(f"✅ The mmgen-pgx tables ({PGX_TABLES_VERSION}) are in `{os.path.abspath(PGX_TABLES_DIR)}`."), kind="success"))
        except Exception as _err:
            _out.append(error_panel(_err))
    if btn_get_hla.value:
        _dir = os.path.join(os.environ.get("MMGEN_HLA_DIR") or "mmgen_tools/hla", PGX_HLA_VERSION)
        try:
            os.makedirs(_dir, exist_ok=True)
            with mo.status.spinner(f"Downloading A_gen.fasta, B_gen.fasta and C_gen.fasta of IPD-IMGT/HLA {PGX_HLA_VERSION} (68 MB)..."):
                for _f in PGX_HLA_FILES:
                    _p = os.path.join(_dir, _f)
                    if not os.path.exists(_p):
                        _url.urlretrieve(PGX_HLA_URL + _f, _p + ".part")
                        with open(_p + ".part", "rb") as _fh:
                            _got = _hashlib.sha256(_fh.read()).hexdigest()
                        if _got != PGX_HLA_SHA256[_f]:
                            os.remove(_p + ".part")
                            raise RuntimeError(f"the downloaded `{_f}` is not the one of IPD-IMGT/HLA {PGX_HLA_VERSION} "
                                               f"(SHA-256 {_got[:12]}..., expected {PGX_HLA_SHA256[_f][:12]}...)")
                        os.replace(_p + ".part", _p)
            _out.append(mo.callout(mo.md(
                f"✅ The HLA reference sequences of IPD-IMGT/HLA {PGX_HLA_VERSION} are in `{os.path.abspath(_dir)}`. They are the "
                "IPD-IMGT/HLA Database's files (Barker et al., Nucleic Acids Research 2023), licensed CC BY-ND: used as they are, "
                "never changed or shared by MMGEN."),
                kind="success"))
        except Exception as _err:
            _out.append(error_panel(f"the download failed ({_err}). Download {', '.join(PGX_HLA_FILES)} from {PGX_HLA_URL} "
                                    f"(IPD-IMGT/HLA {PGX_HLA_VERSION}) into `{_dir}/` by hand."))
    pgx_downloaded = btn_get_pharmcat.value + btn_get_pgx_tables.value + btn_get_hla.value
    mo.vstack(_out) if _out else None
    return PHARMCAT_DIR, pgx_downloaded, pgx_hla_dir, pgx_tables_dir


@app.cell
def pgx_page_example(example_picker, get_page, mo, pgx_engine):
    mo.stop(get_page() != "pgx")
    pgx_example = example_picker("pgx", pgx_engine.value)
    pgx_example
    return (pgx_example,)


@app.cell
def pgx_page_form(NO_FILE, PHARMCAT_DIR, btn_refresh_pgx, example_for, files_under, get_page, label_of, make_form, mo, os, pgx_downloaded, pgx_engine, pgx_example, pgx_folder):
    mo.stop(get_page() != "pgx")
    btn_refresh_pgx, pgx_downloaded  # re-list the folder after a refresh or a download
    _mine = pgx_engine.value.startswith("MMGEN")
    _ex, _about = example_for("pgx", pgx_engine.value, pgx_example.value)
    _folder = pgx_folder.value
    _vcfs = files_under(_folder, "vcf")
    if pgx_engine.value.startswith("A group"):
        _cbeds = {"next to each VCF (x.callable.bed, or Genome in a Bottle's *_benchmark*.bed), when there is one": "auto",
                  "(none: positions without a record are not read)": "none"}
        for _b in files_under(_folder, "bed_any"):
            if _b != NO_FILE:
                _cbeds[f"{_b} (for every VCF)"] = _b
        _gpanels = {"(none: the plain mean of the CPIC population groups)": "none"}
        for _t in files_under(_folder, "table"):
            if _t != NO_FILE and "panel" in os.path.basename(_t).lower():
                _gpanels[_t] = _t
        form_pgx = make_form({
            "vcfs": mo.ui.multiselect(options=[_x for _x in _vcfs if _x != NO_FILE], value=_ex.get("vcfs", []),
                                      label="People's VCFs (GRCh38; one per person, or multi-sample):", full_width=True),
            "callable": mo.ui.dropdown(options=_cbeds, value=label_of(_cbeds, _ex.get("callable", "auto"), next(iter(_cbeds))),
                                       label="Callable regions (BED):", full_width=True),
            "samples": mo.ui.text(value=_ex.get("samples", ""), placeholder="blank: every sample of every VCF",
                                  label="Only these samples (names):", full_width=True),
            "drugs": mo.ui.text_area(value=_ex.get("drugs", ""), placeholder="codeine, ondansetron, celecoxib, warfarin, clopidogrel, omeprazole, sevoflurane...",
                                     label="Medicines of the protocol (commas or lines; blank: every medicine with a recommendation):", full_width=True, rows=3),
            "ancestry_panel": mo.ui.dropdown(options=_gpanels, value=label_of(_gpanels, _ex.get("ancestry_panel", "none"), "(none: the plain mean of the CPIC population groups)"),
                                             label="Estimate each person's ancestry with a reference panel, for NAT2 (optional, screen 17):", full_width=True),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: pgx_cohort/ next to the first VCF",
                                 label="Output folder:", full_width=True),
        }, "💊 Call the pharmacogenes of everyone")
    else:
        _none_bed = ("(none: positions without a record are not read)" if _mine
                     else "(none: a position without a record counts as reference)")
        _beds = {_none_bed: "none"}
        for _b in files_under(_folder, "bed_any"):
            if _b != NO_FILE:
                _beds[_b] = _b
        _fields = {
            "vcf": mo.ui.dropdown(options=_vcfs, value=_ex.get("vcf"), label="The person's VCF (GRCh38):", full_width=True),
            "sample": mo.ui.text(value=_ex.get("sample", ""), placeholder="blank: the first sample", label="Sample, in a multi-sample VCF:"),
            "callable": mo.ui.dropdown(options=_beds, value=label_of(_beds, _ex.get("callable", "none"), _none_bed),
                                       label="Callable regions (BED): inside them, a position without a record is the reference:", full_width=True),
        }
        if _mine:
            _bams = {"(none: only the genes read from the VCF)": "none"}
            for _b in files_under(_folder, "bam"):
                if _b != NO_FILE:
                    _bams[_b] = _b
            _refs = {"(none)": "none"}
            for _f in files_under(_folder, "fasta"):
                if _f != NO_FILE:
                    _refs[_f] = _f
            _fields["bam"] = mo.ui.dropdown(options=_bams, value=label_of(_bams, _ex.get("bam", "none"), "(none: only the genes read from the VCF)"),
                                            label="The person's reads (BAM with its .bai, GRCh38) for CYP2D6 and HLA:", full_width=True)
            _fields["reference"] = mo.ui.dropdown(options=_refs, value="(none)",
                                                  label="GRCh38 genome (FASTA), to left-align the indels of the VCF (optional):", full_width=True)
            _panels = {"(none: the plain mean of the CPIC population groups)": "none"}
            for _t in files_under(_folder, "table"):
                if _t != NO_FILE and "panel" in os.path.basename(_t).lower():
                    _panels[_t] = _t
            _fields["ancestry"] = mo.ui.text(value=_ex.get("ancestry", ""), placeholder="e.g. EUR 0.55, AFR 0.30, AMR 0.15",
                                             label="The person's ancestry, for NAT2 (optional; superpopulations or CPIC groups):", full_width=True)
            _fields["ancestry_panel"] = mo.ui.dropdown(options=_panels, value=label_of(_panels, _ex.get("ancestry_panel", "none"), "(none: the plain mean of the CPIC population groups)"),
                                                       label="... or estimate it from the VCF with a reference panel (screen 17):", full_width=True)
            _submit = "💊 Call the pharmacogenes"
        else:
            _jars = [_j for _j in files_under(_folder, "jar") if "pharmcat" in _j.lower()]
            if os.path.isdir(PHARMCAT_DIR):
                _jars += [f"{PHARMCAT_DIR}/{_n}" for _n in sorted(os.listdir(PHARMCAT_DIR)) if _n.endswith("-all.jar")]
            _jars = list(dict.fromkeys(_jars)) or ["(PharmCAT not found: download it above)"]
            _outside = {"(none)": "none"}
            for _t in files_under(_folder, "table"):
                if _t != NO_FILE:
                    _outside[_t] = _t
            _fields["jar"] = mo.ui.dropdown(options=_jars, value=_jars[-1], label="PharmCAT (pharmcat-<version>-all.jar, with its positions file next to it):", full_width=True)
            _fields["outside"] = mo.ui.dropdown(options=_outside, value=label_of(_outside, _ex.get("outside", "none"), "(none)"),
                                                label="Outside calls (TSV), e.g. CYP2D6 or HLA-B from MMGEN's engine or a laboratory:", full_width=True)
            _submit = "💊 Run PharmCAT"
        _fields["output"] = mo.ui.text(value=_ex.get("output", ""), placeholder=f"blank: <VCF name>_{'pgx' if _mine else 'pharmcat'}/ next to the VCF",
                                       label="Output folder:", full_width=True)
        form_pgx = make_form(_fields, _submit)
    mo.vstack([_x for _x in (_about, form_pgx) if _x is not None])
    return (form_pgx,)


@app.cell
async def pgx_page_backend(asyncio, beside, chosen, code_panel, error_panel, facts, figs, figure_panel, form_pgx, get_page, mmg, mo, os, pgx_engine, pgx_hla_dir, pgx_tables_dir, python_call, re, split_name, warn):
    mo.stop(get_page() != "pgx")
    mo.stop(pgx_engine.value.startswith("A group"))  # the cell below runs that one
    mo.stop(form_pgx.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_pgx.value
    _mine = pgx_engine.value.startswith("MMGEN")
    mo.stop(not chosen(_v["vcf"]), warn("Select the person's VCF."))
    _kw = {}
    if (_v["sample"] or "").strip():
        _kw["sample"] = _v["sample"].strip()
    if _v["callable"] != "none":
        _kw["callable_bed"] = _v["callable"]
    if _mine:
        mo.stop("bam" not in _v, warn("The form is for the other engine: submit it again."))
        _tables = pgx_tables_dir()
        mo.stop(_tables is None, warn("MMGEN's engine needs its tables: click **Get the mmgen-pgx tables** above."))
        _kw["engine"] = "mmgen-pgx"
        _kw["tables"] = _tables
        if _v["bam"] != "none":
            _hla = pgx_hla_dir()
            mo.stop(_hla is None, warn("With the reads (BAM), the HLA genes need their reference sequences: click **Download the HLA "
                                       "reference sequences** above (68 MB, once)."))
            os.environ["MMGEN_HLA_DIR"] = _hla  # where the engine reads them
            _kw["bam"] = _v["bam"]
        if _v["reference"] != "none":
            _kw["reference"] = _v["reference"]
        _anc_text = (_v.get("ancestry") or "").strip()
        if _anc_text:
            _anc = {}
            for _part in re.split(r"[,;\n]", _anc_text):
                if not _part.strip():
                    continue
                _m2 = re.match(r"\s*(.+?)\s*[:=]?\s*([0-9]*\.?[0-9]+)\s*%?\s*$", _part)
                mo.stop(_m2 is None, warn(f"Ancestry: `{_part.strip()}` is not a group and a proportion (e.g. `EUR 0.6`)."))
                _anc[_m2.group(1)] = float(_m2.group(2))
            _kw["ancestry"] = _anc
        elif _v.get("ancestry_panel", "none") != "none":
            _kw["ancestry_panel"] = _v["ancestry_panel"]
        _what, _secs = "Calling the pharmacogenes", "a few seconds"
    else:
        mo.stop("jar" not in _v, warn("The form is for the other engine: submit it again."))
        mo.stop(not _v["jar"].endswith(".jar"), warn("PharmCAT is needed: click **Download PharmCAT** above, or put pharmcat-<version>-all.jar and its positions file in the data folder."))
        _kw["engine"] = "pharmcat"
        _kw["jar"] = _v["jar"]
        if _v["callable"] == "none":
            _kw["absent_is_ref"] = True
        if _v["outside"] != "none":
            _kw["outside_calls"] = _v["outside"]
        _what, _secs = "Preparing the input and running PharmCAT", "about 15 seconds"
    _out = (_v["output"] or "").strip() or beside(_v["vcf"], split_name(_v["vcf"])[0] + ("_pgx" if _mine else "_pharmcat"))
    try:
        with mo.status.spinner(f"{_what} ({_secs})..."):
            _r = await asyncio.to_thread(mmg.pgx_run, _v["vcf"], _out, **_kw)
        _genes = _r["genes"]
        _called = [_g for _g in _genes if _g["called"]]
        _rows = [
            ("Engine", f"{_r['engine']} {_r['engine_version']}; tables: {_r['data_version']}"),
            ("Sample", f"`{_r['sample']}`"),
            ("Genes called", f"`{len(_called)}` of `{len(_genes)}`"),
            ("Drug recommendations", f"`{len(_r['drugs'])}`"),
        ]
        _anc_used = _r.get("x_ancestry")
        if _anc_used:
            _rows.append(("Ancestry used for NAT2", ", ".join(f"{_k} {_x:.2f}" for _k, _x in _anc_used["proportions"].items() if _x >= 0.005)
                          + f" ({_anc_used['source']})"))
        _p = (_r.get("x_pharmcat") or {}).get("prepared")
        if _p:
            _read = _p["with_variant"] + _p["reference_record"] + _p["reference_filled"]
            _rows.insert(2, ("Positions read", f"`{_read:,}` of `{_p['positions']:,}` (`{_p['with_variant']:,}` with a variant)"))
        _items = [mo.callout(mo.md("### 💊 Pharmacogenomics results\n" + facts(_rows)), kind="success")]
        _has_cn = any(_g.get("copy_number") is not None for _g in _genes)
        # mmgen-pgx picks the highest-scoring diplotype (as PharmCAT's Named Allele Matcher) and
        # lists the others that fit the positions read; a warning names the ones that would
        # change the phenotype or a recommendation.
        _has_alt = any(_g.get("x_alternatives") for _g in _called)
        _table = []
        for _g in _called:
            _row = {"gene": _g["gene"], "diplotype": _g["diplotype"], "phenotype": _g["phenotype"] or "",
                    "activity score": "" if _g["activity_score"] is None else f"{_g['activity_score']:g}"}
            if _has_cn:
                _row["copy number"] = "" if _g.get("copy_number") is None else f"{_g['copy_number']:g}"
            if _has_alt:
                _alts = _g.get("x_alternatives") or []
                _row["also possible"] = "; ".join(f"{_a['diplotype']} ({_a['phenotype'] or '-'})" for _a in _alts[:3]) + (
                    f"; and {len(_alts) - 3} more" if len(_alts) > 3 else "")
            _row["from"] = _g["source"]
            _row["positions read"] = _g["positions_read"]
            _table.append(_row)
        _items.append(mo.ui.table(_table, selection=None, pagination=True, page_size=25, label="Genes called"))
        _uncalled = [_g for _g in _genes if not _g["called"]]
        if _uncalled:
            _items.append(mo.ui.table([{"gene": _g["gene"], "why it was not called": _g["not_called_reason"]} for _g in _uncalled],
                                      selection=None, pagination=True, page_size=10, label="Genes not called, and why"))
        if _r["traits"]:
            _items.append(mo.ui.table([{"trait": _t["trait"], "genotype": _t["genotype"], "what it suggests": _t["interpretation"],
                                        "evidence": _t["evidence"]} for _t in _r["traits"]],
                                      selection=None, label="Traits (for teaching)"))
        if _r["drugs"]:
            _items.append(mo.ui.table([{"drug": _d["drug"], "genes": ", ".join(_d["genes"]), "source": _d["source"],
                                        "strength": _d["strength"] or "", "recommendation": _d["recommendation"]} for _d in _r["drugs"]],
                                      selection=None, pagination=True, page_size=10, label="Guideline recommendations"))
        for _w in _r["warnings"]:
            _items.append(warn(_w))
        _oc = (_r.get("x_files") or {}).get("pharmcat_outside_calls")
        if _oc and os.path.exists(_oc):
            _items.append(mo.callout(mo.md(
                f"PharmCAT can use the calls this engine made from the reads: run it with the outside-call file `{_oc}` "
                "(the **PharmCAT** engine above, field *Outside calls*)."), kind="info"))
        _how = {}
        for _g in _called:
            if _g["evidence"]:
                _how[f"{_g['gene']}: {_g['diplotype']}"] = mo.ui.table(
                    [{"position": f"{_e['chrom']}:{_e['pos']:,}", "rsID": _e["rsid"] or "", "genotype": _e["genotype"],
                      "in callable regions": "" if _e["in_callable"] is None else ("yes" if _e["in_callable"] else "no")}
                     for _e in _g["evidence"]], selection=None, pagination=True, page_size=10)
        if _how:
            _items.append(mo.md("#### How each call was made: the positions read"))
            _items.append(mo.accordion(_how))
        # CYP2D6 from the reads: the copy number along the gene and its pseudogene, against diploid control regions.
        _tabdir_cn = _kw.get("tables")
        if _kw.get("bam") and _tabdir_cn and hasattr(mmg, "depth_profile") and os.path.exists(os.path.join(_tabdir_cn, "bam_genes.json")):
            try:
                import json as _json_cn
                with open(os.path.join(_tabdir_cn, "bam_genes.json"), encoding="utf-8") as _fh:
                    _cnv = _json_cn.load(_fh)["cyp2d6"]["cnv_regions"]
                _w, _p = _cnv["target_whole"], _cnv["paralog"]
                _lo, _hi = min(_w["start"], _p["start"]) - 3000, max(_w["end"], _p["end"]) + 3000
                with mo.status.spinner("Measuring the depth over CYP2D6, CYP2D7 and the control regions..."):
                    _loc = await asyncio.to_thread(mmg.depth_profile, _kw["bam"], f"{_w['chrom']}:{_lo}-{_hi}", 500)
                    _ctl = await asyncio.to_thread(mmg.depth_profile, _kw["bam"],
                                                   [f"{_c['chrom']}:{_c['start']}-{_c['end']}" for _c in _cnv["controls"]], 1000)
                _means = sorted(_r0["mean"] for _r0 in _ctl["regions"] if _r0.get("reads"))
                if _means:
                    _med = _means[len(_means) // 2] if len(_means) % 2 else (_means[len(_means) // 2 - 1] + _means[len(_means) // 2]) / 2
                    _bins = [tuple(_b0) for _b0 in _loc["regions"][0]["bins"]]
                    _items.append(figure_panel(
                        figs.copy_number_profile(_bins, _med, features=[
                            ("CYP2D6", _w["start"], _w["end"]), ("exon 9", _cnv["target_ex9"]["start"], _cnv["target_ex9"]["end"]),
                            ("CYP2D7 (pseudogene)", _p["start"], _p["end"])],
                            title=f"CYP2D6 copies along the gene (2 × depth / median of {len(_means)} control regions, {_med:.0f}×)"),
                        "cyp2d6_copy_number",
                        "About 2 copies all along CYP2D6: no deletion or duplication. The whole gene at 1 copy: one copy deleted "
                        "(*5); at 3 or more: a duplication (×N). Only part of the gene changed — the exon 9 end apart from the "
                        "5' end — is a hybrid with CYP2D7 (*36, *68). The CYP2D7 pseudogene is the comparison next to it."))
            except Exception as _cerr:
                _items.append(mo.md(f"*The CYP2D6 copy-number plot could not be drawn: {_cerr}*"))
        # Genes with alternatives: the positions that tell the candidates apart, and what was read there.
        _amb = [_g for _g in _called if _g.get("x_alternatives")]
        _tabdir = _kw.get("tables")
        if _amb and _tabdir and os.path.exists(os.path.join(_tabdir, "genes.json")):
            import json as _json_pgx
            with open(os.path.join(_tabdir, "genes.json"), encoding="utf-8") as _fh:
                _gtab = {_x["gene"]: _x for _x in _json_pgx.load(_fh)["genes"]}
            _mat = {}
            for _g in _amb:
                if _g["gene"] in _gtab:
                    _mat[f"{_g['gene']}: {_g['diplotype']} or " + ", ".join(_a["diplotype"] for _a in _g["x_alternatives"][:4])] = figure_panel(
                        figs.pgx_evidence(_g, _gtab[_g["gene"]]), f"pgx_evidence_{_g['gene']}")
            if _mat:
                _items.append(mo.md("#### Why some calls have alternatives: the evidence matrix\n\n"
                                    "*Rows: the positions that define the candidate alleles; a filled cell is a base the allele "
                                    "requires. The last column is what was read. Unphased variants, or positions not read, leave "
                                    "more than one pair of alleles possible.*"))
                _items.append(mo.accordion(_mat))
        _items.append(mo.callout(mo.md(
            "**Reading the result.** A phenotype such as *poor metabolizer* says how fast an enzyme is expected to work, not what to "
            "prescribe. **Any change to a treatment is a decision for the prescriber, with a laboratory-confirmed result.**"), kind="warn"))
        _html = _r.get("report_html")
        if _html and os.path.exists(_html):
            with open(_html, encoding="utf-8") as _fh:
                _items.append(mo.accordion({"The engine's full report (HTML)": mo.iframe(_fh.read(), height="900px")}))
        _items.append(code_panel(f"r = {python_call('pgx_run', _v['vcf'], _out, **_kw)}",
                                 'for g in r["genes"]:\n    print(g["gene"], g["diplotype"] if g["called"] else "not called: " + g["not_called_reason"])'))
        _panel = mo.vstack(_items)
    except Exception as _err:
        _panel = error_panel(_err)
    _panel
    return


@app.cell
async def pgx_group_backend(asyncio, bed_next_to, beside, chosen, code_panel, error_panel, facts, form_pgx, get_page, mmg, mo, pgx_engine, pgx_tables_dir, python_call, re, warn, words):
    mo.stop(get_page() != "pgx" or not pgx_engine.value.startswith("A group"))
    mo.stop(form_pgx.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_pgx.value
    mo.stop("vcfs" not in _v, warn("The form is for one person: submit it again."))
    _vcfs = [_p for _p in _v["vcfs"] if chosen(_p)]
    mo.stop(not _vcfs, warn("Select the people's VCFs: one per person, or a multi-sample VCF."))
    _tables = pgx_tables_dir()
    mo.stop(_tables is None, warn("MMGEN's engine needs its tables: click **Get the mmgen-pgx tables** above."))
    _kw = {"tables": _tables}
    if _v["callable"] == "auto":
        _found = [bed_next_to(_p) for _p in _vcfs]
        if any(_found):
            _kw["callable_beds"] = [_b or "" for _b in _found]
    elif _v["callable"] != "none":
        _kw["callable_beds"] = _v["callable"]
    if words(_v["samples"]):
        _kw["samples"] = words(_v["samples"])
    _drugs = [_d.strip() for _d in re.split(r"[,;\n]", _v["drugs"] or "") if _d.strip()]
    if _drugs:
        _kw["drugs"] = _drugs
    if _v.get("ancestry_panel", "none") != "none":
        _kw["ancestry_panel"] = _v["ancestry_panel"]
    _out = (_v["output"] or "").strip() or beside(_vcfs[0], "pgx_cohort")
    try:
        with mo.status.spinner(f"Calling the pharmacogenes of {len(_vcfs)} VCF(s) (about a second per person)..."):
            _r = await asyncio.to_thread(mmg.pgx_cohort, _vcfs, _out, **_kw)
        _people = [_p["person"] for _p in _r["people"]]
        _meds = _r["medicines"]
        _n_rec = sum(1 for _m in _meds if _m["status"] == "recommendation")
        _items = [mo.callout(mo.md("### 💊 Pharmacogenomics of a group\n" + facts([
            ("Engine", f"{_r['engine']} {_r['engine_version']}; tables: {_r['data_version']}"),
            ("People", f"`{len(_people)}`: " + ", ".join(f"`{_p}`" for _p in _people[:12]) + (" ..." if len(_people) > 12 else "")),
            ("Medicines asked for", f"`{len(_r['medicines_requested'])}`" if _r["medicines_requested"] else "every medicine with a recommendation"),
            ("Recommendations", f"`{_n_rec}`"),
            ("Tables", f"`{_r['files']['medicines']}` and `{_r['files']['genes']}`"),
        ])), kind="success")]
        if _r["medicines_without_guideline"]:
            _items.append(mo.callout(mo.md(
                "**No CPIC guideline in the tables for:** " + ", ".join(f"`{_d}`" for _d in _r["medicines_without_guideline"])
                + ". Their effect does not depend on a gene with a guideline (or the name is spelt differently): genotype "
                "cannot guide them, which is also worth knowing when an adverse event is reviewed."), kind="info"))
        _mark = {"recommendation": "✅", "gene not called": "⚪ gene not called", "no recommendation": "➖ no recommendation",
                 "no guideline": "∅ no guideline"}
        if _r["medicines_requested"]:
            # One row per medicine, one column per person: what each person's genotype says.
            _grid = []
            for _d in _r["medicines_requested"]:
                _row = {"medicine": _d}
                for _p in _people:
                    _ms = [_m for _m in _meds if _m["person"] == _p and _m["medicine"] == _d]
                    if not _ms:
                        _row[_p] = ""
                    elif _ms[0]["status"] == "recommendation":
                        _row[_p] = "✅ " + "; ".join(dict.fromkeys(str(_m["phenotype"] or "") for _m in _ms))
                    else:
                        _row[_p] = _mark.get(_ms[0]["status"], _ms[0]["status"])
                _grid.append(_row)
            _items.append(mo.ui.table(_grid, selection=None, pagination=True, page_size=20,
                                      label="Each medicine for each person (✅: there is a recommendation for this genotype)"))
        _items.append(mo.ui.table([{
            "person": _m["person"], "medicine": _m["medicine"], "status": _mark.get(_m["status"], _m["status"]),
            "genes": ", ".join(_m["genes"]), "phenotype": _m["phenotype"] or "", "recommendation or reason": _m["recommendation"] or "",
            "strength": _m["strength"] or "",
        } for _m in _meds], selection=None, pagination=True, page_size=10, label="The recommendations, and why some are missing"))
        _relevant = sorted({_g for _m in _meds for _g in _m["genes"]}) if _r["medicines_requested"] else None
        _genes_grid = {}
        for _g in _r["genes"]:
            if _relevant is not None and _g["gene"] not in _relevant:
                continue
            _cell = (f"{_g['diplotype']} ({_g['phenotype'] or '-'})" if _g["called"] else "not called")
            if _g["also_possible_changing"]:
                _cell += " ⚠️ also possible: " + "; ".join(_g["also_possible_changing"][:2])
            _genes_grid.setdefault(_g["gene"], {"gene": _g["gene"]})[_g["person"]] = _cell
        if _genes_grid:
            _items.append(mo.ui.table(list(_genes_grid.values()), selection=None, pagination=True, page_size=20,
                                      label="Each gene for each person" + (" (the genes of these medicines)" if _relevant is not None else "")))
        _anc_rows = [{"person": _p["person"], **{_k: round(_x, 3) for _k, _x in (_p.get("ancestry") or {}).items()}}
                     for _p in _r["people"] if _p.get("ancestry")]
        if _anc_rows:
            _items.append(mo.ui.table(_anc_rows, selection=None, label="Each person's ancestry (it weights the CPIC population groups for NAT2)"))
        if _r["warnings"]:
            _items.append(mo.accordion({f"Warnings ({len(_r['warnings'])})": mo.md("\n".join(f"* {_w}" for _w in _r["warnings"][:60]))}))
        _items.append(mo.callout(mo.md(
            "**Reading the result.** CYP2D6 (codeine, tramadol, ondansetron, tricyclics) needs the person's reads (a BAM) or a "
            "laboratory result: from a VCF alone it is *not called*. A recommendation says what CPIC advises for a genotype, not "
            "what to prescribe: **any change to a treatment is a decision for the prescriber, with a laboratory-confirmed result.**"),
            kind="warn"))
        _items.append(code_panel(f"r = {python_call('pgx_cohort', _vcfs, _out, **_kw)}",
                                 'for m in r["medicines"]:\n    print(m["person"], m["medicine"], m["status"], m["phenotype"] or "")'))
        _panel = mo.vstack(_items)
    except Exception as _err:
        _panel = error_panel(_err)
    _panel
    return


@app.cell
def lookup_page_ui(get_page, human_notice, mo, set_page):
    mo.stop(get_page() != "lookup")
    btn_back_lookup = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_lookup = mo.ui.button(label="🔄 Refresh the file list")
    lookup_task = mo.ui.radio(options=["Look up variants in relatives (cascade testing)", "Allele counts at a position (BAM)"],
                              value="Look up variants in relatives (cascade testing)", inline=True, label="**What to do:**")
    lookup_folder = mo.ui.text(value=".", label="Data folder (its subfolders are searched too):", full_width=True)
    _ui = mo.vstack([
        btn_back_lookup,
        mo.md("<br>\n\n# 🔎 21. Variant Lookup & Allele Counts"),
        human_notice(),
        mo.callout(mo.md("Once a variant is known in a family, relatives are tested for that one variant (**cascade testing**), usually by Sanger sequencing or qPCR. Here: the genotype of each person at given variants, from their VCFs; and, from a BAM, **the bases the reads show at one position** (depth and allele fractions) — what one looks at to confirm a variant."), kind="info"),
        lookup_folder,
        mo.hstack([lookup_task, btn_refresh_lookup], justify="space-between", align="center"),
    ])
    _ui
    return btn_back_lookup, btn_refresh_lookup, lookup_folder, lookup_task


@app.cell
def lookup_page_example(example_picker, get_page, lookup_task, mo):
    mo.stop(get_page() != "lookup")
    lookup_example = example_picker("lookup", lookup_task.value)
    lookup_example
    return (lookup_example,)


@app.cell
def lookup_page_form(NO_FILE, btn_refresh_lookup, example_for, files_under, get_page, label_of, lookup_example, lookup_folder, lookup_task, make_form, mo):
    mo.stop(get_page() != "lookup")
    btn_refresh_lookup  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("lookup", lookup_task.value, lookup_example.value)
    _folder = lookup_folder.value
    if lookup_task.value.startswith("Look up"):
        _vcfs = files_under(_folder, "vcf")
        _cv = {"(none: coordinates only)": "none"}
        for _f in _vcfs:
            if _f != NO_FILE and "clinvar" in _f.lower():
                _cv[_f] = _f
        form_lookup = make_form({
            "variants": mo.ui.text_area(value=_ex.get("variants", ""), placeholder="chr4:186274193:G:T\nrs121965063",
                                        label="Variants, one per line (chrom:pos:ref:alt; or an rsID / ClinVar ID with ClinVar below):", full_width=True),
            "vcfs": mo.ui.multiselect(options=_vcfs, value=_ex.get("vcfs", []), label="The relatives' VCFs:", full_width=True),
            "use_beds": mo.ui.checkbox(value=_ex.get("use_beds", True), label="Use the callable-regions BED next to each VCF (`x.callable.bed`, or the Genome in a Bottle `*_benchmark*.bed`)"),
            "clinvar": mo.ui.dropdown(options=_cv, value=label_of(_cv, _ex.get("clinvar", "none"), "(none: coordinates only)"), label="ClinVar VCF, to resolve rsIDs and show the classification:", full_width=True),
        }, "🔎 Look up")
    else:
        form_lookup = make_form({
            "bam": mo.ui.dropdown(options=files_under(_folder, "sam"), value=_ex.get("bam"), label="Aligned reads (BAM with .bai, or SAM):", full_width=True),
            "position": mo.ui.text(value=_ex.get("position", ""), placeholder="chr1:25000", label="Position (1-based):"),
            "alt": mo.ui.text(value=_ex.get("alt", ""), placeholder="blank: the most frequent other allele; or T, C>T, CAG>C", label="Allele to measure:"),
            "min_mapq": mo.ui.number(start=0, stop=255, step=1, value=_ex.get("min_mapq", 0), label="Minimum mapping quality:"),
            "min_baseq": mo.ui.number(start=0, stop=93, step=1, value=_ex.get("min_baseq", 13), label="Minimum base quality:"),
            "reference": mo.ui.dropdown(options={"(none: mismatches from the reads' MD tags)": "none",
                                                 **{_f: _f for _f in files_under(_folder, "fasta") if _f != NO_FILE}},
                                        value=label_of({"(none: mismatches from the reads' MD tags)": "none",
                                                        **{_f: _f for _f in files_under(_folder, "fasta") if _f != NO_FILE}},
                                                       _ex.get("reference", "none"), "(none: mismatches from the reads' MD tags)"),
                                        label="Reference genome (FASTA), for the browser view:", full_width=True),
            "annotation": mo.ui.dropdown(options={"(none)": "none", **{_f: _f for _f in files_under(_folder, "gff") if _f != NO_FILE}},
                                         value=label_of({"(none)": "none", **{_f: _f for _f in files_under(_folder, "gff") if _f != NO_FILE}},
                                                        _ex.get("annotation", "none"), "(none)"),
                                         label="Genes (GFF/GTF), drawn above the reads:", full_width=True),
        }, "🧮 Count")
    mo.vstack([_x for _x in (_about, form_lookup) if _x is not None])
    return (form_lookup,)


@app.cell
async def lookup_page_backend(asyncio, bed_next_to, chosen, code_panel, error_panel, facts, figs, figure_panel, form_lookup, get_page, lookup_task, mmg, mo, python_call, warn):
    mo.stop(get_page() != "lookup")
    mo.stop(form_lookup.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_lookup.value

    if lookup_task.value.startswith("Look up"):
        _vars = [_l.strip() for _l in (_v["variants"] or "").splitlines() if _l.strip()]
        _vcfs = [_p for _p in _v["vcfs"] if chosen(_p)]
        mo.stop(not _vars, warn("Type at least one variant."))
        mo.stop(not _vcfs, warn("Select at least one VCF."))
        _kw = {}
        if _v["use_beds"]:
            _beds = [bed_next_to(_p) for _p in _vcfs]
            if any(_beds):
                _kw["callable_beds"] = _beds
        if _v["clinvar"] != "none":
            _kw["clinvar"] = _v["clinvar"]
        try:
            with mo.status.spinner("Looking the variants up..."):
                _r = await asyncio.to_thread(mmg.variant_lookup, _vars, _vcfs, **_kw)
            _meaning = {"present": "carries it", "other_allele": "another allele here", "no_call": "no call", "reference": "reference (called)",
                        "absent_callable": "does not carry it (reliable region)", "not_callable": "not read (outside the reliable region)",
                        "absent": "not in the file"}
            _table = []
            for _c in _r["calls"]:
                _x = _r["variants"][_c["variant"]]
                _table.append({
                    "variant": _x["query"] if _x["query"].count(":") == 3 else f"{_x['query']} ({_x['chrom']}:{_x['pos']} {_x['ref']}>{_x['alt']})",
                    "person": _r["people"][_c["person"]],
                    "genotype": _c.get("genotype") or "",
                    "zygosity": _c.get("zygosity") or "",
                    "result": _meaning.get(_c["status"], _c["status"]),
                })
            _info = []
            for _x in _r["variants"]:
                _desc = f"* `{_x['query']}`"
                if _x.get("in_clinvar"):
                    _desc += f": {', '.join(_x.get('genes') or [])} — ClinVar {(_x.get('significance') or '').replace('_', ' ')} ({'★' * (_x.get('stars') or 0)})"
                _info.append(_desc)
            _items = [mo.callout(mo.md("### 🔎 Lookup\n" + "\n".join(_info) + (("\n\nNot resolved: " + ", ".join(f"`{_u}`" for _u in _r["unresolved"])) if _r["unresolved"] else "")), kind="success")]
            _items.append(mo.ui.table(_table, selection=None, pagination=True, page_size=20, label="Each person at each variant"))
            if not _kw.get("callable_beds"):
                _items.append(mo.md("*Without a callable-regions BED, a variant absent from a VCF may be reference or simply not read.*"))
            _items.append(code_panel(f"r = {python_call('variant_lookup', _vars, _vcfs, **_kw)}"))
            _panel = mo.vstack(_items)
        except Exception as _err:
            _panel = error_panel(_err)
    else:
        mo.stop(not chosen(_v["bam"]), warn("Select a BAM (or SAM) file."))
        mo.stop(":" not in (_v["position"] or ""), warn("Type the position as chrom:pos, e.g. chr1:25000."))
        _kw = {}
        if (_v["alt"] or "").strip():
            _kw["alt"] = _v["alt"].strip()
        if int(_v["min_mapq"] or 0):
            _kw["min_mapq"] = int(_v["min_mapq"])
        if int(_v["min_baseq"] or 0) != 13:
            _kw["min_baseq"] = int(_v["min_baseq"] or 0)
        try:
            with mo.status.spinner("Counting the bases..."):
                _r = await asyncio.to_thread(mmg.pileup, _v["bam"], _v["position"].strip(), **_kw)
            _rows = [
                ("Reads over the position", f"`{_r['reads']:,}`; left out: flags `{_r['filtered_flag']}`, mapping quality `{_r['filtered_mapq']}`, base quality `{_r['filtered_baseq']}`"),
                ("Depth", f"`{_r['depth']:,}`"),
                ("Bases", f"A `{_r['a']}` · C `{_r['c']}` · G `{_r['g']}` · T `{_r['t']}` · N `{_r['n']}` · deletion `{_r['deletions']}`"),
                ("Reference allele", f"`{_r['ref_allele']}`: `{_r['ref_count']}` reads ({100 * (_r['ref_fraction'] or 0):.1f}%)" if _r.get("ref_allele") else "unknown"),
            ]
            if _r.get("alt"):
                _rows.append(("Allele measured", f"`{_r['alt']}`: `{_r['alt_count']}` reads ({100 * (_r['alt_fraction'] or 0):.1f}%)"))
            _items = [mo.callout(mo.md(f"### 🧮 {_v['position']}\n" + facts(_rows)), kind="success" if _r["depth"] else "warn")]
            if _r["alleles"]:
                _items.append(mo.ui.table([{
                    "allele": _a.get("allele", ""), "reads": _a.get("count", 0), "forward": _a.get("forward", 0),
                    "reverse": _a.get("reverse", 0), "fraction": f"{100 * _a.get('fraction', 0):.1f}%",
                } for _a in _r["alleles"]], selection=None, label="Alleles"))
            _items.append(mo.md("*About 50% of the reads with an allele is what a heterozygous variant looks like, about 100% a homozygous one; a few percent is usually noise — or, in a tumor, a low-fraction somatic variant.*"))
            _code = [f"r = {python_call('pileup', _v['bam'], _v['position'].strip(), **_kw)}"]
            # The reads themselves, as a genome browser shows them.
            if hasattr(mmg, "region_reads"):
                _kw_v = {"flank": 60}
                if _v["reference"] != "none":
                    _kw_v["reference"] = _v["reference"]
                if int(_v["min_mapq"] or 0):
                    _kw_v["min_mapq"] = int(_v["min_mapq"])
                if int(_v["min_baseq"] or 0) != 13:
                    _kw_v["min_baseq"] = int(_v["min_baseq"] or 0)
                _view = await asyncio.to_thread(mmg.region_reads, _v["bam"], _v["position"].strip(), **_kw_v)
                _feats = []
                if _v["annotation"] != "none":
                    import gzip as _gz_ann
                    _reg = _view["region"]
                    _chrom_alias = {_reg["chrom"], _reg["chrom"][3:] if _reg["chrom"].startswith("chr") else "chr" + _reg["chrom"]}
                    _op = _gz_ann.open if _v["annotation"].endswith(".gz") else open
                    with _op(_v["annotation"], "rt", encoding="utf-8", errors="replace") as _fh:
                        for _line in _fh:
                            if _line.startswith("#"):
                                continue
                            _f = _line.rstrip("\n").split("\t")
                            if len(_f) < 9 or _f[0] not in _chrom_alias or _f[2] not in ("gene", "exon", "CDS"):
                                continue
                            _fs, _fe = int(_f[3]), int(_f[4])
                            if _fe < _reg["start"] or _fs > _reg["end"]:
                                continue
                            _attrs = dict(_x.split("=", 1) for _x in _f[8].split(";") if "=" in _x)
                            _feats.append((_attrs.get("Name") or _attrs.get("gene") or _attrs.get("ID") or _f[2], _fs, _fe, _f[6], _f[2]))
                _pos = int(_v["position"].strip().split(":")[1].replace(",", "").split("-")[0])
                _sb = figs.strand_and_position(_view, _pos, (_v["alt"] or "").strip().upper()[-1:] or None)
                _items.append(figure_panel(figs.genome_browser(_view, features=_feats), "reads_at_position",
                                           "Hover a read, a base or a coverage bar for its details. Reads with mapping quality 0 are faint; "
                                           "so are mismatches in bases of quality below 20."))
                if _sb.get("alt"):
                    _rows_sb = [("Allele checked", f"`{_sb['ref']}>{_sb['alt']}` at {_pos:,}")]
                    if _sb.get("alt_forward") is not None:
                        _rows_sb.append(("Reads with it: forward / reverse", f"`{_sb['alt_forward']}` / `{_sb['alt_reverse']}`; with the reference: `{_sb['ref_forward']}` / `{_sb['ref_reverse']}`"))
                    if _sb.get("fisher_p") is not None:
                        _rows_sb.append(("Strand bias", f"Fisher's p = `{_sb['fisher_p']:.3g}`, SOR = `{_sb['sor']:.2f}` (GATK filters SOR > 3 for SNVs)"))
                    if _sb.get("alt_reads_drawn"):
                        _rows_sb.append(("Where in the reads", f"median `{_sb['alt_median_distance_to_end']:.0f}` bases from the nearest end; "
                                                               f"`{100 * _sb['alt_within_10_of_end']:.0f}%` within 10 bases of an end ({_sb['alt_reads_drawn']} reads drawn)"))
                    _bad = (_sb.get("sor") or 0) > 3 or (_sb.get("alt_within_10_of_end") or 0) > 0.5
                    _items.append(mo.callout(mo.md("### 🔬 Is it real? The usual checks\n" + facts(_rows_sb) +
                                                   ("\n\n⚠️ **One of the checks points to an artefact**: confirm it independently." if _bad else "")),
                                             kind="warn" if _bad else "info"))
                _code.append(f"view = {python_call('region_reads', _v['bam'], _v['position'].strip(), **_kw_v)}")
            _items.append(code_panel(*_code))
            _panel = mo.vstack(_items)
        except Exception as _err:
            _panel = error_panel(_err)
    _panel
    return


@app.cell
def assoc_page_ui(get_page, human_notice, mo, set_page):
    mo.stop(get_page() != "assoc")
    btn_back_assoc = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_assoc = mo.ui.button(label="🔄 Refresh the file list")
    assoc_task = mo.ui.radio(options=["Test variants against an outcome", "Study power: how many people?"],
                             value="Test variants against an outcome", inline=True, label="**What to do:**")
    assoc_folder = mo.ui.text(value=".", label="Data folder (its subfolders are searched too):", full_width=True)
    _ui = mo.vstack([
        btn_back_assoc,
        mo.md("<br>\n\n# 📐 30. Genotype & Outcome: Association and Power"),
        human_notice(),
        mo.callout(mo.md(
            "Does a variant go with an outcome in a group of people — a better recovery, a side effect, a measurement? "
            "Each allele is tested against the outcome, **adjusted for covariates** (age, sex, center, ancestry...): "
            "linear regression for a number, logistic regression for yes/no, and **Firth's method** when a small or "
            "separated sample breaks the plain fit. Many variants mean many tests: the **FDR** and, with permutations, "
            "a **family-wise p** say how many hits chance alone would give.\n\n"
            "**Study power** answers the question to ask first: with this many people, which effect could the study "
            "detect at all? An association is not causation. Relatives inflate the results (check with screen 16), and "
            "population structure is the classic source of false hits (screen 17 gives the ancestry proportions to "
            "use as covariates)."), kind="info"),
        assoc_folder,
        mo.hstack([assoc_task, btn_refresh_assoc], justify="space-between", align="center"),
    ])
    _ui
    return assoc_folder, assoc_task, btn_back_assoc, btn_refresh_assoc


@app.cell
def assoc_page_example(assoc_task, example_picker, get_page, mo):
    mo.stop(get_page() != "assoc")
    assoc_example = example_picker("assoc", assoc_task.value)
    assoc_example
    return (assoc_example,)


@app.cell
def assoc_page_form(assoc_example, assoc_folder, assoc_task, btn_refresh_assoc, example_for, files_under, get_page, label_of, make_form, mo):
    mo.stop(get_page() != "assoc")
    btn_refresh_assoc  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("assoc", assoc_task.value, assoc_example.value)
    _folder = assoc_folder.value
    _models = {"Additive: 0, 1 or 2 copies": "additive", "Dominant: any copy": "dominant", "Recessive: two copies": "recessive"}
    if assoc_task.value.startswith("Test"):
        _types = {"Automatic (two values: yes/no)": "auto", "Yes/no (binary)": "binary", "A number (quantitative)": "quantitative"}
        _firth = {"When the plain fit fails or separates (recommended)": "auto", "Always": "always", "Never": "never"}
        form_assoc = make_form({
            "vcfs": mo.ui.multiselect(options=files_under(_folder, "vcf"), value=_ex.get("vcfs", []),
                                      label="Genotypes: one VCF with everyone, or one VCF per person:", full_width=True),
            "phenotypes": mo.ui.dropdown(options=files_under(_folder, "table"), value=_ex.get("phenotypes"),
                                         label="Phenotype table (a header; one row per person, named as in the VCF):", full_width=True),
            "outcome": mo.ui.text(value=_ex.get("outcome", ""), placeholder="recovery_score", label="Outcome column:", full_width=True),
            "outcome_type": mo.ui.dropdown(options=_types, value=label_of(_types, _ex.get("outcome_type", "auto"), "Automatic (two values: yes/no)"),
                                           label="Kind of outcome:", full_width=True),
            "case_value": mo.ui.text(value=_ex.get("case_value", ""), placeholder="blank: automatic (e.g. yes, 1, case)",
                                     label="Yes/no outcome: the value that counts as the event:", full_width=True),
            "covariates": mo.ui.text(value=_ex.get("covariates", ""), placeholder="age, sex, center, anc_afr, anc_eas",
                                     label="Covariates to adjust for (columns; text columns become categories):", full_width=True),
            "model": mo.ui.dropdown(options=_models, value=label_of(_models, _ex.get("model", "additive"), "Additive: 0, 1 or 2 copies"),
                                    label="Genetic model:", full_width=True),
            "variants": mo.ui.text_area(value=_ex.get("variants", ""), placeholder="blank: every variant that passes the filters\nrs12190259\nchr6:154039662",
                                        label="Only these variants (rsID, chr:pos or chr:pos:ref:alt, one per line):", full_width=True),
            "regions": mo.ui.text(value=_ex.get("regions", ""), placeholder="blank: everywhere; e.g. chr6:154000000-154100000",
                                  label="Only these regions (chr:start-end, separated by commas):", full_width=True),
            "min_maf": mo.ui.number(start=0, stop=0.5, step=0.005, value=_ex.get("min_maf", 0.01), label="Minimum minor-allele frequency:"),
            "min_call_rate": mo.ui.number(start=0, stop=1, step=0.05, value=_ex.get("min_call_rate", 0.9), label="Minimum call rate:"),
            "firth": mo.ui.dropdown(options=_firth, value=label_of(_firth, _ex.get("firth", "auto"), "When the plain fit fails or separates (recommended)"),
                                    label="Firth's penalized logistic regression (yes/no outcomes):", full_width=True),
            "permutations": mo.ui.number(start=0, stop=100000, step=100, value=_ex.get("permutations", 0),
                                         label="Permutations (0: none; 1,000 gives an empirical and a family-wise p):"),
            "absent_is_ref": mo.ui.checkbox(value=_ex.get("absent_is_ref", False),
                                            label="A variant missing from a person's VCF means reference (only for files that list every call of a region that was read)"),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: <phenotype table>.association.tsv",
                                 label="Results table (TSV):", full_width=True),
        }, "📐 Test")
    else:
        _kinds = {"Yes/no: the effect is an odds ratio per copy": "binary", "A number: the effect is in standard deviations per copy": "quantitative"}
        form_assoc = make_form({
            "outcome": mo.ui.dropdown(options=_kinds, value=label_of(_kinds, _ex.get("outcome", "binary"), "Yes/no: the effect is an odds ratio per copy"),
                                      label="Kind of outcome:", full_width=True),
            "n": mo.ui.number(start=2, stop=10_000_000, step=1, value=_ex.get("n", 100), label="People in the study:"),
            "allele_freq": mo.ui.number(start=0.001, stop=0.5, step=0.001, value=_ex.get("allele_freq", 0.25),
                                        label="Allele frequency (yes/no: among the people without the event):"),
            "effect": mo.ui.number(start=0, stop=100, step=0.05, value=_ex.get("effect", 2.0),
                                   label="Effect to detect (odds ratio, or SD per copy; 0: only the smallest detectable effect):"),
            "case_fraction": mo.ui.number(start=0.01, stop=0.99, step=0.01, value=_ex.get("case_fraction", 0.5),
                                          label="Yes/no: the fraction of people with the event:"),
            "alpha": mo.ui.text(value=str(_ex.get("alpha", 0.05)), label="Significance level (α; e.g. 0.05, or 5e-8 genome-wide):"),
            "tests": mo.ui.number(start=1, stop=100_000_000, step=1, value=_ex.get("tests", 1),
                                  label="Variants tested (Bonferroni: α divided by this):"),
            "model": mo.ui.dropdown(options=_models, value=label_of(_models, _ex.get("model", "additive"), "Additive: 0, 1 or 2 copies"),
                                    label="Genetic model:", full_width=True),
            "covariates": mo.ui.number(start=0, stop=100, step=1, value=_ex.get("covariates", 0), label="Covariates in the model:"),
            "target_power": mo.ui.number(start=0.05, stop=0.99, step=0.05, value=_ex.get("target_power", 0.8), label="Target power:"),
            "simulations": mo.ui.number(start=0, stop=20000, step=500, value=_ex.get("simulations", 0),
                                        label="Simulated studies (0: the formula only; 2,000 checks it in small studies):"),
        }, "🧮 Compute the power")
    mo.vstack([_x for _x in (_about, form_assoc) if _x is not None])
    return (form_assoc,)


@app.cell
async def assoc_page_backend(asyncio, assoc_task, chosen, code_panel, error_panel, facts, figs, figure_panel, form_assoc, get_page, mmg, mo, output_name, python_call, scatter_chart, warn, words):
    mo.stop(get_page() != "assoc")
    mo.stop(form_assoc.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_assoc.value

    def _num(x, digits=3):
        if x is None:
            return "—"
        if x != 0 and (abs(x) < 1e-3 or abs(x) >= 1e5):
            return f"{x:.2e}"
        return f"{x:.{digits}g}"

    if assoc_task.value.startswith("Test"):
        _vcfs = [_p for _p in _v["vcfs"] if chosen(_p)]
        mo.stop(not _vcfs, warn("Select the genotype VCF (or one VCF per person)."))
        mo.stop(not chosen(_v["phenotypes"]), warn("Select the phenotype table."))
        mo.stop(not (_v["outcome"] or "").strip(), warn("Type the name of the outcome column."))
        _kw = {}
        if words(_v["covariates"]):
            _kw["covariates"] = words(_v["covariates"])
        if _v["outcome_type"] != "auto":
            _kw["outcome_type"] = _v["outcome_type"]
        if (_v["case_value"] or "").strip():
            _kw["case_value"] = _v["case_value"].strip()
        if _v["model"] != "additive":
            _kw["model"] = _v["model"]
        if words(_v["variants"]):
            _kw["ids"] = words(_v["variants"])
        if words(_v["regions"]):
            _kw["regions"] = words(_v["regions"])
        if abs(float(_v["min_maf"] or 0) - 0.01) > 1e-12:
            _kw["min_maf"] = float(_v["min_maf"] or 0)
        if abs(float(_v["min_call_rate"] or 0) - 0.9) > 1e-12:
            _kw["min_call_rate"] = float(_v["min_call_rate"] or 0)
        if _v["firth"] != "auto":
            _kw["firth"] = _v["firth"]
        if int(_v["permutations"] or 0):
            _kw["permutations"] = int(_v["permutations"])
        if _v["absent_is_ref"]:
            _kw["absent_is_ref"] = True
        _kw["output"] = output_name(_v["output"], _v["phenotypes"], "association", ".tsv")
        _vcf_arg = _vcfs[0] if len(_vcfs) == 1 else _vcfs
        try:
            with mo.status.spinner("Testing the variants..." + (" (with permutations: this takes longer)" if _kw.get("permutations") else "")):
                _r = await asyncio.to_thread(mmg.genotype_association, _vcf_arg, _v["phenotypes"], _v["outcome"].strip(), **_kw)
            _pe, _va = _r["people"], _r["variants"]
            _binary = _r["outcome_type"] == "binary"
            _dropped = []
            for _key, _what in (("not_in_table", "genotyped but not in the table"), ("not_in_vcf", "in the table but not genotyped"),
                                ("missing_outcome", "without the outcome")):
                if _pe.get(_key):
                    _dropped.append(f"{_pe[_key]} {_what}")
            for _col, _n in (_pe.get("missing_covariates") or {}).items():
                if _n:
                    _dropped.append(f"{_n} without `{_col}`")
            _rows_f = [
                ("People analysed", f"`{_pe['analysed']}`" + (f" (left out: {'; '.join(_dropped)})" if _dropped else "")),
                ("Outcome", f"`{_r['outcome']}`: yes/no — `{_r.get('event_value')}` is the event; {_pe.get('cases')} with it, {_pe.get('controls')} without"
                 if _binary else f"`{_r['outcome']}`: a number (linear regression)"),
                ("Adjusted for", ", ".join(f"`{_c}`" for _c in _r.get("covariates") or []) or "nothing"),
                ("Alleles tested", f"`{_va['tested']}`" + (" — methods: " + ", ".join(f"{_m} {_n}" for _m, _n in (_va.get("methods") or {}).items()) if _va.get("methods") else "")),
                ("With FDR < 0.05", f"`{_va.get('fdr_below_0_05', 0)}`"),
            ]
            if _r.get("permutations"):
                _rows_f.append(("Permutations", f"`{_r['permutations']:,}` (seed {_r.get('seed')}): {_r.get('permutation_scheme', '')}"))
            _skipped = {_k: _n for _k, _n in (_va.get("skipped") or {}).items() if _n}
            if _skipped:
                _rows_f.append(("Not tested", "; ".join(f"{_n} {_k.replace('_', ' ')}" for _k, _n in _skipped.items())))
            if _va.get("not_found"):
                _rows_f.append(("Not found in the VCFs", ", ".join(f"`{_x}`" for _x in _va["not_found"][:20])))
            _pw = _r.get("power") or {}
            if _pw.get("detectable_effect"):
                _scale = "odds ratio" if _binary else "SD of the outcome per copy"
                _rows_f.append(("What this study can detect",
                                f"for an allele at frequency `{_num(_pw.get('maf'))}`, an effect of at least `{_num(_pw['detectable_effect'])}` "
                                f"({_scale}) for {100 * (_pw.get('target_power') or 0.8):.0f}% power at α = `{_num(_pw.get('alpha_per_test'))}` per test"))
            _items = [mo.callout(mo.md("### 📐 Association\n" + facts(_rows_f)), kind="success" if _va["tested"] else "warn")]
            if _r.get("warnings"):
                _items.append(mo.callout(mo.md("**Notes and warnings:**\n\n" + "\n".join(f"* {_w}" for _w in _r["warnings"])), kind="warn"))
            _rows = sorted(_r.get("rows") or [], key=lambda _x: (_x.get("p") is None, _x.get("p") or 1.0))
            _perm = bool(_kw.get("permutations"))
            _table = []
            for _x in _rows[:500]:
                _row = {
                    "variant": _x["id"], "position": f"{_x['chrom']}:{_x['pos']} {_x['ref']}>{_x['alt']}",
                    "people": _x.get("n"), "ALT freq.": _num(_x.get("alt_freq")),
                }
                if _binary:
                    _row["odds ratio [95% CI]"] = f"{_num(_x.get('or'))} [{_num(_x.get('or_low'))}, {_num(_x.get('or_high'))}]"
                else:
                    _row["effect per copy [95% CI]"] = f"{_num(_x.get('beta'))} [{_num(_x.get('ci_low'))}, {_num(_x.get('ci_high'))}]"
                _row["p"] = _num(_x.get("p"))
                _row["FDR"] = _num(_x.get("fdr"))
                if _perm:
                    _row["p (permutation)"] = _num(_x.get("p_permutation"))
                    _row["p (family-wise)"] = _num(_x.get("p_fwer"))
                _row["method"] = _x.get("method", "")
                _row["note"] = _x.get("note") or ""
                _table.append(_row)
            if _table:
                _items.append(mo.ui.table(_table, selection=None, pagination=True, page_size=15,
                                          label="Each allele tested, from the smallest p" + (" (the first 500)" if len(_rows) > 500 else "")))
                _lam = figs.genomic_inflation([_x.get("p") for _x in _rows])
                _qq = figure_panel(figs.qq_plot([_x.get("p") for _x in _rows]), "association_qq",
                                   "Points along the dashed line: the test behaves as it should under no effect. The whole cloud "
                                   "lifted above it (λGC well above 1) means inflation — population structure, relatives, or a "
                                   "test that does not fit the data; with few tests, λGC is only a rough guide.")
                _top = _rows[:15]
                if _binary:
                    _fr = figs.forest([(_x["id"], _x.get("or"), _x.get("or_low"), _x.get("or_high"), _x.get("method")) for _x in _top],
                                      title="Odds ratio of the alleles with the smallest p (95% CI)", null=1.0, log=True,
                                      x_label="odds ratio per copy")
                else:
                    _fr = figs.forest([(_x["id"], _x.get("beta"), _x.get("ci_low"), _x.get("ci_high"), _x.get("method")) for _x in _top],
                                      title="Effect of the alleles with the smallest p (95% CI)", x_label="effect per copy (units of the outcome)")
                _fp = figure_panel(_fr, "association_forest", "An interval that crosses the dashed line is compatible with no effect; "
                                   "a wide one means few carriers.")
                _items.append(mo.hstack([_x for _x in (_qq, _fp) if _x is not None], widths="equal", wrap=True))
                if len({str(_x["chrom"]) for _x in _rows}) > 1 or len(_rows) > 50:
                    _items.append(figure_panel(figs.manhattan(_rows), "association_manhattan",
                                               "Each allele at its chromosome and position; orange: FDR < 0.05; dashed: the Bonferroni threshold."))
                if _lam is not None and _lam > 1.1 and len(_rows) >= 100:
                    _items.append(warn(f"λGC = {_lam:.2f}: the p-values are inflated across the board. Adjust for ancestry "
                                       "(screen 17 gives the proportions) and check for relatives (screen 16) before trusting the hits."))
            _items.append(mo.accordion({"📄 The engine's full report": mo.plain_text(_r.get("text", ""))}))
            if _r.get("output"):
                _items.append(mo.md(f"The table was written to `{_r['output']}`."))
            _items.append(code_panel(f"r = {python_call('genotype_association', _vcf_arg, _v['phenotypes'], _v['outcome'].strip(), **_kw)}",
                                     "print(r['text'])"))
            _panel = mo.vstack(_items)
        except Exception as _err:
            _panel = error_panel(_err)
    else:
        try:
            _alpha = float((_v["alpha"] or "").strip())
        except ValueError:
            _alpha = -1.0
        mo.stop(not 0 < _alpha < 1, warn("Type the significance level as a number between 0 and 1, e.g. 0.05 or 5e-8."))
        _kw = {"n": int(_v["n"] or 0)}
        if float(_v["effect"] or 0) > 0:
            _kw["effect"] = float(_v["effect"])
        if _v["outcome"] != "binary":
            _kw["outcome"] = _v["outcome"]
        if _v["outcome"] == "binary" and abs(float(_v["case_fraction"] or 0) - 0.5) > 1e-12:
            _kw["case_fraction"] = float(_v["case_fraction"])
        if abs(_alpha - 0.05) > 1e-15:
            _kw["alpha"] = _alpha
        if int(_v["tests"] or 1) != 1:
            _kw["tests"] = int(_v["tests"])
        if _v["model"] != "additive":
            _kw["model"] = _v["model"]
        if int(_v["covariates"] or 0):
            _kw["covariates"] = int(_v["covariates"])
        if abs(float(_v["target_power"] or 0) - 0.8) > 1e-12:
            _kw["target_power"] = float(_v["target_power"])
        if int(_v["simulations"] or 0):
            _kw["simulations"] = int(_v["simulations"])
        _af = float(_v["allele_freq"] or 0)
        try:
            with mo.status.spinner("Computing the power..." + (" (simulating studies)" if _kw.get("simulations") else "")):
                _r = await asyncio.to_thread(mmg.association_power, _af, **_kw)
            _binary = _r["outcome_type"] == "binary"
            _scale = "odds ratio" if _binary else "SD per copy"
            _rows_f = []
            if _r.get("power") is not None:
                _sim = ""
                if _r.get("simulated_power") is not None:
                    _sim = f" (simulated: {100 * _r['simulated_power']:.1f}% ± {100 * (_r.get('simulated_se') or 0):.1f})"
                _rows_f.append(("Power", f"`{100 * _r['power']:.1f}%`" + _sim))
            if _r.get("n_needed"):
                _rows_f.append((f"People needed for {100 * _r['target_power']:.0f}% power", f"`{_r['n_needed']:,}`"))
            if _r.get("detectable_effect"):
                _rows_f.append((f"Smallest effect detectable with {_r['n']:,} people", f"`{_num(_r['detectable_effect'])}` ({_scale})"))
            elif _r.get("n"):
                _rows_f.append((f"Smallest effect detectable with {_r['n']:,} people", "none in the range searched"))
            _rows_f.append(("α per test", f"`{_num(_r.get('alpha_per_test'))}`"))
            if _binary and _r.get("expected_events") is not None:
                _rows_f.append(("Expected", f"{_num(_r['expected_events'])} people with the event, {_num(_r.get('expected_carriers'))} carriers"))
            _ok = (_r.get("power") or 0) >= (_r.get("target_power") or 0.8)
            _items = [mo.callout(mo.md("### 🧮 Study power\n" + facts(_rows_f)), kind="success" if _ok else "warn")]
            if _r.get("warnings"):
                _items.append(mo.callout(mo.md("\n".join(f"* {_w}" for _w in _r["warnings"])), kind="warn"))
            _curve = _r.get("curve") or {}
            _cn, _cp = _curve.get("n") or [], _curve.get("power") or []
            if len(_cn) >= 2:
                _target = _r.get("target_power") or 0.8
                _pts = [(float(_a), float(_b), "up" if _b >= _target else "down", f"{_a:,} people: {100 * _b:.1f}%")
                        for _a, _b in zip(_cn, _cp)]
                _chart = scatter_chart(_pts, "people in the study", "power", "Power against the size of the study",
                                       hlines=(_target,), vlines=(_r["n"],) if _r.get("n") else (),
                                       names={"up": "at or above the target power", "down": "below it"})
                if _chart is not None:
                    _items.append(_chart)
                    _items.append(mo.md("*The dashed lines mark the target power and the size typed.*"))
            _items.append(mo.accordion({"📄 The engine's full report": mo.plain_text(_r.get("text", ""))}))
            _items.append(code_panel(f"r = {python_call('association_power', _af, **_kw)}", "print(r['text'])"))
            _panel = mo.vstack(_items)
        except Exception as _err:
            _panel = error_panel(_err)
    _panel
    return


@app.cell
def rnaseq_page_ui(get_page, mo, set_page):
    mo.stop(get_page() != "rnaseq")
    btn_back_rnaseq = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_rnaseq = mo.ui.button(label="🔄 Refresh the file list")
    rnaseq_task = mo.ui.radio(options=["Quantify from FASTQ (pseudoalignment)", "Build a transcript index", "See how one read is pseudoaligned",
                                       "Count reads per gene", "Differential expression", "Test gene sets"], value="Count reads per gene",
                              inline=True, label="**What to do:**")
    rnaseq_folder = mo.ui.text(value=".", label="Data folder (its subfolders are searched too):", full_width=True)
    _ui = mo.vstack([
        btn_back_rnaseq,
        mo.md("<br>\n\n# 📊 22. RNA-seq: Gene Counts & Differential Expression"),
        mo.callout(mo.md(
            "Which genes a treatment switches on or off. **From the reads (FASTQ)**, with no aligner: each read is "
            "**pseudoaligned** against an index of the transcripts (the transcripts that hold all its k-mers) and the abundance of "
            "every transcript is estimated by the EM of **kallisto**; the gene counts go straight to differential expression. "
            "**Counting**: each RNA-seq read aligned to the genome (BAM) is assigned "
            "to the gene whose exons it falls in, as featureCounts and htseq-count do. **Differential expression**: the counts of "
            "two conditions are compared gene by gene in the manner of **DESeq2**: samples normalized by their sequencing depth, "
            "the gene's variability estimated and shrunk toward the trend of genes like it, a negative-binomial model and a Wald "
            "test, with the false-discovery rate (padj) for thousands of genes tested at once. Add a blocking factor (for "
            "example the donor) to compare the conditions within each one. **Gene sets**: does a group of genes (a process, a "
            "cell type) move together more than the rest? Each set of a GMT file, or the built-in **spinal cord injury** sets "
            "(glial scar, inflammation, axon regeneration, myelin, laminin receptors...), is tested against the result table, "
            "as limma's wilcoxGST does. Experimental, for teaching and research."), kind="info"),
        rnaseq_folder,
        mo.hstack([rnaseq_task, btn_refresh_rnaseq], justify="space-between", align="center"),
    ])
    _ui
    return btn_back_rnaseq, btn_refresh_rnaseq, rnaseq_folder, rnaseq_task


@app.cell
def rnaseq_page_example(example_picker, get_page, mo, rnaseq_task):
    mo.stop(get_page() != "rnaseq")
    rnaseq_example = example_picker("rnaseq", rnaseq_task.value)
    rnaseq_example
    return (rnaseq_example,)


@app.cell
def rnaseq_page_form(NO_FILE, btn_refresh_rnaseq, example_for, files_under, get_page, label_of, make_form, mo, rnaseq_example, rnaseq_folder, rnaseq_task):
    mo.stop(get_page() != "rnaseq")
    btn_refresh_rnaseq  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("rnaseq", rnaseq_task.value, rnaseq_example.value)
    _folder = rnaseq_folder.value
    if rnaseq_task.value == "Build a transcript index":
        _maps = {"(from the FASTA headers)": "none", **{_f: _f for _f in files_under(_folder, "gff", "table") if _f != NO_FILE}}
        form_rnaseq = make_form({
            "fasta": mo.ui.dropdown(options=files_under(_folder, "fasta"), value=_ex.get("fasta"),
                                    label="Transcript sequences (cDNA FASTA: Ensembl, GENCODE or RefSeq):", full_width=True),
            "gene_map": mo.ui.dropdown(options=_maps, value=label_of(_maps, _ex.get("gene_map", "none"), "(from the FASTA headers)"),
                                       label="Gene of each transcript (a GTF, GFF3 or transcript/gene table; or the headers):", full_width=True),
            "k": mo.ui.number(start=7, stop=31, step=2, value=_ex.get("k", 31), label="k-mer length (odd; 31 for reads of 50 bases or more):"),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: <fasta>.idx next to it", label="Index file:", full_width=True),
        }, "📇 Build the index")
    elif rnaseq_task.value == "Quantify from FASTQ (pseudoalignment)":
        _conds = {"(none)": "none", **{_f: _f for _f in files_under(_folder, "table") if _f != NO_FILE}}
        form_rnaseq = make_form({
            "index": mo.ui.dropdown(options=files_under(_folder, "tx_index"), value=_ex.get("index"),
                                    label="Transcript index (.idx: the task Build a transcript index makes it):", full_width=True),
            "fastqs": mo.ui.multiselect(options=files_under(_folder, "fastq"), value=_ex.get("fastqs", []),
                                        label="Reads (FASTQ; x_R1 and x_R2 are paired automatically):", full_width=True),
            "conditions": mo.ui.dropdown(options=_conds, value=label_of(_conds, _ex.get("conditions", "none"), "(none)"),
                                         label="Sample table (sample, condition), for differential expression next:", full_width=True),
            "fragment_length": mo.ui.number(start=0, stop=1000, step=1, value=_ex.get("fragment_length", 0),
                                            label="Fragment length, mean (single-end reads need it; for pairs, 0 measures it):"),
            "fragment_sd": mo.ui.number(start=0, stop=500, step=1, value=_ex.get("fragment_sd", 0), label="Fragment length, standard deviation:"),
            "max_reads": mo.ui.number(start=0, stop=10**10, step=1000, value=_ex.get("max_reads", 0),
                                      label="Only the first reads of each sample (0: all; for a quick look):"),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: a folder quant next to the index", label="Output folder:", full_width=True),
        }, "📊 Quantify")
    elif rnaseq_task.value == "See how one read is pseudoaligned":
        form_rnaseq = make_form({
            "index": mo.ui.dropdown(options=files_under(_folder, "tx_index"), value=_ex.get("index"), label="Transcript index (.idx):", full_width=True),
            "read": mo.ui.text_area(value=_ex.get("read", ""), placeholder="ACGT...", label="A read (its sequence):", full_width=True),
            "read2": mo.ui.text_area(value=_ex.get("read2", ""), placeholder="optional", label="Its mate, for a pair:", full_width=True),
        }, "🔍 Show")
    elif rnaseq_task.value == "Count reads per gene":
        _attrs = {"default: gene_id, else the parent gene": "default", "gene_id": "gene_id", "gene": "gene", "Name": "Name",
                  "GeneID (RefSeq Dbxref)": "GeneID"}
        _strands = {"unstranded": "unstranded", "forward (read 1 on the gene's strand)": "forward", "reverse (dUTP, most kits)": "reverse"}
        form_rnaseq = make_form({
            "bams": mo.ui.multiselect(options=files_under(_folder, "sam"), value=_ex.get("bams", []), label="Aligned reads (BAM):", full_width=True),
            "annotation": mo.ui.dropdown(options=files_under(_folder, "gff"), value=_ex.get("annotation"), label="Annotation (GFF3 / GTF):", full_width=True),
            "feature": mo.ui.text(value=_ex.get("feature", "exon"), label="Feature counted:"),
            "attribute": mo.ui.dropdown(options=_attrs, value=label_of(_attrs, _ex.get("attribute", "default"), "default: gene_id, else the parent gene"), label="Genes from the attribute:", full_width=True),
            "strandedness": mo.ui.dropdown(options=_strands, value=label_of(_strands, _ex.get("strandedness", "unstranded"), "unstranded"), label="Library strandedness:", full_width=True),
            "min_mapq": mo.ui.number(start=0, stop=255, step=1, value=_ex.get("min_mapq", 10), label="Minimum mapping quality:"),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: counts.tsv next to the first BAM", label="Count table:", full_width=True),
        }, "📊 Count")
    elif rnaseq_task.value == "Test gene sets":
        _tables = [_f for _f in files_under(_folder, "table") if _f != NO_FILE] or [NO_FILE]
        _builtin = "Spinal cord injury: 11 teaching sets (built in)"
        _sets = {_builtin: "spinal_cord_injury", **{_f: _f for _f in files_under(_folder, "geneset", "table") if _f != NO_FILE}}
        _stats = {"stat (the Wald statistic)": "stat", "log2FoldChange": "log2FoldChange"}
        form_rnaseq = make_form({
            "results": mo.ui.dropdown(options=_tables, value=_ex.get("results"),
                                      label="Differential-expression table (gene, name, ..., stat, padj: this screen's or DESeq2's):", full_width=True),
            "sets": mo.ui.dropdown(options=_sets, value=label_of(_sets, _ex.get("sets", "spinal_cord_injury"), _builtin),
                                   label="Gene sets (a GMT file, a set/gene table, or built in):", full_width=True),
            "statistic": mo.ui.dropdown(options=_stats, value=label_of(_stats, _ex.get("statistic", "stat"), "stat (the Wald statistic)"),
                                        label="Rank the genes by:", full_width=True),
            "require_padj": mo.ui.checkbox(value=_ex.get("require_padj", True),
                                           label="Rank only the genes with an adjusted p-value (those that passed the expression filter; recommended)"),
            "alpha": mo.ui.number(start=0.001, stop=0.5, step=0.01, value=_ex.get("alpha", 0.1), label="False-discovery rate threshold:"),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: <results>.sets.tsv next to it", label="Result table:", full_width=True),
        }, "📊 Test the sets")
    else:
        _tables = [_f for _f in files_under(_folder, "table") if _f != NO_FILE] or [NO_FILE]
        _names = {"(none: gene IDs only)": "none", **{_f: _f for _f in _tables if _f != NO_FILE}}
        form_rnaseq = make_form({
            "counts": mo.ui.dropdown(options=_tables, value=_ex.get("counts"), label="Count table (gene + one column per sample):", full_width=True),
            "samples": mo.ui.dropdown(options=_tables, value=_ex.get("samples"), label="Sample table (sample + one column per factor):", full_width=True),
            "condition": mo.ui.text(value=_ex.get("condition", "condition"), label="Condition column:"),
            "reference": mo.ui.text(value=_ex.get("reference", ""), placeholder="control", label="Reference level (the baseline):"),
            "treated": mo.ui.text(value=_ex.get("treated", ""), placeholder="blank: the other level", label="Level compared with it:"),
            "factors": mo.ui.text(value=_ex.get("factors", ""), placeholder="cell_line", label="Blocking factors (columns), if any:", full_width=True),
            "alpha": mo.ui.number(start=0.001, stop=0.5, step=0.01, value=_ex.get("alpha", 0.1), label="False-discovery rate (padj) threshold:"),
            "names": mo.ui.dropdown(options=_names, value=label_of(_names, _ex.get("names", "none"), "(none: gene IDs only)"), label="Gene names (gene_info or a two-column table):", full_width=True),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: <counts>.de.tsv next to it", label="Result table:", full_width=True),
        }, "📊 Compare")
    mo.vstack([_x for _x in (_about, form_rnaseq) if _x is not None])
    return (form_rnaseq,)


@app.cell
async def rnaseq_page_backend(asyncio, beside, chosen, code_panel, error_panel, facts, figs, figure_panel, form_rnaseq, get_page, mmg, mo, np, python_call, rnaseq_task, scatter_chart, split_name, warn, words):
    mo.stop(get_page() != "rnaseq")
    # the tasks that start from the reads are run by rnaseq_fastq_backend
    mo.stop(rnaseq_task.value in ("Build a transcript index", "Quantify from FASTQ (pseudoalignment)", "See how one read is pseudoaligned"))
    mo.stop(form_rnaseq.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    import math as _m
    _v = form_rnaseq.value
    try:
        if rnaseq_task.value == "Count reads per gene":
            _bams = [_p for _p in _v["bams"] if chosen(_p)]
            mo.stop(not _bams, warn("Select the BAM files, one per sample."))
            mo.stop(not chosen(_v["annotation"]), warn("Select the annotation (GFF3 or GTF) of the same assembly."))
            _out = (_v["output"] or "").strip() or beside(_bams[0], "counts.tsv")
            _kw = {}
            if (_v["feature"] or "exon").strip() != "exon":
                _kw["feature_type"] = _v["feature"].strip()
            if _v["attribute"] != "default":
                _kw["attribute"] = _v["attribute"]
            if _v["strandedness"] != "unstranded":
                _kw["strandedness"] = _v["strandedness"]
            if int(_v["min_mapq"]) != 10:
                _kw["min_mapq"] = int(_v["min_mapq"])
            with mo.status.spinner(f"Counting the reads of {len(_bams)} sample(s)..."):
                _r = await asyncio.to_thread(mmg.count_reads, _bams, _v["annotation"], _out, **_kw)
            _items = [mo.callout(mo.md("### 📊 Reads counted\n" + facts([
                ("Genes", f"`{len(_r['genes']):,}` from `{_r['features']:,}` {_r['feature_type']} features"),
                ("Count table", f"`{_out}`"),
            ])), kind="success")]
            _items.append(mo.ui.table([{
                "sample": _s["sample"], "fragments": _s["fragments"], "assigned": _s["assigned"],
                "assigned %": f"{100 * _s['assigned'] / _s['fragments']:.1f}%" if _s["fragments"] else "",
                "no feature": _s["no_feature"], "ambiguous": _s["ambiguous"], "multimapping": _s["multimapping"], "low quality": _s["low_mapq"],
            } for _s in _r["samples"]], selection=None, label="Assignment per sample"))
            for _w in _r["warnings"][:5]:
                _items.append(warn(_w))
            _items.append(code_panel(f"r = {python_call('count_reads', _bams, _v['annotation'], _out, **_kw)}", 'print(r["text"])'))
        elif rnaseq_task.value == "Test gene sets":
            mo.stop(not chosen(_v["results"]), warn("Select a differential-expression table (the one the task above writes, or DESeq2's)."))
            _out = (_v["output"] or "").strip() or beside(_v["results"], split_name(_v["results"])[0] + ".sets.tsv")
            _kw = {"output": _out}
            if _v["sets"] != "spinal_cord_injury":
                _kw["gene_sets"] = _v["sets"]
            if _v["statistic"] != "stat":
                _kw["statistic"] = _v["statistic"]
            if not _v["require_padj"]:
                _kw["require_padj"] = False
            if float(_v["alpha"]) != 0.1:
                _kw["alpha"] = float(_v["alpha"])
            with mo.status.spinner("Ranking the genes and testing each set..."):
                _r = await asyncio.to_thread(mmg.gene_set_test, _v["results"], **_kw)
            _alpha = _r["alpha"]
            _sets = sorted(_r["sets"], key=lambda s: (s["pvalue"] is None, s["pvalue"] if s["pvalue"] is not None else 1.0))
            _dir = lambda s: ("–" if s["fdr"] is None or s["fdr"] >= _alpha else ("⬆️ up" if s["z"] > 0 else "⬇️ down"))
            _items = [mo.callout(mo.md("### 📊 Gene sets\n" + facts([
                ("Sets", f"`{len(_r['sets'])}` ({_r['gene_sets']})"),
                ("Genes ranked", f"`{_r['genes_ranked']:,}`, by `{_r['statistic']}`"),
                (f"Moving one way at FDR < {_alpha:g}", f"`{_r['significant']}`"),
                (f"Changing either way at FDR < {_alpha:g}", f"`{_r['significant_mixed']}`"),
                ("Result table", f"`{_out}`"),
            ])), kind="success")]
            _items.append(mo.ui.table([{
                "set": _s["set"], "genes found": f"{_s['found']} of {_s['size']}", "direction": _dir(_s),
                "z": "" if _s["z"] is None else round(_s["z"], 2),
                "FDR": "too few genes" if _s["fdr"] is None else f"{_s['fdr']:.3g}",
                "FDR, either way": "" if _s["fdr_mixed"] is None else f"{_s['fdr_mixed']:.3g}",
                "mean log2 fold change": "" if _s["mean_log2_fold_change"] is None else round(_s["mean_log2_fold_change"], 2),
                "genes up / down": f"{_s['up']} / {_s['down']}", "what it is": _s["description"],
            } for _s in _sets], selection=None, pagination=True, page_size=12, label="Each set against the other genes, most significant first"))
            _pts = [(_s["z"], -_m.log10(max(_s["fdr"], 1e-300)),
                     "ns" if _s["fdr"] >= _alpha else ("up" if _s["z"] > 0 else "down"), f"{_s['set']}: z {_s['z']:.2f}, FDR {_s['fdr']:.2g}")
                    for _s in _sets if _s["z"] is not None and _s["fdr"] is not None]
            _chart = scatter_chart(_pts, "z (rank of the set's genes against the rest)", "−log10 FDR", "Gene sets",
                                   labels=[(_p[0], _p[1], _p[3].split(":")[0]) for _p in _pts[:10]],
                                   hlines=[-_m.log10(_alpha)], vlines=[0])
            if _chart is not None:
                _items.append(_chart)
            _items.append(mo.accordion({f"{_s['set']} ({_s['found']} of {_s['size']} genes)": mo.md(
                "Found: " + (", ".join(f"`{_g}`" for _g in _s["genes"]) or "none")
                + ("\n\nNot in the table, or not tested: " + ", ".join(f"`{_g}`" for _g in _s["missing"]) if _s["missing"] else ""))
                for _s in _sets}))
            _items.append(mo.callout(mo.md(
                "**Reading it.** *z* above 0: the set's genes rank higher (more up) than the other genes; below 0, lower. The test is "
                "**competitive** and, like limma's wilcoxGST, treats genes as independent: genes of one process often move together, "
                "which makes these p-values smaller than they should be. Read them as a guide to where to look, then check the genes "
                "themselves. The built-in sets are marker genes cited in reviews, for teaching and exploratory work: not validated "
                "signatures."), kind="info"))
            _items.append(code_panel(f"r = {python_call('gene_set_test', _v['results'], **_kw)}", 'print(r["text"])'))
        else:
            mo.stop(not chosen(_v["counts"]) or not chosen(_v["samples"]), warn("Select the count table and the sample table."))
            mo.stop(not (_v["reference"] or "").strip(), warn("Type the reference level of the condition (e.g. control, untrt)."))
            _out = (_v["output"] or "").strip() or beside(_v["counts"], split_name(_v["counts"])[0] + ".de.tsv")
            _kw = {"output": _out}
            if (_v["treated"] or "").strip():
                _kw["treated"] = _v["treated"].strip()
            if words(_v["factors"]):
                _kw["factors"] = words(_v["factors"])
            if float(_v["alpha"]) != 0.1:
                _kw["alpha"] = float(_v["alpha"])
            if _v["names"] != "none":
                _kw["gene_names"] = _v["names"]
            _cond, _ref = _v["condition"].strip(), _v["reference"].strip()
            with mo.status.spinner("Fitting the model for every gene..."):
                _r = await asyncio.to_thread(mmg.differential_expression, _v["counts"], _v["samples"], _cond, _ref, **_kw)
            _p = _r["plot"]
            _alpha = _r["alpha"]
            _n = len(_p["gene"])
            _label = lambda i: (_p["name"][i] if _p.get("name") and _p["name"][i] else _p["gene"][i])
            _kind = lambda i: ("ns" if _p["padj"][i] is None or _p["padj"][i] >= _alpha else ("up" if _p["log2_fold_change"][i] > 0 else "down"))
            _sig = sorted([_i for _i in range(_n) if _p["padj"][_i] is not None and _p["padj"][_i] < _alpha], key=lambda _i: _p["padj"][_i])
            _vol = [(_p["log2_fold_change"][_i], -_m.log10(max(_p["padj"][_i], 1e-300)), _kind(_i),
                     f"{_label(_i)}: log2FC {_p['log2_fold_change'][_i]:.2f}, padj {_p['padj'][_i]:.2g}")
                    for _i in range(_n) if _p["padj"][_i] is not None and _p["log2_fold_change"][_i] is not None]
            _ma = [(_m.log10(_p["base_mean"][_i]), _p["log2_fold_change"][_i], _kind(_i), f"{_label(_i)}: mean {_p['base_mean'][_i]:.0f}")
                   for _i in range(_n) if _p["base_mean"][_i] and _p["log2_fold_change"][_i] is not None]
            _top = [(_p["log2_fold_change"][_i], -_m.log10(max(_p["padj"][_i], 1e-300)), _label(_i)) for _i in _sig[:10]]
            _items = [mo.callout(mo.md("### 📊 Differential expression\n" + facts([
                ("Comparison", f"`{_r['treated']}` against `{_r['reference']}` in `{_r['condition']}`; design `{_r['design']}`"),
                ("Genes tested", f"`{_r['tested']:,}`"),
                (f"Changed at padj < {_alpha:g}", f"`{_r['up']:,}` up, `{_r['down']:,}` down"),
                ("Result table", f"`{_out}`"),
            ])), kind="success")]
            _items.append(mo.ui.table([{
                "gene": _label(_i), "base mean": round(_p["base_mean"][_i], 1), "log2 fold change": round(_p["log2_fold_change"][_i], 3),
                "padj": f"{_p['padj'][_i]:.3g}",
            } for _i in _sig[:500]], selection=None, pagination=True, page_size=10, label=f"The changed genes, most significant first ({len(_sig):,})"))
            _charts = [_c for _c in (
                scatter_chart(_vol, "log2 fold change", "−log10 padj", "Volcano plot", labels=_top, hlines=[-_m.log10(_alpha)], vlines=[0]),
                scatter_chart(_ma, "log10 mean of normalized counts", "log2 fold change", "MA plot", hlines=[0]),
            ) if _c is not None]
            if _charts:
                _items.append(mo.hstack(_charts, widths="equal", wrap=True))
            # Diagnostics: are the samples what the design says (PCA, distances), and is the test calibrated (p-values)?
            _diag = []
            try:
                import gzip as _gz
                _opener = _gz.open if str(_v["counts"]).endswith(".gz") else open
                with _opener(_v["counts"], "rt", encoding="utf-8") as _fh:
                    _head = _fh.readline().rstrip("\n").split("\t")
                    _mat = [_line.rstrip("\n").split("\t") for _line in _fh if _line.strip()]
                _cols = [_head.index(_s0) for _s0 in _r["samples"]]
                _X = np.array([[float(_row[_c]) for _c in _cols] for _row in _mat], dtype=float)
                _X = np.log2(_X / np.array(_r["size_factors"])[None, :] + 1.0)
                _keep = _X.mean(axis=1) > 0
                _X = _X[_keep]
                _var = _X.var(axis=1)
                _sel = _X[np.argsort(-_var)[:500]]
                _C = _sel - _sel.mean(axis=1, keepdims=True)
                _U, _S, _Vt = np.linalg.svd(_C, full_matrices=False)
                _pv = _S ** 2 / max(float((_S ** 2).sum()), 1e-12)
                _sampletab = {}
                with open(_v["samples"], encoding="utf-8") as _fh:
                    _sh = _fh.readline().rstrip("\n").split("\t")
                    _ci = _sh.index(_cond) if _cond in _sh else 1
                    for _line in _fh:
                        _f = _line.rstrip("\n").split("\t")
                        if len(_f) > _ci:
                            _sampletab[_f[0]] = _f[_ci]
                _pts = [(float(_Vt[0][_j] * _S[0]), float(_Vt[1][_j] * _S[1]) if len(_S) > 1 else 0.0, _sampletab.get(_s0, "—"), _s0)
                        for _j, _s0 in enumerate(_r["samples"])]
                _diag.append(figure_panel(
                    figs.group_scatter(_pts, f"PC1 ({100 * _pv[0]:.0f}% of the variance)",
                                       f"PC2 ({100 * _pv[1]:.0f}%)" if len(_pv) > 1 else "PC2", "Samples: PCA of the 500 most variable genes",
                                       label_points=len(_pts) <= 16),
                    "rnaseq_sample_pca",
                    "Samples of the same condition should sit together. A sample among the other group may be swapped or "
                    "mislabeled; a split along PC1 that follows a batch, a donor or a cell line, not the condition, is a factor to "
                    "add to the design."))
                _D = np.sqrt(((_sel[:, :, None] - _sel[:, None, :]) ** 2).sum(axis=0))
                _order = sorted(range(len(_r["samples"])), key=lambda _j: (_sampletab.get(_r["samples"][_j], ""), _r["samples"][_j]))
                _labs = [f"{_r['samples'][_j]} ({_sampletab.get(_r['samples'][_j], '—')})" for _j in _order]
                _diag.append(figure_panel(
                    figs.heatmap([[float(_D[_a, _b]) for _b in _order] for _a in _order], _labs, _labs,
                                 title="Distances between samples", value_label="Euclidean distance (log2 normalized counts)", reverse=True),
                    "rnaseq_sample_distances", "Darker is closer. Blocks along the diagonal should match the conditions."))
            except Exception as _derr:
                _diag.append(mo.md(f"*The sample diagnostics could not be drawn: {_derr}*"))
            try:
                with open(_out, encoding="utf-8") as _fh:
                    _th = _fh.readline().rstrip("\n").split("\t")
                    _pi = _th.index("pvalue")
                    _pvals = []
                    for _line in _fh:
                        _f = _line.rstrip("\n").split("\t")
                        if len(_f) > _pi and _f[_pi] not in ("", "NA", "nan", "None"):
                            _pvals.append(float(_f[_pi]))
                _diag.append(figure_panel(
                    figs.histogram(_pvals, bins=20, title="p-values of the genes tested", x_label="p-value", lo=0, hi=1),
                    "rnaseq_pvalues",
                    "Expected: flat, with a peak at the left for the genes that change. A hump in the middle or a peak at the "
                    "right means the test does not fit the data (outliers, a missing factor, too few replicates)."))
            except Exception:
                pass
            if _diag:
                _items.append(mo.hstack(_diag, widths="equal", wrap=True))
            for _w in _r["warnings"][:5]:
                _items.append(warn(_w))
            _args = {_k: _x for _k, _x in _kw.items() if _k != "output"}
            _items.append(code_panel(f"r = {python_call('differential_expression', _v['counts'], _v['samples'], _cond, _ref, output=_out, **_args)}",
                                     'print(r["text"])'))
        _panel = mo.vstack(_items)
    except Exception as _err:
        _panel = error_panel(_err)
    _panel
    return


@app.cell
async def rnaseq_fastq_backend(asyncio, beside, chosen, code_panel, error_panel, facts, form_rnaseq, get_page, mmg, mo, os, pair_fastqs, python_call, rnaseq_task, split_name, warn):
    # The three tasks that start from the reads (FASTQ) rather than from a BAM.
    mo.stop(get_page() != "rnaseq" or rnaseq_task.value not in (
        "Build a transcript index", "Quantify from FASTQ (pseudoalignment)", "See how one read is pseudoaligned"))
    mo.stop(form_rnaseq.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_rnaseq.value

    def _size(n):
        for _unit in ("bytes", "KB", "MB", "GB"):
            if n < 1024 or _unit == "GB":
                return f"{n:,.0f} {_unit}" if _unit == "bytes" else f"{n:,.1f} {_unit}"
            n /= 1024

    try:
        if rnaseq_task.value == "Build a transcript index":
            mo.stop(not chosen(_v["fasta"]), warn("Select the FASTA of the transcript sequences."))
            _out = (_v["output"] or "").strip() or beside(_v["fasta"], split_name(_v["fasta"])[0] + ".idx")
            _kw = {}
            if int(_v["k"] or 31) != 31:
                _kw["k"] = int(_v["k"])
            if _v["gene_map"] != "none":
                _kw["gene_map"] = _v["gene_map"]
            with mo.status.spinner("Building the index: every k-mer of the transcripts, on both strands..."):
                _r = await asyncio.to_thread(mmg.build_transcript_index, _v["fasta"], _out, **_kw)
            _rows = [
                ("Index", f"`{_r['index']}` ({_size(_r['bytes'])}, {_r['passes']} pass(es), {_r['seconds']:.1f} s)"),
                ("Transcripts", f"`{_r['transcripts']:,}` of `{_r['genes']:,}` genes (gene of each: {_r.get('gene_map') or 'the FASTA headers'})"),
                ("k-mers", f"`{_r['kmers']:,}` (k = {_r['k']}), in `{_r['equivalence_classes']:,}` equivalence classes — the sets of transcripts a k-mer can come from"),
            ]
            if _r.get("transcripts_without_gene"):
                _rows.append(("Transcripts without a gene", f"`{_r['transcripts_without_gene']:,}` (each is its own gene)"))
            if _r.get("short_transcripts"):
                _rows.append(("Shorter than k", f"`{_r['short_transcripts']:,}` (no k-mer: they get no reads)"))
            _items = [mo.callout(mo.md("### 📇 Transcript index\n" + facts(_rows)), kind="success")]
            if _r.get("warnings"):
                _items.append(mo.callout(mo.md("\n".join(f"* {_w}" for _w in _r["warnings"])), kind="warn"))
            _items.append(mo.md("*Next: the task **Quantify from FASTQ** with this index (click **Refresh the file list** to see it).*"))
            _items.append(code_panel(f"r = {python_call('build_transcript_index', _v['fasta'], _out, **_kw)}"))
            _panel = mo.vstack(_items)
        elif rnaseq_task.value == "Quantify from FASTQ (pseudoalignment)":
            mo.stop(not chosen(_v["index"]), warn("Select the transcript index (.idx): build it with the task Build a transcript index."))
            _fq = [_p for _p in _v["fastqs"] if chosen(_p)]
            mo.stop(not _fq, warn("Select the FASTQ files of the samples."))
            _groups = pair_fastqs(_fq)
            _single = any(_r2 is None for _, _, _r2 in _groups)
            mo.stop(_single and not float(_v["fragment_length"] or 0),
                    warn("Single-end reads need the fragment length (mean and standard deviation, e.g. 200 and 20): "
                         "it cannot be measured from one read."))
            _cond = {}
            if _v["conditions"] != "none":
                import csv as _csv
                with open(_v["conditions"], encoding="utf-8", newline="") as _fh:
                    _first = _fh.readline()
                _delim = "\t" if "\t" in _first else ","
                with open(_v["conditions"], encoding="utf-8", newline="") as _fh:
                    _table = [_x for _x in _csv.reader(_fh, delimiter=_delim) if _x and not _x[0].startswith("#")]
                if _table:
                    _head = [_h.strip().lower() for _h in _table[0]]
                    _col = _head.index("condition") if "condition" in _head else 1
                    for _x in _table[1:]:
                        if len(_x) > _col:
                            _cond[_x[0].strip()] = _x[_col].strip()
            _samples = []
            for _name, _r1, _r2 in _groups:
                _s = {"name": _name, "r1": _r1}
                if _r2:
                    _s["r2"] = _r2
                if _name in _cond:
                    _s["condition"] = _cond[_name]
                _samples.append(_s)
            _out = (_v["output"] or "").strip() or beside(_v["index"], "quant")
            _kw = {}
            if float(_v["fragment_length"] or 0):
                _kw["fragment_length"] = float(_v["fragment_length"])
            if float(_v["fragment_sd"] or 0):
                _kw["fragment_sd"] = float(_v["fragment_sd"])
            if int(_v["max_reads"] or 0):
                _kw["max_reads"] = int(_v["max_reads"])
            with mo.status.spinner(f"Pseudoaligning the reads of {len(_samples)} sample(s) and estimating the abundances..."):
                _r = await asyncio.to_thread(mmg.quantify_samples, _v["index"], _samples, _out, **_kw)
            _table_s = []
            for _s in _r["samples"]:
                _table_s.append({
                    "sample": _s["sample"], "condition": _s.get("condition") or "", "reads": "pairs" if _s["paired"] else "single",
                    "fragments": f"{_s['reads']:,}", "pseudoaligned": f"{_s['pseudoaligned_pct']:.1f}%", "unique": f"{_s['unique_pct']:.1f}%",
                    "fragment length": f"{_s['fragment_length_mean']:.0f} ± {_s['fragment_length_sd']:.0f} ({_s['fragment_length_source']})",
                    "EM rounds": _s["em_rounds"],
                })
            _files = _r.get("files") or {}
            _rows = [
                ("Samples", f"`{len(_r['samples'])}`, against `{_r['transcripts']:,}` transcripts of `{_r['genes']:,}` genes"),
                ("Gene counts", f"`{_files.get('gene_counts', '')}` (estimated counts, rounded: the input of differential expression)"),
                ("Also written", ", ".join(f"`{os.path.basename(_p)}`" for _k, _p in _files.items() if _k != "gene_counts" and _p) +
                 f" and one folder per sample (`abundance.tsv`, as kallisto's) in `{_out}`"),
            ]
            _items = [mo.callout(mo.md("### 📊 Quantified from the reads\n" + facts(_rows)), kind="success")]
            if _r.get("warnings"):
                _items.append(mo.callout(mo.md("\n".join(f"* {_w}" for _w in _r["warnings"])), kind="warn"))
            _items.append(mo.ui.table(_table_s, selection=None, label="Each sample"))
            _ids, _names, _tpm = _r.get("gene_ids") or [], _r.get("gene_names") or [], _r.get("gene_tpm") or []
            if _ids and _tpm:
                _mean = [sum(_col[_g] for _col in _tpm) / len(_tpm) for _g in range(len(_ids))]
                _top = sorted(range(len(_ids)), key=lambda _g: -_mean[_g])[:15]
                _items.append(mo.ui.table([
                    dict({"gene": _ids[_g], "name": _names[_g] if _g < len(_names) else ""},
                         **{_s["sample"]: round(_tpm[_j][_g], 1) for _j, _s in enumerate(_r["samples"])})
                    for _g in _top], selection=None, label="The 15 most expressed genes (TPM)"))
            if _files.get("samples"):
                _items.append(mo.md(f"*Next: the task **Differential expression**, with `{_files['gene_counts']}` as the count table, "
                                    f"`{_files['samples']}` as the sample table and `{_files.get('gene_names', '')}` for the gene names "
                                    "(click **Refresh the file list** first).*"))
            else:
                _items.append(mo.md("*Next: the task **Differential expression** with the gene counts and a sample table (sample, "
                                    "condition); choose a sample table here and it is written for you.*"))
            _items.append(mo.accordion({"📄 The engine's full report": mo.plain_text(_r.get("text", ""))}))
            _items.append(code_panel(f"samples = {_samples!r}",
                                     f"r = mmg.quantify_samples({_v['index']!r}, samples, {_out!r}" + "".join(f", {_k}={_x!r}" for _k, _x in _kw.items()) + ")"))
            _panel = mo.vstack(_items)
        else:
            mo.stop(not chosen(_v["index"]), warn("Select the transcript index (.idx)."))
            _read = "".join((_v["read"] or "").split()).upper()
            mo.stop(not _read, warn("Type (or paste) the sequence of a read."))
            _read2 = "".join((_v["read2"] or "").split()).upper() or None
            _kw = {"read2": _read2} if _read2 else {}
            _r = await asyncio.to_thread(mmg.explain_pseudoalignment, _v["index"], _read, **_kw)
            _meaning = {"unique": "compatible with one transcript", "multi": "compatible with several transcripts (the EM shares it)",
                        "conflicting": "its k-mers point to transcripts with nothing in common: not counted",
                        "no_kmer": "none of its k-mers is in the index: not counted", "too_short": "shorter than k: not counted"}
            _rows = [
                ("Result", f"**{_r['status']}** — {_meaning.get(_r['status'], '')}"),
                ("Transcripts", ", ".join(f"`{_t}`" for _t in _r["transcripts"]) or "none"),
                ("Genes", ", ".join(f"`{_g}`" for _g in _r["genes"]) or "none"),
            ]
            _items = [mo.callout(mo.md(f"### 🔍 One read, k = {_r['k']}\n" + facts(_rows)), kind="success" if _r["transcripts"] else "warn")]
            for _i, _m in enumerate(_r["mates"]):
                _items.append(mo.ui.table([{
                    "position": _x["position"], "k-mer": _x["kmer"], "in the index": "yes" if _x["found"] else "no",
                    "class": "" if _x["class"] is None else _x["class"], "transcripts": ", ".join(_x["transcripts"]),
                } for _x in _m["kmers"]], selection=None, pagination=True, page_size=12,
                    label=f"{'Read' if _i == 0 else 'Mate'} ({_m['length']} bases): each k-mer and the transcripts that contain it"))
            _items.append(mo.md("*The read goes to the transcripts **every** k-mer found agrees on (the intersection); k-mers not in "
                                "the index (a sequencing error, a variant) are skipped, as kallisto does.*"))
            _items.append(code_panel(f"r = {python_call('explain_pseudoalignment', _v['index'], _read, **_kw)}"))
            _panel = mo.vstack(_items)
    except Exception as _err:
        _panel = error_panel(_err)
    _panel
    return


@app.cell
def crispr_page_ui(get_page, mo, set_page):
    mo.stop(get_page() != "crispr")
    btn_back_crispr = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_crispr = mo.ui.button(label="🔄 Refresh the file list")
    crispr_task = mo.ui.radio(options=["Count guides in reads", "Rank genes (MAGeCK RRA)"], value="Count guides in reads", inline=True, label="**What to do:**")
    crispr_folder = mo.ui.text(value=".", label="Data folder (its subfolders are searched too):", full_width=True)
    _ui = mo.vstack([
        btn_back_crispr,
        mo.md("<br>\n\n# ✂️ 23. CRISPR Screens (MAGeCK)"),
        mo.callout(mo.md(
            "In a pooled CRISPR knockout screen, each cell gets one guide RNA that switches one gene off; the cells then grow, "
            "with or without a drug. Sequencing the guides before and after tells which knockouts were lost (**negative "
            "selection**: genes the cells need, or that protect them from the drug) and which were enriched (**positive "
            "selection**: genes whose loss gives resistance). **Counting** finds each read's guide in the library, as `mageck "
            "count` does; **ranking** normalizes the counts, scores every guide and combines the guides of each gene by robust "
            "rank aggregation (RRA), as `mageck test` does — here with the same results, line by line. Experimental, for "
            "teaching and research."), kind="info"),
        crispr_folder,
        mo.hstack([crispr_task, btn_refresh_crispr], justify="space-between", align="center"),
    ])
    _ui
    return btn_back_crispr, btn_refresh_crispr, crispr_folder, crispr_task


@app.cell
def crispr_page_example(crispr_task, example_picker, get_page, mo):
    mo.stop(get_page() != "crispr")
    crispr_example = example_picker("crispr", crispr_task.value)
    crispr_example
    return (crispr_example,)


@app.cell
def crispr_page_form(btn_refresh_crispr, crispr_example, crispr_folder, crispr_task, example_for, files_under, get_page, label_of, make_form, mo):
    mo.stop(get_page() != "crispr")
    btn_refresh_crispr  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("crispr", crispr_task.value, crispr_example.value)
    _folder = crispr_folder.value
    _tables = files_under(_folder, "table")
    if crispr_task.value == "Count guides in reads":
        _mm = {"exact only (as MAGeCK)": 0, "also one mismatch, when only one guide is that close": 1}
        form_crispr = make_form({
            "library": mo.ui.dropdown(options=_tables, value=_ex.get("library"), label="Guide library (id, sequence, gene):", full_width=True),
            "fastqs": mo.ui.multiselect(options=files_under(_folder, "fastq"), value=_ex.get("fastqs", []), label="Reads, one FASTQ per sample:", full_width=True),
            "labels": mo.ui.text(value=_ex.get("labels", ""), placeholder="blank: sample1, sample2...", label="Sample names, in the same order:", full_width=True),
            "mismatches": mo.ui.dropdown(options=_mm, value=label_of(_mm, _ex.get("mismatches", 0), "exact only (as MAGeCK)"), label="Matching:", full_width=True),
            "reverse": mo.ui.checkbox(value=_ex.get("reverse", False), label="The reads carry the guide's reverse complement"),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: guide_counts.txt next to the library", label="Count table:", full_width=True),
        }, "✂️ Count")
    else:
        _norms = {"median ratio (MAGeCK's default)": "median", "total counts": "total", "none": "none"}
        form_crispr = make_form({
            "counts": mo.ui.dropdown(options=_tables, value=_ex.get("counts"), label="Count table (sgRNA, Gene, one column per sample):", full_width=True),
            "treatment": mo.ui.text(value=_ex.get("treatment", ""), placeholder="HL60.final, KBM7.final", label="Treatment samples (after the screen):", full_width=True),
            "control": mo.ui.text(value=_ex.get("control", ""), placeholder="blank: every other sample", label="Control samples (before):", full_width=True),
            "norm": mo.ui.dropdown(options=_norms, value=label_of(_norms, _ex.get("norm", "median"), "median ratio (MAGeCK's default)"), label="Normalization:", full_width=True),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: show here only", label="Also write MAGeCK's tables, with this prefix:", full_width=True),
        }, "✂️ Rank the genes")
    mo.vstack([_x for _x in (_about, form_crispr) if _x is not None])
    return (form_crispr,)


@app.cell
async def crispr_page_backend(asyncio, beside, chosen, code_panel, crispr_task, error_panel, facts, form_crispr, get_page, mmg, mo, python_call, scatter_chart, warn, words):
    mo.stop(get_page() != "crispr")
    mo.stop(form_crispr.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_crispr.value
    try:
        if crispr_task.value == "Count guides in reads":
            _fq = [_p for _p in _v["fastqs"] if chosen(_p)]
            mo.stop(not chosen(_v["library"]), warn("Select the guide library."))
            mo.stop(not _fq, warn("Select the FASTQ files, one per sample."))
            _out = (_v["output"] or "").strip() or beside(_v["library"], "guide_counts.txt")
            _kw = {"output": _out}
            if words(_v["labels"]):
                _kw["labels"] = words(_v["labels"])
            if int(_v["mismatches"]):
                _kw["mismatches"] = 1
            if _v["reverse"]:
                _kw["reverse_complement"] = True
            with mo.status.spinner("Counting the guides..."):
                _r = await asyncio.to_thread(mmg.crispr_count, _v["library"], _fq, **_kw)
            _items = [mo.callout(mo.md("### ✂️ Guides counted\n" + facts([
                ("Guides in the library", f"`{len(_r['sgrna']):,}`"), ("Count table", f"`{_out}`")])), kind="success"),
                mo.ui.table([{
                    "sample": _s["label"], "reads": _s["reads"], "mapped": _s["mapped"], "mapped %": f"{_s['mapped_pct']:.1f}%",
                    "guides with no read": _s["zero_guides"], "Gini index": round(_s["gini"], 3),
                    "guide offset": ", ".join(str(_o) for _f in _s["files"] for _o in _f["offsets"]),
                } for _s in _r["samples"]], selection=None, label="Per sample"),
                mo.md("*A Gini index near 0 means even guide counts; a high one, a few guides taking most reads (after strong selection, or a poor library).*")]
            for _w in _r["warnings"][:5]:
                _items.append(warn(_w))
            _items.append(code_panel(f"r = {python_call('crispr_count', _v['library'], _fq, **_kw)}"))
        else:
            mo.stop(not chosen(_v["counts"]), warn("Select the count table."))
            mo.stop(not words(_v["treatment"]), warn("Type the treatment samples (column names of the count table)."))
            _kw = {}
            if words(_v["control"]):
                _kw["control"] = words(_v["control"])
            if _v["norm"] != "median":
                _kw["norm_method"] = _v["norm"]
            if (_v["output"] or "").strip():
                _kw["output_prefix"] = _v["output"].strip()
            with mo.status.spinner("Scoring the guides and ranking the genes (with permutations)..."):
                _r = await asyncio.to_thread(mmg.crispr_test, _v["counts"], words(_v["treatment"]), **_kw)
            _items = [mo.callout(mo.md("### ✂️ Genes ranked\n" + facts([
                ("Treatment against control", f"{', '.join(_r['treatment'])} against {', '.join(_r['control'])}"),
                ("Guides / genes", f"`{_r['guides_used']:,}` of `{_r['guides_total']:,}` guides, `{_r['genes_total']:,}` genes"),
                ("Size factors", ", ".join(f"{_n} {_f:.3f}" for _n, _f in _r["size_factors"])),
                ("Normalization", _r["norm_used"]),
            ])), kind="success")]
            for _side, _title, _kind in (("neg", "Negative selection (depleted: genes the cells need)", "down"),
                                         ("pos", "Positive selection (enriched)", "up")):
                _rp = _r["rank_plot"][_side]
                _rows = [{"rank": _rp["rank"][_i], "gene": _rp["gene"][_i], "score": f"{_rp['score'][_i]:.3g}",
                          "FDR": f"{_rp['fdr'][_i]:.3g}", "LFC": round(_rp["lfc"][_i], 3)} for _i in range(min(200, len(_rp["gene"])))]
                _items.append(mo.ui.table(_rows, selection=None, pagination=True, page_size=10, label=_title))
                _pts = [(_rp["rank"][_i], _rp["minus_log10_score"][_i], _kind if _rp["fdr"][_i] < 0.25 else "ns",
                         f"{_rp['gene'][_i]}: FDR {_rp['fdr'][_i]:.3g}") for _i in range(len(_rp["gene"]))]
                _lab = [(_rp["rank"][_i], _rp["minus_log10_score"][_i], _rp["gene"][_i]) for _i in range(min(5, len(_rp["gene"])))]
                _chart = scatter_chart(_pts, "gene rank", "−log10 RRA score", f"{_title.split(' (')[0]}: genes by rank (coloured: FDR < 0.25)", labels=_lab)
                if _chart is not None:
                    _items.append(_chart)
            for _w in _r["warnings"][:5]:
                _items.append(warn(_w))
            _items.append(code_panel(f"r = {python_call('crispr_test', _v['counts'], words(_v['treatment']), **_kw)}",
                                     'print(r["top_negative"][:10])'))
        _panel = mo.vstack(_items)
    except Exception as _err:
        _panel = error_panel(_err)
    _panel
    return


@app.cell
def meta_page_ui(get_page, mo, set_page):
    mo.stop(get_page() != "meta")
    btn_back_meta = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_meta = mo.ui.button(label="🔄 Refresh the file list")
    meta_task = mo.ui.radio(options=["Which species, and how much", "Classify the reads", "Sketch reference genomes"],
                            value="Which species, and how much", inline=True, label="**What to do:**")
    meta_folder = mo.ui.text(value=".", label="Data folder (its subfolders are searched too):", full_width=True)
    _ui = mo.vstack([
        btn_back_meta,
        mo.md("<br>\n\n# 🦠 24. Metagenomics: Species & Abundance"),
        mo.callout(mo.md(
            "Which organisms are in a sample sequenced as a whole (shotgun metagenomics: stool, soil, a culture, a clinical "
            "sample), and in what proportions. Each reference genome is reduced to a **sketch** — a fixed fraction of its "
            "k-mers, chosen by their hash, as sourmash does (with the same hashes) — and the sample's sketch is covered by the "
            "references, the best match first (`gather`); or each read is assigned to the reference it shares most k-mers "
            "with (as Kraken does). Human (host) reads are counted apart. Only the organisms of the references can be found: "
            "a species missing from them shows as unassigned. Experimental, for teaching and research: not a clinical "
            "diagnostic."), kind="info"),
        meta_folder,
        mo.hstack([meta_task, btn_refresh_meta], justify="space-between", align="center"),
    ])
    _ui
    return btn_back_meta, btn_refresh_meta, meta_folder, meta_task


@app.cell
def meta_page_example(example_picker, get_page, meta_task, mo):
    mo.stop(get_page() != "meta")
    meta_example = example_picker("meta", meta_task.value)
    meta_example
    return (meta_example,)


@app.cell
def meta_page_form(NO_FILE, btn_refresh_meta, example_for, files_under, get_page, label_of, make_form, meta_example, meta_folder, meta_task, mo):
    mo.stop(get_page() != "meta")
    btn_refresh_meta  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("meta", meta_task.value, meta_example.value)
    _folder = meta_folder.value
    _seqs = files_under(_folder, "fastq", "fasta")
    _sigs = files_under(_folder, "sig")
    if meta_task.value == "Which species, and how much":
        form_meta = make_form({
            "samples": mo.ui.multiselect(options=_seqs, value=_ex.get("samples", []), label="Sample (both files of a pair):", full_width=True),
            "databases": mo.ui.multiselect(options=_sigs, value=_ex.get("databases", []), label="Reference sketches (.sig):", full_width=True),
            "threshold": mo.ui.number(start=0, stop=10**8, step=10_000, value=_ex.get("threshold", 50_000), label="Report a match only above this many bases of overlap:"),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: shown here only", label="Also save the table (TSV):", full_width=True),
        }, "🦠 Identify")
    elif meta_task.value == "Classify the reads":
        _mates = {"(none: single-end)": "none", **{_f: _f for _f in _seqs if _f != NO_FILE}}
        form_meta = make_form({
            "input": mo.ui.dropdown(options=_seqs, value=_ex.get("input"), label="Reads (R1):", full_width=True),
            "input2": mo.ui.dropdown(options=_mates, value=label_of(_mates, _ex.get("input2", "none"), "(none: single-end)"), label="Mate reads (R2):", full_width=True),
            "databases": mo.ui.multiselect(options=_sigs, value=_ex.get("databases", []), label="Reference sketches (a dense one, e.g. scaled 10-20):", full_width=True),
            "min_hits": mo.ui.number(start=1, stop=1000, step=1, value=_ex.get("min_hits", 2), label="Minimum shared sketch k-mers to classify a read:"),
            "host_free": mo.ui.text(value=_ex.get("host_free", ""), placeholder="blank: do not write them", label="Write the reads that are not host to (R1):", full_width=True),
        }, "🦠 Classify")
    else:
        form_meta = make_form({
            "genomes": mo.ui.multiselect(options=files_under(_folder, "fasta"), value=_ex.get("genomes", []), label="Genomes (one FASTA each):", full_width=True),
            "host": mo.ui.checkbox(value=_ex.get("host", False), label="They are host genomes (e.g. human): their reads are counted apart"),
            "scaled": mo.ui.number(start=1, stop=100_000, step=10, value=_ex.get("scaled", 1000), label="Scaled (keep one k-mer hash in this many):"),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="references.sig", label="Sketch database:", full_width=True),
        }, "🦠 Sketch")
    mo.vstack([_x for _x in (_about, form_meta) if _x is not None])
    return (form_meta,)


@app.cell
async def meta_page_backend(asyncio, chosen, code_panel, error_panel, facts, form_meta, get_page, meta_task, mmg, mo, os, python_call, warn):
    mo.stop(get_page() != "meta")
    mo.stop(form_meta.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_meta.value

    def _bars(rows):
        """Rows of (label, percent, grey?) as horizontal bars in HTML."""
        _top = max([_p for _, _p, _ in rows] + [1e-9])
        _html = []
        for _lab, _p, _grey in rows:
            _w = max(0.5, 100 * _p / _top)
            _c = "#898781" if _grey else "#2a78d6"
            _html.append(f'<div style="display:grid; grid-template-columns: minmax(120px, 34%) 1fr 64px; gap: 8px; align-items: center; '
                         f'font: 13px system-ui, sans-serif; margin: 3px 0;"><span>{_lab}</span>'
                         f'<span style="background:{_c}; height: 12px; width: {_w:.1f}%; border-radius: 0 4px 4px 0;"></span>'
                         f'<span style="text-align: right;">{_p:.2f}%</span></div>')
        return mo.Html("".join(_html))

    try:
        if meta_task.value == "Which species, and how much":
            _s = [_p for _p in _v["samples"] if chosen(_p)]
            _d = [_p for _p in _v["databases"] if chosen(_p)]
            mo.stop(not _s or not _d, warn("Select the sample and at least one sketch database (make one with **Sketch reference genomes**)."))
            _kw = {}
            if int(_v["threshold"]) != 50_000:
                _kw["threshold_bp"] = int(_v["threshold"])
            if (_v["output"] or "").strip():
                _kw["output"] = _v["output"].strip()
            with mo.status.spinner("Sketching the sample and matching it..."):
                _r = await asyncio.to_thread(mmg.gather, _s, _d, **_kw)
            _m = _r["matches"]
            _rows = [(_x["name"] + (" (host)" if _x["host"] else ""), 100 * (_x["relative_abundance"] or _x["f_unique_weighted"] or 0), _x["host"]) for _x in _m]
            _items = [mo.callout(mo.md("### 🦠 Species found\n" + facts([
                ("Sample", f"`{_r['sequences']:,}` {'reads' if _r['reads'] else 'sequences'}, `{_r['bases']:,}` bases, `{_r['sample_hashes']:,}` sketch k-mers"),
                ("Matches", f"`{len(_m)}` of `{_r['references']}` references"),
                ("Unassigned", f"{100 * (_r['unassigned_weighted'] or 0):.1f}% of the sample's k-mers (weighted by abundance)"),
            ])), kind="success")]
            if _m:
                _items.append(mo.md("**Relative abundance** (host reads apart):"))
                _items.append(_bars(_rows))
                _items.append(mo.ui.table([{
                    "rank": _x["rank"], "reference": _x["name"], "relative abundance": f"{100 * _x['relative_abundance']:.2f}%" if _x["relative_abundance"] is not None else "host",
                    "containment": f"{_x['containment']:.3f}", "ANI": f"{_x['ani']:.4f}", "overlap (bp)": _x["unique_intersect_bp"],
                } for _x in _m], selection=None, label="Matches"))
            _items.append(code_panel(f"r = {python_call('gather', _s, _d, **_kw)}"))
        elif meta_task.value == "Classify the reads":
            mo.stop(not chosen(_v["input"]), warn("Select the reads."))
            _d = [_p for _p in _v["databases"] if chosen(_p)]
            mo.stop(not _d, warn("Select at least one sketch database."))
            _kw = {}
            if _v["input2"] != "none":
                _kw["input2"] = _v["input2"]
            if int(_v["min_hits"]) != 2:
                _kw["min_hits"] = int(_v["min_hits"])
            if (_v["host_free"] or "").strip():
                _kw["host_free"] = _v["host_free"].strip()
            with mo.status.spinner("Classifying the reads..."):
                _r = await asyncio.to_thread(mmg.classify_reads, _v["input"], _d, **_kw)
            _items = [mo.callout(mo.md("### 🦠 Reads classified\n" + facts([
                (("Pairs" if _r["paired"] else "Reads"), f"`{_r['reads']:,}`: `{_r['classified']:,}` classified, `{_r['host']:,}` host, `{_r['unclassified']:,}` unclassified"),
            ])), kind="success")]
            _rows = [(_x["name"] + (" (host)" if _x["host"] else ""), 100 * (_x["fraction_of_classified"] or 0), _x["host"]) for _x in _r["rows"] if _x["reads"]]
            if _rows:
                _items.append(mo.md("**Share of the classified reads:**"))
                _items.append(_bars(sorted(_rows, key=lambda _t: -_t[1])))
            for _w in _r["warnings"][:5]:
                _items.append(warn(_w))
            _items.append(code_panel(f"r = {python_call('classify_reads', _v['input'], _d, **_kw)}"))
        else:
            _g = [_p for _p in _v["genomes"] if chosen(_p)]
            mo.stop(not _g, warn("Select the genomes."))
            _out = (_v["output"] or "").strip() or "references.sig"
            _kw = {}
            if _v["host"]:
                _kw["host"] = [True] * len(_g)
            if int(_v["scaled"]) != 1000:
                _kw["scaled"] = int(_v["scaled"])
            with mo.status.spinner("Sketching..."):
                _r = await asyncio.to_thread(mmg.sketch_references, _g, _out, **_kw)
            _items = [mo.callout(mo.md("### 🦠 Sketches written\n" + facts([("Database", f"`{os.path.abspath(_out)}`"), ("k / scaled", f"{_r['k']} / {_r['scaled']}")])), kind="success"),
                      mo.ui.table([{"genome": _x["name"], "bases": _x["bases"], "sketch k-mers": _x["hashes"], "host": "yes" if _x["host"] else ""}
                                   for _x in _r["genomes"]], selection=None, label="Genomes"),
                      code_panel(f"r = {python_call('sketch_references', _g, _out, **_kw)}")]
        _panel = mo.vstack(_items)
    except Exception as _err:
        _panel = error_panel(_err)
    _panel
    return


@app.cell
def amr_page_ui(example_picker, get_page, mo, set_page):
    mo.stop(get_page() != "amr")
    btn_back_amr = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_amr = mo.ui.button(label="🔄 Refresh the file list")
    amr_folder = mo.ui.text(value=".", label="Data folder (its subfolders are searched too):", full_width=True)
    amr_example = example_picker("amr")
    _ui = mo.vstack([
        btn_back_amr,
        mo.md("<br>\n\n# 💉 25. Antimicrobial Resistance Genes"),
        mo.callout(mo.md(
            "Which known resistance genes a bacterial isolate carries — carbapenemases such as KPC and NDM, mecA of MRSA, "
            "ESBLs, aminoglycoside and other genes — as hospital surveillance looks for them. The isolate's assembly (contigs) "
            "or reads are searched for every gene of the **NCBI AMRFinderPlus** catalog (public domain; download it from "
            "ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/), with the best allele per "
            "locus, its identity and its coverage. A gene's presence does not always mean resistance, and resistance can come "
            "from mutations this search does not check: **it does not replace susceptibility testing.** Experimental, for "
            "teaching and research."), kind="info"),
        amr_folder,
        mo.hstack([amr_example, btn_refresh_amr], justify="space-between", align="center"),
    ])
    _ui
    return amr_example, amr_folder, btn_back_amr, btn_refresh_amr


@app.cell
def amr_page_form(NO_FILE, amr_example, amr_folder, btn_refresh_amr, example_for, files_under, get_page, label_of, make_form, mo):
    mo.stop(get_page() != "amr")
    btn_refresh_amr  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("amr", None, amr_example.value)
    _folder = amr_folder.value
    _seqs = files_under(_folder, "fasta", "fastq")
    _catalogs = [_f for _f in files_under(_folder, "fasta") if "amr" in _f.lower()] or [NO_FILE]
    _tables = {"(none: no drug class)": "none", **{_f: _f for _f in files_under(_folder, "table") if "catalog" in _f.lower()}}
    _modes = {"automatic: FASTA = assembly, FASTQ = reads": "auto", "assembly (contigs)": "assembly", "reads": "reads"}
    form_amr = make_form({
        "inputs": mo.ui.multiselect(options=_seqs, value=_ex.get("inputs", []), label="Isolate (assembly, or both files of reads):", full_width=True),
        "catalog": mo.ui.dropdown(options=_catalogs, value=_ex.get("catalog", _catalogs[0] if _catalogs[0] != NO_FILE else None), label="Catalog (AMR_CDS.fa):", full_width=True),
        "table": mo.ui.dropdown(options=_tables, value=label_of(_tables, _ex.get("table", next((_v for _v in _tables.values() if _v != "none"), "none")), "(none: no drug class)"),
                                label="Catalog table (ReferenceGeneCatalog.txt), for the drug class:", full_width=True),
        "min_identity": mo.ui.number(start=0.5, stop=1.0, step=0.01, value=_ex.get("min_identity", 0.9), label="Minimum identity:"),
        "min_coverage": mo.ui.number(start=0.1, stop=1.0, step=0.05, value=_ex.get("min_coverage", 0.9), label="Minimum coverage of the gene:"),
        "plus": mo.ui.checkbox(value=_ex.get("plus", True), label="Also the 'plus' genes (virulence, stress, biocides)"),
        "mode": mo.ui.dropdown(options=_modes, value=label_of(_modes, _ex.get("mode", "auto"), "automatic: FASTA = assembly, FASTQ = reads"), label="Input:", full_width=True),
        "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: shown here only", label="Also save the table (TSV):", full_width=True),
    }, "💉 Search")
    mo.vstack([_x for _x in (_about, form_amr) if _x is not None])
    return (form_amr,)


@app.cell
async def amr_page_backend(asyncio, chosen, code_panel, error_panel, facts, form_amr, get_page, mmg, mo, python_call, warn):
    mo.stop(get_page() != "amr")
    mo.stop(form_amr.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_amr.value
    _in = [_p for _p in _v["inputs"] if chosen(_p)]
    mo.stop(not _in, warn("Select the isolate's assembly or reads."))
    mo.stop(not chosen(_v["catalog"]), warn("Select the catalog (AMR_CDS.fa of AMRFinderPlus)."))
    _kw = {}
    if _v["table"] != "none":
        _kw["catalog_table"] = _v["table"]
    if float(_v["min_identity"]) != 0.9:
        _kw["min_identity"] = float(_v["min_identity"])
    if float(_v["min_coverage"]) != 0.9:
        _kw["min_coverage"] = float(_v["min_coverage"])
    if not _v["plus"]:
        _kw["include_plus"] = False
    if _v["mode"] != "auto":
        _kw["mode"] = _v["mode"]
    if (_v["output"] or "").strip():
        _kw["output"] = _v["output"].strip()
    try:
        with mo.status.spinner("Searching the catalog's genes..."):
            _r = await asyncio.to_thread(mmg.find_resistance_genes, _in, _v["catalog"], **_kw)
        _hits = _r["hits"]
        _core = [_h for _h in _hits if _h["scope"] == "core"]
        _classes = {}
        for _h in _core:
            _classes.setdefault(_h["class"] or "other", []).append(_h["gene"])
        _items = [mo.callout(mo.md("### 💉 Resistance genes\n" + facts([
            ("Searched", f"`{_r['sequences']:,}` sequence(s), `{_r['bases']:,}` bases, against `{_r['catalog_genes']:,}` catalog genes ({_r['mode']})"),
            ("Found", f"`{len(_core)}` core resistance gene(s)" + (f", and `{len(_hits) - len(_core)}` 'plus' gene(s)" if len(_hits) > len(_core) else "")),
        ]) + ("\n\n**By drug class:** " + "; ".join(f"**{_c}**: {', '.join(sorted(set(_g)))}" for _c, _g in sorted(_classes.items())) if _classes else "")),
                             kind="success" if _hits else "warn")]
        if _hits:
            _items.append(mo.ui.table([{
                "gene": _h["gene"], "product": _h["product"], "class": _h["class"], "subclass": _h["subclass"], "scope": _h["scope"],
                "identity": f"{100 * _h['identity']:.2f}%", "coverage": f"{100 * _h['coverage']:.1f}%",
                "where": f"{_h['contig']}:{_h['start']}-{_h['end']} ({_h['strand']})", "status": _h["status"],
            } for _h in sorted(_hits, key=lambda _h: (_h["scope"] != "core", _h["class"] or "", _h["gene"]))],
                selection=None, pagination=True, page_size=15, label="Genes found"))
        _items.append(code_panel(f"r = {python_call('find_resistance_genes', _in, _v['catalog'], **_kw)}"))
        _panel = mo.vstack(_items)
    except Exception as _err:
        _panel = error_panel(_err)
    _panel
    return


@app.cell
def outbreak_page_ui(get_page, mo, set_page):
    mo.stop(get_page() != "outbreak")
    btn_back_outbreak = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_outbreak = mo.ui.button(label="🔄 Refresh the file list")
    outbreak_task = mo.ui.radio(options=["Distances and a tree", "Mutations of viral genomes", "In-silico PCR"],
                                value="Distances and a tree", inline=True, label="**What to do:**")
    outbreak_folder = mo.ui.text(value=".", label="Data folder (its subfolders are searched too):", full_width=True)
    _ui = mo.vstack([
        btn_back_outbreak,
        mo.md("<br>\n\n# 🌳 26. Outbreaks: Trees, Viral Mutations & PCR"),
        mo.callout(mo.md(
            "Tools of genomic epidemiology. **Distances and a tree**: how close isolates are — each genome reduced to a "
            "k-mer sketch, as Mash and sourmash do (the same hashes as sourmash), compared pair by pair (Jaccard, Mash "
            "distance, ANI) — and the tree they make (neighbor joining or UPGMA): isolates of one outbreak cluster tightly. "
            "**Mutations of viral genomes**: each genome aligned to its reference, its changes named as Nextclade names "
            "them, with the amino-acid changes per gene (`S:D614G`...). **In-silico PCR**: where the primers of an assay "
            "bind in genomes or reads, with their mismatches — what a PCR panel would detect, and which variants could "
            "escape it. A k-mer tree is a quick look, not a SNP phylogeny: transmission chains need a core-genome SNP "
            "analysis. Experimental, for teaching and research: not a clinical or public-health report."), kind="info"),
        outbreak_folder,
        mo.hstack([outbreak_task, btn_refresh_outbreak], justify="space-between", align="center"),
    ])
    _ui
    return btn_back_outbreak, btn_refresh_outbreak, outbreak_folder, outbreak_task


@app.cell
def outbreak_page_example(example_picker, get_page, mo, outbreak_task):
    mo.stop(get_page() != "outbreak")
    outbreak_example = example_picker("outbreak", outbreak_task.value)
    outbreak_example
    return (outbreak_example,)


@app.cell
def outbreak_page_form(NO_FILE, btn_refresh_outbreak, example_for, files_under, get_page, label_of, make_form, mo,
                       outbreak_example, outbreak_folder, outbreak_task):
    mo.stop(get_page() != "outbreak")
    btn_refresh_outbreak  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("outbreak", outbreak_task.value, outbreak_example.value)
    _folder = outbreak_folder.value
    if outbreak_task.value == "Distances and a tree":
        _methods = {"neighbor joining (assumes no molecular clock)": "nj", "UPGMA (assumes every lineage evolves at the same rate)": "upgma"}
        form_outbreak = make_form({
            "genomes": mo.ui.multiselect(options=files_under(_folder, "fasta", "sig"), value=_ex.get("genomes", []),
                                         label="Genomes (FASTA, one per file) or sketch databases (.sig):", full_width=True),
            "records": mo.ui.checkbox(value=_ex.get("records", False), label="Each FASTA record is a genome (a multi-FASTA of viral genomes)"),
            "method": mo.ui.dropdown(options=_methods, value=label_of(_methods, _ex.get("method", "nj"), "neighbor joining (assumes no molecular clock)"),
                                     label="Tree:", full_width=True),
            "outgroup": mo.ui.text(value=_ex.get("outgroup", ""), placeholder="blank: root at the midpoint",
                                   label="Root on this genome (an outgroup, e.g. the oldest lineage):", full_width=True),
            "scaled": mo.ui.number(start=0, stop=100_000, step=10, value=_ex.get("scaled", 0),
                                   label="Scaled (keep one k-mer in this many; 0 = chosen from the genome sizes):"),
            "tree_output": mo.ui.text(value=_ex.get("tree_output", ""), placeholder="blank: shown here only", label="Also save the tree (Newick, .nwk):", full_width=True),
            "matrix_output": mo.ui.text(value=_ex.get("matrix_output", ""), placeholder="blank: shown here only",
                                        label="Also save the distance matrix (.tsv, or .phy for PHYLIP):", full_width=True),
        }, "🌳 Build the tree")
    elif outbreak_task.value == "Mutations of viral genomes":
        _fastas = files_under(_folder, "fasta")
        _gffs = {"(none: nucleotide changes only)": "none", **{_f: _f for _f in files_under(_folder, "gff") if _f != NO_FILE}}
        form_outbreak = make_form({
            "reference": mo.ui.dropdown(options=_fastas, value=_ex.get("reference"), label="Reference genome (FASTA):", full_width=True),
            "annotation": mo.ui.dropdown(options=_gffs, value=label_of(_gffs, _ex.get("annotation", "none"), "(none: nucleotide changes only)"),
                                         label="Its annotation (GFF3), for the amino-acid changes:", full_width=True),
            "genomes": mo.ui.multiselect(options=_fastas, value=_ex.get("genomes", []), label="Genomes (FASTA; every record is a genome):", full_width=True),
            "split": mo.ui.checkbox(value=_ex.get("split", True), label="Report ORF1ab as ORF1a and ORF1b (Nextclade's names)"),
            "summary": mo.ui.text(value=_ex.get("summary", ""), placeholder="blank: shown here only", label="Also save one line per genome (TSV):", full_width=True),
            "output_dir": mo.ui.text(value=_ex.get("output_dir", ""), placeholder="blank: not saved",
                                     label="Also save every change of every genome (one TSV each) in this folder:", full_width=True),
        }, "🧬 Find the mutations")
    else:
        form_outbreak = make_form({
            "panel": mo.ui.dropdown(options=files_under(_folder, "table"), value=_ex.get("panel"),
                                    label="Primer panel (TSV: target, forward, reverse, probe, size, source):", full_width=True),
            "inputs": mo.ui.multiselect(options=files_under(_folder, "fasta", "fastq"), value=_ex.get("inputs", []),
                                        label="Genomes (FASTA) or reads (FASTQ):", full_width=True),
            "records": mo.ui.checkbox(value=_ex.get("records", False), label="Each FASTA record is a genome (a multi-FASTA of viral genomes)"),
            "max_mm": mo.ui.number(start=0, stop=8, step=1, value=_ex.get("max_mm", 3), label="Mismatches allowed in a primer:"),
            "three": mo.ui.number(start=0, stop=10, step=1, value=_ex.get("three", 3), label="Bases at the 3' end of a primer that must match:"),
            "amplicons_output": mo.ui.text(value=_ex.get("amplicons_output", ""), placeholder="blank: shown here only",
                                           label="Also save the products (TSV):", full_width=True),
        }, "🧪 Run the PCR")
    mo.vstack([_x for _x in (_about, form_outbreak) if _x is not None])
    return (form_outbreak,)


@app.cell
async def outbreak_page_backend(asyncio, chosen, code_panel, error_panel, facts, figs, figure_panel, form_outbreak, get_page, mmg, mo,
                                outbreak_task, python_call, tree_chart, warn):
    mo.stop(get_page() != "outbreak")
    mo.stop(form_outbreak.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_outbreak.value
    try:
        if outbreak_task.value == "Distances and a tree":
            _g = [_p for _p in _v["genomes"] if chosen(_p)]
            mo.stop(not _g, warn("Select the genomes (FASTA) or a sketch database (.sig)."))
            _sig = [_p for _p in _g if _p.lower().endswith(".sig")]
            mo.stop(bool(_sig) and len(_sig) != len(_g), warn("Give either genomes (FASTA) or sketch databases (.sig), not both."))
            _kw = {}
            if _v["method"] != "nj":
                _kw["method"] = _v["method"]
            if (_v["outgroup"] or "").strip():
                _kw["outgroup"] = _v["outgroup"].strip()
            if int(_v["scaled"] or 0):
                _kw["scaled"] = int(_v["scaled"])
            for _k in ("tree_output", "matrix_output"):
                if (_v[_k] or "").strip():
                    _kw[_k] = _v[_k].strip()
            if _sig:
                _fn = "sketch_distances"
                with mo.status.spinner("Comparing every pair of sketches..."):
                    _r = await asyncio.to_thread(mmg.sketch_distances, _g, **_kw)
            else:
                _fn = "genome_distances"
                if _v["records"]:
                    _kw["records_as_genomes"] = True
                with mo.status.spinner("Sketching the genomes and comparing every pair..."):
                    _r = await asyncio.to_thread(mmg.genome_distances, _g, **_kw)
            _names, _d, _ani = _r["names"], _r["distance"], _r["ani"]
            _n = len(_names)
            _pairs = [(_i, _j) for _i in range(_n) for _j in range(_i + 1, _n)]
            _close = min(_pairs, key=lambda _t: _d[_t[0]][_t[1]]) if _pairs else None
            _method = "UPGMA" if _r["tree"]["method"] == "upgma" else "Neighbor-joining"
            _root = (f"rooted on `{_kw['outgroup']}`" if "outgroup" in _kw else
                     "rooted at the midpoint" if _method != "UPGMA" else "rooted by the method")
            _items = [mo.callout(mo.md("### 🌳 Tree built\n" + facts([
                ("Genomes", f"`{_n}`"),
                ("Sketches", f"k = {_r['k']}, scaled = {_r['scaled']}" + (" (every k-mer)" if _r["scaled"] == 1 else "")),
                ("Tree", f"{_method}, {_root}"),
                ("Closest pair", f"`{_names[_close[0]]}` and `{_names[_close[1]]}`: ANI {100 * _ani[_close[0]][_close[1]]:.2f}%"
                 if _close else "—"),
            ])), kind="success")]
            _chart = tree_chart(_r["tree"], f"{_method} tree of {_n} genomes")
            if _chart is not None:
                _items.append(_chart)
            _rows = []
            for _i in range(_n):
                _j = min((_j for _j in range(_n) if _j != _i), key=lambda _j: _d[_i][_j], default=None)
                if _j is not None:
                    _rows.append({"genome": _names[_i], "nearest": _names[_j], "ANI": f"{100 * _ani[_i][_j]:.2f}%",
                                  "Mash distance": round(_d[_i][_j], 5), "Jaccard": round(_r["jaccard"][_i][_j], 4)})
            _items.append(mo.ui.table(_rows, selection=None, pagination=True, page_size=25, label="The nearest genome of each"))
            if _n <= 12:
                _items.append(mo.ui.table([{"genome": _names[_i], **{_names[_j]: f"{100 * _ani[_i][_j]:.2f}" for _j in range(_n)}} for _i in range(_n)],
                                          selection=None, label="ANI (%) of every pair"))
            _items.append(mo.accordion({"The tree in Newick format": mo.plain_text(_r["tree"]["newick"])}))
            for _w in _r["warnings"][:5]:
                _items.append(warn(_w))
            _items.append(code_panel(f"r = {python_call(_fn, _g, **_kw)}"))
        elif outbreak_task.value == "Mutations of viral genomes":
            mo.stop(not chosen(_v["reference"]), warn("Select the reference genome."))
            _gs = [_p for _p in _v["genomes"] if chosen(_p)]
            mo.stop(not _gs, warn("Select the genomes."))
            mo.stop(_v["reference"] in _gs, warn("The reference cannot be one of the genomes."))
            _kw = {}
            if _v["annotation"] != "none":
                _kw["annotation"] = _v["annotation"]
            if not _v["split"]:
                _kw["split_orf1ab"] = False
            if (_v["summary"] or "").strip():
                _kw["summary_output"] = _v["summary"].strip()
            if (_v["output_dir"] or "").strip():
                _kw["output_dir"] = _v["output_dir"].strip()
            with mo.status.spinner("Aligning each genome to the reference..."):
                _r = await asyncio.to_thread(mmg.viral_mutations, _v["reference"], _gs, **_kw)
            _G = _r["genomes"]
            _ok = [_x for _x in _G if _x["aligned"]]
            _items = [mo.callout(mo.md("### 🧬 Mutations found\n" + facts([
                ("Genomes", f"`{len(_ok)}` of `{len(_G)}` aligned to `{_r['reference_sequences'][0]['id'] if _r['reference_sequences'] else _v['reference']}`"),
                ("Genes", f"`{len(_r['genes'])}` coding sequences in the annotation" if _r["genes"] else "no annotation: nucleotide changes only"),
            ])), kind="success" if _ok else "warn")]

            def _aa(_x):
                return [_c["label"] for _c in _x["aa_substitutions"] + _x["aa_deletions"]] + [_c["label"] for _c in _x["aa_insertions"]]

            if _ok:
                _items.append(figure_panel(
                    figs.viral_genomes(_r), "viral_mutations",
                    "Each genome against the reference, as Nextclade draws it. A dark grey block is a stretch the sequencing did "
                    "not read (Ns, an amplicon dropout), not a deletion: a real deletion is red and usually shared by a lineage."))

            _items.append(mo.ui.table([{
                "genome": _x["name"],
                "coverage": f"{100 * _x['coverage']:.1f}%" if _x["aligned"] else "not aligned",
                "substitutions": len(_x["substitutions"]), "deletions": len(_x["deletions"]), "insertions": len(_x["insertions"]),
                "missing (N) bases": sum(_b - _a + 1 for _a, _b in _x["missing"]),
                "amino-acid changes": len(_aa(_x)),
                "changes": ", ".join(_aa(_x)[:14]) + (" …" if len(_aa(_x)) > 14 else ""),
            } for _x in _G], selection=None, pagination=True, page_size=25, label="Each genome"))
            if _r["genes"] and len(_ok) > 1:
                _count = {}
                for _x in _ok:
                    for _c in set(_aa(_x)):
                        _count[_c] = _count.get(_c, 0) + 1
                _shared = sorted(((_c, _k) for _c, _k in _count.items() if _k > 1), key=lambda _t: (-_t[1], _t[0]))
                if _shared:
                    _items.append(mo.md("**Amino-acid changes shared by several genomes** — changes that define a lineage are "
                                        "carried by all its descendants:"))
                    _items.append(mo.ui.table([{"change": _c, "genomes": _k, "share": f"{100 * _k / len(_ok):.0f}%",
                                                "in": ", ".join(_x["name"] for _x in _ok if _c in _aa(_x))[:160]}
                                               for _c, _k in _shared], selection=None, pagination=True, page_size=20,
                                              label="Shared changes"))
            _bad = [_x for _x in _G if _x["problem"]]
            for _x in _bad[:5]:
                _items.append(warn(f"`{_x['name']}`: {_x['problem']}"))
            for _w in _r["warnings"][:5]:
                _items.append(warn(_w))
            _items.append(code_panel(f"r = {python_call('viral_mutations', _v['reference'], _gs, **_kw)}"))
        else:
            mo.stop(not chosen(_v["panel"]), warn("Select the primer panel."))
            _ins = [_p for _p in _v["inputs"] if chosen(_p)]
            mo.stop(not _ins, warn("Select the genomes or reads."))
            _kw = {}
            if _v["records"]:
                _kw["records_as_genomes"] = True
            if int(_v["max_mm"]) != 3:
                _kw["max_mismatches"] = int(_v["max_mm"])
            if int(_v["three"]) != 3:
                _kw["three_prime_exact"] = int(_v["three"])
            if (_v["amplicons_output"] or "").strip():
                _kw["amplicons_output"] = _v["amplicons_output"].strip()
            with mo.status.spinner("Looking for the primer sites..."):
                _r = await asyncio.to_thread(mmg.in_silico_pcr, _v["panel"], _ins, **_kw)
            _t = _r["targets"]
            _pos = sum(1 for _row in _r["rows"] for _p in _row["present"] if _p)
            _items = [mo.callout(mo.md("### 🧪 PCR done\n" + facts([
                ("Assays", ", ".join(f"`{_x}`" for _x in _t)),
                ("Inputs", f"`{len(_r['rows'])}` genome(s) or sample(s)"),
                ("Positive", f"`{_pos}` of `{len(_r['rows']) * len(_t)}` (input, assay) pairs"),
                ("Rules", f"up to {_r['max_mismatches']} mismatches per primer, none in the last {_r['three_prime_exact']} bases at 3'"),
            ])), kind="success")]
            _items.append(mo.ui.table([{"genome or sample": _row["name"],
                                        **{_x: (f"✓ {_c}" if _p else "—") for _x, _p, _c in zip(_t, _row["present"], _row["counts"])}}
                                       for _row in _r["rows"]], selection=None, pagination=True, page_size=25,
                                      label="Detected (✓ with the number of products, or of reads with a site)"))
            if _r["amplicons"]:
                _items.append(mo.ui.table([{
                    "genome": _a["genome"], "assay": _a["target"], "where": f"{_a['sequence_id']}:{_a['start']:,}-{_a['end']:,} ({_a['strand']})",
                    "size": _a["size"], "mismatches F / R": f"{_a['forward_mismatches']} / {_a['reverse_mismatches']}",
                    "probe mismatches": "—" if _a["probe_mismatches"] is None else _a["probe_mismatches"],
                    "detected": "yes" if _a["detected"] else ("size out of range" if not _a["in_range"] else "probe does not bind"),
                } for _a in _r["amplicons"]], selection=None, pagination=True, page_size=15, label="Products"))
            if _r["read_counts"]:
                _items.append(mo.ui.table([{"sample": _c["sample"], "assay": _c["target"], "reads": _c["reads"], "with a site": _c["with_site"],
                                            "whole product": _c["amplicon"], "present": "yes" if _c["present"] else ""}
                                           for _c in _r["read_counts"]], selection=None, label="Reads"))
            for _w in _r["warnings"][:5]:
                _items.append(warn(_w))
            _items.append(code_panel(f"r = {python_call('in_silico_pcr', _v['panel'], _ins, **_kw)}"))
        _panel = mo.vstack(_items)
    except Exception as _err:
        _panel = error_panel(_err)
    _panel
    return


@app.cell
def amplicon_page_ui(get_page, mo, set_page):
    mo.stop(get_page() != "amplicon")
    btn_back_amplicon = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_amplicon = mo.ui.button(label="🔄 Refresh the file list")
    amplicon_task = mo.ui.radio(options=["Amplicon variants (ASVs)", "Taxonomy", "Diversity", "Train a classifier"],
                                value="Amplicon variants (ASVs)", inline=True, label="**What to do:**")
    amplicon_folder = mo.ui.text(value=".", label="Data folder (its subfolders are searched too):", full_width=True)
    _ui = mo.vstack([
        btn_back_amplicon,
        mo.md("<br>\n\n# 🌿 27. 16S Amplicons: Taxonomy & Diversity"),
        mo.callout(mo.md(
            "Who is in a microbial community, from the reads of one gene (16S rRNA) amplified from every bacterium in the "
            "sample. **Amplicon variants**: each read pair merged, reads with errors dropped, the rest denoised into exact "
            "sequence variants (ASVs, as UNOISE3 / vsearch) and chimeras removed. **Taxonomy**: each ASV named by the "
            "naive Bayesian classifier of the RDP (as DADA2 and QIIME 2 use it). **Diversity**: within each sample "
            "(alpha: richness, Shannon, Simpson...) and between samples (beta: Bray-Curtis or Jaccard, principal "
            "coordinates, PERMANOVA between groups). 16S tells genera, rarely species, and counts copies of a gene, not "
            "cells. Experimental, for teaching and research."), kind="info"),
        amplicon_folder,
        mo.hstack([amplicon_task, btn_refresh_amplicon], justify="space-between", align="center"),
    ])
    _ui
    return amplicon_folder, amplicon_task, btn_back_amplicon, btn_refresh_amplicon


@app.cell
def amplicon_page_example(amplicon_task, example_picker, get_page, mo):
    mo.stop(get_page() != "amplicon")
    amplicon_example = example_picker("amplicon", amplicon_task.value)
    amplicon_example
    return (amplicon_example,)


@app.cell
def amplicon_page_form(NO_FILE, amplicon_example, amplicon_folder, amplicon_task, btn_refresh_amplicon, example_for, files_under,
                       get_page, label_of, make_form, mo):
    mo.stop(get_page() != "amplicon")
    btn_refresh_amplicon  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("amplicon", amplicon_task.value, amplicon_example.value)
    _folder = amplicon_folder.value
    _tables = [_f for _f in files_under(_folder, "table") if _f != NO_FILE]
    if amplicon_task.value == "Amplicon variants (ASVs)":
        _sheets = {"(none: the FASTQ files below)": "none", **{_f: _f for _f in files_under(_folder, "sheet") if _f != NO_FILE}}
        form_amplicon = make_form({
            "sheet": mo.ui.dropdown(options=_sheets, value=label_of(_sheets, _ex.get("sheet", "none"), "(none: the FASTQ files below)"),
                                    label="Sample sheet (name, R1, R2 per line; mothur's stability.files works):", full_width=True),
            "fastqs": mo.ui.multiselect(options=files_under(_folder, "fastq"), value=_ex.get("fastqs", []),
                                        label="Or the FASTQ files, paired by their _R1 / _R2 names:", full_width=True),
            "primer_f": mo.ui.text(value=_ex.get("primer_f", ""), placeholder="blank: the reads hold no primer (as in the mothur SOP)",
                                   label="Forward primer to remove (5' to 3', IUPAC codes allowed):", full_width=True),
            "primer_r": mo.ui.text(value=_ex.get("primer_r", ""), placeholder="blank: none", label="Reverse primer to remove:", full_width=True),
            "max_ee": mo.ui.number(start=0.1, stop=10, step=0.1, value=_ex.get("max_ee", 1.0), label="Most expected errors in a merged read:"),
            "min_size": mo.ui.number(start=1, stop=1000, step=1, value=_ex.get("min_size", 8),
                                     label="Fewest copies for a sequence to become an ASV (UNOISE3's minsize):"),
            "chimeras": mo.ui.checkbox(value=_ex.get("chimeras", True), label="Remove chimeras (as vsearch --uchime3_denovo)"),
            "out": mo.ui.text(value=_ex.get("out", ""), placeholder="blank: next to the sample sheet (or the first FASTQ)",
                              label="Folder for asvs.fasta, asv_table.tsv and asv_tracking.tsv:", full_width=True),
        }, "🌿 Find the ASVs")
    elif amplicon_task.value == "Taxonomy":
        _models = [_f for _f in files_under(_folder, "model", "fasta") if _f != NO_FILE] or [NO_FILE]
        _tab = {"(none: no composition chart)": "none", **{_f: _f for _f in _tables}}
        _conf = {"0.8 (the RDP's default)": 0.8, "0.5 (DADA2's minBoot 50, for short reads)": 0.5}
        _ranks = ["Phylum", "Class", "Order", "Family", "Genus"]
        form_amplicon = make_form({
            "fasta": mo.ui.dropdown(options=files_under(_folder, "fasta"), value=_ex.get("fasta"), label="ASVs (FASTA):", full_width=True),
            "model": mo.ui.dropdown(options=_models, value=_ex.get("model"),
                                    label="Classifier (a .model from 'Train a classifier', or a reference FASTA):", full_width=True),
            "table": mo.ui.dropdown(options=_tab, value=label_of(_tab, _ex.get("table", "none"), "(none: no composition chart)"),
                                    label="Feature table (ASV x sample), for the composition of each sample:", full_width=True),
            "confidence": mo.ui.dropdown(options=_conf, value=label_of(_conf, _ex.get("confidence", 0.8), "0.8 (the RDP's default)"),
                                         label="Confidence needed at each rank:", full_width=True),
            "rank": mo.ui.dropdown(options=_ranks, value=_ex.get("rank", "Family"), label="Rank of the composition chart:", full_width=True),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: taxonomy.tsv next to the ASVs", label="Taxonomy table:", full_width=True),
        }, "🌿 Classify")
    elif amplicon_task.value == "Diversity":
        _meta = {"(none: no groups)": "none", **{_f: _f for _f in _tables}}
        _metrics = {"Bray-Curtis (abundances)": "braycurtis", "Jaccard (presence / absence)": "jaccard"}
        form_amplicon = make_form({
            "table": mo.ui.dropdown(options=_tables or [NO_FILE], value=_ex.get("table"), label="Feature table (ASV x sample):", full_width=True),
            "metadata": mo.ui.dropdown(options=_meta, value=label_of(_meta, _ex.get("metadata", "none"), "(none: no groups)"),
                                       label="Sample groups (a table: sample, then columns such as time):", full_width=True),
            "column": mo.ui.text(value=_ex.get("column", ""), placeholder="blank: the second column", label="Group column:", full_width=True),
            "metric": mo.ui.dropdown(options=_metrics, value=label_of(_metrics, _ex.get("metric", "braycurtis"), "Bray-Curtis (abundances)"),
                                     label="Distance between samples:", full_width=True),
            "rarefy": mo.ui.number(start=0, stop=10**7, step=100, value=_ex.get("rarefy", 0),
                                   label="Rarefy every sample to this many reads first (0 = no):"),
        }, "🌿 Measure")
    else:
        form_amplicon = make_form({
            "reference": mo.ui.dropdown(options=files_under(_folder, "fasta"), value=_ex.get("reference"),
                                        label="Reference (FASTA with the lineage in each header, as the DADA2 files of RDP and GTDB):", full_width=True),
            "primer_f": mo.ui.text(value=_ex.get("primer_f", ""), placeholder="blank: whole sequences",
                                   label="Train on the amplified region only: forward primer:", full_width=True),
            "primer_r": mo.ui.text(value=_ex.get("primer_r", ""), placeholder="blank: whole sequences", label="Reverse primer:", full_width=True),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: classifier.model next to the reference", label="Model file:", full_width=True),
        }, "🌿 Train")
    mo.vstack([_x for _x in (_about, form_amplicon) if _x is not None])
    return (form_amplicon,)


@app.cell
async def amplicon_page_backend(amplicon_task, asyncio, beside, chosen, code_panel, error_panel, facts, figs, figure_panel, form_amplicon,
                                get_page, mmg, mo, os, python_call, warn):
    mo.stop(get_page() != "amplicon")
    mo.stop(form_amplicon.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_amplicon.value
    # The categorical colours of the charts (validated: the FASTQ report and dataviz palette), grey for "Other".
    _PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
    _GREY = "#898781"

    def _legend(names, colours):
        return ('<div style="display:flex; flex-wrap:wrap; gap: 6px 16px; font: 12px system-ui, sans-serif; margin: 6px 0;">'
                + "".join(f'<span><span style="display:inline-block; width:10px; height:10px; border-radius:2px; background:{_c}; '
                          f'margin-right:5px; vertical-align:-1px"></span>{_n}</span>' for _n, _c in zip(names, colours)) + "</div>")

    def _stacked(samples, taxa, fractions):
        """One 100% bar per sample, a segment per taxon, 2 px gaps; legend above."""
        _cols = [(_GREY if _t in ("Other", "Unassigned") or _i >= len(_PAL) else _PAL[_i]) for _i, _t in enumerate(taxa)]
        _rows = []
        for _j, _s in enumerate(samples):
            _seg = []
            for _i, _t in enumerate(taxa):
                _f = fractions[_i][_j]
                if _f > 0:
                    _seg.append(f'<span title="{_t}: {100 * _f:.1f}%" style="flex: {_f} 0 0; background:{_cols[_i]}; height: 14px;'
                                f' margin-right: 2px; border-radius: 2px;"></span>')
            _rows.append(f'<div style="display:grid; grid-template-columns: minmax(70px, 18%) 1fr; gap: 8px; align-items:center; '
                         f'font: 12px system-ui, sans-serif; margin: 3px 0;"><span>{_s}</span><span style="display:flex;">{"".join(_seg)}</span></div>')
        return mo.Html(_legend(taxa, _cols) + "".join(_rows))

    def _scatter(points, xl, yl, title):
        """PCoA: (x, y, group, label) points, one colour per group, labels on hover."""
        import math
        if not points:
            return None
        _groups = list(dict.fromkeys(_p[2] for _p in points))
        _col = {_g: (_PAL[_i] if _i < len(_PAL) else _GREY) for _i, _g in enumerate(_groups)}
        _xs, _ys = [_p[0] for _p in points], [_p[1] for _p in points]
        _x0, _x1, _y0, _y1 = min(_xs), max(_xs), min(_ys), max(_ys)
        _dx, _dy = (_x1 - _x0) or 1, (_y1 - _y0) or 1
        _x0, _x1, _y0, _y1 = _x0 - 0.08 * _dx, _x1 + 0.08 * _dx, _y0 - 0.08 * _dy, _y1 + 0.08 * _dy
        W, H, L, R, T, B = 560, 420, 50, 16, 16, 40
        X = lambda v: L + (v - _x0) / (_x1 - _x0) * (W - L - R)
        Y = lambda v: T + (H - T - B) - (v - _y0) / (_y1 - _y0) * (H - T - B)
        _parts = [f'<line x1="{L}" x2="{W - R}" y1="{Y(0):.1f}" y2="{Y(0):.1f}" stroke="currentColor" stroke-opacity="0.15"/>',
                  f'<line y1="{T}" y2="{H - B}" x1="{X(0):.1f}" x2="{X(0):.1f}" stroke="currentColor" stroke-opacity="0.15"/>']
        for _x, _y, _g, _lab in points:
            _parts.append(f'<circle cx="{X(_x):.1f}" cy="{Y(_y):.1f}" r="6" fill="{_col[_g]}" stroke="var(--background, #fff)" stroke-width="2">'
                          f'<title>{_lab} ({_g})</title></circle>')
        _parts.append(f'<text x="{L + (W - L - R) / 2:.1f}" y="{H - 8}" text-anchor="middle" fill="currentColor" fill-opacity="0.8">{xl}</text>')
        _parts.append(f'<text transform="translate(14,{T + (H - T - B) / 2:.1f}) rotate(-90)" text-anchor="middle" fill="currentColor" fill-opacity="0.8">{yl}</text>')
        _svg = (f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width:{W}px; font: 12px system-ui, sans-serif;" role="img" aria-label="{title}">'
                + "".join(_parts) + "</svg>")
        return mo.vstack([mo.md(f"**{title}**"), mo.Html(_legend(_groups, [_col[_g] for _g in _groups]) + _svg)])

    def _curves(curves, group_of, title):
        """Rarefaction curves: observed ASVs against reads drawn, one line per sample, coloured by group."""
        if not curves:
            return None
        _groups = list(dict.fromkeys(group_of.get(_c["sample"], "—") for _c in curves))
        _col = {_g: (_PAL[_i] if _i < len(_PAL) and _g != "—" else _GREY) for _i, _g in enumerate(_groups)}
        _xm = max(max(_c["depths"] or [1]) for _c in curves) or 1
        _ym = max(max(_c["observed"] or [1]) for _c in curves) or 1
        W, H, L, R, T, B = 640, 320, 54, 16, 14, 40
        X = lambda v: L + v / _xm * (W - L - R)
        Y = lambda v: T + (H - T - B) - v / _ym * (H - T - B)
        _parts = []
        for _f in (0.25, 0.5, 0.75, 1.0):
            _parts.append(f'<line x1="{L}" x2="{W - R}" y1="{Y(_f * _ym):.1f}" y2="{Y(_f * _ym):.1f}" stroke="currentColor" stroke-opacity="0.12"/>'
                          f'<text x="{L - 6}" y="{Y(_f * _ym) + 4:.1f}" text-anchor="end" fill="currentColor" fill-opacity="0.7">{_f * _ym:.0f}</text>')
        for _c in curves:
            _pts = " ".join(f"{X(_d):.1f},{Y(_o):.1f}" for _d, _o in zip(_c["depths"], _c["observed"]))
            _g = group_of.get(_c["sample"], "—")
            _parts.append(f'<polyline points="{_pts}" fill="none" stroke="{_col[_g]}" stroke-width="2" stroke-opacity="0.85">'
                          f'<title>{_c["sample"]} ({_g})</title></polyline>')
        _parts.append(f'<text x="{L + (W - L - R) / 2:.1f}" y="{H - 8}" text-anchor="middle" fill="currentColor" fill-opacity="0.8">reads drawn (max {_xm:,})</text>')
        _parts.append(f'<text transform="translate(14,{T + (H - T - B) / 2:.1f}) rotate(-90)" text-anchor="middle" fill="currentColor" fill-opacity="0.8">ASVs seen</text>')
        _svg = (f'<svg viewBox="0 0 {W} {H}" width="100%" style="max-width:{W}px; font: 12px system-ui, sans-serif;" role="img" aria-label="{title}">'
                + "".join(_parts) + "</svg>")
        return mo.vstack([mo.md(f"**{title}**"), mo.Html(_legend(_groups, [_col[_g] for _g in _groups]) + _svg)])

    try:
        if amplicon_task.value == "Amplicon variants (ASVs)":
            _fq = [_p for _p in _v["fastqs"] if chosen(_p)]
            mo.stop(_v["sheet"] == "none" and not _fq, warn("Select a sample sheet or the FASTQ files."))
            _samples = _v["sheet"] if _v["sheet"] != "none" else _fq
            _first = _v["sheet"] if _v["sheet"] != "none" else _fq[0]
            _dir = (_v["out"] or "").strip() or os.path.dirname(beside(_first, "x")) or "."
            os.makedirs(_dir, exist_ok=True)
            _fa, _tab, _trk = (os.path.join(_dir, _n) for _n in ("asvs.fasta", "asv_table.tsv", "asv_tracking.tsv"))
            _kw = {"tracking": _trk}
            for _k, _name in (("primer_f", "primer_forward"), ("primer_r", "primer_reverse")):
                if (_v[_k] or "").strip():
                    _kw[_name] = _v[_k].strip()
            if float(_v["max_ee"]) != 1.0:
                _kw["max_ee"] = float(_v["max_ee"])
            if int(_v["min_size"]) != 8:
                _kw["min_size"] = int(_v["min_size"])
            if not _v["chimeras"]:
                _kw["remove_chimeras"] = False
            with mo.status.spinner("Merging, filtering, denoising and removing chimeras..."):
                _r = await asyncio.to_thread(mmg.amplicon_asvs, _samples, _fa, _tab, **_kw)
            _S = _r["samples"]
            _in = sum(_x["input"] for _x in _S)
            _kept = sum(_x["nonchimeric"] for _x in _S)
            _items = [mo.callout(mo.md("### 🌿 ASVs found\n" + facts([
                ("Samples", f"`{len(_S)}`, `{_in:,}` read pairs"),
                ("ASVs", f"`{len(_r['asvs']):,}` (from `{_r['uniques']:,}` distinct reads; `{len(_r['chimeras'])}` chimeras removed)"),
                ("Reads kept", f"`{_kept:,}` ({100 * _kept / max(_in, 1):.1f}%)"),
                ("Written", f"`{_fa}`, `{_tab}`, `{_trk}`"),
            ])), kind="success")]
            _items.append(mo.ui.table([{"sample": _x["sample"], "pairs": _x["input"], "merged": _x["merged"], "filtered": _x["filtered"],
                                        "denoised": _x["denoised"], "non-chimeric": _x["nonchimeric"],
                                        "kept": f"{100 * _x['nonchimeric'] / max(_x['input'], 1):.1f}%"} for _x in _S],
                                      selection=None, pagination=True, page_size=25, label="Reads kept at each step"))
            _items.append(mo.ui.table([{"ASV": _a["name"], "length": len(_a["sequence"]), "reads": _a["reads"], "sequence": _a["sequence"][:60] + "…"}
                                       for _a in _r["asvs"][:200]], selection=None, pagination=True, page_size=10, label="ASVs, the most abundant first"))
            for _w in _r["warnings"][:5]:
                _items.append(warn(_w))
            _items.append(code_panel(f"r = {python_call('amplicon_asvs', _samples, _fa, _tab, **_kw)}"))
        elif amplicon_task.value == "Taxonomy":
            mo.stop(not chosen(_v["fasta"]) or not chosen(_v["model"]), warn("Select the ASVs and the classifier."))
            _out = (_v["output"] or "").strip() or beside(_v["fasta"], "taxonomy.tsv")
            _kw = {}
            if float(_v["confidence"]) != 0.8:
                _kw["min_confidence"] = float(_v["confidence"])
            with mo.status.spinner("Classifying each ASV (100 bootstrap trials)..."):
                _r = await asyncio.to_thread(mmg.classify_asvs, _v["fasta"], _v["model"], _out, **_kw)
            _A = _r["assignments"]
            _ranks = _r["ranks"]
            _items = [mo.callout(mo.md("### 🌿 ASVs classified\n" + facts([
                ("ASVs", f"`{len(_A):,}`, confidence ≥ {_r['min_confidence']}"),
                ("Named down to", ", ".join(f"{_k}: {_n}" for _k, _n in zip(_ranks, _r["assigned_per_rank"]))),
                ("Taxonomy table", f"`{_out}`"),
            ])), kind="success")]
            _items.append(mo.ui.table([{"ASV": _a["id"], **{_k: (_a["lineage"][_i] if _i < len(_a["lineage"]) else "") for _i, _k in enumerate(_ranks)},
                                        "confidence (last rank)": f"{_a['confidence'][len(_a['lineage']) - 1]:.2f}" if _a["lineage"] else ""}
                                       for _a in _A], selection=None, pagination=True, page_size=15, label="Taxonomy"))
            _code = [f"r = {python_call('classify_asvs', _v['fasta'], _v['model'], _out, **_kw)}"]
            # Reads per ASV (from the feature table, when given) for the sunburst; else each ASV counts once.
            _reads = {}
            if _v["table"] != "none":
                try:
                    with open(_v["table"], encoding="utf-8") as _fh:
                        _fh.readline()
                        for _line in _fh:
                            _f = _line.rstrip("\n").split("\t")
                            _reads[_f[0]] = sum(float(_z) for _z in _f[1:] if _z.strip())
                except (OSError, ValueError):
                    _reads = {}
            _skip = 1 if _ranks and _ranks[0].lower() in ("kingdom", "domain") else 0
            _lin = [(tuple(_a["lineage"][_skip:]), _reads.get(_a["id"], 0 if _reads else 1)) for _a in _A]
            _items.append(figure_panel(
                figs.sunburst(_lin, title="Taxonomy" + (" (share of the reads)" if _reads else " (share of the ASVs: no table given)"),
                              ranks=[_k.lower() for _k in _ranks[_skip:]]),
                "taxonomy_sunburst",
                "From the centre outwards, one ring per rank; an arc's angle is its share. 'unclassified' is where the "
                "classifier's confidence fell below the threshold: the reference has no close relative at that rank."))
            if _v["table"] != "none":
                _c = await asyncio.to_thread(mmg.taxa_composition, _v["table"], _out, rank=_v["rank"], top=8)
                _items.append(mo.md(f"**Composition of each sample, by {_c['rank'].lower()}** (the 8 most abundant; the rest as Other):"))
                _items.append(_stacked(_c["samples"], _c["taxa"], _c["fractions"]))
                _code.append(f"c = {python_call('taxa_composition', _v['table'], _out, rank=_v['rank'], top=8)}")
            _items.append(code_panel(*_code))
        elif amplicon_task.value == "Diversity":
            mo.stop(not chosen(_v["table"]), warn("Select the feature table."))
            _kw_a, _kw_b = {}, {}
            if int(_v["rarefy"] or 0):
                _kw_a["rarefy_depth"] = _kw_b["rarefy_depth"] = int(_v["rarefy"])
            if _v["metric"] != "braycurtis":
                _kw_b["metric"] = _v["metric"]
            if _v["metadata"] != "none":
                _kw_b["metadata"] = _v["metadata"]
                if (_v["column"] or "").strip():
                    _kw_b["group_column"] = _v["column"].strip()
            with mo.status.spinner("Alpha and beta diversity, principal coordinates, PERMANOVA..."):
                _a = await asyncio.to_thread(mmg.alpha_diversity, _v["table"], **_kw_a)
                _b = await asyncio.to_thread(mmg.beta_diversity, _v["table"], **_kw_b)
            _group = {}
            if _v["metadata"] != "none":
                with open(_v["metadata"], encoding="utf-8") as _fh:
                    _rows = [_l.rstrip("\n").split("\t") for _l in _fh if _l.strip()]
                _ci = _rows[0].index(_kw_b["group_column"]) if _kw_b.get("group_column") in _rows[0] else 1
                _group = {_r0[0]: _r0[_ci] for _r0 in _rows[1:] if len(_r0) > _ci}
            _pm = _b["permanova"]
            _facts = [("Samples", f"`{len(_b['samples'])}`" + (f" (rarefied to {_b['rarefied_to']:,} reads)" if _b["rarefied_to"] else "")),
                      ("Distance", {"braycurtis": "Bray-Curtis", "jaccard": "Jaccard"}.get(_b["metric"], _b["metric"]))]
            if _pm:
                _facts.append(("PERMANOVA", f"{_pm['column']}: " + ", ".join(f"{_g} ({_n})" for _g, _n in zip(_pm["groups"], _pm["sizes"]))
                               + f" — pseudo-F = {_pm['f']:.2f}, R² = {_pm['r2']:.2f}, p = {_pm['p_value']:.3g} ({_pm['permutations']} permutations)"))
            _items = [mo.callout(mo.md("### 🌿 Diversity measured\n" + facts(_facts)), kind="success")]
            _items.append(mo.ui.table([{"sample": _x["sample"], "group": _group.get(_x["sample"], ""), "reads": _x["depth"], "observed": _x["observed"],
                                        "Chao1": round(_x["chao1"], 1), "Shannon": round(_x["shannon"], 3), "Simpson": round(_x["simpson"], 3),
                                        "Pielou": round(_x["pielou"], 3), "Good's coverage": f"{100 * _x['goods_coverage']:.1f}%"}
                                       for _x in _a["rows"]], selection=None, pagination=True, page_size=25, label="Alpha diversity"))
            _ch = _curves(_a["curves"], _group, "Rarefaction curves: a curve that levels off has been sequenced deeply enough")
            if _ch is not None:
                _items.append(_ch)
            _pc = _b["pcoa"]
            if _pc["coordinates"] and len(_pc["coordinates"][0]) >= 2:
                _pe = _pc["proportion_explained"]
                _items.append(figure_panel(
                    figs.group_scatter([(_c[0], _c[1], _group.get(_s, "—"), _s) for _s, _c in zip(_b["samples"], _pc["coordinates"])],
                                       f"PCo1 ({100 * _pe[0]:.1f}%)", f"PCo2 ({100 * _pe[1]:.1f}%)",
                                       "Principal coordinates: samples close together have similar communities"),
                    "beta_pcoa",
                    "Each ellipse holds about 95% of its group if the points were normally spread (groups of three or more). "
                    "Separate ellipses with little overlap go with a significant PERMANOVA; a group much more spread than "
                    "another can also make PERMANOVA significant (it tests location and dispersion together)."))
            _items.append(code_panel(f"a = {python_call('alpha_diversity', _v['table'], **_kw_a)}",
                                     f"b = {python_call('beta_diversity', _v['table'], **_kw_b)}"))
        else:
            mo.stop(not chosen(_v["reference"]), warn("Select the reference FASTA."))
            _out = (_v["output"] or "").strip() or beside(_v["reference"], "classifier.model")
            _kw = {}
            for _k, _name in (("primer_f", "primer_forward"), ("primer_r", "primer_reverse")):
                if (_v[_k] or "").strip():
                    _kw[_name] = _v[_k].strip()
            mo.stop(len(_kw) == 1, warn("Give both primers, or neither."))
            with mo.status.spinner("Training (8-mer word probabilities per genus)..."):
                _r = await asyncio.to_thread(mmg.train_classifier, _v["reference"], _out, **_kw)
            _items = [mo.callout(mo.md("### 🌿 Classifier trained\n" + facts([
                ("Reference", f"`{_r['sequences']:,}` sequences" + (f", `{_r['without_region']:,}` without the primers' region (left out)" if _r.get("without_region") else "")),
                ("Ranks", ", ".join(_r["ranks"])),
                ("Genera (classes)", f"`{_r['classes']:,}`"),
                ("Model", f"`{_out}` ({_r['model_bytes'] / 1e6:.1f} MB), in {_r['seconds']:.1f} s"),
            ])), kind="success"),
                code_panel(f"r = {python_call('train_classifier', _v['reference'], _out, **_kw)}")]
        _panel = mo.vstack(_items)
    except Exception as _err:
        _panel = error_panel(_err)
    _panel
    return


@app.cell
def tmb_page_ui(example_picker, get_page, human_notice, mo, set_page):
    mo.stop(get_page() != "tmb")
    btn_back_tmb = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_tmb = mo.ui.button(label="🔄 Refresh the file list")
    tmb_folder = mo.ui.text(value=".", label="Data folder (its subfolders are searched too):", full_width=True)
    tmb_example = example_picker("tmb", None)
    _ui = mo.vstack([
        btn_back_tmb,
        mo.md("<br>\n\n# 🎯 28. Tumor Mutational Burden (TMB)"),
        human_notice(),
        mo.callout(mo.md(
            "How many somatic mutations a tumor carries per megabase of the DNA that was sequenced — one of the markers "
            "looked at before an immunotherapy. Give the tumor's somatic VCF (from a tumor/normal caller) and the BED of the "
            "region the test covers (a panel, an exome or the callable genome). Numbers from different panels, exomes and "
            "genomes are not directly comparable: each test calibrates its own cut-off (the FDA approval of "
            "pembrolizumab for TMB-high tumors used 10 mutations/Mb on one panel)."), kind="info"),
        tmb_folder,
        mo.hstack([tmb_example, btn_refresh_tmb], justify="space-between", align="end"),
    ])
    _ui
    return btn_back_tmb, btn_refresh_tmb, tmb_example, tmb_folder


@app.cell
def tmb_page_form(NO_FILE, btn_refresh_tmb, example_for, files_under, get_page, label_of, make_form, mo, tmb_example, tmb_folder):
    mo.stop(get_page() != "tmb")
    btn_refresh_tmb  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("tmb", None, tmb_example.value)
    _folder = tmb_folder.value
    _gffs = {"(none)": "none", **{_f: _f for _f in files_under(_folder, "gff") if _f != NO_FILE}}
    _refs = {"(none)": "none", **{_f: _f for _f in files_under(_folder, "fasta") if _f != NO_FILE}}
    form_tmb = make_form({
        "vcf": mo.ui.dropdown(options=files_under(_folder, "vcf"), value=_ex.get("vcf"), label="Somatic variants of the tumor (VCF):", full_width=True),
        "bed": mo.ui.dropdown(options=files_under(_folder, "bed_any"), value=_ex.get("bed"), label="Region the test covers (BED):", full_width=True),
        "sample": mo.ui.text(value=_ex.get("sample", ""), placeholder="blank: the tumor sample, or the only one", label="Tumor sample (a column of the VCF):", full_width=True),
        "pass_only": mo.ui.checkbox(value=_ex.get("pass_only", True), label="Only variants that passed the caller's filters (FILTER = PASS)"),
        "min_depth": mo.ui.number(start=0, stop=10_000, step=1, value=_ex.get("min_depth", 10), label="Least read depth:"),
        "min_vaf": mo.ui.number(start=0, stop=1, step=0.01, value=_ex.get("min_vaf", 0.05), label="Least variant allele fraction:"),
        "annotation": mo.ui.dropdown(options=_gffs, value=label_of(_gffs, _ex.get("annotation", "none"), "(none)"),
                                     label="Gene annotation (GFF3), to count coding variants only:", full_width=True),
        "coding": mo.ui.checkbox(value=_ex.get("coding", False), label="Count only variants in coding sequences (needs the annotation)"),
        "reference": mo.ui.dropdown(options=_refs, value=label_of(_refs, _ex.get("reference", "none"), "(none)"),
                                    label="Genome (FASTA), to count non-synonymous variants only:", full_width=True),
        "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: shown here only", label="Also save every variant and whether it counted (TSV):", full_width=True),
    }, "🎯 Measure the TMB")
    mo.vstack([_x for _x in (_about, form_tmb) if _x is not None])
    return (form_tmb,)


@app.cell
async def tmb_page_backend(asyncio, chosen, code_panel, error_panel, facts, form_tmb, get_page, mmg, mo, python_call, warn):
    mo.stop(get_page() != "tmb")
    mo.stop(form_tmb.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_tmb.value
    try:
        mo.stop(not chosen(_v["vcf"]) or not chosen(_v["bed"]), warn("Select the tumor's VCF and the BED of the region."))
        mo.stop(_v["coding"] and _v["annotation"] == "none", warn("Counting coding variants only needs the gene annotation (GFF3)."))
        _kw = {}
        if (_v["sample"] or "").strip():
            _kw["tumor_sample"] = _v["sample"].strip()
        if not _v["pass_only"]:
            _kw["pass_only"] = False
        if int(_v["min_depth"]) != 10:
            _kw["min_depth"] = int(_v["min_depth"])
        if float(_v["min_vaf"]) != 0.05:
            _kw["min_vaf"] = float(_v["min_vaf"])
        if _v["annotation"] != "none":
            _kw["annotation"] = _v["annotation"]
            if _v["coding"]:
                _kw["coding_only"] = True
        if _v["reference"] != "none":
            mo.stop(_v["annotation"] == "none", warn("Counting non-synonymous variants needs the gene annotation (GFF3) too."))
            _kw["reference"] = _v["reference"]
            _kw["nonsynonymous_only"] = True
            _kw["coding_only"] = True
        if (_v["output"] or "").strip():
            _kw["output"] = _v["output"].strip()
        with mo.status.spinner("Counting the somatic variants in the region..."):
            _r = await asyncio.to_thread(mmg.tumor_mutational_burden, _v["vcf"], _v["bed"], **_kw)
        _what = ("non-synonymous coding" if _kw.get("nonsynonymous_only") else "coding" if _kw.get("coding_only") else "somatic")
        _tile = (f'<div style="font: 15px system-ui, sans-serif; margin: 8px 0;"><div style="font-size: 2.6em; font-weight: 700; line-height: 1.1;">'
                 f'{_r["tmb"]:.1f}</div><div>{_what} mutations per megabase'
                 f' ({_r["counted"]:,} over {_r["megabases"]:.2f} Mb)</div></div>')
        _items = [mo.Html(_tile), mo.callout(mo.md("### 🎯 TMB measured\n" + facts([
            ("Substitutions", f"`{_r['counted_snv'] + _r['counted_mnv']:,}` ({_r['tmb_snv']:.2f} per Mb)"),
            ("Insertions and deletions", f"`{_r['counted_indel']:,}` ({_r['tmb_indel']:.2f} per Mb)"),
            ("Region", f"{_r['bed_bp'] / 1e6:.2f} Mb in the BED" + (f"; {_r['coding_bp'] / 1e6:.2f} Mb of it coding" if _r["coding_bp"] else "")),
            ("Sample", f"`{_r['sample']}`" if _r["sample"] else "—"),
        ])), kind="success")]
        _why = [("not PASS", _r["excluded_filter"]), ("depth below the minimum", _r["excluded_depth"]),
                ("allele fraction below the minimum", _r["excluded_vaf"]), ("outside the region", _r["outside_region"]),
                ("not in a coding sequence", _r["not_coding"]), ("synonymous", _r["synonymous"]), ("symbolic allele", _r["symbolic"])]
        _items.append(mo.ui.table([{"not counted because": _w, "alleles": _n} for _w, _n in _why if _n]
                                  or [{"not counted because": "(every allele was counted)", "alleles": 0}], selection=None, label="Left out"))
        for _w in _r["warnings"][:5]:
            _items.append(warn(_w))
        _items.append(code_panel(f"r = {python_call('tumor_mutational_burden', _v['vcf'], _v['bed'], **_kw)}"))
        _panel = mo.vstack(_items)
    except Exception as _err:
        _panel = error_panel(_err)
    _panel
    return


@app.cell
def contam_page_ui(get_page, mo, set_page):
    mo.stop(get_page() != "contam")
    btn_back_contam = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_contam = mo.ui.button(label="🔄 Refresh the file list")
    contam_task = mo.ui.radio(options=["Screen samples for adventitious agents", "Build a screening database"],
                              value="Screen samples for adventitious agents", inline=True, label="**What to do:**")
    contam_folder = mo.ui.text(value=".", label="Data folder (its subfolders are searched too):", full_width=True)
    _ui = mo.vstack([
        btn_back_contam,
        mo.md("<br>\n\n# 🛡️ 31. Contaminant Screening: Adventitious Agents"),
        mo.callout(mo.md(
            "**⚠️ For research and teaching only.** This screen is not a validated method: releasing a batch of a "
            "biological product (GMP) requires a validated method run with validated software. ICH Q5A(R2) admits "
            "next-generation sequencing for adventitious viruses, but the validation is the laboratory's. A detection "
            "must be confirmed by an independent method (PCR, culture)."), kind="danger"),
        mo.callout(mo.md(
            "Is anything that should not be there in a biological material — a raw material, a cell bank, a culture "
            "harvest, a tissue extract? Its reads are compared with a **database** of the agents to look for (viruses, "
            "bacteria, mycoplasmas, fungi), the **host** (human, bovine, CHO...) and the **process controls**. Host reads "
            "are set aside first; each other read goes to the reference it matches best; and an organism is called "
            "**detected** only with evidence: its own reads, spread over its genome (**breadth**: a few reads piled on one "
            "region — a vector element, a gene shared with other organisms — are not a detection). An organism also present "
            "in a **negative control** (a reagent or extraction blank) is **background** unless the sample is 10 times "
            "above it. **Spike-ins** give the recovery and the limit of detection. The report records the SHA-256 of every "
            "file and of the database, every parameter, the date and the engine version."), kind="info"),
        contam_folder,
        mo.hstack([contam_task, btn_refresh_contam], justify="space-between", align="center"),
    ])
    _ui
    return btn_back_contam, btn_refresh_contam, contam_folder, contam_task


@app.cell
def contam_page_example(contam_task, example_picker, get_page, mo):
    mo.stop(get_page() != "contam")
    contam_example = example_picker("contam", contam_task.value)
    contam_example
    return (contam_example,)


@app.cell
def contam_page_form(NO_FILE, btn_refresh_contam, contam_example, contam_folder, contam_task, example_for, files_under, get_page, label_of, make_form, mo):
    mo.stop(get_page() != "contam")
    btn_refresh_contam  # a click re-runs this cell, which lists the folder again
    _ex, _about = example_for("contam", contam_task.value, contam_example.value)
    _folder = contam_folder.value
    if contam_task.value.startswith("Screen"):
        _spikes = {"(none)": "none", **{_f: _f for _f in files_under(_folder, "table") if _f != NO_FILE}}
        form_contam = make_form({
            "database": mo.ui.dropdown(options=files_under(_folder, "screen_db"), value=_ex.get("database"),
                                       label="Screening database (.mmgscreen: the task Build a screening database makes it):", full_width=True),
            "fastqs": mo.ui.multiselect(options=files_under(_folder, "fastq"), value=_ex.get("fastqs", []),
                                        label="Reads of the samples (FASTQ; x_R1 and x_R2 are paired automatically):", full_width=True),
            "negative_controls": mo.ui.text(value=_ex.get("negative_controls", ""), placeholder="reagent_blank",
                                            label="Negative controls (sample names: a sample is named after its files, without _R1/_R2):", full_width=True),
            "spike_ins": mo.ui.dropdown(options=_spikes, value=label_of(_spikes, _ex.get("spike_ins", "none"), "(none)"),
                                        label="Spike-ins and process controls (a table: sample, organism, amount, unit):", full_width=True),
            "min_reads": mo.ui.number(start=1, stop=1000, step=1, value=_ex.get("min_reads", 3), label="A detection needs at least this many own reads:"),
            "min_windows": mo.ui.number(start=1, stop=1000, step=1, value=_ex.get("min_windows", 2), label="... over at least this many genome windows:"),
            "min_dispersion": mo.ui.number(start=0, stop=1, step=0.05, value=_ex.get("min_dispersion", 0.5),
                                           label="... spread at least this much (windows hit / windows expected at random; 1 = fully spread):"),
            "min_blank_ratio": mo.ui.number(start=1, stop=1000, step=1, value=_ex.get("min_blank_ratio", 10),
                                            label="Above background when the sample is this many times the blank (reads per million):"),
            "report": mo.ui.text(value=_ex.get("report", ""), placeholder="blank: screen next to the database (screen.html, .tsv, .json)",
                                 label="Report name:", full_width=True),
        }, "🛡️ Screen")
    else:
        form_contam = make_form({
            "table": mo.ui.dropdown(options=files_under(_folder, "table"), value=_ex.get("table"),
                                    label="References: a table with the columns name, kind and file (kind: host, virus, bacterium, "
                                          "mycoplasma, fungus, other or spike; files relative to the table, several joined by ;):",
                                    full_width=True),
            "memory_mb": mo.ui.number(start=64, stop=262144, step=64, value=_ex.get("memory_mb", 1024),
                                      label="Memory for the host index, MB (less memory: more passes, never a failure):"),
            "output": mo.ui.text(value=_ex.get("output", ""), placeholder="blank: <table>.mmgscreen next to it", label="Database file:", full_width=True),
        }, "🧱 Build the database")
    mo.vstack([_x for _x in (_about, form_contam) if _x is not None])
    return (form_contam,)


@app.cell
async def contam_page_backend(asyncio, beside, chosen, code_panel, contam_task, error_panel, facts, figs, figure_panel, form_contam, get_page, mmg, mo, os, pair_fastqs, python_call, split_name, warn, words):
    mo.stop(get_page() != "contam")
    mo.stop(form_contam.value is None, mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_contam.value

    def _size(n):
        for _unit in ("bytes", "KB", "MB", "GB"):
            if n < 1024 or _unit == "GB":
                return f"{n:,.0f} {_unit}" if _unit == "bytes" else f"{n:,.1f} {_unit}"
            n /= 1024

    try:
        if contam_task.value.startswith("Screen"):
            mo.stop(not chosen(_v["database"]), warn("Select the screening database (.mmgscreen): build it with the task Build a screening database."))
            _fq = [_p for _p in _v["fastqs"] if chosen(_p)]
            mo.stop(not _fq, warn("Select the FASTQ files of the samples."))
            _samples = [(_n, _r1, _r2) if _r2 else (_n, _r1) for _n, _r1, _r2 in pair_fastqs(_fq)]
            _names = [_s[0] for _s in _samples]
            _neg = words(_v["negative_controls"])
            _unknown = [_x for _x in _neg if _x not in _names]
            mo.stop(bool(_unknown), warn(f"These negative controls are not among the samples: {', '.join(_unknown)}. "
                                         f"The samples are: {', '.join(_names)}."))
            _prefix = (_v["report"] or "").strip() or beside(_v["database"], "screen")
            for _ext in (".html", ".tsv", ".json"):
                if _prefix.lower().endswith(_ext):
                    _prefix = _prefix[: -len(_ext)]
            _kw = {}
            if _neg:
                _kw["negative_controls"] = _neg
            if _v["spike_ins"] != "none":
                _kw["spike_ins"] = _v["spike_ins"]
            for _key, _default in (("min_reads", 3), ("min_windows", 2)):
                if int(_v[_key] or 0) != _default:
                    _kw[_key] = int(_v[_key] or 0)
            for _key, _default in (("min_dispersion", 0.5), ("min_blank_ratio", 10.0)):
                if abs(float(_v[_key] or 0) - _default) > 1e-12:
                    _kw[_key] = float(_v[_key] or 0)
            _kw.update(output_html=_prefix + ".html", output_tsv=_prefix + ".tsv", output_json=_prefix + ".json")
            with mo.status.spinner(f"Screening {len(_samples)} sample(s): host reads apart, then each read to its best match..."):
                _r = await asyncio.to_thread(mmg.contaminant_screen, _samples, _v["database"], **_kw)
            _symbol = {"detected": "🔴 detected", "background": "⚪ background (also in the blank)",
                       "below_threshold": "🟡 below threshold", "not_detected": "— not detected"}
            _lines, _found = [], False
            for _s in _r["samples"]:
                _det = _s.get("detected") or []
                _found = _found or (bool(_det) and not _s.get("negative_control"))
                _lines.append((f"{_s['name']}" + (" (negative control)" if _s.get("negative_control") else ""),
                               f"{_s['fragments']:,} {'pairs' if _s.get('paired') else 'reads'}, host {100 * (_s.get('host_fraction') or 0):.2f}%; "
                               + ("detected: " + ", ".join(f"**{_d}**" for _d in _det) if _det else "nothing detected")))
            _db = _r.get("database") or {}
            _lines.append(("Database", f"`{os.path.basename(_db.get('path', _v['database']))}` (SHA-256 `{(_db.get('sha256') or '')[:16]}…`, "
                                       f"{len(_db.get('references') or [])} references)"))
            _items = [mo.callout(mo.md("### 🛡️ Screen\n" + facts(_lines)), kind="warn" if _found else "success")]
            if _r.get("warnings"):
                _items.append(mo.callout(mo.md("\n".join(f"* {_w}" for _w in _r["warnings"])), kind="warn"))
            for _s in _r["samples"]:
                _table = []
                for _o in _s["organisms"]:
                    if not (_o.get("reads") or _o.get("unique_reads") or _o.get("role") == "control" or _o.get("call") != "not_detected"):
                        continue
                    _call = _symbol.get(_o["call"], _o["call"])
                    if _o.get("role") == "control":
                        _call = "✅ recovered" if _o["call"] == "detected" else "⚠️ control " + _o["call"].replace("_", " ")
                    _table.append({
                        "organism": _o["name"], "kind": _o["kind"], "call": _call,
                        "reads (own)": f"{_o.get('reads', 0):,} ({_o.get('unique_reads', 0):,})",
                        "per million": f"{_o.get('rpm', 0):,.1f}",
                        "genome windows hit": f"{_o.get('windows', 0)} of {_o.get('genome_windows', 0)}",
                        "dispersion": "" if _o.get("dispersion") is None else f"{_o['dispersion']:.2f}",
                        "× blank": "" if _o.get("blank_ratio") is None else f"{_o['blank_ratio']:.3g}",
                        "recovery": "" if _o.get("recovery") is None else f"{100 * _o['recovery']:.0f}%",
                        "why": "; ".join(_o.get("reasons") or []),
                    })
                if _table:
                    _items.append(mo.ui.table(_table, selection=None, pagination=True, page_size=12,
                                              label=f"{_s['name']}: the organisms with reads"))
                    _items.append(figure_panel(
                        figs.genome_windows(_s["organisms"], title=f"{_s['name']}: where the reads fall on each genome"),
                        f"screen_{_s['name']}_windows",
                        "One cell per genome window. Reads of a real infection spread over the genome; a few lit windows in "
                        "one place are a vector element, an integrated fragment or a gene shared with other organisms — the "
                        "breadth rule turns those into 'below threshold'."))
                else:
                    _items.append(mo.md(f"*{_s['name']}: no read matched any agent or control of the database.*"))
            if _r.get("spike_results"):
                _items.append(mo.ui.table([{
                    "sample": _x["sample"], "organism": _x["organism"], "spiked": f"{_x['amount']:g} {_x['unit']}",
                    "reads": _x["reads"], "recovery": "" if _x.get("recovery") is None else f"{100 * _x['recovery']:.0f}%",
                    "call": _symbol.get(_x["call"], _x["call"]),
                } for _x in _r["spike_results"]], selection=None, label="Spike-ins: what was added, and what was found"))
            if _r.get("lod"):
                _items.append(mo.ui.table([{
                    "organism": _x["organism"], "unit": _x["unit"], "lowest level detected": _x.get("lowest_detected"),
                    "limit of detection": _x.get("lod"), "levels": len(_x.get("levels") or []),
                } for _x in _r["lod"]], selection=None, label="Limit of detection, per spiked organism (with replicates at several levels)"))
            _items.append(mo.md(f"Written: `{_kw['output_html']}` (the report to keep), `{_kw['output_tsv']}`, `{_kw['output_json']}`."))
            try:
                with open(_kw["output_html"], encoding="utf-8") as _fh:
                    _items.append(mo.accordion({"The engine's full report (HTML)": mo.iframe(_fh.read(), height="900px")}))
            except OSError:
                pass
            _items.append(code_panel(f"r = {python_call('contaminant_screen', _samples, _v['database'], **_kw)}", "print(r['text'])"))
            _panel = mo.vstack(_items)
        else:
            mo.stop(not chosen(_v["table"]), warn("Select the table of the references (name, kind, file)."))
            import csv as _csv_refs
            with open(_v["table"], encoding="utf-8", newline="") as _fh:
                _first = _fh.readline()
            with open(_v["table"], encoding="utf-8", newline="") as _fh:
                _rows = [_x for _x in _csv_refs.DictReader(_fh, delimiter="\t" if "\t" in _first else ",")]
            _rows = [{(_k or "").strip().lower(): (_x or "").strip() for _k, _x in _row.items()} for _row in _rows]
            mo.stop(not _rows or not {"kind", "file"} <= set(_rows[0]),
                    warn("The table needs a header with the columns name, kind and file (one genome per line)."))
            _base = os.path.dirname(_v["table"])
            _files = [[_f.strip() if os.path.isabs(_f.strip()) else (f"{_base}/{_f.strip()}" if _base else _f.strip())
                       for _f in _row["file"].split(";") if _f.strip()] for _row in _rows]
            _kinds = [_row["kind"] for _row in _rows]
            _names = [_row.get("name", "") or None for _row in _rows]
            _out = (_v["output"] or "").strip() or beside(_v["table"], split_name(_v["table"])[0] + ".mmgscreen")
            _kw = {"kinds": _kinds}
            if any(_names):
                _kw["names"] = _names
            if int(_v["memory_mb"] or 1024) != 1024:
                _kw["memory_mb"] = int(_v["memory_mb"])
            _refs_arg = [_f[0] if len(_f) == 1 else _f for _f in _files]
            with mo.status.spinner(f"Indexing {len(_files)} reference genome(s)..."):
                _r = await asyncio.to_thread(mmg.build_screen_database, _refs_arg, _out, **_kw)
            _rows_f = [
                ("Database", f"`{_r['path']}` ({_size(_r['file_bytes'])}, SHA-256 `{_r['sha256'][:16]}…`, {_r['seconds']:.1f} s)"),
                ("Agents and controls", f"`{_r['target_minimizers']:,}` minimizers (k = {_r['k']}), with their genome windows"),
                ("Host", f"`{_r['host_minimizers']:,}` minimizers, in {_r.get('host_passes', 1)} pass(es)"),
            ]
            _items = [mo.callout(mo.md("### 🧱 Screening database\n" + facts(_rows_f)), kind="success")]
            if _r.get("warnings"):
                _items.append(mo.callout(mo.md("\n".join(f"* {_w}" for _w in _r["warnings"])), kind="warn"))
            _items.append(mo.ui.table([{
                "reference": _x["name"], "kind": _x["kind"], "bases": f"{_x['bases']:,}", "sequences": _x["sequences"],
                "minimizers": f"{_x['minimizers']:,}", "genome windows": _x.get("windows", ""), "SHA-256": (_x.get("sha256") or [""])[0][:12] + "…"
                if isinstance(_x.get("sha256"), list) else str(_x.get("sha256", ""))[:12] + "…",
            } for _x in _r["references"]], selection=None, label="The references"))
            _items.append(mo.md("*Next: the task **Screen samples for adventitious agents** with this database (click **Refresh the file list**).*"))
            _items.append(code_panel(f"r = {python_call('build_screen_database', _refs_arg, _out, **_kw)}"))
            _panel = mo.vstack(_items)
    except Exception as _err:
        _panel = error_panel(_err)
    _panel
    return


@app.cell
def ai_helpers(error_panel, mo, os, warn):
    # What the AI screens share. The AI layer comes with the engine (an older
    # engine has none); the Anthropic SDK is needed only when a question is sent.
    try:
        from multimerge_genomics import ai as ai_lib
        from multimerge_genomics.ai import assistant as _assistant
        from multimerge_genomics.ai import knowledge as _knowledge
        from multimerge_genomics.ai import _catalog
        AI_MODELS = dict(_assistant.MODELS)
        AI_DISCLAIMER = _knowledge.DISCLAIMER
        AI_HUMAN_PAGES = set(_catalog._HUMAN_PAGES)
        AI_SCREENS = dict(_catalog.SCREEN_BY_PAGE)
    except ImportError:
        ai_lib, AI_MODELS, AI_DISCLAIMER, AI_HUMAN_PAGES, AI_SCREENS = None, {"Claude Opus 5.5 (default)": "claude-opus-5-5"}, "", set(), {}
    AI_EFFORTS = {"Low: the fastest and cheapest": "low", "Medium (default)": "medium", "High: more careful": "high",
                  "Extra high": "xhigh", "Max: the most careful and the slowest": "max"}
    AI_LEVELS = {"Student: clear words and the idea behind each step": "student",
                 "Researcher: concise, with the method's assumptions and limits": "researcher"}
    AI_SENT = ("**What is sent to Anthropic** (the company behind Claude): the name of the tool, its parameters with "
               "file names only (no folders), a short summary of the result (long lists cut, sequences replaced by their "
               "length) and your question; on screen 29, also the names of the files in the folder you give. Never a "
               "file, a genome or a read. The API is paid by the owner of the key.")

    def ai_blocked(settings, human=False):
        """Why the AI cannot be asked now, as a panel; None when it can."""
        if ai_lib is None:
            return warn("This engine has no AI layer: install a newer `multimerge-genomics`.")
        if not settings["consent"]:
            return warn("The AI is off. On screen 29, task *Settings and the session record*, agree to send result "
                        "summaries to Anthropic and save.")
        if human and not settings["human"]:
            return warn("This is a human-genome result: sending its summary needs the second agreement on screen 29 "
                        "(*human-genome results*).")
        return None

    def ai_assistant(settings):
        return ai_lib.Assistant(api_key=settings["key"] or None, model=settings["model"], effort=settings["effort"])

    def ai_is_human(call):
        return call.page in AI_HUMAN_PAGES or any(t["name"] == call.name and t["human"] for t in ai_lib.catalog())

    def ai_streamer(head):
        """An on_text callback that shows the answer as it grows, a few times a second."""
        import time
        last = [0.0]

        async def show(text):
            now = time.monotonic()
            if now - last[0] >= 0.25:
                last[0] = now
                mo.output.replace(mo.vstack([x for x in (head, mo.md(text + " ▍")) if x is not None]))
        return show

    def ai_footer(ans):
        u = ans.usage or {}
        bits = [f"model `{ans.model or '?'}`"]
        if ans.fallback_used:
            bits.append("answered by the fallback model")
        if u:
            bits.append(f"{u.get('input_tokens', 0):,} tokens in (+{u.get('cache_read_input_tokens', 0):,} read from the cache), "
                        f"{u.get('output_tokens', 0):,} out")
        bits.append(f"{ans.seconds:.1f} s")
        return mo.md("<small>" + " · ".join(bits) + ". The AI can be wrong: check what it says against the result and "
                     "the README.</small>")

    def ai_answer_panel(head, ans, human=False):
        items = [head] if head is not None else []
        if ans.refused:
            items.append(warn("The AI declined to answer this" + (f" ({ans.refusal_category})" if ans.refusal_category else "")
                              + ". Rephrase the question: for instance, ask about the method rather than about a person."))
            if ans.text:
                items.append(mo.md(ans.text))
        else:
            items.append(mo.md(ans.text or "*(an empty answer)*"))
        if ans.stop_reason == "max_tokens":
            items.append(warn("The answer was cut at its length limit."))
        if human and AI_DISCLAIMER:
            items.append(mo.callout(mo.md("**⚠️ " + AI_DISCLAIMER + "**"), kind="danger"))
        items.append(ai_footer(ans))
        return mo.vstack(items)

    def ai_error_text(err):
        name = type(err).__name__
        known = {
            "AuthenticationError": "the API key was refused: check it on screen 29 (or ANTHROPIC_API_KEY).",
            "PermissionDeniedError": "this API key may not use this model: pick another on screen 29.",
            "RateLimitError": "too many requests, or no credit left: wait a moment, or check the account at console.anthropic.com.",
            "APIConnectionError": "the Claude API could not be reached (no internet, or a proxy or an antivirus in the way).",
            "APITimeoutError": "the Claude API took too long to answer: try again, or pick a lower effort on screen 29.",
            "NotFoundError": "the model was not found: pick another on screen 29.",
        }
        return known.get(name, str(err) if name in ("AIUnavailable", "ConsentRequired") else f"{name}: {err}")

    def ai_error(err):
        if type(err).__name__ in ("AIUnavailable", "ConsentRequired"):
            return warn(str(err))
        return error_panel(ai_error_text(err))

    def folder_files(folder, limit=200):
        """The file names under a folder (relative, three levels down): what the AI sees of it."""
        out = []
        root = os.path.abspath(folder or ".")
        for base, dirs, names in os.walk(root):
            dirs[:] = sorted(d for d in dirs if not d.startswith((".", "_")) and d not in ("target", "__pycache__", "node_modules"))
            if os.path.relpath(base, root).count(os.sep) >= 2:
                dirs[:] = []
            out += [os.path.relpath(os.path.join(base, n), root).replace(os.sep, "/") for n in sorted(names)]
            if len(out) >= limit:
                break
        return out[:limit]

    return (AI_DISCLAIMER, AI_EFFORTS, AI_LEVELS, AI_MODELS, AI_SCREENS, AI_SENT, ai_answer_panel, ai_assistant, ai_blocked,
            ai_error, ai_error_text, ai_is_human, ai_lib, ai_streamer, folder_files)


@app.cell
def ai_page_ui(get_page, mo, set_page):
    mo.stop(get_page() != "ai")
    btn_back_ai = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_ai = mo.ui.button(label="🔄 Refresh")
    ai_task = mo.ui.radio(
        options=["Describe your problem", "Let the AI run a task", "Settings and the session record"],
        value="Describe your problem", inline=True, label="**What to do:**",
    )
    _ui = mo.vstack([
        btn_back_ai,
        mo.md("<br>\n\n# 🤖 29. AI Assistant"),
        mo.callout(mo.md(
            "An optional assistant built on Claude, Anthropic's AI. **Describe your problem** in your own words: it "
            "explains the problem in detail, lays out a path through MMGEN's screens and writes the Python script that "
            "follows it. **Let the AI run a task**: it runs MMGEN's tools on a folder for you, round by round, and nothing "
            "runs until you approve it. At the bottom of the other screens, **Ask the AI about this result** explains the "
            "last result of the screen. It needs an Anthropic API key (the API is paid) and your agreement to send result "
            "summaries: see *Settings and the session record*, which also saves everything you ran as a Python script."),
            kind="info"),
        mo.callout(mo.md(
            "**The AI can be wrong.** It is a teaching aid: check what it says against the results, the README and a "
            "teacher or a colleague. On human genomes it never replaces a physician or a genetic counselor."), kind="warn"),
        mo.hstack([ai_task, btn_refresh_ai], justify="space-between", align="center"),
    ])
    _ui
    return ai_task, btn_back_ai, btn_refresh_ai


@app.cell
def ai_page_status(AI_MODELS, ai_lib, ai_settings, get_ai_rev, get_page, label_of, mo, os):
    mo.stop(get_page() != "ai")
    get_ai_rev()  # redrawn when the settings are saved
    if ai_lib is None:
        _text = "⚪ **AI: not available** · this engine has no AI layer."
    else:
        try:
            import anthropic as _sdk  # noqa: F401
            _has_sdk = True
        except ImportError:
            _has_sdk = False
        _s = ai_settings
        _key = ("typed in this session (kept in memory only)" if _s["key"]
                else "from ANTHROPIC_API_KEY" if os.environ.get("ANTHROPIC_API_KEY") else "none")
        _on = _s["consent"] and _has_sdk and _key != "none"
        _text = (f"{'🟢' if _on else '⚪'} **AI: {'on' if _on else 'off'}** · key: {_key} · model: "
                 f"{label_of(AI_MODELS, _s['model'], _s['model'])} · effort: {_s['effort']} · human-genome results: "
                 f"{'may be sent' if _s['human'] else 'not sent'}")
        if not _has_sdk:
            _text += ' · the Anthropic SDK is not installed: `pip install "multimerge-genomics[ai]"` (or `pip install anthropic`)'
    mo.md(_text)
    return


@app.cell
def ai_page_example(ai_task, example_picker, get_page, mo):
    mo.stop(get_page() != "ai")
    ai_example = example_picker("ai", ai_task.value)
    ai_example if ai_task.value != "Settings and the session record" else None
    return (ai_example,)


@app.cell
def ai_page_form(AI_EFFORTS, AI_LEVELS, AI_MODELS, AI_SENT, ai_example, ai_settings, ai_task, example_for, get_page, label_of, make_form, mo):
    mo.stop(get_page() != "ai")
    _ex, _about = example_for("ai", ai_task.value, ai_example.value)
    _s = ai_settings
    _note = None
    if ai_task.value == "Describe your problem":
        form_ai = make_form({
            "problem": mo.ui.text_area(value=_ex.get("problem", ""), rows=7, full_width=True,
                                       label="Your problem, in your own words (what you have, what you want to find out):",
                                       placeholder="e.g. I have the variants of a family in one VCF and want to know whether the child is really theirs..."),
            "folder": mo.ui.text(value=_ex.get("folder", "."), full_width=True, label="The folder with your files (the AI sees their names only):"),
            "level": mo.ui.dropdown(options=AI_LEVELS, value=label_of(AI_LEVELS, _ex.get("level", _s["level"]), next(iter(AI_LEVELS))),
                                    label="Explain it for:", full_width=True),
        }, "🤖 Ask the AI")
    elif ai_task.value == "Let the AI run a task":
        _note = mo.md("The agent proposes calls to MMGEN's tools; **you approve each round before it runs**. It reads and "
                      "writes only inside the working folder, and every call it runs goes into the session record.")
        form_ai = make_form({
            "task": mo.ui.text_area(value=_ex.get("task", ""), rows=5, full_width=True, label="The task:",
                                    placeholder="e.g. Validate every VCF in the folder and write sorted, indexed copies of the valid ones."),
            "folder": mo.ui.text(value=_ex.get("folder", "."), full_width=True, label="Working folder (the agent reads and writes only inside it):"),
            "human": mo.ui.checkbox(value=False, label="Let it use the human-genome tools (screens 16-21 and 28)"),
            "rounds": mo.ui.number(start=1, stop=20, step=1, value=_ex.get("rounds", 8), label="At most this many rounds:"),
        }, "🤖 Give the task to the AI")
    else:
        _note = mo.md(AI_SENT)
        form_ai = make_form({
            "key": mo.ui.text(kind="password", value=_s["key"], full_width=True, placeholder="sk-ant-...",
                              label="Anthropic API key (blank: ANTHROPIC_API_KEY of the environment). Kept in memory, never written:"),
            "model": mo.ui.dropdown(options=AI_MODELS, value=label_of(AI_MODELS, _s["model"], next(iter(AI_MODELS))), label="Model:", full_width=True),
            "effort": mo.ui.dropdown(options=AI_EFFORTS, value=label_of(AI_EFFORTS, _s["effort"], "Medium (default)"),
                                     label="Effort (more is slower and costs more):", full_width=True),
            "level": mo.ui.dropdown(options=AI_LEVELS, value=label_of(AI_LEVELS, _s["level"], next(iter(AI_LEVELS))),
                                    label="Explain for (the default of the other screens):", full_width=True),
            "consent": mo.ui.checkbox(value=_s["consent"], label="I agree to send result summaries and my questions to Anthropic to get answers"),
            "human": mo.ui.checkbox(value=_s["human"], label="Also summaries of human-genome results (screens 16-21 and 28), only of people who agreed to it"),
        }, "💾 Save the settings")
    mo.vstack([_x for _x in (_about, _note, form_ai) if _x is not None])
    return (form_ai,)


@app.cell
async def ai_page_backend(ai_answer_panel, ai_assistant, ai_blocked, ai_error, ai_lib, ai_settings, ai_streamer, ai_task, asyncio, folder_files, form_ai, get_page, mmg, mo, os, set_agent, set_ai_rev, warn):
    mo.stop(get_page() != "ai")
    _task = ai_task.value
    mo.stop(form_ai.value is None, None if _task == "Settings and the session record" else mo.md("⏳ *Fill in the form above and submit it.*"))
    _v = form_ai.value
    if _task == "Settings and the session record":
        ai_settings.update(key=(_v["key"] or "").strip(), model=_v["model"], effort=_v["effort"], level=_v["level"],
                           consent=bool(_v["consent"]), human=bool(_v["consent"] and _v["human"]))
        set_ai_rev(lambda n: n + 1)
        _panel = mo.callout(mo.md("✅ Saved for this session. " + (
            "The AI is on: the other screens now offer *Ask the AI about this result* at the bottom." if _v["consent"]
            else "The AI stays off until you agree to send result summaries.")), kind="success")
    elif _task == "Describe your problem":
        mo.stop(not (_v["problem"] or "").strip(), warn("Type your problem first."))
        mo.stop(not os.path.isdir(_v["folder"] or "."), warn(f"The folder `{_v['folder']}` does not exist."))
        _blocked = ai_blocked(ai_settings)
        mo.stop(_blocked is not None, _blocked)
        _files = folder_files(_v["folder"] or ".")
        _head = mo.md(f"**🤖 The AI's answer** (with the names of {len(_files)} file(s) of `{_v['folder'] or '.'}`):")
        try:
            _ans = await ai_assistant(ai_settings).plan_stream(_v["problem"], _files, _v["level"], on_text=ai_streamer(_head))
            _items = [ai_answer_panel(_head, _ans)]
            for _w in _ans.warnings:
                _items.append(warn("Check the script before running it: " + _w))
            _blocks = ai_lib.assistant.python_blocks(_ans.text)
            _downloads = [mo.download(data=_ans.text.encode("utf-8"), filename="mmgen_plan.md", mimetype="text/markdown",
                                      label="⬇️ The answer (Markdown)")]
            if _blocks:
                _downloads.insert(0, mo.download(data=_blocks[-1].encode("utf-8"), filename="mmgen_plan.py",
                                                 mimetype="text/x-python", label="⬇️ The script"))
            _items.append(mo.hstack(_downloads, justify="start"))
            _panel = mo.vstack(_items)
        except Exception as _err:
            _panel = ai_error(_err)
    else:
        mo.stop(not (_v["task"] or "").strip(), warn("Type the task first."))
        mo.stop(not os.path.isdir(_v["folder"] or "."), warn(f"The working folder `{_v['folder']}` does not exist."))
        _blocked = ai_blocked(ai_settings, human=bool(_v["human"]))
        mo.stop(_blocked is not None, _blocked)
        try:
            _agent = ai_lib.Agent(_v["task"], _v["folder"] or ".", assistant=ai_assistant(ai_settings),
                                  allow_human=bool(_v["human"]), max_rounds=int(_v["rounds"]), engine=mmg)
            with mo.status.spinner("The AI is reading the task and the folder..."):
                await asyncio.to_thread(_agent.step)
            set_agent((_agent, 0))
            _panel = mo.callout(mo.md("🤖 The task went to the AI. Its plan and the calls it proposes are below: "
                                      "**nothing runs until you approve each round.**"), kind="info")
        except Exception as _err:
            _panel = ai_error(_err)
    _panel
    return


@app.cell
def ai_page_agent_view(ai_task, error_panel, get_agent, get_page, mo, warn):
    mo.stop(get_page() != "ai" or ai_task.value != "Let the AI run a task")
    _held = get_agent()
    mo.stop(_held is None)
    _agent = _held[0]
    btn_agent_run = mo.ui.run_button(kind="success", label="✅ Run these calls" if _agent.pending else "🔁 Ask the AI again")
    btn_agent_decline = mo.ui.run_button(kind="danger", label="✋ Decline them")
    _items = [mo.md(f"### 🤖 The agent at work\n**Task:** {_agent.task}  \n**Working folder:** `{_agent.root}`")]
    for _i, _r in enumerate(_agent.rounds, 1):
        _parts = [f"**Round {_i}**"]
        if _r.text.strip():
            _parts.append(_r.text.strip())
        if _r.proposals:
            _code = "\n".join(
                f"mmg.{_p.name}(" + ", ".join(f"{_k}={_x!r}" for _k, _x in _p.input.items()) + ")"
                + (f"  # ⛔ will not run: {_p.problem}" if _p.problem else "") for _p in _r.proposals)
            _parts.append("Proposed calls:\n```python\n" + _code + "\n```")
        _items.append(mo.callout(mo.md("\n\n".join(_parts)), kind="neutral"))
        if _r.results:
            _items.append(mo.accordion({
                ("✅ " if _ok else "❌ ") + f"{_name}: what went back to the AI": mo.plain_text(_text[:6000])
                for _name, _ok, _text in _r.results}))
    if getattr(_agent, "ui_error", None):
        _items.append(error_panel(_agent.ui_error))
    if _agent.done:
        _items.append(mo.callout(mo.md("### The agent's report\n\n" + (_agent.final_text or "")),
                                 kind="warn" if _agent.refused else "success"))
        _items.append(mo.download(data=_agent.transcript().encode("utf-8"), filename="mmgen_agent_transcript.md",
                                  mimetype="text/markdown", label="⬇️ The whole run (Markdown)"))
        _items.append(mo.md(f"<small>{_agent.usage['input_tokens']:,} tokens in, {_agent.usage['output_tokens']:,} out, "
                            f"over {len(_agent.rounds)} round(s). The AI can be wrong: check its report against the results.</small>"))
    elif _agent.pending:
        _items.append(mo.md("**Approve this round?** The calls run on this computer, in the working folder; a summary of "
                            "what they return goes back to the AI."))
        _items.append(mo.hstack([btn_agent_run, btn_agent_decline], justify="start"))
    else:
        _items.append(warn("The last request to the AI did not go through."))
        _items.append(btn_agent_run)
    mo.vstack(_items)
    return btn_agent_decline, btn_agent_run


@app.cell
async def ai_page_agent_act(ai_error_text, asyncio, btn_agent_decline, btn_agent_run, get_agent, mo, set_agent):
    mo.stop(not (btn_agent_run.value or btn_agent_decline.value))
    _agent, _n = get_agent()
    _agent.ui_error = None
    try:
        if _agent.pending:
            if btn_agent_run.value:
                with mo.status.spinner("Running the approved calls..."):
                    await asyncio.to_thread(_agent.run_pending)
            else:
                _agent.decline()
        if not _agent.done:
            with mo.status.spinner("The AI is looking at the results..."):
                await asyncio.to_thread(_agent.step)
    except Exception as _err:
        # kept on the agent: this cell's own output is cleared when the view is redrawn
        _agent.ui_error = ai_error_text(_err)
    set_agent((_agent, _n + 1))
    return


@app.cell
def ai_page_session(AI_SCREENS, ai_task, btn_refresh_ai, get_page, mmg, mo):
    mo.stop(get_page() != "ai" or ai_task.value != "Settings and the session record")
    btn_refresh_ai  # a click redraws the record
    mo.stop(not hasattr(mmg, "script"))
    _rows = []
    for _c in reversed(list(mmg.calls)):
        _screen = AI_SCREENS.get(_c.page)
        _rows.append({"time": _c.when, "screen": f"{_screen[0]}. {_screen[2]}" if _screen else "—", "tool": _c.name,
                      "seconds": round(_c.seconds, 2), "status": "failed: " + _c.error[:120] if _c.error else "ok"})
    _items = [mo.md("### 📜 The session record\nEvery tool run in this session, on any screen and by the agent, newest "
                    "first. Saved as a Python script it does the same again: the record of what you did, for your "
                    "notebook or the methods section of a report. Press **Refresh** to see the latest calls.")]
    if _rows:
        _items += [mo.ui.table(_rows, selection=None, label=f"{len(_rows)} call(s)"),
                   mo.download(data=mmg.script().encode("utf-8"), filename="mmgen_session.py", mimetype="text/x-python",
                               label="⬇️ The session as a Python script")]
    else:
        _items.append(mo.md("*Nothing has been run in this session yet.*"))
    mo.vstack(_items)
    return


@app.cell
def ai_panel_ui(AI_LEVELS, AI_SENT, ai_lib, ai_settings, get_ai_rev, get_page, label_of, make_form, mmg, mo, set_page):
    # "Ask the AI about this result", at the bottom of every screen with results.
    get_ai_rev()  # redrawn when the AI settings change
    mo.stop(get_page() in ("home", "ai", "guide") or ai_lib is None or not hasattr(mmg, "last"))
    btn_ai_setup = mo.ui.button(label="🤖 Open the AI settings (screen 29)", on_change=lambda _: set_page("ai"))
    form_ai_ask = make_form({
        "question": mo.ui.text_area(rows=3, full_width=True, label="Your question (blank: explain the result):",
                                    placeholder="e.g. Why were so many lines skipped? Is this number good or bad?"),
        "level": mo.ui.dropdown(options=AI_LEVELS, value=label_of(AI_LEVELS, ai_settings["level"], next(iter(AI_LEVELS))),
                                label="Explain it for:", full_width=True),
    }, "🤖 Ask the AI")
    _body = [mo.md("The AI explains **the last result of this screen**: what the tool did, what the numbers mean, whether "
                   "they look as expected, what could be wrong and what to do next.\n\n" + AI_SENT)]
    if ai_settings["consent"]:
        _body.append(form_ai_ask)
    else:
        _body += [mo.md("*The AI is off: turn it on in the settings of screen 29 (it needs an Anthropic API key and your agreement).*"),
                  btn_ai_setup]
    mo.accordion({"🤖 Ask the AI about this result (optional)": mo.vstack(_body)})
    return btn_ai_setup, form_ai_ask


@app.cell
async def ai_panel_answer(ai_answer_panel, ai_assistant, ai_blocked, ai_error, ai_is_human, ai_settings, ai_streamer, form_ai_ask, get_page, mmg, mo, warn):
    mo.stop(form_ai_ask.value is None)
    _call = mmg.last(get_page())
    mo.stop(_call is None, warn("Run a task on this screen first: the AI explains the last result of the screen."))
    _human = ai_is_human(_call)
    _blocked = ai_blocked(ai_settings, human=_human)
    mo.stop(_blocked is not None, _blocked)
    _q = (form_ai_ask.value["question"] or "").strip()
    _head = mo.md(f"**🤖 About the result of `mmg.{_call.name}`** (run at {_call.when})" + (f"  \n**Your question:** {_q}" if _q else ""))
    try:
        _ans = await ai_assistant(ai_settings).explain_stream(_call.python(), _call.result, _q, form_ai_ask.value["level"],
                                                               on_text=ai_streamer(_head))
        _panel = ai_answer_panel(_head, _ans, human=_human)
    except Exception as _err:
        _panel = ai_error(_err)
    _panel
    return


if __name__ == "__main__":
    app.run()
