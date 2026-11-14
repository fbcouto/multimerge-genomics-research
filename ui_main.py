import marimo

__generated_with = "0.23.11"
# app_title is the browser tab title.
app = marimo.App(width="full", app_title="MMGEN Suite")


@app.cell
def imports():
    import marimo as mo
    import polars as pl
    import numpy as np
    import glob
    import os
    import asyncio
    import multimerge_genomics as mmg
    return asyncio, glob, mmg, mo, np, os, pl


@app.cell
def global_state(mo):
    # allow_self_loops: the cell that calls the setter is re-run too. Without it
    # the menu (which calls set_page) never left the page when a screen opened,
    # and the screen you left through "Back to Menu" stayed on the page.
    get_page, set_page = mo.state("home", allow_self_loops=True)

    # Message of the "Generate simulated patient" button (screen 4). It lives
    # here and not in the screen's own cell: a state created in the very cell
    # that shows it is reborn empty on every run.
    get_patient_msg, set_patient_msg = mo.state("", allow_self_loops=True)
    return get_page, get_patient_msg, set_page, set_patient_msg


@app.cell
def file_helpers(mo, os):
    # Shared by screens 9-14: the files of the folder each screen offers, the
    # default output names, a look at the first records of a result, and the
    # Python call that does the same thing (shown under every result, so each
    # screen is also a recipe for a script).
    NO_FILE = "(no matching file in this folder)"
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
    }
    COMPRESSIONS = (".gz", ".bgz", ".lz4", ".parquet")
    FORMAT_EXTENSIONS = ("vcf", "bed", "gff", "gff3", "gtf", "sam", "bam", "fa", "fasta", "fna", "fq", "fastq", "tsv", "txt")

    def files_with(*kinds):
        """The files of the current folder with an extension of `kinds`."""
        ends = tuple(s for k in kinds for s in FILE_KINDS[k])
        found = sorted(n for n in os.listdir(".") if n.lower().endswith(ends) and os.path.isfile(n))
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

    def output_name(typed, source, tag, ext=None):
        """The name typed in the form, or <source>.<tag><extension> when it was left blank."""
        typed = (typed or "").strip()
        if typed:
            return typed
        stem, source_ext = split_name(source)
        return f"{stem}.{tag}{source_ext if ext is None else ext}"

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
        items = [mo.md("**Abundance spectrum**"), mo.Html(svg)]
        if beyond:
            items.append(mo.md(f"*{beyond:,} distinct k-mer(s) seen more than {x_max:,} times are left out of the chart; the table below has them.*"))
        return mo.vstack(items)

    return (
        chosen,
        clashes,
        code_panel,
        error_panel,
        facts,
        files_with,
        make_form,
        output_name,
        preview_panel,
        python_call,
        real_files,
        size_of,
        spectrum_chart,
        split_name,
        warn,
    )


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
        <p style="color: gray; font-size: 0.85em; margin-top: 40px;">Engine {_engine["version"]} · running on: {_where} · {_engine["threads"]} threads</p>
    </div>
    """)

    menu_screen
    return (
        btn_menu_annotation,
        btn_menu_assembly,
        btn_menu_bloom,
        btn_menu_cluster,
        btn_menu_convert,
        btn_menu_fastq,
        btn_menu_kmercount,
        btn_menu_qc,
        btn_menu_scaffold,
        btn_menu_sort,
        btn_menu_validate,
        btn_menu_variant,
        btn_menu_vcf,
        btn_menu_vcftools,
        menu_screen,
    )


@app.cell
def qc_page_ui(get_page, glob, mo, set_page):
    mo.stop(get_page() != "qc")

    btn_back_qc = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))

    _files = glob.glob("*.fastq") + glob.glob("*.fastq.gz") + glob.glob("*.fasta") + glob.glob("*.fasta.gz") + glob.glob("*.fa")
    if not _files: _files = ["No file found"]

    form_qc = mo.ui.form(
        element=mo.ui.dictionary({
            "files": mo.ui.multiselect(options=_files, label="🧬 DNA files (several allowed, .gz too):"),
            "kmer_size": mo.ui.slider(start=21, stop=63, step=2, value=55, label="K-mer size:"),
            "min_quality": mo.ui.slider(start=0, stop=40, step=1, value=20, label="Minimum mean quality (FASTQ only):")
        }),
        bordered=True,
        submit_button_label="🚀 Start batch ingestion"
    )

    ui_qc = mo.vstack([
        btn_back_qc,
        mo.md("<br>\n\n# 🧬 1. Data Ingestion & QC (Batch / GZIP)"),
        mo.callout(mo.md("The engine reads **several files at once** (for example paired-end R1 and R2) and decompresses `.gz` files on the fly, straight into memory, so nothing extra is written to disk."), kind="info"),
        form_qc
    ])

    ui_qc
    return btn_back_qc, form_qc, ui_qc


@app.cell
async def qc_page_backend(asyncio, form_qc, get_page, mmg, mo, pl):
    mo.stop(get_page() != "qc")
    mo.stop(form_qc.value is None, mo.md("⏳ *Waiting for the settings...*"))

    _paths = form_qc.value.get("files", [])
    mo.stop(not _paths or _paths[0].startswith("No file"), mo.callout(mo.md("⚠️ Select at least one file."), kind="warn"))

    _k = form_qc.value["kmer_size"]
    _min_quality = form_qc.value["min_quality"]

    _is_fasta = any(c.endswith(".fasta") or c.endswith(".fa") or c.endswith(".fasta.gz") for c in _paths)

    try:
        with mo.status.spinner(f"Processing {len(_paths)} file(s) in parallel (K={_k})..."):
            await asyncio.sleep(0.1)

            if _is_fasta:
                _kmer_bytes = mmg.read_fasta_batch(_paths, _k)
            else:
                _kmer_bytes = mmg.read_fastq_batch(_paths, _min_quality, _k)

        with mo.status.spinner("Saving the k-mers (Arrow file)..."):
            await asyncio.sleep(0.1)
            _df = pl.DataFrame({"K-mers": [_kmer_bytes]})

            if len(_paths) > 1:
                _base_name = "Multi_Batch"
            else:
                _base_name = _paths[0].split(".")[0]

            _out_name = f"{_base_name}_K{_k}.arrow"
            _df.write_ipc(_out_name)

        _panel = mo.vstack([
            mo.callout(mo.md(f"""
            ### ✅ Batch processed successfully!
            * **Files read:** `{len(_paths)}`
            * **K-mer size:** `{_k}`
            * **Distinct k-mers extracted:** `{len(_kmer_bytes) // 16:,}`
            * **Arrow file written:** `{_out_name}`
            """), kind="success")
        ])

    except Exception as e:
        _panel = mo.callout(mo.md(f"**🔥 ERROR:** `{str(e)}`"), kind="danger")

    _panel
    return


@app.cell
def cluster_page_ui(get_page, glob, mo, set_page):
    mo.stop(get_page() != "cluster")

    # A global variable (no leading "_"), like on the other screens: marimo only
    # keeps a component alive when it is bound to a global variable. As
    # `_btn_back`, this button did nothing when clicked.
    btn_back_cluster = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))

    _files = glob.glob("*.arrow") + glob.glob("*.txt") + glob.glob("*.fasta")
    if not _files: _files = ["No file found"]

    form_samples = mo.ui.form(
        element=mo.ui.dictionary({
            "sample_a": mo.ui.dropdown(options=_files, label="📂 Sample A (.txt, one 31-mer per line):"),
            "sample_b": mo.ui.dropdown(options=_files, label="📂 Sample B (.txt, one 31-mer per line):")
        }), bordered=True, submit_button_label="🚀 1. Load and merge the k-mers"
    )

    _ui = mo.vstack([btn_back_cluster, mo.md("<br>\n\n# 🧬 Taxonomic Clustering"), form_samples])
    _ui
    return btn_back_cluster, form_samples


@app.cell
async def cluster_page_backend(asyncio, form_samples, get_page, mmg, mo, np, pl):
    mo.stop(get_page() != "cluster")
    mo.stop(form_samples.value is None, mo.md("🧬 *Waiting for the two samples...*"))

    _p_a, _p_b = form_samples.value["sample_a"], form_samples.value["sample_b"]
    mo.stop(not _p_a or not _p_b or _p_a.startswith("No file"), mo.callout(mo.md("⚠️ Select both samples."), kind="warn"))
    # An .arrow next to a .txt went to the text reader, which read the binary as lines.
    mo.stop(
        _p_a.endswith(".arrow") != _p_b.endswith(".arrow"),
        mo.callout(mo.md("⚠️ Both samples must be of the same kind: two `.arrow` files or two text files."), kind="warn")
    )

    with mo.status.spinner("Reading the files, sorting and merging the k-mers..."):
        await asyncio.sleep(0.1)

        if _p_a.endswith(".arrow") and _p_b.endswith(".arrow"):
            _col_a = pl.read_ipc(_p_a).get_column("K-mers")
            _col_b = pl.read_ipc(_p_b).get_column("K-mers")
            # This module works with 64-bit k-mers (k=31), one per row. The
            # .arrow files Module 1 writes today are something else: a single
            # row holding every 128-bit k-mer as binary. On that column
            # .to_numpy().tobytes() returns the bytes of the pointer, and the
            # screen carried on with "2 k-mers" of garbage and no error at all.
            mo.stop(
                not (_col_a.dtype.is_integer() and _col_b.dtype.is_integer()),
                mo.callout(mo.md("⚠️ This module clusters k-mers of 31 bases (64 bits). The `.arrow` files written by Module 1 hold 128-bit k-mers and are not accepted here yet. Use `.txt` files with one 31-base k-mer per line."), kind="warn")
            )
            _bytes_a = _col_a.to_numpy().tobytes()
            _bytes_b = _col_b.to_numpy().tobytes()
            _merged_bytes = mmg.sort_merge_bytes(_bytes_a, _bytes_b)
        else:
            _merged_bytes = mmg.read_and_merge_kmer_files(_p_a, _p_b)

        _merged = np.frombuffer(_merged_bytes, dtype=np.uint64)

    with mo.status.spinner("Computing MinHash signatures on the CPU..."):
        await asyncio.sleep(0.1)
        _signature_bytes = mmg.minhash_signatures(_merged_bytes, 50)
        signature_matrix = np.frombuffer(_signature_bytes, dtype=np.uint64).reshape(-1, 50)

        df_kmers = pl.DataFrame({
            "K-mer (u64)": _merged,
            "DNA_sequence": mmg.decode_kmers(_merged_bytes, 31)
        })

    _panel = mo.md(f"✅ Engine finished. **{df_kmers.height:,}** k-mers ready for clustering.")
    _panel
    return df_kmers, signature_matrix


@app.cell
def cluster_page_hdbscan(df_kmers, get_page, mo):
    mo.stop(get_page() != "cluster")

    form_cluster = mo.ui.form(
        element=mo.ui.dictionary({"sample_percent": mo.ui.slider(start=1, stop=100, step=1, value=5, label="HDBSCAN sample (%):")}),
        bordered=True, submit_button_label="🧬 2. Cluster the k-mers"
    )

    form_cluster
    return (form_cluster,)


@app.cell
async def cluster_page_result(asyncio, df_kmers, form_cluster, get_page, mo, np, pl, signature_matrix):
    from sklearn.cluster import HDBSCAN
    from sklearn.decomposition import PCA

    mo.stop(get_page() != "cluster")
    mo.stop(form_cluster.value is None)

    _fraction = form_cluster.value["sample_percent"] / 100.0
    _total_rows = len(df_kmers)

    # With very little data (a 6 k-mer test, say) use 100% of it.
    if _total_rows < 1000:
        _sample_size = _total_rows
    else:
        _sample_size = int(_total_rows * _fraction)

    if _sample_size < 4:
        mo.stop(True, mo.md("⚠️ **Error:** not enough data. The model needs at least 4 k-mers to compute the 3-D projection (PCA)."))

    if _sample_size > 100000:
        mo.stop(True, mo.md(f"""
        ⚠️ **CPU bottleneck warning:** you are trying to cluster **{_sample_size:,}** k-mers. HDBSCAN in Python (CPU) is not built for that volume and may take hours.
        👉 **What to do:** lower the slider and click Cluster again.
        """))

    with mo.status.spinner(f"Sampling {_sample_size:,} points and computing their geometry..."):
        await asyncio.sleep(0.1)

        _indices = np.arange(_total_rows)
        # Only slice when the sample is smaller than the total.
        if _sample_size < _total_rows:
            _indices = np.random.choice(_indices, _sample_size, replace=False)

        _df_sample = df_kmers[_indices]
        _embeddings = signature_matrix[_indices].astype(np.float64)

        _reduced = PCA(n_components=3).fit_transform(_embeddings)

        # HDBSCAN wants groups of at least 50; adjust for small tests.
        _min_cluster = min(50, max(2, _sample_size // 2))
        _labels = HDBSCAN(min_cluster_size=_min_cluster, n_jobs=-1, copy=False).fit_predict(_reduced)

        _df_clustered = _df_sample.with_columns(pl.Series("Cluster_ID", _labels))
        _total_clusters = len(set(_labels)) - (1 if -1 in _labels else 0)

    _table = mo.ui.table(_df_clustered, selection="single", pagination=True)

    _final_panel = mo.vstack([
        mo.callout(mo.md(f"**Result:** {_total_clusters} cluster(s) found. (K-mers left ungrouped / noise: {sum(_labels == -1)})"), kind="success"),
        _table
    ])

    _final_panel
    return HDBSCAN, PCA


@app.cell
def assembly_page_ui(get_page, glob, mo, set_page):
    mo.stop(get_page() != "assembly")
    btn_back_assembly = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))

    # The .arrow files written by Module 1.
    _files = glob.glob("*.arrow")
    if not _files: _files = ["No file found (.arrow)"]

    form_assembly = mo.ui.form(
        element=mo.ui.dictionary({
            "kmer_file": mo.ui.dropdown(options=_files, label="🧬 K-mer file (.arrow):")
        }),
        bordered=True,
        submit_button_label="🧩 Build the De Bruijn graph and assemble"
    )

    ui_assembly = mo.vstack([
        btn_back_assembly,
        mo.md("<br>\n\n# 🧩 Genome Assembly (De Bruijn)"),
        mo.callout(mo.md("This module takes the k-mers from Module 1, splits each one into a prefix (K-1) and a suffix (K-1), and maps the overlaps in memory, in Rust. The unbranched paths of that graph are the contigs."), kind="info"),
        form_assembly
    ])

    ui_assembly
    return btn_back_assembly, form_assembly, ui_assembly


@app.cell
async def assembly_page_backend(asyncio, form_assembly, get_page, mmg, mo, os, pl):
    mo.stop(get_page() != "assembly")
    mo.stop(form_assembly.value is None, mo.md("🧩 *Waiting for the assembly button...*"))

    _path = form_assembly.value.get("kmer_file")
    mo.stop(not _path or _path.startswith("No file"), mo.callout(mo.md("⚠️ Select a valid .arrow file."), kind="warn"))

    try:
        # K comes from the file name (file_K55.arrow -> 55).
        _k_text = _path.split("_K")[-1].split(".arrow")[0]
        _k_size = int(_k_text)

        with mo.status.spinner(f"1/3 Reading the k-mers (K={_k_size})..."):
            await asyncio.sleep(0.1)
            _df = pl.read_ipc(_path)
            _kmer_bytes = _df.get_column("K-mers").item(0)

        with mo.status.spinner("2/3 Rust: building the graph and walking its paths..."):
            await asyncio.sleep(0.1)
            _nodes, _edges = mmg.build_debruijn_graph(_kmer_bytes, _k_size)
            contig_list = mmg.assemble_debruijn(_kmer_bytes, _k_size)

        with mo.status.spinner("3/3 Saving to disk and preparing the table..."):
            await asyncio.sleep(0.1)

            _out_name = f"MMGEN_Contigs_K{_k_size}.fasta"
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

    except Exception as e:
        _panel = mo.callout(mo.md(f"**🔥 ERROR:** `{str(e)}`"), kind="danger")

    _panel
    return


@app.cell
def variant_page_ui(get_page, get_patient_msg, glob, mo, os, set_page, set_patient_msg):
    mo.stop(get_page() != "variant")

    def _back(_):
        set_patient_msg("")
        set_page("home")

    btn_back_variant = mo.ui.button(label="⬅️ Back to Menu", on_change=_back)

    def generate_patient():
        # The simulated patient is made from a file with exactly this name.
        _base_name = "ecoli.fasta"

        if not os.path.exists(_base_name):
            set_patient_msg(f"❌ Error: the file '{_base_name}' is not in this folder!")
            return

        with open(_base_name, "r") as f: lines = f.readlines()
        seq = "".join([l.strip() for l in lines if not l.startswith(">")])
        seq_list = list(seq)

        if len(seq_list) > 30000:
            seq_list[5000] = 'T' if seq_list[5000] == 'A' else 'A'
            seq_list[15000] = 'G' if seq_list[15000] == 'C' else 'C'
            seq_list[25000] = 'C' if seq_list[25000] == 'G' else 'G'

        mutated = "".join(seq_list)
        with open("mutated_patient.fasta", "w") as f:
            f.write(">Mutated_Patient_Ecoli\n")
            for i in range(0, len(mutated), 80): f.write(mutated[i:i+80] + "\n")

        set_patient_msg("✅ Done! 'mutated_patient.fasta' was written to disk with 3 mutations.")

    btn_generate = mo.ui.button(label="💉 Generate simulated patient", on_click=lambda _: generate_patient())

    _files = glob.glob("*.arrow")
    if not _files: _files = ["No file found (.arrow)"]

    form_variant = mo.ui.form(
        element=mo.ui.dictionary({
            "reference": mo.ui.dropdown(options=_files, label="🧬 1. Healthy reference genome (.arrow):"),
            "patient": mo.ui.dropdown(options=_files, label="🩸 2. Patient sample (.arrow):"),
            "filter_mb": mo.ui.slider(start=1, stop=128, step=1, value=16, label="Bloom filter size (MB):"),
            "num_hashes": mo.ui.slider(start=1, stop=10, step=1, value=5, label="Hashes per k-mer:")
        }),
        bordered=True,
        submit_button_label="🔍 3. Run the search"
    )

    ui_variant = mo.vstack([
        btn_back_variant,
        mo.md("<br>\n\n# ⚠️ Variant Calling (SNPs)"),
        mo.callout(mo.md("To find mutations we compare a patient against a reference. The button below creates a simulated patient from `ecoli.fasta`, with 3 planted mutations; run it through Module 1 to get its `.arrow` file."), kind="info"),
        btn_generate,
        # The button's feedback (an empty one gave "****", which Markdown draws
        # as a horizontal rule).
        mo.md(f"**{get_patient_msg()}**" if get_patient_msg() else ""),
        form_variant
    ])

    ui_variant
    return btn_back_variant, btn_generate, form_variant, generate_patient, ui_variant


@app.cell
async def variant_page_backend(asyncio, form_variant, get_page, mmg, mo, pl):
    mo.stop(get_page() != "variant")
    mo.stop(form_variant.value is None, mo.md("🧪 *Waiting for the settings...*"))

    _ref = form_variant.value.get("reference")
    _pat = form_variant.value.get("patient")

    mo.stop(not _ref or not _pat or _ref.startswith("No file"), mo.callout(mo.md("⚠️ Select both files."), kind="warn"))
    mo.stop(_ref == _pat, mo.callout(mo.md("⚠️ The reference and the patient cannot be the same file!"), kind="warn"))

    _k_text = _ref.split("_K")[-1].split(".arrow")[0]
    _k_size = int(_k_text)
    _mb = form_variant.value["filter_mb"]
    _hashes = form_variant.value["num_hashes"]
    # Where the Bloom filter really runs on this machine.
    _where = "GPU" if mmg.engine_info()["gpu"] else "CPU"

    try:
        with mo.status.spinner("1/3 Loading the k-mers..."):
            await asyncio.sleep(0.1)
            _ref_bytes = pl.read_ipc(_ref).get_column("K-mers").item(0)
            _pat_bytes = pl.read_ipc(_pat).get_column("K-mers").item(0)

        with mo.status.spinner(f"2/3 Building the reference Bloom filter on the {_where}..."):
            await asyncio.sleep(0.1)
            _bloom_bytes = mmg.build_bloom_filter(_ref_bytes, _mb, _hashes)

        with mo.status.spinner(f"3/3 Looking for k-mers found only in the patient ({_where})..."):
            await asyncio.sleep(0.1)
            _mutation_bytes = mmg.extract_mutations(_pat_bytes, _bloom_bytes, _hashes)
            _mutated_seqs = mmg.decode_kmers_128(_mutation_bytes, _k_size)

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
def annotation_page_ui(get_page, mo, set_page):
    mo.stop(get_page() != "annotation")
    btn_back_annotation = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))

    form_annotation = mo.ui.form(
        element=mo.ui.dictionary({
            "dna": mo.ui.text_area(
                label="🧬 Paste the mutated DNA from Module 4 (or a whole contig):",
                # One of the mutated sequences found in the E. coli example.
                value="ATTTACCCTTGCTTTCACATAATCTGCCCGCCGATTTTGCTGCGTTGCGTAAATA"
            )
        }),
        bordered=True,
        submit_button_label="⚙️ Translate to protein"
    )

    ui_annotation = mo.vstack([
        btn_back_annotation,
        mo.md("<br>\n\n# 📘 Functional Annotation (Translation)"),
        mo.callout(mo.md("DNA alone is only an instruction manual. To see whether a mutation matters, we translate the code three letters at a time (codons) and check whether the **protein** changed. Paste a sequence below to try it."), kind="info"),
        form_annotation
    ])

    ui_annotation
    return btn_back_annotation, form_annotation, ui_annotation


@app.cell
async def annotation_page_backend(asyncio, form_annotation, get_page, mmg, mo):
    mo.stop(get_page() != "annotation")
    mo.stop(form_annotation.value is None, mo.md("🧬 *Waiting for the DNA to translate...*"))

    _dna = form_annotation.value.get("dna", "").strip().upper()
    mo.stop(not _dna, mo.callout(mo.md("⚠️ Enter a valid DNA sequence."), kind="warn"))

    try:
        with mo.status.spinner("Translating DNA into protein in Rust (Rayon, multi-core)..."):
            await asyncio.sleep(0.1)
            # Rust takes a list of sequences and returns a list of proteins.
            _results = mmg.translate_proteins([_dna])
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
def bloom_page_ui(get_page, glob, mo, set_page):
    mo.stop(get_page() != "bloom")
    btn_back_bloom = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))

    _files = glob.glob("*.arrow")
    _options = _files if _files else ["No file found (.arrow)"]

    form_bloom = mo.ui.form(
        element=mo.ui.dictionary({
            "kmer_file": mo.ui.dropdown(options=_options, label="🧬 Genome (.arrow):"),
            "filter_mb": mo.ui.slider(start=1, stop=128, step=1, value=16, label="Filter size (MB):"),
            "num_hashes": mo.ui.slider(start=1, stop=10, step=1, value=5, label="Number of hashes:")
        }),
        bordered=True,
        submit_button_label="🔍 Build the Bloom filter"
    )

    _ui = mo.vstack([
        btn_back_bloom,
        mo.md("<br>\n\n# 🔍 Probabilistic Bloom Filters"),
        form_bloom
    ])
    _ui
    return btn_back_bloom, form_bloom


@app.cell
async def bloom_page_backend(asyncio, form_bloom, get_page, mmg, mo, np, pl):
    mo.stop(get_page() != "bloom")
    mo.stop(form_bloom.value is None, mo.md("🔍 *Set the parameters above and click 'Build'...*"))

    _path = form_bloom.value.get("kmer_file")
    mo.stop(not _path or _path.startswith("No file"), mo.md("⚠️ Select a valid .arrow file."))

    _filter_mb = form_bloom.value["filter_mb"]
    _num_hashes = form_bloom.value["num_hashes"]
    # Where the filter really runs on this machine.
    _where = "GPU" if mmg.engine_info()["gpu"] else "CPU"

    with mo.status.spinner(f"🔥 Building the filter on the {_where}..."):
        await asyncio.sleep(0.1)
        _df = pl.read_ipc(_path)
        _kmer_bytes = _df.get_column("K-mers").item(0)
        _total_kmers = len(_kmer_bytes) // 16

        _bloom_bytes = mmg.build_bloom_filter(_kmer_bytes, _filter_mb, _num_hashes)

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
def scaffold_page_ui(get_page, glob, mo, set_page):
    mo.stop(get_page() != "scaffold")
    btn_back_scaffold = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))

    _files = glob.glob("*.fasta") + glob.glob("*.fa")
    if not _files: _files = ["No file found (FASTA)"]

    form_scaffold = mo.ui.form(
        element=mo.ui.dictionary({
            "contigs_fasta": mo.ui.dropdown(options=_files, label="🧩 Contigs file (FASTA):"),
            # With K=55 in the assembly, the contigs overlap by exactly K-1 (54).
            "overlap": mo.ui.slider(start=21, stop=100, step=1, value=54, label="Exact overlap at the ends (bases):")
        }),
        bordered=True,
        submit_button_label="🏗️ Close the genome (scaffolding)"
    )

    ui_scaffold = mo.vstack([
        btn_back_scaffold,
        mo.md("<br>\n\n# 🧬 Genome Scaffolding"),
        mo.callout(mo.md("The De Bruijn graph produced contigs. Scaffolding takes those blocks and stitches them together wherever their ends overlap exactly, trying to rebuild the original chromosome, circular ones included."), kind="info"),
        form_scaffold
    ])

    ui_scaffold
    return btn_back_scaffold, form_scaffold, ui_scaffold


@app.cell
async def scaffold_page_backend(asyncio, form_scaffold, get_page, mmg, mo):
    mo.stop(get_page() != "scaffold")
    mo.stop(form_scaffold.value is None, mo.md("🏗️ *Waiting for a FASTA file...*"))

    _path = form_scaffold.value.get("contigs_fasta")
    _overlap = form_scaffold.value.get("overlap")

    mo.stop(not _path or _path.startswith("No file"), mo.callout(mo.md("⚠️ Select a valid FASTA file produced by Module 3."), kind="warn"))

    try:
        with mo.status.spinner("Reading the contigs..."):
            await asyncio.sleep(0.1)
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
            await asyncio.sleep(0.1)
            scaffolds, circular = mmg.build_scaffolds(contigs, _overlap)

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

    except Exception as e:
        _panel = mo.callout(mo.md(f"**🔥 ERROR:** `{str(e)}`"), kind="danger")

    _panel
    return


@app.cell
def vcf_page_ui(get_page, glob, mo, set_page):
    mo.stop(get_page() != "vcf")
    btn_back_vcf = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))

    _files = sorted(glob.glob("*.vcf") + glob.glob("*.vcf.gz"))
    if not _files:
        _files = ["No file found (put the VCF files in the project folder)"]

    form_vcf = mo.ui.form(
        element=mo.ui.dictionary({
            "files": mo.ui.multiselect(options=_files, label="🗄️ Files (VCF, VCF.GZ):"),
            "sort_by": mo.ui.dropdown(options=["POS", "QUAL"], value="POS", label="📊 Sort by:"),
            "output_name": mo.ui.text(value="Merged_Sorted.vcf", label="💾 Output file name:")
        }),
        bordered=True,
        submit_button_label="🚀 Start the Multimerge engine"
    )

    ui_vcf = mo.vstack([
        btn_back_vcf,
        mo.md("<br>\n\n# 🗄️ 8. VCF Multimerge & Sorting"),
        mo.callout(mo.md("Bulk sorting of genetic variants. Select one or more `.vcf` or `.vcf.gz` files: the engine merges them all into a single output, with a stable sort, for files of any size (what does not fit in memory is sorted in runs on disk)."), kind="info"),
        form_vcf
    ])

    ui_vcf
    return btn_back_vcf, form_vcf, ui_vcf


@app.cell
async def vcf_page_backend(asyncio, form_vcf, get_page, mmg, mo, os):
    mo.stop(get_page() != "vcf")
    mo.stop(form_vcf.value is None, mo.md("⏳ *Waiting for the settings...*"))

    _paths = list(form_vcf.value.get("files", []))
    _sort_by = form_vcf.value.get("sort_by")
    _output = (form_vcf.value.get("output_name") or "").strip() or "Merged_Sorted.vcf"

    mo.stop(not _paths or _paths[0].startswith("No file"), mo.callout(mo.md("⚠️ Select at least one file."), kind="warn"))
    # The result of an earlier run shows up in the list (it is a .vcf in the same
    # folder): reading and writing the same file would destroy the data.
    mo.stop(
        os.path.abspath(_output) in [os.path.abspath(_c) for _c in _paths],
        mo.callout(mo.md("⚠️ The output file cannot be one of the input files. Choose another name."), kind="warn")
    )

    try:
        with mo.status.spinner(f"Merging and sorting {len(_paths)} file(s) in the Multimerge engine..."):
            await asyncio.sleep(0.1)

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
        value="merge: join overlapping or nearby intervals (bedtools merge)", label="Operation:",
    )

    _ui = mo.vstack([
        btn_back_sort,
        mo.md("<br>\n\n# 📏 9. Sorting, Indexing & Intervals"),
        mo.callout(mo.md("Sort and merge VCF, BED, GFF/GTF and SAM/BAM files of any size (what does not fit in memory is sorted in runs on disk); build the `.tbi` / `.bai` index of a sorted file and read only the region you ask for; or run the interval operations of bedtools on BED files, which do not need to be sorted. Results are written to this folder."), kind="info"),
        mo.hstack([sort_task, btn_refresh_sort], justify="space-between", align="center"),
    ])
    _ui
    return btn_back_sort, btn_refresh_sort, interval_op, sort_task


@app.cell
def sort_page_form(btn_refresh_sort, files_with, get_page, interval_op, make_form, mo, real_files, sort_task):
    mo.stop(get_page() != "sort")
    btn_refresh_sort  # a click re-runs this cell, which lists the folder again
    _task = sort_task.value
    _above = []

    if _task == "Sort & merge":
        _orders = {
            "default: the header for VCF and SAM, lexical for BED and GFF": "default",
            "lexical: chr1, chr10, chr2 (bedtools)": "lexical",
            "natural: chr1, chr2, chr10": "natural",
            "header: the ##contig / @SQ lines (bcftools, samtools)": "header",
        }
        for _fai in real_files("genome"):
            _orders[f"the order of {_fai}"] = _fai
        form_sort = make_form({
            "files": mo.ui.multiselect(options=files_with("vcf", "bed", "gff", "sam"), label="Files (several are merged into one output):"),
            "by": mo.ui.dropdown(options={"position": "position", "quality (VCF, highest first)": "quality", "read name (SAM)": "name"},
                                 value="position", label="Sort by:"),
            "order": mo.ui.dropdown(options=_orders, value="default: the header for VCF and SAM, lexical for BED and GFF", label="Chromosome order:"),
            "memory": mo.ui.number(start=0, stop=1_048_576, step=64, value=0, label="Memory budget in MB (0 = half of the free RAM):"),
            "index": mo.ui.checkbox(value=False, label="Also write the index (.tbi for a .gz output, .bai for a .bam output)"),
            "output": mo.ui.text(placeholder="blank: <first file>.sorted.<same extension>", label="Output file:", full_width=True),
        }, "🚀 Sort")
    elif _task == "Index & query a region":
        form_sort = make_form({
            "file": mo.ui.dropdown(options=files_with("indexed"), label="Sorted file (a .gz written as BGZF, or a .bam):"),
            "regions": mo.ui.text(placeholder="chr1:1000000-2000000, chr2", label="Region(s), 1-based, end included (blank: only build the index):", full_width=True),
            "output": mo.ui.text(placeholder="blank: show the records here", label="Save the records to:", full_width=True),
        }, "🔎 Index and query")
    else:
        _op = interval_op.value
        _bed = files_with("bed")
        _fields = {"a": mo.ui.multiselect(options=_bed, label="A: BED file(s):")}
        if _op in ("intersect", "subtract", "closest", "coverage"):
            _fields["b"] = mo.ui.multiselect(options=_bed, label="B: BED file(s):")
        if _op == "intersect":
            _fields["mode"] = mo.ui.dropdown(options={
                "regions: the part of each A that overlaps B": "regions",
                "unique: the lines of A that touch B (-u)": "unique",
                "absent: the lines of A that do not touch B (-v)": "absent",
            }, value="regions: the part of each A that overlaps B", label="Write:")
        if _op == "subtract":
            _fields["whole"] = mo.ui.checkbox(value=False, label="Drop every A that B touches, instead of clipping it (-A)")
        if _op == "merge":
            _fields["distance"] = mo.ui.number(start=0, stop=10**9, step=1, value=0, label="Also join intervals up to this many bases apart (-d):")
        if _op in ("complement", "slop"):
            _fields["genome"] = mo.ui.dropdown(options=files_with("genome"), label="Genome file (.fai, or chromosome sizes: name and length per line):")
        if _op == "slop":
            _fields["left"] = mo.ui.number(start=0, stop=10**9, step=1, value=100, label="Bases to add before the start:")
            _fields["right"] = mo.ui.number(start=0, stop=10**9, step=1, value=100, label="Bases to add after the end:")
        _fields["output"] = mo.ui.text(placeholder=f"blank: <first A file>.{_op}.<same extension>", label="Output file:", full_width=True)
        form_sort = make_form(_fields, f"🚀 Run {_op}")
        _above = [interval_op]

    mo.vstack(_above + [form_sort])
    return (form_sort,)


@app.cell
async def sort_page_backend(asyncio, chosen, clashes, code_panel, error_panel, facts, form_sort, get_page, interval_op, mmg, mo, os, output_name, preview_panel, python_call, size_of, sort_task, split_name, warn):
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
            _panel = mo.vstack([
                mo.callout(mo.md("### ✅ Sorted\n" + facts(_rows)), kind="success"),
                preview_panel(_out),
                code_panel(f"records, runs, skipped = {python_call('sort_genomic', _paths, _out, **_kw)}"),
            ])
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
def validate_page_form(btn_refresh_validate, files_with, get_page, make_form, mo, real_files, validate_task):
    mo.stop(get_page() != "validate")
    btn_refresh_validate  # a click re-runs this cell, which lists the folder again
    _task = validate_task.value

    if _task == "Validate a file":
        form_validate = make_form({
            "file": mo.ui.dropdown(options=files_with("vcf", "bed", "gff", "sam", "fasta", "fastq"), label="File:"),
            "max_lines": mo.ui.number(start=0, stop=10**12, step=1000, value=0, label="Check only the first N lines (0 = the whole file):"),
        }, "🩺 Validate")
    elif _task.startswith("Rename"):
        _maps = {"(none)": "none"}
        for _m in real_files("mapping"):
            _maps[_m] = _m
        form_validate = make_form({
            "file": mo.ui.dropdown(options=files_with("vcf", "bed", "gff", "sam"), label="File:"),
            "style": mo.ui.dropdown(options={"chr1, chrX, chrM (UCSC)": "chr", "1, X, MT (Ensembl, NCBI)": "no_chr"},
                                    value="chr1, chrX, chrM (UCSC)", label="Rename to:"),
            "mapping": mo.ui.dropdown(options=_maps, value="(none)", label="Mapping file for other names (two columns: name, new name):"),
            "output": mo.ui.text(placeholder="blank: <file>.chr.<same extension> (or .nochr)", label="Output file:", full_width=True),
        }, "🔁 Rename")
    else:
        _refs = files_with("vcf", "sam", "gff", "genome")
        form_validate = make_form({
            "a": mo.ui.dropdown(options=_refs, label="File A:"),
            "b": mo.ui.dropdown(options=_refs, label="File B:"),
        }, "⚖️ Compare")

    form_validate
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
def vcftools_page_form(btn_refresh_vcftools, files_with, get_page, make_form, mo, real_files, vcftools_task):
    mo.stop(get_page() != "vcftools")
    btn_refresh_vcftools  # a click re-runs this cell, which lists the folder again
    _task = vcftools_task.value
    _vcfs = files_with("vcf")
    _above = []

    if _task == "Statistics":
        form_vcftools = make_form({"file": mo.ui.dropdown(options=_vcfs, label="VCF file:")}, "📊 Count")
    elif _task == "Filter":
        _above = [mo.accordion({"How to write an expression": mo.md(
            "Fields: `QUAL`, `POS`, `N_ALT` (numbers); `CHROM`, `ID`, `REF`, `ALT`, `FILTER`, `TYPE` (text); `INFO/KEY` "
            "(alone, it tests that the key is present). Operators: `==` `!=` `<` `<=` `>` `>=`, `&&` (and), `||` (or), `!` (not) "
            "and parentheses. A comparison with a missing value is false.\n\n"
            "`TYPE` is `\"snp\"`, `\"mnp\"`, `\"ins\"`, `\"del\"`, `\"indel\"`, `\"other\"` or `\"ref\"`, and works as in bcftools: "
            "`TYPE == \"snp\"` keeps a site when **every** ALT allele is a SNP, `TYPE ~ \"snp\"` when **at least one** is. "
            "`FILTER ~ \"q10\"` keeps a site whose filters include q10; `FILTER == \"q10\"` asks for exactly q10.\n\n"
            "Examples: `QUAL >= 30 && FILTER == \"PASS\"` · `TYPE == \"snp\" && INFO/DP > 10` · `!(INFO/AF < 0.01)`"
        )})]
        form_vcftools = make_form({
            "file": mo.ui.dropdown(options=_vcfs, label="VCF file:"),
            "expression": mo.ui.text(placeholder='QUAL >= 30 && FILTER == "PASS"', label="Expression (blank: no expression):", full_width=True),
            "regions": mo.ui.text(placeholder="chr1:1000000-2000000, chrX", label="Regions (blank: the whole genome):", full_width=True),
            "samples": mo.ui.text(placeholder="sample1, sample2", label="Samples to keep, in this order (blank: all of them):", full_width=True),
            "output": mo.ui.text(placeholder="blank: <file>.filtered.<same extension>", label="Output file:", full_width=True),
        }, "🔍 Filter")
    else:
        _refs = {"(none: split and trim only, no left-alignment)": "none"}
        for _f in real_files("fasta"):
            _refs[_f] = _f
        form_vcftools = make_form({
            "file": mo.ui.dropdown(options=_vcfs, label="VCF file:"),
            "split": mo.ui.checkbox(value=True, label="Split multiallelic sites into one record per allele"),
            "reference": mo.ui.dropdown(options=_refs, value="(none: split and trim only, no left-alignment)",
                                        label="Reference genome (FASTA), to check REF and left-align the indels:"),
            "output": mo.ui.text(placeholder="blank: <file>.norm.<same extension>", label="Output file:", full_width=True),
        }, "🧹 Normalize")

    mo.vstack(_above + [form_vcftools])
    return (form_vcftools,)


@app.cell
async def vcftools_page_backend(asyncio, chosen, clashes, code_panel, error_panel, facts, form_vcftools, get_page, mmg, mo, os, output_name, preview_panel, python_call, size_of, vcftools_task, warn):
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
            _per_sequence = mo.ui.table([{"sequence": _s, "records": _n} for _s, _n in _r["sequences"]],
                                        selection=None, pagination=True, page_size=10, label="Records per sequence")
            _panel = mo.vstack([
                mo.callout(mo.md(f"### 📊 {_path}\n" + facts(_rows)), kind="success"),
                _per_sequence,
                code_panel(f'print({python_call("vcf_stats", _path)}["text"])'),
            ])
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
            _panel = mo.vstack([
                mo.callout(mo.md("### ✅ Filtered\n" + facts(_rows)), kind="success"),
                preview_panel(_out),
                code_panel(python_call("vcf_filter", _path, _out, **_kw)),
            ])
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
            _items = [mo.callout(mo.md("### ✅ Normalized\n" + facts(_rows)), kind="success")]
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
def fastq_page_form(btn_refresh_fastq, fastq_task, files_with, get_page, make_form, mo, real_files):
    mo.stop(get_page() != "fastq")
    btn_refresh_fastq  # a click re-runs this cell, which lists the folder again
    _fastqs = files_with("fastq")

    if fastq_task.value == "Quality report":
        form_fastq = make_form({
            "file": mo.ui.dropdown(options=_fastqs, label="FASTQ file:"),
            "max_reads": mo.ui.number(start=0, stop=10**12, step=10000, value=0, label="Look only at the first N reads (0 = the whole file):"),
            "output": mo.ui.text(placeholder="blank: <file>.report.html", label="Report page:", full_width=True),
        }, "📈 Measure")
    else:
        _second = {"(single-end: no second file)": "none"}
        for _f in real_files("fastq"):
            _second[_f] = _f
        form_fastq = make_form({
            "file": mo.ui.dropdown(options=_fastqs, label="Reads (R1):"),
            "file2": mo.ui.dropdown(options=_second, value="(single-end: no second file)", label="Mate reads (R2), for paired-end data:"),
            "adapter": mo.ui.dropdown(options={
                "auto: the known adapter found in the first reads": "auto",
                "none: no adapter trimming": "none",
                "Illumina Universal: AGATCGGAAGAGC": "AGATCGGAAGAGC",
                "Nextera Transposase: CTGTCTCTTATACACATCT": "CTGTCTCTTATACACATCT",
                "Illumina Small RNA 3': TGGAATTCTCGGGTGCCAAGG": "TGGAATTCTCGGGTGCCAAGG",
                "BGI/MGI: AAGTCGGAGGCCAAGCGGTC": "AAGTCGGAGGCCAAGCGGTC",
                "another one: type it below": "custom",
            }, value="auto: the known adapter found in the first reads", label="3' adapter:"),
            "custom": mo.ui.text(placeholder="ACGT... (only with 'another one')", label="Adapter sequence:", full_width=True),
            "quality": mo.ui.slider(start=0, stop=40, step=1, value=20, show_value=True, label="Trim the 3' end below this Phred quality (0 = off):"),
            "min_length": mo.ui.number(start=0, stop=100000, step=1, value=20, label="Drop reads shorter than (bases):"),
            "trim_n": mo.ui.checkbox(value=False, label="Remove N bases from both ends"),
            "output": mo.ui.text(placeholder="blank: <R1 file>.trimmed.<same extension>", label="Output for R1:", full_width=True),
            "output2": mo.ui.text(placeholder="blank: <R2 file>.trimmed.<same extension>", label="Output for R2:", full_width=True),
        }, "✂️ Trim")

    form_fastq
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
def kmercount_page_ui(get_page, mo, set_page):
    mo.stop(get_page() != "kmercount")
    btn_back_kmercount = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_kmercount = mo.ui.button(label="🔄 Refresh the file list")
    _ui = mo.vstack([
        btn_back_kmercount,
        mo.md("<br>\n\n# 🔢 13. K-mer Counting"),
        mo.callout(mo.md("Count every k-mer (k from 1 to 63) of FASTA and FASTQ files of any size: what does not fit in memory is sorted in runs on disk and merged. Canonical by default, so a k-mer and its reverse complement count as one and reads from the two strands add up. The abundance spectrum shows how many distinct k-mers were seen once, twice, and so on: in sequencing reads, the k-mers seen once are mostly errors, and the peak sits at the coverage."), kind="info"),
        btn_refresh_kmercount,
    ])
    _ui
    return btn_back_kmercount, btn_refresh_kmercount


@app.cell
def kmercount_page_form(btn_refresh_kmercount, files_with, get_page, make_form, mo):
    mo.stop(get_page() != "kmercount")
    btn_refresh_kmercount  # a click re-runs this cell, which lists the folder again
    form_kmercount = make_form({
        "files": mo.ui.multiselect(options=files_with("fasta", "fastq"), label="FASTA / FASTQ files:"),
        "k": mo.ui.number(start=1, stop=63, step=1, value=31, label="k (1 to 63):"),
        "canonical": mo.ui.checkbox(value=True, label="Canonical: count a k-mer and its reverse complement as one"),
        "min_count": mo.ui.number(start=1, stop=10**9, step=1, value=1, label="Leave out of the table the k-mers seen fewer times than:"),
        "output": mo.ui.text(placeholder="blank: the summary only, no table (for example kmers.tsv.gz)", label="Table of counts:", full_width=True),
    }, "🔢 Count")
    form_kmercount
    return (form_kmercount,)


@app.cell
async def kmercount_page_backend(asyncio, chosen, clashes, code_panel, error_panel, facts, form_kmercount, get_page, mmg, mo, os, python_call, size_of, spectrum_chart, warn):
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
        _chart = spectrum_chart(_r["histogram"])
        if _chart is not None:
            _items.append(_chart)
            _items.append(mo.ui.table([{"times seen": _c, "distinct k-mers": _n} for _c, _n in _r["histogram"]],
                                      selection=None, pagination=True, page_size=10, label="The spectrum as a table"))
        _items.append(code_panel(f"result = {python_call('count_kmers', _paths, _out, **_kw)}", 'print(result["text"])'))
        _panel = mo.vstack(_items)
    except Exception as _err:
        _panel = error_panel(_err)

    _panel
    return


@app.cell
def convert_page_ui(get_page, mo, set_page):
    mo.stop(get_page() != "convert")
    btn_back_convert = mo.ui.button(label="⬅️ Back to Menu", on_change=lambda _: set_page("home"))
    btn_refresh_convert = mo.ui.button(label="🔄 Refresh the file list")
    _ui = mo.vstack([
        btn_back_convert,
        mo.md("<br>\n\n# 📦 14. Formats & Conversion"),
        mo.callout(mo.md("Change only the container of a VCF, BED, GFF/GTF or SAM file, keeping every line and its order: plain text, gzip (written as BGZF, which can also be indexed), LZ4, BAM (from and to SAM) and Parquet, a table of typed columns that pandas, polars, DuckDB or Spark open directly, with the header kept in the file, so it converts back to the same text."), kind="info"),
        btn_refresh_convert,
    ])
    _ui
    return btn_back_convert, btn_refresh_convert


@app.cell
def convert_page_form(btn_refresh_convert, files_with, get_page, make_form, mo):
    mo.stop(get_page() != "convert")
    btn_refresh_convert  # a click re-runs this cell, which lists the folder again
    form_convert = make_form({
        "file": mo.ui.dropdown(options=files_with("vcf", "bed", "gff", "sam", "parquet"), label="File:"),
        "target": mo.ui.dropdown(options={
            "plain text": "text",
            "gzip, written as BGZF (.gz)": ".gz",
            "LZ4 (.lz4)": ".lz4",
            "BAM (.bam), for SAM": ".bam",
            "Parquet (.parquet)": ".parquet",
        }, value="gzip, written as BGZF (.gz)", label="Convert to:"),
        "numbers": mo.ui.radio(options={
            "as numbers (filter and sum without conversion)": "number",
            "as text, exactly as written (a byte-for-byte round trip)": "text",
        }, value="as numbers (filter and sum without conversion)", label="QUAL (VCF) and score (GFF) in a Parquet output:"),
        "output": mo.ui.text(placeholder="blank: the same name with the new extension", label="Output file:", full_width=True),
    }, "📦 Convert")
    form_convert
    return (form_convert,)


@app.cell
async def convert_page_backend(asyncio, chosen, clashes, code_panel, error_panel, facts, form_convert, get_page, mmg, mo, os, preview_panel, python_call, size_of, split_name, warn):
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
        _out = _stem + {"text": f".{_fmt}", ".bam": ".bam"}.get(_target, f".{_fmt}{_target}")
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


if __name__ == "__main__":
    app.run()
