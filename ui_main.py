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

    # What this installation really has: the GPU modules fall back to the CPU.
    _engine = mmg.engine_info()
    _where = "NVIDIA GPU (CUDA)" if _engine["gpu"] else "CPU"

    menu_screen = mo.md(f"""
    <div style="text-align: center; margin-top: 40px; font-family: sans-serif;">
        <h1 style="font-size: 3em; margin-bottom: 0px; text-align: center;">🧬 MMGEN</h1>
        <h3 style="color: gray; margin-top: 5px;">High-Performance Genomics Framework</h3>
        <p style="margin-bottom: 40px;">Rust engine with optional NVIDIA CUDA acceleration. Pick a module to start:</p>

        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 15px; max-width: 850px; margin: 0 auto;">
            {btn_menu_qc}
            {btn_menu_cluster}
            {btn_menu_assembly}
            {btn_menu_variant}
            {btn_menu_annotation}
            {btn_menu_bloom}
            {btn_menu_scaffold}
            {btn_menu_vcf}
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
        btn_menu_qc,
        btn_menu_scaffold,
        btn_menu_variant,
        btn_menu_vcf,
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


if __name__ == "__main__":
    app.run()
