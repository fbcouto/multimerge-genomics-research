# 🧬 MMGEN: MultiMerge Genomics Engine Suite

![MMGEN Main Menu](assets/menu.png)

A high-performance metagenomic binning, assembly, and variant calling engine designed to handle millions of DNA sequences (K-mers) with extreme efficiency. Built on a cutting-edge hybrid architecture using **Rust**, **CUDA (C++)**, and **Python**, MMGEN eliminates traditional computational bottlenecks through Zero-Copy memory transfers and hardware-accelerated algorithms.

Available as a native Rust library and a fully reactive Web UI.

---

## 🚀 The 8 Core Modules

The MMGEN Suite is divided into 8 specialized modules, seamlessly integrated via a reactive UI:

### 🔬 1. Fast Data Ingestion & Quality Control (QC)
* **Batch GZIP On-the-Fly:** Utilizes `flate2` and `Rayon` to read and decompress multiple `.fastq.gz` or `.fasta.gz` files simultaneously directly into RAM, bypassing disk bottlenecks.
* **128-bit Dynamic K-mers:** Employs a custom `Kmer128` memory foundation, allowing for ultra-long K-mers (up to `k=63`) mapped perfectly into 16-byte boundaries.
* **Phred Quality Trimming:** Evaluates base quality scores on the fly, eliminating sequencing noise before converting valid DNA sequences.

### 🧬 2. Spatial Taxonomy & Metagenomic Clustering
* **MinHash Avalanche Hashing:** Computes multi-dimensional spatial signatures (MinHash) natively in Rust.
* **HDBSCAN Binning:** Utilizes CPU power to calculate spatial distances, grouping K-mers into distinct biological clusters (Species/Strains).

### 🧩 3. *De Novo* Genome Assembly (De Bruijn Graphs)
* **RAM-Optimized Graph Construction:** Slices purified K-mers into Prefixes and Suffixes, mapping edges and nodes entirely in Rust's memory.
* **Intelligent Bubble Popping (Tourbus Algorithm):** Native graph intelligence that detects diamond-shaped paths caused by diploid SNPs or reading errors, collapsing the weakest branch to prevent Contig fragmentation.
* **Contig Extraction:** Safely traverses the De Bruijn graph to reconstruct massive linear DNA fragments without infinite loops.

### ⚠️ 4. Hardware-Accelerated Variant Calling (SNPs)
* **Zero-Copy Mutation Detection:** Compares a patient's sample against a healthy reference genome entirely within the GPU's VRAM.
* **Bloom Filter Querying:** Employs CUDA to isolate mutated sequences, SNPs, and biological errors in milliseconds.

### 📘 5. Functional Annotation & Translation
* **Multi-Core Protein Factory:** Uses the `Rayon` work-stealing scheduler to translate DNA sequences into Amino Acid chains (Proteins) utilizing 100% of the CPU cores.
* **Universal Codon Mapping:** Highly optimized genetic dictionary mapping triplets to proteins (including Stop Codons).

### 🔍 6. Probabilistic Bloom Filters (cuSBF)
* **Massive GPU Hashing:** Deploys thousands of CUDA threads to build probabilistic data structures directly in the VRAM.
* **Atomic Memory Operations:** Uses `atomicOr` instructions to prevent thread collisions, enabling instant sequence-membership queries.

### 🖇️ 7. Genome Scaffolding & Circularization
* **Greedy O(N) Overlap-Layout-Consensus:** Scans millions of Contigs and perfectly aligns their suffix-prefix overlaps using memory-efficient Rust HashMaps.
* **Circular DNA Detection:** Mathematically identifies and safely closes plasmids and circular chromosomes without triggering infinite processing loops, outputting an export-ready FASTA file.

### 🗄️ 8. VCF Multimerge & Sorting (Legacy Engine Integration)
* **Massive Variant Sorting:** Native integration of a highly optimized, adaptive Multimerge Rust CLI engine to process millions of genetic variants.
* **Polars Parquet Bypass:** Automatically detects `.parquet` formats and routes them through a native Polars Zero-Copy engine for lightning-fast schema inference and sorting.
* **Smart Allele Collapsing:** Identifies and collapses local diploid redundancies directly on the VCF fields.
---

## ⚙️ Zero-Copy Memory Architecture
MMGEN eliminates the *Out-of-Memory (OOM)* crashes common in Python-based bioinformatics tools. Data processed in the Rust/CUDA backend is passed to the Python frontend (Polars) via **Apache Arrow** memory pointers—meaning zero serialization overhead, zero RAM duplication, and instant DataFrame rendering.

---

## 💻 Installation & Usage

### Prerequisites
* Rust Toolchain (`cargo`, `rustc`)
* Python 3.10+
* NVIDIA CUDA Toolkit 11.8+ (Required for GPU acceleration features)
* Visual Studio 2022 Build Tools (For Windows users)

### Build for Development
```bash
# Create and activate a virtual environment
python -m venv .venv
.venv\Scripts\activate

# Install dependencies
pip install maturin marimo polars scikit-learn numpy pyarrow

# Build the Rust/CUDA extensions natively
cargo build --release --bin multimerge_cli
maturin develop --release

```

### Launch the Reactive UI

```bash
marimo run ui_main.py

```

```
