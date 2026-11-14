# 🧬 MMGEN: MultiMerge Genomics Engine Suite

![MMGEN main menu](assets/menu.png)

MMGEN is a free, open-source genomics toolkit written in **Rust**, with optional **CUDA** acceleration and a **Python** API. It has two families of tools:

* **Everyday file tools that work at any size** — validate a file and say what is wrong with it, fix `chr1` versus `1`, tell GRCh37 from GRCh38; sort, index and query by region; merge, intersect, subtract, complement, find the closest and measure the coverage of intervals; count, filter and normalize variants; measure the quality of FASTQ reads (with an HTML report) and trim them; count k-mers. On VCF, BED, GFF/GTF, SAM/BAM, FASTA and FASTQ files that do not fit in memory, read and written as plain text, gzip (BGZF), BAM or Parquet.
* **K-mer modules** — ingestion and quality filtering, De Bruijn assembly, mutation detection with Bloom filters, translation and scaffolding, running on an NVIDIA GPU when there is one and on the CPU when there is not.

You can use it from Python, from a reactive web interface (marimo), or from a Google Colab notebook with nothing installed on your machine.

MMGEN is free (Apache 2.0). If it is useful to you, you can [support its development](#donate) — entirely optional.

---

## ✨ Highlights

* **Any size.** Sorting is bounded by a memory budget, not by the input. Two real dbSNP VCFs totalling 18.7 GB of text (75.2 million variants) are sorted and merged in 242 s within a 1 GB budget on a 4-core laptop.
* **Stable, so merges are deterministic.** Records with the same key keep their input order: file order, then line order.
* **Errors that explain themselves.** The validator gives the line, the reason and the fix. An intersection between a `chr1` file and a `1` file fails with a message that names the mismatch instead of silently returning nothing.
* **The formats you already have, in and out.** Plain text, gzip (written as BGZF, so it can be indexed), `.tbi` and `.bai` indexes for region queries, BAM, and Parquet tables that pandas, polars, DuckDB or Spark read directly. Every tool reads and writes all of them; the file extension decides.
* **Checked against the tools it replaces.** `fastq_trim` writes the same reads as cutadapt (0 differences in 128,800 reads over 84 option combinations); BAM files are read back by an independent reader, using our own `.bai`; Parquet files are read by pyarrow and polars, and theirs by us.
* **No GPU required.** Every GPU path has a CPU twin that returns the same result, bit for bit. With a GPU, the Bloom filter is 3× to 8× faster.
* **One wheel per platform.** The extension is built against Python's stable ABI and serves every Python from 3.10 on.
* **Tested.** 1,669 Rust tests, plus Python tests that drive the API with real data, compare it with independent tools, and drive the interface in a real browser.

What the validator says about one of the example files of the Colab notebook, which has defects planted on purpose:

```text
examples/with_defects.vcf: VCF with 8007 lines and 8000 records -- 4 error(s), 1 warning(s)
  ERROR   line 8: the file mixes the two chromosome naming conventions: `chr1` (line 48) and `1` (line 8)
          fix: normalize_chromosomes(...) rewrites the file in a single naming convention
  ERROR   line 18: has 7 tab-separated column(s) but 8 when split on spaces: the columns must be separated by TABS (a VCF needs at least 8: CHROM POS ID REF ALT QUAL FILTER INFO)
          fix: replace the spaces between the columns with tabs
  ERROR   line 28: POS `12mil` is not an integer
  ERROR   line 38: REF `N?` has the invalid character `?` (only A, C, G, T and N are allowed)
  WARNING line 48: sequence `chr1` is not declared in any `##contig` line of the header
  INFO    line 10: the file is not sorted: sequence `1` shows up again here, after other sequences [7014 occurrences]
          fix: sort_genomic(...) sorts the file, whatever its size
```

---

## 🚀 The Modules

Every module has a screen in the web interface, and every function is also available from Python. On screens 9–14, each result ends with the Python call that does the same thing, ready to paste into a script. The Colab notebook covers modules 9 and 10.

### 🔬 1. Data Ingestion & Quality Control
* **Batch reading, gzip included:** reads several `.fastq`, `.fasta`, `.fastq.gz` or `.fasta.gz` files in parallel and decompresses them on the fly.
* **128-bit k-mers:** k-mers of up to `k=63` packed into 16 bytes each. The output is the sorted set of distinct k-mers.
* **Quality filter:** FASTQ reads whose mean Phred quality is below the threshold are dropped before their k-mers are extracted.

### 🧬 2. Taxonomic Clustering *(experimental)*
* Computes a MinHash signature per k-mer, projects it with PCA and groups the points with HDBSCAN.
* This is the oldest module and the least validated one: it works on 64-bit k-mers (`k=31`) read from text files, and its clustering quality has not been assessed. Treat it as a demonstration.

### 🧩 3. *De Novo* Assembly (De Bruijn Graph)
* Builds the De Bruijn graph of the k-mers in memory and walks its unbranched paths to produce contigs (unitigs), written to a FASTA file.
* Collapses the simplest bubbles (two one-node branches that reconverge). Longer bubbles, tips and coverage are not handled yet.

### ⚠️ 4. Variant Detection by K-mer Difference
* Builds a Bloom filter from the reference k-mers and reports the sample k-mers that are not in it — the k-mers that carry a difference. A single SNP shows up as `k` mutated k-mers.
* Runs on the GPU or on the CPU, with the same result.

### 📘 5. Translation
* Translates DNA into protein with the standard genetic code (reading frame 1), in parallel across sequences. Stop codons are marked.

### 🔍 6. Bloom Filters
* **GPU:** a CUDA kernel sets the bits with `atomicOr`, so thousands of threads can write the same filter safely.
* **CPU:** the same filter, bit for bit, on a machine with no NVIDIA card.

### 🖇️ 7. Scaffolding & Circularization
* Joins contigs whose ends overlap exactly by a fixed number of bases, using hash maps of the ends.
* Detects sequences that close on themselves (plasmids, circular chromosomes) without looping forever.

### 🗄️ 8. VCF Multimerge & Sorting
* Select any number of `.vcf` / `.vcf.gz` files: the engine merges them into a single output sorted by position or by quality — stable, and out of core when the data does not fit in memory.

### 📏 9. Sorting, Indexing & Interval Operations at Any Size
* **`sort_genomic`** sorts VCF, BED, GFF/GTF and SAM/BAM. What fits the memory budget is sorted in memory; what does not is cut into sorted runs on disk (LZ4) and joined by a k-way merge.
* **The chromosome order you choose:** lexical (`bedtools sort`), natural (`chr2` before `chr10`), the order of the header's `##contig` / `@SQ` lines (`bcftools sort`, `samtools sort`), or the order of a `.fai` file.
* **Index and query by region:** `index_file` writes the `.tbi` of a VCF, BED or GFF and the `.bai` of a BAM — the same indexes tabix, samtools and IGV use — and `query_region` reads only the part of the file a region needs (`"chr1:1000000-2000000"`).
* **Interval operations** on inputs that need not be sorted and need not fit in memory: `bed_merge` (`bedtools merge -d`), `bed_intersect` (`bedtools intersect`, with `-u` and `-v`), `bed_subtract` (with `-A`), `bed_complement`, `bed_closest` (`-d`, every tie), `bed_coverage` and `bed_slop`.
* Malformed lines are skipped and counted, never silently sorted as zero.

### 🩺 10. File Validation & Naming Sanity
* **`validate_file`** reads a VCF, BED, GFF/GTF, SAM, FASTA or FASTQ of any size and reports each kind of problem with the **line** of its first occurrence, how many times it repeats and how to fix it.
* **`normalize_chromosomes`** rewrites a whole file, header included, in a single naming convention: UCSC (`chr1`, `chrM`) or Ensembl/NCBI (`1`, `MT`).
* **`compare_references`** compares the sequence lengths two files declare and tells GRCh37 (hg19) from GRCh38 (hg38) before the coordinates of one are used on the other.

### 🧪 11. VCF Statistics, Filtering & Normalization
* **`vcf_stats`**: records per sequence, SNPs, MNPs, insertions, deletions, multiallelic sites, transitions and transversions (Ts/Tv), PASS, missing QUAL and FILTER.
* **`vcf_filter`**: keep the records that match an expression (`QUAL >= 30 && INFO/DP > 10 && FILTER == "PASS"`), fall in some regions, or both; keep some samples, in the order you give. `TYPE` and `FILTER` follow bcftools: `TYPE == "snp"` keeps sites whose every ALT allele is a SNP, `TYPE ~ "snp"` sites with at least one; `FILTER ~ "q10"` keeps sites whose filters include q10. A typo in the expression is reported with the position of the mistake.
* **`vcf_normalize`**: split multiallelic sites into one record per allele (INFO and FORMAT fields declared `Number=A`, `R` or `G` are split with them) and trim every variant to its shortest representation. Give it the genome (`reference="GRCh38.fa"`) and it also checks REF and **left-aligns** every indel, as `bcftools norm -f` does, so files that wrote the same indel at different places of a repeat agree. The output is sorted again when positions move.

### 📈 12. FASTQ Quality Report & Trimming
* **`fastq_stats`** and **`fastq_report`**: quality along the read (median and percentiles), base composition, GC per read, N bases, duplication, read lengths and known adapters, with plain-language warnings that name the fix. The report is one HTML page that opens offline:

  ![FASTQ quality report](assets/fastq_report.png)

* **`fastq_trim`**: quality trimming of both ends, 3' adapter removal (`adapter="auto"` uses the adapter the report found), removal of N at the ends, and length, N and quality filters — single-end or paired-end, with the pairs kept in step. Same results as `cutadapt --no-indels`.

### 🔢 13. K-mer Counting at Any Size
* **`count_kmers`** counts the k-mers (k up to 63) of FASTA and FASTQ files of any size within a memory budget, canonical by default (a k-mer and its reverse complement count together), and returns the abundance spectrum (histogram) with the table of counts.

### 📦 14. Formats and Conversion
* **Every tool reads and writes every container**, chosen by the extension: plain text, `.gz` (read as gzip or BGZF, written as BGZF), `.lz4`, `.bam` (for SAM) and `.parquet` (for VCF, BED, GFF and SAM). **`convert_file`** changes only the container: `convert_file("x.sam", "x.bam")`, `convert_file("x.vcf.gz", "x.vcf.parquet")`.
* **Parquet** keeps one typed column per field (positions as integers, QUAL as a number with `.` as null, one column per sample) and the whole header in the file metadata, so it converts back to the same text. QUAL (and the GFF score) can instead be kept as the text it was, for a byte-for-byte round trip: `parquet_options(numbers="text")`. Parquet tables written by pandas, polars or pyarrow are read by their column names (`chrom`, `start`, `end`...).


---

## 💻 Quick Start

### Prerequisites
* Rust toolchain (`cargo`, `rustc`) — the stable channel is enough
* Python 3.10+
* NVIDIA CUDA Toolkit 11.8+ — optional, only for GPU acceleration
* On Windows: Visual Studio 2022 Build Tools

Pre-built wheels are not published yet; for now, build from source.

### Build and install
```bash
# from a clone of this repository
python -m venv .venv
.venv\Scripts\activate            # Windows
# source .venv/bin/activate       # Linux / macOS

pip install maturin
maturin develop --release         # builds the extension and installs it into the venv

# check the installation: any machine, no data needed, a few seconds
python trials/python/smoke.py
```

To build a wheel instead: `maturin build --release --locked -o dist`, then `pip install dist/*.whl`.

**No GPU? It still works.** Without the CUDA toolkit the build succeeds in CPU mode, and the GPU modules give the same results on the CPU. A build made with CUDA also runs on a machine with no usable NVIDIA card. `engine_info()` tells which build you have and whether a GPU is usable right now.

**Building with CUDA on Windows:** `nvcc` needs `cl.exe`. Either build from the "x64 Native Tools Command Prompt for VS 2022", or load that environment first with `call "...\VC\Auxiliary\Build\vcvars64.bat"`.

| Build variable | Effect |
| :--- | :--- |
| `MULTIMERGE_CPU_ONLY=1` | Do not try CUDA at all |
| `MULTIMERGE_CUDA_ARCH=sm_61` | Target architecture passed to `nvcc` (default: the toolkit's own default) |
| `MULTIMERGE_REQUIRE_CUDA=1` | Fail the build instead of falling back to CPU mode |

### Use it from Python
```python
import multimerge_genomics as mmg

print(mmg.engine_info())
# {'version': '0.1.0', 'gpu_compiled': True, 'gpu': True, 'gpu_memory_mb': [3371, 4095], 'threads': 8}

# What is wrong with this file?
report = mmg.validate_file("sample.vcf.gz")
print(report["text"])

# chr1 or 1? Same reference genome?
mmg.normalize_chromosomes("sample.vcf.gz", "sample.chr.vcf.gz", style="chr")
print(mmg.compare_references("sample.vcf.gz", "cohort.vcf.gz")["verdict"])

# Sort and merge files of any size, within a memory budget; index the result
records, runs, skipped = mmg.sort_genomic(["a.vcf.gz", "b.vcf.gz"], "merged.vcf.gz", memory_mb=1024, index=True)
region = mmg.query_region("merged.vcf.gz", "chr1:1000000-2000000")   # reads only that part

# Variants: summary, filter, normalize
print(mmg.vcf_stats("merged.vcf.gz")["text"])
mmg.vcf_filter("merged.vcf.gz", "pass.vcf.gz", expression='FILTER == "PASS" && QUAL >= 30')
mmg.vcf_normalize("pass.vcf.gz", "pass.norm.vcf.gz", reference="genome.fa")   # split, trim, left-align

# Intervals (inputs are lists of files)
n_read, n_written = mmg.bed_merge(["peaks.bed"], "peaks.merged.bed", distance=100)
n_lines = mmg.bed_intersect(["peaks.bed"], ["genes.bed"], "peaks_in_genes.bed", mode="regions")
mmg.bed_closest(["peaks.bed"], ["genes.bed"], "nearest_gene.bed")

# Reads: quality report, then trimming of a pair of files
mmg.fastq_report("sample_R1.fastq.gz", "sample_R1.report.html")
mmg.fastq_trim("sample_R1.fastq.gz", "trimmed_R1.fastq.gz", input2="sample_R2.fastq.gz", output2="trimmed_R2.fastq.gz",
               adapter="auto", quality_cutoff=20, min_length=30)
print(mmg.count_kmers(["trimmed_R1.fastq.gz", "trimmed_R2.fastq.gz"], "kmers.tsv.gz", k=31)["text"])

# Change the container: BAM, Parquet
mmg.convert_file("aligned.sam", "aligned.bam")
mmg.convert_file("pass.norm.vcf.gz", "pass.norm.vcf.parquet")   # pandas.read_parquet() opens it
```

Every function has a docstring: `help(mmg.sort_genomic)`.

### Launch the web interface
```bash
pip install marimo polars numpy scikit-learn pyarrow
marimo run ui_main.py
```

The interface lists the files of the folder it is started from, and writes its results there. The menu has two groups: the k-mer modules (1–7) and the file tools (8–14). On screens 9–14 you pick a task, fill in a form and get the result with a look at the first records of the output; leave an output name blank and one is made from the input name. **Refresh the file list** picks up files written since the screen was opened.

![Variant detection screen](assets/variant_calling.png)

### Google Colab
[`colab/MultiMerge_Genomics_Colab.ipynb`](colab/MultiMerge_Genomics_Colab.ipynb) runs the validation, naming, sorting and interval tools in the browser: one form per task, with the code one double-click away. It creates its own example files, some with defects planted on purpose. Until a wheel is published, its first cell asks for a Linux wheel (`.whl`) of the package, which you can build with `maturin build` on any Linux machine.

#### The web interface on Colab

The marimo interface runs on Colab too, displayed inside the notebook. With the package installed and `ui_main.py` uploaded to the session (put your data files in the same folder), run this in a cell:

```python
!pip install -q marimo

import subprocess, time, urllib.request
subprocess.Popen(["python", "-m", "marimo", "run", "ui_main.py", "--headless", "--port", "2718", "--no-token"], cwd="/content")
for _ in range(60):  # wait until the server answers
    try:
        urllib.request.urlopen("http://127.0.0.1:2718/", timeout=2)
        break
    except OSError:
        time.sleep(1)

from google.colab import output
output.serve_kernel_port_as_iframe(2718, height=900)
```

> **Blank output? Allow third-party cookies for Colab.** Colab shows the interface in an iframe that depends on them. When the browser blocks them, the cell runs without any error and the output simply stays empty. In Chrome, click the crossed-out eye icon at the right end of the address bar and allow third-party cookies for the site; the page reloads, the session stays connected, and running the cell again shows the menu.

`--no-token` is fine here: the port is reachable only through your own Colab session.

---

## 📚 Python API

| Function | What it does |
| :--- | :--- |
| `engine_info()` | Version, whether CUDA was compiled in, whether a GPU is usable now, GPU memory, threads |
| `validate_file(path, format=None, max_lines=None)` | Validates VCF, BED, GFF/GTF, SAM, FASTA or FASTQ. Returns a dictionary: `valid`, `errors`, `warnings`, `lines`, `records`, `complete`, `sorted`, `naming`, `problems`, `text` |
| `normalize_chromosomes(input, output, style="chr", format=None, mapping=None)` | Rewrites a file with `chr1` (`"chr"`, `"ucsc"`) or `1` (`"no_chr"`, `"ensembl"`, `"ncbi"`) names |
| `file_reference(path)` | What a file declares about its reference: `assembly`, `declared`, `sequences` |
| `compare_references(a, b)` | Same reference genome? `compatible` (True / False / None) and a one-sentence `verdict` |
| `sort_genomic(paths, output, format=None, by="position", chrom_order=None, memory_mb=None, temp_dir=None, index=False)` | Sorts and merges VCF, BED, GFF/GTF or SAM/BAM of any size; `index=True` also writes the `.tbi` / `.bai`. Returns `(records, runs_on_disk, skipped_lines)` |
| `index_file(path, format=None)` | Writes the `.tbi` of a BGZF VCF/BED/GFF or the `.bai` of a BAM. Returns `index`, `kind`, `sequences`, `records`, `skipped` |
| `query_region(path, regions, output=None, header=True)` | The records of an indexed file that overlap `"chr1:1000-2000"` (or a list of regions): a list of lines, or written to `output` |
| `convert_file(input, output, format=None)` | Same lines, another container: text, `.gz` (BGZF), `.lz4`, `.bam`, `.parquet`. Returns `lines`, `records` |
| `parquet_options(numbers=None)` | How every `.parquet` is written from then on: QUAL and score as numbers (`"number"`, the default) or exactly as written (`"text"`). Returns the settings |
| `bed_merge(paths, output, distance=0, ...)` | Merges overlapping or nearby intervals. Returns `(read, written)` |
| `bed_intersect(a, b, output, mode="regions", ...)` | Intersects A with B; `mode` is `"regions"`, `"unique"` (`-u`) or `"absent"` (`-v`). Returns the number of lines |
| `bed_subtract(a, b, output, remove_whole=False, ...)` | What is left of A once B is removed (`bedtools subtract`, `-A`). Returns the number of lines |
| `bed_complement(paths, genome, output, ...)` | The stretches of the genome no interval covers; `genome` is a `.fai` or chrom.sizes file |
| `bed_closest(a, b, output, ...)` | The nearest B of each A and the distance (`bedtools closest -d`), every tie |
| `bed_coverage(a, b, output, ...)` | For each A: how many B overlap it, bases covered, length, covered fraction (`bedtools coverage`) |
| `bed_slop(paths, genome, output, both=None, left=0, right=0)` | Grows every interval, within the sequence ends (`bedtools slop`) |
| `vcf_stats(path)` | Records, SNPs, MNPs, indels, multiallelic sites, Ts/Tv, PASS... and `text`, the report |
| `vcf_filter(input, output, expression=None, regions=None, samples=None)` | Keeps the records that match an expression and/or fall in regions; keeps some samples |
| `vcf_normalize(input, output, split_multiallelic=True, reference=None, ...)` | Splits multiallelic sites and trims alleles to their shortest form; with `reference` (a FASTA), checks REF and left-aligns indels; sorts again when needed |
| `multimerge_vcf(paths, sort_field, output_path, memory_mb=None, temp_dir=None)` | Merges VCFs sorted by `"POS"` or `"QUAL"` |
| `fastq_stats(path, max_reads=None)` / `fastq_report(path, output, max_reads=None)` | Read quality, composition, GC, N, duplication, lengths, adapters, warnings; `fastq_report` also writes the HTML page |
| `fastq_trim(input, output, input2=None, output2=None, adapter=None, ...)` | Quality, adapter and N trimming and the length/N/quality filters, single- or paired-end; `adapter="auto"` detects it |
| `count_kmers(paths, output=None, k=31, canonical=True, min_count=1, memory_mb=None, temp_dir=None)` | K-mer counts of FASTA/FASTQ of any size and the abundance histogram |
| `read_fasta_batch(paths, k)` / `read_fastq_batch(paths, min_quality, k)` | Distinct k-mers of the files, sorted, as bytes (16 per k-mer) |
| `build_debruijn_graph(kmer_bytes, k)` / `assemble_debruijn(kmer_bytes, k)` | Graph size `(nodes, edges)` / the contigs, longest first |
| `build_bloom_filter(kmer_bytes, size_mb, num_hashes, device="auto")` | Bloom filter of the k-mers, as bytes |
| `extract_mutations(sample_bytes, bloom_bytes, num_hashes, device="auto")` | The sample k-mers that are not in the filter |
| `translate_proteins(dna_sequences)` | One protein per DNA sequence |
| `build_scaffolds(sequences, overlap)` | `(scaffolds, circular)` |
| `decode_kmers_128(kmer_bytes, k)` | K-mers back to DNA strings |

* `by`: `"position"` (every format), `"quality"` (VCF, descending) or `"name"` (SAM, QNAME in natural order).
* `chrom_order`: `"lexical"`, `"natural"`, `"header"`, or the path of a `.fai` file. Default: the header for VCF and SAM, lexical for BED and GFF.
* `device`: `"auto"` (the fastest path this machine has), `"gpu"` (GPU only; raises an error if it cannot run there) or `"cpu"`.
* Inputs and outputs: plain text, `.gz` (gzip or BGZF; written as BGZF), `.lz4`, `.bam` and `.parquet`, by extension. A file name such as `x.vcf.parquet` or `x.bed.gz` also tells the format.
* `regions`: `"chr1"`, `"chr1:1000"` or `"chr1:1000-2000"`, 1-based with the end included, as in samtools and tabix.
* The GIL is released during file operations, so a user interface stays responsive.
* Four more functions serve the legacy 64-bit path of module 2: `sort_merge_bytes`, `read_and_merge_kmer_files`, `minhash_signatures` and `decode_kmers`.

---

## ⚡ GPU: what it speeds up, and what it does not

Bloom filter, 128-bit k-mers, 64 MB filter, 5 hashes; best of 3 runs.

| Hardware | k-mers | Build: CPU → GPU | Query: CPU → GPU |
| :--- | ---: | :--- | :--- |
| A100 vs 12 virtual cores | 10 M | 485 → 128 ms (**3.8×**) | 260 → 56 ms (**4.7×**) |
| A100 vs 12 virtual cores | 50 M | 2115 → 272 ms (**7.8×**) | 1276 → 211 ms (**6.0×**) |
| GTX 1050 Ti vs i7-7700HQ | 10 M | 767 → 248 ms (**3.1×**) | 583 → 125 ms (**4.7×**) |

* **Bloom filters gain from the GPU**, and more as the data grows.
* **Merging two sorted blocks does not**: moving the data to the card and back costs more than the merge itself, on both machines. `"auto"` merges on the CPU.
* **Sorting text files (VCF, BED, GFF, SAM) is not a GPU problem**: the measured bottleneck is reading and parsing lines.

## 📊 Measured on real data

4-core laptop (i7-7700HQ), about 3.5 GB of free RAM:

| Input | Size | Result |
| :--- | :--- | :--- |
| dbSNP b151 common, GRCh38 | 9.42 GB, 37.3 M variants | sorted in 164 s with a 1 GB budget (13 runs) |
| The two dbSNP b151 common files | 18.7 GB, 75.2 M variants | sorted and merged in 242 s with a 1 GB budget (25 runs) |
| RefSeq GRCh38.p14 GFF3 | 1.60 GB, 4.93 M records | sorted in 12.8 s with a 256 MB budget (8 runs) |

---

## 🧪 Testing

The tests live in `trials/`.

```bash
# Rust, standalone: needs neither CUDA nor Python
cd trials
cargo run --release                   # 102 tests: sorting, merging, intersection
cargo run --release --bin validate    # 30 tests: validator and chromosome names
cargo run --release --bin bloom       # 100 tests: CPU Bloom filter against a naive reference
cargo run --release --bin formats     # 37 tests: BGZF, tabix and BAI indexes, region queries, BAM
cargo run --release --bin bed_ops     # 308 tests: subtract, complement, closest, coverage, slop
cargo run --release --bin vcf_tools   # 188 tests: VCF statistics, filter expressions, normalization, left-alignment
cargo run --release --bin fastq_qc    # 592 tests: FASTQ statistics and trimming
cargo run --release --bin kmer_count  # 312 tests: k-mer counting, in memory and on disk
# with CUTADAPT=<path of cutadapt> the FASTQ tests also compare every read with cutadapt

# Python, with the extension installed
python trials/python/smoke.py                      # smoke test, no external data
python trials/python/crosscheck.py                 # against pyarrow, polars, biopython, bamnostic, cutadapt
                                                   #   (and samtools, bcftools, bedtools, tabix, jellyfish when present)
python trials/python/test_ui.py ui_main.py         # the interface's cells, without a browser
python trials/python/test_browser.py ui_main.py    # the interface in a real browser (needs playwright)
python trials/python/gpu_test.py 1 10              # GPU against CPU against numpy; sizes in millions
```

Every Rust test compares the complete output with an independent reference written in the test itself, on random inputs from a fixed seed, and runs again with a tiny memory budget so that the on-disk path is exercised too.

## 🚧 Status and known limitations

MMGEN is alpha software (version 0.1.0). The file tools of modules 8–14 are the most mature part; the k-mer modules 1–7 are closer to a teaching and research prototype.

* **The k-mer modules 1–7 work in memory, and their k-mers are strand-specific.** Assembly and variant detection suit same-strand input, such as an assembled reference against a mutated copy of it. `count_kmers` (module 13) is the exception: any size, canonical.
* **Variant detection reports k-mers**, not VCF records.
* **CRAM is not read.** BAM is, and written. **BCF** is not.
* **Indexes:** `.tbi` and `.bai`, which address positions up to 512 Mbp per sequence (enough for human; not for some plant chromosomes). The CSI index is not written.
* **Parquet**, by default, keeps QUAL (VCF) and score (GFF) as numbers, so they come back in their shortest form: `12.50` as `12.5`, `30.0` as `30`. `parquet_options(numbers="text")` keeps them as written.
* **`vcf_normalize` left-aligns only when given the genome** (`reference=`); without it, it splits and trims. Records whose REF does not match the genome are trimmed, not aligned, and counted in `ref_mismatch`. A multiallelic site is aligned only when split. In a split record, a genotype allele that belongs to another record becomes missing (`1/2` gives `1/.` and `./1`).
* **`vcf_filter`** understands a subset of the bcftools expression language: no FORMAT fields, no functions, and `~` only for TYPE and FILTER (no regular expressions). Variant types are decided by allele length, without trimming first: `AC>AT` is an MNP here and a SNP to bcftools.
* **`fastq_trim`** removes 3' adapters without insertions or deletions (`cutadapt --no-indels`); 5' and linked adapters are not supported.
* `bed_merge` writes three columns (chromosome, start, end), like `bedtools merge` without `-c`.
* Reading and parsing lines uses a single core; that is where sorting time goes. Compressing and decompressing BGZF uses every core.
* Module 2 (clustering) is experimental — see above.
* The public API, the messages, the interface and this documentation are in English. Code comments and internal names are still in Portuguese.

## 💬 Feedback and contributions

Suggestions for improvements and requests for new features are welcome: [open an issue](../../issues/new/choose) and pick *Feature request or improvement*, or *Bug report* if something went wrong. To implement something yourself, see [CONTRIBUTING.md](CONTRIBUTING.md); pull requests are welcome.

## 📄 License

[Apache License 2.0](LICENSE).

---

<a id="donate"></a>

## 💛 Support the project (optional)

MMGEN is free and open source, and nothing in it is locked behind a donation. If it saved you time and you would like to help keep it maintained, a contribution is welcome at either address below.

| Currency | Network | Address |
| :--- | :--- | :--- |
| Bitcoin (BTC) | Bitcoin | `12Pea1apDcdMfeC1LgYHxW14QAPeyrL66h` |
| Ether (ETH) | Ethereum (ERC-20) | `0x7e4fc1627ee28bc6622c3e8209c7fc33157f3723` |

Please send only over the network shown for each address: a transfer made over a different network may not arrive, and it cannot be reversed. Copy the address rather than typing it, and compare its first and last characters before confirming.

Other ways to help cost nothing: star the repository, [report a bug](../../issues/new/choose) with a file that reproduces it, or [say which format or tool you miss](../../issues/new/choose).
