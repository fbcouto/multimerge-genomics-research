# Contributing to MMGEN

Suggestions, bug reports and pull requests are all welcome. This page says how to send each one.

## Suggest an improvement or a new feature

[Open an issue](../../issues/new/choose) and pick **Feature request or improvement**. Start from the task: what you were trying to do, with which formats and file sizes, and what got in the way. If another tool already does what you need (bedtools, bcftools, samtools...), name the command — it makes the expected behaviour unambiguous.

## Report a bug

[Open an issue](../../issues/new/choose) and pick **Bug report**. The most useful report has:

* what you ran (the Python call, the interface screen or the notebook cell) and the full error message or validator report;
* a small input file that reproduces it — a few lines are usually enough;
* your system, Python version and the output of `python -c "import multimerge_genomics as m; print(m.engine_info())"`.

A wrong result with no error message is the most serious kind of bug here. Please report it even if you are not sure.

## Implement something yourself

1. For anything larger than a small fix, open an issue first, so the approach can be agreed before you write the code.
2. Fork the repository, create a branch, and build:

   ```bash
   python -m venv .venv
   .venv\Scripts\activate            # Windows
   # source .venv/bin/activate       # Linux / macOS
   pip install maturin
   maturin develop --release
   ```

   No CUDA toolkit is needed: without it the build is CPU-only and every module still works.
3. Run the tests that cover what you changed:

   ```bash
   cd trials
   cargo run --release                   # sorting, merging, intersection
   cargo run --release --bin validate    # validator and chromosome names
   cargo run --release --bin bloom       # CPU Bloom filter
   cargo run --release --bin formats     # BGZF, indexes, region queries, BAM
   cargo run --release --bin bed_ops     # subtract, complement, closest, coverage, slop
   cargo run --release --bin vcf_tools   # VCF statistics, filter, normalization
   cargo run --release --bin fastq_qc    # FASTQ statistics and trimming
   cargo run --release --bin kmer_count  # k-mer counting
   cd ..
   python trials/python/smoke.py                      # the Python API, no external data
   python trials/python/crosscheck.py                 # against independent tools (pyarrow, polars, biopython, bamnostic, cutadapt...)
   python trials/python/test_ui.py ui_main.py         # the interface's cells (needs marimo, polars, numpy)
   ```
4. Open a pull request against `main`, describing what changed and how you tested it.

## Ground rules

* **Everything a user reads is in English**: function and parameter names, dictionary keys, messages, docstrings, the interface and the documentation. Code comments are still mostly in Portuguese; either language is fine there.
* **Sorting stays stable**, and the file tools must keep working on inputs that do not fit in memory.
* **No silent wrong answers.** When something cannot be done, fail with a message that says what is wrong, where, and how to fix it.
* **A GPU path needs a CPU twin** that returns the same result.
* `src/multimerge_core.rs` is a vendored copy of the sorting engine. Please do not change it in a pull request: open an issue instead.
* Add a test for what you fix or add, next to the existing ones in `trials/`.

## License

By contributing you agree that your contribution is licensed under the [Apache License 2.0](LICENSE), the license of the project.
