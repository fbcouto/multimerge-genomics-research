## What this changes

<!-- What the change does and why. Link the issue it answers, if there is one: "Closes #123". -->

## How it was tested

- [ ] `cd trials && cargo run --release`, plus the `--bin` of each module you touched (`validate`, `bloom`, `formats`, `bed_ops`, `vcf_tools`, `fastq_qc`, `kmer_count`)
- [ ] `python trials/python/smoke.py` (and `python trials/python/crosscheck.py` if you touched a file format or a tool that has a reference equivalent)
- [ ] For interface changes: `python trials/python/test_ui.py ui_main.py`
- [ ] Anything a user reads (messages, API names, docstrings, interface) is in English
