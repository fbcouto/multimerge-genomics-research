# MMGEN notebooks for researchers

A starting point, not a catalog: short notebooks for common problems, each one runnable as it is on public demo
data and written to be adapted to your own files. Organize them in folders and add your own (see
[Add your own problem](#add-your-own-problem)).

| Area | Problem | Notebook | Colab | Data used |
| :--- | :--- | :--- | :--- | :--- |
| Files | Two VCFs from different pipelines (chr1 vs 1, build, representation): make them comparable and find the shared variants | [harmonize_two_vcfs.ipynb](files/harmonize_two_vcfs.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fbcouto/multimerge-genomics-research/blob/main/notebooks/files/harmonize_two_vcfs.ipynb) | GIAB HG002 calls and 1000 Genomes frequencies (human pack) |
| Files | A gene list to its regions in an annotation, and only those regions of a large VCF | [extract_genes_from_vcf.ipynb](files/extract_genes_from_vcf.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fbcouto/multimerge-genomics-research/blob/main/notebooks/files/extract_genes_from_vcf.ipynb) | GIAB HG002 calls and RefSeq GFF of panel genes (human pack) |
| Reads | Paired-end FASTQ: quality report, adapter and quality trimming, before and after | [qc_and_trimming.ipynb](reads/qc_and_trimming.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fbcouto/multimerge-genomics-research/blob/main/notebooks/reads/qc_and_trimming.ipynb) | 30,000 read pairs of the ZymoBIOMICS mock community (research pack) |
| Variants | What a VCF holds (types, Ts/Tv, PASS), filter by an expression, split multiallelic records | [summarize_filter_normalize.ipynb](variants/summarize_filter_normalize.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fbcouto/multimerge-genomics-research/blob/main/notebooks/variants/summarize_filter_normalize.ipynb) | GIAB HG002 calls (human pack) |
| RNA-seq | Counts and a sample table to DESeq2-style differential expression, volcano and MA plots, top genes | [differential_expression.ipynb](rna_seq/differential_expression.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fbcouto/multimerge-genomics-research/blob/main/notebooks/rna_seq/differential_expression.ipynb) | airway, GEO GSE52778 (research pack) |
| Outbreaks | Viral genomes against the reference: mutations named as Nextclade, distances, tree | [sars_cov_2_mutations_and_tree.ipynb](outbreaks/sars_cov_2_mutations_and_tree.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fbcouto/multimerge-genomics-research/blob/main/notebooks/outbreaks/sars_cov_2_mutations_and_tree.ipynb) | 20 SARS-CoV-2 genomes and NC_045512.2 (research pack) |
| Microbiology | Bacterial isolates against the AMRFinderPlus catalog: resistance genes and drug classes | [resistance_genes.ipynb](microbiology/resistance_genes.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fbcouto/multimerge-genomics-research/blob/main/notebooks/microbiology/resistance_genes.ipynb) | 4 K. pneumoniae and S. aureus assemblies, AMRFinderPlus catalog (research pack) |
| Human genetics | Sample QC of a cohort: relatedness (KING), genetic sex, trio check | [cohort_sample_qc.ipynb](human/cohort_sample_qc.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/fbcouto/multimerge-genomics-research/blob/main/notebooks/human/cohort_sample_qc.ipynb) | 7 Genome in a Bottle genomes, 1000 Genomes frequencies (human pack) |

Every notebook has the same sections: the problem, what you need and what you get; setup; data (your files or
the demo); numbered steps, each explained (the method and what can go wrong); checks of the result; a methods
paragraph with the version and the parameters used; how to adapt it; references.

## How to run them

**On Google Colab** (nothing to install): click the badge, then *Runtime > Run all*. The first cell installs
MMGEN from the [v2.1.0 release](https://github.com/fbcouto/multimerge-genomics-research/releases/tag/v2.1.0) (about a minute); the data cell downloads the demo pack it needs
once (about 45 MB each).

**On your computer**: install the MMGEN wheel for your system from the [release](https://github.com/fbcouto/multimerge-genomics-research/releases/tag/v2.1.0)
(`python -m pip install <the .whl file>`), plus `python -m pip install jupyter matplotlib pandas`, and open the
notebook with Jupyter in a working folder. The notebooks write the demo files to `data/` and their results to
`results/<notebook name>/`, both in the folder where they run.

**Without Jupyter**: the code cells are plain Python (no `!pip`, no `%` magics), so they also run in order with
`python trials/python/run_notebook.py notebooks/<folder>/<name>.ipynb`. That is how this collection is tested.

## Use your own data

In each notebook, the data cell starts with the variables to set: `MY_VCF = "my_file.vcf.gz"` and so on.
Left as `None`, the demo data is used. With your files, the checks still run, and the demo's expected numbers are
not asserted.

## Add your own problem

```
python notebooks/new_problem.py my_problems/coverage_of_my_panel "Coverage of my gene panel"
```

This creates `notebooks/my_problems/coverage_of_my_panel.ipynb` from the template (`_template/new_problem.ipynb`),
with the title filled in, and the folder if it does not exist. Any folder works (`notebooks/proteomics/...`);
an existing notebook is never overwritten. The Colab badge in the new notebook points to this repository: change
the user and repository in its link if the notebook lives in your fork.

Conventions, so that every notebook reads the same way:

* one problem per notebook, about 15 cells, with the sections of the template;
* inputs at the top of the data cell (`MY_... = None` uses the demo data), results in `results/<name>/`;
* plain Python in the code cells: no shell or IPython magics, plots and tables optional (`plt`, `pd` may be `None`);
* checks print `OK` or `CHECK`; the demo data's expected numbers are asserted only when the demo data is used;
* the methods paragraph prints the MMGEN version, the functions and the parameters actually used;
* a notebook on human genomes carries the notice below, at the top and in its methods text;
* save notebooks without outputs (*Edit > Clear all outputs*, or `jupyter nbconvert --clear-output --inplace`).

To share a notebook with everyone, open a pull request (see `CONTRIBUTING.md`) that adds it and its row in the
table above.

## The demo data

Two packs, attached to the [release](https://github.com/fbcouto/multimerge-genomics-research/releases/tag/v2.1.0); each has a `README.md` with the sources and licenses:

* `mmgen_demo_human_v3.zip`: the Genome in a Bottle (NIST) benchmark calls of seven public genomes (HG001-HG007, volunteers
  who consented to open sharing) cut to chromosome 22 and the genes of some clinical panels, with their callable
  regions; ClinVar records of the same regions; 1000 Genomes allele frequencies; the RefSeq annotation of the
  panel genes; PGS Catalog scores.
* `mmgen_demo_research_v4.zip`: airway RNA-seq counts; MAGeCK CRISPR demos; 20 SARS-CoV-2 genomes with the reference and its
  annotation; the AMRFinderPlus catalog and four bacterial assemblies; ZymoBIOMICS mock-community reads and
  sketches; a PCR primer panel; SEQC2 tumor calls.

> **Notebooks on human genomes.** Experimental software, for teaching and research only: not a medical diagnosis, not a paternity test with legal value and not advice on any treatment or medicine. The authors take no responsibility for decisions made on the basis of these results; a finding must be confirmed by an accredited laboratory and interpreted by a qualified professional. Analyze only genomes of people who agreed to it.

## Cite

Couto F. MMGEN: MultiMerge Genomics (version *the one you used*, printed by `mmg.engine_info()`). Zenodo. https://doi.org/10.5281/zenodo.20417492

Cite also the methods and data each notebook lists in its References.
