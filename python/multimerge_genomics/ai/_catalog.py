"""A machine-readable catalog of the MMGEN tools.

Built at run time from the engine itself (each function's signature and
docstring), plus a little curated metadata: which screen of the web interface
uses it, which area it belongs to, whether it works on human genomes, and
whether it can be called with JSON arguments (the k-mer modules pass raw bytes
between calls and cannot). An AI assistant, the agent and the MCP server all
read the tools from here, so they never drift from the engine.
"""
from __future__ import annotations

import inspect
import re

from .. import _native

# The screens of the web interface: number, page key, title, area, what it does.
SCREENS = [
    (1, "qc", "QC & Ingestion", "K-mer modules", "extract the k-mers of FASTA/FASTQ files into an .arrow file"),
    (2, "cluster", "Taxonomic Clustering", "K-mer modules", "cluster 31-mers by MinHash signatures (experimental)"),
    (3, "assembly", "Genome Assembly (De Bruijn)", "K-mer modules", "assemble contigs from k-mers"),
    (4, "variant", "Variant Calling (SNPs)", "K-mer modules", "k-mers of a patient absent from a reference, through a Bloom filter"),
    (5, "annotation", "Functional Annotation", "K-mer modules", "translate DNA into protein"),
    (6, "bloom", "Bloom Filters", "K-mer modules", "build a Bloom filter of k-mers (GPU when there is one)"),
    (7, "scaffold", "Genome Scaffolding", "K-mer modules", "join contigs by exact end overlaps; detect circular genomes"),
    (8, "vcf", "VCF Multimerge & Sorting", "File tools", "merge several VCFs into one sorted file"),
    (9, "sort", "Sort, Index & Intervals", "File tools", "sort/merge VCF, BED, GFF, SAM/BAM of any size; .tbi/.bai indexes; region queries; bedtools-like interval operations"),
    (10, "validate", "Validation & Chromosome Names", "File tools", "find what is wrong in a file; chr1 vs 1; GRCh37 vs GRCh38"),
    (11, "vcftools", "VCF Statistics, Filtering & Normalization", "File tools", "VCF summary, bcftools-like filter expressions, split and left-align variants"),
    (12, "fastq", "FASTQ Quality Report & Trimming", "File tools", "quality report (HTML) and cutadapt-like trimming"),
    (13, "kmercount", "K-mer Counting", "File tools", "count k-mers of files of any size; abundance spectrum"),
    (14, "convert", "Formats & Conversion", "File tools", "text, BGZF, LZ4, BAM, Parquet"),
    (15, "guide", "Which Test Sees What", "Human genetics", "a guide to genetic tests and their blind spots"),
    (16, "kin", "Relatedness, Trios & Paternity", "Human genetics", "KING kinship, trio consistency, paternity index"),
    (17, "anc", "Ancestry & Polygenic Scores", "Human genetics", "admixture proportions against population panels; PGS Catalog scores"),
    (18, "findings", "Disease Findings & Gene Panels", "Human genetics", "a person's variants listed in ClinVar; gene panels from a GFF"),
    (19, "pgx", "Pharmacogenomics", "Human genetics", "Star alleles, phenotypes and guideline texts from a person's VCF (and BAM, for CYP2D6 and HLA), by MMGEN's engine or PharmCAT"),
    (20, "chrom", "Sex Chromosomes & Homozygosity", "Human genetics", "genetic sex; runs of homozygosity"),
    (21, "lookup", "Variant Lookup & Allele Counts", "Human genetics", "a known variant in relatives; allele counts in a BAM"),
    (22, "rnaseq", "RNA-seq: Gene Counts & Differential Expression", "Research and microbiology", "reads per gene; DESeq2-style differential expression"),
    (23, "crispr", "CRISPR Screens (MAGeCK)", "Research and microbiology", "guide counts; gene ranking by RRA"),
    (24, "meta", "Metagenomics: Species & Abundance", "Research and microbiology", "species and abundance by sketches; read classification; host reads apart"),
    (25, "amr", "Antimicrobial Resistance Genes", "Research and microbiology", "AMRFinderPlus resistance genes in assemblies or reads"),
    (26, "outbreak", "Outbreaks: Trees, Viral Mutations & PCR", "Research and microbiology", "genome distances and trees; viral mutations named as Nextclade; in-silico PCR"),
    (27, "amplicon", "16S Amplicons: Taxonomy & Diversity", "Research and microbiology", "ASVs, RDP classifier, alpha/beta diversity, PCoA, PERMANOVA"),
    (28, "tmb", "Tumor Mutational Burden (TMB)", "Research and microbiology", "somatic variants per megabase of the sequenced region"),
    (29, "ai", "AI Assistant", "AI assistant", "describe a problem and get a path through the tools with a script; an agent that "
     "runs the tools with the person's approval; the session saved as a Python script"),
]
SCREEN_BY_PAGE = {s[1]: s for s in SCREENS}

# Which screen uses each function.
_PAGE_OF = {
    "read_fasta_batch": "qc", "read_fastq_batch": "qc", "test_giant_kmer": "qc",
    "minhash_signatures": "cluster", "sort_merge_bytes": "cluster", "read_and_merge_kmer_files": "cluster", "decode_kmers": "cluster",
    "build_debruijn_graph": "assembly", "assemble_debruijn": "assembly",
    "build_bloom_filter": "variant", "extract_mutations": "variant", "decode_kmers_128": "variant",
    "translate_proteins": "annotation", "build_scaffolds": "scaffold",
    "multimerge_vcf": "vcf",
    "sort_genomic": "sort", "index_file": "sort", "query_region": "sort", "bed_merge": "sort", "bed_intersect": "sort",
    "bed_subtract": "sort", "bed_complement": "sort", "bed_closest": "sort", "bed_coverage": "sort", "bed_slop": "sort",
    "validate_file": "validate", "normalize_chromosomes": "validate", "file_reference": "validate", "compare_references": "validate",
    "vcf_stats": "vcftools", "vcf_filter": "vcftools", "vcf_normalize": "vcftools",
    "fastq_stats": "fastq", "fastq_report": "fastq", "fastq_trim": "fastq",
    "count_kmers": "kmercount", "convert_file": "convert", "parquet_options": "convert",
    "relatedness": "kin", "trio_analysis": "kin",
    "build_panel": "anc", "panel_from_frequencies": "anc", "estimate_ancestry": "anc", "polygenic_score": "anc",
    "clinvar_findings": "findings", "gene_regions": "findings",
    "pharmcat_prepare": "pgx", "pharmcat_run": "pgx", "pgx_run": "pgx", "mmgen_pgx_run": "pgx",
    "genetic_sex": "chrom", "runs_of_homozygosity": "chrom",
    "variant_lookup": "lookup", "pileup": "lookup",
    "count_reads": "rnaseq", "differential_expression": "rnaseq",
    "crispr_count": "crispr", "crispr_test": "crispr",
    "sketch_references": "meta", "gather": "meta", "classify_reads": "meta",
    "find_resistance_genes": "amr",
    "genome_distances": "outbreak", "sketch_distances": "outbreak", "tree_from_matrix": "outbreak", "viral_mutations": "outbreak",
    "in_silico_pcr": "outbreak",
    "amplicon_asvs": "amplicon", "train_classifier": "amplicon", "classify_asvs": "amplicon", "alpha_diversity": "amplicon",
    "beta_diversity": "amplicon", "taxa_composition": "amplicon",
    "tumor_mutational_burden": "tmb",
}
# Functions that pass raw bytes between calls: a person runs them from the
# interface; an agent or an MCP client cannot call them with JSON.
_BYTES = {
    "read_fasta_batch", "read_fastq_batch", "minhash_signatures", "sort_merge_bytes", "read_and_merge_kmer_files", "decode_kmers",
    "build_debruijn_graph", "assemble_debruijn", "build_bloom_filter", "extract_mutations", "decode_kmers_128", "test_giant_kmer",
}
# Human-genome tools: their results describe a person (screens 16-21, 28).
_HUMAN_PAGES = {"kin", "anc", "findings", "pgx", "chrom", "lookup", "tmb"}
_NOT_PERSONAL = {"gene_regions"}  # gene coordinates from an annotation, no person in it

# Parameters that name a file the function writes.
OUTPUT_PARAMS = {
    "output", "output_path", "output_dir", "output_prefix", "output2", "tree_output", "matrix_output", "amplicons_output",
    "reads_output", "summary_output", "pcoa_output", "curves_output", "tracking", "asv_fasta", "feature_table", "host_free",
    "host_free2", "assignments", "details",
}
# Parameters that take one path or a list of them (or of names).
_LISTS = {
    "paths", "vcfs", "genomes", "inputs", "bams", "fastqs", "samples", "databases", "a", "b", "variants", "groups", "af_fields",
    "an_fields", "treatment", "control", "sequences", "dna_sequences", "genes", "regions", "names", "labels", "callable_beds",
    "exclude", "classes", "control_sgrnas", "control_genes", "rank_names", "host", "offsets",
}
_NUMBERS = re.compile(r"^(min_|max_|k$|memory_mb$|seed$|scaled$|threshold|rarefy|curve_|permutations$|bootstrap|block_bp$|window_|"
                      r"distance$|left$|right$|both$|padding$|overlap$|trunc_len|mismatches$|three_prime|alpha$|lfc_|quality_cutoff|"
                      r"phred_offset$|size_mb$|num_hashes$|length$|bootstraps$|curve_steps$|thin_bp$)")


def _schema_for(name: str, default) -> dict:
    if isinstance(default, bool):
        return {"type": "boolean"}
    if isinstance(default, int):
        return {"type": "integer"}
    if isinstance(default, float):
        return {"type": "number"}
    if isinstance(default, str):
        return {"type": "string"}
    if name in _LISTS:
        return {"anyOf": [{"type": "string"}, {"type": "array", "items": {"type": "string"}}]}
    if _NUMBERS.match(name):
        return {"type": "number"}
    return {"type": "string"}


def _summary(doc: str) -> str:
    first = (doc or "").strip().split("\n\n")[0]
    return " ".join(first.split())


def catalog(include_bytes: bool = True) -> list[dict]:
    """Every tool of the engine, with its screen, area, signature, documentation
    and a JSON schema of its parameters."""
    out = []
    for name in sorted(n for n in dir(_native) if not n.startswith("_") and callable(getattr(_native, n))):
        fn = getattr(_native, name)
        if not include_bytes and name in _BYTES:
            continue
        page = _PAGE_OF.get(name)
        screen = SCREEN_BY_PAGE.get(page)
        try:
            sig = inspect.signature(fn)
            params = list(sig.parameters.values())
        except (TypeError, ValueError):
            sig, params = None, []
        props, required, outputs = {}, [], []
        for p in params:
            if p.kind in (p.VAR_POSITIONAL, p.VAR_KEYWORD):
                continue
            default = None if p.default is inspect.Parameter.empty else p.default
            props[p.name] = _schema_for(p.name, default)
            if p.default is inspect.Parameter.empty:
                required.append(p.name)
            if p.name in OUTPUT_PARAMS:
                outputs.append(p.name)
        doc = inspect.getdoc(fn) or ""
        out.append({
            "name": name,
            "screen": screen[0] if screen else None,
            "page": page,
            "area": screen[3] if screen else "General",
            "human": page in _HUMAN_PAGES and name not in _NOT_PERSONAL,
            "json_callable": name not in _BYTES,
            "signature": f"{name}{sig}" if sig else f"{name}(...)",
            "summary": _summary(doc),
            "doc": doc,
            "outputs": outputs,
            "input_schema": {"type": "object", "properties": props, "required": required},
        })
    return out


def tool_definitions(human: bool = False, names=None) -> list[dict]:
    """The tools as Claude API tool definitions (name, description,
    input_schema): only those callable with JSON, and the human-genome ones
    only when `human` is true."""
    tools = []
    for t in catalog(include_bytes=False):
        if names is not None and t["name"] not in names:
            continue
        if t["human"] and not human:
            continue
        desc = t["doc"] if len(t["doc"]) <= 1800 else t["doc"][:1800] + " ..."
        if t["screen"]:
            desc = f"[screen {t['screen']}] " + desc
        tools.append({"name": t["name"], "description": desc, "input_schema": t["input_schema"]})
    return tools


def catalog_text() -> str:
    """The catalog as compact text for a model's context: the screens, then
    every function with its signature and documentation."""
    lines = ["# The screens of the MMGEN web interface", ""]
    for n, page, title, area, what in SCREENS:
        lines.append(f"{n}. {title} ({area}): {what}.")
    lines += ["", "# The functions (Python: `import multimerge_genomics as mmg`, then `mmg.<name>(...)`)", ""]
    for t in catalog():
        where = f"screen {t['screen']}" if t["screen"] else "general"
        flags = []
        if t["human"]:
            flags.append("human genome: experimental, not a diagnosis")
        if not t["json_callable"]:
            flags.append("passes raw bytes between calls")
        lines.append(f"## mmg.{t['signature']}  [{where}{'; ' + '; '.join(flags) if flags else ''}]")
        lines.append(t["doc"].strip() or t["summary"])
        lines.append("")
    return "\n".join(lines)


def page_of(function_name: str):
    return _PAGE_OF.get(function_name)
