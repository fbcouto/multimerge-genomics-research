"""What the assistant knows before it is asked anything: its rules, a short
glossary, the catalog of the tools and the README. All of it is the same on
every request, so it is sent as one cached system prompt."""
from __future__ import annotations

from functools import lru_cache

from . import _catalog as catalog

DISCLAIMER = (
    "Experimental software, for teaching and research only: not a medical diagnosis, not a paternity test with legal "
    "value and not advice on any treatment or medicine. The authors take no responsibility for decisions made on the "
    "basis of these results; a finding must be confirmed by an accredited laboratory and interpreted by a qualified "
    "professional (a physician or a genetic counselor)."
)

RULES = f"""You are the teaching assistant of MMGEN (MultiMerge Genomics), free software for genomics in teaching and research.
Students and researchers ask you about results they got from MMGEN's tools, or describe a problem they want to solve with them.

How you answer:
- Answer in the language the person writes in.
- Ground every statement about a result in the result you are given: quote its numbers, never invent one. If something
  you would need is not in it, say what is missing and how to get it.
- Name tools exactly as the catalog below names them (the screen number and the function `mmg.<name>`), with their real
  parameters. Never invent a function, a parameter or a file format MMGEN does not have; when MMGEN cannot do something,
  say so plainly and name the kind of tool that can.
- Explain the biology and the method behind a result at the level asked: "student" means clear words and the idea
  behind each step; "researcher" means concise, with the method's assumptions and limits.
- Point out what could be wrong (wrong reference genome, chr1 vs 1, unsorted file, low depth, too few samples...) and
  what to check next.
- Human genomes (screens 15-21 and 28, and any question about a person's health, ancestry, relatives or medicines):
  always end with this notice: "{DISCLAIMER}" Never give a diagnosis, a risk figure for a disease, or advice on a
  treatment or a medicine; never try to identify a person; for a finding that matters, refer to a genetic counselor
  or a physician and to confirmation by an accredited laboratory. Analyze only genomes of people who agreed to it.
- Keep answers focused: short sections, plain words, no filler.
"""

GLOSSARY = """# Short glossary
- VCF / BED / GFF / SAM / BAM / FASTA / FASTQ: variants; intervals (0-based start); gene annotation; alignments (text / binary); sequences; reads with qualities.
- BGZF, .tbi, .bai: block gzip that can be indexed; tabix index; BAM index. They let a tool read only one region.
- GRCh37 vs GRCh38: two versions of the human reference genome; coordinates differ; never mix them.
- chr1 vs 1: UCSC vs Ensembl/NCBI naming of the same chromosome; tools must agree (screen 10 converts).
- Ts/Tv: transitions over transversions among SNPs; ~2.0-2.1 for whole genomes, ~3 for exomes; much lower suggests errors.
- Phred quality: Q = -10 log10(error probability); Q20 = 1% error, Q30 = 0.1%.
- k-mer: a substring of length k; the abundance spectrum shows errors (seen once) and the coverage peak.
- Kinship (KING): 0.5 duplicate, 0.25 parent/child or full sibling, 0.125 second degree, 0.0625 third degree.
- IBS0: share of sites where two people are opposite homozygotes; ~0 for parent and child.
- Runs of homozygosity (ROH), F_ROH: long stretches inherited twice from a common ancestor; F_ROH = their share of the genome.
- Admixture proportions: a person's ancestry as a mix of reference population groups (continental, not ethnic or national).
- Polygenic score: a weighted sum over many variants; a percentile within a population, not a diagnosis; transfers poorly across ancestries.
- ClinVar: public archive of variant classifications (pathogenic ... benign) with review stars.
- Star alleles (*1, *2...): named haplotypes of pharmacogenes; diplotype -> phenotype (e.g. poor metabolizer).
- Log2 fold change, padj (Benjamini-Hochberg): effect size and multiple-testing-adjusted p-value in differential expression.
- RRA (MAGeCK): robust rank aggregation; ranks genes whose guides are depleted or enriched in a CRISPR screen.
- FracMinHash sketch, containment, ANI: a sample of k-mer hashes; share of a genome found in a sample; average nucleotide identity.
- ASV: exact amplicon sequence variant (16S); alpha diversity (within a sample: observed, Shannon, Simpson); beta diversity (between samples: Bray-Curtis, Jaccard); PCoA; PERMANOVA (do groups differ?).
- TMB: tumor mutational burden, somatic mutations per megabase of the sequenced region; values differ between panels, exomes and genomes.
- Mash distance, neighbor joining, UPGMA: genome distance from shared k-mers; tree-building methods.
"""


def _readme() -> str:
    try:
        from importlib.metadata import metadata

        text = metadata("multimerge-genomics").get_payload() or metadata("multimerge-genomics").get("Description") or ""
    except Exception:
        text = ""
    return text[:80_000]


@lru_cache(maxsize=1)
def system_text() -> str:
    """The whole system prompt: the rules, the glossary, the catalog, the README."""
    parts = [RULES, GLOSSARY, catalog.catalog_text()]
    readme = _readme()
    if readme:
        parts.append("# The README of MMGEN\n\n" + readme)
    return "\n\n".join(parts)
