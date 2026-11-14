# Builds the human-genome demo data of the web interface's examples: small
# cuts of public data (Genome in a Bottle, ClinVar, RefSeq, PGS Catalog), made
# with the engine itself, packed as one zip.
#
# Usage: python trials/python/make_human_demo.py <human_genomes folder> <RefSeq GFF> <output .zip>
#   e.g. python trials/python/make_human_demo.py human_genomes \
#        refseq_GRCh38p14_genomic.gff.gz mmgen_demo_human_v3.zip
#
# What goes in (all GRCh38):
#  - the seven GIAB genomes (HG001-HG007, benchmark v4.2.1), cut to chr22 plus
#    the regions of the example gene panels (+-5 kb) and of the PharmCAT
#    positions (+-50 bp), with their callable-regions BEDs cut the same way;
#  - chrX/chrY cuts of HG001 (v3.3.2) and HG002 (v5.0q), for genetic sex;
#  - ClinVar records in the same regions; the panel genes from the RefSeq GFF;
#  - the PGS Catalog scoring files that were downloaded.
import gzip
import json
import os
import shutil
import sys
import tempfile
import time
import zipfile

import multimerge_genomics as mmg

VERSION = "3"
PANELS = {
    "Hereditary breast and ovarian cancer": ["BRCA1", "BRCA2", "PALB2", "TP53", "PTEN", "CDH1", "STK11", "ATM", "CHEK2", "BARD1", "BRIP1", "RAD51C", "RAD51D"],
    "Lynch syndrome (hereditary colorectal cancer)": ["MLH1", "MSH2", "MSH6", "PMS2", "EPCAM"],
    "Familial hypercholesterolemia": ["LDLR", "APOB", "PCSK9", "LDLRAP1"],
    "Cardiomyopathy": ["MYH7", "MYBPC3", "TNNT2", "TNNI3", "TPM1", "ACTC1", "MYL2", "MYL3", "LMNA", "TTN", "PLN", "DSP", "PKP2"],
    "Inherited arrhythmia": ["KCNQ1", "KCNH2", "SCN5A", "RYR2", "CASQ2"],
    "Bleeding and clotting": ["F2", "F5", "F8", "F9", "F11", "VWF"],
    "Carrier screening (common recessive conditions)": ["CFTR", "HBB", "HBA1", "HBA2", "HEXA", "GBA1", "ASPA", "FANCC", "BLM", "PAH", "SERPINA1", "GJB2", "MEFV"],
    "Pharmacogenes": ["CYP2C19", "CYP2C9", "CYP2D6", "CYP3A5", "CYP4F2", "DPYD", "TPMT", "NUDT15", "SLCO1B1", "UGT1A1", "VKORC1", "IFNL3", "NAT2", "CYP2B6", "RYR1", "CACNA1S"],
}
PEOPLE = [
    ("HG001", "HG001_NA12878", "noinconsistent"), ("HG002", "HG002_NA24385_son", "noinconsistent"),
    ("HG003", "HG003_NA24149_father", "noinconsistent"), ("HG004", "HG004_NA24143_mother", "noinconsistent"),
    ("HG005", "HG005_NA24631_son", ""), ("HG006", "HG006_NA24694_father", ""), ("HG007", "HG007_NA24695_mother", ""),
]


def main(src, gff, out_zip):
    t0 = time.time()
    tmp = tempfile.mkdtemp(prefix="human_demo_")
    d = os.path.join(tmp, "human")
    os.makedirs(d)

    def p(*n):
        return os.path.join(d, *n)

    # 1. The regions: chr22, the panel genes and the PharmCAT positions.
    genes = sorted({g for v in PANELS.values() for g in v})
    gr = mmg.gene_regions(gff, genes, output=p("panel_genes.bed"), padding=0, naming="chr")
    print(f"panel genes: {len(gr['regions'])} found, not found: {gr['not_found']}")
    regions = [("chr22", 0, 50_818_468)]
    regions += [(r["chrom"], max(0, r["start"] - 5000), r["end"] + 5000) for r in gr["regions"] if r["primary"]]
    positions = None
    for root, _dirs, files in os.walk(os.path.join(os.path.dirname(src), "tools", "pharmcat")):
        for f in files:
            if f.startswith("pharmcat_positions_") and f.endswith(".vcf") and "uniallelic" not in f:
                positions = os.path.join(root, f)
    if positions:
        for line in open(positions, encoding="utf-8"):
            if not line.startswith("#"):
                c, s = line.split("\t")[:2]
                regions.append((c, max(0, int(s) - 50), int(s) + 50))
    # The variants of the polygenic scores (harmonized GRCh38 positions).
    pgs_dir = os.path.join(src, "pgs_catalog")
    pgs_files = sorted(os.path.join(r, f) for r, _d, fs in os.walk(pgs_dir) for f in fs
                       if f.endswith("_hmPOS_GRCh38.txt.gz") or f.endswith("_metadata_scores.csv"))
    pgs_sites = set()
    for f in pgs_files:
        if f.endswith("_hmPOS_GRCh38.txt.gz"):
            with gzip.open(f, "rt", encoding="utf-8") as fh:
                cols = None
                for line in fh:
                    if line.startswith("#"):
                        continue
                    row = line.rstrip("\n").split("\t")
                    if cols is None:
                        cols = row
                        continue
                    rec = dict(zip(cols, row))
                    if rec.get("hm_chr") and rec.get("hm_pos", "").isdigit():
                        pgs_sites.add(("chr" + rec["hm_chr"], int(rec["hm_pos"])))
    regions += [(c, max(0, s - 50), s + 50) for c, s in pgs_sites]
    print(f"{len(pgs_sites)} polygenic-score sites")
    with open(p("demo_regions_unsorted.bed"), "w", newline="\n") as f:
        f.writelines(f"{c}\t{a}\t{b}\n" for c, a, b in regions)
    mmg.bed_merge([p("demo_regions_unsorted.bed")], p("demo_regions.bed"))
    os.remove(p("demo_regions_unsorted.bed"))
    merged = [l.split("\t")[:3] for l in open(p("demo_regions.bed"))]
    query = [f"{c}:{int(a) + 1}-{b.strip()}" for c, a, b in merged]
    autosomes = [q for q in query if not q.startswith(("chrX", "chrY", "chrM"))]  # GIAB v4.2.1 is chr1-22
    print(f"{len(merged)} demo regions")

    # 2. The seven GIAB genomes and their callable regions, cut.
    os.makedirs(p("giab"))
    for sid, folder, kind in PEOPLE:
        base = os.path.join(src, "giab", folder, "v4.2.1", f"{sid}_GRCh38_1_22_v4.2.1_benchmark")
        n = mmg.query_region(base + ".vcf.gz", autosomes, output=p("giab", f"{sid}.vcf.gz"))
        mmg.index_file(p("giab", f"{sid}.vcf.gz"))
        bed = next(b for b in (base + "_noinconsistent.bed", base + ".bed") if os.path.exists(b))
        mmg.bed_intersect([bed], [p("demo_regions.bed")], p("giab", f"{sid}.callable.bed"))
        print(f"{sid}: {n:,} records")

    # 2b. Population allele frequencies (1000 Genomes, GRCh38 sites file), only
    # at the sites of the seven cut genomes: what the paternity index needs.
    sites_dir = os.path.join(src, "1000genomes", "phase3_grch38_sites_20190312")
    sites = next((os.path.join(sites_dir, f) for f in os.listdir(sites_dir) if f.endswith(".vcf.gz")), None) if os.path.isdir(sites_dir) else None
    if sites:
        wanted_sites = set()
        for sid, _f, _k in PEOPLE:
            with gzip.open(p("giab", f"{sid}.vcf.gz"), "rt") as fh:
                for line in fh:
                    if not line.startswith("#"):
                        c, s = line.split("\t", 2)[:2]
                        wanted_sites.add((c.replace("chr", ""), s))
        tmp_af = p("af_all.vcf.gz")
        mmg.query_region(sites, [q[3:] for q in autosomes], output=tmp_af)
        kept = 0
        with gzip.open(tmp_af, "rt") as fin, open(p("af_tmp.vcf"), "w", newline="\n") as fout:
            for line in fin:
                if line.startswith("#"):
                    fout.write(line)
                elif (line.split("\t", 2)[0], line.split("\t", 2)[1]) in wanted_sites:
                    fout.write(line)
                    kept += 1
        mmg.convert_file(p("af_tmp.vcf"), p("af_1000g_demo.vcf.gz"))
        mmg.index_file(p("af_1000g_demo.vcf.gz"))
        os.remove(tmp_af)
        os.remove(p("af_tmp.vcf"))
        print(f"1000 Genomes frequencies: {kept:,} sites")

        # The same frequencies as panels: genome-wide, thinned (ancestry of a
        # whole genome), and at the polygenic-score sites (where a score falls).
        groups = ["AFR", "AMR", "EAS", "EUR", "SAS"]
        r = mmg.panel_from_frequencies([sites], groups, [g + "_AF" for g in groups], p("ancestry_panel_genomewide_1000g.tsv.gz"), thin_bp=20_000)
        print(f"genome-wide ancestry panel: {r['written']:,} sites")
        pgs_regions = sorted({f"{c[3:]}:{s}-{s}" for c, s in pgs_sites})
        mmg.query_region(sites, pgs_regions, output=p("pgs_sites.vcf.gz"))
        r = mmg.panel_from_frequencies([p("pgs_sites.vcf.gz")], groups, [g + "_AF" for g in groups], p("pgs_panel_1000g.tsv.gz"), min_maf=0.0, keep_indels=True)
        os.remove(p("pgs_sites.vcf.gz"))
        print(f"polygenic-score panel: {r['written']:,} sites")

    # 2c. A chr22 ancestry panel with exact frequencies, from the 1000 Genomes
    # 30x genotypes: the unrelated people (the children of the trios left out).
    panel_vcf = os.path.join(src, "1000genomes", "nygc_30x_phased_panel", "1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
    meta = os.path.join(src, "1000genomes", "nygc_30x_metadata")
    if os.path.exists(panel_vcf):
        children = []
        for line in open(os.path.join(meta, "1kGP.3202_samples.pedigree_info.txt"), encoding="utf-8"):
            f = line.split()
            if len(f) >= 3 and f[0] != "sampleID" and (f[1] != "0" or f[2] != "0"):
                children.append(f[0])
        r = mmg.build_panel([panel_vcf], os.path.join(meta, "20130606_g1k_3202_samples_ped_population.txt"),
                            p("ancestry_panel_chr22_1000g.tsv.gz"), "SampleID", "Superpopulation", exclude=children)
        print(f"chr22 ancestry panel: {r['written']:,} sites, {sum(r['samples_per_group'])} samples")

    # 3. Sex chromosomes: HG001 v3.3.2 (chrX) and HG002 v5.0q (chrX + chrY).
    sex_regions = ["chrX:2781480-155701382", "chrY:2781480-56887902"]
    hg001x = os.path.join(src, "giab", "HG001_NA12878", "v3.3.2_with_chrX")
    vcf = next((os.path.join(hg001x, f) for f in os.listdir(hg001x) if f.endswith(".vcf.gz")), None) if os.path.isdir(hg001x) else None
    try:
        if vcf:
            n = mmg.query_region(vcf, sex_regions[:1], output=p("giab", "HG001.chrX.vcf.gz"))
            print(f"HG001 chrX: {n:,} records")
    except RuntimeError as err:
        print(f"HG001 chrX skipped: {err}")
    v5 = os.path.join(src, "giab", "HG002_NA24385_son", "v5.0q", "HG002_GRCh38_v5.0q_smvar.vcf.gz")
    try:
        if os.path.exists(v5):
            # v5.0q: only FILTER "." is a good call (HET1/2, GAP1/2, DIPX, DIPY are not).
            n = mmg.query_region(v5, sex_regions, output=p("giab", "HG002.v5.0q.chrXY_all.vcf.gz"))
            r = mmg.vcf_filter(p("giab", "HG002.v5.0q.chrXY_all.vcf.gz"), p("giab", "HG002.chrXY.vcf.gz"), expression='FILTER == "."')
            os.remove(p("giab", "HG002.v5.0q.chrXY_all.vcf.gz"))
            print(f"HG002 chrX/Y (v5.0q, FILTER .): {r['kept']:,} of {n:,}")
    except RuntimeError as err:
        print(f"HG002 chrX/Y skipped: {err}")

    # 4. ClinVar in the same regions (ClinVar names chromosomes without "chr").
    clinvar = next(os.path.join(src, "clinvar", f) for f in os.listdir(os.path.join(src, "clinvar")) if f.endswith(".vcf.gz") and f.startswith("clinvar_"))
    n = mmg.query_region(clinvar, [q[3:] for q in query], output=p("clinvar_demo.vcf.gz"))
    mmg.index_file(p("clinvar_demo.vcf.gz"))
    print(f"ClinVar: {n:,} records")

    # 5. The panel genes as a small GFF3 (with the accession-to-chromosome lines).
    wanted = {g.upper() for g in genes}
    keep_seq = set()
    lines = []
    with gzip.open(gff, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##gff-version"):
                    lines.append(line)
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            if f[2] == "region" and ";chromosome=" in f[8] and "genome=chromosome" in f[8]:
                lines.append(line)
            elif f[2] == "gene":
                attrs = dict(kv.split("=", 1) for kv in f[8].strip().split(";") if "=" in kv)
                if attrs.get("Name", "").upper() in wanted or attrs.get("gene", "").upper() in wanted:
                    lines.append(line)
                    keep_seq.add(f[0])
    with gzip.open(p("panel_genes.gff3.gz"), "wt", encoding="utf-8", newline="\n") as fh:
        fh.writelines(lines)
    print(f"GFF: {len(lines)} lines")

    # 6. The PGS Catalog scores (harmonized, GRCh38) and the panels.
    os.makedirs(p("pgs"))
    for f in pgs_files:
        shutil.copy(f, p("pgs", os.path.basename(f)))
    print(f"PGS files: {len(pgs_files)}")
    with open(p("panels.json"), "w", encoding="utf-8") as fh:
        json.dump(PANELS, fh, indent=1)

    with open(p("README.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"""# MMGEN human-genome demo data (version {VERSION})

Small cuts of public data, for the examples of the web interface. GRCh38.
Experimental software, for teaching and research only.

* `giab/HG00X.vcf.gz` (+ `.tbi`) and `giab/HG00X.callable.bed`: the Genome in a Bottle
  benchmark small-variant calls v4.2.1 (NIST; public data) of HG001 (NA12878),
  the Ashkenazi trio HG002 (son) / HG003 (father) / HG004 (mother) and the Han
  Chinese trio HG005 (son) / HG006 (father) / HG007 (mother), cut to chr22 and
  to the regions in `demo_regions.bed`. A site absent from a VCF is
  homozygous reference only inside its `.callable.bed`.
* `giab/HG001.chrX.vcf.gz` (v3.3.2) and `giab/HG002.chrXY.vcf.gz` (v5.0q, FILTER "."
  only): the sex chromosomes outside the pseudoautosomal regions.
* `clinvar_demo.vcf.gz`: ClinVar (NCBI; public domain) records in the same regions.
* `af_1000g_demo.vcf.gz`: population allele frequencies of the 1000 Genomes Project
  (GRCh38 sites release of 2019-03-12; INFO AF, AFR_AF, AMR_AF, EAS_AF, EUR_AF,
  SAS_AF, rounded to 2 decimals) at the sites of the cut genomes.
* `ancestry_panel_chr22_1000g.tsv.gz`: allele frequencies of the five 1000 Genomes
  super-populations (AFR, AMR, EAS, EUR, SAS) on chr22, computed from the 30x
  genotypes of the unrelated people; `ancestry_panel_genomewide_1000g.tsv.gz`: the
  same from the sites release, one SNV per 20 kb over the genome (for whole
  genomes); `pgs_panel_1000g.tsv.gz`: the frequencies at the polygenic-score sites.
* `panel_genes.gff3.gz`, `panel_genes.bed`: the genes of the example panels
  (`panels.json`), from the NCBI RefSeq GRCh38.p14 annotation.
* `pgs/`: PGS Catalog scoring files (cite each score's publication).

Source versions and licenses: see the MANIFEST.md of the original downloads.
""")
    with open(p("VERSION"), "w") as fh:
        fh.write(VERSION)

    with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_STORED) as z:
        for root, _dirs, files in os.walk(d):
            for f in files:
                full = os.path.join(root, f)
                z.write(full, os.path.relpath(full, tmp).replace(os.sep, "/"))
    shutil.rmtree(tmp, ignore_errors=True)
    print(f"{out_zip}: {os.path.getsize(out_zip) / 1e6:.1f} MB in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__ if False else "Usage: python make_human_demo.py <human_genomes folder> <RefSeq GFF> <output .zip>")
    main(*sys.argv[1:])
