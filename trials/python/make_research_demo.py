# Builds the research demo data of the web interface's examples (RNA-seq,
# CRISPR screens, microbiology): small public datasets, packed as one zip.
#
# Usage: python trials/python/make_research_demo.py <research_data folder> <output .zip> [RefSeq GFF]
#   e.g. python trials/python/make_research_demo.py research_data mmgen_demo_research_v3.zip
#        refseq_GRCh38p14_genomic.gff.gz
import gzip
import os
import shutil
import sys
import tempfile
import time
import zipfile

VERSION = "3"


def main(src, out_zip, refseq_gff=None):
    t0 = time.time()
    tmp = tempfile.mkdtemp(prefix="research_demo_")
    d = os.path.join(tmp, "research")

    def p(*n):
        return os.path.join(d, *n)

    # RNA-seq: the airway experiment (GSE52778): airway smooth-muscle cells of
    # four donors, untreated or treated with dexamethasone.
    os.makedirs(p("rnaseq"))
    air = os.path.join(src, "rnaseq", "airway")
    shutil.copy(os.path.join(air, "GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz"), p("rnaseq", "airway_counts.tsv.gz"))
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(air, "GSE52778_sample_metadata.tsv"), encoding="utf-8")]
    head = rows[0]
    col = {c: i for i, c in enumerate(head)}
    with open(p("rnaseq", "airway_samples.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("sample\tcell_line\tdex\n")
        for r in rows[1:]:
            if r[col["albuterol"]] == "untrt":
                fh.write(f"{r[col['sample_id']]}\t{r[col['cell_line']]}\t{r[col['dex']]}\n")
    with gzip.open(os.path.join(air, "Homo_sapiens.gene_info.gz"), "rt", encoding="utf-8") as fin, \
            open(p("rnaseq", "airway_gene_names.tsv"), "w", encoding="utf-8", newline="\n") as fout:
        fout.write("gene\tname\n")
        for line in fin:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            fout.write(f"{f[1]}\t{f[2]}\n")

    # CRISPR screens: the MAGeCK demos (BSD 3-clause) and the leukemia screen of
    # its tutorial.
    os.makedirs(p("crispr"))
    mg = os.path.join(src, "crispr", "mageck")
    shutil.copy(os.path.join(mg, "demo", "demo1", "sample.txt"), p("crispr", "demo1_counts.txt"))
    shutil.copy(os.path.join(mg, "demo", "demo2", "library.txt"), p("crispr", "demo2_library.txt"))
    for i in (1, 2):
        with open(os.path.join(mg, "demo", "demo2", f"test{i}.fastq"), "rb") as fin:
            with open(p("crispr", f"demo2_sample{i}.fastq.gz"), "wb") as fout:
                fout.write(gzip.compress(fin.read(), mtime=0))
    with open(os.path.join(mg, "leukemia.new.csv"), "rb") as fin:
        with open(p("crispr", "leukemia_counts.csv.gz"), "wb") as fout:
            fout.write(gzip.compress(fin.read(), mtime=0))
    shutil.copy(os.path.join(mg, "COPYING"), p("crispr", "MAGeCK_COPYING.txt"))

    # Microbiology: reference sketches, mock-community reads, the resistance
    # catalog and four isolates.
    import multimerge_genomics as mmg
    os.makedirs(p("microbio", "isolates"))

    def table(path):
        rows = [l.rstrip("\n").split("\t") for l in open(path, encoding="utf-8") if l.strip()]
        return [dict(zip(rows[0], r)) for r in rows[1:]]

    refs = table(os.path.join(src, "reference_genomes", "genomes.tsv"))
    paths = [[os.path.join(src, r["fasta"])] for r in refs]
    names = [r["organism"] for r in refs]
    host = os.path.join(src, "metagenome", "host_decoy")
    s = mmg.sketch_references(paths + [[os.path.join(host, "GRCh38.p14_chr21.fna.gz"), os.path.join(host, "GRCh38.p14_chrMT.fna.gz")]],
                              p("microbio", "references_k31_s1000.sig"), names=names + ["Homo sapiens (chr21 + chrM, host)"],
                              host=[False] * len(refs) + [True], scaled=1000)
    print(f"sketches: {len(s['genomes'])} genomes")
    zymo = ("Pseudomonas aeruginosa", "Escherichia coli str. K-12", "Salmonella enterica", "Limosilactobacillus fermentum",
            "Enterococcus faecalis", "Staphylococcus aureus subsp. aureus NCTC", "Listeria monocytogenes", "Bacillus subtilis",
            "Saccharomyces cerevisiae", "Cryptococcus neoformans")
    zi = [i for i, n in enumerate(names) if n.startswith(zymo)]
    s = mmg.sketch_references([paths[i] for i in zi] + [[os.path.join(host, "GRCh38.p14_chrMT.fna.gz")]],
                              p("microbio", "zymo_species_k31_s20.sig"), names=[names[i] for i in zi] + ["Homo sapiens (chrM, host)"],
                              host=[False] * len(zi) + [True], scaled=20)
    print(f"Zymo species sketches (scaled 20): {[g['name'] for g in s['genomes']]}")
    mock = os.path.join(src, "metagenome", "mock_community")
    for f in sorted(os.listdir(mock)):
        if f.endswith(".fastq.gz"):
            short = f.replace("_first1M", "_first30k")
            with gzip.open(os.path.join(mock, f), "rt", encoding="utf-8") as fin:
                lines = [next(fin) for _ in range(4 * 30_000)]
            with open(p("microbio", short), "wb") as fout:
                fout.write(gzip.compress("".join(lines).encode(), mtime=0))
    amr = os.path.join(src, "amr", "amrfinderplus", "2026-10-29.1")
    for f in ("AMR_CDS.fa", "ReferenceGeneCatalog.txt"):
        with open(os.path.join(amr, f), "rb") as fin, open(p("microbio", f + ".gz"), "wb") as fout:
            fout.write(gzip.compress(fin.read(), mtime=0))
    outbreak = table(os.path.join(src, "outbreak", "bacteria", "genomes.tsv"))
    wanted = {"KPNIH1": "kpneumoniae_KPNIH1_KPC-3", "ATCC BAA-2146": "kpneumoniae_ATCC_BAA-2146_NDM-1",
              "USA300_FPR3757": "saureus_USA300_MRSA"}
    for r in outbreak:
        for key, short in wanted.items():
            if key in r["strain"] or key in r["organism"]:
                shutil.copy(os.path.join(src, r["fasta"]), p("microbio", "isolates", short + ".fna.gz"))
    mssa = next(r for r in refs if "NCTC 8325" in r["organism"])
    shutil.copy(os.path.join(src, mssa["fasta"]), p("microbio", "isolates", "saureus_NCTC8325_MSSA.fna.gz"))
    print("isolates:", sorted(os.listdir(p("microbio", "isolates"))))

    # Outbreaks and surveillance: 20 SARS-CoV-2 genomes (one per Pango
    # lineage, renamed <lineage>_<WHO label>_<accession> so that the tree and
    # the tables say which is which), the reference and its annotation;
    # sketches of the K. pneumoniae and S. aureus genomes; a PCR panel.
    os.makedirs(p("outbreak"))
    sc = os.path.join(src, "outbreak", "sars_cov_2")
    lineage = {r["accession"]: r for r in table(os.path.join(sc, "sars_cov_2_20_lineages.metadata.tsv"))}
    with open(os.path.join(sc, "sars_cov_2_20_lineages.fasta"), encoding="utf-8") as fin, \
            open(p("outbreak", "sars_cov_2_20_lineages.fasta"), "w", encoding="utf-8", newline="\n") as fout:
        for line in fin:
            if line.startswith(">"):
                acc, _, desc = line[1:].strip().partition(" ")
                r = lineage[acc]
                who = r["who_label"] if r["who_label"] not in ("", "-") else ""
                fout.write(f">{r['pango_lineage']}{'_' + who if who else ''}_{acc} {desc}\n")
            else:
                fout.write(line)
    ref = next(r for r in refs if r["organism"].startswith("Severe acute respiratory syndrome coronavirus 2"))
    shutil.copy(os.path.join(src, ref["fasta"]), p("outbreak", "NC_045512.2.fna.gz"))
    shutil.copy(os.path.join(src, ref["gff3"]), p("outbreak", "NC_045512.2.gff3.gz"))
    bact = table(os.path.join(src, "outbreak", "bacteria", "genomes.tsv"))

    def strain_name(r):
        # "ATCC 700721; MGH 78578" -> "MGH 78578"; a strain named after its ST
        # ("ST307") takes the assembly accession instead.
        strain = r["strain"].split(";")[-1].strip()
        if strain.startswith("ST"):
            strain = r["assembly_accession"]
        return f"{strain} (ST{r['MLST_ST']})"

    for species, short, scaled in (("Klebsiella", "klebsiella_pneumoniae", 100), ("Staphylococcus", "staphylococcus_aureus", 100)):
        rows = [r for r in bact if r["organism"].startswith(species)]
        s = mmg.sketch_references([[os.path.join(src, r["fasta"])] for r in rows],
                                  p("outbreak", f"{short}_{len(rows)}_k31_s{scaled}.sig"),
                                  names=[strain_name(r) for r in rows], scaled=scaled)
        print(f"{species}: {len(s['genomes'])} sketches")

    # 16S amplicons: the mothur MiSeq SOP (Kozich et al. 2013): 19 mouse
    # stool samples (V4, 2 x 250) and the Mock community. To keep the pack
    # small, the first 2,000 pairs of each mouse sample; the Mock whole.
    os.makedirs(p("16s"))
    sop = os.path.join(src, "16s", "mothur_miseq_sop", "MiSeq_SOP")
    sheet = []
    for line in open(os.path.join(sop, "stability.files"), encoding="utf-8"):
        name, r1, r2 = line.split()
        keep = None if name == "Mock" else 2000
        for f in (r1, r2):
            with open(os.path.join(sop, f), encoding="utf-8") as fin:
                lines = fin.readlines() if keep is None else [next(fin) for _ in range(4 * keep)]
            with open(p("16s", f + ".gz"), "wb") as fout:
                fout.write(gzip.compress("".join(lines).encode(), mtime=0))
        sheet.append(f"{name}\t{r1}.gz\t{r2}.gz\n")
    with open(p("16s", "samples.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.writelines(sheet)
    with open(p("16s", "samples_design.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        for line in open(os.path.join(sop, "mouse.time.design"), encoding="utf-8"):
            fh.write(line.replace("\r", ""))
    shutil.copy(os.path.join(sop, "HMP_MOCK.v35.fasta"), p("16s", "HMP_MOCK.v35.fasta"))
    rdp = os.path.join(src, "16s", "reference", "rdp19_dada2", "rdp_19_toGenus_trainset.fa.gz")
    shutil.copy(rdp, p("16s", "rdp_19_toGenus_trainset.fa.gz"))
    t = mmg.train_classifier(rdp, p("16s", "rdp19_v4_classifier.model"),
                             primer_forward="GTGYCAGCMGCCGCGGTAA", primer_reverse="GGACTACNVGGGTWTCTAAT")
    print("16S classifier:", {k: t[k] for k in t if k != "text" and not isinstance(t[k], (list, dict))})
    # The ASVs of these reads, made once here, so that the taxonomy and
    # diversity examples do not depend on running the ASV example first; and
    # the table without the Mock, whose design has no time point.
    here = os.getcwd()
    os.chdir(p("16s"))
    try:
        a = mmg.amplicon_asvs("samples.tsv", "sop_asvs.fasta", "sop_asv_table.tsv")
    finally:
        os.chdir(here)
    print(f"16S: {len(a['asvs'])} ASVs")
    with open(p("16s", "sop_asv_table.tsv"), encoding="utf-8") as fin:
        rows = [l.rstrip("\n").split("\t") for l in fin]
    mock = rows[0].index("Mock")
    with open(p("16s", "sop_asv_table_mouse.tsv"), "w", encoding="utf-8", newline="\n") as fout:
        for i, r in enumerate(rows):
            r = r[:mock] + r[mock + 1:]
            if i == 0 or any(x != "0" for x in r[1:]):
                fout.write("\t".join(r) + "\n")

    # Tumor mutational burden: the SEQC2 somatic truth set of the HCC1395
    # breast-cancer cell line (Fang et al. 2021, Nat Biotechnol 39:1151),
    # chromosome 17: SNVs and indels in its high-confidence region; and the
    # coding sequences (CDS) of chromosome 17 from the RefSeq annotation.
    os.makedirs(p("tumor"))
    seqc = os.path.join(src, "oncology", "seqc2_hcc1395")
    meta, columns, recs = [], None, []
    for f in ("high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz", "high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz"):
        with gzip.open(os.path.join(seqc, f), "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("##"):
                    if line not in meta:
                        meta.append(line)
                elif line.startswith("#"):
                    columns = columns or line
                elif line.startswith("chr17\t"):
                    recs.append(line)
    header = meta + [columns]
    tmp_vcf = p("tumor", "unsorted.vcf")
    with open(tmp_vcf, "w", encoding="utf-8", newline="\n") as fh:
        fh.writelines(header + recs)
    mmg.sort_genomic([tmp_vcf], p("tumor", "HCC1395_somatic_chr17.vcf.gz"), chrom_order="natural")
    os.remove(tmp_vcf)
    with open(os.path.join(seqc, "High-Confidence_Regions_v1.2.bed"), encoding="utf-8") as fin, \
            open(p("tumor", "HCC1395_high_confidence_chr17.bed"), "w", encoding="utf-8", newline="\n") as fout:
        for line in fin:
            if line.startswith("chr17\t"):
                fout.write(line)
    if refseq_gff:
        with gzip.open(refseq_gff, "rt", encoding="utf-8") as fin, \
                gzip.open(p("tumor", "refseq_chr17_cds.gff3.gz"), "wt", encoding="utf-8", newline="\n", compresslevel=9) as fout:
            fout.write("##gff-version 3\n")
            for line in fin:
                if line.startswith("NC_000017.11\t"):
                    kind = line.split("\t", 3)[2]
                    if kind in ("region", "gene", "CDS"):
                        fout.write(line)
    print(f"tumor: {len(recs)} somatic variants on chr17")
    shutil.copy(os.path.join(os.path.dirname(os.path.abspath(__file__)), "pcr_demo_panel.tsv"), p("outbreak", "pcr_demo_panel.tsv"))

    with open(p("README.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"""# MMGEN research demo data (version {VERSION})

Small public datasets for the examples of the web interface.

* `rnaseq/airway_counts.tsv.gz`: raw read counts per gene (NCBI GeneID) of the
  airway experiment, GEO GSE52778 (Himes et al. 2014, PLoS One 9:e99625):
  airway smooth-muscle cells of four donors, untreated, treated with
  dexamethasone, albuterol, or both; counts computed by NCBI (GRCh38.p13).
  `airway_samples.tsv`: the 8 samples without albuterol (cell line, dex:
  trt/untrt). `airway_gene_names.tsv`: GeneID to symbol, from NCBI gene_info.
* `crispr/`: the demo data of MAGeCK 0.5.9.5 (Li et al. 2014, Genome Biology
  15:554; BSD 3-clause license, `MAGeCK_COPYING.txt`): `demo1_counts.txt` (a
  count table), `demo2_library.txt` + `demo2_sample1/2.fastq.gz` (reads to
  count), and `leukemia_counts.csv.gz`, the HL60 / KBM7 screen of MAGeCK's
  tutorial (73,151 guides).
* `microbio/references_k31_s1000.sig`: FracMinHash sketches (k=31, scaled
  1000) of the NCBI RefSeq reference genomes of 25 bacteria, 2 yeasts and 16
  viruses, and of human chr21 + chrM as host; `zymo_species_k31_s20.sig`: the
  ten species of the ZymoBIOMICS standard, denser (scaled 20), for read
  classification. `ERR2984773_zymo_even_*` and `ERR2935805_zymo_log_*`: the
  first 30,000 read pairs of two public runs (ENA) of the ZymoBIOMICS
  microbial community standards (EVEN and LOG).
* `microbio/AMR_CDS.fa.gz`, `ReferenceGeneCatalog.txt.gz`: the NCBI
  AMRFinderPlus resistance-gene catalog, release 2026-10-29.1 (public domain).
  `isolates/`: four public genome assemblies (NCBI RefSeq): K. pneumoniae
  KPNIH1 (blaKPC-3), K. pneumoniae ATCC BAA-2146 (blaNDM-1), S. aureus
  USA300 (MRSA, mecA) and S. aureus NCTC 8325 (no mecA).
* `outbreak/sars_cov_2_20_lineages.fasta`: 20 complete SARS-CoV-2 genomes
  from NCBI Virus (GenBank), one per Pango lineage, from A (January 2020) to
  KP.3.1.1; each record is renamed `<lineage>_<WHO label>_<accession>`.
  `NC_045512.2.fna.gz` / `.gff3.gz`: the reference genome (Wuhan-Hu-1,
  NCBI RefSeq GCF_009858895.2) and its annotation.
* `outbreak/klebsiella_pneumoniae_13_k31_s100.sig`,
  `staphylococcus_aureus_4_k31_s100.sig`: FracMinHash sketches (k=31,
  scaled 100) of public genome assemblies (NCBI RefSeq) of 13 K. pneumoniae
  and 4 S. aureus, named by strain and sequence type.
* `outbreak/pcr_demo_panel.tsv`: published primer pairs, each with its
  source: the CDC 2019-nCoV N1 and N2 assays (Lu et al. 2020, Emerg Infect
  Dis 26:1654), 16S V4 (515F-Y / 806RB) and near-full-length 16S (27F /
  1492R), and mecA (Murakami et al. 1991, J Clin Microbiol 29:2240).
* `16s/`: the mothur MiSeq SOP data (Kozich et al. 2013, Appl Environ
  Microbiol 79:5112; https://mothur.org/wiki/miseq_sop/): 16S V4 reads
  (2 x 250) of 19 mouse stool samples, the first 2,000 read pairs of each,
  and the whole Mock community; `samples.tsv` (name, R1, R2),
  `samples_design.tsv` (early / late after weaning), `HMP_MOCK.v35.fasta`
  (the true sequences of the Mock's strains). `rdp19_v4_classifier.model`:
  the RDP naive Bayesian classifier trained by MMGEN on the V4 region of the
  RDP training set 19 (Cole et al. 2014; DADA2 format by B. Callahan,
  doi:10.5281/zenodo.14168771; CC BY 4.0).
* `tumor/`: the SEQC2 somatic reference call set of the HCC1395 breast-cancer
  cell line (Fang et al. 2021, Nat Biotechnol 39:1151; v1.2.1, NCBI FTP),
  chromosome 17: the high-confidence SNVs and indels and the high-confidence
  region; `refseq_chr17_cds.gff3.gz`: the genes and coding sequences of
  chromosome 17 from the NCBI RefSeq GRCh38.p14 annotation.
""")
    with open(p("VERSION"), "w") as fh:
        fh.write(VERSION)
    with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_DEFLATED) as z:
        for root, _dirs, files in os.walk(d):
            for f in files:
                full = os.path.join(root, f)
                z.write(full, os.path.relpath(full, tmp).replace(os.sep, "/"))
    shutil.rmtree(tmp, ignore_errors=True)
    print(f"{out_zip}: {os.path.getsize(out_zip) / 1e6:.1f} MB in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        sys.exit("Usage: python make_research_demo.py <research_data folder> <output .zip> [RefSeq GRCh38 GFF, for the CDS of chr17]")
    main(*sys.argv[1:])
