"""A small teaching cohort for genotype_association (MMGEN).

`write_demo(folder)` writes, with the same bytes on every run (standard library
only; the .vcf.gz is BGZF with MTIME 0, so index_file can index it):

  cohort.vcf.gz         120 people x 40 SNVs (real GRCh38 positions and alleles of
                        common dbSNP 151 SNVs in five candidate genes: ABCB1, BDNF,
                        COMT, IL6, OPRM1); genotypes SIMULATED
  phenotypes.tsv        122 rows: recovered (yes/no), recovery_score (0-100), age,
                        sex, center (hospital A/B/C) and three ancestry proportions
                        (anc_afr, anc_eur, anc_eas, summing to 1)
  candidate_genes.bed   the five genes (RefSeq GRCh38.p14 gene spans)
  README.md             what the demo shows and the numbers to expect

EVERYTHING IS SIMULATED. The positions and alleles are real, but the genotypes and
outcomes are drawn by this program; the "associations" are planted by it and say
nothing about what these variants do in real people.

What is planted:
  * a true effect: OPRM1 rs12190259 (ALT C) raises recovery_score by 6 points per
    copy and the odds of recovering by e^1.0 = 2.7 per copy, in every ancestry;
  * a confounded variant: BDNF rs2030324 (ALT G) has no effect at all, but its
    frequency differs by ancestry (75% in the AFR-like component, 15% EUR-like,
    30% EAS-like) and ancestry (anc_afr) changes the outcome: unadjusted it looks
    associated; adjusted for ancestry it is not;
  * a rare variant: IL6 rs2069830 (ALT T) carried by 3 people who all recovered:
    complete separation, so the plain logistic fit fails and Firth's method is used;
  * missing data on purpose: 2 genotyped people have no phenotype row, 2 rows have
    no genotypes, 1 person has no age, 2 have no recovery_score, and about 1% of
    the genotypes are missing (./.).

Usage: python -m multimerge_genomics.demos.association <folder>
"""
import math
import os
import random
import struct
import sys
import zlib

DEMO_VERSION = 1

# (gene, chrom, pos, rsID, ref, alt, dbSNP 151 common allele frequency of ALT)
VARIANTS = [
    ("ABCB1", "chr7", 87503600, "rs1055302", "C", "T", 0.2188),
    ("ABCB1", "chr7", 87522729, "rs4148742", "A", "G", 0.1532),
    ("ABCB1", "chr7", 87523474, "rs2141849", "C", "A", 0.1637),
    ("ABCB1", "chr7", 87533700, "rs11760837", "T", "C", 0.1454),
    ("ABCB1", "chr7", 87549310, "rs2235013", "C", "T", 0.4651),
    ("ABCB1", "chr7", 87559294, "rs28381873", "C", "T", 0.1056),
    ("ABCB1", "chr7", 87564281, "rs4148734", "G", "A", 0.1751),
    ("ABCB1", "chr7", 87577592, "rs75232677", "T", "C", 0.1156),
    ("ABCB1", "chr7", 87627286, "rs34800935", "T", "C", 0.252),
    ("ABCB1", "chr7", 87640181, "rs2214101", "T", "C", 0.1052),
    ("ABCB1", "chr7", 87711219, "rs28381724", "G", "A", 0.1516),
    ("BDNF", "chr11", 27656036, "rs11030099", "C", "A", 0.229),
    ("BDNF", "chr11", 27660049, "rs11030102", "C", "G", 0.1142),
    ("BDNF", "chr11", 27673288, "rs11030107", "A", "G", 0.1144),
    ("BDNF", "chr11", 27684012, "rs4633417", "C", "T", 0.4301),
    ("BDNF", "chr11", 27685862, "rs11500197", "G", "A", 0.113),
    ("BDNF", "chr11", 27690834, "rs1013402", "A", "G", 0.2131),
    ("BDNF", "chr11", 27705368, "rs2030324", "A", "G", 0.4341),
    ("BDNF", "chr11", 27720977, "rs58936914", "C", "G", 0.1102),
    ("COMT", "chr22", 19942598, "rs737865", "A", "G", 0.2266),
    ("COMT", "chr22", 19944690, "rs34798238", "C", "T", 0.1054),
    ("COMT", "chr22", 19947625, "rs9332334", "T", "C", 0.1128),
    ("COMT", "chr22", 19951236, "rs112080665", "A", "G", 0.1322),
    ("COMT", "chr22", 19957200, "rs5993888", "C", "T", 0.4499),
    ("COMT", "chr22", 19961490, "rs165722", "C", "T", 0.4279),
    ("COMT", "chr22", 19963748, "rs4680", "G", "A", 0.3692),
    ("COMT", "chr22", 19969258, "rs165599", "G", "A", 0.4908),
    ("IL6", "chr7", 22727518, "rs2069830", "C", "T", 0.02576),
    ("IL6", "chr7", 22728600, "rs1524107", "C", "T", 0.3083),
    ("IL6", "chr7", 22728630, "rs2066992", "G", "T", 0.3081),
    ("IL6", "chr7", 22728953, "rs2069840", "C", "G", 0.1855),
    ("IL6", "chr7", 22730391, "rs2069844", "C", "A", 0.1008),
    ("OPRM1", "chr6", 154011812, "rs1294093", "A", "G", 0.4022),
    ("OPRM1", "chr6", 154036852, "rs12190259", "A", "C", 0.138),
    ("OPRM1", "chr6", 154083178, "rs3778156", "A", "G", 0.1112),
    ("OPRM1", "chr6", 154135053, "rs686633", "C", "T", 0.1138),
    ("OPRM1", "chr6", 154150547, "rs9383692", "A", "G", 0.4952),
    ("OPRM1", "chr6", 154184549, "rs12203560", "A", "G", 0.1607),
    ("OPRM1", "chr6", 154213278, "rs790252", "G", "A", 0.2113),
    ("OPRM1", "chr6", 154246532, "rs9479798", "T", "G", 0.2839),
]
GENES = [  # RefSeq GRCh38.p14 gene spans, BED (0-based start)
    ("chr6", 154010495, 154246867, "OPRM1"),
    ("chr7", 22727199, 22731998, "IL6"),
    ("chr7", 87503016, 87713295, "ABCB1"),
    ("chr11", 27654892, 27722030, "BDNF"),
    ("chr22", 19941771, 19969975, "COMT"),
]
CONTIGS = [("chr6", 170805979), ("chr7", 159345973), ("chr11", 135086622), ("chr22", 50818468)]
CAUSAL = "rs12190259"
CONFOUNDED = "rs2030324"
RARE = "rs2069830"
ORDER = {"chr6": 6, "chr7": 7, "chr11": 11, "chr22": 22}


def _normal(rng):
    u1 = 1.0 - rng.random()
    u2 = rng.random()
    return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)


def _bgzf(data):
    """BGZF (blocked gzip, MTIME 0), any gzip reader opens it; ends with the EOF block."""
    out = bytearray()
    for i in range(0, len(data), 65280):
        bloco = data[i:i + 65280]
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        comp = c.compress(bloco) + c.flush()
        out += struct.pack("<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, len(comp) + 25)
        out += comp
        out += struct.pack("<II", zlib.crc32(bloco) & 0xFFFFFFFF, len(bloco))
    out += bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
    return bytes(out)


def simulate():
    rng = random.Random(2026)
    n = 120
    people = [f"PT{i:03d}" for i in range(1, n + 1)]
    # ancestry: a main component (45% EUR-like, 35% AFR-like, 20% EAS-like) with admixture
    anc = []
    for _ in range(n):
        u = rng.random()
        main = 1 if u < 0.45 else (0 if u < 0.80 else 2)  # 0 AFR, 1 EUR, 2 EAS
        p_main = 0.55 + 0.45 * rng.random()
        rest = 1.0 - p_main
        s = rng.random()
        others = [k for k in range(3) if k != main]
        p = [0.0, 0.0, 0.0]
        p[main] = p_main
        p[others[0]] = rest * s
        p[others[1]] = rest * (1 - s)
        a = [round(x, 3) for x in p]
        a[2] = round(1.0 - a[0] - a[1], 3)  # the three columns sum to 1.000
        anc.append(a)
    age = [max(18, min(90, round(55 + 12 * _normal(rng)))) for _ in range(n)]
    sex = ["F" if rng.random() < 0.5 else "M" for _ in range(n)]
    # hospital
    center = []
    for _ in anc:
        u = rng.random()
        center.append("A" if u < 0.45 else ("B" if u < 0.75 else "C"))
    # genotypes: ALT frequency of each person = sum of ancestry proportion x frequency
    geno = {}
    for (_, _, _, rs, _, _, f) in VARIANTS:
        if rs == RARE:
            geno[rs] = [0] * n
            continue
        if rs == CONFOUNDED:
            fk = (0.75, 0.15, 0.30)
        elif rs == CAUSAL:
            fk = (0.25, 0.25, 0.25)
        else:
            fk = tuple(min(0.95, max(0.02, f * (0.6 + 0.8 * rng.random()))) for _ in range(3))
        col = []
        for a in anc:
            fi = sum(x * y for x, y in zip(a, fk))
            col.append((rng.random() < fi) + (rng.random() < fi))
        geno[rs] = col
    g1 = geno[CAUSAL]
    score, recovered = [], []
    for i in range(n):
        cb = 1.0 if center[i] == "B" else 0.0
        cc = 1.0 if center[i] == "C" else 0.0
        s = 60 + 6 * g1[i] - 0.25 * (age[i] - 55) + 4 * cb - 3 * cc + 25 * (anc[i][0] - 0.35) + 10 * _normal(rng)
        score.append(round(max(0.0, min(100.0, s)), 1))
        e = -0.3 + 1.2 * g1[i] - 0.03 * (age[i] - 55) + 0.3 * cb - 0.2 * cc + 3.0 * (anc[i][0] - 0.35)
        recovered.append("yes" if rng.random() < 1 / (1 + math.exp(-e)) else "no")
    # the rare variant: three people who recovered
    sim = [i for i in range(118) if recovered[i] == "yes" and i not in (16, 43, 89)]
    for k in (3, 11, 29):
        geno[RARE][sim[k]] = 1
    # missing genotypes (~1%), never on the rare variant
    for (_, _, _, rs, _, _, _) in VARIANTS:
        if rs == RARE:
            continue
        for i in range(n):
            if rng.random() < 0.01:
                geno[rs][i] = -1
    return people, anc, age, sex, center, geno, score, recovered


def _vcf_text(people, geno):
    rng = random.Random(7)
    linhas = [
        "##fileformat=VCFv4.2",
        "##source=MMGEN teaching demo (simulated genotypes at real GRCh38 dbSNP 151 positions)",
        "##reference=GRCh38",
    ]
    for c, l in CONTIGS:
        linhas.append(f"##contig=<ID={c},length={l}>")
    linhas.append('##INFO=<ID=GENE,Number=1,Type=String,Description="Candidate gene">')
    linhas.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    linhas.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(people))
    for (gene, chrom, pos, rs, ref, alt, _) in sorted(VARIANTS, key=lambda v: (ORDER[v[1]], v[2])):
        gts = []
        for g in geno[rs]:
            if g < 0:
                gts.append("./.")
            elif g == 1:
                gts.append("0|1" if rng.random() < 0.5 else "1|0")
            else:
                gts.append("1|1" if g == 2 else "0|0")
        linhas.append(f"{chrom}\t{pos}\t{rs}\t{ref}\t{alt}\t.\tPASS\tGENE={gene}\tGT\t" + "\t".join(gts))
    return ("\n".join(linhas) + "\n").encode("ascii")


def _phenotype_text(people, anc, age, sex, center, score, recovered):
    linhas = ["sample\trecovered\trecovery_score\tage\tsex\tcenter\tanc_afr\tanc_eur\tanc_eas"]
    for i, p in enumerate(people[:118]):  # PT119 and PT120: genotyped, no follow-up
        a = "NA" if p == "PT017" else str(age[i])
        s = "NA" if p in ("PT044", "PT090") else f"{score[i]:.1f}"
        linhas.append(f"{p}\t{recovered[i]}\t{s}\t{a}\t{sex[i]}\t{center[i]}\t{anc[i][0]:.3f}\t{anc[i][1]:.3f}\t{anc[i][2]:.3f}")
    # two enrolled people whose DNA sample failed: no genotypes
    linhas.append("PT121\tyes\t71.0\t49\tF\tA\t0.050\t0.900\t0.050")
    linhas.append("PT122\tno\t48.5\t63\tM\tC\t0.100\t0.150\t0.750")
    return ("\n".join(linhas) + "\n").encode("ascii")


README = """# Candidate-gene association: a simulated teaching cohort

**Everything here is simulated.** The 40 SNVs are real common dbSNP 151 SNVs (GRCh38
positions and alleles) in five candidate genes, but the genotypes, the outcomes and
the "associations" were drawn by a program. They say nothing about what these
variants do in real people.

A fictional study followed 122 people after a treatment: did they recover (`recovered`,
yes/no), and how well (`recovery_score`, 0-100)? 120 were genotyped (`cohort.vcf.gz`).
`phenotypes.tsv` also has age, sex, the hospital (`center`: A, B, C) and each person's
ancestry proportions (`anc_afr`, `anc_eur`, `anc_eas`, which sum to 1).

## Try it

```python
import multimerge_genomics as mmg
cov = ["age", "sex", "center", "anc_afr", "anc_eas"]   # two of the three proportions
r = mmg.genotype_association("cohort.vcf.gz", "phenotypes.tsv", outcome="recovery_score", covariates=cov)
print(r["text"])
```

## What it shows

{numbers}

## Lessons

1. **Who was analysed.** Two genotyped people have no phenotype row, two rows have no
   genotypes, one person has no age and two no recovery score: the report says who was
   left out and why, and each test uses the people with a genotype at that variant.
2. **A planted effect** (OPRM1 rs12190259) is found, but look at the confidence interval
   and at the power line: with about 120 people only large effects can be detected, and
   an effect that is just significant is probably overestimated.
3. **Confounding.** BDNF rs2030324 is more frequent in people with more AFR-like
   ancestry, and ancestry changes the outcome here. Without the ancestry covariates the
   variant looks associated; with them it is not. Population structure is the classic
   source of false associations.
4. **Proportions that sum to 1.** Give `anc_afr`, `anc_eur` and `anc_eas` together and
   the analysis refuses: with the intercept, one of them is a combination of the others.
   Leave one out.
5. **Separation.** IL6 rs2069830 is carried by three people who all recovered: the
   plain logistic fit has no finite estimate, so Firth's penalized likelihood is used
   (method `firth`). Three carriers cannot establish anything; the wide interval says so.
6. **Many tests.** 40 variants are tested: look at the FDR column, and at
   `permutations=1000` for family-wise p-values.

Sources of the positions: dbSNP build 151 "common" (GRCh38), RefSeq GRCh38.p14 gene spans.
"""


def write_demo(folder, numbers=None):
    """Writes the demo into `folder` (created if needed). Returns the file paths."""
    os.makedirs(folder, exist_ok=True)
    people, anc, age, sex, center, geno, score, recovered = simulate()
    arquivos = {
        "cohort.vcf.gz": _bgzf(_vcf_text(people, geno)),
        "phenotypes.tsv": _phenotype_text(people, anc, age, sex, center, score, recovered),
        "candidate_genes.bed": ("".join(f"{c}\t{a}\t{b}\t{g}\n" for c, a, b, g in GENES)).encode("ascii"),
        "README.md": README.replace("{numbers}", numbers if numbers is not None else NUMBERS).encode("utf-8"),
        "VERSION": f"{DEMO_VERSION}\n".encode("ascii"),
    }
    caminhos = []
    for nome, dados in arquivos.items():
        p = os.path.join(folder, nome)
        with open(p, "wb") as f:
            f.write(dados)
        caminhos.append(p)
    return caminhos


NUMBERS = """Numbers from MMGEN's genotype_association on these files (additive model, defaults
min_maf=0.01 and min_call_rate=0.9; 40 alleles tested in every run):

* **People.** 120 genotyped, 120 phenotype rows. For `recovery_score`, 115 are analysed:
  PT121 and PT122 have no genotypes, PT119 and PT120 no phenotype row, PT044 and PT090 no
  score, PT017 no age. For `recovered`, 117 (74 yes, 43 no).
* **The planted effect**, `recovery_score` adjusted for age, sex, center, anc_afr, anc_eas:
  OPRM1 rs12190259 (A>C), beta +7.64 points per C allele, 95% CI [4.42, 10.87],
  p = 7.9e-6, FDR = 0.0003: the only allele with FDR < 0.05 (the simulated effect is 6).
  With `permutations=1000`: permutation p = 0.001 and family-wise p = 0.001 (the next
  allele, p = 0.017, has a family-wise p of 0.52).
* **Confounding.** Adjusted only for age, sex and center, BDNF rs2030324 (A>G) is the
  second hit: beta +5.96 [2.04, 9.87], p = 0.0032. Adding anc_afr and anc_eas:
  beta -0.30 [-3.87, 3.28], p = 0.87. Its effect was ancestry all along.
* **A yes/no outcome loses information.** For `recovered` (same covariates) the same
  planted variant gives OR 2.01 [0.92, 4.40], p = 0.074: not detected, although its
  simulated odds ratio is 3.3 per copy. The power line says why: with these 117 people
  an allele at frequency 0.16 needs an odds ratio of at least 4.37 for 80% power at
  0.05/40.
* **Separation.** IL6 rs2069830 (3 carriers, all recovered): method `firth`,
  OR 4.21 [0.27, 633], p = 0.33.
* **All three ancestry proportions** as covariates: refused, "The covariate anc_eas is a
  linear combination of the intercept and the covariates listed before it".
"""

if __name__ == "__main__":
    for p in write_demo(sys.argv[1] if len(sys.argv) > 1 else "association_demo"):
        print(p, os.path.getsize(p))
