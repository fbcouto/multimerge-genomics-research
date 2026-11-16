# `mmgen-pgx` Validation Report (VCF genes)

**Engine:** `mmgen-pgx` (`engine_version: "1.1.0"`, `schema_version: "1"`), MMGEN 2.1.0
**Table bundle:** `cpic-v1.61.0` (`data_version: "CPIC v1.61.0; Open Traits v1.0.0"`, assembly `GRCh38`, CPIC tables `CC0-1.0`)
**Reference schema:** `pgx_result.schema.json` (JSON Schema Draft 2020-12)
**Validation suite:** [trials/python/validate_pgx.py](../trials/python/validate_pgx.py), [trials/python/smoke.py](../trials/python/smoke.py), [trials/src/bin/farmacogenetica.rs](../trials/src/bin/farmacogenetica.rs)

This report replaces the first one (engine 1.0.0). Engine 1.1.0 **chooses the diplotype as PharmCAT does** (section 2): by the Named Allele Matcher score, and for NAT2 by population frequency. Engine 1.0.0 left a gene uncalled when the positions read allowed answers with different effects. Every number below comes from runs of the 2.1.0 wheel (2026-11-16, Windows, CPU only).

---

## 1. Summary

| Criterion | Result |
| :--- | :--- |
| JSON Schema | **0 errors** in the 9 GIAB results, the run without a callable BED and PharmCAT's result through `pgx_run` |
| Determinism | Repeated runs give byte-identical `result.json` |
| CPU time | 0.27 – 0.54 s per GIAB VCF on an idle machine (4.0 s for the 9) |
| PharmCAT 3.4.0 concordance (171 gene × dataset pairs) | **149 `EXACT`, 20 `BOTH_UNCALLED`, 0 `SAFETY_UNCALLED`, 2 `FORMAT_OR_CALL_DIFF`** (engine 1.0.0: 139, 20, 8, 4) |
| Every call with an alternative that changes the phenotype or a recommendation carries a warning | **Yes**: 5 such calls, 5 warnings |
| No silent calls | **Pass** (section 6): nothing read → not called; a key position not read → called with a warning that names it |
| Trio Mendelian consistency | **0 conflicts**; the chosen calls fit in 33 of 34 checks; the other (HG003 CYP2B6, where PharmCAT makes the same call) fits through an alternative the engine listed |
| GeT-RM NA12878 (HG001) | Same as PharmCAT in every gene (CYP2C19 `*1/*2`, NAT2 `*4/*5`); in SLCO1B1 our call is the consensus `*1/*15` where PharmCAT gives two options |
| Rust suite | **110 of 110** with the HG002 data, 94 of 94 without it |
| Python smoke suite | **70 of 70** (CPU wheel), **71 of 71** (CUDA wheel) |

So both acceptance criteria pass:
- "0 `SAFETY_UNCALLED`, CYP2C19 of HG001 = `*1/*2`";
- "GeT-RM ≥ PharmCAT".

What remains different from PharmCAT is in section 5.2.

---

## 2. How the diplotype is chosen

For each gene, every diplotype whose two alleles fit the genotypes read is a candidate. When the VCF has phase sets, the candidates must respect them.

**The score (all genes except NAT2).** Each candidate is scored as PharmCAT's Named Allele Matcher does:

- **an allele is worth the number of positions read that define it**; a position that was not read drops out of the count;
- **the reference allele is defined at every position**, so it is worth every position read (PharmCAT: "the reference allele will always have the maximum score");
- **`*1` is scored as a reference too.** In CPIC's tables `*1` is the reference everywhere except in CYP2C19. There the GRCh38 sequence is `*38`, and `*1` is listed by its one difference (rs3758581). ClinPGx, the source of PharmCAT's definitions, writes CYP2C19 `*1` out at every position. Scored by its one position, `*2/*35` (4 + 3) would beat `*1/*2` (1 + 4) in NA12878; scored as PharmCAT does, `*1/*2` wins (38 against 7 there);
- **a diplotype is worth the sum of its two alleles**, and the highest-scoring one is called (`x_chosen_by: "score"`, `x_score`).

**Population frequency (NAT2).** PharmCAT does not pick NAT2 by score. Its `exemptions.json` carries 29 "unphased priorities" chosen from frequency data: NAT2's reference allele `*1` is the GRCh38 haplotype, rare in every population, and a score that counts every position read for it prefers it over the common `*4`, `*5` and `*6`.

`mmgen-pgx` does the same with its own computation from CPIC's allele frequencies (public domain), which the table bundle now carries per allele and per population group (`pop_freqs`). The rule is marked per gene in the tables (`"diplotype_choice": "frequency"`, NAT2 only):

- the expected frequency of a diplotype is computed in each CPIC group: Latino, European, East Asian, Central/South Asian, Sub-Saharan African. It is 2·f(a)·f(b) when the two alleles differ, and f(a)² when they are the same;
- the mean over the groups is taken, and the most common candidate is called; the score breaks ties;
- the call records it in `x_chosen_by: "frequency"` and `x_frequency`, and the note gives what the score alone would have picked.

In the GIAB sets this gives PharmCAT's NAT2 calls: HG001 `*4/*5` (0.136, against 3.8 × 10⁻⁴ for `*1/*16`) and HG003 `*5/*6` (0.157, against 1.4 × 10⁻⁵ for `*1/*30`).

**Why only NAT2.** With the frequency rule applied to every star-allele gene, the GIAB sets give 145 `EXACT` instead of 149:
- **gain:** it fixes CYP2B6 of HG003 (`*5/*6`, which the trio confirms);
- **loss:** it breaks UGT1A1 in three sets. CPIC gives no frequency for the combined allele `*80+*28`, so the rule picks `*28/*80` and the trio contradicts it.

So the frequency rule stays where PharmCAT uses it.

**What goes with the call.** Ties are all reported, joined by `" OR "`. When tied diplotypes disagree in phenotype or activity score:
- the phenotypes are joined the same way;
- no activity score is given, and no recommendation is matched for the gene;
- a warning says so.

The other candidates go to the gene's `notes` ("Also possible with the positions read: `*2/*35` (Poor Metabolizer)") and to `x_alternatives`. Each entry has the diplotype, phenotype, activity score, score, frequency for NAT2, `changes_phenotype` and `changes_recommendation`. The engine checks each alternative against the guidelines; when one would change the phenotype or a recommendation, `warnings` gets a line such as:

> CYP2C19: `*1/*2` (Intermediate Metabolizer) was chosen as the highest-scoring diplotype, but `*2/*35` (Poor Metabolizer) is also possible with the positions read and would change the phenotype or a recommendation. Phased genotypes would settle it.

**Genes called by variants** have no scoring: DPYD (by activity), and CFTR, RYR1, CACNA1S, MT-RNR1 (by status). When the variants are unphased and cis and trans give different results:
- for DPYD the engine follows CPIC and PharmCAT: the variants are taken to be on different copies, so the lowest-function ones all count (the lowest activity score, `x_chosen_by: "trans"`);
- the cis result goes to `x_alternatives`, with a warning.

**`called: false` is left for two cases:**
- nothing could be read: the sequence is absent, or every position is missing or outside the callable regions;
- no definition fits what was read: an allele outside the curated definitions, heteroplasmy in MT-RNR1, or more than one variant in a haploid gene.

A key position that was not read no longer stops the call. The call uses the positions read, and a note, a warning and `x_missing_key_positions` name the position. The alleles defined only there are not tested, and the warning says so.

---

## 3. Table bundle and licensing

The tables are built by [trials/python/build_pgx_bundle.py](../trials/python/build_pgx_bundle.py) from the official CPIC database dump and GRCh38 reference windows. The dump is `CC0-1.0`, at `https://files.cpicpgx.org/data/database/cpic_db_dump-v1.61.0.sql.gz`. No PharmCAT code or data is copied.

Rebuilt from the dump, `phenotypes.json`, `guidelines.json` and `traits.json` come out **byte-identical** to the earlier bundle. The other two files changed in this round:
- `genes.json`: per-population allele frequencies, and NAT2's `diplotype_choice`;
- `bam_genes.json`: the same frequencies for CYP2D6, and the IPD-IMGT/HLA pin (see the v2 report).

The engine checks the size and SHA-256 of every file against `manifest.json` before each run:

| File | Bytes | SHA-256 (prefix) | Contents |
| :--- | ---: | :--- | :--- |
| `genes.json` | 1,415,814 | `520822874700…` | 19 genes, 1,052 positions (80 key), 1,156 allele definitions; frequencies for 508 alleles |
| `phenotypes.json` | 44,588 | `b7cc096fcb57…` | Diplotype / activity score → phenotype, per gene |
| `guidelines.json` | 1,421,716 | `e43c85b7bce7…` | 2,046 CPIC recommendations, 110 drugs |
| `traits.json` | 9,234 | `4e164c7c2986…` | Lactase persistence, eye colour (teaching only) |
| `bam_genes.json` | 282,517 | `5241d0ab63e0…` | CYP2D6 and HLA-A/B/C (v2) |

The zip of the bundle, `mmgen_pgx_tables_cpic-v1.61.0.zip`, is reproducible (fixed timestamps): SHA-256 `e329b8ec6e8567dc2570cd3823a96401e54135b3a5545305fb05c7265caa1f85`.

### 3.1 Where CPIC's definitions differ from PharmCAT's

Compared allele by allele with the definitions inside `pharmcat-3.4.0-all.jar` (ClinPGx, 2026-07-13), the 19 genes differ in three ways that matter:

1. **CYP2C19 `*1`**: one position in CPIC, every position in PharmCAT. The scoring rule above handles it.
2. **CYP2C19 `*2`**: CPIC lists `R` (A or G) at rs12769205 and rs3758581, PharmCAT `G`. With CPIC's `R`, `*2/*35` fits NA12878's unphased genotypes; with PharmCAT's `G` it does not. So `*2/*35` is reported here as a lower-scoring alternative, and PharmCAT does not mention it.
3. **SLCO1B1 `*45`**: one allele in CPIC, two sub-alleles (`*45.001`, `*45.002`) in PharmCAT, one of which fits NA12878 (section 5.2).

The other differences are of two kinds:
- the same indels written two ways (CPIC `delT`, VCF-style `G`);
- alleles present in only one of the two sources.

---

## 4. Gene scope and calling models

| Gene | Chrom | Model | Lookup | Ploidy | Positions | Key | Alleles | Reference | Choice |
| :--- | :--- | :--- | :--- | :--- | ---: | ---: | ---: | :--- | :--- |
| `ABCG2` | `chr4` | `variant_genotype` | `PHENOTYPE` | diploid | 1 | 1 | 2 | `rs2231142 reference (G)` | score |
| `CACNA1S` | `chr1` | `variant_status` | `PHENOTYPE` | diploid | 2 | 2 | 3 | `Reference` | variants |
| `CFTR` | `chr7` | `variant_status` | `PHENOTYPE` | diploid | 93 | 0 | 104 | `ivacaftor non-responsive CFTR sequence` | variants |
| `CYP2B6` | `chr19` | `star_allele` | `PHENOTYPE` | diploid | 48 | 8 | 47 | `*1` | score |
| `CYP2C19` | `chr10` | `star_allele` | `PHENOTYPE` | diploid | 43 | 8 | 46 | `*38` | score |
| `CYP2C9` | `chr10` | `star_allele` | `ACTIVITY_SCORE` | diploid | 88 | 8 | 94 | `*1` | score |
| `CYP3A5` | `chr7` | `star_allele` | `PHENOTYPE` | diploid | 7 | 3 | 8 | `*1` | score |
| `CYP4F2` | `chr19` | `star_allele` | `PHENOTYPE` | diploid | 20 | 2 | 22 | `*1` | score |
| `DPYD` | `chr1` | `variant_activity` | `ACTIVITY_SCORE` | diploid | 83 | 5 | 84 | `Reference` | variants (trans) |
| `G6PD` | `chrX` | `star_allele` | `PHENOTYPE` | sex_x | 172 | 7 | 186 | `B (reference)` | score |
| `IFNL3` | `chr19` | `variant_genotype` | `ALLELE_STATUS` | diploid | 1 | 1 | 2 | `rs12979860 reference (C)` | score |
| `MT-RNR1` | `chrM` | `variant_status` | `PHENOTYPE` | haploid | 22 | 9 | 25 | `Reference` | variants |
| `NAT2` | `chr8` | `star_allele` | `PHENOTYPE` | diploid | 47 | 9 | 59 | `*1` | **frequency** |
| `NUDT15` | `chr13` | `star_allele` | `PHENOTYPE` | diploid | 20 | 4 | 22 | `*1` | score |
| `RYR1` | `chr19` | `variant_status` | `PHENOTYPE` | diploid | 320 | 0 | 347 | `Reference` | variants |
| `SLCO1B1` | `chr12` | `star_allele` | `PHENOTYPE` | diploid | 35 | 4 | 45 | `*1` | score |
| `TPMT` | `chr6` | `star_allele` | `PHENOTYPE` | diploid | 45 | 4 | 49 | `*1` | score |
| `UGT1A1` | `chr2` | `star_allele` | `PHENOTYPE` | diploid | 4 | 4 | 9 | `*1` | score |
| `VKORC1` | `chr16` | `variant_genotype` | `PHENOTYPE` | diploid | 1 | 1 | 2 | `rs9923231 reference (C)` | score |

---

## 5. GIAB benchmarks against PharmCAT 3.4.0 (9 datasets)

`validate_pgx.py` runs `pharmcat-3.4.0-all.jar` (through `pharmcat_run`, with the callable BED) and `mmgen-pgx` on the 9 GIAB GRCh38 benchmark VCF + BED pairs: 171 gene × dataset pairs. A diplotype is `EXACT` when both engines give the same set, ignoring the order inside a diplotype and among `" OR "` options.

| Dataset | Called | Drugs | `EXACT` | `BOTH_UNCALLED` | `SAFETY_UNCALLED` | `FORMAT_OR_CALL_DIFF` | Calls with alternatives |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `HG001_v421` (chr1–22) | 17 / 19 | 48 | 16 | 2 | 0 | 1 | 4 |
| `HG001_v332` (chr1–X, phased) | 14 / 19 | 45 | 13 | 5 | 0 | 1 | 0 |
| `HG002_v421` | 17 / 19 | 48 | 17 | 2 | 0 | 0 | 0 |
| `HG002_v50q` (with chrX) | 18 / 19 | 81 | 18 | 1 | 0 | 0 | 0 |
| `HG003_v421` | 17 / 19 | 48 | 17 | 2 | 0 | 0 | 3 |
| `HG004_v421` | 17 / 19 | 48 | 17 | 2 | 0 | 0 | 0 |
| `HG005_v421` | 17 / 19 | 48 | 17 | 2 | 0 | 0 | 0 |
| `HG006_v421` | 17 / 19 | 48 | 17 | 2 | 0 | 0 | 0 |
| `HG007_v421` | 17 / 19 | 48 | 17 | 2 | 0 | 0 | 1 |
| **Total** | **151 / 171** | — | **149** | **20** | **0** | **2** | **8** |

`BOTH_UNCALLED` (20):
- MT-RNR1 in all 9 (no chrM);
- G6PD in the 7 chr1–22 releases (no chrX);
- CYP2B6, IFNL3 and VKORC1 in `HG001_v332`, whose BED leaves their positions out.

### 5.1 The 8 `SAFETY_UNCALLED` of engine 1.0.0

| Pair | Engine 1.0.0 | Engine 1.1.0 | PharmCAT |
| :--- | :--- | :--- | :--- |
| `HG001_v421` CYP2C19 | not called (`*1/*2 OR *2/*35`, conflicting phenotypes) | `*1/*2` (Intermediate), alternative `*2/*35` (Poor), warning | `*1/*2` |
| `HG001_v421` UGT1A1 | not called (`*1/*80+*28 OR *28/*80`) | `*1/*80+*28` (Intermediate), alternative `*28/*80` (Indeterminate), warning | `*1/*80+*28` |
| `HG003_v421` UGT1A1 | not called | `*1/*80+*28`, same alternative and warning | `*1/*80+*28` |
| `HG003_v421` NAT2 | not called (`*1/*30 OR *5/*6 OR *16/*34`) | `*5/*6` (Poor) by frequency, alternatives `*16/*34` (Poor), `*1/*30` (Intermediate), warning | `*5/*6` |
| `HG007_v421` NUDT15 | not called (`*1/*3 OR *3/*6`) | `*1/*3` (Intermediate), alternative `*3/*6` (Poor), warning | `*1/*3` |
| `HG001_v332` CYP4F2 | not called (key rs2108622 not read) | `*1/*1`, warning naming rs2108622 | `*1/*1` |
| `HG001_v332` TPMT | not called (keys rs1142345, rs56161402 not read) | `*1/*1`, warning naming both | `*1/*1` |
| `HG001_v332` UGT1A1 | not called (key rs3064744 not read) | `*1/*80 OR *1/*80+*28 OR *1/*80+*37` (tie: Indeterminate OR Intermediate), no recommendation, warnings | the same three |

### 5.2 The 2 `FORMAT_OR_CALL_DIFF`

**`HG001_v421` and `HG001_v332` SLCO1B1**: `*1/*15` (Decreased Function) here, `*1/*15 OR *1/*45` in PharmCAT, same phenotype. PharmCAT splits `*45` into two sub-alleles, one of which fits NA12878; CPIC v1.61.0's single `*45` does not. `*1/*15` is the GeT-RM consensus.

The two NAT2 differences of the first version of this change (`*1/*16` and `*1/*30` where PharmCAT has `*4/*5` and `*5/*6`) are gone with the frequency rule.

### 5.3 The 5 calls whose alternatives would change the phenotype or a recommendation

All 5 carry a warning:

- `HG001_v421` CYP2C19: `*1/*2` (Intermediate); also possible `*2/*35` (Poor)
- `HG001_v421` UGT1A1: `*1/*80+*28` (Intermediate); also possible `*28/*80` (Indeterminate)
- `HG003_v421` NAT2: `*5/*6` (Poor); also possible `*1/*30` (Intermediate)
- `HG003_v421` UGT1A1: `*1/*80+*28` (Intermediate); also possible `*28/*80` (Indeterminate)
- `HG007_v421` NUDT15: `*1/*3` (Intermediate); also possible `*3/*6` (Poor)

---

## 6. Trio Mendelian consistency and "no silent calls"

Both GIAB trios, v4.2.1 calls:
- `PASS`: one allele of the child comes from each parent;
- `PASS_ALT`: the chosen calls do not fit, but they fit when the alternatives the engine listed are included;
- `SKIP`: a member is not called.

| Trio | `PASS` | `PASS_ALT` | `SKIP` | Conflicts |
| :--- | ---: | ---: | ---: | ---: |
| Ashkenazi (`HG002 ← HG003, HG004`) | 16 | 1 | 2 (G6PD, MT-RNR1) | **0** |
| Han Chinese (`HG005 ← HG006, HG007`) | 17 | 0 | 2 (G6PD, MT-RNR1) | **0** |

The `PASS_ALT` is **CYP2B6 of the father, HG003**:
- the engine chose `*1/*7`, as PharmCAT does, but the son `*2/*5` and the mother `*2/*5` need `*5/*6` in the father;
- `*7` carries the `*6` and `*5` variants on one copy, so `*1/*7` and `*5/*6` are the cis and trans readings of the same genotypes, and the reference allele makes `*1/*7` score higher;
- the phenotype is the same (Intermediate), so there is no warning, and `*5/*6` is in the alternatives.

The frequency rule would pick `*5/*6` here, but costs more elsewhere (section 2).

**No silent calls** (`validate_pgx.py`, section 2):

1. **No callable BED** (`HG002_v50q`, `callable_bed=None`): every absent position counts as missing.
   - 12 of 19 genes are not called (16 with engine 1.0.0); TPMT and MT-RNR1 are among them.
   - The genes still called read only a few positions (CYP2C19: 1 of 43). They come out as long ties with every phenotype joined, no recommendation, and warnings for the tie and for each key position not read.
2. **Absent chromosomes** (`HG002_v421`): G6PD ("Sequence chrX is absent from the VCF and callable BED.") and MT-RNR1 are not called.
3. **Key positions removed from VCF and BED** (`HG002_v50q`, CYP2C19 rs4244285 and VKORC1 rs9923231):
   - VKORC1, whose only position is gone, is not called.
   - CYP2C19 is called (`*1/*1`) from its other positions, as PharmCAT would. The result has `x_missing_key_positions: ["rs4244285"]` and the warning "CYP2C19: key position(s) not read (rs4244285 (chr10:94781859)); the alleles defined there were not tested, so the call may miss them."
4. GRCh37 headers and tampered tables (SHA-256) are rejected (Rust suite).

---

## 7. GeT-RM NA12878 / HG001

Consensus genotypes of the CDC/Coriell GeT-RM panels (Pratt et al. 2016; Gaedigk et al. 2022):

| Gene | GeT-RM | `HG001_v421` (unphased) | `HG001_v332` (phased) | PharmCAT, `HG001_v421` |
| :--- | :--- | :--- | :--- | :--- |
| `CYP2B6` | `*1/*1` | **`*1/*1`** | not called (outside the BED) | `*1/*1` |
| `CYP2C9` | `*1/*2` | **`*1/*2`** | **`*1/*2`** | `*1/*2` |
| `CYP2C19` | `*1/*2` | **`*1/*2`** (alternative `*2/*35`, warned) | **`*1/*2`** | `*1/*2` |
| `CYP3A5` | `*3/*3` | **`*3/*3`** | **`*3/*3`** | `*3/*3` |
| `CYP4F2` | `*1/*1` | **`*1/*1`** | `*1/*1` (rs2108622 not read, warned) | `*1/*1` |
| `DPYD` | `*1/(*4)` | `c.1601G>A (*4) AND c.1627A>G (*5)` | `c.1601G>A (*4)/c.1627A>G (*5)` | same as ours |
| `NAT2` | `*4/*5` | **`*4/*5`** (by frequency; alternative `*1/*16`) | **`*4/*5`** | `*4/*5` |
| `SLCO1B1` | `*1/*15` | **`*1/*15`** | **`*1/*15`** | `*1/*15 OR *1/*45` |
| `TPMT` | `*1/*1` | **`*1/*1`** | `*1/*1` (2 keys not read, warned) | `*1/*1` |
| `UGT1A1` | `*60/*28` | `*1/*80+*28` (`*60` is not in CPIC's table) | `*1/*80 OR *1/*80+*28 OR *1/*80+*37` | `*1/*80+*28` |
| `VKORC1` | rs9923231 G/A | **`rs9923231 reference (C)/rs9923231 variant (T)`** (+ strand) | not called (outside the BED) | same |

Pratt 2016 tested DPYD `*4` only; HG001 also carries `*5` (both normal function). Against GeT-RM, `mmgen-pgx` agrees wherever PharmCAT agrees, and in SLCO1B1 it gives the consensus alone.

---

## 8. Tests

- **Rust** ([trials/src/bin/farmacogenetica.rs](../trials/src/bin/farmacogenetica.rs)): **110 of 110** with `PGX_V2_BAM`, `MMGEN_HLA_DIR` and `MMG_HUMAN_GENOMES` set, 94 of 94 without them.
  - The earlier checks still pass: SHA-256 tampering, `.zip` rejection, GRCh37, multi-sample VCFs, left-alignment with a FASTA, haploid chrX, chrM indels, `.vcf` against `.vcf.gz` + `.tbi`.
  - The choice of the diplotype is tested directly:
    - an unphased star-allele gene where the higher score wins and the lower goes to `x_alternatives`, changing the phenotype and the recommendation, with the note and the warning;
    - a CYP2C19-shaped gene (`*38` reference, `*1` at one position, `*2` with IUPAC `R`): unphased gives `*1/*2` over `*2/*35` (6 against 5); phased as `*2/*35`, it gives `*2/*35` alone;
    - two NAT2-shaped genes with the same alleles: with `diplotype_choice: "frequency"` the common `*4/*5` is called over the higher-scoring `*1/*16`, and without it the score picks `*1/*16`;
    - a tie with different activity scores: both are joined, with no activity score and no recommendation;
    - a key position not read: the gene is called, with `x_missing_key_positions` and a warning;
    - the DPYD-style trans assumption: the lowest activity score is called, and cis goes to the alternatives;
    - the warnings and the alternatives appear in `report.html`.
- **Python smoke suite** ([trials/python/smoke.py](../trials/python/smoke.py)): **70 of 70** with the CPU wheel, **71 of 71** with the CUDA wheel.
- **`validate_pgx.py`**: 0 hard failures. It fails if:
  - any call with a phenotype- or recommendation-changing alternative lacks a warning;
  - CYP2C19 of HG001 is not `*1/*2`;
  - a trio conflicts even with the alternatives;
  - a key position is dropped without a warning.

Experimental software for teaching and research: not a medical diagnosis or treatment guidance.
