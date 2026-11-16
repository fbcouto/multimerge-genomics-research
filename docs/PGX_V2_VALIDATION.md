# `mmgen-pgx` v2 Pharmacogenetics Validation (`CPIC v1.61.0` + BAM/CNV/HLA)

## 1. Executive Summary

`mmgen-pgx` v2 extends the verified pure-Rust GRCh38 pharmacogenetics engine (`src/farmacogenetica.rs`, `docs/PGX_V1_VALIDATION.md`) with BAM+VCF structural and star-allele calling in `src/pgx_bam.rs`, adding **`CYP2D6`** (copy-number, deletion/duplication/hybrid SV detection, BAM pileup filling, and read-pair phasing) and Class I **`HLA`** (**`HLA-A`**, **`HLA-B`**, **`HLA-C`**) alongside the 19 VCF-based CPIC genes (`23` total pharmacogenes when `bam` is provided, plus `2` non-clinical genetic trait panels), together with PharmCAT `-po` outside-call export (`outside_calls.tsv`).

Updated 2026-11-16 for engine 1.1.0 (MMGEN 2.1.0), after two changes. The VCF genes are now chosen as PharmCAT chooses them: by the Named Allele Matcher score, and NAT2 by population frequency (see `docs/PGX_V1_VALIDATION.md`, section 2). And **the tables pin the IPD-IMGT/HLA release** used for HLA typing (3.66.0, each FASTA checked by SHA-256; section 2.3). With the 2.1.0 wheels built from that code, the three suites passed with zero failures:

- **Rust integration & unit binary** (`trials/src/bin/farmacogenetica.rs`): `110 passed, 0 failed` with the real data (`PGX_V2_BAM`, `MMGEN_HLA_DIR`, `MMG_HUMAN_GENOMES`), `94 passed` without it: all v1 assertions, the scoring, synthetic BAM CNV/hybrid/HLA/error guards, the IPD-IMGT/HLA pin, and real `HG002` `v4.2.1` + `v5.0q` BAM+VCF integration checks.
- **Python smoke suite** (`trials/python/smoke.py`): `70 of 70 checks passed` with the CPU wheel, `71 of 71` with the CUDA wheel (including synthetic BAM+VCF `mmgen_pgx_run` / `pgx_run` end-to-end execution).
- **GIAB + v2 BAM validation suite** (`trials/python/validate_pgx.py` against `pgx_result.schema.json`, `pharmcat-3.4.0-all.jar`, and the HG002 BAM cut): `0 hard failures` across all 9 GIAB VCF benchmarks, both Mendelian trios, GeT-RM `NA12878`, trait panels, `HG002` v2 `CYP2D6` + `HLA-A/B/C` calling, `PharmCAT -po` round-trip (`callSource="OUTSIDE"`), the offline missing-FASTA guard, and the live download of the pinned IPD-IMGT/HLA release from GitHub (`MMGEN_TEST_IMGT_DOWNLOAD=1`).

CYP2D6 keeps its own resolver (copy number, BAM phasing, IUPAC penalty, frequency tie-break); the Named Allele Matcher scoring applies to it only when the BAM resolver finds no pair and the call falls back to the VCF genes' caller.

---

## 2. Architecture, Bundle Assets & License Separation

### 2.1 CPIC v1.61.0 Bundle (`mmgen_tools/pgx/cpic-v1.61.0/`)

Built via `trials/python/build_pgx_bundle.py` and verified against `mmgen_tools/pgx/cpic-v1.61.0/manifest.json`:

| Bundle File | Bytes | SHA-256 | Role in v2 |
| --- | ---: | --- | --- |
| `genes.json` | 1,415,814 | `520822874700a81623c2dd2c07c643416a36ac5cc96278add65ebd6063a998c5` | 19 VCF-only CPIC gene definitions, with per-population allele frequencies |
| `bam_genes.json` | 282,517 | `5241d0ab63e0e5fd486227e49b9f63b27a8d354dc4f32b5a15d5297e92e6b592` | `CYP2D6` (157 positions, 172 alleles, 6 SVs, 27 duplications, 5 canonical IUPAC defaults, CNV/control regions) + `HLA-A/B/C` metadata with the pinned IPD-IMGT/HLA release |
| `guidelines.json` | 1,421,716 | `e43c85b7bce7810fdb3c2206ab6d93943bf08dd9e25f77ad7133f185717f6d1a` | CPIC gene-drug recommendations including `ACTIVITY_SCORE` (`CYP2D6`) and `ALLELE_STATUS` (`HLA-A`, `HLA-B`) |
| `phenotypes.json` | 44,588 | `b7cc096fcb577b5246911c9ef7a0c9ad0e62a1e6ec078ddcfb6905140639533d` | Diplotype/activity-score/phenotype lookup tables (including `activity_label` such as `>=3.0` / Unicode `≥3.0`) |
| `traits.json` | 9,234 | `4e164c7c298680b3295cac6b5a8008e3df4e34d7156b0deb88c0e6ac7c30cc0c` | Non-clinical trait panels (`lactase_persistence`, `eye_color`) |

### 2.2 BAM Omission vs. BAM Inclusion Semantics

- **`bam=None` (VCF-only mode):** `genes[]` contains only the 19 VCF genes, `outside_calls.tsv` is not written (`x_files.pharmcat_outside_calls` is omitted/null), and an explicit warning (`warnings` / `x_bam_note`) states that `CYP2D6` and Class I `HLA` calling were skipped.
- **`bam=<path>` (v2 BAM+VCF mode):** `genes[]` expands to 23 genes (`19` VCF + `CYP2D6` with `source="vcf+bam"` + `HLA-A`, `HLA-B`, `HLA-C` with `source="bam"`), each annotated with `copy_number` and `source`, rendered in HTML and written to `result.json`.

### 2.3 IPD-IMGT/HLA License Boundary & Runtime Cache (`src/pgx_bam.rs`)

No third-party IPD-IMGT/HLA sequences are bundled in source control, wheels, or bundle archives (their licence is CC BY-ND 4.0); `bam_genes.json` carries only their release and checksums:

| File | IPD-IMGT/HLA | Bytes | SHA-256 | URL (official GitHub mirror, one branch per release) |
| --- | --- | ---: | --- | --- |
| `A_gen.fasta` | 3.66.0 | 20,391,880 | `f265427f646d20da7f26f63f30620f1787c15417dfd731a39e8818026455d931` | `https://raw.githubusercontent.com/ANHIG/IMGTHLA/3660/fasta/A_gen.fasta` |
| `B_gen.fasta` | 3.66.0 | 24,329,147 | `3748309e4aa560f8a0c31a8dfd6f24ddce2f3c9a52bd66db193de5b2f1bcd8f8` | `https://raw.githubusercontent.com/ANHIG/IMGTHLA/3660/fasta/B_gen.fasta` |
| `C_gen.fasta` | 3.66.0 | 23,796,524 | `edec46685befd8d0f47871a6e4f916f4195a8b3d5930ec55ae3fa207e117893a` | `https://raw.githubusercontent.com/ANHIG/IMGTHLA/3660/fasta/C_gen.fasta` |

3.66.0 is the release the HG002 calls below were validated with (the cache used for them is byte-identical to branch `3660`, checked by git blob hash). `build_pgx_bundle.py --imgt-hla-version 3.66.0` writes these fields (with `--hla-dir`, it hashes local copies instead of downloading them).

- At runtime, each FASTA is looked for in `<MMGEN_HLA_DIR>/3.66.0/`, then in `<MMGEN_HLA_DIR>/` (or the default cache directory), and **used only if its size and SHA-256 are the pinned ones**. A file of another release is never used; the error names it and says why.
- If missing and `MMGEN_HLA_OFFLINE` is not `"1"`, `src/pgx_bam.rs` downloads the unmodified file with `curl` from the URL above into `<MMGEN_HLA_DIR>/3.66.0/`, checks the SHA-256, and saves it only if it matches.
- If missing while `MMGEN_HLA_OFFLINE=1`, calling halts cleanly with an explicit error naming `MMGEN_HLA_OFFLINE`, the file, the release and its SHA-256.
- The result says which release was used: `data_version` becomes `CPIC v1.61.0; Open Traits v1.0.0; IPD-IMGT/HLA 3.66.0` when a BAM is given, and each HLA gene's note names the release and the checksum prefix.
- Screen 19 downloads the same release from the same mirror and checks the same SHA-256 sums; `test_ui.py` checks that the screen's constants equal the pins in `bam_genes.json`.

---

## 3. Core Algorithms & Compatibility Fixes (`src/pgx_bam.rs`, `src/farmacogenetica.rs`)

1. **BAM Header & Index Validation (`validate_bam_header`):**
   - Verifies `.bai` indexed random access via `indice::consultar`, checks GRCh38 `@SQ` contig lengths against expected assembly lengths (rejecting GRCh37/non-GRCh38 BAMs), and rejects `SO:queryname` headers.
   - Permits `SO:unsorted` alongside `SO:coordinate` because coordinate-sorted indexed GIAB Novoalign BAMs (specifically `HG002.GRCh38.2x250.bam` and its extracted PGx subset `HG002.GRCh38.2x250.pgx_v2.bam`) carry `@HD VN:1.0 SO:unsorted`.
2. **`CYP2D6` CNV, Hybrid, Pileup & Read-Pair Phasing (`source: "vcf+bam"`):**
   - Measures normalized depth ratios over `target_whole` (`chr22:42126499-42130810`), `target_5p` (`chr22:42128000-42130810`), and `target_ex9` (`chr22:42126500-42126800`) against 7 diploid control regions (`ALB`, `EGFR_3p`, `RPP30`, `GAPDH`, `VDR`, `RPPH1`, `CYP2D8P`), accounting for paralog homology (`chr22:42139576-42144483`).
   - Resolves whole-gene deletions (`CN=0` -> `*5/*5`, `CN=1` -> `*X/*5`), VAF-skewed copy-number gains (`CN >= 3` such as `*1x2/*4`, using `format_activity_score` for `0.25`/`0.75` vs `0.5` formatting and `activity_label` lookup including `≥3.0`), and exon-9 conversion hybrids where `target_ex9` drops relative to `target_5p` (`*36`, `*68`, `*10/*36`, etc.).
   - Merges BAM pileup base observations with VCF genotypes, resolves cis/trans phase across heterozygous sites using paired-end read fragments (`x_phased=true`), and applies `canonical_iupac_defaults` tiebreaking when remaining candidates are ambiguous.
3. **Class I `HLA` Pure-Rust Read Caller (`HLA-A`, `HLA-B`, `HLA-C`, `source: "bam"`):**
   - Parses IMGT/HLA FASTAs (`load_unique_full_imgt`), retaining full-length genomic sequences (`>= 2400` bp for production IMGT FASTAs, `>= 100` bp for minimal synthetic FASTAs), deduplicating identical sequences in numeric allele sort order.
   - Maps paired/unpaired reads extracted from the BAM `HLA-A` (`chr6:29936000-29950000`), `HLA-B` (`chr6:31350000-31362000`), and `HLA-C` (`chr6:31265000-31276000`) windows across forward and reverse-complement orientations using 3 query seed offsets `[min(20, l/4), l/2 - 9, l - 40]`, enforcing `min_ov = min(70, 3*l/4)`, mismatch fraction `<= 5%`, mismatch penalty `-8.0` (`best_ov - 8.0 * best_mm`), unpaired mate penalty `-20.0`, and joint diploid scoring with `-30.0` heterozygosity parsimony.
   - Emits 2-field clinical resolution (`*01:01/*26:01`) and evaluates CPIC `ALLELE_STATUS` risk tokens (`*31:01` for `HLA-A`; `*15:02`, `*57:01`, `*58:01` for `HLA-B`; `HLA-C` has `lookup_method="PHENOTYPE"` and `phenotype=null`).
4. **CPIC Guideline Matching (`match_guidelines`):**
   - Skips `rec.untested_genes` fallback when all required genes for a guideline were actually called in v2, resolves `ACTIVITY_SCORE` via `activity_label` / `format_activity_score`, and matches `;`-delimited `ALLELE_STATUS` phenotypes (e.g., `*31:01 negative`, `*57:01 positive`).
5. **PharmCAT `-po` Outside-Call Export (`outside_calls.tsv`):**
   - When `bam` is supplied, writes headerless 2-column `GENE\tDIPLOTYPE` lines to `outside_calls.tsv` in `output_dir` for called `CYP2D6`, `HLA-A`, and `HLA-B` (omitting `HLA-C`, which PharmCAT `-po` does not accept), exposed via `x_files.pharmcat_outside_calls`.

---

## 4. Synthetic Unit & Integration Verification (`trials/src/bin/farmacogenetica.rs`)

Running `cargo run --release --manifest-path trials\Cargo.toml --bin farmacogenetica` compiled cleanly and passed **110 of 110 checks** with the real data (`94` without it). Key v2 BAM/CNV/HLA assertions verified:

| Test Group | Scenario & Input | Verified Outcome |
| --- | --- | --- |
| BAM omitted (`14`) | Minimal fixture VCF, `bam=None` | Preserves VCF-only gene count (`5` in fixture), omits `x_files.pharmcat_outside_calls`, appends BAM-skip warning |
| `CYP2D6 *2/*4` + HLA alerts (`15`) | Unphased VCF + synthetic BAM (`20x` `chr4 ALB` control, `10x`/copy `CYP2D6`, paired synthetic `HLA-A/B/C`) | Calls `CYP2D6 *2/*4` (`CN=2`, `source="vcf+bam"`, `x_phased=true`), `HLA-A *01:01/*31:01`, `HLA-B *07:02/*57:01`, `HLA-C *04:01/*12:03` (`source="bam"`); triggers `abacavir` and `carbamazepine`; writes `outside_calls.tsv` with `CYP2D6`, `HLA-A`, `HLA-B` and without `HLA-C` |
| `CYP2D6` deletions (`16`) | `0` and `1` synthetic `CYP2D6` copy vs `20x` `ALB` | Calls `*5/*5` (`CN=0`, triggers `codeine` guideline) and `*1/*5` (`CN=1`) |
| `CYP2D6` duplication (`17`) | Heterozygous VCF + `2:1` BAM haplotype copy ratio | Resolves `*1x2/*4` (`CN=3`) via BAM VAF skew |
| Exon-9 hybrid (`18`) | `*10/*36` VCF + BAM `target_ex9` copy drop | Resolves exon-9 hybrid `*10/*36` |
| Input / cache error guards (`19`) | GRCh37 BAM, `SO:queryname` BAM, missing `.bai`, and missing IMGT FASTA with `MMGEN_HLA_OFFLINE=1` | Rejects each cleanly with explicit English diagnostics (`MMGEN_HLA_OFFLINE`, `A_gen.fasta`, etc.) without panics |
| IPD-IMGT/HLA pin (`19b`) | Fixture tables pinned to a release with the synthetic FASTAs' SHA-256, then with a wrong SHA-256 for `B_gen.fasta` | The matching files are used (same HLA calls; `data_version` and the HLA notes name the release); the mismatching one is refused, with `HLA-B`, `B_gen.fasta`, the release and "SHA-256 differs" in the error |
| Real `HG002` v2 BAM (`20`) | `HG002.GRCh38.2x250.pgx_v2.bam` (`PGX_V2_BAM`) + GIAB `HG002` `v4.2.1` & `v5.0q` VCFs (`MMG_HUMAN_GENOMES`) + `cpic-v1.61.0` + IPD-IMGT/HLA 3.66.0 (`MMGEN_HLA_DIR`) | Returns all `23` genes; `CYP2D6=*2/*4 (CN=2)`, `HLA-A=*01:01/*26:01`, `HLA-B=*35:08/*38:01`, `HLA-C=*04:01/*12:03`, exact 3-line `outside_calls.tsv`, `data_version` naming IPD-IMGT/HLA 3.66.0 |

---

## 5. Python Smoke Suite & Full GIAB + HG002 v2 Validation

### 5.1 Package Smoke Suite (`trials/python/smoke.py`)

`python trials\python\smoke.py` (CPU-only wheel) exited with code `0` and reported **`70 of 70 checks passed`**, including end-to-end BAM coordinate sorting (`mmg.sort_genomic`), synthetic `mmgen_pgx_run` BAM+VCF execution (`CYP2D6 *2/*4 CN=2`, `HLA-A/B/C`, drug guidelines, and `outside_calls.tsv`), and wrapper parity via `mmg.pgx_run(..., engine="mmgen-pgx")`.

### 5.2 Nine GIAB v1 VCF Benchmarks (`trials/python/validate_pgx.py`, Section 1)

All 9 GIAB datasets validated against `pgx_result.schema.json` with `0` schema errors. The VCF-only results, their comparison with PharmCAT 3.4.0 (149 `EXACT`, 20 `BOTH_UNCALLED`, 0 `SAFETY_UNCALLED`, 2 `FORMAT_OR_CALL_DIFF`), the trios and GeT-RM are in `docs/PGX_V1_VALIDATION.md`. The runtimes below are from the run before the NAT2 frequency rule, made while a build was compiling alongside; on an idle machine a VCF takes 0.27 – 0.54 s:

| Dataset | Runtime (s, final run, a build compiling alongside) | Genes Called / 19 | Drug Recommendations | Schema |
| --- | ---: | :---: | ---: | --- |
| `HG001_v421` | 1.292 | 17 / 19 | 48 | `OK` |
| `HG001_v332` | 1.202 | 14 / 19 | 45 | `OK` |
| `HG002_v421` | 0.896 | 17 / 19 | 48 | `OK` |
| `HG002_v50q` | 0.534 | 18 / 19 | 81 | `OK` |
| `HG003_v421` | 1.167 | 17 / 19 | 48 | `OK` |
| `HG004_v421` | 1.259 | 17 / 19 | 48 | `OK` |
| `HG005_v421` | 0.980 | 17 / 19 | 48 | `OK` |
| `HG006_v421` | 0.803 | 17 / 19 | 48 | `OK` |
| `HG007_v421` | 0.862 | 17 / 19 | 48 | `OK` |

Sections 2–5 likewise passed with `0` hard failures: negative no-silent-call checks (`callable_bed=None`, missing `chrX`/`chrM`, dropped key SNPs, now called with a warning when other positions were read), Ashkenazi (`HG002/3/4`) and Chinese (`HG005/6/7`) Mendelian trio inheritance checks (`0` conflicts; CYP2B6 of HG003 fits only through a listed alternative, where PharmCAT makes the same call), GeT-RM `NA12878` concordance (the same as PharmCAT in every gene), and `lactase_persistence` / `eye_color` trait panels (including negative `HG002_nobed`).

### 5.3 Real `HG002` v2 BAM+VCF Validation & Consensus Comparison (`validate_pgx.py`, Section 6)

Tested on the HG002 BAM cut `HG002.GRCh38.2x250.pgx_v2.bam` (`32,091,666` B; SHA-256 `ca7d9b790a8350b1612aa1ad472599ff12001d390b841379e9a3055fa1858580`; `212,061` reads, `98.34%` mapped, `65.63x–75.98x` depth across 10 PGx v2 BED windows) paired with both `HG002_v421` and `HG002_v50q` VCFs + callable BEDs:

| Dataset | Runtime (s) | Genes | `CYP2D6` (`vcf+bam`) | `HLA-A` (`bam`) | `HLA-B` (`bam`) | `HLA-C` (`bam`) | Consensus Reference (`Aldy v4.8.3` / `T1K`) |
| --- | ---: | :---: | --- | --- | --- | --- | --- |
| `HG002_v421_v2` | 3.648 | 23 | `*2/*4` (`CN=2`, `AS=1.0`, `Intermediate Metabolizer`, `x_phased=true`) | `*01:01/*26:01` (`*31:01 negative`, `CN=2`) | `*35:08/*38:01` (`*15:02 negative; *57:01 negative; *58:01 negative`, `CN=2`) | `*04:01/*12:03` (`phenotype=null`, `CN=2`) | Exact match at 2-field resolution: `Aldy` `*2/*4` (`*2.001/*4.015`, `cpic_score=1.0`), `T1K` `HLA-A*26:01:01:01/*01:01:01:01`, `HLA-B*38:01:01:01/*35:08:01:01`, `HLA-C*12:03:01:01/*04:01:01:06` |
| `HG002_v50q_v2` | 3.710 | 23 | `*2/*4` (`CN=2`, `AS=1.0`, `Intermediate Metabolizer`, `x_phased=true`) | `*01:01/*26:01` (`*31:01 negative`, `CN=2`) | `*35:08/*38:01` (`*15:02 negative; *57:01 negative; *58:01 negative`, `CN=2`) | `*04:01/*12:03` (`phenotype=null`, `CN=2`) | Identical to `HG002_v421_v2`; also verified bit-for-bit `genes`/`drugs` determinism between `mmg.mmgen_pgx_run` and `mmg.pgx_run(..., engine="mmgen-pgx")` |

### 5.4 PharmCAT `-po` Round-Trip & IMGT/HLA Cache Modes

- **PharmCAT `-po` Round-Trip (`pharmcat-3.4.0-all.jar`, `7.785s`):** Passing `outside_calls.tsv` (`CYP2D6\t*2/*4\nHLA-A\t*01:01/*26:01\nHLA-B\t*35:08/*38:01\n`) via `mmg.pgx_run(..., engine="pharmcat", outside_calls=oc_path)` ingested all three loci cleanly with `callSource="OUTSIDE"` and matching diplotypes (`*2/*4`, `*01:01/*26:01`, `*35:08/*38:01`), with zero hard failures in `pgx_validation_results.json`.
- **Offline Guard:** Under `MMGEN_HLA_OFFLINE=1` with an empty `MMGEN_HLA_DIR`, `mmg.mmgen_pgx_run` raised `RuntimeError` reporting `HLA-A: missing cached IPD-IMGT/HLA file 'A_gen.fasta' of IPD-IMGT/HLA 3.66.0 (SHA-256 f265427f646d...) (MMGEN_HLA_OFFLINE=1)...`.
- **On-demand download of the pinned release (`MMGEN_TEST_IMGT_DOWNLOAD=1`):** with an empty `MMGEN_HLA_DIR`, `mmg.mmgen_pgx_run` downloaded `A_gen.fasta` (`20,391,880` B), `B_gen.fasta` (`24,329,147` B) and `C_gen.fasta` (`23,796,524` B) of IPD-IMGT/HLA 3.66.0 from the GitHub mirror into `<MMGEN_HLA_DIR>/3.66.0/`, checked their SHA-256, and produced the same `HLA-A/B/C` diplotypes (`*01:01/*26:01`, `*35:08/*38:01`, `*04:01/*12:03`); the whole run, downloads included, took **24.3 s**. The earlier download of the same files from `ftp.ebi.ac.uk` (unpinned) took **830.9 s**.
