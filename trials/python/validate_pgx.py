#!/usr/bin/env python3
"""Comprehensive validation suite for mmgen-pgx v1 against Section 7 acceptance criteria:
1. JSON Schema (Draft 2020-12) strict validation, determinism, HTML output, & clear errors.
2. Gene-by-gene comparison vs PharmCAT 3.4.0 across GIAB HG001-HG007 (9 datasets); the
   diplotype is chosen by score, as PharmCAT's Named Allele Matcher does, so a SAFETY_UNCALLED
   (PharmCAT calls, mmgen-pgx does not) is printed with its reason; every call with an
   alternative that changes the phenotype or a recommendation must carry a warning.
3. Negative "no silent calls" safety checks (key SNP not read -> called with a warning; nothing
   read, no callable_bed, missing chrX/M, GRCh37 -> not called).
4. Trio Mendelian inheritance consistency (Ashkenazi HG002/3/4 and Chinese HG005/6/7); a trio
   that fits only through an alternative the engine listed is reported as PASS_ALT.
5. GeT-RM concordance on NA12878 / HG001 (Pratt et al. 2016; Gaedigk et al. 2022).
6. CPU-only execution time check (< 120s per GIAB VCF).
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

import jsonschema
import multimerge_genomics as mmg

REPO = Path(__file__).resolve().parents[2]
# Every outside location comes from the environment, so that the script runs on any machine:
#   MMG_HUMAN_GENOMES     folder with giab/ laid out as the human-genomes download (default ./human_genomes)
#   MMGEN_PGX_TABLES      the mmgen-pgx tables (default mmgen_tools/pgx/cpic-v1.61.0 in the repository)
#   PHARMCAT_JAR          pharmcat-<version>-all.jar: the PharmCAT baseline is computed here with it
#   PGX_PHARMCAT_SUMMARY  or a baseline saved before ({name: {"calls": [...]}}, as pharmcat_run returns)
#   PGX_V2_BAM            the HG002 BAM cut for CYP2D6 and HLA (with its .bai)
#   MMGEN_HLA_DIR         folder with the IPD-IMGT/HLA A_gen.fasta, B_gen.fasta, C_gen.fasta
#   PGX_REPORT_OUT        where the detailed results are saved (default: the temporary folder)
TABLES_DIR = Path(os.environ.get("MMGEN_PGX_TABLES") or REPO / "mmgen_tools" / "pgx" / "cpic-v1.61.0")
SCHEMA_PATH = Path(getattr(mmg, "PGX_SCHEMA_PATH", REPO / "python" / "multimerge_genomics" / "pgx_result.schema.json"))
GIAB_ROOT = Path(os.environ.get("MMG_HUMAN_GENOMES", "human_genomes")) / "giab"
PHARMCAT_SUMMARY = Path(os.environ["PGX_PHARMCAT_SUMMARY"]) if os.environ.get("PGX_PHARMCAT_SUMMARY") else None
REPORT_JSON_OUT = Path(os.environ.get("PGX_REPORT_OUT") or Path(tempfile.gettempdir()) / "pgx_validation_results.json")
PGX_V2_BAM = Path(os.environ.get("PGX_V2_BAM", "pgx_v2/HG002.GRCh38.2x250.pgx_v2.bam"))
PGX_V2_HLA_CACHE = Path(os.environ.get("MMGEN_HLA_DIR", "pgx_v2/hla_cache"))
PHARMCAT_JAR = Path(os.environ.get("PHARMCAT_JAR", "mmgen_tools/pharmcat/pharmcat-3.4.0-all.jar"))

GIAB_CASES: list[tuple[str, Path, Path]] = [
    (
        "HG001_v421",
        GIAB_ROOT / "HG001_NA12878/v4.2.1/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        GIAB_ROOT / "HG001_NA12878/v4.2.1/HG001_GRCh38_1_22_v4.2.1_benchmark.bed",
    ),
    (
        "HG001_v332",
        GIAB_ROOT / "HG001_NA12878/v3.3.2_with_chrX/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz",
        GIAB_ROOT / "HG001_NA12878/v3.3.2_with_chrX/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed",
    ),
    (
        "HG002_v421",
        GIAB_ROOT / "HG002_NA24385_son/v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        GIAB_ROOT / "HG002_NA24385_son/v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
    ),
    (
        "HG002_v50q",
        GIAB_ROOT / "HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz",
        GIAB_ROOT / "HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.benchmark.bed",
    ),
    (
        "HG003_v421",
        GIAB_ROOT / "HG003_NA24149_father/v4.2.1/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        GIAB_ROOT / "HG003_NA24149_father/v4.2.1/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
    ),
    (
        "HG004_v421",
        GIAB_ROOT / "HG004_NA24143_mother/v4.2.1/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        GIAB_ROOT / "HG004_NA24143_mother/v4.2.1/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
    ),
    (
        "HG005_v421",
        GIAB_ROOT / "HG005_NA24631_son/v4.2.1/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        GIAB_ROOT / "HG005_NA24631_son/v4.2.1/HG005_GRCh38_1_22_v4.2.1_benchmark.bed",
    ),
    (
        "HG006_v421",
        GIAB_ROOT / "HG006_NA24694_father/v4.2.1/HG006_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        GIAB_ROOT / "HG006_NA24694_father/v4.2.1/HG006_GRCh38_1_22_v4.2.1_benchmark.bed",
    ),
    (
        "HG007_v421",
        GIAB_ROOT / "HG007_NA24695_mother/v4.2.1/HG007_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        GIAB_ROOT / "HG007_NA24695_mother/v4.2.1/HG007_GRCh38_1_22_v4.2.1_benchmark.bed",
    ),
]


def canon_dip(s: str | None) -> str:
    """Normalize display ordering/labels across mmgen-pgx and PharmCAT for semantic comparison."""
    if not s:
        return "NOT_CALLED"
    alts = [part.strip() for part in s.split(" OR ")]
    norm_alts = []
    for part in alts:
        if "/" in part:
            halves = sorted(h.strip() for h in part.split("/"))
            norm_alts.append("/".join(halves))
        else:
            norm_alts.append(part)
    return " OR ".join(sorted(norm_alts))


def check_schema(schema: dict[str, Any], instance: dict[str, Any], label: str) -> list[str]:
    validator = jsonschema.Draft202012Validator(schema, format_checker=jsonschema.FormatChecker())
    errs = sorted(validator.iter_errors(instance), key=lambda e: list(e.path))
    return [f"{label}: {'/'.join(map(str, e.path))}: {e.message}" for e in errs]


def candidate_allele_pairs(gene_obj: dict[str, Any]) -> list[tuple[str, str]] | None:
    """Return all candidate diploid haplotype pairs if the gene is called in a diploid autosomal locus."""
    if not gene_obj["called"]:
        return None
    if gene_obj.get("x_haplotype_1") is not None and gene_obj.get("x_haplotype_2") is not None:
        return [(gene_obj["x_haplotype_1"], gene_obj["x_haplotype_2"])]
    dip = gene_obj.get("diplotype") or ""
    if " OR " not in dip and len(gene_obj.get("alleles", [])) == 2 and "/" in dip:
        return [(gene_obj["alleles"][0], gene_obj["alleles"][1])]
    pairs = []
    for opt in dip.split(" OR "):
        opt = opt.strip()
        if "/" in opt:
            a, b = opt.split("/", 1)
            pairs.append((a.strip(), b.strip()))
        else:
            # e.g., variant_activity unphased "hapA AND hapB" or homozygous/ambiguous
            return None
    return pairs or None


def alternative_allele_pairs(gene_obj: dict[str, Any]) -> list[tuple[str, str]]:
    """The diploid pairs the engine listed in x_alternatives (lower score, or the other cis/trans)."""
    pairs = []
    for alt in gene_obj.get("x_alternatives") or []:
        dip = alt.get("diplotype") or ""
        if "/" in dip and " OR " not in dip:
            a, b = dip.split("/", 1)
            pairs.append((a.strip(), b.strip()))
    return pairs


def _trio_fit(sp, dp, mp) -> str | None:
    for s1, s2 in sp:
        for d1, d2 in dp:
            for m1, m2 in mp:
                if (s1 in (d1, d2) and s2 in (m1, m2)) or (s2 in (d1, d2) and s1 in (m1, m2)):
                    return f"child {s1}/{s2} <- father {d1}/{d2}, mother {m1}/{m2}"
    return None


def trio_compatible(son: dict[str, Any], dad: dict[str, Any], mom: dict[str, Any]) -> tuple[bool | None, str]:
    """Check whether son alleles can be inherited one from father and one from mother.

    True: the chosen calls fit. "alt": they do not, but the calls plus the alternatives the
    engine listed do (the family shows that the score picked the wrong one). False: nothing fits.
    """
    if not (son["called"] and dad["called"] and mom["called"]):
        return None, "not all trio members called"
    sp = candidate_allele_pairs(son)
    dp = candidate_allele_pairs(dad)
    mp = candidate_allele_pairs(mom)
    if sp is None or dp is None or mp is None:
        return None, "unphased compound 'AND' or non-diploid representation"
    fit = _trio_fit(sp, dp, mp)
    if fit:
        return True, fit
    fit = _trio_fit(sp + alternative_allele_pairs(son), dp + alternative_allele_pairs(dad), mp + alternative_allele_pairs(mom))
    if fit:
        return "alt", (f"chosen calls child={son['diplotype']}, father={dad['diplotype']}, mother={mom['diplotype']} "
                       f"do not fit; with the listed alternatives: {fit}")
    return False, f"child={son['diplotype']} incompatible with dad={dad['diplotype']} & mom={mom['diplotype']}"


def extract_pgx_subset_vcf(src_vcf: Path, dst_vcf: Path, drop_positions: set[tuple[str, int]] = frozenset()) -> None:
    """Extract only PGx chromosome regions using mmgen-pgx tables positions for negative testing."""
    import gzip

    opener = gzip.open if str(src_vcf).endswith(".gz") else open
    genes_data = json.loads((TABLES_DIR / "genes.json").read_text(encoding="utf-8"))
    windows: dict[str, list[tuple[int, int]]] = {}
    for g in genes_data["genes"]:
        c = g["chrom"].replace("chr", "")
        if c == "M":
            c = "MT"
        for p in g["positions"]:
            windows.setdefault(c, []).append((p["pos"] - 100, p["pos"] + 100))

    with opener(src_vcf, "rt", encoding="utf-8", errors="replace") as fin, open(dst_vcf, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            parts = line.split("\t", 5)
            if len(parts) < 5:
                continue
            ch, pos_s = parts[0], parts[1]
            c = ch[3:] if ch.startswith("chr") else ch
            if c == "M":
                c = "MT"
            pos = int(pos_s)
            if (ch, pos) in drop_positions or (c, pos) in drop_positions:
                continue
            if any(s <= pos <= e for s, e in windows.get(c, [])):
                fout.write(line)


def main() -> None:
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    missing = [str(c[1]) for c in GIAB_CASES if not c[1].exists()]
    if missing:
        sys.exit(f"GIAB files not found under {GIAB_ROOT} (set MMG_HUMAN_GENOMES): {missing[:2]}...")
    if PHARMCAT_SUMMARY is not None:
        pc_summary = json.loads(PHARMCAT_SUMMARY.read_text(encoding="utf-8"))
    elif PHARMCAT_JAR.exists():
        print(f"PharmCAT baseline: running {PHARMCAT_JAR.name} on the {len(GIAB_CASES)} GIAB sets...")
        pc_summary = {}
        with tempfile.TemporaryDirectory(prefix="mmgen_pgx_pc_") as pc_tmp:
            for name, vcf_path, bed_path in GIAB_CASES:
                run = mmg.pharmcat_run(str(vcf_path), str(Path(pc_tmp) / name), str(PHARMCAT_JAR), callable_bed=str(bed_path))
                pc_summary[name] = {"calls": run["calls"]}
    else:
        print("PharmCAT baseline: none (set PHARMCAT_JAR or PGX_PHARMCAT_SUMMARY); the comparison is skipped.")
        pc_summary = {name: {"calls": []} for name, _, _ in GIAB_CASES}
    tally: dict[str, int] = {}
    tables_str = str(TABLES_DIR)
    all_failures: list[str] = []
    safety_uncalled: list[str] = []
    changing_alternatives: list[str] = []

    with tempfile.TemporaryDirectory(prefix="mmgen_pgx_val_") as tmp_str:
        tmp = Path(tmp_str)
        giab_results: dict[str, Any] = {}
        timings: dict[str, float] = {}

        # 1 + 2 + 6: Run all 9 GIAB datasets, validate schema, compare to PharmCAT, measure runtime
        print("=" * 80)
        print("SECTION 1: GIAB BENCHMARKS (9 DATASETS) — SCHEMA, TIMING & PHARMCAT COMPARISON")
        print("=" * 80)
        for name, vcf_path, bed_path in GIAB_CASES:
            out_dir = tmp / name
            t0 = time.perf_counter()
            py_res = mmg.mmgen_pgx_run(str(vcf_path), str(out_dir), callable_bed=str(bed_path), tables=tables_str)
            elapsed = time.perf_counter() - t0
            timings[name] = round(elapsed, 3)

            if "no responsibility" not in py_res.get("notice", ""):
                all_failures.append(f"{name}: missing notice on returned Python dict")

            disk_path = Path(py_res["x_files"]["result_json"])
            html_path = Path(py_res["report_html"])
            disk_bytes_1 = disk_path.read_bytes()
            disk_obj = json.loads(disk_bytes_1.decode("utf-8"))

            schema_errs = check_schema(schema, disk_obj, name)
            if schema_errs:
                all_failures.extend(schema_errs)
            if not html_path.exists() or "mmgen-pgx" not in html_path.read_text(encoding="utf-8"):
                all_failures.append(f"{name}: missing or invalid report.html")
            if elapsed >= 120.0:
                all_failures.append(f"{name}: exceeded 120s budget ({elapsed:.2f}s)")

            # Test pgx_run wrapper too on first dataset + determinism check on HG002_v50q
            if name == "HG002_v50q":
                py_res_2 = mmg.pgx_run(
                    str(vcf_path), str(out_dir), callable_bed=str(bed_path), engine="mmgen-pgx", tables=tables_str
                )
                disk_bytes_2 = Path(py_res_2["x_files"]["result_json"]).read_bytes()
                if disk_bytes_1 != disk_bytes_2:
                    all_failures.append(f"{name}: non-deterministic result.json bytes across repeated runs")

            giab_results[name] = disk_obj
            pc_by_gene = {c["gene"]: c for c in pc_summary[name]["calls"]}
            my_by_gene = {g["gene"]: g for g in disk_obj["genes"]}

            called_n = disk_obj["x_summary"]["genes_called"]
            drugs_n = disk_obj["x_summary"]["drugs_with_recommendation"]
            print(f"\n[{name}] time={elapsed:.3f}s  called={called_n}/19  drugs={drugs_n}  schema={'OK' if not schema_errs else f'ERR({len(schema_errs)})'}")
            for gname in sorted(my_by_gene):
                mg = my_by_gene[gname]
                pc = pc_by_gene.get(gname, {})
                pc_dip = pc.get("source_diplotype")
                pc_ph = pc.get("phenotype")
                my_dip = mg["diplotype"]
                my_ph = mg["phenotype"]

                # Determine comparison bucket
                my_called = mg["called"]
                pc_called = bool(pc_dip and pc_dip != "NOT_CALLED")
                exact = my_called and pc_called and canon_dip(my_dip) == canon_dip(pc_dip)
                if not my_called and not pc_called:
                    status = "BOTH_UNCALLED"
                elif exact:
                    status = "EXACT"
                elif not my_called and pc_called:
                    status = "SAFETY_UNCALLED"
                elif my_called and not pc_called:
                    status = "MMGEN_CALLED_ONLY"
                else:
                    status = "FORMAT_OR_CALL_DIFF"
                tally[status] = tally.get(status, 0) + 1
                if status != "EXACT" or name in ("HG001_v421", "HG002_v50q"):
                    print(
                        f"  {gname:9s} | {status:18s} | mmgen={str(my_dip):34s} [{str(my_ph):25s}]"
                        f" | pc={str(pc_dip):30s} [{str(pc_ph)}]"
                    )
                if status == "SAFETY_UNCALLED":
                    print(f"  {'':9s}   reason: {mg['not_called_reason']}")
                    safety_uncalled.append(f"{name} {gname}: {mg['not_called_reason']}")
                changing = [a for a in mg.get("x_alternatives") or [] if a["changes_phenotype"] or a["changes_recommendation"]]
                if mg["called"] and changing:
                    alts = "; ".join(f"{a['diplotype']} ({a['phenotype']})" for a in changing[:4])
                    changing_alternatives.append(f"{name} {gname}: {my_dip} ({my_ph}); also possible {alts}")
                    if not any(w.startswith(f"{gname}: ") and "also possible" in w for w in disk_obj["warnings"]):
                        all_failures.append(f"{name} {gname}: an alternative changes the phenotype/recommendation but no warning")

            if name == "HG001_v421" and my_by_gene["CYP2C19"]["diplotype"] != "*1/*2":
                all_failures.append(f"HG001_v421 CYP2C19 = {my_by_gene['CYP2C19']['diplotype']}, expected *1/*2 (GeT-RM, PharmCAT)")

        print("\nmmgen-pgx x PharmCAT, gene by gene over the GIAB sets: "
              + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())))
        print("  (SAFETY_UNCALLED = PharmCAT calls the gene and mmgen-pgx leaves it uncalled, with its reason)")
        print(f"Calls with an alternative that would change the phenotype or a recommendation: {len(changing_alternatives)}")
        for line in changing_alternatives:
            print("  -", line)

        # 3: Negative safety checks
        print("\n" + "=" * 80)
        print("SECTION 2: NEGATIVE 'NO SILENT CALL' SAFETY CHECKS")
        print("=" * 80)

        # 3a. No callable_bed on variant-only GIAB VCF -> genes absent from VCF must NOT silently call *1
        out_nobed = tmp / "neg_nobed"
        vcf_hg002, _ = GIAB_CASES[3][1], GIAB_CASES[3][2]
        nobed_res = mmg.mmgen_pgx_run(str(vcf_hg002), str(out_nobed), callable_bed=None, tables=tables_str)
        nobed_disk = json.loads(Path(nobed_res["x_files"]["result_json"]).read_text(encoding="utf-8"))
        nobed_schema_errs = check_schema(schema, nobed_disk, "HG002_nobed")
        all_failures.extend(nobed_schema_errs)
        nobed_by_g = {g["gene"]: g for g in nobed_res["genes"]}
        uncalled_nobed = [g for g, o in nobed_by_g.items() if not o["called"]]
        print(f"  callable_bed=None on HG002_v50q: {len(uncalled_nobed)}/19 uncalled ({', '.join(uncalled_nobed[:8])}...)")
        if len(uncalled_nobed) == 0:
            all_failures.append("callable_bed=None failed to guard sparse VCF positions")
        # Specifically VKORC1 / TPMT / ABCG2 / CACNA1S / RYR1 without variants shouldn't be callable without BED
        for must_be_uncalled in ("TPMT", "MT-RNR1"):
            if nobed_by_g[must_be_uncalled]["called"]:
                all_failures.append(f"callable_bed=None silently called {must_be_uncalled}")

        # 3b. Missing chromosome check on HG002_v421 (1_22 only -> G6PD chrX & MT-RNR1 chrM uncalled)
        hg002_421_g = {g["gene"]: g for g in giab_results["HG002_v421"]["genes"]}
        for g_miss_seq in ("G6PD", "MT-RNR1"):
            if hg002_421_g[g_miss_seq]["called"]:
                all_failures.append(f"HG002_v421 missing {g_miss_seq} chromosome still called!")
        print(f"  Missing chrX/chrM on HG002_v421: G6PD reason={hg002_421_g['G6PD']['not_called_reason']!r}")
        print(f"                                   MT-RNR1 reason={hg002_421_g['MT-RNR1']['not_called_reason']!r}")

        # 3c. Delete key SNPs from an extracted HG002_v50q PGx subset (CYP2C19 rs4244285 chr10:94781859, VKORC1 rs9923231 chr16:31096368)
        sub_drop = tmp / "hg002_drop_key.vcf"
        extract_pgx_subset_vcf(vcf_hg002, sub_drop, drop_positions={("chr10", 94781859), ("chr16", 31096368)})
        # Remove key coords from a tiny callable bed too so they're truly missing (not callable_ref!)
        sub_bed = tmp / "hg002_drop_key.bed"
        with open(GIAB_CASES[3][2], encoding="utf-8") as fb_in, open(sub_bed, "w", encoding="utf-8") as fb_out:
            for line in fb_in:
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                ch, s, e = parts[0], int(parts[1]), int(parts[2])
                if ch == "chr10" and s < 94781859 <= e:
                    if s < 94781858:
                        fb_out.write(f"chr10\t{s}\t94781858\n")
                    if 94781860 < e:
                        fb_out.write(f"chr10\t94781860\t{e}\n")
                elif ch == "chr16" and s < 31096368 <= e:
                    if s < 31096367:
                        fb_out.write(f"chr16\t{s}\t31096367\n")
                    if 31096369 < e:
                        fb_out.write(f"chr16\t31096369\t{e}\n")
                elif ch in ("chr1", "chr2", "chr4", "chr6", "chr7", "chr8", "chr10", "chr12", "chr13", "chr16", "chr19", "chr22", "chrX", "chrM"):
                    fb_out.write(line)
        drop_res = mmg.mmgen_pgx_run(str(sub_drop), str(tmp / "out_drop"), callable_bed=str(sub_bed), tables=tables_str)
        drop_by_g = {g["gene"]: g for g in drop_res["genes"]}
        # VKORC1 has one position: nothing read -> not called. CYP2C19 is called from the other
        # positions (as PharmCAT does), and must say, in a warning and in x_missing_key_positions,
        # that rs4244285 was not read: a call, but never a silent one.
        if drop_by_g["VKORC1"]["called"]:
            all_failures.append("Dropped key SNP: VKORC1 (its only position) still called")
        cyp = drop_by_g["CYP2C19"]
        cyp_warned = any(w.startswith("CYP2C19: key position(s) not read") and "rs4244285" in w for w in drop_res["warnings"])
        if cyp["called"] and not (cyp_warned and "rs4244285" in (cyp.get("x_missing_key_positions") or [])):
            all_failures.append("Dropped key SNP: CYP2C19 called without saying that rs4244285 was not read")
        for dropped_gene in ("CYP2C19", "VKORC1"):
            g = drop_by_g[dropped_gene]
            print(f"  Dropped key SNP {dropped_gene}: called={g['called']}, diplotype={g['diplotype']}, "
                  f"missing keys={g.get('x_missing_key_positions')}, reason={g['not_called_reason']}")
        print(f"  CYP2C19 warning present: {cyp_warned}")

        # 4: Trio Mendelian Consistency
        print("\n" + "=" * 80)
        print("SECTION 3: TRIO MENDELIAN CONSISTENCY (ASHKENAZI HG002/3/4 & CHINESE HG005/6/7)")
        print("=" * 80)
        trios = [
            ("Ashkenazi_Trio", "HG002_v421", "HG003_v421", "HG004_v421"),
            ("Chinese_Trio", "HG005_v421", "HG006_v421", "HG007_v421"),
        ]
        trio_summary: dict[str, Any] = {}
        trio_by_alternative: list[str] = []
        for tname, s_k, d_k, m_k in trios:
            sg = {g["gene"]: g for g in giab_results[s_k]["genes"]}
            dg = {g["gene"]: g for g in giab_results[d_k]["genes"]}
            mg = {g["gene"]: g for g in giab_results[m_k]["genes"]}
            trio_rows = []
            for gsym in sorted(sg):
                ok, detail = trio_compatible(sg[gsym], dg[gsym], mg[gsym])
                trio_rows.append({"gene": gsym, "compatible": ok, "detail": detail})
                if ok is False:
                    all_failures.append(f"{tname} {gsym} Mendelian conflict: {detail}")
                elif ok == "alt":
                    trio_by_alternative.append(f"{tname} {gsym}: {detail}")
                tag = {True: "PASS", None: "SKIP", "alt": "PASS_ALT"}.get(ok, "CONFLICT")
                print(f"  [{tname}] {gsym:9s} | {tag:8s} | {detail}")
            trio_summary[tname] = trio_rows
        if trio_by_alternative:
            print("  PASS_ALT = the chosen calls contradict the family; an alternative the engine listed fits:")
            for line in trio_by_alternative:
                print("    -", line)

        # 5: GeT-RM NA12878 Concordance (Pratt 2016 + Gaedigk 2022)
        print("\n" + "=" * 80)
        print("SECTION 4: GeT-RM NA12878 / HG001 CONCORDANCE")
        print("=" * 80)
        getrm_na12878 = {
            "CYP2B6": {"getrm": "*1/*1", "source": "Pratt 2016"},
            "CYP2C9": {"getrm": "*1/*2", "source": "Pratt 2016; Gaedigk 2022"},
            "CYP2C19": {"getrm": "*1/*2", "source": "Pratt 2016; Gaedigk 2022"},
            "CYP3A5": {"getrm": "*3/*3", "source": "Pratt 2016"},
            "CYP4F2": {"getrm": "*1/*1", "source": "Pratt 2016"},
            "DPYD": {"getrm": "*1/(*4)", "source": "Pratt 2016"},
            "NAT2": {"getrm": "*4/*5", "source": "Pratt 2016"},
            "SLCO1B1": {"getrm": "*1/*15", "source": "Pratt 2016"},
            "TPMT": {"getrm": "*1/*1", "source": "Pratt 2016"},
            "UGT1A1": {"getrm": "*60/*28", "source": "Pratt 2016"},
            "VKORC1": {"getrm": "rs9923231 GA (-1639G>A het)", "source": "Pratt 2016"},
        }
        h1_421 = {g["gene"]: g for g in giab_results["HG001_v421"]["genes"]}
        h1_332 = {g["gene"]: g for g in giab_results["HG001_v332"]["genes"]}
        for gsym, meta in getrm_na12878.items():
            g4 = h1_421[gsym]
            g3 = h1_332[gsym]
            print(
                f"  {gsym:8s} | GeT-RM={meta['getrm']:28s} | v4.2.1={str(g4['diplotype']):32s} "
                f"| v3.3.2={str(g3['diplotype'])}"
            )

        # 6: Traits validation (schema, no-silent-calls, GIAB genotypes, trio Mendelian inheritance)
        print("\n" + "=" * 80)
        print("SECTION 5: GENETIC TRAITS (LACTASE + EYE COLOR 6-SNP PANEL)")
        print("=" * 80)
        for name in [c[0] for c in GIAB_CASES]:
            tr_map = {t["x_id"]: t for t in giab_results[name]["traits"]}
            if set(tr_map.keys()) != {"lactase_persistence", "eye_color"}:
                all_failures.append(f"{name}: expected 2 traits, found {list(tr_map.keys())}")
                continue
            lp = tr_map["lactase_persistence"]
            ec = tr_map["eye_color"]
            print(f"  {name:11s} | Lactase: {lp['interpretation']} [{lp['genotype']}]")
            print(f"              | Eye    : {ec['interpretation']} [{ec['genotype']}]")
        # Mendelian trio check per trait SNP across Ashkenazi & Chinese trios
        for tname, s_k, d_k, m_k in trios:
            for tid in ("lactase_persistence", "eye_color"):
                s_ev = {e["rsid"]: e["genotype"] for t in giab_results[s_k]["traits"] if t["x_id"] == tid for e in t["x_evidence"]}
                d_ev = {e["rsid"]: e["genotype"] for t in giab_results[d_k]["traits"] if t["x_id"] == tid for e in t["x_evidence"]}
                m_ev = {e["rsid"]: e["genotype"] for t in giab_results[m_k]["traits"] if t["x_id"] == tid for e in t["x_evidence"]}
                for rsid in sorted(s_ev):
                    sg, dg, mg = s_ev[rsid], d_ev[rsid], m_ev[rsid]
                    if "." in sg or "." in dg or "." in mg:
                        continue
                    sa = sg.replace("|", "/").split("/")
                    da = set(dg.replace("|", "/").split("/"))
                    ma = set(mg.replace("|", "/").split("/"))
                    ok = (sa[0] in da and sa[1] in ma) or (sa[1] in da and sa[0] in ma)
                    if not ok:
                        all_failures.append(f"{tname} {tid} {rsid} Mendelian conflict: child={sg} father={dg} mother={mg}")
        # Negative no-silent-call check on nobed for Lactase (absent from variant-only VCF)
        nobed_tr = {t["x_id"]: t for t in nobed_res["traits"]}
        if nobed_tr["lactase_persistence"]["x_called"]:
            all_failures.append("HG002_nobed: lactase_persistence called without callable_bed!")
        print(f"  HG002_nobed | Lactase x_called={nobed_tr['lactase_persistence']['x_called']} ({nobed_tr['lactase_persistence']['interpretation']})")

        print("\n" + "=" * 80)
        print("SECTION 6: PGX v2 BAM+VCF (CYP2D6, CLASS I HLA, PHARMCAT -po & ON-DEMAND IMGT)")
        print("=" * 80)
        v2_summary: dict[str, Any] = {}
        if PGX_V2_BAM.exists() and PGX_V2_HLA_CACHE.exists():
            saved_dir = os.environ.get("MMGEN_HLA_DIR")
            saved_off = os.environ.get("MMGEN_HLA_OFFLINE")
            try:
                os.environ["MMGEN_HLA_DIR"] = str(PGX_V2_HLA_CACHE)
                os.environ["MMGEN_HLA_OFFLINE"] = "1"
                hg002_cases = [c for c in GIAB_CASES if c[0] in ("HG002_v50q", "HG002_v421") and c[1].exists()]
                for ds_name, ds_vcf, ds_bed in hg002_cases:
                    t0 = time.perf_counter()
                    v2_out = tmp / f"{ds_name}_v2"
                    v2_res = mmg.mmgen_pgx_run(
                        str(ds_vcf),
                        str(v2_out),
                        callable_bed=str(ds_bed),
                        tables=tables_str,
                        bam=str(PGX_V2_BAM),
                    )
                    dt_v2 = time.perf_counter() - t0
                    disk_v2 = json.loads((v2_out / "result.json").read_text(encoding="utf-8"))
                    all_failures.extend(check_schema(schema, disk_v2, f"{ds_name}_v2"))
                    if len(disk_v2["genes"]) != 23:
                        all_failures.append(f"{ds_name}_v2: expected 23 genes (19 VCF + 4 BAM), got {len(disk_v2['genes'])}")
                    by_g2 = {g["gene"]: g for g in disk_v2["genes"]}
                    cyp = by_g2.get("CYP2D6", {})
                    hla_a = by_g2.get("HLA-A", {})
                    hla_b = by_g2.get("HLA-B", {})
                    hla_c = by_g2.get("HLA-C", {})
                    expected_v2 = [
                        ("CYP2D6", cyp, "*2/*4", "vcf+bam", 2),
                        ("HLA-A", hla_a, "*01:01/*26:01", "bam", 2),
                        ("HLA-B", hla_b, "*35:08/*38:01", "bam", 2),
                        ("HLA-C", hla_c, "*04:01/*12:03", "bam", 2),
                    ]
                    for gsym, gobj, want_dip, want_src, want_cn in expected_v2:
                        if gobj.get("diplotype") != want_dip or gobj.get("source") != want_src or gobj.get("copy_number") != want_cn:
                            all_failures.append(
                                f"{ds_name}_v2 {gsym}: got ({gobj.get('diplotype')}, src={gobj.get('source')}, cn={gobj.get('copy_number')}), want ({want_dip}, {want_src}, {want_cn})"
                            )
                    oc_path = v2_res["x_files"].get("pharmcat_outside_calls", "")
                    oc_txt = Path(oc_path).read_text(encoding="utf-8") if oc_path and Path(oc_path).exists() else ""
                    if oc_txt != "CYP2D6\t*2/*4\nHLA-A\t*01:01/*26:01\nHLA-B\t*35:08/*38:01\n":
                        all_failures.append(f"{ds_name}_v2 outside_calls.tsv mismatch: {oc_txt!r}")

                    print(
                        f"  {ds_name}_v2 | {dt_v2:.3f}s | genes={len(disk_v2['genes'])} | "
                        f"CYP2D6={cyp.get('diplotype')} (CN={cyp.get('copy_number')}, AS={cyp.get('activity_score')}) | "
                        f"HLA-A={hla_a.get('diplotype')} | HLA-B={hla_b.get('diplotype')} | HLA-C={hla_c.get('diplotype')}"
                    )
                    v2_summary[ds_name] = {
                        "runtime_sec": round(dt_v2, 3),
                        "genes_count": len(disk_v2["genes"]),
                        "CYP2D6": {k: cyp.get(k) for k in ("diplotype", "phenotype", "activity_score", "copy_number", "source", "x_phased")},
                        "HLA-A": {k: hla_a.get(k) for k in ("diplotype", "phenotype", "copy_number", "source")},
                        "HLA-B": {k: hla_b.get(k) for k in ("diplotype", "phenotype", "copy_number", "source")},
                        "HLA-C": {k: hla_c.get(k) for k in ("diplotype", "phenotype", "copy_number", "source")},
                        "outside_calls_tsv": oc_txt,
                    }

                    if ds_name == "HG002_v50q":
                        v2_wrap = mmg.pgx_run(
                            str(ds_vcf),
                            str(tmp / f"{ds_name}_v2_wrap"),
                            callable_bed=str(ds_bed),
                            engine="mmgen-pgx",
                            tables=tables_str,
                            bam=str(PGX_V2_BAM),
                        )
                        if v2_wrap["genes"] != v2_res["genes"] or v2_wrap["drugs"] != v2_res["drugs"]:
                            all_failures.append("HG002_v50q_v2: non-deterministic genes/drugs between mmgen_pgx_run and pgx_run")

                        if PHARMCAT_JAR.exists() and oc_path:
                            t_pc = time.perf_counter()
                            pc_out = tmp / "hg002_pharmcat_po"
                            pc_res = mmg.pgx_run(
                                str(ds_vcf),
                                str(pc_out),
                                callable_bed=str(ds_bed),
                                engine="pharmcat",
                                jar=str(PHARMCAT_JAR),
                                outside_calls=oc_path,
                            )
                            dt_pc = time.perf_counter() - t_pc
                            # pgx_run gives PharmCAT's result in the same schema as mmgen-pgx's.
                            all_failures.extend(check_schema(schema, pc_res, "PharmCAT via pgx_run"))
                            pc_common = {g["gene"]: g for g in pc_res["genes"]}
                            if pc_common.get("CYP2D6", {}).get("source") != "outside" or pc_common.get("CYP2D6", {}).get("diplotype") != "*2/*4":
                                all_failures.append(f"PharmCAT via pgx_run: CYP2D6 {pc_common.get('CYP2D6')}")
                            pc_rpt = json.loads(Path(pc_res["x_pharmcat"]["files"]["report_json"]).read_text(encoding="utf-8"))
                            pc_po_calls: dict[str, Any] = {}
                            for gsym, want_dip in [("CYP2D6", "*2/*4"), ("HLA-A", "*01:01/*26:01"), ("HLA-B", "*35:08/*38:01")]:
                                gpc = pc_rpt.get("genes", {}).get("CPIC", {}).get(gsym) or pc_rpt.get("genes", {}).get(gsym, {})
                                csrc = gpc.get("callSource")
                                sdips = [x.get("label") for x in gpc.get("sourceDiplotypes", [])]
                                pc_po_calls[gsym] = {"callSource": csrc, "diplotypes": sdips}
                                if csrc != "OUTSIDE" or want_dip not in sdips:
                                    all_failures.append(f"PharmCAT -po {gsym}: got callSource={csrc}, diplotypes={sdips}, want OUTSIDE / {want_dip}")
                            v2_summary["pharmcat_po_roundtrip"] = {"runtime_sec": round(dt_pc, 3), "genes": pc_po_calls}
                            print(f"  PharmCAT -po round-trip | {dt_pc:.3f}s | {pc_po_calls}")

                if hg002_cases:
                    _, ref_vcf, ref_bed = hg002_cases[0]
                    # Offline missing-cache guard check
                    empty_hla_dir = tmp / "hla_empty_guard"
                    empty_hla_dir.mkdir(parents=True, exist_ok=True)
                    os.environ["MMGEN_HLA_DIR"] = str(empty_hla_dir)
                    os.environ["MMGEN_HLA_OFFLINE"] = "1"
                    try:
                        mmg.mmgen_pgx_run(
                            str(ref_vcf),
                            str(tmp / "hg002_v2_offline_err"),
                            callable_bed=str(ref_bed),
                            tables=tables_str,
                            bam=str(PGX_V2_BAM),
                        )
                        all_failures.append("MMGEN_HLA_OFFLINE=1 with missing FASTAs did not raise error")
                    except RuntimeError as exc:
                        if "MMGEN_HLA_OFFLINE" not in str(exc) or "A_gen.fasta" not in str(exc):
                            all_failures.append(f"Unexpected offline HLA error: {exc}")
                        else:
                            print("  Offline HLA guard OK:", str(exc)[:110])

                    # On-demand IMGT/HLA FASTA download (MMGEN_TEST_IMGT_DOWNLOAD=1): the engine fetches
                    # the release the tables pin (bam_genes.json) from the GitHub mirror into
                    # <MMGEN_HLA_DIR>/<release>/ and checks each file's SHA-256 against the tables.
                    if os.environ.get("MMGEN_TEST_IMGT_DOWNLOAD") == "1":
                        hla_pins = {h["fasta_filename"]: h for h in json.loads((TABLES_DIR / "bam_genes.json").read_text(encoding="utf-8"))["hla_genes"]}
                        hla_release = hla_pins["A_gen.fasta"].get("imgt_hla_version") or ""
                        dl_hla_dir = tmp / "hla_ondemand_dl"
                        os.environ["MMGEN_HLA_DIR"] = str(dl_hla_dir)
                        os.environ.pop("MMGEN_HLA_OFFLINE", None)
                        t_dl = time.perf_counter()
                        dl_res = mmg.mmgen_pgx_run(
                            str(ref_vcf),
                            str(tmp / "hg002_v2_dl"),
                            callable_bed=str(ref_bed),
                            tables=tables_str,
                            bam=str(PGX_V2_BAM),
                        )
                        dt_dl = time.perf_counter() - t_dl
                        dl_dir = dl_hla_dir / hla_release
                        dl_files = {fn: (dl_dir / fn).stat().st_size for fn in hla_pins if (dl_dir / fn).exists()}
                        dl_g = {g["gene"]: g["diplotype"] for g in dl_res["genes"]}
                        if len(dl_files) != 3 or any(sz < 100_000 for sz in dl_files.values()):
                            all_failures.append(f"On-demand IMGT download incomplete: {dl_files}")
                        for fn, pin in hla_pins.items():
                            if (dl_dir / fn).exists() and pin.get("fasta_sha256") and hashlib.sha256((dl_dir / fn).read_bytes()).hexdigest() != pin["fasta_sha256"]:
                                all_failures.append(f"On-demand IMGT download: {fn} is not IPD-IMGT/HLA {hla_release}")
                        if f"IPD-IMGT/HLA {hla_release}" not in dl_res["data_version"]:
                            all_failures.append(f"data_version does not name IPD-IMGT/HLA {hla_release}: {dl_res['data_version']}")
                        if dl_g.get("HLA-A") != "*01:01/*26:01" or dl_g.get("HLA-B") != "*35:08/*38:01" or dl_g.get("HLA-C") != "*04:01/*12:03":
                            all_failures.append(f"On-demand IMGT download calls mismatch: {dl_g}")
                        v2_summary["imgt_ondemand_download"] = {"runtime_sec": round(dt_dl, 3), "release": hla_release, "files_bytes": dl_files}
                        print(f"  On-demand IMGT download | {dt_dl:.3f}s | {dl_files}")
                    else:
                        print("  NOTE skipping the live IPD-IMGT/HLA download from GitHub (set MMGEN_TEST_IMGT_DOWNLOAD=1 to run it)")
            finally:
                if saved_dir is None:
                    os.environ.pop("MMGEN_HLA_DIR", None)
                else:
                    os.environ["MMGEN_HLA_DIR"] = saved_dir
                if saved_off is None:
                    os.environ.pop("MMGEN_HLA_OFFLINE", None)
                else:
                    os.environ["MMGEN_HLA_OFFLINE"] = saved_off
        else:
            print("  NOTE skipping HG002 BAM v2 check (local BAM / HLA cache not present)")

        out_summary = {
            "timings_sec": timings,
            "giab_results": {
                k: {
                    "summary": v["x_summary"],
                    "traits": v["traits"],
                    "genes": {g["gene"]: {fld: g[fld] for fld in ("called", "diplotype", "alleles", "phenotype", "activity_score", "not_called_reason", "x_phased", "positions_read", "x_positions_total", "notes") if fld in g} for g in v["genes"]},
                    "pharmcat": {c["gene"]: {"source_diplotype": c.get("source_diplotype"), "phenotype": c.get("phenotype"), "activity_score": c.get("activity_score")} for c in pc_summary[k]["calls"]},
                }
                for k, v in giab_results.items()
            },
            "pharmcat_comparison": tally,
            "safety_uncalled": safety_uncalled,
            "calls_with_changing_alternatives": changing_alternatives,
            "trio_summary": trio_summary,
            "trio_resolved_by_alternative": trio_by_alternative,
            "v2_bam": v2_summary,
            "failures": all_failures,
        }
        REPORT_JSON_OUT.write_text(json.dumps(out_summary, indent=2), encoding="utf-8")

    print("\n" + "=" * 80)
    if all_failures:
        print(f"VALIDATION FAILED ({len(all_failures)} failure(s)):")
        for f in all_failures:
            print("  -", f)
        sys.exit(1)
    else:
        print("VALIDATION SUCCEEDED: 0 hard failures! Saved detailed results to", REPORT_JSON_OUT)


if __name__ == "__main__":
    main()
