












































#![allow(dead_code)]
//! Standalone Rust tests for `src/farmacogenetica.rs` (`mmgen-pgx` v1).
//!
//! Usage:
//!   cargo run --release --manifest-path trials/Cargo.toml --bin farmacogenetica

#[path = "../../../src/multimerge_core.rs"]
pub mod multimerge_core;
#[path = "../../../src/bgzf.rs"]
#[allow(dead_code)]
mod bgzf;
#[path = "../../../src/bam.rs"]
#[allow(dead_code)]
mod bam;
#[path = "../../../src/externo.rs"]
#[allow(dead_code)]
mod externo;
#[path = "../../../src/nomes.rs"]
#[allow(dead_code)]
mod nomes;
#[path = "../../../src/indice.rs"]
#[allow(dead_code)]
mod indice;
#[path = "../../../src/referencia.rs"]
#[allow(dead_code)]
mod referencia;
#[path = "../../../src/farmacogenetica.rs"]
mod farmacogenetica;

use farmacogenetica::*;
use std::io::Write;
use std::path::{Path, PathBuf};

struct Harness {
    ok: usize,
    fail: usize,
    dir: PathBuf,
}

impl Harness {
    fn new(name: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("mmgen_pgx_trial_{name}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Self { ok: 0, fail: 0, dir }
    }

    fn p(&self, rel: &str) -> String {
        self.dir.join(rel).to_string_lossy().to_string()
    }

    fn check(&mut self, cond: bool, what: &str) -> bool {
        if cond {
            self.ok += 1;
        } else {
            self.fail += 1;
            println!("FAIL  {what}");
        }
        cond
    }

    fn eq_str(&mut self, got: &str, want: &str, what: &str) -> bool {
        let ok = got == want;
        if !ok {
            println!("FAIL  {what}\n  got:  {got:?}\n  want: {want:?}");
        }
        self.check(ok, what)
    }

    fn expect_err<T>(&mut self, res: Result<T, String>, needles: &[&str], what: &str) {
        match res {
            Ok(_) => {
                println!("FAIL  {what}: expected error, succeeded");
                self.check(false, what);
            }
            Err(msg) => {
                let missing: Vec<&&str> = needles.iter().filter(|n| !msg.contains(**n)).collect();
                if !missing.is_empty() {
                    println!("FAIL  {what}: missing {missing:?} in error:\n  {msg}");
                }
                self.check(missing.is_empty(), what);
            }
        }
    }
}

impl Drop for Harness {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

fn write_text(path: &str, body: &str) {
    if let Some(parent) = Path::new(path).parent() {
        std::fs::create_dir_all(parent).unwrap();
    }
    let mut w = externo::Escritor::cria(path).unwrap();
    w.write_all(body.as_bytes()).unwrap();
    w.termina().unwrap();
}

fn copy_dir(src: &Path, dst: &Path) {
    std::fs::create_dir_all(dst).unwrap();
    for entry in std::fs::read_dir(src).unwrap() {
        let e = entry.unwrap();
        let p = e.path();
        if p.is_file() {
            std::fs::copy(&p, dst.join(p.file_name().unwrap())).unwrap();
        }
    }
}

fn gene_by_name<'a>(v: &'a serde_json::Value, sym: &str) -> &'a serde_json::Value {
    v["genes"]
        .as_array()
        .unwrap()
        .iter()
        .find(|g| g["gene"].as_str() == Some(sym))
        .unwrap_or_else(|| panic!("missing gene {sym}"))
}

/// Copies the fixture tables and adds MINI2C19, shaped like CYP2C19 in the CPIC tables: the
/// reference sequence is *38, *1 differs from it at one position (rs92 G), *2 carries IUPAC
/// R at rs90 and rs92, *35 is rs90 G + rs92 G. Then rewrites the sizes and SHA-256 sums.
fn add_mini_cyp2c19(src: &Path, dst: &Path) {
    use serde_json::{json, Value};
    copy_dir(src, dst);
    let read_json = |name: &str| -> Value { serde_json::from_slice(&std::fs::read(dst.join(name)).unwrap()).unwrap() };
    let pos = |id: usize, p: u64, rs: &str, r: &str, a: &str| {
        json!({"id": id, "chrom": "chr10", "pos": p, "cpic_pos": p, "rsid": rs, "hgvs": null, "ref_cpic": r,
               "grch38_ref": r, "is_key": true,
               "alt_variants": [{"cpic_value": r, "vcf_pos": p, "vcf_ref": r, "vcf_alt": r, "is_ref": true},
                                {"cpic_value": a, "vcf_pos": p, "vcf_ref": r, "vcf_alt": a, "is_ref": false}]})
    };
    let allele = |name: &str, is_ref: bool, func: &str, vars: Value| {
        json!({"name": name, "is_reference": is_ref, "clinical_function": func, "activity_value": null,
               "activity_label": null, "variants": vars})
    };
    let mut genes = read_json("genes.json");
    genes["genes"].as_array_mut().unwrap().push(json!({
        "gene": "MINI2C19", "chrom": "chr10", "model": "star_allele", "ploidy": "diploid",
        "lookup_method": "PHENOTYPE", "reference_allele": "*38",
        "positions": [pos(0, 1000, "rs90", "A", "G"), pos(1, 1100, "rs91", "G", "A"), pos(2, 1200, "rs92", "A", "G")],
        "alleles": [
            allele("*38", true, "Normal function", json!([])),
            allele("*1", false, "Normal function", json!([{"pos_id": 2, "accepted": ["G"]}])),
            allele("*2", false, "No function", json!([{"pos_id": 0, "accepted": ["A", "G"]}, {"pos_id": 1, "accepted": ["A"]}, {"pos_id": 2, "accepted": ["A", "G"]}])),
            allele("*35", false, "No function", json!([{"pos_id": 0, "accepted": ["G"]}, {"pos_id": 2, "accepted": ["G"]}])),
        ],
    }));
    // MININAT2 / MININAT2B: NAT2-shaped. *1 is the GRCh38 haplotype (rare), *4 = rs95 G,
    // *5 = rs96 A, *16 = both. Het at both, unphased: *1/*16 scores 2 + 2 = 4, *4/*5 1 + 1 = 2,
    // but *4/*5 is far more common. MININAT2 is chosen by frequency, MININAT2B (same alleles,
    // no diplotype_choice) by score.
    let freqs = |eur: f64, eas: f64| json!({"European": eur, "East Asian": eas});
    let nat_allele = |name: &str, is_ref: bool, vars: Value, f: Value| {
        json!({"name": name, "is_reference": is_ref, "clinical_function": "Normal function", "activity_value": null,
               "activity_label": null, "pop_freqs": f, "variants": vars})
    };
    for (gname, chrom, choice) in [("MININAT2", "chr8", true), ("MININAT2B", "chr9", false)] {
        let npos = |id: usize, p: u64, rs: &str, r: &str, a: &str| {
            let mut v = pos(id, p, rs, r, a);
            v["chrom"] = json!(chrom);
            v
        };
        let mut g = json!({
            "gene": gname, "chrom": chrom, "model": "star_allele", "ploidy": "diploid",
            "lookup_method": "PHENOTYPE", "reference_allele": "*1",
            "positions": [npos(0, 2000, "rs95", "A", "G"), npos(1, 2100, "rs96", "G", "A")],
            "alleles": [
                nat_allele("*1", true, json!([]), freqs(0.005, 0.006)),
                nat_allele("*4", false, json!([{"pos_id": 0, "accepted": ["G"]}]), freqs(0.22, 0.51)),
                nat_allele("*5", false, json!([{"pos_id": 1, "accepted": ["A"]}]), freqs(0.42, 0.04)),
                nat_allele("*16", false, json!([{"pos_id": 0, "accepted": ["G"]}, {"pos_id": 1, "accepted": ["A"]}]), freqs(0.03, 0.0)),
            ],
        });
        if choice {
            g["diplotype_choice"] = json!("frequency");
        }
        genes["genes"].as_array_mut().unwrap().push(g);
    }
    let mut phen = read_json("phenotypes.json");
    for gname in ["MININAT2", "MININAT2B"] {
        phen["genes"][gname] = json!({"lookup_method": "PHENOTYPE", "rules": [
            {"counts": {"Normal function": 2}, "phenotype": "Rapid Metabolizer", "activity_score": null}]});
    }
    phen["genes"]["MINI2C19"] = json!({"lookup_method": "PHENOTYPE", "rules": [
        {"counts": {"Normal function": 2}, "phenotype": "Normal Metabolizer", "activity_score": null},
        {"counts": {"Normal function": 1, "No function": 1}, "phenotype": "Intermediate Metabolizer", "activity_score": null},
        {"counts": {"No function": 2}, "phenotype": "Poor Metabolizer", "activity_score": null}]});
    let mut guide = read_json("guidelines.json");
    for (ph, rec) in [("Intermediate Metabolizer", "Use minidrug at the standard dose."), ("Poor Metabolizer", "Avoid minidrug.")] {
        guide["recommendations"].as_array_mut().unwrap().push(json!({
            "drug": "minidrug", "guideline_id": "synth-g9", "guideline_name": "Synthetic CYP2C19-like Guideline",
            "genes": ["MINI2C19"], "untested_genes": [], "lookup": {"MINI2C19": ph}, "population": "general",
            "classification": "Strong", "recommendation": rec, "url": null, "source": "Synthetic CPIC-style Fixture"}));
    }
    let mut manifest = read_json("manifest.json");
    for (name, val) in [("genes.json", genes), ("phenotypes.json", phen), ("guidelines.json", guide)] {
        let bytes = serde_json::to_vec_pretty(&val).unwrap();
        std::fs::write(dst.join(name), &bytes).unwrap();
        manifest["files"][name]["sha256"] = json!(sha256_hex(&bytes));
        manifest["files"][name]["bytes"] = json!(bytes.len());
    }
    std::fs::write(dst.join("manifest.json"), serde_json::to_vec_pretty(&manifest).unwrap()).unwrap();
}

fn main() {
    // Read before the synthetic HLA tests point MMGEN_HLA_DIR at their own files (section 20).
    let real_hla_dir = std::env::var("MMGEN_HLA_DIR").ok();
    let mut h = Harness::new("suite");
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let repo_root = manifest_dir.parent().unwrap().to_path_buf();
    let fixture_dir = manifest_dir.join("fixtures").join("pgx_tables");
    let fixture = fixture_dir.to_string_lossy().to_string();

    // 1. FIPS 180-4 SHA-256 NIST vectors
    h.eq_str(
        &sha256_hex(b""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
        "sha256 empty",
    );
    h.eq_str(
        &sha256_hex(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
        "sha256 abc",
    );
    h.eq_str(
        &sha256_hex(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
        "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1",
        "sha256 448-bit",
    );

    // 2. Fixture bundle load + tampered checksum + zip error
    let b_ok = load_bundle(Some(&fixture));
    h.check(b_ok.is_ok(), "load synthetic fixture bundle");

    let fake_zip = h.p("tables.zip");
    std::fs::write(&fake_zip, b"PK\x03\x04").unwrap();
    h.expect_err(load_bundle(Some(&fake_zip)), &[".zip", "manifest.json"], "reject zip path");

    let tampered_dir = h.dir.join("tampered_tables");
    copy_dir(&fixture_dir, &tampered_dir);
    let mut g_bytes = std::fs::read(tampered_dir.join("genes.json")).unwrap();
    // Flip a whitespace/character without changing JSON validity or change size: mutate last newline -> space (same byte len!)
    let n = g_bytes.len();
    g_bytes[n - 1] = b' ';
    std::fs::write(tampered_dir.join("genes.json"), &g_bytes).unwrap();
    h.expect_err(
        load_bundle(Some(&tampered_dir.to_string_lossy())),
        &["SHA-256", "genes.json"],
        "detect tampered genes.json SHA-256",
    );

    // 3. Header/assembly & sample validation errors
    let grch37_vcf = h.p("grch37.vcf");
    write_text(
        &grch37_vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=249250621>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n",
    );
    h.expect_err(
        run(
            &grch37_vcf,
            &h.p("out_37"),
            &PgxOptions { tables: Some(fixture.clone()), ..Default::default() },
        ),
        &["GRCh37", "GRCh38"],
        "reject GRCh37 VCF",
    );

    let multi_vcf = h.p("multi.vcf");
    write_text(
        &multi_vcf,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n",
    );
    h.expect_err(
        run(
            &multi_vcf,
            &h.p("out_multi1"),
            &PgxOptions { tables: Some(fixture.clone()), ..Default::default() },
        ),
        &["2 samples", "sample="],
        "require sample= on multi-sample VCF",
    );
    h.expect_err(
        run(
            &multi_vcf,
            &h.p("out_multi2"),
            &PgxOptions {
                sample: Some("NO_SUCH".into()),
                tables: Some(fixture.clone()),
                ..Default::default()
            },
        ),
        &["NO_SUCH", "S1"],
        "error on unknown sample=",
    );

    // 4. Unphased ambiguity -> called: false vs phased trans/cis -> called: true
    let unphased_vcf = h.p("unphased.vcf");
    write_text(
        &unphased_vcf,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
         chr1\t100\trs1000\tA\tG\t50\tPASS\t.\tGT\t0/1\n\
         chr1\t120\trs1001\tC\tT\t50\tPASS\t.\tGT\t0/1\n\
         chr2\t200\trs2000\tC\tT\t50\tPASS\t.\tGT\t0/1\n\
         chr3\t300\trs3000\tG\tA\t50\tPASS\t.\tGT\t0/1\n\
         chr3\t320\trs3001\tT\tC\t50\tPASS\t.\tGT\t0/1\n\
         chrX\t400\trs4000\tG\tA\t50\tPASS\t.\tGT\t1\n\
         chrM\t500\trs5000\tA\tG\t50\tPASS\t.\tGT\t1\n\
         chrM\t520\t.\tC\t.\t50\tPASS\t.\tGT\t0\n",
    );
    let r_unph = run(
        &unphased_vcf,
        &h.p("out_unph"),
        &PgxOptions { tables: Some(fixture.clone()), ..Default::default() },
    )
    .unwrap();
    let v_unph: serde_json::Value = serde_json::from_str(&r_unph.json_text).unwrap();
    h.check(v_unph["engine"].as_str() == Some("mmgen-pgx"), "engine == mmgen-pgx");
    // SYNG_STAR: *1/*4 (the reference counts both positions read + *4 two) = 4 beats *2/*3 = 1 + 1.
    let star = gene_by_name(&v_unph, "SYNG_STAR");
    h.eq_str(star["diplotype"].as_str().unwrap_or(""), "*1/*4", "SYNG_STAR unphased: highest score *1/*4 is called");
    h.check(
        star["called"].as_bool() == Some(true)
            && star["x_score"].as_u64() == Some(4)
            && star["phenotype"].as_str() == Some("Intermediate Metabolizer")
            && star["activity_score"].as_f64() == Some(1.0),
        "SYNG_STAR *1/*4: score 4, Intermediate Metabolizer, AS 1.0",
    );
    let star_alts = star["x_alternatives"].as_array().cloned().unwrap_or_default();
    h.check(
        star_alts.len() == 1
            && star_alts[0]["diplotype"].as_str() == Some("*2/*3")
            && star_alts[0]["score"].as_u64() == Some(2)
            && star_alts[0]["phenotype"].as_str() == Some("Poor Metabolizer")
            && star_alts[0]["changes_phenotype"].as_bool() == Some(true)
            && star_alts[0]["changes_recommendation"].as_bool() == Some(true),
        "SYNG_STAR x_alternatives: *2/*3 (score 2, Poor) changes phenotype and recommendation",
    );
    let star_notes: Vec<&str> = star["notes"].as_array().unwrap().iter().filter_map(|n| n.as_str()).collect();
    h.check(
        star_notes.iter().any(|n| *n == "Also possible with the positions read: *2/*3 (Poor Metabolizer, activity score 0.5)"),
        "SYNG_STAR note names the alternative with its phenotype",
    );
    let unph_warnings: Vec<&str> = v_unph["warnings"].as_array().unwrap().iter().filter_map(|w| w.as_str()).collect();
    h.check(
        unph_warnings.iter().any(|w| w.starts_with("SYNG_STAR: *1/*4 (Intermediate Metabolizer, activity score 1.0) was chosen as the highest-scoring diplotype")
            && w.contains("*2/*3 (Poor Metabolizer, activity score 0.5) is also possible")
            && w.contains("would change the phenotype or a recommendation")),
        "warning when the alternative changes phenotype and recommendation",
    );
    // SYNG_ACT (variant model): cis gives AS 1.0, trans 0.5 -> trans, as CPIC/PharmCAT do for DPYD.
    let act = gene_by_name(&v_unph, "SYNG_ACT");
    h.check(
        act["called"].as_bool() == Some(true)
            && act["diplotype"].as_str() == Some("varA/varB")
            && act["activity_score"].as_f64() == Some(0.5)
            && act["phenotype"].as_str() == Some("Poor Metabolizer"),
        "SYNG_ACT unphased: trans (lowest activity score 0.5) is called",
    );
    let act_alts = act["x_alternatives"].as_array().cloned().unwrap_or_default();
    h.check(
        act_alts.len() == 1
            && act_alts[0]["diplotype"].as_str() == Some("Reference/varA+varB")
            && act_alts[0]["activity_score"].as_f64() == Some(1.0)
            && act_alts[0]["score"].is_null(),
        "SYNG_ACT alternative: cis Reference/varA+varB (AS 1.0), no score",
    );
    h.check(
        act["notes"].as_array().unwrap().iter().any(|n| n.as_str().is_some_and(|s| s.contains("in trans"))),
        "SYNG_ACT note says trans was assumed",
    );
    h.check(
        unph_warnings.iter().any(|w| w.starts_with("SYNG_ACT: varA/varB (Poor Metabolizer, activity score 0.5) was chosen assuming the variants are in trans")),
        "SYNG_ACT warning names the trans assumption",
    );
    let html_unph = std::fs::read_to_string(&r_unph.report_html_path).unwrap_or_default();
    h.check(
        html_unph.contains("<h2>Warnings</h2>") && html_unph.contains("Also possible with the positions read"),
        "report.html lists the warnings and the alternatives",
    );
    h.eq_str(
        gene_by_name(&v_unph, "SYNG_VG")["diplotype"].as_str().unwrap_or(""),
        "rs2000 reference (C)/rs2000 variant (T)",
        "SYNG_VG variant_genotype het",
    );
    h.eq_str(
        gene_by_name(&v_unph, "SYNG_X")["diplotype"].as_str().unwrap_or(""),
        "A-",
        "SYNG_X male haploid call",
    );
    h.eq_str(
        gene_by_name(&v_unph, "SYNG_MITO")["diplotype"].as_str().unwrap_or(""),
        "500A>G",
        "SYNG_MITO haploid call",
    );
    let drugs_unph: Vec<&str> = v_unph["drugs"].as_array().unwrap().iter().filter_map(|d| d["drug"].as_str()).collect();
    h.check(drugs_unph.contains(&"syndrug_x") && drugs_unph.contains(&"syndrug_mito"), "matched x + mito drugs");
    h.check(drugs_unph.contains(&"warfarin"), "warfarin given now that SYNG_ACT is called (trans assumed)");
    h.check(
        v_unph["drugs"].as_array().unwrap().iter().any(|d| d["drug"] == "syndrug_a" && d["recommendation"] == "Use reduced starting dose of syndrug_a."),
        "syndrug_a follows the chosen SYNG_STAR call (AS 1.0)",
    );

    // 5. Phased trans -> SYNG_STAR *2/*3 (AS 0.5), SYNG_ACT varA/varB (AS 0.5), triggers warfarin & syndrug_a
    // Also omit SYNG_MITO pos 1 (is_key: false) and keep callable_bed for chrX ref female diploid check.
    let phased_vcf = h.p("phased.vcf");
    let phased_body = "\
        ##fileformat=VCFv4.2\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
        chr1\t100\trs1000\tA\tG\t50\tPASS\t.\tGT:PS\t1|0:10\n\
        chr1\t120\trs1001\tC\tA\t50\tPASS\t.\tGT:PS\t0|1:10\n\
        chr2\t200\trs2000\tC\tT\t50\tPASS\t.\tGT\t0/0\n\
        chr3\t300\trs3000\tG\tA\t50\tPASS\t.\tGT:PS\t1|0:30\n\
        chr3\t320\trs3001\tT\tC\t50\tPASS\t.\tGT:PS\t0|1:30\n\
        chrX\t400\trs4000\tG\tA\t50\tPASS\t.\tGT\t0/1\n\
        chrM\t500\trs5000\tA\t.\t50\tPASS\t.\tGT\t0\n\
        chrM\t520\t.\tC\tCC\t50\tPASS\t.\tGT\t1\n";
    write_text(&phased_vcf, phased_body);
    let r_ph1 = run(
        &phased_vcf,
        &h.p("out_ph1"),
        &PgxOptions { tables: Some(fixture.clone()), ..Default::default() },
    )
    .unwrap();
    let v_ph1: serde_json::Value = serde_json::from_str(&r_ph1.json_text).unwrap();
    h.eq_str(
        gene_by_name(&v_ph1, "SYNG_STAR")["diplotype"].as_str().unwrap_or(""),
        "*2/*3",
        "SYNG_STAR phased trans *2/*3 (IUPAC A accepted by *3)",
    );
    h.check(
        (gene_by_name(&v_ph1, "SYNG_STAR")["activity_score"].as_f64().unwrap_or(-1.0) - 0.5).abs() < 1e-9,
        "SYNG_STAR phased trans AS 0.5",
    );
    h.eq_str(
        gene_by_name(&v_ph1, "SYNG_ACT")["diplotype"].as_str().unwrap_or(""),
        "varA/varB",
        "SYNG_ACT phased trans varA/varB",
    );
    h.eq_str(
        gene_by_name(&v_ph1, "SYNG_X")["diplotype"].as_str().unwrap_or(""),
        "B (reference)/A-",
        "SYNG_X female diploid B (reference)/A-",
    );
    h.eq_str(
        gene_by_name(&v_ph1, "SYNG_MITO")["diplotype"].as_str().unwrap_or(""),
        "520insC",
        "SYNG_MITO indel 520insC",
    );
    let drugs_ph1: Vec<&str> = v_ph1["drugs"].as_array().unwrap().iter().filter_map(|d| d["drug"].as_str()).collect();
    h.check(drugs_ph1.contains(&"warfarin") && drugs_ph1.contains(&"syndrug_a"), "phased calls trigger warfarin + syndrug_a");

    // 6. Plain .vcf vs BGZF .vcf.gz + .tbi + determinism check
    let phased_gz = h.p("phased.vcf.gz");
    write_text(&phased_gz, phased_body);
    indice::indexar(&phased_gz, None).unwrap();
    let r_gz1 = run(
        &phased_gz,
        &h.p("out_ph1"),
        &PgxOptions { tables: Some(fixture.clone()), ..Default::default() },
    )
    .unwrap();
    h.eq_str(&r_gz1.json_text, &r_ph1.json_text, "plain VCF vs .vcf.gz+.tbi byte-for-byte identity");
    let r_gz2 = run(
        &phased_gz,
        &h.p("out_ph1"),
        &PgxOptions { tables: Some(fixture.clone()), ..Default::default() },
    )
    .unwrap();
    h.eq_str(&r_gz2.json_text, &r_gz1.json_text, "repeated run determinism");

    // 7. Key position missing -> called: false; Non-key missing with 0/0 ref at key -> called: true via phantom check
    let key_gate_vcf = h.p("key_gate.vcf");
    write_text(
        &key_gate_vcf,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
         chr1\t100\trs1000\tA\t.\t50\tPASS\t.\tGT\t0/0\n\
         chr2\t200\trs2000\tC\t.\t50\tPASS\t.\tGT\t./.\n",
    );
    let r_kg = run(
        &key_gate_vcf,
        &h.p("out_kg"),
        &PgxOptions { tables: Some(fixture.clone()), ..Default::default() },
    )
    .unwrap();
    let v_kg: serde_json::Value = serde_json::from_str(&r_kg.json_text).unwrap();
    h.eq_str(
        gene_by_name(&v_kg, "SYNG_STAR")["diplotype"].as_str().unwrap_or(""),
        "*1/*1",
        "SYNG_STAR missing non-key pos calls *1/*1 (phantom eliminated)",
    );
    h.check(
        !gene_by_name(&v_kg, "SYNG_VG")["called"].as_bool().unwrap(),
        "SYNG_VG missing key pos not called",
    );
    h.check(
        !gene_by_name(&v_kg, "SYNG_ACT")["called"].as_bool().unwrap(),
        "SYNG_ACT absent from sparse VCF + no callable_bed not called",
    );

    // 7b. A key position not read no longer stops the call: note, warning, x_missing_key_positions.
    //     Here *1/*3 and *1/*4 tie (score 2) with different activity scores -> no recommendation.
    let key_miss_vcf = h.p("key_miss.vcf");
    write_text(
        &key_miss_vcf,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
         chr1\t100\trs1000\tA\tG\t50\tPASS\t.\tGT\t./.\n\
         chr1\t120\trs1001\tC\tT\t50\tPASS\t.\tGT\t0/1\n",
    );
    let r_km = run(&key_miss_vcf, &h.p("out_km"), &PgxOptions { tables: Some(fixture.clone()), ..Default::default() }).unwrap();
    let v_km: serde_json::Value = serde_json::from_str(&r_km.json_text).unwrap();
    let km = gene_by_name(&v_km, "SYNG_STAR");
    h.check(
        km["called"].as_bool() == Some(true)
            && km["diplotype"].as_str() == Some("*1/*3 OR *1/*4")
            && km["phenotype"].as_str() == Some("Normal Metabolizer OR Intermediate Metabolizer")
            && km["activity_score"].is_null()
            && km["x_score"].as_u64() == Some(2),
        "key position not read: tie *1/*3 OR *1/*4 called, both phenotypes, no activity score",
    );
    h.check(
        km["x_missing_key_positions"].as_array().map(|a| a.len() == 1 && a[0] == "rs1000") == Some(true),
        "x_missing_key_positions lists rs1000",
    );
    let km_warn: Vec<&str> = v_km["warnings"].as_array().unwrap().iter().filter_map(|w| w.as_str()).collect();
    h.check(
        km_warn.iter().any(|w| w.starts_with("SYNG_STAR: key position(s) not read (rs1000 (chr1:100))") && w.contains("callable_bed="))
            && km_warn.iter().any(|w| w.starts_with("SYNG_STAR: tied diplotypes give different phenotypes")),
        "warnings: the key position not read and the tie",
    );
    h.check(
        !v_km["drugs"].as_array().unwrap().iter().any(|d| d["drug"] == "syndrug_a"),
        "no syndrug_a recommendation while the tie is unresolved",
    );

    // 7c. CYP2C19-like gene: *38 is the reference and *1 differs from it at one position.
    //     *1 is scored as the reference (every position read), so *1/*2 (3 + 3) beats *2/*35
    //     (3 + 2), as in NA12878; scoring *1 by its one position would pick *2/*35 (1 + 3 < 5).
    let mini_dir = h.dir.join("mini_tables");
    add_mini_cyp2c19(&fixture_dir, &mini_dir);
    let mini = mini_dir.to_string_lossy().to_string();
    h.check(load_bundle(Some(&mini)).is_ok(), "mini CYP2C19 bundle loads (checksums rewritten)");
    let mini_unph = h.p("mini_unph.vcf");
    write_text(
        &mini_unph,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
         chr10\t1000\trs90\tA\tG\t50\tPASS\t.\tGT\t0/1\n\
         chr10\t1100\trs91\tG\tA\t50\tPASS\t.\tGT\t0/1\n\
         chr10\t1200\trs92\tA\tG\t50\tPASS\t.\tGT\t1/1\n",
    );
    let v_mu: serde_json::Value =
        serde_json::from_str(&run(&mini_unph, &h.p("out_mini_u"), &PgxOptions { tables: Some(mini.clone()), ..Default::default() }).unwrap().json_text).unwrap();
    let mu = gene_by_name(&v_mu, "MINI2C19");
    h.check(
        mu["diplotype"].as_str() == Some("*1/*2") && mu["x_score"].as_u64() == Some(6) && mu["phenotype"].as_str() == Some("Intermediate Metabolizer"),
        "mini CYP2C19 unphased: *1/*2 (score 6) called, Intermediate Metabolizer",
    );
    let mu_alts = mu["x_alternatives"].as_array().cloned().unwrap_or_default();
    h.check(
        mu_alts.len() == 1
            && mu_alts[0]["diplotype"].as_str() == Some("*2/*35")
            && mu_alts[0]["score"].as_u64() == Some(5)
            && mu_alts[0]["changes_phenotype"].as_bool() == Some(true)
            && mu_alts[0]["changes_recommendation"].as_bool() == Some(true),
        "mini CYP2C19: *2/*35 (score 5, Poor) is the alternative and would change the recommendation",
    );
    h.check(
        v_mu["drugs"].as_array().unwrap().iter().any(|d| d["drug"] == "minidrug" && d["recommendation"] == "Use minidrug at the standard dose."),
        "mini CYP2C19: recommendation follows *1/*2",
    );
    // Phase decides when it is there: the same genotypes phased as *2/*35 give *2/*35 alone.
    let mini_ph = h.p("mini_ph.vcf");
    write_text(
        &mini_ph,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
         chr10\t1000\trs90\tA\tG\t50\tPASS\t.\tGT:PS\t0|1:7\n\
         chr10\t1100\trs91\tG\tA\t50\tPASS\t.\tGT:PS\t1|0:7\n\
         chr10\t1200\trs92\tA\tG\t50\tPASS\t.\tGT\t1/1\n",
    );
    let v_mp: serde_json::Value =
        serde_json::from_str(&run(&mini_ph, &h.p("out_mini_p"), &PgxOptions { tables: Some(mini.clone()), ..Default::default() }).unwrap().json_text).unwrap();
    let mp = gene_by_name(&v_mp, "MINI2C19");
    h.check(
        mp["diplotype"].as_str() == Some("*2/*35") && mp["phenotype"].as_str() == Some("Poor Metabolizer") && mp.get("x_alternatives").is_none(),
        "mini CYP2C19 phased: *2/*35 alone (phase beats the score)",
    );

    // 7d. NAT2-shaped genes: chosen by population frequency (as PharmCAT does for NAT2) or,
    //     without diplotype_choice, by score. *4/*5: mean of 2·0.22·0.42 and 2·0.51·0.04.
    let nat_vcf = h.p("mini_nat2.vcf");
    write_text(
        &nat_vcf,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
         chr8\t2000\trs95\tA\tG\t50\tPASS\t.\tGT\t0/1\n\
         chr8\t2100\trs96\tG\tA\t50\tPASS\t.\tGT\t0/1\n\
         chr9\t2000\trs95\tA\tG\t50\tPASS\t.\tGT\t0/1\n\
         chr9\t2100\trs96\tG\tA\t50\tPASS\t.\tGT\t0/1\n",
    );
    let v_nat: serde_json::Value =
        serde_json::from_str(&run(&nat_vcf, &h.p("out_mini_nat"), &PgxOptions { tables: Some(mini.clone()), ..Default::default() }).unwrap().json_text).unwrap();
    let nf = gene_by_name(&v_nat, "MININAT2");
    let want_freq = (2.0 * 0.22 * 0.42 + 2.0 * 0.51 * 0.04) / 2.0;
    h.check(
        nf["diplotype"].as_str() == Some("*4/*5")
            && nf["x_chosen_by"].as_str() == Some("frequency")
            && nf["x_frequency"].as_f64().is_some_and(|f| (f - want_freq).abs() < 1e-12)
            && nf["x_score"].as_u64() == Some(2),
        "NAT2-shaped, frequency rule: *4/*5 (the common one) is called over *1/*16 (higher score)",
    );
    let nf_alts = nf["x_alternatives"].as_array().cloned().unwrap_or_default();
    h.check(
        nf_alts.len() == 1
            && nf_alts[0]["diplotype"].as_str() == Some("*1/*16")
            && nf_alts[0]["score"].as_u64() == Some(4)
            && nf_alts[0]["frequency"].as_f64().is_some_and(|f| f > 0.0 && f < 0.001),
        "NAT2-shaped: *1/*16 kept as the alternative, with its score and frequency",
    );
    h.check(
        nf["notes"].as_array().unwrap().iter().any(|n| n.as_str().is_some_and(|s| {
            s.starts_with("Chosen by population frequency, as PharmCAT does for MININAT2: *4/*5")
                && s.ends_with("By score alone, *1/*16 (score 4) would be chosen.")
        })),
        "NAT2-shaped: the note says the call is by frequency and what the score would pick",
    );
    let ns = gene_by_name(&v_nat, "MININAT2B");
    h.check(
        ns["diplotype"].as_str() == Some("*1/*16") && ns["x_chosen_by"].as_str() == Some("score") && ns.get("x_frequency").is_none(),
        "same gene without diplotype_choice: the score picks *1/*16",
    );

    // 8. Callable BED fills homozygous reference on absent VCF positions
    let cal_bed = h.p("callable.bed");
    write_text(
        &cal_bed,
        "chr1\t90\t150\n\
         chr2\t190\t220\n\
         chr3\t290\t360\n\
         chrX\t390\t420\n\
         chrM\t490\t540\n",
    );
    let sparse_vcf = h.p("sparse.vcf");
    write_text(
        &sparse_vcf,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
         chr1\t100\trs1000\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
    );
    let r_cal = run(
        &sparse_vcf,
        &h.p("out_cal"),
        &PgxOptions {
            callable_bed: Some(cal_bed),
            tables: Some(fixture.clone()),
            ..Default::default()
        },
    )
    .unwrap();
    let v_cal: serde_json::Value = serde_json::from_str(&r_cal.json_text).unwrap();
    h.eq_str(gene_by_name(&v_cal, "SYNG_STAR")["diplotype"].as_str().unwrap_or(""), "*1/*2", "callable_bed + sparse variant SYNG_STAR *1/*2");
    h.eq_str(gene_by_name(&v_cal, "SYNG_ACT")["diplotype"].as_str().unwrap_or(""), "Reference/Reference", "callable_bed fills SYNG_ACT Reference/Reference");
    h.check(Path::new(&r_cal.report_html_path).exists(), "report.html written");

    // 9. Optional FASTA reference= enables left-alignment of right-shifted repeat indel
    let fa_path = h.p("ref38.fa");
    // chr3 1-based positions 348..355 = TTCAAAAC (350=C, 351..354=AAAA, 355=C)
    let mut seq3 = vec![b'N'; 347];
    seq3.extend_from_slice(b"TTCAAAAC");
    write_text(&fa_path, &format!(">chr3\n{}\n", String::from_utf8(seq3).unwrap()));

    let shifted_vcf = h.p("shifted.vcf");
    write_text(
        &shifted_vcf,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
         chr3\t300\trs3000\tG\t.\t50\tPASS\t.\tGT\t0/0\n\
         chr3\t320\trs3001\tT\t.\t50\tPASS\t.\tGT\t0/0\n\
         chr3\t353\t.\tAA\tA\t50\tPASS\t.\tGT\t0/1\n",
    );
    let r_nofa = run(
        &shifted_vcf,
        &h.p("out_nofa"),
        &PgxOptions { tables: Some(fixture.clone()), ..Default::default() },
    )
    .unwrap();
    let v_nofa: serde_json::Value = serde_json::from_str(&r_nofa.json_text).unwrap();
    h.eq_str(
        gene_by_name(&v_nofa, "SYNG_ACT")["diplotype"].as_str().unwrap_or(""),
        "Reference/Reference",
        "without FASTA, right-shifted 353 AA->A does not hit 350",
    );
    let r_withfa = run(
        &shifted_vcf,
        &h.p("out_withfa"),
        &PgxOptions {
            reference: Some(fa_path),
            tables: Some(fixture.clone()),
            ..Default::default()
        },
    )
    .unwrap();
    let v_withfa: serde_json::Value = serde_json::from_str(&r_withfa.json_text).unwrap();
    h.eq_str(
        gene_by_name(&v_withfa, "SYNG_ACT")["diplotype"].as_str().unwrap_or(""),
        "Reference/varDel",
        "with FASTA, right-shifted 353 AA->A left-aligns to 350 CA->C (varDel)",
    );

    // 10. Traits evaluation: no silent calls without BED, callable_bed REF fill, and VCF alleles
    let trait_missing_res = run(
        &phased_vcf,
        &h.p("out_trait_miss"),
        &PgxOptions { tables: Some(fixture.clone()), ..Default::default() },
    )
    .unwrap();
    let v_tmiss: serde_json::Value = serde_json::from_str(&trait_missing_res.json_text).unwrap();
    let traits_miss = v_tmiss["traits"].as_array().unwrap();
    h.check(traits_miss.len() == 2, "fixture returns 2 traits");
    h.check(
        traits_miss.iter().all(|tr| !tr["x_called"].as_bool().unwrap() && tr["interpretation"].as_str().unwrap().starts_with("Not called")),
        "traits without VCF records and without callable_bed are never silently called",
    );

    let traits_vcf = h.p("traits.vcf");
    write_text(
        &traits_vcf,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
         chr2\t135851076\trs4988235\tG\tA\t50\tPASS\t.\tGT\t0/1\n\
         chr2\t135859184\trs182549\tC\tT\t50\tPASS\t.\tGT\t0/1\n\
         chr15\t28120472\trs12913832\tA\tG\t50\tPASS\t.\tGT\t0/1\n\
         chr15\t27985172\trs1800407\tC\tT\t50\tPASS\t.\tGT\t0/1\n\
         chr14\t92307319\trs12896399\tG\tT\t50\tPASS\t.\tGT\t0/0\n\
         chr5\t33951588\trs16891982\tC\tG\t50\tPASS\t.\tGT\t1/1\n\
         chr11\t89277878\trs1393350\tG\tA\t50\tPASS\t.\tGT\t0/0\n\
         chr6\t396321\trs12203592\tC\tT\t50\tPASS\t.\tGT\t0/1\n",
    );
    let r_tr = run(
        &traits_vcf,
        &h.p("out_traits"),
        &PgxOptions { tables: Some(fixture.clone()), ..Default::default() },
    )
    .unwrap();
    let v_tr: serde_json::Value = serde_json::from_str(&r_tr.json_text).unwrap();
    let by_tid: std::collections::HashMap<String, &serde_json::Value> = v_tr["traits"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| (x["x_id"].as_str().unwrap().to_string(), x))
        .collect();
    h.check(
        by_tid["lactase_persistence"]["x_called"].as_bool() == Some(true)
            && by_tid["lactase_persistence"]["interpretation"].as_str().unwrap().contains("Lactase persistence likely"),
        "lactase rs4988235 het calls Lactase persistence likely",
    );
    h.check(
        by_tid["eye_color"]["x_called"].as_bool() == Some(true)
            && by_tid["eye_color"]["x_positions_read"].as_u64() == Some(6)
            && by_tid["eye_color"]["interpretation"].as_str().unwrap().contains("intermediate (green / hazel)"),
        "eye color A/G + OCA2*T + 4 modifiers calls intermediate (green / hazel)",
    );

    let traits_bed = h.p("traits_ref.bed");
    write_text(
        &traits_bed,
        "chr2\t135851000\t135859300\n\
         chr5\t33951500\t33951700\n\
         chr6\t396300\t396400\n\
         chr11\t89277800\t89277950\n\
         chr14\t92307250\t92307400\n\
         chr15\t27985100\t28120550\n",
    );
    let r_tr_bed = run(
        &phased_vcf,
        &h.p("out_tr_bed"),
        &PgxOptions { callable_bed: Some(traits_bed), tables: Some(fixture.clone()), ..Default::default() },
    )
    .unwrap();
    let v_trb: serde_json::Value = serde_json::from_str(&r_tr_bed.json_text).unwrap();
    let by_tid_b: std::collections::HashMap<String, &serde_json::Value> = v_trb["traits"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| (x["x_id"].as_str().unwrap().to_string(), x))
        .collect();
    h.check(
        by_tid_b["lactase_persistence"]["interpretation"].as_str().unwrap().contains("Lactase non-persistence")
            && by_tid_b["eye_color"]["interpretation"].as_str().unwrap().contains("brown / dark"),
        "callable_bed REF fill calls Lactase non-persistence + brown/dark eye color",
    );

    // 11. Real CPIC v1.61.0 bundle check (if installed under mmgen_tools/pgx/cpic-v1.61.0)
    let real_bundle_dir = repo_root.join("mmgen_tools").join("pgx").join("cpic-v1.61.0");
    if real_bundle_dir.join("manifest.json").exists() {
        let b = load_bundle(Some(&real_bundle_dir.to_string_lossy())).unwrap();
        h.check(b.genes.genes.len() == 19, "real cpic-v1.61.0 bundle has 19 genes");
        h.check(!b.guidelines.recommendations.is_empty(), "real cpic-v1.61.0 bundle has guidelines");
        h.check(b.traits.traits.len() == 2, "real cpic-v1.61.0 bundle has 2 traits");
    } else {
        println!("NOTE  skipping real bundle check ({} not present)", real_bundle_dir.display());
    }

    run_v2_bam_tests(&mut h, &fixture, &repo_root, real_hla_dir);

    println!("farmacogenetica tests: {} passed, {} failed", h.ok, h.fail);
    if h.fail > 0 {
        std::process::exit(1);
    }
}

fn write_sorted_bam(path: &str, lines: &[String]) {
    let mut w = bam::BamEscritor::novo(std::fs::File::create(path).unwrap());
    for l in lines {
        w.write_all(l.as_bytes()).unwrap();
        w.write_all(b"\n").unwrap();
    }
    w.termina().unwrap();
    indice::indexar(path, None).unwrap();
}

fn synth_dna_seq(seed: u64, len: usize) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut x = seed ^ 0x9E37_79B9_7F4A_7C15;
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        out.push(bases[((x >> 32) & 3) as usize]);
    }
    String::from_utf8(out).unwrap()
}

fn sam_line(qname: &str, flag: u16, rname: &str, pos: u64, seq: &str) -> String {
    let cigar = format!("{}M", seq.len());
    let qual = "I".repeat(seq.len());
    format!("{qname}\t{flag}\t{rname}\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}")
}

fn write_synth_hla_fastas(dir: &Path) -> [Vec<String>; 3] {
    std::fs::create_dir_all(dir).unwrap();
    let specs: [(&str, [(&str, u64); 3]); 3] = [
        (
            "A_gen.fasta",
            [
                ("A*01:01:01:01", 101),
                ("A*31:01:02:01", 202),
                ("A*02:01:01:01", 303),
            ],
        ),
        (
            "B_gen.fasta",
            [
                ("B*07:02:01:01", 404),
                ("B*57:01:01:01", 505),
                ("B*35:01:01:01", 606),
            ],
        ),
        (
            "C_gen.fasta",
            [
                ("C*04:01:01:01", 707),
                ("C*12:03:01:01", 808),
                ("C*07:01:01:01", 909),
            ],
        ),
    ];
    let mut seqs: [Vec<String>; 3] = Default::default();
    for (idx, (fname, alleles)) in specs.iter().enumerate() {
        let mut fa = String::new();
        for (i, (al_name, seed)) in alleles.iter().enumerate() {
            let s = synth_dna_seq(*seed, 240);
            fa.push_str(&format!(">HLA:HLA00{idx}{i} {al_name} 240 bp\n{s}\n"));
            seqs[idx].push(s);
        }
        write_text(&dir.join(fname).to_string_lossy(), &fa);
    }
    seqs
}

struct SynthBamOpts<'a> {
    ex9_h1: usize,
    ex9_h2: usize,
    body_h1: usize,
    body_h2: usize,
    h1_bases: [u8; 7],
    h2_bases: [u8; 7],
    hla_seqs: Option<[(&'a str, &'a str); 3]>,
}

fn build_v2_synth_bam(path: &str, opts: &SynthBamOpts<'_>) {
    let mut records: Vec<(&str, u64, String)> = Vec::new();
    let mut lines = vec![
        "@HD\tVN:1.6\tSO:coordinate".to_string(),
        "@SQ\tSN:chr4\tLN:190214555".to_string(),
        "@SQ\tSN:chr6\tLN:170805979".to_string(),
        "@SQ\tSN:chr22\tLN:50818468".to_string(),
    ];

    // chr4 ALB diploid control: 73404287..73421482 (17,196 bp) -> 20x diploid baseline
    let ctrl_seq = "A".repeat(1000);
    let mut pos = 73404287u64;
    let mut cidx = 0usize;
    while pos <= 73420483 {
        let len = (73421482 - pos + 1).min(1000) as usize;
        let s = &ctrl_seq[..len];
        for copy in 0..20 {
            records.push(("chr4", pos, sam_line(&format!("ctrl_{cidx}_{copy}"), 0, "chr4", pos, s)));
        }
        pos += 1000;
        cidx += 1;
    }

    // Optional HLA paired reads (4 pairs per selected allele)
    if let Some(hla) = opts.hla_seqs {
        let loci = [
            ("hlaA", 29936100u64, hla[0].0, hla[0].1),
            ("hlaC", 31265100u64, hla[2].0, hla[2].1),
            ("hlaB", 31350100u64, hla[1].0, hla[1].1),
        ];
        for (tag, base_pos, seq_a, seq_b) in loci {
            for (hi, seq) in [(1usize, seq_a), (2usize, seq_b)] {
                for p in 0..4usize {
                    let qname = format!("{tag}_{hi}_{p}");
                    let r1 = &seq[0..100];
                    let r2 = &seq[140..240];
                    records.push(("chr6", base_pos + p as u64, sam_line(&qname, 0x41, "chr6", base_pos + p as u64, r1)));
                    records.push(("chr6", base_pos + 140 + p as u64, sam_line(&qname, 0x81, "chr6", base_pos + 140 + p as u64, r2)));
                }
            }
        }
    }

    // chr22 CYP2D6 reads: 10x depth per copy -> exact CN estimates over ex9, mid, and 5p
    let make_hap_reads = |records: &mut Vec<(&'static str, u64, String)>,
                          tag: &str,
                          ex9_copies: usize,
                          body_copies: usize,
                          hb: [u8; 7]| {
        for c in 0..ex9_copies {
            for r in 0..10 {
                let mut ex9 = vec![b'A'; 301];
                ex9[42126611 - 42126500] = hb[0];
                let qn = format!("{tag}_e_{c}_{r}");
                records.push(("chr22", 42126500, sam_line(&qn, 0x41, "chr22", 42126500, &String::from_utf8(ex9).unwrap())));
            }
        }
        for c in 0..body_copies {
            for r in 0..10 {
                let qn = format!("{tag}_b_{c}_{r}");
                let mut mid = vec![b'A'; 1199];
                mid[42127941 - 42126801] = hb[1];
                records.push(("chr22", 42126801, sam_line(&qn, 0x41, "chr22", 42126801, &String::from_utf8(mid).unwrap())));

                let mut fivep = vec![b'A'; 2811];
                fivep[42128945 - 42128000] = hb[2];
                fivep[42129130 - 42128000] = hb[3];
                fivep[42129809 - 42128000] = hb[4];
                fivep[42129819 - 42128000] = hb[5];
                fivep[42130692 - 42128000] = hb[6];
                records.push(("chr22", 42128000, sam_line(&qn, 0x81, "chr22", 42128000, &String::from_utf8(fivep).unwrap())));
            }
        }
    };
    make_hap_reads(&mut records, "h1", opts.ex9_h1, opts.body_h1, opts.h1_bases);
    make_hap_reads(&mut records, "h2", opts.ex9_h2, opts.body_h2, opts.h2_bases);

    let contig_rank = |c: &str| match c {
        "chr4" => 0u8,
        "chr6" => 1u8,
        _ => 2u8,
    };
    records.sort_by(|a, b| contig_rank(a.0).cmp(&contig_rank(b.0)).then_with(|| a.1.cmp(&b.1)));
    for (_, _, rec) in records {
        lines.push(rec);
    }
    write_sorted_bam(path, &lines);
}

struct EnvGuard {
    saved_dir: Option<std::ffi::OsString>,
    saved_off: Option<std::ffi::OsString>,
}
impl EnvGuard {
    fn new(hla_dir: Option<&Path>, offline: bool) -> Self {
        let g = Self {
            saved_dir: std::env::var_os("MMGEN_HLA_DIR"),
            saved_off: std::env::var_os("MMGEN_HLA_OFFLINE"),
        };
        if let Some(d) = hla_dir {
            std::env::set_var("MMGEN_HLA_DIR", d);
        } else {
            std::env::remove_var("MMGEN_HLA_DIR");
        }
        if offline {
            std::env::set_var("MMGEN_HLA_OFFLINE", "1");
        } else {
            std::env::remove_var("MMGEN_HLA_OFFLINE");
        }
        g
    }
}
impl Drop for EnvGuard {
    fn drop(&mut self) {
        match &self.saved_dir {
            Some(v) => std::env::set_var("MMGEN_HLA_DIR", v),
            None => std::env::remove_var("MMGEN_HLA_DIR"),
        }
        match &self.saved_off {
            Some(v) => std::env::set_var("MMGEN_HLA_OFFLINE", v),
            None => std::env::remove_var("MMGEN_HLA_OFFLINE"),
        }
    }
}

fn run_v2_bam_tests(h: &mut Harness, fixture: &str, repo_root: &Path, real_hla_dir: Option<String>) {
    let fixture_dir = PathBuf::from(fixture);
    let hla_dir = h.dir.join("hla_synth");
    let hla_seqs = write_synth_hla_fastas(&hla_dir);
    let _env = EnvGuard::new(Some(&hla_dir), true);

    // Base minimal VCF on chr1 so sample name exists and BAM-omitted returns 19 genes + warning
    let base_vcf = h.p("v2_base.vcf");
    write_text(
        &base_vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr22,length=50818468>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA_V2\n",
    );

    // 14. BAM omitted keeps 19 genes, skips outside_calls.tsv, and emits warning
    let r_nobam = run(
        &base_vcf,
        &h.p("out_v2_nobam"),
        &PgxOptions { tables: Some(fixture.to_string()), ..Default::default() },
    )
    .unwrap();
    let v_nobam: serde_json::Value = serde_json::from_str(&r_nobam.json_text).unwrap();
    h.check(
        v_nobam["genes"].as_array().map(|a| a.len()) == Some(5),
        "BAM omitted preserves VCF-only gene count in fixture",
    );
    h.check(
        v_nobam["x_files"]["pharmcat_outside_calls"].is_null(),
        "BAM omitted skips pharmcat_outside_calls",
    );
    h.check(
        v_nobam["warnings"].as_array().is_some_and(|w| !w.is_empty()),
        "BAM omitted adds BAM skip warning",
    );

    // Star alleles base arrays: [42126611, 42127941, 42128945, 42129130, 42129809, 42129819, 42130692]
    let star1: [u8; 7] = [b'C', b'G', b'C', b'C', b'T', b'G', b'G'];
    let star2: [u8; 7] = [b'G', b'A', b'C', b'C', b'T', b'G', b'G'];
    let star4: [u8; 7] = [b'C', b'G', b'T', b'C', b'T', b'G', b'A'];
    let star10: [u8; 7] = [b'G', b'G', b'C', b'G', b'T', b'G', b'A'];

    // 15. CYP2D6 *2/*4 CN=2 (VCF unphased + BAM pileup + BAM fragment phasing) + HLA alerts + outside_calls.tsv
    let vcf_2_4 = h.p("v2_2_4.vcf");
    write_text(
        &vcf_2_4,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr22,length=50818468>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA_V2\n\
         chr22\t42126611\t.\tC\tG\t50\tPASS\t.\tGT\t0/1\n\
         chr22\t42127941\t.\tG\tA\t50\tPASS\t.\tGT\t0/1\n\
         chr22\t42128945\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n\
         chr22\t42130692\t.\tG\tA\t50\tPASS\t.\tGT\t0/1\n",
    );
    let bam_2_4 = h.p("v2_2_4.bam");
    build_v2_synth_bam(
        &bam_2_4,
        &SynthBamOpts {
            ex9_h1: 1,
            ex9_h2: 1,
            body_h1: 1,
            body_h2: 1,
            h1_bases: star2,
            h2_bases: star4,
            hla_seqs: Some([
                (&hla_seqs[0][0], &hla_seqs[0][1]), // A*01:01 / A*31:01
                (&hla_seqs[1][0], &hla_seqs[1][1]), // B*07:02 / B*57:01
                (&hla_seqs[2][0], &hla_seqs[2][1]), // C*04:01 / C*12:03
            ]),
        },
    );
    let r_2_4 = run(
        &vcf_2_4,
        &h.p("out_v2_2_4"),
        &PgxOptions { bam: Some(bam_2_4.clone()), tables: Some(fixture.to_string()), ..Default::default() },
    )
    .unwrap();
    let v_2_4: serde_json::Value = serde_json::from_str(&r_2_4.json_text).unwrap();
    let g_cyp = gene_by_name(&v_2_4, "CYP2D6");
    h.eq_str(g_cyp["diplotype"].as_str().unwrap_or(""), "*2/*4", "CYP2D6 *2/*4 CN=2 diplotype");
    h.eq_str(g_cyp["source"].as_str().unwrap_or(""), "vcf+bam", "CYP2D6 source is vcf+bam");
    h.check(g_cyp["copy_number"].as_f64() == Some(2.0), "CYP2D6 CN=2");
    h.check(g_cyp["x_phased"].as_bool() == Some(true), "CYP2D6 BAM phasing set x_phased=true");
    h.eq_str(gene_by_name(&v_2_4, "HLA-A")["diplotype"].as_str().unwrap_or(""), "*01:01/*31:01", "HLA-A *01:01/*31:01");
    h.eq_str(gene_by_name(&v_2_4, "HLA-B")["diplotype"].as_str().unwrap_or(""), "*07:02/*57:01", "HLA-B *07:02/*57:01");
    h.eq_str(gene_by_name(&v_2_4, "HLA-C")["diplotype"].as_str().unwrap_or(""), "*04:01/*12:03", "HLA-C *04:01/*12:03");
    h.eq_str(gene_by_name(&v_2_4, "HLA-C")["source"].as_str().unwrap_or(""), "bam", "HLA-C source is bam");

    let drugs_arr = v_2_4["drugs"].as_array().unwrap();
    h.check(drugs_arr.iter().any(|d| d["drug"] == "abacavir"), "HLA-B*57:01 triggers abacavir guideline");
    h.check(drugs_arr.iter().any(|d| d["drug"] == "carbamazepine"), "HLA-A*31:01 + HLA-B*15:02 negative triggers carbamazepine guideline");
    let oc_path = v_2_4["x_files"]["pharmcat_outside_calls"].as_str().unwrap_or("");
    h.check(!oc_path.is_empty() && Path::new(oc_path).exists(), "outside_calls.tsv created and in x_files");
    let oc_txt = std::fs::read_to_string(oc_path).unwrap_or_default();
    h.check(
        oc_txt.contains("CYP2D6\t*2/*4") && oc_txt.contains("HLA-A\t*01:01/*31:01") && oc_txt.contains("HLA-B\t*07:02/*57:01") && !oc_txt.contains("HLA-C"),
        "outside_calls.tsv exports CYP2D6/HLA-A/HLA-B and omits HLA-C",
    );

    // 16. CYP2D6 *5/*5 CN=0 and *1/*5 CN=1 deletions
    let hla_neg = Some([
        (hla_seqs[0][0].as_str(), hla_seqs[0][2].as_str()),
        (hla_seqs[1][0].as_str(), hla_seqs[1][2].as_str()),
        (hla_seqs[2][0].as_str(), hla_seqs[2][2].as_str()),
    ]);
    let bam_0 = h.p("v2_5_5.bam");
    build_v2_synth_bam(
        &bam_0,
        &SynthBamOpts { ex9_h1: 0, ex9_h2: 0, body_h1: 0, body_h2: 0, h1_bases: star1, h2_bases: star1, hla_seqs: hla_neg },
    );
    let r_0 = run(&base_vcf, &h.p("out_v2_0"), &PgxOptions { bam: Some(bam_0), tables: Some(fixture.to_string()), ..Default::default() }).unwrap();
    let v_0: serde_json::Value = serde_json::from_str(&r_0.json_text).unwrap();
    h.eq_str(gene_by_name(&v_0, "CYP2D6")["diplotype"].as_str().unwrap_or(""), "*5/*5", "CYP2D6 *5/*5 deletion");
    h.check(gene_by_name(&v_0, "CYP2D6")["copy_number"].as_f64() == Some(0.0), "CYP2D6 *5/*5 CN=0");
    h.check(v_0["drugs"].as_array().unwrap().iter().any(|d| d["drug"] == "codeine"), "*5/*5 triggers codeine guideline");

    let bam_1 = h.p("v2_1_5.bam");
    build_v2_synth_bam(
        &bam_1,
        &SynthBamOpts { ex9_h1: 1, ex9_h2: 0, body_h1: 1, body_h2: 0, h1_bases: star1, h2_bases: star1, hla_seqs: hla_neg },
    );
    let r_1 = run(&base_vcf, &h.p("out_v2_1"), &PgxOptions { bam: Some(bam_1), tables: Some(fixture.to_string()), ..Default::default() }).unwrap();
    let v_1: serde_json::Value = serde_json::from_str(&r_1.json_text).unwrap();
    h.eq_str(gene_by_name(&v_1, "CYP2D6")["diplotype"].as_str().unwrap_or(""), "*1/*5", "CYP2D6 *1/*5 hemizygous deletion");
    h.check(gene_by_name(&v_1, "CYP2D6")["copy_number"].as_f64() == Some(1.0), "CYP2D6 *1/*5 CN=1");

    // 17. CYP2D6 *1x2/*4 CN=3 duplication via BAM VAF skew (2:1)
    let bam_3 = h.p("v2_dup3.bam");
    build_v2_synth_bam(
        &bam_3,
        &SynthBamOpts { ex9_h1: 2, ex9_h2: 1, body_h1: 2, body_h2: 1, h1_bases: star1, h2_bases: star4, hla_seqs: hla_neg },
    );
    let vcf_1_4 = h.p("v2_1_4.vcf");
    write_text(
        &vcf_1_4,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr22,length=50818468>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA_V2\n\
         chr22\t42128945\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n\
         chr22\t42130692\t.\tG\tA\t50\tPASS\t.\tGT\t0/1\n",
    );
    let r_3 = run(&vcf_1_4, &h.p("out_v2_3"), &PgxOptions { bam: Some(bam_3), tables: Some(fixture.to_string()), ..Default::default() }).unwrap();
    let v_3: serde_json::Value = serde_json::from_str(&r_3.json_text).unwrap();
    h.eq_str(gene_by_name(&v_3, "CYP2D6")["diplotype"].as_str().unwrap_or(""), "*1x2/*4", "CYP2D6 *1x2/*4 CN=3 VAF skew");
    h.check(gene_by_name(&v_3, "CYP2D6")["copy_number"].as_f64() == Some(3.0), "CYP2D6 *1x2/*4 CN=3");

    // 18. CYP2D6 *10/*36 exon-9 hybrid drop
    let vcf_10_36 = h.p("v2_10_36.vcf");
    write_text(
        &vcf_10_36,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr22,length=50818468>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA_V2\n\
         chr22\t42126611\t.\tC\tG\t50\tPASS\t.\tGT\t1/1\n\
         chr22\t42129130\t.\tC\tG\t50\tPASS\t.\tGT\t1/1\n\
         chr22\t42130692\t.\tG\tA\t50\tPASS\t.\tGT\t1/1\n",
    );
    let bam_ex9 = h.p("v2_10_36.bam");
    build_v2_synth_bam(
        &bam_ex9,
        &SynthBamOpts { ex9_h1: 1, ex9_h2: 0, body_h1: 1, body_h2: 1, h1_bases: star10, h2_bases: star10, hla_seqs: hla_neg },
    );
    let r_ex9 = run(&vcf_10_36, &h.p("out_v2_ex9"), &PgxOptions { bam: Some(bam_ex9), tables: Some(fixture.to_string()), ..Default::default() }).unwrap();
    let v_ex9: serde_json::Value = serde_json::from_str(&r_ex9.json_text).unwrap();
    h.eq_str(gene_by_name(&v_ex9, "CYP2D6")["diplotype"].as_str().unwrap_or(""), "*10/*36", "CYP2D6 *10/*36 exon-9 hybrid");

    // 19. Negative error paths: GRCh37 BAM, missing .bai, MMGEN_HLA_OFFLINE missing cache
    let bad37_bam = h.p("bad37.bam");
    write_sorted_bam(&bad37_bam, &["@HD\tVN:1.6\tSO:coordinate".into(), "@SQ\tSN:chr1\tLN:249250621".into()]);
    h.expect_err(
        run(&base_vcf, &h.p("out_bad37"), &PgxOptions { bam: Some(bad37_bam), tables: Some(fixture.to_string()), ..Default::default() }),
        &["GRCh38", "249250621"],
        "reject GRCh37 BAM header",
    );
    let bad_qn_bam = h.p("bad_qn.bam");
    write_sorted_bam(&bad_qn_bam, &["@HD\tVN:1.6\tSO:queryname".into(), "@SQ\tSN:chr1\tLN:248956422".into()]);
    h.expect_err(
        run(&base_vcf, &h.p("out_bad_qn"), &PgxOptions { bam: Some(bad_qn_bam), tables: Some(fixture.to_string()), ..Default::default() }),
        &["SO:queryname", "coordinate"],
        "reject SO:queryname BAM header",
    );
    let nobai_bam = h.p("nobai.bam");
    std::fs::copy(&bam_2_4, &nobai_bam).unwrap();
    h.expect_err(
        run(&base_vcf, &h.p("out_nobai"), &PgxOptions { bam: Some(nobai_bam), tables: Some(fixture.to_string()), ..Default::default() }),
        &[".bai", "index"],
        "reject BAM without .bai index",
    );
    let empty_hla = h.dir.join("hla_empty");
    std::fs::create_dir_all(&empty_hla).unwrap();
    {
        let _off = EnvGuard::new(Some(&empty_hla), true);
        h.expect_err(
            run(&base_vcf, &h.p("out_offline"), &PgxOptions { bam: Some(bam_2_4.clone()), tables: Some(fixture.to_string()), ..Default::default() }),
            &["MMGEN_HLA_OFFLINE", "A_gen.fasta"],
            "reject missing IMGT FASTA under MMGEN_HLA_OFFLINE=1",
        );
    }

    // 19b. The tables pin the IPD-IMGT/HLA release (bam_genes.json): a cached file whose
    //      SHA-256 is not the pinned one is never used, and the release is named in the result.
    {
        let a_sha = sha256_hex(&std::fs::read(hla_dir.join("A_gen.fasta")).unwrap());
        let b_sha = sha256_hex(&std::fs::read(hla_dir.join("B_gen.fasta")).unwrap());
        let c_sha = sha256_hex(&std::fs::read(hla_dir.join("C_gen.fasta")).unwrap());
        let pinned_ok = h.dir.join("pinned_ok");
        pin_fixture_hla(&fixture_dir, &pinned_ok, "9.99.0", [&a_sha, &b_sha, &c_sha]);
        let _pin_env = EnvGuard::new(Some(&hla_dir), true);
        let r_pin = run(
            &vcf_2_4,
            &h.p("out_pinned_ok"),
            &PgxOptions { bam: Some(bam_2_4.clone()), tables: Some(pinned_ok.to_string_lossy().to_string()), ..Default::default() },
        );
        h.check(r_pin.is_ok(), "pinned HLA release: the cached files with the pinned SHA-256 are used");
        if let Ok(r_pin) = r_pin {
            let v_pin: serde_json::Value = serde_json::from_str(&r_pin.json_text).unwrap();
            h.check(
                v_pin["data_version"].as_str().is_some_and(|s| s.ends_with("; IPD-IMGT/HLA 9.99.0"))
                    && gene_by_name(&v_pin, "HLA-A")["notes"][0].as_str().is_some_and(|s| s.contains("IPD-IMGT/HLA 9.99.0 (SHA-256")),
                "pinned HLA release named in data_version and in the HLA note",
            );
            h.eq_str(gene_by_name(&v_pin, "HLA-B")["diplotype"].as_str().unwrap_or(""), "*07:02/*57:01", "pinned HLA release: same HLA-B call");
        }
        let pinned_bad = h.dir.join("pinned_bad");
        let zeros = "0".repeat(64);
        pin_fixture_hla(&fixture_dir, &pinned_bad, "9.99.0", [&a_sha, &zeros, &c_sha]);
        h.expect_err(
            run(
                &vcf_2_4,
                &h.p("out_pinned_bad"),
                &PgxOptions { bam: Some(bam_2_4.clone()), tables: Some(pinned_bad.to_string_lossy().to_string()), ..Default::default() },
            ),
            &["HLA-B", "B_gen.fasta", "IPD-IMGT/HLA 9.99.0", "SHA-256 differs", "MMGEN_HLA_OFFLINE"],
            "pinned HLA release: a cached file of another release is refused, with the reason",
        );
    }

    // 20. Real HG002 BAM + VCF + real CPIC bundle integration check (when the data are given):
    //     PGX_V2_BAM, MMGEN_HLA_DIR (IPD-IMGT/HLA FASTAs) and MMG_HUMAN_GENOMES (giab/...), as
    //     trials/python/validate_pgx.py; MMGEN_HLA_DIR is read before the synthetic tests replace it.
    let real_bundle = repo_root.join("mmgen_tools").join("pgx").join("cpic-v1.61.0");
    let hg002_bam = PathBuf::from(std::env::var("PGX_V2_BAM").unwrap_or_default());
    let hg002_hla = PathBuf::from(real_hla_dir.unwrap_or_default());
    let giab = PathBuf::from(std::env::var("MMG_HUMAN_GENOMES").unwrap_or_default()).join("giab").join("HG002_NA24385_son");
    let hg002_v50q = giab.join("v5.0q").join("HG002_GRCh38_v5.0q_smvar.vcf.gz");
    let hg002_v421 = giab.join("v4.2.1").join("HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz");
    if real_bundle.exists() && hg002_bam.is_file() && hg002_hla.is_dir() && hg002_v50q.exists() {
        let _real_env = EnvGuard::new(Some(&hg002_hla), true);
        let mut vcfs = vec![("v50q", hg002_v50q.as_path())];
        if hg002_v421.exists() {
            vcfs.push(("v421", hg002_v421.as_path()));
        }
        for (tag, vcf_p) in vcfs {
            let r_hg = run(
                &vcf_p.to_string_lossy(),
                &h.p(&format!("out_hg002_{tag}")),
                &PgxOptions {
                    bam: Some(hg002_bam.to_string_lossy().to_string()),
                    tables: Some(real_bundle.to_string_lossy().to_string()),
                    ..Default::default()
                },
            )
            .unwrap();
            let v_hg: serde_json::Value = serde_json::from_str(&r_hg.json_text).unwrap();
            h.check(v_hg["genes"].as_array().map(|a| a.len()) == Some(23), &format!("HG002 {tag} v2 returns all 23 genes (19 VCF + 4 BAM)"));
            h.eq_str(gene_by_name(&v_hg, "CYP2D6")["diplotype"].as_str().unwrap_or(""), "*2/*4", &format!("HG002 {tag} CYP2D6 consensus *2/*4"));
            h.check(gene_by_name(&v_hg, "CYP2D6")["copy_number"].as_f64() == Some(2.0), &format!("HG002 {tag} CYP2D6 CN=2"));
            h.eq_str(gene_by_name(&v_hg, "HLA-A")["diplotype"].as_str().unwrap_or(""), "*01:01/*26:01", &format!("HG002 {tag} HLA-A consensus *01:01/*26:01"));
            h.eq_str(gene_by_name(&v_hg, "HLA-B")["diplotype"].as_str().unwrap_or(""), "*35:08/*38:01", &format!("HG002 {tag} HLA-B consensus *35:08/*38:01"));
            h.eq_str(gene_by_name(&v_hg, "HLA-C")["diplotype"].as_str().unwrap_or(""), "*04:01/*12:03", &format!("HG002 {tag} HLA-C consensus *04:01/*12:03"));
            h.check(
                v_hg["data_version"].as_str().is_some_and(|s| s.contains("IPD-IMGT/HLA 3.66.0")),
                &format!("HG002 {tag} data_version names IPD-IMGT/HLA 3.66.0 (the release the tables pin)"),
            );
            let oc_hg = std::fs::read_to_string(v_hg["x_files"]["pharmcat_outside_calls"].as_str().unwrap()).unwrap();
            h.check(
                oc_hg == "CYP2D6\t*2/*4\nHLA-A\t*01:01/*26:01\nHLA-B\t*35:08/*38:01\n",
                &format!("HG002 {tag} outside_calls.tsv exact 3-line TSV"),
            );
        }
    } else {
        println!("NOTE  skipping real HG002 v2 BAM check (set PGX_V2_BAM, MMGEN_HLA_DIR and MMG_HUMAN_GENOMES to run it)");
    }
}

/// Copies the fixture tables and pins an IPD-IMGT/HLA release in bam_genes.json: version and
/// the SHA-256 of A_gen.fasta, B_gen.fasta and C_gen.fasta (in that order); rewrites the manifest.
fn pin_fixture_hla(src: &Path, dst: &Path, version: &str, shas: [&str; 3]) {
    use serde_json::{json, Value};
    copy_dir(src, dst);
    let mut bam: Value = serde_json::from_slice(&std::fs::read(dst.join("bam_genes.json")).unwrap()).unwrap();
    for hg in bam["hla_genes"].as_array_mut().unwrap() {
        let idx = match hg["fasta_filename"].as_str().unwrap() {
            "A_gen.fasta" => 0,
            "B_gen.fasta" => 1,
            _ => 2,
        };
        hg["imgt_hla_version"] = json!(version);
        hg["fasta_sha256"] = json!(shas[idx]);
    }
    let bytes = serde_json::to_vec_pretty(&bam).unwrap();
    std::fs::write(dst.join("bam_genes.json"), &bytes).unwrap();
    let mut manifest: Value = serde_json::from_slice(&std::fs::read(dst.join("manifest.json")).unwrap()).unwrap();
    manifest["files"]["bam_genes.json"]["sha256"] = json!(sha256_hex(&bytes));
    manifest["files"]["bam_genes.json"]["bytes"] = json!(bytes.len());
    std::fs::write(dst.join("manifest.json"), serde_json::to_vec_pretty(&manifest).unwrap()).unwrap();
}
