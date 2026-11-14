#![allow(dead_code)]
// Tests of src/clinico.rs: ClinVar findings, gene regions, variant lookup
// (cascade testing) and the pileup at one position.
//
// Every result is compared with something this file computes on its own:
//  - findings: a brute-force join of structured records kept by the
//    generator (a hash table on chromosome, position, REF, ALT after an
//    independent trimming), with the filters applied by hand;
//  - review stars and classes: tables written from the ClinVar documentation;
//  - gene regions: coordinates worked out by hand for a small GFF and a GTF;
//  - lookup: the status of each person at each variant decided from the
//    generated records, through the plain file and through the tabix index;
//  - pileup: the generator records, for every read, what it shows at each
//    reference position while it builds the CIGAR; the brute force counts
//    from that, never from the CIGAR string.
//
// Usage: cargo run --release --bin clinico
//        cargo run --release --bin clinico -- real <clinvar.vcf.gz> <refseq.gff.gz>
//            <son.vcf.gz> <son.bed> <father.vcf.gz> <father.bed> <mother.vcf.gz> <mother.bed>
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
#[path = "../../../src/clinico.rs"]
mod clinico;

use clinico::*;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::Debug;
use std::io::Write;
use std::path::PathBuf;

// ---------------------------------------------------------------- harness

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
    fn chance(&mut self, pct: u64) -> bool {
        self.below(100) < pct
    }
    fn pick<'a, X>(&mut self, v: &'a [X]) -> &'a X {
        &v[self.below(v.len() as u64) as usize]
    }
    fn shuffle<X>(&mut self, v: &mut [X]) {
        for i in (1..v.len()).rev() {
            let j = self.below(i as u64 + 1) as usize;
            v.swap(i, j);
        }
    }
}

struct T {
    ok: usize,
    falhas: usize,
    dir: PathBuf,
}

impl T {
    fn p(&self, n: &str) -> String {
        self.dir.join(n).to_string_lossy().to_string()
    }
    fn check(&mut self, cond: bool, what: &str) -> bool {
        if cond {
            self.ok += 1;
        } else {
            self.falhas += 1;
            println!("FAIL  {what}");
        }
        cond
    }
    fn eq<A: PartialEq + Debug>(&mut self, got: &A, want: &A, what: &str) -> bool {
        let ok = got == want;
        if !ok {
            println!("FAIL  {what}\n   got:  {got:?}\n   want: {want:?}");
        }
        self.check(ok, what)
    }
    /// Two lists, reporting only the first difference.
    fn lista(&mut self, got: &[String], want: &[String], what: &str) -> bool {
        for (i, (g, w)) in got.iter().zip(want).enumerate() {
            if g != w {
                println!("FAIL  {what}: first difference at item {i}\n   got:  {g}\n   want: {w}");
                return self.check(false, what);
            }
        }
        if got.len() != want.len() {
            println!(
                "FAIL  {what}: {} items, expected {}; first extra: {:?}",
                got.len(),
                want.len(),
                if got.len() > want.len() { got.get(want.len()) } else { want.get(got.len()) }
            );
            return self.check(false, what);
        }
        self.check(true, what)
    }
    fn erro<X>(&mut self, r: Result<X, String>, agulhas: &[&str], what: &str) {
        match r {
            Ok(_) => {
                println!("FAIL  {what}: should have failed");
                self.check(false, what);
            }
            Err(e) => {
                let falta: Vec<&&str> = agulhas.iter().filter(|a| !e.contains(**a)).collect();
                if !falta.is_empty() {
                    println!("FAIL  {what}: message lacks {falta:?}\n   message: {e}");
                }
                self.check(falta.is_empty(), what);
            }
        }
    }
    fn ok_ou<X>(&mut self, r: Result<X, String>, what: &str) -> Option<X> {
        match r {
            Ok(x) => Some(x),
            Err(e) => {
                println!("FAIL  {what}: {e}");
                self.check(false, what);
                None
            }
        }
    }
}

fn grava(p: &str, linhas: &[String]) {
    let mut w = externo::Escritor::cria(p).unwrap();
    for l in linhas {
        w.write_all(l.as_bytes()).unwrap();
        w.write_all(b"\n").unwrap();
    }
    w.termina().unwrap();
}

fn le(p: &str) -> Vec<String> {
    let mut r = externo::abre_leitura(p).unwrap();
    let mut s = String::new();
    std::io::Read::read_to_string(&mut r, &mut s).unwrap();
    s.lines().map(str::to_string).collect()
}

// ---------------------------------------------------------------- genome and alleles

const BASES: &[u8] = b"ACGT";

struct Crom {
    nome: &'static str,
    seq: Vec<u8>,
}

fn genoma(rng: &mut Rng) -> Vec<Crom> {
    [("1", 30_000), ("2", 30_000), ("3", 20_000), ("10", 20_000), ("22", 20_000), ("X", 20_000), ("Y", 8_000), ("MT", 16_569)]
        .iter()
        .map(|&(n, l)| Crom { nome: n, seq: (0..l).map(|_| BASES[rng.below(4) as usize]).collect() })
        .collect()
}

/// The canonical name, computed here independently: chr1 -> 1, chrM/M -> MT.
fn canon(n: &str) -> String {
    let s = n.strip_prefix("chr").unwrap_or(n);
    if s == "M" {
        "MT".to_string()
    } else {
        s.to_string()
    }
}

fn com_chr(c: &str) -> String {
    if c == "MT" {
        "chrM".to_string()
    } else {
        format!("chr{c}")
    }
}

/// Position of a canonical chromosome in the order 1..22, X, Y, MT.
fn ordem(c: &str) -> u32 {
    match c {
        "X" => 1001,
        "Y" => 1002,
        "MT" => 1003,
        s => s.parse().unwrap_or(u32::MAX),
    }
}

fn outra_base(rng: &mut Rng, b: u8) -> u8 {
    loop {
        let x = BASES[rng.below(4) as usize];
        if x != b {
            return x;
        }
    }
}

/// An allele anchored at `pos` on the reference: SNV, deletion, insertion,
/// MNP or a complex change (first and last bases differ).
fn alelo(rng: &mut Rng, seq: &[u8], pos: u64) -> (String, String) {
    let i = pos as usize - 1;
    let s = |a: usize, b: usize| String::from_utf8(seq[a..b].to_vec()).unwrap();
    match rng.below(20) {
        0..=10 => {
            let a = outra_base(rng, seq[i]);
            (s(i, i + 1), (a as char).to_string())
        }
        11..=14 => {
            let d = 1 + rng.below(6) as usize;
            (s(i, i + 1 + d), s(i, i + 1))
        }
        15..=17 => {
            let n = 1 + rng.below(5) as usize;
            let ins: String = (0..n).map(|_| BASES[rng.below(4) as usize] as char).collect();
            (s(i, i + 1), format!("{}{ins}", seq[i] as char))
        }
        18 => {
            let m = 2 + rng.below(2) as usize;
            let r = s(i, i + m);
            let mut a: Vec<u8> = (0..m).map(|_| BASES[rng.below(4) as usize]).collect();
            a[0] = outra_base(rng, seq[i]);
            a[m - 1] = outra_base(rng, seq[i + m - 1]);
            (r, String::from_utf8(a).unwrap())
        }
        _ => {
            let r = s(i, i + 3);
            let a = vec![outra_base(rng, seq[i]), outra_base(rng, seq[i + 2])];
            (r, String::from_utf8(a).unwrap())
        }
    }
}

/// Trimming written independently of the module: shared trailing bases,
/// then shared leading bases (POS moves), at least one base left on each side.
fn apara(pos: u64, r: &str, a: &str) -> Option<(u64, String, String)> {
    let ok = |s: &str| !s.is_empty() && s.chars().all(|c| "ACGTNacgtn".contains(c));
    if !ok(r) || !ok(a) {
        return None;
    }
    let mut r: Vec<char> = r.to_uppercase().chars().collect();
    let mut a: Vec<char> = a.to_uppercase().chars().collect();
    while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
        r.pop();
        a.pop();
    }
    let mut p = pos;
    while r.len() > 1 && a.len() > 1 && r[0] == a[0] {
        r.remove(0);
        a.remove(0);
        p += 1;
    }
    if r == a {
        return None;
    }
    Some((p, r.iter().collect(), a.iter().collect()))
}

/// Zygosity of ALT `k` in a GT, written independently of the module.
fn zig_esperada(gt: Option<&str>, k: u32) -> &'static str {
    let Some(gt) = gt else { return "unknown" };
    let alelos: Vec<&str> = gt.split(|c| c == '/' || c == '|').collect();
    let n = alelos.len();
    let ks = k.to_string();
    let copias = alelos.iter().filter(|x| **x == ks).count();
    if copias > 0 {
        return if n == 1 {
            "hemizygous"
        } else if copias == n {
            "hom"
        } else {
            "het"
        };
    }
    if alelos.iter().any(|x| *x != "0" && *x != ".") {
        return "other_alt";
    }
    if alelos.iter().any(|x| *x == ".") {
        "no_call"
    } else {
        "hom_ref"
    }
}

fn portado(z: &str) -> bool {
    matches!(z, "het" | "hom" | "hemizygous" | "unknown")
}

// ---------------------------------------------------------------- ClinVar: tables and synthetic file

/// Review status -> stars, from https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
const REVS: &[(&str, u8)] = &[
    ("practice_guideline", 4),
    ("reviewed_by_expert_panel", 3),
    ("criteria_provided,_multiple_submitters,_no_conflicts", 2),
    ("criteria_provided,_single_submitter", 1),
    ("criteria_provided,_conflicting_classifications", 1),
    ("no_assertion_criteria_provided", 0),
    ("no_classification_provided", 0),
    ("no_classification_for_the_single_variant", 0),
    ("no_classifications_from_unflagged_records", 0),
    ("some_future_status", 0),
];

const P: &str = "pathogenic_or_likely_pathogenic";
const B: &str = "benign_or_likely_benign";

/// CLNSIG values seen in the ClinVar release of 2026-11-08, with the class
/// tags and the category expected for each.
const SIGS: &[(&str, &[&str], &str)] = &[
    ("Pathogenic", &["pathogenic"], P),
    ("Likely_pathogenic", &["likely_pathogenic"], P),
    ("Pathogenic/Likely_pathogenic", &["pathogenic", "likely_pathogenic"], P),
    ("Pathogenic,_low_penetrance", &["pathogenic_low_penetrance"], P),
    ("Likely_pathogenic,_low_penetrance", &["likely_pathogenic_low_penetrance"], P),
    ("Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance", &["pathogenic", "likely_pathogenic", "pathogenic_low_penetrance"], P),
    ("Pathogenic|risk_factor", &["pathogenic", "risk_factor"], P),
    ("Likely_pathogenic/Likely_risk_allele", &["likely_pathogenic", "likely_risk_allele"], P),
    ("Uncertain_significance", &["uncertain"], "uncertain"),
    ("VUS-high", &["uncertain"], "uncertain"),
    ("Uncertain_significance/Uncertain_risk_allele", &["uncertain", "uncertain_risk_allele"], "uncertain"),
    ("Conflicting_classifications_of_pathogenicity", &["conflicting"], "conflicting"),
    ("Conflicting_classifications_of_pathogenicity|risk_factor", &["conflicting", "risk_factor"], "conflicting"),
    ("Likely_benign", &["likely_benign"], B),
    ("Benign", &["benign"], B),
    ("Benign/Likely_benign", &["benign", "likely_benign"], B),
    ("Benign|drug_response", &["benign", "drug_response"], B),
    ("drug_response", &["drug_response"], "other"),
    ("risk_factor", &["risk_factor"], "other"),
    ("Established_risk_allele", &["established_risk_allele"], "other"),
    ("association_not_found", &["association_not_found"], "other"),
    ("protective", &["protective"], "other"),
    ("Affects", &["affects"], "other"),
    ("confers_sensitivity", &["confers_sensitivity"], "other"),
    ("other", &["other"], "other"),
    ("not_provided", &["not_provided"], "not_classified"),
    ("no_classification_for_the_single_variant", &["not_classified"], "not_classified"),
    ("no_classifications_from_unflagged_records", &["not_classified"], "not_classified"),
];

const CONDICOES: &[&str] = &[
    "Hereditary_breast_ovarian_cancer_syndrome",
    "not_provided",
    "Breast-ovarian_cancer,_familial,_susceptibility_to,_1",
    "Cystic_fibrosis",
    "not_specified",
    "Lynch_syndrome_1",
];
const CONSEQ: &[&str] = &["SO:0001583|missense_variant", "SO:0001627|intron_variant", "SO:0001587|stop_gained", "SO:0001819|synonymous_variant"];
const ONCS: &[&str] = &["Oncogenic", "Likely_oncogenic", "Uncertain_significance"];
const SCIS: &[&str] = &["Tier_I_-_Strong", "Tier_II_-_Potential"];

#[derive(Clone, Debug)]
struct CvRec {
    crom: String,
    pos: u64,
    ref_: String,
    alt: String,
    vid: String,
    aid: String,
    rs: Option<String>,
    genes: Vec<String>,
    sig: Option<usize>,
    rev: Option<usize>,
    conds: Vec<String>,
    mc: Vec<String>,
    tipo: String,
    onc: Option<(String, usize)>,
    sci: Option<(String, usize)>,
    linha: usize,
}

impl CvRec {
    fn info(&self) -> String {
        let mut v = vec![format!("ALLELEID={}", self.aid)];
        if !self.conds.is_empty() {
            v.push(format!("CLNDN={}", self.conds.join("|")));
        }
        if let Some(r) = self.rev {
            v.push(format!("CLNREVSTAT={}", REVS[r].0));
        }
        match self.sig {
            Some(s) => v.push(format!("CLNSIG={}", SIGS[s].0)),
            None => v.push("CLNSIGINCL=424771:Pathogenic".to_string()),
        }
        v.push(format!("CLNVC={}", self.tipo));
        if !self.genes.is_empty() {
            v.push(format!("GENEINFO={}", self.genes.iter().enumerate().map(|(i, g)| format!("{g}:{}", 100 + i)).collect::<Vec<_>>().join("|")));
        }
        if !self.mc.is_empty() {
            v.push(format!("MC={}", self.mc.join(",")));
        }
        if let Some((o, r)) = &self.onc {
            v.push(format!("ONC={o}"));
            v.push(format!("ONCREVSTAT={}", REVS[*r].0));
        }
        v.push("ORIGIN=1".to_string());
        if let Some(r) = &self.rs {
            v.push(format!("RS={r}"));
        }
        if let Some((o, r)) = &self.sci {
            v.push(format!("SCI={o}"));
            v.push(format!("SCIREVSTAT={}", REVS[*r].0));
        }
        v.join(";")
    }

    fn primario(&self) -> bool {
        !self.crom.starts_with("NT_")
    }

    fn estrelas(&self) -> u8 {
        self.rev.map(|r| REVS[r].1).unwrap_or(0)
    }

    fn classes(&self) -> Vec<&'static str> {
        match self.sig {
            Some(s) => SIGS[s].1.to_vec(),
            None => vec!["not_classified"],
        }
    }

    fn categoria(&self) -> &'static str {
        self.sig.map(|s| SIGS[s].2).unwrap_or("not_classified")
    }

    fn significancia(&self) -> String {
        self.sig.map(|s| SIGS[s].0.replace('_', " ")).unwrap_or_else(|| "not classified".to_string())
    }
}

fn gera_clinvar(rng: &mut Rng, g: &[Crom]) -> Vec<CvRec> {
    let mut pesos: Vec<usize> = Vec::new();
    for (i, s) in SIGS.iter().enumerate() {
        let w = if s.2 == P { 6 } else if s.0 == "Uncertain_significance" || s.0 == "Likely_benign" { 6 } else { 1 };
        pesos.extend(std::iter::repeat(i).take(w));
    }
    let mut v = Vec::new();
    let mut vid = 1000u64;
    let anota = |rng: &mut Rng, crom: &str, pos: u64, r: String, a: String, rs: Option<String>, vid: u64| -> CvRec {
        let tipo = if r.len() == 1 && a.len() == 1 {
            "single_nucleotide_variant"
        } else if a.len() == 1 && r.len() > 1 {
            "Deletion"
        } else if r.len() == 1 {
            "Insertion"
        } else if r.len() == a.len() {
            "MNV"
        } else {
            "Indel"
        };
        let mut genes = vec![format!("G{crom}x{}", pos / 3000)];
        if rng.chance(10) {
            genes.push(format!("AS-G{crom}x{}", pos / 3000));
        }
        if rng.chance(3) {
            genes.clear();
        }
        CvRec {
            crom: crom.to_string(),
            pos,
            ref_: r,
            alt: a,
            vid: vid.to_string(),
            aid: (vid + 7_000_000).to_string(),
            rs,
            genes,
            sig: if rng.chance(5) { None } else { Some(*rng.pick(&pesos)) },
            rev: if rng.chance(3) { None } else { Some(rng.below(REVS.len() as u64) as usize) },
            conds: (0..rng.below(3)).map(|_| rng.pick(CONDICOES).to_string()).collect(),
            mc: (0..rng.below(3)).map(|_| rng.pick(CONSEQ).to_string()).collect(),
            tipo: tipo.to_string(),
            onc: if rng.chance(5) { Some((rng.pick(ONCS).to_string(), rng.below(REVS.len() as u64) as usize)) } else { None },
            sci: if rng.chance(3) { Some((rng.pick(SCIS).to_string(), rng.below(REVS.len() as u64) as usize)) } else { None },
            linha: 0,
        }
    };
    for c in g {
        let mut pos = 2u64;
        loop {
            pos += 1 + rng.below(30);
            if pos + 20 > c.seq.len() as u64 {
                break;
            }
            let n = if rng.chance(15) { 2 + rng.below(2) } else { 1 };
            let rs = if rng.chance(60) { Some((rng.below(900_000_000) + 1).to_string()) } else { None };
            let mut vistos: Vec<(String, String)> = Vec::new();
            for _ in 0..n {
                let (r, a) = alelo(rng, &c.seq, pos);
                if vistos.contains(&(r.clone(), a.clone())) {
                    continue;
                }
                vistos.push((r.clone(), a.clone()));
                vid += 1 + rng.below(5);
                v.push(anota(rng, c.nome, pos, r, a, rs.clone(), vid));
            }
        }
    }
    for i in 0..5u64 {
        vid += 1;
        v.push(anota(rng, "NT_187633.1", 100 + 10 * i, "A".into(), "G".into(), None, vid));
    }
    v
}

const CV_CABECALHO: &[&str] = &[
    "##fileformat=VCFv4.1",
    "##fileDate=2026-11-08",
    "##source=ClinVar",
    "##reference=GRCh38",
    "##ID=<Description=\"ClinVar Variation ID\">",
    "##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description=\"the ClinVar Allele ID\">",
    "##INFO=<ID=CLNDN,Number=.,Type=String,Description=\"ClinVar's preferred disease name\">",
    "##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description=\"ClinVar review status of germline classification\">",
    "##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"Aggregate germline classification\">",
    "##INFO=<ID=CLNSIGINCL,Number=.,Type=String,Description=\"Germline classification for a haplotype or genotype that includes this variant\">",
    "##INFO=<ID=CLNVC,Number=1,Type=String,Description=\"Variant type\">",
    "##INFO=<ID=GENEINFO,Number=1,Type=String,Description=\"Gene(s) for the variant\">",
    "##INFO=<ID=MC,Number=.,Type=String,Description=\"molecular consequence\">",
    "##INFO=<ID=ONC,Number=.,Type=String,Description=\"Aggregate oncogenicity classification\">",
    "##INFO=<ID=ONCREVSTAT,Number=.,Type=String,Description=\"review status of oncogenicity classification\">",
    "##INFO=<ID=ORIGIN,Number=.,Type=String,Description=\"Allele origin\">",
    "##INFO=<ID=RS,Number=.,Type=String,Description=\"dbSNP ID\">",
    "##INFO=<ID=SCI,Number=.,Type=String,Description=\"Aggregate somatic clinical impact\">",
    "##INFO=<ID=SCIREVSTAT,Number=.,Type=String,Description=\"review status of somatic clinical impact\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
];

fn grava_clinvar(p: &str, cv: &mut [CvRec]) {
    let mut l: Vec<String> = CV_CABECALHO.iter().map(|s| s.to_string()).collect();
    for r in cv.iter_mut() {
        l.push(format!("{}\t{}\t{}\t{}\t{}\t.\t.\t{}", r.crom, r.pos, r.vid, r.ref_, r.alt, r.info()));
        r.linha = l.len();
    }
    grava(p, &l);
}

// ---------------------------------------------------------------- person records

#[derive(Clone, Debug)]
struct PRec {
    /// Canonical name for the primary chromosomes; the name itself otherwise.
    crom: String,
    outra: bool,
    pos: u64,
    ref_: String,
    alts: Vec<String>,
    info: String,
    /// GT of each sample, and the full sample field.
    gts: Vec<String>,
    campos: Vec<String>,
    format: String,
    linha: usize,
}

impl PRec {
    fn texto(&self, nome: &str) -> String {
        let mut s = format!("{nome}\t{}\t.\t{}\t{}\t50\tPASS\t{}", self.pos, self.ref_, self.alts.join(","), self.info);
        if !self.campos.is_empty() {
            s.push('\t');
            s.push_str(&self.format);
            for c in &self.campos {
                s.push('\t');
                s.push_str(c);
            }
        }
        s
    }
    fn gt(&self, col: Option<usize>) -> Option<&str> {
        col.filter(|_| !self.campos.is_empty()).map(|j| self.gts[j].as_str())
    }
    fn fim(&self) -> u64 {
        let mut s1 = self.pos + (self.ref_.len() as u64).max(1);
        if let Some(e) = self.info.split(';').find_map(|kv| kv.strip_prefix("END=")).and_then(|e| e.parse::<u64>().ok()) {
            s1 = s1.max(e + 1);
        }
        s1
    }
    /// The value of a FORMAT key for sample `j`, "." as missing.
    fn valor(&self, j: Option<usize>, chave: &str) -> Option<String> {
        let j = j?;
        let i = self.format.split(':').position(|f| f == chave)?;
        let v = self.campos.get(j)?.split(':').nth(i)?;
        if v.is_empty() || v == "." {
            None
        } else {
            Some(v.to_string())
        }
    }
}

fn escolhe_gt(rng: &mut Rng, n_alts: usize, haploide: bool) -> String {
    let um: &[&str] = &["0/1", "0/1", "0/1", "0|1", "1|0", "1/1", "1/1", "0/0", "./.", "0/.", "1/."];
    let dois: &[&str] = &["1/2", "1/2", "0/1", "0/2", "2/2", "1/1", "./.", "2|1", "1/.", "0/0"];
    if haploide && rng.chance(40) {
        let h: &[&str] = if n_alts == 1 { &["1", "1", "0", "."] } else { &["1", "2", "0"] };
        return rng.pick(h).to_string();
    }
    rng.pick(if n_alts == 1 { um } else { dois }).to_string()
}

fn campos_amostra(rng: &mut Rng, gts: &[String], n_alts: usize) -> (String, Vec<String>) {
    if rng.chance(20) {
        return ("GT".to_string(), gts.to_vec());
    }
    let campos = gts
        .iter()
        .map(|gt| {
            let dp = if rng.chance(5) { ".".to_string() } else { (5 + rng.below(200)).to_string() };
            let ad: Vec<String> = (0..=n_alts).map(|_| rng.below(80).to_string()).collect();
            format!("{gt}:{dp}:{}:{}", rng.below(99), ad.join(","))
        })
        .collect();
    ("GT:DP:GQ:AD".to_string(), campos)
}

/// Records of a person (or a family) on every chromosome of `g`: alleles
/// taken from ClinVar (`prob_cv` percent of them) and new ones, some combined
/// into multiallelic records, some padded with reference bases at either end.
fn gera_registros(rng: &mut Rng, g: &[Crom], cv: &[CvRec], prob_cv: u64, n_amostras: usize, evita: &[(&str, u64, u64)]) -> Vec<PRec> {
    let mut todos = Vec::new();
    for c in g {
        let mut alelos: BTreeMap<u64, Vec<(String, String)>> = BTreeMap::new();
        for r in cv.iter().filter(|r| r.crom == c.nome) {
            if rng.chance(prob_cv) {
                alelos.entry(r.pos).or_default().push((r.ref_.clone(), r.alt.clone()));
            }
        }
        for _ in 0..c.seq.len() / 120 {
            let pos = 2 + rng.below(c.seq.len() as u64 - 30);
            let a = alelo(rng, &c.seq, pos);
            alelos.entry(pos).or_default().push(a);
        }
        let haploide = matches!(c.nome, "X" | "Y" | "MT");
        let mut regs: Vec<PRec> = Vec::new();
        for (pos, mut v) in alelos {
            if evita.iter().any(|&(n, a, b)| n == c.nome && pos + 10 >= a && pos < b) {
                continue;
            }
            v.sort();
            v.dedup();
            rng.shuffle(&mut v);
            while !v.is_empty() {
                let n = if v.len() >= 2 && rng.chance(50) { 2 } else { 1 };
                let parte: Vec<(String, String)> = v.drain(..n).collect();
                let rr = parte.iter().map(|x| x.0.clone()).max_by_key(|x| x.len()).unwrap();
                let mut alts: Vec<String> = parte.iter().map(|(r, a)| format!("{a}{}", &rr[r.len()..])).collect();
                let mut rr = rr;
                let mut p = pos;
                if rng.chance(15) {
                    let k = 1 + rng.below(2) as usize;
                    let ini = p as usize - 1 + rr.len();
                    if ini + k <= c.seq.len() {
                        let extra = String::from_utf8(c.seq[ini..ini + k].to_vec()).unwrap();
                        rr.push_str(&extra);
                        for a in &mut alts {
                            a.push_str(&extra);
                        }
                    }
                }
                if rng.chance(15) && p >= 2 {
                    p -= 1;
                    let b = c.seq[p as usize - 1] as char;
                    rr.insert(0, b);
                    for a in &mut alts {
                        a.insert(0, b);
                    }
                }
                let gts: Vec<String> = (0..n_amostras).map(|_| escolhe_gt(rng, alts.len(), haploide)).collect();
                let (format, campos) = campos_amostra(rng, &gts, alts.len());
                regs.push(PRec { crom: c.nome.to_string(), outra: false, pos: p, ref_: rr, alts, info: "DP=40".into(), gts, campos, format, linha: 0 });
            }
        }
        // A record without ALT, and a symbolic deletion.
        for (alt, info) in [(".", "DP=3".to_string()), ("<DEL>", "SVTYPE=DEL;END=0".to_string())] {
            let pos = 100 + rng.below(c.seq.len() as u64 - 500);
            if evita.iter().any(|&(n, a, b)| n == c.nome && pos + 300 >= a && pos < b) {
                continue;
            }
            let info = if alt == "<DEL>" { format!("SVTYPE=DEL;END={}", pos + 40) } else { info };
            let gts: Vec<String> = (0..n_amostras).map(|_| if alt == "." { "0/0".to_string() } else { "0/1".to_string() }).collect();
            let (format, campos) = campos_amostra(rng, &gts, 1);
            regs.push(PRec {
                crom: c.nome.to_string(),
                outra: false,
                pos,
                ref_: (c.seq[pos as usize - 1] as char).to_string(),
                alts: vec![alt.to_string()],
                info,
                gts,
                campos,
                format,
                linha: 0,
            });
        }
        regs.sort_by_key(|r| r.pos);
        todos.extend(regs);
    }
    todos
}

fn registro_outra(rng: &mut Rng, nome: &str, pos: u64, r: &str, a: &str, n_amostras: usize) -> PRec {
    let gts: Vec<String> = (0..n_amostras).map(|_| "0/1".to_string()).collect();
    let (format, campos) = campos_amostra(rng, &gts, 1);
    PRec { crom: nome.to_string(), outra: true, pos, ref_: r.into(), alts: vec![a.into()], info: ".".into(), gts, campos, format, linha: 0 }
}

/// Writes the records in the order of `ordem` (canonical names; others by
/// name), with `chr` names or not; returns them in file order with the line
/// numbers set.
fn grava_pessoa(p: &str, regs: &[PRec], ordem: &[&str], chr: bool, amostras: &[&str], montagem: &str) -> Vec<PRec> {
    let mut l: Vec<String> = vec!["##fileformat=VCFv4.2".into()];
    l.push(format!("##contig=<ID={},length={}>", if chr { "chr1" } else { "1" }, if montagem == "GRCh37" { 249250621 } else { 248956422 }));
    l.push("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">".into());
    l.push("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".into());
    let mut cab = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string();
    if !amostras.is_empty() {
        cab.push_str("\tFORMAT\t");
        cab.push_str(&amostras.join("\t"));
    }
    l.push(cab);
    let mut saida = Vec::new();
    for o in ordem {
        for r in regs.iter().filter(|r| r.crom == *o) {
            let mut r = r.clone();
            if amostras.is_empty() {
                r.campos.clear();
                r.gts.clear();
            }
            let nome = if r.outra { r.crom.clone() } else if chr { com_chr(&r.crom) } else { r.crom.clone() };
            l.push(r.texto(&nome));
            r.linha = l.len();
            saida.push(r);
        }
    }
    grava(p, &l);
    saida
}

// ---------------------------------------------------------------- 1. stars and classes

fn testa_estrelas(t: &mut T) {
    for &(s, e) in REVS.iter().filter(|r| r.0 != "some_future_status") {
        t.eq(&review_stars(s), &Some(e), &format!("review stars of {s}"));
    }
    for (s, e) in [
        ("criteria provided, single submitter", Some(1u8)),
        ("Criteria_Provided,_Multiple_Submitters,_No_Conflicts", Some(2)),
        ("criteria_provided,_multiple_submitters", Some(2)),
        ("criteria_provided,_conflicting_interpretations", Some(1)),
        ("no_assertion_provided", Some(0)),
        ("no_interpretation_for_the_single_variant", Some(0)),
        ("reviewed by expert panel", Some(3)),
        ("some_future_status", None),
        ("", None),
    ] {
        t.eq(&review_stars(s), &e, &format!("review stars of \"{s}\""));
    }
}

fn testa_classes(t: &mut T) {
    for &(s, tags, cat) in SIGS {
        let got = significance_classes(Some(s));
        t.eq(&got, &tags.to_vec(), &format!("classes of {s}"));
        t.eq(&significance_category(&got), &cat, &format!("category of {s}"));
    }
    t.eq(&significance_classes(None), &vec!["not_classified"], "classes without CLNSIG");
    t.eq(&significance_category(&["not_classified"]), &"not_classified", "category without CLNSIG");
}

// ---------------------------------------------------------------- 2. findings against a brute-force join

struct Filtro {
    classes: Option<Vec<&'static str>>,
    min: u8,
    genes: Vec<String>,
    regioes: Option<Vec<(String, u64, u64)>>,
}

impl Filtro {
    fn todos() -> Filtro {
        Filtro { classes: None, min: 0, genes: vec![], regioes: None }
    }
    fn passa(&self, c: &CvRec, pos: u64, r: &str) -> bool {
        if let Some(cl) = &self.classes {
            if !c.classes().iter().any(|x| cl.contains(x)) {
                return false;
            }
        }
        if c.estrelas() < self.min {
            return false;
        }
        if !self.genes.is_empty() && !c.genes.iter().any(|g| self.genes.iter().any(|x| x.eq_ignore_ascii_case(g))) {
            return false;
        }
        if let Some(rs) = &self.regioes {
            let (a, b) = (pos - 1, pos - 1 + r.len() as u64);
            if !rs.iter().any(|(n, x, y)| *n == c.crom && *x < b && a < *y) {
                return false;
            }
        }
        true
    }
}

fn opcoes(f: &Filtro, classes: &[&str], regioes: Option<&str>, saida: Option<&str>) -> FindingsOptions {
    FindingsOptions {
        classes: classes.iter().map(|s| s.to_string()).collect(),
        min_stars: f.min,
        genes: f.genes.clone(),
        regions: regioes.map(str::to_string),
        output: saida.map(str::to_string),
    }
}

#[derive(Default, Debug, PartialEq)]
struct Contadores {
    records: usize,
    alleles: usize,
    not_carried: usize,
    outside: usize,
    other: usize,
    not_comparable: usize,
    matched: usize,
    by_category: Vec<(String, usize)>,
    by_significance: Vec<(String, usize)>,
    plp: [usize; 4],
}

fn linha_achado_esperado(c: &CvRec, pr: &PRec, nome: &str, chave: &(u64, String, String), zig: &str, gt: &str) -> String {
    let conds: Vec<String> = c.conds.iter().map(|x| x.replace('_', " ")).collect();
    let mc: Vec<String> = c.mc.iter().map(|x| x.split('|').nth(1).unwrap().to_string()).collect();
    format!(
        "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{:?}|{}|{:?}|{:?}|{:?}|{:?}|{}|{}|{}|{}|{}",
        nome,
        chave.0,
        chave.1,
        chave.2,
        c.vid,
        zig,
        gt,
        c.estrelas(),
        c.significancia(),
        c.genes.join(","),
        conds.join(";"),
        mc.join(","),
        c.rs.as_ref().map(|r| format!("rs{r}")),
        c.tipo.replace('_', " "),
        c.onc.as_ref().map(|o| o.0.replace('_', " ")),
        c.onc.as_ref().map(|o| REVS[o.1].1),
        c.sci.as_ref().map(|o| o.0.replace('_', " ")),
        c.sci.as_ref().map(|o| REVS[o.1].1),
        c.rev.map(|r| REVS[r].0.replace('_', " ")).unwrap_or_default(),
        c.categoria(),
        c.classes().join(","),
        pr.linha,
        c.linha
    )
}

fn linha_achado(f: &Finding) -> String {
    format!(
        "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{:?}|{}|{:?}|{:?}|{:?}|{:?}|{}|{}|{}|{}|{}",
        f.chrom,
        f.pos,
        f.ref_allele,
        f.alt_allele,
        f.variation_id,
        f.zygosity,
        f.genotype,
        f.stars,
        f.significance,
        f.genes.join(","),
        f.conditions.join(";"),
        f.consequences.join(","),
        f.rsid,
        f.variant_type.clone().unwrap_or_default(),
        f.oncogenicity,
        f.oncogenicity_stars,
        f.somatic_impact,
        f.somatic_impact_stars,
        f.review_status,
        f.category,
        f.classes.join(","),
        f.person_line,
        f.clinvar_line
    )
}

/// The findings and counters the module must give, from the generated
/// records: a hash join, in the order 1..22, X, Y, MT.
fn esperado(cv: &[CvRec], regs: &[PRec], col: Option<usize>, chr: bool, f: &Filtro, chamaveis: Option<&[(String, u64, u64)]>) -> (Vec<String>, Contadores) {
    let mut mapa: HashMap<(String, u64, String, String), Vec<&CvRec>> = HashMap::new();
    for c in cv.iter().filter(|c| c.primario()) {
        let (p, r, a) = apara(c.pos, &c.ref_, &c.alt).unwrap();
        mapa.entry((c.crom.clone(), p, r, a)).or_default().push(c);
    }
    let mut k = Contadores { records: regs.len(), ..Default::default() };
    let mut achados: Vec<(u32, u64, String, String, u64, String)> = Vec::new();
    let mut por_cat: HashMap<&str, usize> = HashMap::new();
    let mut por_sig: HashMap<String, usize> = HashMap::new();
    for pr in regs {
        if pr.outra {
            k.other += 1;
            continue;
        }
        let gt = pr.gt(col);
        for (i, alt) in pr.alts.iter().enumerate() {
            if alt == "." {
                continue;
            }
            let Some(chave) = apara(pr.pos, &pr.ref_, alt) else {
                k.not_comparable += 1;
                continue;
            };
            let z = if col.is_none() { "unknown" } else { zig_esperada(gt, i as u32 + 1) };
            if !portado(z) {
                k.not_carried += 1;
                continue;
            }
            if let Some(b) = chamaveis {
                let (a, f) = (chave.0 - 1, chave.0 - 1 + chave.1.len() as u64);
                if !b.iter().any(|(n, x, y)| *n == pr.crom && *x <= a && f <= *y) {
                    k.outside += 1;
                    continue;
                }
            }
            k.alleles += 1;
            let Some(cs) = mapa.get(&(pr.crom.clone(), chave.0, chave.1.clone(), chave.2.clone())) else { continue };
            for c in cs {
                k.matched += 1;
                *por_cat.entry(c.categoria()).or_default() += 1;
                *por_sig.entry(c.significancia()).or_default() += 1;
                if !f.passa(c, chave.0, &chave.1) {
                    continue;
                }
                if c.categoria() == P {
                    k.plp[match z {
                        "het" => 0,
                        "hom" => 1,
                        "hemizygous" => 2,
                        _ => 3,
                    }] += 1;
                }
                let nome = if chr { com_chr(&pr.crom) } else { pr.crom.clone() };
                let linha = linha_achado_esperado(c, pr, &nome, &chave, z, gt.unwrap_or(""));
                achados.push((ordem(&pr.crom), chave.0, chave.1.clone(), chave.2.clone(), c.vid.parse().unwrap(), linha));
            }
        }
    }
    achados.sort();
    k.by_category = CATEGORIES.iter().map(|c| (c.to_string(), por_cat.get(c).copied().unwrap_or(0))).collect();
    let mut s: Vec<(String, usize)> = por_sig.into_iter().collect();
    s.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    k.by_significance = s;
    (achados.into_iter().map(|x| x.5).collect(), k)
}

fn contadores(r: &FindingsReport) -> Contadores {
    Contadores {
        records: r.person_records,
        alleles: r.person_alleles,
        not_carried: r.person_not_carried,
        outside: r.person_outside_callable,
        other: r.person_other_sequences,
        not_comparable: r.person_not_comparable,
        matched: r.matched,
        by_category: r.by_category.clone(),
        by_significance: r.by_significance.clone(),
        plp: [r.plp_heterozygous, r.plp_homozygous, r.plp_hemizygous, r.plp_unknown_zygosity],
    }
}

#[allow(clippy::too_many_arguments)]
fn confere_achados(t: &mut T, nome: &str, pessoa: &Person, cv_path: &str, cv: &[CvRec], regs: &[PRec], col: Option<usize>, chr: bool, f: &Filtro, opt: &FindingsOptions) -> Option<FindingsReport> {
    confere_achados_em(t, nome, pessoa, cv_path, cv, regs, col, chr, f, opt, None)
}

#[allow(clippy::too_many_arguments)]
fn confere_achados_em(t: &mut T, nome: &str, pessoa: &Person, cv_path: &str, cv: &[CvRec], regs: &[PRec], col: Option<usize>, chr: bool, f: &Filtro, opt: &FindingsOptions, chamaveis: Option<&[(String, u64, u64)]>) -> Option<FindingsReport> {
    let r = t.ok_ou(clinvar_findings(pessoa, cv_path, opt), &format!("findings, {nome}"))?;
    let (esp, k) = esperado(cv, regs, col, chr, f, chamaveis);
    let got: Vec<String> = r.findings.iter().map(linha_achado).collect();
    t.lista(&got, &esp, &format!("findings, {nome}: equal to the brute-force join ({} findings)", esp.len()));
    t.eq(&contadores(&r), &k, &format!("findings, {nome}: counters"));
    t.check(r.clinvar_records == cv.len() && r.clinvar_other_sequences == 5, &format!("findings, {nome}: ClinVar records counted"));
    t.eq(&r.unknown_review_status, &if r.matched > 0 && cv.iter().any(|c| c.rev == Some(REVS.len() - 1)) { r.unknown_review_status.clone() } else { vec![] }, "unknown review statuses");
    for x in &r.unknown_review_status {
        t.check(x == "some_future_status", &format!("findings, {nome}: only the made-up review status is unknown ({x})"));
    }
    let json = r.json();
    t.check(json.starts_with('{') && json.contains("\"findings\":[") && json.contains("not a diagnosis"), &format!("findings, {nome}: JSON"));
    Some(r)
}

fn testa_achados(t: &mut T, rng: &mut Rng, g: &[Crom], cv: &[CvRec], cv_path: &str) {
    let regs = gera_registros(rng, g, cv, 35, 3, &[]);
    let n_por = regs.len();
    let mut com_outras = regs.clone();
    // Records on other sequences: one exactly like a ClinVar record on NT_187633.1.
    com_outras.push(registro_outra(rng, "chr1_KI270706v1_random", 500, "A", "G", 3));
    com_outras.push(registro_outra(rng, "chrUn_KI270742v1", 77, "C", "T", 3));
    let nt = cv.iter().find(|c| !c.primario()).unwrap();
    com_outras.push(registro_outra(rng, "NT_187633.1", nt.pos, &nt.ref_, &nt.alt, 3));
    t.check(n_por > 1000, "enough person records generated");

    // A: chr names, natural order (chrM before chrX, as sort_genomic natural
    // writes them), others in between; three samples, reading KID; plain text.
    let pa = t.p("pessoa_a.vcf");
    let ordem_a = ["1", "chr1_KI270706v1_random", "2", "3", "10", "22", "MT", "chrUn_KI270742v1", "X", "Y", "NT_187633.1"];
    let regs_a = grava_pessoa(&pa, &com_outras, &ordem_a, true, &["MOM", "KID", "DAD"], "GRCh38");
    let pessoa_a = Person { vcf: pa.clone(), sample: Some("KID".into()), callable: None };
    let f = Filtro::todos();
    if let Some(r) = confere_achados(t, "chr names, MT before X, sample KID of 3", &pessoa_a, cv_path, cv, &regs_a, Some(1), true, &f, &opcoes(&f, &[], None, None)) {
        t.eq(&r.sample, &Some("KID".to_string()), "findings: the sample read");
        t.check(r.findings.iter().any(|x| x.zygosity == "hemizygous") && r.findings.iter().any(|x| x.zygosity == "hom"), "findings: the data has hom and hemizygous matches");
        t.check(r.findings.iter().any(|x| x.chrom == "chrM"), "findings: mitochondrial matches are found");
        t.eq(&r.assembly, &Some("GRCh38".to_string()), "findings: assembly from the headers");
    }
    // The same file, another sample.
    let pessoa_mom = Person { vcf: pa.clone(), sample: Some("MOM".into()), callable: None };
    confere_achados(t, "sample MOM", &pessoa_mom, cv_path, cv, &regs_a, Some(0), true, &f, &opcoes(&f, &[], None, None));

    // B: no chr, MT last, BGZF, one sample; P/LP with at least one star, TSV.
    let pb = t.p("pessoa_b.vcf.gz");
    let ordem_b = ["1", "2", "3", "10", "22", "X", "Y", "MT", "NT_187633.1", "chrUn_KI270742v1"];
    let so_um: Vec<PRec> = com_outras.iter().map(|r| PRec { gts: vec![r.gts[2].clone()], campos: vec![r.campos[2].clone()], ..r.clone() }).collect();
    let regs_b = grava_pessoa(&pb, &so_um, &ordem_b, false, &["HG002"], "GRCh38");
    let pessoa_b = Person { vcf: pb.clone(), sample: None, callable: None };
    let f = Filtro { classes: Some(PLP_TAGS.to_vec()), min: 1, ..Filtro::todos() };
    let tsv = t.p("achados.tsv");
    if let Some(r) = confere_achados(t, "no chr, MT last, BGZF, plp >= 1 star", &pessoa_b, cv_path, cv, &regs_b, Some(0), false, &f, &opcoes(&f, &["plp"], None, Some(&tsv))) {
        let linhas = le(&tsv);
        t.eq(&r.written, &Some(r.findings.len()), "TSV: lines written");
        let esp: Vec<String> = std::iter::once(linhas.first().cloned().unwrap_or_default()).chain(r.findings.iter().map(|f| f.tsv())).collect();
        t.lista(&linhas, &esp, "TSV: header and one line per finding");
        t.check(linhas.first().is_some_and(|l| l.starts_with("#chrom\tpos\tref\talt\tgenes")), "TSV: header");
        let ok = linhas.iter().skip(1).zip(&r.findings).all(|(l, f)| {
            let c: Vec<&str> = l.split('\t').collect();
            c.len() == 20 && c[0] == f.chrom && c[1] == f.pos.to_string() && c[2] == f.ref_allele && c[3] == f.alt_allele && c[14] == f.variation_id
        });
        t.check(ok, "TSV: 20 columns; chrom, pos, ref, alt, VariationID in place");
        t.check(r.plp_heterozygous + r.plp_homozygous + r.plp_hemizygous > 0, "findings: P/LP matches by zygosity");
        t.check(r.findings.iter().all(|f| f.stars >= 1 && f.category == P), "findings: filters applied");
    }
    // Class names given in other spellings select the same.
    let f2 = Filtro { classes: Some(PLP_TAGS.to_vec()), min: 1, ..Filtro::todos() };
    confere_achados(t, "class names with spaces", &pessoa_b, cv_path, cv, &regs_b, Some(0), false, &f2, &opcoes(&f2, &["Pathogenic", "likely pathogenic", "Pathogenic, low penetrance", "LIKELY-PATHOGENIC-LOW-PENETRANCE"], None, None));

    // Uncertain and conflicting, two stars, a few genes (any case).
    let genes: Vec<String> = cv.iter().filter(|c| c.crom == "2" || c.crom == "X").flat_map(|c| c.genes.clone()).take(4000).collect::<HashSet<_>>().into_iter().take(6).map(|g| g.to_lowercase()).collect();
    let f = Filtro { classes: Some(vec!["uncertain", "conflicting"]), min: 2, genes: genes.clone(), regioes: None };
    confere_achados(t, "uncertain + conflicting, >= 2 stars, 6 genes", &pessoa_a, cv_path, cv, &regs_a, Some(1), true, &f, &opcoes(&f, &["vus", "conflicting"], None, None));

    // A panel as BED (chr names, while the person file has none).
    let mut regioes: Vec<(String, u64, u64)> = Vec::new();
    let mut bed: Vec<String> = vec!["track name=panel".into()];
    for c in g {
        let mut p = 0u64;
        while p + 3000 < c.seq.len() as u64 {
            p += rng.below(2500);
            let fim = p + 1 + rng.below(1500);
            regioes.push((c.nome.to_string(), p, fim));
            bed.push(format!("{}\t{p}\t{fim}\tgene{}", com_chr(c.nome), regioes.len()));
            p = fim;
        }
    }
    let bed_p = t.p("painel.bed");
    grava(&bed_p, &bed);
    let f = Filtro { classes: Some(vec!["pathogenic", "likely_benign"]), min: 0, genes: vec![], regioes: Some(regioes) };
    confere_achados(t, "panel BED, pathogenic + likely_benign", &pessoa_b, cv_path, cv, &regs_b, Some(0), false, &f, &opcoes(&f, &["pathogenic", "likely_benign"], Some(&bed_p), None));

    // Callable regions (as GIAB gives): only alleles entirely inside count.
    let cham = regioes_chamaveis(rng, g);
    let cham_p = t.p("chamaveis.bed");
    grava_bed(&cham_p, &cham, true);
    let pessoa_c = Person { vcf: pb.clone(), sample: None, callable: Some(cham_p.clone()) };
    let f = Filtro::todos();
    if let Some(r) = confere_achados_em(t, "callable BED", &pessoa_c, cv_path, cv, &regs_b, Some(0), false, &f, &opcoes(&f, &[], None, None), Some(&cham)) {
        t.check(r.person_outside_callable > 100 && r.callable.as_deref() == Some(cham_p.as_str()), "findings: alleles outside the callable regions left out and counted");
    }

    // Sites-only VCF: every allele listed counts, zygosity unknown.
    let ps = t.p("sitios.vcf");
    let regs_s = grava_pessoa(&ps, &regs, &ordem_b, false, &[], "GRCh38");
    let f = Filtro::todos();
    if let Some(r) = confere_achados(t, "sites-only VCF", &Person { vcf: ps.clone(), sample: None, callable: None }, cv_path, cv, &regs_s, None, false, &f, &opcoes(&f, &[], None, None)) {
        t.check(r.sample.is_none() && r.findings.iter().all(|f| f.zygosity == "unknown"), "sites-only VCF: zygosity unknown");
    }

    // ---- errors
    let um = |p: &str, s: Option<&str>| Person { vcf: p.to_string(), sample: s.map(str::to_string), callable: None };
    let todos = FindingsOptions::default();
    t.erro(clinvar_findings(&um(&pa, None), cv_path, &todos), &["3 samples", "MOM, KID, DAD", "sample="], "error: several samples and none chosen");
    t.erro(clinvar_findings(&um(&pa, Some("kid")), cv_path, &todos), &["\"kid\"", "MOM, KID, DAD", "case-sensitive"], "error: sample not in the file");
    t.erro(clinvar_findings(&um(&ps, Some("X")), cv_path, &todos), &["sites-only"], "error: sample asked of a sites-only VCF");
    t.erro(clinvar_findings(&um(cv_path, None), &pa, &todos), &["does not look like the ClinVar VCF", "CLNSIG"], "error: person and ClinVar swapped");
    let nao_existe = t.p("nao_existe.vcf");
    t.erro(clinvar_findings(&um(&nao_existe, None), cv_path, &todos), &["nao_existe.vcf"], "error: missing file named");
    let o = FindingsOptions { classes: vec!["pathogenic".into(), "maybe_bad".into()], ..Default::default() };
    t.erro(clinvar_findings(&um(&pb, None), cv_path, &o), &["maybe_bad", "likely_pathogenic", "plp"], "error: unknown class");
    let o = FindingsOptions { min_stars: 5, ..Default::default() };
    t.erro(clinvar_findings(&um(&pb, None), cv_path, &o), &["0 to 4"], "error: five stars");
    let o = FindingsOptions { output: Some(pb.clone()), ..Default::default() };
    t.erro(clinvar_findings(&um(&pb, None), cv_path, &o), &["input file", "destroy"], "error: output over the input");

    // Lexical order (1, 10, 2): the first record of chromosome 2 is refused.
    let pl = t.p("lexica.vcf");
    let regs_l = grava_pessoa(&pl, &regs, &["1", "10", "2", "22", "3", "MT", "X", "Y"], false, &["S"], "GRCh38");
    let primeiro_2 = regs_l.iter().find(|r| r.crom == "2").unwrap();
    t.erro(
        clinvar_findings(&um(&pl, None), cv_path, &todos),
        &[&format!("line {}", primeiro_2.linha), "2:", "comes after 10:", "sort_genomic", "chrom_order=\"natural\"", "lexica.vcf"],
        "error: lexical chromosome order, with the line and the way out",
    );
    // Two records of chromosome 3 swapped.
    let mut trocados = regs.clone();
    let i3: Vec<usize> = trocados.iter().enumerate().filter(|(_, r)| r.crom == "3").map(|(i, _)| i).collect();
    let (a, b) = (i3[40], i3[41]);
    if trocados[a].pos != trocados[b].pos {
        trocados.swap(a, b);
    } else {
        trocados.swap(a, i3[45]);
    }
    let pt = t.p("trocados.vcf");
    let regs_t = grava_pessoa(&pt, &trocados, &ordem_b, false, &["S"], "GRCh38");
    let ruim = regs_t.windows(2).find(|w| w[0].crom == "3" && w[1].crom == "3" && w[1].pos < w[0].pos).map(|w| w[1].linha).unwrap();
    t.erro(clinvar_findings(&um(&pt, None), cv_path, &todos), &[&format!("line {ruim}"), "sort_genomic"], "error: positions out of order, with the line");
    // ClinVar out of order.
    let mut cv_ruim: Vec<CvRec> = cv.to_vec();
    let iy = cv_ruim.iter().position(|c| c.crom == "Y").unwrap();
    let ix = cv_ruim.iter().position(|c| c.crom == "X").unwrap();
    let y = cv_ruim.remove(iy);
    cv_ruim.insert(ix, y);
    let cvr = t.p("clinvar_ruim.vcf");
    grava_clinvar(&cvr, &mut cv_ruim);
    t.erro(clinvar_findings(&um(&pb, None), &cvr, &todos), &["clinvar_ruim.vcf", "comes after Y", "sort_genomic"], "error: ClinVar out of order");
    // GRCh37 person against GRCh38 ClinVar.
    let p37 = t.p("grch37.vcf");
    grava_pessoa(&p37, &regs, &ordem_b, false, &["S"], "GRCh37");
    t.erro(clinvar_findings(&um(&p37, None), cv_path, &todos), &["GRCh37", "GRCh38", "grch37.vcf"], "error: assemblies differ");
}

const PLP_TAGS: [&str; 4] = ["pathogenic", "likely_pathogenic", "pathogenic_low_penetrance", "likely_pathogenic_low_penetrance"];

// ---------------------------------------------------------------- 3. gene regions

const GFF: &str = "##gff-version 3
##sequence-region NC_000017.11 1 83257441
##sequence-region NC_000013.11 1 114364328
##sequence-region NT_187600.1 1 200000
##sequence-region NC_000001.11 1 2000
NC_000017.11\tRefSeq\tregion\t1\t83257441\t.\t+\t.\tID=NC_000017.11:1..83257441;Dbxref=taxon:9606;Name=17;chromosome=17;gbkey=Src;genome=chromosome;mol_type=genomic DNA
NC_000017.11\tBestRefSeq\tgene\t43044295\t43170327\t.\t-\t.\tID=gene-BRCA1;Dbxref=GeneID:672,HGNC:HGNC:1100;Name=BRCA1;gbkey=Gene;gene=BRCA1;gene_biotype=protein_coding;gene_synonym=BRCAI,BRCC1,FANCS
NC_000017.11\tBestRefSeq\tmRNA\t43044295\t43125483\t.\t-\t.\tID=rna-NM_007294.4;Parent=gene-BRCA1;Name=NM_007294.4;gene=BRCA1
NC_000017.11\tBestRefSeq\texon\t43124017\t43125483\t.\t-\t.\tID=exon-1;Parent=rna-NM_007294.4;Name=BRCA1;gene=BRCA1
NC_000017.11\tBestRefSeq\tgene\t7668421\t7687490\t.\t-\t.\tID=gene-TP53;Dbxref=GeneID:7157;Name=TP53;gene=TP53;gene_biotype=protein_coding
NC_000017.11\tRefSeq\tregion\t5000\t6000\t.\t+\t.\tID=id-repeat;gbkey=repeat_region
NC_000013.11\tRefSeq\tregion\t1\t114364328\t.\t+\t.\tID=NC_000013.11:1..114364328;Name=13;chromosome=13;gbkey=Src;genome=chromosome
NC_000013.11\tBestRefSeq\tgene\t32315077\t32400268\t.\t+\t.\tID=gene-BRCA2;Dbxref=GeneID:675;Name=BRCA2;gene=BRCA2;gene_biotype=protein_coding
NT_187600.1\tRefSeq\tregion\t1\t200000\t.\t+\t.\tID=NT_187600.1:1..200000;Name=17;chromosome=17;gbkey=Src;genome=genomic;map=unlocalized
NT_187600.1\tBestRefSeq\tgene\t100\t5000\t.\t-\t.\tID=gene-BRCA1-2;Dbxref=GeneID:672;Name=BRCA1;gene=BRCA1;gene_biotype=protein_coding
NT_187601.1\tRefSeq\tregion\t1\t90000\t.\t+\t.\tID=NT_187601.1:1..90000;Name=Un;gbkey=Src;genome=genomic
NT_187601.1\tGnomon\tgene\t500\t900\t.\t+\t.\tID=gene-ALTONLY;Name=ALTONLY;gene=ALTONLY;gene_biotype=lncRNA
NC_000023.11\tRefSeq\tregion\t1\t156040895\t.\t+\t.\tID=NC_000023.11:1..156040895;Name=X;chromosome=X;gbkey=Src;genome=chromosome
NC_000023.11\tBestRefSeq\tgene\t624344\t659411\t.\t+\t.\tID=gene-SHOX;Dbxref=GeneID:6473;Name=SHOX;gene=SHOX;gene_biotype=protein_coding
NC_000024.10\tRefSeq\tregion\t1\t57227415\t.\t+\t.\tID=NC_000024.10:1..57227415;Name=Y;chromosome=Y;gbkey=Src;genome=chromosome
NC_000024.10\tBestRefSeq\tgene\t624344\t659411\t.\t+\t.\tID=gene-SHOX-2;Dbxref=GeneID:6473;Name=SHOX;gene=SHOX;gene_biotype=protein_coding
NC_012920.1\tRefSeq\tregion\t1\t16569\t.\t+\t.\tID=NC_012920.1:1..16569;Is_circular=true;Name=MT;gbkey=Src;genome=mitochondrion
NC_012920.1\tRefSeq\tgene\t5904\t7445\t.\t+\t.\tID=gene-COX1;Dbxref=GeneID:4512;Name=COX1;gene=COX1;gene_biotype=protein_coding;gene_synonym=MTCO1
NC_000001.11\tBestRefSeq\tgene\t50\t300\t.\t+\t.\tID=gene-EDGE1;Name=EDGE1;gene=EDGE1;gene_biotype=protein_coding
NC_000001.11\tBestRefSeq\tpseudogene\t1900\t1990\t.\t-\t.\tID=gene-PSG1;Name=PSG1;gene=PSG1;gene_biotype=pseudogene;pseudo=true
NC_000001.11\tBestRefSeq\tncRNA_gene\t1000\t1100\t.\t+\t.\tID=gene-NCX;Name=NCX%2C1;gene=NCX%2C1
###
NC_000002.12\tBestRefSeq\tgene\t1\t10\t.\t+\t.\tID=gene-WEIRD;Name=Weird;gene=WEIRD
";

fn testa_genes(t: &mut T) {
    let gff = t.p("anotacao.gff3.gz");
    grava(&gff, &GFF.lines().map(str::to_string).collect::<Vec<_>>());
    let pedidos: Vec<String> =
        ["BRCA1", "brca2", "TP53", "SHOX", "ALTONLY", "COX1", "EDGE1", "PSG1", "NCX,1", "BRCC1", "NOPE", "BRCA1", "weird"].map(String::from).to_vec();
    let bed = t.p("painel_genes.bed");
    let Some(r) = t.ok_ou(gene_regions(&gff, &pedidos, 100, None, Some(&bed)), "gene regions") else { return };
    let esperado = [
        "1\t0\t400\tEDGE1\t0\t+",
        "1\t899\t1200\tNCX,1\t0\t+",
        "1\t1799\t2000\tPSG1\t0\t-",
        "2\t0\t110\tWeird\t0\t+",
        "13\t32314976\t32400368\tBRCA2\t0\t+",
        "17\t7668320\t7687590\tTP53\t0\t-",
        "17\t43044194\t43170427\tBRCA1\t0\t-",
        "X\t624243\t659511\tSHOX\t0\t+",
        "Y\t624243\t659511\tSHOX\t0\t+",
        "MT\t5803\t7545\tCOX1\t0\t+",
        "NT_187601.1\t399\t1000\tALTONLY\t0\t+",
    ]
    .map(String::from)
    .to_vec();
    t.lista(&le(&bed), &esperado, "gene regions: BED worked out by hand (padding, clipping, accession mapping, sort)");
    let got: Vec<String> = r.regions.iter().map(|x| format!("{}\t{}\t{}\t{}\t0\t{}", x.chrom, x.start, x.end, x.gene, x.strand)).collect();
    t.lista(&got, &esperado, "gene regions: returned regions equal the BED");
    t.eq(&r.written, &Some(11), "gene regions: lines written");
    t.eq(&r.not_found, &vec!["BRCC1".to_string(), "NOPE".to_string()], "gene regions: not found");
    t.eq(&r.synonyms, &vec![("BRCC1".to_string(), vec!["BRCA1".to_string()])], "gene regions: synonym pointed out");
    t.eq(
        &r.several_sequences,
        &vec![
            ("BRCA1".to_string(), vec!["17 (NC_000017.11)".to_string(), "NT_187600.1".to_string()]),
            ("SHOX".to_string(), vec!["X (NC_000023.11)".to_string(), "Y (NC_000024.10)".to_string()]),
        ],
        "gene regions: genes on several sequences",
    );
    t.eq(&r.sequences_mapped, &7, "gene regions: sequences mapped by name");
    let brca1 = r.regions.iter().find(|x| x.gene == "BRCA1").unwrap();
    t.check(
        brca1.sequence == "NC_000017.11" && brca1.primary && brca1.gene_id.as_deref() == Some("672") && brca1.biotype.as_deref() == Some("protein_coding") && brca1.query == "BRCA1",
        "gene regions: BRCA1 details",
    );
    let alt = r.regions.iter().find(|x| x.gene == "ALTONLY").unwrap();
    t.check(!alt.primary && alt.chrom == "NT_187601.1", "gene regions: a gene only on an unplaced sequence is kept, not primary");
    t.check(r.regions.iter().any(|x| x.gene == "PSG1" && x.feature == "pseudogene") && r.regions.iter().any(|x| x.feature == "ncRNA_gene"), "gene regions: pseudogene and ncRNA_gene features");
    let w = r.regions.iter().find(|x| x.gene == "Weird").unwrap();
    t.check(w.query == "weird", "gene regions: query kept as asked");
    t.check(r.json().contains("\"not_found\":[\"BRCC1\",\"NOPE\"]"), "gene regions: JSON");

    // chr naming; padding 0.
    if let Some(r) = t.ok_ou(gene_regions(&gff, &["BRCA1".into(), "COX1".into(), "ALTONLY".into(), "EDGE1".into()], 0, Some("chr"), None), "gene regions, chr naming") {
        let got: Vec<String> = r.regions.iter().map(|x| format!("{}:{}-{}", x.chrom, x.start, x.end)).collect();
        t.eq(&got, &["chr1:49-300", "chr17:43044294-43170327", "chrM:5903-7445", "NT_187601.1:499-900"].map(String::from).to_vec(), "gene regions: chr naming, no padding");
        t.eq(&r.written, &None, "gene regions: nothing written without output");
    }

    // GTF (GENCODE-like, chr names) and Ensembl-like GFF3.
    let gtf = t.p("genes.gtf.gz");
    grava(
        &gtf,
        &[
            "#!genome-build GRCh38",
            "chr7\tHAVANA\tgene\t117480025\t117668665\t.\t+\t.\tgene_id \"ENSG00000001626.16\"; gene_type \"protein_coding\"; gene_name \"CFTR\"; level 2;",
            "chr7\tHAVANA\ttranscript\t117480025\t117668665\t.\t+\t.\tgene_id \"ENSG00000001626.16\"; transcript_id \"ENST1\"; gene_name \"CFTR\";",
            "chr11\tHAVANA\tgene\t5225464\t5227071\t.\t-\t.\tgene_id \"ENSG00000244734.4\"; gene_type \"protein_coding\"; gene_name \"HBB\"; level 1;",
        ]
        .map(String::from),
    );
    if let Some(r) = t.ok_ou(gene_regions(&gtf, &["cftr".into(), "HBB".into()], 0, Some("no_chr"), None), "gene regions, GTF") {
        let got: Vec<String> = r.regions.iter().map(|x| format!("{}:{}-{} {} {:?} {:?}", x.chrom, x.start, x.end, x.gene, x.gene_id, x.biotype)).collect();
        t.eq(
            &got,
            &[
                "7:117480024-117668665 CFTR Some(\"ENSG00000001626.16\") Some(\"protein_coding\")",
                "11:5225463-5227071 HBB Some(\"ENSG00000244734.4\") Some(\"protein_coding\")",
            ]
            .map(String::from)
            .to_vec(),
            "gene regions: GTF attributes, no_chr naming",
        );
    }
    let ens = t.p("ensembl.gff3");
    grava(&ens, &["##gff-version 3", "7\tensembl_havana\tgene\t117480025\t117668665\t.\t+\t.\tID=gene:ENSG00000001626;Name=CFTR;biotype=protein_coding;gene_id=ENSG00000001626"].map(String::from));
    if let Some(r) = t.ok_ou(gene_regions(&ens, &["CFTR".into()], 10, Some("chr"), None), "gene regions, Ensembl GFF3") {
        t.check(
            r.regions.len() == 1 && r.regions[0].chrom == "chr7" && r.regions[0].start == 117480014 && r.regions[0].gene_id.as_deref() == Some("ENSG00000001626"),
            "gene regions: Ensembl GFF3 (gene_id, biotype), chr naming of a plain chromosome",
        );
    }
    t.erro(gene_regions(&gff, &["BRCA1".into()], 0, Some("ucsc_maybe"), None), &["ucsc_maybe", "no_chr"], "gene regions: unknown naming");
    t.erro(gene_regions(&gff, &[" ".into()], 0, None, None), &["No gene"], "gene regions: no gene");
    let exons = t.p("so_exons.gff3");
    grava(&exons, &["##gff-version 3", "1\tx\texon\t1\t10\t.\t+\t.\tName=BRCA1"].map(String::from));
    t.erro(gene_regions(&exons, &["BRCA1".into()], 0, None, None), &["no gene features"], "gene regions: a GFF without gene features");
    t.erro(gene_regions(&gff, &["BRCA1".into()], 0, None, Some(&gff)), &["destroy"], "gene regions: output over the input");
}

// ---------------------------------------------------------------- 4. lookup

struct Pessoa {
    regs: Vec<PRec>,
    col: Option<usize>,
    chr: bool,
    bed: Option<Vec<(String, u64, u64)>>,
}

fn status_esperado(v: &(String, u64, String, String), p: &Pessoa, com_linha: bool) -> String {
    let vfim = v.1 + v.2.len() as u64;
    let mut melhor: Option<(u8, bool, String)> = None;
    let mut outros: Vec<String> = Vec::new();
    for r in p.regs.iter().filter(|r| !r.outra && r.crom == v.0) {
        if !(v.1 < r.fim() && r.pos < vfim) {
            continue;
        }
        let nome = if p.chr { com_chr(&r.crom) } else { r.crom.clone() };
        let alelos: Vec<Option<(u64, String, String)>> = r.alts.iter().map(|a| apara(r.pos, &r.ref_, a)).collect();
        let exato = alelos.iter().position(|x| x.as_ref().is_some_and(|x| x.0 == v.1 && x.1 == v.2 && x.2 == v.3)).map(|i| i as u32 + 1);
        let gt = r.gt(p.col);
        let gtx = gt.unwrap_or("");
        let alts_gt: Vec<u32> = gt.map(|g| g.split(|c| c == '/' || c == '|').filter_map(|x| x.parse::<u32>().ok()).filter(|x| *x > 0).collect()).unwrap_or_default();
        let carrega = |k: u32| gt.is_none() || alts_gt.contains(&k);
        let mut aqui: Vec<String> = Vec::new();
        for (i, a) in alelos.iter().enumerate() {
            let k = i as u32 + 1;
            if Some(k) == exato || !carrega(k) || r.alts[i] == "." {
                continue;
            }
            let rot = if gtx.is_empty() { "listed" } else { gtx };
            aqui.push(match a {
                Some((p_, rr, aa)) => format!("{nome}:{p_}:{rr}:{aa} ({rot})"),
                None => format!("{nome}:{}:{}:{} ({rot})", r.pos, r.ref_, r.alts[i]),
            });
        }
        let estado: Option<(u8, &str, String)> = match exato {
            Some(k) => {
                let z = if gt.is_none() { "unknown" } else { zig_esperada(gt, k) };
                Some(match z {
                    "het" | "hom" | "hemizygous" | "unknown" => (5, "present", z.to_string()),
                    "other_alt" => (4, "other_allele", z.to_string()),
                    "no_call" => (3, "no_call", z.to_string()),
                    _ => (2, "reference", z.to_string()),
                })
            }
            None => {
                let todos_ref = gt.is_some_and(|g| g.split(|c| c == '/' || c == '|').all(|x| x == "0"));
                let falta = gt.is_some_and(|g| g.split(|c| c == '/' || c == '|').any(|x| x == "."));
                if !aqui.is_empty() {
                    Some((4, "other_allele", "other_alt".to_string()))
                } else if gt.is_none() {
                    None
                } else if falta && alts_gt.is_empty() {
                    Some((3, "no_call", "no_call".to_string()))
                } else if todos_ref {
                    Some((2, "reference", "hom_ref".to_string()))
                } else {
                    None
                }
            }
        };
        for o in aqui {
            if !outros.contains(&o) {
                outros.push(o);
            }
        }
        if let Some((pr, st, z)) = estado {
            let troca = match &melhor {
                None => true,
                Some((mp, me, _)) => pr > *mp || (pr == *mp && exato.is_some() && !*me),
            };
            if troca {
                let dp = r.valor(p.col.filter(|_| gt.is_some()), "DP").unwrap_or_default();
                let ad_alt = match exato {
                    Some(k) => r.valor(p.col.filter(|_| gt.is_some()), "AD").and_then(|ad| ad.split(',').nth(k as usize).map(str::to_string)).unwrap_or_default(),
                    None => String::new(),
                };
                let linha = if com_linha { r.linha.to_string() } else { "-".into() };
                melhor = Some((pr, exato.is_some(), format!("{st} {gtx} {z} dp={dp} ad={ad_alt} {nome} line={linha}")));
            }
        }
    }
    let base = match melhor {
        Some((_, _, s)) => s,
        None => match &p.bed {
            None => "absent".to_string(),
            Some(b) => {
                let (a, f) = (v.1 - 1, v.1 - 1 + v.2.len() as u64);
                if b.iter().any(|(n, x, y)| *n == v.0 && *x <= a && f <= *y) {
                    "absent_callable".to_string()
                } else {
                    "not_callable".to_string()
                }
            }
        },
    };
    let cb = p.bed.as_ref().map(|b| {
        let (a, f) = (v.1 - 1, v.1 - 1 + v.2.len() as u64);
        b.iter().any(|(n, x, y)| *n == v.0 && *x <= a && f <= *y)
    });
    format!("{base} | {} | {cb:?}", outros.join(", "))
}

fn status_obtido(c: &Call, com_linha: bool) -> String {
    let base = if c.chrom.is_none() {
        c.status.clone()
    } else {
        format!(
            "{} {} {} dp={} ad={} {} line={}",
            c.status,
            c.genotype,
            c.zygosity,
            c.depth.map(|d| d.to_string()).unwrap_or_default(),
            c.alt_depth.map(|d| d.to_string()).unwrap_or_default(),
            c.chrom.clone().unwrap_or_default(),
            if com_linha { c.line.map(|l| l.to_string()).unwrap_or_default() } else { "-".into() }
        )
    };
    format!("{base} | {} | {:?}", c.other_alleles.join(", "), c.in_callable)
}

/// Random callable regions covering most of each chromosome, merged.
fn regioes_chamaveis(rng: &mut Rng, g: &[Crom]) -> Vec<(String, u64, u64)> {
    let mut v = Vec::new();
    for c in g {
        let mut p = 0u64;
        while p < c.seq.len() as u64 {
            let on = 100 + rng.below(2000);
            v.push((c.nome.to_string(), p, (p + on).min(c.seq.len() as u64)));
            p += on + 1 + rng.below(700);
        }
    }
    v
}

fn grava_bed(p: &str, r: &[(String, u64, u64)], chr: bool) {
    let mut l = vec!["# callable".to_string()];
    for (n, a, b) in r {
        // Split each interval in two touching pieces: the module must merge them.
        let m = (a + b) / 2;
        let nome = if chr { com_chr(n) } else { n.clone() };
        l.push(format!("{nome}\t{a}\t{m}"));
        l.push(format!("{nome}\t{m}\t{b}"));
    }
    grava(p, &l);
}

fn testa_consulta(t: &mut T, rng: &mut Rng, g: &[Crom], cv: &[CvRec], cv_path: &str) {
    // A family file (three samples, chr names) with a reference block
    // (gVCF style, END=) on chromosome 3, and an aunt (one sample, no chr).
    let mut fam = gera_registros(rng, g, cv, 30, 3, &[("3", 15_000, 15_300)]);
    let gts: Vec<String> = vec!["0/0".into(), "0/0".into(), "./.".into()];
    fam.push(PRec { crom: "3".into(), outra: false, pos: 15_000, ref_: "A".into(), alts: vec!["<NON_REF>".into()], info: "END=15200".into(), campos: gts.clone(), gts, format: "GT".into(), linha: 0 });
    fam.sort_by(|a, b| ordem(&a.crom).cmp(&ordem(&b.crom)).then(a.pos.cmp(&b.pos)));
    let tia = gera_registros(rng, g, cv, 30, 1, &[]);
    let ordem_nat = ["1", "2", "3", "10", "22", "X", "Y", "MT"];
    let (fp, tp) = (t.p("familia.vcf"), t.p("tia.vcf"));
    let fam_w = grava_pessoa(&fp, &fam, &ordem_nat, true, &["KID", "MOM", "DAD"], "GRCh38");
    let tia_w = grava_pessoa(&tp, &tia, &ordem_nat, false, &["AUNT"], "GRCh38");
    let (bk, bm, ba) = (regioes_chamaveis(rng, g), regioes_chamaveis(rng, g), regioes_chamaveis(rng, g));
    let (bkp, bmp, bap) = (t.p("kid.bed"), t.p("mom.bed.gz"), t.p("aunt.bed"));
    grava_bed(&bkp, &bk, true);
    grava_bed(&bmp, &bm, false);
    grava_bed(&bap, &ba, false);
    let pessoas = vec![
        Person { vcf: fp.clone(), sample: Some("KID".into()), callable: Some(bkp.clone()) },
        Person { vcf: fp.clone(), sample: Some("MOM".into()), callable: Some(bmp.clone()) },
        Person { vcf: fp.clone(), sample: Some("DAD".into()), callable: None },
        Person { vcf: tp.clone(), sample: None, callable: Some(bap.clone()) },
    ];
    let esp_pessoas = [
        Pessoa { regs: fam_w.clone(), col: Some(0), chr: true, bed: Some(bk.clone()) },
        Pessoa { regs: fam_w.clone(), col: Some(1), chr: true, bed: Some(bm.clone()) },
        Pessoa { regs: fam_w.clone(), col: Some(2), chr: true, bed: None },
        Pessoa { regs: tia_w.clone(), col: Some(0), chr: false, bed: Some(ba.clone()) },
    ];

    // The queries and the variants each must give.
    let mut consultas: Vec<String> = Vec::new();
    let mut esp_var: Vec<(String, String, u64, String, String)> = Vec::new(); // query, chrom, pos, ref, alt
    let comparaveis = |regs: &[PRec]| -> Vec<(usize, usize)> {
        regs.iter().enumerate().flat_map(|(i, r)| (0..r.alts.len()).filter(move |&k| apara(r.pos, &r.ref_, &r.alts[k]).is_some()).map(move |k| (i, k))).collect()
    };
    let cf = comparaveis(&fam_w);
    for j in 0..14 {
        let (i, k) = *rng.pick(&cf);
        let r = &fam_w[i];
        let (p, rr, aa) = apara(r.pos, &r.ref_, &r.alts[k]).unwrap();
        let (txt, crom) = if j % 3 == 0 {
            (format!("{}-{}-{}-{}", r.crom, r.pos, r.ref_, r.alts[k]), r.crom.clone())
        } else {
            (format!("{}:{}:{}:{}", com_chr(&r.crom), r.pos, r.ref_.to_lowercase(), r.alts[k]), com_chr(&r.crom))
        };
        consultas.push(txt.clone());
        esp_var.push((txt, crom, p, rr, aa));
    }
    let ct = comparaveis(&tia_w);
    for _ in 0..4 {
        let (i, k) = *rng.pick(&ct);
        let r = &tia_w[i];
        let (p, rr, aa) = apara(r.pos, &r.ref_, &r.alts[k]).unwrap();
        let txt = format!("{}:{}:{}:{}", r.crom, r.pos, r.ref_, r.alts[k]);
        consultas.push(txt.clone());
        esp_var.push((txt, r.crom.clone(), p, rr, aa));
    }
    for (c, p) in [("3", 15_050u64), ("3", 15_120), ("1", 7_001), ("2", 12_345), ("22", 19_000), ("MT", 333), ("X", 4_444), ("Y", 7_000)] {
        let crom = g.iter().find(|x| x.nome == c).unwrap();
        let b = crom.seq[p as usize - 1];
        let txt = format!("{c}:{p}:{}:{}", b as char, outra_base(rng, b) as char);
        let pr = apara(p, &txt.split(':').nth(2).unwrap().to_string(), txt.split(':').nth(3).unwrap()).unwrap();
        consultas.push(txt.clone());
        esp_var.push((txt, c.to_string(), pr.0, pr.1, pr.2));
    }
    // rsIDs of ClinVar records the family carries, and two VariationIDs.
    let chaves_fam: HashSet<(String, u64, String, String)> = fam_w
        .iter()
        .flat_map(|r| r.alts.iter().filter_map(move |a| apara(r.pos, &r.ref_, a).map(|(p, x, y)| (r.crom.clone(), p, x, y))))
        .collect();
    let com_rs: Vec<&CvRec> = cv.iter().filter(|c| c.rs.is_some() && c.primario() && chaves_fam.contains(&(c.crom.clone(), c.pos, c.ref_.clone(), c.alt.clone()))).collect();
    let mut rs_usados: Vec<String> = Vec::new();
    for _ in 0..5 {
        let c = rng.pick(&com_rs);
        let rs = c.rs.clone().unwrap();
        if rs_usados.contains(&rs) {
            continue;
        }
        rs_usados.push(rs.clone());
        let txt = format!("rs{rs}");
        consultas.push(txt.clone());
        let mut vistos: Vec<(String, u64, String, String)> = Vec::new();
        for x in cv.iter().filter(|x| x.rs.as_deref() == Some(rs.as_str())) {
            let (p, rr, aa) = apara(x.pos, &x.ref_, &x.alt).unwrap();
            let k = (x.crom.clone(), p, rr.clone(), aa.clone());
            if vistos.contains(&k) {
                continue;
            }
            vistos.push(k);
            esp_var.push((txt.clone(), x.crom.clone(), p, rr, aa));
        }
    }
    let primarios: Vec<&CvRec> = cv.iter().filter(|c| c.primario()).collect();
    let vids: Vec<&CvRec> = vec![com_rs[0], *rng.pick(&primarios)];
    for (j, c) in vids.iter().enumerate() {
        let txt = if j == 0 { c.vid.clone() } else { format!("VCV{:0>9}.3", c.vid) };
        consultas.push(txt.clone());
        let (p, rr, aa) = apara(c.pos, &c.ref_, &c.alt).unwrap();
        esp_var.push((txt, c.crom.clone(), p, rr, aa));
    }
    consultas.push("rs1".into());

    // ClinVar annotations of each expected variant.
    let mut cv_por_chave: HashMap<(String, u64, String, String), &CvRec> = HashMap::new();
    for c in cv {
        let (p, rr, aa) = apara(c.pos, &c.ref_, &c.alt).unwrap();
        cv_por_chave.insert((canon(&c.crom), p, rr, aa), c);
    }
    let anot = |crom: &str, p: u64, rr: &str, aa: &str| -> String {
        match cv_por_chave.get(&(canon(crom), p, rr.to_string(), aa.to_string())) {
            Some(c) => format!("{} {:?} {} {}", c.vid, c.rs.as_ref().map(|r| format!("rs{r}")), c.significancia(), c.estrelas()),
            None => "-".into(),
        }
    };
    let esp_vars: Vec<String> = esp_var.iter().map(|(q, c, p, r, a)| format!("{q} => {c}:{p}:{r}:{a} [{}]", anot(c, *p, r, a))).collect();
    let linha_var = |v: &LookupVariant| -> String {
        let a = if v.in_clinvar {
            format!("{} {:?} {} {}", v.variation_id.clone().unwrap_or_default(), v.rsid, v.significance.clone().unwrap_or_default(), v.stars.unwrap_or(9))
        } else {
            "-".into()
        };
        format!("{} => {}:{}:{}:{} [{a}]", v.query, v.chrom, v.pos, v.ref_allele, v.alt_allele)
    };

    let esperado_calls = |vars: &[LookupVariant], com_linha: bool| -> Vec<String> {
        let mut v = Vec::new();
        for (vi, x) in vars.iter().enumerate() {
            let chave = (canon(&x.chrom), x.pos, x.ref_allele.clone(), x.alt_allele.clone());
            for (pi, p) in esp_pessoas.iter().enumerate() {
                v.push(format!("{vi}/{pi} {}", status_esperado(&chave, p, com_linha)));
            }
        }
        v
    };
    let obtido_calls = |r: &Lookup, com_linha: bool| -> Vec<String> {
        r.calls.iter().map(|c| format!("{}/{} {}", c.variant, c.person, status_obtido(c, com_linha))).collect()
    };

    // 1. Plain files, read whole; ClinVar read whole (rsIDs).
    let Some(r) = t.ok_ou(variant_lookup(&consultas, &pessoas, Some(cv_path)), "lookup") else { return };
    t.lista(&r.variants.iter().map(linha_var).collect::<Vec<_>>(), &esp_vars, "lookup: variants resolved (coordinates, rsIDs, VariationIDs) and annotated");
    t.eq(&r.unresolved, &vec!["rs1".to_string()], "lookup: unknown rsID reported");
    t.eq(&r.people, &["KID", "MOM", "DAD", "AUNT"].map(String::from).to_vec(), "lookup: people labels");
    t.eq(&r.indexed, &vec![false; 4], "lookup: files read whole");
    let esp = esperado_calls(&r.variants, true);
    t.lista(&obtido_calls(&r, true), &esp, &format!("lookup: {} calls equal the brute force", esp.len()));
    let estados: HashSet<&str> = r.calls.iter().map(|c| c.status.as_str()).collect();
    for s in ["present", "reference", "no_call", "other_allele", "absent", "absent_callable", "not_callable"] {
        t.check(estados.contains(s), &format!("lookup: the data exercises status {s}"));
    }
    t.check(r.json().contains("\"calls\":[") && r.json().contains("absent_callable"), "lookup: JSON and text");

    // 2. BGZF files with tabix indexes: the same calls (without line numbers).
    let (fz, tz) = (t.p("familia.vcf.gz"), t.p("tia.vcf.gz"));
    grava(&fz, &le(&fp));
    grava(&tz, &le(&tp));
    for p in [&fz, &tz] {
        t.ok_ou(indice::indexar(p, None), "lookup: indexing the person files");
    }
    let mut pz = pessoas.clone();
    pz[0].vcf = fz.clone();
    pz[1].vcf = fz.clone();
    pz[2].vcf = fz.clone();
    pz[3].vcf = tz.clone();
    if let Some(r2) = t.ok_ou(variant_lookup(&consultas, &pz, Some(cv_path)), "lookup through the index") {
        t.eq(&r2.indexed, &vec![true; 4], "lookup: files read through the index");
        t.lista(&obtido_calls(&r2, false), &esperado_calls(&r2.variants, false), "lookup through the index: calls equal the brute force");
        t.check(r2.calls.iter().all(|c| c.line.is_none()), "lookup through the index: no line numbers");
    }

    // 3. Coordinates only: ClinVar through its index gives the same annotations.
    let so_coord: Vec<String> = consultas.iter().filter(|q| q.contains(':') || q.contains('-')).cloned().collect();
    let cvi = t.p("clinvar_indexado.vcf.gz");
    std::fs::copy(cv_path, &cvi).unwrap();
    t.ok_ou(indice::indexar(&cvi, None), "lookup: indexing ClinVar");
    let a = t.ok_ou(variant_lookup(&so_coord, &pessoas, Some(cv_path)), "lookup, coordinates, ClinVar read whole");
    let b = t.ok_ou(variant_lookup(&so_coord, &pessoas, Some(&cvi)), "lookup, coordinates, ClinVar through its index");
    if let (Some(a), Some(b)) = (a, b) {
        t.lista(&b.variants.iter().map(linha_var).collect::<Vec<_>>(), &a.variants.iter().map(linha_var).collect::<Vec<_>>(), "lookup: ClinVar through its index annotates the same");
        t.check(a.variants.iter().filter(|v| v.in_clinvar).count() >= 10, "lookup: coordinate queries found in ClinVar");
    }
    // Without ClinVar: coordinates only, nothing annotated.
    if let Some(r) = t.ok_ou(variant_lookup(&so_coord, &pessoas[..1], None), "lookup without ClinVar") {
        t.check(r.variants.iter().all(|v| !v.in_clinvar) && r.calls.len() == so_coord.len(), "lookup without ClinVar: one call per variant");
    }

    // Errors.
    t.erro(variant_lookup(&["rs123".into()], &pessoas, None), &["rs123", "clinvar="], "lookup error: rsID without ClinVar");
    t.erro(variant_lookup(&["17:abc:T:C".into()], &pessoas, None), &["17:abc:T:C", "not understood", "chrom:pos:ref:alt"], "lookup error: malformed variant");
    t.erro(variant_lookup(&["1:100:A:A".into()], &pessoas, None), &["same"], "lookup error: REF equal to ALT");
    t.erro(variant_lookup(&["1:100:A:G".into()], &[], None), &["No person"], "lookup error: no person");
    t.erro(variant_lookup(&[], &pessoas, None), &["No variant"], "lookup error: no variant");
    let sem = vec![Person { vcf: fp.clone(), sample: Some("UNCLE".into()), callable: None }];
    t.erro(variant_lookup(&["1:100:A:G".into()], &sem, None), &["UNCLE", "KID, MOM, DAD"], "lookup error: unknown sample");
    t.erro(variant_lookup(&["rs5".into()], &pessoas, Some(&fp)), &["does not look like the ClinVar VCF"], "lookup error: ClinVar that is not ClinVar");
}

// ---------------------------------------------------------------- 5. pileup

#[derive(Clone, Copy, Debug, PartialEq)]
enum ObsSim {
    Base(u8, u8),
    /// Index of the query base after the deletion (its quality decides).
    Del(usize),
    Pulo,
}

struct ReadSim {
    nome: String,
    flag: u16,
    crom: usize,
    pos: u64,
    mapq: u8,
    cigar: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
    seq_estrela: bool,
    qual_estrela: bool,
    obs: BTreeMap<u64, ObsSim>,
    ins_depois: HashMap<u64, String>,
    del_depois: HashMap<u64, u32>,
}

impl ReadSim {
    fn qualidade(&self, o: ObsSim) -> u8 {
        match o {
            ObsSim::Base(_, q) => q,
            ObsSim::Del(i) => {
                if self.seq_estrela {
                    0
                } else if self.qual_estrela {
                    255
                } else {
                    self.qual.get(i).map(|q| q - 33).unwrap_or(0)
                }
            }
            ObsSim::Pulo => 0,
        }
    }
    fn sam(&self, refs: &[(&str, usize)]) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}",
            self.nome,
            self.flag,
            refs[self.crom].0,
            self.pos,
            self.mapq,
            self.cigar,
            if self.seq_estrela { "*".to_string() } else { String::from_utf8(self.seq.clone()).unwrap() },
            if self.qual_estrela || self.seq_estrela { "*".to_string() } else { String::from_utf8(self.qual.clone()).unwrap() }
        )
    }
}

fn gera_read(rng: &mut Rng, nome: String, crom: usize, refseq: &[u8]) -> ReadSim {
    let inicio = 1 + rng.below(refseq.len() as u64 - 400);
    let mut ops: Vec<(u32, u8)> = Vec::new();
    let mut seq: Vec<u8> = Vec::new();
    let mut qual: Vec<u8> = Vec::new();
    let mut obs: BTreeMap<u64, ObsSim> = BTreeMap::new();
    let mut r = inicio;
    let base_aleat = |rng: &mut Rng| BASES[rng.below(4) as usize];
    let add_query = |rng: &mut Rng, seq: &mut Vec<u8>, qual: &mut Vec<u8>, b: u8| {
        seq.push(b);
        qual.push(33 + rng.below(42) as u8);
    };
    if rng.chance(10) {
        ops.push((1 + rng.below(5) as u32, b'H'));
    }
    if rng.chance(25) {
        let n = 1 + rng.below(10) as u32;
        ops.push((n, b'S'));
        for _ in 0..n {
            let b = base_aleat(rng);
            add_query(rng, &mut seq, &mut qual, b);
        }
    }
    if rng.chance(5) {
        let n = 1 + rng.below(3) as u32;
        ops.push((n, b'I'));
        for _ in 0..n {
            let b = base_aleat(rng);
            add_query(rng, &mut seq, &mut qual, b);
        }
    }
    let alvo = 30 + rng.below(150);
    let mut anterior = b'S';
    let mut usado = 0u64;
    loop {
        let op = if anterior == b'M' || anterior == b'=' || anterior == b'X' {
            if usado >= alvo {
                break;
            }
            match rng.below(100) {
                0..=9 => b'I',
                10..=19 => b'D',
                20..=23 => b'N',
                24..=27 => b'p', // P then I
                28..=30 => b'q', // P then M
                31..=33 => b'd', // D D
                34..=43 => b'=',
                44..=50 => b'X',
                _ => b'M',
            }
        } else if anterior == b'D' {
            if rng.chance(15) {
                b'I'
            } else {
                b'M'
            }
        } else {
            b'M'
        };
        let mut faz = |op: u8, rng: &mut Rng, ops: &mut Vec<(u32, u8)>, r: &mut u64, usado: &mut u64| {
            let n = match op {
                b'M' | b'=' | b'X' => 1 + rng.below(40) as u32,
                b'I' => 1 + rng.below(4) as u32,
                b'D' => 1 + rng.below(5) as u32,
                b'N' => 5 + rng.below(60) as u32,
                _ => 1 + rng.below(2) as u32,
            };
            ops.push((n, op));
            for _ in 0..n {
                match op {
                    b'M' | b'=' | b'X' => {
                        let rb = refseq[*r as usize - 1];
                        let b = if rng.chance(1) { b'N' } else if rng.chance(8) { outra_base(rng, rb) } else { rb };
                        let b = if rng.chance(3) { b.to_ascii_lowercase() } else { b };
                        seq.push(b);
                        let q = 33 + rng.below(42) as u8;
                        qual.push(q);
                        obs.insert(*r, ObsSim::Base(b.to_ascii_uppercase(), q - 33));
                        *r += 1;
                        *usado += 1;
                    }
                    b'I' => {
                        seq.push(BASES[rng.below(4) as usize]);
                        qual.push(33 + rng.below(42) as u8);
                    }
                    b'D' => {
                        obs.insert(*r, ObsSim::Del(seq.len()));
                        *r += 1;
                        *usado += 1;
                    }
                    b'N' => {
                        obs.insert(*r, ObsSim::Pulo);
                        *r += 1;
                        *usado += 1;
                    }
                    _ => {}
                }
            }
        };
        match op {
            b'p' => {
                faz(b'P', rng, &mut ops, &mut r, &mut usado);
                faz(b'I', rng, &mut ops, &mut r, &mut usado);
                anterior = b'I';
            }
            b'q' => {
                faz(b'P', rng, &mut ops, &mut r, &mut usado);
                faz(b'M', rng, &mut ops, &mut r, &mut usado);
                anterior = b'M';
            }
            b'd' => {
                faz(b'D', rng, &mut ops, &mut r, &mut usado);
                faz(b'D', rng, &mut ops, &mut r, &mut usado);
                anterior = b'D';
            }
            op => {
                faz(op, rng, &mut ops, &mut r, &mut usado);
                anterior = op;
            }
        }
    }
    if rng.chance(20) {
        let n = 1 + rng.below(8) as u32;
        ops.push((n, b'S'));
        for _ in 0..n {
            seq.push(base_aleat(rng));
            qual.push(33 + rng.below(42) as u8);
        }
    }
    if rng.chance(8) {
        ops.push((1 + rng.below(4) as u32, b'H'));
    }
    // Indels after each reference-consuming operation, by htslib's rule.
    let mut ins_depois = HashMap::new();
    let mut del_depois = HashMap::new();
    let (mut rr, mut q) = (inicio, 0usize);
    for i in 0..ops.len() {
        let (n, op) = ops[i];
        match op {
            b'M' | b'=' | b'X' => {
                rr += u64::from(n);
                q += n as usize;
            }
            b'D' | b'N' => rr += u64::from(n),
            b'I' | b'S' => {
                q += n as usize;
                continue;
            }
            _ => continue,
        }
        if op == b'N' {
            continue;
        }
        let fim = rr - 1;
        let Some(&(n2, op2)) = ops.get(i + 1) else { continue };
        if op2 == b'D' && op != b'D' {
            let mut total = n2;
            for &(l, o) in &ops[i + 2..] {
                if o == b'D' {
                    total += l;
                } else {
                    break;
                }
            }
            del_depois.insert(fim, total);
        } else if op2 == b'I' {
            ins_depois.insert(fim, String::from_utf8(seq[q..q + n2 as usize].to_ascii_uppercase()).unwrap());
        } else if op2 == b'P' {
            let mut s = String::new();
            let mut qq = q;
            for &(l, o) in &ops[i + 2..] {
                match o {
                    b'I' => {
                        s.push_str(std::str::from_utf8(&seq[qq..qq + l as usize].to_ascii_uppercase()).unwrap());
                        qq += l as usize;
                    }
                    b'D' | b'M' | b'N' | b'=' | b'X' => break,
                    _ => {}
                }
            }
            if !s.is_empty() {
                ins_depois.insert(fim, s);
            }
        }
    }
    let cigar: String = ops.iter().map(|(n, o)| format!("{n}{}", *o as char)).collect();
    let flag = match rng.below(100) {
        0..=4 => 1024,
        5..=7 => 256,
        8..=10 => 2048,
        11..=12 => 512,
        _ => *rng.pick(&[0u16, 16, 99, 147, 83, 163]),
    };
    let mut rd = ReadSim {
        nome,
        flag,
        crom,
        pos: inicio,
        mapq: if rng.chance(5) { 255 } else { rng.below(61) as u8 },
        cigar,
        seq,
        qual,
        seq_estrela: false,
        qual_estrela: false,
        obs,
        ins_depois,
        del_depois,
    };
    if rng.chance(5) {
        rd.qual_estrela = true;
        for o in rd.obs.values_mut() {
            if let ObsSim::Base(b, _) = o {
                *o = ObsSim::Base(*b, 255);
            }
        }
    } else if rng.chance(1) {
        rd.seq_estrela = true;
        for o in rd.obs.values_mut() {
            if let ObsSim::Base(_, _) = o {
                *o = ObsSim::Base(b'N', 0);
            }
        }
        for s in rd.ins_depois.values_mut() {
            *s = "N".repeat(s.len());
        }
    }
    rd
}

fn pileup_esperado(reads: &[ReadSim], crom: usize, x: u64, o: &PileupOptions) -> String {
    let (mut lidos, mut fflag, mut fmapq, mut fbq, mut pulos) = (0, 0, 0, 0, 0);
    let mut bases: BTreeMap<u8, (usize, usize)> = BTreeMap::new();
    let mut ins: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    let mut dels: BTreeMap<u32, (usize, usize)> = BTreeMap::new();
    for r in reads.iter().filter(|r| r.crom == crom) {
        let Some(&ob) = r.obs.get(&x) else { continue };
        lidos += 1;
        if r.flag & o.exclude_flags != 0 {
            fflag += 1;
            continue;
        }
        if r.mapq < o.min_mapq {
            fmapq += 1;
            continue;
        }
        if ob == ObsSim::Pulo {
            pulos += 1;
            continue;
        }
        if r.qualidade(ob) < o.min_baseq {
            fbq += 1;
            continue;
        }
        let f = usize::from(r.flag & 16 == 0);
        let b = match ob {
            ObsSim::Base(b, _) if b"ACGT".contains(&b) => b,
            ObsSim::Base(..) => b'N',
            _ => b'*',
        };
        let e = bases.entry(b).or_default();
        e.0 += 1;
        e.1 += f;
        if let Some(s) = r.ins_depois.get(&x) {
            let e = ins.entry(s.clone()).or_default();
            e.0 += 1;
            e.1 += f;
        }
        if let Some(d) = r.del_depois.get(&x) {
            let e = dels.entry(*d).or_default();
            e.0 += 1;
            e.1 += f;
        }
    }
    let prof: usize = bases.values().map(|v| v.0).sum();
    let mut alelos: Vec<String> = b"ACGTN*".iter().map(|b| {
        let (n, f) = bases.get(b).copied().unwrap_or((0, 0));
        format!("{}:{n}:{f}:{}", *b as char, n - f)
    }).collect();
    let mut iv: Vec<(&String, &(usize, usize))> = ins.iter().collect();
    iv.sort_by(|a, b| b.1 .0.cmp(&a.1 .0).then_with(|| a.0.cmp(b.0)));
    alelos.extend(iv.iter().map(|(s, (n, f))| format!("+{s}:{n}:{f}:{}", n - f)));
    let mut dv: Vec<(&u32, &(usize, usize))> = dels.iter().collect();
    dv.sort_by(|a, b| b.1 .0.cmp(&a.1 .0).then_with(|| a.0.cmp(b.0)));
    alelos.extend(dv.iter().map(|(d, (n, f))| format!("-{d}:{n}:{f}:{}", n - f)));
    let ni: usize = ins.values().map(|v| v.0).sum();
    let nd: usize = dels.values().map(|v| v.0).sum();
    format!("{x} reads={lidos} flag={fflag} mapq={fmapq} bq={fbq} skips={pulos} depth={prof} ins={ni} dels={nd} {}", alelos.join(" "))
}

fn pileup_obtido(p: &Pileup) -> String {
    let al: Vec<String> = p.alleles.iter().map(|a| format!("{}:{}:{}:{}", a.allele, a.count, a.forward, a.reverse)).collect();
    let soma = p.a + p.c + p.g + p.t + p.n + p.deletions;
    let frac_ok = p.alleles.iter().all(|a| (a.fraction - if p.depth == 0 { 0.0 } else { a.count as f64 / p.depth as f64 }).abs() < 1e-12);
    format!(
        "{} reads={} flag={} mapq={} bq={} skips={} depth={} ins={} dels={} {}{}",
        p.pos,
        p.reads,
        p.filtered_flag,
        p.filtered_mapq,
        p.filtered_baseq,
        p.ref_skips,
        p.depth,
        p.insertions_after,
        p.deletions_after,
        al.join(" "),
        if soma == p.depth && frac_ok { "" } else { " (counts do not add up)" }
    )
}

fn testa_pileup(t: &mut T, rng: &mut Rng) {
    let refs: Vec<(&str, usize)> = vec![("chr7", 3000), ("chr8", 2000)];
    let seqs: Vec<Vec<u8>> = refs.iter().map(|(_, l)| (0..*l).map(|_| BASES[rng.below(4) as usize]).collect()).collect();
    let mut reads: Vec<ReadSim> = (0..2600).map(|i| {
        let c = if i % 5 == 0 { 1 } else { 0 };
        gera_read(rng, format!("r{i}"), c, &seqs[c])
    }).collect();
    reads.sort_by_key(|r| (r.crom, r.pos));
    let mut linhas = vec!["@HD\tVN:1.6\tSO:coordinate".to_string()];
    for (n, l) in &refs {
        linhas.push(format!("@SQ\tSN:{n}\tLN:{l}"));
    }
    let n_cab = linhas.len();
    linhas.extend(reads.iter().map(|r| r.sam(&refs)));
    // Two reads that never span anything: unmapped but placed, CIGAR `*`.
    linhas.insert(n_cab + 10, format!("u1\t4\tchr7\t{}\t0\t*\t*\t0\t0\tACGT\tIIII", reads[9].pos));
    let sam = t.p("leituras.sam");
    let bam_p = t.p("leituras.bam");
    let bam_sem = t.p("sem_indice.bam");
    grava(&sam, &linhas);
    grava(&bam_p, &linhas);
    grava(&bam_sem, &linhas);
    t.ok_ou(indice::indexar(&bam_p, None), "pileup: indexing the BAM");
    t.check(t.dir.join("leituras.bam.bai").exists(), "pileup: .bai written");

    let conjuntos = [
        PileupOptions::default(),
        PileupOptions { min_mapq: 0, min_baseq: 0, exclude_flags: 0, alt: None },
        PileupOptions { min_mapq: 20, min_baseq: 25, exclude_flags: 0x400, alt: None },
    ];
    let mut total = 0usize;
    let mut tipos = [0usize; 4]; // positions with insertions, deletions after, deletions over, skips
    for (oi, o) in conjuntos.iter().enumerate() {
        for (c, (nome, len)) in refs.iter().enumerate() {
            let (mut got, mut esp) = (Vec::new(), Vec::new());
            let mut x = 1 + oi as u64;
            while x <= *len as u64 {
                match pileup(&bam_p, nome, x, o) {
                    Ok(p) => {
                        if oi == 1 {
                            tipos[0] += usize::from(p.insertions_after > 0);
                            tipos[1] += usize::from(p.deletions_after > 0);
                            tipos[2] += usize::from(p.deletions > 0);
                            tipos[3] += usize::from(p.ref_skips > 0);
                        }
                        t.check(p.indexed, "pileup: read through the index");
                        got.push(pileup_obtido(&p));
                    }
                    Err(e) => got.push(format!("error {e}")),
                }
                esp.push(pileup_esperado(&reads, c, x, o));
                x += 3;
                total += 1;
            }
            t.lista(&got, &esp, &format!("pileup, BAM + index, {nome}, options {oi}: {} positions equal the brute force", esp.len()));
        }
    }
    t.check(tipos.iter().all(|&n| n > 20), &format!("pileup: positions with insertions, deletions after, deletions over, skips: {tipos:?}"));
    t.check(total > 3000, "pileup: number of positions compared");
    // SAM text and a BAM without index give the same.
    for (rotulo, caminho) in [("SAM text", &sam), ("BAM without index", &bam_sem)] {
        let (mut got, mut esp) = (Vec::new(), Vec::new());
        for x in (5..3000u64).step_by(41) {
            let o = &conjuntos[(x % 3) as usize];
            match pileup(caminho, "chr7", x, o) {
                Ok(p) => {
                    t.check(!p.indexed, "pileup: read whole");
                    got.push(pileup_obtido(&p));
                }
                Err(e) => got.push(format!("error {e}")),
            }
            esp.push(pileup_esperado(&reads, 0, x, o));
        }
        t.lista(&got, &esp, &format!("pileup, {rotulo}: equal to the brute force"));
    }
    // Naming: "7" finds chr7.
    if let (Some(a), Some(b)) = (t.ok_ou(pileup(&bam_p, "7", 1500, &conjuntos[0]), "pileup 7"), t.ok_ou(pileup(&sam, "7", 1500, &conjuntos[0]), "pileup SAM 7")) {
        t.check(a.chrom == "chr7" && b.chrom == "chr7" && pileup_obtido(&a) == pileup_obtido(&b), "pileup: chromosome named without chr");
    }

    // ALT and REF fractions, in each accepted form, against the brute force.
    let mut testados = [0usize; 5];
    for x in 1..3000u64 {
        let esp = |pred: &dyn Fn(&ReadSim, ObsSim) -> bool| -> (usize, usize) {
            let o = &conjuntos[1];
            let mut prof = 0;
            let mut n = 0;
            for r in reads.iter().filter(|r| r.crom == 0) {
                let Some(&ob) = r.obs.get(&x) else { continue };
                if r.flag & o.exclude_flags != 0 || r.mapq < o.min_mapq || ob == ObsSim::Pulo || r.qualidade(ob) < o.min_baseq {
                    continue;
                }
                prof += 1;
                n += usize::from(pred(r, ob));
            }
            (n, prof)
        };
        let rb = seqs[0][x as usize - 1] as char;
        let mut casos: Vec<(String, (usize, usize), usize)> = Vec::new();
        if let Some(s) = reads.iter().filter(|r| r.crom == 0).find_map(|r| r.ins_depois.get(&x)) {
            let s = s.clone();
            casos.push((format!("+{s}"), esp(&|r, _| r.ins_depois.get(&x) == Some(&s)), 0));
            casos.push((format!("{rb}>{rb}{}", s.to_lowercase()), esp(&|r, _| r.ins_depois.get(&x) == Some(&s)), 1));
        }
        if let Some(d) = reads.iter().filter(|r| r.crom == 0).find_map(|r| r.del_depois.get(&x)).copied() {
            casos.push((format!("-{d}"), esp(&|r, _| r.del_depois.get(&x) == Some(&d)), 2));
            let del: String = (0..d as usize).map(|i| seqs[0].get(x as usize + i).copied().unwrap_or(b'A') as char).collect();
            casos.push((format!("{rb}{del}>{rb}"), esp(&|r, _| r.del_depois.get(&x) == Some(&d)), 2));
        }
        if x % 50 == 0 {
            casos.push(("*".into(), esp(&|_, o| matches!(o, ObsSim::Del(_))), 3));
            casos.push(("t".into(), esp(&|_, o| o == ObsSim::Base(b'T', 0) || matches!(o, ObsSim::Base(b'T', _))), 4));
            let outra = if rb == 'G' { b'C' } else { b'G' };
            casos.push((format!("{rb}>{}", outra as char), esp(&|_, o| matches!(o, ObsSim::Base(b, _) if b == outra)), 4));
        }
        for (alt, (n, prof), tipo) in casos {
            if testados[tipo] > 60 {
                continue;
            }
            let o = PileupOptions { alt: Some(alt.clone()), ..conjuntos[1].clone() };
            match pileup(&bam_p, "chr7", x, &o) {
                Ok(p) => {
                    let f = if prof == 0 { 0.0 } else { n as f64 / prof as f64 };
                    let ok = p.alt_count == Some(n) && p.depth == prof && p.alt_fraction.is_some_and(|y| (y - f).abs() < 1e-12);
                    if !ok {
                        println!("   alt {alt} at {x}: got {:?} of {}, expected {n} of {prof}", p.alt_count, p.depth);
                    }
                    t.check(ok, &format!("pileup: ALT {alt} at chr7:{x}"));
                    if alt.contains('>') {
                        let r_esp = esp(&|_, o| matches!(o, ObsSim::Base(b, _) if b == rb as u8)).0;
                        t.check(p.ref_count == Some(r_esp) && p.ref_allele.as_deref() == Some(&rb.to_string()[..]), &format!("pileup: REF of {alt} at chr7:{x}"));
                    }
                    testados[tipo] += 1;
                }
                Err(e) => {
                    t.check(false, &format!("pileup: ALT {alt} at chr7:{x}: {e}"));
                }
            }
        }
    }
    t.check(testados.iter().all(|&n| n > 5), &format!("pileup: ALT forms tested {testados:?}"));
    if let Some(p) = t.ok_ou(pileup(&bam_p, "chr7", 1000, &PileupOptions { alt: Some("A>T".into()), ..Default::default() }), "pileup JSON") {
        t.check(p.json().contains("\"alleles\":[{\"allele\":\"A\"") && p.json().contains("\"ref\":\"A\""), "pileup: JSON");
    }

    // Errors.
    t.erro(pileup(&bam_p, "chr9", 10, &conjuntos[0]), &["chr9", "chr7, chr8", "leituras.bam"], "pileup error: unknown sequence (BAM)");
    t.erro(pileup(&sam, "chr9", 10, &conjuntos[0]), &["chr9", "chr7, chr8"], "pileup error: unknown sequence (SAM)");
    t.erro(pileup(&bam_p, "chr7", 0, &conjuntos[0]), &["1-based"], "pileup error: position 0");
    for (a, agulha) in [("XYZ", "not understood"), ("AC>GT", "MNP"), ("+12", "not understood"), ("-0", "not understood")] {
        let o = PileupOptions { alt: Some(a.into()), ..Default::default() };
        t.erro(pileup(&bam_p, "chr7", 10, &o), &[agulha], &format!("pileup error: ALT {a}"));
    }
    let ruim = t.p("ruim.sam");
    grava(&ruim, &["@SQ\tSN:chr7\tLN:3000", "q1\t0\tchr7\t5\t60\t10Q\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII"].map(String::from));
    t.erro(pileup(&ruim, "chr7", 6, &conjuntos[0]), &["q1", "10Q", "ruim.sam"], "pileup error: invalid CIGAR");
}

// ---------------------------------------------------------------- real data

fn real(a: &[String]) {
    if a.len() < 8 {
        eprintln!("usage: clinico real <clinvar.vcf.gz> <refseq.gff.gz> <son.vcf.gz> <son.bed> <father.vcf.gz> <father.bed> <mother.vcf.gz> <mother.bed>");
        std::process::exit(2);
    }
    let out = std::env::temp_dir().join(format!("mmgen_clinico_real_{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let o = |n: &str| out.join(n).to_string_lossy().to_string();
    let (cv, gff) = (&a[0], &a[1]);
    let filho = Person { vcf: a[2].clone(), sample: None, callable: Some(a[3].clone()) };
    let pai = Person { vcf: a[4].clone(), sample: None, callable: Some(a[5].clone()) };
    let mae = Person { vcf: a[6].clone(), sample: None, callable: Some(a[7].clone()) };

    let t0 = std::time::Instant::now();
    let opt = FindingsOptions { classes: vec!["plp".into()], min_stars: 1, output: Some(o("findings_plp.tsv")), ..Default::default() };
    let r = clinvar_findings(&filho, cv, &opt).unwrap_or_else(|e| panic!("{e}"));
    println!("=== ClinVar findings, P/LP with >= 1 star ({:.1} s)\n{}", t0.elapsed().as_secs_f64(), r.resumo());
    println!("by significance (all matched alleles):");
    for (s, n) in &r.by_significance {
        println!("  {n:>8}  {s}");
    }
    println!(
        "person records {}, alleles checked {}, not carried {}, other sequences {}, not comparable {}, skipped {}; ClinVar records {}, other sequences {}",
        r.person_records, r.person_alleles, r.person_not_carried, r.person_other_sequences, r.person_not_comparable, r.person_skipped, r.clinvar_records, r.clinvar_other_sequences
    );

    let painel: Vec<String> =
        ["BRCA1", "BRCA2", "PALB2", "TP53", "MLH1", "MSH2", "MSH6", "PMS2", "LDLR", "APOB", "PCSK9", "CFTR", "HBB"].map(String::from).to_vec();
    let t1 = std::time::Instant::now();
    let bed = o("panel.bed");
    let g = gene_regions(gff, &painel, 0, Some("chr"), Some(&bed)).unwrap_or_else(|e| panic!("{e}"));
    println!("\n=== Panel from the GFF ({:.1} s)\n{}", t1.elapsed().as_secs_f64(), g.resumo());

    let t2 = std::time::Instant::now();
    let opt = FindingsOptions { regions: Some(bed.clone()), ..Default::default() };
    let rp = clinvar_findings(&filho, cv, &opt).unwrap_or_else(|e| panic!("{e}"));
    println!("\n=== ClinVar matches inside the panel, all classes ({:.1} s)", t2.elapsed().as_secs_f64());
    println!("  listed {}; by category among them:", rp.findings.len());
    let mut cats: BTreeMap<String, usize> = BTreeMap::new();
    for f in &rp.findings {
        *cats.entry(f.category.clone()).or_default() += 1;
    }
    for (c, n) in cats {
        println!("    {c}: {n}");
    }
    for f in rp.findings.iter().filter(|f| f.category != "benign_or_likely_benign").take(30) {
        println!("    {}:{} {}>{} {} {} {} ({} stars)", f.chrom, f.pos, f.ref_allele, f.alt_allele, f.genes.join(","), f.zygosity, f.significance, f.stars);
    }

    // HG002's P/LP findings and the panel matches that are not benign, by
    // coordinates (ClinVar read through its index), then one by rsID (ClinVar
    // read whole).
    let mut consultas: Vec<String> = r
        .findings
        .iter()
        .chain(rp.findings.iter().filter(|f| f.category != "benign_or_likely_benign"))
        .map(|f| format!("{}:{}:{}:{}", f.chrom, f.pos, f.ref_allele, f.alt_allele))
        .collect();
    consultas.dedup();
    let trio = [filho.clone(), pai.clone(), mae.clone()];
    let t3 = std::time::Instant::now();
    let l = variant_lookup(&consultas, &trio, Some(cv)).unwrap_or_else(|e| panic!("{e}"));
    println!("\n=== Lookup of {} of HG002's matches in the trio ({:.1} s, indexed {:?})\n{}", consultas.len(), t3.elapsed().as_secs_f64(), l.indexed, l.resumo());
    if let Some(rs) = r.findings.iter().find_map(|f| f.rsid.clone()) {
        let t4 = std::time::Instant::now();
        let l = variant_lookup(&[rs.clone()], &trio, Some(cv)).unwrap_or_else(|e| panic!("{e}"));
        println!("=== Lookup of {rs} through ClinVar ({:.1} s)\n{}", t4.elapsed().as_secs_f64(), l.resumo());
    }
    let _ = std::fs::remove_dir_all(&out);
}

// ---------------------------------------------------------------- main

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(String::as_str) == Some("real") {
        real(&args[2..]);
        return;
    }
    let dir = std::env::temp_dir().join(format!("mmgen_clinico_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut t = T { ok: 0, falhas: 0, dir: dir.clone() };
    let inicio = std::time::Instant::now();

    testa_estrelas(&mut t);
    testa_classes(&mut t);
    let mut rng = Rng(20261114);
    let g = genoma(&mut rng);
    let mut cv = gera_clinvar(&mut rng, &g);
    let cv_path = t.p("clinvar.vcf.gz");
    grava_clinvar(&cv_path, &mut cv);
    println!("synthetic ClinVar: {} records", cv.len());
    testa_achados(&mut t, &mut rng, &g, &cv, &cv_path);
    testa_genes(&mut t);
    testa_consulta(&mut t, &mut rng, &g, &cv, &cv_path);
    testa_pileup(&mut t, &mut rng);

    let _ = std::fs::remove_dir_all(&dir);
    println!("{} of {} checks passed ({:.1} s)", t.ok, t.ok + t.falhas, inicio.elapsed().as_secs_f64());
    if t.falhas > 0 {
        std::process::exit(1);
    }
}
