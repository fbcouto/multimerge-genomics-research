#![allow(dead_code)]
// Trials of src/ancestralidade.rs: frequency panels, supervised admixture and
// polygenic scores.
//
// Every result is compared with something computed here on its own:
//  - panel: the generator keeps each record as plain strings, and the panel
//    expected from them is counted by brute force (one sample at a time, one
//    allele at a time), following the rules documented in the module;
//  - ancestry: individuals are simulated with known proportions under the
//    admixture model; the EM must recover them within a tolerance that
//    shrinks with more sites, must reach the maximum found by a grid search
//    of a log-likelihood written here, and its bootstrap standard error must
//    agree with the Fisher-information one;
//  - polygenic score: every variant is planted with a known fate (genotyped,
//    homozygous reference, imputed, ambiguous, flipped, mismatched...) and
//    the score is the direct sum; the group expectations are summed here and
//    the percentile is integrated numerically.
//
// Usage: cargo run --release --bin ancestralidade
//        cargo run --release --bin ancestralidade -- real   (the real-data run;
//        see `real` at the end for the paths)
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
#[path = "../../../src/referencia.rs"]
#[allow(dead_code)]
mod referencia;
#[path = "../../../src/ancestralidade.rs"]
mod ancestralidade;

use ancestralidade::*;
use std::collections::{HashMap, HashSet};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::time::Instant;

// ---------------------------------------------------------------------
//  CHECKS
// ---------------------------------------------------------------------

struct T {
    dir: PathBuf,
    total: usize,
    ok: usize,
    falhas: Vec<String>,
}

impl T {
    fn p(&self, n: &str) -> String {
        self.dir.join(n).to_string_lossy().to_string()
    }

    /// Writes text lines; a `.gz` name gets BGZF.
    fn write(&self, n: &str, linhas: &[String]) -> String {
        let p = self.p(n);
        let mut w = externo::Escritor::cria(&p).unwrap();
        for l in linhas {
            w.linha(l).unwrap();
        }
        w.termina().unwrap();
        p
    }

    fn check(&mut self, what: &str, cond: bool, detalhe: impl FnOnce() -> String) {
        self.total += 1;
        if cond {
            self.ok += 1;
        } else {
            let m = format!("FAILED {what}: {}", detalhe());
            println!("{m}");
            self.falhas.push(m);
        }
    }

    fn near(&mut self, what: &str, got: f64, want: f64, tol: f64) {
        self.check(what, (got - want).abs() <= tol || (got.is_nan() && want.is_nan()), || format!("got {got}, expected {want} (tolerance {tol})"));
    }

    fn eq<X: PartialEq + std::fmt::Debug>(&mut self, what: &str, got: X, want: X) {
        let ok = got == want;
        self.check(what, ok, || format!("got {got:?}, expected {want:?}"));
    }

    /// Must fail with a message holding every needle.
    fn fails<X: std::fmt::Debug>(&mut self, what: &str, r: Result<X, String>, needles: &[&str]) {
        match r {
            Ok(x) => self.check(what, false, || format!("should have failed, returned {x:?}")),
            Err(e) => {
                let falta: Vec<&&str> = needles.iter().filter(|n| !e.contains(**n)).collect();
                self.check(what, falta.is_empty(), || format!("the message should mention {falta:?}; it is: {e}"));
            }
        }
    }
}

fn le_texto(p: &str) -> Vec<String> {
    let mut r = externo::abre_leitura(p).unwrap();
    let mut s = String::new();
    r.read_to_string(&mut s).unwrap();
    s.lines().map(str::to_string).collect()
}

// ---------------------------------------------------------------------
//  RANDOM NUMBERS
// ---------------------------------------------------------------------

struct Rng(u64);

impl Rng {
    fn u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn f(&mut self) -> f64 {
        (self.u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn below(&mut self, n: usize) -> usize {
        (self.u64() % n as u64) as usize
    }
    fn chance(&mut self, p: f64) -> bool {
        self.f() < p
    }
    fn normal(&mut self) -> f64 {
        let u = self.f().max(1e-300);
        let v = self.f();
        (-2.0 * u.ln()).sqrt() * (2.0 * std::f64::consts::PI * v).cos()
    }
    fn base(&mut self) -> char {
        ['A', 'C', 'G', 'T'][self.below(4)]
    }
    /// Two different bases that are not complementary (not A/T, not C/G).
    fn par_seguro(&mut self) -> (char, char) {
        loop {
            let (a, b) = (self.base(), self.base());
            if a != b && comp(a) != b {
                return (a, b);
            }
        }
    }
}

fn comp(c: char) -> char {
    match c {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        x => x,
    }
}

// ---------------------------------------------------------------------
//  1. PANEL FROM GENOTYPES, AGAINST A BRUTE-FORCE COUNT
// ---------------------------------------------------------------------

#[derive(Clone)]
struct Rec {
    crom: String,
    pos: u64,
    ref_: String,
    alt: String,
    filtro: String,
    formato: String,
    /// GT of each sample, as written.
    gts: Vec<String>,
    info: String,
}

fn linha_vcf(r: &Rec) -> String {
    let mut l = format!("{}\t{}\t.\t{}\t{}\t50\t{}\t{}", r.crom, r.pos, r.ref_, r.alt, r.filtro, if r.info.is_empty() { "." } else { &r.info });
    if !r.gts.is_empty() {
        l.push('\t');
        l.push_str(&r.formato);
        for g in &r.gts {
            l.push('\t');
            l.push_str(&match r.formato.as_str() {
                "GT" => g.clone(),
                "GT:DP" => format!("{g}:12"),
                "DP:GT" => format!("9:{g}"),
                f => panic!("format {f}"),
            });
        }
    }
    l
}

fn escreve_vcf(t: &T, nome: &str, amostras: &[String], recs: &[Rec], extra_cab: &[String]) -> String {
    let mut l = vec!["##fileformat=VCFv4.2".to_string()];
    l.extend(extra_cab.iter().cloned());
    let mut c = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string();
    if !amostras.is_empty() {
        c.push_str("\tFORMAT\t");
        c.push_str(&amostras.join("\t"));
    }
    l.push(c);
    l.extend(recs.iter().map(linha_vcf));
    let p = t.p(nome);
    let mut w = externo::Escritor::cria(&p).unwrap();
    for x in &l {
        w.linha(x).unwrap();
    }
    w.termina().unwrap();
    p
}

/// A GT string for a diploid (or now and then haploid / missing) call.
fn gt_aleatorio(r: &mut Rng, p_alt: f64, n_alts: usize) -> String {
    let x = r.f();
    if x < 0.01 {
        return "./.".to_string();
    }
    if x < 0.015 {
        return "0/.".to_string();
    }
    if x < 0.02 {
        return ".".to_string();
    }
    let alelo = |r: &mut Rng| -> String {
        if r.chance(p_alt) {
            (1 + r.below(n_alts)).to_string()
        } else {
            "0".to_string()
        }
    };
    if x < 0.03 {
        return alelo(r);
    }
    let sep = if r.chance(0.5) { "|" } else { "/" };
    format!("{}{sep}{}", alelo(r), alelo(r))
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum TipoT {
    Snv,
    Indel,
    Outro,
    Multi,
}

fn tipo_t(r: &str, a: &str) -> TipoT {
    if a.contains(',') {
        return TipoT::Multi;
    }
    let ok = |s: &str| !s.is_empty() && s.chars().all(|c| "ACGT".contains(c.to_ascii_uppercase()));
    if !ok(r) || !ok(a) {
        TipoT::Outro
    } else if r.len() == 1 && a.len() == 1 {
        if r.eq_ignore_ascii_case(a) {
            TipoT::Outro
        } else {
            TipoT::Snv
        }
    } else if r.len() != a.len() {
        TipoT::Indel
    } else {
        TipoT::Outro
    }
}

#[derive(Default, Debug, PartialEq, Clone)]
struct Contas {
    written: usize,
    not_snv: usize,
    multiallelic: usize,
    filtered: usize,
    no_data: usize,
    below_maf: usize,
    thinned: usize,
}

/// Brute force: the rows the panel must have, as (chrom, pos, ref, alt,
/// per-group (alt, an)).
#[allow(clippy::type_complexity)]
fn painel_esperado(
    recs: &[Rec],
    grupo_de_amostra: &[Option<usize>],
    k: usize,
    o: &PanelOptions,
) -> (Vec<(String, u64, String, String, Vec<(u64, u64)>)>, Contas) {
    let mut c = Contas::default();
    let mut linhas = Vec::new();
    let mut i = 0;
    let mut ultimo: HashMap<String, u64> = HashMap::new();
    while i < recs.len() {
        let mut j = i;
        while j < recs.len() && recs[j].crom == recs[i].crom && recs[j].pos == recs[i].pos {
            j += 1;
        }
        let mut cand = Vec::new();
        for r in &recs[i..j] {
            match tipo_t(&r.ref_, &r.alt) {
                TipoT::Multi => c.multiallelic += 1,
                TipoT::Outro => c.not_snv += 1,
                TipoT::Indel if !o.keep_indels => c.not_snv += 1,
                _ if o.pass_only && r.filtro != "PASS" && r.filtro != "." => c.filtered += 1,
                _ => cand.push(r),
            }
        }
        i = j;
        if cand.len() > 1 {
            c.multiallelic += cand.len();
            continue;
        }
        let Some(r) = cand.pop() else { continue };
        let mut cont = vec![(0u64, 0u64); k];
        for (s, g) in r.gts.iter().enumerate() {
            let Some(grupo) = grupo_de_amostra[s] else { continue };
            // `gts` holds the bare GT; `linha_vcf` adds the other subfields.
            let gt = g.as_str();
            if !gt.chars().all(|ch| ch.is_ascii_digit() || ".|/".contains(ch)) {
                continue;
            }
            for a in gt.split(['/', '|']) {
                if a == "0" {
                    cont[grupo].1 += 1;
                } else if a == "1" {
                    cont[grupo].1 += 1;
                    cont[grupo].0 += 1;
                }
            }
        }
        if cont.iter().any(|x| x.1 == 0) {
            c.no_data += 1;
            continue;
        }
        if !cont.iter().any(|&(x, n)| {
            let p = x as f64 / n as f64;
            p.min(1.0 - p) >= o.min_maf - 1e-12
        }) {
            c.below_maf += 1;
            continue;
        }
        if let Some(&u) = ultimo.get(&r.crom) {
            if o.thin_bp > 0 && r.pos < u + o.thin_bp {
                c.thinned += 1;
                continue;
            }
        }
        ultimo.insert(r.crom.clone(), r.pos);
        c.written += 1;
        linhas.push((r.crom.clone(), r.pos, r.ref_.to_uppercase(), r.alt.to_uppercase(), cont));
    }
    (linhas, c)
}

/// Simulated reference VCF records over some chromosomes.
fn simula_registros(r: &mut Rng, cromos: &[&str], por_crom: usize, freqs_grupo: &[f64], grupo_de: &[Option<usize>]) -> Vec<Rec> {
    let mut out = Vec::new();
    for c in cromos {
        let mut pos = 1000u64;
        for _ in 0..por_crom {
            if !r.chance(0.08) || out.is_empty() {
                pos += 1 + r.below(400) as u64;
            }
            let x = r.f();
            let (ref_, alt): (String, String) = if x < 0.75 {
                let (a, b) = (r.base(), r.base());
                let b = if a == b { comp(a) } else { b };
                if r.chance(0.05) {
                    (a.to_ascii_lowercase().to_string(), b.to_string())
                } else {
                    (a.to_string(), b.to_string())
                }
            } else if x < 0.83 {
                if r.chance(0.5) {
                    ("A".into(), "AT".into())
                } else {
                    ("AC".into(), "A".into())
                }
            } else if x < 0.86 {
                ("AC".into(), "GT".into())
            } else if x < 0.89 {
                ("C".into(), "G,T".into())
            } else if x < 0.91 {
                ("G".into(), "<DEL>".into())
            } else if x < 0.92 {
                ("N".into(), "A".into())
            } else if x < 0.93 {
                ("A".into(), "A".into())
            } else {
                let (a, b) = r.par_seguro();
                (a.to_string(), b.to_string())
            };
            let filtro = if r.chance(0.05) { "LowQual" } else if r.chance(0.03) { "." } else { "PASS" };
            let formato = ["GT", "GT", "GT", "GT:DP", "DP:GT"][r.below(5)];
            let n_alts = alt.split(',').count();
            // Group frequencies of this site: some sites monomorphic, some rare.
            let base_p = if r.chance(0.1) { 0.0 } else if r.chance(0.2) { r.f() * 0.04 } else { r.f() };
            let ps: Vec<f64> = freqs_grupo.iter().map(|d| (base_p + d * r.normal()).clamp(0.0, 1.0)).collect();
            let gts = grupo_de
                .iter()
                .map(|g| {
                    let p = ps[g.unwrap_or(0) % ps.len()];
                    if r.chance(0.002) {
                        "0/2".to_string()
                    } else if r.chance(0.001) {
                        "0/x".to_string()
                    } else {
                        gt_aleatorio(r, p, n_alts)
                    }
                })
                .collect();
            out.push(Rec { crom: c.to_string(), pos, ref_, alt, filtro: filtro.to_string(), formato: formato.to_string(), gts, info: String::new() });
        }
    }
    out
}

/// Data rows of a panel file: (chrom, pos, ref, alt, AF..., AN...).
fn linhas_painel(p: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let l = le_texto(p);
    let mut cab = Vec::new();
    let mut dados = Vec::new();
    for x in l {
        if x.starts_with("##") {
            continue;
        }
        if cab.is_empty() {
            cab = x.split('\t').map(str::to_string).collect();
        } else {
            dados.push(x.split('\t').map(str::to_string).collect());
        }
    }
    (cab, dados)
}

fn compara_painel(t: &mut T, what: &str, saida: &str, esperado: &[(String, u64, String, String, Vec<(u64, u64)>)], grupos: &[&str]) {
    let (cab, dados) = linhas_painel(saida);
    let mut cab_esp = vec!["chrom".to_string(), "pos".into(), "ref".into(), "alt".into()];
    cab_esp.extend(grupos.iter().map(|g| format!("AF_{g}")));
    cab_esp.extend(grupos.iter().map(|g| format!("AN_{g}")));
    t.eq(&format!("{what}: header line"), cab, cab_esp);
    t.eq(&format!("{what}: number of sites"), dados.len(), esperado.len());
    let k = grupos.len();
    let mut ruins = 0;
    let mut primeiro = String::new();
    for (d, e) in dados.iter().zip(esperado) {
        let mut ok = d[0] == e.0 && d[1] == e.1.to_string() && d[2] == e.2 && d[3] == e.3;
        for g in 0..k {
            let (x, n) = e.4[g];
            let af: f64 = d[4 + g].parse().unwrap_or(f64::NAN);
            ok &= (af - x as f64 / n as f64).abs() < 6e-7;
            ok &= d[4 + k + g] == n.to_string();
        }
        if !ok {
            ruins += 1;
            if primeiro.is_empty() {
                primeiro = format!("got {d:?}, expected {e:?}");
            }
        }
    }
    t.check(&format!("{what}: every site's frequencies and allele counts"), ruins == 0, || format!("{ruins} wrong rows; first: {primeiro}"));
}

fn testa_painel(t: &mut T, r: &mut Rng) {
    let grupos = ["AFR", "EAS", "EUR"];
    let mut amostras: Vec<String> = Vec::new();
    let mut grupo_de: Vec<Option<usize>> = Vec::new();
    let mut tabela = vec!["SampleID\tPopulation\tSuperpopulation".to_string()];
    for (g, nome) in grupos.iter().enumerate() {
        for i in 0..15 {
            let s = format!("{nome}{i:02}");
            tabela.push(format!("{s}\tP{g}\t{nome}"));
            amostras.push(s);
            grupo_de.push(Some(g));
        }
    }
    // Samples without a group, and table samples absent from the VCF.
    for i in 0..3 {
        amostras.push(format!("LOOSE{i}"));
        grupo_de.push(None);
        tabela.push(format!("GHOST{i}\tPX\tEUR"));
    }
    let excluir: Vec<String> = vec!["AFR03".into(), "EUR10".into(), "EAS00".into(), "EAS14".into(), "NOT_THERE".into()];
    let grupo_ef: Vec<Option<usize>> = amostras.iter().zip(&grupo_de).map(|(s, g)| if excluir.contains(s) { None } else { *g }).collect();
    // The sample order in the VCF is shuffled against the table.
    let mut ordem: Vec<usize> = (0..amostras.len()).collect();
    for i in (1..ordem.len()).rev() {
        let j = r.below(i + 1);
        ordem.swap(i, j);
    }
    let amostras_vcf: Vec<String> = ordem.iter().map(|&i| amostras[i].clone()).collect();
    let grupo_vcf: Vec<Option<usize>> = ordem.iter().map(|&i| grupo_de[i]).collect();
    let grupo_ef_vcf: Vec<Option<usize>> = ordem.iter().map(|&i| grupo_ef[i]).collect();

    let recs = simula_registros(r, &["chr1", "chr2", "chr10"], 1500, &[0.0, 0.15, 0.3], &grupo_vcf);
    let tab = t.write("pop.tsv", &tabela);
    // Also the 1000 Genomes layout: separated by spaces.
    let tab_esp = t.write("pop_spaces.txt", &tabela.iter().map(|l| l.replace('\t', " ")).collect::<Vec<_>>());
    let vcf_todo = escreve_vcf(t, "ref_all.vcf.gz", &amostras_vcf, &recs, &[]);
    let por_crom: Vec<String> = ["chr10", "chr2", "chr1"]
        .iter()
        .map(|c| {
            let rs: Vec<Rec> = recs.iter().filter(|x| x.crom == *c).cloned().collect();
            escreve_vcf(t, &format!("ref_{c}.vcf"), &amostras_vcf, &rs, &[])
        })
        .collect();

    let casos = [
        PanelOptions::default(),
        PanelOptions { min_maf: 0.0, ..Default::default() },
        PanelOptions { min_maf: 0.1, thin_bp: 700, ..Default::default() },
        PanelOptions { keep_indels: true, pass_only: false, min_maf: 0.02, thin_bp: 0 },
    ];
    for (ci, o) in casos.iter().enumerate() {
        let (esp, contas) = painel_esperado(&recs, &grupo_ef_vcf, 3, o);
        let what = format!("panel case {ci} (min_maf {}, thin {}, indels {}, pass_only {})", o.min_maf, o.thin_bp, o.keep_indels, o.pass_only);
        let saida = t.p(&format!("panel{ci}.tsv.gz"));
        match build_panel(&[vcf_todo.clone()], &tab, "SampleID", "Superpopulation", &excluir, &saida, o) {
            Err(e) => t.check(&what, false, || e),
            Ok(s) => {
                compara_painel(t, &what, &saida, &esp, &grupos);
                let got = Contas {
                    written: s.written,
                    not_snv: s.not_snv,
                    multiallelic: s.multiallelic,
                    filtered: s.filtered,
                    no_data: s.no_data,
                    below_maf: s.below_maf,
                    thinned: s.thinned,
                };
                t.eq(&format!("{what}: summary counts"), got, contas.clone());
                t.eq(&format!("{what}: records"), s.records, recs.len());
                t.eq(&format!("{what}: samples per group"), s.samples_per_group.clone(), vec![14, 13, 14]);
                t.eq(
                    &format!("{what}: sample bookkeeping"),
                    (s.samples_in_vcf, s.samples_excluded, s.samples_without_group, s.table_samples_not_in_vcf),
                    (48, 4, 3, 3),
                );
                t.eq(&format!("{what}: groups"), s.groups.clone(), vec!["AFR".to_string(), "EAS".into(), "EUR".into()]);
                let por_seq: Vec<(String, usize)> = ["chr1", "chr2", "chr10"]
                    .iter()
                    .map(|c| (c.to_string(), esp.iter().filter(|x| x.0 == *c).count()))
                    .filter(|x| x.1 > 0)
                    .collect();
                t.eq(&format!("{what}: sites per chromosome"), s.sequences.clone(), por_seq);
                // The JSON and the text summary say the same numbers.
                t.check(&format!("{what}: json"), s.json().contains(&format!("\"written\":{}", s.written)), || s.json());
                t.check(&format!("{what}: text summary"), s.resumo().contains(&format!("{} site(s) written", s.written)), || s.resumo());
            }
        }
        if ci == 0 {
            // One file per chromosome, given out of order (chr10, chr2, chr1),
            // and the space-separated table: the same panel, in order 1, 2, 10.
            let saida2 = t.p("panel0_split.tsv");
            match build_panel(&por_crom, &tab_esp, "SampleID", "Superpopulation", &excluir, &saida2, o) {
                Err(e) => t.check("panel from one file per chromosome", false, || e),
                Ok(_) => {
                    let a: Vec<Vec<String>> = linhas_painel(&saida).1;
                    let b: Vec<Vec<String>> = linhas_painel(&saida2).1;
                    t.check("panel from one file per chromosome, given out of order: same rows, natural order", a == b, || format!("{} vs {} rows", a.len(), b.len()));
                }
            }
        }
    }

    // ----- panel from INFO frequencies
    let mut recs_f: Vec<Rec> = Vec::new();
    let mut esp_f: Vec<(String, u64, String, String, Vec<f64>, Vec<String>)> = Vec::new();
    let mut c_f = Contas::default();
    let mut pos = 100u64;
    let mut ultimo_kept: Option<u64> = None;
    let o_f = PanelOptions { min_maf: 0.05, thin_bp: 50, ..Default::default() };
    for i in 0..3000 {
        pos += 1 + r.below(60) as u64;
        let (a, b) = r.par_seguro();
        let multi = i % 97 == 5;
        let alt = if multi { format!("{b},{}", comp(b)) } else { b.to_string() };
        let afs: Vec<Option<f64>> = (0..3).map(|_| if r.chance(0.02) { None } else { Some((r.f() * 1e6).round() / 1e6 * if r.chance(0.3) { 0.04 } else { 1.0 }) }).collect();
        let ans: Vec<u64> = (0..3).map(|_| 100 + r.below(900) as u64).collect();
        let mut info = format!("AF=0.5;DP={}", r.below(100));
        for (g, nome) in ["AFR", "EAS", "EUR"].iter().enumerate() {
            match afs[g] {
                Some(x) => info.push_str(&format!(";{nome}_AF={x}")),
                None => {
                    if r.chance(0.5) {
                        info.push_str(&format!(";{nome}_AF=."))
                    }
                }
            }
            info.push_str(&format!(";AN_{nome}={}", ans[g]));
        }
        recs_f.push(Rec { crom: "22".into(), pos, ref_: a.to_string(), alt, filtro: "PASS".into(), formato: String::new(), gts: Vec::new(), info });
        if multi {
            c_f.multiallelic += 1;
        } else if afs.iter().any(Option::is_none) {
            c_f.no_data += 1;
        } else if !afs.iter().any(|p| p.unwrap().min(1.0 - p.unwrap()) >= 0.05 - 1e-12) {
            c_f.below_maf += 1;
        } else if ultimo_kept.is_some_and(|u| pos < u + 50) {
            c_f.thinned += 1;
        } else {
            ultimo_kept = Some(pos);
            c_f.written += 1;
            esp_f.push(("22".into(), pos, a.to_string(), b.to_string(), afs.iter().map(|x| x.unwrap()).collect(), ans.iter().map(|x| x.to_string()).collect()));
        }
    }
    let cab_info: Vec<String> = ["AF", "DP", "AFR_AF", "EAS_AF", "EUR_AF", "AN_AFR", "AN_EAS", "AN_EUR"]
        .iter()
        .map(|i| format!("##INFO=<ID={i},Number=A,Type=Float,Description=\"x\">"))
        .collect();
    let sites = escreve_vcf(t, "sites.vcf.gz", &[], &recs_f, &cab_info);
    let g3: Vec<String> = vec!["AFR".into(), "EAS".into(), "EUR".into()];
    let afc: Vec<String> = vec!["AFR_AF".into(), "EAS_AF".into(), "EUR_AF".into()];
    let anc: Vec<String> = vec!["AN_AFR".into(), "AN_EAS".into(), "AN_EUR".into()];
    let saida_f = t.p("panel_freq.tsv");
    match panel_from_frequencies(&[sites.clone()], &g3, &afc, &anc, &saida_f, &o_f) {
        Err(e) => t.check("panel from INFO frequencies", false, || e),
        Ok(s) => {
            let (_, dados) = linhas_painel(&saida_f);
            let ok = dados.len() == esp_f.len()
                && dados.iter().zip(&esp_f).all(|(d, e)| {
                    d[0] == e.0
                        && d[1] == e.1.to_string()
                        && d[2] == e.2
                        && d[3] == e.3
                        && (0..3).all(|g| (d[4 + g].parse::<f64>().unwrap() - e.4[g]).abs() < 6e-7 && d[7 + g] == e.5[g])
                });
            t.check("panel from INFO frequencies: every row", ok, || format!("{} rows, expected {}", dados.len(), esp_f.len()));
            let got = Contas { written: s.written, multiallelic: s.multiallelic, no_data: s.no_data, below_maf: s.below_maf, thinned: s.thinned, ..Default::default() };
            t.eq("panel from INFO frequencies: counts", got, c_f.clone());
        }
    }
    // Without the AN fields: the AN columns are `.`.
    let saida_f2 = t.p("panel_freq_noan.tsv");
    match panel_from_frequencies(&[sites.clone()], &g3, &afc, &[], &saida_f2, &o_f) {
        Err(e) => t.check("panel from INFO frequencies without AN", false, || e),
        Ok(_) => {
            let (_, dados) = linhas_painel(&saida_f2);
            t.check("panel from INFO frequencies without AN: AN written as `.`", dados.iter().all(|d| d[7..10] == [".", ".", "."]) && dados.len() == esp_f.len(), || format!("{:?}", dados.first()));
        }
    }
    t.fails(
        "panel from INFO: a field the header does not declare",
        panel_from_frequencies(&[sites.clone()], &g3, &["AFR_AF".into(), "EAS_AF".into(), "SAS_AF".into()], &[], &t.p("x.tsv"), &o_f),
        &["SAS_AF", "not declared", "EUR_AF"],
    );
    t.check("panel from INFO: no output left after the error", !Path::new(&t.p("x.tsv")).exists(), || "x.tsv exists".into());
    let mut ruim = recs_f[..20].to_vec();
    ruim[7].info = "AFR_AF=0.3;EAS_AF=abc;EUR_AF=0.1".into();
    let sites_ruim = escreve_vcf(t, "sites_bad.vcf", &[], &ruim, &[]);
    // ##fileformat, #CHROM, then record 7 is on line 10.
    t.fails("panel from INFO: a value that is not a frequency", panel_from_frequencies(&[sites_ruim], &g3, &afc, &[], &t.p("y.tsv"), &o_f), &["line 10", "EAS_AF", "abc"]);
    t.fails("panel from INFO: mismatched lists", panel_from_frequencies(&[sites.clone()], &g3, &afc[..2], &[], &t.p("y.tsv"), &o_f), &["one INFO frequency field per group"]);

    // ----- errors of the genotype builder
    t.fails(
        "panel: table without the named column",
        build_panel(&[vcf_todo.clone()], &tab, "SampleID", "superpop", &[], &t.p("z.tsv"), &PanelOptions::default()),
        &["superpop", "Superpopulation", "SampleID"],
    );
    t.fails(
        "panel: output is the input",
        build_panel(&[vcf_todo.clone()], &tab, "SampleID", "Superpopulation", &[], &vcf_todo, &PanelOptions::default()),
        &["also one of the input files"],
    );
    t.fails("panel: Parquet output refused", build_panel(&[vcf_todo.clone()], &tab, "SampleID", "Superpopulation", &[], &t.p("z.parquet"), &PanelOptions::default()), &[".tsv"]);
    t.fails(
        "panel: a group left without samples",
        build_panel(&[vcf_todo.clone()], &tab, "SampleID", "Superpopulation", &amostras[15..30].to_vec(), &t.p("z.tsv"), &PanelOptions::default()),
        &["EAS", "no sample of group"],
    );
    t.fails("panel: min_maf out of range", build_panel(&[vcf_todo.clone()], &tab, "SampleID", "Superpopulation", &[], &t.p("z.tsv"), &PanelOptions { min_maf: 0.7, ..Default::default() }), &["min_maf"]);
    // Out of order: a position going back, a chromosome coming back.
    let mut fora = recs[..200].to_vec();
    fora.swap(50, 120);
    let v_fora = escreve_vcf(t, "ref_unsorted.vcf", &amostras_vcf, &fora, &[]);
    t.fails(
        "panel: unsorted VCF",
        build_panel(&[v_fora], &tab, "SampleID", "Superpopulation", &[], &t.p("z.tsv"), &PanelOptions::default()),
        &["sort_genomic", "chrom_order=\"natural\"", "line", "comes after"],
    );
    t.check("panel: no output left after an error", !Path::new(&t.p("z.tsv")).exists(), || "z.tsv was left".into());
    let um = recs.iter().position(|x| x.crom == "chr2").unwrap();
    let mut volta: Vec<Rec> = recs[um - 5..um + 5].to_vec();
    volta.push(Rec { pos: 1, ..recs[0].clone() });
    let v_volta = escreve_vcf(t, "ref_back.vcf", &amostras_vcf, &volta, &[]);
    t.fails("panel: a chromosome coming back", build_panel(&[v_volta], &tab, "SampleID", "Superpopulation", &[], &t.p("z.tsv"), &PanelOptions::default()), &["chr1", "appears again", "sort_genomic"]);
    let sem_amostras = escreve_vcf(t, "nosamples.vcf", &[], &recs_f[..10], &[]);
    t.fails("panel: a sites-only VCF for the genotype builder", build_panel(&[sem_amostras], &tab, "SampleID", "Superpopulation", &[], &t.p("z.tsv"), &PanelOptions::default()), &["panel_from_frequencies"]);
}

// ---------------------------------------------------------------------
//  2. ANCESTRY
// ---------------------------------------------------------------------

/// A panel written directly, with known frequencies.
struct SitioT {
    crom: String,
    pos: u64,
    ref_: char,
    alt: char,
    /// As written in the panel (6 decimals).
    p: Vec<f64>,
    an: u64,
}

fn escreve_painel(t: &T, nome: &str, grupos: &[&str], s: &[SitioT]) -> String {
    let mut l = vec!["##fileformat=MMGEN-panel-1".to_string()];
    let mut c = "chrom\tpos\tref\talt".to_string();
    for g in grupos {
        c.push_str(&format!("\tAF_{g}"));
    }
    for g in grupos {
        c.push_str(&format!("\tAN_{g}"));
    }
    l.push(c);
    for x in s {
        let mut r = format!("{}\t{}\t{}\t{}", x.crom, x.pos, x.ref_, x.alt);
        for p in &x.p {
            r.push_str(&format!("\t{p}"));
        }
        for _ in grupos {
            r.push_str(&format!("\t{}", x.an));
        }
        l.push(r);
    }
    t.write(nome, &l)
}

/// Frequencies of K groups around a common one: a structured population.
fn sitios(r: &mut Rng, cromos: &[&str], por_crom: usize, k: usize, passo: u64, dif: f64, an: u64) -> Vec<SitioT> {
    let mut v = Vec::new();
    for c in cromos {
        let mut pos = 10_000u64;
        for _ in 0..por_crom {
            pos += 1 + r.below(passo as usize * 2) as u64;
            let base = 0.05 + 0.9 * r.f();
            let p: Vec<f64> = (0..k).map(|_| ((base + dif * r.normal()).clamp(0.01, 0.99) * 1e6).round() / 1e6).collect();
            let (a, b) = r.par_seguro();
            v.push(SitioT { crom: c.to_string(), pos, ref_: a, alt: b, p, an });
        }
    }
    v
}

fn adj(p: f64, an: u64) -> f64 {
    (p * an as f64 + 0.5) / (an as f64 + 1.0)
}

/// A genotype drawn under the admixture model: each allele picks a group by q.
fn genotipo(r: &mut Rng, q: &[f64], p: &[f64]) -> u8 {
    let mut x = 0;
    for _ in 0..2 {
        let mut u = r.f();
        let mut g = 0;
        while g + 1 < q.len() && u >= q[g] {
            u -= q[g];
            g += 1;
        }
        if r.chance(p[g]) {
            x += 1;
        }
    }
    x
}

/// The log-likelihood written here: (ALT count, adjusted frequencies) per site.
fn ll_t(q: &[f64], dados: &[(u8, Vec<f64>)]) -> f64 {
    let mut ll = 0.0;
    for (x, p) in dados {
        let f: f64 = q.iter().zip(p).map(|(a, b)| a * b).sum();
        ll += match x {
            0 => 2.0 * (1.0 - f).ln(),
            1 => std::f64::consts::LN_2 + f.ln() + (1.0 - f).ln(),
            _ => 2.0 * f.ln(),
        };
    }
    ll
}

/// Maximum on a grid over the simplex (K = 2 or 3), coarse then fine.
fn grade(dados: &[(u8, Vec<f64>)], k: usize) -> (Vec<f64>, f64) {
    let mut melhor = (vec![1.0 / k as f64; k], f64::NEG_INFINITY);
    let avalia = |q: Vec<f64>, m: &mut (Vec<f64>, f64)| {
        if q.iter().any(|x| *x < -1e-12) {
            return;
        }
        let q: Vec<f64> = q.iter().map(|x| x.max(0.0)).collect();
        let l = ll_t(&q, dados);
        if l > m.1 {
            *m = (q, l);
        }
    };
    let passo = 0.02;
    let n = (1.0 / passo) as usize;
    for i in 0..=n {
        if k == 2 {
            let a = i as f64 * passo;
            avalia(vec![a, 1.0 - a], &mut melhor);
        } else {
            for j in 0..=(n - i) {
                let (a, b) = (i as f64 * passo, j as f64 * passo);
                avalia(vec![a, b, 1.0 - a - b], &mut melhor);
            }
        }
    }
    let centro = melhor.0.clone();
    let fino = 0.0005;
    let m = 60i64;
    for i in -m..=m {
        if k == 2 {
            let a = centro[0] + i as f64 * fino;
            avalia(vec![a, 1.0 - a], &mut melhor);
        } else {
            for j in -m..=m {
                let (a, b) = (centro[0] + i as f64 * fino, centro[1] + j as f64 * fino);
                avalia(vec![a, b, 1.0 - a - b], &mut melhor);
            }
        }
    }
    melhor
}

struct PessoaSim {
    nome: String,
    q: Vec<f64>,
    /// ALT count per site; `None` = missing in the person's file.
    g: Vec<Option<u8>>,
}

/// A person's VCF in the GIAB style: only non-reference genotypes are listed
/// (plus explicit missing calls); `callable` decides which absent sites are
/// homozygous reference. Several samples share the file.
fn escreve_pessoas(t: &T, nome: &str, s: &[SitioT], pessoas: &[PessoaSim], prefixo: &str) -> String {
    let nomes: Vec<String> = pessoas.iter().map(|p| p.nome.clone()).collect();
    let mut recs = Vec::new();
    for (j, x) in s.iter().enumerate() {
        let alguem = pessoas.iter().any(|p| p.g[j] != Some(0));
        if !alguem {
            continue;
        }
        let gts = pessoas
            .iter()
            .map(|p| match p.g[j] {
                None => "./.".to_string(),
                Some(0) => "0/0".to_string(),
                Some(1) => if j % 2 == 0 { "0|1" } else { "1|0" }.to_string(),
                Some(_) => "1/1".to_string(),
            })
            .collect();
        recs.push(Rec {
            crom: format!("{prefixo}{}", x.crom),
            pos: x.pos,
            ref_: x.ref_.to_string(),
            alt: x.alt.to_string(),
            filtro: "PASS".into(),
            formato: "GT:DP".into(),
            gts,
            info: String::new(),
        });
    }
    escreve_vcf(t, nome, &nomes, &recs, &[])
}

fn simula_pessoa(r: &mut Rng, nome: &str, q: &[f64], s: &[SitioT], falta: f64) -> PessoaSim {
    let g = s.iter().map(|x| if r.chance(falta) { None } else { Some(genotipo(r, q, &x.p)) }).collect();
    PessoaSim { nome: nome.to_string(), q: q.to_vec(), g }
}

fn dados_t(p: &PessoaSim, s: &[SitioT]) -> Vec<(u8, Vec<f64>)> {
    p.g.iter().zip(s).filter_map(|(g, x)| g.map(|g| (g, x.p.iter().map(|&f| adj(f, x.an)).collect()))).collect()
}

fn testa_ancestralidade(t: &mut T, r: &mut Rng) {
    let grupos = ["AFR", "EUR", "EAS"];
    let qs: Vec<(&str, Vec<f64>)> = vec![
        ("pureAFR", vec![1.0, 0.0, 0.0]),
        ("pureEUR", vec![0.0, 1.0, 0.0]),
        ("pureEAS", vec![0.0, 0.0, 1.0]),
        ("half", vec![0.5, 0.5, 0.0]),
        ("three", vec![0.6, 0.3, 0.1]),
        ("even", vec![1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]),
        ("mostlyEAS", vec![0.05, 0.15, 0.8]),
    ];
    let mut erros_por_tamanho: Vec<(usize, f64)> = Vec::new();
    for (ti, &por_crom) in [400usize, 4000, 20000].iter().enumerate() {
        let s = sitios(r, &["1", "2", "3"], por_crom, 3, 300, 0.18, 1_000_000);
        let painel = escreve_painel(t, &format!("anc_panel{ti}.tsv"), &grupos, &s);
        let pessoas: Vec<PessoaSim> = qs.iter().map(|(n, q)| simula_pessoa(r, n, q, &s, 0.03)).collect();
        // chr-prefixed person against an unprefixed panel.
        let vcf = escreve_pessoas(t, &format!("anc_people{ti}.vcf.gz"), &s, &pessoas, "chr");
        let persons: Vec<Person> = pessoas.iter().map(|p| Person { vcf: vcf.clone(), sample: Some(p.nome.clone()), callable: None }).collect();
        let o = AncestryOptions { bootstrap: if ti == 1 { 60 } else { 0 }, block_bp: 200_000, ..Default::default() };
        let tempo = Instant::now();
        let res = match estimate_ancestry(&painel, &persons, &o) {
            Ok(x) => x,
            Err(e) => {
                t.check(&format!("ancestry, {} sites", 3 * por_crom), false, || e);
                continue;
            }
        };
        let segundos = tempo.elapsed().as_secs_f64();
        let j = 3 * por_crom;
        let tol = 4.0 / (j as f64).sqrt();
        let mut soma_erro = 0.0;
        for (p, a) in pessoas.iter().zip(&res) {
            let what = format!("ancestry {} with {j} sites", p.nome);
            let erro = p.q.iter().zip(&a.proportions).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max);
            soma_erro += erro;
            t.check(&format!("{what}: proportions within {tol:.3} of the truth"), erro <= tol, || format!("got {:?}, truth {:?}", a.proportions, p.q));
            t.near(&format!("{what}: proportions sum to 1"), a.proportions.iter().sum(), 1.0, 1e-9);
            let dados = dados_t(p, &s);
            t.eq(&format!("{what}: sites used"), a.sites_used, dados.len());
            t.eq(&format!("{what}: sites genotyped + assumed"), a.sites_genotyped + a.sites_assumed_ref, dados.len());
            t.eq(&format!("{what}: missing genotypes"), a.sites_missing_genotype, p.g.iter().filter(|g| g.is_none()).count());
            let ll_aqui = ll_t(&a.proportions, &dados);
            t.near(&format!("{what}: the reported log-likelihood is the one at the proportions"), a.log_likelihood, ll_aqui, 1e-6 * ll_aqui.abs() + 1e-6);
            for g in 0..3 {
                let mut e = vec![0.0; 3];
                e[g] = 1.0;
                let l = ll_t(&e, &dados);
                t.near(&format!("{what}: single-group log-likelihood of {}", grupos[g]), a.single_group_log_likelihoods[g], l, 1e-6 * l.abs() + 1e-6);
            }
            let melhor = (0..3).max_by(|&x, &y| ll_t(&{ let mut e = vec![0.0; 3]; e[x] = 1.0; e }, &dados).partial_cmp(&ll_t(&{ let mut e = vec![0.0; 3]; e[y] = 1.0; e }, &dados)).unwrap()).unwrap();
            t.eq(&format!("{what}: best single group"), a.best_group.clone(), grupos[melhor].to_string());
            t.check(&format!("{what}: converged"), a.converged, || format!("{} iterations", a.iterations));
            if ti < 2 {
                // The maximum: a grid search of the likelihood written here.
                // (The module keeps frequencies in f32: both likelihoods are
                // computed here, in f64, from the same frequencies.)
                let (qg, lg) = grade(&dados, 3);
                let d = qg.iter().zip(&a.proportions).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max);
                t.check(&format!("{what}: EM reaches the grid maximum"), ll_aqui >= lg - 1e-6 && d < 0.0015, || {
                    format!("EM {:?} ll {ll_aqui}, grid {qg:?} ll {lg}", a.proportions)
                });
            }
            if ti == 1 {
                t.eq(&format!("{what}: bootstrap replicates"), a.bootstrap_replicates, 60);
                let se_max = a.standard_errors.iter().cloned().fold(0.0, f64::max);
                t.check(&format!("{what}: standard errors present and small"), se_max > 0.0 && se_max < 0.05, || format!("{:?}", a.standard_errors));
            }
        }
        erros_por_tamanho.push((j, soma_erro / pessoas.len() as f64));
        println!("  ancestry: {} persons x {j} sites in {segundos:.2} s; mean max error {:.4}", pessoas.len(), soma_erro / pessoas.len() as f64);
    }
    if erros_por_tamanho.len() == 3 {
        t.check("ancestry: the error shrinks with more sites", erros_por_tamanho[2].1 < erros_por_tamanho[0].1, || format!("{erros_por_tamanho:?}"));
    }

    // ----- the bootstrap standard error against the Fisher information (K = 2)
    {
        let cromos: Vec<String> = (1..=10).map(|c| c.to_string()).collect();
        let cr: Vec<&str> = cromos.iter().map(String::as_str).collect();
        let s = sitios(r, &cr, 2500, 2, 2000, 0.2, 1_000_000);
        let painel = escreve_painel(t, "fisher_panel.tsv", &["A", "B"], &s);
        let p = simula_pessoa(r, "mix", &[0.7, 0.3], &s, 0.0);
        let vcf = escreve_pessoas(t, "fisher_person.vcf", &s, std::slice::from_ref(&p), "");
        let o = AncestryOptions { bootstrap: 200, block_bp: 500_000, seed: 7, ..Default::default() };
        match estimate_ancestry(&painel, &[Person { vcf: vcf.clone(), sample: None, callable: None }], &o) {
            Err(e) => t.check("bootstrap vs Fisher", false, || e),
            Ok(res) => {
                let a = &res[0];
                let q = a.proportions[0];
                let mut info = 0.0;
                for (g, x) in p.g.iter().zip(&s) {
                    if g.is_none() {
                        continue;
                    }
                    let (p1, p2) = (adj(x.p[0], x.an), adj(x.p[1], x.an));
                    let f = q * p1 + (1.0 - q) * p2;
                    info += 2.0 * (p1 - p2).powi(2) / (f * (1.0 - f));
                }
                let se_f = 1.0 / info.sqrt();
                let razao = a.standard_errors[0] / se_f;
                println!("  bootstrap SE {:.5} vs Fisher SE {se_f:.5} ({} blocks)", a.standard_errors[0], a.blocks);
                t.check("bootstrap standard error agrees with the Fisher information (K=2, 0.7/0.3)", (0.75..1.33).contains(&razao), || {
                    format!("bootstrap {} vs Fisher {se_f} (ratio {razao})", a.standard_errors[0])
                });
                t.near("bootstrap: the two standard errors of K=2 are equal", a.standard_errors[0], a.standard_errors[1], 1e-9);
                // 10 chromosomes of about 5 Mb, blocks of 500 kb.
                t.check("bootstrap: blocks counted", (90..=110).contains(&a.blocks), || format!("{} blocks", a.blocks));
                t.near("K=2 proportion near 0.7", q, 0.7, 4.0 * se_f);
                // Same seed, same standard errors; another seed, close.
                let res2 = estimate_ancestry(&painel, &[Person { vcf: vcf.clone(), sample: None, callable: None }], &o).unwrap();
                t.eq("bootstrap: reproducible with the same seed", res2[0].standard_errors.clone(), a.standard_errors.clone());
            }
        }
    }

    // ----- callable regions, absent sites, and the reasons a site is not used
    {
        let s = sitios(r, &["5", "6"], 3000, 3, 200, 0.18, 400);
        let painel = escreve_painel(t, "rules_panel.tsv.gz", &grupos, &s);
        let q = [0.2, 0.7, 0.1];
        let g: Vec<u8> = s.iter().map(|x| genotipo(r, &q, &x.p)).collect();
        // The callable BED covers chromosome 5 but a hole in the middle, and
        // nothing of chromosome 6 past its first half.
        let n5 = s.iter().filter(|x| x.crom == "5").count();
        let buraco = (s[n5 / 3].pos, s[n5 / 2].pos);
        let meio6 = s[n5 + (s.len() - n5) / 2].pos;
        let bed = t.write(
            "callable.bed",
            &[
                "track name=callable".to_string(),
                format!("chr5\t0\t{}", buraco.0 - 1),
                format!("chr5\t{}\t{}", buraco.1, 400_000_000u64),
                format!("chr6\t0\t{meio6}"),
                // Overlapping and unsorted intervals are merged.
                format!("chr6\t5\t{}", meio6 / 2),
            ],
        );
        let callable = |x: &SitioT| match x.crom.as_str() {
            "5" => x.pos < buraco.0 || x.pos > buraco.1,
            _ => x.pos <= meio6,
        };
        // Records: non-reference genotypes, plus planted cases.
        let mut recs = Vec::new();
        let mut esperado_uso = 0usize;
        let mut esp = HashMap::<&str, usize>::new();
        let mut dados = Vec::new();
        let mut j = 0;
        while j < s.len() {
            let x = &s[j];
            let mut base = Rec {
                crom: format!("chr{}", x.crom),
                pos: x.pos,
                ref_: x.ref_.to_string(),
                alt: x.alt.to_string(),
                filtro: "PASS".into(),
                formato: "GT".into(),
                gts: vec![],
                info: String::new(),
            };
            let caso = j % 50;
            let pf: Vec<f64> = x.p.iter().map(|&f| adj(f, x.an)).collect();
            // With a BED nothing outside it is used, records included: every
            // planted record outside the BED is "not callable".
            let cb = callable(x);
            let motivo = |k: &'static str, esp: &mut HashMap<&str, usize>| *esp.entry(if cb { k } else { "notcallable" }).or_default() += 1;
            match caso {
                // Another ALT at the site (the pair is never complementary,
                // so the complement of ALT is neither REF nor ALT): missing.
                7 => {
                    base.alt = comp(x.alt).to_string();
                    base.gts = vec!["0/1".into()];
                    motivo("other", &mut esp);
                    recs.push(base);
                }
                // The VCF's REF differs from the panel's.
                13 => {
                    base.ref_ = comp(x.ref_).to_string();
                    base.alt = x.ref_.to_string();
                    base.gts = vec!["0/1".into()];
                    motivo("refdiff", &mut esp);
                    recs.push(base);
                }
                // Explicit missing call.
                19 => {
                    base.gts = vec!["./.".into()];
                    motivo("missing", &mut esp);
                    recs.push(base);
                }
                // A filtered record.
                23 => {
                    base.filtro = "LowGQ".into();
                    base.gts = vec!["1/1".into()];
                    motivo("filtered", &mut esp);
                    recs.push(base);
                }
                // A deletion of the person, starting before the site, covers it
                // (REF of 3 bases from pos - 2: covers pos - 1 and pos only).
                // The site is absent: the deletion decides, inside the BED or not.
                29 if x.pos > s[j - 1].pos + 3 => {
                    recs.push(Rec { pos: x.pos - 2, ref_: "ACG".into(), alt: "A".into(), gts: vec!["0/1".into()], ..base.clone() });
                    *esp.entry("deletion").or_default() += 1;
                }
                // Haploid call (chrX-like): one allele.
                31 => {
                    let h = if r.chance(0.5) { "1" } else { "0" };
                    base.gts = vec![h.into()];
                    recs.push(base);
                    if cb {
                        esperado_uso += 1;
                        dados.push((if h == "1" { 3 } else { 4 }, pf)); // here: 3 haploid ALT, 4 haploid REF
                    } else {
                        motivo("", &mut esp);
                    }
                }
                _ => {
                    if g[j] > 0 {
                        base.gts = vec![if g[j] == 1 { "0/1" } else { "1/1" }.into()];
                        recs.push(base);
                    }
                    if cb {
                        esperado_uso += 1;
                        dados.push((g[j], pf));
                    } else {
                        motivo("", &mut esp);
                    }
                }
            }
            j += 1;
        }
        let vcf = escreve_vcf(t, "rules_person.vcf.gz", &["HG999".into()], &recs, &[]);
        let o = AncestryOptions { bootstrap: 0, ..Default::default() };
        let pessoa = Person { vcf: vcf.clone(), sample: None, callable: Some(bed.clone()) };
        match estimate_ancestry(&painel, &[pessoa], &o) {
            Err(e) => t.check("ancestry with a callable BED", false, || e),
            Ok(res) => {
                let a = &res[0];
                t.eq("callable BED: sites used", a.sites_used, esperado_uso);
                let got = (a.sites_other_allele, a.sites_ref_mismatch, a.sites_missing_genotype, a.sites_filtered, a.sites_in_deletion, a.sites_not_callable);
                let e = |k: &str| esp.get(k).copied().unwrap_or(0);
                t.eq("callable BED: reasons a site was not used", got, (e("other"), e("refdiff"), e("missing"), e("filtered"), e("deletion"), e("notcallable")));
                // The likelihood, with the haploid calls (codes 3, 4 here).
                let ll = |q: &[f64]| -> f64 {
                    dados
                        .iter()
                        .map(|(x, p)| {
                            let f: f64 = q.iter().zip(p).map(|(a, b)| a * b).sum();
                            match x {
                                3 => f.ln(),
                                4 => (1.0 - f).ln(),
                                0 => 2.0 * (1.0 - f).ln(),
                                1 => std::f64::consts::LN_2 + f.ln() + (1.0 - f).ln(),
                                _ => 2.0 * f.ln(),
                            }
                        })
                        .sum()
                };
                t.near("callable BED: log-likelihood (haploid calls included)", a.log_likelihood, ll(&a.proportions), 1e-6 * a.log_likelihood.abs());
                let mut ok = true;
                for i in 0..3 {
                    for dlt in [-1e-3, 1e-3] {
                        let mut qq = a.proportions.clone();
                        qq[i] += dlt;
                        qq[(i + 1) % 3] -= dlt;
                        if qq.iter().all(|&x| x >= 0.0) && ll(&qq) > a.log_likelihood + 1e-7 {
                            ok = false;
                        }
                    }
                }
                t.check("callable BED: no nearby point has a higher likelihood", ok, || format!("{:?}", a.proportions));
                t.check("callable BED: proportions near the truth", a.proportions.iter().zip(&q).all(|(x, y)| (x - y).abs() < 0.08), || format!("{:?}", a.proportions));
                t.check("callable BED: the summary says what was left out", a.resumo().contains("inside a deletion"), || a.resumo());
                t.check("callable BED: json", a.json().contains("\"sites_in_deletion\""), || a.json());
            }
        }
        // Without the BED: absent sites are 0/0 with absent_is_ref, missing without.
        let sem = Person { vcf: vcf.clone(), sample: Some("HG999".into()), callable: None };
        let a1 = estimate_ancestry(&painel, &[sem.clone()], &o).unwrap();
        let a2 = estimate_ancestry(&painel, &[sem.clone()], &AncestryOptions { absent_is_ref: false, ..o.clone() }).unwrap();
        let presentes = recs.iter().map(|x| (x.crom.clone(), x.pos)).collect::<HashSet<_>>();
        let nao_listados = s.iter().filter(|x| !presentes.contains(&(format!("chr{}", x.crom), x.pos))).count();
        t.eq("no BED, absent_is_ref: every absent site is 0/0", a1[0].sites_assumed_ref, nao_listados - e_count(&esp, "deletion"));
        t.eq("no BED, absent_is_ref off: absent sites are not callable", (a2[0].sites_assumed_ref, a2[0].sites_not_callable), (0, nao_listados - e_count(&esp, "deletion")));
        t.check("no BED, absent_is_ref: a note says so", a1[0].warnings.iter().any(|w| w.contains("absent_is_ref")), || format!("{:?}", a1[0].warnings));
        // pass_only off: the filtered records are used.
        let a3 = estimate_ancestry(&painel, &[Person { callable: Some(bed.clone()), ..sem.clone() }], &AncestryOptions { pass_only: false, ..o.clone() }).unwrap();
        t.eq("pass_only off: filtered records count", (a3[0].sites_filtered, a3[0].sites_used), (0, esperado_uso + e_count(&esp, "filtered")));
    }

    // ----- end to end: a panel built from simulated reference genotypes, and
    // held-out samples of the same VCF.
    {
        let k = 3;
        let por_grupo = 60;
        let s = sitios(r, &["1", "2"], 4000, k, 400, 0.2, 0);
        let mut amostras = Vec::new();
        let mut tabela = vec!["sample\tsuper".to_string()];
        let mut gts_por_amostra: Vec<Vec<u8>> = Vec::new();
        let mut verdade: Vec<Vec<f64>> = Vec::new();
        for (g, nome) in grupos.iter().enumerate() {
            for i in 0..por_grupo {
                let a = format!("{nome}_{i}");
                tabela.push(format!("{a}\t{nome}"));
                let mut q = vec![0.0; k];
                q[g] = 1.0;
                gts_por_amostra.push(s.iter().map(|x| genotipo(r, &q, &x.p)).collect());
                verdade.push(q);
                amostras.push(a);
            }
        }
        // Admixed held-out samples (in the table as AMR, and excluded).
        let mistos = [vec![0.5, 0.5, 0.0], vec![0.2, 0.3, 0.5], vec![0.0, 0.9, 0.1]];
        for (i, q) in mistos.iter().enumerate() {
            let a = format!("MIX_{i}");
            gts_por_amostra.push(s.iter().map(|x| genotipo(r, q, &x.p)).collect());
            verdade.push(q.clone());
            amostras.push(a);
        }
        let recs: Vec<Rec> = s
            .iter()
            .enumerate()
            .map(|(j, x)| Rec {
                crom: x.crom.clone(),
                pos: x.pos,
                ref_: x.ref_.to_string(),
                alt: x.alt.to_string(),
                filtro: "PASS".into(),
                formato: "GT".into(),
                gts: gts_por_amostra.iter().map(|g| ["0|0", "0|1", "1|1"][g[j] as usize].to_string()).collect(),
                info: String::new(),
            })
            .collect();
        let vcf = escreve_vcf(t, "e2e_ref.vcf.gz", &amostras, &recs, &[]);
        let tab = t.write("e2e_pop.tsv", &tabela);
        let fora: Vec<String> = ["AFR_0", "AFR_1", "EUR_0", "EUR_1", "EAS_0", "EAS_1"].iter().map(|x| x.to_string()).collect();
        let painel = t.p("e2e_panel.tsv.gz");
        let s1 = build_panel(&[vcf.clone()], &tab, "sample", "super", &fora, &painel, &PanelOptions { min_maf: 0.01, ..Default::default() });
        match s1 {
            Err(e) => t.check("end to end: panel", false, || e),
            Ok(sm) => {
                t.eq("end to end: samples per group (2 held out of each)", sm.samples_per_group.clone(), vec![58, 58, 58]);
                let mut testadas: Vec<usize> = (0..amostras.len()).filter(|&i| fora.contains(&amostras[i]) || amostras[i].starts_with("MIX")).collect();
                testadas.sort();
                let persons: Vec<Person> = testadas.iter().map(|&i| Person { vcf: vcf.clone(), sample: Some(amostras[i].clone()), callable: None }).collect();
                let tempo = Instant::now();
                match estimate_ancestry(&painel, &persons, &AncestryOptions { bootstrap: 30, block_bp: 300_000, ..Default::default() }) {
                    Err(e) => t.check("end to end: ancestry", false, || e),
                    Ok(res) => {
                        println!("  end to end: {} held-out samples in {:.2} s", res.len(), tempo.elapsed().as_secs_f64());
                        // Groups come out in alphabetical order: AFR, EAS, EUR.
                        let ordem = [0usize, 2, 1];
                        for (a, &i) in res.iter().zip(&testadas) {
                            let v = &verdade[i];
                            let v_alf: Vec<f64> = ordem.iter().map(|&g| v[g]).collect();
                            let erro = a.proportions.iter().zip(&v_alf).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max);
                            t.check(&format!("end to end: {} within 0.06 of the truth", amostras[i]), erro < 0.06, || format!("got {:?}, truth {v_alf:?}", a.proportions));
                            t.eq(&format!("end to end: {} sample name", amostras[i]), a.sample.clone(), amostras[i].clone());
                        }
                        // One person at a time gives the same as all together.
                        let um = estimate_ancestry(&painel, &persons[3..4], &AncestryOptions { bootstrap: 30, block_bp: 300_000, ..Default::default() }).unwrap();
                        t.eq("end to end: one person alone = the same person among others", um[0].proportions.clone(), res[3].proportions.clone());
                        t.eq("end to end: same bootstrap alone and among others", um[0].standard_errors.clone(), res[3].standard_errors.clone());
                    }
                }
            }
        }
    }

    // ----- errors
    {
        let s = sitios(r, &["1", "2"], 300, 2, 100, 0.2, 100);
        let painel = escreve_painel(t, "err_panel.tsv", &["A", "B"], &s);
        let p = simula_pessoa(r, "x", &[0.5, 0.5], &s, 0.0);
        let p2 = simula_pessoa(r, "y", &[0.5, 0.5], &s, 0.0);
        let bom = escreve_pessoas(t, "err_ok.vcf", &s, &[p, p2], "");
        let linhas = le_texto(&bom);
        let o = AncestryOptions { bootstrap: 0, ..Default::default() };
        let pr = |v: &str, s: Option<&str>| vec![Person { vcf: v.to_string(), sample: s.map(str::to_string), callable: None }];
        // Two positions swapped.
        let mut l = linhas.clone();
        let i = l.iter().position(|x| !x.starts_with('#')).unwrap() + 10;
        l.swap(i, i + 1);
        let v = t.write("err_unsorted.vcf", &l);
        t.fails("ancestry: unsorted person VCF", estimate_ancestry(&painel, &pr(&v, Some("x")), &o), &["sort_genomic", "chrom_order=\"natural\"", "comes after", &format!("line {}", i + 2)]);
        // Chromosome 1 again after chromosome 2.
        let mut l = linhas.clone();
        let ult = l.last().unwrap().clone();
        let mut campos: Vec<&str> = ult.split('\t').collect();
        campos[0] = "1";
        let nova = campos.join("\t");
        l.push(nova);
        let v = t.write("err_back.vcf", &l);
        t.fails("ancestry: chromosome in two blocks", estimate_ancestry(&painel, &pr(&v, Some("x")), &o), &["appears again", "sort_genomic"]);
        t.fails("ancestry: sample not in the VCF", estimate_ancestry(&painel, &pr(&bom, Some("HG002")), &o), &["HG002", "x, y"]);
        t.fails("ancestry: a multi-sample VCF without a sample name", estimate_ancestry(&painel, &pr(&bom, None), &o), &["2 samples", "name the sample"]);
        let ruim = t.write("err_panel_bad.tsv", &["chrom\tposition\tref\talt\tAF_A".to_string(), "1\t5\tA\tC\t0.1".to_string()]);
        t.fails("ancestry: a panel without the pos column", estimate_ancestry(&ruim, &pr(&bom, Some("x")), &o), &["chrom, pos, ref, alt", "build_panel"]);
        let ruim2 = t.write("err_panel_bad2.tsv", &["chrom\tpos\tref\talt\tAF_A\tAF_B".to_string(), "1\t5\tA\tC\t0.1\t1.7".to_string()]);
        t.fails("ancestry: a frequency above 1", estimate_ancestry(&ruim2, &pr(&bom, Some("x")), &o), &["line 2", "AF_B", "1.7"]);
        let bed_ruim = t.write("err.bed", &["chr1\tten\t20".to_string()]);
        t.fails(
            "ancestry: a malformed BED",
            estimate_ancestry(&painel, &[Person { vcf: bom.clone(), sample: Some("x".into()), callable: Some(bed_ruim) }], &o),
            &["line 1", "BED"],
        );
        // Another genome build: nothing in common.
        let outro = sitios(r, &["9"], 50, 2, 100, 0.2, 100);
        let painel9 = escreve_painel(t, "err_panel9.tsv", &["A", "B"], &outro);
        t.fails("ancestry: no site in common", estimate_ancestry(&painel9, &pr(&bom, Some("x")), &AncestryOptions { absent_is_ref: false, ..o.clone() }), &["none of the 50 panel SNVs", "genome build"]);
        t.fails("ancestry: no person", estimate_ancestry(&painel, &[], &o), &["No person"]);
        t.fails("ancestry: block_bp 0", estimate_ancestry(&painel, &pr(&bom, Some("x")), &AncestryOptions { block_bp: 0, ..o.clone() }), &["block_bp"]);
        t.fails("ancestry: missing VCF", estimate_ancestry(&painel, &pr(&t.p("nope.vcf"), Some("x")), &o), &["nope.vcf"]);
    }
}

fn e_count(m: &HashMap<&str, usize>, k: &str) -> usize {
    m.get(k).copied().unwrap_or(0)
}

/// The genotypes of the single sample of a person VCF written by
/// `escreve_pessoas`, per panel site (absent = 0/0).
fn dados_pessoa_de(s: &[SitioT], vcf: &str) -> Vec<Option<u8>> {
    let mut m: HashMap<(String, u64), Option<u8>> = HashMap::new();
    for l in le_texto(vcf) {
        if l.starts_with('#') {
            continue;
        }
        let c: Vec<&str> = l.split('\t').collect();
        let gt = c[9].split(':').next().unwrap();
        let g = match gt {
            "./." => None,
            x => Some(x.split(['/', '|']).filter(|a| *a == "1").count() as u8),
        };
        m.insert((c[0].to_string(), c[1].parse().unwrap()), g);
    }
    s.iter().map(|x| m.get(&(x.crom.clone(), x.pos)).copied().unwrap_or(Some(0))).collect()
}

// ---------------------------------------------------------------------
//  3. POLYGENIC SCORE
// ---------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Cat {
    GenoAlt,
    GenoRef,
    AbsentCallablePanel,
    AbsentNotCallablePanel,
    AbsentNotCallableNoPanel,
    AbsentCallableNoPanel,
    Ambiguous,
    Flip,
    Mismatch,
    NoPosition,
    Unsupported,
    OtherUnknown,
    MissingGt,
    Filtered,
    Indel,
    InDeletion,
}

const CATS: [Cat; 16] = [
    Cat::GenoAlt,
    Cat::GenoRef,
    Cat::AbsentCallablePanel,
    Cat::AbsentNotCallablePanel,
    Cat::AbsentNotCallableNoPanel,
    Cat::AbsentCallableNoPanel,
    Cat::Ambiguous,
    Cat::Flip,
    Cat::Mismatch,
    Cat::NoPosition,
    Cat::Unsupported,
    Cat::OtherUnknown,
    Cat::MissingGt,
    Cat::Filtered,
    Cat::Indel,
    Cat::InDeletion,
];

/// The expected fate of one variant for one sample.
#[derive(Clone, Debug, PartialEq)]
enum Fim {
    Usado { dosagem: f64, observado: bool, flip: bool, p: Option<Vec<f64>> },
    Fora(&'static str),
}

struct CasoPgs {
    linhas_escore: Vec<String>,
    linhas_painel: Vec<String>,
    recs: Vec<Rec>,
    bed: Vec<String>,
    fasta: String,
    /// Per variant row, per sample: the fate (None for rows that never reach
    /// the person: no position, unsupported, duplicate).
    fins: Vec<Vec<Option<Fim>>>,
    pesos: Vec<f64>,
    contas: HashMap<&'static str, usize>,
}

/// Builds the scoring file, panel, person VCF (two samples), BED and FASTA,
/// with each variant's fate known. `keep_amb`, `genoma`, `imputa` change the
/// expectation, not the files.
fn caso_pgs(r: &mut Rng, n: usize, keep_amb: bool, com_genoma: bool, imputa: &dyn Fn(&[f64], &[u64]) -> f64) -> CasoPgs {
    let mut c = CasoPgs {
        linhas_escore: vec![
            "###PGS CATALOG SCORING FILE - see https://www.pgscatalog.org/downloads/#dl_ftp_scoring for additional information".into(),
            "#format_version=2.0".into(),
            "##POLYGENIC SCORE (PGS) INFORMATION".into(),
            "#pgs_id=PGS999999".into(),
            "#pgs_name=TEST_SCORE".into(),
            "#trait_reported=Simulated trait".into(),
            "#genome_build=GRCh37".into(),
            "#weight_type=beta".into(),
            "##HARMONIZATION DETAILS".into(),
            "#HmPOS_build=GRCh38".into(),
            "rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\thm_source\thm_rsID\thm_chr\thm_pos\thm_inferOtherAllele".into(),
        ],
        linhas_painel: vec!["chrom\tpos\tref\talt\tAF_G1\tAF_G2\tAN_G1\tAN_G2".into()],
        recs: Vec::new(),
        bed: Vec::new(),
        fasta: String::new(),
        fins: Vec::new(),
        pesos: Vec::new(),
        contas: HashMap::new(),
    };
    // The reference sequence: random, with the REF of every variant in place.
    let tamanho = 200 * (n + 10);
    let mut genoma: Vec<char> = (0..tamanho).map(|_| r.base()).collect();
    let mut pos = 0u64;
    for i in 0..n {
        pos += 200;
        let cat = CATS[r.below(CATS.len())];
        let w = ((r.normal() * 0.2) * 1e4).round() / 1e4;
        let (mut rb, mut ab) = r.par_seguro();
        if cat == Cat::Ambiguous {
            rb = ['A', 'C', 'G', 'T'][r.below(4)];
            ab = comp(rb);
        }
        genoma[pos as usize - 1] = rb;
        let (ref_s, alt_s) = if cat == Cat::Indel { (rb.to_string(), format!("{rb}TT")) } else { (rb.to_string(), ab.to_string()) };
        let af: Vec<f64> = (0..2).map(|_| ((0.02 + 0.96 * r.f()) * 1e6).round() / 1e6).collect();
        let an: Vec<u64> = vec![200 + r.below(800) as u64, 100 + r.below(400) as u64];
        let no_painel = match cat {
            Cat::AbsentCallablePanel | Cat::AbsentNotCallablePanel => true,
            Cat::AbsentNotCallableNoPanel | Cat::AbsentCallableNoPanel | Cat::NoPosition | Cat::Unsupported => false,
            _ => r.chance(0.6),
        };
        if no_painel {
            c.linhas_painel.push(format!("chr1\t{pos}\t{ref_s}\t{alt_s}\t{}\t{}\t{}\t{}", af[0], af[1], an[0], an[1]));
        }
        // Effect allele: ALT or REF.
        let efeito_alt = r.chance(0.5);
        let (mut e, mut o) = if efeito_alt { (alt_s.clone(), ref_s.clone()) } else { (ref_s.clone(), alt_s.clone()) };
        let mut infer = String::new();
        match cat {
            Cat::Flip => {
                e = e.chars().map(comp).collect();
                o = o.chars().map(comp).collect();
            }
            Cat::Mismatch => {
                e = format!("{rb}A");
                o = rb.to_string();
            }
            Cat::OtherUnknown => {
                infer = format!("{o}/{}", if o == "A" { "C" } else { "A" });
                o = String::new();
            }
            Cat::Unsupported => e = "HLA-A*01:01".into(),
            _ => {}
        }
        let p_efeito: Vec<f64> = af.iter().map(|&p| if efeito_alt { p } else { 1.0 - p }).collect();
        // The person's records: two samples with their own genotypes.
        let gts: Vec<u8> = (0..2).map(|_| r.below(3) as u8).collect();
        let txt_gt = |g: u8| ["0/0", "0|1", "1/1"][g as usize].to_string();
        let callable = !matches!(cat, Cat::AbsentNotCallablePanel | Cat::AbsentNotCallableNoPanel);
        if callable {
            c.bed.push(format!("chr1\t{}\t{}", pos - 1, pos));
        }
        let rec = |gts: Vec<String>, filtro: &str| Rec {
            crom: "chr1".into(),
            pos,
            ref_: ref_s.clone(),
            alt: alt_s.clone(),
            filtro: filtro.into(),
            formato: "GT".into(),
            gts,
            info: String::new(),
        };
        let mut del_carrega = [false, false];
        match cat {
            Cat::GenoAlt | Cat::GenoRef | Cat::Ambiguous | Cat::Flip | Cat::Mismatch | Cat::OtherUnknown | Cat::Indel => {
                c.recs.push(rec(gts.iter().map(|&g| txt_gt(g)).collect(), "PASS"));
            }
            Cat::MissingGt => c.recs.push(rec(vec!["./.".into(), txt_gt(gts[1])], "PASS")),
            Cat::Filtered => c.recs.push(rec(gts.iter().map(|&g| txt_gt(g)).collect(), "LowQual")),
            Cat::InDeletion => {
                // Sample 0 carries a deletion that covers the site; sample 1 does not.
                del_carrega = [true, false];
                c.recs.push(Rec {
                    pos: pos - 2,
                    ref_: format!("{}{}{}{rb}", genoma[pos as usize - 3], genoma[pos as usize - 2], rb),
                    alt: genoma[pos as usize - 3].to_string(),
                    gts: vec!["0/1".into(), "0/0".into()],
                    ..rec(vec![], "PASS")
                });
            }
            _ => {}
        }
        if cat == Cat::InDeletion {
            // The deletion REF must match the genome: fix the two bases before.
            let d = &c.recs.last().unwrap().ref_;
            let b: Vec<char> = d.chars().collect();
            genoma[pos as usize - 3] = b[0];
            genoma[pos as usize - 2] = b[1];
        }
        let hm_pos = if cat == Cat::NoPosition { String::new() } else { pos.to_string() };
        let linha = format!("rs{i}\t1\t{}\t{e}\t{o}\t{w}\tENSEMBL\trs{i}\t1\t{hm_pos}\t{infer}", pos + 777);
        let duplicada = r.chance(0.02) && !matches!(cat, Cat::NoPosition | Cat::Unsupported);
        c.linhas_escore.push(linha.clone());
        c.pesos.push(w);
        // Expected fate per sample.
        let fim = |s: usize| -> Option<Fim> {
            let x_alt = gts[s] as f64;
            let dos_g = if efeito_alt { x_alt } else { 2.0 - x_alt };
            let p = if no_painel { Some(p_efeito.clone()) } else { None };
            let imput = || -> Fim {
                if no_painel {
                    Fim::Usado { dosagem: 2.0 * imputa(&p_efeito, &an), observado: false, flip: false, p: None }
                } else {
                    Fim::Fora("not_callable")
                }
            };
            Some(match cat {
                Cat::NoPosition | Cat::Unsupported => return None,
                Cat::GenoAlt | Cat::GenoRef | Cat::Indel | Cat::OtherUnknown => Fim::Usado { dosagem: dos_g, observado: true, flip: false, p },
                Cat::Flip => Fim::Usado { dosagem: dos_g, observado: true, flip: true, p },
                Cat::Ambiguous => {
                    if keep_amb {
                        Fim::Usado { dosagem: dos_g, observado: true, flip: false, p }
                    } else {
                        Fim::Fora("ambiguous_skipped")
                    }
                }
                Cat::Mismatch => Fim::Fora("mismatched"),
                Cat::AbsentCallablePanel => Fim::Usado { dosagem: if efeito_alt { 0.0 } else { 2.0 }, observado: true, flip: false, p },
                Cat::AbsentCallableNoPanel => {
                    if com_genoma {
                        Fim::Usado { dosagem: if efeito_alt { 0.0 } else { 2.0 }, observado: true, flip: false, p: None }
                    } else {
                        Fim::Fora("ref_unknown")
                    }
                }
                Cat::AbsentNotCallablePanel | Cat::AbsentNotCallableNoPanel => imput(),
                Cat::MissingGt => {
                    if s == 0 {
                        imput()
                    } else {
                        Fim::Usado { dosagem: dos_g, observado: true, flip: false, p }
                    }
                }
                Cat::Filtered => imput(),
                Cat::InDeletion => {
                    if del_carrega[s] {
                        imput()
                    } else if no_painel || com_genoma {
                        // Absent and callable: homozygous reference, the REF
                        // told by the panel or by the genome.
                        Fim::Usado { dosagem: if efeito_alt { 0.0 } else { 2.0 }, observado: true, flip: false, p }
                    } else {
                        Fim::Fora("ref_unknown")
                    }
                }
            })
        };
        c.fins.push(vec![fim(0), fim(1)]);
        *c.contas.entry(match cat {
            Cat::NoPosition => "no_position",
            Cat::Unsupported => "unsupported",
            _ => "ok",
        })
        .or_default() += 1;
        if duplicada {
            c.linhas_escore.push(linha);
            c.pesos.push(w);
            c.fins.push(vec![None, None]);
            *c.contas.entry("duplicates").or_default() += 1;
        }
    }
    let mut fa = String::from(">chr1\n");
    for (i, b) in genoma.iter().enumerate() {
        fa.push(*b);
        if i % 60 == 59 {
            fa.push('\n');
        }
    }
    fa.push('\n');
    c.fasta = fa;
    c
}

fn normal_t(z: f64) -> f64 {
    // Simpson's rule on the density from -12 to z.
    let a = -12.0;
    if z <= a {
        return 0.0;
    }
    let n = 20000;
    let h = (z - a) / n as f64;
    let f = |x: f64| (-0.5 * x * x).exp() / (2.0 * std::f64::consts::PI).sqrt();
    let mut s = f(a) + f(z);
    for i in 1..n {
        s += f(a + i as f64 * h) * if i % 2 == 1 { 4.0 } else { 2.0 };
    }
    s * h / 3.0
}

fn confere_pgs(t: &mut T, what: &str, c: &CasoPgs, res: &[ScoreResult]) {
    t.eq(&format!("{what}: one result per sample"), res.len(), 2);
    for (s, a) in res.iter().enumerate() {
        let w = format!("{what}, sample {s}");
        let mut escore = 0.0;
        let mut observado = 0.0;
        let mut n = HashMap::<&str, usize>::new();
        let mut grupos = [(0usize, 0.0f64, 0.0f64, 0.0f64); 2];
        let mut flips = 0;
        for (i, f) in c.fins.iter().enumerate() {
            let Some(f) = &f[s] else { continue };
            let peso = c.pesos[i];
            match f {
                Fim::Usado { dosagem, observado: obs, flip, p } => {
                    escore += dosagem * peso;
                    if *flip {
                        flips += 1;
                    }
                    if *obs {
                        observado += dosagem * peso;
                        *n.entry("observed").or_default() += 1;
                        if let Some(p) = p {
                            for g in 0..2 {
                                grupos[g].0 += 1;
                                grupos[g].1 += dosagem * peso;
                                grupos[g].2 += 2.0 * p[g] * peso;
                                grupos[g].3 += 2.0 * p[g] * (1.0 - p[g]) * peso * peso;
                            }
                        }
                    } else {
                        *n.entry("imputed").or_default() += 1;
                    }
                }
                Fim::Fora(k) => *n.entry(k).or_default() += 1,
            }
        }
        let e = |k: &str| n.get(k).copied().unwrap_or(0);
        t.near(&format!("{w}: score = direct sum"), a.score, escore, 1e-9);
        t.near(&format!("{w}: observed score"), a.score_observed, observado, 1e-9);
        t.eq(&format!("{w}: genotyped + assumed_ref"), a.genotyped + a.assumed_ref, e("observed"));
        t.eq(
            &format!("{w}: imputed, ambiguous skipped, mismatched, REF unknown, not callable"),
            (a.imputed, a.ambiguous_skipped, a.mismatched, a.ref_unknown, a.not_callable),
            (e("imputed"), e("ambiguous_skipped"), e("mismatched"), e("ref_unknown"), e("not_callable")),
        );
        t.eq(&format!("{w}: flipped"), a.flipped, flips);
        let ct = |k: &str| c.contas.get(k).copied().unwrap_or(0);
        t.eq(&format!("{w}: rows, no position, unsupported, duplicates"), (a.variants, a.no_position, a.unsupported_allele, a.duplicates), (c.pesos.len(), ct("no_position"), ct("unsupported"), ct("duplicates")));
        t.eq(
            &format!("{w}: every usable variant lands in one bucket"),
            a.usable,
            a.genotyped + a.assumed_ref + a.imputed + a.ambiguous_skipped + a.mismatched + a.ref_unknown + a.not_callable,
        );
        t.near(&format!("{w}: coverage"), a.coverage, e("observed") as f64 / c.pesos.len() as f64, 1e-12);
        if !a.groups.is_empty() {
            for g in 0..2 {
                let (nv, ob, m, v) = grupos[g];
                let gr = &a.groups[g];
                t.eq(&format!("{w}: group {g} variants"), gr.variants, nv);
                t.near(&format!("{w}: group {g} observed"), gr.observed, ob, 1e-9);
                t.near(&format!("{w}: group {g} mean = sum 2pw"), gr.mean, m, 1e-6);
                t.near(&format!("{w}: group {g} sd = sqrt(sum 2p(1-p)w^2)"), gr.sd, v.sqrt(), 1e-6);
                let z = (ob - m) / v.sqrt();
                t.near(&format!("{w}: group {g} z"), gr.z, z, 1e-6);
                t.near(&format!("{w}: group {g} percentile = normal integral"), gr.percentile, 100.0 * normal_t(z), 2e-5);
            }
        }
    }
}

fn testa_pgs(t: &mut T, r: &mut Rng) {
    let pooled = |p: &[f64], an: &[u64]| (p[0] * an[0] as f64 + p[1] * an[1] as f64) / (an[0] + an[1]) as f64;
    for (ci, (keep, genoma)) in [(false, false), (true, false), (false, true)].iter().enumerate() {
        let c = caso_pgs(r, 2500, *keep, *genoma, &pooled);
        let escore = t.write(&format!("pgs{ci}.txt"), &c.linhas_escore);
        let painel = t.write(&format!("pgs_panel{ci}.tsv"), &c.linhas_painel);
        let vcf = escreve_vcf(t, &format!("pgs_person{ci}.vcf.gz"), &["P1".into(), "P2".into()], &c.recs, &[]);
        let bed = t.write(&format!("pgs_callable{ci}.bed"), &c.bed);
        let fasta = t.p(&format!("pgs_ref{ci}.fa"));
        std::fs::write(&fasta, &c.fasta).unwrap();
        let persons: Vec<Person> = ["P1", "P2"].iter().map(|s| Person { vcf: vcf.clone(), sample: Some(s.to_string()), callable: Some(bed.clone()) }).collect();
        let det = t.p(&format!("pgs_details{ci}.tsv"));
        let o = ScoreOptions {
            panel: Some(painel.clone()),
            reference: if *genoma { Some(fasta.clone()) } else { None },
            keep_ambiguous: *keep,
            details: Some(det.clone()),
            ..Default::default()
        };
        let what = format!("polygenic score case {ci} (keep_ambiguous {keep}, reference {genoma})");
        let tempo = Instant::now();
        match polygenic_score(&escore, &persons, &o) {
            Err(e) => t.check(&what, false, || e),
            Ok(res) => {
                println!("  {what}: {:.2} s", tempo.elapsed().as_secs_f64());
                confere_pgs(t, &what, &c, &res);
                t.eq(&format!("{what}: metadata"), (res[0].pgs_id.clone(), res[0].genome_build.clone(), res[0].positions_from.clone()), ("PGS999999".to_string(), "GRCh38".to_string(), "hm_pos".to_string()));
                let d = le_texto(&det);
                let usados = res.iter().map(|a| a.genotyped + a.assumed_ref + a.imputed + a.ambiguous_skipped + a.mismatched + a.ref_unknown + a.not_callable).sum::<usize>();
                t.eq(&format!("{what}: details table, one row per variant and sample"), d.len() - 1, usados);
                let soma: f64 = d[1..].iter().filter_map(|l| l.split('\t').nth(9).and_then(|x| x.parse::<f64>().ok())).sum();
                t.near(&format!("{what}: details contributions add up to the scores"), soma, res[0].score + res[1].score, 1e-6);
                t.check(&format!("{what}: summary text"), res[0].resumo().contains("PGS999999") && res[0].resumo().contains("percentile"), || res[0].resumo());
                t.check(&format!("{what}: json"), res[0].json().contains("\"coverage\"") && res[0].json().contains("\"groups\":[{"), || res[0].json());
            }
        }
        if ci == 0 {
            // Imputation from one group, and from proportions.
            let g2 = |p: &[f64], _: &[u64]| p[1];
            let mut r2 = Rng(4242);
            let c2 = caso_pgs(&mut r2, 800, false, false, &g2);
            let e2 = t.write("pgs_g2.txt", &c2.linhas_escore);
            let p2 = t.write("pgs_g2_panel.tsv", &c2.linhas_painel);
            let v2 = escreve_vcf(t, "pgs_g2.vcf", &["P1".into(), "P2".into()], &c2.recs, &[]);
            let b2 = t.write("pgs_g2.bed", &c2.bed);
            let ps: Vec<Person> = ["P1", "P2"].iter().map(|s| Person { vcf: v2.clone(), sample: Some(s.to_string()), callable: Some(b2.clone()) }).collect();
            match polygenic_score(&e2, &ps, &ScoreOptions { panel: Some(p2.clone()), impute: Impute::Group("G2".into()), ..Default::default() }) {
                Ok(res) => confere_pgs(t, "polygenic score imputed from group G2", &c2, &res),
                Err(e) => t.check("polygenic score imputed from group G2", false, || e),
            }
            let mut r3 = Rng(4242);
            let mix = |p: &[f64], _: &[u64]| 0.25 * p[0] + 0.75 * p[1];
            let c3 = caso_pgs(&mut r3, 800, false, false, &mix);
            match polygenic_score(&e2, &ps, &ScoreOptions { panel: Some(p2.clone()), impute: Impute::Proportions(vec![1.0, 3.0]), ..Default::default() }) {
                Ok(res) => confere_pgs(t, "polygenic score imputed from proportions 1:3", &c3, &res),
                Err(e) => t.check("polygenic score imputed from proportions", false, || e),
            }
            t.fails("polygenic score: unknown imputation group", polygenic_score(&e2, &ps, &ScoreOptions { panel: Some(p2.clone()), impute: Impute::Group("EUR".into()), ..Default::default() }), &["EUR", "G1, G2"]);
            t.fails(
                "polygenic score: proportions of the wrong length",
                polygenic_score(&e2, &ps, &ScoreOptions { panel: Some(p2.clone()), impute: Impute::Proportions(vec![1.0]), ..Default::default() }),
                &["one per panel group"],
            );
            // Without a panel: nothing imputed, no group, a note.
            match polygenic_score(&escore, &persons[..1], &ScoreOptions::default()) {
                Err(e) => t.check("polygenic score without panel", false, || e),
                Ok(res) => {
                    let a = &res[0];
                    t.check("polygenic score without panel: nothing imputed, no group", a.imputed == 0 && a.groups.is_empty(), || format!("{} imputed", a.imputed));
                    t.check("polygenic score without panel: a note on the variants left out", a.warnings.iter().any(|w| w.contains("give a panel")), || format!("{:?}", a.warnings));
                }
            }
        }
    }

    // ----- unsupported layouts and malformed files
    let base: Vec<String> = vec!["#pgs_id=PGS1".into(), "rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight".into(), "rs1\t1\t100\tA\tG\t0.1".into()];
    let p = Person { vcf: t.p("pgs_person0.vcf.gz"), sample: Some("P1".into()), callable: None };
    let ps = vec![p];
    let o = ScoreOptions::default();
    let f = t.write("bad_dosage.txt", &["#pgs_id=PGS2".into(), "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\tdosage_0_weight\tdosage_1_weight\tdosage_2_weight".into(), "1\t5\tA\tG\t0\t0\t0.2\t0.5".into()]);
    t.fails("polygenic score: dosage-specific weights refused", polygenic_score(&f, &ps, &o), &["dosage-specific", "not supported"]);
    let f = t.write("bad_haplo.txt", &["chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\tis_haplotype".into(), "1\t5\tA\tG\t0.1\tFalse".into(), "1\t9\tA\tG\t0.1\tTrue".into()]);
    t.fails("polygenic score: haplotype weights refused", polygenic_score(&f, &ps, &o), &["line 3", "is_haplotype", "not supported"]);
    let f = t.write("bad_rsid.txt", &["rsID\teffect_allele\tother_allele\teffect_weight".into(), "rs1\tA\tG\t0.1".into()]);
    t.fails("polygenic score: rsIDs only", polygenic_score(&f, &ps, &o), &["no variant positions", "_hmPOS_GRCh38"]);
    let f = t.write("bad_weight.txt", &[base[0].clone(), base[1].clone(), "rs1\t1\t100\tA\tG\tbig".into()]);
    t.fails("polygenic score: a weight that is not a number", polygenic_score(&f, &ps, &o), &["line 3", "big"]);
    let f = t.write("bad_cols.txt", &["chr_name\tchr_position\teffect_allele\tweight".into()]);
    t.fails("polygenic score: no effect_weight", polygenic_score(&f, &ps, &o), &["effect_weight", "weight"]);
    let f = t.write("ok_unharmonized.txt", &base);
    match polygenic_score(&f, &ps, &o) {
        Ok(res) => t.check("polygenic score: author positions give a note on the build", res[0].warnings.iter().any(|w| w.contains("genome_build")), || format!("{:?}", res[0].warnings)),
        Err(e) => t.check("polygenic score: author positions", false, || e),
    }
    let f2 = t.write("ok_or.txt", &["#weight_type=OR".into(), base[1].clone(), base[2].clone()]);
    match polygenic_score(&f2, &ps, &o) {
        Ok(res) => t.check("polygenic score: weight_type OR gives a note", res[0].warnings.iter().any(|w| w.contains("log(OR)")), || format!("{:?}", res[0].warnings)),
        Err(e) => t.check("polygenic score: weight_type OR", false, || e),
    }
    t.fails("polygenic score: details over the input", polygenic_score(&f, &ps, &ScoreOptions { details: Some(f.clone()), ..Default::default() }), &["also one of the input files"]);
}

// ---------------------------------------------------------------------
//  REAL DATA
// ---------------------------------------------------------------------

/// Arguments: real [chr22 VCF]. The data folder (1000genomes/, giab/,
/// pgs_catalog/, as MANIFEST.md of the downloads lays them out) is
/// MMG_HUMAN_GENOMES (default: ./human_genomes); the run is skipped for files
/// that are missing or still growing.
fn real(args: &[String]) {
    let base = std::env::var("MMG_HUMAN_GENOMES").unwrap_or_else(|_| "human_genomes".to_string());
    let kg = format!("{base}/1000genomes");
    let vcf22 = args.first().cloned().unwrap_or_else(|| format!("{kg}/nygc_30x_phased_panel/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"));
    let pop = format!("{kg}/nygc_30x_metadata/20130606_g1k_3202_samples_ped_population.txt");
    let ped = format!("{kg}/nygc_30x_metadata/1kGP.3202_samples.pedigree_info.txt");
    let out = std::env::var("MMG_OUT").unwrap_or_else(|_| std::env::temp_dir().join("mmg_ancestry_real").to_string_lossy().to_string());
    std::fs::create_dir_all(&out).unwrap();
    let estavel = |p: &str| -> bool {
        let a = std::fs::metadata(p).map(|m| m.len()).unwrap_or(0);
        std::thread::sleep(std::time::Duration::from_secs(3));
        let b = std::fs::metadata(p).map(|m| m.len()).unwrap_or(0);
        a > 0 && a == b
    };
    let painel = format!("{out}/panel_chr22_5sp.tsv.gz");
    let pessoas_1kg: Vec<(String, String, String)>;
    {
        // Children of trios out (they are relatives); 5 samples of each
        // population held out.
        let mut filhos: HashSet<String> = HashSet::new();
        for l in std::fs::read_to_string(&ped).unwrap().lines().skip(1) {
            let c: Vec<&str> = l.split_whitespace().collect();
            if c.len() >= 3 && (c[1] != "0" || c[2] != "0") {
                filhos.insert(c[0].to_string());
            }
        }
        let mut por_pop: HashMap<String, Vec<(String, String)>> = HashMap::new();
        for l in std::fs::read_to_string(&pop).unwrap().lines().skip(1) {
            let c: Vec<&str> = l.split_whitespace().collect();
            if !filhos.contains(c[1]) {
                por_pop.entry(c[5].to_string()).or_default().push((c[1].to_string(), c[6].to_string()));
            }
        }
        let mut pops: Vec<&String> = por_pop.keys().collect();
        pops.sort();
        let mut fora = Vec::new();
        for p in pops {
            for (s, sp) in por_pop[p].iter().take(5) {
                fora.push((s.clone(), p.clone(), sp.clone()));
            }
        }
        pessoas_1kg = fora;
        if Path::new(&painel).exists() && std::env::var("MMG_REBUILD").is_err() {
            println!("panel already built: {painel}");
        } else {
            if !estavel(&vcf22) {
                println!("{vcf22} is missing or still growing: real run skipped");
                return;
            }
            let mut excluir: Vec<String> = filhos.iter().cloned().collect();
            excluir.extend(pessoas_1kg.iter().map(|x| x.0.clone()));
            let t0 = Instant::now();
            match build_panel(&[vcf22.clone()], &pop, "SampleID", "Superpopulation", &excluir, &painel, &PanelOptions::default()) {
                Ok(s) => println!("panel in {:.1} s\n{}", t0.elapsed().as_secs_f64(), s.resumo()),
                Err(e) => {
                    println!("panel failed: {e}");
                    return;
                }
            }
        }
    }
    // Held-out 1000 Genomes samples.
    let persons: Vec<Person> = pessoas_1kg.iter().map(|x| Person { vcf: vcf22.clone(), sample: Some(x.0.clone()), callable: None }).collect();
    let o = AncestryOptions { block_bp: 2_000_000, absent_is_ref: false, ..Default::default() };
    // MMG_ONLY_PGS: skip the chr22 runs (the scores need the genome-wide
    // ancestry only).
    let so_pgs = std::env::var("MMG_ONLY_PGS").is_ok();
    let t0 = Instant::now();
    match if so_pgs { Err("skipped (MMG_ONLY_PGS)".to_string()) } else { estimate_ancestry(&painel, &persons, &o) } {
        Ok(res) => {
            println!("ancestry of {} held-out samples in {:.1} s (one pass over the chr22 VCF)", res.len(), t0.elapsed().as_secs_f64());
            println!("sample\tpop\tsuperpop\t{}\tbest\tsites", res[0].groups.join("\t"));
            for (a, x) in res.iter().zip(&pessoas_1kg) {
                let q: Vec<String> = a.proportions.iter().zip(&a.standard_errors).map(|(p, e)| format!("{p:.3}+/-{e:.3}")).collect();
                println!("{}\t{}\t{}\t{}\t{}\t{}", a.sample, x.1, x.2, q.join("\t"), a.best_group, a.sites_used);
            }
        }
        Err(e) => println!("ancestry failed: {e}"),
    }
    // GIAB genomes, restricted to chr22 by the panel.
    let giab = format!("{base}/giab");
    let mut gs = Vec::new();
    for (d, h) in [
        ("HG001_NA12878", "HG001"),
        ("HG002_NA24385_son", "HG002"),
        ("HG003_NA24149_father", "HG003"),
        ("HG004_NA24143_mother", "HG004"),
        ("HG005_NA24631_son", "HG005"),
        ("HG006_NA24694_father", "HG006"),
        ("HG007_NA24695_mother", "HG007"),
    ] {
        let v = format!("{giab}/{d}/v4.2.1/{h}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz");
        let b1 = format!("{giab}/{d}/v4.2.1/{h}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed");
        let b2 = format!("{giab}/{d}/v4.2.1/{h}_GRCh38_1_22_v4.2.1_benchmark.bed");
        let bed = if Path::new(&b1).exists() { b1 } else { b2 };
        if Path::new(&v).exists() && Path::new(&bed).exists() {
            gs.push(Person { vcf: v, sample: None, callable: Some(bed) });
        }
    }
    let t0 = Instant::now();
    match if so_pgs { Err("skipped (MMG_ONLY_PGS)".to_string()) } else { estimate_ancestry(&painel, &gs, &o) } {
        Ok(res) => {
            println!("GIAB ancestry on chr22 in {:.1} s", t0.elapsed().as_secs_f64());
            for a in &res {
                print!("{}", a.resumo());
            }
        }
        Err(e) => println!("GIAB ancestry failed: {e}"),
    }
    // Genome-wide ancestry of the GIAB genomes from the phase 3 frequencies
    // (2 decimals), one site per 20 kb with MAF >= 0.05 in some group.
    let sites = format!("{kg}/phase3_grch38_sites_20190312/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz");
    let g: Vec<String> = ["AFR", "AMR", "EAS", "EUR", "SAS"].iter().map(|x| x.to_string()).collect();
    let f: Vec<String> = g.iter().map(|x| format!("{x}_AF")).collect();
    let painel_gw_anc = format!("{out}/panel_gw_anc_5sp.tsv.gz");
    if !Path::new(&painel_gw_anc).exists() && estavel(&sites) {
        let t0 = Instant::now();
        match panel_from_frequencies(&[sites.clone()], &g, &f, &[], &painel_gw_anc, &PanelOptions { thin_bp: 20_000, ..Default::default() }) {
            Ok(s) => println!("genome-wide ancestry panel in {:.1} s\n{}", t0.elapsed().as_secs_f64(), s.resumo()),
            Err(e) => println!("genome-wide ancestry panel failed: {e}"),
        }
    }
    let mut q_hg002: Option<Vec<f64>> = None;
    if Path::new(&painel_gw_anc).exists() {
        let t0 = Instant::now();
        match estimate_ancestry(&painel_gw_anc, &gs, &AncestryOptions { block_bp: 5_000_000, ..Default::default() }) {
            Ok(res) => {
                println!("GIAB ancestry genome-wide (phase 3 frequencies) in {:.1} s", t0.elapsed().as_secs_f64());
                for a in &res {
                    print!("{}", a.resumo());
                }
                q_hg002 = res.iter().find(|a| a.sample == "HG002").map(|a| a.proportions.clone());
            }
            Err(e) => println!("GIAB genome-wide ancestry failed: {e}"),
        }
    }
    // Polygenic score on HG002, with a genome-wide panel if one is there.
    let painel_gw = format!("{out}/panel_gw_5sp.tsv.gz");
    if !Path::new(&painel_gw).exists() && estavel(&sites) {
        let t0 = Instant::now();
        match panel_from_frequencies(&[sites.clone()], &g, &f, &[], &painel_gw, &PanelOptions { min_maf: 0.0, keep_indels: true, ..Default::default() }) {
            Ok(s) => println!("genome-wide panel in {:.1} s\n{}", t0.elapsed().as_secs_f64(), s.resumo()),
            Err(e) => println!("genome-wide panel failed: {e}"),
        }
    }
    let hg002 = gs.iter().find(|p| p.vcf.contains("HG002")).cloned();
    if let Some(p) = hg002 {
        for id in ["PGS000004", "PGS000662", "PGS000832", "PGS003438"] {
            let f = format!("{base}/pgs_catalog/{id}/{id}_hmPOS_GRCh38.txt.gz");
            if !Path::new(&f).exists() {
                continue;
            }
            let o = ScoreOptions {
                panel: if Path::new(&painel_gw).exists() { Some(painel_gw.clone()) } else { None },
                details: Some(format!("{out}/{id}_HG002_details.tsv")),
                ..Default::default()
            };
            let t0 = Instant::now();
            match polygenic_score(&f, &[p.clone()], &o) {
                Ok(r) => print!("{:.1} s\n{}", t0.elapsed().as_secs_f64(), r[0].resumo()),
                Err(e) => println!("{id} failed: {e}"),
            }
            // Keeping the strand-ambiguous SNVs (harmonized files are on the
            // forward strand of GRCh38).
            let o3 = ScoreOptions { keep_ambiguous: true, details: None, ..o.clone() };
            if let Ok(r) = polygenic_score(&f, &[p.clone()], &o3) {
                let a = &r[0];
                let z: Vec<String> = a.groups.iter().map(|g| format!("{} z {:.2}", g.group, g.z)).collect();
                println!("  {id} keeping A/T and C/G SNVs: score {:.6}, coverage {:.1}%; {}", a.score, 100.0 * a.coverage, z.join(", "));
            }
            // The same, imputing with HG002's own genome-wide proportions.
            if let (Some(q), Some(_)) = (&q_hg002, &o.panel) {
                let o2 = ScoreOptions { impute: Impute::Proportions(q.clone()), details: None, ..o.clone() };
                if let Ok(r) = polygenic_score(&f, &[p.clone()], &o2) {
                    let a = &r[0];
                    println!("  {id} imputed with HG002's proportions: score {:.6} ({} imputed; observed part {:.6})", a.score, a.imputed, a.score_observed);
                }
            }
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(String::as_str) == Some("real") {
        real(&args[2..]);
        return;
    }
    let dir = std::env::temp_dir().join(format!("trial_ancestry_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut t = T { dir: dir.clone(), total: 0, ok: 0, falhas: Vec::new() };
    let mut r = Rng(20261114);
    let t0 = Instant::now();
    testa_painel(&mut t, &mut r);
    println!("  panels: {} checks so far, {:.1} s", t.total, t0.elapsed().as_secs_f64());
    testa_ancestralidade(&mut t, &mut r);
    println!("  ancestry: {} checks so far, {:.1} s", t.total, t0.elapsed().as_secs_f64());
    testa_pgs(&mut t, &mut r);
    let _ = std::fs::remove_dir_all(&dir);
    println!("{} of {} checks passed ({:.1} s)", t.ok, t.total, t0.elapsed().as_secs_f64());
    if t.ok != t.total {
        std::process::exit(1);
    }
}
