#![allow(dead_code)]
// Tests of src/expressao.rs: reads per gene (count_reads) and differential
// expression (the DESeq2 method).
//
// Every result is compared with something this file computes on its own:
//  - counting: the generator records, for every read, the reference bases
//    its CIGAR aligns while it builds the CIGAR; the brute force paints the
//    exons of each gene on an array of positions and applies the documented
//    rules to the structured records, never to the SAM text. A small
//    hand-made file has its expected counts worked out by hand;
//  - special functions: published values and closed forms;
//  - size factors, dispersions, the GLM: hand-worked cases and an
//    independent likelihood (lgamma differences as sums of logs, determinant
//    by elimination) maximized on a fine grid;
//  - the whole method: negative-binomial counts simulated with known fold
//    changes (recovery, false discoveries, calibration under the null).
//
// Usage: cargo run --release --bin expressao
//        cargo run --release --bin expressao -- airway <counts> <samples> <gene_info> <out.tsv> [prior_variance]
//        cargo run --release --bin expressao -- de <counts> <samples> <condition> <reference> <factors|-> <out_prefix> [names|-] [prior_variance|-] [treated|-]
//        cargo run --release --bin expressao -- dump <dir>          (the simulated data sets, for other tools)
//        cargo run --release --bin expressao -- gff <annotation> <a sequence name> [attribute]
//        cargo run --release --bin expressao -- bench <fragments> <dir>
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
#[path = "../../../src/expressao.rs"]
mod expressao;

use expressao::*;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::io::Write;
use std::path::PathBuf;

// ---------------------------------------------------------------- harness

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
    /// |got - want| <= tol * max(|want|, piso)
    fn perto(&mut self, got: f64, want: f64, tol: f64, what: &str) -> bool {
        let ok = (got.is_nan() && want.is_nan()) || (got - want).abs() <= tol * want.abs().max(1e-300);
        if !ok {
            println!("FAIL  {what}\n   got:  {got:e}\n   want: {want:e} (relative tolerance {tol:e})");
        }
        self.check(ok, what)
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

// ---------------------------------------------------------------- random numbers

struct Rng(u64);

impl Rng {
    fn u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn unif(&mut self) -> f64 {
        (self.u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn below(&mut self, n: u64) -> u64 {
        self.u64() % n.max(1)
    }
    fn chance(&mut self, p: f64) -> bool {
        self.unif() < p
    }
    fn normal(&mut self) -> f64 {
        let u1 = 1.0 - self.unif();
        let u2 = self.unif();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
    /// Gamma(shape, scale): Marsaglia and Tsang.
    fn gamma(&mut self, shape: f64, scale: f64) -> f64 {
        if shape < 1.0 {
            let u = 1.0 - self.unif();
            return self.gamma(shape + 1.0, scale) * u.powf(1.0 / shape);
        }
        let d = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let x = self.normal();
            let v = (1.0 + c * x).powi(3);
            if v <= 0.0 {
                continue;
            }
            let u = 1.0 - self.unif();
            if u < 1.0 - 0.0331 * x.powi(4) || u.ln() < 0.5 * x * x + d * (1.0 - v + v.ln()) {
                return d * v * scale;
            }
        }
    }
    /// Poisson: multiplication for small means, PTRS (Hoermann) otherwise.
    fn poisson(&mut self, lambda: f64) -> u64 {
        if lambda <= 0.0 {
            return 0;
        }
        if lambda < 10.0 {
            let l = (-lambda).exp();
            let (mut k, mut p) = (0u64, 1.0);
            loop {
                p *= self.unif();
                if p <= l {
                    return k;
                }
                k += 1;
            }
        }
        let slam = lambda.sqrt();
        let loglam = lambda.ln();
        let b = 0.931 + 2.53 * slam;
        let a = -0.059 + 0.02483 * b;
        let invalpha = 1.1239 + 1.1328 / (b - 3.4);
        let vr = 0.9277 - 3.6224 / (b - 2.0);
        loop {
            let u = self.unif() - 0.5;
            let v = self.unif();
            let us = 0.5 - u.abs();
            let k = ((2.0 * a / us + b) * u + lambda + 0.43).floor();
            if us >= 0.07 && v <= vr {
                return k as u64;
            }
            if k < 0.0 || (us < 0.013 && v > us) {
                continue;
            }
            if v.ln() + invalpha.ln() - (a / (us * us) + b).ln() <= -lambda + k * loglam - ln_gamma(k + 1.0) {
                return k as u64;
            }
        }
    }
    /// Negative binomial with mean mu and dispersion alpha (Gamma-Poisson).
    fn nb(&mut self, mu: f64, alpha: f64) -> u64 {
        if alpha <= 0.0 {
            return self.poisson(mu);
        }
        let l = self.gamma(1.0 / alpha, mu * alpha);
        self.poisson(l)
    }
}

// ---------------------------------------------------------------- a JSON checker

fn json_valido(s: &str) -> bool {
    fn ws(b: &[u8], i: &mut usize) {
        while *i < b.len() && matches!(b[*i], b' ' | b'\n' | b'\t' | b'\r') {
            *i += 1;
        }
    }
    fn texto(b: &[u8], i: &mut usize) -> bool {
        *i += 1;
        while *i < b.len() {
            match b[*i] {
                b'"' => {
                    *i += 1;
                    return true;
                }
                b'\\' => *i += 2,
                c if c < 0x20 => return false,
                _ => *i += 1,
            }
        }
        false
    }
    fn valor(b: &[u8], i: &mut usize) -> bool {
        ws(b, i);
        if *i >= b.len() {
            return false;
        }
        match b[*i] {
            b'{' => {
                *i += 1;
                ws(b, i);
                if b.get(*i) == Some(&b'}') {
                    *i += 1;
                    return true;
                }
                loop {
                    ws(b, i);
                    if b.get(*i) != Some(&b'"') || !texto(b, i) {
                        return false;
                    }
                    ws(b, i);
                    if b.get(*i) != Some(&b':') {
                        return false;
                    }
                    *i += 1;
                    if !valor(b, i) {
                        return false;
                    }
                    ws(b, i);
                    match b.get(*i) {
                        Some(b',') => *i += 1,
                        Some(b'}') => {
                            *i += 1;
                            return true;
                        }
                        _ => return false,
                    }
                }
            }
            b'[' => {
                *i += 1;
                ws(b, i);
                if b.get(*i) == Some(&b']') {
                    *i += 1;
                    return true;
                }
                loop {
                    if !valor(b, i) {
                        return false;
                    }
                    ws(b, i);
                    match b.get(*i) {
                        Some(b',') => *i += 1,
                        Some(b']') => {
                            *i += 1;
                            return true;
                        }
                        _ => return false,
                    }
                }
            }
            b'"' => texto(b, i),
            b't' => {
                *i += 4;
                &b[*i - 4..*i] == b"true"
            }
            b'f' => {
                *i += 5;
                &b[*i - 5..*i] == b"false"
            }
            b'n' => {
                *i += 4;
                &b[*i - 4..*i] == b"null"
            }
            _ => {
                let ini = *i;
                while *i < b.len() && matches!(b[*i], b'-' | b'+' | b'.' | b'e' | b'E' | b'0'..=b'9') {
                    *i += 1;
                }
                std::str::from_utf8(&b[ini..*i]).ok().and_then(|t| t.parse::<f64>().ok()).is_some()
            }
        }
    }
    let b = s.as_bytes();
    let mut i = 0;
    let ok = valor(b, &mut i);
    ws(b, &mut i);
    ok && i == b.len()
}

// =====================================================================
//  COUNTING: a synthetic genome, annotation and alignments
// =====================================================================

/// Names in the annotation and in the alignments (chr1 vs 1 on purpose).
const CROM_ANOT: [&str; 2] = ["1", "2"];
const CROM_BAM: [&str; 4] = ["chr1", "chr2", "chr3", "chrUn_KI270742v1"];
const COMPS: [u64; 4] = [60_000, 40_000, 20_000, 5_000];

#[derive(Clone, Debug)]
struct GeneS {
    nome: String,
    geneid: u64,
    crom: usize,
    fita: u8,
    transcritos: Vec<Vec<(u64, u64)>>,
}

struct Anot {
    genes: Vec<GeneS>,
    /// Per chromosome, per position: the (gene, strand) whose exons cover it.
    pint: Vec<Vec<Vec<(u32, u8)>>>,
}

fn gera_anotacao(rng: &mut Rng) -> Anot {
    let mut genes: Vec<GeneS> = Vec::new();
    for crom in 0..2 {
        let comp = COMPS[crom];
        let mut pos = 500;
        let mut ultimo: Option<(u64, u64)> = None;
        while pos + 9000 < comp {
            let span = 1500 + rng.below(6000);
            let mut ini = pos;
            if let Some((a, b)) = ultimo {
                if rng.chance(0.35) {
                    ini = a + rng.below(b - a);
                }
            }
            let fim = (ini + span).min(comp - 100);
            let fita = match rng.below(10) {
                0 => b'.',
                1..=4 => b'-',
                _ => b'+',
            };
            let mut exons = Vec::new();
            let mut r = ini;
            while r + 50 < fim {
                let e = (r + 50 + rng.below(350)).min(fim);
                exons.push((r, e));
                r = e + 100 + rng.below(1400);
            }
            if exons.is_empty() {
                exons.push((ini, ini + 50));
            }
            let mut ts = vec![exons.clone()];
            for _ in 1..1 + rng.below(3) {
                let mut t: Vec<(u64, u64)> = exons.iter().copied().filter(|_| rng.chance(0.7)).collect();
                if t.is_empty() {
                    t.push(exons[rng.below(exons.len() as u64) as usize]);
                }
                let u = t.len() - 1;
                let (a, b) = t[u];
                let nb = if rng.chance(0.5) { b + rng.below(30) } else { b.saturating_sub(rng.below(30)).max(a + 10) };
                t[u] = (a, nb.min(comp - 1));
                ts.push(t);
            }
            let g = genes.len();
            genes.push(GeneS { nome: format!("G{g}"), geneid: 1000 + g as u64, crom, fita, transcritos: ts });
            ultimo = Some((ini, fim));
            pos = fim + 200 + rng.below(3000);
        }
    }
    let mut pint: Vec<Vec<Vec<(u32, u8)>>> = COMPS.iter().map(|&c| vec![Vec::new(); c as usize]).collect();
    for (g, gene) in genes.iter().enumerate() {
        for t in &gene.transcritos {
            for &(a, b) in t {
                for x in a..b {
                    let v = &mut pint[gene.crom][x as usize];
                    if !v.contains(&(g as u32, gene.fita)) {
                        v.push((g as u32, gene.fita));
                    }
                }
            }
        }
    }
    Anot { genes, pint }
}

fn fita_txt(f: u8) -> char {
    f as char
}

/// GFF3 like Ensembl/RefSeq: gene -> mRNA -> exon through Parent, with CDS
/// lines to ignore. `refseq`: sequences named by accession, with `region`
/// lines, and the GeneID only on the gene lines.
fn escreve_gff(a: &Anot, refseq: bool) -> Vec<String> {
    let seq = |c: usize| -> String {
        if refseq {
            ["NC_000001.11", "NC_000002.12"][c].to_string()
        } else {
            CROM_ANOT[c].to_string()
        }
    };
    let mut l = vec!["##gff-version 3".to_string()];
    for c in 0..2 {
        l.push(format!("##sequence-region {} 1 {}", seq(c), COMPS[c]));
        if refseq {
            l.push(format!(
                "{}\tRefSeq\tregion\t1\t{}\t.\t+\t.\tID={}:1..{};Dbxref=taxon:9606;Name={};chromosome={};gbkey=Src;genome=chromosome",
                seq(c),
                COMPS[c],
                seq(c),
                COMPS[c],
                CROM_ANOT[c],
                CROM_ANOT[c]
            ));
        }
    }
    for g in &a.genes {
        let ini = g.transcritos.iter().flatten().map(|e| e.0).min().unwrap();
        let fim = g.transcritos.iter().flatten().map(|e| e.1).max().unwrap();
        let f = fita_txt(g.fita);
        l.push(format!(
            "{}\ttest\tgene\t{}\t{}\t.\t{f}\t.\tID=gene-{};Name=SYM{};Dbxref=GeneID:{}",
            seq(g.crom),
            ini + 1,
            fim,
            g.nome,
            g.nome,
            g.geneid
        ));
        for (k, t) in g.transcritos.iter().enumerate() {
            let (ti, tf) = (t[0].0, t[t.len() - 1].1);
            l.push(format!("{}\ttest\tmRNA\t{}\t{}\t.\t{f}\t.\tID=rna-{}-{k};Parent=gene-{}", seq(g.crom), ti + 1, tf, g.nome, g.nome));
            for (e, &(x, y)) in t.iter().enumerate() {
                l.push(format!(
                    "{}\ttest\texon\t{}\t{}\t.\t{f}\t.\tID=exon-{}-{k}-{e};Parent=rna-{}-{k};gene=SYM{}",
                    seq(g.crom),
                    x + 1,
                    y,
                    g.nome,
                    g.nome,
                    g.nome
                ));
            }
            l.push(format!("{}\ttest\tCDS\t{}\t{}\t.\t{f}\t0\tID=cds-{}-{k};Parent=rna-{}-{k}", seq(g.crom), ti + 1, tf, g.nome, g.nome));
        }
    }
    l
}

fn escreve_gtf(a: &Anot) -> Vec<String> {
    let mut l = vec!["#!genome-build test".to_string()];
    for g in &a.genes {
        let f = fita_txt(g.fita);
        for (k, t) in g.transcritos.iter().enumerate() {
            let (ti, tf) = (t[0].0, t[t.len() - 1].1);
            l.push(format!(
                "{}\ttest\ttranscript\t{}\t{}\t.\t{f}\t.\tgene_id \"{}\"; transcript_id \"{}.t{k}\";",
                CROM_ANOT[g.crom],
                ti + 1,
                tf,
                g.nome,
                g.nome
            ));
            for &(x, y) in t {
                l.push(format!(
                    "{}\ttest\texon\t{}\t{}\t.\t{f}\t.\tgene_id \"{}\"; transcript_id \"{}.t{k}\"; gene_name \"SYM{}\";",
                    CROM_ANOT[g.crom],
                    x + 1,
                    y,
                    g.nome,
                    g.nome,
                    g.nome
                ));
            }
        }
    }
    l
}

#[derive(Clone, Debug)]
struct Rec {
    nome: String,
    flag: u16,
    crom: Option<usize>,
    pos: u64,
    cigar: String,
    blocos: Vec<(u64, u64)>,
    mapq: u8,
    nh: Option<u32>,
    mate: Option<(usize, u64)>,
}

impl Rec {
    fn sam(&self) -> String {
        let rname = self.crom.map(|c| CROM_BAM[c]).unwrap_or("*");
        let (rnext, pnext) = match (self.mate, self.crom) {
            (Some((mc, mp)), Some(c)) if mc == c => ("=".to_string(), mp),
            (Some((mc, mp)), _) => (CROM_BAM[mc].to_string(), mp),
            _ => ("*".to_string(), 0),
        };
        let nh = self.nh.map(|n| format!("\tNH:i:{n}")).unwrap_or_default();
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t*\t*{}",
            self.nome, self.flag, rname, self.pos, self.mapq, self.cigar, rnext, pnext, nh
        )
    }
}

/// A CIGAR covering `segs` (genomic pieces separated by introns, as N),
/// with clips, insertions, deletions and =/X runs; returns it with the
/// reference intervals of the aligned bases.
fn monta_cigar(rng: &mut Rng, segs: &[(u64, u64)]) -> (String, Vec<(u64, u64)>) {
    let mut ops: Vec<(u64, char)> = Vec::new();
    let mut blocos = Vec::new();
    if rng.chance(0.15) {
        ops.push((1 + rng.below(8), 'H'));
    }
    if rng.chance(0.25) {
        ops.push((1 + rng.below(12), 'S'));
    }
    let alinhada = |rng: &mut Rng| match rng.below(20) {
        0..=16 => 'M',
        17 | 18 => '=',
        _ => 'X',
    };
    for (k, &(a, b)) in segs.iter().enumerate() {
        if k > 0 {
            ops.push((a - segs[k - 1].1, 'N'));
        }
        let mut r = a;
        let mut primeiro = true;
        while r < b {
            let l = (b - r).min(4 + rng.below(50));
            let ultimo = r + l >= b;
            let op = if primeiro || ultimo {
                alinhada(rng)
            } else if rng.chance(0.15) {
                'D'
            } else {
                alinhada(rng)
            };
            if !primeiro && rng.chance(0.08) {
                ops.push((1 + rng.below(4), 'I'));
            }
            ops.push((l, op));
            if op != 'D' {
                blocos.push((r, r + l));
            }
            r += l;
            primeiro = false;
        }
    }
    if rng.chance(0.25) {
        ops.push((1 + rng.below(12), 'S'));
    }
    if rng.chance(0.15) {
        ops.push((1 + rng.below(8), 'H'));
    }
    (ops.iter().map(|(n, o)| format!("{n}{o}")).collect(), blocos)
}

/// Transcript coordinates [off, off + len) as genomic pieces.
fn no_transcrito(t: &[(u64, u64)], off: u64, len: u64) -> Vec<(u64, u64)> {
    let mut segs = Vec::new();
    let mut acc = 0;
    for &(a, b) in t {
        let l = b - a;
        let s = off.max(acc);
        let e = (off + len).min(acc + l);
        if s < e {
            segs.push((a + s - acc, a + e - acc));
        }
        acc += l;
    }
    segs
}

fn mapq(rng: &mut Rng) -> u8 {
    match rng.below(20) {
        0..=9 => 60,
        10 | 11 => 255,
        12 | 13 => 30,
        14 | 15 => 10,
        16 => 3,
        17 => 1,
        _ => 0,
    }
}

fn gera_leituras(rng: &mut Rng, a: &Anot, n: usize) -> Vec<Rec> {
    let mut recs: Vec<Rec> = Vec::new();
    let genomico = |rng: &mut Rng, len: u64| -> (usize, u64) {
        let c = match rng.below(20) {
            0..=10 => 0,
            11..=16 => 1,
            17 | 18 => 2,
            _ => 3,
        };
        (c, rng.below(COMPS[c] - len - 1))
    };
    for f in 0..n {
        let nome = format!("f{f}");
        let par = rng.chance(0.6);
        let (l1, l2) = (50 + rng.below(60), 50 + rng.below(60));
        let frag = 160 + rng.below(250);
        // where the fragment comes from
        let mut s1: Vec<(u64, u64)>;
        let mut s2: Vec<(u64, u64)>;
        let c1: usize;
        let mut c2: usize;
        let reversa: bool;
        let g = &a.genes[rng.below(a.genes.len() as u64) as usize];
        let t = &g.transcritos[rng.below(g.transcritos.len() as u64) as usize];
        let tl: u64 = t.iter().map(|e| e.1 - e.0).sum();
        if rng.chance(0.65) && tl >= frag {
            let off = rng.below(tl - frag + 1);
            s1 = no_transcrito(t, off, l1);
            s2 = no_transcrito(t, off + frag - l2, l2);
            c1 = g.crom;
            c2 = g.crom;
            let na_fita = rng.chance(0.85);
            reversa = match g.fita {
                b'+' => !na_fita,
                b'-' => na_fita,
                _ => rng.chance(0.5),
            };
        } else {
            let (c, s) = genomico(rng, frag + 10);
            s1 = vec![(s, s + l1)];
            s2 = vec![(s + frag - l2, s + frag)];
            c1 = c;
            c2 = c;
            reversa = rng.chance(0.5);
        }
        if par && rng.chance(0.03) {
            let (c, s) = genomico(rng, l2 + 10);
            c2 = c;
            s2 = vec![(s, s + l2)];
        }
        if !par && reversa && rng.chance(0.5) {
            // a single read can come from either end
            s1 = s2.clone();
        }
        let multi = rng.chance(0.06);
        let nh = if multi {
            Some(2 + rng.below(2) as u32)
        } else if rng.chance(0.5) {
            Some(1)
        } else {
            None
        };
        let mut lista: Vec<Rec> = Vec::new();
        let (cig1, b1) = monta_cigar(rng, &s1);
        let mut r1 = Rec {
            nome: nome.clone(),
            flag: 0,
            crom: Some(c1),
            pos: s1[0].0 + 1,
            cigar: cig1,
            blocos: b1,
            mapq: mapq(rng),
            nh,
            mate: None,
        };
        if !par {
            if reversa {
                r1.flag |= 0x10;
            }
            lista.push(r1);
        } else {
            let (cig2, b2) = monta_cigar(rng, &s2);
            let mut r2 = Rec {
                nome: nome.clone(),
                flag: 0x1 | 0x2 | 0x80,
                crom: Some(c2),
                pos: s2[0].0 + 1,
                cigar: cig2,
                blocos: b2,
                mapq: mapq(rng),
                nh,
                mate: None,
            };
            r1.flag |= 0x1 | 0x2 | 0x40;
            if reversa {
                r1.flag |= 0x10;
                r2.flag |= 0x20;
            } else {
                r1.flag |= 0x20;
                r2.flag |= 0x10;
            }
            r1.mate = Some((c2, r2.pos));
            r2.mate = Some((c1, r1.pos));
            lista.push(r1);
            lista.push(r2);
        }
        // unmapped
        if rng.chance(0.05) {
            if !par {
                let r = &mut lista[0];
                r.flag |= 0x4;
                r.cigar = "*".into();
                r.blocos.clear();
                if rng.chance(0.5) {
                    r.crom = None;
                    r.pos = 0;
                }
            } else if rng.chance(0.5) {
                // both
                for r in lista.iter_mut() {
                    r.flag = (r.flag & !(0x2 | 0x10 | 0x20)) | 0x4 | 0x8;
                    r.cigar = "*".into();
                    r.blocos.clear();
                    r.crom = None;
                    r.pos = 0;
                    r.mate = None;
                }
            } else {
                // one mate unmapped, placed at the other's position
                let k = rng.below(2) as usize;
                let (oc, op) = (lista[1 - k].crom, lista[1 - k].pos);
                lista[k].flag = (lista[k].flag & !0x2) | 0x4;
                lista[k].cigar = "*".into();
                lista[k].blocos.clear();
                lista[k].crom = oc;
                lista[k].pos = op;
                lista[1 - k].flag = (lista[1 - k].flag & !0x2) | 0x8;
            }
        }
        // whole-fragment flags
        if rng.chance(0.04) {
            if par && rng.chance(0.3) {
                lista[1].flag |= 0x400;
            } else {
                for r in lista.iter_mut() {
                    r.flag |= 0x400;
                }
            }
        }
        if rng.chance(0.02) {
            for r in lista.iter_mut() {
                r.flag |= 0x200;
            }
        }
        // secondary and supplementary alignments
        let base: Vec<Rec> = lista.clone();
        if multi {
            for r in &base {
                if r.flag & 0x4 != 0 {
                    continue;
                }
                let (c, s) = genomico(rng, 200);
                let (cig, bl) = monta_cigar(rng, &[(s, s + 60)]);
                let mut x = r.clone();
                x.flag = (r.flag & !(0x400 | 0x200)) | 0x100;
                x.crom = Some(c);
                x.pos = s + 1;
                x.cigar = cig;
                x.blocos = bl;
                lista.push(x);
            }
        }
        if rng.chance(0.03) && base[0].flag & 0x4 == 0 {
            let (c, s) = genomico(rng, 200);
            let (cig, bl) = monta_cigar(rng, &[(s, s + 40)]);
            let mut x = base[0].clone();
            x.flag = base[0].flag | 0x800;
            x.crom = Some(c);
            x.pos = s + 1;
            x.cigar = cig;
            x.blocos = bl;
            lista.push(x);
        }
        // a mate missing from the file
        if par && base.len() == 2 && base[1].flag & 0x4 == 0 && base[0].flag & 0x8 == 0 && rng.chance(0.03) {
            lista.remove(1);
        }
        recs.extend(lista);
    }
    recs
}

#[derive(Clone, Debug, Default, PartialEq)]
struct Esperado {
    counts: Vec<u64>,
    records: usize,
    fragments: usize,
    assigned: usize,
    ambiguous: usize,
    no_feature: usize,
    multimapping: usize,
    low_mapq: usize,
    unmapped: usize,
    secondary: usize,
    supplementary: usize,
    duplicate: usize,
    qc_fail: usize,
    pairs: usize,
    mate_missing: usize,
    off_annotation: usize,
}

fn do_resultado(s: &SampleCounts, c: &[u64]) -> Esperado {
    Esperado {
        counts: c.to_vec(),
        records: s.records,
        fragments: s.fragments,
        assigned: s.assigned,
        ambiguous: s.ambiguous,
        no_feature: s.no_feature,
        multimapping: s.multimapping,
        low_mapq: s.low_mapq,
        unmapped: s.unmapped,
        secondary: s.secondary,
        supplementary: s.supplementary,
        duplicate: s.duplicate,
        qc_fail: s.qc_fail,
        pairs: s.pairs,
        mate_missing: s.mate_missing,
        off_annotation: s.off_annotation,
    }
}

/// The documented rules, applied to the structured records.
fn forca_bruta(recs: &[Rec], a: &Anot, fita: Strandedness, min_mapq: u8) -> Esperado {
    let mut e = Esperado { counts: vec![0; a.genes.len()], ..Default::default() };
    let mut frags: Vec<Vec<&Rec>> = Vec::new();
    let mut por_nome: HashMap<&str, usize> = HashMap::new();
    for r in recs {
        e.records += 1;
        if r.flag & 0x100 != 0 {
            e.secondary += 1;
            continue;
        }
        if r.flag & 0x800 != 0 {
            e.supplementary += 1;
            continue;
        }
        if r.flag & 0x200 != 0 {
            e.qc_fail += 1;
            continue;
        }
        if r.flag & 0x400 != 0 {
            e.duplicate += 1;
            continue;
        }
        let par = r.flag & 0x1 != 0;
        if r.flag & 0x4 != 0 {
            if par && r.flag & 0x8 == 0 {
                continue;
            }
            if !par || r.flag & 0x80 == 0 {
                e.unmapped += 1;
                e.fragments += 1;
            }
            continue;
        }
        if r.crom.map_or(true, |c| c >= 2) {
            e.off_annotation += 1;
        }
        if !par || r.flag & 0x8 != 0 {
            frags.push(vec![r]);
        } else {
            match por_nome.get(r.nome.as_str()) {
                Some(&i) => frags[i].push(r),
                None => {
                    por_nome.insert(&r.nome, frags.len());
                    frags.push(vec![r]);
                }
            }
        }
    }
    for f in &frags {
        e.fragments += 1;
        let par = f[0].flag & 0x1 != 0 && f[0].flag & 0x8 == 0;
        if f.len() == 2 {
            e.pairs += 1;
        } else if par {
            e.mate_missing += 1;
        }
        if f.iter().any(|r| r.nh.map_or(false, |n| n > 1)) {
            e.multimapping += 1;
            continue;
        }
        if f.iter().any(|r| r.mapq < min_mapq) {
            e.low_mapq += 1;
            continue;
        }
        let mut genes: HashSet<u32> = HashSet::new();
        for r in f {
            let c = r.crom.unwrap();
            if c >= 2 {
                continue;
            }
            let rev = (r.flag & 0x10 != 0) ^ (r.flag & 0x1 != 0 && r.flag & 0x80 != 0);
            let frag = if rev { b'-' } else { b'+' };
            let pedida = match fita {
                Strandedness::Unstranded => None,
                Strandedness::Forward => Some(frag),
                Strandedness::Reverse => Some(if frag == b'+' { b'-' } else { b'+' }),
            };
            for &(x, y) in &r.blocos {
                for p in x..y {
                    for &(g, s) in &a.pint[c][p as usize] {
                        if pedida.map_or(true, |q| q == s || s == b'.') {
                            genes.insert(g);
                        }
                    }
                }
            }
        }
        match genes.len() {
            0 => e.no_feature += 1,
            1 => {
                e.assigned += 1;
                e.counts[*genes.iter().next().unwrap() as usize] += 1;
            }
            _ => e.ambiguous += 1,
        }
    }
    e
}

fn cabecalho_sam(ordenado: bool) -> Vec<String> {
    let mut l = vec![format!("@HD\tVN:1.6\tSO:{}", if ordenado { "coordinate" } else { "unsorted" })];
    for (n, c) in CROM_BAM.iter().zip(COMPS) {
        l.push(format!("@SQ\tSN:{n}\tLN:{c}"));
    }
    l
}

fn ordena_recs(recs: &mut [Rec]) {
    recs.sort_by_key(|r| (r.crom.unwrap_or(usize::MAX), r.pos));
}

fn testa_contagem_sintetica(t: &mut T) {
    let mut rng = Rng(20261114);
    let a = gera_anotacao(&mut rng);
    let n_genes = a.genes.len();
    println!("synthetic annotation: {} genes, {} transcripts", n_genes, a.genes.iter().map(|g| g.transcritos.len()).sum::<usize>());
    let gff = t.p("anot.gff3");
    let gtf = t.p("anot.gtf.gz");
    let refseq = t.p("anot_refseq.gff");
    grava(&gff, &escreve_gff(&a, false));
    grava(&gtf, &escreve_gtf(&a));
    grava(&refseq, &escreve_gff(&a, true));

    // three samples
    let mut amostras: Vec<Vec<Rec>> = Vec::new();
    for k in 0..3 {
        let mut r = gera_leituras(&mut rng, &a, 6000 + 1500 * k);
        ordena_recs(&mut r);
        amostras.push(r);
    }
    let mut caminhos = Vec::new();
    for (k, r) in amostras.iter().enumerate() {
        let mut l = cabecalho_sam(true);
        l.extend(r.iter().map(Rec::sam));
        let p = t.p(&format!("s{k}.sam"));
        grava(&p, &l);
        caminhos.push(p);
    }
    let esperados: Vec<Esperado> = amostras.iter().map(|r| forca_bruta(r, &a, Strandedness::Unstranded, 10)).collect();
    println!(
        "synthetic reads: {} records; sample 0: {} assigned, {} ambiguous, {} no feature, {} pairs",
        amostras.iter().map(Vec::len).sum::<usize>(),
        esperados[0].assigned,
        esperados[0].ambiguous,
        esperados[0].no_feature,
        esperados[0].pairs
    );
    t.check(esperados[0].assigned > 1000 && esperados[0].ambiguous > 50 && esperados[0].no_feature > 300, "the synthetic data exercise every category");
    t.check(
        esperados[0].multimapping > 50 && esperados[0].low_mapq > 50 && esperados[0].unmapped > 50 && esperados[0].mate_missing > 10,
        "the synthetic data exercise the filters"
    );

    let opt = CountOptions::default();
    let matriz = t.p("counts.tsv");
    let resumo = t.p("counts.summary.tsv");
    let opt_r = CountOptions { summary_output: Some(resumo.clone()), ..CountOptions::default() };
    let Some(r) = t.ok_ou(count_reads(&caminhos, &gff, Some(&matriz), &opt_r), "count_reads, GFF3, three SAM files") else { return };
    t.eq(&r.genes, &a.genes.iter().map(|g| g.nome.clone()).collect::<Vec<_>>(), "GFF3 automatic attribute: genes named by the top-level ID without 'gene-', in annotation order");
    for k in 0..3 {
        t.eq(&do_resultado(&r.samples[k], &r.counts[k]), &esperados[k], &format!("GFF3, sample {k}: counts and summary equal the brute force"));
    }
    t.eq(&r.samples.iter().map(|s| s.sample.clone()).collect::<Vec<_>>(), &vec!["s0".to_string(), "s1".into(), "s2".into()], "sample names from the file names");
    // gene lengths
    let comps: Vec<u64> = a
        .genes
        .iter()
        .map(|g| {
            let mut v: HashSet<u64> = HashSet::new();
            for tr in &g.transcritos {
                for &(x, y) in tr {
                    v.extend(x..y);
                }
            }
            v.len() as u64
        })
        .collect();
    t.eq(&r.lengths, &comps, "gene lengths are the union of their exons");
    t.eq(&r.features, &a.genes.iter().map(|g| g.transcritos.iter().map(Vec::len).sum::<usize>()).sum::<usize>(), "every exon line is a feature");
    // the matrix file
    let linhas = le(&matriz);
    t.eq(&linhas[0], &"gene\ts0\ts1\ts2".to_string(), "matrix header");
    let mut esperado_txt = Vec::new();
    for (g, gene) in a.genes.iter().enumerate() {
        esperado_txt.push(format!("{}\t{}\t{}\t{}", gene.nome, esperados[0].counts[g], esperados[1].counts[g], esperados[2].counts[g]));
    }
    t.eq(&linhas[1..].to_vec(), &esperado_txt, "matrix rows");
    let lr = le(&resumo);
    t.eq(&lr[0], &"status\ts0\ts1\ts2".to_string(), "summary header");
    t.eq(
        &lr[1],
        &format!("assigned\t{}\t{}\t{}", esperados[0].assigned, esperados[1].assigned, esperados[2].assigned),
        "summary: assigned row",
    );
    t.check(json_valido(&r.json()), "GeneCounts::json is valid JSON");
    t.check(r.resumo().contains("assigned"), "GeneCounts::resumo");

    // GTF: the same counts
    if let Some(r2) = t.ok_ou(count_reads(&caminhos, &gtf, None, &opt), "count_reads, GTF") {
        t.eq(&r2.counts, &r.counts, "GTF (gene_id) gives the same counts as the GFF3");
        t.eq(&r2.genes, &r.genes, "GTF: same genes, same order");
        t.eq(&r2.lengths, &r.lengths, "GTF: same lengths");
    }
    // explicit attributes on the GFF3
    let o = CountOptions { attribute: Some("gene".into()), ..CountOptions::default() };
    if let Some(r2) = t.ok_ou(count_reads(&caminhos[..1], &gff, None, &o), "attribute gene") {
        t.eq(&r2.genes, &a.genes.iter().map(|g| format!("SYM{}", g.nome)).collect::<Vec<_>>(), "attribute \"gene\" read on the exons");
        t.eq(&r2.counts[0], &esperados[0].counts, "attribute \"gene\": same counts");
    }
    // RefSeq-style: accessions, region lines, GeneID through the parents
    let o = CountOptions { attribute: Some("GeneID".into()), ..CountOptions::default() };
    if let Some(r2) = t.ok_ou(count_reads(&caminhos[..1], &refseq, None, &o), "RefSeq accessions and GeneID") {
        t.eq(&r2.genes, &a.genes.iter().map(|g| g.geneid.to_string()).collect::<Vec<_>>(), "GeneID from Dbxref on the gene lines, inherited by the exons");
        t.eq(&r2.counts[0], &esperados[0].counts, "NC_000001.11 matches chr1: same counts");
    }
    // BAM, gzip, shuffled order
    let bam = t.p("s0.bam");
    {
        let mut l = cabecalho_sam(true);
        l.extend(amostras[0].iter().map(Rec::sam));
        grava(&bam, &l);
    }
    if let Some(r2) = t.ok_ou(count_reads(&[bam.clone()], &gff, None, &opt), "BAM input") {
        t.eq(&do_resultado(&r2.samples[0], &r2.counts[0]), &esperados[0], "BAM: same as the SAM and the brute force");
    }
    let mut emb = amostras[0].clone();
    for i in (1..emb.len()).rev() {
        let j = rng.below(i as u64 + 1) as usize;
        emb.swap(i, j);
    }
    let gz = t.p("s0_shuffled.sam.gz");
    {
        let mut l = cabecalho_sam(false);
        l.extend(emb.iter().map(Rec::sam));
        grava(&gz, &l);
    }
    if let Some(r2) = t.ok_ou(count_reads(&[gz.clone()], &gff, None, &opt), "shuffled .sam.gz") {
        t.eq(&do_resultado(&r2.samples[0], &r2.counts[0]), &esperados[0], "any record order gives the same counts (mates matched by name)");
    }
    // strandedness and mapping quality
    for (fita, nome) in [(Strandedness::Forward, "forward"), (Strandedness::Reverse, "reverse")] {
        let o = CountOptions { strandedness: fita, ..CountOptions::default() };
        if let Some(r2) = t.ok_ou(count_reads(&caminhos[..1], &gtf, None, &o), nome) {
            let e = forca_bruta(&amostras[0], &a, fita, 10);
            t.eq(&do_resultado(&r2.samples[0], &r2.counts[0]), &e, &format!("{nome} strandedness equals the brute force"));
            if fita == Strandedness::Forward {
                t.check(e.assigned > esperados[0].assigned / 2, "the simulated forward protocol assigns most reads in forward mode");
            } else {
                t.check(e.assigned < esperados[0].assigned / 3, "few reads are assigned in reverse mode");
            }
        }
    }
    for q in [0u8, 30] {
        let o = CountOptions { min_mapq: q, ..CountOptions::default() };
        if let Some(r2) = t.ok_ou(count_reads(&caminhos[1..2], &gff, None, &o), "min_mapq") {
            let e = forca_bruta(&amostras[1], &a, Strandedness::Unstranded, q);
            t.eq(&do_resultado(&r2.samples[0], &r2.counts[0]), &e, &format!("min_mapq {q} equals the brute force"));
        }
    }
    // sample names, errors
    let o = CountOptions { sample_names: Some(vec!["A".into(), "B".into()]), ..CountOptions::default() };
    if let Some(r2) = t.ok_ou(count_reads(&caminhos[..2], &gff, None, &o), "sample names") {
        t.eq(&r2.samples[1].sample, &"B".to_string(), "sample names given");
    }
    t.erro(count_reads(&caminhos[..2], &gff, None, &CountOptions { sample_names: Some(vec!["A".into()]), ..CountOptions::default() }), &["1 sample name(s) for 2"], "one name for two files");
    t.erro(count_reads(&caminhos[..2], &gff, None, &CountOptions { sample_names: Some(vec!["A".into(), "A".into()]), ..CountOptions::default() }), &["same sample name"], "duplicate sample names");
    t.erro(count_reads(&[caminhos[0].clone(), caminhos[0].clone()], &gff, None, &opt), &["same sample name"], "the same file twice");
    t.erro(count_reads(&caminhos[..1], &gff, Some(&gff), &opt), &["itself"], "output over the annotation");
    t.erro(count_reads(&caminhos[..1], &t.p("nada.gtf"), None, &opt), &["Error opening", "nada.gtf"], "missing annotation");
    t.erro(count_reads(&[t.p("nada.bam")], &gff, None, &opt), &["Error opening", "nada.bam"], "missing alignments");
    t.erro(count_reads(&[], &gff, None, &opt), &["No alignment file"], "no file");
    t.erro(count_reads(&caminhos[..1], &gff, None, &CountOptions { feature_type: "CDSX".into(), ..CountOptions::default() }), &["no feature of type \"CDSX\"", "Types present", "exon ("], "unknown feature type");
    t.erro(count_reads(&caminhos[..1], &gtf, None, &CountOptions { attribute: Some("gene_biotype".into()), ..CountOptions::default() }), &["has a gene", "gene_biotype"], "attribute absent everywhere");
}

/// A small file whose counts are worked out by hand.
fn testa_contagem_a_mao(t: &mut T) {
    let gtf = t.p("mao.gtf");
    grava(
        &gtf,
        &[
            "1\tt\texon\t1001\t1200\t.\t+\t.\tgene_id \"A\";".to_string(),
            "1\tt\texon\t2001\t2200\t.\t+\t.\tgene_id \"A\";".into(),
            "1\tt\texon\t3001\t3200\t.\t+\t.\tgene_id \"A\";".into(),
            "1\tt\texon\t2101\t2300\t.\t-\t.\tgene_id \"B\";".into(),
            "1\tt\texon\t1501\t1600\t.\t+\t.\tgene_id \"C\";".into(),
            "1\tt\texon\t5001\t5100\t.\t.\t.\tgene_id \"D\";".into(),
            "2\tt\texon\t101\t300\t.\t+\t.\tgene_id \"E\";".into(),
            "2\tt\tgene\t101\t300\t.\t+\t.\tgene_id \"E\";".into(),
        ],
    );
    let r = |n: &str, f: u16, c: &str, p: u64, q: u8, cig: &str, rn: &str, pn: u64, tag: &str| -> String {
        format!("{n}\t{f}\t{c}\t{p}\t{q}\t{cig}\t{rn}\t{pn}\t0\t*\t*{tag}")
    };
    let mut sam = vec![
        "@HD\tVN:1.6\tSO:unsorted".to_string(),
        "@SQ\tSN:chr1\tLN:10000".into(),
        "@SQ\tSN:chr2\tLN:5000".into(),
        "@SQ\tSN:chr3\tLN:1000".into(),
    ];
    sam.extend([
        r("r1", 0, "chr1", 1051, 60, "50M", "*", 0, ""),                         // A
        r("r2", 0, "chr1", 2151, 60, "30M", "*", 0, ""),                         // A and B: ambiguous
        r("r3", 0, "chr1", 1171, 60, "30M770N20M", "*", 0, ""),                  // A; the N jumps over C
        r("r4", 0, "chr1", 1551, 60, "20M", "*", 0, "\tNH:i:1"),                 // C
        r("r5", 0, "chr1", 1491, 60, "5M200D5M", "*", 0, ""),                    // the D covers C: no feature
        r("r6", 0, "chr1", 4001, 60, "50M", "*", 0, ""),                         // intergenic
        r("r7", 0x1 | 0x2 | 0x40 | 0x20, "chr1", 1101, 60, "50M", "=", 3101, ""), // pair A + A
        r("r7", 0x1 | 0x2 | 0x80 | 0x10, "chr1", 3101, 60, "50M", "=", 1101, ""),
        r("r8", 0x1 | 0x2 | 0x40 | 0x20, "chr1", 1101, 60, "50M", "=", 1521, ""), // pair A + C: ambiguous
        r("r8", 0x1 | 0x2 | 0x80 | 0x10, "chr1", 1521, 60, "50M", "=", 1101, ""),
        r("r9", 0x1 | 0x2 | 0x40 | 0x20, "chr1", 1101, 60, "50M", "=", 4001, ""), // pair A + nothing: A
        r("r9", 0x1 | 0x2 | 0x80 | 0x10, "chr1", 4001, 60, "50M", "=", 1101, ""),
        r("r10", 0, "chr1", 1101, 60, "50M", "*", 0, "\tNH:i:2"),                // multimapping
        r("r11", 0, "chr1", 1101, 5, "50M", "*", 0, ""),                         // low quality
        r("r12", 0, "chr2", 151, 60, "50M", "*", 0, ""),                         // E
        r("r13", 0, "chr3", 100, 60, "50M", "*", 0, ""),                         // not annotated
        r("r14", 0x100, "chr1", 1101, 60, "50M", "*", 0, ""),
        r("r15", 0x800, "chr1", 1101, 60, "50M", "*", 0, ""),
        r("r16", 0x400, "chr1", 1101, 60, "50M", "*", 0, ""),
        r("r17", 0x200, "chr1", 1101, 60, "50M", "*", 0, ""),
        r("r18", 0x4, "*", 0, 0, "*", "*", 0, ""),
        r("r19", 0x1 | 0x4 | 0x8 | 0x40, "*", 0, 0, "*", "*", 0, ""),
        r("r19", 0x1 | 0x4 | 0x8 | 0x80, "*", 0, 0, "*", "*", 0, ""),
        r("r20", 0x1 | 0x8 | 0x40, "chr1", 1101, 60, "50M", "=", 1101, ""),     // mate unmapped: A
        r("r20", 0x1 | 0x4 | 0x80, "chr1", 1101, 0, "*", "=", 1101, ""),
        r("r21", 0, "chr1", 5051, 60, "20M", "*", 0, ""),                        // D, no strand
        r("r22", 0x10, "chr1", 1101, 60, "50M", "*", 0, ""),                     // reverse read on A
        r("r23", 0x1 | 0x2 | 0x40, "chr1", 3051, 60, "50M", "=", 3500, ""),      // mate missing: A
    ]);
    let p = t.p("mao.sam");
    grava(&p, &sam);
    let base = |t: &mut T, fita: Strandedness| -> Option<GeneCounts> {
        t.ok_ou(count_reads(&[p.clone()], &gtf, None, &CountOptions { strandedness: fita, ..CountOptions::default() }), "hand-made file")
    };
    if let Some(r) = base(t, Strandedness::Unstranded) {
        t.eq(&r.genes, &vec!["A".to_string(), "B".into(), "C".into(), "D".into(), "E".into()], "hand-made: genes");
        t.eq(&r.counts[0], &vec![7, 0, 1, 1, 1], "hand-made, unstranded: A 7, B 0, C 1, D 1, E 1");
        let s = &r.samples[0];
        t.eq(
            &(s.assigned, s.ambiguous, s.no_feature, s.multimapping, s.low_mapq, s.unmapped, s.fragments),
            &(10, 2, 3, 1, 1, 2, 19),
            "hand-made: assigned, ambiguous, no feature, multimapping, low quality, unmapped, fragments",
        );
        t.eq(
            &(s.secondary, s.supplementary, s.duplicate, s.qc_fail, s.pairs, s.mate_missing, s.off_annotation, s.records),
            &(1, 1, 1, 1, 3, 1, 1, 28),
            "hand-made: skipped records, pairs, mates missing, off-annotation, records",
        );
        t.eq(&r.lengths, &vec![600, 200, 100, 100, 200], "hand-made: gene lengths");
    }
    if let Some(r) = base(t, Strandedness::Forward) {
        t.eq(&r.counts[0], &vec![7, 0, 1, 1, 1], "hand-made, forward: r2 goes to A, r22 to none");
        t.eq(&(r.samples[0].ambiguous, r.samples[0].no_feature), &(1, 4), "hand-made, forward: 1 ambiguous, 4 no feature");
    }
    if let Some(r) = base(t, Strandedness::Reverse) {
        t.eq(&r.counts[0], &vec![1, 1, 0, 1, 0], "hand-made, reverse: r22 to A, r2 to B, D unstranded");
        t.eq(&(r.samples[0].assigned, r.samples[0].ambiguous, r.samples[0].no_feature), &(3, 0, 12), "hand-made, reverse: summary");
    }
    // the same genes as an Ensembl GFF3: gene_id only on the gene line, two
    // levels above the exons, which carry an exon Name
    let ens = t.p("mao_ensembl.gff3");
    let mut l = vec!["##gff-version 3".to_string()];
    for (g, sq, f, exons) in [
        ("A", "1", '+', vec![(1001, 1200), (2001, 2200), (3001, 3200)]),
        ("B", "1", '-', vec![(2101, 2300)]),
        ("C", "1", '+', vec![(1501, 1600)]),
        ("D", "1", '.', vec![(5001, 5100)]),
        ("E", "2", '+', vec![(101, 300)]),
    ] {
        let (ini, fim) = (exons[0].0, exons[exons.len() - 1].1);
        l.push(format!("{sq}\tens\tgene\t{ini}\t{fim}\t.\t{f}\t.\tID=gene:ENSG{g};gene_id=ENSG{g};Name=SYM{g}"));
        l.push(format!("{sq}\tens\tmRNA\t{ini}\t{fim}\t.\t{f}\t.\tID=transcript:ENST{g};Parent=gene:ENSG{g};transcript_id=ENST{g}"));
        for (k, (a, b)) in exons.iter().enumerate() {
            l.push(format!("{sq}\tens\texon\t{a}\t{b}\t.\t{f}\t.\tParent=transcript:ENST{g};Name=ENSE{g}{k}"));
            l.push(format!("{sq}\tens\tCDS\t{a}\t{b}\t.\t{f}\t0\tID=CDS:ENSP{g};Parent=transcript:ENST{g}"));
        }
    }
    grava(&ens, &l);
    if let Some(r) = t.ok_ou(count_reads(&[p.clone()], &ens, None, &CountOptions::default()), "Ensembl GFF3") {
        t.eq(&r.genes, &vec!["ENSGA".to_string(), "ENSGB".into(), "ENSGC".into(), "ENSGD".into(), "ENSGE".into()], "Ensembl GFF3: gene_id found two levels up");
        t.eq(&r.counts[0], &vec![7, 0, 1, 1, 1], "Ensembl GFF3: the same counts as the GTF");
    }
    if let Some(r) = t.ok_ou(count_reads(&[p.clone()], &ens, None, &CountOptions { attribute: Some("Name".into()), ..CountOptions::default() }), "attribute Name") {
        t.eq(&r.genes[..2].to_vec(), &vec!["ENSEA0".to_string(), "ENSEA1".into()], "an explicit attribute is read on the feature first (here the exon's own Name)");
    }

    // as BAM
    let pb = t.p("mao.bam");
    grava(&pb, &sam);
    if let Some(r) = t.ok_ou(count_reads(&[pb], &gtf, None, &CountOptions::default()), "hand-made BAM") {
        t.eq(&r.counts[0], &vec![7, 0, 1, 1, 1], "hand-made BAM: same counts");
    }

    // errors in the inputs
    let ruim = t.p("ruim.gtf");
    grava(&ruim, &["1\tt\texon\t1001\t1200\t.\t+\t.\tgene_id \"A\";".to_string(), "1\tt\texon\t2001".into()]);
    t.erro(count_reads(&[p.clone()], &ruim, None, &CountOptions::default()), &["ruim.gtf, line 2", "9"], "annotation line with 4 columns");
    let ruim = t.p("ruim2.gtf");
    grava(&ruim, &["1\tt\texon\t1001\tx\t.\t+\t.\tgene_id \"A\";".to_string()]);
    t.erro(count_reads(&[p.clone()], &ruim, None, &CountOptions::default()), &["line 1", "whole numbers"], "coordinate not a number");
    let ruim = t.p("ruim3.gtf");
    grava(&ruim, &["1\tt\texon\t1001\t1200\t.\t+\t.\ttranscript_id \"x\";".to_string()]);
    t.erro(count_reads(&[p.clone()], &ruim, None, &CountOptions { attribute: Some("gene_id".into()), ..CountOptions::default() }), &["has a gene", "gene_id", "line 1"], "exon without gene_id");
    let outro = t.p("outro.sam");
    grava(&outro, &["@SQ\tSN:contig7\tLN:1000".to_string(), "x\t0\tcontig7\t1\t60\t10M\t*\t0\t0\t*\t*".into()]);
    t.erro(count_reads(&[outro], &gtf, None, &CountOptions::default()), &["None of the 1 sequences", "contig7", "same genome"], "sequence names that match nothing");
    let curto = t.p("curto.sam");
    grava(&curto, &["@SQ\tSN:chr1\tLN:10000".to_string(), "x\t0\tchr1\t1\t60".into()]);
    t.erro(count_reads(&[curto], &gtf, None, &CountOptions::default()), &["curto.sam, line 2", "at least 11"], "SAM line with 5 columns");
    let cig = t.p("cigar.sam");
    grava(&cig, &["@SQ\tSN:chr1\tLN:10000".to_string(), "x\t0\tchr1\t1\t60\t10Z\t*\t0\t0\t*\t*".into()]);
    t.erro(count_reads(&[cig], &gtf, None, &CountOptions::default()), &["line 2", "invalid CIGAR \"10Z\""], "invalid CIGAR");
    t.erro(Strandedness::from_text("sideways"), &["Unknown strandedness", "reverse"], "unknown strandedness");
    t.eq(&Strandedness::from_text("2").ok(), &Some(Strandedness::Reverse), "strandedness 2 is reverse");
    t.eq(&Strandedness::from_text("yes").ok(), &Some(Strandedness::Forward), "strandedness yes is forward");
}

// =====================================================================
//  SPECIAL FUNCTIONS AND SMALL PIECES
// =====================================================================

fn testa_funcoes(t: &mut T) {
    t.perto(ln_gamma(0.5), 0.5723649429247001, 1e-14, "lnGamma(1/2) = ln sqrt(pi)");
    t.check(ln_gamma(1.0).abs() < 1e-14 && ln_gamma(2.0).abs() < 1e-14, "lnGamma(1) = lnGamma(2) = 0");
    t.perto(ln_gamma(10.0), 12.801827480081469, 1e-14, "lnGamma(10) = ln 9!");
    t.perto(ln_gamma(101.0), 363.73937555556347, 1e-14, "lnGamma(101) = ln 100!");
    t.perto(ln_gamma(1e-5), 11.512919692895824, 1e-12, "lnGamma(1e-5)");
    t.perto(ln_gamma(1e10), 220258509288.81058, 1e-14, "lnGamma(1e10) (mpmath)");
    t.perto(digamma(1.0), -0.5772156649015329, 1e-14, "digamma(1) = -Euler's gamma");
    t.perto(digamma(0.5), -1.9635100260214235, 1e-14, "digamma(1/2)");
    t.perto(digamma(10.0), 2.251752589066721, 1e-14, "digamma(10)");
    t.perto(trigamma(1.0), 1.6449340668482264, 1e-14, "trigamma(1) = pi^2/6");
    t.perto(trigamma(0.5), 4.934802200544679, 1e-14, "trigamma(1/2) = pi^2/2");
    t.perto(trigamma(1.5), 0.9348022005446793, 1e-13, "trigamma(3/2) = pi^2/2 - 4");
    t.perto(normal_cdf(0.0), 0.5, 1e-15, "Phi(0)");
    t.perto(normal_cdf(1.0), 0.8413447460685429, 1e-15, "Phi(1)");
    t.perto(normal_cdf(-1.96), 0.024997895148220435, 1e-14, "Phi(-1.96)");
    t.perto(normal_cdf(-5.0), 2.866515718791939e-07, 1e-13, "Phi(-5)");
    t.perto(normal_cdf(-10.0), 7.619853024160527e-24, 1e-13, "Phi(-10)");
    t.perto(normal_cdf(-20.0), 2.7536241186062337e-89, 1e-13, "Phi(-20)");
    t.perto(normal_cdf(-35.0), 1.1249107064724062e-268, 1e-13, "Phi(-35) (mpmath)");
    t.perto(two_sided_p(1.959963984540054), 0.05, 1e-13, "two-sided p of 1.96 is 0.05");
    t.perto(two_sided_p(-3.0), 0.0026997960632601866, 1e-13, "two-sided p of -3");
    t.eq(&two_sided_p(0.0), &1.0, "two-sided p of 0 is 1");
    for d2 in [1.0, 2.0, 3.0, 6.0, 10.0, 25.0] {
        let exato = d2 / 2.0 * (0.01f64.powf(-2.0 / d2) - 1.0);
        t.perto(f_quantile(0.99, 2.0, d2), exato, 1e-10, &format!("qf(0.99, 2, {d2}) closed form"));
    }
    t.perto(f_quantile(0.99, 5.0, 3.0) * f_quantile(0.01, 3.0, 5.0), 1.0, 1e-9, "qf(q, a, b) = 1 / qf(1 - q, b, a)");
    t.perto(f_quantile(0.99, 1.0, 10.0), 10.04428927339658, 1e-9, "qf(0.99, 1, 10) = qt(0.995, 10)^2");
    t.perto(f_quantile(0.99, 5.0, 3.0), 28.237080837755045, 1e-9, "qf(0.99, 5, 3) (scipy)");
    t.perto(f_cdf(f_quantile(0.95, 4.0, 7.0), 4.0, 7.0), 0.95, 1e-12, "F CDF inverts the quantile");
    t.eq(&median(&[3.0, 1.0, 2.0]), &2.0, "median of 3");
    t.eq(&median(&[4.0, 1.0, 3.0, 2.0]), &2.5, "median of 4 averages the middle two");
    t.perto(mad(&[1.0, 2.0, 3.0, 4.0, 100.0]), 1.4826, 1e-15, "mad as R: 1.4826 x median absolute deviation");
    // BH, by hand and against the definition
    let p = [0.01, 0.04, 0.03, 0.005, f64::NAN];
    let a = benjamini_hochberg(&p);
    t.perto(a[0], 0.02, 1e-15, "BH by hand (1)");
    t.perto(a[1], 0.04, 1e-15, "BH by hand (2)");
    t.perto(a[2], 0.04, 1e-15, "BH by hand (3)");
    t.perto(a[3], 0.02, 1e-15, "BH by hand (4)");
    t.check(a[4].is_nan(), "BH keeps NA and leaves it out of n");
    let mut rng = Rng(7);
    let p: Vec<f64> = (0..300).map(|i| if i % 37 == 0 { f64::NAN } else if i % 11 == 0 { 0.25 } else { rng.unif().powi(3) }).collect();
    let a = benjamini_hochberg(&p);
    let validos: Vec<f64> = p.iter().copied().filter(|x| !x.is_nan()).collect();
    let n = validos.len() as f64;
    let mut ok = true;
    for (i, &pi) in p.iter().enumerate() {
        if pi.is_nan() {
            ok &= a[i].is_nan();
            continue;
        }
        // min over p_j >= p_i of n p_j / rank_j, rank = number of p <= p_j
        let mut m = f64::INFINITY;
        for &pj in &validos {
            if pj >= pi {
                let rank = validos.iter().filter(|&&x| x <= pj).count() as f64;
                m = m.min(n * pj / rank);
            }
        }
        ok &= (a[i] - m.min(1.0)).abs() < 1e-12;
    }
    t.check(ok, "BH equals its definition, with ties and NA");
    // lowess
    let x: Vec<f64> = (0..50).map(|i| i as f64 / 7.0).collect();
    let y: Vec<f64> = x.iter().map(|v| 2.0 * v + 1.0).collect();
    let f = lowess(&x, &y, 0.2);
    t.check(f.iter().zip(&y).all(|(a, b)| (a - b).abs() < 1e-10), "lowess reproduces a straight line");
    // R's lowess stops the robustness iterations when the median residual is
    // zero, so the outlier test needs noise
    let mut rng = Rng(31);
    let limpo: Vec<f64> = y.iter().map(|v| v + 0.05 * rng.normal()).collect();
    let mut y2 = limpo.clone();
    y2[20] += 50.0;
    let f = lowess(&x, &y2, 0.4);
    let f0 = lowess(&x, &limpo, 0.4);
    let pior = f.iter().zip(&f0).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    t.check(pior < 0.05, &format!("lowess's robustness iterations ignore one outlier (largest change {pior:.4})"));
    // with a 10-point window the first pass spreads the outlier over every
    // neighbour, all robustness weights of the window become 0 and the raw
    // value is kept: what R's clowess and statsmodels do too
    let f = lowess(&x, &y2, 0.2);
    t.check((f[20] - y2[20]).abs() < 1e-9, "lowess with a 10-point window keeps the outlier, as R and statsmodels do");
    let x4: Vec<f64> = (0..200).map(|i| i as f64 / 40.0).collect();
    let y3: Vec<f64> = x4.iter().map(|v| (v * 1.3).sin()).collect();
    let f = lowess(&x4, &y3, 0.1);
    t.check(f.iter().zip(&y3).all(|(a, b)| (a - b).abs() < 0.01), "lowess follows a smooth curve");
    lowess_cars(t);
    t.eq(&num(f64::NAN), &"NA".to_string(), "NA written as R does");
    t.eq(&num(1e-300), &"1e-300".to_string(), "tiny numbers in scientific notation");
    t.eq(&num(0.25), &"0.25".to_string(), "ordinary numbers as they are");
}

/// R's `lowess(cars)` (f = 2/3, iter = 3), the example of R's help page.
fn lowess_cars(t: &mut T) {
    let speed = [
        4., 4., 7., 7., 8., 9., 10., 10., 10., 11., 11., 12., 12., 12., 12., 13., 13., 13., 13., 14., 14., 14., 14., 15., 15., 15., 16., 16., 17., 17., 17., 18., 18., 18., 18., 19., 19., 19., 20., 20., 20., 20., 20., 22., 23., 24., 24., 24., 24., 25.,
    ];
    let dist = [
        2., 10., 4., 22., 16., 10., 18., 26., 34., 17., 28., 14., 20., 24., 28., 26., 34., 34., 46., 26., 36., 60., 80., 20., 26., 54., 32., 40., 32., 40., 50., 42., 56., 76., 84., 36., 46., 68., 32., 48., 52., 56., 64., 66., 54., 70., 92., 93., 120., 85.,
    ];
    let r_y = [
        4.965459, 4.965459, 13.124495, 13.124495, 15.858633, 18.579691, 21.280313, 21.280313, 21.280313, 24.129277, 24.129277, 27.119549, 27.119549, 27.119549, 27.119549, 30.027276, 30.027276, 30.027276, 30.027276, 32.962506, 32.962506,
        32.962506, 32.962506, 36.757728, 36.757728, 36.757728, 40.435075, 40.435075, 43.463492, 43.463492, 43.463492, 46.885479, 46.885479, 46.885479, 46.885479, 50.793152, 50.793152, 50.793152, 56.491224, 56.491224, 56.491224, 56.491224,
        56.491224, 67.585824, 73.079695, 78.643164, 78.643164, 78.643164, 78.643164, 84.328698,
    ];
    let f = lowess(&speed, &dist, 2.0 / 3.0);
    let pior = f.iter().zip(&r_y).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    if pior > 1e-5 {
        println!("   lowess(cars): {f:?}");
    }
    t.check(pior < 1e-5, "lowess(cars) equals R's output (6 decimals)");
}

// =====================================================================
//  DESeq2 STEPS
// =====================================================================

/// ln Gamma(y + s) - ln Gamma(s) as a sum of logs (y a whole number).
fn lgd(y: u64, s: f64) -> f64 {
    (0..y).map(|k| (s + k as f64).ln()).sum()
}

fn det(mut a: Vec<f64>, p: usize) -> f64 {
    let mut d = 1.0;
    for c in 0..p {
        let piv = (c..p).max_by(|&i, &j| a[i * p + c].abs().total_cmp(&a[j * p + c].abs())).unwrap();
        if a[piv * p + c] == 0.0 {
            return 0.0;
        }
        if piv != c {
            for k in 0..p {
                a.swap(piv * p + k, c * p + k);
            }
            d = -d;
        }
        d *= a[c * p + c];
        for i in c + 1..p {
            let f = a[i * p + c] / a[c * p + c];
            for k in c..p {
                a[i * p + k] -= f * a[c * p + k];
            }
        }
    }
    d
}

/// The Cox-Reid adjusted log likelihood of a dispersion (plus the prior),
/// written independently of the engine.
fn lp_indep(alfa: f64, y: &[u64], mu: &[f64], x: &[f64], p: usize, prior: Option<(f64, f64)>) -> f64 {
    let s = 1.0 / alfa;
    let mut ll = 0.0;
    let mut b = vec![0.0; p * p];
    for j in 0..y.len() {
        ll += lgd(y[j], s) - y[j] as f64 * (mu[j] + s).ln() - s * (1.0 + mu[j] * alfa).ln();
        let w = mu[j] / (1.0 + mu[j] * alfa);
        for r in 0..p {
            for c in 0..p {
                b[r * p + c] += w * x[j * p + r] * x[j * p + c];
            }
        }
    }
    let mut v = ll - 0.5 * det(b, p).ln();
    if let Some((m, var)) = prior {
        v -= 0.5 * (alfa.ln() - m).powi(2) / var;
    }
    v
}

fn melhor_na_grade(y: &[u64], mu: &[f64], x: &[f64], p: usize, prior: Option<(f64, f64)>) -> f64 {
    let mut best: (f64, f64) = (0.0, f64::NEG_INFINITY);
    let mut a: f64 = -14.0;
    while a <= 3.0 {
        let v = lp_indep(a.exp(), y, mu, x, p, prior);
        if v > best.1 {
            best = (a, v);
        }
        a += 0.02;
    }
    let centro = best.0;
    let mut a = centro - 0.04;
    while a <= centro + 0.04 {
        let v = lp_indep(a.exp(), y, mu, x, p, prior);
        if v > best.1 {
            best = (a, v);
        }
        a += 0.0005;
    }
    // refine
    let (mut lo, mut hi) = (best.0 - 0.001, best.0 + 0.001);
    for _ in 0..60 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if lp_indep(m1.exp(), y, mu, x, p, prior) < lp_indep(m2.exp(), y, mu, x, p, prior) {
            lo = m1;
        } else {
            hi = m2;
        }
    }
    (0.5 * (lo + hi)).exp()
}

fn testa_passos(t: &mut T) {
    // size factors by hand
    let c = [10.0, 20.0, 40.0, 8.0, 8.0, 8.0, 2.0, 6.0, 18.0];
    let sf = size_factors(&c, 3).unwrap();
    t.perto(sf[0], 0.5, 1e-14, "size factor by hand (sample 1: median of 0.5, 1, 1/3)");
    t.perto(sf[1], 1.0, 1e-14, "size factor by hand (sample 2)");
    t.perto(sf[2], 2.0, 1e-14, "size factor by hand (sample 3: median of 2, 1, 3)");
    let c = [10.0, 20.0, 40.0, 8.0, 8.0, 8.0, 2.0, 6.0, 18.0, 0.0, 5.0, 9.0, 3.0, 3.0, 12.0];
    // genes 1, 2, 3, 5 (gene 4 has a zero); sample 1 ratios: 0.5, 1, 1/3, (3/sqrt3(108))
    let lg5 = (3.0f64 * 3.0 * 12.0).ln() / 3.0;
    let r5 = 3.0f64.ln() - lg5;
    let mut v = vec![0.5f64.ln(), 0.0, (1.0f64 / 3.0).ln(), r5];
    v.sort_by(|a, b| a.total_cmp(b));
    let sf = size_factors(&c, 3).unwrap();
    t.perto(sf[0], ((v[1] + v[2]) / 2.0).exp(), 1e-13, "size factor with an even number of genes: mean of the middle two log ratios; genes with a zero left out");
    t.erro(size_factors(&[0.0, 1.0, 1.0, 0.0], 2), &["zero"], "size factors when every gene has a zero");
    let pc = size_factors_poscounts(&[0.0, 4.0, 8.0, 2.0, 0.0, 6.0], 3);
    t.check(pc.iter().map(|s| s.ln()).sum::<f64>().abs() < 1e-12, "poscounts: geometric mean 1");

    // gene-wise and MAP dispersions against the grid
    let mut rng = Rng(99);
    let casos: Vec<(usize, usize, f64, f64)> = vec![(6, 2, 0.1, 50.0), (4, 1, 0.3, 20.0), (8, 2, 1.0, 5.0), (10, 3, 0.05, 300.0), (6, 2, 0.02, 1000.0)];
    let mut n_ok = 0;
    let mut n_map = 0;
    let mut total = 0;
    for &(m, p, alfa, media) in &casos {
        for _ in 0..8 {
            let mut x = vec![0.0; m * p];
            for j in 0..m {
                x[j * p] = 1.0;
                if p >= 2 {
                    x[j * p + 1] = if j % 2 == 0 { 1.0 } else { 0.0 };
                }
                if p >= 3 {
                    x[j * p + 2] = if j % 3 == 0 { 1.0 } else { 0.0 };
                }
            }
            let mu: Vec<f64> = (0..m).map(|j| media * (1.0 + 0.3 * x[j * p + if p > 1 { 1 } else { 0 }])).collect();
            let y: Vec<u64> = mu.iter().map(|&u| rng.nb(u, alfa)).collect();
            let yf: Vec<f64> = y.iter().map(|&v| v as f64).collect();
            let best = melhor_na_grade(&y, &mu, &x, p, None);
            if best < 1e-5 {
                continue;
            }
            total += 1;
            let got = cox_reid_dispersion(&yf, &mu, &x, p, 0.5);
            let lp_g = lp_indep(got, &y, &mu, &x, p, None);
            let lp_b = lp_indep(best, &y, &mu, &x, p, None);
            if (got.ln() - best.ln()).abs() < 0.02 || lp_g >= lp_b - 1e-5 {
                n_ok += 1;
            } else {
                println!("   gene-wise dispersion: y {y:?} mu {mu:?}: got {got}, grid {best}");
            }
            let tr = alfa * 1.3;
            let bm = melhor_na_grade(&y, &mu, &x, p, Some((tr.ln(), 0.5)));
            let gm = map_dispersion(&yf, &mu, &x, p, got, tr, 0.5);
            if (gm.ln() - bm.ln()).abs() < 0.02 || lp_indep(gm, &y, &mu, &x, p, Some((tr.ln(), 0.5))) >= lp_indep(bm, &y, &mu, &x, p, Some((tr.ln(), 0.5))) - 1e-5 {
                n_map += 1;
            } else {
                println!("   MAP dispersion: y {y:?}: got {gm}, grid {bm}");
            }
        }
    }
    t.eq(&n_ok, &total, "gene-wise (Cox-Reid) dispersions equal the grid maximum of an independent likelihood");
    t.eq(&n_map, &total, "MAP dispersions equal the grid maximum of likelihood + prior");
    let y = [10.0, 10.0, 10.0, 10.0];
    let x = [1.0, 1.0, 1.0, 1.0];
    t.check(cox_reid_dispersion(&y, &[10.0; 4], &x, 1, 0.1) < 1e-4, "no extra-Poisson variation: dispersion at the minimum");

    // the trend
    let medias: Vec<f64> = (0..400).map(|i| 10f64.powf(i as f64 / 100.0)).collect();
    let disps: Vec<f64> = medias.iter().map(|m| 0.05 + 3.0 / m).collect();
    if let Some(tr) = t.ok_ou(dispersion_trend(&medias, &disps), "trend without noise") {
        t.perto(tr.asymptotic, 0.05, 1e-6, "trend without noise: asymptotic dispersion");
        t.perto(tr.extra_poisson, 3.0, 1e-6, "trend without noise: extra-Poisson term");
    }
    let medias: Vec<f64> = (0..6000).map(|_| 10f64.powf(rng.unif() * 4.0)).collect();
    let disps: Vec<f64> = medias.iter().map(|m| (0.08 + 1.5 / m) * rng.gamma(4.0, 0.25)).collect();
    if let Some(tr) = t.ok_ou(dispersion_trend(&medias, &disps), "trend with noise") {
        t.perto(tr.asymptotic, 0.08, 0.08, "trend with gamma noise: asymptotic dispersion within 8%");
        t.perto(tr.extra_poisson, 1.5, 0.08, "trend with gamma noise: extra-Poisson term within 8%");
        t.perto(tr.at(100.0), 0.08 + 1.5 / 100.0, 0.08, "trend at mean 100");
    }
}

// =====================================================================
//  THE WHOLE METHOD
// =====================================================================

struct Sim {
    counts: Vec<f64>,
    genes: Vec<String>,
    lfc: Vec<f64>,
    sf: Vec<f64>,
}

/// NB counts: `nivel[j]` the condition level of sample j (0 the reference),
/// `efeito[l][g]` the log2 fold change of level l for gene g, `bloco[j]`
/// (optional) a per-sample, per-gene log effect.
fn simula(rng: &mut Rng, n_genes: usize, nivel: &[usize], efeito: &[Vec<f64>], bloco: Option<&[Vec<f64>]>, a0: f64, a1: f64) -> Sim {
    let n = nivel.len();
    let sf: Vec<f64> = (0..n).map(|_| 0.5 + 1.5 * rng.unif()).collect();
    let mut counts = Vec::with_capacity(n * n_genes);
    let mut lfc = Vec::new();
    for g in 0..n_genes {
        let media = (5.0f64.ln() + 3.5 * rng.unif() * 2.3).exp();
        let disp = (a0 + a1 / media) * (0.5 * rng.normal()).exp();
        lfc.push(if efeito.len() > 1 { efeito[1][g] } else { 0.0 });
        for j in 0..n {
            let mut mu = sf[j] * media * 2f64.powf(efeito[nivel[j]][g]);
            if let Some(b) = bloco {
                mu *= b[j][g].exp();
            }
            counts.push(rng.nb(mu, disp) as f64);
        }
    }
    Sim { counts, genes: (0..n_genes).map(|g| format!("g{g}")).collect(), lfc, sf }
}

fn corr(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let cov: f64 = a.iter().zip(b).map(|(x, y)| (x - ma) * (y - mb)).sum();
    let va: f64 = a.iter().map(|x| (x - ma).powi(2)).sum();
    let vb: f64 = b.iter().map(|y| (y - mb).powi(2)).sum();
    cov / (va * vb).sqrt()
}

fn niveis_txt(nivel: &[usize], nomes: &[&str]) -> Vec<String> {
    nivel.iter().map(|&l| nomes[l].to_string()).collect()
}

fn amostras(n: usize) -> Vec<String> {
    (0..n).map(|j| format!("s{j}")).collect()
}

fn testa_simulacao(t: &mut T) {
    let mut rng = Rng(424242);
    let n_genes = 3000;
    // 4 vs 4, 10% DE
    let nivel = [0, 0, 0, 0, 1, 1, 1, 1];
    let ef: Vec<f64> = (0..n_genes).map(|g| if g % 10 == 0 { (1.0 + 2.0 * rng.unif()) * if rng.chance(0.5) { 1.0 } else { -1.0 } } else { 0.0 }).collect();
    let s = simula(&mut rng, n_genes, &nivel, &[vec![0.0; n_genes], ef.clone()], None, 0.05, 2.0);
    let opt = DeOptions { condition: "cond".into(), reference: "ctl".into(), ..DeOptions::default() };
    let Some(r) = t.ok_ou(deseq(&s.counts, &s.genes, &amostras(8), &niveis_txt(&nivel, &["ctl", "trt"]), &[], &opt), "simulated 4 vs 4") else { return };
    let razao: Vec<f64> = r.size_factors.iter().zip(&s.sf).map(|(a, b)| a / b).collect();
    let mr = razao.iter().sum::<f64>() / 8.0;
    t.check(razao.iter().all(|x| (x / mr - 1.0).abs() < 0.04), "size factors proportional to the true ones (within 4%)");
    t.eq(&r.treated, &"trt".to_string(), "the treated level is the only other level");
    t.eq(&r.design, &"~ cond".to_string(), "design formula");
    t.eq(&r.prior_method, &"mad".to_string(), "6 residual df: prior variance by the MAD");
    // The simulated dispersions scatter log-normally around the curve, and
    // DESeq2's fit leaves out the gene-wise estimates near zero, frequent at
    // low counts with 6 df: the 1/mean term comes out inflated, by design of
    // the method (pydeseq2 gives the same on these data).
    println!("4 vs 4: trend {:.4} + {:.3}/mean (simulated 0.05 + 2/mean)", r.trend.asymptotic, r.trend.extra_poisson);
    t.perto(r.trend.asymptotic, 0.05, 0.35, "trend: asymptotic dispersion near the simulated 0.05");
    t.check(r.trend.extra_poisson > 1.0 && r.trend.extra_poisson < 5.0, "trend: extra-Poisson term of the simulated order (2)");
    let de: Vec<usize> = (0..n_genes).filter(|&g| ef[g] != 0.0 && r.genes[g].base_mean > 20.0).collect();
    let est: Vec<f64> = de.iter().map(|&g| r.genes[g].log2_fold_change).collect();
    let ver: Vec<f64> = de.iter().map(|&g| ef[g]).collect();
    let c = corr(&est, &ver);
    println!("4 vs 4: LFC correlation {c:.3} over {} DE genes; {} up, {} down", de.len(), r.up, r.down);
    t.check(c > 0.95, "log2 fold changes recovered (correlation > 0.95 for mean > 20)");
    let sig: Vec<usize> = (0..n_genes).filter(|&g| r.genes[g].padj < 0.1).collect();
    let falsos = sig.iter().filter(|&&g| ef[g] == 0.0).count();
    let verdadeiros = sig.len() - falsos;
    let total_de = ef.iter().filter(|e| **e != 0.0).count();
    println!("4 vs 4: {} discoveries, {} false ({:.1}%), power {:.1}%", sig.len(), falsos, 100.0 * falsos as f64 / sig.len().max(1) as f64, 100.0 * verdadeiros as f64 / total_de as f64);
    t.check((falsos as f64) < 0.15 * sig.len() as f64, "false discovery proportion below 15% at padj < 0.1");
    t.check(verdadeiros as f64 > 0.6 * total_de as f64, "power above 60%");
    t.eq(&(r.up + r.down), &sig.len(), "summary up + down = discoveries");
    t.check(json_valido(&r.json()), "DeResult::json is valid JSON");
    t.check(r.resumo().contains("LFC > 0 (up)"), "DeResult::resumo");
    // the same run twice: same numbers (the threads do not change anything)
    if let Ok(r2) = deseq(&s.counts, &s.genes, &amostras(8), &niveis_txt(&nivel, &["ctl", "trt"]), &[], &opt) {
        t.check(r2.genes.iter().zip(&r.genes).all(|(a, b)| a.pvalue.to_bits() == b.pvalue.to_bits() && a.dispersion.to_bits() == b.dispersion.to_bits()), "deterministic: two runs give the same bits");
    }
    // filtering off: plain BH over every p-value
    let o2 = DeOptions { independent_filtering: false, ..opt.clone() };
    if let Ok(r2) = deseq(&s.counts, &s.genes, &amostras(8), &niveis_txt(&nivel, &["ctl", "trt"]), &[], &o2) {
        let p: Vec<f64> = r2.genes.iter().map(|g| g.pvalue).collect();
        let bh = benjamini_hochberg(&p);
        t.check(r2.genes.iter().zip(&bh).all(|(g, b)| g.padj.to_bits() == b.to_bits()), "without independent filtering padj is BH over all p-values");
        t.eq(&r2.low_counts, &0, "without filtering, no low-count gene");
    }
    // filtering on: the threshold and the subset
    let limiar = r.filter_threshold;
    let filtrados: Vec<&DeGene> = r.genes.iter().filter(|g| g.filtered).collect();
    t.check(filtrados.iter().all(|g| g.base_mean < limiar && g.padj.is_nan()), "filtered genes are below the mean-count threshold");
    let ficam: Vec<usize> = (0..n_genes).filter(|&g| r.genes[g].base_mean >= limiar).collect();
    let bh = benjamini_hochberg(&ficam.iter().map(|&g| r.genes[g].pvalue).collect::<Vec<_>>());
    t.check(ficam.iter().zip(&bh).all(|(&g, b)| r.genes[g].padj.to_bits() == b.to_bits()), "padj of the kept genes is BH over the kept genes");

    // null data: calibration
    let s0 = simula(&mut rng, n_genes, &nivel, &[vec![0.0; n_genes], vec![0.0; n_genes]], None, 0.05, 2.0);
    if let Some(r0) = t.ok_ou(deseq(&s0.counts, &s0.genes, &amostras(8), &niveis_txt(&nivel, &["ctl", "trt"]), &[], &opt), "null data") {
        let pv: Vec<f64> = r0.genes.iter().filter(|g| g.base_mean > 10.0).map(|g| g.pvalue).filter(|p| !p.is_nan()).collect();
        let f05 = pv.iter().filter(|p| **p < 0.05).count() as f64 / pv.len() as f64;
        let n_sig = r0.genes.iter().filter(|g| g.padj < 0.1).count();
        println!("null: {:.3} of p-values below 0.05 (mean > 10), {} discoveries", f05, n_sig);
        t.check(f05 > 0.02 && f05 < 0.07, "under the null about 5% of p-values fall below 0.05");
        for g in r0.genes.iter().filter(|g| g.padj < 0.1) {
            println!("   null discovery {}: baseMean {:.1}, LFC {:.2}, disp gene {:.3} trend {:.3} final {:.3}, p {:.2e}", g.gene, g.base_mean, g.log2_fold_change, g.disp_gene, g.disp_trend, g.dispersion, g.pvalue);
        }
        t.check(n_sig <= 10, "under the null few discoveries (dispersion outliers shrunk to the trend)");
    }

    // a blocking factor: 6 pairs with strong pair effects
    let nivel = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1];
    let ef: Vec<f64> = (0..n_genes).map(|g| if g % 10 == 3 { (0.6 + 1.0 * rng.unif()) * if rng.chance(0.5) { 1.0 } else { -1.0 } } else { 0.0 }).collect();
    let par_ef: Vec<Vec<f64>> = (0..6).map(|_| (0..n_genes).map(|_| 0.8 * rng.normal()).collect()).collect();
    let bloco: Vec<Vec<f64>> = (0..12).map(|j| par_ef[j / 2].clone()).collect();
    let s = simula(&mut rng, n_genes, &nivel, &[vec![0.0; n_genes], ef.clone()], Some(&bloco), 0.02, 1.0);
    let pares: Vec<String> = (0..12).map(|j| format!("p{}", j / 2)).collect();
    let cond = niveis_txt(&nivel, &["ctl", "trt"]);
    let com = deseq(&s.counts, &s.genes, &amostras(12), &cond, &[("pair".into(), pares.clone())], &opt);
    let sem = deseq(&s.counts, &s.genes, &amostras(12), &cond, &[], &opt);
    if let (Some(com), Some(sem)) = (t.ok_ou(com, "paired design"), t.ok_ou(sem, "unpaired design")) {
        let acertos = |r: &DeResult| (0..n_genes).filter(|&g| r.genes[g].padj < 0.1 && ef[g] != 0.0).count();
        let falsos = |r: &DeResult| (0..n_genes).filter(|&g| r.genes[g].padj < 0.1 && ef[g] == 0.0).count();
        println!("paired: {} true / {} false with ~ pair + cond; {} / {} with ~ cond", acertos(&com), falsos(&com), acertos(&sem), falsos(&sem));
        t.eq(&com.design, &"~ pair + cond".to_string(), "design with a blocking factor");
        t.eq(&com.coefficients.len(), &7, "intercept, 5 pair indicators, the condition");
        t.check(acertos(&com) as f64 > 1.3 * acertos(&sem) as f64, "the blocking factor gains power");
        t.check((falsos(&com) as f64) < 0.15 * (acertos(&com) + falsos(&com)) as f64, "paired design: false discoveries below 15%");
        let de: Vec<usize> = (0..n_genes).filter(|&g| ef[g] != 0.0 && com.genes[g].base_mean > 20.0).collect();
        let c = corr(&de.iter().map(|&g| com.genes[g].log2_fold_change).collect::<Vec<_>>(), &de.iter().map(|&g| ef[g]).collect::<Vec<_>>());
        t.check(c > 0.9, "paired design: log2 fold changes recovered");
    }

    // three levels: B and C against A
    let nivel = [0, 0, 0, 1, 1, 1, 2, 2, 2];
    let efb: Vec<f64> = (0..n_genes).map(|g| if g % 10 == 1 { 2.0 } else { 0.0 }).collect();
    let efc: Vec<f64> = (0..n_genes).map(|g| if g % 10 == 2 { -2.0 } else { 0.0 }).collect();
    let s = simula(&mut rng, n_genes, &nivel, &[vec![0.0; n_genes], efb.clone(), efc.clone()], None, 0.05, 1.0);
    let cond = niveis_txt(&nivel, &["A", "B", "C"]);
    t.erro(deseq(&s.counts, &s.genes, &amostras(9), &cond, &[], &DeOptions { condition: "g".into(), reference: "A".into(), ..DeOptions::default() }), &["has 3 levels", "A, B, C"], "three levels and no treated level");
    let o3 = DeOptions { condition: "g".into(), reference: "A".into(), treated: Some("C".into()), ..DeOptions::default() };
    if let Some(r) = t.ok_ou(deseq(&s.counts, &s.genes, &amostras(9), &cond, &[], &o3), "three levels, C vs A") {
        let sig_c = (0..n_genes).filter(|&g| r.genes[g].padj < 0.1 && efc[g] != 0.0 && r.genes[g].log2_fold_change < 0.0).count();
        let sig_b = (0..n_genes).filter(|&g| r.genes[g].padj < 0.1 && efb[g] != 0.0).count();
        t.check(sig_c > 200 && sig_b < 15, "C vs A finds the C genes, not the B genes");
        t.eq(&r.coefficients, &vec!["Intercept".to_string(), "g_B_vs_A".into(), "g_C_vs_A".into()], "coefficients of a three-level factor");
    }

    // Cook's distance
    let mut c: Vec<f64> = Vec::new();
    let mut nomes = Vec::new();
    for g in 0..300 {
        let base = 50.0 + (g % 17) as f64 * 20.0;
        for j in 0..8 {
            c.push((base * (1.0 + 0.05 * ((g * 7 + j * 3) % 5) as f64)).round());
        }
        nomes.push(format!("h{g}"));
    }
    c.extend([100.0, 105.0, 98.0, 2000.0, 100.0, 95.0, 110.0, 102.0]);
    nomes.push("alto".into());
    c.extend([1000.0, 1100.0, 1050.0, 3.0, 1000.0, 990.0, 1020.0, 1060.0]);
    nomes.push("baixo".into());
    let cond = niveis_txt(&[0, 0, 0, 0, 1, 1, 1, 1], &["ctl", "trt"]);
    if let Some(r) = t.ok_ou(deseq(&c, &nomes, &amostras(8), &cond, &[], &opt), "Cook's distances") {
        let alto = &r.genes[300];
        let baixo = &r.genes[301];
        let corte = f_quantile(0.99, 2.0, 6.0);
        t.check(alto.outlier && alto.pvalue.is_nan() && alto.max_cooks > corte, "a high outlier count is flagged and its p-value removed");
        t.check(!baixo.outlier && !baixo.pvalue.is_nan(), "an outlier below three other counts is not flagged (DESeq2's rule)");
        t.check(r.outliers >= 1, "outliers counted in the summary");
        let o = DeOptions { cooks_filter: false, ..opt.clone() };
        if let Ok(r2) = deseq(&c, &nomes, &amostras(8), &cond, &[], &o) {
            t.check(!r2.genes[300].outlier && !r2.genes[300].pvalue.is_nan(), "cooks_filter off keeps the p-value");
        }
    }

    // a GLM worked by hand: size factors exactly 1, two groups
    let mut c: Vec<f64> = Vec::new();
    let mut nomes = Vec::new();
    for g in 0..200 {
        for _ in 0..6 {
            c.push(30.0 + g as f64);
        }
        nomes.push(format!("h{g}"));
    }
    for g in 0..20 {
        let v = 40.0 + 7.0 * g as f64;
        c.extend([v, v * 1.5, v * 0.7, v * 1.1, v * 0.9, v * 1.3].map(|x: f64| x.round()));
        nomes.push(format!("v{g}"));
    }
    c.extend([10.0, 20.0, 30.0, 40.0, 50.0, 60.0]);
    nomes.push("G".into());
    let cond = niveis_txt(&[0, 0, 0, 1, 1, 1], &["ctl", "trt"]);
    if let Some(r) = t.ok_ou(deseq(&c, &nomes, &amostras(6), &cond, &[], &opt), "GLM by hand") {
        t.check(r.size_factors.iter().all(|s| (s - 1.0).abs() < 1e-12), "size factors are exactly 1");
        let g = &r.genes[220];
        t.perto(g.log2_fold_change, (50.0f64 / 20.0).log2(), 1e-6, "MLE log2 fold change = log2 of the ratio of group means");
        t.perto(g.base_mean, 35.0, 1e-12, "baseMean = mean of normalized counts");
        let a = g.dispersion;
        let w = |m: f64| m / (1.0 + a * m);
        let se = std::f64::consts::LOG2_E * (1.0 / (3.0 * w(20.0)) + 1.0 / (3.0 * w(50.0))).sqrt();
        t.perto(g.lfc_se, se, 1e-5, "lfcSE from the Fisher information of the NB GLM");
        t.perto(g.stat, g.log2_fold_change / g.lfc_se, 1e-14, "Wald statistic = LFC / SE");
        t.perto(g.pvalue, two_sided_p(g.stat), 1e-14, "p-value = two-sided normal");
    }

    // genes with zeros only, poscounts
    let mut c = s.counts[..9 * 50].to_vec();
    c.extend([0.0; 9]);
    let mut g2: Vec<String> = s.genes[..50].to_vec();
    g2.push("vazio".into());
    let cond = niveis_txt(&[0, 0, 0, 1, 1, 1, 2, 2, 2], &["A", "B", "C"]);
    if let Some(r) = t.ok_ou(deseq(&c, &g2, &amostras(9), &cond, &[], &o3), "a gene with zeros only") {
        let z = &r.genes[50];
        t.check(z.all_zero && z.pvalue.is_nan() && z.padj.is_nan() && z.log2_fold_change.is_nan() && z.base_mean == 0.0, "all-zero gene: baseMean 0, everything else NA");
        t.eq(&r.all_zero, &1, "all-zero genes counted");
    }
    let mut c = s.counts[..9 * 300].to_vec();
    for g in 0..300 {
        c[g * 9 + g % 9] = 0.0;
    }
    if let Some(r) = t.ok_ou(deseq(&c, &s.genes[..300].to_vec(), &amostras(9), &cond, &[], &o3), "every gene has a zero") {
        t.eq(&r.size_factor_method, &"poscounts".to_string(), "every gene with a zero: poscounts size factors");
    }
}

// =====================================================================
//  FILES AND ERRORS
// =====================================================================

fn testa_arquivos(t: &mut T) {
    let mut rng = Rng(5150);
    let n_genes = 800;
    let nivel = [0, 1, 0, 1, 0, 1, 0, 1];
    let ef: Vec<f64> = (0..n_genes).map(|g| if g % 8 == 0 { 2.0 } else { 0.0 }).collect();
    let mut s = simula(&mut rng, n_genes, &nivel, &[vec![0.0; n_genes], ef], None, 0.05, 1.5);
    // a gene with no read: NA everywhere in the results
    s.counts.extend([0.0; 8]);
    s.genes.push("zeros".into());
    let n_genes = n_genes + 1;
    let nomes_a = ["ctl_1", "trt_1", "ctl_2", "trt_2", "ctl_3", "trt_3", "ctl_4", "trt_4"];
    let lote = ["b1", "b1", "b1", "b1", "b2", "b2", "b2", "b2"];
    // featureCounts layout, with an extra sample column, htseq rows and a comment
    let mut l = vec![
        "# Program:featureCounts v2.0.6; Command:\"featureCounts\" \"-a\" \"x.gtf\"".to_string(),
        format!("Geneid\tChr\tStart\tEnd\tStrand\tLength\t{}\textra", nomes_a.join("\t")),
    ];
    for g in 0..n_genes {
        let v: Vec<String> = (0..8).map(|j| (s.counts[g * 8 + j] as u64).to_string()).collect();
        l.push(format!("{}\t1\t{}\t{}\t+\t1000\t{}\t7", s.genes[g], 1000 * g + 1, 1000 * g + 500, v.join("\t")));
    }
    l.push("__no_feature\t.\t.\t.\t.\t.\t5\t5\t5\t5\t5\t5\t5\t5\t5".into());
    let mtz = t.p("counts.tsv.gz");
    grava(&mtz, &l);
    let tab = t.p("samples.tsv");
    let mut ls = vec!["sample\tcondition\tbatch".to_string()];
    for j in 0..8 {
        ls.push(format!("{}\t{}\t{}", nomes_a[j], if nivel[j] == 0 { "ctl" } else { "trt" }, lote[j]));
    }
    grava(&tab, &ls);
    let csv = t.p("samples.csv");
    let lc: Vec<String> = ls.iter().map(|x| x.replace('\t', ",")).collect();
    grava(&csv, &lc);
    let nomes_genes = t.p("names.tsv");
    let mut ln = vec!["gene_id\tname".to_string()];
    for g in 0..n_genes - 6 {
        ln.push(format!("{}\tSYM{g}", s.genes[g]));
    }
    grava(&nomes_genes, &ln);
    let opt = DeOptions {
        condition: "condition".into(),
        reference: "ctl".into(),
        factors: vec!["batch".into()],
        gene_names: Some(nomes_genes.clone()),
        ..DeOptions::default()
    };
    let saida = t.p("de.tsv");
    let Some(r) = t.ok_ou(differential_expression(&mtz, &tab, Some(&saida), &opt), "differential_expression from files") else { return };
    t.eq(&r.ignored_columns, &vec!["Chr".to_string(), "Start".into(), "End".into(), "Strand".into(), "Length".into(), "extra".into()], "columns not in the sample table are ignored");
    t.eq(&r.skipped_rows, &1, "htseq-count's __ rows skipped");
    t.eq(&r.names_missing, &6, "genes missing from the name table");
    t.eq(&r.all_zero, &1, "the gene without reads is counted");
    t.eq(&r.design, &"~ batch + condition".to_string(), "design from the files");
    let cond = niveis_txt(&nivel, &["ctl", "trt"]);
    let mem = deseq(&s.counts, &s.genes, &nomes_a.iter().map(|x| x.to_string()).collect::<Vec<_>>(), &cond, &[("batch".into(), lote.iter().map(|x| x.to_string()).collect())], &opt).unwrap();
    t.check(mem.genes.iter().zip(&r.genes).all(|(a, b)| a.pvalue.to_bits() == b.pvalue.to_bits() && a.padj.to_bits() == b.padj.to_bits()), "files give the same result as the matrix in memory");
    let linhas = le(&saida);
    t.eq(&linhas[0], &"gene\tname\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj".to_string(), "results header");
    t.eq(&linhas.len(), &(n_genes + 1), "one row per gene");
    let ordem = r.order();
    let mut ok = true;
    let mut ultimo = -1.0f64;
    let mut na = false;
    for (k, &i) in ordem.iter().enumerate() {
        let g = &r.genes[i as usize];
        let c: Vec<&str> = linhas[k + 1].split('\t').collect();
        ok &= c[0] == g.gene && c[1] == g.name && c[7] == num(g.padj) && c[6] == num(g.pvalue) && c[3] == num(g.log2_fold_change);
        if g.padj.is_nan() {
            na = true;
        } else {
            ok &= !na && g.padj >= ultimo;
            ultimo = g.padj;
        }
        ok &= c[2].parse::<f64>().map_or(false, |v| v == g.base_mean);
    }
    t.check(ok, "results table sorted by padj, NA last, numbers exact");
    t.check(linhas.iter().any(|l| l.ends_with("\tNA")), "NA written for missing padj");
    t.check(linhas[1..].iter().filter(|l| l.split('\t').nth(1).map_or(false, |n| n.starts_with('g'))).count() == 5, "a gene without name keeps its id as name");
    if let Some(r2) = t.ok_ou(differential_expression(&mtz, &csv, None, &DeOptions { gene_names: None, ..opt.clone() }), "comma-separated sample table") {
        t.check(r2.genes.iter().zip(&r.genes).all(|(a, b)| a.padj.to_bits() == b.padj.to_bits()), "CSV sample table: same result");
    }
    // a plain matrix written by count_reads' layout, CSV
    let mtz2 = t.p("counts.csv");
    let mut l2 = vec![format!("gene,{}", nomes_a.join(","))];
    for g in 0..n_genes {
        l2.push(format!("{},{}", s.genes[g], (0..8).map(|j| (s.counts[g * 8 + j] as u64).to_string()).collect::<Vec<_>>().join(",")));
    }
    grava(&mtz2, &l2);
    if let Some(r2) = t.ok_ou(differential_expression(&mtz2, &tab, None, &DeOptions { gene_names: None, ..opt.clone() }), "CSV matrix") {
        t.check(r2.genes.iter().zip(&r.genes).all(|(a, b)| a.padj.to_bits() == b.padj.to_bits()), "CSV matrix: same result");
    }
    // NCBI gene_info
    let gi = t.p("gene_info.tsv");
    grava(&gi, &["#tax_id\tGeneID\tSymbol\tLocusTag".to_string(), format!("9606\t{}\tFIRST\t-", s.genes[0])]);
    if let Some(r2) = t.ok_ou(differential_expression(&mtz2, &tab, None, &DeOptions { gene_names: Some(gi), ..opt.clone() }), "gene_info names") {
        t.eq(&r2.genes[0].name, &"FIRST".to_string(), "names from NCBI gene_info (GeneID -> Symbol)");
        t.eq(&r2.names_missing, &(n_genes - 1), "gene_info: the others keep their id");
    }
    // only the compared levels
    let tab3 = t.p("samples3.tsv");
    let mut l3 = ls.clone();
    l3.push("other_1\tdrug\tb1".into());
    grava(&tab3, &l3);
    let mtz3 = t.p("counts3.tsv");
    let mut l4 = vec![format!("gene\t{}\tother_1", nomes_a.join("\t"))];
    for g in 0..n_genes {
        l4.push(format!("{}\t{}\t{}", s.genes[g], (0..8).map(|j| (s.counts[g * 8 + j] as u64).to_string()).collect::<Vec<_>>().join("\t"), 10 + g % 7));
    }
    grava(&mtz3, &l4);
    let o = DeOptions { treated: Some("trt".into()), only_compared_levels: true, gene_names: None, ..opt.clone() };
    if let Some(r2) = t.ok_ou(differential_expression(&mtz3, &tab3, None, &o), "only the compared levels") {
        t.eq(&r2.samples.len(), &8, "only_compared_levels leaves the third level out");
        t.check(r2.genes.iter().zip(&r.genes).all(|(a, b)| a.padj.to_bits() == b.padj.to_bits()), "only_compared_levels: same as without the extra sample");
    }

    // errors
    let o = DeOptions { gene_names: None, ..opt.clone() };
    t.erro(differential_expression(&mtz, &tab, None, &DeOptions { condition: "treatment".into(), ..o.clone() }), &["samples.tsv has no column \"treatment\"", "sample, condition, batch"], "condition column missing");
    t.erro(differential_expression(&mtz, &tab, None, &DeOptions { reference: "placebo".into(), ..o.clone() }), &["\"placebo\" is not a level", "ctl, trt"], "reference level absent");
    t.erro(differential_expression(&mtz, &tab, None, &DeOptions { reference: String::new(), ..o.clone() }), &["reference level"], "no reference level");
    t.erro(differential_expression(&mtz, &tab, Some(&mtz), &o), &["itself"], "output over the counts");
    t.erro(differential_expression(&mtz3, &tab3, None, &DeOptions { only_compared_levels: true, ..o.clone() }), &["only_compared_levels needs the treated level"], "only_compared_levels without treated");
    let ruim = t.p("ruim_samples.tsv");
    grava(&ruim, &["sample\tcondition".to_string(), "ctl_1\tctl".into(), "nobody\ttrt".into()]);
    t.erro(differential_expression(&mtz, &ruim, None, &DeOptions { factors: vec![], ..o.clone() }), &["Sample \"nobody\"", "is not a column of", "ctl_1"], "sample absent from the matrix");
    grava(&ruim, &["sample\tcondition".to_string(), "ctl_1\tctl".into(), "ctl_1\ttrt".into()]);
    t.erro(differential_expression(&mtz, &ruim, None, &DeOptions { factors: vec![], ..o.clone() }), &["appears twice", "lines 2 and 3"], "duplicate sample");
    grava(&ruim, &["sample\tcondition\tbatch".to_string(), "ctl_1\tctl\tb1".into(), "trt_1\ttrt\tb2".into(), "ctl_2\tctl\tb1".into(), "trt_2\ttrt\tb2".into()]);
    t.erro(differential_expression(&mtz, &ruim, None, &o), &["not full rank", "~ batch + condition"], "blocking factor confounded with the condition");
    grava(&ruim, &["sample\tcondition".to_string(), "ctl_1\tctl".into(), "trt_1\ttrt".into()]);
    t.erro(differential_expression(&mtz, &ruim, None, &DeOptions { factors: vec![], ..o.clone() }), &["no residual degrees of freedom", "replicates"], "no replicates");
    let m_ruim = t.p("ruim_counts.tsv");
    grava(&m_ruim, &[format!("gene\t{}", nomes_a.join("\t")), "a\t1\t2\t3\t4\t5\t6\t7\t8".into(), "b\t1\t2\t3.5\t4\t5\t6\t7\t8".into()]);
    t.erro(differential_expression(&m_ruim, &tab, None, &o), &["ruim_counts.tsv, line 3, column 4 (sample ctl_2)", "3.5", "whole number", "rounded"], "non-integer count");
    grava(&m_ruim, &[format!("gene\t{}", nomes_a.join("\t")), "a\t1\t2\t3\t4\t5\t6\t-7\t8".into()]);
    t.erro(differential_expression(&m_ruim, &tab, None, &o), &["line 2", "-7"], "negative count");
    grava(&m_ruim, &[format!("gene\t{}", nomes_a.join("\t")), "a\t1\t2\t3\t4\t5\t6\tNA\t8".into()]);
    t.erro(differential_expression(&m_ruim, &tab, None, &o), &["line 2", "\"NA\" is not a number"], "NA count");
    grava(&m_ruim, &[format!("gene\t{}", nomes_a.join("\t")), "a\t1\t2\t3\t4\t5\t6\t7\t8".into(), "a\t1\t2\t3\t4\t5\t6\t7\t8".into()]);
    t.erro(differential_expression(&m_ruim, &tab, None, &o), &["gene \"a\" appears twice", "lines 2 and 3"], "duplicate gene");
    grava(&m_ruim, &[format!("gene\t{}", nomes_a.join("\t")), "a\t1\t2\t3\t4\t5".into()]);
    t.erro(differential_expression(&m_ruim, &tab, None, &o), &["line 2", "column(s)"], "short row");
    grava(&m_ruim, &[format!("gene\t{}", nomes_a.join("\t")), "a\t0\t0\t0\t0\t0\t0\t0\t0".into()]);
    t.erro(differential_expression(&m_ruim, &tab, None, &o), &["every count is zero"], "all counts zero");
    grava(&m_ruim, &["# nothing".to_string()]);
    t.erro(differential_expression(&m_ruim, &tab, None, &o), &["is empty"], "empty matrix");
    t.erro(differential_expression(&t.p("nada.tsv"), &tab, None, &o), &["Error opening", "nada.tsv"], "missing matrix");
    t.erro(differential_expression(&mtz, &tab, None, &DeOptions { alpha: 1.5, ..o.clone() }), &["alpha"], "alpha out of range");
    t.erro(differential_expression(&mtz, &tab, None, &DeOptions { prior_variance: Some(-1.0), ..o.clone() }), &["prior_variance"], "negative prior variance");
    t.erro(differential_expression(&mtz, &tab, None, &DeOptions { factors: vec!["condition".into()], ..o.clone() }), &["is the condition itself"], "condition listed as a factor");
}

// =====================================================================
//  REAL DATA
// =====================================================================

fn real(args: &[String]) {
    if args.len() < 4 {
        eprintln!("usage: expressao airway <counts> <samples> <gene_info> <out.tsv> [prior_variance]");
        std::process::exit(2);
    }
    let inicio = std::time::Instant::now();
    let opt = DeOptions {
        condition: "dex".into(),
        reference: "untrt".into(),
        factors: vec!["cell_line".into()],
        gene_names: Some(args[2].clone()),
        prior_variance: args.get(4).and_then(|v| v.parse().ok()),
        ..DeOptions::default()
    };
    match differential_expression(&args[0], &args[1], Some(&args[3]), &opt) {
        Ok(r) => {
            print!("{}", r.resumo());
            println!("  log dispersion variance {:.4}; {:.1} s", r.log_dispersion_variance, inicio.elapsed().as_secs_f64());
            println!("  top genes:");
            for &i in r.order().iter().take(15) {
                let g = &r.genes[i as usize];
                println!("    {:<10} baseMean {:>9.1}  log2FC {:>6.2}  padj {:.2e}", g.name, g.base_mean, g.log2_fold_change, g.padj);
            }
            println!("  known dexamethasone genes:");
            for nome in ["FKBP5", "TSC22D3", "DUSP1", "KLF15", "CRISPLD2", "ZBTB16", "PER1"] {
                if let Some(g) = r.genes.iter().find(|g| g.name == nome) {
                    println!("    {:<10} baseMean {:>9.1}  log2FC {:>6.2}  padj {:.2e}", g.name, g.base_mean, g.log2_fold_change, g.padj);
                }
            }
        }
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(1);
        }
    }
}

/// `de <counts> <samples> <condition> <reference> <factors|-> <out_prefix>
/// [gene_names|-] [prior_variance|-] [treated|-]`: writes <out>.tsv (the
/// results table), <out>_full.tsv (with dispersions and flags) and
/// <out>.json.
fn modo_de(args: &[String]) {
    if args.len() < 6 {
        eprintln!("usage: expressao de <counts> <samples> <condition> <reference> <factors|-> <out_prefix> [gene_names|-] [prior_variance|-] [treated|-]");
        std::process::exit(2);
    }
    let opcional = |i: usize| args.get(i).filter(|s| s.as_str() != "-").cloned();
    let opt = DeOptions {
        condition: args[2].clone(),
        reference: args[3].clone(),
        factors: if args[4] == "-" { vec![] } else { args[4].split(',').map(str::to_string).collect() },
        gene_names: opcional(6),
        prior_variance: opcional(7).and_then(|v| v.parse().ok()),
        treated: opcional(8),
        ..DeOptions::default()
    };
    let inicio = std::time::Instant::now();
    let saida = format!("{}.tsv", args[5]);
    let r = match differential_expression(&args[0], &args[1], Some(&saida), &opt) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(1);
        }
    };
    let s = inicio.elapsed().as_secs_f64();
    let mut l = vec!["gene\tname\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\tdisp_gene\tdisp_trend\tdispersion\tmax_cooks\tall_zero\toutlier\tfiltered\tnot_converged".to_string()];
    for g in &r.genes {
        l.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            g.gene,
            g.name,
            num(g.base_mean),
            num(g.log2_fold_change),
            num(g.lfc_se),
            num(g.stat),
            num(g.pvalue),
            num(g.padj),
            num(g.disp_gene),
            num(g.disp_trend),
            num(g.dispersion),
            num(g.max_cooks),
            g.all_zero,
            g.outlier,
            g.filtered,
            g.not_converged
        ));
    }
    grava(&format!("{}_full.tsv", args[5]), &l);
    std::fs::write(format!("{}.json", args[5]), r.json()).unwrap();
    print!("{}", r.resumo());
    println!("  log dispersion variance {:.4}; {:.2} s", r.log_dispersion_variance, s);
}

/// `dump <dir>`: the simulated data sets of the tests, as files, for
/// comparison with other implementations.
fn modo_dump(dir: &str) {
    std::fs::create_dir_all(dir).unwrap();
    let escreve = |nome: &str, s: &Sim, nivel: &[usize], extra: Option<(&str, Vec<String>)>| {
        let n = nivel.len();
        let mut l = vec![format!("gene\t{}", amostras(n).join("\t"))];
        for (g, gene) in s.genes.iter().enumerate() {
            l.push(format!("{gene}\t{}", (0..n).map(|j| (s.counts[g * n + j] as u64).to_string()).collect::<Vec<_>>().join("\t")));
        }
        grava(&format!("{dir}/{nome}_counts.tsv"), &l);
        let mut a = vec![format!("sample\tcond{}", extra.as_ref().map(|e| format!("\t{}", e.0)).unwrap_or_default())];
        for j in 0..n {
            a.push(format!("s{j}\t{}{}", if nivel[j] == 0 { "ctl" } else { "trt" }, extra.as_ref().map(|e| format!("\t{}", e.1[j])).unwrap_or_default()));
        }
        grava(&format!("{dir}/{nome}_samples.tsv"), &a);
        let lf: Vec<String> = std::iter::once("gene\ttrue_lfc".to_string()).chain(s.genes.iter().zip(&s.lfc).map(|(g, v)| format!("{g}\t{v}"))).collect();
        grava(&format!("{dir}/{nome}_truth.tsv"), &lf);
    };
    // the same draws as testa_simulacao
    let mut rng = Rng(424242);
    let n_genes = 3000;
    let nivel = [0, 0, 0, 0, 1, 1, 1, 1];
    let ef: Vec<f64> = (0..n_genes).map(|g| if g % 10 == 0 { (1.0 + 2.0 * rng.unif()) * if rng.chance(0.5) { 1.0 } else { -1.0 } } else { 0.0 }).collect();
    let s = simula(&mut rng, n_genes, &nivel, &[vec![0.0; n_genes], ef], None, 0.05, 2.0);
    escreve("sim4v4", &s, &nivel, None);
    let s0 = simula(&mut rng, n_genes, &nivel, &[vec![0.0; n_genes], vec![0.0; n_genes]], None, 0.05, 2.0);
    escreve("simnull", &s0, &nivel, None);
    let nivel = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1];
    let ef: Vec<f64> = (0..n_genes).map(|g| if g % 10 == 3 { (0.6 + 1.0 * rng.unif()) * if rng.chance(0.5) { 1.0 } else { -1.0 } } else { 0.0 }).collect();
    let par_ef: Vec<Vec<f64>> = (0..6).map(|_| (0..n_genes).map(|_| 0.8 * rng.normal()).collect()).collect();
    let bloco: Vec<Vec<f64>> = (0..12).map(|j| par_ef[j / 2].clone()).collect();
    let s = simula(&mut rng, n_genes, &nivel, &[vec![0.0; n_genes], ef], Some(&bloco), 0.02, 1.0);
    escreve("simpaired", &s, &nivel, Some(("pair", (0..12).map(|j| format!("p{}", j / 2)).collect())));
    println!("written to {dir}");
}

/// `gff <annotation> <first sequence name> [attribute]`: builds the index of
/// a real annotation (through count_reads on an alignment file with no read)
/// and reports genes, features and time.
fn modo_gff(args: &[String]) {
    let dir = std::env::temp_dir().join(format!("mmgen_expressao_gff_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let sam = dir.join("vazio.sam").to_string_lossy().to_string();
    grava(&sam, &[format!("@SQ\tSN:{}\tLN:1000", args[1])]);
    let opt = CountOptions { attribute: args.get(2).cloned(), ..CountOptions::default() };
    let inicio = std::time::Instant::now();
    match count_reads(&[sam], &args[0], None, &opt) {
        Ok(r) => {
            println!(
                "{} genes, {} features on {} sequences, gene = {}; {:.1} s",
                r.genes.len(),
                r.features,
                r.sequences,
                r.attribute,
                inicio.elapsed().as_secs_f64()
            );
            println!("first genes: {:?}", &r.genes[..r.genes.len().min(8)]);
            for w in &r.warnings {
                println!("note: {w}");
            }
            let ids: Vec<String> = r.genes.clone();
            std::fs::write(dir.parent().unwrap().join("mmgen_expressao_gff_genes.txt"), ids.join("\n")).unwrap();
        }
        Err(e) => println!("error: {e}"),
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// `bench <fragments> <dir>`: a large synthetic BAM, counted twice (BAM and
/// SAM), with the time.
fn modo_bench(n: usize, dir: &str) {
    std::fs::create_dir_all(dir).unwrap();
    let mut rng = Rng(77);
    let a = gera_anotacao(&mut rng);
    let gtf = format!("{dir}/bench.gtf");
    grava(&gtf, &escreve_gtf(&a));
    let mut recs = Vec::new();
    let lote = 200_000;
    let mut feitos = 0;
    while feitos < n {
        let k = lote.min(n - feitos);
        let mut r = gera_leituras(&mut rng, &a, k);
        for x in r.iter_mut() {
            x.nome = format!("b{feitos}_{}", x.nome);
        }
        recs.extend(r);
        feitos += k;
    }
    ordena_recs(&mut recs);
    let mut l = cabecalho_sam(true);
    l.extend(recs.iter().map(Rec::sam));
    let bam = format!("{dir}/bench.bam");
    let sam = format!("{dir}/bench.sam");
    grava(&bam, &l);
    grava(&sam, &l);
    let e = forca_bruta(&recs, &a, Strandedness::Unstranded, 10);
    drop(recs);
    drop(l);
    for p in [&bam, &sam] {
        let inicio = std::time::Instant::now();
        let r = count_reads(&[p.clone()], &gtf, None, &CountOptions::default()).unwrap();
        let s = inicio.elapsed().as_secs_f64();
        let mb = std::fs::metadata(p).unwrap().len() as f64 / 1e6;
        println!(
            "{p}: {:.0} MB, {} records in {:.2} s ({:.2} M records/s); equals brute force: {}",
            mb,
            r.samples[0].records,
            s,
            r.samples[0].records as f64 / s / 1e6,
            do_resultado(&r.samples[0], &r.counts[0]) == e
        );
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("airway") => return real(&args[2..]),
        Some("de") => return modo_de(&args[2..]),
        Some("dump") if args.len() > 2 => return modo_dump(&args[2]),
        Some("gff") if args.len() > 3 => return modo_gff(&args[2..]),
        Some("bench") if args.len() > 3 => return modo_bench(args[2].parse().unwrap(), &args[3]),
        _ => {}
    }
    let dir = std::env::temp_dir().join(format!("mmgen_expressao_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut t = T { ok: 0, falhas: 0, dir: dir.clone() };
    let inicio = std::time::Instant::now();

    testa_funcoes(&mut t);
    testa_contagem_a_mao(&mut t);
    testa_contagem_sintetica(&mut t);
    testa_passos(&mut t);
    testa_simulacao(&mut t);
    testa_arquivos(&mut t);

    let _ = std::fs::remove_dir_all(&dir);
    println!("{} of {} checks passed ({:.1} s)", t.ok, t.ok + t.falhas, inicio.elapsed().as_secs_f64());
    if t.falhas > 0 {
        std::process::exit(1);
    }
}
