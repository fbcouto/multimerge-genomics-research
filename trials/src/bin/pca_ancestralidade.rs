// Trials of src/pca_ancestralidade.rs: the PCA of people against the groups
// of a frequency panel.
//
// Every result is compared with something computed here on its own:
//  - the eigensolver against a cyclic Jacobi written here, and against the
//    definition (A v = lambda v, V'V = I, trace);
//  - the PCA of a reference against the Gram matrix built and diagonalised
//    here; the least-squares projection against an exact answer when sites
//    are missing (and against the shrinkage of plain mean imputation);
//  - the whole tool on populations simulated with a known F_ST
//    (Balding-Nichols): each person must fall in the cloud of the group it
//    was drawn from, admixed people at the mixture of the centres, a person
//    with most sites missing not pulled toward 0, the distance between the
//    centres equal to the one computed from the panel's frequencies, the
//    callable-BED reading exact, the projection equal to the one computed
//    here from the sites table, the same JSON for the same seed.
//
// Usage: cargo run --release --bin pca_ancestralidade
//        cargo run --release --bin pca_ancestralidade -- real
//          (the 1000 Genomes and GIAB run; see `real` at the end for the
//          environment variables)
#[path = "../../../src/multimerge_core.rs"]
#[allow(dead_code)]
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
#[allow(dead_code)]
mod ancestralidade;
#[path = "../../../src/pca_ancestralidade.rs"]
mod pca_ancestralidade;

use pca_ancestralidade::*;
use serde_json::Value;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Read};
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
        self.check(what, (got - want).abs() <= tol, || format!("got {got}, expected {want} (tolerance {tol})"));
    }

    fn eq<X: PartialEq + std::fmt::Debug>(&mut self, what: &str, got: X, want: X) {
        let ok = got == want;
        self.check(what, ok, || format!("got {got:?}, expected {want:?}"));
    }

    fn fails<X: std::fmt::Debug>(&mut self, what: &str, r: Result<X, String>, needles: &[&str]) {
        match r {
            Ok(x) => {
                let s = format!("{x:?}");
                self.check(what, false, || format!("should have failed, returned {}", &s[..s.len().min(200)]))
            }
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
    /// Gamma(a, 1), Marsaglia and Tsang.
    fn gamma(&mut self, a: f64) -> f64 {
        if a < 1.0 {
            let u = self.f().max(1e-300);
            return self.gamma(a + 1.0) * u.powf(1.0 / a);
        }
        let d = a - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let x = self.normal();
            let v = (1.0 + c * x).powi(3);
            if v <= 0.0 {
                continue;
            }
            let u = self.f().max(1e-300);
            if u.ln() < 0.5 * x * x + d - d * v + d * v.ln() {
                return d * v;
            }
        }
    }
    fn beta(&mut self, a: f64, b: f64) -> f64 {
        let x = self.gamma(a);
        let y = self.gamma(b);
        x / (x + y)
    }
    fn par_bases(&mut self) -> (char, char) {
        let b = ['A', 'C', 'G', 'T'];
        let i = self.below(4);
        let j = (i + 1 + self.below(3)) % 4;
        (b[i], b[j])
    }
}

// ---------------------------------------------------------------------
//  SIMULATED PANELS AND PEOPLE
// ---------------------------------------------------------------------

struct Painel {
    grupos: Vec<String>,
    crom: Vec<String>,
    pos: Vec<u64>,
    ra: Vec<(char, char)>,
    /// Per site, per group, as written (6 decimals).
    p: Vec<Vec<f64>>,
}

/// Balding-Nichols: an ancestral frequency p0 ~ U(0.05, 0.95) per site, and
/// each group's frequency ~ Beta(p0 (1-F)/F, (1-p0) (1-F)/F), so that F is
/// the F_ST of each group against the ancestral population.
fn painel_bn(r: &mut Rng, grupos: &[&str], cromos: &[&str], por_crom: usize, passo: usize, f: f64) -> Painel {
    let mut p = Painel { grupos: grupos.iter().map(|s| s.to_string()).collect(), crom: Vec::new(), pos: Vec::new(), ra: Vec::new(), p: Vec::new() };
    for c in cromos {
        let mut pos = 10_000u64;
        for _ in 0..por_crom {
            pos += 1 + r.below(passo) as u64;
            let p0 = 0.05 + 0.9 * r.f();
            let fs: Vec<f64> = grupos
                .iter()
                .map(|_| {
                    let x = r.beta(p0 * (1.0 - f) / f, (1.0 - p0) * (1.0 - f) / f);
                    (x.clamp(0.0, 1.0) * 1e6).round() / 1e6
                })
                .collect();
            p.crom.push(c.to_string());
            p.pos.push(pos);
            p.ra.push(r.par_bases());
            p.p.push(fs);
        }
    }
    p
}

fn escreve_painel(t: &T, nome: &str, p: &Painel) -> String {
    let mut l = vec!["##fileformat=MMGEN-panel-1".to_string()];
    let mut c = "chrom\tpos\tref\talt".to_string();
    for g in &p.grupos {
        c.push_str(&format!("\tAF_{g}"));
    }
    for g in &p.grupos {
        c.push_str(&format!("\tAN_{g}"));
    }
    l.push(c);
    for j in 0..p.pos.len() {
        let mut s = format!("{}\t{}\t{}\t{}", p.crom[j], p.pos[j], p.ra[j].0, p.ra[j].1);
        for x in &p.p[j] {
            s.push_str(&format!("\t{x}"));
        }
        for _ in &p.grupos {
            s.push_str("\t2000");
        }
        l.push(s);
    }
    t.write(nome, &l)
}

/// ALT count per site of a person drawn under the admixture model (each
/// allele picks a group by q); `None` = missing with probability `falta`.
fn pessoa(r: &mut Rng, p: &Painel, q: &[f64], falta: f64) -> Vec<Option<u8>> {
    (0..p.pos.len())
        .map(|j| {
            let mut x = 0u8;
            for _ in 0..2 {
                let mut u = r.f();
                let mut g = 0;
                while g + 1 < q.len() && u >= q[g] {
                    u -= q[g];
                    g += 1;
                }
                if r.chance(p.p[j][g]) {
                    x += 1;
                }
            }
            if falta > 0.0 && r.chance(falta) {
                None
            } else {
                Some(x)
            }
        })
        .collect()
}

/// A VCF of several people. `so_variantes`: only records where someone is
/// not 0/0 (the GIAB style; absent = 0/0 inside the callable BED).
fn escreve_vcf(t: &T, nome: &str, prefixo: &str, p: &Painel, nomes: &[String], gts: &[Vec<Option<u8>>], so_variantes: bool) -> String {
    let mut l = vec!["##fileformat=VCFv4.2".to_string(), "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string()];
    l.push(format!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", nomes.join("\t")));
    for j in 0..p.pos.len() {
        if so_variantes && gts.iter().all(|g| g[j] == Some(0)) {
            continue;
        }
        let campos: Vec<&str> = gts
            .iter()
            .map(|g| match g[j] {
                None => "./.",
                Some(0) => "0/0",
                Some(1) => {
                    if j % 2 == 0 {
                        "0|1"
                    } else {
                        "1/0"
                    }
                }
                Some(_) => "1/1",
            })
            .collect();
        l.push(format!("{prefixo}{}\t{}\t.\t{}\t{}\t50\tPASS\t.\tGT\t{}", p.crom[j], p.pos[j], p.ra[j].0, p.ra[j].1, campos.join("\t")));
    }
    t.write(nome, &l)
}

fn pcs(v: &Value) -> Vec<f64> {
    v.as_array().map(|a| a.iter().map(|x| x.as_f64().unwrap_or(f64::NAN)).collect()).unwrap_or_default()
}

fn dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum::<f64>().sqrt()
}

/// The sites table: (chrom, pos, AF per group, mean, sd, loadings).
struct TabelaSitios {
    crom: Vec<String>,
    pos: Vec<u64>,
    ra: Vec<(String, String)>,
    af: Vec<Vec<f64>>,
    media: Vec<f64>,
    dp: Vec<f64>,
    v: Vec<Vec<f64>>,
}

fn le_sitios(p: &str, k: usize) -> TabelaSitios {
    let mut t = TabelaSitios { crom: Vec::new(), pos: Vec::new(), ra: Vec::new(), af: Vec::new(), media: Vec::new(), dp: Vec::new(), v: Vec::new() };
    let linhas = le_texto(p);
    let mut ng = 0;
    for l in &linhas {
        if l.starts_with("##") {
            continue;
        }
        let c: Vec<&str> = l.split('\t').collect();
        if c[0] == "chrom" {
            ng = c.iter().filter(|x| x.starts_with("AF_")).count();
            continue;
        }
        t.crom.push(c[0].to_string());
        t.pos.push(c[1].parse().unwrap());
        t.ra.push((c[2].to_string(), c[3].to_string()));
        t.af.push((0..ng).map(|g| c[4 + g].parse().unwrap()).collect());
        t.media.push(c[4 + ng].parse().unwrap());
        t.dp.push(c[5 + ng].parse().unwrap());
        t.v.push((0..k).map(|q| c[6 + ng + q].parse().unwrap()).collect());
    }
    t
}

impl TabelaSitios {
    /// A PcaFit built from the table, to project with the published axes.
    fn ajuste(&self, k: usize) -> PcaFit {
        PcaFit {
            individuals: 0,
            sites: self.pos.len(),
            components: k,
            mean: self.media.clone(),
            sd: self.dp.clone(),
            eigenvalues: Vec::new(),
            total_variance: 0.0,
            loadings: self.v.iter().flat_map(|x| x.iter().copied()).collect(),
            scores: Vec::new(),
        }
    }
}

// ---------------------------------------------------------------------
//  LINEAR ALGEBRA WRITTEN HERE
// ---------------------------------------------------------------------

/// Cyclic Jacobi: the eigenvalues, ascending.
fn jacobi(a: &[f64], n: usize) -> Vec<f64> {
    let mut m = a.to_vec();
    let frob: f64 = m.iter().map(|x| x * x).sum();
    for _ in 0..100 {
        let mut off = 0.0;
        for p in 0..n {
            for q in 0..n {
                if p != q {
                    off += m[p * n + q] * m[p * n + q];
                }
            }
        }
        if off <= 1e-26 * frob.max(1e-300) {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                let apq = m[p * n + q];
                if apq.abs() < 1e-300 {
                    continue;
                }
                let theta = (m[q * n + q] - m[p * n + p]) / (2.0 * apq);
                let sinal = if theta >= 0.0 { 1.0 } else { -1.0 };
                let tt = sinal / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (tt * tt + 1.0).sqrt();
                let s = tt * c;
                for k in 0..n {
                    let (akp, akq) = (m[k * n + p], m[k * n + q]);
                    m[k * n + p] = c * akp - s * akq;
                    m[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (m[p * n + k], m[q * n + k]);
                    m[p * n + k] = c * apk - s * aqk;
                    m[q * n + k] = s * apk + c * aqk;
                }
            }
        }
    }
    let mut d: Vec<f64> = (0..n).map(|i| m[i * n + i]).collect();
    d.sort_by(|a, b| a.partial_cmp(b).unwrap());
    d
}

/// Procrustes correlation of two configurations (rows are points): both
/// centred and scaled to unit size, then the best rotation or reflection;
/// the sum of the singular values of A'B (1 = identical shapes).
fn procrustes(a: &[Vec<f64>], b: &[Vec<f64>]) -> f64 {
    let n = a.len();
    let d = a[0].len();
    let centra = |x: &[Vec<f64>]| -> Vec<Vec<f64>> {
        let m: Vec<f64> = (0..d).map(|c| x.iter().map(|r| r[c]).sum::<f64>() / n as f64).collect();
        let y: Vec<Vec<f64>> = x.iter().map(|r| (0..d).map(|c| r[c] - m[c]).collect()).collect();
        let s: f64 = y.iter().flat_map(|r| r.iter()).map(|v| v * v).sum::<f64>().sqrt();
        y.into_iter().map(|r| r.into_iter().map(|v| v / s).collect()).collect()
    };
    let (a, b) = (centra(a), centra(b));
    let mut mtm = vec![0.0; d * d];
    let mut m = vec![0.0; d * d];
    for i in 0..n {
        for p in 0..d {
            for q in 0..d {
                m[p * d + q] += a[i][p] * b[i][q];
            }
        }
    }
    for p in 0..d {
        for q in 0..d {
            mtm[p * d + q] = (0..d).map(|r| m[r * d + p] * m[r * d + q]).sum();
        }
    }
    let (ev, _) = symmetric_eigen(&mtm, d).unwrap();
    ev.iter().map(|x| x.max(0.0).sqrt()).sum()
}

/// The mixture of the centres closest to a point: least squares with the
/// weights >= 0 summing to 1 (projected gradient on the simplex).
fn mistura(centros: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    let k = centros.len();
    let d = x.len();
    let mut q = vec![1.0 / k as f64; k];
    let mut l = 0.0;
    for a in centros {
        for b in centros {
            l += a.iter().zip(b).map(|(u, v)| u * v).sum::<f64>().abs();
        }
    }
    let passo = 1.0 / l.max(1e-12);
    for _ in 0..20000 {
        let r: Vec<f64> = (0..d).map(|c| (0..k).map(|g| q[g] * centros[g][c]).sum::<f64>() - x[c]).collect();
        let grad: Vec<f64> = (0..k).map(|g| (0..d).map(|c| centros[g][c] * r[c]).sum::<f64>()).collect();
        let y: Vec<f64> = (0..k).map(|g| q[g] - passo * grad[g]).collect();
        // Projection onto the simplex (sort-based).
        let mut u = y.clone();
        u.sort_by(|a, b| b.partial_cmp(a).unwrap());
        let mut acum = 0.0;
        let mut theta = 0.0;
        for (i, ui) in u.iter().enumerate() {
            acum += ui;
            let t = (acum - 1.0) / (i as f64 + 1.0);
            if ui - t > 0.0 {
                theta = t;
            }
        }
        q = y.iter().map(|v| (v - theta).max(0.0)).collect();
    }
    q
}

fn correlacao(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        sab += (x - ma) * (y - mb);
        saa += (x - ma) * (x - ma);
        sbb += (y - mb) * (y - mb);
    }
    sab / (saa * sbb).sqrt()
}

// ---------------------------------------------------------------------
//  1. THE EIGENSOLVER
// ---------------------------------------------------------------------

fn confere_autos(t: &mut T, what: &str, a: &[f64], n: usize) {
    let (d, v) = match symmetric_eigen(a, n) {
        Ok(x) => x,
        Err(e) => {
            t.check(what, false, || e);
            return;
        }
    };
    let norma = a.iter().map(|x| x.abs()).fold(0.0f64, f64::max).max(1e-300);
    let mut res: f64 = 0.0;
    let mut orto: f64 = 0.0;
    for c in 0..n {
        let vc = &v[c * n..c * n + n];
        for i in 0..n {
            let av: f64 = (0..n).map(|j| a[i * n + j] * vc[j]).sum();
            res = res.max((av - d[c] * vc[i]).abs());
        }
        for e in 0..n {
            let ve = &v[e * n..e * n + n];
            let p: f64 = vc.iter().zip(ve).map(|(x, y)| x * y).sum();
            orto = orto.max((p - if c == e { 1.0 } else { 0.0 }).abs());
        }
    }
    t.check(&format!("{what}: A v = lambda v"), res <= 1e-10 * norma * n as f64, || format!("residual {res:e}, |A| {norma:e}"));
    t.check(&format!("{what}: eigenvectors orthonormal"), orto <= 1e-10 * n as f64, || format!("max |V'V - I| {orto:e}"));
    t.check(&format!("{what}: ascending"), d.windows(2).all(|w| w[0] <= w[1]), || format!("{d:?}"));
    let traco: f64 = (0..n).map(|i| a[i * n + i]).sum();
    t.near(&format!("{what}: sum of eigenvalues = trace"), d.iter().sum(), traco, 1e-9 * norma * n as f64);
    let j = jacobi(a, n);
    let dif = d.iter().zip(&j).map(|(x, y)| (x - y).abs()).fold(0.0f64, f64::max);
    t.check(&format!("{what}: eigenvalues = Jacobi's"), dif <= 1e-9 * norma * n as f64, || format!("max difference {dif:e}: {:?} vs {:?}", &d[..n.min(5)], &j[..n.min(5)]));
}

fn testa_autos(t: &mut T, r: &mut Rng) {
    for &n in &[1usize, 2, 3, 5, 10, 40, 120] {
        let mut a = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..=i {
                let x = r.normal();
                a[i * n + j] = x;
                a[j * n + i] = x;
            }
        }
        confere_autos(t, &format!("random symmetric {n}x{n}"), &a, n);
    }
    // A Gram matrix (positive semidefinite, rank 7): what the PCA feeds it.
    let (n, m) = (60, 7);
    let x: Vec<f64> = (0..n * m).map(|_| r.normal()).collect();
    let mut g = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            g[i * n + j] = (0..m).map(|c| x[i * m + c] * x[j * m + c]).sum();
        }
    }
    confere_autos(t, "Gram matrix of rank 7", &g, n);
    let (d, _) = symmetric_eigen(&g, n).unwrap();
    t.check("Gram of rank 7: 53 zero eigenvalues", d[..53].iter().all(|x| x.abs() < 1e-9 * d[n - 1]), || format!("{:?}", &d[50..56]));
    // Special matrices.
    confere_autos(t, "zero 4x4", &[0.0; 16], 4);
    let mut diag = vec![0.0; 25];
    for (i, v) in [3.0, 1.0, 3.0, -2.0, 1.0].iter().enumerate() {
        diag[i * 5 + i] = *v;
    }
    confere_autos(t, "diagonal with repeats", &diag, 5);
    let (d, _) = symmetric_eigen(&diag, 5).unwrap();
    let dif = d.iter().zip([-2.0, 1.0, 1.0, 3.0, 3.0]).map(|(a, b)| (a - b).abs()).fold(0.0f64, f64::max);
    t.check("diagonal: eigenvalues -2, 1, 1, 3, 3", dif < 1e-14, || format!("{d:?}"));
    let u = [1.0, 2.0, -1.0, 0.5, 3.0, 0.0];
    let mut a = vec![0.0; 36];
    for i in 0..6 {
        for j in 0..6 {
            a[i * 6 + j] = u[i] * u[j] + if i == j { 1.0 } else { 0.0 };
        }
    }
    confere_autos(t, "identity + rank one", &a, 6);
    let (d, _) = symmetric_eigen(&a, 6).unwrap();
    t.near("identity + rank one: largest = 1 + |u|^2", d[5], 1.0 + u.iter().map(|x| x * x).sum::<f64>(), 1e-12);
    t.check("identity + rank one: the others are 1", d[..5].iter().all(|x| (x - 1.0).abs() < 1e-12), || format!("{d:?}"));
    let (d, v) = symmetric_eigen(&[2.0, 1.0, 1.0, 2.0], 2).unwrap();
    t.check("2x2: eigenvalues 1 and 3", (d[0] - 1.0).abs() < 1e-14 && (d[1] - 3.0).abs() < 1e-14, || format!("{d:?}"));
    t.near("2x2: eigenvector of 3 is (1,1)/sqrt 2", (v[2] * v[3]).abs(), 0.5, 1e-12);
    t.fails("wrong length", symmetric_eigen(&[1.0, 2.0, 3.0], 2), &["3 values", "2 x 2"]);
    t.fails("not finite", symmetric_eigen(&[1.0, f64::NAN, f64::NAN, 1.0], 2), &["not finite"]);
}

// ---------------------------------------------------------------------
//  2. THE PCA AND THE PROJECTION
// ---------------------------------------------------------------------

fn testa_pca(t: &mut T, r: &mut Rng) {
    // 40 individuals of two groups, 300 sites, a few missing genotypes.
    let (n, m, k) = (40usize, 300usize, 3usize);
    let fa: Vec<(f64, f64)> = (0..m).map(|_| (0.05 + 0.9 * r.f(), 0.05 + 0.9 * r.f())).collect();
    let mut g = vec![0u8; n * m];
    for i in 0..n {
        for j in 0..m {
            let p = if i < 20 { fa[j].0 } else { fa[j].1 };
            g[i * m + j] = u8::from(r.chance(p)) + u8::from(r.chance(p));
            if r.chance(0.02) {
                g[i * m + j] = 9;
            }
        }
    }
    let fit = pca_fit(&g, n, m, k).unwrap();
    // The same computed here.
    let mut x = vec![0.0; n * m];
    let mut medias = vec![0.0; m];
    for j in 0..m {
        let vistos: Vec<f64> = (0..n).filter(|&i| g[i * m + j] <= 2).map(|i| f64::from(g[i * m + j])).collect();
        let soma: f64 = vistos.iter().sum();
        let media = soma / vistos.len() as f64;
        let p = (1.0 + soma) / (2.0 + 2.0 * vistos.len() as f64);
        let sd = (2.0 * p * (1.0 - p)).sqrt();
        medias[j] = media;
        for i in 0..n {
            x[i * m + j] = if g[i * m + j] <= 2 { (f64::from(g[i * m + j]) - media) / sd } else { 0.0 };
        }
    }
    let dm = medias.iter().zip(&fit.mean).map(|(a, b)| (a - b).abs()).fold(0.0f64, f64::max);
    t.check("PCA: site means over the observed genotypes", dm < 1e-12, || format!("max difference {dm:e}"));
    let mut gram = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            gram[i * n + j] = (0..m).map(|c| x[i * m + c] * x[j * m + c]).sum::<f64>() / m as f64;
        }
    }
    let ev = jacobi(&gram, n);
    for c in 0..k {
        t.check(&format!("PCA: eigenvalue {} = Jacobi's of the Gram matrix built here", c + 1), (fit.eigenvalues[c] - ev[n - 1 - c]).abs() < 1e-5 * ev[n - 1], || {
            format!("{} vs {} (f32 storage)", fit.eigenvalues[c], ev[n - 1 - c])
        });
    }
    t.near("PCA: total variance = trace", fit.total_variance, (0..n).map(|i| gram[i * n + i]).sum(), 1e-4 * n as f64);
    // Loadings: unit length, orthogonal; scores = X v / sqrt(m).
    let mut orto: f64 = 0.0;
    for a in 0..k {
        for b in 0..k {
            let p: f64 = (0..m).map(|j| fit.loadings[j * k + a] * fit.loadings[j * k + b]).sum();
            orto = orto.max((p - if a == b { 1.0 } else { 0.0 }).abs());
        }
    }
    t.check("PCA: axes of unit length and orthogonal", orto < 1e-5, || format!("max |V'V - I| {orto:e}"));
    let mut dif: f64 = 0.0;
    for i in 0..n {
        for c in 0..k {
            let s: f64 = (0..m).map(|j| x[i * m + j] * fit.loadings[j * k + c]).sum::<f64>() / (m as f64).sqrt();
            dif = dif.max((s - fit.scores[i * k + c]).abs());
        }
    }
    t.check("PCA: scores = X v / sqrt(sites)", dif < 1e-5, || format!("max difference {dif:e}"));
    // Projecting a fitted individual gives back its score (missing ones count as the mean).
    let mut dif: f64 = 0.0;
    for i in 0..n {
        let linha: Vec<u8> = g[i * m..i * m + m].iter().map(|&x| if x == 9 { 255 } else { x }).collect();
        let y: Vec<f64> = (0..m).map(|j| if linha[j] <= 2 { x[i * m + j] } else { 0.0 }).collect();
        let (s, _) = fit.project_standardized(&y).unwrap();
        for c in 0..k {
            dif = dif.max((s[c] - fit.scores[i * k + c]).abs());
        }
    }
    t.check("PCA: a fitted individual projects onto its own score", dif < 1e-5, || format!("max difference {dif:e}"));
    // Group 1 vs group 2 on PC1.
    let m1: f64 = (0..20).map(|i| fit.scores[i * k]).sum::<f64>() / 20.0;
    let m2: f64 = (20..40).map(|i| fit.scores[i * k]).sum::<f64>() / 20.0;
    t.check("PCA: PC1 separates the two groups", (m1 - m2).abs() > 0.5, || format!("means {m1}, {m2}"));

    // Least squares with missing sites: a vector exactly on the axes is
    // recovered whatever is missing; mean imputation alone shrinks it.
    let b = [0.7, -0.4, 0.25];
    let y: Vec<f64> = (0..m).map(|j| (0..k).map(|c| b[c] * fit.loadings[j * k + c]).sum::<f64>() * (m as f64).sqrt()).collect();
    let (s, vistos) = fit.project_standardized(&y).unwrap();
    t.check("projection, nothing missing: exact", (0..k).all(|c| (s[c] - b[c]).abs() < 1e-9), || format!("{s:?}"));
    t.eq("projection: sites observed", vistos, m);
    let mut ym = y.clone();
    let mut faltam = 0;
    for v in ym.iter_mut() {
        if r.chance(0.6) {
            *v = f64::NAN;
            faltam += 1;
        }
    }
    let (s, vistos) = fit.project_standardized(&ym).unwrap();
    t.eq("projection, 60% missing: sites observed", vistos, m - faltam);
    t.check("projection, 60% missing: still exact (least squares on the observed sites)", (0..k).all(|c| (s[c] - b[c]).abs() < 1e-9), || format!("{s:?} vs {b:?}"));
    let ingenua: Vec<f64> = (0..k)
        .map(|c| (0..m).filter(|&j| !ym[j].is_nan()).map(|j| ym[j] * fit.loadings[j * k + c]).sum::<f64>() / (m as f64).sqrt())
        .collect();
    let razao = ingenua.iter().zip(&b).map(|(x, y)| x * y).sum::<f64>() / b.iter().map(|y| y * y).sum::<f64>();
    t.check("projection, 60% missing: plain mean imputation would shrink it to about 40%", (0.25..0.55).contains(&razao), || format!("ratio {razao}"));
    let quase: Vec<f64> = (0..m).map(|j| if j < 3 { y[j] } else { f64::NAN }).collect();
    t.eq("projection: too few sites -> None", fit.project_standardized(&quase).is_none(), true);
    let gq: Vec<u8> = (0..m).map(|j| if j % 3 == 0 { g[j] } else { 255 }).collect();
    let (_, vq) = fit.project(&gq).unwrap();
    t.eq("project(): 255 is missing", vq, gq.iter().filter(|&&x| x <= 2).count());

    t.fails("PCA: 0 components", pca_fit(&g, n, m, 0), &["at least 1"]);
    t.fails("PCA: too few individuals", pca_fit(&g[..3 * m], 3, m, 3), &["4 individuals"]);
    t.fails("PCA: wrong size", pca_fit(&g[..10], n, m, 2), &["genotypes for"]);
    // Two identical individuals and nothing else: one direction only.
    let dois: Vec<u8> = g[..m].iter().chain(g[..m].iter()).chain(g[m..2 * m].iter()).copied().collect();
    t.fails("PCA: more components than the data hold", pca_fit(&dois, 3, m, 2), &["no variance"]);
}

// ---------------------------------------------------------------------
//  3. THE WHOLE TOOL ON SIMULATED POPULATIONS
// ---------------------------------------------------------------------

fn roda(panel: &str, vcfs: &[String], o: &PcaOptions) -> Result<Value, String> {
    let s = ancestry_pca(panel, vcfs, o)?;
    serde_json::from_str(&s).map_err(|e| format!("the JSON does not parse: {e}"))
}

fn testa_ferramenta(t: &mut T, r: &mut Rng) {
    let grupos = ["A", "B", "C"];
    let f = 0.08;
    let painel = painel_bn(r, &grupos, &["1", "2", "3"], 6000, 3000, f);
    let caminho = escreve_painel(t, "bn3.tsv.gz", &painel);
    // Realised F_ST (Hudson) between the pairs, for the record.
    for (a, b) in [(0, 1), (0, 2), (1, 2)] {
        let (mut num, mut den) = (0.0, 0.0);
        for p in &painel.p {
            num += (p[a] - p[b]).powi(2);
            den += p[a] * (1.0 - p[b]) + p[b] * (1.0 - p[a]);
        }
        println!("    simulated panel: F_ST (Hudson) {}-{} = {:.4}", grupos[a], grupos[b], num / den);
    }

    let mut nomes = Vec::new();
    let mut gts = Vec::new();
    let mut verdade = Vec::new();
    for (g, nome) in grupos.iter().enumerate() {
        for i in 0..5 {
            let mut q = vec![0.0; 3];
            q[g] = 1.0;
            nomes.push(format!("{nome}{}", i + 1));
            gts.push(pessoa(r, &painel, &q, 0.0));
            verdade.push(q);
        }
    }
    for (nome, q) in [("AB1", [0.5, 0.5, 0.0]), ("AB2", [0.5, 0.5, 0.0]), ("BC1", [0.0, 0.5, 0.5]), ("ABC", [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0])] {
        nomes.push(nome.to_string());
        gts.push(pessoa(r, &painel, &q, 0.0));
        verdade.push(q.to_vec());
    }
    nomes.push("A_missing".to_string());
    gts.push(pessoa(r, &painel, &[1.0, 0.0, 0.0], 0.7));
    verdade.push(vec![1.0, 0.0, 0.0]);
    // "chr" names in the VCF, none in the panel.
    let vcf = escreve_vcf(t, "people.vcf.gz", "chr", &painel, &nomes, &gts, false);
    let tabela = t.p("pca_sites.tsv");
    let o = PcaOptions { max_sites: 3000, per_group: 100, components: 2, seed: 7, sites_output: Some(tabela.clone()), ..Default::default() };
    // The pure and admixed people (not the one with sites missing): sites
    // observed by all of them are plenty.
    let sem_falta: Vec<String> = nomes.iter().filter(|n| *n != "A_missing").cloned().collect();
    let o1 = PcaOptions { samples: Some(sem_falta.clone()), ..o.clone() };
    let t0 = Instant::now();
    let j = match roda(&caminho, &[vcf.clone()], &o1) {
        Ok(j) => j,
        Err(e) => {
            t.check("simulated: runs", false, || e);
            return;
        }
    };
    println!("    simulated, {} people x {} sites x 300 reference: {:.2} s", sem_falta.len(), j["sites_used"], t0.elapsed().as_secs_f64());
    for chave in ["groups", "sites_used", "variance_explained", "reference", "centroids", "people", "warnings", "text", "shrinkage", "candidate_sites", "min_people_per_site", "spacing_bp"] {
        t.check(&format!("JSON has `{chave}`"), !j[chave].is_null(), String::new);
    }
    t.eq("groups", j["groups"].clone(), serde_json::json!(["A", "B", "C"]));
    let m = j["sites_used"].as_u64().unwrap() as usize;
    t.check("sites used: at most max_sites, and close to it", (2700..=3000).contains(&m), || format!("{m}"));
    t.eq("every site observed in all the people", j["min_people_per_site"].as_u64(), Some(sem_falta.len() as u64));
    let ve = pcs(&j["variance_explained"]);
    t.check("variance explained: PC1 >= PC2, both structure (well above the noise, ~0.006)", ve[0] >= ve[1] && ve[1] > 0.02, || format!("{ve:?}"));
    let reference = j["reference"].as_array().unwrap();
    t.eq("reference: 3 x 100 simulated individuals", reference.len(), 300);
    t.eq("reference: 100 per group, in group order", reference.iter().filter(|x| x["group"] == "B").count(), 100);
    t.eq("reference: 2 components each", pcs(&reference[150]["pcs"]).len(), 2);
    let centros: Vec<Vec<f64>> = j["centroids"].as_array().unwrap().iter().map(|c| pcs(&c["pcs"])).collect();
    let dps: Vec<Vec<f64>> = j["centroids"].as_array().unwrap().iter().map(|c| pcs(&c["sd"])).collect();
    // Clouds centred on their centroids.
    for g in 0..3 {
        for c in 0..2 {
            let v: Vec<f64> = reference[g * 100..g * 100 + 100].iter().map(|x| pcs(&x["pcs"])[c]).collect();
            let media = v.iter().sum::<f64>() / 100.0;
            t.check(&format!("cloud {} PC{}: centred on the centroid", grupos[g], c + 1), (media - centros[g][c]).abs() < 4.0 * dps[g][c] / 10.0 + 1e-4, || {
                format!("cloud mean {media}, centroid {}, sd {}", centros[g][c], dps[g][c])
            });
        }
    }
    let pessoas = j["people"].as_array().unwrap();
    t.eq("people: one entry each, in order", pessoas.iter().map(|p| p["name"].as_str().unwrap().to_string()).collect::<Vec<_>>(), sem_falta.clone());
    let pos_de = |nome: &str| -> Vec<f64> { pcs(&pessoas.iter().find(|p| p["name"] == nome).unwrap()["pcs"]) };
    let mut razoes = Vec::new();
    for (i, nome) in nomes.iter().enumerate().take(15) {
        let p = pessoas.iter().find(|p| p["name"] == nome.as_str()).unwrap();
        let g = i / 5;
        t.eq(&format!("{nome}: closest cloud is its group"), p["closest_group"].as_str(), Some(grupos[g]));
        let d = p["distance"].as_f64().unwrap();
        t.check(&format!("{nome}: inside the cloud (< 4 spreads)"), d < 4.0, || format!("distance {d}"));
        let x = pcs(&p["pcs"]);
        razoes.push(x.iter().zip(&centros[g]).map(|(a, b)| a * b).sum::<f64>() / centros[g].iter().map(|b| b * b).sum::<f64>());
    }
    let rm = razoes.iter().sum::<f64>() / razoes.len() as f64;
    t.check("people are not pulled toward 0 compared with the clouds' centres (mean ratio ~1)", (0.95..1.05).contains(&rm), || format!("mean ratio {rm}: {razoes:?}"));
    let sh = pcs(&j["shrinkage"]);
    t.check("shrinkage of a projection: between 0.5 and 1, PC1 the least", sh.iter().all(|x| (0.5..=1.02).contains(x)) && sh[0] >= sh[1] - 0.02, || format!("{sh:?}"));
    // Admixed people at the mixture of the centres.
    let dab = dist(&centros[0], &centros[1]);
    for (nome, q) in [("AB1", [0.5, 0.5, 0.0]), ("AB2", [0.5, 0.5, 0.0]), ("BC1", [0.0, 0.5, 0.5]), ("ABC", [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0])] {
        let esperado: Vec<f64> = (0..2).map(|c| (0..3).map(|g| q[g] * centros[g][c]).sum()).collect();
        let d = dist(&pos_de(nome), &esperado);
        t.check(&format!("{nome}: at the mixture of the centres (within 15% of the A-B distance)"), d < 0.15 * dab, || {
            format!("{:?} vs expected {esperado:?}: {d} (A-B {dab})", pos_de(nome))
        });
    }
    // The distance between centres, against the one computed from the
    // panel's frequencies at the sites used (about 4 F_ST). The axes of a
    // finite simulated sample catch a share cos^2 of the true between-group
    // plane: with sites m, individuals n and a spike lambda (in noise units),
    // cos^2 = (1 - g / (lambda - 1)^2) / (1 + g / (lambda - 1)), g = m / n
    // (the spiked covariance model). Here lambda - 1 ~ m F / 2 per axis.
    let tab = le_sitios(&tabela, 2);
    t.eq("sites table: one row per site used", tab.pos.len(), m);
    let razao_dist = |cs: &[Vec<f64>], tab: &TabelaSitios| -> Vec<f64> {
        let m = tab.pos.len();
        [(0usize, 1usize), (0, 2), (1, 2)]
            .iter()
            .map(|&(a, b)| {
                let teoria: f64 = (0..m).map(|jj| ((2.0 * tab.af[jj][a] - 2.0 * tab.af[jj][b]) / tab.dp[jj]).powi(2)).sum::<f64>() / m as f64;
                dist(&cs[a], &cs[b]).powi(2) / teoria
            })
            .collect()
    };
    let r100 = razao_dist(&centros, &tab);
    let lambda = m as f64 * f / 2.0;
    let gama = m as f64 / 300.0;
    let cos2 = (1.0 - gama / (lambda * lambda)) / (1.0 + gama / lambda);
    println!("    squared distances between centres / from the frequencies, 100 per group: {r100:.4?} (spiked model: ~{cos2:.3})");
    t.check("distances between centres: the share the spiked model predicts (+/- 0.04)", r100.iter().all(|x| (x - cos2).abs() < 0.04), || format!("{r100:?} vs {cos2}"));
    let o400 = PcaOptions { per_group: 400, sites_output: Some(t.p("pca_sites400.tsv")), ..o1.clone() };
    let j400 = roda(&caminho, &[vcf.clone()], &o400).unwrap();
    let c400: Vec<Vec<f64>> = j400["centroids"].as_array().unwrap().iter().map(|c| pcs(&c["pcs"])).collect();
    let r400 = razao_dist(&c400, &le_sitios(&t.p("pca_sites400.tsv"), 2));
    let cos2_400 = (1.0 - gama / 4.0 / (lambda * lambda)) / (1.0 + gama / 4.0 / lambda);
    println!("    the same with 400 per group: {r400:.4?} (spiked model: ~{cos2_400:.3})");
    t.check("distances between centres, 400 per group: closer to the whole (+/- 0.03 of the model, and above 100 per group)", r400.iter().zip(&r100).all(|(x, y)| (x - cos2_400).abs() < 0.03 && x > y), || {
        format!("{r400:?} vs {cos2_400}; 100 per group {r100:?}")
    });
    // The projection redone here from the table.
    let ajuste = tab.ajuste(2);
    let idx: HashMap<(String, u64), usize> = (0..painel.pos.len()).map(|jj| (painel.crom[jj].clone(), painel.pos[jj])).zip(0..).collect();
    for nome in ["A1", "B3", "AB1"] {
        let i = nomes.iter().position(|n| n == nome).unwrap();
        let g: Vec<u8> = (0..m).map(|jj| gts[i][idx[&(tab.crom[jj].clone(), tab.pos[jj])]].unwrap()).collect();
        let mut s = vec![0.0; 2];
        for jj in 0..m {
            for c in 0..2 {
                s[c] += (f64::from(g[jj]) - tab.media[jj]) / tab.dp[jj] * tab.v[jj][c];
            }
        }
        let s: Vec<f64> = s.iter().map(|x| x / (m as f64).sqrt()).collect();
        let d = dist(&s, &pos_de(nome));
        t.check(&format!("{nome}: projection = sum((g - mean) / sd x axis) / sqrt(sites), from the table"), d < 2e-5, || format!("{s:?} vs {:?}", pos_de(nome)));
        let (s2, _) = ajuste.project(&g).unwrap();
        t.check(&format!("{nome}: the same through PcaFit::project with the published axes"), dist(&s2, &pos_de(nome)) < 2e-5, || format!("{s2:?}"));
    }
    let mut orto: f64 = 0.0;
    for a in 0..2 {
        for b in 0..2 {
            let p: f64 = (0..m).map(|jj| tab.v[jj][a] * tab.v[jj][b]).sum();
            orto = orto.max((p - if a == b { 1.0 } else { 0.0 }).abs());
        }
    }
    t.check("sites table: the axes have unit length and are orthogonal", orto < 1e-4, || format!("{orto:e}"));
    let mut tab_ord = true;
    for w in 0..m.saturating_sub(1) {
        if tab.crom[w] == tab.crom[w + 1] && tab.pos[w + 1] < tab.pos[w] + j["spacing_bp"].as_u64().unwrap() {
            tab_ord = false;
        }
    }
    t.check("sites: sorted and at least spacing_bp apart", tab_ord, String::new);
    let por_crom: Vec<usize> = ["1", "2", "3"].iter().map(|c| tab.crom.iter().filter(|x| x == c).count()).collect();
    t.check("sites: spread over the three chromosomes", por_crom.iter().all(|&x| x > m / 4), || format!("{por_crom:?}"));
    t.check("text: says the clouds are simulated", j["text"].as_str().unwrap().contains("SIMULATED from the panel's allele frequencies"), String::new);

    // Same seed: same JSON; another seed: another reference, same answers.
    let a1 = ancestry_pca(&caminho, &[vcf.clone()], &o1).unwrap();
    let a2 = ancestry_pca(&caminho, &[vcf.clone()], &o1).unwrap();
    t.check("deterministic: the same JSON twice", a1 == a2, String::new);
    let o2 = PcaOptions { seed: 8, sites_output: None, ..o1.clone() };
    let j2 = roda(&caminho, &[vcf.clone()], &o2).unwrap();
    t.check("another seed: another reference", j2["reference"] != j["reference"], String::new);
    let mesmos = j2["people"].as_array().unwrap().iter().zip(pessoas).all(|(a, b)| a["closest_group"] == b["closest_group"]);
    t.check("another seed: the same closest groups", mesmos, String::new);
    let c2: Vec<Vec<f64>> = j2["centroids"].as_array().unwrap().iter().map(|c| pcs(&c["pcs"])).collect();
    let dc = (0..3).map(|g| dist(&c2[g], &centros[g])).fold(0.0f64, f64::max);
    t.check("another seed: centroids move by little (< 10% of the A-B distance)", dc < 0.1 * dab, || format!("{dc}"));

    // A person with 70% of the genotypes missing, with only one other: the
    // sites observed in both are too few, so sites observed in one are let
    // in, and the person is placed on the 30% observed.
    let o3 = PcaOptions {
        samples: Some(vec!["A1".to_string(), "A_missing".to_string()]),
        max_sites: 20_000,
        min_distance_bp: 0,
        sites_output: None,
        ..o.clone()
    };
    let j3 = roda(&caminho, &[vcf.clone()], &o3).unwrap();
    t.eq("missing: sites observed in at least 1 person", j3["min_people_per_site"].as_u64(), Some(1));
    let pm = &j3["people"][1];
    let mf = pm["missing_fraction"].as_f64().unwrap();
    t.check("missing: about 70% missing", (0.67..0.73).contains(&mf), || format!("{mf}"));
    t.eq("missing: still closest to A", pm["closest_group"].as_str(), Some("A"));
    let c3: Vec<Vec<f64>> = j3["centroids"].as_array().unwrap().iter().map(|c| pcs(&c["pcs"])).collect();
    let x = pcs(&pm["pcs"]);
    let razao = x.iter().zip(&c3[0]).map(|(a, b)| a * b).sum::<f64>() / c3[0].iter().map(|b| b * b).sum::<f64>();
    t.check("missing: not pulled toward 0 (ratio to the centre ~1)", (0.85..1.15).contains(&razao), || format!("ratio {razao}"));
    let d = pm["distance"].as_f64().unwrap();
    t.check("missing: near its centre (less precise: < 6 spreads)", d < 6.0, || format!("{d}"));
    let avisos = j3["warnings"].to_string();
    t.check("missing: the notes say so", avisos.contains("A_missing: 70% of the sites are missing") && avisos.contains("observed in at least 1"), || avisos.clone());

    // GIAB style: only the non-reference genotypes listed, and a callable
    // BED (without "chr") covering chromosome 1 and half of chromosome 2.
    let b_pessoa = pessoa(r, &painel, &[0.0, 1.0, 0.0], 0.0);
    let vcf_giab = escreve_vcf(t, "giab_like.vcf", "chr", &painel, &["HGB".to_string()], &[b_pessoa.clone()], true);
    let meio2 = painel.pos[6000 + 3000];
    let bed = t.write("giab_like.bed", &["1\t0\t100000000".to_string(), format!("2\t0\t{meio2}")]);
    let o4 = PcaOptions { samples: Some(vec!["".to_string(), "A1".to_string()]), callable_beds: Some(vec![bed.clone(), "".to_string()]), sites_output: Some(t.p("giab_sites.tsv")), ..o.clone() };
    let j4 = roda(&caminho, &[vcf_giab.clone(), vcf.clone()], &o4).unwrap();
    let p4 = &j4["people"][0];
    t.eq("GIAB style: the sample's name", p4["name"].as_str(), Some("HGB"));
    t.eq("GIAB style: closest to B", p4["closest_group"].as_str(), Some("B"));
    let tab4 = le_sitios(&t.p("giab_sites.tsv"), 2);
    let dentro = (0..tab4.pos.len()).filter(|&jj| tab4.crom[jj] == "1" || (tab4.crom[jj] == "2" && tab4.pos[jj] <= meio2)).count();
    t.eq("GIAB style: observed = the sites inside the BED (absent = 0/0 there)", p4["sites_observed"].as_u64(), Some(dentro as u64));
    t.eq("GIAB style: A1 (full VCF) observed everywhere", j4["people"][1]["sites_observed"].as_u64(), Some(tab4.pos.len() as u64));
    // The GIAB-style genotypes match: project here with the table.
    let idx4: Vec<usize> = (0..tab4.pos.len()).map(|jj| idx[&(tab4.crom[jj].clone(), tab4.pos[jj])]).collect();
    let g4: Vec<u8> = (0..tab4.pos.len())
        .map(|jj| if tab4.crom[jj] == "1" || (tab4.crom[jj] == "2" && tab4.pos[jj] <= meio2) { b_pessoa[idx4[jj]].unwrap() } else { 255 })
        .collect();
    let (s4, _) = tab4.ajuste(2).project(&g4).unwrap();
    t.check("GIAB style: the position = least squares redone here on the BED's sites", dist(&s4, &pcs(&p4["pcs"])) < 2e-5, || format!("{s4:?} vs {}", p4["pcs"]));

    // Several VCFs, a single-sample one without a name.
    let vcf_a = escreve_vcf(t, "one_a.vcf", "", &painel, &["Solo".to_string()], &[gts[0].clone()], false);
    let j5 = roda(&caminho, &[vcf_a.clone(), vcf_giab.clone()], &PcaOptions { samples: None, callable_beds: Some(vec!["".to_string(), bed.clone()]), ..o1.clone() }).unwrap();
    let nomes5: Vec<&str> = j5["people"].as_array().unwrap().iter().map(|p| p["name"].as_str().unwrap()).collect();
    t.eq("several VCFs: each file's only sample", nomes5, vec!["Solo", "HGB"]);
    t.eq("several VCFs: Solo is A", j5["people"][0]["closest_group"].as_str(), Some("A"));
    // All the samples of one VCF when none is named.
    let j6 = roda(&caminho, &[vcf.clone()], &PcaOptions { samples: None, sites_output: None, ..o1.clone() }).unwrap();
    t.eq("one VCF, no names: every sample", j6["people"].as_array().unwrap().len(), nomes.len());

    // Components beyond groups - 1 are noise, and the notes say so.
    let j7 = roda(&caminho, &[vcf.clone()], &PcaOptions { components: 4, sites_output: None, ..o1.clone() }).unwrap();
    let ve7 = pcs(&j7["variance_explained"]);
    t.check("4 components of 3 groups: PC3 and PC4 are noise (far below PC2)", ve7[2] < 0.3 * ve7[1] && ve7[3] <= ve7[2], || format!("{ve7:?}"));
    t.check("4 components of 3 groups: a note", j7["warnings"].to_string().contains("differ along only 2 direction(s)"), || j7["warnings"].to_string());
    t.eq("4 components: the same closest groups", j7["people"].as_array().unwrap().iter().map(|p| p["closest_group"].clone()).collect::<Vec<_>>(), pessoas.iter().map(|p| p["closest_group"].clone()).collect::<Vec<_>>());

    // Errors.
    let ruim = |o: PcaOptions| ancestry_pca(&caminho, &[vcf.clone()], &o);
    t.fails("no panel", ancestry_pca(&t.p("none.tsv"), &[vcf.clone()], &o1), &["Error opening"]);
    let um = t.write("one_group.tsv", &["chrom\tpos\tref\talt\tAF_X".to_string(), "1\t100\tA\tG\t0.3".to_string()]);
    t.fails("single group", ancestry_pca(&um, &[vcf.clone()], &o1), &["single group"]);
    t.fails("unknown sample", ruim(PcaOptions { samples: Some(vec!["Nobody".to_string()]), ..o1.clone() }), &["no sample named `Nobody`"]);
    t.fails("reference too large", ruim(PcaOptions { per_group: 1000, ..o1.clone() }), &["at most 2000", "per_group"]);
    t.fails("0 components", ruim(PcaOptions { components: 0, ..o1.clone() }), &["components is 0"]);
    t.fails("max_sites too small", ruim(PcaOptions { max_sites: 2, ..o1.clone() }), &["max_sites is 2"]);
    t.fails("min_maf out of range", ruim(PcaOptions { min_maf: 0.5, ..o1.clone() }), &["min_maf"]);
    t.fails("sites table over the panel", ruim(PcaOptions { sites_output: Some(caminho.clone()), ..o1.clone() }), &["also one of the input files"]);
    t.fails("callable BEDs: wrong count", ruim(PcaOptions { callable_beds: Some(vec![bed.clone(), bed.clone()]), ..o1.clone() }), &["callable_beds", "one per person"]);
    t.fails("several VCFs: wrong number of names", ancestry_pca(&caminho, &[vcf.clone(), vcf_a.clone()], &PcaOptions { samples: Some(vec!["A1".to_string()]), ..o1.clone() }), &["one per VCF"]);
    t.fails("no VCF", ancestry_pca(&caminho, &[], &o1), &["No VCF given"]);
    let vcf_longe = escreve_vcf(t, "elsewhere.vcf", "Z", &painel, &["Far".to_string()], &[gts[0].clone()], false);
    t.fails("no shared chromosome", ancestry_pca(&caminho, &[vcf_longe], &PcaOptions { absent_is_ref: false, ..o.clone() }), &["None of the", "same genome build"]);
    let alto = t.write("high_maf.tsv", &["chrom\tpos\tref\talt\tAF_A\tAF_B".to_string(), "1\t100\tA\tG\t0.01\t0.02".to_string()]);
    t.fails("nothing above the MAF", ancestry_pca(&alto, &[vcf.clone()], &o1), &["no biallelic SNV", "below the MAF"]);
    // A VCF out of order.
    let mut desordem = le_texto(&vcf_a);
    let n = desordem.len();
    desordem.swap(n - 1, n - 5);
    let vcf_des = t.write("unsorted.vcf", &desordem);
    t.fails("unsorted VCF", ancestry_pca(&caminho, &[vcf_des], &PcaOptions { samples: None, ..o1.clone() }), &["must be sorted by position"]);
}

// ---------------------------------------------------------------------
//  4. TIME: 5 people x 5,000 sites x 500 reference individuals
// ---------------------------------------------------------------------

fn testa_tempo(t: &mut T, r: &mut Rng) {
    let grupos = ["AFR", "AMR", "EAS", "EUR", "SAS"];
    let painel = painel_bn(r, &grupos, &["1", "2"], 25_000, 1500, 0.05);
    let caminho = escreve_painel(t, "bn5.tsv.gz", &painel);
    let mut nomes = Vec::new();
    let mut gts = Vec::new();
    for (g, n) in grupos.iter().enumerate() {
        let mut q = vec![0.0; 5];
        q[g] = 1.0;
        nomes.push(format!("{n}_person"));
        gts.push(pessoa(r, &painel, &q, 0.0));
    }
    let vcf = escreve_vcf(t, "people5.vcf.gz", "", &painel, &nomes, &gts, false);
    let o = PcaOptions { max_sites: 5000, per_group: 100, components: 4, ..Default::default() };
    let _ = ancestry_pca(&caminho, &[vcf.clone()], &o);
    let mut tempos = Vec::new();
    let mut j = Value::Null;
    for _ in 0..3 {
        let t0 = Instant::now();
        j = roda(&caminho, &[vcf.clone()], &o).unwrap();
        tempos.push(t0.elapsed().as_secs_f64());
    }
    tempos.sort_by(|a, b| a.partial_cmp(b).unwrap());
    println!(
        "    5 people x {} sites x 500 reference individuals, 4 components (panel of {} SNVs): best of 3 {:.2} s",
        j["sites_used"], painel.pos.len(), tempos[0]
    );
    let usados = j["sites_used"].as_u64().unwrap();
    t.check("time run: 4950 to 5000 sites (the even spacing lands just under max_sites)", (4950..=5000).contains(&usados), || format!("{usados}"));
    let certos = j["people"].as_array().unwrap().iter().zip(grupos).filter(|(p, g)| p["closest_group"] == *g).count();
    t.eq("time run: every person in its own cloud", certos, 5);
    t.check("time run: under 10 s", tempos[0] < 10.0, || format!("{tempos:?}"));
    // The core alone.
    let (n, m) = (500usize, 5000usize);
    let g: Vec<u8> = (0..n * m).map(|_| r.below(3) as u8).collect();
    let t0 = Instant::now();
    let fit = pca_fit(&g, n, m, 4).unwrap();
    let tf = t0.elapsed().as_secs_f64();
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let x = r.normal();
            a[i * n + j] = x;
            a[j * n + i] = x;
        }
    }
    let t0 = Instant::now();
    symmetric_eigen(&a, n).unwrap();
    let te = t0.elapsed().as_secs_f64();
    let t0 = Instant::now();
    let (dd, _) = symmetric_eigen(&vec![1.0; n * n], n).unwrap();
    let t1 = t0.elapsed().as_secs_f64();
    println!("    pca_fit 500 x 5000: {tf:.2} s; a 500 x 500 eigendecomposition {te:.2} s (all-ones matrix {t1:.2} s); top eigenvalue of noise {:.3}", fit.eigenvalues[0]);
    t.near("all-ones 500x500: eigenvalue 500", dd[n - 1], 500.0, 1e-9);
    t.check("all-ones 500x500: the others 0", dd[..n - 1].iter().all(|x| x.abs() < 1e-10), || format!("{:?}", &dd[..3]));
    t.check("all-ones 500x500: not slower than a random matrix (no iterating on negligible blocks)", t1 < 2.0 * te + 0.05, || format!("{t1} s vs {te} s"));
}

// ---------------------------------------------------------------------
//  REAL DATA (not part of the pass/fail count)
// ---------------------------------------------------------------------

/// Reads the dosage of `amostras` at the positions of `sitios` (one
/// chromosome) from a VCF; 255 = missing or not found.
fn le_genotipos(vcf: &str, sitios: &TabelaSitios, amostras: &[String]) -> Vec<Vec<u8>> {
    let alvo: HashMap<u64, usize> = sitios.pos.iter().enumerate().map(|(i, &p)| (p, i)).collect();
    let m = sitios.pos.len();
    let mut out = vec![vec![255u8; m]; amostras.len()];
    let mut r = externo::abre_leitura(vcf).unwrap();
    let mut col: Vec<usize> = Vec::new();
    let mut linha = String::new();
    let mut achados = 0;
    loop {
        linha.clear();
        if r.read_line(&mut linha).unwrap_or(0) == 0 {
            break;
        }
        if linha.starts_with("##") {
            continue;
        }
        if linha.starts_with("#CHROM") {
            let nomes: Vec<&str> = linha.trim_end().split('\t').skip(9).collect();
            let idx: HashMap<&str, usize> = nomes.iter().enumerate().map(|(i, n)| (*n, i)).collect();
            col = amostras.iter().map(|a| idx[a.as_str()]).collect();
            continue;
        }
        let mut it = linha.splitn(3, '\t');
        let _ = it.next();
        let pos: u64 = match it.next().and_then(|x| x.parse().ok()) {
            Some(p) => p,
            None => continue,
        };
        let Some(&j) = alvo.get(&pos) else { continue };
        let c: Vec<&str> = linha.trim_end().split('\t').collect();
        if c[3] != sitios.ra[j].0 || c[4] != sitios.ra[j].1 {
            continue;
        }
        achados += 1;
        for (s, &ci) in col.iter().enumerate() {
            let gt = c[9 + ci].split(':').next().unwrap_or(".");
            let b = gt.as_bytes();
            out[s][j] = if b.len() == 3 && b[0].is_ascii_digit() && b[2].is_ascii_digit() { (b[0] - b'0') + (b[2] - b'0') } else { 255 };
        }
    }
    println!("    read {achados} of {m} sites for {} samples from {vcf}", amostras.len());
    out
}

/// The 1000 Genomes and GIAB run. Environment:
///   MMG_HUMAN_GENOMES  the human_genomes folder (1000genomes/...)
///   MMG_HUMAN_DEMO     the unzipped human demo pack (the folder holding human/)
///   MMG_OUT            where to write the panel and the tables
fn real() {
    let base = std::env::var("MMG_HUMAN_GENOMES").unwrap_or_else(|_| "human_genomes".to_string());
    let kg = format!("{base}/1000genomes");
    let vcf22 = format!("{kg}/nygc_30x_phased_panel/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz");
    let pop = format!("{kg}/nygc_30x_metadata/20130606_g1k_3202_samples_ped_population.txt");
    let ped = format!("{kg}/nygc_30x_metadata/1kGP.3202_samples.pedigree_info.txt");
    let out = std::env::var("MMG_OUT").unwrap_or_else(|_| std::env::temp_dir().join("mmg_pca_real").to_string_lossy().to_string());
    std::fs::create_dir_all(&out).unwrap();
    let k = 4;

    // ---- held-out people and the panel without them (children out too)
    let mut filhos: HashSet<String> = HashSet::new();
    for l in std::fs::read_to_string(&ped).unwrap().lines().skip(1) {
        let c: Vec<&str> = l.split_whitespace().collect();
        if c.len() >= 3 && (c[1] != "0" || c[2] != "0") {
            filhos.insert(c[0].to_string());
        }
    }
    let mut por_pop: Vec<(String, Vec<(String, String)>)> = Vec::new();
    let mut superpop: HashMap<String, String> = HashMap::new();
    for l in std::fs::read_to_string(&pop).unwrap().lines().skip(1) {
        let c: Vec<&str> = l.split_whitespace().collect();
        superpop.insert(c[1].to_string(), c[6].to_string());
        if filhos.contains(c[1]) {
            continue;
        }
        match por_pop.iter_mut().find(|(p, _)| p == c[5]) {
            Some((_, v)) => v.push((c[1].to_string(), c[6].to_string())),
            None => por_pop.push((c[5].to_string(), vec![(c[1].to_string(), c[6].to_string())])),
        }
    }
    por_pop.sort_by(|a, b| a.0.cmp(&b.0));
    let mut fora: Vec<(String, String, String)> = Vec::new();
    for (p, v) in &por_pop {
        for (s, sp) in v.iter().take(5) {
            fora.push((s.clone(), p.clone(), sp.clone()));
        }
    }
    let painel = format!("{out}/panel_chr22_heldout.tsv.gz");
    if !Path::new(&painel).exists() {
        let mut excluir: Vec<String> = filhos.iter().cloned().collect();
        excluir.extend(fora.iter().map(|x| x.0.clone()));
        let t0 = Instant::now();
        match ancestralidade::build_panel(&[vcf22.clone()], &pop, "SampleID", "Superpopulation", &excluir, &painel, &ancestralidade::PanelOptions::default()) {
            Ok(s) => println!("panel in {:.1} s\n{}", t0.elapsed().as_secs_f64(), s.resumo()),
            Err(e) => {
                println!("panel failed: {e}");
                return;
            }
        }
    }

    // ---- the held-out people projected
    let tabela = format!("{out}/sites_1kg_heldout.tsv");
    let o = PcaOptions {
        samples: Some(fora.iter().map(|x| x.0.clone()).collect()),
        absent_is_ref: false,
        sites_output: Some(tabela.clone()),
        components: k,
        ..Default::default()
    };
    let t0 = Instant::now();
    let j: Value = match ancestry_pca(&painel, &[vcf22.clone()], &o) {
        Ok(s) => {
            std::fs::write(format!("{out}/pca_heldout.json"), &s).unwrap();
            serde_json::from_str(&s).unwrap()
        }
        Err(e) => {
            println!("PCA failed: {e}");
            return;
        }
    };
    println!("PCA of {} held-out 1000 Genomes people in {:.1} s (one pass over the chr22 VCF)", fora.len(), t0.elapsed().as_secs_f64());
    let grupos: Vec<String> = j["groups"].as_array().unwrap().iter().map(|x| x.as_str().unwrap().to_string()).collect();
    println!("groups {grupos:?}; sites {}; spacing {} bp; variance explained {}; shrinkage {}", j["sites_used"], j["spacing_bp"], j["variance_explained"], j["shrinkage"]);
    let pessoas = j["people"].as_array().unwrap();
    let mut certo: HashMap<String, (usize, usize)> = HashMap::new();
    let mut amr: HashMap<String, usize> = HashMap::new();
    let mut dists: HashMap<String, Vec<f64>> = HashMap::new();
    for (p, (s, popu, sp)) in pessoas.iter().zip(&fora) {
        assert_eq!(p["name"].as_str().unwrap(), s);
        let cg = p["closest_group"].as_str().unwrap().to_string();
        let e = certo.entry(sp.clone()).or_default();
        e.1 += 1;
        if cg == *sp {
            e.0 += 1;
        } else {
            println!("    {s} ({popu}, {sp}) -> closest {cg}, distance {:.2}", p["distance"].as_f64().unwrap());
        }
        if sp == "AMR" {
            *amr.entry(format!("{popu}->{cg}")).or_default() += 1;
        }
        dists.entry(sp.clone()).or_default().push(p["distance"].as_f64().unwrap());
    }
    let mut sps: Vec<&String> = certo.keys().collect();
    sps.sort();
    for sp in sps {
        let (a, n) = certo[sp];
        let mut d = dists[sp].clone();
        d.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!("  {sp}: {a}/{n} closest to their own group; distance median {:.2}, max {:.2}", d[d.len() / 2], d[d.len() - 1]);
    }
    let mut amr_v: Vec<_> = amr.into_iter().collect();
    amr_v.sort();
    println!("  AMR breakdown: {amr_v:?}");

    // ---- the PCA position against the bar chart (estimate_ancestry)
    let persons: Vec<ancestralidade::Person> =
        fora.iter().map(|x| ancestralidade::Person { vcf: vcf22.clone(), sample: Some(x.0.clone()), callable: None }).collect();
    let oa = ancestralidade::AncestryOptions { bootstrap: 0, absent_is_ref: false, block_bp: 2_000_000, ..Default::default() };
    let t0 = Instant::now();
    match ancestralidade::estimate_ancestry(&painel, &persons, &oa) {
        Ok(res) => {
            println!("estimate_ancestry of the same {} people in {:.1} s", res.len(), t0.elapsed().as_secs_f64());
            let cs: Vec<Vec<f64>> = j["centroids"].as_array().unwrap().iter().map(|c| pcs(&c["pcs"])).collect();
            let mut qp: Vec<Vec<f64>> = Vec::new();
            for p in pessoas {
                qp.push(mistura(&cs, &pcs(&p["pcs"])));
            }
            for (g, nome) in grupos.iter().enumerate() {
                let a: Vec<f64> = res.iter().map(|x| x.proportions[g]).collect();
                let b: Vec<f64> = qp.iter().map(|x| x[g]).collect();
                let mad = a.iter().zip(&b).map(|(x, y)| (x - y).abs()).sum::<f64>() / a.len() as f64;
                println!("  {nome}: proportion from the bar chart vs mixture of the PCA centres: r = {:.3}, mean |difference| {:.3}", correlacao(&a, &b), mad);
            }
            for nome in ["NA19625", "HG03714"] {
                if let Some(i) = fora.iter().position(|x| x.0 == nome) {
                    let a: Vec<String> = res[i].proportions.iter().map(|x| format!("{x:.2}")).collect();
                    let b: Vec<String> = qp[i].iter().map(|x| format!("{x:.2}")).collect();
                    println!("    {nome}: bar chart {a:?}; PCA mixture {b:?}");
                }
            }
        }
        Err(e) => println!("estimate_ancestry failed: {e}"),
    }

    // ---- real 1000 Genomes individuals on the same sites
    let tab = le_sitios(&tabela, k);
    let mut r = Rng(20261116);
    let fora_set: HashSet<&String> = fora.iter().map(|x| &x.0).collect();
    let mut reais: Vec<(String, String)> = Vec::new();
    for g in &grupos {
        let mut cand: Vec<String> = por_pop.iter().flat_map(|(_, v)| v.iter()).filter(|(s, sp)| sp == g && !fora_set.contains(s)).map(|(s, _)| s.clone()).collect();
        // Deterministic shuffle, then 100.
        for i in (1..cand.len()).rev() {
            let j = r.below(i + 1);
            cand.swap(i, j);
        }
        for s in cand.into_iter().take(100) {
            reais.push((s, g.clone()));
        }
    }
    let mut todas: Vec<String> = reais.iter().map(|x| x.0.clone()).collect();
    todas.extend(fora.iter().map(|x| x.0.clone()));
    let t0 = Instant::now();
    let gen = le_genotipos(&vcf22, &tab, &todas);
    println!("    ({:.1} s)", t0.elapsed().as_secs_f64());
    let m = tab.pos.len();
    let nr = reais.len();
    let mat: Vec<u8> = gen[..nr].iter().flat_map(|v| v.iter().copied()).collect();
    let fit_real = pca_fit(&mat, nr, m, k).unwrap();
    let ve_real: Vec<f64> = fit_real.eigenvalues.iter().map(|x| x / fit_real.total_variance).collect();
    println!("real 1000G PCA ({nr} individuals x {m} sites): variance explained {:?}", ve_real.iter().map(|x| format!("{:.4}", x)).collect::<Vec<_>>());
    // Centres: simulated (projection of 2p) and real (in-sample means).
    let c_sim: Vec<Vec<f64>> = j["centroids"].as_array().unwrap().iter().map(|c| pcs(&c["pcs"])).collect();
    let sd_sim: Vec<Vec<f64>> = j["centroids"].as_array().unwrap().iter().map(|c| pcs(&c["sd"])).collect();
    let mut c_real = vec![vec![0.0; k]; grupos.len()];
    for (i, (_, g)) in reais.iter().enumerate() {
        let gi = grupos.iter().position(|x| x == g).unwrap();
        for c in 0..k {
            c_real[gi][c] += fit_real.scores[i * k + c] / 100.0;
        }
    }
    for c in 0..k {
        let a: Vec<f64> = c_sim.iter().map(|v| v[c]).collect();
        let row: Vec<String> = (0..k)
            .map(|d| {
                let b: Vec<f64> = c_real.iter().map(|v| v[d]).collect();
                format!("{:+.3}", correlacao(&a, &b))
            })
            .collect();
        println!("  correlation of the group centres: simulated PC{} vs real PC1..PC{k}: {}", c + 1, row.join(" "));
    }
    let a2: Vec<Vec<f64>> = c_sim.iter().map(|v| v[..2].to_vec()).collect();
    let b2: Vec<Vec<f64>> = c_real.iter().map(|v| v[..2].to_vec()).collect();
    println!("  Procrustes correlation of the 5 group centres, PC1-PC2: {:.4}; PC1-PC{k}: {:.4}", procrustes(&a2, &b2), procrustes(&c_sim, &c_real));
    // The held-out people on both sets of axes.
    let sim_h: Vec<Vec<f64>> = pessoas.iter().map(|p| pcs(&p["pcs"])).collect();
    let real_h: Vec<Vec<f64>> = gen[nr..].iter().map(|g| fit_real.project(g).unwrap().0).collect();
    let s2: Vec<Vec<f64>> = sim_h.iter().map(|v| v[..2].to_vec()).collect();
    let r2: Vec<Vec<f64>> = real_h.iter().map(|v| v[..2].to_vec()).collect();
    println!(
        "  Procrustes correlation of the {} held-out people, simulated axes vs real-individual axes: PC1-PC2 {:.4}; PC1-PC{k} {:.4}",
        sim_h.len(),
        procrustes(&s2, &r2),
        procrustes(&sim_h, &real_h)
    );
    // Our own projection from the table = the JSON.
    let ajuste = tab.ajuste(k);
    let dmax = gen[nr..].iter().zip(&sim_h).map(|(g, s)| dist(&ajuste.project(g).unwrap().0, s)).fold(0.0f64, f64::max);
    println!("  held-out positions recomputed from the sites table: max difference {dmax:.2e}");
    // Real individuals projected onto the simulated axes, against the clouds.
    let proj_reais: Vec<Vec<f64>> = gen[..nr].iter().map(|g| ajuste.project(g).unwrap().0).collect();
    for (gi, g) in grupos.iter().enumerate() {
        let v: Vec<&Vec<f64>> = proj_reais.iter().zip(&reais).filter(|(_, x)| &x.1 == g).map(|(p, _)| p).collect();
        let mut linha = format!("  real {g} on the simulated axes:");
        for c in 0..k {
            let media = v.iter().map(|p| p[c]).sum::<f64>() / v.len() as f64;
            let sd = (v.iter().map(|p| (p[c] - media).powi(2)).sum::<f64>() / (v.len() as f64 - 1.0)).sqrt();
            linha.push_str(&format!(" PC{} centre {:+.4} (simulated {:+.4}, {:+.1} sd), sd x{:.2};", c + 1, media, c_sim[gi][c], (media - c_sim[gi][c]) / sd_sim[gi][c], sd / sd_sim[gi][c]));
        }
        println!("{linha}");
    }
    // Nearest centre of the real reference individuals (projected) as a sanity check.
    let mut acerto = 0;
    for (p, (_, g)) in proj_reais.iter().zip(&reais) {
        let pooled: Vec<f64> = (0..k).map(|c| (sd_sim.iter().map(|s| s[c] * s[c]).sum::<f64>() / grupos.len() as f64).sqrt()).collect();
        let d: Vec<f64> = c_sim.iter().map(|cg| (0..k).map(|c| ((p[c] - cg[c]) / pooled[c]).powi(2)).sum::<f64>().sqrt()).collect();
        let best = (0..grupos.len()).fold(0, |a, b| if d[b] < d[a] { b } else { a });
        if &grupos[best] == g {
            acerto += 1;
        }
    }
    println!("  real reference individuals (in the panel) closest to their own simulated cloud: {acerto}/{nr}");
    // Files for the Python cross-check.
    let mut w = String::new();
    for (i, (s, g)) in reais.iter().enumerate() {
        w.push_str(&format!("{s}\t{g}\t{}\n", gen[i].iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\t")));
    }
    std::fs::write(format!("{out}/real_reference_genotypes.tsv"), w).unwrap();
    let mut w = String::new();
    for (i, (s, p, g)) in fora.iter().enumerate() {
        w.push_str(&format!("{s}\t{p}\t{g}\t{}\n", gen[nr + i].iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\t")));
    }
    std::fs::write(format!("{out}/heldout_genotypes.tsv"), w).unwrap();
    let mut w = format!("eigenvalues\t{}\ntotal\t{}\n", fit_real.eigenvalues.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\t"), fit_real.total_variance);
    for (i, h) in real_h.iter().enumerate() {
        w.push_str(&format!("heldout\t{}\t{}\n", fora[i].0, h.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\t")));
    }
    std::fs::write(format!("{out}/rust_real_fit.tsv"), w).unwrap();

    // ---- GIAB, with the demo pack's panel
    if let Ok(demo) = std::env::var("MMG_HUMAN_DEMO") {
        let h = format!("{demo}/human");
        let ids = ["HG001", "HG002", "HG003", "HG004", "HG005", "HG006", "HG007"];
        let vcfs: Vec<String> = ids.iter().map(|i| format!("{h}/giab/{i}.vcf.gz")).collect();
        let beds: Vec<String> = ids.iter().map(|i| format!("{h}/giab/{i}.callable.bed")).collect();
        for (nome_painel, max) in [("ancestry_panel_chr22_1000g.tsv.gz", 5000usize), ("ancestry_panel_genomewide_1000g.tsv.gz", 5000)] {
            let o = PcaOptions { callable_beds: Some(beds.clone()), max_sites: max, components: k, sites_output: Some(format!("{out}/giab_{nome_painel}.sites.tsv")), ..Default::default() };
            let t0 = Instant::now();
            match ancestry_pca(&format!("{h}/{nome_painel}"), &vcfs, &o) {
                Ok(s) => {
                    let el = t0.elapsed().as_secs_f64();
                    std::fs::write(format!("{out}/giab_{nome_painel}.json"), &s).unwrap();
                    let j: Value = serde_json::from_str(&s).unwrap();
                    println!("GIAB x {nome_painel}: {:.2} s for 7 people x {} sites x 500 reference", el, j["sites_used"]);
                    print!("{}", j["text"].as_str().unwrap());
                }
                Err(e) => println!("GIAB x {nome_painel} failed: {e}"),
            }
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(String::as_str) == Some("real") {
        real();
        return;
    }
    let dir = std::env::temp_dir().join(format!("trial_pca_ancestry_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut t = T { dir: dir.clone(), total: 0, ok: 0, falhas: Vec::new() };
    let mut r = Rng(20261116);
    let t0 = Instant::now();
    testa_autos(&mut t, &mut r);
    println!("  eigensolver: {} checks so far, {:.1} s", t.total, t0.elapsed().as_secs_f64());
    testa_pca(&mut t, &mut r);
    println!("  PCA and projection: {} checks so far, {:.1} s", t.total, t0.elapsed().as_secs_f64());
    testa_ferramenta(&mut t, &mut r);
    println!("  the tool on simulated populations: {} checks so far, {:.1} s", t.total, t0.elapsed().as_secs_f64());
    testa_tempo(&mut t, &mut r);
    let _ = std::fs::remove_dir_all(&dir);
    println!("pca_ancestralidade tests: {} passed, {} failed ({:.1} s)", t.ok, t.total - t.ok, t0.elapsed().as_secs_f64());
    if t.ok != t.total {
        std::process::exit(1);
    }
}
