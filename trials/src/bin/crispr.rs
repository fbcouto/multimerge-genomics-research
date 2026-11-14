#![allow(dead_code)]
// Tests of src/crispr.rs: guide counting and the MAGeCK-style screen test.
//
//  - counting is compared with a brute force written here: every window of
//    every read against every guide, as strings, for exact matches and for
//    one mismatch (guides that collide at one mismatch included);
//  - the statistics are compared with hand-worked values, with exact
//    formulas (the beta CDF against the binomial sum, the random generator
//    against the Park-Miller check value) and with the numbers MAGeCK's own
//    Python functions give on the same small inputs;
//  - the screen test runs on simulated screens: planted essential (depleted)
//    and enriched genes must be found, and a screen without effects must
//    give calibrated gene p-values.
//
// With the MAGeCK data folder as first argument it also counts the demo2
// FASTQs and tests demo1 and the leukemia table, checking the numbers MAGeCK
// printed for them and the report MAGeCK ships in demo5; with a folder of
// MAGeCK 0.5.9.5 outputs as second argument it compares guide by guide and
// gene by gene.
//
// Usage: cargo run --release --bin crispr [-- <mageck data dir> [<mageck outputs dir>]]
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
#[path = "../../../src/crispr.rs"]
mod crispr;

use crispr::*;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

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
    /// Relative difference (absolute when `want` is 0) at most `tol`.
    fn perto(&mut self, got: f64, want: f64, tol: f64, what: &str) -> bool {
        let d = if want == 0.0 { got.abs() } else { ((got - want) / want).abs() };
        let ok = d <= tol || (got.is_nan() && want.is_nan());
        if !ok {
            println!("FAIL  {what}\n   got:  {got:e}\n   want: {want:e} (relative difference {d:e}, tolerance {tol:e})");
        }
        self.check(ok, what)
    }
    fn erro<X: Debug>(&mut self, r: Result<X, String>, agulhas: &[&str], what: &str) {
        match r {
            Ok(x) => {
                println!("FAIL  {what}: should have failed, got {:?}", format!("{x:?}").chars().take(200).collect::<String>());
                self.check(false, what);
            }
            Err(e) => {
                let falta: Vec<&&str> = agulhas.iter().filter(|a| !e.contains(**a)).collect();
                if !falta.is_empty() {
                    println!("FAIL  {what}: the message lacks {falta:?}\n   message: {e}");
                }
                self.check(falta.is_empty(), what);
            }
        }
    }
}

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
    /// Uniform in (0, 1).
    fn u(&mut self) -> f64 {
        (self.next() as f64 + 0.5) / (1u64 << 31) as f64
    }
    fn normal(&mut self) -> f64 {
        (-2.0 * self.u().ln()).sqrt() * (2.0 * std::f64::consts::PI * self.u()).cos()
    }
    /// Marsaglia-Tsang.
    fn gamma(&mut self, k: f64, theta: f64) -> f64 {
        if k < 1.0 {
            return self.gamma(k + 1.0, theta) * self.u().powf(1.0 / k);
        }
        let d = k - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let x = self.normal();
            let v = (1.0 + c * x).powi(3);
            if v <= 0.0 {
                continue;
            }
            let u = self.u();
            if u.ln() < 0.5 * x * x + d - d * v + d * v.ln() {
                return d * v * theta;
            }
        }
    }
    fn poisson(&mut self, l: f64) -> u64 {
        if l < 30.0 {
            let lim = (-l).exp();
            let (mut k, mut p) = (0u64, 1.0);
            loop {
                p *= self.u();
                if p <= lim {
                    return k;
                }
                k += 1;
            }
        }
        (l + l.sqrt() * self.normal()).round().max(0.0) as u64
    }
    /// Negative binomial with this mean and dispersion (var = m + phi m^2).
    fn nb(&mut self, m: f64, phi: f64) -> u64 {
        let l = self.gamma(1.0 / phi, m * phi);
        self.poisson(l)
    }
    fn dna(&mut self, n: usize) -> String {
        (0..n).map(|_| b"ACGT"[self.below(4) as usize] as char).collect()
    }
}

fn outra_base(b: u8, r: &mut Rng) -> u8 {
    loop {
        let x = b"ACGT"[r.below(4) as usize];
        if x != b.to_ascii_uppercase() {
            return x;
        }
    }
}

fn rc(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}

// ---------------------------------------------------------------- formats

fn testa_formatos(t: &mut T) {
    // Python's "{:.5g}" and str(): the number formats of MAGeCK's tables.
    let g5 = [
        (0.00012345, "0.00012345"),
        (123456.0, "1.2346e+05"),
        (1.0, "1"),
        (-0.86379123, "-0.86379"),
        (4.6563e-69, "4.6563e-69"),
        (99999.5, "1e+05"),
        (0.0, "0"),
        (1e-5, "1e-05"),
        (2.5e-5, "2.5e-05"),
        (1234.5678, "1234.6"),
        (0.5, "0.5"),
        (-1e-7, "-1e-07"),
        (17.52512345, "17.525"),
        (1e16, "1e+16"),
        (3.68e-06, "3.68e-06"),
    ];
    for (x, s) in g5 {
        t.eq(&py_g(x, 5), &s.to_string(), &format!("{{:.5g}} of {x}"));
    }
    t.eq(&py_g(0.00012345, 4), &"0.0001234".to_string(), "{:.4g} of 0.00012345");
    t.eq(&py_g(1234.5678, 4), &"1235".to_string(), "{:.4g} of 1234.5678");
    let rep = [
        (3.6799e-06, "3.6799e-06"),
        (0.000123, "0.000123"),
        (1.0, "1.0"),
        (1e-5, "1e-05"),
        (100.0, "100.0"),
        (5.489e-05, "5.489e-05"),
        (0.005489, "0.005489"),
        (123456789012345680.0, "1.2345678901234568e+17"),
        (0.1, "0.1"),
        (123456.0, "123456.0"),
        (17.52512345, "17.52512345"),
    ];
    for (x, s) in rep {
        t.eq(&py_repr(x), &s.to_string(), &format!("str() of {x}"));
    }
}

// ---------------------------------------------------------------- normal distribution

fn testa_normal(t: &mut T) {
    // Python's math.erfc at the same points.
    let py = [
        (0.0, 1.0),
        (0.1, 0.8875370839817152),
        (0.5, 0.4795001221869534),
        (0.99, 0.1614919304446302),
        (1.0, 0.1572992070502851),
        (1.5, 0.033894853524689274),
        (2.0, 0.004677734981047266),
        (3.0, 2.2090496998585434e-05),
        (5.0, 1.537459794428035e-12),
        (10.0, 2.088487583762545e-45),
        (20.0, 5.395865611607901e-176),
        (26.0, 5.663192408856143e-296),
        (-0.3, 1.3286267594591274),
        (-2.0, 1.9953222650189528),
    ];
    for (x, e) in py {
        t.perto(erfc(x), e, 2e-14, &format!("erfc({x}) against Python"));
    }
    // continuity where the series hands over to the continued fraction
    let (a, b) = (erfc(1.0 - 1e-12), erfc(1.0 + 1e-12));
    t.check((a - b).abs() < 1e-12, "erfc is continuous at x = 1");
    // the quantile inverts the CDF, in both tails
    let mut pior = 0f64;
    for k in 1..=200 {
        let lp = -(k as f64) * 0.35; // ln p from -0.35 to -70
        let x = qnorm_log(lp);
        pior = pior.max(((norm_cdf(x)).ln() - lp).abs() / lp.abs());
        let q = 1.0 - (-0.001 * k as f64).exp(); // p near 1
        let y = qnorm_log((1.0 - q).ln());
        pior = pior.max((norm_sf(y) - q).abs() / q);
    }
    t.check(pior < 1e-12, &format!("qnorm inverts the normal CDF in both tails (worst relative error {pior:e})"));
    t.perto(qnorm_log(0.5f64.ln()), 0.0, 1e-15, "qnorm(0.5) = 0");
    t.perto(qnorm_log(0.975f64.ln()), 1.959963984540054, 1e-14, "qnorm(0.975) = 1.95996...");
    // MAGeCK's interpolated table (QNormConverter.get_qnorm with islog=True)
    let mageck: [(f64, Option<f64>); 9] = [
        (-24.0, Some(-6.5092982359543)),
        (-23.9995, Some(-6.509223119514205)),
        (-10.1234, Some(-3.943625657502798)),
        (-0.7, Some(-0.00855947858248099)),
        (-0.0015, Some(2.9844289533586355)),
        (-0.001, Some(3.090380786916421)),
        (-0.0005, None),
        (-25.0, None),
        (-3.0, Some(-1.6469217205277105)),
    ];
    for (lp, e) in mageck {
        let g = crispr::internals::qnorm_table(lp);
        match (g, e) {
            (Some(g), Some(e)) => {
                t.perto(g, e, 1e-12, &format!("MAGeCK's qnorm table at ln p = {lp}"));
            }
            _ => {
                t.eq(&g.is_some(), &e.is_some(), &format!("MAGeCK's qnorm table at ln p = {lp}: outside the table"));
            }
        }
    }
}

// ---------------------------------------------------------------- RRA math

/// I_x(a, b) for integer a, b: P(Binomial(a + b - 1, x) >= a), summed exactly.
fn beta_binomial(a: usize, b: usize, x: f64) -> f64 {
    let n = a + b - 1;
    let mut s = 0.0;
    for j in a..=n {
        // C(n, j) in logs
        let lc: f64 = (1..=n).map(|i| (i as f64).ln()).sum::<f64>()
            - (1..=j).map(|i| (i as f64).ln()).sum::<f64>()
            - (1..=(n - j)).map(|i| (i as f64).ln()).sum::<f64>();
        s += (lc + j as f64 * x.ln() + (n - j) as f64 * (1.0 - x).ln()).exp();
    }
    s
}

/// The lo-value by its definition, independently: min over the sorted
/// percentiles kept of P(the (i+1)-th smallest of n uniforms <= x).
fn lo_ingenuo(pct: &[f64], alpha: f64) -> (f64, usize) {
    let mut v = pct.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    let (mut lo, mut bons) = (1.0f64, 0);
    for (i, &x) in v.iter().enumerate() {
        if x > alpha && i > 0 {
            break;
        }
        if x <= alpha {
            bons += 1;
        }
        lo = lo.min(beta_binomial(i + 1, n - i, x));
    }
    (lo, bons)
}

fn bh_ingenuo(p: &[f64]) -> Vec<f64> {
    let n = p.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| p[a].partial_cmp(&p[b]).unwrap());
    let mut out = vec![0.0; n];
    for (r, &i) in idx.iter().enumerate() {
        // min over the ranks at or after r of p * n / rank
        let m = (r..n).map(|s| p[idx[s]] * n as f64 / (s + 1) as f64).fold(f64::INFINITY, f64::min);
        out[i] = m.min(1.0);
    }
    out
}

fn testa_rra_matematica(t: &mut T, r: &mut Rng) {
    // beta CDF (RRA's AS 63 with its LogGamma) against the exact binomial sum
    let mut pior = 0f64;
    for _ in 0..400 {
        let a = 1 + r.below(30) as usize;
        let b = 1 + r.below(30) as usize;
        let x = r.u();
        pior = pior.max((beta_cdf(a as f64, b as f64, x) - beta_binomial(a, b, x)).abs());
    }
    t.check(pior < 1e-9, &format!("beta CDF = binomial sum on 400 random (a, b, x) (worst difference {pior:e})"));
    t.perto(beta_cdf(1.0, 3.0, 0.01), 1.0 - 0.99f64.powi(3), 1e-10, "Beta(1, 3) CDF at 0.01 = 1 - 0.99^3");

    // lo-value, worked by hand
    let (lo, bons) = lo_value(&[0.01, 0.2, 0.5], 0.25, -1.0);
    t.perto(lo, 0.029701, 1e-9, "lo of [0.01, 0.2, 0.5], alpha 0.25: Beta(1,3) at 0.01 = 0.029701");
    t.eq(&bons, &2, "two good guides below alpha 0.25");
    let (lo, bons) = lo_value(&[0.6, 0.7], 0.25, -1.0);
    t.perto(lo, 0.84, 1e-10, "lo of [0.6, 0.7]: the first percentile always counts (1 - 0.4^2)");
    t.eq(&bons, &0, "no good guide when all are above alpha");
    let (lo, _) = lo_value(&[0.2, 0.21, 0.22], 0.25, -1.0);
    t.perto(lo, beta_binomial(3, 1, 0.22).min(beta_binomial(2, 2, 0.21)).min(beta_binomial(1, 3, 0.2)), 1e-10, "lo takes the minimum over all kept guides");
    let mut pior = 0f64;
    let mut bons_ok = true;
    for _ in 0..300 {
        let n = 1 + r.below(12) as usize;
        let alpha = 0.05 + 0.45 * r.u();
        let mut v: Vec<f64> = (0..n).map(|_| r.u()).collect();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let (a, ba) = lo_value(&v, alpha, -1.0);
        let (b, bb) = lo_ingenuo(&v, alpha);
        pior = pior.max((a - b).abs());
        bons_ok &= ba == bb;
    }
    t.check(pior < 1e-9 && bons_ok, &format!("lo-value = its definition on 300 random genes (worst {pior:e})"));

    // percentiles: mid-rank over 2n, ties sharing the middle
    let lista = [1.0, 2.0, 2.0, 3.0];
    t.eq(&crispr::internals::percentile(1.0, &lista), &0.125, "percentile of the smallest of 4: 0.5/4");
    t.eq(&crispr::internals::percentile(2.0, &lista), &0.5, "two tied values share the middle rank");
    t.eq(&crispr::internals::percentile(3.0, &lista), &0.875, "percentile of the largest of 4: 3.5/4");
    // permutation p-value: lo below every random value -> 0.5/M
    let aleat: Vec<f64> = (1..=10).map(|i| i as f64 / 10.0).collect();
    t.eq(&crispr::internals::perm_p(0.01, &aleat), &0.05, "p of a lo below all 10 random lo-values = 0.5/10");
    t.eq(&crispr::internals::perm_p(0.5, &aleat), &0.45, "p of a lo equal to the 5th random value = 4.5/10");

    // Park & Miller: from seed 1, 10,000 steps reach 399268537 (rngs.c's CHECK)
    let mut s = 1u64;
    for _ in 0..10_000 {
        crispr::internals::park_miller(&mut s);
    }
    t.eq(&s, &399_268_537, "Park-Miller generator: 10,000 steps from 1 give the check value 399268537");
    let mut ok = true;
    for _ in 0..50 {
        let n = r.below(100_000);
        let semente = 1 + r.below(2_000_000_000);
        let mut a = semente;
        for _ in 0..n {
            crispr::internals::park_miller(&mut a);
        }
        ok &= a == crispr::internals::park_miller_jump(semente, n);
    }
    t.check(ok, "jumping the generator ahead = stepping it");
    // random lo-values computed in parallel blocks = one sequential stream
    for (tam, total, nulo) in [(3usize, 10_000usize, false), (1, 5000, false), (5, 9000, true)] {
        let ctrl: Vec<f64> = (0..37).map(|i| (i as f64 + 0.5) / 37.0).collect();
        let mut s1 = 123456u64;
        let par = crispr::internals::random_lo(&mut s1, tam, total, 0.2, if nulo { Some(&ctrl) } else { None });
        let mut s2 = 123456u64;
        let mut seq: Vec<f64> = (0..total)
            .map(|_| {
                let mut v: Vec<f64> = (0..tam)
                    .map(|_| {
                        let u = crispr::internals::park_miller(&mut s2);
                        if nulo {
                            ctrl[((ctrl.len() as f64 * u) as usize).min(ctrl.len() - 1)]
                        } else {
                            u
                        }
                    })
                    .collect();
                v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                lo_value(&v, 0.2, -1.0).0
            })
            .collect();
        seq.sort_by(|a, b| a.partial_cmp(b).unwrap());
        t.check(par == seq && s1 == s2, &format!("{total} random lo-values of size {tam}{}: parallel blocks = RRA's sequential stream", if nulo { " from control percentiles" } else { "" }));
    }

    // Benjamini-Hochberg: worked case, MAGeCK's pFDR gives [0.04, 0.0533, 0.0533, 0.5]
    let p = [0.01, 0.04, 0.03, 0.5];
    t.eq(&bh_fdr(&p), &vec![0.04, 0.05333333333333333, 0.05333333333333333, 0.5], "BH of [0.01, 0.04, 0.03, 0.5] = MAGeCK's pFDR");
    let q = crispr::internals::bh_rra(&p);
    t.check(q.iter().zip([0.04, 0.0533333333333333, 0.0533333333333333, 0.5]).all(|(a, b)| (a - b).abs() < 1e-15), "RRA's BH of the same list");
    let mut ok = true;
    for _ in 0..100 {
        let n = 1 + r.below(60) as usize;
        let v: Vec<f64> = (0..n).map(|_| if r.chance(20) { (r.below(5) as f64) / 5.0 } else { r.u() }).collect();
        let e = bh_ingenuo(&v);
        ok &= bh_fdr(&v).iter().zip(&e).all(|(a, b)| (a - b).abs() < 1e-12);
        ok &= crispr::internals::bh_rra(&v).iter().zip(&e).all(|(a, b)| (a - b).abs() < 1e-12);
    }
    t.check(ok, "both BH versions = min over later ranks of p n / rank, on 100 random lists with ties");
}

// ---------------------------------------------------------------- normalization and model

fn testa_normalizacao(t: &mut T) {
    // Factors MAGeCK's normalizeCounts gives on the same tables.
    let cols = vec![vec![10.0, 30.0, 0.0, 40.0, 7.0], vec![20.0, 30.0, 5.0, 10.0, 7.0]];
    let (f, m, _) = crispr::internals::size_factors(&cols, NormMethod::Median).unwrap();
    t.eq(&m, &"median".to_string(), "median normalization kept");
    t.perto(f[0], 1.1428571428571426, 1e-14, "median factor of sample 1 = MAGeCK (8/7)");
    t.perto(f[1], 1.0333333333333332, 1e-14, "median factor of sample 2 = MAGeCK (31/30)");
    let (f, _, _) = crispr::internals::size_factors(&cols, NormMethod::Total).unwrap();
    t.perto(f[0], 0.9137931034482759, 1e-14, "total factor of sample 1 = MAGeCK");
    t.perto(f[1], 1.1041666666666667, 1e-14, "total factor of sample 2 = MAGeCK");
    let z = vec![vec![0.0, 0.0, 0.0, 3.0, 4.0], vec![5.0, 6.0, 7.0, 8.0, 0.0]];
    let (f, m, av) = crispr::internals::size_factors(&z, NormMethod::Median).unwrap();
    t.eq(&m, &"total".to_string(), "60% zero-count guides: median normalization falls back to total, as in MAGeCK");
    t.perto(f[0], 2.357142857142857, 1e-14, "fallback factor of sample 1 = MAGeCK");
    t.perto(f[1], 0.6346153846153846, 1e-14, "fallback factor of sample 2 = MAGeCK");
    t.check(av.iter().any(|a| a.contains("switched to total")), "the fallback is reported");
    let (f, m, _) = crispr::internals::size_factors(&cols, NormMethod::None).unwrap();
    t.eq(&(f, m), &(vec![1.0, 1.0], "none".to_string()), "no normalization");
    t.erro(crispr::internals::size_factors(&[vec![1.0, 2.0], vec![0.0, 0.0]], NormMethod::Total), &["s1", "no reads"], "a sample without reads cannot be normalized");

    // Gini index, as mageckcount_gini on ln(count + 1)
    t.perto(gini_index(&[0, 0, 10, 100]), 0.7720528662014617, 1e-14, "Gini of [0, 0, 10, 100] = MAGeCK");
    t.check(gini_index(&[5, 5, 5, 5]).abs() < 1e-15, "Gini of an even library is 0");

    // Median as mmedian
    t.eq(&median(&[3.0, 1.0, 2.0]), &2.0, "median of 3 values");
    t.eq(&median(&[4.0, 1.0, 3.0, 2.0]), &2.5, "median of 4 values: mean of the two middle ones");
    t.eq(&median(&[]), &0.0, "median of nothing is 0, as mmedian");

    // weighted least squares on an exact line, and the variance model on the
    // table MAGeCK's modelmeanvar fitted to (1.8276983668133715, 0)
    t.eq(&least_squares(&[1.0, 2.0, 3.0, 4.0], &[5.0, 8.0, 11.0, 14.0], &[1.0, 2.0, 0.5, 3.0]), &(3.0, 2.0), "least squares on y = 2 + 3x");
    let tabela = modelo_tabela();
    let (k, b, n, fora) = crispr::internals::variance_model(&tabela);
    t.perto(k, 1.8276983668133715, 1e-12, "variance model slope k = MAGeCK's modelmeanvar");
    t.eq(&b, &0.0, "variance model intercept clamped to 0, as MAGeCK");
    t.check(n > 10 && fora == 0, "variance model fitted on the guides whose variance exceeds the mean");
    // an exact power law in log2(x + 1) space is recovered
    let mut vals = Vec::new();
    for i in 0..40 {
        let m = 20.0 + 15.0 * i as f64;
        let alvo = 2f64.powf(0.7) * (m + 1.0).powf(1.3) - 1.0; // var - mean
        // [m - d, m, m + d]: median m, variance 2 d^2 / (3 - 1) = d^2
        let d = (m + alvo).sqrt();
        vals.push(vec![m - d, m, m + d]);
    }
    let (k, b, _, _) = crispr::internals::variance_model(&vals);
    t.check((k - 1.3).abs() < 1e-9 && (b - 0.7).abs() < 1e-9, &format!("variance model recovers k = 1.3, b = 0.7 (got {k}, {b})"));
}

/// The 40 x 3 table of refs2.py (random.seed(5)).
fn modelo_tabela() -> Vec<Vec<f64>> {
    let v: [[f64; 3]; 40] = [
        [243.60769311715217, 545.7651584812367, 313.33633897040244],
        [285.91672490288056, 774.3011151288174, 398.62081495501985],
        [268.5259048550242, 713.6627599670707, 371.2497944962657],
        [16.653588548259588, 28.37091950748832, 19.357588000389296],
        [380.0742592608495, 778.2526270918847, 471.9615749141653],
        [435.630819049648, 501.99549130625905, 450.94574341655823],
        [219.64385635806002, 295.67358629920784, 237.18917865217105],
        [226.9558306453402, 431.51427424327005, 274.16162532178555],
        [10.744357145422912, 13.982079516785516, 11.49152384650659],
        [103.93805078722366, 274.69617246975275, 143.3437711754996],
        [365.6460606108069, 445.3275583751526, 384.03409855642514],
        [382.95283216317176, 455.03752274269436, 399.5877607584462],
        [298.8316808507282, 349.99672023072594, 310.6389976307277],
        [4.341776106506934, 11.001159069413696, 5.878556790254649],
        [101.6553025103041, 132.09959875101424, 108.68090933508336],
        [362.7146742867654, 919.9110308090585, 491.298448868833],
        [105.45699924589354, 290.7029206641788, 148.2060580347286],
        [216.6218093396363, 456.2283095957679, 271.9156170910513],
        [76.33954330201905, 206.45358100158876, 106.36585969422745],
        [246.28676122922712, 682.1377595768943, 346.86776084791956],
        [407.29849356841834, 581.0809198401519, 447.40213040035684],
        [174.6387474718262, 214.28991840676002, 183.78901768758016],
        [75.61532502964768, 82.14617919927177, 77.12244522263785],
        [126.27781594192167, 247.15588467742475, 154.17275488088393],
        [5.317153140598615, 11.19961392266653, 6.674644090306595],
        [156.24103847807234, 225.6519732648629, 172.25894650579323],
        [351.01078342374944, 607.351992038802, 410.16644695030004],
        [138.0288906707818, 238.94657647134983, 161.3175873939898],
        [347.7609508497564, 373.978789899729, 353.8112213997501],
        [484.3289996018093, 498.8252309383536, 487.67428375639645],
        [280.8083231341999, 693.9492457007814, 376.14853603418027],
        [10.648297839487572, 24.927203921988536, 13.94343001237241],
        [153.93461209198412, 294.01699370751044, 186.26131554172096],
        [9.360716079129281, 9.937419464124236, 9.493801475666578],
        [67.45999392360962, 184.87232190963223, 94.55514653576869],
        [79.08961282225962, 179.57360717620378, 102.27822690393904],
        [333.71357294439156, 903.3987507220444, 465.179383200773],
        [156.79243506196923, 237.72420545881323, 175.46899746124092],
        [203.13039009982202, 470.050771720669, 264.7274012430944],
        [45.35488935731696, 102.2571062863434, 58.48617018709229],
    ];
    v.iter().map(|r| r.to_vec()).collect()
}

// ---------------------------------------------------------------- library

fn testa_biblioteca(t: &mut T) {
    std::fs::write(t.p("lib.tsv"), "g1\tACGTACGTAC\tA\ng2\tTTTTGGGGCC\tA\ng3\tacgtaaaaac\tB\n").unwrap();
    let l = read_library(&t.p("lib.tsv"), false).unwrap();
    t.eq(&(l.ids.len(), l.header, l.separator), &(3, false, '\t'), "TSV library without header");
    t.eq(&l.sequences[2], &"ACGTAAAAAC".to_string(), "sequences are upper-cased");
    t.eq(&l.gene_count(), &2, "two genes");
    std::fs::write(t.p("lib.csv"), "\u{feff}sgRNAID,Seq,gene\r\n\"g1\",\"ACGTACGTAC\",\"A\"\r\ng2,TTTTGGGGCCA,A,extra\r\n\r\ng3,ACGTAAAAAC,B\r\n").unwrap();
    let l = read_library(&t.p("lib.csv"), false).unwrap();
    t.eq(&(l.ids.clone(), l.header, l.separator), &(vec!["g1".to_string(), "g2".to_string(), "g3".to_string()], true, ','), "CSV with BOM, header, quotes, CRLF, blank line, extra column");
    t.eq(&l.lengths, &vec![(11, 1), (10, 2)], "guide lengths, longest first");
    let l = read_library(&t.p("lib.csv"), true).unwrap();
    t.eq(&l.sequences[0], &"GTACGTACGT".to_string(), "reverse complement of the library");
    std::fs::write(t.p("dup.tsv"), "g1\tAAAACCCCGG\tA\ng1\tAAAACCCCGT\tA\ng2\tAAAACCCCGG\tB\ng3\tCCCCAAAAGG\tB\n").unwrap();
    let l = read_library(&t.p("dup.tsv"), false).unwrap();
    t.eq(&(l.ids.clone(), l.duplicate_ids, l.duplicate_sequences), &(vec!["g1".to_string(), "g3".to_string()], 1, 1), "a repeated id and a repeated sequence are skipped and counted, as in MAGeCK");
    t.check(l.resumo().contains("repeated") && l.json().contains("\"duplicate_sequences\":1"), "duplicates in the summary and the JSON");

    std::fs::write(t.p("dois.tsv"), "g1\tACGT\n").unwrap();
    t.erro(read_library(&t.p("dois.tsv"), false), &["dois.tsv", "line 1", "3 -- guide id, guide sequence, gene"], "library with 2 columns");
    std::fs::write(t.p("n.tsv"), "g1\tACGTACGT\tA\ng2\tACGNACGT\tB\n").unwrap();
    t.erro(read_library(&t.p("n.tsv"), false), &["n.tsv", "line 2", "ACGNACGT", "column order"], "guide with an N");
    std::fs::write(t.p("ordem.tsv"), "g1\tA1CF\tACGTACGTAC\n").unwrap();
    let r = read_library(&t.p("ordem.tsv"), false);
    t.erro(r, &["no guides found"], "columns in the wrong order: the only line is taken as a header, and nothing is left");
    std::fs::write(t.p("semgene.tsv"), "g1\tACGTACGT\t\n").unwrap();
    t.erro(read_library(&t.p("semgene.tsv"), false), &["line 1", "no gene"], "guide without a gene");
    std::fs::write(t.p("longo.tsv"), format!("g1\t{}\tA\n", "A".repeat(65))).unwrap();
    t.erro(read_library(&t.p("longo.tsv"), false), &["65 nt", "64"], "guide longer than 64 nt");
    std::fs::write(t.p("vazio.tsv"), "\n\n").unwrap();
    t.erro(read_library(&t.p("vazio.tsv"), false), &["vazio.tsv", "no guides"], "empty library");
    std::fs::write(t.p("espaco.txt"), "g1 ACGTACGT A\n").unwrap();
    t.erro(read_library(&t.p("espaco.txt"), false), &["no tab or comma"], "space-separated library");
    t.erro(read_library(&t.p("nao_existe.tsv"), false), &["nao_existe.tsv"], "missing library file");
}

// ---------------------------------------------------------------- counting

/// What the brute force says about one read.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Ref {
    Nada,
    Exato(usize),
    Troca(usize),
    Ambiguo,
}

fn bruto(seqs: &[String], read: &str, offsets: &[usize], mm: bool) -> Ref {
    let r: Vec<u8> = read.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let mut comps: Vec<usize> = seqs.iter().map(|s| s.len()).collect::<HashSet<_>>().into_iter().collect();
    comps.sort_by(|a, b| b.cmp(a));
    let acgt = |b: u8| matches!(b, b'A' | b'C' | b'G' | b'T');
    for &o in offsets {
        for &l in &comps {
            if o + l > r.len() {
                continue;
            }
            let w = &r[o..o + l];
            if !w.iter().all(|&b| acgt(b)) {
                continue;
            }
            for (g, s) in seqs.iter().enumerate() {
                if s.as_bytes() == w {
                    return Ref::Exato(g);
                }
            }
        }
    }
    if !mm {
        return Ref::Nada;
    }
    let mut c: HashSet<usize> = HashSet::new();
    for &o in offsets {
        for &l in &comps {
            if o + l > r.len() {
                continue;
            }
            let w = &r[o..o + l];
            for (g, s) in seqs.iter().enumerate() {
                if s.len() == l && s.bytes().zip(w).filter(|(a, b)| a != *b || !acgt(**b)).count() <= 1 {
                    c.insert(g);
                }
            }
        }
    }
    match c.len() {
        0 => Ref::Nada,
        1 => Ref::Troca(*c.iter().next().unwrap()),
        _ => Ref::Ambiguo,
    }
}

/// A library with guides that collide at one mismatch, guides two
/// mismatches apart, two lengths, and a 19-mer that is the prefix of a
/// 20-mer.
fn biblioteca_teste(r: &mut Rng) -> Vec<(String, String, String)> {
    let mut v: Vec<(String, String, String)> = Vec::new();
    let mut vistos = HashSet::new();
    let mut add = |v: &mut Vec<(String, String, String)>, s: String, gene: String| {
        if vistos.insert(s.clone()) {
            let id = format!("sg{}", v.len());
            v.push((id, s, gene));
        }
    };
    for i in 0..220 {
        add(&mut v, r.dna(20), format!("G{}", i / 4));
    }
    for i in 0..40 {
        add(&mut v, r.dna(19), format!("H{}", i / 4));
    }
    for i in 0..12 {
        // one mismatch apart
        let a = r.dna(20);
        let mut b = a.clone().into_bytes();
        let p = r.below(20) as usize;
        b[p] = outra_base(b[p], r);
        add(&mut v, a, format!("C{i}"));
        add(&mut v, String::from_utf8(b).unwrap(), format!("C{i}"));
    }
    for i in 0..12 {
        // two mismatches apart
        let a = r.dna(20);
        let mut b = a.clone().into_bytes();
        let p = r.below(10) as usize;
        let q = 10 + r.below(10) as usize;
        b[p] = outra_base(b[p], r);
        b[q] = outra_base(b[q], r);
        add(&mut v, a, format!("D{i}"));
        add(&mut v, String::from_utf8(b).unwrap(), format!("D{i}"));
    }
    for i in 0..5 {
        let a = r.dna(20);
        add(&mut v, a[..19].to_string(), format!("P{i}"));
        add(&mut v, a, format!("P{i}"));
    }
    v
}

/// A read of `len` bases (or shorter), of a random kind.
fn leitura(r: &mut Rng, seqs: &[String], offsets: &[usize], len: usize, inverte: bool) -> String {
    let g = &seqs[r.below(seqs.len() as u64) as usize];
    let mut guia: Vec<u8> = g.clone().into_bytes();
    let tipo = r.below(100);
    match tipo {
        0..=44 => {}
        45..=59 => {
            let p = r.below(guia.len() as u64) as usize;
            guia[p] = outra_base(guia[p], r);
        }
        60..=64 => {
            let p = r.below(guia.len() as u64) as usize;
            guia[p] = b'N';
        }
        65..=69 => {
            for _ in 0..2 {
                let p = r.below(guia.len() as u64) as usize;
                guia[p] = outra_base(guia[p], r);
            }
        }
        70..=74 => guia = guia.iter().map(|b| b.to_ascii_lowercase()).collect(),
        _ => {}
    }
    let o = if tipo >= 90 { r.below(len as u64) as usize } else { offsets[r.below(offsets.len() as u64) as usize] };
    let mut s: Vec<u8> = r.dna(len).into_bytes();
    if (75..=84).contains(&tipo) {
        // random read
    } else {
        for (i, &b) in guia.iter().enumerate() {
            if o + i < s.len() {
                s[o + i] = b;
            }
        }
    }
    if (85..=89).contains(&tipo) {
        s.truncate(o + guia.len() / 2);
    }
    let s = String::from_utf8(s).unwrap();
    if inverte {
        rc(&s)
    } else {
        s
    }
}

fn grava_fastq(caminho: &str, reads: &[String]) {
    let mut w = externo::Escritor::cria(caminho).unwrap();
    for (i, s) in reads.iter().enumerate() {
        write!(w, "@r{i} x\n{s}\n+\n{}\n", "I".repeat(s.len())).unwrap();
    }
    w.termina().unwrap();
}

fn grava_gzip_comum(caminho: &str, reads: &[String]) {
    let f = std::fs::File::create(caminho).unwrap();
    let mut w = flate2::write::GzEncoder::new(f, flate2::Compression::fast());
    for (i, s) in reads.iter().enumerate() {
        write!(w, "@r{i}\r\n{s}\r\n+\r\n{}\r\n", "I".repeat(s.len())).unwrap();
    }
    w.finish().unwrap();
}

fn gini_ref(c: &[u64]) -> f64 {
    let mut x: Vec<f64> = c.iter().map(|&v| (v as f64 + 1.0).ln()).collect();
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = x.len() as f64;
    let gs: f64 = x.iter().enumerate().map(|(i, v)| (i as f64 + 1.0) * v).sum();
    let ys: f64 = x.iter().sum();
    1.0 - 2.0 * (n - gs / if ys == 0.0 { 1.0 } else { ys }) / (n - 1.0)
}

fn testa_contagem(t: &mut T, r: &mut Rng) {
    let lib = biblioteca_teste(r);
    let texto: String = lib.iter().map(|(i, s, g)| format!("{i}\t{s}\t{g}\n")).collect();
    std::fs::write(t.p("biblio.txt"), &texto).unwrap();
    let seqs: Vec<String> = lib.iter().map(|x| x.1.clone()).collect();
    let seqs_rc: Vec<String> = seqs.iter().map(|s| rc(s)).collect();

    // (offsets in the reads, offsets given to the counter, mismatches, reverse complement)
    let casos: [(&[usize], &[usize], u8, bool); 6] = [
        (&[3], &[3], 0, false),
        (&[3], &[3], 1, false),
        (&[0, 1, 2], &[0, 1, 2], 0, false),
        (&[0, 1, 2], &[2, 0, 1], 1, false),
        (&[5], &[5], 1, true),
        (&[4], &[3], 1, false),
    ];
    for (k, (no_read, dado, mm, inv)) in casos.iter().enumerate() {
        let reads: Vec<String> = (0..6000).map(|_| leitura(r, &seqs, no_read, 32, *inv)).collect();
        let arq = t.p(&format!("c{k}.fastq"));
        grava_fastq(&arq, &reads);
        let opts = CountOptions { offsets: Some(dado.to_vec()), reverse_complement: *inv, mismatches: *mm };
        let res = crispr_count(&t.p("biblio.txt"), &[vec![arq.clone()]], &[], None, &opts).unwrap();
        // brute force: the guides as the reads carry them
        let alvo: Vec<String> = if *inv { seqs_rc.clone() } else { seqs.clone() };
        let mut cont = vec![0u64; seqs.len()];
        let (mut ex, mut tr, mut am) = (0u64, 0u64, 0u64);
        for s in &reads {
            match bruto(&alvo, s, dado, *mm == 1) {
                Ref::Exato(g) => {
                    cont[g] += 1;
                    ex += 1;
                }
                Ref::Troca(g) => {
                    cont[g] += 1;
                    tr += 1;
                }
                Ref::Ambiguo => am += 1,
                Ref::Nada => {}
            }
        }
        let s = &res.samples[0];
        let nome = format!("reads at offsets {no_read:?}, counted at {dado:?}, {mm} mismatch(es){}", if *inv { ", reverse complement" } else { "" });
        t.eq(&res.counts[0], &cont, &format!("{nome}: counts per guide = brute force"));
        t.eq(&(s.reads, s.mapped_exact, s.mapped_mismatch, s.ambiguous), &(6000, ex, tr, am), &format!("{nome}: exact / one mismatch / ambiguous = brute force"));
        t.eq(&s.zero_guides, &cont.iter().filter(|&&c| c == 0).count(), &format!("{nome}: zero-count guides"));
        t.perto(s.gini, gini_ref(&cont), 1e-12, &format!("{nome}: Gini index"));
        if *mm == 1 && k == 1 {
            t.check(am > 0 && tr > 0, "the synthetic reads include one-mismatch reads and reads ambiguous between colliding guides");
        }
    }

    // collisions worked by hand: two guides one mismatch apart
    std::fs::write(t.p("colide.tsv"), "a\tAAAAAAAAAA\tX\nb\tAAAAAAAAAC\tX\nc\tCCCCCCCCCC\tY\nd\tGGGGGGGGGG\tZ\ne\tGGGGGTTGGG\tZ\n").unwrap();
    let reads: Vec<String> = [
        "AAAAAAAAAA", // exact a
        "AAAAAAAAAC", // exact b
        "AAAAAAAAAG", // one from a and from b: ambiguous
        "AAAAATAAAA", // one from a only
        "CCCCCNCCCC", // N as the mismatch: c
        "CCCCCNCCCN", // two N: nothing
        "GGGGGTGGGG", // one from d and one from e: ambiguous
        "TTTTTTTTTT", // nothing
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    grava_fastq(&t.p("colide.fastq"), &reads);
    let o1 = CountOptions { offsets: Some(vec![0]), reverse_complement: false, mismatches: 1 };
    let res = crispr_count(&t.p("colide.tsv"), &[vec![t.p("colide.fastq")]], &["S".to_string()], None, &o1).unwrap();
    t.eq(&res.counts[0], &vec![2, 1, 1, 0, 0], "hand-made collisions: a 2 (exact + one mismatch), b 1, c 1 (N), d and e 0 (ambiguous)");
    t.eq(&(res.samples[0].mapped_exact, res.samples[0].mapped_mismatch, res.samples[0].ambiguous), &(2, 2, 2), "hand-made collisions: 2 exact, 2 one mismatch, 2 ambiguous");
    let o0 = CountOptions { offsets: Some(vec![0]), ..Default::default() };
    let res = crispr_count(&t.p("colide.tsv"), &[vec![t.p("colide.fastq")]], &[], None, &o0).unwrap();
    t.eq(&res.counts[0], &vec![1, 1, 0, 0, 0], "exact only: just the two exact reads");

    // compressed input, replicates summed, labels, the table written
    let reads: Vec<String> = (0..3000).map(|_| leitura(r, &seqs, &[6], 30, false)).collect();
    grava_fastq(&t.p("rep1.fastq.gz"), &reads[..1500]);
    grava_gzip_comum(&t.p("rep2.fq.gz"), &reads[1500..]);
    grava_fastq(&t.p("todos.fastq"), &reads);
    let o = CountOptions { offsets: Some(vec![6]), ..Default::default() };
    let a = crispr_count(&t.p("biblio.txt"), &[vec![t.p("rep1.fastq.gz"), t.p("rep2.fq.gz")], vec![t.p("todos.fastq")]], &["reps".to_string(), "um".to_string()], Some(&t.p("tabela.txt")), &o).unwrap();
    t.eq(&a.counts[0], &a.counts[1], "BGZF + plain gzip with CRLF, as two replicates = the same reads in one text file");
    t.eq(&a.samples[0].files.len(), &2, "two files in the replicate sample");
    let tab = read_count_table(&t.p("tabela.txt")).unwrap();
    t.eq(&tab.samples, &vec!["reps".to_string(), "um".to_string()], "count table header: sgRNA, Gene, the labels");
    t.eq(&tab.ids, &lib.iter().map(|x| x.0.clone()).collect::<Vec<_>>(), "count table in library order");
    t.check(tab.counts[1].iter().zip(&a.counts[1]).all(|(x, &y)| *x == y as f64), "count table holds the counts");
    t.check(a.countsummary().starts_with("File\tLabel\tReads\tMapped\tPercentage") && a.countsummary().lines().count() == 3, "MAGeCK's countsummary table");
    t.check(a.json().contains("\"labels\":[\"reps\",\"um\"]") && a.resumo().contains("reps:"), "count JSON and summary");

    // offset found from the reads, as MAGeCK's AUTO
    let reads: Vec<String> = (0..4000).map(|_| leitura(r, &seqs, &[7], 36, false)).collect();
    grava_fastq(&t.p("auto.fastq"), &reads);
    let res = crispr_count(&t.p("biblio.txt"), &[vec![t.p("auto.fastq")]], &[], None, &CountOptions::default()).unwrap();
    t.eq(&res.samples[0].files[0].offsets, &vec![7], "offset 7 found from the reads");
    t.check(res.samples[0].files[0].detected && res.samples[0].files[0].offset_probe[0].0 == 7, "the probe is reported");
    let dado = crispr_count(&t.p("biblio.txt"), &[vec![t.p("auto.fastq")]], &[], None, &CountOptions { offsets: Some(vec![7]), ..Default::default() }).unwrap();
    t.eq(&res.counts, &dado.counts, "found offset counts the same as the offset given");
    // staggered: offsets 2, 3 and 4 in equal parts are all kept
    let reads: Vec<String> = (0..6000).map(|_| leitura(r, &seqs, &[2, 3, 4], 36, false)).collect();
    grava_fastq(&t.p("escada.fastq"), &reads);
    let res = crispr_count(&t.p("biblio.txt"), &[vec![t.p("escada.fastq")]], &[], None, &CountOptions::default()).unwrap();
    let mut o = res.samples[0].files[0].offsets.clone();
    o.sort();
    t.eq(&o, &vec![2, 3, 4], "staggered offsets 2, 3, 4 found");
    // MAGeCK's selection rule, worked by hand
    t.eq(&crispr::internals::auto_offsets(&[(0, 900), (1, 50), (2, 50)], 1000), &vec![0], "offset rule: 90% at one offset takes only it");
    t.eq(&crispr::internals::auto_offsets(&[(3, 400), (4, 350), (5, 250)], 1000), &vec![3, 4, 5], "offset rule: 40/35/25% takes all three");
    t.eq(&crispr::internals::auto_offsets(&[(3, 700), (4, 100), (9, 5)], 805), &vec![3], "offset rule: a third of the previous is needed");
    t.eq(&crispr::internals::auto_offsets(&[(2, 500), (3, 400), (7, 100)], 1000), &vec![2, 3], "offset rule: stops once 75% is covered");
    // random reads: no offset can be found
    let lixo: Vec<String> = (0..500).map(|_| r.dna(30)).collect();
    grava_fastq(&t.p("lixo.fastq"), &lixo);
    t.erro(crispr_count(&t.p("biblio.txt"), &[vec![t.p("lixo.fastq")]], &[], None, &CountOptions::default()), &["lixo.fastq", "cannot find where the guide is", "reverse_complement"], "no guide in the reads: the offset cannot be found");

    // malformed input
    std::fs::write(t.p("fasta.fq"), ">r1\nACGT\n").unwrap();
    t.erro(crispr_count(&t.p("biblio.txt"), &[vec![t.p("fasta.fq")]], &[], None, &o0), &["fasta.fq", "line 1", "FASTA"], "FASTA given as FASTQ");
    std::fs::write(t.p("corte.fq"), "@r1\nACGTACGT\n+\nIIII\n").unwrap();
    t.erro(crispr_count(&t.p("biblio.txt"), &[vec![t.p("corte.fq")]], &[], None, &o0), &["corte.fq", "line 4", "8 bases but 4 quality"], "quality shorter than the sequence");
    std::fs::write(t.p("meio.fq"), "@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nACGT\n").unwrap();
    t.erro(crispr_count(&t.p("biblio.txt"), &[vec![t.p("meio.fq")]], &[], None, &o0), &["meio.fq", "ends at line 6", "record 2"], "file cut in the middle of a record");
    t.erro(crispr_count(&t.p("biblio.txt"), &[vec![t.p("todos.fastq")]], &["a".to_string(), "b".to_string()], None, &o0), &["2 label(s) for 1 sample"], "more labels than samples");
    t.erro(crispr_count(&t.p("biblio.txt"), &[vec![t.p("todos.fastq")], vec![t.p("todos.fastq")]], &["a".to_string(), "a".to_string()], None, &o0), &["used twice"], "repeated label");
    t.erro(crispr_count(&t.p("biblio.txt"), &[vec![t.p("todos.fastq")]], &[], Some(&t.p("todos.fastq")), &o0), &["destroy"], "output = input");
    t.erro(crispr_count(&t.p("biblio.txt"), &[vec![t.p("todos.fastq")]], &[], None, &CountOptions { mismatches: 2, ..Default::default() }), &["mismatches=2"], "two mismatches are not offered");
    t.erro(crispr_count(&t.p("biblio.txt"), &[vec![]], &[], None, &o0), &["no FASTQ"], "sample without files");
    t.erro(crispr_count(&t.p("biblio.txt"), &[vec![t.p("nada.fastq")]], &[], None, &o0), &["nada.fastq"], "missing FASTQ");
    t.check(!Path::new(&t.p("x_saida.txt")).exists(), "no output left behind");
}

// ---------------------------------------------------------------- count table

fn testa_tabela(t: &mut T) {
    std::fs::write(t.p("t1.txt"), "sgRNA\tGene\tA\tB\ns1\tG1\t10\t20\ns2\tG1\t0\t5.5\ns1\tG2\t1\t1\n").unwrap();
    let c = read_count_table(&t.p("t1.txt")).unwrap();
    t.eq(&(c.samples.clone(), c.ids.clone(), c.duplicate_ids), &(vec!["A".to_string(), "B".to_string()], vec!["s1".to_string(), "s2".to_string()], 1), "TSV count table; repeated id skipped as in MAGeCK");
    t.eq(&c.counts, &vec![vec![10.0, 0.0], vec![20.0, 5.5]], "counts per sample, fractional counts accepted");
    std::fs::write(t.p("t2.csv"), "s1,G1,10,20\ns2,G2,3,4\n").unwrap();
    let c = read_count_table(&t.p("t2.csv")).unwrap();
    t.eq(&(c.header, c.samples.clone()), &(false, vec!["sample1".to_string(), "sample2".to_string()]), "CSV without header: samples sample1, sample2");
    std::fs::write(t.p("t3.txt"), "sgRNA\tGene\tA\ns1\tG1\tdez\n").unwrap();
    t.erro(read_count_table(&t.p("t3.txt")), &["t3.txt", "line 2", "`dez`", "column `A`"], "non-numeric count");
    std::fs::write(t.p("t4.txt"), "sgRNA\tGene\tA\ns1\tG1\t-3\n").unwrap();
    t.erro(read_count_table(&t.p("t4.txt")), &["line 2", ">= 0"], "negative count");
    std::fs::write(t.p("t5.txt"), "sgRNA\tGene\tA\tB\ns1\tG1\t3\n").unwrap();
    t.erro(read_count_table(&t.p("t5.txt")), &["line 2", "3 columns", "first line has 4"], "ragged line");
    std::fs::write(t.p("t6.txt"), "sgRNA\tGene\tA\ns1\t\t3\n").unwrap();
    t.erro(read_count_table(&t.p("t6.txt")), &["line 2", "no gene"], "empty gene");
    std::fs::write(t.p("t7.txt"), "sgRNA\tGene\tA\tA\ns1\tG\t3\t4\n").unwrap();
    t.erro(read_count_table(&t.p("t7.txt")), &["`A` names two columns"], "repeated sample label");
    std::fs::write(t.p("t8.txt"), "sgRNA\tGene\n").unwrap();
    t.erro(read_count_table(&t.p("t8.txt")), &["at least 3"], "table without count columns");
}

// ---------------------------------------------------------------- simulated screens

struct Tela {
    tabela: CountTable,
    essenciais: HashSet<String>,
    enriquecidos: HashSet<String>,
}

/// A screen of `ng` genes of 4 guides; NB counts (dispersion `phi`) with
/// lognormal guide abundances; `ess` depleted and `enr` enriched genes with
/// three of their four guides at +/- `efeito` (log2).
#[allow(clippy::too_many_arguments)]
fn simula(r: &mut Rng, ng: usize, ess: usize, enr: usize, efeito: f64, phi: f64, nc: usize, nt: usize) -> Tela {
    let mut ids = Vec::new();
    let mut genes = Vec::new();
    let mut mu = Vec::new();
    let mut ef = Vec::new();
    let mut essenciais = HashSet::new();
    let mut enriquecidos = HashSet::new();
    for g in 0..ng {
        let nome = format!("GENE{g:05}");
        let e = if g < ess {
            essenciais.insert(nome.clone());
            -efeito
        } else if g < ess + enr {
            enriquecidos.insert(nome.clone());
            efeito
        } else {
            0.0
        };
        for k in 0..4 {
            ids.push(format!("{nome}_{k}"));
            genes.push(nome.clone());
            mu.push((6.0 + 0.6 * r.normal()).exp());
            ef.push(if k == 3 { 0.0 } else { e });
        }
    }
    let mut samples = Vec::new();
    let mut counts = Vec::new();
    for s in 0..nc + nt {
        let trat = s >= nc;
        samples.push(if trat { format!("T{}", s - nc + 1) } else { format!("C{}", s + 1) });
        let prof = 0.7 + 0.7 * r.u();
        counts.push(mu.iter().zip(&ef).map(|(&m, &e)| r.nb(prof * m * if trat { 2f64.powf(e) } else { 1.0 }, phi) as f64).collect());
    }
    Tela { tabela: CountTable { path: "simulated".to_string(), samples, ids, genes, counts, duplicate_ids: 0, header: true }, essenciais, enriquecidos }
}

fn rotulos(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| s.to_string()).collect()
}

fn testa_tela(t: &mut T, r: &mut Rng) {
    let o = TestOptions::default();
    // planted essential and enriched genes
    let tela = simula(r, 600, 25, 15, 2.0, 0.05, 2, 2);
    let res = crispr_test_table(&tela.tabela, &rotulos(&["T1", "T2"]), &rotulos(&["C1", "C2"]), Some(&t.p("sim")), &o).unwrap();
    let neg: Vec<&GeneStat> = res.genes.iter().filter(|g| g.neg_fdr < 0.1).collect();
    let rec = res.genes.iter().filter(|g| tela.essenciais.contains(&g.gene) && g.neg_fdr < 0.1).count();
    let fp = neg.iter().filter(|g| !tela.essenciais.contains(&g.gene)).count();
    t.check(rec >= 23, &format!("planted depleted genes found at FDR < 0.1: {rec} of 25"));
    t.check(fp <= 2, &format!("negative hits at FDR < 0.1 that were not planted: {fp}"));
    let topo: Vec<String> = res.top(20, false).iter().map(|g| g.gene.clone()).collect();
    t.check(topo.iter().all(|g| tela.essenciais.contains(g)), "the 20 most depleted genes are all planted essential genes");
    let rec = res.genes.iter().filter(|g| tela.enriquecidos.contains(&g.gene) && g.pos_fdr < 0.1).count();
    let fp = res.genes.iter().filter(|g| !tela.enriquecidos.contains(&g.gene) && g.pos_fdr < 0.1).count();
    t.check(rec >= 13, &format!("planted enriched genes found at FDR < 0.1: {rec} of 15"));
    t.check(fp <= 2, &format!("positive hits at FDR < 0.1 that were not planted: {fp}"));
    println!(
        "Simulated screen (600 genes x 4 guides, 2 vs 2, NB dispersion 0.05, effect 2 log2 on 3 of 4 guides): depleted {} of 25 found at FDR < 0.1 with {} false; enriched {} of 15 with {} false",
        res.genes.iter().filter(|g| tela.essenciais.contains(&g.gene) && g.neg_fdr < 0.1).count(),
        neg.iter().filter(|g| !tela.essenciais.contains(&g.gene)).count(),
        rec,
        fp
    );
    t.check(res.genes.iter().filter(|g| tela.essenciais.contains(&g.gene)).all(|g| g.neg_lfc < -1.0), "gene LFC of the depleted genes is strongly negative");

    // internal consistency of the result
    let n = res.genes.len();
    t.eq(&n, &600, "every gene in the gene table");
    let mut ranks: Vec<usize> = res.genes.iter().map(|g| g.neg_rank).collect();
    t.check(ranks.windows(2).all(|w| w[0] < w[1]) && ranks[0] == 1 && ranks[n - 1] == n, "gene table ordered by negative rank 1..n");
    ranks = res.genes.iter().map(|g| g.pos_rank).collect();
    ranks.sort();
    t.check(ranks == (1..=n).collect::<Vec<_>>(), "positive ranks are 1..n");
    t.check(res.genes.windows(2).all(|w| w[0].neg_score <= w[1].neg_score), "negative rank follows the RRA score");
    let mut por_score: Vec<&GeneStat> = res.genes.iter().collect();
    por_score.sort_by(|a, b| a.neg_p.partial_cmp(&b.neg_p).unwrap());
    t.check(por_score.windows(2).all(|w| w[0].neg_fdr <= w[1].neg_fdr + 1e-15), "FDR grows with the p-value");
    t.check(res.genes.iter().all(|g| g.neg_p > 0.0 && g.neg_p <= 1.0 && g.neg_good <= g.guides), "p-values in (0, 1]; good guides <= guides");
    let lfc_de: HashMap<&str, Vec<f64>> = res.guides.iter().fold(HashMap::new(), |mut m, g| {
        m.entry(g.gene.as_str()).or_default().push(g.lfc);
        m
    });
    t.check(res.genes.iter().all(|g| (g.neg_lfc - median(&lfc_de[g.gene.as_str()])).abs() < 1e-12 && g.neg_lfc == g.pos_lfc), "gene LFC = median of its guides' LFC");
    let g = &res.guides;
    t.check(g.windows(2).all(|w| w[0].score >= w[1].score), "guide table ordered by |score|");
    t.check(g.iter().all(|x| x.p_twosided == 2.0 * x.p_low.min(x.p_high)), "two-sided p = 2 min(p.low, p.high)");
    t.check(g.iter().all(|x| (x.lfc - ((x.treatment_mean + 1.0).log2() - (x.control_mean + 1.0).log2())).abs() < 1e-12), "guide LFC = log2(treat + 1) - log2(control + 1)");
    t.check(g.iter().all(|x| x.high_in_treatment == (x.treatment_mean > x.control_mean)), "high_in_treatment");
    let p2: Vec<f64> = g.iter().map(|x| x.p_twosided).collect();
    t.check(bh_fdr(&p2).iter().zip(g).all(|(a, x)| *a == x.fdr), "guide FDR = BH of the two-sided p-values");
    t.check(g.iter().all(|x| x.adj_var >= x.control_mean && (x.adj_var - (x.control_mean.powf(res.model_k) * 2f64.powf(res.model_b) + x.control_mean)).abs() <= 1e-9 * x.adj_var), "adjusted variance follows the fitted model");
    t.check(res.alpha_neg >= 0.05 && res.alpha_neg <= 0.5, "negative alpha within [0.05, 0.5]");
    t.eq(&res.random_genes, &(600 * (100_000 / 600 + 2)), "random genes per size: RRA's rule (100,000 / genes + 1 passes, + 1)");
    // files
    let gs = std::fs::read_to_string(t.p("sim.gene_summary.txt")).unwrap();
    let ss = std::fs::read_to_string(t.p("sim.sgrna_summary.txt")).unwrap();
    t.check(gs.starts_with("id\tnum\tneg|score\tneg|p-value\tneg|fdr\tneg|rank\tneg|goodsgrna\tneg|lfc\tpos|score") && gs.lines().count() == 601, "gene_summary.txt: MAGeCK's header, one line per gene");
    t.check(ss.starts_with("sgrna\tGene\tcontrol_count\ttreatment_count\tcontrol_mean\ttreat_mean\tLFC") && ss.lines().count() == 2401, "sgrna_summary.txt: MAGeCK's header, one line per guide");
    t.eq(&res.outputs.len(), &2, "both outputs listed");
    let j = res.json();
    t.check(j.contains("\"rank_plot\":{\"neg\":{\"gene\":[") && j.contains("\"top_negative\":[") && j.contains("\"size_factors\":[[\"C1\","), "JSON: rank plot, top genes, size factors");
    t.check(res.resumo().contains("negative (depleted) selection") && res.resumo().contains("size factors C1"), "text summary");
    let jc = res.json_compact();
    t.check(jc.contains("\"guides\":null") && jc.contains("\"genes\":{\"id\":[") && jc.len() * 3 < j.len(), "compact JSON: everything but the guide table");
    // determinism
    let de_novo = crispr_test_table(&tela.tabela, &rotulos(&["T1", "T2"]), &rotulos(&["C1", "C2"]), None, &o).unwrap();
    t.check(de_novo.genes == res.genes && de_novo.guides == res.guides, "the same input gives the same result, bit for bit");
    let outra = crispr_test_table(&tela.tabela, &rotulos(&["T1", "T2"]), &rotulos(&["C1", "C2"]), None, &TestOptions { seed: 77, ..o.clone() }).unwrap();
    t.check(outra.genes.iter().zip(&res.genes).any(|(a, b)| a.neg_p != b.neg_p), "another seed gives other permutation p-values");
    // control by default = the rest; columns by number
    let resto = crispr_test_table(&tela.tabela, &rotulos(&["2", "3"]), &[], None, &o).unwrap();
    t.check(resto.control == rotulos(&["C1", "C2"]) && resto.genes == res.genes, "control omitted = the other samples; samples by column number");

    // calibration: no effect at all
    let nulo = simula(r, 2000, 0, 0, 0.0, 0.05, 2, 2);
    let res = crispr_test_table(&nulo.tabela, &rotulos(&["T1", "T2"]), &rotulos(&["C1", "C2"]), None, &o).unwrap();
    let n = res.genes.len() as f64;
    for (lado, ps) in [("negative", res.genes.iter().map(|g| g.neg_p).collect::<Vec<f64>>()), ("positive", res.genes.iter().map(|g| g.pos_p).collect())] {
        let f05 = ps.iter().filter(|&&p| p < 0.05).count() as f64 / n;
        let f20 = ps.iter().filter(|&&p| p < 0.2).count() as f64 / n;
        println!("Null screen (2,000 genes, no effect), {lado} selection: {:.1}% of the genes at p < 0.05, {:.1}% at p < 0.2", 100.0 * f05, 100.0 * f20);
        t.check((0.03..=0.07).contains(&f05), &format!("null screen, {lado} selection: {:.1}% of the genes at p < 0.05 (expected ~5%)", 100.0 * f05));
        t.check((0.16..=0.24).contains(&f20), &format!("null screen, {lado} selection: {:.1}% of the genes at p < 0.2 (expected ~20%)", 100.0 * f20));
    }
    let falsos = res.genes.iter().filter(|g| g.neg_fdr < 0.1 || g.pos_fdr < 0.1).count();
    println!("Null screen: {falsos} gene(s) at FDR < 0.1 in either direction");
    t.check(falsos <= 2, &format!("null screen: {falsos} gene(s) at FDR < 0.1"));

    // options
    let tela = simula(r, 300, 15, 0, 2.0, 0.05, 3, 3);
    let mut tab = tela.tabela.clone();
    // a gene with no reads anywhere, and a gene of 4 non-targeting controls
    for s in 0..6 {
        for k in 0..4 {
            tab.counts[s][4 * 299 + k] = 0.0;
        }
    }
    for k in 0..4 {
        tab.genes[4 * 298 + k] = "NonTargeting".to_string();
    }
    let trat = rotulos(&["T1", "T2", "T3"]);
    let ctrl = rotulos(&["C1", "C2", "C3"]);
    let base = crispr_test_table(&tab, &trat, &ctrl, None, &o).unwrap();
    t.eq(&(base.guides_used, base.genes.len()), &(1196, 299), "remove_zero=both: the 4 guides without reads leave, and their gene");
    let todos = crispr_test_table(&tab, &trat, &ctrl, None, &TestOptions { remove_zero: RemoveZero::None, ..o.clone() }).unwrap();
    t.eq(&(todos.guides_used, todos.genes.len()), &(1200, 300), "remove_zero=none keeps every guide");
    let ctl: Vec<String> = (0..4).map(|k| tab.ids[4 * 298 + k].clone()).collect();
    let com_ctrl = crispr_test_table(&tab, &trat, &ctrl, None, &TestOptions { control_sgrnas: ctl.clone(), ..o.clone() }).unwrap();
    t.check(com_ctrl.control_guides == 4 && com_ctrl.genes.iter().all(|g| g.gene != "NonTargeting") && com_ctrl.genes.len() == 298, "control sgRNAs give the null; their gene leaves the gene table");
    let por_gene = crispr_test_table(&tab, &trat, &ctrl, None, &TestOptions { control_genes: vec!["NonTargeting".to_string()], ..o.clone() }).unwrap();
    t.check(por_gene.genes == com_ctrl.genes, "control_genes = the same guides as control_sgrnas");
    let norm_ctrl = crispr_test_table(&tab, &trat, &ctrl, None, &TestOptions { control_sgrnas: ctl.clone(), norm_method: NormMethod::Control, ..o.clone() }).unwrap();
    let cols: Vec<Vec<f64>> = [0usize, 1, 2, 3, 4, 5].iter().map(|&s| (0..4).map(|k| tab.counts[s][4 * 298 + k]).collect()).collect();
    let (f, _, _) = crispr::internals::size_factors(&cols, NormMethod::Median).unwrap();
    t.check(norm_ctrl.norm_used == "control" && norm_ctrl.size_factors.iter().zip(&f).all(|(a, b)| (a.1 - b).abs() < 1e-12), "control normalization = median ratio over the control guides only");
    let tot = crispr_test_table(&tab, &trat, &ctrl, None, &TestOptions { norm_method: NormMethod::Total, ..o.clone() }).unwrap();
    t.eq(&tot.norm_used, &"total".to_string(), "total normalization");
    let pos = crispr_test_table(&tab, &trat, &ctrl, None, &TestOptions { sort_by_positive: true, ..o.clone() }).unwrap();
    t.check(pos.genes.windows(2).all(|w| w[0].pos_rank < w[1].pos_rank), "gene table sorted by positive selection");
    let sb = crispr_test_table(&tab, &trat, &ctrl, None, &TestOptions { gene_lfc: GeneLfc::SecondBest, ..o.clone() }).unwrap();
    let mut ok = true;
    for g in &sb.genes {
        let mut v: Vec<f64> = sb.guides.iter().filter(|x| x.gene == g.gene && x.used).map(|x| x.lfc).collect();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let neg = if v.len() > 1 { v[1] } else { 0.0 };
        let posv = if v.len() > 1 { v[v.len() - 2] } else { 0.0 };
        ok &= g.neg_lfc == neg && g.pos_lfc == posv;
    }
    t.check(ok, "gene_lfc=secondbest: the second most depleted / enriched guide");
    let mean = crispr_test_table(&tab, &trat, &ctrl, None, &TestOptions { gene_lfc: GeneLfc::Mean, ..o.clone() }).unwrap();
    t.check(mean.genes.iter().all(|g| {
        let v: Vec<f64> = mean.guides.iter().filter(|x| x.gene == g.gene && x.used).map(|x| x.lfc).collect();
        (g.neg_lfc - v.iter().sum::<f64>() / v.len() as f64).abs() < 1e-12
    }), "gene_lfc=mean");
    let perm = crispr_test_table(&tab, &trat, &ctrl, None, &TestOptions { permutations: Some(20), ..o.clone() }).unwrap();
    t.eq(&perm.random_genes, &(299 * 21), "permutations=20: 21 passes over the genes");
    let skip = crispr_test_table(&tab, &trat, &ctrl, None, &TestOptions { skip_genes: vec!["GENE00000".to_string()], ..o.clone() }).unwrap();
    t.check(skip.genes.iter().find(|g| g.gene == "GENE00000").map(|g| g.neg_p == 1.0).unwrap_or(false), "a skipped gene gets p = 1");
    // one control sample: the variance comes from all samples
    let um = crispr_test_table(&tab, &trat, &rotulos(&["C1"]), None, &o).unwrap();
    t.check(um.genes.iter().filter(|g| g.neg_fdr < 0.1).count() >= 10, "one control sample still finds the depleted genes");

    // errors
    t.erro(crispr_test_table(&tab, &rotulos(&["T9"]), &ctrl, None, &o), &["no sample `T9`", "C1, C2, C3, T1, T2, T3"], "unknown sample");
    t.erro(crispr_test_table(&tab, &trat, &rotulos(&["T1"]), None, &o), &["`T1` is both treatment and control"], "sample in both groups");
    t.erro(crispr_test_table(&tab, &[], &ctrl, None, &o), &["No treatment"], "no treatment");
    t.erro(crispr_test_table(&tab, &trat, &ctrl, None, &TestOptions { control_sgrnas: vec!["nada".to_string()], ..o.clone() }), &["Only 0 of the 1 control sgRNA"], "control sgRNAs not in the table");
    let mut mistura = ctl.clone();
    mistura.push(tab.ids[0].clone());
    t.erro(crispr_test_table(&tab, &trat, &ctrl, None, &TestOptions { control_sgrnas: mistura, ..o.clone() }), &["GENE00000", "both control and non-control"], "gene with control and non-control guides");
    t.erro(crispr_test_table(&tab, &trat, &ctrl, None, &TestOptions { norm_method: NormMethod::Control, ..o.clone() }), &["needs control_sgrnas"], "control normalization without controls");
    let mut zero = tab.clone();
    for s in 0..3 {
        zero.counts[s].iter_mut().for_each(|x| *x = 0.0);
    }
    t.erro(crispr_test_table(&zero, &trat, &ctrl, None, &o), &["no reads"], "control samples without reads");
    t.erro(crispr_test_table(&tab, &trat, &ctrl, None, &TestOptions { gene_test_fdr_threshold: 0.0, ..o.clone() }), &["gene_test_fdr_threshold"], "threshold 0");
    t.erro(crispr_test_table(&tab, &trat, &ctrl, None, &TestOptions { seed: 0, ..o.clone() }), &["seed=0"], "seed 0");
    t.erro(crispr_test(&t.p("nao_ha.txt"), &trat, &ctrl, None, &o), &["nao_ha.txt"], "missing count table");
    t.check(NormMethod::from_text("Median") == Ok(NormMethod::Median) && NormMethod::from_text("x").is_err(), "normalization names");
    t.check(RemoveZero::from_text("any") == Ok(RemoveZero::Any) && GeneLfc::from_text("alphamedian") == Ok(GeneLfc::AlphaMedian), "option names");

    // count -> table file -> test
    let lib: String = (0..120).map(|i| format!("g{i}\t{}\tGENE{}\n", { let mut q = Rng(i as u64 + 99); q.dna(20) }, i / 4)).collect();
    std::fs::write(t.p("rt.tsv"), &lib).unwrap();
    let seqs: Vec<String> = lib.lines().map(|l| l.split('\t').nth(1).unwrap().to_string()).collect();
    for (k, nome) in ["c1", "c2", "t1", "t2"].iter().enumerate() {
        let reads: Vec<String> = (0..20000)
            .map(|_| {
                let mut g = r.below(120) as usize;
                // genes 0..4 depleted in treatment
                if k >= 2 && g < 20 && r.chance(80) {
                    g = 20 + r.below(100) as usize;
                }
                format!("{}{}", seqs[g], r.dna(5))
            })
            .collect();
        grava_fastq(&t.p(&format!("{nome}.fq")), &reads);
    }
    let c = crispr_count(&t.p("rt.tsv"), &[vec![t.p("c1.fq")], vec![t.p("c2.fq")], vec![t.p("t1.fq")], vec![t.p("t2.fq")]], &rotulos(&["c1", "c2", "t1", "t2"]), Some(&t.p("rt.count.txt")), &CountOptions::default()).unwrap();
    t.check(c.samples.iter().all(|s| s.mapped == 20000 && s.files[0].offsets == vec![0]), "count of the round trip: every read mapped at offset 0");
    let res = crispr_test(&t.p("rt.count.txt"), &rotulos(&["t1", "t2"]), &rotulos(&["c1", "c2"]), None, &o).unwrap();
    let topo: HashSet<String> = res.top(5, false).iter().map(|g| g.gene.clone()).collect();
    t.eq(&topo, &(0..5).map(|i| format!("GENE{i}")).collect::<HashSet<_>>(), "count -> test: the 5 depleted genes on top");
}

// ---------------------------------------------------------------- real data

fn le_tsv(p: &str) -> Vec<Vec<String>> {
    std::fs::read_to_string(p).unwrap().lines().map(|l| l.split('\t').map(str::to_string).collect()).collect()
}

/// Spearman correlation (average ranks for ties).
fn spearman(a: &[f64], b: &[f64]) -> f64 {
    fn postos(v: &[f64]) -> Vec<f64> {
        let mut i: Vec<usize> = (0..v.len()).collect();
        i.sort_by(|&x, &y| v[x].partial_cmp(&v[y]).unwrap());
        let mut r = vec![0.0; v.len()];
        let mut k = 0;
        while k < i.len() {
            let mut j = k;
            while j + 1 < i.len() && v[i[j + 1]] == v[i[k]] {
                j += 1;
            }
            for m in k..=j {
                r[i[m]] = (k + j) as f64 / 2.0 + 1.0;
            }
            k = j + 1;
        }
        r
    }
    let (x, y) = (postos(a), postos(b));
    let n = x.len() as f64;
    let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
    let c: f64 = x.iter().zip(&y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sx: f64 = x.iter().map(|a| (a - mx) * (a - mx)).sum();
    let sy: f64 = y.iter().map(|b| (b - my) * (b - my)).sum();
    c / (sx * sy).sqrt()
}

/// Compares our gene and guide tables with MAGeCK's files.
fn compara_mageck(t: &mut T, nome: &str, res: &TestResult, dir: &Path, prefixo: &str) {
    let gp = dir.join(format!("{prefixo}.gene_summary.txt"));
    let sp = dir.join(format!("{prefixo}.sgrna_summary.txt"));
    if !gp.exists() || !sp.exists() {
        println!("  ({nome}: no MAGeCK output in {}, comparison skipped)", dir.display());
        return;
    }
    let mg = le_tsv(&gp.to_string_lossy());
    let ms = le_tsv(&sp.to_string_lossy());
    // guide table: same order, every column within the 5 digits MAGeCK prints
    let nossa: Vec<Vec<String>> = res.sgrna_summary().lines().map(|l| l.split('\t').map(str::to_string).collect()).collect();
    let mut iguais = 0;
    let mut pior = (0f64, String::new());
    let por_id: HashMap<&str, &Vec<String>> = ms[1..].iter().map(|l| (l[0].as_str(), l)).collect();
    for l in &nossa[1..] {
        let m = por_id[l[0].as_str()];
        if l == m {
            iguais += 1;
        }
        for c in 4..14 {
            let (a, b): (f64, f64) = (l[c].parse().unwrap(), m[c].parse().unwrap());
            let d = if b == 0.0 { a.abs() } else { ((a - b) / b).abs() };
            if d > pior.0 {
                pior = (d, format!("{} column {}: {} vs {}", l[0], ms[0][c], l[c], m[c]));
            }
        }
    }
    let ordem_igual = nossa[1..].iter().zip(&ms[1..]).filter(|(a, b)| a[0] == b[0]).count();
    println!(
        "  {nome}: guides identical to MAGeCK's line {iguais} of {} (same position {ordem_igual}); largest difference {:.1e} ({})",
        ms.len() - 1,
        pior.0,
        pior.1
    );
    t.check(pior.0 < 2e-4, &format!("{nome}: every number of the guide table within MAGeCK's printed precision"));
    // gene table
    let mgd: HashMap<&str, &Vec<String>> = mg[1..].iter().map(|l| (l[0].as_str(), l)).collect();
    t.eq(&res.genes.len(), &(mg.len() - 1), &format!("{nome}: same number of genes as MAGeCK"));
    let ns: Vec<String> = res.gene_summary().lines().map(str::to_string).collect();
    let linhas_iguais = ns[1..].iter().zip(&mg[1..]).filter(|(a, b)| a.split('\t').collect::<Vec<_>>() == b.iter().map(String::as_str).collect::<Vec<_>>()).count();
    let mut a = Vec::new();
    let mut b = Vec::new();
    let mut fa = Vec::new();
    let mut fb = Vec::new();
    let mut pa = Vec::new();
    let mut pb = Vec::new();
    for g in &res.genes {
        let m = mgd[g.gene.as_str()];
        a.push(g.neg_score);
        b.push(m[2].parse::<f64>().unwrap());
        pa.push(g.pos_score);
        pb.push(m[8].parse::<f64>().unwrap());
        fa.push(g.neg_fdr);
        fb.push(m[4].parse::<f64>().unwrap());
    }
    // Ranks up to ties: genes with exactly the same lo-value form a group, and
    // RRA orders a group with an unstable quicksort (we keep the order of
    // first appearance). Within each group the set of ranks must be the same.
    let ranks_empatados = |pos: bool| -> (bool, usize) {
        let mut ordem: Vec<&GeneStat> = res.genes.iter().collect();
        ordem.sort_by_key(|g| if pos { g.pos_rank } else { g.neg_rank });
        let mut ok = true;
        let mut empatados = 0;
        let mut i = 0;
        while i < ordem.len() {
            let lo = |g: &GeneStat| if pos { g.pos_score } else { g.neg_score };
            let mut j = i;
            while j + 1 < ordem.len() && lo(ordem[j + 1]) == lo(ordem[i]) {
                j += 1;
            }
            let mut nossos: Vec<usize> = ordem[i..=j].iter().map(|g| if pos { g.pos_rank } else { g.neg_rank }).collect();
            let mut deles: Vec<usize> = ordem[i..=j].iter().map(|g| mgd[g.gene.as_str()][if pos { 11 } else { 5 }].parse().unwrap()).collect();
            nossos.sort();
            deles.sort();
            ok &= nossos == deles;
            if j > i {
                empatados += j - i + 1;
            }
            i = j + 1;
        }
        (ok, empatados)
    };
    let (rn, en) = ranks_empatados(false);
    let (rp, ep) = ranks_empatados(true);
    // every column but the two ranks
    let sem_rank = |l: &[&str]| -> Vec<String> { l.iter().enumerate().filter(|(i, _)| *i != 5 && *i != 11).map(|(_, x)| x.to_string()).collect() };
    let por_gene_nosso: HashMap<String, Vec<String>> = ns[1..].iter().map(|l| {
        let c: Vec<&str> = l.split('\t').collect();
        (c[0].to_string(), sem_rank(&c))
    }).collect();
    let resto_igual = mg[1..].iter().filter(|l| {
        let c: Vec<&str> = l.iter().map(String::as_str).collect();
        por_gene_nosso.get(c[0]) == Some(&sem_rank(&c))
    }).count();
    println!("  {nome}: gene lines identical apart from the ranks {resto_igual} of {}; ranks equal up to the order of tied genes: {} (neg, {en} genes in ties), {} (pos, {ep} in ties)", mg.len() - 1, rn, rp);
    t.check(resto_igual == mg.len() - 1, &format!("{nome}: every gene line identical to MAGeCK's apart from the rank columns"));
    t.check(rn && rp, &format!("{nome}: gene ranks equal to MAGeCK's up to the order of genes with the same RRA score"));
    let rho_n = spearman(&a, &b);
    let rho_p = spearman(&pa, &pb);
    let topo = |pos: bool| -> Vec<String> { res.top(20, pos).iter().map(|g| g.gene.clone()).collect() };
    let topo_m = |col: usize| -> Vec<String> {
        let mut v: Vec<&Vec<String>> = mg[1..].iter().filter(|l| !res.skipped_genes.contains(&l[0])).collect();
        v.sort_by_key(|l| l[col].parse::<usize>().unwrap());
        v.iter().take(20).map(|l| l[0].clone()).collect()
    };
    let (tn, tmn) = (topo(false), topo_m(5));
    let (tp, tmp) = (topo(true), topo_m(11));
    let comum_n = tn.iter().filter(|g| tmn.contains(g)).count();
    let comum_p = tp.iter().filter(|g| tmp.contains(g)).count();
    let fdr_concorda = |lim: f64| fa.iter().zip(&fb).filter(|(x, y)| (**x < lim) == (**y < lim)).count();
    let dif_fdr = fa.iter().zip(&fb).map(|(x, y)| (x - y).abs()).fold(0f64, f64::max);
    println!(
        "  {nome}: gene lines identical {linhas_iguais} of {}; Spearman of RRA scores {rho_n:.6} (neg) {rho_p:.6} (pos); top-20 in common {comum_n} (neg) {comum_p} (pos); FDR < 0.05 / 0.25 agree for {} / {} genes; largest FDR difference {dif_fdr:.2e}",
        mg.len() - 1,
        fdr_concorda(0.05),
        fdr_concorda(0.25)
    );
    t.check(rho_n > 0.999 && rho_p > 0.999, &format!("{nome}: gene RRA scores rank like MAGeCK's (Spearman > 0.999)"));
    t.check(fdr_concorda(0.05) == fa.len() && fdr_concorda(0.25) == fa.len(), &format!("{nome}: the same genes pass FDR 0.05 and 0.25 as in MAGeCK"));
}

fn dados_reais(t: &mut T, dir: &str, ref_dir: Option<&str>) {
    let d = Path::new(dir);
    let demo = d.join("demo");
    let p = |x: &Path| x.to_string_lossy().to_string();
    println!("\nReal data in {dir}:");

    // demo2: counting (MAGeCK printed: 2500 reads; 1453 / 1471 mapped; 1276 / 1199 zero; Gini 0.5267 / 0.4931)
    let c = crispr_count(
        &p(&demo.join("demo2/library.txt")),
        &[vec![p(&demo.join("demo2/test1.fastq"))], vec![p(&demo.join("demo2/test2.fastq"))]],
        &rotulos(&["L1", "CTRL"]),
        Some(&t.p("demo2.count.txt")),
        &CountOptions::default(),
    )
    .unwrap();
    print!("{}", c.resumo());
    let s: Vec<(u64, u64, usize, String)> = c.samples.iter().map(|x| (x.reads, x.mapped, x.zero_guides, py_g(x.gini, 4))).collect();
    t.eq(&s, &vec![(2500, 1453, 1276, "0.5267".to_string()), (2500, 1471, 1199, "0.4931".to_string())], "demo2: reads, mapped, zero-count guides and Gini = MAGeCK's countsummary");
    t.check(c.samples.iter().all(|x| x.files[0].offsets == vec![0]), "demo2: offset 0 found, as MAGeCK");
    t.eq(&c.library.ids.len(), &2550, "demo2: 2,550 guides in the library");
    let um = crispr_count(&p(&demo.join("demo2/library.txt")), &[vec![p(&demo.join("demo2/test1.fastq"))]], &[], None, &CountOptions { mismatches: 1, ..Default::default() }).unwrap();
    println!("  demo2 L1 with one mismatch: {} exact + {} with one mismatch, {} ambiguous", um.samples[0].mapped_exact, um.samples[0].mapped_mismatch, um.samples[0].ambiguous);
    t.check(um.samples[0].mapped_exact == 1453 && um.samples[0].mapped_mismatch > 0, "demo2: one mismatch adds reads to the exact ones");
    if let Some(r) = ref_dir {
        let mc = Path::new(r).join("demo2_count/demo.count.txt");
        if mc.exists() {
            let m = read_count_table(&p(&mc)).unwrap();
            let nossa = read_count_table(&t.p("demo2.count.txt")).unwrap();
            let de: HashMap<&str, (f64, f64)> = m.ids.iter().enumerate().map(|(i, id)| (id.as_str(), (m.counts[0][i], m.counts[1][i]))).collect();
            let iguais = nossa.ids.iter().enumerate().filter(|(i, id)| de.get(id.as_str()) == Some(&(nossa.counts[0][*i], nossa.counts[1][*i]))).count();
            t.check(iguais == 2550 && m.ids.len() == 2550, &format!("demo2: count of every guide = MAGeCK's demo.count.txt ({iguais} of {})", m.ids.len()));
        }
    }

    // demo1: the test
    let tab = p(&demo.join("demo1/sample.txt"));
    let trat = rotulos(&["HL60.final", "KBM7.final"]);
    let ctrl = rotulos(&["HL60.initial", "KBM7.initial"]);
    let o = TestOptions::default();
    let t0 = Instant::now();
    let r1 = crispr_test(&tab, &trat, &ctrl, Some(&t.p("demo1")), &o).unwrap();
    print!("{}", r1.resumo());
    println!("  ({:.2} s)", t0.elapsed().as_secs_f64());
    let f: Vec<f64> = r1.size_factors.iter().map(|x| x.1).collect();
    let mf = [1.469870985815957, 1.7759474888175155, 0.6308897693810006, 0.5721813161032618];
    t.check(f.iter().zip(&mf).all(|(a, b)| ((a - b) / b).abs() < 1e-13), "demo1: size factors = MAGeCK's log (1.4699, 1.7759, 0.6309, 0.5722)");
    t.perto(r1.alpha_neg, 0.3563563563563564, 1e-15, "demo1: negative alpha = MAGeCK's RRA -p 0.3563563563563564");
    t.perto(r1.alpha_pos, 0.36936936936936937, 1e-15, "demo1: positive alpha = MAGeCK's RRA -p 0.36936936936936937");
    // MAGeCK's own report for demo1 (demo5/demo_summary.Rnw, control normalization)
    let ctl = read_name_list(&p(&demo.join("demo5/control_sgrna.txt"))).unwrap();
    let rc = crispr_test(&tab, &trat, &ctrl, Some(&t.p("demo1_control")), &TestOptions { norm_method: NormMethod::Control, control_sgrnas: ctl.clone(), ..o.clone() }).unwrap();
    let f: Vec<f64> = rc.size_factors.iter().map(|x| x.1).collect();
    let rnw = [1.582701653316349, 1.7345173601833033, 0.5375843168437389, 0.5105538115034376];
    t.check(f.iter().zip(&rnw).all(|(a, b)| ((a - b) / b).abs() < 1e-13), "demo1, control normalization: size factors = MAGeCK's demo5 report (1.5827, 1.7345, 0.5376, 0.5106)");
    let neg: Vec<String> = rc.top(10, false).iter().map(|g| g.gene.clone()).collect();
    let pos: Vec<String> = rc.top(10, true).iter().map(|g| g.gene.clone()).collect();
    let rnw_neg = rotulos(&["ACTR8", "ACIN1", "AHCY", "AATF", "ABT1", "AHCTF1", "ACLY", "ADIRF", "AGBL5", "ACTR3"]);
    let rnw_pos = rotulos(&["ACRC", "ABCB8", "ADRA1A", "AGTPBP1", "AHRR", "ADCK1", "ACTR1A", "ADNP2", "ADK", "ADRBK1"]);
    println!("  demo1, control normalization, top 10 negative: {}\n    MAGeCK's demo5 report:                    {}", neg.join(" "), rnw_neg.join(" "));
    println!("  demo1, control normalization, top 10 positive: {}\n    MAGeCK's demo5 report:                    {}", pos.join(" "), rnw_pos.join(" "));
    let em_comum = |a: &[String], b: &[String]| a.iter().filter(|g| b.contains(g)).count();
    t.check(em_comum(&neg, &rnw_neg) >= 9, &format!("demo1: top 10 negative genes vs MAGeCK's shipped report ({} of 10 in common)", em_comum(&neg, &rnw_neg)));
    t.check(em_comum(&pos, &rnw_pos) >= 9, &format!("demo1: top 10 positive genes vs MAGeCK's shipped report ({} of 10 in common)", em_comum(&pos, &rnw_pos)));
    // the normalized counts of ACTR8's first guide in that report
    let actr8 = rc.guides.iter().find(|g| g.sgrna == "ACTR8_m53916047").unwrap();
    let alvo = [2674.7657941046295, 2294.76646752251, 1851.4403872098367, 812.8016679134727];
    let nossos: Vec<f64> = actr8.control_counts.iter().chain(&actr8.treatment_counts).copied().collect();
    t.check(nossos.iter().zip(&alvo).all(|(a, b)| ((a - b) / b).abs() < 1e-13), "demo1: normalized counts of ACTR8_m53916047 = MAGeCK's demo5 report");

    // other modes of demo1, and demo2's test (one control, one treatment)
    let r_tot = crispr_test(&tab, &trat, &ctrl, None, &TestOptions { norm_method: NormMethod::Total, ..o.clone() }).unwrap();
    let r_um = crispr_test(&tab, &trat, &rotulos(&["HL60.initial"]), None, &o).unwrap();
    let r_d2 = crispr_test(&t.p("demo2.count.txt"), &rotulos(&["L1"]), &rotulos(&["CTRL"]), None, &o).unwrap();
    t.eq(&r_d2.norm_used, &"total".to_string(), "demo2 test: too many zero-count guides, median falls back to total (as MAGeCK)");

    // the leukemia table: 73,151 guides
    let lk = p(&d.join("leukemia.new.csv"));
    let t0 = Instant::now();
    let rl = crispr_test(&lk, &trat, &ctrl, Some(&t.p("leuk")), &o).unwrap();
    print!("{}", rl.resumo());
    println!("  ({:.2} s)", t0.elapsed().as_secs_f64());
    let f: Vec<f64> = rl.size_factors.iter().map(|x| x.1).collect();
    let ml = [1.4168444424340858, 1.694643621525505, 0.6479653176768615, 0.5966574481207096];
    t.check(f.iter().zip(&ml).all(|(a, b)| ((a - b) / b).abs() < 1e-13), "leukemia: size factors = MAGeCK's log");
    t.eq(&rl.guides_total, &73151, "leukemia: 73,151 guides");
    // Brunello (77,441 guides): reads with the guide after the U6 promoter
    // tail, as sequenced; the offset is found and every read counted
    let bru = d.join("libraries/broadgpp-brunello-library-corrected.txt");
    if bru.exists() {
        let lib = read_library(&p(&bru), false).unwrap();
        t.eq(&(lib.ids.len(), lib.gene_count(), lib.header, lib.duplicate_sequences), &(77441, 19113, true, 0), "Brunello: 77,441 guides of 19,113 genes, header detected");
        let mut rr = Rng(5);
        let mut esperado = vec![0u64; lib.ids.len()];
        let reads: Vec<String> = (0..400_000)
            .map(|_| {
                let g = rr.below(lib.ids.len() as u64) as usize;
                esperado[g] += 1;
                format!("TTGTGGAAAGGACGAAACACCG{}GTTTTAGAGCTAGAAATAGC", lib.sequences[g])
            })
            .collect();
        grava_fastq(&t.p("brunello.fastq.gz"), &reads);
        let t0 = Instant::now();
        let c = crispr_count(&p(&bru), &[vec![t.p("brunello.fastq.gz")]], &[], None, &CountOptions { mismatches: 1, ..Default::default() }).unwrap();
        println!("  Brunello: 400,000 reads counted with one mismatch allowed in {:.2} s (offset {:?})", t0.elapsed().as_secs_f64(), c.samples[0].files[0].offsets);
        t.check(c.samples[0].files[0].offsets == vec![22] && c.counts[0] == esperado, "Brunello: offset 22 found, every read on its guide");
    }
    // MMG_CRISPR_DUMP=<dir>: keep the JSON of each result, to check it with a JSON parser
    if let Ok(dump) = std::env::var("MMG_CRISPR_DUMP") {
        for (nome, j) in [("count_demo2.json", c.json()), ("test_demo1.json", r1.json()), ("test_demo1_control.json", rc.json()), ("test_leukemia.json", rl.json()), ("test_leukemia_compact.json", rl.json_compact()), ("library_demo2.json", c.library.json())] {
            std::fs::write(Path::new(&dump).join(nome), j).unwrap();
        }
    }
    let topo: Vec<String> = rl.top(3, false).iter().map(|g| g.gene.clone()).collect();
    t.eq(&topo, &rotulos(&["RPL3", "RPL8", "RPS2"]), "leukemia: top 3 depleted genes = MAGeCK's (RPL3, RPL8, RPS2)");
    let na = rl.genes.iter().find(|g| g.gene == "NA").unwrap();
    t.check(na.pos_rank == 1 && na.pos_p == 1.0 && rl.top(1, true)[0].gene != "NA", "leukemia: gene NA stays in the table (positive rank 1, p = 1, as in MAGeCK) but out of the top list");

    if let Some(r) = ref_dir {
        let r = Path::new(r);
        println!("\nComparison with MAGeCK 0.5.9.5 outputs in {}:", r.display());
        compara_mageck(t, "demo1 (median)", &r1, &r.join("demo1"), "demo");
        compara_mageck(t, "demo1 (control sgRNAs)", &rc, &r.join("demo1_control"), "demo");
        compara_mageck(t, "demo1 (total)", &r_tot, &r.join("demo1_total"), "demo");
        compara_mageck(t, "demo1 (one control sample)", &r_um, &r.join("demo1_onecontrol"), "demo");
        // on MAGeCK's own count table, whose rows are in the order the reads
        // came (ours are in library order), so tied guides and genes keep the
        // same order
        let mc = r.join("demo2_count/demo.count.txt");
        if mc.exists() {
            let r_d2m = crispr_test(&p(&mc), &rotulos(&["L1"]), &rotulos(&["CTRL"]), None, &o).unwrap();
            t.eq(&r_d2m.genes.len(), &r_d2.genes.len(), "demo2 test: the same genes from our count table and from MAGeCK's");
            compara_mageck(t, "demo2 test (1 vs 1)", &r_d2m, &r.join("demo2_test"), "demo");
        }
        compara_mageck(t, "leukemia", &rl, &r.join("leukemia"), "leuk");
    }
}

fn main() {
    let dir = std::env::temp_dir().join(format!("test_crispr_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut t = T { ok: 0, falhas: 0, dir };
    // MMG_CRISPR_SEED=<n> draws other synthetic data (the checks must hold for any seed)
    let semente = std::env::var("MMG_CRISPR_SEED").ok().and_then(|s| s.parse().ok()).unwrap_or(20261114u64);
    let mut r = Rng(semente);
    let t0 = Instant::now();

    testa_formatos(&mut t);
    testa_normal(&mut t);
    testa_rra_matematica(&mut t, &mut r);
    testa_normalizacao(&mut t);
    testa_biblioteca(&mut t);
    testa_contagem(&mut t, &mut r);
    testa_tabela(&mut t);
    testa_tela(&mut t, &mut r);

    let args: Vec<String> = std::env::args().skip(1).collect();
    if let Some(d) = args.first() {
        dados_reais(&mut t, d, args.get(1).map(String::as_str));
    }

    let _ = std::fs::remove_dir_all(&t.dir);
    println!("{} of {} checks passed ({:.1} s)", t.ok, t.ok + t.falhas, t0.elapsed().as_secs_f64());
    if t.falhas > 0 {
        std::process::exit(1);
    }
}
