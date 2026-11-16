#![allow(dead_code)]
// Tests of src/associacao.rs: genotype-outcome association (linear, logistic,
// Firth), permutations (per-variant and family-wise), FDR, and power.
//
// Every result is compared with something computed independently:
//  - special functions: values from scipy (embedded);
//  - OLS and logistic regression: this file's own implementations (normal
//    equations / Newton with Gauss-Jordan inverses), written apart from the
//    module's Gram-Schmidt and Frisch-Waugh path; and values from statsmodels
//    0.15 on deterministic data sets (embedded; `dump` writes them);
//  - Firth: the modified score at the estimate, the penalized likelihood
//    computed here, and values from logistf (via firthlogist's tests) and an
//    independent scipy maximization (embedded);
//  - the files: hand-counted genotypes, a split into several VCFs that must
//    give the same rows, indexed against streamed reading;
//  - permutations: uniformity under the null and family-wise error rates in
//    simulation;
//  - power: simulation.
//
// Usage: cargo run --release --bin associacao            (MMG_VERBOSE=1 prints every check)
//        cargo run --release --bin associacao -- run vcfs=a.vcf,b.vcf pheno=p.tsv outcome=y [covariates=a,b] [categorical=c]
//              [type=auto] [model=additive] [firth=auto] [permutations=0] [seed=1] [ids=..] [regions=..]
//              [min_maf=0.01] [min_call_rate=0.9] [absent_is_ref=0] [case_value=..] [output=x.tsv] [json=out.json]
//        cargo run --release --bin associacao -- power type=binary freq=0.25 [effect=2] [n=120] [case_fraction=0.5]
//              [alpha=0.05] [tests=1] [model=additive] [covariates=0] [target=0.8] [simulations=0] [seed=1]
//        cargo run --release --bin associacao -- fit <csv: y,x1,x2...> <logit|firth>
//        cargo run --release --bin associacao -- dump <dir>
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
#[path = "../../../src/associacao.rs"]
mod associacao;

use associacao::*;
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
            if std::env::var_os("MMG_VERBOSE").is_some() {
                println!("ok    {what}");
            }
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
    /// |got - want| <= tol * max(|want|, floor)
    fn perto(&mut self, got: f64, want: f64, tol: f64, what: &str) -> bool {
        self.perto_piso(got, want, tol, 1e-300, what)
    }
    fn perto_piso(&mut self, got: f64, want: f64, tol: f64, piso: f64, what: &str) -> bool {
        let ok = (got.is_nan() && want.is_nan())
            || (got.is_infinite() && want.is_infinite() && got.signum() == want.signum())
            || (got - want).abs() <= tol * want.abs().max(piso);
        if !ok {
            println!("FAIL  {what}\n   got:  {got:e}\n   want: {want:e} (relative tolerance {tol:e})");
        }
        self.check(ok, what)
    }
    fn erro_contem(&mut self, r: Result<AssocReport, String>, trecho: &str, what: &str) {
        match r {
            Ok(_) => {
                self.check(false, &format!("{what}: expected an error with \"{trecho}\", got Ok"));
            }
            Err(e) => {
                let ok = e.contains(trecho);
                if !ok {
                    println!("FAIL  {what}: error does not contain \"{trecho}\":\n   {e}");
                }
                self.check(ok, what);
            }
        }
    }
}

// ---------------------------------------------------------------- independent numerics

/// Gauss-Jordan inverse with partial pivoting (row-major p x p).
fn inversa(a: &[f64], p: usize) -> Option<Vec<f64>> {
    let mut m = a.to_vec();
    let mut inv = vec![0.0; p * p];
    for i in 0..p {
        inv[i * p + i] = 1.0;
    }
    for c in 0..p {
        let piv = (c..p).max_by(|&i, &j| m[i * p + c].abs().total_cmp(&m[j * p + c].abs()))?;
        if m[piv * p + c].abs() < 1e-14 {
            return None;
        }
        for k in 0..p {
            m.swap(c * p + k, piv * p + k);
            inv.swap(c * p + k, piv * p + k);
        }
        let d = m[c * p + c];
        for k in 0..p {
            m[c * p + k] /= d;
            inv[c * p + k] /= d;
        }
        for i in 0..p {
            if i != c {
                let f = m[i * p + c];
                if f != 0.0 {
                    for k in 0..p {
                        m[i * p + k] -= f * m[c * p + k];
                        inv[i * p + k] -= f * inv[c * p + k];
                    }
                }
            }
        }
    }
    Some(inv)
}

/// OLS by the normal equations: (beta, se, t, df) of every coefficient.
/// `x`: rows (people) of p values (the intercept included).
fn ols(x: &[Vec<f64>], y: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>, f64) {
    let p = x[0].len();
    let n = x.len();
    let mut xtx = vec![0.0; p * p];
    let mut xty = vec![0.0; p];
    for (r, yi) in x.iter().zip(y) {
        for i in 0..p {
            xty[i] += r[i] * yi;
            for j in 0..p {
                xtx[i * p + j] += r[i] * r[j];
            }
        }
    }
    let inv = inversa(&xtx, p).expect("singular OLS design");
    let beta: Vec<f64> = (0..p).map(|i| (0..p).map(|j| inv[i * p + j] * xty[j]).sum()).collect();
    let rss: f64 = x.iter().zip(y).map(|(r, yi)| (yi - (0..p).map(|j| r[j] * beta[j]).sum::<f64>()).powi(2)).sum();
    let df = (n - p) as f64;
    let s2 = rss / df;
    let se: Vec<f64> = (0..p).map(|i| (s2 * inv[i * p + i]).sqrt()).collect();
    let t = beta.iter().zip(&se).map(|(b, s)| b / s).collect();
    (beta, se, t, df)
}

/// Logistic regression by plain Newton: (beta, se, loglik).
fn logit_indep(x: &[Vec<f64>], y: &[f64]) -> (Vec<f64>, Vec<f64>, f64) {
    let p = x[0].len();
    let mut beta = vec![0.0; p];
    let ll = |b: &[f64]| -> f64 {
        x.iter()
            .zip(y)
            .map(|(r, yi)| {
                let e: f64 = (0..p).map(|j| r[j] * b[j]).sum();
                yi * e - (1.0 + e.exp()).ln()
            })
            .sum()
    };
    let mut inv = vec![0.0; p * p];
    for _ in 0..200 {
        let mut h = vec![0.0; p * p];
        let mut g = vec![0.0; p];
        for (r, yi) in x.iter().zip(y) {
            let e: f64 = (0..p).map(|j| r[j] * beta[j]).sum();
            let pi = 1.0 / (1.0 + (-e).exp());
            let w = pi * (1.0 - pi);
            for i in 0..p {
                g[i] += r[i] * (yi - pi);
                for j in 0..p {
                    h[i * p + j] += w * r[i] * r[j];
                }
            }
        }
        inv = inversa(&h, p).unwrap();
        let passo: Vec<f64> = (0..p).map(|i| (0..p).map(|j| inv[i * p + j] * g[j]).sum()).collect();
        for i in 0..p {
            beta[i] += passo[i];
        }
        if passo.iter().all(|s| s.abs() < 1e-13) {
            break;
        }
    }
    let se = (0..p).map(|i| inv[i * p + i].sqrt()).collect();
    (beta.clone(), se, ll(&beta))
}

/// Firth's penalized log-likelihood l(b) + 1/2 log det(X'WX), computed here.
fn ll_firth(x: &[Vec<f64>], y: &[f64], b: &[f64]) -> f64 {
    let p = x[0].len();
    let mut inf = vec![0.0; p * p];
    let mut l = 0.0;
    for (r, yi) in x.iter().zip(y) {
        let e: f64 = (0..p).map(|j| r[j] * b[j]).sum();
        let pi = 1.0 / (1.0 + (-e).exp());
        l += yi * e - (1.0 + e.exp()).ln();
        for i in 0..p {
            for j in 0..p {
                inf[i * p + j] += pi * (1.0 - pi) * r[i] * r[j];
            }
        }
    }
    // determinant by elimination
    let mut m = inf;
    let mut det = 1.0;
    for c in 0..p {
        let piv = (c..p).max_by(|&i, &j| m[i * p + c].abs().total_cmp(&m[j * p + c].abs())).unwrap();
        if piv != c {
            for k in 0..p {
                m.swap(c * p + k, piv * p + k);
            }
            det = -det;
        }
        let d = m[c * p + c];
        det *= d;
        for i in c + 1..p {
            let f = m[i * p + c] / d;
            for k in c..p {
                m[i * p + k] -= f * m[c * p + k];
            }
        }
    }
    l + 0.5 * det.ln()
}

/// The modified score of Firth's method at b (hat values by an explicit inverse).
fn score_firth(x: &[Vec<f64>], y: &[f64], b: &[f64]) -> Vec<f64> {
    let p = x[0].len();
    let mut inf = vec![0.0; p * p];
    let pis: Vec<f64> = x
        .iter()
        .map(|r| {
            let e: f64 = (0..p).map(|j| r[j] * b[j]).sum();
            1.0 / (1.0 + (-e).exp())
        })
        .collect();
    for (r, pi) in x.iter().zip(&pis) {
        for i in 0..p {
            for j in 0..p {
                inf[i * p + j] += pi * (1.0 - pi) * r[i] * r[j];
            }
        }
    }
    let inv = inversa(&inf, p).unwrap();
    let mut u = vec![0.0; p];
    for ((r, yi), pi) in x.iter().zip(y).zip(&pis) {
        let mut h = 0.0;
        for i in 0..p {
            for j in 0..p {
                h += r[i] * inv[i * p + j] * r[j];
            }
        }
        h *= pi * (1.0 - pi);
        for i in 0..p {
            u[i] += r[i] * (yi - pi + h * (0.5 - pi));
        }
    }
    u
}

fn cols_de(x: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let p = x[0].len();
    (0..p).map(|j| x.iter().map(|r| r[j]).collect()).collect()
}

// ---------------------------------------------------------------- data sets

#[derive(Clone)]
struct Variante {
    chrom: String,
    pos: u64,
    id: String,
    refe: String,
    alts: Vec<String>,
    /// GT text per person.
    gts: Vec<String>,
}

fn escreve_vcf(caminho: &str, amostras: &[String], vars: &[Variante], cols: Option<&[usize]>) {
    let mut f = std::io::BufWriter::new(std::fs::File::create(caminho).unwrap());
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(f, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">").unwrap();
    let idx: Vec<usize> = cols.map(|c| c.to_vec()).unwrap_or_else(|| (0..amostras.len()).collect());
    let nomes: Vec<&str> = idx.iter().map(|&i| amostras[i].as_str()).collect();
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", nomes.join("\t")).unwrap();
    for v in vars {
        let g: Vec<&str> = idx.iter().map(|&i| v.gts[i].as_str()).collect();
        writeln!(f, "{}\t{}\t{}\t{}\t{}\t.\tPASS\t.\tGT\t{}", v.chrom, v.pos, v.id, v.refe, v.alts.join(","), g.join("\t")).unwrap();
    }
}

fn escreve_texto(caminho: &str, texto: &str) {
    std::fs::write(caminho, texto).unwrap();
}

/// A cohort: genotypes (copies of ALT, -1 missing), covariates and outcomes.
struct Coorte {
    nomes: Vec<String>,
    idade: Vec<f64>,
    sexo: Vec<f64>,
    centro: Vec<String>,
    pc1: Vec<f64>,
    /// genotypes, variant-major: g[v][i] in {-1, 0, 1, 2}
    g: Vec<Vec<i32>>,
    yq: Vec<f64>,
    yb: Vec<f64>,
}

/// Additive genotypes at frequency f.
fn genotipo(rng: &mut Rng, f: f64) -> i32 {
    (rng.unif() < f) as i32 + (rng.unif() < f) as i32
}

fn gt_texto(g: i32, rng: &mut Rng) -> String {
    match g {
        -1 => "./.".into(),
        0 => "0/0".into(),
        1 => {
            if rng.unif() < 0.5 {
                "0|1".into()
            } else {
                "1|0".into()
            }
        }
        _ => "1/1".into(),
    }
}

/// n people, nv variants on chr1 (and chr2 for the second half); variant 0
/// has a real effect on both outcomes; variant 1's frequency follows pc1
/// (confounded: pc1 changes the outcomes); 4% missing genotypes.
fn coorte(seed: u64, n: usize, nv: usize, efeito_q: f64, efeito_b: f64, falta: f64) -> Coorte {
    let mut rng = Rng::new(seed);
    let nomes: Vec<String> = (0..n).map(|i| format!("P{:03}", i + 1)).collect();
    let idade: Vec<f64> = (0..n).map(|_| ((50.0 + 10.0 * rng.normal()) * 10.0).round() / 10.0).collect();
    let sexo: Vec<f64> = (0..n).map(|_| (rng.unif() < 0.5) as u8 as f64).collect();
    let centro: Vec<String> = (0..n).map(|_| ["A", "B", "C"][rng.below(3)].to_string()).collect();
    let pc1: Vec<f64> = (0..n).map(|_| (rng.normal() * 1000.0).round() / 1000.0).collect();
    let mut g = Vec::new();
    for v in 0..nv {
        let base = 0.1 + 0.4 * ((v * 7919) % 13) as f64 / 13.0;
        let col: Vec<i32> = (0..n)
            .map(|i| {
                let f = if v == 1 { (0.3 + 0.2 * pc1[i]).clamp(0.02, 0.98) } else { base };
                if rng.unif() < falta {
                    -1
                } else {
                    genotipo(&mut rng, f)
                }
            })
            .collect();
        g.push(col);
    }
    let efc = |c: &str| match c {
        "A" => 0.0,
        "B" => 0.4,
        _ => -0.3,
    };
    let yq: Vec<f64> = (0..n)
        .map(|i| {
            let x = g[0][i].max(0) as f64;
            let v = 1.0 + 0.02 * (idade[i] - 50.0) + 0.3 * sexo[i] + efc(&centro[i]) + 0.5 * pc1[i] + efeito_q * x + rng.normal();
            (v * 1e4).round() / 1e4
        })
        .collect();
    let yb: Vec<f64> = (0..n)
        .map(|i| {
            let x = g[0][i].max(0) as f64;
            let e = -0.5 + 0.03 * (idade[i] - 50.0) + 0.4 * sexo[i] + efc(&centro[i]) + 0.6 * pc1[i] + efeito_b * x;
            (rng.unif() < 1.0 / (1.0 + (-e).exp())) as u8 as f64
        })
        .collect();
    Coorte { nomes, idade, sexo, centro, pc1, g, yq, yb }
}

fn variantes_de(c: &Coorte, seed: u64) -> Vec<Variante> {
    let mut rng = Rng::new(seed ^ 77);
    let nv = c.g.len();
    (0..nv)
        .map(|v| {
            let (chrom, pos) = if v < nv / 2 { ("chr1", 1000 + 100 * v as u64) } else { ("chr2", 5000 + 100 * v as u64) };
            Variante {
                chrom: chrom.into(),
                pos,
                id: format!("rs{}", 100 + v),
                refe: "A".into(),
                alts: vec!["G".into()],
                gts: c.g[v].iter().map(|&x| gt_texto(x, &mut rng)).collect(),
            }
        })
        .collect()
}

fn escreve_tabela_coorte(caminho: &str, c: &Coorte, faltas: &[(usize, &str)]) {
    let mut s = String::from("sample\tyq\tyb\tage\tsex\tcenter\tpc1\n");
    for i in 0..c.nomes.len() {
        let mut campos = vec![
            c.nomes[i].clone(),
            format!("{}", c.yq[i]),
            format!("{}", c.yb[i] as u8),
            format!("{}", c.idade[i]),
            format!("{}", c.sexo[i] as u8),
            c.centro[i].clone(),
            format!("{}", c.pc1[i]),
        ];
        for (k, col) in faltas {
            if *k == i {
                let j = ["sample", "yq", "yb", "age", "sex", "center", "pc1"].iter().position(|h| h == col).unwrap();
                campos[j] = "NA".into();
            }
        }
        s.push_str(&campos.join("\t"));
        s.push('\n');
    }
    escreve_texto(caminho, &s);
}

/// The design rows of the people analysed (complete phenotypes) with a
/// genotype at variant v; the coded genotype last.
fn desenho(c: &Coorte, v: usize, excluir: &[usize], codif: &dyn Fn(i32) -> f64, binario: bool) -> (Vec<Vec<f64>>, Vec<f64>) {
    let mut x = Vec::new();
    let mut y = Vec::new();
    for i in 0..c.nomes.len() {
        if excluir.contains(&i) || c.g[v][i] < 0 {
            continue;
        }
        x.push(vec![
            1.0,
            c.idade[i],
            c.sexo[i],
            (c.centro[i] == "B") as u8 as f64,
            (c.centro[i] == "C") as u8 as f64,
            c.pc1[i],
            codif(c.g[v][i]),
        ]);
        y.push(if binario { c.yb[i] } else { c.yq[i] });
    }
    (x, y)
}

fn opcoes(outcome: &str) -> AssocOptions {
    AssocOptions {
        outcome: outcome.into(),
        covariates: vec!["age".into(), "sex".into(), "center".into(), "pc1".into()],
        min_maf: 0.0,
        min_call_rate: 0.0,
        ..Default::default()
    }
}

// ---------------------------------------------------------------- tests

fn testa_funcoes(t: &mut T) {
    for (x, want) in [
        (0.5, 0.3085375387259869),
        (1.96, 0.024997895148220435),
        (5.0, 2.866515718791933e-07),
        (10.0, 7.61985302416047e-24),
        (30.0, 4.906713927147908e-198),
        (-2.5, 0.9937903346742238),
    ] {
        t.perto(normal_upper(x), want, 1e-13, &format!("normal_upper({x})"));
    }
    for (p, want) in [
        (1e-300, -37.0470962993612),
        (1e-10, -6.361340902404056),
        (0.025, -1.9599639845400545),
        (0.3, -0.5244005127080409),
        (0.975, 1.959963984540054),
        (0.999999, 4.753424308817087),
    ] {
        t.perto(normal_quantile(p), want, 1e-13, &format!("normal_quantile({p})"));
    }
    t.eq(&normal_quantile(0.5), &0.0, "normal_quantile(0.5)");
    for (tt, df, want) in [
        (2.0, 5.0, 0.10193947882985835),
        (2.228138851986274, 10.0, 0.050000000000000044),
        (10.0, 3.0, 0.0021283990584141503),
        (50.0, 100.0, 1.4472163679761465e-72),
        (0.1, 1.0, 0.936548965138893),
        (40.0, 2000.0, 1.4278618533562838e-257),
        (-3.3, 17.0, 0.004230097403725666),
    ] {
        t.perto(t_two_sided(tt, df), want, 1e-10, &format!("t_two_sided({tt}, {df})"));
    }
    for (df, a, want) in [
        (1.0, 0.05, 12.706204736174705),
        (5.0, 0.05, 2.5705818356363155),
        (30.0, 1e-06, 6.119075620384216),
        (1000.0, 0.05, 1.9623390808264087),
        (3.0, 0.2, 1.6377443536962089),
    ] {
        t.perto(t_critical(df, a), want, 1e-11, &format!("t_critical({df}, {a})"));
    }
    for (x, want) in [(0.5, 0.5723649429247), (3.7, 1.428072326665388), (50.5, 146.51925549072064), (100000.0, 1051287.7089736569), (0.01, 4.599479878042022)] {
        t.perto(ln_gamma(x), want, 1e-13, &format!("ln_gamma({x})"));
    }
    t.check(ln_gamma(1.0).abs() < 1e-14, "ln_gamma(1) = 0");
    for (a, b, x, lo, up) in [
        (2.5, 0.5, 0.3, 0.018927124071945658, 0.9810728759280543),
        (50.0, 0.5, 0.01, 7.998227417904836e-102, 1.0),
        (0.5, 0.5, 0.9, 0.7951672353008665, 0.20483276469913342),
        (3.0, 7.0, 0.25, 0.399322509765625, 0.600677490234375),
    ] {
        let (l, u) = beta_reg(a, b, x, 1.0 - x);
        t.perto(l, lo, 1e-11, &format!("beta_reg lower ({a}, {b}, {x})"));
        t.perto(u, up, 1e-11, &format!("beta_reg upper ({a}, {b}, {x})"));
    }
    for (df, d, c, want) in [
        (3.0, 2.0, 3.1824463052837086, 0.2887524164017148),
        (10.0, 3.0, 2.2281388519862753, 0.7714014353512886),
        (58.0, 2.8, 2.001717484145236, 0.7862286390766343),
        (1.0, 5.0, 12.706204736174705, 0.3051601027351557),
        (500.0, 4.5, 5.5361563774683615, 0.15435886060092893),
        (200.0, 6.0, 5.6677397680687145, 0.6280337609782534),
        (40.0, 7.0, 6.69489551093783, 0.6109042969167199),
        (118.0, 1.2, 3.3070141321536544, 0.020021864837709093),
    ] {
        t.perto(noncentral_t_power(df, d, c), want, 1e-8, &format!("noncentral_t_power({df}, {d}, {c})"));
    }
    for (x, want) in [(3.841458820694124, 0.04999999999999989), (100.0, 1.5239706048320995e-23), (0.5, 0.47950012218695337), (1e-6, 0.9992021155721779)] {
        t.perto(chi2_1_upper(x), want, 1e-12, &format!("chi2_1_upper({x})"));
    }
    // BH against hand values (R: p.adjust(c(0.01,0.04,0.03,0.2,NA,0.04), "BH"))
    let bh = benjamini_hochberg(&[0.01, 0.04, 0.03, 0.2, f64::NAN, 0.04]);
    let want = [0.05, 0.05, 0.05, 0.2, f64::NAN, 0.05];
    for (i, (g, w)) in bh.iter().zip(want).enumerate() {
        t.perto(*g, w, 1e-14, &format!("BH[{i}]"));
    }
    let bh = benjamini_hochberg(&[0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205]);
    let want = [0.008, 0.032, 0.0672, 0.0672, 0.0672, 0.08, 0.0845714285714286, 0.205];
    for (i, (g, w)) in bh.iter().zip(want).enumerate() {
        t.perto(*g, w, 1e-12, &format!("BH (8 values)[{i}]"));
    }
}

fn testa_logistica(t: &mut T) {
    // a data set with two covariates and a genotype
    let mut rng = Rng::new(11);
    let n = 300;
    let mut x = Vec::new();
    let mut y = Vec::new();
    for _ in 0..n {
        let a = rng.normal();
        let s = (rng.unif() < 0.4) as u8 as f64;
        let g = genotipo(&mut rng, 0.3) as f64;
        let e = -0.3 + 0.5 * a - 0.4 * s + 0.6 * g;
        x.push(vec![1.0, a, s, g]);
        y.push((rng.unif() < 1.0 / (1.0 + (-e).exp())) as u8 as f64);
    }
    let cols = cols_de(&x);
    let f = logistic_fit(&cols, &y, None);
    let (b, se, ll) = logit_indep(&x, &y);
    t.check(f.converged && !f.separated, "logistic: converged, not separated");
    for j in 0..4 {
        t.perto(f.beta[j], b[j], 1e-9, &format!("logistic beta[{j}] = independent Newton"));
        t.perto(f.se[j], se[j], 1e-9, &format!("logistic se[{j}] = independent Newton"));
    }
    t.perto(f.loglik, ll, 1e-12, "logistic loglik = independent");

    // complete separation: every carrier has the event
    let mut x2 = Vec::new();
    let mut y2 = Vec::new();
    for i in 0..40 {
        let g = (i % 5 == 0) as u8 as f64;
        let yy = if g == 1.0 { 1.0 } else { (i % 3 == 0) as u8 as f64 };
        x2.push(vec![1.0, (i as f64) / 10.0, g]);
        y2.push(yy);
    }
    let f2 = logistic_fit(&cols_de(&x2), &y2, None);
    t.check(f2.separated, "logistic: quasi-complete separation is detected");

    // Firth on the separated data: finite estimate, zero modified score,
    // the maximum of the penalized likelihood
    let ff = firth_fit(&cols_de(&x2), &y2, None, None).expect("firth fit");
    t.check(ff.converged, "firth: converged under separation");
    t.check(ff.beta.iter().all(|b| b.is_finite() && b.abs() < 20.0), "firth: finite estimates under separation");
    let u = score_firth(&x2, &y2, &ff.beta);
    t.check(u.iter().all(|v| v.abs() < 1e-7), &format!("firth: modified score is zero at the estimate ({u:?})"));
    t.perto(ff.loglik_pen, ll_firth(&x2, &y2, &ff.beta), 1e-12, "firth: penalized loglik = independent");
    let mut melhor = true;
    for j in 0..3 {
        for d in [-1e-3, 1e-3] {
            let mut b = ff.beta.clone();
            b[j] += d;
            if ll_firth(&x2, &y2, &b) > ff.loglik_pen + 1e-12 {
                melhor = false;
            }
        }
    }
    t.check(melhor, "firth: the estimate is a maximum of the penalized likelihood");
    // the restricted fit maximizes over the free coefficients
    let rf = firth_fit(&cols_de(&x2), &y2, Some((2, 0.0)), None).unwrap();
    t.eq(&rf.beta[2], &0.0, "firth restricted: the fixed coefficient stays");
    let u = score_firth(&x2, &y2, &rf.beta);
    t.check(u[0].abs() < 1e-7 && u[1].abs() < 1e-7, &format!("firth restricted: free scores are zero ({u:?})"));
    let (s, p) = firth_lrt(&cols_de(&x2), &y2, &ff, 2);
    t.perto(s, 2.0 * (ff.loglik_pen - ll_firth(&x2, &y2, &rf.beta)), 1e-9, "firth LRT statistic = 2 (l* - l*_0)");
    t.perto(p, chi2_1_upper(s), 1e-14, "firth LRT p");
    // the profile CI ends where 2 (l* - l*_profile) = 3.84
    let (lo, hi) = firth_profile_ci(&cols_de(&x2), &y2, &ff, 2, CHI2_1_95);
    for (b, nome) in [(lo, "lower"), (hi, "upper")] {
        let r = firth_fit(&cols_de(&x2), &y2, Some((2, b)), None).unwrap();
        t.perto(2.0 * (ff.loglik_pen - r.loglik_pen), CHI2_1_95, 1e-7, &format!("firth profile CI {nome} bound at the 3.84 contour"));
    }
    t.check(lo < ff.beta[2] && ff.beta[2] < hi, "firth profile CI contains the estimate");

    // endometrial-like quasi-separation with logistf's reference values is
    // checked in Python (the data are not redistributed here); a
    // deterministic data set with values from firthlogist and scipy:
    let (xf, yf) = dados_firth();
    let ff = firth_fit(&cols_de(&xf), &yf, None, None).unwrap();
    for (j, (b, se)) in REF_FIRTH_BETA.iter().zip(REF_FIRTH_SE).enumerate() {
        t.perto(ff.beta[j], *b, 1e-6, &format!("firth beta[{j}] = firthlogist/scipy"));
        t.perto(ff.se[j], se, 1e-6, &format!("firth se[{j}] = firthlogist/scipy"));
    }
    let j = xf[0].len() - 1;
    let (_, p) = firth_lrt(&cols_de(&xf), &yf, &ff, j);
    t.perto(p, REF_FIRTH_P_LRT, 1e-5, "firth penalized LRT p = scipy (full-design penalty, as logistf)");
    let (lo, hi) = firth_profile_ci(&cols_de(&xf), &yf, &ff, j, CHI2_1_95);
    t.perto(lo, REF_FIRTH_CI.0, 1e-5, "firth profile CI lower = scipy");
    t.perto(hi, REF_FIRTH_CI.1, 1e-5, "firth profile CI upper = scipy");
}

/// A small deterministic data set with quasi-complete separation in the last
/// column (all 6 carriers have the event): y, intercept, a, b, g.
fn dados_firth() -> (Vec<Vec<f64>>, Vec<f64>) {
    let mut rng = Rng::new(2024);
    let mut x = Vec::new();
    let mut y = Vec::new();
    for i in 0..30 {
        let a = (rng.normal() * 100.0).round() / 100.0;
        let b = (rng.unif() < 0.5) as u8 as f64;
        let g = (i % 5 == 0) as u8 as f64;
        let e = -0.4 + 0.8 * a + 0.5 * b;
        let yy = if g == 1.0 { 1.0 } else { (rng.unif() < 1.0 / (1.0 + (-e).exp())) as u8 as f64 };
        x.push(vec![1.0, a, b, g]);
        y.push(yy);
    }
    (x, y)
}

// An independent maximization of the penalized likelihood with scipy (BFGS
// and Nelder-Mead on l + 1/2 log det I; restricted fits with the full-design
// penalty; profile bounds by Brent's method), on `dump`'s firth.csv. The
// firthlogist 0.5.0 code gives the same estimates (to 1e-8) and CI, and a
// penalized-LRT p of 0.021949: it computes the hat values of the restricted fit
// without the fixed column, a small departure from logistf. On logistf's own
// data (endometrial, sex2; values from firthlogist's tests) this module's
// coefficients and profile CIs agree with logistf to 3e-7, and its
// penalized-LRT p-values with scipy to 3e-14 (endometrial NV: p = 0.009124).
const REF_FIRTH_BETA: [f64; 4] = [-2.678990480229634, 1.0642355223130322, 3.224109594744089, 3.0160764257165265];
const REF_FIRTH_SE: [f64; 4] = [1.4163405532685045, 0.6336003633575521, 1.4898814313063873, 1.866097704294576];
const REF_FIRTH_P_LRT: f64 = 0.021965390984583655;
const REF_FIRTH_CI: (f64, f64) = (0.363167384755351, 8.54544464737418);

fn linha<'a>(r: &'a AssocReport, id: &str) -> &'a AssocRow {
    r.rows.iter().find(|x| x.id == id).unwrap_or_else(|| panic!("row {id} missing"))
}

fn testa_ponta_a_ponta(t: &mut T) {
    let c = coorte(5, 160, 8, 0.7, 1.2, 0.04);
    let vars = variantes_de(&c, 5);
    let vcf = t.p("coorte.vcf");
    escreve_vcf(&vcf, &c.nomes, &vars, None);
    // two people without age, one without the outcome, two not in the VCF
    let tab = t.p("coorte.tsv");
    escreve_tabela_coorte(&tab, &c, &[(3, "age"), (10, "age"), (20, "yq")]);
    {
        let mut s = std::fs::read_to_string(&tab).unwrap();
        s.push_str("X001\t1.0\t1\t40\t0\tA\t0.1\nX002\t2.0\t0\t41\t1\tB\t0.2\n");
        escreve_texto(&tab, &s);
    }
    let excl = [3usize, 10, 20];

    // --- quantitative ---
    let r = genotype_association(&[vcf.clone()], &tab, &opcoes("yq")).unwrap();
    t.eq(&r.outcome_type, &"quantitative".to_string(), "quantitative outcome detected");
    t.eq(&r.analysed, &157, "people analysed (160 - 2 without age - 1 without outcome)");
    t.eq(&r.not_in_vcf, &2, "people in the table and not in the VCF");
    t.eq(&r.missing_outcome, &1, "missing outcome counted");
    t.eq(&r.missing_covariates, &vec![("age".to_string(), 2), ("sex".to_string(), 0), ("center".to_string(), 0), ("pc1".to_string(), 0)], "missing covariates counted");
    t.eq(&r.design, &vec!["intercept".to_string(), "age".into(), "sex".into(), "center=B".into(), "center=C".into(), "pc1".into()], "design columns (center: reference A)");
    t.eq(&r.tested, &8, "8 variants tested");
    for v in 0..8 {
        let (x, y) = desenho(&c, v, &excl, &|g| g as f64, false);
        let (b, se, tt, df) = ols(&x, &y);
        let row = linha(&r, &format!("rs{}", 100 + v));
        t.eq(&row.n, &x.len(), &format!("rs{} n", 100 + v));
        t.perto(row.beta, b[6], 1e-9, &format!("rs{} OLS beta = normal equations", 100 + v));
        t.perto(row.se, se[6], 1e-9, &format!("rs{} OLS se", 100 + v));
        t.perto(row.stat, tt[6], 1e-9, &format!("rs{} OLS t", 100 + v));
        t.perto(row.p, t_two_sided(tt[6], df), 1e-8, &format!("rs{} OLS p", 100 + v));
        let tc = t_critical(df, 0.05);
        t.perto(row.ci_low, b[6] - tc * se[6], 1e-9, &format!("rs{} OLS CI low", 100 + v));
        let freq = x.iter().map(|r| r[6]).sum::<f64>() / (2.0 * x.len() as f64);
        t.perto(row.alt_freq, freq, 1e-12, &format!("rs{} alt freq", 100 + v));
        t.perto(row.call_rate, x.len() as f64 / 157.0, 1e-12, &format!("rs{} call rate", 100 + v));
    }
    let p100 = linha(&r, "rs100").p;
    t.check(p100 < 1e-3, &format!("the planted effect (rs100) is found (p {p100:e})"));
    // dominant and recessive
    for (m, nome) in [(GeneticModel::Dominant, "dominant"), (GeneticModel::Recessive, "recessive")] {
        let mut o = opcoes("yq");
        o.model = m;
        let r = genotype_association(&[vcf.clone()], &tab, &o).unwrap();
        for v in 0..8usize {
            let cod = |g: i32| -> f64 {
                match m {
                    GeneticModel::Dominant => (g >= 1) as u8 as f64,
                    _ => (g == 2) as u8 as f64,
                }
            };
            let (x, y) = desenho(&c, v, &excl, &cod, false);
            let k = x.iter().filter(|r| r[6] == 1.0).count();
            if k == 0 || k == x.len() {
                t.check(!r.rows.iter().any(|z| z.id == format!("rs{}", 100 + v)), &format!("{nome} rs{}: no variation, not tested", 100 + v));
                continue;
            }
            let (b, se, ..) = ols(&x, &y);
            let row = linha(&r, &format!("rs{}", 100 + v));
            t.perto(row.beta, b[6], 1e-9, &format!("{nome} rs{} beta", 100 + v));
            t.perto(row.se, se[6], 1e-9, &format!("{nome} rs{} se", 100 + v));
        }
    }

    // --- binary ---
    let r = genotype_association(&[vcf.clone()], &tab, &opcoes("yb")).unwrap();
    t.eq(&r.outcome_type, &"binary".to_string(), "binary outcome detected");
    t.eq(&r.event_value, &Some("1".to_string()), "event = 1 (0/1 coding)");
    t.eq(&r.analysed, &158, "binary: 158 analysed (the yq gap does not count)");
    for v in 0..8 {
        let (x, y) = desenho(&c, v, &[3, 10], &|g| g as f64, true);
        let (b, se, ll1) = logit_indep(&x, &y);
        let x0: Vec<Vec<f64>> = x.iter().map(|r| r[..6].to_vec()).collect();
        let (_, _, ll0) = logit_indep(&x0, &y);
        let row = linha(&r, &format!("rs{}", 100 + v));
        t.eq(&row.method, &"logistic".to_string(), &format!("rs{} method", 100 + v));
        t.perto(row.beta, b[6], 1e-8, &format!("rs{} logistic beta = independent Newton", 100 + v));
        t.perto(row.se, se[6], 1e-8, &format!("rs{} logistic se", 100 + v));
        t.perto(row.p_wald, p_two_sided_z(b[6] / se[6]), 1e-7, &format!("rs{} logistic Wald p", 100 + v));
        t.perto(row.p_lrt, chi2_1_upper(2.0 * (ll1 - ll0)), 1e-7, &format!("rs{} logistic LRT p", 100 + v));
        t.perto(row.odds_ratio, b[6].exp(), 1e-8, &format!("rs{} OR", 100 + v));
        let casos = y.iter().filter(|v| **v == 1.0).count();
        t.eq(&row.n_cases, &casos, &format!("rs{} cases", 100 + v));
    }
    let p100 = linha(&r, "rs100").p;
    t.check(p100 < 0.01, &format!("binary: the planted effect is found (p {p100:e})"));
    // without the confounder, rs101 (frequency follows pc1) looks associated
    let mut o = opcoes("yq");
    o.covariates = vec!["age".into(), "sex".into(), "center".into()];
    let sem = genotype_association(&[vcf.clone()], &tab, &o).unwrap();
    let com = genotype_association(&[vcf.clone()], &tab, &opcoes("yq")).unwrap();
    t.check(linha(&sem, "rs101").p < 0.01 && linha(&com, "rs101").p > 0.05, &format!(
        "confounding: rs101 p {} without pc1, {} with it",
        linha(&sem, "rs101").p,
        linha(&com, "rs101").p
    ));

    // --- output file and JSON ---
    let mut o = opcoes("yb");
    o.output = Some(t.p("coorte_res.tsv"));
    let r = genotype_association(&[vcf.clone()], &tab, &o).unwrap();
    let txt = std::fs::read_to_string(t.p("coorte_res.tsv")).unwrap();
    let linhas: Vec<&str> = txt.lines().collect();
    t.eq(&linhas.len(), &9, "TSV: header + 8 rows");
    t.check(linhas[0].starts_with("id\tchrom\tpos\tref\talt\tn\talt_freq\tcall_rate\tbeta\tse\tci_low\tci_high\tor\tor_low\tor_high\tp\tmethod\tp_permutation\tp_fwer\tfdr"), "TSV header");
    let js: serde_json::Value = serde_json::from_str(&r.json()).expect("valid JSON");
    t.eq(&js["rows"].as_array().map(|a| a.len()), &Some(8), "JSON rows");
    t.eq(&js["people"]["analysed"].as_u64(), &Some(158), "JSON people.analysed");
    t.check(js["text"].as_str().map_or(false, |s| s.contains("not causation")), "JSON text");
    t.check(js["power"]["detectable_effect"].as_f64().map_or(false, |x| x > 1.0), "JSON power.detectable_effect is an odds ratio > 1");
    let p_rs = js["rows"][0]["p"].as_f64().unwrap();
    t.perto(p_rs, r.rows[0].p, 1e-12, "JSON p round-trips");
    // FDR is BH over the tested p-values
    let bh = benjamini_hochberg(&r.rows.iter().map(|x| x.p).collect::<Vec<_>>());
    t.check(r.rows.iter().zip(&bh).all(|(x, b)| (x.fdr - b).abs() < 1e-15), "fdr = BH of p");

    // --- several VCFs give the same rows ---
    let grupos: [Vec<usize>; 3] = [(0..50).collect(), (50..51).collect(), (51..160).collect()];
    let mut arqs = Vec::new();
    for (k, gr) in grupos.iter().enumerate() {
        let p = t.p(&format!("parte{k}.vcf"));
        // the middle file names chromosomes without "chr"
        let vs: Vec<Variante> = if k == 1 {
            vars.iter().map(|v| Variante { chrom: v.chrom.trim_start_matches("chr").to_string(), ..v.clone() }).collect()
        } else {
            vars.clone()
        };
        escreve_vcf(&p, &c.nomes, &vs, Some(gr));
        arqs.push(p);
    }
    let r1 = genotype_association(&[vcf.clone()], &tab, &opcoes("yq")).unwrap();
    let r3 = genotype_association(&arqs, &tab, &opcoes("yq")).unwrap();
    t.eq(&r3.tested, &r1.tested, "split VCFs: same number of tests");
    let mesmas = r1.rows.iter().zip(&r3.rows).all(|(a, b)| a.id == b.id && a.n == b.n && (a.beta - b.beta).abs() < 1e-12 && (a.p - b.p).abs() <= 1e-12 * a.p);
    t.check(mesmas, "split VCFs (one with `1` instead of `chr1`): identical rows");
    // variant-only files: records where every person of the file is 0/0 are left out
    let mut arqs2 = Vec::new();
    for (k, gr) in grupos.iter().enumerate() {
        let p = t.p(&format!("so_variantes{k}.vcf"));
        let vs: Vec<Variante> = vars.iter().filter(|v| gr.iter().any(|&i| v.gts[i] != "0/0")).cloned().collect();
        escreve_vcf(&p, &c.nomes, &vs, Some(gr));
        arqs2.push(p);
    }
    let mut o = opcoes("yq");
    o.absent_is_ref = true;
    let ra = genotype_association(&arqs2, &tab, &o).unwrap();
    let mesmas = r1.rows.iter().zip(&ra.rows).all(|(a, b)| a.id == b.id && a.n == b.n && (a.beta - b.beta).abs() < 1e-12);
    t.check(mesmas, "variant-only VCFs with absent_is_ref: identical rows");
    let rf = genotype_association(&arqs2, &tab, &opcoes("yq")).unwrap();
    let n_falta: usize = rf.rows.iter().map(|x| x.n).sum();
    let n_cheio: usize = r1.rows.iter().map(|x| x.n).sum();
    t.check(n_falta <= n_cheio, "variant-only VCFs without absent_is_ref: absent = missing");
    t.check(rf.warnings.iter().any(|w| w.contains("counts as missing")), "the absent-as-missing rule is stated");
    // the person in file 1 (P051): the second file has no record where P051 is 0/0
    let ausente: Vec<usize> = (0..8).filter(|&v| c.g[v][50] == 0).collect();
    if let Some(&v) = ausente.first() {
        let a = linha(&rf, &format!("rs{}", 100 + v));
        let b = linha(&r1, &format!("rs{}", 100 + v));
        t.eq(&a.n, &(b.n - 1), "a person whose file lacks the record is missing there");
    }

    // --- indexed reading of regions ---
    let gz = t.p("coorte.vcf.gz");
    externo::converter(&vcf, &gz, None).unwrap();
    indice::indexar(&gz, None).unwrap();
    let mut o = opcoes("yq");
    o.regions = vec!["chr1:1100-1250".into(), "2:5600-5700".into()];
    let ri = genotype_association(&[gz.clone()], &tab, &o).unwrap();
    let rs = genotype_association(&[vcf.clone()], &tab, &o).unwrap();
    t.eq(&ri.reading, &"indexed".to_string(), "regions with a .tbi: indexed reading");
    t.eq(&rs.reading, &"streamed".to_string(), "regions without a .tbi: streamed");
    let ids: Vec<String> = ri.rows.iter().map(|x| x.id.clone()).collect();
    t.eq(&ids, &vec!["rs101".to_string(), "rs102".into(), "rs106".into(), "rs107".into()], "regions select by POS (chr1 1100-1250, chr2 5600-5700)");
    let mesmas = ri.rows.iter().zip(&rs.rows).all(|(a, b)| a.id == b.id && (a.beta - b.beta).abs() < 1e-12 && a.n == b.n);
    t.check(mesmas && ri.rows.len() == rs.rows.len(), "indexed and streamed give the same rows");
    // several indexed files
    let mut gzs = Vec::new();
    for (k, a) in arqs.iter().enumerate() {
        let g = t.p(&format!("parte{k}.vcf.gz"));
        externo::converter(a, &g, None).unwrap();
        indice::indexar(&g, None).unwrap();
        gzs.push(g);
    }
    let rm = genotype_association(&gzs, &tab, &o).unwrap();
    t.eq(&rm.reading, &"indexed".to_string(), "several indexed files");
    let mesmas = ri.rows.iter().zip(&rm.rows).all(|(a, b)| a.id == b.id && (a.beta - b.beta).abs() < 1e-12 && a.n == b.n);
    t.check(mesmas && ri.rows.len() == rm.rows.len(), "several indexed files: same rows");

    // --- IDs ---
    let mut o = opcoes("yq");
    o.ids = vec!["rs103".into(), "RS104".into(), "chr2:5500".into(), "1:1000:A:G".into(), "chr2-5700-A-G".into(), "rs999".into(), "chr1:1000:A:T".into()];
    let r = genotype_association(&[vcf.clone()], &tab, &o).unwrap();
    let mut ids: Vec<String> = r.rows.iter().map(|x| x.id.clone()).collect();
    ids.sort();
    t.eq(&ids, &vec!["rs100".to_string(), "rs103".into(), "rs104".into(), "rs105".into(), "rs107".into()], "IDs: rsID (any case), chr:pos, chr:pos:ref:alt, chr-pos-ref-alt");
    t.eq(&r.not_found, &vec!["rs999".to_string(), "chr1:1000:A:T".into()], "IDs not found are listed");
    let e = genotype_association(&[vcf.clone()], &tab, &AssocOptions { ids: vec!["chr1:100-200".into()], ..opcoes("yq") });
    t.erro_contem(e, "looks like a region", "a region given as an ID");
}

fn testa_multialelico_e_casos(t: &mut T) {
    // 12 people, by hand
    let nomes: Vec<String> = (1..=12).map(|i| format!("S{i}")).collect();
    let gts1 = ["0/0", "0/1", "0/2", "1/2", "2/2", "1/1", "./.", "0/1", "0|2", "1|1", "0/0", "2/1"];
    let v = vec![
        Variante { chrom: "chr3".into(), pos: 500, id: "rsM".into(), refe: "C".into(), alts: vec!["T".into(), "G".into()], gts: gts1.iter().map(|s| s.to_string()).collect() },
        Variante { chrom: "chr3".into(), pos: 600, id: ".".into(), refe: "A".into(), alts: vec!["<DEL>".into()], gts: vec!["0/1".into(); 12] },
        Variante { chrom: "chrX".into(), pos: 700, id: "rsX".into(), refe: "A".into(), alts: vec!["G".into()], gts: ["0", "1", "0/1", "1", "0", "0/0", "1/1", "1", "0", "0/1", "1", "0"].iter().map(|s| s.to_string()).collect() },
    ];
    let vcf = t.p("multi.vcf");
    escreve_vcf(&vcf, &nomes, &v, None);
    let mut s = String::from("id,y,grp\n");
    let ys = [1.2, 0.4, 2.2, 1.9, 3.1, 0.5, 1.0, 1.4, 2.6, 0.9, 1.1, 2.8];
    for (i, n) in nomes.iter().enumerate() {
        s.push_str(&format!("{n},{},{}\n", ys[i], if i % 2 == 0 { "a" } else { "b" }));
    }
    let tab = t.p("multi.csv");
    escreve_texto(&tab, &s);
    let o = AssocOptions { outcome: "y".into(), covariates: vec!["grp".into()], min_maf: 0.0, min_call_rate: 0.0, ..Default::default() };
    let r = genotype_association(&[vcf.clone()], &tab, &o).unwrap();
    t.eq(&r.rows.len(), &3, "multiallelic: one row per ALT (T, G) + chrX; <DEL> not tested");
    t.eq(&r.symbolic_alleles, &1, "symbolic allele counted");
    // copies of T: 0,1,0,1,0,2,-,1,0,2,0,1 ; of G: 0,0,1,1,2,0,-,0,1,0,0,1
    let dt = [0, 1, 0, 1, 0, 2, -1, 1, 0, 2, 0, 1];
    let dg = [0, 0, 1, 1, 2, 0, -1, 0, 1, 0, 0, 1];
    for (alt, d) in [("T", dt), ("G", dg)] {
        let row = r.rows.iter().find(|x| x.alt == alt).unwrap();
        t.eq(&row.id, &"rsM".to_string(), &format!("multiallelic {alt}: the record's ID"));
        t.eq(&row.n, &11, &format!("multiallelic {alt}: n"));
        let soma: i32 = d.iter().filter(|x| **x >= 0).sum();
        t.perto(row.alt_freq, soma as f64 / 22.0, 1e-12, &format!("multiallelic {alt}: frequency of that ALT"));
        let mut x = Vec::new();
        let mut y = Vec::new();
        for i in 0..12 {
            if d[i] >= 0 {
                x.push(vec![1.0, (i % 2 == 1) as u8 as f64, d[i] as f64]);
                y.push(ys[i]);
            }
        }
        let (b, se, ..) = ols(&x, &y);
        t.perto(row.beta, b[2], 1e-10, &format!("multiallelic {alt}: beta of that ALT's dosage"));
        t.perto(row.se, se[2], 1e-10, &format!("multiallelic {alt}: se"));
        t.eq(&row.carriers, &d.iter().filter(|x| **x >= 1).count(), &format!("multiallelic {alt}: carriers"));
    }
    // haploid calls: 0/1 copies
    let row = r.rows.iter().find(|x| x.id == "rsX").unwrap();
    let dx = [0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 2.0, 1.0, 0.0, 1.0, 1.0, 0.0];
    let ploidia = [1.0, 1.0, 2.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 2.0, 1.0, 1.0];
    t.perto(row.alt_freq, dx.iter().sum::<f64>() / ploidia.iter().sum::<f64>(), 1e-12, "haploid calls: frequency over alleles called");
    t.eq(&r.haploid_calls, &8, "haploid calls counted");
    t.check(r.warnings.iter().any(|w| w.contains("haploid")), "haploid warning");

    // --- errors, each with a useful message ---
    let e = genotype_association(&[vcf.clone()], &tab, &AssocOptions { outcome: "yy".into(), ..o.clone() });
    t.erro_contem(e, "has no column \"yy\"", "missing outcome column");
    let e = genotype_association(&[vcf.clone(), vcf.clone()], &tab, &o);
    t.erro_contem(e, "each person must be in one VCF", "a sample in two VCFs");
    let mut s2 = String::from("id,y,a,b,c\n");
    for (i, n) in nomes.iter().enumerate() {
        let a = (i as f64) / 20.0;
        let b = 0.3 + (i % 3) as f64 / 10.0;
        s2.push_str(&format!("{n},{},{a},{b},{}\n", ys[i], 1.0 - a - b));
    }
    let tab2 = t.p("colinear.csv");
    escreve_texto(&tab2, &s2);
    let e = genotype_association(&[vcf.clone()], &tab2, &AssocOptions { covariates: vec!["a".into(), "b".into(), "c".into()], ..o.clone() });
    t.erro_contem(e, "proportions that sum to 1", "proportions summing to 1 are refused, by name");
    let tab3 = t.p("tres.csv");
    escreve_texto(&tab3, &nomes.iter().enumerate().map(|(i, n)| format!("{n},{}\n", ["x", "y", "z"][i % 3])).fold("id,o\n".to_string(), |a, b| a + &b));
    let e = genotype_association(&[vcf.clone()], &tab3, &AssocOptions { outcome: "o".into(), ..Default::default() });
    t.erro_contem(e, "3 distinct values", "a text outcome with three values");
    let tab4 = t.p("plink.txt");
    escreve_texto(&tab4, &nomes.iter().enumerate().map(|(i, n)| format!("{n} {n} {}\n", [1, 2, -9][i % 3])).fold("FID IID PHENO\n".to_string(), |a, b| a + &b));
    let e = genotype_association(&[vcf.clone()], &tab4, &AssocOptions { outcome: "PHENO".into(), ..Default::default() });
    t.erro_contem(e, "-9 means missing", "PLINK -9 hint");
    let tab5 = t.p("outros.csv");
    escreve_texto(&tab5, "id,y\nZ1,1\nZ2,2\n");
    let e = genotype_association(&[vcf.clone()], &tab5, &AssocOptions { outcome: "y".into(), ..Default::default() });
    t.erro_contem(e, "No sample of the VCF", "no sample in common");
    let e = genotype_association(&[vcf.clone()], &tab, &AssocOptions { output: Some(vcf.clone()), ..o.clone() });
    t.erro_contem(e, "would destroy the input", "output = input refused");
    let e = genotype_association(&[vcf.clone()], &tab, &AssocOptions { case_value: Some("2".into()), outcome_type: OutcomeType::Binary, ..o.clone() });
    t.erro_contem(e, "declared binary", "a quantitative outcome declared binary");
    // unsorted second file
    let desord = t.p("desord.vcf");
    let mut vv = v.clone();
    vv.swap(0, 2);
    escreve_vcf(&desord, &nomes, &vv, Some(&[0, 1, 2, 3, 4, 5]));
    let outro = t.p("multi_b.vcf");
    escreve_vcf(&outro, &nomes, &v, Some(&[6, 7, 8, 9, 10, 11]));
    let e = genotype_association(&[desord.clone(), outro.clone()], &tab, &o);
    t.erro_contem(e, "sort_genomic", "several VCFs: an unsorted file is refused with the fix");
    // plink-style table matched by IID; text binary outcome with a positive word
    let mut s6 = String::from("FID IID status sex\n");
    for (i, n) in nomes.iter().enumerate() {
        s6.push_str(&format!("F{i} {n} {} {}\n", if [0, 2, 3, 4, 8, 11].contains(&i) { "case" } else { "control" }, i % 2));
    }
    let tab6 = t.p("plink6.txt");
    escreve_texto(&tab6, &s6);
    let r = genotype_association(&[vcf.clone()], &tab6, &AssocOptions { outcome: "status".into(), min_maf: 0.0, min_call_rate: 0.0, ..Default::default() }).unwrap();
    t.eq(&r.event_value, &Some("case".to_string()), "text outcome: \"case\" is the event");
    t.eq(&r.analysed, &12, "FID/IID table: matched by IID");
    t.check(r.warnings.iter().any(|w| w.contains("IID")), "IID rule stated");
    t.check(r.warnings.iter().any(|w| w.contains("rarer outcome")), "small-sample warning for 6 events");

    // pass_only: a failed record gives missing genotypes
    let pf = t.p("filtro.vcf");
    let txt = std::fs::read_to_string(&vcf).unwrap().replace("chrX\t700\trsX\tA\tG\t.\tPASS", "chrX\t700\trsX\tA\tG\t.\tLowQual");
    escreve_texto(&pf, &txt);
    let r = genotype_association(&[pf.clone()], &tab, &AssocOptions { pass_only: true, ..o.clone() }).unwrap();
    t.check(!r.rows.iter().any(|x| x.id == "rsX"), "pass_only: a LowQual record is not tested (genotypes missing)");
}

fn testa_separacao(t: &mut T) {
    // 40 people; the 5 carriers of rsS all have the event
    let n = 40;
    let nomes: Vec<String> = (0..n).map(|i| format!("Q{i}")).collect();
    let mut rng = Rng::new(99);
    let y: Vec<u8> = (0..n).map(|i| if i % 8 == 0 { 1 } else { (rng.unif() < 0.35) as u8 }).collect();
    let g: Vec<&str> = (0..n).map(|i| if i % 8 == 0 { "0/1" } else { "0/0" }).collect();
    let g2: Vec<String> = (0..n).map(|_| ["0/0", "0/1", "1/1"][genotipo(&mut rng, 0.4) as usize].to_string()).collect();
    let v = vec![
        Variante { chrom: "chr5".into(), pos: 10, id: "rsS".into(), refe: "G".into(), alts: vec!["A".into()], gts: g.iter().map(|s| s.to_string()).collect() },
        Variante { chrom: "chr5".into(), pos: 20, id: "rsN".into(), refe: "G".into(), alts: vec!["A".into()], gts: g2 },
    ];
    let vcf = t.p("sep.vcf");
    escreve_vcf(&vcf, &nomes, &v, None);
    let tab = t.p("sep.tsv");
    let ages: Vec<f64> = (0..n).map(|_| (40.0 + 10.0 * rng.normal()).round()).collect();
    escreve_texto(&tab, &(0..n).map(|i| format!("{}\t{}\t{}\n", nomes[i], y[i], ages[i])).fold("s\ty\tage\n".to_string(), |a, b| a + &b));
    let o = AssocOptions { outcome: "y".into(), covariates: vec!["age".into()], min_maf: 0.0, ..Default::default() };
    let r = genotype_association(&[vcf.clone()], &tab, &o).unwrap();
    let s = linha(&r, "rsS");
    t.eq(&s.method, &"firth".to_string(), "separation: Firth used (auto)");
    t.check(s.note.contains("separation") || s.note.contains("did not converge"), &format!("separation noted ({})", s.note));
    t.check(s.beta.is_finite() && s.ci_low.is_finite() && s.ci_high.is_finite() && s.ci_method == "profile", "Firth: finite estimate and profile CI");
    t.eq(&linha(&r, "rsN").method, &"logistic".to_string(), "no separation: plain logistic");
    // the Firth row equals the direct fit
    let mut x = Vec::new();
    for i in 0..n {
        x.push(vec![1.0, ages[i], if i % 8 == 0 { 1.0 } else { 0.0 }]);
    }
    let yf: Vec<f64> = y.iter().map(|v| *v as f64).collect();
    let ff = firth_fit(&cols_de(&x), &yf, None, None).unwrap();
    t.perto(s.beta, ff.beta[2], 1e-10, "Firth row = direct fit");
    let rn = genotype_association(&[vcf.clone()], &tab, &AssocOptions { firth: FirthPolicy::Never, ..o.clone() }).unwrap();
    let sn = linha(&rn, "rsS");
    t.check(sn.method == "logistic" && sn.beta.is_nan() && sn.p_lrt < 0.05, &format!("firth=never: no estimate, only the LRT p ({})", sn.p_lrt));
    let ra = genotype_association(&[vcf.clone()], &tab, &AssocOptions { firth: FirthPolicy::Always, ..o.clone() }).unwrap();
    t.check(ra.rows.iter().all(|x| x.method == "firth"), "firth=always");
}

/// A null data set for the permutation tests. `confundido`: the genotype of
/// each variant follows a covariate c that changes the outcome (confounding;
/// the null holds given c); otherwise genotypes are unrelated to c.
fn dados_nulos(seed: u64, n: usize, nv: usize, ld: bool, dir: &T, binario: bool, confundido: bool) -> (String, String) {
    let fr = |c: f64| if confundido { (0.3 + 0.15 * c).clamp(0.05, 0.95) } else { 0.3 };
    let mut rng = Rng::new(seed);
    let nomes: Vec<String> = (0..n).map(|i| format!("N{i}")).collect();
    let c: Vec<f64> = (0..n).map(|_| (rng.normal() * 1000.0).round() / 1000.0).collect();
    let base: Vec<i32> = (0..n).map(|i| genotipo(&mut rng, fr(c[i]))).collect();
    let mut vars = Vec::new();
    for v in 0..nv {
        let gts: Vec<String> = (0..n)
            .map(|i| {
                let g = if ld && rng.unif() < 0.85 { base[i] } else { genotipo(&mut rng, fr(c[i])) };
                ["0/0", "0/1", "1/1"][g as usize].to_string()
            })
            .collect();
        vars.push(Variante { chrom: "chr7".into(), pos: 100 + v as u64, id: format!("v{v}"), refe: "C".into(), alts: vec!["T".into()], gts });
    }
    let vcf = dir.p(&format!("nulo_{seed}.vcf"));
    escreve_vcf(&vcf, &nomes, &vars, None);
    let mut s = String::from("s\ty\tc\n");
    for i in 0..n {
        let y = if binario {
            let e = -0.3 + 1.0 * c[i];
            format!("{}", (rng.unif() < 1.0 / (1.0 + (-e).exp())) as u8)
        } else {
            format!("{:.4}", 0.8 * c[i] + rng.normal())
        };
        s.push_str(&format!("{}\t{y}\t{}\n", nomes[i], c[i]));
    }
    let tab = dir.p(&format!("nulo_{seed}.tsv"));
    escreve_texto(&tab, &s);
    (vcf, tab)
}

fn testa_permutacoes(t: &mut T) {
    // determinism and range
    let (vcf, tab) = dados_nulos(1, 60, 5, false, t, false, true);
    let o = AssocOptions { outcome: "y".into(), covariates: vec!["c".into()], permutations: 199, seed: 7, min_maf: 0.0, ..Default::default() };
    let a = genotype_association(&[vcf.clone()], &tab, &o).unwrap();
    let b = genotype_association(&[vcf.clone()], &tab, &o).unwrap();
    t.check(a.rows.iter().zip(&b.rows).all(|(x, y)| x.p_permutation == y.p_permutation && x.p_fwer == y.p_fwer), "permutations: same seed, same p");
    t.check(a.rows.iter().all(|x| x.p_permutation >= 1.0 / 200.0 && x.p_permutation <= 1.0 && x.p_fwer >= x.p_permutation - 1e-12), "p_perm in [1/(B+1), 1]; p_fwer >= p_perm");
    t.check(a.permutation_scheme.contains("Freedman-Lane"), "linear with covariates: Freedman-Lane");
    // the observed statistic is the row's t: the permutation p of a strong effect is minimal
    let c = coorte(5, 400, 4, 0.8, 1.5, 0.0);
    let vars = variantes_de(&c, 5);
    let vcf2 = t.p("forte.vcf");
    escreve_vcf(&vcf2, &c.nomes, &vars, None);
    let tab2 = t.p("forte.tsv");
    escreve_tabela_coorte(&tab2, &c, &[]);
    for (y, nome) in [("yq", "linear"), ("yb", "logistic")] {
        let mut o = opcoes(y);
        o.permutations = 299;
        let r = genotype_association(&[vcf2.clone()], &tab2, &o).unwrap();
        let x = linha(&r, "rs100");
        t.check(x.p < 1e-4 && x.p_permutation <= 2.0 / 300.0 && x.p_fwer <= 2.0 / 300.0, &format!("{nome}: a strong effect has the smallest permutation p ({}, {})", x.p_permutation, x.p_fwer));
    }

    // uniform under a confounded null: linear (Freedman-Lane), logistic
    // without covariates (exact) and with (parametric bootstrap)
    for (binario, cov, nome) in [(false, true, "linear, Freedman-Lane"), (true, false, "logistic, label permutation"), (true, true, "logistic, null-model bootstrap")] {
        let mut ps = Vec::new();
        let mut fw = 0usize;
        let mut ingenuo = 0usize;
        let mut assint = Vec::new();
        let d = 300u64;
        for s in 0..d {
            let (vcf, tab) = dados_nulos(1000 + s + if binario { 5000 } else { 0 }, 80, 8, true, t, binario, cov);
            let o = AssocOptions {
                outcome: "y".into(),
                covariates: if cov { vec!["c".into()] } else { vec![] },
                permutations: 199,
                seed: s,
                min_maf: 0.0,
                ..Default::default()
            };
            let r = genotype_association(&[vcf], &tab, &o).unwrap();
            for x in &r.rows {
                ps.push(x.p_permutation);
                assint.push(x.p);
            }
            if r.rows.iter().map(|x| x.p_fwer).fold(1.0, f64::min) <= 0.05 {
                fw += 1;
            }
            if r.rows.iter().map(|x| x.p).fold(1.0, f64::min) <= 0.05 {
                ingenuo += 1;
            }
        }
        let m = ps.len() as f64;
        for lim in [0.05, 0.1, 0.25, 0.5] {
            let f = ps.iter().filter(|p| **p <= lim + 1e-12).count() as f64 / m;
            // the variants of a data set are correlated (LD): an effective
            // sample of about 3 per data set
            let se = (lim * (1.0 - lim) / (3.0 * d as f64)).sqrt();
            t.check((f - lim).abs() <= 3.5 * se + 0.005, &format!("{nome}: P(p_perm <= {lim}) = {f:.4} (null)"));
        }
        let f = fw as f64 / d as f64;
        let fi = ingenuo as f64 / d as f64;
        let se = (0.05f64 * 0.95 / d as f64).sqrt();
        t.check((f - 0.05).abs() <= 3.5 * se, &format!("{nome}: family-wise error {f:.3} at 0.05 (min p unadjusted: {fi:.3})"));
        t.check(fi > f, &format!("{nome}: the unadjusted minimum p over 8 variants rejects more often ({fi:.3} > {f:.3})"));
        if cov {
            let fa = assint.iter().filter(|p| **p <= 0.05).count() as f64 / assint.len() as f64;
            t.check((fa - 0.05).abs() < 0.03, &format!("{nome}: the asymptotic p adjusted for the confounder is calibrated too ({fa:.3})"));
        }
    }
}

fn testa_poder(t: &mut T) {
    // analytic against simulated
    let casos: Vec<(PowerOptions, usize, f64, &str)> = vec![
        (PowerOptions { outcome: OutcomeType::Quantitative, allele_freq: 0.3, ..Default::default() }, 60, 0.35, "quantitative additive n=60"),
        (PowerOptions { outcome: OutcomeType::Quantitative, allele_freq: 0.2, covariates: 3, model: GeneticModel::Dominant, ..Default::default() }, 40, 0.6, "quantitative dominant n=40, 3 covariates"),
        (PowerOptions { outcome: OutcomeType::Quantitative, allele_freq: 0.4, model: GeneticModel::Recessive, alpha: 0.05, tests: 20, ..Default::default() }, 150, 0.6, "quantitative recessive, 20 tests"),
        (PowerOptions { outcome: OutcomeType::Binary, allele_freq: 0.25, case_fraction: 0.4, ..Default::default() }, 1000, 1.35, "binary additive n=1000"),
        (PowerOptions { outcome: OutcomeType::Binary, allele_freq: 0.15, case_fraction: 0.5, model: GeneticModel::Dominant, ..Default::default() }, 600, 1.6, "binary dominant n=600"),
        (PowerOptions { outcome: OutcomeType::Binary, allele_freq: 0.35, case_fraction: 0.3, model: GeneticModel::Recessive, tests: 10, ..Default::default() }, 1500, 2.0, "binary recessive n=1500, 10 tests"),
    ];
    for (o, n, e, nome) in casos {
        let o = PowerOptions { effect: Some(e), n: Some(n), seed: 3, ..o };
        let r = association_power(&o).unwrap();
        let reps = 4000;
        let (s, sw) = simulated_power(&o, n, e, reps);
        let alvo = if o.outcome == OutcomeType::Binary { sw } else { s };
        let se = (r.power * (1.0 - r.power) / reps as f64).sqrt();
        t.check((alvo - r.power).abs() <= 4.0 * se + 0.012, &format!("{nome}: analytic power {:.4} vs simulated {:.4} (+-{:.4})", r.power, alvo, se));
        // the people needed
        if let Some(k) = r.n_needed {
            let a = association_power(&PowerOptions { n: Some(k), ..o.clone() }).unwrap().power;
            let b = association_power(&PowerOptions { n: Some(k - 1), ..o.clone() }).map(|x| x.power).unwrap_or(0.0);
            t.check(a >= 0.8 && b < 0.8, &format!("{nome}: n_needed {k}: power {a:.4} there, {b:.4} with one fewer"));
        }
        // the detectable effect gives the target power
        let d = r.detectable_effect;
        let pd = association_power(&PowerOptions { effect: Some(d), ..o.clone() }).unwrap().power;
        t.perto_piso(pd, 0.8, 1e-6, 1.0, &format!("{nome}: detectable effect {d:.4} has power 0.8"));
    }
    // a known value: Wald power of the additive logistic test by the expected
    // information computed here independently (enumerating the 6 cells)
    {
        let (p, or_, fi, n) = (0.3f64, 1.5f64, 0.4f64, 500.0);
        let q = 1.0 - p;
        let d0 = [q * q, 2.0 * p * q, p * p];
        let z: f64 = (0..3).map(|x| d0[x] * or_.powi(x as i32)).sum();
        let d1: Vec<f64> = (0..3).map(|x| d0[x] * or_.powi(x as i32) / z).collect();
        let mut inf = [0.0; 4];
        for x in 0..3 {
            for (caso, w) in [(1.0, fi * d1[x]), (0.0, (1.0 - fi) * d0[x])] {
                let _ = caso;
                let px = fi * d1[x] + (1.0 - fi) * d0[x];
                let pi = fi * d1[x] / px;
                let xv = x as f64;
                inf[0] += w * pi * (1.0 - pi);
                inf[1] += w * pi * (1.0 - pi) * xv;
                inf[3] += w * pi * (1.0 - pi) * xv * xv;
            }
        }
        let var = inf[0] / (inf[0] * inf[3] - inf[1] * inf[1]);
        let lam = n * or_.ln().powi(2) / var;
        let want = normal_upper(Z975 - lam.sqrt()) + normal_upper(Z975 + lam.sqrt());
        let r = association_power(&PowerOptions { outcome: OutcomeType::Binary, allele_freq: p, effect: Some(or_), n: Some(500), case_fraction: fi, ..Default::default() }).unwrap();
        t.perto(r.power, want, 1e-12, "binary power = expected information enumerated by hand");
    }
    // five patients: nothing can be detected
    let r = association_power(&PowerOptions { outcome: OutcomeType::Binary, allele_freq: 0.25, effect: Some(2.0), n: Some(5), case_fraction: 0.4, simulations: 2000, ..Default::default() }).unwrap();
    t.check(r.power < 0.1, &format!("5 patients, OR 2: power {:.4}", r.power));
    t.check(r.detectable_effect.is_nan() || r.detectable_effect > 10.0, "5 patients: no realistic odds ratio is detectable");
    t.check(r.text().contains("cannot be expected to detect"), "5 patients: the text says so");
    t.check(r.simulated_power < 0.1, &format!("5 patients: simulated power {:.4}", r.simulated_power));
    t.check(r.n_needed.map_or(false, |k| k > 100 && k < 400), &format!("OR 2 at 25%: people needed {:?}", r.n_needed));
    let js: serde_json::Value = serde_json::from_str(&r.json()).expect("power JSON");
    t.check(js["curve"]["n"].as_array().map_or(false, |a| a.len() > 10), "power JSON: curve");
    // errors
    t.check(association_power(&PowerOptions { allele_freq: 1.2, n: Some(10), ..Default::default() }).is_err(), "power: frequency out of range");
    t.check(association_power(&PowerOptions { outcome: OutcomeType::Quantitative, allele_freq: 0.3, effect: Some(3.0), n: Some(10), ..Default::default() }).is_err(), "power: beta explaining all the variance is refused");
    t.check(association_power(&PowerOptions::default()).is_err(), "power: needs an effect or a number of people");
}

// ---------------------------------------------------------------- statsmodels references

fn testa_referencias(t: &mut T) {
    // deterministic data sets written by `dump`; values from statsmodels
    // (scratchpad/agent_assoc/sm_refs.py); empty until filled.
    let c = coorte(31, 200, 6, 0.3, 0.7, 0.05);
    let vars = variantes_de(&c, 31);
    let vcf = t.p("ref.vcf");
    escreve_vcf(&vcf, &c.nomes, &vars, None);
    let tab = t.p("ref.tsv");
    escreve_tabela_coorte(&tab, &c, &[(7, "pc1"), (8, "center")]);
    let rq = genotype_association(&[vcf.clone()], &tab, &opcoes("yq")).unwrap();
    let rb = genotype_association(&[vcf.clone()], &tab, &opcoes("yb")).unwrap();
    for (id, b, se, p) in REF_OLS {
        let r = linha(&rq, id);
        t.perto(r.beta, *b, 1e-6, &format!("statsmodels OLS {id} beta"));
        t.perto(r.se, *se, 1e-6, &format!("statsmodels OLS {id} se"));
        t.perto(r.p, *p, 1e-6, &format!("statsmodels OLS {id} p"));
    }
    for (id, b, se, pw, plr) in REF_LOGIT {
        let r = linha(&rb, id);
        t.perto(r.beta, *b, 1e-6, &format!("statsmodels Logit {id} beta"));
        t.perto(r.se, *se, 1e-6, &format!("statsmodels Logit {id} se"));
        t.perto(r.p_wald, *pw, 1e-6, &format!("statsmodels Logit {id} Wald p"));
        t.perto(r.p_lrt, *plr, 1e-6, &format!("statsmodels Logit {id} LRT p"));
    }
}

// statsmodels 0.15.0 (OLS; Logit with method="newton", tol=1e-14; LRT
// against the covariates-only Logit), from scratchpad/agent_assoc/refs.py.
const REF_OLS: &[(&str, f64, f64, f64)] = &[
    ("rs100", 0.17346554985341894, 0.18117512869777944, 0.3396188897420376),
    ("rs101", 0.11771250246258225, 0.1469587394540674, 0.42418753724770025),
    ("rs102", 0.04503521714108874, 0.13553940027675287, 0.7400745570011064),
    ("rs103", -0.2838730844748523, 0.11354208748369102, 0.01328826303099758),
    ("rs104", -0.00873779080916505, 0.11365616991290978, 0.9388045914867669),
    ("rs105", 0.00695356032348869, 0.10936138654039826, 0.9493714094861708),
];
const REF_LOGIT: &[(&str, f64, f64, f64, f64)] = &[
    ("rs100", 0.45837616603697323, 0.39490631327895437, 0.24575526025751904, 0.2402512728737995),
    ("rs101", -0.08514134924008894, 0.3053583207848121, 0.7803795916172771, 0.7803724222655626),
    ("rs102", 0.18439321323063262, 0.28644684463641745, 0.5197532815261052, 0.5194135788322612),
    ("rs103", -0.25015274927095366, 0.25001348262238887, 0.3170410098673103, 0.314754936552725),
    ("rs104", -0.16498087038199774, 0.24370916938115464, 0.49843260665189226, 0.49734779167096843),
    ("rs105", 0.13947801740437477, 0.23375477197420133, 0.5507175626321476, 0.5500738574840894),
];

fn dump(dir: &str) {
    std::fs::create_dir_all(dir).unwrap();
    let t = T { ok: 0, falhas: 0, dir: PathBuf::from(dir) };
    let c = coorte(31, 200, 6, 0.3, 0.7, 0.05);
    let vars = variantes_de(&c, 31);
    escreve_vcf(&t.p("ref.vcf"), &c.nomes, &vars, None);
    escreve_tabela_coorte(&t.p("ref.tsv"), &c, &[(7, "pc1"), (8, "center")]);
    let (x, y) = dados_firth();
    let mut s = String::from("y,a,b,g\n");
    for (r, yi) in x.iter().zip(&y) {
        s.push_str(&format!("{},{},{},{}\n", yi, r[1], r[2], r[3]));
    }
    escreve_texto(&t.p("firth.csv"), &s);
    let c = coorte(5, 160, 4, 0.8, 1.5, 0.0);
    escreve_vcf(&t.p("forte.vcf"), &c.nomes, &variantes_de(&c, 5), None);
    escreve_tabela_coorte(&t.p("forte.tsv"), &c, &[]);
    println!("written to {dir}");
}

// ---------------------------------------------------------------- CLI modes

fn chaves(args: &[String]) -> std::collections::HashMap<String, String> {
    args.iter().filter_map(|a| a.split_once('=')).map(|(k, v)| (k.to_string(), v.to_string())).collect()
}

fn lista(s: Option<&String>) -> Vec<String> {
    s.map(|v| v.split(',').filter(|x| !x.is_empty()).map(str::to_string).collect()).unwrap_or_default()
}

fn modo_run(args: &[String]) -> Result<(), String> {
    let k = chaves(args);
    let o = AssocOptions {
        outcome: k.get("outcome").cloned().unwrap_or_default(),
        covariates: lista(k.get("covariates")),
        categorical: lista(k.get("categorical")),
        sample_column: k.get("sample_column").cloned(),
        outcome_type: OutcomeType::from_text(k.get("type").map_or("auto", |s| s))?,
        case_value: k.get("case_value").cloned(),
        model: GeneticModel::from_text(k.get("model").map_or("additive", |s| s))?,
        ids: lista(k.get("ids")),
        regions: lista(k.get("regions")),
        min_maf: k.get("min_maf").map_or(Ok(0.01), |s| s.parse()).map_err(|e| format!("{e}"))?,
        min_call_rate: k.get("min_call_rate").map_or(Ok(0.9), |s| s.parse()).map_err(|e| format!("{e}"))?,
        firth: FirthPolicy::from_text(k.get("firth").map_or("auto", |s| s))?,
        permutations: k.get("permutations").map_or(Ok(0), |s| s.parse()).map_err(|e| format!("{e}"))?,
        seed: k.get("seed").map_or(Ok(1), |s| s.parse()).map_err(|e| format!("{e}"))?,
        absent_is_ref: k.get("absent_is_ref").map_or(false, |s| s == "1" || s == "true"),
        pass_only: k.get("pass_only").map_or(false, |s| s == "1" || s == "true"),
        output: k.get("output").cloned(),
    };
    let r = genotype_association(&lista(k.get("vcfs")), k.get("pheno").map_or("", |s| s), &o)?;
    match k.get("json") {
        Some(p) => std::fs::write(p, r.json()).map_err(|e| e.to_string())?,
        None => println!("{}", r.json()),
    }
    eprint!("{}", r.text());
    Ok(())
}

fn modo_power(args: &[String]) -> Result<(), String> {
    let k = chaves(args);
    let num = |c: &str, d: f64| -> Result<f64, String> { k.get(c).map_or(Ok(d), |s| s.parse::<f64>().map_err(|e| format!("{c}: {e}"))) };
    let o = PowerOptions {
        outcome: OutcomeType::from_text(k.get("type").map_or("binary", |s| s))?,
        allele_freq: num("freq", 0.2)?,
        effect: k.get("effect").map(|s| s.parse().unwrap()),
        n: k.get("n").map(|s| s.parse().unwrap()),
        case_fraction: num("case_fraction", 0.5)?,
        alpha: num("alpha", 0.05)?,
        tests: num("tests", 1.0)? as usize,
        model: GeneticModel::from_text(k.get("model").map_or("additive", |s| s))?,
        covariates: num("covariates", 0.0)? as usize,
        target_power: num("target", 0.8)?,
        simulations: num("simulations", 0.0)? as usize,
        seed: num("seed", 1.0)? as u64,
    };
    let r = association_power(&o)?;
    println!("{}", r.json());
    eprint!("{}", r.text());
    Ok(())
}

/// `fit <csv> <logit|firth>`: y in the first column, the intercept added.
fn modo_fit(csv: &str, metodo: &str) -> Result<(), String> {
    let txt = std::fs::read_to_string(csv).map_err(|e| e.to_string())?;
    let mut x = Vec::new();
    let mut y = Vec::new();
    for l in txt.lines().skip(1).filter(|l| !l.trim().is_empty()) {
        let v: Vec<f64> = l.split(',').map(|s| s.trim().parse::<f64>().unwrap()).collect();
        y.push(v[0]);
        let mut r = vec![1.0];
        r.extend_from_slice(&v[1..]);
        x.push(r);
    }
    let cols = cols_de(&x);
    let p = cols.len();
    let fmt = |v: &[f64]| format!("[{}]", v.iter().map(|x| format!("{x:e}")).collect::<Vec<_>>().join(","));
    if metodo == "logit" {
        let f = logistic_fit(&cols, &y, None);
        println!("{{\"beta\":{},\"se\":{},\"loglik\":{:e},\"converged\":{},\"separated\":{}}}", fmt(&f.beta), fmt(&f.se), f.loglik, f.converged, f.separated);
    } else {
        let f = firth_fit(&cols, &y, None, None).ok_or("firth failed")?;
        let mut pl = Vec::new();
        let mut lo = Vec::new();
        let mut hi = Vec::new();
        for j in 0..p {
            pl.push(firth_lrt(&cols, &y, &f, j).1);
            let (a, b) = firth_profile_ci(&cols, &y, &f, j, CHI2_1_95);
            lo.push(a);
            hi.push(b);
        }
        println!(
            "{{\"beta\":{},\"se\":{},\"loglik_pen\":{:e},\"converged\":{},\"p_lrt\":{},\"ci_low\":{},\"ci_high\":{}}}",
            fmt(&f.beta),
            fmt(&f.se),
            f.loglik_pen,
            f.converged,
            fmt(&pl),
            fmt(&lo),
            fmt(&hi)
        );
    }
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let r = match args.get(1).map(String::as_str) {
        Some("run") => Some(modo_run(&args[2..])),
        Some("power") => Some(modo_power(&args[2..])),
        Some("fit") => Some(modo_fit(&args[2], &args[3])),
        Some("dump") => {
            dump(&args[2]);
            Some(Ok(()))
        }
        _ => None,
    };
    if let Some(r) = r {
        if let Err(e) = r {
            eprintln!("error: {e}");
            std::process::exit(1);
        }
        return;
    }
    let dir = std::env::temp_dir().join(format!("mmg_associacao_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut t = T { ok: 0, falhas: 0, dir: dir.clone() };
    let inicio = std::time::Instant::now();
    testa_funcoes(&mut t);
    testa_logistica(&mut t);
    testa_ponta_a_ponta(&mut t);
    testa_multialelico_e_casos(&mut t);
    testa_separacao(&mut t);
    testa_referencias(&mut t);
    testa_poder(&mut t);
    testa_permutacoes(&mut t);
    let _ = std::fs::remove_dir_all(&dir);
    println!("associacao tests: {} passed, {} failed ({:.1} s)", t.ok, t.falhas, inicio.elapsed().as_secs_f64());
    if t.falhas > 0 {
        std::process::exit(1);
    }
}
