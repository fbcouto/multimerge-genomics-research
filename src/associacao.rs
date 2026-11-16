// =====================================================================
//  ASSOCIATION BETWEEN GENOTYPES AND AN OUTCOME IN A COHORT, AND ITS POWER
//
//  Two tools:
//
//   - `genotype_association`  one test per ALT allele of each selected
//                             variant, against an outcome from a phenotype
//                             table, adjusted for covariates:
//                               * quantitative outcome: linear regression
//                                 (ordinary least squares, t test of the
//                                 genotype coefficient);
//                               * binary outcome: logistic regression (IRLS;
//                                 Wald and likelihood-ratio tests), with
//                                 Firth's penalized likelihood when the plain
//                                 fit separates or fails (or always);
//                             additive, dominant or recessive coding;
//                             Benjamini-Hochberg FDR across the variants;
//                             optional permutations: an empirical p per
//                             variant and a family-wise (single-step max-T /
//                             Westfall-Young) adjusted p;
//   - `association_power`     the power of that test for a number of people,
//                             an allele frequency and an effect, the number of
//                             people needed for a target power and the
//                             smallest effect detectable with the people one
//                             has. This is what tells a group with 5 patients
//                             that nothing can be detected.
//
//  GENOTYPES come from one multi-sample VCF or several VCFs (one or more
//  samples each; text, .gz -- gzip or BGZF --, .lz4, .bam is not a VCF).
//  People are matched by sample name between the VCFs and the phenotype
//  table. Each ALT allele of a record is one test, with the dosage of THAT
//  allele (a 1/2 genotype has one copy of ALT 1 and one of ALT 2), as plink2
//  does with split multiallelic variants.
//
//  SEVERAL VCFs are joined by position (a k-way merge; every file sorted by
//  position, chromosomes in natural order -- 1, 2, ..., 22, X, Y; `chr1` and
//  `1` are the same). At a position, a person whose file has a record with
//  the same REF that lists the allele gets the dosage from that record; a
//  record with the same REF that does not list it (another ALT) means zero
//  copies; no record with that REF means NO INFORMATION: the genotype is
//  missing, unless `absent_is_ref` -- there is no callable-region (BED) logic
//  here, and the result says so. When several files are merged, the
//  mitochondria and sequences that are not primary chromosomes are skipped
//  (their order differs between references) and counted.
//
//  SELECTION: every variant (with `min_maf` and `min_call_rate`), or a list
//  of IDs (rsID, `chr:pos`, `chr:pos:ref:alt`, `chr-pos-ref-alt`), or regions
//  (`chr:start-end`, 1-based, the variant's POS inside). IDs and regions add
//  up. With regions (or positional IDs only) and a `.tbi` next to every
//  (BGZF) file, only the indexed blocks are read.
//
//  PEOPLE: complete cases. The analysis set is everyone in a VCF and in the
//  table with the outcome and every covariate; each test then uses the
//  people of that set whose genotype is known. The report counts who was
//  dropped and why, per run and per variant.
//
//  THE TESTS
//   - Linear: by Frisch-Waugh-Lovell -- the genotype and the outcome are
//     residualized on the covariates (a QR by Gram-Schmidt with
//     re-orthogonalization), which gives exactly the OLS coefficient, its
//     standard error, t and the t-based CI of the full model, at O(n k) per
//     variant.
//   - Logistic: Newton-Raphson (IRLS) with step-halving on the
//     log-likelihood. Wald z, its CI, and the likelihood-ratio test against
//     the covariates-only model fitted on the same people. The main `p` is
//     the LRT, better behaved than the Wald test in small samples (the Wald
//     test loses power for large effects: Hauck-Donner); `p_wald` is given
//     too (it is what plink2 and R's summary() print).
//   - Separation (some genotype class with no events, or only events): the
//     MLE is infinite. It shows as a fit that does not converge, or fitted
//     probabilities within e^-25 of 0 or 1.
//   - Firth (logistf, Heinze & Schemper 2002): maximizes l(b) + 1/2 log|I(b)|
//     by the modified score U* = X'(y - p + h(1/2 - p)), h the hat values, with
//     step-halving on the penalized likelihood. p: penalized likelihood ratio
//     with the genotype fixed at 0 and the penalty of the full model (as
//     logistf); CI: profile penalized likelihood (as logistf); Wald SE from
//     I(b)^-1.
//
//  PERMUTATIONS (optional). The statistic of each replicate is the t of the
//  genotype (linear) or the score statistic of the genotype given the
//  covariates (logistic; it needs only the covariates-only fit, and stays
//  finite under separation). The same replicate is applied to every variant,
//  so the correlation between variants (LD) is kept, which the max-T needs.
//   - Linear, Freedman-Lane (1983): residuals of the covariates-only model
//     are permuted and added back to its fitted values; exact without
//     covariates, and the scheme with the best properties with them
//     (Anderson & Legendre 1999; Winkler et al. 2014).
//   - Logistic without covariates: the outcome labels are permuted (exact:
//     the number of events is kept).
//   - Logistic with covariates: there is no exact permutation that keeps the
//     covariate-outcome relation. Each replicate draws the outcome from the
//     covariates-only logistic model fitted to the data (a parametric
//     bootstrap under the null, i.e. a Monte Carlo test; Davison & Hinkley
//     1997, ch. 4) -- valid when that null model is right, and it keeps the
//     genotype-covariate relation (confounding by ancestry, for instance),
//     which permuting genotypes would break.
//   Empirical p = (r + 1) / (B + 1) (Phipson & Smyth 2010), r the replicates
//   at least as extreme. Family-wise p (single-step max-T, Westfall & Young
//   1993): the fraction of replicates whose smallest p, over all variants, is
//   at most the variant's p -- the p of the permutation statistic, so that
//   variants with different numbers of people compare on one scale.
//
//  NUMERICS: the normal (Cody's algorithm, as in R), its quantile (Wichura's
//  AS241), log-gamma (Lanczos), the regularized incomplete beta (continued
//  fraction, both tails, so tiny p-values keep their digits), Student's t,
//  and the non-central t power (by integration over the chi distribution).
//  Only std and rayon.
//
//  Software for teaching and research: an association is not causation, and
//  nothing here is a diagnosis.
// =====================================================================

use crate::bgzf;
use crate::externo::{abre_leitura, natural_cmp, Escritor};
use crate::indice::{consultar, le_regiao, le_tbi};
use crate::nomes::{canonico, json_texto, parece_cromossomo};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

// ---------------------------------------------------------------------
//  OPTIONS
// ---------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutcomeType {
    /// Two distinct values: binary; otherwise quantitative.
    Auto,
    Binary,
    Quantitative,
}

impl OutcomeType {
    pub fn from_text(s: &str) -> Result<OutcomeType, String> {
        match s.trim().to_ascii_lowercase().as_str() {
            "" | "auto" => Ok(OutcomeType::Auto),
            "binary" | "logistic" | "case-control" | "case_control" | "dichotomous" => Ok(OutcomeType::Binary),
            "quantitative" | "linear" | "continuous" | "numeric" => Ok(OutcomeType::Quantitative),
            _ => Err(format!("Unknown outcome type \"{s}\": use \"auto\", \"binary\" or \"quantitative\"")),
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            OutcomeType::Auto => "auto",
            OutcomeType::Binary => "binary",
            OutcomeType::Quantitative => "quantitative",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GeneticModel {
    /// 0, 1, 2 copies of the allele.
    Additive,
    /// Carrier (1 or 2 copies) against non-carrier.
    Dominant,
    /// Two copies (or one, in a haploid call) against the rest.
    Recessive,
}

impl GeneticModel {
    pub fn from_text(s: &str) -> Result<GeneticModel, String> {
        match s.trim().to_ascii_lowercase().as_str() {
            "" | "additive" | "add" | "allelic" | "log-additive" | "multiplicative" => Ok(GeneticModel::Additive),
            "dominant" | "dom" | "carrier" => Ok(GeneticModel::Dominant),
            "recessive" | "rec" => Ok(GeneticModel::Recessive),
            _ => Err(format!("Unknown genetic model \"{s}\": use \"additive\", \"dominant\" or \"recessive\"")),
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            GeneticModel::Additive => "additive",
            GeneticModel::Dominant => "dominant",
            GeneticModel::Recessive => "recessive",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FirthPolicy {
    /// Plain logistic regression; Firth when it separates or fails.
    Auto,
    Always,
    Never,
}

impl FirthPolicy {
    pub fn from_text(s: &str) -> Result<FirthPolicy, String> {
        match s.trim().to_ascii_lowercase().as_str() {
            "" | "auto" | "fallback" | "firth-fallback" => Ok(FirthPolicy::Auto),
            "always" | "yes" | "true" | "firth" => Ok(FirthPolicy::Always),
            "never" | "no" | "false" | "off" => Ok(FirthPolicy::Never),
            _ => Err(format!("Unknown firth policy \"{s}\": use \"auto\" (when the plain fit separates or fails), \"always\" or \"never\"")),
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            FirthPolicy::Auto => "auto",
            FirthPolicy::Always => "always",
            FirthPolicy::Never => "never",
        }
    }
}

#[derive(Clone, Debug)]
pub struct AssocOptions {
    /// The outcome column of the phenotype table.
    pub outcome: String,
    /// Covariate columns. Numbers are used as numbers; anything else is a
    /// category (indicators, the first level in sorted order is the reference).
    pub covariates: Vec<String>,
    /// Covariates to treat as categories even when their values are numbers
    /// (batch 1, 2, 3).
    pub categorical: Vec<String>,
    /// The sample column; `None`: the first column (or `IID` in a plink file).
    pub sample_column: Option<String>,
    pub outcome_type: OutcomeType,
    /// The value of a binary outcome that counts as the event (case).
    pub case_value: Option<String>,
    pub model: GeneticModel,
    /// rsIDs, `chr:pos`, `chr:pos:ref:alt` or `chr-pos-ref-alt`.
    pub ids: Vec<String>,
    /// `chr:start-end` (1-based, end included) or a whole `chr`.
    pub regions: Vec<String>,
    /// Minimum minor allele frequency among the people analysed. Default 0.01.
    pub min_maf: f64,
    /// Minimum fraction of the people analysed with a genotype. Default 0.9.
    pub min_call_rate: f64,
    pub firth: FirthPolicy,
    /// Replicates for the empirical and family-wise p-values (0: none).
    pub permutations: usize,
    pub seed: u64,
    /// A person whose VCF has no record at a position is homozygous for the
    /// reference there (single-sample VCFs that list only variant sites).
    /// Default false: missing.
    pub absent_is_ref: bool,
    /// Only records with FILTER PASS or `.`; the others make the genotype missing.
    pub pass_only: bool,
    /// The results table (TSV; `.gz` gives BGZF).
    pub output: Option<String>,
}

impl Default for AssocOptions {
    fn default() -> Self {
        AssocOptions {
            outcome: String::new(),
            covariates: Vec::new(),
            categorical: Vec::new(),
            sample_column: None,
            outcome_type: OutcomeType::Auto,
            case_value: None,
            model: GeneticModel::Additive,
            ids: Vec::new(),
            regions: Vec::new(),
            min_maf: 0.01,
            min_call_rate: 0.9,
            firth: FirthPolicy::Auto,
            permutations: 0,
            seed: 1,
            absent_is_ref: false,
            pass_only: false,
            output: None,
        }
    }
}

// ---------------------------------------------------------------------
//  SMALL HELPERS: NUMBERS AS TEXT, JSON
// ---------------------------------------------------------------------

/// A number for a table: shortest exact form, scientific notation for very
/// small or very large values, `NA` when missing.
pub fn num(x: f64) -> String {
    if x.is_nan() {
        return "NA".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "Inf".into() } else { "-Inf".into() };
    }
    let a = x.abs();
    if a == 0.0 || (1e-4..1e15).contains(&a) {
        format!("{x}")
    } else {
        format!("{x:e}")
    }
}

fn jnum(x: f64) -> String {
    if x.is_finite() {
        num(x)
    } else {
        "null".into()
    }
}

fn jtextos(v: &[String]) -> String {
    format!("[{}]", v.iter().map(|s| json_texto(s)).collect::<Vec<_>>().join(","))
}

fn jopt(s: Option<&str>) -> String {
    s.map(json_texto).unwrap_or_else(|| "null".into())
}

/// A p-value for a sentence.
fn fmt_p(p: f64) -> String {
    if p.is_nan() {
        "NA".into()
    } else if p < 1e-4 {
        format!("{p:.2e}")
    } else {
        format!("{p:.4}")
    }
}

fn fmt_g(x: f64, casas: usize) -> String {
    if !x.is_finite() {
        return if x.is_nan() { "NA".into() } else if x > 0.0 { "Inf".into() } else { "-Inf".into() };
    }
    if x != 0.0 && (x.abs() < 1e-3 || x.abs() >= 1e6) {
        format!("{x:.casas$e}")
    } else {
        format!("{x:.casas$}")
    }
}

/// 12345 -> "12,345".
fn milhar(n: usize) -> String {
    let s = n.to_string();
    let mut o = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            o.push(',');
        }
        o.push(c);
    }
    o
}

fn lista_curta(v: &[String], n: usize) -> String {
    let mut s: Vec<String> = v.iter().take(n).cloned().collect();
    if v.len() > n {
        s.push(format!("and {} more", v.len() - n));
    }
    s.join(", ")
}

/// Removes a half-written output when the operation ends in error.
struct Saida {
    caminho: String,
    armada: bool,
}

impl Drop for Saida {
    fn drop(&mut self) {
        if self.armada {
            let _ = std::fs::remove_file(&self.caminho);
        }
    }
}

fn mesmo_arquivo(a: &str, b: &str) -> bool {
    a == b || matches!((std::fs::canonicalize(a), std::fs::canonicalize(b)), (Ok(x), Ok(y)) if x == y)
}

// =====================================================================
//  1) SPECIAL FUNCTIONS
// =====================================================================

/// 1.96: the 97.5% quantile of the standard normal.
pub const Z975: f64 = 1.959_963_984_540_054;
/// The 95% quantile of a chi-square with 1 degree of freedom (= Z975^2).
pub const CHI2_1_95: f64 = 3.841_458_820_694_124;

/// (lower, upper) tails of the standard normal at x: Cody's algorithm, as in
/// R's pnorm (relative accuracy ~1e-16 down to 1e-300).
fn normal_caudas(x: f64) -> (f64, f64) {
    const A: [f64; 5] = [2.2352520354606839287, 161.02823106855587881, 1067.6894854603709582, 18154.981253343561249, 0.065682337918207449113];
    const B: [f64; 4] = [47.20258190468824187, 976.09855173777669322, 10260.932208618978205, 45507.789335026729956];
    const C: [f64; 9] = [
        0.39894151208813466764,
        8.8831497943883759412,
        93.506656132177855979,
        597.27027639480026226,
        2494.5375852903726711,
        6848.1904505362823326,
        11602.651437647350124,
        9842.7148383839780218,
        1.0765576773720192317e-8,
    ];
    const D: [f64; 8] = [
        22.266688044328115691,
        235.38790178262499861,
        1519.377599407554805,
        6485.558298266760755,
        18615.571640885098091,
        34900.952721145977266,
        38912.003286093271411,
        19685.429676859990727,
    ];
    const P: [f64; 6] = [0.21589853405795699, 0.1274011611602473639, 0.022235277870649807, 0.001421619193227893466, 2.9112874951168792e-5, 0.02307344176494017303];
    const Q: [f64; 5] = [1.28426009614491121, 0.468238212480865118, 0.0659881378689285515, 0.00378239633202758244, 7.29751555083966205e-5];
    const UM_RAIZ_2PI: f64 = 0.398942280401432677939946059934;
    if x.is_nan() {
        return (f64::NAN, f64::NAN);
    }
    let y = x.abs();
    let com_del = |v: f64, temp: f64| -> (f64, f64) {
        let xsq = (v * 16.0).trunc() / 16.0;
        let del = (v - xsq) * (v + xsq);
        let cum = (-xsq * xsq * 0.5).exp() * (-del * 0.5).exp() * temp;
        if x > 0.0 {
            (1.0 - cum, cum)
        } else {
            (cum, 1.0 - cum)
        }
    };
    if y <= 0.67448975 {
        let (mut xnum, mut xden) = (0.0, 0.0);
        if y > 1.1102230246251565e-16 {
            let xsq = x * x;
            xnum = A[4] * xsq;
            xden = xsq;
            for i in 0..3 {
                xnum = (xnum + A[i]) * xsq;
                xden = (xden + B[i]) * xsq;
            }
        }
        let temp = x * (xnum + A[3]) / (xden + B[3]);
        (0.5 + temp, 0.5 - temp)
    } else if y <= 5.656854249492380195206754896838 {
        let mut xnum = C[8] * y;
        let mut xden = y;
        for i in 0..7 {
            xnum = (xnum + C[i]) * y;
            xden = (xden + D[i]) * y;
        }
        let temp = (xnum + C[7]) / (xden + D[7]);
        com_del(y, temp)
    } else if y < 50.0 {
        let xsq = 1.0 / (x * x);
        let mut xnum = P[5] * xsq;
        let mut xden = xsq;
        for i in 0..4 {
            xnum = (xnum + P[i]) * xsq;
            xden = (xden + Q[i]) * xsq;
        }
        let mut temp = xsq * (xnum + P[4]) / (xden + Q[4]);
        temp = (UM_RAIZ_2PI - temp) / y;
        com_del(y, temp)
    } else if x > 0.0 {
        (1.0, 0.0)
    } else {
        (0.0, 1.0)
    }
}

/// P(Z <= x).
pub fn normal_lower(x: f64) -> f64 {
    normal_caudas(x).0
}

/// P(Z > x), accurate far in the tail.
pub fn normal_upper(x: f64) -> f64 {
    normal_caudas(x).1
}

/// Two-sided p of a z statistic.
pub fn p_two_sided_z(z: f64) -> f64 {
    if z.is_nan() {
        return f64::NAN;
    }
    (2.0 * normal_upper(z.abs())).min(1.0)
}

/// P(X > x) for a chi-square with 1 degree of freedom (= the two-sided p of
/// sqrt(x)).
pub fn chi2_1_upper(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x <= 0.0 {
        return 1.0;
    }
    (2.0 * normal_upper(x.sqrt())).min(1.0)
}

/// The p quantile of the standard normal (Wichura's AS241, PPND16; relative
/// accuracy ~1e-16). For an upper quantile use -normal_quantile(1 - q) with q
/// small, i.e. `-normal_quantile(alpha / 2)`, which keeps the digits.
pub fn normal_quantile(p: f64) -> f64 {
    if !(p > 0.0 && p < 1.0) {
        return if p == 0.0 {
            f64::NEG_INFINITY
        } else if p == 1.0 {
            f64::INFINITY
        } else {
            f64::NAN
        };
    }
    let q = p - 0.5;
    if q.abs() <= 0.425 {
        let r = 0.180625 - q * q;
        return q
            * (((((((2509.0809287301226727 * r + 33430.575583588128105) * r + 67265.770927008700853) * r + 45921.953931549871457) * r
                + 13731.693765509461125)
                * r
                + 1971.5909503065514427)
                * r
                + 133.14166789178437745)
                * r
                + 3.387132872796366608)
            / (((((((5226.495278852545925 * r + 28729.085735721942674) * r + 39307.89580009271061) * r + 21213.794301586595867) * r
                + 5394.1960214247511077)
                * r
                + 687.1870074920579083)
                * r
                + 42.313330701600911252)
                * r
                + 1.0);
    }
    let mut r = if q < 0.0 { p } else { 1.0 - p };
    r = (-r.ln()).sqrt();
    let val = if r <= 5.0 {
        r -= 1.6;
        (((((((7.7454501427834140764e-4 * r + 0.0227238449892691845833) * r + 0.24178072517745061177) * r + 1.27045825245236838258) * r
            + 3.64784832476320460504)
            * r
            + 5.7694972214606914055)
            * r
            + 4.6303378461565452959)
            * r
            + 1.42343711074968357734)
            / (((((((1.05075007164441684324e-9 * r + 5.475938084995344946e-4) * r + 0.0151986665636164571966) * r + 0.14810397642748007459)
                * r
                + 0.68976733498510000455)
                * r
                + 1.6763848301838038494)
                * r
                + 2.05319162663775882187)
                * r
                + 1.0)
    } else {
        r -= 5.0;
        (((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r + 0.0012426609473880784386) * r + 0.026532189526576123093)
            * r
            + 0.29656057182850489123)
            * r
            + 1.7848265399172913358)
            * r
            + 5.4637849111641143699)
            * r
            + 6.6579046435011037772)
            / (((((((2.04426310338993978564e-15 * r + 1.4215117583164458887e-7) * r + 1.8463183175100546818e-5) * r
                + 7.868691311456132591e-4)
                * r
                + 0.0148753612908506148525)
                * r
                + 0.13692988092273580531)
                * r
                + 0.59983220655588793769)
                * r
                + 1.0)
    };
    if q < 0.0 {
        -val
    } else {
        val
    }
}

/// ln Gamma(x), x > 0 (Lanczos, g = 7; relative error ~1e-15).
pub fn ln_gamma(x: f64) -> f64 {
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if !(x > 0.0) {
        return f64::NAN;
    }
    if x < 0.5 {
        // reflection
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = C[0];
    let t = x + 7.5;
    for (i, c) in C.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.918_938_533_204_672_8 + (x + 0.5) * t.ln() - t + a.ln()
}

/// Continued fraction of the incomplete beta function (Lentz).
fn beta_cf(a: f64, b: f64, x: f64) -> f64 {
    let tiny = 1e-300;
    let (qab, qap, qam) = (a + b, a + 1.0, a - 1.0);
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < tiny {
        d = tiny;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..20_000 {
        let m = m as f64;
        let m2 = 2.0 * m;
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < tiny {
            d = tiny;
        }
        c = 1.0 + aa / c;
        if c.abs() < tiny {
            c = tiny;
        }
        d = 1.0 / d;
        h *= d * c;
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < tiny {
            d = tiny;
        }
        c = 1.0 + aa / c;
        if c.abs() < tiny {
            c = tiny;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-16 {
            break;
        }
    }
    h
}

/// The regularized incomplete beta I_x(a, b) and its complement 1 - I_x(a, b),
/// each computed directly where it is small. `y` must be 1 - x (passed apart
/// so that it keeps its digits when x is near 1).
pub fn beta_reg(a: f64, b: f64, x: f64, y: f64) -> (f64, f64) {
    if x <= 0.0 {
        return (0.0, 1.0);
    }
    if y <= 0.0 {
        return (1.0, 0.0);
    }
    let ln_bt = ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * y.ln();
    if x < (a + 1.0) / (a + b + 2.0) {
        let lo = (ln_bt.exp() * beta_cf(a, b, x) / a).clamp(0.0, 1.0);
        (lo, 1.0 - lo)
    } else {
        let up = (ln_bt.exp() * beta_cf(b, a, y) / b).clamp(0.0, 1.0);
        (1.0 - up, up)
    }
}

/// Two-sided p of a t statistic with `df` degrees of freedom.
pub fn t_two_sided(t: f64, df: f64) -> f64 {
    if t.is_nan() || !(df > 0.0) {
        return f64::NAN;
    }
    if t.is_infinite() {
        return 0.0;
    }
    let t2 = t * t;
    beta_reg(df / 2.0, 0.5, df / (df + t2), t2 / (df + t2)).0
}

/// The t for which the two-sided p is `alpha` (R's qt(1 - alpha/2, df)).
pub fn t_critical(df: f64, alpha: f64) -> f64 {
    if !(df > 0.0) || !(alpha > 0.0 && alpha < 1.0) {
        return f64::NAN;
    }
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    while t_two_sided(hi, df) > alpha && hi < 1e150 {
        lo = hi;
        hi *= 2.0;
    }
    for _ in 0..300 {
        let mid = 0.5 * (lo + hi);
        if t_two_sided(mid, df) > alpha {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo <= 1e-15 * hi {
            break;
        }
    }
    0.5 * (lo + hi)
}

/// P(|T| > crit) for a non-central t with `df` degrees of freedom and
/// non-centrality `delta`: E_S[Q(crit S - delta) + Q(crit S + delta)], with
/// S = sqrt(V / df), V chi-square(df), integrated by Simpson's rule over the
/// density of S (accurate to ~1e-9).
pub fn noncentral_t_power(df: f64, delta: f64, crit: f64) -> f64 {
    if !(df > 0.0) || delta.is_nan() || crit.is_nan() {
        return f64::NAN;
    }
    let d = delta.abs();
    if df > 1e7 {
        return (normal_upper(crit - d) + normal_upper(crit + d)).min(1.0);
    }
    // density of S: 2 df s f_V(df s^2)
    let k = df / 2.0;
    let lconst = std::f64::consts::LN_2 + df.ln() + (k - 1.0) * df.ln() - k * std::f64::consts::LN_2 - ln_gamma(k);
    let dens = |s: f64| -> f64 {
        if s <= 0.0 {
            return if df == 1.0 { (lconst).exp() } else { 0.0 };
        }
        (lconst + (df - 1.0) * s.ln() - df * s * s / 2.0).exp()
    };
    let sd = 1.0 / (2.0 * df).sqrt();
    let lo = (1.0 - 14.0 * sd).max(0.0);
    let hi = 1.0 + 14.0 * sd + if df < 30.0 { 6.0 } else { 0.0 };
    let n = 6000usize;
    let h = (hi - lo) / n as f64;
    let mut soma = 0.0;
    for i in 0..=n {
        let s = lo + i as f64 * h;
        let f = dens(s) * (normal_upper(crit * s - d) + normal_upper(crit * s + d));
        let w = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        soma += w * f;
    }
    (soma * h / 3.0).clamp(0.0, 1.0)
}

/// Benjamini-Hochberg adjusted p-values (R's p.adjust(p, "BH")). NaN entries
/// stay NaN and do not count as tests.
pub fn benjamini_hochberg(p: &[f64]) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..p.len()).filter(|&i| !p[i].is_nan()).collect();
    let n = idx.len();
    let mut r = vec![f64::NAN; p.len()];
    if n == 0 {
        return r;
    }
    // decreasing p; ties keep the input order (stable sort)
    idx.sort_by(|a, b| p[*b].total_cmp(&p[*a]));
    let mut cmin = f64::INFINITY;
    for (k, &i) in idx.iter().enumerate() {
        let rank = (n - k) as f64;
        let v = (n as f64 / rank * p[i]).min(cmin);
        cmin = v;
        r[i] = v.min(1.0);
    }
    r
}

// ---------------------------------------------------------------------
//  RANDOM NUMBERS (xoshiro256**, seeded by splitmix64): one stream per
//  replicate, so that results do not depend on the number of threads.
// ---------------------------------------------------------------------

fn splitmix(x: &mut u64) -> u64 {
    *x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[derive(Clone, Debug)]
pub struct Rng {
    s: [u64; 4],
}

impl Rng {
    pub fn new(seed: u64) -> Rng {
        let mut x = seed;
        Rng { s: [splitmix(&mut x), splitmix(&mut x), splitmix(&mut x), splitmix(&mut x)] }
    }
    /// The k-th independent stream of `seed`.
    pub fn stream(seed: u64, k: u64) -> Rng {
        let mut x = k ^ 0x5851_F42D_4C95_7F2D;
        let m = splitmix(&mut x);
        Rng::new(seed.wrapping_mul(0x2545_F491_4F6C_DD1D) ^ m)
    }
    pub fn next_u64(&mut self) -> u64 {
        let r = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        r
    }
    /// Uniform in [0, 1).
    pub fn unif(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    }
    /// Uniform integer in [0, n).
    pub fn below(&mut self, n: usize) -> usize {
        ((self.next_u64() as u128 * n as u128) >> 64) as usize
    }
    pub fn normal(&mut self) -> f64 {
        let u1 = 1.0 - self.unif();
        let u2 = self.unif();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
    pub fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            let j = self.below(i + 1);
            v.swap(i, j);
        }
    }
}

// =====================================================================
//  2) LINEAR ALGEBRA (small dense, column-major designs)
// =====================================================================

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn axpy(a: f64, x: &[f64], y: &mut [f64]) {
    for (yi, xi) in y.iter_mut().zip(x) {
        *yi += a * xi;
    }
}

/// Thin QR by Gram-Schmidt with re-orthogonalization ("twice is enough").
/// A column whose part orthogonal to the earlier ones is below `tol` times
/// its norm is numerically a combination of them: it is dropped.
#[derive(Clone, Debug)]
struct Qr {
    q: Vec<Vec<f64>>,
    /// k x k upper triangular, row-major, for the kept columns.
    r: Vec<f64>,
    kept: Vec<usize>,
    dropped: Vec<usize>,
}

fn qr_gs(cols: &[Vec<f64>], tol: f64) -> Qr {
    let mut q: Vec<Vec<f64>> = Vec::with_capacity(cols.len());
    let mut rcols: Vec<Vec<f64>> = Vec::with_capacity(cols.len());
    let (mut kept, mut dropped) = (Vec::new(), Vec::new());
    for (j, c) in cols.iter().enumerate() {
        let n0 = dot(c, c).sqrt();
        if !(n0 > 0.0) || !n0.is_finite() {
            dropped.push(j);
            continue;
        }
        let mut v = c.clone();
        let mut coef = vec![0.0; q.len()];
        for _ in 0..2 {
            for (i, qi) in q.iter().enumerate() {
                let d = dot(qi, &v);
                axpy(-d, qi, &mut v);
                coef[i] += d;
            }
        }
        let nv = dot(&v, &v).sqrt();
        if !(nv > tol * n0) {
            dropped.push(j);
            continue;
        }
        for x in v.iter_mut() {
            *x /= nv;
        }
        coef.push(nv);
        q.push(v);
        rcols.push(coef);
        kept.push(j);
    }
    let k = q.len();
    let mut r = vec![0.0; k * k];
    for (jj, col) in rcols.iter().enumerate() {
        for (i, &v) in col.iter().enumerate() {
            r[i * k + jj] = v;
        }
    }
    Qr { q, r, kept, dropped }
}

impl Qr {
    fn k(&self) -> usize {
        self.q.len()
    }
    /// v minus its projection on the column space.
    fn residual(&self, v: &[f64]) -> Vec<f64> {
        let mut r = v.to_vec();
        for _ in 0..2 {
            for qi in &self.q {
                let d = dot(qi, &r);
                axpy(-d, qi, &mut r);
            }
        }
        r
    }
    /// (R'R)^-1 g.
    fn solve_rtr(&self, g: &[f64]) -> Vec<f64> {
        let k = self.k();
        // R' u = g
        let mut u = vec![0.0; k];
        for i in 0..k {
            let mut s = g[i];
            for j in 0..i {
                s -= self.r[j * k + i] * u[j];
            }
            u[i] = s / self.r[i * k + i];
        }
        // R x = u
        let mut x = vec![0.0; k];
        for i in (0..k).rev() {
            let mut s = u[i];
            for j in i + 1..k {
                s -= self.r[i * k + j] * x[j];
            }
            x[i] = s / self.r[i * k + i];
        }
        x
    }
    /// (R'R)^-1, k x k row-major.
    fn inv_rtr(&self) -> Vec<f64> {
        let k = self.k();
        // R^-1 (upper)
        let mut ri = vec![0.0; k * k];
        for c in 0..k {
            for i in (0..=c).rev() {
                let mut s = if i == c { 1.0 } else { 0.0 };
                for j in i + 1..=c {
                    s -= self.r[i * k + j] * ri[j * k + c];
                }
                ri[i * k + c] = s / self.r[i * k + i];
            }
        }
        let mut out = vec![0.0; k * k];
        for i in 0..k {
            for j in 0..k {
                let mut s = 0.0;
                for m in i.max(j)..k {
                    s += ri[i * k + m] * ri[j * k + m];
                }
                out[i * k + j] = s;
            }
        }
        out
    }
    /// R'R (the cross-product matrix), k x k row-major.
    fn rtr(&self) -> Vec<f64> {
        let k = self.k();
        let mut out = vec![0.0; k * k];
        for i in 0..k {
            for j in 0..k {
                let mut s = 0.0;
                for m in 0..=i.min(j) {
                    s += self.r[m * k + i] * self.r[m * k + j];
                }
                out[i * k + j] = s;
            }
        }
        out
    }
    fn ln_det_rtr_half(&self) -> f64 {
        let k = self.k();
        (0..k).map(|i| self.r[i * k + i].abs().ln()).sum()
    }
}

/// Solves A x = b for a symmetric positive-definite A (m x m, row-major).
fn chol_solve(a: &[f64], m: usize, b: &[f64]) -> Option<Vec<f64>> {
    let mut l = vec![0.0; m * m];
    for i in 0..m {
        for j in 0..=i {
            let mut s = a[i * m + j];
            for k in 0..j {
                s -= l[i * m + k] * l[j * m + k];
            }
            if i == j {
                if !(s > 0.0) {
                    return None;
                }
                l[i * m + i] = s.sqrt();
            } else {
                l[i * m + j] = s / l[j * m + j];
            }
        }
    }
    let mut z = vec![0.0; m];
    for i in 0..m {
        let mut s = b[i];
        for k in 0..i {
            s -= l[i * m + k] * z[k];
        }
        z[i] = s / l[i * m + i];
    }
    let mut x = vec![0.0; m];
    for i in (0..m).rev() {
        let mut s = z[i];
        for k in i + 1..m {
            s -= l[k * m + i] * x[k];
        }
        x[i] = s / l[i * m + i];
    }
    Some(x)
}

fn eta_de(cols: &[Vec<f64>], beta: &[f64], n: usize) -> Vec<f64> {
    let mut eta = vec![0.0; n];
    for (c, b) in cols.iter().zip(beta) {
        if *b != 0.0 {
            axpy(*b, c, &mut eta);
        }
    }
    eta
}

// =====================================================================
//  3) MODELS
// =====================================================================

/// log(1 + e^x) without overflow.
fn softplus(x: f64) -> f64 {
    if x > 0.0 {
        x + (-x).exp().ln_1p()
    } else {
        x.exp().ln_1p()
    }
}

fn sigmoid(x: f64) -> f64 {
    if x >= 0.0 {
        1.0 / (1.0 + (-x).exp())
    } else {
        let e = x.exp();
        e / (1.0 + e)
    }
}

/// p (1 - p) at the linear predictor x, accurate in the tails.
fn peso_logit(x: f64) -> f64 {
    let e = (-x.abs()).exp();
    e / ((1.0 + e) * (1.0 + e))
}

fn loglik_logit(eta: &[f64], y: &[f64]) -> f64 {
    eta.iter().zip(y).map(|(e, yi)| yi * e - softplus(*e)).sum()
}

/// A logistic regression fitted by maximum likelihood.
#[derive(Clone, Debug)]
pub struct LogitFit {
    pub beta: Vec<f64>,
    pub se: Vec<f64>,
    pub loglik: f64,
    pub converged: bool,
    /// The fit did not converge or some fitted probability is within e^-25 of
    /// 0 or 1: (quasi-)complete separation, the MLE is infinite.
    pub separated: bool,
    pub iterations: usize,
    pub fitted: Vec<f64>,
}

/// Logistic regression by Newton-Raphson (IRLS) with step-halving. `cols`
/// are the design columns (the first is the intercept); every column must be
/// linearly independent of the others.
pub fn logistic_fit(cols: &[Vec<f64>], y: &[f64], start: Option<&[f64]>) -> LogitFit {
    let n = y.len();
    let p = cols.len();
    let mut beta = match start {
        Some(s) if s.len() == p => s.to_vec(),
        _ => {
            let mut b = vec![0.0; p];
            let m = (y.iter().sum::<f64>() + 0.5) / (n as f64 + 1.0);
            b[0] = (m / (1.0 - m)).ln();
            b
        }
    };
    let mut eta = eta_de(cols, &beta, n);
    let mut ll = loglik_logit(&eta, y);
    let mut converged = false;
    let mut iters = 0;
    let mut singular = false;
    for it in 1..=100 {
        iters = it;
        let sw: Vec<f64> = eta.iter().map(|e| peso_logit(*e).sqrt()).collect();
        let colsw: Vec<Vec<f64>> = cols.iter().map(|c| c.iter().zip(&sw).map(|(x, w)| x * w).collect()).collect();
        let qr = qr_gs(&colsw, 1e-12);
        if !qr.dropped.is_empty() {
            singular = true;
            break;
        }
        let resid: Vec<f64> = eta.iter().zip(y).map(|(e, yi)| yi - sigmoid(*e)).collect();
        let g: Vec<f64> = cols.iter().map(|c| dot(c, &resid)).collect();
        let step = qr.solve_rtr(&g);
        let mut t = 1.0;
        let mut aceito = false;
        let (mut nb, mut ne, mut nll) = (beta.clone(), eta.clone(), ll);
        for _ in 0..40 {
            nb = beta.iter().zip(&step).map(|(b, s)| b + t * s).collect();
            ne = eta_de(cols, &nb, n);
            nll = loglik_logit(&ne, y);
            if nll.is_finite() && nll >= ll - 1e-10 * (ll.abs() + 1.0) {
                aceito = true;
                break;
            }
            t *= 0.5;
        }
        if !aceito {
            break;
        }
        let dmax = step.iter().fold(0.0f64, |a, s| a.max((t * s).abs()));
        let dll = nll - ll;
        beta = nb;
        eta = ne;
        ll = nll;
        if dmax < 1e-10 || (dll.abs() < 1e-14 * (ll.abs() + 1.0) && dmax < 1e-7) {
            converged = true;
            break;
        }
    }
    let sw: Vec<f64> = eta.iter().map(|e| peso_logit(*e).sqrt()).collect();
    let colsw: Vec<Vec<f64>> = cols.iter().map(|c| c.iter().zip(&sw).map(|(x, w)| x * w).collect()).collect();
    let qr = qr_gs(&colsw, 1e-12);
    let se = if qr.dropped.is_empty() {
        let c = qr.inv_rtr();
        (0..p).map(|i| c[i * p + i].sqrt()).collect()
    } else {
        vec![f64::NAN; p]
    };
    let eta_max = eta.iter().fold(0.0f64, |a, e| a.max(e.abs()));
    let fitted = eta.iter().map(|e| sigmoid(*e)).collect();
    LogitFit {
        separated: !converged || singular || eta_max > 25.0,
        converged: converged && !singular,
        beta,
        se,
        loglik: ll,
        iterations: iters,
        fitted,
    }
}

/// A Firth (bias-reduced, Jeffreys-penalized) logistic regression.
#[derive(Clone, Debug)]
pub struct FirthFit {
    pub beta: Vec<f64>,
    pub se: Vec<f64>,
    /// l(b) + 1/2 log |I(b)|.
    pub loglik_pen: f64,
    pub converged: bool,
    pub iterations: usize,
}

struct EstadoFirth {
    eta: Vec<f64>,
    qr: Qr,
    llpen: f64,
}

fn avalia_firth(cols: &[Vec<f64>], y: &[f64], beta: &[f64]) -> Option<EstadoFirth> {
    let n = y.len();
    let eta = eta_de(cols, beta, n);
    if eta.iter().any(|e| !e.is_finite() || e.abs() > 700.0) {
        return None;
    }
    let sw: Vec<f64> = eta.iter().map(|e| peso_logit(*e).sqrt()).collect();
    let colsw: Vec<Vec<f64>> = cols.iter().map(|c| c.iter().zip(&sw).map(|(x, w)| x * w).collect()).collect();
    let qr = qr_gs(&colsw, 1e-13);
    if !qr.dropped.is_empty() {
        return None;
    }
    let llpen = loglik_logit(&eta, y) + qr.ln_det_rtr_half();
    if !llpen.is_finite() {
        return None;
    }
    Some(EstadoFirth { eta, qr, llpen })
}

/// Firth's penalized logistic regression (as R's logistf): Newton steps on
/// the modified score with the Fisher information, at most 5 per coefficient,
/// halved until the penalized likelihood does not decrease. `fixed` holds one
/// coefficient at a value (the restricted fits of the penalized likelihood
/// ratio test and of the profile CI); the penalty is always that of the full
/// design, as in logistf.
pub fn firth_fit(cols: &[Vec<f64>], y: &[f64], fixed: Option<(usize, f64)>, start: Option<&[f64]>) -> Option<FirthFit> {
    let n = y.len();
    let p = cols.len();
    let mut beta = match start {
        Some(s) if s.len() == p => s.to_vec(),
        _ => {
            let mut b = vec![0.0; p];
            let m = (y.iter().sum::<f64>() + 0.5) / (n as f64 + 1.0);
            b[0] = (m / (1.0 - m)).ln();
            b
        }
    };
    if let Some((j, v)) = fixed {
        beta[j] = v;
    }
    let mut st = match avalia_firth(cols, y, &beta) {
        Some(s) => s,
        None => {
            // a start too far out: begin at zero
            beta = vec![0.0; p];
            if let Some((j, v)) = fixed {
                beta[j] = v;
            }
            avalia_firth(cols, y, &beta)?
        }
    };
    let livres: Vec<usize> = (0..p).filter(|&k| fixed.map_or(true, |(j, _)| j != k)).collect();
    let mut converged = false;
    let mut iters = 0;
    for it in 1..=300 {
        iters = it;
        let h: Vec<f64> = (0..n).map(|i| st.qr.q.iter().map(|qk| qk[i] * qk[i]).sum()).collect();
        let resid: Vec<f64> = (0..n)
            .map(|i| {
                let pi = sigmoid(st.eta[i]);
                y[i] - pi + h[i] * (0.5 - pi)
            })
            .collect();
        let u: Vec<f64> = cols.iter().map(|c| dot(c, &resid)).collect();
        let mut step = vec![0.0; p];
        if livres.len() == p {
            step = st.qr.solve_rtr(&u);
        } else if !livres.is_empty() {
            let inf = st.qr.rtr();
            let m = livres.len();
            let mut a = vec![0.0; m * m];
            for (a_i, &i) in livres.iter().enumerate() {
                for (a_j, &j) in livres.iter().enumerate() {
                    a[a_i * m + a_j] = inf[i * p + j];
                }
            }
            let uf: Vec<f64> = livres.iter().map(|&i| u[i]).collect();
            let sf = chol_solve(&a, m, &uf)?;
            for (a_i, &i) in livres.iter().enumerate() {
                step[i] = sf[a_i];
            }
        }
        let umax = livres.iter().fold(0.0f64, |a, &i| a.max(u[i].abs()));
        let smax = step.iter().fold(0.0f64, |a, s| a.max(s.abs()));
        if smax > 5.0 {
            for s in step.iter_mut() {
                *s *= 5.0 / smax;
            }
        }
        if smax < 1e-11 && umax < 1e-8 {
            converged = true;
            break;
        }
        let mut t = 1.0;
        let mut novo = None;
        for _ in 0..50 {
            let nb: Vec<f64> = beta.iter().zip(&step).map(|(b, s)| b + t * s).collect();
            if let Some(ns) = avalia_firth(cols, y, &nb) {
                if ns.llpen >= st.llpen - 1e-11 * (st.llpen.abs() + 1.0) {
                    novo = Some((nb, ns));
                    break;
                }
            }
            t *= 0.5;
        }
        let Some((nb, ns)) = novo else {
            // no step improves: at the maximum within rounding
            converged = umax < 1e-6;
            break;
        };
        let dmax = step.iter().fold(0.0f64, |a, s| a.max((t * s).abs()));
        beta = nb;
        st = ns;
        if dmax < 1e-11 {
            converged = true;
            break;
        }
    }
    let c = st.qr.inv_rtr();
    let se = (0..p).map(|i| c[i * p + i].sqrt()).collect();
    Some(FirthFit { beta, se, loglik_pen: st.llpen, converged, iterations: iters })
}

/// Penalized likelihood ratio test of coefficient `j` (fixed at 0, the
/// penalty of the full design): (statistic, p).
pub fn firth_lrt(cols: &[Vec<f64>], y: &[f64], fit: &FirthFit, j: usize) -> (f64, f64) {
    let mut inicio = fit.beta.clone();
    inicio[j] = 0.0;
    match firth_fit(cols, y, Some((j, 0.0)), Some(&inicio)) {
        Some(r) => {
            let s = (2.0 * (fit.loglik_pen - r.loglik_pen)).max(0.0);
            (s, chi2_1_upper(s))
        }
        None => (f64::NAN, f64::NAN),
    }
}

/// Profile penalized-likelihood confidence interval of coefficient `j` (as
/// logistf): the values b where 2 (l*(b_hat) - max l* with b_j = b) equals the
/// chi-square quantile `q` (3.84 for 95%). An infinite bound is reported as
/// +-Inf.
pub fn firth_profile_ci(cols: &[Vec<f64>], y: &[f64], fit: &FirthFit, j: usize, q: f64) -> (f64, f64) {
    let bj = fit.beta[j];
    let mut quente = fit.beta.clone();
    let f = |b: f64, quente: &mut Vec<f64>| -> f64 {
        let mut ini = quente.clone();
        ini[j] = b;
        match firth_fit(cols, y, Some((j, b)), Some(&ini)) {
            Some(r) => {
                *quente = r.beta.clone();
                2.0 * (fit.loglik_pen - r.loglik_pen) - q
            }
            None => f64::INFINITY,
        }
    };
    let mut lim = [f64::NAN; 2];
    for (k, lado) in [-1.0f64, 1.0].iter().enumerate() {
        quente.clone_from(&fit.beta);
        let mut passo = if fit.se[j].is_finite() && fit.se[j] > 0.0 { fit.se[j] } else { 1.0 };
        let (mut a, mut fa) = (bj, -q);
        let mut b = bj + lado * passo;
        let mut fb = f(b, &mut quente);
        let mut dobras = 0;
        while fb < 0.0 && dobras < 40 {
            a = b;
            fa = fb;
            passo *= 2.0;
            b = bj + lado * passo;
            fb = f(b, &mut quente);
            dobras += 1;
            if passo > 1e3 {
                break;
            }
        }
        if fb < 0.0 {
            lim[k] = lado * f64::INFINITY;
            continue;
        }
        // Illinois (regula falsi) between a (f < 0) and b (f >= 0)
        let mut lado_ant = 0i8;
        for _ in 0..200 {
            let c = if fb.is_finite() { (a * fb - b * fa) / (fb - fa) } else { 0.5 * (a + b) };
            let c = if c.is_finite() && (c - a) * (c - b) < 0.0 { c } else { 0.5 * (a + b) };
            let fc = f(c, &mut quente);
            if fc < 0.0 {
                a = c;
                fa = fc;
                if lado_ant == -1 {
                    fb *= 0.5;
                }
                lado_ant = -1;
            } else {
                b = c;
                fb = fc;
                if lado_ant == 1 {
                    fa *= 0.5;
                }
                lado_ant = 1;
            }
            if (b - a).abs() <= 1e-10 * (1.0 + a.abs()) || fc.abs() < 1e-11 {
                break;
            }
        }
        lim[k] = if fa.abs() < fb.abs() { a } else { b };
    }
    (lim[0], lim[1])
}

// =====================================================================
//  4) THE PHENOTYPE TABLE
// =====================================================================

fn faltante(v: &str) -> bool {
    matches!(v.trim().to_ascii_lowercase().as_str(), "" | "na" | "n/a" | "nan" | "." | "null" | "none" | "missing")
}

fn divide(t: &str, sep: Option<char>) -> Vec<String> {
    match sep {
        Some(s) => t.split(s).map(|c| c.trim().trim_matches('"').trim().to_string()).collect(),
        None => t.split_whitespace().map(|c| c.trim_matches('"').to_string()).collect(),
    }
}

struct Tabela {
    caminho: String,
    cab: Vec<String>,
    linhas: Vec<(usize, Vec<String>)>,
    col_amostra: usize,
    nota: Option<String>,
}

impl Tabela {
    fn coluna(&self, nome: &str) -> Result<usize, String> {
        self.cab
            .iter()
            .position(|c| c == nome)
            .or_else(|| self.cab.iter().position(|c| c.eq_ignore_ascii_case(nome)))
            .ok_or_else(|| format!("The phenotype table {} has no column \"{nome}\"; its columns are: {}", self.caminho, self.cab.join(", ")))
    }
}

fn le_tabela(caminho: &str, sample_column: Option<&str>) -> Result<Tabela, String> {
    let mut r = abre_leitura(caminho)?;
    let mut buf = String::new();
    let mut n = 0usize;
    let mut cab: Option<(Vec<String>, Option<char>)> = None;
    let mut linhas = Vec::new();
    loop {
        buf.clear();
        let lidos = r.read_line(&mut buf).map_err(|e| format!("Error reading {caminho} after line {n}: {e}"))?;
        if lidos == 0 {
            break;
        }
        n += 1;
        let t = buf.trim_end_matches(['\n', '\r']).trim_start_matches('\u{feff}');
        if t.trim().is_empty() || t.starts_with("##") || t.starts_with("# ") {
            continue;
        }
        match &cab {
            None => {
                let t = t.strip_prefix('#').unwrap_or(t);
                let sep = if t.contains('\t') {
                    Some('\t')
                } else if t.contains(',') {
                    Some(',')
                } else if t.contains(';') {
                    Some(';')
                } else {
                    None
                };
                cab = Some((divide(t, sep), sep));
            }
            Some((c, sep)) => {
                let v = divide(t, *sep);
                if v.len() != c.len() {
                    return Err(format!(
                        "{caminho}, line {n}: {} column(s), but the header has {}. Every line of the phenotype table needs one value per column (write NA for a missing value).",
                        v.len(),
                        c.len()
                    ));
                }
                linhas.push((n, v));
            }
        }
    }
    let Some((cab, _)) = cab else {
        return Err(format!(
            "{caminho} is empty: a phenotype table has a header line (sample, outcome, covariates...) and one line per person, separated by tabs or commas"
        ));
    };
    if cab.len() < 2 {
        return Err(format!(
            "{caminho}: the header has a single column ({}); a phenotype table has the sample name and at least the outcome, separated by tabs or commas",
            cab.join("")
        ));
    }
    let mut t = Tabela { caminho: caminho.to_string(), cab, linhas, col_amostra: 0, nota: None };
    match sample_column {
        Some(s) if !s.trim().is_empty() => t.col_amostra = t.coluna(s.trim())?,
        _ => {
            if t.cab[0].eq_ignore_ascii_case("FID") {
                if let Some(i) = t.cab.iter().position(|c| c.eq_ignore_ascii_case("IID")) {
                    t.col_amostra = i;
                    t.nota = Some("the table has FID and IID columns (plink style): people are matched by IID".into());
                }
            }
        }
    }
    Ok(t)
}

// =====================================================================
//  5) READING VCF RECORDS
// =====================================================================

/// A genotype: allele indices; `ploidia` 0 is missing.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Gt {
    a: [u16; 2],
    ploidia: u8,
}

/// `0/1`, `1|2`, `1` (haploid). A `.` allele, more than two alleles or junk:
/// missing.
fn le_gt(s: &str) -> Gt {
    let mut alelos = [0u16; 2];
    let mut n = 0usize;
    let mut atual: Option<u32> = None;
    for c in s.bytes() {
        match c {
            b'0'..=b'9' => {
                let v = atual.unwrap_or(0) * 10 + u32::from(c - b'0');
                if v > u32::from(u16::MAX) {
                    return Gt::default();
                }
                atual = Some(v);
            }
            b'/' | b'|' => {
                let Some(v) = atual else { return Gt::default() };
                if n == 2 {
                    return Gt::default();
                }
                alelos[n] = v as u16;
                n += 1;
                atual = None;
            }
            _ => return Gt::default(),
        }
    }
    let Some(v) = atual else { return Gt::default() };
    if n == 2 {
        return Gt::default();
    }
    alelos[n] = v as u16;
    Gt { a: alelos, ploidia: (n + 1) as u8 }
}

/// Genotype code for one allele: -1 missing; 0, 1, 2 copies in a diploid
/// call; 3 / 4: 0 / 1 copy in a haploid call.
fn codigo(g: Gt, alelo: u16) -> i8 {
    match g.ploidia {
        2 => (g.a[0] == alelo) as i8 + (g.a[1] == alelo) as i8,
        1 => 3 + (g.a[0] == alelo) as i8,
        _ => -1,
    }
}

fn copias(c: i8) -> f64 {
    match c {
        0..=2 => c as f64,
        3 => 0.0,
        4 => 1.0,
        _ => f64::NAN,
    }
}

fn ploidia_de(c: i8) -> f64 {
    match c {
        0..=2 => 2.0,
        3 | 4 => 1.0,
        _ => 0.0,
    }
}

fn codifica(c: i8, m: GeneticModel) -> f64 {
    match m {
        GeneticModel::Additive => copias(c),
        GeneticModel::Dominant => (copias(c) >= 1.0) as u8 as f64,
        GeneticModel::Recessive => (c == 2 || c == 4) as u8 as f64,
    }
}

/// The people of the analysis one VCF holds: (sample column, 0-based among
/// the samples; person index), by column.
struct Arquivo {
    servidos: Vec<(usize, u32)>,
}

/// One record: the line, with the fixed fields located; the genotypes are
/// read only when the record is used.
struct Registro {
    arq: usize,
    linha: String,
    canon: String,
    pos: u64,
    id: (usize, usize),
    refe: String,
    alts: Vec<String>,
    passa: bool,
    formato: Option<(usize, usize)>,
    amostras_ini: Option<usize>,
}

impl Registro {
    fn bruto(&self) -> &str {
        self.linha.split('\t').next().unwrap_or("")
    }
    fn id(&self) -> &str {
        &self.linha[self.id.0..self.id.1]
    }
    fn ids(&self) -> impl Iterator<Item = &str> {
        self.id().split(';').filter(|s| !s.is_empty() && *s != ".")
    }
    /// One genotype per served person of the file, in `servidos` order.
    fn gts(&self, arq: &Arquivo) -> Vec<Gt> {
        let mut out = vec![Gt::default(); arq.servidos.len()];
        let (Some((f0, f1)), Some(mut p)) = (self.formato, self.amostras_ini) else {
            return out;
        };
        let Some(gi) = self.linha[f0..f1].split(':').position(|k| k == "GT") else {
            return out;
        };
        let b = self.linha.as_bytes();
        let mut col = 0usize;
        for (k, &(alvo, _)) in arq.servidos.iter().enumerate() {
            while col < alvo && p < b.len() {
                if b[p] == b'\t' {
                    col += 1;
                }
                p += 1;
            }
            if col < alvo {
                break;
            }
            let mut fim = p;
            while fim < b.len() && b[fim] != b'\t' {
                fim += 1;
            }
            out[k] = self.linha[p..fim].split(':').nth(gi).map_or(Gt::default(), le_gt);
        }
        out
    }
}

/// Splits the fixed columns. `None`: not a usable record (counted).
fn analisa_linha(linha: String, arq: usize, cache: &mut (String, String)) -> Option<Registro> {
    let mut lim = [(0usize, 0usize); 9];
    let mut ncol = 0usize;
    let mut ini = 0usize;
    while ncol < 9 {
        match linha[ini..].find('\t') {
            Some(k) => {
                lim[ncol] = (ini, ini + k);
                ini += k + 1;
                ncol += 1;
            }
            None => {
                lim[ncol] = (ini, linha.len());
                ncol += 1;
                break;
            }
        }
    }
    if ncol < 8 {
        return None;
    }
    let pos: u64 = linha[lim[1].0..lim[1].1].trim().parse().ok()?;
    let bruto = &linha[lim[0].0..lim[0].1];
    if bruto != cache.0 {
        cache.0 = bruto.to_string();
        cache.1 = canonico(bruto);
    }
    let refe = linha[lim[3].0..lim[3].1].to_ascii_uppercase();
    let alts: Vec<String> = linha[lim[4].0..lim[4].1].split(',').map(|a| a.to_ascii_uppercase()).collect();
    let filtro = &linha[lim[6].0..lim[6].1];
    let amostras_ini = if ncol == 9 && lim[8].1 < linha.len() { Some(lim[8].1 + 1) } else { None };
    Some(Registro {
        arq,
        canon: cache.1.clone(),
        pos,
        id: lim[2],
        refe,
        alts,
        passa: filtro == "PASS" || filtro == ".",
        formato: if ncol == 9 { Some(lim[8]) } else { None },
        amostras_ini,
        linha,
    })
}

/// A usable ALT allele: bases (ACGTN), not `.`, `*`, `<DEL>` or a breakend.
fn alelo_util(a: &str) -> bool {
    !a.is_empty() && a.bytes().all(|c| matches!(c, b'A' | b'C' | b'G' | b'T' | b'N'))
}

/// Reads a VCF header: the sample names, and the reader positioned at the
/// first record (returned apart).
fn abre_vcf(caminho: &str) -> Result<(Box<dyn BufRead>, Vec<String>, usize), String> {
    let mut r = abre_leitura(caminho)?;
    let mut buf = String::new();
    let mut n = 0usize;
    loop {
        buf.clear();
        let lidos = r.read_line(&mut buf).map_err(|e| format!("Error reading {caminho} after line {n}: {e}"))?;
        if lidos == 0 {
            return Err(format!(
                "{caminho}: no `#CHROM` header line was found. This is not a VCF (or it is empty); validate_file(...) tells what the file is."
            ));
        }
        n += 1;
        let t = buf.trim_end_matches(['\n', '\r']);
        if t.starts_with("##") || t.is_empty() {
            continue;
        }
        if t.starts_with("#CHROM") {
            let amostras: Vec<String> = t.split('\t').skip(9).map(str::to_string).collect();
            if amostras.is_empty() {
                return Err(format!(
                    "{caminho} has no sample columns (a sites-only VCF): an association needs each person's genotype (GT)."
                ));
            }
            return Ok((r, amostras, n));
        }
        if t.starts_with('#') {
            continue;
        }
        return Err(format!(
            "{caminho}, line {n}: a record comes before the `#CHROM` header line. This is not a VCF, or its header was cut off."
        ));
    }
}

// ---------------------------------------------------------------------
//  Record sources: a stream, or the lines of an indexed window
// ---------------------------------------------------------------------

enum Fonte {
    Fluxo(Box<dyn BufRead>),
    Linhas(std::vec::IntoIter<String>),
}

struct Cursor {
    arq: usize,
    caminho: String,
    fonte: Fonte,
    n_linha: usize,
    cab: Option<Registro>,
    cache: (String, String),
    /// Several files: natural chromosome order is checked and sequences that
    /// are not primary chromosomes (and the mitochondria) are skipped.
    multiplos: bool,
    ultimo: Option<(String, u64)>,
    nome_ultimo: String,
    registros: u64,
    malformados: u64,
    pulados_seq: u64,
}

impl Cursor {
    fn le(&mut self, buf: &mut String) -> Result<bool, String> {
        match &mut self.fonte {
            Fonte::Fluxo(r) => {
                buf.clear();
                let n = r.read_line(buf).map_err(|e| {
                    format!(
                        "Error reading {} after line {}: {e}. The file may be truncated or not a VCF.",
                        self.caminho, self.n_linha
                    )
                })?;
                if n == 0 {
                    return Ok(false);
                }
                self.n_linha += 1;
                while buf.ends_with('\n') || buf.ends_with('\r') {
                    buf.pop();
                }
                Ok(true)
            }
            Fonte::Linhas(it) => match it.next() {
                Some(l) => {
                    *buf = l;
                    while buf.ends_with('\n') || buf.ends_with('\r') {
                        buf.pop();
                    }
                    Ok(true)
                }
                None => Ok(false),
            },
        }
    }

    fn avanca(&mut self) -> Result<(), String> {
        loop {
            let mut l = String::new();
            if !self.le(&mut l)? {
                self.cab = None;
                return Ok(());
            }
            if l.is_empty() || l.starts_with('#') {
                continue;
            }
            self.registros += 1;
            let Some(r) = analisa_linha(l, self.arq, &mut self.cache) else {
                self.malformados += 1;
                continue;
            };
            if self.multiplos {
                if r.canon == "MT" || r.canon == "M" || !parece_cromossomo(&r.canon) {
                    self.pulados_seq += 1;
                    continue;
                }
                if let Some((c, p)) = &self.ultimo {
                    let ordem = if *c == r.canon { r.pos.cmp(p) } else { natural_cmp(&r.canon, c) };
                    if ordem == std::cmp::Ordering::Less {
                        let onde = if *c == r.canon {
                            format!("position {} on `{}` comes after position {p}", r.pos, r.bruto())
                        } else {
                            format!("chromosome `{}` comes after `{}`", r.bruto(), self.nome_ultimo)
                        };
                        return Err(format!(
                            "{}, line {}: {onde}. Several VCFs are joined by position, so each must be sorted by position with the chromosomes in natural order (1, 2, ..., 10, ..., 22, X, Y; `chr1` and `1` are the same). Sort it first: sort_genomic([\"{}\"], \"sorted.vcf.gz\", chrom_order=\"natural\"); or merge the files into one multi-sample VCF.",
                            self.caminho, self.n_linha, self.caminho
                        ));
                    }
                }
                if self.ultimo.as_ref().map_or(true, |(c, _)| *c != r.canon) {
                    self.nome_ultimo = r.bruto().to_string();
                }
                self.ultimo = Some((r.canon.clone(), r.pos));
            }
            self.cab = Some(r);
            return Ok(());
        }
    }
}

// =====================================================================
//  6) WHICH VARIANTS
// =====================================================================

enum Pedido {
    Rs(String),
    Pos(String, u64),
    Exato(String, u64, String, String),
}

fn le_pedido(s: &str) -> Result<Pedido, String> {
    let t = s.trim();
    if t.is_empty() {
        return Err("An empty variant ID in `ids`".into());
    }
    let pos = |v: &str| -> Option<u64> { v.replace(',', "").trim().parse::<u64>().ok().filter(|p| *p > 0) };
    let alelos_ok = |r: &str, a: &str| alelo_util(&r.to_ascii_uppercase()) && alelo_util(&a.to_ascii_uppercase());
    if t.contains(':') {
        let p: Vec<&str> = t.split(':').collect();
        match p.len() {
            2 => {
                if let Some(x) = pos(p[1]) {
                    return Ok(Pedido::Pos(canonico(p[0].trim()), x));
                }
                if p[1].contains('-') {
                    return Err(format!("\"{t}\" looks like a region: pass it in `regions`, not in `ids`"));
                }
            }
            3 => {
                if let (Some(x), Some((r, a))) = (pos(p[1]), p[2].split_once('>')) {
                    if alelos_ok(r, a) {
                        return Ok(Pedido::Exato(canonico(p[0].trim()), x, r.trim().to_ascii_uppercase(), a.trim().to_ascii_uppercase()));
                    }
                }
            }
            4 => {
                if let Some(x) = pos(p[1]) {
                    if alelos_ok(p[2].trim(), p[3].trim()) {
                        return Ok(Pedido::Exato(canonico(p[0].trim()), x, p[2].trim().to_ascii_uppercase(), p[3].trim().to_ascii_uppercase()));
                    }
                }
            }
            _ => {}
        }
        return Err(format!(
            "Variant \"{t}\" is not understood: use an rsID (rs429358), \"chr:pos\" (all alleles at that position), \"chr:pos:ref:alt\" or \"chr-pos-ref-alt\""
        ));
    }
    let p: Vec<&str> = t.split('-').collect();
    if p.len() == 4 {
        if let Some(x) = pos(p[1]) {
            if alelos_ok(p[2], p[3]) {
                return Ok(Pedido::Exato(canonico(p[0]), x, p[2].to_ascii_uppercase(), p[3].to_ascii_uppercase()));
            }
        }
    }
    Ok(Pedido::Rs(t.to_ascii_lowercase()))
}

struct Selecao {
    todas: bool,
    /// canonical chromosome -> merged 1-based closed intervals
    regioes: HashMap<String, Vec<(u64, u64)>>,
    rs: HashMap<String, Vec<usize>>,
    pos: HashMap<(String, u64), Vec<usize>>,
    exatos: HashMap<(String, u64, String, String), Vec<usize>>,
    pos_exatos: HashSet<(String, u64)>,
    pedidos: Vec<String>,
    achados: Vec<bool>,
}

impl Selecao {
    fn nova(ids: &[String], regioes: &[String]) -> Result<Selecao, String> {
        let mut s = Selecao {
            todas: ids.is_empty() && regioes.is_empty(),
            regioes: HashMap::new(),
            rs: HashMap::new(),
            pos: HashMap::new(),
            exatos: HashMap::new(),
            pos_exatos: HashSet::new(),
            pedidos: Vec::new(),
            achados: Vec::new(),
        };
        for (i, t) in ids.iter().enumerate() {
            s.pedidos.push(t.trim().to_string());
            s.achados.push(false);
            match le_pedido(t)? {
                Pedido::Rs(r) => s.rs.entry(r).or_default().push(i),
                Pedido::Pos(c, p) => s.pos.entry((c, p)).or_default().push(i),
                Pedido::Exato(c, p, r, a) => {
                    s.pos_exatos.insert((c.clone(), p));
                    s.exatos.entry((c, p, r, a)).or_default().push(i)
                }
            }
        }
        for r in regioes {
            let (nome, ini0, fim) = le_regiao(r, &|_| false)?;
            s.regioes.entry(canonico(&nome)).or_default().push((ini0 as u64 + 1, fim as u64));
        }
        for v in s.regioes.values_mut() {
            v.sort_unstable();
            let mut u: Vec<(u64, u64)> = Vec::with_capacity(v.len());
            for &(a, b) in v.iter() {
                match u.last_mut() {
                    Some(l) if a <= l.1.saturating_add(1) => l.1 = l.1.max(b),
                    _ => u.push((a, b)),
                }
            }
            *v = u;
        }
        Ok(s)
    }

    fn em_regiao(&self, canon: &str, pos: u64) -> bool {
        let Some(v) = self.regioes.get(canon) else { return false };
        let i = v.partition_point(|iv| iv.1 < pos);
        i < v.len() && v[i].0 <= pos
    }

    /// Positions to look at, for the indexed reading: the regions and the
    /// positional IDs. `None` when rsIDs are asked (they can be anywhere).
    fn alvos(&self) -> Option<Vec<(String, u64, u64)>> {
        if !self.rs.is_empty() || self.todas {
            return None;
        }
        let mut v: Vec<(String, u64, u64)> = Vec::new();
        for (c, iv) in &self.regioes {
            for &(a, b) in iv {
                v.push((c.clone(), a, b));
            }
        }
        for (c, p) in self.pos.keys() {
            v.push((c.clone(), *p, *p));
        }
        for (c, p, _, _) in self.exatos.keys() {
            v.push((c.clone(), *p, *p));
        }
        v.sort_by(|a, b| natural_cmp(&a.0, &b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
        let mut u: Vec<(String, u64, u64)> = Vec::new();
        for x in v {
            match u.last_mut() {
                Some(l) if l.0 == x.0 && x.1 <= l.2.saturating_add(1) => l.2 = l.2.max(x.2),
                _ => u.push(x),
            }
        }
        Some(u)
    }
}

// =====================================================================
//  7) THE ANALYSIS SET: OUTCOME AND DESIGN
// =====================================================================

struct Analise {
    nomes: Vec<String>,
    y: Vec<f64>,
    /// Design columns, the intercept first.
    z: Vec<Vec<f64>>,
    z_nomes: Vec<String>,
    binario: bool,
    evento: Option<String>,
    outro: Option<String>,
    covariaveis: Vec<String>,
    // counts
    linhas_tabela: usize,
    amostras_vcf: usize,
    fora_da_tabela: usize,
    fora_do_vcf: Vec<String>,
    falta_desfecho: usize,
    falta_cov: Vec<(String, usize)>,
    avisos: Vec<String>,
}

fn e_numero(v: &str) -> Option<f64> {
    v.trim().parse::<f64>().ok().filter(|x| x.is_finite())
}

fn monta_analise(tab: &Tabela, amostras_vcf: &HashMap<String, usize>, opt: &AssocOptions) -> Result<Analise, String> {
    let caminho = &tab.caminho;
    let ia = tab.col_amostra;
    let iy = tab.coluna(opt.outcome.trim())?;
    if iy == ia {
        return Err(format!("The outcome column \"{}\" is the sample column itself", opt.outcome));
    }
    let mut ic = Vec::new();
    for c in &opt.covariates {
        let i = tab.coluna(c.trim())?;
        if i == iy {
            return Err(format!("The covariate \"{c}\" is the outcome itself"));
        }
        if i == ia {
            return Err(format!("The covariate \"{c}\" is the sample column"));
        }
        if ic.iter().any(|(_, j)| *j == i) {
            return Err(format!("The covariate \"{c}\" is listed twice"));
        }
        ic.push((tab.cab[i].clone(), i));
    }
    for c in &opt.categorical {
        if !ic.iter().any(|(n, _)| n.eq_ignore_ascii_case(c.trim())) {
            return Err(format!("\"{c}\" is listed in `categorical` but not in `covariates`"));
        }
    }
    let mut vistos: HashMap<&str, usize> = HashMap::new();
    let mut an = Analise {
        nomes: Vec::new(),
        y: Vec::new(),
        z: Vec::new(),
        z_nomes: Vec::new(),
        binario: false,
        evento: None,
        outro: None,
        covariaveis: ic.iter().map(|(n, _)| n.clone()).collect(),
        linhas_tabela: tab.linhas.len(),
        amostras_vcf: amostras_vcf.len(),
        fora_da_tabela: 0,
        fora_do_vcf: Vec::new(),
        falta_desfecho: 0,
        falta_cov: ic.iter().map(|(n, _)| (n.clone(), 0)).collect(),
        avisos: Vec::new(),
    };
    if let Some(n) = &tab.nota {
        an.avisos.push(n.clone());
    }
    let mut na_tabela = HashSet::new();
    let mut sem_cov = 0usize;
    let mut linhas_ok: Vec<&Vec<String>> = Vec::new();
    for (ln, v) in &tab.linhas {
        let s = v[ia].as_str();
        if s.is_empty() || faltante(s) {
            return Err(format!("{caminho}, line {ln}: no sample name in column \"{}\"", tab.cab[ia]));
        }
        if let Some(l0) = vistos.insert(s, *ln) {
            return Err(format!("{caminho}: sample \"{s}\" appears twice (lines {l0} and {ln}); each person must have one line"));
        }
        na_tabela.insert(s.to_string());
        if !amostras_vcf.contains_key(s) {
            an.fora_do_vcf.push(s.to_string());
            continue;
        }
        if faltante(&v[iy]) {
            an.falta_desfecho += 1;
            continue;
        }
        let mut falta = false;
        for (k, (_, i)) in ic.iter().enumerate() {
            if faltante(&v[*i]) {
                an.falta_cov[k].1 += 1;
                falta = true;
            }
        }
        if falta {
            sem_cov += 1;
            continue;
        }
        linhas_ok.push(v);
    }
    an.fora_da_tabela = amostras_vcf.keys().filter(|s| !na_tabela.contains(*s)).count();
    if na_tabela.iter().all(|s| !amostras_vcf.contains_key(s)) {
        let mut ex_vcf: Vec<String> = amostras_vcf.keys().cloned().collect();
        ex_vcf.sort();
        let mut ex_tab: Vec<String> = tab.linhas.iter().map(|(_, v)| v[ia].clone()).collect();
        ex_tab.truncate(5);
        return Err(format!(
            "No sample of the VCF(s) is in the phenotype table {caminho}: people are matched by name. The VCF has {}; the table's column \"{}\" has {}.",
            lista_curta(&ex_vcf, 5),
            tab.cab[ia],
            lista_curta(&ex_tab, 5)
        ));
    }
    if linhas_ok.is_empty() {
        return Err(format!(
            "Nobody is left to analyse: of the people in both the VCF(s) and {caminho}, {} have no value of the outcome and {} lack some covariate.",
            an.falta_desfecho,
            sem_cov
        ));
    }
    an.nomes = linhas_ok.iter().map(|v| v[ia].clone()).collect();
    let n = an.nomes.len();

    // --- the outcome ---
    let valores: Vec<&str> = linhas_ok.iter().map(|v| v[iy].as_str()).collect();
    let mut distintos: Vec<&str> = valores.clone();
    distintos.sort_unstable();
    distintos.dedup();
    let todos_num = valores.iter().all(|v| e_numero(v).is_some());
    let binario = match opt.outcome_type {
        OutcomeType::Binary => {
            if distintos.len() != 2 {
                return Err(format!(
                    "The outcome \"{}\" was declared binary but has {} distinct value(s) among the people analysed: {}",
                    opt.outcome,
                    distintos.len(),
                    lista_curta(&distintos.iter().map(|s| s.to_string()).collect::<Vec<_>>(), 6)
                ));
            }
            true
        }
        OutcomeType::Quantitative => {
            if !todos_num {
                let mau = valores.iter().find(|v| e_numero(v).is_none()).unwrap();
                return Err(format!("The outcome \"{}\" was declared quantitative but has the value \"{mau}\", which is not a number", opt.outcome));
            }
            false
        }
        OutcomeType::Auto => {
            if distintos.len() == 2 {
                true
            } else {
                if todos_num && distintos.len() == 3 {
                    let nums: HashSet<i64> = distintos.iter().filter_map(|v| e_numero(v)).map(|x| x as i64).collect();
                    if nums.contains(&-9) && (nums.is_subset(&HashSet::from([-9, 1, 2])) || nums.is_subset(&HashSet::from([-9, 0, 1]))) {
                        return Err(format!(
                            "The outcome \"{}\" has the values {}: in PLINK files -9 means missing. Write NA for missing values (or declare the outcome quantitative if -9 is a real value).",
                            opt.outcome,
                            distintos.join(", ")
                        ));
                    }
                }
                if !todos_num {
                    return Err(format!(
                        "The outcome \"{}\" has {} distinct values among the people analysed and not all are numbers ({}): a binary outcome has exactly two values, a quantitative one only numbers.",
                        opt.outcome,
                        distintos.len(),
                        lista_curta(&distintos.iter().map(|s| s.to_string()).collect::<Vec<_>>(), 6)
                    ));
                }
                if distintos.len() < 2 {
                    return Err(format!("The outcome \"{}\" has a single value among the people analysed: there is nothing to compare", opt.outcome));
                }
                false
            }
        }
    };
    an.binario = binario;
    if binario {
        let (a, b) = (distintos[0], distintos[1]);
        let evento: &str = if let Some(cv) = opt.case_value.as_deref().map(str::trim).filter(|s| !s.is_empty()) {
            let igual = |v: &str| v == cv || (e_numero(v).is_some() && e_numero(v) == e_numero(cv));
            if igual(a) {
                a
            } else if igual(b) {
                b
            } else {
                return Err(format!("case_value \"{cv}\" is not one of the outcome's two values ({a}, {b})"));
            }
        } else if let (Some(x), Some(z)) = (e_numero(a), e_numero(b)) {
            let hi = if x < z { b } else { a };
            let (xl, xh) = (x.min(z), x.max(z));
            an.avisos.push(if xl == 0.0 && xh == 1.0 {
                format!("event = {hi} (0/1 coding)")
            } else if xl == 1.0 && xh == 2.0 {
                format!("event = {hi} (PLINK coding: 1 = control, 2 = case)")
            } else {
                format!("event = {hi}, the larger value; set case_value to choose")
            });
            hi
        } else {
            const POSITIVOS: [&str; 13] =
                ["case", "cases", "yes", "y", "true", "t", "affected", "event", "positive", "pos", "present", "1", "disease"];
            let pa = POSITIVOS.contains(&a.to_ascii_lowercase().as_str());
            let pb = POSITIVOS.contains(&b.to_ascii_lowercase().as_str());
            if pa != pb {
                let e = if pa { a } else { b };
                an.avisos.push(format!("event = \"{e}\"; set case_value to choose the other"));
                e
            } else {
                an.avisos.push(format!("event = \"{b}\", the second value in sorted order; set case_value to choose"));
                b
            }
        };
        let outro = if evento == a { b } else { a };
        an.evento = Some(evento.to_string());
        an.outro = Some(outro.to_string());
        an.y = valores.iter().map(|v| if *v == evento { 1.0 } else { 0.0 }).collect();
    } else {
        an.y = valores.iter().map(|v| e_numero(v).unwrap()).collect();
    }

    // --- the design ---
    an.z.push(vec![1.0; n]);
    an.z_nomes.push("intercept".into());
    for (nome, i) in &ic {
        let vals: Vec<&str> = linhas_ok.iter().map(|v| v[*i].as_str()).collect();
        let forcar = opt.categorical.iter().any(|c| c.trim().eq_ignore_ascii_case(nome));
        let todos: Vec<&str> = tab.linhas.iter().map(|(_, v)| v[*i].as_str()).filter(|v| !faltante(v)).collect();
        let numerico = !forcar && todos.iter().all(|v| e_numero(v).is_some());
        if numerico {
            let col: Vec<f64> = vals.iter().map(|v| e_numero(v).unwrap()).collect();
            let m = col[0];
            if col.iter().all(|x| *x == m) {
                an.avisos.push(format!("the covariate {nome} has a single value ({m}) among the people analysed and was left out"));
                continue;
            }
            an.z.push(col);
            an.z_nomes.push(nome.clone());
        } else {
            let mut niveis: Vec<&str> = vals.clone();
            niveis.sort_unstable();
            niveis.dedup();
            if niveis.len() < 2 {
                an.avisos.push(format!(
                    "the covariate {nome} has a single level ({}) among the people analysed and was left out",
                    niveis.first().copied().unwrap_or("")
                ));
                continue;
            }
            if niveis.len() > 50 && niveis.len() * 2 > n {
                return Err(format!(
                    "The covariate {nome} has {} different text values among {n} people: it looks like an identifier, not a category. Leave it out, or fix the values if it should be numeric.",
                    niveis.len()
                ));
            }
            for l in &niveis[1..] {
                an.z.push(vals.iter().map(|v| if v == l { 1.0 } else { 0.0 }).collect());
                an.z_nomes.push(format!("{nome}={l}"));
            }
            an.avisos.push(format!("{nome} is a category: reference level \"{}\", {} indicator(s)", niveis[0], niveis.len() - 1));
        }
    }
    let qr = qr_gs(&an.z, 1e-9);
    if let Some(&j) = qr.dropped.first() {
        let nome = &an.z_nomes[j];
        return Err(format!(
            "The covariate {nome} is a linear combination of the intercept and the covariates listed before it (among the {n} people analysed), so its effect cannot be told apart. Typical causes: proportions that sum to 1 (ancestry fractions: leave one out), a category with the same information as another covariate, or two copies of one column. Remove it."
        ));
    }
    if n <= an.z.len() + 1 {
        return Err(format!(
            "{n} people and {} parameters (intercept, covariates and the genotype): a regression needs more people than parameters. Use fewer covariates or more people.",
            an.z.len() + 1
        ));
    }
    Ok(an)
}

// =====================================================================
//  8) THE TESTS
// =====================================================================

/// One tested variant (one ALT allele).
#[derive(Clone, Debug, Default)]
pub struct AssocRow {
    pub id: String,
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt: String,
    /// People used (analysis set with a genotype).
    pub n: usize,
    /// People of the analysis set without a genotype here.
    pub n_missing: usize,
    pub alt_freq: f64,
    pub call_rate: f64,
    /// People with at least one copy of the allele.
    pub carriers: usize,
    pub hom_alt: usize,
    pub n_cases: usize,
    pub n_controls: usize,
    pub alt_freq_cases: f64,
    pub alt_freq_controls: f64,
    pub beta: f64,
    pub se: f64,
    pub ci_low: f64,
    pub ci_high: f64,
    pub odds_ratio: f64,
    pub or_low: f64,
    pub or_high: f64,
    /// t (linear) or Wald z (logistic, Firth).
    pub stat: f64,
    pub df: f64,
    pub p: f64,
    pub p_wald: f64,
    pub p_lrt: f64,
    /// "linear", "logistic" or "firth".
    pub method: String,
    /// "t", "wald" or "profile".
    pub ci_method: String,
    pub p_permutation: f64,
    pub p_fwer: f64,
    pub fdr: f64,
    pub note: String,
}

/// The people with a genotype at a variant and their design (shared by every
/// variant with the same missing genotypes).
struct Grupo {
    s: Vec<u32>,
    /// The covariate columns (the intercept first) restricted to s, the
    /// dependent ones dropped.
    z: Vec<Vec<f64>>,
    z_nomes_fora: Vec<String>,
    qr: Qr,
    y: Vec<f64>,
    y_res: Vec<f64>,
    yy: f64,
    /// The covariates-only logistic fit on these people.
    nulo: Option<LogitFit>,
}

fn prepara_grupo(an: &Analise, faltam: &[u32]) -> Grupo {
    let n = an.y.len();
    let mut s = Vec::with_capacity(n - faltam.len());
    let mut k = 0usize;
    for i in 0..n as u32 {
        if k < faltam.len() && faltam[k] == i {
            k += 1;
        } else {
            s.push(i);
        }
    }
    let z_todas: Vec<Vec<f64>> = an.z.iter().map(|c| s.iter().map(|&i| c[i as usize]).collect()).collect();
    let qr = qr_gs(&z_todas, 1e-9);
    let z_nomes_fora: Vec<String> = qr.dropped.iter().map(|&j| an.z_nomes[j].clone()).collect();
    let z: Vec<Vec<f64>> = qr.kept.iter().map(|&j| z_todas[j].clone()).collect();
    let y: Vec<f64> = s.iter().map(|&i| an.y[i as usize]).collect();
    let y_res = qr.residual(&y);
    let yy = dot(&y_res, &y_res);
    let nulo = if an.binario {
        let soma: f64 = y.iter().sum();
        if soma > 0.0 && soma < y.len() as f64 {
            Some(logistic_fit(&z, &y, None))
        } else {
            None
        }
    } else {
        None
    };
    Grupo { s, z, z_nomes_fora, qr, y, y_res, yy, nulo }
}

/// What makes a variant untestable.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum Motivo {
    Chamada,
    Monomorfico,
    Maf,
    SemVariacao,
    Colinear,
    Poucos,
    SemEvento,
    Ajuste,
}

impl Motivo {
    fn codigo(self) -> &'static str {
        match self {
            Motivo::Chamada => "low_call_rate",
            Motivo::Monomorfico => "monomorphic",
            Motivo::Maf => "low_maf",
            Motivo::SemVariacao => "no_variation_in_model",
            Motivo::Colinear => "collinear_with_covariates",
            Motivo::Poucos => "too_few_people",
            Motivo::SemEvento => "single_outcome_value",
            Motivo::Ajuste => "fit_failed",
        }
    }
    fn texto(self) -> &'static str {
        match self {
            Motivo::Chamada => "call rate below min_call_rate",
            Motivo::Monomorfico => "no copy of the allele (or every allele is it) among the people analysed",
            Motivo::Maf => "minor allele frequency below min_maf",
            Motivo::SemVariacao => "every person has the same genotype under the genetic model",
            Motivo::Colinear => "the genotype is a combination of the covariates among these people",
            Motivo::Poucos => "fewer people with a genotype than parameters plus one",
            Motivo::SemEvento => "everyone with a genotype has the same outcome",
            Motivo::Ajuste => "the regression could not be fitted",
        }
    }
}

struct Candidato {
    id: String,
    crom: String,
    pos: u64,
    refe: String,
    alt: String,
    cod: Vec<i8>,
    pedido: bool,
}

/// Runs one test. The row, or why the variant could not be tested.
fn testa(c: &Candidato, g: &Grupo, an: &Analise, opt: &AssocOptions, falta: usize) -> Result<AssocRow, Motivo> {
    let n_a = an.y.len();
    let mut r = AssocRow {
        id: c.id.clone(),
        chrom: c.crom.clone(),
        pos: c.pos,
        ref_allele: c.refe.clone(),
        alt: c.alt.clone(),
        n: g.s.len(),
        n_missing: falta,
        call_rate: g.s.len() as f64 / n_a as f64,
        p_permutation: f64::NAN,
        p_fwer: f64::NAN,
        fdr: f64::NAN,
        alt_freq_cases: f64::NAN,
        alt_freq_controls: f64::NAN,
        p_wald: f64::NAN,
        p_lrt: f64::NAN,
        odds_ratio: f64::NAN,
        or_low: f64::NAN,
        or_high: f64::NAN,
        df: f64::NAN,
        ..Default::default()
    };
    let (mut alelos, mut total) = (0.0, 0.0);
    let (mut ac, mut tc, mut ao, mut to) = (0.0, 0.0, 0.0, 0.0);
    for &i in &g.s {
        let cd = c.cod[i as usize];
        let (x, pl) = (copias(cd), ploidia_de(cd));
        alelos += x;
        total += pl;
        if x >= 1.0 {
            r.carriers += 1;
        }
        if cd == 2 || cd == 4 {
            r.hom_alt += 1;
        }
        if an.binario {
            if an.y[i as usize] == 1.0 {
                ac += x;
                tc += pl;
                r.n_cases += 1;
            } else {
                ao += x;
                to += pl;
                r.n_controls += 1;
            }
        }
    }
    r.alt_freq = if total > 0.0 { alelos / total } else { f64::NAN };
    if an.binario {
        r.alt_freq_cases = if tc > 0.0 { ac / tc } else { f64::NAN };
        r.alt_freq_controls = if to > 0.0 { ao / to } else { f64::NAN };
    }
    if r.call_rate < opt.min_call_rate {
        return Err(Motivo::Chamada);
    }
    if !(r.alt_freq > 0.0 && r.alt_freq < 1.0) {
        return Err(Motivo::Monomorfico);
    }
    let maf = r.alt_freq.min(1.0 - r.alt_freq);
    if maf < opt.min_maf {
        return Err(Motivo::Maf);
    }
    let x: Vec<f64> = g.s.iter().map(|&i| codifica(c.cod[i as usize], opt.model)).collect();
    let x0 = x[0];
    if x.iter().all(|v| *v == x0) {
        return Err(Motivo::SemVariacao);
    }
    let k = g.z.len();
    let ns = g.s.len();
    if ns < k + 2 {
        return Err(Motivo::Poucos);
    }
    let x_res = g.qr.residual(&x);
    let xx = dot(&x_res, &x_res);
    let media = x.iter().sum::<f64>() / ns as f64;
    let xc: f64 = x.iter().map(|v| (v - media) * (v - media)).sum();
    if !(xx > 1e-9 * xc) {
        return Err(Motivo::Colinear);
    }
    let mut notas: Vec<String> = Vec::new();
    if !g.z_nomes_fora.is_empty() {
        notas.push(format!("covariate(s) with no variation among these people left out: {}", g.z_nomes_fora.join(", ")));
    }
    let n1 = x.iter().filter(|v| **v != 0.0).count();
    let n0 = ns - n1;
    if n1.min(n0) < 5 {
        if n1 < n0 {
            notas.push(format!("only {n1} people with the coded genotype"));
        } else {
            notas.push(format!("only {n0} people without the coded genotype"));
        }
    }
    if !an.binario {
        let df = (ns - k - 1) as f64;
        let c_ = dot(&x_res, &g.y_res);
        let beta = c_ / xx;
        let rss = (g.yy - c_ * c_ / xx).max(0.0);
        let s2 = rss / df;
        let se = (s2 / xx).sqrt();
        let t = beta / se;
        r.beta = beta;
        r.se = se;
        r.stat = t;
        r.df = df;
        r.p = t_two_sided(t, df);
        let tc = t_critical(df, 0.05);
        r.ci_low = beta - tc * se;
        r.ci_high = beta + tc * se;
        r.method = "linear".into();
        r.ci_method = "t".into();
        if ns < 30 {
            notas.push(format!("{ns} people: the t test assumes normal residuals"));
        }
    } else {
        let eventos = g.y.iter().sum::<f64>() as usize;
        if eventos == 0 || eventos == ns {
            return Err(Motivo::SemEvento);
        }
        if eventos.min(ns - eventos) < 10 {
            notas.push(format!("only {} people with the rarer outcome", eventos.min(ns - eventos)));
        }
        let mut cols = g.z.clone();
        cols.push(x);
        let j = cols.len() - 1;
        let mut usa_firth = opt.firth == FirthPolicy::Always;
        let mut plano: Option<LogitFit> = None;
        if !usa_firth {
            let f = logistic_fit(&cols, &g.y, None);
            if f.separated || !f.converged {
                if opt.firth == FirthPolicy::Auto {
                    usa_firth = true;
                    notas.push(if f.converged { "separation: Firth".into() } else { "plain fit did not converge: Firth".into() });
                } else {
                    notas.push("separation (infinite estimate): only the likelihood-ratio p is reported; use firth".into());
                }
            }
            plano = Some(f);
        }
        if usa_firth {
            let Some(f) = firth_fit(&cols, &g.y, None, None) else {
                return Err(Motivo::Ajuste);
            };
            if !f.converged {
                notas.push("Firth fit did not converge fully".into());
            }
            let (_, p) = firth_lrt(&cols, &g.y, &f, j);
            let (lo, hi) = firth_profile_ci(&cols, &g.y, &f, j, CHI2_1_95);
            r.beta = f.beta[j];
            r.se = f.se[j];
            r.stat = r.beta / r.se;
            r.p_wald = p_two_sided_z(r.stat);
            r.p_lrt = p;
            r.p = p;
            r.ci_low = lo;
            r.ci_high = hi;
            r.method = "firth".into();
            r.ci_method = "profile".into();
        } else {
            let f = plano.unwrap();
            let l0 = g.nulo.as_ref().map(|z| z.loglik).unwrap_or(f64::NAN);
            let lr = (2.0 * (f.loglik - l0)).max(0.0);
            r.p_lrt = chi2_1_upper(lr);
            r.p = r.p_lrt;
            r.method = "logistic".into();
            if f.separated {
                r.beta = f64::NAN;
                r.se = f64::NAN;
                r.stat = f64::NAN;
                r.ci_low = f64::NAN;
                r.ci_high = f64::NAN;
                r.ci_method = "none".into();
            } else {
                r.beta = f.beta[j];
                r.se = f.se[j];
                r.stat = r.beta / r.se;
                r.p_wald = p_two_sided_z(r.stat);
                r.ci_low = r.beta - Z975 * r.se;
                r.ci_high = r.beta + Z975 * r.se;
                r.ci_method = "wald".into();
            }
        }
        r.odds_ratio = r.beta.exp();
        r.or_low = r.ci_low.exp();
        r.or_high = r.ci_high.exp();
    }
    r.note = notas.join("; ");
    Ok(r)
}

// =====================================================================
//  9) PERMUTATIONS
// =====================================================================

/// The statistic of every variant of a group in one replicate.
struct GrupoPerm {
    s: Vec<u32>,
    z: Vec<Vec<f64>>,
    qr: Qr,
    /// Index in the rows, coded genotype on s, its residual (linear), |.|^2.
    vars: Vec<(usize, Vec<f64>, Vec<f64>, f64)>,
    df: f64,
    nulo_beta: Option<Vec<f64>>,
}

/// Statistic of each variant (t^2, or the score chi-square) given an outcome
/// vector on the analysis set.
fn estatisticas(g: &GrupoPerm, y_a: &[f64], binario: bool, saida: &mut Vec<(usize, f64)>) {
    let y: Vec<f64> = g.s.iter().map(|&i| y_a[i as usize]).collect();
    if !binario {
        let r = g.qr.residual(&y);
        let rr = dot(&r, &r);
        for (v, _, xres, xx) in &g.vars {
            let c = dot(xres, &r);
            let den = xx * rr - c * c;
            let s = if den > 0.0 { c * c * g.df / den } else { f64::INFINITY };
            saida.push((*v, s));
        }
        return;
    }
    let soma: f64 = y.iter().sum();
    if soma <= 0.0 || soma >= y.len() as f64 {
        for (v, ..) in &g.vars {
            saida.push((*v, 0.0));
        }
        return;
    }
    let nulo = logistic_fit(&g.z, &y, g.nulo_beta.as_deref());
    let sw: Vec<f64> = nulo.fitted.iter().map(|p| (p * (1.0 - p)).sqrt()).collect();
    let zw: Vec<Vec<f64>> = g.z.iter().map(|c| c.iter().zip(&sw).map(|(a, b)| a * b).collect()).collect();
    let qw = qr_gs(&zw, 1e-12);
    let resid: Vec<f64> = y.iter().zip(&nulo.fitted).map(|(a, b)| a - b).collect();
    for (v, x, ..) in &g.vars {
        let u = dot(x, &resid);
        let xw: Vec<f64> = x.iter().zip(&sw).map(|(a, b)| a * b).collect();
        let xr = qw.residual(&xw);
        let vv = dot(&xr, &xr);
        saida.push((*v, if vv > 1e-12 * dot(&xw, &xw) { u * u / vv } else { 0.0 }));
    }
}

fn p_da_estatistica(s: f64, df: f64, binario: bool) -> f64 {
    if binario {
        chi2_1_upper(s)
    } else {
        t_two_sided(s.max(0.0).sqrt(), df)
    }
}

struct Permutacoes {
    p_emp: Vec<f64>,
    p_fwer: Vec<f64>,
    esquema: String,
}

fn permuta(an: &Analise, linhas: &[AssocRow], cods: &[Vec<i8>], opt: &AssocOptions) -> Result<Permutacoes, String> {
    let b = opt.permutations;
    let nv = linhas.len();
    let n_a = an.y.len();
    // groups by missing genotypes
    let mut mapa: HashMap<Vec<u32>, usize> = HashMap::new();
    let mut grupos: Vec<GrupoPerm> = Vec::new();
    for (v, cod) in cods.iter().enumerate() {
        let faltam: Vec<u32> = (0..n_a as u32).filter(|&i| cod[i as usize] < 0).collect();
        let gi = match mapa.get(&faltam) {
            Some(&g) => g,
            None => {
                let g = prepara_grupo(an, &faltam);
                let df = (g.s.len() - g.z.len() - 1) as f64;
                grupos.push(GrupoPerm {
                    s: g.s,
                    z: g.z,
                    qr: g.qr,
                    vars: Vec::new(),
                    df,
                    nulo_beta: g.nulo.map(|f| f.beta),
                });
                mapa.insert(faltam, grupos.len() - 1);
                grupos.len() - 1
            }
        };
        let g = &mut grupos[gi];
        let x: Vec<f64> = g.s.iter().map(|&i| codifica(cod[i as usize], opt.model)).collect();
        let xr = g.qr.residual(&x);
        let xx = dot(&xr, &xr);
        g.vars.push((v, x, xr, xx));
    }
    // budget
    let mut custo = 0f64;
    for g in &grupos {
        let n = g.s.len() as f64;
        let k = g.z.len() as f64;
        custo += n * k * 4.0 + g.vars.len() as f64 * n * if an.binario { 2.0 * k + 4.0 } else { 2.0 };
        if an.binario {
            custo += 8.0 * n * k * k;
        }
    }
    custo *= b as f64;
    // ~1e9 of these per second on 8 threads: 1e12 is about 20 minutes
    if custo > 1e12 {
        return Err(format!(
            "{} permutations of {} variants with {} people would take a long time (about {:.0e} operations, tens of minutes or more). Use fewer permutations, or permute only the variants of interest (ids or regions).",
            b, nv, n_a, custo
        ));
    }
    // observed
    let mut obs = vec![f64::NAN; nv];
    let mut obs_p = vec![f64::NAN; nv];
    let mut tmp = Vec::new();
    for g in &grupos {
        tmp.clear();
        estatisticas(g, &an.y, an.binario, &mut tmp);
        for &(v, s) in &tmp {
            obs[v] = s;
            obs_p[v] = p_da_estatistica(s, g.df, an.binario);
        }
    }
    // the null model on the analysis set
    let semcov = an.z.len() == 1;
    let qa = qr_gs(&an.z, 1e-9);
    let (ajuste_a, resid_a): (Vec<f64>, Vec<f64>) = if !an.binario {
        let r = qa.residual(&an.y);
        (an.y.iter().zip(&r).map(|(y, e)| y - e).collect(), r)
    } else {
        let f = logistic_fit(&an.z, &an.y, None);
        (f.fitted, Vec::new())
    };
    let esquema = if !an.binario {
        if semcov {
            "outcome permuted (exact)".to_string()
        } else {
            "Freedman-Lane: residuals of the covariates-only linear model permuted and added to its fitted values".to_string()
        }
    } else if semcov {
        "outcome labels permuted (exact; the number of events is kept)".to_string()
    } else {
        "parametric bootstrap under the null: each replicate draws the outcome from the covariates-only logistic model fitted to the data (no exact permutation keeps the covariate-outcome relation)".to_string()
    };
    let tol = 1e-9;
    let (contagens, mut minps) = (0..b as u64)
        .into_par_iter()
        .fold(
            || (vec![0u32; nv], Vec::<f64>::new(), Vec::<(usize, f64)>::new()),
            |(mut cont, mut mins, mut buf), rep| {
                let mut rng = Rng::stream(opt.seed, rep);
                let ystar: Vec<f64> = if !an.binario {
                    let mut idx: Vec<usize> = (0..n_a).collect();
                    rng.shuffle(&mut idx);
                    (0..n_a).map(|i| ajuste_a[i] + resid_a[idx[i]]).collect()
                } else if semcov {
                    let mut v = an.y.clone();
                    rng.shuffle(&mut v);
                    v
                } else {
                    ajuste_a.iter().map(|p| if rng.unif() < *p { 1.0 } else { 0.0 }).collect()
                };
                let mut minp = 1.0f64;
                for g in &grupos {
                    buf.clear();
                    estatisticas(g, &ystar, an.binario, &mut buf);
                    let mut maior = f64::NEG_INFINITY;
                    for &(v, s) in &buf {
                        if s >= obs[v] * (1.0 - tol) - 1e-12 {
                            cont[v] += 1;
                        }
                        if s > maior {
                            maior = s;
                        }
                    }
                    if maior > f64::NEG_INFINITY {
                        minp = minp.min(p_da_estatistica(maior, g.df, an.binario));
                    }
                }
                mins.push(minp);
                (cont, mins, buf)
            },
        )
        .map(|(c, m, _)| (c, m))
        .reduce(
            || (vec![0u32; nv], Vec::new()),
            |(mut a, mut ma), (c, mc)| {
                for (x, y) in a.iter_mut().zip(c) {
                    *x += y;
                }
                ma.extend(mc);
                (a, ma)
            },
        );
    minps.sort_by(|a, b| a.total_cmp(b));
    let bf = b as f64;
    let p_emp: Vec<f64> = contagens.iter().map(|&c| (c as f64 + 1.0) / (bf + 1.0)).collect();
    let p_fwer: Vec<f64> = obs_p
        .iter()
        .map(|&p| {
            let lim = p * (1.0 + tol) + 1e-300;
            let c = minps.partition_point(|m| *m <= lim);
            (c as f64 + 1.0) / (bf + 1.0)
        })
        .collect();
    Ok(Permutacoes { p_emp, p_fwer, esquema })
}

// =====================================================================
//  10) THE REPORT
// =====================================================================

#[derive(Clone, Debug, Default)]
pub struct AssocReport {
    pub vcfs: Vec<String>,
    pub phenotypes: String,
    pub outcome: String,
    /// "binary" or "quantitative".
    pub outcome_type: String,
    pub event_value: Option<String>,
    pub other_value: Option<String>,
    pub model: String,
    pub covariates: Vec<String>,
    /// The design columns besides the genotype (intercept, covariates, indicators).
    pub design: Vec<String>,
    pub firth: String,
    pub permutations: usize,
    pub permutation_scheme: String,
    pub seed: u64,
    pub absent_is_ref: bool,
    /// "indexed" (only the .tbi blocks of the regions) or "streamed".
    pub reading: String,
    // people
    pub vcf_samples: usize,
    pub table_rows: usize,
    pub analysed: usize,
    pub not_in_table: usize,
    pub not_in_vcf: usize,
    pub not_in_vcf_examples: Vec<String>,
    pub missing_outcome: usize,
    pub missing_covariates: Vec<(String, usize)>,
    pub cases: usize,
    pub controls: usize,
    // variants
    pub records_read: u64,
    pub malformed_records: u64,
    pub skipped_sequences: u64,
    pub symbolic_alleles: usize,
    pub alleles_considered: usize,
    pub tested: usize,
    pub skipped: Vec<(String, usize)>,
    pub skipped_examples: Vec<(String, String)>,
    pub not_found: Vec<String>,
    pub haploid_calls: usize,
    pub rows: Vec<AssocRow>,
    // checks and power
    pub residual_skewness: f64,
    pub residual_excess_kurtosis: f64,
    pub power_maf: f64,
    pub power_alpha: f64,
    /// The smallest odds ratio (binary) or beta in SD units (quantitative)
    /// detectable with 80% power at `power_maf` with these people.
    pub detectable_effect: f64,
    pub output: Option<String>,
    pub warnings: Vec<String>,
}

impl AssocReport {
    fn linha_json(r: &AssocRow) -> String {
        format!(
            "{{\"id\":{},\"chrom\":{},\"pos\":{},\"ref\":{},\"alt\":{},\"n\":{},\"n_missing\":{},\"alt_freq\":{},\"call_rate\":{},\"carriers\":{},\"hom_alt\":{},\"n_cases\":{},\"n_controls\":{},\"alt_freq_cases\":{},\"alt_freq_controls\":{},\"beta\":{},\"se\":{},\"ci_low\":{},\"ci_high\":{},\"or\":{},\"or_low\":{},\"or_high\":{},\"stat\":{},\"df\":{},\"p\":{},\"p_wald\":{},\"p_lrt\":{},\"method\":{},\"ci_method\":{},\"p_permutation\":{},\"p_fwer\":{},\"fdr\":{},\"note\":{}}}",
            json_texto(&r.id),
            json_texto(&r.chrom),
            r.pos,
            json_texto(&r.ref_allele),
            json_texto(&r.alt),
            r.n,
            r.n_missing,
            jnum(r.alt_freq),
            jnum(r.call_rate),
            r.carriers,
            r.hom_alt,
            r.n_cases,
            r.n_controls,
            jnum(r.alt_freq_cases),
            jnum(r.alt_freq_controls),
            jnum(r.beta),
            jnum(r.se),
            jnum(r.ci_low),
            jnum(r.ci_high),
            jnum(r.odds_ratio),
            jnum(r.or_low),
            jnum(r.or_high),
            jnum(r.stat),
            jnum(r.df),
            jnum(r.p),
            jnum(r.p_wald),
            jnum(r.p_lrt),
            json_texto(&r.method),
            json_texto(&r.ci_method),
            jnum(r.p_permutation),
            jnum(r.p_fwer),
            jnum(r.fdr),
            json_texto(&r.note)
        )
    }

    /// The rows by p (smallest first), NaN last, ties in genomic order.
    pub fn order_by_p(&self) -> Vec<usize> {
        let mut idx: Vec<usize> = (0..self.rows.len()).collect();
        idx.sort_by(|a, b| {
            let (x, y) = (self.rows[*a].p, self.rows[*b].p);
            match (x.is_nan(), y.is_nan()) {
                (false, false) => x.total_cmp(&y),
                (true, false) => std::cmp::Ordering::Greater,
                (false, true) => std::cmp::Ordering::Less,
                (true, true) => std::cmp::Ordering::Equal,
            }
        });
        idx
    }

    pub fn json(&self) -> String {
        const MAX_LINHAS: usize = 50_000;
        let (linhas, truncado): (Vec<usize>, bool) = if self.rows.len() > MAX_LINHAS {
            let mut v: Vec<usize> = self.order_by_p().into_iter().take(MAX_LINHAS).collect();
            v.sort_unstable();
            (v, true)
        } else {
            ((0..self.rows.len()).collect(), false)
        };
        let rows = linhas.iter().map(|&i| Self::linha_json(&self.rows[i])).collect::<Vec<_>>().join(",");
        let top = self
            .order_by_p()
            .into_iter()
            .take(10)
            .map(|i| {
                let r = &self.rows[i];
                format!("{{\"id\":{},\"p\":{},\"fdr\":{}}}", json_texto(&r.id), jnum(r.p), jnum(r.fdr))
            })
            .collect::<Vec<_>>()
            .join(",");
        let pares = |v: &[(String, usize)]| -> String {
            format!("{{{}}}", v.iter().map(|(k, n)| format!("{}:{n}", json_texto(k))).collect::<Vec<_>>().join(","))
        };
        let exemplos = format!(
            "[{}]",
            self.skipped_examples
                .iter()
                .map(|(id, m)| format!("{{\"id\":{},\"reason\":{}}}", json_texto(id), json_texto(m)))
                .collect::<Vec<_>>()
                .join(",")
        );
        let metodos = {
            let mut m: Vec<(String, usize)> = Vec::new();
            for r in &self.rows {
                match m.iter_mut().find(|(k, _)| *k == r.method) {
                    Some(x) => x.1 += 1,
                    None => m.push((r.method.clone(), 1)),
                }
            }
            pares(&m)
        };
        let signif = self.rows.iter().filter(|r| r.fdr < 0.05).count();
        format!(
            "{{\"vcfs\":{},\"phenotypes\":{},\"outcome\":{},\"outcome_type\":{},\"event_value\":{},\"other_value\":{},\"model\":{},\"covariates\":{},\"design\":{},\"firth\":{},\"permutations\":{},\"permutation_scheme\":{},\"seed\":{},\"absent_is_ref\":{},\"reading\":{},\"people\":{{\"vcf_samples\":{},\"table_rows\":{},\"analysed\":{},\"not_in_table\":{},\"not_in_vcf\":{},\"not_in_vcf_examples\":{},\"missing_outcome\":{},\"missing_covariates\":{},\"cases\":{},\"controls\":{}}},\"variants\":{{\"records_read\":{},\"malformed_records\":{},\"skipped_sequences\":{},\"symbolic_alleles\":{},\"alleles_considered\":{},\"tested\":{},\"skipped\":{},\"skipped_examples\":{},\"not_found\":{},\"haploid_calls\":{},\"methods\":{},\"fdr_below_0_05\":{}}},\"residuals\":{{\"skewness\":{},\"excess_kurtosis\":{}}},\"power\":{{\"maf\":{},\"alpha_per_test\":{},\"target_power\":0.8,\"detectable_effect\":{},\"effect_scale\":{}}},\"top\":[{}],\"rows_truncated\":{},\"rows\":[{}],\"output\":{},\"warnings\":{},\"text\":{}}}",
            jtextos(&self.vcfs),
            json_texto(&self.phenotypes),
            json_texto(&self.outcome),
            json_texto(&self.outcome_type),
            jopt(self.event_value.as_deref()),
            jopt(self.other_value.as_deref()),
            json_texto(&self.model),
            jtextos(&self.covariates),
            jtextos(&self.design),
            json_texto(&self.firth),
            self.permutations,
            json_texto(&self.permutation_scheme),
            self.seed,
            self.absent_is_ref,
            json_texto(&self.reading),
            self.vcf_samples,
            self.table_rows,
            self.analysed,
            self.not_in_table,
            self.not_in_vcf,
            jtextos(&self.not_in_vcf_examples),
            self.missing_outcome,
            pares(&self.missing_covariates),
            self.cases,
            self.controls,
            self.records_read,
            self.malformed_records,
            self.skipped_sequences,
            self.symbolic_alleles,
            self.alleles_considered,
            self.tested,
            pares(&self.skipped),
            exemplos,
            jtextos(&self.not_found),
            self.haploid_calls,
            metodos,
            signif,
            jnum(self.residual_skewness),
            jnum(self.residual_excess_kurtosis),
            jnum(self.power_maf),
            jnum(self.power_alpha),
            jnum(self.detectable_effect),
            json_texto(if self.outcome_type == "binary" { "odds_ratio" } else { "beta_sd" }),
            top,
            truncado,
            rows,
            jopt(self.output.as_deref()),
            jtextos(&self.warnings),
            json_texto(&self.text())
        )
    }

    pub fn text(&self) -> String {
        let bin = self.outcome_type == "binary";
        let mut o = String::new();
        let ajuste = if self.covariates.is_empty() { "no covariates".to_string() } else { format!("adjusted for {}", self.covariates.join(", ")) };
        if bin {
            o.push_str(&format!(
                "Association of {} (binary: {} with \"{}\", {} without) with genotypes; {} model, {}.\n",
                self.outcome,
                self.cases,
                self.event_value.as_deref().unwrap_or("?"),
                self.controls,
                self.model,
                ajuste
            ));
        } else {
            o.push_str(&format!("Association of {} (quantitative) with genotypes; {} model, {}.\n", self.outcome, self.model, ajuste));
        }
        let mut pessoas = format!("  people: {} analysed", milhar(self.analysed));
        let mut fora = Vec::new();
        if self.not_in_vcf > 0 {
            fora.push(format!("{} of the table not in any VCF", self.not_in_vcf));
        }
        if self.not_in_table > 0 {
            fora.push(format!("{} VCF sample(s) not in the table", self.not_in_table));
        }
        if self.missing_outcome > 0 {
            fora.push(format!("{} without the outcome", self.missing_outcome));
        }
        for (c, k) in &self.missing_covariates {
            if *k > 0 {
                fora.push(format!("{k} without {c}"));
            }
        }
        if !fora.is_empty() {
            pessoas.push_str(&format!(" (left out: {})", fora.join("; ")));
        }
        o.push_str(&pessoas);
        o.push_str("; each test uses those with a genotype at the variant.\n");
        let mut m: Vec<(String, usize)> = Vec::new();
        for r in &self.rows {
            match m.iter_mut().find(|(k, _)| *k == r.method) {
                Some(x) => x.1 += 1,
                None => m.push((r.method.clone(), 1)),
            }
        }
        let metodos = m.iter().map(|(k, n)| format!("{k} {n}")).collect::<Vec<_>>().join(", ");
        o.push_str(&format!(
            "  variants: {} allele(s) tested{}",
            milhar(self.tested),
            if metodos.is_empty() { String::new() } else { format!(" ({metodos})") }
        ));
        let pul: Vec<String> = self.skipped.iter().filter(|(_, k)| *k > 0).map(|(c, k)| format!("{k} {c}")).collect();
        if !pul.is_empty() {
            o.push_str(&format!("; not tested: {}", pul.join(", ")));
        }
        o.push_str(".\n");
        if !self.not_found.is_empty() {
            o.push_str(&format!("  not found in the VCF(s): {}\n", lista_curta(&self.not_found, 10)));
        }
        if bin {
            o.push_str("  p: likelihood-ratio test (Firth: penalized likelihood ratio, as logistf); p_wald also given; OR = odds ratio per copy (additive) or for the coded genotype.\n");
        } else {
            o.push_str("  p: t test of the genotype coefficient (OLS); beta in outcome units per copy (additive) or for the coded genotype.\n");
        }
        let ordem = self.order_by_p();
        if !ordem.is_empty() {
            o.push_str("  smallest p-values:\n");
            for &i in ordem.iter().take(5) {
                let r = &self.rows[i];
                let efeito = if bin {
                    format!("OR {} [{}, {}]", fmt_g(r.odds_ratio, 2), fmt_g(r.or_low, 2), fmt_g(r.or_high, 2))
                } else {
                    format!("beta {} [{}, {}]", fmt_g(r.beta, 3), fmt_g(r.ci_low, 3), fmt_g(r.ci_high, 3))
                };
                let mut extra = format!("FDR {}", fmt_p(r.fdr));
                if !r.p_permutation.is_nan() {
                    extra.push_str(&format!(", permutation p {}, family-wise p {}", fmt_p(r.p_permutation), fmt_p(r.p_fwer)));
                }
                o.push_str(&format!(
                    "    {} ({}:{} {}>{}), n {}, alt freq {:.3}: {}, p {} ({}), {}{}\n",
                    r.id,
                    r.chrom,
                    r.pos,
                    r.ref_allele,
                    r.alt,
                    r.n,
                    r.alt_freq,
                    efeito,
                    fmt_p(r.p),
                    r.method,
                    extra,
                    if r.note.is_empty() { String::new() } else { format!(" -- {}", r.note) }
                ));
            }
            let sig = self.rows.iter().filter(|r| r.fdr < 0.05).count();
            o.push_str(&format!("  {} allele(s) with FDR < 0.05 (Benjamini-Hochberg over the {} tests).\n", sig, self.tested));
        }
        if self.permutations > 0 {
            o.push_str(&format!(
                "  permutations: {} replicates (seed {}); {}. Empirical p = (r+1)/(B+1); family-wise p by single-step max-T over all tested alleles.\n",
                self.permutations, self.seed, self.permutation_scheme
            ));
        }
        if self.detectable_effect.is_finite() {
            if bin {
                o.push_str(&format!(
                    "  power: with these people, an allele at frequency {:.2} needs an odds ratio of at least {} to be found with 80% power at alpha {} (0.05 over {} tests).\n",
                    self.power_maf,
                    fmt_g(self.detectable_effect, 2),
                    fmt_p(self.power_alpha),
                    self.tested.max(1)
                ));
            } else {
                o.push_str(&format!(
                    "  power: with these people, an allele at frequency {:.2} needs an effect of at least {} SD of the outcome per copy to be found with 80% power at alpha {} (0.05 over {} tests).\n",
                    self.power_maf,
                    fmt_g(self.detectable_effect, 3),
                    fmt_p(self.power_alpha),
                    self.tested.max(1)
                ));
            }
        } else if self.tested > 0 && self.power_maf.is_finite() {
            o.push_str("  power: with these people no effect of a realistic size reaches 80% power (association_power(...) gives the details).\n");
        }
        if let Some(p) = &self.output {
            o.push_str(&format!("  results written to {p}\n"));
        }
        for w in &self.warnings {
            o.push_str(&format!("  note: {w}\n"));
        }
        o.push_str("  An association is not causation, and this is software for teaching and research, not a diagnosis.\n");
        o
    }

    fn grava(&mut self, saida: &str) -> Result<(), String> {
        let mut w = Escritor::cria(saida)?;
        let mut guarda = Saida { caminho: saida.to_string(), armada: true };
        let erro = |e: String| format!("Error writing {saida}: {e}");
        w.linha("id\tchrom\tpos\tref\talt\tn\talt_freq\tcall_rate\tbeta\tse\tci_low\tci_high\tor\tor_low\tor_high\tp\tmethod\tp_permutation\tp_fwer\tfdr\tn_missing\tcarriers\thom_alt\tn_cases\tn_controls\talt_freq_cases\talt_freq_controls\tstat\tdf\tp_wald\tp_lrt\tci_method\tnote")
            .map_err(erro)?;
        for r in &self.rows {
            w.linha(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.id,
                r.chrom,
                r.pos,
                r.ref_allele,
                r.alt,
                r.n,
                num(r.alt_freq),
                num(r.call_rate),
                num(r.beta),
                num(r.se),
                num(r.ci_low),
                num(r.ci_high),
                num(r.odds_ratio),
                num(r.or_low),
                num(r.or_high),
                num(r.p),
                r.method,
                num(r.p_permutation),
                num(r.p_fwer),
                num(r.fdr),
                r.n_missing,
                r.carriers,
                r.hom_alt,
                if self.outcome_type == "binary" { r.n_cases.to_string() } else { "NA".into() },
                if self.outcome_type == "binary" { r.n_controls.to_string() } else { "NA".into() },
                num(r.alt_freq_cases),
                num(r.alt_freq_controls),
                num(r.stat),
                num(r.df),
                num(r.p_wald),
                num(r.p_lrt),
                r.ci_method,
                if r.note.is_empty() { "." } else { r.note.as_str() }
            ))
            .map_err(erro)?;
        }
        w.termina().map_err(erro)?;
        guarda.armada = false;
        self.output = Some(saida.to_string());
        Ok(())
    }
}

// =====================================================================
//  11) THE DRIVER
// =====================================================================

struct Motor<'a> {
    an: &'a Analise,
    arqs: &'a [Arquivo],
    opt: &'a AssocOptions,
    sel: Selecao,
    lote: Vec<Candidato>,
    linhas: Vec<AssocRow>,
    cods: Vec<Vec<i8>>,
    grupos: HashMap<Vec<u32>, (Grupo, u64)>,
    relogio: u64,
    pulados: HashMap<Motivo, usize>,
    exemplos: Vec<(String, String)>,
    simbolicos: usize,
    alelos: usize,
    haploides: usize,
}

const LOTE: usize = 512;

impl<'a> Motor<'a> {
    /// The records of one position (from one or several files).
    fn posicao(&mut self, regs: &[Registro]) -> Result<(), String> {
        let canon = regs[0].canon.clone();
        let pos = regs[0].pos;
        let sel = &self.sel;
        let pos_ok = sel.todas || sel.em_regiao(&canon, pos) || sel.pos.contains_key(&(canon.clone(), pos));
        let tem_exato = !sel.pos_exatos.is_empty() && sel.pos_exatos.contains(&(canon.clone(), pos));
        let tem_rs = !sel.rs.is_empty() && regs.iter().any(|r| r.ids().any(|i| sel.rs.contains_key(&i.to_ascii_lowercase())));
        if !pos_ok && !tem_exato && !tem_rs {
            return Ok(());
        }
        if let Some(v) = self.sel.pos.get(&(canon.clone(), pos)) {
            for &i in v {
                self.sel.achados[i] = true;
            }
        }
        // the alleles at this position
        let mut chaves: Vec<(String, String)> = Vec::new();
        for r in regs {
            for a in &r.alts {
                if !alelo_util(a) || a == &r.refe {
                    if a != "." && a != &r.refe {
                        self.simbolicos += 1;
                    }
                    continue;
                }
                if !chaves.iter().any(|(x, y)| *x == r.refe && y == a) {
                    chaves.push((r.refe.clone(), a.clone()));
                }
            }
        }
        // which are selected
        let mut escolhidas: Vec<(String, String, bool)> = Vec::new();
        for (rf, al) in chaves {
            let mut pedido = false;
            let mut ok = pos_ok;
            if let Some(v) = self.sel.exatos.get(&(canon.clone(), pos, rf.clone(), al.clone())) {
                ok = true;
                pedido = true;
                for &i in v.clone().iter() {
                    self.sel.achados[i] = true;
                }
            }
            if !self.sel.rs.is_empty() {
                for r in regs.iter().filter(|r| r.refe == rf && r.alts.iter().any(|a| *a == al)) {
                    for i in r.ids() {
                        if let Some(v) = self.sel.rs.get(&i.to_ascii_lowercase()) {
                            ok = true;
                            pedido = true;
                            for &k in v.clone().iter() {
                                self.sel.achados[k] = true;
                            }
                        }
                    }
                }
            }
            if self.sel.pos.contains_key(&(canon.clone(), pos)) {
                pedido = true;
            }
            if ok {
                escolhidas.push((rf, al, pedido));
            }
        }
        if escolhidas.is_empty() {
            return Ok(());
        }
        // genotypes of every record at this position, once
        let gts: Vec<Vec<Gt>> = regs.iter().map(|r| r.gts(&self.arqs[r.arq])).collect();
        let n_a = self.an.y.len();
        let padrao: i8 = if self.opt.absent_is_ref { 0 } else { -1 };
        for (rf, al, pedido) in escolhidas {
            self.alelos += 1;
            let mut cod = vec![padrao; n_a];
            for (fi, arq) in self.arqs.iter().enumerate() {
                if arq.servidos.is_empty() {
                    continue;
                }
                let mut listado: Option<(usize, u16)> = None;
                let mut mesma_ref: Option<usize> = None;
                for (k, r) in regs.iter().enumerate() {
                    if r.arq != fi || r.refe != rf {
                        continue;
                    }
                    if mesma_ref.is_none() {
                        mesma_ref = Some(k);
                    }
                    if listado.is_none() {
                        if let Some(i) = r.alts.iter().position(|a| *a == al) {
                            listado = Some((k, (i + 1) as u16));
                        }
                    }
                }
                if let Some((k, idx)) = listado {
                    let ok = regs[k].passa || !self.opt.pass_only;
                    for (s, &(_, pessoa)) in arq.servidos.iter().enumerate() {
                        cod[pessoa as usize] = if ok { codigo(gts[k][s], idx) } else { -1 };
                    }
                } else if let Some(k) = mesma_ref {
                    let ok = regs[k].passa || !self.opt.pass_only;
                    for (s, &(_, pessoa)) in arq.servidos.iter().enumerate() {
                        let g = gts[k][s];
                        cod[pessoa as usize] = if !ok {
                            -1
                        } else {
                            match g.ploidia {
                                2 => 0,
                                1 => 3,
                                _ => -1,
                            }
                        };
                    }
                }
            }
            let id = regs
                .iter()
                .filter(|r| r.refe == rf && r.alts.iter().any(|a| *a == al))
                .flat_map(|r| r.ids())
                .next()
                .map(str::to_string)
                .unwrap_or_else(|| format!("{}:{}:{}:{}", regs[0].bruto(), pos, rf, al));
            self.lote.push(Candidato { id, crom: regs[0].bruto().to_string(), pos, refe: rf, alt: al, cod, pedido });
        }
        if self.lote.len() >= LOTE {
            self.testa_lote();
        }
        Ok(())
    }

    fn testa_lote(&mut self) {
        if self.lote.is_empty() {
            return;
        }
        let lote = std::mem::take(&mut self.lote);
        let n_a = self.an.y.len();
        // groups, prepared before the parallel tests
        let padroes: Vec<Vec<u32>> = lote.iter().map(|c| (0..n_a as u32).filter(|&i| c.cod[i as usize] < 0).collect()).collect();
        for p in &padroes {
            self.relogio += 1;
            let t = self.relogio;
            if let Some(g) = self.grupos.get_mut(p) {
                g.1 = t;
            } else {
                let g = prepara_grupo(self.an, p);
                self.grupos.insert(p.clone(), (g, t));
            }
        }
        let grupos = &self.grupos;
        let (an, opt) = (self.an, self.opt);
        let res: Vec<Result<AssocRow, Motivo>> = lote
            .par_iter()
            .zip(padroes.par_iter())
            .map(|(c, p)| testa(c, &grupos[p].0, an, opt, p.len()))
            .collect();
        for (c, r) in lote.into_iter().zip(res) {
            match r {
                Ok(row) => {
                    self.haploides += c.cod.iter().filter(|x| **x >= 3).count();
                    self.linhas.push(row);
                    if self.opt.permutations > 0 {
                        self.cods.push(c.cod);
                    }
                }
                Err(m) => {
                    *self.pulados.entry(m).or_insert(0) += 1;
                    if c.pedido || self.exemplos.len() < 20 {
                        if self.exemplos.len() < 200 {
                            self.exemplos.push((c.id.clone(), m.texto().to_string()));
                        }
                    }
                }
            }
        }
        // keep the cache small: the most recently used groups
        if self.grupos.len() > 64 {
            let mut v: Vec<(Vec<u32>, u64)> = self.grupos.iter().map(|(k, g)| (k.clone(), g.1)).collect();
            v.sort_by(|a, b| b.1.cmp(&a.1));
            for (k, _) in v.into_iter().skip(32) {
                self.grupos.remove(&k);
            }
        }
    }
}

/// Reads records from cursors, joining by position, and hands each position
/// to the motor.
fn junta(cursores: &mut [Cursor], motor: &mut Motor) -> Result<(), String> {
    for c in cursores.iter_mut() {
        c.avanca()?;
    }
    loop {
        // the smallest position among the heads
        let mut menor: Option<(String, u64)> = None;
        for c in cursores.iter() {
            if let Some(r) = &c.cab {
                let melhor = match &menor {
                    None => true,
                    Some((mc, mp)) => {
                        let o = if *mc == r.canon { r.pos.cmp(mp) } else { natural_cmp(&r.canon, mc) };
                        o == std::cmp::Ordering::Less
                    }
                };
                if melhor {
                    menor = Some((r.canon.clone(), r.pos));
                }
            }
        }
        let Some((mc, mp)) = menor else { break };
        // every record of that position: with one file, the consecutive ones
        // (a single file need not be sorted)
        let mut regs: Vec<Registro> = Vec::new();
        for c in cursores.iter_mut() {
            while matches!(&c.cab, Some(r) if r.canon == mc && r.pos == mp) {
                regs.push(c.cab.take().unwrap());
                c.avanca()?;
            }
        }
        motor.posicao(&regs)?;
    }
    Ok(())
}

fn indexavel(caminho: &str) -> bool {
    let c = caminho.to_ascii_lowercase();
    (c.ends_with(".gz") || c.ends_with(".bgz")) && std::path::Path::new(&format!("{caminho}.tbi")).exists() && bgzf::e_bgzf(caminho)
}

/// Genotype-outcome association. `vcfs`: one multi-sample VCF or several;
/// `phenotypes`: a table (tab, comma or whitespace separated) with a header.
pub fn genotype_association(vcfs: &[String], phenotypes: &str, opt: &AssocOptions) -> Result<AssocReport, String> {
    if vcfs.is_empty() {
        return Err("No VCF given: pass one multi-sample VCF or several VCFs".into());
    }
    if opt.outcome.trim().is_empty() {
        return Err("Name the outcome column of the phenotype table (outcome=\"...\")".into());
    }
    if !(0.0..=0.5).contains(&opt.min_maf) {
        return Err(format!("min_maf must be between 0 and 0.5, not {}", opt.min_maf));
    }
    if !(0.0..=1.0).contains(&opt.min_call_rate) {
        return Err(format!("min_call_rate must be between 0 and 1, not {}", opt.min_call_rate));
    }
    if let Some(o) = &opt.output {
        let entradas: Vec<&str> = vcfs.iter().map(String::as_str).chain(std::iter::once(phenotypes)).collect();
        for e in entradas {
            if mesmo_arquivo(e, o) {
                return Err(format!("The output `{o}` is the input file `{e}` itself; writing it would destroy the input. Choose a different output path."));
            }
        }
    }
    let sel = Selecao::nova(&opt.ids, &opt.regions)?;
    let tab = le_tabela(phenotypes, opt.sample_column.as_deref())?;

    // the samples of every VCF
    let mut cabecalhos = Vec::new();
    let mut amostra_arq: HashMap<String, usize> = HashMap::new();
    for (i, v) in vcfs.iter().enumerate() {
        let (r, amostras, n) = abre_vcf(v)?;
        for a in &amostras {
            if let Some(j) = amostra_arq.insert(a.clone(), i) {
                return Err(if j == i {
                    format!("Sample \"{a}\" appears twice in {v}")
                } else {
                    format!(
                        "Sample \"{a}\" is in both {} and {v}: each person must be in one VCF (rename one, or merge the files)",
                        vcfs[j]
                    )
                });
            }
        }
        cabecalhos.push((r, amostras, n));
    }
    let an = monta_analise(&tab, &amostra_arq, opt)?;
    let pessoa: HashMap<&str, u32> = an.nomes.iter().enumerate().map(|(i, s)| (s.as_str(), i as u32)).collect();
    let arqs: Vec<Arquivo> = cabecalhos
        .iter()
        .map(|(_, amostras, _)| Arquivo {
            servidos: amostras.iter().enumerate().filter_map(|(c, a)| pessoa.get(a.as_str()).map(|&p| (c, p))).collect(),
        })
        .collect();

    let mut motor = Motor {
        an: &an,
        arqs: &arqs,
        opt,
        sel,
        lote: Vec::new(),
        linhas: Vec::new(),
        cods: Vec::new(),
        grupos: HashMap::new(),
        relogio: 0,
        pulados: HashMap::new(),
        exemplos: Vec::new(),
        simbolicos: 0,
        alelos: 0,
        haploides: 0,
    };
    let mut avisos: Vec<String> = an.avisos.clone();
    let alvos = motor.sel.alvos();
    let usa_indice = alvos.is_some() && vcfs.iter().all(|v| indexavel(v));
    let (mut lidos, mut malformados, mut pulados_seq) = (0u64, 0u64, 0u64);
    let multiplos = vcfs.len() > 1;
    if usa_indice {
        drop(cabecalhos);
        let alvos = alvos.unwrap();
        // the names and extent of each sequence in each file's index
        let mut nomes_arq: Vec<HashMap<String, (String, u64)>> = Vec::new();
        for v in vcfs {
            let (idx, _, nomes) = le_tbi(&format!("{v}.tbi"))?;
            let mut m = HashMap::new();
            for (i, nome) in nomes.iter().enumerate() {
                let lim = idx.refs.get(i).map_or(0, |r| (r.linear.len() as u64 + 1) << 14);
                m.insert(canonico(nome), (nome.clone(), lim));
            }
            nomes_arq.push(m);
        }
        // windows of at most 200 kb, each a list of intervals
        const JANELA: u64 = 200_000;
        let mut janelas: Vec<(String, Vec<(u64, u64)>)> = Vec::new();
        for (c, a, b) in alvos {
            let lim = nomes_arq.iter().filter_map(|m| m.get(&c).map(|x| x.1)).max().unwrap_or(0);
            if lim == 0 {
                continue;
            }
            let b = b.min(lim);
            let mut ini = a;
            while ini <= b {
                let fim = (ini + JANELA - 1).min(b);
                match janelas.last_mut() {
                    Some((jc, iv)) if *jc == c && fim - iv[0].0 < JANELA => iv.push((ini, fim)),
                    _ => janelas.push((c.clone(), vec![(ini, fim)])),
                }
                ini = fim + 1;
            }
        }
        for (c, iv) in &janelas {
            let mut cursores = Vec::new();
            for (fi, v) in vcfs.iter().enumerate() {
                let mut linhas: Vec<String> = Vec::new();
                if let Some((nome, _)) = nomes_arq[fi].get(c) {
                    let regioes: Vec<String> = iv.iter().map(|(a, b)| format!("{nome}:{a}-{b}")).collect();
                    consultar(v, &regioes, &mut |cab, l| {
                        if !cab {
                            let p = l.split('\t').nth(1).and_then(|x| x.trim().parse::<u64>().ok()).unwrap_or(0);
                            let k = iv.partition_point(|x| x.1 < p);
                            if k < iv.len() && iv[k].0 <= p {
                                linhas.push(l.to_string());
                            }
                        }
                        Ok(())
                    })?;
                }
                cursores.push(Cursor {
                    arq: fi,
                    caminho: v.clone(),
                    fonte: Fonte::Linhas(linhas.into_iter()),
                    n_linha: 0,
                    cab: None,
                    cache: (String::new(), String::new()),
                    multiplos: false,
                    ultimo: None,
                    nome_ultimo: String::new(),
                    registros: 0,
                    malformados: 0,
                    pulados_seq: 0,
                });
            }
            junta(&mut cursores, &mut motor)?;
            for c in &cursores {
                lidos += c.registros;
                malformados += c.malformados;
            }
        }
    } else {
        if alvos.is_some() {
            let sem: Vec<&String> = vcfs.iter().filter(|v| !indexavel(v)).collect();
            avisos.push(format!(
                "no .tbi index next to {}: the whole file was read; index_file(...) on a sorted BGZF file makes region queries fast",
                lista_curta(&sem.iter().map(|s| s.to_string()).collect::<Vec<_>>(), 3)
            ));
        }
        let mut cursores: Vec<Cursor> = cabecalhos
            .into_iter()
            .enumerate()
            .map(|(fi, (r, _, n))| Cursor {
                arq: fi,
                caminho: vcfs[fi].clone(),
                fonte: Fonte::Fluxo(r),
                n_linha: n,
                cab: None,
                cache: (String::new(), String::new()),
                multiplos,
                ultimo: None,
                nome_ultimo: String::new(),
                registros: 0,
                malformados: 0,
                pulados_seq: 0,
            })
            .collect();
        junta(&mut cursores, &mut motor)?;
        for c in &cursores {
            lidos += c.registros;
            malformados += c.malformados;
            pulados_seq += c.pulados_seq;
        }
    }
    motor.testa_lote();

    let Motor { linhas, cods, pulados, exemplos, simbolicos, alelos, haploides, sel, .. } = motor;
    let mut linhas = linhas;
    // FDR
    let ps: Vec<f64> = linhas.iter().map(|r| r.p).collect();
    let fdr = benjamini_hochberg(&ps);
    for (r, f) in linhas.iter_mut().zip(fdr) {
        r.fdr = f;
    }
    // permutations
    let mut esquema = String::new();
    if opt.permutations > 0 && !linhas.is_empty() {
        let p = permuta(&an, &linhas, &cods, opt)?;
        for (i, r) in linhas.iter_mut().enumerate() {
            r.p_permutation = p.p_emp[i];
            r.p_fwer = p.p_fwer[i];
        }
        esquema = p.esquema;
    }

    let mut rep = AssocReport {
        vcfs: vcfs.to_vec(),
        phenotypes: phenotypes.to_string(),
        outcome: opt.outcome.trim().to_string(),
        outcome_type: if an.binario { "binary".into() } else { "quantitative".into() },
        event_value: an.evento.clone(),
        other_value: an.outro.clone(),
        model: opt.model.name().into(),
        covariates: an.covariaveis.clone(),
        design: an.z_nomes.clone(),
        firth: opt.firth.name().into(),
        permutations: opt.permutations,
        permutation_scheme: esquema,
        seed: opt.seed,
        absent_is_ref: opt.absent_is_ref,
        reading: if usa_indice { "indexed".into() } else { "streamed".into() },
        vcf_samples: an.amostras_vcf,
        table_rows: an.linhas_tabela,
        analysed: an.y.len(),
        not_in_table: an.fora_da_tabela,
        not_in_vcf: an.fora_do_vcf.len(),
        not_in_vcf_examples: an.fora_do_vcf.iter().take(10).cloned().collect(),
        missing_outcome: an.falta_desfecho,
        missing_covariates: an.falta_cov.clone(),
        cases: if an.binario { an.y.iter().filter(|v| **v == 1.0).count() } else { 0 },
        controls: if an.binario { an.y.iter().filter(|v| **v == 0.0).count() } else { 0 },
        records_read: lidos,
        malformed_records: malformados,
        skipped_sequences: pulados_seq,
        symbolic_alleles: simbolicos,
        alleles_considered: alelos,
        tested: linhas.len(),
        skipped: Vec::new(),
        skipped_examples: exemplos,
        not_found: sel.pedidos.iter().zip(&sel.achados).filter(|(_, a)| !**a).map(|(p, _)| p.clone()).collect(),
        haploid_calls: haploides,
        rows: linhas,
        residual_skewness: f64::NAN,
        residual_excess_kurtosis: f64::NAN,
        power_maf: f64::NAN,
        power_alpha: f64::NAN,
        detectable_effect: f64::NAN,
        output: None,
        warnings: Vec::new(),
    };
    for m in [Motivo::Chamada, Motivo::Monomorfico, Motivo::Maf, Motivo::SemVariacao, Motivo::Colinear, Motivo::Poucos, Motivo::SemEvento, Motivo::Ajuste] {
        if let Some(&k) = pulados.get(&m) {
            rep.skipped.push((m.codigo().to_string(), k));
        }
    }

    // --- checks and honest warnings ---
    let n = an.y.len();
    let k = an.z.len();
    if !an.binario {
        let qa = qr_gs(&an.z, 1e-9);
        let e = qa.residual(&an.y);
        let m = e.iter().sum::<f64>() / n as f64;
        let m2 = e.iter().map(|x| (x - m).powi(2)).sum::<f64>() / n as f64;
        let m3 = e.iter().map(|x| (x - m).powi(3)).sum::<f64>() / n as f64;
        let m4 = e.iter().map(|x| (x - m).powi(4)).sum::<f64>() / n as f64;
        if m2 > 0.0 {
            rep.residual_skewness = m3 / m2.powf(1.5);
            rep.residual_excess_kurtosis = m4 / (m2 * m2) - 3.0;
            if rep.residual_skewness.abs() > 1.0 || rep.residual_excess_kurtosis > 3.0 {
                avisos.push(format!(
                    "the residuals of the outcome given the covariates are far from normal (skewness {:.2}, excess kurtosis {:.2}): with few people the t test can be off; transform the outcome (log, rank) or use permutations",
                    rep.residual_skewness, rep.residual_excess_kurtosis
                ));
            }
        }
    }
    if an.binario {
        let raros = rep.cases.min(rep.controls);
        if raros < 10 * (k + 1) {
            avisos.push(format!(
                "only {raros} people with the rarer outcome for {} parameters (rule of thumb: at least 10 per parameter): logistic p-values and confidence intervals are approximate; firth=\"always\" and permutations are safer",
                k + 1
            ));
        }
    }
    if n < 30 {
        avisos.push(format!("{n} people: the p-values rest on large-sample approximations and the study can detect only very large effects (see the power line)"));
    }
    if multiplos {
        if opt.absent_is_ref {
            avisos.push("a variant absent from a person's VCF was taken as homozygous reference (absent_is_ref): right only where that person's VCF could have called it".into());
        } else {
            avisos.push("a variant absent from a person's VCF counts as missing for that person (there is no information on which regions each VCF covers); single-sample VCFs that list only variant sites need absent_is_ref=True".into());
            let baixa = rep.skipped.iter().find(|(c, _)| c == "low_call_rate").map_or(0, |x| x.1);
            if baixa > 0 && baixa * 2 >= rep.alleles_considered {
                avisos.push(format!(
                    "{baixa} of {} alleles were left out for a low call rate: with single-sample VCFs of variant sites only, people without the variant have no record; absent_is_ref=True treats them as homozygous reference",
                    rep.alleles_considered
                ));
            }
        }
        if pulados_seq > 0 {
            avisos.push(format!(
                "{} record(s) on the mitochondria or on sequences that are not primary chromosomes were skipped: their order differs between references, so several VCFs cannot be joined there",
                milhar(pulados_seq as usize)
            ));
        }
    }
    if rep.haploid_calls > 0 {
        avisos.push(format!(
            "{} haploid genotype call(s) (e.g. chrX in men) counted as 0 or 1 copy; include sex as a covariate",
            milhar(rep.haploid_calls)
        ));
    }
    if rep.symbolic_alleles > 0 {
        avisos.push(format!("{} symbolic or spanning-deletion allele(s) (<DEL>, *) were not tested", rep.symbolic_alleles));
    }
    if malformados > 0 {
        avisos.push(format!("{malformados} malformed record(s) skipped"));
    }
    if !rep.not_found.is_empty() {
        avisos.push(format!("{} requested variant(s) not found in the VCF(s)", rep.not_found.len()));
    }
    let pedidos_pulados: Vec<String> = rep.skipped_examples.iter().take(5).map(|(i, m)| format!("{i}: {m}")).collect();
    if !opt.ids.is_empty() && !pedidos_pulados.is_empty() {
        avisos.push(format!("not tested: {}", pedidos_pulados.join("; ")));
    }
    if opt.permutations > 0 && opt.permutations < 1000 {
        avisos.push(format!("with {} permutations the smallest possible empirical p is {:.2e}", opt.permutations, 1.0 / (opt.permutations as f64 + 1.0)));
    }
    // the power line
    if rep.tested > 0 {
        let mut mafs: Vec<f64> = rep.rows.iter().map(|r| r.alt_freq.min(1.0 - r.alt_freq)).collect();
        mafs.sort_by(|a, b| a.total_cmp(b));
        let maf = mafs[mafs.len() / 2].max(0.01);
        let alfa = 0.05 / rep.tested as f64;
        rep.power_maf = maf;
        rep.power_alpha = alfa;
        let po = PowerOptions {
            outcome: if an.binario { OutcomeType::Binary } else { OutcomeType::Quantitative },
            allele_freq: maf,
            effect: None,
            n: Some(n),
            case_fraction: if an.binario { rep.cases as f64 / n as f64 } else { 0.5 },
            alpha: 0.05,
            tests: rep.tested,
            model: opt.model,
            covariates: k - 1,
            target_power: 0.8,
            simulations: 0,
            seed: opt.seed,
        };
        rep.detectable_effect = efeito_detectavel(&po, n).unwrap_or(f64::NAN);
    }
    rep.warnings = avisos;
    if let Some(o) = &opt.output {
        rep.grava(o)?;
    }
    Ok(rep)
}

// =====================================================================
//  12) POWER
// =====================================================================

#[derive(Clone, Debug)]
pub struct PowerOptions {
    /// Binary or quantitative (Auto is refused).
    pub outcome: OutcomeType,
    /// Frequency of the tested allele: in the people without the event
    /// (binary), or in the population (quantitative).
    pub allele_freq: f64,
    /// Odds ratio per copy / for the coded genotype (binary), or beta in SD
    /// units of the outcome (quantitative). `None`: only the people needed is
    /// not computed; the detectable effect is.
    pub effect: Option<f64>,
    /// People in the study. `None`: only the people needed.
    pub n: Option<usize>,
    /// Fraction of people with the event (binary).
    pub case_fraction: f64,
    pub alpha: f64,
    /// Number of tests (Bonferroni): the test level is alpha / tests.
    pub tests: usize,
    pub model: GeneticModel,
    /// Covariates in the model (they take degrees of freedom).
    pub covariates: usize,
    pub target_power: f64,
    /// Simulated data sets for an empirical power (0: none).
    pub simulations: usize,
    pub seed: u64,
}

impl Default for PowerOptions {
    fn default() -> Self {
        PowerOptions {
            outcome: OutcomeType::Binary,
            allele_freq: 0.2,
            effect: None,
            n: None,
            case_fraction: 0.5,
            alpha: 0.05,
            tests: 1,
            model: GeneticModel::Additive,
            covariates: 0,
            target_power: 0.8,
            simulations: 0,
            seed: 1,
        }
    }
}

/// The coded genotype's distribution under Hardy-Weinberg at allele
/// frequency p: (value, probability).
fn distribuicao(m: GeneticModel, p: f64) -> Vec<(f64, f64)> {
    let q = 1.0 - p;
    match m {
        GeneticModel::Additive => vec![(0.0, q * q), (1.0, 2.0 * p * q), (2.0, p * p)],
        GeneticModel::Dominant => vec![(0.0, q * q), (1.0, 1.0 - q * q)],
        GeneticModel::Recessive => vec![(0.0, 1.0 - p * p), (1.0, p * p)],
    }
}

/// The variance of the coded genotype in the population.
fn variancia_codigo(m: GeneticModel, p: f64) -> f64 {
    let d = distribuicao(m, p);
    let mu: f64 = d.iter().map(|(x, w)| x * w).sum();
    d.iter().map(|(x, w)| w * (x - mu).powi(2)).sum()
}

/// Non-centrality per person of the Wald test of a logistic regression of a
/// case-control sample (fraction `fi` with the event; the genotype in those
/// without the event is `dist`, in those with it the same tilted by the odds
/// ratio): beta^2 / [I^-1]_22, I the expected information per person.
fn lambda_por_pessoa(dist: &[(f64, f64)], beta: f64, fi: f64) -> f64 {
    let z: f64 = dist.iter().map(|(x, w)| w * (beta * x).exp()).sum();
    let (mut i11, mut i12, mut i22) = (0.0, 0.0, 0.0);
    for (x, w0) in dist {
        let w1 = w0 * (beta * x).exp() / z;
        let px = fi * w1 + (1.0 - fi) * w0;
        if px <= 0.0 {
            continue;
        }
        let pi = fi * w1 / px;
        let v = px * pi * (1.0 - pi);
        i11 += v;
        i12 += v * x;
        i22 += v * x * x;
    }
    let det = i11 * i22 - i12 * i12;
    if !(det > 0.0) {
        return 0.0;
    }
    beta * beta * det / i11
}

fn poder_de_lambda(lambda: f64, alfa: f64) -> f64 {
    let z = -normal_quantile(alfa / 2.0);
    let s = lambda.max(0.0).sqrt();
    (normal_upper(z - s) + normal_upper(z + s)).min(1.0)
}

/// Analytic power: binary, the Wald test of the logistic regression (large
/// sample); quantitative, the t test (non-central t, exact under normality
/// for a given genotype variance).
fn poder(o: &PowerOptions, n: usize, efeito: f64) -> f64 {
    let alfa = o.alpha / o.tests.max(1) as f64;
    match o.outcome {
        OutcomeType::Quantitative => {
            let v = variancia_codigo(o.model, o.allele_freq);
            let q2 = efeito * efeito * v;
            if q2 >= 1.0 {
                return f64::NAN;
            }
            let df = n as f64 - 2.0 - o.covariates as f64;
            if df < 1.0 {
                return f64::NAN;
            }
            let delta = (efeito * efeito * v * (n as f64 - 1.0) / (1.0 - q2)).sqrt();
            noncentral_t_power(df, delta, t_critical(df, alfa))
        }
        _ => {
            let d = distribuicao(o.model, o.allele_freq);
            let l1 = lambda_por_pessoa(&d, efeito.ln(), o.case_fraction);
            poder_de_lambda(n as f64 * l1, alfa)
        }
    }
}

fn n_minimo(o: &PowerOptions, n: usize) -> usize {
    match o.outcome {
        OutcomeType::Quantitative => (o.covariates + 3).max(n),
        _ => n.max(2),
    }
}

/// The smallest number of people with power >= target (None: above 10^9).
fn pessoas_necessarias(o: &PowerOptions, efeito: f64) -> Option<usize> {
    let lo0 = n_minimo(o, 2);
    let alvo = o.target_power;
    let ok = |n: usize| poder(o, n, efeito) >= alvo;
    if ok(lo0) {
        return Some(lo0);
    }
    let mut hi = lo0 * 2;
    while !ok(hi) {
        if hi > 1_000_000_000 {
            return None;
        }
        hi *= 2;
    }
    let mut lo = hi / 2;
    while hi - lo > 1 {
        let m = lo + (hi - lo) / 2;
        if ok(m) {
            hi = m;
        } else {
            lo = m;
        }
    }
    Some(hi)
}

/// The smallest effect (odds ratio > 1, or beta > 0 in SD units) with power
/// >= target with n people. None: none in range (OR up to 1000; beta up to
/// where the variant would explain 99% of the variance).
fn efeito_detectavel(o: &PowerOptions, n: usize) -> Option<f64> {
    let alvo = o.target_power;
    match o.outcome {
        OutcomeType::Quantitative => {
            let v = variancia_codigo(o.model, o.allele_freq);
            if !(v > 0.0) {
                return None;
            }
            let bmax = (0.99 / v).sqrt();
            if !(poder(o, n, bmax) >= alvo) {
                return None;
            }
            let (mut lo, mut hi) = (0.0, bmax);
            for _ in 0..80 {
                let m = 0.5 * (lo + hi);
                if poder(o, n, m) >= alvo {
                    hi = m;
                } else {
                    lo = m;
                }
            }
            Some(hi)
        }
        _ => {
            let lmax = 1000f64.ln();
            if !(poder(o, n, lmax.exp()) >= alvo) {
                return None;
            }
            let (mut lo, mut hi) = (0.0, lmax);
            for _ in 0..80 {
                let m = 0.5 * (lo + hi);
                if poder(o, n, m.exp()) >= alvo {
                    hi = m;
                } else {
                    lo = m;
                }
            }
            Some(hi.exp())
        }
    }
}

/// Simulated power: binary, a case-control sample tested as
/// `genotype_association` does by default (logistic likelihood-ratio test,
/// Firth's penalized one when the plain fit separates); quantitative, the
/// OLS t test with `covariates` unrelated normal covariates. Also returns the
/// rejection rate of the Wald test (binary), the analytic formula's test.
pub fn simulated_power(o: &PowerOptions, n: usize, efeito: f64, reps: usize) -> (f64, f64) {
    let alfa = o.alpha / o.tests.max(1) as f64;
    let dist = distribuicao(o.model, o.allele_freq);
    let (rej, rej_wald) = (0..reps as u64)
        .into_par_iter()
        .map(|rep| {
            let mut rng = Rng::stream(o.seed, rep);
            let sorteia = |rng: &mut Rng, d: &[(f64, f64)]| -> f64 {
                let u = rng.unif();
                let mut acc = 0.0;
                for (x, w) in d {
                    acc += w;
                    if u < acc {
                        return *x;
                    }
                }
                d.last().unwrap().0
            };
            match o.outcome {
                OutcomeType::Quantitative => {
                    let v = variancia_codigo(o.model, o.allele_freq);
                    let mu: f64 = dist.iter().map(|(x, w)| x * w).sum();
                    let sd_e = (1.0 - efeito * efeito * v).max(0.0).sqrt();
                    let x: Vec<f64> = (0..n).map(|_| sorteia(&mut rng, &dist)).collect();
                    let y: Vec<f64> = x.iter().map(|xi| efeito * (xi - mu) + sd_e * rng.normal()).collect();
                    let mut z = vec![vec![1.0; n]];
                    for _ in 0..o.covariates {
                        z.push((0..n).map(|_| rng.normal()).collect());
                    }
                    let qr = qr_gs(&z, 1e-9);
                    let xr = qr.residual(&x);
                    let yr = qr.residual(&y);
                    let xx = dot(&xr, &xr);
                    if !(xx > 0.0) {
                        return (0u32, 0u32);
                    }
                    let c = dot(&xr, &yr);
                    let df = (n - qr.k() - 1) as f64;
                    let rss = (dot(&yr, &yr) - c * c / xx).max(0.0);
                    let t = (c / xx) / (rss / df / xx).sqrt();
                    let r = (t_two_sided(t, df) < alfa) as u32;
                    (r, r)
                }
                _ => {
                    let n1 = ((n as f64) * o.case_fraction).round() as usize;
                    let n1 = n1.clamp(1, n - 1);
                    let zt: f64 = dist.iter().map(|(x, w)| w * (efeito.ln() * x).exp()).sum();
                    let d1: Vec<(f64, f64)> = dist.iter().map(|(x, w)| (*x, w * (efeito.ln() * x).exp() / zt)).collect();
                    let mut x = Vec::with_capacity(n);
                    let mut y = Vec::with_capacity(n);
                    for i in 0..n {
                        if i < n1 {
                            x.push(sorteia(&mut rng, &d1));
                            y.push(1.0);
                        } else {
                            x.push(sorteia(&mut rng, &dist));
                            y.push(0.0);
                        }
                    }
                    let x0 = x[0];
                    if x.iter().all(|v| *v == x0) {
                        return (0, 0);
                    }
                    let cols = vec![vec![1.0; n], x];
                    let f = logistic_fit(&cols, &y, None);
                    let m = y.iter().sum::<f64>() / n as f64;
                    let l0 = y.iter().map(|v| if *v == 1.0 { m.ln() } else { (1.0 - m).ln() }).sum::<f64>();
                    let wald = if !f.separated && f.converged { (p_two_sided_z(f.beta[1] / f.se[1]) < alfa) as u32 } else { 0 };
                    let p = if !f.separated && f.converged {
                        chi2_1_upper(2.0 * (f.loglik - l0))
                    } else {
                        match firth_fit(&cols, &y, None, None) {
                            Some(ff) => firth_lrt(&cols, &y, &ff, 1).1,
                            None => 1.0,
                        }
                    };
                    ((p < alfa) as u32, wald)
                }
            }
        })
        .reduce(|| (0, 0), |a, b| (a.0 + b.0, a.1 + b.1));
    (rej as f64 / reps as f64, rej_wald as f64 / reps as f64)
}

#[derive(Clone, Debug, Default)]
pub struct PowerReport {
    pub outcome_type: String,
    pub model: String,
    pub allele_freq: f64,
    pub effect: f64,
    pub n: Option<usize>,
    pub case_fraction: f64,
    pub alpha: f64,
    pub tests: usize,
    pub alpha_per_test: f64,
    pub covariates: usize,
    pub target_power: f64,
    pub power: f64,
    pub n_needed: Option<usize>,
    pub detectable_effect: f64,
    pub simulations: usize,
    pub simulated_power: f64,
    pub simulated_se: f64,
    pub simulated_power_wald: f64,
    pub expected_events: f64,
    pub expected_carriers: f64,
    pub curve_n: Vec<usize>,
    pub curve_power: Vec<f64>,
    pub method: String,
    pub warnings: Vec<String>,
}

impl PowerReport {
    pub fn json(&self) -> String {
        format!(
            "{{\"outcome_type\":{},\"model\":{},\"allele_freq\":{},\"effect\":{},\"effect_scale\":{},\"n\":{},\"case_fraction\":{},\"alpha\":{},\"tests\":{},\"alpha_per_test\":{},\"covariates\":{},\"target_power\":{},\"power\":{},\"n_needed\":{},\"detectable_effect\":{},\"simulations\":{},\"simulated_power\":{},\"simulated_se\":{},\"simulated_power_wald\":{},\"expected_events\":{},\"expected_carriers\":{},\"curve\":{{\"n\":[{}],\"power\":[{}]}},\"method\":{},\"warnings\":{},\"text\":{}}}",
            json_texto(&self.outcome_type),
            json_texto(&self.model),
            jnum(self.allele_freq),
            jnum(self.effect),
            json_texto(if self.outcome_type == "binary" { "odds_ratio" } else { "beta_sd" }),
            self.n.map_or("null".into(), |x| x.to_string()),
            if self.outcome_type == "binary" { jnum(self.case_fraction) } else { "null".into() },
            jnum(self.alpha),
            self.tests,
            jnum(self.alpha_per_test),
            self.covariates,
            jnum(self.target_power),
            jnum(self.power),
            self.n_needed.map_or("null".into(), |x| x.to_string()),
            jnum(self.detectable_effect),
            self.simulations,
            jnum(self.simulated_power),
            jnum(self.simulated_se),
            jnum(self.simulated_power_wald),
            jnum(self.expected_events),
            jnum(self.expected_carriers),
            self.curve_n.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","),
            self.curve_power.iter().map(|x| jnum(*x)).collect::<Vec<_>>().join(","),
            json_texto(&self.method),
            jtextos(&self.warnings),
            json_texto(&self.text())
        )
    }

    pub fn text(&self) -> String {
        let bin = self.outcome_type == "binary";
        let mut o = format!(
            "Power of a {} association test ({} model), alpha {} / {} test(s) = {} per test.\n",
            self.outcome_type,
            self.model,
            self.alpha,
            self.tests,
            fmt_p(self.alpha_per_test)
        );
        if bin {
            o.push_str(&format!(
                "  allele frequency {} (in people without the event); {}% of the people have the event",
                self.allele_freq,
                fmt_g(100.0 * self.case_fraction, 1)
            ));
        } else {
            o.push_str(&format!("  allele frequency {}", self.allele_freq));
        }
        if self.effect.is_finite() {
            o.push_str(&format!("; {} {}", if bin { "odds ratio" } else { "effect (SD per copy)" }, self.effect));
        }
        o.push_str(".\n");
        if let Some(n) = self.n {
            if self.power.is_finite() {
                o.push_str(&format!("  with {n} people: power {:.3} ({})", self.power, self.method));
                if self.simulations > 0 {
                    o.push_str(&format!(
                        "; simulated over {} data sets: {:.3} +- {:.3}",
                        self.simulations, self.simulated_power, self.simulated_se
                    ));
                }
                o.push_str(".\n");
            }
            if self.detectable_effect.is_finite() {
                o.push_str(&format!(
                    "  smallest {} detectable with {:.0}% power with {n} people: {}.\n",
                    if bin { "odds ratio" } else { "effect (SD per copy)" },
                    100.0 * self.target_power,
                    fmt_g(self.detectable_effect, 3)
                ));
            } else {
                o.push_str(&format!(
                    "  with {n} people no {} reaches {:.0}% power.\n",
                    if bin { "odds ratio up to 1000" } else { "effect" },
                    100.0 * self.target_power
                ));
            }
            let enorme = if bin { !(self.detectable_effect <= 10.0) } else { !(self.detectable_effect <= 1.0) };
            if enorme {
                o.push_str(&format!(
                    "  That is far beyond the effects common variants have (most replicated associations are odds ratios below 1.5, or below 0.1 SD per allele): a study of {n} people cannot be expected to detect a real association, and a \"significant\" result would most likely be a false positive or a large overestimate.\n"
                ));
            }
        }
        if self.effect.is_finite() {
            match self.n_needed {
                Some(k) => o.push_str(&format!("  people needed for {:.0}% power: {}.\n", 100.0 * self.target_power, milhar(k))),
                None => o.push_str("  no study size up to 10^9 reaches that power for this effect.\n"),
            }
        }
        for w in &self.warnings {
            o.push_str(&format!("  note: {w}\n"));
        }
        o
    }
}

/// Power of an association test, the people needed and the smallest
/// detectable effect.
pub fn association_power(o: &PowerOptions) -> Result<PowerReport, String> {
    let bin = match o.outcome {
        OutcomeType::Binary => true,
        OutcomeType::Quantitative => false,
        OutcomeType::Auto => return Err("Say whether the outcome is \"binary\" or \"quantitative\"".into()),
    };
    if !(o.allele_freq > 0.0 && o.allele_freq < 1.0) {
        return Err(format!("allele_freq must be between 0 and 1 (exclusive), not {}", o.allele_freq));
    }
    if !(o.alpha > 0.0 && o.alpha < 1.0) {
        return Err(format!("alpha must be between 0 and 1, not {}", o.alpha));
    }
    if o.tests == 0 {
        return Err("tests must be at least 1".into());
    }
    if !(o.target_power > 0.0 && o.target_power < 1.0) {
        return Err(format!("target_power must be between 0 and 1, not {}", o.target_power));
    }
    if bin && !(o.case_fraction > 0.0 && o.case_fraction < 1.0) {
        return Err(format!("case_fraction must be between 0 and 1 (exclusive), not {}", o.case_fraction));
    }
    if o.effect.is_none() && o.n.is_none() {
        return Err("Give the effect (odds ratio or beta), the number of people, or both".into());
    }
    if let Some(e) = o.effect {
        if bin && !(e > 0.0 && e.is_finite()) {
            return Err(format!("The odds ratio must be positive, not {e}"));
        }
        if !bin {
            let v = variancia_codigo(o.model, o.allele_freq);
            if !(e.is_finite()) || e * e * v >= 1.0 {
                return Err(format!(
                    "An effect of {e} SD per copy at frequency {} would explain {:.0}% or more of the outcome's variance: beta is in units of the outcome's standard deviation",
                    o.allele_freq,
                    100.0 * e * e * v
                ));
            }
        }
    }
    if let Some(n) = o.n {
        let minimo = if bin { 4 } else { o.covariates + 3 };
        if n < minimo {
            return Err(format!("{n} people: a test with {} covariate(s) needs at least {minimo}", o.covariates));
        }
    }
    // an OR below 1 is the same problem as 1/OR for the other allele only
    // approximately; it is computed as given
    let mut r = PowerReport {
        outcome_type: if bin { "binary".into() } else { "quantitative".into() },
        model: o.model.name().into(),
        allele_freq: o.allele_freq,
        effect: o.effect.unwrap_or(f64::NAN),
        n: o.n,
        case_fraction: o.case_fraction,
        alpha: o.alpha,
        tests: o.tests,
        alpha_per_test: o.alpha / o.tests as f64,
        covariates: o.covariates,
        target_power: o.target_power,
        power: f64::NAN,
        n_needed: None,
        detectable_effect: f64::NAN,
        simulations: o.simulations,
        simulated_power: f64::NAN,
        simulated_se: f64::NAN,
        simulated_power_wald: f64::NAN,
        expected_events: f64::NAN,
        expected_carriers: f64::NAN,
        method: if bin {
            "Wald test of the logistic regression, large-sample".into()
        } else {
            "t test, non-central t".into()
        },
        ..Default::default()
    };
    if let Some(e) = o.effect {
        if let Some(n) = o.n {
            r.power = poder(o, n, e);
            if o.simulations > 0 {
                let (s, sw) = simulated_power(o, n, e, o.simulations);
                r.simulated_power = s;
                r.simulated_power_wald = if bin { sw } else { f64::NAN };
                r.simulated_se = (s * (1.0 - s) / o.simulations as f64).sqrt();
            }
        }
        r.n_needed = pessoas_necessarias(o, e);
    }
    if let Some(n) = o.n {
        r.detectable_effect = efeito_detectavel(o, n).unwrap_or(f64::NAN);
        let d = distribuicao(o.model, o.allele_freq);
        let portador: f64 = d.iter().filter(|(x, _)| *x > 0.0).map(|(_, w)| w).sum();
        if bin {
            r.expected_events = n as f64 * o.case_fraction;
            r.expected_carriers = n as f64 * portador;
            let eventos = (n as f64 * o.case_fraction.min(1.0 - o.case_fraction)).round();
            if eventos < 10.0 || r.expected_carriers < 10.0 {
                r.warnings.push(format!(
                    "about {} people with the rarer outcome and {:.1} carriers expected: the large-sample power is unreliable here; simulations=1000 gives a better estimate",
                    eventos, r.expected_carriers
                ));
            }
        } else {
            r.expected_carriers = n as f64 * portador;
            if r.expected_carriers < 5.0 {
                r.warnings.push(format!(
                    "{:.1} carriers expected: the t test then rests on a handful of people",
                    r.expected_carriers
                ));
            }
            if n < 100 {
                r.warnings.push(
                    "below ~100 people the analytic power, which uses the expected genotype variance, can be a few points too high (0.03 at 40 people in the tests); simulations=2000 checks it".into(),
                );
            }
        }
    }
    if !bin {
        r.warnings.push("beta is the change of the outcome, in standard deviations, per copy of the allele (or for the coded genotype)".into());
    } else {
        r.warnings.push("covariates are assumed unrelated to the genotype; a covariate strongly related to both lowers the power".into());
    }
    // the curve
    let alvo_n = r.n_needed.unwrap_or(0).max(o.n.unwrap_or(0)).max(20);
    let lo = n_minimo(o, 4) as f64;
    let hi = ((2 * alvo_n) as f64).max(lo * 4.0).min(1e8);
    if let Some(e) = o.effect {
        let mut vistos = HashSet::new();
        for i in 0..30 {
            let n = (lo * (hi / lo).powf(i as f64 / 29.0)).round() as usize;
            if vistos.insert(n) {
                r.curve_n.push(n);
                r.curve_power.push(poder(o, n, e));
            }
        }
    }
    Ok(r)
}
