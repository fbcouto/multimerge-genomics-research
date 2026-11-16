// =====================================================================
//  ANCESTRY PCA: PEOPLE PROJECTED ONTO THE GROUPS OF A FREQUENCY PANEL
//
//  The plot researchers draw with the 1000 Genomes: the reference
//  individuals as coloured clouds on the first principal components, and
//  the people analysed projected onto the same axes. A panel written by
//  `build_panel` / `panel_from_frequencies` holds allele frequencies per
//  group, not individuals, so the reference individuals are SIMULATED from
//  the frequencies. In order:
//
//   1. SITES. The panel's biallelic SNVs whose mean frequency over the
//      groups has a minor-allele frequency of at least `min_maf`, observed
//      in the people (genotyped, or homozygous reference inside a callable
//      BED -- the same reading as `estimate_ancestry`). Sites observed in
//      every person are preferred; when they are too few, sites observed in
//      fewer people are let in, most observed first. Then at most
//      `max_sites`, spread evenly along the genome: the smallest spacing (at
//      least `min_distance_bp`) that keeps no more than `max_sites`, taken
//      left to right on each chromosome.
//   2. REFERENCE. `per_group` individuals per group, each genotype drawn
//      Binomial(2, p_group) independently at each site. Independent sites:
//      no linkage, and no structure inside a group (the 1000 Genomes EUR
//      mixes Finns and Italians; a simulated EUR does not). The clouds show
//      the drift between the groups, not real people.
//   3. AXES. Each site is centred on its mean in the reference and divided
//      by sqrt(2 p (1 - p)), p = (1 + ALT alleles) / (2 + 2 n), as
//      EIGENSOFT does. The PCA is the eigendecomposition of the n x n Gram
//      matrix X X' / m of the reference alone (Householder reduction to a
//      tridiagonal matrix, then the implicit QL algorithm: no library).
//   4. PROJECTION. A person's standardized genotypes are projected onto the
//      axes by least squares on the sites the person has (smartpca's
//      `lsqproject`). With every site observed it is the ordinary
//      projection; with sites missing it equals mean imputation followed by
//      a rescaling by the share of each axis the observed sites carry, so a
//      person with half the sites missing is not pulled toward 0 -- the
//      position is only less precise.
//   5. SHRINKAGE. With many more sites than reference individuals, people
//      PROJECTED onto PCs lie closer to 0 than the individuals the PCs were
//      computed from, because each of those carries its own noise into the
//      axes (Lee, Zou & Wright 2010; the reason for LASER/TRACE and OADP).
//      Since the reference is simulated, the cure is free: the axes come
//      from one simulated sample, and the clouds returned are a SECOND,
//      independent sample, projected exactly as the people are. People and
//      clouds are on the same footing. `shrinkage` says, per component, how
//      much projection pulls the group centres in (1 = not at all), for
//      information only: nothing needs correcting.
//
//  The centroid of a group is the projection of its expected genotypes
//  (2 x frequency): the centre of its cloud without sampling noise. An
//  admixed person falls, on average, at the mixture of the centroids
//  (projection is linear). With K groups the reference has K - 1 dimensions
//  of real structure: components beyond K - 1 describe only the simulated
//  sampling noise.
//
//  Scale: a score is (projection onto a unit axis) / sqrt(sites), so it does
//  not depend on how many sites or reference individuals were used. The
//  squared distance between two centroids, summed over the components that
//  hold the structure, approaches the mean over the sites of
//  (2 p_a - 2 p_b)^2 / (2 p (1 - p)) -- about 4 F_ST between the two groups
//  -- times the share of the true between-group directions the axes of a
//  finite simulated sample catch: ~0.94 with 100 per group and F_ST 0.08 on
//  3,000 sites, ~0.985 with 400 per group (the spiked covariance model
//  predicts it; trials/src/bin/pca_ancestralidade.rs checks it).
//
//  Memory: the panel sites that pass the MAF filter (about 30 bytes each),
//  one byte per such site per person, and the reference (per_group x groups
//  x sites, a byte and a 4-byte float each).
// =====================================================================

use crate::ancestralidade::Person;
use crate::externo::Escritor;
use crate::multimerge_core::multi_merge_sort_stable_by;
use crate::nomes::json_texto;
use rayon::prelude::*;
use std::collections::HashMap;

// The panel, the BEDs and the people's genotypes are read by the helpers of
// ancestralidade.rs, exactly as `estimate_ancestry` reads them.
use crate::ancestralidade::{Canon, FALTA, Linhas, Motivo, Obs, Regioes, abre_painel, classifica_snv, compara_crom, divide, duas_casas, e_acgt, junta, le_regioes_das_pessoas, linha_painel, por_vcf, splitmix, txt};

// =====================================================================
//  SYMMETRIC EIGENDECOMPOSITION
// =====================================================================

/// Eigenvalues (ascending) and eigenvectors of a symmetric n x n matrix
/// given row by row (only its lower triangle is read, so a matrix that is
/// symmetric up to rounding is fine). Eigenvector `c` is
/// `vectors[c * n .. c * n + n]`, of unit length.
///
/// Householder reduction to a tridiagonal matrix (`tred2`) and the implicit
/// QL algorithm (`tql2`) of Bowdler, Martin, Reinsch and Wilkinson
/// (Handbook for Automatic Computation, vol. II; EISPACK), as transcribed in
/// the public-domain JAMA library. The matrix is kept column by column, so
/// every O(n^3) loop walks contiguous memory.
pub fn symmetric_eigen(a: &[f64], n: usize) -> Result<(Vec<f64>, Vec<f64>), String> {
    if a.len() != n * n {
        return Err(format!("symmetric_eigen: {} values for a {n} x {n} matrix.", a.len()));
    }
    if n == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    if a.iter().any(|x| !x.is_finite()) {
        return Err("symmetric_eigen: the matrix has a value that is not finite.".to_string());
    }
    // Column-major copy of the symmetrised lower triangle: v[j * n + i] = V[i][j].
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let x = a[i * n + j];
            v[j * n + i] = x;
            v[i * n + j] = x;
        }
    }
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];
    let norma = a.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    tred2(&mut v, &mut d, &mut e, n, norma);
    tql2(&mut v, &mut d, &mut e, n)?;
    Ok((d, v))
}

/// Householder tridiagonalisation; `v` (column-major) ends with the
/// orthogonal transformation, `d` the diagonal, `e` the subdiagonal.
/// `norma`: the largest |entry| of the matrix. A row whose part left of the
/// diagonal is below eps x norma is taken as already reduced (JAMA tests
/// for exactly 0): otherwise a rank-one matrix sends the reduction down to
/// denormal numbers, ten times slower. The change is backward stable.
fn tred2(v: &mut [f64], d: &mut [f64], e: &mut [f64], n: usize, norma: f64) {
    let at = |i: usize, j: usize| j * n + i;
    for j in 0..n {
        d[j] = v[at(n - 1, j)];
    }
    for i in (1..n).rev() {
        let mut scale = 0.0;
        let mut h = 0.0;
        for x in &d[..i] {
            scale += x.abs();
        }
        if scale <= f64::EPSILON * norma {
            e[i] = d[i - 1];
            for j in 0..i {
                d[j] = v[at(i - 1, j)];
                v[at(i, j)] = 0.0;
                v[at(j, i)] = 0.0;
            }
        } else {
            for x in &mut d[..i] {
                *x /= scale;
                h += *x * *x;
            }
            let mut f = d[i - 1];
            let mut g = h.sqrt();
            if f > 0.0 {
                g = -g;
            }
            e[i] = scale * g;
            h -= f * g;
            d[i - 1] = f - g;
            for x in &mut e[..i] {
                *x = 0.0;
            }
            for j in 0..i {
                f = d[j];
                v[at(j, i)] = f;
                g = e[j] + v[at(j, j)] * f;
                for k in j + 1..i {
                    g += v[at(k, j)] * d[k];
                    e[k] += v[at(k, j)] * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for j in 0..i {
                e[j] /= h;
                f += e[j] * d[j];
            }
            let hh = f / (h + h);
            for j in 0..i {
                e[j] -= hh * d[j];
            }
            for j in 0..i {
                f = d[j];
                g = e[j];
                for k in j..i {
                    v[at(k, j)] -= f * e[k] + g * d[k];
                }
                d[j] = v[at(i - 1, j)];
                v[at(i, j)] = 0.0;
            }
        }
        d[i] = h;
    }
    // Accumulate the transformations.
    for i in 0..n - 1 {
        v[at(n - 1, i)] = v[at(i, i)];
        v[at(i, i)] = 1.0;
        let h = d[i + 1];
        if h != 0.0 {
            for k in 0..=i {
                d[k] = v[at(k, i + 1)] / h;
            }
            for j in 0..=i {
                let mut g = 0.0;
                for k in 0..=i {
                    g += v[at(k, i + 1)] * v[at(k, j)];
                }
                for k in 0..=i {
                    v[at(k, j)] -= g * d[k];
                }
            }
        }
        for k in 0..=i {
            v[at(k, i + 1)] = 0.0;
        }
    }
    for j in 0..n {
        d[j] = v[at(n - 1, j)];
        v[at(n - 1, j)] = 0.0;
    }
    v[at(n - 1, n - 1)] = 1.0;
    e[0] = 0.0;
}

/// Implicit QL on the tridiagonal matrix (d, e), accumulating the rotations
/// in `v`; ends with the eigenvalues in `d`, sorted ascending.
fn tql2(v: &mut [f64], d: &mut [f64], e: &mut [f64], n: usize) -> Result<(), String> {
    let at = |i: usize, j: usize| j * n + i;
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;
    let mut f = 0.0;
    // Deflation threshold: relative to the whole tridiagonal matrix, not to
    // the running maximum of JAMA/EISPACK -- with that, a matrix whose first
    // rows are negligible (a rank-one matrix) is iterated on numbers near
    // underflow, 50 times slower. Dropping an element below eps x |T| is
    // backward stable.
    let tst1: f64 = (0..n).map(|i| d[i].abs() + e[i].abs()).fold(0.0, f64::max);
    let eps = f64::EPSILON;
    for l in 0..n {
        let mut m = l;
        while m < n - 1 && e[m].abs() > eps * tst1 {
            m += 1;
        }
        if m > l {
            let mut iteracoes = 0;
            loop {
                iteracoes += 1;
                if iteracoes > 100 {
                    return Err(format!("symmetric_eigen: the QL iteration did not converge for eigenvalue {l} of {n}."));
                }
                let mut g = d[l];
                let mut p = (d[l + 1] - g) / (2.0 * e[l]);
                let mut r = p.hypot(1.0);
                if p < 0.0 {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                let dl1 = d[l + 1];
                let mut h = g - d[l];
                for x in &mut d[l + 2..n] {
                    *x -= h;
                }
                f += h;
                p = d[m];
                let mut c = 1.0;
                let mut c2 = c;
                let mut c3 = c;
                let el1 = e[l + 1];
                let mut s = 0.0;
                let mut s2 = 0.0;
                for i in (l..m).rev() {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = p.hypot(e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    let (esq, dir) = v.split_at_mut(at(0, i + 1));
                    let col_i = &mut esq[at(0, i)..at(0, i) + n];
                    let col_i1 = &mut dir[..n];
                    for (vi, vi1) in col_i.iter_mut().zip(col_i1.iter_mut()) {
                        let t = *vi1;
                        *vi1 = s * *vi + c * t;
                        *vi = c * *vi - s * t;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
                if e[l].abs() <= eps * tst1 {
                    break;
                }
            }
        }
        d[l] += f;
        e[l] = 0.0;
    }
    // Sort ascending, with the vectors.
    for i in 0..n.saturating_sub(1) {
        let mut k = i;
        let mut p = d[i];
        for (j, &x) in d.iter().enumerate().skip(i + 1) {
            if x < p {
                k = j;
                p = x;
            }
        }
        if k != i {
            d[k] = d[i];
            d[i] = p;
            for r in 0..n {
                v.swap(at(r, i), at(r, k));
            }
        }
    }
    Ok(())
}

// =====================================================================
//  PCA OF A REFERENCE, AND PROJECTION
// =====================================================================

/// A PCA fitted on reference genotypes. Scores are (projection onto a unit
/// axis) / sqrt(sites); see the module notes.
#[derive(Clone, Debug)]
pub struct PcaFit {
    pub individuals: usize,
    pub sites: usize,
    pub components: usize,
    /// Per site: the mean genotype of the reference (0..2 ALT alleles).
    pub mean: Vec<f64>,
    /// Per site: sqrt(2 p (1 - p)), p = (1 + ALT alleles) / (2 + 2 n).
    pub sd: Vec<f64>,
    /// The `components` largest eigenvalues of X X' / sites, descending.
    pub eigenvalues: Vec<f64>,
    /// Sum of all the eigenvalues (the trace): the total variance.
    pub total_variance: f64,
    /// Unit-length axes, site by site: `loadings[j * components + c]`.
    pub loadings: Vec<f64>,
    /// Scores of the reference individuals: `scores[i * components + c]`.
    pub scores: Vec<f64>,
}

/// Fits a PCA on `individuals` x `sites` genotypes (row by row; 0, 1 or 2
/// ALT alleles, anything else missing). Sites are standardized as EIGENSOFT
/// does (missing genotypes count as the mean); the axes are the eigenvectors
/// of the Gram matrix X X' / sites.
pub fn pca_fit(genotypes: &[u8], individuals: usize, sites: usize, components: usize) -> Result<PcaFit, String> {
    let (n, m, k) = (individuals, sites, components);
    if genotypes.len() != n * m {
        return Err(format!("pca_fit: {} genotypes for {n} individuals x {m} sites.", genotypes.len()));
    }
    if k == 0 {
        return Err("components must be at least 1.".to_string());
    }
    if n < k + 1 || m < k + 1 {
        return Err(format!("{k} component(s) need at least {} individuals and {} sites ({n} and {m} given).", k + 1, k + 1));
    }
    let mut mean = vec![0.0; m];
    let mut sd = vec![1.0; m];
    let mut soma = vec![0u64; m];
    let mut vistos = vec![0u64; m];
    for linha in genotypes.chunks(m) {
        for (j, &g) in linha.iter().enumerate() {
            if g <= 2 {
                soma[j] += u64::from(g);
                vistos[j] += 1;
            }
        }
    }
    for j in 0..m {
        if vistos[j] > 0 {
            mean[j] = soma[j] as f64 / vistos[j] as f64;
            let p = (1.0 + soma[j] as f64) / (2.0 + 2.0 * vistos[j] as f64);
            sd[j] = (2.0 * p * (1.0 - p)).sqrt();
        }
    }
    let mut x = vec![0f32; n * m];
    x.par_chunks_mut(m).zip(genotypes.par_chunks(m)).for_each(|(xl, gl)| {
        for j in 0..m {
            xl[j] = if gl[j] <= 2 { ((f64::from(gl[j]) - mean[j]) / sd[j]) as f32 } else { 0.0 };
        }
    });
    let gram = gram(&x, n, m);
    let (d, vetores) = symmetric_eigen(&gram, n)?;
    let total: f64 = d.iter().sum();
    let mut eigenvalues = Vec::with_capacity(k);
    let mut scores = vec![0.0; n * k];
    let mut loadings = vec![0.0; m * k];
    for c in 0..k {
        let idx = n - 1 - c;
        let mu = d[idx];
        if mu.is_nan() || mu <= 1e-12 * total.max(1e-300) {
            return Err(format!(
                "Component {} has no variance: the reference has only {c} independent direction(s) at these sites. Ask for fewer components or use more sites.",
                c + 1
            ));
        }
        eigenvalues.push(mu);
        let u = &vetores[idx * n..idx * n + n];
        for i in 0..n {
            scores[i * k + c] = mu.sqrt() * u[i];
        }
    }
    // Loadings: v_c = X' u_c / sqrt(m mu_c), site by site, in parallel.
    let fatores: Vec<f64> = eigenvalues.iter().map(|mu| 1.0 / (m as f64 * mu).sqrt()).collect();
    loadings.par_chunks_mut(k).enumerate().for_each(|(j, lj)| {
        for c in 0..k {
            let idx = n - 1 - c;
            let u = &vetores[idx * n..idx * n + n];
            let mut s = 0.0;
            for i in 0..n {
                s += f64::from(x[i * m + j]) * u[i];
            }
            lj[c] = s * fatores[c];
        }
    });
    Ok(PcaFit { individuals: n, sites: m, components: k, mean, sd, eigenvalues, total_variance: total, loadings, scores })
}

/// X X' / m of a row-major n x m matrix, computed row by row in parallel.
fn gram(x: &[f32], n: usize, m: usize) -> Vec<f64> {
    let linhas: Vec<Vec<f64>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let xi = &x[i * m..i * m + m];
            (0..=i)
                .map(|j| {
                    let xj = &x[j * m..j * m + m];
                    let mut acc = [0.0f64; 4];
                    let mut t = 0;
                    while t + 4 <= m {
                        for q in 0..4 {
                            acc[q] += f64::from(xi[t + q]) * f64::from(xj[t + q]);
                        }
                        t += 4;
                    }
                    let mut s = acc[0] + acc[1] + acc[2] + acc[3];
                    while t < m {
                        s += f64::from(xi[t]) * f64::from(xj[t]);
                        t += 1;
                    }
                    s / m as f64
                })
                .collect()
        })
        .collect();
    let mut g = vec![0.0; n * n];
    for (i, l) in linhas.iter().enumerate() {
        for (j, &v) in l.iter().enumerate() {
            g[i * n + j] = v;
            g[j * n + i] = v;
        }
    }
    g
}

impl PcaFit {
    /// Projects standardized values (NaN = missing) by least squares on the
    /// observed sites: the scores b minimising sum over observed j of
    /// (y_j - sum_c b_c v_jc)^2, divided by sqrt(sites). With every site
    /// observed this is the plain projection. `None` when fewer than
    /// `components + 1` sites are observed or the system is singular.
    /// Returns the scores and the number of sites observed.
    pub fn project_standardized(&self, y: &[f64]) -> Option<(Vec<f64>, usize)> {
        let k = self.components;
        if y.len() != self.sites {
            return None;
        }
        let mut a = vec![0.0; k * k];
        let mut b = vec![0.0; k];
        let mut vistos = 0usize;
        for (j, &yj) in y.iter().enumerate() {
            if yj.is_nan() {
                continue;
            }
            vistos += 1;
            let v = &self.loadings[j * k..j * k + k];
            for c in 0..k {
                b[c] += v[c] * yj;
                for e in 0..=c {
                    a[c * k + e] += v[c] * v[e];
                }
            }
        }
        if vistos < k + 1 {
            return None;
        }
        for c in 0..k {
            for e in 0..c {
                a[e * k + c] = a[c * k + e];
            }
        }
        let sol = resolve_spd(a, b, k)?;
        let escala = 1.0 / (self.sites as f64).sqrt();
        Some((sol.into_iter().map(|x| x * escala).collect(), vistos))
    }

    /// Projects genotypes (0, 1, 2 ALT alleles; anything else is missing).
    pub fn project(&self, genotypes: &[u8]) -> Option<(Vec<f64>, usize)> {
        if genotypes.len() != self.sites {
            return None;
        }
        let y: Vec<f64> = genotypes
            .iter()
            .enumerate()
            .map(|(j, &g)| if g <= 2 { (f64::from(g) - self.mean[j]) / self.sd[j] } else { f64::NAN })
            .collect();
        self.project_standardized(&y)
    }

    /// Turns axis `c` around (its sign is arbitrary).
    fn inverte(&mut self, c: usize) {
        let k = self.components;
        for j in 0..self.sites {
            self.loadings[j * k + c] = -self.loadings[j * k + c];
        }
        for i in 0..self.individuals {
            self.scores[i * k + c] = -self.scores[i * k + c];
        }
    }
}

/// Solves A x = b for a small symmetric positive definite A (Cholesky).
/// `None` when A is (numerically) singular.
fn resolve_spd(mut a: Vec<f64>, mut b: Vec<f64>, k: usize) -> Option<Vec<f64>> {
    let escala = (0..k).map(|c| a[c * k + c]).fold(0.0f64, f64::max);
    if escala <= 0.0 || !escala.is_finite() {
        return None;
    }
    for c in 0..k {
        let mut s = a[c * k + c];
        for e in 0..c {
            s -= a[c * k + e] * a[c * k + e];
        }
        if s <= 1e-10 * escala {
            return None;
        }
        let l = s.sqrt();
        a[c * k + c] = l;
        for r in c + 1..k {
            let mut t = a[r * k + c];
            for e in 0..c {
                t -= a[r * k + e] * a[c * k + e];
            }
            a[r * k + c] = t / l;
        }
    }
    for c in 0..k {
        let mut t = b[c];
        for e in 0..c {
            t -= a[c * k + e] * b[e];
        }
        b[c] = t / a[c * k + c];
    }
    for c in (0..k).rev() {
        let mut t = b[c];
        for e in c + 1..k {
            t -= a[e * k + c] * b[e];
        }
        b[c] = t / a[c * k + c];
    }
    Some(b)
}

// =====================================================================
//  OPTIONS AND THE PANEL
// =====================================================================

/// Options of `ancestry_pca`.
#[derive(Clone, Debug, PartialEq)]
pub struct PcaOptions {
    /// One VCF: the samples to analyse (None: every sample of the file).
    /// Several VCFs: one name per VCF ("" = the file's only sample).
    pub samples: Option<Vec<String>>,
    /// Callable-regions BEDs, as in `estimate_ancestry`: None, one for
    /// everybody, or one per person ("" = none for that person).
    pub callable_beds: Option<Vec<String>>,
    /// Without a callable BED: a panel site absent from the VCF is
    /// homozygous reference (on chromosomes the VCF has records for).
    pub absent_is_ref: bool,
    /// A record whose FILTER is not PASS or `.` makes the site missing.
    pub pass_only: bool,
    /// At most this many sites, spread evenly along the genome.
    pub max_sites: usize,
    /// Minor-allele frequency of the mean of the groups' frequencies.
    pub min_maf: f64,
    /// Minimum distance between two chosen sites on a chromosome.
    pub min_distance_bp: u64,
    /// Simulated reference individuals per group.
    pub per_group: usize,
    pub components: usize,
    pub seed: u64,
    /// Optional TSV of the sites used: frequencies, mean, sd and loadings.
    pub sites_output: Option<String>,
}

impl Default for PcaOptions {
    fn default() -> Self {
        PcaOptions {
            samples: None,
            callable_beds: None,
            absent_is_ref: true,
            pass_only: true,
            max_sites: 5000,
            min_maf: 0.05,
            min_distance_bp: 1000,
            per_group: 100,
            components: 4,
            seed: 1,
            sites_output: None,
        }
    }
}

/// Simulated reference individuals at most (the eigendecomposition is
/// O(n^3): 500 take a fraction of a second, 2000 some tens of seconds).
const MAX_REFERENCIA: usize = 2000;
/// Reference individuals x sites at most (a byte and a 4-byte float each).
const MAX_CELULAS: usize = 50_000_000;
/// People x candidate sites at most (one byte each).
const MAX_CODIGOS: usize = 1_000_000_000;

/// The SNV sites of a panel that pass the MAF filter, with the frequencies
/// as written.
struct PainelPca {
    grupos: Vec<String>,
    /// Canonical chromosome names, in order (1..22, X, Y, MT, others).
    cromos: Vec<String>,
    /// The name the panel writes for each.
    nomes_crom: Vec<String>,
    inicio: Vec<usize>,
    pos: Vec<u64>,
    ra: Vec<[u8; 2]>,
    /// K frequencies per site.
    p: Vec<f32>,
    linhas: usize,
    snvs: usize,
    nao_snv: usize,
    sem_dados: usize,
    abaixo_maf: usize,
    duplicados: usize,
    arredondado: bool,
}

impl PainelPca {
    fn le(caminho: &str, min_maf: f64) -> Result<PainelPca, String> {
        let (mut l, c) = abre_painel(caminho)?;
        let k = c.grupos.len();
        struct Tmp {
            pos: Vec<u64>,
            ra: Vec<[u8; 2]>,
            p: Vec<f32>,
            nome: String,
        }
        let mut indice: HashMap<String, usize> = HashMap::new();
        let mut tmp: Vec<(String, Tmp)> = Vec::new();
        let mut canon = Canon::default();
        let mut atual = usize::MAX;
        let mut chave_atual = String::new();
        let mut sp = Vec::new();
        let (mut linhas, mut nao_snv, mut sem_dados, mut abaixo_maf, mut snvs) = (0usize, 0usize, 0usize, 0usize, 0usize);
        let mut arredondado = true;
        while l.proxima()? {
            if l.buf.is_empty() || l.buf[0] == b'#' {
                continue;
            }
            linhas += 1;
            divide(&l.buf, &mut sp);
            let Some(r) = linha_painel(caminho, l.n, &l.buf, &sp, &c)? else {
                sem_dados += 1;
                continue;
            };
            if r.ref_.len() != 1 || r.alt.len() != 1 || !e_acgt(r.ref_[0]) || !e_acgt(r.alt[0]) || r.ref_ == r.alt {
                nao_snv += 1;
                continue;
            }
            snvs += 1;
            let media = r.af.iter().sum::<f64>() / k as f64;
            if media.min(1.0 - media) < min_maf {
                abaixo_maf += 1;
                continue;
            }
            let chave = canon.de(r.crom);
            if atual == usize::MAX || chave != chave_atual {
                chave_atual = chave.to_string();
                atual = match indice.get(&chave_atual) {
                    Some(&i) => i,
                    None => {
                        indice.insert(chave_atual.clone(), tmp.len());
                        tmp.push((chave_atual.clone(), Tmp { pos: Vec::new(), ra: Vec::new(), p: Vec::new(), nome: txt(r.crom) }));
                        tmp.len() - 1
                    }
                };
            }
            arredondado &= r.af.iter().all(|&x| duas_casas(x));
            let t = &mut tmp[atual].1;
            t.pos.push(r.pos);
            t.ra.push([r.ref_[0], r.alt[0]]);
            t.p.extend(r.af.iter().map(|&x| x as f32));
        }
        tmp.sort_by(|a, b| compara_crom(&a.0, &b.0));
        let total: usize = tmp.iter().map(|(_, t)| t.pos.len()).sum();
        let mut p = PainelPca {
            grupos: c.grupos,
            cromos: Vec::new(),
            nomes_crom: Vec::new(),
            inicio: vec![0],
            pos: Vec::with_capacity(total),
            ra: Vec::with_capacity(total),
            p: Vec::with_capacity(total * k),
            linhas,
            snvs,
            nao_snv,
            sem_dados,
            abaixo_maf,
            duplicados: 0,
            arredondado: arredondado && total > 0,
        };
        for (nome, t) in tmp {
            let mut ix: Vec<(u64, u32)> = t.pos.iter().enumerate().map(|(i, &x)| (x, i as u32)).collect();
            if !t.pos.windows(2).all(|w| w[0] <= w[1]) {
                multi_merge_sort_stable_by(&mut ix, |a, b| a.0.cmp(&b.0));
            }
            let antes = p.pos.len();
            for &(x, i) in &ix {
                let i = i as usize;
                let dup = (antes..p.pos.len()).rev().take_while(|&j| p.pos[j] == x).any(|j| p.ra[j] == t.ra[i]);
                if dup {
                    p.duplicados += 1;
                    continue;
                }
                p.pos.push(x);
                p.ra.push(t.ra[i]);
                p.p.extend_from_slice(&t.p[i * k..i * k + k]);
            }
            p.cromos.push(nome);
            p.nomes_crom.push(t.nome);
            p.inicio.push(p.pos.len());
        }
        Ok(p)
    }

    fn n(&self) -> usize {
        self.pos.len()
    }
}

// =====================================================================
//  CHOOSING THE SITES
// =====================================================================

/// Left to right over the candidates `ok` accepts, keeping a site when it
/// is at least `gap` bases after the last one kept on its chromosome. For a
/// given gap this keeps the most sites any choice can keep, so the count
/// never grows with the gap nor with a stricter `ok` (the searches below
/// rely on it).
fn espaca(pa: &PainelPca, ok: &dyn Fn(usize) -> bool, gap: u64, mut guarda: Option<&mut Vec<usize>>) -> usize {
    let mut n = 0;
    for c in 0..pa.cromos.len() {
        let mut ultimo: Option<u64> = None;
        for j in pa.inicio[c]..pa.inicio[c + 1] {
            if !ok(j) {
                continue;
            }
            let x = pa.pos[j];
            if matches!(ultimo, Some(u) if x - u < gap) {
                continue;
            }
            ultimo = Some(x);
            n += 1;
            if let Some(v) = guarda.as_mut() {
                v.push(j);
            }
        }
    }
    n
}

/// The sites to use: the most-observed threshold that still offers
/// `max_sites`, then the smallest spacing that keeps no more than that.
/// Returns (sites, people threshold, spacing).
fn escolhe_sitios(pa: &PainelPca, obs: &[u32], pessoas: u32, max: usize, min_gap: u64) -> (Vec<usize>, u32, u64) {
    let conta = |t: u32, gap: u64| espaca(pa, &|j| obs[j] >= t, gap, None);
    let t = if conta(pessoas, min_gap) >= max || pessoas <= 1 {
        pessoas.max(1)
    } else if conta(1, min_gap) < max {
        1
    } else {
        // conta(lo) >= max > conta(hi)
        let (mut lo, mut hi) = (1u32, pessoas);
        while hi - lo > 1 {
            let meio = lo + (hi - lo) / 2;
            if conta(meio, min_gap) >= max {
                lo = meio;
            } else {
                hi = meio;
            }
        }
        lo
    };
    let gap = if conta(t, min_gap) <= max {
        min_gap
    } else {
        // conta(lo) > max >= conta(hi), or hi so large that each chromosome keeps one site.
        let maior = pa.pos.iter().copied().max().unwrap_or(0);
        let (mut lo, mut hi) = (min_gap, maior + 1);
        while hi - lo > 1 {
            let meio = lo + (hi - lo) / 2;
            if conta(t, meio) <= max {
                hi = meio;
            } else {
                lo = meio;
            }
        }
        hi
    };
    let mut sel = Vec::new();
    espaca(pa, &|j| obs[j] >= t, gap, Some(&mut sel));
    if sel.len() > max {
        // More chromosomes than max_sites: one in so many, evenly.
        let n = sel.len();
        sel = (0..max).map(|i| sel[i * n / max]).collect();
    }
    (sel, t, gap)
}

// =====================================================================
//  SIMULATED REFERENCE
// =====================================================================

/// `per_group` individuals per group, row by row (group by group), each
/// genotype Binomial(2, p) drawn independently. Every individual has its own
/// random stream (seed, set, row), so the result does not depend on the
/// number of threads.
fn simula(freqs: &[f32], k: usize, m: usize, por_grupo: usize, semente: u64, conjunto: u64) -> Vec<u8> {
    let n = k * por_grupo;
    let mut g = vec![0u8; n * m];
    g.par_chunks_mut(m).enumerate().for_each(|(r, linha)| {
        let grupo = r / por_grupo;
        let mut s = semente;
        let a = splitmix(&mut s);
        let mut s = a ^ conjunto.wrapping_mul(0xD1B5_4A32_D192_ED03) ^ (r as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        splitmix(&mut s);
        for (j, x) in linha.iter_mut().enumerate() {
            let p = f64::from(freqs[j * k + grupo]);
            let u1 = (splitmix(&mut s) >> 11) as f64 / (1u64 << 53) as f64;
            let u2 = (splitmix(&mut s) >> 11) as f64 / (1u64 << 53) as f64;
            *x = u8::from(u1 < p) + u8::from(u2 < p);
        }
    });
    g
}

// =====================================================================
//  ANCESTRY PCA
// =====================================================================

#[derive(Clone, Default)]
struct Contas {
    genotipados: usize,
    assumidos: usize,
    nao_chamaveis: usize,
    sob_indel: usize,
    outro_alelo: usize,
    ref_diferente: usize,
    sem_gt: usize,
    filtrados: usize,
}

impl Contas {
    fn conta(&mut self, m: Motivo) {
        match m {
            Motivo::Genotipado => self.genotipados += 1,
            Motivo::Assumido => self.assumidos += 1,
            Motivo::NaoChamavel => self.nao_chamaveis += 1,
            Motivo::SobIndel => self.sob_indel += 1,
            Motivo::OutroAlelo => self.outro_alelo += 1,
            Motivo::RefDiferente => self.ref_diferente += 1,
            Motivo::SemGt => self.sem_gt += 1,
            Motivo::Filtrado => self.filtrados += 1,
        }
    }
}

/// Sample names of a VCF's #CHROM line.
fn amostras_do_vcf(vcf: &str) -> Result<Vec<String>, String> {
    let mut l = Linhas::abre(vcf)?;
    while l.proxima()? {
        if l.buf.starts_with(b"#CHROM") {
            let v: Vec<String> = l.buf.split(|&b| b == b'\t').skip(9).map(txt).collect();
            if v.is_empty() {
                return Err(format!("{vcf} has no sample columns (a sites-only VCF): a person needs genotypes."));
            }
            return Ok(v);
        }
        if !l.buf.is_empty() && l.buf[0] != b'#' {
            break;
        }
    }
    Err(format!("{vcf}: no `#CHROM` header line: is this a VCF file?"))
}

fn pessoas_de(vcfs: &[String], o: &PcaOptions) -> Result<Vec<Person>, String> {
    if vcfs.is_empty() {
        return Err("No VCF given: give the people's VCF file(s).".to_string());
    }
    let nao_vazio = |s: &String| if s.trim().is_empty() { None } else { Some(s.clone()) };
    let base: Vec<(String, Option<String>)> = if vcfs.len() == 1 {
        let nomes = match &o.samples {
            Some(s) if !s.is_empty() => s.clone(),
            _ => amostras_do_vcf(&vcfs[0])?,
        };
        nomes.iter().map(|s| (vcfs[0].clone(), nao_vazio(s))).collect()
    } else {
        match &o.samples {
            None => vcfs.iter().map(|v| (v.clone(), None)).collect(),
            Some(s) if s.len() == vcfs.len() => vcfs.iter().zip(s).map(|(v, s)| (v.clone(), nao_vazio(s))).collect(),
            Some(s) => {
                return Err(format!(
                    "`samples` has {} name(s): with several VCFs give one per VCF ({}), \"\" for a file's only sample, or none.",
                    s.len(),
                    vcfs.len()
                ))
            }
        }
    };
    let beds: Vec<Option<String>> = match &o.callable_beds {
        None => vec![None; base.len()],
        Some(b) if b.is_empty() => vec![None; base.len()],
        Some(b) if b.len() == 1 => vec![nao_vazio(&b[0]); base.len()],
        Some(b) if b.len() == base.len() => b.iter().map(nao_vazio).collect(),
        Some(b) => {
            return Err(format!(
                "`callable_beds` has {} file(s): give one for everybody or one per person ({}).",
                b.len(),
                base.len()
            ))
        }
    };
    Ok(base.into_iter().zip(beds).map(|((vcf, sample), callable)| Person { vcf, sample, callable }).collect())
}

/// A PCA of the people against the groups of a frequency panel; see the
/// module notes. `vcfs`, `samples` and `callable_beds` name the people as
/// in `estimate_ancestry`. Returns the result as JSON (see `PcaResult`).
pub fn ancestry_pca(panel: &str, vcfs: &[String], opt: &PcaOptions) -> Result<String, String> {
    let pessoas = pessoas_de(vcfs, opt)?;
    Ok(ancestry_pca_persons(panel, &pessoas, opt)?.json())
}

/// One person's place on the axes.
#[derive(Clone, Debug, PartialEq)]
pub struct PcaPerson {
    pub name: String,
    pub vcf: String,
    /// Empty when too few sites were observed to place the person.
    pub pcs: Vec<f64>,
    pub sites_observed: usize,
    pub missing_fraction: f64,
    /// The group whose cloud centre is closest, and the distance in units of
    /// the clouds' spread (per component, pooled over the groups). The
    /// simulated clouds are tighter than real populations, so real people
    /// sit several spreads away (median 4-7 in the 1000 Genomes).
    pub closest_group: String,
    pub distance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PcaResult {
    pub panel: String,
    pub groups: Vec<String>,
    /// Panel SNVs that passed the MAF filter.
    pub candidate_sites: usize,
    pub sites_used: usize,
    /// A site was used only when observed in at least this many people.
    pub min_people_per_site: usize,
    pub min_maf: f64,
    pub spacing_bp: u64,
    pub per_group: usize,
    pub components: usize,
    pub seed: u64,
    pub variance_explained: Vec<f64>,
    /// Projected / fitted group centres, per component (1 = no shrinkage).
    pub shrinkage: Vec<f64>,
    /// The simulated individuals drawn as clouds: (group index, scores).
    pub reference: Vec<(usize, Vec<f64>)>,
    /// Per group: the projection of 2 x frequency, and the cloud's sd.
    pub centroids: Vec<(Vec<f64>, Vec<f64>)>,
    pub people: Vec<PcaPerson>,
    pub sites_output: Option<String>,
    pub warnings: Vec<String>,
}

/// A number for the JSON: 6 decimals, without trailing zeros.
fn jx(x: f64) -> String {
    if !x.is_finite() {
        return "null".to_string();
    }
    let mut s = format!("{x:.6}");
    if s.contains('.') {
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
    }
    if s == "-0" {
        s = "0".to_string();
    }
    s
}

fn jv(v: &[f64]) -> String {
    format!("[{}]", v.iter().map(|&x| jx(x)).collect::<Vec<_>>().join(","))
}

fn jt(v: &[String]) -> String {
    format!("[{}]", v.iter().map(|s| json_texto(s)).collect::<Vec<_>>().join(","))
}

impl PcaResult {
    pub fn json(&self) -> String {
        let reference: Vec<String> =
            self.reference.iter().map(|(g, pcs)| format!("{{\"group\":{},\"pcs\":{}}}", json_texto(&self.groups[*g]), jv(pcs))).collect();
        let centroids: Vec<String> = self
            .centroids
            .iter()
            .enumerate()
            .map(|(g, (c, sd))| format!("{{\"group\":{},\"pcs\":{},\"sd\":{}}}", json_texto(&self.groups[g]), jv(c), jv(sd)))
            .collect();
        let people: Vec<String> = self
            .people
            .iter()
            .map(|p| {
                format!(
                    "{{\"name\":{},\"vcf\":{},\"pcs\":{},\"sites_observed\":{},\"missing_fraction\":{},\"closest_group\":{},\"distance\":{}}}",
                    json_texto(&p.name),
                    json_texto(&p.vcf),
                    if p.pcs.is_empty() { "null".to_string() } else { jv(&p.pcs) },
                    p.sites_observed,
                    jx(p.missing_fraction),
                    if p.pcs.is_empty() { "null".to_string() } else { json_texto(&p.closest_group) },
                    if p.pcs.is_empty() { "null".to_string() } else { jx(p.distance) }
                )
            })
            .collect();
        format!(
            "{{\"panel\":{},\"groups\":{},\"candidate_sites\":{},\"sites_used\":{},\"min_people_per_site\":{},\"spacing_bp\":{},\"per_group\":{},\"components\":{},\"seed\":{},\"variance_explained\":{},\"shrinkage\":{},\"reference\":[{}],\"centroids\":[{}],\"people\":[{}],\"sites_output\":{},\"warnings\":{},\"text\":{}}}",
            json_texto(&self.panel),
            jt(&self.groups),
            self.candidate_sites,
            self.sites_used,
            self.min_people_per_site,
            self.spacing_bp,
            self.per_group,
            self.components,
            self.seed,
            jv(&self.variance_explained),
            jv(&self.shrinkage),
            reference.join(","),
            centroids.join(","),
            people.join(","),
            self.sites_output.as_deref().map(json_texto).unwrap_or_else(|| "null".to_string()),
            jt(&self.warnings),
            json_texto(&self.resumo())
        )
    }

    pub fn resumo(&self) -> String {
        let k = self.components;
        let mut o = format!(
            "PCA of {} person(s) against the {} groups of the panel ({}): {} sites, chosen evenly (at least {} bp apart) among the {} panel SNVs whose mean frequency over the groups has a minor-allele frequency of at least {}, each observed in at least {} of the {} people.\n",
            self.people.len(),
            self.groups.len(),
            self.groups.join(", "),
            self.sites_used,
            self.spacing_bp,
            self.candidate_sites,
            self.min_maf,
            self.min_people_per_site,
            self.people.len()
        );
        o.push_str(&format!(
            "  Reference clouds: {} individuals per group SIMULATED from the panel's allele frequencies (two independent draws per site). They show the drift between the groups, not real people: they ignore linkage and the structure inside each group, so real people spread more than the clouds.\n",
            self.per_group
        ));
        let ve: Vec<String> = self.variance_explained.iter().enumerate().map(|(c, v)| format!("PC{} {:.1}%", c + 1, 100.0 * v)).collect();
        o.push_str(&format!("  Variance of the reference explained: {}\n", ve.join(", ")));
        let sh: Vec<String> = self.shrinkage.iter().enumerate().map(|(c, v)| format!("PC{} {:.3}", c + 1, v)).collect();
        o.push_str(&format!(
            "  The axes come from one simulated sample; the clouds are a second, independent sample projected exactly as the people are, so both are on the same footing (projection pulls points toward 0 by: {}).\n",
            sh.join(", ")
        ));
        for (g, (c, _)) in self.centroids.iter().enumerate() {
            let v: Vec<String> = c.iter().map(|x| format!("{x:.4}")).collect();
            o.push_str(&format!("  centre of {:<8} {}\n", self.groups[g], v.join("  ")));
        }
        for p in &self.people {
            if p.pcs.is_empty() {
                o.push_str(&format!("  {}: not placed ({} sites observed)\n", p.name, p.sites_observed));
                continue;
            }
            let v: Vec<String> = p.pcs.iter().enumerate().map(|(c, x)| format!("PC{} {x:.4}", c + 1)).collect();
            o.push_str(&format!(
                "  {}: {}; closest cloud {} ({:.1} spreads from its centre); {} of {} sites observed ({:.1}% missing)\n",
                p.name,
                v.join(", "),
                p.closest_group,
                p.distance,
                p.sites_observed,
                self.sites_used,
                100.0 * p.missing_fraction
            ));
        }
        if self.people.iter().any(|p| !p.pcs.is_empty()) {
            o.push_str("  Distances are in units of the simulated clouds' spread. Real populations spread 2 to 4 times more (linkage, structure inside the group, admixture): in the 1000 Genomes (chr22, 5,000 sites) people sit a median of 4 to 7 spreads from their own group's centre. Compare people with each other, not with 1.
");
        }
        if k > self.groups.len().saturating_sub(1) {
            o.push_str(&format!(
                "  With {} groups only the first {} component(s) hold differences between groups; the rest show sampling noise.\n",
                self.groups.len(),
                self.groups.len() - 1
            ));
        }
        for w in &self.warnings {
            o.push_str(&format!("  note: {w}\n"));
        }
        o.push_str("  Positions describe genetic similarity to the reference groups of the panel at these sites; the closest cloud is not an ancestry call.\n");
        o
    }
}

/// `ancestry_pca` over a list of persons.
pub fn ancestry_pca_persons(panel: &str, persons: &[Person], o: &PcaOptions) -> Result<PcaResult, String> {
    if persons.is_empty() {
        return Err("No person given: give at least one (VCF, sample, callable regions).".to_string());
    }
    if o.components == 0 || o.components > 20 {
        return Err(format!("components is {}: give 1 to 20 (4 is usual; the plot uses the first two).", o.components));
    }
    if o.per_group < 2 {
        return Err(format!("per_group is {}: give at least 2 simulated individuals per group (100 is usual).", o.per_group));
    }
    if o.max_sites < o.components + 2 {
        return Err(format!("max_sites is {}: give at least {} (thousands are usual).", o.max_sites, o.components + 2));
    }
    if o.min_maf.is_nan() || !(0.0..0.5).contains(&o.min_maf) {
        return Err(format!("min_maf is {}: give a frequency from 0 to below 0.5 (0.05 is usual).", o.min_maf));
    }
    if let Some(s) = &o.sites_output {
        if s.trim().is_empty() {
            return Err("sites_output is empty: give a file name such as pca_sites.tsv, or none.".to_string());
        }
        let c = s.to_ascii_lowercase();
        if c.ends_with(".parquet") || c.ends_with(".bam") || c.ends_with(".vcf") || c.ends_with(".vcf.gz") {
            return Err(format!("`{s}`: the sites table is written as tab-separated text. Give it a .tsv or .tsv.gz name."));
        }
        let mut entradas: Vec<&str> = vec![panel];
        entradas.extend(persons.iter().map(|p| p.vcf.as_str()));
        entradas.extend(persons.iter().filter_map(|p| p.callable.as_deref()));
        for e in entradas {
            let mesmo = e == s || matches!((std::fs::canonicalize(e), std::fs::canonicalize(s)), (Ok(x), Ok(y)) if x == y);
            if mesmo {
                return Err(format!("sites_output `{s}` is also one of the input files: choose a different path."));
            }
        }
    }
    let pa = PainelPca::le(panel, o.min_maf)?;
    let k = pa.grupos.len();
    if k < 2 {
        return Err(format!("{panel} has a single group ({}): a PCA against reference groups needs at least two.", pa.grupos[0]));
    }
    if pa.n() == 0 {
        return Err(format!(
            "{panel}: no biallelic SNV with a minor-allele frequency of at least {} (mean of the groups) in the panel ({} row(s) read: {} SNVs, {} of them below the MAF; {} not SNVs; {} with a frequency missing).",
            o.min_maf, pa.linhas, pa.snvs, pa.abaixo_maf, pa.nao_snv, pa.sem_dados
        ));
    }
    let n_ref = k * o.per_group;
    if n_ref > MAX_REFERENCIA {
        return Err(format!(
            "{} groups x per_group {} = {n_ref} simulated individuals: at most {MAX_REFERENCIA} (the PCA grows with the cube). Lower per_group to {} or less.",
            k,
            o.per_group,
            MAX_REFERENCIA / k
        ));
    }
    if n_ref.saturating_mul(o.max_sites) > MAX_CELULAS {
        return Err(format!(
            "{n_ref} simulated individuals x max_sites {} is too large for memory: keep the product under {MAX_CELULAS} (e.g. max_sites {}).",
            o.max_sites,
            MAX_CELULAS / n_ref
        ));
    }
    let n_cand = pa.n();
    if persons.len().saturating_mul(n_cand) > MAX_CODIGOS {
        return Err(format!(
            "{} people x {n_cand} panel sites is too large for memory: thin the panel (build_panel(..., thin_bp=20000)), raise min_maf, or analyse fewer people at a time.",
            persons.len()
        ));
    }

    // ---- the people's genotypes at every candidate site
    let regioes = le_regioes_das_pessoas(persons)?;
    let alvos: HashMap<String, &[u64]> =
        pa.cromos.iter().enumerate().map(|(i, c)| (c.clone(), &pa.pos[pa.inicio[i]..pa.inicio[i + 1]])).collect();
    let offset: HashMap<&str, usize> = pa.cromos.iter().enumerate().map(|(i, c)| (c.as_str(), pa.inicio[i])).collect();
    let mut codigos: Vec<Vec<u8>> = vec![Vec::new(); persons.len()];
    let mut contas: Vec<Contas> = vec![Contas::default(); persons.len()];
    let mut nomes: Vec<String> = vec![String::new(); persons.len()];
    for (vcf, idxs) in por_vcf(persons) {
        let pedidas: Vec<Option<String>> = idxs.iter().map(|&i| persons[i].sample.clone()).collect();
        let regs: Vec<Option<&Regioes>> = idxs.iter().map(|&i| persons[i].callable.as_ref().map(|b| &regioes[b])).collect();
        let mut g: Vec<Vec<u8>> = idxs.iter().map(|_| vec![FALTA; n_cand]).collect();
        let mut c: Vec<Contas> = vec![Contas::default(); idxs.len()];
        let (mut crom_atual, mut off) = (String::new(), 0usize);
        let mut f = |crom: &str, t: usize, s: usize, obs: &Obs<'_>| {
            if crom != crom_atual {
                crom_atual = crom.to_string();
                off = offset[crom];
            }
            let j = off + t;
            let [pr, palt] = pa.ra[j];
            let (codigo, motivo) = classifica_snv(obs, s, pr, palt, o.pass_only);
            // Haploid calls (a man's X) count as homozygous.
            g[s][j] = match codigo {
                3 => 0,
                4 => 2,
                x => x,
            };
            c[s].conta(motivo);
        };
        let r = junta(&vcf, &pedidas, &regs, o.absent_is_ref, &alvos, &mut f)?;
        for (s, &i) in idxs.iter().enumerate() {
            codigos[i] = std::mem::take(&mut g[s]);
            contas[i] = c[s].clone();
            nomes[i] = r.nomes[s].clone();
        }
    }
    let mut obs = vec![0u32; n_cand];
    for cod in &codigos {
        for (j, &x) in cod.iter().enumerate() {
            if x != FALTA {
                obs[j] += 1;
            }
        }
    }
    if obs.iter().all(|&x| x == 0) {
        let c = &contas[0];
        return Err(format!(
            "None of the {n_cand} panel SNVs was observed in any person ({}: {} not callable, {} with another allele, {} with a different REF, {} missing or filtered). Check that the VCFs and the panel are on the same genome build and share chromosomes; for a VCF that lists only variants, give the callable-regions BED or keep absent_is_ref on.",
            nomes[0],
            c.nao_chamaveis,
            c.outro_alelo,
            c.ref_diferente,
            c.sem_gt + c.filtrados + c.sob_indel
        ));
    }

    // ---- the sites
    let (mut sel, limiar, gap) = escolhe_sitios(&pa, &obs, persons.len() as u32, o.max_sites, o.min_distance_bp);
    let freqs_de = |sel: &[usize]| -> Vec<f32> { sel.iter().flat_map(|&j| pa.p[j * k..j * k + k].iter().copied()).collect() };
    let mut freqs = freqs_de(&sel);
    let mut ref1 = simula(&freqs, k, sel.len(), o.per_group, o.seed, 1);
    // A site the first sample happens to have monomorphic carries nothing.
    let m0 = sel.len();
    let variavel: Vec<bool> = (0..m0)
        .map(|j| {
            let x0 = ref1[j];
            (1..n_ref).any(|i| ref1[i * m0 + j] != x0)
        })
        .collect();
    let monomorficos = variavel.iter().filter(|v| !**v).count();
    if monomorficos > 0 {
        let manter: Vec<usize> = (0..m0).filter(|&j| variavel[j]).collect();
        let mut novo = Vec::with_capacity(n_ref * manter.len());
        for i in 0..n_ref {
            novo.extend(manter.iter().map(|&j| ref1[i * m0 + j]));
        }
        ref1 = novo;
        sel = manter.iter().map(|&j| sel[j]).collect();
        freqs = freqs_de(&sel);
    }
    let m = sel.len();
    let kc = o.components;
    let minimo = (10 * kc).max(50);
    if m < minimo {
        return Err(format!(
            "Only {m} usable site(s) are shared by the panel and the people (at least {minimo} are needed for {kc} component(s)): {n_cand} panel SNVs pass the MAF filter, and {} of them are observed in someone. Check the genome build and the chromosome names, give the callable-regions BEDs, or lower min_maf.",
            obs.iter().filter(|&&x| x > 0).count()
        ));
    }
    let mut fit = pca_fit(&ref1, n_ref, m, kc)?;
    drop(ref1);

    // ---- centres: the projection of the expected genotypes
    let centro = |fit: &PcaFit, g: usize| -> Vec<f64> {
        let y: Vec<f64> = (0..m).map(|j| (2.0 * f64::from(freqs[j * k + g]) - fit.mean[j]) / fit.sd[j]).collect();
        fit.project_standardized(&y).map(|r| r.0).unwrap_or_else(|| vec![f64::NAN; kc])
    };
    let mut centros: Vec<Vec<f64>> = (0..k).map(|g| centro(&fit, g)).collect();
    // Sign: the group farthest from 0 on an axis is on its positive side.
    for c in 0..kc {
        let g = (0..k).fold(0, |a, g| if centros[g][c].abs() > centros[a][c].abs() { g } else { a });
        if centros[g][c] < 0.0 {
            fit.inverte(c);
            for cg in centros.iter_mut() {
                cg[c] = -cg[c];
            }
        }
    }
    // Shrinkage: projected centres against the fitted sample's own centres.
    let mut ajustados = vec![vec![0.0; kc]; k];
    for i in 0..n_ref {
        for c in 0..kc {
            ajustados[i / o.per_group][c] += fit.scores[i * kc + c] / o.per_group as f64;
        }
    }
    let shrinkage: Vec<f64> = (0..kc)
        .map(|c| {
            let num: f64 = (0..k).map(|g| centros[g][c] * ajustados[g][c]).sum();
            let den: f64 = (0..k).map(|g| ajustados[g][c] * ajustados[g][c]).sum();
            if den > 0.0 {
                num / den
            } else {
                f64::NAN
            }
        })
        .collect();

    // ---- the clouds: a second, independent sample, projected
    let ref2 = simula(&freqs, k, m, o.per_group, o.seed, 2);
    let nuvem: Vec<Vec<f64>> = ref2
        .par_chunks(m)
        .map(|linha| fit.project(linha).map(|r| r.0).unwrap_or_else(|| vec![f64::NAN; kc]))
        .collect();
    drop(ref2);
    let mut dp = vec![vec![0.0; kc]; k];
    for g in 0..k {
        for c in 0..kc {
            let v: Vec<f64> = (0..o.per_group).map(|i| nuvem[g * o.per_group + i][c]).collect();
            let media = v.iter().sum::<f64>() / v.len() as f64;
            dp[g][c] = (v.iter().map(|x| (x - media) * (x - media)).sum::<f64>() / (v.len() as f64 - 1.0)).sqrt();
        }
    }
    let dp_comum: Vec<f64> = (0..kc).map(|c| ((0..k).map(|g| dp[g][c] * dp[g][c]).sum::<f64>() / k as f64).sqrt().max(1e-12)).collect();

    // ---- the people
    let min_pessoa = (5 * kc).max(20);
    let pessoas_pca: Vec<PcaPerson> = persons
        .par_iter()
        .enumerate()
        .map(|(i, p)| {
            let g: Vec<u8> = sel.iter().map(|&j| codigos[i][j]).collect();
            let vistos = g.iter().filter(|&&x| x != FALTA).count();
            let proj = if vistos >= min_pessoa { fit.project(&g) } else { None };
            let (pcs, grupo, dist) = match proj {
                Some((t, _)) => {
                    let d: Vec<f64> = centros
                        .iter()
                        .map(|cg| (0..kc).map(|c| ((t[c] - cg[c]) / dp_comum[c]).powi(2)).sum::<f64>().sqrt())
                        .collect();
                    let melhor = (0..k).fold(0, |a, g| if d[g] < d[a] { g } else { a });
                    (t, pa.grupos[melhor].clone(), d[melhor])
                }
                None => (Vec::new(), String::new(), f64::NAN),
            };
            PcaPerson {
                name: nomes[i].clone(),
                vcf: p.vcf.clone(),
                pcs,
                sites_observed: vistos,
                missing_fraction: 1.0 - vistos as f64 / m as f64,
                closest_group: grupo,
                distance: dist,
            }
        })
        .collect();
    if pessoas_pca.iter().all(|p| p.pcs.is_empty()) {
        return Err(format!(
            "No person could be placed: each needs at least {min_pessoa} of the {m} sites observed. Check the genome build and the chromosome names, or give the callable-regions BEDs."
        ));
    }

    // ---- notes
    let mut avisos = Vec::new();
    if m * 10 < o.max_sites * 9 {
        avisos.push(format!(
            "{m} sites were usable, fewer than max_sites ({}): the positions are less precise (a panel with more sites, or people observed at more of them, helps)",
            o.max_sites
        ));
    }
    if (limiar as usize) < persons.len() {
        avisos.push(format!(
            "too few sites are observed in all {} people: sites observed in at least {limiar} were used, so people differ in the sites they are placed with",
            persons.len()
        ));
    }
    if kc > k - 1 {
        avisos.push(format!(
            "{kc} components were asked for, but {k} groups differ along only {} direction(s): PC{} onward show only the simulated sampling noise",
            k - 1,
            k
        ));
    }
    if monomorficos > 0 {
        avisos.push(format!("{monomorficos} chosen site(s) came out with a single genotype in the simulated reference and were dropped"));
    }
    for (i, p) in pessoas_pca.iter().enumerate() {
        if p.pcs.is_empty() {
            avisos.push(format!("{}: only {} of the {m} sites observed (at least {min_pessoa} are needed): not placed", p.name, p.sites_observed));
        } else if p.missing_fraction > 0.5 {
            avisos.push(format!(
                "{}: {:.0}% of the sites are missing; the position is projected on the observed sites only (not pulled toward 0) but is less precise",
                p.name,
                100.0 * p.missing_fraction
            ));
        }
        let c = &contas[i];
        let vistos = c.genotipados + c.ref_diferente + c.outro_alelo;
        if vistos > 100 && c.ref_diferente * 20 > vistos {
            avisos.push(format!(
                "{}: the REF of the VCF differs from the panel's at {} of {} positions found in both: are they on the same genome build (GRCh37 vs GRCh38)?",
                p.name, c.ref_diferente, vistos
            ));
        }
        if persons[i].callable.is_none() && o.absent_is_ref && c.assumidos > 0 && c.genotipados * 2 < c.assumidos {
            avisos.push(format!(
                "{}: {} panel sites absent from the VCF were taken as homozygous reference (absent_is_ref, no callable-regions BED): right for a whole-genome VCF, wrong where the genome was not assessed",
                p.name, c.assumidos
            ));
        }
    }
    if pa.arredondado {
        avisos.push("the panel's frequencies have at most 2 decimals: the clouds of groups where an allele is rare are drawn roughly".to_string());
    }
    if pa.duplicados > 0 || pa.sem_dados > 0 {
        avisos.push(format!("panel rows left out: {} duplicate site(s), {} with a frequency missing", pa.duplicados, pa.sem_dados));
    }

    // ---- the sites table
    if let Some(saida) = &o.sites_output {
        let mut w = Escritor::cria(saida)?;
        let r = (|| -> Result<(), String> {
            w.linha(&format!(
                "##mmgen_ancestry_pca panel={panel} per_group={} components={kc} seed={} sites={m}",
                o.per_group, o.seed
            ))?;
            w.linha("##mean and sd standardize a genotype (ALT alleles): (g - mean) / sd; PCk is the unit-length axis; a score is sum((g - mean) / sd * PCk) / sqrt(sites)")?;
            let mut cab = "chrom\tpos\tref\talt".to_string();
            for g in &pa.grupos {
                cab.push_str(&format!("\tAF_{g}"));
            }
            cab.push_str("\tmean\tsd");
            for c in 0..kc {
                cab.push_str(&format!("\tPC{}", c + 1));
            }
            w.linha(&cab)?;
            for (jj, &j) in sel.iter().enumerate() {
                let c = pa.inicio.partition_point(|&x| x <= j) - 1;
                let mut l = format!("{}\t{}\t{}\t{}", pa.nomes_crom[c], pa.pos[j], pa.ra[j][0] as char, pa.ra[j][1] as char);
                for g in 0..k {
                    l.push_str(&format!("\t{}", jx(f64::from(pa.p[j * k + g]))));
                }
                l.push_str(&format!("\t{}\t{}", jx(fit.mean[jj]), jx(fit.sd[jj])));
                for c in 0..kc {
                    l.push_str(&format!("\t{:.6e}", fit.loadings[jj * kc + c]));
                }
                w.linha(&l)?;
            }
            Ok(())
        })();
        let fim = w.termina();
        if r.is_err() || fim.is_err() {
            let _ = std::fs::remove_file(saida);
            r?;
            fim?;
        }
    }

    let reference: Vec<(usize, Vec<f64>)> = nuvem.into_iter().enumerate().map(|(i, v)| (i / o.per_group, v)).collect();
    Ok(PcaResult {
        panel: panel.to_string(),
        groups: pa.grupos.clone(),
        candidate_sites: n_cand,
        sites_used: m,
        min_people_per_site: limiar as usize,
        min_maf: o.min_maf,
        spacing_bp: gap,
        per_group: o.per_group,
        components: kc,
        seed: o.seed,
        variance_explained: fit.eigenvalues.iter().map(|mu| mu / fit.total_variance).collect(),
        shrinkage,
        reference,
        centroids: centros.into_iter().zip(dp).collect(),
        people: pessoas_pca,
        sites_output: o.sites_output.clone(),
        warnings: avisos,
    })
}
