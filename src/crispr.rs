// =====================================================================
//  POOLED CRISPR KNOCKOUT SCREENS: COUNTING GUIDES AND RANKING GENES
//
//  The two steps of a pooled CRISPR screen analysis, done the way MAGeCK
//  does them (Li et al. 2014, Genome Biology 15:554; MAGeCK 0.5.9.5):
//
//   crispr_count  reads one FASTQ per sample (several = technical
//                 replicates, summed) and counts the reads that carry each
//                 guide of the library -- `mageck count`.
//   crispr_test   compares treatment with control samples of a count table
//                 -- `mageck test` with robust rank aggregation (RRA):
//                 normalization, a mean-variance model fitted across guides,
//                 a p-value and a log2 fold change per guide, then alpha-RRA
//                 per gene for negative and positive selection, with
//                 permutation p-values and Benjamini-Hochberg FDR.
//
//  FIDELITY. The test reproduces MAGeCK step by step, quirks included, so
//  the two give the same numbers (compared in trials/src/bin/crispr.rs):
//
//   - the "mean" of a group of samples is the MEDIAN (MAGeCK >= 0.5.7);
//   - median normalization: size factor = 1 / upper median of
//     count / exp(mean(ln(count + 1))), over guides with any read; it falls
//     back to total-count normalization when a sample has more than 45%
//     zero-count guides or a zero median ratio;
//   - the variance model is var = mean + mean^k * 2^b, with (k, b) a least
//     squares line of log2(var - mean + 1) on log2(mean + 1) weighted by the
//     mean, over the guides whose variance exceeds the mean (high-count
//     outliers left out); k >= 1 and b >= 0;
//   - p.low / p.high come from a normal distribution truncated at zero; the
//     guide score is qnorm(p.low), read from MAGeCK's table of ln p in steps
//     of 0.001 and interpolated linearly. Outside ln p in [-24, -0.001] it
//     is the plain z-score, as in MAGeCK (so the score is not monotone in
//     p.low at p.low ~ 4e-11; this is MAGeCK's behaviour, kept on purpose);
//   - RRA: percentiles by mid-rank among the guides that passed the zero
//     filter; lo-value = min over the sorted percentiles below alpha (the
//     first one always counts) of the beta CDF; alpha = fraction of guides
//     with p <= 0.25 (clamped to [0.05, 0.5] for negative selection);
//     permutation per gene size with the same Park-Miller generator and
//     seed as MAGeCK's RRA (123456), so the permutation p-values match too.
//
//  On MAGeCK's demo1 (four normalization modes) every line of the two
//  summary files comes out identical to MAGeCK's. Where order can differ:
//  genes with exactly the same lo-value keep their order of first appearance
//  here (RRA sorts them with an unstable quicksort), and `crispr_count`
//  writes the count table in library order (MAGeCK in the order the reads
//  arrive), which reorders tied guides downstream.
//
//  ANY SIZE. Counting streams the FASTQ files (text, .gz/BGZF, .lz4); the
//  memory is the library (one counter per guide and sample) plus a batch of
//  reads. The test holds the count table, one row per guide.
//
//  Every error names the file (and the line) and says how to fix it.
// =====================================================================

use crate::externo::{abre_leitura, Escritor};
use crate::multimerge_core::multi_merge_sort_stable_by;
use crate::nomes::json_texto;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{BufRead, Write};
use std::sync::OnceLock;

/// Longest guide accepted (two bits per base in a u128).
pub const MAX_GUIDE_LENGTH: usize = 64;
/// Reads used to find the guide offset in each FASTQ, as in MAGeCK.
pub const OFFSET_PROBE_READS: usize = 100_000;
/// Reads matched together (in parallel) after the probe.
const LOTE: usize = 262_144;

// ---------------------------------------------------------------------
//  SMALL HELPERS
// ---------------------------------------------------------------------

/// A, C, G, T (either case) -> 0..3; anything else -> 4.
const fn tabela_bases() -> [u8; 256] {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
}
static BASE: [u8; 256] = tabela_bases();

fn json_num(x: f64) -> String {
    if !x.is_finite() {
        return "null".to_string();
    }
    let a = x.abs();
    if a != 0.0 && !(1e-5..1e16).contains(&a) {
        format!("{x:e}")
    } else {
        format!("{x}")
    }
}

fn json_lista<T>(v: &[T], f: impl Fn(&T) -> String) -> String {
    let mut o = String::with_capacity(v.len() * 8 + 2);
    o.push('[');
    for (i, x) in v.iter().enumerate() {
        if i > 0 {
            o.push(',');
        }
        o.push_str(&f(x));
    }
    o.push(']');
    o
}

fn json_textos(v: &[String]) -> String {
    json_lista(v, |s| json_texto(s))
}

fn json_nums(v: &[f64]) -> String {
    json_lista(v, |x| json_num(*x))
}

/// Builds a JSON object field by field.
struct Objeto {
    s: String,
    algum: bool,
}

impl Objeto {
    fn novo() -> Objeto {
        Objeto { s: String::from("{"), algum: false }
    }
    fn campo(&mut self, k: &str, v: impl AsRef<str>) -> &mut Objeto {
        if self.algum {
            self.s.push(',');
        }
        self.algum = true;
        self.s.push_str(&json_texto(k));
        self.s.push(':');
        self.s.push_str(v.as_ref());
        self
    }
    fn fim(&mut self) -> String {
        self.s.push('}');
        std::mem::take(&mut self.s)
    }
}

/// Python's `"{:.Ng}".format(x)` -- the number format of MAGeCK's tables.
pub fn py_g(x: f64, digitos: usize) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".to_string() } else { "-inf".to_string() };
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0".to_string() } else { "0".to_string() };
    }
    let p = digitos.max(1);
    let e = format!("{:.*e}", p - 1, x);
    let (m, ex) = e.split_once('e').unwrap_or((&e, "0"));
    let ex: i32 = ex.parse().unwrap_or(0);
    if ex < -4 || ex >= p as i32 {
        format!("{}e{}{:02}", sem_zeros(m), if ex < 0 { '-' } else { '+' }, ex.abs())
    } else {
        sem_zeros(&format!("{:.*}", (p as i32 - 1 - ex) as usize, x))
    }
}

fn sem_zeros(s: &str) -> String {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s.to_string()
    }
}

/// Python's `str(x)` of a float: shortest digits that read back the same.
pub fn py_repr(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".to_string() } else { "-inf".to_string() };
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0.0".to_string() } else { "0.0".to_string() };
    }
    let e = format!("{x:e}");
    let (m, ex) = e.split_once('e').unwrap_or((&e, "0"));
    let ex: i32 = ex.parse().unwrap_or(0);
    let sinal = if m.starts_with('-') { "-" } else { "" };
    let d: String = m.chars().filter(|c| c.is_ascii_digit()).collect();
    if (-4..16).contains(&ex) {
        if ex >= 0 {
            let ex = ex as usize;
            if d.len() > ex + 1 {
                format!("{sinal}{}.{}", &d[..ex + 1], &d[ex + 1..])
            } else {
                format!("{sinal}{}{}.0", d, "0".repeat(ex + 1 - d.len()))
            }
        } else {
            format!("{sinal}0.{}{}", "0".repeat((-ex - 1) as usize), d)
        }
    } else {
        let mm = if d.len() > 1 { format!("{}.{}", &d[..1], &d[1..]) } else { d.clone() };
        format!("{sinal}{mm}e{}{:02}", if ex < 0 { '-' } else { '+' }, ex.abs())
    }
}

/// What MAGeCK writes for a number RRA printed with `%10.4e` and Python read
/// back: five significant digits, in Python's `str`.
fn fmt_e4(x: f64) -> String {
    py_repr(format!("{x:.4e}").parse().unwrap_or(x))
}

/// The same for a number RRA printed with `%f` (six decimals).
fn fmt_f6(x: f64) -> String {
    py_repr(format!("{x:.6}").parse().unwrap_or(x))
}

/// A piece of the file, for a message.
fn mostra(s: &str, max: usize) -> String {
    if s.chars().count() > max {
        format!("{}...", s.chars().take(max).collect::<String>())
    } else {
        s.to_string()
    }
}

fn mesmo_arquivo(a: &str, b: &str) -> bool {
    a == b || matches!((std::fs::canonicalize(a), std::fs::canonicalize(b)), (Ok(x), Ok(y)) if x == y)
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

/// Line reader that remembers the path and the line number for messages.
struct Linhas {
    caminho: String,
    r: Box<dyn BufRead>,
    buf: String,
    linha: usize,
}

impl Linhas {
    fn abre(caminho: &str) -> Result<Linhas, String> {
        Ok(Linhas { caminho: caminho.to_string(), r: abre_leitura(caminho)?, buf: String::new(), linha: 0 })
    }

    fn avanca(&mut self) -> Result<bool, String> {
        self.buf.clear();
        match self.r.read_line(&mut self.buf) {
            Ok(0) => Ok(false),
            Ok(_) => {
                self.linha += 1;
                Ok(true)
            }
            Err(e) => Err(format!(
                "Error reading {} after line {}: {e}. The file may be truncated, or not a text table (a spreadsheet .xlsx must be saved as .csv or as tab-separated .txt first).",
                self.caminho, self.linha
            )),
        }
    }

    /// The current line without the line ending and without a UTF-8 BOM.
    fn texto(&self) -> &str {
        self.buf.trim_end_matches(['\n', '\r']).trim_start_matches('\u{feff}')
    }
}

/// The fields of a table line: trimmed, without surrounding double quotes.
fn campos(linha: &str, sep: char) -> Vec<&str> {
    linha
        .split(sep)
        .map(|c| {
            let c = c.trim();
            c.strip_prefix('"').and_then(|c| c.strip_suffix('"')).unwrap_or(c)
        })
        .collect()
}

/// Tab if the line has one, else comma if it has one.
fn separador(linha: &str) -> Option<char> {
    if linha.contains('\t') {
        Some('\t')
    } else if linha.contains(',') {
        Some(',')
    } else {
        None
    }
}

fn so_acgt(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
}

fn complemento_reverso(s: &str) -> String {
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

fn ordena_f64(v: &mut [f64]) {
    multi_merge_sort_stable_by(v, |a: &f64, b: &f64| a.partial_cmp(b).unwrap_or(Ordering::Equal));
}

/// Indices 0..n sorted by `chave`, ascending or descending; ties keep their
/// order (Python's `sorted`, with or without `reverse=True`).
fn ordem_por(chave: &[f64], crescente: bool) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..chave.len()).collect();
    if crescente {
        multi_merge_sort_stable_by(&mut idx, |&a: &usize, &b: &usize| chave[a].partial_cmp(&chave[b]).unwrap_or(Ordering::Equal));
    } else {
        multi_merge_sort_stable_by(&mut idx, |&a: &usize, &b: &usize| chave[b].partial_cmp(&chave[a]).unwrap_or(Ordering::Equal));
    }
    idx
}

/// MAGeCK's `mmedian`: the middle value, or the mean of the two middle
/// values; 0 for an empty list.
pub fn median(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let mut s = v.to_vec();
    ordena_f64(&mut s);
    let n = s.len();
    let i = (n - 1) / 2;
    if n % 2 == 1 {
        s[i]
    } else {
        (s[i] + s[i + 1]) / 2.0
    }
}

/// MAGeCK's `getVars`: squared deviations from `centro` (the median) over
/// n - 1 (or 1).
fn variancia(v: &[f64], centro: f64) -> f64 {
    let s: f64 = v.iter().map(|&x| (x - centro) * (x - centro)).sum();
    s / ((v.len() as f64) - 1.0).max(1.0)
}

/// `math.log(x, 2)` as CPython computes it: ln(x) / ln(2).
fn log2_py(x: f64) -> f64 {
    x.ln() / 2f64.ln()
}

/// A list of names, one per line (first column of a tab- or comma-separated
/// file); blank lines are skipped. For the control sgRNA and control gene
/// lists of `crispr_test`.
pub fn read_name_list(caminho: &str) -> Result<Vec<String>, String> {
    let mut l = Linhas::abre(caminho)?;
    let mut v = Vec::new();
    while l.avanca()? {
        let t = l.texto();
        let nome = t.split(['\t', ',']).next().unwrap_or("").trim();
        if !nome.is_empty() {
            v.push(nome.to_string());
        }
    }
    if v.is_empty() {
        return Err(format!("{caminho}: the list is empty. Write one name per line."));
    }
    Ok(v)
}

// =====================================================================
//  PART 1: THE LIBRARY
// =====================================================================

/// A guide library: id, guide (spacer) sequence and gene of each guide.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Library {
    pub path: String,
    pub ids: Vec<String>,
    /// Upper case; reverse-complemented when `reverse_complement` is on.
    pub sequences: Vec<String>,
    pub genes: Vec<String>,
    /// The first line was a header (its second column is not a DNA
    /// sequence), and was skipped.
    pub header: bool,
    /// `\t` or `,`.
    pub separator: char,
    pub reverse_complement: bool,
    /// Lines whose guide id was already seen: skipped, as in MAGeCK.
    pub duplicate_ids: usize,
    /// Guides whose sequence was already in the library: skipped (the first
    /// one keeps all the reads), as in MAGeCK.
    pub duplicate_sequences: usize,
    /// Up to 5 of the skipped duplicates, as `id (line N): ...`.
    pub duplicate_examples: Vec<String>,
    /// (guide length, guides of that length), longest first.
    pub lengths: Vec<(usize, usize)>,
}

impl Library {
    /// Distinct genes.
    pub fn gene_count(&self) -> usize {
        self.genes.iter().collect::<HashSet<_>>().len()
    }

    pub fn json(&self) -> String {
        Objeto::novo()
            .campo("path", json_texto(&self.path))
            .campo("guides", self.ids.len().to_string())
            .campo("genes", self.gene_count().to_string())
            .campo("header", self.header.to_string())
            .campo("separator", json_texto(&self.separator.to_string()))
            .campo("reverse_complement", self.reverse_complement.to_string())
            .campo("duplicate_ids", self.duplicate_ids.to_string())
            .campo("duplicate_sequences", self.duplicate_sequences.to_string())
            .campo("duplicate_examples", json_textos(&self.duplicate_examples))
            .campo("lengths", json_lista(&self.lengths, |(l, n)| format!("[{l},{n}]")))
            .fim()
    }

    pub fn resumo(&self) -> String {
        let comp: Vec<String> = self.lengths.iter().map(|(l, n)| format!("{n} of {l} nt")).collect();
        let mut o = format!(
            "{}: {} guides, {} genes ({}){}{}\n",
            self.path,
            self.ids.len(),
            self.gene_count(),
            comp.join(", "),
            if self.header { "; header line skipped" } else { "" },
            if self.reverse_complement { "; sequences reverse-complemented" } else { "" }
        );
        if self.duplicate_ids > 0 || self.duplicate_sequences > 0 {
            o.push_str(&format!(
                "  skipped: {} line(s) with a repeated guide id, {} guide(s) with a repeated sequence (the first guide with that sequence gets the reads)\n",
                self.duplicate_ids, self.duplicate_sequences
            ));
            for e in &self.duplicate_examples {
                o.push_str(&format!("    {e}\n"));
            }
        }
        o
    }
}

/// Reads a guide library: three columns -- guide id, guide sequence, gene --
/// separated by tabs or commas (detected), with or without a header line
/// (detected: the first line is a header when its second column is not a
/// DNA sequence). Extra columns are ignored. Text, `.gz` or `.lz4`.
///
/// As in MAGeCK, a repeated guide id or a repeated sequence is skipped and
/// counted. `reverse_complement` reverse-complements the library sequences,
/// for reads that carry the guide on the other strand.
pub fn read_library(caminho: &str, reverse_complement: bool) -> Result<Library, String> {
    let mut l = Linhas::abre(caminho)?;
    let mut lib = Library { path: caminho.to_string(), reverse_complement, ..Default::default() };
    let mut sep: Option<char> = None;
    let mut primeira = true;
    let mut ids: HashSet<String> = HashSet::new();
    let mut por_seq: HashMap<String, usize> = HashMap::new();
    while l.avanca()? {
        let t = l.texto();
        if t.trim().is_empty() {
            continue;
        }
        let s = match sep {
            Some(s) => s,
            None => match separador(t) {
                Some(s) => {
                    sep = Some(s);
                    s
                }
                None => {
                    return Err(format!(
                        "{caminho}, line {}: no tab or comma in `{}`. A library has three columns -- guide id, guide sequence, gene -- separated by tabs (.txt) or commas (.csv).",
                        l.linha,
                        mostra(t, 60)
                    ))
                }
            },
        };
        let c = campos(t, s);
        if primeira {
            primeira = false;
            if c.len() >= 2 && !so_acgt(c[1]) {
                lib.header = true;
                continue;
            }
        }
        if c.len() < 3 {
            return Err(format!(
                "{caminho}, line {}: {} column(s) in `{}`; a library needs 3 -- guide id, guide sequence, gene, in this order.",
                l.linha,
                c.len(),
                mostra(t, 60)
            ));
        }
        let (id, seq, gene) = (c[0], c[1], c[2]);
        if id.is_empty() {
            return Err(format!("{caminho}, line {}: the guide id (first column) is empty. Every guide needs a name.", l.linha));
        }
        if !so_acgt(seq) {
            return Err(format!(
                "{caminho}, line {}: the sequence `{}` of guide `{id}` has characters other than A, C, G and T. The second column must be the guide (spacer) sequence: check the column order (id, sequence, gene).",
                l.linha,
                mostra(seq, 40)
            ));
        }
        if seq.len() > MAX_GUIDE_LENGTH {
            return Err(format!(
                "{caminho}, line {}: guide `{id}` is {} nt long; guides up to {MAX_GUIDE_LENGTH} nt are supported. Put only the spacer in the second column, without the scaffold.",
                l.linha,
                seq.len()
            ));
        }
        if gene.is_empty() {
            return Err(format!(
                "{caminho}, line {}: guide `{id}` has no gene (third column). Give every guide a gene; non-targeting controls usually get a name such as `NonTargeting`.",
                l.linha
            ));
        }
        if !ids.insert(id.to_string()) {
            lib.duplicate_ids += 1;
            if lib.duplicate_examples.len() < 5 {
                lib.duplicate_examples.push(format!("{id} (line {}): the guide id appears again", l.linha));
            }
            continue;
        }
        let mut sq = seq.to_ascii_uppercase();
        if reverse_complement {
            sq = complemento_reverso(&sq);
        }
        if let Some(&k) = por_seq.get(&sq) {
            lib.duplicate_sequences += 1;
            if lib.duplicate_examples.len() < 5 {
                lib.duplicate_examples.push(format!("{id} (line {}): same sequence as {}", l.linha, lib.ids[k]));
            }
            continue;
        }
        por_seq.insert(sq.clone(), lib.ids.len());
        lib.ids.push(id.to_string());
        lib.sequences.push(sq);
        lib.genes.push(gene.to_string());
    }
    if lib.ids.is_empty() {
        return Err(format!(
            "{caminho}: no guides found. A library has one guide per line: guide id, guide sequence, gene (tab- or comma-separated)."
        ));
    }
    lib.separator = sep.unwrap_or('\t');
    let mut comp: BTreeMap<usize, usize> = BTreeMap::new();
    for s in &lib.sequences {
        *comp.entry(s.len()).or_insert(0) += 1;
    }
    lib.lengths = comp.into_iter().rev().collect();
    Ok(lib)
}

// =====================================================================
//  PART 2: MATCHING READS TO GUIDES
// =====================================================================

#[derive(Default)]
struct Mistura(u64);

impl Hasher for Mistura {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(5) ^ b as u64).wrapping_mul(0x517c_c1b7_2722_0a95);
        }
    }
    fn write_u128(&mut self, x: u128) {
        let mut h = (x as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ ((x >> 64) as u64).rotate_left(32).wrapping_mul(0xC2B2_AE3D_27D4_EB4F);
        h ^= h >> 29;
        h = h.wrapping_mul(0xBF58_476D_1CE4_E5B9);
        h ^= h >> 32;
        self.0 = h;
    }
}

type Mapa = HashMap<u128, u32, BuildHasherDefault<Mistura>>;

/// Two or more guides share this key.
const AMBIGUO: u32 = u32::MAX;

/// Two bits per base; how many non-ACGT bases and where the last one is
/// (encoded as A).
fn codifica(s: &[u8]) -> (u128, u32, usize) {
    let mut c = 0u128;
    let (mut n, mut pos) = (0u32, 0usize);
    for (i, &b) in s.iter().enumerate() {
        let v = BASE[b as usize];
        c <<= 2;
        if v == 4 {
            n += 1;
            pos = i;
        } else {
            c |= v as u128;
        }
    }
    (c, n, pos)
}

fn mascara(c: u128, len: usize, i: usize) -> u128 {
    c & !(3u128 << (2 * (len - 1 - i)))
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Acerto {
    Nada,
    Exato(u32),
    Troca(u32),
    Ambiguo,
}

/// Exact lookup per guide length, and -- for one mismatch -- one map per
/// position with that base masked: a window one mismatch away from a guide
/// meets it under the mask of the position where they differ.
struct Indice {
    /// Guide lengths, longest first (MAGeCK tries the longest first).
    classes: Vec<usize>,
    exatos: Vec<Mapa>,
    vizinhos: Vec<Vec<Mapa>>,
}

impl Indice {
    fn novo(lib: &Library, mismatches: u8) -> Indice {
        let classes: Vec<usize> = lib.lengths.iter().map(|&(l, _)| l).collect();
        let mut exatos: Vec<Mapa> = classes.iter().map(|_| Mapa::default()).collect();
        let mut vizinhos: Vec<Vec<Mapa>> = classes
            .iter()
            .map(|&l| if mismatches > 0 { (0..l).map(|_| Mapa::default()).collect() } else { Vec::new() })
            .collect();
        for (g, s) in lib.sequences.iter().enumerate() {
            let ci = classes.iter().position(|&l| l == s.len()).unwrap_or(0);
            let (c, _, _) = codifica(s.as_bytes());
            exatos[ci].insert(c, g as u32);
            if mismatches > 0 {
                for (i, mapa) in vizinhos[ci].iter_mut().enumerate() {
                    let e = mapa.entry(mascara(c, s.len(), i)).or_insert(g as u32);
                    if *e != g as u32 {
                        *e = AMBIGUO;
                    }
                }
            }
        }
        Indice { classes, exatos, vizinhos }
    }

    /// First the exact match, trying the offsets in order and, at each
    /// offset, the longest guides first (MAGeCK's order); windows with an N
    /// are skipped. Then, with `troca`, one mismatch (an N counts as the
    /// mismatch): the read goes to a guide only when that guide is the only
    /// one within one mismatch over all the windows tried.
    fn casa(&self, seq: &[u8], offsets: &[usize], troca: bool) -> Acerto {
        for &o in offsets {
            for (ci, &len) in self.classes.iter().enumerate() {
                if o + len > seq.len() {
                    continue;
                }
                let (c, n, _) = codifica(&seq[o..o + len]);
                if n > 0 {
                    continue;
                }
                if let Some(&g) = self.exatos[ci].get(&c) {
                    return Acerto::Exato(g);
                }
            }
        }
        if !troca {
            return Acerto::Nada;
        }
        let mut achado: Option<u32> = None;
        for &o in offsets {
            for (ci, &len) in self.classes.iter().enumerate() {
                if o + len > seq.len() {
                    continue;
                }
                let (c, n, pos) = codifica(&seq[o..o + len]);
                if n > 1 {
                    continue;
                }
                let (de, ate) = if n == 1 { (pos, pos + 1) } else { (0, len) };
                for i in de..ate {
                    if let Some(&g) = self.vizinhos[ci][i].get(&mascara(c, len, i)) {
                        if g == AMBIGUO {
                            return Acerto::Ambiguo;
                        }
                        match achado {
                            None => achado = Some(g),
                            Some(h) if h != g => return Acerto::Ambiguo,
                            _ => {}
                        }
                    }
                }
            }
        }
        match achado {
            Some(g) => Acerto::Troca(g),
            None => Acerto::Nada,
        }
    }

    /// The first offset (from 0) where the read has a library guide,
    /// exactly -- MAGeCK's `--trim-5 AUTO` probe.
    fn primeiro_offset(&self, seq: &[u8]) -> Option<usize> {
        let menor = *self.classes.last()?;
        if seq.len() < menor {
            return None;
        }
        for o in 0..=(seq.len() - menor) {
            for (ci, &len) in self.classes.iter().enumerate() {
                if seq.len() - o < len {
                    continue;
                }
                let (c, n, _) = codifica(&seq[o..o + len]);
                if n > 0 {
                    continue;
                }
                if self.exatos[ci].contains_key(&c) {
                    return Some(o);
                }
            }
        }
        None
    }
}

// ---------------------------------------------------------------------
//  FASTQ, streamed in batches
// ---------------------------------------------------------------------

struct Fastq {
    r: Box<dyn BufRead>,
    caminho: String,
    linha: u64,
    registros: u64,
    cab: Vec<u8>,
    seq: Vec<u8>,
    sep: Vec<u8>,
    qual: Vec<u8>,
}

impl Fastq {
    fn abre(caminho: &str) -> Result<Fastq, String> {
        Ok(Fastq {
            r: abre_leitura(caminho)?,
            caminho: caminho.to_string(),
            linha: 0,
            registros: 0,
            cab: Vec::new(),
            seq: Vec::new(),
            sep: Vec::new(),
            qual: Vec::new(),
        })
    }

    fn le(&mut self, qual: u8) -> Result<bool, String> {
        let buf = match qual {
            0 => &mut self.cab,
            1 => &mut self.seq,
            2 => &mut self.sep,
            _ => &mut self.qual,
        };
        buf.clear();
        let n = self.r.read_until(b'\n', buf).map_err(|e| {
            format!(
                "Error reading {} after line {}: {e}. The file may be corrupted or incomplete (an interrupted download leaves this behind); fetch it again.",
                self.caminho, self.linha
            )
        })?;
        if n == 0 {
            return Ok(false);
        }
        self.linha += 1;
        if buf.last() == Some(&b'\n') {
            buf.pop();
        }
        if buf.last() == Some(&b'\r') {
            buf.pop();
        }
        Ok(true)
    }

    fn truncado(&self, faltou: &str) -> String {
        format!(
            "{}: the file ends at line {}, in the middle of record {}: expected {faltou} next. An interrupted download or copy usually leaves this behind; fetch the file again.",
            self.caminho,
            self.linha,
            self.registros + 1
        )
    }

    /// The next record into `seq`. `false` at the end of the file.
    fn proximo(&mut self) -> Result<bool, String> {
        loop {
            if !self.le(0)? {
                return Ok(false);
            }
            if !self.cab.is_empty() {
                break;
            }
        }
        if self.cab[0] != b'@' {
            let dica = if self.cab[0] == b'>' {
                " The file looks like FASTA; counting needs the FASTQ reads."
            } else {
                " A FASTQ record has exactly 4 lines: @name, sequence, +, qualities."
            };
            return Err(format!(
                "{}, line {}: expected the header of record {}, which starts with `@`, but the line is `{}`.{dica}",
                self.caminho,
                self.linha,
                self.registros + 1,
                mostra(&String::from_utf8_lossy(&self.cab), 30)
            ));
        }
        if !self.le(1)? {
            return Err(self.truncado("the sequence line"));
        }
        if !self.le(2)? {
            return Err(self.truncado("the `+` separator line"));
        }
        if self.sep.first() != Some(&b'+') {
            return Err(format!(
                "{}, line {}: expected the `+` separator of record {}, but the line is `{}`. Every record must have exactly 4 lines: a sequence wrapped over several lines has to be joined into one first.",
                self.caminho,
                self.linha,
                self.registros + 1,
                mostra(&String::from_utf8_lossy(&self.sep), 30)
            ));
        }
        if !self.le(3)? && !self.seq.is_empty() {
            return Err(self.truncado("the quality line"));
        }
        if self.qual.len() != self.seq.len() {
            return Err(format!(
                "{}, line {}: record {} has {} bases but {} quality values. The file was probably cut or edited: fetch it again.",
                self.caminho,
                self.linha,
                self.registros + 1,
                self.seq.len(),
                self.qual.len()
            ));
        }
        self.registros += 1;
        Ok(true)
    }
}

/// A batch of read sequences, one after the other.
#[derive(Default)]
struct Lote {
    dados: Vec<u8>,
    fins: Vec<usize>,
}

impl Lote {
    /// Reads up to `max` records. `true` when the file ended.
    fn enche(&mut self, f: &mut Fastq, max: usize) -> Result<bool, String> {
        self.dados.clear();
        self.fins.clear();
        while self.fins.len() < max {
            if !f.proximo()? {
                return Ok(true);
            }
            self.dados.extend_from_slice(&f.seq);
            self.fins.push(self.dados.len());
        }
        Ok(false)
    }

    fn leitura(&self, k: usize) -> &[u8] {
        let ini = if k == 0 { 0 } else { self.fins[k - 1] };
        &self.dados[ini..self.fins[k]]
    }

    fn len(&self) -> usize {
        self.fins.len()
    }
}

/// MAGeCK's choice of offsets from the probe: the most frequent first, and
/// more while each one holds at least a third of the previous one's reads,
/// at least 5% of the matched reads, and the ones taken hold no more than
/// 75% together.
fn escolhe_offsets(contagem: &[(usize, u64)], total: u64) -> Vec<usize> {
    let mut v: Vec<(usize, u64)> = contagem.to_vec();
    multi_merge_sort_stable_by(&mut v, |a: &(usize, u64), b: &(usize, u64)| b.1.cmp(&a.1));
    let mut escolhidos = Vec::new();
    let mut soma = 0.0;
    let mut anterior = 0.01;
    for &(o, n) in &v {
        let f = n as f64 / total as f64;
        if anterior > f * 3.0 {
            break;
        }
        escolhidos.push(o);
        soma += f;
        if soma > 0.75 {
            break;
        }
        if f < 0.05 {
            break;
        }
        anterior = f;
    }
    escolhidos
}

// =====================================================================
//  PART 3: COUNTING
// =====================================================================

/// How `crispr_count` finds the guide in each read. The default is MAGeCK's:
/// offset found from the reads, library as given, exact match.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CountOptions {
    /// Where the guide starts in the read (0-based; MAGeCK's `--trim-5`).
    /// Several offsets are tried in order (staggered primers). `None`: found
    /// in each file from its first 100,000 reads, as MAGeCK's `AUTO` does.
    pub offsets: Option<Vec<usize>>,
    /// Reverse-complement the library sequences (MAGeCK's
    /// `--reverse-complement`), for reads that carry the guide on the other
    /// strand.
    pub reverse_complement: bool,
    /// 0: exact match only (MAGeCK). 1: a read with no exact match also
    /// counts for a guide one mismatch away (an N is a mismatch), when that
    /// guide is the only one that close.
    pub mismatches: u8,
}

/// One FASTQ file of a sample.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FileCount {
    pub path: String,
    pub reads: u64,
    pub mapped: u64,
    /// The offsets used, in the order they were tried.
    pub offsets: Vec<usize>,
    /// The offsets were found from the first reads (not given).
    pub detected: bool,
    /// The probe: (offset, reads whose first exact guide starts there), most
    /// frequent first. Empty when the offsets were given.
    pub offset_probe: Vec<(usize, u64)>,
    /// Reads looked at by the probe.
    pub probe_reads: u64,
}

/// One sample: its files and its totals.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SampleCount {
    pub label: String,
    pub files: Vec<FileCount>,
    pub reads: u64,
    /// Reads assigned to a guide (exact + one mismatch).
    pub mapped: u64,
    pub mapped_exact: u64,
    pub mapped_mismatch: u64,
    /// Reads one mismatch away from two or more guides: not counted.
    pub ambiguous: u64,
    /// 100 * mapped / reads.
    pub mapped_pct: f64,
    /// Library guides with no read.
    pub zero_guides: usize,
    /// MAGeCK's Gini index of ln(count + 1) over the library guides: 0 is an
    /// even library, 1 all reads on one guide.
    pub gini: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CountResult {
    pub library: Library,
    pub samples: Vec<SampleCount>,
    /// `counts[sample][guide]`, guides in library order.
    pub counts: Vec<Vec<u64>>,
    pub mismatches: u8,
    pub output: Option<String>,
    pub warnings: Vec<String>,
}

impl CountResult {
    pub fn json(&self) -> String {
        let amostras = json_lista(&self.samples, |s| {
            let arquivos = json_lista(&s.files, |f| {
                Objeto::novo()
                    .campo("path", json_texto(&f.path))
                    .campo("reads", f.reads.to_string())
                    .campo("mapped", f.mapped.to_string())
                    .campo("offsets", json_lista(&f.offsets, |o| o.to_string()))
                    .campo("detected", f.detected.to_string())
                    .campo("probe_reads", f.probe_reads.to_string())
                    .campo("offset_probe", json_lista(&f.offset_probe, |(o, n)| format!("[{o},{n}]")))
                    .fim()
            });
            Objeto::novo()
                .campo("label", json_texto(&s.label))
                .campo("files", arquivos)
                .campo("reads", s.reads.to_string())
                .campo("mapped", s.mapped.to_string())
                .campo("mapped_exact", s.mapped_exact.to_string())
                .campo("mapped_mismatch", s.mapped_mismatch.to_string())
                .campo("ambiguous", s.ambiguous.to_string())
                .campo("mapped_pct", json_num(s.mapped_pct))
                .campo("zero_guides", s.zero_guides.to_string())
                .campo("gini", json_num(s.gini))
                .fim()
        });
        let rotulos: Vec<String> = self.samples.iter().map(|s| s.label.clone()).collect();
        Objeto::novo()
            .campo("library", self.library.json())
            .campo("samples", amostras)
            .campo("labels", json_textos(&rotulos))
            .campo("sgrna", json_textos(&self.library.ids))
            .campo("gene", json_textos(&self.library.genes))
            .campo("counts", json_lista(&self.counts, |c| json_lista(c, |x| x.to_string())))
            .campo("mismatches", self.mismatches.to_string())
            .campo("output", self.output.as_deref().map(json_texto).unwrap_or_else(|| "null".to_string()))
            .campo("warnings", json_textos(&self.warnings))
            .fim()
    }

    pub fn resumo(&self) -> String {
        let mut o = format!(
            "Counted {} sample(s) against {} guides ({} genes) of {}\n",
            self.samples.len(),
            self.library.ids.len(),
            self.library.gene_count(),
            self.library.path
        );
        for s in &self.samples {
            let offs: Vec<String> = s
                .files
                .iter()
                .map(|f| {
                    let l: Vec<String> = f.offsets.iter().map(|x| x.to_string()).collect();
                    format!("{}{}", l.join(","), if f.detected { " (detected)" } else { "" })
                })
                .collect();
            o.push_str(&format!(
                "  {}: {} reads, {} matched a guide ({:.2}%), {} guide(s) with no read, Gini {:.4}; offset {}\n",
                s.label,
                s.reads,
                s.mapped,
                s.mapped_pct,
                s.zero_guides,
                s.gini,
                offs.join(" / ")
            ));
            if self.mismatches > 0 {
                o.push_str(&format!(
                    "    {} exact, {} with one mismatch, {} one mismatch from two or more guides (not counted)\n",
                    s.mapped_exact, s.mapped_mismatch, s.ambiguous
                ));
            }
        }
        if let Some(p) = &self.output {
            o.push_str(&format!("  count table written to {p}\n"));
        }
        for a in &self.warnings {
            o.push_str(&format!("  warning: {a}\n"));
        }
        o
    }

    /// The table MAGeCK writes as `*.countsummary.txt` (its first 8 columns).
    pub fn countsummary(&self) -> String {
        let mut o = String::from("File\tLabel\tReads\tMapped\tPercentage\tTotalsgRNAs\tZerocounts\tGiniIndex\n");
        for s in &self.samples {
            let arquivos: Vec<&str> = s.files.iter().map(|f| f.path.as_str()).collect();
            o.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                arquivos.join(","),
                s.label,
                s.reads,
                s.mapped,
                py_g(s.mapped as f64 / (s.reads as f64).max(1.0), 4),
                self.library.ids.len(),
                s.zero_guides,
                py_g(s.gini, 4)
            ));
        }
        o
    }
}

/// MAGeCK's Gini index: over ln(count + 1) of every library guide,
/// 1 - 2 (n - sum((i + 1) x_i) / sum(x)) / (n - 1), x sorted.
pub fn gini_index(counts: &[u64]) -> f64 {
    let mut xs: Vec<f64> = counts.iter().map(|&c| (c as f64 + 1.0).ln()).collect();
    ordena_f64(&mut xs);
    let n = xs.len();
    if n < 2 {
        return 0.0;
    }
    let gs: f64 = xs.iter().enumerate().map(|(i, &x)| (i as f64 + 1.0) * x).sum();
    let mut ys: f64 = xs.iter().sum();
    if ys == 0.0 {
        ys = 1.0;
    }
    1.0 - 2.0 * (n as f64 - gs / ys) / (n as f64 - 1.0)
}

/// Counts the reads of each guide in each sample (`mageck count`).
///
/// `fastq[i]` holds the file(s) of sample i -- several files are technical
/// replicates and are summed. `labels` names the samples (empty: `sample1`,
/// `sample2`...). With `output`, writes the count table (`sgRNA`, `Gene`, one
/// column per sample, guides in library order), the format `crispr_test` and
/// `mageck test` read.
pub fn crispr_count(
    library: &str,
    fastq: &[Vec<String>],
    labels: &[String],
    output: Option<&str>,
    opts: &CountOptions,
) -> Result<CountResult, String> {
    if fastq.is_empty() {
        return Err("No FASTQ files given: pass one file (or a list of files) per sample.".to_string());
    }
    if let Some((i, _)) = fastq.iter().enumerate().find(|(_, f)| f.is_empty()) {
        return Err(format!("Sample {} has no FASTQ file: give at least one file per sample.", i + 1));
    }
    if opts.mismatches > 1 {
        return Err(format!(
            "mismatches={} is not supported: use 0 (exact match, as MAGeCK) or 1 (one mismatch when only one guide is that close).",
            opts.mismatches
        ));
    }
    if let Some(o) = &opts.offsets {
        if o.is_empty() {
            return Err("offsets is an empty list: give at least one offset, or None to find it from the reads.".to_string());
        }
    }
    let rotulos: Vec<String> = if labels.is_empty() {
        (1..=fastq.len()).map(|i| format!("sample{i}")).collect()
    } else {
        if labels.len() != fastq.len() {
            return Err(format!(
                "{} label(s) for {} sample(s): give one label per sample (or none, for sample1, sample2, ...).",
                labels.len(),
                fastq.len()
            ));
        }
        let mut vistos = HashSet::new();
        for r in labels {
            if r.trim().is_empty() || r.contains(['\t', '\n', '\r', ',']) {
                return Err(format!("The sample label `{r}` is empty or has a tab, comma or line break: it becomes a column name of the count table."));
            }
            if !vistos.insert(r.as_str()) {
                return Err(format!("The sample label `{r}` is used twice: every sample needs its own label."));
            }
        }
        labels.to_vec()
    };
    if let Some(saida) = output {
        let mut entradas: Vec<&str> = vec![library];
        entradas.extend(fastq.iter().flatten().map(String::as_str));
        if let Some(e) = entradas.iter().find(|e| mesmo_arquivo(e, saida)) {
            return Err(format!("The output `{saida}` is the input `{e}`: writing would destroy it. Choose another output path."));
        }
    }
    let lib = read_library(library, opts.reverse_complement)?;
    let indice = Indice::novo(&lib, opts.mismatches);
    let mut res = CountResult { mismatches: opts.mismatches, ..Default::default() };

    for (arquivos, rotulo) in fastq.iter().zip(&rotulos) {
        let mut col = vec![0u64; lib.ids.len()];
        let mut sc = SampleCount { label: rotulo.clone(), ..Default::default() };
        for caminho in arquivos {
            let fc = conta_arquivo(caminho, &indice, opts, &mut col, &mut sc, &mut res.warnings)?;
            sc.files.push(fc);
        }
        sc.mapped = sc.mapped_exact + sc.mapped_mismatch;
        sc.mapped_pct = if sc.reads > 0 { 100.0 * sc.mapped as f64 / sc.reads as f64 } else { 0.0 };
        sc.zero_guides = col.iter().filter(|&&c| c == 0).count();
        sc.gini = gini_index(&col);
        if sc.reads > 0 && sc.mapped_pct < 50.0 {
            res.warnings.push(format!(
                "sample {}: only {:.1}% of the reads matched a guide. Check the offset, try reverse_complement, or allow one mismatch.",
                sc.label, sc.mapped_pct
            ));
        }
        res.samples.push(sc);
        res.counts.push(col);
    }
    if lib.duplicate_sequences > 0 {
        res.warnings.push(format!(
            "{} guide(s) of the library repeat the sequence of another guide and were left out (the first guide with each sequence gets the reads).",
            lib.duplicate_sequences
        ));
    }
    res.library = lib;

    if let Some(saida) = output {
        let mut g = Saida { caminho: saida.to_string(), armada: false };
        let mut w = Escritor::cria(saida)?;
        g.armada = true;
        let erro = |e: std::io::Error| format!("Error writing {saida}: {e}");
        let mut cab = String::from("sgRNA\tGene");
        for s in &res.samples {
            cab.push('\t');
            cab.push_str(&s.label);
        }
        w.linha(&cab).map_err(|e| format!("Error writing {saida}: {e}"))?;
        let mut linha = String::new();
        for (k, id) in res.library.ids.iter().enumerate() {
            linha.clear();
            linha.push_str(id);
            linha.push('\t');
            linha.push_str(&res.library.genes[k]);
            for c in &res.counts {
                linha.push('\t');
                linha.push_str(&c[k].to_string());
            }
            linha.push('\n');
            w.write_all(linha.as_bytes()).map_err(erro)?;
        }
        w.termina().map_err(|e| format!("Error writing {saida}: {e}"))?;
        g.armada = false;
        res.output = Some(saida.to_string());
    }
    Ok(res)
}

fn conta_arquivo(
    caminho: &str,
    indice: &Indice,
    opts: &CountOptions,
    col: &mut [u64],
    sc: &mut SampleCount,
    avisos: &mut Vec<String>,
) -> Result<FileCount, String> {
    let mut f = Fastq::abre(caminho)?;
    let mut fc = FileCount { path: caminho.to_string(), ..Default::default() };
    let mut lote = Lote::default();
    let mut fim = lote.enche(&mut f, OFFSET_PROBE_READS)?;
    let offsets: Vec<usize> = match &opts.offsets {
        Some(o) => o.clone(),
        None => {
            let n = lote.len();
            if n == 0 {
                return Err(format!("{caminho}: the file has no reads, so the guide offset cannot be found. Check the file, or give the offset."));
            }
            let achados: Vec<Option<usize>> =
                (0..n).into_par_iter().map(|k| indice.primeiro_offset(lote.leitura(k))).collect();
            // Offsets in order of first appearance, as MAGeCK's dict keeps them.
            let mut ordem: Vec<(usize, u64)> = Vec::new();
            let mut onde: HashMap<usize, usize> = HashMap::new();
            for o in achados.into_iter().flatten() {
                let i = *onde.entry(o).or_insert_with(|| {
                    ordem.push((o, 0));
                    ordem.len() - 1
                });
                ordem[i].1 += 1;
            }
            let total: u64 = ordem.iter().map(|x| x.1).sum();
            let frac = total as f64 / n as f64;
            if frac < 0.001 {
                return Err(format!(
                    "{caminho}: cannot find where the guide is in the reads -- only {total} of the first {n} reads ({:.3}%) contain a library guide at any position. Check that the FASTQ belongs to this library; if the reads carry the guide reverse-complemented, use reverse_complement=True; or give the offset.",
                    100.0 * frac
                ));
            }
            if frac < 0.5 {
                avisos.push(format!(
                    "{caminho}: only {:.1}% of the first {n} reads contain a library guide; the offset was found from those.",
                    100.0 * frac
                ));
            }
            let escolhidos = escolhe_offsets(&ordem, total);
            multi_merge_sort_stable_by(&mut ordem, |a: &(usize, u64), b: &(usize, u64)| b.1.cmp(&a.1));
            fc.offset_probe = ordem;
            fc.probe_reads = n as u64;
            fc.detected = true;
            escolhidos
        }
    };
    let troca = opts.mismatches > 0;
    loop {
        let acertos: Vec<Acerto> =
            (0..lote.len()).into_par_iter().map(|k| indice.casa(lote.leitura(k), &offsets, troca)).collect();
        fc.reads += lote.len() as u64;
        for a in acertos {
            match a {
                Acerto::Exato(g) => {
                    col[g as usize] += 1;
                    sc.mapped_exact += 1;
                    fc.mapped += 1;
                }
                Acerto::Troca(g) => {
                    col[g as usize] += 1;
                    sc.mapped_mismatch += 1;
                    fc.mapped += 1;
                }
                Acerto::Ambiguo => sc.ambiguous += 1,
                Acerto::Nada => {}
            }
        }
        if fim {
            break;
        }
        fim = lote.enche(&mut f, LOTE)?;
    }
    if fc.reads == 0 {
        avisos.push(format!("{caminho}: the file has no reads."));
    }
    sc.reads += fc.reads;
    fc.offsets = offsets;
    Ok(fc)
}

// =====================================================================
//  PART 4: THE COUNT TABLE
// =====================================================================

/// A count table: `sgRNA`, `Gene`, then one column of counts per sample
/// (MAGeCK's format; tab or comma).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CountTable {
    pub path: String,
    pub samples: Vec<String>,
    pub ids: Vec<String>,
    pub genes: Vec<String>,
    /// `counts[sample][guide]`.
    pub counts: Vec<Vec<f64>>,
    /// Lines whose guide id was already seen: skipped, as in MAGeCK.
    pub duplicate_ids: usize,
    /// The first line held the sample labels.
    pub header: bool,
}

impl CountTable {
    pub fn json(&self) -> String {
        Objeto::novo()
            .campo("path", json_texto(&self.path))
            .campo("samples", json_textos(&self.samples))
            .campo("guides", self.ids.len().to_string())
            .campo("genes", self.genes.iter().collect::<HashSet<_>>().len().to_string())
            .campo("duplicate_ids", self.duplicate_ids.to_string())
            .campo("header", self.header.to_string())
            .fim()
    }

    pub fn resumo(&self) -> String {
        format!(
            "{}: {} guides, {} genes, {} sample(s): {}{}\n",
            self.path,
            self.ids.len(),
            self.genes.iter().collect::<HashSet<_>>().len(),
            self.samples.len(),
            self.samples.join(", "),
            if self.duplicate_ids > 0 { format!("; {} repeated guide id(s) skipped", self.duplicate_ids) } else { String::new() }
        )
    }
}

fn numero(s: &str) -> Option<f64> {
    s.parse::<f64>().ok().filter(|x| x.is_finite())
}

/// Reads a count table: guide id, gene, one column per sample; tab- or
/// comma-separated (detected); the first line is the header when its count
/// columns are not numbers (otherwise the samples are `sample1`...).
pub fn read_count_table(caminho: &str) -> Result<CountTable, String> {
    let mut l = Linhas::abre(caminho)?;
    let mut t = CountTable { path: caminho.to_string(), ..Default::default() };
    let mut sep: Option<char> = None;
    let mut largura = 0usize;
    let mut ids: HashSet<String> = HashSet::new();
    let mut linhas: Vec<Vec<f64>> = Vec::new();
    while l.avanca()? {
        let txt = l.texto();
        if txt.trim().is_empty() {
            continue;
        }
        let s = match sep {
            Some(s) => s,
            None => match separador(txt) {
                Some(s) => {
                    sep = Some(s);
                    s
                }
                None => {
                    return Err(format!(
                        "{caminho}, line {}: no tab or comma in `{}`. A count table has the columns sgRNA, Gene and one column of counts per sample, separated by tabs or commas.",
                        l.linha,
                        mostra(txt, 60)
                    ))
                }
            },
        };
        let c = campos(txt, s);
        if largura == 0 {
            if c.len() < 3 {
                return Err(format!(
                    "{caminho}, line {}: {} column(s); a count table needs at least 3 -- sgRNA, Gene, and one column of counts per sample.",
                    l.linha,
                    c.len()
                ));
            }
            largura = c.len();
            if !c[2..].iter().all(|x| numero(x).is_some()) {
                t.header = true;
                let mut vistos = HashSet::new();
                for r in &c[2..] {
                    if r.is_empty() {
                        return Err(format!("{caminho}, line 1: a sample column has no name in the header."));
                    }
                    if !vistos.insert(*r) {
                        return Err(format!(
                            "{caminho}, line {}: the sample `{r}` names two columns. Rename one: samples are chosen by name.",
                            l.linha
                        ));
                    }
                }
                t.samples = c[2..].iter().map(|x| x.to_string()).collect();
                continue;
            }
            t.samples = (1..=largura - 2).map(|i| format!("sample{i}")).collect();
        }
        if c.len() != largura {
            return Err(format!(
                "{caminho}, line {}: {} columns, but the first line has {largura}. Every line needs the guide id, the gene and one count per sample.",
                l.linha,
                c.len()
            ));
        }
        if c[0].is_empty() {
            return Err(format!("{caminho}, line {}: the guide id (first column) is empty.", l.linha));
        }
        if c[1].is_empty() {
            return Err(format!(
                "{caminho}, line {}: guide `{}` has no gene (second column). Give every guide a gene; MAGeCK refuses empty genes too.",
                l.linha, c[0]
            ));
        }
        let mut v = Vec::with_capacity(largura - 2);
        for (j, x) in c[2..].iter().enumerate() {
            match numero(x) {
                Some(n) if n >= 0.0 => v.push(n),
                _ => {
                    return Err(format!(
                        "{caminho}, line {}: `{}` in column `{}` (guide `{}`) is not a read count. Counts must be numbers >= 0.",
                        l.linha,
                        mostra(x, 30),
                        t.samples[j],
                        c[0]
                    ))
                }
            }
        }
        if !ids.insert(c[0].to_string()) {
            t.duplicate_ids += 1;
            continue;
        }
        t.ids.push(c[0].to_string());
        t.genes.push(c[1].to_string());
        linhas.push(v);
    }
    if t.ids.is_empty() {
        return Err(format!("{caminho}: no guides found. A count table has one line per guide: sgRNA, Gene, one count per sample."));
    }
    t.counts = (0..t.samples.len()).map(|s| linhas.iter().map(|v| v[s]).collect()).collect();
    Ok(t)
}

// =====================================================================
//  PART 5: THE NORMAL DISTRIBUTION
// =====================================================================

const SQRT_PI: f64 = 1.772_453_850_905_516;

/// The complementary error function, relative error ~1e-15: a series for
/// x < 1, a continued fraction (modified Lentz) above; exp(-x^2) is taken in
/// two parts so the rounding of x^2 does not spoil the far tail.
pub fn erfc(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 {
        return 2.0 - erfc(-x);
    }
    if x > 27.3 {
        return 0.0;
    }
    // exp(-x^2) = exp(-z^2) exp((z - x)(z + x)), z = x with the low 32 bits
    // of the mantissa cleared, so z^2 is exact.
    let z = f64::from_bits(x.to_bits() & 0xFFFF_FFFF_0000_0000);
    let e = (-z * z).exp() * ((z - x) * (z + x)).exp();
    if x < 1.0 {
        // erf(x) = 2x/sqrt(pi) exp(-x^2) sum (2x^2)^n / (1 3 5 ... (2n+1))
        let x2 = x * x;
        let (mut termo, mut soma, mut n) = (1.0f64, 1.0f64, 0u32);
        loop {
            n += 1;
            termo *= 2.0 * x2 / (2 * n + 1) as f64;
            soma += termo;
            if termo < soma * 1e-17 || n > 200 {
                break;
            }
        }
        return 1.0 - 2.0 * x / SQRT_PI * e * soma;
    }
    // erfc(x) = exp(-x^2)/sqrt(pi) / f, f = x + (1/2)/(x + 1/(x + (3/2)/(x + ...)))
    let minusculo = 1e-300;
    let mut f = x;
    let mut c = f;
    let mut d = 0.0f64;
    for n in 1..5000 {
        let a = n as f64 / 2.0;
        d = x + a * d;
        if d.abs() < minusculo {
            d = minusculo;
        }
        c = x + a / c;
        if c.abs() < minusculo {
            c = minusculo;
        }
        d = 1.0 / d;
        let delta = c * d;
        f *= delta;
        if (delta - 1.0).abs() < 1e-16 {
            break;
        }
    }
    e / SQRT_PI / f
}

/// P(Z <= x), as MAGeCK writes it: erfc(-x / sqrt 2) / 2.
pub fn norm_cdf(x: f64) -> f64 {
    erfc(-x / std::f64::consts::SQRT_2) / 2.0
}

/// P(Z > x): erfc(x / sqrt 2) / 2.
pub fn norm_sf(x: f64) -> f64 {
    erfc(x / std::f64::consts::SQRT_2) / 2.0
}

/// Acklam's rational approximation of the normal quantile (the same
/// coefficients as `normalInv` in MAGeCK's RRA), for 0 < p < 1.
fn acklam(p: f64) -> f64 {
    let a = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.357751867269, -30.66479806614716, 2.506628277459239];
    let b = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572];
    let c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783];
    let d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416];
    let (pb, pa) = (0.02425, 1.0 - 0.02425);
    if p < pb {
        let q = (-2.0 * p.ln()).sqrt();
        (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    } else if p <= pa {
        let q = p - 0.5;
        let r = q * q;
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
            / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    }
}

/// The normal quantile of a lower-tail probability p <= 0.5 (x <= 0):
/// Acklam, then two Halley steps.
fn qnorm_baixo(p: f64) -> f64 {
    let mut x = acklam(p);
    for _ in 0..2 {
        let e = norm_cdf(x) - p;
        let u = e * (2.0 * std::f64::consts::PI).sqrt() * (x * x / 2.0).exp();
        x -= u / (1.0 + x * u / 2.0);
    }
    x
}

/// The normal quantile of exp(ln_p), precise at both ends.
pub fn qnorm_log(ln_p: f64) -> f64 {
    if ln_p < 0.5f64.ln() {
        qnorm_baixo(ln_p.exp())
    } else {
        -qnorm_baixo(-ln_p.exp_m1())
    }
}

/// Rows of MAGeCK's `qnorm_table.txt`: ln p from -24.000 to -0.001.
const QN: usize = 24000;

/// The point i of MAGeCK's `qnorm_table.txt`: ln p = -24 + i/1000, as the
/// file's decimal text reads.
fn grade_logp(i: usize) -> f64 {
    (i as f64 - 24000.0) / 1000.0
}

fn grade_q() -> &'static [f64] {
    static T: OnceLock<Vec<f64>> = OnceLock::new();
    T.get_or_init(|| (0..QN).into_par_iter().map(|i| qnorm_log(grade_logp(i))).collect())
}

/// MAGeCK's `QNormConverter.get_qnorm(ln_p, islog=True)`: linear
/// interpolation in the table; `None` outside ln p in [-24, -0.001].
fn qnorm_mageck(lp: f64) -> Option<f64> {
    let t = grade_q();
    let (qmin, qmax) = (grade_logp(0), grade_logp(QN - 1));
    let passo = grade_logp(1) - grade_logp(0);
    let i = (((lp - qmin) / passo).floor() as i64).clamp(0, QN as i64 - 2) as usize;
    let frac = ((lp - grade_logp(i)) / passo).clamp(0.0, 1.0);
    let v = t[i] * (1.0 - frac) + t[i + 1] * frac;
    if lp >= qmin && lp <= qmax {
        Some(v)
    } else {
        None
    }
}

// =====================================================================
//  PART 6: ROBUST RANK AGGREGATION (MAGeCK's RRA program, ported)
// =====================================================================

/// `LogGamma` of RRA's math_api.cpp.
fn log_gamma(x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    let mut y = x;
    let f;
    if x < 7.0 {
        let mut ff = 1.0;
        let mut z = y;
        while z < 7.0 {
            ff *= z;
            z += 1.0;
        }
        y = z;
        f = -ff.ln();
    } else {
        f = 0.0;
    }
    let z = 1.0 / y / y;
    f + (y - 0.5) * y.ln() - y
        + 0.918938533204673
        + (((-0.000595238095238 * z + 0.000793650793651) * z - 0.002777777777778) * z + 0.083333333333333) / y
}

/// The regularized incomplete beta function I_x(p, q), algorithm AS 63 as
/// RRA has it (`betain`); `beta` = ln B(p, q).
fn betain(x: f64, p: f64, q: f64, beta: f64) -> f64 {
    let acu = 0.1e-14;
    if p <= 0.0 || q <= 0.0 || !(0.0..=1.0).contains(&x) {
        return x;
    }
    if x == 0.0 || x == 1.0 {
        return x;
    }
    let mut psq = p + q;
    let mut cx = 1.0 - x;
    let (xx, pp, qq, indx);
    if p < psq * x {
        xx = cx;
        cx = x;
        pp = q;
        qq = p;
        indx = true;
    } else {
        xx = x;
        pp = p;
        qq = q;
        indx = false;
    }
    let mut term = 1.0;
    let mut ai = 1.0;
    let mut value = 1.0;
    let mut ns = (qq + cx * psq) as i64;
    let mut rx = xx / cx;
    let mut temp = qq - ai;
    if ns == 0 {
        rx = xx;
    }
    loop {
        term = term * temp * rx / (pp + ai);
        value += term;
        temp = term.abs();
        if temp <= acu && temp <= acu * value {
            value = value * (pp * xx.ln() + (qq - 1.0) * cx.ln() - beta).exp() / pp;
            if indx {
                value = 1.0 - value;
            }
            break;
        }
        ai += 1.0;
        ns -= 1;
        if 0 <= ns {
            temp = qq - ai;
            if ns == 0 {
                rx = xx;
            }
        } else {
            temp = psq;
            psq += 1.0;
        }
    }
    value
}

/// ln B(a, b) as RRA computes it.
fn log_beta(a: f64, b: f64) -> f64 {
    log_gamma(a) + log_gamma(b) - log_gamma(a + b)
}

/// The beta CDF, as RRA computes it (`BetaNoncentralCdf` with lambda 0).
pub fn beta_cdf(a: f64, b: f64, x: f64) -> f64 {
    betain(x, a, b, log_beta(a, b))
}

/// RRA's `ComputeLoValue`: `pct` sorted ascending. Returns (lo-value, good
/// guides = how many percentiles fall in [min_p, max_p] before the first one
/// outside, the first one excepted).
pub fn lo_value(pct: &[f64], max_p: f64, min_p: f64) -> (f64, usize) {
    let n = pct.len();
    if n == 0 {
        return (1.0, 0);
    }
    let mut lo = 1.0f64;
    let mut bons = 0;
    for (i, &x) in pct.iter().enumerate() {
        if x > max_p || x < min_p {
            if i > 0 {
                break;
            }
        } else {
            bons += 1;
        }
        let f = beta_cdf((i + 1) as f64, (n - i) as f64, x);
        if f < lo {
            lo = f;
        }
    }
    (lo, bons)
}

/// RRA's `bTreeSearchingF` on an ascending array: the index of `v` if it is
/// outside, else of the nearer of its two neighbours.
fn busca_btree(v: f64, a: &[f64]) -> usize {
    let (mut lo, mut hi) = (0usize, a.len() - 1);
    loop {
        if v <= a[lo] {
            return lo;
        }
        if v >= a[hi] {
            return hi;
        }
        if hi - lo <= 1 {
            return if (a[hi] - v).abs() > (a[lo] - v).abs() { lo } else { hi };
        }
        let mid = (lo + hi) / 2;
        if v >= a[mid] {
            lo = mid;
        } else {
            hi = mid;
        }
    }
}

/// RRA's mid-rank of `v` in the ascending array `a`, over 2n: the percentile
/// of an item, and the permutation p-value of a lo-value.
fn meio_posto(v: f64, a: &[f64]) -> f64 {
    (busca_btree(v - 0.000000001, a) + busca_btree(v + 0.000000001, a) + 1) as f64
}

const PM_M: u64 = 2_147_483_647;
const PM_A: u64 = 48_271;

/// One step of RRA's random generator (Park & Miller, rngs.c): seed * 48271
/// mod 2^31 - 1, and seed / (2^31 - 1). The modulo of a Mersenne prime is
/// folded instead of divided.
#[inline]
fn pm_proximo(s: &mut u64) -> f64 {
    let x = *s * PM_A;
    let mut r = (x & PM_M) + (x >> 31);
    if r >= PM_M {
        r -= PM_M;
    }
    *s = r;
    r as f64 / PM_M as f64
}

/// The generator `n` steps ahead: seed * 48271^n mod 2^31 - 1.
fn pm_salta(s: u64, mut n: u64) -> u64 {
    let (mut r, mut b) = (s % PM_M, PM_A);
    while n > 0 {
        if n & 1 == 1 {
            r = r * b % PM_M;
        }
        b = b * b % PM_M;
        n >>= 1;
    }
    r
}

/// `total` lo-values of random genes of `s` guides, sorted -- RRA's
/// `ComputePermutationbyNumSG` for one size. The uniforms come from the
/// generator in RRA's order; blocks are computed in parallel by jumping the
/// generator ahead, which gives exactly the same numbers.
///
/// Three shortcuts, none of which changes a bit of the result: only the
/// percentiles at or below alpha (and the smallest one) reach the lo-value,
/// so only those are sorted; ln B(i + 1, s - i) is computed once per size;
/// and since the lo-value keeps only the minimum, a beta CDF that cannot go
/// below the current minimum is not computed. That is decided with
/// thresholds: F_i is increasing in x, so x >= q_i(L), where F_i(q_i(L)) is
/// at least L (1 + 1e-6), means F_i(x) > L > the current minimum.
fn lo_aleatorios(semente: &mut u64, s: usize, total: usize, alpha: f64, nulo: Option<&[f64]>) -> Vec<f64> {
    const BLOCO: usize = 4096;
    const NIVEIS: [f64; 12] = [1e-4, 1e-3, 5e-3, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 0.9];
    let inicio = *semente;
    let blocos = total.div_ceil(BLOCO);
    let lb: Vec<f64> = (0..s).map(|i| log_beta((i + 1) as f64, (s - i) as f64)).collect();
    let limiar: Vec<[f64; 12]> = (0..s)
        .into_par_iter()
        .map(|i| {
            let (a, b) = ((i + 1) as f64, (s - i) as f64);
            let mut q = [1.0f64; 12];
            for (k, &nivel) in NIVEIS.iter().enumerate() {
                let alvo = nivel * (1.0 + 1e-6);
                let (mut baixo, mut alto) = (0.0f64, 1.0f64);
                for _ in 0..64 {
                    let m = 0.5 * (baixo + alto);
                    if betain(m, a, b, lb[i]) >= alvo {
                        alto = m;
                    } else {
                        baixo = m;
                    }
                }
                q[k] = alto;
            }
            q
        })
        .collect();
    let mut v: Vec<f64> = (0..blocos)
        .into_par_iter()
        .flat_map_iter(|b| {
            let mut st = pm_salta(inicio, (b * BLOCO * s) as u64);
            let n = BLOCO.min(total - b * BLOCO);
            let mut abaixo: Vec<f64> = Vec::with_capacity(s);
            let mut out = Vec::with_capacity(n);
            for _ in 0..n {
                abaixo.clear();
                let mut menor = f64::INFINITY;
                for _ in 0..s {
                    let u = pm_proximo(&mut st);
                    let x = match nulo {
                        Some(c) => c[((c.len() as f64 * u) as usize).min(c.len() - 1)],
                        None => u,
                    };
                    if x <= alpha {
                        abaixo.push(x);
                    }
                    if x < menor {
                        menor = x;
                    }
                }
                let lo = if abaixo.is_empty() {
                    betain(menor, 1.0, s as f64, lb[0]).min(1.0)
                } else {
                    ordena_f64(&mut abaixo);
                    let mut lo = 1.0f64;
                    for (i, &x) in abaixo.iter().enumerate() {
                        let k = NIVEIS.partition_point(|&v| v <= lo);
                        if k < NIVEIS.len() && x >= limiar[i][k] {
                            continue;
                        }
                        let f = betain(x, (i + 1) as f64, (s - i) as f64, lb[i]);
                        if f < lo {
                            lo = f;
                        }
                    }
                    lo
                };
                out.push(lo);
            }
            out
        })
        .collect();
    *semente = pm_salta(inicio, (total * s) as u64);
    ordena_f64(&mut v);
    v
}

/// Benjamini-Hochberg as RRA computes it (gene FDR).
fn bh_rra(p: &[f64]) -> Vec<f64> {
    let n = p.len();
    if n == 0 {
        return Vec::new();
    }
    let idx = ordem_por(p, true);
    let mut q: Vec<f64> = idx.iter().enumerate().map(|(i, &k)| p[k] / (i as f64 + 1.0) * n as f64).collect();
    if q[n - 1] > 1.0 {
        q[n - 1] = 1.0;
    }
    for i in (0..n - 1).rev() {
        if q[i] > q[i + 1] {
            q[i] = q[i + 1];
        }
    }
    let mut out = vec![0.0; n];
    for (i, &k) in idx.iter().enumerate() {
        out[k] = q[i];
    }
    out
}

/// Benjamini-Hochberg as MAGeCK's Python computes it (`pFDR`, guide FDR).
pub fn bh_fdr(p: &[f64]) -> Vec<f64> {
    let n = p.len();
    if n == 0 {
        return Vec::new();
    }
    let mut idx = ordem_por(p, true);
    idx.reverse();
    let mut v: Vec<f64> = idx.iter().enumerate().map(|(i, &k)| (n as f64 / (n - i) as f64) * p[k]).collect();
    for i in 0..n - 1 {
        if v[i] < v[i + 1] {
            v[i + 1] = v[i];
        }
    }
    let mut out = vec![0.0; n];
    for (i, &k) in idx.iter().enumerate() {
        out[k] = v[i];
    }
    out
}

/// One gene of an RRA run.
#[derive(Clone, Debug)]
struct GeneRra {
    nome: String,
    itens: usize,
    lo: f64,
    p: f64,
    fdr: f64,
    bons: usize,
}

struct Rra {
    /// In RRA's output order (lo-value ascending; ties in order of first
    /// appearance), without the genes made only of control guides.
    genes: Vec<GeneRra>,
    passes: usize,
    aleatorios: usize,
}

/// MAGeCK's RRA program on one ranked list. `ordem`: the guides sorted by
/// `valor` (ascending: most selected first); only `valido` guides enter.
#[allow(clippy::too_many_arguments)]
fn rra(
    valor: &[f64],
    ordem: &[usize],
    valido: &[bool],
    genes: &[String],
    ids: &[String],
    controle: &[bool],
    alpha: f64,
    o: &TestOptions,
) -> Result<Rra, String> {
    let itens: Vec<usize> = ordem.iter().copied().filter(|&g| valido[g]).collect();
    let n = itens.len();
    let mut lista: Vec<f64> = itens.iter().map(|&g| valor[g]).collect();
    ordena_f64(&mut lista);
    let usa_controle = controle.iter().any(|&c| c);

    let mut onde: HashMap<&str, usize> = HashMap::new();
    let mut grupos: Vec<(&str, Vec<f64>, usize)> = Vec::new();
    let mut pct_controle: Vec<(&str, f64)> = Vec::new();
    for &g in &itens {
        let pct = meio_posto(valor[g], &lista) / (n * 2) as f64;
        let i = *onde.entry(genes[g].as_str()).or_insert_with(|| {
            grupos.push((genes[g].as_str(), Vec::new(), 0));
            grupos.len() - 1
        });
        grupos[i].1.push(pct);
        if usa_controle && controle[g] {
            grupos[i].2 += 1;
            pct_controle.push((ids[g].as_str(), pct));
        }
    }
    let ng = grupos.len();
    let pular: HashSet<&str> = o.skip_genes.iter().map(String::as_str).collect();
    let mut res: Vec<GeneRra> = Vec::with_capacity(ng);
    let mut pulado: Vec<bool> = Vec::with_capacity(ng);
    let mut so_controle: Vec<bool> = Vec::with_capacity(ng);
    for (nome, pct, nc) in grupos.iter_mut() {
        ordena_f64(pct);
        let (lo, bons) = lo_value(pct, alpha, -1.0);
        let sc = usa_controle && *nc == pct.len();
        so_controle.push(sc);
        pulado.push(sc || pular.contains(nome));
        res.push(GeneRra { nome: nome.to_string(), itens: pct.len(), lo, p: 1.0, fdr: 1.0, bons });
    }

    // Permutation, gene size by gene size.
    let passes = match o.permutations {
        Some(p) => p,
        None if 100 * ng < 100_000 => 100_000 / ng + 1,
        None => 100,
    };
    let total = ng * (passes + 1);
    let nulo: Option<Vec<f64>> = if usa_controle {
        // RRA's control list is a std::map: names in byte order.
        multi_merge_sort_stable_by(&mut pct_controle, |a: &(&str, f64), b: &(&str, f64)| a.0.cmp(b.0));
        if pct_controle.is_empty() {
            return Err("None of the control sgRNAs is left after the zero-count filter, so they cannot give the null distribution. Use remove_zero=\"none\" or check the control list.".to_string());
        }
        Some(pct_controle.iter().map(|x| x.1).collect())
    } else {
        None
    };
    let mut tamanhos: Vec<usize> = (0..ng).filter(|&i| !pulado[i]).map(|i| res[i].itens).collect::<HashSet<_>>().into_iter().collect();
    multi_merge_sort_stable_by(&mut tamanhos, |a: &usize, b: &usize| a.cmp(b));
    let mut semente = (o.seed as u64) % PM_M;
    let mut ultimo: Option<Vec<f64>> = None;
    for s in tamanhos {
        if s <= o.max_sgrna_per_gene_permutation || ultimo.is_none() {
            ultimo = Some(lo_aleatorios(&mut semente, s, total, alpha, nulo.as_deref()));
        }
        let r = ultimo.as_ref().unwrap();
        for i in 0..ng {
            if res[i].itens == s && !pulado[i] {
                res[i].p = meio_posto(res[i].lo, r) / 2.0 / r.len() as f64;
            }
        }
    }
    let ps: Vec<f64> = res.iter().map(|g| g.p).collect();
    for (g, q) in res.iter_mut().zip(bh_rra(&ps)) {
        g.fdr = q;
    }
    let los: Vec<f64> = res.iter().map(|g| g.lo).collect();
    let ordem_lo = ordem_por(&los, true);
    let genes_saida: Vec<GeneRra> = ordem_lo.into_iter().filter(|&i| !so_controle[i]).map(|i| res[i].clone()).collect();
    Ok(Rra { genes: genes_saida, passes, aleatorios: total })
}

// =====================================================================
//  PART 7: THE SCREEN TEST
// =====================================================================

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NormMethod {
    /// No normalization.
    None,
    /// Median ratio (MAGeCK's default).
    Median,
    /// Total read count.
    Total,
    /// Median ratio over the control guides only.
    Control,
}

impl NormMethod {
    pub fn from_text(s: &str) -> Result<NormMethod, String> {
        match s.to_ascii_lowercase().as_str() {
            "none" => Ok(NormMethod::None),
            "median" => Ok(NormMethod::Median),
            "total" => Ok(NormMethod::Total),
            "control" => Ok(NormMethod::Control),
            _ => Err(format!("Unknown normalization `{s}`: use \"median\" (default), \"total\", \"control\" or \"none\".")),
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            NormMethod::None => "none",
            NormMethod::Median => "median",
            NormMethod::Total => "total",
            NormMethod::Control => "control",
        }
    }
}

/// Which guides leave the test because they have no reads (MAGeCK's
/// `--remove-zero`): their median normalized count is not above the
/// threshold in...
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RemoveZero {
    /// ... nothing: every guide stays.
    None,
    /// ... the control group.
    Control,
    /// ... the treatment group.
    Treatment,
    /// ... both groups (default).
    Both,
    /// ... either group.
    Any,
}

impl RemoveZero {
    pub fn from_text(s: &str) -> Result<RemoveZero, String> {
        match s.to_ascii_lowercase().as_str() {
            "none" => Ok(RemoveZero::None),
            "control" => Ok(RemoveZero::Control),
            "treatment" => Ok(RemoveZero::Treatment),
            "both" => Ok(RemoveZero::Both),
            "any" => Ok(RemoveZero::Any),
            _ => Err(format!("Unknown remove_zero `{s}`: use \"both\" (default), \"none\", \"control\", \"treatment\" or \"any\".")),
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            RemoveZero::None => "none",
            RemoveZero::Control => "control",
            RemoveZero::Treatment => "treatment",
            RemoveZero::Both => "both",
            RemoveZero::Any => "any",
        }
    }
}

/// How the gene log2 fold change is taken from its guides (MAGeCK's
/// `--gene-lfc-method`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GeneLfc {
    /// Median of the guides (default).
    Median,
    Mean,
    /// Median of the guides ranked before the alpha cut-off.
    AlphaMedian,
    AlphaMean,
    /// The second strongest guide.
    SecondBest,
}

impl GeneLfc {
    pub fn from_text(s: &str) -> Result<GeneLfc, String> {
        match s.to_ascii_lowercase().as_str() {
            "median" => Ok(GeneLfc::Median),
            "mean" => Ok(GeneLfc::Mean),
            "alphamedian" => Ok(GeneLfc::AlphaMedian),
            "alphamean" => Ok(GeneLfc::AlphaMean),
            "secondbest" => Ok(GeneLfc::SecondBest),
            _ => Err(format!(
                "Unknown gene_lfc `{s}`: use \"median\" (default), \"mean\", \"alphamedian\", \"alphamean\" or \"secondbest\"."
            )),
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            GeneLfc::Median => "median",
            GeneLfc::Mean => "mean",
            GeneLfc::AlphaMedian => "alphamedian",
            GeneLfc::AlphaMean => "alphamean",
            GeneLfc::SecondBest => "secondbest",
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct TestOptions {
    pub norm_method: NormMethod,
    /// Negative-control guide ids: they give the null distribution of RRA,
    /// and the factors with `NormMethod::Control`. Genes made only of
    /// control guides are left out of the gene table.
    pub control_sgrnas: Vec<String>,
    /// Or genes whose guides are the negative controls.
    pub control_genes: Vec<String>,
    pub remove_zero: RemoveZero,
    pub remove_zero_threshold: f64,
    /// p-value threshold that sets the alpha of RRA (default 0.25).
    pub gene_test_fdr_threshold: f64,
    pub gene_lfc: GeneLfc,
    /// Rounds of permutation. `None`: 100, raised so that at least 100,000
    /// random genes are drawn (RRA's rule).
    pub permutations: Option<usize>,
    /// Gene sizes above this reuse the null of the largest size permuted
    /// (RRA's `--max-sgrnapergene-permutation`, 100).
    pub max_sgrna_per_gene_permutation: usize,
    /// Genes left out of the permutation (p = 1): by default `NA` and `na`.
    pub skip_genes: Vec<String>,
    /// Order of the gene table: by negative (default) or positive selection.
    pub sort_by_positive: bool,
    /// Seed of the permutation generator (RRA's is 123456).
    pub seed: u32,
}

impl Default for TestOptions {
    fn default() -> Self {
        TestOptions {
            norm_method: NormMethod::Median,
            control_sgrnas: Vec::new(),
            control_genes: Vec::new(),
            remove_zero: RemoveZero::Both,
            remove_zero_threshold: 0.0,
            gene_test_fdr_threshold: 0.25,
            gene_lfc: GeneLfc::Median,
            permutations: None,
            max_sgrna_per_gene_permutation: 100,
            skip_genes: vec!["NA".to_string(), "na".to_string()],
            sort_by_positive: false,
            seed: 123456,
        }
    }
}

/// One guide of the test: a line of `*.sgrna_summary.txt`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GuideStat {
    pub sgrna: String,
    pub gene: String,
    /// Normalized counts, control samples in the order given.
    pub control_counts: Vec<f64>,
    pub treatment_counts: Vec<f64>,
    /// Median of the control counts (zero raised to the smallest non-zero
    /// median, as in MAGeCK).
    pub control_mean: f64,
    /// Median of the treatment counts.
    pub treatment_mean: f64,
    /// log2(treatment_mean + 1) - log2(control_mean + 1).
    pub lfc: f64,
    pub control_var: f64,
    /// Variance from the mean-variance model.
    pub adj_var: f64,
    /// |guide score|: the RRA ranking value (MAGeCK's `score` column).
    pub score: f64,
    /// The signed score: negative when depleted.
    pub theta: f64,
    pub p_low: f64,
    pub p_high: f64,
    pub p_twosided: f64,
    /// BH over p_twosided of all the guides.
    pub fdr: f64,
    pub high_in_treatment: bool,
    /// Passed the zero-count filter and entered RRA.
    pub used: bool,
}

/// One gene: a line of `*.gene_summary.txt`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GeneStat {
    pub gene: String,
    /// Guides of the gene that entered RRA.
    pub guides: usize,
    /// RRA lo-value (MAGeCK's `neg|score`): smaller = more depleted.
    pub neg_score: f64,
    pub neg_p: f64,
    pub neg_fdr: f64,
    /// 1 = most depleted.
    pub neg_rank: usize,
    /// Guides below the alpha cut-off.
    pub neg_good: usize,
    pub neg_lfc: f64,
    pub pos_score: f64,
    pub pos_p: f64,
    pub pos_fdr: f64,
    pub pos_rank: usize,
    pub pos_good: usize,
    pub pos_lfc: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TestResult {
    pub count_table: String,
    pub treatment: Vec<String>,
    pub control: Vec<String>,
    pub norm_method: String,
    /// What was applied: differs from `norm_method` when median
    /// normalization fell back to total counts.
    pub norm_used: String,
    /// (sample, factor), control samples then treatment samples.
    pub size_factors: Vec<(String, f64)>,
    /// The variance model var = mean + mean^k * 2^b.
    pub model_k: f64,
    pub model_b: f64,
    /// Guides that entered the fit, and guides left out as high-count outliers.
    pub model_guides: usize,
    pub model_outliers: usize,
    pub guides_total: usize,
    pub guides_used: usize,
    pub genes_total: usize,
    /// Control guides found in the table.
    pub control_guides: usize,
    /// RRA alpha cut-offs.
    pub alpha_neg: f64,
    pub alpha_pos: f64,
    pub permutations: usize,
    /// Random genes drawn per gene size.
    pub random_genes: usize,
    /// In `*.sgrna_summary.txt` order: by |score|, strongest first.
    pub guides: Vec<GuideStat>,
    /// In `*.gene_summary.txt` order: by negative rank (or positive).
    pub genes: Vec<GeneStat>,
    /// Genes left out of the permutation (`TestOptions::skip_genes`): they
    /// stay in the tables, with p = 1, but not in the "top" lists.
    pub skipped_genes: Vec<String>,
    pub outputs: Vec<String>,
    pub warnings: Vec<String>,
}

impl TestResult {
    /// Genes with FDR below each threshold (0.01, 0.05, 0.25): (negative, positive).
    pub fn hits(&self) -> [(usize, usize); 3] {
        let mut h = [(0, 0); 3];
        for (k, t) in [0.01, 0.05, 0.25].iter().enumerate() {
            h[k] = (
                self.genes.iter().filter(|g| g.neg_fdr < *t).count(),
                self.genes.iter().filter(|g| g.pos_fdr < *t).count(),
            );
        }
        h
    }

    /// The `n` best genes in one direction, leaving out the skipped genes
    /// (`NA` by default), which have no p-value.
    pub fn top(&self, n: usize, positive: bool) -> Vec<&GeneStat> {
        let mut v = self.por_rank(positive);
        v.retain(|g| !self.skipped_genes.contains(&g.gene));
        v.truncate(n);
        v
    }

    /// Every gene, by rank in one direction.
    fn por_rank(&self, positive: bool) -> Vec<&GeneStat> {
        let mut v: Vec<&GeneStat> = self.genes.iter().collect();
        if positive {
            multi_merge_sort_stable_by(&mut v, |a: &&GeneStat, b: &&GeneStat| a.pos_rank.cmp(&b.pos_rank));
        } else {
            multi_merge_sort_stable_by(&mut v, |a: &&GeneStat, b: &&GeneStat| a.neg_rank.cmp(&b.neg_rank));
        }
        v
    }

    /// What a rank plot needs, for one direction: every gene from rank 1 on,
    /// with its RRA score, -log10 of it, LFC and FDR.
    fn grafico(&self, positivo: bool) -> String {
        let v = self.por_rank(positivo);
        let nomes: Vec<String> = v.iter().map(|g| g.gene.clone()).collect();
        let sel = |g: &GeneStat| if positivo { (g.pos_rank, g.pos_score, g.pos_lfc, g.pos_fdr) } else { (g.neg_rank, g.neg_score, g.neg_lfc, g.neg_fdr) };
        Objeto::novo()
            .campo("gene", json_textos(&nomes))
            .campo("rank", json_lista(&v, |g| sel(g).0.to_string()))
            .campo("score", json_lista(&v, |g| json_num(sel(g).1)))
            .campo("minus_log10_score", json_lista(&v, |g| json_num(-sel(g).1.max(f64::MIN_POSITIVE).log10())))
            .campo("lfc", json_lista(&v, |g| json_num(sel(g).2)))
            .campo("fdr", json_lista(&v, |g| json_num(sel(g).3)))
            .fim()
    }

    /// Everything, the per-guide table included (for a genome-wide library
    /// that is tens of MB: see `json_compact`).
    pub fn json(&self) -> String {
        self.json_com(true)
    }

    /// Everything but the per-guide table (`guides`): the summary, the gene
    /// table and the rank plot -- what a screen of results needs.
    pub fn json_compact(&self) -> String {
        self.json_com(false)
    }

    fn json_com(&self, com_guias: bool) -> String {
        let g = &self.guides;
        let guias = if !com_guias {
            "null".to_string()
        } else {
            Objeto::novo()
                .campo("sgrna", json_lista(g, |x| json_texto(&x.sgrna)))
                .campo("gene", json_lista(g, |x| json_texto(&x.gene)))
                .campo("control_counts", json_lista(g, |x| json_nums(&x.control_counts)))
                .campo("treatment_counts", json_lista(g, |x| json_nums(&x.treatment_counts)))
                .campo("control_mean", json_lista(g, |x| json_num(x.control_mean)))
                .campo("treatment_mean", json_lista(g, |x| json_num(x.treatment_mean)))
                .campo("lfc", json_lista(g, |x| json_num(x.lfc)))
                .campo("control_var", json_lista(g, |x| json_num(x.control_var)))
                .campo("adj_var", json_lista(g, |x| json_num(x.adj_var)))
                .campo("score", json_lista(g, |x| json_num(x.score)))
                .campo("theta", json_lista(g, |x| json_num(x.theta)))
                .campo("p_low", json_lista(g, |x| json_num(x.p_low)))
                .campo("p_high", json_lista(g, |x| json_num(x.p_high)))
                .campo("p_twosided", json_lista(g, |x| json_num(x.p_twosided)))
                .campo("fdr", json_lista(g, |x| json_num(x.fdr)))
                .campo("high_in_treatment", json_lista(g, |x| x.high_in_treatment.to_string()))
                .campo("used", json_lista(g, |x| x.used.to_string()))
                .fim()
        };
        let n = &self.genes;
        let genes = Objeto::novo()
            .campo("id", json_lista(n, |x| json_texto(&x.gene)))
            .campo("num", json_lista(n, |x| x.guides.to_string()))
            .campo("neg_score", json_lista(n, |x| json_num(x.neg_score)))
            .campo("neg_p", json_lista(n, |x| json_num(x.neg_p)))
            .campo("neg_fdr", json_lista(n, |x| json_num(x.neg_fdr)))
            .campo("neg_rank", json_lista(n, |x| x.neg_rank.to_string()))
            .campo("neg_goodsgrna", json_lista(n, |x| x.neg_good.to_string()))
            .campo("neg_lfc", json_lista(n, |x| json_num(x.neg_lfc)))
            .campo("pos_score", json_lista(n, |x| json_num(x.pos_score)))
            .campo("pos_p", json_lista(n, |x| json_num(x.pos_p)))
            .campo("pos_fdr", json_lista(n, |x| json_num(x.pos_fdr)))
            .campo("pos_rank", json_lista(n, |x| x.pos_rank.to_string()))
            .campo("pos_goodsgrna", json_lista(n, |x| x.pos_good.to_string()))
            .campo("pos_lfc", json_lista(n, |x| json_num(x.pos_lfc)))
            .fim();
        let h = self.hits();
        let top = |pos: bool| json_lista(&self.top(10, pos), |g| json_texto(&g.gene));
        Objeto::novo()
            .campo("count_table", json_texto(&self.count_table))
            .campo("treatment", json_textos(&self.treatment))
            .campo("control", json_textos(&self.control))
            .campo("norm_method", json_texto(&self.norm_method))
            .campo("norm_used", json_texto(&self.norm_used))
            .campo("size_factors", json_lista(&self.size_factors, |(s, f)| format!("[{},{}]", json_texto(s), json_num(*f))))
            .campo("model_k", json_num(self.model_k))
            .campo("model_b", json_num(self.model_b))
            .campo("model_guides", self.model_guides.to_string())
            .campo("model_outliers", self.model_outliers.to_string())
            .campo("guides_total", self.guides_total.to_string())
            .campo("guides_used", self.guides_used.to_string())
            .campo("genes_total", self.genes_total.to_string())
            .campo("control_guides", self.control_guides.to_string())
            .campo("alpha_neg", json_num(self.alpha_neg))
            .campo("alpha_pos", json_num(self.alpha_pos))
            .campo("permutations", self.permutations.to_string())
            .campo("random_genes", self.random_genes.to_string())
            .campo("neg_fdr_below", format!("{{\"0.01\":{},\"0.05\":{},\"0.25\":{}}}", h[0].0, h[1].0, h[2].0))
            .campo("pos_fdr_below", format!("{{\"0.01\":{},\"0.05\":{},\"0.25\":{}}}", h[0].1, h[1].1, h[2].1))
            .campo("skipped_genes", json_textos(&self.skipped_genes))
            .campo("top_negative", top(false))
            .campo("top_positive", top(true))
            .campo("rank_plot", format!("{{\"neg\":{},\"pos\":{}}}", self.grafico(false), self.grafico(true)))
            .campo("genes", genes)
            .campo("guides", guias)
            .campo("outputs", json_textos(&self.outputs))
            .campo("warnings", json_textos(&self.warnings))
            .fim()
    }

    pub fn resumo(&self) -> String {
        let mut o = format!(
            "CRISPR screen test: {} vs {} ({})\n",
            self.treatment.join(","),
            self.control.join(","),
            self.count_table
        );
        o.push_str(&format!(
            "  {} guides, {} entered the test ({} left out: no reads), {} genes\n",
            self.guides_total,
            self.guides_used,
            self.guides_total - self.guides_used,
            self.genes_total
        ));
        let f: Vec<String> = self.size_factors.iter().map(|(s, x)| format!("{s} {x:.4}")).collect();
        let modo = if self.norm_used != self.norm_method {
            format!("{} (asked: {})", self.norm_used, self.norm_method)
        } else {
            self.norm_used.clone()
        };
        o.push_str(&format!("  normalization: {modo}; size factors {}\n", f.join(", ")));
        o.push_str(&format!(
            "  variance model: var = mean + mean^{:.4} * 2^{:.4} (fitted on {} guides, {} high-count outliers left out)\n",
            self.model_k, self.model_b, self.model_guides, self.model_outliers
        ));
        o.push_str(&format!(
            "  RRA: alpha {:.4} (negative), {:.4} (positive); {} random genes per gene size\n",
            self.alpha_neg, self.alpha_pos, self.random_genes
        ));
        if self.control_guides > 0 {
            o.push_str(&format!("  null distribution from {} control guides\n", self.control_guides));
        }
        let h = self.hits();
        for (pos, nome) in [(false, "negative (depleted)"), (true, "positive (enriched)")] {
            let (a, b, c) = if pos { (h[0].1, h[1].1, h[2].1) } else { (h[0].0, h[1].0, h[2].0) };
            let top: Vec<String> = self
                .top(5, pos)
                .iter()
                .map(|g| format!("{} (FDR {})", g.gene, py_g(if pos { g.pos_fdr } else { g.neg_fdr }, 3)))
                .collect();
            o.push_str(&format!(
                "  {nome} selection: {a} gene(s) at FDR < 0.01, {b} < 0.05, {c} < 0.25; top: {}\n",
                top.join(", ")
            ));
        }
        for p in &self.outputs {
            o.push_str(&format!("  written: {p}\n"));
        }
        for a in &self.warnings {
            o.push_str(&format!("  warning: {a}\n"));
        }
        o
    }

    /// The text of `*.sgrna_summary.txt`, as MAGeCK writes it.
    pub fn sgrna_summary(&self) -> String {
        let mut o = String::from(
            "sgrna\tGene\tcontrol_count\ttreatment_count\tcontrol_mean\ttreat_mean\tLFC\tcontrol_var\tadj_var\tscore\tp.low\tp.high\tp.twosided\tFDR\thigh_in_treatment\n",
        );
        for g in &self.guides {
            let cc: Vec<String> = g.control_counts.iter().map(|&x| py_g(x, 5)).collect();
            let tc: Vec<String> = g.treatment_counts.iter().map(|&x| py_g(x, 5)).collect();
            o.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                g.sgrna,
                g.gene,
                cc.join("/"),
                tc.join("/"),
                py_g(g.control_mean, 5),
                py_g(g.treatment_mean, 5),
                py_g(g.lfc, 5),
                py_g(g.control_var, 5),
                py_g(g.adj_var, 5),
                py_g(g.score, 5),
                py_g(g.p_low, 5),
                py_g(g.p_high, 5),
                py_g(g.p_twosided, 5),
                py_g(g.fdr, 5),
                if g.high_in_treatment { "True" } else { "False" }
            ));
        }
        o
    }

    /// The text of `*.gene_summary.txt`, as MAGeCK writes it.
    pub fn gene_summary(&self) -> String {
        let mut o = String::from(
            "id\tnum\tneg|score\tneg|p-value\tneg|fdr\tneg|rank\tneg|goodsgrna\tneg|lfc\tpos|score\tpos|p-value\tpos|fdr\tpos|rank\tpos|goodsgrna\tpos|lfc\n",
        );
        for g in &self.genes {
            o.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                g.gene,
                g.guides,
                fmt_e4(g.neg_score),
                fmt_e4(g.neg_p),
                fmt_f6(g.neg_fdr),
                g.neg_rank,
                g.neg_good,
                py_g(g.neg_lfc, 5),
                fmt_e4(g.pos_score),
                fmt_e4(g.pos_p),
                fmt_f6(g.pos_fdr),
                g.pos_rank,
                g.pos_good,
                py_g(g.pos_lfc, 5)
            ));
        }
        o
    }
}

/// Sample columns by label, or by 0-based index when no label matches.
fn resolve_amostras(spec: &[String], t: &CountTable, papel: &str) -> Result<Vec<usize>, String> {
    let mut v = Vec::new();
    for s in spec {
        let s = s.trim();
        let i = match t.samples.iter().position(|x| x == s) {
            Some(i) => i,
            None => match s.parse::<usize>() {
                Ok(i) if i < t.samples.len() => i,
                _ => {
                    return Err(format!(
                        "{}: no sample `{s}` (given as {papel}). The samples are: {}. Use these names, or 0-based column numbers among them.",
                        t.path,
                        t.samples.join(", ")
                    ))
                }
            },
        };
        if v.contains(&i) {
            return Err(format!("The {papel} sample `{}` is given twice.", t.samples[i]));
        }
        v.push(i);
    }
    Ok(v)
}

/// MAGeCK's median-ratio factor of each column, over the rows `linhas`.
fn fator_mediana(colunas: &[Vec<f64>], linhas: &[usize]) -> Vec<f64> {
    let n = colunas.len() as f64;
    let gm: Vec<(usize, f64)> = linhas
        .iter()
        .filter(|&&g| colunas.iter().map(|c| c[g]).sum::<f64>() > 0.0)
        .map(|&g| (g, (colunas.iter().map(|c| (c[g] + 1.0).ln()).sum::<f64>() / n).exp()))
        .collect();
    colunas
        .iter()
        .map(|c| {
            if gm.is_empty() {
                return 0.0;
            }
            let mut r: Vec<f64> = gm.iter().map(|&(g, m)| c[g] / m).collect();
            ordena_f64(&mut r);
            let x = r[r.len() / 2];
            if x > 0.0 {
                1.0 / x
            } else {
                0.0
            }
        })
        .collect()
}

/// MAGeCK's `normalizeCounts` factors. Returns (factors, method applied).
fn fatores(
    colunas: &[Vec<f64>],
    linhas: &[usize],
    metodo: NormMethod,
    rotulos: &[String],
    avisos: &mut Vec<String>,
) -> Result<(Vec<f64>, &'static str), String> {
    let n = colunas.len();
    let m = linhas.len();
    let total = || -> Result<Vec<f64>, String> {
        let somas: Vec<f64> = colunas.iter().map(|c| linhas.iter().map(|&g| c[g]).sum::<f64>()).collect();
        if let Some(i) = somas.iter().position(|&s| s <= 0.0) {
            return Err(format!(
                "Sample `{}` has no reads in the guides used for normalization, so its counts cannot be normalized. Check the sample, or leave it out.",
                rotulos[i]
            ));
        }
        let media = somas.iter().sum::<f64>() / n as f64;
        Ok(somas.iter().map(|s| media / s).collect())
    };
    let (f, usado): (Vec<f64>, &'static str) = match metodo {
        NormMethod::None => (vec![1.0; n], "none"),
        NormMethod::Total => (total()?, "total"),
        NormMethod::Median | NormMethod::Control => {
            let med = fator_mediana(colunas, linhas);
            let mut usa_total = false;
            for (s, c) in colunas.iter().enumerate() {
                if med[s] == 0.0 {
                    avisos.push(format!(
                        "sample {}: the median count ratio is zero, so median normalization is impossible; switched to total read count normalization.",
                        rotulos[s]
                    ));
                    usa_total = true;
                }
                let zf = linhas.iter().filter(|&&g| c[g] == 0.0).count() as f64 / m as f64;
                if zf > 0.45 {
                    avisos.push(format!(
                        "sample {}: {:.1}% of the guides have no reads, too many for median normalization; switched to total read count normalization.",
                        rotulos[s],
                        100.0 * zf
                    ));
                    usa_total = true;
                } else if zf > 0.3 {
                    avisos.push(format!(
                        "sample {}: {:.1}% of the guides have no reads; median normalization may be unstable (consider norm_method=\"total\").",
                        rotulos[s],
                        100.0 * zf
                    ));
                }
            }
            if usa_total {
                (total()?, "total")
            } else if metodo == NormMethod::Control {
                (med, "control")
            } else {
                (med, "median")
            }
        }
    };
    if f.iter().any(|&x| !(0.2..=5.0).contains(&x)) {
        avisos.push(
            "some size factors are below 0.2 or above 5: the samples differ a lot in depth or composition; check the normalization.".to_string(),
        );
    }
    Ok((f, usado))
}

/// MAGeCK's `leastsquare` with weights: y = a + b x; returns (b, a).
pub fn least_squares(x: &[f64], y: &[f64], w: &[f64]) -> (f64, f64) {
    let nw: f64 = w.iter().sum();
    let sy: f64 = y.iter().zip(w).map(|(y, w)| y * w).sum();
    let sx: f64 = x.iter().zip(w).map(|(x, w)| x * w).sum();
    let sx2: f64 = x.iter().zip(w).map(|(x, w)| x * x * w).sum();
    let sxy: f64 = x.iter().zip(y).zip(w).map(|((x, y), w)| x * y * w).sum();
    let den = nw * sx2 - sx * sx;
    let a = (sy * sx2 - sx * sxy) / den;
    let b = (nw * sxy - sx * sy) / den;
    (b, a)
}

struct Modelo {
    k: f64,
    b: f64,
    guias: usize,
    fora: usize,
}

/// MAGeCK's `modelmeanvar(..., 'linear')`. `valores[g]`: the counts of guide
/// g in the samples used for the variance.
fn modelo_variancia(valores: &[Vec<f64>], avisos: &mut Vec<String>) -> Modelo {
    let media: Vec<f64> = valores.iter().map(|v| median(v)).collect();
    let var: Vec<f64> = valores.iter().zip(&media).map(|(v, &m)| variancia(v, m)).collect();
    let corte = if media.len() > 10 {
        let mm = median(&media);
        mm + 4.0 * variancia(&media, mm).powf(0.5)
    } else {
        media.iter().cloned().fold(f64::NEG_INFINITY, f64::max) * 2.0
    };
    let fora = media.iter().filter(|&&x| x > corte).count();
    let mut mg = Vec::new();
    let mut vg = Vec::new();
    for (&m, &v) in media.iter().zip(&var) {
        if m < v && m < corte {
            mg.push(m);
            vg.push(v - m);
        }
    }
    if mg.len() < 10 {
        avisos.push(
            "the counts vary too little between the replicates to fit the variance model (if they are technical replicates, merge them into one sample); all guides were used."
                .to_string(),
        );
        mg = media.clone();
        vg = var.iter().zip(&media).map(|(&v, &m)| if v - m > 0.01 { v - m } else { 0.01 }).collect();
    }
    let x: Vec<f64> = mg.iter().map(|&m| log2_py(m + 1.0)).collect();
    let y: Vec<f64> = vg.iter().map(|&v| log2_py(v + 1.0)).collect();
    let (mut k, mut b) = least_squares(&x, &y, &mg);
    if !k.is_finite() || !b.is_finite() {
        avisos.push("the variance model could not be fitted (too few distinct counts); var = 2 * mean was used.".to_string());
        k = 1.0;
        b = 0.0;
    }
    if k < 1.0 {
        k = 1.0;
    }
    if b < 0.0 {
        b = 0.0;
    }
    Modelo { k, b, guias: mg.len(), fora }
}

/// MAGeCK's `getNormalPValue`: a normal distribution truncated at zero.
fn p_normal(m0: &[f64], v0: &[f64], m1: &[f64], baixo: bool) -> Result<Vec<f64>, String> {
    let min1 = m1.iter().cloned().filter(|&x| x > 0.0).fold(f64::INFINITY, f64::min);
    if !min1.is_finite() {
        return Err("Every guide has a treatment median of zero: there is nothing to test. Check the treatment samples.".to_string());
    }
    Ok((0..m0.len())
        .map(|i| {
            let a = if m1[i] > min1 { m1[i] } else { min1 };
            let th = (a - m0[i]) / v0[i].sqrt();
            let th0 = (0.0 - m0[i]) / v0[i].sqrt();
            let tp = if baixo { norm_cdf(th) } else { norm_sf(th) };
            let tp0 = norm_cdf(th0);
            if baixo {
                (tp - tp0) / (1.0 - tp0)
            } else {
                tp / (1.0 - tp0)
            }
        })
        .collect())
}

/// MAGeCK's `calculate_gene_lfc`.
fn lfc_genes(lfc: &[f64], ordem: &[usize], valido: &[bool], genes: &[String], n_corte: usize, metodo: GeneLfc, positivo: bool) -> HashMap<String, f64> {
    let mut por: HashMap<&str, Vec<f64>> = HashMap::new();
    let alfa = matches!(metodo, GeneLfc::AlphaMean | GeneLfc::AlphaMedian);
    for (ni, &g) in ordem.iter().enumerate() {
        if !valido[g] {
            continue;
        }
        let e = por.entry(genes[g].as_str()).or_default();
        // MAGeCK counts the position from 1 and keeps it while <= the cut
        if !alfa || ni < n_corte {
            e.push(lfc[g]);
        }
    }
    por.into_iter()
        .map(|(gene, v)| {
            let x = match metodo {
                GeneLfc::Median | GeneLfc::AlphaMedian => median(&v),
                GeneLfc::Mean | GeneLfc::AlphaMean => {
                    if v.is_empty() {
                        0.0
                    } else {
                        v.iter().sum::<f64>() / v.len() as f64
                    }
                }
                GeneLfc::SecondBest => {
                    let mut s = v.clone();
                    ordena_f64(&mut s);
                    if positivo {
                        s.reverse();
                    }
                    if s.len() > 1 {
                        s[1]
                    } else {
                        0.0
                    }
                }
            };
            (gene.to_string(), x)
        })
        .collect()
}

/// Tests a count table file (`mageck test -k`): treatment vs control.
///
/// `treatment` and `control` are sample names (or 0-based column numbers);
/// an empty `control` means every other sample. With `output_prefix`,
/// writes `<prefix>.sgrna_summary.txt` and `<prefix>.gene_summary.txt` in
/// MAGeCK's format.
pub fn crispr_test(
    count_table: &str,
    treatment: &[String],
    control: &[String],
    output_prefix: Option<&str>,
    opts: &TestOptions,
) -> Result<TestResult, String> {
    let t = read_count_table(count_table)?;
    crispr_test_table(&t, treatment, control, output_prefix, opts)
}

/// `crispr_test` on a table already in memory.
pub fn crispr_test_table(
    t: &CountTable,
    treatment: &[String],
    control: &[String],
    output_prefix: Option<&str>,
    o: &TestOptions,
) -> Result<TestResult, String> {
    if !(o.gene_test_fdr_threshold > 0.0 && o.gene_test_fdr_threshold <= 1.0) {
        return Err(format!("gene_test_fdr_threshold={} must be in (0, 1]; MAGeCK's default is 0.25.", o.gene_test_fdr_threshold));
    }
    if o.remove_zero_threshold.is_nan() || o.remove_zero_threshold < 0.0 {
        return Err(format!("remove_zero_threshold={} must be >= 0.", o.remove_zero_threshold));
    }
    if o.permutations == Some(0) {
        return Err("permutations=0: give at least 1 round, or None for RRA's rule (at least 100,000 random genes).".to_string());
    }
    if o.max_sgrna_per_gene_permutation == 0 {
        return Err("max_sgrna_per_gene_permutation must be at least 1 (RRA's default is 100).".to_string());
    }
    if o.seed == 0 || o.seed as u64 >= PM_M {
        return Err(format!("seed={} must be between 1 and 2147483646 (RRA's is 123456).", o.seed));
    }
    let trat = resolve_amostras(treatment, t, "treatment")?;
    if trat.is_empty() {
        return Err("No treatment sample given: name at least one column of the count table.".to_string());
    }
    let mut ctrl = resolve_amostras(control, t, "control")?;
    if ctrl.is_empty() {
        ctrl = (0..t.samples.len()).filter(|i| !trat.contains(i)).collect();
        if ctrl.is_empty() {
            return Err("No control sample: every column of the table is a treatment. Name the control samples.".to_string());
        }
    }
    if let Some(&i) = trat.iter().find(|i| ctrl.contains(i)) {
        return Err(format!("Sample `{}` is both treatment and control: a sample can only be in one group.", t.samples[i]));
    }
    if let Some(p) = output_prefix {
        for suf in [".sgrna_summary.txt", ".gene_summary.txt"] {
            if mesmo_arquivo(&format!("{p}{suf}"), &t.path) {
                return Err(format!("The output `{p}{suf}` is the count table itself. Choose another output prefix."));
            }
        }
    }
    let ng = t.ids.len();
    let mut avisos: Vec<String> = Vec::new();

    // Control guides.
    if !o.control_sgrnas.is_empty() && !o.control_genes.is_empty() {
        return Err("Give control_sgrnas or control_genes, not both.".to_string());
    }
    let mut controle = vec![false; ng];
    if !o.control_genes.is_empty() {
        let cg: HashSet<&str> = o.control_genes.iter().map(String::as_str).collect();
        for (c, gene) in controle.iter_mut().zip(&t.genes) {
            *c = cg.contains(gene.as_str());
        }
        let n = controle.iter().filter(|&&c| c).count();
        if n < 2 {
            return Err(format!(
                "Only {n} guide(s) of the table belong to the control genes ({}). At least 2 control guides are needed: check the gene names.",
                mostra(&o.control_genes.join(", "), 80)
            ));
        }
    } else if !o.control_sgrnas.is_empty() {
        let cs: HashSet<&str> = o.control_sgrnas.iter().map(String::as_str).collect();
        for (c, id) in controle.iter_mut().zip(&t.ids) {
            *c = cs.contains(id.as_str());
        }
        let n = controle.iter().filter(|&&c| c).count();
        if n < 2 {
            return Err(format!(
                "Only {n} of the {} control sgRNA ids are in the count table. At least 2 are needed: check that the ids match the first column of the table.",
                cs.len()
            ));
        }
        let mut por_gene: HashMap<&str, (usize, usize)> = HashMap::new();
        for (gene, &c) in t.genes.iter().zip(&controle) {
            let e = por_gene.entry(gene.as_str()).or_insert((0, 0));
            e.0 += 1;
            if c {
                e.1 += 1;
            }
        }
        let mut mistos: Vec<&str> = por_gene.iter().filter(|(_, &(n, c))| c > 0 && c < n).map(|(g, _)| *g).collect();
        mistos.sort_unstable();
        if let Some(g) = mistos.first() {
            return Err(format!(
                "Gene `{g}` has both control and non-control guides. A gene's guides must be all controls or none (MAGeCK refuses this too): fix the control list."
            ));
        }
    }
    let n_controle = controle.iter().filter(|&&c| c).count();
    if o.norm_method == NormMethod::Control && n_controle == 0 {
        return Err("norm_method=\"control\" needs control_sgrnas or control_genes.".to_string());
    }

    // Normalization of the selected samples: control, then treatment.
    let sel: Vec<usize> = ctrl.iter().chain(trat.iter()).copied().collect();
    let rotulos: Vec<String> = sel.iter().map(|&s| t.samples[s].clone()).collect();
    let colunas: Vec<Vec<f64>> = sel.iter().map(|&s| t.counts[s].clone()).collect();
    let linhas: Vec<usize> = if o.norm_method == NormMethod::Control {
        (0..ng).filter(|&g| controle[g]).collect()
    } else {
        (0..ng).collect()
    };
    let (fat, usado) = fatores(&colunas, &linhas, o.norm_method, &rotulos, &mut avisos)?;
    let norm: Vec<Vec<f64>> = colunas.iter().zip(&fat).map(|(c, &f)| c.iter().map(|&x| f * x).collect()).collect();
    let (nc, nt) = (ctrl.len(), trat.len());

    // Guide statistics.
    let por_guia = |cols: &[Vec<f64>], g: usize| -> Vec<f64> { cols.iter().map(|c| c[g]).collect() };
    let cols_var: &[Vec<f64>] = if nc > 1 { &norm[..nc] } else { &norm[..] };
    let vals_var: Vec<Vec<f64>> = (0..ng).map(|g| por_guia(cols_var, g)).collect();
    let modelo = modelo_variancia(&vals_var, &mut avisos);
    let vals_c: Vec<Vec<f64>> = (0..ng).map(|g| por_guia(&norm[..nc], g)).collect();
    let vals_t: Vec<Vec<f64>> = (0..ng).map(|g| por_guia(&norm[nc..], g)).collect();
    let cmed: Vec<f64> = vals_c.iter().map(|v| median(v)).collect();
    let valido1: Vec<bool> = cmed.iter().map(|&x| x > o.remove_zero_threshold).collect();
    let cmin = cmed.iter().cloned().filter(|&x| x > 0.0).fold(f64::INFINITY, f64::min);
    if !cmin.is_finite() {
        return Err("Every guide has a control median of zero: the control samples have no reads. Check the control samples.".to_string());
    }
    let cmean: Vec<f64> = cmed.iter().map(|&x| if x > cmin { x } else { cmin }).collect();
    let var_c: Vec<f64> = vals_var.iter().map(|v| variancia(v, median(v))).collect();
    let var_t: Vec<f64> = vals_t.iter().map(|v| variancia(v, median(v))).collect();
    let cvar: Vec<f64> = if nc > 1 || (nc == 1 && nt == 1) {
        var_c.clone()
    } else {
        let fc = (nc as f64 - 1.0) / (nc as f64 - 1.0 + nt as f64 - 1.0);
        let ft = 1.0 - fc;
        var_c.iter().zip(&var_t).map(|(&a, &b)| a * fc + b * ft).collect()
    };
    let adjvar: Vec<f64> = cmean.iter().map(|&z| z.powf(modelo.k) * 2f64.powf(modelo.b) + z).collect();
    let tmean: Vec<f64> = vals_t.iter().map(|v| median(v)).collect();
    let valido2: Vec<bool> = tmean.iter().map(|&x| x > o.remove_zero_threshold).collect();
    let valido: Vec<bool> = (0..ng)
        .map(|g| match o.remove_zero {
            RemoveZero::None => true,
            RemoveZero::Control => valido1[g],
            RemoveZero::Treatment => valido2[g],
            RemoveZero::Any => valido1[g] && valido2[g],
            RemoveZero::Both => valido1[g] || valido2[g],
        })
        .collect();
    let n_validos = valido.iter().filter(|&&v| v).count();
    if n_validos == 0 {
        return Err(format!(
            "No guide is left after removing the guides without reads (remove_zero=\"{}\"): check the samples, or use remove_zero=\"none\".",
            o.remove_zero.name()
        ));
    }
    let p_low = p_normal(&cmean, &adjvar, &tmean, true)?;
    let p_high = p_normal(&cmean, &adjvar, &tmean, false)?;
    let theta0: Vec<f64> = (0..ng).map(|g| (tmean[g] - cmean[g]) / adjvar[g].sqrt()).collect();
    let menor = p_low.iter().map(|&x| x * 0.1).filter(|&x| x > 0.0).fold(f64::INFINITY, f64::min);
    let menor = if menor.is_finite() { menor } else { f64::MIN_POSITIVE };
    let theta: Vec<f64> = (0..ng)
        .map(|g| {
            if tmean[g] < cmean[g] {
                qnorm_mageck(p_low[g].max(menor).ln()).unwrap_or(theta0[g])
            } else {
                theta0[g]
            }
        })
        .collect();
    let abs_theta: Vec<f64> = theta.iter().map(|x| x.abs()).collect();
    let p2: Vec<f64> = (0..ng).map(|g| if p_low[g] < p_high[g] { 2.0 * p_low[g] } else { 2.0 * p_high[g] }).collect();
    let fdr = bh_fdr(&p2);
    let lfc: Vec<f64> = (0..ng).map(|g| log2_py(tmean[g] + 1.0) - log2_py(cmean[g] + 1.0)).collect();

    // Gene test, both directions.
    let limiar = o.gene_test_fdr_threshold;
    let ordem_neg = ordem_por(&theta, true);
    let n_low = p_low.iter().filter(|&&p| p <= limiar).count();
    let n_low_v = (0..ng).filter(|&g| p_low[g] <= limiar && valido[g]).count();
    let alpha_neg = (n_low_v as f64 / n_validos as f64).clamp(0.05, 0.5);
    let lfc_neg = lfc_genes(&lfc, &ordem_neg, &valido, &t.genes, n_low, o.gene_lfc, false);
    let rra_neg = rra(&theta, &ordem_neg, &valido, &t.genes, &t.ids, &controle, alpha_neg, o)?;

    let menos_theta: Vec<f64> = theta.iter().map(|x| -1.0 * x).collect();
    let ordem_pos = ordem_por(&menos_theta, true);
    let n_high = p_high.iter().filter(|&&p| p <= limiar).count();
    let n_high_v = (0..ng).filter(|&g| p_high[g] <= limiar && valido[g]).count();
    let alpha_pos = if n_high > 0 { n_high_v as f64 / n_validos as f64 } else { 0.01 };
    let lfc_pos = lfc_genes(&lfc, &ordem_pos, &valido, &t.genes, n_high, o.gene_lfc, true);
    let rra_pos = rra(&menos_theta, &ordem_pos, &valido, &t.genes, &t.ids, &controle, alpha_pos, o)?;

    let mut pos_de: HashMap<&str, (usize, &GeneRra)> = HashMap::new();
    for (i, g) in rra_pos.genes.iter().enumerate() {
        pos_de.insert(g.nome.as_str(), (i + 1, g));
    }
    let mut genes: Vec<GeneStat> = Vec::with_capacity(rra_neg.genes.len());
    for (i, g) in rra_neg.genes.iter().enumerate() {
        let mut e = GeneStat {
            gene: g.nome.clone(),
            guides: g.itens,
            neg_score: g.lo,
            neg_p: g.p,
            neg_fdr: g.fdr,
            neg_rank: i + 1,
            neg_good: g.bons,
            neg_lfc: lfc_neg.get(&g.nome).copied().unwrap_or(0.0),
            ..Default::default()
        };
        if let Some(&(r, p)) = pos_de.get(g.nome.as_str()) {
            e.pos_score = p.lo;
            e.pos_p = p.p;
            e.pos_fdr = p.fdr;
            e.pos_rank = r;
            e.pos_good = p.bons;
            e.pos_lfc = lfc_pos.get(&p.nome).copied().unwrap_or(0.0);
        }
        genes.push(e);
    }
    if o.sort_by_positive {
        let mut idx: Vec<usize> = (0..genes.len()).collect();
        multi_merge_sort_stable_by(&mut idx, |&a: &usize, &b: &usize| genes[a].pos_rank.cmp(&genes[b].pos_rank));
        let mut velhos: Vec<Option<GeneStat>> = genes.into_iter().map(Some).collect();
        genes = idx.into_iter().map(|i| velhos[i].take().unwrap_or_default()).collect();
    }

    let ordem_guias = ordem_por(&abs_theta, false);
    let guides: Vec<GuideStat> = ordem_guias
        .iter()
        .map(|&g| GuideStat {
            sgrna: t.ids[g].clone(),
            gene: t.genes[g].clone(),
            control_counts: vals_c[g].clone(),
            treatment_counts: vals_t[g].clone(),
            control_mean: cmean[g],
            treatment_mean: tmean[g],
            lfc: lfc[g],
            control_var: cvar[g],
            adj_var: adjvar[g],
            score: abs_theta[g],
            theta: theta[g],
            p_low: p_low[g],
            p_high: p_high[g],
            p_twosided: p2[g],
            fdr: fdr[g],
            high_in_treatment: tmean[g] > cmean[g],
            used: valido[g],
        })
        .collect();

    let mut r = TestResult {
        count_table: t.path.clone(),
        treatment: trat.iter().map(|&i| t.samples[i].clone()).collect(),
        control: ctrl.iter().map(|&i| t.samples[i].clone()).collect(),
        norm_method: o.norm_method.name().to_string(),
        norm_used: usado.to_string(),
        size_factors: rotulos.into_iter().zip(fat).collect(),
        model_k: modelo.k,
        model_b: modelo.b,
        model_guides: modelo.guias,
        model_outliers: modelo.fora,
        guides_total: ng,
        guides_used: n_validos,
        genes_total: genes.len(),
        control_guides: n_controle,
        alpha_neg,
        alpha_pos,
        permutations: rra_neg.passes,
        random_genes: rra_neg.aleatorios,
        guides,
        genes,
        outputs: Vec::new(),
        skipped_genes: o.skip_genes.clone(),
        warnings: avisos,
    };
    if let Some(p) = output_prefix {
        let a = format!("{p}.sgrna_summary.txt");
        let b = format!("{p}.gene_summary.txt");
        grava_texto(&a, &r.sgrna_summary())?;
        grava_texto(&b, &r.gene_summary())?;
        r.outputs = vec![a, b];
    }
    Ok(r)
}

fn grava_texto(caminho: &str, texto: &str) -> Result<(), String> {
    let mut g = Saida { caminho: caminho.to_string(), armada: false };
    let mut w = Escritor::cria(caminho)?;
    g.armada = true;
    w.write_all(texto.as_bytes()).map_err(|e| format!("Error writing {caminho}: {e}"))?;
    w.termina().map_err(|e| format!("Error writing {caminho}: {e}"))?;
    g.armada = false;
    Ok(())
}

// Internal pieces the trials check against hand-worked values.
#[doc(hidden)]
pub mod internals {
    pub fn percentile(v: f64, sorted: &[f64]) -> f64 {
        super::meio_posto(v, sorted) / (sorted.len() * 2) as f64
    }
    pub fn perm_p(lo: f64, sorted_random: &[f64]) -> f64 {
        super::meio_posto(lo, sorted_random) / 2.0 / sorted_random.len() as f64
    }
    pub fn park_miller(seed: &mut u64) -> f64 {
        super::pm_proximo(seed)
    }
    pub fn park_miller_jump(seed: u64, n: u64) -> u64 {
        super::pm_salta(seed, n)
    }
    pub fn random_lo(seed: &mut u64, size: usize, total: usize, alpha: f64, null: Option<&[f64]>) -> Vec<f64> {
        super::lo_aleatorios(seed, size, total, alpha, null)
    }
    pub fn qnorm_table(ln_p: f64) -> Option<f64> {
        super::qnorm_mageck(ln_p)
    }
    pub fn bh_rra(p: &[f64]) -> Vec<f64> {
        super::bh_rra(p)
    }
    /// (factors, method applied, warnings).
    pub fn size_factors(cols: &[Vec<f64>], method: super::NormMethod) -> Result<(Vec<f64>, String, Vec<String>), String> {
        let linhas: Vec<usize> = (0..cols.first().map(|c| c.len()).unwrap_or(0)).collect();
        let rot: Vec<String> = (0..cols.len()).map(|i| format!("s{i}")).collect();
        let mut av = Vec::new();
        let (f, u) = super::fatores(cols, &linhas, method, &rot, &mut av)?;
        Ok((f, u.to_string(), av))
    }
    /// (k, b, guides fitted, outliers).
    pub fn variance_model(values: &[Vec<f64>]) -> (f64, f64, usize, usize) {
        let mut av = Vec::new();
        let m = super::modelo_variancia(values, &mut av);
        (m.k, m.b, m.guias, m.fora)
    }
    pub fn auto_offsets(probe: &[(usize, u64)], total: u64) -> Vec<usize> {
        super::escolhe_offsets(probe, total)
    }
}
