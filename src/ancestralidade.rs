// =====================================================================
//  ANCESTRY: FREQUENCY PANELS, SUPERVISED ADMIXTURE, POLYGENIC SCORES
//
//  Three tools over the genotypes of a VCF, for teaching and research:
//
//   1. PANEL (`build_panel`, `panel_from_frequencies`). For each biallelic
//      SNV, the frequency of the ALT allele in each group (population) and
//      how many alleles were called in each group. Built from a multi-sample
//      VCF plus a sample -> group table (the 1000 Genomes superpopulations,
//      say), or from a sites VCF whose INFO already holds the frequencies of
//      each group (AFR_AF, EUR_AF...). Written as a tab-separated table.
//
//   2. ANCESTRY (`estimate_ancestry`). Supervised admixture: the panel's
//      frequencies are FIXED (the P matrix of ADMIXTURE), and the proportions
//      q of one person maximise the binomial likelihood of the person's
//      genotypes. EM with the update of ADMIXTURE / fastNGSadmix for fixed P,
//      accelerated by SQUAREM, then Newton steps on the simplex (with an
//      active set for the groups at zero) to land exactly on the maximum.
//      With P fixed the log-likelihood is concave in q, so the maximum found
//      is the global one. Standard errors come from a block bootstrap (blocks
//      of a few Mb resampled with replacement).
//
//   3. POLYGENIC SCORE (`polygenic_score`). A PGS Catalog scoring file
//      applied to a person: the sum of effect-allele dosage x weight. With a
//      panel, the variants the person's file cannot call are filled with
//      2 x frequency, and each group gets the expected mean and variance of
//      the score (Hardy-Weinberg, independent variants), a z-score and a
//      normal percentile.
//
//  A PERSON is a VCF, optionally a sample in it, and optionally a BED of
//  callable regions. GIAB benchmark files list only the variants: a site
//  absent from the VCF is homozygous reference only where the BED says the
//  genome was assessed. With a BED, nothing outside it is used, not even a
//  record (the benchmark is not reliable there). Without a BED,
//  `absent_is_ref` decides (and only on chromosomes the VCF has records for).
//
//  JOIN. The panel (or the score) is held in memory, grouped by chromosome
//  and sorted by position; the person's VCF is read once, line by line, and
//  merged against it. Chromosome names are compared in canonical form
//  (chr1 == 1, chrM == MT). The VCF must be sorted by position with each
//  chromosome in a single block; the chromosome order itself is free, so a
//  VCF in header order (1..22, X, Y, M) and one sorted by
//  `sort_genomic(..., chrom_order="natural")` (M before X) both work. The
//  panel is written in the order 1..22, X, Y, M, then the rest.
//
//  Memory: the panel (about 30 bytes per site with five groups) plus one byte
//  per panel site per person; the person's VCF is streamed.
//
//  Nothing here makes medical statements: proportions describe similarity to
//  the reference groups of the panel, and a polygenic score is a weighted
//  sum compared with those groups.
// =====================================================================

use crate::externo::{abre_leitura, natural_cmp, Escritor};
use crate::multimerge_core::multi_merge_sort_stable_by;
use crate::nomes::{canonico, json_texto};
use crate::referencia::Genoma;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

// ---------------------------------------------------------------------
//  READING
// ---------------------------------------------------------------------

/// Line reader over bytes (no UTF-8 check on multi-gigabyte genotype
/// files) that remembers the path and the line number for messages.
struct Linhas {
    caminho: String,
    r: Box<dyn BufRead>,
    buf: Vec<u8>,
    n: usize,
}

impl Linhas {
    fn abre(caminho: &str) -> Result<Linhas, String> {
        Ok(Linhas { caminho: caminho.to_string(), r: abre_leitura(caminho)?, buf: Vec::with_capacity(1 << 16), n: 0 })
    }

    /// Reads the next line into `buf`, without its line ending. `false` at
    /// the end of the file.
    fn le_em(&mut self, buf: &mut Vec<u8>) -> Result<bool, String> {
        buf.clear();
        match self.r.read_until(b'\n', buf) {
            Ok(0) => Ok(false),
            Ok(_) => {
                self.n += 1;
                while matches!(buf.last(), Some(b'\n' | b'\r')) {
                    buf.pop();
                }
                Ok(true)
            }
            Err(e) => Err(format!(
                "Error reading {} after line {}: {e}. The file may be truncated, still being written, or not text (BCF is not read: convert it to VCF first).",
                self.caminho, self.n
            )),
        }
    }

    fn proxima(&mut self) -> Result<bool, String> {
        let mut b = std::mem::take(&mut self.buf);
        let r = self.le_em(&mut b);
        self.buf = b;
        r
    }
}

fn txt(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

/// The start of a line, for messages.
fn trecho(b: &[u8]) -> String {
    let mut t = txt(&b[..b.len().min(60)]);
    if b.len() > 60 {
        t.push_str("...");
    }
    t.replace('\t', " ")
}

fn inteiro(b: &[u8]) -> Option<u64> {
    if b.is_empty() || b.len() > 18 {
        return None;
    }
    let mut v = 0u64;
    for &c in b {
        if !c.is_ascii_digit() {
            return None;
        }
        v = v * 10 + u64::from(c - b'0');
    }
    Some(v)
}

/// A finite number (`parse::<f64>` alone accepts "inf" and "nan").
fn real(b: &[u8]) -> Option<f64> {
    std::str::from_utf8(b).ok()?.trim().parse::<f64>().ok().filter(|x| x.is_finite())
}

fn e_acgt(c: u8) -> bool {
    matches!(c, b'A' | b'C' | b'G' | b'T')
}

/// Non-empty and only A, C, G, T (already in upper case).
fn so_acgt(b: &[u8]) -> bool {
    !b.is_empty() && b.iter().all(|&c| e_acgt(c))
}

fn comp_base(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        x => x,
    }
}

/// A/T or C/G: the two alleles of an SNV read the same on both strands.
fn complementares(a: &[u8], b: &[u8]) -> bool {
    a.len() == 1 && b.len() == 1 && comp_base(a[0]) == b[0]
}

fn passa_filtro(f: &[u8]) -> bool {
    f == b"PASS" || f == b"." || f.is_empty()
}

/// Splits on tabs into reusable (start, end) pairs.
fn divide(l: &[u8], out: &mut Vec<(usize, usize)>) {
    out.clear();
    let mut ini = 0;
    for (i, &c) in l.iter().enumerate() {
        if c == b'\t' {
            out.push((ini, i));
            ini = i + 1;
        }
    }
    out.push((ini, l.len()));
}

/// The eight fixed columns of a VCF record, borrowed from the line, and the
/// rest (FORMAT and the samples, still joined by tabs; empty when absent).
struct Fixas<'a> {
    crom: &'a [u8],
    pos: u64,
    ref_: &'a [u8],
    alt: &'a [u8],
    filtro: &'a [u8],
    info: &'a [u8],
    resto: &'a [u8],
}

/// `None`: not a record (fewer than 8 columns, or POS is not a number).
fn fixas(l: &[u8]) -> Option<Fixas<'_>> {
    let mut c: [&[u8]; 8] = [&[]; 8];
    let mut ini = 0usize;
    let mut resto: &[u8] = &[];
    for (k, campo) in c.iter_mut().enumerate() {
        match l[ini..].iter().position(|&b| b == b'\t') {
            Some(p) => {
                *campo = &l[ini..ini + p];
                ini += p + 1;
                if k == 7 {
                    resto = &l[ini..];
                }
            }
            None => {
                if k < 7 {
                    return None;
                }
                *campo = &l[ini..];
            }
        }
    }
    if c[0].is_empty() {
        return None;
    }
    let pos = inteiro(c[1])?;
    Some(Fixas { crom: c[0], pos, ref_: c[3], alt: c[4], filtro: c[6], info: c[7], resto })
}

// ---------------------------------------------------------------------
//  CHROMOSOMES
// ---------------------------------------------------------------------

/// Order of canonical names: 1..22 by number, X, Y, MT, then the rest in
/// natural order.
fn classe_crom(c: &str) -> (u8, u64) {
    if !c.is_empty() && c.len() <= 6 && c.bytes().all(|b| b.is_ascii_digit()) {
        return (0, c.parse().unwrap_or(u64::MAX));
    }
    match c {
        "X" => (1, 0),
        "Y" => (2, 0),
        "MT" => (3, 0),
        _ => (4, 0),
    }
}

fn compara_crom(a: &str, b: &str) -> Ordering {
    classe_crom(a).cmp(&classe_crom(b)).then_with(|| natural_cmp(a, b))
}

/// Canonical name of the last chromosome seen: inside a chromosome the name
/// is not converted again.
#[derive(Default)]
struct Canon {
    cru: Vec<u8>,
    chave: String,
    vazio: bool,
}

impl Canon {
    fn de(&mut self, cru: &[u8]) -> &str {
        if !self.vazio || cru != self.cru.as_slice() {
            self.cru.clear();
            self.cru.extend_from_slice(cru);
            self.chave = canonico(&String::from_utf8_lossy(cru));
            self.vazio = true;
        }
        &self.chave
    }
}

fn erro_ordem(caminho: &str, linha: usize, o_que: &str) -> String {
    format!(
        "{caminho}, line {linha}: {o_que}. The file must be sorted by position, with each chromosome in a single block. Sort it first -- sort_genomic([\"{caminho}\"], \"sorted.vcf.gz\", chrom_order=\"natural\") -- and use the sorted file."
    )
}

// ---------------------------------------------------------------------
//  OUTPUT FILES
// ---------------------------------------------------------------------

fn mesmo_arquivo(a: &str, b: &str) -> bool {
    a == b || matches!((std::fs::canonicalize(a), std::fs::canonicalize(b)), (Ok(x), Ok(y)) if x == y)
}

/// Refuses outputs that would destroy an input or that are not text tables.
fn confere_saida(saida: &str, entradas: &[&str]) -> Result<(), String> {
    if saida.trim().is_empty() {
        return Err("The output path is empty: give a file name such as panel.tsv.gz.".to_string());
    }
    for e in entradas {
        if mesmo_arquivo(e, saida) {
            return Err(format!(
                "The output `{saida}` is also one of the input files. Writing would destroy the input while it is being read: choose a different output path."
            ));
        }
    }
    let c = saida.to_ascii_lowercase();
    if c.ends_with(".parquet") || c.ends_with(".bam") {
        return Err(format!(
            "`{saida}`: this table is written as tab-separated text. Give it a .tsv, .tsv.gz or .txt name (convert it afterwards if needed)."
        ));
    }
    Ok(())
}

/// Deletes a half-written output when the operation ends in error.
struct Guarda {
    caminho: String,
    armada: bool,
}

impl Drop for Guarda {
    fn drop(&mut self) {
        if self.armada {
            let _ = std::fs::remove_file(&self.caminho);
        }
    }
}

/// A frequency with up to 6 decimals, without trailing zeros.
fn fmt_freq(x: f64) -> String {
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

fn jn(x: f64) -> String {
    if x.is_finite() {
        format!("{x}")
    } else {
        "null".to_string()
    }
}

fn jlista_txt(v: &[String]) -> String {
    format!("[{}]", v.iter().map(|s| json_texto(s)).collect::<Vec<_>>().join(","))
}

fn jlista_num(v: &[f64]) -> String {
    format!("[{}]", v.iter().map(|&x| jn(x)).collect::<Vec<_>>().join(","))
}

fn jlista_us(v: &[usize]) -> String {
    format!("[{}]", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","))
}

// =====================================================================
//  1. PANEL
// =====================================================================

/// Options of the panel builders.
#[derive(Clone, Debug, PartialEq)]
pub struct PanelOptions {
    /// Keep a site when its minor-allele frequency is at least this in at
    /// least one group. 0 keeps every biallelic site.
    pub min_maf: f64,
    /// Keep at most one site per this many bases (the first one that passes
    /// the other filters). 0 keeps them all.
    pub thin_bp: u64,
    /// Also keep biallelic indels (useful for a panel meant for polygenic
    /// scores; ancestry uses only the SNVs).
    pub keep_indels: bool,
    /// Leave out records whose FILTER is not PASS or `.`.
    pub pass_only: bool,
}

impl Default for PanelOptions {
    fn default() -> Self {
        PanelOptions { min_maf: 0.05, thin_bp: 0, keep_indels: false, pass_only: true }
    }
}

/// What a panel builder did.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PanelSummary {
    /// "genotypes" or "frequencies".
    pub source: String,
    pub groups: Vec<String>,
    /// Samples counted in each group (genotype source, first VCF); empty for
    /// the frequency source.
    pub samples_per_group: Vec<usize>,
    /// Sample columns of the (first) VCF.
    pub samples_in_vcf: usize,
    /// Samples of the VCF left out because they were in `exclude`.
    pub samples_excluded: usize,
    /// Samples of the VCF that are not in the table (not counted).
    pub samples_without_group: usize,
    /// Samples of the table that are not in the VCF.
    pub table_samples_not_in_vcf: usize,
    pub records: usize,
    /// Indels (unless kept), MNPs, symbolic and other non-SNV records.
    pub not_snv: usize,
    /// Records with several ALT alleles, and records at a position that has
    /// more than one variant record (a multiallelic site split in lines).
    pub multiallelic: usize,
    /// FILTER not PASS (with `pass_only`).
    pub filtered: usize,
    /// Sites where some group has no called allele (or no frequency).
    pub no_data: usize,
    pub below_maf: usize,
    pub thinned: usize,
    pub written: usize,
    /// Lines that are not records (fewer than 8 columns, POS not a number).
    pub skipped: usize,
    /// Sites written per chromosome, in output order.
    pub sequences: Vec<(String, usize)>,
}

impl PanelSummary {
    pub fn json(&self) -> String {
        let seqs: Vec<String> = self.sequences.iter().map(|(n, q)| format!("[{},{q}]", json_texto(n))).collect();
        format!(
            "{{\"source\":{},\"groups\":{},\"samples_per_group\":{},\"samples_in_vcf\":{},\"samples_excluded\":{},\"samples_without_group\":{},\"table_samples_not_in_vcf\":{},\"records\":{},\"not_snv\":{},\"multiallelic\":{},\"filtered\":{},\"no_data\":{},\"below_maf\":{},\"thinned\":{},\"written\":{},\"skipped\":{},\"sequences\":[{}]}}",
            json_texto(&self.source),
            jlista_txt(&self.groups),
            jlista_us(&self.samples_per_group),
            self.samples_in_vcf,
            self.samples_excluded,
            self.samples_without_group,
            self.table_samples_not_in_vcf,
            self.records,
            self.not_snv,
            self.multiallelic,
            self.filtered,
            self.no_data,
            self.below_maf,
            self.thinned,
            self.written,
            self.skipped,
            seqs.join(",")
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!("Panel of {} group(s) from {}: {} site(s) written of {} record(s)\n", self.groups.len(), self.source, self.written, self.records);
        if !self.samples_per_group.is_empty() {
            let g: Vec<String> = self.groups.iter().zip(&self.samples_per_group).map(|(g, n)| format!("{g} {n}")).collect();
            o.push_str(&format!("  samples per group: {}\n", g.join(", ")));
            o.push_str(&format!(
                "  VCF samples: {} ({} excluded, {} not in the table); {} table sample(s) not in the VCF\n",
                self.samples_in_vcf, self.samples_excluded, self.samples_without_group, self.table_samples_not_in_vcf
            ));
        }
        o.push_str(&format!(
            "  left out: {} not SNV, {} multiallelic, {} filtered, {} without data in some group, {} below the MAF, {} thinned, {} malformed line(s)\n",
            self.not_snv, self.multiallelic, self.filtered, self.no_data, self.below_maf, self.thinned, self.skipped
        ));
        if !self.sequences.is_empty() {
            let l: Vec<String> = self.sequences.iter().take(25).map(|(n, q)| format!("{n} {q}")).collect();
            o.push_str(&format!("  sites per chromosome: {}\n", l.join(", ")));
        }
        o
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Tipo {
    Snv,
    Indel,
    NaoSnv,
    Multi,
}

fn tipo_variante(ref_: &[u8], alt: &[u8]) -> Tipo {
    if alt.contains(&b',') {
        return Tipo::Multi;
    }
    let ok = |s: &[u8]| !s.is_empty() && s.iter().all(|c| e_acgt(c.to_ascii_uppercase()));
    if !ok(ref_) || !ok(alt) {
        return Tipo::NaoSnv;
    }
    if ref_.len() == 1 && alt.len() == 1 {
        if ref_.eq_ignore_ascii_case(alt) {
            Tipo::NaoSnv
        } else {
            Tipo::Snv
        }
    } else if ref_.len() != alt.len() {
        Tipo::Indel
    } else {
        Tipo::NaoSnv
    }
}

/// What a source line gives, decided in parallel.
enum Analise {
    /// Not a record.
    Ruim,
    /// A malformed value that stops the build (message without the line).
    Erro(String),
    Registro { crom: Vec<u8>, pos: u64, v: Veredito },
}

enum Veredito {
    NaoSnv,
    Multi,
    Filtrado,
    Candidato { ref_: Vec<u8>, alt: Vec<u8>, af: Vec<f64>, an: Vec<Option<u64>>, sem_dados: bool },
}

/// Type and FILTER of a record; `None` lets the source compute the values.
fn veredito_inicial(fx: &Fixas<'_>, o: &PanelOptions) -> Option<Veredito> {
    match tipo_variante(fx.ref_, fx.alt) {
        Tipo::Multi => Some(Veredito::Multi),
        Tipo::NaoSnv => Some(Veredito::NaoSnv),
        Tipo::Indel if !o.keep_indels => Some(Veredito::NaoSnv),
        _ if o.pass_only && !passa_filtro(fx.filtro) => Some(Veredito::Filtrado),
        _ => None,
    }
}

/// A source of panel values: genotypes + table, or INFO frequencies.
trait Fonte: Sync {
    /// The `#CHROM` line of each file (sample names) and the INFO IDs its
    /// header declares.
    fn cabecalho(&mut self, caminho: &str, amostras: &[String], infos: &[String]) -> Result<(), String>;
    fn analisa(&self, linha: &[u8], o: &PanelOptions) -> Analise;
}

const SEM_GRUPO: u16 = u16::MAX;

struct FonteGenotipos {
    k: usize,
    grupo_de: HashMap<String, u16>,
    excluir: HashSet<String>,
    /// Group of each sample column of the current file.
    grupo_col: Vec<u16>,
    arquivos: usize,
    resumo: PanelSummary,
    grupos: Vec<String>,
}

/// Counts the called alleles of a GT and how many are ALT (allele 1).
/// Each allele counts on its own (`0/.` is one called REF allele), as AN
/// does in bcftools; an allele index above 1 is not counted.
fn conta_alelos(gt: &[u8], alt: &mut u64, an: &mut u64) {
    if !gt.iter().all(|&c| c.is_ascii_digit() || matches!(c, b'.' | b'/' | b'|')) {
        return;
    }
    for a in gt.split(|&c| c == b'/' || c == b'|') {
        match a {
            b"0" => *an += 1,
            b"1" => {
                *an += 1;
                *alt += 1
            }
            _ => {}
        }
    }
}

impl Fonte for FonteGenotipos {
    fn cabecalho(&mut self, caminho: &str, amostras: &[String], _infos: &[String]) -> Result<(), String> {
        if amostras.is_empty() {
            return Err(format!(
                "{caminho} has no sample columns (a sites-only VCF): a panel from genotypes needs them. For a VCF with per-group frequencies in INFO, use panel_from_frequencies."
            ));
        }
        let mut por_grupo = vec![0usize; self.k];
        let (mut excl, mut sem) = (0usize, 0usize);
        let mut vistos: HashSet<&str> = HashSet::new();
        self.grupo_col = amostras
            .iter()
            .map(|s| {
                vistos.insert(s.as_str());
                if self.excluir.contains(s) {
                    excl += 1;
                    return SEM_GRUPO;
                }
                match self.grupo_de.get(s) {
                    Some(&g) => {
                        por_grupo[g as usize] += 1;
                        g
                    }
                    None => {
                        sem += 1;
                        SEM_GRUPO
                    }
                }
            })
            .collect();
        if let Some(g) = por_grupo.iter().position(|&n| n == 0) {
            return Err(format!(
                "{caminho}: no sample of group `{}` is in this VCF (after the exclusions). Check that the sample names of the table match the #CHROM line, or remove the group from the table.",
                self.grupos[g]
            ));
        }
        if self.arquivos == 0 {
            self.resumo.samples_per_group = por_grupo;
            self.resumo.samples_in_vcf = amostras.len();
            self.resumo.samples_excluded = excl;
            self.resumo.samples_without_group = sem;
            self.resumo.table_samples_not_in_vcf = self.grupo_de.keys().filter(|s| !vistos.contains(s.as_str())).count();
        }
        self.arquivos += 1;
        Ok(())
    }

    fn analisa(&self, l: &[u8], o: &PanelOptions) -> Analise {
        let Some(fx) = fixas(l) else { return Analise::Ruim };
        if fx.resto.is_empty() {
            return Analise::Ruim;
        }
        let crom = fx.crom.to_vec();
        if let Some(v) = veredito_inicial(&fx, o) {
            return Analise::Registro { crom, pos: fx.pos, v };
        }
        let (mut alt, mut an) = (vec![0u64; self.k], vec![0u64; self.k]);
        let mut campos = fx.resto.split(|&b| b == b'\t');
        let formato = campos.next().unwrap_or(&[]);
        if let Some(ig) = formato.split(|&b| b == b':').position(|x| x == b"GT") {
            for (i, campo) in campos.enumerate() {
                let g = match self.grupo_col.get(i) {
                    Some(&g) if g != SEM_GRUPO => g as usize,
                    Some(_) => continue,
                    None => break,
                };
                let gt = if ig == 0 {
                    campo.split(|&b| b == b':').next().unwrap_or(&[])
                } else {
                    match campo.split(|&b| b == b':').nth(ig) {
                        Some(x) => x,
                        None => continue,
                    }
                };
                conta_alelos(gt, &mut alt[g], &mut an[g]);
            }
        }
        let sem_dados = an.contains(&0);
        let af = alt.iter().zip(&an).map(|(&x, &n)| if n > 0 { x as f64 / n as f64 } else { f64::NAN }).collect();
        Analise::Registro {
            crom,
            pos: fx.pos,
            v: Veredito::Candidato {
                ref_: fx.ref_.to_ascii_uppercase(),
                alt: fx.alt.to_ascii_uppercase(),
                af,
                an: an.into_iter().map(Some).collect(),
                sem_dados,
            },
        }
    }
}

struct FonteFrequencias {
    grupos: Vec<String>,
    af: Vec<Vec<u8>>,
    an: Vec<Vec<u8>>,
}

impl Fonte for FonteFrequencias {
    fn cabecalho(&mut self, caminho: &str, _amostras: &[String], infos: &[String]) -> Result<(), String> {
        if infos.is_empty() {
            return Ok(());
        }
        for campo in self.af.iter().chain(&self.an) {
            let c = txt(campo);
            if !infos.contains(&c) {
                let mut lista: Vec<&str> = infos.iter().map(String::as_str).take(30).collect();
                if infos.len() > 30 {
                    lista.push("...");
                }
                return Err(format!(
                    "{caminho}: the INFO field `{c}` is not declared in the header (the ##INFO lines declare: {}). Name the per-group fields exactly as they appear in the file.",
                    lista.join(", ")
                ));
            }
        }
        Ok(())
    }

    fn analisa(&self, l: &[u8], o: &PanelOptions) -> Analise {
        let Some(fx) = fixas(l) else { return Analise::Ruim };
        let crom = fx.crom.to_vec();
        if let Some(v) = veredito_inicial(&fx, o) {
            return Analise::Registro { crom, pos: fx.pos, v };
        }
        let k = self.grupos.len();
        let mut af = vec![f64::NAN; k];
        let mut an: Vec<Option<u64>> = vec![None; k];
        for kv in fx.info.split(|&b| b == b';') {
            let (chave, valor) = match kv.iter().position(|&b| b == b'=') {
                Some(p) => (&kv[..p], &kv[p + 1..]),
                None => (kv, &[][..]),
            };
            for g in 0..k {
                if chave == self.af[g].as_slice() && valor != b"." {
                    match real(valor) {
                        Some(x) if (0.0..=1.0).contains(&x) => af[g] = x,
                        _ => {
                            return Analise::Erro(format!(
                                "the INFO field {} is `{}`, not a frequency between 0 and 1",
                                txt(chave),
                                trecho(valor)
                            ))
                        }
                    }
                }
                if self.an.get(g).map(|x| x.as_slice()) == Some(chave) && valor != b"." {
                    match inteiro(valor) {
                        Some(n) => an[g] = Some(n),
                        None => {
                            return Analise::Erro(format!(
                                "the INFO field {} is `{}`, not a whole number of alleles",
                                txt(chave),
                                trecho(valor)
                            ))
                        }
                    }
                }
            }
        }
        let sem_dados = af.iter().any(|x| x.is_nan()) || an.contains(&Some(0));
        Analise::Registro {
            crom,
            pos: fx.pos,
            v: Veredito::Candidato {
                ref_: fx.ref_.to_ascii_uppercase(),
                alt: fx.alt.to_ascii_uppercase(),
                af,
                an,
                sem_dados,
            },
        }
    }
}

/// Canonical name of the first record of a VCF (`None`: no record).
fn primeiro_cromossomo(caminho: &str) -> Result<Option<String>, String> {
    let mut l = Linhas::abre(caminho)?;
    while l.proxima()? {
        if l.buf.is_empty() || l.buf[0] == b'#' {
            continue;
        }
        if let Some(fx) = fixas(&l.buf) {
            return Ok(Some(canonico(&txt(fx.crom))));
        }
    }
    Ok(None)
}

/// The sequential end of the panel build: order checks, the one-record-per-
/// position rule, MAF, thinning and writing.
struct Montagem {
    w: Option<Escritor>,
    min_maf: f64,
    thin: u64,
    canon: Canon,
    vistos: HashSet<String>,
    crom: Option<String>,
    crom_cru: String,
    pos_atual: Option<u64>,
    grupo: Vec<Veredito>,
    ultimo_mantido: Option<u64>,
    r: PanelSummary,
}

impl Montagem {
    fn registro(&mut self, caminho: &str, linha: usize, crom: &[u8], pos: u64, v: Veredito) -> Result<(), String> {
        self.r.records += 1;
        let chave = self.canon.de(crom).to_string();
        if self.crom.as_deref() != Some(chave.as_str()) {
            self.fecha_posicao()?;
            if let Some(c) = self.crom.take() {
                self.vistos.insert(c);
            }
            if self.vistos.contains(&chave) {
                return Err(erro_ordem(caminho, linha, &format!("chromosome {} appears again after other chromosomes", txt(crom))));
            }
            self.crom = Some(chave);
            self.crom_cru = txt(crom);
            self.pos_atual = None;
            self.ultimo_mantido = None;
        }
        if let Some(p) = self.pos_atual {
            if pos < p {
                return Err(erro_ordem(caminho, linha, &format!("position {pos} comes after position {p} on {}", txt(crom))));
            }
            if pos != p {
                self.fecha_posicao()?;
            }
        }
        self.pos_atual = Some(pos);
        match v {
            Veredito::NaoSnv => self.r.not_snv += 1,
            Veredito::Multi => self.r.multiallelic += 1,
            Veredito::Filtrado => self.r.filtered += 1,
            c @ Veredito::Candidato { .. } => self.grupo.push(c),
        }
        Ok(())
    }

    fn fecha_posicao(&mut self) -> Result<(), String> {
        if self.grupo.len() > 1 {
            self.r.multiallelic += self.grupo.len();
            self.grupo.clear();
            return Ok(());
        }
        let Some(Veredito::Candidato { ref_, alt, af, an, sem_dados }) = self.grupo.pop() else {
            return Ok(());
        };
        let pos = self.pos_atual.unwrap_or(0);
        if sem_dados {
            self.r.no_data += 1;
            return Ok(());
        }
        if !af.iter().any(|&p| p.min(1.0 - p) >= self.min_maf - 1e-12) {
            self.r.below_maf += 1;
            return Ok(());
        }
        if let Some(u) = self.ultimo_mantido {
            if self.thin > 0 && pos < u + self.thin {
                self.r.thinned += 1;
                return Ok(());
            }
        }
        self.ultimo_mantido = Some(pos);
        let mut l = format!("{}\t{pos}\t{}\t{}", self.crom_cru, txt(&ref_), txt(&alt));
        for p in &af {
            l.push('\t');
            l.push_str(&fmt_freq(*p));
        }
        for n in &an {
            l.push('\t');
            match n {
                Some(n) => l.push_str(&n.to_string()),
                None => l.push('.'),
            }
        }
        if let Some(w) = self.w.as_mut() {
            w.linha(&l).map_err(|e| format!("Error writing the panel: {e}"))?;
        }
        self.r.written += 1;
        match self.r.sequences.last_mut() {
            Some((n, q)) if *n == self.crom_cru => *q += 1,
            _ => self.r.sequences.push((self.crom_cru.clone(), 1)),
        }
        Ok(())
    }
}

/// Lines read and analysed per parallel batch.
const LOTE_LINHAS: usize = 2048;
const LOTE_BYTES: usize = 48 << 20;

fn monta_painel<F: Fonte>(
    vcfs: &[String],
    fonte: &mut F,
    grupos: &[String],
    output: &str,
    o: &PanelOptions,
    cabecalho_extra: &[String],
    mut r: PanelSummary,
) -> Result<PanelSummary, String> {
    if vcfs.is_empty() {
        return Err("No input VCF: give at least one file.".to_string());
    }
    if !(0.0..=0.5).contains(&o.min_maf) {
        return Err(format!("min_maf is {}: it must be between 0 and 0.5.", o.min_maf));
    }
    let entradas: Vec<&str> = vcfs.iter().map(String::as_str).collect();
    confere_saida(output, &entradas)?;
    for g in grupos {
        if g.is_empty() || g.contains(['\t', ' ', ',', '\n']) {
            return Err(format!("The group name `{g}` is empty or has a space, comma or tab: use short names such as EUR or AFR."));
        }
    }
    // The files in chromosome order (one file per chromosome, as 1000
    // Genomes ships them, then gives a panel in order 1..22, X, Y, M).
    let mut ordem: Vec<(Option<String>, usize)> = Vec::with_capacity(vcfs.len());
    for (i, v) in vcfs.iter().enumerate() {
        ordem.push((primeiro_cromossomo(v)?, i));
    }
    ordem.sort_by(|a, b| match (&a.0, &b.0) {
        (Some(x), Some(y)) => compara_crom(x, y),
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
        (None, None) => Ordering::Equal,
    });

    let mut guarda = Guarda { caminho: output.to_string(), armada: true };
    let mut w = Escritor::cria(output)?;
    let mut cab = vec![
        "##fileformat=MMGEN-panel-1".to_string(),
        format!("##source={}", r.source),
        format!("##inputs={}", ordem.iter().map(|(_, i)| vcfs[*i].as_str()).collect::<Vec<_>>().join(",")),
        format!("##groups={}", grupos.join(",")),
    ];
    cab.extend(cabecalho_extra.iter().cloned());
    cab.push(format!("##min_maf={}", o.min_maf));
    cab.push(format!("##thin_bp={}", o.thin_bp));
    cab.push(format!("##keep_indels={}", o.keep_indels));
    cab.push(format!("##pass_only={}", o.pass_only));
    let mut colunas = vec!["chrom".to_string(), "pos".to_string(), "ref".to_string(), "alt".to_string()];
    colunas.extend(grupos.iter().map(|g| format!("AF_{g}")));
    colunas.extend(grupos.iter().map(|g| format!("AN_{g}")));
    cab.push(colunas.join("\t"));
    for l in &cab {
        w.linha(l).map_err(|e| format!("Error writing {output}: {e}"))?;
    }
    r.groups = grupos.to_vec();
    let mut m = Montagem {
        w: Some(w),
        min_maf: o.min_maf,
        thin: o.thin_bp,
        canon: Canon::default(),
        vistos: HashSet::new(),
        crom: None,
        crom_cru: String::new(),
        pos_atual: None,
        grupo: Vec::new(),
        ultimo_mantido: None,
        r,
    };

    let mut lote: Vec<Vec<u8>> = Vec::new();
    let mut numeros: Vec<usize> = Vec::new();
    for (_, i) in &ordem {
        let caminho = vcfs[*i].as_str();
        let mut l = Linhas::abre(caminho)?;
        let mut infos: Vec<String> = Vec::new();
        loop {
            if !l.proxima()? {
                return Err(format!("{caminho}: no `#CHROM` header line: is this a VCF file?"));
            }
            let t = &l.buf;
            if t.starts_with(b"##INFO=<ID=") {
                let r = &t[11..];
                let fim = r.iter().position(|&b| b == b',' || b == b'>').unwrap_or(r.len());
                infos.push(txt(&r[..fim]));
                continue;
            }
            if t.starts_with(b"#CHROM") {
                let amostras: Vec<String> = t.split(|&b| b == b'\t').skip(9).map(txt).collect();
                fonte.cabecalho(caminho, &amostras, &infos)?;
                break;
            }
            if t.is_empty() || t[0] == b'#' {
                continue;
            }
            return Err(format!("{caminho}, line {}: a record before the `#CHROM` header line: is this a VCF file?", l.n));
        }
        let fonte_ref: &F = fonte;
        loop {
            let mut usados = 0usize;
            let mut bytes = 0usize;
            numeros.clear();
            while usados < LOTE_LINHAS && bytes < LOTE_BYTES {
                if usados == lote.len() {
                    lote.push(Vec::new());
                }
                if !l.le_em(&mut lote[usados])? {
                    break;
                }
                if lote[usados].is_empty() || lote[usados][0] == b'#' {
                    continue;
                }
                bytes += lote[usados].len();
                numeros.push(l.n);
                usados += 1;
            }
            if usados == 0 {
                break;
            }
            let analises: Vec<Analise> = lote[..usados].par_iter().map(|x| fonte_ref.analisa(x, o)).collect();
            for (a, &n) in analises.into_iter().zip(&numeros) {
                match a {
                    Analise::Ruim => m.r.skipped += 1,
                    Analise::Erro(e) => return Err(format!("{caminho}, line {n}: {e}.")),
                    Analise::Registro { crom, pos, v } => m.registro(caminho, n, &crom, pos, v)?,
                }
            }
        }
    }
    m.fecha_posicao()?;
    let w = m.w.take().expect("the writer is there until the end");
    w.termina().map_err(|e| format!("Error writing {output}: {e}"))?;
    guarda.armada = false;
    Ok(m.r)
}

/// Reads a sample -> group table: a header line, then one sample per line.
/// Tab-separated, or separated by spaces when the header has no tab (as the
/// 1000 Genomes `..._ped_population.txt`). Groups come out in alphabetical
/// order.
fn le_tabela(caminho: &str, col_amostra: &str, col_grupo: &str) -> Result<(Vec<String>, HashMap<String, u16>), String> {
    let mut l = Linhas::abre(caminho)?;
    let mut cab: Option<(bool, usize, usize)> = None;
    let mut pares: Vec<(String, String)> = Vec::new();
    let mut ja: HashMap<String, String> = HashMap::new();
    while l.proxima()? {
        let t = txt(&l.buf);
        if t.trim().is_empty() {
            continue;
        }
        let Some((tab, ia, ig)) = cab else {
            let tab = t.contains('\t');
            let nomes: Vec<&str> = if tab { t.split('\t').map(str::trim).collect() } else { t.split_whitespace().collect() };
            let acha = |n: &str| {
                nomes.iter().position(|x| *x == n).or_else(|| nomes.iter().position(|x| x.trim_start_matches('#').eq_ignore_ascii_case(n)))
            };
            match (acha(col_amostra), acha(col_grupo)) {
                (Some(a), Some(g)) => cab = Some((tab, a, g)),
                _ => {
                    return Err(format!(
                        "{caminho}: the header line has no column `{}` (it has: {}). Name the sample column and the group column as written in the first line.",
                        if acha(col_amostra).is_none() { col_amostra } else { col_grupo },
                        nomes.join(", ")
                    ))
                }
            }
            continue;
        };
        let v: Vec<&str> = if tab { t.split('\t').map(str::trim).collect() } else { t.split_whitespace().collect() };
        let (Some(a), Some(g)) = (v.get(ia), v.get(ig)) else {
            return Err(format!("{caminho}, line {}: fewer columns than the header line.", l.n));
        };
        if a.is_empty() || g.is_empty() {
            continue;
        }
        if let Some(antes) = ja.get(*a) {
            if antes != g {
                return Err(format!("{caminho}, line {}: sample {a} is in group {g} here and in group {antes} before.", l.n));
            }
            continue;
        }
        ja.insert(a.to_string(), g.to_string());
        pares.push((a.to_string(), g.to_string()));
    }
    if pares.is_empty() {
        return Err(format!("{caminho}: no sample in the table (only a header line, or an empty file)."));
    }
    let mut grupos: Vec<String> = pares.iter().map(|(_, g)| g.clone()).collect::<HashSet<_>>().into_iter().collect();
    grupos.sort();
    if grupos.len() > 1000 {
        return Err(format!("{caminho}: {} different groups: is `{col_grupo}` the group column?", grupos.len()));
    }
    let idx: HashMap<&str, u16> = grupos.iter().enumerate().map(|(i, g)| (g.as_str(), i as u16)).collect();
    let mapa = pares.iter().map(|(a, g)| (a.clone(), idx[g.as_str()])).collect();
    Ok((grupos, mapa))
}

/// Builds a frequency panel from the genotypes of a multi-sample VCF (or one
/// VCF per chromosome) and a sample -> group table.
///
/// For each biallelic SNV (and biallelic indel with `keep_indels`) the ALT
/// frequency and the number of called alleles of each group are written.
/// Samples in `exclude` (a hold-out set, relatives) and samples absent from
/// the table are not counted. A position with more than one variant record
/// (a multiallelic site split in several lines) is left out.
pub fn build_panel(
    vcfs: &[String],
    samples_table: &str,
    sample_column: &str,
    group_column: &str,
    exclude: &[String],
    output: &str,
    opts: &PanelOptions,
) -> Result<PanelSummary, String> {
    let (grupos, grupo_de) = le_tabela(samples_table, sample_column, group_column)?;
    let mut entradas = vcfs.to_vec();
    entradas.push(samples_table.to_string());
    confere_saida(output, &entradas.iter().map(String::as_str).collect::<Vec<_>>())?;
    let mut f = FonteGenotipos {
        k: grupos.len(),
        grupo_de,
        excluir: exclude.iter().cloned().collect(),
        grupo_col: Vec::new(),
        arquivos: 0,
        resumo: PanelSummary::default(),
        grupos: grupos.clone(),
    };
    let r = PanelSummary { source: "genotypes".to_string(), ..Default::default() };
    // The sample counts are known only after the first header: they go in
    // the summary, and the table and exclusions in the header of the panel.
    let extra = vec![
        format!("##samples_table={samples_table} (sample column {sample_column}, group column {group_column})"),
        format!("##excluded_samples={}", exclude.len()),
    ];
    let mut r = monta_painel(vcfs, &mut f, &grupos, output, opts, &extra, r)?;
    r.samples_per_group = f.resumo.samples_per_group;
    r.samples_in_vcf = f.resumo.samples_in_vcf;
    r.samples_excluded = f.resumo.samples_excluded;
    r.samples_without_group = f.resumo.samples_without_group;
    r.table_samples_not_in_vcf = f.resumo.table_samples_not_in_vcf;
    Ok(r)
}

/// Builds a frequency panel from a sites VCF whose INFO holds the frequency
/// of each group (`groups[i]` is read from the INFO field `af_fields[i]`,
/// e.g. EUR from EUR_AF). `an_fields`, when not empty, names the INFO field
/// with the number of called alleles of each group; without it the AN
/// columns are written as `.`.
pub fn panel_from_frequencies(
    vcfs: &[String],
    groups: &[String],
    af_fields: &[String],
    an_fields: &[String],
    output: &str,
    opts: &PanelOptions,
) -> Result<PanelSummary, String> {
    if groups.is_empty() || groups.len() != af_fields.len() {
        return Err(format!(
            "Give one INFO frequency field per group: {} group(s) and {} field(s).",
            groups.len(),
            af_fields.len()
        ));
    }
    if !an_fields.is_empty() && an_fields.len() != groups.len() {
        return Err(format!(
            "Give one INFO allele-number field per group, or none: {} group(s) and {} field(s).",
            groups.len(),
            an_fields.len()
        ));
    }
    let mut vistos = HashSet::new();
    for g in groups {
        if !vistos.insert(g) {
            return Err(format!("The group `{g}` is named twice."));
        }
    }
    let mut f = FonteFrequencias {
        grupos: groups.to_vec(),
        af: af_fields.iter().map(|s| s.as_bytes().to_vec()).collect(),
        an: an_fields.iter().map(|s| s.as_bytes().to_vec()).collect(),
    };
    let extra = vec![format!(
        "##info_fields={}",
        groups.iter().zip(af_fields).map(|(g, a)| format!("{g}:{a}")).collect::<Vec<_>>().join(",")
    )];
    let r = PanelSummary { source: "frequencies".to_string(), ..Default::default() };
    monta_painel(vcfs, &mut f, groups, output, opts, &extra, r)
}

// ---------------------------------------------------------------------
//  READING A PANEL
// ---------------------------------------------------------------------

struct ColunasPainel {
    crom: usize,
    pos: usize,
    ref_: usize,
    alt: usize,
    af: Vec<usize>,
    an: Vec<Option<usize>>,
    grupos: Vec<String>,
}

fn colunas_painel(caminho: &str, linha: &[u8]) -> Result<ColunasPainel, String> {
    let nomes: Vec<String> = linha.split(|&b| b == b'\t').map(|c| txt(c).trim().trim_start_matches('#').to_string()).collect();
    let acha = |n: &str| nomes.iter().position(|x| x.eq_ignore_ascii_case(n));
    let (Some(crom), Some(pos), Some(ref_), Some(alt)) = (acha("chrom"), acha("pos"), acha("ref"), acha("alt")) else {
        return Err(format!(
            "{caminho}: the header line must name the columns chrom, pos, ref, alt and one AF_<group> per group (it is: `{}`). Is this a panel written by build_panel?",
            trecho(linha)
        ));
    };
    let mut af = Vec::new();
    let mut grupos = Vec::new();
    for (i, n) in nomes.iter().enumerate() {
        if let Some(g) = n.strip_prefix("AF_") {
            if g.is_empty() {
                return Err(format!("{caminho}: a column is named `AF_` without a group name."));
            }
            af.push(i);
            grupos.push(g.to_string());
        }
    }
    if af.is_empty() {
        return Err(format!("{caminho}: no AF_<group> column in the header line `{}`.", trecho(linha)));
    }
    let an = grupos.iter().map(|g| acha(&format!("AN_{g}"))).collect();
    Ok(ColunasPainel { crom, pos, ref_, alt, af, an, grupos })
}

/// Reads the header of a panel and returns the reader at the first row.
fn abre_painel(caminho: &str) -> Result<(Linhas, ColunasPainel), String> {
    let mut l = Linhas::abre(caminho)?;
    loop {
        if !l.proxima()? {
            return Err(format!("{caminho}: empty panel (no header line)."));
        }
        if l.buf.is_empty() || l.buf.starts_with(b"##") {
            continue;
        }
        let c = colunas_painel(caminho, &l.buf)?;
        return Ok((l, c));
    }
}

/// One row of a panel, as written.
struct LinhaPainel<'a> {
    crom: &'a [u8],
    pos: u64,
    ref_: Vec<u8>,
    alt: Vec<u8>,
    af: Vec<f64>,
    an: Vec<Option<u64>>,
}

/// `Ok(None)`: a group has no frequency (`.`): the row is skipped.
fn linha_painel<'a>(caminho: &str, n: usize, l: &'a [u8], sp: &[(usize, usize)], c: &ColunasPainel) -> Result<Option<LinhaPainel<'a>>, String> {
    let max = [c.crom, c.pos, c.ref_, c.alt].into_iter().chain(c.af.iter().copied()).chain(c.an.iter().flatten().copied()).max().unwrap_or(0);
    if sp.len() <= max {
        return Err(format!("{caminho}, line {n}: {} column(s), fewer than the header line names.", sp.len()));
    }
    let campo = move |i: usize| -> &'a [u8] { &l[sp[i].0..sp[i].1] };
    let pos = inteiro(campo(c.pos)).ok_or_else(|| format!("{caminho}, line {n}: pos `{}` is not a whole number.", trecho(campo(c.pos))))?;
    let mut af = Vec::with_capacity(c.af.len());
    for (g, &i) in c.af.iter().enumerate() {
        let v = campo(i);
        if v == b"." || v.eq_ignore_ascii_case(b"NA") || v.is_empty() {
            return Ok(None);
        }
        match real(v) {
            Some(x) if (0.0..=1.0).contains(&x) => af.push(x),
            _ => {
                return Err(format!(
                    "{caminho}, line {n}: AF_{} is `{}`, not a frequency between 0 and 1.",
                    c.grupos[g],
                    trecho(v)
                ))
            }
        }
    }
    let mut an = Vec::with_capacity(c.an.len());
    for (g, i) in c.an.iter().enumerate() {
        an.push(match i {
            None => None,
            Some(i) => {
                let v = campo(*i);
                if v == b"." || v.is_empty() || v.eq_ignore_ascii_case(b"NA") {
                    None
                } else {
                    Some(inteiro(v).ok_or_else(|| {
                        format!("{caminho}, line {n}: AN_{} is `{}`, not a whole number.", c.grupos[g], trecho(v))
                    })?)
                }
            }
        });
    }
    Ok(Some(LinhaPainel {
        crom: campo(c.crom),
        pos,
        ref_: campo(c.ref_).to_ascii_uppercase(),
        alt: campo(c.alt).to_ascii_uppercase(),
        af,
        an,
    }))
}

/// Frequency used by the likelihood: never 0 or 1, or a person carrying an
/// allele never seen in a group would rule the group out entirely. With the
/// allele count known, (x + 1/2) / (n + 1) (x = ALT alleles seen); without
/// it, the frequency is kept within [0.001, 0.999].
fn freq_ajustada(p: f64, an: Option<u64>) -> f64 {
    match an {
        Some(n) if n > 0 => (p * n as f64 + 0.5) / (n as f64 + 1.0),
        _ => p.clamp(1e-3, 1.0 - 1e-3),
    }
}

/// The SNV sites of a panel, in memory, for ancestry.
struct PainelSnv {
    grupos: Vec<String>,
    /// Canonical chromosome names, in order (1..22, X, Y, MT, others).
    cromos: Vec<String>,
    /// Offset of each chromosome in the arrays below; one more at the end.
    inicio: Vec<usize>,
    pos: Vec<u64>,
    /// REF and ALT.
    ra: Vec<[u8; 2]>,
    /// K adjusted frequencies per site.
    p: Vec<f32>,
    linhas: usize,
    nao_snv: usize,
    sem_dados: usize,
    duplicados: usize,
    /// Every frequency has at most 2 decimals.
    arredondado: bool,
}

impl PainelSnv {
    fn le(caminho: &str) -> Result<PainelSnv, String> {
        let (mut l, c) = abre_painel(caminho)?;
        let k = c.grupos.len();
        struct Tmp {
            pos: Vec<u64>,
            ra: Vec<[u8; 2]>,
            p: Vec<f32>,
        }
        let mut indice: HashMap<String, usize> = HashMap::new();
        let mut tmp: Vec<(String, Tmp)> = Vec::new();
        let mut canon = Canon::default();
        let mut atual = usize::MAX;
        let mut chave_atual = String::new();
        let mut sp = Vec::new();
        let (mut linhas, mut nao_snv, mut sem_dados) = (0usize, 0usize, 0usize);
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
            let chave = canon.de(r.crom);
            if atual == usize::MAX || chave != chave_atual {
                chave_atual = chave.to_string();
                atual = match indice.get(&chave_atual) {
                    Some(&i) => i,
                    None => {
                        indice.insert(chave_atual.clone(), tmp.len());
                        tmp.push((chave_atual.clone(), Tmp { pos: Vec::new(), ra: Vec::new(), p: Vec::new() }));
                        tmp.len() - 1
                    }
                };
            }
            arredondado &= r.af.iter().all(|&x| duas_casas(x));
            let t = &mut tmp[atual].1;
            t.pos.push(r.pos);
            t.ra.push([r.ref_[0], r.alt[0]]);
            for g in 0..k {
                t.p.push(freq_ajustada(r.af[g], r.an[g]) as f32);
            }
        }
        tmp.sort_by(|a, b| compara_crom(&a.0, &b.0));
        let total: usize = tmp.iter().map(|(_, t)| t.pos.len()).sum();
        let mut p = PainelSnv {
            grupos: c.grupos,
            cromos: Vec::new(),
            inicio: vec![0],
            pos: Vec::with_capacity(total),
            ra: Vec::with_capacity(total),
            p: Vec::with_capacity(total * k),
            linhas,
            nao_snv,
            sem_dados,
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
                // The same site twice: the first one stays.
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
            p.inicio.push(p.pos.len());
        }
        Ok(p)
    }

    fn n(&self) -> usize {
        self.pos.len()
    }
}

// ---------------------------------------------------------------------
//  CALLABLE REGIONS
// ---------------------------------------------------------------------

/// A BED of callable regions, merged, by canonical chromosome name.
struct Regioes {
    mapa: HashMap<String, Vec<(u64, u64)>>,
}

impl Regioes {
    fn le(caminho: &str) -> Result<Regioes, String> {
        let mut l = Linhas::abre(caminho)?;
        let mut mapa: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        let mut canon = Canon::default();
        while l.proxima()? {
            let t = &l.buf;
            if t.is_empty() || t[0] == b'#' || t.starts_with(b"track") || t.starts_with(b"browser") {
                continue;
            }
            let mut c = t.split(|&b| b == b'\t');
            let crom = c.next().unwrap_or(&[]);
            let (Some(a), Some(b)) = (c.next().and_then(inteiro), c.next().and_then(inteiro)) else {
                return Err(format!(
                    "{caminho}, line {}: not a BED line (chromosome, start and end separated by tabs): `{}`. Give the callable regions as a BED file.",
                    l.n,
                    trecho(t)
                ));
            };
            if b < a {
                return Err(format!("{caminho}, line {}: the end ({b}) is before the start ({a}).", l.n));
            }
            mapa.entry(canon.de(crom).to_string()).or_default().push((a, b));
        }
        if mapa.is_empty() {
            return Err(format!("{caminho}: no region in this BED file."));
        }
        for v in mapa.values_mut() {
            multi_merge_sort_stable_by(v, |a, b| a.cmp(b));
            let mut out: Vec<(u64, u64)> = Vec::with_capacity(v.len());
            for &(a, b) in v.iter() {
                if b <= a {
                    continue;
                }
                if let Some(u) = out.last_mut() {
                    if a <= u.1 {
                        u.1 = u.1.max(b);
                        continue;
                    }
                }
                out.push((a, b));
            }
            *v = out;
        }
        Ok(Regioes { mapa })
    }

    /// Is the 1-based position inside a region?
    fn contem(&self, crom: &str, pos1: u64) -> bool {
        let Some(v) = self.mapa.get(crom) else { return false };
        if pos1 == 0 {
            return false;
        }
        let p = pos1 - 1;
        let i = v.partition_point(|&(_, e)| e <= p);
        i < v.len() && v[i].0 <= p
    }
}

// ---------------------------------------------------------------------
//  THE PERSON'S GENOTYPES: STREAMING JOIN
// ---------------------------------------------------------------------

/// One person: a VCF, a sample in it (`None`: the only sample), and
/// optionally a BED of callable regions (inside it, a site absent from the
/// VCF is homozygous reference; outside it, nothing is used).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Person {
    pub vcf: String,
    pub sample: Option<String>,
    pub callable: Option<String>,
}

/// A genotype: allele indices; `n` = 0 missing, 1 haploid, 2 diploid.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Gt {
    a: [u8; 2],
    n: u8,
}

const SEM_GT: Gt = Gt { a: [0, 0], n: 0 };

impl Gt {
    fn alelos(&self) -> &[u8] {
        &self.a[..self.n as usize]
    }
}

/// The GT subfield `ig` of a sample column. A partly missing call (`0/.`) is
/// missing; more than two alleles is not read (missing).
fn le_gt(campo: &[u8], ig: usize) -> Gt {
    let gt = match campo.split(|&b| b == b':').nth(ig) {
        Some(x) => x,
        None => return SEM_GT,
    };
    let mut g = SEM_GT;
    let mut v: u32 = 0;
    let (mut tem, mut falta) = (false, false);
    for &c in gt.iter().chain(std::iter::once(&b'/')) {
        match c {
            b'0'..=b'9' => {
                v = v * 10 + u32::from(c - b'0');
                tem = true;
                if v > 254 {
                    return SEM_GT;
                }
            }
            b'/' | b'|' => {
                if falta || !tem || g.n == 2 {
                    return SEM_GT;
                }
                g.a[g.n as usize] = v as u8;
                g.n += 1;
                v = 0;
                tem = false;
            }
            b'.' => falta = true,
            _ => return SEM_GT,
        }
    }
    g
}

/// A record of the person's VCF at a position of interest.
struct RegPessoa {
    ref_: Vec<u8>,
    alts: Vec<Vec<u8>>,
    filtro_ok: bool,
    /// One per selected sample.
    gts: Vec<Gt>,
}

impl RegPessoa {
    fn alelo(&self, i: u8) -> Option<&[u8]> {
        if i == 0 {
            Some(&self.ref_)
        } else {
            self.alts.get(i as usize - 1).map(Vec::as_slice)
        }
    }

    fn indice(&self, a: &[u8]) -> Option<u8> {
        if self.ref_ == a {
            Some(0)
        } else {
            self.alts.iter().position(|x| x == a).map(|i| (i + 1) as u8)
        }
    }
}

/// What the person's file says at a target position, for one sample.
enum Obs<'a> {
    /// The records at this position (the sample's GT is `gts[sample]`).
    Registros(&'a [RegPessoa]),
    /// No record, inside the callable regions (or `absent_is_ref`):
    /// homozygous reference.
    Ausente,
    /// No record, outside the callable regions.
    NaoChamavel,
    /// No record, but a longer variant of the sample (a deletion) starting
    /// before this position covers it.
    SobIndel,
}

fn resolve_amostras(vcf: &str, nomes: &[String], pedidas: &[Option<String>]) -> Result<Vec<usize>, String> {
    let mut idx: HashMap<&str, usize> = HashMap::new();
    for (i, n) in nomes.iter().enumerate() {
        idx.entry(n.as_str()).or_insert(i);
    }
    let exemplo = |n: usize| -> String {
        let mut v: Vec<&str> = nomes.iter().take(n).map(String::as_str).collect();
        if nomes.len() > n {
            v.push("...");
        }
        v.join(", ")
    };
    pedidas
        .iter()
        .map(|p| match p {
            Some(n) => idx.get(n.as_str()).copied().ok_or_else(|| {
                format!(
                    "{vcf}: no sample named `{n}` in the #CHROM line (it has {} sample(s): {}). Check the spelling.",
                    nomes.len(),
                    exemplo(5)
                )
            }),
            None => match nomes.len() {
                0 => Err(format!("{vcf} has no sample columns (a sites-only VCF): a person needs genotypes.")),
                1 => Ok(0),
                n => Err(format!("{vcf} has {n} samples ({}): name the sample to use.", exemplo(3))),
            },
        })
        .collect()
}

/// GTs of the selected samples; `ordem` is (column, selected index) sorted
/// by column.
fn le_gts(resto: &[u8], ordem: &[(usize, usize)], ns: usize) -> Vec<Gt> {
    let mut gts = vec![SEM_GT; ns];
    let mut campos = resto.split(|&b| b == b'\t');
    let Some(formato) = campos.next() else { return gts };
    let Some(ig) = formato.split(|&b| b == b':').position(|x| x == b"GT") else { return gts };
    let mut k = 0;
    for (c, campo) in campos.enumerate() {
        while k < ordem.len() && ordem[k].0 == c {
            gts[ordem[k].1] = le_gt(campo, ig);
            k += 1;
        }
        if k >= ordem.len() {
            break;
        }
    }
    gts
}

type Retorno<'f> = &'f mut dyn FnMut(&str, usize, usize, &Obs<'_>);

struct Juncao<'a, 'f> {
    alvos: &'a HashMap<String, &'a [u64]>,
    regioes: &'a [Option<&'a Regioes>],
    absent_is_ref: bool,
    f: Retorno<'f>,
    chave: String,
    posicoes: &'a [u64],
    cursor: usize,
    /// Last base covered by a longer allele of each sample.
    del_fim: Vec<u64>,
    grupo: Vec<RegPessoa>,
    pos_grupo: u64,
}

impl Juncao<'_, '_> {
    fn inicia(&mut self, chave: &str) {
        self.chave = chave.to_string();
        self.posicoes = self.alvos.get(chave).copied().unwrap_or(&[]);
        self.cursor = 0;
        self.del_fim.iter_mut().for_each(|x| *x = 0);
        self.grupo.clear();
    }

    /// Targets before `limite` have no record: absent.
    fn ausentes_ate(&mut self, limite: u64, presente: bool) {
        while self.cursor < self.posicoes.len() && self.posicoes[self.cursor] < limite {
            let p = self.posicoes[self.cursor];
            for s in 0..self.del_fim.len() {
                let obs = if self.del_fim[s] >= p {
                    Obs::SobIndel
                } else if match self.regioes[s] {
                    Some(r) => r.contem(&self.chave, p),
                    None => self.absent_is_ref && presente,
                } {
                    Obs::Ausente
                } else {
                    Obs::NaoChamavel
                };
                (self.f)(&self.chave, self.cursor, s, &obs);
            }
            self.cursor += 1;
        }
    }

    /// Hands the records of the current position to its targets, then
    /// notes the bases the longer alleles cover.
    fn fecha_grupo(&mut self) {
        if self.grupo.is_empty() {
            return;
        }
        let p = self.pos_grupo;
        while self.cursor < self.posicoes.len() && self.posicoes[self.cursor] == p {
            for s in 0..self.del_fim.len() {
                // With a BED, only its regions are trusted: a record outside
                // them (GIAB files have a few) is not used.
                let fora = matches!(self.regioes[s], Some(r) if !r.contem(&self.chave, p));
                let obs = if fora { Obs::NaoChamavel } else { Obs::Registros(&self.grupo) };
                (self.f)(&self.chave, self.cursor, s, &obs);
            }
            self.cursor += 1;
        }
        for r in &self.grupo {
            marca_cobertura(&mut self.del_fim, p, r.ref_.len(), &r.gts);
        }
        self.grupo.clear();
    }
}

/// A record whose REF spans several bases covers the following positions
/// for every sample that carries a non-REF allele (or has no call).
fn marca_cobertura(del_fim: &mut [u64], pos: u64, len_ref: usize, gts: &[Gt]) {
    if len_ref <= 1 {
        return;
    }
    let fim = pos + len_ref as u64 - 1;
    for (s, g) in gts.iter().enumerate() {
        if g.n == 0 || g.alelos().iter().any(|&a| a != 0) {
            del_fim[s] = del_fim[s].max(fim);
        }
    }
}

struct ResumoJuncao {
    nomes: Vec<String>,
}

/// Streams `vcf` against the targets (sorted positions per canonical
/// chromosome) and calls `f(chromosome, target index, sample, observation)`
/// exactly once per target and selected sample.
fn junta(
    vcf: &str,
    pedidas: &[Option<String>],
    regioes: &[Option<&Regioes>],
    absent_is_ref: bool,
    alvos: &HashMap<String, &[u64]>,
    f: Retorno<'_>,
) -> Result<ResumoJuncao, String> {
    let mut l = Linhas::abre(vcf)?;
    let nomes_vcf: Vec<String> = loop {
        if !l.proxima()? {
            return Err(format!("{vcf}: no `#CHROM` header line: is this a VCF file?"));
        }
        let t = &l.buf;
        if t.is_empty() || t.starts_with(b"##") {
            continue;
        }
        if t.starts_with(b"#CHROM") {
            break t.split(|&b| b == b'\t').skip(9).map(txt).collect();
        }
        if t[0] == b'#' {
            continue;
        }
        return Err(format!("{vcf}, line {}: a record before the `#CHROM` header line: is this a VCF file?", l.n));
    };
    let colunas = resolve_amostras(vcf, &nomes_vcf, pedidas)?;
    let ns = colunas.len();
    let mut ordem: Vec<(usize, usize)> = colunas.iter().enumerate().map(|(s, &c)| (c, s)).collect();
    multi_merge_sort_stable_by(&mut ordem, |a, b| a.0.cmp(&b.0));
    let resumo = ResumoJuncao { nomes: colunas.iter().map(|&c| nomes_vcf[c].clone()).collect() };

    let mut j = Juncao {
        alvos,
        regioes,
        absent_is_ref,
        f,
        chave: String::new(),
        posicoes: &[],
        cursor: 0,
        del_fim: vec![0; ns],
        grupo: Vec::new(),
        pos_grupo: 0,
    };
    let mut vistos: HashSet<String> = HashSet::new();
    let mut canon = Canon::default();
    let mut tem_crom = false;
    let mut pos_atual: Option<u64> = None;
    let mut grupo_alvo = false;
    while l.proxima()? {
        if l.buf.is_empty() || l.buf[0] == b'#' {
            continue;
        }
        let Some(fx) = fixas(&l.buf) else { continue };
        let chave = canon.de(fx.crom);
        if !tem_crom || chave != j.chave {
            if tem_crom {
                j.fecha_grupo();
                j.ausentes_ate(u64::MAX, true);
                vistos.insert(j.chave.clone());
            }
            if vistos.contains(chave) {
                return Err(erro_ordem(vcf, l.n, &format!("chromosome {} appears again after other chromosomes", txt(fx.crom))));
            }
            j.inicia(chave);
            tem_crom = true;
            pos_atual = None;
        }
        if let Some(p) = pos_atual {
            if fx.pos < p {
                return Err(erro_ordem(vcf, l.n, &format!("position {} comes after position {p} on {}", fx.pos, txt(fx.crom))));
            }
        }
        if pos_atual != Some(fx.pos) {
            j.fecha_grupo();
            j.ausentes_ate(fx.pos, true);
            grupo_alvo = j.cursor < j.posicoes.len() && j.posicoes[j.cursor] == fx.pos;
            j.pos_grupo = fx.pos;
            pos_atual = Some(fx.pos);
        }
        if grupo_alvo {
            let alts = if fx.alt == b"." { Vec::new() } else { fx.alt.split(|&b| b == b',').map(|a| a.to_ascii_uppercase()).collect() };
            j.grupo.push(RegPessoa {
                ref_: fx.ref_.to_ascii_uppercase(),
                alts,
                filtro_ok: passa_filtro(fx.filtro),
                gts: le_gts(fx.resto, &ordem, ns),
            });
        } else if fx.ref_.len() > 1 && j.cursor < j.posicoes.len() {
            // Only targets still ahead on this chromosome can be covered.
            let gts = le_gts(fx.resto, &ordem, ns);
            marca_cobertura(&mut j.del_fim, fx.pos, fx.ref_.len(), &gts);
        }
    }
    if tem_crom {
        j.fecha_grupo();
        j.ausentes_ate(u64::MAX, true);
        vistos.insert(j.chave.clone());
    }
    let mut faltam: Vec<&String> = alvos.keys().filter(|c| !vistos.contains(c.as_str())).collect();
    faltam.sort_by(|a, b| compara_crom(a, b));
    for c in faltam {
        j.inicia(c);
        j.ausentes_ate(u64::MAX, false);
    }
    Ok(resumo)
}

/// Persons grouped by VCF (each VCF is read once), in order of appearance.
fn por_vcf(persons: &[Person]) -> Vec<(String, Vec<usize>)> {
    let mut out: Vec<(String, Vec<usize>)> = Vec::new();
    for (i, p) in persons.iter().enumerate() {
        match out.iter_mut().find(|(v, _)| *v == p.vcf) {
            Some((_, l)) => l.push(i),
            None => out.push((p.vcf.clone(), vec![i])),
        }
    }
    out
}

fn le_regioes_das_pessoas(persons: &[Person]) -> Result<HashMap<String, Regioes>, String> {
    let mut m = HashMap::new();
    for p in persons {
        if let Some(b) = &p.callable {
            if !m.contains_key(b) {
                m.insert(b.clone(), Regioes::le(b)?);
            }
        }
    }
    Ok(m)
}

// =====================================================================
//  2. ANCESTRY
// =====================================================================

#[derive(Clone, Debug, PartialEq)]
pub struct AncestryOptions {
    /// Without a callable BED: a panel site absent from the VCF is
    /// homozygous reference (on chromosomes the VCF has records for).
    pub absent_is_ref: bool,
    /// A record whose FILTER is not PASS or `.` makes the site missing.
    pub pass_only: bool,
    /// Bootstrap replicates for the standard errors (0: none).
    pub bootstrap: usize,
    /// Block length of the block bootstrap, in bases.
    pub block_bp: u64,
    pub seed: u64,
    /// Maximum EM updates of the warm-up (at most 30 are used) and Newton
    /// steps of the finish (at most 200).
    pub max_iterations: usize,
    /// Stop when a Newton step would raise the log-likelihood by less than
    /// this.
    pub tolerance: f64,
}

impl Default for AncestryOptions {
    fn default() -> Self {
        AncestryOptions {
            absent_is_ref: true,
            pass_only: true,
            bootstrap: 100,
            block_bp: 5_000_000,
            seed: 1,
            max_iterations: 10_000,
            tolerance: 1e-7,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct AncestryResult {
    pub vcf: String,
    pub sample: String,
    pub groups: Vec<String>,
    /// Proportion of each group (sums to 1).
    pub proportions: Vec<f64>,
    /// Block-bootstrap standard errors (zeros without bootstrap).
    pub standard_errors: Vec<f64>,
    /// SNV sites of the panel.
    pub panel_sites: usize,
    /// Sites with a genotype: genotyped + assumed_ref.
    pub sites_used: usize,
    pub sites_genotyped: usize,
    /// Absent from the VCF inside the callable regions: 0/0.
    pub sites_assumed_ref: usize,
    /// Outside the callable regions (records there are not used either), or
    /// absent from a VCF without BED on a chromosome it does not have.
    pub sites_not_callable: usize,
    /// Absent, but inside a deletion (or other longer allele) of the person.
    pub sites_in_deletion: usize,
    /// The person has another ALT, or another kind of variant, there.
    pub sites_other_allele: usize,
    /// The VCF's REF differs from the panel's REF (a different genome build?).
    pub sites_ref_mismatch: usize,
    /// GT missing (`./.`).
    pub sites_missing_genotype: usize,
    /// FILTER not PASS.
    pub sites_filtered: usize,
    /// At the proportions (binomial, with its constant).
    pub log_likelihood: f64,
    pub iterations: usize,
    pub converged: bool,
    /// The group that alone explains the genotypes best.
    pub best_group: String,
    /// Log-likelihood of each group taken alone (proportion 1).
    pub single_group_log_likelihoods: Vec<f64>,
    /// Blocks with data, and bootstrap replicates done.
    pub blocks: usize,
    pub bootstrap_replicates: usize,
    pub warnings: Vec<String>,
}

impl AncestryResult {
    pub fn json(&self) -> String {
        format!(
            "{{\"vcf\":{},\"sample\":{},\"groups\":{},\"proportions\":{},\"standard_errors\":{},\"panel_sites\":{},\"sites_used\":{},\"sites_genotyped\":{},\"sites_assumed_ref\":{},\"sites_not_callable\":{},\"sites_in_deletion\":{},\"sites_other_allele\":{},\"sites_ref_mismatch\":{},\"sites_missing_genotype\":{},\"sites_filtered\":{},\"log_likelihood\":{},\"iterations\":{},\"converged\":{},\"best_group\":{},\"single_group_log_likelihoods\":{},\"blocks\":{},\"bootstrap_replicates\":{},\"warnings\":{}}}",
            json_texto(&self.vcf),
            json_texto(&self.sample),
            jlista_txt(&self.groups),
            jlista_num(&self.proportions),
            jlista_num(&self.standard_errors),
            self.panel_sites,
            self.sites_used,
            self.sites_genotyped,
            self.sites_assumed_ref,
            self.sites_not_callable,
            self.sites_in_deletion,
            self.sites_other_allele,
            self.sites_ref_mismatch,
            self.sites_missing_genotype,
            self.sites_filtered,
            jn(self.log_likelihood),
            self.iterations,
            self.converged,
            json_texto(&self.best_group),
            jlista_num(&self.single_group_log_likelihoods),
            self.blocks,
            self.bootstrap_replicates,
            jlista_txt(&self.warnings)
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!(
            "{} ({}): {} of {} panel sites used ({} genotyped, {} assumed homozygous reference)\n",
            self.sample, self.vcf, self.sites_used, self.panel_sites, self.sites_genotyped, self.sites_assumed_ref
        );
        for (i, g) in self.groups.iter().enumerate() {
            let se = if self.bootstrap_replicates > 0 { format!(" +/- {:.4}", self.standard_errors[i]) } else { String::new() };
            o.push_str(&format!("  {g:<10} {:.4}{se}\n", self.proportions[i]));
        }
        o.push_str(&format!(
            "  best single group: {}; log-likelihood {:.2} after {} update(s) (EM, then Newton){}\n",
            self.best_group,
            self.log_likelihood,
            self.iterations,
            if self.converged { "" } else { " (not converged)" }
        ));
        let fora = self.sites_not_callable + self.sites_in_deletion + self.sites_other_allele + self.sites_ref_mismatch + self.sites_missing_genotype + self.sites_filtered;
        if fora > 0 {
            o.push_str(&format!(
                "  not used: {} not callable, {} inside a deletion, {} other allele, {} REF differs from the panel, {} missing GT, {} filtered\n",
                self.sites_not_callable, self.sites_in_deletion, self.sites_other_allele, self.sites_ref_mismatch, self.sites_missing_genotype, self.sites_filtered
            ));
        }
        if self.bootstrap_replicates > 0 {
            o.push_str(&format!("  standard errors: {} bootstrap replicates over {} blocks\n", self.bootstrap_replicates, self.blocks));
        }
        for w in &self.warnings {
            o.push_str(&format!("  note: {w}\n"));
        }
        o.push_str("  The proportions measure similarity to the reference groups of the panel, under the admixture model.\n");
        o
    }
}

/// Genotype codes: 0..=2 ALT alleles of a diploid call; 3 and 4: haploid REF
/// and ALT; FALTA: missing.
const FALTA: u8 = u8::MAX;
/// (ALT alleles, called alleles, ln of the binomial coefficient).
const DECODIFICA: [(f64, f64, f64); 5] =
    [(0.0, 2.0, 0.0), (1.0, 2.0, std::f64::consts::LN_2), (2.0, 2.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0)];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Motivo {
    Genotipado,
    Assumido,
    NaoChamavel,
    SobIndel,
    OutroAlelo,
    RefDiferente,
    SemGt,
    Filtrado,
}

/// The genotype code of a panel SNV (REF `pr`, ALT `pa`) for one sample.
fn classifica_snv(obs: &Obs<'_>, s: usize, pr: u8, pa: u8, pass_only: bool) -> (u8, Motivo) {
    match obs {
        Obs::Ausente => (0, Motivo::Assumido),
        Obs::NaoChamavel => (FALTA, Motivo::NaoChamavel),
        Obs::SobIndel => (FALTA, Motivo::SobIndel),
        Obs::Registros(regs) => {
            if pass_only && regs.iter().any(|r| !r.filtro_ok) {
                return (FALTA, Motivo::Filtrado);
            }
            let com_ref = |r: &&RegPessoa| r.ref_.len() == 1 && r.ref_[0] == pr;
            let r = regs
                .iter()
                .filter(com_ref)
                .find(|r| r.alts.iter().any(|a| a.len() == 1 && a[0] == pa))
                .or_else(|| regs.iter().find(com_ref));
            let Some(r) = r else {
                return if regs.iter().any(|r| r.ref_.len() == 1 && r.ref_[0] != pr) {
                    (FALTA, Motivo::RefDiferente)
                } else {
                    (FALTA, Motivo::OutroAlelo)
                };
            };
            let g = r.gts[s];
            if g.n == 0 {
                return (FALTA, Motivo::SemGt);
            }
            let mut x = 0u8;
            for &a in g.alelos() {
                if a == 0 {
                    continue;
                }
                match r.alelo(a) {
                    Some(seq) if seq.len() == 1 && seq[0] == pa => x += 1,
                    _ => return (FALTA, Motivo::OutroAlelo),
                }
            }
            (if g.n == 2 { x } else { 3 + x }, Motivo::Genotipado)
        }
    }
}

/// The data of one person for the EM.
struct DadosEm<'a> {
    k: usize,
    p: &'a [f32],
    g: &'a [u8],
    bloco: &'a [u32],
}

/// One EM update from `q` into `novo`; returns the log-likelihood AT `q`.
/// `peso`: multiplicity of each block (bootstrap), or all ones.
fn passo_em(d: &DadosEm<'_>, q: &[f64], peso: Option<&[u32]>, novo: &mut [f64]) -> f64 {
    let k = d.k;
    novo.iter_mut().for_each(|x| *x = 0.0);
    let (mut total, mut ll) = (0.0f64, 0.0f64);
    for (j, &c) in d.g.iter().enumerate() {
        if c == FALTA {
            continue;
        }
        let w = match peso {
            Some(w) => {
                let x = w[d.bloco[j] as usize];
                if x == 0 {
                    continue;
                }
                f64::from(x)
            }
            None => 1.0,
        };
        let (x, n, lc) = DECODIFICA[c as usize];
        let p = &d.p[j * k..j * k + k];
        let mut f = 0.0;
        for i in 0..k {
            f += q[i] * f64::from(p[i]);
        }
        let a = x / f;
        let b = (n - x) / (1.0 - f);
        for i in 0..k {
            let pi = f64::from(p[i]);
            novo[i] += w * q[i] * (pi * a + (1.0 - pi) * b);
        }
        total += w * n;
        let mut l = lc;
        if x > 0.0 {
            l += x * f.ln();
        }
        if n > x {
            l += (n - x) * (1.0 - f).ln();
        }
        ll += w * l;
    }
    if total > 0.0 {
        novo.iter_mut().for_each(|x| *x /= total);
    }
    ll
}

/// Floor of an extrapolated proportion: above zero (EM cannot leave zero),
/// small enough that a group absent from the person costs nothing.
const QMIN: f64 = 1e-12;

fn projeta(q: &mut [f64]) {
    for x in q.iter_mut() {
        if x.is_nan() || *x <= QMIN {
            *x = QMIN;
        }
    }
    let s: f64 = q.iter().sum();
    q.iter_mut().for_each(|x| *x /= s);
}

struct Ajuste {
    q: Vec<f64>,
    ll: f64,
    iteracoes: usize,
    convergiu: bool,
}

/// Log-likelihood, gradient and Hessian in q (the Hessian is negative
/// semi-definite: the log-likelihood is concave in q).
fn derivadas(d: &DadosEm<'_>, q: &[f64], peso: Option<&[u32]>, g: &mut [f64], h: &mut [f64]) -> f64 {
    let k = d.k;
    g.iter_mut().for_each(|x| *x = 0.0);
    h.iter_mut().for_each(|x| *x = 0.0);
    let mut ll = 0.0;
    let mut p = vec![0.0f64; k];
    for (j, &c) in d.g.iter().enumerate() {
        if c == FALTA {
            continue;
        }
        let w = match peso {
            Some(w) => {
                let x = w[d.bloco[j] as usize];
                if x == 0 {
                    continue;
                }
                f64::from(x)
            }
            None => 1.0,
        };
        let (x, n, lc) = DECODIFICA[c as usize];
        let mut f = 0.0;
        for i in 0..k {
            p[i] = f64::from(d.p[j * k + i]);
            f += q[i] * p[i];
        }
        let a = x / f - (n - x) / (1.0 - f);
        let b = x / (f * f) + (n - x) / ((1.0 - f) * (1.0 - f));
        for i in 0..k {
            g[i] += w * a * p[i];
            let wb = w * b * p[i];
            for l in i..k {
                h[i * k + l] -= wb * p[l];
            }
        }
        let mut t = lc;
        if x > 0.0 {
            t += x * f.ln();
        }
        if n > x {
            t += (n - x) * (1.0 - f).ln();
        }
        ll += w * t;
    }
    for i in 0..k {
        for l in 0..i {
            h[i * k + l] = h[l * k + i];
        }
    }
    ll
}

/// Solves `a x = b` (n x n, row-major) by Gaussian elimination with partial
/// pivoting. `None` when singular.
fn resolve_sistema(mut a: Vec<f64>, mut b: Vec<f64>, n: usize) -> Option<Vec<f64>> {
    for c in 0..n {
        let piv = (c..n).max_by(|&x, &y| a[x * n + c].abs().total_cmp(&a[y * n + c].abs()))?;
        if a[piv * n + c].abs() < 1e-300 {
            return None;
        }
        if piv != c {
            for l in 0..n {
                a.swap(c * n + l, piv * n + l);
            }
            b.swap(c, piv);
        }
        for r in c + 1..n {
            let f = a[r * n + c] / a[c * n + c];
            if f != 0.0 {
                for l in c..n {
                    a[r * n + l] -= f * a[c * n + l];
                }
                b[r] -= f * b[c];
            }
        }
    }
    let mut x = vec![0.0; n];
    for c in (0..n).rev() {
        let mut s = b[c];
        for l in c + 1..n {
            s -= a[c * n + l] * x[l];
        }
        x[c] = s / a[c * n + c];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// Finishes the maximisation with Newton steps on the simplex, with an
/// active set for the proportions at zero (the KKT conditions: a group left
/// at zero has a gradient no larger than the common gradient of the others).
/// EM alone creeps towards a boundary by a near-constant factor per update;
/// this ends at the exact maximum. Returns (log-likelihood, steps, converged).
fn polimento_newton(d: &DadosEm<'_>, q: &mut [f64], peso: Option<&[u32]>, max_passos: usize, tol: f64) -> (f64, usize, bool) {
    let k = d.k;
    let mut livre: Vec<bool> = q.iter().map(|&x| x > 1e-7).collect();
    for i in 0..k {
        if !livre[i] {
            q[i] = 0.0;
        }
    }
    let s: f64 = q.iter().sum();
    q.iter_mut().for_each(|x| *x /= s);
    let (mut g, mut h) = (vec![0.0; k], vec![0.0; k * k]);
    let mut lixo = vec![0.0; k];
    let mut ll = derivadas(d, q, peso, &mut g, &mut h);
    for passo in 1..=max_passos {
        let f: Vec<usize> = (0..k).filter(|&i| livre[i]).collect();
        let m = f.len();
        // [H_FF  -1] [dq]   [-g_F]
        // [ 1ᵀ    0] [mu] = [  0 ]
        let mut a = vec![0.0; (m + 1) * (m + 1)];
        let mut b = vec![0.0; m + 1];
        for (r, &i) in f.iter().enumerate() {
            for (c, &l) in f.iter().enumerate() {
                a[r * (m + 1) + c] = h[i * k + l];
            }
            a[r * (m + 1) + r] -= 1e-9 * (1.0 + h[i * k + i].abs());
            a[r * (m + 1) + m] = -1.0;
            a[m * (m + 1) + r] = 1.0;
            b[r] = -g[i];
        }
        let Some(sol) = resolve_sistema(a, b, m + 1) else {
            return (ll, passo, false);
        };
        let mu = sol[m];
        let mut dq = vec![0.0; k];
        for (r, &i) in f.iter().enumerate() {
            dq[i] = sol[r];
        }
        let ganho: f64 = (0..k).map(|i| g[i] * dq[i]).sum::<f64>()
            + 0.5 * (0..k).map(|i| (0..k).map(|l| dq[i] * h[i * k + l] * dq[l]).sum::<f64>()).sum::<f64>();
        if ganho.abs() < tol || dq.iter().all(|x| x.abs() < 1e-13) {
            // Optimal on the free set: release a group at zero whose gradient
            // beats the common one, or stop.
            let solto = (0..k).filter(|&i| !livre[i] && g[i] > mu + 1e-7 * (1.0 + mu.abs())).max_by(|&x, &y| g[x].total_cmp(&g[y]));
            match solto {
                Some(i) => {
                    livre[i] = true;
                    continue;
                }
                None => return (ll, passo, true),
            }
        }
        // The largest step that keeps every proportion >= 0.
        let mut t = 1.0f64;
        let mut bate: Option<usize> = None;
        for &i in &f {
            if dq[i] < 0.0 && q[i] + t * dq[i] < 0.0 {
                t = q[i] / -dq[i];
                bate = Some(i);
            }
        }
        let mut novo = vec![0.0; k];
        let mut aceito = false;
        for _ in 0..60 {
            for i in 0..k {
                novo[i] = (q[i] + t * dq[i]).max(0.0);
            }
            if let Some(i) = bate {
                if t >= q[i] / -dq[i] * (1.0 - 1e-12) {
                    novo[i] = 0.0;
                }
            }
            let s: f64 = novo.iter().sum();
            novo.iter_mut().for_each(|x| *x /= s);
            let l = passo_em(d, &novo, peso, &mut lixo);
            if l.is_finite() && l >= ll - 1e-12 * ll.abs().max(1.0) {
                aceito = true;
                break;
            }
            t *= 0.5;
            bate = None;
        }
        if !aceito {
            return (ll, passo, true);
        }
        q.copy_from_slice(&novo);
        for i in 0..k {
            if q[i] <= 0.0 {
                livre[i] = false;
                q[i] = 0.0;
            }
        }
        ll = derivadas(d, q, peso, &mut g, &mut h);
    }
    (ll, max_passos, false)
}

/// EM updates before the Newton finish. EM is robust far from the maximum
/// but slow near it, above all along the ridge of two similar groups (EUR
/// and SAS); Newton is fast near it.
const AQUECIMENTO_EM: usize = 30;

/// The maximum likelihood proportions. EM (the ADMIXTURE / fastNGSadmix
/// update for fixed P) accelerated by SQUAREM (Varadhan & Roland 2008,
/// scheme S3) brings q near the maximum -- a step that would lower the
/// likelihood falls back to plain EM updates -- and Newton steps on the
/// simplex finish it. `aquecimento`: at most this many EM updates first (0:
/// Newton only, for a start already near the maximum, as a bootstrap
/// replicate started at the estimate).
fn ajusta(d: &DadosEm<'_>, q0: &[f64], peso: Option<&[u32]>, aquecimento: usize, max_iter: usize, tol: f64) -> Ajuste {
    let k = d.k;
    let mut q = q0.to_vec();
    let (mut q1, mut q2, mut qx, mut q3) = (vec![0.0; k], vec![0.0; k], vec![0.0; k], vec![0.0; k]);
    if k == 1 {
        let ll = passo_em(d, &q, peso, &mut q1);
        return Ajuste { q, ll, iteracoes: 1, convergiu: true };
    }
    let newton_max = max_iter.clamp(1, 200);
    if aquecimento == 0 {
        let (ll, passos, ok) = polimento_newton(d, &mut q, peso, newton_max, tol);
        return Ajuste { q, ll, iteracoes: passos, convergiu: ok };
    }
    let mut it = 0usize;
    loop {
        let ll0 = passo_em(d, &q, peso, &mut q1);
        let ll1 = passo_em(d, &q1, peso, &mut q2);
        it += 2;
        // Near the maximum (or at the end of the warm-up): the Newton finish.
        if (ll1 - ll0).abs() < 1e-3 || it >= aquecimento.min(max_iter) {
            q.copy_from_slice(&q1);
            let (ll, passos, ok) = polimento_newton(d, &mut q, peso, newton_max, tol);
            let (q, ll) = if ll >= ll1 { (q, ll) } else { (q1.clone(), ll1) };
            return Ajuste { q, ll, iteracoes: it + passos, convergiu: ok };
        }
        let (mut sr, mut sv) = (0.0, 0.0);
        for i in 0..k {
            let r = q1[i] - q[i];
            let v = q2[i] - 2.0 * q1[i] + q[i];
            sr += r * r;
            sv += v * v;
        }
        let mut alfa = if sv > 0.0 { -(sr / sv).sqrt() } else { -1.0 };
        if alfa > -1.0 {
            alfa = -1.0;
        }
        for i in 0..k {
            let r = q1[i] - q[i];
            let v = q2[i] - 2.0 * q1[i] + q[i];
            qx[i] = q[i] - 2.0 * alfa * r + alfa * alfa * v;
        }
        projeta(&mut qx);
        let llx = passo_em(d, &qx, peso, &mut q3);
        it += 1;
        if llx.is_finite() && llx >= ll1 {
            q.copy_from_slice(&q3);
        } else {
            q.copy_from_slice(&q2);
        }
        projeta(&mut q);
    }
}

fn splitmix(s: &mut u64) -> u64 {
    *s = s.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *s;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[derive(Clone, Default)]
struct ContaSitios {
    genotipados: usize,
    assumidos: usize,
    nao_chamaveis: usize,
    sob_indel: usize,
    outro_alelo: usize,
    ref_diferente: usize,
    sem_gt: usize,
    filtrados: usize,
}

impl ContaSitios {
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

/// Estimates the ancestry proportions of each person against the groups of
/// a panel written by `build_panel` / `panel_from_frequencies`.
///
/// Persons that share a VCF are read in one pass over it. Only the panel's
/// biallelic SNVs are used. Results come in the order of `persons`.
pub fn estimate_ancestry(panel: &str, persons: &[Person], opts: &AncestryOptions) -> Result<Vec<AncestryResult>, String> {
    if persons.is_empty() {
        return Err("No person given: give at least one (VCF, sample, callable regions).".to_string());
    }
    if opts.block_bp == 0 {
        return Err("block_bp is 0: give the bootstrap block length in bases (e.g. 5000000).".to_string());
    }
    if opts.tolerance.is_nan() || opts.tolerance <= 0.0 || opts.max_iterations == 0 {
        return Err("tolerance must be above 0 and max_iterations at least 1.".to_string());
    }
    if opts.bootstrap > 100_000 {
        return Err(format!("{} bootstrap replicates is too many: 100 to 1000 is plenty.", opts.bootstrap));
    }
    let painel = PainelSnv::le(panel)?;
    if painel.n() == 0 {
        return Err(format!(
            "{panel}: no biallelic SNV in the panel ({} row(s) read, {} of them not SNVs): ancestry uses the SNV rows.",
            painel.linhas, painel.nao_snv
        ));
    }
    let k = painel.grupos.len();
    let n = painel.n();
    let regioes = le_regioes_das_pessoas(persons)?;

    // Bootstrap blocks: consecutive ids along each chromosome.
    let mut bloco = vec![0u32; n];
    let mut nb = 0u32;
    for c in 0..painel.cromos.len() {
        let mut ultimo: Option<u64> = None;
        let (ini, fim) = (painel.inicio[c], painel.inicio[c + 1]);
        for (id, &pos) in bloco[ini..fim].iter_mut().zip(&painel.pos[ini..fim]) {
            let b = pos / opts.block_bp;
            if ultimo != Some(b) {
                nb += 1;
                ultimo = Some(b);
            }
            *id = nb - 1;
        }
    }

    let alvos: HashMap<String, &[u64]> = painel
        .cromos
        .iter()
        .enumerate()
        .map(|(i, c)| (c.clone(), &painel.pos[painel.inicio[i]..painel.inicio[i + 1]]))
        .collect();
    let offset: HashMap<&str, usize> = painel.cromos.iter().enumerate().map(|(i, c)| (c.as_str(), painel.inicio[i])).collect();

    let mut genos: Vec<Vec<u8>> = vec![Vec::new(); persons.len()];
    let mut contas: Vec<ContaSitios> = vec![ContaSitios::default(); persons.len()];
    let mut nomes: Vec<String> = vec![String::new(); persons.len()];
    for (vcf, idxs) in por_vcf(persons) {
        let pedidas: Vec<Option<String>> = idxs.iter().map(|&i| persons[i].sample.clone()).collect();
        let regs: Vec<Option<&Regioes>> = idxs.iter().map(|&i| persons[i].callable.as_ref().map(|b| &regioes[b])).collect();
        let mut g: Vec<Vec<u8>> = idxs.iter().map(|_| vec![FALTA; n]).collect();
        let mut c: Vec<ContaSitios> = vec![ContaSitios::default(); idxs.len()];
        let (mut crom_atual, mut off) = (String::new(), 0usize);
        let mut f = |crom: &str, t: usize, s: usize, obs: &Obs<'_>| {
            if crom != crom_atual {
                crom_atual = crom.to_string();
                off = offset[crom];
            }
            let j = off + t;
            let [pr, pa] = painel.ra[j];
            let (codigo, motivo) = classifica_snv(obs, s, pr, pa, opts.pass_only);
            g[s][j] = codigo;
            c[s].conta(motivo);
        };
        let r = junta(&vcf, &pedidas, &regs, opts.absent_is_ref, &alvos, &mut f)?;
        for (s, &i) in idxs.iter().enumerate() {
            genos[i] = std::mem::take(&mut g[s]);
            contas[i] = c[s].clone();
            nomes[i] = r.nomes[s].clone();
        }
    }

    let resultados: Vec<Result<AncestryResult, String>> = (0..persons.len())
        .into_par_iter()
        .map(|i| {
            let d = DadosEm { k, p: &painel.p, g: &genos[i], bloco: &bloco };
            ancestralidade_de(&d, &persons[i], &nomes[i], &contas[i], &painel, nb as usize, opts)
        })
        .collect();
    resultados.into_iter().collect()
}

fn ancestralidade_de(
    d: &DadosEm<'_>,
    pessoa: &Person,
    nome: &str,
    c: &ContaSitios,
    painel: &PainelSnv,
    nb: usize,
    o: &AncestryOptions,
) -> Result<AncestryResult, String> {
    let k = d.k;
    let usados = c.genotipados + c.assumidos;
    if usados == 0 {
        return Err(format!(
            "{} (sample {nome}): none of the {} panel SNVs could be used ({} not callable, {} with another allele, {} with a different REF, {} missing). Check that the VCF and the panel are on the same genome build and share chromosomes; for a VCF that lists only variants, give the callable-regions BED or keep absent_is_ref on.",
            pessoa.vcf,
            painel.n(),
            c.nao_chamaveis,
            c.outro_alelo,
            c.ref_diferente,
            c.sem_gt + c.filtrados + c.sob_indel
        ));
    }
    let uniforme = vec![1.0 / k as f64; k];
    let a = ajusta(d, &uniforme, None, AQUECIMENTO_EM, o.max_iterations, o.tolerance);
    let mut lixo = vec![0.0; k];
    let isolados: Vec<f64> = (0..k)
        .map(|g| {
            let mut e = vec![0.0; k];
            e[g] = 1.0;
            passo_em(d, &e, None, &mut lixo)
        })
        .collect();
    let melhor = (0..k).fold(0, |m, g| if isolados[g] > isolados[m] { g } else { m });

    // Blocks that hold data for this person.
    let mut presentes: Vec<u32> = Vec::new();
    for (j, &g) in d.g.iter().enumerate() {
        if g != FALTA && presentes.last() != Some(&d.bloco[j]) {
            presentes.push(d.bloco[j]);
        }
    }
    let mut erros = vec![0.0; k];
    let mut reps = 0;
    // The replicates of a person depend on the seed and on the sample name --
    // not on the place in the list, nor on where the file is: the same person
    // gets the same standard errors alone or among others, here or elsewhere.
    let mut semente = 0xcbf2_9ce4_8422_2325u64;
    for b in nome.bytes() {
        semente = (semente ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3);
    }
    if o.bootstrap > 0 && presentes.len() >= 2 {
        let qs: Vec<Vec<f64>> = (0..o.bootstrap)
            .into_par_iter()
            .map(|r| {
                let mut s = o.seed ^ semente ^ (r as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15);
                let mut w = vec![0u32; nb];
                let m = presentes.len() as u64;
                for _ in 0..presentes.len() {
                    w[presentes[(splitmix(&mut s) % m) as usize] as usize] += 1;
                }
                // A replicate is near the estimate: Newton from there.
                ajusta(d, &a.q, Some(&w), 0, o.max_iterations, o.tolerance).q
            })
            .collect();
        reps = qs.len();
        for g in 0..k {
            let media = qs.iter().map(|q| q[g]).sum::<f64>() / reps as f64;
            let var = qs.iter().map(|q| (q[g] - media).powi(2)).sum::<f64>() / (reps as f64 - 1.0).max(1.0);
            erros[g] = var.sqrt();
        }
    }

    let mut avisos = Vec::new();
    if usados < 1000 {
        avisos.push(format!("only {usados} panel sites were usable: the proportions are rough"));
    }
    if reps > 0 && presentes.len() < 20 {
        avisos.push(format!(
            "the standard errors come from only {} blocks of {} bp: they are rough (use more chromosomes or a smaller block_bp)",
            presentes.len(),
            o.block_bp
        ));
    }
    if o.bootstrap > 0 && presentes.len() < 2 {
        avisos.push("the data fall in a single bootstrap block: no standard errors".to_string());
    }
    if !a.convergiu {
        avisos.push(format!("the EM stopped at max_iterations ({}) before converging", o.max_iterations));
    }
    let vistos = c.genotipados + c.ref_diferente + c.outro_alelo;
    if vistos > 100 && c.ref_diferente * 20 > vistos {
        avisos.push(format!(
            "the REF of the VCF differs from the panel's at {} of {} positions found in both: are they on the same genome build (GRCh37 vs GRCh38)?",
            c.ref_diferente, vistos
        ));
    }
    if painel.arredondado {
        avisos.push("the panel's frequencies have at most 2 decimals: groups where an allele is rare are poorly described (a panel built from genotypes is more precise)".to_string());
    }
    if painel.duplicados > 0 || painel.sem_dados > 0 {
        avisos.push(format!(
            "panel rows left out: {} duplicate site(s), {} with a frequency missing",
            painel.duplicados, painel.sem_dados
        ));
    }
    if pessoa.callable.is_none() && o.absent_is_ref && c.assumidos > 0 {
        avisos.push(format!(
            "{} panel sites absent from the VCF were taken as homozygous reference (absent_is_ref, no callable-regions BED): right for a whole-genome VCF, wrong where the genome was not assessed",
            c.assumidos
        ));
    }
    Ok(AncestryResult {
        vcf: pessoa.vcf.clone(),
        sample: nome.to_string(),
        groups: painel.grupos.clone(),
        proportions: a.q,
        standard_errors: erros,
        panel_sites: painel.n(),
        sites_used: usados,
        sites_genotyped: c.genotipados,
        sites_assumed_ref: c.assumidos,
        sites_not_callable: c.nao_chamaveis,
        sites_in_deletion: c.sob_indel,
        sites_other_allele: c.outro_alelo,
        sites_ref_mismatch: c.ref_diferente,
        sites_missing_genotype: c.sem_gt,
        sites_filtered: c.filtrados,
        log_likelihood: a.ll,
        iterations: a.iteracoes,
        converged: a.convergiu,
        best_group: painel.grupos[melhor].clone(),
        single_group_log_likelihoods: isolados,
        blocks: presentes.len(),
        bootstrap_replicates: reps,
        warnings: avisos,
    })
}

// =====================================================================
//  3. POLYGENIC SCORE
// =====================================================================

/// Allele frequency used for score variants the person's file cannot call.
#[derive(Clone, Debug, PartialEq)]
pub enum Impute {
    /// Pooled over the panel's groups (weighted by called alleles when known).
    Pooled,
    /// The frequency of one group of the panel.
    Group(String),
    /// A mixture of the groups (e.g. the person's ancestry proportions).
    Proportions(Vec<f64>),
}

#[derive(Clone, Debug, PartialEq)]
pub struct ScoreOptions {
    /// Frequency panel: fills variants the person's file cannot call, tells
    /// which allele is REF at absent sites, and gives each group's expected
    /// mean and variance.
    pub panel: Option<String>,
    /// Reference FASTA: tells which allele is REF at absent SNV sites the
    /// panel does not have.
    pub reference: Option<String>,
    pub absent_is_ref: bool,
    pub pass_only: bool,
    /// Keep A/T and C/G SNVs (the strand cannot be checked on them).
    pub keep_ambiguous: bool,
    /// Match an SNV whose alleles are the complement of the file's.
    pub allow_flips: bool,
    pub impute: Impute,
    /// Optional per-variant table (TSV): what was done with each variant.
    pub details: Option<String>,
}

impl Default for ScoreOptions {
    fn default() -> Self {
        ScoreOptions {
            panel: None,
            reference: None,
            absent_is_ref: true,
            pass_only: true,
            keep_ambiguous: false,
            allow_flips: true,
            impute: Impute::Pooled,
            details: None,
        }
    }
}

/// The score against one group of the panel.
#[derive(Clone, Debug, PartialEq)]
pub struct GroupScore {
    pub group: String,
    /// Score variants observed in the person AND present in the panel: the
    /// mean, sd and z below are over these only.
    pub variants: usize,
    /// The person's score over those variants.
    pub observed: f64,
    /// Sum of 2 p w (expected score in the group).
    pub mean: f64,
    /// Square root of the sum of 2 p (1 - p) w^2.
    pub sd: f64,
    pub z: f64,
    /// 100 x normal CDF of z.
    pub percentile: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ScoreResult {
    pub vcf: String,
    pub sample: String,
    pub pgs_id: String,
    pub pgs_name: String,
    pub trait_reported: String,
    /// Build of the positions used (HmPOS_build for hm_pos, else genome_build).
    pub genome_build: String,
    /// "hm_pos" (harmonized) or "chr_position" (as the authors gave them).
    pub positions_from: String,
    /// Data rows of the scoring file.
    pub variants: usize,
    /// With a position and A/C/G/T alleles, not duplicated.
    pub usable: usize,
    pub no_position: usize,
    pub unsupported_allele: usize,
    pub duplicates: usize,
    /// Found in the VCF with matching alleles.
    pub genotyped: usize,
    /// Absent from the VCF but callable: homozygous reference.
    pub assumed_ref: usize,
    /// Not callable; 2 x panel frequency.
    pub imputed: usize,
    /// A/T or C/G pairs met (kept or not).
    pub ambiguous: usize,
    pub ambiguous_skipped: usize,
    /// Matched on the opposite strand.
    pub flipped: usize,
    /// The person (or the reference genome) has other alleles there.
    pub mismatched: usize,
    /// Homozygous reference, but nothing tells which allele is REF (give the
    /// panel or a reference FASTA).
    pub ref_unknown: usize,
    /// Not callable and not in the panel: left out.
    pub not_callable: usize,
    /// genotyped + assumed_ref + imputed.
    pub score: f64,
    /// genotyped + assumed_ref only.
    pub score_observed: f64,
    /// (genotyped + assumed_ref) / variants.
    pub coverage: f64,
    /// (genotyped + assumed_ref + imputed) / variants.
    pub coverage_with_imputed: f64,
    pub groups: Vec<GroupScore>,
    pub warnings: Vec<String>,
}

impl ScoreResult {
    pub fn json(&self) -> String {
        let grupos: Vec<String> = self
            .groups
            .iter()
            .map(|g| {
                format!(
                    "{{\"group\":{},\"variants\":{},\"observed\":{},\"mean\":{},\"sd\":{},\"z\":{},\"percentile\":{}}}",
                    json_texto(&g.group),
                    g.variants,
                    jn(g.observed),
                    jn(g.mean),
                    jn(g.sd),
                    jn(g.z),
                    jn(g.percentile)
                )
            })
            .collect();
        format!(
            "{{\"vcf\":{},\"sample\":{},\"pgs_id\":{},\"pgs_name\":{},\"trait_reported\":{},\"genome_build\":{},\"positions_from\":{},\"variants\":{},\"usable\":{},\"no_position\":{},\"unsupported_allele\":{},\"duplicates\":{},\"genotyped\":{},\"assumed_ref\":{},\"imputed\":{},\"ambiguous\":{},\"ambiguous_skipped\":{},\"flipped\":{},\"mismatched\":{},\"ref_unknown\":{},\"not_callable\":{},\"score\":{},\"score_observed\":{},\"coverage\":{},\"coverage_with_imputed\":{},\"groups\":[{}],\"warnings\":{}}}",
            json_texto(&self.vcf),
            json_texto(&self.sample),
            json_texto(&self.pgs_id),
            json_texto(&self.pgs_name),
            json_texto(&self.trait_reported),
            json_texto(&self.genome_build),
            json_texto(&self.positions_from),
            self.variants,
            self.usable,
            self.no_position,
            self.unsupported_allele,
            self.duplicates,
            self.genotyped,
            self.assumed_ref,
            self.imputed,
            self.ambiguous,
            self.ambiguous_skipped,
            self.flipped,
            self.mismatched,
            self.ref_unknown,
            self.not_callable,
            jn(self.score),
            jn(self.score_observed),
            jn(self.coverage),
            jn(self.coverage_with_imputed),
            grupos.join(","),
            jlista_txt(&self.warnings)
        )
    }

    pub fn resumo(&self) -> String {
        let id = if self.pgs_id.is_empty() { "score".to_string() } else { self.pgs_id.clone() };
        let mut o = format!(
            "{id}{} for {} ({}): {:.6} over {} variant(s) ({:.6} without the imputed ones)\n",
            if self.trait_reported.is_empty() { String::new() } else { format!(" ({})", self.trait_reported) },
            self.sample,
            self.vcf,
            self.score,
            self.genotyped + self.assumed_ref + self.imputed,
            self.score_observed
        );
        o.push_str(&format!(
            "  {} variant(s) in the file: {} genotyped, {} assumed homozygous reference, {} imputed from the panel; coverage {:.1}% ({:.1}% with imputed)\n",
            self.variants,
            self.genotyped,
            self.assumed_ref,
            self.imputed,
            100.0 * self.coverage,
            100.0 * self.coverage_with_imputed
        ));
        o.push_str(&format!(
            "  left out: {} strand-ambiguous, {} other alleles, {} REF unknown, {} not callable, {} without position, {} unsupported allele, {} duplicate(s); {} matched on the other strand\n",
            self.ambiguous_skipped, self.mismatched, self.ref_unknown, self.not_callable, self.no_position, self.unsupported_allele, self.duplicates, self.flipped
        ));
        for g in &self.groups {
            o.push_str(&format!(
                "  vs {:<8} {} variant(s): observed {:.4}, expected {:.4} +/- {:.4}, z {:.2}, percentile {:.1}\n",
                g.group, g.variants, g.observed, g.mean, g.sd, g.z, g.percentile
            ));
        }
        for w in &self.warnings {
            o.push_str(&format!("  note: {w}\n"));
        }
        o.push_str("  A polygenic score is a weighted sum of allele counts, for research and teaching; the z-scores compare it with the panel's groups under Hardy-Weinberg and independent variants.\n");
        o
    }
}

struct VarEscore {
    pos: u64,
    efeito: Vec<u8>,
    /// Possible other alleles; empty when the file does not say.
    outros: Vec<Vec<u8>>,
    peso: f64,
}

struct Escore {
    meta: HashMap<String, String>,
    por_crom: HashMap<String, Vec<VarEscore>>,
    linhas: usize,
    sem_posicao: usize,
    alelo_nao_suportado: usize,
    duplicados: usize,
    fonte: &'static str,
    avisos: Vec<String>,
}

fn verdadeiro(v: &[u8]) -> bool {
    matches!(v.to_ascii_lowercase().as_slice(), b"true" | b"t" | b"1" | b"yes")
}

fn le_escore(caminho: &str) -> Result<Escore, String> {
    let mut l = Linhas::abre(caminho)?;
    let mut meta: HashMap<String, String> = HashMap::new();
    let cab: Vec<String> = loop {
        if !l.proxima()? {
            return Err(format!("{caminho}: no header line with the column names: is this a PGS Catalog scoring file?"));
        }
        let t = txt(&l.buf);
        if t.trim().is_empty() {
            continue;
        }
        if let Some(m) = t.strip_prefix('#') {
            let m = m.trim_start_matches('#');
            if let Some((a, b)) = m.split_once('=') {
                meta.insert(a.trim().to_string(), b.trim().to_string());
            }
            continue;
        }
        break t.split('\t').map(|s| s.trim().to_string()).collect();
    };
    let acha = |n: &str| cab.iter().position(|x| x.eq_ignore_ascii_case(n));
    for d in ["dosage_0_weight", "dosage_1_weight", "dosage_2_weight"] {
        if acha(d).is_some() {
            return Err(format!(
                "{caminho} has dosage-specific weights ({d}): the score is not additive in the allele dosage and is not supported. Choose a score with one effect_weight per variant."
            ));
        }
    }
    let (Some(ie), Some(iw)) = (acha("effect_allele"), acha("effect_weight")) else {
        return Err(format!(
            "{caminho}: the scoring file needs the columns effect_allele and effect_weight (it has: {}). Is this a PGS Catalog scoring file?",
            cab.join(", ")
        ));
    };
    let io = acha("other_allele").or_else(|| acha("reference_allele"));
    let iinf = acha("hm_inferOtherAllele");
    let (ic, ip, fonte) = match (acha("hm_chr"), acha("hm_pos"), acha("chr_name"), acha("chr_position")) {
        (Some(c), Some(p), _, _) => (c, p, "hm_pos"),
        (_, _, Some(c), Some(p)) => (c, p, "chr_position"),
        _ => {
            return Err(format!(
                "{caminho} has no variant positions (only rsIDs?). Download the harmonized file from the PGS Catalog (the one ending in _hmPOS_GRCh38.txt.gz, or _hmPOS_GRCh37 for a GRCh37 VCF): it has hm_chr and hm_pos."
            ))
        }
    };
    let recusas: Vec<(usize, &str)> = ["is_haplotype", "is_diplotype", "is_interaction", "is_dominant", "is_recessive"]
        .iter()
        .filter_map(|n| acha(n).map(|i| (i, *n)))
        .collect();
    let mut avisos = Vec::new();
    if let Some(wt) = meta.get("weight_type") {
        if wt.eq_ignore_ascii_case("OR") || wt.eq_ignore_ascii_case("HR") {
            avisos.push(format!(
                "weight_type is {wt}: if effect_weight holds odds or hazard ratios rather than their logarithms, the sum is not meaningful (the PGS Catalog asks for log(OR)/log(HR) in effect_weight)"
            ));
        }
    }
    if fonte == "chr_position" {
        avisos.push(format!(
            "positions are the ones the authors reported (genome_build={}): check that the VCF uses the same build, or use the harmonized _hmPOS_ file",
            meta.get("genome_build").map(String::as_str).unwrap_or("not stated")
        ));
    }

    let mut e = Escore { meta, por_crom: HashMap::new(), linhas: 0, sem_posicao: 0, alelo_nao_suportado: 0, duplicados: 0, fonte, avisos };
    let mut canon = Canon::default();
    let mut sp = Vec::new();
    while l.proxima()? {
        if l.buf.is_empty() || l.buf[0] == b'#' {
            continue;
        }
        e.linhas += 1;
        divide(&l.buf, &mut sp);
        let (buf, sp) = (&l.buf[..], &sp[..]);
        let campo = move |i: usize| -> &[u8] {
            match sp.get(i) {
                Some(&(a, b)) => buf[a..b].trim_ascii(),
                None => &[],
            }
        };
        for &(i, nome) in &recusas {
            if verdadeiro(campo(i)) {
                return Err(format!(
                    "{caminho}, line {}: {nome} is True. Haplotype, diplotype, interaction and dominant/recessive weights are not supported: only additive per-allele weights are."
                , l.n));
            }
        }
        let peso = real(campo(iw)).ok_or_else(|| format!("{caminho}, line {}: effect_weight `{}` is not a number.", l.n, trecho(campo(iw))))?;
        let (c, p) = (campo(ic), campo(ip));
        let Some(pos) = inteiro(p).filter(|_| !c.is_empty()) else {
            e.sem_posicao += 1;
            continue;
        };
        let efeito = campo(ie).to_ascii_uppercase();
        if !so_acgt(&efeito) {
            e.alelo_nao_suportado += 1;
            continue;
        }
        let mut outros: Vec<Vec<u8>> = Vec::new();
        let o = io.map(campo).unwrap_or(&[]);
        let fonte_outro = if !o.is_empty() { o } else { iinf.map(campo).unwrap_or(&[]) };
        for a in fonte_outro.split(|&b| b == b'/' || b == b',') {
            let a = a.trim_ascii().to_ascii_uppercase();
            if so_acgt(&a) && a != efeito && !outros.contains(&a) {
                outros.push(a);
            }
        }
        let chave = canon.de(c).to_string();
        e.por_crom.entry(chave).or_default().push(VarEscore { pos, efeito, outros, peso });
    }
    for v in e.por_crom.values_mut() {
        let mut ix: Vec<(u64, u32)> = v.iter().enumerate().map(|(i, x)| (x.pos, i as u32)).collect();
        multi_merge_sort_stable_by(&mut ix, |a, b| a.0.cmp(&b.0));
        let mut velho: Vec<Option<VarEscore>> = std::mem::take(v).into_iter().map(Some).collect();
        for (_, i) in ix {
            let x = velho[i as usize].take().expect("each index once");
            let dup = v.iter().rev().take_while(|y| y.pos == x.pos).any(|y| y.efeito == x.efeito && y.outros == x.outros);
            if dup {
                e.duplicados += 1;
            } else {
                v.push(x);
            }
        }
    }
    Ok(e)
}

/// A frequency written with at most 2 decimals (0.04, 0.5, 1).
fn duas_casas(x: f64) -> bool {
    ((x * 100.0) - (x * 100.0).round()).abs() < 1e-7
}

type PainelPgs = HashMap<(String, u64), Vec<LinhaPainelDona>>;

/// The panel rows at the score's positions, and whether all their
/// frequencies have at most 2 decimals.
fn painel_nas_posicoes(caminho: &str, quer: &HashMap<String, Vec<u64>>) -> Result<(Vec<String>, PainelPgs, bool), String> {
    let (mut l, c) = abre_painel(caminho)?;
    let mut out: PainelPgs = HashMap::new();
    let mut arredondado = true;
    let mut canon = Canon::default();
    let mut sp = Vec::new();
    while l.proxima()? {
        if l.buf.is_empty() || l.buf[0] == b'#' {
            continue;
        }
        divide(&l.buf, &mut sp);
        // The position first: most rows are not in the score.
        let (a, b) = sp.get(c.crom).copied().unwrap_or((0, 0));
        let chave = canon.de(&l.buf[a..b]);
        let Some(ps) = quer.get(chave) else { continue };
        let (pa, pb) = sp.get(c.pos).copied().unwrap_or((0, 0));
        let Some(pos) = inteiro(&l.buf[pa..pb]) else {
            return Err(format!("{caminho}, line {}: pos `{}` is not a whole number.", l.n, trecho(&l.buf[pa..pb])));
        };
        if ps.binary_search(&pos).is_err() {
            continue;
        }
        let Some(r) = linha_painel(caminho, l.n, &l.buf, &sp, &c)? else { continue };
        arredondado &= r.af.iter().all(|&x| duas_casas(x));
        out.entry((chave.to_string(), pos)).or_default().push(LinhaPainelDona { ref_: r.ref_, alt: r.alt, af: r.af, an: r.an });
    }
    let arredondado = arredondado && !out.is_empty();
    Ok((c.grupos, out, arredondado))
}

struct LinhaPainelDona {
    ref_: Vec<u8>,
    alt: Vec<u8>,
    af: Vec<f64>,
    an: Vec<Option<u64>>,
}

/// A panel row that has the variant's alleles: (row, the effect allele is
/// the ALT, matched on the other strand).
fn casa_painel<'a>(v: &VarEscore, linhas: &'a [LinhaPainelDona], flips: bool) -> Option<(&'a LinhaPainelDona, bool, bool)> {
    let casa = |e: &[u8], outros: &[Vec<u8>], r: &LinhaPainelDona| -> Option<bool> {
        if e == r.alt.as_slice() && (outros.is_empty() || outros.contains(&r.ref_)) {
            Some(true)
        } else if e == r.ref_.as_slice() && (outros.is_empty() || outros.contains(&r.alt)) {
            Some(false)
        } else {
            None
        }
    };
    for r in linhas {
        if let Some(alt) = casa(&v.efeito, &v.outros, r) {
            return Some((r, alt, false));
        }
    }
    if flips && v.efeito.len() == 1 && v.outros.iter().all(|o| o.len() == 1) {
        let e = vec![comp_base(v.efeito[0])];
        let outros: Vec<Vec<u8>> = v.outros.iter().map(|o| vec![comp_base(o[0])]).collect();
        for r in linhas {
            if r.ref_.len() == 1 && r.alt.len() == 1 && !complementares(&r.ref_, &r.alt) {
                if let Some(alt) = casa(&e, &outros, r) {
                    return Some((r, alt, true));
                }
            }
        }
    }
    None
}

/// Frequency of the effect allele in each group, from a matched row.
fn freqs_efeito(r: &LinhaPainelDona, efeito_alt: bool) -> Vec<f64> {
    r.af.iter().map(|&p| if efeito_alt { p } else { 1.0 - p }).collect()
}

/// Normal CDF, from erfc (Numerical Recipes' Chebyshev fit, relative error
/// below 1.2e-7).
fn normal_cdf(z: f64) -> f64 {
    let x = -z / std::f64::consts::SQRT_2;
    let a = x.abs();
    let t = 1.0 / (1.0 + 0.5 * a);
    let y = t
        * (-a * a - 1.265_512_23
            + t * (1.000_023_68
                + t * (0.374_091_96
                    + t * (0.096_784_18
                        + t * (-0.186_288_06
                            + t * (0.278_868_07 + t * (-1.135_203_98 + t * (1.488_515_87 + t * (-0.822_152_23 + t * 0.170_872_77)))))))))
            .exp();
    let erfc = if x >= 0.0 { y } else { 2.0 - y };
    0.5 * erfc
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Estado {
    Genotipado,
    Assumido,
    Imputado,
    Ambiguo,
    Diferente,
    RefDesconhecida,
    NaoChamavel,
}

impl Estado {
    fn nome(self) -> &'static str {
        match self {
            Estado::Genotipado => "genotyped",
            Estado::Assumido => "assumed_ref",
            Estado::Imputado => "imputed",
            Estado::Ambiguo => "ambiguous_skipped",
            Estado::Diferente => "mismatched",
            Estado::RefDesconhecida => "ref_unknown",
            Estado::NaoChamavel => "not_callable",
        }
    }
}

/// What happened to one score variant for one person.
struct Desfecho {
    estado: Estado,
    dosagem: f64,
    flip: bool,
    ambiguo: bool,
    /// Frequencies of the effect allele per group, when in the panel.
    p: Option<Vec<f64>>,
}

#[derive(Default)]
struct Acumula {
    n: [usize; 7],
    ambiguos: usize,
    flips: usize,
    /// Variants used on X or Y (scored as diploid).
    xy: usize,
    escore: f64,
    observado: f64,
    /// Per group: variants, observed, mean, variance.
    grupos: Vec<(usize, f64, f64, f64)>,
}

struct ContextoEscore<'a> {
    painel: Option<&'a HashMap<(String, u64), Vec<LinhaPainelDona>>>,
    genoma: Option<Genoma>,
    pesos_imputacao: Option<Vec<f64>>,
    o: &'a ScoreOptions,
}

impl ContextoEscore<'_> {
    fn linhas(&self, crom: &str, pos: u64) -> &[LinhaPainelDona] {
        self.painel.and_then(|p| p.get(&(crom.to_string(), pos))).map(Vec::as_slice).unwrap_or(&[])
    }

    fn p_imputacao(&self, r: &LinhaPainelDona, efeito_alt: bool) -> f64 {
        let p = freqs_efeito(r, efeito_alt);
        match &self.pesos_imputacao {
            Some(w) => p.iter().zip(w).map(|(a, b)| a * b).sum(),
            None => {
                let soma: u64 = r.an.iter().map(|n| n.unwrap_or(0)).sum();
                if r.an.iter().all(Option::is_some) && soma > 0 {
                    p.iter().zip(&r.an).map(|(a, n)| a * n.unwrap_or(0) as f64).sum::<f64>() / soma as f64
                } else {
                    p.iter().sum::<f64>() / p.len() as f64
                }
            }
        }
    }

    /// Not callable: 2 x frequency from the panel, or left out.
    fn imputa(&self, v: &VarEscore, crom: &str) -> Desfecho {
        match casa_painel(v, self.linhas(crom, v.pos), self.o.allow_flips) {
            Some((r, alt, flip)) => {
                let outro = if alt { &r.ref_ } else { &r.alt };
                let efeito = if alt { &r.alt } else { &r.ref_ };
                let ambiguo = complementares(efeito, outro);
                if ambiguo && !self.o.keep_ambiguous {
                    return Desfecho { estado: Estado::Ambiguo, dosagem: 0.0, flip, ambiguo, p: None };
                }
                Desfecho { estado: Estado::Imputado, dosagem: 2.0 * self.p_imputacao(r, alt), flip, ambiguo, p: None }
            }
            None => Desfecho { estado: Estado::NaoChamavel, dosagem: 0.0, flip: false, ambiguo: false, p: None },
        }
    }

    fn resolve(&mut self, v: &VarEscore, crom: &str, s: usize, obs: &Obs<'_>) -> Result<Desfecho, String> {
        let o = self.o;
        // Ambiguous from the scoring file alone (both alleles known).
        if !o.keep_ambiguous && !v.outros.is_empty() && v.outros.iter().all(|x| complementares(&v.efeito, x)) {
            return Ok(Desfecho { estado: Estado::Ambiguo, dosagem: 0.0, flip: false, ambiguo: true, p: None });
        }
        let painel = casa_painel(v, self.linhas(crom, v.pos), o.allow_flips);
        let p_grupos = painel.map(|(r, alt, _)| freqs_efeito(r, alt));
        let mut d = match obs {
            Obs::NaoChamavel | Obs::SobIndel => return Ok(self.imputa(v, crom)),
            Obs::Registros(regs) => {
                if o.pass_only && regs.iter().any(|r| !r.filtro_ok) {
                    return Ok(self.imputa(v, crom));
                }
                // Direct, then on the other strand.
                let mut achado: Option<(&RegPessoa, u8, Vec<u8>, bool)> = None;
                let tentativas: Vec<(Vec<u8>, Vec<Vec<u8>>, bool)> = {
                    let mut t = vec![(v.efeito.clone(), v.outros.clone(), false)];
                    if o.allow_flips && v.efeito.len() == 1 && v.outros.iter().all(|x| x.len() == 1) {
                        t.push((vec![comp_base(v.efeito[0])], v.outros.iter().map(|x| vec![comp_base(x[0])]).collect(), true));
                    }
                    t
                };
                'fora: for (e, outros, flip) in &tentativas {
                    for r in regs.iter() {
                        let Some(ie) = r.indice(e) else { continue };
                        let outro = if outros.is_empty() {
                            if ie == 0 {
                                r.alts.first().cloned()
                            } else {
                                Some(r.ref_.clone())
                            }
                        } else {
                            outros.iter().find(|x| r.indice(x).is_some()).cloned()
                        };
                        let Some(outro) = outro else { continue };
                        if *flip && complementares(e, &outro) {
                            continue;
                        }
                        achado = Some((r, ie, outro, *flip));
                        break 'fora;
                    }
                }
                let Some((r, ie, outro, flip)) = achado else {
                    return Ok(Desfecho { estado: Estado::Diferente, dosagem: 0.0, flip: false, ambiguo: false, p: None });
                };
                let efeito = r.alelo(ie).unwrap_or(&[]).to_vec();
                let ambiguo = complementares(&efeito, &outro);
                if ambiguo && !o.keep_ambiguous {
                    return Ok(Desfecho { estado: Estado::Ambiguo, dosagem: 0.0, flip, ambiguo, p: None });
                }
                let g = r.gts[s];
                if g.n == 0 {
                    return Ok(self.imputa(v, crom));
                }
                let x = g.alelos().iter().filter(|&&a| a == ie).count() as f64;
                Desfecho { estado: Estado::Genotipado, dosagem: x, flip, ambiguo, p: None }
            }
            Obs::Ausente => {
                if let Some((r, alt, flip)) = painel {
                    let outro = if alt { &r.ref_ } else { &r.alt };
                    let efeito = if alt { &r.alt } else { &r.ref_ };
                    let ambiguo = complementares(efeito, outro);
                    if ambiguo && !o.keep_ambiguous {
                        return Ok(Desfecho { estado: Estado::Ambiguo, dosagem: 0.0, flip, ambiguo, p: None });
                    }
                    Desfecho { estado: Estado::Assumido, dosagem: if alt { 0.0 } else { 2.0 }, flip, ambiguo, p: None }
                } else if let (Some(g), true) = (self.genoma.as_mut(), v.efeito.len() == 1 && v.outros.iter().all(|x| x.len() == 1)) {
                    let Some(nome) = g.nome(crom) else {
                        return Ok(Desfecho { estado: Estado::RefDesconhecida, dosagem: 0.0, flip: false, ambiguo: false, p: None });
                    };
                    let b = match g.trecho(&nome, v.pos - 1, v.pos)? {
                        [b] => *b,
                        _ => return Ok(Desfecho { estado: Estado::RefDesconhecida, dosagem: 0.0, flip: false, ambiguo: false, p: None }),
                    };
                    let e = v.efeito[0];
                    let tem_outro = |x: u8| v.outros.iter().any(|o| o[0] == x);
                    let (dos, flip, outro) = if e == b {
                        (2.0, false, v.outros.first().map(|o| o[0]))
                    } else if v.outros.is_empty() || tem_outro(b) {
                        (0.0, false, Some(b))
                    } else if o.allow_flips && comp_base(e) == b {
                        (2.0, true, v.outros.first().map(|o| comp_base(o[0])))
                    } else if o.allow_flips && v.outros.iter().any(|x| comp_base(x[0]) == b) {
                        (0.0, true, Some(b))
                    } else {
                        return Ok(Desfecho { estado: Estado::Diferente, dosagem: 0.0, flip: false, ambiguo: false, p: None });
                    };
                    let ef = if flip { comp_base(e) } else { e };
                    let ambiguo = outro.is_some_and(|x| comp_base(ef) == x);
                    if ambiguo && !o.keep_ambiguous {
                        return Ok(Desfecho { estado: Estado::Ambiguo, dosagem: 0.0, flip, ambiguo, p: None });
                    }
                    Desfecho { estado: Estado::Assumido, dosagem: dos, flip, ambiguo, p: None }
                } else {
                    return Ok(Desfecho { estado: Estado::RefDesconhecida, dosagem: 0.0, flip: false, ambiguo: false, p: None });
                }
            }
        };
        d.p = p_grupos;
        Ok(d)
    }
}

/// Applies a PGS Catalog scoring file to each person.
///
/// The positions are `hm_chr`/`hm_pos` when the file has them (harmonized
/// files), else `chr_name`/`chr_position`. The effect allele can be the REF
/// or the ALT of the VCF. Persons that share a VCF are read in one pass.
pub fn polygenic_score(scoring_file: &str, persons: &[Person], opts: &ScoreOptions) -> Result<Vec<ScoreResult>, String> {
    if persons.is_empty() {
        return Err("No person given: give at least one (VCF, sample, callable regions).".to_string());
    }
    let escore = le_escore(scoring_file)?;
    let usaveis: usize = escore.por_crom.values().map(Vec::len).sum();
    let posicoes: HashMap<String, Vec<u64>> = escore.por_crom.iter().map(|(c, v)| (c.clone(), v.iter().map(|x| x.pos).collect())).collect();
    let (grupos, linhas_painel, arredondado) = match &opts.panel {
        Some(p) => {
            let (g, l, a) = painel_nas_posicoes(p, &posicoes)?;
            (g, Some(l), a)
        }
        None => (Vec::new(), None, false),
    };
    let pesos_imputacao = match &opts.impute {
        Impute::Pooled => None,
        Impute::Group(nome) => {
            if opts.panel.is_none() {
                None
            } else {
                let i = grupos.iter().position(|g| g == nome).ok_or_else(|| {
                    format!("The imputation group `{nome}` is not in the panel (it has: {}).", grupos.join(", "))
                })?;
                let mut w = vec![0.0; grupos.len()];
                w[i] = 1.0;
                Some(w)
            }
        }
        Impute::Proportions(q) => {
            if opts.panel.is_none() {
                None
            } else {
                if q.len() != grupos.len() || q.iter().any(|x| x.is_nan() || *x < 0.0) || q.iter().sum::<f64>() <= 0.0 {
                    return Err(format!(
                        "The imputation proportions must be {} non-negative numbers, one per panel group ({}), and not all zero.",
                        grupos.len(),
                        grupos.join(", ")
                    ));
                }
                let s: f64 = q.iter().sum();
                Some(q.iter().map(|x| x / s).collect())
            }
        }
    };
    let genoma = match &opts.reference {
        Some(r) => Some(Genoma::abre(r, None)?),
        None => None,
    };
    let regioes = le_regioes_das_pessoas(persons)?;
    let mut entradas: Vec<&str> = vec![scoring_file];
    entradas.extend(persons.iter().map(|p| p.vcf.as_str()));
    let mut detalhes = match &opts.details {
        Some(c) => {
            confere_saida(c, &entradas)?;
            let mut w = Escritor::cria(c)?;
            w.linha("sample\tchrom\tpos\teffect_allele\tother_allele\teffect_weight\tstatus\tflipped\tdosage\tcontribution")
                .map_err(|e| format!("Error writing {c}: {e}"))?;
            Some((w, Guarda { caminho: c.clone(), armada: true }))
        }
        None => None,
    };

    let mut ctx = ContextoEscore { painel: linhas_painel.as_ref(), genoma, pesos_imputacao, o: opts };
    let alvos: HashMap<String, &[u64]> = posicoes.iter().map(|(c, v)| (c.clone(), v.as_slice())).collect();
    let k = grupos.len();
    let mut acums: Vec<Acumula> = (0..persons.len()).map(|_| Acumula { grupos: vec![(0, 0.0, 0.0, 0.0); k], ..Default::default() }).collect();
    let mut nomes: Vec<String> = vec![String::new(); persons.len()];

    for (vcf, idxs) in por_vcf(persons) {
        let pedidas: Vec<Option<String>> = idxs.iter().map(|&i| persons[i].sample.clone()).collect();
        let regs: Vec<Option<&Regioes>> = idxs.iter().map(|&i| persons[i].callable.as_ref().map(|b| &regioes[b])).collect();
        let mut ac: Vec<Acumula> = idxs.iter().map(|_| Acumula { grupos: vec![(0, 0.0, 0.0, 0.0); k], ..Default::default() }).collect();
        let mut erro: Option<String> = None;
        let mut linhas_det: Vec<(usize, String)> = Vec::new();
        {
            let mut f = |crom: &str, t: usize, s: usize, obs: &Obs<'_>| {
                if erro.is_some() {
                    return;
                }
                let v = &escore.por_crom[crom][t];
                let d = match ctx.resolve(v, crom, s, obs) {
                    Ok(d) => d,
                    Err(e) => {
                        erro = Some(e);
                        return;
                    }
                };
                let a = &mut ac[s];
                let i = d.estado as usize;
                a.n[i] += 1;
                if d.ambiguo {
                    a.ambiguos += 1;
                }
                let contrib = d.dosagem * v.peso;
                if matches!(d.estado, Estado::Genotipado | Estado::Assumido | Estado::Imputado) && (crom == "X" || crom == "Y") {
                    a.xy += 1;
                }
                match d.estado {
                    Estado::Genotipado | Estado::Assumido => {
                        if d.flip {
                            a.flips += 1;
                        }
                        a.escore += contrib;
                        a.observado += contrib;
                        if let Some(p) = &d.p {
                            for (g, &pg) in p.iter().enumerate() {
                                let x = &mut a.grupos[g];
                                x.0 += 1;
                                x.1 += contrib;
                                x.2 += 2.0 * pg * v.peso;
                                x.3 += 2.0 * pg * (1.0 - pg) * v.peso * v.peso;
                            }
                        }
                    }
                    Estado::Imputado => {
                        if d.flip {
                            a.flips += 1;
                        }
                        a.escore += contrib;
                    }
                    _ => {}
                }
                if detalhes.is_some() {
                    let usado = matches!(d.estado, Estado::Genotipado | Estado::Assumido | Estado::Imputado);
                    let outros: Vec<String> = v.outros.iter().map(|x| txt(x)).collect();
                    linhas_det.push((
                        s,
                        format!(
                            "{crom}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            v.pos,
                            txt(&v.efeito),
                            if outros.is_empty() { ".".to_string() } else { outros.join("/") },
                            v.peso,
                            d.estado.nome(),
                            d.flip,
                            if usado { format!("{}", d.dosagem) } else { ".".to_string() },
                            if usado { format!("{contrib}") } else { ".".to_string() }
                        ),
                    ));
                }
            };
            let r = junta(&vcf, &pedidas, &regs, opts.absent_is_ref, &alvos, &mut f)?;
            for (s, &i) in idxs.iter().enumerate() {
                nomes[i] = r.nomes[s].clone();
            }
        }
        if let Some(e) = erro {
            return Err(e);
        }
        if let Some((w, _)) = detalhes.as_mut() {
            for (s, l) in &linhas_det {
                w.linha(&format!("{}\t{l}", nomes[idxs[*s]])).map_err(|e| format!("Error writing the details table: {e}"))?;
            }
        }
        for (s, &i) in idxs.iter().enumerate() {
            acums[i] = std::mem::take(&mut ac[s]);
        }
    }
    if let Some((w, mut g)) = detalhes.take() {
        w.termina().map_err(|e| format!("Error writing the details table: {e}"))?;
        g.armada = false;
    }

    let meta = |k: &str| escore.meta.get(k).cloned().unwrap_or_default();
    let build = if escore.fonte == "hm_pos" { meta("HmPOS_build") } else { meta("genome_build") };
    Ok(persons
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let a = &acums[i];
            let n = |e: Estado| a.n[e as usize];
            let observados = n(Estado::Genotipado) + n(Estado::Assumido);
            let mut avisos = escore.avisos.clone();
            if opts.panel.is_none() && n(Estado::NaoChamavel) > 0 {
                avisos.push(format!(
                    "{} variant(s) could not be called in this file and were left out (give a panel to fill them with 2 x frequency)",
                    n(Estado::NaoChamavel)
                ));
            }
            if n(Estado::RefDesconhecida) > 0 {
                avisos.push(format!(
                    "{} variant(s) are absent from the VCF (homozygous reference) but neither the panel nor a reference FASTA tells which allele is REF: left out",
                    n(Estado::RefDesconhecida)
                ));
            }
            if arredondado && !grupos.is_empty() {
                avisos.push(
                    "the panel's frequencies at the score's variants have at most 2 decimals (as in the 1000 Genomes phase 3 sites file): imputed dosages and the groups' expected means and variances are approximate, and poor for rare variants".to_string(),
                );
            }
            if a.xy > 0 {
                avisos.push(format!(
                    "{} variant(s) on chromosome X or Y were scored like autosomal ones (dosage 0 to 2, or 2 x frequency when imputed): for a male genome their contribution is not comparable",
                    a.xy
                ));
            }
            if escore.linhas > 0 && observados * 2 < escore.linhas {
                avisos.push(format!(
                    "fewer than half of the score's variants were observed ({observados} of {}): check the genome build and the chromosome names",
                    escore.linhas
                ));
            }
            let grupos_r: Vec<GroupScore> = grupos
                .iter()
                .enumerate()
                .map(|(g, nome)| {
                    let (nv, obs, media, var) = a.grupos[g];
                    let sd = var.sqrt();
                    let z = if sd > 0.0 { (obs - media) / sd } else { f64::NAN };
                    GroupScore {
                        group: nome.clone(),
                        variants: nv,
                        observed: obs,
                        mean: media,
                        sd,
                        z,
                        percentile: if z.is_finite() { 100.0 * normal_cdf(z) } else { f64::NAN },
                    }
                })
                .collect();
            let total = escore.linhas.max(1) as f64;
            ScoreResult {
                vcf: p.vcf.clone(),
                sample: nomes[i].clone(),
                pgs_id: meta("pgs_id"),
                pgs_name: meta("pgs_name"),
                trait_reported: meta("trait_reported"),
                genome_build: build.clone(),
                positions_from: escore.fonte.to_string(),
                variants: escore.linhas,
                usable: usaveis,
                no_position: escore.sem_posicao,
                unsupported_allele: escore.alelo_nao_suportado,
                duplicates: escore.duplicados,
                genotyped: n(Estado::Genotipado),
                assumed_ref: n(Estado::Assumido),
                imputed: n(Estado::Imputado),
                ambiguous: a.ambiguos,
                ambiguous_skipped: n(Estado::Ambiguo),
                flipped: a.flips,
                mismatched: n(Estado::Diferente),
                ref_unknown: n(Estado::RefDesconhecida),
                not_callable: n(Estado::NaoChamavel),
                score: a.escore,
                score_observed: a.observado,
                coverage: observados as f64 / total,
                coverage_with_imputed: (observados + n(Estado::Imputado)) as f64 / total,
                groups: grupos_r,
                warnings: avisos,
            }
        })
        .collect())
}
