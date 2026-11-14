// =====================================================================
//  RNA-SEQ: READS PER GENE AND DIFFERENTIAL EXPRESSION
//
//  Two tools:
//
//   - `count_reads`             reads (or read pairs) per gene, from BAM/SAM
//                               alignments and a GFF3/GTF annotation, in the
//                               "union" mode of htseq-count and featureCounts;
//   - `differential_expression` the DESeq2 method on a count matrix: size
//                               factors by the median of ratios, gene-wise
//                               dispersions (Cox-Reid adjusted), a parametric
//                               trend, empirical-Bayes shrinkage toward it, a
//                               negative-binomial GLM, the Wald test, Cook's
//                               distances, independent filtering and the
//                               Benjamini-Hochberg adjustment.
//
//  COUNTING. A read is assigned to a gene when the reference bases its CIGAR
//  aligns (M, = and X; an N -- an intron -- and a D are not aligned bases)
//  overlap the features (exons) of exactly one gene; two or more genes make it
//  ambiguous, none makes it "no feature". The two mates of a pair are one
//  fragment: the union of the genes of both. Multimappers (NH > 1) and reads
//  below the mapping-quality threshold are counted apart and never assigned;
//  secondary, supplementary, duplicate and QC-fail records are skipped.
//
//  The alignments are streamed: what stays in memory is the annotation index
//  (a few tens of MB for a human annotation) and the mates waiting for their
//  pair (a few per read length in a coordinate-sorted file; any order works,
//  only the memory changes).
//
//  DIFFERENTIAL EXPRESSION follows the R code of DESeq2 1.4x, step by step
//  (the comments name the DESeq2 function each step reproduces). The count
//  matrix is held in memory (genes x samples x 8 bytes: 60,000 genes and 100
//  samples are 48 MB). Differences from DESeq2 are listed in `DeResult`'s
//  documentation and in the report text.
//
//  The special functions (log-gamma, digamma, trigamma, the normal and F
//  distributions, the incomplete gamma and beta functions) and lowess are
//  implemented here: no dependency beyond std and rayon.
// =====================================================================

use crate::externo::{abre_leitura, Escritor};
use crate::multimerge_core::multi_merge_sort_stable_by;
use crate::nomes::{canonico, json_texto};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::BufRead;

// ---------------------------------------------------------------------
//  READING, WRITING, JSON
// ---------------------------------------------------------------------

/// Line reader that remembers the path and the line number for messages.
struct Leitor {
    caminho: String,
    r: Box<dyn BufRead>,
    buf: String,
    linha: usize,
}

impl Leitor {
    fn abre(caminho: &str) -> Result<Leitor, String> {
        Ok(Leitor { caminho: caminho.to_string(), r: abre_leitura(caminho)?, buf: String::new(), linha: 0 })
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
                "Error reading {} after line {}: {e}. The file may be truncated, or binary in a format that is not read here (CRAM is not; convert it to BAM).",
                self.caminho, self.linha
            )),
        }
    }

    fn texto(&self) -> &str {
        self.buf.trim_end_matches(['\n', '\r'])
    }
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

fn recusa_mesmo_arquivo(entrada: &str, saida: &str) -> Result<(), String> {
    if mesmo_arquivo(entrada, saida) {
        return Err(format!(
            "The output `{saida}` is the input file `{entrada}` itself; writing it would destroy the input. Choose a different output path."
        ));
    }
    Ok(())
}

/// A number for a table: shortest exact form, scientific notation for very
/// small or very large values, `NA` for a missing one (as R writes it).
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

fn jnums(v: &[f64]) -> String {
    format!("[{}]", v.iter().map(|x| jnum(*x)).collect::<Vec<_>>().join(","))
}

fn jinteiros<T: std::fmt::Display>(v: &[T]) -> String {
    format!("[{}]", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","))
}

fn jtextos(v: &[String]) -> String {
    format!("[{}]", v.iter().map(|s| json_texto(s)).collect::<Vec<_>>().join(","))
}

fn jopt(s: Option<&str>) -> String {
    s.map(json_texto).unwrap_or_else(|| "null".into())
}

fn percent(x: usize, de: usize) -> String {
    if de == 0 {
        "-".into()
    } else {
        format!("{:.1}%", 100.0 * x as f64 / de as f64)
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

// =====================================================================
//  1) READS PER GENE
// =====================================================================

/// Which strand of the gene the reads come from (the library protocol).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strandedness {
    /// Either strand counts (featureCounts -s 0, htseq-count -s no).
    Unstranded,
    /// The read (read 1 of a pair) is on the gene's strand (-s 1, -s yes).
    Forward,
    /// The read (read 1 of a pair) is on the opposite strand: dUTP / TruSeq
    /// stranded libraries (-s 2, -s reverse).
    Reverse,
}

impl Strandedness {
    pub fn from_text(s: &str) -> Result<Strandedness, String> {
        match s.trim().to_ascii_lowercase().as_str() {
            "" | "unstranded" | "no" | "none" | "0" | "false" => Ok(Strandedness::Unstranded),
            "forward" | "yes" | "stranded" | "sense" | "1" | "fr-secondstrand" => Ok(Strandedness::Forward),
            "reverse" | "antisense" | "2" | "fr-firststrand" => Ok(Strandedness::Reverse),
            _ => Err(format!(
                "Unknown strandedness \"{s}\": use \"unstranded\" (0), \"forward\" (1: read 1 on the gene's strand) or \"reverse\" (2: read 1 on the opposite strand, as in dUTP/TruSeq stranded kits)"
            )),
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Strandedness::Unstranded => "unstranded",
            Strandedness::Forward => "forward",
            Strandedness::Reverse => "reverse",
        }
    }
}

#[derive(Clone, Debug)]
pub struct CountOptions {
    /// Column 3 of the annotation lines that are counted. Default "exon".
    pub feature_type: String,
    /// The attribute that names the gene of each feature. It is read on the
    /// feature itself or, when absent there, on its nearest ancestor through
    /// `Parent=` (GFF3). "GeneID" also reads `Dbxref=GeneID:672` (NCBI).
    /// `None`: `gene_id` (feature or ancestor), else the ID of the top-level
    /// ancestor without a "gene-"/"gene:" prefix (RefSeq: the symbol;
    /// Ensembl GFF3: the gene's gene_id), else, for a feature with no
    /// parent, `gene`/`gene_name`/`Name`.
    pub attribute: Option<String>,
    pub strandedness: Strandedness,
    /// Reads (fragments: both mates) below this mapping quality are counted
    /// as low quality, not assigned. Default 10, as htseq-count.
    pub min_mapq: u8,
    /// One name per input file. `None`: the file name without extensions.
    pub sample_names: Option<Vec<String>>,
    /// Where to write the assignment summary as a table (status x sample),
    /// like featureCounts' `.summary`.
    pub summary_output: Option<String>,
}

impl Default for CountOptions {
    fn default() -> Self {
        CountOptions {
            feature_type: "exon".into(),
            attribute: None,
            strandedness: Strandedness::Unstranded,
            min_mapq: 10,
            sample_names: None,
            summary_output: None,
        }
    }
}

/// What happened to the alignments of one file. The first six counters are
/// FRAGMENTS (a read, or a pair of mates counted once) and add up to
/// `fragments`; the skipped ones are alignment RECORDS.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SampleCounts {
    pub sample: String,
    pub path: String,
    /// Alignment lines read (header excluded).
    pub records: usize,
    pub fragments: usize,
    pub assigned: usize,
    /// Overlapping the features of two or more genes.
    pub ambiguous: usize,
    /// Overlapping no feature (introns, intergenic, sequences absent from
    /// the annotation).
    pub no_feature: usize,
    /// NH > 1 in either mate.
    pub multimapping: usize,
    /// Mapping quality below `min_mapq` in either mate.
    pub low_mapq: usize,
    /// Unmapped reads, and pairs with both mates unmapped.
    pub unmapped: usize,
    pub secondary: usize,
    pub supplementary: usize,
    pub duplicate: usize,
    pub qc_fail: usize,
    /// Fragments made of two mates seen together.
    pub pairs: usize,
    /// Paired reads whose mate never came (filtered, or missing from the
    /// file): counted as single reads.
    pub mate_missing: usize,
    /// Mapped records on sequences the annotation does not have.
    pub off_annotation: usize,
}

#[derive(Clone, Debug, Default)]
pub struct GeneCounts {
    pub annotation: String,
    pub feature_type: String,
    /// The gene attribute as used ("gene_id", or the automatic rule).
    pub attribute: String,
    pub strandedness: String,
    pub min_mapq: u8,
    /// Feature lines of the annotation used (one per exon line).
    pub features: usize,
    /// Annotation sequences with features.
    pub sequences: usize,
    /// Genes, in the order they first appear in the annotation.
    pub genes: Vec<String>,
    /// Length of each gene: the union of its features, in bp (what TPM and
    /// FPKM divide by).
    pub lengths: Vec<u64>,
    pub samples: Vec<SampleCounts>,
    /// `counts[sample][gene]`.
    pub counts: Vec<Vec<u64>>,
    pub output: Option<String>,
    pub summary_output: Option<String>,
    pub warnings: Vec<String>,
}

const STATUS: [&str; 12] = [
    "assigned",
    "ambiguous",
    "no_feature",
    "multimapping",
    "low_mapq",
    "unmapped",
    "secondary",
    "supplementary",
    "duplicate",
    "qc_fail",
    "mate_missing",
    "off_annotation",
];

impl SampleCounts {
    fn status(&self) -> [usize; 12] {
        [
            self.assigned,
            self.ambiguous,
            self.no_feature,
            self.multimapping,
            self.low_mapq,
            self.unmapped,
            self.secondary,
            self.supplementary,
            self.duplicate,
            self.qc_fail,
            self.mate_missing,
            self.off_annotation,
        ]
    }

    pub fn json(&self) -> String {
        format!(
            "{{\"sample\":{},\"path\":{},\"records\":{},\"fragments\":{},\"assigned\":{},\"ambiguous\":{},\"no_feature\":{},\"multimapping\":{},\"low_mapq\":{},\"unmapped\":{},\"secondary\":{},\"supplementary\":{},\"duplicate\":{},\"qc_fail\":{},\"pairs\":{},\"mate_missing\":{},\"off_annotation\":{}}}",
            json_texto(&self.sample),
            json_texto(&self.path),
            self.records,
            self.fragments,
            self.assigned,
            self.ambiguous,
            self.no_feature,
            self.multimapping,
            self.low_mapq,
            self.unmapped,
            self.secondary,
            self.supplementary,
            self.duplicate,
            self.qc_fail,
            self.pairs,
            self.mate_missing,
            self.off_annotation
        )
    }
}

impl GeneCounts {
    pub fn json(&self) -> String {
        format!(
            "{{\"annotation\":{},\"feature_type\":{},\"attribute\":{},\"strandedness\":{},\"min_mapq\":{},\"features\":{},\"sequences\":{},\"genes\":{},\"lengths\":{},\"samples\":[{}],\"counts\":[{}],\"output\":{},\"summary_output\":{},\"warnings\":{},\"text\":{}}}",
            json_texto(&self.annotation),
            json_texto(&self.feature_type),
            json_texto(&self.attribute),
            json_texto(&self.strandedness),
            self.min_mapq,
            self.features,
            self.sequences,
            jtextos(&self.genes),
            jinteiros(&self.lengths),
            self.samples.iter().map(SampleCounts::json).collect::<Vec<_>>().join(","),
            self.counts.iter().map(|c| jinteiros(c)).collect::<Vec<_>>().join(","),
            jopt(self.output.as_deref()),
            jopt(self.summary_output.as_deref()),
            jtextos(&self.warnings),
            json_texto(&self.resumo())
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!(
            "Reads per gene: {} genes from {} ({} {} features on {} sequences; gene = {}), {}, min mapping quality {}\n",
            milhar(self.genes.len()),
            self.annotation,
            milhar(self.features),
            self.feature_type,
            self.sequences,
            self.attribute,
            self.strandedness,
            self.min_mapq
        );
        for s in &self.samples {
            o.push_str(&format!(
                "  {}: {} fragments ({} pairs) -- assigned {} ({}), ambiguous {}, no feature {}, multimapping {}, low quality {}, unmapped {}\n",
                s.sample,
                milhar(s.fragments),
                milhar(s.pairs),
                milhar(s.assigned),
                percent(s.assigned, s.fragments),
                milhar(s.ambiguous),
                milhar(s.no_feature),
                milhar(s.multimapping),
                milhar(s.low_mapq),
                milhar(s.unmapped)
            ));
            let pulados = s.secondary + s.supplementary + s.duplicate + s.qc_fail;
            if pulados > 0 || s.mate_missing > 0 {
                o.push_str(&format!(
                    "      records skipped: {} secondary, {} supplementary, {} duplicate, {} QC fail; {} mates without their pair\n",
                    milhar(s.secondary),
                    milhar(s.supplementary),
                    milhar(s.duplicate),
                    milhar(s.qc_fail),
                    milhar(s.mate_missing)
                ));
            }
        }
        if let Some(p) = &self.output {
            o.push_str(&format!("  count matrix written to {p}\n"));
        }
        for w in &self.warnings {
            o.push_str(&format!("  note: {w}\n"));
        }
        o
    }
}

// ---------------------------------------------------------------------
//  The annotation
// ---------------------------------------------------------------------

const NENHUM: u32 = u32::MAX;

/// Strings stored once, referred to by number.
#[derive(Default)]
struct Interna {
    mapa: HashMap<String, u32>,
    nomes: Vec<String>,
}

impl Interna {
    fn id(&mut self, s: &str) -> u32 {
        if let Some(&i) = self.mapa.get(s) {
            return i;
        }
        let i = self.nomes.len() as u32;
        self.nomes.push(s.to_string());
        self.mapa.insert(s.to_string(), i);
        i
    }
}

/// Decodes the `%XX` escapes of GFF3 values.
fn decodifica(v: &str) -> String {
    if !v.contains('%') {
        return v.to_string();
    }
    let b = v.as_bytes();
    let hex = |c: u8| (c as char).to_digit(16).map(|d| d as u8);
    let mut o: Vec<u8> = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'%' && i + 2 < b.len() {
            if let (Some(x), Some(y)) = (hex(b[i + 1]), hex(b[i + 2])) {
                o.push(x * 16 + y);
                i += 3;
                continue;
            }
        }
        o.push(b[i]);
        i += 1;
    }
    String::from_utf8_lossy(&o).to_string()
}

/// The attributes of column 9, GFF3 (`key=value;`) or GTF (`key "value";`).
fn atributos(col9: &str) -> Vec<(&str, &str)> {
    let mut v = Vec::new();
    for parte in col9.split(';') {
        let p = parte.trim();
        if p.is_empty() {
            continue;
        }
        let (eq, sp) = (p.find('='), p.find([' ', '\t']));
        match (eq, sp) {
            (Some(e), s) if s.map_or(true, |s| e < s) => v.push((p[..e].trim(), p[e + 1..].trim())),
            (_, Some(s)) => v.push((p[..s].trim(), p[s + 1..].trim().trim_matches('"'))),
            _ => {}
        }
    }
    v
}

/// The value of attribute `nome` on one line; "GeneID" also reads the NCBI
/// cross-reference (`Dbxref=GeneID:672,...` or `db_xref "GeneID:672"`).
fn valor_de<'a>(at: &[(&'a str, &'a str)], nome: &str) -> Option<&'a str> {
    if let Some((_, v)) = at.iter().find(|(k, _)| *k == nome).or_else(|| at.iter().find(|(k, _)| k.eq_ignore_ascii_case(nome))) {
        return Some(v);
    }
    if nome.eq_ignore_ascii_case("GeneID") {
        for (k, v) in at {
            if k.eq_ignore_ascii_case("Dbxref") || k.eq_ignore_ascii_case("db_xref") {
                for x in v.split(',') {
                    if let Some(id) = x.trim().strip_prefix("GeneID:") {
                        return Some(id);
                    }
                }
            }
        }
    }
    None
}

/// The chromosome of a RefSeq chromosome accession (NC_000001.11 -> 1,
/// NC_000023 -> X, NC_012920 -> MT), any assembly version.
fn cromossomo_refseq(s: &str) -> Option<String> {
    if s.starts_with("NC_012920") || s.starts_with("NC_001807") {
        return Some("MT".to_string());
    }
    let r = s.strip_prefix("NC_0000")?;
    let n: u32 = r.get(..2)?.parse().ok()?;
    if !r[2..].is_empty() && !r[2..].starts_with('.') {
        return None;
    }
    match n {
        1..=22 => Some(n.to_string()),
        23 => Some("X".to_string()),
        24 => Some("Y".to_string()),
        _ => None,
    }
}

/// The key two names of the same sequence share: chr1 == 1 ==
/// NC_000001.11, chrM == MT.
fn chave_seq(nome: &str) -> String {
    match cromossomo_refseq(nome) {
        Some(c) => canonico(&c),
        None => canonico(nome),
    }
}

/// A stretch of one sequence covered by the same set of genes.
#[derive(Clone, Copy, Debug)]
struct Segmento {
    ini: u64,
    fim: u64,
    de: u32,
    ate: u32,
}

/// The genes a read touches, as far as the union rule needs to know.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Conjunto {
    Vazio,
    Um(u32),
    Varios,
}

impl Conjunto {
    fn soma(&mut self, g: u32) {
        *self = match *self {
            Conjunto::Vazio => Conjunto::Um(g),
            Conjunto::Um(h) if h == g => Conjunto::Um(h),
            _ => Conjunto::Varios,
        }
    }
    fn une(self, o: Conjunto) -> Conjunto {
        let mut r = self;
        match o {
            Conjunto::Vazio => {}
            Conjunto::Um(g) => r.soma(g),
            Conjunto::Varios => r = Conjunto::Varios,
        }
        r
    }
}

struct Anotacao {
    caminho: String,
    genes: Vec<String>,
    comprimentos: Vec<u64>,
    seqs: Vec<String>,
    por_nome: HashMap<String, usize>,
    por_chave: HashMap<String, usize>,
    segs: Vec<Vec<Segmento>>,
    /// (gene, strand) of the segments, `Segmento::de..ate`.
    membros: Vec<(u32, u8)>,
    features: usize,
    atributo: String,
    avisos: Vec<String>,
}

impl Anotacao {
    fn seq_de(&self, nome: &str) -> Option<usize> {
        self.por_nome.get(nome).or_else(|| self.por_chave.get(&chave_seq(nome))).copied()
    }

    /// Adds to `c` the genes whose features overlap [ini, fim) on sequence
    /// `seq` (0-based, end excluded); with `fita`, only genes on that strand
    /// (and genes without strand).
    fn junta(&self, seq: usize, ini: u64, fim: u64, fita: Option<u8>, c: &mut Conjunto) {
        let segs = &self.segs[seq];
        let mut i = segs.partition_point(|s| s.fim <= ini);
        while i < segs.len() && segs[i].ini < fim {
            let s = segs[i];
            for &(g, f) in &self.membros[s.de as usize..s.ate as usize] {
                if fita.map_or(true, |x| f == x || f == b'.') {
                    c.soma(g);
                    if *c == Conjunto::Varios {
                        return;
                    }
                }
            }
            i += 1;
        }
    }
}

/// Feature types that are never the parent of another feature (unless they
/// are the type being counted): their lines are skipped without reading the
/// attributes. In a RefSeq annotation the CDS lines alone are 1.6 million.
const FOLHAS: [&str; 11] = [
    "CDS",
    "start_codon",
    "stop_codon",
    "five_prime_UTR",
    "three_prime_UTR",
    "UTR",
    "cDNA_match",
    "match",
    "match_part",
    "Selenocysteine",
    "stop_codon_redefined_as_selenocysteine",
];

/// One feature line as read, before its gene is resolved.
struct Bruto {
    seq: u32,
    ini: u64,
    fim: u64,
    fita: u8,
    valor: u32,
    pai: u32,
    alt: u32,
    linha: usize,
}

fn le_anotacao(caminho: &str, tipo: &str, atributo: Option<&str>) -> Result<Anotacao, String> {
    let auto = atributo.map_or(true, |a| a.trim().is_empty());
    let chave_valor: String = if auto { "gene_id".into() } else { atributo.unwrap().trim().to_string() };
    let mut pool = Interna::default();
    let mut seqs = Interna::default();
    let mut ids: HashMap<u32, (u32, u32)> = HashMap::new();
    let mut aliases: Vec<(String, u32)> = Vec::new();
    let mut brutos: Vec<Bruto> = Vec::new();
    let mut tipos: HashMap<String, usize> = HashMap::new();
    let mut l = Leitor::abre(caminho)?;
    while l.avanca()? {
        let t = l.texto();
        if t.is_empty() {
            continue;
        }
        if t.starts_with("##FASTA") {
            break;
        }
        if t.starts_with('#') {
            continue;
        }
        let c: Vec<&str> = t.splitn(9, '\t').collect();
        if c.len() < 9 {
            return Err(format!(
                "{caminho}, line {}: {} tab-separated column(s); a GFF3/GTF line has 9 (sequence, source, type, start, end, score, strand, frame, attributes). Is this a GFF3 or GTF annotation, with tabs between the columns?",
                l.linha,
                c.len()
            ));
        }
        let e_alvo = c[2].eq_ignore_ascii_case(tipo);
        if !e_alvo {
            if tipos.len() < 1000 || tipos.contains_key(c[2]) {
                *tipos.entry(c[2].to_string()).or_default() += 1;
            }
            // features that are never the parent of another: nothing to keep
            if FOLHAS.iter().any(|f| c[2].eq_ignore_ascii_case(f)) {
                continue;
            }
        }
        let at = atributos(c[8]);
        if !e_alvo {
            if c[2] == "region" && c[3] == "1" {
                if let Some(cr) = valor_de(&at, "chromosome") {
                    aliases.push((cr.to_string(), seqs.id(c[0])));
                }
            }
            if let Some(id) = valor_de(&at, "ID") {
                let id = pool.id(id);
                if !ids.contains_key(&id) {
                    let pai = valor_de(&at, "Parent").map(|p| pool.id(p.split(',').next().unwrap_or(p))).unwrap_or(NENHUM);
                    let valor = valor_de(&at, &chave_valor).map(|v| pool.id(&decodifica(v))).unwrap_or(NENHUM);
                    ids.insert(id, (pai, valor));
                }
            }
            continue;
        }
        let (Ok(ini), Ok(fim)) = (c[3].trim().parse::<u64>(), c[4].trim().parse::<u64>()) else {
            return Err(format!("{caminho}, line {}: start \"{}\" and end \"{}\" must be whole numbers", l.linha, c[3], c[4]));
        };
        if ini == 0 || fim < ini {
            return Err(format!(
                "{caminho}, line {}: a feature from {ini} to {fim}; GFF/GTF coordinates are 1-based with start <= end",
                l.linha
            ));
        }
        let fita = match c[6] {
            "+" => b'+',
            "-" => b'-',
            _ => b'.',
        };
        let valor = valor_de(&at, &chave_valor).map(|v| pool.id(&decodifica(v))).unwrap_or(NENHUM);
        let pai = valor_de(&at, "Parent").map(|p| pool.id(p.split(',').next().unwrap_or(p))).unwrap_or(NENHUM);
        // the name fallback is used only for a feature with no gene_id and no parent
        let alt = if auto && valor == NENHUM && pai == NENHUM {
            ["gene", "gene_name", "Name"].iter().find_map(|k| valor_de(&at, k)).map(|v| pool.id(&decodifica(v))).unwrap_or(NENHUM)
        } else {
            NENHUM
        };
        brutos.push(Bruto { seq: seqs.id(c[0]), ini: ini - 1, fim, fita, valor, pai, alt, linha: l.linha });
    }
    if brutos.is_empty() {
        let mut v: Vec<(String, usize)> = tipos.into_iter().collect();
        v.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        let lista: Vec<String> = v.iter().take(12).map(|(t, n)| format!("{t} ({n})")).collect();
        return Err(if lista.is_empty() {
            format!("{caminho} has no annotation lines: is it a GFF3 or GTF file?")
        } else {
            format!(
                "{caminho} has no feature of type \"{tipo}\" (column 3). Types present: {}. Choose one of them as the feature type.",
                lista.join(", ")
            )
        });
    }

    // the gene of each feature
    let mut cache_raiz: HashMap<u32, u32> = HashMap::new();
    let mut gene_de: HashMap<u32, u32> = HashMap::new();
    let mut genes: Vec<String> = Vec::new();
    let mut exons: Vec<(u32, u64, u64, u32, u8)> = Vec::with_capacity(brutos.len());
    let mut sem_gene = 0usize;
    let mut primeira_sem: Option<usize> = None;
    for b in &brutos {
        let mut rotulo = b.valor;
        if rotulo == NENHUM {
            let (mut p, mut raiz, mut passos) = (b.pai, NENHUM, 0);
            while p != NENHUM && passos < 64 {
                raiz = p;
                match ids.get(&p) {
                    Some(&(pp, v)) => {
                        if v != NENHUM {
                            rotulo = v;
                            break;
                        }
                        p = pp;
                    }
                    None => break,
                }
                passos += 1;
            }
            if rotulo == NENHUM && auto {
                if raiz != NENHUM {
                    rotulo = *cache_raiz.entry(raiz).or_insert_with(|| {
                        let n = pool.nomes[raiz as usize].clone();
                        let curto = n.strip_prefix("gene-").or_else(|| n.strip_prefix("gene:")).unwrap_or(&n).to_string();
                        pool.id(&curto)
                    });
                } else {
                    rotulo = b.alt;
                }
            }
        }
        if rotulo == NENHUM {
            sem_gene += 1;
            primeira_sem.get_or_insert(b.linha);
            continue;
        }
        let g = *gene_de.entry(rotulo).or_insert_with(|| {
            genes.push(pool.nomes[rotulo as usize].clone());
            (genes.len() - 1) as u32
        });
        exons.push((b.seq, b.ini, b.fim, g, b.fita));
    }
    let atributo_txt = if auto {
        "gene_id (else the top-level parent's ID, else gene/Name)".to_string()
    } else {
        chave_valor.clone()
    };
    if exons.is_empty() {
        return Err(format!(
            "{caminho}: none of the {} \"{tipo}\" features has a gene: no `{chave_valor}` attribute on the feature or on its parents (first such feature at line {}). Choose the attribute that names the gene (gene_id, gene, gene_name, GeneID...).",
            brutos.len(),
            primeira_sem.unwrap_or(0)
        ));
    }
    let mut avisos = Vec::new();
    if sem_gene > 0 {
        avisos.push(format!(
            "{} {tipo} feature(s) of {caminho} have no gene ({atributo_txt}) and were left out; the first is at line {}",
            sem_gene,
            primeira_sem.unwrap_or(0)
        ));
    }
    drop(brutos);
    drop(ids);

    // gene lengths: the union of each gene's features
    let mut por_gene = exons.clone();
    multi_merge_sort_stable_by(&mut por_gene, |a, b| (a.3, a.0, a.1).cmp(&(b.3, b.0, b.1)));
    let mut comprimentos = vec![0u64; genes.len()];
    let mut i = 0;
    while i < por_gene.len() {
        let (s, g) = (por_gene[i].0, por_gene[i].3);
        let (mut a, mut z) = (por_gene[i].1, por_gene[i].2);
        let mut j = i + 1;
        while j < por_gene.len() && por_gene[j].3 == g && por_gene[j].0 == s {
            if por_gene[j].1 > z {
                comprimentos[g as usize] += z - a;
                a = por_gene[j].1;
                z = por_gene[j].2;
            } else {
                z = z.max(por_gene[j].2);
            }
            j += 1;
        }
        comprimentos[g as usize] += z - a;
        i = j;
    }
    drop(por_gene);

    // segments of constant gene set, per sequence
    multi_merge_sort_stable_by(&mut exons, |a, b| (a.0, a.1, a.2).cmp(&(b.0, b.1, b.2)));
    let n_seq = seqs.nomes.len();
    let mut segs: Vec<Vec<Segmento>> = vec![Vec::new(); n_seq];
    let mut membros: Vec<(u32, u8)> = Vec::new();
    let mut i = 0;
    while i < exons.len() {
        let s = exons[i].0;
        let mut j = i;
        while j < exons.len() && exons[j].0 == s {
            j += 1;
        }
        let bloco = &exons[i..j];
        let mut limites: Vec<u64> = bloco.iter().flat_map(|e| [e.1, e.2]).collect();
        limites.sort_unstable();
        limites.dedup();
        let mut ativos: Vec<(u64, u32, u8)> = Vec::new();
        let mut k = 0;
        let lista = &mut segs[s as usize];
        for w in 0..limites.len().saturating_sub(1) {
            let (b0, b1) = (limites[w], limites[w + 1]);
            while k < bloco.len() && bloco[k].1 == b0 {
                ativos.push((bloco[k].2, bloco[k].3, bloco[k].4));
                k += 1;
            }
            ativos.retain(|a| a.0 > b0);
            if ativos.is_empty() {
                continue;
            }
            let mut m: Vec<(u32, u8)> = ativos.iter().map(|a| (a.1, a.2)).collect();
            m.sort_unstable();
            m.dedup();
            if let Some(u) = lista.last_mut() {
                if u.fim == b0 && membros[u.de as usize..u.ate as usize] == m[..] {
                    u.fim = b1;
                    continue;
                }
            }
            let de = membros.len() as u32;
            membros.extend_from_slice(&m);
            lista.push(Segmento { ini: b0, fim: b1, de, ate: membros.len() as u32 });
        }
        i = j;
    }
    let mut por_nome = HashMap::new();
    let mut por_chave = HashMap::new();
    for (i, n) in seqs.nomes.iter().enumerate() {
        por_nome.insert(n.clone(), i);
    }
    for (i, n) in seqs.nomes.iter().enumerate() {
        por_chave.entry(chave_seq(n)).or_insert(i);
    }
    for (cr, i) in aliases {
        por_chave.entry(chave_seq(&cr)).or_insert(i as usize);
    }
    let features = exons.len();
    Ok(Anotacao {
        caminho: caminho.to_string(),
        genes,
        comprimentos,
        seqs: seqs.nomes,
        por_nome,
        por_chave,
        segs,
        membros,
        features,
        atributo: atributo_txt,
        avisos,
    })
}

// ---------------------------------------------------------------------
//  Counting one file
// ---------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
struct Estado {
    genes: Conjunto,
    multi: bool,
    baixa: bool,
}

fn finaliza(e: Estado, s: &mut SampleCounts, cont: &mut [u64]) {
    s.fragments += 1;
    if e.multi {
        s.multimapping += 1;
    } else if e.baixa {
        s.low_mapq += 1;
    } else {
        match e.genes {
            Conjunto::Vazio => s.no_feature += 1,
            Conjunto::Varios => s.ambiguous += 1,
            Conjunto::Um(g) => {
                s.assigned += 1;
                cont[g as usize] += 1;
            }
        }
    }
}

/// The aligned blocks of a CIGAR (M, = and X runs; N and D end a block), as
/// 0-based half-open reference intervals. `false`: the CIGAR is invalid.
fn blocos_cigar(pos1: u64, cigar: &str, blocos: &mut Vec<(u64, u64)>) -> bool {
    blocos.clear();
    let mut r = pos1 - 1;
    let mut n: u64 = 0;
    let mut tem = false;
    let mut aberto: Option<u64> = None;
    for b in cigar.bytes() {
        if b.is_ascii_digit() {
            n = n * 10 + u64::from(b - b'0');
            if n > 1 << 32 {
                return false;
            }
            tem = true;
            continue;
        }
        if !tem {
            return false;
        }
        match b {
            b'M' | b'=' | b'X' => {
                if aberto.is_none() {
                    aberto = Some(r);
                }
                r += n;
            }
            b'D' | b'N' => {
                if let Some(a) = aberto.take() {
                    if r > a {
                        blocos.push((a, r));
                    }
                }
                r += n;
            }
            b'I' | b'S' | b'H' | b'P' => {}
            _ => return false,
        }
        n = 0;
        tem = false;
    }
    if tem {
        return false;
    }
    if let Some(a) = aberto {
        if r > a {
            blocos.push((a, r));
        }
    }
    true
}

fn tag_nh(opcionais: &str) -> Option<u64> {
    opcionais.split('\t').find_map(|c| c.strip_prefix("NH:i:")).and_then(|v| v.trim().parse().ok())
}

fn lista_curta(v: &[String]) -> String {
    let mut x: Vec<&str> = v.iter().take(6).map(String::as_str).collect();
    if v.len() > 6 {
        x.push("...");
    }
    x.join(", ")
}

fn conta_arquivo(caminho: &str, amostra: &str, an: &Anotacao, opt: &CountOptions) -> Result<(SampleCounts, Vec<u64>), String> {
    let mut l = Leitor::abre(caminho)?;
    let mut s = SampleCounts { sample: amostra.to_string(), path: caminho.to_string(), ..Default::default() };
    let mut cont = vec![0u64; an.genes.len()];
    let mut pendentes: HashMap<String, Estado> = HashMap::new();
    let mut sq: Vec<String> = Vec::new();
    let mut conferido = false;
    let mut cache: (String, Option<usize>) = (String::new(), None);
    let mut blocos: Vec<(u64, u64)> = Vec::new();
    let confere = |sq: &[String]| -> Result<(), String> {
        if !sq.is_empty() && !sq.iter().any(|n| an.seq_de(n).is_some()) {
            return Err(format!(
                "None of the {} sequences in the header of {caminho} ({}) is in the annotation {} (its sequences: {}). The alignments and the annotation must be of the same genome; if only the naming differs beyond chr1/1 and RefSeq accessions, rename one of them (normalize_chromosomes).",
                sq.len(),
                lista_curta(sq),
                an.caminho,
                lista_curta(&an.seqs)
            ));
        }
        Ok(())
    };
    while l.avanca()? {
        let t = l.texto();
        if t.is_empty() {
            continue;
        }
        if t.starts_with('@') {
            if let Some(sn) = t.strip_prefix("@SQ").and_then(|r| r.split('\t').find_map(|c| c.strip_prefix("SN:"))) {
                sq.push(sn.to_string());
            }
            continue;
        }
        if !conferido {
            confere(&sq)?;
            conferido = true;
        }
        let c: Vec<&str> = t.splitn(12, '\t').collect();
        if c.len() < 11 {
            return Err(format!(
                "{caminho}, line {}: {} tab-separated column(s); a SAM alignment line has at least 11. Is this a SAM or BAM file?",
                l.linha,
                c.len()
            ));
        }
        let Ok(flag) = c[1].parse::<u16>() else {
            return Err(format!("{caminho}, line {}: FLAG \"{}\" is not a number between 0 and 65535", l.linha, c[1]));
        };
        s.records += 1;
        if flag & 0x100 != 0 {
            s.secondary += 1;
            continue;
        }
        if flag & 0x800 != 0 {
            s.supplementary += 1;
            continue;
        }
        if flag & 0x200 != 0 {
            s.qc_fail += 1;
            continue;
        }
        if flag & 0x400 != 0 {
            s.duplicate += 1;
            continue;
        }
        let par = flag & 0x1 != 0;
        let mate_fora = flag & 0x8 != 0;
        let pos: u64 = c[3].parse().map_err(|_| format!("{caminho}, line {}: POS \"{}\" is not a whole number", l.linha, c[3]))?;
        if flag & 0x4 != 0 || c[2] == "*" || pos == 0 || c[5] == "*" {
            if par && !mate_fora {
                // the mapped mate stands for the fragment
                continue;
            }
            if !par || flag & 0x80 == 0 {
                s.unmapped += 1;
                s.fragments += 1;
            }
            continue;
        }
        let mapq: u32 = c[4].parse().map_err(|_| format!("{caminho}, line {}: MAPQ \"{}\" is not a whole number", l.linha, c[4]))?;
        let multi = c.get(11).and_then(|o| tag_nh(o)).map_or(false, |nh| nh > 1);
        let baixa = mapq < u32::from(opt.min_mapq);
        if !blocos_cigar(pos, c[5], &mut blocos) {
            return Err(format!(
                "{caminho}, line {}: read {} has an invalid CIGAR \"{}\" (operations are M I D N S H P = X, each after a length)",
                l.linha, c[0], c[5]
            ));
        }
        if c[2] != cache.0 {
            cache = (c[2].to_string(), an.seq_de(c[2]));
        }
        let mut genes = Conjunto::Vazio;
        match cache.1 {
            None => s.off_annotation += 1,
            Some(sq_i) if !multi && !baixa => {
                let reversa = (flag & 0x10 != 0) ^ (par && flag & 0x80 != 0);
                let frag = if reversa { b'-' } else { b'+' };
                let fita = match opt.strandedness {
                    Strandedness::Unstranded => None,
                    Strandedness::Forward => Some(frag),
                    Strandedness::Reverse => Some(if frag == b'+' { b'-' } else { b'+' }),
                };
                for &(a, z) in &blocos {
                    an.junta(sq_i, a, z, fita, &mut genes);
                }
            }
            Some(_) => {}
        }
        let e = Estado { genes, multi, baixa };
        if par && !mate_fora {
            match pendentes.remove(c[0]) {
                Some(o) => {
                    s.pairs += 1;
                    finaliza(Estado { genes: o.genes.une(e.genes), multi: o.multi || e.multi, baixa: o.baixa || e.baixa }, &mut s, &mut cont);
                }
                None => {
                    pendentes.insert(c[0].to_string(), e);
                }
            }
        } else {
            finaliza(e, &mut s, &mut cont);
        }
    }
    if !conferido {
        confere(&sq)?;
    }
    for (_, e) in pendentes.drain() {
        s.mate_missing += 1;
        finaliza(e, &mut s, &mut cont);
    }
    Ok((s, cont))
}

fn nome_de_amostra(caminho: &str) -> String {
    let base = std::path::Path::new(caminho).file_name().map(|f| f.to_string_lossy().to_string()).unwrap_or_else(|| caminho.to_string());
    let mut b = base.as_str();
    for suf in [".gz", ".bgz", ".lz4"] {
        if b.to_ascii_lowercase().ends_with(suf) {
            b = &b[..b.len() - suf.len()];
        }
    }
    for suf in [".bam", ".sam", ".cram"] {
        if b.to_ascii_lowercase().ends_with(suf) {
            b = &b[..b.len() - suf.len()];
            break;
        }
    }
    b.to_string()
}

/// Reads per gene, one column per alignment file (BAM, SAM, `.sam.gz`).
///
/// The annotation is a GFF3 or GTF file (`.gz` accepted); features of type
/// `opt.feature_type` are grouped into genes by `opt.attribute` (see
/// `CountOptions`). Sequence names are matched exactly or by their canonical
/// form (chr1 == 1, chrM == MT, NC_000001.11 == 1).
///
/// Union mode (htseq-count `-m union`, featureCounts' default): a fragment is
/// assigned when its aligned bases overlap features of exactly one gene. The
/// two mates of a pair (FLAG 0x1) are one fragment, matched by name, in any
/// file order; a mate whose pair never comes is counted alone. A fragment is
/// multimapping when either mate has NH > 1, low quality when either mate is
/// below `min_mapq` (htseq-count's rule; featureCounts asks only one). In a
/// stranded library the fragment's strand is read 1's (read 2 flipped);
/// features without strand ('.') match both.
///
/// `output` (if any) gets the matrix: `gene` and one column per sample, the
/// genes in annotation order. The files are counted in parallel.
pub fn count_reads(paths: &[String], annotation: &str, output: Option<&str>, opt: &CountOptions) -> Result<GeneCounts, String> {
    if paths.is_empty() {
        return Err("No alignment file given: pass one BAM or SAM file per sample".into());
    }
    let nomes: Vec<String> = match &opt.sample_names {
        Some(v) => {
            if v.len() != paths.len() {
                return Err(format!(
                    "{} sample name(s) for {} alignment file(s): give one name per file, in the same order",
                    v.len(),
                    paths.len()
                ));
            }
            v.iter().map(|s| s.trim().to_string()).collect()
        }
        None => paths.iter().map(|p| nome_de_amostra(p)).collect(),
    };
    for (i, n) in nomes.iter().enumerate() {
        if n.is_empty() || n.contains(['\t', '\n']) {
            return Err(format!("Sample name \"{n}\" (file {}) is empty or has a tab or line break", paths[i]));
        }
        if let Some(j) = nomes[..i].iter().position(|m| m == n) {
            return Err(format!(
                "Two files give the same sample name \"{n}\" ({} and {}): pass sample_names with one distinct name per file",
                paths[j], paths[i]
            ));
        }
    }
    for (i, p) in paths.iter().enumerate() {
        if let Some(j) = paths[..i].iter().position(|q| mesmo_arquivo(p, q)) {
            return Err(format!("The same file is given twice ({} and {}): each file is one sample", paths[j], p));
        }
    }
    for o in [output, opt.summary_output.as_deref()].into_iter().flatten() {
        recusa_mesmo_arquivo(annotation, o)?;
        for p in paths {
            recusa_mesmo_arquivo(p, o)?;
        }
    }
    let tipo = if opt.feature_type.trim().is_empty() { "exon" } else { opt.feature_type.trim() };
    let an = le_anotacao(annotation, tipo, opt.attribute.as_deref())?;
    let res: Vec<Result<(SampleCounts, Vec<u64>), String>> =
        paths.par_iter().zip(nomes.par_iter()).map(|(p, n)| conta_arquivo(p, n, &an, opt)).collect();
    let mut r = GeneCounts {
        annotation: annotation.to_string(),
        feature_type: tipo.to_string(),
        attribute: an.atributo.clone(),
        strandedness: opt.strandedness.as_str().to_string(),
        min_mapq: opt.min_mapq,
        features: an.features,
        sequences: an.segs.iter().filter(|s| !s.is_empty()).count(),
        genes: an.genes.clone(),
        lengths: an.comprimentos.clone(),
        warnings: an.avisos.clone(),
        ..Default::default()
    };
    for x in res {
        let (s, c) = x?;
        if s.off_annotation > 0 {
            r.warnings.push(format!(
                "{}: {} mapped record(s) on sequences absent from the annotation, counted as no feature",
                s.sample,
                milhar(s.off_annotation)
            ));
        }
        if s.records > 0 && s.assigned == 0 {
            r.warnings.push(format!(
                "{}: no fragment was assigned to a gene. Check the strandedness (try \"unstranded\"), the feature type and that the annotation is of the genome the reads were aligned to",
                s.sample
            ));
        }
        r.samples.push(s);
        r.counts.push(c);
    }
    if let Some(saida) = output {
        let mut w = Escritor::cria(saida)?;
        let mut guarda = Saida { caminho: saida.to_string(), armada: true };
        let erro = |e: String| format!("Error writing {saida}: {e}");
        w.linha(&format!("gene\t{}", nomes.join("\t"))).map_err(erro)?;
        let mut linha = String::new();
        for (g, nome) in r.genes.iter().enumerate() {
            linha.clear();
            linha.push_str(nome);
            for c in &r.counts {
                linha.push('\t');
                linha.push_str(&c[g].to_string());
            }
            w.linha(&linha).map_err(erro)?;
        }
        w.termina().map_err(erro)?;
        guarda.armada = false;
        r.output = Some(saida.to_string());
    }
    if let Some(saida) = &opt.summary_output {
        let mut w = Escritor::cria(saida)?;
        let mut guarda = Saida { caminho: saida.clone(), armada: true };
        let erro = |e: String| format!("Error writing {saida}: {e}");
        w.linha(&format!("status\t{}", nomes.join("\t"))).map_err(erro)?;
        let st: Vec<[usize; 12]> = r.samples.iter().map(SampleCounts::status).collect();
        for (k, nome) in STATUS.iter().enumerate() {
            let v: Vec<String> = st.iter().map(|s| s[k].to_string()).collect();
            w.linha(&format!("{nome}\t{}", v.join("\t"))).map_err(erro)?;
        }
        w.termina().map_err(erro)?;
        guarda.armada = false;
        r.summary_output = Some(saida.clone());
    }
    Ok(r)
}

// =====================================================================
//  2) SPECIAL FUNCTIONS AND SMALL NUMERICS
// =====================================================================

const LN_RAIZ_2PI: f64 = 0.918_938_533_204_672_8;

/// lnGamma(x) - [(x - 1/2) ln x - x + ln(2 pi)/2], for x >= 10 (Stirling).
fn stirling(x: f64) -> f64 {
    let z = 1.0 / (x * x);
    (1.0 / 12.0 - z * (1.0 / 360.0 - z * (1.0 / 1260.0 - z * (1.0 / 1680.0 - z * (1.0 / 1188.0 - z * (691.0 / 360360.0 - z / 156.0))))))
        / x
}

/// ln Gamma(x), x > 0 (relative error ~1e-15).
pub fn ln_gamma(x: f64) -> f64 {
    if !(x > 0.0) {
        return f64::NAN;
    }
    let (mut x, mut ajuste) = (x, 0.0);
    if x < 10.0 {
        let mut prod = 1.0;
        while x < 10.0 {
            prod *= x;
            x += 1.0;
        }
        ajuste = -prod.ln();
    }
    ajuste + (x - 0.5) * x.ln() - x + LN_RAIZ_2PI + stirling(x)
}

/// ln Gamma(s + y) - ln Gamma(s) - y ln s, accurate when s is huge (a
/// dispersion near zero) where the plain difference loses every digit.
fn ln_gamma_razao(y: f64, s: f64) -> f64 {
    if y == 0.0 {
        return 0.0;
    }
    if s >= 10.0 {
        let r = (y / s).ln_1p();
        (s + y - 0.5) * r - y + (stirling(s + y) - stirling(s))
    } else {
        ln_gamma(s + y) - ln_gamma(s) - y * s.ln()
    }
}

/// psi(x) - ln x + 1/(2x), for x >= 10.
fn serie_psi(x: f64) -> f64 {
    let z = 1.0 / (x * x);
    z * (1.0 / 12.0 - z * (1.0 / 120.0 - z * (1.0 / 252.0 - z * (1.0 / 240.0 - z * (1.0 / 132.0 - z * (691.0 / 32760.0 - z / 12.0))))))
}

/// The digamma function psi(x) = d/dx ln Gamma(x), x > 0.
pub fn digamma(x: f64) -> f64 {
    if !(x > 0.0) {
        return f64::NAN;
    }
    let (mut x, mut acc) = (x, 0.0);
    while x < 10.0 {
        acc -= 1.0 / x;
        x += 1.0;
    }
    acc + x.ln() - 0.5 / x - serie_psi(x)
}

/// psi(s + y) - psi(s), accurate for huge s.
fn digamma_dif(y: f64, s: f64) -> f64 {
    if y == 0.0 {
        return 0.0;
    }
    if s >= 10.0 {
        let t = s + y;
        (y / s).ln_1p() + 0.5 * y / (s * t) - (serie_psi(t) - serie_psi(s))
    } else {
        digamma(s + y) - digamma(s)
    }
}

/// The trigamma function psi'(x), x > 0.
pub fn trigamma(x: f64) -> f64 {
    if !(x > 0.0) {
        return f64::NAN;
    }
    let (mut x, mut acc) = (x, 0.0);
    while x < 10.0 {
        acc += 1.0 / (x * x);
        x += 1.0;
    }
    let z = 1.0 / (x * x);
    acc + 1.0 / x
        + 0.5 * z
        + z / x * (1.0 / 6.0 - z * (1.0 / 30.0 - z * (1.0 / 42.0 - z * (1.0 / 30.0 - z * (5.0 / 66.0 - z * (691.0 / 2730.0 - z * 7.0 / 6.0))))))
}

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
        com_del(x.abs(), temp)
    } else if x > 0.0 {
        (1.0, 0.0)
    } else {
        (0.0, 1.0)
    }
}

/// P(Z <= x) for a standard normal Z.
pub fn normal_cdf(x: f64) -> f64 {
    normal_caudas(x).0
}

/// Two-sided p-value of a Wald statistic: 2 P(Z > |z|).
pub fn two_sided_p(z: f64) -> f64 {
    if z.is_nan() {
        return f64::NAN;
    }
    (2.0 * normal_caudas(-z.abs()).0).min(1.0)
}

/// Regularized lower incomplete gamma P(a, x).
fn gama_inc(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x.is_infinite() {
        return 1.0;
    }
    let ln_pre = -x + a * x.ln() - ln_gamma(a);
    if x < a + 1.0 {
        let (mut ap, mut soma, mut del) = (a, 1.0 / a, 1.0 / a);
        for _ in 0..10_000 {
            ap += 1.0;
            del *= x / ap;
            soma += del;
            if del.abs() < soma.abs() * 1e-16 {
                break;
            }
        }
        (soma.ln() + ln_pre).exp().min(1.0)
    } else {
        // continued fraction for Q (Lentz)
        let tiny = 1e-300;
        let mut b = x + 1.0 - a;
        let mut c = 1.0 / tiny;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..10_000 {
            let an = -(i as f64) * (i as f64 - a);
            b += 2.0;
            d = an * d + b;
            if d.abs() < tiny {
                d = tiny;
            }
            c = b + an / c;
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
        (1.0 - (ln_pre + h.ln()).exp()).max(0.0)
    }
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
    for m in 1..10_000 {
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

/// Regularized incomplete beta I_x(a, b).
fn beta_inc(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let ln_bt = ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln();
    if x < (a + 1.0) / (a + b + 2.0) {
        (ln_bt.exp() * beta_cf(a, b, x) / a).clamp(0.0, 1.0)
    } else {
        (1.0 - ln_bt.exp() * beta_cf(b, a, 1.0 - x) / b).clamp(0.0, 1.0)
    }
}

/// P(F <= f) for an F distribution with (d1, d2) degrees of freedom.
pub fn f_cdf(f: f64, d1: f64, d2: f64) -> f64 {
    if f <= 0.0 {
        return 0.0;
    }
    beta_inc(d1 / 2.0, d2 / 2.0, d1 * f / (d1 * f + d2))
}

/// The q quantile of the F distribution (R's qf(q, d1, d2)).
pub fn f_quantile(q: f64, d1: f64, d2: f64) -> f64 {
    if !(q > 0.0 && q < 1.0) || !(d1 > 0.0) || !(d2 > 0.0) {
        return f64::NAN;
    }
    let mut hi = 1.0;
    while f_cdf(hi, d1, d2) < q && hi < 1e300 {
        hi *= 2.0;
    }
    let mut lo = 0.0;
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f_cdf(mid, d1, d2) < q {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo <= 1e-14 * hi {
            break;
        }
    }
    0.5 * (lo + hi)
}

/// Median as R computes it (mean of the two middle values when even).
/// NaN values are ignored; NaN for an empty input.
pub fn median(v: &[f64]) -> f64 {
    let mut w: Vec<f64> = v.iter().copied().filter(|x| !x.is_nan()).collect();
    if w.is_empty() {
        return f64::NAN;
    }
    w.sort_by(|a, b| a.total_cmp(b));
    let n = w.len();
    if n % 2 == 1 {
        w[n / 2]
    } else {
        0.5 * (w[n / 2 - 1] + w[n / 2])
    }
}

/// R's mad(): 1.4826 x median absolute deviation from the median.
pub fn mad(v: &[f64]) -> f64 {
    let m = median(v);
    let d: Vec<f64> = v.iter().map(|x| (x - m).abs()).collect();
    1.4826 * median(&d)
}

/// R's quantile(x, p) of type 7, on sorted data.
fn quantil7(ordenados: &[f64], p: f64) -> f64 {
    let n = ordenados.len();
    if n == 0 {
        return f64::NAN;
    }
    let h = (n - 1) as f64 * p;
    let lo = h.floor() as usize;
    let hi = (lo + 1).min(n - 1);
    ordenados[lo] + (h - lo as f64) * (ordenados[hi] - ordenados[lo])
}

/// R's mean(x, trim = t).
fn media_aparada(v: &[f64], t: f64) -> f64 {
    let n = v.len();
    if n == 0 {
        return f64::NAN;
    }
    if t >= 0.5 {
        return median(v);
    }
    let mut w = v.to_vec();
    w.sort_by(|a, b| a.total_cmp(b));
    let lo = (n as f64 * t).floor() as usize;
    let hi = n - lo;
    w[lo..hi].iter().sum::<f64>() / (hi - lo) as f64
}

/// Benjamini-Hochberg adjusted p-values (R's p.adjust(p, "BH")). NaN entries
/// stay NaN and do not count in the number of tests.
pub fn benjamini_hochberg(p: &[f64]) -> Vec<f64> {
    let mut idx: Vec<u32> = (0..p.len() as u32).filter(|&i| !p[i as usize].is_nan()).collect();
    let n = idx.len();
    let mut r = vec![f64::NAN; p.len()];
    if n == 0 {
        return r;
    }
    // decreasing p
    multi_merge_sort_stable_by(&mut idx, |a, b| p[*b as usize].total_cmp(&p[*a as usize]));
    let mut cmin = f64::INFINITY;
    for (k, &i) in idx.iter().enumerate() {
        let rank = (n - k) as f64;
        let v = (n as f64 / rank * p[i as usize]).min(cmin);
        cmin = v;
        r[i as usize] = v.min(1.0);
    }
    r
}

/// R's lowess(x, y, f, iter = 3, delta = 0.01 * diff(range(x))), on x
/// already in increasing order. Returns the fitted values.
pub fn lowess(x: &[f64], y: &[f64], f: f64) -> Vec<f64> {
    let n = x.len();
    let mut ys = vec![0.0; n];
    if n < 2 {
        if n == 1 {
            ys[0] = y[0];
        }
        return ys;
    }
    let delta = 0.01 * (x[n - 1] - x[0]);
    let nsteps = 3;
    let mut rw = vec![0.0; n];
    let mut res = vec![0.0; n];
    let ns = 2.max(n.min((f * n as f64 + 1e-7) as usize));
    let fcube = |v: f64| v * v * v;
    // `lowest` of R's clowess, 0-based: fitted value at xs from points
    // nleft..=nright (and ties to the right).
    let lowest = |xs: f64, nleft: usize, nright: usize, userw: bool, rw: &[f64], w: &mut [f64]| -> Option<f64> {
        let range = x[n - 1] - x[0];
        let h = (xs - x[nleft]).max(x[nright] - xs);
        let h9 = 0.999 * h;
        let h1 = 0.001 * h;
        let mut a = 0.0;
        let mut j = nleft;
        while j < n {
            w[j] = 0.0;
            let r = (x[j] - xs).abs();
            if r <= h9 {
                if r <= h1 {
                    w[j] = 1.0;
                } else {
                    w[j] = fcube(1.0 - fcube(r / h));
                }
                if userw {
                    w[j] *= rw[j];
                }
                a += w[j];
            } else if x[j] > xs {
                break;
            }
            j += 1;
        }
        let nrt = j - 1;
        if a <= 0.0 {
            return None;
        }
        for k in nleft..=nrt {
            w[k] /= a;
        }
        if h > 0.0 {
            let mut a = 0.0;
            for k in nleft..=nrt {
                a += w[k] * x[k];
            }
            let mut b = xs - a;
            let mut c = 0.0;
            for k in nleft..=nrt {
                c += w[k] * (x[k] - a) * (x[k] - a);
            }
            if c.sqrt() > 0.001 * range {
                b /= c;
                for k in nleft..=nrt {
                    w[k] *= b * (x[k] - a) + 1.0;
                }
            }
        }
        let mut ysv = 0.0;
        for k in nleft..=nrt {
            ysv += w[k] * y[k];
        }
        Some(ysv)
    };
    let mut iter = 1;
    while iter <= nsteps + 1 {
        let mut nleft = 0usize;
        let mut nright = ns - 1;
        let mut last: isize = -1;
        let mut i = 0usize;
        loop {
            if nright < n - 1 {
                let d1 = x[i] - x[nleft];
                let d2 = x[nright + 1] - x[i];
                if d1 > d2 {
                    nleft += 1;
                    nright += 1;
                    continue;
                }
            }
            let mut wbuf = std::mem::take(&mut res);
            match lowest(x[i], nleft, nright, iter > 1, &rw, &mut wbuf) {
                Some(v) => ys[i] = v,
                None => ys[i] = y[i],
            }
            res = wbuf;
            if last < i as isize - 1 {
                let l = last as usize;
                let denom = x[i] - x[l];
                for j in l + 1..i {
                    let alpha = (x[j] - x[l]) / denom;
                    ys[j] = alpha * ys[i] + (1.0 - alpha) * ys[l];
                }
            }
            last = i as isize;
            let lu = last as usize;
            let cut = x[lu] + delta;
            let mut ii = lu + 1;
            let mut lu2 = lu;
            while ii < n {
                if x[ii] > cut {
                    break;
                }
                if x[ii] == x[lu2] {
                    ys[ii] = ys[lu2];
                    lu2 = ii;
                }
                ii += 1;
            }
            last = lu2 as isize;
            i = (lu2 + 1).max(ii - 1);
            if lu2 >= n - 1 {
                break;
            }
        }
        for k in 0..n {
            res[k] = y[k] - ys[k];
        }
        let sc = res.iter().map(|r| r.abs()).sum::<f64>() / n as f64;
        if iter > nsteps {
            break;
        }
        for k in 0..n {
            rw[k] = res[k].abs();
        }
        let mut ord = rw.clone();
        ord.sort_by(|a, b| a.total_cmp(b));
        let m1 = n / 2;
        let cmad = if n % 2 == 0 { 3.0 * (ord[m1] + ord[n - m1 - 1]) } else { 6.0 * ord[m1] };
        if cmad < 1e-7 * sc {
            break;
        }
        let c9 = 0.999 * cmad;
        let c1 = 0.001 * cmad;
        for k in 0..n {
            let r = res[k].abs();
            rw[k] = if r <= c1 {
                1.0
            } else if r <= c9 {
                let u = r / cmad;
                (1.0 - u * u) * (1.0 - u * u)
            } else {
                0.0
            };
        }
        iter += 1;
    }
    ys
}

// --- small dense linear algebra (p x p, row-major) ---------------------

fn cholesky(a: &[f64], p: usize) -> Option<Vec<f64>> {
    let mut l = vec![0.0; p * p];
    for i in 0..p {
        for j in 0..=i {
            let mut s = a[i * p + j];
            for k in 0..j {
                s -= l[i * p + k] * l[j * p + k];
            }
            if i == j {
                if !(s > 0.0) || !s.is_finite() {
                    return None;
                }
                l[i * p + i] = s.sqrt();
            } else {
                l[i * p + j] = s / l[j * p + j];
            }
        }
    }
    Some(l)
}

fn chol_resolve(l: &[f64], p: usize, b: &[f64]) -> Vec<f64> {
    let mut z = b.to_vec();
    for i in 0..p {
        let mut s = z[i];
        for k in 0..i {
            s -= l[i * p + k] * z[k];
        }
        z[i] = s / l[i * p + i];
    }
    for i in (0..p).rev() {
        let mut s = z[i];
        for k in i + 1..p {
            s -= l[k * p + i] * z[k];
        }
        z[i] = s / l[i * p + i];
    }
    z
}

fn chol_inversa(l: &[f64], p: usize) -> Vec<f64> {
    let mut inv = vec![0.0; p * p];
    let mut e = vec![0.0; p];
    for j in 0..p {
        e.iter_mut().for_each(|v| *v = 0.0);
        e[j] = 1.0;
        let c = chol_resolve(l, p, &e);
        for i in 0..p {
            inv[i * p + j] = c[i];
        }
    }
    inv
}

/// X^T diag(w) X, X m x p row-major.
fn xtwx(x: &[f64], m: usize, p: usize, w: &[f64]) -> Vec<f64> {
    let mut b = vec![0.0; p * p];
    for j in 0..m {
        let xj = &x[j * p..(j + 1) * p];
        let wj = w[j];
        for r in 0..p {
            let v = wj * xj[r];
            for c in 0..=r {
                b[r * p + c] += v * xj[c];
            }
        }
    }
    for r in 0..p {
        for c in 0..r {
            b[c * p + r] = b[r * p + c];
        }
    }
    b
}

// =====================================================================
//  3) DIFFERENTIAL EXPRESSION (the DESeq2 method)
// =====================================================================

const MIN_DISP: f64 = 1e-8;
const MIN_MU: f64 = 0.5;
const KAPPA_0: f64 = 1.0;
const DISP_TOL: f64 = 1e-6;
const MAXIT: usize = 100;
const BETA_TOL: f64 = 1e-8;
/// Ridge on the coefficients, 1e-6 on the log2 scale (DESeq2's default when
/// there is no beta prior), here on the natural-log scale.
const LAMBDA_LOG2: f64 = 1e-6;

/// The model matrix and what it means.
#[derive(Clone, Debug)]
struct Desenho {
    x: Vec<f64>,
    m: usize,
    p: usize,
    nomes: Vec<String>,
    coef: usize,
    formula: String,
}

/// DESeq2's estimateSizeFactorsForMatrix: for each sample, the median over
/// the genes with no zero count of count / geometric mean of the gene.
/// `counts` is gene-major (`counts[g * n + j]`). Errors when every gene has a
/// zero somewhere; `size_factors_poscounts` handles that case.
pub fn size_factors(counts: &[f64], n: usize) -> Result<Vec<f64>, String> {
    let g = counts.len() / n.max(1);
    let mut lg = vec![f64::NAN; g];
    let mut algum = false;
    for i in 0..g {
        let l = &counts[i * n..(i + 1) * n];
        if l.iter().all(|&c| c > 0.0) {
            lg[i] = l.iter().map(|c| c.ln()).sum::<f64>() / n as f64;
            algum = true;
        }
    }
    if !algum {
        return Err("every gene has a zero count in at least one sample, so the geometric means are all zero".into());
    }
    Ok((0..n)
        .map(|j| {
            let r: Vec<f64> = (0..g).filter(|&i| lg[i].is_finite()).map(|i| counts[i * n + j].ln() - lg[i]).collect();
            median(&r).exp()
        })
        .collect())
}

/// DESeq2's type = "poscounts": the geometric mean over the positive counts
/// (divided by all n), the ratios of positive counts only, size factors
/// scaled to a geometric mean of 1. For data where every gene has a zero.
pub fn size_factors_poscounts(counts: &[f64], n: usize) -> Vec<f64> {
    let g = counts.len() / n.max(1);
    let lg: Vec<f64> = (0..g)
        .map(|i| {
            let l = &counts[i * n..(i + 1) * n];
            if l.iter().all(|&c| c == 0.0) {
                f64::NEG_INFINITY
            } else {
                l.iter().filter(|&&c| c > 0.0).map(|c| c.ln()).sum::<f64>() / n as f64
            }
        })
        .collect();
    let sf: Vec<f64> = (0..n)
        .map(|j| {
            let r: Vec<f64> =
                (0..g).filter(|&i| lg[i].is_finite() && counts[i * n + j] > 0.0).map(|i| counts[i * n + j].ln() - lg[i]).collect();
            median(&r).exp()
        })
        .collect();
    let gm = (sf.iter().map(|s| s.ln()).sum::<f64>() / n as f64).exp();
    sf.iter().map(|s| s / gm).collect()
}

// --- the dispersion likelihood (DESeq2.cpp: log_posterior, dlog_posterior)

fn log_post(a: f64, y: &[f64], mu: &[f64], x: &[f64], p: usize, prior: Option<(f64, f64)>) -> f64 {
    let alfa = a.exp();
    let s = 1.0 / alfa;
    let m = y.len();
    let mut ll = 0.0;
    let mut w = vec![0.0; m];
    for j in 0..m {
        let ma = (mu[j] * alfa).ln_1p();
        ll += ln_gamma_razao(y[j], s) - (y[j] + s) * ma;
        w[j] = 1.0 / (1.0 / mu[j] + alfa);
    }
    let b = xtwx(x, m, p, &w);
    let cr = match cholesky(&b, p) {
        Some(l) => -(0..p).map(|i| l[i * p + i].ln()).sum::<f64>(),
        None => f64::NEG_INFINITY,
    };
    let pr = match prior {
        Some((media, var)) => -0.5 * (a - media) * (a - media) / var,
        None => 0.0,
    };
    ll + pr + cr
}

fn dlog_post(a: f64, y: &[f64], mu: &[f64], x: &[f64], p: usize, prior: Option<(f64, f64)>) -> f64 {
    let alfa = a.exp();
    let s = 1.0 / alfa;
    let m = y.len();
    let mut soma = 0.0;
    let mut w = vec![0.0; m];
    let mut dw = vec![0.0; m];
    for j in 0..m {
        let ma = mu[j] * alfa;
        soma += -digamma_dif(y[j], s) + ma.ln_1p() - ma / (1.0 + ma) + y[j] * alfa / (1.0 + ma);
        w[j] = 1.0 / (1.0 / mu[j] + alfa);
        dw[j] = -w[j] * w[j];
    }
    let ll = soma / (alfa * alfa);
    let b = xtwx(x, m, p, &w);
    let db = xtwx(x, m, p, &dw);
    let cr = match cholesky(&b, p) {
        Some(l) => {
            let inv = chol_inversa(&l, p);
            let mut tr = 0.0;
            for i in 0..p {
                for k in 0..p {
                    tr += inv[i * p + k] * db[k * p + i];
                }
            }
            -0.5 * tr
        }
        None => 0.0,
    };
    let pr = match prior {
        Some((media, var)) => -(a - media) / var,
        None => 0.0,
    };
    (ll + cr) * alfa + pr
}

struct AjusteDisp {
    a: f64,
    iter: usize,
    lp_inicial: f64,
    lp_final: f64,
}

/// DESeq2.cpp fitDisp: line search (Armijo rule) on log alpha.
fn fit_disp(y: &[f64], mu: &[f64], x: &[f64], p: usize, a0: f64, prior: Option<(f64, f64)>) -> AjusteDisp {
    let min_log_alpha = (MIN_DISP / 10.0).ln();
    let epsilon = 1.0e-4;
    let mut a = a0;
    let mut lp = log_post(a, y, mu, x, p, prior);
    let mut dlp = dlog_post(a, y, mu, x, p, prior);
    let mut kappa = KAPPA_0;
    let lp_inicial = lp;
    let mut iter = 0;
    let mut aceitos = 0;
    for _ in 0..MAXIT {
        iter += 1;
        let a_prop = a + kappa * dlp;
        if a_prop < -30.0 {
            kappa = (-30.0 - a) / dlp;
        }
        if a_prop > 10.0 {
            kappa = (10.0 - a) / dlp;
        }
        let theta_kappa = -log_post(a + kappa * dlp, y, mu, x, p, prior);
        let theta_hat_kappa = -lp - kappa * epsilon * dlp * dlp;
        if theta_kappa <= theta_hat_kappa {
            aceitos += 1;
            a += kappa * dlp;
            let lpnew = log_post(a, y, mu, x, p, prior);
            let change = lpnew - lp;
            if change < DISP_TOL {
                lp = lpnew;
                break;
            }
            if a < min_log_alpha {
                break;
            }
            lp = lpnew;
            dlp = dlog_post(a, y, mu, x, p, prior);
            kappa = (kappa * 1.1).min(KAPPA_0);
            if aceitos % 5 == 0 {
                kappa /= 2.0;
            }
        } else {
            kappa /= 2.0;
        }
    }
    AjusteDisp { a, iter, lp_inicial, lp_final: lp }
}

/// DESeq2.cpp fitDispGrid: the best of 20 points of log alpha between
/// log(1e-8) and log(max(10, m)), then of 20 around it.
fn fit_disp_grade(y: &[f64], mu: &[f64], x: &[f64], p: usize, prior: Option<(f64, f64)>) -> f64 {
    let m = y.len();
    let (lo, hi) = (MIN_DISP.ln(), (m.max(10) as f64).ln());
    let grade: Vec<f64> = (0..20).map(|i| lo + (hi - lo) * i as f64 / 19.0).collect();
    let melhor = |g: &[f64]| -> f64 {
        let mut b = (g[0], f64::NEG_INFINITY);
        for &a in g {
            let v = log_post(a, y, mu, x, p, prior);
            if v > b.1 {
                b = (a, v);
            }
        }
        b.0
    };
    let a_hat = melhor(&grade);
    let delta = grade[1] - grade[0];
    let fina: Vec<f64> = (0..20).map(|i| (a_hat - delta) + 2.0 * delta * i as f64 / 19.0).collect();
    melhor(&fina)
}

/// The gene-wise dispersion of one gene (DESeq2's
/// estimateDispersionsGeneEst for one row): Cox-Reid adjusted maximum
/// likelihood with the means `mu` fixed, from the starting value `inicio`.
/// `x` is the m x p model matrix (row-major).
pub fn cox_reid_dispersion(y: &[f64], mu: &[f64], x: &[f64], p: usize, inicio: f64) -> f64 {
    let max_disp = (y.len().max(10)) as f64;
    let ini = inicio.clamp(MIN_DISP, max_disp);
    let r = fit_disp(y, mu, x, p, ini.ln(), None);
    let mut alfa = r.a.exp().min(max_disp);
    if r.lp_final < r.lp_inicial + r.lp_inicial.abs() / 1e6 {
        alfa = ini;
    }
    let conv = r.iter < MAXIT && r.iter != 1;
    if !conv && alfa > MIN_DISP * 10.0 {
        alfa = fit_disp_grade(y, mu, x, p, None).exp();
    }
    alfa.clamp(MIN_DISP, max_disp)
}

/// The maximum a posteriori dispersion of one gene (DESeq2's
/// estimateDispersionsMAP for one row, before the outlier rule): the log
/// dispersion has a normal prior with mean ln(`trend`) and variance
/// `prior_var`.
pub fn map_dispersion(y: &[f64], mu: &[f64], x: &[f64], p: usize, inicio: f64, trend: f64, prior_var: f64) -> f64 {
    let max_disp = (y.len().max(10)) as f64;
    let prior = Some((trend.ln(), prior_var));
    let r = fit_disp(y, mu, x, p, inicio.ln(), prior);
    let mut alfa = r.a.exp();
    if r.iter >= MAXIT {
        alfa = fit_disp_grade(y, mu, x, p, prior).exp();
    }
    alfa.clamp(MIN_DISP, max_disp)
}

/// ln P(Y = y) for a negative binomial of mean mu and dispersion alfa
/// (variance mu + alfa mu^2).
fn ln_nb(y: f64, mu: f64, alfa: f64) -> f64 {
    let s = 1.0 / alfa;
    let ma = mu * alfa;
    let mut v = ln_gamma_razao(y, s) + y * s.ln() - ln_gamma(y + 1.0) - s * ma.ln_1p();
    if y > 0.0 {
        v += y * (ma.ln() - ma.ln_1p());
    }
    v
}

struct AjusteBeta {
    /// Natural-log scale.
    beta: Vec<f64>,
    /// Variances of the coefficients (natural-log scale).
    var: Vec<f64>,
    /// Diagonal of the hat matrix.
    hat: Vec<f64>,
    /// sf * exp(X beta), not clamped.
    mu: Vec<f64>,
    convergiu: bool,
}

/// Least squares of log(normalized + 0.1) on X: DESeq2's starting values.
fn beta_inicial(q: &[f64], x: &[f64], m: usize, p: usize, lxtx: &[f64]) -> Vec<f64> {
    let mut b = vec![0.0; p];
    for j in 0..m {
        let v = (q[j] + 0.1).ln();
        for k in 0..p {
            b[k] += x[j * p + k] * v;
        }
    }
    chol_resolve(lxtx, p, &b)
}

fn mu_de(beta: &[f64], sf: &[f64], x: &[f64], p: usize) -> Vec<f64> {
    sf.iter()
        .enumerate()
        .map(|(j, s)| {
            let eta: f64 = (0..p).map(|k| x[j * p + k] * beta[k]).sum();
            s * eta.exp()
        })
        .collect()
}

fn finaliza_beta(beta: Vec<f64>, y: &[f64], sf: &[f64], x: &[f64], p: usize, alfa: f64, lambda: &[f64], convergiu: bool) -> AjusteBeta {
    let m = y.len();
    let mu = mu_de(&beta, sf, x, p);
    let w: Vec<f64> = mu.iter().map(|&u| {
        let u = u.max(MIN_MU);
        u / (1.0 + alfa * u)
    }).collect();
    let xw = xtwx(x, m, p, &w);
    let mut r = xw.clone();
    for k in 0..p {
        r[k * p + k] += lambda[k];
    }
    let (var, hat) = match cholesky(&r, p) {
        Some(l) => {
            let inv = chol_inversa(&l, p);
            // sigma = inv * xw * inv
            let mut t = vec![0.0; p * p];
            for i in 0..p {
                for k in 0..p {
                    let mut s = 0.0;
                    for z in 0..p {
                        s += inv[i * p + z] * xw[z * p + k];
                    }
                    t[i * p + k] = s;
                }
            }
            let var: Vec<f64> = (0..p).map(|i| (0..p).map(|z| t[i * p + z] * inv[z * p + i]).sum()).collect();
            let hat: Vec<f64> = (0..m)
                .map(|j| {
                    let xj = &x[j * p..(j + 1) * p];
                    let mut h = 0.0;
                    for a in 0..p {
                        for b in 0..p {
                            h += xj[a] * inv[a * p + b] * xj[b];
                        }
                    }
                    h * w[j]
                })
                .collect();
            (var, hat)
        }
        None => (vec![f64::NAN; p], vec![f64::NAN; m]),
    };
    AjusteBeta { beta, var, hat, mu, convergiu }
}

/// Maximizes the penalized NB log-likelihood by Fisher scoring with step
/// halving, inside |beta| <= 30 on the log2 scale: the fallback for genes
/// where the IRLS does not converge (DESeq2 uses optim's L-BFGS-B there).
fn otimiza_beta(y: &[f64], sf: &[f64], x: &[f64], p: usize, alfa: f64, lambda: &[f64], inicio: Vec<f64>) -> (Vec<f64>, bool) {
    let m = y.len();
    let lim = 30.0 * std::f64::consts::LN_2;
    let obj = |b: &[f64]| -> f64 {
        let mu = mu_de(b, sf, x, p);
        let mut v: f64 = (0..m).map(|j| ln_nb(y[j], mu[j], alfa)).sum();
        for k in 0..p {
            v -= 0.5 * lambda[k] * b[k] * b[k];
        }
        if v.is_finite() {
            v
        } else {
            -1e300
        }
    };
    let mut b: Vec<f64> = inicio.iter().map(|v| if v.is_finite() { v.clamp(-lim, lim) } else { 0.0 }).collect();
    let mut f = obj(&b);
    for _ in 0..1000 {
        let mu = mu_de(&b, sf, x, p);
        let mut g = vec![0.0; p];
        let mut w = vec![0.0; m];
        for j in 0..m {
            let r = (y[j] - mu[j]) / (1.0 + alfa * mu[j]);
            for k in 0..p {
                g[k] += x[j * p + k] * r;
            }
            w[j] = mu[j] / (1.0 + alfa * mu[j]);
        }
        for k in 0..p {
            g[k] -= lambda[k] * b[k];
        }
        let mut info = xtwx(x, m, p, &w);
        for k in 0..p {
            info[k * p + k] += lambda[k] + 1e-12;
        }
        let d = match cholesky(&info, p) {
            Some(l) => chol_resolve(&l, p, &g),
            None => g.clone(),
        };
        let mut passo = 1.0;
        let mut aceito = None;
        for _ in 0..60 {
            let nb: Vec<f64> = (0..p).map(|k| (b[k] + passo * d[k]).clamp(-lim, lim)).collect();
            let nf = obj(&nb);
            if nf >= f {
                aceito = Some((nb, nf));
                break;
            }
            passo /= 2.0;
        }
        match aceito {
            Some((nb, nf)) => {
                let mudou = (nf - f).abs();
                b = nb;
                let velho = f;
                f = nf;
                if mudou < 1e-10 * (velho.abs() + 0.1) {
                    return (b, true);
                }
            }
            None => return (b, true),
        }
    }
    (b, false)
}

/// DESeq2.cpp fitBeta: IRLS with a small ridge, for one gene; the optim
/// fallback when it does not converge, is unstable or gives a non-positive
/// variance.
fn fit_beta(y: &[f64], sf: &[f64], x: &[f64], p: usize, alfa: f64, beta0: Vec<f64>, lambda: &[f64]) -> AjusteBeta {
    let m = y.len();
    let large = 30.0;
    let mut beta = beta0.clone();
    let mut mu_hat: Vec<f64> = mu_de(&beta, sf, x, p).into_iter().map(|u| u.max(MIN_MU)).collect();
    let mut dev_old = 0.0;
    let mut convergiu = false;
    let mut w = vec![0.0; m];
    let mut ok = true;
    for t in 0..MAXIT {
        for j in 0..m {
            w[j] = mu_hat[j] / (1.0 + alfa * mu_hat[j]);
        }
        let mut a = xtwx(x, m, p, &w);
        for k in 0..p {
            a[k * p + k] += lambda[k];
        }
        let mut rhs = vec![0.0; p];
        for j in 0..m {
            let z = (mu_hat[j] / sf[j]).ln() + (y[j] - mu_hat[j]) / mu_hat[j];
            for k in 0..p {
                rhs[k] += x[j * p + k] * w[j] * z;
            }
        }
        let Some(l) = cholesky(&a, p) else {
            ok = false;
            break;
        };
        beta = chol_resolve(&l, p, &rhs);
        if beta.iter().any(|b| b.abs() > large || !b.is_finite()) {
            ok = false;
            break;
        }
        mu_hat = mu_de(&beta, sf, x, p).into_iter().map(|u| u.max(MIN_MU)).collect();
        let dev: f64 = (0..m).map(|j| -2.0 * ln_nb(y[j], mu_hat[j], alfa)).sum();
        let conv_test = (dev - dev_old).abs() / (dev.abs() + 0.1);
        if conv_test.is_nan() {
            ok = false;
            break;
        }
        if t > 0 && conv_test < BETA_TOL {
            convergiu = true;
            break;
        }
        dev_old = dev;
    }
    let ajuste = finaliza_beta(beta.clone(), y, sf, x, p, alfa, lambda, convergiu && ok);
    let estavel = ajuste.beta.iter().all(|b| b.is_finite());
    let var_pos = ajuste.var.iter().all(|v| *v > 0.0);
    if convergiu && ok && estavel && var_pos {
        return ajuste;
    }
    let inicio = if ok && estavel { beta } else { beta0 };
    let (b, conv) = otimiza_beta(y, sf, x, p, alfa, lambda, inicio);
    finaliza_beta(b, y, sf, x, p, alfa, lambda, conv)
}

// --- the dispersion trend (DESeq2's parametricDispersionFit) ------------

/// The fitted mean-dispersion relation: alpha(mean) = asymptotic +
/// extra_poisson / mean.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DispersionTrend {
    pub asymptotic: f64,
    pub extra_poisson: f64,
    /// "parametric" or "mean" (the fallback when the parametric fit fails).
    pub kind: String,
}

impl DispersionTrend {
    pub fn at(&self, mean: f64) -> f64 {
        self.asymptotic + self.extra_poisson / mean
    }
}

/// R's glm(y ~ x, family = Gamma(link = "identity"), start = start), by
/// IRLS: (coefficients, converged).
fn glm_gama_identidade(y: &[f64], x: &[f64], inicio: (f64, f64)) -> Option<((f64, f64), bool)> {
    let n = y.len();
    if n < 2 {
        return None;
    }
    let dev = |c: (f64, f64)| -> f64 {
        let mut d = 0.0;
        for i in 0..n {
            let mu = c.0 + c.1 * x[i];
            d += -2.0 * ((y[i] / mu).ln() - (y[i] - mu) / mu);
        }
        d
    };
    let valido = |c: (f64, f64)| (0..n).all(|i| c.0 + c.1 * x[i] > 0.0);
    let mut c = inicio;
    if !valido(c) {
        return None;
    }
    let mut dev_old = dev(c);
    for _ in 0..25 {
        // weighted least squares of y on (1, x) with weights 1/mu^2
        let (mut s0, mut s1, mut s2, mut t0, mut t1) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for i in 0..n {
            let mu = c.0 + c.1 * x[i];
            let w = 1.0 / (mu * mu);
            s0 += w;
            s1 += w * x[i];
            s2 += w * x[i] * x[i];
            t0 += w * y[i];
            t1 += w * x[i] * y[i];
        }
        let det = s0 * s2 - s1 * s1;
        if !(det.abs() > 0.0) {
            return None;
        }
        let mut novo = ((s2 * t0 - s1 * t1) / det, (s0 * t1 - s1 * t0) / det);
        let mut d = dev(novo);
        let mut k = 0;
        while (!valido(novo) || !d.is_finite()) && k < 25 {
            novo = ((novo.0 + c.0) / 2.0, (novo.1 + c.1) / 2.0);
            d = dev(novo);
            k += 1;
        }
        if !valido(novo) || !d.is_finite() {
            return None;
        }
        c = novo;
        if (d - dev_old).abs() / (d.abs() + 0.1) < 1e-8 {
            return Some((c, true));
        }
        dev_old = d;
    }
    Some((c, false))
}

/// DESeq2's parametricDispersionFit: iteratively, a gamma-family GLM of the
/// dispersions on 1/mean (identity link), leaving out the genes more than
/// 15x above or 1e4x below the current curve.
pub fn dispersion_trend(means: &[f64], disps: &[f64]) -> Result<DispersionTrend, String> {
    let mut coefs = (0.1, 1.0);
    let mut iter = 0;
    loop {
        let mut ys = Vec::new();
        let mut xs = Vec::new();
        for i in 0..means.len() {
            let r = disps[i] / (coefs.0 + coefs.1 / means[i]);
            if r > 1e-4 && r < 15.0 {
                ys.push(disps[i]);
                xs.push(1.0 / means[i]);
            }
        }
        let Some((novo, conv)) = glm_gama_identidade(&ys, &xs, coefs) else {
            return Err("the gamma GLM of the dispersion trend could not be fitted".into());
        };
        let velho = coefs;
        coefs = novo;
        if !(coefs.0 > 0.0 && coefs.1 > 0.0) {
            return Err("parametric dispersion fit failed (a coefficient is not positive)".into());
        }
        let mudanca = (coefs.0 / velho.0).ln().powi(2) + (coefs.1 / velho.1).ln().powi(2);
        if mudanca < 1e-6 && conv {
            break;
        }
        iter += 1;
        if iter > 10 {
            return Err("the parametric dispersion fit did not converge".into());
        }
    }
    Ok(DispersionTrend { asymptotic: coefs.0, extra_poisson: coefs.1, kind: "parametric".into() })
}

// --- the prior variance for df <= 3 -------------------------------------

/// P(Y <= b) for Y = ln(X / df) + N(0, v), X ~ chi-squared(df): the integral
/// over the normal of the chi-squared CDF (Simpson, 801 points).
fn cdf_log_chi2_normal(b: f64, df: f64, v: f64) -> f64 {
    let cdf_w = |w: f64| -> f64 {
        let t = df * w.min(700.0).exp() / 2.0;
        gama_inc(df / 2.0, t)
    };
    if v <= 0.0 {
        return cdf_w(b);
    }
    let sd = v.sqrt();
    let n = 800;
    let (lo, hi) = (-8.0, 8.0);
    let h = (hi - lo) / n as f64;
    let mut s = 0.0;
    for i in 0..=n {
        let z = lo + h * i as f64;
        let phi = (-0.5 * z * z).exp() * 0.398942280401432677939946059934;
        let f = phi * cdf_w(b - sd * z);
        let c = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        s += c * f;
    }
    s * h / 3.0
}

/// DESeq2's estimateDispersionsPriorVar when the residual degrees of freedom
/// are 1 to 3: the variance v of a normal that, added to the sampling
/// distribution of a log dispersion (ln(chi2_df / df)), best matches the
/// histogram of the observed log residuals (Kullback-Leibler over 40 bins of
/// width 0.5 between -10 and 10). DESeq2 draws 10,000 random values per
/// candidate v and smooths the curve with loess; here the bin probabilities
/// are computed exactly, and the floor added to the densities is the one a
/// 10,000-draw histogram has. The result is at least 0.25.
fn variancia_prior_kl(residuos: &[f64], df: f64) -> f64 {
    let largura = 0.5;
    let brks: Vec<f64> = (-20..=20).map(|i| i as f64 / 2.0).collect();
    let obs: Vec<f64> = residuos.iter().copied().filter(|r| *r > -10.0 && *r < 10.0).collect();
    let mut cont = vec![0usize; 40];
    for r in &obs {
        let k = (((r + 10.0) / largura).ceil() as isize - 1).clamp(0, 39) as usize;
        cont[k] += 1;
    }
    let n_obs = obs.len().max(1) as f64;
    let dens_obs: Vec<f64> = cont.iter().map(|&c| c as f64 / (n_obs * largura)).collect();
    let min_obs = dens_obs.iter().copied().filter(|d| *d > 0.0).fold(f64::INFINITY, f64::min);
    let piso = min_obs.min(1.0 / (1e4 * largura));
    let kl = |v: f64| -> f64 {
        let c: Vec<f64> = brks.iter().map(|&b| cdf_log_chi2_normal(b, df, v)).collect();
        let total = c[40] - c[0];
        let mut s = 0.0;
        for k in 0..40 {
            if dens_obs[k] > 0.0 {
                let dm = ((c[k + 1] - c[k]) / total / largura).max(0.0);
                s += dens_obs[k] * ((dens_obs[k] + piso).ln() - (dm + piso).ln());
            }
        }
        s
    };
    let grade: Vec<f64> = (0..200).map(|i| 8.0 * i as f64 / 199.0).collect();
    let valores: Vec<f64> = grade.par_iter().map(|&v| kl(v)).collect();
    let mut im = 0;
    for i in 1..200 {
        if valores[i] < valores[im] {
            im = i;
        }
    }
    // golden-section search between the neighbours of the best grid point
    let (mut a, mut b) = (grade[im.saturating_sub(1)], grade[(im + 1).min(199)]);
    let g = 0.5 * (5f64.sqrt() - 1.0);
    let mut c = b - g * (b - a);
    let mut d = a + g * (b - a);
    let (mut fc, mut fd) = (kl(c), kl(d));
    for _ in 0..40 {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - g * (b - a);
            fc = kl(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + g * (b - a);
            fd = kl(d);
        }
    }
    let mut melhor = 0.5 * (a + b);
    if valores[im] < kl(melhor) {
        melhor = grade[im];
    }
    melhor.max(0.25)
}

// --- the model -----------------------------------------------------------

fn ordena_niveis(v: &[String]) -> Vec<String> {
    let mut n: Vec<String> = v.to_vec();
    n.sort();
    n.dedup();
    n
}

/// Builds `~ factor1 + ... + condition`: an intercept, one indicator per
/// non-reference level of each blocking factor (reference: the first level
/// in sorted order, as R does) and of the condition.
fn monta_desenho(
    condicao: &str,
    niveis: &[String],
    referencia: &str,
    tratado: &str,
    fatores: &[(String, Vec<String>)],
    avisos: &mut Vec<String>,
) -> Result<Desenho, String> {
    let m = niveis.len();
    let mut colunas: Vec<(String, Vec<f64>)> = vec![("Intercept".into(), vec![1.0; m])];
    let mut termos: Vec<String> = Vec::new();
    for (nome, v) in fatores {
        let lv = ordena_niveis(v);
        if lv.len() < 2 {
            avisos.push(format!("the factor {nome} has a single level ({}) among these samples and was left out of the design", lv.join("")));
            continue;
        }
        termos.push(nome.clone());
        for l in &lv[1..] {
            colunas.push((format!("{nome}_{l}_vs_{}", lv[0]), v.iter().map(|x| if x == l { 1.0 } else { 0.0 }).collect()));
        }
    }
    termos.push(condicao.to_string());
    let lv = ordena_niveis(niveis);
    let mut coef = 0;
    for l in lv.iter().filter(|l| *l != referencia) {
        if l == tratado {
            coef = colunas.len();
        }
        colunas.push((format!("{condicao}_{l}_vs_{referencia}"), niveis.iter().map(|x| if x == l { 1.0 } else { 0.0 }).collect()));
    }
    let p = colunas.len();
    let mut x = vec![0.0; m * p];
    for j in 0..m {
        for (k, c) in colunas.iter().enumerate() {
            x[j * p + k] = c.1[j];
        }
    }
    let formula = format!("~ {}", termos.join(" + "));
    let xtx = xtwx(&x, m, p, &vec![1.0; m]);
    let ok = cholesky(&xtx, p).map_or(false, |l| {
        let d: Vec<f64> = (0..p).map(|i| l[i * p + i]).collect();
        let maior = d.iter().fold(0.0f64, |a, b| a.max(*b));
        d.iter().all(|v| *v > 1e-8 * maior)
    });
    if !ok {
        return Err(format!(
            "The design {formula} cannot be fitted: its model matrix is not full rank -- some factor is a combination of the others (for example a blocking factor whose levels follow the condition's, so that every level of it has a single level of {condicao}). Remove that factor or check the sample table."
        ));
    }
    if m <= p {
        return Err(format!(
            "{m} samples and {p} coefficients in the design {formula}: with no residual degrees of freedom the dispersions cannot be estimated. Differential expression needs replicates (more samples than coefficients)."
        ));
    }
    Ok(Desenho { x, m, p, nomes: colunas.into_iter().map(|c| c.0).collect(), coef, formula })
}

#[derive(Clone, Debug)]
pub struct DeOptions {
    /// The column of the sample table that is compared.
    pub condition: String,
    /// Its reference (control) level: log2 fold changes are treated / reference.
    pub reference: String,
    /// The level compared with the reference. `None`: the only other level.
    pub treated: Option<String>,
    /// Blocking factors (columns of the sample table), e.g. the cell line or
    /// the batch of a paired design.
    pub factors: Vec<String>,
    /// Use only the samples at the reference and treated levels (otherwise
    /// all samples enter the fit and the other levels get their own
    /// coefficients, as in DESeq2 with a multi-level factor).
    pub only_compared_levels: bool,
    /// Significance level: the target of the independent filtering and of
    /// the summary. Default 0.1, as DESeq2.
    pub alpha: f64,
    /// |log2 fold change| above which a significant gene counts as up/down
    /// in the summary (not a test threshold). Default 0.
    pub lfc_threshold: f64,
    /// DESeq2's independent filtering by mean count. Default true.
    pub independent_filtering: bool,
    /// Flag genes with an outlier count (Cook's distance) and leave their
    /// p-value out. Default true.
    pub cooks_filter: bool,
    /// The variance of the prior on the log dispersions; `None` estimates it.
    pub prior_variance: Option<f64>,
    /// A table of gene names: NCBI's gene_info (GeneID -> Symbol) or two
    /// columns (gene id, name). Adds a `name` column.
    pub gene_names: Option<String>,
}

impl Default for DeOptions {
    fn default() -> Self {
        DeOptions {
            condition: "condition".into(),
            reference: String::new(),
            treated: None,
            factors: Vec::new(),
            only_compared_levels: false,
            alpha: 0.1,
            lfc_threshold: 0.0,
            independent_filtering: true,
            cooks_filter: true,
            prior_variance: None,
            gene_names: None,
        }
    }
}

/// One gene of the result. NaN where DESeq2 writes NA.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DeGene {
    pub gene: String,
    /// From `gene_names`; the gene id itself when the table does not have it.
    pub name: String,
    pub base_mean: f64,
    pub log2_fold_change: f64,
    pub lfc_se: f64,
    pub stat: f64,
    pub pvalue: f64,
    pub padj: f64,
    /// Gene-wise, trend and final dispersions (for a dispersion plot).
    pub disp_gene: f64,
    pub disp_trend: f64,
    pub dispersion: f64,
    pub max_cooks: f64,
    /// Every count is zero: nothing is tested.
    pub all_zero: bool,
    /// An outlier count (Cook's distance): p-value left out.
    pub outlier: bool,
    /// Removed by the independent filtering (low mean count).
    pub filtered: bool,
    /// The GLM did not converge (the estimate is still reported).
    pub not_converged: bool,
}

/// The result of `differential_expression`.
///
/// Differences from DESeq2, all documented in the text summary too:
///  - with 1 to 3 residual degrees of freedom, the prior variance of the log
///    dispersions is matched to the observed residuals with exact densities
///    instead of DESeq2's 10,000 random draws smoothed by loess
///    (`prior_method` = "distribution_matching");
///  - the dispersion trend is parametric; if that fit fails, the mean of the
///    gene-wise estimates is used (DESeq2 would switch to a local fit);
///  - genes with a Cook's outlier get their p-value removed, as DESeq2 does
///    for cells with fewer than 7 replicates; DESeq2's replacement of
///    outliers and refit in cells with 7 or more replicates is not done;
///  - no log fold change shrinkage (lfcShrink); the fold changes are the MLE,
///    as in DESeq2's `results()`.
#[derive(Clone, Debug, Default)]
pub struct DeResult {
    pub counts: String,
    pub sample_table: String,
    pub samples: Vec<String>,
    pub size_factors: Vec<f64>,
    /// "ratio" (median of ratios) or "poscounts" (when every gene has a zero).
    pub size_factor_method: String,
    pub design: String,
    pub coefficients: Vec<String>,
    pub condition: String,
    pub reference: String,
    pub treated: String,
    /// Genes in the input order.
    pub genes: Vec<DeGene>,
    pub trend: DispersionTrend,
    pub log_dispersion_variance: f64,
    pub prior_variance: f64,
    /// "mad" (residual df > 3), "distribution_matching" (1 to 3) or "given".
    pub prior_method: String,
    pub residual_df: usize,
    pub alpha: f64,
    pub lfc_threshold: f64,
    pub tested: usize,
    pub all_zero: usize,
    pub up: usize,
    pub down: usize,
    pub outliers: usize,
    pub low_counts: usize,
    /// Genes with a mean normalized count below this were filtered (0: none).
    pub filter_threshold: f64,
    /// The quantile of the mean counts the threshold sits at.
    pub filter_quantile: f64,
    pub not_converged: usize,
    pub names_missing: usize,
    /// Columns of the count matrix not in the sample table (not used).
    pub ignored_columns: Vec<String>,
    /// Rows of the matrix that are not genes (htseq-count's `__no_feature`...).
    pub skipped_rows: usize,
    pub output: Option<String>,
    pub warnings: Vec<String>,
}

impl DeResult {
    /// The genes in the order of the output table: by adjusted p-value,
    /// missing values last, ties in input order (R's `order(res$padj)`).
    pub fn order(&self) -> Vec<u32> {
        let mut idx: Vec<u32> = (0..self.genes.len() as u32).collect();
        let g = &self.genes;
        multi_merge_sort_stable_by(&mut idx, |a, b| {
            let (x, y) = (g[*a as usize].padj, g[*b as usize].padj);
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
        let g = &self.genes;
        let col = |f: &dyn Fn(&DeGene) -> f64| -> String { jnums(&g.iter().map(f).collect::<Vec<_>>()) };
        format!(
            "{{\"counts\":{},\"sample_table\":{},\"samples\":{},\"size_factors\":{},\"size_factor_method\":{},\"design\":{},\"coefficients\":{},\"condition\":{},\"reference\":{},\"treated\":{},\"trend\":{{\"asymptotic\":{},\"extra_poisson\":{},\"kind\":{}}},\"log_dispersion_variance\":{},\"prior_variance\":{},\"prior_method\":{},\"residual_df\":{},\"alpha\":{},\"lfc_threshold\":{},\"genes_total\":{},\"tested\":{},\"all_zero\":{},\"up\":{},\"down\":{},\"outliers\":{},\"low_counts\":{},\"filter_threshold\":{},\"filter_quantile\":{},\"not_converged\":{},\"names_missing\":{},\"ignored_columns\":{},\"skipped_rows\":{},\"output\":{},\"warnings\":{},\"plot\":{{\"gene\":{},\"name\":{},\"base_mean\":{},\"log2_fold_change\":{},\"lfc_se\":{},\"stat\":{},\"pvalue\":{},\"padj\":{},\"disp_gene\":{},\"disp_trend\":{},\"dispersion\":{},\"outlier\":{},\"filtered\":{}}},\"text\":{}}}",
            json_texto(&self.counts),
            json_texto(&self.sample_table),
            jtextos(&self.samples),
            jnums(&self.size_factors),
            json_texto(&self.size_factor_method),
            json_texto(&self.design),
            jtextos(&self.coefficients),
            json_texto(&self.condition),
            json_texto(&self.reference),
            json_texto(&self.treated),
            jnum(self.trend.asymptotic),
            jnum(self.trend.extra_poisson),
            json_texto(&self.trend.kind),
            jnum(self.log_dispersion_variance),
            jnum(self.prior_variance),
            json_texto(&self.prior_method),
            self.residual_df,
            jnum(self.alpha),
            jnum(self.lfc_threshold),
            g.len(),
            self.tested,
            self.all_zero,
            self.up,
            self.down,
            self.outliers,
            self.low_counts,
            jnum(self.filter_threshold),
            jnum(self.filter_quantile),
            self.not_converged,
            self.names_missing,
            jtextos(&self.ignored_columns),
            self.skipped_rows,
            jopt(self.output.as_deref()),
            jtextos(&self.warnings),
            jtextos(&g.iter().map(|x| x.gene.clone()).collect::<Vec<_>>()),
            jtextos(&g.iter().map(|x| x.name.clone()).collect::<Vec<_>>()),
            col(&|x| x.base_mean),
            col(&|x| x.log2_fold_change),
            col(&|x| x.lfc_se),
            col(&|x| x.stat),
            col(&|x| x.pvalue),
            col(&|x| x.padj),
            col(&|x| x.disp_gene),
            col(&|x| x.disp_trend),
            col(&|x| x.dispersion),
            format!("[{}]", g.iter().map(|x| x.outlier.to_string()).collect::<Vec<_>>().join(",")),
            format!("[{}]", g.iter().map(|x| x.filtered.to_string()).collect::<Vec<_>>().join(",")),
            json_texto(&self.resumo())
        )
    }

    pub fn resumo(&self) -> String {
        let nz = self.genes.len() - self.all_zero;
        let mut o = format!(
            "Differential expression (DESeq2 method): {} {} vs {}, design {}\n",
            self.condition, self.treated, self.reference, self.design
        );
        o.push_str(&format!(
            "  {} samples; size factors ({}) {}\n",
            self.samples.len(),
            self.size_factor_method,
            self.samples.iter().zip(&self.size_factors).map(|(s, f)| format!("{s} {f:.3}")).collect::<Vec<_>>().join(", ")
        ));
        o.push_str(&format!(
            "  dispersion trend ({}): {:.4} + {:.4} / mean; prior variance of log dispersions {:.3} ({}, residual df {})\n",
            self.trend.kind, self.trend.asymptotic, self.trend.extra_poisson, self.prior_variance, self.prior_method, self.residual_df
        ));
        o.push_str(&format!("  out of {} genes with nonzero total read count, adjusted p-value < {}:\n", milhar(nz), self.alpha));
        let lim = if self.lfc_threshold > 0.0 { format!(" > {}", self.lfc_threshold) } else { " > 0".into() };
        let liminf = if self.lfc_threshold > 0.0 { format!(" < -{}", self.lfc_threshold) } else { " < 0".into() };
        o.push_str(&format!("    LFC{lim} (up):   {} ({})\n", milhar(self.up), percent(self.up, nz)));
        o.push_str(&format!("    LFC{liminf} (down): {} ({})\n", milhar(self.down), percent(self.down, nz)));
        o.push_str(&format!("    outliers (Cook's distance): {} ({})\n", milhar(self.outliers), percent(self.outliers, nz)));
        o.push_str(&format!(
            "    low counts (independent filtering): {} ({}); mean count < {}\n",
            milhar(self.low_counts),
            percent(self.low_counts, nz),
            num((self.filter_threshold * 100.0).round() / 100.0)
        ));
        if self.not_converged > 0 {
            o.push_str(&format!("  {} gene(s) whose GLM did not converge (estimates still reported)\n", self.not_converged));
        }
        if let Some(p) = &self.output {
            o.push_str(&format!("  results written to {p}, sorted by adjusted p-value\n"));
        }
        for w in &self.warnings {
            o.push_str(&format!("  note: {w}\n"));
        }
        o
    }
}

/// DESeq2 on a count matrix held in memory.
///
/// `counts` is gene-major (`counts[g * samples.len() + j]`), whole numbers.
/// `condition` gives each sample's level; `factors` are the blocking factors
/// (name, the level of each sample).
pub fn deseq(
    counts: &[f64],
    genes: &[String],
    samples: &[String],
    condition: &[String],
    factors: &[(String, Vec<String>)],
    opt: &DeOptions,
) -> Result<DeResult, String> {
    let n = samples.len();
    let n_genes = genes.len();
    if counts.len() != n * n_genes {
        return Err(format!("the count matrix has {} values for {} genes x {} samples", counts.len(), n_genes, n));
    }
    if condition.len() != n || factors.iter().any(|f| f.1.len() != n) {
        return Err("the condition and every factor need one level per sample".into());
    }
    if n_genes == 0 {
        return Err("the count matrix has no gene".into());
    }
    let lv = ordena_niveis(condition);
    let referencia = opt.reference.trim();
    if referencia.is_empty() {
        return Err(format!("Give the reference level of {} (one of: {})", opt.condition, lv.join(", ")));
    }
    if !lv.iter().any(|l| l == referencia) {
        return Err(format!(
            "The reference level \"{referencia}\" is not a level of {} in these samples; its levels are: {}",
            opt.condition,
            lv.join(", ")
        ));
    }
    let tratado = match opt.treated.as_deref().map(str::trim).filter(|t| !t.is_empty()) {
        Some(t) => {
            if !lv.iter().any(|l| l == t) {
                return Err(format!("The level \"{t}\" is not a level of {}; its levels are: {}", opt.condition, lv.join(", ")));
            }
            if t == referencia {
                return Err(format!("The treated level and the reference are the same (\"{t}\"): choose two different levels"));
            }
            t.to_string()
        }
        None => {
            let outros: Vec<&String> = lv.iter().filter(|l| *l != referencia).collect();
            match outros.len() {
                0 => {
                    return Err(format!(
                        "{} has a single level (\"{referencia}\") in these samples: there is nothing to compare it with",
                        opt.condition
                    ))
                }
                1 => outros[0].clone(),
                _ => {
                    return Err(format!(
                        "{} has {} levels ({}): say which one is compared with \"{referencia}\" (treated level)",
                        opt.condition,
                        lv.len(),
                        lv.join(", ")
                    ))
                }
            }
        }
    };
    let mut avisos: Vec<String> = Vec::new();
    let des = monta_desenho(&opt.condition, condition, referencia, &tratado, factors, &mut avisos)?;
    let (m, p) = (des.m, des.p);
    for (i, v) in counts.iter().enumerate() {
        if !(v.is_finite() && *v >= 0.0 && v.fract() == 0.0) {
            return Err(format!(
                "gene {}, sample {}: count {v} is not a non-negative whole number",
                genes[i / n],
                samples[i % n]
            ));
        }
    }
    if counts.iter().all(|c| *c == 0.0) {
        return Err("every count is zero".into());
    }
    if let Some(j) = (0..n).find(|&j| (0..n_genes).all(|g| counts[g * n + j] == 0.0)) {
        return Err(format!("sample {} has no read in any gene: remove it from the sample table", samples[j]));
    }

    // size factors
    let (sf, metodo_sf) = match size_factors(counts, n) {
        Ok(s) => (s, "ratio"),
        Err(_) => {
            avisos.push("every gene has a zero in some sample, so the size factors were computed from the positive counts only (DESeq2's \"poscounts\")".into());
            (size_factors_poscounts(counts, n), "poscounts")
        }
    };
    if sf.iter().any(|s| !(s.is_finite() && *s > 0.0)) {
        return Err("a sample has no counts in the genes used for the size factors: remove empty samples".into());
    }
    let max_disp = (m.max(10)) as f64;
    let lambda: Vec<f64> = vec![LAMBDA_LOG2 / (std::f64::consts::LN_2 * std::f64::consts::LN_2); p];

    // per gene: normalized counts, mean, variance
    let norm = |g: usize| -> Vec<f64> { (0..n).map(|j| counts[g * n + j] / sf[j]).collect() };
    let mut base_mean = vec![0.0; n_genes];
    let mut base_var = vec![0.0; n_genes];
    let mut zero = vec![false; n_genes];
    for g in 0..n_genes {
        let q = norm(g);
        let mu = q.iter().sum::<f64>() / n as f64;
        base_mean[g] = mu;
        base_var[g] = q.iter().map(|v| (v - mu) * (v - mu)).sum::<f64>() / (n as f64 - 1.0);
        zero[g] = counts[g * n..(g + 1) * n].iter().all(|c| *c == 0.0);
    }
    let nz: Vec<usize> = (0..n_genes).filter(|&g| !zero[g]).collect();

    // hat matrix of the linear model and (X^T X)^-1
    let xtx = xtwx(&des.x, m, p, &vec![1.0; m]);
    let lxtx = cholesky(&xtx, p).ok_or("the model matrix is not full rank")?;
    let ixtx = chol_inversa(&lxtx, p);
    let mut hat_lin = vec![0.0; m * m];
    for i in 0..m {
        for k in 0..m {
            let mut s = 0.0;
            for a in 0..p {
                for b in 0..p {
                    s += des.x[i * p + a] * ixtx[a * p + b] * des.x[k * p + b];
                }
            }
            hat_lin[i * m + k] = s;
        }
    }
    let grupos = {
        let mut linhas: Vec<Vec<u64>> = (0..m).map(|j| des.x[j * p..(j + 1) * p].iter().map(|v| v.to_bits()).collect()).collect();
        linhas.sort();
        linhas.dedup();
        linhas.len()
    };
    let linear_mu = grupos == p;
    let xim = sf.iter().map(|s| 1.0 / s).sum::<f64>() / n as f64;

    // 1. gene-wise dispersions (estimateDispersionsGeneEst)
    struct Gene1 {
        disp: f64,
        mu: Vec<f64>,
    }
    let x = &des.x;
    let est1: Vec<Gene1> = nz
        .par_iter()
        .map(|&g| {
            let y = &counts[g * n..(g + 1) * n];
            let q = norm(g);
            let mut mu_lin = vec![0.0; m];
            for i in 0..m {
                mu_lin[i] = (0..m).map(|k| hat_lin[i * m + k] * q[k]).sum();
            }
            let mut rough = 0.0;
            for j in 0..m {
                let u = mu_lin[j].max(1.0);
                rough += ((q[j] - u) * (q[j] - u) - u) / (u * u);
            }
            let rough = (rough / (m - p) as f64).max(0.0);
            let bm = base_mean[g];
            let momentos = (base_var[g] - xim * bm) / (bm * bm);
            let ini = rough.min(momentos).max(MIN_DISP).min(max_disp);
            let ini = if ini.is_nan() { MIN_DISP } else { ini };
            let mu: Vec<f64> = if linear_mu {
                (0..m).map(|j| (mu_lin[j] * sf[j]).max(MIN_MU)).collect()
            } else {
                let b0 = beta_inicial(&q, x, m, p, &lxtx);
                let a = fit_beta(y, &sf, x, p, ini, b0, &lambda);
                a.mu.iter().map(|u| u.max(MIN_MU)).collect()
            };
            let disp = cox_reid_dispersion(y, &mu, x, p, ini);
            Gene1 { disp, mu }
        })
        .collect();

    // 2. the trend (estimateDispersionsFit)
    let usa: Vec<usize> = (0..nz.len()).filter(|&i| est1[i].disp > 100.0 * MIN_DISP).collect();
    if usa.is_empty() {
        return Err("all gene-wise dispersion estimates are within 2 orders of magnitude of the minimum value (1e-8): the counts look like Poisson (technical replicates?) and the DESeq2 method does not apply".into());
    }
    let medias: Vec<f64> = usa.iter().map(|&i| base_mean[nz[i]]).collect();
    let disps: Vec<f64> = usa.iter().map(|&i| est1[i].disp).collect();
    let trend = match dispersion_trend(&medias, &disps) {
        Ok(t) => t,
        Err(e) => {
            let v: Vec<f64> = est1.iter().map(|e| e.disp).filter(|d| *d > 10.0 * MIN_DISP).collect();
            avisos.push(format!("{e}; the dispersion trend is the mean of the gene-wise estimates instead (DESeq2 would use a local fit)"));
            DispersionTrend { asymptotic: media_aparada(&v, 0.001), extra_poisson: 0.0, kind: "mean".into() }
        }
    };
    let disp_fit: Vec<f64> = nz.iter().map(|&g| trend.at(base_mean[g])).collect();
    let acima: Vec<usize> = (0..nz.len()).filter(|&i| est1[i].disp >= 100.0 * MIN_DISP).collect();
    let res_log: Vec<f64> = acima.iter().map(|&i| est1[i].disp.ln() - disp_fit[i].ln()).collect();
    let var_log = mad(&res_log).powi(2);

    // 3. the prior variance (estimateDispersionsPriorVar)
    let df = m - p;
    let (prior_var, metodo_prior) = match opt.prior_variance {
        Some(v) => {
            if !(v > 0.0 && v.is_finite()) {
                return Err(format!("prior_variance must be a positive number, not {v}"));
            }
            (v, "given")
        }
        None if df <= 3 => (variancia_prior_kl(&res_log, df as f64), "distribution_matching"),
        None => ((var_log - trigamma(df as f64 / 2.0)).max(0.25), "mad"),
    };

    // 4. MAP dispersions (estimateDispersionsMAP)
    let limite_outlier = 2.0 * var_log.sqrt();
    let disp_final: Vec<f64> = (0..nz.len())
        .into_par_iter()
        .map(|i| {
            let g = nz[i];
            let y = &counts[g * n..(g + 1) * n];
            let dg = est1[i].disp;
            let df_ = disp_fit[i];
            let ini = if dg > 0.1 * df_ { dg } else { df_ };
            let d = map_dispersion(y, &est1[i].mu, x, p, ini, df_, prior_var);
            if dg.ln() > df_.ln() + limite_outlier {
                dg
            } else {
                d
            }
        })
        .collect();

    // 5. the GLM and the Wald test (nbinomWaldTest)
    let ajustes: Vec<AjusteBeta> = (0..nz.len())
        .into_par_iter()
        .map(|i| {
            let g = nz[i];
            let y = &counts[g * n..(g + 1) * n];
            let q = norm(g);
            let b0 = beta_inicial(&q, x, m, p, &lxtx);
            fit_beta(y, &sf, x, p, disp_final[i], b0, &lambda)
        })
        .collect();

    // 6. Cook's distances (calculateCooksDistance, recordMaxCooks)
    let linhas: Vec<Vec<u64>> = (0..m).map(|j| x[j * p..(j + 1) * p].iter().map(|v| v.to_bits()).collect()).collect();
    let tres_ou_mais: Vec<bool> = (0..m).map(|j| linhas.iter().filter(|l| **l == linhas[j]).count() >= 3).collect();
    let faz_cooks = opt.cooks_filter && m > p && tres_ou_mais.iter().any(|b| *b);
    let corte_cooks = f_quantile(0.99, p as f64, (m - p) as f64);
    let celulas: Vec<usize> = (0..m).map(|j| linhas.iter().position(|l| *l == linhas[j]).unwrap()).collect();
    let max_cooks: Vec<(f64, bool)> = if faz_cooks {
        (0..nz.len())
            .into_par_iter()
            .map(|i| {
                let g = nz[i];
                let q = norm(g);
                let y = &counts[g * n..(g + 1) * n];
                // robustMethodOfMomentsDisp
                let v = {
                    let mut maior = f64::NEG_INFINITY;
                    let mut cel: Vec<usize> = (0..m).filter(|&j| tres_ou_mais[j]).map(|j| celulas[j]).collect();
                    cel.sort_unstable();
                    cel.dedup();
                    for c in cel {
                        let membros: Vec<f64> = (0..m).filter(|&j| tres_ou_mais[j] && celulas[j] == c).map(|j| q[j]).collect();
                        let k = membros.len();
                        let (trim, escala) = if k <= 3 {
                            (1.0 / 3.0, 2.04)
                        } else if k <= 23 {
                            (1.0 / 4.0, 1.86)
                        } else {
                            (1.0 / 8.0, 1.51)
                        };
                        let media = media_aparada(&membros, trim);
                        let sq: Vec<f64> = membros.iter().map(|v| (v - media) * (v - media)).collect();
                        maior = maior.max(escala * media_aparada(&sq, trim));
                    }
                    maior
                };
                let media = q.iter().sum::<f64>() / m as f64;
                let alfa = ((v - media) / (media * media)).max(0.04);
                let a = &ajustes[i];
                let cooks: Vec<f64> = (0..m)
                    .map(|j| {
                        let mu = a.mu[j];
                        let var = mu + alfa * mu * mu;
                        let pr = (y[j] - mu) * (y[j] - mu) / var;
                        let h = a.hat[j];
                        pr / p as f64 * h / ((1.0 - h) * (1.0 - h))
                    })
                    .collect();
                let mx = (0..m).filter(|&j| tres_ou_mais[j]).map(|j| cooks[j]).fold(f64::NEG_INFINITY, f64::max);
                let mut outlier = mx > corte_cooks;
                if outlier {
                    let jmax = (0..m).fold(0, |b, j| if cooks[j] > cooks[b] { j } else { b });
                    let c_out = y[jmax];
                    if y.iter().filter(|c| **c > c_out).count() >= 3 {
                        outlier = false;
                    }
                }
                (mx, outlier)
            })
            .collect()
    } else {
        vec![(f64::NAN, false); nz.len()]
    };

    // the result per gene
    let log2e = std::f64::consts::LOG2_E;
    let mut res_genes: Vec<DeGene> = (0..n_genes)
        .map(|g| DeGene {
            gene: genes[g].clone(),
            name: genes[g].clone(),
            base_mean: base_mean[g],
            log2_fold_change: f64::NAN,
            lfc_se: f64::NAN,
            stat: f64::NAN,
            pvalue: f64::NAN,
            padj: f64::NAN,
            disp_gene: f64::NAN,
            disp_trend: f64::NAN,
            dispersion: f64::NAN,
            max_cooks: f64::NAN,
            all_zero: zero[g],
            ..Default::default()
        })
        .collect();
    let mut nao_conv = 0;
    for (i, &g) in nz.iter().enumerate() {
        let a = &ajustes[i];
        let r = &mut res_genes[g];
        r.disp_gene = est1[i].disp;
        r.disp_trend = disp_fit[i];
        r.dispersion = disp_final[i];
        r.log2_fold_change = a.beta[des.coef] * log2e;
        r.lfc_se = a.var[des.coef].max(0.0).sqrt() * log2e;
        r.stat = r.log2_fold_change / r.lfc_se;
        r.pvalue = two_sided_p(r.stat);
        r.max_cooks = max_cooks[i].0;
        r.not_converged = !a.convergiu;
        if !a.convergiu {
            nao_conv += 1;
        }
        if max_cooks[i].1 {
            r.outlier = true;
            r.pvalue = f64::NAN;
        }
    }

    // 7. independent filtering and BH (pValueAdjustment)
    let pv: Vec<f64> = res_genes.iter().map(|r| r.pvalue).collect();
    let mut limiar = 0.0;
    let mut quantil = 0.0;
    let padj: Vec<f64> = if opt.independent_filtering {
        let filtro: Vec<f64> = base_mean.clone();
        let mut ord = filtro.clone();
        ord.sort_by(|a, b| a.total_cmp(b));
        let lower = filtro.iter().filter(|v| **v == 0.0).count() as f64 / n_genes as f64;
        let upper = if lower < 0.95 { 0.95 } else { 1.0 };
        let theta: Vec<f64> = (0..50).map(|i| lower + (upper - lower) * i as f64 / 49.0).collect();
        let cortes: Vec<f64> = theta.iter().map(|t| quantil7(&ord, *t)).collect();
        let ajustados: Vec<Vec<f64>> = cortes
            .par_iter()
            .map(|c| {
                let usa: Vec<usize> = (0..n_genes).filter(|&g| filtro[g] >= *c).collect();
                let sub: Vec<f64> = usa.iter().map(|&g| pv[g]).collect();
                let adj = benjamini_hochberg(&sub);
                let mut v = vec![f64::NAN; n_genes];
                for (k, &g) in usa.iter().enumerate() {
                    v[g] = adj[k];
                }
                v
            })
            .collect();
        let num_rej: Vec<f64> = ajustados.iter().map(|v| v.iter().filter(|x| **x < opt.alpha).count() as f64).collect();
        let lo = lowess(&theta, &num_rej, 1.0 / 5.0);
        let maxr = num_rej.iter().fold(0.0f64, |a, b| a.max(*b));
        let j = if maxr <= 10.0 {
            0
        } else {
            let resid: Vec<f64> = if num_rej.iter().all(|v| *v == 0.0) {
                vec![0.0]
            } else {
                (0..50).filter(|&k| num_rej[k] > 0.0).map(|k| num_rej[k] - lo[k]).collect()
            };
            let rms = (resid.iter().map(|r| r * r).sum::<f64>() / resid.len() as f64).sqrt();
            let thresh = lo.iter().fold(f64::NEG_INFINITY, |a, b| a.max(*b)) - rms;
            (0..50).find(|&k| num_rej[k] > thresh).unwrap_or(0)
        };
        limiar = cortes[j];
        quantil = theta[j];
        ajustados[j].clone()
    } else {
        benjamini_hochberg(&pv)
    };
    let mut r = DeResult {
        samples: samples.to_vec(),
        size_factors: sf,
        size_factor_method: metodo_sf.into(),
        design: des.formula.clone(),
        coefficients: des.nomes.clone(),
        condition: opt.condition.clone(),
        reference: referencia.to_string(),
        treated: tratado,
        trend,
        log_dispersion_variance: var_log,
        prior_variance: prior_var,
        prior_method: metodo_prior.into(),
        residual_df: df,
        alpha: opt.alpha,
        lfc_threshold: opt.lfc_threshold,
        filter_threshold: limiar,
        filter_quantile: quantil,
        not_converged: nao_conv,
        warnings: avisos,
        ..Default::default()
    };
    for (g, rg) in res_genes.iter_mut().enumerate() {
        rg.padj = padj[g];
        if rg.all_zero {
            r.all_zero += 1;
            continue;
        }
        if rg.outlier {
            r.outliers += 1;
        } else if !rg.pvalue.is_nan() && rg.padj.is_nan() {
            rg.filtered = true;
            r.low_counts += 1;
        }
        if !rg.pvalue.is_nan() {
            r.tested += 1;
        }
        if rg.padj < opt.alpha {
            if rg.log2_fold_change > opt.lfc_threshold {
                r.up += 1;
            } else if rg.log2_fold_change < -opt.lfc_threshold {
                r.down += 1;
            }
        }
    }
    r.genes = res_genes;
    Ok(r)
}

// --- the files -----------------------------------------------------------

/// Splits a table line by tab, or by comma when the file is comma-separated;
/// surrounding double quotes are removed.
fn campos(t: &str, sep: char) -> Vec<String> {
    t.split(sep).map(|c| c.trim().trim_matches('"').to_string()).collect()
}

fn separador(cabecalho: &str) -> char {
    if cabecalho.contains('\t') {
        '\t'
    } else if cabecalho.contains(',') {
        ','
    } else {
        '\t'
    }
}

struct Amostras {
    nomes: Vec<String>,
    condicao: Vec<String>,
    fatores: Vec<(String, Vec<String>)>,
}

fn le_amostras(caminho: &str, opt: &DeOptions, avisos: &mut Vec<String>) -> Result<Amostras, String> {
    let mut l = Leitor::abre(caminho)?;
    let mut cab: Option<(Vec<String>, char)> = None;
    let mut linhas: Vec<(usize, Vec<String>)> = Vec::new();
    while l.avanca()? {
        let t = l.texto();
        if t.trim().is_empty() || t.starts_with('#') {
            continue;
        }
        match &cab {
            None => {
                let s = separador(t);
                cab = Some((campos(t, s), s));
            }
            Some((_, s)) => linhas.push((l.linha, campos(t, *s))),
        }
    }
    let Some((cab, _)) = cab else {
        return Err(format!("{caminho} is empty: a sample table has a header line and one line per sample (sample, condition, factors...)"));
    };
    let col = |nome: &str| -> Result<usize, String> {
        cab.iter().position(|c| c == nome).or_else(|| cab.iter().position(|c| c.eq_ignore_ascii_case(nome))).ok_or_else(|| {
            format!("{caminho} has no column \"{nome}\"; its columns are: {}", cab.join(", "))
        })
    };
    if cab.len() < 2 {
        return Err(format!(
            "{caminho}: the header has a single column; a sample table has the sample name in the first column and the condition (and factors) in others, separated by tabs or commas"
        ));
    }
    let ic = col(&opt.condition)?;
    let mut ifs = Vec::new();
    for f in &opt.factors {
        let i = col(f)?;
        if i == ic {
            return Err(format!("the factor {f} is the condition itself: list only the other factors"));
        }
        ifs.push((f.clone(), i));
    }
    let mut a = Amostras { nomes: Vec::new(), condicao: Vec::new(), fatores: ifs.iter().map(|(f, _)| (f.clone(), Vec::new())).collect() };
    let mut vistos: HashMap<String, usize> = HashMap::new();
    let tratado = opt.treated.as_deref().map(str::trim).filter(|t| !t.is_empty());
    if opt.only_compared_levels && tratado.is_none() {
        return Err("only_compared_levels needs the treated level: name the level compared with the reference".into());
    }
    let mut fora = 0;
    for (ln, c) in linhas {
        let pega = |i: usize, nome: &str| -> Result<String, String> {
            match c.get(i) {
                Some(v) if !v.is_empty() && v != "NA" => Ok(v.clone()),
                _ => Err(format!("{caminho}, line {ln}: no value for {nome}")),
            }
        };
        let s = pega(0, "the sample name (first column)")?;
        if let Some(l0) = vistos.insert(s.clone(), ln) {
            return Err(format!("{caminho}: sample \"{s}\" appears twice (lines {l0} and {ln})"));
        }
        let cond = pega(ic, &opt.condition)?;
        if opt.only_compared_levels {
            let e_ref = cond == opt.reference.trim();
            let e_trat = tratado.map_or(false, |t| cond == t);
            if !(e_ref || e_trat) {
                fora += 1;
                continue;
            }
        }
        a.nomes.push(s);
        a.condicao.push(cond);
        for (k, (f, i)) in ifs.iter().enumerate() {
            a.fatores[k].1.push(pega(*i, f)?);
        }
    }
    if fora > 0 {
        avisos.push(format!("{fora} sample(s) at other levels of {} were left out (only_compared_levels)", opt.condition));
    }
    if a.nomes.len() < 2 {
        return Err(format!("{caminho} lists {} sample(s): differential expression needs samples at two levels, with replicates", a.nomes.len()));
    }
    Ok(a)
}

struct Matriz {
    genes: Vec<String>,
    counts: Vec<f64>,
    ignoradas: Vec<String>,
    puladas: usize,
}

fn le_matriz(caminho: &str, amostras: &[String]) -> Result<Matriz, String> {
    let mut l = Leitor::abre(caminho)?;
    let mut cab: Option<(Vec<String>, char, usize)> = None;
    let mut idx: Vec<usize> = Vec::new();
    let mut m = Matriz { genes: Vec::new(), counts: Vec::new(), ignoradas: Vec::new(), puladas: 0 };
    let mut vistos: HashMap<String, usize> = HashMap::new();
    let n = amostras.len();
    while l.avanca()? {
        let t = l.texto();
        if t.trim().is_empty() || t.starts_with('#') {
            continue;
        }
        match &cab {
            None => {
                let s = separador(t);
                let c = campos(t, s);
                for a in amostras {
                    match c.iter().skip(1).position(|x| x == a) {
                        Some(i) => idx.push(i + 1),
                        None => {
                            return Err(format!(
                                "Sample \"{a}\" of the sample table is not a column of {caminho}; its columns are: {}",
                                lista_curta(&c[1..].to_vec())
                            ))
                        }
                    }
                }
                m.ignoradas = c.iter().skip(1).filter(|x| !amostras.contains(x)).cloned().collect();
                cab = Some((c, s, l.linha));
            }
            Some((cabec, s, _)) => {
                let c: Vec<&str> = t.split(*s).collect();
                let gene = c[0].trim().trim_matches('"');
                if gene.starts_with("__") {
                    m.puladas += 1;
                    continue;
                }
                if gene.is_empty() {
                    return Err(format!("{caminho}, line {}: no gene name in the first column", l.linha));
                }
                if let Some(l0) = vistos.insert(gene.to_string(), l.linha) {
                    return Err(format!("{caminho}: gene \"{gene}\" appears twice (lines {l0} and {}); each row must be a distinct gene", l.linha));
                }
                for (k, &i) in idx.iter().enumerate() {
                    let Some(v) = c.get(i) else {
                        return Err(format!(
                            "{caminho}, line {}: {} column(s), but the header has {}",
                            l.linha,
                            c.len(),
                            cabec.len()
                        ));
                    };
                    let v = v.trim().trim_matches('"');
                    let x: f64 = v.parse().map_err(|_| {
                        format!("{caminho}, line {}, column {} (sample {}): \"{v}\" is not a number", l.linha, i + 1, amostras[k])
                    })?;
                    if !(x >= 0.0) || x.fract() != 0.0 || !x.is_finite() {
                        return Err(format!(
                            "{caminho}, line {}, column {} (sample {}): {v} is not a non-negative whole number. The DESeq2 method models raw read counts; estimated counts (salmon, RSEM, kallisto) must be rounded to whole numbers first.",
                            l.linha,
                            i + 1,
                            amostras[k]
                        ));
                    }
                    m.counts.push(x);
                }
                m.genes.push(gene.to_string());
            }
        }
    }
    if cab.is_none() {
        return Err(format!("{caminho} is empty: a count matrix has a header (gene, then one column per sample) and one line per gene"));
    }
    if m.genes.is_empty() {
        return Err(format!("{caminho} has a header and no gene"));
    }
    debug_assert_eq!(m.counts.len(), m.genes.len() * n);
    Ok(m)
}

/// Gene id -> name, from NCBI's gene_info (GeneID -> Symbol) or a table of
/// two columns (id, name).
fn le_nomes(caminho: &str) -> Result<HashMap<String, String>, String> {
    let mut l = Leitor::abre(caminho)?;
    let mut mapa = HashMap::new();
    let mut gene_info = false;
    let mut primeira = true;
    while l.avanca()? {
        let t = l.texto();
        if primeira {
            primeira = false;
            if t.starts_with("#tax_id") || t.starts_with("#Format: tax_id") {
                gene_info = true;
                continue;
            }
        }
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let c: Vec<&str> = if t.contains('\t') { t.split('\t').collect() } else { t.split(',').collect() };
        let (id, nome) = if gene_info {
            match (c.get(1), c.get(2)) {
                (Some(a), Some(b)) => (*a, *b),
                _ => continue,
            }
        } else {
            match (c.first(), c.get(1)) {
                (Some(a), Some(b)) => (*a, *b),
                _ => {
                    return Err(format!(
                        "{caminho}, line {}: a table of gene names needs two columns (gene id, name), tab- or comma-separated",
                        l.linha
                    ))
                }
            }
        };
        let (id, nome) = (id.trim().trim_matches('"'), nome.trim().trim_matches('"'));
        if !id.is_empty() && !nome.is_empty() && nome != "-" {
            mapa.entry(id.to_string()).or_insert_with(|| nome.to_string());
        }
    }
    if mapa.is_empty() {
        return Err(format!("{caminho} gave no gene name: use NCBI's gene_info or a two-column table (gene id, name)"));
    }
    Ok(mapa)
}

/// The DESeq2 method from files: `counts` is a count matrix (tab- or
/// comma-separated, `.gz` accepted; first column the gene; featureCounts
/// and htseq-count tables work: columns not in the sample table are ignored,
/// `#` lines and `__` rows skipped); `samples` a table whose first column
/// names the samples (the matrix's column names) and whose other columns
/// hold the condition and the blocking factors.
///
/// `output` gets gene, (name,) baseMean, log2FoldChange, lfcSE, stat,
/// pvalue, padj, sorted by padj (NA last), as DESeq2's results table.
pub fn differential_expression(counts: &str, samples: &str, output: Option<&str>, opt: &DeOptions) -> Result<DeResult, String> {
    if let Some(o) = output {
        recusa_mesmo_arquivo(counts, o)?;
        recusa_mesmo_arquivo(samples, o)?;
    }
    if !(opt.alpha > 0.0 && opt.alpha < 1.0) {
        return Err(format!("alpha must be between 0 and 1, not {}", opt.alpha));
    }
    let mut avisos = Vec::new();
    let a = le_amostras(samples, opt, &mut avisos)?;
    let mtz = le_matriz(counts, &a.nomes)?;
    let nomes = match &opt.gene_names {
        Some(p) if !p.trim().is_empty() => Some(le_nomes(p)?),
        _ => None,
    };
    let mut r = deseq(&mtz.counts, &mtz.genes, &a.nomes, &a.condicao, &a.fatores, opt)?;
    r.counts = counts.to_string();
    r.sample_table = samples.to_string();
    r.ignored_columns = mtz.ignoradas;
    r.skipped_rows = mtz.puladas;
    avisos.append(&mut r.warnings);
    r.warnings = avisos;
    if let Some(mapa) = &nomes {
        for g in r.genes.iter_mut() {
            match mapa.get(&g.gene) {
                Some(nm) => g.name = nm.clone(),
                None => r.names_missing += 1,
            }
        }
        if r.names_missing > 0 {
            r.warnings.push(format!(
                "{} gene(s) are not in the name table (retired or unknown ids): their id is used as the name",
                milhar(r.names_missing)
            ));
        }
    }
    if let Some(saida) = output {
        let mut w = Escritor::cria(saida)?;
        let mut guarda = Saida { caminho: saida.to_string(), armada: true };
        let erro = |e: String| format!("Error writing {saida}: {e}");
        let com_nome = nomes.is_some();
        w.linha(if com_nome {
            "gene\tname\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj"
        } else {
            "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj"
        })
        .map_err(erro)?;
        for i in r.order() {
            let g = &r.genes[i as usize];
            let nome = if com_nome { format!("\t{}", g.name) } else { String::new() };
            w.linha(&format!(
                "{}{nome}\t{}\t{}\t{}\t{}\t{}\t{}",
                g.gene,
                num(g.base_mean),
                num(g.log2_fold_change),
                num(g.lfc_se),
                num(g.stat),
                num(g.pvalue),
                num(g.padj)
            ))
            .map_err(erro)?;
        }
        w.termina().map_err(erro)?;
        guarda.armada = false;
        r.output = Some(saida.to_string());
    }
    Ok(r)
}
