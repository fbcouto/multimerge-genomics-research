// =====================================================================
//  CLINICAL VARIANT TOOLS, FOR TEACHING AND RESEARCH
//
//  Four tools that sit between a person's variant calls and the public
//  databases:
//
//   - `clinvar_findings`  joins a person's VCF with the ClinVar VCF and lists
//                         the alleles the person carries that ClinVar has a
//                         record for: genes, classification, review stars,
//                         conditions, molecular consequence, zygosity;
//   - `gene_regions`      the coordinates of named genes in a GFF3 or GTF
//                         annotation, written as BED (a "virtual panel");
//   - `variant_lookup`    the genotype of several people (files and/or
//                         samples) at a list of variants -- the lookup done
//                         in cascade testing of relatives;
//   - `pileup`            the bases that the reads of a BAM (or SAM) show at
//                         one reference position, with filters -- the numbers
//                         a Sanger trace is usually compared with.
//
//  None of this is a diagnosis. The tools report database matches and read
//  counts, and the texts say what was found, never what it means for the
//  health of anyone.
//
//  House rules, as in the other modules:
//
//   - ANY SIZE. The person's VCF (~5 million records) and ClinVar (~4.6
//     million) are joined by a streaming merge that holds a few records at a
//     time; the GFF is streamed; only BED files of regions are loaded whole
//     (a callable-regions BED of a genome is a few MB in memory).
//   - Chromosome names are compared by their canonical form: chr1 == 1,
//     chrM == M == MT.
//   - Variants are compared after the same trimming on both sides: a
//     multiallelic record is split into one allele per ALT, then the bases
//     shared at the end and then at the start are removed (POS moves forward
//     with the start). There is NO left-alignment: an indel inside a repeat
//     only matches if both files place it the same way. ClinVar is
//     left-aligned; a person's VCF can be left-aligned with
//     vcf_normalize(..., reference="genome.fa").
//   - Every message says what is wrong, in which file (and line) and how to
//     fix it.
// =====================================================================

use crate::bam;
use crate::bgzf::BgzfReader;
use crate::externo::{abre_leitura, natural_cmp, Escritor, Formato};
use crate::indice;
use crate::nomes::{canonico, converte_nome, convencao, json_texto, Convencao, Estilo};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::rc::Rc;

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

    /// Reads the next line. `false` at the end of the file.
    fn avanca(&mut self) -> Result<bool, String> {
        self.buf.clear();
        match self.r.read_line(&mut self.buf) {
            Ok(0) => Ok(false),
            Ok(_) => {
                self.linha += 1;
                Ok(true)
            }
            Err(e) => Err(format!(
                "Error reading {} after line {}: {e}. The file may be truncated, binary (BCF is not read; convert it to VCF first) or not UTF-8 text.",
                self.caminho, self.linha
            )),
        }
    }

    /// The current line without its line ending (LF or CRLF).
    fn texto(&self) -> &str {
        self.buf.trim_end_matches(['\n', '\r'])
    }
}

/// Removes a half-written output when the operation ends in error. Armed
/// only after the output was created, so a file that could not be created is
/// never touched.
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

fn jopt(s: Option<&str>) -> String {
    s.map(json_texto).unwrap_or_else(|| "null".to_string())
}

fn jlista(v: &[String]) -> String {
    format!("[{}]", v.iter().map(|s| json_texto(s)).collect::<Vec<_>>().join(","))
}

fn jnum(x: f64) -> String {
    if x.is_finite() {
        format!("{x}")
    } else {
        "null".to_string()
    }
}

fn jnum_opt<T: std::fmt::Display>(x: Option<T>) -> String {
    x.map(|v| v.to_string()).unwrap_or_else(|| "null".to_string())
}

fn jpares(v: &[(String, usize)]) -> String {
    format!("[{}]", v.iter().map(|(n, q)| format!("[{},{q}]", json_texto(n))).collect::<Vec<_>>().join(","))
}

fn percent(x: usize, de: usize) -> String {
    if de == 0 {
        "-".to_string()
    } else {
        format!("{:.1}%", 100.0 * x as f64 / de as f64)
    }
}

// ---------------------------------------------------------------------
//  CHROMOSOMES
// ---------------------------------------------------------------------

const RANK_X: u32 = 1001;
const RANK_Y: u32 = 1002;
const RANK_MT: u32 = 1003;
const RANK_OUTRA: u32 = u32::MAX;

/// Where a sequence goes in the join: the primary chromosomes in the order
/// 1, 2, ..., 22, X, Y; the mitochondrion apart (sorters disagree on where it
/// goes, so it is joined in memory: a few thousand records); everything else
/// (alt contigs, unplaced scaffolds, decoys) is skipped.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Lugar {
    Primario(u32),
    Mito,
    Outra,
}

fn lugar(nome: &str) -> Lugar {
    let c = canonico(nome);
    match c.as_str() {
        "X" => Lugar::Primario(RANK_X),
        "Y" => Lugar::Primario(RANK_Y),
        "MT" => Lugar::Mito,
        s if !s.is_empty() && s.len() <= 3 && s.bytes().all(|b| b.is_ascii_digit()) => match s.parse::<u32>() {
            Ok(n) if n >= 1 => Lugar::Primario(n),
            _ => Lugar::Outra,
        },
        _ => Lugar::Outra,
    }
}

/// Sort key of a sequence name in the outputs: 1..22, X, Y, MT, then the rest.
fn chave_crom(nome: &str) -> u32 {
    match lugar(nome) {
        Lugar::Primario(r) => r,
        Lugar::Mito => RANK_MT,
        Lugar::Outra => RANK_OUTRA,
    }
}

fn compara_crom(a: &str, b: &str) -> Ordering {
    chave_crom(a).cmp(&chave_crom(b)).then_with(|| natural_cmp(a, b))
}

/// Which assembly a header line points to, when it says: `##reference=` or
/// the length of chromosome 1 in a `##contig` line.
fn montagem_da_linha(l: &str) -> Option<&'static str> {
    if let Some(r) = l.strip_prefix("##reference=").or_else(|| l.strip_prefix("##assembly=")) {
        let r = r.to_ascii_lowercase();
        if r.contains("grch38") || r.contains("hg38") {
            return Some("GRCh38");
        }
        if r.contains("grch37") || r.contains("hg19") || r.contains("b37") || r.contains("hs37") {
            return Some("GRCh37");
        }
        return None;
    }
    let r = l.strip_prefix("##contig=<")?;
    let id = r.split([',', '>']).find_map(|c| c.strip_prefix("ID="))?;
    if canonico(id) != "1" {
        return None;
    }
    match r.split([',', '>']).find_map(|c| c.strip_prefix("length="))? {
        "248956422" => Some("GRCh38"),
        "249250621" => Some("GRCh37"),
        _ => None,
    }
}

fn confere_montagens(a: &str, ma: Option<&str>, b: &str, mb: Option<&str>) -> Result<(), String> {
    if let (Some(x), Some(y)) = (ma, mb) {
        if x != y {
            return Err(format!(
                "{a} is on {x} and {b} is on {y}. The coordinates of the two assemblies differ, so nothing would match correctly. Use files on the same assembly (ClinVar publishes both: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/ and vcf_GRCh38/)."
            ));
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------
//  VCF RECORDS, ALLELES, GENOTYPES
// ---------------------------------------------------------------------

/// The columns of a VCF record, borrowed from the line.
struct Campos<'a> {
    chrom: &'a str,
    pos: u64,
    id: &'a str,
    ref_: &'a str,
    alt: &'a str,
    qual: &'a str,
    filter: &'a str,
    info: &'a str,
    format: Option<&'a str>,
    /// The sample columns, still joined by tabs.
    amostras: Option<&'a str>,
}

/// `None`: fewer than 8 columns, or POS is not a whole number.
fn campos(t: &str) -> Option<Campos<'_>> {
    let mut c = t.splitn(10, '\t');
    let chrom = c.next()?;
    let pos = c.next()?.trim().parse().ok()?;
    let id = c.next()?;
    let ref_ = c.next()?;
    let alt = c.next()?;
    let qual = c.next()?;
    let filter = c.next()?;
    let info = c.next()?;
    if chrom.is_empty() {
        return None;
    }
    Some(Campos { chrom, pos, id, ref_, alt, qual, filter, info, format: c.next(), amostras: c.next() })
}

fn info_valor<'a>(info: &'a str, chave: &str) -> Option<&'a str> {
    info.split(';').find_map(|kv| kv.strip_prefix(chave).and_then(|r| r.strip_prefix('=')))
}

/// The value of FORMAT key `chave` for the sample in column `col` (0 = the
/// first sample), as written (a missing GT is `.`, which is a no call).
/// `None` when the key or the value is not there.
fn campo_amostra<'a>(format: &str, amostras: &'a str, col: usize, chave: &str) -> Option<&'a str> {
    let i = format.split(':').position(|f| f == chave)?;
    let v = amostras.split('\t').nth(col)?.split(':').nth(i)?;
    if v.is_empty() {
        None
    } else {
        Some(v)
    }
}

fn valor_formato<'a>(c: &Campos<'a>, col: Option<usize>, chave: &str) -> Option<&'a str> {
    match (col, c.format, c.amostras) {
        (Some(j), Some(f), Some(s)) => campo_amostra(f, s, j, chave),
        _ => None,
    }
}

/// Like `valor_formato`, with `.` (missing) as None: for DP, GQ, AD.
fn valor_presente<'a>(c: &Campos<'a>, col: Option<usize>, chave: &str) -> Option<&'a str> {
    valor_formato(c, col, chave).filter(|v| *v != ".")
}

/// Non-empty and made only of A, C, G, T, N in either case.
fn so_bases(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}

/// The minimal form of one allele, in upper case: the bases shared at the end
/// are removed, then those shared at the start (moving POS forward), always
/// keeping at least one base on each side. `None`: not a sequence change
/// (symbolic `<DEL>`, `*`, breakend, IUPAC codes) or ALT equal to REF.
fn minimo(pos: u64, r: &str, a: &str) -> Option<(u64, String, String)> {
    if !so_bases(r) || !so_bases(a) {
        return None;
    }
    let (r, a) = (r.as_bytes(), a.as_bytes());
    let (mut fr, mut fa) = (r.len(), a.len());
    while fr > 1 && fa > 1 && r[fr - 1].eq_ignore_ascii_case(&a[fa - 1]) {
        fr -= 1;
        fa -= 1;
    }
    let mut i = 0usize;
    while fr - i > 1 && fa - i > 1 && r[i].eq_ignore_ascii_case(&a[i]) {
        i += 1;
    }
    let rr = String::from_utf8(r[i..fr].to_ascii_uppercase()).ok()?;
    let aa = String::from_utf8(a[i..fa].to_ascii_uppercase()).ok()?;
    if rr == aa {
        return None;
    }
    Some((pos + i as u64, rr, aa))
}

/// How many copies of one ALT allele a genotype holds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Zygosity {
    /// Some, but not all, of the called alleles (`0/1`, `1/2`, `1/.`).
    Het,
    /// Every allele of a diploid (or polyploid) call (`1/1`).
    Hom,
    /// The only allele of a haploid call (`1`, as on chrY or chrX of a male).
    Hemizygous,
    /// Every allele is the reference (`0/0`, `0`).
    HomRef,
    /// The call holds other ALT alleles of the record but not this one.
    OtherAlt,
    /// Not called (`./.`, or `0/.`).
    NoCall,
    /// No genotype to read: a sites-only VCF, or no GT in FORMAT.
    Unknown,
}

impl Zygosity {
    pub fn as_str(self) -> &'static str {
        match self {
            Zygosity::Het => "het",
            Zygosity::Hom => "hom",
            Zygosity::Hemizygous => "hemizygous",
            Zygosity::HomRef => "hom_ref",
            Zygosity::OtherAlt => "other_alt",
            Zygosity::NoCall => "no_call",
            Zygosity::Unknown => "unknown",
        }
    }

    /// Does the person carry the allele (or, without a genotype, list it)?
    fn carrega(self) -> bool {
        matches!(self, Zygosity::Het | Zygosity::Hom | Zygosity::Hemizygous | Zygosity::Unknown)
    }
}

/// The zygosity of ALT number `k` (1-based) in a GT value.
fn zigosidade(gt: Option<&str>, k: u32) -> Zygosity {
    let Some(gt) = gt else { return Zygosity::Unknown };
    let (mut ploidia, mut copias, mut faltam, mut outro) = (0usize, 0usize, 0usize, false);
    for a in gt.split(['/', '|']) {
        ploidia += 1;
        match a.trim().parse::<u32>() {
            Ok(x) if x == k => copias += 1,
            Ok(0) => {}
            Ok(_) => outro = true,
            Err(_) => faltam += 1,
        }
    }
    if ploidia == 0 {
        Zygosity::Unknown
    } else if copias > 0 {
        if ploidia == 1 {
            Zygosity::Hemizygous
        } else if copias == ploidia {
            Zygosity::Hom
        } else {
            Zygosity::Het
        }
    } else if outro {
        Zygosity::OtherAlt
    } else if faltam > 0 {
        Zygosity::NoCall
    } else {
        Zygosity::HomRef
    }
}

/// The sample column to read. `pedida`: the sample asked for.
fn escolhe_amostra(caminho: &str, nomes: &[String], pedida: Option<&str>) -> Result<Option<usize>, String> {
    let lista = || {
        let mut s: Vec<&str> = nomes.iter().take(10).map(String::as_str).collect();
        if nomes.len() > 10 {
            s.push("...");
        }
        s.join(", ")
    };
    match pedida {
        Some(p) => {
            if nomes.is_empty() {
                return Err(format!(
                    "{caminho} has no sample columns (a sites-only VCF), so sample \"{p}\" cannot be read. Leave the sample empty to take every listed allele."
                ));
            }
            nomes.iter().position(|n| n == p).map(Some).ok_or_else(|| {
                format!("Sample \"{p}\" is not in {caminho}. Its samples are: {}. Names are case-sensitive.", lista())
            })
        }
        None => match nomes.len() {
            0 => Ok(None),
            1 => Ok(Some(0)),
            n => Err(format!("{caminho} has {n} samples ({}); say which one with sample=\"...\".", lista())),
        },
    }
}

// ---------------------------------------------------------------------
//  CLINVAR: CLASSIFICATIONS AND REVIEW STARS
// ---------------------------------------------------------------------

/// Review stars of a ClinVar review status (CLNREVSTAT, ONCREVSTAT or
/// SCIREVSTAT), as on https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/:
///
///   4  practice guideline
///   3  reviewed by expert panel
///   2  criteria provided, multiple submitters, no conflicts
///      (and, for somatic classifications, "criteria provided, multiple submitters")
///   1  criteria provided, single submitter
///      criteria provided, conflicting classifications (formerly "interpretations")
///   0  no assertion criteria provided; no classification provided; no
///      classification for the single variant; no classifications from
///      unflagged records (and the older "no assertion provided", "no
///      interpretation for the single variant")
///
/// Underscores, commas and case are ignored. `None`: a value not in this list.
pub fn review_stars(status: &str) -> Option<u8> {
    let s: String = status.to_ascii_lowercase().replace(['_', ','], " ");
    let s = s.split_whitespace().collect::<Vec<_>>().join(" ");
    Some(match s.as_str() {
        "practice guideline" => 4,
        "reviewed by expert panel" => 3,
        "criteria provided multiple submitters no conflicts" | "criteria provided multiple submitters" => 2,
        "criteria provided single submitter"
        | "criteria provided conflicting classifications"
        | "criteria provided conflicting interpretations" => 1,
        "no assertion criteria provided"
        | "no classification provided"
        | "no assertion provided"
        | "no classification for the single variant"
        | "no classifications from unflagged records"
        | "no interpretation for the single variant"
        | "no assertion for the individual variant" => 0,
        _ => return None,
    })
}

/// The class names a caller can ask for, one per ClinVar term (the VUS-low,
/// -mid and -high terms count as `uncertain`).
pub const CLASSES: [&str; 21] = [
    "pathogenic",
    "likely_pathogenic",
    "pathogenic_low_penetrance",
    "likely_pathogenic_low_penetrance",
    "uncertain",
    "conflicting",
    "likely_benign",
    "benign",
    "drug_response",
    "risk_factor",
    "established_risk_allele",
    "likely_risk_allele",
    "uncertain_risk_allele",
    "association",
    "association_not_found",
    "protective",
    "affects",
    "confers_sensitivity",
    "other",
    "not_provided",
    "not_classified",
];

const PLP: [&str; 4] = ["pathogenic", "likely_pathogenic", "pathogenic_low_penetrance", "likely_pathogenic_low_penetrance"];

/// The categories of the summary, in the order shown.
pub const CATEGORIES: [&str; 6] = [
    "pathogenic_or_likely_pathogenic",
    "conflicting",
    "uncertain",
    "benign_or_likely_benign",
    "other",
    "not_classified",
];

fn etiqueta(componente: &str) -> &'static str {
    let c = componente.trim().to_ascii_lowercase().replace(",_", "_").replace(", ", "_").replace(' ', "_");
    match c.as_str() {
        "pathogenic" => "pathogenic",
        "likely_pathogenic" => "likely_pathogenic",
        "pathogenic_low_penetrance" => "pathogenic_low_penetrance",
        "likely_pathogenic_low_penetrance" => "likely_pathogenic_low_penetrance",
        "uncertain_significance" | "vus-high" | "vus-mid" | "vus-low" | "vus" => "uncertain",
        "conflicting_classifications_of_pathogenicity"
        | "conflicting_interpretations_of_pathogenicity"
        | "conflicting_data_from_submitters" => "conflicting",
        "likely_benign" => "likely_benign",
        "benign" => "benign",
        "drug_response" => "drug_response",
        "risk_factor" => "risk_factor",
        "established_risk_allele" => "established_risk_allele",
        "likely_risk_allele" => "likely_risk_allele",
        "uncertain_risk_allele" => "uncertain_risk_allele",
        "association" => "association",
        "association_not_found" => "association_not_found",
        "protective" => "protective",
        "affects" => "affects",
        "confers_sensitivity" => "confers_sensitivity",
        "not_provided" => "not_provided",
        "no_classification_for_the_single_variant"
        | "no_classifications_from_unflagged_records"
        | "no_interpretation_for_the_single_variant" => "not_classified",
        _ => "other",
    }
}

/// The class tags of a CLNSIG value: every term between `/` and `|`
/// ("Pathogenic/Likely_pathogenic|risk_factor" has three). No CLNSIG (an
/// "included" variant, or a somatic-only record) is `not_classified`.
pub fn significance_classes(clnsig: Option<&str>) -> Vec<&'static str> {
    let mut v: Vec<&'static str> = Vec::new();
    if let Some(s) = clnsig {
        for parte in s.split(['|', '/']) {
            if parte.trim().is_empty() {
                continue;
            }
            let t = etiqueta(parte);
            if !v.contains(&t) {
                v.push(t);
            }
        }
    }
    if v.is_empty() {
        v.push("not_classified");
    }
    v
}

/// One category per record, the first that applies: P/LP (any of the four
/// pathogenic terms), conflicting, uncertain, benign or likely benign,
/// not classified (only `not_provided` / no classification), other.
pub fn significance_category(classes: &[&str]) -> &'static str {
    if classes.iter().any(|c| PLP.contains(c)) {
        "pathogenic_or_likely_pathogenic"
    } else if classes.contains(&"conflicting") {
        "conflicting"
    } else if classes.contains(&"uncertain") {
        "uncertain"
    } else if classes.iter().any(|c| *c == "benign" || *c == "likely_benign") {
        "benign_or_likely_benign"
    } else if classes.iter().all(|c| *c == "not_classified" || *c == "not_provided") {
        "not_classified"
    } else {
        "other"
    }
}

/// The classes a caller asked for. Empty: every class. Accepts the names in
/// `CLASSES`, in any case, with spaces or hyphens for underscores, and the
/// groups `plp` (the four pathogenic terms), `blb` (benign, likely benign),
/// `vus` / `uncertain_significance`, `risk` (risk factor and risk alleles)
/// and `all`.
fn classes_pedidas(v: &[String]) -> Result<Option<HashSet<&'static str>>, String> {
    if v.is_empty() {
        return Ok(None);
    }
    let mut s: HashSet<&'static str> = HashSet::new();
    for t in v {
        let n = t.trim().to_ascii_lowercase().replace([' ', '-'], "_").replace(",_", "_");
        match n.as_str() {
            "" => continue,
            "all" => return Ok(None),
            "plp" | "p/lp" | "p_lp" | "pathogenic_or_likely_pathogenic" | "pathogenic/likely_pathogenic" => s.extend(PLP),
            "blb" | "b/lb" | "b_lb" | "benign_or_likely_benign" | "benign/likely_benign" => {
                s.extend(["benign", "likely_benign"])
            }
            "vus" | "uncertain_significance" => {
                s.insert("uncertain");
            }
            "conflicting_classifications_of_pathogenicity" | "conflicting_interpretations_of_pathogenicity" => {
                s.insert("conflicting");
            }
            "risk" => s.extend(["risk_factor", "established_risk_allele", "likely_risk_allele", "uncertain_risk_allele"]),
            x => match CLASSES.iter().find(|c| **c == x) {
                Some(c) => {
                    s.insert(*c);
                }
                None => {
                    return Err(format!(
                        "Unknown ClinVar class \"{t}\". Use any of: {}; or the groups plp, blb, vus, risk, all.",
                        CLASSES.join(", ")
                    ))
                }
            },
        }
    }
    Ok(Some(s))
}

fn espacos(s: &str) -> String {
    s.replace('_', " ")
}

/// What a ClinVar record says, from its ID and INFO.
#[derive(Clone, Debug, Default)]
struct ClinVar {
    variation_id: String,
    allele_id: Option<String>,
    rsids: Vec<String>,
    genes: Vec<String>,
    significance: String,
    classes: Vec<&'static str>,
    category: &'static str,
    review_status: String,
    review_raw: Option<String>,
    stars: u8,
    conditions: Vec<String>,
    consequences: Vec<String>,
    variant_type: Option<String>,
    oncogenicity: Option<String>,
    oncogenicity_stars: Option<u8>,
    somatic_impact: Option<String>,
    somatic_impact_stars: Option<u8>,
}

fn le_clinvar(id: &str, info: &str) -> ClinVar {
    let mut cv = ClinVar { variation_id: id.to_string(), ..Default::default() };
    let (mut sig, mut onc_rev, mut sci_rev) = (None, None, None);
    for kv in info.split(';') {
        let (k, v) = kv.split_once('=').unwrap_or((kv, ""));
        match k {
            "ALLELEID" => cv.allele_id = Some(v.to_string()),
            "RS" => {
                cv.rsids = v.split(',').filter(|x| !x.is_empty() && *x != ".").map(|x| format!("rs{x}")).collect();
            }
            "GENEINFO" => {
                cv.genes =
                    v.split('|').filter_map(|g| g.split(':').next()).filter(|g| !g.is_empty()).map(str::to_string).collect();
            }
            "CLNSIG" => sig = Some(v),
            "CLNREVSTAT" => cv.review_raw = Some(v.to_string()),
            "CLNDN" => cv.conditions = v.split('|').filter(|x| !x.is_empty()).map(espacos).collect(),
            "MC" => {
                cv.consequences = v
                    .split(',')
                    .filter_map(|m| m.split_once('|').map(|x| x.1))
                    .filter(|x| !x.is_empty())
                    .map(str::to_string)
                    .collect();
            }
            "CLNVC" => cv.variant_type = Some(espacos(v)),
            "ONC" => cv.oncogenicity = Some(espacos(v)),
            "ONCREVSTAT" => onc_rev = Some(v),
            "SCI" => cv.somatic_impact = Some(espacos(v)),
            "SCIREVSTAT" => sci_rev = Some(v),
            _ => {}
        }
    }
    cv.significance = sig.map(espacos).unwrap_or_else(|| "not classified".to_string());
    cv.classes = significance_classes(sig);
    cv.category = significance_category(&cv.classes);
    cv.review_status = cv.review_raw.as_deref().map(espacos).unwrap_or_default();
    cv.stars = cv.review_raw.as_deref().and_then(review_stars).unwrap_or(0);
    // An unknown status counts as 0 stars, as for the germline classification.
    cv.oncogenicity_stars = onc_rev.map(|r| review_stars(r).unwrap_or(0));
    cv.somatic_impact_stars = sci_rev.map(|r| review_stars(r).unwrap_or(0));
    cv
}

// ---------------------------------------------------------------------
//  REGIONS (BED), LOADED WHOLE
// ---------------------------------------------------------------------

/// Regions of a BED file per canonical sequence name: 0-based, half-open,
/// sorted and merged (touching intervals too).
struct Intervalos {
    por_crom: HashMap<String, Vec<(u64, u64)>>,
}

impl Intervalos {
    /// `so`: keep only these canonical sequence names (the rest is not even
    /// stored).
    fn le(caminho: &str, so: Option<&HashSet<String>>) -> Result<Intervalos, String> {
        let mut l = Leitor::abre(caminho)?;
        let mut m: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        let (mut ultimo, mut ultimo_canon) = (String::new(), String::new());
        while l.avanca()? {
            let t = l.texto();
            if t.trim().is_empty() || t.starts_with('#') || t.starts_with("track") || t.starts_with("browser") {
                continue;
            }
            let mut c = t.split('\t');
            let (Some(nome), Some(a), Some(b)) = (c.next(), c.next(), c.next()) else {
                return Err(format!(
                    "{caminho}, line {}: a BED line needs at least 3 tab-separated columns (chrom, start, end)",
                    l.linha
                ));
            };
            let (Ok(a), Ok(b)) = (a.trim().parse::<u64>(), b.trim().parse::<u64>()) else {
                return Err(format!("{caminho}, line {}: start and end must be whole numbers (0-based BED coordinates)", l.linha));
            };
            if b < a {
                return Err(format!("{caminho}, line {}: the end ({b}) is before the start ({a})", l.linha));
            }
            if nome != ultimo {
                ultimo = nome.to_string();
                ultimo_canon = canonico(nome);
            }
            if so.is_some_and(|f| !f.contains(&ultimo_canon)) {
                continue;
            }
            m.entry(ultimo_canon.clone()).or_default().push((a, b));
        }
        for v in m.values_mut() {
            v.sort_unstable();
            let mut unidos: Vec<(u64, u64)> = Vec::with_capacity(v.len());
            for &(a, b) in v.iter() {
                match unidos.last_mut() {
                    Some(u) if a <= u.1 => u.1 = u.1.max(b),
                    _ => unidos.push((a, b)),
                }
            }
            *v = unidos;
        }
        Ok(Intervalos { por_crom: m })
    }

    /// Is `[ini, fim)` entirely inside one region?
    fn contem(&self, canon: &str, ini: u64, fim: u64) -> bool {
        let Some(v) = self.por_crom.get(canon) else { return false };
        let i = v.partition_point(|&(_, e)| e <= ini);
        i < v.len() && v[i].0 <= ini && fim <= v[i].1
    }

    /// Does `[ini, fim)` overlap a region? (An empty span counts as one base.)
    fn toca(&self, canon: &str, ini: u64, fim: u64) -> bool {
        let Some(v) = self.por_crom.get(canon) else { return false };
        let fim = fim.max(ini + 1);
        let i = v.partition_point(|&(_, e)| e <= ini);
        i < v.len() && v[i].0 < fim
    }
}

// ---------------------------------------------------------------------
//  SORTED STREAM OF ALLELES
//
//  Splitting and trimming can move an allele forward (never backward), so a
//  file sorted by record POS is not always sorted by allele position. The
//  alleles wait in a small heap until the reader has passed their position:
//  an allele at P is final once a record with POS > P has been read, because
//  every later allele is at or after its record's POS.
// ---------------------------------------------------------------------

struct Reg {
    linha: String,
    numero: usize,
    /// Person files: the sample column read.
    amostra: Option<usize>,
}

struct Alelo {
    pos: u64,
    ref_: String,
    alt: String,
    zig: Zygosity,
    /// Order of arrival, to release equal positions in file order.
    seq: u64,
    reg: Rc<Reg>,
}

struct Pendente(Alelo);

impl PartialEq for Pendente {
    fn eq(&self, o: &Self) -> bool {
        self.cmp(o) == Ordering::Equal
    }
}
impl Eq for Pendente {}
impl PartialOrd for Pendente {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Pendente {
    /// Reversed: the heap of the standard library is a max-heap.
    fn cmp(&self, o: &Self) -> Ordering {
        (o.0.pos, o.0.seq).cmp(&(self.0.pos, self.0.seq))
    }
}

struct Fluxo {
    l: Leitor,
    /// A person's VCF (alleles not carried by the sample are dropped) or the
    /// ClinVar VCF (every allele, and the header must declare CLNSIG).
    pessoa: bool,
    pedida: Option<String>,
    amostra: Option<usize>,
    amostra_nome: Option<String>,
    decidida: bool,
    viu_clnsig: bool,
    inicio_dados: bool,
    montagem: Option<&'static str>,
    heap: BinaryHeap<Pendente>,
    rank: Option<u32>,
    crom_atual: String,
    ultima_pos: u64,
    segurada: Option<(String, usize)>,
    fim: bool,
    seq: u64,
    cache_nome: String,
    cache_canon: String,
    cache_lugar: Lugar,
    mito: Vec<Alelo>,
    /// Person files: only alleles inside these regions are compared.
    chamaveis: Option<Intervalos>,
    registros: usize,
    alelos: usize,
    nao_portados: usize,
    fora: usize,
    outras: usize,
    ignoradas: usize,
    nao_comparaveis: usize,
}

impl Fluxo {
    fn novo(caminho: &str, pessoa: bool, pedida: Option<String>) -> Result<Fluxo, String> {
        Ok(Fluxo {
            l: Leitor::abre(caminho)?,
            pessoa,
            pedida,
            amostra: None,
            amostra_nome: None,
            decidida: false,
            viu_clnsig: false,
            inicio_dados: false,
            montagem: None,
            heap: BinaryHeap::new(),
            rank: None,
            crom_atual: String::new(),
            ultima_pos: 0,
            segurada: None,
            fim: false,
            seq: 0,
            cache_nome: String::new(),
            cache_canon: String::new(),
            cache_lugar: Lugar::Outra,
            mito: Vec::new(),
            chamaveis: None,
            registros: 0,
            alelos: 0,
            nao_portados: 0,
            fora: 0,
            outras: 0,
            ignoradas: 0,
            nao_comparaveis: 0,
        })
    }

    /// The next allele on 1..22, X, Y, in order, with the rank of its
    /// chromosome. `None` at the end of the file.
    fn proximo(&mut self) -> Result<Option<(u32, Alelo)>, String> {
        loop {
            if let Some(topo) = self.heap.peek() {
                if self.fim || self.segurada.is_some() || topo.0.pos < self.ultima_pos {
                    let a = self.heap.pop().expect("peeked").0;
                    return Ok(Some((self.rank.unwrap_or(0), a)));
                }
            } else if let Some((linha, n)) = self.segurada.take() {
                self.processa(&linha, n)?;
                continue;
            } else if self.fim {
                return Ok(None);
            }
            if !self.l.avanca()? {
                self.fim = true;
                continue;
            }
            let n = self.l.linha;
            let buf = std::mem::take(&mut self.l.buf);
            let r = self.processa(buf.trim_end_matches(['\n', '\r']), n);
            self.l.buf = buf;
            r?;
        }
    }

    /// Reads the whole header (and nothing else), so the sample and the
    /// assembly are known before the join starts.
    fn le_cabecalho(&mut self) -> Result<(), String> {
        while self.segurada.is_none() && !self.fim {
            if !self.l.avanca()? {
                self.fim = true;
                break;
            }
            let n = self.l.linha;
            let buf = std::mem::take(&mut self.l.buf);
            let t = buf.trim_end_matches(['\n', '\r']);
            let r = if t.starts_with('#') || t.is_empty() {
                self.processa(t, n)
            } else {
                // The first record: kept to be processed by `proximo`.
                self.segurada = Some((t.to_string(), n));
                Ok(())
            };
            self.l.buf = buf;
            r?;
        }
        Ok(())
    }

    fn cabecalho(&mut self, t: &str) -> Result<(), String> {
        if self.montagem.is_none() {
            self.montagem = montagem_da_linha(t);
        }
        if t.starts_with("##INFO=<ID=CLNSIG,") {
            self.viu_clnsig = true;
        }
        if t.starts_with("#CHROM") && self.pessoa && !self.decidida {
            let nomes: Vec<String> = t.split('\t').skip(9).map(str::to_string).collect();
            self.amostra = escolhe_amostra(&self.l.caminho, &nomes, self.pedida.as_deref())?;
            self.amostra_nome = self.amostra.map(|i| nomes[i].clone());
            self.decidida = true;
        }
        Ok(())
    }

    fn primeiro_registro(&mut self) -> Result<(), String> {
        if self.pessoa && !self.decidida {
            // No #CHROM line: the samples have no names.
            if let Some(p) = &self.pedida {
                return Err(format!(
                    "{} has no #CHROM header line, so its samples have no names and sample \"{p}\" cannot be found. Add the VCF header.",
                    self.l.caminho
                ));
            }
            self.amostra = Some(0);
            self.decidida = true;
        }
        if !self.pessoa && !self.viu_clnsig {
            return Err(format!(
                "{} does not look like the ClinVar VCF: its header does not declare INFO/CLNSIG. Download clinvar.vcf.gz for your assembly from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/ (or vcf_GRCh37/) and check that the person's VCF and ClinVar were not passed in each other's place.",
                self.l.caminho
            ));
        }
        Ok(())
    }

    fn fora_de_ordem(&self, c: &Campos, n: usize) -> String {
        format!(
            "{}, line {n}: {}:{} comes after {}:{}. This tool reads the files side by side and needs them sorted by position, with the chromosomes in the order 1, 2, ..., 22, X, Y (chrM/MT may be anywhere; other sequences are skipped). Sort the file first: sort_genomic([\"{}\"], \"sorted.vcf.gz\", chrom_order=\"natural\")",
            self.l.caminho, c.chrom, c.pos, self.crom_atual, self.ultima_pos, self.l.caminho
        )
    }

    fn processa(&mut self, texto: &str, n: usize) -> Result<(), String> {
        if texto.is_empty() {
            return Ok(());
        }
        if texto.starts_with('#') {
            return self.cabecalho(texto);
        }
        let Some(c) = campos(texto) else {
            self.ignoradas += 1;
            return Ok(());
        };
        if !self.inicio_dados {
            self.inicio_dados = true;
            self.primeiro_registro()?;
        }
        let lugar_ = if c.chrom == self.cache_nome {
            self.cache_lugar
        } else {
            let l = lugar(c.chrom);
            self.cache_nome = c.chrom.to_string();
            self.cache_canon = canonico(c.chrom);
            self.cache_lugar = l;
            l
        };
        match lugar_ {
            Lugar::Outra => {
                self.registros += 1;
                self.outras += 1;
            }
            Lugar::Mito => {
                self.registros += 1;
                self.empurra(&c, texto, n, true);
            }
            Lugar::Primario(r) => {
                match self.rank {
                    Some(a) if r < a => return Err(self.fora_de_ordem(&c, n)),
                    Some(a) if r == a => {
                        if c.pos < self.ultima_pos {
                            return Err(self.fora_de_ordem(&c, n));
                        }
                    }
                    _ => {
                        // A new chromosome: the alleles of the previous one
                        // leave first.
                        if !self.heap.is_empty() {
                            self.segurada = Some((texto.to_string(), n));
                            return Ok(());
                        }
                        self.rank = Some(r);
                        self.crom_atual = c.chrom.to_string();
                    }
                }
                self.ultima_pos = c.pos;
                self.registros += 1;
                self.empurra(&c, texto, n, false);
            }
        }
        Ok(())
    }

    fn empurra(&mut self, c: &Campos, texto: &str, n: usize, mito: bool) {
        let gt = if self.pessoa { valor_formato(c, self.amostra, "GT") } else { None };
        let mut reg: Option<Rc<Reg>> = None;
        for (i, alt) in c.alt.split(',').enumerate() {
            let k = i as u32 + 1;
            let Some((pos, r, a)) = minimo(c.pos, c.ref_, alt) else {
                if alt != "." && !alt.is_empty() {
                    self.nao_comparaveis += 1;
                }
                continue;
            };
            let zig = if self.pessoa && self.amostra.is_some() { zigosidade(gt, k) } else { Zygosity::Unknown };
            if !zig.carrega() {
                self.nao_portados += 1;
                continue;
            }
            if let Some(b) = &self.chamaveis {
                if !b.contem(&self.cache_canon, pos - 1, pos - 1 + r.len() as u64) {
                    self.fora += 1;
                    continue;
                }
            }
            self.alelos += 1;
            self.seq += 1;
            let amostra = self.amostra;
            let reg = reg.get_or_insert_with(|| Rc::new(Reg { linha: texto.to_string(), numero: n, amostra })).clone();
            let al = Alelo { pos, ref_: r, alt: a, zig, seq: self.seq, reg };
            if mito {
                self.mito.push(al);
            } else {
                self.heap.push(Pendente(al));
            }
        }
    }
}

/// Joins two sorted streams by (chromosome, position, REF, ALT) and calls
/// `casou(person allele, ClinVar allele)` for every match. Both files are read
/// to the end (the mitochondrion is joined last, in memory), so the counters
/// cover the whole files.
fn junta(p: &mut Fluxo, c: &mut Fluxo, casou: &mut dyn FnMut(&Alelo, &Alelo) -> Result<(), String>) -> Result<(), String> {
    let mut a = p.proximo()?;
    let mut b = c.proximo()?;
    let (mut ga, mut gb): (Vec<Alelo>, Vec<Alelo>) = (Vec::new(), Vec::new());
    loop {
        let (ka, kb) = match (&a, &b) {
            (Some((ra, x)), Some((rb, y))) => ((*ra, x.pos), (*rb, y.pos)),
            _ => break,
        };
        match ka.cmp(&kb) {
            Ordering::Less => a = p.proximo()?,
            Ordering::Greater => b = c.proximo()?,
            Ordering::Equal => {
                ga.clear();
                gb.clear();
                while let Some((r, x)) = a.take() {
                    if (r, x.pos) != ka {
                        a = Some((r, x));
                        break;
                    }
                    ga.push(x);
                    a = p.proximo()?;
                }
                while let Some((r, y)) = b.take() {
                    if (r, y.pos) != kb {
                        b = Some((r, y));
                        break;
                    }
                    gb.push(y);
                    b = c.proximo()?;
                }
                for x in &ga {
                    for y in &gb {
                        if x.ref_ == y.ref_ && x.alt == y.alt {
                            casou(x, y)?;
                        }
                    }
                }
            }
        }
    }
    while p.proximo()?.is_some() {}
    while c.proximo()?.is_some() {}
    let mut mapa: HashMap<(u64, &str, &str), Vec<&Alelo>> = HashMap::new();
    for y in &c.mito {
        mapa.entry((y.pos, y.ref_.as_str(), y.alt.as_str())).or_default().push(y);
    }
    for x in &p.mito {
        if let Some(v) = mapa.get(&(x.pos, x.ref_.as_str(), x.alt.as_str())) {
            for y in v {
                casou(x, y)?;
            }
        }
    }
    Ok(())
}

// =====================================================================
//  1) CLINVAR FINDINGS
// =====================================================================

/// A person: a VCF, the sample to read (needed when the file has several),
/// and optionally a BED of the regions where the calls are reliable (the
/// "callable" or benchmark regions; GIAB files have calls outside them).
/// `clinvar_findings` then checks only the alleles entirely inside them;
/// `variant_lookup` uses them to tell "reference" from "not callable" and
/// flags each call with `in_callable`.
#[derive(Clone, Debug, Default)]
pub struct Person {
    pub vcf: String,
    pub sample: Option<String>,
    pub callable: Option<String>,
}

#[derive(Clone, Debug, Default)]
pub struct FindingsOptions {
    /// ClinVar classes to keep (see `CLASSES` and the groups of
    /// `classes_pedidas`). Empty: all.
    pub classes: Vec<String>,
    /// Minimum review stars (0 to 4) of the germline classification.
    pub min_stars: u8,
    /// Keep only variants in these genes (GENEINFO symbols, any case). Empty: all.
    pub genes: Vec<String>,
    /// Keep only variants that overlap the regions of this BED (a panel).
    pub regions: Option<String>,
    /// Also write the findings as TSV (`.gz` gives BGZF).
    pub output: Option<String>,
}

/// One allele the person carries that matches a ClinVar record.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Finding {
    /// As named in the person's VCF.
    pub chrom: String,
    /// Minimal form (shared bases trimmed), 1-based: the form ClinVar uses.
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    /// GENEINFO symbols.
    pub genes: Vec<String>,
    pub rsid: Option<String>,
    pub variation_id: String,
    pub allele_id: Option<String>,
    /// The GT of the sample, as written.
    pub genotype: String,
    /// het, hom, hemizygous, or unknown (no genotype in the file).
    pub zygosity: String,
    pub depth: Option<u64>,
    pub genotype_quality: Option<f64>,
    /// CLNSIG with underscores turned into spaces.
    pub significance: String,
    /// The class tags of `significance` (see `CLASSES`).
    pub classes: Vec<String>,
    /// One of `CATEGORIES`.
    pub category: String,
    pub review_status: String,
    pub stars: u8,
    /// CLNDN, one entry per condition, underscores turned into spaces.
    pub conditions: Vec<String>,
    /// MC: Sequence Ontology terms (missense_variant, ...).
    pub consequences: Vec<String>,
    pub variant_type: Option<String>,
    /// Somatic classifications (ONC, SCI) and their stars, when present.
    pub oncogenicity: Option<String>,
    pub oncogenicity_stars: Option<u8>,
    pub somatic_impact: Option<String>,
    pub somatic_impact_stars: Option<u8>,
    /// Line numbers of the two records, for checking by hand.
    pub person_line: usize,
    pub clinvar_line: usize,
}

const CABECALHO_TSV: &str = "#chrom\tpos\tref\talt\tgenes\tzygosity\tgenotype\tsignificance\tstars\treview_status\tconditions\tconsequences\tvariant_type\trsid\tvariation_id\tallele_id\tdepth\tgenotype_quality\toncogenicity\tsomatic_impact";

impl Finding {
    pub fn json(&self) -> String {
        format!(
            "{{\"chrom\":{},\"pos\":{},\"ref\":{},\"alt\":{},\"genes\":{},\"rsid\":{},\"variation_id\":{},\"allele_id\":{},\"genotype\":{},\"zygosity\":{},\"depth\":{},\"genotype_quality\":{},\"significance\":{},\"classes\":{},\"category\":{},\"review_status\":{},\"stars\":{},\"conditions\":{},\"consequences\":{},\"variant_type\":{},\"oncogenicity\":{},\"oncogenicity_stars\":{},\"somatic_impact\":{},\"somatic_impact_stars\":{},\"person_line\":{},\"clinvar_line\":{}}}",
            json_texto(&self.chrom),
            self.pos,
            json_texto(&self.ref_allele),
            json_texto(&self.alt_allele),
            jlista(&self.genes),
            jopt(self.rsid.as_deref()),
            json_texto(&self.variation_id),
            jopt(self.allele_id.as_deref()),
            json_texto(&self.genotype),
            json_texto(&self.zygosity),
            jnum_opt(self.depth),
            self.genotype_quality.map(jnum).unwrap_or_else(|| "null".into()),
            json_texto(&self.significance),
            jlista(&self.classes),
            json_texto(&self.category),
            json_texto(&self.review_status),
            self.stars,
            jlista(&self.conditions),
            jlista(&self.consequences),
            jopt(self.variant_type.as_deref()),
            jopt(self.oncogenicity.as_deref()),
            jnum_opt(self.oncogenicity_stars),
            jopt(self.somatic_impact.as_deref()),
            jnum_opt(self.somatic_impact_stars),
            self.person_line,
            self.clinvar_line
        )
    }

    pub fn tsv(&self) -> String {
        let o = |s: Option<&str>| s.filter(|x| !x.is_empty()).unwrap_or(".").to_string();
        let l = |v: &[String], sep: &str| if v.is_empty() { ".".to_string() } else { v.join(sep) };
        [
            self.chrom.clone(),
            self.pos.to_string(),
            self.ref_allele.clone(),
            self.alt_allele.clone(),
            l(&self.genes, ","),
            self.zygosity.clone(),
            o(Some(self.genotype.as_str())),
            self.significance.clone(),
            self.stars.to_string(),
            o(Some(self.review_status.as_str())),
            l(&self.conditions, "; "),
            l(&self.consequences, ","),
            o(self.variant_type.as_deref()),
            o(self.rsid.as_deref()),
            self.variation_id.clone(),
            o(self.allele_id.as_deref()),
            self.depth.map(|d| d.to_string()).unwrap_or_else(|| ".".into()),
            self.genotype_quality.map(|d| d.to_string()).unwrap_or_else(|| ".".into()),
            o(self.oncogenicity.as_deref()),
            o(self.somatic_impact.as_deref()),
        ]
        .join("\t")
    }
}

/// What `clinvar_findings` found.
#[derive(Clone, Debug, Default)]
pub struct FindingsReport {
    pub person: String,
    pub clinvar: String,
    /// The sample read (None: a sites-only VCF).
    pub sample: Option<String>,
    pub person_records: usize,
    /// Alleles the sample carries (or that a sites-only file lists) on
    /// 1..22, X, Y, MT: the person sites checked against ClinVar.
    pub person_alleles: usize,
    /// Alleles of the records that the sample does not carry (0/0, no call,
    /// another ALT).
    pub person_not_carried: usize,
    /// Carried alleles left out because they are not entirely inside the
    /// person's callable regions (when `Person::callable` is given).
    pub person_outside_callable: usize,
    pub callable: Option<String>,
    /// Records on other sequences (alt contigs, unplaced), not checked.
    pub person_other_sequences: usize,
    /// Symbolic, `*` and breakend alleles, not checked.
    pub person_not_comparable: usize,
    /// Lines with fewer than 8 columns or a non-numeric POS.
    pub person_skipped: usize,
    pub clinvar_records: usize,
    pub clinvar_alleles: usize,
    pub clinvar_other_sequences: usize,
    pub clinvar_skipped: usize,
    /// Carried alleles found in ClinVar, before the filters.
    pub matched: usize,
    /// Over `matched`: CLNSIG as written (underscores as spaces), most common first.
    pub by_significance: Vec<(String, usize)>,
    /// Over `matched`: the six `CATEGORIES`, in that order.
    pub by_category: Vec<(String, usize)>,
    /// The matches that pass the filters, in genome order.
    pub findings: Vec<Finding>,
    /// Among `findings`, the pathogenic / likely pathogenic ones by zygosity.
    pub plp_heterozygous: usize,
    pub plp_homozygous: usize,
    pub plp_hemizygous: usize,
    pub plp_unknown_zygosity: usize,
    /// CLNREVSTAT values not in the star table (counted as 0 stars).
    pub unknown_review_status: Vec<String>,
    /// Lines written to the TSV (findings, without the header).
    pub written: Option<usize>,
    pub filters: String,
    pub assembly: Option<String>,
}

impl FindingsReport {
    pub fn json(&self) -> String {
        let achados: Vec<String> = self.findings.iter().map(Finding::json).collect();
        format!(
            "{{\"person\":{},\"clinvar\":{},\"sample\":{},\"callable\":{},\"person_records\":{},\"person_alleles\":{},\"person_not_carried\":{},\"person_outside_callable\":{},\"person_other_sequences\":{},\"person_not_comparable\":{},\"person_skipped\":{},\"clinvar_records\":{},\"clinvar_alleles\":{},\"clinvar_other_sequences\":{},\"clinvar_skipped\":{},\"matched\":{},\"by_significance\":{},\"by_category\":{},\"findings\":[{}],\"plp_heterozygous\":{},\"plp_homozygous\":{},\"plp_hemizygous\":{},\"plp_unknown_zygosity\":{},\"unknown_review_status\":{},\"written\":{},\"filters\":{},\"assembly\":{},\"text\":{}}}",
            json_texto(&self.person),
            json_texto(&self.clinvar),
            jopt(self.sample.as_deref()),
            jopt(self.callable.as_deref()),
            self.person_records,
            self.person_alleles,
            self.person_not_carried,
            self.person_outside_callable,
            self.person_other_sequences,
            self.person_not_comparable,
            self.person_skipped,
            self.clinvar_records,
            self.clinvar_alleles,
            self.clinvar_other_sequences,
            self.clinvar_skipped,
            self.matched,
            jpares(&self.by_significance),
            jpares(&self.by_category),
            achados.join(","),
            self.plp_heterozygous,
            self.plp_homozygous,
            self.plp_hemizygous,
            self.plp_unknown_zygosity,
            jlista(&self.unknown_review_status),
            jnum_opt(self.written),
            json_texto(&self.filters),
            jopt(self.assembly.as_deref()),
            json_texto(&self.resumo())
        )
    }

    /// A short, factual text report.
    pub fn resumo(&self) -> String {
        let quem = self.sample.as_deref().unwrap_or("(sites-only VCF: every listed allele)");
        let mut o = format!("ClinVar matches for {quem} ({})\n", self.person);
        o.push_str(&format!(
            "  ClinVar: {} ({} records{})\n",
            self.clinvar,
            self.clinvar_records,
            self.assembly.as_deref().map(|a| format!(", {a}")).unwrap_or_default()
        ));
        o.push_str(&format!(
            "  person: {} records; {} carried alleles on 1-22, X, Y, MT checked; {} alleles not carried (0/0, no call, other ALT); {} records on other sequences and {} symbolic alleles not checked\n",
            self.person_records, self.person_alleles, self.person_not_carried, self.person_other_sequences, self.person_not_comparable
        ));
        if let Some(b) = &self.callable {
            o.push_str(&format!(
                "  only calls inside the callable regions of {b} were checked; {} carried alleles outside them were left out\n",
                self.person_outside_callable
            ));
        }
        let cat: Vec<String> = self.by_category.iter().map(|(n, q)| format!("{} {q}", n.replace('_', " "))).collect();
        o.push_str(&format!("  found in ClinVar: {} alleles ({})\n", self.matched, cat.join(", ")));
        o.push_str(&format!("  filters: {}\n", self.filters));
        o.push_str(&format!("  listed: {}\n", self.findings.len()));
        let plp = self.plp_heterozygous + self.plp_homozygous + self.plp_hemizygous + self.plp_unknown_zygosity;
        if plp > 0 {
            o.push_str(&format!(
                "  classified pathogenic or likely pathogenic: {} heterozygous (one copy), {} homozygous (two copies), {} hemizygous, {} without genotype\n",
                self.plp_heterozygous, self.plp_homozygous, self.plp_hemizygous, self.plp_unknown_zygosity
            ));
            o.push_str("  One copy of a variant linked to a recessive condition is what is called a carrier state; ClinVar does not record the mode of inheritance.\n");
        }
        const MOSTRA: usize = 40;
        for f in self.findings.iter().take(MOSTRA) {
            let genes = if f.genes.is_empty() { "-".to_string() } else { f.genes.join(",") };
            let cond = f.conditions.iter().take(2).cloned().collect::<Vec<_>>().join("; ");
            o.push_str(&format!(
                "    {}:{} {}>{}  {genes}  {} ({})  {}, {} star(s)  {}{}\n",
                f.chrom,
                f.pos,
                f.ref_allele,
                f.alt_allele,
                f.zygosity,
                if f.genotype.is_empty() { "-" } else { f.genotype.as_str() },
                f.significance,
                f.stars,
                cond,
                if f.conditions.len() > 2 { "; ..." } else { "" }
            ));
        }
        if self.findings.len() > MOSTRA {
            o.push_str(&format!("    ... and {} more\n", self.findings.len() - MOSTRA));
        }
        if !self.unknown_review_status.is_empty() {
            o.push_str(&format!(
                "  review status values not in the star table (counted as 0 stars): {}\n",
                self.unknown_review_status.join(", ")
            ));
        }
        o.push_str("  These are matches against a public database, for research and teaching; they are not a diagnosis.\n");
        o
    }
}

/// Joins a person's VCF with the ClinVar VCF and lists, for the alleles the
/// person carries, what ClinVar says about them.
///
/// Both files must be sorted by position with the chromosomes in the order
/// 1..22, X, Y (chrM/MT anywhere); `sort_genomic(..., chrom_order="natural")`
/// gives that. Names are compared canonically (chr1 == 1). Records on other
/// sequences (alt contigs, NT_/NW_ scaffolds: ClinVar has 51) are skipped and
/// counted. Variants match on (chromosome, position, REF, ALT) after
/// splitting and trimming; see the top of this file about left-alignment.
/// With `person.callable`, only alleles entirely inside those regions count.
pub fn clinvar_findings(person: &Person, clinvar: &str, opt: &FindingsOptions) -> Result<FindingsReport, String> {
    let classes = classes_pedidas(&opt.classes)?;
    if opt.min_stars > 4 {
        return Err(format!("min_stars is {}, but ClinVar review stars go from 0 to 4", opt.min_stars));
    }
    let genes: HashSet<String> =
        opt.genes.iter().map(|g| g.trim().to_ascii_lowercase()).filter(|g| !g.is_empty()).collect();
    let regioes = match &opt.regions {
        Some(r) => Some(Intervalos::le(r, None)?),
        None => None,
    };
    if let Some(o) = &opt.output {
        recusa_mesmo_arquivo(&person.vcf, o)?;
        recusa_mesmo_arquivo(clinvar, o)?;
        if let Some(r) = &opt.regions {
            recusa_mesmo_arquivo(r, o)?;
        }
    }

    let mut p = Fluxo::novo(&person.vcf, true, person.sample.clone())?;
    if let Some(b) = &person.callable {
        p.chamaveis = Some(Intervalos::le(b, None)?);
    }
    let mut c = Fluxo::novo(clinvar, false, None)?;
    p.le_cabecalho()?;
    c.le_cabecalho()?;
    confere_montagens(&person.vcf, p.montagem, clinvar, c.montagem)?;

    let mut rel = FindingsReport { person: person.vcf.clone(), clinvar: clinvar.to_string(), ..Default::default() };
    let mut por_sig: HashMap<String, usize> = HashMap::new();
    let mut por_cat: HashMap<&'static str, usize> = HashMap::new();
    let mut desconhecidos: BTreeSet<String> = BTreeSet::new();
    let mut achados: Vec<Finding> = Vec::new();
    let mut matched = 0usize;

    let mut casou = |x: &Alelo, y: &Alelo| -> Result<(), String> {
        let pc = campos(&x.reg.linha).ok_or("internal error: a person record could not be read again")?;
        let cc = campos(&y.reg.linha).ok_or("internal error: a ClinVar record could not be read again")?;
        let cv = le_clinvar(cc.id, cc.info);
        matched += 1;
        *por_sig.entry(cv.significance.clone()).or_default() += 1;
        *por_cat.entry(cv.category).or_default() += 1;
        if let Some(r) = &cv.review_raw {
            if review_stars(r).is_none() {
                desconhecidos.insert(r.clone());
            }
        }
        let passa = classes.as_ref().map_or(true, |cl| cv.classes.iter().any(|t| cl.contains(t)))
            && cv.stars >= opt.min_stars
            && (genes.is_empty() || cv.genes.iter().any(|g| genes.contains(&g.to_ascii_lowercase())))
            && regioes.as_ref().map_or(true, |r| {
                r.toca(&canonico(pc.chrom), x.pos - 1, x.pos - 1 + x.ref_.len() as u64)
            });
        if !passa {
            return Ok(());
        }
        let col = x.reg.amostra;
        achados.push(Finding {
            chrom: pc.chrom.to_string(),
            pos: x.pos,
            ref_allele: x.ref_.clone(),
            alt_allele: x.alt.clone(),
            genes: cv.genes.clone(),
            rsid: cv.rsids.first().cloned(),
            variation_id: cv.variation_id.clone(),
            allele_id: cv.allele_id.clone(),
            genotype: valor_formato(&pc, col, "GT").unwrap_or("").to_string(),
            zygosity: x.zig.as_str().to_string(),
            depth: valor_presente(&pc, col, "DP").and_then(|v| v.parse().ok()),
            genotype_quality: valor_presente(&pc, col, "GQ").and_then(|v| v.parse::<f64>().ok()).filter(|v| v.is_finite()),
            significance: cv.significance.clone(),
            classes: cv.classes.iter().map(|s| s.to_string()).collect(),
            category: cv.category.to_string(),
            review_status: cv.review_status.clone(),
            stars: cv.stars,
            conditions: cv.conditions.clone(),
            consequences: cv.consequences.clone(),
            variant_type: cv.variant_type.clone(),
            oncogenicity: cv.oncogenicity.clone(),
            oncogenicity_stars: cv.oncogenicity_stars,
            somatic_impact: cv.somatic_impact.clone(),
            somatic_impact_stars: cv.somatic_impact_stars,
            person_line: x.reg.numero,
            clinvar_line: y.reg.numero,
        });
        Ok(())
    };
    junta(&mut p, &mut c, &mut casou)?;

    achados.sort_by(|a, b| {
        compara_crom(&a.chrom, &b.chrom)
            .then(a.pos.cmp(&b.pos))
            .then_with(|| a.ref_allele.cmp(&b.ref_allele))
            .then_with(|| a.alt_allele.cmp(&b.alt_allele))
            .then_with(|| natural_cmp(&a.variation_id, &b.variation_id))
    });
    for f in &achados {
        if f.category == "pathogenic_or_likely_pathogenic" {
            match f.zygosity.as_str() {
                "het" => rel.plp_heterozygous += 1,
                "hom" => rel.plp_homozygous += 1,
                "hemizygous" => rel.plp_hemizygous += 1,
                _ => rel.plp_unknown_zygosity += 1,
            }
        }
    }
    let mut sig: Vec<(String, usize)> = por_sig.into_iter().collect();
    sig.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    rel.by_significance = sig;
    rel.by_category = CATEGORIES.iter().map(|c| (c.to_string(), por_cat.get(c).copied().unwrap_or(0))).collect();
    rel.unknown_review_status = desconhecidos.into_iter().collect();
    rel.matched = matched;
    rel.sample = p.amostra_nome.clone();
    rel.person_records = p.registros;
    rel.person_alleles = p.alelos;
    rel.person_not_carried = p.nao_portados;
    rel.person_outside_callable = p.fora;
    rel.callable = person.callable.clone();
    rel.person_other_sequences = p.outras;
    rel.person_not_comparable = p.nao_comparaveis;
    rel.person_skipped = p.ignoradas;
    rel.clinvar_records = c.registros;
    rel.clinvar_alleles = c.alelos;
    rel.clinvar_other_sequences = c.outras;
    rel.clinvar_skipped = c.ignoradas;
    rel.assembly = c.montagem.or(p.montagem).map(str::to_string);
    rel.filters = {
        let mut f = vec![match &classes {
            None => "all classes".to_string(),
            Some(s) => {
                let mut v: Vec<&str> = s.iter().copied().collect();
                v.sort_by_key(|c| CLASSES.iter().position(|x| x == c));
                format!("classes {}", v.join(", "))
            }
        }];
        f.push(format!("at least {} star(s)", opt.min_stars));
        if !genes.is_empty() {
            f.push(format!("genes {}", opt.genes.join(", ")));
        }
        if let Some(r) = &opt.regions {
            f.push(format!("regions of {r}"));
        }
        f.join("; ")
    };
    rel.findings = achados;

    if let Some(saida) = &opt.output {
        let mut w = Escritor::cria(saida)?;
        let mut guarda = Saida { caminho: saida.clone(), armada: true };
        let erro = |e: String| format!("Error writing {saida}: {e}");
        w.linha(CABECALHO_TSV).map_err(erro)?;
        for f in &rel.findings {
            w.linha(&f.tsv()).map_err(erro)?;
        }
        w.termina().map_err(erro)?;
        guarda.armada = false;
        rel.written = Some(rel.findings.len());
    }
    Ok(rel)
}

// =====================================================================
//  2) GENE REGIONS FROM AN ANNOTATION
// =====================================================================

/// One gene found in the annotation.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GeneRegion {
    /// The gene as asked.
    pub query: String,
    /// The name in the annotation (Name=, gene_name or gene=).
    pub gene: String,
    /// The sequence name written to the BED: the chromosome when the
    /// annotation's own `region` features (or a RefSeq NC_ accession) say
    /// which one, in the naming asked for; otherwise the sequence as written.
    pub chrom: String,
    /// The sequence as written in the annotation (NC_000017.11, chr17, 17).
    pub sequence: String,
    /// On a primary chromosome (not an alt locus, patch or unplaced scaffold).
    pub primary: bool,
    /// BED coordinates (0-based, end excluded), padding included.
    pub start: u64,
    pub end: u64,
    pub strand: String,
    pub gene_id: Option<String>,
    pub biotype: Option<String>,
    pub feature: String,
    pub line: usize,
}

#[derive(Clone, Debug, Default)]
pub struct GeneRegions {
    pub annotation: String,
    pub regions: Vec<GeneRegion>,
    /// Asked for and not found.
    pub not_found: Vec<String>,
    /// Genes found on more than one sequence: (gene, the sequences). Only the
    /// primary chromosomes are kept when there is one (two for genes of the
    /// pseudoautosomal regions, on X and Y).
    pub several_sequences: Vec<(String, Vec<String>)>,
    /// Not found by name, but listed as a synonym (gene_synonym=) of these
    /// genes.
    pub synonyms: Vec<(String, Vec<String>)>,
    pub padding: u64,
    /// Sequences whose name was mapped (NC_000017.11 -> 17).
    pub sequences_mapped: usize,
    pub written: Option<usize>,
}

impl GeneRegions {
    pub fn json(&self) -> String {
        let reg: Vec<String> = self
            .regions
            .iter()
            .map(|r| {
                format!(
                    "{{\"query\":{},\"gene\":{},\"chrom\":{},\"sequence\":{},\"primary\":{},\"start\":{},\"end\":{},\"strand\":{},\"gene_id\":{},\"biotype\":{},\"feature\":{},\"line\":{}}}",
                    json_texto(&r.query),
                    json_texto(&r.gene),
                    json_texto(&r.chrom),
                    json_texto(&r.sequence),
                    r.primary,
                    r.start,
                    r.end,
                    json_texto(&r.strand),
                    jopt(r.gene_id.as_deref()),
                    jopt(r.biotype.as_deref()),
                    json_texto(&r.feature),
                    r.line
                )
            })
            .collect();
        let lista = |v: &[(String, Vec<String>)]| -> String {
            format!("[{}]", v.iter().map(|(g, s)| format!("[{},{}]", json_texto(g), jlista(s))).collect::<Vec<_>>().join(","))
        };
        format!(
            "{{\"annotation\":{},\"regions\":[{}],\"not_found\":{},\"several_sequences\":{},\"synonyms\":{},\"padding\":{},\"sequences_mapped\":{},\"written\":{},\"text\":{}}}",
            json_texto(&self.annotation),
            reg.join(","),
            jlista(&self.not_found),
            lista(&self.several_sequences),
            lista(&self.synonyms),
            self.padding,
            self.sequences_mapped,
            jnum_opt(self.written),
            json_texto(&self.resumo())
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!("Gene regions from {} (padding {} bp)\n", self.annotation, self.padding);
        for r in &self.regions {
            o.push_str(&format!(
                "  {:<10} {}:{}-{} ({}, {} bp){}{}\n",
                r.gene,
                r.chrom,
                r.start + 1,
                r.end,
                r.strand,
                r.end - r.start,
                if r.chrom != r.sequence { format!(", {}", r.sequence) } else { String::new() },
                if r.primary { "" } else { ", not a primary chromosome" }
            ));
        }
        for (g, s) in &self.several_sequences {
            o.push_str(&format!("  {g} is on several sequences: {}\n", s.join(", ")));
        }
        if !self.not_found.is_empty() {
            o.push_str(&format!("  not found: {}\n", self.not_found.join(", ")));
        }
        for (g, s) in &self.synonyms {
            o.push_str(&format!("  {g} is listed as a synonym of {}\n", s.join(", ")));
        }
        o
    }
}

#[derive(Clone, Debug, Default)]
struct InfoSeq {
    genome: Option<String>,
    cromossomo: Option<String>,
    nome: Option<String>,
    comprimento: Option<u64>,
    tem_regiao: bool,
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
fn atributos(col9: &str) -> Vec<(String, String)> {
    let mut v = Vec::new();
    for parte in col9.split(';') {
        let p = parte.trim();
        if p.is_empty() {
            continue;
        }
        let (eq, sp) = (p.find('='), p.find(' '));
        match (eq, sp) {
            (Some(e), s) if s.map_or(true, |s| e < s) => {
                v.push((p[..e].trim().to_string(), decodifica(p[e + 1..].trim())));
            }
            (_, Some(s)) => v.push((p[..s].trim().to_string(), p[s + 1..].trim().trim_matches('"').to_string())),
            _ => {}
        }
    }
    v
}

fn attr<'a>(at: &'a [(String, String)], chaves: &[&str]) -> Option<&'a str> {
    chaves.iter().find_map(|k| at.iter().find(|(c, _)| c.eq_ignore_ascii_case(k)).map(|(_, v)| v.as_str()))
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

/// (name, primary chromosome?) of a sequence of the annotation.
fn nome_da_sequencia(seqid: &str, info: Option<&InfoSeq>) -> (String, bool) {
    if let Some(i) = info.filter(|i| i.tem_regiao) {
        match i.genome.as_deref() {
            Some("chromosome") => {
                if let Some(c) = i.cromossomo.as_deref().or(i.nome.as_deref()) {
                    return (c.to_string(), true);
                }
            }
            Some("mitochondrion") => return ("MT".to_string(), true),
            Some(_) => return (seqid.to_string(), false),
            None => {}
        }
    }
    if let Some(c) = cromossomo_refseq(seqid) {
        return (c, true);
    }
    (seqid.to_string(), convencao(seqid) != Convencao::Outra)
}

struct Acerto {
    gene: String,
    seqid: String,
    ini: u64,
    fim: u64,
    fita: String,
    gene_id: Option<String>,
    biotipo: Option<String>,
    tipo: String,
    linha: usize,
}

/// The coordinates of the genes `genes` in a GFF3 or GTF annotation
/// (`.gz` accepted), from the gene features (type `gene`, `pseudogene`,
/// `ncRNA_gene`...: any type ending in "gene"), matching Name=, gene= or
/// gene_name= without regard to case.
///
/// RefSeq names its sequences by accession (NC_000017.11): they are mapped to
/// chromosome names by the annotation's own `region` features (chromosome=,
/// genome=chromosome / mitochondrion), or by the accession itself. Alt loci,
/// patches and unplaced scaffolds keep their accession and are not primary.
///
/// `padding` widens each region on both sides (clipped at 0 and at the
/// sequence length when the annotation gives it). `naming`: None keeps the
/// names as mapped (1, X, MT for RefSeq), "chr" gives chr1/chrM, "no_chr"
/// gives 1/MT. The BED written to `output` (if any) is sorted, with the gene
/// name in column 4 and the strand in column 6.
pub fn gene_regions(
    gff: &str,
    genes: &[String],
    padding: u64,
    naming: Option<&str>,
    output: Option<&str>,
) -> Result<GeneRegions, String> {
    let estilo = match naming.map(str::trim) {
        None | Some("") => None,
        Some(s) => Some(Estilo::de_texto(s).ok_or_else(|| {
            format!("Unknown naming \"{s}\": use \"chr\" (chr1, chrM), \"no_chr\" (1, MT), or leave it empty to keep the annotation's names")
        })?),
    };
    let mut pedidos: Vec<String> = Vec::new();
    let mut indice: HashMap<String, usize> = HashMap::new();
    for g in genes {
        let t = g.trim();
        if t.is_empty() {
            continue;
        }
        let k = t.to_ascii_lowercase();
        if !indice.contains_key(&k) {
            indice.insert(k, pedidos.len());
            pedidos.push(t.to_string());
        }
    }
    if pedidos.is_empty() {
        return Err("No gene given: pass gene symbols such as [\"BRCA1\", \"BRCA2\"]".into());
    }
    if let Some(o) = output {
        recusa_mesmo_arquivo(gff, o)?;
    }

    let mut seqs: HashMap<String, InfoSeq> = HashMap::new();
    let mut acertos: Vec<Vec<Acerto>> = (0..pedidos.len()).map(|_| Vec::new()).collect();
    let mut sinonimos: Vec<Vec<String>> = vec![Vec::new(); pedidos.len()];
    let mut l = Leitor::abre(gff)?;
    let mut viu_gene = false;
    while l.avanca()? {
        let t = l.texto();
        if t.is_empty() {
            continue;
        }
        if let Some(r) = t.strip_prefix("##sequence-region") {
            let c: Vec<&str> = r.split_whitespace().collect();
            if c.len() >= 3 {
                if let Ok(fim) = c[2].parse::<u64>() {
                    seqs.entry(c[0].to_string()).or_default().comprimento = Some(fim);
                }
            }
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
            continue;
        }
        let tipo = c[2];
        if tipo == "region" {
            let at = atributos(c[8]);
            let src = attr(&at, &["gbkey"]) == Some("Src");
            let inteira = c[3] == "1" && (attr(&at, &["chromosome"]).is_some() || attr(&at, &["genome"]).is_some());
            if src || inteira {
                let e = seqs.entry(c[0].to_string()).or_default();
                if !e.tem_regiao {
                    e.tem_regiao = true;
                    e.genome = attr(&at, &["genome"]).map(str::to_string);
                    e.cromossomo = attr(&at, &["chromosome"]).map(str::to_string);
                    e.nome = attr(&at, &["Name"]).map(str::to_string);
                    if let Ok(fim) = c[4].parse::<u64>() {
                        e.comprimento = Some(fim);
                    }
                }
            }
            continue;
        }
        if !tipo.to_ascii_lowercase().ends_with("gene") {
            continue;
        }
        viu_gene = true;
        let at = atributos(c[8]);
        let mut casados: Vec<usize> = Vec::new();
        let mut sin: Vec<usize> = Vec::new();
        for (k, v) in &at {
            let kl = k.to_ascii_lowercase();
            if kl == "name" || kl == "gene" || kl == "gene_name" {
                if let Some(&i) = indice.get(&v.to_ascii_lowercase()) {
                    if !casados.contains(&i) {
                        casados.push(i);
                    }
                }
            } else if kl == "gene_synonym" {
                for s in v.split(',') {
                    if let Some(&i) = indice.get(&s.trim().to_ascii_lowercase()) {
                        sin.push(i);
                    }
                }
            }
        }
        if casados.is_empty() && sin.is_empty() {
            continue;
        }
        let oficial = attr(&at, &["Name", "gene_name", "gene"]).unwrap_or("").to_string();
        for i in sin {
            if !casados.contains(&i) && !sinonimos[i].contains(&oficial) {
                sinonimos[i].push(oficial.clone());
            }
        }
        if casados.is_empty() {
            continue;
        }
        let (Ok(ini), Ok(fim)) = (c[3].parse::<u64>(), c[4].parse::<u64>()) else {
            return Err(format!("{gff}, line {}: start and end of a gene feature must be whole numbers", l.linha));
        };
        if ini == 0 || fim < ini {
            return Err(format!("{gff}, line {}: a gene feature from {ini} to {fim} (GFF coordinates are 1-based, start <= end)", l.linha));
        }
        let gene_id = attr(&at, &["Dbxref"])
            .and_then(|d| d.split(',').find_map(|x| x.strip_prefix("GeneID:")))
            .or_else(|| attr(&at, &["gene_id", "ID"]))
            .map(str::to_string);
        let biotipo = attr(&at, &["gene_biotype", "gene_type", "biotype"]).map(str::to_string);
        for i in casados {
            acertos[i].push(Acerto {
                gene: oficial.clone(),
                seqid: c[0].to_string(),
                ini,
                fim,
                fita: c[6].to_string(),
                gene_id: gene_id.clone(),
                biotipo: biotipo.clone(),
                tipo: tipo.to_string(),
                linha: l.linha,
            });
        }
    }
    if !viu_gene && l.linha > 0 && acertos.iter().all(Vec::is_empty) {
        return Err(format!(
            "{gff} has no gene features (column 3 `gene`, `pseudogene`, `ncRNA_gene`...): is it a GFF3 or GTF gene annotation?"
        ));
    }

    let mut rel = GeneRegions { annotation: gff.to_string(), padding, ..Default::default() };
    let mut mapeadas: HashSet<String> = HashSet::new();
    for (i, q) in pedidos.iter().enumerate() {
        let lista = &acertos[i];
        if lista.is_empty() {
            rel.not_found.push(q.clone());
            if !sinonimos[i].is_empty() {
                rel.synonyms.push((q.clone(), sinonimos[i].clone()));
            }
            continue;
        }
        let nomes: Vec<(String, bool)> = lista.iter().map(|a| nome_da_sequencia(&a.seqid, seqs.get(&a.seqid))).collect();
        let mut distintas: Vec<String> = Vec::new();
        for (a, (n, _)) in lista.iter().zip(&nomes) {
            let d = if *n != a.seqid { format!("{n} ({})", a.seqid) } else { a.seqid.clone() };
            if !distintas.contains(&d) {
                distintas.push(d);
            }
        }
        if distintas.len() > 1 {
            rel.several_sequences.push((q.clone(), distintas));
        }
        let algum_primario = nomes.iter().any(|(_, p)| *p);
        for (a, (n, primario)) in lista.iter().zip(&nomes) {
            if algum_primario && !primario {
                continue;
            }
            if *n != a.seqid {
                mapeadas.insert(a.seqid.clone());
            }
            let chrom = match (estilo, primario) {
                (Some(e), true) => converte_nome(n, e).unwrap_or_else(|| n.clone()),
                _ => n.clone(),
            };
            let mut fim = a.fim + padding;
            if let Some(t) = seqs.get(&a.seqid).and_then(|s| s.comprimento) {
                fim = fim.min(t);
            }
            rel.regions.push(GeneRegion {
                query: q.clone(),
                gene: if a.gene.is_empty() { q.clone() } else { a.gene.clone() },
                chrom,
                sequence: a.seqid.clone(),
                primary: *primario,
                start: (a.ini - 1).saturating_sub(padding),
                end: fim,
                strand: a.fita.clone(),
                gene_id: a.gene_id.clone(),
                biotype: a.biotipo.clone(),
                feature: a.tipo.clone(),
                line: a.linha,
            });
        }
    }
    rel.sequences_mapped = mapeadas.len();
    rel.regions.sort_by(|a, b| {
        compara_crom(&a.chrom, &b.chrom).then(a.start.cmp(&b.start)).then(a.end.cmp(&b.end)).then_with(|| a.gene.cmp(&b.gene))
    });

    if let Some(saida) = output {
        let mut w = Escritor::cria_para(saida, Formato::Bed)?;
        let mut guarda = Saida { caminho: saida.to_string(), armada: true };
        let erro = |e: String| format!("Error writing {saida}: {e}");
        for r in &rel.regions {
            let fita = if r.strand == "+" || r.strand == "-" { r.strand.as_str() } else { "." };
            w.linha(&format!("{}\t{}\t{}\t{}\t0\t{fita}", r.chrom, r.start, r.end, r.gene)).map_err(erro)?;
        }
        w.termina().map_err(erro)?;
        guarda.armada = false;
        rel.written = Some(rel.regions.len());
    }
    Ok(rel)
}

// =====================================================================
//  3) VARIANT LOOKUP ACROSS PEOPLE (CASCADE TESTING)
// =====================================================================

/// One variant of the lookup.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LookupVariant {
    /// The text it came from (a coordinate, an rsID or a VariationID).
    pub query: String,
    pub chrom: String,
    /// Minimal form, 1-based.
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    /// ClinVar's data, when a ClinVar VCF was given and has the variant.
    pub in_clinvar: bool,
    pub variation_id: Option<String>,
    pub rsid: Option<String>,
    pub genes: Vec<String>,
    pub significance: Option<String>,
    pub stars: Option<u8>,
    pub conditions: Vec<String>,
}

/// The state of one person at one variant.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Call {
    /// Index into `Lookup::variants` and `Lookup::people`.
    pub variant: usize,
    pub person: usize,
    /// present: the person carries the allele;
    /// other_allele: at this site the person carries another allele (listed
    ///   in `other_alleles`) and not this one;
    /// no_call: the record of the site has no genotype for the person;
    /// reference: a record covering the site calls the person 0/0;
    /// absent_callable: no record, and the site is inside the person's
    ///   callable regions -- so reference;
    /// not_callable: no record, and the site is outside the callable regions
    ///   -- nothing can be said;
    /// absent: no record, and no callable regions were given -- reference
    ///   only if the site was callable.
    pub status: String,
    /// The GT of the record that decided the status ("" when none).
    pub genotype: String,
    /// het, hom, hemizygous, hom_ref, other_alt, no_call, unknown; "" when
    /// there is no record.
    pub zygosity: String,
    pub depth: Option<u64>,
    pub genotype_quality: Option<f64>,
    /// AD as written (reference first, then each ALT).
    pub allele_depths: Option<String>,
    /// AD of this ALT, when the record has it.
    pub alt_depth: Option<u64>,
    pub qual: Option<String>,
    pub filter: Option<String>,
    /// The sequence name in this person's file.
    pub chrom: Option<String>,
    /// Line of the record (None when it was read through the index).
    pub line: Option<usize>,
    /// Other alleles carried at the site, as "chrom:pos:ref:alt (GT)".
    pub other_alleles: Vec<String>,
    /// Is the variant entirely inside the person's callable regions? None
    /// when no callable BED was given. A call outside them (GIAB files have
    /// some) is reported, but is not part of the benchmark.
    pub in_callable: Option<bool>,
}

#[derive(Clone, Debug, Default)]
pub struct Lookup {
    pub variants: Vec<LookupVariant>,
    /// rsIDs and VariationIDs not found in ClinVar.
    pub unresolved: Vec<String>,
    /// One label per person: the sample, or the file name.
    pub people: Vec<String>,
    /// One call per (variant, person), variant-major.
    pub calls: Vec<Call>,
    /// Per person: was the file read through its tabix index?
    pub indexed: Vec<bool>,
}

impl Lookup {
    pub fn json(&self) -> String {
        let vs: Vec<String> = self
            .variants
            .iter()
            .map(|v| {
                format!(
                    "{{\"query\":{},\"chrom\":{},\"pos\":{},\"ref\":{},\"alt\":{},\"in_clinvar\":{},\"variation_id\":{},\"rsid\":{},\"genes\":{},\"significance\":{},\"stars\":{},\"conditions\":{}}}",
                    json_texto(&v.query),
                    json_texto(&v.chrom),
                    v.pos,
                    json_texto(&v.ref_allele),
                    json_texto(&v.alt_allele),
                    v.in_clinvar,
                    jopt(v.variation_id.as_deref()),
                    jopt(v.rsid.as_deref()),
                    jlista(&v.genes),
                    jopt(v.significance.as_deref()),
                    jnum_opt(v.stars),
                    jlista(&v.conditions)
                )
            })
            .collect();
        let cs: Vec<String> = self
            .calls
            .iter()
            .map(|c| {
                format!(
                    "{{\"variant\":{},\"person\":{},\"status\":{},\"genotype\":{},\"zygosity\":{},\"depth\":{},\"genotype_quality\":{},\"allele_depths\":{},\"alt_depth\":{},\"qual\":{},\"filter\":{},\"chrom\":{},\"line\":{},\"other_alleles\":{},\"in_callable\":{}}}",
                    c.variant,
                    c.person,
                    json_texto(&c.status),
                    json_texto(&c.genotype),
                    json_texto(&c.zygosity),
                    jnum_opt(c.depth),
                    c.genotype_quality.map(jnum).unwrap_or_else(|| "null".into()),
                    jopt(c.allele_depths.as_deref()),
                    jnum_opt(c.alt_depth),
                    jopt(c.qual.as_deref()),
                    jopt(c.filter.as_deref()),
                    jopt(c.chrom.as_deref()),
                    jnum_opt(c.line),
                    jlista(&c.other_alleles),
                    jnum_opt(c.in_callable)
                )
            })
            .collect();
        format!(
            "{{\"variants\":[{}],\"unresolved\":{},\"people\":{},\"calls\":[{}],\"indexed\":[{}],\"text\":{}}}",
            vs.join(","),
            jlista(&self.unresolved),
            jlista(&self.people),
            cs.join(","),
            self.indexed.iter().map(|b| b.to_string()).collect::<Vec<_>>().join(","),
            json_texto(&self.resumo())
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!("Variant lookup: {} variant(s) in {} person(s)\n", self.variants.len(), self.people.len());
        let largura = self.people.iter().map(|p| p.len()).max().unwrap_or(4).max(4);
        for (i, v) in self.variants.iter().enumerate() {
            let mut extra: Vec<String> = Vec::new();
            if !v.genes.is_empty() {
                extra.push(v.genes.join(","));
            }
            if let Some(s) = &v.significance {
                extra.push(format!("{s}, {} star(s)", v.stars.unwrap_or(0)));
            }
            if let Some(r) = &v.rsid {
                extra.push(r.clone());
            }
            if let Some(id) = &v.variation_id {
                extra.push(format!("VariationID {id}"));
            }
            o.push_str(&format!(
                "  {}:{}:{}:{}{}{}\n",
                v.chrom,
                v.pos,
                v.ref_allele,
                v.alt_allele,
                if v.query != format!("{}:{}:{}:{}", v.chrom, v.pos, v.ref_allele, v.alt_allele) { format!(" [{}]", v.query) } else { String::new() },
                if extra.is_empty() { String::new() } else { format!("  ({})", extra.join("; ")) }
            ));
            for c in self.calls.iter().filter(|c| c.variant == i) {
                let mut det: Vec<String> = Vec::new();
                if !c.genotype.is_empty() {
                    det.push(format!("GT {}", c.genotype));
                }
                if !c.zygosity.is_empty() && c.status == "present" {
                    det.push(c.zygosity.clone());
                }
                if let Some(d) = c.depth {
                    det.push(format!("DP {d}"));
                }
                if let Some(a) = &c.allele_depths {
                    det.push(format!("AD {a}"));
                }
                if !c.other_alleles.is_empty() {
                    det.push(format!("other: {}", c.other_alleles.join(", ")));
                }
                o.push_str(&format!(
                    "      {:<largura$}  {:<15} {}\n",
                    self.people[c.person],
                    c.status,
                    det.join("  ")
                ));
            }
        }
        if !self.unresolved.is_empty() {
            o.push_str(&format!("  not found in ClinVar: {}\n", self.unresolved.join(", ")));
        }
        o.push_str("  absent_callable: no call, inside the callable regions (so reference); not_callable: no call, outside them; absent: no call and no callable regions given.\n");
        o
    }
}

enum Pedido {
    Coord { chrom: String, pos: u64, ref_: String, alt: String },
    Rs(String),
    Vid(String),
}

const FORMAS_VARIANTE: &str = "Use chrom:pos:ref:alt (1-based, e.g. 17:43045712:T:C or chr17-43045712-T-C), an rsID (rs80357906) or a ClinVar VariationID (17661 or VCV000017661)";

fn le_pedido(t: &str) -> Result<Pedido, String> {
    let s = t.trim();
    let digitos = |x: &str| !x.is_empty() && x.bytes().all(|b| b.is_ascii_digit());
    if let Some(r) = s.strip_prefix("rs").or_else(|| s.strip_prefix("RS")).or_else(|| s.strip_prefix("Rs")) {
        if digitos(r) {
            return Ok(Pedido::Rs(r.to_string()));
        }
    }
    let vcv = s.strip_prefix("VCV").or_else(|| s.strip_prefix("vcv")).map(|r| r.split('.').next().unwrap_or(r));
    if let Some(r) = vcv.or(if digitos(s) { Some(s) } else { None }) {
        if digitos(r) {
            let n = r.trim_start_matches('0');
            return Ok(Pedido::Vid(if n.is_empty() { "0".into() } else { n.to_string() }));
        }
    }
    let sep = if s.contains(':') { ':' } else { '-' };
    let partes: Vec<&str> = s.rsplitn(4, sep).collect();
    if partes.len() == 4 {
        let (alt, ref_, pos, chrom) = (partes[0].trim(), partes[1].trim(), partes[2].trim(), partes[3].trim());
        let pos_ok: Option<u64> = pos.replace(',', "").parse().ok().filter(|p| *p >= 1);
        if let Some(p) = pos_ok {
            if !chrom.is_empty() && so_bases(ref_) && so_bases(alt) {
                return Ok(Pedido::Coord {
                    chrom: chrom.to_string(),
                    pos: p,
                    ref_: ref_.to_ascii_uppercase(),
                    alt: alt.to_ascii_uppercase(),
                });
            }
        }
    }
    Err(format!("Variant \"{t}\" is not understood. {FORMAS_VARIANTE}."))
}

/// Reads the records of a VCF and hands each line (header lines too) to
/// `recebe(line, line number)`. With `regioes` (canonical name -> 1-based
/// inclusive spans) and a `.tbi` next to the file, only the records touching
/// the regions are read, through the index, and the line number is None.
/// Returns whether the index was used.
fn varre_vcf(
    caminho: &str,
    regioes: Option<&HashMap<String, Vec<(u64, u64)>>>,
    recebe: &mut dyn FnMut(&str, Option<usize>) -> Result<(), String>,
) -> Result<bool, String> {
    if let Some(reg) = regioes {
        let tbi = format!("{caminho}.tbi");
        if Path::new(&tbi).exists() {
            let (_, _, nomes) = indice::le_tbi(&tbi)?;
            let mut textos: Vec<String> = Vec::new();
            for nome in &nomes {
                if let Some(v) = reg.get(&canonico(nome)) {
                    for &(a, b) in v {
                        textos.push(format!("{nome}:{a}-{b}"));
                    }
                }
            }
            if textos.is_empty() {
                // None of the sequences asked for is in the file: only the header.
                let mut l = Leitor::abre(caminho)?;
                while l.avanca()? {
                    let t = l.texto();
                    if !t.starts_with('#') {
                        break;
                    }
                    recebe(t, Some(l.linha))?;
                }
            } else {
                indice::consultar(caminho, &textos, &mut |_, linha| recebe(linha, None))?;
            }
            return Ok(true);
        }
    }
    let mut l = Leitor::abre(caminho)?;
    while l.avanca()? {
        let n = l.linha;
        recebe(l.texto(), Some(n))?;
    }
    Ok(false)
}

/// Above this many variants a person's file is read whole instead of through
/// its index.
const LIMITE_INDICE: usize = 10_000;

fn anota(v: &mut LookupVariant, cv: &ClinVar) {
    v.in_clinvar = true;
    v.variation_id = Some(cv.variation_id.clone());
    v.rsid = cv.rsids.first().cloned();
    v.genes = cv.genes.clone();
    v.significance = Some(cv.significance.clone());
    v.stars = Some(cv.stars);
    v.conditions = cv.conditions.clone();
}

/// Resolves rsIDs and VariationIDs through ClinVar and annotates the
/// coordinate queries with what ClinVar says.
fn resolve_clinvar(cv: &str, pedidos: &[Pedido], queries: &[String], por_pedido: &mut [Vec<LookupVariant>]) -> Result<(), String> {
    let mut vids: HashMap<String, Vec<usize>> = HashMap::new();
    let mut rs: HashMap<String, Vec<usize>> = HashMap::new();
    let mut coords: HashMap<String, HashMap<u64, Vec<usize>>> = HashMap::new();
    let mut regioes: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    for (i, p) in pedidos.iter().enumerate() {
        match p {
            Pedido::Rs(n) => rs.entry(n.clone()).or_default().push(i),
            Pedido::Vid(n) => vids.entry(n.clone()).or_default().push(i),
            Pedido::Coord { chrom, .. } => {
                let v = &por_pedido[i][0];
                let canon = canonico(chrom);
                coords.entry(canon.clone()).or_default().entry(v.pos).or_default().push(i);
                regioes.entry(canon).or_default().push((v.pos, v.pos + v.ref_allele.len() as u64 - 1));
            }
        }
    }
    let so_coordenadas = vids.is_empty() && rs.is_empty();
    let mut viu_clnsig = false;
    let mut viu_registro = false;
    let mut cache = (String::new(), String::new());
    varre_vcf(cv, if so_coordenadas { Some(&regioes) } else { None }, &mut |linha, _| {
        if linha.starts_with('#') {
            if linha.starts_with("##INFO=<ID=CLNSIG,") {
                viu_clnsig = true;
            }
            return Ok(());
        }
        let Some(c) = campos(linha) else { return Ok(()) };
        if !viu_registro {
            viu_registro = true;
            if !viu_clnsig {
                return Err(format!(
                    "{cv} does not look like the ClinVar VCF: its header does not declare INFO/CLNSIG. Download clinvar.vcf.gz from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/ (or vcf_GRCh37/)."
                ));
            }
        }
        let mut alvos: Vec<usize> = Vec::new();
        if let Some(v) = vids.get(c.id) {
            alvos.extend(v);
        }
        if !rs.is_empty() {
            if let Some(r) = info_valor(c.info, "RS") {
                for x in r.split(',') {
                    if let Some(v) = rs.get(x) {
                        alvos.extend(v);
                    }
                }
            }
        }
        if c.chrom != cache.0 {
            cache = (c.chrom.to_string(), canonico(c.chrom));
        }
        let por_pos = coords.get(&cache.1);
        if alvos.is_empty() && por_pos.is_none() {
            return Ok(());
        }
        let info = le_clinvar(c.id, c.info);
        for alt in c.alt.split(',') {
            let Some((np, nr, na)) = minimo(c.pos, c.ref_, alt) else { continue };
            for &i in &alvos {
                if por_pedido[i].iter().any(|v| v.pos == np && v.ref_allele == nr && v.alt_allele == na && canonico(&v.chrom) == cache.1) {
                    continue;
                }
                let mut v = LookupVariant {
                    query: queries[i].clone(),
                    chrom: c.chrom.to_string(),
                    pos: np,
                    ref_allele: nr.clone(),
                    alt_allele: na.clone(),
                    ..Default::default()
                };
                anota(&mut v, &info);
                por_pedido[i].push(v);
            }
            if let Some(is) = por_pos.and_then(|pp| pp.get(&np)) {
                for &i in is {
                    let v = &mut por_pedido[i][0];
                    if v.ref_allele == nr && v.alt_allele == na {
                        anota(v, &info);
                    }
                }
            }
        }
        Ok(())
    })?;
    if !viu_registro && !viu_clnsig {
        return Err(format!("{cv} does not look like the ClinVar VCF: no INFO/CLNSIG declaration and no record was read."));
    }
    Ok(())
}

#[derive(Clone)]
struct Estado {
    prioridade: u8,
    exato: bool,
    call: Call,
}

/// The genotype facts of one call: the ALT numbers it holds, whether every
/// allele is the reference, whether something is missing.
fn le_gt(gt: Option<&str>) -> (Vec<u32>, bool, bool) {
    let Some(gt) = gt else { return (Vec::new(), false, false) };
    let (mut alts, mut todos_ref, mut falta) = (Vec::new(), true, false);
    for a in gt.split(['/', '|']) {
        match a.trim().parse::<u32>() {
            Ok(0) => {}
            Ok(x) => {
                todos_ref = false;
                if !alts.contains(&x) {
                    alts.push(x);
                }
            }
            Err(_) => {
                todos_ref = false;
                falta = true;
            }
        }
    }
    (alts, todos_ref, falta)
}

/// Looks up each variant in each person.
///
/// `queries`: "chrom:pos:ref:alt" (1-based; `-` also separates), rsIDs or
/// ClinVar VariationIDs; the last two need `clinvar`, through which they are
/// resolved (an rsID may give several alleles). With `clinvar`, coordinate
/// queries are also annotated with what ClinVar says.
///
/// `people`: several files and/or several samples of one file (a file shared
/// by several people is read once). A file with a `.tbi` index is read
/// through it when there are at most 10,000 variants; otherwise it is read
/// whole, so any VCF works, sorted or not.
pub fn variant_lookup(queries: &[String], people: &[Person], clinvar: Option<&str>) -> Result<Lookup, String> {
    if queries.iter().all(|q| q.trim().is_empty()) {
        return Err(format!("No variant given. {FORMAS_VARIANTE}."));
    }
    if people.is_empty() {
        return Err("No person given: pass at least one VCF".into());
    }
    let queries: Vec<String> = queries.iter().filter(|q| !q.trim().is_empty()).map(|q| q.trim().to_string()).collect();
    let pedidos: Vec<Pedido> = queries.iter().map(|q| le_pedido(q)).collect::<Result<_, _>>()?;
    if clinvar.is_none() {
        if let Some(i) = pedidos.iter().position(|p| !matches!(p, Pedido::Coord { .. })) {
            return Err(format!(
                "\"{}\" is an rsID or a ClinVar VariationID; these are resolved through ClinVar: pass clinvar=\"clinvar.vcf.gz\", or give the variant as chrom:pos:ref:alt",
                queries[i]
            ));
        }
    }
    let mut por_pedido: Vec<Vec<LookupVariant>> = vec![Vec::new(); pedidos.len()];
    for (i, p) in pedidos.iter().enumerate() {
        if let Pedido::Coord { chrom, pos, ref_, alt } = p {
            let (np, nr, na) = minimo(*pos, ref_, alt)
                .ok_or_else(|| format!("Variant \"{}\": REF and ALT are the same, so there is nothing to look up", queries[i]))?;
            por_pedido[i].push(LookupVariant {
                query: queries[i].clone(),
                chrom: chrom.clone(),
                pos: np,
                ref_allele: nr,
                alt_allele: na,
                ..Default::default()
            });
        }
    }
    let mut rel = Lookup::default();
    if let Some(cv) = clinvar {
        resolve_clinvar(cv, &pedidos, &queries, &mut por_pedido)?;
        for (i, p) in pedidos.iter().enumerate() {
            if !matches!(p, Pedido::Coord { .. }) && por_pedido[i].is_empty() {
                rel.unresolved.push(queries[i].clone());
            }
        }
    }
    let variantes: Vec<LookupVariant> = por_pedido.into_iter().flatten().collect();

    // The variants by sequence, sorted by start, for the overlap test.
    let mut por_crom: HashMap<String, Vec<(u64, u64, usize)>> = HashMap::new();
    let mut maior: HashMap<String, u64> = HashMap::new();
    for (i, v) in variantes.iter().enumerate() {
        let canon = canonico(&v.chrom);
        let fim = v.pos + v.ref_allele.len() as u64;
        por_crom.entry(canon.clone()).or_default().push((v.pos, fim, i));
        let m = maior.entry(canon).or_insert(0);
        *m = (*m).max(fim - v.pos);
    }
    for l in por_crom.values_mut() {
        l.sort_unstable();
    }
    let regioes: HashMap<String, Vec<(u64, u64)>> =
        por_crom.iter().map(|(c, l)| (c.clone(), l.iter().map(|&(a, b, _)| (a, b - 1)).collect())).collect();

    let mut estados: Vec<Vec<Option<Estado>>> = vec![vec![None; people.len()]; variantes.len()];
    rel.people = vec![String::new(); people.len()];
    rel.indexed = vec![false; people.len()];

    // One pass per distinct file.
    let mut grupos: Vec<(String, Vec<usize>)> = Vec::new();
    for (i, p) in people.iter().enumerate() {
        match grupos.iter_mut().find(|g| g.0 == p.vcf) {
            Some(g) => g.1.push(i),
            None => grupos.push((p.vcf.clone(), vec![i])),
        }
    }
    for (caminho, membros) in &grupos {
        let mut cols: Vec<Option<usize>> = vec![None; membros.len()];
        let mut decidido = false;
        let mut cache = (String::new(), String::new());
        let usar_indice = variantes.len() <= LIMITE_INDICE;
        let mut rotulos: Vec<String> = vec![String::new(); membros.len()];
        let indexado = varre_vcf(caminho, if usar_indice { Some(&regioes) } else { None }, &mut |linha, numero| {
            if linha.starts_with("#CHROM") {
                let nomes: Vec<String> = linha.split('\t').skip(9).map(str::to_string).collect();
                for (mi, &pi) in membros.iter().enumerate() {
                    cols[mi] = escolhe_amostra(caminho, &nomes, people[pi].sample.as_deref())?;
                    rotulos[mi] = match cols[mi] {
                        Some(j) => nomes[j].clone(),
                        None => Path::new(caminho).file_name().map(|f| f.to_string_lossy().to_string()).unwrap_or_else(|| caminho.clone()),
                    };
                }
                decidido = true;
                return Ok(());
            }
            if linha.starts_with('#') || linha.is_empty() {
                return Ok(());
            }
            let Some(c) = campos(linha) else { return Ok(()) };
            if !decidido {
                for (mi, &pi) in membros.iter().enumerate() {
                    if let Some(s) = &people[pi].sample {
                        return Err(format!("{caminho} has no #CHROM header line, so sample \"{s}\" cannot be found. Add the VCF header."));
                    }
                    cols[mi] = Some(0);
                    rotulos[mi] = Path::new(caminho).file_name().map(|f| f.to_string_lossy().to_string()).unwrap_or_default();
                }
                decidido = true;
            }
            if c.chrom != cache.0 {
                cache = (c.chrom.to_string(), canonico(c.chrom));
            }
            let Some(lista) = por_crom.get(&cache.1) else { return Ok(()) };
            let l = maior[&cache.1];
            let s0 = c.pos;
            let mut s1 = c.pos + (c.ref_.len() as u64).max(1);
            if let Some(e) = info_valor(c.info, "END").and_then(|e| e.parse::<u64>().ok()) {
                s1 = s1.max(e + 1);
            }
            let a = lista.partition_point(|x| x.0 + l <= s0);
            let b = lista.partition_point(|x| x.0 < s1);
            if a >= b {
                return Ok(());
            }
            let alelos: Vec<Option<(u64, String, String)>> = c.alt.split(',').map(|alt| minimo(c.pos, c.ref_, alt)).collect();
            let texto_alelo = |k: usize, gt: &str| -> String {
                match &alelos[k] {
                    Some((p, r, a)) => format!("{}:{p}:{r}:{a} ({gt})", c.chrom),
                    None => format!("{}:{}:{}:{} ({gt})", c.chrom, c.pos, c.ref_, c.alt.split(',').nth(k).unwrap_or(".")),
                }
            };
            for &(ini, fim, vi) in &lista[a..b] {
                if !(ini < s1 && s0 < fim) {
                    continue;
                }
                let v = &variantes[vi];
                let exato = alelos
                    .iter()
                    .position(|x| matches!(x, Some((p, r, aa)) if *p == v.pos && *r == v.ref_allele && *aa == v.alt_allele))
                    .map(|i| i as u32 + 1);
                for (mi, &pi) in membros.iter().enumerate() {
                    let col = cols[mi];
                    let gt = valor_formato(&c, col, "GT");
                    let sem_genotipo = col.is_none() || gt.is_none();
                    let gt_txt = gt.unwrap_or("");
                    let (alts, todos_ref, falta) = le_gt(gt);
                    let carrega = |k: u32| sem_genotipo || alts.contains(&k);
                    let outros: Vec<String> = (0..alelos.len())
                        .filter(|&i| Some(i as u32 + 1) != exato && carrega(i as u32 + 1) && c.alt.split(',').nth(i) != Some("."))
                        .map(|i| texto_alelo(i, if gt_txt.is_empty() { "listed" } else { gt_txt }))
                        .collect();
                    let (status, zig, prioridade) = match exato {
                        Some(k) => {
                            let z = if sem_genotipo { Zygosity::Unknown } else { zigosidade(gt, k) };
                            match z {
                                Zygosity::Het | Zygosity::Hom | Zygosity::Hemizygous | Zygosity::Unknown => ("present", z, 5u8),
                                Zygosity::OtherAlt => ("other_allele", z, 4),
                                Zygosity::NoCall => ("no_call", z, 3),
                                Zygosity::HomRef => ("reference", z, 2),
                            }
                        }
                        None => {
                            if !outros.is_empty() {
                                ("other_allele", Zygosity::OtherAlt, 4)
                            } else if sem_genotipo {
                                continue;
                            } else if falta && alts.is_empty() {
                                ("no_call", Zygosity::NoCall, 3)
                            } else if todos_ref {
                                ("reference", Zygosity::HomRef, 2)
                            } else {
                                continue;
                            }
                        }
                    };
                    let ad = valor_presente(&c, col, "AD");
                    let novo = Call {
                        variant: vi,
                        person: pi,
                        status: status.to_string(),
                        genotype: gt_txt.to_string(),
                        zygosity: zig.as_str().to_string(),
                        depth: valor_presente(&c, col, "DP").and_then(|x| x.parse().ok()),
                        genotype_quality: valor_presente(&c, col, "GQ").and_then(|x| x.parse::<f64>().ok()).filter(|x| x.is_finite()),
                        allele_depths: ad.map(str::to_string),
                        alt_depth: match (exato, ad) {
                            (Some(k), Some(ad)) => ad.split(',').nth(k as usize).and_then(|x| x.parse().ok()),
                            _ => None,
                        },
                        qual: Some(c.qual.to_string()),
                        filter: Some(c.filter.to_string()),
                        chrom: Some(c.chrom.to_string()),
                        line: numero,
                        other_alleles: Vec::new(),
                        in_callable: None,
                    };
                    let e = &mut estados[vi][pi];
                    let mut acumulados = e.as_ref().map(|x| x.call.other_alleles.clone()).unwrap_or_default();
                    for o in outros {
                        if !acumulados.contains(&o) {
                            acumulados.push(o);
                        }
                    }
                    let troca = match e {
                        None => true,
                        Some(x) => prioridade > x.prioridade || (prioridade == x.prioridade && exato.is_some() && !x.exato),
                    };
                    if troca {
                        *e = Some(Estado { prioridade, exato: exato.is_some(), call: novo });
                    }
                    if let Some(x) = e {
                        x.call.other_alleles = acumulados;
                    }
                }
            }
            Ok(())
        })?;
        if !decidido {
            for (mi, &pi) in membros.iter().enumerate() {
                if let Some(s) = &people[pi].sample {
                    return Err(format!("{caminho} has no #CHROM header line, so sample \"{s}\" cannot be found. Add the VCF header."));
                }
                rotulos[mi] = Path::new(caminho).file_name().map(|f| f.to_string_lossy().to_string()).unwrap_or_default();
            }
        }
        for (mi, &pi) in membros.iter().enumerate() {
            rel.people[pi] = rotulos[mi].clone();
            rel.indexed[pi] = indexado;
        }
    }

    // No record: the callable regions decide.
    let cromossomos: HashSet<String> = por_crom.keys().cloned().collect();
    let mut bed: Vec<Option<Intervalos>> = Vec::with_capacity(people.len());
    for p in people {
        bed.push(match &p.callable {
            Some(b) => Some(Intervalos::le(b, Some(&cromossomos))?),
            None => None,
        });
    }
    for (vi, v) in variantes.iter().enumerate() {
        let canon = canonico(&v.chrom);
        for pi in 0..people.len() {
            let dentro = bed[pi].as_ref().map(|b| b.contem(&canon, v.pos - 1, v.pos - 1 + v.ref_allele.len() as u64));
            let mut call = match estados[vi][pi].take() {
                Some(e) => e.call,
                None => {
                    let status = match dentro {
                        None => "absent",
                        Some(true) => "absent_callable",
                        Some(false) => "not_callable",
                    };
                    Call { variant: vi, person: pi, status: status.to_string(), ..Default::default() }
                }
            };
            call.in_callable = dentro;
            rel.calls.push(call);
        }
    }
    rel.variants = variantes;
    Ok(rel)
}

// =====================================================================
//  4) PILEUP AT ONE POSITION
// =====================================================================

#[derive(Clone, Debug)]
pub struct PileupOptions {
    /// Reads with a lower mapping quality are left out.
    pub min_mapq: u8,
    /// Bases with a lower base quality are left out (13, as samtools).
    pub min_baseq: u8,
    /// Reads with any of these FLAG bits are left out. Default 0xF04:
    /// unmapped (0x4), secondary (0x100), QC fail (0x200), duplicate (0x400),
    /// supplementary (0x800).
    pub exclude_flags: u16,
    /// An allele to report the fraction of: a base ("T"), "REF>ALT" in VCF
    /// form anchored at this position ("C>T", "CAG>C", "C>CTT"), "+SEQ" (an
    /// insertion after this position), "-N" or "-SEQ" (a deletion of N bases
    /// after it) or "*" (reads with a deletion over this position).
    pub alt: Option<String>,
}

impl Default for PileupOptions {
    fn default() -> Self {
        PileupOptions { min_mapq: 0, min_baseq: 13, exclude_flags: 0xF04, alt: None }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AlleleCount {
    /// A, C, G, T, N, `*` (deletion over the position), `+SEQ` (insertion
    /// after it), `-N` (deletion of N bases after it).
    pub allele: String,
    pub count: usize,
    pub forward: usize,
    pub reverse: usize,
    /// count / depth.
    pub fraction: f64,
}

#[derive(Clone, Debug, Default)]
pub struct Pileup {
    pub path: String,
    /// The sequence name in the file.
    pub chrom: String,
    pub pos: u64,
    /// Read through the `.bai` index.
    pub indexed: bool,
    /// Alignments that span the position, before the filters.
    pub reads: usize,
    pub filtered_flag: usize,
    pub filtered_mapq: usize,
    pub filtered_baseq: usize,
    /// Reads that skip the position (N in the CIGAR: an intron), after the
    /// FLAG and mapping-quality filters. Not part of the depth.
    pub ref_skips: usize,
    /// Reads with a base or a deletion at the position, after the filters.
    pub depth: usize,
    pub a: usize,
    pub c: usize,
    pub g: usize,
    pub t: usize,
    pub n: usize,
    pub deletions: usize,
    /// Reads of `depth` with an insertion right after the position.
    pub insertions_after: usize,
    /// Reads of `depth` with a deletion starting right after the position.
    pub deletions_after: usize,
    pub alleles: Vec<AlleleCount>,
    pub alt: Option<String>,
    pub alt_count: Option<usize>,
    pub alt_fraction: Option<f64>,
    pub ref_allele: Option<String>,
    pub ref_count: Option<usize>,
    pub ref_fraction: Option<f64>,
    pub options: String,
}

impl Pileup {
    pub fn json(&self) -> String {
        let al: Vec<String> = self
            .alleles
            .iter()
            .map(|a| {
                format!(
                    "{{\"allele\":{},\"count\":{},\"forward\":{},\"reverse\":{},\"fraction\":{}}}",
                    json_texto(&a.allele),
                    a.count,
                    a.forward,
                    a.reverse,
                    jnum(a.fraction)
                )
            })
            .collect();
        format!(
            "{{\"path\":{},\"chrom\":{},\"pos\":{},\"indexed\":{},\"reads\":{},\"filtered_flag\":{},\"filtered_mapq\":{},\"filtered_baseq\":{},\"ref_skips\":{},\"depth\":{},\"a\":{},\"c\":{},\"g\":{},\"t\":{},\"n\":{},\"deletions\":{},\"insertions_after\":{},\"deletions_after\":{},\"alleles\":[{}],\"alt\":{},\"alt_count\":{},\"alt_fraction\":{},\"ref\":{},\"ref_count\":{},\"ref_fraction\":{},\"options\":{},\"text\":{}}}",
            json_texto(&self.path),
            json_texto(&self.chrom),
            self.pos,
            self.indexed,
            self.reads,
            self.filtered_flag,
            self.filtered_mapq,
            self.filtered_baseq,
            self.ref_skips,
            self.depth,
            self.a,
            self.c,
            self.g,
            self.t,
            self.n,
            self.deletions,
            self.insertions_after,
            self.deletions_after,
            al.join(","),
            jopt(self.alt.as_deref()),
            jnum_opt(self.alt_count),
            self.alt_fraction.map(jnum).unwrap_or_else(|| "null".into()),
            jopt(self.ref_allele.as_deref()),
            jnum_opt(self.ref_count),
            self.ref_fraction.map(jnum).unwrap_or_else(|| "null".into()),
            json_texto(&self.options),
            json_texto(&self.resumo())
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!("Pileup of {} at {}:{} ({})\n", self.path, self.chrom, self.pos, self.options);
        o.push_str(&format!(
            "  {} read(s) span the position; left out: {} by FLAG, {} by mapping quality, {} by base quality; {} skip it (N in CIGAR)\n",
            self.reads, self.filtered_flag, self.filtered_mapq, self.filtered_baseq, self.ref_skips
        ));
        o.push_str(&format!("  depth {}:", self.depth));
        for a in &self.alleles {
            if a.count > 0 || matches!(a.allele.as_str(), "A" | "C" | "G" | "T") {
                o.push_str(&format!("  {} {} ({}; +{} -{})", a.allele, a.count, percent(a.count, self.depth), a.forward, a.reverse));
            }
        }
        o.push('\n');
        if let (Some(r), Some(n)) = (&self.ref_allele, self.ref_count) {
            o.push_str(&format!("  REF {r}: {n} of {} ({})\n", self.depth, percent(n, self.depth)));
        }
        if let (Some(a), Some(n)) = (&self.alt, self.alt_count) {
            o.push_str(&format!("  ALT {a}: {n} of {} ({})\n", self.depth, percent(n, self.depth)));
        }
        o
    }
}

/// What the alt option asks for.
#[derive(Clone, Debug, PartialEq)]
enum Alvo {
    Base(u8),
    DelSobre,
    Ins(String),
    DelDepois(u32),
}

fn le_alvo(s: &str) -> Result<(Alvo, Option<String>, String), String> {
    let t = s.trim().to_ascii_uppercase();
    let erro = || {
        format!(
            "Allele \"{s}\" is not understood. Use a base (T), REF>ALT anchored at the position (C>T, CAG>C, C>CTT), +SEQ (insertion after the position), -N or -SEQ (deletion after it) or * (deletion over it)"
        )
    };
    let base = |b: u8| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N');
    if t.len() == 1 && base(t.as_bytes()[0]) {
        return Ok((Alvo::Base(t.as_bytes()[0]), None, t));
    }
    if t == "*" {
        return Ok((Alvo::DelSobre, None, t));
    }
    if let Some(r) = t.strip_prefix('+') {
        if so_bases(r) {
            return Ok((Alvo::Ins(r.to_string()), None, t));
        }
        return Err(erro());
    }
    if let Some(r) = t.strip_prefix('-') {
        if let Ok(n) = r.parse::<u32>() {
            if n > 0 {
                return Ok((Alvo::DelDepois(n), None, t));
            }
        }
        if so_bases(r) {
            return Ok((Alvo::DelDepois(r.len() as u32), None, format!("-{}", r.len())));
        }
        return Err(erro());
    }
    if let Some((r, a)) = t.split_once('>') {
        let (r, a) = (r.trim(), a.trim());
        if !so_bases(r) || !so_bases(a) || r == a {
            return Err(erro());
        }
        let ref1 = r[..1].to_string();
        if r.len() == 1 && a.len() == 1 {
            return Ok((Alvo::Base(a.as_bytes()[0]), Some(ref1), a.to_string()));
        }
        if r.len() == 1 && a.len() > 1 && a.as_bytes()[0] == r.as_bytes()[0] {
            return Ok((Alvo::Ins(a[1..].to_string()), Some(ref1), format!("+{}", &a[1..])));
        }
        if a.len() == 1 && r.len() > 1 && a.as_bytes()[0] == r.as_bytes()[0] {
            return Ok((Alvo::DelDepois(r.len() as u32 - 1), Some(ref1), format!("-{}", r.len() - 1)));
        }
        return Err(format!(
            "Allele \"{s}\": a pileup counts one position, so only SNVs and indels anchored at it (first base shared, as VCF writes them) can be counted; split an MNP or a complex change into its parts"
        ));
    }
    Err(erro())
}

fn le_cigar(c: &str) -> Option<Vec<(u32, u8)>> {
    let mut v = Vec::new();
    let mut n: u64 = 0;
    let mut tem = false;
    for b in c.bytes() {
        if b.is_ascii_digit() {
            n = n * 10 + u64::from(b - b'0');
            if n > u64::from(u32::MAX >> 4) {
                return None;
            }
            tem = true;
        } else {
            if !tem || !b"MIDNSHP=X".contains(&b) {
                return None;
            }
            v.push((n as u32, b));
            n = 0;
            tem = false;
        }
    }
    if tem || v.is_empty() {
        return None;
    }
    Some(v)
}

enum Obs {
    Base(u8, u8),
    Del(u8),
    Pulo,
}

struct Observacao {
    obs: Obs,
    ins: Option<String>,
    del_depois: Option<u32>,
}

/// The insertion or deletion that follows the operation `i`, read at its
/// last reference base, as htslib's pileup reports it: the next operation is
/// D (and the current one is not) -- the run of D operations; the next is I
/// -- that insertion; the next is P -- the I operations up to the next
/// reference-consuming one. `q`: query index of the base after the current one.
fn indel_depois(cigar: &[(u32, u8)], i: usize, op: u8, seq: &[u8], q: usize) -> (Option<String>, Option<u32>) {
    let Some(&(l2, op2)) = cigar.get(i + 1) else { return (None, None) };
    let pedaco = |a: usize, n: usize| -> String {
        (a..a + n).map(|j| seq.get(j).map(|b| b.to_ascii_uppercase() as char).unwrap_or('N')).collect()
    };
    match op2 {
        b'D' if op != b'D' => {
            let mut total = l2;
            for &(l, o) in &cigar[i + 2..] {
                if o == b'D' {
                    total += l;
                } else {
                    break;
                }
            }
            (None, Some(total))
        }
        b'I' => (Some(pedaco(q, l2 as usize)), None),
        b'P' => {
            let mut s = String::new();
            let mut qq = q;
            for &(l, o) in &cigar[i + 2..] {
                match o {
                    b'I' => {
                        s.push_str(&pedaco(qq, l as usize));
                        qq += l as usize;
                    }
                    b'D' | b'M' | b'N' | b'=' | b'X' => break,
                    _ => {}
                }
            }
            if s.is_empty() {
                (None, None)
            } else {
                (Some(s), None)
            }
        }
        _ => (None, None),
    }
}

/// What one alignment shows at reference position `alvo` (1-based). `None`:
/// its alignment does not span the position.
fn observa(pos1: u64, cigar: &[(u32, u8)], seq: &[u8], qual: Option<&[u8]>, alvo: u64) -> Option<Observacao> {
    let qualidade = |i: usize| -> u8 {
        if seq.is_empty() {
            return 0;
        }
        match qual {
            None => 255,
            Some(q) => q.get(i).map(|c| c.saturating_sub(33)).unwrap_or(0),
        }
    };
    let mut r = pos1;
    let mut q = 0usize;
    for (i, &(len, op)) in cigar.iter().enumerate() {
        let l = u64::from(len);
        match op {
            b'M' | b'=' | b'X' => {
                if alvo >= r && alvo < r + l {
                    let qi = q + (alvo - r) as usize;
                    let b = seq.get(qi).map(|b| b.to_ascii_uppercase()).unwrap_or(b'N');
                    let (ins, del_depois) = if alvo == r + l - 1 { indel_depois(cigar, i, op, seq, qi + 1) } else { (None, None) };
                    return Some(Observacao { obs: Obs::Base(b, qualidade(qi)), ins, del_depois });
                }
                r += l;
                q += len as usize;
            }
            b'D' => {
                if alvo >= r && alvo < r + l {
                    let (ins, del_depois) = if alvo == r + l - 1 { indel_depois(cigar, i, op, seq, q) } else { (None, None) };
                    return Some(Observacao { obs: Obs::Del(qualidade(q)), ins, del_depois });
                }
                r += l;
            }
            b'N' => {
                if alvo >= r && alvo < r + l {
                    return Some(Observacao { obs: Obs::Pulo, ins: None, del_depois: None });
                }
                r += l;
            }
            b'I' | b'S' => q += len as usize,
            _ => {}
        }
        if r > alvo {
            return None;
        }
    }
    None
}

#[derive(Default)]
struct Contagem {
    reads: usize,
    flag: usize,
    mapq: usize,
    baseq: usize,
    pulos: usize,
    /// A C G T N * : total and forward strand.
    bases: [(usize, usize); 6],
    ins: BTreeMap<String, (usize, usize)>,
    dels: BTreeMap<u32, (usize, usize)>,
    ins_total: usize,
    dels_total: usize,
}

impl Contagem {
    fn leitura(&mut self, linha: &str, alvo: u64, opt: &PileupOptions, caminho: &str) -> Result<(), String> {
        let c: Vec<&str> = linha.splitn(12, '\t').collect();
        if c.len() < 11 {
            return Ok(());
        }
        let (Ok(flag), Ok(pos1)) = (c[1].parse::<u16>(), c[3].parse::<u64>()) else { return Ok(()) };
        if c[5] == "*" || pos1 == 0 {
            return Ok(());
        }
        let cigar = le_cigar(c[5]).ok_or_else(|| format!("{caminho}: read {} has an invalid CIGAR \"{}\"", c[0], c[5]))?;
        let seq: &[u8] = if c[9] == "*" { &[] } else { c[9].as_bytes() };
        let qual: Option<&[u8]> = if c[10] == "*" { None } else { Some(c[10].as_bytes()) };
        let Some(o) = observa(pos1, &cigar, seq, qual, alvo) else { return Ok(()) };
        self.reads += 1;
        if flag & opt.exclude_flags != 0 {
            self.flag += 1;
            return Ok(());
        }
        let mapq: u32 = c[4].parse().unwrap_or(0);
        if mapq < u32::from(opt.min_mapq) {
            self.mapq += 1;
            return Ok(());
        }
        let (i, q) = match o.obs {
            Obs::Pulo => {
                self.pulos += 1;
                return Ok(());
            }
            Obs::Base(b, q) => (
                match b {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => 4,
                },
                q,
            ),
            Obs::Del(q) => (5, q),
        };
        if q < opt.min_baseq {
            self.baseq += 1;
            return Ok(());
        }
        let frente = usize::from(flag & 0x10 == 0);
        self.bases[i].0 += 1;
        self.bases[i].1 += frente;
        if let Some(s) = o.ins {
            let e = self.ins.entry(s).or_default();
            e.0 += 1;
            e.1 += frente;
            self.ins_total += 1;
        }
        if let Some(d) = o.del_depois {
            let e = self.dels.entry(d).or_default();
            e.0 += 1;
            e.1 += frente;
            self.dels_total += 1;
        }
        Ok(())
    }
}

fn acha_bai(caminho: &str) -> Option<String> {
    let ao_lado = format!("{caminho}.bai");
    if Path::new(&ao_lado).exists() {
        return Some(ao_lado);
    }
    let trocado = format!("{}.bai", &caminho[..caminho.len() - 4]);
    if Path::new(&trocado).exists() {
        return Some(trocado);
    }
    None
}

fn lista_nomes(nomes: &[String]) -> String {
    let mut v: Vec<&str> = nomes.iter().take(8).map(String::as_str).collect();
    if nomes.len() > 8 {
        v.push("...");
    }
    v.join(", ")
}

/// The name, among `nomes`, of the sequence `pedido` (exactly, or as the
/// same chromosome under the other naming).
fn resolve_nome(pedido: &str, nomes: &[String], caminho: &str) -> Result<String, String> {
    if let Some(n) = nomes.iter().find(|n| *n == pedido) {
        return Ok(n.clone());
    }
    let c = canonico(pedido);
    nomes.iter().find(|n| canonico(n) == c).cloned().ok_or_else(|| {
        format!("Sequence \"{pedido}\" is not in {caminho}. Its sequences are named like: {}", lista_nomes(nomes))
    })
}

/// Counts what the reads show at `chrom:pos` (1-based): each base, deletions
/// over the position, and insertions and deletions right after it, with the
/// FLAG, mapping-quality and base-quality filters of `opt`.
///
/// A coordinate-sorted BAM with its `.bai` is read through the index; a BAM
/// without index, a SAM or a `.sam.gz` is read whole. The rules follow
/// `samtools mpileup` (htslib's pileup): a deletion is filtered by the quality
/// of the base after it; a read with QUAL `*` passes any base-quality filter;
/// reads skipping the position (N) are counted apart and are not in the
/// depth. Comparable output: `samtools mpileup -A -B -x -d 0 -Q <min_baseq>
/// -q <min_mapq> --ff <exclude_flags>` (count the A/C/G/T/N and `*` of the
/// read-bases column, not `<` and `>`). Overlapping mates are both counted.
pub fn pileup(path: &str, chrom: &str, pos: u64, opt: &PileupOptions) -> Result<Pileup, String> {
    if pos == 0 {
        return Err("Positions are 1-based: the first base of a sequence is 1".into());
    }
    let alvo = match &opt.alt {
        Some(a) if !a.trim().is_empty() => Some(le_alvo(a)?),
        _ => None,
    };
    let mut cont = Contagem::default();
    let e_bam = path.to_ascii_lowercase().ends_with(".bam");
    let mut nome = chrom.to_string();
    let mut indexado = false;
    if e_bam && acha_bai(path).is_some() {
        let arq = std::fs::File::open(path).map_err(|e| format!("Error opening {path}: {e}"))?;
        let mut r = BgzfReader::novo(BufReader::new(arq));
        let cab = bam::le_cabecalho(&mut r).map_err(|e| format!("{path}: {e}"))?;
        let nomes: Vec<String> = cab.refs.iter().map(|x| x.0.clone()).collect();
        nome = resolve_nome(chrom, &nomes, path)?;
        if pos > indice::POS_MAX as u64 {
            return Err(format!("Position {pos} is beyond what a BAM index addresses (512 Mbp)"));
        }
        let regiao = format!("{nome}:{pos}-{pos}");
        indice::consultar(path, &[regiao], &mut |cab, l| if cab { Ok(()) } else { cont.leitura(l, pos, opt, path) })?;
        indexado = true;
    } else {
        let mut l = Leitor::abre(path)?;
        let mut sq: Vec<String> = Vec::new();
        let mut alvo_nome: Option<String> = None;
        let canon = canonico(chrom);
        let mut cache = (String::new(), false);
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
            if alvo_nome.is_none() && !sq.is_empty() {
                alvo_nome = Some(resolve_nome(chrom, &sq, path)?);
            }
            let rname = t.split('\t').nth(2).unwrap_or("");
            if rname != cache.0 {
                let casa = match &alvo_nome {
                    Some(n) => rname == n,
                    None => rname == chrom || canonico(rname) == canon,
                };
                cache = (rname.to_string(), casa);
                if casa {
                    nome = rname.to_string();
                }
            }
            if cache.1 {
                cont.leitura(t, pos, opt, path)?;
            }
        }
        match alvo_nome {
            Some(n) => nome = n,
            // A header with @SQ and no record: the name must still exist.
            None if !sq.is_empty() => nome = resolve_nome(chrom, &sq, path)?,
            None => {}
        }
    }

    let mut p = Pileup {
        path: path.to_string(),
        chrom: nome,
        pos,
        indexed: indexado,
        reads: cont.reads,
        filtered_flag: cont.flag,
        filtered_mapq: cont.mapq,
        filtered_baseq: cont.baseq,
        ref_skips: cont.pulos,
        ..Default::default()
    };
    p.depth = cont.bases.iter().map(|b| b.0).sum();
    [p.a, p.c, p.g, p.t, p.n, p.deletions] = cont.bases.map(|b| b.0);
    p.insertions_after = cont.ins_total;
    p.deletions_after = cont.dels_total;
    let profundidade = p.depth;
    let fr = move |n: usize| if profundidade == 0 { 0.0 } else { n as f64 / profundidade as f64 };
    let mut alelos: Vec<AlleleCount> = ["A", "C", "G", "T", "N", "*"]
        .iter()
        .zip(cont.bases.iter())
        .map(|(n, &(t, f))| AlleleCount { allele: n.to_string(), count: t, forward: f, reverse: t - f, fraction: fr(t) })
        .collect();
    let mut ins: Vec<(&String, &(usize, usize))> = cont.ins.iter().collect();
    ins.sort_by(|a, b| b.1 .0.cmp(&a.1 .0).then_with(|| a.0.cmp(b.0)));
    for (s, &(t, f)) in ins {
        alelos.push(AlleleCount { allele: format!("+{s}"), count: t, forward: f, reverse: t - f, fraction: fr(t) });
    }
    let mut dels: Vec<(&u32, &(usize, usize))> = cont.dels.iter().collect();
    dels.sort_by(|a, b| b.1 .0.cmp(&a.1 .0).then_with(|| a.0.cmp(b.0)));
    for (d, &(t, f)) in dels {
        alelos.push(AlleleCount { allele: format!("-{d}"), count: t, forward: f, reverse: t - f, fraction: fr(t) });
    }
    if let Some((a, r, rotulo)) = alvo {
        let n = match &a {
            Alvo::Base(b) => cont.bases[match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            }]
            .0,
            Alvo::DelSobre => cont.bases[5].0,
            Alvo::Ins(s) => cont.ins.get(s).map(|x| x.0).unwrap_or(0),
            Alvo::DelDepois(d) => cont.dels.get(d).map(|x| x.0).unwrap_or(0),
        };
        p.alt = Some(rotulo);
        p.alt_count = Some(n);
        p.alt_fraction = Some(fr(n));
        if let Some(r) = r {
            let i = match r.as_bytes()[0] {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            };
            p.ref_count = Some(cont.bases[i].0);
            p.ref_fraction = Some(fr(cont.bases[i].0));
            p.ref_allele = Some(r);
        }
    }
    p.alleles = alelos;
    p.options = format!(
        "min mapping quality {}, min base quality {}, FLAG filter 0x{:X}",
        opt.min_mapq, opt.min_baseq, opt.exclude_flags
    );
    Ok(p)
}
