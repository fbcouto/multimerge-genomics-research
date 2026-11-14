// =====================================================================
//  VCF TOOLS: STATISTICS, FILTER AND NORMALIZATION
//
//  The three things people do to a VCF before anything else:
//
//   - `vcf_stats`      what is in the file: how many variants of each type,
//                      Ts/Tv, FILTER and QUAL, records per sequence;
//   - `vcf_filter`     keep the records that satisfy an expression, fall in
//                      a region, and/or keep only some samples;
//   - `vcf_normalize`  split multiallelic sites and trim every variant to its
//                      minimal representation.
//
//  House rules:
//
//   - ANY SIZE. Every function is a streaming pass that holds one line at a
//     time. Normalization can move POS forward and break the order of the
//     file; when that happens the output is rebuilt by the out-of-core stable
//     sorter (`externo.rs`), never by loading the file.
//   - Input is plain text, `.gz` (BGZF included) or `.lz4`; the extension of
//     the output chooses its compression. CRLF is accepted, LF is written.
//   - Header lines are copied. Each operation adds one `##MMGEN_...` line
//     just before `#CHROM` recording what was done.
//   - A line with fewer than 8 columns or a non-numeric POS is not a record:
//     it is left out and COUNTED. Blank lines are dropped silently.
//   - Every message says what is wrong, where, and how to fix it.
// =====================================================================

use crate::externo::{abre_leitura, ordenar_arquivos, Criterio, DirTemp, Escritor, Formato, Opcoes};
use crate::nomes::json_texto;
use crate::referencia::Genoma;
use std::collections::{HashMap, VecDeque};
use std::io::BufRead;

// ---------------------------------------------------------------------
//  READING AND WRITING
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

/// Removes a half-written output when the operation ends in error -- but only
/// an output this operation created: a file that already existed is never
/// deleted.
struct Saida<'a> {
    caminho: &'a str,
    armada: bool,
}

impl Drop for Saida<'_> {
    fn drop(&mut self) {
        if self.armada {
            let _ = std::fs::remove_file(self.caminho);
        }
    }
}

fn mesmo_arquivo(a: &str, b: &str) -> bool {
    a == b || matches!((std::fs::canonicalize(a), std::fs::canonicalize(b)), (Ok(x), Ok(y)) if x == y)
}

fn recusa_mesmo_arquivo(entrada: &str, saida: &str) -> Result<(), String> {
    if mesmo_arquivo(entrada, saida) {
        return Err(format!(
            "The output `{saida}` is the input file itself. Writing would destroy the input while it is being read: choose a different output path."
        ));
    }
    Ok(())
}

fn json_num(x: f64) -> String {
    if x.is_finite() {
        format!("{x}")
    } else {
        "null".to_string()
    }
}

// ---------------------------------------------------------------------
//  RECORDS AND ALLELES
// ---------------------------------------------------------------------

/// The eight fixed columns of a record, borrowed from the line.
struct Registro<'a> {
    chrom: &'a str,
    pos_txt: &'a str,
    pos: u64,
    id: &'a str,
    ref_: &'a str,
    alt: &'a str,
    qual: &'a str,
    filter: &'a str,
    info: &'a str,
    /// FORMAT and the sample columns, still joined by tabs.
    resto: Option<&'a str>,
}

/// `None`: the line is not a record (fewer than 8 columns, or POS is not a
/// non-negative integer).
fn registro(t: &str) -> Option<Registro<'_>> {
    let mut c = t.splitn(9, '\t');
    let chrom = c.next()?;
    let pos_txt = c.next()?;
    let pos = pos_txt.parse().ok()?;
    Some(Registro {
        chrom,
        pos_txt,
        pos,
        id: c.next()?,
        ref_: c.next()?,
        alt: c.next()?,
        qual: c.next()?,
        filter: c.next()?,
        info: c.next()?,
        resto: c.next(),
    })
}

/// A finite number, or nothing. (`parse::<f64>` alone accepts "inf" and "nan".)
fn numero(s: &str) -> Option<f64> {
    s.parse::<f64>().ok().filter(|x| x.is_finite())
}

/// Non-empty and made only of A, C, G, T, N in either case.
fn so_bases(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}

fn sem_alt(alt: &str) -> bool {
    alt == "." || alt.is_empty()
}

fn n_alt(alt: &str) -> usize {
    if sem_alt(alt) {
        0
    } else {
        alt.bytes().filter(|b| *b == b',').count() + 1
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TipoAlelo {
    Snp,
    Mnp,
    Ins,
    Del,
    Outro,
}

/// The type of one ALT allele against REF, by LENGTH (case ignored):
/// SNP = both one base and different; MNP = same length above one;
/// insertion = ALT longer; deletion = ALT shorter.
///
/// `Outro` is everything that is not a plain sequence change: symbolic
/// (`<DEL>`), breakend (`G]17:198982]`), `*`, IUPAC codes, a REF that is not
/// made of ACGTN, and an ALT identical to REF (which is not a variant).
fn tipo_alelo(ref_: &str, alt: &str) -> TipoAlelo {
    if !so_bases(ref_) || !so_bases(alt) {
        TipoAlelo::Outro
    } else if alt.len() > ref_.len() {
        TipoAlelo::Ins
    } else if alt.len() < ref_.len() {
        TipoAlelo::Del
    } else if alt.eq_ignore_ascii_case(ref_) {
        TipoAlelo::Outro
    } else if alt.len() == 1 {
        TipoAlelo::Snp
    } else {
        TipoAlelo::Mnp
    }
}

/// `Some(true)`: transition (A<->G, C<->T). `Some(false)`: transversion.
/// `None`: N is involved, and the change is neither.
fn transicao(r: u8, a: u8) -> Option<bool> {
    match (r.to_ascii_uppercase(), a.to_ascii_uppercase()) {
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C') => Some(true),
        (b'N', _) | (_, b'N') => None,
        _ => Some(false),
    }
}

// =====================================================================
//  A) STATISTICS
// =====================================================================

/// What `vcf_stats` found. Type counts are per SITE: a site counts once for
/// each type present among its ALT alleles, so a multiallelic site with a SNP
/// and an insertion is in both `snps` and `insertions`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct VcfStats {
    pub records: usize,
    /// Columns after FORMAT in the `#CHROM` line.
    pub samples: usize,
    /// (sequence name, records), in order of first appearance.
    pub sequences: Vec<(String, usize)>,
    pub snps: usize,
    pub mnps: usize,
    pub insertions: usize,
    pub deletions: usize,
    /// Symbolic, breakend, `*`, non-ACGTN alleles, ALT equal to REF.
    pub other: usize,
    /// ALT is `.`.
    pub no_alt: usize,
    /// Sites with more than one ALT allele.
    pub multiallelic: usize,
    /// Per SNP allele. A SNP to or from `N` is in neither count.
    pub transitions: usize,
    pub transversions: usize,
    /// transitions / transversions; 0.0 when there are no transversions.
    pub ts_tv: f64,
    /// FILTER is `PASS`.
    pub pass: usize,
    /// FILTER is `.`.
    pub filter_missing: usize,
    /// QUAL is `.`.
    pub qual_missing: usize,
    /// Mean over the records with a numeric QUAL; 0.0 when there are none.
    pub mean_qual: f64,
    /// Lines with fewer than 8 columns or a non-numeric POS, left out.
    pub skipped: usize,
}

impl VcfStats {
    pub fn json(&self) -> String {
        let seqs: Vec<String> = self.sequences.iter().map(|(n, q)| format!("[{},{q}]", json_texto(n))).collect();
        format!(
            "{{\"records\":{},\"samples\":{},\"sequences\":[{}],\"snps\":{},\"mnps\":{},\"insertions\":{},\"deletions\":{},\"other\":{},\"no_alt\":{},\"multiallelic\":{},\"transitions\":{},\"transversions\":{},\"ts_tv\":{},\"pass\":{},\"filter_missing\":{},\"qual_missing\":{},\"mean_qual\":{},\"skipped\":{}}}",
            self.records,
            self.samples,
            seqs.join(","),
            self.snps,
            self.mnps,
            self.insertions,
            self.deletions,
            self.other,
            self.no_alt,
            self.multiallelic,
            self.transitions,
            self.transversions,
            json_num(self.ts_tv),
            self.pass,
            self.filter_missing,
            self.qual_missing,
            json_num(self.mean_qual),
            self.skipped
        )
    }

    /// The statistics as a short text report.
    pub fn resumo(&self) -> String {
        let mut o = format!(
            "{} record(s), {} sample(s), {} sequence(s)\n",
            self.records,
            self.samples,
            self.sequences.len()
        );
        o.push_str(&format!(
            "  sites by type: {} SNP, {} MNP, {} insertion, {} deletion, {} other, {} without ALT; {} multiallelic\n",
            self.snps, self.mnps, self.insertions, self.deletions, self.other, self.no_alt, self.multiallelic
        ));
        let razao = if self.transversions > 0 { format!("{:.3}", self.ts_tv) } else { "n/a".to_string() };
        o.push_str(&format!(
            "  SNP alleles: {} transition(s), {} transversion(s), Ts/Tv {razao}\n",
            self.transitions, self.transversions
        ));
        o.push_str(&format!(
            "  FILTER: {} PASS, {} missing; QUAL: mean {:.2}, {} missing\n",
            self.pass, self.filter_missing, self.mean_qual, self.qual_missing
        ));
        if !self.sequences.is_empty() {
            const MOSTRA: usize = 12;
            let lista: Vec<String> = self.sequences.iter().take(MOSTRA).map(|(n, q)| format!("{n} {q}")).collect();
            let mais = self.sequences.len().saturating_sub(MOSTRA);
            let resto = if mais > 0 { format!(", and {mais} more") } else { String::new() };
            o.push_str(&format!("  records per sequence: {}{resto}\n", lista.join(", ")));
        }
        if self.skipped > 0 {
            o.push_str(&format!(
                "  {} line(s) left out: fewer than 8 tab-separated columns or a POS that is not a number (validate_file(...) shows the lines)\n",
                self.skipped
            ));
        }
        o
    }
}

/// Counts what is in a VCF, in one streaming pass.
///
/// Allele types are decided by length, as documented in `VcfStats`; unlike
/// bcftools, alleles are not trimmed first, so `AC -> AT` is an MNP here.
pub fn vcf_stats(caminho: &str) -> Result<VcfStats, String> {
    let mut l = Leitor::abre(caminho)?;
    let mut s = VcfStats::default();
    let mut indice: HashMap<String, usize> = HashMap::new();
    // Index, in `s.sequences`, of the sequence of the previous record: inside
    // a chromosome the name is not looked up again.
    let mut atual = usize::MAX;
    let (mut soma_qual, mut com_qual) = (0f64, 0usize);

    while l.avanca()? {
        let t = l.texto();
        if t.is_empty() {
            continue;
        }
        if t.starts_with('#') {
            if t.starts_with("#CHROM") {
                s.samples = t.split('\t').count().saturating_sub(9);
            }
            continue;
        }
        let Some(r) = registro(t) else {
            s.skipped += 1;
            continue;
        };
        s.records += 1;

        if atual == usize::MAX || s.sequences[atual].0 != r.chrom {
            atual = match indice.get(r.chrom) {
                Some(&i) => i,
                None => {
                    s.sequences.push((r.chrom.to_string(), 0));
                    indice.insert(r.chrom.to_string(), s.sequences.len() - 1);
                    s.sequences.len() - 1
                }
            };
        }
        s.sequences[atual].1 += 1;

        if sem_alt(r.alt) {
            s.no_alt += 1;
        } else {
            let (mut snp, mut mnp, mut ins, mut del, mut outro) = (false, false, false, false, false);
            let mut alelos = 0usize;
            for a in r.alt.split(',') {
                alelos += 1;
                match tipo_alelo(r.ref_, a) {
                    TipoAlelo::Snp => {
                        snp = true;
                        match transicao(r.ref_.as_bytes()[0], a.as_bytes()[0]) {
                            Some(true) => s.transitions += 1,
                            Some(false) => s.transversions += 1,
                            None => {}
                        }
                    }
                    TipoAlelo::Mnp => mnp = true,
                    TipoAlelo::Ins => ins = true,
                    TipoAlelo::Del => del = true,
                    TipoAlelo::Outro => outro = true,
                }
            }
            s.snps += snp as usize;
            s.mnps += mnp as usize;
            s.insertions += ins as usize;
            s.deletions += del as usize;
            s.other += outro as usize;
            if alelos > 1 {
                s.multiallelic += 1;
            }
        }

        match r.filter {
            "PASS" => s.pass += 1,
            "." => s.filter_missing += 1,
            _ => {}
        }
        if r.qual == "." {
            s.qual_missing += 1;
        } else if let Some(q) = numero(r.qual) {
            soma_qual += q;
            com_qual += 1;
        }
    }

    if s.transversions > 0 {
        s.ts_tv = s.transitions as f64 / s.transversions as f64;
    }
    if com_qual > 0 {
        s.mean_qual = soma_qual / com_qual as f64;
    }
    Ok(s)
}

// =====================================================================
//  B) FILTER
// =====================================================================

// ---------------------------------------------------------------------
//  The expression language: tokens
// ---------------------------------------------------------------------

const AJUDA_EXPRESSAO: &str = "Fields: QUAL, POS, N_ALT (numbers); CHROM, ID, REF, ALT, FILTER, TYPE (text, in double quotes); INFO/KEY (number or text). Operators: == != < <= > >=, ~ and !~ (TYPE and FILTER), && || ! and parentheses. As in bcftools, TYPE==\"snp\" needs every ALT allele to be a SNP and TYPE~\"snp\" at least one; FILTER~\"q10\" means q10 is among the filters. Example: QUAL>=30 && TYPE==\"snp\" && INFO/DP>10";

const VALORES_TYPE: &str = "\"snp\", \"mnp\", \"ins\", \"del\", \"indel\", \"other\" or \"ref\"";

/// Longest expression accepted, in tokens, and deepest nesting of `(` and `!`.
/// The limits keep a pathological expression from exhausting the stack.
const MAX_TOKENS: usize = 2000;
const MAX_NIVEL: usize = 100;

/// A syntax error: what was expected, the position (in characters, from 1),
/// the expression with a mark under the position, and the language in a line.
fn erro_expressao(expr: &str, pos: usize, msg: &str) -> String {
    let visivel: String = expr.chars().map(|c| if c.is_control() { ' ' } else { c }).collect();
    format!(
        "Invalid filter expression at position {pos}: {msg}\n  {visivel}\n  {}^\n{AJUDA_EXPRESSAO}",
        " ".repeat(pos.saturating_sub(1))
    )
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Cmp {
    Ig,
    Dif,
    Men,
    MenIg,
    Mai,
    MaiIg,
    /// `~`: TYPE de pelo menos um alelo; FILTER entre os do registro.
    Tem,
    /// `!~`: o contrario de `~`.
    NaoTem,
}

impl Cmp {
    fn simbolo(self) -> &'static str {
        match self {
            Cmp::Ig => "==",
            Cmp::Dif => "!=",
            Cmp::Men => "<",
            Cmp::MenIg => "<=",
            Cmp::Mai => ">",
            Cmp::MaiIg => ">=",
            Cmp::Tem => "~",
            Cmp::NaoTem => "!~",
        }
    }

    /// The operator that says the same with the operands swapped: `30 < QUAL`
    /// is `QUAL > 30`.
    fn invertido(self) -> Cmp {
        match self {
            Cmp::Men => Cmp::Mai,
            Cmp::MenIg => Cmp::MaiIg,
            Cmp::Mai => Cmp::Men,
            Cmp::MaiIg => Cmp::MenIg,
            igual => igual,
        }
    }

    fn vale(self, a: f64, b: f64) -> bool {
        match self {
            Cmp::Ig => a == b,
            Cmp::Dif => a != b,
            Cmp::Men => a < b,
            Cmp::MenIg => a <= b,
            Cmp::Mai => a > b,
            Cmp::MaiIg => a >= b,
            // Nunca chega aqui: `monta` recusa `~` com numeros.
            Cmp::Tem | Cmp::NaoTem => false,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
enum Campo {
    Qual,
    Pos,
    NAlt,
    Chrom,
    Id,
    Ref,
    Alt,
    Filter,
    Tipo,
    Info(String),
}

impl Campo {
    fn de_nome(nome: &str) -> Option<Campo> {
        Some(match nome {
            "QUAL" => Campo::Qual,
            "POS" => Campo::Pos,
            "N_ALT" => Campo::NAlt,
            "CHROM" => Campo::Chrom,
            "ID" => Campo::Id,
            "REF" => Campo::Ref,
            "ALT" => Campo::Alt,
            "FILTER" => Campo::Filter,
            "TYPE" => Campo::Tipo,
            _ => return None,
        })
    }

    fn nome(&self) -> String {
        match self {
            Campo::Qual => "QUAL".to_string(),
            Campo::Pos => "POS".to_string(),
            Campo::NAlt => "N_ALT".to_string(),
            Campo::Chrom => "CHROM".to_string(),
            Campo::Id => "ID".to_string(),
            Campo::Ref => "REF".to_string(),
            Campo::Alt => "ALT".to_string(),
            Campo::Filter => "FILTER".to_string(),
            Campo::Tipo => "TYPE".to_string(),
            Campo::Info(k) => format!("INFO/{k}"),
        }
    }

    fn numerico(&self) -> bool {
        matches!(self, Campo::Qual | Campo::Pos | Campo::NAlt)
    }
}

#[derive(Clone, Debug, PartialEq)]
enum Tok {
    Campo(Campo),
    Num(f64),
    Txt(String),
    Cmp(Cmp),
    E,
    Ou,
    Nao,
    Abre,
    Fecha,
    Fim,
}

#[derive(Clone, Debug)]
struct Token {
    tok: Tok,
    /// Position of the first character, from 1.
    pos: usize,
    /// The token as written, for messages.
    fonte: String,
}

fn descreve(t: &Token) -> String {
    if t.tok == Tok::Fim {
        "the end of the expression".to_string()
    } else {
        format!("`{}`", t.fonte)
    }
}

fn tokeniza(expr: &str) -> Result<Vec<Token>, String> {
    let c: Vec<char> = expr.chars().collect();
    let n = c.len();
    let erro = |pos: usize, msg: &str| erro_expressao(expr, pos, msg);
    let mut v: Vec<Token> = Vec::new();
    let mut i = 0usize;

    while i < n {
        if c[i].is_whitespace() {
            i += 1;
            continue;
        }
        let ini = i;
        let prox = c.get(i + 1).copied();
        let tok = match c[i] {
            '(' => {
                i += 1;
                Tok::Abre
            }
            ')' => {
                i += 1;
                Tok::Fecha
            }
            '&' => {
                i += if prox == Some('&') { 2 } else { 1 };
                Tok::E
            }
            '|' => {
                i += if prox == Some('|') { 2 } else { 1 };
                Tok::Ou
            }
            '=' => {
                i += if prox == Some('=') { 2 } else { 1 };
                Tok::Cmp(Cmp::Ig)
            }
            '!' if prox == Some('=') => {
                i += 2;
                Tok::Cmp(Cmp::Dif)
            }
            '~' => {
                i += 1;
                Tok::Cmp(Cmp::Tem)
            }
            '!' if prox == Some('~') => {
                i += 2;
                Tok::Cmp(Cmp::NaoTem)
            }
            '!' => {
                i += 1;
                Tok::Nao
            }
            '<' => {
                i += if prox == Some('=') { 2 } else { 1 };
                Tok::Cmp(if prox == Some('=') { Cmp::MenIg } else { Cmp::Men })
            }
            '>' => {
                i += if prox == Some('=') { 2 } else { 1 };
                Tok::Cmp(if prox == Some('=') { Cmp::MaiIg } else { Cmp::Mai })
            }
            '\'' => {
                return Err(erro(ini + 1, "text is written in double quotes, not single quotes: for example FILTER==\"PASS\""));
            }
            '"' => {
                let mut s = String::new();
                i += 1;
                loop {
                    match c.get(i) {
                        None => {
                            return Err(erro(
                                ini + 1,
                                "the text opened by this double quote is never closed: add the closing `\"`",
                            ));
                        }
                        Some('"') => {
                            i += 1;
                            break;
                        }
                        // `\"` and `\\` put a quote or a backslash in the text.
                        Some('\\') if matches!(c.get(i + 1), Some('"') | Some('\\')) => {
                            s.push(c[i + 1]);
                            i += 2;
                        }
                        Some(&x) => {
                            s.push(x);
                            i += 1;
                        }
                    }
                }
                Tok::Txt(s)
            }
            '0'..='9' | '.' | '+' | '-' => {
                let mut j = i;
                if c[j] == '+' || c[j] == '-' {
                    j += 1;
                }
                let d0 = j;
                while j < n && c[j].is_ascii_digit() {
                    j += 1;
                }
                let mut digitos = j - d0;
                if j < n && c[j] == '.' {
                    let f0 = j + 1;
                    let mut k = f0;
                    while k < n && c[k].is_ascii_digit() {
                        k += 1;
                    }
                    digitos += k - f0;
                    j = k;
                }
                if digitos == 0 {
                    return Err(erro(
                        ini + 1,
                        &format!(
                            "`{}` does not start a number here; write the digits (for example -1.5), or \".\" in double quotes to test for a missing value",
                            c[ini]
                        ),
                    ));
                }
                if j < n && (c[j] == 'e' || c[j] == 'E') {
                    let mut k = j + 1;
                    if k < n && (c[k] == '+' || c[k] == '-') {
                        k += 1;
                    }
                    let e0 = k;
                    while k < n && c[k].is_ascii_digit() {
                        k += 1;
                    }
                    if k > e0 {
                        j = k;
                    }
                }
                let colado = j < n && (c[j].is_alphanumeric() || c[j] == '_' || c[j] == '.');
                while j < n && colado && (c[j].is_alphanumeric() || c[j] == '_' || c[j] == '.') {
                    j += 1;
                }
                let txt: String = c[i..j].iter().collect();
                match numero(&txt) {
                    Some(x) if !colado => {
                        i = j;
                        Tok::Num(x)
                    }
                    _ => {
                        return Err(erro(
                            ini + 1,
                            &format!("`{txt}` is not a number; if it is text, put it in double quotes"),
                        ));
                    }
                }
            }
            x if x.is_ascii_alphabetic() || x == '_' => {
                let mut j = i;
                while j < n && (c[j].is_ascii_alphanumeric() || c[j] == '_') {
                    j += 1;
                }
                let nome: String = c[i..j].iter().collect();
                if nome == "INFO" && c.get(j) == Some(&'/') {
                    let k0 = j + 1;
                    let mut k = k0;
                    while k < n && (c[k].is_ascii_alphanumeric() || c[k] == '_' || c[k] == '.') {
                        k += 1;
                    }
                    if k == k0 {
                        return Err(erro(
                            k0 + 1,
                            "expected the name of an INFO key right after `INFO/`, for example INFO/DP",
                        ));
                    }
                    i = k;
                    Tok::Campo(Campo::Info(c[k0..k].iter().collect()))
                } else if let Some(campo) = Campo::de_nome(&nome) {
                    i = j;
                    Tok::Campo(campo)
                } else {
                    let maiusculo = nome.to_ascii_uppercase();
                    let dica = if Campo::de_nome(&maiusculo).is_some() {
                        format!("field names are upper case: write {maiusculo}")
                    } else if (maiusculo == "FMT" || maiusculo == "FORMAT") && c.get(j) == Some(&'/') {
                        "per-sample (FORMAT) fields cannot be used in this filter".to_string()
                    } else if maiusculo == "INFO" {
                        "an INFO key is written INFO/KEY, for example INFO/DP".to_string()
                    } else {
                        format!("if it is an INFO key write INFO/{nome}; if it is text, put it in double quotes")
                    };
                    return Err(erro(ini + 1, &format!("unknown field `{nome}`: {dica}")));
                }
            }
            x => {
                return Err(erro(ini + 1, &format!("unexpected character `{x}`")));
            }
        };
        v.push(Token { tok, pos: ini + 1, fonte: c[ini..i].iter().collect() });
        if v.len() > MAX_TOKENS {
            return Err(erro(
                ini + 1,
                &format!("the expression is too long (more than {MAX_TOKENS} terms); split the filter into several runs"),
            ));
        }
    }
    v.push(Token { tok: Tok::Fim, pos: n + 1, fonte: String::new() });
    Ok(v)
}

// ---------------------------------------------------------------------
//  The expression language: syntax tree
// ---------------------------------------------------------------------

/// The value of TYPE in a comparison.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Classe {
    Snp,
    Mnp,
    Ins,
    Del,
    Indel,
    Outro,
    /// No ALT allele (ALT is `.`).
    Ref,
}

impl Classe {
    fn de_texto(s: &str) -> Option<Classe> {
        Some(match s {
            "snp" => Classe::Snp,
            "mnp" => Classe::Mnp,
            "ins" => Classe::Ins,
            "del" => Classe::Del,
            "indel" => Classe::Indel,
            "other" => Classe::Outro,
            "ref" => Classe::Ref,
            _ => return None,
        })
    }

    fn inclui(self, t: TipoAlelo) -> bool {
        match self {
            Classe::Snp => t == TipoAlelo::Snp,
            Classe::Mnp => t == TipoAlelo::Mnp,
            Classe::Ins => t == TipoAlelo::Ins,
            Classe::Del => t == TipoAlelo::Del,
            Classe::Indel => t == TipoAlelo::Ins || t == TipoAlelo::Del,
            Classe::Outro => t == TipoAlelo::Outro,
            Classe::Ref => false,
        }
    }
}

/// How TYPE is compared, as in bcftools.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ModoTipo {
    /// `==`: every ALT allele is of the type.
    Todos,
    /// `!=`: not every ALT allele is of the type.
    NemTodos,
    /// `~`: at least one ALT allele is of the type.
    Algum,
    /// `!~`: no ALT allele is of the type.
    Nenhum,
}

#[derive(Clone, Debug)]
enum No {
    E(Box<No>, Box<No>),
    Ou(Box<No>, Box<No>),
    Nao(Box<No>),
    /// Numeric field against a number.
    Num(Campo, Cmp, f64),
    /// Text field against text: `true` for ==, `false` for !=.
    Txt(Campo, bool, String),
    /// TYPE against a class.
    Tipo(ModoTipo, Classe),
    /// `FILTER ~ "A;B"` (`false`): A and B are among the filters of the
    /// record; `FILTER !~ "A;B"` (`true`): not both.
    FiltroTem(bool, Vec<String>),
    /// `FIELD == "."` (`true`: the value is missing) or `FIELD != "."`.
    Ausente(Campo, bool),
    /// `INFO/KEY` alone: the key or flag is present.
    Presente(String),
}

/// The value of an INFO key. `None`: the key is absent. `Some(None)`: present
/// without a value (a flag). `Some(Some(v))`: `KEY=v`.
fn info_valor<'a>(info: &'a str, chave: &str) -> Option<Option<&'a str>> {
    for item in info.split(';') {
        match item.split_once('=') {
            Some((k, v)) if k == chave => return Some(Some(v)),
            None if item == chave => return Some(None),
            _ => {}
        }
    }
    None
}

fn ausente(v: &str) -> bool {
    v == "." || v.is_empty()
}

impl No {
    fn avalia(&self, r: &Registro) -> bool {
        match self {
            No::E(a, b) => a.avalia(r) && b.avalia(r),
            No::Ou(a, b) => a.avalia(r) || b.avalia(r),
            No::Nao(a) => !a.avalia(r),
            // A missing value satisfies no comparison, `!=` included.
            No::Num(campo, op, v) => match campo {
                Campo::Qual => numero(r.qual).is_some_and(|q| op.vale(q, *v)),
                Campo::Pos => op.vale(r.pos as f64, *v),
                Campo::NAlt => op.vale(n_alt(r.alt) as f64, *v),
                // Several values: true when ANY of them satisfies the
                // comparison. A value that is not a number satisfies none.
                Campo::Info(k) => match info_valor(r.info, k) {
                    Some(Some(valor)) => valor.split(',').any(|x| numero(x).is_some_and(|x| op.vale(x, *v))),
                    _ => false,
                },
                _ => false,
            },
            No::Txt(campo, igual, s) => match campo {
                Campo::Info(k) => match info_valor(r.info, k) {
                    Some(Some(valor)) => valor.split(',').any(|x| !ausente(x) && (x == s) == *igual),
                    _ => false,
                },
                _ => {
                    let v = match campo {
                        Campo::Chrom => r.chrom,
                        Campo::Id => r.id,
                        Campo::Ref => r.ref_,
                        Campo::Alt => r.alt,
                        Campo::Filter => r.filter,
                        _ => ".",
                    };
                    !ausente(v) && (v == s) == *igual
                }
            },
            // As in bcftools: `==` asks that every ALT allele be of the type,
            // `~` that at least one be; `!=` and `!~` are their negations.
            // A record without ALT is of type "ref" and of no other.
            No::Tipo(modo, classe) => {
                let (todos, algum) = if *classe == Classe::Ref || sem_alt(r.alt) {
                    let v = *classe == Classe::Ref && sem_alt(r.alt);
                    (v, v)
                } else {
                    r.alt.split(',').fold((true, false), |(t, a), alelo| {
                        let x = classe.inclui(tipo_alelo(r.ref_, alelo));
                        (t && x, a || x)
                    })
                };
                match modo {
                    ModoTipo::Todos => todos,
                    ModoTipo::NemTodos => !todos,
                    ModoTipo::Algum => algum,
                    ModoTipo::Nenhum => !algum,
                }
            }
            No::FiltroTem(nega, pedidos) => {
                let tem = !ausente(r.filter) && pedidos.iter().all(|p| r.filter.split(';').any(|f| f == p));
                tem != *nega
            }
            No::Ausente(campo, falta) => {
                let f = match campo {
                    Campo::Qual => ausente(r.qual),
                    Campo::Pos | Campo::NAlt => false,
                    Campo::Chrom => ausente(r.chrom),
                    Campo::Id => ausente(r.id),
                    Campo::Ref => ausente(r.ref_),
                    Campo::Alt | Campo::Tipo => sem_alt(r.alt),
                    Campo::Filter => ausente(r.filter),
                    Campo::Info(k) => match info_valor(r.info, k) {
                        None => true,
                        Some(None) => false,
                        Some(Some(valor)) => valor.split(',').all(ausente),
                    },
                };
                f == *falta
            }
            No::Presente(k) => info_valor(r.info, k).is_some(),
        }
    }
}

// ---------------------------------------------------------------------
//  The expression language: recursive-descent parser
//
//      expression := and { ("||" | "|") and }
//      and        := not { ("&&" | "&") not }
//      not        := "!" not | primary
//      primary    := "(" expression ")" | comparison
//      comparison := operand [ operator operand ]
//      operand    := field | number | string
// ---------------------------------------------------------------------

struct Analisador<'a> {
    expr: &'a str,
    t: Vec<Token>,
    i: usize,
    nivel: usize,
}

impl Analisador<'_> {
    fn erro(&self, pos: usize, msg: &str) -> String {
        erro_expressao(self.expr, pos, msg)
    }

    /// The token under the cursor. The list always ends in `Fim`, and the
    /// cursor never moves past it.
    fn atual(&self) -> &Token {
        &self.t[self.i.min(self.t.len() - 1)]
    }

    fn desce(&mut self) -> Result<(), String> {
        self.nivel += 1;
        if self.nivel > MAX_NIVEL {
            return Err(self.erro(
                self.atual().pos,
                &format!("more than {MAX_NIVEL} nested `(` or `!`; simplify the expression"),
            ));
        }
        Ok(())
    }

    fn ou(&mut self) -> Result<No, String> {
        let mut e = self.e()?;
        while self.atual().tok == Tok::Ou {
            self.i += 1;
            e = No::Ou(Box::new(e), Box::new(self.e()?));
        }
        Ok(e)
    }

    fn e(&mut self) -> Result<No, String> {
        let mut e = self.nao()?;
        while self.atual().tok == Tok::E {
            self.i += 1;
            e = No::E(Box::new(e), Box::new(self.nao()?));
        }
        Ok(e)
    }

    fn nao(&mut self) -> Result<No, String> {
        if self.atual().tok == Tok::Nao {
            self.desce()?;
            self.i += 1;
            let dentro = self.nao()?;
            self.nivel -= 1;
            return Ok(No::Nao(Box::new(dentro)));
        }
        self.primario()
    }

    fn primario(&mut self) -> Result<No, String> {
        let pos = self.atual().pos;
        match self.atual().tok {
            Tok::Abre => {
                self.desce()?;
                self.i += 1;
                let dentro = self.ou()?;
                self.nivel -= 1;
                if self.atual().tok != Tok::Fecha {
                    return Err(self.erro(
                        self.atual().pos,
                        &format!(
                            "expected `)` to close the parenthesis opened at position {pos}, found {}",
                            descreve(self.atual())
                        ),
                    ));
                }
                self.i += 1;
                Ok(dentro)
            }
            Tok::Campo(_) | Tok::Num(_) | Tok::Txt(_) => self.comparacao(),
            _ => Err(self.erro(
                pos,
                &format!(
                    "expected a condition (a field such as QUAL or INFO/DP, `!` or `(`), found {}",
                    descreve(self.atual())
                ),
            )),
        }
    }

    fn comparacao(&mut self) -> Result<No, String> {
        let esq = self.atual().clone();
        self.i += 1;
        let Tok::Cmp(op) = self.atual().tok else {
            return match esq.tok {
                Tok::Campo(Campo::Info(k)) => Ok(No::Presente(k)),
                Tok::Campo(c) => Err(self.erro(
                    self.atual().pos,
                    &format!(
                        "expected a comparison operator (== != < <= > >= ~ !~) after `{}`, found {}; only an INFO key can stand alone, as a test that the key is present",
                        c.nome(),
                        descreve(self.atual())
                    ),
                )),
                _ => Err(self.erro(
                    esq.pos,
                    &format!("`{}` alone is not a condition: compare it with a field, for example QUAL>30", esq.fonte),
                )),
            };
        };
        let pos_op = self.atual().pos;
        self.i += 1;
        let dir = self.atual().clone();
        if !matches!(dir.tok, Tok::Campo(_) | Tok::Num(_) | Tok::Txt(_)) {
            return Err(self.erro(
                dir.pos,
                &format!(
                    "expected a number, text in double quotes or a field after `{}`, found {}",
                    op.simbolo(),
                    descreve(&dir)
                ),
            ));
        }
        self.i += 1;
        match (esq.tok, dir.tok) {
            (Tok::Campo(_), Tok::Campo(_)) => Err(self.erro(
                pos_op,
                "comparing two fields is not supported: one side must be a number or text in double quotes",
            )),
            (Tok::Campo(c), literal) => self.monta(c, op, op, literal, pos_op, dir.pos, &dir.fonte),
            (literal, Tok::Campo(c)) => self.monta(c, op.invertido(), op, literal, pos_op, esq.pos, &esq.fonte),
            _ => Err(self.erro(
                esq.pos,
                "a comparison needs a field on one side: QUAL, POS, N_ALT, CHROM, ID, REF, ALT, FILTER, TYPE or INFO/KEY",
            )),
        }
    }

    /// Checks the types and builds the comparison node. `escrito` is the
    /// operator as the user wrote it (`op` is already swapped when the
    /// literal came first).
    #[allow(clippy::too_many_arguments)]
    fn monta(
        &self,
        campo: Campo,
        op: Cmp,
        escrito: Cmp,
        literal: Tok,
        pos_op: usize,
        pos_lit: usize,
        fonte_lit: &str,
    ) -> Result<No, String> {
        let so_igualdade = |o_que: &str| match op {
            Cmp::Ig => Ok(true),
            Cmp::Dif => Ok(false),
            _ => Err(self.erro(
                pos_op,
                &format!("{o_que} can only be compared with == or != (found `{}`)", escrito.simbolo()),
            )),
        };
        let contem = matches!(op, Cmp::Tem | Cmp::NaoTem);
        match (campo, literal) {
            (campo, Tok::Txt(s)) if s == "." => Ok(No::Ausente(campo, so_igualdade("a missing value (\".\")")?)),
            (Campo::Tipo, Tok::Txt(s)) => match Classe::de_texto(&s) {
                Some(classe) => Ok(No::Tipo(
                    match op {
                        Cmp::Ig => ModoTipo::Todos,
                        Cmp::Dif => ModoTipo::NemTodos,
                        Cmp::Tem => ModoTipo::Algum,
                        Cmp::NaoTem => ModoTipo::Nenhum,
                        _ => {
                            return Err(self.erro(
                                pos_op,
                                &format!(
                                    "TYPE can only be compared with == (every ALT allele), ~ (at least one), != or !~ (found `{}`)",
                                    escrito.simbolo()
                                ),
                            ))
                        }
                    },
                    classe,
                )),
                None => Err(self.erro(
                    pos_lit,
                    &format!("unknown variant type \"{s}\": TYPE is compared with {VALORES_TYPE} (lower case)"),
                )),
            },
            (Campo::Tipo, _) => Err(self.erro(
                pos_lit,
                &format!("TYPE is compared with {VALORES_TYPE}, in double quotes"),
            )),
            (Campo::Filter, Tok::Txt(s)) if contem => {
                Ok(No::FiltroTem(op == Cmp::NaoTem, s.split(';').filter(|f| !f.is_empty()).map(str::to_string).collect()))
            }
            (campo, _) if contem => Err(self.erro(
                pos_op,
                &format!(
                    "`{}` works only with TYPE (TYPE~\"indel\": at least one ALT allele is an indel) and FILTER (FILTER~\"q10\": q10 is among the filters); {} takes == or != (regular expressions are not supported)",
                    escrito.simbolo(),
                    campo.nome()
                ),
            )),
            (campo, Tok::Num(v)) if campo.numerico() || matches!(campo, Campo::Info(_)) => Ok(No::Num(campo, op, v)),
            (campo, Tok::Num(_)) => Err(self.erro(
                pos_lit,
                &format!(
                    "{0} holds text: put the value in double quotes, {0}==\"{fonte_lit}\"",
                    campo.nome()
                ),
            )),
            (campo, Tok::Txt(_)) if campo.numerico() => Err(self.erro(
                pos_lit,
                &format!(
                    "{} is a number: compare it with a number without quotes (or with \".\" to test for a missing value)",
                    campo.nome()
                ),
            )),
            (campo, Tok::Txt(s)) => Ok(No::Txt(campo, so_igualdade("text")?, s)),
            (_, _) => Err(self.erro(pos_lit, "expected a number or text in double quotes")),
        }
    }
}

/// A filter expression, compiled once and evaluated on every record.
///
/// A deliberately small subset of the `-i` expressions of bcftools:
///
/// - fields: `QUAL`, `POS`, `N_ALT` (numbers); `CHROM`, `ID`, `REF`, `ALT`,
///   `FILTER` (text, compared with the whole column, exactly and
///   case-sensitively); `TYPE` (`"snp"`, `"mnp"`, `"ins"`, `"del"`, `"indel"`,
///   `"other"`, `"ref"`); `INFO/KEY`;
/// - comparison: `== != < <= > >=` (`=` is `==`); text takes only `==` and `!=`;
/// - TYPE, as in bcftools: `TYPE=="snp"` is true when EVERY ALT allele is a
///   SNP, `TYPE~"snp"` when AT LEAST ONE is; `!=` and `!~` are the negations.
///   A record without ALT is `"ref"` and nothing else;
/// - `FILTER~"A"` (and `FILTER~"A;B"`): the filters named are among those of
///   the record, as in bcftools; `FILTER!~"A"` is the negation (true for a
///   missing FILTER). `~` is not a regular expression and works only with
///   TYPE and FILTER;
/// - logic: `&& || !` and parentheses (`&` and `|` are synonyms); `!` binds
///   tighter than `&&`, which binds tighter than `||`;
/// - a missing value (`.`, or an INFO key that is absent) satisfies no
///   comparison, `!=` included. `FIELD=="."` tests that the value is missing,
///   `FIELD!="."` that it is there; `INFO/KEY` alone tests that the key is
///   present;
/// - an INFO value with several comma-separated entries satisfies a
///   comparison when ANY entry does (TYPE is the exception: see above).
#[derive(Clone, Debug)]
pub struct Expressao {
    raiz: No,
}

impl Expressao {
    /// Compiles the expression. The error says what was expected and at which
    /// character position (from 1), and shows the expression.
    pub fn compila(texto: &str) -> Result<Expressao, String> {
        let t = tokeniza(texto)?;
        let mut a = Analisador { expr: texto, t, i: 0, nivel: 0 };
        if a.atual().tok == Tok::Fim {
            return Err(a.erro(a.atual().pos, "the expression is empty: expected a condition such as QUAL>30"));
        }
        let raiz = a.ou()?;
        if a.atual().tok != Tok::Fim {
            let sobra = if a.atual().tok == Tok::Fecha { " (this `)` has no matching `(`)" } else { "" };
            return Err(a.erro(
                a.atual().pos,
                &format!("expected `&&`, `||` or the end of the expression, found {}{sobra}", descreve(a.atual())),
            ));
        }
        Ok(Expressao { raiz })
    }

    fn aceita(&self, r: &Registro) -> bool {
        self.raiz.avalia(r)
    }
}

// ---------------------------------------------------------------------
//  Regions
// ---------------------------------------------------------------------

const FORMAS_REGIAO: &str = "Accepted forms: `chr1` (the whole sequence), `chr1:1000-2000` (1-based, both ends included) and `chr1:1500` (one position); numbers may use `,` as thousands separator, and a sequence name that contains `:` goes in braces, `{HLA-A*01:01}:1-500`";

/// One region as (sequence, start, end), 1-based and inclusive.
fn analisa_regiao(texto: &str) -> Result<(String, u64, u64), String> {
    let erro = |motivo: &str| format!("Invalid region `{texto}`: {motivo}. {FORMAS_REGIAO}.");
    let s = texto.trim();
    if s.is_empty() {
        return Err(erro("it is empty"));
    }
    let (nome, faixa) = if let Some(r) = s.strip_prefix('{') {
        let Some(fim) = r.find('}') else {
            return Err(erro("the `{` is never closed"));
        };
        match &r[fim + 1..] {
            "" => (&r[..fim], None),
            depois => match depois.strip_prefix(':') {
                Some(f) => (&r[..fim], Some(f)),
                None => return Err(erro("after `}` there must be nothing, or `:` and the positions")),
            },
        }
    } else {
        match s.rsplit_once(':') {
            Some((n, f)) => (n, Some(f)),
            None => (s, None),
        }
    };
    if nome.is_empty() {
        return Err(erro("the sequence name is missing"));
    }
    let num = |x: &str| -> Option<u64> {
        let limpo: String = x.trim().chars().filter(|c| *c != ',').collect();
        if limpo.is_empty() || !limpo.bytes().all(|b| b.is_ascii_digit()) {
            return None;
        }
        limpo.parse().ok()
    };
    let (ini, fim) = match faixa {
        None => (0, u64::MAX),
        Some(f) => match f.split_once('-') {
            None => {
                let Some(p) = num(f) else {
                    return Err(erro(&format!("`{f}` is not a position")));
                };
                (p, p)
            }
            Some((a, b)) => {
                if a.trim().is_empty() || b.trim().is_empty() {
                    return Err(erro("open ranges (`start-` or `-end`) are not supported: give both ends"));
                }
                let (Some(x), Some(y)) = (num(a), num(b)) else {
                    return Err(erro(&format!("`{a}` and `{b}` must both be positions (whole numbers)")));
                };
                if x > y {
                    return Err(erro(&format!("the start ({x}) is greater than the end ({y})")));
                }
                (x, y)
            }
        },
    };
    Ok((nome.to_string(), ini, fim))
}

/// The regions of a filter: per sequence, disjoint intervals in order, so a
/// record is tested with one binary search however many regions were given.
struct Regioes {
    por_nome: HashMap<String, Vec<(u64, u64)>>,
    /// Each region in canonical form, in the order given, for the header.
    textos: Vec<String>,
}

impl Regioes {
    fn nova(textos: &[String]) -> Result<Regioes, String> {
        let mut r = Regioes { por_nome: HashMap::new(), textos: Vec::new() };
        for t in textos {
            let (nome, ini, fim) = analisa_regiao(t)?;
            r.textos.push(if fim == u64::MAX {
                nome.clone()
            } else if ini == fim {
                format!("{nome}:{ini}")
            } else {
                format!("{nome}:{ini}-{fim}")
            });
            r.por_nome.entry(nome).or_default().push((ini, fim));
        }
        // "At least one region" only needs the union: merge what overlaps or
        // touches.
        for v in r.por_nome.values_mut() {
            v.sort_unstable();
            let mut unidos: Vec<(u64, u64)> = Vec::with_capacity(v.len());
            for &(ini, fim) in v.iter() {
                match unidos.last_mut() {
                    Some(u) if ini <= u.1.saturating_add(1) => u.1 = u.1.max(fim),
                    _ => unidos.push((ini, fim)),
                }
            }
            *v = unidos;
        }
        Ok(r)
    }

    /// Does the span `ini..=fim` on `chrom` overlap any region?
    fn contem(&self, chrom: &str, ini: u64, fim: u64) -> bool {
        let Some(v) = self.por_nome.get(chrom) else {
            return false;
        };
        // The intervals are disjoint and in order, so their ends are in
        // order too: the first one that does not end before `ini` is the only
        // candidate.
        let i = v.partition_point(|&(_, f)| f < ini);
        i < v.len() && v[i].0 <= fim
    }
}

// ---------------------------------------------------------------------
//  The filter
// ---------------------------------------------------------------------

/// What `vcf_filter` did.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FilterSummary {
    /// Records read.
    pub lidos: usize,
    /// Records written.
    pub mantidos: usize,
    /// Sample columns in the input and in the output.
    pub amostras_antes: usize,
    pub amostras_depois: usize,
    /// Lines with fewer than 8 columns or a non-numeric POS, left out (they
    /// are not in `lidos`).
    pub ignoradas: usize,
}

impl FilterSummary {
    pub fn json(&self) -> String {
        format!(
            "{{\"read\":{},\"kept\":{},\"samples_before\":{},\"samples_after\":{},\"skipped\":{}}}",
            self.lidos, self.mantidos, self.amostras_antes, self.amostras_depois, self.ignoradas
        )
    }
}

/// `s` in double quotes, with `"` and `\` escaped: a value that a VCF header
/// parser reads back whole, whatever it contains.
fn entre_aspas(s: &str) -> String {
    let mut o = String::with_capacity(s.len() + 2);
    o.push('"');
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            c if c.is_control() => o.push(' '),
            c => o.push(c),
        }
    }
    o.push('"');
    o
}

/// Where each requested sample is among the sample columns of the file.
fn localiza_amostras(pedidas: &[String], nomes: &[&str], entrada: &str) -> Result<Vec<usize>, String> {
    let mut idx = Vec::with_capacity(pedidas.len());
    let mut desconhecidas: Vec<&str> = Vec::new();
    for (i, p) in pedidas.iter().enumerate() {
        if pedidas[..i].contains(p) {
            return Err(format!(
                "Sample `{p}` is listed more than once; each sample can be selected only once (the output would have two columns with the same name)."
            ));
        }
        match nomes.iter().position(|n| n == p) {
            Some(j) => idx.push(j),
            None => desconhecidas.push(p),
        }
    }
    if desconhecidas.is_empty() {
        return Ok(idx);
    }
    const MOSTRA: usize = 200;
    let disponiveis = if nomes.is_empty() {
        "the file has no sample columns".to_string()
    } else {
        let mais = nomes.len().saturating_sub(MOSTRA);
        format!(
            "samples available ({}): {}{}",
            nomes.len(),
            nomes.iter().take(MOSTRA).copied().collect::<Vec<_>>().join(", "),
            if mais > 0 { format!(", ... and {mais} more") } else { String::new() }
        )
    };
    Err(format!(
        "Sample(s) not found in {entrada}: {}. Names are case-sensitive; {disponiveis}.",
        desconhecidas.iter().take(MOSTRA).map(|d| format!("`{d}`")).collect::<Vec<_>>().join(", ")
    ))
}

/// Column `j` of `resto` (0 is FORMAT), given the positions of its tabs.
fn coluna<'a>(resto: &'a str, cortes: &[usize], j: usize) -> Option<&'a str> {
    if j > cortes.len() {
        return None;
    }
    let ini = if j == 0 { 0 } else { cortes[j - 1] + 1 };
    let fim = if j < cortes.len() { cortes[j] } else { resto.len() };
    Some(&resto[ini..fim])
}

/// Filters a VCF: keeps the records that satisfy the expression (if given)
/// AND overlap at least one of the regions (if any is given); then, if
/// `amostras` is given, keeps only those sample columns, in the order given.
///
/// - `expressao`: see `Expressao`. A blank expression is no expression.
/// - `regioes`: `name`, `name:start-end` or `name:pos`, 1-based and inclusive.
///   A record matches when its span `POS ..= POS + len(REF) - 1` overlaps the
///   region. (INFO/END of symbolic alleles is not considered.)
/// - `amostras`: `None` keeps every column. `Some(&[])` writes sites only: no
///   FORMAT and no sample columns. An unknown name is an error.
///
/// One line is added just before `#CHROM` recording what was done:
/// `##MMGEN_filter=<expression="...";regions="...";samples="...">`, with only
/// the parts that were used. The values are in double quotes because an
/// expression contains `>` and `,`, which would end an unquoted header value.
pub fn vcf_filter(
    entrada: &str,
    saida: &str,
    expressao: Option<&str>,
    regioes: &[String],
    amostras: Option<&[String]>,
) -> Result<FilterSummary, String> {
    recusa_mesmo_arquivo(entrada, saida)?;
    // Everything that can be wrong in the request is checked before the
    // output is created.
    let texto_expr = expressao.map(str::trim).filter(|e| !e.is_empty());
    let expr = match texto_expr {
        Some(e) => Some(Expressao::compila(e)?),
        None => None,
    };
    let reg = if regioes.is_empty() { None } else { Some(Regioes::nova(regioes)?) };

    let mut marca: Vec<String> = Vec::new();
    if let Some(e) = texto_expr {
        marca.push(format!("expression={}", entre_aspas(e)));
    }
    if let Some(g) = &reg {
        marca.push(format!("regions={}", entre_aspas(&g.textos.join(","))));
    }
    if let Some(a) = amostras {
        marca.push(format!("samples={}", entre_aspas(&a.join(","))));
    }
    let marca = format!("##MMGEN_filter=<{}>", marca.join(";"));

    let mut l = Leitor::abre(entrada)?;
    let existia = std::path::Path::new(saida).exists();
    let mut guarda = Saida { caminho: saida, armada: false };
    let mut w = Escritor::cria_para(saida, Formato::Vcf)?;
    guarda.armada = !existia;

    let mut s = FilterSummary::default();
    let mut viu_chrom = false;
    // Positions, among the sample columns, of the samples to keep.
    let mut escolhidas: Vec<usize> = Vec::new();
    let mut cortes: Vec<usize> = Vec::new();
    let mut nova = String::new();

    while l.avanca()? {
        let t = l.texto();
        if t.is_empty() {
            continue;
        }
        if t.starts_with('#') {
            if viu_chrom || !t.starts_with("#CHROM") {
                w.linha(t)?;
                continue;
            }
            viu_chrom = true;
            w.linha(&marca)?;
            let nomes: Vec<&str> = t.split('\t').skip(9).collect();
            s.amostras_antes = nomes.len();
            match amostras {
                None => {
                    s.amostras_depois = nomes.len();
                    w.linha(t)?;
                }
                Some(pedidas) => {
                    escolhidas = localiza_amostras(pedidas, &nomes, entrada)?;
                    s.amostras_depois = pedidas.len();
                    nova.clear();
                    for (i, c) in t.split('\t').take(8).enumerate() {
                        if i > 0 {
                            nova.push('\t');
                        }
                        nova.push_str(c);
                    }
                    if !pedidas.is_empty() {
                        nova.push_str("\tFORMAT");
                        for p in pedidas {
                            nova.push('\t');
                            nova.push_str(p);
                        }
                    }
                    w.linha(&nova)?;
                }
            }
            continue;
        }
        if !viu_chrom && amostras.is_some_and(|a| !a.is_empty()) {
            return Err(format!(
                "{entrada} has no `#CHROM` header line before its first record (line {}), so its sample columns have no names and cannot be selected. Add the header line, or pass an empty sample list to keep only the sites.",
                l.linha
            ));
        }
        let Some(r) = registro(t) else {
            s.ignoradas += 1;
            continue;
        };
        s.lidos += 1;
        if expr.as_ref().is_some_and(|e| !e.aceita(&r)) {
            continue;
        }
        if let Some(g) = &reg {
            let fim = r.pos.saturating_add(r.ref_.len().max(1) as u64 - 1);
            if !g.contem(r.chrom, r.pos, fim) {
                continue;
            }
        }
        s.mantidos += 1;
        let Some(pedidas) = amostras else {
            w.linha(t)?;
            continue;
        };
        // The eight fixed columns are a prefix of the line.
        let fixas = match r.resto {
            Some(resto) => t.len() - resto.len() - 1,
            None => t.len(),
        };
        nova.clear();
        nova.push_str(&t[..fixas]);
        if !pedidas.is_empty() {
            let resto = r.resto.unwrap_or("");
            cortes.clear();
            cortes.extend(resto.bytes().enumerate().filter(|(_, b)| *b == b'\t').map(|(i, _)| i));
            nova.push('\t');
            nova.push_str(coluna(resto, &cortes, 0).filter(|f| !f.is_empty()).unwrap_or("."));
            for &j in &escolhidas {
                nova.push('\t');
                // A record with fewer sample columns than the header gets `.`.
                nova.push_str(coluna(resto, &cortes, j + 1).unwrap_or("."));
            }
        }
        w.linha(&nova)?;
    }
    w.termina()?;
    guarda.armada = false;
    Ok(s)
}

// =====================================================================
//  C) NORMALIZATION
// =====================================================================

/// What `vcf_normalize` did.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NormSummary {
    /// Records read.
    pub lidos: usize,
    /// Records written (more than `lidos` when sites were split).
    pub gravados: usize,
    /// Input records that were split into one record per ALT allele.
    pub separados: usize,
    /// Output records whose REF and ALT were trimmed (got shorter).
    pub aparados: usize,
    /// Output records whose POS changed (trimming moves it forward,
    /// left-alignment back).
    pub pos_alteradas: usize,
    /// Output records moved to the left by the left-alignment (only with a
    /// reference).
    pub alinhados: usize,
    /// Records whose REF is not the reference at their position: neither
    /// checked nor aligned, only trimmed (only with a reference).
    pub ref_divergente: usize,
    /// Records on a sequence the reference does not have, under either naming
    /// (`chr1` / `1`): only trimmed (only with a reference).
    pub fora_da_referencia: usize,
    /// The output had to be sorted again (some POS changed).
    pub reordenado: bool,
    /// Lines with fewer than 8 columns or a non-numeric POS, left out (they
    /// are not in `lidos`).
    pub ignoradas: usize,
    /// Runs the sort wrote to disk. Zero when it was not needed or fit in
    /// memory.
    pub lotes: usize,
}

impl NormSummary {
    pub fn json(&self) -> String {
        format!(
            "{{\"read\":{},\"written\":{},\"split\":{},\"trimmed\":{},\"position_changed\":{},\"left_aligned\":{},\"ref_mismatch\":{},\"not_in_reference\":{},\"resorted\":{},\"skipped\":{},\"runs\":{}}}",
            self.lidos,
            self.gravados,
            self.separados,
            self.aparados,
            self.pos_alteradas,
            self.alinhados,
            self.ref_divergente,
            self.fora_da_referencia,
            self.reordenado,
            self.ignoradas,
            self.lotes
        )
    }
}

/// How the values of an INFO or FORMAT key follow the alleles.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Numero {
    /// One value per ALT allele.
    A,
    /// One value per allele, REF first.
    R,
    /// One value per genotype.
    G,
}

/// `##INFO=<ID=AF,Number=A,...>` gives `("AF", Numero::A)`. `None` when the
/// line is not such a declaration or its Number is not A, R or G.
fn declaracao(linha: &str, prefixo: &str) -> Option<(String, Numero)> {
    let c = linha.strip_prefix(prefixo)?;
    let c = c.strip_suffix('>').unwrap_or(c);
    let b = c.as_bytes();
    let (mut id, mut num): (Option<&str>, Option<&str>) = (None, None);
    let (mut ini, mut i, mut aspas) = (0usize, 0usize, false);
    // `key=value` pairs separated by commas; a comma inside a quoted
    // Description does not separate anything.
    loop {
        if i >= b.len() || (b[i] == b',' && !aspas) {
            let fim = i.min(b.len());
            if let Some((k, v)) = c[ini..fim].split_once('=') {
                match k.trim() {
                    "ID" if id.is_none() => id = Some(v.trim()),
                    "Number" if num.is_none() => num = Some(v.trim()),
                    _ => {}
                }
            }
            if i >= b.len() {
                break;
            }
            ini = i + 1;
        } else if b[i] == b'"' {
            aspas = !aspas;
        } else if b[i] == b'\\' && aspas {
            i += 1;
        }
        i += 1;
    }
    let numero = match num? {
        "A" => Numero::A,
        "R" => Numero::R,
        "G" => Numero::G,
        _ => return None,
    };
    Some((id?.to_string(), numero))
}

/// Appends to `sai` the values of `valor` that belong to ALT allele `k`
/// (from 1) of `n`. `false`, with nothing appended, when the number of values
/// is not what `num` implies: the caller then copies the value unchanged.
fn reparte(valor: &str, num: Numero, k: usize, n: usize, sai: &mut String) -> bool {
    let quantos = valor.bytes().filter(|b| *b == b',').count() + 1;
    let item = |j: usize| valor.split(',').nth(j).unwrap_or(".");
    match num {
        Numero::A if quantos == n => sai.push_str(item(k - 1)),
        Numero::R if quantos == n + 1 => {
            sai.push_str(item(0));
            sai.push(',');
            sai.push_str(item(k));
        }
        // Diploid genotypes in VCF order: 0/0, 0/1, 1/1, 0/2, 1/2, 2/2, ...
        // Genotype a/b (a <= b) is at b(b+1)/2 + a.
        Numero::G if quantos == (n + 1) * (n + 2) / 2 => {
            let base = k * (k + 1) / 2;
            for (i, j) in [0, base, base + k].into_iter().enumerate() {
                if i > 0 {
                    sai.push(',');
                }
                sai.push_str(item(j));
            }
        }
        _ => return false,
    }
    true
}

/// The INFO column of the record that keeps only ALT allele `k` of `n`.
fn info_do_alelo(info: &str, decl: &HashMap<String, Numero>, k: usize, n: usize, sai: &mut String) {
    for (i, item) in info.split(';').enumerate() {
        if i > 0 {
            sai.push(';');
        }
        match item.split_once('=') {
            Some((chave, valor)) => match decl.get(chave) {
                Some(&num) => {
                    sai.push_str(chave);
                    sai.push('=');
                    if !reparte(valor, num, k, n, sai) {
                        sai.push_str(valor);
                    }
                }
                None => sai.push_str(item),
            },
            None => sai.push_str(item),
        }
    }
}

/// GT for the record that keeps only ALT allele `k`: `k` becomes 1, 0 stays
/// 0, any other ALT allele becomes `.` -- that sample carries an allele this
/// record no longer describes, and calling it reference would be a false
/// statement. The separators (`/`, `|`) are kept.
fn recodifica_gt(gt: &str, k: usize, sai: &mut String) {
    let alelo = |a: &str, sai: &mut String| match a.parse::<usize>() {
        Ok(0) => sai.push('0'),
        Ok(x) if x == k => sai.push('1'),
        Ok(_) => sai.push('.'),
        Err(_) => sai.push_str(a),
    };
    let mut ini = 0usize;
    for (i, b) in gt.bytes().enumerate() {
        if b == b'/' || b == b'|' {
            alelo(&gt[ini..i], sai);
            sai.push(b as char);
            ini = i + 1;
        }
    }
    alelo(&gt[ini..], sai);
}

enum Acao {
    Copia,
    Gt,
    Reparte(Numero),
}

/// FORMAT and the sample columns of the record that keeps only ALT allele `k`.
fn amostras_do_alelo(
    resto: &str,
    decl: &HashMap<String, Numero>,
    k: usize,
    n: usize,
    acoes: &mut Vec<Acao>,
    sai: &mut String,
) {
    let mut cols = resto.split('\t');
    let formato = cols.next().unwrap_or("");
    acoes.clear();
    acoes.extend(formato.split(':').map(|chave| match decl.get(chave) {
        _ if chave == "GT" => Acao::Gt,
        Some(&num) => Acao::Reparte(num),
        None => Acao::Copia,
    }));
    let so_copia = acoes.iter().all(|a| matches!(a, Acao::Copia));
    sai.push_str(formato);
    for amostra in cols {
        sai.push('\t');
        if so_copia {
            sai.push_str(amostra);
            continue;
        }
        for (i, campo) in amostra.split(':').enumerate() {
            if i > 0 {
                sai.push(':');
            }
            match acoes.get(i) {
                Some(Acao::Gt) => recodifica_gt(campo, k, sai),
                Some(Acao::Reparte(num)) => {
                    if !reparte(campo, *num, k, n, sai) {
                        sai.push_str(campo);
                    }
                }
                _ => sai.push_str(campo),
            }
        }
    }
}

/// How many bases to drop from the start and from the end of both alleles to
/// reach the minimal representation: first the shared end, then the shared
/// start, always leaving at least one base in each.
fn apara(ref_: &[u8], alt: &[u8]) -> (usize, usize) {
    let (mut r, mut a) = (ref_.len(), alt.len());
    while r > 1 && a > 1 && ref_[r - 1].eq_ignore_ascii_case(&alt[a - 1]) {
        r -= 1;
        a -= 1;
    }
    let mut esq = 0usize;
    while r - esq > 1 && a - esq > 1 && ref_[esq].eq_ignore_ascii_case(&alt[esq]) {
        esq += 1;
    }
    (esq, ref_.len() - r)
}

/// What `alinha` found.
enum Alinhamento {
    /// The variant in its leftmost, shortest form: POS (1-based), REF, ALT.
    Feito(u64, Vec<u8>, Vec<u8>),
    /// REF is not the reference at POS.
    RefDiferente,
    /// The reference has no such sequence.
    ForaDaReferencia,
}

/// Left-aligns and trims one variant with a single ALT allele against the
/// reference, by the algorithm of vt normalize and bcftools norm (Tan,
/// Abecasis and Kang, Bioinformatics 2015): while both alleles end with the
/// same base, drop it; when that leaves an allele empty, put the reference
/// base to the left in front of both and move POS one back; at the end, drop
/// the bases both share at the start, leaving at least one in each. An indel
/// inside a repeat thus ends at the leftmost position that gives the same
/// sequence. At position 1 there is nothing to the left, and it stops.
/// `pos` is 1-based and at least 1; both alleles are made of ACGTN.
fn alinha(g: &mut Genoma, chrom: &str, pos: u64, ref_: &str, alt: &str) -> Result<Alinhamento, String> {
    let Some(nome) = g.nome(chrom) else {
        return Ok(Alinhamento::ForaDaReferencia);
    };
    let p0 = pos - 1;
    let seg = g.trecho(&nome, p0, p0 + ref_.len() as u64)?;
    if seg.len() != ref_.len() || !seg.eq_ignore_ascii_case(ref_.as_bytes()) {
        return Ok(Alinhamento::RefDiferente);
    }
    let mut r: VecDeque<u8> = ref_.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let mut a: VecDeque<u8> = alt.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let mut p = p0;
    loop {
        let mut mudou = false;
        if r.back().is_some() && r.back() == a.back() {
            if (r.len() == 1 || a.len() == 1) && p == 0 {
                break;
            }
            r.pop_back();
            a.pop_back();
            mudou = true;
        }
        if r.is_empty() || a.is_empty() {
            p -= 1;
            let b = g.trecho(&nome, p, p + 1)?.first().copied().unwrap_or(b'N');
            r.push_front(b);
            a.push_front(b);
            mudou = true;
        }
        if !mudou {
            break;
        }
    }
    let mut esq = 0usize;
    while r.len() - esq >= 2 && a.len() - esq >= 2 && r[esq] == a[esq] {
        esq += 1;
    }
    Ok(Alinhamento::Feito(p + esq as u64 + 1, r.range(esq..).copied().collect(), a.range(esq..).copied().collect()))
}

/// Writes one record: `alt` in the ALT column and `cauda` (QUAL onwards)
/// after it, trimmed when it has a single plain ALT allele -- and, with a
/// reference, left-aligned. `original` is the input line when nothing else
/// in it changed, to be copied as it is.
#[allow(clippy::too_many_arguments)]
fn grava_registro(
    w: &mut Escritor,
    s: &mut NormSummary,
    linha: &mut String,
    r: &Registro,
    alt: &str,
    cauda: &str,
    original: Option<&str>,
    genoma: Option<&mut Genoma>,
) -> Result<(), String> {
    s.gravados += 1;
    // Symbolic alleles, breakends, `*`, `.` and lists of alleles are not
    // plain bases and are left alone.
    let bases = so_bases(r.ref_) && so_bases(alt);
    // An ALT equal to REF is not a variant: there is nothing to align.
    if let (Some(g), true) = (genoma, bases && r.pos >= 1 && !r.ref_.eq_ignore_ascii_case(alt)) {
        match alinha(g, r.chrom, r.pos, r.ref_, alt)? {
            Alinhamento::ForaDaReferencia => s.fora_da_referencia += 1,
            Alinhamento::RefDiferente => s.ref_divergente += 1,
            Alinhamento::Feito(pos, nr, na) => {
                let mesmo = pos == r.pos && nr.eq_ignore_ascii_case(r.ref_.as_bytes()) && na.eq_ignore_ascii_case(alt.as_bytes());
                if mesmo {
                    if let Some(t) = original {
                        return w.linha(t);
                    }
                    return monta_registro(w, linha, r, r.pos_txt, r.ref_, alt, cauda);
                }
                if nr.len() + na.len() < r.ref_.len() + alt.len() {
                    s.aparados += 1;
                }
                if pos < r.pos {
                    s.alinhados += 1;
                }
                if pos != r.pos {
                    s.pos_alteradas += 1;
                }
                // Made of the alleles and of reference bases: ASCII.
                let (nr, na) = (String::from_utf8_lossy(&nr), String::from_utf8_lossy(&na));
                return monta_registro(w, linha, r, &pos.to_string(), &nr, &na, cauda);
            }
        }
    }
    let (esq, dir) = if r.ref_.len() > 1 && alt.len() > 1 && bases {
        apara(r.ref_.as_bytes(), alt.as_bytes())
    } else {
        (0, 0)
    };
    if esq + dir == 0 {
        if let Some(t) = original {
            return w.linha(t);
        }
    } else {
        s.aparados += 1;
        if esq > 0 {
            s.pos_alteradas += 1;
        }
    }
    let pos = if esq == 0 { r.pos_txt.to_string() } else { r.pos.saturating_add(esq as u64).to_string() };
    monta_registro(w, linha, r, &pos, &r.ref_[esq..r.ref_.len() - dir], &alt[esq..alt.len() - dir], cauda)
}

fn monta_registro(
    w: &mut Escritor,
    linha: &mut String,
    r: &Registro,
    pos: &str,
    ref_: &str,
    alt: &str,
    cauda: &str,
) -> Result<(), String> {
    linha.clear();
    linha.push_str(r.chrom);
    linha.push('\t');
    linha.push_str(pos);
    linha.push('\t');
    linha.push_str(r.id);
    linha.push('\t');
    linha.push_str(ref_);
    linha.push('\t');
    linha.push_str(alt);
    linha.push('\t');
    linha.push_str(cauda);
    w.linha(linha)
}

/// `.gz`, `.lz4` or nothing: the compression `Escritor` picks for this path.
fn compressao(caminho: &str) -> &'static str {
    // Every container the writer knows: the temporary file must be in the
    // same one as the output, because it may become the output by a rename.
    let c = caminho.to_ascii_lowercase();
    [".gz", ".bgz", ".lz4", ".parquet", ".bam"].into_iter().find(|s| c.ends_with(s)).unwrap_or("")
}

/// Puts the finished temporary file at its destination: a rename when both
/// are on the same volume, a copy otherwise (the temporary copy is removed
/// with its directory).
fn move_arquivo(de: &str, para: &str) -> Result<(), String> {
    if std::fs::rename(de, para).is_ok() {
        return Ok(());
    }
    std::fs::copy(de, para).map(|_| ()).map_err(|e| format!("Error writing {para}: {e}"))
}

/// Normalizes a VCF: optionally splits multiallelic sites, then trims every
/// variant to its minimal representation.
///
/// 1. SPLITTING (only with `separar_multialelicos`). A record with n > 1 ALT
///    alleles becomes n records, one per allele, in ALT order. INFO and
///    FORMAT keys declared in the header as `Number=A`, `Number=R` or
///    `Number=G` (diploid) keep only the values of that allele; every other
///    key -- and any key whose number of values is not what its Number
///    implies -- is copied unchanged. GT is recoded: the allele of the record
///    becomes 1, the reference stays 0, and any OTHER ALT allele becomes `.`
///    (missing), so `1/2` gives `1/.` and `./1`.
/// 2. TRIMMING, on records with exactly one ALT allele made of A, C, G, T, N
///    (symbolic alleles, breakends, `*` and `.` are left alone): bases shared
///    at the end of REF and ALT are dropped, then bases shared at the start,
///    adding 1 to POS for each -- always leaving one base in each allele.
///    Without splitting, a multiallelic record is not trimmed.
/// 3. LEFT-ALIGNMENT, only when `referencia` (a FASTA: plain, `.gz` or
///    `.lz4`; a `.fai` next to a plain one is used) is given, on the same
///    records as the trimming, SNPs included: REF is checked against the
///    reference, and an indel inside a repeat moves to its leftmost position
///    (see `alinha`) -- what `bcftools norm -f` does. A record whose REF
///    differs from the reference, or whose sequence the reference lacks
///    under either naming (`chr1` / `1`), is only trimmed, and counted.
///    Without a reference, two files that describe the same indel at
///    different positions of a repeat still differ after normalization.
///
/// Trimming and alignment move POS, so the output may stop being sorted. The
/// normalized records go to a temporary file (in `op.dir_temp`, or the system
/// temporary directory); if any POS changed, the output is that file sorted
/// by position with the out-of-core stable sorter, using `op` as given
/// (memory budget and chromosome order); otherwise it is that file as it is,
/// in the input order. Any file size.
///
/// One line is added just before `#CHROM`:
/// `##MMGEN_normalize=<split_multiallelic=true|false>`, with
/// `,reference=NAME` (the file name, without its directory) when aligned.
pub fn vcf_normalize(
    entrada: &str,
    saida: &str,
    separar_multialelicos: bool,
    referencia: Option<&str>,
    op: &Opcoes,
) -> Result<NormSummary, String> {
    recusa_mesmo_arquivo(entrada, saida)?;
    let mut l = Leitor::abre(entrada)?;
    let mut genoma = match referencia {
        Some(c) => Some(Genoma::abre(c, op.dir_temp.as_deref())?),
        None => None,
    };
    let marca = match referencia {
        Some(c) => {
            let nome = std::path::Path::new(c).file_name().map_or_else(|| c.to_string(), |n| n.to_string_lossy().to_string());
            format!("##MMGEN_normalize=<split_multiallelic={separar_multialelicos},reference={nome}>")
        }
        None => format!("##MMGEN_normalize=<split_multiallelic={separar_multialelicos}>"),
    };
    // The temporary file has the compression of the output: when no sorting
    // is needed it becomes the output by a rename.
    let mut temp = DirTemp::novo(op.dir_temp.as_deref())?;
    let provisorio = temp.arquivo(&format!("normalized.vcf{}", compressao(saida)));
    let mut w = Escritor::cria_para(&provisorio, Formato::Vcf)?;

    let mut s = NormSummary::default();
    let mut decl_info: HashMap<String, Numero> = HashMap::new();
    let mut decl_formato: HashMap<String, Numero> = HashMap::new();
    let mut viu_chrom = false;
    let mut acoes: Vec<Acao> = Vec::new();
    let mut cauda = String::new();
    let mut linha = String::new();

    while l.avanca()? {
        let t = l.texto();
        if t.is_empty() {
            continue;
        }
        if t.starts_with('#') {
            if let Some((id, num)) = declaracao(t, "##INFO=<") {
                decl_info.insert(id, num);
            } else if let Some((id, num)) = declaracao(t, "##FORMAT=<") {
                decl_formato.insert(id, num);
            } else if !viu_chrom && t.starts_with("#CHROM") {
                viu_chrom = true;
                w.linha(&marca)?;
            }
            w.linha(t)?;
            continue;
        }
        let Some(r) = registro(t) else {
            s.ignoradas += 1;
            continue;
        };
        s.lidos += 1;
        let n = n_alt(r.alt);
        if n > 1 && separar_multialelicos {
            s.separados += 1;
            for (k, alelo) in r.alt.split(',').enumerate() {
                cauda.clear();
                cauda.push_str(r.qual);
                cauda.push('\t');
                cauda.push_str(r.filter);
                cauda.push('\t');
                info_do_alelo(r.info, &decl_info, k + 1, n, &mut cauda);
                if let Some(resto) = r.resto {
                    cauda.push('\t');
                    amostras_do_alelo(resto, &decl_formato, k + 1, n, &mut acoes, &mut cauda);
                }
                grava_registro(&mut w, &mut s, &mut linha, &r, alelo, &cauda, None, genoma.as_mut())?;
            }
        } else {
            // QUAL onwards, exactly as in the input. Only a single ALT allele
            // is aligned: a multiallelic record that was not split stays.
            let ini = r.chrom.len() + r.pos_txt.len() + r.id.len() + r.ref_.len() + r.alt.len() + 5;
            let g = if n == 1 { genoma.as_mut() } else { None };
            grava_registro(&mut w, &mut s, &mut linha, &r, r.alt, &t[ini..], Some(t), g)?;
        }
    }
    w.termina()?;
    // The input is closed before the output is produced.
    drop(l);

    let mut guarda = Saida { caminho: saida, armada: !std::path::Path::new(saida).exists() };
    if s.pos_alteradas > 0 {
        let r = ordenar_arquivos(std::slice::from_ref(&provisorio), saida, Formato::Vcf, Criterio::Posicao, op)?;
        s.reordenado = true;
        s.lotes = r.lotes;
    } else {
        move_arquivo(&provisorio, saida)?;
    }
    guarda.armada = false;
    Ok(s)
}
