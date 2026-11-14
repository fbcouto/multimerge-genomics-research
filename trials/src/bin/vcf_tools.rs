#![allow(dead_code)]
// Tests of src/vcf_tools.rs: vcf_stats, vcf_filter and vcf_normalize.
//
// Every result is compared with something this file computes on its own:
//  - statistics: the generator knows the type of each allele it builds, and
//    the tallies come from that, not from reading the alleles back;
//  - filter: each expression has a plain Rust closure over the generated
//    record (no second parser);
//  - normalization: an independent split-and-trim over structured records,
//    plus properties checked on the output itself (one ALT per record, minimal
//    alleles, same haplotype against a random reference, sorted output).
// Usage: cargo run --release --bin vcf_tools [-- real_file.vcf.gz ...]
#[path = "../../../src/multimerge_core.rs"]
pub mod multimerge_core;
#[path = "../../../src/bgzf.rs"]
#[allow(dead_code)]
mod bgzf;
#[path = "../../../src/bam.rs"]
#[allow(dead_code)]
mod bam;
#[path = "../../../src/externo.rs"]
mod externo;
#[path = "../../../src/nomes.rs"]
mod nomes;
#[path = "../../../src/referencia.rs"]
mod referencia;
#[path = "../../../src/vcf_tools.rs"]
mod vcf_tools;

use externo::*;
use std::io::Write;
use std::path::PathBuf;
use vcf_tools::*;

// ---------------------------------------------------------------- helpers

fn lcg(s: &mut u64) -> u64 {
    *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    *s >> 11
}

fn escolhe<'a>(s: &mut u64, v: &[&'a str]) -> &'a str {
    v[(lcg(s) % v.len() as u64) as usize]
}

const BASES: [u8; 4] = *b"ACGT";

fn base(s: &mut u64) -> u8 {
    BASES[(lcg(s) % 4) as usize]
}

/// A base different from `b`, whatever the case of `b`.
fn outra(b: u8, s: &mut u64) -> u8 {
    loop {
        let x = base(s);
        if x != b.to_ascii_uppercase() {
            return x;
        }
    }
}

/// Lower case now and then: the tools must ignore the case of the bases.
fn caixa(b: u8, s: &mut u64) -> u8 {
    if lcg(s) % 6 == 0 {
        b.to_ascii_lowercase()
    } else {
        b
    }
}

fn seq(n: usize, s: &mut u64) -> String {
    (0..n).map(|_| caixa(base(s), s) as char).collect()
}

fn muda_caixa(t: &str, s: &mut u64) -> String {
    t.bytes().map(|b| caixa(b, s) as char).collect()
}

fn lista(n: usize, s: &mut u64) -> String {
    (0..n).map(|_| (lcg(s) % 100).to_string()).collect::<Vec<_>>().join(",")
}

struct Amb {
    dir: PathBuf,
    runs: PathBuf,
}

impl Amb {
    fn p(&self, n: &str) -> String {
        self.dir.join(n).to_string_lossy().to_string()
    }
    fn op(&self, mem: usize, ordem: OrdemCromossomos) -> Opcoes {
        Opcoes { memoria_bytes: mem, dir_temp: Some(self.runs.clone()), ordem }
    }
    fn sem_temporarios(&self) -> bool {
        std::fs::read_dir(&self.runs).unwrap().next().is_none()
    }
}

fn le(p: &str) -> Vec<String> {
    let mut r = abre_leitura(p).unwrap();
    let mut s = String::new();
    std::io::Read::read_to_string(&mut r, &mut s).unwrap();
    assert!(!s.contains('\r'), "{p}: the output must have LF line endings only");
    s.lines().map(str::to_string).collect()
}

/// Writes the lines with the compression the extension asks for.
fn grava(p: &str, linhas: &[String], fim: &str) {
    let mut w = Escritor::cria(p).unwrap();
    for l in linhas {
        w.write_all(l.as_bytes()).unwrap();
        w.write_all(fim.as_bytes()).unwrap();
    }
    w.termina().unwrap();
}

fn dados(linhas: &[String]) -> Vec<String> {
    linhas.iter().filter(|l| !l.starts_with('#')).cloned().collect()
}

fn cab(linhas: &[String]) -> Vec<String> {
    linhas.iter().filter(|l| l.starts_with('#')).cloned().collect()
}

/// Shows the first difference between two lists of lines.
fn iguais(got: &[String], esp: &[String], contexto: &str) {
    for (i, (g, e)) in got.iter().zip(esp.iter()).enumerate() {
        assert!(g == e, "{contexto}: first difference at record {i}\n  got:      {g}\n  expected: {e}");
    }
    assert!(got.len() == esp.len(), "{contexto}: got {} lines, expected {}", got.len(), esp.len());
}

// ---------------------------------------------------------------- a small JSON reader

#[derive(Debug, Clone, PartialEq)]
enum J {
    Nulo,
    Bool(bool),
    Num(f64),
    Txt(String),
    Lista(Vec<J>),
    Obj(Vec<(String, J)>),
}

struct LeJson<'a> {
    b: &'a [u8],
    i: usize,
}

impl LeJson<'_> {
    fn brancos(&mut self) {
        while self.i < self.b.len() && self.b[self.i].is_ascii_whitespace() {
            self.i += 1;
        }
    }
    fn espera(&mut self, c: u8) {
        self.brancos();
        assert!(self.i < self.b.len() && self.b[self.i] == c, "JSON: expected `{}` at byte {}", c as char, self.i);
        self.i += 1;
    }
    fn palavra(&mut self, p: &str, v: J) -> J {
        assert!(self.b[self.i..].starts_with(p.as_bytes()), "JSON: bad literal at byte {}", self.i);
        self.i += p.len();
        v
    }
    fn valor(&mut self) -> J {
        self.brancos();
        assert!(self.i < self.b.len(), "JSON: unexpected end");
        match self.b[self.i] {
            b'{' => {
                self.i += 1;
                let mut v = Vec::new();
                self.brancos();
                if self.b[self.i] == b'}' {
                    self.i += 1;
                    return J::Obj(v);
                }
                loop {
                    let k = match self.valor() {
                        J::Txt(k) => k,
                        x => panic!("JSON: key is not a string: {x:?}"),
                    };
                    self.espera(b':');
                    v.push((k, self.valor()));
                    self.brancos();
                    if self.b[self.i] == b'}' {
                        self.i += 1;
                        return J::Obj(v);
                    }
                    self.espera(b',');
                }
            }
            b'[' => {
                self.i += 1;
                let mut v = Vec::new();
                self.brancos();
                if self.b[self.i] == b']' {
                    self.i += 1;
                    return J::Lista(v);
                }
                loop {
                    v.push(self.valor());
                    self.brancos();
                    if self.b[self.i] == b']' {
                        self.i += 1;
                        return J::Lista(v);
                    }
                    self.espera(b',');
                }
            }
            b'"' => {
                self.i += 1;
                let mut o: Vec<u8> = Vec::new();
                loop {
                    let c = self.b[self.i];
                    self.i += 1;
                    match c {
                        b'"' => return J::Txt(String::from_utf8(o).expect("JSON: text is not UTF-8")),
                        b'\\' => {
                            let e = self.b[self.i];
                            self.i += 1;
                            match e {
                                b'n' => o.push(b'\n'),
                                b'r' => o.push(b'\r'),
                                b't' => o.push(b'\t'),
                                b'u' => {
                                    let h = std::str::from_utf8(&self.b[self.i..self.i + 4]).unwrap();
                                    let ch = char::from_u32(u32::from_str_radix(h, 16).unwrap()).unwrap();
                                    o.extend(ch.to_string().bytes());
                                    self.i += 4;
                                }
                                b'"' | b'\\' | b'/' => o.push(e),
                                x => panic!("JSON: bad escape \\{}", x as char),
                            }
                        }
                        c => {
                            assert!(c >= 0x20, "JSON: raw control character in a string");
                            o.push(c);
                        }
                    }
                }
            }
            b't' => self.palavra("true", J::Bool(true)),
            b'f' => self.palavra("false", J::Bool(false)),
            b'n' => self.palavra("null", J::Nulo),
            _ => {
                let ini = self.i;
                while self.i < self.b.len() && b"+-0123456789.eE".contains(&self.b[self.i]) {
                    self.i += 1;
                }
                let t = std::str::from_utf8(&self.b[ini..self.i]).unwrap();
                J::Num(t.parse().unwrap_or_else(|_| panic!("JSON: `{t}` is not a number (byte {ini})")))
            }
        }
    }
}

fn json(s: &str) -> J {
    let mut p = LeJson { b: s.as_bytes(), i: 0 };
    let v = p.valor();
    p.brancos();
    assert!(p.i == s.len(), "JSON: text left after the value: {s}");
    v
}

impl J {
    fn chaves(&self) -> Vec<&str> {
        match self {
            J::Obj(v) => v.iter().map(|(k, _)| k.as_str()).collect(),
            x => panic!("JSON: not an object: {x:?}"),
        }
    }
    fn campo(&self, k: &str) -> &J {
        match self {
            J::Obj(v) => &v.iter().find(|(x, _)| x == k).unwrap_or_else(|| panic!("JSON: missing key `{k}`")).1,
            x => panic!("JSON: not an object: {x:?}"),
        }
    }
    fn num(&self) -> f64 {
        match self {
            J::Num(x) => *x,
            x => panic!("JSON: not a number: {x:?}"),
        }
    }
}

// ---------------------------------------------------------------- records for statistics and filter

/// One ALT allele, with the type it was BUILT to have:
/// 's' SNP, 'm' MNP, 'i' insertion, 'd' deletion, 'o' other.
#[derive(Clone, Debug)]
struct Alelo {
    txt: String,
    tipo: char,
    /// For a SNP: transition (true) or transversion (false); `None` with N.
    ts: Option<bool>,
}

#[derive(Clone, Debug)]
struct Reg {
    chrom: String,
    pos: u64,
    id: String,
    rf: String,
    alts: Vec<Alelo>,
    qual_txt: String,
    qual: Option<f64>,
    filter: String,
    /// (key, value); `None` is a flag.
    info: Vec<(String, Option<String>)>,
    amostras: Vec<String>,
}

impl Reg {
    fn alt_txt(&self) -> String {
        if self.alts.is_empty() {
            ".".to_string()
        } else {
            self.alts.iter().map(|a| a.txt.as_str()).collect::<Vec<_>>().join(",")
        }
    }
    fn info_txt(&self) -> String {
        if self.info.is_empty() {
            return ".".to_string();
        }
        let itens: Vec<String> = self
            .info
            .iter()
            .map(|(k, v)| match v {
                Some(v) => format!("{k}={v}"),
                None => k.clone(),
            })
            .collect();
        itens.join(";")
    }
    fn fixas(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.chrom,
            self.pos,
            self.id,
            self.rf,
            self.alt_txt(),
            self.qual_txt,
            self.filter,
            self.info_txt()
        )
    }
    fn linha(&self) -> String {
        let mut l = self.fixas();
        if !self.amostras.is_empty() {
            l.push_str("\tGT:DP");
            for a in &self.amostras {
                l.push('\t');
                l.push_str(a);
            }
        }
        l
    }
    fn tem(&self, k: &str) -> bool {
        self.info.iter().any(|(x, _)| x == k)
    }
    /// The comma-separated entries of a key that has a value.
    fn itens(&self, k: &str) -> Vec<&str> {
        match self.info.iter().find(|(x, _)| x == k) {
            Some((_, Some(v))) => v.split(',').collect(),
            _ => Vec::new(),
        }
    }
    /// Is there an entry that is a number and satisfies `f`?
    fn algum_num(&self, k: &str, f: impl Fn(f64) -> bool) -> bool {
        self.itens(k).iter().filter_map(|x| x.parse::<f64>().ok()).any(f)
    }
    /// Is there an entry that is not missing and satisfies `f`?
    fn algum_txt(&self, k: &str, f: impl Fn(&str) -> bool) -> bool {
        self.itens(k).iter().filter(|x| **x != ".").any(|x| f(x))
    }
    fn algum(&self, t: char) -> bool {
        self.alts.iter().any(|a| a.tipo == t)
    }
    /// Every ALT allele is of one of the types in `t` (and there is one).
    fn todos(&self, t: &str) -> bool {
        !self.alts.is_empty() && self.alts.iter().all(|a| t.contains(a.tipo))
    }
    fn filtros(&self) -> Vec<&str> {
        if self.filter == "." {
            Vec::new()
        } else {
            self.filter.split(';').collect()
        }
    }
    fn q(&self, f: impl Fn(f64) -> bool) -> bool {
        self.qual.is_some_and(f)
    }
    fn fim(&self) -> u64 {
        self.pos + self.rf.len() as u64 - 1
    }
}

fn purina(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'G')
}

fn gera_alelo(rf: &str, s: &mut u64) -> Alelo {
    let r = rf.as_bytes();
    let sorte = lcg(s) % 100;
    let alelo = |txt: String, tipo: char, ts: Option<bool>| Alelo { txt, tipo, ts };
    if sorte < 8 {
        // symbolic, breakends, `*`, IUPAC codes: not a plain sequence change
        let t = escolhe(s, &["<DEL>", "<INS:ME:ALU>", "G]17:198982]", "]13:123456]T", "*", "R", "AY", "<NON_REF>", ".A"]);
        return alelo(t.to_string(), 'o', None);
    }
    if sorte < 10 {
        // the same bases as REF, in the other case: not a variant
        let t: String = rf
            .chars()
            .map(|c| if c.is_ascii_uppercase() { c.to_ascii_lowercase() } else { c.to_ascii_uppercase() })
            .collect();
        return alelo(t, 'o', None);
    }
    if r.len() == 1 {
        if sorte >= 75 {
            return alelo(seq(2 + (lcg(s) % 4) as usize, s), 'i', None);
        }
        if r[0].eq_ignore_ascii_case(&b'N') {
            return alelo((caixa(base(s), s) as char).to_string(), 's', None);
        }
        if sorte < 13 {
            return alelo("N".to_string(), 's', None);
        }
        let a = outra(r[0], s);
        return alelo((caixa(a, s) as char).to_string(), 's', Some(purina(a) == purina(r[0])));
    }
    match sorte % 3 {
        0 => alelo(r.iter().map(|&b| caixa(outra(b, s), s) as char).collect(), 'm', None),
        1 => alelo(seq(1 + (lcg(s) % (r.len() as u64 - 1)) as usize, s), 'd', None),
        _ => alelo(seq(r.len() + 1 + (lcg(s) % 3) as usize, s), 'i', None),
    }
}

fn gera_reg(s: &mut u64, n_amostras: usize) -> Reg {
    let chrom = escolhe(s, &["chr1", "chr1", "chr2", "chrX"]).to_string();
    let pos = 1 + lcg(s) % 3000;
    let rf = if lcg(s) % 40 == 0 { "N".to_string() } else { seq([1, 1, 1, 1, 2, 3, 4][(lcg(s) % 7) as usize], s) };
    let n_alts = match lcg(s) % 25 {
        0 => 0,
        1..=5 => 2 + (lcg(s) % 3) as usize,
        _ => 1,
    };
    let alts: Vec<Alelo> = (0..n_alts).map(|_| gera_alelo(&rf, s)).collect();
    let (qual_txt, qual) = match lcg(s) % 10 {
        0 => (".".to_string(), None),
        1 => ("50".to_string(), Some(50.0)),
        2..=5 => {
            let q = lcg(s) % 100;
            (q.to_string(), Some(q as f64))
        }
        _ => {
            let t = format!("{}.{}", lcg(s) % 100, lcg(s) % 100);
            let v: f64 = t.parse().unwrap();
            (t, Some(v))
        }
    };
    let filter = escolhe(s, &["PASS", "PASS", "PASS", ".", "q10", "q10;s50"]).to_string();
    let id = if lcg(s) % 3 == 0 { ".".to_string() } else { format!("rs{}", lcg(s) % 40) };

    let mut info: Vec<(String, Option<String>)> = Vec::new();
    let mut poe = |k: &str, v: Option<String>| info.push((k.to_string(), v));
    match lcg(s) % 20 {
        0 => poe("DP", Some(".".to_string())),
        1..=3 => {}
        _ => poe("DP", Some((lcg(s) % 60).to_string())),
    }
    if n_alts > 0 && lcg(s) % 7 != 0 {
        let v: Vec<&str> = (0..n_alts).map(|_| escolhe(s, &["0.01", "0.1", "0.25", "0.5", "0.75"])).collect();
        poe("AF", Some(v.join(",")));
    }
    if lcg(s) % 3 == 0 {
        poe("DB", None);
    }
    if lcg(s) % 2 == 0 {
        let v: Vec<&str> =
            (0..1 + lcg(s) % 3).map(|_| escolhe(s, &["missense", "synonymous", "intron", "stop_gained"])).collect();
        poe("CSQ", Some(v.join(",")));
    }
    if lcg(s) % 5 == 0 {
        let v: Vec<&str> = (0..1 + lcg(s) % 4).map(|_| escolhe(s, &["3", ".", "abc", "9", "5", "12.5"])).collect();
        poe("MX", Some(v.join(",")));
    }
    if lcg(s) % 5 < 2 {
        poe("BQ", Some(escolhe(s, &["-3.5", "-1", "0", "2", "4.25"]).to_string()));
    }
    if lcg(s) % 10 == 0 {
        poe("1000G", None);
    }
    if lcg(s) % 12 == 0 {
        poe("NT", Some("a\"b".to_string()));
    }
    let amostras = (0..n_amostras).map(|_| format!("{}/{}:{}", lcg(s) % 2, lcg(s) % 2, lcg(s) % 50)).collect();
    Reg { chrom, pos, id, rf, alts, qual_txt, qual, filter, info, amostras }
}

/// A record written by hand: one ALT allele of a known type.
fn plantado(chrom: &str, pos: u64, rf: &str, alt: &str, tipo: char, n_amostras: usize) -> Reg {
    Reg {
        chrom: chrom.to_string(),
        pos,
        id: format!("planted{pos}"),
        rf: rf.to_string(),
        alts: vec![Alelo { txt: alt.to_string(), tipo, ts: None }],
        qual_txt: "30".to_string(),
        qual: Some(30.0),
        filter: "PASS".to_string(),
        info: vec![("DP".to_string(), Some("10".to_string()))],
        amostras: (0..n_amostras).map(|i| format!("0/1:{i}")).collect(),
    }
}

fn cabecalho(amostras: &[&str]) -> Vec<String> {
    let mut c: Vec<String> = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=5000>",
        "##contig=<ID=chr2,length=5000>",
        "##contig=<ID=chrX,length=5000>",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
    ]
    .iter()
    .map(|x| x.to_string())
    .collect();
    let mut col = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string();
    if !amostras.is_empty() {
        col.push_str("\tFORMAT");
        for a in amostras {
            col.push('\t');
            col.push_str(a);
        }
    }
    c.push(col);
    c
}

/// Three lines that are not records (and one blank), to be left out.
const TORTAS: [&str; 3] = ["chr1\tabc\t.\tA\tT\t1\tPASS\t.", "chr1\t5\t.\tA", ""];

// ---------------------------------------------------------------- A) statistics

fn stats_esperado(regs: &[Reg], n_amostras: usize, ignoradas: usize) -> VcfStats {
    let mut e = VcfStats { records: regs.len(), samples: n_amostras, skipped: ignoradas, ..VcfStats::default() };
    let (mut soma, mut com_qual) = (0f64, 0usize);
    for r in regs {
        match e.sequences.iter_mut().find(|(n, _)| *n == r.chrom) {
            Some(x) => x.1 += 1,
            None => e.sequences.push((r.chrom.clone(), 1)),
        }
        e.no_alt += r.alts.is_empty() as usize;
        e.multiallelic += (r.alts.len() > 1) as usize;
        e.snps += r.algum('s') as usize;
        e.mnps += r.algum('m') as usize;
        e.insertions += r.algum('i') as usize;
        e.deletions += r.algum('d') as usize;
        e.other += r.algum('o') as usize;
        for a in &r.alts {
            match a.ts {
                Some(true) => e.transitions += 1,
                Some(false) => e.transversions += 1,
                None => {}
            }
        }
        e.pass += (r.filter == "PASS") as usize;
        e.filter_missing += (r.filter == ".") as usize;
        match r.qual {
            Some(q) => {
                soma += q;
                com_qual += 1;
            }
            None => e.qual_missing += 1,
        }
    }
    if e.transversions > 0 {
        e.ts_tv = e.transitions as f64 / e.transversions as f64;
    }
    if com_qual > 0 {
        e.mean_qual = soma / com_qual as f64;
    }
    e
}

fn confere_stats(got: &VcfStats, esp: &VcfStats, contexto: &str) {
    assert!((got.ts_tv - esp.ts_tv).abs() < 1e-12, "{contexto}: ts_tv {} != {}", got.ts_tv, esp.ts_tv);
    assert!((got.mean_qual - esp.mean_qual).abs() < 1e-9, "{contexto}: mean_qual {} != {}", got.mean_qual, esp.mean_qual);
    let mut g = got.clone();
    g.ts_tv = esp.ts_tv;
    g.mean_qual = esp.mean_qual;
    assert_eq!(&g, esp, "{contexto}");
    // the JSON says the same, with the keys in English
    let j = json(&got.json());
    assert_eq!(
        j.chaves(),
        [
            "records", "samples", "sequences", "snps", "mnps", "insertions", "deletions", "other", "no_alt",
            "multiallelic", "transitions", "transversions", "ts_tv", "pass", "filter_missing", "qual_missing",
            "mean_qual", "skipped"
        ],
        "{contexto}"
    );
    for (k, v) in [
        ("records", got.records),
        ("samples", got.samples),
        ("snps", got.snps),
        ("mnps", got.mnps),
        ("insertions", got.insertions),
        ("deletions", got.deletions),
        ("other", got.other),
        ("no_alt", got.no_alt),
        ("multiallelic", got.multiallelic),
        ("transitions", got.transitions),
        ("transversions", got.transversions),
        ("pass", got.pass),
        ("filter_missing", got.filter_missing),
        ("qual_missing", got.qual_missing),
        ("skipped", got.skipped),
    ] {
        assert_eq!(j.campo(k).num(), v as f64, "{contexto}: JSON key {k}");
    }
    assert!((j.campo("ts_tv").num() - got.ts_tv).abs() < 1e-12 && (j.campo("mean_qual").num() - got.mean_qual).abs() < 1e-9);
    let seqs: Vec<J> =
        got.sequences.iter().map(|(n, q)| J::Lista(vec![J::Txt(n.clone()), J::Num(*q as f64)])).collect();
    assert_eq!(j.campo("sequences"), &J::Lista(seqs), "{contexto}");
    let texto = got.resumo();
    assert!(texto.starts_with(&format!("{} record(s), {} sample(s)", got.records, got.samples)), "{texto}");
}

fn testa_stats(amb: &Amb, s: &mut u64) -> usize {
    let mut ensaios = 0;

    // random files: every allele type, several chromosomes, with and without samples
    for (n, n_amostras, nome, fim) in [
        (0usize, 0usize, "st0.vcf", "\n"),
        (1, 3, "st1.vcf", "\n"),
        (5000, 3, "st2.vcf", "\r\n"),
        (20_000, 0, "st3.vcf.gz", "\n"),
        (3000, 2, "st4.vcf.lz4", "\r\n"),
    ] {
        let nomes_amostras: Vec<String> = (0..n_amostras).map(|i| format!("S{}", i + 1)).collect();
        let regs: Vec<Reg> = (0..n).map(|_| gera_reg(s, n_amostras)).collect();
        let mut linhas = cabecalho(&nomes_amostras.iter().map(String::as_str).collect::<Vec<_>>());
        let mut ignoradas = 0;
        for (i, r) in regs.iter().enumerate() {
            linhas.push(r.linha());
            if i % 400 == 7 {
                linhas.extend(TORTAS.iter().map(|x| x.to_string()));
                ignoradas += 2;
            }
        }
        grava(&amb.p(nome), &linhas, fim);
        let got = vcf_stats(&amb.p(nome)).unwrap();
        confere_stats(&got, &stats_esperado(&regs, n_amostras, ignoradas), nome);
        if n >= 3000 {
            // the random data really covers every category
            assert!(
                got.snps > 0 && got.mnps > 0 && got.insertions > 0 && got.deletions > 0 && got.other > 0,
                "{nome}: {}",
                got.resumo()
            );
            assert!(got.no_alt > 0 && got.multiallelic > 0 && got.transitions > 0 && got.transversions > 0);
            assert!(got.qual_missing > 0 && got.filter_missing > 0 && got.sequences.len() == 3);
        }
        ensaios += 1;
    }

    // the decisions the specification left open, one per line, counted by hand
    let fixo = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB",
        "chr2\t1\t.\tA\tA\t10\tPASS\t.",           // ALT equal to REF: other
        "chr2\t2\t.\tAC\tac\t20\tPASS\t.",         // the same, ignoring case: other
        "chr2\t3\t.\tN\tA\t.\t.\t.",               // SNP from N: a SNP, neither Ts nor Tv
        "chr2\t4\t.\tA\tN\t.\tq10\t.",             // SNP to N: the same
        "chr2\t5\t.\tR\tA\t30\tPASS\t.",           // REF that is not ACGTN: other
        "chr1\t6\t.\tA\tG,T,<DEL>,AT,*\t40\tPASS\t.", // SNP (1 Ts + 1 Tv), other, insertion; multiallelic
        "chr1\t7\t.\tAT\tA,ATT,GC\t50\tPASS\t.",   // deletion, insertion, MNP; multiallelic
        "chr1\t8\t.\ta\tg\tabc\tPASS\t.",          // lower case: SNP, transition; QUAL neither number nor `.`
        "chr1\t9\t.\t.\t.\tinf\tPASS\t.",          // no ALT; "inf" is not a QUAL
        "chr2\t10\t.\tC\tT,T\t60\ts50\t.",         // two SNP alleles on one site: 1 SNP site, 2 transitions
        "chr3\t11\t.\tAC\tAT\t70\tPASS\t.",        // same length above one: MNP (no trimming first)
        "chr3\t12\t.\tACG\tAT\t80\tPASS\t.",       // ALT shorter: deletion
    ];
    std::fs::write(amb.p("fixo.vcf"), fixo.join("\n")).unwrap(); // no line ending on the last line
    let got = vcf_stats(&amb.p("fixo.vcf")).unwrap();
    let esp = VcfStats {
        records: 12,
        samples: 2,
        sequences: vec![("chr2".to_string(), 6), ("chr1".to_string(), 4), ("chr3".to_string(), 2)],
        snps: 5,
        mnps: 2,
        insertions: 2,
        deletions: 2,
        other: 4,
        no_alt: 1,
        multiallelic: 3,
        transitions: 4,
        transversions: 1,
        ts_tv: 4.0,
        pass: 9,
        filter_missing: 1,
        qual_missing: 2,
        mean_qual: 45.0,
        skipped: 0,
    };
    confere_stats(&got, &esp, "fixo.vcf");
    ensaios += 1;

    // empty file, and errors
    std::fs::write(amb.p("vazio.vcf"), "").unwrap();
    confere_stats(&vcf_stats(&amb.p("vazio.vcf")).unwrap(), &VcfStats::default(), "empty file");
    let e = vcf_stats(&amb.p("does_not_exist.vcf")).unwrap_err();
    assert!(e.contains("does_not_exist.vcf"), "{e}");
    ensaios += 1;
    ensaios
}

// ---------------------------------------------------------------- B) filter

struct Caso {
    expr: &'static str,
    f: Box<dyn Fn(&Reg) -> bool>,
    /// The expected selection may be empty or everything.
    livre: bool,
}

fn caso(expr: &'static str, f: impl Fn(&Reg) -> bool + 'static) -> Caso {
    Caso { expr, f: Box::new(f), livre: false }
}

fn caso_livre(expr: &'static str, f: impl Fn(&Reg) -> bool + 'static) -> Caso {
    Caso { expr, f: Box::new(f), livre: true }
}

/// Every feature of the language, each with the same condition written
/// directly in Rust over the generated record.
fn casos() -> Vec<Caso> {
    vec![
        // QUAL: every operator, the literal on either side, missing values
        caso("QUAL>30", |r| r.q(|q| q > 30.0)),
        caso("QUAL>=30 && QUAL<=60", |r| r.q(|q| (30.0..=60.0).contains(&q))),
        caso("QUAL<10.5", |r| r.q(|q| q < 10.5)),
        caso("QUAL==50", |r| r.q(|q| q == 50.0)),
        caso("QUAL = 50.0", |r| r.q(|q| q == 50.0)),
        caso("QUAL!=50", |r| r.q(|q| q != 50.0)), // a missing QUAL is NOT "different from 50"
        caso("30<QUAL", |r| r.q(|q| q > 30.0)),
        caso("60 >= QUAL", |r| r.q(|q| q <= 60.0)),
        caso("QUAL==\".\"", |r| r.qual.is_none()),
        caso("QUAL!=\".\"", |r| r.qual.is_some()),
        caso("!(QUAL!=50)", |r| !r.q(|q| q != 50.0)), // true for 50 AND for a missing QUAL
        // POS
        caso("POS>=1000 & POS<2000", |r| r.pos >= 1000 && r.pos < 2000),
        caso("POS<=150", |r| r.pos <= 150),
        caso("POS != 150 && POS < 152 && POS > 147", |r| r.pos != 150 && r.pos < 152 && r.pos > 147),
        // text columns: exact, case-sensitive, whole column
        caso("CHROM==\"chr2\"", |r| r.chrom == "chr2"),
        caso("CHROM!=\"chr2\" && CHROM!=\"chrX\"", |r| r.chrom != "chr2" && r.chrom != "chrX"),
        caso("\"chrX\" == CHROM", |r| r.chrom == "chrX"),
        caso("ID==\".\"", |r| r.id == "."),
        caso("ID!=\".\"", |r| r.id != "."),
        caso("ID==\"rs7\"", |r| r.id == "rs7"),
        caso("ID!=\"rs7\"", |r| r.id != "." && r.id != "rs7"),
        caso("REF==\"A\"", |r| r.rf == "A"),
        caso("REF!=\"A\"", |r| r.rf != "A"),
        caso("REF==\"a\"", |r| r.rf == "a"),
        caso("ALT==\"T\"", |r| r.alt_txt() == "T"),
        caso("ALT!=\"T\"", |r| !r.alts.is_empty() && r.alt_txt() != "T"),
        caso("ALT==\".\"", |r| r.alts.is_empty()),
        caso("FILTER==\"PASS\"", |r| r.filter == "PASS"),
        caso("FILTER!=\"PASS\"", |r| r.filter != "." && r.filter != "PASS"),
        caso("FILTER=\"q10;s50\"", |r| r.filter == "q10;s50"),
        caso("FILTER==\".\"", |r| r.filter == "."),
        // FILTER ~: the filters named are among those of the record (bcftools)
        caso("FILTER~\"q10\"", |r| r.filtros().contains(&"q10")),
        caso("FILTER~\"s50\"", |r| r.filtros().contains(&"s50")),
        caso("FILTER!~\"q10\"", |r| !r.filtros().contains(&"q10")), // a missing FILTER passes
        caso("FILTER~\"s50;q10\"", |r| r.filtros().contains(&"q10") && r.filtros().contains(&"s50")),
        caso("FILTER~\"PASS\" || FILTER==\".\"", |r| r.filter == "PASS" || r.filter == "."),
        // N_ALT
        caso("N_ALT>1", |r| r.alts.len() > 1),
        caso("N_ALT==0", |r| r.alts.is_empty()),
        caso("N_ALT=1 | N_ALT=3", |r| r.alts.len() == 1 || r.alts.len() == 3),
        // TYPE, as in bcftools: == every ALT allele is of the type, ~ at least
        // one; != and !~ are the negations; no ALT is "ref" and nothing else
        caso("TYPE==\"snp\"", |r| r.todos("s")),
        caso("TYPE=\"mnp\"", |r| r.todos("m")),
        caso("TYPE==\"ins\"", |r| r.todos("i")),
        caso("TYPE==\"del\"", |r| r.todos("d")),
        caso("TYPE==\"indel\"", |r| r.todos("id")),
        caso("TYPE==\"other\"", |r| r.todos("o")),
        caso("TYPE!=\"snp\"", |r| !r.todos("s")),
        caso("TYPE~\"snp\"", |r| r.algum('s')),
        caso("TYPE~\"mnp\"", |r| r.algum('m')),
        caso("TYPE~\"ins\"", |r| r.algum('i')),
        caso("TYPE~\"del\"", |r| r.algum('d')),
        caso("TYPE~\"indel\"", |r| r.algum('i') || r.algum('d')),
        caso("TYPE~\"other\"", |r| r.algum('o')),
        caso("TYPE!~\"snp\"", |r| !r.algum('s')),
        caso("!(TYPE~\"snp\")", |r| !r.algum('s')),
        caso("TYPE~\"snp\" && TYPE~\"indel\"", |r| r.algum('s') && (r.algum('i') || r.algum('d'))),
        caso("TYPE~\"snp\" && TYPE!=\"snp\"", |r| r.algum('s') && !r.todos("s")),
        caso("TYPE==\"ref\"", |r| r.alts.is_empty()),
        caso("TYPE!=\"ref\"", |r| !r.alts.is_empty()),
        caso("TYPE~\"ref\"", |r| r.alts.is_empty()),
        caso("TYPE==\".\"", |r| r.alts.is_empty()),
        caso("\"indel\" ~ TYPE", |r| r.algum('i') || r.algum('d')),
        // INFO: numbers, several values, missing, presence
        caso("INFO/DP>20", |r| r.algum_num("DP", |x| x > 20.0)),
        caso("INFO/DP<=20", |r| r.algum_num("DP", |x| x <= 20.0)),
        caso("INFO/DP!=20", |r| r.algum_num("DP", |x| x != 20.0)),
        caso("INFO/DP==\".\"", |r| !r.tem("DP") || r.itens("DP") == ["."]),
        caso("INFO/DP!=\".\"", |r| r.tem("DP") && r.itens("DP") != ["."]),
        caso("INFO/DP", |r| r.tem("DP")),
        caso("INFO/AF>0.5", |r| r.algum_num("AF", |x| x > 0.5)),
        caso("INFO/AF!=0.5", |r| r.algum_num("AF", |x| x != 0.5)),
        caso("INFO/AF==0.25", |r| r.algum_num("AF", |x| x == 0.25)),
        caso("INFO/AF<1e-1", |r| r.algum_num("AF", |x| x < 0.1)),
        caso("INFO/AF >= .5", |r| r.algum_num("AF", |x| x >= 0.5)),
        caso("0.1 > INFO/AF", |r| r.algum_num("AF", |x| x < 0.1)),
        caso("INFO/CSQ==\"missense\"", |r| r.algum_txt("CSQ", |x| x == "missense")),
        caso("INFO/CSQ!=\"missense\"", |r| r.algum_txt("CSQ", |x| x != "missense")),
        caso("INFO/DB", |r| r.tem("DB")),
        caso("!INFO/DB", |r| !r.tem("DB")),
        caso("!!INFO/DB", |r| r.tem("DB")),
        caso_livre("INFO/DB==1", |_| false), // a flag has no value to compare
        caso("INFO/1000G && QUAL>20", |r| r.tem("1000G") && r.q(|q| q > 20.0)),
        caso("INFO/MX>5", |r| r.algum_num("MX", |x| x > 5.0)), // `.` and `abc` satisfy nothing
        caso("INFO/MX!=3", |r| r.algum_num("MX", |x| x != 3.0)),
        caso("INFO/MX==\"abc\"", |r| r.algum_txt("MX", |x| x == "abc")),
        caso("INFO/MX!=\"abc\"", |r| r.algum_txt("MX", |x| x != "abc")),
        caso("INFO/BQ>-1.5", |r| r.algum_num("BQ", |x| x > -1.5)),
        caso("INFO/BQ>=+2", |r| r.algum_num("BQ", |x| x >= 2.0)),
        caso("INFO/BQ<-1", |r| r.algum_num("BQ", |x| x < -1.0)),
        caso("INFO/NT==\"a\\\"b\"", |r| r.algum_txt("NT", |x| x == "a\"b")),
        caso_livre("INFO/NOPE>1 || INFO/NOPE", |_| false),
        // precedence: ! before &&, && before ||; parentheses
        caso("QUAL>50 || POS<500 && CHROM==\"chr1\"", |r| r.q(|q| q > 50.0) || (r.pos < 500 && r.chrom == "chr1")),
        caso("(QUAL>50 || POS<500) && CHROM==\"chr1\"", |r| (r.q(|q| q > 50.0) || r.pos < 500) && r.chrom == "chr1"),
        caso("POS<500 && CHROM==\"chr1\" || QUAL>50", |r| (r.pos < 500 && r.chrom == "chr1") || r.q(|q| q > 50.0)),
        caso("!QUAL>50 && POS<500", |r| !r.q(|q| q > 50.0) && r.pos < 500),
        caso("!(QUAL>50 && POS<500)", |r| !(r.q(|q| q > 50.0) && r.pos < 500)),
        caso("!(QUAL>50 || INFO/DB)", |r| !(r.q(|q| q > 50.0) || r.tem("DB"))),
        caso("QUAL>50 | INFO/DB & !TYPE=\"snp\"", |r| r.q(|q| q > 50.0) || (r.tem("DB") && !r.todos("s"))),
        caso("  QUAL>30&&POS<1000  ", |r| r.q(|q| q > 30.0) && r.pos < 1000),
        caso("((QUAL > 30))", |r| r.q(|q| q > 30.0)),
        caso("QUAL>30 && (FILTER==\"PASS\" || FILTER==\"q10\") && !(INFO/DP<10)", |r| {
            r.q(|q| q > 30.0) && (r.filter == "PASS" || r.filter == "q10") && !r.algum_num("DP", |x| x < 10.0)
        }),
    ]
}

/// (expression, character position the error must name)
const ERRADAS: [(&str, usize); 36] = [
    ("QUAL >", 7),                  // nothing after the operator
    ("QUAL > > 3", 8),              // an operator where a value should be
    ("(QUAL>3", 8),                 // `)` missing
    ("QUAL>3)", 7),                 // `)` without `(`
    ("FOO>3", 1),                   // unknown field
    ("CHROM<\"chr2\"", 6),          // text only takes == and !=
    ("QUAL==\"high\"", 7),          // a number compared with text
    ("POS>1 && \"abc", 10),         // text never closed
    ("QUAL>3 &&", 10),              // nothing after &&
    ("QUAL>3 POS<5", 8),            // && or || missing
    ("INFO/ >3", 6),                // INFO/ without a key
    ("QUAL # 3", 6),                // unexpected character
    ("TYPE==\"snv\"", 7),           // unknown variant type
    ("POS", 4),                     // a field alone
    ("QUAL>POS", 5),                // two fields
    ("3>2", 1),                     // two literals
    ("QUAL ~ 3", 6),                // ~ only with TYPE and FILTER
    ("CHROM!~\"chr\"", 6),          // no regular expressions
    ("CHROM==1", 8),                // text field compared with a number
    ("INFO/DP>\"x\"", 8),           // ordering text
    ("qual>3", 1),                  // field names are upper case
    ("QUAL>3 && ()", 12),           // empty parentheses
    ("!", 2),
    ("ID=='x'", 5),                 // single quotes
    ("QUAL>30abc", 6),              // not a number
    ("QUAL>1e999", 6),              // not a finite number
    ("TYPE<\"snp\"", 5),
    ("ID==\"é\" && FOO>1", 12),     // the position counts characters, not bytes
    ("POS>=", 6),
    ("QUAL>3 || || POS<2", 11),
    ("\"PASS\"", 1),                // a literal alone
    ("QUAL<\".\"", 5),              // a missing value only takes == and !=
    ("N_ALT==\"two\"", 8),
    ("DP>10", 1),                   // an INFO key without INFO/
    ("FMT/DP>10", 1),               // per-sample fields
    ("TYPE==5", 7),
];

fn testa_filtro(amb: &Amb, s: &mut u64) -> usize {
    let mut ensaios = 0;
    let amostras = ["S1", "S2", "S3"];

    // ---- the data: random records plus a few planted by hand for the regions
    let mut regs: Vec<Reg> = (0..4000).map(|_| gera_reg(s, 3)).collect();
    let longo = "ACGTACGTACGTACGTACGTACGTACGTAC"; // 30 bases
    for (chrom, pos, rf, alt, tipo) in [
        ("chr1", 80, longo, "A", 'd'),   // spans 80..109: reaches chr1:100-200 through its REF
        ("chr1", 60, longo, "A", 'd'),   // spans 60..89: does not
        ("chr1", 150, "A", "G", 's'),
        ("chr1", 148, "ACGT", "A", 'd'), // spans 148..151: covers position 150
        ("chr1", 151, "A", "G", 's'),
        ("chr1", 200, "AC", "A", 'd'),   // starts on the last base of chr1:100-200
        ("chr1", 201, "A", "G", 's'),
        ("chr2", 999, "AC", "A", 'd'),   // spans 999..1000
        ("chr2", 2001, "A", "G", 's'),
    ] {
        regs.push(plantado(chrom, pos, rf, alt, tipo, 3));
    }
    let mut linhas = cabecalho(&amostras);
    let mut ignoradas = 0;
    for (i, r) in regs.iter().enumerate() {
        linhas.push(r.linha());
        if i % 500 == 3 {
            linhas.extend(TORTAS.iter().map(|x| x.to_string()));
            ignoradas += 2;
        }
    }
    let entrada = amb.p("f.vcf");
    grava(&entrada, &linhas, "\n");
    let cab_entrada = cabecalho(&amostras);
    let saida = amb.p("f_out.vcf");

    // ---- expressions against the closures
    let lista = casos();
    assert!(lista.len() >= 20);
    for c in &lista {
        let r = vcf_filter(&entrada, &saida, Some(c.expr), &[], None)
            .unwrap_or_else(|e| panic!("expression `{}` was rejected:\n{e}", c.expr));
        let esp: Vec<String> = regs.iter().filter(|r| (c.f)(r)).map(Reg::linha).collect();
        if !c.livre {
            assert!(
                !esp.is_empty() && esp.len() < regs.len(),
                "`{}` selects {} of {} records: the test data does not exercise it",
                c.expr,
                esp.len(),
                regs.len()
            );
        }
        let got = le(&saida);
        iguais(&dados(&got), &esp, &format!("expression `{}`", c.expr));
        assert_eq!(
            r,
            FilterSummary { lidos: regs.len(), mantidos: esp.len(), amostras_antes: 3, amostras_depois: 3, ignoradas },
            "`{}`",
            c.expr
        );
        ensaios += 1;
    }

    // ---- the header: copied, with one line before #CHROM saying what was done
    {
        let expr = "QUAL>30 && FILTER==\"PASS\"";
        vcf_filter(&entrada, &saida, Some(expr), &[], None).unwrap();
        let got = cab(&le(&saida));
        let mut esp = cab_entrada.clone();
        esp.insert(esp.len() - 1, "##MMGEN_filter=<expression=\"QUAL>30 && FILTER==\\\"PASS\\\"\">".to_string());
        iguais(&got, &esp, "header after an expression");
        // a blank expression is no expression: everything is kept
        let r = vcf_filter(&entrada, &saida, Some("   "), &[], None).unwrap();
        assert_eq!((r.lidos, r.mantidos), (regs.len(), regs.len()));
        let got = le(&saida);
        assert_eq!(got[cab_entrada.len() - 1], "##MMGEN_filter=<>");
        iguais(&dados(&got), &regs.iter().map(Reg::linha).collect::<Vec<_>>(), "no expression");
        let j = json(&r.json());
        assert_eq!(j.chaves(), ["read", "kept", "samples_before", "samples_after", "skipped"]);
        assert_eq!(
            (j.campo("read").num(), j.campo("kept").num(), j.campo("samples_before").num(), j.campo("skipped").num()),
            (regs.len() as f64, regs.len() as f64, 3.0, ignoradas as f64)
        );
        ensaios += 1;
    }

    // ---- malformed expressions: an error that names the position and shows the expression
    assert!(ERRADAS.len() >= 8);
    let nao_criada = amb.p("never_created.vcf");
    for (expr, pos) in ERRADAS {
        let e = match vcf_filter(&entrada, &nao_criada, Some(expr), &[], None) {
            Err(e) => e,
            Ok(_) => panic!("`{expr}` should be a syntax error"),
        };
        assert!(e.contains(&format!("at position {pos}:")), "`{expr}`: expected position {pos} in:\n{e}");
        assert!(e.contains(&format!("\n  {expr}\n  {}^\n", " ".repeat(pos - 1))), "`{expr}`: expression and mark missing in:\n{e}");
        assert!(Expressao::compila(expr).is_err());
        assert!(!std::path::Path::new(&nao_criada).exists(), "`{expr}`: an output was created");
        ensaios += 1;
    }
    // hints that make the message useful
    let dica = |expr: &str, trecho: &str| {
        let e = Expressao::compila(expr).unwrap_err();
        assert!(e.contains(trecho), "`{expr}`: expected `{trecho}` in:\n{e}");
    };
    dica("qual>3", "write QUAL");
    dica("DP>10", "INFO/DP");
    dica("TYPE==\"snv\"", "\"indel\"");
    dica("TYPE<\"snp\"", "~ (at least one)");
    dica("CHROM~\"chr\"", "regular expressions are not supported");
    dica("FILTER~\".\"", "== or !=");
    dica("CHROM==1", "CHROM==\"1\"");
    dica("(QUAL>3", "opened at position 1");
    dica("QUAL>3)", "no matching `(`");
    dica("", "empty");
    // nesting and length are bounded, with an error instead of a stack overflow
    dica(&format!("{}QUAL>1{}", "(".repeat(300), ")".repeat(300)), "nested");
    dica(&format!("{}INFO/DB", "!".repeat(300)), "nested");
    dica(&vec!["QUAL>=0"; 700].join(" && "), "too long");
    let comprida = vec!["QUAL>=0"; 400].join(" && ");
    let r = vcf_filter(&entrada, &saida, Some(&comprida), &[], None).unwrap();
    assert_eq!(r.mantidos, regs.iter().filter(|r| r.qual.is_some()).count());
    ensaios += 1;

    // ---- regions
    type Pred = Box<dyn Fn(&Reg) -> bool>;
    let em = |c: &'static str, a: u64, b: u64| -> Pred { Box::new(move |r: &Reg| r.chrom == c && r.pos <= b && r.fim() >= a) };
    let regioes: Vec<(Vec<&str>, Pred)> = vec![
        (vec!["chr1"], Box::new(|r: &Reg| r.chrom == "chr1") as Pred),
        (vec!["chr1:100-200"], em("chr1", 100, 200)),
        (vec!["chr1:150"], em("chr1", 150, 150)),
        (vec!["chr2:1,000-2,000"], em("chr2", 1000, 2000)),
        (vec!["  chrX:2,500-2600 "], em("chrX", 2500, 2600)),
        (vec!["{chr1}:100-200"], em("chr1", 100, 200)),
        (vec!["{chrX}"], Box::new(|r: &Reg| r.chrom == "chrX") as Pred),
        (vec!["chr1:100-200", "chrX", "chr2:1000-2000"], {
            let (a, b) = (em("chr1", 100, 200), em("chr2", 1000, 2000));
            Box::new(move |r: &Reg| a(r) || r.chrom == "chrX" || b(r)) as Pred
        }),
        // overlapping, touching and repeated regions are one union
        (vec!["chr1:300-400", "chr1:100-200", "chr1:150-310", "chr1:401-500", "chr1:100-200"], em("chr1", 100, 500)),
        // regions far apart on one sequence: the gap is not selected
        (vec!["chr1:2000-2100", "chr1:10-20", "chr1:900-950"], {
            let (a, b, c) = (em("chr1", 10, 20), em("chr1", 900, 950), em("chr1", 2000, 2100));
            Box::new(move |r: &Reg| a(r) || b(r) || c(r)) as Pred
        }),
    ];
    for (textos, f) in &regioes {
        let lista: Vec<String> = textos.iter().map(|x| x.to_string()).collect();
        let r = vcf_filter(&entrada, &saida, None, &lista, None).unwrap();
        let esp: Vec<String> = regs.iter().filter(|r| f(r)).map(Reg::linha).collect();
        assert!(!esp.is_empty() && esp.len() < regs.len(), "{textos:?} selects {} records", esp.len());
        iguais(&dados(&le(&saida)), &esp, &format!("regions {textos:?}"));
        assert_eq!((r.lidos, r.mantidos), (regs.len(), esp.len()));
        ensaios += 1;
    }
    {
        // the overlap through a long REF, spelled out
        vcf_filter(&entrada, &saida, None, &["chr1:100-200".to_string()], None).unwrap();
        let got = dados(&le(&saida));
        let tem = |id: &str| got.iter().any(|l| l.split('\t').nth(2) == Some(id));
        assert!(tem("planted80"), "POS 80 with a 30-base REF reaches position 109");
        assert!(tem("planted200") && !tem("planted60") && !tem("planted201"));
        vcf_filter(&entrada, &saida, None, &["chr1:150".to_string()], None).unwrap();
        let got = dados(&le(&saida));
        let tem = |id: &str| got.iter().any(|l| l.split('\t').nth(2) == Some(id));
        assert!(tem("planted150") && tem("planted148") && !tem("planted151"));
        // a sequence that is not in the file selects nothing, without error
        let r = vcf_filter(&entrada, &saida, None, &["chrNope:1-100".to_string()], None).unwrap();
        assert_eq!(r.mantidos, 0);
        assert_eq!(le(&saida)[cab_entrada.len() - 1], "##MMGEN_filter=<regions=\"chrNope:1-100\">");
        // expression AND region
        let r = vcf_filter(&entrada, &saida, Some("TYPE==\"snp\""), &["chr2:1,000-2,000".to_string(), "chr1".to_string()], None)
            .unwrap();
        let esp: Vec<String> = regs
            .iter()
            .filter(|r| r.todos("s") && (r.chrom == "chr1" || (r.chrom == "chr2" && r.pos <= 2000 && r.fim() >= 1000)))
            .map(Reg::linha)
            .collect();
        let got = le(&saida);
        iguais(&dados(&got), &esp, "expression and regions");
        assert_eq!(r.mantidos, esp.len());
        assert_eq!(
            got[cab_entrada.len() - 1],
            "##MMGEN_filter=<expression=\"TYPE==\\\"snp\\\"\";regions=\"chr2:1000-2000,chr1\">"
        );
        ensaios += 1;
    }
    for ruim in ["", "chr1:", "chr1:abc", "chr1:200-100", "chr1:100-", "chr1:-100", ":100", "chr1:1-2-3", "{chr1", "{chr1}x"] {
        let e = match vcf_filter(&entrada, &nao_criada, None, &[ruim.to_string()], None) {
            Err(e) => e,
            Ok(_) => panic!("region `{ruim}` should be an error"),
        };
        assert!(e.contains(&format!("`{ruim}`")) && e.contains("Accepted forms") && e.contains("chr1:1000-2000"), "{e}");
        assert!(!std::path::Path::new(&nao_criada).exists());
    }
    ensaios += 1;

    // ---- samples: subset, order, sites only
    let nomes = |v: &[&str]| -> Vec<String> { v.iter().map(|x| x.to_string()).collect() };
    {
        let pedidas = nomes(&["S3", "S1"]);
        let r = vcf_filter(&entrada, &saida, None, &[], Some(&pedidas)).unwrap();
        let got = le(&saida);
        let mut esp_cab = cabecalho(&["S3", "S1"]);
        esp_cab.insert(esp_cab.len() - 1, "##MMGEN_filter=<samples=\"S3,S1\">".to_string());
        iguais(&cab(&got), &esp_cab, "header with two samples");
        let esp: Vec<String> =
            regs.iter().map(|r| format!("{}\tGT:DP\t{}\t{}", r.fixas(), r.amostras[2], r.amostras[0])).collect();
        iguais(&dados(&got), &esp, "samples S3,S1");
        assert_eq!(
            r,
            FilterSummary { lidos: regs.len(), mantidos: regs.len(), amostras_antes: 3, amostras_depois: 2, ignoradas }
        );
        ensaios += 1;

        // sites only: no FORMAT, no samples
        let r = vcf_filter(&entrada, &saida, None, &[], Some(&[])).unwrap();
        let got = le(&saida);
        let mut esp_cab = cabecalho(&[]);
        esp_cab.insert(esp_cab.len() - 1, "##MMGEN_filter=<samples=\"\">".to_string());
        iguais(&cab(&got), &esp_cab, "header of sites only");
        iguais(&dados(&got), &regs.iter().map(Reg::fixas).collect::<Vec<_>>(), "sites only");
        assert_eq!((r.amostras_antes, r.amostras_depois), (3, 0));
        ensaios += 1;

        // everything at once, compressed in and out, CRLF in
        let (gz, lz4) = (amb.p("f_crlf.vcf.gz"), amb.p("f_out.vcf.lz4"));
        grava(&gz, &linhas, "\r\n");
        let pedidas = nomes(&["S2"]);
        let r = vcf_filter(&gz, &lz4, Some("QUAL>=20 && !INFO/DB"), &["chr1".to_string()], Some(&pedidas)).unwrap();
        let esp: Vec<String> = regs
            .iter()
            .filter(|r| r.q(|q| q >= 20.0) && !r.tem("DB") && r.chrom == "chr1")
            .map(|r| format!("{}\tGT:DP\t{}", r.fixas(), r.amostras[1]))
            .collect();
        let got = le(&lz4);
        iguais(&dados(&got), &esp, "expression, region and sample; gz in, lz4 out");
        assert_eq!(
            got[cab_entrada.len() - 1],
            "##MMGEN_filter=<expression=\"QUAL>=20 && !INFO/DB\";regions=\"chr1\";samples=\"S2\">"
        );
        assert_eq!(got[cab_entrada.len()], "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2");
        assert_eq!(
            r,
            FilterSummary { lidos: regs.len(), mantidos: esp.len(), amostras_antes: 3, amostras_depois: 1, ignoradas }
        );
        ensaios += 1;
    }
    {
        // unknown sample: an error that lists the names there are, and no output left behind
        let e = vcf_filter(&entrada, &nao_criada, None, &[], Some(&nomes(&["S1", "nope", "s2"]))).unwrap_err();
        assert!(e.contains("`nope`") && e.contains("`s2`") && e.contains("S1, S2, S3") && e.contains("f.vcf"), "{e}");
        assert!(!std::path::Path::new(&nao_criada).exists(), "the half-written output was left behind");
        let e = vcf_filter(&entrada, &nao_criada, None, &[], Some(&nomes(&["S1", "S1"]))).unwrap_err();
        assert!(e.contains("`S1`") && e.contains("more than once"), "{e}");
        // an output that already existed is not deleted by a failed run
        std::fs::write(amb.p("ja_existia.vcf"), "x\n").unwrap();
        assert!(vcf_filter(&entrada, &amb.p("ja_existia.vcf"), None, &[], Some(&nomes(&["nope"]))).is_err());
        assert!(std::path::Path::new(&amb.p("ja_existia.vcf")).exists());
        // a file with no samples, and one with no #CHROM line
        let sem = amb.p("sem_amostras.vcf");
        grava(&sem, &[cabecalho(&[]), vec![regs[0].fixas()]].concat(), "\n");
        let e = vcf_filter(&sem, &nao_criada, None, &[], Some(&nomes(&["S1"]))).unwrap_err();
        assert!(e.contains("no sample columns"), "{e}");
        let r = vcf_filter(&sem, &saida, None, &[], Some(&[])).unwrap();
        assert_eq!((r.amostras_antes, r.amostras_depois, r.mantidos), (0, 0, 1));
        let nu = amb.p("sem_chrom.vcf");
        grava(&nu, &[regs[0].linha(), regs[1].linha()], "\n");
        let e = vcf_filter(&nu, &nao_criada, None, &[], Some(&nomes(&["S1"]))).unwrap_err();
        assert!(e.contains("#CHROM") && e.contains("line 1"), "{e}");
        let r = vcf_filter(&nu, &saida, Some("POS>0"), &[], Some(&[])).unwrap();
        assert_eq!(r.mantidos, 2);
        iguais(&le(&saida), &[regs[0].fixas(), regs[1].fixas()], "no header at all, sites only");
        ensaios += 1;
    }
    {
        // empty file, header only, input == output, missing input
        std::fs::write(amb.p("vazio.vcf"), "").unwrap();
        let r = vcf_filter(&amb.p("vazio.vcf"), &saida, Some("QUAL>1"), &["chr1".to_string()], Some(&[])).unwrap();
        assert_eq!(r, FilterSummary::default());
        assert!(le(&saida).is_empty());
        grava(&amb.p("so_cab.vcf"), &cab_entrada, "\r\n");
        let r = vcf_filter(&amb.p("so_cab.vcf"), &saida, Some("QUAL>1"), &[], None).unwrap();
        assert_eq!((r.lidos, r.mantidos, r.amostras_antes), (0, 0, 3));
        assert_eq!(le(&saida).len(), cab_entrada.len() + 1);
        let e = vcf_filter(&entrada, &entrada, Some("QUAL>1"), &[], None).unwrap_err();
        assert!(e.contains("input file itself"), "{e}");
        assert_eq!(le(&entrada).len(), linhas.len(), "the input was damaged");
        assert!(vcf_filter(&amb.p("does_not_exist.vcf"), &nao_criada, None, &[], None).is_err());
        assert!(!std::path::Path::new(&nao_criada).exists());
        ensaios += 1;
    }
    ensaios
}

// ---------------------------------------------------------------- C) normalization

#[derive(Clone, Debug, PartialEq)]
struct NReg {
    chrom: String,
    pos: u64,
    id: String,
    rf: String,
    /// Empty: ALT is `.`.
    alts: Vec<String>,
    qual: String,
    filter: String,
    info: Vec<(String, Option<String>)>,
    /// Empty: the file has no FORMAT column.
    formato: Vec<String>,
    amostras: Vec<Vec<String>>,
}

impl NReg {
    fn linha(&self) -> String {
        let alt = if self.alts.is_empty() { ".".to_string() } else { self.alts.join(",") };
        let info = if self.info.is_empty() {
            ".".to_string()
        } else {
            let v: Vec<String> = self
                .info
                .iter()
                .map(|(k, v)| match v {
                    Some(v) => format!("{k}={v}"),
                    None => k.clone(),
                })
                .collect();
            v.join(";")
        };
        let mut l = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.chrom, self.pos, self.id, self.rf, alt, self.qual, self.filter, info
        );
        if !self.formato.is_empty() {
            l.push('\t');
            l.push_str(&self.formato.join(":"));
            for a in &self.amostras {
                l.push('\t');
                l.push_str(&a.join(":"));
            }
        }
        l
    }
}

const CROMS_N: [&str; 3] = ["chr1", "chr2", "chr10"];

fn cabecalho_norm(n_amostras: usize) -> Vec<String> {
    let mut c: Vec<String> = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=4000>",
        "##contig=<ID=chr2,length=4000>",
        "##contig=<ID=chr10,length=4000>",
        "##INFO=<ID=AI,Number=A,Type=Integer,Description=\"Index of the ALT allele\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency, one per ALT\">",
        "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Depths\">",
        "##INFO=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">",
        // a Description that tries to pass for other attributes
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth, not Number=A, ID=AF\">",
        "##INFO=<ID=XA,Number=.,Type=String,Description=\"Per ALT only in FORMAT\">",
        "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"In dbSNP\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depths\">",
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
        // the attributes in another order
        "##FORMAT=<Number=A,Type=Integer,ID=XA,Description=\"Per ALT, \\\"ID\\\" after Number\">",
    ]
    .iter()
    .map(|x| x.to_string())
    .collect();
    let mut col = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string();
    if n_amostras > 0 {
        col.push_str("\tFORMAT");
        for i in 0..n_amostras {
            col.push_str(&format!("\tS{}", i + 1));
        }
    }
    c.push(col);
    c
}

/// What the header above declares, written again by hand.
fn num_info(chave: &str) -> char {
    match chave {
        "AI" | "AF" => 'A',
        "AD" => 'R',
        "PL" => 'G',
        _ => '.',
    }
}

fn num_formato(chave: &str) -> char {
    match chave {
        "XA" => 'A',
        "AD" => 'R',
        "PL" => 'G',
        _ => '.',
    }
}

/// The values of ALT allele `k` (from 1) of `n`, by the rule of the
/// specification; unchanged when the count is not the one the Number implies.
fn parte(valor: &str, tipo: char, k: usize, n: usize) -> String {
    let v: Vec<&str> = valor.split(',').collect();
    // position of diploid genotype a/b (a <= b) in VCF order
    let g = |a: usize, b: usize| b * (b + 1) / 2 + a;
    match tipo {
        'A' if v.len() == n => v[k - 1].to_string(),
        'R' if v.len() == n + 1 => format!("{},{}", v[0], v[k]),
        'G' if v.len() == (n + 1) * (n + 2) / 2 => format!("{},{},{}", v[g(0, 0)], v[g(0, k)], v[g(k, k)]),
        _ => valor.to_string(),
    }
}

/// GT for the record of ALT allele `k`: k -> 1, 0 -> 0, another ALT -> `.`.
fn gt_esperado(gt: &str, k: usize) -> String {
    let mut o = String::new();
    let mut atual = String::new();
    for c in gt.chars().chain(std::iter::once('\n')) {
        if c == '/' || c == '|' || c == '\n' {
            o.push(if atual == "0" {
                '0'
            } else if atual == k.to_string() {
                '1'
            } else {
                '.'
            });
            atual.clear();
            if c != '\n' {
                o.push(c);
            }
        } else {
            atual.push(c);
        }
    }
    o
}

fn separa(r: &NReg, k: usize) -> NReg {
    let n = r.alts.len();
    let mut o = r.clone();
    o.alts = vec![r.alts[k - 1].clone()];
    for (chave, valor) in o.info.iter_mut() {
        if let Some(v) = valor {
            *v = parte(v, num_info(chave), k, n);
        }
    }
    for amostra in o.amostras.iter_mut() {
        for (i, campo) in amostra.iter_mut().enumerate() {
            let chave = &r.formato[i];
            *campo = if chave == "GT" { gt_esperado(campo, k) } else { parte(campo, num_formato(chave), k, n) };
        }
    }
    o
}

fn puro(x: &str) -> bool {
    !x.is_empty() && x.chars().all(|c| "ACGTNacgtn".contains(c))
}

/// Trims in place; returns (alleles trimmed, POS changed).
fn apara_esperado(r: &mut NReg) -> (bool, bool) {
    if r.alts.len() != 1 || !puro(&r.rf) || !puro(&r.alts[0]) {
        return (false, false);
    }
    let mut a: Vec<char> = r.rf.chars().collect();
    let mut b: Vec<char> = r.alts[0].chars().collect();
    let (tamanho, pos) = (a.len(), r.pos);
    while a.len() > 1 && b.len() > 1 && a[a.len() - 1].eq_ignore_ascii_case(&b[b.len() - 1]) {
        a.pop();
        b.pop();
    }
    while a.len() > 1 && b.len() > 1 && a[0].eq_ignore_ascii_case(&b[0]) {
        a.remove(0);
        b.remove(0);
        r.pos += 1;
    }
    let aparado = a.len() != tamanho;
    r.rf = a.into_iter().collect();
    r.alts[0] = b.into_iter().collect();
    (aparado, r.pos != pos)
}

struct Esperado {
    /// The normalized records in the order they are produced (not sorted).
    fluxo: Vec<NReg>,
    separados: usize,
    aparados: usize,
    pos_alteradas: usize,
}

fn normaliza_esperado(regs: &[NReg], separar: bool) -> Esperado {
    let mut e = Esperado { fluxo: Vec::new(), separados: 0, aparados: 0, pos_alteradas: 0 };
    for r in regs {
        let partes: Vec<NReg> = if separar && r.alts.len() > 1 {
            e.separados += 1;
            (1..=r.alts.len()).map(|k| separa(r, k)).collect()
        } else {
            vec![r.clone()]
        };
        for mut x in partes {
            let (aparado, moveu) = apara_esperado(&mut x);
            e.aparados += aparado as usize;
            e.pos_alteradas += moveu as usize;
            e.fluxo.push(x);
        }
    }
    e
}

fn posto(ordem: &OrdemCromossomos, c: &str) -> usize {
    let lexica = ["chr1", "chr10", "chr2"];
    match ordem {
        OrdemCromossomos::Lexica => lexica.iter().position(|x| *x == c).unwrap(),
        _ => CROMS_N.iter().position(|x| *x == c).unwrap(),
    }
}

/// An ALT allele derived from REF, so that shared bases at both ends are common.
fn gera_alt_norm(rf: &str, s: &mut u64) -> String {
    let mut a: Vec<u8> = rf.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let n = a.len();
    match lcg(s) % 20 {
        0 => return "<DEL>".to_string(),
        1 => return "*".to_string(),
        // one base changed, anywhere: a SNP with context around it
        2..=6 => {
            let i = (lcg(s) % n as u64) as usize;
            a[i] = outra(a[i], s);
        }
        // a stretch removed, leaving at least one base
        7..=9 if n >= 2 => {
            let ini = (lcg(s) % n as u64) as usize;
            let max = (n - 1).min(n - ini);
            let q = 1 + (lcg(s) % max as u64) as usize;
            a.drain(ini..ini + q);
        }
        // bases inserted, anywhere
        7..=13 => {
            let i = (lcg(s) % (n as u64 + 1)) as usize;
            for _ in 0..1 + lcg(s) % 3 {
                a.insert(i, base(s));
            }
        }
        // two bases changed
        14..=15 => {
            for _ in 0..2 {
                let i = (lcg(s) % n as u64) as usize;
                a[i] = outra(a[i], s);
            }
        }
        // the same as REF
        16 => {}
        _ => a = (0..1 + lcg(s) % 6).map(|_| base(s)).collect(),
    }
    muda_caixa(std::str::from_utf8(&a).unwrap(), s)
}

fn gera_gt(n: usize, s: &mut u64) -> String {
    let ploidia = match lcg(s) % 14 {
        0 => 1,
        1 => 3,
        _ => 2,
    };
    let mut g = String::new();
    for i in 0..ploidia {
        if i > 0 {
            g.push(if lcg(s) % 2 == 0 { '/' } else { '|' });
        }
        if lcg(s) % 10 == 0 {
            g.push('.');
        } else {
            g.push_str(&(lcg(s) % (n as u64 + 1)).to_string());
        }
    }
    g
}

/// `l_ref_max` = 1 gives records that no trimming can move.
fn gera_nreg(s: &mut u64, chrom: &str, genoma: &[u8], n_amostras: usize, l_ref_max: u64) -> NReg {
    let l_ref = 1 + (lcg(s) % l_ref_max) as usize;
    let pos = 1 + lcg(s) % (genoma.len() as u64 - 10);
    let i = pos as usize - 1;
    let rf = muda_caixa(std::str::from_utf8(&genoma[i..i + l_ref]).unwrap(), s);
    let n = match lcg(s) % 12 {
        0 => 0,
        1..=4 => 2 + (lcg(s) % 3) as usize,
        _ => 1,
    };
    let alts: Vec<String> = (0..n).map(|_| gera_alt_norm(&rf, s)).collect();
    let g = (n + 1) * (n + 2) / 2;

    let mut info: Vec<(String, Option<String>)> = Vec::new();
    let mut poe = |k: &str, v: Option<String>| info.push((k.to_string(), v));
    if n > 0 {
        poe("AI", Some((1..=n).map(|k| k.to_string()).collect::<Vec<_>>().join(",")));
        match lcg(s) % 10 {
            0 => poe("AF", Some(".".to_string())),
            1 => poe("AF", Some(lista(n + 1, s))), // one value too many: copied as it is
            2 => {}
            _ => poe("AF", Some(lista(n, s))),
        }
        if lcg(s) % 2 == 0 {
            poe("AD", Some(lista(n + 1, s)));
        }
        if lcg(s) % 3 == 0 {
            poe("PL", Some(lista(g, s)));
        }
        if lcg(s) % 4 == 0 {
            poe("XA", Some(lista(n, s)));
        }
    }
    if lcg(s) % 2 == 0 {
        poe("DP", Some((lcg(s) % 90).to_string()));
    }
    if lcg(s) % 4 == 0 {
        poe("DB", None);
    }
    if lcg(s) % 5 == 0 {
        poe("UNK", Some("1,2,3".to_string())); // not declared in the header
    }

    let formato: Vec<String> = if n_amostras == 0 {
        Vec::new()
    } else {
        escolhe(s, &["GT:AD:PL:DP:XA", "GT:AD:PL:DP:XA", "GT:DP", "DP:GT:XA", "AD:PL", "GT"])
            .split(':')
            .map(str::to_string)
            .collect()
    };
    let amostras: Vec<Vec<String>> = (0..n_amostras)
        .map(|_| {
            if lcg(s) % 15 == 0 {
                return vec![".".to_string()];
            }
            let mut campos: Vec<String> = formato
                .iter()
                .map(|k| match k.as_str() {
                    "GT" => gera_gt(n, s),
                    "AD" if lcg(s) % 10 == 0 => ".".to_string(),
                    "AD" => lista(n + 1, s),
                    "PL" => match lcg(s) % 12 {
                        0 => ".".to_string(),
                        1 => lista(n + 1, s), // one per allele, as in a haploid call: copied as it is
                        _ => lista(g, s),
                    },
                    "XA" if n == 0 => ".".to_string(),
                    "XA" => lista(n, s),
                    _ => (lcg(s) % 99).to_string(),
                })
                .collect();
            if lcg(s) % 8 == 0 {
                // trailing fields may be dropped
                campos.truncate(1 + (lcg(s) % campos.len() as u64) as usize);
            }
            campos
        })
        .collect();

    NReg {
        chrom: chrom.to_string(),
        pos,
        id: String::new(),
        rf,
        alts,
        qual: (lcg(s) % 90).to_string(),
        filter: escolhe(s, &["PASS", "PASS", "q10", "."]).to_string(),
        info,
        formato,
        amostras,
    }
}

/// The sequence that results from applying the variant to the reference.
fn aplica(genoma: &[u8], pos: u64, rf: &str, alt: &str, contexto: &str) -> Vec<u8> {
    let i = pos as usize - 1;
    assert!(
        i + rf.len() <= genoma.len() && genoma[i..i + rf.len()].eq_ignore_ascii_case(rf.as_bytes()),
        "{contexto}: REF `{rf}` at POS {pos} is not what the reference has there"
    );
    let mut v = genoma[..i].to_vec();
    v.extend(alt.bytes().map(|b| b.to_ascii_uppercase()));
    v.extend_from_slice(&genoma[i + rf.len()..]);
    v
}

struct Dados {
    genomas: Vec<Vec<u8>>,
    regs: Vec<NReg>,
    cab: Vec<String>,
}

/// A sorted VCF over random reference sequences.
fn gera_dados(s: &mut u64, n: usize, n_amostras: usize, l_ref_max: u64) -> Dados {
    let genomas: Vec<Vec<u8>> = CROMS_N.iter().map(|_| (0..4000).map(|_| base(s)).collect()).collect();
    let mut regs: Vec<NReg> = Vec::new();
    for (c, chrom) in CROMS_N.iter().enumerate() {
        let mut deste: Vec<NReg> = (0..n / 3).map(|_| gera_nreg(s, chrom, &genomas[c], n_amostras, l_ref_max)).collect();
        deste.sort_by_key(|r| r.pos);
        regs.extend(deste);
    }
    for (i, r) in regs.iter_mut().enumerate() {
        r.id = format!("r{i}");
    }
    Dados { genomas, regs, cab: cabecalho_norm(n_amostras) }
}

/// Runs vcf_normalize and checks the result in every way this file knows.
#[allow(clippy::too_many_arguments)]
fn confere_normaliza(
    amb: &Amb,
    d: &Dados,
    entrada: &str,
    saida: &str,
    separar: bool,
    mem: usize,
    ordem: OrdemCromossomos,
    contexto: &str,
) -> NormSummary {
    let r = vcf_normalize(entrada, saida, separar, None, &amb.op(mem, ordem.clone())).unwrap_or_else(|e| panic!("{contexto}: {e}"));
    assert!(amb.sem_temporarios(), "{contexto}: a temporary file was left behind");
    let got = le(saida);

    // the header: the input header plus one line before #CHROM
    let mut esp_cab = d.cab.clone();
    esp_cab.insert(esp_cab.len() - 1, format!("##MMGEN_normalize=<split_multiallelic={separar}>"));
    iguais(&cab(&got), &esp_cab, &format!("{contexto}: header"));
    let got = dados(&got);

    // the records: exactly the independent normalization, sorted (stably) only if a POS changed
    let e = normaliza_esperado(&d.regs, separar);
    let fluxo: Vec<String> = e.fluxo.iter().map(NReg::linha).collect();
    let mut ordenado = e.fluxo.clone();
    if e.pos_alteradas > 0 {
        ordenado.sort_by_key(|x| (posto(&ordem, &x.chrom), x.pos));
    }
    iguais(&got, &ordenado.iter().map(NReg::linha).collect::<Vec<_>>(), &format!("{contexto}: records"));
    assert_eq!(
        (r.lidos, r.gravados, r.separados, r.aparados, r.pos_alteradas, r.reordenado, r.ignoradas),
        (d.regs.len(), fluxo.len(), e.separados, e.aparados, e.pos_alteradas, e.pos_alteradas > 0, 0),
        "{contexto}: summary"
    );

    // (v) sorted, and the same multiset of records as the unsorted stream
    if r.reordenado {
        let chave = |l: &str| {
            let mut c = l.split('\t');
            let crom = posto(&ordem, c.next().unwrap());
            (crom, c.next().unwrap().parse::<u64>().unwrap())
        };
        for w in got.windows(2) {
            assert!(chave(&w[0]) <= chave(&w[1]), "{contexto}: out of order:\n  {}\n  {}", w[0], w[1]);
        }
        let (mut a, mut b) = (got.clone(), fluxo.clone());
        a.sort();
        b.sort();
        assert!(a == b, "{contexto}: the sorted output is not the same multiset of records");
    }

    // properties, read from the output itself
    let mut conferidos = 0usize;
    for l in &got {
        let c: Vec<&str> = l.split('\t').collect();
        let (pos, id, rf, alt, info) = (c[1].parse::<u64>().unwrap(), c[2], c[3], c[4], c[7]);
        let orig = &d.regs[id[1..].parse::<usize>().unwrap()];
        // (i) one ALT per record after splitting
        if separar {
            assert!(!alt.contains(','), "{contexto}: still multiallelic: {l}");
        }
        if alt.contains(',') {
            assert_eq!(l, &orig.linha(), "{contexto}: a multiallelic record must be untouched without splitting");
            continue;
        }
        if !puro(rf) || !puro(alt) {
            continue;
        }
        // (ii) minimal alleles
        if rf.len() > 1 && alt.len() > 1 {
            let (r, a) = (rf.as_bytes(), alt.as_bytes());
            assert!(!r[r.len() - 1].eq_ignore_ascii_case(&a[a.len() - 1]), "{contexto}: shared last base: {l}");
            assert!(!r[0].eq_ignore_ascii_case(&a[0]), "{contexto}: shared first base: {l}");
        }
        // (iii) the variant still means the same thing
        let k = if orig.alts.len() == 1 {
            1
        } else {
            let ai = info.split(';').find_map(|x| x.strip_prefix("AI=")).unwrap_or_else(|| panic!("no AI in {l}"));
            ai.parse::<usize>().unwrap_or_else(|_| panic!("{contexto}: AI was not split in {l}"))
        };
        let genoma = &d.genomas[CROMS_N.iter().position(|x| *x == c[0]).unwrap()];
        let antes = aplica(genoma, orig.pos, &orig.rf, &orig.alts[k - 1], contexto);
        let depois = aplica(genoma, pos, rf, alt, contexto);
        assert!(antes == depois, "{contexto}: the variant changed meaning:\n  before: {}\n  after:  {l}", orig.linha());
        conferidos += 1;
    }
    assert!(d.regs.is_empty() || conferidos > got.len() / 4, "{contexto}: too few records checked against the reference");

    let j = json(&r.json());
    assert_eq!(
        j.chaves(),
        ["read", "written", "split", "trimmed", "position_changed", "left_aligned", "ref_mismatch", "not_in_reference", "resorted", "skipped", "runs"]
    );
    assert_eq!(
        (j.campo("read").num(), j.campo("written").num(), j.campo("split").num(), j.campo("trimmed").num()),
        (r.lidos as f64, r.gravados as f64, r.separados as f64, r.aparados as f64)
    );
    assert_eq!(
        (j.campo("position_changed").num(), j.campo("resorted"), j.campo("runs").num()),
        (r.pos_alteradas as f64, &J::Bool(r.reordenado), r.lotes as f64)
    );

    // normalizing what is already normalized changes nothing
    let de_novo = amb.p("again.vcf.lz4");
    let r2 = vcf_normalize(saida, &de_novo, separar, None, &amb.op(mem, ordem)).unwrap();
    assert_eq!(
        (r2.lidos, r2.gravados, r2.separados, r2.aparados, r2.pos_alteradas, r2.reordenado),
        (got.len(), got.len(), 0, 0, 0, false),
        "{contexto}: second pass"
    );
    iguais(&dados(&le(&de_novo)), &got, &format!("{contexto}: second pass"));
    r
}

/// A real file, if one is given on the command line: the three tools run on
/// it and are checked against each other.
fn arquivo_real(amb: &Amb, arq: &str) {
    let relogio = std::time::Instant::now();
    let st = vcf_stats(arq).unwrap();
    print!("\n[{:.1} s] {arq}\n{}", relogio.elapsed().as_secs_f64(), st.resumo());

    // the filter and the statistics count the same sites
    let saida = amb.p("real_filtered.vcf.lz4");
    for (expr, esperado) in [
        ("TYPE~\"snp\"", st.snps),
        ("TYPE==\"indel\" && N_ALT==1", usize::MAX),
        ("N_ALT>1", st.multiallelic),
        ("FILTER==\"PASS\"", st.pass),
        ("QUAL==\".\"", st.qual_missing),
        ("ALT==\".\"", st.no_alt),
    ] {
        let relogio = std::time::Instant::now();
        let f = vcf_filter(arq, &saida, Some(expr), &[], None).unwrap();
        println!("[{:.1} s] filter {expr}: kept {} of {}", relogio.elapsed().as_secs_f64(), f.mantidos, f.lidos);
        assert_eq!(f.lidos, st.records);
        if esperado != usize::MAX {
            assert_eq!(f.mantidos, esperado, "`{expr}` and vcf_stats disagree");
        }
    }

    // after splitting there is no multiallelic site left, and a second pass finds nothing to do
    let normal = amb.p("real_normalized.vcf.lz4");
    let relogio = std::time::Instant::now();
    let n = vcf_normalize(arq, &normal, true, None, &amb.op(256 << 20, OrdemCromossomos::Cabecalho)).unwrap();
    println!("[{:.1} s] normalize: {}", relogio.elapsed().as_secs_f64(), n.json());
    let depois = vcf_stats(&normal).unwrap();
    print!("{}", depois.resumo());
    assert_eq!((depois.records, depois.multiallelic, n.lidos), (n.gravados, 0, st.records));
    // trimming can only reveal SNPs (`AC -> AT` becomes `C -> T`), never lose one
    assert!(depois.transitions >= st.transitions && depois.transversions >= st.transversions);
    let relogio = std::time::Instant::now();
    let n2 = vcf_normalize(&normal, &amb.p("real_again.vcf.lz4"), true, None, &amb.op(256 << 20, OrdemCromossomos::Cabecalho)).unwrap();
    println!("[{:.1} s] normalize again: {}", relogio.elapsed().as_secs_f64(), n2.json());
    assert_eq!((n2.lidos, n2.separados, n2.aparados, n2.reordenado), (n.gravados, 0, 0, false));
    assert!(amb.sem_temporarios());
}

fn testa_normaliza(amb: &Amb, s: &mut u64) -> usize {
    let mut ensaios = 0;
    let natural = || OrdemCromossomos::Natural;

    // ---- a file written by hand, with the expected output written by hand
    {
        let cabecalho = [
            "##fileformat=VCFv4.2",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency, per ALT\">",
            "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Depths\">",
            "##INFO=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depths\">",
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
        ];
        let am = "GT:AD:PL\t0/1:5,4:10,0,20\t1|1:0,9:90,9,0\t./.:.:.";
        let entrada = [
            format!("chr1\t100\ta\tACGT\tACT\t9\tPASS\tDP=3\t{am}"),
            format!("chr1\t200\tb\tGCACA\tGCA\t9\tPASS\tDP=3\t{am}"),
            "chr1\t300\tc\tACG\tATG,A,ACGT\t50\tq10\tAF=0.1,0.2,0.3;AD=10,5,3,1;PL=0,1,2,3,4,5,6,7,8,9;DP=7;DB\tGT:AD:PL\t1/2:10,5,3,1:0,1,2,3,4,5,6,7,8,9\t2|3:1,2,3,4:9,8,7,6,5,4,3,2,1,0\t0/3:7,0,0,6:.".to_string(),
            "chr1\t301\td\tT\tG\t9\tPASS\t.\tGT:AD:PL\t0/0:9,0:0,9,90\t0/1:5,4:10,0,20\t1/1:0,9:90,9,0".to_string(),
            "chr1\t400\te\ta\t<DEL>,aT\t9\tPASS\tAF=0.5,0.25\tGT\t1/2\t0|2\t2".to_string(),
            "chr1\t500\tf\tacg\taTg\t9\tPASS\t.\tGT\t0/1\t0/1\t0/1".to_string(),
            "chr1\t600\tg\tA\t.\t9\tPASS\t.\tGT\t0/0\t0/0\t0/0".to_string(),
        ];
        let mut linhas: Vec<String> = cabecalho.iter().map(|x| x.to_string()).collect();
        linhas.extend(entrada.iter().cloned());
        linhas.push("chr1\tnot_a_number\th\tA\tT\t9\tPASS\t.".to_string());
        grava(&amb.p("mao.vcf"), &linhas, "\r\n");

        // splitting: 1/2 becomes 1/. and ./1; Number=A, R and G keep the values of the allele
        let esperado = [
            format!("chr1\t101\ta\tCG\tC\t9\tPASS\tDP=3\t{am}"),
            format!("chr1\t200\tb\tGCA\tG\t9\tPASS\tDP=3\t{am}"),
            "chr1\t300\tc\tACG\tA\t50\tq10\tAF=0.2;AD=10,3;PL=0,3,5;DP=7;DB\tGT:AD:PL\t./1:10,3:0,3,5\t1|.:1,3:9,6,4\t0/.:7,0:.".to_string(),
            "chr1\t301\tc\tC\tT\t50\tq10\tAF=0.1;AD=10,5;PL=0,1,2;DP=7;DB\tGT:AD:PL\t1/.:10,5:0,1,2\t.|.:1,2:9,8,7\t0/.:7,0:.".to_string(),
            "chr1\t301\td\tT\tG\t9\tPASS\t.\tGT:AD:PL\t0/0:9,0:0,9,90\t0/1:5,4:10,0,20\t1/1:0,9:90,9,0".to_string(),
            "chr1\t302\tc\tG\tGT\t50\tq10\tAF=0.3;AD=10,1;PL=0,6,9;DP=7;DB\tGT:AD:PL\t./.:10,1:0,6,9\t.|1:1,4:9,3,0\t0/1:7,6:.".to_string(),
            "chr1\t400\te\ta\t<DEL>\t9\tPASS\tAF=0.5\tGT\t1/.\t0|.\t.".to_string(),
            "chr1\t400\te\ta\taT\t9\tPASS\tAF=0.25\tGT\t./1\t0|1\t1".to_string(),
            "chr1\t501\tf\tc\tT\t9\tPASS\t.\tGT\t0/1\t0/1\t0/1".to_string(),
            "chr1\t600\tg\tA\t.\t9\tPASS\t.\tGT\t0/0\t0/0\t0/0".to_string(),
        ];
        let r = vcf_normalize(&amb.p("mao.vcf"), &amb.p("mao_out.vcf"), true, None, &amb.op(usize::MAX, natural())).unwrap();
        let got = le(&amb.p("mao_out.vcf"));
        iguais(&dados(&got), &esperado, "hand-written file, splitting");
        assert_eq!(got[cabecalho.len() - 1], "##MMGEN_normalize=<split_multiallelic=true>");
        assert_eq!(got[cabecalho.len()], cabecalho[cabecalho.len() - 1]);
        assert_eq!(
            r,
            NormSummary {
                lidos: 7,
                gravados: 10,
                separados: 2,
                aparados: 5,
                pos_alteradas: 4,
                reordenado: true,
                ignoradas: 1,
                lotes: 0,
                ..NormSummary::default()
            }
        );
        ensaios += 1;

        // without splitting: the multiallelic records stay as they are
        let esperado = [
            format!("chr1\t101\ta\tCG\tC\t9\tPASS\tDP=3\t{am}"),
            format!("chr1\t200\tb\tGCA\tG\t9\tPASS\tDP=3\t{am}"),
            entrada[2].clone(),
            entrada[3].clone(),
            entrada[4].clone(),
            "chr1\t501\tf\tc\tT\t9\tPASS\t.\tGT\t0/1\t0/1\t0/1".to_string(),
            entrada[6].clone(),
        ];
        let r = vcf_normalize(&amb.p("mao.vcf"), &amb.p("mao_out.vcf.gz"), false, None, &amb.op(usize::MAX, natural())).unwrap();
        let got = le(&amb.p("mao_out.vcf.gz"));
        iguais(&dados(&got), &esperado, "hand-written file, no splitting");
        assert_eq!(got[cabecalho.len() - 1], "##MMGEN_normalize=<split_multiallelic=false>");
        assert_eq!(
            r,
            NormSummary {
                lidos: 7,
                gravados: 7,
                separados: 0,
                aparados: 3,
                pos_alteradas: 2,
                reordenado: true,
                ignoradas: 1,
                lotes: 0,
                ..NormSummary::default()
            }
        );
        ensaios += 1;

        // trimming that moves no POS: no sorting, and the input order is kept even if it was not sorted
        let fora_de_ordem = [
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
            "chr2\t50\tx\tT\tC\t1\tPASS\t.".to_string(),
            "chr1\t200\tb\tGCACA\tGCA\t9\tPASS\tDP=3".to_string(),
            "chr1\t10\ty\tAC\tA,ACC\t9\tPASS\t.".to_string(),
        ];
        grava(&amb.p("fora.vcf"), &fora_de_ordem, "\n");
        for saida in ["fora_out.vcf", "fora_out.vcf.gz", "fora_out.vcf.lz4"] {
            let r = vcf_normalize(&amb.p("fora.vcf"), &amb.p(saida), true, None, &amb.op(usize::MAX, natural())).unwrap();
            let esp = [
                "##MMGEN_normalize=<split_multiallelic=true>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr2\t50\tx\tT\tC\t1\tPASS\t.",
                "chr1\t200\tb\tGCA\tG\t9\tPASS\tDP=3",
                "chr1\t10\ty\tAC\tA\t9\tPASS\t.",
                "chr1\t10\ty\tA\tAC\t9\tPASS\t.", // AC -> ACC loses the shared last base; POS stays
            ];
            assert_eq!(le(&amb.p(saida)), esp, "{saida}");
            assert_eq!(
                r,
                NormSummary {
                    lidos: 3,
                    gravados: 4,
                    separados: 1,
                    aparados: 2,
                    pos_alteradas: 0,
                    reordenado: false,
                    ignoradas: 0,
                    lotes: 0,
                    ..NormSummary::default()
                }
            );
            assert!(amb.sem_temporarios());
        }
        ensaios += 1;
    }

    // ---- random records over a random reference: every property, several configurations
    let d = gera_dados(s, 3000, 3, 6);
    let mut linhas = d.cab.clone();
    linhas.extend(d.regs.iter().map(NReg::linha));
    grava(&amb.p("n.vcf"), &linhas, "\n");
    grava(&amb.p("n_crlf.vcf.gz"), &linhas, "\r\n");
    for (entrada, saida, separar, ordem) in [
        ("n.vcf", "n_out.vcf", true, OrdemCromossomos::Natural),
        ("n.vcf", "n_out.vcf", false, OrdemCromossomos::Natural),
        ("n.vcf", "n_out.vcf.gz", true, OrdemCromossomos::Cabecalho), // (vi) a .gz output read back
        ("n_crlf.vcf.gz", "n_out.vcf.lz4", true, OrdemCromossomos::Lexica),
        ("n_crlf.vcf.gz", "n_out.vcf.gz", false, OrdemCromossomos::Lexica),
    ] {
        let contexto = format!("{entrada} -> {saida}, split={separar}, {ordem:?}");
        let r = confere_normaliza(amb, &d, &amb.p(entrada), &amb.p(saida), separar, usize::MAX, ordem, &contexto);
        assert!(r.reordenado && r.lotes == 0 && r.aparados > r.pos_alteradas && r.pos_alteradas > 100, "{contexto}: {r:?}");
        assert!(separar == (r.separados > 0) && separar == (r.gravados > r.lidos), "{contexto}: {r:?}");
        ensaios += 1;
    }

    // ---- no samples at all
    let d = gera_dados(s, 1500, 0, 6);
    let mut linhas = d.cab.clone();
    linhas.extend(d.regs.iter().map(NReg::linha));
    grava(&amb.p("n0.vcf.lz4"), &linhas, "\n");
    for separar in [true, false] {
        confere_normaliza(amb, &d, &amb.p("n0.vcf.lz4"), &amb.p("n0_out.vcf"), separar, usize::MAX, natural(), "sites only");
        ensaios += 1;
    }

    // ---- nothing to trim: no POS changes, so no sorting, whatever the output compression
    let d = gera_dados(s, 3000, 2, 1);
    let mut linhas = d.cab.clone();
    linhas.extend(d.regs.iter().map(NReg::linha));
    grava(&amb.p("n1.vcf"), &linhas, "\n");
    for (saida, separar) in [("n1_out.vcf", true), ("n1_out.vcf.gz", true), ("n1_out.vcf.lz4", false)] {
        let r = confere_normaliza(amb, &d, &amb.p("n1.vcf"), &amb.p(saida), separar, 1, natural(), saida);
        assert!(!r.reordenado && r.pos_alteradas == 0 && r.aparados == 0 && r.lotes == 0, "{saida}: {r:?}");
        ensaios += 1;
    }

    // ---- (v) a file larger than the memory budget: the sort goes through runs on disk
    let d = gera_dados(s, 30_000, 3, 6);
    let mut linhas = d.cab.clone();
    linhas.extend(d.regs.iter().map(NReg::linha));
    grava(&amb.p("ng.vcf.gz"), &linhas, "\n");
    for (saida, separar, ordem) in [
        ("ng_out.vcf", true, OrdemCromossomos::Natural),
        ("ng_out.vcf.gz", true, OrdemCromossomos::Lexica),
        ("ng_out.vcf.lz4", false, OrdemCromossomos::Cabecalho),
    ] {
        let contexto = format!("tiny memory budget -> {saida}, split={separar}, {ordem:?}");
        let r = confere_normaliza(amb, &d, &amb.p("ng.vcf.gz"), &amb.p(saida), separar, 1, ordem, &contexto);
        assert!(r.reordenado && r.lotes >= 2, "{contexto}: the sort should have used runs on disk: {r:?}");
        ensaios += 1;
    }

    // ---- empty file, header only, input == output, missing input
    {
        std::fs::write(amb.p("vazio.vcf"), "").unwrap();
        let r = vcf_normalize(&amb.p("vazio.vcf"), &amb.p("vazio_out.vcf"), true, None, &amb.op(1, natural())).unwrap();
        assert_eq!(r, NormSummary::default());
        assert!(le(&amb.p("vazio_out.vcf")).is_empty());
        let cab_n = cabecalho_norm(2);
        grava(&amb.p("so_cab.vcf"), &cab_n, "\r\n");
        let r = vcf_normalize(&amb.p("so_cab.vcf"), &amb.p("so_cab_out.vcf.gz"), true, None, &amb.op(1, natural())).unwrap();
        assert_eq!(r, NormSummary::default());
        let mut esp = cab_n.clone();
        esp.insert(esp.len() - 1, "##MMGEN_normalize=<split_multiallelic=true>".to_string());
        iguais(&le(&amb.p("so_cab_out.vcf.gz")), &esp, "header only");
        let e = vcf_normalize(&amb.p("n.vcf"), &amb.p("n.vcf"), true, None, &amb.op(1, natural())).unwrap_err();
        assert!(e.contains("input file itself"), "{e}");
        let e = vcf_normalize(&amb.p("does_not_exist.vcf"), &amb.p("x.vcf"), true, None, &amb.op(1, natural())).unwrap_err();
        assert!(e.contains("does_not_exist.vcf") && !std::path::Path::new(&amb.p("x.vcf")).exists(), "{e}");
        assert!(amb.sem_temporarios());
        ensaios += 1;
    }
    ensaios
}

// ---------------------------------------------------------------- D) left-alignment against a reference

/// A sequence full of repeats (runs of one base, of two, of three), where an
/// indel has several equivalent positions.
fn genoma_com_repeticoes(n: usize, s: &mut u64) -> Vec<u8> {
    let mut g = Vec::with_capacity(n + 64);
    while g.len() < n {
        match lcg(s) % 4 {
            0 => {
                let b = base(s);
                g.extend(std::iter::repeat(b).take(2 + (lcg(s) % 10) as usize));
            }
            1 => {
                let u = [base(s), base(s)];
                for _ in 0..2 + lcg(s) % 6 {
                    g.extend_from_slice(&u);
                }
            }
            2 => {
                let u = [base(s), base(s), base(s)];
                for _ in 0..2 + lcg(s) % 4 {
                    g.extend_from_slice(&u);
                }
            }
            _ => {
                for _ in 0..1 + lcg(s) % 8 {
                    g.push(base(s));
                }
            }
        }
    }
    g.truncate(n);
    g
}

const JANELA_BRUTA: usize = 400;

fn txt(b: &[u8]) -> String {
    String::from_utf8(b.to_vec()).unwrap()
}

/// The leftmost form of the deletion of `g[p..p+l]`, by comparing haplotypes
/// (no shifting rule): the first position, from the left, whose deletion
/// gives the same sequence. (POS, REF, ALT), with the anchor base on the left.
fn delecao_esperada(g: &[u8], p: usize, l: usize) -> (u64, String, String) {
    let (w0, w1) = (p - JANELA_BRUTA, (p + l + JANELA_BRUTA).min(g.len()));
    let hap = |q: usize| [&g[w0..q], &g[q + l..w1]].concat();
    let alvo = hap(p);
    let q = (w0 + 1..=p).find(|&q| hap(q) == alvo).unwrap();
    assert!(q > w0 + 1, "a repeat longer than the brute-force window");
    (q as u64, txt(&g[q - 1..q + l]), txt(&g[q - 1..q]))
}

/// The leftmost form of the insertion of `ins` before `g[gap]`, the same way.
fn insercao_esperada(g: &[u8], gap: usize, ins: &[u8]) -> (u64, String, String) {
    let (w0, w1, l) = (gap - JANELA_BRUTA, (gap + JANELA_BRUTA).min(g.len()), ins.len());
    let alvo = [&g[w0..gap], ins, &g[gap..w1]].concat();
    let q = (w0 + 1..=gap)
        .find(|&q| [&g[w0..q], &alvo[q - w0..q - w0 + l], &g[q..w1]].concat() == alvo)
        .unwrap();
    assert!(q > w0 + 1, "a repeat longer than the brute-force window");
    (q as u64, txt(&g[q - 1..q]), txt(&[&g[q - 1..q], &alvo[q - w0..q - w0 + l]].concat()))
}

/// What a record must become.
#[derive(Clone)]
enum Espera {
    /// Exactly this POS, REF and ALT (upper case).
    Exato(u64, String, String),
    /// The input line, untouched.
    Igual,
    /// Something with the same meaning against the reference.
    MesmoHaplotipo,
}

struct AReg {
    chrom: String,
    pos: u64,
    id: String,
    rf: String,
    alt: String,
    espera: Espera,
    /// 'a' aligned, 'r' REF differs from the reference, 'f' sequence not in it.
    classe: char,
}

impl AReg {
    fn linha(&self) -> String {
        format!("{}\t{}\t{}\t{}\t{}\t30\tPASS\tDP=7", self.chrom, self.pos, self.id, self.rf, self.alt)
    }
}

/// Writes a FASTA with `largura` bases per line (random lengths when 0) and
/// returns its `.fai` (meaningful only for fixed-width lines).
fn grava_fasta(p: &str, seqs: &[(&str, &[u8])], largura: usize, fim: &str, s: &mut u64) -> String {
    let mut texto = String::new();
    let mut fai = String::new();
    for (nome, g) in seqs {
        texto.push_str(&format!(">{nome} a description\n").replace('\n', fim));
        fai.push_str(&format!("{nome}\t{}\t{}\t{largura}\t{}\n", g.len(), texto.len(), largura + fim.len()));
        let mut i = 0;
        while i < g.len() {
            let n = if largura == 0 { 10 + (lcg(s) % 70) as usize } else { largura };
            let pedaco = &g[i..(i + n).min(g.len())];
            // soft-masked stretches, as in the real genomes
            if lcg(s) % 5 == 0 {
                texto.push_str(&txt(pedaco).to_ascii_lowercase());
            } else {
                texto.push_str(&txt(pedaco));
            }
            texto.push_str(fim);
            i += n;
        }
    }
    let mut w = Escritor::cria(p).unwrap();
    w.write_all(texto.as_bytes()).unwrap();
    w.termina().unwrap();
    fai
}

fn testa_alinhamento(amb: &Amb, s: &mut u64) -> usize {
    let mut ensaios = 0;
    // The FASTA names them chr1 and 2; the VCF calls them 1 and chr2.
    let g1 = genoma_com_repeticoes(60_000, s);
    let g2 = genoma_com_repeticoes(20_000, s);
    let mut regs: Vec<AReg> = Vec::new();
    for (chrom, g) in [("1", &g1), ("chr2", &g2)] {
        for _ in 0..1500 {
            let p = 1000 + (lcg(s) % (g.len() as u64 - 2000)) as usize;
            let (a, b) = ((lcg(s) % 3) as usize, (lcg(s) % 3) as usize);
            let reg = |pos: u64, rf: String, alt: String, espera: Espera| AReg {
                chrom: chrom.to_string(),
                pos,
                id: String::new(),
                rf,
                alt,
                espera,
                classe: 'a',
            };
            regs.push(match lcg(s) % 10 {
                // deletion, with a bases of context on the left and b on the right
                0..=3 => {
                    let l = 1 + (lcg(s) % 6) as usize;
                    let (pos, rf, alt) = delecao_esperada(g, p, l);
                    reg(
                        (p - a) as u64,
                        txt(&g[p - 1 - a..p + l + b]),
                        txt(&[&g[p - 1 - a..p], &g[p + l..p + l + b]].concat()),
                        Espera::Exato(pos, rf, alt),
                    )
                }
                // insertion: mostly a copy of what follows or precedes it, which can shift
                4..=6 => {
                    let l = 1 + (lcg(s) % 6) as usize;
                    let ins: Vec<u8> = match lcg(s) % 4 {
                        0 | 1 => g[p..p + l].to_vec(),
                        2 => g[p - l..p].to_vec(),
                        _ => (0..l).map(|_| base(s)).collect(),
                    };
                    let (pos, rf, alt) = insercao_esperada(g, p, &ins);
                    reg(
                        (p - a) as u64,
                        txt(&g[p - 1 - a..p + b]),
                        txt(&[&g[p - 1 - a..p], &ins[..], &g[p..p + b]].concat()),
                        Espera::Exato(pos, rf, alt),
                    )
                }
                // SNP and MNP: nothing moves
                7 => reg(p as u64 + 1, txt(&g[p..p + 1]), (outra(g[p], s) as char).to_string(), Espera::Igual),
                8 => {
                    let alt: String = [outra(g[p], s) as char, outra(g[p + 1], s) as char].iter().collect();
                    reg(p as u64 + 1, txt(&g[p..p + 2]), alt, Espera::Igual)
                }
                // complex: length and bases change; only the meaning is checked
                _ => {
                    let alt: String = (0..1 + lcg(s) % 2).map(|_| base(s) as char).collect();
                    reg(p as u64 + 1, txt(&g[p..p + 3]), alt, Espera::MesmoHaplotipo)
                }
            });
        }
    }
    // REF that is not the reference (a SNP, and an indel that is only trimmed), and a sequence it lacks
    let p = 5000usize;
    let errado = (outra(g1[p], s) as char).to_string();
    for (chrom, pos, rf, alt, espera, classe) in [
        ("1", p as u64 + 1, errado.clone(), "N".to_string(), Espera::Igual, 'r'),
        ("1", p as u64 + 1, format!("{errado}TTTT"), format!("{errado}T"), Espera::Exato(p as u64 + 1, format!("{errado}TTT"), errado.clone()), 'r'),
        ("chrZ", 100, "CAAA".to_string(), "CA".to_string(), Espera::Exato(100, "CAA".to_string(), "C".to_string()), 'f'),
        ("chrZ", 200, "A".to_string(), "G".to_string(), Espera::Igual, 'f'),
    ] {
        regs.push(AReg { chrom: chrom.to_string(), pos, id: String::new(), rf, alt, espera, classe });
    }
    for (i, r) in regs.iter_mut().enumerate() {
        r.id = format!("v{i}");
    }
    let mut linhas: Vec<String> = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=1,length=60000>",
        "##contig=<ID=chr2,length=20000>",
        "##contig=<ID=chrZ,length=1000>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    .iter()
    .map(|x| x.to_string())
    .collect();
    linhas.extend(regs.iter().map(AReg::linha));
    let entrada = amb.p("al.vcf");
    grava(&entrada, &linhas, "\n");

    // what the counters must say
    let mut esp = NormSummary { lidos: regs.len(), gravados: regs.len(), reordenado: true, ..NormSummary::default() };
    for r in &regs {
        match r.classe {
            'r' => esp.ref_divergente += 1,
            'f' => esp.fora_da_referencia += 1,
            _ => {}
        }
        if let Espera::Exato(pos, rf, alt) = &r.espera {
            esp.aparados += usize::from(rf.len() + alt.len() < r.rf.len() + r.alt.len());
            esp.alinhados += usize::from(*pos < r.pos);
            esp.pos_alteradas += usize::from(*pos != r.pos);
        }
    }

    let confere = |saida: &str, r: &NormSummary, contexto: &str| {
        let got = le(saida);
        assert!(
            got.iter().any(|l| l.starts_with("##MMGEN_normalize=<split_multiallelic=true,reference=")),
            "{contexto}: the header does not record the reference"
        );
        let got = dados(&got);
        assert_eq!(got.len(), regs.len(), "{contexto}");
        let posto = |c: &str| ["1", "chr2", "chrZ"].iter().position(|x| *x == c).unwrap();
        let mut anterior = (0usize, 0u64);
        let mut mesmo_hap = 0;
        // the complex variants are checked by meaning; their counts come from the output
        let mut esp = esp.clone();
        for l in &got {
            let c: Vec<&str> = l.split('\t').collect();
            let pos: u64 = c[1].parse().unwrap();
            let chave = (posto(c[0]), pos);
            assert!(chave >= anterior, "{contexto}: out of order at {l}");
            anterior = chave;
            let orig = &regs[c[2][1..].parse::<usize>().unwrap()];
            match &orig.espera {
                Espera::Igual => assert_eq!(l, &orig.linha(), "{contexto}"),
                Espera::Exato(p, rf, alt) => assert!(
                    pos == *p && c[3].eq_ignore_ascii_case(rf) && c[4].eq_ignore_ascii_case(alt),
                    "{contexto}: {}\n  got:      {l}\n  expected: POS {p} REF {rf} ALT {alt}",
                    orig.linha()
                ),
                Espera::MesmoHaplotipo => {
                    let g = if c[0] == "1" { &g1 } else { &g2 };
                    assert_eq!(aplica(g, orig.pos, &orig.rf, &orig.alt, contexto), aplica(g, pos, c[3], c[4], contexto), "{contexto}: {l}");
                    mesmo_hap += 1;
                    esp.aparados += usize::from(c[3].len() + c[4].len() < orig.rf.len() + orig.alt.len());
                    esp.alinhados += usize::from(pos < orig.pos);
                    esp.pos_alteradas += usize::from(pos != orig.pos);
                }
            }
        }
        assert!(mesmo_hap > 100, "{contexto}");
        assert_eq!(r, &NormSummary { lotes: r.lotes, ..esp }, "{contexto}: summary");
    };

    // ---- the same result whatever the shape of the reference file
    let seqs: [(&str, &[u8]); 2] = [("chr1", &g1), ("2", &g2)];
    let fasta = amb.p("ref.fa");
    let fai = grava_fasta(&fasta, &seqs, 60, "\n", s);
    let saida = amb.p("al_out.vcf");
    let mut referencia_lida = Vec::new();
    for (forma, caminho) in [
        ("plain, no .fai", fasta.clone()),
        ("plain, with .fai", {
            std::fs::write(format!("{fasta}.fai"), &fai).unwrap();
            fasta.clone()
        }),
        ("CRLF, 70 per line, no .fai", {
            let c = amb.p("ref_crlf.fa");
            grava_fasta(&c, &seqs, 70, "\r\n", s);
            c
        }),
        ("lines of random length", {
            let c = amb.p("ref_irregular.fa");
            grava_fasta(&c, &seqs, 0, "\n", s);
            c
        }),
        ("BGZF", {
            let c = amb.p("ref.fa.gz");
            grava_fasta(&c, &seqs, 60, "\n", s);
            c
        }),
        ("a wrong .fai older than the FASTA, ignored", {
            let c = amb.p("ref_old_fai.fa");
            grava_fasta(&c, &seqs, 80, "\n", s);
            std::fs::write(format!("{c}.fai"), &fai).unwrap(); // the .fai of the 60-per-line file
            let antes = std::time::SystemTime::now() - std::time::Duration::from_secs(3600);
            std::fs::File::options().write(true).open(format!("{c}.fai")).unwrap().set_modified(antes).unwrap();
            c
        }),
    ] {
        let r = vcf_normalize(&entrada, &saida, true, Some(&caminho), &amb.op(usize::MAX, OrdemCromossomos::Cabecalho))
            .unwrap_or_else(|e| panic!("{forma}: {e}"));
        confere(&saida, &r, forma);
        assert!(amb.sem_temporarios(), "{forma}: a temporary file was left behind");
        referencia_lida.push(dados(&le(&saida)));
        ensaios += 1;
    }
    assert!(referencia_lida.windows(2).all(|w| w[0] == w[1]));
    assert!(esp.alinhados > 300 && esp.aparados > 300, "{esp:?}");

    // ---- the smallest memory budget, and a compressed output
    let r = vcf_normalize(&entrada, &amb.p("al_out.vcf.gz"), true, Some(&fasta), &amb.op(1, OrdemCromossomos::Cabecalho)).unwrap();
    confere(&amb.p("al_out.vcf.gz"), &r, "smallest budget, .gz output");
    ensaios += 1;

    // ---- aligning what is aligned changes nothing
    let r2 = vcf_normalize(&saida, &amb.p("al_again.vcf"), true, Some(&fasta), &amb.op(usize::MAX, OrdemCromossomos::Cabecalho)).unwrap();
    assert_eq!(
        (r2.aparados, r2.alinhados, r2.pos_alteradas, r2.reordenado, r2.ref_divergente, r2.fora_da_referencia),
        (0, 0, 0, false, esp.ref_divergente, esp.fora_da_referencia),
        "second pass"
    );
    iguais(&dados(&le(&amb.p("al_again.vcf"))), &referencia_lida[0], "second pass");
    ensaios += 1;

    // ---- by hand: a run of A, the first position of a sequence, a split site, lower case
    {
        //          1234567890123456789
        let g = b"AAAACGTAAAAACGTCACACAGGT";
        let fa_um = amb.p("mao_um.fa");
        std::fs::write(&fa_um, format!(">s\n{}\n", txt(g))).unwrap();
        let vcf = [
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "s\t11\ta\tAA\tA\t.\t.\t.",        // one A of AAAAA (8-12) deleted, written at the right end
            "s\t2\tb\tAA\tA\t.\t.\t.",         // one A of AAAA (1-4): it cannot go before position 1
            "s\t20\tc\tC\tCAC\t.\t.\t.",       // AC inserted in CACACA (16-21)
            "s\t11\td\taa\ta,aaa\t.\t.\t.",     // a split site, in lower case
            "s\t9\te\tA\tG\t.\t.\t.",          // a SNP: untouched
        ];
        std::fs::write(amb.p("mao_al.vcf"), vcf.join("\n")).unwrap();
        let r = vcf_normalize(&amb.p("mao_al.vcf"), &amb.p("mao_al_out.vcf"), true, Some(&fa_um), &amb.op(usize::MAX, OrdemCromossomos::Natural)).unwrap();
        let esp = [
            "s\t1\tb\tAA\tA\t.\t.\t.",
            "s\t7\ta\tTA\tT\t.\t.\t.",
            "s\t7\td\tTA\tT\t.\t.\t.",
            "s\t7\td\tT\tTA\t.\t.\t.",
            "s\t9\te\tA\tG\t.\t.\t.",
            "s\t15\tc\tT\tTCA\t.\t.\t.",
        ];
        iguais(&dados(&le(&amb.p("mao_al_out.vcf"))), &esp.map(str::to_string), "by hand");
        assert_eq!((r.lidos, r.gravados, r.separados, r.alinhados, r.pos_alteradas), (5, 6, 1, 5, 5), "{r:?}");
        // the multiallelic site is not aligned when it is not split
        let r = vcf_normalize(&amb.p("mao_al.vcf"), &amb.p("mao_al_out.vcf"), false, Some(&fa_um), &amb.op(usize::MAX, OrdemCromossomos::Natural)).unwrap();
        assert!(dados(&le(&amb.p("mao_al_out.vcf"))).contains(&vcf[4].to_string()));
        assert_eq!((r.gravados, r.alinhados), (5, 3), "{r:?}");
        ensaios += 1;

        // errors: no such file, not a FASTA, a name twice, the reference is never written
        let e = vcf_normalize(&amb.p("mao_al.vcf"), &amb.p("x.vcf"), true, Some(&amb.p("nope.fa")), &amb.op(1, OrdemCromossomos::Natural)).unwrap_err();
        assert!(e.contains("nope.fa"), "{e}");
        let e = vcf_normalize(&amb.p("mao_al.vcf"), &amb.p("x.vcf"), true, Some(&amb.p("mao_al.vcf")), &amb.op(1, OrdemCromossomos::Natural)).unwrap_err();
        assert!(e.contains("not a FASTA") || e.contains("no `>` sequence"), "{e}");
        std::fs::write(amb.p("dup.fa"), ">s\nACGT\n>s\nACGT\n").unwrap();
        let e = vcf_normalize(&amb.p("mao_al.vcf"), &amb.p("x.vcf"), true, Some(&amb.p("dup.fa")), &amb.op(1, OrdemCromossomos::Natural)).unwrap_err();
        assert!(e.contains("appears twice"), "{e}");
        assert!(!std::path::Path::new(&amb.p("x.vcf")).exists() && !std::path::Path::new(&format!("{fa_um}.fai")).exists());
        assert!(amb.sem_temporarios());
        ensaios += 1;
    }
    ensaios
}

fn main() {
    let dir = std::env::temp_dir().join(format!("test_vcf_tools_{}", std::process::id()));
    let amb = Amb { runs: dir.join("runs"), dir };
    std::fs::create_dir_all(&amb.runs).unwrap();
    let mut s = 7u64;

    let mut ensaios = testa_stats(&amb, &mut s);
    ensaios += testa_filtro(&amb, &mut s);
    ensaios += testa_normaliza(&amb, &mut s);
    ensaios += testa_alinhamento(&amb, &mut s);

    assert!(amb.sem_temporarios(), "a temporary file was left behind");
    println!("{ensaios} VCF tool tests passed");

    // real files, if given: cargo run --release --bin vcf_tools -- file.vcf.gz
    for arq in std::env::args().skip(1) {
        arquivo_real(&amb, &arq);
    }
    let _ = std::fs::remove_dir_all(&amb.dir);
}
