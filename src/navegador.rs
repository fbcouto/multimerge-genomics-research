// =====================================================================
//  NAVEGADOR DE GENOMA: OS DADOS DE UMA VISTA COMO A DO IGV
//
//  A interface desenha (em SVG, no Python) as leituras de um BAM sobre um
//  trecho do genoma; aqui so se produzem os dados, em JSON. E a vista que um
//  pesquisador abre para decidir se uma variante candidata e real ou um
//  artefato: diferencas so na ponta das leituras, so numa fita, em bases
//  recortadas (soft clip) ou de baixa qualidade.
//
//   - `region_reads`    as leituras de uma janela: blocos alinhados (os
//                       intervalos entre eles mostram delecoes e introns),
//                       diferencas para a referencia, insercoes, delecoes,
//                       recortes, o par, a linha em que cada uma e desenhada
//                       (empacotamento guloso por inicio, como o IGV); e a
//                       cobertura por posicao, contada com as MESMAS regras do
//                       `clinico::pileup` (que bate com o `samtools mpileup`);
//   - `depth_profile`   profundidade media por faixa (bin) em varias regioes,
//                       para graficos de numero de copias (CYP2D6 x CYP2D7
//                       contra regioes de controle);
//   - `fetch_sequence`  as bases de um trecho do genoma de referencia.
//
//  Regras da casa: BAM com `.bai` e lido pelo indice (so os blocos da
//  janela); sem indice, ou SAM, o arquivo e lido inteiro, e o resultado diz
//  isso. `chr1` x `1` (e `chrM` x `MT`) casam sozinhos, no BAM e no FASTA.
//  Toda mensagem diz o que esta errado, em que arquivo, e como consertar.
//
//  Regras da cobertura (as do `pileup`, conferidas posicao a posicao no
//  ensaio): leitura sem posicao, sem CIGAR ou com a FLAG 0x4 nao conta (o
//  htslib tambem a descarta antes de qualquer filtro); a FLAG e a qualidade
//  de mapeamento filtram a leitura; uma base abaixo de `min_baseq` fica de
//  fora; uma delecao e filtrada pela qualidade da base que vem depois dela;
//  QUAL `*` passa em qualquer filtro; SEQ `*` vale qualidade 0; um salto `N`
//  (intron) nao entra na profundidade; a coluna `ins` conta as leituras com
//  uma insercao logo depois da posicao (ou depois de um `P`), como o pileup.
//  Mates sobrepostos contam os dois.
// =====================================================================

use crate::bam;
use crate::bgzf::BgzfReader;
use crate::externo::abre_leitura;
use crate::indice;
use crate::nomes::{canonico, json_texto};
use crate::referencia::Genoma;
use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeSet, BinaryHeap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Largest window `region_reads` draws (bases). Wider: `depth_profile`.
pub const MAX_WINDOW: u64 = 100_000;
/// Longest stretch `fetch_sequence` returns (bases).
pub const MAX_SEQUENCE: u64 = 10_000_000;
/// Most bins `depth_profile` makes for one region.
pub const MAX_BINS: u64 = 2_000_000;
/// Reads left out by `depth_profile`: unmapped, secondary, QC fail and
/// duplicate (0x704), as `samtools depth` and mosdepth. Supplementary
/// alignments count (the pieces of a split read are real coverage).
pub const DEPTH_EXCLUDE_FLAGS: u16 = 0x704;

// ---------------------------------------------------------------------
//  OPCOES E JSON
// ---------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct ReadsOptions {
    /// FASTA the reads were aligned to. Given: mismatches are computed
    /// against it (and the window's bases are returned). Not given: from the
    /// MD tag of each read, when it has one.
    pub reference: Option<String>,
    /// Most reads returned for drawing (default 500). More than that: a fixed
    /// pseudo-random sample (by a hash of the read name, so both mates of a
    /// pair stay together and the same call gives the same reads). The
    /// coverage always counts every read that passes the filters.
    pub max_reads: usize,
    /// Reads with a lower mapping quality are left out (default 0).
    pub min_mapq: u8,
    /// Bases with a lower base quality are left out of the coverage
    /// (default 13, as samtools mpileup). The reads still show them.
    pub min_baseq: u8,
    /// Reads with any of these FLAG bits are left out. Default 0xF04:
    /// unmapped, secondary, QC fail, duplicate, supplementary.
    pub exclude_flags: u16,
    /// For a region given as one position ("chr1:1000"): bases on each side
    /// (default 60).
    pub flank: u64,
}

impl Default for ReadsOptions {
    fn default() -> Self {
        ReadsOptions { reference: None, max_reads: 500, min_mapq: 0, min_baseq: 13, exclude_flags: 0xF04, flank: 60 }
    }
}

fn jopt(s: Option<&str>) -> String {
    s.map(json_texto).unwrap_or_else(|| "null".to_string())
}

fn jlista(v: &[String]) -> String {
    format!("[{}]", v.iter().map(|s| json_texto(s)).collect::<Vec<_>>().join(","))
}

/// A number with at most `casas` decimals, as JSON.
fn jdec(x: f64, casas: i32) -> String {
    if !x.is_finite() {
        return "null".to_string();
    }
    let f = 10f64.powi(casas);
    let y = (x * f).round() / f;
    // -0 -> 0
    format!("{}", if y == 0.0 { 0.0 } else { y })
}

/// 1234567 -> "1,234,567".
fn milhar(n: u64) -> String {
    let s = n.to_string();
    let primeiro = match s.len() % 3 {
        0 => 3,
        r => r,
    };
    let mut o = s[..primeiro.min(s.len())].to_string();
    for grupo in s.as_bytes()[primeiro.min(s.len())..].chunks(3) {
        o.push(',');
        o.push_str(std::str::from_utf8(grupo).unwrap_or(""));
    }
    o
}

fn nome_curto(caminho: &str) -> &str {
    caminho.rsplit(['/', '\\']).next().unwrap_or(caminho)
}

// ---------------------------------------------------------------------
//  O ARQUIVO DE ALINHAMENTOS
// ---------------------------------------------------------------------

/// O `.bai` de um BAM: `x.bam.bai` ou `x.bai`, como o samtools aceita.
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

struct Fonte {
    caminho: String,
    e_bam: bool,
    bai: Option<String>,
    /// As sequencias do cabecalho: nome e comprimento (0 = sem LN).
    seqs: Vec<(String, u64)>,
}

impl Fonte {
    fn abre(caminho: &str) -> Result<Fonte, String> {
        if !Path::new(caminho).is_file() {
            return Err(format!("{caminho} does not exist (or is not a file). Give the path of a BAM, SAM or .sam.gz file of aligned reads."));
        }
        let baixo = caminho.to_ascii_lowercase();
        if baixo.ends_with(".cram") {
            return Err(format!(
                "{caminho}: CRAM is not read here. Convert it to BAM first, e.g. `samtools view -b -T genome.fa -o reads.bam {}` and `samtools index reads.bam`.",
                nome_curto(caminho)
            ));
        }
        if baixo.ends_with(".bam") {
            let arq = File::open(caminho).map_err(|e| format!("Error opening {caminho}: {e}"))?;
            let mut r = BgzfReader::novo(BufReader::new(arq));
            let cab = bam::le_cabecalho(&mut r).map_err(|e| format!("{caminho}: {e}"))?;
            return Ok(Fonte {
                caminho: caminho.to_string(),
                e_bam: true,
                bai: acha_bai(caminho),
                seqs: cab.refs.iter().map(|(n, l)| (n.clone(), u64::from(*l))).collect(),
            });
        }
        // SAM em texto (ou .sam.gz): as linhas @SQ do topo.
        let mut r = abre_leitura(caminho)?;
        let mut seqs = Vec::new();
        let mut linha = String::new();
        loop {
            linha.clear();
            let n = r.read_line(&mut linha).map_err(|e| erro_leitura(caminho, &e))?;
            if n == 0 {
                break;
            }
            let t = linha.trim_end_matches(['\n', '\r']);
            if t.is_empty() {
                continue;
            }
            if !t.starts_with('@') {
                break;
            }
            if let Some(resto) = t.strip_prefix("@SQ\t") {
                let sn = resto.split('\t').find_map(|c| c.strip_prefix("SN:"));
                let ln = resto.split('\t').find_map(|c| c.strip_prefix("LN:")).and_then(|v| v.parse::<u64>().ok());
                if let Some(sn) = sn {
                    seqs.push((sn.to_string(), ln.unwrap_or(0)));
                }
            }
        }
        Ok(Fonte { caminho: caminho.to_string(), e_bam: false, bai: None, seqs })
    }

    fn indexado(&self) -> bool {
        self.e_bam && self.bai.is_some()
    }

    /// O nome no arquivo da sequencia `pedido` (igual, ou o mesmo
    /// cromossomo na outra convencao) e o comprimento (0 = desconhecido).
    /// `None`: o arquivo nao declara sequencias (SAM sem @SQ).
    fn resolve(&self, pedido: &str) -> Result<Option<(String, u64)>, String> {
        if self.seqs.is_empty() {
            return Ok(None);
        }
        if let Some(s) = self.seqs.iter().find(|s| s.0 == pedido) {
            return Ok(Some(s.clone()));
        }
        let c = canonico(pedido);
        match self.seqs.iter().find(|s| canonico(&s.0) == c) {
            Some(s) => Ok(Some(s.clone())),
            None => {
                let mut v: Vec<&str> = self.seqs.iter().take(8).map(|s| s.0.as_str()).collect();
                if self.seqs.len() > 8 {
                    v.push("...");
                }
                Err(format!(
                    "Sequence \"{pedido}\" is not in {}. Its sequences are named like: {}. Check the name, or the genome build of the file.",
                    self.caminho,
                    v.join(", ")
                ))
            }
        }
    }

    /// Algum nome do arquivo e exatamente `t` (ou o mesmo cromossomo)?
    fn conhece(&self, t: &str) -> bool {
        self.seqs.iter().any(|s| s.0 == t)
    }

    /// Avisos sobre o modo de leitura.
    fn avisos(&self, avisos: &mut Vec<String>) {
        if self.e_bam {
            match &self.bai {
                None => avisos.push(format!(
                    "No index (.bai) next to {}: the whole file was read to find the reads of the region. Build it with index_file(\"{}\") (the BAM must be sorted by position) to make this fast.",
                    self.caminho, self.caminho
                )),
                Some(b) => {
                    let velho = match (std::fs::metadata(b).and_then(|m| m.modified()), std::fs::metadata(&self.caminho).and_then(|m| m.modified())) {
                        (Ok(ti), Ok(tb)) => ti < tb,
                        _ => false,
                    };
                    if velho {
                        avisos.push(format!(
                            "The index {b} is older than {}: if the BAM changed after it was indexed, the reads shown can be wrong. Index it again with index_file(\"{}\").",
                            self.caminho, self.caminho
                        ));
                    }
                }
            }
        } else {
            avisos.push(format!(
                "{} is SAM text, which has no index: the whole file was read. For fast views, write it as a sorted, indexed BAM: sort_genomic([\"{}\"], \"reads.bam\", index=True).",
                self.caminho, self.caminho
            ));
        }
    }

    /// Entrega cada linha de alinhamento que pode tocar `[ini, fim]` (base 1)
    /// da sequencia `nome`. Com o `.bai`, so os blocos que o indice aponta;
    /// sem ele, o arquivo inteiro, filtrado pelo nome (`exato`: o nome do
    /// cabecalho; senao, o mesmo cromossomo em qualquer convencao).
    fn linhas(&self, nome: &str, exato: bool, ini: u64, fim: u64, recebe: &mut dyn FnMut(&str) -> Result<(), String>) -> Result<(), String> {
        if self.indexado() {
            if fim > indice::POS_MAX as u64 {
                return Err(format!("Position {fim} is beyond what a BAM index addresses (512 Mbp per sequence)"));
            }
            let regiao = format!("{nome}:{ini}-{fim}");
            indice::consultar(&self.caminho, &[regiao], &mut |cab, l| if cab { Ok(()) } else { recebe(l) })?;
            return Ok(());
        }
        let canon = canonico(nome);
        let mut r = abre_leitura(&self.caminho)?;
        let mut linha = String::new();
        let mut cache: (String, bool) = (String::new(), false);
        loop {
            linha.clear();
            let n = r.read_line(&mut linha).map_err(|e| erro_leitura(&self.caminho, &e))?;
            if n == 0 {
                break;
            }
            let t = linha.trim_end_matches(['\n', '\r']);
            if t.is_empty() || t.starts_with('@') {
                continue;
            }
            let rname = t.split('\t').nth(2).unwrap_or("");
            if rname != cache.0 {
                let casa = rname == nome || (!exato && canonico(rname) == canon);
                cache = (rname.to_string(), casa);
            }
            if cache.1 {
                recebe(t)?;
            }
        }
        Ok(())
    }
}

fn erro_leitura(caminho: &str, e: &std::io::Error) -> String {
    format!("Error reading {caminho}: {e}. The file may be truncated, or not SAM/BAM (CRAM is not read; convert it to BAM).")
}

// ---------------------------------------------------------------------
//  UM ALINHAMENTO
// ---------------------------------------------------------------------

/// CIGAR em texto -> operacoes. `None`: invalido. (O mesmo que o
/// `clinico::le_cigar`, que e privado.)
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

/// Ha uma insercao logo depois da operacao `i`, pela regra do pileup do
/// htslib (e do `clinico::pileup`)? A seguinte e I; ou e P e ha um I com
/// bases antes da proxima operacao que consome a referencia.
fn insercao_depois(cigar: &[(u32, u8)], i: usize) -> bool {
    match cigar.get(i + 1) {
        Some(&(_, b'I')) => true,
        Some(&(_, b'P')) => cigar[i + 2..]
            .iter()
            .take_while(|&&(_, o)| !matches!(o, b'D' | b'M' | b'N' | b'=' | b'X'))
            .any(|&(l, o)| o == b'I' && l > 0),
        _ => false,
    }
}

struct Leitura<'a> {
    nome: &'a str,
    flag: u16,
    pos: u64,
    mapq: u32,
    cigar_txt: &'a str,
    cigar: Vec<(u32, u8)>,
    rname: &'a str,
    rnext: &'a str,
    pnext: u64,
    tlen: i64,
    seq: &'a [u8],
    qual: Option<&'a [u8]>,
    tags: &'a str,
    /// Ultima base alinhada (base 1, incluida).
    fim: u64,
}

impl<'a> Leitura<'a> {
    /// `None`: a linha nao e um alinhamento com posicao (colunas de menos,
    /// FLAG ou POS nao numericos, POS 0, CIGAR `*`, FLAG 0x4).
    fn le(linha: &'a str, caminho: &str) -> Result<Option<Leitura<'a>>, String> {
        let c: Vec<&str> = linha.splitn(12, '\t').collect();
        if c.len() < 11 {
            return Ok(None);
        }
        let (Ok(flag), Ok(pos)) = (c[1].parse::<u16>(), c[3].parse::<u64>()) else { return Ok(None) };
        if c[5] == "*" || pos == 0 || flag & 4 != 0 {
            return Ok(None);
        }
        let cigar = le_cigar(c[5]).ok_or_else(|| {
            format!("{caminho}: read {} has an invalid CIGAR \"{}\". Fix or remove that line (e.g. with samtools view).", c[0], c[5])
        })?;
        let compr: u64 = cigar.iter().filter(|(_, o)| matches!(o, b'M' | b'D' | b'N' | b'=' | b'X')).map(|&(l, _)| u64::from(l)).sum();
        Ok(Some(Leitura {
            nome: c[0],
            flag,
            pos,
            mapq: c[4].parse().unwrap_or(0),
            cigar_txt: c[5],
            cigar,
            rname: c[2],
            rnext: c[6],
            pnext: c[7].parse().unwrap_or(0),
            tlen: c[8].parse().unwrap_or(0),
            seq: if c[9] == "*" { &[] } else { c[9].as_bytes() },
            qual: if c[10] == "*" { None } else { Some(c[10].as_bytes()) },
            tags: c.get(11).copied().unwrap_or(""),
            fim: pos + compr.max(1) - 1,
        }))
    }

    /// Qualidade da base `i` como o pileup a ve: SEQ `*` = 0, QUAL `*` = 255.
    fn qualidade(&self, i: usize) -> u8 {
        if self.seq.is_empty() {
            return 0;
        }
        match self.qual {
            None => 255,
            Some(q) => q.get(i).map(|c| c.saturating_sub(33)).unwrap_or(0),
        }
    }

    /// Qualidade para mostrar (`None` sem QUAL).
    fn qual_mostrada(&self, i: usize) -> Option<u8> {
        if self.seq.is_empty() {
            return None;
        }
        self.qual.and_then(|q| q.get(i)).map(|c| c.saturating_sub(33))
    }

    fn base(&self, i: usize) -> u8 {
        self.seq.get(i).map(|b| b.to_ascii_uppercase()).unwrap_or(b'N')
    }

    fn md(&self) -> Option<&'a str> {
        self.tags.split('\t').find_map(|t| t.strip_prefix("MD:Z:"))
    }

    fn sobrepoe(&self, ini: u64, fim: u64) -> bool {
        self.pos <= fim && self.fim >= ini
    }
}

fn classe(b: u8) -> usize {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

// ---------------------------------------------------------------------
//  COBERTURA (as regras do pileup)
// ---------------------------------------------------------------------

/// Por posicao: A C G T N del ins.
type Contagem = [u32; 7];

struct Cobertura {
    ini: u64,
    tot: Vec<Contagem>,
    frente: Vec<Contagem>,
}

impl Cobertura {
    fn nova(ini: u64, fim: u64) -> Cobertura {
        let n = (fim - ini + 1) as usize;
        Cobertura { ini, tot: vec![[0; 7]; n], frente: vec![[0; 7]; n] }
    }

    fn soma(&mut self, p: u64, k: usize, frente: bool) {
        let s = (p - self.ini) as usize;
        self.tot[s][k] += 1;
        if frente {
            self.frente[s][k] += 1;
        }
    }

    /// Conta o que a leitura mostra em cada posicao da janela -- o mesmo
    /// que `clinico::pileup` conta numa posicao, numa passada so pelo CIGAR.
    fn conta(&mut self, l: &Leitura, min_baseq: u8) {
        let (ini, fim) = (self.ini, self.ini + self.tot.len() as u64 - 1);
        let frente = l.flag & 0x10 == 0;
        let (mut r, mut q) = (l.pos, 0usize);
        for (i, &(len, op)) in l.cigar.iter().enumerate() {
            let n = u64::from(len);
            match op {
                b'M' | b'=' | b'X' => {
                    if n > 0 && r <= fim && r + n > ini {
                        let ultima = r + n - 1;
                        let ins = insercao_depois(&l.cigar, i);
                        for p in r.max(ini)..=ultima.min(fim) {
                            let qi = q + (p - r) as usize;
                            if l.qualidade(qi) < min_baseq {
                                continue;
                            }
                            self.soma(p, classe(l.base(qi)), frente);
                            if p == ultima && ins {
                                self.soma(p, 6, frente);
                            }
                        }
                    }
                    r += n;
                    q += len as usize;
                }
                b'D' => {
                    // A delecao vale a qualidade da base depois dela.
                    if n > 0 && r <= fim && r + n > ini && l.qualidade(q) >= min_baseq {
                        let ultima = r + n - 1;
                        let ins = insercao_depois(&l.cigar, i);
                        for p in r.max(ini)..=ultima.min(fim) {
                            self.soma(p, 5, frente);
                            if p == ultima && ins {
                                self.soma(p, 6, frente);
                            }
                        }
                    }
                    r += n;
                }
                b'N' => r += n,
                b'I' | b'S' => q += len as usize,
                _ => {}
            }
            if r > fim {
                break;
            }
        }
    }
}

// ---------------------------------------------------------------------
//  O QUE SE DESENHA DE CADA LEITURA
// ---------------------------------------------------------------------

/// Uma posicao alinhada (referencia, indice na leitura) ou uma delecao
/// (inicio, comprimento), na ordem do CIGAR: o que o MD descreve.
#[derive(Clone, Copy)]
enum Evento {
    Al(u64, usize),
    Del(u64, u32),
}

fn eventos(l: &Leitura) -> Vec<Evento> {
    let mut v = Vec::new();
    let (mut r, mut q) = (l.pos, 0usize);
    for &(len, op) in &l.cigar {
        let n = u64::from(len);
        match op {
            b'M' | b'=' | b'X' => {
                for k in 0..n {
                    v.push(Evento::Al(r + k, q + k as usize));
                }
                r += n;
                q += len as usize;
            }
            b'D' => {
                if len > 0 {
                    v.push(Evento::Del(r, len));
                }
                r += n;
            }
            b'N' => r += n,
            b'I' | b'S' => q += len as usize,
            _ => {}
        }
    }
    v
}

/// O que o MD diz dentro da janela: as posicoes alinhadas (com a base da
/// referencia quando e uma diferenca) e as bases deletadas.
#[derive(Default)]
struct Md {
    /// (posicao, indice na leitura, base da referencia se difere).
    alinhadas: Vec<(u64, usize, Option<u8>)>,
    deletadas: Vec<(u64, u8)>,
}

/// Le o MD contra o CIGAR. `None`: o MD nao descreve este CIGAR.
fn le_md(md: &str, ev: &[Evento], ini: u64, fim: u64) -> Option<Md> {
    let b = md.as_bytes();
    let (mut i, mut e) = (0usize, 0usize);
    let mut out = Md::default();
    let dentro = |p: u64| p >= ini && p <= fim;
    while i < b.len() {
        if b[i].is_ascii_digit() {
            let mut n: u64 = 0;
            while i < b.len() && b[i].is_ascii_digit() {
                n = n * 10 + u64::from(b[i] - b'0');
                if n > u64::from(u32::MAX) {
                    return None;
                }
                i += 1;
            }
            for _ in 0..n {
                match ev.get(e) {
                    Some(&Evento::Al(p, qi)) => {
                        if dentro(p) {
                            out.alinhadas.push((p, qi, None));
                        }
                    }
                    _ => return None,
                }
                e += 1;
            }
        } else if b[i] == b'^' {
            i += 1;
            let s = i;
            while i < b.len() && b[i].is_ascii_alphabetic() {
                i += 1;
            }
            // Um grupo `^` cobre uma ou mais operacoes D seguidas (o calmd
            // escreve "^AB0^CD" para 2D2D; outros, "^ABCD"). Letras que
            // sobram sao uma diferenca logo depois, sem o `0` que deveria
            // separa-las ("^AC" + "T").
            let m = i - s;
            let mut usadas = 0usize;
            while let Some(&Evento::Del(p, len)) = ev.get(e) {
                if usadas + len as usize > m {
                    break;
                }
                for k in 0..len as usize {
                    if dentro(p + k as u64) {
                        out.deletadas.push((p + k as u64, b[s + usadas + k].to_ascii_uppercase()));
                    }
                }
                usadas += len as usize;
                e += 1;
            }
            if usadas == 0 {
                return None;
            }
            i = s + usadas;
        } else if b[i].is_ascii_alphabetic() {
            match ev.get(e) {
                Some(&Evento::Al(p, qi)) => {
                    if dentro(p) {
                        out.alinhadas.push((p, qi, Some(b[i].to_ascii_uppercase())));
                    }
                }
                _ => return None,
            }
            e += 1;
            i += 1;
        } else {
            return None;
        }
    }
    if e != ev.len() {
        return None;
    }
    Some(out)
}

#[derive(Clone, Copy, PartialEq)]
enum OrigemDif {
    Referencia,
    Md,
    /// Sem referencia e sem MD (ou MD que nao casa com o CIGAR).
    Nenhuma,
}

struct Par {
    chrom: String,
    pos: Option<u64>,
    reversa: bool,
    nao_mapeado: bool,
    proprio: bool,
    tlen: i64,
    numero: Option<u8>,
}

struct Vista {
    nome: String,
    flag: u16,
    mapq: u32,
    inicio: u64,
    fim: u64,
    cigar: String,
    /// soft esquerdo, soft direito, hard esquerdo, hard direito.
    recortes: [u32; 4],
    seq_recortes: (String, String),
    blocos: Vec<(u64, u64)>,
    difs: Vec<(u64, u8, Option<u8>)>,
    insercoes: Vec<(u64, u32, String)>,
    delecoes: Vec<(u64, u32)>,
    saltos: Vec<(u64, u32)>,
    par: Option<Par>,
    media_q: Option<f64>,
    tem_seq: bool,
    origem: OrigemDif,
    tem_md: bool,
    md_invalido: bool,
    /// Posicoes da janela em que o MD e o FASTA discordam.
    discordancias: usize,
    /// Bases da referencia que o MD revela na janela (sem FASTA).
    ref_md: Vec<(u64, u8)>,
    linha: usize,
}

impl Vista {
    fn json(&self) -> String {
        let blocos: Vec<String> = self.blocos.iter().map(|(a, b)| format!("[{a},{b}]")).collect();
        let difs: Vec<String> = self
            .difs
            .iter()
            .map(|(p, b, q)| format!("[{p},{},{}]", json_texto(&(*b as char).to_string()), q.map(|q| q.to_string()).unwrap_or_else(|| "null".into())))
            .collect();
        let ins: Vec<String> = self.insercoes.iter().map(|(p, n, s)| format!("[{p},{n},{}]", json_texto(s))).collect();
        let dels: Vec<String> = self.delecoes.iter().map(|(p, n)| format!("[{p},{n}]")).collect();
        let saltos: Vec<String> = self.saltos.iter().map(|(p, n)| format!("[{p},{n}]")).collect();
        let par = match &self.par {
            None => "null".to_string(),
            Some(m) => format!(
                "{{\"chrom\":{},\"start\":{},\"strand\":\"{}\",\"unmapped\":{},\"proper_pair\":{},\"insert_size\":{},\"read_number\":{}}}",
                json_texto(&m.chrom),
                m.pos.map(|p| p.to_string()).unwrap_or_else(|| "null".into()),
                if m.reversa { '-' } else { '+' },
                m.nao_mapeado,
                m.proprio,
                m.tlen,
                m.numero.map(|n| n.to_string()).unwrap_or_else(|| "null".into())
            ),
        };
        format!(
            "{{\"name\":{},\"flag\":{},\"strand\":\"{}\",\"mapq\":{},\"start\":{},\"end\":{},\"row\":{},\"cigar\":{},\"clip_left\":{},\"clip_right\":{},\"hard_clip_left\":{},\"hard_clip_right\":{},\"clip_left_seq\":{},\"clip_right_seq\":{},\"blocks\":[{}],\"mismatches\":[{}],\"insertions\":[{}],\"deletions\":[{}],\"skips\":[{}],\"mate\":{},\"mean_qual\":{},\"has_sequence\":{},\"mismatch_source\":\"{}\"}}",
            json_texto(&self.nome),
            self.flag,
            if self.flag & 0x10 != 0 { '-' } else { '+' },
            self.mapq,
            self.inicio,
            self.fim,
            self.linha,
            json_texto(&self.cigar),
            self.recortes[0],
            self.recortes[1],
            self.recortes[2],
            self.recortes[3],
            json_texto(&self.seq_recortes.0),
            json_texto(&self.seq_recortes.1),
            blocos.join(","),
            difs.join(","),
            ins.join(","),
            dels.join(","),
            saltos.join(","),
            par,
            self.media_q.map(|m| jdec(m, 2)).unwrap_or_else(|| "null".into()),
            self.tem_seq,
            match self.origem {
                OrigemDif::Referencia => "reference",
                OrigemDif::Md => "MD",
                OrigemDif::Nenhuma => "none",
            }
        )
    }
}

/// `refj`: as bases da janela no FASTA (posicao `ini` em `refj[0]`), quando
/// ha FASTA com esta sequencia.
fn monta_vista(l: &Leitura, ini: u64, fim: u64, refj: Option<&[u8]>) -> Vista {
    let ops = &l.cigar;
    // Recortes nas pontas.
    let mut recortes = [0u32; 4];
    for &(n, o) in ops.iter() {
        match o {
            b'H' => recortes[2] += n,
            b'S' => recortes[0] += n,
            _ => break,
        }
    }
    for &(n, o) in ops.iter().rev() {
        match o {
            b'H' => recortes[3] += n,
            b'S' => recortes[1] += n,
            _ => break,
        }
    }
    let tem_seq = !l.seq.is_empty();
    let pedaco = |a: usize, n: usize| -> String { (a..a + n).map(|j| l.base(j) as char).collect() };
    let seq_recortes = if tem_seq {
        let (e, d) = (recortes[0] as usize, recortes[1] as usize);
        (pedaco(0, e.min(l.seq.len())), pedaco(l.seq.len().saturating_sub(d), d.min(l.seq.len())))
    } else {
        (String::new(), String::new())
    };

    // Blocos, insercoes, delecoes e saltos, num passo pelo CIGAR.
    let mut blocos: Vec<(u64, u64)> = Vec::new();
    let mut insercoes: Vec<(u64, u32, String)> = Vec::new();
    let mut delecoes: Vec<(u64, u32)> = Vec::new();
    let mut saltos: Vec<(u64, u32)> = Vec::new();
    let (mut r, mut q) = (l.pos, 0usize);
    // Desde a ultima insercao / delecao, nada consumiu a referencia (ou so P).
    let (mut junta_ins, mut junta_del) = (false, false);
    for &(len, op) in ops {
        let n = u64::from(len);
        match op {
            b'M' | b'=' | b'X' => {
                if n > 0 {
                    match blocos.last_mut() {
                        Some(b) if b.1 + 1 == r => b.1 = r + n - 1,
                        _ => blocos.push((r, r + n - 1)),
                    }
                    junta_ins = false;
                    junta_del = false;
                }
                r += n;
                q += len as usize;
            }
            b'I' => {
                if len > 0 {
                    let s = if tem_seq { pedaco(q, len as usize) } else { "N".repeat(len as usize) };
                    match insercoes.last_mut() {
                        Some(x) if junta_ins && x.0 == r - 1 => {
                            x.1 += len;
                            x.2.push_str(&s);
                        }
                        _ => insercoes.push((r - 1, len, s)),
                    }
                    junta_ins = true;
                    junta_del = false;
                }
                q += len as usize;
            }
            b'D' => {
                if n > 0 {
                    match delecoes.last_mut() {
                        Some(x) if junta_del && x.0 + u64::from(x.1) == r => x.1 += len,
                        _ => delecoes.push((r, len)),
                    }
                    junta_del = true;
                    junta_ins = false;
                }
                r += n;
            }
            b'N' => {
                if n > 0 {
                    saltos.push((r, len));
                    junta_del = false;
                    junta_ins = false;
                }
                r += n;
            }
            b'S' => q += len as usize,
            _ => {}
        }
    }
    insercoes.retain(|x| x.0 + 1 >= ini && x.0 <= fim);
    delecoes.retain(|x| x.0 <= fim && x.0 + u64::from(x.1) > ini);
    saltos.retain(|x| x.0 <= fim && x.0 + u64::from(x.1) > ini);

    // Diferencas: contra o FASTA, ou pelo MD.
    let md_txt = l.md();
    let md = md_txt.and_then(|m| le_md(m, &eventos(l), ini, fim));
    let md_invalido = md_txt.is_some() && md.is_none();
    let mut difs: Vec<(u64, u8, Option<u8>)> = Vec::new();
    let mut discordancias = 0usize;
    let mut ref_md: Vec<(u64, u8)> = Vec::new();
    let origem;
    if let Some(rj) = refj {
        origem = OrigemDif::Referencia;
        if tem_seq {
            let (mut r, mut q) = (l.pos, 0usize);
            for &(len, op) in ops {
                let n = u64::from(len);
                match op {
                    b'M' | b'=' | b'X' => {
                        if n > 0 && r <= fim && r + n > ini {
                            for p in r.max(ini)..=(r + n - 1).min(fim) {
                                let qi = q + (p - r) as usize;
                                let b = l.base(qi);
                                let rb = rj.get((p - ini) as usize).copied().unwrap_or(b'N');
                                if b != b'=' && b != rb {
                                    difs.push((p, b, l.qual_mostrada(qi)));
                                }
                            }
                        }
                        r += n;
                        q += len as usize;
                    }
                    b'D' | b'N' => r += n,
                    b'I' | b'S' => q += len as usize,
                    _ => {}
                }
                if r > fim {
                    break;
                }
            }
            // O MD, quando ha, tem de dizer o mesmo (onde o FASTA tem A/C/G/T).
            if let Some(m) = &md {
                for &(p, qi, rb) in &m.alinhadas {
                    let fb = rj.get((p - ini) as usize).copied().unwrap_or(b'N');
                    if !matches!(fb, b'A' | b'C' | b'G' | b'T') {
                        continue;
                    }
                    let b = l.base(qi);
                    let pelo_md = rb.is_some();
                    let pelo_fasta = b != b'=' && b != fb;
                    if pelo_md != pelo_fasta || rb.is_some_and(|x| x != fb) {
                        discordancias += 1;
                    }
                }
                for &(p, rb) in &m.deletadas {
                    let fb = rj.get((p - ini) as usize).copied().unwrap_or(b'N');
                    if matches!(fb, b'A' | b'C' | b'G' | b'T') && fb != rb {
                        discordancias += 1;
                    }
                }
            }
        }
    } else if let Some(m) = &md {
        origem = OrigemDif::Md;
        for &(p, qi, rb) in &m.alinhadas {
            match rb {
                Some(x) => {
                    if tem_seq {
                        difs.push((p, l.base(qi), l.qual_mostrada(qi)));
                    }
                    ref_md.push((p, x));
                }
                None => {
                    let b = l.base(qi);
                    if tem_seq && matches!(b, b'A' | b'C' | b'G' | b'T') {
                        ref_md.push((p, b));
                    }
                }
            }
        }
        ref_md.extend(m.deletadas.iter().copied());
    } else {
        origem = OrigemDif::Nenhuma;
    }

    let par = if l.flag & 1 != 0 {
        Some(Par {
            chrom: match l.rnext {
                "=" => l.rname.to_string(),
                "*" => String::new(),
                x => x.to_string(),
            },
            pos: if l.pnext > 0 { Some(l.pnext) } else { None },
            reversa: l.flag & 0x20 != 0,
            nao_mapeado: l.flag & 0x8 != 0,
            proprio: l.flag & 0x2 != 0,
            tlen: l.tlen,
            numero: match (l.flag & 0x40 != 0, l.flag & 0x80 != 0) {
                (true, false) => Some(1),
                (false, true) => Some(2),
                _ => None,
            },
        })
    } else {
        None
    };
    let media_q = match l.qual {
        Some(q) if tem_seq && !q.is_empty() => Some(q.iter().map(|&c| f64::from(c.saturating_sub(33))).sum::<f64>() / q.len() as f64),
        _ => None,
    };
    Vista {
        nome: l.nome.to_string(),
        flag: l.flag,
        mapq: l.mapq,
        inicio: l.pos,
        fim: l.fim,
        cigar: l.cigar_txt.to_string(),
        recortes,
        seq_recortes,
        blocos,
        difs,
        insercoes,
        delecoes,
        saltos,
        par,
        media_q,
        tem_seq,
        origem,
        tem_md: md_txt.is_some(),
        md_invalido,
        discordancias,
        ref_md,
        linha: 0,
    }
}

// ---------------------------------------------------------------------
//  AMOSTRA DETERMINISTICA E EMPACOTAMENTO
// ---------------------------------------------------------------------

/// FNV-1a + o finalizador do splitmix64: o mesmo numero em toda maquina.
fn espalha(s: &str) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for b in s.bytes() {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100000001b3);
    }
    h = (h ^ (h >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    h = (h ^ (h >> 27)).wrapping_mul(0x94d049bb133111eb);
    h ^ (h >> 31)
}

/// Chave da amostra: ficam as `max_reads` menores. O nome primeiro, para os
/// dois mates ficarem juntos; o numero de chegada desempata linhas iguais.
type Chave = (u64, String, u16, u64, u64);

struct Candidato {
    chave: Chave,
    vista: Vista,
}

impl PartialEq for Candidato {
    fn eq(&self, o: &Self) -> bool {
        self.chave == o.chave
    }
}
impl Eq for Candidato {}
impl PartialOrd for Candidato {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Candidato {
    fn cmp(&self, o: &Self) -> Ordering {
        self.chave.cmp(&o.chave)
    }
}

/// Distribui as leituras em linhas sem sobreposicao, como o IGV: em ordem
/// de inicio, cada uma vai para a linha livre de menor numero; uma linha esta
/// livre quando a ultima leitura dela termina pelo menos `1 + folga` bases
/// antes. Da o menor numero de linhas possivel (o guloso por inicio e otimo
/// para intervalos). `vistas` deve estar em ordem de inicio. Devolve o numero
/// de linhas.
fn empacota(vistas: &mut [Vista], folga: u64) -> usize {
    let mut ocupadas: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new();
    let mut livres: BTreeSet<usize> = BTreeSet::new();
    let mut linhas = 0usize;
    for v in vistas.iter_mut() {
        while let Some(&Reverse((fim, linha))) = ocupadas.peek() {
            if fim + folga < v.inicio {
                ocupadas.pop();
                livres.insert(linha);
            } else {
                break;
            }
        }
        let linha = match livres.pop_first() {
            Some(l) => l,
            None => {
                linhas += 1;
                linhas - 1
            }
        };
        v.linha = linha;
        ocupadas.push(Reverse((v.fim, linha)));
    }
    linhas
}

// ---------------------------------------------------------------------
//  A JANELA
// ---------------------------------------------------------------------

struct Janela {
    /// O nome na fonte (ou o pedido, num SAM sem @SQ).
    nome: String,
    /// O nome casa exatamente (ha @SQ) ou pela forma canonica.
    exato: bool,
    ini: u64,
    fim: u64,
    foco: Option<u64>,
    comprimento: Option<u64>,
}

const FORMAS_REGIAO: &str = "Use \"chr1:1000-1200\" (1-based, end included) or \"chr1:1000\" (a window around that position)";

fn janela(texto: &str, fonte: &Fonte, flanco: u64) -> Result<Janela, String> {
    let t = texto.trim();
    let (pedido, ini0, fim) = indice::le_regiao(t, &|n| fonte.conhece(n))?;
    if fonte.conhece(t) || !t.contains(':') {
        return Err(format!("Region \"{texto}\" names a whole sequence; the read view needs a stretch of it. {FORMAS_REGIAO}."));
    }
    let faixa = &t[pedido.len() + 1..];
    let (nome, exato, comprimento) = match fonte.resolve(&pedido)? {
        Some((n, l)) => (n, true, if l > 0 { Some(l) } else { None }),
        None => (pedido.clone(), false, None),
    };
    let (mut ini, mut fim, foco) = if faixa.contains('-') {
        (ini0 as u64 + 1, fim as u64, None)
    } else {
        let p = ini0 as u64 + 1;
        (p.saturating_sub(flanco).max(1), p.saturating_add(flanco), Some(p))
    };
    if let Some(len) = comprimento {
        if ini > len {
            return Err(format!(
                "Region \"{texto}\" starts after the end of {nome}, which has {} bases in {}.",
                milhar(len),
                fonte.caminho
            ));
        }
        fim = fim.min(len);
    } else if fim >= indice::POS_MAX as u64 {
        return Err(format!("Region \"{texto}\" has no end, and {} does not declare the length of {nome}. {FORMAS_REGIAO}.", fonte.caminho));
    }
    ini = ini.max(1);
    let largura = fim - ini + 1;
    if largura > MAX_WINDOW {
        return Err(format!(
            "Region \"{texto}\" spans {} bases; the read view shows at most {} at a time. Choose a narrower region, or use depth_profile for the coverage of large regions.",
            milhar(largura),
            milhar(MAX_WINDOW)
        ));
    }
    Ok(Janela { nome, exato, ini, fim, foco, comprimento })
}

/// As reads of a BAM (or SAM) over a window of the genome, ready to draw as
/// an IGV-like view; JSON. See the module header for the coverage rules.
///
/// `region`: "chr1:1000-1200" (1-based, end included), or "chr1:1000" for a
/// window of `opt.flank` bases on each side of that position (`focus`).
pub fn region_reads(alignments: &str, region: &str, opt: &ReadsOptions) -> Result<String, String> {
    let fonte = Fonte::abre(alignments)?;
    let j = janela(region, &fonte, opt.flank)?;
    let (ini, fim) = (j.ini, j.fim);
    let mut avisos: Vec<String> = Vec::new();
    fonte.avisos(&mut avisos);

    // O genoma de referencia da janela.
    let mut ref_fasta: Option<Vec<u8>> = None;
    if let Some(fa) = opt.reference.as_deref().filter(|f| !f.trim().is_empty()) {
        let mut g = Genoma::abre(fa, None)?;
        match g.nome(&j.nome) {
            None => avisos.push(format!(
                "Sequence {} is not in the reference {fa} (also tried it with and without \"chr\"): mismatches are taken from the MD tags of the reads instead. Check that the FASTA is the genome the reads were aligned to.",
                j.nome
            )),
            Some(nf) => {
                let len_fa = g.comprimento(&nf).unwrap_or(0);
                if let Some(len_bam) = j.comprimento {
                    if len_fa != len_bam {
                        avisos.push(format!(
                            "In the reference {fa}, {nf} has {} bases; the header of {} says {} has {}. They are probably different builds of the genome: mismatches may be wrong. Use the FASTA the reads were aligned to.",
                            milhar(len_fa),
                            fonte.caminho,
                            j.nome,
                            milhar(len_bam)
                        ));
                    }
                }
                let mut bases = g.trecho(&nf, ini - 1, fim)?.to_vec();
                if (bases.len() as u64) < fim - ini + 1 {
                    avisos.push(format!(
                        "The window goes past the end of {nf} in the reference {fa} ({} bases): the bases beyond it are shown as N.",
                        milhar(len_fa)
                    ));
                    bases.resize((fim - ini + 1) as usize, b'N');
                }
                ref_fasta = Some(bases);
            }
        }
    }

    // Uma passada pelas leituras da janela.
    let mut cob = Cobertura::nova(ini, fim);
    let mut amostra: BinaryHeap<Candidato> = BinaryHeap::new();
    let (mut total, mut f_flag, mut f_mapq, mut chegada) = (0u64, 0u64, 0u64, 0u64);
    let mut nome_achado: Option<String> = None;
    let caminho = fonte.caminho.clone();
    fonte.linhas(&j.nome, j.exato, ini, fim, &mut |linha| {
        let Some(l) = Leitura::le(linha, &caminho)? else { return Ok(()) };
        if !l.sobrepoe(ini, fim) {
            return Ok(());
        }
        if nome_achado.is_none() {
            nome_achado = Some(l.rname.to_string());
        }
        if l.flag & opt.exclude_flags != 0 {
            f_flag += 1;
            return Ok(());
        }
        if l.mapq < u32::from(opt.min_mapq) {
            f_mapq += 1;
            return Ok(());
        }
        total += 1;
        chegada += 1;
        cob.conta(&l, opt.min_baseq);
        if opt.max_reads == 0 {
            return Ok(());
        }
        let chave: Chave = (espalha(l.nome), l.nome.to_string(), l.flag, l.pos, chegada);
        if amostra.len() >= opt.max_reads {
            if let Some(maior) = amostra.peek() {
                if chave >= maior.chave {
                    return Ok(());
                }
            }
        }
        let vista = monta_vista(&l, ini, fim, ref_fasta.as_deref());
        amostra.push(Candidato { chave, vista });
        if amostra.len() > opt.max_reads {
            amostra.pop();
        }
        Ok(())
    })?;
    let nome_final = if j.exato { j.nome.clone() } else { nome_achado.unwrap_or_else(|| j.nome.clone()) };

    let mut vistas: Vec<Vista> = amostra.into_iter().map(|c| c.vista).collect();
    vistas.sort_by(|a, b| a.inicio.cmp(&b.inicio).then(a.fim.cmp(&b.fim)).then_with(|| a.nome.cmp(&b.nome)).then(a.flag.cmp(&b.flag)));
    let linhas = empacota(&mut vistas, 1);
    let mostradas = vistas.len();
    let amostrado = (mostradas as u64) < total;
    if amostrado {
        avisos.push(format!(
            "{} reads pass the filters; a fixed sample of {} is drawn (max_reads). The coverage counts all {}.",
            milhar(total),
            milhar(mostradas as u64),
            milhar(total)
        ));
    }

    // De onde vem a referencia e as diferencas.
    let mut ref_janela: Option<Vec<u8>> = ref_fasta.clone();
    let origem_ref: Option<&str>;
    let fonte_difs: &str;
    if ref_fasta.is_some() {
        origem_ref = Some("FASTA");
        fonte_difs = "reference";
        let discordantes = vistas.iter().filter(|v| v.discordancias > 0).count();
        let posicoes: usize = vistas.iter().map(|v| v.discordancias).sum();
        if posicoes > 0 {
            avisos.push(format!(
                "At {posicoes} base(s) of {discordantes} read(s) drawn, the MD tag and the reference FASTA disagree: the reads were probably aligned to another version of the genome than {}, or the MD tags are stale (samtools calmd recomputes them).",
                opt.reference.as_deref().unwrap_or("")
            ));
        }
    } else {
        let com_md = vistas.iter().filter(|v| v.origem == OrigemDif::Md).count();
        if com_md > 0 {
            fonte_difs = "MD";
            origem_ref = Some("MD");
            let mut r = vec![b'N'; (fim - ini + 1) as usize];
            for v in &vistas {
                for &(p, b) in &v.ref_md {
                    let s = (p - ini) as usize;
                    if r[s] == b'N' {
                        r[s] = b;
                    }
                }
            }
            ref_janela = Some(r);
            let sem = vistas.iter().filter(|v| !v.tem_md).count();
            if sem > 0 {
                avisos.push(format!(
                    "{sem} of the {mostradas} reads drawn have no MD tag: their mismatches are not shown. Pass reference=\"genome.fa\" to compute them for every read."
                ));
            }
        } else if mostradas > 0 {
            fonte_difs = "none";
            origem_ref = None;
            avisos.push(
                "No reference FASTA was given and the reads have no usable MD tag: mismatches cannot be shown. Pass reference=\"genome.fa\" (the FASTA the reads were aligned to).".to_string(),
            );
        } else {
            fonte_difs = "none";
            origem_ref = None;
        }
    }
    let invalidos: Vec<&Vista> = vistas.iter().filter(|v| v.md_invalido).collect();
    if !invalidos.is_empty() {
        avisos.push(format!(
            "{} read(s) drawn have an MD tag that does not fit their CIGAR (e.g. {}): {}. Recompute the tags with samtools calmd.",
            invalidos.len(),
            invalidos[0].nome,
            if ref_fasta.is_some() { "the FASTA was used for them" } else { "their mismatches are not shown" }
        ));
    }
    let sem_seq = vistas.iter().filter(|v| !v.tem_seq).count();
    if sem_seq > 0 {
        avisos.push(format!("{sem_seq} read(s) drawn carry no sequence (SEQ is *): they are drawn without mismatches."));
    }

    // Cobertura.
    let mut cob_json: Vec<String> = Vec::with_capacity(cob.tot.len());
    let mut fr_json: Vec<String> = Vec::with_capacity(cob.tot.len());
    let (mut maxp, mut soma) = (0u32, 0u64);
    for (s, (t, f)) in cob.tot.iter().zip(&cob.frente).enumerate() {
        let p = ini + s as u64;
        let prof: u32 = t[..6].iter().sum();
        maxp = maxp.max(prof);
        soma += u64::from(prof);
        cob_json.push(format!("[{p},{prof},{},{},{},{},{},{}]", t[0], t[1], t[2], t[3], t[5], t[6]));
        fr_json.push(format!("[{p},{},{},{},{},{},{}]", f[0], f[1], f[2], f[3], f[5], f[6]));
    }
    let largura = fim - ini + 1;
    let media = soma as f64 / largura as f64;

    let opcoes = format!(
        "min mapping quality {}, min base quality {} (coverage), FLAG filter 0x{:X}, max reads {}",
        opt.min_mapq, opt.min_baseq, opt.exclude_flags, opt.max_reads
    );
    let mut texto = format!(
        "Reads of {} in {}:{}-{} ({} bases{})\n",
        fonte.caminho,
        nome_final,
        milhar(ini),
        milhar(fim),
        milhar(largura),
        j.foco.map(|p| format!(", centred on {}", milhar(p))).unwrap_or_default()
    );
    texto.push_str(&format!(
        "  {} read(s) pass the filters ({} left out by FLAG, {} by mapping quality); {} drawn in {} row(s)\n",
        milhar(total),
        milhar(f_flag),
        milhar(f_mapq),
        milhar(mostradas as u64),
        linhas
    ));
    texto.push_str(&format!("  depth: mean {}, max {} ({})\n", jdec(media, 1), maxp, opcoes));
    texto.push_str(&format!(
        "  mismatches: {}\n",
        match fonte_difs {
            "reference" => format!("against {}", opt.reference.as_deref().unwrap_or("")),
            "MD" => "from the MD tags".to_string(),
            _ => "not shown".to_string(),
        }
    ));
    for a in &avisos {
        texto.push_str(&format!("  warning: {a}\n"));
    }

    let leituras: Vec<String> = vistas.iter().map(Vista::json).collect();
    Ok(format!(
        "{{\"path\":{},\"region\":{{\"chrom\":{},\"start\":{},\"end\":{},\"focus\":{},\"length\":{},\"sequence_length\":{},\"requested\":{}}},\"indexed\":{},\"reference\":{},\"reference_source\":{},\"mismatch_source\":\"{}\",\"coverage\":[{}],\"coverage_forward\":[{}],\"mean_depth\":{},\"max_depth\":{},\"reads\":[{}],\"rows\":{},\"reads_total\":{},\"reads_shown\":{},\"downsampled\":{},\"filtered_flag\":{},\"filtered_mapq\":{},\"options\":{},\"warnings\":{},\"text\":{}}}",
        json_texto(&fonte.caminho),
        json_texto(&nome_final),
        ini,
        fim,
        j.foco.map(|p| p.to_string()).unwrap_or_else(|| "null".into()),
        largura,
        j.comprimento.map(|p| p.to_string()).unwrap_or_else(|| "null".into()),
        json_texto(region.trim()),
        fonte.indexado(),
        ref_janela.as_ref().map(|r| json_texto(&String::from_utf8_lossy(r))).unwrap_or_else(|| "null".into()),
        jopt(origem_ref),
        fonte_difs,
        cob_json.join(","),
        fr_json.join(","),
        jdec(media, 3),
        maxp,
        leituras.join(","),
        linhas,
        total,
        mostradas,
        amostrado,
        f_flag,
        f_mapq,
        json_texto(&opcoes),
        jlista(&avisos),
        json_texto(&texto)
    ))
}

// ---------------------------------------------------------------------
//  PERFIL DE PROFUNDIDADE
// ---------------------------------------------------------------------

struct RegiaoPerfil {
    texto: String,
    nome: String,
    exato: bool,
    ini: u64,
    fim: u64,
    bin: u64,
    /// Bases alinhadas por bin.
    somas: Vec<u64>,
    leituras: u64,
}

impl RegiaoPerfil {
    /// Soma as bases alinhadas (M, =, X) da leitura que caem na regiao.
    fn soma(&mut self, l: &Leitura) {
        let mut r = l.pos;
        let mut tocou = false;
        for &(len, op) in &l.cigar {
            let n = u64::from(len);
            match op {
                b'M' | b'=' | b'X' => {
                    if n > 0 {
                        let (mut x, y) = (r.max(self.ini), (r + n - 1).min(self.fim));
                        while x <= y {
                            let bi = (x - self.ini) / self.bin;
                            let fim_bin = self.ini + (bi + 1) * self.bin - 1;
                            let z = y.min(fim_bin);
                            self.somas[bi as usize] += z - x + 1;
                            tocou = true;
                            x = z + 1;
                        }
                    }
                    r += n;
                }
                b'D' | b'N' => r += n,
                _ => {}
            }
            if r > self.fim {
                break;
            }
        }
        if tocou {
            self.leituras += 1;
        }
    }

    fn medias(&self) -> Vec<(u64, u64, f64)> {
        self.somas
            .iter()
            .enumerate()
            .map(|(i, &s)| {
                let a = self.ini + i as u64 * self.bin;
                let b = (a + self.bin - 1).min(self.fim);
                (a, b, s as f64 / (b - a + 1) as f64)
            })
            .collect()
    }
}

fn mediana(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

fn passa_profundidade(l: &Leitura, min_mapq: u8) -> bool {
    l.flag & DEPTH_EXCLUDE_FLAGS == 0 && l.mapq >= u32::from(min_mapq)
}

/// Mean depth per bin of `bin` bases over each region, counting aligned
/// bases only (M, =, X: no deletions, no N skips, no clipped bases), for
/// copy-number plots. Reads with FLAG 0x704 (unmapped, secondary, QC fail,
/// duplicate) or a mapping quality below `min_mapq` are left out; no
/// base-quality filter; overlapping mates both count (as `samtools depth`
/// without -s). Regions: "chr22:42126499-42130881", "chr22:42126499" (one
/// base) or "chr22" (the whole sequence). JSON.
pub fn depth_profile(alignments: &str, regions: &[String], bin: u32, min_mapq: u8) -> Result<String, String> {
    if regions.is_empty() {
        return Err("No region given. Use \"chr22:42126499-42130881\" (1-based, end included), or a whole sequence (\"chr22\")".into());
    }
    if bin == 0 {
        return Err("bin must be at least 1 base".into());
    }
    let fonte = Fonte::abre(alignments)?;
    let mut avisos: Vec<String> = Vec::new();
    fonte.avisos(&mut avisos);
    let bin = u64::from(bin);
    let mut regs: Vec<RegiaoPerfil> = Vec::with_capacity(regions.len());
    for texto in regions {
        let t = texto.trim();
        let (pedido, ini0, fim) = indice::le_regiao(t, &|n| fonte.conhece(n))?;
        let (nome, exato, comprimento) = match fonte.resolve(&pedido)? {
            Some((n, l)) => (n, true, if l > 0 { Some(l) } else { None }),
            None => (pedido.clone(), false, None),
        };
        let ini = ini0 as u64 + 1;
        let mut fim = fim as u64;
        match comprimento {
            Some(len) => {
                if ini > len {
                    return Err(format!("Region \"{t}\" starts after the end of {nome}, which has {} bases in {}.", milhar(len), fonte.caminho));
                }
                fim = fim.min(len);
            }
            None if fim >= indice::POS_MAX as u64 => {
                return Err(format!(
                    "Region \"{t}\" has no end, and {} does not declare the length of {nome}. Give it as \"{nome}:start-end\".",
                    fonte.caminho
                ));
            }
            None => {}
        }
        let n_bins = (fim - ini) / bin + 1;
        if n_bins > MAX_BINS {
            return Err(format!(
                "Region \"{t}\" would make {} bins of {bin} bases (at most {} per region). Use a larger bin.",
                milhar(n_bins),
                milhar(MAX_BINS)
            ));
        }
        regs.push(RegiaoPerfil { texto: t.to_string(), nome, exato, ini, fim, bin, somas: vec![0; n_bins as usize], leituras: 0 });
    }

    let caminho = fonte.caminho.clone();
    if fonte.indexado() {
        for reg in regs.iter_mut() {
            let (nome, exato, ini, fim) = (reg.nome.clone(), reg.exato, reg.ini, reg.fim);
            fonte.linhas(&nome, exato, ini, fim, &mut |linha| {
                if let Some(l) = Leitura::le(linha, &caminho)? {
                    if l.sobrepoe(ini, fim) && passa_profundidade(&l, min_mapq) {
                        reg.soma(&l);
                    }
                }
                Ok(())
            })?;
        }
    } else {
        // Uma passada so: cada leitura vai para as regioes que toca.
        // Por sequencia, as regioes em ordem de inicio, com o maior fim ate
        // cada uma (as regioes podem se sobrepor).
        let mut por_nome: Vec<(String, Vec<usize>, Vec<u64>)> = Vec::new();
        for (i, r) in regs.iter().enumerate() {
            let chave = if r.exato { r.nome.clone() } else { canonico(&r.nome) };
            match por_nome.iter_mut().find(|x| x.0 == chave) {
                Some(x) => x.1.push(i),
                None => por_nome.push((chave, vec![i], Vec::new())),
            }
        }
        for (_, idx, maxfim) in por_nome.iter_mut() {
            idx.sort_by_key(|&i| regs[i].ini);
            let mut m = 0u64;
            *maxfim = idx.iter().map(|&i| {
                m = m.max(regs[i].fim);
                m
            }).collect();
        }
        let exato = !fonte.seqs.is_empty();
        let mut r = abre_leitura(&fonte.caminho)?;
        let mut linha = String::new();
        let mut cache: (String, Option<usize>) = (String::new(), None);
        loop {
            linha.clear();
            let n = r.read_line(&mut linha).map_err(|e| erro_leitura(&fonte.caminho, &e))?;
            if n == 0 {
                break;
            }
            let t = linha.trim_end_matches(['\n', '\r']);
            if t.is_empty() || t.starts_with('@') {
                continue;
            }
            let rname = t.split('\t').nth(2).unwrap_or("");
            if rname != cache.0 {
                let chave = if exato { rname.to_string() } else { canonico(rname) };
                cache = (rname.to_string(), por_nome.iter().position(|x| x.0 == chave));
            }
            let Some(g) = cache.1 else { continue };
            let Some(l) = Leitura::le(t, &caminho)? else { continue };
            if !passa_profundidade(&l, min_mapq) {
                continue;
            }
            let (_, idx, maxfim) = &por_nome[g];
            // Regioes com inicio <= fim da leitura; de tras para a frente
            // enquanto o maior fim ainda alcanca o inicio da leitura.
            let k = idx.partition_point(|&i| regs[i].ini <= l.fim);
            for j in (0..k).rev() {
                if maxfim[j] < l.pos {
                    break;
                }
                let i = idx[j];
                if l.sobrepoe(regs[i].ini, regs[i].fim) {
                    regs[i].soma(&l);
                }
            }
        }
    }

    let mut partes: Vec<String> = Vec::with_capacity(regs.len());
    let mut texto = format!(
        "Depth profile of {} (bins of {} bases; min mapping quality {min_mapq}; FLAG filter 0x{:X}; aligned bases only)\n",
        fonte.caminho,
        milhar(bin),
        DEPTH_EXCLUDE_FLAGS
    );
    for reg in &regs {
        let medias = reg.medias();
        let total: u64 = reg.somas.iter().sum();
        let media = total as f64 / (reg.fim - reg.ini + 1) as f64;
        let mut ms: Vec<f64> = medias.iter().map(|m| m.2).collect();
        let med = mediana(&mut ms);
        let bins: Vec<String> = medias.iter().map(|(a, b, m)| format!("[{a},{b},{}]", jdec(*m, 4))).collect();
        texto.push_str(&format!(
            "  {} ({}:{}-{}): mean depth {}, median of the bins {}, {} read(s), {} bin(s)\n",
            reg.texto,
            reg.nome,
            milhar(reg.ini),
            milhar(reg.fim),
            jdec(media, 2),
            jdec(med, 2),
            milhar(reg.leituras),
            medias.len()
        ));
        partes.push(format!(
            "{{\"region\":{},\"chrom\":{},\"start\":{},\"end\":{},\"bins\":[{}],\"mean\":{},\"median\":{},\"reads\":{}}}",
            json_texto(&reg.texto),
            json_texto(&reg.nome),
            reg.ini,
            reg.fim,
            bins.join(","),
            jdec(media, 4),
            jdec(med, 4),
            reg.leituras
        ));
    }
    for a in &avisos {
        texto.push_str(&format!("  warning: {a}\n"));
    }
    Ok(format!(
        "{{\"path\":{},\"bin\":{},\"min_mapq\":{},\"exclude_flags\":{},\"indexed\":{},\"regions\":[{}],\"warnings\":{},\"text\":{}}}",
        json_texto(&fonte.caminho),
        bin,
        min_mapq,
        DEPTH_EXCLUDE_FLAGS,
        fonte.indexado(),
        partes.join(","),
        jlista(&avisos),
        json_texto(&texto)
    ))
}

// ---------------------------------------------------------------------
//  BASES DO GENOMA DE REFERENCIA
// ---------------------------------------------------------------------

/// The bases (uppercase) of `region` in a FASTA: "chr1:1000-1200" (1-based,
/// end included), "chr1:1000" (one base), "chr1:1000-" (to the end) or
/// "chr1" (the whole sequence; at most `MAX_SEQUENCE` bases). With or without
/// a `.fai` next to it, plain or compressed; `chr1` and `1` are matched.
pub fn fetch_sequence(fasta: &str, region: &str) -> Result<String, String> {
    if !Path::new(fasta).is_file() {
        return Err(format!("{fasta} does not exist (or is not a file). Give the path of a FASTA file."));
    }
    let mut g = Genoma::abre(fasta, None)?;
    let t = region.trim();
    let inteiro = !t.is_empty() && g.nome(t).is_some();
    let (pedido, ini0, fim) = indice::le_regiao(t, &|n| inteiro && n == t)?;
    let nome = g.nome(&pedido).ok_or_else(|| {
        format!("Sequence \"{pedido}\" is not in {fasta} (the name was also tried with and without \"chr\"). Check the name, or the genome build of the FASTA.")
    })?;
    let len = g.comprimento(&nome).unwrap_or(0);
    let ini = ini0 as u64;
    let aberto = fim >= indice::POS_MAX;
    let fim = fim as u64;
    if ini >= len {
        return Err(format!("Region \"{t}\" starts after the end of {nome}, which has {} bases in {fasta}.", milhar(len)));
    }
    if fim > len && !aberto {
        return Err(format!(
            "Region \"{t}\" ends after the end of {nome}, which has {} bases in {fasta}. Use \"{nome}:{}-{}\", or an open end (\"{nome}:{}-\").",
            milhar(len),
            ini + 1,
            len,
            ini + 1
        ));
    }
    let fim = fim.min(len);
    if fim - ini > MAX_SEQUENCE {
        return Err(format!(
            "Region \"{t}\" has {} bases; at most {} are returned at a time. Ask for a narrower region.",
            milhar(fim - ini),
            milhar(MAX_SEQUENCE)
        ));
    }
    let b = g.trecho(&nome, ini, fim)?;
    Ok(String::from_utf8_lossy(b).into_owned())
}
