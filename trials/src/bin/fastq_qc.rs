#![allow(dead_code)]
// Ensaios de src/fastq_qc.rs: o relatorio de qualidade e o aparo de FASTQ.
// Cada resultado e comparado com uma referencia independente, escrita aqui,
// sobre os registros em memoria: contagens por forca bruta (percentis por
// ordenacao), a transcricao direta dos lacos de aparo por qualidade, e uma
// busca ingenua de adaptador que tenta todas as posicoes e conta as trocas.
//
// Com a variavel de ambiente CUTADAPT apontando para o executavel do cutadapt,
// roda tambem o ensaio DIFERENCIAL: os mesmos arquivos passam por `fastq_trim`
// e pelo `cutadapt --no-indels`, e as saidas sao comparadas leitura por leitura.
//
// Uso: cargo run --release --bin fastq_qc [-- arquivo.fastq ...]
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
#[path = "../../../src/fastq_qc.rs"]
mod fastq_qc;

use fastq_qc::*;
use std::collections::HashSet;
use std::io::{Read, Write};

fn lcg(s: &mut u64) -> u64 {
    *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    *s >> 11
}

fn rnd(s: &mut u64, n: usize) -> usize {
    (lcg(s) % n as u64) as usize
}

#[derive(Clone, Debug, PartialEq)]
struct Rec {
    nome: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

fn rec(nome: &str, seq: &str, qual: &str) -> Rec {
    Rec { nome: nome.to_string(), seq: seq.as_bytes().to_vec(), qual: qual.as_bytes().to_vec() }
}

/// Um registro com qualidade 40 (`I`) em todas as bases.
fn reci(nome: &str, seq: &str) -> Rec {
    rec(nome, seq, &"I".repeat(seq.len()))
}

fn texto(recs: &[Rec], eol: &str) -> Vec<u8> {
    let mut o = Vec::new();
    for r in recs {
        o.push(b'@');
        o.extend_from_slice(r.nome.as_bytes());
        o.extend_from_slice(eol.as_bytes());
        o.extend_from_slice(&r.seq);
        o.extend_from_slice(eol.as_bytes());
        o.push(b'+');
        o.extend_from_slice(eol.as_bytes());
        o.extend_from_slice(&r.qual);
        o.extend_from_slice(eol.as_bytes());
    }
    o
}

/// Le um arquivo inteiro (texto, .gz ou .lz4) pelo leitor do projeto.
fn le(p: &str) -> Vec<u8> {
    let mut r = externo::abre_leitura(p).unwrap();
    let mut v = Vec::new();
    r.read_to_end(&mut v).unwrap();
    v
}

/// Os registros de um FASTQ ja gravado (4 linhas cada, fim de linha LF).
fn registros(bytes: &[u8]) -> Vec<Rec> {
    let mut linhas: Vec<&[u8]> = bytes.split(|b| *b == b'\n').collect();
    assert_eq!(linhas.pop(), Some(&b""[..]), "the output must end with a line feed");
    assert_eq!(linhas.len() % 4, 0, "the output must have 4 lines per record");
    linhas
        .chunks(4)
        .map(|c| {
            assert_eq!(c[0].first(), Some(&b'@'));
            assert_eq!(c[2], b"+", "the third line must be `+` alone");
            Rec { nome: String::from_utf8_lossy(&c[0][1..]).to_string(), seq: c[1].to_vec(), qual: c[3].to_vec() }
        })
        .collect()
}

fn bases(s: &mut u64, n: usize) -> Vec<u8> {
    (0..n).map(|_| b"ACGT"[rnd(s, 4)]).collect()
}

fn perto(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.0)
}

// ---------------------------------------------------------------------
//  Gerador de leituras para o relatorio
// ---------------------------------------------------------------------

struct Plano {
    n: usize,
    max_len: usize,
    zero: bool,
    qmin: u8,
    qmax: u8,
    por_mil_n: usize,
    minusculas: bool,
    adaptadores: bool,
    duplicadas: bool,
    cauda_ruim: bool,
}

fn gera(s: &mut u64, pl: &Plano) -> Vec<Rec> {
    let mut recs: Vec<Rec> = Vec::with_capacity(pl.n);
    for i in 0..pl.n {
        let seq = if pl.duplicadas && i > 10 && rnd(s, 100) < 10 {
            // leitura duplicada: a sequencia exata de uma anterior
            recs[rnd(s, i)].seq.clone()
        } else {
            let len = match rnd(s, 20) {
                0 if pl.zero => 0,
                1..=5 => 1 + rnd(s, pl.max_len),
                _ => pl.max_len,
            };
            let mut q = bases(s, len);
            for b in q.iter_mut() {
                if rnd(s, 1000) < pl.por_mil_n {
                    *b = b'N';
                }
            }
            if pl.adaptadores {
                // Illumina em 20%, Nextera em 4%, Small RNA em 0,25% (abaixo do limiar)
                let plantado = if i % 5 == 0 {
                    Some(0)
                } else if i % 25 == 1 {
                    Some(1)
                } else if i % 400 == 2 {
                    Some(2)
                } else {
                    None
                };
                if let Some(a) = plantado {
                    let p = rnd(s, len + 1);
                    q.truncate(p);
                    q.extend_from_slice(KNOWN_ADAPTERS[a].1.as_bytes());
                    let extra = rnd(s, 10);
                    q.extend(bases(s, extra));
                }
            }
            if pl.minusculas {
                match rnd(s, 10) {
                    0 => q.make_ascii_lowercase(),
                    1 => {
                        for b in q.iter_mut() {
                            if rnd(s, 2) == 0 {
                                *b = b.to_ascii_lowercase();
                            }
                        }
                    }
                    _ => {}
                }
            }
            q
        };
        let len = seq.len();
        // Com `cauda_ruim`, o ultimo terco da leitura tem as 12 piores
        // qualidades e o resto as 21 melhores; sem ela, a faixa inteira.
        let qual = (0..len)
            .map(|j| {
                let (de, ate) = match (pl.cauda_ruim, j * 3 >= len * 2) {
                    (false, _) => (pl.qmin, pl.qmax),
                    (true, true) => (pl.qmin, pl.qmin + 11),
                    (true, false) => (pl.qmax - 20, pl.qmax),
                };
                de + rnd(s, (ate - de) as usize + 1) as u8
            })
            .collect();
        recs.push(Rec { nome: format!("r{i} len={len}"), seq, qual });
    }
    recs
}

// ---------------------------------------------------------------------
//  Contagem independente, sobre os registros em memoria
// ---------------------------------------------------------------------

struct Conta {
    reads: u64,
    bases: u64,
    min_len: usize,
    max_len: usize,
    mean_len: f64,
    hist_len: Vec<(usize, u64)>,
    offset: u8,
    faixa_q: Option<(u8, u8)>,
    pos: Vec<PositionStats>,
    hist_q: Vec<u64>,
    mean_q: f64,
    gc_percent: f64,
    hist_gc: Vec<u64>,
    n_bases: u64,
    n_percent: f64,
    com_n: u64,
    dup_percent: f64,
    dup_sample: u64,
    adaptadores: Vec<(u64, f64)>,
    provavel: Option<(&'static str, &'static str)>,
}

fn conta(recs: &[Rec]) -> Conta {
    let reads = recs.len() as u64;
    let mut lens: Vec<usize> = recs.iter().map(|r| r.seq.len()).collect();
    lens.sort();
    let bases: u64 = lens.iter().map(|&l| l as u64).sum();
    let mut hist_len: Vec<(usize, u64)> = Vec::new();
    for &l in &lens {
        match hist_len.last_mut() {
            Some(u) if u.0 == l => u.1 += 1,
            _ => hist_len.push((l, 1)),
        }
    }
    // Phred+64 so com o menor caractere em 64 (`@`) ou mais E o maior acima
    // de 75 (`K`); em qualquer outro caso, Phred+33.
    let todos = || recs.iter().flat_map(|r| r.qual.iter().copied());
    let offset: u8 = match (todos().min(), todos().max()) {
        (Some(menor), Some(maior)) if menor >= 64 && maior > 75 => 64,
        _ => 33,
    };
    let q = |c: u8| (c - offset) as u64;
    let max_len = lens.last().copied().unwrap_or(0);

    let mut pos = Vec::new();
    for j in 0..max_len.min(1000) {
        let mut vals: Vec<u64> = recs.iter().filter(|r| r.seq.len() > j).map(|r| q(r.qual[j])).collect();
        vals.sort();
        let n = vals.len();
        // posto mais proximo: o valor de posto ceil(p*n/100) na lista ordenada
        let pc = |p: f64| -> u8 { vals[((p * n as f64 / 100.0).ceil() as usize).max(1) - 1] as u8 };
        let letra = |alvo: u8| recs.iter().filter(|r| r.seq.len() > j && r.seq[j].to_ascii_uppercase() == alvo).count() as u64;
        let (a, c, g, t) = (letra(b'A'), letra(b'C'), letra(b'G'), letra(b'T'));
        pos.push(PositionStats {
            bases: n as u64,
            mean_quality: vals.iter().sum::<u64>() as f64 / n as f64,
            q10: pc(10.0),
            q25: pc(25.0),
            q50: pc(50.0),
            q75: pc(75.0),
            q90: pc(90.0),
            a,
            c,
            g,
            t,
            n: n as u64 - a - c - g - t,
        });
    }

    let mut hist_q = vec![0u64; 94];
    let mut hist_gc = vec![0u64; 101];
    let (mut soma_q, mut gc_tot, mut acgt_tot, mut n_bases, mut com_n) = (0u64, 0u64, 0u64, 0u64, 0u64);
    for r in recs {
        let s: u64 = r.qual.iter().map(|&c| q(c)).sum();
        soma_q += s;
        if !r.seq.is_empty() {
            hist_q[(s as f64 / r.seq.len() as f64).floor() as usize] += 1;
        }
        let m = r.seq.to_ascii_uppercase();
        let gc = m.iter().filter(|b| matches!(b, b'G' | b'C')).count();
        let acgt = m.iter().filter(|b| matches!(b, b'A' | b'C' | b'G' | b'T')).count();
        if acgt > 0 {
            hist_gc[(100.0 * gc as f64 / acgt as f64).round() as usize] += 1;
        }
        gc_tot += gc as u64;
        acgt_tot += acgt as u64;
        n_bases += (m.len() - acgt) as u64;
        if m.len() > acgt {
            com_n += 1;
        }
    }

    // duplicacao: a sequencia exata, sem hash
    let amostra = recs.len().min(1_000_000);
    let mut vistos: HashSet<&[u8]> = HashSet::new();
    let repetidas = recs[..amostra].iter().filter(|r| !vistos.insert(&r.seq)).count();

    let pct = |a: u64, b: u64| if b == 0 { 0.0 } else { 100.0 * a as f64 / b as f64 };
    let adaptadores: Vec<(u64, f64)> = KNOWN_ADAPTERS
        .iter()
        .map(|(_, a)| {
            let c = recs.iter().filter(|r| String::from_utf8_lossy(&r.seq).to_ascii_uppercase().contains(&a[..12])).count() as u64;
            (c, pct(c, reads))
        })
        .collect();
    let mut provavel = None;
    let mut maior = -1.0;
    for (i, (_, p)) in adaptadores.iter().enumerate() {
        if *p > maior {
            maior = *p;
            provavel = Some(KNOWN_ADAPTERS[i]);
        }
    }
    if maior < 0.5 {
        provavel = None;
    }

    Conta {
        reads,
        bases,
        min_len: lens.first().copied().unwrap_or(0),
        max_len,
        mean_len: if reads == 0 { 0.0 } else { bases as f64 / reads as f64 },
        hist_len,
        offset,
        faixa_q: recs.iter().flat_map(|r| r.qual.iter().copied()).min().zip(recs.iter().flat_map(|r| r.qual.iter().copied()).max()),
        pos,
        hist_q,
        mean_q: if bases == 0 { 0.0 } else { soma_q as f64 / bases as f64 },
        gc_percent: pct(gc_tot, acgt_tot),
        hist_gc,
        n_bases,
        n_percent: pct(n_bases, bases),
        com_n,
        dup_percent: pct(repetidas as u64, amostra as u64),
        dup_sample: amostra as u64,
        adaptadores,
        provavel,
    }
}

fn confere(st: &FastqStats, c: &Conta, rot: &str) {
    assert_eq!((st.reads, st.bases), (c.reads, c.bases), "{rot}: reads, bases");
    assert_eq!((st.min_length, st.max_length), (c.min_len, c.max_len), "{rot}: lengths");
    assert!(perto(st.mean_length, c.mean_len), "{rot}: mean_length {} x {}", st.mean_length, c.mean_len);
    assert_eq!(st.length_histogram, c.hist_len, "{rot}: length_histogram");
    assert_eq!(st.phred_offset, c.offset, "{rot}: phred_offset");
    assert_eq!(st.quality_char_range, c.faixa_q, "{rot}: quality_char_range");
    assert_eq!(st.positions_reported, c.pos.len(), "{rot}: positions_reported");
    assert_eq!(st.positions.len(), c.pos.len(), "{rot}: positions");
    for (j, (a, b)) in st.positions.iter().zip(&c.pos).enumerate() {
        assert!(perto(a.mean_quality, b.mean_quality), "{rot}: cycle {}: mean quality {} x {}", j + 1, a.mean_quality, b.mean_quality);
        let sem_media = |p: &PositionStats| (p.bases, p.q10, p.q25, p.q50, p.q75, p.q90, p.a, p.c, p.g, p.t, p.n);
        assert_eq!(sem_media(a), sem_media(b), "{rot}: cycle {}", j + 1);
    }
    assert_eq!(st.read_quality_histogram, c.hist_q, "{rot}: read_quality_histogram");
    assert!(perto(st.mean_quality, c.mean_q), "{rot}: mean_quality {} x {}", st.mean_quality, c.mean_q);
    assert!(perto(st.gc_percent, c.gc_percent), "{rot}: gc_percent {} x {}", st.gc_percent, c.gc_percent);
    assert_eq!(st.gc_histogram, c.hist_gc, "{rot}: gc_histogram");
    assert_eq!((st.n_bases, st.reads_with_n), (c.n_bases, c.com_n), "{rot}: n_bases, reads_with_n");
    assert!(perto(st.n_percent, c.n_percent), "{rot}: n_percent");
    assert_eq!(st.duplication_sample, c.dup_sample, "{rot}: duplication_sample");
    assert!(perto(st.duplication_percent, c.dup_percent), "{rot}: duplication {} x {}", st.duplication_percent, c.dup_percent);
    assert_eq!(st.adapters.len(), 4, "{rot}: adapters");
    for (i, (a, b)) in st.adapters.iter().zip(&c.adaptadores).enumerate() {
        assert_eq!((a.name, a.sequence), KNOWN_ADAPTERS[i], "{rot}: adapter {i}");
        assert_eq!(a.reads, b.0, "{rot}: reads with adapter {}", a.name);
        assert!(perto(a.percent, b.1), "{rot}: percent of adapter {}", a.name);
    }
    assert_eq!(st.likely_adapter, c.provavel, "{rot}: likely_adapter");
}

// ---------------------------------------------------------------------
//  Um analisador de JSON minimo, para conferir que json() e JSON
// ---------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq)]
enum J {
    Nulo,
    Bool(bool),
    Num(f64),
    Txt(String),
    Lista(Vec<J>),
    Obj(Vec<(String, J)>),
}

struct Analisador<'a> {
    b: &'a [u8],
    i: usize,
}

impl Analisador<'_> {
    fn brancos(&mut self) {
        while self.i < self.b.len() && self.b[self.i].is_ascii_whitespace() {
            self.i += 1;
        }
    }
    fn come(&mut self, c: u8) {
        assert_eq!(self.b.get(self.i), Some(&c), "JSON: expected `{}` at byte {}", c as char, self.i);
        self.i += 1;
    }
    fn palavra(&mut self, p: &str, v: J) -> J {
        assert!(self.b[self.i..].starts_with(p.as_bytes()), "JSON: bad literal at byte {}", self.i);
        self.i += p.len();
        v
    }
    fn texto(&mut self) -> String {
        self.come(b'"');
        let mut o = Vec::new();
        loop {
            let c = *self.b.get(self.i).expect("JSON: unterminated string");
            self.i += 1;
            match c {
                b'"' => break,
                b'\\' => {
                    let e = self.b[self.i];
                    self.i += 1;
                    match e {
                        b'n' => o.push(b'\n'),
                        b'r' => o.push(b'\r'),
                        b't' => o.push(b'\t'),
                        b'u' => {
                            let h = std::str::from_utf8(&self.b[self.i..self.i + 4]).unwrap();
                            let c = u32::from_str_radix(h, 16).expect("JSON: bad \\u escape");
                            o.extend_from_slice(char::from_u32(c).unwrap_or('?').to_string().as_bytes());
                            self.i += 4;
                        }
                        b'"' | b'\\' | b'/' => o.push(e),
                        _ => panic!("JSON: bad escape `\\{}`", e as char),
                    }
                }
                c if c < 0x20 => panic!("JSON: raw control character in a string"),
                c => o.push(c),
            }
        }
        String::from_utf8(o).expect("JSON: string is not UTF-8")
    }
    fn valor(&mut self) -> J {
        self.brancos();
        match self.b[self.i] {
            b'{' => {
                self.i += 1;
                let mut campos = Vec::new();
                self.brancos();
                if self.b[self.i] == b'}' {
                    self.i += 1;
                    return J::Obj(campos);
                }
                loop {
                    self.brancos();
                    let k = self.texto();
                    self.brancos();
                    self.come(b':');
                    campos.push((k, self.valor()));
                    self.brancos();
                    if self.b[self.i] == b',' {
                        self.i += 1;
                    } else {
                        self.come(b'}');
                        return J::Obj(campos);
                    }
                }
            }
            b'[' => {
                self.i += 1;
                let mut itens = Vec::new();
                self.brancos();
                if self.b[self.i] == b']' {
                    self.i += 1;
                    return J::Lista(itens);
                }
                loop {
                    itens.push(self.valor());
                    self.brancos();
                    if self.b[self.i] == b',' {
                        self.i += 1;
                    } else {
                        self.come(b']');
                        return J::Lista(itens);
                    }
                }
            }
            b'"' => J::Txt(self.texto()),
            b't' => self.palavra("true", J::Bool(true)),
            b'f' => self.palavra("false", J::Bool(false)),
            b'n' => self.palavra("null", J::Nulo),
            _ => {
                let ini = self.i;
                while self.i < self.b.len() && matches!(self.b[self.i], b'0'..=b'9' | b'-' | b'+' | b'.' | b'e' | b'E') {
                    self.i += 1;
                }
                let t = std::str::from_utf8(&self.b[ini..self.i]).unwrap();
                J::Num(t.parse().unwrap_or_else(|_| panic!("JSON: `{t}` is not a number (byte {ini})")))
            }
        }
    }
}

fn json(s: &str) -> J {
    let mut a = Analisador { b: s.as_bytes(), i: 0 };
    let v = a.valor();
    a.brancos();
    assert_eq!(a.i, s.len(), "JSON: text after the value");
    v
}

impl J {
    fn campo(&self, k: &str) -> &J {
        match self {
            J::Obj(c) => &c.iter().find(|(n, _)| n == k).unwrap_or_else(|| panic!("JSON: missing key `{k}`")).1,
            _ => panic!("JSON: not an object"),
        }
    }
    fn chaves(&self) -> Vec<&str> {
        match self {
            J::Obj(c) => c.iter().map(|(n, _)| n.as_str()).collect(),
            _ => panic!("JSON: not an object"),
        }
    }
    fn num(&self) -> f64 {
        match self {
            J::Num(x) => *x,
            _ => panic!("JSON: not a number: {self:?}"),
        }
    }
    fn lista(&self) -> &Vec<J> {
        match self {
            J::Lista(l) => l,
            _ => panic!("JSON: not a list"),
        }
    }
}

// ---------------------------------------------------------------------
//  Referencias do aparo
// ---------------------------------------------------------------------

/// Transcricao direta dos dois lacos do enunciado (`q` ja em valores Phred).
fn ref_qualidade(q: &[i64], corte5: i64, corte3: i64) -> (usize, usize) {
    let len = q.len();
    let mut stop = len;
    if corte3 != 0 {
        let mut s = 0;
        let mut max_qual = 0;
        for i in (0..len).rev() {
            s += corte3 - q[i];
            if s < 0 {
                break;
            }
            if s > max_qual {
                max_qual = s;
                stop = i;
            }
        }
    }
    let mut start = 0;
    if corte5 != 0 {
        let mut s = 0;
        let mut max_qual = 0;
        for i in 0..len {
            s += corte5 - q[i];
            if s < 0 {
                break;
            }
            if s > max_qual {
                max_qual = s;
                start = i + 1;
            }
        }
    }
    if start >= stop {
        (0, 0)
    } else {
        (start, stop)
    }
}

/// A ponta 3' de outro jeito, como o BWA a descreve: as somas parciais de
/// (corte - q) do fim para o comeco, recalculadas a cada posicao; o corte e
/// onde a soma e maxima, sem passar de uma soma negativa.
fn ref_qualidade_3_bruta(q: &[i64], corte: i64) -> usize {
    let len = q.len();
    let soma = |i: usize| -> i64 { (i..len).map(|j| corte - q[j]).sum() };
    let mut melhor = (0i64, len);
    for i in (0..len).rev() {
        if soma(i) < 0 {
            break;
        }
        if soma(i) > melhor.0 {
            melhor = (soma(i), i);
        }
    }
    melhor.1
}

/// Busca ingenua do adaptador, com a regra do cutadapt escrita de outro
/// jeito: tenta TODAS as posicoes, conta TODAS as trocas, guarda os casamentos
/// aceitaveis e so depois escolhe entre eles.
///
/// Aceitavel: sobreposicao >= min_overlap (limitada ao adaptador) e trocas <=
/// taxa * (sobreposicao - N do adaptador nela). Pontuacao: acertos - trocas.
fn ref_adaptador(ad: &str, leitura: &[u8], min_overlap: usize, taxa: f64) -> Option<usize> {
    ref_adaptador_regra(ad, leitura, min_overlap, taxa, Regra::Cutadapt)
}

/// A regra do cutadapt e duas regras mais simples, que NAO sao a dele. As
/// duas ultimas so servem para medir, no ensaio diferencial, em quantas
/// leituras a diferenca apareceria: e o que mostra que o ensaio a pegaria.
#[derive(Clone, Copy, PartialEq)]
enum Regra {
    Cutadapt,
    /// No passo 2, compara o comeco do proprio casamento parcial (e nao o do
    /// prefixo que o cutadapt ainda acompanha) com a ocorrencia inteira.
    ComecoDoParcial,
    /// A posicao aceitavel mais a esquerda, sem olhar a pontuacao.
    MaisAEsquerda,
}

fn ref_adaptador_regra(ad: &str, leitura: &[u8], min_overlap: usize, taxa: f64, regra: Regra) -> Option<usize> {
    let a = ad.to_ascii_uppercase().into_bytes();
    let (m, n) = (a.len(), leitura.len());
    let curinga = a.contains(&b'N');
    let minimo = min_overlap.min(m);
    // trocas entre as `len` primeiras bases do adaptador e a leitura a partir de `p`
    let trocas = |p: usize, len: usize| -> usize {
        (0..len)
            .filter(|&i| {
                let (x, y) = (a[i], leitura[p + i].to_ascii_uppercase());
                !(x == b'N' || x == y || (curinga && x == b'T' && y == b'U'))
            })
            .count()
    };
    // (posicao, trocas, pontuacao) de cada casamento aceitavel, da esquerda para a direita
    let mut aceitos: Vec<(usize, usize, i64)> = Vec::new();
    for p in 0..n {
        let len = m.min(n - p);
        let t = trocas(p, len);
        let sem_n = a[..len].iter().filter(|c| **c != b'N').count();
        if len >= minimo && t as f64 <= sem_n as f64 * taxa {
            aceitos.push((p, t, len as i64 - 2 * t as i64));
        }
    }
    if regra == Regra::MaisAEsquerda {
        return aceitos.first().map(|c| c.0);
    }
    // (posicao, pontuacao, e uma ocorrencia inteira?)
    let mut melhor: Option<(usize, i64, bool)> = None;
    // 1. Ocorrencias inteiras: fica a primeira; uma posterior a substitui se
    //    comeca a ate meio adaptador dela e pontua mais; uma exata encerra.
    for &(p, t, pontos) in aceitos.iter().filter(|c| c.0 + m <= n) {
        let troca = match melhor {
            None => true,
            Some((onde, quanto, _)) => p <= onde + m / 2 && pontos > quanto,
        };
        if troca {
            if t == 0 {
                return Some(p);
            }
            melhor = Some((p, pontos, true));
        }
    }
    // 2. Casamentos que chegam ao fim da leitura, do mais longo para o mais
    //    curto. Contra uma ocorrencia inteira, o cutadapt nao compara o comeco
    //    do casamento parcial, e sim o comeco do prefixo mais longo do
    //    adaptador que ele ainda acompanha na ultima coluna: o mais longo que
    //    termina na penultima base com ate floor(taxa * m) trocas, mais uma base.
    let teto = (taxa * m as f64) as usize;
    let vivo = if n < 2 { teto } else { (0..=m.min(n - 1)).rev().find(|&len| trocas(n - 1 - len, len) <= teto).unwrap_or(0) };
    let acompanhada = n.saturating_sub(m.min(vivo + 1));
    for &(p, _, pontos) in aceitos.iter().filter(|c| c.0 + m >= n) {
        let troca = match melhor {
            None => true,
            Some((_, quanto, false)) => pontos > quanto,
            Some((onde, quanto, true)) => pontos > quanto && (if regra == Regra::ComecoDoParcial { p } else { acompanhada }) <= onde + m / 2,
        };
        if troca {
            melhor = Some((p, pontos, false));
        }
    }
    melhor.map(|b| b.0)
}

struct Aparada {
    seq: Vec<u8>,
    qual: Vec<u8>,
    rq: u64,
    ra: u64,
    rn: u64,
    achou: bool,
}

fn e_n(b: &u8) -> bool {
    b.eq_ignore_ascii_case(&b'N')
}

fn ref_apara(r: &Rec, o: &TrimOptions, ad: Option<&str>) -> Aparada {
    let q: Vec<i64> = r.qual.iter().map(|&c| c as i64 - o.phred_offset as i64).collect();
    let (ini, fim) = ref_qualidade(&q, o.quality_cutoff_5 as i64, o.quality_cutoff as i64);
    let mut seq = r.seq[ini..fim].to_vec();
    let mut qual = r.qual[ini..fim].to_vec();
    let rq = (r.seq.len() - seq.len()) as u64;
    let (mut ra, mut achou) = (0u64, false);
    if let Some(p) = ad.and_then(|a| ref_adaptador(a, &seq, o.min_overlap, o.max_error_rate)) {
        ra = (seq.len() - p) as u64;
        seq.truncate(p);
        qual.truncate(p);
        achou = true;
    }
    let mut rn = 0u64;
    if o.trim_n {
        let antes = seq.len();
        // so o N maiusculo, como o --trim-n do cutadapt
        let na_frente = seq.iter().take_while(|b| **b == b'N').count();
        seq.drain(..na_frente);
        qual.drain(..na_frente);
        while seq.last() == Some(&b'N') {
            seq.pop();
            qual.pop();
        }
        rn = (antes - seq.len()) as u64;
    }
    Aparada { seq, qual, rq, ra, rn, achou }
}

fn ref_motivos(a: &Aparada, o: &TrimOptions) -> [bool; 3] {
    let curta = a.seq.len() < o.min_length;
    let ns = a.seq.iter().filter(|b| e_n(b)).count();
    let muitos = o.max_n.map(|m| ns > m).unwrap_or(false);
    let baixa = match o.min_mean_quality {
        None => false,
        Some(q) => {
            a.seq.is_empty()
                || a.qual.iter().map(|&c| c as f64 - o.phred_offset as f64).sum::<f64>() / (a.qual.len() as f64) < q as f64
        }
    };
    [curta, muitos, baixa]
}

fn ref_conta(s: &mut TrimSummary, r: &Rec, a: &Aparada) {
    s.bases_in += r.seq.len() as u64;
    s.bases_removed_quality += a.rq;
    s.bases_removed_adapter += a.ra;
    s.bases_removed_n += a.rn;
}

fn ref_descarta(s: &mut TrimSummary, m: [bool; 3]) -> bool {
    if m[0] {
        s.discarded_too_short += 1;
    } else if m[1] {
        s.discarded_too_many_n += 1;
    } else if m[2] {
        s.discarded_low_quality += 1;
    }
    m[0] || m[1] || m[2]
}

fn ref_simples(recs: &[Rec], o: &TrimOptions) -> (Vec<u8>, TrimSummary) {
    let mut s = TrimSummary::default();
    let mut fica = Vec::new();
    for r in recs {
        let a = ref_apara(r, o, o.adapter.as_deref());
        s.reads_in += 1;
        ref_conta(&mut s, r, &a);
        s.reads_with_adapter += a.achou as u64;
        if ref_descarta(&mut s, ref_motivos(&a, o)) {
            continue;
        }
        s.reads_out += 1;
        s.bases_out += a.seq.len() as u64;
        fica.push(Rec { nome: r.nome.clone(), seq: a.seq, qual: a.qual });
    }
    (texto(&fica, "\n"), s)
}

fn ref_pareado(r1: &[Rec], r2: &[Rec], o: &TrimOptions) -> (Vec<u8>, Vec<u8>, TrimSummary) {
    let mut s = TrimSummary { paired: true, ..TrimSummary::default() };
    let (mut f1, mut f2) = (Vec::new(), Vec::new());
    for (x, y) in r1.iter().zip(r2) {
        let a = ref_apara(x, o, o.adapter.as_deref());
        let b = ref_apara(y, o, o.adapter2.as_deref().or(o.adapter.as_deref()));
        s.reads_in += 1;
        ref_conta(&mut s, x, &a);
        ref_conta(&mut s, y, &b);
        s.reads_with_adapter += a.achou as u64;
        s.reads_with_adapter2 += b.achou as u64;
        let (m, n) = (ref_motivos(&a, o), ref_motivos(&b, o));
        if ref_descarta(&mut s, [m[0] || n[0], m[1] || n[1], m[2] || n[2]]) {
            continue;
        }
        s.reads_out += 1;
        s.bases_out += (a.seq.len() + b.seq.len()) as u64;
        f1.push(Rec { nome: x.nome.clone(), seq: a.seq, qual: a.qual });
        f2.push(Rec { nome: y.nome.clone(), seq: b.seq, qual: b.qual });
    }
    (texto(&f1, "\n"), texto(&f2, "\n"), s)
}

/// Troca `n` bases de `seq` por outras, em posicoes distintas.
fn planta_trocas(s: &mut u64, seq: &mut [u8], n: usize) {
    let mut feitas: Vec<usize> = Vec::new();
    while feitas.len() < n.min(seq.len()) {
        let p = rnd(s, seq.len());
        if !feitas.contains(&p) {
            let outras: Vec<u8> = b"ACGT".iter().copied().filter(|b| *b != seq[p].to_ascii_uppercase()).collect();
            seq[p] = outras[rnd(s, outras.len())];
            feitas.push(p);
        }
    }
}

/// Uma leitura para os ensaios de adaptador. `tipo`: 0 sem adaptador; 1
/// adaptador inteiro no meio, com 0 a 2 trocas; 2 comeco do adaptador na
/// ponta 3', com 0 a 2 trocas; 3 e 4 como 1 e 2, com N espalhado; 5 leitura
/// curta; 6 como 1, em minusculas.
fn leitura_com_adaptador(s: &mut u64, ad: &[u8], tipo: usize) -> Vec<u8> {
    let trocas = rnd(s, 3);
    let mut r = match tipo {
        0 => {
            let n = rnd(s, 80);
            bases(s, n)
        }
        1 | 3 | 6 => {
            let n = rnd(s, 50);
            let mut r = bases(s, n);
            let mut a = ad.to_vec();
            planta_trocas(s, &mut a, trocas);
            r.extend(a);
            let n = rnd(s, 20);
            r.extend(bases(s, n));
            r
        }
        2 | 4 => {
            let n = rnd(s, 50);
            let mut r = bases(s, n);
            let mut a = ad[..1 + rnd(s, ad.len())].to_vec();
            let quantas = trocas.min(a.len() / 3);
            planta_trocas(s, &mut a, quantas);
            r.extend(a);
            r
        }
        _ => {
            let n = rnd(s, 6);
            let mut r = bases(s, n);
            let k = rnd(s, ad.len().min(4) + 1);
            r.extend_from_slice(&ad[..k]);
            r
        }
    };
    // O N do adaptador chega a leitura como uma base qualquer.
    for b in r.iter_mut() {
        if *b == b'N' {
            *b = b"ACGT"[rnd(s, 4)];
        }
    }
    if tipo == 3 || tipo == 4 {
        for b in r.iter_mut() {
            if rnd(s, 12) == 0 {
                *b = b'N';
            }
        }
    }
    if tipo == 6 {
        r.make_ascii_lowercase();
    }
    r
}

const ILLUMINA: &str = "AGATCGGAAGAGC";

/// Leituras para o aparo: adaptador em parte delas, N nas pontas e no meio,
/// qualidade que cai no fim, algumas em minusculas, algumas vazias.
fn gera_aparo(s: &mut u64, n: usize, ads: &[&str], base: u8, rotulo: &str) -> Vec<Rec> {
    (0..n)
        .map(|i| {
            let mut seq = match rnd(s, 10) {
                0..=3 => {
                    let (qual, tipo) = (rnd(s, ads.len()), 1 + rnd(s, 4));
                    leitura_com_adaptador(s, ads[qual].as_bytes(), tipo)
                }
                4 => {
                    let k = rnd(s, 4);
                    bases(s, k)
                }
                _ => {
                    let k = rnd(s, 100);
                    bases(s, k)
                }
            };
            match rnd(s, 8) {
                0 => {
                    let k = 1 + rnd(s, 4);
                    seq.splice(0..0, std::iter::repeat(b'N').take(k));
                }
                1 => {
                    let k = 1 + rnd(s, 4);
                    seq.extend(std::iter::repeat(b'N').take(k));
                }
                2 => {
                    for b in seq.iter_mut() {
                        if rnd(s, 10) == 0 {
                            *b = b'N';
                        }
                    }
                }
                3 => seq.make_ascii_lowercase(),
                _ => {}
            }
            let len = seq.len();
            let (frente, cauda) = (rnd(s, 3) == 0, rnd(s, 2) == 0);
            let qual = (0..len)
                .map(|j| {
                    let ruim = (frente && j < 4) || (cauda && j * 4 >= len * 3);
                    base + if ruim { 2 + rnd(s, 12) as u8 } else { 15 + rnd(s, 26) as u8 }
                })
                .collect();
            Rec { nome: format!("{rotulo}{i}"), seq, qual }
        })
        .collect()
}

fn opcoes_aleatorias(s: &mut u64, ads: &[&str], pareado: bool) -> TrimOptions {
    TrimOptions {
        adapter: if rnd(s, 4) == 0 { None } else { Some(ads[rnd(s, ads.len())].to_string()) },
        adapter2: if pareado && rnd(s, 2) == 0 { Some(ads[rnd(s, ads.len())].to_string()) } else { None },
        min_overlap: [1, 3, 5][rnd(s, 3)],
        max_error_rate: [0.0, 0.1, 0.2][rnd(s, 3)],
        quality_cutoff: [0, 10, 20, 30][rnd(s, 4)],
        quality_cutoff_5: [0, 0, 15, 25][rnd(s, 4)],
        trim_n: rnd(s, 2) == 0,
        min_length: [0, 0, 10, 30][rnd(s, 4)],
        max_n: [None, None, Some(0), Some(2)][rnd(s, 4)],
        min_mean_quality: [None, None, Some(15), Some(25)][rnd(s, 4)],
        phred_offset: 33,
    }
}

// ---------------------------------------------------------------------
//  Ensaio diferencial contra o cutadapt de verdade (variavel CUTADAPT)
// ---------------------------------------------------------------------

const NEXTERA: &str = "CTGTCTCTTATACACATCT";

/// Adaptadores do ensaio diferencial: os comuns, repetitivos (poli-A de 10 a
/// 30, poli-G, repeticoes de 2 e 3 bases), com N no meio, no comeco e no fim,
/// e curtos.
const ADS_DIF: [&str; 16] = [
    ILLUMINA,
    "AAAAAAAAAA",
    "ACGNNACGTTGCA",
    "AAAAAAAAAAAAAAAAAAAA",
    NEXTERA,
    "GGGGGGGGGG",
    "NNACGTACGTAC",
    "ACACACACACAC",
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    "TGC",
    "AGATCGGAAGAGCNNNNNNATCTCGTAT",
    "AACAACAACAAC",
    "ACGTACGTNN",
    "AAAAAAAAAAAAAAA",
    "GATTACA",
    "TTTTTTTTTTNNTTTTTTTT",
];

struct CasoDif {
    o: TrimOptions,
    pareado: bool,
}

/// As combinacoes de opcoes do ensaio diferencial. Os ciclos tem periodos
/// diferentes, de modo que cada adaptador passa pelas tres taxas de erro e
/// por varias sobreposicoes, cortes de qualidade e filtros.
fn casos_dif() -> Vec<CasoDif> {
    let mut v = Vec::new();
    for c in 0..56usize {
        v.push(CasoDif {
            pareado: false,
            o: TrimOptions {
                adapter: if c % 11 == 10 { None } else { Some(ADS_DIF[c % 16].to_string()) },
                adapter2: None,
                max_error_rate: [0.0, 0.1, 0.2][c % 3],
                min_overlap: [1, 3, 5, 10][(c / 3) % 4],
                quality_cutoff: [0, 20, 20, 10, 0, 30, 15][c % 7],
                quality_cutoff_5: [0, 0, 15, 0, 25, 0, 20][c % 7],
                trim_n: (c / 2) % 2 == 1,
                min_length: [0, 0, 1, 10, 30][c % 5],
                max_n: [None, Some(0), None, Some(2), Some(1)][(c / 2) % 5],
                min_mean_quality: None,
                phred_offset: if c % 9 == 8 { 64 } else { 33 },
            },
        });
    }
    for c in 0..28usize {
        let so_o_segundo = c % 7 == 6;
        v.push(CasoDif {
            pareado: true,
            o: TrimOptions {
                adapter: if so_o_segundo { None } else { Some(ADS_DIF[(c * 3) % 16].to_string()) },
                adapter2: if c % 3 == 0 && !so_o_segundo { None } else { Some(ADS_DIF[(c * 5 + 1) % 16].to_string()) },
                max_error_rate: [0.1, 0.2, 0.0][c % 3],
                min_overlap: [3, 1, 10, 5][c % 4],
                quality_cutoff: [0, 20, 15, 0, 25][c % 5],
                quality_cutoff_5: [0, 0, 10, 20, 0][c % 5],
                trim_n: c % 2 == 0,
                min_length: [0, 5, 0, 20][(c / 2) % 4],
                max_n: [None, Some(1), Some(0), None][(c / 3) % 4],
                min_mean_quality: None,
                phred_offset: if c % 10 == 9 { 64 } else { 33 },
            },
        });
    }
    v
}

/// Os argumentos do cutadapt equivalentes as nossas opcoes.
fn args_cutadapt(o: &TrimOptions, pareado: bool) -> Vec<String> {
    let mut a: Vec<String> = ["--no-indels", "-j", "1"].iter().map(|x| x.to_string()).collect();
    let mut par = |chave: &str, valor: String| a.extend([chave.to_string(), valor]);
    if let Some(ad) = &o.adapter {
        par("-a", ad.clone());
    }
    if pareado {
        // adapter2 ausente: a leitura 2 usa `adapter`, isto e, `-a X -A X`
        if let Some(ad) = o.adapter2.as_ref().or(o.adapter.as_ref()) {
            par("-A", ad.clone());
        }
    }
    par("-e", format!("{}", o.max_error_rate));
    par("-O", o.min_overlap.to_string());
    if o.quality_cutoff_5 > 0 {
        par("-q", format!("{},{}", o.quality_cutoff_5, o.quality_cutoff));
    } else if o.quality_cutoff > 0 {
        par("-q", o.quality_cutoff.to_string());
    }
    if o.min_length > 0 {
        par("-m", o.min_length.to_string());
    }
    if let Some(m) = o.max_n {
        par("--max-n", m.to_string());
    }
    if o.phred_offset == 64 {
        a.push("--quality-base=64".to_string());
    }
    if o.trim_n {
        a.push("--trim-n".to_string());
    }
    a
}

/// Uma leitura do ensaio diferencial, feita para provocar o adaptador `ad`:
/// adaptador inteiro ou so o comeco, com 0 a 2 trocas, uma ou duas copias,
/// caudas repetitivas com interrupcoes, leituras de baixa complexidade (que
/// dao muitos casamentos sobrepostos), curtas, com N, com minusculas.
fn leitura_dif(s: &mut u64, ad: &[u8]) -> Vec<u8> {
    let comum = *b"ACGT".iter().max_by_key(|b| ad.iter().filter(|x| x == b).count()).unwrap();
    let sorteada = b"ACGT"[rnd(s, 4)];
    let outra = ad.iter().copied().find(|b| *b != comum && *b != b'N').unwrap_or(sorteada);
    let copia = |s: &mut u64, ate: usize| -> Vec<u8> {
        let mut a = ad[..ate].to_vec();
        let k = rnd(s, 3).min(ate / 3);
        planta_trocas(s, &mut a, k);
        a
    };
    let trecho = |s: &mut u64, max: usize| -> Vec<u8> {
        let k = rnd(s, max + 1);
        bases(s, k)
    };
    let mut r = match rnd(s, 14) {
        0 => trecho(s, 90),
        1 | 2 => [trecho(s, 60), copia(s, ad.len()), trecho(s, 25)].concat(),
        3 | 4 => {
            let ate = 1 + rnd(s, ad.len());
            [trecho(s, 60), copia(s, ate)].concat()
        }
        5 => [copia(s, ad.len()), trecho(s, 30)].concat(),
        6 => {
            let ate = rnd(s, ad.len().min(5) + 1);
            [trecho(s, 6), ad[..ate].to_vec()].concat()
        }
        7 => [trecho(s, 30), copia(s, ad.len()), trecho(s, 6), copia(s, ad.len()), trecho(s, 10)].concat(),
        8 | 9 => {
            let mut r = trecho(s, 40);
            for _ in 0..3 + rnd(s, 43) {
                r.push(if rnd(s, 8) == 0 { b"ACGT"[rnd(s, 4)] } else { comum });
            }
            if rnd(s, 2) == 0 {
                r.extend(trecho(s, 5));
            }
            r
        }
        10 | 11 => (0..rnd(s, 70)).map(|_| if rnd(s, 3) == 0 { outra } else { comum }).collect(),
        _ => {
            // Cauda quase pura, pouco maior que o adaptador, com 1 a 3
            // interrupcoes longe do fim: ocorrencias inteiras com trocas
            // seguidas de um casamento parcial exato. E onde o cutadapt
            // decide entre a ocorrencia inteira e o parcial.
            let mut r = trecho(s, 30);
            let cauda = ad.len() + 2 + rnd(s, ad.len() + 4);
            let inicio = r.len();
            r.extend(std::iter::repeat(comum).take(cauda));
            for _ in 0..1 + rnd(s, 3) {
                let onde = inicio + rnd(s, cauda * 2 / 3 + 1);
                r[onde] = outra;
            }
            r
        }
    };
    // o N do adaptador chega a leitura como uma base qualquer
    for b in r.iter_mut() {
        if *b == b'N' {
            *b = b"ACGT"[rnd(s, 4)];
        }
    }
    match rnd(s, 12) {
        0 => {
            for b in r.iter_mut() {
                if rnd(s, 10) == 0 {
                    *b = b'N';
                }
            }
        }
        1 | 2 => {
            // N nas pontas; as vezes minusculo, que o --trim-n nao tira
            let letra = if rnd(s, 3) == 0 { b'n' } else { b'N' };
            let k = 1 + rnd(s, 4);
            if rnd(s, 2) == 0 {
                r.splice(0..0, std::iter::repeat(letra).take(k));
            } else {
                r.extend(std::iter::repeat(letra).take(k));
            }
            if rnd(s, 3) == 0 {
                r.extend(std::iter::repeat(b'N').take(k));
            }
        }
        3 => r.make_ascii_lowercase(),
        4 => {
            for b in r.iter_mut() {
                if rnd(s, 2) == 0 {
                    *b = b.to_ascii_lowercase();
                }
            }
        }
        5 => {
            // RNA e codigos de ambiguidade: U no lugar de T, e letras IUPAC
            for b in r.iter_mut() {
                if b.eq_ignore_ascii_case(&b'T') && rnd(s, 2) == 0 {
                    *b = if b.is_ascii_lowercase() { b'u' } else { b'U' };
                } else if rnd(s, 15) == 0 {
                    *b = b"RYSWKMBDHV"[rnd(s, 10)];
                }
            }
        }
        _ => {}
    }
    r
}

/// Qualidades de varios perfis: boas, cauda ruim, comeco ruim, tudo ruim,
/// ao acaso, e as quatro faixas dos instrumentos recentes.
fn qualidades_dif(s: &mut u64, len: usize, base: u8) -> Vec<u8> {
    let perfil = rnd(s, 6);
    (0..len)
        .map(|j| {
            base + match perfil {
                0 => 20 + rnd(s, 22),
                1 => if j * 3 >= len * 2 { 2 + rnd(s, 13) } else { 20 + rnd(s, 22) },
                2 => if j * 4 < len { 2 + rnd(s, 13) } else { 20 + rnd(s, 22) },
                3 => rnd(s, 16),
                4 => rnd(s, 42),
                _ => [2, 12, 23, 37][rnd(s, 4)],
            } as u8
        })
        .collect()
}

fn leituras_dif(s: &mut u64, n: usize, ad: &str, base: u8, nome: impl Fn(usize) -> String) -> Vec<Rec> {
    (0..n)
        .map(|i| {
            let seq = leitura_dif(s, ad.as_bytes());
            let qual = qualidades_dif(s, seq.len(), base);
            Rec { nome: nome(i), seq, qual }
        })
        .collect()
}

/// Como chamar o cutadapt. O lancador `cutadapt.exe` do Windows pode levar
/// dezenas de segundos para subir a cada chamada; se ao lado dele estiver o
/// Python do mesmo ambiente e `python -m cutadapt` responder, e por ali que se
/// chama -- e o mesmo programa.
fn comando_cutadapt(caminho: &str) -> (std::path::PathBuf, Vec<String>) {
    let exe = std::path::PathBuf::from(caminho);
    if let Some(pasta) = exe.parent() {
        for nome in ["python.exe", "python", "python3"] {
            let py = pasta.join(nome);
            let responde = py.is_file()
                && std::process::Command::new(&py).args(["-m", "cutadapt", "--version"]).output().map(|o| o.status.success()).unwrap_or(false);
            if responde {
                return (py, vec!["-m".to_string(), "cutadapt".to_string()]);
            }
        }
    }
    (exe, Vec::new())
}

fn roda_cutadapt(cmd: &(std::path::PathBuf, Vec<String>), pasta: &std::path::Path, args: &[String]) -> String {
    let saida = std::process::Command::new(&cmd.0)
        .args(&cmd.1)
        .args(args)
        .current_dir(pasta)
        .output()
        .unwrap_or_else(|e| panic!("could not run cutadapt ({}): {e}", cmd.0.display()));
    assert!(
        saida.status.success(),
        "cutadapt failed with the arguments {}:\n{}\n{}",
        args.join(" "),
        String::from_utf8_lossy(&saida.stdout),
        String::from_utf8_lossy(&saida.stderr)
    );
    String::from_utf8_lossy(&saida.stdout).to_string()
}

/// Compara, leitura por leitura e na ordem da entrada, o que nos gravamos com
/// o que o cutadapt gravou: nome ate o primeiro espaco, sequencia, qualidade.
fn compara_dif(rot: &str, entrada: &[Rec], nosso: &[Rec], deles: &[Rec], o: &TrimOptions, args: &[String]) {
    let chave = |r: &Rec| r.nome.split_whitespace().next().unwrap_or("").to_string();
    let mostra = |r: Option<&Rec>| match r {
        Some(r) => format!("seq `{}` qual `{}`", String::from_utf8_lossy(&r.seq), String::from_utf8_lossy(&r.qual)),
        None => "(read not written)".to_string(),
    };
    let (mut i, mut j) = (0, 0);
    for e in entrada {
        let k = chave(e);
        let a = nosso.get(i).filter(|r| chave(r) == k);
        let b = deles.get(j).filter(|r| chave(r) == k);
        if a.map(|r| (&r.seq, &r.qual)) != b.map(|r| (&r.seq, &r.qual)) {
            panic!(
                "DIFFERENCE from cutadapt in {rot}\n  options:  {o:?}\n  cutadapt: {}\n  read `{}`: {}\n  ours:     {}\n  cutadapt: {}",
                args.join(" "),
                e.nome,
                mostra(Some(e)),
                mostra(a),
                mostra(b)
            );
        }
        i += a.is_some() as usize;
        j += b.is_some() as usize;
    }
    assert!(i == nosso.len() && j == deles.len(), "{rot}: the outputs are not in the order of the input");
}

/// Compara os contadores do nosso resumo com o relatorio JSON do cutadapt.
fn contadores_dif(rot: &str, r: &TrimSummary, rel: &J, o: &TrimOptions, args: &[String]) {
    let num = |v: &J| match v {
        J::Nulo => 0u64,
        x => x.num() as u64,
    };
    // bases tiradas pelo adaptador: soma de (comprimento removido x leituras)
    let removidas = |chave: &str| -> u64 {
        match rel.campo(chave) {
            J::Nulo => 0,
            lista => lista
                .lista()
                .iter()
                .map(|a| match a.campo("three_prime_end") {
                    J::Nulo => 0,
                    t => t
                        .campo("trimmed_lengths")
                        .lista()
                        .iter()
                        .map(|x| num(x.campo("len")) * x.campo("counts").lista().iter().map(&num).sum::<u64>())
                        .sum::<u64>(),
                })
                .sum(),
        }
    };
    let (rc, bp) = (rel.campo("read_counts"), rel.campo("basepair_counts"));
    let deles = [
        ("reads_in", num(rc.campo("input"))),
        ("reads_out", num(rc.campo("output"))),
        ("discarded_too_short", num(rc.campo("filtered").campo("too_short"))),
        ("discarded_too_many_n", num(rc.campo("filtered").campo("too_many_n"))),
        ("reads_with_adapter", num(rc.campo("read1_with_adapter"))),
        ("reads_with_adapter2", num(rc.campo("read2_with_adapter"))),
        ("bases_in", num(bp.campo("input"))),
        ("bases_out", num(bp.campo("output"))),
        ("bases_removed_quality", num(bp.campo("quality_trimmed"))),
        ("bases_removed_adapter", removidas("adapters_read1") + removidas("adapters_read2")),
    ];
    let nossos = [
        r.reads_in,
        r.reads_out,
        r.discarded_too_short,
        r.discarded_too_many_n,
        r.reads_with_adapter,
        r.reads_with_adapter2,
        r.bases_in,
        r.bases_out,
        r.bases_removed_quality,
        r.bases_removed_adapter,
    ];
    for ((nome, d), n) in deles.iter().zip(nossos) {
        assert!(*d == n, "DIFFERENCE from cutadapt in {rot}: counter {nome} is {n} here and {d} in cutadapt's report\n  options:  {o:?}\n  cutadapt: {}", args.join(" "));
    }
}

/// O ensaio diferencial inteiro. Devolve quantos casos passaram.
fn diferencial(pasta: &std::path::Path, caminho: &str, s: &mut u64) -> u32 {
    let cmd = comando_cutadapt(caminho);
    let versao = roda_cutadapt(&cmd, pasta, &["--version".to_string()]);
    println!("differential test against cutadapt {} ({} {})", versao.trim(), cmd.0.display(), cmd.1.join(" "));
    let p = |n: &str| pasta.join(n).to_string_lossy().to_string();
    let casos = casos_dif();
    let (mut leituras, mut simples, mut pareados) = (0u64, 0u32, 0u32);
    let mut soma = TrimSummary::default();
    // Em quantas leituras uma regra mais simples daria outro corte: (a mais a
    // esquerda, o comeco do parcial). Mede a sensibilidade do ensaio.
    let mut outra_regra = (0u64, 0u64);
    let mut sensibilidade = |entrada: &[Rec], ad: Option<&str>, o: &TrimOptions| {
        let Some(ad) = ad else { return };
        for r in entrada {
            let q: Vec<i64> = r.qual.iter().map(|&c| c as i64 - o.phred_offset as i64).collect();
            let (ini, fim) = ref_qualidade(&q, o.quality_cutoff_5 as i64, o.quality_cutoff as i64);
            let certo = ref_adaptador(ad, &r.seq[ini..fim], o.min_overlap, o.max_error_rate);
            outra_regra.0 += (ref_adaptador_regra(ad, &r.seq[ini..fim], o.min_overlap, o.max_error_rate, Regra::MaisAEsquerda) != certo) as u64;
            outra_regra.1 += (ref_adaptador_regra(ad, &r.seq[ini..fim], o.min_overlap, o.max_error_rate, Regra::ComecoDoParcial) != certo) as u64;
        }
    };
    for (c, caso) in casos.iter().enumerate() {
        let o = &caso.o;
        let rot = format!("differential case {c} ({})", if caso.pareado { "paired" } else { "single-end" });
        let mut args = args_cutadapt(o, caso.pareado);
        let ad1 = o.adapter.as_deref().unwrap_or(ILLUMINA);
        let r = if caso.pareado {
            let ad2 = o.adapter2.as_deref().or(o.adapter.as_deref()).unwrap_or(ILLUMINA);
            let estilo = c % 2;
            let e1 = leituras_dif(s, 800, ad1, o.phred_offset, |i| if estilo == 0 { format!("d{i}/1") } else { format!("d{i} 1:N:0:ACGT") });
            let e2 = leituras_dif(s, 800, ad2, o.phred_offset, |i| if estilo == 0 { format!("d{i}/2") } else { format!("d{i} 2:N:0:ACGT") });
            std::fs::write(p("dif_e1.fastq"), texto(&e1, "\n")).unwrap();
            std::fs::write(p("dif_e2.fastq"), texto(&e2, "\n")).unwrap();
            let r = fastq_trim(&p("dif_e1.fastq"), &p("dif_n1.fastq"), Some(&p("dif_e2.fastq")), Some(&p("dif_n2.fastq")), o).unwrap();
            args.extend(["-o", "dif_c1.fastq", "-p", "dif_c2.fastq", "--json", "dif.json", "dif_e1.fastq", "dif_e2.fastq"].iter().map(|x| x.to_string()));
            roda_cutadapt(&cmd, pasta, &args);
            compara_dif(&format!("{rot}, read 1"), &e1, &registros(&le(&p("dif_n1.fastq"))), &registros(&le(&p("dif_c1.fastq"))), o, &args);
            compara_dif(&format!("{rot}, read 2"), &e2, &registros(&le(&p("dif_n2.fastq"))), &registros(&le(&p("dif_c2.fastq"))), o, &args);
            sensibilidade(&e1, o.adapter.as_deref(), o);
            sensibilidade(&e2, o.adapter2.as_deref().or(o.adapter.as_deref()), o);
            leituras += 1600;
            pareados += 1;
            r
        } else {
            let e1 = leituras_dif(s, 1500, ad1, o.phred_offset, |i| format!("d{i} x"));
            std::fs::write(p("dif_e1.fastq"), texto(&e1, "\n")).unwrap();
            let r = fastq_trim(&p("dif_e1.fastq"), &p("dif_n1.fastq"), None, None, o).unwrap();
            args.extend(["-o", "dif_c1.fastq", "--json", "dif.json", "dif_e1.fastq"].iter().map(|x| x.to_string()));
            roda_cutadapt(&cmd, pasta, &args);
            compara_dif(&rot, &e1, &registros(&le(&p("dif_n1.fastq"))), &registros(&le(&p("dif_c1.fastq"))), o, &args);
            sensibilidade(&e1, o.adapter.as_deref(), o);
            leituras += 1500;
            simples += 1;
            r
        };
        contadores_dif(&rot, &r, &json(&String::from_utf8_lossy(&std::fs::read(p("dif.json")).unwrap())), o, &args);
        soma.reads_with_adapter += r.reads_with_adapter + r.reads_with_adapter2;
        soma.bases_removed_adapter += r.bases_removed_adapter;
        soma.bases_removed_quality += r.bases_removed_quality;
        soma.bases_removed_n += r.bases_removed_n;
        soma.reads_out += r.reads_out;
        soma.discarded_too_short += r.discarded_too_short;
        soma.discarded_too_many_n += r.discarded_too_many_n;
    }
    println!(
        "cutadapt differential: {} option combinations ({simples} single-end, {pareados} paired), {leituras} reads compared, 0 differences",
        casos.len()
    );
    println!(
        "  in those runs: adapter found in {} reads ({} bases removed), {} bases removed by quality, {} N removed, {} reads or pairs written, {} dropped as too short, {} for too many N",
        soma.reads_with_adapter, soma.bases_removed_adapter, soma.bases_removed_quality, soma.bases_removed_n, soma.reads_out, soma.discarded_too_short, soma.discarded_too_many_n
    );
    println!(
        "  sensitivity: the leftmost acceptable position would give another cut in {} of those reads; comparing the start of a partial match instead of cutadapt's tracked prefix, in {}",
        outra_regra.0, outra_regra.1
    );
    assert!(outra_regra.0 > 1000, "the differential reads do not exercise the choice between overlapping matches");
    casos.len() as u32
}

fn nome_par(nome: &str) -> &str {
    let n = nome.split_whitespace().next().unwrap_or("");
    n.strip_suffix("/1").or_else(|| n.strip_suffix("/2")).unwrap_or(n)
}

fn main() {
    let dir = std::env::temp_dir().join(format!("teste_fastq_qc_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = |n: &str| dir.join(n).to_string_lossy().to_string();
    let grava = |n: &str, t: &[u8]| {
        std::fs::write(p(n), t).unwrap();
        p(n)
    };
    let erro = |r: Result<FastqStats, String>| r.expect_err("an error was expected");
    let mut s = 20261113u64;
    let mut n = 0u32;

    // ================================================================
    //  PARTE 1: fastq_stats
    // ================================================================

    // ---------------------------------------- arquivo variado, contra a contagem independente
    let variado = gera(
        &mut s,
        &Plano { n: 4000, max_len: 150, zero: true, qmin: 35, qmax: 73, por_mil_n: 15, minusculas: true, adaptadores: true, duplicadas: true, cauda_ruim: true },
    );
    let arq = grava("variado.fastq", &texto(&variado, "\n"));
    let st = fastq_stats(&arq, None).unwrap();
    let c_variado = conta(&variado);
    confere(&st, &c_variado, "varied");
    assert!(st.complete && st.path == arq);
    assert_eq!(st.phred_offset, 33);
    // o adaptador plantado alonga algumas leituras alem das 150 bases
    assert!(st.max_length > 150 && st.max_length < 180 && st.positions_reported == st.max_length);
    assert_eq!((st.read_quality_histogram.len(), st.gc_histogram.len()), (94, 101));
    // o que foi plantado esta la: zeros, N, minusculas, duplicadas, adaptadores
    assert!(st.min_length == 0 && st.reads_with_n > 500 && st.duplication_percent > 8.0, "{}", st.resumo());
    assert!(variado.iter().any(|r| r.seq.iter().any(|b| b.is_ascii_lowercase())));
    assert!(st.adapters[0].percent > 15.0 && st.adapters[1].percent > 3.0 && st.adapters[2].reads > 0 && st.adapters[2].percent < 0.5);
    assert_eq!(st.likely_adapter, Some(("Illumina Universal", ILLUMINA)));
    n += 1;

    // ---------------------------------------- o JSON e JSON, com as chaves em ingles
    let j = json(&st.json());
    assert_eq!(
        j.chaves(),
        [
            "path", "reads", "bases", "complete", "min_length", "max_length", "mean_length", "length_histogram", "phred_offset",
            "quality_char_range", "mean_quality", "read_quality_histogram", "positions_reported", "per_position", "gc_percent", "gc_histogram", "n_bases",
            "n_percent", "reads_with_n", "duplication_percent", "duplication_sample", "adapters", "likely_adapter", "warnings"
        ]
    );
    assert_eq!(j.campo("path"), &J::Txt(arq.clone()));
    assert_eq!((j.campo("reads").num(), j.campo("bases").num()), (st.reads as f64, st.bases as f64));
    assert_eq!(j.campo("complete"), &J::Bool(true));
    assert_eq!(j.campo("phred_offset").num(), 33.0);
    assert_eq!(j.campo("quality_char_range"), &J::Lista(vec![J::Num(35.0), J::Num(73.0)]));
    assert!(perto(j.campo("mean_quality").num(), st.mean_quality) && perto(j.campo("gc_percent").num(), st.gc_percent));
    let hl = j.campo("length_histogram").lista();
    assert_eq!(hl.len(), st.length_histogram.len());
    assert_eq!(hl[0], J::Lista(vec![J::Num(st.length_histogram[0].0 as f64), J::Num(st.length_histogram[0].1 as f64)]));
    assert_eq!((j.campo("read_quality_histogram").lista().len(), j.campo("gc_histogram").lista().len()), (94, 101));
    let pp = j.campo("per_position");
    assert_eq!(pp.chaves(), ["bases", "mean_quality", "q10", "q25", "q50", "q75", "q90", "a", "c", "g", "t", "n"]);
    for k in pp.chaves() {
        assert_eq!(pp.campo(k).lista().len(), st.positions_reported, "per_position.{k}");
    }
    assert_eq!(pp.campo("q50").lista()[149].num(), st.positions[149].q50 as f64);
    assert_eq!(pp.campo("n").lista()[7].num(), st.positions[7].n as f64);
    let ads = j.campo("adapters").lista();
    assert_eq!(ads.len(), 4);
    assert_eq!(ads[0].chaves(), ["name", "sequence", "reads", "percent"]);
    assert_eq!(ads[2].campo("name"), &J::Txt("Illumina Small RNA 3'".to_string()));
    assert_eq!(j.campo("likely_adapter").campo("sequence"), &J::Txt(ILLUMINA.to_string()));
    assert_eq!(j.campo("warnings").lista().len(), st.warnings().len());
    n += 1;

    // ---------------------------------------- o relatorio em texto e os avisos, cada um com o conserto
    let r = st.resumo();
    assert!(r.starts_with(&format!("{arq}: FASTQ with 4000 reads and {} bases\n", st.bases)), "{r}");
    assert!(r.contains("Phred+33") && r.contains("read length") && r.contains("GC content") && r.contains("duplication"), "{r}");
    assert!(r.contains("Adapter found: Illumina Universal") && r.contains("fastq_trim(adapter=\"AGATCGGAAGAGC\")"), "{r}");
    assert!(r.contains("Low quality tail: from cycle") && r.contains("fastq_trim(quality_cutoff=20)"), "{r}");
    assert!(r.contains("Many undetermined bases") && r.contains("fastq_trim(trim_n=True)") && r.contains("fastq_trim(max_n=0)"), "{r}");
    assert!(!r.contains("Phred+64") && !r.contains("no problems found"), "{r}");
    assert!(r.lines().count() >= 12 && r.lines().count() <= 16, "{r}");
    n += 1;

    // ---------------------------------------- arquivo limpo: nenhum aviso
    let limpo = gera(
        &mut s,
        &Plano { n: 1500, max_len: 100, zero: false, qmin: 58, qmax: 73, por_mil_n: 0, minusculas: false, adaptadores: false, duplicadas: false, cauda_ruim: false },
    );
    let st = fastq_stats(&grava("limpo.fq", &texto(&limpo, "\n")), None).unwrap();
    confere(&st, &conta(&limpo), "clean");
    assert!(st.warnings().is_empty() && st.likely_adapter.is_none(), "{}", st.resumo());
    assert!(st.resumo().contains("no problems found"));
    assert_eq!(json(&st.json()).campo("likely_adapter"), &J::Nulo);
    assert_eq!((st.n_bases, st.reads_with_n), (0, 0));
    n += 1;

    // ---------------------------------------- Phred+64, e a fronteira da deteccao
    let antigo = gera(
        &mut s,
        &Plano { n: 1200, max_len: 80, zero: true, qmin: 66, qmax: 104, por_mil_n: 5, minusculas: true, adaptadores: true, duplicadas: true, cauda_ruim: false },
    );
    let st = fastq_stats(&grava("phred64.fastq", &texto(&antigo, "\n")), None).unwrap();
    let c = conta(&antigo);
    assert_eq!((st.phred_offset, c.offset), (64, 64));
    confere(&st, &c, "Phred+64");
    assert!(st.mean_quality > 15.0 && st.mean_quality < 25.0, "{}", st.mean_quality);
    let r = st.resumo();
    assert!(r.contains("quality scale     Phred+64\n") && r.contains("fastq_trim(adapter=\"AGATCGGAAGAGC\")"), "{r}");
    // a escala e o primeiro aviso, e diz o que passar ao aparo
    let w = st.warnings();
    assert!(w[0].contains("old Phred+64 scale") && w[0].contains("from `B` to `h`") && w[0].contains("pass phred_offset=64 to fastq_trim"), "{}", w[0]);
    n += 1;
    // Qualidade constante `I`, como nas leituras simuladas (e o caso do
    // amostra_paciente.fastq do projeto): e Phred+33, qualidade 40, sem aviso
    // nenhum -- e nao Phred+64 de qualidade 9.
    let simulado: Vec<Rec> = (0..50).map(|i| reci(&format!("sim{i}"), &String::from_utf8(bases(&mut s, 80)).unwrap())).collect();
    let st = fastq_stats(&grava("simulado.fastq", &texto(&simulado, "\n")), None).unwrap();
    confere(&st, &conta(&simulado), "constant quality");
    assert_eq!((st.phred_offset, st.quality_char_range, st.positions[0].q50), (33, Some((b'I', b'I')), 40));
    assert!(perto(st.mean_quality, 40.0) && st.read_quality_histogram[40] == 50);
    assert!(st.warnings().is_empty() && st.resumo().contains("quality scale     Phred+33\n") && st.resumo().contains("no problems found"), "{}", st.resumo());
    assert_eq!(json(&st.json()).campo("quality_char_range"), &J::Lista(vec![J::Num(73.0), J::Num(73.0)]));
    n += 1;
    // A fronteira: o menor caractere tem de ser `@` (64) ou mais, e o maior
    // tem de passar de `K` (75).
    for (menor, maior, esperado) in [
        (64u8, 104u8, 64u8),
        (63, 104, 33),
        (59, 104, 33),
        (58, 104, 33),
        (33, 104, 33),
        (64, 76, 64),
        (64, 75, 33),
        (70, 75, 33),
        (75, 75, 33),
        (76, 76, 64),
        (66, 126, 64),
    ] {
        let recs = vec![
            Rec { nome: "a".into(), seq: b"ACGT".to_vec(), qual: vec![maior; 4] },
            Rec { nome: "b".into(), seq: b"ACG".to_vec(), qual: vec![maior, menor, maior] },
        ];
        let st = fastq_stats(&grava("fronteira.fastq", &texto(&recs, "\n")), None).unwrap();
        assert_eq!(st.phred_offset, esperado, "quality characters from {menor} to {maior}");
        confere(&st, &conta(&recs), "offset boundary");
        assert_eq!((st.positions[1].q10, st.positions[1].q90, st.quality_char_range), (menor - esperado, maior - esperado, Some((menor, maior))));
        assert_eq!(st.warnings().len(), usize::from(esperado == 64) + usize::from(maior - esperado < 20), "{}", st.resumo());
        n += 1;
    }

    // ---------------------------------------- leituras de mais de 1000 bases
    let mut longas = gera(
        &mut s,
        &Plano { n: 60, max_len: 2500, zero: false, qmin: 40, qmax: 70, por_mil_n: 10, minusculas: true, adaptadores: false, duplicadas: false, cauda_ruim: false },
    );
    longas.push(Rec { nome: "curta".into(), seq: b"ACGTN".to_vec(), qual: b"IIII#".to_vec() });
    let st = fastq_stats(&grava("longas.fastq", &texto(&longas, "\n")), None).unwrap();
    let c = conta(&longas);
    confere(&st, &c, "long reads");
    assert_eq!((st.positions_reported, st.positions.len(), st.max_length, MAX_POS), (1000, 1000, 2500, 1000));
    // os totais contam as bases alem do ciclo 1000; as tabelas por ciclo, nao
    let nas_tabelas: u64 = st.positions.iter().map(|x| x.bases).sum();
    assert!(st.bases > nas_tabelas && nas_tabelas == longas.iter().map(|r| r.seq.len().min(1000) as u64).sum::<u64>());
    assert!(st.resumo().contains("only the first 1000 cycles are reported"));
    n += 1;

    // ---------------------------------------- .gz, .lz4 e CRLF dao o mesmo relatorio
    let st0 = fastq_stats(&p("variado.fastq"), None).unwrap();
    let mut gz = flate2::write::GzEncoder::new(std::fs::File::create(p("variado.fastq.gz")).unwrap(), flate2::Compression::fast());
    gz.write_all(&texto(&variado, "\n")).unwrap();
    gz.finish().unwrap();
    let mut lz = externo::Escritor::cria(&p("variado.fq.lz4")).unwrap();
    lz.write_all(&texto(&variado, "\n")).unwrap();
    lz.termina().unwrap();
    grava("variado_crlf.fastq", &texto(&variado, "\r\n"));
    for arq in ["variado.fastq.gz", "variado.fq.lz4", "variado_crlf.fastq"] {
        let st = fastq_stats(&p(arq), None).unwrap();
        confere(&st, &c_variado, arq);
        assert_eq!(st.json().replace(&nomes::json_texto(&p(arq)), "X"), st0.json().replace(&nomes::json_texto(&p("variado.fastq")), "X"), "{arq}");
        n += 1;
    }

    // ---------------------------------------- max_reads
    for k in [0usize, 1, 37, 3999, 4000, 4001] {
        let st = fastq_stats(&p("variado.fastq"), Some(k as u64)).unwrap();
        let ate = k.min(variado.len());
        confere(&st, &conta(&variado[..ate]), "max_reads");
        assert_eq!(st.complete, k >= 4000, "max_reads={k}");
        assert_eq!(st.resumo().contains("only the first"), k < 4000);
        n += 1;
    }
    // linhas vazias depois do ultimo registro nao sao "mais conteudo"
    let mut t = texto(&variado[..5], "\n");
    t.extend_from_slice(b"\n\n");
    assert!(fastq_stats(&grava("cinco.fastq", &t), Some(5)).unwrap().complete);
    assert!(!fastq_stats(&p("cinco.fastq"), Some(4)).unwrap().complete);
    n += 1;

    // ---------------------------------------- arquivo vazio: zero leituras, sem erro
    for (nome, conteudo) in [("vazio.fastq", &b""[..]), ("so_linhas.fastq", &b"\n\r\n\n"[..])] {
        for limite in [None, Some(0), Some(10)] {
            let st = fastq_stats(&grava(nome, conteudo), limite).unwrap();
            confere(&st, &conta(&[]), "empty");
            assert!(st.complete && st.reads == 0 && st.positions.is_empty() && st.length_histogram.is_empty());
            assert_eq!((st.mean_length, st.mean_quality, st.gc_percent, st.n_percent, st.duplication_percent), (0.0, 0.0, 0.0, 0.0, 0.0));
            assert!(st.warnings().is_empty() && st.resumo().contains("0 reads") && st.quality_char_range.is_none());
            assert_eq!(json(&st.json()).campo("quality_char_range"), &J::Nulo);
        }
        n += 1;
    }

    // ---------------------------------------- casos feitos a mao
    // linhas vazias entre registros; leitura de tamanho zero; ultimo registro sem fim de linha
    let st = fastq_stats(&grava("mao.fastq", b"@r1 x\nACGTN\n+r1\nIIII#\n\n\n@r2\n\n+\n\n@r3\nggcc\n+\n!!!~"), None).unwrap();
    assert_eq!((st.reads, st.bases, st.min_length, st.max_length), (3, 9, 0, 5));
    assert_eq!(st.length_histogram, vec![(0, 1), (4, 1), (5, 1)]);
    assert_eq!((st.n_bases, st.reads_with_n), (1, 1));
    assert!(perto(st.gc_percent, 100.0 * 6.0 / 8.0));
    // GC por leitura: r1 = 2 de 4 -> 50; r3 = 4 de 4 -> 100; r2 nao tem base
    assert_eq!((st.gc_histogram[50], st.gc_histogram[100], st.gc_histogram.iter().sum::<u64>()), (1, 1, 2));
    // qualidade media por leitura: r1 = (40*4+2)/5 = 32,4 -> 32; r3 = (0+0+0+93)/4 = 23,25 -> 23
    assert_eq!((st.read_quality_histogram[32], st.read_quality_histogram[23], st.read_quality_histogram.iter().sum::<u64>()), (1, 1, 2));
    assert!(perto(st.mean_quality, (162.0 + 93.0) / 9.0));
    let p1 = &st.positions[0];
    assert_eq!((p1.bases, p1.a, p1.c, p1.g, p1.t, p1.n, p1.q10, p1.q50, p1.q90), (2, 1, 0, 1, 0, 0, 0, 0, 40));
    let p4 = &st.positions[3];
    assert_eq!((p4.bases, p4.c, p4.t, p4.q10, p4.q25, p4.q50, p4.q75, p4.q90), (2, 1, 1, 40, 40, 40, 93, 93));
    assert_eq!((st.positions[4].bases, st.positions[4].n, st.positions[4].q50), (1, 1, 2));
    n += 1;
    // percentis exatos: qualidades 1..=10 num ciclo -> postos 1, 3, 5, 8, 9
    let dez: Vec<Rec> = (1..=10u8).rev().map(|q| Rec { nome: format!("q{q}"), seq: b"A".to_vec(), qual: vec![33 + q] }).collect();
    let st = fastq_stats(&grava("dez.fastq", &texto(&dez, "\n")), None).unwrap();
    let x = &st.positions[0];
    assert_eq!((x.q10, x.q25, x.q50, x.q75, x.q90), (1, 3, 5, 8, 9));
    assert!(perto(x.mean_quality, 5.5));
    n += 1;
    // duplicacao: 3 de 6 leituras (b, d, f) repetem uma sequencia ja vista; maiusculas e minusculas diferem
    let dup = [rec("a", "ACGT", "IIII"), rec("b", "ACGT", "####"), rec("c", "acgt", "IIII"), rec("d", "ACGT", "IIII"), rec("e", "", ""), rec("f", "", "")];
    let st = fastq_stats(&grava("dup.fastq", &texto(&dup, "\n")), None).unwrap();
    assert_eq!(st.duplication_sample, 6);
    assert!(perto(st.duplication_percent, 100.0 * 3.0 / 6.0), "{}", st.duplication_percent);
    confere(&st, &conta(&dup), "duplicates");
    n += 1;
    // adaptador: so as 12 primeiras bases, exatas, em qualquer lugar, sem diferenca de caixa
    let ad = [
        reci("a", "TTTTAGATCGGAAGAGTTTT"),     // 12 bases do Illumina: conta
        reci("b", "TTTTAGATCGGAAGATTTTT"),     // 11 bases: nao conta
        reci("c", "ttttagatcggaagagcacac"),    // minusculas: conta
        reci("d", "AGATCGNAAGAGC"),            // um N no meio: nao conta
        reci("e", "CTGTCTCTTATAAGATCGGAAGAG"), // Nextera (12) e Illumina (12)
    ];
    let st = fastq_stats(&grava("ad.fastq", &texto(&ad, "\n")), None).unwrap();
    assert_eq!(st.adapters.iter().map(|a| a.reads).collect::<Vec<_>>(), [3, 1, 0, 0]);
    assert!(perto(st.adapters[0].percent, 60.0) && perto(st.adapters[1].percent, 20.0));
    assert_eq!(st.likely_adapter, Some(KNOWN_ADAPTERS[0]));
    confere(&st, &conta(&ad), "adapters");
    n += 1;

    // ---------------------------------------- a amostra de duplicacao para em DUP_SAMPLE
    {
        let total = DUP_SAMPLE as usize + 300;
        let arq = p("milhao.fastq");
        let mut w = std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(&arq).unwrap());
        let mut codigos: Vec<u64> = Vec::with_capacity(total);
        let mut vistos: HashSet<u64> = HashSet::new();
        let (mut repetidas, mut linha) = (0u64, [0u8; 20]);
        let mut qualidade = vec![b'5'; 19];
        qualidade.extend_from_slice(b"I\n");
        for i in 0..total {
            // depois da amostra, so repeticoes: nao podem entrar na conta
            let cod = if i >= DUP_SAMPLE as usize || (i > 0 && rnd(&mut s, 20) == 0) { codigos[rnd(&mut s, i)] } else { lcg(&mut s) & ((1 << 40) - 1) };
            codigos.push(cod);
            if i < DUP_SAMPLE as usize && !vistos.insert(cod) {
                repetidas += 1;
            }
            for (k, b) in linha.iter_mut().enumerate() {
                *b = b"ACGT"[((cod >> (2 * k)) & 3) as usize];
            }
            w.write_all(b"@d\n").unwrap();
            w.write_all(&linha).unwrap();
            w.write_all(b"\n+\n").unwrap();
            w.write_all(&qualidade).unwrap();
        }
        w.flush().unwrap();
        drop(w);
        let st = fastq_stats(&arq, None).unwrap();
        assert_eq!((st.reads, st.bases, st.duplication_sample), (total as u64, 20 * total as u64, DUP_SAMPLE));
        assert!(perto(st.duplication_percent, 100.0 * repetidas as f64 / DUP_SAMPLE as f64), "{} x {repetidas}", st.duplication_percent);
        assert!(repetidas > 40_000 && repetidas < 60_000);
        assert_eq!((st.min_length, st.max_length, st.length_histogram.clone()), (20, 20, vec![(20, total as u64)]));
        assert_eq!(st.positions[19].q50, 40);
        std::fs::remove_file(&arq).unwrap();
        n += 1;
    }

    // ---------------------------------------- FASTQ invalido: erro com a linha e o que se esperava
    let bom = "@r1\nACGT\n+\nIIII\n";
    for (nome, conteudo, linha, trechos) in [
        ("sem_arroba.fastq", "r1\nACGT\n+\nIIII\n".to_string(), "line 1:", vec!["must start with `@`", "record 1", "`r1`"]),
        ("sem_arroba2.fastq", format!("{bom}{bom}ACGT\n+\nIIII\n"), "line 9:", vec!["must start with `@`", "record 3", "`ACGT`"]),
        ("sem_mais.fastq", format!("{bom}@r2\nACGT\n-\nIIII\n"), "line 7:", vec!["must start with `+`", "record 2", "`r2`", "`-`"]),
        ("multilinha.fastq", format!("{bom}@r2\nACGT\nACGT\n+\nIIIIIIII\n"), "line 7:", vec!["must start with `+`", "wrapped over several lines"]),
        ("tamanhos.fastq", format!("{bom}@r2\nACGT\n+\nIII\n"), "line 8:", vec!["record 2", "`r2`", "4 bases but 3 quality values"]),
        ("corta1.fastq", format!("{bom}@r2\n"), "line 5,", vec!["ends at", "middle of record 2", "`r2`", "the sequence line"]),
        ("corta2.fastq", format!("{bom}@r2\nAC\n"), "line 6,", vec!["middle of record 2", "the `+` separator line", "fetch the file again"]),
        ("corta3.fastq", format!("{bom}{bom}@r3 x\nAC\n+\n"), "line 11,", vec!["middle of record 3", "`r3`", "the quality line"]),
        ("corta4.fastq", format!("{bom}@r2\nAC\n+"), "line 7,", vec!["middle of record 2", "the quality line"]),
        ("qual_ruim.fastq", format!("{bom}@r2\nACGT\n+\nII I\n"), "line 8:", vec!["record 2", "column 3", "ASCII code 32"]),
        ("qual_ruim2.fastq", "@r1\nAC\n+\nI\u{7f}\n".to_string(), "line 4:", vec!["column 2", "ASCII code 127"]),
    ] {
        let e = erro(fastq_stats(&grava(nome, conteudo.as_bytes()), None));
        assert!(e.contains(&p(nome)) && e.contains(linha), "{nome}: {e}");
        for t in trechos {
            assert!(e.contains(t), "{nome}: missing `{t}` in: {e}");
        }
        // com CRLF a linha apontada e a mesma
        let e2 = erro(fastq_stats(&grava(nome, conteudo.replace('\n', "\r\n").as_bytes()), None));
        assert_eq!(e2, e, "{nome} with CRLF");
        n += 1;
    }
    // o que parece truncado e nao e: leitura de tamanho zero no fim, sem a ultima linha vazia
    let st = fastq_stats(&grava("zero_no_fim.fastq", format!("{bom}@r2\n\n+\n").as_bytes()), None).unwrap();
    assert_eq!((st.reads, st.bases, st.min_length), (2, 4, 0));
    // max_reads para antes do defeito: nao ha erro
    assert_eq!(fastq_stats(&p("tamanhos.fastq"), Some(1)).unwrap().reads, 1);
    assert!(erro(fastq_stats(&p("nao_existe.fastq"), None)).contains("nao_existe.fastq"));
    n += 1;

    // ================================================================
    //  PARTE 2: aparo
    // ================================================================

    // ---------------------------------------- aparo por qualidade: casos feitos a mao
    let ph = |v: &[u8]| -> Vec<u8> { v.iter().map(|q| q + 33).collect() };
    for (q, c5, c3, esperado) in [
        (vec![30, 30, 30, 30, 10, 10], 0, 20, (0, 4)),
        (vec![30, 30, 10, 30, 10, 10], 0, 20, (0, 4)),
        (vec![10, 30, 30, 2, 30], 0, 20, (0, 5)),  // base ruim no meio, ponta boa: nada sai
        (vec![5, 5, 5], 0, 20, (0, 0)),
        (vec![10, 10, 30, 30, 30, 30], 20, 0, (2, 6)),
        (vec![10, 10, 10, 10], 20, 20, (0, 0)),    // os cortes se cruzam
        (vec![5, 30, 30, 30, 5], 20, 20, (1, 4)),
        (vec![5, 30, 30, 30, 5], 0, 0, (0, 5)),    // corte zero desliga
        (vec![5, 30, 30, 30, 5], 20, 0, (1, 5)),
        (vec![20, 20, 20], 20, 20, (0, 3)),        // qualidade igual ao corte fica
        (vec![40, 40, 3, 40, 3], 0, 10, (0, 4)),
        (vec![40, 40, 3, 3, 40, 3, 3], 0, 10, (0, 5)),  // BWA nao volta ate o primeiro trecho ruim
        (vec![], 20, 20, (0, 0)),
    ] {
        assert_eq!(quality_trim_bounds(&ph(&q), c5, c3, 33), esperado, "q={q:?} 5'={c5} 3'={c3}");
        n += 1;
    }
    // a escala muda o resultado: 'B' e 2 em Phred+64 e 33 em Phred+33
    assert_eq!(quality_trim_bounds(b"hhhhBB", 0, 20, 64), (0, 4));
    assert_eq!(quality_trim_bounds(b"hhhhBB", 0, 20, 33), (0, 6));
    n += 1;

    // ---------------------------------------- aparo por qualidade: aleatorio contra a transcricao dos lacos
    let (mut cortou5, mut cortou3, mut tudo) = (0, 0, 0);
    for i in 0..60_000 {
        let base = if i % 3 == 0 { 64u8 } else { 33 };
        let len = rnd(&mut s, 70);
        let jeito = rnd(&mut s, 4);
        let qual: Vec<u8> = (0..len)
            .map(|k| {
                base + match jeito {
                    0 => rnd(&mut s, 42),
                    1 => if k * 2 < len { 20 + rnd(&mut s, 21) } else { rnd(&mut s, 25) },
                    2 => if k * 3 < len { rnd(&mut s, 20) } else { 15 + rnd(&mut s, 27) },
                    _ => [2, 12, 23, 37][rnd(&mut s, 4)],
                } as u8
            })
            .collect();
        let (c5, c3) = (if rnd(&mut s, 3) == 0 { 0 } else { rnd(&mut s, 46) as u8 }, if rnd(&mut s, 5) == 0 { 0 } else { rnd(&mut s, 46) as u8 });
        let q: Vec<i64> = qual.iter().map(|&c| c as i64 - base as i64).collect();
        let esperado = ref_qualidade(&q, c5 as i64, c3 as i64);
        let obtido = quality_trim_bounds(&qual, c5, c3, base);
        assert_eq!(obtido, esperado, "qualities {q:?}, cutoffs 5'={c5} 3'={c3}");
        if c5 == 0 && c3 > 0 {
            assert_eq!(obtido, (0, ref_qualidade_3_bruta(&q, c3 as i64)), "{q:?} 3'={c3}");
        }
        cortou5 += (obtido.0 > 0) as u32;
        cortou3 += (obtido.1 < len && obtido != (0, 0)) as u32;
        tudo += (obtido == (0, 0) && len > 0) as u32;
    }
    // os tres desfechos aparecem: so a ponta 5', so a 3', e nada sobrando
    assert!(cortou5 > 5_000 && cortou3 > 5_000 && tudo > 5_000, "{cortou5} {cortou3} {tudo}");
    n += 1;

    // ---------------------------------------- adaptador: casos feitos a mao
    let acha = |ad: &str, min: usize, taxa: f64, leitura: &str| -> Option<usize> {
        let r = Adapter::new(ad, min, taxa).unwrap().find(leitura.as_bytes());
        assert_eq!(r, ref_adaptador(ad, leitura.as_bytes(), min, taxa), "adapter {ad}, read {leitura}");
        r
    };
    for (ad, min, taxa, leitura, esperado) in [
        (ILLUMINA, 3, 0.1, "ACGTACGTAGATCGGAAGAGCTTTT", Some(8)),   // inteiro, no meio
        (ILLUMINA, 3, 0.1, "ACGTACGTAGATCGGAAGAGC", Some(8)),       // inteiro, no fim
        (ILLUMINA, 3, 0.1, "AGATCGGAAGAGCACGT", Some(0)),           // no comeco: nao sobra nada
        (ILLUMINA, 3, 0.1, "ACGTACGTAGA", Some(8)),                 // 3 bases na ponta 3'
        (ILLUMINA, 3, 0.1, "ACGTACGTAG", None),                     // 2 bases: abaixo de min_overlap
        (ILLUMINA, 2, 0.1, "ACGTACGTAG", Some(8)),
        (ILLUMINA, 1, 0.1, "CCCCCCCA", Some(7)),
        (ILLUMINA, 3, 0.1, "ACGTAGATCGGTAGAGC", Some(4)),           // 1 troca em 13: floor(1,3) = 1
        (ILLUMINA, 3, 0.0, "ACGTAGATCGGTAGAGC", None),              // sem tolerancia
        (ILLUMINA, 3, 0.1, "ACGTAGTTCGGTAGAGC", None),              // 2 trocas em 13
        (ILLUMINA, 3, 0.2, "ACGTAGTTCGGTAGAGC", Some(4)),           // floor(2,6) = 2
        (ILLUMINA, 3, 0.1, "ACGTAGATCGGNAGAGC", Some(4)),           // N na leitura e UMA troca
        (ILLUMINA, 3, 0.0, "ACGTAGATCGGNAGAGC", None),
        (ILLUMINA, 3, 0.1, "acgtacgtagatcggaagagc", Some(8)),       // minusculas
        (ILLUMINA, 3, 0.1, "TTAGATCGGAAGAGCGGAGATCGGAAGAGC", Some(2)), // duas copias: a da esquerda
        (ILLUMINA, 3, 0.1, "AGATC", Some(0)),                       // leitura mais curta que o adaptador
        (ILLUMINA, 3, 0.1, "AG", None),                             // leitura mais curta que min_overlap
        (ILLUMINA, 3, 0.1, "", None),
        (ILLUMINA, 3, 0.1, "CCCCCCCCCCCCCCCCCCCC", None),
        ("AGNTC", 5, 0.0, "TTTTAGGTCAA", Some(4)),                  // N do adaptador casa com qualquer base
        ("AGNTC", 5, 0.0, "TTTTAGNTC", Some(4)),                    // ... inclusive com N
        ("AGNTC", 5, 0.0, "TTTTANGTC", None),                       // N da leitura contra G: troca
        ("agatcggaagagc", 3, 0.1, "ACGTACGTAGATCGGAAGAGC", Some(8)), // adaptador em minusculas
        ("TGC", 10, 0.0, "AAATGCAAA", Some(3)),                     // min_overlap maior que o adaptador: vale o adaptador inteiro
        ("TGC", 10, 0.0, "AAAAAATG", None),
        // Adaptador repetitivo: a posicao 4 ja e aceitavel (`GT` como duas
        // trocas), mas as ocorrencias seguintes, sobrepostas, pontuam mais, e a
        // da posicao 6 e exata: o corte fica no comeco verdadeiro da cauda.
        ("AAAAAAAAAA", 3, 0.2, "CCCCGTAAAAAAAAAAAA", Some(6)),
        ("AAAAAAAAAA", 3, 0.2, "CCCCGTAAAAAAAAA", Some(6)),          // o mesmo com a cauda mais curta que o adaptador
        // Ocorrencia com 1 troca na posicao 4 e outra, exata, depois do `T`:
        // a exata vence se comeca a ate meio adaptador (5) da que ja ficou...
        ("AAAAAAAAAA", 3, 0.2, "CCCCAAAATAAAAAAAAAAAA", Some(9)),
        // ... e perde se comeca mais longe: fica a ocorrencia da esquerda.
        ("AAAAAAAAAA", 3, 0.2, "CCCCAAAAATAAAAAAAAAAAA", Some(4)),
        ("GGGGGGGGGG", 3, 0.1, "ACGTACGTGGGGGGGGGGGGGG", Some(8)),
        // N no adaptador nao conta no comprimento: 14 bases com 6 N valem 8, e
        // floor(0,1 * 8) = 0 trocas; com taxa 0,125, floor(1,0) = 1.
        ("TGCANNNNNNGACT", 3, 0.1, "CCCCTGCTGGGGGGGACT", None),
        ("TGCANNNNNNGACT", 3, 0.125, "CCCCTGCTGGGGGGGACT", Some(4)),
        ("TGCANNNNNNGACT", 3, 0.1, "CCCCTGCAGGGGGGGACT", Some(4)),
        // N no comeco do adaptador: com min_overlap 2, os dois N casam com as
        // duas ultimas bases de qualquer leitura
        ("NNACGT", 3, 0.0, "TTTTTTTT", None),
        ("NNACGT", 2, 0.0, "TTTTTTTT", Some(6)),
    ] {
        assert_eq!(acha(ad, min, taxa, leitura), esperado, "adapter {ad}, min_overlap {min}, rate {taxa}, read {leitura}");
        n += 1;
    }
    // sobreposicao de 10 bases com 1 troca: floor(0,1 * 10) = 1, aceita; com 9 bases, floor(0,9) = 0
    assert_eq!(acha(ILLUMINA, 3, 0.1, "ACGTAGATCGGTAG"), Some(4));
    assert_eq!(acha(ILLUMINA, 3, 0.1, "ACGTAGATCGGTA"), None);
    assert_eq!(acha(ILLUMINA, 3, 0.1, "ACGTAGATCGGAA"), Some(4));
    n += 1;

    // ---------------------------------------- adaptador: aleatorio contra a busca ingenua
    let mut ads_teste: Vec<String> = [ILLUMINA, "CTGTCTCTTATACACATCT", "ACGNNACGTTGCA", "NNACGTACGTAC", "AAAAAAAAAAAAAAA", "ACACACACACAC", "TGC", "GATTACA", "A"]
        .iter()
        .map(|a| a.to_string())
        .collect();
    for _ in 0..6 {
        let k = 4 + rnd(&mut s, 30);
        ads_teste.push(String::from_utf8(bases(&mut s, k)).unwrap());
    }
    let (mut casos, mut achados, mut nao_achados, mut no_fim) = (0u32, 0u32, 0u32, 0u32);
    for ad in &ads_teste {
        for taxa in [0.0, 0.05, 0.1, 0.2, 0.34] {
            for min in [1usize, 3, 5, 8, 40] {
                let a = Adapter::new(ad, min, taxa).unwrap();
                for i in 0..350 {
                    let leitura = leitura_com_adaptador(&mut s, ad.as_bytes(), i % 7);
                    let (obtido, esperado) = (a.find(&leitura), ref_adaptador(ad, &leitura, min, taxa));
                    assert_eq!(obtido, esperado, "adapter {ad}, min_overlap {min}, rate {taxa}, read {}", String::from_utf8_lossy(&leitura));
                    casos += 1;
                    match obtido {
                        Some(pos) => {
                            achados += 1;
                            no_fim += (pos + ad.len() > leitura.len()) as u32;
                        }
                        None => nao_achados += 1,
                    }
                }
                n += 1;
            }
        }
    }
    assert!(casos == 15 * 25 * 350 && achados > 30_000 && nao_achados > 15_000 && no_fim > 5_000, "{casos} {achados} {nao_achados} {no_fim}");

    // ---------------------------------------- adaptador e opcoes invalidos
    for (ad, min, taxa, trecho) in [
        ("AGXT", 3, 0.1, "invalid character `X`"),
        ("auto", 3, 0.1, "invalid character `U`"),
        ("ACG-T", 3, 0.1, "invalid character `-`"),
        ("", 3, 0.1, "is empty"),
        ("NNNN", 3, 0.1, "only N"),
        (ILLUMINA, 0, 0.1, "min_overlap is 0"),
        (ILLUMINA, 3, 1.0, "max_error_rate is 1"),
        (ILLUMINA, 3, -0.1, "max_error_rate is -0.1"),
        (ILLUMINA, 3, f64::NAN, "max_error_rate is NaN"),
    ] {
        let e = Adapter::new(ad, min, taxa).expect_err("an error was expected");
        assert!(e.contains(trecho), "{ad}: {e}");
        n += 1;
    }
    let padrao = TrimOptions::default();
    assert_eq!((padrao.adapter.clone(), padrao.adapter2.clone(), padrao.min_overlap, padrao.max_error_rate), (None, None, 3, 0.1));
    assert_eq!((padrao.quality_cutoff, padrao.quality_cutoff_5, padrao.trim_n, padrao.min_length), (0, 0, false, 0));
    assert_eq!((padrao.max_n, padrao.min_mean_quality, padrao.phred_offset), (None, None, 33));
    n += 1;

    // ---------------------------------------- fastq_trim: sem nenhuma opcao, copia (com LF e `+` sozinho)
    let apara = |entrada: &str, saida: &str, o: &TrimOptions| fastq_trim(&p(entrada), &p(saida), None, None, o).unwrap();
    let s0 = apara("variado.fastq", "copia.fastq", &padrao);
    assert_eq!(le(&p("copia.fastq")), texto(&variado, "\n"));
    assert_eq!((s0.paired, s0.reads_in, s0.reads_out, s0.bases_in, s0.bases_out), (false, 4000, 4000, st0.bases, st0.bases));
    assert_eq!(s0, TrimSummary { reads_in: 4000, reads_out: 4000, bases_in: st0.bases, bases_out: st0.bases, ..TrimSummary::default() });
    let j = json(&s0.json());
    assert_eq!(
        j.chaves(),
        [
            "paired", "reads_in", "reads_out", "bases_in", "bases_out", "reads_with_adapter", "reads_with_adapter2", "bases_removed_quality",
            "bases_removed_adapter", "bases_removed_n", "discarded_too_short", "discarded_too_many_n", "discarded_low_quality"
        ]
    );
    assert_eq!((j.campo("paired"), j.campo("reads_out").num(), j.campo("bases_in").num()), (&J::Bool(false), 4000.0, st0.bases as f64));
    n += 1;
    // CRLF e `+nome` na entrada; LF e `+` sozinho na saida
    grava("crlf.fastq", b"@r1 a b\r\nACGT\r\n+r1 a b\r\nIIII\r\n@r2\r\n\r\n+\r\n\r\n");
    apara("crlf.fastq", "crlf_saida.fastq", &padrao);
    assert_eq!(le(&p("crlf_saida.fastq")), b"@r1 a b\nACGT\n+\nIIII\n@r2\n\n+\n\n");
    n += 1;

    // ---------------------------------------- a ordem dos passos
    // (a) qualidade ANTES do adaptador: a cauda ruim leva 10 das 13 bases do
    //     adaptador; as 3 que sobram ainda sao achadas (min_overlap 3).
    grava("ordem_a.fastq", &texto(&[rec("a", "ACGTACGTACAGATCGGAAGAGC", "IIIIIIIIIIIII##########")], "\n"));
    let o = TrimOptions { adapter: Some(ILLUMINA.into()), quality_cutoff: 20, ..TrimOptions::default() };
    let r = apara("ordem_a.fastq", "ordem_a_saida.fastq", &o);
    assert_eq!(le(&p("ordem_a_saida.fastq")), b"@a\nACGTACGTAC\n+\nIIIIIIIIII\n");
    assert_eq!((r.bases_removed_quality, r.bases_removed_adapter, r.reads_with_adapter, r.bases_in, r.bases_out), (10, 3, 1, 23, 10));
    //     com min_overlap 4 as 3 bases que sobraram nao bastam: ficam na leitura
    let r = apara("ordem_a.fastq", "ordem_a4.fastq", &TrimOptions { min_overlap: 4, ..o.clone() });
    assert_eq!(le(&p("ordem_a4.fastq")), b"@a\nACGTACGTACAGA\n+\nIIIIIIIIIIIII\n");
    assert_eq!((r.bases_removed_quality, r.bases_removed_adapter, r.reads_with_adapter), (10, 0, 0));
    // (b) adaptador ANTES do aparo de N: os N so ficam na ponta depois que o adaptador sai
    grava("ordem_b.fastq", &texto(&[reci("b", "ACGTNNAGATCGGAAGAGC")], "\n"));
    let o = TrimOptions { adapter: Some(ILLUMINA.into()), trim_n: true, ..TrimOptions::default() };
    let r = apara("ordem_b.fastq", "ordem_b_saida.fastq", &o);
    assert_eq!(le(&p("ordem_b_saida.fastq")), b"@b\nACGT\n+\nIIII\n");
    assert_eq!((r.bases_removed_adapter, r.bases_removed_n, r.bases_removed_quality), (13, 2, 0));
    // (c) qualidade ANTES do aparo de N, e os filtros por ultimo, sobre o que sobrou:
    //     a leitura inteira tem 4 N e qualidade media baixa; aparada, tem 0 N e media 40
    grava("ordem_c.fastq", &texto(&[rec("c", "NNACGTACGTNNACGT", "IIIIIIIIII######")], "\n"));
    let o = TrimOptions { quality_cutoff: 20, trim_n: true, max_n: Some(0), min_mean_quality: Some(35), min_length: 8, ..TrimOptions::default() };
    let r = apara("ordem_c.fastq", "ordem_c_saida.fastq", &o);
    assert_eq!(le(&p("ordem_c_saida.fastq")), b"@c\nACGTACGT\n+\nIIIIIIII\n");
    assert_eq!((r.bases_removed_quality, r.bases_removed_n, r.reads_out, r.bases_out), (6, 2, 1, 8));
    n += 3;

    // ---------------------------------------- trim_n
    let ns = [
        rec("a", "NNACGTNN", "IIIIIIII"),
        rec("b", "ACNNGT", "IIIIII"),   // N no meio fica
        rec("c", "NNNN", "IIII"),       // so N: sobra uma leitura vazia, que e gravada
        rec("d", "nnACGTn", "IIIIIII"), // n minusculo nao e aparado (como no cutadapt)
        rec("e", "", ""),
        rec("f", "ACGT", "IIII"),
    ];
    grava("ns.fastq", &texto(&ns, "\n"));
    let r = apara("ns.fastq", "ns_saida.fastq", &TrimOptions { trim_n: true, ..TrimOptions::default() });
    assert_eq!(le(&p("ns_saida.fastq")), b"@a\nACGT\n+\nIIII\n@b\nACNNGT\n+\nIIIIII\n@c\n\n+\n\n@d\nnnACGTn\n+\nIIIIIII\n@e\n\n+\n\n@f\nACGT\n+\nIIII\n");
    assert_eq!((r.bases_removed_n, r.reads_in, r.reads_out, r.bases_in, r.bases_out), (8, 6, 6, 29, 21));
    // ... mas conta para max_n: com max_n 0 saem `b` (N no meio) e `d` (os tres n)
    let r = apara("ns.fastq", "ns_max.fastq", &TrimOptions { trim_n: true, max_n: Some(0), ..TrimOptions::default() });
    assert_eq!(registros(&le(&p("ns_max.fastq"))).iter().map(|x| x.nome.clone()).collect::<Vec<_>>(), ["a", "c", "e", "f"]);
    assert_eq!((r.discarded_too_many_n, r.reads_out), (2, 4));
    // sem trim_n nada muda
    assert_eq!(apara("ns.fastq", "ns_igual.fastq", &padrao).bases_removed_n, 0);
    assert_eq!(le(&p("ns_igual.fastq")), texto(&ns, "\n"));
    // uma leitura vazia e reprovada por min_mean_quality, qualquer que seja o limite
    let r = apara("ns.fastq", "ns_q.fastq", &TrimOptions { trim_n: true, min_mean_quality: Some(0), ..TrimOptions::default() });
    assert_eq!(le(&p("ns_q.fastq")), b"@a\nACGT\n+\nIIII\n@b\nACNNGT\n+\nIIIIII\n@d\nnnACGTn\n+\nIIIIIII\n@f\nACGT\n+\nIIII\n");
    assert_eq!((r.discarded_low_quality, r.discarded_too_short, r.discarded_too_many_n, r.reads_out), (2, 0, 0, 4));
    n += 3;

    // ---------------------------------------- cada filtro, e a precedencia entre eles
    let filtros = [
        rec("curta", "ACG", "III"),
        rec("muitos_n", "ACNNGT", "IIIIII"),
        rec("baixa", "ACGTAC", "######"),
        rec("curta_n_baixa", "NN", "##"),        // os tres motivos: conta como curta
        rec("n_baixa", "ACNNGT", "######"),      // dois motivos: conta como N demais
        rec("boa", "ACGTAC", "IIIIII"),
        rec("no_limite", "ACGTNC", "555555"),    // 1 N (<= 1) e media exatamente 20: fica
        rec("quase", "ACGTAC", "555554"),        // media 19,83: sai
        rec("cinco", "ACGTA", "IIIII"),          // comprimento exatamente 5: fica
    ];
    grava("filtros.fastq", &texto(&filtros, "\n"));
    let o = TrimOptions { min_length: 5, max_n: Some(1), min_mean_quality: Some(20), ..TrimOptions::default() };
    let r = apara("filtros.fastq", "filtros_saida.fastq", &o);
    assert_eq!(registros(&le(&p("filtros_saida.fastq"))).iter().map(|x| x.nome.clone()).collect::<Vec<_>>(), ["boa", "no_limite", "cinco"]);
    assert_eq!((r.discarded_too_short, r.discarded_too_many_n, r.discarded_low_quality, r.reads_in, r.reads_out), (2, 2, 2, 9, 3));
    assert_eq!((r.bases_in, r.bases_out), (46, 17));
    // um filtro de cada vez
    for (o, ficam, descartes) in [
        (TrimOptions { min_length: 5, ..TrimOptions::default() }, 7, (2, 0, 0)),
        (TrimOptions { max_n: Some(1), ..TrimOptions::default() }, 6, (0, 3, 0)),
        (TrimOptions { max_n: Some(0), ..TrimOptions::default() }, 5, (0, 4, 0)),
        (TrimOptions { min_mean_quality: Some(20), ..TrimOptions::default() }, 5, (0, 0, 4)),
        (TrimOptions { min_mean_quality: Some(41), ..TrimOptions::default() }, 0, (0, 0, 9)),
    ] {
        let r = apara("filtros.fastq", "filtros_um.fastq", &o);
        assert_eq!((r.reads_out, (r.discarded_too_short, r.discarded_too_many_n, r.discarded_low_quality)), (ficam, descartes), "{o:?}");
        assert_eq!(registros(&le(&p("filtros_um.fastq"))).len() as u64, ficam);
        n += 1;
    }
    n += 1;

    // ---------------------------------------- leitura simples, aleatorio, contra a referencia
    let ads_aparo = [ILLUMINA, "CTGTCTCTTATACACATCT", "ACGNNACGTTGCA", "AAAAAAAAAAAA", "TGC"];
    let mut soma = TrimSummary::default();
    for rodada in 0..48 {
        let base = if rodada % 6 == 5 { 64 } else { 33 };
        let recs = gera_aparo(&mut s, 400, &ads_aparo, base, "s");
        let mut o = opcoes_aleatorias(&mut s, &ads_aparo, false);
        o.phred_offset = base;
        // entrada e saida comprimidas em algumas rodadas
        let (entrada, saida) = match rodada % 4 {
            0 => ("ale.fastq.gz", "ale_saida.fastq.gz"),
            1 => ("ale.fastq", "ale_saida.fq.lz4"),
            _ => ("ale.fastq", "ale_saida.fastq"),
        };
        let mut w = externo::Escritor::cria(&p(entrada)).unwrap();
        w.write_all(&texto(&recs, if rodada % 5 == 0 { "\r\n" } else { "\n" })).unwrap();
        w.termina().unwrap();
        let r = apara(entrada, saida, &o);
        let (esperado, resumo) = ref_simples(&recs, &o);
        assert_eq!(r, resumo, "round {rodada}: {o:?}");
        assert!(le(&p(saida)) == esperado, "round {rodada}: the output differs from the reference; {o:?}");
        assert_eq!(registros(&esperado).len() as u64, r.reads_out);
        // o que entrou = o que saiu + o que cada passo tirou + as leituras descartadas
        assert!(r.bases_in >= r.bases_out + r.bases_removed_quality + r.bases_removed_adapter + r.bases_removed_n);
        assert_eq!(r.reads_in, r.reads_out + r.discarded_too_short + r.discarded_too_many_n + r.discarded_low_quality);
        for (a, b) in [
            (&mut soma.reads_with_adapter, r.reads_with_adapter),
            (&mut soma.bases_removed_quality, r.bases_removed_quality),
            (&mut soma.bases_removed_adapter, r.bases_removed_adapter),
            (&mut soma.bases_removed_n, r.bases_removed_n),
            (&mut soma.discarded_too_short, r.discarded_too_short),
            (&mut soma.discarded_too_many_n, r.discarded_too_many_n),
            (&mut soma.discarded_low_quality, r.discarded_low_quality),
            (&mut soma.reads_out, r.reads_out),
        ] {
            *a += b;
        }
        n += 1;
    }
    // todos os caminhos foram exercitados
    assert!(soma.reads_with_adapter > 2000 && soma.bases_removed_quality > 20_000 && soma.bases_removed_adapter > 20_000, "{soma:?}");
    assert!(soma.bases_removed_n > 500 && soma.discarded_too_short > 500 && soma.discarded_too_many_n > 300, "{soma:?}");
    assert!(soma.discarded_low_quality > 200 && soma.reads_out > 4000, "{soma:?}");

    // ---------------------------------------- .gz de saida, lido de volta, igual ao texto
    let o = TrimOptions { adapter: Some(ILLUMINA.into()), quality_cutoff: 20, trim_n: true, min_length: 20, ..TrimOptions::default() };
    let r1 = apara("variado.fastq", "aparado.fastq", &o);
    let r2 = apara("variado.fastq.gz", "aparado.fastq.gz", &o);
    let r3 = apara("variado_crlf.fastq", "aparado.fq.lz4", &o);
    assert!(r1 == r2 && r1 == r3 && r1 == ref_simples(&variado, &o).1);
    let plano = std::fs::read(p("aparado.fastq")).unwrap();
    assert!(plano == ref_simples(&variado, &o).0 && !plano.contains(&b'\r'));
    assert!(le(&p("aparado.fastq.gz")) == plano && le(&p("aparado.fq.lz4")) == plano);
    let bruto = std::fs::read(p("aparado.fastq.gz")).unwrap();
    assert!(bruto.starts_with(&[0x1f, 0x8b]) && bruto.len() < plano.len() * 3 / 4, "the .gz output must be compressed");
    assert!(r1.reads_with_adapter > 200 && r1.bases_removed_quality > 50_000 && r1.reads_out < 4000 && r1.discarded_too_short > 0, "{r1:?}");
    // o que o relatorio achou, o aparo tira: no arquivo aparado o adaptador some
    let depois = fastq_stats(&p("aparado.fastq.gz"), None).unwrap();
    assert!(depois.adapters[0].reads == 0 && depois.likely_adapter != Some(KNOWN_ADAPTERS[0]) && depois.min_length >= 20, "{}", depois.resumo());
    assert_eq!((depois.reads, depois.bases), (r1.reads_out, r1.bases_out));
    n += 2;

    // ---------------------------------------- entrada vazia
    grava("vazio.fastq", b"");
    let r = apara("vazio.fastq", "vazio_saida.fastq", &o);
    assert_eq!(r, TrimSummary::default());
    assert_eq!(std::fs::read(p("vazio_saida.fastq")).unwrap(), b"");
    let r = fastq_trim(&p("vazio.fastq"), &p("v1.fastq.gz"), Some(&p("so_linhas.fastq")), Some(&p("v2.fastq")), &o).unwrap();
    assert_eq!(r, TrimSummary { paired: true, ..TrimSummary::default() });
    assert!(le(&p("v1.fastq.gz")).is_empty() && le(&p("v2.fastq")).is_empty());
    n += 2;

    // ---------------------------------------- pareado, aleatorio, contra a referencia
    let mut soma = TrimSummary::default();
    for rodada in 0..30 {
        let r1 = gera_aparo(&mut s, 300, &ads_aparo, 33, "x");
        let mut r2 = gera_aparo(&mut s, 300, &ads_aparo, 33, "y");
        // os dois estilos de nome de par
        let (mut r1, estilo) = (r1, rodada % 3);
        for (i, (a, b)) in r1.iter_mut().zip(r2.iter_mut()).enumerate() {
            (a.nome, b.nome) = match estilo {
                0 => (format!("par{i}/1"), format!("par{i}/2")),
                1 => (format!("par{i} 1:N:0:ACGT"), format!("par{i} 2:N:0:ACGT")),
                _ => (format!("par{i}"), format!("par{i}\tsegunda")),
            };
        }
        let o = opcoes_aleatorias(&mut s, &ads_aparo, true);
        grava("p1.fastq", &texto(&r1, "\n"));
        grava("p2.fastq", &texto(&r2, if rodada % 4 == 0 { "\r\n" } else { "\n" }));
        let (s1, s2) = if rodada % 2 == 0 { ("p1_saida.fastq", "p2_saida.fastq.gz") } else { ("p1_saida.fq.lz4", "p2_saida.fastq") };
        let r = fastq_trim(&p("p1.fastq"), &p(s1), Some(&p("p2.fastq")), Some(&p(s2)), &o).unwrap();
        let (e1, e2, resumo) = ref_pareado(&r1, &r2, &o);
        assert_eq!(r, resumo, "round {rodada}: {o:?}");
        let (o1, o2) = (le(&p(s1)), le(&p(s2)));
        assert!(o1 == e1 && o2 == e2, "round {rodada}: the paired output differs from the reference; {o:?}");
        // as duas saidas tem o mesmo numero de registros, e os nomes continuam em par
        let (g1, g2) = (registros(&o1), registros(&o2));
        assert_eq!((g1.len() as u64, g2.len() as u64), (r.reads_out, r.reads_out));
        assert!(g1.iter().zip(&g2).all(|(a, b)| nome_par(&a.nome) == nome_par(&b.nome) && a.nome != b.nome));
        assert!(r.paired && r.reads_in == 300);
        assert_eq!(r.reads_in, r.reads_out + r.discarded_too_short + r.discarded_too_many_n + r.discarded_low_quality);
        soma.reads_with_adapter += r.reads_with_adapter;
        soma.reads_with_adapter2 += r.reads_with_adapter2;
        soma.reads_out += r.reads_out;
        soma.discarded_too_short += r.discarded_too_short;
        n += 1;
    }
    assert!(soma.reads_with_adapter > 1000 && soma.reads_with_adapter2 > 1000 && soma.reads_out > 2000 && soma.discarded_too_short > 500, "{soma:?}");

    // ---------------------------------------- pareado, a mao: o par sai quando QUALQUER leitura e reprovada
    let m1 = [reci("a/1", "ACGTACGT"), reci("b/1", "ACGTACGT"), reci("c/1", "ACG"), reci("d/1", "ACGTACGTAGATCGG"), reci("e/1", "ACNNACGT")];
    let m2 = [reci("a/2", "TTTTGGGG"), reci("b/2", "TT"), reci("c/2", "TTTTGGGG"), reci("d/2", "TTTTGGGGCTGTCTC"), reci("e/2", "TT")];
    grava("m1.fastq", &texto(&m1, "\n"));
    grava("m2.fastq", &texto(&m2, "\n"));
    let o = TrimOptions { adapter: Some(ILLUMINA.into()), adapter2: Some("CTGTCTCTTATACACATCT".into()), min_length: 5, max_n: Some(0), ..TrimOptions::default() };
    let r = fastq_trim(&p("m1.fastq"), &p("m1_saida.fastq"), Some(&p("m2.fastq")), Some(&p("m2_saida.fastq")), &o).unwrap();
    assert_eq!(le(&p("m1_saida.fastq")), b"@a/1\nACGTACGT\n+\nIIIIIIII\n@d/1\nACGTACGT\n+\nIIIIIIII\n");
    assert_eq!(le(&p("m2_saida.fastq")), b"@a/2\nTTTTGGGG\n+\nIIIIIIII\n@d/2\nTTTTGGGG\n+\nIIIIIIII\n");
    // b: so a leitura 2 e curta; c: so a leitura 1; e: a 1 tem N demais e a 2 e curta -> conta como curta
    assert_eq!(
        r,
        TrimSummary {
            paired: true,
            reads_in: 5,
            reads_out: 2,
            bases_in: 42 + 35,
            bases_out: 32,
            reads_with_adapter: 1,
            reads_with_adapter2: 1,
            bases_removed_adapter: 14,
            discarded_too_short: 3,
            ..TrimSummary::default()
        }
    );
    // sem adapter2 a leitura 2 usa `adapter`: o Nextera fica; so com adapter2, so a leitura 2 e aparada
    let r = fastq_trim(&p("m1.fastq"), &p("m1_b.fastq"), Some(&p("m2.fastq")), Some(&p("m2_b.fastq")), &TrimOptions { adapter2: None, ..o.clone() }).unwrap();
    assert_eq!((r.reads_with_adapter, r.reads_with_adapter2, r.bases_removed_adapter), (1, 0, 7));
    assert!(le(&p("m2_b.fastq")).ends_with(&texto(&m2[3..4], "\n")));
    let r = fastq_trim(&p("m1.fastq"), &p("m1_c.fastq"), Some(&p("m2.fastq")), Some(&p("m2_c.fastq")), &TrimOptions { adapter: None, ..o.clone() }).unwrap();
    assert_eq!((r.reads_with_adapter, r.reads_with_adapter2, r.bases_removed_adapter), (0, 1, 7));
    assert!(le(&p("m1_c.fastq")).ends_with(&texto(&m1[3..4], "\n")));
    n += 3;

    // ---------------------------------------- pareado: arquivos fora de passo ou truncados
    let par = |a: &str, b: &str| fastq_trim(&p(a), &p("e1.fastq"), Some(&p(b)), Some(&p("e2.fastq")), &TrimOptions::default());
    let mut fora = m2.to_vec();
    fora.swap(2, 3);
    grava("fora.fastq", &texto(&fora, "\n"));
    let e = par("m1.fastq", "fora.fastq").expect_err("out of step");
    assert!(e.contains("out of step at record 3") && e.contains("read 1 is `c/1`") && e.contains("read 2 is `d/2`"), "{e}");
    assert!(e.contains("line 9") && e.contains(&p("m1.fastq")) && e.contains("incomplete: do not use it"), "{e}");
    // os dois registros que estavam em passo foram gravados antes do erro
    assert_eq!(registros(&std::fs::read(p("e1.fastq")).unwrap()).len(), 2);
    grava("curto.fastq", &texto(&m2[..3], "\n"));
    let e = par("m1.fastq", "curto.fastq").expect_err("truncated pair");
    assert!(e.contains("out of step at record 4") && e.contains("`d/1`") && e.contains(&p("curto.fastq")) && e.contains("ends after 3 records"), "{e}");
    assert!(e.contains("read 2: none"), "{e}");
    grava("curto1.fastq", &texto(&m1[..1], "\n"));
    let e = par("curto1.fastq", "m2.fastq").expect_err("truncated pair");
    assert!(e.contains("out of step at record 2") && e.contains("`b/2`") && e.contains(&p("curto1.fastq")) && e.contains("ends after 1 records"), "{e}");
    assert!(e.contains("read 1: none"), "{e}");
    // nomes diferentes so depois do espaco, ou so no /1 e /2, estao em passo
    grava("n1.fastq", b"@x/1 primeira\nAC\n+\nII\n@y\tz\nAC\n+\nII\n");
    grava("n2.fastq", b"@x/2 segunda\nAC\n+\nII\n@y w\nAC\n+\nII\n");
    assert_eq!(par("n1.fastq", "n2.fastq").unwrap().reads_out, 2);
    grava("n3.fastq", b"@x/2 segunda\nAC\n+\nII\n@yy w\nAC\n+\nII\n");
    let e = par("n1.fastq", "n3.fastq").expect_err("out of step");
    assert!(e.contains("record 2") && e.contains("`y`") && e.contains("`yy`"), "{e}");
    n += 5;

    // ---------------------------------------- chamadas invalidas
    let chamada = |e2: Option<&str>, s2: Option<&str>, o: &TrimOptions| {
        fastq_trim(&p("m1.fastq"), &p("c1.fastq"), e2.map(&p).as_deref(), s2.map(&p).as_deref(), o).expect_err("an error was expected")
    };
    assert!(chamada(Some("m2.fastq"), None, &padrao).contains("without a second output file"));
    assert!(chamada(None, Some("c2.fastq"), &padrao).contains("without a second input file"));
    assert!(chamada(None, None, &TrimOptions { adapter2: Some("ACGT".into()), ..TrimOptions::default() }).contains("adapter2 was given without a second input file"));
    assert!(chamada(None, None, &TrimOptions { adapter: Some("auto".into()), ..TrimOptions::default() }).contains("invalid character `U`"));
    assert!(chamada(Some("m2.fastq"), Some("c2.fastq"), &TrimOptions { adapter: Some(ILLUMINA.into()), adapter2: Some("AC!".into()), ..TrimOptions::default() }).contains("invalid character `!`"));
    assert!(chamada(None, None, &TrimOptions { adapter: Some(ILLUMINA.into()), max_error_rate: 1.5, ..TrimOptions::default() }).contains("max_error_rate is 1.5"));
    assert!(chamada(None, None, &TrimOptions { adapter: Some(ILLUMINA.into()), min_overlap: 0, ..TrimOptions::default() }).contains("min_overlap is 0"));
    assert!(chamada(None, None, &TrimOptions { phred_offset: 50, ..TrimOptions::default() }).contains("phred_offset is 50"));
    // nenhuma dessas chamadas chegou a criar a saida
    assert!(!std::path::Path::new(&p("c1.fastq")).exists() && !std::path::Path::new(&p("c2.fastq")).exists());
    assert!(fastq_trim(&p("nao_existe.fastq"), &p("c1.fastq"), None, None, &padrao).expect_err("missing input").contains("nao_existe.fastq"));
    assert!(!std::path::Path::new(&p("c1.fastq")).exists());
    n += 1;
    // a saida nao pode ser um arquivo de entrada, nem as duas saidas o mesmo arquivo
    let antes = std::fs::read(p("m1.fastq")).unwrap();
    let e = fastq_trim(&p("m1.fastq"), &p("m1.fastq"), None, None, &padrao).expect_err("same file");
    assert!(e.contains("is the input file"), "{e}");
    let e = fastq_trim(&p("m1.fastq"), &p("c1.fastq"), Some(&p("m2.fastq")), Some(&p("m1.fastq")), &padrao).expect_err("same file");
    assert!(e.contains("is the input file"), "{e}");
    let e = fastq_trim(&p("m1.fastq"), &p("c1.fastq"), Some(&p("m2.fastq")), Some(&p("c1.fastq")), &padrao).expect_err("same file");
    assert!(e.contains("two outputs are the same file"), "{e}");
    assert_eq!(std::fs::read(p("m1.fastq")).unwrap(), antes, "the input must be intact");
    assert!(!std::path::Path::new(&p("c1.fastq")).exists());
    n += 1;

    // ---------------------------------------- FASTQ invalido no aparo: o mesmo erro, com a linha
    let e = fastq_trim(&p("tamanhos.fastq"), &p("inv.fastq"), None, None, &padrao).expect_err("malformed");
    assert!(e.contains("line 8:") && e.contains("4 bases but 3 quality values") && e.contains("incomplete: do not use it"), "{e}");
    assert_eq!(std::fs::read(p("inv.fastq")).unwrap(), b"@r1\nACGT\n+\nIIII\n", "what came before the defect was written");
    let e = fastq_trim(&p("corta2.fastq"), &p("inv.fastq.gz"), None, None, &padrao).expect_err("truncated");
    assert!(e.contains("line 6,") && e.contains("the `+` separator line"), "{e}");
    grava("bons.fastq", b"@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n");
    let e = fastq_trim(&p("bons.fastq"), &p("inv1.fastq"), Some(&p("sem_mais.fastq")), Some(&p("inv2.fastq")), &padrao).expect_err("malformed");
    assert!(e.contains(&p("sem_mais.fastq")) && e.contains("line 7:") && e.contains("must start with `+`"), "{e}");
    n += 3;

    // ================================================================
    //  PARTE 3: diferencial contra o cutadapt, se CUTADAPT estiver definida
    // ================================================================
    if let Some(caminho) = std::env::var("CUTADAPT").ok().filter(|c| !c.trim().is_empty()) {
        n += diferencial(&dir, caminho.trim(), &mut s);
    }

    println!("{n} FASTQ QC tests passed");

    // ---------------------------------------- arquivos reais, se informados
    for arq in std::env::args().skip(1) {
        let t = std::time::Instant::now();
        let st = fastq_stats(&arq, None).unwrap();
        print!("\n[{:.1} s] {}", t.elapsed().as_secs_f64(), st.resumo());
        json(&st.json());
        let o = TrimOptions {
            adapter: st.likely_adapter.map(|a| a.1.to_string()),
            quality_cutoff: 20,
            trim_n: true,
            min_length: 20,
            phred_offset: st.phred_offset,
            ..TrimOptions::default()
        };
        let t = std::time::Instant::now();
        let r = fastq_trim(&arq, &p("real_aparado.fastq.gz"), None, None, &o).unwrap();
        println!("[{:.1} s] fastq_trim(quality_cutoff=20, trim_n, min_length=20, adapter={:?}): {}", t.elapsed().as_secs_f64(), o.adapter, r.json());
        let depois = fastq_stats(&p("real_aparado.fastq.gz"), None).unwrap();
        assert_eq!((depois.reads, depois.bases), (r.reads_out, r.bases_out));
    }
    let _ = std::fs::remove_dir_all(&dir);
}
