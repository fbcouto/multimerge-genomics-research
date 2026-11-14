// =====================================================================
//  CONTROLE DE QUALIDADE DE FASTQ: RELATORIO, APARO E FILTRO
//
//  Duas ferramentas em fluxo para leituras brutas de sequenciador, de
//  qualquer tamanho:
//
//   fastq_stats   uma passada que diz o que ha no arquivo -- comprimentos,
//                 qualidade por ciclo, composicao de bases, GC, Ns,
//                 duplicacao, adaptadores -- e o que fazer a respeito.
//   fastq_trim    aparo por qualidade, remocao de adaptador 3', aparo de N e
//                 filtros, em leitura simples ou pareada. O resultado e o do
//                 cutadapt (conferido contra o cutadapt 5.2, leitura por
//                 leitura, em trials/src/bin/fastq_qc.rs).
//
//  A memoria nao depende do tamanho do arquivo: o relatorio guarda um
//  histograma por ciclo (so os MAX_POS primeiros), um por comprimento e os
//  hashes das DUP_SAMPLE primeiras leituras; o aparo guarda um registro (ou
//  um par) de cada vez.
//
//  O que nao e FASTQ valido -- cabecalho sem `@`, separador sem `+`,
//  sequencia e qualidade de tamanhos diferentes, arquivo que acaba no meio
//  de um registro -- e ERRO, com a linha e o que se esperava: seguir adiante
//  daria numeros errados em silencio.
// =====================================================================

use crate::externo::{abre_leitura, Escritor};
use crate::nomes::json_texto;
use std::collections::{BTreeMap, HashSet};
use std::fmt::Display;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

/// Quantos ciclos (posicoes na leitura) entram nas tabelas por posicao.
/// Leituras mais longas contam inteiras nos totais.
pub const MAX_POS: usize = 1000;
/// Quantas leituras, do comeco do arquivo, entram na medida de duplicacao.
pub const DUP_SAMPLE: u64 = 1_000_000;
/// Quantas bases do comeco de cada adaptador conhecido sao procuradas.
pub const ADAPTER_PROBE: usize = 12;
/// Adaptadores conhecidos: (nome, sequencia). A sequencia inteira e o que se
/// passa a `fastq_trim` em `adapter`.
pub const KNOWN_ADAPTERS: [(&str, &str); 4] = [
    ("Illumina Universal", "AGATCGGAAGAGC"),
    ("Nextera Transposase", "CTGTCTCTTATACACATCT"),
    ("Illumina Small RNA 3'", "TGGAATTCTCGGGTGCCAAGG"),
    ("BGI/MGI", "AAGTCGGAGGCCAAGCGGTC"),
];

/// Caracteres de qualidade validos: ASCII 33 a 126.
const FAIXAS_Q: usize = 94;
const MASCARA_SONDA: u32 = (1 << (2 * ADAPTER_PROBE)) - 1;

/// A, C, G, T (maiuscula ou minuscula) -> 0..3; qualquer outra coisa -> 4 (N).
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

// ---------------------------------------------------------------------
//  LEITURA DE REGISTROS
// ---------------------------------------------------------------------

/// Um trecho do arquivo, para mostrar numa mensagem.
fn mostra(b: &[u8], max: usize) -> String {
    let s = String::from_utf8_lossy(b);
    if s.chars().count() > max {
        format!("{}...", s.chars().take(max).collect::<String>())
    } else {
        s.into_owned()
    }
}

/// O nome da leitura: o cabecalho sem o `@`, ate o primeiro espaco.
fn nome_do_registro(cab: &[u8]) -> &[u8] {
    let s = cab.get(1..).unwrap_or(&[]);
    let fim = s.iter().position(|c| c.is_ascii_whitespace()).unwrap_or(s.len());
    &s[..fim]
}

/// O nome que as duas leituras de um par tem em comum: sem `/1` ou `/2`.
fn nome_do_par(cab: &[u8]) -> &[u8] {
    let n = nome_do_registro(cab);
    if n.ends_with(b"/1") || n.ends_with(b"/2") {
        &n[..n.len() - 2]
    } else {
        n
    }
}

/// Le uma linha, sem o fim de linha (LF ou CRLF). `false` no fim do arquivo.
fn le_linha(r: &mut dyn BufRead, buf: &mut Vec<u8>, caminho: &str, linha: &mut u64) -> Result<bool, String> {
    buf.clear();
    let n = r.read_until(b'\n', buf).map_err(|e| {
        format!("Error reading {caminho} after line {linha}: {e}. The file may be corrupted or incomplete (an interrupted download leaves this behind); fetch it again")
    })?;
    if n == 0 {
        return Ok(false);
    }
    *linha += 1;
    if buf.last() == Some(&b'\n') {
        buf.pop();
    }
    if buf.last() == Some(&b'\r') {
        buf.pop();
    }
    Ok(true)
}

/// Leitor de FASTQ de 4 linhas por registro, com os buffers reaproveitados.
struct Leitor {
    r: Box<dyn BufRead>,
    caminho: String,
    /// Linhas lidas ate aqui.
    linha: u64,
    /// Registros completos lidos ate aqui.
    registros: u64,
    cab: Vec<u8>,
    seq: Vec<u8>,
    sep: Vec<u8>,
    qual: Vec<u8>,
}

impl Leitor {
    fn abre(caminho: &str) -> Result<Leitor, String> {
        Ok(Leitor {
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

    fn nome(&self) -> String {
        mostra(nome_do_registro(&self.cab), 60)
    }

    fn truncado(&self, k: u64, faltou: &str) -> String {
        format!(
            "{}: the file ends at line {}, in the middle of record {k} (`{}`): expected {faltou} next. An interrupted download or copy usually leaves this behind; fetch the file again",
            self.caminho,
            self.linha,
            self.nome()
        )
    }

    /// Le o proximo registro para `cab`, `seq` e `qual`. `false` no fim do
    /// arquivo. Linhas vazias ENTRE registros sao puladas; dentro de um
    /// registro, uma linha vazia e conteudo (leitura de tamanho zero).
    fn proximo(&mut self) -> Result<bool, String> {
        loop {
            if !le_linha(&mut *self.r, &mut self.cab, &self.caminho, &mut self.linha)? {
                return Ok(false);
            }
            if !self.cab.is_empty() {
                break;
            }
        }
        let k = self.registros + 1;
        if self.cab[0] != b'@' {
            return Err(format!(
                "{}, line {}: expected the header of record {k}, which must start with `@`, but the line is `{}`. A FASTQ record has exactly 4 lines (@name, sequence, +, qualities); validate_file(...) lists everything that is wrong with the file",
                self.caminho,
                self.linha,
                mostra(&self.cab, 30)
            ));
        }
        if !le_linha(&mut *self.r, &mut self.seq, &self.caminho, &mut self.linha)? {
            return Err(self.truncado(k, "the sequence line"));
        }
        if !le_linha(&mut *self.r, &mut self.sep, &self.caminho, &mut self.linha)? {
            return Err(self.truncado(k, "the `+` separator line"));
        }
        if self.sep.first() != Some(&b'+') {
            return Err(format!(
                "{}, line {}: expected the separator line of record {k} (`{}`), which must start with `+`, but the line is `{}`. Every record must have exactly 4 lines: a sequence wrapped over several lines has to be joined into one first",
                self.caminho,
                self.linha,
                self.nome(),
                mostra(&self.sep, 30)
            ));
        }
        // Uma leitura de tamanho zero no fim de um arquivo sem a ultima linha
        // vazia: a linha de qualidade que falta e a propria linha vazia.
        if !le_linha(&mut *self.r, &mut self.qual, &self.caminho, &mut self.linha)? && !self.seq.is_empty() {
            return Err(self.truncado(k, "the quality line"));
        }
        if self.qual.len() != self.seq.len() {
            return Err(format!(
                "{}, line {}: record {k} (`{}`) has {} bases but {} quality values; expected one quality character per base. The file was probably cut or edited: fetch it again",
                self.caminho,
                self.linha,
                self.nome(),
                self.seq.len(),
                self.qual.len()
            ));
        }
        self.registros = k;
        Ok(true)
    }

    /// Ainda ha alguma linha nao vazia?
    fn resta_conteudo(&mut self) -> Result<bool, String> {
        loop {
            if !le_linha(&mut *self.r, &mut self.cab, &self.caminho, &mut self.linha)? {
                return Ok(false);
            }
            if !self.cab.is_empty() {
                return Ok(true);
            }
        }
    }
}

// ---------------------------------------------------------------------
//  PARTE 1: RELATORIO DE QUALIDADE
// ---------------------------------------------------------------------

/// Um ciclo (posicao na leitura, a partir de 1).
#[derive(Clone, Debug, PartialEq)]
pub struct PositionStats {
    /// Bases vistas nesta posicao (leituras com pelo menos este comprimento).
    pub bases: u64,
    pub mean_quality: f64,
    /// Percentis da qualidade, pelo posto mais proximo: o menor valor com
    /// pelo menos p% das bases nele ou abaixo dele.
    pub q10: u8,
    pub q25: u8,
    pub q50: u8,
    pub q75: u8,
    pub q90: u8,
    pub a: u64,
    pub c: u64,
    pub g: u64,
    pub t: u64,
    /// N e qualquer outra coisa que nao seja A, C, G ou T.
    pub n: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AdapterHit {
    pub name: &'static str,
    pub sequence: &'static str,
    /// Leituras que contem as ADAPTER_PROBE primeiras bases, sem erro.
    pub reads: u64,
    pub percent: f64,
}

#[derive(Clone, Debug)]
pub struct FastqStats {
    pub path: String,
    pub reads: u64,
    pub bases: u64,
    /// `false` quando `max_reads` parou a passada antes do fim do arquivo.
    pub complete: bool,
    pub min_length: usize,
    pub max_length: usize,
    pub mean_length: f64,
    /// (comprimento, leituras), em ordem de comprimento.
    pub length_histogram: Vec<(usize, u64)>,
    /// 33 ou 64. Todas as qualidades abaixo ja estao nesta escala.
    pub phred_offset: u8,
    /// (menor, maior) caractere de qualidade visto, em codigo ASCII: a
    /// evidencia de que a escala foi tirada. `None` sem nenhuma base.
    pub quality_char_range: Option<(u8, u8)>,
    /// Qualidade media, ponderada pelas bases.
    pub mean_quality: f64,
    /// 94 contagens: leituras pela parte inteira da sua qualidade media.
    /// Leituras de tamanho zero nao entram.
    pub read_quality_histogram: Vec<u64>,
    /// min(max_length, MAX_POS): quantos ciclos ha em `positions`.
    pub positions_reported: usize,
    pub positions: Vec<PositionStats>,
    /// (G+C) sobre (A+C+G+T), em porcentagem.
    pub gc_percent: f64,
    /// 101 contagens: leituras pelo GC arredondado. Leituras sem nenhuma
    /// base A/C/G/T nao entram.
    pub gc_histogram: Vec<u64>,
    pub n_bases: u64,
    pub n_percent: f64,
    pub reads_with_n: u64,
    /// Das `duplication_sample` primeiras leituras, a porcentagem cuja
    /// sequencia exata ja tinha aparecido.
    pub duplication_percent: f64,
    pub duplication_sample: u64,
    pub adapters: Vec<AdapterHit>,
    /// (nome, sequencia) do adaptador mais frequente, se em 0,5% ou mais.
    pub likely_adapter: Option<(&'static str, &'static str)>,
}

/// Hash de 64 bits de uma sequencia: FNV-1a seguido do finalizador do
/// MurmurHash3 (fmix64), que espalha os bits que o FNV deixa fracos.
fn hash64(b: &[u8]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &c in b {
        h ^= u64::from(c);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51_afd7_ed55_8ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    h ^ (h >> 33)
}

/// O conjunto guarda hashes ja espalhados: nao ha por que espalhar de novo.
#[derive(Default)]
struct Identidade(u64);

impl Hasher for Identidade {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0 << 8) | u64::from(b);
        }
    }
    fn write_u64(&mut self, v: u64) {
        self.0 = v;
    }
}

fn razao(a: u64, b: u64) -> f64 {
    if b == 0 {
        0.0
    } else {
        a as f64 / b as f64
    }
}

/// Percentil pelo posto mais proximo, a partir de um histograma.
fn percentil(h: &[u64; FAIXAS_Q], n: u64, p: u64) -> usize {
    let posto = ((u128::from(p) * u128::from(n) + 99) / 100).max(1);
    let mut acumulado = 0u128;
    for (v, &c) in h.iter().enumerate() {
        acumulado += u128::from(c);
        if acumulado >= posto {
            return v;
        }
    }
    0
}

/// Tudo e acumulado sobre o caractere BRUTO de qualidade: a escala (33 ou 64)
/// so se conhece no fim, pelo menor caractere visto, e e descontada la.
struct Acumulador {
    reads: u64,
    bases: u64,
    min_len: usize,
    max_len: usize,
    comprimentos: BTreeMap<usize, u64>,
    min_char: u8,
    max_char: u8,
    soma_bruta: u64,
    pos_q: Vec<[u64; FAIXAS_Q]>,
    pos_b: Vec<[u64; 5]>,
    media_leitura: [u64; FAIXAS_Q],
    total: [u64; 5],
    gc: [u64; 101],
    com_n: u64,
    vistos: HashSet<u64, BuildHasherDefault<Identidade>>,
    dup_examinadas: u64,
    dup_repetidas: u64,
    sondas: [u32; 4],
    com_adaptador: [u64; 4],
}

impl Acumulador {
    fn novo() -> Acumulador {
        let mut sondas = [0u32; 4];
        for (s, (_, seq)) in sondas.iter_mut().zip(KNOWN_ADAPTERS) {
            for &b in seq.as_bytes().iter().take(ADAPTER_PROBE) {
                *s = (*s << 2) | u32::from(BASE[b as usize] & 3);
            }
        }
        Acumulador {
            reads: 0,
            bases: 0,
            min_len: 0,
            max_len: 0,
            comprimentos: BTreeMap::new(),
            min_char: u8::MAX,
            max_char: 0,
            soma_bruta: 0,
            pos_q: Vec::new(),
            pos_b: Vec::new(),
            media_leitura: [0; FAIXAS_Q],
            total: [0; 5],
            gc: [0; 101],
            com_n: 0,
            vistos: HashSet::default(),
            dup_examinadas: 0,
            dup_repetidas: 0,
            sondas,
            com_adaptador: [0; 4],
        }
    }

    /// Soma uma leitura. `Err((coluna, byte))`: caractere de qualidade fora
    /// de ASCII 33 a 126.
    fn soma(&mut self, seq: &[u8], qual: &[u8]) -> Result<(), (usize, u8)> {
        let n = seq.len();
        let mut soma_q = 0u64;
        let (mut min_c, mut max_c) = (self.min_char, self.max_char);
        for (i, &c) in qual.iter().enumerate() {
            if !(33..=126).contains(&c) {
                return Err((i + 1, c));
            }
            soma_q += u64::from(c);
            min_c = min_c.min(c);
            max_c = max_c.max(c);
        }
        self.min_char = min_c;
        self.max_char = max_c;
        self.soma_bruta += soma_q;
        self.reads += 1;
        self.bases += n as u64;
        *self.comprimentos.entry(n).or_insert(0) += 1;
        if self.reads == 1 || n < self.min_len {
            self.min_len = n;
        }
        self.max_len = self.max_len.max(n);

        let lim = n.min(MAX_POS);
        if self.pos_q.len() < lim {
            self.pos_q.resize(lim, [0; FAIXAS_Q]);
            self.pos_b.resize(lim, [0; 5]);
        }
        let mut comp = [0u64; 5];
        let mut achou = [false; 4];
        // As 12 ultimas bases, 2 bits cada: acha os adaptadores numa passada.
        let (mut codigo, mut validas) = (0u32, 0usize);
        for (i, (&s, &q)) in seq.iter().zip(qual).enumerate() {
            let b = BASE[s as usize] as usize;
            comp[b] += 1;
            if i < lim {
                self.pos_b[i][b] += 1;
                self.pos_q[i][(q - 33) as usize] += 1;
            }
            if b == 4 {
                validas = 0;
            } else {
                codigo = ((codigo << 2) | b as u32) & MASCARA_SONDA;
                validas += 1;
                if validas >= ADAPTER_PROBE {
                    for (a, sonda) in achou.iter_mut().zip(&self.sondas) {
                        *a |= codigo == *sonda;
                    }
                }
            }
        }
        for (t, c) in self.total.iter_mut().zip(&comp) {
            *t += c;
        }
        for (t, a) in self.com_adaptador.iter_mut().zip(&achou) {
            *t += u64::from(*a);
        }
        if n > 0 {
            // floor(media) nao depende da escala: desconta-se um inteiro.
            self.media_leitura[(soma_q / n as u64) as usize - 33] += 1;
        }
        let acgt = comp[0] + comp[1] + comp[2] + comp[3];
        if acgt > 0 {
            let gc = comp[1] + comp[2];
            self.gc[((200 * gc + acgt) / (2 * acgt)) as usize] += 1;
        }
        if comp[4] > 0 {
            self.com_n += 1;
        }
        if self.dup_examinadas < DUP_SAMPLE {
            self.dup_examinadas += 1;
            if !self.vistos.insert(hash64(seq)) {
                self.dup_repetidas += 1;
            }
        }
        Ok(())
    }

    fn fecha(self, caminho: &str, complete: bool) -> FastqStats {
        // Phred+64 so quando NENHUM caractere fica abaixo de `@` (64) E algum
        // passa de `K` (75). Sem a segunda condicao, um arquivo em que toda
        // base tem qualidade 31 ou mais em Phred+33 -- leitura simulada com
        // `I` constante, por exemplo -- sairia como Phred+64 de qualidade
        // pessima: nenhum dado real de Phred+64 tem todas as bases entre 0 e 11.
        let offset: u8 = if self.min_char != u8::MAX && self.min_char >= 64 && self.max_char > b'K' { 64 } else { 33 };
        let desloc = (offset - 33) as usize;

        let positions: Vec<PositionStats> = self
            .pos_q
            .iter()
            .zip(&self.pos_b)
            .map(|(h, b)| {
                let n: u64 = h.iter().sum();
                let soma: u64 = h.iter().enumerate().map(|(v, &c)| v.saturating_sub(desloc) as u64 * c).sum();
                let p = |pc: u64| percentil(h, n, pc).saturating_sub(desloc) as u8;
                PositionStats {
                    bases: n,
                    mean_quality: razao(soma, n),
                    q10: p(10),
                    q25: p(25),
                    q50: p(50),
                    q75: p(75),
                    q90: p(90),
                    a: b[0],
                    c: b[1],
                    g: b[2],
                    t: b[3],
                    n: b[4],
                }
            })
            .collect();

        let mut read_quality_histogram = vec![0u64; FAIXAS_Q];
        for v in desloc..FAIXAS_Q {
            read_quality_histogram[v - desloc] = self.media_leitura[v];
        }

        let adapters: Vec<AdapterHit> = KNOWN_ADAPTERS
            .iter()
            .zip(&self.com_adaptador)
            .map(|(&(name, sequence), &reads)| AdapterHit { name, sequence, reads, percent: 100.0 * razao(reads, self.reads) })
            .collect();
        let mut provavel: Option<&AdapterHit> = None;
        for a in &adapters {
            if !provavel.is_some_and(|p| a.percent <= p.percent) {
                provavel = Some(a);
            }
        }
        let likely_adapter = provavel.filter(|a| a.percent >= 0.5).map(|a| (a.name, a.sequence));

        let acgt = self.total[0] + self.total[1] + self.total[2] + self.total[3];
        FastqStats {
            path: caminho.to_string(),
            reads: self.reads,
            bases: self.bases,
            complete,
            min_length: self.min_len,
            max_length: self.max_len,
            mean_length: razao(self.bases, self.reads),
            length_histogram: self.comprimentos.into_iter().collect(),
            phred_offset: offset,
            quality_char_range: if self.min_char == u8::MAX { None } else { Some((self.min_char, self.max_char)) },
            mean_quality: if self.bases == 0 { 0.0 } else { razao(self.soma_bruta, self.bases) - f64::from(offset) },
            read_quality_histogram,
            positions_reported: positions.len(),
            positions,
            gc_percent: 100.0 * razao(self.total[1] + self.total[2], acgt),
            gc_histogram: self.gc.to_vec(),
            n_bases: self.total[4],
            n_percent: 100.0 * razao(self.total[4], self.bases),
            reads_with_n: self.com_n,
            duplication_percent: 100.0 * razao(self.dup_repetidas, self.dup_examinadas),
            duplication_sample: self.dup_examinadas,
            adapters,
            likely_adapter,
        }
    }
}

/// Uma passada em fluxo sobre um FASTQ (texto, `.gz` ou `.lz4`), parando
/// depois de `max_reads` registros, se dado.
///
/// So devolve `Err` quando o arquivo nao abre, nao le ou nao e FASTQ valido;
/// um arquivo vazio tem zero leituras e nao e erro.
pub fn fastq_stats(caminho: &str, max_reads: Option<u64>) -> Result<FastqStats, String> {
    let mut leitor = Leitor::abre(caminho)?;
    let mut a = Acumulador::novo();
    let mut complete = true;
    loop {
        if max_reads.is_some_and(|m| a.reads >= m) {
            complete = !leitor.resta_conteudo()?;
            break;
        }
        if !leitor.proximo()? {
            break;
        }
        if let Err((coluna, c)) = a.soma(&leitor.seq, &leitor.qual) {
            return Err(format!(
                "{caminho}, line {}: the quality line of record {} (`{}`) has a character outside the valid range at column {coluna} (ASCII code {c}; qualities go from `!`, code 33, to `~`, code 126). The file is damaged or is not FASTQ: fetch it again or check it with validate_file(...)",
                leitor.linha,
                leitor.registros,
                leitor.nome()
            ));
        }
    }
    Ok(a.fecha(caminho, complete))
}

fn junta<T: Display>(itens: impl Iterator<Item = T>) -> String {
    let mut o = String::new();
    for (i, x) in itens.enumerate() {
        if i > 0 {
            o.push(',');
        }
        o.push_str(&x.to_string());
    }
    o
}

/// Um numero de ponto flutuante como numero JSON.
fn num(x: f64) -> String {
    if x.is_finite() {
        format!("{x}")
    } else {
        "0".to_string()
    }
}

impl FastqStats {
    /// A cauda de baixa qualidade: (primeiro ciclo, ultimo ciclo) do trecho
    /// final em que a mediana fica abaixo de 20. So contam os ciclos a que
    /// pelo menos 5% das leituras chegam: uma unica leitura comprida nao
    /// define a cauda.
    fn cauda_ruim(&self) -> Option<(usize, usize)> {
        let ultimo = self.positions.iter().rposition(|p| p.bases.saturating_mul(20) >= self.reads)?;
        let mut primeiro = ultimo + 1;
        while primeiro > 0 && self.positions[primeiro - 1].q50 < 20 {
            primeiro -= 1;
        }
        if primeiro > ultimo {
            None
        } else {
            Some((primeiro + 1, ultimo + 1))
        }
    }

    /// O que merece atencao, em linguagem de gente, cada item com o conserto.
    /// A escala vem primeiro: os outros avisos dependem dela.
    pub fn warnings(&self) -> Vec<String> {
        let mut w = Vec::new();
        if let (64, Some((menor, maior))) = (self.phred_offset, self.quality_char_range) {
            w.push(format!(
                "The qualities are in the old Phred+64 scale (Illumina 1.3 to 1.7): the quality characters go from `{}` to `{}`, none below `@`. The numbers in this report already account for it, but current tools assume Phred+33 and would read inflated qualities: pass phred_offset=64 to fastq_trim(...).",
                char::from(menor),
                char::from(maior)
            ));
        }
        let escala = if self.phred_offset == 64 { ", phred_offset=64" } else { "" };
        if let Some((nome, seq)) = self.likely_adapter {
            let pct = self.adapters.iter().find(|a| a.name == nome).map(|a| a.percent).unwrap_or(0.0);
            w.push(format!(
                "Adapter found: {nome} in {pct:.1}% of the reads. Adapter bases are not part of the sample and disturb alignment and assembly: trim it with fastq_trim(adapter=\"{seq}\")."
            ));
        }
        if let Some((primeiro, ultimo)) = self.cauda_ruim() {
            if primeiro == 1 {
                w.push(format!(
                    "Low quality throughout: the median quality is below 20 at every cycle up to {ultimo}. Keep only what is usable with fastq_trim(quality_cutoff=20, min_length=30{escala}), and consider repeating the run."
                ));
            } else {
                w.push(format!(
                    "Low quality tail: from cycle {primeiro} to cycle {ultimo} the median quality is below 20 (it is {} at cycle {ultimo}). Trim it with fastq_trim(quality_cutoff=20{escala}).",
                    self.positions[ultimo - 1].q50
                ));
            }
        }
        let mut pior: Option<(usize, f64)> = None;
        for (i, p) in self.positions.iter().enumerate() {
            if p.bases > 0 && p.bases.saturating_mul(20) >= self.reads {
                let pct = 100.0 * razao(p.n, p.bases);
                if !pior.is_some_and(|x| pct <= x.1) {
                    pior = Some((i + 1, pct));
                }
            }
        }
        if self.n_percent >= 1.0 || pior.is_some_and(|x| x.1 >= 5.0) {
            let onde = pior.map(|(c, p)| format!("; the worst cycle is {c}, with {p:.1}% N")).unwrap_or_default();
            w.push(format!(
                "Many undetermined bases: N is {:.2}% of all bases and {:.1}% of the reads have at least one{onde}. Cut the Ns at the ends with fastq_trim(trim_n=True) and drop the reads that still have them with fastq_trim(max_n=0).",
                self.n_percent,
                100.0 * razao(self.reads_with_n, self.reads)
            ));
        }
        w
    }

    /// O relatorio em texto, pronto para mostrar.
    pub fn resumo(&self) -> String {
        let pct = |a: u64, b: u64| 100.0 * razao(a, b);
        let parcial = if self.complete { String::new() } else { format!(" (only the first {} reads were examined)", self.reads) };
        let mut o = format!("{}: FASTQ with {} reads and {} bases{parcial}\n", self.path, self.reads, self.bases);
        o.push_str(&format!(
            "  read length       min {}, mean {:.1}, max {}\n",
            self.min_length, self.mean_length, self.max_length
        ));
        o.push_str(&format!("  quality scale     Phred+{}\n", self.phred_offset));
        let com_media: u64 = self.read_quality_histogram.iter().sum();
        let de30: u64 = self.read_quality_histogram.iter().skip(30).sum();
        let ate20: u64 = self.read_quality_histogram.iter().take(20).sum();
        o.push_str(&format!(
            "  mean quality      {:.1} ({:.1}% of the reads average 30 or more, {:.1}% below 20)\n",
            self.mean_quality,
            pct(de30, com_media),
            pct(ate20, com_media)
        ));
        if let Some(fim) = self.positions.len().checked_sub(1) {
            let meio = fim / 2;
            let resto = if self.max_length > MAX_POS { format!(" (only the first {MAX_POS} cycles are reported)") } else { String::new() };
            o.push_str(&format!(
                "  quality by cycle  median {} at cycle 1, {} at cycle {}, {} at cycle {}{resto}\n",
                self.positions[0].q50,
                self.positions[meio].q50,
                meio + 1,
                self.positions[fim].q50,
                fim + 1
            ));
        }
        o.push_str(&format!("  GC content        {:.1}%\n", self.gc_percent));
        o.push_str(&format!(
            "  undetermined (N)  {} bases ({:.2}%), in {} reads ({:.1}%)\n",
            self.n_bases,
            self.n_percent,
            self.reads_with_n,
            pct(self.reads_with_n, self.reads)
        ));
        o.push_str(&format!(
            "  duplication       {:.1}% of the first {} reads repeat an earlier sequence\n",
            self.duplication_percent, self.duplication_sample
        ));
        let ads: Vec<String> = self.adapters.iter().map(|a| format!("{} {:.2}%", a.name, a.percent)).collect();
        o.push_str(&format!("  adapters          {}\n", ads.join(", ")));
        let avisos = self.warnings();
        if avisos.is_empty() {
            o.push_str("  no problems found\n");
        } else {
            o.push_str("  warnings:\n");
            for a in &avisos {
                o.push_str(&format!("    - {a}\n"));
            }
        }
        o
    }

    pub fn json(&self) -> String {
        let p = &self.positions;
        let por_posicao = format!(
            "{{\"bases\":[{}],\"mean_quality\":[{}],\"q10\":[{}],\"q25\":[{}],\"q50\":[{}],\"q75\":[{}],\"q90\":[{}],\"a\":[{}],\"c\":[{}],\"g\":[{}],\"t\":[{}],\"n\":[{}]}}",
            junta(p.iter().map(|x| x.bases)),
            junta(p.iter().map(|x| num(x.mean_quality))),
            junta(p.iter().map(|x| x.q10)),
            junta(p.iter().map(|x| x.q25)),
            junta(p.iter().map(|x| x.q50)),
            junta(p.iter().map(|x| x.q75)),
            junta(p.iter().map(|x| x.q90)),
            junta(p.iter().map(|x| x.a)),
            junta(p.iter().map(|x| x.c)),
            junta(p.iter().map(|x| x.g)),
            junta(p.iter().map(|x| x.t)),
            junta(p.iter().map(|x| x.n))
        );
        let adaptadores = junta(self.adapters.iter().map(|a| {
            format!(
                "{{\"name\":{},\"sequence\":{},\"reads\":{},\"percent\":{}}}",
                json_texto(a.name),
                json_texto(a.sequence),
                a.reads,
                num(a.percent)
            )
        }));
        let provavel = match self.likely_adapter {
            Some((n, s)) => format!("{{\"name\":{},\"sequence\":{}}}", json_texto(n), json_texto(s)),
            None => "null".to_string(),
        };
        format!(
            "{{\"path\":{},\"reads\":{},\"bases\":{},\"complete\":{},\"min_length\":{},\"max_length\":{},\"mean_length\":{},\"length_histogram\":[{}],\"phred_offset\":{},\"quality_char_range\":{},\"mean_quality\":{},\"read_quality_histogram\":[{}],\"positions_reported\":{},\"per_position\":{},\"gc_percent\":{},\"gc_histogram\":[{}],\"n_bases\":{},\"n_percent\":{},\"reads_with_n\":{},\"duplication_percent\":{},\"duplication_sample\":{},\"adapters\":[{}],\"likely_adapter\":{},\"warnings\":[{}]}}",
            json_texto(&self.path),
            self.reads,
            self.bases,
            self.complete,
            self.min_length,
            self.max_length,
            num(self.mean_length),
            junta(self.length_histogram.iter().map(|(l, c)| format!("[{l},{c}]"))),
            self.phred_offset,
            match self.quality_char_range {
                Some((menor, maior)) => format!("[{menor},{maior}]"),
                None => "null".to_string(),
            },
            num(self.mean_quality),
            junta(self.read_quality_histogram.iter()),
            self.positions_reported,
            por_posicao,
            num(self.gc_percent),
            junta(self.gc_histogram.iter()),
            self.n_bases,
            num(self.n_percent),
            self.reads_with_n,
            num(self.duplication_percent),
            self.duplication_sample,
            adaptadores,
            provavel,
            junta(self.warnings().iter().map(|a| json_texto(a)))
        )
    }
}

// ---------------------------------------------------------------------
//  PARTE 2: APARO E FILTRO
// ---------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct TrimOptions {
    /// Adaptador 3' da leitura 1 (`cutadapt -a`). `None`: nao procura.
    pub adapter: Option<String>,
    /// Adaptador 3' da leitura 2 (`cutadapt -A`). `None`: usa `adapter`.
    pub adapter2: Option<String>,
    /// Sobreposicao minima entre adaptador e leitura (`-O`).
    pub min_overlap: usize,
    /// Fracao de erros tolerada sobre o trecho casado (`-e`), de 0 a menos de 1.
    pub max_error_rate: f64,
    /// Corte de qualidade na ponta 3' (`-q`). Zero desliga.
    pub quality_cutoff: u8,
    /// Corte de qualidade na ponta 5' (`-q 5',3'`). Zero desliga.
    pub quality_cutoff_5: u8,
    /// Tira os N (maiusculos) das duas pontas (`--trim-n`).
    pub trim_n: bool,
    /// Descarta a leitura mais curta que isto, depois do aparo (`-m`).
    pub min_length: usize,
    /// Descarta a leitura com mais N que isto (`--max-n`).
    pub max_n: Option<usize>,
    /// Descarta a leitura cuja qualidade media fica abaixo disto.
    pub min_mean_quality: Option<u8>,
    /// 33 ou 64.
    pub phred_offset: u8,
}

impl Default for TrimOptions {
    fn default() -> TrimOptions {
        TrimOptions {
            adapter: None,
            adapter2: None,
            min_overlap: 3,
            max_error_rate: 0.1,
            quality_cutoff: 0,
            quality_cutoff_5: 0,
            trim_n: false,
            min_length: 0,
            max_n: None,
            min_mean_quality: None,
            phred_offset: 33,
        }
    }
}

/// O que o aparo fez. Em modo pareado, `reads_in`, `reads_out` e os
/// descartes contam PARES; as bases somam as duas leituras.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct TrimSummary {
    pub paired: bool,
    pub reads_in: u64,
    pub reads_out: u64,
    pub bases_in: u64,
    pub bases_out: u64,
    /// Leituras 1 em que o adaptador foi achado (inclusive as descartadas depois).
    pub reads_with_adapter: u64,
    /// O mesmo para as leituras 2. Zero em leitura simples.
    pub reads_with_adapter2: u64,
    pub bases_removed_quality: u64,
    pub bases_removed_adapter: u64,
    pub bases_removed_n: u64,
    /// Cada leitura (ou par) descartada conta no primeiro motivo que vale,
    /// nesta ordem.
    pub discarded_too_short: u64,
    pub discarded_too_many_n: u64,
    pub discarded_low_quality: u64,
}

impl TrimSummary {
    pub fn json(&self) -> String {
        format!(
            "{{\"paired\":{},\"reads_in\":{},\"reads_out\":{},\"bases_in\":{},\"bases_out\":{},\"reads_with_adapter\":{},\"reads_with_adapter2\":{},\"bases_removed_quality\":{},\"bases_removed_adapter\":{},\"bases_removed_n\":{},\"discarded_too_short\":{},\"discarded_too_many_n\":{},\"discarded_low_quality\":{}}}",
            self.paired,
            self.reads_in,
            self.reads_out,
            self.bases_in,
            self.bases_out,
            self.reads_with_adapter,
            self.reads_with_adapter2,
            self.bases_removed_quality,
            self.bases_removed_adapter,
            self.bases_removed_n,
            self.discarded_too_short,
            self.discarded_too_many_n,
            self.discarded_low_quality
        )
    }
}

/// Aparo por qualidade, o algoritmo do BWA como o cutadapt o escreve (`-q`).
/// Devolve o trecho `[inicio, fim)` que fica. As duas pontas sao calculadas
/// sobre a leitura original; se os cortes se cruzam, nao sobra nada. Corte
/// zero desliga a ponta.
pub fn quality_trim_bounds(qual: &[u8], cutoff_5: u8, cutoff_3: u8, phred_offset: u8) -> (usize, usize) {
    let n = qual.len();
    let base = i64::from(phred_offset);
    let (mut inicio, mut fim) = (0usize, n);
    if cutoff_5 > 0 {
        let (mut s, mut maximo) = (0i64, 0i64);
        for (i, &c) in qual.iter().enumerate() {
            s += i64::from(cutoff_5) - (i64::from(c) - base);
            if s < 0 {
                break;
            }
            if s > maximo {
                maximo = s;
                inicio = i + 1;
            }
        }
    }
    if cutoff_3 > 0 {
        let (mut s, mut maximo) = (0i64, 0i64);
        for (i, &c) in qual.iter().enumerate().rev() {
            s += i64::from(cutoff_3) - (i64::from(c) - base);
            if s < 0 {
                break;
            }
            if s > maximo {
                maximo = s;
                fim = i;
            }
        }
    }
    if inicio >= fim {
        (0, 0)
    } else {
        (inicio, fim)
    }
}

/// Um adaptador 3' pronto para procurar, com o resultado do
/// `cutadapt -a SEQ -e TAXA -O MIN --no-indels`.
#[derive(Clone, Debug)]
pub struct Adapter {
    /// Em maiusculas: A, C, G, T e N.
    seq: Vec<u8>,
    /// Tem N: o adaptador e comparado com curingas.
    curinga: bool,
    /// Ja limitada ao comprimento do adaptador.
    min_overlap: usize,
    /// floor(taxa * comprimento): o teto de trocas de qualquer casamento.
    teto: usize,
    /// permitidos[L]: trocas toleradas quando as L primeiras bases do
    /// adaptador se sobrepoem a leitura = floor(taxa * (L - quantas dessas L
    /// bases sao N)). O N casa com tudo e nao conta no comprimento.
    permitidos: Vec<usize>,
}

impl Adapter {
    pub fn new(sequence: &str, min_overlap: usize, max_error_rate: f64) -> Result<Adapter, String> {
        let seq: Vec<u8> = sequence.trim().bytes().map(|b| b.to_ascii_uppercase()).collect();
        if seq.is_empty() {
            return Err("The adapter sequence is empty: pass the bases of the adapter (A, C, G, T, N), or leave the adapter out to skip adapter removal".to_string());
        }
        if let Some(&c) = seq.iter().find(|c| !matches!(c, b'A' | b'C' | b'G' | b'T' | b'N')) {
            return Err(format!(
                "The adapter `{}` has the invalid character `{}`: an adapter is a sequence of A, C, G, T and N only (N matches any base). To find out which adapter a file has, run fastq_stats(...) and pass the sequence in its `likely_adapter`",
                mostra(sequence.as_bytes(), 60),
                char::from(c).escape_default()
            ));
        }
        if seq.iter().all(|&c| c == b'N') {
            return Err(format!(
                "The adapter `{}` has only N: it would match everywhere and remove every read entirely. Pass the real bases of the adapter",
                mostra(sequence.as_bytes(), 60)
            ));
        }
        if !(0.0..1.0).contains(&max_error_rate) {
            return Err(format!(
                "max_error_rate is {max_error_rate}, but it must be at least 0 and below 1: it is the fraction of the matched bases that may differ (0.1 allows one mismatch every 10 bases)"
            ));
        }
        if min_overlap == 0 {
            return Err("min_overlap is 0, but it must be at least 1: it is the smallest number of adapter bases that have to match at the end of a read (3 is the usual value)".to_string());
        }
        let m = seq.len();
        let mut permitidos = Vec::with_capacity(m + 1);
        permitidos.push(0);
        let mut sem_n = 0usize;
        for &c in &seq {
            sem_n += usize::from(c != b'N');
            permitidos.push((sem_n as f64 * max_error_rate) as usize);
        }
        Ok(Adapter {
            curinga: sem_n < m,
            // Um adaptador mais curto que a sobreposicao minima ainda vale
            // inteiro (o cutadapt faz o mesmo); sem isto ele nunca casaria.
            min_overlap: min_overlap.min(m),
            teto: (max_error_rate * m as f64) as usize,
            permitidos,
            seq,
        })
    }

    /// Trocas entre as `len` primeiras bases do adaptador e a leitura a partir
    /// de `p`. Para de contar ao passar de `limite`. N no adaptador casa com
    /// qualquer coisa; N na leitura e uma troca; maiusculas e minusculas sao a
    /// mesma base.
    #[inline]
    fn trocas(&self, read: &[u8], p: usize, len: usize, limite: usize) -> usize {
        let mut t = 0usize;
        for (&a, &r) in self.seq[..len].iter().zip(&read[p..p + len]) {
            if a == b'N' {
                continue;
            }
            let r = r.to_ascii_uppercase();
            // Com curingas o cutadapt le o U da leitura como T; sem eles
            // compara as letras.
            if a != r && !(self.curinga && a == b'T' && r == b'U') {
                t += 1;
                if t > limite {
                    break;
                }
            }
        }
        t
    }

    /// Onde comeca, na leitura, o prefixo mais longo do adaptador que o
    /// cutadapt ainda acompanha na ultima coluna da sua matriz: ele segue um
    /// prefixo enquanto as trocas nao passam de `teto`, e mais uma base.
    ///
    /// So serve a um detalhe de `find`: ao comparar um casamento parcial com
    /// uma ocorrencia inteira ja achada, o cutadapt 5.2 nao olha onde o
    /// parcial comeca, e sim este ponto (a variavel `origin` de
    /// `Aligner.locate` nao e atualizada naquele laco). Reproduzido como e,
    /// para o resultado ser o mesmo.
    fn origem_acompanhada(&self, read: &[u8]) -> usize {
        let (m, n) = (self.seq.len(), read.len());
        let acompanhado = if n < 2 {
            m.min(self.teto + 1)
        } else {
            let mut vivo = 0;
            for len in (1..=m.min(n - 1)).rev() {
                if self.trocas(read, n - 1 - len, len, self.teto) <= self.teto {
                    vivo = len;
                    break;
                }
            }
            m.min(vivo + 1)
        };
        n.saturating_sub(acompanhado)
    }

    /// Onde cortar a leitura: o que o `cutadapt -a ... --no-indels` responde.
    ///
    /// Um casamento e ACEITAVEL quando a sobreposicao L entre adaptador e
    /// leitura e pelo menos `min_overlap` e as trocas nao passam de
    /// floor(taxa * (L - N do adaptador nessas L bases)). Sem insercoes nem
    /// delecoes. A PONTUACAO de um casamento e acertos menos trocas.
    ///
    /// 1. Ocorrencias INTEIRAS do adaptador, da esquerda para a direita. Fica
    ///    a primeira aceitavel; uma posterior so a substitui se comeca a ate
    ///    meio adaptador de distancia dela e tem pontuacao maior. Uma
    ///    ocorrencia sem nenhuma troca, quando aceita, encerra a busca. E o
    ///    que faz um adaptador repetitivo (poli-A) parar no comeco exato da
    ///    cauda, em vez de levar junto as bases toleradas como troca.
    /// 2. O adaptador passando do fim da leitura (casamento parcial na ponta
    ///    3'), do mais longo para o mais curto. Sem ocorrencia inteira, fica o
    ///    primeiro aceitavel, e um mais curto o substitui se tem pontuacao
    ///    maior. Havendo ocorrencia inteira, um parcial so a substitui com
    ///    pontuacao maior e perto dela (ver `origem_acompanhada`).
    pub fn find(&self, read: &[u8]) -> Option<usize> {
        let (m, n) = (self.seq.len(), read.len());
        let meio = m / 2;
        // (pontuacao, posicao, e uma ocorrencia inteira achada no passo 1?)
        let mut melhor: Option<(i64, usize, bool)> = None;
        if n >= m {
            let limite = self.permitidos[m];
            for p in 0..=n - m {
                let t = self.trocas(read, p, m, limite);
                if t > limite {
                    continue;
                }
                let pontos = m as i64 - 2 * t as i64;
                let troca = match melhor {
                    None => true,
                    Some((s, o, _)) => p <= o + meio && pontos > s,
                };
                if troca {
                    if t == 0 {
                        return Some(p);
                    }
                    melhor = Some((pontos, p, true));
                }
            }
        }
        let mut acompanhada: Option<usize> = None;
        for len in (self.min_overlap..=m.min(n)).rev() {
            let limite = self.permitidos[len];
            let p = n - len;
            let t = self.trocas(read, p, len, limite);
            if t > limite {
                continue;
            }
            let pontos = len as i64 - 2 * t as i64;
            let troca = match melhor {
                None => true,
                Some((s, _, false)) => pontos > s,
                Some((s, o, true)) => pontos > s && *acompanhada.get_or_insert_with(|| self.origem_acompanhada(read)) <= o + meio,
            };
            if troca {
                melhor = Some((pontos, p, false));
            }
        }
        melhor.map(|(_, p, _)| p)
    }
}

/// N para o filtro `max_n`: maiusculo ou minusculo, como no cutadapt.
#[inline]
fn e_n(b: u8) -> bool {
    b == b'N' || b == b'n'
}

/// O que sobrou de uma leitura depois dos tres passos de aparo.
struct Aparo {
    inicio: usize,
    fim: usize,
    por_qualidade: usize,
    por_adaptador: usize,
    por_n: usize,
    achou: bool,
}

/// Os passos de aparo, NA ORDEM DO CUTADAPT: qualidade, adaptador, N.
fn apara(seq: &[u8], qual: &[u8], op: &TrimOptions, adaptador: Option<&Adapter>) -> Aparo {
    let (mut inicio, mut fim) = quality_trim_bounds(qual, op.quality_cutoff_5, op.quality_cutoff, op.phred_offset);
    let por_qualidade = seq.len() - (fim - inicio);
    let (mut por_adaptador, mut achou) = (0, false);
    if let Some(p) = adaptador.and_then(|a| a.find(&seq[inicio..fim])) {
        por_adaptador = fim - inicio - p;
        fim = inicio + p;
        achou = true;
    }
    let mut por_n = 0;
    if op.trim_n {
        let antes = fim - inicio;
        // So o N maiusculo, como o `--trim-n` do cutadapt.
        while inicio < fim && seq[inicio] == b'N' {
            inicio += 1;
        }
        while fim > inicio && seq[fim - 1] == b'N' {
            fim -= 1;
        }
        por_n = antes - (fim - inicio);
    }
    Aparo { inicio, fim, por_qualidade, por_adaptador, por_n, achou }
}

/// Os motivos de descarte que valem para uma leitura ja aparada: curta
/// demais, N demais, qualidade media baixa.
fn motivos(seq: &[u8], qual: &[u8], op: &TrimOptions) -> [bool; 3] {
    let curta = seq.len() < op.min_length;
    let muitos_n = op.max_n.is_some_and(|m| seq.iter().filter(|b| e_n(**b)).count() > m);
    let baixa = op.min_mean_quality.is_some_and(|q| {
        let soma: i64 = qual.iter().map(|&c| i64::from(c) - i64::from(op.phred_offset)).sum();
        qual.is_empty() || soma < i64::from(q) * qual.len() as i64
    });
    [curta, muitos_n, baixa]
}

fn grava(w: &mut Escritor, caminho: &str, cab: &[u8], seq: &[u8], qual: &[u8]) -> Result<(), String> {
    w.write_all(cab)
        .and_then(|_| w.write_all(b"\n"))
        .and_then(|_| w.write_all(seq))
        .and_then(|_| w.write_all(b"\n+\n"))
        .and_then(|_| w.write_all(qual))
        .and_then(|_| w.write_all(b"\n"))
        .map_err(|e| format!("Error writing {caminho}: {e}. Check the free space on the disk"))
}

/// O caminho sem atalhos, exista o arquivo ou nao.
fn caminho_absoluto(p: &str) -> PathBuf {
    let caminho = Path::new(p);
    if let Ok(c) = std::fs::canonicalize(caminho) {
        return c;
    }
    let pai = caminho.parent().filter(|d| !d.as_os_str().is_empty()).unwrap_or(Path::new("."));
    match (std::fs::canonicalize(pai), caminho.file_name()) {
        (Ok(d), Some(f)) => d.join(f),
        _ => caminho.to_path_buf(),
    }
}

fn mesmo_arquivo(a: &str, b: &str) -> bool {
    let (x, y) = (caminho_absoluto(a), caminho_absoluto(b));
    if cfg!(windows) {
        x.to_string_lossy().eq_ignore_ascii_case(&y.to_string_lossy())
    } else {
        x == y
    }
}

/// Soma ao resumo o que o aparo fez numa leitura.
fn conta_aparo(s: &mut TrimSummary, comprimento: usize, a: &Aparo) {
    s.bases_in += comprimento as u64;
    s.bases_removed_quality += a.por_qualidade as u64;
    s.bases_removed_adapter += a.por_adaptador as u64;
    s.bases_removed_n += a.por_n as u64;
}

/// Conta o descarte no primeiro motivo que vale. `false`: a leitura fica.
fn descarta(s: &mut TrimSummary, m: [bool; 3]) -> bool {
    if m[0] {
        s.discarded_too_short += 1;
    } else if m[1] {
        s.discarded_too_many_n += 1;
    } else if m[2] {
        s.discarded_low_quality += 1;
    } else {
        return false;
    }
    true
}

struct Aparador<'a> {
    op: &'a TrimOptions,
    ad1: Option<Adapter>,
    ad2: Option<Adapter>,
}

impl Aparador<'_> {
    fn simples(&self, l: &mut Leitor, w: &mut Escritor, saida: &str, s: &mut TrimSummary) -> Result<(), String> {
        while l.proximo()? {
            let a = apara(&l.seq, &l.qual, self.op, self.ad1.as_ref());
            s.reads_in += 1;
            conta_aparo(s, l.seq.len(), &a);
            s.reads_with_adapter += u64::from(a.achou);
            let (seq, qual) = (&l.seq[a.inicio..a.fim], &l.qual[a.inicio..a.fim]);
            if descarta(s, motivos(seq, qual, self.op)) {
                continue;
            }
            grava(w, saida, &l.cab, seq, qual)?;
            s.reads_out += 1;
            s.bases_out += seq.len() as u64;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn pareado(
        &self,
        l1: &mut Leitor,
        l2: &mut Leitor,
        w1: &mut Escritor,
        w2: &mut Escritor,
        saida1: &str,
        saida2: &str,
        s: &mut TrimSummary,
    ) -> Result<(), String> {
        loop {
            let (tem1, tem2) = (l1.proximo()?, l2.proximo()?);
            let k = s.reads_in + 1;
            match (tem1, tem2) {
                (false, false) => return Ok(()),
                (true, true) => {}
                (true, false) => {
                    return Err(format!(
                        "The two files are out of step at record {k}: {} has read `{}` there, but {} ends after {} records (read 2: none, end of file). Both files must have the same reads in the same order; one of them is probably incomplete or was filtered on its own",
                        l1.caminho,
                        l1.nome(),
                        l2.caminho,
                        l2.registros
                    ));
                }
                (false, true) => {
                    return Err(format!(
                        "The two files are out of step at record {k}: {} has read `{}` there, but {} ends after {} records (read 1: none, end of file). Both files must have the same reads in the same order; one of them is probably incomplete or was filtered on its own",
                        l2.caminho,
                        l2.nome(),
                        l1.caminho,
                        l1.registros
                    ));
                }
            }
            if nome_do_par(&l1.cab) != nome_do_par(&l2.cab) {
                return Err(format!(
                    "The two files are out of step at record {k} (line {} of {}): read 1 is `{}` and read 2 is `{}`. Both files must list the same reads in the same order; use the two files exactly as they came from the sequencer, or re-pair them first",
                    l1.linha.saturating_sub(3),
                    l1.caminho,
                    l1.nome(),
                    l2.nome()
                ));
            }
            let a1 = apara(&l1.seq, &l1.qual, self.op, self.ad1.as_ref());
            let a2 = apara(&l2.seq, &l2.qual, self.op, self.ad2.as_ref());
            s.reads_in += 1;
            conta_aparo(s, l1.seq.len(), &a1);
            conta_aparo(s, l2.seq.len(), &a2);
            s.reads_with_adapter += u64::from(a1.achou);
            s.reads_with_adapter2 += u64::from(a2.achou);
            let (seq1, qual1) = (&l1.seq[a1.inicio..a1.fim], &l1.qual[a1.inicio..a1.fim]);
            let (seq2, qual2) = (&l2.seq[a2.inicio..a2.fim], &l2.qual[a2.inicio..a2.fim]);
            // O par sai quando QUALQUER das duas leituras e reprovada
            // (`--pair-filter=any`), e conta no primeiro motivo que vale para
            // uma delas.
            let (m1, m2) = (motivos(seq1, qual1, self.op), motivos(seq2, qual2, self.op));
            if descarta(s, [m1[0] || m2[0], m1[1] || m2[1], m1[2] || m2[2]]) {
                continue;
            }
            grava(w1, saida1, &l1.cab, seq1, qual1)?;
            grava(w2, saida2, &l2.cab, seq2, qual2)?;
            s.reads_out += 1;
            s.bases_out += (seq1.len() + seq2.len()) as u64;
        }
    }
}

/// Apara e filtra um FASTQ (`entrada2` ausente) ou um par de FASTQ.
///
/// Cada leitura passa, nesta ordem (a do cutadapt): 1. aparo por qualidade
/// (`quality_cutoff`, `quality_cutoff_5`); 2. remocao do adaptador 3'
/// (`adapter` na leitura 1; `adapter2`, ou `adapter` se ele faltar, na
/// leitura 2); 3. aparo de N nas pontas (`trim_n`); 4. filtros (`min_length`,
/// `max_n`, `min_mean_quality`). Em modo pareado o par e descartado quando
/// qualquer das duas leituras e reprovada, e as duas saidas ficam com os
/// mesmos registros, na mesma ordem.
///
/// Entrada em texto, `.gz` ou `.lz4`, com fim de linha LF ou CRLF; a extensao
/// de cada saida escolhe a compressao, e o fim de linha e sempre LF. Uma
/// leitura aparada ate o tamanho zero que passa nos filtros e gravada com as
/// linhas de sequencia e de qualidade vazias.
pub fn fastq_trim(
    entrada1: &str,
    saida1: &str,
    entrada2: Option<&str>,
    saida2: Option<&str>,
    opt: &TrimOptions,
) -> Result<TrimSummary, String> {
    if opt.phred_offset != 33 && opt.phred_offset != 64 {
        return Err(format!(
            "phred_offset is {}, but it must be 33 (the current standard) or 64 (Illumina 1.3 to 1.7); fastq_stats(...) tells which one a file uses",
            opt.phred_offset
        ));
    }
    let saida2 = match (entrada2, saida2) {
        (Some(_), None) => {
            return Err("A second input file was given without a second output file: in paired-end mode both outputs are needed, one for each read of the pair".to_string());
        }
        (None, Some(s)) => {
            return Err(format!(
                "A second output file (`{s}`) was given without a second input file: pass both files of the pair, or leave the second output out for single-end reads"
            ));
        }
        (_, s) => s,
    };
    if entrada2.is_none() && opt.adapter2.is_some() {
        return Err("adapter2 was given without a second input file: adapter2 is the adapter of read 2 of a pair. For single-end reads pass the adapter in `adapter`".to_string());
    }
    let adaptador = |s: &Option<String>| -> Result<Option<Adapter>, String> {
        s.as_deref().map(|a| Adapter::new(a, opt.min_overlap, opt.max_error_rate)).transpose()
    };
    let ad1 = adaptador(&opt.adapter)?;
    let ad2 = match (entrada2, &opt.adapter2) {
        (None, _) => None,
        (Some(_), Some(_)) => adaptador(&opt.adapter2)?,
        (Some(_), None) => ad1.clone(),
    };

    // Uma saida que e um dos arquivos de entrada o apagaria antes de ser lido.
    let entradas: Vec<&str> = std::iter::once(entrada1).chain(entrada2).collect();
    let saidas: Vec<&str> = std::iter::once(saida1).chain(saida2).collect();
    for s in &saidas {
        if let Some(e) = entradas.iter().find(|e| mesmo_arquivo(s, e)) {
            return Err(format!(
                "The output `{s}` is the input file `{e}` itself: writing it would destroy the input before it is read. Choose another name for the output"
            ));
        }
    }
    if let Some(s2) = saida2 {
        if mesmo_arquivo(saida1, s2) {
            return Err(format!(
                "The two outputs are the same file (`{saida1}`): each read of the pair needs its own output file"
            ));
        }
    }

    // As entradas abrem antes de as saidas serem criadas: um caminho de
    // entrada errado nao deixa um arquivo de saida vazio para tras.
    let mut l1 = Leitor::abre(entrada1)?;
    let mut l2 = entrada2.map(Leitor::abre).transpose()?;
    let mut w1 = Escritor::cria(saida1)?;
    let mut w2 = saida2.map(Escritor::cria).transpose()?;

    let aparador = Aparador { op: opt, ad1, ad2 };
    let mut s = TrimSummary { paired: l2.is_some(), ..TrimSummary::default() };
    let feito = match (l2.as_mut(), w2.as_mut(), saida2) {
        (Some(l2), Some(w2), Some(s2)) => aparador.pareado(&mut l1, l2, &mut w1, w2, saida1, s2, &mut s),
        _ => aparador.simples(&mut l1, &mut w1, saida1, &mut s),
    };
    // Um erro no meio deixa a saida pela metade: quem ler a mensagem tem de
    // saber que ela nao serve.
    feito.map_err(|e| format!("{e}. The output written so far is incomplete: do not use it"))?;
    w1.termina().map_err(|e| format!("Error writing {saida1}: {e}"))?;
    if let (Some(w2), Some(s2)) = (w2, saida2) {
        w2.termina().map_err(|e| format!("Error writing {s2}: {e}"))?;
    }
    Ok(s)
}
