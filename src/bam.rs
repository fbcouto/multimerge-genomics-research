// =====================================================================
//  BAM -- o SAM em binario, dentro de um BGZF
//
//  Aqui o BAM e tratado como um ENVELOPE do SAM: o leitor entrega as mesmas
//  linhas de texto que um `.sam` teria (cabecalho e alinhamentos), e o
//  escritor recebe linhas de SAM e grava o binario. Com isso tudo o que o
//  projeto faz com SAM -- validar, ordenar, normalizar nomes -- vale para BAM
//  sem que cada ferramenta precise saber que o formato existe.
//
//  O que o binario guarda e o texto nao: o indice numerico da sequencia de
//  referencia (por isso todo RNAME precisa estar declarado num @SQ), o `bin`
//  do esquema de indexacao e os tipos dos campos opcionais. Um inteiro e
//  gravado no menor tipo em que cabe, como o htslib faz.
// =====================================================================

use crate::bgzf::{BgzfParLeitor, BgzfWriter};
use std::collections::HashMap;
use std::io::{self, Read, Write};

fn erro(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// O cabecalho de um BAM: o texto (o mesmo de um SAM) e a lista binaria de
/// sequencias de referencia, com nome e comprimento.
#[derive(Clone, Debug, Default)]
pub struct Cabecalho {
    pub texto: String,
    pub refs: Vec<(String, u32)>,
}

/// O necessario de um alinhamento para indexar ou filtrar por regiao.
#[derive(Clone, Copy, Debug)]
pub struct Info {
    pub tid: i32,
    /// Inicio, base 0. -1 quando nao ha posicao.
    pub pos: i64,
    /// Fim exclusivo sobre a referencia. Sem CIGAR, ou nao mapeado: `pos + 1`.
    pub fim: i64,
    pub flag: u16,
}

// ---------------------------------------------------------------------
//  O BIN (esquema de indexacao do UCSC: 14 bits no nivel mais fino, 5 niveis)
// ---------------------------------------------------------------------

/// O menor bin que contem `[beg, end)`. Para um alinhamento sem posicao
/// (`beg = -1, end = 0`) da 4680, o valor que o samtools grava.
pub fn reg2bin(beg: i64, end: i64) -> u32 {
    let end = end - 1;
    if beg >> 14 == end >> 14 {
        return (4681 + (beg >> 14)) as u32;
    }
    if beg >> 17 == end >> 17 {
        return (585 + (beg >> 17)) as u32;
    }
    if beg >> 20 == end >> 20 {
        return (73 + (beg >> 20)) as u32;
    }
    if beg >> 23 == end >> 23 {
        return (9 + (beg >> 23)) as u32;
    }
    if beg >> 26 == end >> 26 {
        return (1 + (beg >> 26)) as u32;
    }
    0
}

// ---------------------------------------------------------------------
//  LEITURA DO BINARIO
// ---------------------------------------------------------------------

/// Le `buf.len()` bytes. `false` se a entrada acabou antes do primeiro;
/// acabar no meio e erro.
pub fn le_exato_ou_fim<R: Read>(r: &mut R, buf: &mut [u8]) -> io::Result<bool> {
    let mut lidos = 0;
    while lidos < buf.len() {
        match r.read(&mut buf[lidos..]) {
            Ok(0) if lidos == 0 => return Ok(false),
            Ok(0) => return Err(erro("truncated BAM: the file ends in the middle of a record".into())),
            Ok(n) => lidos += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        }
    }
    Ok(true)
}

fn le_i32<R: Read>(r: &mut R, o_que: &str) -> io::Result<i32> {
    let mut b = [0u8; 4];
    if !le_exato_ou_fim(r, &mut b)? {
        return Err(erro(format!("truncated BAM: the file ends before {o_que}")));
    }
    Ok(i32::from_le_bytes(b))
}

/// Le `n` bytes sem reservar a memoria de antemao: um tamanho corrompido nao
/// vira uma alocacao gigante.
fn le_bytes<R: Read>(r: &mut R, n: usize, o_que: &str) -> io::Result<Vec<u8>> {
    let mut v = Vec::new();
    r.take(n as u64).read_to_end(&mut v)?;
    if v.len() != n {
        return Err(erro(format!("truncated BAM: the file ends in the middle of {o_que}")));
    }
    Ok(v)
}

pub fn le_cabecalho<R: Read>(r: &mut R) -> io::Result<Cabecalho> {
    let mut magico = [0u8; 4];
    if !le_exato_ou_fim(r, &mut magico)? || &magico != b"BAM\x01" {
        return Err(erro("not a BAM file: the BAM signature is missing at the start".into()));
    }
    let l_texto = le_i32(r, "the header")?;
    if l_texto < 0 {
        return Err(erro("damaged BAM header: negative text length".into()));
    }
    let mut texto = le_bytes(r, l_texto as usize, "the header text")?;
    while texto.last() == Some(&0) {
        texto.pop();
    }
    let n_ref = le_i32(r, "the reference list")?;
    if n_ref < 0 {
        return Err(erro("damaged BAM header: negative number of reference sequences".into()));
    }
    let mut refs = Vec::new();
    for _ in 0..n_ref {
        let l_nome = le_i32(r, "a reference name")?;
        if l_nome < 1 {
            return Err(erro("damaged BAM header: a reference sequence has no name".into()));
        }
        let mut nome = le_bytes(r, l_nome as usize, "a reference name")?;
        while nome.last() == Some(&0) {
            nome.pop();
        }
        let l_ref = le_i32(r, "a reference length")?;
        refs.push((String::from_utf8_lossy(&nome).to_string(), l_ref as u32));
    }
    Ok(Cabecalho { texto: String::from_utf8_lossy(&texto).to_string(), refs })
}

/// Le o proximo alinhamento (sem o campo de tamanho) para `rec`. `false` no
/// fim do arquivo.
pub fn le_registro<R: Read>(r: &mut R, rec: &mut Vec<u8>) -> io::Result<bool> {
    let mut b = [0u8; 4];
    if !le_exato_ou_fim(r, &mut b)? {
        return Ok(false);
    }
    let tamanho = i32::from_le_bytes(b);
    if tamanho < 32 {
        return Err(erro(format!("damaged BAM: a record declares an impossible size ({tamanho})")));
    }
    rec.clear();
    r.take(tamanho as u64).read_to_end(rec)?;
    if rec.len() != tamanho as usize {
        return Err(erro("truncated BAM: the file ends in the middle of a record".into()));
    }
    Ok(true)
}

#[inline]
fn i32_em(rec: &[u8], p: usize) -> i32 {
    i32::from_le_bytes([rec[p], rec[p + 1], rec[p + 2], rec[p + 3]])
}

#[inline]
fn u16_em(rec: &[u8], p: usize) -> u16 {
    u16::from_le_bytes([rec[p], rec[p + 1]])
}

const OPS_CIGAR: &[u8; 9] = b"MIDNSHP=X";

/// Comprimento sobre a referencia de um CIGAR binario: M, D, N, = e X.
fn comprimento_ref(cigar: &[u8]) -> i64 {
    cigar
        .chunks_exact(4)
        .map(|c| {
            let v = u32::from_le_bytes([c[0], c[1], c[2], c[3]]);
            match v & 0xf {
                0 | 2 | 3 | 7 | 8 => i64::from(v >> 4),
                _ => 0,
            }
        })
        .sum()
}

pub fn info(rec: &[u8]) -> Result<Info, String> {
    if rec.len() < 32 {
        return Err("damaged BAM record: shorter than its fixed part".into());
    }
    let (tid, pos, flag) = (i32_em(rec, 0), i64::from(i32_em(rec, 4)), u16_em(rec, 14));
    let ini = 32 + rec[8] as usize;
    let fim_cigar = ini + 4 * u16_em(rec, 12) as usize;
    if fim_cigar > rec.len() {
        return Err("damaged BAM record: shorter than its fields say".into());
    }
    let mut l = if flag & 4 != 0 { 0 } else { comprimento_ref(&rec[ini..fim_cigar]) };
    if l == 0 {
        l = 1;
    }
    Ok(Info { tid, pos, fim: pos + l, flag })
}

/// Tamanho, em bytes, de um elemento de campo opcional do tipo dado.
fn tamanho_tipo(t: u8) -> Option<usize> {
    match t {
        b'A' | b'c' | b'C' => Some(1),
        b's' | b'S' => Some(2),
        b'i' | b'I' | b'f' => Some(4),
        b'd' => Some(8),
        _ => None,
    }
}

/// Onde termina o campo opcional que comeca em `p` (tag, tipo e valor).
fn fim_do_campo(aux: &[u8], p: usize) -> Result<usize, String> {
    let curto = || "damaged BAM record: an optional field runs past the end of the record".to_string();
    if p + 3 > aux.len() {
        return Err(curto());
    }
    let v = p + 3;
    let fim = match aux[p + 2] {
        b'Z' | b'H' => v + aux[v..].iter().position(|&b| b == 0).ok_or_else(curto)? + 1,
        b'B' => {
            if v + 5 > aux.len() {
                return Err(curto());
            }
            let n = i32_em(aux, v + 1).max(0) as usize;
            let t = tamanho_tipo(aux[v]).filter(|_| aux[v] != b'A').ok_or_else(|| {
                format!("damaged BAM record: unknown array type '{}' in field {}", aux[v] as char, String::from_utf8_lossy(&aux[p..p + 2]))
            })?;
            v + 5 + n * t
        }
        t => {
            v + tamanho_tipo(t).ok_or_else(|| {
                format!("damaged BAM record: unknown type '{}' in field {}", t as char, String::from_utf8_lossy(&aux[p..p + 2]))
            })?
        }
    };
    if fim > aux.len() {
        return Err(curto());
    }
    Ok(fim)
}

/// Escreve um numero como o `%g` do C (6 algarismos significativos), que e
/// como o samtools imprime os campos de ponto flutuante.
fn escreve_g(x: f64, s: &mut Vec<u8>) {
    if x.is_nan() {
        s.extend_from_slice(b"nan");
        return;
    }
    if x.is_infinite() {
        s.extend_from_slice(if x < 0.0 { b"-inf" } else { b"inf" });
        return;
    }
    if x == 0.0 {
        s.extend_from_slice(if x.is_sign_negative() { b"-0" } else { b"0" });
        return;
    }
    let cientifico = format!("{x:.5e}");
    let (mantissa, expoente) = cientifico.split_once('e').unwrap_or((&cientifico, "0"));
    let expoente: i32 = expoente.parse().unwrap_or(0);
    let enxuga = |t: &str| -> String {
        if t.contains('.') {
            t.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            t.to_string()
        }
    };
    if (-4..6).contains(&expoente) {
        let casas = (5 - expoente).max(0) as usize;
        s.extend_from_slice(enxuga(&format!("{x:.casas$}")).as_bytes());
    } else {
        s.extend_from_slice(enxuga(mantissa).as_bytes());
        s.extend_from_slice(format!("e{}{:02}", if expoente < 0 { '-' } else { '+' }, expoente.abs()).as_bytes());
    }
}

fn escreve_int(v: i64, s: &mut Vec<u8>) {
    s.extend_from_slice(v.to_string().as_bytes());
}

/// Escreve o valor de um elemento de tipo numerico.
fn escreve_numero(t: u8, b: &[u8], s: &mut Vec<u8>) {
    match t {
        b'c' => escreve_int(i64::from(b[0] as i8), s),
        b'C' => escreve_int(i64::from(b[0]), s),
        b's' => escreve_int(i64::from(i16::from_le_bytes([b[0], b[1]])), s),
        b'S' => escreve_int(i64::from(u16::from_le_bytes([b[0], b[1]])), s),
        b'i' => escreve_int(i64::from(i32::from_le_bytes([b[0], b[1], b[2], b[3]])), s),
        b'I' => escreve_int(i64::from(u32::from_le_bytes([b[0], b[1], b[2], b[3]])), s),
        b'f' => escreve_g(f64::from(f32::from_le_bytes([b[0], b[1], b[2], b[3]])), s),
        _ => escreve_g(f64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]), s),
    }
}

fn escreve_cigar(cigar: &[u8], s: &mut Vec<u8>) -> Result<(), String> {
    for c in cigar.chunks_exact(4) {
        let v = u32::from_le_bytes([c[0], c[1], c[2], c[3]]);
        let op = *OPS_CIGAR.get((v & 0xf) as usize).ok_or("damaged BAM record: unknown CIGAR operation")?;
        escreve_int(i64::from(v >> 4), s);
        s.push(op);
    }
    Ok(())
}

fn nome_da_ref<'a>(refs: &'a [(String, u32)], id: i32) -> Result<&'a str, String> {
    refs.get(id as usize)
        .map(|r| r.0.as_str())
        .ok_or_else(|| format!("damaged BAM record: reference sequence number {id} is not in the header (it has {})", refs.len()))
}

/// Acrescenta a `s` a linha de SAM do alinhamento `rec`, sem o fim de linha.
pub fn decodifica(rec: &[u8], refs: &[(String, u32)], s: &mut Vec<u8>) -> Result<(), String> {
    if rec.len() < 32 {
        return Err("damaged BAM record: shorter than its fixed part".into());
    }
    let (refid, pos) = (i32_em(rec, 0), i32_em(rec, 4));
    let (l_nome, mapq) = (rec[8] as usize, rec[9]);
    let (n_cigar, flag) = (u16_em(rec, 12) as usize, u16_em(rec, 14));
    let l_seq = i32_em(rec, 16);
    let (prox_refid, prox_pos, tlen) = (i32_em(rec, 20), i32_em(rec, 24), i32_em(rec, 28));
    if l_seq < 0 {
        return Err("damaged BAM record: negative sequence length".into());
    }
    let l_seq = l_seq as usize;
    let fim_nome = 32 + l_nome;
    let fim_cigar = fim_nome + 4 * n_cigar;
    let fim_seq = fim_cigar + (l_seq + 1) / 2;
    let fim_qual = fim_seq + l_seq;
    if fim_qual > rec.len() {
        return Err("damaged BAM record: shorter than its fields say".into());
    }
    let aux = &rec[fim_qual..];

    // CIGAR com mais de 65535 operacoes: o campo fixo traz um substituto
    // (`<l_seq>S<comprimento>N`) e o CIGAR de verdade vem no campo CG.
    let mut cigar = &rec[fim_nome..fim_cigar];
    let mut cg: Option<(usize, usize)> = None;
    if n_cigar > 0 && refid >= 0 && pos >= 0 {
        let primeiro = u32::from_le_bytes([cigar[0], cigar[1], cigar[2], cigar[3]]);
        if primeiro & 0xf == 4 && (primeiro >> 4) as usize == l_seq {
            let mut p = 0;
            while p < aux.len() {
                let f = fim_do_campo(aux, p)?;
                if &aux[p..p + 2] == b"CG" && aux[p + 2] == b'B' && (aux[p + 3] == b'I' || aux[p + 3] == b'i') {
                    let n = i32_em(aux, p + 4).max(0) as usize;
                    if n >= n_cigar {
                        cigar = &aux[p + 8..p + 8 + 4 * n];
                        cg = Some((p, f));
                    }
                    break;
                }
                p = f;
            }
        }
    }

    // QNAME FLAG RNAME POS MAPQ
    let nome = &rec[32..fim_nome];
    let nome = &nome[..nome.iter().position(|&b| b == 0).unwrap_or(nome.len())];
    s.extend_from_slice(if nome.is_empty() { b"*" } else { nome });
    s.push(b'\t');
    escreve_int(i64::from(flag), s);
    s.push(b'\t');
    s.extend_from_slice(if refid < 0 { "*" } else { nome_da_ref(refs, refid)? }.as_bytes());
    s.push(b'\t');
    escreve_int(i64::from(pos) + 1, s);
    s.push(b'\t');
    escreve_int(i64::from(mapq), s);
    s.push(b'\t');
    // CIGAR
    if cigar.is_empty() {
        s.push(b'*');
    } else {
        escreve_cigar(cigar, s)?;
    }
    s.push(b'\t');
    // RNEXT PNEXT TLEN
    if prox_refid < 0 {
        s.push(b'*');
    } else if prox_refid == refid {
        s.push(b'=');
    } else {
        s.extend_from_slice(nome_da_ref(refs, prox_refid)?.as_bytes());
    }
    s.push(b'\t');
    escreve_int(i64::from(prox_pos) + 1, s);
    s.push(b'\t');
    escreve_int(i64::from(tlen), s);
    s.push(b'\t');
    // SEQ QUAL
    if l_seq == 0 {
        s.extend_from_slice(b"*\t*");
    } else {
        let seq = &rec[fim_cigar..fim_seq];
        for i in 0..l_seq {
            let b = seq[i / 2];
            s.push(b"=ACMGRSVTWYHKDBN"[(if i % 2 == 0 { b >> 4 } else { b & 0xf }) as usize]);
        }
        s.push(b'\t');
        let qual = &rec[fim_seq..fim_qual];
        if qual[0] == 0xff {
            s.push(b'*');
        } else {
            s.extend(qual.iter().map(|&q| q.saturating_add(33)));
        }
    }
    // Campos opcionais
    let mut p = 0;
    while p < aux.len() {
        let f = fim_do_campo(aux, p)?;
        if cg == Some((p, f)) {
            p = f;
            continue;
        }
        s.push(b'\t');
        s.extend_from_slice(&aux[p..p + 2]);
        s.push(b':');
        let (tipo, v) = (aux[p + 2], &aux[p + 3..f]);
        match tipo {
            b'A' => {
                s.extend_from_slice(b"A:");
                s.push(v[0]);
            }
            b'Z' | b'H' => {
                s.push(tipo);
                s.push(b':');
                s.extend_from_slice(&v[..v.len() - 1]);
            }
            b'B' => {
                s.extend_from_slice(b"B:");
                s.push(v[0]);
                let t = tamanho_tipo(v[0]).unwrap_or(1);
                for e in v[5..].chunks_exact(t) {
                    s.push(b',');
                    escreve_numero(v[0], e, s);
                }
            }
            b'f' | b'd' => {
                s.push(tipo);
                s.push(b':');
                escreve_numero(tipo, v, s);
            }
            _ => {
                s.extend_from_slice(b"i:");
                escreve_numero(tipo, v, s);
            }
        }
        p = f;
    }
    Ok(())
}

// ---------------------------------------------------------------------
//  ESCRITA DO BINARIO
// ---------------------------------------------------------------------

#[inline]
fn codigo_base(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'=' => 0,
        b'A' => 1,
        b'C' => 2,
        b'M' => 3,
        b'G' => 4,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
        b'T' => 8,
        b'W' => 9,
        b'Y' => 10,
        b'H' => 11,
        b'K' => 12,
        b'D' => 13,
        b'B' => 14,
        _ => 15,
    }
}

/// Grava um inteiro no menor tipo em que cabe, com a letra do tipo antes.
fn grava_menor_int(v: i64, s: &mut Vec<u8>) -> Result<(), String> {
    if v < 0 {
        if v >= i64::from(i8::MIN) {
            s.push(b'c');
            s.push(v as i8 as u8);
        } else if v >= i64::from(i16::MIN) {
            s.push(b's');
            s.extend_from_slice(&(v as i16).to_le_bytes());
        } else if v >= i64::from(i32::MIN) {
            s.push(b'i');
            s.extend_from_slice(&(v as i32).to_le_bytes());
        } else {
            return Err(format!("integer {v} does not fit in 32 bits"));
        }
    } else if v <= i64::from(u8::MAX) {
        s.push(b'C');
        s.push(v as u8);
    } else if v <= i64::from(u16::MAX) {
        s.push(b'S');
        s.extend_from_slice(&(v as u16).to_le_bytes());
    } else if v <= i64::from(u32::MAX) {
        s.push(b'I');
        s.extend_from_slice(&(v as u32).to_le_bytes());
    } else {
        return Err(format!("integer {v} does not fit in 32 bits"));
    }
    Ok(())
}

/// Grava um elemento de vetor do tipo `t`, lido do texto `v`.
fn grava_elemento(t: u8, v: &str, s: &mut Vec<u8>) -> Result<(), String> {
    let fora = || format!("value '{v}' does not fit the array type '{}'", t as char);
    match t {
        b'c' => s.push(v.parse::<i8>().map_err(|_| fora())? as u8),
        b'C' => s.push(v.parse::<u8>().map_err(|_| fora())?),
        b's' => s.extend_from_slice(&v.parse::<i16>().map_err(|_| fora())?.to_le_bytes()),
        b'S' => s.extend_from_slice(&v.parse::<u16>().map_err(|_| fora())?.to_le_bytes()),
        b'i' => s.extend_from_slice(&v.parse::<i32>().map_err(|_| fora())?.to_le_bytes()),
        b'I' => s.extend_from_slice(&v.parse::<u32>().map_err(|_| fora())?.to_le_bytes()),
        b'f' => s.extend_from_slice(&v.parse::<f32>().map_err(|_| fora())?.to_le_bytes()),
        _ => return Err(format!("unknown array type '{}'", t as char)),
    }
    Ok(())
}

/// Converte linhas de SAM em registros BAM. Guarda o numero de cada
/// sequencia declarada no cabecalho.
pub struct Codificador {
    ids: HashMap<String, i32>,
}

impl Codificador {
    pub fn novo(refs: &[(String, u32)]) -> Codificador {
        Codificador { ids: refs.iter().enumerate().map(|(i, r)| (r.0.clone(), i as i32)).collect() }
    }

    fn id(&self, nome: &str) -> Result<i32, String> {
        if nome == "*" {
            return Ok(-1);
        }
        self.ids.get(nome).copied().ok_or_else(|| {
            format!("sequence '{nome}' is not declared in the header: BAM needs an @SQ line (with SN and LN) for every sequence the alignments use")
        })
    }

    /// Acrescenta a `s` o registro BAM da linha (com o campo de tamanho).
    pub fn codifica(&self, linha: &str, s: &mut Vec<u8>) -> Result<(), String> {
        let mut c = linha.split('\t');
        let mut campo = |nome: &str| c.next().ok_or_else(|| format!("the {nome} column is missing (an alignment has 11 mandatory columns)"));
        let qname = campo("QNAME")?;
        let flag: u16 = campo("FLAG")?.parse().map_err(|_| "FLAG is not a number between 0 and 65535".to_string())?;
        let rname = campo("RNAME")?;
        let pos: i64 = campo("POS")?.parse().map_err(|_| "POS is not a number".to_string())?;
        let mapq: u8 = campo("MAPQ")?.parse().map_err(|_| "MAPQ is not a number between 0 and 255".to_string())?;
        let cigar_txt = campo("CIGAR")?;
        let rnext = campo("RNEXT")?;
        let pnext: i64 = campo("PNEXT")?.parse().map_err(|_| "PNEXT is not a number".to_string())?;
        let tlen: i32 = campo("TLEN")?.parse().map_err(|_| "TLEN is not a number that fits in 32 bits".to_string())?;
        let seq = campo("SEQ")?;
        let qual = campo("QUAL")?;

        if qname.is_empty() || qname.len() > 254 {
            return Err(format!("QNAME has {} characters; BAM allows 1 to 254", qname.len()));
        }
        if !(0..=i64::from(i32::MAX)).contains(&pos) || !(0..=i64::from(i32::MAX)).contains(&pnext) {
            return Err("POS or PNEXT is outside the range BAM can store (0 to 2147483647)".into());
        }
        let refid = self.id(rname)?;
        let prox_refid = if rnext == "=" { refid } else { self.id(rnext)? };

        // CIGAR
        let mut cigar: Vec<u32> = Vec::new();
        if cigar_txt != "*" {
            let mut n: u64 = 0;
            let mut tem_digito = false;
            for b in cigar_txt.bytes() {
                if b.is_ascii_digit() {
                    n = n * 10 + u64::from(b - b'0');
                    tem_digito = true;
                    if n >= 1 << 28 {
                        return Err(format!("CIGAR '{}' has an operation longer than BAM can store", resumo(cigar_txt)));
                    }
                } else {
                    let op = OPS_CIGAR.iter().position(|&o| o == b).filter(|_| tem_digito);
                    let Some(op) = op else {
                        return Err(format!("CIGAR '{}' is not valid", resumo(cigar_txt)));
                    };
                    cigar.push((n as u32) << 4 | op as u32);
                    n = 0;
                    tem_digito = false;
                }
            }
            if tem_digito {
                return Err(format!("CIGAR '{}' is not valid", resumo(cigar_txt)));
            }
        }
        let cigar_bytes: Vec<u8> = cigar.iter().flat_map(|v| v.to_le_bytes()).collect();
        let l_ref = if flag & 4 != 0 { 0 } else { comprimento_ref(&cigar_bytes) };
        let pos0 = pos - 1;
        let bin = reg2bin(pos0, pos0 + l_ref.max(1));

        let l_seq = if seq == "*" { 0 } else { seq.len() };
        if qual != "*" && qual.len() != l_seq {
            return Err(format!("SEQ has {l_seq} bases and QUAL has {} values; they must have the same length", qual.len()));
        }

        let inicio = s.len();
        s.extend_from_slice(&[0u8; 4]);
        s.extend_from_slice(&refid.to_le_bytes());
        s.extend_from_slice(&(pos0 as i32).to_le_bytes());
        s.push(qname.len() as u8 + 1);
        s.push(mapq);
        s.extend_from_slice(&(bin as u16).to_le_bytes());
        let longo = cigar.len() > 0xffff;
        s.extend_from_slice(&(if longo { 2u16 } else { cigar.len() as u16 }).to_le_bytes());
        s.extend_from_slice(&flag.to_le_bytes());
        s.extend_from_slice(&(l_seq as i32).to_le_bytes());
        s.extend_from_slice(&prox_refid.to_le_bytes());
        s.extend_from_slice(&((pnext - 1) as i32).to_le_bytes());
        s.extend_from_slice(&tlen.to_le_bytes());
        s.extend_from_slice(qname.as_bytes());
        s.push(0);
        if longo {
            // O substituto do CIGAR longo: tudo recortado, e um salto do
            // comprimento certo sobre a referencia.
            s.extend_from_slice(&((l_seq as u32) << 4 | 4).to_le_bytes());
            s.extend_from_slice(&((comprimento_ref(&cigar_bytes) as u32) << 4 | 3).to_le_bytes());
        } else {
            s.extend_from_slice(&cigar_bytes);
        }
        let sb = seq.as_bytes();
        for par in sb[..l_seq].chunks(2) {
            s.push(codigo_base(par[0]) << 4 | if par.len() > 1 { codigo_base(par[1]) } else { 0 });
        }
        if qual == "*" {
            s.extend(std::iter::repeat(0xffu8).take(l_seq));
        } else {
            for q in qual.bytes() {
                if !(33..=126).contains(&q) {
                    return Err("QUAL has a character outside the printable range".into());
                }
                s.push(q - 33);
            }
        }

        // Campos opcionais: TAG:TIPO:VALOR
        for campo in c {
            let b = campo.as_bytes();
            if b.len() < 5 || b[2] != b':' || b[4] != b':' {
                return Err(format!("optional field '{}' is not in the TAG:TYPE:VALUE form", resumo(campo)));
            }
            let v = &campo[5..];
            let invalido = || format!("optional field '{}' has a value that does not match its type", resumo(campo));
            s.extend_from_slice(&b[..2]);
            match b[3] {
                b'A' => {
                    if v.len() != 1 {
                        return Err(invalido());
                    }
                    s.push(b'A');
                    s.push(v.as_bytes()[0]);
                }
                b'i' => grava_menor_int(v.parse().map_err(|_| invalido())?, s).map_err(|e| format!("optional field '{}': {e}", resumo(campo)))?,
                b'f' => {
                    s.push(b'f');
                    s.extend_from_slice(&v.parse::<f32>().map_err(|_| invalido())?.to_le_bytes());
                }
                b'Z' | b'H' => {
                    s.push(b[3]);
                    s.extend_from_slice(v.as_bytes());
                    s.push(0);
                }
                b'B' => {
                    let mut partes = v.split(',');
                    let t = partes.next().filter(|t| t.len() == 1).ok_or_else(invalido)?.as_bytes()[0];
                    s.push(b'B');
                    s.push(t);
                    let onde = s.len();
                    s.extend_from_slice(&[0u8; 4]);
                    let mut n = 0i32;
                    for e in partes {
                        grava_elemento(t, e, s).map_err(|e| format!("optional field '{}': {e}", resumo(campo)))?;
                        n += 1;
                    }
                    s[onde..onde + 4].copy_from_slice(&n.to_le_bytes());
                }
                t => return Err(format!("optional field '{}' has the unknown type '{}'", resumo(campo), t as char)),
            }
        }
        if longo {
            s.extend_from_slice(b"CGBI");
            s.extend_from_slice(&(cigar.len() as i32).to_le_bytes());
            s.extend_from_slice(&cigar_bytes);
        }
        let tamanho = (s.len() - inicio - 4) as i32;
        s[inicio..inicio + 4].copy_from_slice(&tamanho.to_le_bytes());
        Ok(())
    }
}

/// Corta um texto longo para caber numa mensagem.
fn resumo(t: &str) -> String {
    if t.len() <= 60 {
        t.to_string()
    } else {
        format!("{}...", t.chars().take(57).collect::<String>())
    }
}

/// As sequencias declaradas nas linhas @SQ de um cabecalho de SAM.
pub fn refs_do_texto(texto: &str) -> Result<Vec<(String, u32)>, String> {
    let mut refs = Vec::new();
    for l in texto.lines().filter(|l| l.starts_with("@SQ\t")) {
        let nome = l.split('\t').find_map(|c| c.strip_prefix("SN:"));
        let tam = l.split('\t').find_map(|c| c.strip_prefix("LN:")).and_then(|t| t.parse::<u32>().ok());
        match (nome, tam) {
            (Some(n), Some(t)) => refs.push((n.to_string(), t)),
            _ => return Err(format!("header line '{}' needs both SN (name) and LN (length) to be written as BAM", resumo(l))),
        }
    }
    Ok(refs)
}

/// Grava o cabecalho binario: assinatura, texto e lista de sequencias.
pub fn grava_cabecalho<W: Write>(w: &mut W, cab: &Cabecalho) -> io::Result<()> {
    w.write_all(b"BAM\x01")?;
    w.write_all(&(cab.texto.len() as i32).to_le_bytes())?;
    w.write_all(cab.texto.as_bytes())?;
    w.write_all(&(cab.refs.len() as i32).to_le_bytes())?;
    for (nome, tam) in &cab.refs {
        w.write_all(&(nome.len() as i32 + 1).to_le_bytes())?;
        w.write_all(nome.as_bytes())?;
        w.write_all(&[0])?;
        w.write_all(&(*tam as i32).to_le_bytes())?;
    }
    Ok(())
}

/// Recebe SAM em texto e grava BAM. As linhas `@` do comeco formam o
/// cabecalho; o binario so comeca a sair quando chega o primeiro alinhamento
/// (ou no `termina`, se nao chegar nenhum).
pub struct BamEscritor<W: Write> {
    w: BgzfWriter<W>,
    /// Bytes recebidos que ainda nao fecharam uma linha.
    resto: Vec<u8>,
    cabecalho: String,
    cod: Option<Codificador>,
    linha_n: u64,
    buf: Vec<u8>,
}

impl<W: Write> BamEscritor<W> {
    pub fn novo(w: W) -> BamEscritor<W> {
        BamEscritor {
            w: BgzfWriter::novo(w, crate::bgzf::NIVEL_PADRAO),
            resto: Vec::new(),
            cabecalho: String::new(),
            cod: None,
            linha_n: 0,
            buf: Vec::new(),
        }
    }

    fn fecha_cabecalho(&mut self) -> io::Result<()> {
        let refs = refs_do_texto(&self.cabecalho).map_err(erro)?;
        grava_cabecalho(&mut self.w, &Cabecalho { texto: std::mem::take(&mut self.cabecalho), refs: refs.clone() })?;
        self.cod = Some(Codificador::novo(&refs));
        Ok(())
    }

    fn linha(&mut self, linha: &[u8]) -> io::Result<()> {
        self.linha_n += 1;
        let linha = match linha.last() {
            Some(b'\r') => &linha[..linha.len() - 1],
            _ => linha,
        };
        if linha.is_empty() {
            return Ok(());
        }
        let n = self.linha_n;
        let texto = std::str::from_utf8(linha).map_err(|_| erro(format!("line {n}: not valid text (UTF-8)")))?;
        if self.cod.is_none() {
            if texto.starts_with('@') {
                self.cabecalho.push_str(texto);
                self.cabecalho.push('\n');
                return Ok(());
            }
            self.fecha_cabecalho()?;
        }
        if texto.starts_with('@') {
            return Err(erro(format!("line {n}: a header line after the first alignment cannot be written to BAM; put every @ line at the top")));
        }
        self.buf.clear();
        match self.cod.as_ref() {
            Some(c) => c.codifica(texto, &mut self.buf).map_err(|e| erro(format!("line {n}: {e}")))?,
            None => return Err(erro("BAM header was not written".into())),
        }
        self.w.write_all(&self.buf)
    }

    /// Grava o que falta e fecha o BGZF.
    pub fn termina(mut self) -> io::Result<W> {
        if !self.resto.is_empty() {
            let resto = std::mem::take(&mut self.resto);
            self.linha(&resto)?;
        }
        if self.cod.is_none() {
            self.fecha_cabecalho()?;
        }
        self.w.termina()
    }
}

impl<W: Write> Write for BamEscritor<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let mut inicio = 0;
        while let Some(p) = buf[inicio..].iter().position(|&b| b == b'\n') {
            if self.resto.is_empty() {
                self.linha(&buf[inicio..inicio + p])?;
            } else {
                self.resto.extend_from_slice(&buf[inicio..inicio + p]);
                let completa = std::mem::take(&mut self.resto);
                self.linha(&completa)?;
            }
            inicio += p + 1;
        }
        self.resto.extend_from_slice(&buf[inicio..]);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// O texto do cabecalho como um SAM o teria: se o binario declara sequencias
/// e o texto nao traz nenhum @SQ, eles sao gerados (depois do @HD).
pub fn texto_do_cabecalho(cab: &Cabecalho) -> String {
    let mut linhas: Vec<String> = cab.texto.lines().filter(|l| !l.is_empty()).map(str::to_string).collect();
    if !cab.refs.is_empty() && !linhas.iter().any(|l| l.starts_with("@SQ")) {
        let onde = usize::from(linhas.first().is_some_and(|l| l.starts_with("@HD")));
        for (i, (nome, tam)) in cab.refs.iter().enumerate() {
            linhas.insert(onde + i, format!("@SQ\tSN:{nome}\tLN:{tam}"));
        }
    }
    let mut t = linhas.join("\n");
    if !t.is_empty() {
        t.push('\n');
    }
    t
}

/// Le BAM e entrega SAM em texto: o cabecalho e depois uma linha por
/// alinhamento.
pub struct BamLeitor<R: Read> {
    r: BgzfParLeitor<R>,
    refs: Vec<(String, u32)>,
    /// Texto ja convertido e ainda nao entregue.
    texto: Vec<u8>,
    pos: usize,
    rec: Vec<u8>,
    n: u64,
    acabou: bool,
}

impl<R: Read> BamLeitor<R> {
    pub fn novo(r: R) -> io::Result<BamLeitor<R>> {
        let mut r = BgzfParLeitor::novo(r);
        let cab = le_cabecalho(&mut r)?;
        Ok(BamLeitor {
            texto: texto_do_cabecalho(&cab).into_bytes(),
            r,
            refs: cab.refs,
            pos: 0,
            rec: Vec::new(),
            n: 0,
            acabou: false,
        })
    }

    /// Converte mais um punhado de alinhamentos. `false` quando nao ha mais.
    fn enche(&mut self) -> io::Result<bool> {
        self.texto.clear();
        self.pos = 0;
        while !self.acabou && self.texto.len() < (256 << 10) {
            if !le_registro(&mut self.r, &mut self.rec)? {
                self.acabou = true;
                break;
            }
            self.n += 1;
            let n = self.n;
            decodifica(&self.rec, &self.refs, &mut self.texto).map_err(|e| erro(format!("alignment {n}: {e}")))?;
            self.texto.push(b'\n');
        }
        Ok(!self.texto.is_empty())
    }
}

impl<R: Read> Read for BamLeitor<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.pos >= self.texto.len() && !self.enche()? {
            return Ok(0);
        }
        let n = buf.len().min(self.texto.len() - self.pos);
        buf[..n].copy_from_slice(&self.texto[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}
