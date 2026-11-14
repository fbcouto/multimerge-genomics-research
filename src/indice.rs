// =====================================================================
//  INDICES: tabix (.tbi) para VCF, BED e GFF em BGZF, e BAI (.bai) para BAM
//
//  Os dois indices tem a mesma estrutura -- o esquema de bins do UCSC mais um
//  indice linear de janelas de 16 kb -- e respondem a mesma pergunta: em que
//  pedacos do arquivo comprimido podem estar os registros que tocam uma
//  regiao. A consulta descomprime so esses pedacos.
//
//  A construcao segue o htslib passo a passo (hts_idx_push / hts_idx_finish),
//  para que o indice sirva ao tabix, ao samtools, ao bcftools e ao IGV:
//
//   * cada registro entra no menor bin que o contem; registros consecutivos
//     do mesmo bin formam um pedaco (offset virtual de inicio e de fim);
//   * um bin cujos pedacos cabem em 64 KiB de arquivo e fundido no bin pai;
//   * pedacos que comecam no bloco em que o anterior termina sao unidos;
//   * o indice linear guarda, por janela, onde comeca o primeiro registro que
//     a toca.
//
//  O arquivo precisa estar ordenado por posicao, com cada sequencia num bloco
//  so. O limite do formato e 2^29 bases por sequencia (512 Mbp).
// =====================================================================

use crate::bam;
use crate::bgzf::{self, BgzfReader, BgzfWriter};
use crate::externo::Formato;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};

/// Pseudo-bin do htslib: nao indexa nada, guarda estatisticas da sequencia.
const BIN_META: u32 = 37450;
const SEM_BIN: u32 = u32::MAX;
const SEM_OFFSET: u64 = u64::MAX;
/// Maior coordenada que um .tbi ou .bai representa.
pub const POS_MAX: i64 = 1 << 29;
/// Primeiro bin de cada nivel, do mais largo (0) ao mais fino (5).
const PRIMEIRO_DO_NIVEL: [u32; 7] = [0, 1, 9, 73, 585, 4681, 37449];

type Pedaco = (u64, u64);

/// O indice de uma sequencia.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IndiceRef {
    pub bins: BTreeMap<u32, Vec<Pedaco>>,
    pub linear: Vec<u64>,
    /// Inicio e fim, no arquivo, do bloco de registros da sequencia, e quantos
    /// sao mapeados e nao mapeados (o pseudo-bin 37450).
    pub meta: Option<[u64; 4]>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Indice {
    pub refs: Vec<IndiceRef>,
    /// Registros sem coordenada (no BAM: nao mapeados, no fim do arquivo).
    pub sem_coordenada: u64,
}

/// A configuracao de um .tbi: que colunas trazem sequencia, inicio e fim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tabix {
    /// 0 generico, 2 VCF; `0x10000` somado quando o inicio e base 0 (BED).
    pub formato: i32,
    pub col_seq: i32,
    pub col_ini: i32,
    pub col_fim: i32,
    pub meta: u8,
    pub pular: i32,
}

impl Tabix {
    pub fn de(formato: Formato) -> Result<Tabix, String> {
        match formato {
            Formato::Vcf => Ok(Tabix { formato: 2, col_seq: 1, col_ini: 2, col_fim: 0, meta: b'#', pular: 0 }),
            Formato::Bed => Ok(Tabix { formato: 0x10000, col_seq: 1, col_ini: 2, col_fim: 3, meta: b'#', pular: 0 }),
            Formato::Gff => Ok(Tabix { formato: 0, col_seq: 1, col_ini: 4, col_fim: 5, meta: b'#', pular: 0 }),
            Formato::Sam => Err("SAM text is not indexed here: write it as BAM (a .bam output) and index the BAM".into()),
        }
    }

    /// Sequencia e intervalo `[inicio, fim)` em base 0 de uma linha de dados.
    /// `None`: a linha nao e um registro valido.
    pub fn intervalo<'a>(&self, linha: &'a str) -> Option<(&'a str, i64, i64)> {
        let (mut seq, mut ini, mut fim) = (None, -1i64, -1i64);
        for (i, campo) in linha.split('\t').enumerate() {
            let col = i as i32 + 1;
            if col == self.col_seq {
                seq = Some(campo);
            } else if col == self.col_ini {
                ini = campo.trim_end().parse().ok()?;
                if self.col_ini <= self.col_fim || self.col_fim == 0 {
                    fim = ini;
                }
                if self.formato & 0x10000 == 0 {
                    ini -= 1;
                } else {
                    fim += 1;
                }
                if ini < 0 {
                    ini = 0;
                }
                if fim < 1 {
                    fim = 1;
                }
            } else if self.formato & 0xffff == 0 {
                if col == self.col_fim {
                    fim = campo.trim_end().parse().ok()?;
                }
            } else if self.formato & 0xffff == 2 {
                if col == 4 {
                    if !campo.is_empty() && ini >= 0 {
                        fim = ini + campo.len() as i64;
                    }
                } else if col == 8 {
                    // INFO/END estende a variante (delecoes simbolicas, gVCF).
                    let end = campo.split(';').find_map(|c| c.strip_prefix("END="));
                    if let Some(e) = end.and_then(|e| e.parse::<i64>().ok()) {
                        if e > ini {
                            fim = e;
                        }
                    }
                }
            }
        }
        let seq = seq?;
        if seq.is_empty() || ini < 0 || fim < 0 {
            return None;
        }
        Some((seq, ini, fim.max(ini)))
    }
}

// ---------------------------------------------------------------------
//  CONSTRUCAO
// ---------------------------------------------------------------------

pub struct Construtor {
    idx: Indice,
    visto: Vec<bool>,
    ultimo_tid: i32,
    ultimo_bin: u32,
    guardado_tid: i32,
    guardado_bin: u32,
    guardado_off: u64,
    ultimo_off: u64,
    off_ini: u64,
    n_mapeados: u64,
    n_nao_mapeados: u64,
    ultima_coord: i64,
}

fn nivel_do_bin(bin: u32) -> usize {
    (0..6).rev().find(|&l| bin >= PRIMEIRO_DO_NIVEL[l]).unwrap_or(0)
}

impl Construtor {
    /// `offset0`: offset virtual do primeiro registro (logo depois do
    /// cabecalho).
    pub fn novo(n_refs: usize, offset0: u64) -> Construtor {
        Construtor {
            idx: Indice { refs: vec![IndiceRef::default(); n_refs], sem_coordenada: 0 },
            visto: vec![false; n_refs],
            ultimo_tid: -1,
            ultimo_bin: SEM_BIN,
            guardado_tid: -1,
            guardado_bin: SEM_BIN,
            guardado_off: offset0,
            ultimo_off: offset0,
            off_ini: offset0,
            n_mapeados: 0,
            n_nao_mapeados: 0,
            ultima_coord: 0,
        }
    }

    fn poe_no_bin(&mut self, tid: i32, bin: u32, ini: u64, fim: u64) {
        if tid >= 0 {
            self.idx.refs[tid as usize].bins.entry(bin).or_default().push((ini, fim));
        }
    }

    /// Registra um registro: sequencia `tid` (negativa: sem coordenada),
    /// intervalo `[ini, fim)` em base 0 e o offset virtual logo DEPOIS dele.
    pub fn empurra(&mut self, tid: i32, mut ini: i64, mut fim: i64, offset: u64, mapeado: bool) -> Result<(), String> {
        if tid < 0 {
            ini = -1;
            fim = 0;
        } else if ini > POS_MAX || fim > POS_MAX {
            return Err(format!(
                "a record at position {} is beyond the {POS_MAX} bases (512 Mbp) per sequence that this kind of index can address",
                ini + 1
            ));
        }
        if tid >= 0 && tid as usize >= self.idx.refs.len() {
            self.idx.refs.resize(tid as usize + 1, IndiceRef::default());
            self.visto.resize(tid as usize + 1, false);
        }
        if self.ultimo_tid != tid {
            if tid >= 0 && self.idx.sem_coordenada > 0 {
                return Err("records with a position come after records without one: the file is not sorted by position".into());
            }
            if tid >= 0 && self.visto[tid as usize] {
                return Err("the records of one sequence are split in more than one stretch of the file: it is not sorted by position".into());
            }
            self.ultimo_tid = tid;
            self.ultimo_bin = SEM_BIN;
        } else if tid >= 0 && self.ultima_coord > ini {
            return Err(format!("position {} comes after position {} on the same sequence: the file is not sorted by position", ini + 1, self.ultima_coord + 1));
        }
        if fim < ini {
            return Err(format!("a record ends (at {fim}) before it starts (at {})", ini + 1));
        }
        if tid >= 0 {
            self.visto[tid as usize] = true;
            ini = ini.max(0);
            if fim <= 0 {
                fim = 1;
            }
            // Indice linear: o primeiro registro que toca cada janela de 16 kb.
            let linear = &mut self.idx.refs[tid as usize].linear;
            let (j0, j1) = ((ini >> 14) as usize, ((fim - 1) >> 14) as usize);
            if linear.len() < j1 + 1 {
                linear.resize(j1 + 1, SEM_OFFSET);
            }
            for v in &mut linear[j0..=j1] {
                if *v == SEM_OFFSET {
                    *v = self.ultimo_off;
                }
            }
        } else {
            self.idx.sem_coordenada += 1;
        }
        let bin = bam::reg2bin(ini, fim);
        if self.ultimo_bin != bin {
            if self.guardado_bin != SEM_BIN {
                self.poe_no_bin(self.guardado_tid, self.guardado_bin, self.guardado_off, self.ultimo_off);
            }
            if self.ultimo_bin == SEM_BIN && self.guardado_bin != SEM_BIN {
                // Mudou de sequencia: fecha as estatisticas da anterior.
                if self.guardado_tid >= 0 {
                    self.idx.refs[self.guardado_tid as usize].meta =
                        Some([self.off_ini, self.ultimo_off, self.n_mapeados, self.n_nao_mapeados]);
                }
                self.n_mapeados = 0;
                self.n_nao_mapeados = 0;
                self.off_ini = self.ultimo_off;
            }
            self.guardado_off = self.ultimo_off;
            self.guardado_bin = bin;
            self.ultimo_bin = bin;
            self.guardado_tid = tid;
        }
        if mapeado {
            self.n_mapeados += 1;
        } else {
            self.n_nao_mapeados += 1;
        }
        self.ultimo_off = offset;
        self.ultima_coord = ini;
        Ok(())
    }

    /// Fecha o indice. `offset_final`: offset virtual do fim dos registros.
    pub fn termina(mut self, offset_final: u64) -> Indice {
        if self.guardado_tid >= 0 {
            self.poe_no_bin(self.guardado_tid, self.guardado_bin, self.guardado_off, offset_final);
            self.idx.refs[self.guardado_tid as usize].meta =
                Some([self.off_ini, offset_final, self.n_mapeados, self.n_nao_mapeados]);
        }
        for r in &mut self.idx.refs {
            // Janela sem registro proprio herda o offset da janela seguinte.
            for i in (0..r.linear.len().saturating_sub(1)).rev() {
                if r.linear[i] == SEM_OFFSET {
                    r.linear[i] = r.linear[i + 1];
                }
            }
            // Bin pequeno demais (menos de 64 KiB de arquivo) sobe para o pai.
            for nivel in (1..=5usize).rev() {
                let bins: Vec<u32> =
                    r.bins.range(PRIMEIRO_DO_NIVEL[nivel]..PRIMEIRO_DO_NIVEL[nivel + 1]).map(|(&b, _)| b).collect();
                for bin in bins {
                    let Some(mut pedacos) = r.bins.remove(&bin) else { continue };
                    pedacos.sort_by_key(|p| p.0);
                    let largura = (pedacos[pedacos.len() - 1].1 >> 16) - (pedacos[0].0 >> 16);
                    if largura < 0x10000 {
                        r.bins.entry((bin - 1) >> 3).or_default().extend(pedacos);
                    } else {
                        r.bins.insert(bin, pedacos);
                    }
                }
            }
            debug_assert!(r.bins.keys().all(|&b| nivel_do_bin(b) <= 5));
            // Une os pedacos que comecam no bloco em que o anterior termina.
            for pedacos in r.bins.values_mut() {
                pedacos.sort_by_key(|p| p.0);
                let mut m = 0;
                for l in 1..pedacos.len() {
                    if pedacos[m].1 >> 16 >= pedacos[l].0 >> 16 {
                        if pedacos[m].1 < pedacos[l].1 {
                            pedacos[m].1 = pedacos[l].1;
                        }
                    } else {
                        m += 1;
                        pedacos[m] = pedacos[l];
                    }
                }
                pedacos.truncate(m + 1);
            }
        }
        self.idx
    }
}

// ---------------------------------------------------------------------
//  GRAVAR E LER
// ---------------------------------------------------------------------

fn serializa_refs(idx: &Indice, s: &mut Vec<u8>) {
    for r in &idx.refs {
        let n_bins = r.bins.len() + usize::from(r.meta.is_some());
        s.extend_from_slice(&(n_bins as i32).to_le_bytes());
        for (&bin, pedacos) in &r.bins {
            s.extend_from_slice(&bin.to_le_bytes());
            s.extend_from_slice(&(pedacos.len() as i32).to_le_bytes());
            for p in pedacos {
                s.extend_from_slice(&p.0.to_le_bytes());
                s.extend_from_slice(&p.1.to_le_bytes());
            }
        }
        if let Some(m) = r.meta {
            s.extend_from_slice(&BIN_META.to_le_bytes());
            s.extend_from_slice(&2i32.to_le_bytes());
            for v in m {
                s.extend_from_slice(&v.to_le_bytes());
            }
        }
        s.extend_from_slice(&(r.linear.len() as i32).to_le_bytes());
        for v in &r.linear {
            s.extend_from_slice(&v.to_le_bytes());
        }
    }
    s.extend_from_slice(&idx.sem_coordenada.to_le_bytes());
}

pub fn grava_bai(idx: &Indice, caminho: &str) -> Result<(), String> {
    let mut s = Vec::new();
    s.extend_from_slice(b"BAI\x01");
    s.extend_from_slice(&(idx.refs.len() as i32).to_le_bytes());
    serializa_refs(idx, &mut s);
    std::fs::write(caminho, s).map_err(|e| format!("Error writing {caminho}: {e}"))
}

pub fn grava_tbi(idx: &Indice, conf: &Tabix, nomes: &[String], caminho: &str) -> Result<(), String> {
    let mut s = Vec::new();
    s.extend_from_slice(b"TBI\x01");
    s.extend_from_slice(&(idx.refs.len() as i32).to_le_bytes());
    for v in [conf.formato, conf.col_seq, conf.col_ini, conf.col_fim, i32::from(conf.meta), conf.pular] {
        s.extend_from_slice(&v.to_le_bytes());
    }
    let l_nomes: usize = nomes.iter().map(|n| n.len() + 1).sum();
    s.extend_from_slice(&(l_nomes as i32).to_le_bytes());
    for n in nomes {
        s.extend_from_slice(n.as_bytes());
        s.push(0);
    }
    serializa_refs(idx, &mut s);
    let arq = File::create(caminho).map_err(|e| format!("Error creating {caminho}: {e}"))?;
    let mut w = BgzfWriter::novo(arq, bgzf::NIVEL_PADRAO);
    w.write_all(&s).and_then(|_| w.termina().map(|_| ())).map_err(|e| format!("Error writing {caminho}: {e}"))
}

/// Cursor de leitura sobre o indice ja em memoria.
struct Cursor<'a> {
    b: &'a [u8],
    p: usize,
    caminho: &'a str,
}

impl Cursor<'_> {
    fn bytes(&mut self, n: usize) -> Result<&[u8], String> {
        if self.p + n > self.b.len() {
            return Err(format!("{}: the index is truncated or damaged; build it again", self.caminho));
        }
        self.p += n;
        Ok(&self.b[self.p - n..self.p])
    }
    fn i32(&mut self) -> Result<i32, String> {
        let b = self.bytes(4)?;
        Ok(i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }
    fn u32(&mut self) -> Result<u32, String> {
        Ok(self.i32()? as u32)
    }
    fn u64(&mut self) -> Result<u64, String> {
        let b = self.bytes(8)?;
        Ok(u64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]))
    }
    fn tamanho(&mut self) -> Result<usize, String> {
        let n = self.i32()?;
        if n < 0 {
            return Err(format!("{}: the index is damaged; build it again", self.caminho));
        }
        Ok(n as usize)
    }

    fn refs(&mut self, n_refs: usize) -> Result<Indice, String> {
        let mut idx = Indice::default();
        for _ in 0..n_refs {
            let mut r = IndiceRef::default();
            for _ in 0..self.tamanho()? {
                let bin = self.u32()?;
                let n = self.tamanho()?;
                let mut pedacos = Vec::with_capacity(n.min(1 << 16));
                for _ in 0..n {
                    pedacos.push((self.u64()?, self.u64()?));
                }
                if bin == BIN_META {
                    if pedacos.len() == 2 {
                        r.meta = Some([pedacos[0].0, pedacos[0].1, pedacos[1].0, pedacos[1].1]);
                    }
                } else {
                    r.bins.insert(bin, pedacos);
                }
            }
            let n = self.tamanho()?;
            r.linear.reserve(n.min(1 << 16));
            for _ in 0..n {
                r.linear.push(self.u64()?);
            }
            idx.refs.push(r);
        }
        // O total de registros sem coordenada e opcional no fim do arquivo.
        if self.p + 8 <= self.b.len() {
            idx.sem_coordenada = self.u64()?;
        }
        Ok(idx)
    }
}

pub fn le_bai(caminho: &str) -> Result<Indice, String> {
    let b = std::fs::read(caminho).map_err(|e| format!("Error opening {caminho}: {e}"))?;
    let mut c = Cursor { b: &b, p: 0, caminho };
    if c.bytes(4)? != b"BAI\x01" {
        return Err(format!("{caminho}: not a BAI index"));
    }
    let n = c.tamanho()?;
    c.refs(n)
}

pub fn le_tbi(caminho: &str) -> Result<(Indice, Tabix, Vec<String>), String> {
    let arq = File::open(caminho).map_err(|e| format!("Error opening {caminho}: {e}"))?;
    let mut b = Vec::new();
    BgzfReader::novo(arq).read_to_end(&mut b).map_err(|e| format!("{caminho}: {e}"))?;
    let mut c = Cursor { b: &b, p: 0, caminho };
    if c.bytes(4)? != b"TBI\x01" {
        return Err(format!("{caminho}: not a tabix index"));
    }
    let n = c.tamanho()?;
    let conf = Tabix {
        formato: c.i32()?,
        col_seq: c.i32()?,
        col_ini: c.i32()?,
        col_fim: c.i32()?,
        meta: c.i32()? as u8,
        pular: c.i32()?,
    };
    let l_nomes = c.tamanho()?;
    let nomes: Vec<String> = c
        .bytes(l_nomes)?
        .split(|&x| x == 0)
        .filter(|n| !n.is_empty())
        .map(|n| String::from_utf8_lossy(n).to_string())
        .collect();
    let idx = c.refs(n)?;
    Ok((idx, conf, nomes))
}

// ---------------------------------------------------------------------
//  INDEXAR UM ARQUIVO
// ---------------------------------------------------------------------

/// O que a indexacao fez.
#[derive(Clone, Debug, Default)]
pub struct ResumoIndice {
    /// Caminho do indice gravado.
    pub indice: String,
    /// "bai" ou "tbi".
    pub tipo: &'static str,
    pub sequencias: usize,
    pub registros: u64,
    /// Linhas que nao eram cabecalho nem registro valido.
    pub ignoradas: u64,
}

fn e_bam(caminho: &str) -> bool {
    caminho.to_ascii_lowercase().ends_with(".bam")
}

fn abre_bgzf(caminho: &str) -> Result<BgzfReader<BufReader<File>>, String> {
    let arq = File::open(caminho).map_err(|e| format!("Error opening {caminho}: {e}"))?;
    Ok(BgzfReader::novo(BufReader::with_capacity(1 << 20, arq)))
}

const COMO_ORDENAR: &str = "sort it by position first (sort_genomic with by=\"position\") and index the sorted file";

fn indexa_bam(caminho: &str) -> Result<ResumoIndice, String> {
    let mut r = abre_bgzf(caminho)?;
    let cab = bam::le_cabecalho(&mut r).map_err(|e| format!("{caminho}: {e}"))?;
    let mut c = Construtor::novo(cab.refs.len(), r.pos_virtual());
    let mut rec = Vec::new();
    let mut n = 0u64;
    while bam::le_registro(&mut r, &mut rec).map_err(|e| format!("{caminho}: {e}"))? {
        n += 1;
        let i = bam::info(&rec).map_err(|e| format!("{caminho}, alignment {n}: {e}"))?;
        if i.tid >= cab.refs.len() as i32 {
            return Err(format!("{caminho}, alignment {n}: reference sequence number {} is not in the header", i.tid));
        }
        c.empurra(i.tid, i.pos, i.fim, r.pos_virtual(), i.flag & 4 == 0)
            .map_err(|e| format!("{caminho}, alignment {n}: {e}; {COMO_ORDENAR}"))?;
    }
    let idx = c.termina(r.pos_virtual());
    let indice = format!("{caminho}.bai");
    grava_bai(&idx, &indice)?;
    Ok(ResumoIndice { indice, tipo: "bai", sequencias: idx.refs.iter().filter(|r| r.meta.is_some()).count(), registros: n, ignoradas: 0 })
}

/// Le uma linha (sem o fim de linha) para `linha`. `false` no fim do arquivo.
fn le_linha<R: BufRead>(r: &mut R, bruto: &mut Vec<u8>, caminho: &str) -> Result<bool, String> {
    bruto.clear();
    if r.read_until(b'\n', bruto).map_err(|e| format!("{caminho}: {e}"))? == 0 {
        return Ok(false);
    }
    while matches!(bruto.last(), Some(b'\n' | b'\r')) {
        bruto.pop();
    }
    Ok(true)
}

fn e_meta(conf: &Tabix, linha: &[u8]) -> bool {
    linha.first() == Some(&conf.meta)
        || (conf.formato & 0x10000 != 0 && (linha.starts_with(b"track") || linha.starts_with(b"browser")))
}

fn indexa_texto(caminho: &str, formato: Formato) -> Result<ResumoIndice, String> {
    let conf = Tabix::de(formato)?;
    if !bgzf::e_bgzf(caminho) {
        return Err(format!(
            "{caminho} is not BGZF-compressed, so it cannot be indexed: write it again with a .gz output from this library (which is BGZF), or compress it with bgzip"
        ));
    }
    let mut r = abre_bgzf(caminho)?;
    let mut nomes: Vec<String> = Vec::new();
    let mut ids: HashMap<String, i32> = HashMap::new();
    let mut c: Option<Construtor> = None;
    let mut offset0 = 0u64;
    let (mut n, mut ignoradas, mut n_linha) = (0u64, 0u64, 0u64);
    let mut bruto = Vec::new();
    while le_linha(&mut r, &mut bruto, caminho)? {
        n_linha += 1;
        let depois = r.pos_virtual();
        if bruto.is_empty() || e_meta(&conf, &bruto) {
            if c.is_none() {
                offset0 = depois;
            }
            continue;
        }
        let c = c.get_or_insert_with(|| Construtor::novo(0, offset0));
        let Some((seq, ini, fim)) = std::str::from_utf8(&bruto).ok().and_then(|l| conf.intervalo(l)) else {
            ignoradas += 1;
            continue;
        };
        let tid = match ids.get(seq) {
            Some(&t) => t,
            None => {
                let t = nomes.len() as i32;
                nomes.push(seq.to_string());
                ids.insert(seq.to_string(), t);
                t
            }
        };
        c.empurra(tid, ini, fim, depois, true).map_err(|e| format!("{caminho}, line {n_linha}: {e}; {COMO_ORDENAR}"))?;
        n += 1;
    }
    let idx = c.unwrap_or_else(|| Construtor::novo(0, offset0)).termina(r.pos_virtual());
    let indice = format!("{caminho}.tbi");
    grava_tbi(&idx, &conf, &nomes, &indice)?;
    Ok(ResumoIndice { indice, tipo: "tbi", sequencias: nomes.len(), registros: n, ignoradas })
}

/// Constroi o indice de `caminho`: `.bai` para um BAM, `.tbi` para VCF, BED ou
/// GFF comprimido em BGZF. O indice fica ao lado do arquivo, com a extensao
/// acrescentada (`x.vcf.gz.tbi`, `x.bam.bai`).
pub fn indexar(caminho: &str, formato: Option<Formato>) -> Result<ResumoIndice, String> {
    if e_bam(caminho) {
        return indexa_bam(caminho);
    }
    let formato = formato.or_else(|| Formato::pelo_nome(caminho)).ok_or_else(|| {
        format!("Cannot tell the format of {caminho} from its name; pass format=\"vcf\", \"bed\" or \"gff\"")
    })?;
    indexa_texto(caminho, formato)
}

// ---------------------------------------------------------------------
//  CONSULTA POR REGIAO
// ---------------------------------------------------------------------

/// Os bins que podem conter registros que tocam `[ini, fim)`.
fn bins_da_regiao(ini: i64, fim: i64) -> Vec<u32> {
    let fim = fim - 1;
    let mut v = vec![0u32];
    for (desloca, primeiro) in [(26, 1i64), (23, 9), (20, 73), (17, 585), (14, 4681)] {
        v.extend(((primeiro + (ini >> desloca))..=(primeiro + (fim >> desloca))).map(|b| b as u32));
    }
    v
}

impl IndiceRef {
    /// Os pedacos do arquivo a ler para a regiao `[ini, fim)`, em ordem e sem
    /// sobreposicao.
    pub fn pedacos(&self, ini: i64, fim: i64) -> Vec<Pedaco> {
        // Nada que toque a regiao comeca antes do primeiro registro que toca
        // a janela onde ela comeca.
        let janela = (ini >> 14) as usize;
        let minimo = self.linear.get(janela).or(self.linear.last()).copied().unwrap_or(0);
        let mut ps: Vec<Pedaco> = bins_da_regiao(ini, fim)
            .iter()
            .filter_map(|b| self.bins.get(b))
            .flatten()
            .filter(|p| p.1 > minimo)
            .map(|p| (p.0.max(minimo), p.1))
            .collect();
        ps.sort();
        let mut unidos: Vec<Pedaco> = Vec::with_capacity(ps.len());
        for p in ps {
            match unidos.last_mut() {
                Some(u) if p.0 <= u.1 => u.1 = u.1.max(p.1),
                _ => unidos.push(p),
            }
        }
        unidos
    }
}

/// Le uma regiao como o samtools e o tabix: `chr1`, `chr1:1000` ou
/// `chr1:1000-2000` (base 1, fim incluido; virgulas nos numeros sao aceitas).
/// Devolve `(nome, inicio base 0, fim exclusivo)`. `conhece` diz se um nome
/// existe no arquivo: e o que decide se um `:` faz parte do nome.
pub fn le_regiao(texto: &str, conhece: &dyn Fn(&str) -> bool) -> Result<(String, i64, i64), String> {
    let t = texto.trim();
    if t.is_empty() {
        return Err("Empty region. Use \"chr1\", \"chr1:1000\" or \"chr1:1000-2000\"".into());
    }
    if conhece(t) {
        return Ok((t.to_string(), 0, POS_MAX));
    }
    let Some((nome, faixa)) = t.rsplit_once(':') else {
        return Ok((t.to_string(), 0, POS_MAX));
    };
    let errada = || format!("Region \"{texto}\" is not valid. Use \"chr1\", \"chr1:1000\" or \"chr1:1000-2000\" (1-based, end included)");
    let numero = |s: &str| s.replace(',', "").trim().parse::<i64>().map_err(|_| errada());
    let (ini, fim) = match faixa.split_once('-') {
        Some((a, b)) if b.trim().is_empty() => (numero(a)?, POS_MAX),
        Some((a, b)) => (numero(a)?, numero(b)?),
        None => {
            let p = numero(faixa)?;
            (p, p)
        }
    };
    if nome.is_empty() || ini < 1 || fim < ini {
        return Err(errada());
    }
    Ok((nome.to_string(), ini - 1, fim.min(POS_MAX)))
}

/// Junta as regioes por sequencia, na ordem do arquivo, e une as que se tocam:
/// cada registro sai uma vez so, e a saida fica na ordem do arquivo.
fn agrupa_regioes(regioes: &[(i32, i64, i64)]) -> Vec<(i32, i64, i64)> {
    let mut r = regioes.to_vec();
    r.sort();
    let mut unidas: Vec<(i32, i64, i64)> = Vec::with_capacity(r.len());
    for x in r {
        match unidas.last_mut() {
            Some(u) if u.0 == x.0 && x.1 <= u.2 => u.2 = u.2.max(x.2),
            _ => unidas.push(x),
        }
    }
    unidas
}

fn acha_indice(caminho: &str, extensao: &str) -> Result<String, String> {
    let ao_lado = format!("{caminho}.{extensao}");
    if std::path::Path::new(&ao_lado).exists() {
        return Ok(ao_lado);
    }
    // O samtools tambem aceita `x.bai` para `x.bam`.
    if extensao == "bai" {
        let trocado = format!("{}.bai", &caminho[..caminho.len() - 4]);
        if std::path::Path::new(&trocado).exists() {
            return Ok(trocado);
        }
    }
    Err(format!("No index found for {caminho} (expected {ao_lado}). Build it with index_file(\"{caminho}\")"))
}

fn sequencia_ausente(nome: &str, nomes: &[&str], caminho: &str) -> String {
    let exemplos: Vec<&str> = nomes.iter().take(5).copied().collect();
    format!(
        "Sequence \"{nome}\" is not in {caminho}. Its sequences are named like: {}{}. If the names differ only by the \"chr\" prefix, use normalize_chromosomes first",
        exemplos.join(", "),
        if nomes.len() > 5 { ", ..." } else { "" }
    )
}

/// Entrega, em ordem, os registros de `caminho` que tocam alguma das regioes.
/// `recebe(true, linha)` para cada linha de cabecalho, `recebe(false, linha)`
/// para cada registro. Devolve quantos registros entregou.
///
/// O arquivo precisa do indice ao lado (`indexar`). So os pedacos apontados
/// pelo indice sao descomprimidos.
pub fn consultar(
    caminho: &str,
    regioes: &[String],
    recebe: &mut dyn FnMut(bool, &str) -> Result<(), String>,
) -> Result<u64, String> {
    if regioes.is_empty() {
        return Err("No region given. Use \"chr1\", \"chr1:1000\" or \"chr1:1000-2000\"".into());
    }
    if e_bam(caminho) {
        consulta_bam(caminho, regioes, recebe)
    } else {
        consulta_texto(caminho, regioes, recebe)
    }
}

fn resolve_regioes(regioes: &[String], nomes: &[&str], caminho: &str) -> Result<Vec<(i32, i64, i64)>, String> {
    let ids: HashMap<&str, i32> = nomes.iter().enumerate().map(|(i, n)| (*n, i as i32)).collect();
    let mut v = Vec::with_capacity(regioes.len());
    for texto in regioes {
        let (nome, ini, fim) = le_regiao(texto, &|n| ids.contains_key(n))?;
        let tid = *ids.get(nome.as_str()).ok_or_else(|| sequencia_ausente(&nome, nomes, caminho))?;
        v.push((tid, ini, fim));
    }
    Ok(agrupa_regioes(&v))
}

fn consulta_texto(
    caminho: &str,
    regioes: &[String],
    recebe: &mut dyn FnMut(bool, &str) -> Result<(), String>,
) -> Result<u64, String> {
    let (idx, conf, nomes) = le_tbi(&acha_indice(caminho, "tbi")?)?;
    let nomes_ref: Vec<&str> = nomes.iter().map(String::as_str).collect();
    let alvo = resolve_regioes(regioes, &nomes_ref, caminho)?;
    let mut r = abre_bgzf(caminho)?;
    let mut bruto = Vec::new();

    // O cabecalho: as linhas de comentario do topo.
    while le_linha(&mut r, &mut bruto, caminho)? {
        if !bruto.is_empty() && !e_meta(&conf, &bruto) {
            break;
        }
        if !bruto.is_empty() {
            recebe(true, &String::from_utf8_lossy(&bruto))?;
        }
    }

    let mut total = 0u64;
    for (tid, ini, fim) in alvo {
        let nome = nomes_ref[tid as usize];
        let Some(iref) = idx.refs.get(tid as usize) else { continue };
        'regiao: for (p0, p1) in iref.pedacos(ini, fim) {
            r.busca(p0).map_err(|e| format!("{caminho}: {e}"))?;
            while r.pos_virtual() < p1 && le_linha(&mut r, &mut bruto, caminho)? {
                if bruto.is_empty() || e_meta(&conf, &bruto) {
                    continue;
                }
                let Ok(linha) = std::str::from_utf8(&bruto) else { continue };
                let Some((seq, a, b)) = conf.intervalo(linha) else { continue };
                if seq != nome {
                    continue;
                }
                if a >= fim {
                    break 'regiao;
                }
                // Um registro de comprimento zero conta como uma base.
                if b.max(a + 1) > ini {
                    recebe(false, linha)?;
                    total += 1;
                }
            }
        }
    }
    Ok(total)
}

fn consulta_bam(
    caminho: &str,
    regioes: &[String],
    recebe: &mut dyn FnMut(bool, &str) -> Result<(), String>,
) -> Result<u64, String> {
    let idx = le_bai(&acha_indice(caminho, "bai")?)?;
    let mut r = abre_bgzf(caminho)?;
    let cab = bam::le_cabecalho(&mut r).map_err(|e| format!("{caminho}: {e}"))?;
    let nomes_ref: Vec<&str> = cab.refs.iter().map(|x| x.0.as_str()).collect();
    let alvo = resolve_regioes(regioes, &nomes_ref, caminho)?;
    for l in bam::texto_do_cabecalho(&cab).lines() {
        recebe(true, l)?;
    }
    let (mut rec, mut texto) = (Vec::new(), Vec::new());
    let mut total = 0u64;
    for (tid, ini, fim) in alvo {
        let Some(iref) = idx.refs.get(tid as usize) else { continue };
        'regiao: for (p0, p1) in iref.pedacos(ini, fim) {
            r.busca(p0).map_err(|e| format!("{caminho}: {e}"))?;
            while r.pos_virtual() < p1 && bam::le_registro(&mut r, &mut rec).map_err(|e| format!("{caminho}: {e}"))? {
                let i = bam::info(&rec).map_err(|e| format!("{caminho}: {e}"))?;
                if i.tid != tid {
                    continue;
                }
                if i.pos >= fim {
                    break 'regiao;
                }
                if i.fim > ini {
                    texto.clear();
                    bam::decodifica(&rec, &cab.refs, &mut texto).map_err(|e| format!("{caminho}: {e}"))?;
                    recebe(false, &String::from_utf8_lossy(&texto))?;
                    total += 1;
                }
            }
        }
    }
    Ok(total)
}
