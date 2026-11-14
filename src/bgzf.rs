// =====================================================================
//  BGZF -- gzip em blocos, com acesso aleatorio
//
//  Um arquivo BGZF e uma sequencia de membros gzip pequenos (no maximo
//  64 KiB cada), cada um com o proprio tamanho anotado num campo extra do
//  cabecalho. Qualquer leitor de gzip o le como um `.gz` comum; quem conhece
//  o formato pula direto para um bloco sem descomprimir os anteriores. E o
//  que o `bgzip` grava e o que o `tabix`, o `samtools` e o `bcftools` exigem
//  para indexar.
//
//  A posicao dentro do arquivo e um OFFSET VIRTUAL de 64 bits:
//
//      (offset do bloco no arquivo comprimido) << 16 | (offset dentro do bloco)
//
//  Como os blocos sao independentes, tambem da para comprimir e descomprimir
//  varios ao mesmo tempo -- e o que o escritor e o leitor sequencial fazem.
// =====================================================================

use flate2::{Compress, Compression, Crc, Decompress, FlushCompress, FlushDecompress, Status};
use rayon::prelude::*;
use std::io::{self, BufRead, Read, Seek, SeekFrom, Write};

/// Dado por bloco. O formato admite 64 KiB; o htslib usa 0xff00 para que o
/// bloco comprimido caiba em 64 KiB mesmo quando o dado nao comprime.
pub const DADO_POR_BLOCO: usize = 0xff00;
/// Tamanho maximo de um bloco inteiro: cabecalho, dado comprimido e rodape.
const BLOCO_MAX: usize = 0x10000;
const CABECALHO: usize = 18;
const RODAPE: usize = 8;
/// Blocos comprimidos (ou descomprimidos) de uma vez, em paralelo: ~4 MB.
const BLOCOS_POR_LOTE: usize = 64;
/// Nivel padrao de compressao: o mesmo do `bgzip`.
pub const NIVEL_PADRAO: u32 = 6;

/// O bloco vazio que marca o fim de um arquivo BGZF. Sem ele o `samtools`
/// avisa que o arquivo pode estar truncado.
pub const MARCA_DE_FIM: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00,
    0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

fn erro(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

// ---------------------------------------------------------------------
//  ESCRITA
// ---------------------------------------------------------------------

/// Comprime `dado` (no maximo `DADO_POR_BLOCO` bytes) num bloco BGZF completo,
/// acrescentado a `saida`.
fn comprime_bloco(c: &mut Compress, dado: &[u8], saida: &mut Vec<u8>) {
    let inicio = saida.len();
    saida.extend_from_slice(&[0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0, 0, 0]);
    let limite = BLOCO_MAX - CABECALHO - RODAPE;
    let p = inicio + CABECALHO;
    saida.resize(p + limite, 0);
    c.reset();
    let n = match c.compress(dado, &mut saida[p..], FlushCompress::Finish) {
        Ok(Status::StreamEnd) => c.total_out() as usize,
        // Nao coube: o dado nao comprime. Vai um bloco deflate "armazenado",
        // que custa 5 bytes alem do dado e sempre cabe.
        _ => {
            let len = dado.len() as u16;
            saida[p] = 1;
            saida[p + 1..p + 3].copy_from_slice(&len.to_le_bytes());
            saida[p + 3..p + 5].copy_from_slice(&(!len).to_le_bytes());
            saida[p + 5..p + 5 + dado.len()].copy_from_slice(dado);
            5 + dado.len()
        }
    };
    saida.truncate(p + n);
    let mut crc = Crc::new();
    crc.update(dado);
    saida.extend_from_slice(&crc.sum().to_le_bytes());
    saida.extend_from_slice(&(dado.len() as u32).to_le_bytes());
    let total = saida.len() - inicio;
    saida[inicio + 16..inicio + 18].copy_from_slice(&((total - 1) as u16).to_le_bytes());
}

/// Escritor BGZF. Junta o dado em blocos e comprime varios blocos ao mesmo
/// tempo. `termina` grava o resto e a marca de fim.
pub struct BgzfWriter<W: Write> {
    w: Option<W>,
    nivel: Compression,
    pendente: Vec<u8>,
}

impl<W: Write> BgzfWriter<W> {
    pub fn novo(w: W, nivel: u32) -> BgzfWriter<W> {
        BgzfWriter {
            w: Some(w),
            nivel: Compression::new(nivel.min(9)),
            pendente: Vec::with_capacity(BLOCOS_POR_LOTE * DADO_POR_BLOCO + (64 << 10)),
        }
    }

    /// Comprime e grava os blocos cheios; com `tudo`, tambem o bloco parcial.
    fn despeja(&mut self, tudo: bool) -> io::Result<()> {
        let cheios = self.pendente.len() / DADO_POR_BLOCO * DADO_POR_BLOCO;
        let ate = if tudo { self.pendente.len() } else { cheios };
        if ate == 0 {
            return Ok(());
        }
        let nivel = self.nivel;
        let blocos: Vec<Vec<u8>> = self.pendente[..ate]
            .par_chunks(DADO_POR_BLOCO)
            .map_init(
                || Compress::new(nivel, false),
                |c, dado| {
                    let mut s = Vec::with_capacity(dado.len() / 3 + 64);
                    comprime_bloco(c, dado, &mut s);
                    s
                },
            )
            .collect();
        let w = self.w.as_mut().ok_or_else(|| erro("BGZF output already closed".into()))?;
        for b in &blocos {
            w.write_all(b)?;
        }
        self.pendente.drain(..ate);
        Ok(())
    }

    /// Grava o que falta e a marca de fim, e devolve o escritor de baixo.
    pub fn termina(mut self) -> io::Result<W> {
        self.fecha()?;
        self.w.take().ok_or_else(|| erro("BGZF output already closed".into()))
    }

    fn fecha(&mut self) -> io::Result<()> {
        self.despeja(true)?;
        if let Some(w) = self.w.as_mut() {
            w.write_all(&MARCA_DE_FIM)?;
            w.flush()?;
        }
        Ok(())
    }
}

impl<W: Write> Write for BgzfWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.pendente.extend_from_slice(buf);
        if self.pendente.len() >= BLOCOS_POR_LOTE * DADO_POR_BLOCO {
            self.despeja(false)?;
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.despeja(true)?;
        match self.w.as_mut() {
            Some(w) => w.flush(),
            None => Ok(()),
        }
    }
}

impl<W: Write> Drop for BgzfWriter<W> {
    /// Quem esquece o `termina` ainda fica com um arquivo inteiro.
    fn drop(&mut self) {
        if self.w.is_some() {
            let _ = self.fecha();
        }
    }
}

// ---------------------------------------------------------------------
//  LEITURA
// ---------------------------------------------------------------------

/// Le um bloco inteiro, ainda comprimido, para `bloco`. Devolve `false` no fim
/// limpo do arquivo. `onde` e o offset do bloco, so para as mensagens.
fn le_bloco_cru<R: Read>(r: &mut R, bloco: &mut Vec<u8>, onde: u64) -> io::Result<bool> {
    bloco.clear();
    bloco.resize(12, 0);
    let mut lidos = 0;
    while lidos < 12 {
        match r.read(&mut bloco[lidos..12]) {
            Ok(0) => break,
            Ok(n) => lidos += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        }
    }
    if lidos == 0 {
        bloco.clear();
        return Ok(false);
    }
    if lidos < 12 {
        return Err(erro(format!("truncated BGZF block at byte {onde}: the file ends in the middle of a block header")));
    }
    if bloco[0] != 0x1f || bloco[1] != 0x8b {
        return Err(erro(format!("not a gzip block at byte {onde}: the file is damaged or is not BGZF")));
    }
    if bloco[2] != 8 || bloco[3] & 4 == 0 {
        return Err(erro(format!(
            "the gzip member at byte {onde} is not a BGZF block: the file was compressed with plain gzip, which cannot be indexed or read in blocks; recompress it as BGZF (write it again with a .gz output, or use bgzip)"
        )));
    }
    let xlen = u16::from_le_bytes([bloco[10], bloco[11]]) as usize;
    bloco.resize(12 + xlen, 0);
    r.read_exact(&mut bloco[12..]).map_err(|_| erro(format!("truncated BGZF block at byte {onde}")))?;
    // Procura o subcampo "BC", que traz o tamanho do bloco menos 1.
    let mut total = None;
    let mut i = 12;
    while i + 4 <= 12 + xlen {
        let slen = u16::from_le_bytes([bloco[i + 2], bloco[i + 3]]) as usize;
        if bloco[i] == b'B' && bloco[i + 1] == b'C' && slen == 2 && i + 6 <= 12 + xlen {
            total = Some(u16::from_le_bytes([bloco[i + 4], bloco[i + 5]]) as usize + 1);
            break;
        }
        i += 4 + slen;
    }
    let Some(total) = total else {
        return Err(erro(format!(
            "the gzip member at byte {onde} is not a BGZF block: the file was compressed with plain gzip, which cannot be indexed or read in blocks; recompress it as BGZF (write it again with a .gz output, or use bgzip)"
        )));
    };
    if total < 12 + xlen + RODAPE {
        return Err(erro(format!("damaged BGZF block at byte {onde}: impossible block size {total}")));
    }
    let ja = bloco.len();
    bloco.resize(total, 0);
    r.read_exact(&mut bloco[ja..])
        .map_err(|_| erro(format!("truncated BGZF block at byte {onde}: the file ends in the middle of a block")))?;
    Ok(true)
}

/// Descomprime um bloco lido por `le_bloco_cru` para `dado`, conferindo o
/// tamanho e o CRC.
fn descomprime_bloco(inflador: &mut Decompress, bloco: &[u8], dado: &mut Vec<u8>, onde: u64) -> io::Result<()> {
    let xlen = u16::from_le_bytes([bloco[10], bloco[11]]) as usize;
    let fim = bloco.len() - RODAPE;
    let crc_esperado = u32::from_le_bytes([bloco[fim], bloco[fim + 1], bloco[fim + 2], bloco[fim + 3]]);
    let tamanho = u32::from_le_bytes([bloco[fim + 4], bloco[fim + 5], bloco[fim + 6], bloco[fim + 7]]) as usize;
    if tamanho > BLOCO_MAX {
        return Err(erro(format!("damaged BGZF block at byte {onde}: declared size {tamanho} is over 64 KiB")));
    }
    dado.clear();
    dado.resize(tamanho, 0);
    inflador.reset(false);
    let status = inflador
        .decompress(&bloco[12 + xlen..fim], dado, FlushDecompress::Finish)
        .map_err(|e| erro(format!("damaged BGZF block at byte {onde}: {e}")))?;
    if inflador.total_out() as usize != tamanho || (tamanho > 0 && status != Status::StreamEnd) {
        return Err(erro(format!("damaged BGZF block at byte {onde}: the data does not have the declared size")));
    }
    let mut crc = Crc::new();
    crc.update(dado);
    if crc.sum() != crc_esperado {
        return Err(erro(format!("damaged BGZF block at byte {onde}: checksum mismatch")));
    }
    Ok(())
}

/// Leitor BGZF com acesso aleatorio: `busca` salta para um offset virtual e
/// `pos_virtual` diz onde se esta. Um bloco por vez, numa thread so -- e o
/// leitor das consultas por regiao e da construcao dos indices.
pub struct BgzfReader<R: Read> {
    r: R,
    /// Bloco corrente, descomprimido.
    dado: Vec<u8>,
    pos: usize,
    /// Offset, no arquivo comprimido, do bloco corrente e do seguinte.
    inicio: u64,
    proximo: u64,
    cru: Vec<u8>,
    inflador: Decompress,
}

impl<R: Read> BgzfReader<R> {
    pub fn novo(r: R) -> BgzfReader<R> {
        BgzfReader {
            r,
            dado: Vec::new(),
            pos: 0,
            inicio: 0,
            proximo: 0,
            cru: Vec::new(),
            inflador: Decompress::new(false),
        }
    }

    /// Carrega o proximo bloco que tenha dado. `false` no fim do arquivo.
    fn carrega(&mut self) -> io::Result<bool> {
        loop {
            if !le_bloco_cru(&mut self.r, &mut self.cru, self.proximo)? {
                self.dado.clear();
                self.pos = 0;
                return Ok(false);
            }
            descomprime_bloco(&mut self.inflador, &self.cru, &mut self.dado, self.proximo)?;
            self.inicio = self.proximo;
            self.proximo += self.cru.len() as u64;
            self.pos = 0;
            if !self.dado.is_empty() {
                return Ok(true);
            }
        }
    }

    /// Offset virtual do proximo byte a ler. No fim de um bloco aponta para o
    /// comeco do seguinte, como o htslib.
    pub fn pos_virtual(&self) -> u64 {
        if self.pos < self.dado.len() {
            (self.inicio << 16) | self.pos as u64
        } else {
            self.proximo << 16
        }
    }

}

impl<R: Read + Seek> BgzfReader<R> {
    /// Salta para um offset virtual.
    pub fn busca(&mut self, voff: u64) -> io::Result<()> {
        let (bloco, dentro) = (voff >> 16, (voff & 0xffff) as usize);
        self.r.seek(SeekFrom::Start(bloco))?;
        self.proximo = bloco;
        self.inicio = bloco;
        self.dado.clear();
        self.pos = 0;
        if dentro > 0 {
            if !self.carrega()? || self.inicio != bloco || dentro > self.dado.len() {
                return Err(erro(format!("invalid BGZF offset {voff}: the index does not match this file; build the index again")));
            }
            self.pos = dentro;
        }
        Ok(())
    }
}

impl<R: Read> Read for BgzfReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.pos >= self.dado.len() && !self.carrega()? {
            return Ok(0);
        }
        let n = buf.len().min(self.dado.len() - self.pos);
        buf[..n].copy_from_slice(&self.dado[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

impl<R: Read> BufRead for BgzfReader<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.pos >= self.dado.len() && !self.carrega()? {
            return Ok(&[]);
        }
        Ok(&self.dado[self.pos..])
    }

    fn consume(&mut self, n: usize) {
        self.pos = (self.pos + n).min(self.dado.len());
    }
}

/// Leitor BGZF so para leitura em sequencia: le um lote de blocos e os
/// descomprime em paralelo. Sem acesso aleatorio.
pub struct BgzfParLeitor<R: Read> {
    r: R,
    /// Blocos do lote corrente, ja descomprimidos, na ordem do arquivo.
    lote: Vec<Vec<u8>>,
    qual: usize,
    pos: usize,
    lido: u64,
    acabou: bool,
}

impl<R: Read> BgzfParLeitor<R> {
    pub fn novo(r: R) -> BgzfParLeitor<R> {
        BgzfParLeitor { r, lote: Vec::new(), qual: 0, pos: 0, lido: 0, acabou: false }
    }

    /// Garante que haja dado no bloco corrente. `false` no fim do arquivo.
    fn enche(&mut self) -> io::Result<bool> {
        loop {
            while self.qual < self.lote.len() {
                if self.pos < self.lote[self.qual].len() {
                    return Ok(true);
                }
                self.qual += 1;
                self.pos = 0;
            }
            if self.acabou {
                return Ok(false);
            }
            // Le os blocos comprimidos em sequencia (barato) e descomprime
            // todos de uma vez (caro, e o que se divide entre os nucleos).
            let mut crus: Vec<(u64, Vec<u8>)> = Vec::with_capacity(4 * BLOCOS_POR_LOTE);
            while crus.len() < 4 * BLOCOS_POR_LOTE {
                let mut b = Vec::new();
                if !le_bloco_cru(&mut self.r, &mut b, self.lido)? {
                    self.acabou = true;
                    break;
                }
                let onde = self.lido;
                self.lido += b.len() as u64;
                crus.push((onde, b));
            }
            self.lote = crus
                .par_iter()
                .map_init(
                    || Decompress::new(false),
                    |inflador, (onde, cru)| {
                        let mut dado = Vec::new();
                        descomprime_bloco(inflador, cru, &mut dado, *onde).map(|_| dado)
                    },
                )
                .collect::<io::Result<Vec<Vec<u8>>>>()?;
            self.qual = 0;
            self.pos = 0;
        }
    }
}

impl<R: Read> Read for BgzfParLeitor<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if !self.enche()? {
            return Ok(0);
        }
        let bloco = &self.lote[self.qual];
        let n = buf.len().min(bloco.len() - self.pos);
        buf[..n].copy_from_slice(&bloco[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

impl<R: Read> BufRead for BgzfParLeitor<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if !self.enche()? {
            return Ok(&[]);
        }
        Ok(&self.lote[self.qual][self.pos..])
    }

    fn consume(&mut self, n: usize) {
        self.pos += n;
    }
}

/// O arquivo comeca com um bloco BGZF? (Um `.gz` comum responde `false`.)
pub fn e_bgzf(caminho: &str) -> bool {
    let Ok(mut f) = std::fs::File::open(caminho) else {
        return false;
    };
    let mut c = [0u8; 18];
    if f.read_exact(&mut c).is_err() {
        return false;
    }
    c[0] == 0x1f && c[1] == 0x8b && c[2] == 8 && c[3] & 4 != 0 && c[12] == b'B' && c[13] == b'C'
}
