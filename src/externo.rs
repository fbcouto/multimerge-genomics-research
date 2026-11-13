// =====================================================================
//  ORDENACAO DE ARQUIVOS GENOMICOS DE QUALQUER TAMANHO
//
//  VCF, BED, GFF/GTF e SAM. Um arquivo que cabe na RAM e ordenado em memoria;
//  um que nao cabe e cortado em LOTES: cada lote e ordenado em memoria e
//  gravado em disco, comprimido; no fim, um merge k-vias junta os lotes
//  direto no destino.
//
//      entrada(s) --> [lote 0] [lote 1] ... [lote k]  --> merge k-vias --> destino
//                      ordenados pelo Multimerge_Stable
//
//  A memoria usada e limitada pelo orcamento, nao pelo tamanho da entrada.
//
//  ESTAVEL DE PONTA A PONTA. Dois registros com a mesma chave saem na ordem em
//  que entraram -- a ordem dos arquivos, e dentro de cada arquivo a ordem das
//  linhas. Isso depende de duas coisas:
//
//   1. cada lote e ordenado pelo motor estavel;
//   2. no merge, o empate entre lotes e decidido pelo lote MAIS ANTIGO.
//
//  A versao de maio de 2026 deste motor fazia (1) e nao (2): o heap comparava
//  so a chave, e linhas iguais vindas de lotes diferentes saiam em ordem
//  arbitraria. Num VCF isso e comum -- varias variantes na mesma posicao,
//  milhares com a mesma QUAL.
//
//  O destino e uma interface, nao um arquivo: a fusao de intervalos BED
//  (`intervalos.rs`) consome a saida ordenada em fluxo, sem gravar o arquivo
//  ordenado intermediario.
// =====================================================================

use crate::multimerge_core::multi_merge_sort_stable_by;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use lz4_flex::frame::{FrameDecoder, FrameEncoder};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

// ---------------------------------------------------------------------
//  O QUE ORDENAR E COMO
// ---------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Formato {
    Vcf,
    Bed,
    /// GFF3 e GTF: as colunas de sequencia, inicio e fim sao as mesmas.
    Gff,
    Sam,
}

impl Formato {
    pub fn de_texto(s: &str) -> Option<Formato> {
        match s.to_ascii_lowercase().as_str() {
            "vcf" => Some(Formato::Vcf),
            "bed" => Some(Formato::Bed),
            "gff" | "gff3" | "gtf" => Some(Formato::Gff),
            "sam" => Some(Formato::Sam),
            _ => None,
        }
    }

    /// Deduz o formato pela extensao, ignorando `.gz` e `.lz4`.
    pub fn pelo_nome(caminho: &str) -> Option<Formato> {
        let c = caminho.to_ascii_lowercase();
        let c = c.strip_suffix(".gz").or_else(|| c.strip_suffix(".lz4")).unwrap_or(&c);
        Formato::de_texto(c.rsplit('.').next()?)
    }

    fn e_cabecalho(self, linha: &str) -> bool {
        match self {
            Formato::Vcf | Formato::Gff => linha.starts_with('#'),
            Formato::Bed => {
                linha.starts_with('#') || linha.starts_with("track") || linha.starts_with("browser")
            }
            Formato::Sam => linha.starts_with('@'),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Criterio {
    /// Sequencia (cromossomo) e depois coordenada, crescente. VCF: CHROM, POS.
    /// BED: sequencia, inicio, fim. GFF: sequencia, inicio e, no empate, o
    /// mais longo primeiro (pai antes do filho). SAM: RNAME, POS, nao mapeados
    /// no fim.
    Posicao,
    /// So VCF: QUAL, da maior para a menor. QUAL ausente (".") vale zero.
    Qualidade,
    /// So SAM: QNAME, em ordem natural (`read2` antes de `read10`).
    Nome,
}

impl Criterio {
    pub fn de_texto(s: &str) -> Option<Criterio> {
        match s.to_ascii_lowercase().as_str() {
            "position" | "pos" | "coordinate" => Some(Criterio::Posicao),
            "quality" | "qual" => Some(Criterio::Qualidade),
            "name" | "qname" => Some(Criterio::Nome),
            _ => None,
        }
    }
}

/// Em que ordem os cromossomos saem.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum OrdemCromossomos {
    /// Byte a byte: chr1, chr10, chr11, chr2. E a do `bedtools sort` e a do
    /// `sort -k1,1`.
    Lexica,
    /// Numeros como numeros: chr1, chr2, chr10, chrM, chrX, chrY.
    Natural,
    /// A ordem em que as sequencias aparecem no cabecalho do primeiro arquivo
    /// (`##contig` no VCF, `@SQ` no SAM). E a do `bcftools sort` e a do
    /// `samtools sort`. Sequencia fora do cabecalho vai para o fim, em ordem
    /// natural.
    Cabecalho,
    /// Uma lista explicita -- tipicamente a primeira coluna de um `.fai` ou de
    /// um arquivo de tamanhos de cromossomo. O resto vai para o fim.
    Lista(Vec<String>),
}

impl OrdemCromossomos {
    /// "lexica", "natural", "cabecalho", ou o caminho de um arquivo cuja
    /// primeira coluna da a ordem (`.fai`, `chrom.sizes`).
    pub fn de_texto(s: &str) -> Result<OrdemCromossomos, String> {
        match s.to_ascii_lowercase().as_str() {
            "lexical" | "lexicographic" => Ok(OrdemCromossomos::Lexica),
            "natural" => Ok(OrdemCromossomos::Natural),
            "header" => Ok(OrdemCromossomos::Cabecalho),
            _ => {
                let txt = std::fs::read_to_string(s).map_err(|e| {
                    format!("Unknown chromosome order (use \"lexical\", \"natural\", \"header\" or the path of a .fai file), and '{s}' is not a readable file: {e}")
                })?;
                Ok(OrdemCromossomos::Lista(
                    txt.lines()
                        .filter(|l| !l.is_empty() && !l.starts_with('#'))
                        .filter_map(|l| l.split(['\t', ' ']).next())
                        .map(str::to_string)
                        .collect(),
                ))
            }
        }
    }
}

pub struct Opcoes {
    /// Orcamento de memoria para os registros, em bytes.
    pub memoria_bytes: usize,
    /// Onde gravar os lotes temporarios. `None` usa o diretorio temporario do
    /// sistema. Precisa de espaco para a entrada comprimida.
    pub dir_temp: Option<PathBuf>,
    pub ordem: OrdemCromossomos,
}

/// O que a ordenacao fez.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct Resumo {
    /// Registros entregues ao destino.
    pub registros: usize,
    /// Lotes gravados em disco. Zero quando tudo coube na memoria.
    pub lotes: usize,
    /// Linhas que nao eram cabecalho nem registro valido (colunas de menos,
    /// coordenada nao numerica) e ficaram de fora.
    pub ignoradas: usize,
}

/// Quem recebe a saida ordenada.
pub trait Destino {
    fn cabecalho(&mut self, linha: &str) -> Result<(), String>;
    fn registro(&mut self, linha: &str) -> Result<(), String>;
    fn fim(&mut self) -> Result<(), String>;
}

// ---------------------------------------------------------------------
//  ORDEM DOS CROMOSSOMOS
// ---------------------------------------------------------------------

/// Compara tratando sequencias de digitos como numeros.
pub fn natural_cmp(a: &str, b: &str) -> Ordering {
    let (x, y) = (a.as_bytes(), b.as_bytes());
    let (mut i, mut j) = (0usize, 0usize);
    while i < x.len() && j < y.len() {
        if x[i].is_ascii_digit() && y[j].is_ascii_digit() {
            let (si, sj) = (i, j);
            while i < x.len() && x[i].is_ascii_digit() {
                i += 1;
            }
            while j < y.len() && y[j].is_ascii_digit() {
                j += 1;
            }
            let (mut da, mut db) = (&x[si..i], &y[sj..j]);
            while da.len() > 1 && da[0] == b'0' {
                da = &da[1..];
            }
            while db.len() > 1 && db[0] == b'0' {
                db = &db[1..];
            }
            let o = da.len().cmp(&db.len()).then_with(|| da.cmp(db));
            if o != Ordering::Equal {
                return o;
            }
        } else {
            if x[i] != y[j] {
                return x[i].cmp(&y[j]);
            }
            i += 1;
            j += 1;
        }
    }
    // Um e prefixo do outro, ou so diferem em zeros a esquerda: a comparacao
    // byte a byte fecha a ordem total.
    (x.len() - i).cmp(&(y.len() - j)).then_with(|| a.cmp(b))
}

/// A regra de ordem dos nomes de sequencia, sem estado de execucao: serve a
/// ordenacao e a qualquer varredura que precise casar dois arquivos ordenados.
#[derive(Clone, Debug)]
pub struct OrdemNomes {
    natural: bool,
    /// Posicao fixa de cada nome (cabecalho ou lista). Vazio nas ordens lexica
    /// e natural.
    fixos: HashMap<String, u32>,
}

impl OrdemNomes {
    pub fn nova(ordem: &OrdemCromossomos) -> OrdemNomes {
        let mut o = OrdemNomes { natural: !matches!(ordem, OrdemCromossomos::Lexica), fixos: HashMap::new() };
        if let OrdemCromossomos::Lista(nomes) = ordem {
            for n in nomes {
                o.fixa(n);
            }
        }
        o
    }

    fn fixa(&mut self, nome: &str) {
        if !self.fixos.contains_key(nome) {
            let r = self.fixos.len() as u32;
            self.fixos.insert(nome.to_string(), r);
        }
    }

    fn posicao(&self, nome: &str) -> u32 {
        self.fixos.get(nome).copied().unwrap_or(u32::MAX)
    }

    pub fn compara(&self, a: &str, b: &str) -> Ordering {
        if a == b {
            return Ordering::Equal;
        }
        // "*" e o RNAME de um alinhamento nao mapeado: sempre por ultimo.
        let (ea, eb) = (a == "*", b == "*");
        if ea != eb {
            return ea.cmp(&eb);
        }
        self.posicao(a).cmp(&self.posicao(b)).then_with(|| {
            if self.natural {
                natural_cmp(a, b)
            } else {
                a.cmp(b)
            }
        })
    }
}

/// Nomes de sequencia vistos na execucao, cada um com um numero. A chave de um
/// registro guarda o numero; a comparacao so olha os nomes quando os numeros
/// diferem, o que e raro dentro de um cromossomo.
struct Cromossomos {
    ordem: OrdemNomes,
    ids: HashMap<String, u32>,
    nomes: Vec<String>,
}

impl Cromossomos {
    fn id(&mut self, nome: &str) -> u32 {
        if let Some(&i) = self.ids.get(nome) {
            return i;
        }
        let i = self.nomes.len() as u32;
        self.nomes.push(nome.to_string());
        self.ids.insert(nome.to_string(), i);
        i
    }

    #[inline]
    fn compara(&self, a: u32, b: u32) -> Ordering {
        if a == b {
            return Ordering::Equal;
        }
        self.ordem.compara(&self.nomes[a as usize], &self.nomes[b as usize])
    }
}

// ---------------------------------------------------------------------
//  CHAVES
// ---------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
struct Chave {
    crom: u32,
    a: u64,
    b: u64,
    /// Trecho da linha comparado como texto (criterio `Nome`).
    t0: u32,
    t1: u32,
}

/// Leva um f32 a um inteiro que preserva a ordem de `f32::total_cmp`.
#[inline]
fn ordem_total(q: f32) -> u32 {
    let b = q.to_bits();
    if b >> 31 == 1 {
        !b
    } else {
        b | 0x8000_0000
    }
}

struct Plano {
    formato: Formato,
    criterio: Criterio,
    crom: Cromossomos,
}

impl Plano {
    /// Extrai a chave. `None`: a linha nao e um registro valido.
    fn chave(&mut self, linha: &str) -> Option<Chave> {
        let mut cols = linha.split('\t');
        let (crom, a, b, t0, t1): (&str, u64, u64, u32, u32) = match (self.formato, self.criterio) {
            (Formato::Vcf, Criterio::Posicao) => {
                let crom = cols.next()?;
                let pos = cols.next()?.parse().ok()?;
                // Uma variante tem pelo menos CHROM POS ID REF ALT QUAL.
                cols.nth(3)?;
                (crom, pos, 0, 0, 0)
            }
            (Formato::Vcf, Criterio::Qualidade) => {
                cols.next()?;
                cols.next()?;
                let qual: f32 = cols.nth(3)?.parse().unwrap_or(0.0);
                // Decrescente: a maior qualidade vira a menor chave.
                ("", u64::from(!ordem_total(qual)), 0, 0, 0)
            }
            (Formato::Bed, Criterio::Posicao) => {
                let crom = cols.next()?;
                (crom, cols.next()?.parse().ok()?, cols.next()?.parse().ok()?, 0, 0)
            }
            (Formato::Gff, Criterio::Posicao) => {
                let crom = cols.next()?;
                let inicio: u64 = cols.nth(2)?.parse().ok()?;
                let fim: u64 = cols.next()?.parse().ok()?;
                // PAI ANTES DO FILHO. Um gene e seu primeiro exon comecam na
                // mesma base; quem contem o outro e o mais longo. Com o fim em
                // ordem DECRESCENTE no desempate, gene sai antes de mRNA antes
                // de exon -- a convencao que os leitores de GFF esperam. Gene e
                // mRNA de mesmas coordenadas ficam na ordem do arquivo, porque
                // a ordenacao e estavel.
                (crom, inicio, u64::MAX - fim, 0, 0)
            }
            (Formato::Sam, Criterio::Posicao) => {
                cols.next()?;
                cols.next()?;
                let crom = cols.next()?;
                (crom, cols.next()?.parse().ok()?, 0, 0, 0)
            }
            (Formato::Sam, Criterio::Nome) => {
                let nome = cols.next()?;
                // Um alinhamento tem pelo menos QNAME FLAG RNAME POS.
                cols.nth(2)?;
                ("", 0, 0, 0, nome.len() as u32)
            }
            _ => return None,
        };
        Some(Chave { crom: self.crom.id(crom), a, b, t0, t1 })
    }

    #[inline]
    fn compara(&self, ka: &Chave, la: &str, kb: &Chave, lb: &str) -> Ordering {
        match self.criterio {
            Criterio::Posicao => self
                .crom
                .compara(ka.crom, kb.crom)
                .then_with(|| ka.a.cmp(&kb.a))
                .then_with(|| ka.b.cmp(&kb.b)),
            Criterio::Qualidade => ka.a.cmp(&kb.a),
            Criterio::Nome => {
                natural_cmp(&la[ka.t0 as usize..ka.t1 as usize], &lb[kb.t0 as usize..kb.t1 as usize])
            }
        }
    }

    fn valida(formato: Formato, criterio: Criterio) -> Result<(), String> {
        match (formato, criterio) {
            (_, Criterio::Posicao) | (Formato::Vcf, Criterio::Qualidade) | (Formato::Sam, Criterio::Nome) => Ok(()),
            (f, c) => {
                let chave = match c {
                    Criterio::Posicao => "position",
                    Criterio::Qualidade => "quality",
                    Criterio::Nome => "name",
                };
                let nome = match f {
                    Formato::Vcf => "VCF",
                    Formato::Bed => "BED",
                    Formato::Gff => "GFF",
                    Formato::Sam => "SAM",
                };
                Err(format!("Sort key \"{chave}\" is not available for {nome} files"))
            }
        }
    }

    /// Nome de sequencia declarado numa linha de cabecalho, se houver.
    fn sequencia_do_cabecalho(&self, linha: &str) -> Option<String> {
        match self.formato {
            Formato::Vcf => {
                let r = linha.strip_prefix("##contig=<")?;
                let r = &r[r.find("ID=")? + 3..];
                Some(r[..r.find([',', '>']).unwrap_or(r.len())].to_string())
            }
            Formato::Sam => {
                linha.strip_prefix("@SQ")?;
                linha.split('\t').find_map(|c| c.strip_prefix("SN:")).map(str::to_string)
            }
            _ => None,
        }
    }
}

// ---------------------------------------------------------------------
//  ENTRADA E SAIDA
// ---------------------------------------------------------------------

/// Abre para leitura: `.gz` (inclusive BGZF, que e gzip de varios membros),
/// `.lz4` ou texto puro.
pub fn abre_leitura(caminho: &str) -> Result<Box<dyn BufRead>, String> {
    let file = File::open(caminho).map_err(|e| format!("Error opening {caminho}: {e}"))?;
    let c = caminho.to_ascii_lowercase();
    Ok(if c.ends_with(".gz") {
        Box::new(BufReader::with_capacity(16 << 20, MultiGzDecoder::new(file)))
    } else if c.ends_with(".lz4") {
        Box::new(BufReader::with_capacity(4 << 20, FrameDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(16 << 20, file))
    })
}

/// Escritor que escolhe a compressao pela extensao: `.gz`, `.lz4` ou nenhuma.
pub enum Escritor {
    Puro(BufWriter<File>),
    Gz(GzEncoder<BufWriter<File>>),
    Lz4(FrameEncoder<BufWriter<File>>),
}

impl Escritor {
    pub fn cria(caminho: &str) -> Result<Escritor, String> {
        let file = File::create(caminho).map_err(|e| format!("Error creating {caminho}: {e}"))?;
        let w = BufWriter::with_capacity(16 << 20, file);
        let c = caminho.to_ascii_lowercase();
        Ok(if c.ends_with(".gz") {
            Escritor::Gz(GzEncoder::new(w, flate2::Compression::fast()))
        } else if c.ends_with(".lz4") {
            Escritor::Lz4(FrameEncoder::new(w))
        } else {
            Escritor::Puro(w)
        })
    }

    pub fn linha(&mut self, linha: &str) -> Result<(), String> {
        self.write_all(linha.as_bytes())
            .and_then(|_| self.write_all(b"\n"))
            .map_err(|e| e.to_string())
    }

    /// Fecha o fluxo de compressao e descarrega. Sem isto um `.gz` ou `.lz4`
    /// fica truncado.
    pub fn termina(self) -> Result<(), String> {
        match self {
            Escritor::Puro(mut w) => w.flush(),
            Escritor::Gz(e) => e.finish().and_then(|mut w| w.flush()),
            Escritor::Lz4(e) => e
                .finish()
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))
                .and_then(|mut w| w.flush()),
        }
        .map_err(|e| e.to_string())
    }
}

impl Write for Escritor {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Escritor::Puro(w) => w.write(buf),
            Escritor::Gz(w) => w.write(buf),
            Escritor::Lz4(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Escritor::Puro(w) => w.flush(),
            Escritor::Gz(w) => w.flush(),
            Escritor::Lz4(w) => w.flush(),
        }
    }
}

/// Destino que grava cabecalho e registros num arquivo.
pub struct DestinoArquivo {
    w: Option<Escritor>,
}

impl DestinoArquivo {
    pub fn cria(caminho: &str) -> Result<DestinoArquivo, String> {
        Ok(DestinoArquivo { w: Some(Escritor::cria(caminho)?) })
    }
}

impl Destino for DestinoArquivo {
    fn cabecalho(&mut self, linha: &str) -> Result<(), String> {
        self.registro(linha)
    }
    fn registro(&mut self, linha: &str) -> Result<(), String> {
        self.w.as_mut().ok_or("output already closed")?.linha(linha)
    }
    fn fim(&mut self) -> Result<(), String> {
        self.w.take().ok_or("output already closed")?.termina()
    }
}

/// Diretorio temporario proprio de uma operacao. Some junto com ela, inclusive
/// quando a operacao termina em erro -- e nunca apaga nada que nao criou.
pub struct DirTemp {
    caminho: PathBuf,
    arquivos: Vec<PathBuf>,
}

impl DirTemp {
    pub fn novo(base: Option<&Path>) -> Result<DirTemp, String> {
        let base = base.map(Path::to_path_buf).unwrap_or_else(std::env::temp_dir);
        let carimbo = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        let caminho = base.join(format!("mmgen_{}_{}", std::process::id(), carimbo));
        std::fs::create_dir_all(&caminho).map_err(|e| format!("Error creating {}: {e}", caminho.display()))?;
        Ok(DirTemp { caminho, arquivos: Vec::new() })
    }

    /// Reserva um nome de arquivo dentro do diretorio. O arquivo e apagado
    /// quando o diretorio sai de escopo, exista ele ou nao.
    pub fn arquivo(&mut self, nome: &str) -> String {
        let p = self.caminho.join(nome);
        self.arquivos.push(p.clone());
        p.to_string_lossy().to_string()
    }
}

impl Drop for DirTemp {
    fn drop(&mut self) {
        for a in &self.arquivos {
            let _ = std::fs::remove_file(a);
        }
        let _ = std::fs::remove_dir(&self.caminho);
    }
}

// ---------------------------------------------------------------------
//  LOTES
// ---------------------------------------------------------------------

struct Registro {
    chave: Chave,
    linha: String,
}

impl Registro {
    /// Memoria que este registro ocupa, aproximadamente: o texto, a
    /// estrutura, a referencia usada na ordenacao e a folga do alocador.
    fn bytes(&self) -> usize {
        self.linha.capacity() + std::mem::size_of::<Registro>() + 8 + 32
    }
}

/// Ordena o lote com o motor estavel. Ordena REFERENCIAS: mover 8 bytes por
/// registro em vez de um `String` e uma chave.
fn ordena_lote<'a>(lote: &'a [Registro], plano: &Plano) -> Vec<&'a Registro> {
    let mut ordem: Vec<&Registro> = lote.iter().collect();
    multi_merge_sort_stable_by(&mut ordem, |x, y| plano.compara(&x.chave, &x.linha, &y.chave, &y.linha));
    ordem
}

fn grava_lote(lote: &[Registro], plano: &Plano, temp: &mut DirTemp, lotes: &mut Vec<String>) -> Result<(), String> {
    let ordem = ordena_lote(lote, plano);
    // O nome entra no diretorio temporario antes da escrita: se ela falhar no
    // meio, o arquivo parcial e apagado junto com o resto.
    let caminho = temp.arquivo(&format!("lote_{}.lz4", lotes.len()));
    let mut w = Escritor::cria(&caminho)?;
    lotes.push(caminho);
    for r in ordem {
        w.linha(&r.linha)?;
    }
    w.termina()
}

/// Merge k-vias dos lotes. Um heap de indices de lote; a cabeca de cada lote
/// (chave e linha) mora fora do heap, em buffers reaproveitados a cada leitura.
struct Fusao<'p> {
    plano: &'p mut Plano,
    leitores: Vec<Box<dyn BufRead>>,
    linhas: Vec<String>,
    chaves: Vec<Chave>,
    heap: Vec<usize>,
}

impl Fusao<'_> {
    /// O lote `i` sai antes do lote `j`? No empate ganha o MAIS ANTIGO: e isso
    /// que mantem o merge estavel.
    #[inline]
    fn antes(&self, i: usize, j: usize) -> bool {
        match self.plano.compara(&self.chaves[i], &self.linhas[i], &self.chaves[j], &self.linhas[j]) {
            Ordering::Less => true,
            Ordering::Greater => false,
            Ordering::Equal => i < j,
        }
    }

    fn desce(&mut self, mut p: usize) {
        let n = self.heap.len();
        loop {
            let (e, d) = (2 * p + 1, 2 * p + 2);
            let mut m = p;
            if e < n && self.antes(self.heap[e], self.heap[m]) {
                m = e;
            }
            if d < n && self.antes(self.heap[d], self.heap[m]) {
                m = d;
            }
            if m == p {
                return;
            }
            self.heap.swap(p, m);
            p = m;
        }
    }

    /// Le o proximo registro do lote `i`. `false` no fim do lote.
    fn avanca(&mut self, i: usize) -> Result<bool, String> {
        let linha = &mut self.linhas[i];
        linha.clear();
        if self.leitores[i].read_line(linha).map_err(|e| e.to_string())? == 0 {
            return Ok(false);
        }
        if linha.ends_with('\n') {
            linha.pop();
        }
        self.chaves[i] = self.plano.chave(&self.linhas[i]).ok_or("corrupted temporary run")?;
        Ok(true)
    }
}

fn funde_lotes(lotes: &[String], plano: &mut Plano, destino: &mut dyn Destino) -> Result<usize, String> {
    let k = lotes.len();
    let mut f = Fusao {
        plano,
        leitores: Vec::with_capacity(k),
        linhas: vec![String::new(); k],
        chaves: vec![Chave { crom: 0, a: 0, b: 0, t0: 0, t1: 0 }; k],
        heap: Vec::with_capacity(k),
    };
    for (i, caminho) in lotes.iter().enumerate() {
        f.leitores.push(abre_leitura(caminho)?);
        if f.avanca(i)? {
            f.heap.push(i);
        }
    }
    for p in (0..f.heap.len() / 2).rev() {
        f.desce(p);
    }

    let mut total = 0usize;
    while let Some(&i) = f.heap.first() {
        destino.registro(&f.linhas[i])?;
        total += 1;
        if !f.avanca(i)? {
            let ultimo = f.heap.pop().unwrap();
            if f.heap.is_empty() {
                break;
            }
            f.heap[0] = ultimo;
        }
        f.desce(0);
    }
    Ok(total)
}

// ---------------------------------------------------------------------
//  CABECALHO
// ---------------------------------------------------------------------

const CABECALHO_VCF_MINIMO: [&str; 3] = [
    "##fileformat=VCFv4.2",
    "##source=MMGEN_Multimerge_Engine",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
];

/// Ajusta o cabecalho ao que a saida passou a ser. Um SAM ordenado tem que
/// dizer isso na linha `@HD`, ou as ferramentas que o leem depois ordenam de
/// novo -- ou, pior, confiam numa ordem que nao existe mais.
fn ajusta_cabecalho(mut linhas: Vec<String>, formato: Formato, criterio: Criterio) -> Vec<String> {
    match formato {
        Formato::Vcf if linhas.is_empty() => CABECALHO_VCF_MINIMO.iter().map(|s| s.to_string()).collect(),
        Formato::Sam => {
            let so = if criterio == Criterio::Nome { "queryname" } else { "coordinate" };
            if linhas.first().is_some_and(|l| l.starts_with("@HD")) {
                let mut campos: Vec<String> =
                    linhas[0].split('\t').filter(|c| !c.starts_with("SO:")).map(str::to_string).collect();
                campos.push(format!("SO:{so}"));
                linhas[0] = campos.join("\t");
            } else {
                linhas.insert(0, format!("@HD\tVN:1.6\tSO:{so}"));
            }
            linhas
        }
        _ => linhas,
    }
}

// ---------------------------------------------------------------------
//  ENTRADA PRINCIPAL
// ---------------------------------------------------------------------

/// Ordena um ou mais arquivos e entrega o resultado, em fluxo, a `destino`.
///
/// O cabecalho do PRIMEIRO arquivo e preservado; os dos demais sao
/// descartados. Um VCF sem cabecalho recebe um cabecalho minimo.
pub fn ordenar(
    caminhos: &[String],
    formato: Formato,
    criterio: Criterio,
    op: &Opcoes,
    destino: &mut dyn Destino,
) -> Result<Resumo, String> {
    Plano::valida(formato, criterio)?;
    let orcamento = op.memoria_bytes.max(1 << 20);
    let por_cabecalho = op.ordem == OrdemCromossomos::Cabecalho;
    let mut plano = Plano {
        formato,
        criterio,
        crom: Cromossomos { ordem: OrdemNomes::nova(&op.ordem), ids: HashMap::new(), nomes: Vec::new() },
    };

    let mut cabecalho: Vec<String> = Vec::new();
    let mut lote: Vec<Registro> = Vec::new();
    let mut usado = 0usize;
    let mut temp: Option<DirTemp> = None;
    let mut lotes: Vec<String> = Vec::new();
    let mut ignoradas = 0usize;
    let mut linha = String::new();

    for (n_arq, caminho) in caminhos.iter().enumerate() {
        let mut leitor = abre_leitura(caminho)?;
        loop {
            linha.clear();
            if leitor.read_line(&mut linha).map_err(|e| format!("Error reading {caminho}: {e}"))? == 0 {
                break;
            }
            let texto = linha.trim_end_matches(['\n', '\r']);
            if texto.is_empty() {
                continue;
            }
            if formato.e_cabecalho(texto) {
                // `###` no GFF3 diz "tudo acima ja esta resolvido". Depois de
                // reordenar as linhas ele nao significa mais nada, e um GFF do
                // Ensembl traz dezenas de milhares.
                if formato == Formato::Gff && texto == "###" {
                    continue;
                }
                if formato == Formato::Gff && texto.starts_with("##FASTA") {
                    return Err(format!(
                        "{caminho}: a GFF with an embedded ##FASTA section cannot be sorted line by line; split the sequence off first"
                    ));
                }
                // So o cabecalho do primeiro arquivo fixa a ordem das
                // sequencias: ele e lido antes de qualquer registro, entao
                // todos os lotes enxergam a mesma ordem.
                if n_arq == 0 {
                    if por_cabecalho {
                        if let Some(seq) = plano.sequencia_do_cabecalho(texto) {
                            plano.crom.ordem.fixa(&seq);
                        }
                    }
                    cabecalho.push(texto.to_string());
                }
                continue;
            }
            let Some(chave) = plano.chave(texto) else {
                ignoradas += 1;
                continue;
            };
            let r = Registro { chave, linha: texto.to_string() };
            usado += r.bytes();
            lote.push(r);

            if usado >= orcamento {
                if temp.is_none() {
                    temp = Some(DirTemp::novo(op.dir_temp.as_deref())?);
                }
                grava_lote(&lote, &plano, temp.as_mut().unwrap(), &mut lotes)?;
                lote.clear();
                usado = 0;
            }
        }
    }

    for l in ajusta_cabecalho(cabecalho, formato, criterio) {
        destino.cabecalho(&l)?;
    }

    let registros = match temp.as_mut() {
        // Tudo coube: ordena e entrega direto, sem tocar no disco temporario.
        None => {
            let ordem = ordena_lote(&lote, &plano);
            for r in &ordem {
                destino.registro(&r.linha)?;
            }
            ordem.len()
        }
        Some(temp) => {
            if !lote.is_empty() {
                grava_lote(&lote, &plano, temp, &mut lotes)?;
            }
            // Solta o ultimo lote antes do merge: a partir daqui a memoria e
            // so um buffer de leitura por lote.
            lote = Vec::new();
            funde_lotes(&lotes, &mut plano, destino)?
        }
    };
    drop(lote);
    destino.fim()?;
    Ok(Resumo { registros, lotes: lotes.len(), ignoradas })
}

/// Ordena um ou mais arquivos para `saida`. A extensao de `saida` escolhe a
/// compressao: `.gz`, `.lz4` ou nenhuma.
pub fn ordenar_arquivos(
    caminhos: &[String],
    saida: &str,
    formato: Formato,
    criterio: Criterio,
    op: &Opcoes,
) -> Result<Resumo, String> {
    let mut destino = DestinoArquivo::cria(saida)?;
    ordenar(caminhos, formato, criterio, op, &mut destino)
}
