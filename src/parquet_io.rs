// =====================================================================
//  PARQUET como envelope de VCF, BED, GFF e SAM
//
//  Um `.parquet` aqui e so mais uma forma de guardar as mesmas linhas: o
//  escritor recebe o texto que iria para um `.vcf` (ou `.bed`, `.gff`,
//  `.sam`) e grava colunas tipadas; o leitor devolve as linhas. Com isso
//  toda ferramenta do projeto le e grava Parquet sem saber que ele existe.
//
//  O que vai para o arquivo:
//
//   * uma coluna por campo do formato, com os nomes do proprio formato
//     (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT e uma coluna por
//     amostra; chrom, start, end...). Coordenadas sao inteiros de 64 bits;
//     QUAL do VCF e score do GFF sao numeros, com `.` virando nulo -- ou
//     texto, tal como estava escrito, se o usuario pedir
//     (`define_numeros_como_texto`): ai a volta ao texto e byte a byte;
//   * o cabecalho, inteiro, nos metadados do arquivo (`mmgen.header`), e o
//     formato em `mmgen.format`. E o que permite voltar ao texto sem perda.
//
//  A escrita e em grupos de linhas, e a leitura tambem: a memoria nao cresce
//  com o tamanho do arquivo. O Parquet escrito por outra ferramenta (pandas,
//  DuckDB, Spark) e lido pelos nomes das colunas.
// =====================================================================

use crate::externo::Formato;
use polars::prelude::{DataType, Series};
use polars_arrow::array::{Array, MutablePrimitiveArray, MutableUtf8Array, PrimitiveArray, Utf8Array};
use polars_arrow::chunk::Chunk;
use polars_arrow::datatypes::{ArrowDataType, ArrowSchema, Field, PhysicalType};
use polars_arrow::types::PrimitiveType;
use polars_parquet::read::{infer_schema, read_metadata, FileReader};
use polars_parquet::write::{
    transverse, CompressionOptions, Encoding, FileWriter, KeyValue, RowGroupIterator, Version, WriteOptions,
};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::sync::atomic::{AtomicBool, Ordering};

/// QUAL (VCF) e score (GFF) gravados como texto em vez de numero. Vale para
/// todo Parquet criado depois, no processo inteiro (a API Python expoe como
/// `parquet_options(numbers=...)`).
static NUMEROS_COMO_TEXTO: AtomicBool = AtomicBool::new(false);

pub fn define_numeros_como_texto(sim: bool) {
    NUMEROS_COMO_TEXTO.store(sim, Ordering::Relaxed);
}

pub fn numeros_como_texto() -> bool {
    NUMEROS_COMO_TEXTO.load(Ordering::Relaxed)
}

const CHAVE_FORMATO: &str = "mmgen.format";
const CHAVE_CABECALHO: &str = "mmgen.header";
/// Um grupo de linhas fecha com este numero de registros ou este tanto de
/// texto, o que vier primeiro.
const LINHAS_POR_GRUPO: usize = 1 << 17;
const BYTES_POR_GRUPO: usize = 64 << 20;

fn erro(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

fn nome_formato(f: Formato) -> &'static str {
    match f {
        Formato::Vcf => "vcf",
        Formato::Bed => "bed",
        Formato::Gff => "gff",
        Formato::Sam => "sam",
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Tipo {
    Texto,
    Inteiro,
    Real,
}

const VCF: [(&str, Tipo); 8] = [
    ("CHROM", Tipo::Texto),
    ("POS", Tipo::Inteiro),
    ("ID", Tipo::Texto),
    ("REF", Tipo::Texto),
    ("ALT", Tipo::Texto),
    ("QUAL", Tipo::Real),
    ("FILTER", Tipo::Texto),
    ("INFO", Tipo::Texto),
];
const BED: [(&str, Tipo); 12] = [
    ("chrom", Tipo::Texto),
    ("start", Tipo::Inteiro),
    ("end", Tipo::Inteiro),
    ("name", Tipo::Texto),
    ("score", Tipo::Texto),
    ("strand", Tipo::Texto),
    ("thickStart", Tipo::Texto),
    ("thickEnd", Tipo::Texto),
    ("itemRgb", Tipo::Texto),
    ("blockCount", Tipo::Texto),
    ("blockSizes", Tipo::Texto),
    ("blockStarts", Tipo::Texto),
];
const GFF: [(&str, Tipo); 9] = [
    ("seqid", Tipo::Texto),
    ("source", Tipo::Texto),
    ("type", Tipo::Texto),
    ("start", Tipo::Inteiro),
    ("end", Tipo::Inteiro),
    ("score", Tipo::Real),
    ("strand", Tipo::Texto),
    ("phase", Tipo::Texto),
    ("attributes", Tipo::Texto),
];
const SAM: [(&str, Tipo); 12] = [
    ("QNAME", Tipo::Texto),
    ("FLAG", Tipo::Inteiro),
    ("RNAME", Tipo::Texto),
    ("POS", Tipo::Inteiro),
    ("MAPQ", Tipo::Inteiro),
    ("CIGAR", Tipo::Texto),
    ("RNEXT", Tipo::Texto),
    ("PNEXT", Tipo::Inteiro),
    ("TLEN", Tipo::Inteiro),
    ("SEQ", Tipo::Texto),
    ("QUAL", Tipo::Texto),
    ("TAGS", Tipo::Texto),
];

fn colunas_fixas(f: Formato) -> &'static [(&'static str, Tipo)] {
    match f {
        Formato::Vcf => &VCF,
        Formato::Bed => &BED,
        Formato::Gff => &GFF,
        Formato::Sam => &SAM,
    }
}

/// Quantas colunas uma linha precisa ter para ser gravada.
fn obrigatorias(f: Formato) -> usize {
    match f {
        Formato::Vcf => 5,
        Formato::Bed => 3,
        Formato::Gff => 5,
        Formato::Sam => 11,
    }
}

/// Quantas colunas a linha tem ao voltar ao texto, no minimo: as que faltarem
/// saem como `.`.
fn minimo_na_leitura(f: Formato) -> usize {
    match f {
        Formato::Vcf => 8,
        Formato::Bed => 3,
        Formato::Gff => 9,
        Formato::Sam => 11,
    }
}

fn e_cabecalho(f: Option<Formato>, linha: &str) -> bool {
    match f {
        Some(Formato::Vcf | Formato::Gff) => linha.starts_with('#'),
        Some(Formato::Bed) => linha.starts_with('#') || linha.starts_with("track") || linha.starts_with("browser"),
        Some(Formato::Sam) => linha.starts_with('@'),
        None => linha.starts_with('#') || linha.starts_with('@') || linha.starts_with("track ") || linha.starts_with("browser "),
    }
}

fn e_numero(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

/// Descobre o formato pelo que chegou: o cabecalho, se houver, e a primeira
/// linha de dados.
fn fareja(cabecalho: &[String], primeira: &str) -> Option<Formato> {
    if cabecalho.iter().any(|l| l.starts_with("##fileformat=VCF") || l.starts_with("#CHROM\t")) {
        return Some(Formato::Vcf);
    }
    if cabecalho.iter().any(|l| l.starts_with("##gff-version") || l.starts_with("##gtf")) {
        return Some(Formato::Gff);
    }
    if cabecalho.iter().any(|l| ["@HD\t", "@SQ\t", "@RG\t", "@PG\t", "@CO\t"].iter().any(|p| l.starts_with(p))) {
        return Some(Formato::Sam);
    }
    let c: Vec<&str> = primeira.split('\t').collect();
    if c.len() >= 11 && e_numero(c[1]) && e_numero(c[3]) && e_numero(c[4]) && (c[5] == "*" || c[5].bytes().all(|b| b.is_ascii_digit() || b"MIDNSHP=X".contains(&b))) {
        return Some(Formato::Sam);
    }
    if c.len() == 9 && e_numero(c[3]) && e_numero(c[4]) && !e_numero(c[2]) && ["+", "-", ".", "?"].contains(&c[6]) {
        return Some(Formato::Gff);
    }
    if c.len() >= 3 && e_numero(c[1]) && e_numero(c[2]) {
        return Some(Formato::Bed);
    }
    if c.len() >= 5 && e_numero(c[1]) {
        return Some(Formato::Vcf);
    }
    None
}

/// O formato de um nome como `x.vcf.parquet`.
pub fn formato_pelo_nome(caminho: &str) -> Option<Formato> {
    let c = caminho.to_ascii_lowercase();
    let c = c.strip_suffix(".parquet").unwrap_or(&c);
    Formato::de_texto(c.rsplit('.').next()?)
}

// ---------------------------------------------------------------------
//  ESCRITA
// ---------------------------------------------------------------------

enum Coluna {
    Texto(MutableUtf8Array<i64>),
    Inteiro(MutablePrimitiveArray<i64>),
    Real(MutablePrimitiveArray<f64>),
}

impl Coluna {
    fn nova(t: Tipo) -> Coluna {
        match t {
            Tipo::Texto => Coluna::Texto(MutableUtf8Array::new()),
            Tipo::Inteiro => Coluna::Inteiro(MutablePrimitiveArray::new()),
            Tipo::Real => Coluna::Real(MutablePrimitiveArray::new()),
        }
    }

    fn nulo(&mut self) {
        match self {
            Coluna::Texto(a) => a.push::<&str>(None),
            Coluna::Inteiro(a) => a.push(None),
            Coluna::Real(a) => a.push(None),
        }
    }

    /// Guarda o valor de um campo. `Err`: o texto nao e do tipo da coluna.
    fn poe(&mut self, v: &str) -> Result<(), ()> {
        match self {
            Coluna::Texto(a) => a.push(Some(v)),
            Coluna::Inteiro(a) => a.push(Some(v.parse::<i64>().map_err(|_| ())?)),
            Coluna::Real(a) => {
                if v == "." || v.is_empty() {
                    a.push(None)
                } else {
                    a.push(Some(v.parse::<f64>().map_err(|_| ())?))
                }
            }
        }
        Ok(())
    }

    fn fecha(&mut self) -> Box<dyn Array> {
        match self {
            Coluna::Texto(a) => {
                let pronto: Utf8Array<i64> = std::mem::take(a).into();
                pronto.boxed()
            }
            Coluna::Inteiro(a) => {
                let pronto: PrimitiveArray<i64> = std::mem::take(a).into();
                pronto.boxed()
            }
            Coluna::Real(a) => {
                let pronto: PrimitiveArray<f64> = std::mem::take(a).into();
                pronto.boxed()
            }
        }
    }
}

/// A mesma escolha de codificacao do Polars: dicionario para texto e inteiros
/// (com queda para o simples quando ha valores distintos demais).
fn codificacao(tipo: &ArrowDataType) -> Encoding {
    match tipo.to_physical_type() {
        PhysicalType::Dictionary(_) | PhysicalType::LargeBinary | PhysicalType::LargeUtf8 => Encoding::RleDictionary,
        PhysicalType::Primitive(PrimitiveType::Float16 | PrimitiveType::Float32 | PrimitiveType::Float64) => Encoding::Plain,
        PhysicalType::Primitive(_) => Encoding::RleDictionary,
        _ => Encoding::Plain,
    }
}

/// O arquivo ja aberto, depois que o formato e as colunas ficaram conhecidos.
struct Aberto {
    formato: Formato,
    w: FileWriter<BufWriter<File>>,
    esquema: ArrowSchema,
    codificacoes: Vec<Vec<Encoding>>,
    opcoes: WriteOptions,
    nomes: Vec<String>,
    colunas: Vec<Coluna>,
    no_grupo: usize,
    bytes_no_grupo: usize,
}

/// Recebe as linhas de um VCF, BED, GFF ou SAM e grava Parquet.
pub struct ParquetEscritor {
    caminho: String,
    formato: Option<Formato>,
    /// Tabela derivada de um BED (a saida de `closest` ou `coverage`): so as
    /// tres primeiras colunas tem nome de BED; as demais sao `column_N`.
    tabela: bool,
    /// O valor de `numeros_como_texto()` quando o arquivo foi criado.
    reais_como_texto: bool,
    cabecalho: Vec<String>,
    resto: Vec<u8>,
    aberto: Option<Aberto>,
    linha_n: u64,
    registros: u64,
}

fn nomes_unicos(nomes: Vec<String>) -> Vec<String> {
    let mut vistos = std::collections::HashSet::new();
    nomes
        .into_iter()
        .map(|n| {
            let mut nome = n.clone();
            let mut k = 2;
            while !vistos.insert(nome.to_ascii_lowercase()) {
                nome = format!("{n}_{k}");
                k += 1;
            }
            nome
        })
        .collect()
}

impl ParquetEscritor {
    /// `formato`: quando quem chama sabe o que esta gravando. Sem ele, vale o
    /// nome (`x.vcf.parquet`) e, por fim, o proprio conteudo.
    pub fn cria(caminho: &str, formato: Option<Formato>, tabela: bool) -> Result<ParquetEscritor, String> {
        // Cria o arquivo ja agora: um caminho invalido falha antes do trabalho.
        File::create(caminho).map_err(|e| format!("Error creating {caminho}: {e}"))?;
        Ok(ParquetEscritor {
            caminho: caminho.to_string(),
            formato: formato.or_else(|| formato_pelo_nome(caminho)),
            tabela,
            reais_como_texto: numeros_como_texto(),
            cabecalho: Vec::new(),
            resto: Vec::new(),
            aberto: None,
            linha_n: 0,
            registros: 0,
        })
    }

    /// Fixa formato e colunas e abre o arquivo. `primeira`: a primeira linha
    /// de dados, ou `None` se nao ha nenhuma.
    fn abre(&mut self, primeira: Option<&str>) -> io::Result<()> {
        let formato = match self.formato.or_else(|| fareja(&self.cabecalho, primeira.unwrap_or(""))) {
            Some(f) => f,
            None => {
                return Err(erro(format!(
                    "Cannot tell which format is being written to {}: name the output like out.vcf.parquet, out.bed.parquet, out.gff.parquet or out.sam.parquet",
                    self.caminho
                )))
            }
        };
        self.formato = Some(formato);
        let fixas = colunas_fixas(formato);
        let n_primeira = primeira.map_or(0, |l| l.split('\t').count());
        let mut nomes: Vec<String> = Vec::new();
        let mut tipos: Vec<Tipo> = Vec::new();
        match formato {
            Formato::Vcf => {
                // As colunas alem das 8 fixas vem da linha #CHROM (FORMAT e os
                // nomes das amostras) ou, sem ela, do que a primeira linha tem.
                let titulo = self.cabecalho.iter().rev().find(|l| l.starts_with("#CHROM"));
                let extras: Vec<String> = match titulo {
                    Some(t) => t.split('\t').skip(8).map(str::to_string).collect(),
                    None => (8..n_primeira).map(|i| if i == 8 { "FORMAT".to_string() } else { format!("sample_{}", i - 8) }).collect(),
                };
                nomes.extend(fixas.iter().map(|c| c.0.to_string()));
                tipos.extend(fixas.iter().map(|c| c.1));
                tipos.extend(extras.iter().map(|_| Tipo::Texto));
                nomes.extend(extras);
            }
            Formato::Bed => {
                let n = n_primeira.max(3);
                for i in 0..n {
                    match fixas.get(i).filter(|_| i < 3 || !self.tabela) {
                        Some(c) => {
                            nomes.push(c.0.to_string());
                            tipos.push(c.1);
                        }
                        None => {
                            nomes.push(format!("column_{}", i + 1));
                            tipos.push(Tipo::Texto);
                        }
                    }
                }
            }
            Formato::Gff | Formato::Sam => {
                nomes.extend(fixas.iter().map(|c| c.0.to_string()));
                tipos.extend(fixas.iter().map(|c| c.1));
            }
        }
        if self.reais_como_texto {
            for t in tipos.iter_mut().filter(|t| **t == Tipo::Real) {
                *t = Tipo::Texto;
            }
        }
        let nomes = nomes_unicos(nomes);
        let campos: Vec<Field> = nomes
            .iter()
            .zip(&tipos)
            .map(|(n, t)| {
                let tipo = match t {
                    Tipo::Texto => ArrowDataType::LargeUtf8,
                    Tipo::Inteiro => ArrowDataType::Int64,
                    Tipo::Real => ArrowDataType::Float64,
                };
                Field::new(n.clone(), tipo, true)
            })
            .collect();
        let esquema = ArrowSchema::from(campos);
        let opcoes = WriteOptions {
            write_statistics: true,
            version: Version::V2,
            compression: CompressionOptions::Zstd(None),
            data_pagesize_limit: None,
        };
        let codificacoes = esquema.fields.iter().map(|f| transverse(&f.data_type, codificacao)).collect();
        let arq = File::create(&self.caminho).map_err(|e| erro(format!("Error creating {}: {e}", self.caminho)))?;
        let w = FileWriter::try_new(BufWriter::with_capacity(1 << 20, arq), esquema.clone(), opcoes)
            .map_err(|e| erro(format!("{}: {e}", self.caminho)))?;
        self.aberto = Some(Aberto {
            formato,
            w,
            esquema,
            codificacoes,
            opcoes,
            nomes,
            colunas: tipos.iter().map(|t| Coluna::nova(*t)).collect(),
            no_grupo: 0,
            bytes_no_grupo: 0,
        });
        Ok(())
    }

    fn linha(&mut self, bruta: &[u8]) -> io::Result<()> {
        self.linha_n += 1;
        let n = self.linha_n;
        let bruta = match bruta.last() {
            Some(b'\r') => &bruta[..bruta.len() - 1],
            _ => bruta,
        };
        if bruta.is_empty() {
            return Ok(());
        }
        let linha = std::str::from_utf8(bruta).map_err(|_| erro(format!("line {n}: not valid text (UTF-8)")))?;
        if e_cabecalho(self.formato, linha) {
            self.cabecalho.push(linha.to_string());
            return Ok(());
        }
        if self.aberto.is_none() {
            self.abre(Some(linha))?;
        }
        let Some(a) = self.aberto.as_mut() else {
            return Err(erro("Parquet output was not opened".into()));
        };
        let total = a.colunas.len();
        let minimo = obrigatorias(a.formato);
        // No SAM os campos opcionais ficam todos na ultima coluna.
        let mut campos = if a.formato == Formato::Sam { linha.splitn(total, '\t') } else { linha.splitn(usize::MAX, '\t') };
        let mut i = 0;
        for campo in campos.by_ref() {
            if i >= total {
                return Err(erro(format!(
                    "line {n}: {} columns, but the {} has {total}; a Parquet file needs the same columns in every row",
                    linha.split('\t').count(),
                    if a.formato == Formato::Vcf { "header (or the first record)" } else { "first record" }
                )));
            }
            if a.colunas[i].poe(campo).is_err() {
                let visto: String = campo.chars().take(40).collect();
                return Err(erro(format!("line {n}: column {} ({}) should be a number, found '{visto}'", i + 1, a.nomes[i])));
            }
            i += 1;
        }
        if i < minimo {
            return Err(erro(format!("line {n}: {i} columns, but a record needs at least {minimo}")));
        }
        for c in &mut a.colunas[i..] {
            c.nulo();
        }
        a.no_grupo += 1;
        a.bytes_no_grupo += linha.len();
        self.registros += 1;
        if a.no_grupo >= LINHAS_POR_GRUPO || a.bytes_no_grupo >= BYTES_POR_GRUPO {
            self.fecha_grupo()?;
        }
        Ok(())
    }

    fn fecha_grupo(&mut self) -> io::Result<()> {
        let Some(a) = self.aberto.as_mut() else { return Ok(()) };
        if a.no_grupo == 0 {
            return Ok(());
        }
        let arrays: Vec<Box<dyn Array>> = a.colunas.iter_mut().map(Coluna::fecha).collect();
        let falha = |e: polars::error::PolarsError| erro(format!("Error writing Parquet: {e}"));
        let grupos = RowGroupIterator::try_new(std::iter::once(Ok(Chunk::new(arrays))), &a.esquema, a.opcoes, a.codificacoes.clone())
            .map_err(falha)?;
        for g in grupos {
            a.w.write(g.map_err(falha)?).map_err(falha)?;
        }
        a.no_grupo = 0;
        a.bytes_no_grupo = 0;
        Ok(())
    }

    /// Grava o que falta e o rodape, com o cabecalho nos metadados.
    pub fn termina(mut self) -> io::Result<u64> {
        if !self.resto.is_empty() {
            let resto = std::mem::take(&mut self.resto);
            self.linha(&resto)?;
        }
        if self.aberto.is_none() {
            self.abre(None)?;
        }
        self.fecha_grupo()?;
        let cabecalho = self.cabecalho.join("\n");
        let Some(a) = self.aberto.as_mut() else {
            return Err(erro("Parquet output was not opened".into()));
        };
        let metadados = vec![
            KeyValue::new(CHAVE_FORMATO.to_string(), nome_formato(a.formato).to_string()),
            KeyValue::new(CHAVE_CABECALHO.to_string(), cabecalho),
        ];
        a.w.end(Some(metadados)).map_err(|e| erro(format!("Error writing Parquet: {e}")))?;
        Ok(self.registros)
    }
}

impl Write for ParquetEscritor {
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

// ---------------------------------------------------------------------
//  LEITURA
// ---------------------------------------------------------------------

/// Uma coluna de um grupo, ja levada a um de dois tipos: texto, ou numero de
/// ponto flutuante (que precisa ser impresso com cuidado).
enum Lida {
    Texto(Utf8Array<i64>),
    Real(PrimitiveArray<f64>),
}

/// Imprime um numero sem perder nada e sem enfeite: 30 e nao 30.0.
fn escreve_real(v: f64, s: &mut Vec<u8>) {
    let a = v.abs();
    let t = if v.is_finite() && a != 0.0 && !(1e-5..1e16).contains(&a) { format!("{v:e}") } else { format!("{v}") };
    s.extend_from_slice(t.as_bytes());
}

impl Lida {
    fn de(nome: &str, array: &Box<dyn Array>, caminho: &str) -> io::Result<Lida> {
        // Caminho curto: e o que este mesmo modulo grava.
        if let Some(a) = array.as_any().downcast_ref::<Utf8Array<i64>>() {
            return Ok(Lida::Texto(a.clone()));
        }
        if let Some(a) = array.as_any().downcast_ref::<PrimitiveArray<f64>>() {
            return Ok(Lida::Real(a.clone()));
        }
        // Qualquer outro tipo (inteiros, dicionario, texto de 32 bits...)
        // passa pelo Polars.
        let falha = |e: polars::error::PolarsError| erro(format!("{caminho}: cannot read column '{nome}': {e}"));
        let serie = Series::try_from((nome, array.clone())).map_err(falha)?;
        if serie.dtype().is_float() {
            let s = serie.cast(&DataType::Float64).map_err(falha)?.rechunk();
            let a = s.f64().map_err(falha)?.downcast_iter().next().cloned();
            return Ok(Lida::Real(a.unwrap_or_else(|| PrimitiveArray::<f64>::from_vec(Vec::new()))));
        }
        let s = serie.cast(&DataType::String).map_err(falha)?.rechunk();
        let a = s.str().map_err(falha)?.downcast_iter().next().cloned();
        Ok(Lida::Texto(a.unwrap_or_else(|| Utf8Array::<i64>::from_slice::<&str, _>([]))))
    }

    fn nulo(&self, i: usize) -> bool {
        match self {
            Lida::Texto(a) => a.is_null(i),
            Lida::Real(a) => a.is_null(i),
        }
    }

    fn escreve(&self, i: usize, s: &mut Vec<u8>) {
        match self {
            Lida::Texto(a) => s.extend_from_slice(a.value(i).as_bytes()),
            Lida::Real(a) => escreve_real(a.value(i), s),
        }
    }
}

/// Nomes aceitos, em minusculas, para cada coluna fixa quando o arquivo veio
/// de outra ferramenta.
fn apelidos(formato: Formato, i: usize) -> &'static [&'static str] {
    match (formato, i) {
        (Formato::Vcf, 0) => &["chrom", "#chrom", "chr", "chromosome", "contig"],
        (Formato::Vcf, 1) => &["pos", "position"],
        (Formato::Vcf, 3) => &["ref", "reference"],
        (Formato::Vcf, 4) => &["alt", "alternate"],
        (Formato::Bed, 0) => &["chrom", "#chrom", "chr", "chromosome", "seqname", "seqnames", "contig"],
        (Formato::Bed, 1) => &["start", "chromstart"],
        (Formato::Bed, 2) => &["end", "chromend", "stop"],
        (Formato::Gff, 0) => &["seqid", "seqname", "seqnames", "chrom", "chr"],
        (Formato::Gff, 2) => &["type", "feature"],
        (Formato::Gff, 5) => &["score"],
        (Formato::Gff, 7) => &["phase", "frame"],
        (Formato::Gff, 8) => &["attributes", "attribute", "group"],
        _ => &[],
    }
}

/// O tipo em que uma coluna e decodificada. O decodificador do Polars 0.36 so
/// sabe entregar texto com offsets de 64 bits; o pyarrow e o pandas gravam os
/// de 32 e, muitas vezes, em dicionario. O proprio Polars faz esta mesma troca
/// antes de ler (`enlarge_data_type`).
fn alarga(tipo: ArrowDataType) -> ArrowDataType {
    match tipo {
        ArrowDataType::Utf8 => ArrowDataType::LargeUtf8,
        ArrowDataType::Binary => ArrowDataType::LargeBinary,
        ArrowDataType::Dictionary(_, valores, _) => alarga(*valores),
        outro => outro,
    }
}

/// O esquema Arrow do arquivo. Quem grava com uma versao do Arrow mais nova que
/// a deste leitor (tipos `Utf8View`, por exemplo) deixa nos metadados um
/// esquema que ele nao entende; nesse caso vale o que o proprio Parquet diz
/// das colunas.
fn esquema_do_arquivo(meta: &polars_parquet::read::FileMetaData) -> Result<ArrowSchema, String> {
    let mut esquema = match infer_schema(meta) {
        Ok(e) => e,
        Err(_) => {
            let mut sem = meta.clone();
            sem.key_value_metadata =
                sem.key_value_metadata.map(|kv| kv.into_iter().filter(|k| k.key != "ARROW:schema").collect());
            infer_schema(&sem).map_err(|e| e.to_string())?
        }
    };
    for f in &mut esquema.fields {
        f.data_type = alarga(f.data_type.clone());
    }
    Ok(esquema)
}

/// Le Parquet e entrega as linhas de texto do formato.
pub struct ParquetLeitor {
    leitor: FileReader<BufReader<File>>,
    caminho: String,
    formato: Formato,
    nomes: Vec<String>,
    /// Para cada coluna da linha de texto, a coluna do arquivo (ou nenhuma).
    mapa: Vec<Option<usize>>,
    texto: Vec<u8>,
    pos: usize,
    linha_n: u64,
    acabou: bool,
}

impl ParquetLeitor {
    pub fn abre(caminho: &str) -> Result<ParquetLeitor, String> {
        let arq = File::open(caminho).map_err(|e| format!("Error opening {caminho}: {e}"))?;
        let mut r = BufReader::with_capacity(1 << 20, arq);
        let meta = read_metadata(&mut r).map_err(|e| format!("{caminho}: not a readable Parquet file: {e}"))?;
        let esquema = esquema_do_arquivo(&meta).map_err(|e| format!("{caminho}: {e}"))?;
        let valor = |chave: &str| -> Option<String> {
            meta.key_value_metadata.as_ref()?.iter().find(|k| k.key == chave).and_then(|k| k.value.clone())
        };
        let nomes: Vec<String> = esquema.fields.iter().map(|f| f.name.clone()).collect();
        let minusculos: Vec<String> = nomes.iter().map(|n| n.to_ascii_lowercase()).collect();
        let tem = |n: &str| minusculos.iter().any(|m| m == n);
        let acha = |formato: Formato, i: usize| -> Option<usize> {
            let padrao = colunas_fixas(formato)[i].0.to_ascii_lowercase();
            minusculos.iter().position(|m| *m == padrao || apelidos(formato, i).contains(&m.as_str()))
        };

        let proprio = valor(CHAVE_FORMATO).and_then(|f| Formato::de_texto(&f));
        let formato = proprio
            .or_else(|| formato_pelo_nome(caminho))
            .or_else(|| {
                if acha(Formato::Vcf, 0).is_some() && tem("pos") && tem("ref") && tem("alt") {
                    Some(Formato::Vcf)
                } else if tem("qname") && tem("flag") && tem("rname") {
                    Some(Formato::Sam)
                } else if tem("attributes") || tem("attribute") {
                    Some(Formato::Gff)
                } else if acha(Formato::Bed, 0).is_some() && acha(Formato::Bed, 1).is_some() && acha(Formato::Bed, 2).is_some() {
                    Some(Formato::Bed)
                } else {
                    None
                }
            })
            .ok_or_else(|| {
                format!(
                    "{caminho}: cannot tell which genomic format these columns are ({}). Expected CHROM, POS, REF, ALT (VCF); chrom, start, end (BED); seqid, start, end, attributes (GFF); or QNAME, FLAG, RNAME... (SAM)",
                    nomes.join(", ")
                )
            })?;

        // Que coluna do arquivo alimenta cada coluna da linha.
        let fixas = colunas_fixas(formato);
        let mut mapa: Vec<Option<usize>> = Vec::new();
        if proprio.is_some() {
            // Gravado por este modulo: as colunas estao na ordem da linha.
            mapa.extend((0..nomes.len()).map(Some));
        } else {
            let n_fixas = if formato == Formato::Bed { 3 } else { fixas.len() };
            for i in 0..n_fixas {
                mapa.push(acha(formato, i));
            }
            let usadas: Vec<usize> = mapa.iter().flatten().copied().collect();
            let sobram = (0..nomes.len()).filter(|i| !usadas.contains(i));
            match formato {
                // FORMAT primeiro, depois as amostras na ordem do arquivo.
                Formato::Vcf => {
                    if let Some(f) = minusculos.iter().position(|m| m == "format") {
                        mapa.push(Some(f));
                        mapa.extend(sobram.filter(|&i| i != f).map(Some));
                    }
                }
                Formato::Bed => mapa.extend(sobram.map(Some)),
                _ => {}
            }
            for i in 0..obrigatorias(formato).min(n_fixas) {
                if mapa[i].is_none() {
                    return Err(format!(
                        "{caminho}: column {} is missing; a {} table needs it. Columns found: {}",
                        fixas[i].0,
                        nome_formato(formato).to_ascii_uppercase(),
                        nomes.join(", ")
                    ));
                }
            }
        }

        // O cabecalho: o guardado nos metadados ou, sem ele, o minimo.
        let mut texto = Vec::new();
        match valor(CHAVE_CABECALHO) {
            Some(c) => {
                if !c.is_empty() {
                    texto.extend_from_slice(c.as_bytes());
                    texto.push(b'\n');
                }
            }
            None => match formato {
                Formato::Vcf => {
                    texto.extend_from_slice(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
                    for m in mapa.iter().skip(8).flatten() {
                        texto.push(b'\t');
                        texto.extend_from_slice(nomes[*m].as_bytes());
                    }
                    texto.push(b'\n');
                }
                Formato::Gff => texto.extend_from_slice(b"##gff-version 3\n"),
                _ => {}
            },
        }

        let leitor = FileReader::new(r, meta.row_groups, esquema, Some(1 << 16), None, None);
        Ok(ParquetLeitor { leitor, caminho: caminho.to_string(), formato, nomes, mapa, texto, pos: 0, linha_n: 0, acabou: false })
    }

    /// O formato das linhas que este arquivo entrega.
    pub fn formato(&self) -> Formato {
        self.formato
    }

    /// Converte o proximo grupo de linhas. `false` quando nao ha mais.
    fn enche(&mut self) -> io::Result<bool> {
        self.texto.clear();
        self.pos = 0;
        while !self.acabou && self.texto.is_empty() {
            let Some(grupo) = self.leitor.next() else {
                self.acabou = true;
                break;
            };
            let grupo = grupo.map_err(|e| erro(format!("{}: {e}", self.caminho)))?;
            let arrays = grupo.arrays();
            let mut colunas: Vec<Option<Lida>> = Vec::with_capacity(self.mapa.len());
            for m in &self.mapa {
                colunas.push(match m {
                    Some(i) => Some(Lida::de(&self.nomes[*i], &arrays[*i], &self.caminho)?),
                    None => None,
                });
            }
            let minimo = minimo_na_leitura(self.formato);
            let obrig = obrigatorias(self.formato).min(colunas.len());
            for linha in 0..grupo.len() {
                self.linha_n += 1;
                // A linha vai ate a ultima coluna com valor, e nunca fica
                // menor que o minimo do formato.
                let ultimo = (0..colunas.len()).rev().find(|&c| colunas[c].as_ref().is_some_and(|l| !l.nulo(linha)));
                let n = ultimo.map_or(0, |u| u + 1).max(minimo.min(colunas.len()));
                for (c, coluna) in colunas.iter().enumerate().take(n) {
                    if c > 0 {
                        self.texto.push(b'\t');
                    }
                    match coluna {
                        Some(l) if !l.nulo(linha) => l.escreve(linha, &mut self.texto),
                        _ if c < obrig => {
                            return Err(erro(format!(
                                "{}, row {}: column {} is empty; a record cannot be written without it",
                                self.caminho,
                                self.linha_n,
                                colunas_fixas(self.formato)[c].0
                            )));
                        }
                        // Os campos opcionais do SAM, ausentes, nao deixam nada.
                        _ if self.formato == Formato::Sam && c == 11 => {
                            self.texto.pop();
                        }
                        _ => self.texto.push(b'.'),
                    }
                }
                self.texto.push(b'\n');
            }
        }
        Ok(!self.texto.is_empty())
    }
}

impl Read for ParquetLeitor {
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

/// O formato das linhas guardadas num `.parquet`, sem ler os dados.
pub fn formato_do_arquivo(caminho: &str) -> Result<Formato, String> {
    Ok(ParquetLeitor::abre(caminho)?.formato())
}
