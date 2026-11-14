// =====================================================================
//  PARENTESCO, TRIOS, CORRIDAS DE HOMOZIGOSE E SEXO GENETICO
//
//  Quatro analises sobre os genotipos (GT) de um ou mais VCFs:
//
//   - `relatedness`           todos os pares de pessoas: HETHET, IBS0, IBS1,
//                             IBS2, o parentesco KING-robust e a relacao
//                             inferida (duplicata/gemeo MZ, pai-filho, irmaos,
//                             2o e 3o grau, nao aparentados);
//   - `trio_analysis`         filho, pai e (opcional) mae: erros mendelianos
//                             por tipo, homozigotos opostos (exclusoes),
//                             candidatos a mutacao de novo e o indice de
//                             paternidade combinado (CPI);
//   - `runs_of_homozygosity`  trechos sem nenhum heterozigoto (ROH) e F_ROH;
//   - `genetic_sex`           heterozigose do X fora das PARs contra a dos
//                             autossomos, mais as chamadas no Y.
//
//  Regras da casa:
//
//   - EM FLUXO. Os VCFs sao lidos uma linha por vez e juntados por um merge
//     k-vias em (cromossomo, posicao). A memoria nao depende do tamanho dos
//     VCFs: so dos BEDs de regioes chamaveis (que ficam na memoria, um vetor
//     de intervalos por cromossomo) e do numero de pares.
//   - NOMES. `chr1` e `1` sao o mesmo cromossomo, `chrM` e `MT` tambem
//     (`nomes::canonico`). A ordem esperada e a natural (1, 2, ..., 22, X, Y),
//     a mesma do `sort_genomic(..., chrom_order="natural")`. Sequencias que nao
//     sao cromossomos primarios (scaffolds, alt, decoys) e a mitocondria sao
//     puladas -- e contadas -- antes de qualquer conferencia de ordem: a ordem
//     delas varia de uma referencia para outra (o hg19 poe o chrM primeiro).
//   - UMA PESSOA = (arquivo, amostra, BED de regioes chamaveis). Um VCF com
//     varias amostras da varias pessoas, lidas numa passada so.
//   - CHAVE = (cromossomo, posicao, REF, ALT), REF e ALT sem diferenciar
//     maiusculas. Registros multialelicos e alelos simbolicos nao viram chave.
//   - SITIO AUSENTE do arquivo de uma pessoa: homozigoto da referencia so
//     DENTRO do BED dela (e assim que os benchmarks do GIAB devem ser lidos);
//     fora do BED a pessoa nao tem dado, nem onde o arquivo tem registro. Sem
//     BED, decide `absent_is_ref`. Registro na mesma posicao com outro ALT:
//     genotipo faltante para aquela chave.
//   - Parentesco e paternidade usam so SNVs bialelicos dos autossomos; as
//     contagens mendelianas usam todo sitio bialelico (SNV, indel, MNP).
//   - Nenhum texto faz afirmacao medica ou legal: e software de ensino e
//     pesquisa, e os relatorios dizem isso.
// =====================================================================

use crate::externo::{abre_leitura, natural_cmp, Escritor, Formato};
use crate::nomes::{canonico, json_texto, parece_cromossomo};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::io::BufRead;

// ---------------------------------------------------------------------
//  CONSTANTES
// ---------------------------------------------------------------------

/// Limiares de parentesco do KING (Manichaikul et al. 2010, Bioinformatics
/// 26:2867): 2^-1.5, 2^-2.5, 2^-3.5 e 2^-4.5, os pontos medios em escala log
/// entre os parentescos esperados de cada grau (1/2, 1/4, 1/8, 1/16, 1/32).
pub const KING_DUPLICATE: f64 = 0.354;
pub const KING_FIRST: f64 = 0.177;
pub const KING_SECOND: f64 = 0.0884;
pub const KING_THIRD: f64 = 0.0442;

/// Pai-filho x irmaos (os dois tem parentesco 1/4). Pai e filho compartilham
/// um alelo em TODO sitio, entao IBS0 (homozigotos opostos) so aparece por
/// erro de genotipagem; irmaos nao compartilham nada em 1/4 do genoma.
///
/// A razao IBS0/HETHET nao depende de quantos sitios ha nem do painel de
/// forma grosseira: para dois nao aparentados vale 1/2 (sum 2p^2q^2 sobre
/// sum 4p^2q^2); para irmaos vale S/(2(S+T)), com S = sum p^2q^2 e
/// T = sum pq, que da ~0.07 no espectro neutro (S/T = 1/6). Medido nos dados
/// reais (S e T tirados de pares nao aparentados): 0.075 no GIAB (HG003 x
/// HG004, genoma todo) e 0.069 no chr22 do 1000G (NA12891 x NA12892). Para
/// pai-filho e da ordem da taxa de erro: 0 a 3e-4 nos pares do GIAB, ate
/// 0.016 no chr22 do 1000G (NA12878 x NA12892). O limiar fica no meio, em
/// escala log, dos dois lados.
pub const PO_IBS0_RATIO: f64 = 0.03;

/// Pai e filho compartilham exatamente um haplotipo em todo ponto: o
/// parentesco e 1/4 (0.242 a 0.273 nos pares reais medidos, inclusive num
/// cromossomo so). Bem acima disso o par compartilha os dois haplotipos em
/// algum trecho (IBD2), o que irmaos fazem e pai e filho nao (fora de
/// familias consanguineas). Importa com pouco genoma: no chr22 sozinho, dois
/// irmaos sem trecho IBD0 tem IBS0 ~ 0 e parentesco 0.30-0.36.
pub const PO_MAX_KINSHIP: f64 = 0.30;

/// Com menos sitios informativos que isto (o denominador do KING: sitios em
/// que alguem do par e heterozigoto), a relacao nao e inferida.
pub const MIN_KINSHIP_SITES: u64 = 100;

/// Menos sitios informativos para exclusao que isto: veredito inconclusivo.
pub const MIN_INFORMATIVE_SITES: u64 = 100;

/// Frequencias alelicas sao limitadas a [MIN_AF, 1 - MIN_AF] no CPI: um
/// alelo que o painel nao viu nao pode valer um indice infinito.
pub const MIN_AF: f64 = 0.001;

/// Sexo genetico: razao (fracao de heterozigotos no X fora das PARs) /
/// (fracao nos autossomos). Abaixo de ONE_X, um X; acima de TWO_X, dois.
pub const ONE_X_RATIO: f64 = 0.2;
pub const TWO_X_RATIO: f64 = 0.5;
/// Chamadas homozigotas/haploides nao-referencia no Y (fora das PARs) a
/// partir das quais se considera que ha um Y.
pub const MIN_Y_CALLS: u64 = 20;
/// Chamadas nao-referencia no X fora das PARs abaixo disto: sem dado de X.
pub const MIN_X_CALLS: u64 = 100;
/// Fracao de heterozigotos tipica nos autossomos, usada quando o arquivo nao
/// tem autossomos suficientes para medir.
const AUTOSOME_HET_FRACTION_FALLBACK: f64 = 0.6;

/// Tamanho dos autossomos 1..22 do GRCh38 (soma 2.875.001.522 pb).
pub const GRCH38_AUTOSOMES: [u64; 22] = [
    248_956_422, 242_193_529, 198_295_559, 190_214_555, 181_538_259, 170_805_979, 159_345_973, 145_138_636,
    138_394_717, 133_797_422, 135_086_622, 133_275_309, 114_364_328, 107_043_718, 101_991_189, 90_338_345,
    83_257_441, 80_373_285, 58_617_616, 64_444_167, 46_709_983, 50_818_468,
];
/// Tamanho dos autossomos 1..22 do GRCh37 (soma 2.881.033.286 pb).
pub const GRCH37_AUTOSOMES: [u64; 22] = [
    249_250_621, 243_199_373, 198_022_430, 191_154_276, 180_915_260, 171_115_067, 159_138_663, 146_364_022,
    141_213_431, 135_534_747, 135_006_516, 133_851_895, 115_169_878, 107_349_540, 102_531_392, 90_354_753,
    81_195_210, 78_077_248, 59_128_983, 63_025_520, 48_129_895, 51_304_566,
];

pub fn grch38_autosome_length() -> u64 {
    GRCH38_AUTOSOMES.iter().sum()
}

pub fn grch37_autosome_length() -> u64 {
    GRCH37_AUTOSOMES.iter().sum()
}

/// Regioes pseudoautossomicas (1-based, fim incluido).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Pars {
    x: [(u64, u64); 2],
    y: [(u64, u64); 2],
}

const PAR_GRCH38: Pars = Pars {
    x: [(10_001, 2_781_479), (155_701_383, 156_030_895)],
    y: [(10_001, 2_781_479), (56_887_903, 57_217_415)],
};
const PAR_GRCH37: Pars = Pars {
    x: [(60_001, 2_699_520), (154_931_044, 155_260_560)],
    y: [(10_001, 2_649_520), (59_034_050, 59_363_566)],
};
const CHRX_GRCH38: u64 = 156_040_895;
const CHRX_GRCH37: u64 = 155_270_560;

fn na_par(faixas: &[(u64, u64); 2], pos: u64) -> bool {
    faixas.iter().any(|&(a, b)| pos >= a && pos <= b)
}

// ---------------------------------------------------------------------
//  TIPOS PUBLICOS: PESSOAS E OPCOES
// ---------------------------------------------------------------------

/// One person: a VCF, which of its samples, and optionally the BED of the
/// regions where that person's calls can be trusted.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Person {
    pub path: String,
    /// The sample column. `None`: the file's only sample (an error if the
    /// file has several).
    pub sample: Option<String>,
    /// Callable regions. Inside them, a site absent from the file is
    /// homozygous reference; outside them the person has no data, even where
    /// the file has a record.
    pub callable_bed: Option<String>,
    /// Name in the reports. `None`: the sample name.
    pub label: Option<String>,
}

impl Person {
    pub fn new(path: &str) -> Person {
        Person { path: path.to_string(), ..Default::default() }
    }
    pub fn sample(mut self, s: &str) -> Person {
        self.sample = Some(s.to_string());
        self
    }
    pub fn bed(mut self, b: &str) -> Person {
        self.callable_bed = Some(b.to_string());
        self
    }
    pub fn label(mut self, l: &str) -> Person {
        self.label = Some(l.to_string());
        self
    }
}

/// How the files are joined.
#[derive(Clone, Debug, PartialEq)]
pub struct JoinOptions {
    /// For a person WITHOUT a callable BED: is a site absent from the file
    /// homozygous reference (true, the default: single-sample VCFs list only
    /// variants) or missing (false: the file is not trusted to be complete)?
    pub absent_is_ref: bool,
    /// Use only records whose FILTER is `PASS` or `.`; a person's genotype
    /// in any other record counts as missing. Default false.
    pub pass_only: bool,
}

impl Default for JoinOptions {
    fn default() -> Self {
        JoinOptions { absent_is_ref: true, pass_only: false }
    }
}

// ---------------------------------------------------------------------
//  UTILIDADES
// ---------------------------------------------------------------------

fn json_num(x: f64) -> String {
    if x.is_finite() {
        format!("{x}")
    } else {
        "null".to_string()
    }
}

fn json_opt(x: Option<f64>) -> String {
    x.map(json_num).unwrap_or_else(|| "null".to_string())
}

fn json_lista_u64(v: &[u64]) -> String {
    let t: Vec<String> = v.iter().map(|x| x.to_string()).collect();
    format!("[{}]", t.join(","))
}

/// 12345678 -> "12,345,678"
fn milhar(n: u64) -> String {
    let s = n.to_string();
    let mut o = String::with_capacity(s.len() + s.len() / 3);
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            o.push(',');
        }
        o.push(c);
    }
    o
}

fn fmt_opt(x: Option<f64>, casas: usize) -> String {
    match x {
        Some(v) if v.is_finite() => format!("{v:.casas$}"),
        _ => "n/a".to_string(),
    }
}

fn erro_leitura(caminho: &str, linha: usize, e: std::io::Error) -> String {
    format!(
        "Error reading {caminho} after line {linha}: {e}. The file may be truncated, binary (BCF is not read; convert it to VCF first) or not UTF-8 text."
    )
}

/// Non-empty and made only of A, C, G, T, N in either case.
fn so_bases(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}

fn e_acgt(s: &str) -> bool {
    s.bytes().all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
}

fn caminho_canonico(p: &str) -> String {
    std::fs::canonicalize(p).map(|c| c.to_string_lossy().to_string()).unwrap_or_else(|_| p.to_string())
}

/// Removes a half-written output when the operation ends in error.
struct Saida {
    caminho: Option<String>,
}

impl Saida {
    fn desarma(&mut self) {
        self.caminho = None;
    }
}

impl Drop for Saida {
    fn drop(&mut self) {
        if let Some(c) = &self.caminho {
            let _ = std::fs::remove_file(c);
        }
    }
}

fn recusa_sobrescrever_entrada(saida: &str, entradas: &[&str]) -> Result<(), String> {
    let s = caminho_canonico(saida);
    for e in entradas {
        if *e == saida || caminho_canonico(e) == s {
            return Err(format!(
                "The output `{saida}` is one of the input files. Writing would destroy it while it is being read: choose a different output path."
            ));
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------
//  GENOTIPOS
// ---------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
enum Gt {
    #[default]
    Falta,
    Hap(u16),
    Dip(u16, u16),
}

/// Reads a GT value: `0/1`, `1|0`, `1` (haploid). Anything with a `.`
/// allele, more than two alleles or junk is missing.
fn le_gt(s: &str) -> Gt {
    let mut alelos = [0u16; 2];
    let mut n = 0usize;
    let mut atual: Option<u32> = None;
    for c in s.bytes() {
        match c {
            b'0'..=b'9' => {
                let v = atual.unwrap_or(0) * 10 + u32::from(c - b'0');
                if v > u32::from(u16::MAX) {
                    return Gt::Falta;
                }
                atual = Some(v);
            }
            b'/' | b'|' => {
                let Some(v) = atual else { return Gt::Falta };
                if n == 2 {
                    return Gt::Falta;
                }
                alelos[n] = v as u16;
                n += 1;
                atual = None;
            }
            _ => return Gt::Falta,
        }
    }
    let Some(v) = atual else { return Gt::Falta };
    if n == 2 {
        return Gt::Falta;
    }
    alelos[n] = v as u16;
    match n + 1 {
        1 => Gt::Hap(alelos[0]),
        _ => Gt::Dip(alelos[0], alelos[1]),
    }
}

/// Number of ALT alleles in a diploid biallelic genotype; -1 when missing,
/// haploid or with an allele index above 1.
fn codigo(g: Gt) -> i8 {
    match g {
        Gt::Dip(a, b) if a <= 1 && b <= 1 => (a + b) as i8,
        _ => -1,
    }
}

// ---------------------------------------------------------------------
//  CROMOSSOMOS
// ---------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TipoSeq {
    Autossomo,
    X,
    Y,
}

/// O papel de um nome canonico. `None`: nao e cromossomo primario (scaffold,
/// contig alternativo, decoy) ou e a mitocondria -- fica de fora.
fn tipo_seq(canon: &str) -> Option<TipoSeq> {
    match canon {
        "X" => Some(TipoSeq::X),
        "Y" => Some(TipoSeq::Y),
        "MT" | "M" | "W" | "Z" => None,
        c if parece_cromossomo(c) => Some(TipoSeq::Autossomo),
        _ => None,
    }
}

#[derive(Default)]
struct Cromos {
    ids: HashMap<String, u32>,
    nomes: Vec<String>,
    tipos: Vec<TipoSeq>,
}

impl Cromos {
    fn id(&mut self, canon: &str, tipo: TipoSeq) -> u32 {
        if let Some(&i) = self.ids.get(canon) {
            return i;
        }
        let i = self.nomes.len() as u32;
        self.ids.insert(canon.to_string(), i);
        self.nomes.push(canon.to_string());
        self.tipos.push(tipo);
        i
    }

    fn compara(&self, a: u32, b: u32) -> Ordering {
        if a == b {
            Ordering::Equal
        } else {
            natural_cmp(&self.nomes[a as usize], &self.nomes[b as usize])
        }
    }

    fn tipo(&self, id: u32) -> TipoSeq {
        self.tipos[id as usize]
    }
}

fn compara_pos(c: &Cromos, a: (u32, u64), b: (u32, u64)) -> Ordering {
    c.compara(a.0, b.0).then(a.1.cmp(&b.1))
}

// ---------------------------------------------------------------------
//  REGIOES CHAMAVEIS (BED)
// ---------------------------------------------------------------------

struct Chamavel {
    /// Intervalos 0-based semiabertos, ordenados e fundidos, por cromossomo.
    por_crom: HashMap<u32, Vec<(u64, u64)>>,
    crom: u32,
    idx: usize,
    ultimo: u64,
    valido: bool,
}

impl Chamavel {
    fn carrega(caminho: &str, cromos: &mut Cromos) -> Result<Chamavel, String> {
        let mut r = abre_leitura(caminho)?;
        let mut por_crom: HashMap<u32, Vec<(u64, u64)>> = HashMap::new();
        let mut l = String::new();
        let mut n = 0usize;
        let mut cache: (String, Option<u32>) = (String::new(), None);
        loop {
            l.clear();
            if r.read_line(&mut l).map_err(|e| erro_leitura(caminho, n, e))? == 0 {
                break;
            }
            n += 1;
            let t = l.trim_end_matches(['\n', '\r']);
            if t.is_empty() || t.starts_with('#') || t.starts_with("track") || t.starts_with("browser") {
                continue;
            }
            let mut c = t.split('\t');
            let (Some(nome), Some(a), Some(b)) = (c.next(), c.next(), c.next()) else {
                return Err(format!(
                    "{caminho}, line {n}: a BED line needs at least three tab-separated columns (chromosome, start, end)."
                ));
            };
            let (Ok(a), Ok(b)) = (a.trim().parse::<u64>(), b.trim().parse::<u64>()) else {
                return Err(format!("{caminho}, line {n}: start and end must be whole numbers (0-based start, end excluded)."));
            };
            if b < a {
                return Err(format!("{caminho}, line {n}: the end ({b}) is before the start ({a})."));
            }
            if nome != cache.0 {
                let canon = canonico(nome);
                cache = (nome.to_string(), tipo_seq(&canon).map(|tp| cromos.id(&canon, tp)));
            }
            if let Some(id) = cache.1 {
                if b > a {
                    por_crom.entry(id).or_default().push((a, b));
                }
            }
        }
        for v in por_crom.values_mut() {
            v.sort_unstable();
            let mut fundidos: Vec<(u64, u64)> = Vec::with_capacity(v.len());
            for &(a, b) in v.iter() {
                match fundidos.last_mut() {
                    Some(u) if a <= u.1 => u.1 = u.1.max(b),
                    _ => fundidos.push((a, b)),
                }
            }
            *v = fundidos;
        }
        if por_crom.is_empty() {
            return Err(format!(
                "{caminho}: no region on a primary chromosome. A callable-regions BED lists `chromosome start end` (for example `chr1 10000 20000`)."
            ));
        }
        Ok(Chamavel { por_crom, crom: 0, idx: 0, ultimo: 0, valido: false })
    }

    /// Is the 1-based `pos` inside a region? Fast when the positions of a
    /// chromosome come in order (a cursor walks the intervals).
    fn contem(&mut self, crom: u32, pos: u64) -> bool {
        let Some(v) = self.por_crom.get(&crom) else {
            return false;
        };
        let p0 = pos.saturating_sub(1);
        if !self.valido || crom != self.crom || p0 < self.ultimo {
            self.crom = crom;
            self.valido = true;
            self.idx = v.partition_point(|iv| iv.1 <= p0);
        }
        self.ultimo = p0;
        while self.idx < v.len() && v[self.idx].1 <= p0 {
            self.idx += 1;
        }
        self.idx < v.len() && v[self.idx].0 <= p0
    }

    /// The longest stretch outside the regions between two positions inside
    /// them (1-based, a <= b).
    fn maior_lacuna(&self, crom: u32, a: u64, b: u64) -> u64 {
        let Some(v) = self.por_crom.get(&crom) else {
            return b.saturating_sub(a);
        };
        let ia = v.partition_point(|iv| iv.1 <= a.saturating_sub(1));
        let ib = v.partition_point(|iv| iv.1 <= b.saturating_sub(1));
        let mut maior = 0;
        for k in ia..ib.min(v.len().saturating_sub(1)) {
            maior = maior.max(v[k + 1].0.saturating_sub(v[k].1));
        }
        maior
    }
}

// ---------------------------------------------------------------------
//  LEITURA DE UM VCF EM FLUXO
// ---------------------------------------------------------------------

#[derive(Default)]
struct Reg {
    linha: String,
    n_linha: usize,
    crom: u32,
    pos: u64,
    r: (usize, usize),
    a: (usize, usize),
    info: (usize, usize),
    /// Fim da coluna FORMAT (ou da INFO, se nao houver FORMAT).
    fim_fixas: usize,
    passa: bool,
    /// Um por pessoa servida por este fluxo (na ordem dos `slots`).
    gts: Vec<Gt>,
    cols: Vec<(usize, usize)>,
}

impl Reg {
    fn refs(&self) -> &str {
        &self.linha[self.r.0..self.r.1]
    }
    fn alts(&self) -> &str {
        &self.linha[self.a.0..self.a.1]
    }
    fn info(&self) -> &str {
        &self.linha[self.info.0..self.info.1]
    }
    fn crom_bruto(&self) -> &str {
        self.linha.split('\t').next().unwrap_or("")
    }
}

/// Um VCF aberto: cabecalho lido, um registro a frente (`cab`), conferindo a
/// ordem a cada linha.
struct Fluxo {
    caminho: String,
    r: Box<dyn BufRead>,
    n_linha: usize,
    primeira: Option<(String, usize)>,
    amostras: Vec<String>,
    /// Linhas de cabecalho `##`, na ordem.
    cabecalho: Vec<String>,
    contig_x: Option<u64>,
    contig_y: bool,
    /// Coluna (0-based na linha) de cada pessoa servida; `ordem`: os slots
    /// em ordem de coluna, para extrair tudo numa varredura da linha.
    colunas: Vec<usize>,
    ordem: Vec<usize>,
    crom_atual: Option<u32>,
    nome_atual: String,
    pos_atual: u64,
    bruto_ant: String,
    id_ant: Option<u32>,
    formato_ant: String,
    gt_ant: Option<usize>,
    cab: Option<Reg>,
    livres: Vec<Reg>,
    registros: u64,
    outras_seqs: u64,
    malformados: u64,
}

impl Fluxo {
    fn abre(caminho: &str) -> Result<Fluxo, String> {
        let r = abre_leitura(caminho)?;
        let mut f = Fluxo {
            caminho: caminho.to_string(),
            r,
            n_linha: 0,
            primeira: None,
            amostras: Vec::new(),
            cabecalho: Vec::new(),
            contig_x: None,
            contig_y: false,
            colunas: Vec::new(),
            ordem: Vec::new(),
            crom_atual: None,
            nome_atual: String::new(),
            pos_atual: 0,
            bruto_ant: String::from("\u{0}"),
            id_ant: None,
            formato_ant: String::from("\u{0}"),
            gt_ant: None,
            cab: None,
            livres: Vec::new(),
            registros: 0,
            outras_seqs: 0,
            malformados: 0,
        };
        let mut l = String::new();
        let mut viu_chrom = false;
        while f.le(&mut l)? {
            if l.is_empty() {
                continue;
            }
            if l.starts_with("##") {
                if let Some((id, len)) = contig(&l) {
                    match canonico(&id).as_str() {
                        "X" => f.contig_x = Some(len.unwrap_or(0)),
                        "Y" => f.contig_y = true,
                        _ => {}
                    }
                }
                f.cabecalho.push(l.clone());
                continue;
            }
            if l.starts_with('#') {
                if l.starts_with("#CHROM") {
                    f.amostras = l.split('\t').skip(9).map(str::to_string).collect();
                    viu_chrom = true;
                }
                continue;
            }
            if !viu_chrom {
                return Err(format!(
                    "{caminho}, line {}: a record comes before the `#CHROM` header line. This is not a VCF, or its header was cut off.",
                    f.n_linha
                ));
            }
            f.primeira = Some((std::mem::take(&mut l), f.n_linha));
            break;
        }
        if !viu_chrom {
            return Err(format!(
                "{caminho}: no `#CHROM` header line was found. This is not a VCF (or it is empty); validate_file(...) tells what the file is."
            ));
        }
        Ok(f)
    }

    /// Reads one line into `buf`, without its line ending. `false` at the end.
    fn le(&mut self, buf: &mut String) -> Result<bool, String> {
        buf.clear();
        let n = self.r.read_line(buf).map_err(|e| erro_leitura(&self.caminho, self.n_linha, e))?;
        if n == 0 {
            return Ok(false);
        }
        self.n_linha += 1;
        while buf.ends_with('\n') || buf.ends_with('\r') {
            buf.pop();
        }
        Ok(true)
    }

    /// Puts the next valid record on a primary chromosome in `cab` (`None` at
    /// the end of the file).
    fn avanca(&mut self, cromos: &mut Cromos) -> Result<(), String> {
        loop {
            let mut reg = self.livres.pop().unwrap_or_default();
            if let Some((linha, n)) = self.primeira.take() {
                reg.linha = linha;
                reg.n_linha = n;
            } else {
                if !self.le(&mut reg.linha)? {
                    self.livres.push(reg);
                    self.cab = None;
                    return Ok(());
                }
                reg.n_linha = self.n_linha;
            }
            if reg.linha.is_empty() || reg.linha.starts_with('#') {
                self.livres.push(reg);
                continue;
            }
            self.registros += 1;
            if self.analisa(&mut reg, cromos)? {
                self.cab = Some(reg);
                return Ok(());
            }
            self.livres.push(reg);
        }
    }

    fn analisa(&mut self, reg: &mut Reg, cromos: &mut Cromos) -> Result<bool, String> {
        let linha = reg.linha.as_str();
        let mut lim = [(0usize, 0usize); 9];
        let mut ncol = 0usize;
        let mut ini = 0usize;
        while ncol < 9 {
            match linha[ini..].find('\t') {
                Some(k) => {
                    lim[ncol] = (ini, ini + k);
                    ini += k + 1;
                    ncol += 1;
                }
                None => {
                    lim[ncol] = (ini, linha.len());
                    ncol += 1;
                    break;
                }
            }
        }
        if ncol < 8 {
            self.malformados += 1;
            return Ok(false);
        }
        let amostras_ini = if ncol == 9 && lim[8].1 < linha.len() { Some(lim[8].1 + 1) } else { None };

        let bruto = &linha[lim[0].0..lim[0].1];
        if bruto != self.bruto_ant {
            self.bruto_ant.clear();
            self.bruto_ant.push_str(bruto);
            let canon = canonico(bruto);
            self.id_ant = tipo_seq(&canon).map(|t| cromos.id(&canon, t));
        }
        let Some(id) = self.id_ant else {
            self.outras_seqs += 1;
            return Ok(false);
        };
        let Ok(pos) = linha[lim[1].0..lim[1].1].parse::<u64>() else {
            self.malformados += 1;
            return Ok(false);
        };

        match self.crom_atual {
            Some(c) if c == id => {
                if pos < self.pos_atual {
                    return Err(format!(
                        "{}, line {}: position {pos} on `{bruto}` comes after position {}: the file is not sorted by position. Sort it first: sort_genomic([\"{}\"], \"sorted.vcf.gz\", chrom_order=\"natural\").",
                        self.caminho, reg.n_linha, self.pos_atual, self.caminho
                    ));
                }
            }
            Some(c) => {
                if cromos.compara(id, c) != Ordering::Greater {
                    return Err(format!(
                        "{}, line {}: chromosome `{bruto}` comes after `{}`, but these analyses read files in natural chromosome order (1, 2, ..., 10, ..., 22, X, Y; `chr1` and `1` are the same) and every chromosome in one block. Sort the file first: sort_genomic([\"{}\"], \"sorted.vcf.gz\", chrom_order=\"natural\").",
                        self.caminho, reg.n_linha, self.nome_atual, self.caminho
                    ));
                }
            }
            None => {}
        }
        if self.crom_atual != Some(id) {
            self.nome_atual.clear();
            self.nome_atual.push_str(bruto);
        }
        self.crom_atual = Some(id);
        self.pos_atual = pos;

        reg.crom = id;
        reg.pos = pos;
        reg.r = lim[3];
        reg.a = lim[4];
        reg.info = lim[7];
        let filtro = &linha[lim[6].0..lim[6].1];
        reg.passa = filtro == "PASS" || filtro == ".";
        reg.fim_fixas = if ncol == 9 { lim[8].1 } else { lim[7].1 };
        reg.gts.clear();
        reg.cols.clear();
        let ns = self.colunas.len();
        if ns > 0 {
            reg.gts.resize(ns, Gt::Falta);
            reg.cols.resize(ns, (0, 0));
            let gi = if ncol == 9 {
                let fmt = &linha[lim[8].0..lim[8].1];
                if fmt != self.formato_ant {
                    self.formato_ant.clear();
                    self.formato_ant.push_str(fmt);
                    self.gt_ant = fmt.split(':').position(|k| k == "GT");
                }
                self.gt_ant
            } else {
                None
            };
            if let (Some(mut p), Some(gi)) = (amostras_ini, gi) {
                // Um laco de bytes, nao `find` por coluna: com milhares de
                // amostras de poucos bytes cada (`0|0`), o custo de preparar
                // cada busca dominava a leitura.
                let b = linha.as_bytes();
                let mut col = 9usize;
                for &s in &self.ordem {
                    let alvo = self.colunas[s];
                    while col < alvo && p < b.len() {
                        if b[p] == b'\t' {
                            col += 1;
                        }
                        p += 1;
                    }
                    if col < alvo {
                        break;
                    }
                    let mut fim = p;
                    while fim < b.len() && b[fim] != b'\t' {
                        fim += 1;
                    }
                    reg.cols[s] = (p, fim);
                    reg.gts[s] = linha[p..fim].split(':').nth(gi).map_or(Gt::Falta, le_gt);
                }
            }
        }
        Ok(true)
    }
}

/// `##contig=<ID=chrX,length=156040895>` -> ("chrX", Some(156040895)).
fn contig(linha: &str) -> Option<(String, Option<u64>)> {
    let r = linha.strip_prefix("##contig=<")?;
    let r = r.strip_suffix('>').unwrap_or(r);
    let mut id = None;
    let mut len = None;
    for parte in r.split(',') {
        if let Some(v) = parte.strip_prefix("ID=") {
            id = Some(v.to_string());
        } else if let Some(v) = parte.strip_prefix("length=") {
            len = v.parse().ok();
        }
    }
    Some((id?, len))
}

fn lista_amostras(a: &[String]) -> String {
    const N: usize = 8;
    let v: Vec<&str> = a.iter().take(N).map(String::as_str).collect();
    if a.len() > N {
        format!("Samples in the file: {}, and {} more (vcf_samples(...) lists them all).", v.join(", "), a.len() - N)
    } else {
        format!("Samples in the file: {}.", v.join(", "))
    }
}

fn resolve_amostra(f: &Fluxo, p: &Person) -> Result<usize, String> {
    if f.amostras.is_empty() {
        return Err(format!(
            "{} has no sample columns (it is a sites-only VCF): these analyses need each person's genotype (GT).",
            f.caminho
        ));
    }
    match &p.sample {
        Some(s) => f
            .amostras
            .iter()
            .position(|a| a == s)
            .ok_or_else(|| format!("Sample `{s}` is not in {}. {}", f.caminho, lista_amostras(&f.amostras))),
        None if f.amostras.len() == 1 => Ok(0),
        None => Err(format!(
            "{} has {} samples: say which one is the person (sample=\"...\"). {}",
            f.caminho,
            f.amostras.len(),
            lista_amostras(&f.amostras)
        )),
    }
}

/// The sample names of a VCF, from its `#CHROM` line.
pub fn vcf_samples(path: &str) -> Result<Vec<String>, String> {
    Ok(Fluxo::abre(path)?.amostras)
}

// ---------------------------------------------------------------------
//  O MERGE K-VIAS: SITIOS COM O GENOTIPO DE CADA PESSOA
// ---------------------------------------------------------------------

struct Membro {
    rotulo: String,
    caminho: String,
    amostra: String,
    fluxo: usize,
    slot: usize,
    bed: Option<Chamavel>,
    dentro: Option<bool>,
}

#[derive(Clone, Copy, Debug)]
struct Chave {
    fluxo: usize,
    reg: usize,
    snv: bool,
    af: Option<f64>,
}

/// What the scan read and left out.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ScanCounts {
    /// Data lines read from the persons' files.
    pub records: u64,
    /// Records on sequences that are not primary chromosomes (scaffolds,
    /// alternative contigs, decoys) or on the mitochondrion: left out.
    pub other_sequences: u64,
    /// Lines with fewer than 8 columns or a POS that is not a number.
    pub malformed: u64,
    /// Records with more than one ALT allele: not used as sites (a person
    /// with such a record is missing at the other sites of that position).
    pub multiallelic: u64,
    /// Records with a symbolic or non-ACGTN allele, or ALT equal to REF.
    pub not_plain: u64,
    /// Distinct biallelic sites (chromosome, position, REF, ALT) built.
    pub sites: u64,
    /// Data lines read from the allele-frequency file.
    pub af_records: u64,
}

impl ScanCounts {
    fn json(&self) -> String {
        format!(
            "{{\"records\":{},\"other_sequences\":{},\"malformed\":{},\"multiallelic\":{},\"not_plain\":{},\"sites\":{},\"af_records\":{}}}",
            self.records, self.other_sequences, self.malformed, self.multiallelic, self.not_plain, self.sites, self.af_records
        )
    }

    fn nota(&self) -> String {
        let mut o = Vec::new();
        if self.other_sequences > 0 {
            o.push(format!("{} record(s) on other sequences (scaffolds, mitochondrion) left out", milhar(self.other_sequences)));
        }
        if self.malformed > 0 {
            o.push(format!(
                "{} malformed line(s) left out (validate_file(...) shows them)",
                milhar(self.malformed)
            ));
        }
        if self.multiallelic > 0 {
            o.push(format!("{} multiallelic record(s) not used (vcf_normalize(...) splits them)", milhar(self.multiallelic)));
        }
        if o.is_empty() {
            String::new()
        } else {
            format!("  note: {}\n", o.join("; "))
        }
    }
}

struct Varredura {
    cromos: Cromos,
    fluxos: Vec<Fluxo>,
    membros: Vec<Membro>,
    grupo: Vec<Vec<Reg>>,
    af: Option<usize>,
    campo_af: String,
    absent_is_ref: bool,
    pass_only: bool,
    crom: u32,
    tipo: TipoSeq,
    pos: u64,
    chaves: Vec<Chave>,
    /// chaves.len() x membros.len(): -1 faltante, 0/1/2 copias do ALT.
    cods: Vec<i8>,
    /// Indice, em `grupo[fluxo do membro]`, do registro de onde veio o codigo.
    regs: Vec<Option<usize>>,
    cont: ScanCounts,
}

impl Varredura {
    fn nova(persons: &[Person], af: Option<&str>, campo_af: &str, opts: &JoinOptions) -> Result<Varredura, String> {
        let mut cromos = Cromos::default();
        let mut fluxos: Vec<Fluxo> = Vec::new();
        let mut chaves_caminho: Vec<String> = Vec::new();
        let mut membros: Vec<Membro> = Vec::new();
        for (i, p) in persons.iter().enumerate() {
            if p.path.is_empty() {
                return Err(format!("Person {} has an empty path: give the VCF file of each person.", i + 1));
            }
            let chave = caminho_canonico(&p.path);
            let s = match chaves_caminho.iter().position(|c| *c == chave) {
                Some(s) => s,
                None => {
                    fluxos.push(Fluxo::abre(&p.path)?);
                    chaves_caminho.push(chave);
                    fluxos.len() - 1
                }
            };
            let f = &mut fluxos[s];
            let a = resolve_amostra(f, p)?;
            let coluna = 9 + a;
            let slot = match f.colunas.iter().position(|&c| c == coluna) {
                Some(k) => k,
                None => {
                    f.colunas.push(coluna);
                    f.colunas.len() - 1
                }
            };
            let amostra = f.amostras[a].clone();
            let bed = match &p.callable_bed {
                Some(b) => Some(Chamavel::carrega(b, &mut cromos)?),
                None => None,
            };
            let base = p.label.clone().unwrap_or_else(|| amostra.clone());
            let mut rotulo = base.clone();
            let mut k = 2;
            while membros.iter().any(|m| m.rotulo == rotulo) {
                rotulo = format!("{base}_{k}");
                k += 1;
            }
            membros.push(Membro { rotulo, caminho: p.path.clone(), amostra, fluxo: s, slot, bed, dentro: None });
        }
        let af_idx = match af {
            Some(a) => {
                if campo_af.is_empty() {
                    return Err("Name the INFO field that holds the allele frequency (for example \"AF\" or \"EUR_AF\").".to_string());
                }
                fluxos.push(Fluxo::abre(a)?);
                Some(fluxos.len() - 1)
            }
            None => None,
        };
        for f in fluxos.iter_mut() {
            let cols = f.colunas.clone();
            f.ordem = (0..cols.len()).collect();
            f.ordem.sort_by_key(|&s| cols[s]);
            f.avanca(&mut cromos)?;
        }
        let n = fluxos.len();
        Ok(Varredura {
            cromos,
            fluxos,
            membros,
            grupo: (0..n).map(|_| Vec::new()).collect(),
            af: af_idx,
            campo_af: campo_af.to_string(),
            absent_is_ref: opts.absent_is_ref,
            pass_only: opts.pass_only,
            crom: 0,
            tipo: TipoSeq::Autossomo,
            pos: 0,
            chaves: Vec::new(),
            cods: Vec::new(),
            regs: Vec::new(),
            cont: ScanCounts::default(),
        })
    }

    /// Moves to the next position that has at least one biallelic site.
    /// `false` when every person's file has ended.
    fn proxima(&mut self) -> Result<bool, String> {
        loop {
            for (s, g) in self.grupo.iter_mut().enumerate() {
                self.fluxos[s].livres.append(g);
            }
            let mut menor: Option<(u32, u64)> = None;
            for (s, f) in self.fluxos.iter().enumerate() {
                if Some(s) == self.af {
                    continue;
                }
                if let Some(r) = &f.cab {
                    if menor.is_none_or(|m| compara_pos(&self.cromos, (r.crom, r.pos), m) == Ordering::Less) {
                        menor = Some((r.crom, r.pos));
                    }
                }
            }
            let Some((c, p)) = menor else {
                return Ok(false);
            };
            for s in 0..self.fluxos.len() {
                if Some(s) == self.af {
                    while let Some(r) = &self.fluxos[s].cab {
                        if compara_pos(&self.cromos, (r.crom, r.pos), (c, p)) != Ordering::Less {
                            break;
                        }
                        let r = self.fluxos[s].cab.take().unwrap();
                        self.fluxos[s].livres.push(r);
                        self.fluxos[s].avanca(&mut self.cromos)?;
                    }
                }
                while let Some(r) = &self.fluxos[s].cab {
                    if r.crom != c || r.pos != p {
                        break;
                    }
                    let r = self.fluxos[s].cab.take().unwrap();
                    self.grupo[s].push(r);
                    self.fluxos[s].avanca(&mut self.cromos)?;
                }
            }

            self.crom = c;
            self.tipo = self.cromos.tipo(c);
            self.pos = p;
            for m in self.membros.iter_mut() {
                m.dentro = m.bed.as_mut().map(|b| b.contem(c, p));
            }

            // As chaves desta posicao
            self.chaves.clear();
            let grupo = &self.grupo;
            for (s, g) in grupo.iter().enumerate() {
                if Some(s) == self.af {
                    continue;
                }
                for (k, r) in g.iter().enumerate() {
                    let (rf, al) = (r.refs(), r.alts());
                    if al.contains(',') {
                        self.cont.multiallelic += 1;
                        continue;
                    }
                    if !so_bases(rf) || !so_bases(al) || rf.eq_ignore_ascii_case(al) {
                        self.cont.not_plain += 1;
                        continue;
                    }
                    let repetida = self.chaves.iter().any(|ch| {
                        let o = &grupo[ch.fluxo][ch.reg];
                        o.refs().eq_ignore_ascii_case(rf) && o.alts().eq_ignore_ascii_case(al)
                    });
                    if repetida {
                        continue;
                    }
                    let snv = rf.len() == 1 && al.len() == 1 && e_acgt(rf) && e_acgt(al);
                    self.chaves.push(Chave { fluxo: s, reg: k, snv, af: None });
                }
            }
            if self.chaves.is_empty() {
                continue;
            }
            self.cont.sites += self.chaves.len() as u64;

            // O genotipo de cada pessoa em cada chave
            self.cods.clear();
            self.regs.clear();
            for ch in &self.chaves {
                let modelo = &grupo[ch.fluxo][ch.reg];
                let (rf, al) = (modelo.refs(), modelo.alts());
                for m in &self.membros {
                    let g = &grupo[m.fluxo];
                    let achado = g
                        .iter()
                        .position(|r| r.refs().eq_ignore_ascii_case(rf) && r.alts().eq_ignore_ascii_case(al));
                    let cod = match achado {
                        Some(i) => {
                            let r = &g[i];
                            if (self.pass_only && !r.passa) || m.dentro == Some(false) {
                                -1
                            } else {
                                codigo(r.gts[m.slot])
                            }
                        }
                        None if !g.is_empty() => -1,
                        None => match m.dentro {
                            Some(true) => 0,
                            Some(false) => -1,
                            None => {
                                if self.absent_is_ref {
                                    0
                                } else {
                                    -1
                                }
                            }
                        },
                    };
                    self.cods.push(cod);
                    self.regs.push(achado);
                }
            }
            if let Some(sa) = self.af {
                for ch in self.chaves.iter_mut() {
                    let modelo = &grupo[ch.fluxo][ch.reg];
                    ch.af = busca_af(&grupo[sa], modelo.refs(), modelo.alts(), &self.campo_af);
                }
            }
            return Ok(true);
        }
    }

    /// Marca o autossomo corrente como visto (para contar quantos ha).
    fn marca_autossomo(&self, vistos: &mut Vec<u32>) {
        if vistos.last() != Some(&self.crom) && !vistos.contains(&self.crom) {
            vistos.push(self.crom);
        }
    }

    /// Fecha as contagens dos fluxos.
    fn contagens(&self) -> ScanCounts {
        let mut c = self.cont.clone();
        for (s, f) in self.fluxos.iter().enumerate() {
            if Some(s) == self.af {
                c.af_records += f.registros;
            } else {
                c.records += f.registros;
                c.other_sequences += f.outras_seqs;
                c.malformed += f.malformados;
            }
        }
        c
    }
}

/// The allele frequency of (REF, ALT) in the records of the frequency file at
/// this position: the INFO value of `campo` for that ALT allele (a
/// multiallelic record carries one value per ALT).
fn busca_af(regs: &[Reg], rf: &str, al: &str, campo: &str) -> Option<f64> {
    for r in regs {
        if !r.refs().eq_ignore_ascii_case(rf) {
            continue;
        }
        let Some(i) = r.alts().split(',').position(|a| a.eq_ignore_ascii_case(al)) else {
            continue;
        };
        for kv in r.info().split(';') {
            if let Some(v) = kv.strip_prefix(campo).and_then(|x| x.strip_prefix('=')) {
                return v
                    .split(',')
                    .nth(i)
                    .and_then(|x| x.parse::<f64>().ok())
                    .filter(|x| x.is_finite() && (0.0..=1.0).contains(x));
            }
        }
    }
    None
}

// =====================================================================
//  1) PARENTESCO DE TODOS OS PARES
// =====================================================================

/// One person in a report.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PersonSummary {
    pub label: String,
    pub path: String,
    pub sample: String,
    /// Biallelic autosomal SNV sites where this person has a genotype.
    pub called: u64,
    pub hom_ref: u64,
    pub het: u64,
    pub hom_alt: u64,
    pub missing: u64,
}

impl PersonSummary {
    fn json(&self) -> String {
        format!(
            "{{\"label\":{},\"path\":{},\"sample\":{},\"called\":{},\"hom_ref\":{},\"het\":{},\"hom_alt\":{},\"missing\":{}}}",
            json_texto(&self.label),
            json_texto(&self.path),
            json_texto(&self.sample),
            self.called,
            self.hom_ref,
            self.het,
            self.hom_alt,
            self.missing
        )
    }
}

/// The comparison of two persons over the biallelic autosomal SNVs where
/// both have a genotype.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PairResult {
    pub a: String,
    pub b: String,
    /// Sites where both have a genotype.
    pub sites: u64,
    /// Both heterozygous.
    pub hethet: u64,
    /// Opposite homozygotes (0/0 against 1/1): no allele in common.
    pub ibs0: u64,
    /// `a` heterozygous, `b` homozygous.
    pub het1_hom2: u64,
    /// `b` heterozygous, `a` homozygous.
    pub het2_hom1: u64,
    /// One allele in common: het1_hom2 + het2_hom1.
    pub ibs1: u64,
    /// The same genotype.
    pub ibs2: u64,
    /// KING-robust kinship, between-family form -- the value plink2
    /// --make-king-table reports, and the one the relationship is read from:
    /// 0.5 - (4 IBS0 + HET1_HOM2 + HET2_HOM1) / (4 min(het_a, het_b)), with
    /// het_x = HETHET + the sites where only x is heterozygous. `None` when
    /// neither is heterozygous anywhere.
    pub kinship: Option<f64>,
    /// KING-robust, within-family (symmetric) form:
    /// (HETHET - 2 IBS0) / (HET1_HOM2 + HET2_HOM1 + 2 HETHET). The same as
    /// `kinship` when both have as many heterozygous sites; higher when they
    /// differ, as for two people of different populations.
    pub kinship_within: Option<f64>,
    /// IBS0 / HETHET: ~0 for parent and child, ~0.5 for unrelated persons.
    pub ibs0_ratio: Option<f64>,
    /// "duplicate_or_mz_twin", "parent_offspring", "full_siblings",
    /// "second_degree", "third_degree", "unrelated" or "insufficient_data".
    pub relationship: &'static str,
}

impl PairResult {
    fn de(a: &str, b: &str, c: &[u64; 6]) -> PairResult {
        let [n, hh, i0, h12, h21, i2] = *c;
        let denom = h12 + h21 + 2 * hh;
        let kinship_within = (denom > 0).then(|| (hh as f64 - 2.0 * i0 as f64) / denom as f64);
        let menor_het = (hh + h12).min(hh + h21);
        let kinship =
            (menor_het > 0).then(|| 0.5 - (4.0 * i0 as f64 + h12 as f64 + h21 as f64) / (4.0 * menor_het as f64));
        let ibs0_ratio = (hh > 0).then(|| i0 as f64 / hh as f64);
        PairResult {
            a: a.to_string(),
            b: b.to_string(),
            sites: n,
            hethet: hh,
            ibs0: i0,
            het1_hom2: h12,
            het2_hom1: h21,
            ibs1: h12 + h21,
            ibs2: i2,
            kinship,
            kinship_within,
            ibs0_ratio,
            relationship: classifica(kinship, ibs0_ratio, denom),
        }
    }

    pub fn json(&self) -> String {
        format!(
            "{{\"a\":{},\"b\":{},\"sites\":{},\"hethet\":{},\"ibs0\":{},\"het1_hom2\":{},\"het2_hom1\":{},\"ibs1\":{},\"ibs2\":{},\"kinship\":{},\"kinship_within\":{},\"ibs0_ratio\":{},\"relationship\":\"{}\"}}",
            json_texto(&self.a),
            json_texto(&self.b),
            self.sites,
            self.hethet,
            self.ibs0,
            self.het1_hom2,
            self.het2_hom1,
            self.ibs1,
            self.ibs2,
            json_opt(self.kinship),
            json_opt(self.kinship_within),
            json_opt(self.ibs0_ratio),
            self.relationship
        )
    }
}

/// The relationship class from the kinship (KING thresholds) and, in the
/// first-degree range, IBS0/HETHET.
pub fn classify_relationship(kinship: Option<f64>, ibs0_ratio: Option<f64>, informative_sites: u64) -> &'static str {
    classifica(kinship, ibs0_ratio, informative_sites)
}

fn classifica(kinship: Option<f64>, razao: Option<f64>, denom: u64) -> &'static str {
    let Some(k) = kinship else { return "insufficient_data" };
    if denom < MIN_KINSHIP_SITES {
        "insufficient_data"
    } else if k > KING_DUPLICATE {
        "duplicate_or_mz_twin"
    } else if k >= KING_FIRST {
        match razao {
            Some(r) if r >= PO_IBS0_RATIO => "full_siblings",
            _ if k > PO_MAX_KINSHIP => "full_siblings",
            _ => "parent_offspring",
        }
    } else if k >= KING_SECOND {
        "second_degree"
    } else if k >= KING_THIRD {
        "third_degree"
    } else {
        "unrelated"
    }
}

fn nome_relacao(r: &str) -> &'static str {
    match r {
        "duplicate_or_mz_twin" => "duplicate or identical twin",
        "parent_offspring" => "parent and child",
        "full_siblings" => "full siblings",
        "second_degree" => "second degree (half-siblings, grandparent, uncle/aunt)",
        "third_degree" => "third degree (first cousins, ...)",
        "unrelated" => "unrelated",
        _ => "insufficient data",
    }
}

const PALAVRAS: usize = 16;

/// Contagens de todos os pares, com os genotipos guardados em bits: blocos
/// de 64 sitios por palavra, tres mascaras por pessoa (0/0, 0/1, 1/1).
struct Pares {
    #[allow(dead_code)]
    p: usize,
    bits: Vec<[u64; 3]>,
    k: usize,
    lista: Vec<(u32, u32)>,
    /// sitios, hethet, ibs0, het1_hom2, het2_hom1, ibs2
    cont: Vec<[u64; 6]>,
    /// hom_ref, het, hom_alt, faltante
    pessoa: Vec<[u64; 4]>,
}

impl Pares {
    fn novo(p: usize) -> Pares {
        let mut lista = Vec::with_capacity(p * p.saturating_sub(1) / 2);
        for i in 0..p {
            for j in i + 1..p {
                lista.push((i as u32, j as u32));
            }
        }
        let n = lista.len();
        Pares { p, bits: vec![[0; 3]; p * PALAVRAS], k: 0, lista, cont: vec![[0; 6]; n], pessoa: vec![[0; 4]; p] }
    }

    fn adiciona(&mut self, cods: &[i8]) {
        let w = self.k / 64;
        let b = 1u64 << (self.k % 64);
        for (i, &c) in cods.iter().enumerate() {
            if (0..=2).contains(&c) {
                self.bits[i * PALAVRAS + w][c as usize] |= b;
                self.pessoa[i][c as usize] += 1;
            } else {
                self.pessoa[i][3] += 1;
            }
        }
        self.k += 1;
        if self.k == PALAVRAS * 64 {
            self.descarrega();
        }
    }

    fn descarrega(&mut self) {
        if self.k == 0 {
            return;
        }
        let nw = self.k.div_ceil(64);
        let bits = &self.bits;
        let conta = |&(i, j): &(u32, u32), c: &mut [u64; 6]| {
            let (bi, bj) = (i as usize * PALAVRAS, j as usize * PALAVRAS);
            for w in 0..nw {
                let [ri, hi, ai] = bits[bi + w];
                let [rj, hj, aj] = bits[bj + w];
                let hh = hi & hj;
                c[0] += u64::from(((ri | hi | ai) & (rj | hj | aj)).count_ones());
                c[1] += u64::from(hh.count_ones());
                c[2] += u64::from(((ri & aj) | (ai & rj)).count_ones());
                c[3] += u64::from((hi & (rj | aj)).count_ones());
                c[4] += u64::from((hj & (ri | ai)).count_ones());
                c[5] += u64::from(((ri & rj) | hh | (ai & aj)).count_ones());
            }
        };
        if self.lista.len() >= 256 {
            self.lista.par_iter().zip(self.cont.par_iter_mut()).for_each(|(p, c)| conta(p, c));
        } else {
            for (p, c) in self.lista.iter().zip(self.cont.iter_mut()) {
                conta(p, c);
            }
        }
        for b in self.bits.iter_mut() {
            *b = [0; 3];
        }
        self.k = 0;
    }
}

/// What `relatedness` found.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RelatednessReport {
    pub persons: Vec<PersonSummary>,
    /// Every pair, in the order of the persons: (0,1), (0,2), ..., (1,2), ...
    pub pairs: Vec<PairResult>,
    /// Biallelic autosomal SNV sites used (in the union of the files).
    pub snv_sites: u64,
    /// Autosomes with at least one site. With few, relationships are
    /// tentative (see `FEW_AUTOSOMES`).
    pub autosomes: u64,
    pub scan: ScanCounts,
}

/// Below this many autosomes the reports warn that the classes are tentative:
/// relatives share whole chromosomes by chance (two siblings can share a
/// haplotype along all of chr22, and look like parent and child there).
pub const FEW_AUTOSOMES: u64 = 10;

fn aviso_poucos(autossomos: u64) -> String {
    if autossomos >= FEW_AUTOSOMES {
        return String::new();
    }
    format!(
        "  WARNING: the data cover only {autossomos} autosome(s). On so little of the genome relatives share whole chromosomes by chance: full siblings can share a haplotype along all of a chromosome and look like parent and child (or, sharing both, like twins), and unrelated pairs scatter more. Read the relationships as tentative and use genome-wide data to conclude.\n"
    )
}

impl RelatednessReport {
    fn aviso_cromossomos(&self) -> String {
        aviso_poucos(self.autosomes)
    }

    pub fn json(&self) -> String {
        let p: Vec<String> = self.persons.iter().map(PersonSummary::json).collect();
        let q: Vec<String> = self.pairs.iter().map(PairResult::json).collect();
        format!(
            "{{\"persons\":[{}],\"pairs\":[{}],\"snv_sites\":{},\"autosomes\":{},\"scan\":{},\"thresholds\":{{\"duplicate\":{},\"first_degree\":{},\"second_degree\":{},\"third_degree\":{},\"parent_offspring_ibs0_ratio\":{},\"parent_offspring_max_kinship\":{}}}}}",
            p.join(","),
            q.join(","),
            self.snv_sites,
            self.autosomes,
            self.scan.json(),
            KING_DUPLICATE,
            KING_FIRST,
            KING_SECOND,
            KING_THIRD,
            PO_IBS0_RATIO,
            PO_MAX_KINSHIP
        )
    }

    /// The report as text, ready to print.
    pub fn resumo(&self) -> String {
        let mut o = format!(
            "Relatedness of {} person(s) over {} biallelic autosomal SNV site(s), KING-robust kinship\n",
            self.persons.len(),
            milhar(self.snv_sites)
        );
        for p in &self.persons {
            o.push_str(&format!(
                "  {}: {} called ({} het, {} hom-alt), {} missing\n",
                p.label,
                milhar(p.called),
                milhar(p.het),
                milhar(p.hom_alt),
                milhar(p.missing)
            ));
        }
        let largura = self.pairs.iter().map(|p| p.a.len() + p.b.len() + 3).max().unwrap_or(4).max(4);
        o.push_str(&format!(
            "  {:<largura$}  {:>12}  {:>9}  {:>11}  {:>8}  {:>8}  relationship\n",
            "pair", "sites", "IBS0", "IBS0/HETHET", "kinship", "within"
        ));
        for p in &self.pairs {
            let prop = if p.sites > 0 { format!("{:.2e}", p.ibs0 as f64 / p.sites as f64) } else { "n/a".to_string() };
            o.push_str(&format!(
                "  {:<largura$}  {:>12}  {:>9}  {:>11}  {:>8}  {:>8}  {}\n",
                format!("{} - {}", p.a, p.b),
                milhar(p.sites),
                prop,
                fmt_opt(p.ibs0_ratio, 4),
                fmt_opt(p.kinship, 4),
                fmt_opt(p.kinship_within, 4),
                nome_relacao(p.relationship)
            ));
        }
        o.push_str(&self.scan.nota());
        o.push_str(&self.aviso_cromossomos());
        o.push_str("  IBS0 is the share of sites where the two are opposite homozygotes. kinship: KING-robust, the between-family form plink2 --make-king-table reports, 0.5 - (4 IBS0 + HET1_HOM2 + HET2_HOM1) / (4 x the smaller heterozygous count); within: the symmetric form (HETHET - 2 IBS0) / (HET1_HOM2 + HET2_HOM1 + 2 HETHET), higher for two people of different populations.\n");
        o.push_str(&format!(
            "  Classes from the kinship (KING thresholds: duplicate > {KING_DUPLICATE}, first degree >= {KING_FIRST}, second >= {KING_SECOND}, third >= {KING_THIRD}); a first-degree pair is parent and child when IBS0/HETHET < {PO_IBS0_RATIO} and the kinship is at most {PO_MAX_KINSHIP} (above it the pair shares both haplotypes somewhere, as siblings do), otherwise full siblings. These are statistical estimates for teaching and research, not a medical or legal test.\n"
        ));
        o
    }
}

/// Pairwise relatedness of every pair of persons, in one streaming pass over
/// their files (a multi-sample VCF can supply several persons).
///
/// Uses the biallelic SNVs of the autosomes. The time grows with the number
/// of pairs: a few hundred persons are fine; thousands are slow.
pub fn relatedness(persons: &[Person], opts: &JoinOptions) -> Result<RelatednessReport, String> {
    if persons.len() < 2 {
        return Err("relatedness needs at least two persons.".to_string());
    }
    let mut v = Varredura::nova(persons, None, "", opts)?;
    let np = v.membros.len();
    let mut acc = Pares::novo(np);
    let mut sitios = 0u64;
    let mut autossomos: Vec<u32> = Vec::new();
    while v.proxima()? {
        if v.tipo != TipoSeq::Autossomo {
            continue;
        }
        for (k, ch) in v.chaves.iter().enumerate() {
            if ch.snv {
                acc.adiciona(&v.cods[k * np..(k + 1) * np]);
                sitios += 1;
                v.marca_autossomo(&mut autossomos);
            }
        }
    }
    acc.descarrega();
    let persons_out = v
        .membros
        .iter()
        .zip(acc.pessoa.iter())
        .map(|(m, c)| PersonSummary {
            label: m.rotulo.clone(),
            path: m.caminho.clone(),
            sample: m.amostra.clone(),
            called: c[0] + c[1] + c[2],
            hom_ref: c[0],
            het: c[1],
            hom_alt: c[2],
            missing: c[3],
        })
        .collect();
    let pairs = acc
        .lista
        .iter()
        .zip(acc.cont.iter())
        .map(|(&(i, j), c)| PairResult::de(&v.membros[i as usize].rotulo, &v.membros[j as usize].rotulo, c))
        .collect();
    Ok(RelatednessReport {
        persons: persons_out,
        pairs,
        snv_sites: sitios,
        autosomes: autossomos.len() as u64,
        scan: v.contagens(),
    })
}

// =====================================================================
//  2) TRIOS: MENDEL, EXCLUSOES, DE NOVO, CPI
// =====================================================================

/// Options of `trio_analysis`.
#[derive(Clone, Debug, PartialEq)]
pub struct TrioOptions {
    pub join: JoinOptions,
    /// A sorted VCF of population allele frequencies (a sites file is
    /// enough). Needed for the paternity index and for the exclusion rate an
    /// unrelated man would give.
    pub af_path: Option<String>,
    /// The INFO field of `af_path` with the ALT frequency: "AF", "EUR_AF"...
    pub af_field: String,
    /// Genotyping-error (and mutation) term of the paternity index: a parent
    /// transmits the allele it does not carry with this probability, so one
    /// wrong genotype lowers the index instead of zeroing it. Also the error
    /// level assumed for the verdict when there is no mother. Default 0.001.
    pub genotype_error: f64,
    /// Prior probability of paternity. Default 0.5.
    pub prior: f64,
    /// Where to write the de novo candidates as a VCF of the child's records.
    pub de_novo_output: Option<String>,
}

impl Default for TrioOptions {
    fn default() -> Self {
        TrioOptions {
            join: JoinOptions::default(),
            af_path: None,
            af_field: "AF".to_string(),
            genotype_error: 0.001,
            prior: 0.5,
            de_novo_output: None,
        }
    }
}

/// Mendelian error types for a trio, with REF = A and ALT = B. They are the
/// codes 1-8 of plink's --mendel, in that order.
pub const TRIO_ERROR_TYPES: [&str; 8] = [
    "parents_hom_ref_child_het",
    "parents_hom_alt_child_het",
    "father_hom_alt_child_hom_ref",
    "mother_hom_alt_child_hom_ref",
    "parents_hom_alt_child_hom_ref",
    "father_hom_ref_child_hom_alt",
    "mother_hom_ref_child_hom_alt",
    "parents_hom_ref_child_hom_alt",
];

/// Mendelian error types for a duo (father and child only).
pub const DUO_ERROR_TYPES: [&str; 2] = ["father_hom_alt_child_hom_ref", "father_hom_ref_child_hom_alt"];

/// The plink --mendel code (1-8) of a trio genotype combination, or 0 when
/// it is consistent. Genotypes are ALT counts 0, 1, 2.
pub fn mendel_code(father: i8, mother: i8, child: i8) -> u8 {
    match child {
        1 if father == 0 && mother == 0 => 1,
        1 if father == 2 && mother == 2 => 2,
        0 if father == 2 && mother == 2 => 5,
        0 if father == 2 => 3,
        0 if mother == 2 => 4,
        2 if father == 0 && mother == 0 => 8,
        2 if father == 0 => 6,
        2 if mother == 0 => 7,
        _ => 0,
    }
}

/// What `trio_analysis` found. Sites are biallelic and autosomal; the
/// Mendelian counts use every such site (SNV, indel, MNP), the kinship,
/// exclusions and paternity index only SNVs.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TrioReport {
    pub child: String,
    pub father: String,
    pub mother: Option<String>,
    /// Sites where the child, the father (and the mother) all have a genotype.
    pub sites_checked: u64,
    pub consistent: u64,
    pub errors: u64,
    pub snv_checked: u64,
    pub snv_errors: u64,
    /// Indels and MNPs.
    pub indel_checked: u64,
    pub indel_errors: u64,
    /// (type, count): the 8 types of `TRIO_ERROR_TYPES`, or the 2 of
    /// `DUO_ERROR_TYPES` without a mother.
    pub errors_by_type: Vec<(&'static str, u64)>,
    /// Trio: 27 counts, index (father*3 + mother)*3 + child; duo: 9 counts,
    /// index father*3 + child. Genotype = number of ALT alleles.
    pub genotype_table: Vec<u64>,
    pub father_child: PairResult,
    pub mother_child: Option<PairResult>,
    /// SNVs where the child is homozygous and the father has a genotype.
    pub father_informative: u64,
    /// Of those, the father is homozygous for the other allele: he cannot
    /// have transmitted the child's allele (an exclusion, or an error).
    pub father_exclusions: u64,
    pub mother_informative: u64,
    pub mother_exclusions: u64,
    /// Child heterozygous, both parents homozygous reference (with a genotype).
    pub de_novo: u64,
    pub de_novo_snv: u64,
    pub de_novo_indel: u64,
    /// Records written to `de_novo_output` (0 when there is none).
    pub de_novo_written: u64,
    /// father_exclusions / father_informative.
    pub exclusion_rate: f64,
    /// The exclusion rate genotyping errors alone would give.
    pub error_level: f64,
    /// "mother_child" (measured on the mother, when she looks like the
    /// child's parent) or "genotype_error" (the option).
    pub error_level_source: &'static str,
    /// The exclusion rate expected from an unrelated man, from the allele
    /// frequencies (`None` without them).
    pub unrelated_expected_rate: Option<f64>,
    /// "not_excluded", "excluded" or "inconclusive".
    pub verdict: &'static str,
    pub reason: String,
    /// SNVs used in the paternity index (all genotyped and with a frequency).
    pub cpi_sites: u64,
    /// SNVs with every genotype but no frequency in the frequency file (0
    /// when no frequency file was given).
    pub cpi_sites_without_af: u64,
    /// log10 of the combined paternity index; `None` without frequencies.
    pub log10_cpi: Option<f64>,
    /// Probability of paternity from the CPI and the prior.
    pub probability_of_paternity: Option<f64>,
    pub prior: f64,
    pub genotype_error: f64,
    /// Autosomes with at least one SNV used.
    pub autosomes: u64,
    pub scan: ScanCounts,
}

impl TrioReport {
    pub fn json(&self) -> String {
        let tipos: Vec<String> = self.errors_by_type.iter().map(|(n, c)| format!("\"{n}\":{c}")).collect();
        format!(
            "{{\"child\":{},\"father\":{},\"mother\":{},\"sites_checked\":{},\"consistent\":{},\"errors\":{},\"error_rate\":{},\"snv_checked\":{},\"snv_errors\":{},\"indel_checked\":{},\"indel_errors\":{},\"errors_by_type\":{{{}}},\"genotype_table\":{},\"father_child\":{},\"mother_child\":{},\"father_informative\":{},\"father_exclusions\":{},\"mother_informative\":{},\"mother_exclusions\":{},\"de_novo\":{},\"de_novo_snv\":{},\"de_novo_indel\":{},\"de_novo_written\":{},\"exclusion_rate\":{},\"error_level\":{},\"error_level_source\":\"{}\",\"unrelated_expected_rate\":{},\"verdict\":\"{}\",\"reason\":{},\"cpi_sites\":{},\"cpi_sites_without_af\":{},\"log10_cpi\":{},\"probability_of_paternity\":{},\"prior\":{},\"genotype_error\":{},\"autosomes\":{},\"scan\":{}}}",
            json_texto(&self.child),
            json_texto(&self.father),
            self.mother.as_deref().map(json_texto).unwrap_or_else(|| "null".to_string()),
            self.sites_checked,
            self.consistent,
            self.errors,
            json_num(self.errors as f64 / self.sites_checked.max(1) as f64),
            self.snv_checked,
            self.snv_errors,
            self.indel_checked,
            self.indel_errors,
            tipos.join(","),
            json_lista_u64(&self.genotype_table),
            self.father_child.json(),
            self.mother_child.as_ref().map(PairResult::json).unwrap_or_else(|| "null".to_string()),
            self.father_informative,
            self.father_exclusions,
            self.mother_informative,
            self.mother_exclusions,
            self.de_novo,
            self.de_novo_snv,
            self.de_novo_indel,
            self.de_novo_written,
            json_num(self.exclusion_rate),
            json_num(self.error_level),
            self.error_level_source,
            json_opt(self.unrelated_expected_rate),
            self.verdict,
            json_texto(&self.reason),
            self.cpi_sites,
            self.cpi_sites_without_af,
            json_opt(self.log10_cpi),
            json_opt(self.probability_of_paternity),
            json_num(self.prior),
            json_num(self.genotype_error),
            self.autosomes,
            self.scan.json()
        )
    }

    pub fn resumo(&self) -> String {
        let quem = match &self.mother {
            Some(m) => format!("child {}, father {}, mother {}", self.child, self.father, m),
            None => format!("child {}, father {} (no mother: duo)", self.child, self.father),
        };
        let mut o = format!("Trio analysis: {quem}\n");
        o.push_str(&format!(
            "  Mendelian check: {} site(s) with every genotype, {} consistent, {} error(s) ({:.4}%): SNV {} of {}, indel {} of {}\n",
            milhar(self.sites_checked),
            milhar(self.consistent),
            milhar(self.errors),
            100.0 * self.errors as f64 / self.sites_checked.max(1) as f64,
            milhar(self.snv_errors),
            milhar(self.snv_checked),
            milhar(self.indel_errors),
            milhar(self.indel_checked)
        ));
        let nz: Vec<String> =
            self.errors_by_type.iter().filter(|(_, c)| *c > 0).map(|(n, c)| format!("{n} {}", milhar(*c))).collect();
        if !nz.is_empty() {
            o.push_str(&format!("  errors by type: {}\n", nz.join(", ")));
        }
        o.push_str(&format!(
            "  father-child: kinship {}, IBS0/HETHET {} ({}); opposite homozygotes {} of {} informative SNVs\n",
            fmt_opt(self.father_child.kinship, 4),
            fmt_opt(self.father_child.ibs0_ratio, 4),
            nome_relacao(self.father_child.relationship),
            milhar(self.father_exclusions),
            milhar(self.father_informative)
        ));
        if let Some(mc) = &self.mother_child {
            o.push_str(&format!(
                "  mother-child: kinship {}, IBS0/HETHET {} ({}); opposite homozygotes {} of {} informative SNVs\n",
                fmt_opt(mc.kinship, 4),
                fmt_opt(mc.ibs0_ratio, 4),
                nome_relacao(mc.relationship),
                milhar(self.mother_exclusions),
                milhar(self.mother_informative)
            ));
            o.push_str(&format!(
                "  de novo candidates (child het, both parents hom-ref): {} ({} SNV, {} indel){}\n",
                milhar(self.de_novo),
                milhar(self.de_novo_snv),
                milhar(self.de_novo_indel),
                if self.de_novo_written > 0 { format!(", {} written", milhar(self.de_novo_written)) } else { String::new() }
            ));
        }
        o.push_str(&format!("  verdict: {} -- {}\n", self.verdict, self.reason));
        match self.log10_cpi {
            Some(l) => o.push_str(&format!(
                "  combined paternity index: log10 CPI = {l:.2} over {} SNV(s); probability of paternity {} (prior {}); {} SNV(s) had no frequency\n",
                milhar(self.cpi_sites),
                self.probability_of_paternity.map_or("n/a".to_string(), prob_texto),
                self.prior,
                milhar(self.cpi_sites_without_af)
            )),
            None if self.cpi_sites_without_af > 0 && self.cpi_sites == 0 => o.push_str(
                "  combined paternity index: no SNV had a population frequency (check the af_field and that the frequency file uses the same REF/ALT)\n",
            ),
            None => o.push_str("  combined paternity index: not computed (give a population allele-frequency VCF)\n"),
        }
        o.push_str(&self.scan.nota());
        o.push_str(&aviso_poucos(self.autosomes));
        o.push_str("  The index treats SNVs as independent (linked SNVs make it much larger than a forensic panel would); for teaching and research, not a legal or medical paternity test.\n");
        o
    }
}

fn prob_texto(w: f64) -> String {
    if w >= 1.0 - 1e-12 {
        "> 0.999999999999".to_string()
    } else if w > 0.0 && w < 1e-12 {
        "< 1e-12".to_string()
    } else {
        format!("{w:.12}")
    }
}

/// Probabilidade de o pai transmitir o alelo `alelo` (0 = REF, 1 = ALT)
/// tendo `g` copias do ALT, com o termo de erro: com prob. `eps` o alelo
/// transmitido e trocado.
fn transmite(g: i8, alelo: i8, eps: f64) -> f64 {
    let t = f64::from(g) / 2.0;
    let t1 = (1.0 - eps) * t + eps * (1.0 - t);
    if alelo == 1 {
        t1
    } else {
        1.0 - t1
    }
}

/// The paternity index of one biallelic site: P(child | mother, alleged
/// father) / P(child | mother, random man), with Hardy-Weinberg population
/// frequency `p` of ALT and the transmission error `eps`. Without the
/// mother, her allele comes from the population. Genotypes are ALT counts.
pub fn paternity_index(child: i8, father: i8, mother: Option<i8>, p: f64, eps: f64) -> f64 {
    let pop = |a: i8| if a == 1 { p } else { 1.0 - p };
    let (mut x, mut y) = (0.0, 0.0);
    for a in 0..=1i8 {
        let b = child - a;
        if !(0..=1).contains(&b) {
            continue;
        }
        let mae = match mother {
            Some(m) => transmite(m, b, eps),
            None => pop(b),
        };
        x += transmite(father, a, eps) * mae;
        y += pop(a) * mae;
    }
    x / y
}

/// Probability of paternity from log10 CPI and the prior.
pub fn probability_of_paternity(log10_cpi: f64, prior: f64) -> f64 {
    if log10_cpi == f64::NEG_INFINITY {
        return 0.0;
    }
    let z = ((1.0 - prior) / prior).log10() - log10_cpi;
    if z > 300.0 {
        0.0
    } else if z < -300.0 {
        1.0
    } else {
        1.0 / (1.0 + 10f64.powf(z))
    }
}

#[derive(Default, Clone, Copy)]
struct ContaPar([u64; 6]);

impl ContaPar {
    fn adiciona(&mut self, a: i8, b: i8) {
        if a < 0 || b < 0 {
            return;
        }
        let c = &mut self.0;
        c[0] += 1;
        match (a, b) {
            (1, 1) => {
                c[1] += 1;
                c[5] += 1;
            }
            (0, 2) | (2, 0) => c[2] += 1,
            (1, _) => c[3] += 1,
            (_, 1) => c[4] += 1,
            _ => c[5] += 1,
        }
    }
}

/// The verdict from the exclusion counts. `nivel` is the exclusion rate
/// errors alone would give; `u` the one of an unrelated man.
fn veredito(inf: u64, exc: u64, nivel: f64, u: Option<f64>, parece_pai: bool) -> (&'static str, String) {
    if inf < MIN_INFORMATIVE_SITES {
        return (
            "inconclusive",
            format!(
                "only {inf} informative SNV(s) (child homozygous and father genotyped); at least {MIN_INFORMATIVE_SITES} are needed"
            ),
        );
    }
    let taxa = exc as f64 / inf as f64;
    let lambda = nivel * inf as f64;
    let tolera = lambda + 4.0 * lambda.sqrt() + 3.0;
    let base = format!(
        "{} exclusion(s) in {} informative SNVs (rate {:.2e}); genotyping errors alone would give about {:.1} (error level {:.2e})",
        milhar(exc),
        milhar(inf),
        taxa,
        lambda,
        nivel
    );
    let comp_u = match u {
        Some(u) if u > 0.0 => format!(", and an unrelated man about {:.0} (rate {:.2e})", u * inf as f64, u),
        _ => String::new(),
    };
    if exc as f64 <= tolera {
        return ("not_excluded", format!("{base}{comp_u}. The exclusions are within the error level: consistent with paternity, which this does not prove."));
    }
    let longe_do_erro = taxa >= 10.0 * nivel;
    let perto_de_estranho = u.is_none_or(|u| taxa >= 0.1 * u);
    if longe_do_erro && perto_de_estranho {
        if parece_pai {
            return ("inconclusive", format!("{base}{comp_u}. The exclusions are far above the error level, yet the kinship and IBS0/HETHET look like parent and child: check the genotype quality of the two samples."));
        }
        return ("excluded", format!("{base}{comp_u}. The exclusions are far above what errors explain: the man is not consistent with being the biological father."));
    }
    if parece_pai {
        // Pouco acima do nivel de erro, e o KING diz pai-filho: IBS0/HETHET
        // abaixo de PO_IBS0_RATIO ja e uma medida de exclusoes, numa escala
        // que nao depende do painel, e nenhum parente nao-pai fica abaixo dela.
        return ("not_excluded", format!("{base}{comp_u}. The exclusions are somewhat above the assumed error level, but the kinship and IBS0/HETHET are those of parent and child, which no other relative gives: consistent with paternity, which this does not prove."));
    }
    ("inconclusive", format!("{base}{comp_u}. The exclusions are above the error level but not clearly at the level of a non-father: check the data (genotype quality, sample identity)."))
}

/// Trio (or duo) analysis: Mendelian consistency, exclusions, de novo
/// candidates and, with population frequencies, the combined paternity index.
pub fn trio_analysis(child: &Person, father: &Person, mother: Option<&Person>, opts: &TrioOptions) -> Result<TrioReport, String> {
    if !(0.0..0.5).contains(&opts.genotype_error) {
        return Err(format!(
            "genotype_error must be at least 0 and below 0.5 (it is a probability per site); got {}.",
            opts.genotype_error
        ));
    }
    if !(opts.prior > 0.0 && opts.prior < 1.0) {
        return Err(format!("prior must be between 0 and 1, both excluded; got {}.", opts.prior));
    }
    let mut pessoas = vec![child.clone(), father.clone()];
    if let Some(m) = mother {
        pessoas.push(m.clone());
    }
    let trio = mother.is_some();
    if let Some(saida) = &opts.de_novo_output {
        if !trio {
            return Err("De novo candidates need the mother too (child het where BOTH parents are hom-ref): give mother=... or leave de_novo_output out.".to_string());
        }
        let mut ent: Vec<&str> = pessoas.iter().map(|p| p.path.as_str()).collect();
        if let Some(a) = &opts.af_path {
            ent.push(a);
        }
        recusa_sobrescrever_entrada(saida, &ent)?;
    }
    let mut v = Varredura::nova(&pessoas, opts.af_path.as_deref(), &opts.af_field, &opts.join)?;
    let np = v.membros.len();

    // A saida de novo: cabecalho do arquivo do filho, so com a coluna dele.
    let mut guarda = Saida { caminho: None };
    let mut escritor = match &opts.de_novo_output {
        Some(caminho) => {
            let mut w = Escritor::cria_para(caminho, Formato::Vcf)?;
            guarda.caminho = Some(caminho.clone());
            let f = &v.fluxos[v.membros[0].fluxo];
            for l in &f.cabecalho {
                w.linha(l)?;
            }
            w.linha("##MMGEN_trio_de_novo=<Description=\"Child heterozygous where both parents are homozygous reference (with a genotype); candidates, not confirmed mutations\">")?;
            w.linha(&format!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", v.membros[0].amostra))?;
            Some(w)
        }
        None => None,
    };

    let eps = opts.genotype_error;
    let mut tabela = vec![0u64; if trio { 27 } else { 9 }];
    let mut tab_snv = vec![0u64; tabela.len()];
    let (mut kf, mut km) = (ContaPar::default(), ContaPar::default());
    let (mut inf_f, mut exc_f, mut inf_m, mut exc_m) = (0u64, 0u64, 0u64, 0u64);
    let (mut u_soma, mut u_n) = (0f64, 0u64);
    let (mut dn, mut dn_snv, mut escritos) = (0u64, 0u64, 0u64);
    let (mut cpi, mut cpi_n, mut sem_af) = (0f64, 0u64, 0u64);
    let mut autossomos: Vec<u32> = Vec::new();

    while v.proxima()? {
        if v.tipo != TipoSeq::Autossomo {
            continue;
        }
        for (k, ch) in v.chaves.iter().enumerate() {
            let c = v.cods[k * np];
            let f = v.cods[k * np + 1];
            let m = if trio { v.cods[k * np + 2] } else { -1 };
            if trio {
                if c >= 0 && f >= 0 && m >= 0 {
                    let i = ((f * 3 + m) * 3 + c) as usize;
                    tabela[i] += 1;
                    if ch.snv {
                        tab_snv[i] += 1;
                    }
                }
                if c == 1 && f == 0 && m == 0 {
                    dn += 1;
                    if ch.snv {
                        dn_snv += 1;
                    }
                    if let (Some(w), Some(i)) = (escritor.as_mut(), v.regs[k * np]) {
                        let r = &v.grupo[v.membros[0].fluxo][i];
                        let (a, b) = r.cols[v.membros[0].slot];
                        let mut l = String::with_capacity(r.fim_fixas + b - a + 1);
                        l.push_str(&r.linha[..r.fim_fixas]);
                        l.push('\t');
                        l.push_str(&r.linha[a..b]);
                        w.linha(&l)?;
                        escritos += 1;
                    }
                }
            } else if c >= 0 && f >= 0 {
                let i = (f * 3 + c) as usize;
                tabela[i] += 1;
                if ch.snv {
                    tab_snv[i] += 1;
                }
            }
            if !ch.snv {
                continue;
            }
            v.marca_autossomo(&mut autossomos);
            kf.adiciona(f, c);
            km.adiciona(m, c);
            if c == 0 || c == 2 {
                if f >= 0 {
                    inf_f += 1;
                    if f == 2 - c {
                        exc_f += 1;
                    }
                    if let Some(p) = ch.af {
                        u_soma += if c == 0 { p * p } else { (1.0 - p) * (1.0 - p) };
                        u_n += 1;
                    }
                }
                if m >= 0 {
                    inf_m += 1;
                    if m == 2 - c {
                        exc_m += 1;
                    }
                }
            }
            if c >= 0 && f >= 0 && (!trio || m >= 0) {
                match ch.af {
                    Some(p) => {
                        let p = p.clamp(MIN_AF, 1.0 - MIN_AF);
                        let pi = paternity_index(c, f, trio.then_some(m), p, eps);
                        if pi.is_nan() {
                            // eps = 0 e mae incompativel com o filho: o sitio nao informa.
                            continue;
                        }
                        cpi += pi.log10();
                        cpi_n += 1;
                    }
                    None if opts.af_path.is_some() => sem_af += 1,
                    None => {}
                }
            }
        }
    }
    if let Some(w) = escritor.take() {
        w.termina()?;
    }
    guarda.desarma();

    let rotulo = |i: usize| v.membros[i].rotulo.clone();
    let father_child = PairResult::de(&rotulo(1), &rotulo(0), &kf.0);
    let mother_child = trio.then(|| PairResult::de(&rotulo(2), &rotulo(0), &km.0));

    let mut rep = TrioReport {
        child: rotulo(0),
        father: rotulo(1),
        mother: trio.then(|| rotulo(2)),
        genotype_table: tabela.clone(),
        father_informative: inf_f,
        father_exclusions: exc_f,
        mother_informative: if trio { inf_m } else { 0 },
        mother_exclusions: if trio { exc_m } else { 0 },
        de_novo: dn,
        de_novo_snv: dn_snv,
        de_novo_indel: dn - dn_snv,
        de_novo_written: escritos,
        prior: opts.prior,
        genotype_error: eps,
        cpi_sites: cpi_n,
        cpi_sites_without_af: sem_af,
        autosomes: autossomos.len() as u64,
        ..Default::default()
    };
    // Mendel
    let erro = |i: usize| -> u8 {
        if trio {
            let (f, m, c) = ((i / 9) as i8, ((i / 3) % 3) as i8, (i % 3) as i8);
            mendel_code(f, m, c)
        } else {
            let (f, c) = ((i / 3) as i8, (i % 3) as i8);
            if f == 2 - c && c != 1 {
                if c == 0 {
                    1
                } else {
                    2
                }
            } else {
                0
            }
        }
    };
    let mut por_tipo = vec![0u64; if trio { 8 } else { 2 }];
    for (i, (&n, &ns)) in tabela.iter().zip(tab_snv.iter()).enumerate() {
        rep.sites_checked += n;
        rep.snv_checked += ns;
        let e = erro(i);
        if e > 0 {
            rep.errors += n;
            rep.snv_errors += ns;
            por_tipo[e as usize - 1] += n;
        }
    }
    rep.consistent = rep.sites_checked - rep.errors;
    rep.indel_checked = rep.sites_checked - rep.snv_checked;
    rep.indel_errors = rep.errors - rep.snv_errors;
    let nomes: &[&'static str] = if trio { &TRIO_ERROR_TYPES } else { &DUO_ERROR_TYPES };
    rep.errors_by_type = nomes.iter().copied().zip(por_tipo).collect();

    // Veredito
    rep.exclusion_rate = if inf_f > 0 { exc_f as f64 / inf_f as f64 } else { 0.0 };
    let mae_parece_mae = mother_child.as_ref().is_some_and(|p| p.relationship == "parent_offspring");
    let mut nota_mae = String::new();
    if trio && mae_parece_mae && inf_m >= 1000 {
        rep.error_level = (exc_m as f64 + 1.0) / inf_m as f64;
        rep.error_level_source = "mother_child";
    } else {
        rep.error_level = eps;
        rep.error_level_source = "genotype_error";
        if trio && !mae_parece_mae {
            nota_mae = " The mother-child pair does not look like parent and child, so the error level is the genotype_error option, not a measurement.".to_string();
        }
    }
    rep.unrelated_expected_rate = (u_n >= MIN_INFORMATIVE_SITES).then(|| u_soma / u_n as f64);
    let parece_pai = matches!(father_child.relationship, "parent_offspring" | "duplicate_or_mz_twin");
    let (vd, razao) = veredito(inf_f, exc_f, rep.error_level, rep.unrelated_expected_rate, parece_pai);
    rep.verdict = vd;
    rep.reason = format!("{razao}{nota_mae}");
    if opts.af_path.is_some() && cpi_n > 0 {
        rep.log10_cpi = Some(cpi);
        rep.probability_of_paternity = Some(probability_of_paternity(cpi, opts.prior));
    }
    rep.father_child = father_child;
    rep.mother_child = mother_child;
    rep.scan = v.contagens();
    Ok(rep)
}

// =====================================================================
//  3) CORRIDAS DE HOMOZIGOSE
// =====================================================================

/// Options of `runs_of_homozygosity`.
#[derive(Clone, Debug, PartialEq)]
pub struct RohOptions {
    /// At least this many SNV calls in a run. Default 50.
    pub min_sites: u64,
    /// At least this length, in bases (first to last site). Default 1,000,000.
    pub min_length: u64,
    /// A stretch without data longer than this breaks a run: the distance
    /// between two calls, or, with a callable BED, the longest stretch
    /// outside it. Default 1,000,000.
    pub max_gap: u64,
    /// Length of the autosomes for F_ROH. `None`: GRCh38 (2,875,001,522).
    pub autosome_length: Option<u64>,
    pub pass_only: bool,
    /// The sliding-window method of `plink --homozyg`: windows of this many
    /// consecutive SNV calls slide along each stretch with data; a window is
    /// homozygous when it has at most `window_het` heterozygous and
    /// `window_missing` missing calls; a call belongs to a run when at least
    /// `window_threshold` of the windows that cover it are homozygous. A
    /// genotyping error (one het call in a long run) no longer splits the run.
    /// 0 = the strict method: any heterozygous call ends a run. Default 50.
    pub window_snps: u64,
    /// Heterozygous calls a homozygous window may hold. Default 1.
    pub window_het: u64,
    /// Missing calls a homozygous window may hold. Default 5.
    pub window_missing: u64,
    /// Share of the windows covering a call that must be homozygous. Default 0.05.
    pub window_threshold: f64,
    /// Window method: a run with fewer than one call per this many bases on
    /// average is dropped (`plink --homozyg-density`, in bases). 0 = no limit.
    /// Default 50,000.
    pub max_density: u64,
}

impl Default for RohOptions {
    fn default() -> Self {
        RohOptions {
            min_sites: 50,
            min_length: 1_000_000,
            max_gap: 1_000_000,
            autosome_length: None,
            pass_only: false,
            window_snps: 50,
            window_het: 1,
            window_missing: 5,
            window_threshold: 0.05,
            max_density: 50_000,
        }
    }
}

/// The runs of one stretch of calls (one chromosome, no gap over `max_gap`)
/// by the sliding-window method. `chamadas`: (position, 0 hom, 1 het, 2
/// missing), in order.
fn roh_janelas(chamadas: &[(u64, u8)], nome: &str, o: &RohOptions, runs: &mut Vec<RohRun>) {
    let n = chamadas.len();
    let w = o.window_snps as usize;
    if w == 0 || n < w {
        return;
    }
    let mut het = vec![0u64; n + 1];
    let mut falta = vec![0u64; n + 1];
    for (i, &(_, s)) in chamadas.iter().enumerate() {
        het[i + 1] = het[i] + u64::from(s == 1);
        falta[i + 1] = falta[i] + u64::from(s == 2);
    }
    let nj = n - w + 1;
    // acumulado das janelas homozigotas (a janela i cobre as chamadas i..i+w)
    let mut hom = vec![0u64; nj + 1];
    for i in 0..nj {
        let e_hom = het[i + w] - het[i] <= o.window_het && falta[i + w] - falta[i] <= o.window_missing;
        hom[i + 1] = hom[i] + u64::from(e_hom);
    }
    let candidata = |j: usize| {
        let (a, b) = (j.saturating_sub(w - 1), j.min(nj - 1));
        let total = (b + 1 - a) as f64;
        let h = hom[b + 1] - hom[a];
        h > 0 && h as f64 >= o.window_threshold * total
    };
    let mut j = 0;
    while j < n {
        if !candidata(j) {
            j += 1;
            continue;
        }
        let ini = j;
        while j < n && candidata(j) {
            j += 1;
        }
        let (p0, p1, sites) = (chamadas[ini].0, chamadas[j - 1].0, (j - ini) as u64);
        let comprimento = p1 - p0 + 1;
        let densa = o.max_density == 0 || comprimento <= o.max_density.saturating_mul(sites);
        if sites >= o.min_sites && comprimento >= o.min_length && densa {
            runs.push(RohRun { chrom: nome.to_string(), start: p0 - 1, end: p1, sites });
        }
    }
}

/// One run of homozygosity. `start` is 0-based (BED), `end` the 1-based
/// position of its last site.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RohRun {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub sites: u64,
}

impl RohRun {
    pub fn length(&self) -> u64 {
        self.end - self.start
    }
}

/// What `runs_of_homozygosity` found.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RohReport {
    pub person: String,
    pub runs: Vec<RohRun>,
    pub count: u64,
    pub total_length: u64,
    pub longest: u64,
    pub mean_length: f64,
    pub autosome_length: u64,
    /// total_length / autosome_length.
    pub f_roh: f64,
    /// Sum, over the autosomes in the file, of first-to-last used call.
    pub covered_span: u64,
    /// total_length / covered_span: F_ROH over the part of the genome the
    /// file covers (useful when it has one chromosome).
    pub f_roh_covered: f64,
    /// Autosomal SNV calls used (inside the BED, with a genotype).
    pub snv_calls: u64,
    pub het_calls: u64,
    pub hom_calls: u64,
    pub missing: u64,
    pub outside_bed: u64,
    /// Autosomal records that are not SNVs (indels, MNPs, symbolic): ignored.
    pub non_snv: u64,
    pub min_sites: u64,
    pub min_length: u64,
    pub max_gap: u64,
    /// "window" (as plink --homozyg) or "strict" (any het call ends a run).
    pub method: &'static str,
    pub window_snps: u64,
    pub window_het: u64,
    pub window_missing: u64,
    pub window_threshold: f64,
    pub max_density: u64,
    /// Lines written to the BED output.
    pub written: u64,
    pub scan: ScanCounts,
}

impl RohReport {
    pub fn json(&self) -> String {
        let runs: Vec<String> = self
            .runs
            .iter()
            .map(|r| format!("[{},{},{},{}]", json_texto(&r.chrom), r.start, r.end, r.sites))
            .collect();
        format!(
            "{{\"person\":{},\"runs\":[{}],\"count\":{},\"total_length\":{},\"longest\":{},\"mean_length\":{},\"autosome_length\":{},\"f_roh\":{},\"covered_span\":{},\"f_roh_covered\":{},\"snv_calls\":{},\"het_calls\":{},\"hom_calls\":{},\"missing\":{},\"outside_bed\":{},\"non_snv\":{},\"min_sites\":{},\"min_length\":{},\"max_gap\":{},\"method\":\"{}\",\"window_snps\":{},\"window_het\":{},\"window_missing\":{},\"window_threshold\":{},\"max_density\":{},\"written\":{},\"scan\":{}}}",
            json_texto(&self.person),
            runs.join(","),
            self.count,
            self.total_length,
            self.longest,
            json_num(self.mean_length),
            self.autosome_length,
            json_num(self.f_roh),
            self.covered_span,
            json_num(self.f_roh_covered),
            self.snv_calls,
            self.het_calls,
            self.hom_calls,
            self.missing,
            self.outside_bed,
            self.non_snv,
            self.min_sites,
            self.min_length,
            self.max_gap,
            self.method,
            self.window_snps,
            self.window_het,
            self.window_missing,
            json_num(self.window_threshold),
            self.max_density,
            self.written,
            self.scan.json()
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!(
            "Runs of homozygosity of {}: {} run(s) of at least {} SNV calls and {} kb (gaps over {} kb break a run)\n",
            self.person,
            self.count,
            self.min_sites,
            self.min_length / 1000,
            self.max_gap / 1000
        );
        o.push_str(&format!(
            "  total {:.2} Mb, longest {:.2} Mb, mean {:.2} Mb; F_ROH {:.4} over {:.0} Mb of autosomes, {:.4} over the {:.0} Mb the file covers\n",
            self.total_length as f64 / 1e6,
            self.longest as f64 / 1e6,
            self.mean_length / 1e6,
            self.f_roh,
            self.autosome_length as f64 / 1e6,
            self.f_roh_covered,
            self.covered_span as f64 / 1e6
        ));
        o.push_str(&format!(
            "  {} autosomal SNV call(s) used: {} het, {} hom; {} missing, {} outside the callable BED, {} non-SNV record(s) ignored\n",
            milhar(self.snv_calls),
            milhar(self.het_calls),
            milhar(self.hom_calls),
            milhar(self.missing),
            milhar(self.outside_bed),
            milhar(self.non_snv)
        ));
        let mut maiores: Vec<&RohRun> = self.runs.iter().collect();
        maiores.sort_by_key(|r| std::cmp::Reverse(r.length()));
        for r in maiores.iter().take(5) {
            o.push_str(&format!(
                "    {}:{}-{}  {:.2} Mb, {} sites\n",
                r.chrom,
                milhar(r.start + 1),
                milhar(r.end),
                r.length() as f64 / 1e6,
                milhar(r.sites)
            ));
        }
        o.push_str(&self.scan.nota());
        if self.method == "window" {
            o.push_str(&format!(
                "  Method: sliding windows of {} SNV calls, as plink --homozyg: a window is homozygous with at most {} het and {} missing call(s); a call is in a run when at least {:.0}% of the windows over it are; runs sparser than one call per {} kb are dropped. For teaching and research.\n",
                self.window_snps,
                self.window_het,
                self.window_missing,
                100.0 * self.window_threshold,
                self.max_density / 1000
            ));
        } else {
            o.push_str("  Method: strict, a run is a stretch with no heterozygous call (no allele frequencies, no error model): one wrong het call splits it. For teaching and research.\n");
        }
        o
    }
}

struct Corrida {
    crom: u32,
    nome: String,
    ini: u64,
    fim: u64,
    n: u64,
}

/// SNV record: one-base REF, every ALT one base, all A/C/G/T.
fn e_snv(r: &Reg) -> bool {
    let rf = r.refs();
    rf.len() == 1 && e_acgt(rf) && r.alts().split(',').all(|a| a.len() == 1 && e_acgt(a) && !a.eq_ignore_ascii_case(rf))
}

/// Opens one person's file as a single stream (no merge).
fn abre_pessoa(person: &Person, cromos: &mut Cromos) -> Result<(Fluxo, Option<Chamavel>, String), String> {
    let mut f = Fluxo::abre(&person.path)?;
    let a = resolve_amostra(&f, person)?;
    f.colunas.push(9 + a);
    f.ordem = vec![0];
    let bed = match &person.callable_bed {
        Some(b) => Some(Chamavel::carrega(b, cromos)?),
        None => None,
    };
    let rotulo = person.label.clone().unwrap_or_else(|| f.amostras[a].clone());
    f.avanca(cromos)?;
    Ok((f, bed, rotulo))
}

fn contagens_fluxo(f: &Fluxo) -> ScanCounts {
    ScanCounts { records: f.registros, other_sequences: f.outras_seqs, malformed: f.malformados, ..Default::default() }
}

/// Runs of homozygosity of one person: stretches of the autosomes that are
/// homozygous by the sliding-window method of `plink --homozyg` (default), or,
/// with `window_snps = 0`, with no heterozygous SNV call at all; at least
/// `min_sites` SNV calls and `min_length` bases long, broken by stretches
/// without data longer than `max_gap`, and limited to the callable BED when
/// the person has one. With `output_bed`, the runs are written as BED:
/// chromosome, start, end, number of sites.
pub fn runs_of_homozygosity(person: &Person, output_bed: Option<&str>, opts: &RohOptions) -> Result<RohReport, String> {
    if opts.min_sites == 0 {
        return Err("min_sites must be at least 1.".to_string());
    }
    if !(0.0..=1.0).contains(&opts.window_threshold) {
        return Err(format!("window_threshold is {}: give a share between 0 and 1 (plink's default is 0.05).", opts.window_threshold));
    }
    let janela = opts.window_snps > 0;
    if let Some(s) = output_bed {
        let mut ent = vec![person.path.as_str()];
        if let Some(b) = &person.callable_bed {
            ent.push(b);
        }
        recusa_sobrescrever_entrada(s, &ent)?;
    }
    let mut cromos = Cromos::default();
    let (mut f, mut bed, rotulo) = abre_pessoa(person, &mut cromos)?;
    let mut rep = RohReport {
        person: rotulo,
        min_sites: opts.min_sites,
        min_length: opts.min_length,
        max_gap: opts.max_gap,
        autosome_length: opts.autosome_length.unwrap_or_else(grch38_autosome_length),
        method: if janela { "window" } else { "strict" },
        window_snps: opts.window_snps,
        window_het: opts.window_het,
        window_missing: opts.window_missing,
        window_threshold: opts.window_threshold,
        max_density: opts.max_density,
        ..Default::default()
    };
    // Método das janelas: as chamadas de um trecho sem lacuna (um cromossomo,
    // nenhum vão maior que max_gap), processadas quando o trecho fecha.
    let mut trecho: Vec<(u64, u8)> = Vec::new();
    let mut trecho_crom: Option<(u32, String)> = None;
    let mut aberta: Option<Corrida> = None;
    let mut span: Option<(u32, u64, u64)> = None;
    let fecha = |c: Corrida, rep: &mut RohReport| {
        if c.n >= opts.min_sites && c.fim - c.ini + 1 >= opts.min_length {
            rep.runs.push(RohRun { chrom: c.nome, start: c.ini - 1, end: c.fim, sites: c.n });
        }
    };
    while let Some(r) = f.cab.take() {
        let tipo = cromos.tipo(r.crom);
        let mut usa = tipo == TipoSeq::Autossomo;
        if usa && !e_snv(&r) {
            rep.non_snv += 1;
            usa = false;
        }
        let mut het = None;
        if usa {
            if let Some(b) = bed.as_mut() {
                if !b.contem(r.crom, r.pos) {
                    rep.outside_bed += 1;
                    usa = false;
                }
            }
        }
        if usa {
            het = match r.gts[0] {
                _ if opts.pass_only && !r.passa => None,
                Gt::Dip(a, b) => Some(a != b),
                _ => None,
            };
            if het.is_none() {
                rep.missing += 1;
            }
            if janela {
                let estado = match het {
                    Some(false) => 0,
                    Some(true) => 1,
                    None => 2,
                };
                let quebra = match (&trecho_crom, trecho.last()) {
                    (Some((c, _)), Some(&(ultima, _))) if *c == r.crom => {
                        let lacuna = match &bed {
                            Some(b) => b.maior_lacuna(r.crom, ultima, r.pos),
                            None => r.pos - ultima,
                        };
                        lacuna > opts.max_gap
                    }
                    _ => true,
                };
                if quebra {
                    if let Some((_, nome)) = trecho_crom.take() {
                        roh_janelas(&trecho, &nome, opts, &mut rep.runs);
                    }
                    trecho.clear();
                    trecho_crom = Some((r.crom, r.crom_bruto().to_string()));
                }
                trecho.push((r.pos, estado));
            }
        }
        if let Some(het) = het {
            rep.snv_calls += 1;
            match &mut span {
                Some((c, _, fim)) if *c == r.crom => *fim = r.pos,
                _ => {
                    if let Some((_, a, b)) = span {
                        rep.covered_span += b - a + 1;
                    }
                    span = Some((r.crom, r.pos, r.pos));
                }
            }
            if let Some(c) = aberta.take() {
                let quebra = c.crom != r.crom || {
                    let lacuna = match &bed {
                        Some(b) => b.maior_lacuna(r.crom, c.fim, r.pos),
                        None => r.pos - c.fim,
                    };
                    lacuna > opts.max_gap
                };
                if quebra || het {
                    fecha(c, &mut rep);
                } else {
                    aberta = Some(c);
                }
            }
            if het {
                rep.het_calls += 1;
            } else if janela {
                rep.hom_calls += 1;
            } else {
                rep.hom_calls += 1;
                match &mut aberta {
                    Some(c) => {
                        c.fim = r.pos;
                        c.n += 1;
                    }
                    None => {
                        aberta = Some(Corrida { crom: r.crom, nome: r.crom_bruto().to_string(), ini: r.pos, fim: r.pos, n: 1 })
                    }
                }
            }
        }
        f.livres.push(r);
        f.avanca(&mut cromos)?;
    }
    if let Some(c) = aberta.take() {
        fecha(c, &mut rep);
    }
    if let Some((_, nome)) = trecho_crom.take() {
        roh_janelas(&trecho, &nome, opts, &mut rep.runs);
    }
    if let Some((_, a, b)) = span {
        rep.covered_span += b - a + 1;
    }
    rep.count = rep.runs.len() as u64;
    rep.total_length = rep.runs.iter().map(RohRun::length).sum();
    rep.longest = rep.runs.iter().map(RohRun::length).max().unwrap_or(0);
    rep.mean_length = if rep.count > 0 { rep.total_length as f64 / rep.count as f64 } else { 0.0 };
    rep.f_roh = rep.total_length as f64 / rep.autosome_length.max(1) as f64;
    rep.f_roh_covered = if rep.covered_span > 0 { rep.total_length as f64 / rep.covered_span as f64 } else { 0.0 };
    rep.scan = contagens_fluxo(&f);
    if let Some(s) = output_bed {
        let mut guarda = Saida { caminho: Some(s.to_string()) };
        let mut w = Escritor::cria_para(s, Formato::Bed)?;
        for r in &rep.runs {
            w.linha(&format!("{}\t{}\t{}\t{}", r.chrom, r.start, r.end, r.sites))?;
            rep.written += 1;
        }
        w.termina()?;
        guarda.desarma();
    }
    Ok(rep)
}

// =====================================================================
//  4) SEXO GENETICO
// =====================================================================

/// What `genetic_sex` found. Counts are SNV calls with a genotype.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SexReport {
    pub person: String,
    /// "GRCh38" or "GRCh37": which pseudoautosomal regions were used.
    pub build: &'static str,
    /// "argument", "header" (chrX length in `##contig`) or "default".
    pub build_source: &'static str,
    pub has_x_data: bool,
    pub has_y_data: bool,
    pub autosome_het: u64,
    pub autosome_hom_alt: u64,
    pub autosome_het_fraction: Option<f64>,
    /// chrX outside the PARs.
    pub x_het: u64,
    pub x_hom_alt: u64,
    pub x_hom_ref: u64,
    /// Haploid calls (the caller already treated X as a single copy).
    pub x_haploid: u64,
    pub x_par_calls: u64,
    /// x_het / (x_het + x_hom_alt + haploid ALT calls).
    pub x_het_fraction: Option<f64>,
    /// x_het_fraction / autosome_het_fraction: ~1 with two X, ~0 with one.
    pub x_het_ratio: Option<f64>,
    /// chrY outside the PARs: calls with a genotype, non-reference calls,
    /// and heterozygous calls (not expected on a Y).
    pub y_calls: u64,
    pub y_non_ref: u64,
    pub y_het: u64,
    /// "XX", "XY" or "undetermined".
    pub verdict: &'static str,
    /// "x0_pattern", "xxy_pattern", "intermediate_x_heterozygosity",
    /// "haploid_x_calls", "no_x_data", "no_y_data".
    pub flags: Vec<&'static str>,
    pub reason: String,
    pub scan: ScanCounts,
}

impl SexReport {
    pub fn json(&self) -> String {
        let flags: Vec<String> = self.flags.iter().map(|f| format!("\"{f}\"")).collect();
        format!(
            "{{\"person\":{},\"build\":\"{}\",\"build_source\":\"{}\",\"has_x_data\":{},\"has_y_data\":{},\"autosome_het\":{},\"autosome_hom_alt\":{},\"autosome_het_fraction\":{},\"x_het\":{},\"x_hom_alt\":{},\"x_hom_ref\":{},\"x_haploid\":{},\"x_par_calls\":{},\"x_het_fraction\":{},\"x_het_ratio\":{},\"y_calls\":{},\"y_non_ref\":{},\"y_het\":{},\"verdict\":\"{}\",\"flags\":[{}],\"reason\":{},\"scan\":{}}}",
            json_texto(&self.person),
            self.build,
            self.build_source,
            self.has_x_data,
            self.has_y_data,
            self.autosome_het,
            self.autosome_hom_alt,
            json_opt(self.autosome_het_fraction),
            self.x_het,
            self.x_hom_alt,
            self.x_hom_ref,
            self.x_haploid,
            self.x_par_calls,
            json_opt(self.x_het_fraction),
            json_opt(self.x_het_ratio),
            self.y_calls,
            self.y_non_ref,
            self.y_het,
            self.verdict,
            flags.join(","),
            json_texto(&self.reason),
            self.scan.json()
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!("Genetic sex of {}: {} ({})\n", self.person, self.verdict, self.reason);
        o.push_str(&format!(
            "  chrX outside the PARs ({}, from the {}): {} het, {} hom-alt, {} haploid; het fraction {} against {} on the autosomes (ratio {})\n",
            self.build,
            self.build_source,
            milhar(self.x_het),
            milhar(self.x_hom_alt),
            milhar(self.x_haploid),
            fmt_opt(self.x_het_fraction, 3),
            fmt_opt(self.autosome_het_fraction, 3),
            fmt_opt(self.x_het_ratio, 3)
        ));
        o.push_str(&format!(
            "  chrY outside the PARs: {} call(s), {} non-reference, {} heterozygous\n",
            milhar(self.y_calls),
            milhar(self.y_non_ref),
            milhar(self.y_het)
        ));
        if !self.flags.is_empty() {
            o.push_str(&format!("  flags: {}\n", self.flags.join(", ")));
        }
        o.push_str("  Chromosomal sex from genotypes only; it says nothing about gender, and it is for teaching and research, not a clinical test.\n");
        o
    }
}

fn le_build(b: &str) -> Result<&'static str, String> {
    match b.to_ascii_lowercase().as_str() {
        "grch38" | "hg38" | "38" => Ok("GRCh38"),
        "grch37" | "hg19" | "b37" | "37" => Ok("GRCh37"),
        _ => Err(format!("Unknown genome build `{b}`: use \"GRCh38\" (hg38) or \"GRCh37\" (hg19, b37).")),
    }
}

/// Chromosomal sex of one person from the heterozygosity of chrX outside the
/// pseudoautosomal regions, compared with the autosomes, and the calls on
/// chrY. `build`: "GRCh38" or "GRCh37"; `None` reads the chrX length in the
/// header, and uses GRCh38 when the header does not say.
pub fn genetic_sex(person: &Person, build: Option<&str>, pass_only: bool) -> Result<SexReport, String> {
    let mut cromos = Cromos::default();
    let (mut f, mut bed, rotulo) = abre_pessoa(person, &mut cromos)?;
    let (nome_build, fonte) = match build {
        Some(b) => (le_build(b)?, "argument"),
        None => match f.contig_x {
            Some(CHRX_GRCH37) => ("GRCh37", "header"),
            Some(CHRX_GRCH38) => ("GRCh38", "header"),
            _ => ("GRCh38", "default"),
        },
    };
    let pars = if nome_build == "GRCh37" { PAR_GRCH37 } else { PAR_GRCH38 };
    let mut rep = SexReport { person: rotulo, build: nome_build, build_source: fonte, ..Default::default() };
    let mut x_hap_alt = 0u64;
    let mut y_hom = 0u64;
    let mut x_registros = 0u64;
    let mut y_registros = 0u64;
    while let Some(r) = f.cab.take() {
        let tipo = cromos.tipo(r.crom);
        match tipo {
            TipoSeq::X => x_registros += 1,
            TipoSeq::Y => y_registros += 1,
            _ => {}
        }
        let dentro = match bed.as_mut() {
            Some(b) => b.contem(r.crom, r.pos),
            None => true,
        };
        let g = if (pass_only && !r.passa) || !dentro { Gt::Falta } else { r.gts[0] };
        if g != Gt::Falta && e_snv(&r) {
            match tipo {
                TipoSeq::Autossomo => match g {
                    Gt::Dip(a, b) if a != b => rep.autosome_het += 1,
                    Gt::Dip(a, _) if a > 0 => rep.autosome_hom_alt += 1,
                    _ => {}
                },
                TipoSeq::X => {
                    if na_par(&pars.x, r.pos) {
                        rep.x_par_calls += 1;
                    } else {
                        match g {
                            Gt::Dip(a, b) if a != b => rep.x_het += 1,
                            Gt::Dip(a, _) if a > 0 => rep.x_hom_alt += 1,
                            Gt::Dip(_, _) => rep.x_hom_ref += 1,
                            Gt::Hap(a) => {
                                rep.x_haploid += 1;
                                if a > 0 {
                                    x_hap_alt += 1;
                                }
                            }
                            Gt::Falta => {}
                        }
                    }
                }
                TipoSeq::Y => {
                    if !na_par(&pars.y, r.pos) {
                        rep.y_calls += 1;
                        match g {
                            Gt::Dip(a, b) if a != b => {
                                rep.y_het += 1;
                                rep.y_non_ref += 1;
                            }
                            Gt::Dip(a, _) | Gt::Hap(a) if a > 0 => {
                                rep.y_non_ref += 1;
                                y_hom += 1;
                            }
                            _ => {}
                        }
                    }
                }
            }
        }
        f.livres.push(r);
        f.avanca(&mut cromos)?;
    }
    rep.scan = contagens_fluxo(&f);

    let x_nao_ref = rep.x_het + rep.x_hom_alt + x_hap_alt;
    let auto_nao_ref = rep.autosome_het + rep.autosome_hom_alt;
    rep.has_x_data = x_nao_ref >= MIN_X_CALLS;
    rep.has_y_data = rep.y_calls > 0;
    rep.x_het_fraction = (x_nao_ref > 0).then(|| rep.x_het as f64 / x_nao_ref as f64);
    rep.autosome_het_fraction = (auto_nao_ref > 0).then(|| rep.autosome_het as f64 / auto_nao_ref as f64);
    let (ref_auto, nota_auto) = match rep.autosome_het_fraction {
        Some(h) if auto_nao_ref >= 1000 && h > 0.0 => (h, ""),
        _ => (AUTOSOME_HET_FRACTION_FALLBACK, " (too few autosomal calls to measure: a typical 0.6 was assumed)"),
    };
    rep.x_het_ratio = rep.x_het_fraction.map(|h| h / ref_auto);
    let y_presente = y_hom >= MIN_Y_CALLS;
    if !rep.has_y_data {
        rep.flags.push("no_y_data");
    }
    if x_nao_ref > 0 && rep.x_haploid * 2 > x_nao_ref {
        rep.flags.push("haploid_x_calls");
    }
    let texto_y = if y_presente {
        format!("{} homozygous/haploid non-reference calls on chrY", milhar(y_hom))
    } else if rep.has_y_data {
        format!("only {} homozygous/haploid non-reference call(s) on chrY", milhar(y_hom))
    } else if f.contig_y || y_registros > 0 {
        "no chrY call with a genotype".to_string()
    } else {
        "the file has no chrY data".to_string()
    };
    if !rep.has_x_data {
        rep.flags.push("no_x_data");
        rep.verdict = "undetermined";
        rep.reason = if x_registros == 0 {
            format!("the file has no chrX data outside the pseudoautosomal regions; {texto_y}")
        } else {
            format!(
                "only {} non-reference SNV call(s) on chrX outside the pseudoautosomal regions (at least {MIN_X_CALLS} are needed); {texto_y}",
                x_nao_ref
            )
        };
        return Ok(rep);
    }
    let ratio = rep.x_het_ratio.unwrap_or(0.0);
    let texto_x = format!(
        "chrX heterozygosity is {:.2} of the autosomal one{nota_auto}",
        ratio
    );
    if ratio < ONE_X_RATIO {
        if y_presente {
            rep.verdict = "XY";
            rep.reason = format!("{texto_x} (one X) and there are {texto_y}");
        } else {
            rep.verdict = "undetermined";
            rep.flags.push("x0_pattern");
            rep.reason = format!(
                "{texto_x} (one X), and {texto_y}: the pattern of X0, or of a male sample in a file without chrY calls"
            );
        }
    } else if ratio > TWO_X_RATIO {
        if y_presente {
            rep.verdict = "undetermined";
            rep.flags.push("xxy_pattern");
            rep.reason = format!(
                "{texto_x} (two X) and there are {texto_y}: the pattern of XXY, of a mixture of two samples, or of Y reads mapped by mistake"
            );
        } else {
            rep.verdict = "XX";
            rep.reason = format!("{texto_x} (two X) and {texto_y}");
        }
    } else {
        rep.verdict = "undetermined";
        rep.flags.push("intermediate_x_heterozygosity");
        rep.reason = format!(
            "{texto_x}, between one X (< {ONE_X_RATIO}) and two (> {TWO_X_RATIO}): possible mosaicism, contamination or poor chrX calls; {texto_y}"
        );
    }
    Ok(rep)
}
