// =====================================================================
//  VALIDADOR DE ARQUIVOS GENOMICOS
//
//  Uma passada em fluxo que diz, em linguagem de gente e com o numero da
//  linha, o que ha de errado num arquivo -- ANTES de ele entrar numa analise
//  de horas. VCF, BED, GFF/GTF, SAM, FASTA e FASTQ, de qualquer tamanho.
//
//  Regras da casa:
//
//   - Todo problema aponta a LINHA da primeira ocorrencia e diz o que esta
//     errado nela. "POS `abc` nao e um numero inteiro", nao "parse error".
//   - Problemas iguais sao agregados: um arquivo com um milhao de linhas
//     tortas gera uma entrada com a contagem, nao um milhao de entradas.
//   - Quando existe conserto, o problema diz qual funcao conserta.
//   - ERRO e o que faz uma ferramenta falhar ou, pior, dar resultado errado
//     em silencio. AVISO e o que merece um olhar. INFO e so um fato.
// =====================================================================

use crate::externo::abre_leitura;
use crate::nomes::{convencao, json_texto, Convencao};
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tipo {
    Vcf,
    Bed,
    Gff,
    Sam,
    Fasta,
    Fastq,
}

impl Tipo {
    pub fn de_texto(s: &str) -> Option<Tipo> {
        match s.to_ascii_lowercase().as_str() {
            "vcf" => Some(Tipo::Vcf),
            "bed" => Some(Tipo::Bed),
            "gff" | "gff3" | "gtf" => Some(Tipo::Gff),
            "sam" | "bam" => Some(Tipo::Sam),
            "fasta" | "fa" | "fna" | "faa" | "ffn" | "fas" => Some(Tipo::Fasta),
            "fastq" | "fq" => Some(Tipo::Fastq),
            _ => None,
        }
    }

    /// Deduz o tipo pela extensao, ignorando a do envelope: `.gz`, `.bgz`,
    /// `.lz4` e `.parquet`. Um `.bam` e lido como o SAM que ele guarda.
    pub fn pelo_nome(caminho: &str) -> Option<Tipo> {
        let c = caminho.to_ascii_lowercase();
        let c = [".gz", ".bgz", ".lz4", ".parquet"].iter().find_map(|s| c.strip_suffix(s)).unwrap_or(&c);
        Tipo::de_texto(c.rsplit('.').next()?)
    }

    pub fn nome(self) -> &'static str {
        match self {
            Tipo::Vcf => "vcf",
            Tipo::Bed => "bed",
            Tipo::Gff => "gff",
            Tipo::Sam => "sam",
            Tipo::Fasta => "fasta",
            Tipo::Fastq => "fastq",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Gravidade {
    Erro,
    Aviso,
    Info,
}

impl Gravidade {
    pub fn nome(self) -> &'static str {
        match self {
            Gravidade::Erro => "error",
            Gravidade::Aviso => "warning",
            Gravidade::Info => "info",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Problema {
    /// Identificador estavel do tipo de problema, para quem trata por codigo.
    pub codigo: &'static str,
    pub gravidade: Gravidade,
    /// Linha da primeira ocorrencia (a partir de 1). Zero: vale para o arquivo.
    pub linha: usize,
    /// O que esta errado na primeira ocorrencia.
    pub mensagem: String,
    pub ocorrencias: usize,
    /// Como consertar, quando ha conserto.
    pub correcao: Option<&'static str>,
}

#[derive(Clone, Debug)]
pub struct Relatorio {
    pub caminho: String,
    pub tipo: Tipo,
    pub linhas: usize,
    pub registros: usize,
    /// `false` quando a validacao parou no limite de linhas pedido.
    pub completo: bool,
    /// Do mais grave para o menos grave; dentro de cada gravidade, por linha.
    pub problemas: Vec<Problema>,
    /// Nomes de sequencia distintos (formatos posicionais).
    pub sequencias: usize,
    /// "chr", "no_chr", "mixed", "other" ou "-" (formato sem cromossomos).
    pub nomenclatura: &'static str,
    /// O arquivo esta ordenado por sequencia e posicao? `None` em FASTA/FASTQ.
    pub ordenado: Option<bool>,
}

impl Relatorio {
    pub fn conta(&self, g: Gravidade) -> usize {
        self.problemas.iter().filter(|p| p.gravidade == g).count()
    }

    /// Sem nenhum ERRO.
    pub fn valido(&self) -> bool {
        self.conta(Gravidade::Erro) == 0
    }

    /// O relatorio em texto, pronto para mostrar.
    pub fn resumo(&self) -> String {
        let mut o = format!(
            "{}: {} with {} lines and {} records{} -- {} error(s), {} warning(s)\n",
            self.caminho,
            self.tipo.nome().to_ascii_uppercase(),
            self.linhas,
            self.registros,
            if self.completo { "" } else { " (only the beginning was validated)" },
            self.conta(Gravidade::Erro),
            self.conta(Gravidade::Aviso)
        );
        for p in &self.problemas {
            let onde = if p.linha > 0 { format!("line {}: ", p.linha) } else { String::new() };
            let mais = if p.ocorrencias > 1 { format!(" [{} occurrences]", p.ocorrencias) } else { String::new() };
            o.push_str(&format!("  {:<7} {onde}{}{mais}\n", p.gravidade.nome().to_ascii_uppercase(), p.mensagem));
            if let Some(c) = p.correcao {
                o.push_str(&format!("          fix: {c}\n"));
            }
        }
        if self.problemas.is_empty() {
            o.push_str("  no problems found\n");
        }
        o
    }

    pub fn json(&self) -> String {
        let ps: Vec<String> = self
            .problemas
            .iter()
            .map(|p| {
                format!(
                    "{{\"code\":{},\"severity\":{},\"line\":{},\"message\":{},\"occurrences\":{},\"fix\":{}}}",
                    json_texto(p.codigo),
                    json_texto(p.gravidade.nome()),
                    p.linha,
                    json_texto(&p.mensagem),
                    p.ocorrencias,
                    p.correcao.map(json_texto).unwrap_or_else(|| "null".to_string())
                )
            })
            .collect();
        format!(
            "{{\"path\":{},\"format\":{},\"lines\":{},\"records\":{},\"complete\":{},\"valid\":{},\"errors\":{},\"warnings\":{},\"sequences\":{},\"naming\":{},\"sorted\":{},\"problems\":[{}]}}",
            json_texto(&self.caminho),
            json_texto(self.tipo.nome()),
            self.linhas,
            self.registros,
            self.completo,
            self.valido(),
            self.conta(Gravidade::Erro),
            self.conta(Gravidade::Aviso),
            self.sequencias,
            json_texto(self.nomenclatura),
            match self.ordenado {
                Some(v) => v.to_string(),
                None => "null".to_string(),
            },
            ps.join(",")
        )
    }
}

const CONSERTO_ORDENAR: &str = "sort_genomic(...) sorts the file, whatever its size";
const CONSERTO_NOMES: &str = "normalize_chromosomes(...) rewrites the file in a single naming convention";
const CONSERTO_TAB: &str = "replace the spaces between the columns with tabs";

/// Ate quantos identificadores de FASTA guardar para achar repetidos.
const MAX_IDS: usize = 2_000_000;

struct Estado {
    tipo: Tipo,
    linha: usize,
    registros: usize,
    lista: Vec<Problema>,
    onde: HashMap<&'static str, usize>,
    primeira_util: bool,
    /// Linhas vazias vistas e ainda nao seguidas de conteudo (as do fim do
    /// arquivo nao sao problema).
    vazias: (usize, usize),
    crlf: usize,

    // posicional
    vistos: HashSet<String>,
    ultima: String,
    ultimo_inicio: u64,
    ordenado: bool,
    ex_chr: Option<(usize, String)>,
    ex_sem: Option<(usize, String)>,

    // VCF
    viu_chrom: bool,
    colunas_cab: usize,
    contigs: HashSet<String>,
    // SAM
    sq: HashSet<String>,
    so: Option<String>,
    // GFF
    secao_fasta: bool,
    // FASTA
    ids: HashSet<String>,
    bases: u64,
    linha_cab: usize,
    id_atual: String,
    // FASTQ
    fase: u8,
    len_seq: usize,
    min_qual: u8,
}

impl Estado {
    fn registra(
        &mut self,
        codigo: &'static str,
        gravidade: Gravidade,
        linha: usize,
        correcao: Option<&'static str>,
        mensagem: impl FnOnce() -> String,
    ) {
        match self.onde.get(codigo) {
            Some(&i) => self.lista[i].ocorrencias += 1,
            None => {
                self.onde.insert(codigo, self.lista.len());
                self.lista.push(Problema { codigo, gravidade, linha, mensagem: mensagem(), ocorrencias: 1, correcao });
            }
        }
    }

    fn erro(&mut self, codigo: &'static str, correcao: Option<&'static str>, mensagem: impl FnOnce() -> String) {
        let l = self.linha;
        self.registra(codigo, Gravidade::Erro, l, correcao, mensagem);
    }

    fn aviso(&mut self, codigo: &'static str, correcao: Option<&'static str>, mensagem: impl FnOnce() -> String) {
        let l = self.linha;
        self.registra(codigo, Gravidade::Aviso, l, correcao, mensagem);
    }

    /// Acompanha ordem e nomenclatura de um registro posicional.
    fn posicional(&mut self, seq: &str, inicio: u64) {
        match convencao(seq) {
            Convencao::ComChr if self.ex_chr.is_none() => self.ex_chr = Some((self.linha, seq.to_string())),
            Convencao::SemChr if self.ex_sem.is_none() => self.ex_sem = Some((self.linha, seq.to_string())),
            _ => {}
        }
        if seq != self.ultima {
            if !self.vistos.insert(seq.to_string()) {
                self.ordenado = false;
                let l = self.linha;
                self.registra("unsorted", Gravidade::Info, l, Some(CONSERTO_ORDENAR), || {
                    format!("the file is not sorted: sequence `{seq}` shows up again here, after other sequences")
                });
            }
            self.ultima.clear();
            self.ultima.push_str(seq);
        } else if inicio < self.ultimo_inicio {
            self.ordenado = false;
            let (l, ant) = (self.linha, self.ultimo_inicio);
            self.registra("unsorted", Gravidade::Info, l, Some(CONSERTO_ORDENAR), || {
                format!("the file is not sorted: position {inicio} comes after position {ant} on sequence `{seq}`")
            });
        }
        self.ultimo_inicio = inicio;
    }

    /// O registro tem colunas de menos. Diz se a causa provavel e espaco no
    /// lugar de tabulacao -- o erro mais comum de quem editou o arquivo a mao.
    fn colunas_de_menos(&mut self, codigo: &'static str, t: &str, tem: usize, precisa: &str) {
        let por_espaco = t.split_whitespace().count();
        if por_espaco > tem {
            self.erro(codigo, Some(CONSERTO_TAB), || {
                format!("has {tem} tab-separated column(s) but {por_espaco} when split on spaces: the columns must be separated by TABS ({precisa})")
            });
        } else {
            self.erro(codigo, None, || format!("has {tem} column(s); {precisa}"));
        }
    }

    fn vcf(&mut self, t: &str) {
        if self.primeira_util && !t.starts_with("##fileformat=VCF") {
            self.aviso("vcf_missing_fileformat", None, || {
                "the first line should be `##fileformat=VCFv4.x`; some tools reject the file without it".to_string()
            });
        }
        if t.starts_with("##") {
            if let Some(c) = t.strip_prefix("##contig=<") {
                if let Some(p) = c.find("ID=") {
                    let r = &c[p + 3..];
                    self.contigs.insert(r[..r.find([',', '>']).unwrap_or(r.len())].to_string());
                }
            }
            return;
        }
        if t.starts_with("#CHROM") {
            self.viu_chrom = true;
            self.colunas_cab = t.split('\t').count();
            if self.colunas_cab < 8 {
                self.colunas_de_menos("vcf_header_columns", t, self.colunas_cab, "the `#CHROM` line needs 8 columns");
            }
            return;
        }
        if t.starts_with('#') {
            return;
        }
        self.registros += 1;
        if !self.viu_chrom {
            self.erro("vcf_missing_chrom", None, || {
                "there are variants before the column-title line `#CHROM POS ID REF ALT QUAL FILTER INFO`, or that line is missing".to_string()
            });
        }
        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() < 8 {
            self.colunas_de_menos("vcf_columns", t, cols.len(), "a VCF needs at least 8: CHROM POS ID REF ALT QUAL FILTER INFO");
            return;
        }
        if self.colunas_cab >= 8 && cols.len() != self.colunas_cab {
            let cab = self.colunas_cab;
            self.erro("vcf_column_count_mismatch", None, || {
                format!("has {} columns, but the `#CHROM` line declares {cab}", cols.len())
            });
        }
        let pos = match cols[1].parse::<u64>() {
            Ok(p) => Some(p),
            Err(_) => {
                self.erro("vcf_pos", None, || format!("POS `{}` is not an integer", cols[1]));
                None
            }
        };
        if let Some(c) = cols[3].chars().find(|c| !matches!(c.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T' | 'N')) {
            self.erro("vcf_ref", None, || format!("REF `{}` has the invalid character `{c}` (only A, C, G, T and N are allowed)", cols[3]));
        }
        // ALT simbolico (`<DEL>`) e de quebra (`G]17:198982]`) tem gramatica
        // propria; aqui so se confere o alelo simples.
        let alt = cols[4];
        if !alt.contains(['<', '[', ']']) {
            if let Some(c) = alt.chars().find(|c| !matches!(c.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T' | 'N' | ',' | '.' | '*')) {
                self.erro("vcf_alt", None, || format!("ALT `{alt}` has the invalid character `{c}`"));
            }
        }
        if cols[5] != "." && cols[5].parse::<f64>().is_err() {
            self.erro("vcf_qual", None, || format!("QUAL `{}` is neither a number nor `.`", cols[5]));
        }
        if !self.contigs.is_empty() && !self.contigs.contains(cols[0]) {
            self.aviso("vcf_contig", None, || {
                format!("sequence `{}` is not declared in any `##contig` line of the header", cols[0])
            });
        }
        if let Some(p) = pos {
            self.posicional(cols[0], p);
        }
    }

    fn bed(&mut self, t: &str) {
        if t.starts_with('#') || t.starts_with("track") || t.starts_with("browser") {
            return;
        }
        self.registros += 1;
        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() < 3 {
            self.colunas_de_menos("bed_columns", t, cols.len(), "a BED needs at least 3: chromosome, start, end");
            return;
        }
        let (Ok(ini), Ok(fim)) = (cols[1].parse::<u64>(), cols[2].parse::<u64>()) else {
            self.erro("bed_coordinate", None, || {
                format!("start `{}` and end `{}` must be non-negative integers", cols[1], cols[2])
            });
            return;
        };
        if ini > fim {
            self.erro("bed_inverted", None, || format!("the start ({ini}) is greater than the end ({fim})"));
        }
        self.posicional(cols[0], ini);
    }

    fn gff(&mut self, t: &str) {
        if self.secao_fasta {
            return;
        }
        if t.starts_with("##FASTA") {
            self.secao_fasta = true;
            let l = self.linha;
            self.registra("gff_fasta", Gravidade::Info, l, None, || {
                "from here on the file carries the sequence itself (`##FASTA` section); it was not validated, and a GFF like this cannot be sorted line by line".to_string()
            });
            return;
        }
        if t.starts_with('#') {
            return;
        }
        self.registros += 1;
        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() != 9 {
            if cols.len() < 9 {
                self.colunas_de_menos("gff_columns", t, cols.len(), "a GFF has exactly 9");
            } else {
                self.erro("gff_columns", None, || format!("has {} columns; a GFF has exactly 9", cols.len()));
            }
            return;
        }
        let (Ok(ini), Ok(fim)) = (cols[3].parse::<u64>(), cols[4].parse::<u64>()) else {
            self.erro("gff_coordinate", None, || format!("start `{}` and end `{}` must be integers", cols[3], cols[4]));
            return;
        };
        if ini == 0 {
            self.erro("gff_zero_based", None, || {
                "the start is 0, but GFF coordinates start at 1 (BED is the format that starts at 0)".to_string()
            });
        }
        if ini > fim {
            self.erro("gff_inverted", None, || format!("the start ({ini}) is greater than the end ({fim})"));
        }
        if !matches!(cols[6], "+" | "-" | "." | "?") {
            self.erro("gff_strand", None, || format!("strand `{}` is not `+`, `-`, `.` or `?`", cols[6]));
        }
        if cols[2] == "CDS" && !matches!(cols[7], "0" | "1" | "2") {
            self.aviso("gff_phase", None, || format!("a CDS feature needs phase 0, 1 or 2; here it is `{}`", cols[7]));
        }
        self.posicional(cols[0], ini);
    }

    fn sam(&mut self, t: &str) {
        if t.starts_with('@') {
            if self.registros > 0 {
                self.aviso("sam_late_header", None, || "header line (`@`) after the alignments".to_string());
            }
            let campo = |p: &str| t.split('\t').find_map(|c| c.strip_prefix(p));
            if t.starts_with("@HD") {
                self.so = campo("SO:").map(str::to_string);
            } else if t.starts_with("@SQ") {
                if let Some(sn) = campo("SN:") {
                    self.sq.insert(sn.to_string());
                }
            }
            return;
        }
        self.registros += 1;
        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() < 11 {
            self.colunas_de_menos("sam_columns", t, cols.len(), "a SAM alignment has 11 mandatory fields");
            return;
        }
        if cols[1].parse::<u16>().is_err() {
            self.erro("sam_flag", None, || format!("FLAG `{}` is not an integer from 0 to 65535", cols[1]));
        }
        let pos = cols[3].parse::<u64>();
        if pos.is_err() {
            self.erro("sam_pos", None, || format!("POS `{}` is not an integer", cols[3]));
        }
        if cols[4].parse::<u8>().is_err() {
            self.erro("sam_mapq", None, || format!("MAPQ `{}` is not an integer from 0 to 255", cols[4]));
        }
        if cols[2] != "*" && !self.sq.is_empty() && !self.sq.contains(cols[2]) {
            self.erro("sam_rname", None, || {
                format!("reference `{}` is not in any `@SQ` line of the header", cols[2])
            });
        }
        if cols[9] != "*" && cols[10] != "*" && cols[9].len() != cols[10].len() {
            self.erro("sam_quality", None, || {
                format!("the sequence has {} bases and the quality string has {} values", cols[9].len(), cols[10].len())
            });
        }
        if let (Ok(p), true) = (pos, cols[2] != "*") {
            self.posicional(cols[2], p);
        }
    }

    fn fecha_fasta(&mut self) {
        if self.registros > 0 && self.bases == 0 {
            let (l, id) = (self.linha_cab, self.id_atual.clone());
            self.registra("fasta_empty", Gravidade::Aviso, l, None, || format!("sequence `{id}` has no bases"));
        }
    }

    fn fasta(&mut self, t: &str) {
        if let Some(cab) = t.strip_prefix('>') {
            self.fecha_fasta();
            self.registros += 1;
            self.bases = 0;
            self.linha_cab = self.linha;
            let id = cab.split_whitespace().next().unwrap_or("");
            self.id_atual.clear();
            self.id_atual.push_str(id);
            if id.is_empty() {
                self.erro("fasta_missing_id", None, || "the `>` header has no identifier".to_string());
            } else if self.ids.len() < MAX_IDS && !self.ids.insert(id.to_string()) {
                self.erro("fasta_duplicate_id", None, || {
                    format!("identifier `{id}` already appeared earlier; every sequence needs a unique name")
                });
            }
            return;
        }
        if self.registros == 0 {
            self.erro("fasta_missing_header", None, || {
                "there is sequence before the first header; every FASTA record starts with a `>name` line".to_string()
            });
        }
        // Letras (nucleotideos, ambiguidade IUPAC, aminoacidos), `*` de parada
        // e `-`/`.` de lacuna. Digito, espaco e pontuacao nao sao sequencia.
        if let Some((col, c)) = t.chars().enumerate().find(|(_, c)| !(c.is_ascii_alphabetic() || matches!(c, '*' | '-' | '.'))) {
            self.erro("fasta_character", None, || {
                let visivel = if c == ' ' { "a space".to_string() } else { format!("`{c}`") };
                format!("invalid character {visivel} in the sequence (column {})", col + 1)
            });
        }
        self.bases += t.len() as u64;
    }

    fn fastq(&mut self, t: &str) {
        match self.fase {
            0 => {
                self.registros += 1;
                if !t.starts_with('@') {
                    self.erro("fastq_header", None, || {
                        format!(
                            "expected a record header (it starts with `@`), but the line starts with `{}`; from here on the 4-line records are out of step",
                            t.chars().next().unwrap_or(' ')
                        )
                    });
                }
            }
            1 => {
                self.len_seq = t.len();
                if let Some(c) = t.chars().find(|c| !(c.is_ascii_alphabetic() || *c == '.')) {
                    self.erro("fastq_base", None, || format!("invalid character `{c}` in the sequence"));
                }
            }
            2 => {
                if !t.starts_with('+') {
                    self.erro("fastq_separator", None, || {
                        "the third line of every record must start with `+`".to_string()
                    });
                }
            }
            _ => {
                if t.len() != self.len_seq {
                    let n = self.len_seq;
                    self.erro("fastq_length", None, || {
                        format!("the sequence has {n} bases and the quality line has {} values; there must be one per base", t.len())
                    });
                }
                for &b in t.as_bytes() {
                    if !(33..=126).contains(&b) {
                        self.erro("fastq_quality", None, || {
                            "the quality line has a character outside the valid range (ASCII 33 to 126)".to_string()
                        });
                        break;
                    }
                    self.min_qual = self.min_qual.min(b);
                }
            }
        }
        self.fase = (self.fase + 1) % 4;
    }

    fn processa(&mut self, t: &str) {
        if t.is_empty() {
            // No FASTQ uma linha vazia dentro do registro e conteudo (leitura
            // de tamanho zero); entre registros, e so uma linha vazia.
            if self.tipo != Tipo::Fastq || self.fase == 0 {
                if self.vazias.0 == 0 {
                    self.vazias.1 = self.linha;
                }
                self.vazias.0 += 1;
                return;
            }
        } else if self.vazias.0 > 0 {
            let (n, l) = self.vazias;
            self.vazias = (0, 0);
            for _ in 0..n {
                self.registra("blank_line", Gravidade::Aviso, l, None, || {
                    "blank line in the middle of the file; several tools stop at it or count it as a record".to_string()
                });
            }
        }
        match self.tipo {
            Tipo::Vcf => self.vcf(t),
            Tipo::Bed => self.bed(t),
            Tipo::Gff => self.gff(t),
            Tipo::Sam => self.sam(t),
            Tipo::Fasta => self.fasta(t),
            Tipo::Fastq => self.fastq(t),
        }
        self.primeira_util = false;
    }

    fn encerra(mut self, caminho: &str, completo: bool) -> Relatorio {
        let posicional = matches!(self.tipo, Tipo::Vcf | Tipo::Bed | Tipo::Gff | Tipo::Sam);
        if self.linha == 0 {
            self.registra("empty", Gravidade::Erro, 0, None, || "the file is empty".to_string());
        } else if self.registros == 0 {
            self.registra("no_records", Gravidade::Aviso, 0, None, || {
                "the file has only a header: no records".to_string()
            });
        }
        if self.crlf > 0 {
            let (n, total) = (self.crlf, self.linha);
            self.registra("crlf", Gravidade::Aviso, 0, Some("any sort or normalize operation writes its output with Unix line endings"), || {
                if n == total {
                    "every line ends in CRLF (Windows line ending); Linux tools see an extra character in the last column".to_string()
                } else {
                    format!("{n} of {total} lines end in CRLF (Windows line ending) and the others do not: the file mixes both")
                }
            });
            // `registra` contou uma ocorrencia; o numero certo e o de linhas.
            if let Some(&i) = self.onde.get("crlf") {
                self.lista[i].ocorrencias = n;
            }
        }
        match self.tipo {
            Tipo::Fasta if completo => self.fecha_fasta(),
            Tipo::Fastq => {
                if completo && self.fase != 0 {
                    let (l, faltam) = (self.linha, 4 - self.fase);
                    self.registra("fastq_truncated", Gravidade::Erro, l, None, || {
                        format!("the file ends in the middle of a record: {faltam} line(s) missing. An interrupted download or copy usually leaves this behind")
                    });
                }
                if self.registros > 0 && self.min_qual >= 64 && self.min_qual < u8::MAX {
                    self.registra("fastq_phred64", Gravidade::Aviso, 0, None, || {
                        "the qualities seem to be in the old Phred+64 scale (Illumina 1.3 to 1.7); current tools assume Phred+33 and would read inflated qualities".to_string()
                    });
                }
            }
            Tipo::Sam => {
                if self.registros > 0 && self.sq.is_empty() {
                    self.registra("sam_missing_sq", Gravidade::Aviso, 0, None, || {
                        "the header has no `@SQ` line; without them the file cannot be converted to BAM and its reference cannot be checked".to_string()
                    });
                }
                if self.so.as_deref() == Some("coordinate") && !self.ordenado {
                    self.registra("sam_so_false", Gravidade::Erro, 0, Some(CONSERTO_ORDENAR), || {
                        "the header declares `SO:coordinate`, but the alignments are not in coordinate order; anything that trusts the header will index or merge incorrectly".to_string()
                    });
                }
            }
            _ => {}
        }
        let nomenclatura = if !posicional {
            "-"
        } else {
            match (&self.ex_chr, &self.ex_sem) {
                (Some((la, a)), Some((lb, b))) => {
                    let (la, a, lb, b) = (*la, a.clone(), *lb, b.clone());
                    self.registra("mixed_names", Gravidade::Erro, la.min(lb), Some(CONSERTO_NOMES), || {
                        format!("the file mixes the two chromosome naming conventions: `{a}` (line {la}) and `{b}` (line {lb})")
                    });
                    "mixed"
                }
                (Some(_), None) => "chr",
                (None, Some(_)) => "no_chr",
                (None, None) => "other",
            }
        };
        self.lista.sort_by_key(|p| (p.gravidade, p.linha));
        Relatorio {
            caminho: caminho.to_string(),
            tipo: self.tipo,
            linhas: self.linha,
            registros: self.registros,
            completo,
            problemas: self.lista,
            sequencias: self.vistos.len(),
            nomenclatura,
            ordenado: if posicional { Some(self.ordenado) } else { None },
        }
    }
}

/// Valida um arquivo e devolve o relatorio.
///
/// `tipo`: sem ele, vem da extensao. `limite_linhas`: valida so o comeco do
/// arquivo -- util para uma conferencia rapida de um arquivo de dezenas de GB;
/// o relatorio avisa que nao viu tudo.
///
/// Le `.gz` (inclusive BGZF) e `.lz4` direto. So devolve `Err` quando o
/// arquivo nao pode ser aberto ou lido; problema de CONTEUDO vai no relatorio.
pub fn validar_arquivo(caminho: &str, tipo: Option<Tipo>, limite_linhas: Option<usize>) -> Result<Relatorio, String> {
    let tipo = tipo.or_else(|| Tipo::pelo_nome(caminho)).ok_or_else(|| {
        format!("Could not infer the format of `{caminho}` from its extension; pass the format (vcf, bed, gff, sam, fasta, fastq)")
    })?;
    let mut leitor = abre_leitura(caminho)?;
    let mut e = Estado {
        tipo,
        linha: 0,
        registros: 0,
        lista: Vec::new(),
        onde: HashMap::new(),
        primeira_util: true,
        vazias: (0, 0),
        crlf: 0,
        vistos: HashSet::new(),
        ultima: String::new(),
        ultimo_inicio: 0,
        ordenado: true,
        ex_chr: None,
        ex_sem: None,
        viu_chrom: false,
        colunas_cab: 0,
        contigs: HashSet::new(),
        sq: HashSet::new(),
        so: None,
        secao_fasta: false,
        ids: HashSet::new(),
        bases: 0,
        linha_cab: 0,
        id_atual: String::new(),
        fase: 0,
        len_seq: 0,
        min_qual: u8::MAX,
    };
    let mut bruto: Vec<u8> = Vec::new();
    let mut completo = true;
    loop {
        bruto.clear();
        if leitor.read_until(b'\n', &mut bruto).map_err(|x| format!("Error reading {caminho}: {x}"))? == 0 {
            break;
        }
        if limite_linhas.is_some_and(|l| e.linha >= l) {
            completo = false;
            break;
        }
        e.linha += 1;
        let mut fim = bruto.len();
        if fim > 0 && bruto[fim - 1] == b'\n' {
            fim -= 1;
        }
        if fim > 0 && bruto[fim - 1] == b'\r' {
            fim -= 1;
            e.crlf += 1;
        }
        match std::str::from_utf8(&bruto[..fim]) {
            Ok(t) => e.processa(t),
            Err(_) => {
                e.erro("not_text", None, || {
                    "the line has bytes that are not UTF-8 text; the file may be binary (BAM, BCF), corrupted, or in another encoding".to_string()
                });
                let t = String::from_utf8_lossy(&bruto[..fim]).into_owned();
                e.processa(&t);
            }
        }
    }
    Ok(e.encerra(caminho, completo))
}
