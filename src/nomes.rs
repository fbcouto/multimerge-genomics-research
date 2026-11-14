// =====================================================================
//  NOMES DE SEQUENCIA E GENOMA DE REFERENCIA
//
//  Os dois erros silenciosos mais caros da genomica:
//
//   1. NOMENCLATURA. Um arquivo chama o cromossomo de `chr1` (convencao
//      UCSC), o outro de `1` (Ensembl/NCBI). Nada casa, nenhuma ferramenta
//      reclama, e o resultado sai vazio depois de horas.
//   2. REFERENCIA. Um arquivo esta em GRCh37/hg19, o outro em GRCh38/hg38. Os
//      nomes sao iguais, as coordenadas nao: tudo casa, e casa errado.
//
//  Este modulo reconhece a convencao de um nome, converte entre as duas,
//  reescreve um arquivo inteiro na convencao pedida, e compara as referencias
//  de dois arquivos pelos TAMANHOS das sequencias que cada um declara.
// =====================================================================

use crate::externo::{abre_leitura, Escritor, Formato};
use std::collections::HashMap;
use std::io::BufRead;

// ---------------------------------------------------------------------
//  CONVENCOES
// ---------------------------------------------------------------------

/// `s` tem cara de nome de cromossomo primario, sem o prefixo? Numeros
/// (1, 22), sexuais (X, Y, W, Z), mitocondrial (M, MT), romanos (levedura:
/// I ... XVI) e bracos (drosofila: 2L, 3R).
///
/// Scaffolds e contigs nao montados (`GL000192.1`, `Un_KI270742v1`) ficam de
/// fora de proposito: o nome deles muda de forma entre as convencoes e nao
/// basta tirar ou por um prefixo.
pub fn parece_cromossomo(s: &str) -> bool {
    let b = s.as_bytes();
    if b.is_empty() || b.len() > 5 {
        return false;
    }
    if b.iter().all(u8::is_ascii_digit) {
        return b.len() <= 3;
    }
    if matches!(s, "X" | "Y" | "W" | "Z" | "M" | "MT") {
        return true;
    }
    if b.iter().all(|c| matches!(c, b'I' | b'V' | b'X')) {
        return true;
    }
    let (num, braco) = b.split_at(b.len() - 1);
    !num.is_empty() && num.iter().all(u8::is_ascii_digit) && matches!(braco[0], b'L' | b'R')
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Convencao {
    /// `chr1`, `chrX`, `chrM` -- UCSC.
    ComChr,
    /// `1`, `X`, `MT` -- Ensembl, NCBI.
    SemChr,
    /// Qualquer outra coisa: scaffold, acesso RefSeq (`NC_000001.11`), plasmideo.
    Outra,
}

pub fn convencao(nome: &str) -> Convencao {
    match nome.strip_prefix("chr") {
        Some(r) if parece_cromossomo(r) => Convencao::ComChr,
        _ if parece_cromossomo(nome) => Convencao::SemChr,
        _ => Convencao::Outra,
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Estilo {
    ComChr,
    SemChr,
}

impl Estilo {
    pub fn de_texto(s: &str) -> Option<Estilo> {
        match s.to_ascii_lowercase().as_str() {
            "chr" | "with_chr" | "ucsc" => Some(Estilo::ComChr),
            "no_chr" | "ensembl" | "ncbi" => Some(Estilo::SemChr),
            _ => None,
        }
    }
}

/// O nome na convencao pedida. `None`: nenhuma regra se aplica e o nome fica
/// como esta -- ja esta na convencao, ou e um scaffold que so um mapa
/// explicito sabe converter.
///
/// O mitocondrial e o unico irregular: `chrM` de um lado, `MT` do outro.
pub fn converte_nome(nome: &str, estilo: Estilo) -> Option<String> {
    match (convencao(nome), estilo) {
        (Convencao::SemChr, Estilo::ComChr) => {
            Some(if nome == "MT" { "chrM".to_string() } else { format!("chr{nome}") })
        }
        (Convencao::ComChr, Estilo::SemChr) => {
            let r = &nome[3..];
            Some(if r == "M" { "MT".to_string() } else { r.to_string() })
        }
        (Convencao::ComChr, Estilo::ComChr) if nome == "chrMT" => Some("chrM".to_string()),
        (Convencao::SemChr, Estilo::SemChr) if nome == "M" => Some("MT".to_string()),
        _ => None,
    }
}

/// Nome sem convencao, para comparar dois arquivos que seguem convencoes
/// diferentes: `chr1` e `1` viram `1`; `chrM`, `M` e `MT` viram `MT`.
pub fn canonico(nome: &str) -> String {
    converte_nome(nome, Estilo::SemChr).unwrap_or_else(|| nome.to_string())
}

// ---------------------------------------------------------------------
//  NORMALIZAR UM ARQUIVO
// ---------------------------------------------------------------------

#[derive(Clone, Debug, Default)]
pub struct ResumoNomes {
    pub linhas: usize,
    pub registros: usize,
    /// Registros cujo nome de sequencia mudou.
    pub alterados: usize,
    /// (nome original, nome novo), na ordem em que apareceram.
    pub convertidos: Vec<(String, String)>,
    /// Nomes que nenhuma regra soube converter e ficaram como estavam.
    pub mantidos: Vec<String>,
}

struct Tradutor {
    estilo: Estilo,
    mapa: HashMap<String, String>,
    /// Decisao ja tomada para cada nome visto.
    cache: HashMap<String, Option<String>>,
    resumo: ResumoNomes,
}

impl Tradutor {
    /// O nome novo, ou `None` se o nome fica como esta.
    fn novo_nome(&mut self, nome: &str) -> Option<&str> {
        if !self.cache.contains_key(nome) {
            let novo = match self.mapa.get(nome) {
                Some(m) => Some(m.clone()),
                None => converte_nome(nome, self.estilo),
            }
            .filter(|n| n != nome);
            match &novo {
                Some(n) => self.resumo.convertidos.push((nome.to_string(), n.clone())),
                // So vale avisar de quem NAO esta na convencao pedida.
                None => {
                    let ja_esta = matches!(
                        (convencao(nome), self.estilo),
                        (Convencao::ComChr, Estilo::ComChr) | (Convencao::SemChr, Estilo::SemChr)
                    );
                    if !ja_esta && nome != "*" && nome != "=" {
                        self.resumo.mantidos.push(nome.to_string());
                    }
                }
            }
            self.cache.insert(nome.to_string(), novo);
        }
        self.cache[nome].as_deref()
    }
}

fn le_mapa(caminho: &str) -> Result<HashMap<String, String>, String> {
    let txt = std::fs::read_to_string(caminho).map_err(|e| format!("Error reading the name mapping {caminho}: {e}"))?;
    let mut mapa = HashMap::new();
    for (i, l) in txt.lines().enumerate() {
        if l.trim().is_empty() || l.starts_with('#') {
            continue;
        }
        let mut c = l.split(['\t', ' ']).filter(|x| !x.is_empty());
        match (c.next(), c.next()) {
            (Some(de), Some(para)) => {
                mapa.insert(de.to_string(), para.to_string());
            }
            // Mapas publicos trazem linhas com a segunda coluna vazia para os
            // nomes sem equivalente: nao e erro, e so "nao converte".
            (Some(_), None) => {}
            _ => return Err(format!("{caminho}, line {}: expected two columns (original name, new name)", i + 1)),
        }
    }
    Ok(mapa)
}

/// Reescreve `entrada` em `saida` com os nomes de sequencia na convencao
/// pedida. Em fluxo: serve para arquivo de qualquer tamanho.
///
/// Troca o nome nos registros (VCF, BED e GFF: coluna 1; SAM: RNAME e RNEXT) e
/// tambem onde o cabecalho o declara (`##contig` no VCF, `@SQ` no SAM,
/// `##sequence-region` no GFF) -- senao o arquivo sairia inconsistente consigo
/// mesmo.
///
/// `mapa` e um arquivo de duas colunas (nome original, nome novo) que tem
/// precedencia sobre a regra do prefixo. E o que converte scaffolds e acessos
/// RefSeq, cujo nome nao se deduz.
pub fn normalizar_cromossomos(
    entrada: &str,
    saida: &str,
    formato: Formato,
    estilo: Estilo,
    mapa: Option<&str>,
) -> Result<ResumoNomes, String> {
    let mut t = Tradutor {
        estilo,
        mapa: match mapa {
            Some(m) => le_mapa(m)?,
            None => HashMap::new(),
        },
        cache: HashMap::new(),
        resumo: ResumoNomes::default(),
    };
    let mut leitor = abre_leitura(entrada)?;
    let mut w = Escritor::cria_para(saida, formato)?;
    let mut linha = String::new();
    let mut nova = String::new();

    loop {
        linha.clear();
        if leitor.read_line(&mut linha).map_err(|e| format!("Error reading {entrada}: {e}"))? == 0 {
            break;
        }
        let texto = linha.trim_end_matches(['\n', '\r']);
        t.resumo.linhas += 1;
        nova.clear();

        let e_cabecalho = match formato {
            Formato::Vcf | Formato::Gff => texto.starts_with('#'),
            Formato::Bed => texto.starts_with('#') || texto.starts_with("track") || texto.starts_with("browser"),
            Formato::Sam => texto.starts_with('@'),
        };
        if e_cabecalho {
            match formato {
                Formato::Vcf if texto.starts_with("##contig=<") => {
                    if let Some(p) = texto.find("ID=") {
                        let ini = p + 3;
                        let fim = ini + texto[ini..].find([',', '>']).unwrap_or(texto.len() - ini);
                        if let Some(n) = t.novo_nome(&texto[ini..fim]) {
                            nova.push_str(&texto[..ini]);
                            nova.push_str(n);
                            nova.push_str(&texto[fim..]);
                        }
                    }
                }
                Formato::Sam if texto.starts_with("@SQ") => {
                    let mut mudou = false;
                    for (i, campo) in texto.split('\t').enumerate() {
                        if i > 0 {
                            nova.push('\t');
                        }
                        match campo.strip_prefix("SN:").and_then(|n| t.novo_nome(n)) {
                            Some(n) => {
                                nova.push_str("SN:");
                                nova.push_str(n);
                                mudou = true;
                            }
                            None => nova.push_str(campo),
                        }
                    }
                    if !mudou {
                        nova.clear();
                    }
                }
                Formato::Gff if texto.starts_with("##sequence-region") => {
                    let mut partes = texto.split_whitespace();
                    partes.next();
                    if let Some(nome) = partes.next() {
                        if let Some(n) = t.novo_nome(nome) {
                            nova.push_str("##sequence-region ");
                            nova.push_str(n);
                            for p in partes {
                                nova.push(' ');
                                nova.push_str(p);
                            }
                        }
                    }
                }
                _ => {}
            }
            w.linha(if nova.is_empty() { texto } else { &nova })?;
            continue;
        }
        if texto.is_empty() {
            w.linha(texto)?;
            continue;
        }

        t.resumo.registros += 1;
        match formato {
            Formato::Sam => {
                // RNAME e a coluna 3; RNEXT, a 7 (`=` e `*` nao sao nomes).
                let mut mudou = false;
                for (i, campo) in texto.split('\t').enumerate() {
                    if i > 0 {
                        nova.push('\t');
                    }
                    match if i == 2 || i == 6 { t.novo_nome(campo) } else { None } {
                        Some(n) => {
                            nova.push_str(n);
                            mudou = true;
                        }
                        None => nova.push_str(campo),
                    }
                }
                if mudou {
                    t.resumo.alterados += 1;
                } else {
                    nova.clear();
                }
            }
            _ => {
                let (nome, resto) = match texto.split_once('\t') {
                    Some((n, r)) => (n, Some(r)),
                    None => (texto, None),
                };
                if let Some(n) = t.novo_nome(nome) {
                    nova.push_str(n);
                    if let Some(r) = resto {
                        nova.push('\t');
                        nova.push_str(r);
                    }
                    t.resumo.alterados += 1;
                }
            }
        }
        w.linha(if nova.is_empty() { texto } else { &nova })?;
    }
    w.termina()?;
    Ok(t.resumo)
}

// ---------------------------------------------------------------------
//  GENOMA DE REFERENCIA
// ---------------------------------------------------------------------

/// Tamanho do cromossomo 1 humano em cada montagem. E a assinatura mais
/// simples de uma referencia: nenhuma outra sequencia tem esses tamanhos.
const MONTAGENS_HUMANAS: [(u64, &str); 2] = [(249_250_621, "GRCh37 (hg19)"), (248_956_422, "GRCh38 (hg38)")];

#[derive(Clone, Debug, Default)]
pub struct Referencia {
    /// (nome, tamanho) de cada sequencia que o arquivo declara.
    pub sequencias: Vec<(String, u64)>,
    /// O que o proprio arquivo diz sobre a referencia (`##reference=`,
    /// `@SQ AS:`, `#!genome-build`), se disser.
    pub declarada: Option<String>,
    /// A montagem humana reconhecida, pelos tamanhos ou pelo texto declarado.
    pub montagem: Option<&'static str>,
}

fn atributo<'a>(texto: &'a str, chave: &str) -> Option<&'a str> {
    let ini = texto.find(chave)? + chave.len();
    let resto = &texto[ini..];
    Some(&resto[..resto.find([',', '>', '\t']).unwrap_or(resto.len())])
}

/// Le o que um arquivo declara sobre o genoma de referencia.
///
/// Entende o cabecalho de VCF (`##contig`, `##reference`), de SAM (`@SQ`), as
/// diretivas de GFF3 (`##sequence-region`, `#!genome-build`) e arquivos de
/// tamanhos (`.fai`, `chrom.sizes`: nome e tamanho por linha). Em VCF e SAM le
/// so o cabecalho; em GFF percorre o arquivo, porque as diretivas vem
/// espalhadas.
pub fn referencia_do_arquivo(caminho: &str) -> Result<Referencia, String> {
    let mut leitor = abre_leitura(caminho)?;
    let mut r = Referencia::default();
    let mut linha = String::new();
    let mut e_gff = false;
    let mut e_tamanhos = false;
    let mut n = 0usize;

    loop {
        linha.clear();
        if leitor.read_line(&mut linha).map_err(|e| format!("Error reading {caminho}: {e}"))? == 0 {
            break;
        }
        let t = linha.trim_end_matches(['\n', '\r']);
        if t.is_empty() {
            continue;
        }
        n += 1;
        if let Some(c) = t.strip_prefix("##contig=<") {
            if let (Some(id), Some(tam)) = (atributo(c, "ID="), atributo(c, "length=").and_then(|x| x.parse().ok())) {
                r.sequencias.push((id.to_string(), tam));
            }
        } else if let Some(v) = t.strip_prefix("##reference=").or_else(|| t.strip_prefix("##assembly=")) {
            r.declarada.get_or_insert_with(|| v.to_string());
        } else if t.starts_with("@SQ") {
            let campo = |p: &str| t.split('\t').find_map(|c| c.strip_prefix(p));
            if let (Some(sn), Some(ln)) = (campo("SN:"), campo("LN:").and_then(|x| x.parse().ok())) {
                r.sequencias.push((sn.to_string(), ln));
            }
            if let Some(a) = campo("AS:") {
                r.declarada.get_or_insert_with(|| a.to_string());
            }
        } else if let Some(resto) = t.strip_prefix("##sequence-region") {
            let mut p = resto.split_whitespace();
            if let (Some(nome), Some(_), Some(fim)) = (p.next(), p.next(), p.next().and_then(|x| x.parse().ok())) {
                r.sequencias.push((nome.to_string(), fim));
            }
        } else if let Some(v) = t.strip_prefix("#!genome-build ") {
            r.declarada.get_or_insert_with(|| v.to_string());
        } else if t.starts_with("##gff-version") {
            e_gff = true;
        } else if t.starts_with('#') || t.starts_with('@') {
            // outra linha de cabecalho
        } else {
            // Primeira linha que nao e cabecalho.
            if n == 1 || e_tamanhos {
                let mut c = t.split('\t');
                if let (Some(nome), Some(tam)) = (c.next(), c.next().and_then(|x| x.parse::<u64>().ok())) {
                    // `.fai` e `chrom.sizes`: nome e tamanho, sem cabecalho.
                    e_tamanhos = true;
                    r.sequencias.push((nome.to_string(), tam));
                    continue;
                }
            }
            if !e_gff {
                break;
            }
        }
    }

    r.montagem = MONTAGENS_HUMANAS
        .iter()
        .find(|(tam, _)| r.sequencias.iter().any(|(_, t)| t == tam))
        .map(|(_, nome)| *nome)
        .or_else(|| {
            let d = r.declarada.as_deref()?.to_ascii_lowercase();
            if d.contains("grch37") || d.contains("hg19") || d.contains("b37") {
                Some(MONTAGENS_HUMANAS[0].1)
            } else if d.contains("grch38") || d.contains("hg38") {
                Some(MONTAGENS_HUMANAS[1].1)
            } else {
                None
            }
        });
    Ok(r)
}

#[derive(Clone, Debug, Default)]
pub struct Comparacao {
    pub a: Referencia,
    pub b: Referencia,
    /// Sequencias presentes nos dois com o mesmo tamanho.
    pub iguais: usize,
    /// (nome canonico, tamanho em A, tamanho em B) das que diferem.
    pub diferentes: Vec<(String, u64, u64)>,
    pub so_em_a: usize,
    pub so_em_b: usize,
    /// Os dois arquivos usam convencoes de nome diferentes (`chr1` x `1`)?
    pub convencoes_diferentes: bool,
    /// `Some(true)`: mesma referencia. `Some(false)`: referencias diferentes.
    /// `None`: os arquivos nao dizem o bastante para decidir.
    pub compativeis: Option<bool>,
    /// A conclusao, em uma frase.
    pub veredito: String,
}

fn convencao_dominante(r: &Referencia) -> Convencao {
    let conta = |c: Convencao| r.sequencias.iter().filter(|(n, _)| convencao(n) == c).count();
    let (com, sem) = (conta(Convencao::ComChr), conta(Convencao::SemChr));
    if com > sem {
        Convencao::ComChr
    } else if sem > com {
        Convencao::SemChr
    } else {
        Convencao::Outra
    }
}

/// Dois arquivos estao sobre o mesmo genoma de referencia?
///
/// Compara os TAMANHOS das sequencias que cada um declara, casando os nomes
/// sem a convencao (`chr1` e `1` sao a mesma sequencia). Tamanho e a
/// evidencia: `chr1` com 249.250.621 bases e GRCh37, com 248.956.422 e GRCh38,
/// e nenhum cabecalho precisa dizer isso.
pub fn comparar_referencias(caminho_a: &str, caminho_b: &str) -> Result<Comparacao, String> {
    let (a, b) = (referencia_do_arquivo(caminho_a)?, referencia_do_arquivo(caminho_b)?);
    let mut c = Comparacao::default();
    let em_b: HashMap<String, u64> = b.sequencias.iter().map(|(n, t)| (canonico(n), *t)).collect();
    let mut vistos = 0usize;
    for (nome, tam_a) in &a.sequencias {
        let k = canonico(nome);
        match em_b.get(&k) {
            Some(tam_b) if tam_b == tam_a => {
                c.iguais += 1;
                vistos += 1;
            }
            Some(tam_b) => {
                c.diferentes.push((k, *tam_a, *tam_b));
                vistos += 1;
            }
            None => c.so_em_a += 1,
        }
    }
    c.so_em_b = em_b.len().saturating_sub(vistos);
    let (ca, cb) = (convencao_dominante(&a), convencao_dominante(&b));
    c.convencoes_diferentes = ca != Convencao::Outra && cb != Convencao::Outra && ca != cb;

    let nome = |r: &Referencia| r.montagem.map(str::to_string).or_else(|| r.declarada.clone());
    let aviso_nomes = if c.convencoes_diferentes {
        " The names follow different conventions (`chr1` in one, `1` in the other): run `normalize_chromosomes` before combining the two."
    } else {
        ""
    };

    (c.compativeis, c.veredito) = if let Some((seq, ta, tb)) = c.diferentes.first() {
        let montagens = match (nome(&a), nome(&b)) {
            (Some(x), Some(y)) if x != y => format!(" A is {x}; B is {y}."),
            _ => String::new(),
        };
        (
            Some(false),
            format!(
                "DIFFERENT reference genomes: sequence `{seq}` has {ta} bases in A and {tb} in B ({} sequences differ, {} match).{montagens} Coordinates from one file are not valid in the other.",
                c.diferentes.len(),
                c.iguais
            ),
        )
    } else if c.iguais > 0 {
        (
            Some(true),
            format!(
                "Same reference: {} sequences with the same length in both files{}.{aviso_nomes}",
                c.iguais,
                nome(&a).map(|m| format!(" ({m})")).unwrap_or_default()
            ),
        )
    } else {
        match (a.montagem, b.montagem) {
            (Some(x), Some(y)) if x != y => (
                Some(false),
                format!("DIFFERENT reference genomes according to what the files declare: A is {x}; B is {y}. Coordinates from one file are not valid in the other."),
            ),
            (Some(x), Some(_)) => (Some(true), format!("Same reference according to what the files declare: {x}.{aviso_nomes}")),
            _ if a.sequencias.is_empty() || b.sequencias.is_empty() => (
                None,
                "Could not compare: at least one of the files does not declare its sequence lengths (missing `##contig` with `length=` in a VCF, `@SQ` in a SAM or `##sequence-region` in a GFF).".to_string(),
            ),
            _ => (
                None,
                "The two files declare no sequence in common, so the references cannot be compared.".to_string(),
            ),
        }
    };
    c.a = a;
    c.b = b;
    Ok(c)
}

// ---------------------------------------------------------------------
//  SAIDA EM JSON (para a camada Python)
// ---------------------------------------------------------------------

/// `s` como literal de texto JSON, com aspas.
pub fn json_texto(s: &str) -> String {
    let mut o = String::with_capacity(s.len() + 2);
    o.push('"');
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            '\r' => o.push_str("\\r"),
            '\t' => o.push_str("\\t"),
            c if (c as u32) < 0x20 => o.push_str(&format!("\\u{:04x}", c as u32)),
            c => o.push(c),
        }
    }
    o.push('"');
    o
}

fn json_opcional(s: Option<&str>) -> String {
    s.map(json_texto).unwrap_or_else(|| "null".to_string())
}

impl ResumoNomes {
    pub fn json(&self) -> String {
        let pares: Vec<String> =
            self.convertidos.iter().map(|(a, b)| format!("[{},{}]", json_texto(a), json_texto(b))).collect();
        let mantidos: Vec<String> = self.mantidos.iter().map(|n| json_texto(n)).collect();
        format!(
            "{{\"lines\":{},\"records\":{},\"changed\":{},\"converted\":[{}],\"kept\":[{}]}}",
            self.linhas,
            self.registros,
            self.alterados,
            pares.join(","),
            mantidos.join(",")
        )
    }
}

impl Referencia {
    pub fn json(&self) -> String {
        let seqs: Vec<String> = self.sequencias.iter().map(|(n, t)| format!("[{},{t}]", json_texto(n))).collect();
        format!(
            "{{\"assembly\":{},\"declared\":{},\"sequences\":[{}]}}",
            json_opcional(self.montagem),
            json_opcional(self.declarada.as_deref()),
            seqs.join(",")
        )
    }
}

impl Comparacao {
    pub fn json(&self) -> String {
        let dif: Vec<String> =
            self.diferentes.iter().map(|(n, a, b)| format!("[{},{a},{b}]", json_texto(n))).collect();
        format!(
            "{{\"compatible\":{},\"verdict\":{},\"matching\":{},\"different\":[{}],\"only_in_a\":{},\"only_in_b\":{},\"different_naming\":{},\"assembly_a\":{},\"assembly_b\":{},\"declared_a\":{},\"declared_b\":{}}}",
            match self.compativeis {
                Some(v) => v.to_string(),
                None => "null".to_string(),
            },
            json_texto(&self.veredito),
            self.iguais,
            dif.join(","),
            self.so_em_a,
            self.so_em_b,
            self.convencoes_diferentes,
            json_opcional(self.a.montagem),
            json_opcional(self.b.montagem),
            json_opcional(self.a.declarada.as_deref()),
            json_opcional(self.b.declarada.as_deref())
        )
    }
}
