// Preparo do VCF de entrada do PharmCAT e a chamada dele.
//
// O PharmCAT (PharmGKB/CPIC, licenca MPL-2.0, https://pharmcat.org) chama os
// alelos-estrela, os fenotipos e as recomendacoes a partir de um VCF com o
// genotipo da pessoa em cada posicao que ele examina (o arquivo
// `pharmcat_positions_<versao>.vcf` de cada versao). O MMGEN nao chama alelos:
// so prepara essa entrada e le o que o PharmCAT devolve.
//
// Por que preparar: um VCF de chamada de variantes so traz as posicoes onde a
// pessoa difere da referencia, e o PharmCAT trata como "nao lida" toda posicao
// ausente. Aqui cada posicao recebe:
//  - o genotipo da pessoa, quando ela tem um registro ali, com os alelos
//    traduzidos para a lista de alelos da posicao (um alelo que a lista nao tem
//    e acrescentado: o PharmCAT avisa que nao o conhece);
//  - 0/0 quando nao ha registro e a posicao esta inteira dentro da regiao em
//    que as chamadas sao confiaveis (o BED de um GIAB, por exemplo) ou, sem
//    BED, quando `ausente_e_ref`;
//  - ./. no resto, e sempre numa sequencia que nao aparece no arquivo da
//    pessoa (um VCF so do chr22 nao vira "referencia" no genoma inteiro).

use crate::externo;
use crate::indice;
use crate::nomes;
use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::path::Path;
use std::time::Instant;

/// Uma posicao do arquivo de posicoes do PharmCAT.
pub struct Posicao {
    pub cromo: String,
    pub pos: u64,
    pub id: String,
    pub refe: String,
    pub alts: Vec<String>,
    pub info: String,
}

pub struct Posicoes {
    pub cabecalho: Vec<String>,
    pub lista: Vec<Posicao>,
}

/// Le o `pharmcat_positions_<versao>.vcf` (texto, .gz ou BGZF).
pub fn le_posicoes(caminho: &str) -> Result<Posicoes, String> {
    let leitor = externo::abre_leitura(caminho)?;
    let mut cabecalho = Vec::new();
    let mut lista = Vec::new();
    for (i, linha) in leitor.lines().enumerate() {
        let linha = linha.map_err(|e| format!("{caminho}: {e}"))?;
        let l = linha.trim_end_matches(['\r', '\n']);
        if l.starts_with("##") {
            cabecalho.push(l.to_string());
            continue;
        }
        if l.is_empty() || l.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() < 8 {
            return Err(format!(
                "{caminho}, line {}: expected the 8 columns of a VCF, found {}. This must be the \
                 pharmcat_positions_<version>.vcf file of a PharmCAT release.",
                i + 1,
                f.len()
            ));
        }
        let pos = f[1]
            .parse::<u64>()
            .map_err(|_| format!("{caminho}, line {}: POS `{}` is not a number", i + 1, f[1]))?;
        lista.push(Posicao {
            cromo: f[0].to_string(),
            pos,
            id: f[2].to_string(),
            refe: f[3].to_ascii_uppercase(),
            alts: if f[4] == "." { Vec::new() } else { f[4].split(',').map(|s| s.to_ascii_uppercase()).collect() },
            info: f[7].to_string(),
        });
    }
    if lista.is_empty() {
        return Err(format!(
            "{caminho} has no positions. It must be the pharmcat_positions_<version>.vcf file of a PharmCAT release."
        ));
    }
    Ok(Posicoes { cabecalho, lista })
}

/// Regioes em que as chamadas sao confiaveis, por cromossomo (nome canonico),
/// em base 0, meio-abertas, ordenadas e fundidas.
struct Chamavel {
    por_cromo: HashMap<String, Vec<(u64, u64)>>,
}

impl Chamavel {
    fn le(caminho: &str, cromos: &HashSet<String>) -> Result<Chamavel, String> {
        let leitor = externo::abre_leitura(caminho)?;
        let mut por_cromo: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        for (i, linha) in leitor.lines().enumerate() {
            let linha = linha.map_err(|e| format!("{caminho}: {e}"))?;
            let l = linha.trim_end_matches(['\r', '\n']);
            if l.is_empty() || l.starts_with('#') || l.starts_with("track") || l.starts_with("browser") {
                continue;
            }
            let mut f = l.split('\t');
            let (Some(c), Some(a), Some(b)) = (f.next(), f.next(), f.next()) else {
                return Err(format!("{caminho}, line {}: a BED line needs chrom, start and end separated by tabs", i + 1));
            };
            let canon = nomes::canonico(c);
            if !cromos.contains(&canon) {
                continue;
            }
            let (Ok(a), Ok(b)) = (a.parse::<u64>(), b.parse::<u64>()) else {
                return Err(format!("{caminho}, line {}: start and end must be numbers, found `{a}` and `{b}`", i + 1));
            };
            por_cromo.entry(canon).or_default().push((a, b));
        }
        for v in por_cromo.values_mut() {
            v.sort_unstable();
            let mut fundido: Vec<(u64, u64)> = Vec::with_capacity(v.len());
            for &(a, b) in v.iter() {
                match fundido.last_mut() {
                    Some(u) if a <= u.1 => u.1 = u.1.max(b),
                    _ => fundido.push((a, b)),
                }
            }
            *v = fundido;
        }
        Ok(Chamavel { por_cromo })
    }

    /// A faixa [inicio0, fim0) esta inteira dentro de uma regiao?
    fn contem(&self, cromo: &str, inicio0: u64, fim0: u64) -> bool {
        let Some(v) = self.por_cromo.get(cromo) else { return false };
        let i = v.partition_point(|&(a, _)| a <= inicio0);
        i > 0 && v[i - 1].1 >= fim0
    }
}

/// Um registro da pessoa numa posicao de interesse.
struct Registro {
    refe: String,
    alts: Vec<String>,
    gt: String,
}

pub struct Opcoes {
    pub amostra: Option<String>,
    pub chamavel: Option<String>,
    pub ausente_e_ref: bool,
}

/// O que o preparo fez.
#[derive(Default)]
pub struct Preparo {
    pub entrada: String,
    pub saida: String,
    pub amostra: String,
    pub indexado: bool,
    pub posicoes: usize,
    /// A pessoa tem ali um alelo diferente da referencia.
    pub com_variante: usize,
    /// A pessoa tem ali um registro com genotipo de referencia.
    pub registro_referencia: usize,
    /// Sem registro, dentro da regiao confiavel (ou sem BED, com ausente_e_ref): 0/0.
    pub referencia_preenchida: usize,
    /// Sem registro e fora da regiao confiavel: ./.
    pub faltando_fora_da_regiao: usize,
    /// Numa sequencia que o arquivo da pessoa nao tem: ./.
    pub faltando_sem_sequencia: usize,
    /// O registro da pessoa tem genotipo ausente (./.).
    pub faltando_no_registro: usize,
    /// O REF da pessoa nao bate com o da posicao: ./.
    pub ref_divergente: usize,
    /// Alelos da pessoa que a lista da posicao nao tem (acrescentados).
    pub alelos_acrescentados: usize,
    pub sequencias_ausentes: Vec<String>,
    pub notas: Vec<String>,
}

impl Preparo {
    pub fn json(&self) -> String {
        let lista = |v: &[String]| format!("[{}]", v.iter().map(|s| nomes::json_texto(s)).collect::<Vec<_>>().join(","));
        format!(
            concat!(
                "{{\"input\":{},\"output\":{},\"sample\":{},\"indexed\":{},\"positions\":{},",
                "\"with_variant\":{},\"reference_record\":{},\"reference_filled\":{},",
                "\"missing_outside_callable\":{},\"missing_no_sequence\":{},\"missing_in_record\":{},",
                "\"ref_mismatch\":{},\"alleles_added\":{},\"absent_sequences\":{},\"notes\":{}}}"
            ),
            nomes::json_texto(&self.entrada),
            nomes::json_texto(&self.saida),
            nomes::json_texto(&self.amostra),
            self.indexado,
            self.posicoes,
            self.com_variante,
            self.registro_referencia,
            self.referencia_preenchida,
            self.faltando_fora_da_regiao,
            self.faltando_sem_sequencia,
            self.faltando_no_registro,
            self.ref_divergente,
            self.alelos_acrescentados,
            lista(&self.sequencias_ausentes),
            lista(&self.notas),
        )
    }

    pub fn lidas(&self) -> usize {
        self.com_variante + self.registro_referencia + self.referencia_preenchida
    }
}

/// O campo GT da amostra (`./.` se nao houver).
fn gt_da_amostra(f: &[&str], col: usize) -> String {
    let Some(formato) = f.get(8) else { return "./.".into() };
    let Some(i) = formato.split(':').position(|k| k == "GT") else { return "./.".into() };
    f.get(col)
        .and_then(|s| s.split(':').nth(i))
        .map(|s| s.to_string())
        .unwrap_or_else(|| "./.".into())
}

/// Coluna da amostra na linha #CHROM.
fn coluna_da_amostra(caminho: &str, linha: &str, amostra: Option<&str>) -> Result<(usize, String), String> {
    let f: Vec<&str> = linha.trim_end_matches(['\r', '\n']).split('\t').collect();
    if f.len() < 10 {
        return Err(format!(
            "{caminho} has no sample column (a sites-only VCF): PharmCAT needs the genotypes of one person."
        ));
    }
    match amostra {
        None => Ok((9, f[9].to_string())),
        Some(a) => match f[9..].iter().position(|s| *s == a) {
            Some(i) => Ok((9 + i, a.to_string())),
            None => Err(format!(
                "Sample `{a}` is not in {caminho}. Its samples are: {}",
                f[9..].iter().take(20).cloned().collect::<Vec<_>>().join(", ")
                    + if f.len() > 29 { ", ..." } else { "" }
            )),
        },
    }
}

/// Traduz o genotipo de um registro da pessoa para a lista de alelos da
/// posicao. `None` quando o REF dos dois nao se reconcilia.
fn traduz(p: &Posicao, r: &Registro, alts_saida: &mut Vec<String>, acrescentados: &mut usize) -> Option<String> {
    // Os alelos da pessoa escritos sobre o REF da posicao.
    let alelos: Vec<String> = if r.refe == p.refe {
        r.alts.clone()
    } else if r.refe.starts_with(&p.refe) {
        // REF da pessoa mais longo: tira o mesmo sufixo de todos os alelos.
        let sufixo = &r.refe[p.refe.len()..];
        let mut v = Vec::new();
        for a in &r.alts {
            if a == "*" {
                v.push(a.clone());
            } else if a.len() > sufixo.len() && a.ends_with(sufixo) {
                v.push(a[..a.len() - sufixo.len()].to_string());
            } else {
                return None;
            }
        }
        v
    } else if p.refe.starts_with(&r.refe) {
        // REF da posicao mais longo: estende os alelos da pessoa com o resto.
        let sufixo = &p.refe[r.refe.len()..];
        r.alts.iter().map(|a| if a == "*" { a.clone() } else { format!("{a}{sufixo}") }).collect()
    } else {
        return None;
    };
    let mut partes = Vec::new();
    let mut seps = Vec::new();
    let mut atual = String::new();
    for ch in r.gt.chars() {
        if ch == '/' || ch == '|' {
            partes.push(std::mem::take(&mut atual));
            seps.push(ch);
        } else {
            atual.push(ch);
        }
    }
    partes.push(atual);
    let mut saida = String::new();
    for (i, a) in partes.iter().enumerate() {
        let traduzido = if a == "." || a.is_empty() {
            ".".to_string()
        } else if a == "0" {
            "0".to_string()
        } else {
            let k: usize = a.parse().ok()?;
            let alelo = alelos.get(k.checked_sub(1)?)?;
            if *alelo == p.refe {
                "0".to_string()
            } else {
                let j = match alts_saida.iter().position(|x| x == alelo) {
                    Some(j) => j,
                    None => {
                        alts_saida.push(alelo.clone());
                        *acrescentados += 1;
                        alts_saida.len() - 1
                    }
                };
                (j + 1).to_string()
            }
        };
        saida.push_str(&traduzido);
        if let Some(s) = seps.get(i) {
            saida.push(*s);
        }
    }
    Some(saida)
}

/// Escreve o VCF de entrada do PharmCAT para uma pessoa.
pub fn preparar(vcf: &str, posicoes: &str, saida: &str, op: &Opcoes) -> Result<Preparo, String> {
    if !Path::new(vcf).exists() {
        return Err(format!("{vcf} does not exist"));
    }
    let pos = le_posicoes(posicoes)?;
    let cromos: HashSet<String> = pos.lista.iter().map(|p| nomes::canonico(&p.cromo)).collect();
    let mut chave: HashMap<(String, u64), Vec<usize>> = HashMap::new();
    for (i, p) in pos.lista.iter().enumerate() {
        chave.entry((nomes::canonico(&p.cromo), p.pos)).or_default().push(i);
    }
    let chamavel = match &op.chamavel {
        Some(b) => Some(Chamavel::le(b, &cromos)?),
        None => None,
    };

    let mut registros: HashMap<(String, u64), Vec<Registro>> = HashMap::new();
    let mut col: Option<(usize, String)> = None;
    let mut sequencias: HashSet<String> = HashSet::new();
    let tbi = format!("{vcf}.tbi");
    let indexado = Path::new(&tbi).exists();

    // Uma linha de registro: guarda se cair numa posicao de interesse.
    let recebe_registro = |linha: &str, col: &(usize, String), registros: &mut HashMap<(String, u64), Vec<Registro>>,
                               sequencias: &mut HashSet<String>|
     -> Result<(), String> {
        let f: Vec<&str> = linha.trim_end_matches(['\r', '\n']).split('\t').collect();
        if f.len() < 8 {
            return Ok(());
        }
        let canon = nomes::canonico(f[0]);
        sequencias.insert(canon.clone());
        let Ok(p) = f[1].parse::<u64>() else { return Ok(()) };
        let k = (canon, p);
        if chave.contains_key(&k) {
            registros.entry(k).or_default().push(Registro {
                refe: f[3].to_ascii_uppercase(),
                alts: if f[4] == "." { Vec::new() } else { f[4].split(',').map(|s| s.to_ascii_uppercase()).collect() },
                gt: gt_da_amostra(&f, col.0),
            });
        }
        Ok(())
    };

    if indexado {
        let (_, _, nomes_arquivo) = indice::le_tbi(&tbi)?;
        let por_canon: HashMap<String, String> =
            nomes_arquivo.iter().map(|n| (nomes::canonico(n), n.clone())).collect();
        for n in &nomes_arquivo {
            sequencias.insert(nomes::canonico(n));
        }
        let regioes: Vec<String> = pos
            .lista
            .iter()
            .filter_map(|p| {
                por_canon
                    .get(&nomes::canonico(&p.cromo))
                    .map(|n| format!("{n}:{}-{}", p.pos, p.pos + p.refe.len().max(1) as u64 - 1))
            })
            .collect();
        // Sem nenhuma regiao, consultar so le o cabecalho: uma regiao qualquer
        // de uma sequencia do arquivo basta para isso.
        let regioes = if regioes.is_empty() {
            nomes_arquivo.first().map(|n| vec![format!("{n}:1-1")]).unwrap_or_default()
        } else {
            regioes
        };
        if regioes.is_empty() {
            return Err(format!("{vcf}: its index lists no sequence"));
        }
        let mut erro_amostra: Option<String> = None;
        indice::consultar(vcf, &regioes, &mut |cabecalho, linha| {
            if cabecalho {
                if linha.starts_with("#CHROM") {
                    match coluna_da_amostra(vcf, linha, op.amostra.as_deref()) {
                        Ok(c) => col = Some(c),
                        Err(e) => erro_amostra = Some(e),
                    }
                }
                return Ok(());
            }
            match &col {
                Some(c) => {
                    let c = c.clone();
                    recebe_registro(linha, &c, &mut registros, &mut sequencias)
                }
                None => Err(erro_amostra.clone().unwrap_or_else(|| format!("{vcf} has no #CHROM header line"))),
            }
        })?;
        if let Some(e) = erro_amostra {
            return Err(e);
        }
        if col.is_none() {
            return Err(format!("{vcf} has no #CHROM header line"));
        }
    } else {
        let leitor = externo::abre_leitura(vcf)?;
        for linha in leitor.lines() {
            let linha = linha.map_err(|e| format!("{vcf}: {e}"))?;
            if linha.starts_with("##") {
                continue;
            }
            if linha.starts_with("#CHROM") {
                col = Some(coluna_da_amostra(vcf, &linha, op.amostra.as_deref())?);
                continue;
            }
            if linha.is_empty() || linha.starts_with('#') {
                continue;
            }
            let Some(c) = &col else {
                return Err(format!("{vcf}: records before the #CHROM header line"));
            };
            let c = c.clone();
            recebe_registro(&linha, &c, &mut registros, &mut sequencias)?;
        }
        if col.is_none() {
            return Err(format!("{vcf} has no #CHROM header line: is it a VCF?"));
        }
    }
    let (_, nome_amostra) = col.unwrap();

    let mut r = Preparo {
        entrada: vcf.to_string(),
        saida: saida.to_string(),
        amostra: nome_amostra.clone(),
        indexado,
        posicoes: pos.lista.len(),
        ..Default::default()
    };
    let mut ausentes: Vec<String> = Vec::new();
    let mut linhas = Vec::with_capacity(pos.lista.len());
    for p in &pos.lista {
        let canon = nomes::canonico(&p.cromo);
        let mut alts_saida = p.alts.clone();
        let candidatos = registros.get(&(canon.clone(), p.pos));
        let mut gt: Option<String> = None;
        let mut divergente = false;
        if let Some(lista) = candidatos {
            for reg in lista {
                match traduz(p, reg, &mut alts_saida, &mut r.alelos_acrescentados) {
                    Some(g) => {
                        gt = Some(g);
                        break;
                    }
                    None => divergente = true,
                }
            }
        }
        let gt = match gt {
            Some(g) => {
                if g.split(['/', '|']).all(|a| a == ".") {
                    r.faltando_no_registro += 1;
                } else if g.split(['/', '|']).any(|a| a != "0" && a != ".") {
                    r.com_variante += 1;
                } else {
                    r.registro_referencia += 1;
                }
                g
            }
            None if divergente => {
                r.ref_divergente += 1;
                if r.notas.len() < 20 {
                    r.notas.push(format!(
                        "{}:{} {}: the person's record has another REF; the position is left missing",
                        p.cromo, p.pos, p.id
                    ));
                }
                "./.".to_string()
            }
            None => {
                if !sequencias.contains(&canon) && chamavel.as_ref().map_or(true, |c| !c.por_cromo.contains_key(&canon)) {
                    r.faltando_sem_sequencia += 1;
                    if !ausentes.contains(&p.cromo) {
                        ausentes.push(p.cromo.clone());
                    }
                    "./.".to_string()
                } else {
                    let dentro = match &chamavel {
                        Some(c) => c.contem(&canon, p.pos - 1, p.pos - 1 + p.refe.len() as u64),
                        None => op.ausente_e_ref,
                    };
                    if dentro {
                        r.referencia_preenchida += 1;
                        "0/0".to_string()
                    } else {
                        r.faltando_fora_da_regiao += 1;
                        "./.".to_string()
                    }
                }
            }
        };
        let alt = if alts_saida.is_empty() { ".".to_string() } else { alts_saida.join(",") };
        linhas.push(format!("{}\t{}\t{}\t{}\t{}\t.\tPASS\t{}\tGT\t{}", p.cromo, p.pos, p.id, p.refe, alt, p.info, gt));
    }
    r.sequencias_ausentes = ausentes;

    let mut w = externo::Escritor::cria(saida)?;
    for l in &pos.cabecalho {
        w.linha(l)?;
    }
    w.linha(&format!(
        "##MMGEN=<Description=\"pharmcat_prepare: genotypes of {} at the PharmCAT positions; absent records are 0/0 only {}\">",
        nome_amostra,
        if op.chamavel.is_some() { "inside the callable regions" } else if op.ausente_e_ref { "on the sequences present in the file" } else { "never" }
    ))?;
    w.linha(&format!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{nome_amostra}"))?;
    for l in &linhas {
        w.linha(l)?;
    }
    w.termina()?;
    Ok(r)
}

/// Uma linha da tabela de chamadas do PharmCAT (`<base>.report.tsv`), como
/// pares (coluna, valor), com os nomes das colunas em snake_case.
pub struct Chamada {
    pub campos: Vec<(String, String)>,
}

pub struct Rodada {
    pub preparo: Preparo,
    pub versao: String,
    pub chamadas: Vec<Chamada>,
    pub arquivos: Vec<(String, String)>,
    pub avisos_de_posicao: usize,
    pub segundos_preparo: f64,
    pub segundos_pharmcat: f64,
    pub comando: Vec<String>,
}

fn snake(s: &str) -> String {
    let mut out = String::new();
    for ch in s.trim().chars() {
        if ch.is_ascii_alphanumeric() {
            out.push(ch.to_ascii_lowercase());
        } else if !out.ends_with('_') && !out.is_empty() {
            out.push('_');
        }
    }
    out.trim_end_matches('_').to_string()
}

/// Le a tabela de chamadas do PharmCAT.
pub fn le_chamadas(caminho: &str) -> Result<(String, Vec<Chamada>), String> {
    let texto = std::fs::read_to_string(caminho).map_err(|e| format!("Error reading {caminho}: {e}"))?;
    let mut linhas = texto.lines();
    let versao = linhas.next().unwrap_or("").trim().to_string();
    let Some(cab) = linhas.next() else {
        return Err(format!("{caminho} has no header line"));
    };
    let colunas: Vec<String> = cab.split('\t').map(snake).collect();
    let mut chamadas = Vec::new();
    for l in linhas {
        if l.trim().is_empty() {
            continue;
        }
        let valores: Vec<&str> = l.split('\t').collect();
        let campos = colunas
            .iter()
            .enumerate()
            .map(|(i, c)| (c.clone(), valores.get(i).map(|v| v.trim().to_string()).unwrap_or_default()))
            .collect();
        chamadas.push(Chamada { campos });
    }
    Ok((versao, chamadas))
}

/// O arquivo de posicoes que acompanha um `pharmcat-<versao>-all.jar`: o de
/// mesma versao na mesma pasta, ou o mais novo de la.
pub fn posicoes_junto_do_jar(jar: &str) -> Option<String> {
    let p = Path::new(jar);
    let pasta = p.parent().filter(|d| !d.as_os_str().is_empty()).unwrap_or(Path::new("."));
    let nome = p.file_name()?.to_string_lossy().to_string();
    let versao = nome.strip_prefix("pharmcat-").and_then(|s| s.strip_suffix("-all.jar")).map(|s| s.to_string());
    let mut achados: Vec<String> = std::fs::read_dir(pasta)
        .ok()?
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().to_string())
        .filter(|n| n.starts_with("pharmcat_positions_") && n.ends_with(".vcf") && !n.contains("uniallelic"))
        .collect();
    achados.sort();
    if let Some(v) = versao {
        if let Some(n) = achados.iter().find(|n| **n == format!("pharmcat_positions_{v}.vcf")) {
            return Some(pasta.join(n).to_string_lossy().to_string());
        }
    }
    achados.last().map(|n| pasta.join(n).to_string_lossy().to_string())
}

/// Prepara a entrada e roda o PharmCAT (`java -jar`), com os relatorios HTML,
/// JSON e a tabela de chamadas.
#[allow(clippy::too_many_arguments)]
pub fn rodar(
    vcf: &str,
    pasta: &str,
    jar: &str,
    posicoes: Option<&str>,
    op: &Opcoes,
    java: &str,
    fora: Option<&str>,
    base: Option<&str>,
) -> Result<Rodada, String> {
    if !Path::new(jar).exists() {
        return Err(format!(
            "{jar} does not exist. Download pharmcat-<version>-all.jar and pharmcat_positions_<version>.vcf from \
             https://github.com/PharmGKB/PharmCAT/releases and put them in the same folder."
        ));
    }
    let posicoes = match posicoes {
        Some(p) => p.to_string(),
        None => posicoes_junto_do_jar(jar).ok_or_else(|| {
            format!(
                "No pharmcat_positions_<version>.vcf next to {jar}. Download it from the same PharmCAT release \
                 (https://github.com/PharmGKB/PharmCAT/releases) into that folder, or pass positions=."
            )
        })?,
    };
    if let Some(f) = fora {
        if !Path::new(f).exists() {
            return Err(format!("The outside-call file {f} does not exist"));
        }
    }
    std::fs::create_dir_all(pasta).map_err(|e| format!("Error creating {pasta}: {e}"))?;
    let base = match base {
        Some(b) if !b.trim().is_empty() => b.trim().to_string(),
        _ => {
            let n = Path::new(vcf).file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_else(|| "person".into());
            let mut n = n.as_str();
            for ext in [".gz", ".bgz", ".lz4", ".vcf"] {
                n = n.strip_suffix(ext).unwrap_or(n);
            }
            n.to_string()
        }
    };
    let entrada = Path::new(pasta).join(format!("{base}.pharmcat_input.vcf")).to_string_lossy().to_string();
    let t0 = Instant::now();
    let preparo = preparar(vcf, &posicoes, &entrada, op)?;
    let segundos_preparo = t0.elapsed().as_secs_f64();

    // PharmCAT runs inside the output folder (it writes a pharmcat.log where it
    // runs), so every path it gets is absolute.
    let abs = |p: &str| std::path::absolute(Path::new(p)).map(|a| a.to_string_lossy().to_string()).unwrap_or_else(|_| p.to_string());
    let mut comando: Vec<String> = vec![
        java.to_string(),
        "-jar".into(),
        abs(jar),
        "-vcf".into(),
        abs(&entrada),
        "-o".into(),
        abs(pasta),
        "-bf".into(),
        base.clone(),
        "-reporterHtml".into(),
        "-reporterJson".into(),
        "-reporterCallsOnlyTsv".into(),
    ];
    if let Some(f) = fora {
        comando.push("-po".into());
        comando.push(abs(f));
    }
    let t1 = Instant::now();
    let saida = std::process::Command::new(&comando[0]).args(&comando[1..]).current_dir(pasta).output().map_err(|e| {
        if e.kind() == std::io::ErrorKind::NotFound {
            format!(
                "Java was not found (`{java}`). PharmCAT needs Java 17 or newer: install it (for example from \
                 https://adoptium.net) and try again, or pass java= with the full path to the java program."
            )
        } else {
            format!("Error starting PharmCAT with `{java}`: {e}")
        }
    })?;
    let segundos_pharmcat = t1.elapsed().as_secs_f64();
    if !saida.status.success() {
        let cauda = |b: &[u8]| {
            let s = String::from_utf8_lossy(b);
            let s = s.trim();
            if s.len() > 1500 { format!("...{}", &s[s.len() - 1500..]) } else { s.to_string() }
        };
        return Err(format!(
            "PharmCAT stopped with an error ({}). Its messages:\n{}\n{}",
            saida.status,
            cauda(&saida.stdout),
            cauda(&saida.stderr)
        ));
    }
    let arq = |sufixo: &str| Path::new(pasta).join(format!("{base}.{sufixo}")).to_string_lossy().to_string();
    let tsv = arq("report.tsv");
    let (versao, chamadas) = le_chamadas(&tsv)?;
    let avisos = arq("match_warnings.txt");
    let avisos_de_posicao = std::fs::read_to_string(&avisos)
        .map(|t| t.lines().filter(|l| !l.starts_with(char::is_whitespace) && !l.trim().is_empty()).count())
        .unwrap_or(0);
    let mut arquivos = vec![("input_vcf".to_string(), entrada.clone())];
    for (chave, sufixo) in [
        ("report_html", "report.html"),
        ("report_json", "report.json"),
        ("calls_tsv", "report.tsv"),
        ("match_json", "match.json"),
        ("phenotype_json", "phenotype.json"),
        ("match_warnings", "match_warnings.txt"),
    ] {
        let c = arq(sufixo);
        if Path::new(&c).exists() {
            arquivos.push((chave.to_string(), c));
        }
    }
    Ok(Rodada { preparo, versao, chamadas, arquivos, avisos_de_posicao, segundos_preparo, segundos_pharmcat, comando })
}

impl Rodada {
    pub fn json(&self) -> String {
        let chamadas: Vec<String> = self
            .chamadas
            .iter()
            .map(|c| {
                format!(
                    "{{{}}}",
                    c.campos
                        .iter()
                        .map(|(k, v)| format!("{}:{}", nomes::json_texto(k), nomes::json_texto(v)))
                        .collect::<Vec<_>>()
                        .join(",")
                )
            })
            .collect();
        let arquivos: Vec<String> = self
            .arquivos
            .iter()
            .map(|(k, v)| format!("{}:{}", nomes::json_texto(k), nomes::json_texto(v)))
            .collect();
        let comando: Vec<String> = self.comando.iter().map(|s| nomes::json_texto(s)).collect();
        format!(
            "{{\"pharmcat_version\":{},\"prepared\":{},\"calls\":[{}],\"files\":{{{}}},\"position_warnings\":{},\"seconds_prepare\":{:.3},\"seconds_pharmcat\":{:.3},\"command\":[{}]}}",
            nomes::json_texto(&self.versao),
            self.preparo.json(),
            chamadas.join(","),
            arquivos.join(","),
            self.avisos_de_posicao,
            self.segundos_preparo,
            self.segundos_pharmcat,
            comando.join(","),
        )
    }
}
