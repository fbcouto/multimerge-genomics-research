// =====================================================================
//  GENOMA DE REFERENCIA (FASTA) COM ACESSO ALEATORIO
//
//  O alinhamento a esquerda do `vcf_normalize` precisa das bases do genoma
//  em volta de cada indel. O genoma nunca vai inteiro para a memoria:
//
//   * FASTA de texto: as bases sao lidas direto do arquivo, pela conta do
//     `.fai` (o do `samtools faidx`, se estiver ao lado e nao for mais velho
//     que o FASTA; senao, uma passada pelo arquivo monta o mesmo indice em
//     memoria, sem gravar nada ao lado do arquivo do usuario);
//   * `.gz` (gzip ou BGZF), `.lz4`, ou um FASTA de texto com linhas de
//     tamanhos irregulares: uma passada copia so as bases, sem quebras de
//     linha, para um arquivo temporario, e guarda onde cada sequencia comeca.
//
//  As bases sao lidas em janelas de 64 kb, guardadas em maiusculas.
//  `chr1` x `1` (e `chrM` x `MT`): se o nome do VCF nao existe no FASTA,
//  vale o mesmo nome na outra convencao.
// =====================================================================

use crate::externo::{abre_leitura, DirTemp};
use crate::nomes::{converte_nome, Estilo};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

/// Bases por janela lida do arquivo.
const JANELA: u64 = 1 << 16;

/// Onde uma sequencia esta no arquivo, como numa linha do `.fai`.
#[derive(Clone, Copy, Debug)]
struct Entrada {
    comprimento: u64,
    offset: u64,
    bases_por_linha: u64,
    bytes_por_linha: u64,
}

impl Entrada {
    /// Byte do arquivo onde esta a base `p` (0-based).
    fn byte(&self, p: u64) -> u64 {
        self.offset + (p / self.bases_por_linha) * self.bytes_por_linha + p % self.bases_por_linha
    }
}

pub struct Genoma {
    caminho: String,
    arq: File,
    seqs: HashMap<String, Entrada>,
    /// Nome do VCF -> nome no FASTA (ou nenhum).
    nomes: HashMap<String, Option<String>>,
    /// A janela em memoria: sequencia, inicio (0-based) e bases.
    janela: (String, u64, Vec<u8>),
    /// O temporario com as bases, quando foi preciso copia-las.
    _temp: Option<DirTemp>,
}

fn comprimido(caminho: &str) -> bool {
    let c = caminho.to_ascii_lowercase();
    [".gz", ".bgz", ".lz4", ".bz2", ".xz", ".zst"].iter().any(|s| c.ends_with(s))
}

fn nome_da_sequencia(cabecalho: &[u8]) -> String {
    let t = String::from_utf8_lossy(&cabecalho[1..]);
    t.split_whitespace().next().unwrap_or("").to_string()
}

/// Le o `.fai` ao lado de um FASTA de texto, se existir e nao for mais velho
/// que o FASTA.
fn le_fai(caminho: &str) -> Result<Option<HashMap<String, Entrada>>, String> {
    let fai = format!("{caminho}.fai");
    let (Ok(mf), Ok(mi)) = (std::fs::metadata(caminho), std::fs::metadata(&fai)) else {
        return Ok(None);
    };
    if let (Ok(tf), Ok(ti)) = (mf.modified(), mi.modified()) {
        if ti < tf {
            return Ok(None);
        }
    }
    let texto = std::fs::read_to_string(&fai).map_err(|e| format!("Error reading {fai}: {e}"))?;
    let mut seqs = HashMap::new();
    for (i, l) in texto.lines().enumerate() {
        if l.trim().is_empty() {
            continue;
        }
        let c: Vec<&str> = l.split('\t').collect();
        let num = |j: usize| c.get(j).and_then(|v| v.trim().parse::<u64>().ok());
        let (Some(comprimento), Some(offset), Some(bl), Some(bw)) = (num(1), num(2), num(3), num(4)) else {
            return Err(format!(
                "{fai}, line {}: not a .fai line (name, length, offset, bases per line, bytes per line). Rebuild it with `samtools faidx`, or delete it.",
                i + 1
            ));
        };
        let e = if bl == 0 || bw == 0 {
            Entrada { comprimento: 0, offset, bases_por_linha: 1, bytes_por_linha: 1 }
        } else {
            Entrada { comprimento, offset, bases_por_linha: bl, bytes_por_linha: bw }
        };
        seqs.insert(c[0].to_string(), e);
    }
    Ok(Some(seqs))
}

/// Monta em memoria o mesmo indice do `samtools faidx`, numa passada por um
/// FASTA de texto. `None` quando as linhas de alguma sequencia nao tem todas
/// o mesmo tamanho (fora a ultima): ai o indice nao serve.
fn indexa_texto(caminho: &str) -> Result<Option<HashMap<String, Entrada>>, String> {
    let arq = File::open(caminho).map_err(|e| format!("Error opening {caminho}: {e}"))?;
    let mut r = BufReader::with_capacity(16 << 20, arq);
    let mut seqs: HashMap<String, Entrada> = HashMap::new();
    let mut linha: Vec<u8> = Vec::new();
    let mut byte: u64 = 0;
    // A sequencia aberta: nome, entrada, e se a linha curta (a ultima) ja passou.
    let mut atual: Option<(String, Entrada, bool)> = None;
    let fecha = |atual: &mut Option<(String, Entrada, bool)>, seqs: &mut HashMap<String, Entrada>| -> Result<(), String> {
        if let Some((nome, mut e, _)) = atual.take() {
            if e.bases_por_linha == 0 {
                e.bases_por_linha = 1;
                e.bytes_por_linha = 1;
            }
            if seqs.insert(nome.clone(), e).is_some() {
                return Err(format!("{caminho}: the sequence name `{nome}` appears twice"));
            }
        }
        Ok(())
    };
    loop {
        linha.clear();
        let n = r.read_until(b'\n', &mut linha).map_err(|e| format!("Error reading {caminho}: {e}"))?;
        if n == 0 {
            break;
        }
        let inicio = byte;
        byte += n as u64;
        if linha[0] == b'>' {
            fecha(&mut atual, &mut seqs)?;
            atual = Some((
                nome_da_sequencia(&linha),
                Entrada { comprimento: 0, offset: byte, bases_por_linha: 0, bytes_por_linha: 0 },
                false,
            ));
            continue;
        }
        let mut bases = n;
        while bases > 0 && matches!(linha[bases - 1], b'\n' | b'\r') {
            bases -= 1;
        }
        let Some((_, e, curta)) = atual.as_mut() else {
            if bases == 0 {
                continue;
            }
            return Err(format!("{caminho}: sequence before the first `>` header line: this is not a FASTA file"));
        };
        if linha[..bases].iter().any(|b| b.is_ascii_whitespace()) {
            return Ok(None);
        }
        let (bases, n) = (bases as u64, n as u64);
        let com_fim_de_linha = linha.last() == Some(&b'\n');
        if e.bases_por_linha == 0 {
            if bases == 0 {
                // Linha vazia antes de qualquer base.
                return Ok(None);
            }
            e.offset = inicio;
            e.bases_por_linha = bases;
            e.bytes_por_linha = n;
            if !com_fim_de_linha {
                // A unica linha, no fim do arquivo e sem quebra: a largura
                // de uma linha cheia ainda e desconhecida, mas nao ha outra.
                e.bytes_por_linha = bases + 1;
            }
        } else if *curta || bases > e.bases_por_linha {
            // Depois da linha curta nao pode vir outra; e nenhuma e maior que a primeira.
            return Ok(None);
        } else if bases < e.bases_por_linha {
            *curta = true;
        } else if n != e.bytes_por_linha && com_fim_de_linha {
            // Mesma quantidade de bases, quebra de linha diferente (LF x CRLF).
            return Ok(None);
        }
        e.comprimento += bases;
    }
    fecha(&mut atual, &mut seqs)?;
    Ok(Some(seqs))
}

impl Genoma {
    /// Abre o genoma de referencia. `dir_temp`: onde copiar as bases quando
    /// o arquivo nao pode ser lido direto (comprimido, ou com linhas
    /// irregulares).
    pub fn abre(caminho: &str, dir_temp: Option<&Path>) -> Result<Genoma, String> {
        if !comprimido(caminho) {
            let seqs = match le_fai(caminho)? {
                Some(s) => Some(s),
                None => indexa_texto(caminho)?,
            };
            if let Some(seqs) = seqs {
                if seqs.is_empty() {
                    return Err(format!("{caminho}: no `>` sequence in the reference: is this a FASTA file?"));
                }
                let arq = File::open(caminho).map_err(|e| format!("Error opening {caminho}: {e}"))?;
                return Ok(Genoma {
                    caminho: caminho.to_string(),
                    arq,
                    seqs,
                    nomes: HashMap::new(),
                    janela: (String::new(), 0, Vec::new()),
                    _temp: None,
                });
            }
        }
        Genoma::copia(caminho, dir_temp)
    }

    /// Copia so as bases para um temporario: uma sequencia depois da outra,
    /// sem quebras de linha.
    fn copia(caminho: &str, dir_temp: Option<&Path>) -> Result<Genoma, String> {
        let mut temp = DirTemp::novo(dir_temp)?;
        let destino = temp.arquivo("reference.seq");
        let mut w = BufWriter::with_capacity(16 << 20, File::create(&destino).map_err(|e| format!("Error creating {destino}: {e}"))?);
        let mut r = abre_leitura(caminho)?;
        let mut seqs: HashMap<String, Entrada> = HashMap::new();
        let mut linha: Vec<u8> = Vec::new();
        let mut escritos: u64 = 0;
        let mut atual: Option<(String, Entrada)> = None;
        let erro_escrita = |e: std::io::Error| format!("Error writing the temporary copy of {caminho} ({destino}): {e}");
        let mut fecha = |atual: &mut Option<(String, Entrada)>| -> Result<(), String> {
            if let Some((nome, e)) = atual.take() {
                if seqs.insert(nome.clone(), e).is_some() {
                    return Err(format!("{caminho}: the sequence name `{nome}` appears twice"));
                }
            }
            Ok(())
        };
        loop {
            linha.clear();
            let n = r.read_until(b'\n', &mut linha).map_err(|e| format!("Error reading {caminho}: {e}"))?;
            if n == 0 {
                break;
            }
            if linha[0] == b'>' {
                fecha(&mut atual)?;
                atual = Some((
                    nome_da_sequencia(&linha),
                    Entrada { comprimento: 0, offset: escritos, bases_por_linha: u64::MAX, bytes_por_linha: u64::MAX },
                ));
                continue;
            }
            let bases: Vec<u8> = linha.iter().copied().filter(|b| !b.is_ascii_whitespace()).collect();
            if bases.is_empty() {
                continue;
            }
            let Some((_, e)) = atual.as_mut() else {
                return Err(format!("{caminho}: sequence before the first `>` header line: this is not a FASTA file"));
            };
            w.write_all(&bases).map_err(erro_escrita)?;
            escritos += bases.len() as u64;
            e.comprimento += bases.len() as u64;
        }
        fecha(&mut atual)?;
        w.flush().map_err(erro_escrita)?;
        drop(w);
        if seqs.is_empty() {
            return Err(format!("{caminho}: no `>` sequence in the reference: is this a FASTA file?"));
        }
        let arq = File::open(&destino).map_err(|e| format!("Error opening {destino}: {e}"))?;
        Ok(Genoma {
            caminho: caminho.to_string(),
            arq,
            seqs,
            nomes: HashMap::new(),
            janela: (String::new(), 0, Vec::new()),
            _temp: Some(temp),
        })
    }

    /// O nome que `chrom` tem no FASTA: o mesmo, ou o da outra convencao
    /// (`chr1` x `1`, `chrM` x `MT`). `None` se nao esta la.
    pub fn nome(&mut self, chrom: &str) -> Option<String> {
        if let Some(n) = self.nomes.get(chrom) {
            return n.clone();
        }
        let achado = if self.seqs.contains_key(chrom) {
            Some(chrom.to_string())
        } else {
            [Estilo::ComChr, Estilo::SemChr]
                .into_iter()
                .filter_map(|e| converte_nome(chrom, e))
                .chain(if chrom == "chrM" { Some("chrMT".to_string()) } else { None })
                .find(|n| self.seqs.contains_key(n))
        };
        self.nomes.insert(chrom.to_string(), achado.clone());
        achado
    }

    pub fn comprimento(&self, nome: &str) -> Option<u64> {
        self.seqs.get(nome).map(|e| e.comprimento)
    }

    /// As bases `[ini, fim)` (0-based) de `nome`, em maiusculas; menos que o
    /// pedido no fim da sequencia. `nome` e um nome do FASTA (ver `nome`).
    pub fn trecho(&mut self, nome: &str, ini: u64, fim: u64) -> Result<&[u8], String> {
        let Some(&e) = self.seqs.get(nome) else {
            return Ok(&[]);
        };
        let fim = fim.min(e.comprimento);
        if ini >= fim {
            return Ok(&[]);
        }
        let (ref j_nome, j_ini, ref j) = self.janela;
        let cobre = j_nome == nome && j_ini <= ini && fim <= j_ini + j.len() as u64;
        if !cobre {
            // Metade da janela de cada lado: o alinhamento anda para a
            // esquerda, e o VCF ordenado anda para a direita.
            let a = ini.saturating_sub(JANELA / 2);
            let b = fim.max(ini + JANELA / 2).min(e.comprimento);
            let bases = self.le(&e, a, b)?;
            self.janela = (nome.to_string(), a, bases);
        }
        let (_, j_ini, ref j) = self.janela;
        Ok(&j[(ini - j_ini) as usize..(fim - j_ini) as usize])
    }

    fn le(&mut self, e: &Entrada, a: u64, b: u64) -> Result<Vec<u8>, String> {
        let (de, ate) = (e.byte(a), e.byte(b - 1) + 1);
        let mut bruto = vec![0u8; (ate - de) as usize];
        self.arq
            .seek(SeekFrom::Start(de))
            .and_then(|_| self.arq.read_exact(&mut bruto))
            .map_err(|err| format!("Error reading the reference {}: {err}", self.caminho))?;
        bruto.retain(|c| *c != b'\n' && *c != b'\r');
        bruto.make_ascii_uppercase();
        if bruto.len() as u64 != b - a {
            return Err(format!(
                "The .fai index does not match {}: rebuild it with `samtools faidx`, or delete it.",
                self.caminho
            ));
        }
        Ok(bruto)
    }
}
