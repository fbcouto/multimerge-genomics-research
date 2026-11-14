// =====================================================================
//  OPERACOES SOBRE INTERVALOS (BED)
//
//  As duas perguntas mais frequentes sobre intervalos genomicos:
//
//    fundir      "quais regioes este conjunto cobre?"      (bedtools merge)
//    intersecao  "o que do conjunto A cai sobre o B?"      (bedtools intersect)
//
//  As duas sao varreduras lineares sobre dado ORDENADO, e o custo real esta
//  na ordenacao. Aqui ela e a de `externo.rs`: a entrada nao precisa vir
//  ordenada nem caber na memoria.
//
//  Coordenadas BED: inicio em base zero, fim exclusivo. Dois intervalos se
//  sobrepoem quando max(inicios) < min(fins). Intervalos encostados
//  (fim de um == inicio do outro) NAO se sobrepoem, mas sao fundidos.
// =====================================================================

use crate::externo::{
    abre_leitura, ordenar, ordenar_arquivos, Criterio, Destino, DirTemp, Escritor, Formato, Opcoes, OrdemNomes,
};
use std::cmp::Ordering;
use std::io::BufRead;

fn campos(linha: &str) -> Option<(&str, u64, u64)> {
    let mut c = linha.split('\t');
    Some((c.next()?, c.next()?.parse().ok()?, c.next()?.parse().ok()?))
}

// ---------------------------------------------------------------------
//  FUNDIR
// ---------------------------------------------------------------------

/// Recebe intervalos ordenados e emite a uniao. Roda como destino da
/// ordenacao: o BED ordenado nunca chega a existir em disco.
struct Fusor {
    w: Option<Escritor>,
    distancia: u64,
    atual: Option<(String, u64, u64)>,
    saidas: usize,
}

impl Fusor {
    fn descarrega(&mut self) -> Result<(), String> {
        if let Some((c, i, f)) = self.atual.take() {
            self.w.as_mut().ok_or("output already closed")?.linha(&format!("{c}\t{i}\t{f}"))?;
            self.saidas += 1;
        }
        Ok(())
    }
}

impl Destino for Fusor {
    // Linhas `track` e `browser` descrevem o arquivo original, nao a uniao.
    fn cabecalho(&mut self, _linha: &str) -> Result<(), String> {
        Ok(())
    }

    fn registro(&mut self, linha: &str) -> Result<(), String> {
        let (c, i, f) = campos(linha).ok_or("invalid interval in the sorted output")?;
        match &mut self.atual {
            Some((ca, _, fa)) if ca == c && i <= fa.saturating_add(self.distancia) => {
                *fa = (*fa).max(f);
            }
            _ => {
                self.descarrega()?;
                self.atual = Some((c.to_string(), i, f));
            }
        }
        Ok(())
    }

    fn fim(&mut self) -> Result<(), String> {
        self.descarrega()?;
        self.w.take().ok_or("output already closed")?.termina()
    }
}

/// Funde os intervalos de um ou mais BED: sobrepostos, encostados, ou
/// separados por ate `distancia` bases. Grava `cromossomo inicio fim`.
/// Devolve (intervalos lidos, intervalos gravados).
pub fn bed_fundir(caminhos: &[String], saida: &str, distancia: u64, op: &Opcoes) -> Result<(usize, usize), String> {
    let mut fusor = Fusor { w: Some(Escritor::cria_para(saida, Formato::Bed)?), distancia, atual: None, saidas: 0 };
    let r = ordenar(caminhos, Formato::Bed, Criterio::Posicao, op, &mut fusor)?;
    Ok((r.registros, fusor.saidas))
}

// ---------------------------------------------------------------------
//  INTERSECAO
// ---------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ModoIntersecao {
    /// Para cada sobreposicao entre um intervalo de A e um de B, a linha de A
    /// com as coordenadas recortadas ao trecho comum. E o padrao do
    /// `bedtools intersect`.
    Regioes,
    /// Cada linha de A que toca B, uma vez, inteira (`-u`).
    Unicos,
    /// Cada linha de A que NAO toca B (`-v`).
    Ausentes,
}

impl ModoIntersecao {
    pub fn de_texto(s: &str) -> Option<ModoIntersecao> {
        match s.to_ascii_lowercase().as_str() {
            "regions" => Some(ModoIntersecao::Regioes),
            "unique" | "u" => Some(ModoIntersecao::Unicos),
            "absent" | "v" => Some(ModoIntersecao::Ausentes),
            _ => None,
        }
    }
}

/// Leitor de um BED ja ordenado, com um intervalo de antecipacao.
struct Fluxo {
    leitor: Box<dyn BufRead>,
    linha: String,
}

impl Fluxo {
    /// Carrega o proximo intervalo em `self.linha`. `false` no fim.
    fn proximo(&mut self) -> Result<bool, String> {
        self.linha.clear();
        if self.leitor.read_line(&mut self.linha).map_err(|e| e.to_string())? == 0 {
            return Ok(false);
        }
        if self.linha.ends_with('\n') {
            self.linha.pop();
        }
        Ok(true)
    }
}

/// Cruza A com B. As entradas nao precisam vir ordenadas nem caber na
/// memoria: cada lado e ordenado para um arquivo temporario, e uma varredura
/// unica os percorre juntos. A saida sai na ordem de A. Devolve o numero de
/// linhas gravadas.
///
/// A memoria da varredura e o conjunto de intervalos de B que ainda podem
/// tocar o A corrente -- proporcional a profundidade de sobreposicao, nao ao
/// tamanho dos arquivos.
pub fn bed_intersecao(
    a: &[String],
    b: &[String],
    saida: &str,
    modo: ModoIntersecao,
    op: &Opcoes,
) -> Result<usize, String> {
    let mut temp = DirTemp::novo(op.dir_temp.as_deref())?;
    let (arq_a, arq_b) = (temp.arquivo("a_ordenado.bed.lz4"), temp.arquivo("b_ordenado.bed.lz4"));
    ordenar_arquivos(a, &arq_a, Formato::Bed, Criterio::Posicao, op)?;
    ordenar_arquivos(b, &arq_b, Formato::Bed, Criterio::Posicao, op)?;

    let ordem = OrdemNomes::nova(&op.ordem);
    let mut fa = Fluxo { leitor: abre_leitura(&arq_a)?, linha: String::new() };
    let mut fb = Fluxo { leitor: abre_leitura(&arq_b)?, linha: String::new() };
    let mut w = Escritor::cria_para(saida, Formato::Bed)?;

    // Intervalos de B, do cromossomo corrente, que ainda podem tocar um A.
    let mut ativos: Vec<(u64, u64)> = Vec::new();
    let mut crom_ativos = String::new();
    let mut tem_b = fb.proximo()?;
    let mut gravadas = 0usize;
    // Nomes de sequencia vistos de cada lado, para a checagem do fim.
    let mut nomes_a: Vec<String> = Vec::new();
    let mut nomes_b: Vec<String> = Vec::new();

    while fa.proximo()? {
        // Cabecalhos nao chegam aqui: a ordenacao os poe no topo, e o `#`
        // do BED faz `campos` falhar na coordenada.
        let Some((ca, ia, fim_a)) = campos(&fa.linha) else { continue };
        if ca != crom_ativos {
            ativos.clear();
            crom_ativos.clear();
            crom_ativos.push_str(ca);
            nomes_a.push(ca.to_string());
        }
        // Traz de B tudo que comeca antes do fim deste A.
        while tem_b {
            let Some((cb, ib, fim_b)) = campos(&fb.linha) else {
                tem_b = fb.proximo()?;
                continue;
            };
            if nomes_b.last().map(String::as_str) != Some(cb) {
                nomes_b.push(cb.to_string());
            }
            match ordem.compara(cb, ca) {
                Ordering::Less => {}
                Ordering::Equal if ib < fim_a => ativos.push((ib, fim_b)),
                _ => break,
            }
            tem_b = fb.proximo()?;
        }
        // A esta ordenado por inicio: quem termina antes deste A comecar nao
        // toca este nem nenhum dos proximos.
        ativos.retain(|&(_, fim_b)| fim_b > ia);

        let mut tocou = false;
        for &(ib, fim_b) in &ativos {
            let (i, f) = (ia.max(ib), fim_a.min(fim_b));
            if i < f {
                tocou = true;
                if modo != ModoIntersecao::Regioes {
                    break;
                }
                // A linha de A, com as coordenadas recortadas.
                let resto = fa.linha.splitn(4, '\t').nth(3);
                match resto {
                    Some(r) => w.linha(&format!("{ca}\t{i}\t{f}\t{r}"))?,
                    None => w.linha(&format!("{ca}\t{i}\t{f}"))?,
                }
                gravadas += 1;
            }
        }
        let emite = match modo {
            ModoIntersecao::Regioes => false,
            ModoIntersecao::Unicos => tocou,
            ModoIntersecao::Ausentes => !tocou,
        };
        if emite {
            w.linha(&fa.linha)?;
            gravadas += 1;
        }
    }
    w.termina()?;

    // ERRO SILENCIOSO DE NOMENCLATURA. Se A e B nao tem nenhuma sequencia em
    // comum, o resultado acima e "nada se cruza" -- formalmente correto e
    // quase sempre errado: um arquivo chama o cromossomo de `chr1` e o outro
    // de `1`, ou os dois vem de versoes diferentes do genoma de referencia.
    // Devolver um arquivo vazio sem dizer nada e o pior desfecho, entao a
    // operacao falha e explica. (Como B e lido ate passar do ultimo
    // cromossomo de A, qualquer nome em comum ja teria aparecido.)
    if !nomes_a.is_empty() && !nomes_b.is_empty() && !nomes_a.iter().any(|n| nomes_b.contains(n)) {
        let _ = std::fs::remove_file(saida);
        let (ea, eb) = (&nomes_a[0], &nomes_b[0]);
        let dica = if ea.starts_with("chr") != eb.starts_with("chr") {
            "One file uses the `chr` prefix and the other does not: rename the sequences of one of them before intersecting (normalize_chromosomes does it)."
        } else {
            "Check that both files come from the same version of the reference genome (GRCh37/hg19 and GRCh38/hg38, for instance, do not mix)."
        };
        return Err(format!(
            "A and B have no sequence name in common, so nothing was intersected. A uses names like `{ea}`; B uses names like `{eb}`. {dica}"
        ));
    }
    Ok(gravadas)
}
