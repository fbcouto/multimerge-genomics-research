use std::collections::{HashMap, HashSet};
use crate::kmer::Kmer128;

pub fn construir_grafo_estatisticas(kmers: &[u128], k: u8) -> (usize, usize) {
    let mut grafo: HashMap<u128, Vec<u128>> = HashMap::with_capacity(kmers.len());
    for &bits in kmers {
        let kmer = Kmer128 { bits, k };
        let prefixo = kmer.prefix();
        let sufixo = kmer.suffix();
        grafo.entry(prefixo).or_insert_with(Vec::new).push(sufixo);
    }
    let total_nos = grafo.len();
    let total_arestas: usize = grafo.values().map(|v| v.len()).sum();
    (total_nos, total_arestas)
}

pub fn montar_todos_contigs(kmers: &[u128], k: u8) -> Vec<String> {
    println!("🧩 Extracting contigs (u128, K={})...", k);
    
    let mut grafo: HashMap<u128, Vec<u128>> = HashMap::with_capacity(kmers.len());
    let mut in_degree: HashMap<u128, usize> = HashMap::with_capacity(kmers.len());
    
    for &bits in kmers {
        let kmer = Kmer128 { bits, k };
        let prefixo = kmer.prefix();
        let sufixo = kmer.suffix();
        grafo.entry(prefixo).or_insert_with(Vec::new).push(sufixo);
        *in_degree.entry(sufixo).or_insert(0) += 1;
        in_degree.entry(prefixo).or_insert(0);
    }
    
    let mut contigs = Vec::new();
    let mut arestas_visitadas = HashSet::new();
    let mut bolhas_colapsadas = 0;
    
    for (&u, destinos) in &grafo {
        let out_u = destinos.len();
        let in_u = *in_degree.get(&u).unwrap_or(&0);
        
        if !(in_u == 1 && out_u == 1) && out_u > 0 {
            for &v in destinos {
                if arestas_visitadas.insert((u, v)) {
                    extrair_caminho_linear(u, v, &grafo, &in_degree, &mut arestas_visitadas, &mut contigs, &mut bolhas_colapsadas, k);
                }
            }
        }
    }
    
    for (&u, destinos) in &grafo {
        for &v in destinos {
            if arestas_visitadas.insert((u, v)) {
                extrair_caminho_linear(u, v, &grafo, &in_degree, &mut arestas_visitadas, &mut contigs, &mut bolhas_colapsadas, k);
            }
        }
    }
    
    println!("🫧 Done: {} simple bubble(s) collapsed.", bolhas_colapsadas);
    contigs.sort_by(|a, b| b.len().cmp(&a.len()));
    contigs
}

fn extrair_caminho_linear(
    inicio_u: u128,
    inicio_v: u128,
    grafo: &HashMap<u128, Vec<u128>>,
    in_degree: &HashMap<u128, usize>,
    visitados: &mut HashSet<(u128, u128)>,
    contigs: &mut Vec<String>,
    bolhas_colapsadas: &mut u32,
    k: u8
) {
    let mut caminho_atual = vec![inicio_u, inicio_v];
    let mut atual = inicio_v;
    
    while let Some(proximos) = grafo.get(&atual) {
        let out_atual = proximos.len();
        let in_atual = *in_degree.get(&atual).unwrap_or(&0);
        
        if in_atual == 1 && out_atual == 1 {
            let proximo_no = proximos[0];
            if visitados.insert((atual, proximo_no)) {
                caminho_atual.push(proximo_no);
                atual = proximo_no;
                continue;
            } else { break; }
        } 
        else if out_atual == 2 {
            let v1 = proximos[0];
            let v2 = proximos[1];
            
            if let (Some(d1), Some(d2)) = (grafo.get(&v1), grafo.get(&v2)) {
                if d1.len() == 1 && d2.len() == 1 && d1[0] == d2[0] {
                    let no_convergencia = d1[0];
                    visitados.insert((atual, v1));
                    visitados.insert((atual, v2));
                    visitados.insert((v1, no_convergencia));
                    visitados.insert((v2, no_convergencia));
                    
                    caminho_atual.push(v1);
                    caminho_atual.push(no_convergencia);
                    atual = no_convergencia;
                    *bolhas_colapsadas += 1;
                    continue;
                }
            }
            break;
        } else {
            break;
        }
    }
    
    // Decodifica a sequência genética usando a nossa estrutura Kmer128 (Lembrando que os nós têm k-1)
    let mut dna = String::with_capacity(caminho_atual.len() + (k as usize));
    let primeiro_no = Kmer128 { bits: caminho_atual[0], k: k - 1 };
    dna.push_str(&primeiro_no.decode());
    
    for &no in &caminho_atual[1..] {
        let seq_no = Kmer128 { bits: no, k: k - 1 };
        dna.push(seq_no.decode().chars().last().unwrap());
    }
    contigs.push(dna);
}