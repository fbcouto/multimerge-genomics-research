use std::fs::File;
use std::io::{BufRead, BufReader};
use rayon::prelude::*;
use flate2::read::MultiGzDecoder;
use crate::kmer::Kmer128;

/// Processa múltiplos ficheiros FASTA (suporta .gz e texto limpo) em paralelo
pub fn processar_fasta_lote(caminhos: Vec<String>, k: u8) -> Vec<u128> {
    println!("🚀 Starting multi-threaded extraction from {} FASTA file(s)...", caminhos.len());

        let mut todos_kmers: Vec<u128> = caminhos.par_iter().flat_map(|caminho| {
        let file = File::open(caminho).expect(&format!("Error opening file: {}", caminho));
        
        // Inteligência de Descompressão
        let mut reader: Box<dyn BufRead> = if caminho.ends_with(".gz") {
            Box::new(BufReader::with_capacity(1024 * 1024 * 16, MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::with_capacity(1024 * 1024 * 16, file))
        };

        let mut kmers_locais = Vec::new();
        let mut linha = String::new();
        let mut seq_atual = String::new();

        while let Ok(bytes_lidos) = reader.read_line(&mut linha) {
            if bytes_lidos == 0 { break; }
            
            let l = linha.trim();
            if l.starts_with('>') {
                // Quando acha um novo cabeçalho, processa a sequência acumulada anterior
                if seq_atual.len() >= k as usize {
                    for janela in seq_atual.as_bytes().windows(k as usize) {
                        if !janela.contains(&b'N') && !janela.contains(&b'n') {
                            kmers_locais.push(Kmer128::encode(janela).bits);
                        }
                    }
                }
                seq_atual.clear();
            } else {
                seq_atual.push_str(l);
            }
            linha.clear();
        }
        
        // Não esquecer de processar a última sequência do ficheiro!
        if seq_atual.len() >= k as usize {
            for janela in seq_atual.as_bytes().windows(k as usize) {
                if !janela.contains(&b'N') && !janela.contains(&b'n') {
                    kmers_locais.push(Kmer128::encode(janela).bits);
                }
            }
        }
        
        kmers_locais
    }).collect();

    println!("🔄 Sorting {} k-mers from the FASTA batch...", todos_kmers.len());
    crate::multimerge_core::multi_merge_sort_stable(&mut todos_kmers);
    todos_kmers.dedup(); 

    todos_kmers
}