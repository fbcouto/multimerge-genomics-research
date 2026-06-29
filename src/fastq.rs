use std::fs::File;
use std::io::{BufRead, BufReader};
use rayon::prelude::*;
use flate2::read::MultiGzDecoder;
use crate::kmer::Kmer128;

/// Processa múltiplos ficheiros FASTQ (suporta .gz e texto limpo) em paralelo
pub fn processar_fastq_lote(caminhos: Vec<String>, qualidade_minima: u8, k: u8) -> Vec<u128> {
    println!("🚀 Iniciando extração Multi-Thread de {} ficheiro(s) FASTQ...", caminhos.len());

    // 1. O Rayon processa CADA ficheiro numa thread diferente simultaneamente
        let mut todos_kmers: Vec<u128> = caminhos.par_iter().flat_map(|caminho| {
        let file = File::open(caminho).expect(&format!("Erro ao abrir o ficheiro: {}", caminho));
        
        // 2. Inteligência de Descompressão On-The-Fly
        let mut reader: Box<dyn BufRead> = if caminho.ends_with(".gz") {
            Box::new(BufReader::with_capacity(1024 * 1024 * 16, MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::with_capacity(1024 * 1024 * 16, file))
        };

        let mut kmers_locais = Vec::new();
        let mut linha = String::new();
        let mut num_linha = 0;
        let mut seq = String::new();

        // 3. Leitura rápida do bloco FASTQ (4 linhas por sequência)
        while let Ok(bytes_lidos) = reader.read_line(&mut linha) {
            if bytes_lidos == 0 { break; }
            
            let idx = num_linha % 4;
            if idx == 1 { 
                seq = linha.trim().to_string(); // Linha 2: A Sequência de DNA
            } else if idx == 3 { 
                let qual = linha.trim();        // Linha 4: A Qualidade (Phred)
                
                if seq.len() >= k as usize && seq.len() == qual.len() {
                    // Cálculo rápido de qualidade média
                    let mut soma_qual: usize = 0;
                    for &q in qual.as_bytes() {
                        soma_qual += (q as usize).saturating_sub(33); // Tabela Phred+33
                    }
                    let media = soma_qual / qual.len();
                    
                    // Se passar no QC, fatiamos a sequência em K-mers de 128-bits
                    if media >= qualidade_minima as usize {
                        for janela in seq.as_bytes().windows(k as usize) {
                            if !janela.contains(&b'N') && !janela.contains(&b'n') {
                                kmers_locais.push(Kmer128::encode(janela).bits);
                            }
                        }
                    }
                }
            }
            linha.clear();
            num_linha += 1;
        }
        kmers_locais
    }).collect(); // O Rayon junta os resultados de todas as threads aqui!

    // 4. Esmagamento Final: Ordenação Instável In-Place (Ultra rápida)
    println!("🔄 Ordenando {} K-mers do lote e removendo duplicatas...", todos_kmers.len());
    todos_kmers.par_sort_unstable();
    todos_kmers.dedup(); // Remove cópias exatas do DNA para economizar RAM no Grafo

    todos_kmers
}