use rayon::prelude::*;

pub fn translate_dna_to_protein(dna: &str) -> String {
    // Garante que todas as letras estão maiúsculas para o dicionário não falhar
    let seq_upper = dna.to_ascii_uppercase();
    let bytes = seq_upper.as_bytes();
    
    // Otimização: A proteína terá exatamente 1/3 do tamanho do DNA
    let mut protein = String::with_capacity(bytes.len() / 3);

    for chunk in bytes.chunks(3) {
        if chunk.len() < 3 { break; } // Ignora pedaços incompletos no final
        
        // O Dicionário Universal da Vida (Tabela de Códons)
        let aa = match chunk {
            b"ATA" | b"ATC" | b"ATT" => 'I',
            b"ATG" => 'M', // Start Codon
            b"ACA" | b"ACC" | b"ACG" | b"ACT" => 'T',
            b"AAC" | b"AAT" => 'N',
            b"AAA" | b"AAG" => 'K',
            b"AGC" | b"AGT" | b"TCA" | b"TCC" | b"TCG" | b"TCT" => 'S',
            b"AGA" | b"AGG" | b"CGA" | b"CGC" | b"CGG" | b"CGT" => 'R',
            b"CTA" | b"CTC" | b"CTG" | b"CTT" | b"TTA" | b"TTG" => 'L',
            b"CCA" | b"CCC" | b"CCG" | b"CCT" => 'P',
            b"CAC" | b"CAT" => 'H',
            b"CAA" | b"CAG" => 'Q',
            b"GTA" | b"GTC" | b"GTG" | b"GTT" => 'V',
            b"GCA" | b"GCC" | b"GCG" | b"GCT" => 'A',
            b"GAC" | b"GAT" => 'D',
            b"GAA" | b"GAG" => 'E',
            b"GGA" | b"GGC" | b"GGG" | b"GGT" => 'G',
            b"TTC" | b"TTT" => 'F',
            b"TAC" | b"TAT" => 'Y',
            b"TGC" | b"TGT" => 'C',
            b"TGG" => 'W',
            b"TAA" | b"TAG" | b"TGA" => '*', // Stop Codons (Fim da proteína)
            _ => 'X', // X = Ambíguo ou erro de sequenciamento
        };
        protein.push(aa);
    }
    protein
}

/// Recebe milhares de fragmentos de DNA e traduz todos ao mesmo tempo usando a CPU inteira
pub fn translate_batch(sequences: &[String]) -> Vec<String> {
    sequences.par_iter()
        .map(|seq| translate_dna_to_protein(seq))
        .collect()
}