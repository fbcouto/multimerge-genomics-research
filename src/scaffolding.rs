use std::collections::HashMap;

/// Recebe a lista de milhares de contigs e o tamanho do "encaixe" (overlap)
/// Retorna os Super-Contigs montados e a contagem de Genomas Circulares encontrados.
pub fn fast_hash_scaffold(contigs: Vec<String>, k_overlap: usize) -> (Vec<String>, usize) {
    println!("🛠️ RUST: Iniciando Scaffolding de {} peças...", contigs.len());
    
    let mut prefix_map: HashMap<String, usize> = HashMap::with_capacity(contigs.len());
    let mut next_link = vec![usize::MAX; contigs.len()];
    let mut in_degree = vec![0; contigs.len()];

    // 1. Mapear o INÍCIO de cada peça (Prefixos)
    for (i, c) in contigs.iter().enumerate() {
        if c.len() >= k_overlap {
            let prefix = c[..k_overlap].to_string();
            prefix_map.insert(prefix, i);
        }
    }

    // 2. Procurar se o FINAL de uma peça encaixa no INÍCIO de outra
    for (i, c) in contigs.iter().enumerate() {
        if c.len() >= k_overlap {
            let suffix = &c[c.len() - k_overlap..];
            if let Some(&j) = prefix_map.get(suffix) {
                // Não permite ligar a si mesmo no primeiro passo e garante linha única
                if i != j && in_degree[j] == 0 && next_link[i] == usize::MAX {
                    next_link[i] = j;
                    in_degree[j] += 1;
                }
            }
        }
    }

    // 3. Montagem Final (O Andaime)
    let mut scaffolds = Vec::new();
    let mut circulares = 0;
    let mut visited = vec![false; contigs.len()];

    // 3.1. Costurar caminhos lineares
    for i in 0..contigs.len() {
        if in_degree[i] == 0 && !visited[i] {
            let mut curr = i;
            let mut super_contig = contigs[curr].clone();
            visited[curr] = true;

            while next_link[curr] != usize::MAX {
                let nxt = next_link[curr];
                if visited[nxt] { break; } // Segurança máxima contra loops infinitos
                visited[nxt] = true;
                
                // Cola o pedaço, omitindo o trecho que se sobrepõe para não duplicar o DNA!
                super_contig.push_str(&contigs[nxt][k_overlap..]);
                curr = nxt;
            }
            
            // 🚨 TESTE DE DNA CIRCULAR! 🚨
            // Verifica se, depois de todo o caminho longo, a "cauda" encaixa na "cabeça" original
            if super_contig.len() >= k_overlap * 2 {
                let pref = &super_contig[..k_overlap];
                let suff = &super_contig[super_contig.len() - k_overlap..];
                if pref == suff {
                    circulares += 1;
                    // Limpa a duplicação final do círculo perfeito
                    let clean_len = super_contig.len() - k_overlap;
                    super_contig.truncate(clean_len);
                }
            }
            scaffolds.push(super_contig);
        }
    }

    // 3.2. Capturar pequenos anéis biológicos (Plasmídeos) que ficaram isolados
    for i in 0..contigs.len() {
        if !visited[i] {
            let mut curr = i;
            let mut anel = contigs[curr].clone();
            visited[curr] = true;

            while next_link[curr] != usize::MAX {
                let nxt = next_link[curr];
                if visited[nxt] {
                    circulares += 1;
                    let clean_len = anel.len() - k_overlap;
                    anel.truncate(clean_len);
                    break;
                }
                visited[nxt] = true;
                anel.push_str(&contigs[nxt][k_overlap..]);
                curr = nxt;
            }
            scaffolds.push(anel);
        }
    }

    // Ordenar do maior para o menor
    scaffolds.sort_by(|a, b| b.len().cmp(&a.len()));
    (scaffolds, circulares)
}