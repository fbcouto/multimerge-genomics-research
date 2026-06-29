use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use flate2::read::MultiGzDecoder;

/// Estrutura leve para manter a linha inteira na RAM e apenas as chaves necessárias para a ordenação
#[derive(Clone)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: u64,
    pub qual: f32,
    pub raw_line: String,
}

pub fn processar_e_ordenar_vcfs(caminhos: Vec<String>, campo_ordenacao: &str, output_path: &str) -> Result<usize, String> {
    println!("🧬 RUST: Iniciando VCF Multimerge de {} ficheiro(s)...", caminhos.len());

    // 1. LER E EXTRAIR EM PARALELO
    let mut todos_registos: Vec<VcfRecord> = caminhos.par_iter().flat_map(|caminho| {
        let file = File::open(caminho).expect("Erro ao abrir VCF");
        let reader: Box<dyn BufRead> = if caminho.ends_with(".gz") {
            Box::new(BufReader::with_capacity(1024 * 1024 * 16, MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::with_capacity(1024 * 1024 * 16, file))
        };

        let mut registos_locais = Vec::new();
        for linha_result in reader.lines() {
            if let Ok(linha) = linha_result {
                // Ignora os cabeçalhos (mas numa versão final poderíamos guardá-los)
                if linha.starts_with('#') { continue; }

                let colunas: Vec<&str> = linha.split('\t').collect();
                if colunas.len() >= 6 {
                    let chrom = colunas[0].to_string();
                    let pos = colunas[1].parse::<u64>().unwrap_or(0);
                    let qual = colunas[5].parse::<f32>().unwrap_or(0.0);

                    registos_locais.push(VcfRecord {
                        chrom,
                        pos,
                        qual,
                        raw_line: linha,
                    });
                }
            }
        }
        registos_locais
    }).collect();

    let total_variantes = todos_registos.len();
    println!("🔄 RUST: Ordenando {:#?} variantes em memória...", total_variantes);

    // 2. O PODER DO MULTIMERGE (Ordenação Baseada no Campo Escolhido)
    match campo_ordenacao {
        "POS" => {
            // Ordena primeiro por Cromossoma (Alfabético) e depois por Posição (Numérico)
            todos_registos.par_sort_unstable_by(|a, b| {
                a.chrom.cmp(&b.chrom).then_with(|| a.pos.cmp(&b.pos))
            });
        },
        "QUAL" => {
            // Ordena por Qualidade da mutação (do maior para o menor)
            todos_registos.par_sort_unstable_by(|a, b| {
                b.qual.partial_cmp(&a.qual).unwrap_or(std::cmp::Ordering::Equal)
            });
        },
        _ => return Err("Campo de ordenação desconhecido!".to_string()),
    }

    // 3. GRAVAÇÃO SUPER RÁPIDA NO DISCO
    println!("💾 RUST: Escrevendo super-ficheiro VCF em {}...", output_path);
    let out_file = File::create(output_path).map_err(|e| e.to_string())?;
    let mut writer = BufWriter::with_capacity(1024 * 1024 * 32, out_file);

    // Escreve um cabeçalho VCF mínimo
    writeln!(writer, "##fileformat=VCFv4.2").unwrap();
    writeln!(writer, "##source=MMGEN_Multimerge_Engine").unwrap();
    writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT").unwrap();

    for reg in todos_registos {
        writeln!(writer, "{}", reg.raw_line).unwrap();
    }

    Ok(total_variantes)
}