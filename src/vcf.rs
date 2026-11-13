// Ordenacao de VCF pela interface original do projeto ("POS" / "QUAL").
//
// O trabalho esta em `externo.rs`: ordena em memoria quando o dado cabe no
// orcamento e em lotes no disco, com merge k-vias, quando nao cabe. Nos dois
// casos a ordenacao e ESTAVEL -- variantes com a mesma chave saem na ordem em
// que entraram (ordem dos arquivos, depois ordem das linhas).
//
// Estas funcoes mantem a ordem alfabetica de cromossomos que o projeto sempre
// teve (chr1, chr10, chr11, chr2). Para a ordem natural ou a do cabecalho
// `##contig`, que e a do `bcftools sort`, use `externo::ordenar_arquivos`.
use crate::externo::{ordenar_arquivos, Criterio, Formato, OrdemCromossomos, Resumo};
use crate::orcamento::opcoes;

/// Ordena um ou mais VCFs por "POS" (cromossomo e posicao, crescente) ou
/// "QUAL" (qualidade, decrescente) e grava em `output_path`. O orcamento de
/// memoria e calculado a partir da RAM livre. Devolve o numero de variantes.
pub fn processar_e_ordenar_vcfs(
    caminhos: Vec<String>,
    campo_ordenacao: &str,
    output_path: &str,
) -> Result<usize, String> {
    processar_e_ordenar_vcfs_com_orcamento(caminhos, campo_ordenacao, output_path, None, None)
        .map(|r| r.registros)
}

/// Igual a [`processar_e_ordenar_vcfs`], com o orcamento de memoria (em MB) e
/// o diretorio dos lotes temporarios escolhidos por quem chama.
pub fn processar_e_ordenar_vcfs_com_orcamento(
    caminhos: Vec<String>,
    campo_ordenacao: &str,
    output_path: &str,
    memoria_mb: Option<usize>,
    dir_temp: Option<&str>,
) -> Result<Resumo, String> {
    let criterio = match campo_ordenacao {
        "POS" => Criterio::Posicao,
        "QUAL" => Criterio::Qualidade,
        _ => return Err("Unknown sort field (use \"POS\" or \"QUAL\")".to_string()),
    };
    let op = opcoes(memoria_mb, dir_temp, OrdemCromossomos::Lexica);
    println!(
        "🧬 RUST: starting VCF Multimerge of {} file(s), memory budget {} MB...",
        caminhos.len(),
        op.memoria_bytes >> 20
    );
    let resumo = ordenar_arquivos(&caminhos, output_path, Formato::Vcf, criterio, &op)?;
    if resumo.lotes == 0 {
        println!("✅ RUST: {} variants sorted in memory.", resumo.registros);
    } else {
        println!(
            "✅ RUST: {} variants sorted out of core ({} runs on disk, k-way merge).",
            resumo.registros, resumo.lotes
        );
    }
    Ok(resumo)
}
