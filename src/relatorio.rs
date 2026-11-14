// =====================================================================
//  RELATORIO DE QUALIDADE DE FASTQ EM HTML
//
//  Uma pagina so, sem nada de fora: o modelo (`fastq_report.html`) leva os
//  estilos, o codigo que desenha os graficos em SVG e os numeros, que entram
//  como JSON. Abre em qualquer navegador, sem rede, e serve de anexo.
//
//  O que a pagina mostra e o que `fastq_qc::FastqStats` mediu: resumo,
//  avisos, qualidade e composicao ao longo da leitura, qualidade e GC por
//  leitura, comprimentos e adaptadores -- cada grafico com a tabela dos
//  mesmos numeros.
// =====================================================================

use crate::fastq_qc::FastqStats;

const MODELO: &str = include_str!("fastq_report.html");

/// A pagina do relatorio para estas estatisticas.
pub fn fastq_html(stats: &FastqStats) -> String {
    // Os numeros vao dentro de um <script>: um `<` no caminho do arquivo ou
    // num aviso nao pode fechar o bloco nem abrir outra marca.
    let dados = stats.json().replace('<', "\\u003c");
    MODELO.replace("__MMGEN_DATA__", &dados).replace("__MMGEN_VERSION__", env!("CARGO_PKG_VERSION"))
}

/// Grava o relatorio em `saida`.
pub fn grava_fastq_html(stats: &FastqStats, saida: &str) -> Result<(), String> {
    std::fs::write(saida, fastq_html(stats)).map_err(|e| format!("Error writing {saida}: {e}"))
}
