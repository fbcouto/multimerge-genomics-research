use pyo3::prelude::*;
use pyo3::types::PyBytes;
use rayon::prelude::*; 
use std::fs::File;
use std::io::{BufRead, BufReader};

// Declaração ÚNICA de todos os nossos submódulos
mod kmer; 
mod fastq;
mod fasta;
mod debruijn;
mod annotation;
mod scaffolding;
mod vcf;
pub mod multimerge;

// ---------------------------------------------------------
// ROTA 1: O Motor da GPU (CUDA)
// ---------------------------------------------------------
#[cfg(cuda_enabled)]
extern "C" {
    fn launch_merge_path_kernel(a: *const u64, b: *const u64, out: *mut u64, len_a: usize, len_b: usize);
    fn launch_bloom_filter_build(kmers: *const u64, num_kmers: usize, bloom: *mut u8, bloom_bytes: usize, num_hashes: i32);
    fn launch_bloom_filter_query(query_kmers: *const u64, num_queries: usize, bloom: *const u8, bloom_bytes: usize, num_hashes: i32, out_mask: *mut u8);
}

#[cfg(cuda_enabled)]
fn ordenar_e_mesclar_backend(amostra_a: &[u64], amostra_b: &[u64]) -> Vec<u64> {
    println!("🚀 Enviando dados para a VRAM da NVIDIA GPU (CUDA)...");
    let mut resultado = vec![0; amostra_a.len() + amostra_b.len()];
    
    unsafe {
        launch_merge_path_kernel(
            amostra_a.as_ptr(),
            amostra_b.as_ptr(),
            resultado.as_mut_ptr(),
            amostra_a.len(),
            amostra_b.len()
        );
    }
    resultado
}

#[cfg(not(cuda_enabled))] // Fallback CPU automático caso não haja GPU disponível
fn ordenar_e_mesclar_backend(amostra_a: &[u64], amostra_b: &[u64]) -> Vec<u64> {
    println!("⚙️ Executando via CPU Multi-Core (Rayon Fallback)...");
    let mut resultado = Vec::with_capacity(amostra_a.len() + amostra_b.len());
    resultado.extend_from_slice(amostra_a);
    resultado.extend_from_slice(amostra_b);
    resultado.par_sort_unstable();
    resultado
}

// ---------------------------------------------------------
// PIPELINE ARROW: RECEBE DADOS LIMPOS DO MÓDULO 1
// ---------------------------------------------------------
#[pyfunction]
fn ordenar_bytes_gpu<'py>(py: Python<'py>, bytes_a: &[u8], bytes_b: &[u8]) -> PyResult<&'py PyBytes> {
    println!("🚀 Recebendo blocos Arrow IPC para a GPU...");

    let mut kmer_a: Vec<u64> = unsafe {
        std::slice::from_raw_parts(bytes_a.as_ptr() as *const u64, bytes_a.len() / 8).to_vec()
    };
    let mut kmer_b: Vec<u64> = unsafe {
        std::slice::from_raw_parts(bytes_b.as_ptr() as *const u64, bytes_b.len() / 8).to_vec()
    };

    kmer_a.par_sort_unstable();
    kmer_b.par_sort_unstable();

    let resultado = ordenar_e_mesclar_backend(&kmer_a, &kmer_b);

    let byte_slice: &[u8] = unsafe {
        std::slice::from_raw_parts(resultado.as_ptr() as *const u8, resultado.len() * 8)
    };
    Ok(PyBytes::new(py, byte_slice))
}

// ---------------------------------------------------------
// MÓDULO 1: QC E LIMPEZA FASTQ/FASTA EM LOTE (BATCH U128)
// ---------------------------------------------------------
#[pyfunction]
fn ler_fastq_em_lote<'py>(py: Python<'py>, caminhos: Vec<String>, qualidade_minima: u8, k: u8) -> PyResult<&'py PyBytes> {
    let kmers = fastq::processar_fastq_lote(caminhos, qualidade_minima, k);
    
    // MAGIA MONOLÍTICA (128-bits = 16 bytes)
    let byte_slice: &[u8] = unsafe {
        std::slice::from_raw_parts(kmers.as_ptr() as *const u8, kmers.len() * 16) 
    };
    Ok(PyBytes::new(py, byte_slice))
}

#[pyfunction]
fn ler_fasta_em_lote<'py>(py: Python<'py>, caminhos: Vec<String>, k: u8) -> PyResult<&'py PyBytes> {
    let kmers = fasta::processar_fasta_lote(caminhos, k);
    
    let byte_slice: &[u8] = unsafe {
        std::slice::from_raw_parts(kmers.as_ptr() as *const u8, kmers.len() * 16)
    };
    Ok(PyBytes::new(py, byte_slice))
}

// ---------------------------------------------------------
// MÓDULO 3: GRAFOS DE DE BRUIJN (128-BITS)
// ---------------------------------------------------------
#[pyfunction]
fn construir_grafo_debruijn(bytes_kmer: &[u8], k: u8) -> PyResult<(usize, usize)> {
    let kmers: &[u128] = unsafe {
        std::slice::from_raw_parts(bytes_kmer.as_ptr() as *const u128, bytes_kmer.len() / 16)
    };
    Ok(debruijn::construir_grafo_estatisticas(kmers, k))
}

#[pyfunction]
fn montar_genoma_debruijn(bytes_kmer: &[u8], k: u8) -> PyResult<Vec<String>> {
    let kmers: &[u128] = unsafe {
        std::slice::from_raw_parts(bytes_kmer.as_ptr() as *const u128, bytes_kmer.len() / 16)
    };
    Ok(debruijn::montar_todos_contigs(kmers, k))
}

// ---------------------------------------------------------
// MÓDULO 4 e 6: FILTRO BLOOM & CHAMADA DE VARIANTES (GPU)
// ---------------------------------------------------------
#[pyfunction]
fn construir_filtro_bloom_gpu<'py>(py: Python<'py>, bytes_kmer: &[u8], tamanho_mb: usize, num_hashes: i32) -> PyResult<&'py PyBytes> {
    let kmers_128: &[u128] = unsafe {
        std::slice::from_raw_parts(bytes_kmer.as_ptr() as *const u128, bytes_kmer.len() / 16)
    };

    let kmers_gpu: Vec<u64> = kmers_128.par_iter().map(|&k| ((k >> 64) as u64) ^ (k as u64)).collect();

    let bloom_bytes = tamanho_mb * 1024 * 1024;
    let mut bloom_array = vec![0u8; bloom_bytes]; 

    #[cfg(cuda_enabled)]
    unsafe {
        println!("🔥 CUDA: Inserindo {} K-mers (Comprimidos) no Filtro...", kmers_gpu.len());
        launch_bloom_filter_build(kmers_gpu.as_ptr(), kmers_gpu.len(), bloom_array.as_mut_ptr(), bloom_bytes, num_hashes);
    }
    Ok(PyBytes::new(py, &bloom_array))
}

#[pyfunction]
fn extrair_mutacoes_gpu<'py>(py: Python<'py>, amostra_bytes: &[u8], bloom_bytes: &[u8], num_hashes: i32) -> PyResult<&'py PyBytes> {
    let kmers_128: &[u128] = unsafe {
        std::slice::from_raw_parts(amostra_bytes.as_ptr() as *const u128, amostra_bytes.len() / 16)
    };

    let kmers_gpu: Vec<u64> = kmers_128.par_iter().map(|&k| ((k >> 64) as u64) ^ (k as u64)).collect();

    let mut mascara_resultados = vec![0u8; kmers_gpu.len()];

    #[cfg(cuda_enabled)]
    unsafe {
        println!("⚠️ CUDA: Vasculhando {} K-mers em busca de mutações...", kmers_gpu.len());
        launch_bloom_filter_query(kmers_gpu.as_ptr(), kmers_gpu.len(), bloom_bytes.as_ptr(), bloom_bytes.len(), num_hashes, mascara_resultados.as_mut_ptr());
    }

    let kmers_mutantes: Vec<u128> = kmers_128.par_iter()
        .zip(mascara_resultados.par_iter())
        .filter(|(_, &existe_na_ref)| existe_na_ref == 0)
        .map(|(&k, _)| k)
        .collect();

    let bytes_mutantes: &[u8] = unsafe {
        std::slice::from_raw_parts(kmers_mutantes.as_ptr() as *const u8, kmers_mutantes.len() * 16)
    };
    Ok(PyBytes::new(py, bytes_mutantes))
}

// ---------------------------------------------------------
// MÓDULO 5: ANOTAÇÃO FUNCIONAL
// ---------------------------------------------------------
#[pyfunction]
fn traduzir_proteinas(sequencias_dna: Vec<String>) -> PyResult<Vec<String>> {
    println!("⚙️ RAYON: Traduzindo {} sequências de DNA em Proteínas...", sequencias_dna.len());
    let proteinas = annotation::translate_batch(&sequencias_dna);
    Ok(proteinas)
}

// ---------------------------------------------------------
// MÓDULO 7: SCAFFOLDING
// ---------------------------------------------------------
#[pyfunction]
fn montar_scaffolds(sequencias: Vec<String>, overlap_k: usize) -> PyResult<(Vec<String>, usize)> {
    let (scaffolds, circulares) = scaffolding::fast_hash_scaffold(sequencias, overlap_k);
    Ok((scaffolds, circulares))
}

// ---------------------------------------------------------
// DESCODIFICADORES VISUAIS E TESTES
// ---------------------------------------------------------
#[pyfunction]
fn decodificar_lote_128(kmers_bytes: &[u8], k: u8) -> PyResult<Vec<String>> {
    let kmers: &[u128] = unsafe {
        std::slice::from_raw_parts(kmers_bytes.as_ptr() as *const u128, kmers_bytes.len() / 16)
    };
    let seqs: Vec<String> = kmers.par_iter().map(|&bits| {
        let kmer = kmer::Kmer128 { bits, k };
        kmer.decode()
    }).collect();
    Ok(seqs)
}

#[pyfunction]
fn testar_kmer_gigante(sequencia_dna: &str) -> PyResult<String> {
    let kmer_binario = kmer::Kmer128::encode(sequencia_dna.as_bytes());
    let dna_recuperado = kmer_binario.decode();
    Ok(format!("Sucesso! DNA Original: {} | Bits na RAM: {} | DNA Recuperado: {}", sequencia_dna, kmer_binario.bits, dna_recuperado))
}

// ---------------------------------------------------------
// MÓDULOS LEGADOS (MANTIDOS PARA COMPATIBILIDADE U64)
// ---------------------------------------------------------
#[pyfunction]
fn ler_arquivos_e_processar<'py>(py: Python<'py>, caminho_a: &str, caminho_b: &str) -> PyResult<&'py PyBytes> {
    let ler_kmer = |caminho: &str| -> Vec<u64> {
        let file = File::open(caminho).unwrap();
        let reader = BufReader::with_capacity(1024 * 1024 * 16, file);
        reader.lines().filter_map(|l| l.ok()).map(|s| kmer::encode_kmer(s.as_bytes())).filter(|&k| k != u64::MAX).collect()
    };
    let mut kmer_a = ler_kmer(caminho_a);
    let mut kmer_b = ler_kmer(caminho_b);
    kmer_a.par_sort_unstable();
    kmer_b.par_sort_unstable();
    let resultado = ordenar_e_mesclar_backend(&kmer_a, &kmer_b);
    let byte_slice: &[u8] = unsafe { std::slice::from_raw_parts(resultado.as_ptr() as *const u8, resultado.len() * 8) };
    Ok(PyBytes::new(py, byte_slice))
}

#[pyfunction]
fn gerar_assinaturas_lote<'py>(py: Python<'py>, kmers_bytes: &[u8], num_hashes: usize) -> PyResult<&'py PyBytes> {
    let kmers: &[u64] = unsafe { std::slice::from_raw_parts(kmers_bytes.as_ptr() as *const u64, kmers_bytes.len() / 8) };
    let assinaturas: Vec<u64> = kmers.par_iter().flat_map(|&k| kmer::compute_minhash_signature(k, num_hashes)).collect();
    let byte_slice: &[u8] = unsafe { std::slice::from_raw_parts(assinaturas.as_ptr() as *const u8, assinaturas.len() * 8) };
    Ok(PyBytes::new(py, byte_slice))
}

#[pyfunction]
fn decodificar_lote(kmers_bytes: &[u8], length: usize) -> PyResult<Vec<String>> {
    let kmers: &[u64] = unsafe { std::slice::from_raw_parts(kmers_bytes.as_ptr() as *const u64, kmers_bytes.len() / 8) };
    let seqs: Vec<String> = kmers.par_iter().map(|&k| kmer::decode_kmer(k, length as u8)).collect();
    Ok(seqs)
}
// ---------------------------------------------------------
// MÓDULO 8: VCF MULTIMERGE & SORTING
// ---------------------------------------------------------
#[pyfunction]
fn multimerge_vcf(caminhos: Vec<String>, campo_ordenacao: &str, output_path: &str) -> PyResult<usize> {
    match vcf::processar_e_ordenar_vcfs(caminhos, campo_ordenacao, output_path) {
        Ok(total) => Ok(total),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

// ---------------------------------------------------------
// O SUMÁRIO (Exporta todas as funções para o Python)
// ---------------------------------------------------------
#[pymodule]
fn multimerge_genomics(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ler_arquivos_e_processar, m)?)?;
    m.add_function(wrap_pyfunction!(decodificar_lote, m)?)?;
    m.add_function(wrap_pyfunction!(gerar_assinaturas_lote, m)?)?;
    m.add_function(wrap_pyfunction!(ler_fastq_em_lote, m)?)?; 
    m.add_function(wrap_pyfunction!(ler_fasta_em_lote, m)?)?;
    m.add_function(wrap_pyfunction!(ordenar_bytes_gpu, m)?)?;
    m.add_function(wrap_pyfunction!(construir_grafo_debruijn, m)?)?;
    m.add_function(wrap_pyfunction!(montar_genoma_debruijn, m)?)?;
    m.add_function(wrap_pyfunction!(construir_filtro_bloom_gpu, m)?)?;
    m.add_function(wrap_pyfunction!(extrair_mutacoes_gpu, m)?)?;
    m.add_function(wrap_pyfunction!(traduzir_proteinas, m)?)?;
    m.add_function(wrap_pyfunction!(montar_scaffolds, m)?)?;
    m.add_function(wrap_pyfunction!(testar_kmer_gigante, m)?)?;
    m.add_function(wrap_pyfunction!(decodificar_lote_128, m)?)?;
    m.add_function(wrap_pyfunction!(multimerge_vcf, m)?)?;
    Ok(())
}