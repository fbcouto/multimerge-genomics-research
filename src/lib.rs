use pyo3::prelude::*;
use pyo3::types::PyBytes;
use rayon::prelude::*;
use std::borrow::Cow;
use std::fs::File;
use std::io::{BufRead, BufReader};

// Declaração ÚNICA de todos os nossos submódulos
mod kmer;
mod fastq;
mod fasta;
mod debruijn;
mod annotation;
mod scaffolding;
pub mod vcf;
pub mod multimerge;
// O motor de ordenacao (copia de adaptive-parallel-multimerge-sort) e a
// ordenacao de VCF de qualquer tamanho que roda sobre ele.
pub mod multimerge_core;
pub mod externo;
pub mod intervalos;
pub mod orcamento;
pub mod nomes;
pub mod validador;
// Filtro de Bloom na CPU (gemeo do kernel CUDA) e a ponte para a GPU.
pub mod bloom;
pub mod gpu;
// Envelopes: BGZF (o .gz que se indexa), BAM, os indices tabix e BAI com a
// consulta por regiao, e Parquet.
pub mod bam;
pub mod bgzf;
pub mod indice;
#[cfg(feature = "parquet")]
pub mod parquet_io;
// Mais operacoes de intervalos, ferramentas de VCF, qualidade e aparo de
// FASTQ (com o relatorio em HTML) e contagem de k-mers de qualquer tamanho.
pub mod bed_ops;
pub mod fastq_qc;
pub mod kmer_count;
pub mod referencia;
pub mod relatorio;
pub mod vcf_tools;
// Genoma humano (passo 15): achados do ClinVar, regioes de genes, consulta de
// variantes em parentes, pileup; e a entrada do PharmCAT, que chama os alelos
// de farmacogenetica (o MMGEN so prepara o VCF e le o resultado dele).
pub mod clinico;
pub mod farmacogenetica;
pub mod pharmcat;
// Parentesco (KING), trio (Mendel, de novo, paternidade), regioes de
// homozigose e sexo genetico.
pub mod parentesco;
// Paineis de frequencias por populacao, ascendencia (EM supervisionado) e
// escore poligenico (arquivos do PGS Catalog).
pub mod ancestralidade;
// Pesquisa: contagem de leituras por gene e expressao diferencial (DESeq2).
pub mod expressao;
pub mod conjuntos;
// RNA-seq direto do FASTQ: pseudoalinhamento e EM, como o kallisto.
pub mod quantificacao;
// Triagens CRISPR: contagem de guias e teste de genes (MAGeCK, RRA).
pub mod crispr;
// Microbiologia: esbocos FracMinHash (como o sourmash), especies e
// abundancia, classificacao de leituras, genes de resistencia (AMRFinderPlus).
pub mod microbio;
// Surto e vigilancia: distancias e arvores (NJ, UPGMA), mutacoes de um virus
// contra a referencia (como o Nextclade), PCR in silico.
pub mod surto;
// 16S: pares -> ASVs (UNOISE3, quimeras), taxonomia (classificador do RDP),
// diversidade alfa e beta, PCoA, PERMANOVA; e a carga mutacional do tumor.
pub mod amplicon;
// Coorte: associacao genotipo x desfecho (linear, logistica, Firth,
// permutacoes, FDR) e poder do estudo.
pub mod associacao;
// Triagem de agentes adventicios (virus, bacterias, micoplasmas, fungos)
// nas leituras de um material, com hospedeiro, brancos e spike-ins.
pub mod contaminantes;

// A API PYTHON E EM INGLES (2026-11-13): nomes das funcoes, parametros, chaves
// dos dicionarios, valores aceitos, mensagens e a documentacao que o `help()`
// mostra. Os nomes internos do Rust e os comentarios continuam em portugues.

// ---------------------------------------------------------
// ROTA 1: GPU (CUDA) com a CPU de reserva
// ---------------------------------------------------------
// As funcoes que tem kernel CUDA recebem `device`:
//   "auto" (padrao)  o caminho mais rapido que esta maquina tem.
//   "gpu"            so a GPU; levanta erro se ela nao puder (para medir e
//                    para ter certeza de que foi ela que rodou).
//   "cpu"            so a CPU.
// Os dois caminhos devolvem o mesmo resultado, bit a bit. Antes, sem CUDA, o
// filtro de Bloom saia zerado e a consulta dava todo k-mer como mutacao.
#[derive(Clone, Copy, PartialEq)]
enum Dispositivo {
    Auto,
    Cpu,
    Gpu,
}

impl Dispositivo {
    fn de_texto(texto: &str) -> PyResult<Self> {
        match texto.to_ascii_lowercase().as_str() {
            "auto" => Ok(Self::Auto),
            "cpu" => Ok(Self::Cpu),
            "gpu" | "cuda" => Ok(Self::Gpu),
            outro => Err(erro_py(format!("Unknown device: {outro} (use \"auto\", \"gpu\" or \"cpu\")"))),
        }
    }
}

/// Roda `na_gpu`; em "auto", se a GPU nao puder, roda `na_cpu`.
fn despacha<T>(
    onde: Dispositivo,
    tarefa: &str,
    na_gpu: impl FnOnce() -> Result<T, String>,
    na_cpu: impl FnOnce() -> T,
) -> PyResult<T> {
    match onde {
        Dispositivo::Cpu => Ok(na_cpu()),
        Dispositivo::Gpu => na_gpu().map_err(|e| erro_py(format!("{tarefa}: GPU unavailable ({e})"))),
        Dispositivo::Auto => match na_gpu() {
            Ok(v) => Ok(v),
            Err(e) => {
                // Sem placa nenhuma nao ha o que avisar: a CPU e o caminho normal.
                if gpu::disponivel() {
                    eprintln!("⚠️ {tarefa}: the GPU could not handle it ({e}); falling back to the CPU.");
                }
                Ok(na_cpu())
            }
        },
    }
}

/// Le os bytes vindos do Python como palavras de 64 bits. So copia se o
/// ponteiro vier desalinhado (um `bytes` do CPython vem alinhado, mas isso nao
/// e promessa de ninguem). Bytes que sobram no fim sao ignorados.
fn como_u64(bytes: &[u8]) -> Cow<'_, [u64]> {
    // SAFETY: todo padrao de bits e um u64 valido, e `align_to` so poe no meio
    // o trecho corretamente alinhado.
    let (antes, meio, _) = unsafe { bytes.align_to::<u64>() };
    if antes.is_empty() {
        Cow::Borrowed(meio)
    } else {
        Cow::Owned(bytes.chunks_exact(8).map(|c| u64::from_ne_bytes(c.try_into().unwrap())).collect())
    }
}

/// O mesmo para k-mers de 128 bits.
fn como_u128(bytes: &[u8]) -> Cow<'_, [u128]> {
    // SAFETY: como em `como_u64`.
    let (antes, meio, _) = unsafe { bytes.align_to::<u128>() };
    if antes.is_empty() {
        Cow::Borrowed(meio)
    } else {
        Cow::Owned(bytes.chunks_exact(16).map(|c| u128::from_ne_bytes(c.try_into().unwrap())).collect())
    }
}

/// Funde dois blocos ordenados na CPU.
fn mesclar_cpu(amostra_a: &[u64], amostra_b: &[u64]) -> Vec<u64> {
    let mut resultado = Vec::with_capacity(amostra_a.len() + amostra_b.len());
    resultado.extend_from_slice(amostra_a);
    resultado.extend_from_slice(amostra_b);
    // Dois runs ja ordenados: o motor os detecta e faz uma fusao so.
    multimerge_core::multi_merge_sort_stable(&mut resultado);
    resultado
}

fn ordenar_e_mesclar_backend(amostra_a: &[u64], amostra_b: &[u64], onde: Dispositivo) -> PyResult<Vec<u64>> {
    // Em "auto" a fusao fica na CPU. Medido em 2026-11-13 numa A100 com 12
    // nucleos (ensaios/gpu_teste.py), a GPU perde em todo tamanho -- 2 x 1M:
    // +4 ms; 2 x 10M: +25 ms; 2 x 50M: +62 ms -- porque levar o dado ate a
    // placa e trazer de volta custa mais que fundir. A GPU so funde quando
    // pedida por nome. Nao inverter sem medir de novo.
    if onde == Dispositivo::Auto {
        return Ok(mesclar_cpu(amostra_a, amostra_b));
    }
    despacha(onde, "merge", || gpu::mesclar(amostra_a, amostra_b), || mesclar_cpu(amostra_a, amostra_b))
}

// ---------------------------------------------------------
// PIPELINE ARROW: RECEBE DADOS LIMPOS DO MÓDULO 1
// ---------------------------------------------------------
/// Sorts two blocks of 64-bit keys and merges them into a single sorted block.
///
/// `bytes_a` and `bytes_b` are the raw bytes of two arrays of unsigned 64-bit
/// integers (native byte order); the result has the same layout. `device`:
/// "auto" (the merge runs on the CPU, which measured faster), "gpu" (use the
/// CUDA merge-path kernel, raising an error if no GPU can run it) or "cpu".
#[pyfunction]
#[pyo3(signature = (bytes_a, bytes_b, device="auto"))]
fn sort_merge_bytes<'py>(py: Python<'py>, bytes_a: &[u8], bytes_b: &[u8], device: &str) -> PyResult<Bound<'py, PyBytes>> {
    let onde = Dispositivo::de_texto(device)?;
    // Sem o GIL: a interface (marimo) continua respondendo enquanto ordena.
    let resultado = py.allow_threads(|| {
        let mut kmer_a: Vec<u64> = como_u64(bytes_a).into_owned();
        let mut kmer_b: Vec<u64> = como_u64(bytes_b).into_owned();
        multimerge_core::multi_merge_sort_stable(&mut kmer_a);
        multimerge_core::multi_merge_sort_stable(&mut kmer_b);
        ordenar_e_mesclar_backend(&kmer_a, &kmer_b, onde)
    })?;

    let byte_slice: &[u8] = unsafe {
        std::slice::from_raw_parts(resultado.as_ptr() as *const u8, resultado.len() * 8)
    };
    Ok(PyBytes::new(py, byte_slice))
}

// ---------------------------------------------------------
// MÓDULO 1: QC E LIMPEZA FASTQ/FASTA EM LOTE (BATCH U128)
// ---------------------------------------------------------
// Confere as entradas dos leitores de k-mers antes de começar: um caminho
// errado era um `expect` dentro do Rayon (pânico, que no Python chega como
// PanicException e as telas não pegavam), e k = 0 entrava em pânico no
// `windows(0)`. Agora são erros comuns, com o nome do arquivo.
fn confere_entradas_kmer(paths: &[String], k: u8) -> PyResult<()> {
    if k == 0 || k > 63 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!("k must be between 1 and 63 (got {k})")));
    }
    if paths.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err("no input file given"));
    }
    for p in paths {
        if let Err(e) = std::fs::File::open(p) {
            return Err(pyo3::exceptions::PyOSError::new_err(format!("cannot open {p}: {e}")));
        }
    }
    Ok(())
}

/// Reads one or more FASTQ files (plain text or `.gz`) and returns their
/// distinct k-mers, sorted, as raw bytes: 16 bytes (one unsigned 128-bit
/// integer) per k-mer.
///
/// Reads whose mean Phred quality is below `min_quality` are dropped. `k` goes
/// up to 63. The files are read in parallel.
#[pyfunction]
fn read_fastq_batch<'py>(py: Python<'py>, paths: Vec<String>, min_quality: u8, k: u8) -> PyResult<Bound<'py, PyBytes>> {
    confere_entradas_kmer(&paths, k)?;
    let kmers = py.allow_threads(|| fastq::processar_fastq_lote(paths, min_quality, k));

    // MAGIA MONOLÍTICA (128-bits = 16 bytes)
    let byte_slice: &[u8] = unsafe {
        std::slice::from_raw_parts(kmers.as_ptr() as *const u8, kmers.len() * 16)
    };
    Ok(PyBytes::new(py, byte_slice))
}

/// Reads one or more FASTA files (plain text or `.gz`) and returns their
/// distinct k-mers, sorted, as raw bytes: 16 bytes (one unsigned 128-bit
/// integer) per k-mer. Windows containing `N` are skipped. `k` goes up to 63.
#[pyfunction]
fn read_fasta_batch<'py>(py: Python<'py>, paths: Vec<String>, k: u8) -> PyResult<Bound<'py, PyBytes>> {
    confere_entradas_kmer(&paths, k)?;
    let kmers = py.allow_threads(|| fasta::processar_fasta_lote(paths, k));

    let byte_slice: &[u8] = unsafe {
        std::slice::from_raw_parts(kmers.as_ptr() as *const u8, kmers.len() * 16)
    };
    Ok(PyBytes::new(py, byte_slice))
}

// ---------------------------------------------------------
// MÓDULO 3: GRAFOS DE DE BRUIJN (128-BITS)
// ---------------------------------------------------------
/// Builds the De Bruijn graph of a set of k-mers (16 bytes each, as returned
/// by `read_fasta_batch` / `read_fastq_batch`) and returns `(nodes, edges)`.
#[pyfunction]
fn build_debruijn_graph(py: Python<'_>, kmer_bytes: &[u8], k: u8) -> PyResult<(usize, usize)> {
    let kmers = como_u128(kmer_bytes);
    Ok(py.allow_threads(|| debruijn::construir_grafo_estatisticas(&kmers, k)))
}

/// Assembles contigs from a set of k-mers by walking the unbranched paths of
/// the De Bruijn graph. Returns the contigs as strings, longest first.
#[pyfunction]
fn assemble_debruijn(py: Python<'_>, kmer_bytes: &[u8], k: u8) -> PyResult<Vec<String>> {
    let kmers = como_u128(kmer_bytes);
    Ok(py.allow_threads(|| debruijn::montar_todos_contigs(&kmers, k)))
}

// ---------------------------------------------------------
// MÓDULO 4 e 6: FILTRO BLOOM & CHAMADA DE VARIANTES (GPU)
// ---------------------------------------------------------
/// Builds the Bloom filter of the reference k-mers: `size_mb` megabytes,
/// `num_hashes` bits per k-mer. Returns the filter as bytes.
///
/// Runs on the GPU when there is one, otherwise on the CPU; the filter is
/// bit-for-bit identical either way. `device`: "auto", "gpu" (raises an error
/// if no GPU can run it) or "cpu".
#[pyfunction]
#[pyo3(signature = (kmer_bytes, size_mb, num_hashes, device="auto"))]
fn build_bloom_filter<'py>(
    py: Python<'py>,
    kmer_bytes: &[u8],
    size_mb: usize,
    num_hashes: i32,
    device: &str,
) -> PyResult<Bound<'py, PyBytes>> {
    let onde = Dispositivo::de_texto(device)?;
    if size_mb == 0 || num_hashes <= 0 {
        return Err(erro_py("The Bloom filter needs at least 1 MB and 1 hash".to_string()));
    }
    let bloom_array = py.allow_threads(|| {
        let kmers_128 = como_u128(kmer_bytes);
        let kmers_gpu: Vec<u64> = kmers_128.par_iter().map(|&k| bloom::dobra(k)).collect();
        let bloom_bytes = size_mb * 1024 * 1024;
        despacha(
            onde,
            "Bloom filter",
            || gpu::bloom_construir(&kmers_gpu, bloom_bytes, num_hashes),
            || bloom::construir(&kmers_gpu, bloom_bytes, num_hashes),
        )
    })?;
    Ok(PyBytes::new(py, &bloom_array))
}

/// Returns the k-mers of the sample that are NOT in the reference filter (the
/// candidate mutations), in the order they came, as raw bytes (16 per k-mer).
/// `num_hashes` must be the value used to build the filter. `device`: "auto",
/// "gpu" or "cpu".
#[pyfunction]
#[pyo3(signature = (sample_bytes, bloom_bytes, num_hashes, device="auto"))]
fn extract_mutations<'py>(
    py: Python<'py>,
    sample_bytes: &[u8],
    bloom_bytes: &[u8],
    num_hashes: i32,
    device: &str,
) -> PyResult<Bound<'py, PyBytes>> {
    let onde = Dispositivo::de_texto(device)?;
    if bloom_bytes.len() < 4 || num_hashes <= 0 {
        return Err(erro_py("Empty Bloom filter or invalid number of hashes".to_string()));
    }
    let kmers_mutantes: Vec<u128> = py.allow_threads(|| -> PyResult<Vec<u128>> {
        let kmers_128 = como_u128(sample_bytes);
        let kmers_gpu: Vec<u64> = kmers_128.par_iter().map(|&k| bloom::dobra(k)).collect();
        let mascara_resultados = despacha(
            onde,
            "Bloom filter query",
            || gpu::bloom_consultar(&kmers_gpu, bloom_bytes, num_hashes),
            || bloom::consultar(&kmers_gpu, bloom_bytes, num_hashes),
        )?;
        Ok(kmers_128.par_iter()
            .zip(mascara_resultados.par_iter())
            .filter(|(_, &existe_na_ref)| existe_na_ref == 0)
            .map(|(&k, _)| k)
            .collect())
    })?;

    let bytes_mutantes: &[u8] = unsafe {
        std::slice::from_raw_parts(kmers_mutantes.as_ptr() as *const u8, kmers_mutantes.len() * 16)
    };
    Ok(PyBytes::new(py, bytes_mutantes))
}

// ---------------------------------------------------------
// MÓDULO 5: ANOTAÇÃO FUNCIONAL
// ---------------------------------------------------------
/// Translates DNA sequences into proteins (standard genetic code, reading
/// frame 1; `*` marks a stop codon). One protein per input sequence.
#[pyfunction]
fn translate_proteins(py: Python<'_>, dna_sequences: Vec<String>) -> PyResult<Vec<String>> {
    println!("⚙️ RAYON: translating {} DNA sequence(s) into proteins...", dna_sequences.len());
    Ok(py.allow_threads(|| annotation::translate_batch(&dna_sequences)))
}

// ---------------------------------------------------------
// MÓDULO 7: SCAFFOLDING
// ---------------------------------------------------------
/// Joins contigs whose ends overlap exactly by `overlap` bases. Returns
/// `(scaffolds, circular)`: the joined sequences and how many of them closed
/// into a circle.
#[pyfunction]
fn build_scaffolds(py: Python<'_>, sequences: Vec<String>, overlap: usize) -> PyResult<(Vec<String>, usize)> {
    Ok(py.allow_threads(|| scaffolding::fast_hash_scaffold(sequences, overlap)))
}

// ---------------------------------------------------------
// DESCODIFICADORES VISUAIS E TESTES
// ---------------------------------------------------------
/// Decodes 128-bit k-mers (16 bytes each) of size `k` back into DNA strings.
#[pyfunction]
fn decode_kmers_128(py: Python<'_>, kmer_bytes: &[u8], k: u8) -> PyResult<Vec<String>> {
    let kmers = como_u128(kmer_bytes);
    Ok(py.allow_threads(|| {
        kmers.par_iter().map(|&bits| kmer::Kmer128 { bits, k }.decode()).collect()
    }))
}

/// Encodes a DNA string of up to 63 bases into a 128-bit k-mer and decodes it
/// back; returns a one-line report. A quick self-test of the encoding.
#[pyfunction]
fn test_giant_kmer(dna_sequence: &str) -> PyResult<String> {
    let kmer_binario = kmer::Kmer128::encode(dna_sequence.as_bytes());
    let dna_recuperado = kmer_binario.decode();
    Ok(format!("Success! Original DNA: {} | Bits in RAM: {} | Recovered DNA: {}", dna_sequence, kmer_binario.bits, dna_recuperado))
}

// ---------------------------------------------------------
// MÓDULOS LEGADOS (MANTIDOS PARA COMPATIBILIDADE U64)
// ---------------------------------------------------------
/// Legacy 64-bit path. Reads two text files with one k-mer of up to 31 bases
/// per line, sorts each and merges them; returns the merged keys as raw bytes
/// (8 per k-mer).
#[pyfunction]
fn read_and_merge_kmer_files<'py>(py: Python<'py>, path_a: &str, path_b: &str) -> PyResult<Bound<'py, PyBytes>> {
    // Caminho errado vira erro com o nome do arquivo (antes era um `unwrap`,
    // ou seja, um panico).
    let ler_kmer = |caminho: &str| -> PyResult<Vec<u64>> {
        let file = File::open(caminho).map_err(|e| erro_py(format!("Error opening {caminho}: {e}")))?;
        let reader = BufReader::with_capacity(1024 * 1024 * 16, file);
        Ok(reader.lines().filter_map(|l| l.ok()).map(|s| kmer::encode_kmer(s.as_bytes())).filter(|&k| k != u64::MAX).collect())
    };
    let resultado = py.allow_threads(|| -> PyResult<Vec<u64>> {
        let mut kmer_a = ler_kmer(path_a)?;
        let mut kmer_b = ler_kmer(path_b)?;
        multimerge_core::multi_merge_sort_stable(&mut kmer_a);
        multimerge_core::multi_merge_sort_stable(&mut kmer_b);
        ordenar_e_mesclar_backend(&kmer_a, &kmer_b, Dispositivo::Auto)
    })?;
    let byte_slice: &[u8] = unsafe { std::slice::from_raw_parts(resultado.as_ptr() as *const u8, resultado.len() * 8) };
    Ok(PyBytes::new(py, byte_slice))
}

/// Legacy 64-bit path. Computes a MinHash signature of `num_hashes` values for
/// each 64-bit k-mer; returns them as raw bytes (`num_hashes` x 8 per k-mer).
#[pyfunction]
fn minhash_signatures<'py>(py: Python<'py>, kmer_bytes: &[u8], num_hashes: usize) -> PyResult<Bound<'py, PyBytes>> {
    let kmers = como_u64(kmer_bytes);
    let assinaturas: Vec<u64> =
        py.allow_threads(|| kmers.par_iter().flat_map(|&k| kmer::compute_minhash_signature(k, num_hashes)).collect());
    let byte_slice: &[u8] = unsafe { std::slice::from_raw_parts(assinaturas.as_ptr() as *const u8, assinaturas.len() * 8) };
    Ok(PyBytes::new(py, byte_slice))
}

/// Legacy 64-bit path. Decodes 64-bit k-mers (8 bytes each) of `length` bases
/// back into DNA strings.
#[pyfunction]
fn decode_kmers(py: Python<'_>, kmer_bytes: &[u8], length: usize) -> PyResult<Vec<String>> {
    let kmers = como_u64(kmer_bytes);
    Ok(py.allow_threads(|| kmers.par_iter().map(|&k| kmer::decode_kmer(k, length as u8)).collect()))
}
// ---------------------------------------------------------
// MÓDULO 8: VCF MULTIMERGE & SORTING
// ---------------------------------------------------------
/// Sorts one or more VCF files by "POS" or "QUAL" into a single output.
/// Stable, and for files of any size: what does not fit in the memory budget
/// is sorted in runs on disk and joined by a k-way merge.
///
/// `memory_mb` is the budget for the records; without it, half of the free
/// RAM. `temp_dir` is where the runs are written; without it, the system's
/// temporary directory. Returns the number of variants written.
///
/// The GIL is released while sorting, so a user interface does not freeze.
#[pyfunction]
#[pyo3(signature = (paths, sort_field, output_path, memory_mb=None, temp_dir=None))]
fn multimerge_vcf(
    py: Python<'_>,
    paths: Vec<String>,
    sort_field: &str,
    output_path: &str,
    memory_mb: Option<usize>,
    temp_dir: Option<&str>,
) -> PyResult<usize> {
    let resultado = py.allow_threads(|| {
        vcf::processar_e_ordenar_vcfs_com_orcamento(paths, sort_field, output_path, memory_mb, temp_dir)
    });
    match resultado {
        Ok(resumo) => Ok(resumo.registros),
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(e)),
    }
}

// ---------------------------------------------------------
// O SUMÁRIO (Exporta todas as funções para o Python)
// ---------------------------------------------------------
// ---------------------------------------------------------
// ORDENACAO E INTERVALOS DE QUALQUER TAMANHO
// ---------------------------------------------------------

fn erro_py(e: String) -> PyErr {
    pyo3::exceptions::PyRuntimeError::new_err(e)
}

/// Ordem padrao dos cromossomos de cada formato: a que a ferramenta de
/// referencia daquele formato usa. VCF e SAM seguem o cabecalho (bcftools,
/// samtools); BED e GFF seguem a ordem lexicografica (bedtools).
fn ordem_ou_padrao(texto: Option<&str>, formato: externo::Formato) -> Result<externo::OrdemCromossomos, String> {
    match texto {
        Some(t) => externo::OrdemCromossomos::de_texto(t),
        None => Ok(match formato {
            externo::Formato::Vcf | externo::Formato::Sam => externo::OrdemCromossomos::Cabecalho,
            externo::Formato::Bed | externo::Formato::Gff => externo::OrdemCromossomos::Lexica,
        }),
    }
}

/// Sorts one or more VCF, BED, GFF/GTF or SAM files of any size into a single
/// output. Stable: records with the same key come out in the order they went
/// in (file order, then line order).
///
/// - `paths`: a list of input files. Plain text, `.gz` (BGZF included) or
///   `.lz4`.
/// - `output`: the output file; its extension picks the compression.
/// - `format`: "vcf", "bed", "gff", "gtf" or "sam". Without it, taken from the
///   extension of the first file.
/// - `by`: "position" (every format), "quality" (VCF, descending) or "name"
///   (SAM, QNAME in natural order).
/// - `chrom_order`: "lexical", "natural", "header", or the path of a `.fai` /
///   chromosome-sizes file whose first column gives the order. Default: the
///   header for VCF and SAM, lexical for BED and GFF.
/// - `memory_mb`: budget for the records; without it, half of the free RAM.
/// - `temp_dir`: where to write the runs when the input does not fit the
///   budget.
/// - `index`: also build the index of the output (`.tbi` for a `.gz`, `.bai`
///   for a `.bam`), as `index_file` does. Needs `by="position"`.
///
/// Inputs and output can be plain text, `.gz` (BGZF), `.lz4`, `.bam` (for
/// SAM) or `.parquet`.
///
/// Returns `(records, runs_on_disk, skipped_lines)`.
#[pyfunction]
#[pyo3(signature = (paths, output, format=None, by="position", chrom_order=None, memory_mb=None, temp_dir=None, index=false))]
#[allow(clippy::too_many_arguments)]
fn sort_genomic(
    py: Python<'_>,
    paths: Vec<String>,
    output: &str,
    format: Option<&str>,
    by: &str,
    chrom_order: Option<&str>,
    memory_mb: Option<usize>,
    temp_dir: Option<&str>,
    index: bool,
) -> PyResult<(usize, usize, usize)> {
    py.allow_threads(|| {
        let formato = match format {
            Some(f) => externo::Formato::de_texto(f).ok_or_else(|| format!("Unknown format: {f}"))?,
            None => paths
                .first()
                .and_then(|c| externo::Formato::do_arquivo(c))
                .ok_or("Could not infer the format from the file extension; pass `format`")?,
        };
        let criterio = externo::Criterio::de_texto(by)
            .ok_or_else(|| format!("Unknown sort key: {by} (use \"position\", \"quality\" or \"name\")"))?;
        if index {
            // Antes de ordenar: descobrir isto depois de uma hora de trabalho
            // seria tarde.
            let s = output.to_ascii_lowercase();
            if criterio != externo::Criterio::Posicao {
                return Err("index=True needs by=\"position\": an index only works on a file sorted by position".to_string());
            }
            if !(s.ends_with(".gz") || s.ends_with(".bgz") || s.ends_with(".bam")) {
                return Err(format!("index=True needs a .gz (BGZF) or .bam output, and the output is {output}"));
            }
        }
        let op = orcamento::opcoes(memory_mb, temp_dir, ordem_ou_padrao(chrom_order, formato)?);
        let r = externo::ordenar_arquivos(&paths, output, formato, criterio, &op)?;
        if index {
            indice::indexar(output, Some(formato))?;
        }
        Ok((r.registros, r.lotes, r.ignoradas))
    })
    .map_err(erro_py)
}

/// Merges the intervals of one or more BED files: overlapping, touching, or
/// up to `distance` bases apart (like `bedtools merge -d`). The input does not
/// need to be sorted or to fit in memory. Writes `chromosome start end` and
/// returns `(intervals_read, intervals_written)`.
///
/// `paths` is a list of files. `chrom_order`, `memory_mb` and `temp_dir` work
/// as in `sort_genomic`.
#[pyfunction]
#[pyo3(signature = (paths, output, distance=0, chrom_order=None, memory_mb=None, temp_dir=None))]
fn bed_merge(
    py: Python<'_>,
    paths: Vec<String>,
    output: &str,
    distance: u64,
    chrom_order: Option<&str>,
    memory_mb: Option<usize>,
    temp_dir: Option<&str>,
) -> PyResult<(usize, usize)> {
    py.allow_threads(|| {
        let op = orcamento::opcoes(memory_mb, temp_dir, ordem_ou_padrao(chrom_order, externo::Formato::Bed)?);
        intervalos::bed_fundir(&paths, output, distance, &op)
    })
    .map_err(erro_py)
}

/// Intersects the intervals of A with those of B (like `bedtools intersect`).
/// `a` and `b` are lists of BED files.
///
/// `mode`: "regions" (the line of A clipped to each shared stretch, the
/// bedtools default), "unique" (`-u`: lines of A that touch B) or "absent"
/// (`-v`: lines of A that do not touch B). The inputs do not need to be sorted
/// or to fit in memory. Returns the number of lines written.
///
/// Raises an error, instead of writing an empty file, when A and B have no
/// sequence name in common (`chr1` against `1`, or different reference
/// genomes).
#[pyfunction]
#[pyo3(signature = (a, b, output, mode="regions", chrom_order=None, memory_mb=None, temp_dir=None))]
fn bed_intersect(
    py: Python<'_>,
    a: Vec<String>,
    b: Vec<String>,
    output: &str,
    mode: &str,
    chrom_order: Option<&str>,
    memory_mb: Option<usize>,
    temp_dir: Option<&str>,
) -> PyResult<usize> {
    py.allow_threads(|| {
        let modo = intervalos::ModoIntersecao::de_texto(mode)
            .ok_or_else(|| format!("Unknown mode: {mode} (use \"regions\", \"unique\" or \"absent\")"))?;
        let op = orcamento::opcoes(memory_mb, temp_dir, ordem_ou_padrao(chrom_order, externo::Formato::Bed)?);
        intervalos::bed_intersecao(&a, &b, output, modo, &op)
    })
    .map_err(erro_py)
}

// ---------------------------------------------------------
// VALIDACAO, NOMES DE CROMOSSOMO E GENOMA DE REFERENCIA
// ---------------------------------------------------------

/// Converte um texto JSON produzido no lado Rust num objeto Python.
fn json_para_py(py: Python<'_>, texto: &str) -> PyResult<Py<PyAny>> {
    Ok(py.import("json")?.call_method1("loads", (texto,))?.unbind())
}

/// Validates a VCF, BED, GFF/GTF, SAM, FASTA or FASTQ file of any size and
/// says in plain language what is wrong with it: each problem carries the LINE
/// of its first occurrence, the message, how many times it repeats and, when
/// there is one, how to fix it.
///
/// - `format`: "vcf", "bed", "gff", "sam", "fasta" or "fastq". Without it,
///   taken from the extension.
/// - `max_lines`: validate only the beginning of the file (a quick check).
///
/// Returns a dictionary: `valid`, `errors`, `warnings`, `lines`, `records`,
/// `complete`, `sorted`, `naming` ("chr", "no_chr", "mixed", "other" or "-"),
/// `sequences`, `path`, `format`, `problems` (a list of dictionaries with
/// `code`, `severity`, `line`, `message`, `occurrences` and `fix`) and `text`
/// (the report, ready to print).
///
/// Only raises an exception when the file cannot be opened; problems with the
/// CONTENT go in the report.
#[pyfunction]
#[pyo3(signature = (path, format=None, max_lines=None))]
fn validate_file(py: Python<'_>, path: &str, format: Option<&str>, max_lines: Option<usize>) -> PyResult<Py<PyAny>> {
    let relatorio = py
        .allow_threads(|| {
            let tipo = match format {
                Some(f) => Some(validador::Tipo::de_texto(f).ok_or_else(|| format!("Unknown format: {f}"))?),
                None if validador::Tipo::pelo_nome(path).is_none() => externo::Formato::do_arquivo(path).map(|f| match f {
                    externo::Formato::Vcf => validador::Tipo::Vcf,
                    externo::Formato::Bed => validador::Tipo::Bed,
                    externo::Formato::Gff => validador::Tipo::Gff,
                    externo::Formato::Sam => validador::Tipo::Sam,
                }),
                None => None,
            };
            validador::validar_arquivo(path, tipo, max_lines)
        })
        .map_err(erro_py)?;
    let obj = json_para_py(py, &relatorio.json())?;
    obj.bind(py).set_item("text", relatorio.resumo())?;
    Ok(obj)
}

/// Rewrites a file with its chromosome names in a single convention:
/// `style="chr"` (chr1, chrX, chrM -- UCSC) or `"no_chr"` (1, X, MT --
/// Ensembl/NCBI). Changes the records and the header. Streaming: works for a
/// file of any size.
///
/// `mapping` is a two-column file (original name, new name) for what cannot
/// be derived from the prefix: scaffolds and RefSeq accessions.
///
/// Returns a dictionary: `lines`, `records`, `changed`, `converted` (pairs)
/// and `kept` (names that stayed as they were).
#[pyfunction]
#[pyo3(signature = (input, output, style="chr", format=None, mapping=None))]
fn normalize_chromosomes(
    py: Python<'_>,
    input: &str,
    output: &str,
    style: &str,
    format: Option<&str>,
    mapping: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let resumo = py
        .allow_threads(|| {
            let estilo = nomes::Estilo::de_texto(style)
                .ok_or_else(|| format!("Unknown style: {style} (use \"chr\" or \"no_chr\")"))?;
            let formato = match format {
                Some(f) => externo::Formato::de_texto(f).ok_or_else(|| format!("Unknown format: {f}"))?,
                None => externo::Formato::do_arquivo(input)
                    .ok_or("Could not infer the format from the file extension; pass `format`")?,
            };
            nomes::normalizar_cromossomos(input, output, formato, estilo, mapping)
        })
        .map_err(erro_py)?;
    json_para_py(py, &resumo.json())
}

/// What a file declares about its reference genome. Reads VCF, SAM, GFF3,
/// `.fai` and `chrom.sizes`.
///
/// Returns a dictionary: `assembly` (the human assembly when recognized:
/// GRCh37/hg19 or GRCh38/hg38), `declared` (the text the file declares) and
/// `sequences` (pairs of name and length).
#[pyfunction]
fn file_reference(py: Python<'_>, path: &str) -> PyResult<Py<PyAny>> {
    let r = py.allow_threads(|| nomes::referencia_do_arquivo(path)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// Are two files on the same reference genome? Compares the LENGTHS of the
/// sequences each one declares, matching `chr1` with `1`.
///
/// Returns a dictionary with `compatible` (True, False, or None when the
/// files do not say enough), `verdict` (the conclusion in one sentence),
/// `matching`, `different`, `only_in_a`, `only_in_b`, `different_naming`,
/// `assembly_a`, `assembly_b`, `declared_a` and `declared_b`.
#[pyfunction]
fn compare_references(py: Python<'_>, a: &str, b: &str) -> PyResult<Py<PyAny>> {
    let c = py.allow_threads(|| nomes::comparar_referencias(a, b)).map_err(erro_py)?;
    json_para_py(py, &c.json())
}

/// What this installation has. Returns a dictionary: `version`,
/// `gpu_compiled` (the extension carries the CUDA kernel), `gpu` (there is a
/// usable NVIDIA card right now), `gpu_memory_mb` (`[free, total]`, or None)
/// and `threads` (cores the CPU paths will use).
#[pyfunction]
fn engine_info(py: Python<'_>) -> PyResult<Py<PyAny>> {
    let memoria = match gpu::memoria() {
        Some((livre, total)) => format!("[{},{}]", livre >> 20, total >> 20),
        None => "null".to_string(),
    };
    let texto = format!(
        "{{\"version\":\"{}\",\"gpu_compiled\":{},\"gpu\":{},\"gpu_memory_mb\":{},\"threads\":{}}}",
        env!("CARGO_PKG_VERSION"),
        gpu::COMPILADA,
        gpu::disponivel(),
        memoria,
        rayon::current_num_threads()
    );
    json_para_py(py, &texto)
}

// ---------------------------------------------------------
// MAIS OPERACOES DE INTERVALOS (BED)
// ---------------------------------------------------------

/// Removes from each interval of A the stretches covered by B (like
/// `bedtools subtract`) and writes what is left: the line of A with the
/// coordinates of each remaining piece. `a` and `b` are lists of BED files.
///
/// `remove_whole=True` (`bedtools subtract -A`) drops every interval of A
/// that B touches at all, instead of clipping it.
///
/// The inputs do not need to be sorted or to fit in memory. Returns the number
/// of lines written. Raises an error, instead of writing A unchanged, when A
/// and B have no sequence name in common.
#[pyfunction]
#[pyo3(signature = (a, b, output, remove_whole=false, chrom_order=None, memory_mb=None, temp_dir=None))]
#[allow(clippy::too_many_arguments)]
fn bed_subtract(
    py: Python<'_>,
    a: Vec<String>,
    b: Vec<String>,
    output: &str,
    remove_whole: bool,
    chrom_order: Option<&str>,
    memory_mb: Option<usize>,
    temp_dir: Option<&str>,
) -> PyResult<usize> {
    py.allow_threads(|| {
        let op = orcamento::opcoes(memory_mb, temp_dir, ordem_ou_padrao(chrom_order, externo::Formato::Bed)?);
        bed_ops::bed_subtract(&a, &b, output, remove_whole, &op)
    })
    .map_err(erro_py)
}

/// Writes the stretches of the genome that no interval covers (like
/// `bedtools complement`). `paths` is a list of BED files; `genome` is a
/// `.fai` index or a chromosome-sizes file (name and length, tab-separated),
/// and its order is the order of the output.
///
/// Returns the number of lines written. A sequence that is in the BED files
/// and not in the genome file is an error.
#[pyfunction]
#[pyo3(signature = (paths, genome, output, memory_mb=None, temp_dir=None))]
fn bed_complement(
    py: Python<'_>,
    paths: Vec<String>,
    genome: &str,
    output: &str,
    memory_mb: Option<usize>,
    temp_dir: Option<&str>,
) -> PyResult<usize> {
    py.allow_threads(|| {
        let op = orcamento::opcoes(memory_mb, temp_dir, externo::OrdemCromossomos::Lexica);
        bed_ops::bed_complement(&paths, genome, output, &op)
    })
    .map_err(erro_py)
}

/// For each interval of A, finds the nearest interval(s) of B on the same
/// sequence (like `bedtools closest -d`). Writes one line per closest B: the
/// line of A, the line of B and the distance. `a` and `b` are lists of BED
/// files.
///
/// The distance is 0 when the two overlap; otherwise the number of bases
/// between them plus 1, so two intervals that touch are at distance 1. Every
/// B tied at the smallest distance is written. When B has nothing on the
/// sequence, the line of A is followed by `.`, `-1`, `-1` and the distance
/// `-1`.
///
/// Returns the number of lines written. Raises an error when A and B have no
/// sequence name in common.
#[pyfunction]
#[pyo3(signature = (a, b, output, chrom_order=None, memory_mb=None, temp_dir=None))]
fn bed_closest(
    py: Python<'_>,
    a: Vec<String>,
    b: Vec<String>,
    output: &str,
    chrom_order: Option<&str>,
    memory_mb: Option<usize>,
    temp_dir: Option<&str>,
) -> PyResult<usize> {
    py.allow_threads(|| {
        let op = orcamento::opcoes(memory_mb, temp_dir, ordem_ou_padrao(chrom_order, externo::Formato::Bed)?);
        bed_ops::bed_closest(&a, &b, output, &op)
    })
    .map_err(erro_py)
}

/// For each interval of A, how much of it B covers (like `bedtools
/// coverage`). Writes the line of A plus four columns: the number of
/// intervals of B that overlap it, the number of its bases covered by at
/// least one of them, its length, and the covered fraction. `a` and `b` are
/// lists of BED files.
///
/// Returns the number of lines written. Raises an error when A and B have no
/// sequence name in common.
#[pyfunction]
#[pyo3(signature = (a, b, output, chrom_order=None, memory_mb=None, temp_dir=None))]
fn bed_coverage(
    py: Python<'_>,
    a: Vec<String>,
    b: Vec<String>,
    output: &str,
    chrom_order: Option<&str>,
    memory_mb: Option<usize>,
    temp_dir: Option<&str>,
) -> PyResult<usize> {
    py.allow_threads(|| {
        let op = orcamento::opcoes(memory_mb, temp_dir, ordem_ou_padrao(chrom_order, externo::Formato::Bed)?);
        bed_ops::bed_coverage(&a, &b, output, &op)
    })
    .map_err(erro_py)
}

/// Grows every interval (like `bedtools slop`): the start goes down by `left`
/// bases and the end goes up by `right`, never below 0 nor beyond the length
/// of the sequence in `genome` (a `.fai` index or a chromosome-sizes file).
/// `both` sets the two sides at once. The lines keep their order and their
/// other columns. Returns the number of lines written.
#[pyfunction]
#[pyo3(signature = (paths, genome, output, both=None, left=0, right=0))]
fn bed_slop(
    py: Python<'_>,
    paths: Vec<String>,
    genome: &str,
    output: &str,
    both: Option<u64>,
    left: u64,
    right: u64,
) -> PyResult<usize> {
    py.allow_threads(|| {
        let (l, r) = (both.unwrap_or(left), both.unwrap_or(right));
        let op = orcamento::opcoes(Some(1), None, externo::OrdemCromossomos::Lexica);
        bed_ops::bed_slop(&paths, genome, output, l, r, &op)
    })
    .map_err(erro_py)
}

// ---------------------------------------------------------
// VCF: ESTATISTICAS, FILTRO E NORMALIZACAO
// ---------------------------------------------------------

/// Aceita um texto ou uma lista de textos.
fn lista_de_textos(valor: Option<&Bound<'_, PyAny>>, nome: &str) -> PyResult<Vec<String>> {
    let Some(v) = valor else {
        return Ok(Vec::new());
    };
    if let Ok(s) = v.extract::<String>() {
        return Ok(vec![s]);
    }
    v.extract::<Vec<String>>()
        .map_err(|_| pyo3::exceptions::PyTypeError::new_err(format!("`{nome}` must be a string or a list of strings")))
}

/// Counts what a VCF holds, in one pass over a file of any size.
///
/// Returns a dictionary: `records`, `samples`, `sequences` (pairs of name and
/// records), `snps`, `mnps`, `insertions`, `deletions`, `other` (symbolic
/// alleles, breakends, `*`), `no_alt`, `multiallelic`, `transitions`,
/// `transversions`, `ts_tv`, `pass`, `filter_missing`, `qual_missing`,
/// `mean_qual`, `skipped` (lines that are not valid records) and `text` (the
/// report, ready to print). A site counts once for each type present among
/// its ALT alleles.
#[pyfunction]
fn vcf_stats(py: Python<'_>, path: &str) -> PyResult<Py<PyAny>> {
    let r = py.allow_threads(|| vcf_tools::vcf_stats(path)).map_err(erro_py)?;
    let obj = json_para_py(py, &r.json())?;
    obj.bind(py).set_item("text", r.resumo())?;
    Ok(obj)
}

/// Keeps the records of a VCF that satisfy an expression, fall in a region,
/// or both, and optionally keeps only some samples. Streaming: works for a
/// file of any size, and needs no index.
///
/// - `expression`: for example `QUAL >= 30 && INFO/DP > 10`,
///   `TYPE == "snp" && FILTER == "PASS"`, `!(INFO/AF < 0.01)`. Fields: `QUAL`,
///   `POS`, `N_ALT` (numbers); `CHROM`, `ID`, `REF`, `ALT`, `FILTER`, `TYPE`
///   (text; TYPE is "snp", "mnp", "ins", "del", "indel", "other" or "ref");
///   `INFO/KEY` (alone, it tests that the key is present). Operators: `==`,
///   `!=`, `<`, `<=`, `>`, `>=`, `~`, `!~`, `&&`, `||`, `!` and parentheses.
///   A comparison with a missing value is false.
///   TYPE works as in bcftools: `TYPE == "snp"` keeps a site when EVERY ALT
///   allele is a SNP (the default reading), `TYPE ~ "snp"` when AT LEAST ONE
///   is; `!=` and `!~` are their negations. `FILTER ~ "q10"` keeps a site
///   whose FILTER includes q10 (`FILTER == "q10"` asks for exactly q10).
/// - `regions`: a string or a list such as `"chr1"`, `"chr1:1000-2000"`
///   (1-based, end included). A record is kept when it overlaps any of them.
/// - `samples`: names of the sample columns to keep, in this order. An empty
///   list keeps the sites only.
///
/// Returns a dictionary: `read`, `kept`, `samples_before`, `samples_after`
/// and `skipped`. A mistake in the expression raises an error that points at
/// the character where it stopped making sense.
#[pyfunction]
#[pyo3(signature = (input, output, expression=None, regions=None, samples=None))]
fn vcf_filter(
    py: Python<'_>,
    input: &str,
    output: &str,
    expression: Option<&str>,
    regions: Option<&Bound<'_, PyAny>>,
    samples: Option<&Bound<'_, PyAny>>,
) -> PyResult<Py<PyAny>> {
    let regioes = lista_de_textos(regions, "regions")?;
    let amostras = match samples {
        Some(s) => Some(lista_de_textos(Some(s), "samples")?),
        None => None,
    };
    let r = py
        .allow_threads(|| vcf_tools::vcf_filter(input, output, expression, &regioes, amostras.as_deref()))
        .map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// Normalizes a VCF: with `split_multiallelic`, a site with several ALT
/// alleles becomes one record per allele (INFO and FORMAT fields declared
/// with `Number=A`, `R` or `G` are split with it); then every variant is
/// trimmed to its shortest representation (shared bases at the end and at the
/// start are removed, and POS moves with them).
///
/// - `reference`: a FASTA of the genome (plain, `.gz` or `.lz4`; a `.fai`
///   from `samtools faidx` next to a plain one saves a pass over it). With
///   it, REF is checked against the genome and every indel is LEFT-ALIGNED:
///   inside a repeat it moves to the leftmost position that gives the same
///   sequence, as `bcftools norm -f` does, so two files that wrote the same
///   indel at different places of a repeat agree. Records whose REF differs
///   from the genome, or whose sequence the genome lacks (`chr1` and `1` are
///   matched), are only trimmed, and counted. Without it (the default), there
///   is no left-alignment.
///
/// If positions change, the output is sorted again, whatever its size. In a
/// split record, a genotype allele that belongs to another record is written
/// as missing (`1/2` becomes `1/.` and `./1`).
///
/// Returns a dictionary: `read`, `written`, `split`, `trimmed`,
/// `position_changed`, `left_aligned`, `ref_mismatch`, `not_in_reference`,
/// `resorted`, `skipped` and `runs`.
#[pyfunction]
#[pyo3(signature = (input, output, split_multiallelic=true, reference=None, chrom_order=None, memory_mb=None, temp_dir=None))]
#[allow(clippy::too_many_arguments)]
fn vcf_normalize(
    py: Python<'_>,
    input: &str,
    output: &str,
    split_multiallelic: bool,
    reference: Option<&str>,
    chrom_order: Option<&str>,
    memory_mb: Option<usize>,
    temp_dir: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let r = py
        .allow_threads(|| {
            let op = orcamento::opcoes(memory_mb, temp_dir, ordem_ou_padrao(chrom_order, externo::Formato::Vcf)?);
            vcf_tools::vcf_normalize(input, output, split_multiallelic, reference, &op)
        })
        .map_err(erro_py)?;
    json_para_py(py, &r.json())
}

// ---------------------------------------------------------
// FASTQ: QUALIDADE, RELATORIO E APARO
// ---------------------------------------------------------

/// Measures the quality of a FASTQ file of any size in one pass: read
/// lengths, quality at each position of the read, base composition, GC
/// content, N bases, duplicated reads and known adapters.
///
/// `max_reads` looks only at the beginning of the file (a quick check).
///
/// Returns a dictionary: `path`, `reads`, `bases`, `complete`, `min_length`,
/// `max_length`, `mean_length`, `length_histogram`, `phred_offset` (33 or
/// 64, detected), `mean_quality`, `read_quality_histogram`, `per_position`
/// (columns `bases`, `mean_quality`, `q10`, `q25`, `q50`, `q75`, `q90`, `a`,
/// `c`, `g`, `t`, `n`; index 0 is the first base of the read),
/// `positions_reported`, `gc_percent`, `gc_histogram`, `n_bases`, `n_percent`,
/// `reads_with_n`, `duplication_percent`, `duplication_sample`, `adapters`,
/// `likely_adapter`, `warnings` (plain-language notes, each with the fix)
/// and `text` (the report, ready to print). For the same numbers as a page
/// with charts, use `fastq_report`.
#[pyfunction]
#[pyo3(signature = (path, max_reads=None))]
fn fastq_stats(py: Python<'_>, path: &str, max_reads: Option<u64>) -> PyResult<Py<PyAny>> {
    let r = py.allow_threads(|| fastq_qc::fastq_stats(path, max_reads)).map_err(erro_py)?;
    let obj = json_para_py(py, &r.json())?;
    obj.bind(py).set_item("text", r.resumo())?;
    Ok(obj)
}

/// Writes the quality report of a FASTQ file as one self-contained HTML page
/// (`output`), with charts: quality and base composition along the read,
/// quality and GC content per read, read lengths and adapters. The page opens
/// in any browser, offline. Returns the same dictionary as `fastq_stats`.
#[pyfunction]
#[pyo3(signature = (path, output, max_reads=None))]
fn fastq_report(py: Python<'_>, path: &str, output: &str, max_reads: Option<u64>) -> PyResult<Py<PyAny>> {
    let r = py
        .allow_threads(|| {
            let r = fastq_qc::fastq_stats(path, max_reads)?;
            relatorio::grava_fastq_html(&r, output)?;
            Ok::<_, String>(r)
        })
        .map_err(erro_py)?;
    let obj = json_para_py(py, &r.json())?;
    obj.bind(py).set_item("text", r.resumo())?;
    Ok(obj)
}

/// Trims and filters the reads of a FASTQ file, single-end or paired-end,
/// in the order cutadapt uses: low-quality ends, then the 3' adapter, then N
/// at the ends, then the filters.
///
/// - `input`, `output`: the reads and where to write them (`.gz` and `.lz4`
///   by extension). For paired-end data also pass `input2` and `output2`: a
///   pair is kept or dropped as a whole, and the two outputs stay in step.
/// - `adapter`: the 3' adapter of read 1, or `"auto"` to use the known
///   adapter found in the first reads of the file. `adapter2` is the adapter
///   of read 2; without it, read 2 uses `adapter`.
/// - `min_overlap`, `max_error_rate`: how much of the adapter must match, and
///   the fraction of mismatches allowed (no insertions or deletions).
/// - `quality_cutoff` (3' end) and `quality_cutoff_5` (5' end): Phred cutoffs
///   for the quality trimming; 0 turns it off.
/// - `trim_n`: remove N bases from both ends.
/// - `min_length`, `max_n`, `min_mean_quality`: reads shorter than, with more
///   N than, or with a lower mean quality than these are discarded.
/// - `phred_offset`: 33 (the default) or 64.
///
/// Returns a dictionary: `paired`, `reads_in`, `reads_out`, `bases_in`,
/// `bases_out`, `reads_with_adapter`, `reads_with_adapter2`,
/// `bases_removed_quality`, `bases_removed_adapter`, `bases_removed_n`,
/// `discarded_too_short`, `discarded_too_many_n`, `discarded_low_quality`,
/// and `adapter` / `adapter2` (the sequences that were actually used).
#[pyfunction]
#[pyo3(signature = (input, output, input2=None, output2=None, adapter=None, adapter2=None, min_overlap=3,
                    max_error_rate=0.1, quality_cutoff=0, quality_cutoff_5=0, trim_n=false, min_length=0,
                    max_n=None, min_mean_quality=None, phred_offset=33))]
#[allow(clippy::too_many_arguments)]
fn fastq_trim(
    py: Python<'_>,
    input: &str,
    output: &str,
    input2: Option<&str>,
    output2: Option<&str>,
    adapter: Option<&str>,
    adapter2: Option<&str>,
    min_overlap: usize,
    max_error_rate: f64,
    quality_cutoff: u8,
    quality_cutoff_5: u8,
    trim_n: bool,
    min_length: usize,
    max_n: Option<usize>,
    min_mean_quality: Option<u8>,
    phred_offset: u8,
) -> PyResult<Py<PyAny>> {
    let (r, a1, a2) = py
        .allow_threads(|| {
            // "auto": o adaptador conhecido mais frequente no comeco do arquivo.
            let detecta = |arquivo: &str| -> Result<Option<String>, String> {
                Ok(fastq_qc::fastq_stats(arquivo, Some(200_000))?.likely_adapter.map(|a| a.1.to_string()))
            };
            let e_auto = |a: Option<&str>| a.is_some_and(|t| t.eq_ignore_ascii_case("auto"));
            let mut a1 = adapter.map(str::to_string);
            let mut a2 = adapter2.map(str::to_string);
            if e_auto(adapter) {
                a1 = detecta(input)?;
                if let (None, Some(segundo)) = (adapter2, input2) {
                    a2 = detecta(segundo)?;
                }
            }
            if e_auto(adapter2) {
                a2 = match input2 {
                    Some(segundo) => detecta(segundo)?,
                    None => return Err("adapter2=\"auto\" needs input2: there is no second file to look at".to_string()),
                };
            }
            let opt = fastq_qc::TrimOptions {
                adapter: a1.clone(),
                adapter2: a2.clone(),
                min_overlap,
                max_error_rate,
                quality_cutoff,
                quality_cutoff_5,
                trim_n,
                min_length,
                max_n,
                min_mean_quality,
                phred_offset,
            };
            let r = fastq_qc::fastq_trim(input, output, input2, output2, &opt)?;
            // Sem adaptador proprio, a leitura 2 usou o da leitura 1.
            let usado2 = if input2.is_some() { a2.or_else(|| a1.clone()) } else { None };
            Ok((r, a1, usado2))
        })
        .map_err(erro_py)?;
    let obj = json_para_py(py, &r.json())?;
    obj.bind(py).set_item("adapter", a1)?;
    obj.bind(py).set_item("adapter2", a2)?;
    Ok(obj)
}

// ---------------------------------------------------------
// CONTAGEM DE K-MERS DE QUALQUER TAMANHO
// ---------------------------------------------------------

/// Counts the k-mers of FASTA and FASTQ files of any size: what does not fit
/// the memory budget is sorted in runs on disk and merged.
///
/// - `paths`: a list of FASTA and/or FASTQ files (plain, `.gz` or `.lz4`).
/// - `output`: a table with one line per distinct k-mer, `KMER<TAB>count`, in
///   alphabetical order (the extension picks the compression). Without it,
///   only the summary is computed.
/// - `k`: 1 to 63.
/// - `canonical`: count a k-mer and its reverse complement as one (the
///   smaller of the two is written). This is what makes reads from the two
///   strands add up; turn it off only for single-strand data.
/// - `min_count`: leave rarer k-mers out of the table (they still count in
///   the summary).
/// - `memory_mb`, `temp_dir`: as in `sort_genomic`.
///
/// A k-mer never crosses from one sequence to the next, and never contains a
/// base other than A, C, G or T.
///
/// Returns a dictionary: `k`, `canonical`, `sequences`, `bases`,
/// `total_kmers`, `distinct_kmers`, `unique_kmers` (seen once), `max_count`,
/// `written`, `runs` (run files on disk; 0 when everything fit in memory),
/// `histogram` (pairs of count and number of distinct k-mers with it: the
/// abundance spectrum; counts above 10000 share the last entry) and `text`
/// (the summary, ready to print).
#[pyfunction]
#[pyo3(signature = (paths, output=None, k=31, canonical=true, min_count=1, memory_mb=None, temp_dir=None))]
#[allow(clippy::too_many_arguments)]
fn count_kmers(
    py: Python<'_>,
    paths: Vec<String>,
    output: Option<&str>,
    k: u8,
    canonical: bool,
    min_count: u32,
    memory_mb: Option<usize>,
    temp_dir: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let r = py
        .allow_threads(|| {
            let op = orcamento::opcoes(memory_mb, temp_dir, externo::OrdemCromossomos::Lexica);
            let opt = kmer_count::CountOptions { k, canonical, min_count, memoria_bytes: op.memoria_bytes, dir_temp: op.dir_temp };
            kmer_count::count_kmers(&paths, output, &opt)
        })
        .map_err(erro_py)?;
    let obj = json_para_py(py, &r.json())?;
    obj.bind(py).set_item("text", r.resumo())?;
    Ok(obj)
}

// ---------------------------------------------------------
// ENVELOPES: CONVERSAO, INDICE E CONSULTA POR REGIAO
// ---------------------------------------------------------

/// Rewrites a file in another container, keeping every line and its order:
/// plain text, `.gz` (BGZF: ordinary gzip that can also be indexed), `.lz4`,
/// `.bam` (from and to SAM) and `.parquet` (typed columns, with the header
/// kept in the file metadata). The extensions of `input` and `output` say
/// what to do: `convert_file("x.sam", "x.bam")`,
/// `convert_file("x.vcf.gz", "x.vcf.parquet")`, `convert_file("x.bam",
/// "x.sam")`. Streaming: works for a file of any size.
///
/// `format` ("vcf", "bed", "gff" or "sam") is only needed when neither name
/// says it and the output is Parquet.
///
/// Returns a dictionary: `lines` and `records` (lines that are not header).
#[pyfunction]
#[pyo3(signature = (input, output, format=None))]
fn convert_file(py: Python<'_>, input: &str, output: &str, format: Option<&str>) -> PyResult<Py<PyAny>> {
    let r = py
        .allow_threads(|| {
            let formato = match format {
                Some(f) => Some(externo::Formato::de_texto(f).ok_or_else(|| format!("Unknown format: {f}"))?),
                None => None,
            };
            externo::converter(input, output, formato)
        })
        .map_err(erro_py)?;
    json_para_py(py, &format!("{{\"lines\":{},\"records\":{}}}", r.linhas, r.registros))
}

/// How every function writes `.parquet` from now on, in this Python process.
/// Called without arguments, it only returns the current settings.
///
/// - `numbers`: how QUAL (VCF) and score (GFF) are stored.
///   `"number"` (the default) stores them as numbers, so tables can be
///   filtered and summed without conversion (`pl.col("QUAL") > 30`); a
///   missing value (`.`) is stored as null, and on the way back to text the
///   number is written in its shortest form: `12.50` comes back as `12.5`,
///   `30.0` as `30`. `"text"` stores them exactly as written, so converting
///   back gives the original file byte for byte.
///
/// Returns a dictionary with the settings now in force: `numbers`.
#[pyfunction]
#[pyo3(signature = (numbers=None))]
fn parquet_options(py: Python<'_>, numbers: Option<&str>) -> PyResult<Py<PyAny>> {
    #[cfg(feature = "parquet")]
    {
        if let Some(n) = numbers {
            match n {
                "number" => parquet_io::define_numeros_como_texto(false),
                "text" => parquet_io::define_numeros_como_texto(true),
                _ => {
                    return Err(pyo3::exceptions::PyValueError::new_err(format!(
                        "numbers must be \"number\" (QUAL and score stored as numbers, the default) or \"text\" (stored exactly as written), not \"{n}\""
                    )))
                }
            }
        }
        let atual = if parquet_io::numeros_como_texto() { "text" } else { "number" };
        json_para_py(py, &format!("{{\"numbers\":\"{atual}\"}}"))
    }
    #[cfg(not(feature = "parquet"))]
    {
        let _ = (py, numbers);
        Err(erro_py("Parquet support is not compiled into this build".to_string()))
    }
}

/// Builds the index that lets `query_region` (and tabix, samtools, bcftools,
/// IGV) jump straight to a region: a `.tbi` next to a BGZF-compressed VCF,
/// BED or GFF (`x.vcf.gz.tbi`), or a `.bai` next to a BAM (`x.bam.bai`).
///
/// The file must be sorted by position (`sort_genomic`, or
/// `sort_genomic(..., index=True)` to do both at once) and, for text, be a
/// `.gz` written as BGZF: every `.gz` this library writes is.
///
/// Returns a dictionary: `index` (the path written), `kind` ("tbi" or
/// "bai"), `sequences`, `records` and `skipped`.
#[pyfunction]
#[pyo3(signature = (path, format=None))]
fn index_file(py: Python<'_>, path: &str, format: Option<&str>) -> PyResult<Py<PyAny>> {
    let r = py
        .allow_threads(|| {
            let formato = match format {
                Some(f) => Some(externo::Formato::de_texto(f).ok_or_else(|| format!("Unknown format: {f}"))?),
                None => None,
            };
            indice::indexar(path, formato)
        })
        .map_err(erro_py)?;
    json_para_py(
        py,
        &format!(
            "{{\"index\":{},\"kind\":\"{}\",\"sequences\":{},\"records\":{},\"skipped\":{}}}",
            nomes::json_texto(&r.indice),
            r.tipo,
            r.sequencias,
            r.registros,
            r.ignoradas
        ),
    )
}

/// Returns the records of an indexed file that overlap one or more regions,
/// reading only the parts of the file the index points to. Works on a
/// BGZF-compressed VCF, BED or GFF with a `.tbi`, and on a BAM with a `.bai`
/// (see `index_file`).
///
/// - `regions`: a string or a list: `"chr1"`, `"chr1:1000"`,
///   `"chr1:1000-2000"` (1-based, end included, as in samtools and tabix).
///   Overlapping regions are merged: each record comes out once, in file
///   order.
/// - `output`: where to write the records (any container, by extension). The
///   header goes first unless `header=False`.
///
/// With `output`, returns the number of records written. Without it, returns
/// the records as a list of strings (a BAM gives SAM lines).
#[pyfunction]
#[pyo3(signature = (path, regions, output=None, header=true))]
fn query_region(
    py: Python<'_>,
    path: &str,
    regions: &Bound<'_, PyAny>,
    output: Option<&str>,
    header: bool,
) -> PyResult<Py<PyAny>> {
    let regioes = lista_de_textos(Some(regions), "regions")?;
    match output {
        Some(saida) => {
            let n = py
                .allow_threads(|| {
                    let formato = externo::Formato::do_arquivo(path);
                    let mut w = match formato {
                        Some(f) => externo::Escritor::cria_para(saida, f)?,
                        None => externo::Escritor::cria(saida)?,
                    };
                    // BAM e Parquet nao existem sem o cabecalho.
                    let s = saida.to_ascii_lowercase();
                    let com_cabecalho = header || s.ends_with(".bam") || s.ends_with(".parquet");
                    let r = indice::consultar(path, &regioes, &mut |cabecalho, linha| {
                        if cabecalho && !com_cabecalho {
                            Ok(())
                        } else {
                            w.linha(linha)
                        }
                    });
                    let fim = match r {
                        Ok(n) => w.termina().map(|_| n),
                        Err(e) => {
                            drop(w);
                            Err(e)
                        }
                    };
                    if fim.is_err() {
                        let _ = std::fs::remove_file(saida);
                    }
                    fim
                })
                .map_err(erro_py)?;
            Ok(n.into_pyobject(py)?.into_any().unbind())
        }
        None => {
            let linhas = py
                .allow_threads(|| {
                    let mut linhas: Vec<String> = Vec::new();
                    indice::consultar(path, &regioes, &mut |cabecalho, linha| {
                        if !cabecalho {
                            linhas.push(linha.to_string());
                        }
                        Ok(())
                    })?;
                    Ok::<_, String>(linhas)
                })
                .map_err(erro_py)?;
            Ok(linhas.into_pyobject(py)?.into_any().unbind())
        }
    }
}

// ---------------------------------------------------------
// Genoma humano (passo 15). Experimental, para ensino e pesquisa.
// ---------------------------------------------------------

const AVISO: &str = "Experimental software, for teaching and research only: not a medical diagnosis, not a paternity test \
with legal value and not advice on any treatment or medicine. The authors take no responsibility for decisions made on \
the basis of these results; a finding must be confirmed by an accredited laboratory and interpreted by a qualified professional.";

fn opcoes_de_lista(valor: Option<&Bound<'_, PyAny>>, nome: &str, n: usize) -> PyResult<Vec<Option<String>>> {
    let Some(v) = valor else {
        return Ok(vec![None; n]);
    };
    if v.is_none() {
        return Ok(vec![None; n]);
    }
    if let Ok(s) = v.extract::<String>() {
        return Ok(vec![Some(s); n]);
    }
    let lista: Vec<Option<String>> = v
        .extract()
        .map_err(|_| pyo3::exceptions::PyTypeError::new_err(format!("`{nome}` must be a string, a list or None")))?;
    if lista.len() != n {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "`{nome}` has {} item(s), one per file is needed ({n})",
            lista.len()
        )));
    }
    Ok(lista)
}

/// The ClinVar findings of one person: the variants of a VCF that ClinVar
/// lists, with their classification, review stars, conditions and zygosity.
/// Experimental, for teaching and research only: not a medical diagnosis.
///
/// - `clinvar`: the ClinVar VCF (GRCh38 or GRCh37, as the person's file).
/// - `sample`: which sample of a multi-sample VCF (default: the first).
/// - `callable_bed`: the regions where the person's calls are reliable; only
///   variants inside them are compared (as in the Genome in a Bottle files).
/// - `classes`: which classifications to keep: "plp" (pathogenic and likely
///   pathogenic), "vus", "conflicting", "blb" (benign and likely benign),
///   "risk", "all", or ClinVar's own terms. Default: "plp".
/// - `min_stars`: minimum review stars (0-4). Default 1.
/// - `genes`, `regions`: keep only these genes, or a BED of regions (a panel).
/// - `output`: also write the findings as a TSV table.
///
/// Both files must be sorted in natural chromosome order (1..22, X, Y);
/// variants are split per allele and trimmed, not left-aligned: left-align
/// first with `vcf_normalize(reference=...)` if the VCF was not.
#[pyfunction]
#[pyo3(signature = (vcf, clinvar, output=None, sample=None, callable_bed=None, classes=None, min_stars=1, genes=None, regions=None))]
#[allow(clippy::too_many_arguments)]
fn clinvar_findings(
    py: Python<'_>,
    vcf: &str,
    clinvar: &str,
    output: Option<&str>,
    sample: Option<String>,
    callable_bed: Option<String>,
    classes: Option<&Bound<'_, PyAny>>,
    min_stars: u8,
    genes: Option<&Bound<'_, PyAny>>,
    regions: Option<String>,
) -> PyResult<Py<PyAny>> {
    let mut classes = lista_de_textos(classes, "classes")?;
    if classes.is_empty() {
        classes.push("plp".into());
    }
    let genes = lista_de_textos(genes, "genes")?;
    let pessoa = clinico::Person { vcf: vcf.to_string(), sample, callable: callable_bed };
    let op = clinico::FindingsOptions { classes, min_stars, genes, regions, output: output.map(|s| s.to_string()) };
    let r = py.allow_threads(|| clinico::clinvar_findings(&pessoa, clinvar, &op)).map_err(erro_py)?;
    let d = json_para_py(py, &r.json())?;
    d.bind(py).set_item("notice", AVISO)?;
    Ok(d)
}

/// The coordinates of named genes, from a GFF3/GTF annotation (for example
/// the RefSeq GFF, whose accessions NC_000001.11... are mapped to chromosome
/// names). Use them as a gene panel: `vcf_filter(regions=...)`,
/// `clinvar_findings(regions=...)`, `query_region`.
///
/// - `padding`: bases added on each side.
/// - `naming`: "chr" or "no_chr" for the chromosome names (default: as the GFF says).
/// - `output`: also write them as BED.
#[pyfunction]
#[pyo3(signature = (gff, genes, output=None, padding=0, naming=None))]
fn gene_regions(
    py: Python<'_>,
    gff: &str,
    genes: &Bound<'_, PyAny>,
    output: Option<&str>,
    padding: u64,
    naming: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let genes = lista_de_textos(Some(genes), "genes")?;
    let r = py.allow_threads(|| clinico::gene_regions(gff, &genes, padding, naming, output)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// The genotype of each person at each variant: testing relatives for a
/// variant already known in the family (cascade testing), as Sanger
/// sequencing would. Experimental, for teaching and research only.
///
/// - `variants`: "chrom:pos:ref:alt" (1-based), or, with `clinvar`, an rsID
///   or a ClinVar VariationID.
/// - `vcfs`: one or more VCF files; `samples` and `callable_beds` give, per
///   file, the sample and the callable regions (a string applies to all).
#[pyfunction]
#[pyo3(signature = (variants, vcfs, samples=None, callable_beds=None, clinvar=None))]
fn variant_lookup(
    py: Python<'_>,
    variants: &Bound<'_, PyAny>,
    vcfs: &Bound<'_, PyAny>,
    samples: Option<&Bound<'_, PyAny>>,
    callable_beds: Option<&Bound<'_, PyAny>>,
    clinvar: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let variantes = lista_de_textos(Some(variants), "variants")?;
    let arquivos = lista_de_textos(Some(vcfs), "vcfs")?;
    let amostras = opcoes_de_lista(samples, "samples", arquivos.len())?;
    let beds = opcoes_de_lista(callable_beds, "callable_beds", arquivos.len())?;
    let pessoas: Vec<clinico::Person> = arquivos
        .iter()
        .zip(amostras)
        .zip(beds)
        .map(|((v, s), b)| clinico::Person { vcf: v.clone(), sample: s, callable: b })
        .collect();
    let r = py.allow_threads(|| clinico::variant_lookup(&variantes, &pessoas, clinvar)).map_err(erro_py)?;
    let d = json_para_py(py, &r.json())?;
    d.bind(py).set_item("notice", AVISO)?;
    Ok(d)
}

/// The bases the reads of a BAM (or SAM) show at one position: depth, counts
/// of A, C, G, T, N, deletions and insertions, and the fraction of an allele —
/// what one looks at to confirm a variant. Uses the .bai when there is one.
///
/// - `position`: "chr1:12345" (1-based).
/// - `alt`: the allele whose fraction to report: "T", "C>T", "CAG>C",
///   "C>CTT", "+SEQ", "-N" or "*".
/// - `min_mapq`, `min_baseq`: quality filters (as samtools mpileup -q / -Q).
#[pyfunction]
#[pyo3(signature = (path, position, alt=None, min_mapq=0, min_baseq=13))]
fn pileup(py: Python<'_>, path: &str, position: &str, alt: Option<String>, min_mapq: u8, min_baseq: u8) -> PyResult<Py<PyAny>> {
    let (cromo, pos) = position
        .trim()
        .rsplit_once(':')
        .and_then(|(c, p)| p.replace(',', "").parse::<u64>().ok().map(|p| (c.to_string(), p)))
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err(format!("`{position}`: write the position as chrom:pos, e.g. chr1:12345")))?;
    let op = clinico::PileupOptions { min_mapq, min_baseq, alt, ..Default::default() };
    let r = py.allow_threads(|| clinico::pileup(path, &cromo, pos, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// Writes the input VCF of PharmCAT (the pharmacogenomics tool of PharmGKB
/// and CPIC, https://pharmcat.org) from a person's VCF: the genotype at each
/// position PharmCAT examines (`pharmcat_positions_<version>.vcf` of its
/// release). A position without a record becomes 0/0 only inside
/// `callable_bed` (or, without it, when `absent_is_ref`), and stays missing
/// on a sequence the file does not have. MMGEN calls no allele itself.
#[pyfunction]
#[pyo3(signature = (vcf, output, positions, sample=None, callable_bed=None, absent_is_ref=true))]
fn pharmcat_prepare(
    py: Python<'_>,
    vcf: &str,
    output: &str,
    positions: &str,
    sample: Option<String>,
    callable_bed: Option<String>,
    absent_is_ref: bool,
) -> PyResult<Py<PyAny>> {
    let op = pharmcat::Opcoes { amostra: sample, chamavel: callable_bed, ausente_e_ref: absent_is_ref };
    let r = py.allow_threads(|| pharmcat::preparar(vcf, positions, output, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// Runs PharmCAT on a person's VCF: prepares its input (see
/// `pharmcat_prepare`), calls `java -jar <jar>` and returns its calls per gene
/// (diplotype, phenotype, activity score...) and the paths of its reports
/// (HTML with the CPIC/DPWG recommendations, JSON, the calls TSV).
/// Experimental, for teaching and research only: not medical advice; a
/// pharmacogenomic result must be confirmed by an accredited laboratory.
///
/// - `jar`: pharmcat-<version>-all.jar (https://github.com/PharmGKB/PharmCAT/releases);
///   `positions` defaults to the pharmcat_positions_<version>.vcf next to it.
/// - `java`: the java program (Java 17 or newer).
/// - `outside_calls`: a PharmCAT outside-call TSV (for genes such as CYP2D6
///   or HLA-B, which PharmCAT does not call from a VCF).
#[pyfunction]
#[pyo3(signature = (vcf, output_dir, jar, positions=None, sample=None, callable_bed=None, absent_is_ref=true, java="java", outside_calls=None, base_name=None))]
#[allow(clippy::too_many_arguments)]
fn pharmcat_run(
    py: Python<'_>,
    vcf: &str,
    output_dir: &str,
    jar: &str,
    positions: Option<&str>,
    sample: Option<String>,
    callable_bed: Option<String>,
    absent_is_ref: bool,
    java: &str,
    outside_calls: Option<&str>,
    base_name: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let op = pharmcat::Opcoes { amostra: sample, chamavel: callable_bed, ausente_e_ref: absent_is_ref };
    let r = py
        .allow_threads(|| pharmcat::rodar(vcf, output_dir, jar, positions, &op, java, outside_calls, base_name))
        .map_err(erro_py)?;
    let d = json_para_py(py, &r.json())?;
    d.bind(py).set_item("notice", AVISO)?;
    Ok(d)
}

/// Calls pharmacogenetic star alleles, diplotypes, phenotypes, activity scores,
/// and CPIC drug recommendations directly from a GRCh38 VCF using MMGEN's
/// native `mmgen-pgx` engine and external CC0 JSON tables.
///
/// When several diplotypes fit the positions read, the highest-scoring one is
/// called, as PharmCAT's Named Allele Matcher does (each allele counts the
/// positions read that define it; ties are joined with " OR "; for NAT2, as PharmCAT
/// also does, the most common one in CPIC's population frequencies); the others are in
/// the gene's `notes` and `x_alternatives`, and a warning names those that would
/// change the phenotype or a recommendation. A gene is left uncalled only when
/// nothing could be read or no definition fits. With `bam=`, the HLA genes use the
/// IPD-IMGT/HLA release the tables pin (its SHA-256 is checked).
///
/// - `ancestry`: the person's ancestry proportions, e.g. {"EUR": 0.6, "AFR": 0.3,
///   "AMR": 0.1} (1000 Genomes superpopulations or CPIC group names), or
///   `ancestry_panel`: a reference panel to estimate them from the VCF (as
///   `estimate_ancestry`). With either, the frequency choice of NAT2 weights the CPIC
///   population groups by the person's ancestry instead of their plain mean: what
///   matters in admixed populations. The result says which was used (`x_ancestry`).
///
/// Writes `result.json` (matching `pgx_result.schema.json` v1) and `report.html`
/// inside `output_dir`, and returns the parsed result dictionary with `notice`.
#[pyfunction]
#[pyo3(signature = (vcf, output_dir, sample=None, callable_bed=None, bam=None, reference=None, tables=None, ancestry=None, ancestry_panel=None))]
#[allow(clippy::too_many_arguments)]
fn mmgen_pgx_run(
    py: Python<'_>,
    vcf: &str,
    output_dir: &str,
    sample: Option<String>,
    callable_bed: Option<String>,
    bam: Option<String>,
    reference: Option<String>,
    tables: Option<String>,
    ancestry: Option<std::collections::BTreeMap<String, f64>>,
    ancestry_panel: Option<String>,
) -> PyResult<Py<PyAny>> {
    let op = farmacogenetica::PgxOptions { sample, callable_bed, bam, reference, tables, ancestry, ancestry_panel };
    let r = py.allow_threads(|| farmacogenetica::run(vcf, output_dir, &op)).map_err(erro_py)?;
    let d = json_para_py(py, &r.json_text)?;
    d.bind(py).set_item("notice", AVISO)?;
    Ok(d)
}

/// Runs MMGEN's pharmacogenomics engine (`mmgen-pgx`) for a group of people and
/// the medicines of a protocol (e.g. the drugs given around a surgery).
///
/// - `vcfs`: one VCF per person, or multi-sample VCFs (every sample is a person;
///   `samples` keeps only some). `callable_beds`: one BED for all, or one per VCF
///   (None or "" for a VCF without one).
/// - `drugs`: the medicines (a list, or one text with commas or lines); each person
///   gets, for each one, a recommendation, or the reason there is none: a gene not
///   called (CYP2D6 and HLA need the reads), or no CPIC guideline in the tables.
/// - `ancestry_panel`: a reference panel to estimate each person's ancestry, which
///   weights the frequency choice of NAT2 (see `mmgen_pgx_run`).
///
/// Writes one folder per person (result.json, report.html) and two tables,
/// `cohort_genes.tsv` and `cohort_medicines.tsv`, inside `output_dir`; returns
/// `people`, `genes`, `medicines`, `medicines_without_guideline`, `files`,
/// `warnings`, `text` and `notice`.
#[pyfunction]
#[pyo3(signature = (vcfs, output_dir, callable_beds=None, samples=None, drugs=None, reference=None, tables=None, ancestry_panel=None))]
#[allow(clippy::too_many_arguments)]
fn pgx_cohort(
    py: Python<'_>,
    vcfs: &Bound<'_, PyAny>,
    output_dir: &str,
    callable_beds: Option<&Bound<'_, PyAny>>,
    samples: Option<&Bound<'_, PyAny>>,
    drugs: Option<&Bound<'_, PyAny>>,
    reference: Option<String>,
    tables: Option<String>,
    ancestry_panel: Option<String>,
) -> PyResult<Py<PyAny>> {
    let arquivos = lista_de_textos(Some(vcfs), "vcfs")?;
    let beds: Vec<Option<String>> = match callable_beds {
        None => Vec::new(),
        Some(b) => {
            if let Ok(s) = b.extract::<String>() {
                vec![Some(s).filter(|x| !x.is_empty())]
            } else {
                b.extract::<Vec<Option<String>>>()
                    .map_err(|_| pyo3::exceptions::PyTypeError::new_err("`callable_beds` must be a path or a list of paths (None for none)"))?
                    .into_iter()
                    .map(|x| x.filter(|s| !s.is_empty()))
                    .collect()
            }
        }
    };
    let amostras = samples.map(|s| lista_de_textos(Some(s), "samples")).transpose()?;
    let remedios = drugs
        .map(|d| lista_de_textos(Some(d), "drugs"))
        .transpose()?
        .map(|v| v.iter().flat_map(|s| s.split([',', ';', '\n'])).map(|s| s.trim().to_string()).filter(|s| !s.is_empty()).collect::<Vec<_>>());
    let op = farmacogenetica::CohortOptions { callable_beds: beds, samples: amostras, reference, tables, drugs: remedios, ancestry_panel };
    let texto = py.allow_threads(|| farmacogenetica::run_cohort(&arquivos, output_dir, &op)).map_err(erro_py)?;
    let d = json_para_py(py, &texto)?;
    d.bind(py).set_item("notice", AVISO)?;
    Ok(d)
}

/// Unified pharmacogenetics entry point: `engine="mmgen-pgx"` (default, MMGEN's own
/// engine, no Java needed; with `bam=` it also calls CYP2D6 with copy number and
/// HLA-A/B/C) or `engine="pharmcat"` (PharmCAT, Java 17+; `jar=` or PHARMCAT_JAR).
/// Both return the same dictionary (`pgx_result.schema.json`, version 1): `genes`
/// (with `called`, `diplotype`, `phenotype`, `evidence`, and `not_called_reason`
/// when a gene could not be called), `drugs`, `traits`, `warnings`, `report_html`,
/// `text` and `notice`. PharmCAT's own run details are under `x_pharmcat`.
/// `ancestry` / `ancestry_panel` (MMGEN's engine): see `mmgen_pgx_run`.
#[pyfunction]
#[pyo3(signature = (
    vcf,
    output_dir,
    sample=None,
    callable_bed=None,
    bam=None,
    reference=None,
    engine="mmgen-pgx",
    tables=None,
    jar=None,
    positions=None,
    absent_is_ref=false,
    java="java",
    outside_calls=None,
    base_name=None,
    ancestry=None,
    ancestry_panel=None
))]
#[allow(clippy::too_many_arguments)]
fn pgx_run(
    py: Python<'_>,
    vcf: &str,
    output_dir: &str,
    sample: Option<String>,
    callable_bed: Option<String>,
    bam: Option<String>,
    reference: Option<String>,
    engine: &str,
    tables: Option<String>,
    jar: Option<String>,
    positions: Option<String>,
    absent_is_ref: bool,
    java: &str,
    outside_calls: Option<String>,
    base_name: Option<String>,
    ancestry: Option<std::collections::BTreeMap<String, f64>>,
    ancestry_panel: Option<String>,
) -> PyResult<Py<PyAny>> {
    match engine {
        "mmgen-pgx" => mmgen_pgx_run(py, vcf, output_dir, sample, callable_bed, bam, reference, tables, ancestry, ancestry_panel),
        "pharmcat" => {
            let jar_path = jar
                .or_else(|| std::env::var("PHARMCAT_JAR").ok().filter(|s| !s.trim().is_empty()))
                .ok_or_else(|| {
                    pyo3::exceptions::PyValueError::new_err(
                        "engine=\"pharmcat\" requires `jar=<path/to/pharmcat-*-all.jar>` (or PHARMCAT_JAR set in environment).",
                    )
                })?;
            let op = pharmcat::Opcoes { amostra: sample, chamavel: callable_bed, ausente_e_ref: absent_is_ref };
            let r = py
                .allow_threads(|| {
                    pharmcat::rodar(
                        vcf,
                        output_dir,
                        &jar_path,
                        positions.as_deref(),
                        &op,
                        java,
                        outside_calls.as_deref(),
                        base_name.as_deref(),
                    )
                })
                .map_err(erro_py)?;
            let d = json_para_py(py, &r.esquema_comum().map_err(erro_py)?)?;
            d.bind(py).set_item("notice", AVISO)?;
            Ok(d)
        }
        other => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "Unknown PGx engine {other:?}; choose \"mmgen-pgx\" or \"pharmcat\"."
        ))),
    }
}

fn pessoa(path: &str, sample: Option<String>, bed: Option<String>, label: Option<String>) -> parentesco::Person {
    parentesco::Person { path: path.to_string(), sample, callable_bed: bed, label }
}

/// Relatedness of every pair of people: the KING-robust kinship (as plink2
/// --make-king-table), IBS0 and the inferred relationship (duplicate or
/// monozygotic twin, parent-offspring, full siblings, 2nd or 3rd degree,
/// unrelated). Experimental, for teaching and research only.
///
/// - `vcfs`: one VCF per person, or one multi-sample VCF (then `samples`
///   picks the samples; without it, all of them).
/// - `samples`, `callable_beds`, `labels`: per file (or per sample of a
///   single file); a single string applies to all.
/// - `absent_is_ref`: without a callable BED, is a site absent from a file
///   homozygous reference? (True for whole-genome variant calls.)
/// Files must be sorted in natural chromosome order (1..22, X, Y).
#[pyfunction]
#[pyo3(signature = (vcfs, samples=None, callable_beds=None, labels=None, absent_is_ref=true, pass_only=false))]
fn relatedness(
    py: Python<'_>,
    vcfs: &Bound<'_, PyAny>,
    samples: Option<&Bound<'_, PyAny>>,
    callable_beds: Option<&Bound<'_, PyAny>>,
    labels: Option<&Bound<'_, PyAny>>,
    absent_is_ref: bool,
    pass_only: bool,
) -> PyResult<Py<PyAny>> {
    let arquivos = lista_de_textos(Some(vcfs), "vcfs")?;
    let pessoas: Vec<parentesco::Person> = if arquivos.len() == 1 {
        let todas = match samples {
            Some(s) if !s.is_none() => lista_de_textos(Some(s), "samples")?,
            _ => py.allow_threads(|| parentesco::vcf_samples(&arquivos[0])).map_err(erro_py)?,
        };
        let beds = opcoes_de_lista(callable_beds, "callable_beds", todas.len())?;
        let rotulos = opcoes_de_lista(labels, "labels", todas.len())?;
        todas
            .into_iter()
            .zip(beds)
            .zip(rotulos)
            .map(|((s, b), l)| pessoa(&arquivos[0], Some(s), b, l))
            .collect()
    } else {
        let amostras = opcoes_de_lista(samples, "samples", arquivos.len())?;
        let beds = opcoes_de_lista(callable_beds, "callable_beds", arquivos.len())?;
        let rotulos = opcoes_de_lista(labels, "labels", arquivos.len())?;
        arquivos
            .iter()
            .zip(amostras)
            .zip(beds)
            .zip(rotulos)
            .map(|(((v, s), b), l)| pessoa(v, s, b, l))
            .collect()
    };
    let op = parentesco::JoinOptions { absent_is_ref, pass_only };
    let r = py.allow_threads(|| parentesco::relatedness(&pessoas, &op)).map_err(erro_py)?;
    let d = json_para_py(py, &r.json())?;
    d.bind(py).set_item("notice", AVISO)?;
    d.bind(py).set_item("text", r.resumo())?;
    Ok(d)
}

/// Trio analysis: Mendelian consistency of a child with the father and the
/// mother, de novo candidates, and a parentage test: the exclusions (opposite
/// homozygotes) and, given population allele frequencies, the combined
/// paternity index and the probability of paternity. Without `mother`, a duo.
/// Experimental, for teaching and research only: NOT a paternity test with
/// legal value (that needs an accredited laboratory and a chain of custody).
///
/// - `child`, `father`, `mother`: VCF paths; `*_sample` and `*_bed` pick the
///   sample and the callable regions of each.
/// - `allele_frequencies`: a sorted VCF of population frequencies (e.g. the
///   1000 Genomes sites file); `af_field` names its INFO field ("AF", "EUR_AF").
/// - `de_novo_output`: write the de novo candidates as a VCF.
#[pyfunction]
#[pyo3(signature = (child, father, mother=None, child_sample=None, father_sample=None, mother_sample=None,
                    child_bed=None, father_bed=None, mother_bed=None, allele_frequencies=None, af_field="AF",
                    genotype_error=0.001, prior=0.5, de_novo_output=None, absent_is_ref=true, pass_only=false))]
#[allow(clippy::too_many_arguments)]
fn trio_analysis(
    py: Python<'_>,
    child: &str,
    father: &str,
    mother: Option<&str>,
    child_sample: Option<String>,
    father_sample: Option<String>,
    mother_sample: Option<String>,
    child_bed: Option<String>,
    father_bed: Option<String>,
    mother_bed: Option<String>,
    allele_frequencies: Option<String>,
    af_field: &str,
    genotype_error: f64,
    prior: f64,
    de_novo_output: Option<String>,
    absent_is_ref: bool,
    pass_only: bool,
) -> PyResult<Py<PyAny>> {
    let c = pessoa(child, child_sample, child_bed, Some("child".into()));
    let f = pessoa(father, father_sample, father_bed, Some("father".into()));
    let m = mother.map(|p| pessoa(p, mother_sample, mother_bed, Some("mother".into())));
    let op = parentesco::TrioOptions {
        join: parentesco::JoinOptions { absent_is_ref, pass_only },
        af_path: allele_frequencies,
        af_field: af_field.to_string(),
        genotype_error,
        prior,
        de_novo_output,
    };
    let r = py.allow_threads(|| parentesco::trio_analysis(&c, &f, m.as_ref(), &op)).map_err(erro_py)?;
    let d = json_para_py(py, &r.json())?;
    d.bind(py).set_item("notice", AVISO)?;
    d.bind(py).set_item("text", r.resumo())?;
    Ok(d)
}

/// Runs of homozygosity of one person: stretches of at least `min_length`
/// bases and `min_sites` SNV calls (broken by gaps without data longer than
/// `max_gap`), written as BED, with F_ROH (their share of the autosomes).
/// Long runs come from related parents (consanguinity) or, on a single
/// chromosome, from uniparental disomy — as a SNP array shows. Experimental,
/// for teaching and research only.
///
/// By default, the sliding-window method of `plink --homozyg`: windows of
/// `window_snps` calls, each homozygous with at most `window_het` het and
/// `window_missing` missing calls; a call is in a run when at least
/// `window_threshold` of the windows over it are homozygous, so one
/// genotyping error does not split a run; runs sparser than one call per
/// `max_density` bases are dropped. `window_snps=0`: the strict method, where
/// any heterozygous call ends a run.
#[pyfunction]
#[pyo3(signature = (vcf, output=None, sample=None, callable_bed=None, min_sites=50, min_length=1_000_000, max_gap=1_000_000, pass_only=false,
                    window_snps=50, window_het=1, window_missing=5, window_threshold=0.05, max_density=50_000))]
#[allow(clippy::too_many_arguments)]
fn runs_of_homozygosity(
    py: Python<'_>,
    vcf: &str,
    output: Option<&str>,
    sample: Option<String>,
    callable_bed: Option<String>,
    min_sites: u64,
    min_length: u64,
    max_gap: u64,
    pass_only: bool,
    window_snps: u64,
    window_het: u64,
    window_missing: u64,
    window_threshold: f64,
    max_density: u64,
) -> PyResult<Py<PyAny>> {
    let p = pessoa(vcf, sample, callable_bed, None);
    let op = parentesco::RohOptions {
        min_sites,
        min_length,
        max_gap,
        pass_only,
        window_snps,
        window_het,
        window_missing,
        window_threshold,
        max_density,
        ..Default::default()
    };
    let r = py.allow_threads(|| parentesco::runs_of_homozygosity(&p, output, &op)).map_err(erro_py)?;
    let d = json_para_py(py, &r.json())?;
    d.bind(py).set_item("notice", AVISO)?;
    d.bind(py).set_item("text", r.resumo())?;
    Ok(d)
}

/// Genetic sex from a VCF: the heterozygosity of chrX outside the
/// pseudoautosomal regions (compared with the autosomes) and the calls on
/// chrY: XX, XY or undetermined, with flags for patterns that suggest X0 or
/// XXY. The file must have chrX (and chrY) calls. Experimental, for teaching
/// and research only.
#[pyfunction]
#[pyo3(signature = (vcf, sample=None, callable_bed=None, build=None, pass_only=false))]
fn genetic_sex(
    py: Python<'_>,
    vcf: &str,
    sample: Option<String>,
    callable_bed: Option<String>,
    build: Option<&str>,
    pass_only: bool,
) -> PyResult<Py<PyAny>> {
    let p = pessoa(vcf, sample, callable_bed, None);
    let r = py.allow_threads(|| parentesco::genetic_sex(&p, build, pass_only)).map_err(erro_py)?;
    let d = json_para_py(py, &r.json())?;
    d.bind(py).set_item("notice", AVISO)?;
    d.bind(py).set_item("text", r.resumo())?;
    Ok(d)
}

/// The people of `vcfs`: one per file, or, with a single file, one per sample
/// (`samples`, or every sample of the file).
fn pessoas_anc(
    py: Python<'_>,
    vcfs: &Bound<'_, PyAny>,
    samples: Option<&Bound<'_, PyAny>>,
    callable_beds: Option<&Bound<'_, PyAny>>,
) -> PyResult<Vec<ancestralidade::Person>> {
    let arquivos = lista_de_textos(Some(vcfs), "vcfs")?;
    if arquivos.len() == 1 {
        let todas = match samples {
            Some(s) if !s.is_none() => lista_de_textos(Some(s), "samples")?,
            _ => py.allow_threads(|| parentesco::vcf_samples(&arquivos[0])).map_err(erro_py)?,
        };
        let beds = opcoes_de_lista(callable_beds, "callable_beds", todas.len())?;
        Ok(todas
            .into_iter()
            .zip(beds)
            .map(|(s, b)| ancestralidade::Person { vcf: arquivos[0].clone(), sample: Some(s), callable: b })
            .collect())
    } else {
        let amostras = opcoes_de_lista(samples, "samples", arquivos.len())?;
        let beds = opcoes_de_lista(callable_beds, "callable_beds", arquivos.len())?;
        Ok(arquivos
            .iter()
            .zip(amostras)
            .zip(beds)
            .map(|((v, s), b)| ancestralidade::Person { vcf: v.clone(), sample: s, callable: b })
            .collect())
    }
}

fn lista_json_para_py<T>(py: Python<'_>, itens: &[T], json: fn(&T) -> String, resumo: fn(&T) -> String) -> PyResult<Py<PyAny>> {
    let texto = format!("[{}]", itens.iter().map(json).collect::<Vec<_>>().join(","));
    let lista = json_para_py(py, &texto)?;
    for (i, it) in itens.iter().enumerate() {
        let d = lista.bind(py).get_item(i)?;
        d.set_item("text", resumo(it))?;
        d.set_item("notice", AVISO)?;
    }
    Ok(lista)
}

/// Builds a population frequency panel (a sorted TSV: chrom, pos, ref, alt,
/// AF_<group>..., AN_<group>...) from a multi-sample VCF (e.g. the 1000
/// Genomes genotypes) and a table that gives each sample's group. Keeps
/// biallelic SNVs with a minor-allele frequency of at least `min_maf` in some
/// group, optionally one per `thin_bp` bases. `exclude`: samples left out (the
/// relatives, or people held out to test the panel).
#[pyfunction]
#[pyo3(signature = (vcfs, samples_table, output, sample_column, group_column, exclude=None, min_maf=0.05, thin_bp=0, keep_indels=false, pass_only=true))]
#[allow(clippy::too_many_arguments)]
fn build_panel(
    py: Python<'_>,
    vcfs: &Bound<'_, PyAny>,
    samples_table: &str,
    output: &str,
    sample_column: &str,
    group_column: &str,
    exclude: Option<&Bound<'_, PyAny>>,
    min_maf: f64,
    thin_bp: u64,
    keep_indels: bool,
    pass_only: bool,
) -> PyResult<Py<PyAny>> {
    let arquivos = lista_de_textos(Some(vcfs), "vcfs")?;
    let fora = lista_de_textos(exclude, "exclude")?;
    let op = ancestralidade::PanelOptions { min_maf, thin_bp, keep_indels, pass_only };
    let r = py
        .allow_threads(|| ancestralidade::build_panel(&arquivos, samples_table, sample_column, group_column, &fora, output, &op))
        .map_err(erro_py)?;
    let d = json_para_py(py, &r.json())?;
    d.bind(py).set_item("text", r.resumo())?;
    Ok(d)
}

/// Builds a population frequency panel from a sites VCF whose INFO holds each
/// group's frequency, e.g. the 1000 Genomes GRCh38 sites file:
/// `groups=["AFR","AMR","EAS","EUR","SAS"]`, `af_fields=["AFR_AF",...]`.
#[pyfunction]
#[pyo3(signature = (vcfs, groups, af_fields, output, an_fields=None, min_maf=0.05, thin_bp=0, keep_indels=false, pass_only=true))]
#[allow(clippy::too_many_arguments)]
fn panel_from_frequencies(
    py: Python<'_>,
    vcfs: &Bound<'_, PyAny>,
    groups: &Bound<'_, PyAny>,
    af_fields: &Bound<'_, PyAny>,
    output: &str,
    an_fields: Option<&Bound<'_, PyAny>>,
    min_maf: f64,
    thin_bp: u64,
    keep_indels: bool,
    pass_only: bool,
) -> PyResult<Py<PyAny>> {
    let arquivos = lista_de_textos(Some(vcfs), "vcfs")?;
    let grupos = lista_de_textos(Some(groups), "groups")?;
    let afs = lista_de_textos(Some(af_fields), "af_fields")?;
    let ans = lista_de_textos(an_fields, "an_fields")?;
    let op = ancestralidade::PanelOptions { min_maf, thin_bp, keep_indels, pass_only };
    let r = py
        .allow_threads(|| ancestralidade::panel_from_frequencies(&arquivos, &grupos, &afs, &ans, output, &op))
        .map_err(erro_py)?;
    let d = json_para_py(py, &r.json())?;
    d.bind(py).set_item("text", r.resumo())?;
    Ok(d)
}

/// Ancestry proportions of each person against the groups of a frequency
/// panel (`build_panel`, `panel_from_frequencies`): the maximum-likelihood
/// admixture proportions with the panel's frequencies fixed (as ADMIXTURE's
/// supervised mode), standard errors from a block bootstrap, and the best
/// single group. Experimental, for teaching and research only: reference
/// groups are broad continental samples, not ethnic or national identities.
///
/// Returns one dictionary per person.
#[pyfunction]
#[pyo3(signature = (panel, vcfs, samples=None, callable_beds=None, absent_is_ref=true, pass_only=true, bootstrap=100, block_bp=5_000_000, seed=1))]
#[allow(clippy::too_many_arguments)]
fn estimate_ancestry(
    py: Python<'_>,
    panel: &str,
    vcfs: &Bound<'_, PyAny>,
    samples: Option<&Bound<'_, PyAny>>,
    callable_beds: Option<&Bound<'_, PyAny>>,
    absent_is_ref: bool,
    pass_only: bool,
    bootstrap: usize,
    block_bp: u64,
    seed: u64,
) -> PyResult<Py<PyAny>> {
    let pessoas = pessoas_anc(py, vcfs, samples, callable_beds)?;
    let op = ancestralidade::AncestryOptions { absent_is_ref, pass_only, bootstrap, block_bp, seed, ..Default::default() };
    let r = py.allow_threads(|| ancestralidade::estimate_ancestry(panel, &pessoas, &op)).map_err(erro_py)?;
    lista_json_para_py(py, &r, |x| x.json(), |x| x.resumo())
}

/// A polygenic score (a PGS Catalog scoring file, GRCh38 harmonized) for each
/// person: the weighted sum of effect-allele dosages, its coverage, and, with
/// a frequency `panel`, where it falls in each panel group (mean, sd, z,
/// percentile, under Hardy-Weinberg and independent variants). Experimental,
/// for teaching and research only: a polygenic score is a statistical
/// summary, not a diagnosis, and it transfers poorly across ancestries.
///
/// - `impute`: for score variants the person's file cannot call: "pooled"
///   (the panel's pooled frequency), a group name, or a list of proportions.
/// - `keep_ambiguous`: keep A/T and C/G SNVs (their strand cannot be checked).
/// - `details`: write a per-variant table (TSV).
#[pyfunction]
#[pyo3(signature = (scoring_file, vcfs, samples=None, callable_beds=None, panel=None, reference=None, impute=None,
                    keep_ambiguous=false, allow_flips=true, details=None, absent_is_ref=true, pass_only=true))]
#[allow(clippy::too_many_arguments)]
fn polygenic_score(
    py: Python<'_>,
    scoring_file: &str,
    vcfs: &Bound<'_, PyAny>,
    samples: Option<&Bound<'_, PyAny>>,
    callable_beds: Option<&Bound<'_, PyAny>>,
    panel: Option<String>,
    reference: Option<String>,
    impute: Option<&Bound<'_, PyAny>>,
    keep_ambiguous: bool,
    allow_flips: bool,
    details: Option<String>,
    absent_is_ref: bool,
    pass_only: bool,
) -> PyResult<Py<PyAny>> {
    let pessoas = pessoas_anc(py, vcfs, samples, callable_beds)?;
    let imputar = match impute {
        None => ancestralidade::Impute::Pooled,
        Some(v) if v.is_none() => ancestralidade::Impute::Pooled,
        Some(v) => {
            if let Ok(s) = v.extract::<String>() {
                if s.eq_ignore_ascii_case("pooled") {
                    ancestralidade::Impute::Pooled
                } else {
                    ancestralidade::Impute::Group(s)
                }
            } else {
                let p: Vec<f64> = v.extract().map_err(|_| {
                    pyo3::exceptions::PyTypeError::new_err("`impute` must be \"pooled\", a group name or a list of proportions")
                })?;
                ancestralidade::Impute::Proportions(p)
            }
        }
    };
    let op = ancestralidade::ScoreOptions { panel, reference, absent_is_ref, pass_only, keep_ambiguous, allow_flips, impute: imputar, details };
    let r = py.allow_threads(|| ancestralidade::polygenic_score(scoring_file, &pessoas, &op)).map_err(erro_py)?;
    lista_json_para_py(py, &r, |x| x.json(), |x| x.resumo())
}

/// Counts the reads (fragments) of RNA-seq alignments per gene, as
/// featureCounts / htseq-count do in union mode: from coordinate-sorted
/// BAM/SAM files and a GFF3/GTF annotation. Writes a count matrix (gene, one
/// column per sample) when `output` is given, and returns the assignment
/// summary of each sample (assigned, ambiguous, no feature, multimapping...).
///
/// - `feature_type`: the annotation feature counted ("exon").
/// - `attribute`: the attribute that groups features into genes; default:
///   gene_id, else the parent gene; "GeneID" reads Dbxref=GeneID (RefSeq).
/// - `strandedness`: "unstranded", "forward" or "reverse".
#[pyfunction]
#[pyo3(signature = (bams, annotation, output=None, feature_type="exon", attribute=None, strandedness="unstranded", min_mapq=10, sample_names=None, summary_output=None))]
#[allow(clippy::too_many_arguments)]
fn count_reads(
    py: Python<'_>,
    bams: &Bound<'_, PyAny>,
    annotation: &str,
    output: Option<&str>,
    feature_type: &str,
    attribute: Option<String>,
    strandedness: &str,
    min_mapq: u8,
    sample_names: Option<&Bound<'_, PyAny>>,
    summary_output: Option<String>,
) -> PyResult<Py<PyAny>> {
    let arquivos = lista_de_textos(Some(bams), "bams")?;
    let nomes_amostras = lista_de_textos(sample_names, "sample_names")?;
    let op = expressao::CountOptions {
        feature_type: feature_type.to_string(),
        attribute,
        strandedness: expressao::Strandedness::from_text(strandedness).map_err(erro_py)?,
        min_mapq,
        sample_names: if nomes_amostras.is_empty() { None } else { Some(nomes_amostras) },
        summary_output,
    };
    let r = py.allow_threads(|| expressao::count_reads(&arquivos, annotation, output, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// Differential expression in the manner of DESeq2: size factors by the
/// median of ratios, gene-wise dispersions shrunk toward a fitted trend, a
/// negative-binomial GLM and a Wald test for `condition` (treated against
/// `reference`), with optional blocking `factors` (e.g. the cell line),
/// Benjamini-Hochberg adjustment, independent filtering and the Cook's
/// distance outlier flag.
///
/// - `counts`: a count matrix (gene + one column per sample; featureCounts and
///   htseq-count tables work). `samples`: a table with a column per factor.
/// - `gene_names`: NCBI gene_info, or a two-column table (id, name).
/// Returns the summary and, in "plot", the per-gene values for MA and volcano plots.
#[pyfunction]
#[pyo3(signature = (counts, samples, condition, reference, output=None, treated=None, factors=None, alpha=0.1, lfc_threshold=0.0,
                    independent_filtering=true, cooks_filter=true, gene_names=None, only_compared_levels=false))]
#[allow(clippy::too_many_arguments)]
fn differential_expression(
    py: Python<'_>,
    counts: &str,
    samples: &str,
    condition: &str,
    reference: &str,
    output: Option<&str>,
    treated: Option<String>,
    factors: Option<&Bound<'_, PyAny>>,
    alpha: f64,
    lfc_threshold: f64,
    independent_filtering: bool,
    cooks_filter: bool,
    gene_names: Option<String>,
    only_compared_levels: bool,
) -> PyResult<Py<PyAny>> {
    let op = expressao::DeOptions {
        condition: condition.to_string(),
        reference: reference.to_string(),
        treated,
        factors: lista_de_textos(factors, "factors")?,
        only_compared_levels,
        alpha,
        lfc_threshold,
        independent_filtering,
        cooks_filter,
        gene_names,
        ..Default::default()
    };
    let r = py.allow_threads(|| expressao::differential_expression(counts, samples, output, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// Tests gene sets against a differential-expression table (the one
/// `differential_expression` writes, or DESeq2's from R): for each set, the
/// Wilcoxon rank-sum test of its genes' statistic against the other genes tested,
/// as limma's `wilcoxGST` (a competitive test that assumes independent genes):
/// does the set go up or down (`pvalue`, `fdr`), or change in either direction
/// (`pvalue_mixed`, `fdr_mixed`), more than the rest?
///
/// - `gene_sets`: a GMT file (name, description, genes), a two-column table (set,
///   gene), or the built-in "spinal_cord_injury" (reactive astrocytes, glial-scar
///   proteoglycans, fibrotic scar, inflammation, regeneration-associated genes,
///   axonal cytoskeleton, myelin, growth inhibition and RhoA/ROCK, neurotrophins,
///   laminin receptors and chains: marker genes from reviews, for teaching and
///   exploratory analysis, not validated signatures). Gene symbols are matched
///   without case (human GFAP = rodent Gfap).
/// - `statistic`: the column the genes are ranked by ("stat", or "log2FoldChange").
/// - `require_padj`: rank only the genes with an adjusted p-value, those DESeq2
///   kept after independent filtering (as limma users filter lowly expressed genes
///   first); otherwise thousands of barely expressed genes, with statistics near
///   zero, make any set of expressed genes look changed.
#[pyfunction]
#[pyo3(signature = (results, gene_sets="spinal_cord_injury", output=None, statistic="stat", alpha=0.1, min_size=3, require_padj=true))]
#[allow(clippy::too_many_arguments)]
fn gene_set_test(
    py: Python<'_>,
    results: &str,
    gene_sets: &str,
    output: Option<&str>,
    statistic: &str,
    alpha: f64,
    min_size: usize,
    require_padj: bool,
) -> PyResult<Py<PyAny>> {
    let op = conjuntos::SetOptions { statistic: statistic.to_string(), alpha, min_size, require_padj };
    let texto = py.allow_threads(|| conjuntos::gene_set_test(results, gene_sets, output, &op)).map_err(erro_py)?;
    json_para_py(py, &texto)
}

/// Association between genotypes and an outcome in a cohort: one test per ALT
/// allele of each selected variant, adjusted for covariates. Quantitative
/// outcome: linear regression (t test). Binary outcome: logistic regression
/// (likelihood-ratio p, Wald p too), with Firth's penalized likelihood when the
/// plain fit separates or fails (`firth="auto"`), always, or never.
/// Benjamini-Hochberg FDR; with `permutations`, an empirical p per variant and a
/// family-wise (max-T) p.
///
/// - `vcfs`: one multi-sample VCF or several (people matched by sample name;
///   several files must be sorted). A position missing from a person's file is a
///   missing genotype, unless `absent_is_ref=True`.
/// - `phenotypes`: a table (tab, comma or space separated) with a header; the
///   first column (or `IID`, or `sample_column`) names the people. `outcome` is
///   the column tested; `covariates` the columns adjusted for. Numeric
///   covariates are numbers, text ones categories (the first level in sorted
///   order is the reference; `categorical` forces a column to be a category).
/// - `outcome_type`: "auto" (two values: binary), "binary" or "quantitative";
///   `case_value`: which value is the event (binary).
/// - `ids`: rsIDs, "chr:pos", "chr:pos:ref:alt"; `regions`: "chr:start-end".
///   Without either, every variant with MAF >= `min_maf` and call rate >=
///   `min_call_rate`.
/// - `model`: "additive", "dominant" or "recessive". `output`: the table (TSV).
///
/// Returns `people` (used and dropped, and why), `variants` (read, tested,
/// skipped and why, methods), `rows` (one per tested allele: beta, se, ci, or,
/// p, method, p_permutation, p_fwer, fdr...), `top`, `power` (the effect the
/// study could detect with 80% power), `warnings`, `text` and `notice`.
/// Only hard GT calls are read; there is no relatedness correction (run
/// `relatedness` first). An association is not causation.
#[pyfunction]
#[pyo3(signature = (vcfs, phenotypes, outcome, covariates=None, categorical=None, sample_column=None, outcome_type="auto", case_value=None, model="additive", ids=None, regions=None, min_maf=0.01, min_call_rate=0.9, firth="auto", permutations=0, seed=1, absent_is_ref=false, pass_only=false, output=None))]
#[allow(clippy::too_many_arguments)]
fn genotype_association(
    py: Python<'_>,
    vcfs: &Bound<'_, PyAny>,
    phenotypes: &str,
    outcome: &str,
    covariates: Option<&Bound<'_, PyAny>>,
    categorical: Option<&Bound<'_, PyAny>>,
    sample_column: Option<String>,
    outcome_type: &str,
    case_value: Option<String>,
    model: &str,
    ids: Option<&Bound<'_, PyAny>>,
    regions: Option<&Bound<'_, PyAny>>,
    min_maf: f64,
    min_call_rate: f64,
    firth: &str,
    permutations: usize,
    seed: u64,
    absent_is_ref: bool,
    pass_only: bool,
    output: Option<String>,
) -> PyResult<Py<PyAny>> {
    let vcfs = lista_de_textos(Some(vcfs), "vcfs")?;
    let op = associacao::AssocOptions {
        outcome: outcome.to_string(),
        covariates: lista_de_textos(covariates, "covariates")?,
        categorical: lista_de_textos(categorical, "categorical")?,
        sample_column,
        outcome_type: associacao::OutcomeType::from_text(outcome_type).map_err(erro_py)?,
        case_value,
        model: associacao::GeneticModel::from_text(model).map_err(erro_py)?,
        ids: lista_de_textos(ids, "ids")?,
        regions: lista_de_textos(regions, "regions")?,
        min_maf,
        min_call_rate,
        firth: associacao::FirthPolicy::from_text(firth).map_err(erro_py)?,
        permutations,
        seed,
        absent_is_ref,
        pass_only,
        output,
    };
    let r = py.allow_threads(|| associacao::genotype_association(&vcfs, phenotypes, &op)).map_err(erro_py)?;
    let d = json_para_py(py, &r.json())?;
    d.bind(py).set_item("notice", AVISO)?;
    Ok(d)
}

/// Power of a genotype association test, the people needed for
/// `target_power`, and the smallest effect detectable with `n` people.
///
/// - `outcome="binary"`: `effect` is the odds ratio per copy (additive model),
///   `allele_freq` the frequency in people without the event, `case_fraction`
///   the fraction with it (Wald test, large-sample).
/// - `outcome="quantitative"`: `effect` is beta in standard deviations of the
///   outcome per copy (non-central t).
/// - `tests` applies Bonferroni (alpha / tests); `covariates` is how many the
///   model adjusts for. `simulations` adds a simulated power (the default tests
///   of `genotype_association`, with Firth when needed): use it for small
///   studies.
///
/// Returns `power`, `n_needed`, `detectable_effect`, `simulated_power`,
/// `expected_events`, `expected_carriers`, `curve` (n and power), `warnings`
/// and `text`.
#[pyfunction]
#[pyo3(signature = (allele_freq, effect=None, n=None, outcome="binary", case_fraction=0.5, alpha=0.05, tests=1, model="additive", covariates=0, target_power=0.8, simulations=0, seed=1))]
#[allow(clippy::too_many_arguments)]
fn association_power(
    py: Python<'_>,
    allele_freq: f64,
    effect: Option<f64>,
    n: Option<usize>,
    outcome: &str,
    case_fraction: f64,
    alpha: f64,
    tests: usize,
    model: &str,
    covariates: usize,
    target_power: f64,
    simulations: usize,
    seed: u64,
) -> PyResult<Py<PyAny>> {
    let op = associacao::PowerOptions {
        outcome: associacao::OutcomeType::from_text(outcome).map_err(erro_py)?,
        allele_freq,
        effect,
        n,
        case_fraction,
        alpha,
        tests,
        model: associacao::GeneticModel::from_text(model).map_err(erro_py)?,
        covariates,
        target_power,
        simulations,
        seed,
    };
    let r = py.allow_threads(|| associacao::association_power(&op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// Builds the database of an adventitious-agent screen: the genomes of the
/// agents to look for (viruses, bacteria, mycoplasmas, fungi), of the
/// process controls and of the host, in one file. Agents and controls are
/// indexed by minimizers with their genome windows (for the breadth of
/// coverage); the host goes into a compact set (a human genome takes ~0.7
/// GB and is built in passes within `memory_mb`). The SHA-256 of every
/// reference file is kept. For research and teaching: not a validated
/// method (see `notice`).
///
/// - `references`: one entry per genome: a FASTA path, or a list of paths
///   (one genome in several files: chromosomes, segments, plasmids).
/// - `kinds`: one per genome: "virus", "bacterium", "mycoplasma", "fungus",
///   "other", "spike" (a spike-in or process control) or "host".
/// - `names`: the names in the reports (default: the first FASTA header).
/// - `minimizer_window`: one k-mer in about (w + 1) / 2 is kept for the
///   agents; `host_minimizer_window` the same for the host (sparser).
/// - `genome_window`: bases per genome window for the breadth (a small
///   genome gets length / 20, at least 100).
#[pyfunction]
#[pyo3(signature = (references, output, kinds, names=None, k=31, minimizer_window=12, host_minimizer_window=31,
                    genome_window=1000, memory_mb=1024))]
#[allow(clippy::too_many_arguments)]
fn build_screen_database(
    py: Python<'_>,
    references: &Bound<'_, PyAny>,
    output: &str,
    kinds: Vec<String>,
    names: Option<&Bound<'_, PyAny>>,
    k: u8,
    minimizer_window: u16,
    host_minimizer_window: u16,
    genome_window: u32,
    memory_mb: u64,
) -> PyResult<Py<PyAny>> {
    let listas: Vec<Vec<String>> = if let Ok(s) = references.extract::<String>() {
        vec![vec![s]]
    } else if let Ok(v) = references.extract::<Vec<String>>() {
        v.into_iter().map(|f| vec![f]).collect()
    } else {
        references.extract::<Vec<Vec<String>>>().map_err(|_| {
            pyo3::exceptions::PyTypeError::new_err("`references` must be a path, a list of paths, or a list of lists of paths")
        })?
    };
    if kinds.len() != listas.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "`kinds` needs one kind per reference ({} references, {} kinds): virus, bacterium, mycoplasma, fungus, other, spike or host",
            listas.len(),
            kinds.len()
        )));
    }
    let tipos: Vec<contaminantes::AgentKind> =
        kinds.iter().map(|s| contaminantes::AgentKind::from_text(s)).collect::<Result<_, _>>().map_err(erro_py)?;
    let nomes = opcoes_de_lista(names, "names", listas.len())?;
    let refs: Vec<contaminantes::ScreenReference> = listas
        .into_iter()
        .zip(nomes)
        .zip(tipos)
        .map(|((files, name), kind)| contaminantes::ScreenReference { files, name: name.unwrap_or_default(), kind })
        .collect();
    let op = contaminantes::ScreenDbOptions { k, minimizer_window, host_minimizer_window, genome_window_bp: genome_window, memory_mb };
    let r = py.allow_threads(|| contaminantes::build_screen_database(&refs, output, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// What a screen database holds: its references (name, kind, files and
/// their SHA-256, sizes, genome windows), its parameters and its SHA-256.
#[pyfunction]
fn describe_screen_database(py: Python<'_>, path: &str) -> PyResult<Py<PyAny>> {
    let r = py.allow_threads(|| contaminantes::describe_screen_database(path)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// Screens the sequencing reads of a biological material (a raw material,
/// a cell bank, a harvest, a tissue extract) for adventitious agents with a
/// database of `build_screen_database`: host reads set aside, each read (or
/// pair) given to the reference it matches best, and for each organism the
/// evidence -- unique and tied reads, reads per million (RPM), distinct
/// k-mers, the genome windows hit (breadth) and whether they are spread as
/// an infection would spread them -- and a call: "detected", "background"
/// (also in a negative control, not 10x above it), "below_threshold" (with
/// the reasons) or "not_detected". Spike-ins get their recovery, and the
/// run its limit of detection. The report (TSV, JSON, HTML) records the
/// SHA-256 of every input and of the database, all parameters, the date and
/// the engine version.
///
/// For research and teaching: releasing a GMP batch requires a validated
/// method and validated software (ICH Q5A(R2) admits NGS for adventitious
/// virus detection, but the validation is the laboratory's).
///
/// - `samples`: a list of (name, r1) or (name, r1, r2) -- r1 and r2 a path
///   or a list of paths (lanes, in the same order) -- or of dicts with the
///   keys name, r1 and r2.
/// - `negative_controls`: names of the samples that are blanks.
/// - The thresholds of a call: `min_reads` own reads, `min_windows` genome
///   windows, `min_dispersion` (windows hit / windows expected at random for
///   that many reads), `min_unique_fraction` (own reads among those that match
///   it best) and `max_weak_fraction` (reads that only resemble it); a sample
///   is "background" unless its RPM is `min_blank_ratio` times a blank's.
/// - `spike_ins`: a TSV path (columns sample, organism, amount, unit; sample
///   empty or "*" = every sample) or a list of dicts with those keys. unit:
///   "reads", "rpm", or any other (copies, CFU...: recovery in reads per
///   unit).
#[pyfunction]
#[pyo3(signature = (samples, database, negative_controls=None, spike_ins=None, output_tsv=None, output_json=None,
                    output_html=None, min_hits=3, min_fraction=0.2, unique_margin=0.2, host_min_hits=1, min_reads=3,
                    min_windows=2, min_dispersion=0.5, min_unique_fraction=0.25, max_weak_fraction=0.5, min_rpm=0.0,
                    min_blank_ratio=10.0, lod_rate=0.95))]
#[allow(clippy::too_many_arguments)]
fn contaminant_screen(
    py: Python<'_>,
    samples: &Bound<'_, PyAny>,
    database: &str,
    negative_controls: Option<&Bound<'_, PyAny>>,
    spike_ins: Option<&Bound<'_, PyAny>>,
    output_tsv: Option<String>,
    output_json: Option<String>,
    output_html: Option<String>,
    min_hits: u32,
    min_fraction: f64,
    unique_margin: f64,
    host_min_hits: u32,
    min_reads: u64,
    min_windows: u32,
    min_dispersion: f64,
    min_unique_fraction: f64,
    max_weak_fraction: f64,
    min_rpm: f64,
    min_blank_ratio: f64,
    lod_rate: f64,
) -> PyResult<Py<PyAny>> {
    use pyo3::types::{PyDict, PyList};
    let paths = |v: Option<Bound<'_, PyAny>>| -> PyResult<Vec<String>> {
        match v {
            None => Ok(Vec::new()),
            Some(x) if x.is_none() => Ok(Vec::new()),
            Some(x) => lista_de_textos(Some(&x), "r1/r2"),
        }
    };
    let mut amostras = Vec::new();
    let itens: Vec<Bound<'_, PyAny>> = samples
        .extract()
        .map_err(|_| pyo3::exceptions::PyTypeError::new_err("`samples` must be a list of (name, r1[, r2]) or of dicts"))?;
    for item in itens {
        let (name, r1, r2) = if let Ok(d) = item.downcast::<PyDict>() {
            let name: String = d.get_item("name")?.ok_or_else(|| pyo3::exceptions::PyValueError::new_err("a sample dict needs \"name\""))?.extract()?;
            (name, paths(d.get_item("r1")?)?, paths(d.get_item("r2")?)?)
        } else {
            let t: Vec<Bound<'_, PyAny>> = item.extract().map_err(|_| {
                pyo3::exceptions::PyTypeError::new_err("each sample is (name, r1), (name, r1, r2) or a dict with name, r1 and r2")
            })?;
            if t.len() < 2 || t.len() > 3 {
                return Err(pyo3::exceptions::PyValueError::new_err("each sample is (name, r1) or (name, r1, r2)"));
            }
            (t[0].extract()?, paths(Some(t[1].clone()))?, paths(t.get(2).cloned())?)
        };
        amostras.push(contaminantes::ScreenSample { name, r1, r2 });
    }
    let mut op = contaminantes::ScreenOptions {
        negative_controls: lista_de_textos(negative_controls, "negative_controls")?,
        min_hits,
        min_fraction,
        unique_margin,
        host_min_hits,
        min_reads,
        min_windows,
        min_dispersion,
        min_unique_fraction,
        max_weak_fraction,
        min_rpm,
        min_blank_ratio,
        lod_rate,
        output_tsv,
        output_json,
        output_html,
        ..Default::default()
    };
    if let Some(sp) = spike_ins.filter(|x| !x.is_none()) {
        if let Ok(p) = sp.extract::<String>() {
            op.spike_table = Some(p);
        } else {
            for d in sp.downcast::<PyList>().map_err(|_| pyo3::exceptions::PyTypeError::new_err("`spike_ins` must be a TSV path or a list of dicts"))?.iter() {
                let d = d.downcast::<PyDict>().map_err(|_| pyo3::exceptions::PyTypeError::new_err("each spike-in is a dict: organism, amount, and optionally sample and unit"))?;
                let get = |k: &str| d.get_item(k).ok().flatten().filter(|v| !v.is_none());
                op.spike_ins.push(contaminantes::SpikeIn {
                    sample: get("sample").map(|v| v.extract::<String>()).transpose()?.filter(|s| !s.is_empty() && s != "*"),
                    organism: get("organism").ok_or_else(|| pyo3::exceptions::PyValueError::new_err("a spike-in needs \"organism\""))?.extract()?,
                    amount: get("amount").ok_or_else(|| pyo3::exceptions::PyValueError::new_err("a spike-in needs \"amount\""))?.extract()?,
                    unit: get("unit").map(|v| v.extract::<String>()).transpose()?.unwrap_or_else(|| "reads".to_string()),
                });
            }
        }
    }
    let r = py.allow_threads(|| contaminantes::contaminant_screen(&amostras, database, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// Builds the pseudoalignment index of a transcriptome, for
/// `quantify_samples` (the method of kallisto): every k-mer of the
/// transcripts, on both strands, mapped to the set of transcripts that contain
/// it. Reads a FASTA of transcript sequences (cDNA: Ensembl
/// `*.cdna.all.fa.gz`, GENCODE `*.transcripts.fa.gz`, RefSeq `*_rna.fna.gz`),
/// plain or .gz. The index is a file; build it once per transcriptome.
///
/// - `k`: k-mer length, odd, 7 to 31; 31 (kallisto's default) suits reads of
///   50 bp or more.
/// - `gene_map`: a GTF, a GFF3, or a table (transcript, gene[, gene name]);
///   by default the gene of each transcript is read from the FASTA headers
///   (Ensembl `gene:`/`gene_symbol:`, GENCODE `|` fields, RefSeq `(SYMBOL)`),
///   and a transcript without one is its own gene.
/// - `make_unique`: rename repeated transcript names NAME_2, NAME_3... instead
///   of stopping.
/// - `memory_mb`: memory for one building pass (default 1024); less memory
///   means more passes, never a failure.
/// - `clip_poly_a`: leave poly-A tails (10 or more A's at the end) out of the
///   k-mers, as kallisto does.
///
/// Returns the summary: `transcripts`, `genes`, `gene_map`, `kmers`,
/// `equivalence_classes`, `bytes`, `bytes_per_million_kmers`, `passes`,
/// `seconds`, `warnings`, `text`...
#[pyfunction]
#[pyo3(signature = (fasta, output, k=31, gene_map=None, make_unique=false, memory_mb=None, clip_poly_a=true))]
#[allow(clippy::too_many_arguments)]
fn build_transcript_index(
    py: Python<'_>,
    fasta: &str,
    output: &str,
    k: usize,
    gene_map: Option<String>,
    make_unique: bool,
    memory_mb: Option<usize>,
    clip_poly_a: bool,
) -> PyResult<Py<PyAny>> {
    let op = quantificacao::IndexOptions { k, gene_map, make_unique, memory_mb, clip_poly_a };
    let r = py.allow_threads(|| quantificacao::build_transcript_index(fasta, output, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// What a transcript index holds -- k, transcripts, genes, k-mers,
/// equivalence classes, the FASTA and gene map it was built from, its size --
/// read from the header of the file (fast).
#[pyfunction]
fn transcript_index_info(py: Python<'_>, index: &str) -> PyResult<Py<PyAny>> {
    let r = quantificacao::transcript_index_info(index).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// One item of `samples`: a dict {"name", "r1", "r2", "condition"} ("r2" and
/// "condition" optional), or a tuple/list (name, r1[, r2[, condition]]).
fn amostra_quant(item: &Bound<'_, PyAny>, i: usize) -> PyResult<quantificacao::QuantSample> {
    use pyo3::types::PyDict;
    let texto = |v: Option<Bound<'_, PyAny>>, campo: &str| -> PyResult<Option<String>> {
        match v {
            Some(v) if !v.is_none() => v.extract::<String>().map(Some).map_err(|_| {
                pyo3::exceptions::PyTypeError::new_err(format!("sample {}: `{campo}` must be a string", i + 1))
            }),
            _ => Ok(None),
        }
    };
    if let Ok(d) = item.downcast::<PyDict>() {
        Ok(quantificacao::QuantSample {
            name: texto(d.get_item("name")?, "name")?.unwrap_or_default(),
            r1: texto(d.get_item("r1")?, "r1")?.ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!("sample {}: no \"r1\" (the FASTQ file, or the R1 file of a pair)", i + 1))
            })?,
            r2: texto(d.get_item("r2")?, "r2")?,
            condition: texto(d.get_item("condition")?, "condition")?,
        })
    } else {
        let v: Vec<Option<String>> = item.extract().map_err(|_| {
            pyo3::exceptions::PyTypeError::new_err(format!(
                "sample {}: give a dict {{\"name\", \"r1\", \"r2\", \"condition\"}} or a tuple (name, r1[, r2[, condition]])",
                i + 1
            ))
        })?;
        if v.len() < 2 || v.len() > 4 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "sample {}: a tuple (name, r1[, r2[, condition]]) has 2 to 4 items, not {}",
                i + 1,
                v.len()
            )));
        }
        let pega = |j: usize| v.get(j).cloned().flatten().filter(|s| !s.is_empty());
        Ok(quantificacao::QuantSample {
            name: pega(0).unwrap_or_default(),
            r1: pega(1).ok_or_else(|| pyo3::exceptions::PyValueError::new_err(format!("sample {}: no r1 file", i + 1)))?,
            r2: pega(2),
            condition: pega(3),
        })
    }
}

/// Quantifies RNA-seq samples straight from their FASTQ files, without an
/// aligner: pseudoalignment against an index made by `build_transcript_index`
/// and the EM of kallisto. Single-end or paired-end, plain or gzip.
///
/// - `samples`: a list with one item per sample: a dict {"name", "r1", "r2",
///   "condition"} or a tuple (name, r1[, r2[, condition]]). Without "r2" the
///   reads are single-end. A name left empty comes from the file name. When
///   every sample has a condition, `samples.tsv` is written too, ready for
///   `differential_expression`.
/// - `fragment_length`, `fragment_sd`: mean and standard deviation of the
///   fragment (insert) length. Required for single-end reads (kallisto's -l
///   and -s; 200 and 20 are common for Illumina mRNA libraries); for pairs
///   they are estimated from the first `fragment_pairs` pairs that fall on a
///   single transcript, unless given.
/// - `gene_map`: replaces the gene map stored in the index (GTF, GFF3 or a
///   transcript-gene table).
/// - `check_pair_names`: stop when R1 and R2 are out of step (read names
///   differ).
/// - `single_overhang`: single-end reads only. True (default) keeps every
///   compatible transcript; False applies kallisto's default rule (drops the
///   transcripts on which a fragment of the mean length, starting at the read,
///   would run past the end).
/// - `max_reads`: only the first reads (pairs) of each sample, for a quick look.
///
/// Writes in `output_dir`: a folder per sample (`abundance.tsv` as kallisto's,
/// `abundance_genes.tsv`, `run_info.json`) and the tables of all samples:
/// `gene_counts.tsv` (rounded estimated counts: the input of
/// `differential_expression`), `gene_tpm.tsv`, `transcript_counts.tsv`,
/// `transcript_tpm.tsv`, `gene_names.tsv`, `quant_stats.tsv`, `samples.tsv`.
/// Returns, per sample, the reads, pseudoaligned, unique, fragment length and
/// EM rounds; the gene ids, names, counts and TPM of every sample; the files
/// written; warnings; and `text`, the summary.
#[pyfunction]
#[pyo3(signature = (index, samples, output_dir, fragment_length=None, fragment_sd=None, gene_map=None, fragment_pairs=10000,
                    check_pair_names=true, single_overhang=true, max_reads=None))]
#[allow(clippy::too_many_arguments)]
fn quantify_samples(
    py: Python<'_>,
    index: &str,
    samples: &Bound<'_, PyAny>,
    output_dir: &str,
    fragment_length: Option<f64>,
    fragment_sd: Option<f64>,
    gene_map: Option<String>,
    fragment_pairs: usize,
    check_pair_names: bool,
    single_overhang: bool,
    max_reads: Option<u64>,
) -> PyResult<Py<PyAny>> {
    let mut amostras = Vec::new();
    for (i, item) in samples.try_iter()?.enumerate() {
        amostras.push(amostra_quant(&item?, i)?);
    }
    let op = quantificacao::QuantOptions {
        fragment_length,
        fragment_sd,
        fragment_pairs,
        gene_map,
        check_pair_names,
        single_overhang,
        max_reads,
        ..Default::default()
    };
    let r = py.allow_threads(|| quantificacao::quantify_samples(index, &amostras, output_dir, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// How one read (or pair) is pseudoaligned, step by step, for teaching: each
/// k-mer of the read, whether it is in the index and which transcripts
/// contain it, and the intersection -- the transcripts the read is compatible
/// with (`status`: "unique", "multi", "conflicting", "no_kmer" or
/// "too_short"). Loads the index on every call: meant for small indexes.
#[pyfunction]
#[pyo3(signature = (index, read, read2=None))]
fn explain_pseudoalignment(py: Python<'_>, index: &str, read: &str, read2: Option<&str>) -> PyResult<Py<PyAny>> {
    let j = py.allow_threads(|| quantificacao::explain_pseudoalignment(index, read, read2)).map_err(erro_py)?;
    json_para_py(py, &j)
}

/// Counts the guides of a pooled CRISPR screen in FASTQ files, as `mageck
/// count` does: each read's guide is found at its offset (found from the first
/// reads when not given) and matched to the library exactly (or with one
/// mismatch, when only one guide is that close).
///
/// - `library`: guide id, guide sequence, gene (TSV or CSV, header optional).
/// - `fastqs`: one file per sample, or a list of files per sample (technical
///   replicates, summed). `labels`: the sample names.
/// Returns per-sample statistics (reads, mapped %, guides with no read, Gini).
#[pyfunction]
#[pyo3(signature = (library, fastqs, output=None, labels=None, offsets=None, reverse_complement=false, mismatches=0))]
#[allow(clippy::too_many_arguments)]
fn crispr_count(
    py: Python<'_>,
    library: &str,
    fastqs: &Bound<'_, PyAny>,
    output: Option<&str>,
    labels: Option<&Bound<'_, PyAny>>,
    offsets: Option<Vec<usize>>,
    reverse_complement: bool,
    mismatches: u8,
) -> PyResult<Py<PyAny>> {
    let amostras: Vec<Vec<String>> = if let Ok(s) = fastqs.extract::<String>() {
        vec![vec![s]]
    } else if let Ok(v) = fastqs.extract::<Vec<String>>() {
        v.into_iter().map(|f| vec![f]).collect()
    } else {
        fastqs.extract::<Vec<Vec<String>>>().map_err(|_| {
            pyo3::exceptions::PyTypeError::new_err("`fastqs` must be a path, a list of paths, or a list of lists of paths")
        })?
    };
    let rotulos = lista_de_textos(labels, "labels")?;
    let op = crispr::CountOptions { offsets, reverse_complement, mismatches };
    let r = py.allow_threads(|| crispr::crispr_count(library, &amostras, &rotulos, output, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// Tests a CRISPR screen as `mageck test` (RRA) does: normalizes the counts,
/// models each guide's variance from its mean, scores each guide in the
/// treatment against the control, and ranks the genes by robust rank
/// aggregation for negative selection (depleted: essential genes) and positive
/// selection (enriched), with permutation p-values and FDR.
///
/// - `treatment`, `control`: sample names (or 0-based columns) of the count
///   table; without `control`, every other sample.
/// - `norm_method`: "median" (default), "total", "control" (with
///   `control_sgrnas` or `control_genes`) or "none".
/// - `output_prefix`: also write <prefix>.sgrna_summary.txt and
///   <prefix>.gene_summary.txt in MAGeCK's format.
/// Returns the summary, the top genes and, in "rank_plot", the ranked genes.
#[pyfunction]
#[pyo3(signature = (count_table, treatment, control=None, output_prefix=None, norm_method="median", control_sgrnas=None,
                    control_genes=None, gene_lfc="median", permutations=None, seed=123456))]
#[allow(clippy::too_many_arguments)]
fn crispr_test(
    py: Python<'_>,
    count_table: &str,
    treatment: &Bound<'_, PyAny>,
    control: Option<&Bound<'_, PyAny>>,
    output_prefix: Option<&str>,
    norm_method: &str,
    control_sgrnas: Option<String>,
    control_genes: Option<String>,
    gene_lfc: &str,
    permutations: Option<usize>,
    seed: u32,
) -> PyResult<Py<PyAny>> {
    let tratamento = lista_de_textos(Some(treatment), "treatment")?;
    let controle = lista_de_textos(control, "control")?;
    let mut op = crispr::TestOptions {
        norm_method: crispr::NormMethod::from_text(norm_method).map_err(erro_py)?,
        gene_lfc: crispr::GeneLfc::from_text(gene_lfc).map_err(erro_py)?,
        permutations,
        seed,
        ..Default::default()
    };
    if let Some(p) = control_sgrnas {
        op.control_sgrnas = py.allow_threads(|| crispr::read_name_list(&p)).map_err(erro_py)?;
    }
    if let Some(p) = control_genes {
        op.control_genes = py.allow_threads(|| crispr::read_name_list(&p)).map_err(erro_py)?;
    }
    let r = py
        .allow_threads(|| crispr::crispr_test(count_table, &tratamento, &controle, output_prefix, &op))
        .map_err(erro_py)?;
    json_para_py(py, &r.json_compact())
}

/// Sketches reference genomes for species identification: the FracMinHash
/// sketch of each genome (its canonical k-mers whose hash falls below
/// 2^64/scaled), with the same hashes as sourmash, saved as one database file.
///
/// - `genomes`: one entry per genome: a FASTA path, or a list of paths (one
///   genome in several files). `names`: their names (default: the first
///   FASTA header). `host`: True for the genomes that are host (e.g. human).
#[pyfunction]
#[pyo3(signature = (genomes, output, names=None, host=None, k=31, scaled=1000))]
fn sketch_references(
    py: Python<'_>,
    genomes: &Bound<'_, PyAny>,
    output: &str,
    names: Option<&Bound<'_, PyAny>>,
    host: Option<Vec<bool>>,
    k: u8,
    scaled: u64,
) -> PyResult<Py<PyAny>> {
    let listas: Vec<Vec<String>> = if let Ok(s) = genomes.extract::<String>() {
        vec![vec![s]]
    } else if let Ok(v) = genomes.extract::<Vec<String>>() {
        v.into_iter().map(|f| vec![f]).collect()
    } else {
        genomes.extract::<Vec<Vec<String>>>().map_err(|_| {
            pyo3::exceptions::PyTypeError::new_err("`genomes` must be a path, a list of paths, or a list of lists of paths")
        })?
    };
    let nomes_g = opcoes_de_lista(names, "names", listas.len())?;
    let hospede = host.unwrap_or_else(|| vec![false; listas.len()]);
    if hospede.len() != listas.len() {
        return Err(pyo3::exceptions::PyValueError::new_err("`host` needs one True/False per genome"));
    }
    let entradas: Vec<microbio::GenomeInput> = listas
        .into_iter()
        .zip(nomes_g)
        .zip(hospede)
        .map(|((files, name), host)| microbio::GenomeInput { files, name: name.unwrap_or_default(), host })
        .collect();
    let op = microbio::SketchOptions { k, scaled };
    let r = py.allow_threads(|| microbio::sketch_references(&entradas, output, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// Which species are in a sample, and how much of each, as `sourmash gather`
/// does: the sample (reads or an assembly; both mates of a pair) is sketched
/// and covered greedily by the references of the databases
/// (`sketch_references`), the largest overlap first. Each match has its
/// containment, an ANI estimate and its relative abundance; host genomes
/// are reported apart. Experimental, for teaching and research.
#[pyfunction]
#[pyo3(signature = (samples, databases, output=None, scaled=None, threshold_bp=50000, min_abundance=1))]
fn gather(
    py: Python<'_>,
    samples: &Bound<'_, PyAny>,
    databases: &Bound<'_, PyAny>,
    output: Option<String>,
    scaled: Option<u64>,
    threshold_bp: u64,
    min_abundance: u32,
) -> PyResult<Py<PyAny>> {
    let amostras = lista_de_textos(Some(samples), "samples")?;
    let bancos = lista_de_textos(Some(databases), "databases")?;
    let op = microbio::GatherOptions { scaled, threshold_bp, min_abundance, output };
    let r = py.allow_threads(|| microbio::gather(&amostras, &bancos, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// Classifies each read (or pair) to the reference genome it shares the most
/// sketch k-mers with, as Kraken does without a taxonomy: an abundance table
/// (reads per reference), host reads counted apart and, optionally, the
/// host-free reads written out. Use a database sketched with a small
/// `scaled` (e.g. 10) for reads. Experimental, for teaching and research.
#[pyfunction]
#[pyo3(signature = (input, databases, input2=None, output=None, scaled=None, min_hits=2, min_fraction=0.0, only=None,
                    assignments=None, host_free=None, host_free2=None))]
#[allow(clippy::too_many_arguments)]
fn classify_reads(
    py: Python<'_>,
    input: &str,
    databases: &Bound<'_, PyAny>,
    input2: Option<&str>,
    output: Option<String>,
    scaled: Option<u64>,
    min_hits: u32,
    min_fraction: f64,
    only: Option<&Bound<'_, PyAny>>,
    assignments: Option<String>,
    host_free: Option<String>,
    host_free2: Option<String>,
) -> PyResult<Py<PyAny>> {
    let bancos = lista_de_textos(Some(databases), "databases")?;
    let op = microbio::ClassifyOptions {
        scaled,
        min_hits,
        min_fraction,
        only: lista_de_textos(only, "only")?,
        output,
        assignments,
        host_free,
        host_free2,
    };
    let r = py.allow_threads(|| microbio::classify_reads(input, input2, &bancos, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// The known antimicrobial-resistance genes of an isolate: its assembly
/// (contigs) or its reads are searched for each gene of the NCBI
/// AMRFinderPlus catalog (AMR_CDS.fa, with ReferenceGeneCatalog.txt for the
/// drug class), reporting the best allele per locus with its identity and
/// coverage. Experimental, for teaching and research: it does not replace
/// susceptibility testing.
///
/// - `mode`: "auto" (FASTQ = reads, FASTA = assembly), "assembly" or "reads".
/// - `include_plus`: also the "plus" genes (virulence, stress, biocides).
#[pyfunction]
#[pyo3(signature = (inputs, catalog, catalog_table=None, output=None, min_identity=0.9, min_coverage=0.9,
                    min_partial_coverage=0.5, include_plus=true, mode="auto", min_depth=1))]
#[allow(clippy::too_many_arguments)]
fn find_resistance_genes(
    py: Python<'_>,
    inputs: &Bound<'_, PyAny>,
    catalog: &str,
    catalog_table: Option<&str>,
    output: Option<String>,
    min_identity: f64,
    min_coverage: f64,
    min_partial_coverage: f64,
    include_plus: bool,
    mode: &str,
    min_depth: u32,
) -> PyResult<Py<PyAny>> {
    let entradas = lista_de_textos(Some(inputs), "inputs")?;
    let modo = match mode.to_ascii_lowercase().as_str() {
        "auto" => microbio::AmrMode::Auto,
        "assembly" => microbio::AmrMode::Assembly,
        "reads" => microbio::AmrMode::Reads,
        outro => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "mode `{outro}`: use \"auto\", \"assembly\" or \"reads\""
            )))
        }
    };
    let op = microbio::AmrOptions { min_identity, min_coverage, min_partial_coverage, include_plus, mode: modo, min_depth, output };
    let r = py
        .allow_threads(|| microbio::find_resistance_genes(&entradas, catalog, catalog_table, &op))
        .map_err(erro_py)?;
    json_para_py(py, &r.json())
}

// ---------------------------------------------------------
// Surto e vigilancia (src/surto.rs)
// ---------------------------------------------------------

fn opcoes_de_distancia(
    k: u8,
    scaled: u64,
    records_as_genomes: bool,
    method: &str,
    outgroup: Option<String>,
    matrix_output: Option<String>,
    matrix_format: Option<&str>,
    matrix_value: &str,
    tree_output: Option<String>,
) -> PyResult<surto::DistanceOptions> {
    Ok(surto::DistanceOptions {
        k,
        scaled,
        records_as_genomes,
        method: surto::TreeMethod::from_text(method).map_err(erro_py)?,
        outgroup,
        matrix_output,
        matrix_format: matrix_format.map(surto::MatrixFormat::from_text).transpose().map_err(erro_py)?,
        matrix_value: surto::MatrixValue::from_text(matrix_value).map_err(erro_py)?,
        tree_output,
    })
}

/// Distances between genomes and the tree they make, for an outbreak or a
/// collection of isolates: each genome (FASTA; a list of files when it is in
/// several, e.g. chromosome and plasmids) is reduced to a FracMinHash sketch
/// (the same hashes as sourmash), and each pair gets its Jaccard index, the
/// Mash distance and an ANI estimate. The tree is built by neighbor joining
/// (rooted at the midpoint, or on `outgroup`) or UPGMA.
///
/// - `genomes`: a path, a list of paths (one genome each), or a list of lists.
/// - `names`: the names in the matrix and the tree (default: the file names).
/// - `records_as_genomes`: each FASTA record is a genome (a multi-FASTA of
///   viral genomes).
/// - `scaled`: one k-mer in `scaled` is kept; 0 = chosen from the sizes.
/// - `method`: "nj" or "upgma". `matrix_value`: "distance", "ani" or
///   "jaccard"; `matrix_format`: "tsv" or "phylip" (default: by the name).
/// - `matrix_output`, `tree_output` (Newick): files to write.
///
/// The result has `names`, the `jaccard`, `distance` and `ani` matrices and
/// `tree` (Newick, and nodes with x = distance from the root and y = leaf
/// order, ready to draw).
#[pyfunction]
#[pyo3(signature = (genomes, names=None, tree_output=None, matrix_output=None, k=31, scaled=0, records_as_genomes=false,
                    method="nj", outgroup=None, matrix_value="distance", matrix_format=None))]
#[allow(clippy::too_many_arguments)]
fn genome_distances(
    py: Python<'_>,
    genomes: &Bound<'_, PyAny>,
    names: Option<&Bound<'_, PyAny>>,
    tree_output: Option<String>,
    matrix_output: Option<String>,
    k: u8,
    scaled: u64,
    records_as_genomes: bool,
    method: &str,
    outgroup: Option<String>,
    matrix_value: &str,
    matrix_format: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let listas: Vec<Vec<String>> = if let Ok(s) = genomes.extract::<String>() {
        vec![vec![s]]
    } else if let Ok(v) = genomes.extract::<Vec<String>>() {
        v.into_iter().map(|f| vec![f]).collect()
    } else {
        genomes.extract::<Vec<Vec<String>>>().map_err(|_| {
            pyo3::exceptions::PyTypeError::new_err("`genomes` must be a path, a list of paths, or a list of lists of paths")
        })?
    };
    let nomes = opcoes_de_lista(names, "names", listas.len())?;
    let entradas: Vec<surto::GenomeFiles> =
        listas.into_iter().zip(nomes).map(|(files, name)| surto::GenomeFiles { files, name: name.unwrap_or_default() }).collect();
    let op = opcoes_de_distancia(k, scaled, records_as_genomes, method, outgroup, matrix_output, matrix_format, matrix_value, tree_output)?;
    let r = py.allow_threads(|| surto::genome_distances(&entradas, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// The same as `genome_distances`, from sketch databases already made
/// (`sketch_references`, or sourmash signatures of the same k): every genome
/// of the databases goes into the matrix and the tree.
#[pyfunction]
#[pyo3(signature = (databases, tree_output=None, matrix_output=None, scaled=0, method="nj", outgroup=None,
                    matrix_value="distance", matrix_format=None))]
#[allow(clippy::too_many_arguments)]
fn sketch_distances(
    py: Python<'_>,
    databases: &Bound<'_, PyAny>,
    tree_output: Option<String>,
    matrix_output: Option<String>,
    scaled: u64,
    method: &str,
    outgroup: Option<String>,
    matrix_value: &str,
    matrix_format: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let bases = lista_de_textos(Some(databases), "databases")?;
    let op = opcoes_de_distancia(31, scaled, false, method, outgroup, matrix_output, matrix_format, matrix_value, tree_output)?;
    let r = py.allow_threads(|| surto::sketch_distances(&bases, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// A tree from a distance (or similarity) matrix made elsewhere: TSV, CSV
/// (as `sourmash compare --csv` writes it) or PHYLIP. `similarity=True` for a
/// matrix of similarities (1 = identical), which is turned into 1 - value.
#[pyfunction]
#[pyo3(signature = (path, method="nj", similarity=false, outgroup=None, output=None))]
fn tree_from_matrix(
    py: Python<'_>,
    path: &str,
    method: &str,
    similarity: bool,
    outgroup: Option<&str>,
    output: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let m = surto::TreeMethod::from_text(method).map_err(erro_py)?;
    let r = py.allow_threads(|| surto::tree_from_matrix(path, m, similarity, outgroup, output)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// The mutations of viral genomes against their reference, named as
/// Nextclade names them: nucleotide substitutions (`C241T`), deletions
/// (`21765-21770`), insertions (`22204:GAGCCAGAA`) and missing (N) stretches,
/// in reference coordinates; and, with the reference's GFF3, the amino-acid
/// changes they cause in each gene (`S:D614G`, `S:H69-`, `ORF1b:P314L`),
/// with multi-part CDS such as ORF1ab and its ribosomal frameshift.
///
/// - `genomes`: FASTA files; every record is a genome.
/// - `output_dir`: one TSV per genome; `summary_output`: one line per genome.
/// - `split_orf1ab`: report ORF1ab as ORF1a + ORF1b, as Nextclade does.
#[pyfunction]
#[pyo3(signature = (reference, genomes, annotation=None, output_dir=None, summary_output=None, split_orf1ab=true))]
fn viral_mutations(
    py: Python<'_>,
    reference: &str,
    genomes: &Bound<'_, PyAny>,
    annotation: Option<String>,
    output_dir: Option<String>,
    summary_output: Option<String>,
    split_orf1ab: bool,
) -> PyResult<Py<PyAny>> {
    let arquivos = lista_de_textos(Some(genomes), "genomes")?;
    let op = surto::MutationOptions { annotation, output_dir, summary_output, split_orf1ab };
    let r = py.allow_threads(|| surto::viral_mutations(reference, &arquivos, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// In-silico PCR: the primer pairs of a panel (a TSV with `target`,
/// `forward`, `reverse`, optional `probe`, `size` such as "60-90", `source`;
/// IUPAC codes allowed) against genomes (FASTA: the products, on both
/// strands, with their mismatches) and read samples (FASTQ: reads carrying a
/// primer site). A primer may have up to `max_mismatches`, but none in its
/// last `three_prime_exact` bases; a product is detected when its size is in
/// the assay's range and its probe (if any) binds.
///
/// Returns `rows` (the presence/absence matrix, one row per genome or
/// sample), `amplicons` and `read_counts`; writes them as TSV when asked.
#[pyfunction]
#[pyo3(signature = (panel, inputs, amplicons_output=None, matrix_output=None, reads_output=None, max_mismatches=3,
                    three_prime_exact=3, max_total_mismatches=4, max_probe_mismatches=3, max_amplicon=2000,
                    records_as_genomes=false, min_reads=3))]
#[allow(clippy::too_many_arguments)]
fn in_silico_pcr(
    py: Python<'_>,
    panel: &str,
    inputs: &Bound<'_, PyAny>,
    amplicons_output: Option<String>,
    matrix_output: Option<String>,
    reads_output: Option<String>,
    max_mismatches: u32,
    three_prime_exact: u32,
    max_total_mismatches: u32,
    max_probe_mismatches: u32,
    max_amplicon: u32,
    records_as_genomes: bool,
    min_reads: u64,
) -> PyResult<Py<PyAny>> {
    let entradas = lista_de_textos(Some(inputs), "inputs")?;
    let op = surto::PcrOptions {
        max_mismatches,
        three_prime_exact,
        max_total_mismatches,
        max_probe_mismatches,
        max_amplicon,
        records_as_genomes,
        min_reads,
        amplicons_output,
        matrix_output,
        reads_output,
    };
    let r = py.allow_threads(|| surto::in_silico_pcr_panel(panel, &entradas, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

// ---------------------------------------------------------
// 16S e carga mutacional do tumor (src/amplicon.rs)
// ---------------------------------------------------------

fn amostras_16s(samples: &Bound<'_, PyAny>) -> PyResult<Vec<amplicon::AmpliconSample>> {
    if let Ok(s) = samples.extract::<String>() {
        let baixo = s.to_ascii_lowercase();
        let fastq = [".fastq", ".fq", ".fastq.gz", ".fq.gz"].iter().any(|e| baixo.ends_with(e));
        return if fastq { amplicon::pair_fastq_files(&[s]) } else { amplicon::read_sample_sheet(&s) }.map_err(erro_py);
    }
    if let Ok(v) = samples.extract::<Vec<String>>() {
        return amplicon::pair_fastq_files(&v).map_err(erro_py);
    }
    let linhas: Vec<Vec<String>> = samples.extract().map_err(|_| {
        pyo3::exceptions::PyTypeError::new_err(
            "`samples` must be a sample sheet (name, R1, R2 per line, as mothur's stability.files), a list of FASTQ files (paired by their _R1/_R2 names) or a list of [name, R1, R2]",
        )
    })?;
    linhas
        .into_iter()
        .map(|l| match l.as_slice() {
            [n, r1] => Ok(amplicon::AmpliconSample { name: n.clone(), r1: r1.clone(), r2: None }),
            [n, r1, r2] => Ok(amplicon::AmpliconSample { name: n.clone(), r1: r1.clone(), r2: Some(r2.clone()) }),
            _ => Err(pyo3::exceptions::PyValueError::new_err("each sample is [name, R1] or [name, R1, R2]")),
        })
        .collect()
}

/// Amplicon sequence variants (ASVs) from 16S reads, as vsearch / UNOISE3
/// do: each pair merged by its overlap, primers removed (when given),
/// reads with more than `max_ee` expected errors dropped, the reads of all
/// samples dereplicated and denoised (UNOISE3: alpha 2, at least `min_size`
/// copies), chimeras removed (as `--uchime3_denovo`). Writes the ASVs
/// (FASTA, ASV_1 the most abundant) and the feature table (ASV x sample).
///
/// - `samples`: a sample sheet (name, R1, R2 per line; mothur's
///   `stability.files` works as it is), a list of FASTQ files paired by
///   their `_R1`/`_R2` names, or a list of `[name, R1, R2]`.
/// - `tracking`: also write the reads kept at each step, per sample.
#[pyfunction]
#[pyo3(signature = (samples, asv_fasta, feature_table, tracking=None, primer_forward=None, primer_reverse=None,
                    require_primers=true, max_ee=1.0, min_length=1, max_length=0, min_overlap=10, max_diffs=10,
                    trunc_len_r1=0, trunc_len_r2=0, min_size=8, alpha=2.0, remove_chimeras=true))]
#[allow(clippy::too_many_arguments)]
fn amplicon_asvs(
    py: Python<'_>,
    samples: &Bound<'_, PyAny>,
    asv_fasta: &str,
    feature_table: &str,
    tracking: Option<&str>,
    primer_forward: Option<String>,
    primer_reverse: Option<String>,
    require_primers: bool,
    max_ee: Option<f64>,
    min_length: usize,
    max_length: usize,
    min_overlap: usize,
    max_diffs: usize,
    trunc_len_r1: usize,
    trunc_len_r2: usize,
    min_size: u64,
    alpha: f64,
    remove_chimeras: bool,
) -> PyResult<Py<PyAny>> {
    let amostras = amostras_16s(samples)?;
    let mut op = amplicon::AsvOptions { primer_forward, primer_reverse, require_primers, max_ee, min_length, max_length, ..Default::default() };
    op.merge.min_overlap = min_overlap;
    op.merge.max_diffs = max_diffs;
    op.merge.trunc_len_r1 = trunc_len_r1;
    op.merge.trunc_len_r2 = trunc_len_r2;
    op.denoise.min_size = min_size;
    op.denoise.alpha = alpha;
    op.chimera.enabled = remove_chimeras;
    let r = py
        .allow_threads(|| amplicon::amplicon_asvs(&amostras, asv_fasta, feature_table, tracking, &op))
        .map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// Trains the naive Bayesian classifier of the RDP (Wang et al. 2007) on a
/// reference FASTA whose headers carry the lineage (the DADA2-formatted RDP
/// and GTDB files, or GTDB's own `d__...;g__...` headers), and saves the
/// model to reload fast. With both primers, only the amplified region of
/// each reference is used, which suits short reads (V4) better.
#[pyfunction]
#[pyo3(signature = (reference, output=None, k=8, primer_forward=None, primer_reverse=None, rank_names=None, max_depth=None))]
#[allow(clippy::too_many_arguments)]
fn train_classifier(
    py: Python<'_>,
    reference: &str,
    output: Option<&str>,
    k: usize,
    primer_forward: Option<String>,
    primer_reverse: Option<String>,
    rank_names: Option<Vec<String>>,
    max_depth: Option<usize>,
) -> PyResult<Py<PyAny>> {
    let op = amplicon::TrainOptions { k, max_depth, rank_names, primer_forward, primer_reverse };
    let r = py.allow_threads(|| amplicon::train_classifier(reference, output, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// The taxonomy of each sequence (ASVs) by the RDP naive Bayesian
/// classifier: 8-mer words, 100 bootstrap trials of 1/8 of the words, and
/// the deepest rank whose confidence reaches `min_confidence` (0.8 as RDP;
/// 0.5 as DADA2's minBoot for short reads). `model`: a model saved by
/// `train_classifier`, or a reference FASTA (trained on the fly).
/// `species_reference`: a FASTA of species for exact matches (DADA2's
/// addSpecies).
#[pyfunction]
#[pyo3(signature = (sequences, model, output, min_confidence=0.8, bootstraps=100, seed=42, try_reverse_complement=false,
                    species_reference=None))]
#[allow(clippy::too_many_arguments)]
fn classify_asvs(
    py: Python<'_>,
    sequences: &str,
    model: &str,
    output: &str,
    min_confidence: f64,
    bootstraps: u32,
    seed: u64,
    try_reverse_complement: bool,
    species_reference: Option<String>,
) -> PyResult<Py<PyAny>> {
    let op = amplicon::ClassifyOptions { bootstraps, min_confidence, seed, try_reverse_complement, species_reference, ..Default::default() };
    let r = py.allow_threads(|| amplicon::classify_asvs(sequences, model, output, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// Alpha diversity of each sample of a feature table: observed, Chao1,
/// Shannon, Simpson, Pielou's evenness, Good's coverage; optionally after
/// rarefying to `rarefy_depth`, and rarefaction curves (the exact
/// expectation at each depth).
#[pyfunction]
#[pyo3(signature = (table, output=None, rarefy_depth=None, seed=42, curve_steps=20, curve_max_depth=None, curves_output=None))]
#[allow(clippy::too_many_arguments)]
fn alpha_diversity(
    py: Python<'_>,
    table: &str,
    output: Option<&str>,
    rarefy_depth: Option<u64>,
    seed: u64,
    curve_steps: usize,
    curve_max_depth: Option<u64>,
    curves_output: Option<String>,
) -> PyResult<Py<PyAny>> {
    let op = amplicon::AlphaOptions { rarefy_depth, seed, curve_steps, curve_max_depth, curves_output };
    let r = py.allow_threads(|| amplicon::alpha_diversity(table, output, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// Beta diversity between the samples of a feature table: Bray-Curtis,
/// Jaccard, or "kmer" (a distance between the samples' ASV k-mer profiles,
/// with `asv_fasta`), the principal coordinates (PCoA) and, with a
/// `metadata` table and its `group_column`, PERMANOVA (permutations, fixed
/// seed).
#[pyfunction]
#[pyo3(signature = (table, output=None, metric="braycurtis", relative=false, rarefy_depth=None, seed=42, asv_fasta=None,
                    metadata=None, group_column=None, permutations=999, pcoa_output=None))]
#[allow(clippy::too_many_arguments)]
fn beta_diversity(
    py: Python<'_>,
    table: &str,
    output: Option<&str>,
    metric: &str,
    relative: bool,
    rarefy_depth: Option<u64>,
    seed: u64,
    asv_fasta: Option<String>,
    metadata: Option<String>,
    group_column: Option<String>,
    permutations: u32,
    pcoa_output: Option<String>,
) -> PyResult<Py<PyAny>> {
    let op = amplicon::BetaOptions {
        metric: metric.to_string(),
        relative,
        rarefy_depth,
        seed,
        asv_fasta,
        metadata,
        group_column,
        permutations,
        pcoa_output,
    };
    let r = py.allow_threads(|| amplicon::beta_diversity(table, output, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// The relative abundance of each taxon (at `rank`: "Phylum", "Family",
/// "Genus"...) in each sample, from a feature table and the taxonomy of its
/// ASVs (`classify_asvs`), for stacked bars; `top` keeps the most abundant
/// taxa and adds the rest as "Other".
#[pyfunction]
#[pyo3(signature = (table, taxonomy, output=None, rank="Phylum", top=0))]
fn taxa_composition(py: Python<'_>, table: &str, taxonomy: &str, output: Option<&str>, rank: &str, top: usize) -> PyResult<Py<PyAny>> {
    let op = amplicon::CompositionOptions { rank: rank.to_string(), top };
    let r = py.allow_threads(|| amplicon::taxa_composition(table, taxonomy, output, &op)).map_err(erro_py)?;
    json_para_py(py, &r.json())
}

/// Tumor mutational burden: the somatic variants of a tumor VCF counted per
/// megabase of the sequenced region (`bed`: the panel's or exome's callable
/// region). Filters: PASS only, minimum depth and variant allele fraction
/// (from FORMAT AD/DP or AF); optionally only coding variants (`annotation`,
/// a GFF3) and only non-synonymous ones (also `reference`, the genome
/// FASTA). Experimental, for teaching and research: not a clinical report.
#[pyfunction]
#[pyo3(signature = (vcf, bed, output=None, tumor_sample=None, pass_only=true, min_depth=10, min_vaf=0.05, annotation=None,
                    coding_only=false, reference=None, nonsynonymous_only=false))]
#[allow(clippy::too_many_arguments)]
fn tumor_mutational_burden(
    py: Python<'_>,
    vcf: &str,
    bed: &str,
    output: Option<String>,
    tumor_sample: Option<String>,
    pass_only: bool,
    min_depth: u32,
    min_vaf: f64,
    annotation: Option<String>,
    coding_only: bool,
    reference: Option<String>,
    nonsynonymous_only: bool,
) -> PyResult<Py<PyAny>> {
    let op = amplicon::TmbOptions { tumor_sample, pass_only, min_depth, min_vaf, annotation, coding_only, reference, nonsynonymous_only, output };
    let r = py.allow_threads(|| amplicon::tumor_mutational_burden(vcf, bed, &op)).map_err(erro_py)?;
    let d = json_para_py(py, &r.json())?;
    d.bind(py).set_item("notice", AVISO)?;
    Ok(d)
}

/// MultiMerge Genomics: tools for genomic files of any size. Validate,
/// normalize chromosome names, compare reference genomes, sort, index and
/// query by region; interval operations on BED; statistics, filtering and
/// normalization of VCF; quality report and trimming of FASTQ; k-mer
/// counting; plain text, BGZF, BAM and Parquet in and out; and k-mer modules
/// with optional CUDA acceleration. Start with `engine_info()`; every
/// function has a docstring.
//
// O módulo compilado é `multimerge_genomics._native`; o pacote Python
// (python/multimerge_genomics/__init__.py) o reexporta inteiro e acrescenta a
// parte escrita em Python (o assistente de IA e o servidor MCP).
#[pymodule]
#[pyo3(name = "_native")]
fn multimerge_genomics(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_and_merge_kmer_files, m)?)?;
    m.add_function(wrap_pyfunction!(decode_kmers, m)?)?;
    m.add_function(wrap_pyfunction!(minhash_signatures, m)?)?;
    m.add_function(wrap_pyfunction!(read_fastq_batch, m)?)?;
    m.add_function(wrap_pyfunction!(read_fasta_batch, m)?)?;
    m.add_function(wrap_pyfunction!(sort_merge_bytes, m)?)?;
    m.add_function(wrap_pyfunction!(build_debruijn_graph, m)?)?;
    m.add_function(wrap_pyfunction!(assemble_debruijn, m)?)?;
    m.add_function(wrap_pyfunction!(build_bloom_filter, m)?)?;
    m.add_function(wrap_pyfunction!(extract_mutations, m)?)?;
    m.add_function(wrap_pyfunction!(translate_proteins, m)?)?;
    m.add_function(wrap_pyfunction!(build_scaffolds, m)?)?;
    m.add_function(wrap_pyfunction!(test_giant_kmer, m)?)?;
    m.add_function(wrap_pyfunction!(decode_kmers_128, m)?)?;
    m.add_function(wrap_pyfunction!(multimerge_vcf, m)?)?;
    m.add_function(wrap_pyfunction!(sort_genomic, m)?)?;
    m.add_function(wrap_pyfunction!(bed_merge, m)?)?;
    m.add_function(wrap_pyfunction!(bed_intersect, m)?)?;
    m.add_function(wrap_pyfunction!(validate_file, m)?)?;
    m.add_function(wrap_pyfunction!(normalize_chromosomes, m)?)?;
    m.add_function(wrap_pyfunction!(file_reference, m)?)?;
    m.add_function(wrap_pyfunction!(compare_references, m)?)?;
    m.add_function(wrap_pyfunction!(engine_info, m)?)?;
    m.add_function(wrap_pyfunction!(bed_subtract, m)?)?;
    m.add_function(wrap_pyfunction!(bed_complement, m)?)?;
    m.add_function(wrap_pyfunction!(bed_closest, m)?)?;
    m.add_function(wrap_pyfunction!(bed_coverage, m)?)?;
    m.add_function(wrap_pyfunction!(bed_slop, m)?)?;
    m.add_function(wrap_pyfunction!(vcf_stats, m)?)?;
    m.add_function(wrap_pyfunction!(vcf_filter, m)?)?;
    m.add_function(wrap_pyfunction!(vcf_normalize, m)?)?;
    m.add_function(wrap_pyfunction!(fastq_stats, m)?)?;
    m.add_function(wrap_pyfunction!(fastq_report, m)?)?;
    m.add_function(wrap_pyfunction!(fastq_trim, m)?)?;
    m.add_function(wrap_pyfunction!(count_kmers, m)?)?;
    m.add_function(wrap_pyfunction!(convert_file, m)?)?;
    m.add_function(wrap_pyfunction!(parquet_options, m)?)?;
    m.add_function(wrap_pyfunction!(index_file, m)?)?;
    m.add_function(wrap_pyfunction!(query_region, m)?)?;
    m.add_function(wrap_pyfunction!(clinvar_findings, m)?)?;
    m.add_function(wrap_pyfunction!(gene_regions, m)?)?;
    m.add_function(wrap_pyfunction!(variant_lookup, m)?)?;
    m.add_function(wrap_pyfunction!(pileup, m)?)?;
    m.add_function(wrap_pyfunction!(pharmcat_prepare, m)?)?;
    m.add_function(wrap_pyfunction!(pharmcat_run, m)?)?;
    m.add_function(wrap_pyfunction!(mmgen_pgx_run, m)?)?;
    m.add_function(wrap_pyfunction!(pgx_run, m)?)?;
    m.add_function(wrap_pyfunction!(pgx_cohort, m)?)?;
    m.add_function(wrap_pyfunction!(relatedness, m)?)?;
    m.add_function(wrap_pyfunction!(trio_analysis, m)?)?;
    m.add_function(wrap_pyfunction!(runs_of_homozygosity, m)?)?;
    m.add_function(wrap_pyfunction!(genetic_sex, m)?)?;
    m.add_function(wrap_pyfunction!(build_panel, m)?)?;
    m.add_function(wrap_pyfunction!(panel_from_frequencies, m)?)?;
    m.add_function(wrap_pyfunction!(estimate_ancestry, m)?)?;
    m.add_function(wrap_pyfunction!(polygenic_score, m)?)?;
    m.add_function(wrap_pyfunction!(count_reads, m)?)?;
    m.add_function(wrap_pyfunction!(differential_expression, m)?)?;
    m.add_function(wrap_pyfunction!(gene_set_test, m)?)?;
    m.add_function(wrap_pyfunction!(genotype_association, m)?)?;
    m.add_function(wrap_pyfunction!(association_power, m)?)?;
    m.add_function(wrap_pyfunction!(build_screen_database, m)?)?;
    m.add_function(wrap_pyfunction!(describe_screen_database, m)?)?;
    m.add_function(wrap_pyfunction!(contaminant_screen, m)?)?;
    m.add_function(wrap_pyfunction!(build_transcript_index, m)?)?;
    m.add_function(wrap_pyfunction!(transcript_index_info, m)?)?;
    m.add_function(wrap_pyfunction!(quantify_samples, m)?)?;
    m.add_function(wrap_pyfunction!(explain_pseudoalignment, m)?)?;
    m.add_function(wrap_pyfunction!(crispr_count, m)?)?;
    m.add_function(wrap_pyfunction!(crispr_test, m)?)?;
    m.add_function(wrap_pyfunction!(sketch_references, m)?)?;
    m.add_function(wrap_pyfunction!(gather, m)?)?;
    m.add_function(wrap_pyfunction!(classify_reads, m)?)?;
    m.add_function(wrap_pyfunction!(find_resistance_genes, m)?)?;
    m.add_function(wrap_pyfunction!(genome_distances, m)?)?;
    m.add_function(wrap_pyfunction!(sketch_distances, m)?)?;
    m.add_function(wrap_pyfunction!(tree_from_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(viral_mutations, m)?)?;
    m.add_function(wrap_pyfunction!(in_silico_pcr, m)?)?;
    m.add_function(wrap_pyfunction!(amplicon_asvs, m)?)?;
    m.add_function(wrap_pyfunction!(train_classifier, m)?)?;
    m.add_function(wrap_pyfunction!(classify_asvs, m)?)?;
    m.add_function(wrap_pyfunction!(alpha_diversity, m)?)?;
    m.add_function(wrap_pyfunction!(beta_diversity, m)?)?;
    m.add_function(wrap_pyfunction!(taxa_composition, m)?)?;
    m.add_function(wrap_pyfunction!(tumor_mutational_burden, m)?)?;
    Ok(())
}
