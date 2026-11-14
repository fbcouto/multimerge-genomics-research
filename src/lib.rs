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
    let mut kmer_a: Vec<u64> = como_u64(bytes_a).into_owned();
    let mut kmer_b: Vec<u64> = como_u64(bytes_b).into_owned();

    multimerge_core::multi_merge_sort_stable(&mut kmer_a);
    multimerge_core::multi_merge_sort_stable(&mut kmer_b);

    let resultado = ordenar_e_mesclar_backend(&kmer_a, &kmer_b, onde)?;

    let byte_slice: &[u8] = unsafe {
        std::slice::from_raw_parts(resultado.as_ptr() as *const u8, resultado.len() * 8)
    };
    Ok(PyBytes::new(py, byte_slice))
}

// ---------------------------------------------------------
// MÓDULO 1: QC E LIMPEZA FASTQ/FASTA EM LOTE (BATCH U128)
// ---------------------------------------------------------
/// Reads one or more FASTQ files (plain text or `.gz`) and returns their
/// distinct k-mers, sorted, as raw bytes: 16 bytes (one unsigned 128-bit
/// integer) per k-mer.
///
/// Reads whose mean Phred quality is below `min_quality` are dropped. `k` goes
/// up to 63. The files are read in parallel.
#[pyfunction]
fn read_fastq_batch<'py>(py: Python<'py>, paths: Vec<String>, min_quality: u8, k: u8) -> PyResult<Bound<'py, PyBytes>> {
    let kmers = fastq::processar_fastq_lote(paths, min_quality, k);

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
    let kmers = fasta::processar_fasta_lote(paths, k);

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
fn build_debruijn_graph(kmer_bytes: &[u8], k: u8) -> PyResult<(usize, usize)> {
    let kmers = como_u128(kmer_bytes);
    Ok(debruijn::construir_grafo_estatisticas(&kmers, k))
}

/// Assembles contigs from a set of k-mers by walking the unbranched paths of
/// the De Bruijn graph. Returns the contigs as strings, longest first.
#[pyfunction]
fn assemble_debruijn(kmer_bytes: &[u8], k: u8) -> PyResult<Vec<String>> {
    let kmers = como_u128(kmer_bytes);
    Ok(debruijn::montar_todos_contigs(&kmers, k))
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
    let kmers_128 = como_u128(kmer_bytes);
    let kmers_gpu: Vec<u64> = kmers_128.par_iter().map(|&k| bloom::dobra(k)).collect();
    let bloom_bytes = size_mb * 1024 * 1024;

    let bloom_array = despacha(
        onde,
        "Bloom filter",
        || gpu::bloom_construir(&kmers_gpu, bloom_bytes, num_hashes),
        || bloom::construir(&kmers_gpu, bloom_bytes, num_hashes),
    )?;
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
    let kmers_128 = como_u128(sample_bytes);
    let kmers_gpu: Vec<u64> = kmers_128.par_iter().map(|&k| bloom::dobra(k)).collect();

    let mascara_resultados = despacha(
        onde,
        "Bloom filter query",
        || gpu::bloom_consultar(&kmers_gpu, bloom_bytes, num_hashes),
        || bloom::consultar(&kmers_gpu, bloom_bytes, num_hashes),
    )?;

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
/// Translates DNA sequences into proteins (standard genetic code, reading
/// frame 1; `*` marks a stop codon). One protein per input sequence.
#[pyfunction]
fn translate_proteins(dna_sequences: Vec<String>) -> PyResult<Vec<String>> {
    println!("⚙️ RAYON: translating {} DNA sequence(s) into proteins...", dna_sequences.len());
    let proteinas = annotation::translate_batch(&dna_sequences);
    Ok(proteinas)
}

// ---------------------------------------------------------
// MÓDULO 7: SCAFFOLDING
// ---------------------------------------------------------
/// Joins contigs whose ends overlap exactly by `overlap` bases. Returns
/// `(scaffolds, circular)`: the joined sequences and how many of them closed
/// into a circle.
#[pyfunction]
fn build_scaffolds(sequences: Vec<String>, overlap: usize) -> PyResult<(Vec<String>, usize)> {
    let (scaffolds, circulares) = scaffolding::fast_hash_scaffold(sequences, overlap);
    Ok((scaffolds, circulares))
}

// ---------------------------------------------------------
// DESCODIFICADORES VISUAIS E TESTES
// ---------------------------------------------------------
/// Decodes 128-bit k-mers (16 bytes each) of size `k` back into DNA strings.
#[pyfunction]
fn decode_kmers_128(kmer_bytes: &[u8], k: u8) -> PyResult<Vec<String>> {
    let kmers = como_u128(kmer_bytes);
    let seqs: Vec<String> = kmers.par_iter().map(|&bits| {
        let kmer = kmer::Kmer128 { bits, k };
        kmer.decode()
    }).collect();
    Ok(seqs)
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
    let mut kmer_a = ler_kmer(path_a)?;
    let mut kmer_b = ler_kmer(path_b)?;
    multimerge_core::multi_merge_sort_stable(&mut kmer_a);
    multimerge_core::multi_merge_sort_stable(&mut kmer_b);
    let resultado = ordenar_e_mesclar_backend(&kmer_a, &kmer_b, Dispositivo::Auto)?;
    let byte_slice: &[u8] = unsafe { std::slice::from_raw_parts(resultado.as_ptr() as *const u8, resultado.len() * 8) };
    Ok(PyBytes::new(py, byte_slice))
}

/// Legacy 64-bit path. Computes a MinHash signature of `num_hashes` values for
/// each 64-bit k-mer; returns them as raw bytes (`num_hashes` x 8 per k-mer).
#[pyfunction]
fn minhash_signatures<'py>(py: Python<'py>, kmer_bytes: &[u8], num_hashes: usize) -> PyResult<Bound<'py, PyBytes>> {
    let kmers = como_u64(kmer_bytes);
    let assinaturas: Vec<u64> = kmers.par_iter().flat_map(|&k| kmer::compute_minhash_signature(k, num_hashes)).collect();
    let byte_slice: &[u8] = unsafe { std::slice::from_raw_parts(assinaturas.as_ptr() as *const u8, assinaturas.len() * 8) };
    Ok(PyBytes::new(py, byte_slice))
}

/// Legacy 64-bit path. Decodes 64-bit k-mers (8 bytes each) of `length` bases
/// back into DNA strings.
#[pyfunction]
fn decode_kmers(kmer_bytes: &[u8], length: usize) -> PyResult<Vec<String>> {
    let kmers = como_u64(kmer_bytes);
    let seqs: Vec<String> = kmers.par_iter().map(|&k| kmer::decode_kmer(k, length as u8)).collect();
    Ok(seqs)
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

/// MultiMerge Genomics: tools for genomic files of any size. Validate,
/// normalize chromosome names, compare reference genomes, sort, index and
/// query by region; interval operations on BED; statistics, filtering and
/// normalization of VCF; quality report and trimming of FASTQ; k-mer
/// counting; plain text, BGZF, BAM and Parquet in and out; and k-mer modules
/// with optional CUDA acceleration. Start with `engine_info()`; every
/// function has a docstring.
#[pymodule]
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
    Ok(())
}
