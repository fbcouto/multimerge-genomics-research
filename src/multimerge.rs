use rayon::prelude::*;

// --- Generic Utilities ---

fn calculate_minrun(mut n: usize) -> usize {
    let mut r = 0;
    while n >= 64 { 
        r |= n & 1; 
        n >>= 1; 
    }
    n + r
}

pub fn insertion_sort<T: Ord>(arr: &mut [T]) {
    for i in 1..arr.len() {
        let mut j = i;
        while j > 0 && arr[j - 1] > arr[j] {
            arr.swap(j - 1, j);
            j -= 1;
        }
    }
}

fn detect_global_trend<T: Ord + Sync>(arr: &mut [T]) -> bool {
    let n = arr.len();
    if n < 100_000 { return false; } 

    let chunk_size = 32768; 
    let num_blocks = (n + chunk_size - 1) / chunk_size;
    let immutable_arr: &[T] = arr;

    let (ascending, descending) = (0..num_blocks)
        .into_par_iter()
        .map(|i| {
            let start = i * chunk_size;
            let end = std::cmp::min(start + chunk_size + 1, n);
            let overlapping_chunk = &immutable_arr[start..end];
            
            let mut asc = 0;
            let mut desc = 0;
            for j in 1..overlapping_chunk.len() {
                if overlapping_chunk[j - 1] < overlapping_chunk[j] { asc += 1; }
                else if overlapping_chunk[j - 1] > overlapping_chunk[j] { desc += 1; }
            }
            (asc, desc)
        })
        .reduce(|| (0, 0), |a, b| (a.0 + b.0, a.1 + b.1));

    if descending == 0 { return true; }
    if ascending == 0 { arr.reverse(); return true; }
    if descending > 0 && descending < (n / 20) { return false; }
    
    false
}

// --- Adaptive Main Engine ---

pub fn multi_merge_sort_generic<T: Ord + Sync + Send + Clone>(arr: &mut [T]) {
    let n = arr.len();
    
    // Safety check: Prevent panic if the array is empty
    if n == 0 { return; } 

    if n < 1024 {
        insertion_sort(arr);
        return;
    }
    
    if detect_global_trend(arr) { return; }

    let mut is_pure_chaos = false;
    if n > 120 {
        let mid = n / 2;
        let mut direction_changes = 0;
        let mut is_ascending = arr[mid] <= arr[mid + 1];
        
        for i in (mid + 1)..(mid + 100).min(n - 1) {
            let current_direction = arr[i] <= arr[i + 1];
            if current_direction != is_ascending {
                direction_changes += 1;
                is_ascending = current_direction;
            }
        }
        if direction_changes > 15 { is_pure_chaos = true; }
    }

    if !is_pure_chaos {
        let mut buffer = vec![arr[0].clone(); n];
        let num_threads = rayon::current_num_threads();
        let threshold = (n / num_threads).max(1_000_000); 
        parallel_recursive_sort(arr, &mut buffer, threshold);
    } else {
        arr.par_sort_unstable();
    }
}

// --- Processing Core (Path A) ---

fn sequential_timsort_style<T: Ord + Clone>(arr: &mut [T], buffer: &mut [T]) {
    let n = arr.len();
    let minrun = calculate_minrun(n);

    for i in (0..n).step_by(minrun) {
        let end = (i + minrun).min(n);
        insertion_sort(&mut arr[i..end]);
    }

    let mut block_size = minrun;
    while block_size < n {
        for left in (0..n).step_by(block_size * 2) {
            let mid = (left + block_size).min(n);
            let right = (left + block_size * 2).min(n);
            if mid < right {
                stable_merge(&mut arr[left..right], buffer, mid - left);
            }
        }
        block_size *= 2;
    }
}

fn parallel_recursive_sort<T: Ord + Clone + Send>(arr: &mut [T], buffer: &mut [T], threshold: usize) {
    let n = arr.len();
    if n <= threshold {
        sequential_timsort_style(arr, buffer);
        return;
    }

    let mid = n / 2;
    let (arr_left, arr_right) = arr.split_at_mut(mid);
    let (buf_left, buf_right) = buffer.split_at_mut(mid);

    rayon::join(
        || parallel_recursive_sort(arr_left, buf_left, threshold),
        || parallel_recursive_sort(arr_right, buf_right, threshold),
    );

    stable_merge(arr, buffer, mid);
}

fn stable_merge<T: Ord + Clone>(arr: &mut [T], buffer: &mut [T], mid: usize) {
    let n = arr.len();
    buffer[..n].clone_from_slice(&arr[..n]);

    let mut i = 0;
    let mut j = mid;
    let mut k = 0;

    while i < mid && j < n {
        if buffer[i] <= buffer[j] { 
            arr[k] = buffer[i].clone(); 
            i += 1; 
        } else { 
            arr[k] = buffer[j].clone(); 
            j += 1; 
        }
        k += 1;
    }

    if i < mid { 
        arr[k..k + (mid - i)].clone_from_slice(&buffer[i..mid]); 
    } else if j < n { 
        arr[k..k + (n - j)].clone_from_slice(&buffer[j..n]); 
    }
}

use std::path::Path;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

// =====================================================================
// ESTRUTURAS DE DADOS PARA O CLI
// =====================================================================

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ColumnType {
    Chromosome,
    Numeric,
    Float,
    Genotype,
    Text,
}

#[derive(Debug, Clone)]
pub struct Rule {
    pub index: usize,
    pub col_type: ColumnType,
}

// Estrutura que o seu `multi_merge_sort_generic` vai ordenar
#[derive(Clone, Eq, PartialEq)]
pub struct VcfRecord {
    pub raw_line: String,
    pub chrom: String,
    pub num: i64,
    pub float_val: String, // Mantido como string para simplificar ordens parciais, ou parse se preferir
}

// Ensinamos ao Rust como comparar dois VcfRecords para o seu motor os poder ordenar
impl PartialOrd for VcfRecord {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for VcfRecord {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Ordena primeiro por Cromossoma, depois por Número (POS)
        self.chrom.cmp(&other.chrom)
            .then_with(|| self.num.cmp(&other.num))
            .then_with(|| self.float_val.cmp(&other.float_val))
    }
}

// =====================================================================
// FUNÇÕES PÚBLICAS EXPORTADAS PARA O CLI
// =====================================================================

pub fn infer_schema(_path: &Path, partials: &[(usize, Option<ColumnType>)]) -> Vec<Rule> {
    // Mapeia as regras recebidas do terminal para o formato interno
    partials.iter().map(|(idx, col_type)| {
        Rule {
            index: *idx,
            col_type: col_type.clone().unwrap_or(ColumnType::Text),
        }
    }).collect()
}

pub fn process_files(input: &str, output: &str, rules: &[Rule]) {
    println!("🚀 [RUST] Reading and parsing VCF: {}", input);
    
    let file = File::open(input).expect("Error opening the input file");
    let reader = BufReader::with_capacity(1024 * 1024 * 16, file);

    let mut headers = Vec::new();
    let mut records = Vec::new();

    // 1. LEITURA ULTRA-RÁPIDA
    for line_res in reader.lines() {
        if let Ok(line) = line_res {
            if line.starts_with('#') {
                headers.push(line);
                continue;
            }

            let cols: Vec<&str> = line.split('\t').collect();
            let mut chrom = String::new();
            let mut num = 0;
            let mut float_val = String::new();

            // Aplica as regras de ordenação pedidas pelo CLI
            for rule in rules {
                if rule.index < cols.len() {
                    let val = cols[rule.index];
                    match rule.col_type {
                        ColumnType::Chromosome => chrom = val.to_string(),
                        ColumnType::Numeric => num = val.parse::<i64>().unwrap_or(0),
                        ColumnType::Float => float_val = val.to_string(),
                        _ => {}
                    }
                }
            }
            records.push(VcfRecord { raw_line: line, chrom, num, float_val });
        }
    }

    println!("🔄 [HYBRID ENGINE] Sorting {:#?} variants with multi_merge_sort...", records.len());
    
    // 2. Multimerge_Stable, ordenando REFERENCIAS: move 8 bytes por variante
    //    em vez de tres `String`, e variantes com a mesma chave mantem a ordem
    //    do arquivo.
    let mut ordem: Vec<&VcfRecord> = records.iter().collect();
    crate::multimerge_core::multi_merge_sort_stable(&mut ordem);

    println!("💾 [RUST] Writing results to: {}", output);
    let out_file = File::create(output).expect("Error creating the output file");
    let mut writer = BufWriter::with_capacity(1024 * 1024 * 32, out_file);

    // 3. ESCRITA NO DISCO
    for h in headers {
        writeln!(writer, "{}", h).unwrap();
    }
    for r in ordem {
        writeln!(writer, "{}", r.raw_line).unwrap();
    }
}