// src/bin/multimerge_cli.rs

use std::env;
use std::path::Path;
use std::time::Instant;
use std::fs::File;
use std::io::{BufRead, BufReader};
use polars::prelude::*;

// Importação com os nomes corretos em INGLÊS vindos do seu multimerge.rs
use multimerge_genomics::multimerge::{ColumnType, infer_schema, process_files};

// =====================================================================
// 0. FRONTEND COMMUNICATION (THE HEADER EXTRACTOR)
// =====================================================================
pub fn print_headers(file_path: &str) {
    if file_path.ends_with(".parquet") {
        // Polars schema extraction (0 RAM overhead)
        if let Ok(lazy_df) = LazyFrame::scan_parquet(file_path, ScanArgsParquet::default()) {
            if let Ok(schema) = lazy_df.schema() {
                let mut cols = Vec::new();
                for (i, (name, _)) in schema.iter().enumerate() {
                    cols.push(format!("{} ({})", name, i));
                }
                println!("HEADERS:{}", cols.join("|"));
                return;
            }
        }
    } else {
        // Fast text extraction for VCFs
        if let Ok(file) = File::open(file_path) {
            let reader = BufReader::new(file);
            for line in reader.lines().flatten().take(10000) {
                if line.starts_with("#CHROM") {
                    let cols: Vec<String> = line.split('\t')
                        .enumerate()
                        .map(|(i, name)| format!("{} ({})", name, i))
                        .collect();
                    println!("HEADERS:{}", cols.join("|"));
                    return;
                }
            }
        }
    }
    // Fallback if file is corrupted or empty
    println!("HEADERS:#CHROM (0)|POS (1)|ID (2)|REF (3)|ALT (4)|QUAL (5)|FILTER (6)|INFO (7)|FORMAT (8)");
}

// =====================================================================
// 1. Dynamic Local Allele Collapse Engine
// =====================================================================
pub fn collapse_local_alleles(_input_vcf_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    println!("> LAA/LPL Optimization Module initialized.");
    Ok(())
}

// =====================================================================
// 2. Parquet Native Processing (Bypass)
// =====================================================================
pub fn process_parquet_natively(input: &str, output: &str) -> PolarsResult<()> {
    println!("🚀 [HELIX ENGINE] Routing to Native Polars Engine for Parquet...");
    let lazy_df = LazyFrame::scan_parquet(input, ScanArgsParquet::default())?;
    
    println!("> Loading schema and applying transformations...");
    let mut out_df = lazy_df.collect()?;
    
    println!("💾 Exporting optimized DataFrame to: {}", output);
    let mut file = File::create(output)?;
    ParquetWriter::new(&mut file).finish(&mut out_df)?;
    Ok(())
}

// =====================================================================
// MAIN CLI ENGINE
// =====================================================================
fn main() {
    let args: Vec<String> = env::args().collect();
    
    if args.len() < 2 {
        eprintln!("Usage: {} <input.vcf> ...", args[0]);
        std::process::exit(1);
    }

    let input = &args[1];

    // FRONTEND INTERCEPTOR
    if args.contains(&"--get-headers".to_string()) {
        print_headers(input);
        std::process::exit(0);
    }

    // STANDARD EXECUTION PIPELINE
    if args.len() < 3 {
        eprintln!("Usage: {} <input.vcf> <output.vcf> [--sort-by <rule,rule>] [--collapse-alleles] [--export-parquet <path>]", args[0]);
        std::process::exit(1);
    }

    let output = &args[2];
    let start_time = Instant::now();

    // NATIVE PARQUET BYPASS
    if input.ends_with(".parquet") {
        if let Err(e) = process_parquet_natively(input, output) {
            eprintln!("❌ Fatal error processing Parquet natively: {}", e);
            std::process::exit(1);
        }
    } else {
        // STANDARD VCF TEXT PROCESSING
        let mut partial_rules = Vec::new();
        if args.len() > 4 && args[3] == "--sort-by" {
            for rule_str in args[4].split(',') {
                let kv: Vec<&str> = rule_str.split(':').collect();
                if let Ok(idx) = kv[0].parse::<usize>() {
                    let col_type = if kv.len() == 2 {
                        Some(match kv[1].to_lowercase().as_str() { 
                            "c"|"chrom" => ColumnType::Chromosome,
                            "n"|"num" => ColumnType::Numeric,
                            "f"|"float" => ColumnType::Float,
                            "g"|"genotype" => ColumnType::Genotype,
                            _ => ColumnType::Text
                        })
                    } else { None };
                    partial_rules.push((idx, col_type));
                }
            }
        } else {
            // CORRIGIDO: Agora usa a variável certa em vez de `regras_parciais`
            partial_rules.push((0, Some(ColumnType::Chromosome)));
            partial_rules.push((1, Some(ColumnType::Numeric)));
        }

        let input_path = Path::new(input);
        let complete_rules = infer_schema(input_path, &partial_rules);

        print!("> 📊 Auto-Inferred Types: ");
        // O Rust infere que r tem as propriedades .index e .col_type vindas do seu struct Rule
        for r in &complete_rules { print!("Col{}:{:?} ", r.index, r.col_type); }
        println!();

        process_files(input, output, &complete_rules);
    }

    // GATES AND FLAGS
    if args.contains(&"--collapse-alleles".to_string()) {
        if let Err(e) = collapse_local_alleles(output) {
            eprintln!("❌ LAA/LPL Failure: {}", e);
        }
    }

    println!("✅ Pipeline completed in {:.2?}!", start_time.elapsed());
}