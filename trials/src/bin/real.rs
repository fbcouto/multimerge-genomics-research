#![allow(dead_code)]
// Ensaio com VCF real, em dois passos para caber em sessoes curtas:
//
//   cargo run --release --bin real -- sort  <mem_mb> <saida> <entrada>...   ordena com orcamento de memoria
//   cargo run --release --bin real -- check <saida> [entrada]...            confere a saida em fluxo: a
//                                                 chave nao decresce; com as
//                                                 entradas, reconta as variantes
//
// Nada e carregado na memoria em nenhum dos dois.
#[path = "../../../src/multimerge_core.rs"]
pub mod multimerge_core;
#[path = "../../../src/externo.rs"]
mod externo;

use externo::*;
use std::io::BufRead;
use std::time::Instant;

fn variantes(caminho: &str, mut cada: impl FnMut(&str, u64)) -> usize {
    let mut leitor = abre_leitura(caminho).unwrap();
    let mut linha = String::new();
    let mut n = 0usize;
    loop {
        linha.clear();
        if leitor.read_line(&mut linha).unwrap() == 0 {
            return n;
        }
        let t = linha.trim_end();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let mut c = t.split('\t');
        let (Some(crom), Some(pos)) = (c.next(), c.next().and_then(|p| p.parse::<u64>().ok())) else { continue };
        if c.nth(3).is_none() {
            continue;
        }
        cada(crom, pos);
        n += 1;
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let t = Instant::now();
    match a[1].as_str() {
        "sort" => {
            let mem_mb: usize = a[2].parse().unwrap();
            let op = Opcoes { memoria_bytes: mem_mb << 20, dir_temp: None, ordem: OrdemCromossomos::Lexica };
            let r = ordenar_arquivos(&a[4..], &a[3], Formato::Vcf, Criterio::Posicao, &op).unwrap();
            let bytes = std::fs::metadata(&a[3]).unwrap().len();
            println!(
                "sorted in {:.1} s: {} variants, {} runs on disk, {} lines skipped, output of {:.2} GB (budget {mem_mb} MB)",
                t.elapsed().as_secs_f64(),
                r.registros,
                r.lotes,
                r.ignoradas,
                bytes as f64 / 1e9
            );
        }
        "check" => {
            let mut ant: Option<(String, u64)> = None;
            let n = variantes(&a[2], |crom, pos| {
                if let Some((c, p)) = &ant {
                    assert!((c.as_str(), *p) <= (crom, pos), "OUT OF ORDER: ({c}, {p}) > ({crom}, {pos})");
                }
                match &mut ant {
                    Some((c, p)) if c == crom => *p = pos,
                    _ => ant = Some((crom.to_string(), pos)),
                }
            });
            print!("checked in {:.1} s: {n} variants in order", t.elapsed().as_secs_f64());
            if a.len() > 3 {
                let esperado: usize = a[3..].iter().map(|e| variantes(e, |_, _| {})).sum();
                assert_eq!(n, esperado, "the number of variants changed");
                print!("; the recount of the inputs matches ({:.1} s in total)", t.elapsed().as_secs_f64());
            }
            println!();
        }
        _ => eprintln!("usage: real sort <mem_mb> <output> <input>... | real check <output> [input]..."),
    }
}
