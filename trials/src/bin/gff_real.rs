#![allow(dead_code)]
// Ensaio com GFF3 real: ordena e confere em fluxo (sequencia em ordem natural,
// inicio crescente, e no empate o mais longo primeiro).
#[path = "../../../src/multimerge_core.rs"]
pub mod multimerge_core;
#[path = "../../../src/externo.rs"]
mod externo;
use externo::*;
use std::io::BufRead;
use std::time::Instant;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mem_mb: usize = a[1].parse().unwrap();
    let op = Opcoes { memoria_bytes: mem_mb << 20, dir_temp: None, ordem: OrdemCromossomos::Natural };
    let t = Instant::now();
    let r = ordenar_arquivos(&a[3..], &a[2], Formato::Gff, Criterio::Posicao, &op).unwrap();
    println!("sorted in {:.1} s: {} records, {} runs, {} skipped", t.elapsed().as_secs_f64(), r.registros, r.lotes, r.ignoradas);

    let (mut n, mut cab, mut tres) = (0usize, 0usize, 0usize);
    let mut ant: Option<(String, u64, u64)> = None;
    let mut l = String::new();
    let mut leitor = abre_leitura(&a[2]).unwrap();
    loop {
        l.clear();
        if leitor.read_line(&mut l).unwrap() == 0 { break; }
        let t = l.trim_end();
        if t.starts_with("#") { cab += 1; if t == "###" { tres += 1; } assert_eq!(n, 0, "header after a record"); continue; }
        let c: Vec<&str> = t.split("\t").collect();
        let k = (c[0].to_string(), c[3].parse::<u64>().unwrap(), c[4].parse::<u64>().unwrap());
        if let Some(p) = &ant {
            let o = natural_cmp(&p.0, &k.0).then(p.1.cmp(&k.1)).then(k.2.cmp(&p.2));
            assert!(o.is_le(), "OUT OF ORDER at record {n}: {p:?} then {k:?}");
        }
        ant = Some(k);
        n += 1;
    }
    assert_eq!(n, r.registros);
    println!("checked: {n} records in order, {cab} header lines on top, {tres} ### separators in the output");
}
