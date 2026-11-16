#![allow(dead_code)]
// Tests of src/quantificacao.rs: the pseudoalignment index, pseudoalignment,
// the EM and the quantification of FASTQ files (kallisto's method).
//
// Every result is compared with something this file computes on its own:
//  - the index: a naive map (k-mer TEXT, canonical by string comparison with
//    its reverse complement -> the transcripts that contain it), built from
//    the transcripts with no code of the module;
//  - pseudoalignment: the same naive map, intersecting sets read by read;
//    the reads come from a simulator that records the source transcript and
//    the fragment length of every read;
//  - the EM: closed forms and the optimality (KKT) conditions of the
//    likelihood, which is concave; a grid search for a 3-transcript case;
//  - effective lengths: direct sums over the fragment-length distribution;
//  - quantification: the truth of the simulation (fragments drawn from known
//    abundances), and the class counts of the naive pseudoalignment;
//  - the chain with differential expression: planted fold changes;
//  - kallisto's single-end fragment rule: reconstructed here by text search.
// The modes below run the real-data checks: transcripts spliced from a genome
// and its RefSeq GFF3, simulated reads with their truth, and the comparison
// with kallisto 0.46.2's abundance.tsv and pseudoalignment classes.
//
// Usage: cargo run --release --bin quantificacao                      (the tests)
//   -- refseq <genome.fa[.gz]> <annotation.gff[.gz]> <out_dir>         transcripts of a genome + RefSeq GFF3
//   -- synth <out_dir> <genes> <seed>                                    a synthetic transcriptome
//   -- sim <transcripts.fa> <out_dir> <fragments> <seed> [paired=1] [error=0.005] [read_len=100] [frag_mean=200] [frag_sd=20] [zero_fraction=0.2]
//   -- index <transcripts.fa> <out.idx> [k=31] [gene_map|-] [memory_mb]
//   -- quant <index> <out_dir> <name> <r1> [r2|-] [frag_mean frag_sd [filter]]   (filter: kallisto's single-end rule)
//   -- accuracy <index> <truth.tsv> <abundance.tsv>
//   -- compare <index> <ours abundance.tsv> <kallisto abundance.tsv>
//   -- compare-ec <ours ec_counts.tsv> <kallisto pseudoalignments.ec> <kallisto pseudoalignments.tsv>
//   -- dechain <transcripts.fa> <index> <out_dir> <fragments per sample> <seed> [planted genes per direction=15]
//   -- explain <index> <read> [read2]
//   -- demo <folder written by quant_demo.py>
#[path = "../../../src/multimerge_core.rs"]
pub mod multimerge_core;
#[path = "../../../src/bgzf.rs"]
#[allow(dead_code)]
mod bgzf;
#[path = "../../../src/bam.rs"]
#[allow(dead_code)]
mod bam;
#[path = "../../../src/externo.rs"]
#[allow(dead_code)]
mod externo;
#[path = "../../../src/nomes.rs"]
#[allow(dead_code)]
mod nomes;
#[path = "../../../src/expressao.rs"]
#[allow(dead_code)]
mod expressao;
#[path = "../../../src/quantificacao.rs"]
mod quantificacao;

use quantificacao::*;
use std::collections::{BTreeSet, HashMap};
use std::fmt::Debug;
use std::io::{BufRead, Write};
use std::path::PathBuf;
use std::time::Instant;

// ---------------------------------------------------------------- harness

struct T {
    ok: usize,
    falhas: usize,
    dir: PathBuf,
}

impl T {
    fn p(&self, n: &str) -> String {
        self.dir.join(n).to_string_lossy().to_string()
    }
    fn check(&mut self, cond: bool, what: &str) -> bool {
        if cond {
            self.ok += 1;
        } else {
            self.falhas += 1;
            println!("FAIL  {what}");
        }
        cond
    }
    fn eq<A: PartialEq + Debug>(&mut self, got: &A, want: &A, what: &str) -> bool {
        let ok = got == want;
        if !ok {
            println!("FAIL  {what}\n   got:  {got:?}\n   want: {want:?}");
        }
        self.check(ok, what)
    }
    fn perto(&mut self, got: f64, want: f64, tol: f64, what: &str) -> bool {
        let ok = (got - want).abs() <= tol * want.abs().max(1e-300);
        if !ok {
            println!("FAIL  {what}\n   got:  {got:e}\n   want: {want:e} (relative tolerance {tol:e})");
        }
        self.check(ok, what)
    }
    fn maior(&mut self, got: f64, limite: f64, what: &str) -> bool {
        let ok = got > limite;
        if !ok {
            println!("FAIL  {what}: {got} is not above {limite}");
        }
        self.check(ok, what)
    }
    fn menor(&mut self, got: f64, limite: f64, what: &str) -> bool {
        let ok = got < limite;
        if !ok {
            println!("FAIL  {what}: {got} is not below {limite}");
        }
        self.check(ok, what)
    }
    fn erro<X>(&mut self, r: Result<X, String>, agulhas: &[&str], what: &str) {
        match r {
            Ok(_) => {
                println!("FAIL  {what}: should have failed");
                self.check(false, what);
            }
            Err(e) => {
                let falta: Vec<&&str> = agulhas.iter().filter(|a| !e.contains(**a)).collect();
                if !falta.is_empty() {
                    println!("FAIL  {what}: message lacks {falta:?}\n   message: {e}");
                }
                self.check(falta.is_empty(), what);
            }
        }
    }
    fn ok_ou<X>(&mut self, r: Result<X, String>, what: &str) -> Option<X> {
        match r {
            Ok(x) => Some(x),
            Err(e) => {
                println!("FAIL  {what}: {e}");
                self.check(false, what);
                None
            }
        }
    }
}

fn grava(p: &str, texto: &str) {
    std::fs::write(p, texto).unwrap();
}

fn le_tsv(p: &str) -> Vec<Vec<String>> {
    let r = externo::abre_leitura(p).unwrap();
    r.lines().map(|l| l.unwrap().split('\t').map(str::to_string).collect()).collect()
}

// ---------------------------------------------------------------- memory (peak working set)

#[cfg(windows)]
fn pico_mb() -> f64 {
    #[repr(C)]
    struct Pmc {
        cb: u32,
        page_fault_count: u32,
        peak_working_set_size: usize,
        working_set_size: usize,
        quota_peak_paged_pool_usage: usize,
        quota_paged_pool_usage: usize,
        quota_peak_non_paged_pool_usage: usize,
        quota_non_paged_pool_usage: usize,
        pagefile_usage: usize,
        peak_pagefile_usage: usize,
    }
    #[link(name = "kernel32")]
    extern "system" {
        fn GetCurrentProcess() -> isize;
        fn K32GetProcessMemoryInfo(p: isize, c: *mut Pmc, cb: u32) -> i32;
    }
    // SAFETY: a plain struct of integers, filled by the documented Windows call.
    unsafe {
        let mut c: Pmc = std::mem::zeroed();
        c.cb = std::mem::size_of::<Pmc>() as u32;
        if K32GetProcessMemoryInfo(GetCurrentProcess(), &mut c, c.cb) != 0 {
            c.peak_working_set_size as f64 / 1048576.0
        } else {
            0.0
        }
    }
}

#[cfg(not(windows))]
fn pico_mb() -> f64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| s.lines().find(|l| l.starts_with("VmHWM:")).and_then(|l| l.split_whitespace().nth(1).and_then(|x| x.parse::<f64>().ok())))
        .map(|kb| kb / 1024.0)
        .unwrap_or(0.0)
}

// ---------------------------------------------------------------- random numbers

struct Rng(u64);

impl Rng {
    fn u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn unif(&mut self) -> f64 {
        (self.u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn abaixo(&mut self, n: u64) -> u64 {
        ((self.u64() as u128 * n as u128) >> 64) as u64
    }
    fn normal(&mut self) -> f64 {
        let u1 = self.unif().max(1e-300);
        let u2 = self.unif();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[(self.u64() >> 62) as usize]
    }
    fn outra_base(&mut self, b: u8) -> u8 {
        loop {
            let c = self.base();
            if c != b {
                return c;
            }
        }
    }
}

// ---------------------------------------------------------------- naive DNA (no code of the module)

fn rc(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

fn acgt(s: &[u8]) -> bool {
    s.iter().all(|c| matches!(c, b'A' | b'C' | b'G' | b'T'))
}

/// The canonical k-mer by string comparison (A < C < G < T, as the codes).
fn canon(s: &[u8]) -> Vec<u8> {
    let r = rc(s);
    if r.as_slice() < s {
        r
    } else {
        s.to_vec()
    }
}

fn codigo(s: &[u8]) -> u64 {
    s.iter().fold(0u64, |x, c| {
        (x << 2)
            | match c {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                _ => 3,
            }
    })
}

// ---------------------------------------------------------------- synthetic transcriptome

#[derive(Clone)]
struct Tx {
    name: String,
    gene: String,
    symbol: String,
    seq: Vec<u8>,
}

/// Genes of 1 to 7 exons with up to 4 isoforms (exon skipping, alternative
/// last exon); every 10th gene a paralog of the previous one (2%
/// substitutions); a transcript with Ns, one shorter than k, two with the
/// same sequence.
fn transcritoma_sintetico(rng: &mut Rng, genes: usize) -> Vec<Tx> {
    let mut txs = Vec::new();
    let mut exons_de: Vec<Vec<Vec<u8>>> = Vec::new();
    for g in 0..genes {
        let exons: Vec<Vec<u8>> = if g % 10 == 9 {
            exons_de[g - 1].iter().map(|e| e.iter().map(|&b| if rng.unif() < 0.02 { rng.outra_base(b) } else { b }).collect()).collect()
        } else {
            let ne = 1 + rng.abaixo(7) as usize;
            (0..ne).map(|_| (0..(60 + rng.abaixo(340) as usize)).map(|_| rng.base()).collect()).collect()
        };
        let ne = exons.len();
        let n_iso = if ne >= 3 { 1 + rng.abaixo(4) as usize } else { 1 };
        let mut vistas: BTreeSet<Vec<usize>> = BTreeSet::new();
        let mut isos = vec![(0..ne).collect::<Vec<_>>()];
        vistas.insert(isos[0].clone());
        for _ in 1..n_iso {
            for _ in 0..20 {
                let mut v: Vec<usize> = (0..ne).filter(|&i| i == 0 || i == ne - 1 || rng.unif() > 0.4).collect();
                if rng.unif() < 0.3 && v.len() > 2 {
                    v.pop();
                }
                if vistas.insert(v.clone()) {
                    isos.push(v);
                    break;
                }
            }
        }
        for (i, iso) in isos.iter().enumerate() {
            let seq: Vec<u8> = iso.iter().flat_map(|&e| exons[e].iter().copied()).collect();
            txs.push(Tx { name: format!("T{g:04}_{}", i + 1), gene: format!("G{g:04}"), symbol: format!("Sym{g}"), seq });
        }
        exons_de.push(exons);
    }
    // casos especiais
    let mut com_n = txs[3].clone();
    com_n.name = "TN_1".into();
    com_n.gene = "GN".into();
    com_n.symbol = "SymN".into();
    for i in [10usize, 11, 50, 51, 52] {
        if i < com_n.seq.len() {
            com_n.seq[i] = b'N';
        }
    }
    txs.push(com_n);
    txs.push(Tx { name: "TSHORT_1".into(), gene: "GSHORT".into(), symbol: "SymShort".into(), seq: (0..20).map(|_| rng.base()).collect() });
    let gemea_seq: Vec<u8> = (0..800).map(|_| rng.base()).collect();
    txs.push(Tx { name: "TDUP_1".into(), gene: "GDUP".into(), symbol: "SymDup".into(), seq: gemea_seq.clone() });
    txs.push(Tx { name: "TDUP_2".into(), gene: "GDUP".into(), symbol: "SymDup".into(), seq: gemea_seq });
    let mut poli_a: Vec<u8> = (0..400).map(|_| rng.base()).collect();
    poli_a.extend_from_slice(&[b'A'; 25]);
    txs.push(Tx { name: "TPOLYA_1".into(), gene: "GPOLYA".into(), symbol: "SymPolyA".into(), seq: poli_a });
    txs
}

fn escreve_fasta(p: &str, txs: &[Tx]) {
    let mut w = std::io::BufWriter::new(std::fs::File::create(p).unwrap());
    for t in txs {
        writeln!(w, ">{} cdna synthetic gene:{} gene_symbol:{}", t.name, t.gene, t.symbol).unwrap();
        for c in t.seq.chunks(60) {
            w.write_all(c).unwrap();
            w.write_all(b"\n").unwrap();
        }
    }
}

fn le_fasta(p: &str) -> Vec<Tx> {
    let r = externo::abre_leitura(p).unwrap();
    let mut v: Vec<Tx> = Vec::new();
    for l in r.lines() {
        let l = l.unwrap();
        if let Some(h) = l.strip_prefix('>') {
            let (nome, gene, simb) = parse_header(h);
            v.push(Tx { name: nome, gene, symbol: simb, seq: Vec::new() });
        } else if let Some(t) = v.last_mut() {
            t.seq.extend(l.trim().bytes().map(|c| c.to_ascii_uppercase()));
        }
    }
    v
}

// ---------------------------------------------------------------- naive index

struct Ingenuo {
    k: usize,
    mapa: HashMap<Vec<u8>, Vec<u32>>,
}

impl Ingenuo {
    fn novo(txs: &[Tx], k: usize) -> Ingenuo {
        let mut mapa: HashMap<Vec<u8>, Vec<u32>> = HashMap::new();
        for (t, tx) in txs.iter().enumerate() {
            // a cauda poli-A (10 ou mais A no fim) fica fora, como no kallisto
            let mut fim = tx.seq.len();
            if fim >= 10 && tx.seq[fim - 10..].iter().all(|&c| c == b'A') {
                while fim > 0 && tx.seq[fim - 1] == b'A' {
                    fim -= 1;
                }
            }
            let seq = &tx.seq[..fim];
            if seq.len() < k {
                continue;
            }
            for w in seq.windows(k) {
                if !acgt(w) {
                    continue;
                }
                let e = mapa.entry(canon(w)).or_default();
                if e.last() != Some(&(t as u32)) {
                    e.push(t as u32);
                }
            }
        }
        Ingenuo { k, mapa }
    }

    /// 0 too short, 1 no k-mer in the index, 2 conflicting, 3 assigned.
    fn classe(&self, mates: &[&[u8]]) -> (u8, Vec<u32>) {
        let mut valido = false;
        let mut acc: Option<Vec<u32>> = None;
        for m in mates {
            if m.len() < self.k {
                continue;
            }
            for w in m.windows(self.k) {
                if !acgt(w) {
                    continue;
                }
                valido = true;
                if let Some(v) = self.mapa.get(&canon(w)) {
                    acc = Some(match acc {
                        None => v.clone(),
                        Some(a) => a.into_iter().filter(|x| v.contains(x)).collect(),
                    });
                }
            }
        }
        match acc {
            None => (if valido { 1 } else { 0 }, Vec::new()),
            Some(a) if a.is_empty() => (2, Vec::new()),
            Some(a) => (3, a),
        }
    }
}

fn status_de(c: &Classe) -> u8 {
    match c {
        Classe::Curta => 0,
        Classe::SemKmer => 1,
        Classe::Conflito => 2,
        _ => 3,
    }
}

// ---------------------------------------------------------------- read simulation

struct SimOpt {
    n: usize,
    read_len: usize,
    media: f64,
    dp: f64,
    erro: f64,
    paired: bool,
}

struct Sim {
    r1: Vec<Vec<u8>>,
    r2: Vec<Vec<u8>>,
    origem: Vec<u32>,
    frag: Vec<u32>,
    contagem: Vec<u64>,
}

/// Effective length under the simulator's model, by a direct sum:
/// sum_f p(f) (L - f + 1) / sum_f p(f), f = 1..min(L, 999).
fn eff_direto(l: usize, media: f64, dp: f64) -> f64 {
    let (mut num, mut den) = (0.0, 0.0);
    for f in 0..=l.min(999) {
        let p = (-(f as f64 - media).powi(2) / (2.0 * dp * dp)).exp();
        num += p * (l as f64 - f as f64 + 1.0);
        den += p;
    }
    if den > 0.0 {
        num / den
    } else {
        l as f64 + 1.0
    }
}

/// Fragments drawn with probability proportional to theta x effective length;
/// length ~ N(mean, sd) given 1 <= f <= min(L, 999); uniform start; random
/// strand; substitution errors.
fn simula(seqs: &[&[u8]], theta: &[f64], rng: &mut Rng, o: &SimOpt) -> Sim {
    let pesos: Vec<f64> = seqs.iter().zip(theta).map(|(s, &th)| if s.len() < 31 { 0.0 } else { th * eff_direto(s.len(), o.media, o.dp) }).collect();
    let mut cum = Vec::with_capacity(pesos.len());
    let mut acc = 0.0;
    for p in &pesos {
        acc += p;
        cum.push(acc);
    }
    let mut s = Sim { r1: Vec::with_capacity(o.n), r2: Vec::new(), origem: Vec::with_capacity(o.n), frag: Vec::with_capacity(o.n), contagem: vec![0; seqs.len()] };
    let erra = |r: &mut Vec<u8>, rng: &mut Rng| {
        for b in r.iter_mut() {
            if rng.unif() < o.erro {
                *b = rng.outra_base(*b);
            }
        }
    };
    for _ in 0..o.n {
        let u = rng.unif() * acc;
        let t = cum.partition_point(|&c| c <= u).min(seqs.len() - 1);
        let seq = seqs[t];
        let l = seq.len();
        let mut f = 0usize;
        for _ in 0..100 {
            let x = (o.media + o.dp * rng.normal()).round();
            if x >= 1.0 && x as usize <= l && (x as usize) < 1000 {
                f = x as usize;
                break;
            }
        }
        if f == 0 {
            f = l.min(999);
        }
        let ini = rng.abaixo((l - f + 1) as u64) as usize;
        let mut fr = seq[ini..ini + f].to_vec();
        if rng.unif() < 0.5 {
            fr = rc(&fr);
        }
        let n = o.read_len.min(f);
        let mut a = fr[..n].to_vec();
        erra(&mut a, rng);
        s.r1.push(a);
        if o.paired {
            let mut b = rc(&fr)[..n].to_vec();
            erra(&mut b, rng);
            s.r2.push(b);
        }
        s.origem.push(t as u32);
        s.frag.push(f as u32);
        s.contagem[t] += 1;
    }
    s
}

fn escreve_fastq(p: &str, reads: &[Vec<u8>], sufixo: &str) {
    let f = std::fs::File::create(p).unwrap();
    let mut w: Box<dyn Write> = if p.ends_with(".gz") {
        Box::new(std::io::BufWriter::new(flate2::write::GzEncoder::new(f, flate2::Compression::fast())))
    } else {
        Box::new(std::io::BufWriter::new(f))
    };
    let mut q = Vec::new();
    for (i, r) in reads.iter().enumerate() {
        q.clear();
        q.resize(r.len(), b'I');
        writeln!(w, "@r{i}{sufixo} simulated").unwrap();
        w.write_all(r).unwrap();
        w.write_all(b"\n+\n").unwrap();
        w.write_all(&q).unwrap();
        w.write_all(b"\n").unwrap();
    }
    w.flush().unwrap();
}

// ---------------------------------------------------------------- statistics

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        sab += (x - ma) * (y - mb);
        saa += (x - ma).powi(2);
        sbb += (y - mb).powi(2);
    }
    sab / (saa * sbb).sqrt()
}

fn postos(v: &[f64]) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| v[a].total_cmp(&v[b]));
    let mut r = vec![0.0; v.len()];
    let mut i = 0;
    while i < idx.len() {
        let mut j = i;
        while j + 1 < idx.len() && v[idx[j + 1]] == v[idx[i]] {
            j += 1;
        }
        let m = (i + j) as f64 / 2.0 + 1.0;
        for x in &idx[i..=j] {
            r[*x] = m;
        }
        i = j + 1;
    }
    r
}

fn spearman(a: &[f64], b: &[f64]) -> f64 {
    pearson(&postos(a), &postos(b))
}

/// Mean absolute relative difference: |x - y| / (x + y), 0 when both are 0.
fn mard(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| if x + y > 0.0 { (x - y).abs() / (x + y) } else { 0.0 }).sum::<f64>() / a.len() as f64
}

struct Metricas {
    n: usize,
    pearson_log: f64,
    spearman: f64,
    mard: f64,
}

/// Over the items expressed in either vector: Pearson of log(1 + x),
/// Spearman, MARD.
fn metricas(est: &[f64], verdade: &[f64]) -> Metricas {
    let (mut a, mut b) = (Vec::new(), Vec::new());
    for (x, y) in est.iter().zip(verdade) {
        if *x > 0.0 || *y > 0.0 {
            a.push(*x);
            b.push(*y);
        }
    }
    let la: Vec<f64> = a.iter().map(|x| (1.0 + x).ln()).collect();
    let lb: Vec<f64> = b.iter().map(|x| (1.0 + x).ln()).collect();
    Metricas { n: a.len(), pearson_log: pearson(&la, &lb), spearman: spearman(&a, &b), mard: mard(&a, &b) }
}

fn mostra(nome: &str, m: &Metricas) {
    println!("  {nome:<34} n = {:>6}  Pearson(log1p) = {:.4}  Spearman = {:.4}  MARD = {:.4}", m.n, m.pearson_log, m.spearman, m.mard);
}

fn por_gene(v: &[f64], gene_de: &[usize], ng: usize) -> Vec<f64> {
    let mut g = vec![0.0; ng];
    for (t, x) in v.iter().enumerate() {
        g[gene_de[t]] += x;
    }
    g
}

/// gene index of each transcript, by the gene name in the Tx list.
fn genes_de(txs: &[Tx]) -> (Vec<usize>, Vec<String>) {
    let mut pos: HashMap<String, usize> = HashMap::new();
    let mut nomes = Vec::new();
    let v = txs
        .iter()
        .map(|t| {
            *pos.entry(t.gene.clone()).or_insert_with(|| {
                nomes.push(t.gene.clone());
                nomes.len() - 1
            })
        })
        .collect();
    (v, nomes)
}

/// Molecule abundances: log-normal, a fraction of transcripts silent.
fn abundancias(rng: &mut Rng, n: usize, zeros: f64) -> Vec<f64> {
    (0..n).map(|_| if rng.unif() < zeros { 0.0 } else { (rng.normal() * 1.5 + 3.0).exp() }).collect()
}

// ================================================================ the tests

fn testa_kmers(t: &mut T) {
    let mut rng = Rng(11);
    for &k in &[7usize, 15, 21, 31] {
        let mut ok = true;
        for _ in 0..50 {
            let mut s: Vec<u8> = (0..200).map(|_| rng.base()).collect();
            for _ in 0..3 {
                let i = rng.abaixo(200) as usize;
                s[i] = b'N';
            }
            // lowercase and U must read as the same bases
            let mut s2 = s.clone();
            for c in s2.iter_mut().step_by(7) {
                *c = c.to_ascii_lowercase();
            }
            let mut vistos = Vec::new();
            kmers(&s2, k, |pos, f, r| vistos.push((pos, f, r)));
            let mut esperado = Vec::new();
            for (pos, w) in s.windows(k).enumerate() {
                if acgt(w) {
                    esperado.push((pos, codigo(w), codigo(&rc(w))));
                }
            }
            ok &= vistos == esperado;
        }
        t.check(ok, &format!("k = {k}: rolling k-mers, reverse complements, N skipped, lowercase read as uppercase"));
    }
    let mut v = Vec::new();
    kmers(b"ACGUACGUACGUA", 7, |_, f, _| v.push(f));
    t.eq(&v[0], &codigo(b"ACGTACG"), "U reads as T");
    t.eq(&decodifica_kmer(codigo(b"GATTACAGATTACAGATTACAGATTACAGAT"), 31), &"GATTACAGATTACAGATTACAGATTACAGAT".to_string(), "k-mer decoding");
    // a mistura e uma bijecao: sem colisao em 1 milhao de chaves
    let mut hs: Vec<u64> = (0..1_000_000u64).map(|i| mistura(i * 7919)).collect();
    hs.sort_unstable();
    hs.dedup();
    t.eq(&hs.len(), &1_000_000, "the mixing function has no collision in 1,000,000 keys");
    t.eq(&codifica(b"ACGN"), &None, "codifica refuses N");
}

fn testa_cabecalhos(t: &mut T) {
    let casos: [(&str, (&str, &str, &str)); 7] = [
        (
            "ENST00000631435.1 cdna chromosome:GRCh38:CHR_HSCHR7_2_CTG6:142847306:142847317:1 gene:ENSG00000282253.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene gene_symbol:TRBD1 description:T cell receptor beta diversity 1",
            ("ENST00000631435.1", "ENSG00000282253.1", "TRBD1"),
        ),
        (
            "ENST00000456328.2|ENSG00000290825.1|-|-|DDX11L2-202|DDX11L2|1657|lncRNA|",
            ("ENST00000456328.2", "ENSG00000290825.1", "DDX11L2"),
        ),
        ("NM_000014.6 Homo sapiens alpha-2-macroglobulin (A2M), transcript variant 1, mRNA", ("NM_000014.6", "A2M", "A2M")),
        (
            "NM_006516.4 Homo sapiens solute carrier family 2 (facilitated glucose transporter), member 1 (SLC2A1), mRNA",
            ("NM_006516.4", "SLC2A1", "SLC2A1"),
        ),
        ("XR_001737578.2 PREDICTED: Homo sapiens uncharacterized LOC105379484 (LOC105379484), ncRNA", ("XR_001737578.2", "LOC105379484", "LOC105379484")),
        ("tx1 some free text (ABC)", ("tx1", "", "")),
        ("  tx2  ", ("tx2", "", "")),
    ];
    for (h, (a, b, c)) in casos {
        let r = parse_header(h);
        t.eq(&r, &(a.to_string(), b.to_string(), c.to_string()), &format!("header {h:.40}"));
    }
}

fn testa_mapas(t: &mut T) {
    // um FASTA de 4 transcritos sem gene no cabecalho
    let mut rng = Rng(5);
    let seqs: Vec<Vec<u8>> = (0..4).map(|_| (0..300).map(|_| rng.base()).collect()).collect();
    let nomes = ["ENST01.2", "ENST02.1", "NM_9.3", "rna5"];
    let mut fa = String::new();
    for (n, s) in nomes.iter().zip(&seqs) {
        fa.push_str(&format!(">{n}\n{}\n", String::from_utf8_lossy(s)));
    }
    let pfa = t.p("mapas.fa");
    grava(&pfa, &fa);
    let genes_do = |t: &mut T, mapa: Option<String>, what: &str| -> Option<Vec<(String, String)>> {
        let idx = t.p("mapas.idx");
        let r = t.ok_ou(build_transcript_index(&pfa, &idx, &IndexOptions { gene_map: mapa, ..Default::default() }), what)?;
        let i = TranscriptIndex::load(&r.index).unwrap();
        Some((0..4).map(|x| (i.gene_ids[i.gene_of[x] as usize].clone(), i.gene_names[i.gene_of[x] as usize].clone())).collect())
    };
    // sem mapa: cada transcrito e o proprio gene, com aviso
    let idx = t.p("mapas.idx");
    if let Some(r) = t.ok_ou(build_transcript_index(&pfa, &idx, &IndexOptions::default()), "index without any gene") {
        t.eq(&r.genes, &4, "no gene map: each transcript its own gene");
        t.check(r.warnings.iter().any(|w| w.contains("gene map")), "no gene map: a warning says how to give one");
    }
    // GTF do Ensembl: ids sem versao + transcript_version
    let gtf = "#!genome-build GRCh38\n\
1\thavana\tgene\t1\t1000\t.\t+\t.\tgene_id \"ENSG1\"; gene_version \"5\"; gene_name \"AAA\";\n\
1\thavana\ttranscript\t1\t1000\t.\t+\t.\tgene_id \"ENSG1\"; gene_version \"5\"; transcript_id \"ENST01\"; transcript_version \"2\"; gene_name \"AAA\";\n\
1\thavana\texon\t1\t1000\t.\t+\t.\tgene_id \"ENSG1\"; transcript_id \"ENST01\"; transcript_version \"2\"; gene_name \"AAA\";\n\
1\thavana\ttranscript\t1\t1000\t.\t+\t.\tgene_id \"ENSG1\"; transcript_id \"ENST02\"; transcript_version \"1\"; gene_name \"AAA\";\n\
1\tRefSeq\texon\t1\t1000\t.\t+\t.\tgene_id \"BBB\"; transcript_id \"NM_9.3\"; gene \"BBB\";\n";
    let p = t.p("mapa.gtf");
    grava(&p, gtf);
    if let Some(g) = genes_do(t, Some(p.clone()), "GTF gene map") {
        t.eq(
            &g,
            &vec![("ENSG1".into(), "AAA".into()), ("ENSG1".into(), "AAA".into()), ("BBB".into(), "BBB".into()), ("rna5".into(), "rna5".into())],
            "GTF: versioned FASTA names found through transcript_version; the absent one is its own gene",
        );
    }
    // GFF3 do Ensembl e do RefSeq (miRNA sob primary_transcript)
    let gff = "##gff-version 3\n\
1\tensembl\tgene\t1\t1000\t.\t+\t.\tID=gene:ENSG1;Name=AAA;gene_id=ENSG1;version=5\n\
1\tensembl\tmRNA\t1\t1000\t.\t+\t.\tID=transcript:ENST01;Parent=gene:ENSG1;Name=AAA-201;transcript_id=ENST01;version=2\n\
1\tensembl\texon\t1\t1000\t.\t+\t.\tParent=transcript:ENST01;Name=ENSE1\n\
1\tensembl\tlnc_RNA\t1\t1000\t.\t+\t.\tID=transcript:ENST02;Parent=gene:ENSG1;transcript_id=ENST02;version=1\n\
1\tRefSeq\tgene\t1\t1000\t.\t+\t.\tID=gene-BBB;Name=BBB%2Dx;gene=BBB-x\n\
1\tRefSeq\tmRNA\t1\t1000\t.\t+\t.\tID=rna-NM_9.3;Parent=gene-BBB;transcript_id=NM_9.3\n\
1\tRefSeq\texon\t1\t1000\t.\t+\t.\tID=exon-NM_9.3-1;Parent=rna-NM_9.3\n\
1\tRefSeq\tgene\t1\t1000\t.\t+\t.\tID=gene-MIR1;Name=MIR1\n\
1\tRefSeq\tprimary_transcript\t1\t1000\t.\t+\t.\tID=rna-NR_1.1;Parent=gene-MIR1\n\
1\tRefSeq\tmiRNA\t1\t20\t.\t+\t.\tID=rna-rna5;Parent=rna-NR_1.1\n";
    let p = t.p("mapa.gff3");
    grava(&p, gff);
    if let Some(g) = genes_do(t, Some(p.clone()), "GFF3 gene map") {
        t.eq(
            &g,
            &vec![("ENSG1".into(), "AAA".into()), ("ENSG1".into(), "AAA".into()), ("BBB".into(), "BBB-x".into()), ("MIR1".into(), "MIR1".into())],
            "GFF3: Ensembl transcripts by ID and version, RefSeq by ID, a miRNA through its primary transcript to the gene, %-escapes decoded",
        );
    }
    // tabela com cabecalho, separada por virgula, sem a terceira coluna numa linha
    let p = t.p("mapa.csv");
    grava(&p, "transcript_id,gene_id,gene_name\nENST01,G1,One\nENST02,G1,One\nNM_9,G9\nrna5,G5,Five\n");
    if let Some(g) = genes_do(t, Some(p.clone()), "table gene map") {
        t.eq(
            &g,
            &vec![("G1".into(), "One".into()), ("G1".into(), "One".into()), ("G9".into(), "G9".into()), ("G5".into(), "Five".into())],
            "table: header skipped, unversioned map names match versioned FASTA names, missing name = gene id",
        );
    }
    let p = t.p("mapa_ruim.tsv");
    grava(&p, "x1\tg1\nx2\tg2\n");
    t.erro(build_transcript_index(&pfa, &idx, &IndexOptions { gene_map: Some(p), ..Default::default() }), &["None of the 4", "same release"], "a map of other transcripts");
    let p = t.p("mapa_uma_coluna.tsv");
    grava(&p, "abc\ndef\n");
    t.erro(build_transcript_index(&pfa, &idx, &IndexOptions { gene_map: Some(p), ..Default::default() }), &["two columns"], "a one-column map");
}

/// Builds and checks an index against the naive map.
fn confere_indice(t: &mut T, txs: &[Tx], i: &TranscriptIndex, k: usize, rotulo: &str) {
    let ing = Ingenuo::novo(txs, k);
    t.eq(&i.kmers(), &ing.mapa.len(), &format!("{rotulo}: number of distinct k-mers"));
    let mut ruins = 0;
    for (km, ts) in &ing.mapa {
        match i.lookup(codigo(km)) {
            Some(e) if i.ec(e) == ts.as_slice() => {}
            _ => ruins += 1,
        }
    }
    t.eq(&ruins, &0, &format!("{rotulo}: every k-mer finds exactly its transcripts"));
    let mut rng = Rng(77);
    let mut falsos = 0;
    let mut testados = 0;
    while testados < 100_000 {
        let km: Vec<u8> = (0..k).map(|_| rng.base()).collect();
        let c = canon(&km);
        if ing.mapa.contains_key(&c) {
            continue;
        }
        testados += 1;
        if i.lookup(codigo(&c)).is_some() {
            falsos += 1;
        }
    }
    t.eq(&falsos, &0, &format!("{rotulo}: 100,000 absent k-mers are not found"));
    // as classes: as T primeiras sao os transcritos sozinhos; as outras distintas e ordenadas
    let nt = i.transcripts();
    let mut ok = (0..nt).all(|x| i.ec(x as u32) == [x as u32]);
    let mut vistas = BTreeSet::new();
    for e in nt..i.classes() {
        let v = i.ec(e as u32);
        ok &= v.len() >= 2 && v.windows(2).all(|w| w[0] < w[1]) && vistas.insert(v.to_vec());
    }
    t.check(ok, &format!("{rotulo}: classes are single transcripts first, then distinct sorted sets"));
}

fn testa_indice(t: &mut T, txs: &[Tx], pfa: &str) {
    for &k in &[31usize, 21, 15] {
        let p = t.p(&format!("synt_k{k}.idx"));
        let Some(r) = t.ok_ou(build_transcript_index(pfa, &p, &IndexOptions { k, ..Default::default() }), &format!("index k={k}")) else {
            continue;
        };
        t.eq(&r.transcripts, &txs.len(), &format!("k={k}: transcripts"));
        t.eq(&r.short_transcripts, &(if k > 20 { 1 } else { 0 }), &format!("k={k}: the 20-bp transcript is reported as shorter than k (only when k > 20)"));
        t.eq(&r.non_acgt_bases, &5, &format!("k={k}: the 5 Ns are counted"));
        t.eq(&r.poly_a_clipped, &1, &format!("k={k}: one poly-A tail clipped"));
        let i = TranscriptIndex::load(&p).unwrap();
        confere_indice(t, txs, &i, k, &format!("k={k}"));
        t.eq(&r.kmers, &(i.kmers() as u64), &format!("k={k}: the summary counts the k-mers loaded"));
        t.eq(&r.bytes, &std::fs::metadata(&p).unwrap().len(), &format!("k={k}: the summary gives the file size"));
        // genes do cabecalho
        let g = i.gene_of[0] as usize;
        t.eq(&(i.gene_ids[g].as_str(), i.gene_names[g].as_str()), &("G0000", "Sym0"), &format!("k={k}: gene from the Ensembl-style header"));
        if k == 31 {
            let info = transcript_index_info(&p).unwrap();
            for (j, rot) in [(r.json(), "build"), (info.json(), "info")] {
                match serde_json::from_str::<serde_json::Value>(&j) {
                    Ok(v) => t.check(v["kmers"].as_u64() == Some(r.kmers) && v["text"].is_string(), &format!("index summary JSON ({rot}) parses, with its keys")),
                    Err(e) => t.check(false, &format!("index summary JSON ({rot}): {e}")),
                };
            }
            t.eq(&(info.k, info.transcripts, info.kmers, info.equivalence_classes, info.genes), &(31, r.transcripts, r.kmers, r.equivalence_classes, r.genes), "transcript_index_info reads the header");
            for x in [0usize, 5, txs.len() - 1] {
                t.eq(&i.sequence(x).into_bytes(), &txs[x].seq, &format!("transcript {x} stored as it is (Ns included)"));
            }
        }
    }
    // cauda poli-A: fora dos k-mers, mas no comprimento
    let i31 = TranscriptIndex::load(&t.p("synt_k31.idx")).unwrap();
    let npa = txs.iter().position(|x| x.name == "TPOLYA_1").unwrap();
    let sq = &txs[npa].seq;
    let a31 = sq[sq.len() - 31..].to_vec();
    t.check(i31.lookup(codigo(&canon(&a31))).is_none(), "poly-A: the last k-mer, into the A tail, is not in the index");
    let antes = &sq[sq.len() - 25 - 31..sq.len() - 25];
    t.check(i31.lookup(codigo(&canon(antes))).is_some(), "poly-A: the k-mer just before the tail is");
    t.eq(&i31.lengths[npa], &425, "poly-A: the length stays the original one (425)");
    let psem = t.p("synt_nopolya.idx");
    if let Some(r) = t.ok_ou(build_transcript_index(pfa, &psem, &IndexOptions { clip_poly_a: false, ..Default::default() }), "index without poly-A clipping") {
        let i = TranscriptIndex::load(&psem).unwrap();
        t.eq(&i.lookup(codigo(&canon(&a31))).map(|e| i.ec(e).to_vec()), &Some(vec![npa as u32]), "clip_poly_a=False keeps the tail's k-mers");
        t.eq(&r.poly_a_clipped, &0, "clip_poly_a=False clips nothing");
    }
    // passadas e threads nao mudam um byte
    let base = t.p("synt_k31.idx");
    let p2 = t.p("synt_passes.idx");
    if let Some(r) = t.ok_ou(build_transcript_index(pfa, &p2, &IndexOptions { memory_mb: Some(1), ..Default::default() }), "index in many passes") {
        t.check(r.passes > 3, &format!("a 1 MB budget makes several passes ({})", r.passes));
        t.check(std::fs::read(&base).unwrap() == std::fs::read(&p2).unwrap(), "the index in several passes is the same file, byte for byte");
    }
    let p3 = t.p("synt_1thread.idx");
    let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    if let Some(_r) = t.ok_ou(pool.install(|| build_transcript_index(pfa, &p3, &IndexOptions::default())), "index with 1 thread") {
        t.check(std::fs::read(&base).unwrap() == std::fs::read(&p3).unwrap(), "the index built by 1 thread is the same file, byte for byte");
    }
    t.check(!std::path::Path::new(&format!("{base}.part")).exists() && !std::path::Path::new(&format!("{base}.kmers.tmp")).exists(), "no temporary file is left");
    // fasta comprimido = mesmo indice
    let pgz = t.p("synt.fa.gz");
    {
        let mut w = flate2::write::GzEncoder::new(std::fs::File::create(&pgz).unwrap(), flate2::Compression::default());
        w.write_all(&std::fs::read(pfa).unwrap()).unwrap();
        w.finish().unwrap();
    }
    let p4 = t.p("synt_gz.idx");
    if t.ok_ou(build_transcript_index(&pgz, &p4, &IndexOptions::default()), "index from .fa.gz").is_some() {
        let a = TranscriptIndex::load(&base).unwrap();
        let b = TranscriptIndex::load(&p4).unwrap();
        t.check(a.kmers() == b.kmers() && a.classes() == b.classes() && a.names == b.names, "a gzip FASTA gives the same index");
    }
    // arquivos danificados
    let bytes = std::fs::read(&base).unwrap();
    let pt = t.p("truncado.idx");
    std::fs::write(&pt, &bytes[..bytes.len() / 2]).unwrap();
    t.erro(TranscriptIndex::load(&pt), &["truncated", "build it again"], "a truncated index");
    let mut v = bytes.clone();
    v[8] = 9;
    std::fs::write(&pt, &v).unwrap();
    t.erro(TranscriptIndex::load(&pt), &["format version 9", "reads version 1"], "an index of another format version");
    std::fs::write(&pt, b"not an index at all, just text").unwrap();
    t.erro(TranscriptIndex::load(&pt), &["not an MMGEN transcript index"], "a file that is not an index");
    let mut v = bytes.clone();
    let n = v.len();
    v[n - 20] ^= 0xFF;
    std::fs::write(&pt, &v).unwrap();
    let r = TranscriptIndex::load(&pt);
    t.check(r.is_err(), "a damaged byte near the end is noticed");
    t.erro(TranscriptIndex::load(&t.p("nao_existe.idx")), &["Error opening the transcript index"], "a missing index");
}

fn testa_pseudoalinhamento(t: &mut T, txs: &[Tx], idx: &TranscriptIndex) {
    let ing = Ingenuo::novo(txs, 31);
    let mut rng = Rng(2024);
    let seqs: Vec<&[u8]> = txs.iter().map(|x| x.seq.as_slice()).collect();
    let theta = vec![1.0; txs.len()];
    for &(paired, erro) in &[(false, 0.0), (false, 0.02), (true, 0.0), (true, 0.02)] {
        let sim = simula(&seqs, &theta, &mut rng, &SimOpt { n: 20_000, read_len: 75, media: 180.0, dp: 30.0, erro, paired });
        let mut difs = 0;
        let mut contas = [0usize; 4];
        let mut contem_origem = 0;
        let mut atribuidas = 0;
        for i in 0..sim.r1.len() {
            let mates: Vec<&[u8]> = if paired { vec![&sim.r1[i], &sim.r2[i]] } else { vec![&sim.r1[i]] };
            let (st, ts) = ing.classe(&mates);
            let c = idx.pseudoalign(&sim.r1[i], if paired { Some(&sim.r2[i]) } else { None });
            contas[st as usize] += 1;
            if st != status_de(&c) || ts != idx.transcripts_of(&c) {
                difs += 1;
            }
            if st == 3 {
                atribuidas += 1;
                if ts.contains(&sim.origem[i]) {
                    contem_origem += 1;
                }
            }
        }
        let rot = format!("{}, {}% errors", if paired { "pairs" } else { "single reads" }, erro * 100.0);
        t.eq(&difs, &0, &format!("{rot}: 20,000 pseudoalignments equal to the naive intersection"));
        if erro == 0.0 {
            t.eq(&contem_origem, &atribuidas, &format!("{rot}: without errors the source transcript is always in the class"));
        } else {
            t.maior(contem_origem as f64 / atribuidas as f64, 0.995, &format!("{rot}: with errors the source is in the class"));
        }
        println!("  {rot}: assigned {}, no k-mer {}, conflicting {}, too short {}", contas[3], contas[1], contas[2], contas[0]);
    }
    // casos feitos a mao
    let s = &txs[0].seq;
    let leitura = &s[10..110];
    t.eq(&idx.pseudoalign(&leitura[..20], None), &Classe::Curta, "a read shorter than k is too short");
    let tudo_n = vec![b'N'; 80];
    t.eq(&idx.pseudoalign(&tudo_n, None), &Classe::Curta, "a read of Ns has no k-mer");
    let lixo: Vec<u8> = (0..100).map(|_| rng.base()).collect();
    t.eq(&idx.pseudoalign(&lixo, None), &Classe::SemKmer, "a random read has k-mers, none in the index");
    let c1 = idx.pseudoalign(leitura, None);
    let c2 = idx.pseudoalign(&rc(leitura), None);
    t.eq(&c1, &c2, "a read and its reverse complement pseudoalign the same (canonical k-mers)");
    t.eq(&idx.pseudoalign(leitura, Some(&lixo)), &c1, "a mate with no k-mer in the index adds nothing to the pair");
    // um par cujos mates caem em genes sem nada em comum: conflito
    let outro = &txs[txs.len() - 6].seq;
    t.eq(&idx.pseudoalign(leitura, Some(&outro[0..100])), &Classe::Conflito, "mates on unrelated genes: conflicting");
    // quimera dentro de uma leitura
    let mut quim = s[0..60].to_vec();
    quim.extend_from_slice(&outro[0..60]);
    t.eq(&idx.pseudoalign(&quim, None), &Classe::Conflito, "a chimeric read: conflicting");
    // as duas copias identicas
    let ndup = txs.iter().position(|x| x.name == "TDUP_1").unwrap() as u32;
    let d = idx.pseudoalign(&txs[ndup as usize].seq[100..200], None);
    t.eq(&idx.transcripts_of(&d), &vec![ndup, ndup + 1], "identical transcripts share every class");
}

fn testa_em(t: &mut T) {
    // dois transcritos: {A} c1, {B} c2, {A,B} c3, comprimentos efetivos la, lb.
    // Ponto fixo: x = c1 + c3 (x/la) / (x/la + (N-x)/lb), resolvido por bissecao.
    let rigido = EmOptions { min_rounds: 50, max_rounds: 100_000, change: 1e-12, change_limit: 1e-12, alpha_limit: 1e-30 };
    for &(c1, c2, c3, la, lb) in &[(100.0, 50.0, 300.0, 1000.0, 1000.0), (10.0, 200.0, 1000.0, 500.0, 2000.0), (0.0, 30.0, 70.0, 800.0, 300.0)] {
        let n: f64 = c1 + c2 + c3;
        let g = |x: f64| c1 + c3 * (x / la) / (x / la + (n - x) / lb) - x;
        let (mut lo, mut hi) = (c1, c1 + c3);
        for _ in 0..200 {
            let m = 0.5 * (lo + hi);
            if g(m) > 0.0 {
                lo = m;
            } else {
                hi = m;
            }
        }
        let x = 0.5 * (lo + hi);
        let (a, rodadas) = em(2, &[0, 1, 2, 4], &[0, 1, 0, 1], &[c1, c2, c3], &[la, lb], &rigido);
        let ok = (a[0] - x).abs() <= 1e-6 * n && (a[1] - (n - x)).abs() <= 1e-6 * n;
        if !ok {
            println!("   EM {a:?}, closed form {x} {}", n - x);
        }
        t.check(ok, &format!("EM, 2 transcripts ({c1},{c2},{c3}; {la},{lb}): the fixed point solved by bisection ({rodadas} rounds)"));
    }
    // tres transcritos: condicoes KKT (a verossimilhanca e concava) e uma grade
    let ofs = [0usize, 2, 4, 7, 8, 10];
    let dados = [0u32, 1, 1, 2, 0, 1, 2, 0, 0, 2];
    let cont = [120.0, 80.0, 300.0, 40.0, 25.0];
    let eff = [900.0, 1500.0, 400.0];
    let (a, _) = em(3, &ofs, &dados, &cont, &eff, &rigido);
    let n: f64 = cont.iter().sum();
    t.perto(a.iter().sum::<f64>(), n, 1e-9, "EM keeps every read");
    let eta: Vec<f64> = a.iter().map(|x| x / n).collect();
    let mut grad = [0.0; 3];
    let mut logl = 0.0;
    for c in 0..cont.len() {
        let ts = &dados[ofs[c]..ofs[c + 1]];
        let s: f64 = ts.iter().map(|&u| eta[u as usize] / eff[u as usize]).sum();
        logl += cont[c] * s.ln();
        for &u in ts {
            grad[u as usize] += cont[c] / (eff[u as usize] * s);
        }
    }
    let ok = (0..3).all(|u| if eta[u] > 1e-9 { (grad[u] - n).abs() < 1e-5 * n } else { grad[u] <= n * (1.0 + 1e-6) });
    t.check(ok, "EM, 3 transcripts: the optimality conditions of the likelihood hold");
    let ll = |e: [f64; 3]| -> f64 {
        (0..cont.len())
            .map(|c| {
                let s: f64 = dados[ofs[c]..ofs[c + 1]].iter().map(|&u| e[u as usize] / eff[u as usize]).sum();
                if s > 0.0 {
                    cont[c] * s.ln()
                } else {
                    f64::NEG_INFINITY
                }
            })
            .sum()
    };
    let mut melhor = f64::NEG_INFINITY;
    let passo = 0.002;
    let mut x = 0.0;
    while x <= 1.0 {
        let mut y = 0.0;
        while x + y <= 1.0 {
            melhor = melhor.max(ll([x, y, 1.0 - x - y]));
            y += passo;
        }
        x += passo;
    }
    t.check(logl >= melhor - 1e-9, "EM, 3 transcripts: no point of a 0.002 grid has a higher likelihood");
    // regras do kallisto: um transcrito sem leitura propria vai a zero; rodadas >= 51
    // {A,B} 50 reads, {A} 100: B shrinks by 1/3 a round, and below 1e-8 it is set to 0
    let (a, rodadas) = em(2, &[0, 2, 3], &[0, 1, 0], &[50.0, 100.0], &[1000.0, 1000.0], &EmOptions::default());
    t.check(a[1] == 0.0 && (a[0] - 150.0).abs() < 1e-9, "EM: all reads go to the transcript with evidence of its own; the other is set to 0");
    // {A,B} 100, {A} 50: B shrinks by 2/3 a round and is still 5e-8 after 52 rounds: kept, as kallisto keeps it
    let (a, _) = em(2, &[0, 2, 3], &[0, 1, 0], &[100.0, 50.0], &[1000.0, 1000.0], &EmOptions::default());
    t.check(a[1] > 0.0 && a[1] < 1e-6 && (a[0] + a[1] - 150.0).abs() < 1e-9, "EM: a count above kallisto's zeroing limit is kept");
    t.check(rodadas >= 52, &format!("EM: at least 50 rounds, then a final one ({rodadas})"));
    let (_, r1) = em(1, &[0, 1], &[0], &[10.0], &[100.0], &EmOptions::default());
    t.eq(&r1, &52, "EM: kallisto's count of rounds for an easy case (52)");
    // TPM
    let tp = tpm(&[10.0, 20.0, 0.0], &[100.0, 400.0, 50.0]);
    t.perto(tp[0], 1e6 * (0.1 / 0.15), 1e-12, "TPM of a transcript");
    t.perto(tp.iter().sum::<f64>(), 1e6, 1e-12, "TPM add up to a million");
}

fn testa_comprimentos(t: &mut T) {
    let fld = gaussian_fld(200.0, 20.0);
    let ls: Vec<u32> = vec![0, 1, 50, 150, 200, 250, 999, 1000, 5000];
    let e = effective_lengths(&ls, &fld);
    for (i, &l) in ls.iter().enumerate() {
        let d = eff_direto(l as usize, 200.0, 20.0);
        t.perto(e[i], d.max(1.0), 1e-9, &format!("effective length of {l} bp, single-end N(200, 20)"));
    }
    // a distribuicao vazia: comprimento + 1
    let e = effective_lengths(&[100], &vec![0.0; MAX_FRAG_LEN]);
    t.eq(&e[0], &101.0, "no fragment length known: length + 1");
}

/// Quantification of simulated samples.
fn testa_quantificacao(t: &mut T, txs: &[Tx], idx_path: &str) {
    let ing = Ingenuo::novo(txs, 31);
    let mut rng = Rng(99);
    let seqs: Vec<&[u8]> = txs.iter().map(|x| x.seq.as_slice()).collect();
    let theta = abundancias(&mut rng, txs.len(), 0.2);
    let (gene_de, gnomes) = genes_de(txs);
    for &paired in &[true, false] {
        let o = SimOpt { n: 200_000, read_len: 100, media: 200.0, dp: 25.0, erro: 0.005, paired };
        let sim = simula(&seqs, &theta, &mut rng, &o);
        let nome = if paired { "pe" } else { "se" };
        let r1 = t.p(&format!("{nome}_1.fastq.gz"));
        escreve_fastq(&r1, &sim.r1, "/1");
        let r2 = if paired {
            let p = t.p(&format!("{nome}_2.fastq.gz"));
            escreve_fastq(&p, &sim.r2, "/2");
            Some(p)
        } else {
            None
        };
        let mut opt = QuantOptions { keep_class_counts: true, ..Default::default() };
        if !paired {
            opt.fragment_length = Some(200.0);
            opt.fragment_sd = Some(25.0);
        }
        let saida = t.p(&format!("quant_{nome}"));
        let amostra = QuantSample { name: nome.into(), r1: r1.clone(), r2: r2.clone(), condition: None };
        let inicio = Instant::now();
        let Some(r) = t.ok_ou(quantify_samples(idx_path, &[amostra.clone()], &saida, &opt), &format!("quantify {nome}")) else {
            continue;
        };
        let seg = inicio.elapsed().as_secs_f64();
        let s = &r.samples[0];
        match serde_json::from_str::<serde_json::Value>(&r.json()) {
            Ok(v) => {
                let ok = v["samples"][0]["reads"].as_u64() == Some(o.n as u64)
                    && v["gene_counts"][0].as_array().map(|a| a.len()) == Some(r.gene_ids.len())
                    && v["files"]["gene_counts"].is_string()
                    && v["samples"][0]["overhang"].as_u64() == Some(0);
                t.check(ok, &format!("{nome}: the result JSON parses, with its keys"));
            }
            Err(e) => {
                t.check(false, &format!("{nome}: result JSON: {e}"));
            }
        }
        // contagens por classe = as da pseudoalinhamento ingenua
        let mut esperado: HashMap<Vec<u32>, u64> = HashMap::new();
        let mut st = [0u64; 4];
        let mut frags: Vec<i64> = Vec::new();
        let mut frags_verdade: Vec<i64> = Vec::new();
        for i in 0..sim.r1.len() {
            let mates: Vec<&[u8]> = if paired { vec![&sim.r1[i], &sim.r2[i]] } else { vec![&sim.r1[i]] };
            let (c, ts) = ing.classe(&mates);
            st[c as usize] += 1;
            if c == 3 {
                if ts.len() == 1 && frags.len() < 10_000 && paired {
                    if let Some(f) = frag_reconstruido(&ing, &txs[ts[0] as usize].seq, &sim.r1[i], &sim.r2[i]) {
                        frags.push(f);
                        frags_verdade.push(sim.frag[i] as i64);
                    }
                }
                *esperado.entry(ts).or_insert(0) += 1;
            }
        }
        let obtido: HashMap<Vec<u32>, u64> = s.class_counts.iter().cloned().collect();
        t.check(obtido == esperado, &format!("{nome}: reads per equivalence class equal the naive pseudoalignment ({} classes)", esperado.len()));
        t.eq(&(s.reads, s.pseudoaligned, s.no_kmer, s.conflicting, s.too_short), &(o.n as u64, st[3], st[1], st[2], st[0]), &format!("{nome}: read counts"));
        let unicos: u64 = esperado.iter().filter(|(k, _)| k.len() == 1).map(|(_, v)| v).sum();
        t.eq(&s.unique, &unicos, &format!("{nome}: unique = reads on a single transcript"));
        if paired {
            // comprimento de fragmento: os 10.000 primeiros pares unicos que se deixam
            // posicionar, reconstruidos aqui por busca de texto no transcrito
            let ms = |v: &[i64]| {
                let n = v.len() as f64;
                let m = v.iter().map(|&f| f as f64).sum::<f64>() / n;
                (m, (v.iter().map(|&f| (f as f64 - m).powi(2)).sum::<f64>() / n).sqrt())
            };
            let (m, sd) = ms(&frags);
            let (mv, sdv) = ms(&frags_verdade);
            t.perto(s.fragment_length_mean, m, 1e-9, "pe: the fragment length mean is that of the first 10,000 placeable unique pairs");
            t.perto(s.fragment_length_sd, sd, 1e-9, "pe: and so is the standard deviation");
            t.eq(&s.fragment_pairs, &10_000, "pe: 10,000 pairs used");
            let errados = frags.iter().zip(&frags_verdade).filter(|(a, b)| a != b).count();
            println!("  pe: fragment length {m:.2} +/- {sd:.2} (true lengths of the same pairs {mv:.2} +/- {sdv:.2}; {errados} of 10,000 placed on another isoform by a sequencing error)");
            t.check(errados <= 10, &format!("pe: at most 10 of the 10,000 pairs measured on the wrong isoform ({errados})"));
            t.perto(s.fragment_length_mean, mv, 0.002, "pe: the estimated mean is within 0.2% of the true one");
        }
        // exatidao
        let verdade: Vec<f64> = sim.contagem.iter().map(|&c| c as f64).collect();
        let tot_theta: f64 = theta.iter().zip(&seqs).filter(|(_, s)| s.len() >= 31).map(|(x, _)| x).sum();
        let tpm_verdade: Vec<f64> = theta.iter().zip(&seqs).map(|(x, s)| if s.len() >= 31 { x / tot_theta * 1e6 } else { 0.0 }).collect();
        // TPM verdadeiro da amostra: contagem sorteada / comprimento efetivo do simulador
        let effv: Vec<f64> = seqs.iter().map(|x| eff_direto(x.len(), o.media, o.dp)).collect();
        let rho: Vec<f64> = verdade.iter().zip(&effv).map(|(c, e)| c / e).collect();
        let srho: f64 = rho.iter().sum();
        let tpm_amostra: Vec<f64> = rho.iter().map(|r| r / srho * 1e6).collect();
        let mt = metricas(&s.est_counts, &verdade);
        let mtp = metricas(&s.tpm, &tpm_amostra);
        let mtp_pop = metricas(&s.tpm, &tpm_verdade);
        let ng = gnomes.len();
        let mg = metricas(&por_gene(&s.est_counts, &gene_de, ng), &por_gene(&verdade, &gene_de, ng));
        println!("  {nome}: {} reads, {:.1}% pseudoaligned, {:.2} s", o.n, 100.0 * s.pseudoaligned as f64 / o.n as f64, seg);
        mostra("transcripts, est. counts vs truth", &mt);
        mostra("transcripts, TPM vs sampled truth", &mtp);
        mostra("transcripts, TPM vs abundance drawn", &mtp_pop);
        mostra("genes, est. counts vs truth", &mg);
        let (pl, sp, md) = if paired { (0.93, 0.97, 0.13) } else { (0.85, 0.92, 0.20) };
        t.maior(mt.pearson_log, pl, &format!("{nome}: transcript counts, Pearson on log"));
        t.maior(mt.spearman, sp, &format!("{nome}: transcript counts, Spearman"));
        t.menor(mt.mard, md, &format!("{nome}: transcript counts, MARD"));
        t.maior(mtp.spearman, sp, &format!("{nome}: transcript TPM, Spearman"));
        t.menor(mtp.mard, md, &format!("{nome}: transcript TPM, MARD"));
        t.maior(mg.pearson_log, 0.995, &format!("{nome}: gene counts, Pearson on log"));
        t.menor(mg.mard, 0.02, &format!("{nome}: gene counts, MARD"));
        // as tabelas
        let ab = le_tsv(&format!("{saida}/{nome}/abundance.tsv"));
        t.eq(&ab[0].join("\t"), &"target_id\tlength\teff_length\test_counts\ttpm".to_string(), &format!("{nome}: abundance.tsv has kallisto's columns"));
        t.eq(&(ab.len() - 1), &txs.len(), &format!("{nome}: one line per transcript"));
        let soma_tpm: f64 = ab[1..].iter().map(|l| l[4].parse::<f64>().unwrap()).sum();
        let soma_c: f64 = ab[1..].iter().map(|l| l[3].parse::<f64>().unwrap()).sum();
        t.perto(soma_tpm, 1e6, 1e-9, &format!("{nome}: TPM add up to a million"));
        t.perto(soma_c, s.pseudoaligned as f64, 1e-6, &format!("{nome}: estimated counts add up to the pseudoaligned reads"));
        let gc = le_tsv(&format!("{saida}/gene_counts.tsv"));
        t.eq(&gc[0], &vec!["gene".to_string(), nome.to_string()], &format!("{nome}: gene_counts.tsv header"));
        let ok = gc[1..].iter().enumerate().all(|(g, l)| l[1].parse::<u64>().map_or(false, |x| x == s.gene_counts[g].round() as u64));
        t.check(ok, &format!("{nome}: gene_counts.tsv holds the rounded gene counts"));
        let info: String = std::fs::read_to_string(format!("{saida}/{nome}/run_info.json")).unwrap_or_default();
        t.check(info.contains("\"pseudoaligned\":") && info.contains(&format!("\"sample\":\"{nome}\"")), &format!("{nome}: run_info.json"));
        if !paired {
            // a regra do kallisto para leituras simples (sem --single-overhang),
            // reconstruida aqui: em cada transcrito da classe, a leitura e posta pelo
            // primeiro k-mer achado (busca de texto); fica o transcrito em que um
            // fragmento de 200 pb, a partir da leitura, cabe
            let mut esperado_f: HashMap<Vec<u32>, u64> = HashMap::new();
            let mut fora = 0u64;
            for r in &sim.r1 {
                let (c, ts) = ing.classe(&[r]);
                if c != 3 {
                    continue;
                }
                let (o1, km) = r.windows(31).enumerate().find(|(_, w)| acgt(w) && ing.mapa.contains_key(&canon(w))).map(|(o, w)| (o, w.to_vec())).unwrap();
                let kmr = rc(&km);
                let fica: Vec<u32> = ts
                    .iter()
                    .copied()
                    .filter(|&x| {
                        // cai o transcrito em que alguma ocorrencia, na orientacao da
                        // leitura, deixaria o fragmento passar da ponta
                        let tx = &txs[x as usize].seq;
                        let l = tx.len() as i64;
                        !(0..=tx.len() - 31).any(|q| {
                            (tx[q..q + 31] == km[..] && q as i64 - o1 as i64 + 1 + 200 > l)
                                || (tx[q..q + 31] == kmr[..] && q as i64 + 31 + o1 as i64 - 200 < 0)
                        })
                    })
                    .collect();
                if fica.is_empty() {
                    fora += 1;
                } else {
                    *esperado_f.entry(fica).or_insert(0) += 1;
                }
            }
            let mut of = opt.clone();
            of.single_overhang = false;
            if let Some(r) = t.ok_ou(quantify_samples(idx_path, &[amostra.clone()], &t.p("quant_se_filter"), &of), "single-end with kallisto's fragment rule") {
                let obtido: HashMap<Vec<u32>, u64> = r.samples[0].class_counts.iter().cloned().collect();
                t.check(obtido == esperado_f, "se, single_overhang off: reads per class equal kallisto's rule applied here by text search");
                t.eq(&r.samples[0].overhang, &fora, "se, single_overhang off: the reads dropped by the rule");
                let mf = metricas(&r.samples[0].est_counts, &verdade);
                let mgf = metricas(&por_gene(&r.samples[0].est_counts, &gene_de, ng), &por_gene(&verdade, &gene_de, ng));
                println!("  se with kallisto's single-end rule: {fora} reads dropped");
                mostra("transcripts, est. counts vs truth", &mf);
                mostra("genes, est. counts vs truth", &mgf);
            }
        }
        if paired {
            // 1 thread = 8 threads, byte a byte
            let saida1 = t.p("quant_pe_1t");
            let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
            if t.ok_ou(pool.install(|| quantify_samples(idx_path, &[amostra.clone()], &saida1, &opt)), "quantify with 1 thread").is_some() {
                let mut iguais = true;
                for f in ["gene_counts.tsv", "gene_tpm.tsv", "transcript_counts.tsv", "transcript_tpm.tsv", "pe/abundance.tsv", "pe/abundance_genes.tsv"] {
                    iguais &= std::fs::read(format!("{saida}/{f}")).unwrap() == std::fs::read(format!("{saida1}/{f}")).unwrap();
                }
                t.check(iguais, "1 thread and all threads write the same tables, byte for byte");
            }
            // max_reads
            let mut o2 = opt.clone();
            o2.max_reads = Some(1000);
            if let Some(r) = t.ok_ou(quantify_samples(idx_path, &[amostra.clone()], &t.p("quant_pe_max"), &o2), "max_reads") {
                t.eq(&r.samples[0].reads, &1000, "max_reads = 1000 reads 1000 pairs");
            }
            // texto puro = gzip
            let p1 = t.p("pe_plain_1.fastq");
            let p2 = t.p("pe_plain_2.fastq");
            escreve_fastq(&p1, &sim.r1, "/1");
            escreve_fastq(&p2, &sim.r2, "/2");
            let a2 = QuantSample { name: "pe".into(), r1: p1, r2: Some(p2), condition: None };
            if let Some(r) = t.ok_ou(quantify_samples(idx_path, &[a2], &t.p("quant_pe_plain"), &opt), "plain FASTQ") {
                t.check(r.samples[0].est_counts == s.est_counts, "plain and gzip FASTQ give the same estimates");
            }
            // comprimento dado em vez de estimado
            let mut o3 = opt.clone();
            o3.fragment_length = Some(200.0);
            o3.fragment_sd = Some(25.0);
            if let Some(r) = t.ok_ou(quantify_samples(idx_path, &[amostra.clone()], &t.p("quant_pe_given"), &o3), "pairs with a given fragment length") {
                t.eq(&r.samples[0].fragment_length_source, &"given".to_string(), "a given fragment length is used for pairs");
                let m = metricas(&r.samples[0].est_counts, &s.est_counts);
                t.maior(m.pearson_log, 0.999, "estimated and given fragment lengths give nearly the same counts");
            }
            // mapa de genes na quantificacao
            let pm = t.p("quant_mapa.tsv");
            let linhas: String = txs.iter().map(|x| format!("{}\tNEW_{}\n", x.name, x.gene)).collect();
            grava(&pm, &linhas);
            let mut o4 = opt.clone();
            o4.gene_map = Some(pm);
            o4.max_reads = Some(2000);
            if let Some(r) = t.ok_ou(quantify_samples(idx_path, &[amostra.clone()], &t.p("quant_pe_map"), &o4), "gene map at quantification") {
                t.check(r.gene_ids.iter().all(|g| g.starts_with("NEW_")), "a gene map given to quantify_samples replaces the one in the index");
            }
        }
    }
}

/// The fragment length of a pair placed on transcript `tx` by the first k-mer
/// of each mate found in the naive map, searched as text in the transcript:
/// mates on opposite strands; |start of the forward one - end of the reverse
/// one|; (0, 1000).
fn frag_reconstruido(ing: &Ingenuo, tx: &[u8], r1: &[u8], r2: &[u8]) -> Option<i64> {
    let k = ing.k;
    let primeiro = |r: &[u8]| -> Option<(usize, Vec<u8>)> {
        if r.len() < k {
            return None;
        }
        r.windows(k).enumerate().find(|(_, w)| acgt(w) && ing.mapa.contains_key(&canon(w))).map(|(o, w)| (o, w.to_vec()))
    };
    let lugar = |o: usize, km: &[u8]| -> Option<(i64, bool)> {
        let r = rc(km);
        for q in 0..=tx.len().checked_sub(k)? {
            if &tx[q..q + k] == km {
                return Some((q as i64 - o as i64, true));
            }
            if tx[q..q + k] == r[..] {
                return Some((q as i64 + k as i64 + o as i64, false));
            }
        }
        None
    };
    let (o1, k1) = primeiro(r1)?;
    let (o2, k2) = primeiro(r2)?;
    let (a, da) = lugar(o1, &k1)?;
    let (b, db) = lugar(o2, &k2)?;
    let tl = (a - b).abs();
    (da != db && tl > 0 && tl < 1000).then_some(tl)
}

/// 3 + 3 samples with planted fold changes -> gene_counts.tsv ->
/// differential_expression.
fn cadeia_de(txs: &[Tx], idx_path: &str, pasta: &str, fragmentos: usize, semente: u64, plantados: usize, mostrar: bool) -> Result<(usize, usize, usize, usize, usize), String> {
    let mut rng = Rng(semente);
    let seqs: Vec<&[u8]> = txs.iter().map(|x| x.seq.as_slice()).collect();
    let (gene_de, gnomes) = genes_de(txs);
    let ng = gnomes.len();
    let base = abundancias(&mut rng, txs.len(), 0.1);
    // genes expressos: soma de theta acima da mediana dos que tem algo
    let mut expr = vec![0.0; ng];
    for (t, &x) in base.iter().enumerate() {
        if seqs[t].len() >= 200 {
            expr[gene_de[t]] += x;
        }
    }
    let mut candidatos: Vec<usize> = (0..ng).filter(|&g| expr[g] > 0.0).collect();
    candidatos.sort_by(|&a, &b| expr[b].total_cmp(&expr[a]));
    candidatos.truncate(candidatos.len() / 2);
    let mut sobe = BTreeSet::new();
    let mut desce = BTreeSet::new();
    while sobe.len() < plantados || desce.len() < plantados {
        let g = candidatos[rng.abaixo(candidatos.len() as u64) as usize];
        if sobe.contains(&g) || desce.contains(&g) {
            continue;
        }
        if sobe.len() < plantados {
            sobe.insert(g);
        } else {
            desce.insert(g);
        }
    }
    std::fs::create_dir_all(pasta).map_err(|e| e.to_string())?;
    let mut amostras = Vec::new();
    let mut verdade_genes: Vec<Vec<u64>> = Vec::new();
    for s in 0..6 {
        let tratado = s >= 3;
        let ruido_gene: Vec<f64> = (0..ng).map(|_| (0.1 * rng.normal()).exp()).collect();
        let theta: Vec<f64> = (0..txs.len())
            .map(|t| {
                let g = gene_de[t];
                let fc = if tratado && sobe.contains(&g) {
                    4.0
                } else if tratado && desce.contains(&g) {
                    0.25
                } else {
                    1.0
                };
                base[t] * fc * ruido_gene[g]
            })
            .collect();
        let n = (fragmentos as f64 * (0.8 + 0.4 * rng.unif())) as usize;
        let sim = simula(&seqs, &theta, &mut rng, &SimOpt { n, read_len: 100, media: 200.0, dp: 25.0, erro: 0.005, paired: true });
        let nome = format!("{}{}", if tratado { "treated" } else { "control" }, s % 3 + 1);
        let r1 = format!("{pasta}/{nome}_1.fastq.gz");
        let r2 = format!("{pasta}/{nome}_2.fastq.gz");
        escreve_fastq(&r1, &sim.r1, "/1");
        escreve_fastq(&r2, &sim.r2, "/2");
        let mut vg = vec![0u64; ng];
        for (t, &c) in sim.contagem.iter().enumerate() {
            vg[gene_de[t]] += c;
        }
        verdade_genes.push(vg);
        amostras.push(QuantSample { name: nome, r1, r2: Some(r2), condition: Some(if tratado { "treated".into() } else { "control".into() }) });
    }
    let saida = format!("{pasta}/quant");
    let r = quantify_samples(idx_path, &amostras, &saida, &QuantOptions::default())?;
    if mostrar {
        print!("{}", r.resumo());
    }
    let de = expressao::differential_expression(
        &format!("{saida}/gene_counts.tsv"),
        &format!("{saida}/samples.tsv"),
        Some(&format!("{saida}/de.tsv")),
        &expressao::DeOptions { condition: "condition".into(), reference: "control".into(), ..Default::default() },
    )?;
    if mostrar {
        print!("{}", de.resumo());
    }
    let nome_para_g: HashMap<&str, usize> = gnomes.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let (mut ac_sobe, mut ac_desce, mut falsos, mut testados) = (0, 0, 0, 0);
    for g in &de.genes {
        let Some(&gi) = nome_para_g.get(g.gene.as_str()) else { continue };
        if !g.pvalue.is_nan() {
            testados += 1;
        }
        let sig = g.padj < 0.1;
        if sobe.contains(&gi) {
            if sig && g.log2_fold_change > 0.0 {
                ac_sobe += 1;
            }
        } else if desce.contains(&gi) {
            if sig && g.log2_fold_change < 0.0 {
                ac_desce += 1;
            }
        } else if sig {
            falsos += 1;
            if mostrar {
                println!("   false positive {} log2FC {:.2} padj {:.2e}", g.gene, g.log2_fold_change, g.padj);
            }
        }
    }
    if mostrar {
        // a mesma analise sobre as contagens VERDADEIRAS (fragmentos simulados por gene):
        // o que sobra de falso positivo ali e da estatistica, nao da quantificacao
        let pv = format!("{saida}/truth_gene_counts.tsv");
        let mut txt = format!("gene{}\n", amostras.iter().map(|a| format!("\t{}", a.name)).collect::<String>());
        for g in 0..ng {
            txt.push_str(&gnomes[g]);
            for vg in &verdade_genes {
                txt.push_str(&format!("\t{}", vg[g]));
            }
            txt.push('\n');
        }
        std::fs::write(&pv, txt).map_err(|e| e.to_string())?;
        let dv = expressao::differential_expression(
            &pv,
            &format!("{saida}/samples.tsv"),
            None,
            &expressao::DeOptions { condition: "condition".into(), reference: "control".into(), ..Default::default() },
        )?;
        let (mut vs, mut vd, mut vf) = (0, 0, 0);
        let mut fs = Vec::new();
        for g in &dv.genes {
            let Some(&gi) = nome_para_g.get(g.gene.as_str()) else { continue };
            let sig = g.padj < 0.1;
            if sobe.contains(&gi) {
                vs += (sig && g.log2_fold_change > 0.0) as usize;
            } else if desce.contains(&gi) {
                vd += (sig && g.log2_fold_change < 0.0) as usize;
            } else if sig {
                vf += 1;
                fs.push(g.gene.clone());
            }
        }
        println!("   on the TRUE gene counts of the same samples: {vs} up and {vd} down found, {vf} false positives ({})", fs.join(", "));
        let falt: Vec<String> = de
            .genes
            .iter()
            .filter(|g| nome_para_g.get(g.gene.as_str()).map_or(false, |gi| (sobe.contains(gi) || desce.contains(gi)) && !(g.padj < 0.1)))
            .map(|g| format!("{} (log2FC {:.2}, baseMean {:.0}, padj {:.2e})", g.gene, g.log2_fold_change, g.base_mean, g.padj))
            .collect();
        if !falt.is_empty() {
            println!("   planted but not found: {}", falt.join("; "));
        }
    }
    Ok((ac_sobe, ac_desce, falsos, testados, plantados))
}

fn testa_cadeia(t: &mut T, txs: &[Tx], idx_path: &str) {
    let pasta = t.p("de_chain");
    if let Some((s, d, f, testados, p)) = t.ok_ou(cadeia_de(txs, idx_path, &pasta, 150_000, 31337, 10, false), "DE chain") {
        println!("  DE chain (synthetic): planted 10 up + 10 down x4: found {s} up, {d} down; {f} false positives among {} other tested genes", testados - 2 * p);
        t.check(s + d >= 18, &format!("DE chain: at least 18 of the 20 planted genes found with the right sign ({})", s + d));
        t.check(f <= 3, &format!("DE chain: at most 3 false positives ({f})"));
    }
    let samples = le_tsv(&format!("{pasta}/quant/samples.tsv"));
    t.eq(&samples[0], &vec!["sample".to_string(), "condition".to_string()], "samples.tsv is written when every sample has a condition");
}

fn testa_erros(t: &mut T, txs: &[Tx], idx_path: &str) {
    let pfa = t.p("erros.fa");
    t.erro(build_transcript_index(&pfa, &t.p("x.idx"), &IndexOptions::default()), &["Error opening"], "a missing FASTA");
    grava(&pfa, ">a\nACGT\n>a\nACGTACGT\n");
    t.erro(build_transcript_index(&pfa, &t.p("x.idx"), &IndexOptions::default()), &["\"a\" was already used on line 1", "make_unique"], "repeated names");
    grava(&pfa, &format!(">a\n{}\n>a\n{}\n", "ACGT".repeat(20), "TTGA".repeat(20)));
    if let Some(r) = t.ok_ou(build_transcript_index(&pfa, &t.p("x.idx"), &IndexOptions { make_unique: true, ..Default::default() }), "make_unique") {
        let i = TranscriptIndex::load(&r.index).unwrap();
        t.eq(&i.names, &vec!["a".to_string(), "a_2".to_string()], "make_unique renames the repeat a_2");
    }
    grava(&pfa, "@read1\nACGT\n+\nIIII\n");
    t.erro(build_transcript_index(&pfa, &t.p("x.idx"), &IndexOptions::default()), &["looks like a FASTQ"], "reads given as the transcriptome");
    grava(&pfa, ">a\nACGTNNNN\n");
    t.erro(build_transcript_index(&pfa, &t.p("x.idx"), &IndexOptions::default()), &["no k-mer of length 31"], "a FASTA without a single k-mer");
    t.erro(build_transcript_index(&pfa, &t.p("x.idx"), &IndexOptions { k: 30, ..Default::default() }), &["odd number from 7 to 31"], "an even k");
    t.erro(build_transcript_index(&pfa, &pfa, &IndexOptions::default()), &["the input file"], "output = input");
    // quantificacao
    let mut rng = Rng(3);
    let seqs: Vec<&[u8]> = txs.iter().map(|x| x.seq.as_slice()).collect();
    let theta = vec![1.0; txs.len()];
    let sim = simula(&seqs, &theta, &mut rng, &SimOpt { n: 300, read_len: 80, media: 200.0, dp: 20.0, erro: 0.0, paired: true });
    let r1 = t.p("err_1.fastq");
    let r2 = t.p("err_2.fastq");
    escreve_fastq(&r1, &sim.r1, "/1");
    escreve_fastq(&r2, &sim.r2[..299], "/2");
    let saida = t.p("q_err");
    let a = |r1: &str, r2: Option<&str>| QuantSample { name: "s".into(), r1: r1.into(), r2: r2.map(str::to_string), condition: None };
    t.erro(quantify_samples(idx_path, &[a(&r1, Some(&r2))], &saida, &QuantOptions::default()), &["ends after 299 reads", "same reads in the same order"], "R2 shorter than R1");
    // nomes fora de passo
    let mut trocado = sim.r2.clone();
    trocado.swap(10, 11);
    escreve_fastq(&r2, &trocado, "/2");
    let r2b = t.p("err_2b.fastq");
    {
        // nomes trocados: o arquivo R2 com os nomes das leituras 10 e 11 invertidos
        let txt = std::fs::read_to_string(&r2).unwrap();
        let txt = txt.replacen("@r10/2", "@TMP", 1).replacen("@r11/2", "@r10/2", 1).replacen("@TMP", "@r11/2", 1);
        grava(&r2b, &txt);
    }
    t.erro(quantify_samples(idx_path, &[a(&r1, Some(&r2b))], &saida, &QuantOptions::default()), &["Pair 11", "\"r10/1\"", "\"r11/2\"", "same order"], "mates out of step");
    let mut o = QuantOptions::default();
    o.check_pair_names = false;
    t.check(quantify_samples(idx_path, &[a(&r1, Some(&r2b))], &saida, &o).is_ok(), "check_pair_names=False skips the name test");
    // nomes do SRA: com --readids (SRR1.7.1 / SRR1.7.2) e sem (iguais nos dois arquivos)
    for (s1, s2) in [(".1", ".2"), ("", "")] {
        let p1 = t.p("sra_1.fastq");
        let p2 = t.p("sra_2.fastq");
        let nomeia = |reads: &[Vec<u8>], suf: &str| -> String {
            reads.iter().enumerate().map(|(i, r)| format!("@SRR1.{i}{suf} length={}\n{}\n+\n{}\n", r.len(), String::from_utf8_lossy(r), "I".repeat(r.len()))).collect()
        };
        grava(&p1, &nomeia(&sim.r1, s1));
        grava(&p2, &nomeia(&sim.r2, s2));
        t.check(quantify_samples(idx_path, &[a(&p1, Some(&p2))], &saida, &QuantOptions::default()).is_ok(), &format!("mates named SRR1.N{s1} / SRR1.N{s2} are accepted"));
        // e um descompasso continua sendo visto: R2 comeca na leitura 1
        grava(&p2, &nomeia(&sim.r2, s2).replacen(&format!("@SRR1.0{s2} "), &format!("@SRR1.9999{s2} "), 1));
        t.erro(quantify_samples(idx_path, &[a(&p1, Some(&p2))], &saida, &QuantOptions::default()), &["Pair 1 of", "SRR1.9999"], &format!("SRR1.N{s1} / SRR1.N{s2}: a different name is still caught"));
    }
    t.erro(quantify_samples(idx_path, &[a(&r1, None)], &saida, &QuantOptions::default()), &["single-end", "fragment_length and fragment_sd"], "single-end without a fragment length");
    let mut o = QuantOptions::default();
    o.fragment_length = Some(200.0);
    t.erro(quantify_samples(idx_path, &[a(&r1, None)], &saida, &o), &["both fragment_length and fragment_sd"], "a mean without the sd");
    t.erro(quantify_samples(idx_path, &[a(&r1, Some(&r1))], &saida, &QuantOptions::default()), &["R1 and R2 are the same file"], "R1 = R2");
    t.erro(quantify_samples(idx_path, &[a(&r1, Some(&r2)), a(&r1, Some(&r2))], &saida, &QuantOptions::default()), &["both named \"s\""], "two samples with one name");
    let mut b = a(&r1, Some(&r2));
    b.name = "a/b".into();
    t.erro(quantify_samples(idx_path, &[b], &saida, &QuantOptions::default()), &["cannot be used"], "a sample name with a slash");
    t.erro(quantify_samples(idx_path, &[a(&t.p("nada.fq"), None)], &saida, &QuantOptions::default()), &["does not exist"], "a missing FASTQ");
    t.erro(quantify_samples(idx_path, &[], &saida, &QuantOptions::default()), &["No sample given"], "no sample");
    let ruim = t.p("ruim.fastq");
    grava(&ruim, "@r1\nACGTACGTAC\n+\nIIII\n");
    let mut o = QuantOptions::default();
    o.fragment_length = Some(200.0);
    o.fragment_sd = Some(20.0);
    t.erro(quantify_samples(idx_path, &[a(&ruim, None)], &saida, &o), &["line 1", "10 bases and the quality line 4"], "a FASTQ with a short quality line");
    grava(&ruim, "@r1\nACGTACGTAC\n+\nIIIIIIIIII\nACGT\n");
    t.erro(quantify_samples(idx_path, &[a(&ruim, None)], &saida, &o), &["line 5", "starts with '@'"], "a FASTQ with a broken record");
    grava(&ruim, "@r1\nACGTACGTAC\n+\n");
    t.erro(quantify_samples(idx_path, &[a(&ruim, None)], &saida, &o), &["incomplete", "truncated"], "a truncated FASTQ");
    grava(&ruim, ">r1\nACGTACGTAC\n");
    t.erro(quantify_samples(idx_path, &[a(&ruim, None)], &saida, &o), &["this is a FASTA file"], "reads in FASTA");
    t.erro(quantify_samples(&pfa, &[a(&r1, None)], &saida, &o), &["not an MMGEN transcript index"], "a FASTA given as the index");
    // nome automatico
    t.eq(&nome_de_arquivo("dir/ctrl1_R1_001.fastq.gz"), &"ctrl1".to_string(), "sample name from ctrl1_R1_001.fastq.gz");
    t.eq(&nome_de_arquivo("C:\\x\\amostra_1.fq"), &"amostra".to_string(), "sample name from amostra_1.fq");
}

fn testa_explicacao(t: &mut T, txs: &[Tx], idx_path: &str) {
    let leitura = String::from_utf8(txs[0].seq[0..60].to_vec()).unwrap();
    match explain_pseudoalignment(idx_path, &leitura, None) {
        Ok(j) => {
            let v: serde_json::Value = serde_json::from_str(&j).unwrap();
            t.eq(&v["mates"][0]["kmers"].as_array().map(|a| a.len()), &Some(30), "explain: the 30 k-mers of a 60-bp read");
            t.check(v["transcripts"].as_array().map_or(false, |a| a.iter().any(|x| x == "T0000_1")), "explain: the read's transcripts");
            t.check(v["status"] == "unique" || v["status"] == "multi", "explain: the status");
        }
        Err(e) => {
            t.check(false, &format!("explain: {e}"));
        }
    }
}

// ================================================================ modes for real data

/// Transcripts from a genome and its GFF3 (RefSeq: exons whose Parent is a
/// transcript `rna-...`), spliced and reverse-complemented on the - strand.
fn modo_refseq(args: &[String]) {
    let (genoma, gff, pasta) = (&args[0], &args[1], &args[2]);
    std::fs::create_dir_all(pasta).unwrap();
    let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
    {
        let r = externo::abre_leitura(genoma).unwrap();
        let mut nome = String::new();
        for l in r.lines() {
            let l = l.unwrap();
            if let Some(h) = l.strip_prefix('>') {
                nome = h.split_whitespace().next().unwrap_or("").to_string();
                seqs.insert(nome.clone(), Vec::new());
            } else if let Some(s) = seqs.get_mut(&nome) {
                s.extend(l.trim().bytes().map(|c| c.to_ascii_uppercase()));
            }
        }
    }
    println!("genome: {} sequences, {} bp", seqs.len(), seqs.values().map(Vec::len).sum::<usize>());
    // transcrito -> (seqid, fita, exons, gene)
    let mut tx: HashMap<String, (String, char, Vec<(usize, usize)>, String)> = HashMap::new();
    let mut ordem: Vec<String> = Vec::new();
    let r = externo::abre_leitura(gff).unwrap();
    for l in r.lines() {
        let l = l.unwrap();
        if l.starts_with('#') {
            continue;
        }
        let c: Vec<&str> = l.split('\t').collect();
        if c.len() < 9 || c[2] != "exon" || !seqs.contains_key(c[0]) {
            continue;
        }
        let at: HashMap<&str, &str> = c[8].split(';').filter_map(|p| p.split_once('=')).collect();
        let Some(pai) = at.get("Parent") else { continue };
        let Some(id) = pai.strip_prefix("rna-") else { continue };
        let gene = at.get("gene").copied().unwrap_or(id).to_string();
        let (a, b): (usize, usize) = (c[3].parse().unwrap(), c[4].parse().unwrap());
        let e = tx.entry(id.to_string()).or_insert_with(|| {
            ordem.push(id.to_string());
            (c[0].to_string(), c[6].chars().next().unwrap_or('+'), Vec::new(), gene)
        });
        e.2.push((a, b));
    }
    let mut fa = std::io::BufWriter::new(std::fs::File::create(format!("{pasta}/transcripts.fa")).unwrap());
    let mut mapa = std::io::BufWriter::new(std::fs::File::create(format!("{pasta}/gene_map.tsv")).unwrap());
    let (mut n, mut bases) = (0, 0);
    let mut genes = BTreeSet::new();
    for id in &ordem {
        let (sq, fita, exons, gene) = &tx[id];
        let mut ex = exons.clone();
        ex.sort();
        let g = &seqs[sq];
        let mut s: Vec<u8> = Vec::new();
        for (a, b) in ex {
            s.extend_from_slice(&g[a - 1..b.min(g.len())]);
        }
        if *fita == '-' {
            s = rc(&s);
        }
        writeln!(fa, ">{id} gene:{gene} gene_symbol:{gene}").unwrap();
        for c in s.chunks(70) {
            fa.write_all(c).unwrap();
            fa.write_all(b"\n").unwrap();
        }
        writeln!(mapa, "{id}\t{gene}\t{gene}").unwrap();
        n += 1;
        bases += s.len();
        genes.insert(gene.clone());
    }
    println!("{n} transcripts of {} genes, {bases} bases -> {pasta}/transcripts.fa, gene_map.tsv", genes.len());
}

fn modo_sintetico(args: &[String]) {
    let pasta = &args[0];
    let genes: usize = args[1].parse().unwrap();
    let semente: u64 = args[2].parse().unwrap();
    std::fs::create_dir_all(pasta).unwrap();
    let txs = transcritoma_sintetico(&mut Rng(semente), genes);
    escreve_fasta(&format!("{pasta}/transcripts.fa"), &txs);
    let mut m = String::new();
    for t in &txs {
        m.push_str(&format!("{}\t{}\t{}\n", t.name, t.gene, t.symbol));
    }
    grava(&format!("{pasta}/gene_map.tsv"), &m);
    println!("{} transcripts -> {pasta}/transcripts.fa", txs.len());
}

fn modo_sim(args: &[String]) {
    let txs = le_fasta(&args[0]);
    let pasta = &args[1];
    let n: usize = args[2].parse().unwrap();
    let semente: u64 = args[3].parse().unwrap();
    let arg = |i: usize, d: f64| args.get(i).map(|x| x.parse::<f64>().unwrap()).unwrap_or(d);
    let paired = arg(4, 1.0) != 0.0;
    let o = SimOpt { n, paired, erro: arg(5, 0.005), read_len: arg(6, 100.0) as usize, media: arg(7, 200.0), dp: arg(8, 20.0) };
    let zeros = arg(9, 0.2);
    std::fs::create_dir_all(pasta).unwrap();
    let mut rng = Rng(semente);
    let theta = abundancias(&mut rng, txs.len(), zeros);
    let seqs: Vec<&[u8]> = txs.iter().map(|x| x.seq.as_slice()).collect();
    let inicio = Instant::now();
    let sim = simula(&seqs, &theta, &mut rng, &o);
    escreve_fastq(&format!("{pasta}/reads_1.fastq.gz"), &sim.r1, "/1");
    if paired {
        escreve_fastq(&format!("{pasta}/reads_2.fastq.gz"), &sim.r2, "/2");
    }
    let tot: f64 = theta.iter().zip(&seqs).filter(|(_, s)| s.len() >= 31).map(|(x, _)| x).sum();
    let mut tv = String::from("transcript\tgene\tlength\ttrue_tpm\ttrue_count\n");
    for (t, x) in txs.iter().enumerate() {
        let tp = if x.seq.len() >= 31 { theta[t] / tot * 1e6 } else { 0.0 };
        tv.push_str(&format!("{}\t{}\t{}\t{}\t{}\n", x.name, x.gene, x.seq.len(), tp, sim.contagem[t]));
    }
    grava(&format!("{pasta}/truth.tsv"), &tv);
    println!("{n} {} written to {pasta} in {:.1} s", if paired { "pairs" } else { "reads" }, inicio.elapsed().as_secs_f64());
}

fn modo_indice(args: &[String]) {
    let k = args.get(2).map(|x| x.parse().unwrap()).unwrap_or(31);
    let mapa = args.get(3).filter(|x| *x != "-").cloned();
    let memoria = args.get(4).map(|x| x.parse().unwrap());
    let r = build_transcript_index(&args[0], &args[1], &IndexOptions { k, gene_map: mapa, memory_mb: memoria, ..Default::default() });
    match r {
        Ok(r) => {
            print!("{}", r.resumo());
            println!("peak memory (working set): {:.0} MB", pico_mb());
        }
        Err(e) => {
            println!("ERROR: {e}");
            std::process::exit(1);
        }
    }
}

fn modo_quant(args: &[String]) {
    let (idx, pasta, nome, r1) = (&args[0], &args[1], &args[2], &args[3]);
    let r2 = args.get(4).filter(|x| *x != "-").cloned();
    let mut o = QuantOptions { keep_class_counts: true, ..Default::default() };
    if args.len() >= 7 {
        o.fragment_length = Some(args[5].parse().unwrap());
        o.fragment_sd = Some(args[6].parse().unwrap());
    }
    // an 8th argument "filter": kallisto's single-end rule (no --single-overhang)
    if args.get(7).map(String::as_str) == Some("filter") {
        o.single_overhang = false;
    }
    // o tempo de carregar o indice, medido a parte
    let tc = Instant::now();
    match TranscriptIndex::load(idx) {
        Ok(i) => println!("index loaded in {:.2} s ({} k-mers, {} classes, {:.0} MB)", tc.elapsed().as_secs_f64(), i.kmers(), i.classes(), i.bytes as f64 / 1048576.0),
        Err(e) => {
            println!("ERROR: {e}");
            std::process::exit(1);
        }
    }
    let t0 = Instant::now();
    let r = quantify_samples(idx, &[QuantSample { name: nome.clone(), r1: r1.clone(), r2, condition: None }], pasta, &o);
    match r {
        Ok(r) => {
            print!("{}", r.resumo());
            let s = &r.samples[0];
            let mut ec = String::from("transcripts\tcount\n");
            for (ts, c) in &s.class_counts {
                ec.push_str(&format!("{}\t{c}\n", ts.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")));
            }
            grava(&format!("{pasta}/{nome}/ec_counts.tsv"), &ec);
            println!(
                "time {:.2} s ({:.2} s per million reads), peak memory (working set) {:.0} MB",
                t0.elapsed().as_secs_f64(),
                s.seconds / s.reads as f64 * 1e6,
                pico_mb()
            );
        }
        Err(e) => {
            println!("ERROR: {e}");
            std::process::exit(1);
        }
    }
}

fn le_abundancia(p: &str) -> HashMap<String, (f64, f64)> {
    le_tsv(p).into_iter().skip(1).filter(|l| l.len() >= 5).map(|l| (l[0].clone(), (l[3].parse().unwrap(), l[4].parse().unwrap()))).collect()
}

fn modo_exatidao(args: &[String]) {
    let idx = TranscriptIndex::load(&args[0]).unwrap();
    let verdade: HashMap<String, (f64, f64)> =
        le_tsv(&args[1]).into_iter().skip(1).map(|l| (l[0].clone(), (l[4].parse().unwrap(), l[3].parse().unwrap()))).collect();
    let est = le_abundancia(&args[2]);
    let nt = idx.transcripts();
    let ng = idx.gene_ids.len();
    let gene_de: Vec<usize> = idx.gene_of.iter().map(|&g| g as usize).collect();
    let pega = |m: &HashMap<String, (f64, f64)>, i: usize| -> Vec<f64> {
        idx.names.iter().map(|n| m.get(n).map(|x| if i == 0 { x.0 } else { x.1 }).unwrap_or(0.0)).collect()
    };
    let (ec, et) = (pega(&est, 0), pega(&est, 1));
    let (vc, vt) = (pega(&verdade, 0), pega(&verdade, 1));
    println!("{nt} transcripts, {ng} genes");
    mostra("transcripts, est. counts vs truth", &metricas(&ec, &vc));
    mostra("transcripts, TPM vs truth", &metricas(&et, &vt));
    mostra("genes, est. counts vs truth", &metricas(&por_gene(&ec, &gene_de, ng), &por_gene(&vc, &gene_de, ng)));
    mostra("genes, TPM vs truth", &metricas(&por_gene(&et, &gene_de, ng), &por_gene(&vt, &gene_de, ng)));
    // genes com mais de uma isoforma
    let mut n_iso = vec![0; ng];
    for &g in &gene_de {
        n_iso[g] += 1;
    }
    let sel: Vec<usize> = (0..nt).filter(|&t| n_iso[gene_de[t]] > 1).collect();
    let f = |v: &[f64]| -> Vec<f64> { sel.iter().map(|&t| v[t]).collect() };
    mostra("transcripts of multi-isoform genes", &metricas(&f(&ec), &f(&vc)));
}

fn modo_compara(args: &[String]) {
    let idx = TranscriptIndex::load(&args[0]).unwrap();
    let a = le_abundancia(&args[1]);
    let b = le_abundancia(&args[2]);
    let ng = idx.gene_ids.len();
    let gene_de: Vec<usize> = idx.gene_of.iter().map(|&g| g as usize).collect();
    let pega = |m: &HashMap<String, (f64, f64)>, i: usize| -> Vec<f64> {
        idx.names.iter().map(|n| m.get(n).map(|x| if i == 0 { x.0 } else { x.1 }).unwrap_or(0.0)).collect()
    };
    let (ac, at, bc, bt) = (pega(&a, 0), pega(&a, 1), pega(&b, 0), pega(&b, 1));
    mostra("transcripts, est. counts", &metricas(&ac, &bc));
    println!("    Pearson (linear) = {:.6}", pearson(&ac, &bc));
    mostra("transcripts, TPM", &metricas(&at, &bt));
    let (ga, gb) = (por_gene(&ac, &gene_de, ng), por_gene(&bc, &gene_de, ng));
    mostra("genes, est. counts", &metricas(&ga, &gb));
    println!("    Pearson (linear) = {:.6}", pearson(&ga, &gb));
    let dif: Vec<f64> = ac.iter().zip(&bc).map(|(x, y)| (x - y).abs()).collect();
    let mut ord: Vec<usize> = (0..dif.len()).collect();
    ord.sort_by(|&x, &y| dif[y].total_cmp(&dif[x]));
    let maiores: f64 = ac.iter().sum::<f64>();
    println!("  total counts: ours {:.1}, other {:.1}", maiores, bc.iter().sum::<f64>());
    let (mut d1, mut d10) = (0, 0);
    for i in 0..ac.len() {
        if dif[i] > 1.0 {
            d1 += 1;
        }
        if dif[i] > 10.0 {
            d10 += 1;
        }
    }
    println!("  transcripts differing by > 1 read: {d1}; by > 10 reads: {d10}");
    for &i in ord.iter().take(5) {
        println!("    {}: ours {:.2}, other {:.2}", idx.names[i], ac[i], bc[i]);
    }
}

fn modo_compara_ec(args: &[String]) {
    let mut nossas: HashMap<Vec<u32>, u64> = HashMap::new();
    for l in le_tsv(&args[0]).into_iter().skip(1) {
        let ts: Vec<u32> = l[0].split(',').map(|x| x.parse().unwrap()).collect();
        nossas.insert(ts, l[1].parse().unwrap());
    }
    let mut ecs: HashMap<String, Vec<u32>> = HashMap::new();
    for l in le_tsv(&args[1]) {
        let mut ts: Vec<u32> = l[1].split(',').map(|x| x.parse().unwrap()).collect();
        ts.sort();
        ecs.insert(l[0].clone(), ts);
    }
    let mut delas: HashMap<Vec<u32>, u64> = HashMap::new();
    for l in le_tsv(&args[2]) {
        let c: u64 = l[1].parse().unwrap();
        if c > 0 {
            *delas.entry(ecs[&l[0]].clone()).or_insert(0) += c;
        }
    }
    let todas: BTreeSet<&Vec<u32>> = nossas.keys().chain(delas.keys()).collect();
    let (mut iguais, mut dif, mut so_n, mut so_k, mut lidos_dif) = (0, 0, 0, 0, 0u64);
    for k in &todas {
        match (nossas.get(*k), delas.get(*k)) {
            (Some(a), Some(b)) if a == b => iguais += 1,
            (Some(a), Some(b)) => {
                dif += 1;
                lidos_dif += a.abs_diff(*b);
            }
            (Some(a), None) => {
                so_n += 1;
                lidos_dif += a;
            }
            (None, Some(b)) => {
                so_k += 1;
                lidos_dif += b;
            }
            _ => {}
        }
    }
    let tn: u64 = nossas.values().sum();
    let tk: u64 = delas.values().sum();
    println!(
        "classes: {} in all; same count {iguais}, different count {dif}, only ours {so_n}, only kallisto {so_k}; reads in the differences {lidos_dif} (ours {tn} pseudoaligned, kallisto {tk})",
        todas.len()
    );
}

fn modo_dechain(args: &[String]) {
    let txs = le_fasta(&args[0]);
    let idx = &args[1];
    let pasta = &args[2];
    let n: usize = args[3].parse().unwrap();
    let semente: u64 = args[4].parse().unwrap();
    let plantados: usize = args.get(5).map(|x| x.parse().unwrap()).unwrap_or(15);
    let t0 = Instant::now();
    match cadeia_de(&txs, idx, pasta, n, semente, plantados, true) {
        Ok((s, d, f, testados, p)) => println!(
            "planted {p} up + {p} down (x4 / x0.25): found {s} up and {d} down at padj < 0.1; {f} false positives among {} other genes tested ({:.1} s)",
            testados - 2 * p,
            t0.elapsed().as_secs_f64()
        ),
        Err(e) => {
            println!("ERROR: {e}");
            std::process::exit(1);
        }
    }
}

/// The teaching demo of quant_demo.py: index, quantification of the six
/// samples, differential expression, and what was planted.
fn modo_demo(args: &[String]) {
    let pasta = &args[0];
    let p = |n: &str| format!("{pasta}/{n}");
    let r = build_transcript_index(&p("transcripts.fa"), &p("demo.idx"), &IndexOptions::default()).unwrap();
    print!("{}", r.resumo());
    let mut amostras = Vec::new();
    for l in le_tsv(&p("samples.tsv")).into_iter().skip(1) {
        amostras.push(QuantSample {
            name: l[0].clone(),
            r1: p(&format!("{}_R1.fastq.gz", l[0])),
            r2: Some(p(&format!("{}_R2.fastq.gz", l[0]))),
            condition: Some(l[1].clone()),
        });
    }
    let t0 = Instant::now();
    let q = quantify_samples(&p("demo.idx"), &amostras, &p("quant"), &QuantOptions::default()).unwrap();
    print!("{}", q.resumo());
    println!("quantified in {:.2} s", t0.elapsed().as_secs_f64());
    let de = expressao::differential_expression(
        &p("quant/gene_counts.tsv"),
        &p("quant/samples.tsv"),
        Some(&p("quant/de.tsv")),
        &expressao::DeOptions { condition: "condition".into(), reference: "control".into(), gene_names: Some(p("quant/gene_names.tsv")), ..Default::default() },
    )
    .unwrap();
    let verdade = le_tsv(&p("truth.tsv"));
    let mut plantado: HashMap<String, String> = HashMap::new();
    for l in verdade.iter().skip(1) {
        if !l[4].is_empty() {
            plantado.insert(l[1].clone(), l[4].clone());
        }
    }
    let mut achados = Vec::new();
    for g in &de.genes {
        if g.padj < 0.1 {
            achados.push(format!("{} ({}) log2FC {:+.2} padj {:.1e} {}", g.name, g.gene, g.log2_fold_change, g.padj, plantado.get(&g.gene).map(String::as_str).unwrap_or("NOT PLANTED")));
        }
    }
    println!("significant genes (padj < 0.1): {}", achados.len());
    for a in &achados {
        println!("  {a}");
    }
    for (g, k) in &plantado {
        if let Some(x) = de.genes.iter().find(|x| &x.gene == g) {
            println!("  planted {k} {}: log2FC {:+.2}, padj {:.2e}", x.name, x.log2_fold_change, x.padj);
        }
    }
    // a troca de isoformas
    let tt = le_tsv(&p("quant/transcript_tpm.tsv"));
    for l in &tt {
        if l[0].starts_with("DEMO0090") {
            println!("  isoform {}: TPM {}", l[0], l[2..].iter().map(|x| format!("{:.0}", x.parse::<f64>().unwrap_or(0.0))).collect::<Vec<_>>().join(" "));
        }
    }
    // exatidao por gene contra os fragmentos simulados
    let idx = TranscriptIndex::load(&p("demo.idx")).unwrap();
    let pos: HashMap<&str, usize> = idx.names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let gene_de: Vec<usize> = idx.gene_of.iter().map(|&g| g as usize).collect();
    for (si, s) in q.samples.iter().enumerate() {
        let mut v = vec![0.0; idx.transcripts()];
        for l in verdade.iter().skip(1) {
            v[pos[l[0].as_str()]] = l[6 + 2 * si].parse().unwrap();
        }
        let ng = idx.gene_ids.len();
        let mg = metricas(&por_gene(&s.est_counts, &gene_de, ng), &por_gene(&v, &gene_de, ng));
        let mt = metricas(&s.est_counts, &v);
        println!(
            "  {}: genes Pearson(log1p) {:.4} MARD {:.4}; transcripts Pearson(log1p) {:.4} Spearman {:.4} MARD {:.4}",
            s.sample, mg.pearson_log, mg.mard, mt.pearson_log, mt.spearman, mt.mard
        );
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("refseq") => return modo_refseq(&args[2..]),
        Some("synth") => return modo_sintetico(&args[2..]),
        Some("sim") => return modo_sim(&args[2..]),
        Some("index") => return modo_indice(&args[2..]),
        Some("quant") => return modo_quant(&args[2..]),
        Some("accuracy") => return modo_exatidao(&args[2..]),
        Some("compare") => return modo_compara(&args[2..]),
        Some("compare-ec") => return modo_compara_ec(&args[2..]),
        Some("dechain") => return modo_dechain(&args[2..]),
        Some("demo") => return modo_demo(&args[2..]),
        Some("explain") => {
            match explain_pseudoalignment(&args[2], &args[3], args.get(4).map(String::as_str)) {
                Ok(j) => println!("{j}"),
                Err(e) => println!("ERROR: {e}"),
            }
            return;
        }
        _ => {}
    }
    let dir = std::env::temp_dir().join(format!("mmgen_quantificacao_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut t = T { ok: 0, falhas: 0, dir: dir.clone() };
    let inicio = Instant::now();

    testa_kmers(&mut t);
    testa_cabecalhos(&mut t);
    testa_mapas(&mut t);
    let txs = transcritoma_sintetico(&mut Rng(1), 300);
    let pfa = t.p("synt.fa");
    escreve_fasta(&pfa, &txs);
    testa_indice(&mut t, &txs, &pfa);
    let idx_path = t.p("synt_k31.idx");
    let idx = TranscriptIndex::load(&idx_path).unwrap();
    testa_pseudoalinhamento(&mut t, &txs, &idx);
    drop(idx);
    testa_em(&mut t);
    testa_comprimentos(&mut t);
    testa_quantificacao(&mut t, &txs, &idx_path);
    testa_cadeia(&mut t, &txs, &idx_path);
    testa_erros(&mut t, &txs, &idx_path);
    testa_explicacao(&mut t, &txs, &idx_path);

    let _ = std::fs::remove_dir_all(&dir);
    println!("quantificacao tests: {} passed, {} failed ({:.1} s)", t.ok, t.falhas, inicio.elapsed().as_secs_f64());
    if t.falhas > 0 {
        std::process::exit(1);
    }
}

