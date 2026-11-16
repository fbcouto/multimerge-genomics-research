// Tests of src/navegador.rs: the data behind the genome browser view
// (region_reads), the depth profile for copy-number plots (depth_profile)
// and the reference bases of a region (fetch_sequence).
//
// Every result is compared with something worked out independently:
//  - coverage: `clinico::pileup` (validated against `samtools mpileup`) at
//    EVERY position of every window, with three sets of filters: depth, A, C,
//    G, T, deletions, insertions after, and the forward-strand counts;
//  - the drawn reads: a generator that records, while it builds each read,
//    its blocks, mismatches, insertions, deletions, introns (N), clips and
//    mate; mismatches checked against a FASTA and against the MD tag (written
//    as `samtools calmd` writes it, sometimes without the `0` after a
//    deletion, sometimes missing or broken);
//  - indexed BAM, BAM without index, SAM text and SAM.gz give the same view;
//  - row packing: no two reads of a row overlap (1 free base between them),
//    and the number of rows is the minimum (the largest set of reads that
//    pairwise touch, found by a sweep);
//  - depth_profile: per-position sums of the aligned bases of the generator;
//  - fetch_sequence: slices of the generated genome, through plain FASTA,
//    FASTA + .fai, BGZF FASTA and FASTA with irregular lines.
//
// Usage: cargo run --release --bin navegador
// Optional environment, to add checks on real data:
//   NAV_EXAMPLES=<folder>  pileup_demo.bam(+.bai), toy_genome.fa,
//                          alignments.sorted.bam, rnaseq_ctrl_1.bam
//                          (the files the web interface writes in mmgen_examples/)
//   NAV_REAL_BAM=<bam>     a real BAM with its .bai (e.g. the HG002 slice)
#[path = "../../../src/multimerge_core.rs"]
#[allow(dead_code)]
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
#[path = "../../../src/indice.rs"]
#[allow(dead_code)]
mod indice;
#[path = "../../../src/referencia.rs"]
#[allow(dead_code)]
mod referencia;
#[path = "../../../src/clinico.rs"]
#[allow(dead_code)]
mod clinico;
#[path = "../../../src/navegador.rs"]
mod navegador;

use navegador::{depth_profile, fetch_sequence, region_reads, ReadsOptions};
use serde_json::{json, Value};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::Debug;
use std::io::{BufRead, Write};
use std::path::PathBuf;

// ---------------------------------------------------------------- harness

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
    fn chance(&mut self, pct: u64) -> bool {
        self.below(100) < pct
    }
    fn pick<'a, X>(&mut self, v: &'a [X]) -> &'a X {
        &v[self.below(v.len() as u64) as usize]
    }
}

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
    /// Two lists, reporting only the first difference.
    fn lista(&mut self, got: &[String], want: &[String], what: &str) -> bool {
        for (i, (g, w)) in got.iter().zip(want).enumerate() {
            if g != w {
                println!("FAIL  {what}: first difference at item {i}\n   got:  {g}\n   want: {w}");
                return self.check(false, what);
            }
        }
        if got.len() != want.len() {
            println!("FAIL  {what}: {} items, expected {}", got.len(), want.len());
            return self.check(false, what);
        }
        self.check(true, what)
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
    /// Runs and parses a JSON result.
    fn json(&mut self, r: Result<String, String>, what: &str) -> Option<Value> {
        let s = self.ok_ou(r, what)?;
        match serde_json::from_str::<Value>(&s) {
            Ok(v) => Some(v),
            Err(e) => {
                println!("FAIL  {what}: invalid JSON: {e}\n{}", &s[..s.len().min(400)]);
                self.check(false, what);
                None
            }
        }
    }
}

fn grava(p: &str, linhas: &[String]) {
    let mut w = externo::Escritor::cria(p).unwrap();
    for l in linhas {
        w.write_all(l.as_bytes()).unwrap();
        w.write_all(b"\n").unwrap();
    }
    w.termina().unwrap();
}

fn opcoes(reference: Option<&str>, max_reads: usize) -> ReadsOptions {
    ReadsOptions { reference: reference.map(str::to_string), max_reads, ..Default::default() }
}

fn u(v: &Value) -> u64 {
    v.as_u64().unwrap_or(u64::MAX)
}

/// The view without the keys that depend on how the file was read.
fn sem_modo(v: &Value) -> Value {
    let mut v = v.clone();
    if let Some(o) = v.as_object_mut() {
        for k in ["path", "indexed", "warnings", "text"] {
            o.remove(k);
        }
    }
    v
}

fn avisos(v: &Value) -> Vec<String> {
    v["warnings"].as_array().map(|a| a.iter().map(|x| x.as_str().unwrap_or("").to_string()).collect()).unwrap_or_default()
}

// ---------------------------------------------------------------- generator

const BASES: &[u8] = b"ACGT";

fn outra(rng: &mut Rng, b: u8) -> u8 {
    loop {
        let x = BASES[rng.below(4) as usize];
        if x != b {
            return x;
        }
    }
}

#[derive(Clone, Default)]
struct Sim {
    nome: String,
    flag: u16,
    crom: usize,
    pos: u64,
    mapq: u8,
    ops: Vec<(u32, u8)>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    seq_estrela: bool,
    qual_estrela: bool,
    /// MD tag as written (None: no tag) and whether it fits the CIGAR.
    md: Option<String>,
    md_valido: bool,
    rnext: String,
    pnext: u64,
    tlen: i64,
    // Truth, recorded while the read is built.
    fim: u64,
    alinhadas: Vec<(u64, usize)>,
    deletadas: Vec<u64>,
    /// Every aligned base that differs from the reference: (pos, base, qual).
    mism: Vec<(u64, u8, u8)>,
    ins: Vec<(u64, u32, String)>,
    dels: Vec<(u64, u32)>,
    skips: Vec<(u64, u32)>,
    blocks: Vec<(u64, u64)>,
    clips: [u32; 4],
    clip_seqs: (String, String),
}

impl Sim {
    fn cigar(&self) -> String {
        self.ops.iter().map(|(n, o)| format!("{n}{}", *o as char)).collect()
    }
    fn sam(&self, nomes: &[&str]) -> String {
        let seq = if self.seq_estrela { "*".to_string() } else { String::from_utf8(self.seq.clone()).unwrap() };
        let qual = if self.seq_estrela || self.qual_estrela { "*".to_string() } else { String::from_utf8(self.qual.clone()).unwrap() };
        let md = self.md.as_ref().map(|m| format!("\tMD:Z:{m}")).unwrap_or_default();
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNM:i:0{}",
            self.nome,
            self.flag,
            nomes[self.crom],
            self.pos,
            self.mapq,
            self.cigar(),
            self.rnext,
            self.pnext,
            self.tlen,
            seq,
            qual,
            md
        )
    }
    fn tem_seq(&self) -> bool {
        !self.seq_estrela
    }
}

/// Builds a read at `inicio` over `refseq`, recording what it shows.
fn gera_read(rng: &mut Rng, nome: String, crom: usize, refseq: &[u8], inicio: u64) -> Option<Sim> {
    let mut s = Sim { nome, crom, pos: inicio, ..Default::default() };
    let mut r = inicio;
    let mut md = String::new();
    let mut um: u32 = 0;
    let qual = |rng: &mut Rng| 33 + 2 + rng.below(40) as u8;
    let aleat = |rng: &mut Rng, n: u32, s: &mut Sim| -> String {
        let mut t = String::new();
        for _ in 0..n {
            let b = BASES[rng.below(4) as usize];
            let b = if rng.chance(5) { b.to_ascii_lowercase() } else { b };
            s.seq.push(b);
            s.qual.push(qual(rng));
            t.push(b.to_ascii_uppercase() as char);
        }
        t
    };
    if rng.chance(10) {
        let n = 1 + rng.below(5) as u32;
        s.ops.push((n, b'H'));
        s.clips[2] = n;
    }
    if rng.chance(30) {
        let n = 1 + rng.below(10) as u32;
        s.ops.push((n, b'S'));
        s.clip_seqs.0 = aleat(rng, n, &mut s);
        s.clips[0] = n;
    }
    if rng.chance(5) {
        let n = 1 + rng.below(3) as u32;
        s.ops.push((n, b'I'));
        let t = aleat(rng, n, &mut s);
        s.ins.push((inicio - 1, n, t));
    }
    let alvo = 30 + rng.below(170) as u32;
    let (mut usado, mut introns) = (0u32, 0);
    loop {
        let tipo = match rng.below(10) {
            0 => b'=',
            1 => b'X',
            _ => b'M',
        };
        let n = 1 + rng.below(40) as u32;
        if r + u64::from(n) > refseq.len() as u64 {
            return None;
        }
        s.ops.push((n, tipo));
        match s.blocks.last_mut() {
            Some(b) if b.1 + 1 == r => b.1 = r + u64::from(n) - 1,
            _ => s.blocks.push((r, r + u64::from(n) - 1)),
        }
        for _ in 0..n {
            let rb = refseq[r as usize - 1];
            let b = match tipo {
                b'=' => rb,
                b'X' => outra(rng, rb),
                _ => {
                    if rng.chance(1) {
                        b'N'
                    } else if rng.chance(8) {
                        outra(rng, rb)
                    } else {
                        rb
                    }
                }
            };
            let q = qual(rng);
            let qi = s.seq.len();
            s.seq.push(if rng.chance(3) { b.to_ascii_lowercase() } else { b });
            s.qual.push(q);
            s.alinhadas.push((r, qi));
            if b != rb {
                s.mism.push((r, b, q - 33));
                md.push_str(&um.to_string());
                md.push(rb as char);
                um = 0;
            } else {
                um += 1;
            }
            r += 1;
        }
        usado += n;
        if usado >= alvo {
            break;
        }
        let deleta = |rng: &mut Rng, s: &mut Sim, md: &mut String, um: &mut u32, r: &mut u64| -> Option<u32> {
            let n = 1 + rng.below(5) as u32;
            if *r + u64::from(n) + 50 > refseq.len() as u64 {
                return None;
            }
            s.ops.push((n, b'D'));
            md.push_str(&um.to_string());
            md.push('^');
            for k in 0..n {
                md.push(refseq[(*r + u64::from(k)) as usize - 1] as char);
                s.deletadas.push(*r + u64::from(k));
            }
            *um = 0;
            *r += u64::from(n);
            Some(n)
        };
        match rng.below(100) {
            0..=24 => {
                let n = 1 + rng.below(4) as u32;
                s.ops.push((n, b'I'));
                let t = aleat(rng, n, &mut s);
                s.ins.push((r - 1, n, t));
            }
            25..=49 => {
                let p = r;
                let n = deleta(rng, &mut s, &mut md, &mut um, &mut r)?;
                s.dels.push((p, n));
            }
            50..=59 if introns < 2 => {
                let n = 20 + rng.below(280) as u32;
                if r + u64::from(n) + 50 > refseq.len() as u64 {
                    return None;
                }
                s.ops.push((n, b'N'));
                s.skips.push((r, n));
                r += u64::from(n);
                introns += 1;
            }
            60..=64 => {
                // Two D operations in a row: one deletion on screen, two `^` groups in MD.
                let p = r;
                let a = deleta(rng, &mut s, &mut md, &mut um, &mut r)?;
                let b = deleta(rng, &mut s, &mut md, &mut um, &mut r)?;
                s.dels.push((p, a + b));
            }
            65..=69 => {
                let p = r;
                let n = deleta(rng, &mut s, &mut md, &mut um, &mut r)?;
                s.dels.push((p, n));
                let k = 1 + rng.below(3) as u32;
                s.ops.push((k, b'I'));
                let t = aleat(rng, k, &mut s);
                s.ins.push((r - 1, k, t));
            }
            70..=74 => {
                s.ops.push((1, b'P'));
                let k = 1 + rng.below(3) as u32;
                s.ops.push((k, b'I'));
                let t = aleat(rng, k, &mut s);
                s.ins.push((r - 1, k, t));
            }
            _ => {}
        }
    }
    md.push_str(&um.to_string());
    if rng.chance(25) {
        let n = 1 + rng.below(8) as u32;
        s.ops.push((n, b'S'));
        s.clip_seqs.1 = aleat(rng, n, &mut s);
        s.clips[1] = n;
    }
    if rng.chance(8) {
        let n = 1 + rng.below(4) as u32;
        s.ops.push((n, b'H'));
        s.clips[3] = n;
    }
    s.fim = r - 1;
    // The MD tag: as calmd writes it; sometimes without the `0` between a
    // deletion and a mismatch right after it; sometimes missing or broken.
    match rng.below(100) {
        0..=9 => s.md = None,
        10..=12 => {
            s.md = Some(format!("{md}1"));
            s.md_valido = false;
        }
        13..=32 => {
            let b = md.as_bytes();
            let mut o = String::new();
            let mut em_del = false;
            for (i, &c) in b.iter().enumerate() {
                if c == b'^' {
                    em_del = true;
                } else if c.is_ascii_digit() {
                    if em_del && c == b'0' && b.get(i + 1).is_some_and(|x| x.is_ascii_alphabetic()) && (i == 0 || b[i - 1].is_ascii_alphabetic()) {
                        continue;
                    }
                    em_del = false;
                }
                o.push(c as char);
            }
            s.md = Some(o);
            s.md_valido = true;
        }
        _ => {
            s.md = Some(md);
            s.md_valido = true;
        }
    }
    Some(s)
}

const NOMES: [&str; 2] = ["chr7", "chr8"];
const TAMANHOS: [usize; 2] = [3000, 2000];

struct Dados {
    genoma: Vec<Vec<u8>>,
    reads: Vec<Sim>,
    bam: String,
    bam_sem: String,
    sam: String,
    sam_gz: String,
    fasta: String,
}

fn gera_dados(t: &T, rng: &mut Rng) -> Dados {
    let genoma: Vec<Vec<u8>> = TAMANHOS.iter().map(|&l| (0..l).map(|_| BASES[rng.below(4) as usize]).collect()).collect();
    let mut reads: Vec<Sim> = Vec::new();
    let mut i = 0;
    while reads.len() < 2400 {
        i += 1;
        let c = if i % 4 == 0 { 1 } else { 0 };
        let inicio = 1 + rng.below(genoma[c].len() as u64 - 100);
        let Some(mut s) = gera_read(rng, format!("r{i}"), c, &genoma[c], inicio) else { continue };
        s.flag = match rng.below(100) {
            0..=4 => 1024,
            5..=7 => 256,
            8..=10 => 2048,
            11..=12 => 512,
            _ => *rng.pick(&[0u16, 16, 99, 147, 83, 163, 97, 73]),
        };
        s.mapq = if rng.chance(5) { 255 } else { rng.below(61) as u8 };
        if rng.chance(5) {
            s.qual_estrela = true;
        } else if rng.chance(1) {
            s.seq_estrela = true;
        }
        if s.flag & 1 != 0 {
            if s.flag & 8 != 0 {
                s.rnext = "=".into();
                s.pnext = s.pos;
                s.tlen = 0;
            } else {
                s.rnext = if rng.chance(85) { "=".into() } else { NOMES[1 - c].to_string() };
                s.pnext = 1 + rng.below(TAMANHOS[1 - c] as u64 - 10);
                s.tlen = (100 + rng.below(500) as i64) * if rng.chance(50) { 1 } else { -1 };
            }
        } else {
            s.rnext = "*".into();
        }
        // A mate with the same name, nearby.
        let par = s.flag == 99 && rng.chance(70);
        let pos = s.pos;
        reads.push(s);
        if par {
            let inicio = (pos + rng.below(250)).min(genoma[c].len() as u64 - 100);
            if let Some(mut m) = gera_read(rng, format!("r{i}"), c, &genoma[c], inicio) {
                m.flag = 147;
                m.mapq = 60;
                m.rnext = "=".into();
                m.pnext = pos;
                m.tlen = -(m.fim as i64 - pos as i64 + 1);
                reads.push(m);
            }
        }
    }
    reads.sort_by_key(|r| (r.crom, r.pos));
    let mut linhas = vec!["@HD\tVN:1.6\tSO:coordinate".to_string()];
    for (n, l) in NOMES.iter().zip(TAMANHOS) {
        linhas.push(format!("@SQ\tSN:{n}\tLN:{l}"));
    }
    linhas.extend(reads.iter().map(|r| r.sam(&NOMES)));
    let d = Dados {
        bam: t.p("leituras.bam"),
        bam_sem: t.p("sem_indice.bam"),
        sam: t.p("leituras.sam"),
        sam_gz: t.p("leituras.sam.gz"),
        fasta: t.p("genoma.fa"),
        genoma,
        reads,
    };
    grava(&d.bam, &linhas);
    grava(&d.bam_sem, &linhas);
    grava(&d.sam, &linhas);
    grava(&d.sam_gz, &linhas);
    indice::indexar(&d.bam, None).unwrap();
    escreve_fasta(&d.fasta, &d.genoma, &NOMES, 60);
    d
}

fn escreve_fasta(p: &str, genoma: &[Vec<u8>], nomes: &[&str], largura: usize) {
    let mut s = String::new();
    for (g, n) in genoma.iter().zip(nomes) {
        s.push_str(&format!(">{n} test sequence\n"));
        for c in g.chunks(largura) {
            s.push_str(std::str::from_utf8(c).unwrap());
            s.push('\n');
        }
    }
    std::fs::write(p, s).unwrap();
}

// ---------------------------------------------------------------- expected reads

fn passa(s: &Sim, o: &ReadsOptions) -> bool {
    s.flag & 4 == 0 && s.flag & o.exclude_flags == 0 && s.mapq >= o.min_mapq
}

fn esperado(s: &Sim, ws: u64, we: u64, com_ref: bool) -> Value {
    let mism: Vec<Value> = if !s.tem_seq() || (!com_ref && !s.md_valido) {
        Vec::new()
    } else {
        s.mism
            .iter()
            .filter(|m| m.0 >= ws && m.0 <= we)
            .map(|&(p, b, q)| json!([p, (b as char).to_string(), if s.qual_estrela { Value::Null } else { json!(q) }]))
            .collect()
    };
    let ins: Vec<Value> = s
        .ins
        .iter()
        .filter(|x| x.0 + 1 >= ws && x.0 <= we)
        .map(|(p, n, t)| json!([p, n, if s.tem_seq() { t.clone() } else { "N".repeat(*n as usize) }]))
        .collect();
    let toca = |p: u64, n: u32| p <= we && p + u64::from(n) > ws;
    let dels: Vec<Value> = s.dels.iter().filter(|x| toca(x.0, x.1)).map(|(p, n)| json!([p, n])).collect();
    let skips: Vec<Value> = s.skips.iter().filter(|x| toca(x.0, x.1)).map(|(p, n)| json!([p, n])).collect();
    let mate = if s.flag & 1 != 0 {
        json!({
            "chrom": if s.rnext == "=" { NOMES[s.crom].to_string() } else { s.rnext.clone() },
            "start": if s.pnext > 0 { json!(s.pnext) } else { Value::Null },
            "strand": if s.flag & 0x20 != 0 { "-" } else { "+" },
            "unmapped": s.flag & 8 != 0,
            "proper_pair": s.flag & 2 != 0,
            "insert_size": s.tlen,
            "read_number": if s.flag & 0x40 != 0 { json!(1) } else if s.flag & 0x80 != 0 { json!(2) } else { Value::Null },
        })
    } else {
        Value::Null
    };
    let fonte = if com_ref {
        "reference"
    } else if s.md_valido {
        "MD"
    } else {
        "none"
    };
    json!({
        "name": s.nome,
        "flag": s.flag,
        "strand": if s.flag & 16 != 0 { "-" } else { "+" },
        "mapq": s.mapq,
        "start": s.pos,
        "end": s.fim,
        "cigar": s.cigar(),
        "clip_left": s.clips[0],
        "clip_right": s.clips[1],
        "hard_clip_left": s.clips[2],
        "hard_clip_right": s.clips[3],
        "clip_left_seq": if s.tem_seq() { s.clip_seqs.0.clone() } else { String::new() },
        "clip_right_seq": if s.tem_seq() { s.clip_seqs.1.clone() } else { String::new() },
        "blocks": s.blocks.iter().map(|(a, b)| json!([a, b])).collect::<Vec<_>>(),
        "mismatches": mism,
        "insertions": ins,
        "deletions": dels,
        "skips": skips,
        "mate": mate,
        "has_sequence": s.tem_seq(),
        "mismatch_source": fonte,
    })
}

fn media_q(s: &Sim) -> Option<f64> {
    if !s.tem_seq() || s.qual_estrela {
        return None;
    }
    Some(s.qual.iter().map(|&q| f64::from(q - 33)).sum::<f64>() / s.qual.len() as f64)
}

/// Compares every drawn read with the generator's record. Returns the
/// number of reads compared.
fn confere_leituras(t: &mut T, d: &Dados, v: &Value, crom: usize, ws: u64, we: u64, o: &ReadsOptions, com_ref: bool, rotulo: &str) -> usize {
    let mut esp: BTreeMap<(String, u64), (Value, Option<f64>)> = BTreeMap::new();
    for s in d.reads.iter().filter(|s| s.crom == crom && s.pos <= we && s.fim >= ws && passa(s, o)) {
        esp.insert((s.nome.clone(), u64::from(s.flag)), (esperado(s, ws, we, com_ref), media_q(s)));
    }
    let lidas = v["reads"].as_array().cloned().unwrap_or_default();
    let mut ok = t.eq(&lidas.len(), &esp.len(), &format!("{rotulo}: number of reads drawn"));
    let mut n = 0;
    for r in &lidas {
        let chave = (r["name"].as_str().unwrap_or("").to_string(), u(&r["flag"]));
        let Some((e, mq)) = esp.get(&chave) else {
            println!("FAIL  {rotulo}: unexpected read {chave:?}");
            ok = false;
            continue;
        };
        for (k, ev) in e.as_object().unwrap() {
            if &r[k] != ev {
                println!("FAIL  {rotulo}: read {} (flag {}), key `{k}`\n   got:  {}\n   want: {}", chave.0, chave.1, r[k], ev);
                ok = false;
            }
        }
        let mq_ok = match (mq, r["mean_qual"].as_f64()) {
            (None, None) => r["mean_qual"].is_null(),
            (Some(a), Some(b)) => (a - b).abs() <= 0.005 + 1e-9,
            _ => false,
        };
        if !mq_ok {
            println!("FAIL  {rotulo}: read {}, mean_qual {} expected {mq:?}", chave.0, r["mean_qual"]);
            ok = false;
        }
        n += 1;
    }
    t.check(ok, &format!("{rotulo}: {n} reads equal to the generator (blocks, mismatches, insertions, deletions, skips, clips, mate)"));
    n
}

// ---------------------------------------------------------------- packing

fn confere_linhas(t: &mut T, v: &Value, rotulo: &str) {
    let lidas = v["reads"].as_array().cloned().unwrap_or_default();
    let linhas = u(&v["rows"]) as usize;
    let mut por_linha: Vec<Vec<(u64, u64)>> = vec![Vec::new(); linhas];
    let mut ok = true;
    for r in &lidas {
        let l = u(&r["row"]) as usize;
        if l >= linhas {
            ok = false;
            continue;
        }
        por_linha[l].push((u(&r["start"]), u(&r["end"])));
    }
    for (i, pl) in por_linha.iter_mut().enumerate() {
        pl.sort();
        if pl.is_empty() {
            println!("FAIL  {rotulo}: row {i} is empty");
            ok = false;
        }
        for w in pl.windows(2) {
            if w[0].1 + 1 >= w[1].0 {
                println!("FAIL  {rotulo}: row {i}: {:?} and {:?} touch", w[0], w[1]);
                ok = false;
            }
        }
    }
    // The minimum: the most reads (each with its 1 free base after it) over one point.
    let mut ev: Vec<(u64, i32)> = Vec::new();
    for r in &lidas {
        ev.push((u(&r["start"]), 1));
        ev.push((u(&r["end"]) + 2, -1));
    }
    ev.sort();
    let (mut cur, mut max) = (0i32, 0i32);
    for (_, d) in ev {
        cur += d;
        max = max.max(cur);
    }
    t.check(ok, &format!("{rotulo}: no two reads of a row overlap ({} reads, {linhas} rows)", lidas.len()));
    t.eq(&linhas, &(max as usize), &format!("{rotulo}: the number of rows is the minimum"));
}

// ---------------------------------------------------------------- coverage vs pileup

fn confere_pileup(t: &mut T, arquivo: &str, chrom: &str, v: &Value, o: &ReadsOptions, rotulo: &str) -> usize {
    let po = clinico::PileupOptions { min_mapq: o.min_mapq, min_baseq: o.min_baseq, exclude_flags: o.exclude_flags, alt: None };
    let cob = v["coverage"].as_array().cloned().unwrap_or_default();
    let fr = v["coverage_forward"].as_array().cloned().unwrap_or_default();
    t.eq(&cob.len(), &fr.len(), &format!("{rotulo}: coverage and coverage_forward have the same length"));
    let (mut got, mut esp) = (Vec::new(), Vec::new());
    for (c, f) in cob.iter().zip(&fr) {
        let p = u(&c[0]);
        let n: Vec<u64> = c.as_array().unwrap().iter().map(u).collect();
        let m: Vec<u64> = f.as_array().unwrap().iter().map(u).collect();
        got.push(format!("{p} depth={} A={} C={} G={} T={} del={} ins={} | fwd {:?}", n[1], n[2], n[3], n[4], n[5], n[6], n[7], &m[1..]));
        match clinico::pileup(arquivo, chrom, p, &po) {
            Ok(x) => {
                let fw = |a: &str| x.alleles.iter().find(|y| y.allele == a).map(|y| y.forward as u64).unwrap_or(0);
                let ins_f: u64 = x.alleles.iter().filter(|y| y.allele.starts_with('+')).map(|y| y.forward as u64).sum();
                let fwd = vec![fw("A"), fw("C"), fw("G"), fw("T"), fw("*"), ins_f];
                esp.push(format!(
                    "{p} depth={} A={} C={} G={} T={} del={} ins={} | fwd {:?}",
                    x.depth, x.a, x.c, x.g, x.t, x.deletions, x.insertions_after, fwd
                ));
            }
            Err(e) => esp.push(format!("error {e}")),
        }
    }
    t.lista(&got, &esp, &format!("{rotulo}: {} positions equal to pileup", esp.len()));
    esp.len()
}

// ---------------------------------------------------------------- 1. synthetic tests

fn testa_cobertura_e_leituras(t: &mut T, d: &Dados) {
    let conjuntos = [
        ReadsOptions { max_reads: 100_000, ..Default::default() },
        ReadsOptions { max_reads: 100_000, min_mapq: 0, min_baseq: 0, exclude_flags: 0, ..Default::default() },
        ReadsOptions { max_reads: 100_000, min_mapq: 20, min_baseq: 25, exclude_flags: 0x400, ..Default::default() },
    ];
    let mut posicoes = 0;
    let mut leituras = 0;
    let mut tipos = [0u64; 4]; // positions with ins, del, N skip reads, mismatches drawn
    for (oi, o) in conjuntos.iter().enumerate() {
        for (c, nome) in NOMES.iter().enumerate() {
            let mut ini = 1u64;
            while ini <= TAMANHOS[c] as u64 {
                let fim = (ini + 599).min(TAMANHOS[c] as u64);
                let reg = format!("{nome}:{ini}-{fim}");
                let rot = format!("view {reg}, options {oi}");
                if let Some(v) = t.json(region_reads(&d.bam, &reg, o), &rot) {
                    t.check(v["indexed"] == json!(true), &format!("{rot}: read through the index"));
                    posicoes += confere_pileup(t, &d.bam, nome, &v, o, &rot);
                    leituras += confere_leituras(t, d, &v, c, ini, fim, o, false, &rot);
                    confere_linhas(t, &v, &rot);
                    if oi == 1 {
                        for x in v["coverage"].as_array().unwrap() {
                            tipos[0] += u64::from(u(&x[7]) > 0);
                            tipos[1] += u64::from(u(&x[6]) > 0);
                        }
                        for r in v["reads"].as_array().unwrap() {
                            tipos[2] += u64::from(!r["skips"].as_array().unwrap().is_empty());
                            tipos[3] += r["mismatches"].as_array().unwrap().len() as u64;
                        }
                    }
                }
                ini = fim + 1;
            }
        }
    }
    println!("  coverage: {posicoes} positions compared with pileup; {leituras} reads compared with the generator");
    t.check(posicoes == 3 * 5000, "coverage: every position of both sequences, with 3 sets of filters");
    t.check(tipos.iter().all(|&n| n > 50), &format!("coverage: positions with insertions after, deletions; reads with N; mismatches: {tipos:?}"));

    // With the FASTA: mismatches against it, the window's bases; same coverage.
    let o = ReadsOptions { max_reads: 100_000, reference: Some(d.fasta.clone()), ..Default::default() };
    for (c, nome) in NOMES.iter().enumerate() {
        for (ini, fim) in [(1u64, 700u64), (650, 1400), (1300, TAMANHOS[c] as u64)] {
            let reg = format!("{nome}:{ini}-{fim}");
            let rot = format!("view {reg} with the FASTA");
            if let Some(v) = t.json(region_reads(&d.bam, &reg, &o), &rot) {
                confere_leituras(t, d, &v, c, ini, fim, &o, true, &rot);
                let janela = std::str::from_utf8(&d.genoma[c][ini as usize - 1..fim as usize]).unwrap();
                t.eq(&v["reference"].as_str().unwrap_or(""), &janela, &format!("{rot}: the window's bases"));
                t.eq(&v["reference_source"], &json!("FASTA"), &format!("{rot}: reference source"));
                t.eq(&v["mismatch_source"], &json!("reference"), &format!("{rot}: mismatch source"));
                t.check(!avisos(&v).iter().any(|a| a.contains("disagree")), &format!("{rot}: MD and FASTA agree everywhere"));
                let sem_ref = t.json(region_reads(&d.bam, &reg, &opcoes(None, 100_000)), &rot);
                if let Some(s) = sem_ref {
                    t.eq(&v["coverage"], &s["coverage"], &format!("{rot}: the coverage does not depend on the FASTA"));
                }
            }
        }
    }

    // Without the FASTA: the window's bases as the MD tags tell them.
    let reg = "chr7:1-3000";
    if let Some(v) = t.json(region_reads(&d.bam, reg, &opcoes(None, 100_000)), "MD reference") {
        let o = opcoes(None, 100_000);
        let mut r = vec![b'N'; 3000];
        for s in d.reads.iter().filter(|s| s.crom == 0 && passa(s, &o) && s.md_valido) {
            let fora: HashSet<u64> = s.mism.iter().map(|m| m.0).collect();
            for &(p, qi) in &s.alinhadas {
                let b = s.seq[qi].to_ascii_uppercase();
                if (s.tem_seq() && matches!(b, b'A' | b'C' | b'G' | b'T')) || fora.contains(&p) {
                    r[p as usize - 1] = d.genoma[0][p as usize - 1];
                }
            }
            for &p in &s.deletadas {
                r[p as usize - 1] = d.genoma[0][p as usize - 1];
            }
        }
        t.eq(&v["reference"].as_str().unwrap_or(""), &std::str::from_utf8(&r).unwrap(), "MD reference: the window's bases rebuilt from the MD tags");
        t.eq(&v["reference_source"], &json!("MD"), "MD reference: source");
        t.eq(&v["mismatch_source"], &json!("MD"), "MD reference: mismatch source");
        let sem_md = d.reads.iter().filter(|s| s.crom == 0 && passa(s, &o) && s.md.is_none()).count();
        let quebrado = d.reads.iter().filter(|s| s.crom == 0 && passa(s, &o) && s.md.is_some() && !s.md_valido).count();
        let w = avisos(&v);
        t.check(w.iter().any(|a| a.starts_with(&format!("{sem_md} of the "))), &format!("MD reference: warning about the {sem_md} reads without MD"));
        t.check(w.iter().any(|a| a.starts_with(&format!("{quebrado} read(s) drawn have an MD tag that does not fit"))), &format!("MD reference: warning about the {quebrado} broken MD tags"));
    }

    // A FASTA that is not the one the reads were aligned to.
    let mut errado = d.genoma.clone();
    let mut trocadas: HashSet<u64> = HashSet::new();
    let mut rng = Rng(77);
    while trocadas.len() < 30 {
        let p = 1 + rng.below(3000);
        if trocadas.insert(p) {
            let b = errado[0][p as usize - 1];
            errado[0][p as usize - 1] = outra(&mut rng, b);
        }
    }
    let fa_errado = t.p("errado.fa");
    escreve_fasta(&fa_errado, &errado, &NOMES, 70);
    let o = ReadsOptions { max_reads: 100_000, reference: Some(fa_errado.clone()), ..Default::default() };
    if let Some(v) = t.json(region_reads(&d.bam, "chr7:1-3000", &o), "wrong FASTA") {
        let mut n_pos = 0;
        let mut n_reads = 0;
        for s in d.reads.iter().filter(|s| s.crom == 0 && passa(s, &o) && s.md_valido && s.tem_seq()) {
            let k = s.alinhadas.iter().filter(|a| trocadas.contains(&a.0)).count() + s.deletadas.iter().filter(|p| trocadas.contains(p)).count();
            n_pos += k;
            n_reads += usize::from(k > 0);
        }
        let w = avisos(&v);
        t.check(
            w.iter().any(|a| a.starts_with(&format!("At {n_pos} base(s) of {n_reads} read(s) drawn, the MD tag and the reference FASTA disagree"))),
            &format!("wrong FASTA: warns that MD and FASTA disagree at {n_pos} bases of {n_reads} reads"),
        );
        // Mismatches now follow the (wrong) FASTA.
        let mut ok = true;
        for r in v["reads"].as_array().unwrap() {
            for m in r["mismatches"].as_array().unwrap() {
                let p = u(&m[0]);
                if m[1].as_str().unwrap().as_bytes()[0] == errado[0][p as usize - 1] {
                    ok = false;
                }
            }
        }
        t.check(ok, "wrong FASTA: every mismatch differs from that FASTA");
    }
    // Lengths that differ from the BAM header.
    let mut curto = d.genoma.clone();
    curto[1].truncate(1500);
    let fa_curto = t.p("curto.fa");
    escreve_fasta(&fa_curto, &curto, &["7", "8"], 60);
    let o = ReadsOptions { reference: Some(fa_curto.clone()), ..Default::default() };
    if let Some(v) = t.json(region_reads(&d.bam, "chr8:1450-1600", &o), "short FASTA") {
        let w = avisos(&v);
        t.check(w.iter().any(|a| a.contains("has 1,500 bases") && a.contains("says chr8 has 2,000")), "short FASTA: warns that the lengths differ (and matched 8 to chr8)");
        t.check(w.iter().any(|a| a.contains("past the end of 8")), "short FASTA: warns that the window goes past its end");
        let r = v["reference"].as_str().unwrap_or("");
        t.check(r.len() == 151 && r[51..].bytes().all(|b| b == b'N') && r[..51].as_bytes() == &curto[1][1449..1500], "short FASTA: bases beyond its end are N");
    }
    if let Some(v) = t.json(region_reads(&d.bam, "chr7:100-200", &opcoes(Some(&fa_curto), 500)), "FASTA chr names") {
        t.eq(&v["mismatch_source"], &json!("reference"), "FASTA named 7, 8: chr7 of the BAM matched to 7");
    }
    let fa_outro = t.p("outro.fa");
    escreve_fasta(&fa_outro, &d.genoma[..1], &["chr7"], 60);
    if let Some(v) = t.json(region_reads(&d.bam, "chr8:100-200", &opcoes(Some(&fa_outro), 500)), "FASTA without the sequence") {
        t.check(avisos(&v).iter().any(|a| a.contains("chr8 is not in the reference")), "FASTA without the sequence: warns, and falls back to MD");
        t.eq(&v["mismatch_source"], &json!("MD"), "FASTA without the sequence: mismatches from MD");
    }
}

fn testa_modos(t: &mut T, d: &Dados) {
    // Indexed BAM, BAM without index, SAM and SAM.gz: the same view.
    let regioes = ["chr7:1-600", "chr7:1201-2400", "chr8:300-1300", "7:2500", "chr8:1990-2000"];
    for reg in regioes {
        for o in [opcoes(Some(&d.fasta), 100_000), opcoes(None, 37)] {
            let rot = format!("modes {reg} (max_reads {})", o.max_reads);
            let a = t.json(region_reads(&d.bam, reg, &o), &rot);
            let b = t.json(region_reads(&d.bam_sem, reg, &o), &rot);
            let c = t.json(region_reads(&d.sam, reg, &o), &rot);
            let e = t.json(region_reads(&d.sam_gz, reg, &o), &rot);
            if let (Some(a), Some(b), Some(c), Some(e)) = (a, b, c, e) {
                t.check(a["indexed"] == json!(true) && b["indexed"] == json!(false) && c["indexed"] == json!(false), &format!("{rot}: indexed flag"));
                t.check(sem_modo(&a) == sem_modo(&b), &format!("{rot}: BAM without index = indexed BAM"));
                t.check(sem_modo(&a) == sem_modo(&c), &format!("{rot}: SAM = indexed BAM"));
                t.check(sem_modo(&a) == sem_modo(&e), &format!("{rot}: SAM.gz = indexed BAM"));
                t.check(avisos(&b).iter().any(|w| w.contains("No index (.bai)")), &format!("{rot}: BAM without index says so"));
                t.check(avisos(&c).iter().any(|w| w.contains("SAM text")), &format!("{rot}: SAM says it was read whole"));
                t.check(!avisos(&a).iter().any(|w| w.contains("index")), &format!("{rot}: indexed BAM has no index warning"));
            }
        }
    }
    // A point and its flank.
    if let Some(v) = t.json(region_reads(&d.bam, "7:1500", &ReadsOptions { flank: 60, ..Default::default() }), "point") {
        t.eq(&v["region"], &json!({"chrom": "chr7", "start": 1440, "end": 1560, "focus": 1500, "length": 121, "sequence_length": 3000, "requested": "7:1500"}), "point: a window of 60 bases on each side, chr7 matched");
    }
    if let Some(v) = t.json(region_reads(&d.bam, "chr7:10", &ReadsOptions { flank: 60, ..Default::default() }), "point near the start") {
        t.check(u(&v["region"]["start"]) == 1 && u(&v["region"]["end"]) == 70, "point near the start: clamped at 1");
    }
    if let Some(v) = t.json(region_reads(&d.bam, "chr8:1,990-", &ReadsOptions::default()), "open end") {
        t.check(u(&v["region"]["start"]) == 1990 && u(&v["region"]["end"]) == 2000, "open end: to the end of the sequence");
    }
    if let Some(v) = t.json(region_reads(&d.bam, "chr8:1999", &ReadsOptions { flank: 0, ..Default::default() }), "one base") {
        t.check(u(&v["region"]["length"]) == 1 && v["coverage"].as_array().map(|a| a.len()) == Some(1), "flank 0: a window of one base");
    }

    // Downsampling: deterministic, coverage untouched, mates together.
    let reg = "chr7:1-3000";
    let todos = t.json(region_reads(&d.bam, reg, &opcoes(None, 100_000)), "downsampling: all");
    let a = t.json(region_reads(&d.bam, reg, &opcoes(None, 150)), "downsampling: 150");
    let b = t.json(region_reads(&d.bam, reg, &opcoes(None, 150)), "downsampling: 150 again");
    let z = t.json(region_reads(&d.bam, reg, &opcoes(None, 0)), "downsampling: 0");
    if let (Some(todos), Some(a), Some(b), Some(z)) = (todos, a, b, z) {
        let total = u(&todos["reads_total"]);
        t.check(total > 1000 && u(&todos["reads_shown"]) == total && todos["downsampled"] == json!(false), "downsampling: all reads drawn when they fit");
        t.check(u(&a["reads_total"]) == total && u(&a["reads_shown"]) == 150 && a["downsampled"] == json!(true), "downsampling: 150 drawn, the total said");
        t.check(avisos(&a).iter().any(|w| w.contains("a fixed sample of 150 is drawn")), "downsampling: warning");
        t.check(a == b, "downsampling: the same call gives the same reads");
        t.check(a["coverage"] == todos["coverage"] && z["coverage"] == todos["coverage"], "downsampling: the coverage counts every read");
        t.check(u(&z["reads_shown"]) == 0 && u(&z["rows"]) == 0, "max_reads 0: coverage only");
        confere_linhas(t, &a, "downsampling");
        // The sample: every drawn read is a read of the full view, unchanged
        // except for its row; mates stay together (at most one pair cut).
        let mut chaves: HashMap<(String, u64), Value> = HashMap::new();
        for r in todos["reads"].as_array().unwrap() {
            let mut r = r.clone();
            r.as_object_mut().unwrap().remove("row");
            chaves.insert((r["name"].as_str().unwrap().to_string(), u(&r["flag"])), r);
        }
        let mut ok = true;
        let mut nomes: HashMap<String, usize> = HashMap::new();
        for r in a["reads"].as_array().unwrap() {
            let mut r = r.clone();
            r.as_object_mut().unwrap().remove("row");
            let k = (r["name"].as_str().unwrap().to_string(), u(&r["flag"]));
            ok &= chaves.get(&k) == Some(&r);
            *nomes.entry(k.0).or_default() += 1;
        }
        t.check(ok, "downsampling: the drawn reads are reads of the full view");
        let mut pares_total: HashMap<String, usize> = HashMap::new();
        for (n, _) in chaves.keys() {
            *pares_total.entry(n.clone()).or_default() += 1;
        }
        let cortados = nomes.iter().filter(|(n, &k)| k == 1 && pares_total.get(*n) == Some(&2)).count();
        let pares = nomes.values().filter(|&&k| k == 2).count();
        t.check(cortados <= 1 && pares > 5, &format!("downsampling: mates kept together ({pares} pairs drawn whole, {cortados} cut)"));
    }

    // Filters: counts of what was left out.
    let o = ReadsOptions { min_mapq: 30, ..Default::default() };
    if let Some(v) = t.json(region_reads(&d.bam, "chr8:1-2000", &o), "filters") {
        let fora_flag = d.reads.iter().filter(|s| s.crom == 1 && s.flag & 0xF04 != 0).count() as u64;
        let fora_mapq = d.reads.iter().filter(|s| s.crom == 1 && s.flag & 0xF04 == 0 && s.mapq < 30).count() as u64;
        let dentro = d.reads.iter().filter(|s| s.crom == 1 && passa(s, &o)).count() as u64;
        t.check(u(&v["filtered_flag"]) == fora_flag && u(&v["filtered_mapq"]) == fora_mapq && u(&v["reads_total"]) == dentro, &format!("filters: {fora_flag} by FLAG, {fora_mapq} by MAPQ, {dentro} kept"));
    }

    // A SAM without @SQ: names matched in any convention, no length known.
    let sem_sq = t.p("sem_sq.sam");
    let linhas: Vec<String> = d.reads.iter().filter(|s| s.crom == 0).map(|s| s.sam(&NOMES)).collect();
    grava(&sem_sq, &linhas);
    let a = t.json(region_reads(&sem_sq, "7:800-900", &opcoes(None, 100_000)), "SAM without @SQ");
    let b = t.json(region_reads(&d.bam, "chr7:800-900", &opcoes(None, 100_000)), "SAM without @SQ (BAM)");
    if let (Some(a), Some(b)) = (a, b) {
        t.check(a["reads"] == b["reads"] && a["coverage"] == b["coverage"], "SAM without @SQ: 7 matched to chr7, same reads");
        t.check(a["region"]["chrom"] == json!("chr7") && a["region"]["sequence_length"].is_null(), "SAM without @SQ: the name found in the file, no length");
    }
    t.erro(region_reads(&sem_sq, "chr7:800-", &ReadsOptions::default()), &["has no end", "does not declare the length"], "SAM without @SQ: an open end is refused");
}

fn testa_profundidade(t: &mut T, d: &Dados) {
    let regioes: Vec<String> = ["chr7:1-3000", "7:250-1249", "chr7:1000-1500", "chr7:1200-1300", "chr8", "chr8:1990-2000", "chr8:5"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let limites = [(0u8, 1u64, 3000u64), (0, 250, 1249), (0, 1000, 1500), (0, 1200, 1300), (1, 1, 2000), (1, 1990, 2000), (1, 5, 5)];
    for (bin, mq) in [(100u32, 0u8), (7, 30), (1, 0), (5000, 0)] {
        // Brute force: the aligned bases of the generator, position by position.
        let mut prof: Vec<Vec<u64>> = TAMANHOS.iter().map(|&l| vec![0; l + 1]).collect();
        for s in d.reads.iter().filter(|s| s.flag & 0x704 == 0 && s.mapq >= mq) {
            for &(p, _) in &s.alinhadas {
                prof[s.crom][p as usize] += 1;
            }
        }
        let rot = format!("depth_profile bin {bin}, min_mapq {mq}");
        let a = t.json(depth_profile(&d.bam, &regioes, bin, mq), &rot);
        let b = t.json(depth_profile(&d.sam, &regioes, bin, mq), &rot);
        let c = t.json(depth_profile(&d.bam_sem, &regioes, bin, mq), &rot);
        let (Some(a), Some(b), Some(c)) = (a, b, c) else { continue };
        t.check(a["regions"] == b["regions"] && a["regions"] == c["regions"], &format!("{rot}: indexed BAM = SAM = BAM without index"));
        t.check(a["indexed"] == json!(true) && b["indexed"] == json!(false), &format!("{rot}: indexed flag"));
        let regs = a["regions"].as_array().unwrap();
        t.eq(&regs.len(), &regioes.len(), &format!("{rot}: one result per region, in order"));
        let mut ok = true;
        let mut nbins = 0;
        for ((r, &(c, ini, fim)), texto) in regs.iter().zip(&limites).zip(&regioes) {
            ok &= r["region"] == json!(texto) && r["chrom"] == json!(NOMES[c as usize]) && u(&r["start"]) == ini && u(&r["end"]) == fim;
            let mut medias = Vec::new();
            let mut a0 = ini;
            let bins = r["bins"].as_array().unwrap();
            for x in bins {
                let (s, e, m) = (u(&x[0]), u(&x[1]), x[2].as_f64().unwrap());
                let e_esp = (a0 + u64::from(bin) - 1).min(fim);
                let soma: u64 = (s..=e).map(|p| prof[c as usize][p as usize]).sum();
                let esp = soma as f64 / (e - s + 1) as f64;
                if s != a0 || e != e_esp || (m - esp).abs() > 6e-5 {
                    println!("FAIL  {rot}: {texto} bin {s}-{e} mean {m}, expected {a0}-{e_esp} {esp}");
                    ok = false;
                }
                medias.push(esp);
                a0 = e + 1;
                nbins += 1;
            }
            ok &= a0 == fim + 1;
            let soma: u64 = (ini..=fim).map(|p| prof[c as usize][p as usize]).sum();
            let media = soma as f64 / (fim - ini + 1) as f64;
            medias.sort_by(|x, y| x.partial_cmp(y).unwrap());
            let n = medias.len();
            let med = if n % 2 == 1 { medias[n / 2] } else { (medias[n / 2 - 1] + medias[n / 2]) / 2.0 };
            let leituras = d
                .reads
                .iter()
                .filter(|s| s.crom == c as usize && s.flag & 0x704 == 0 && s.mapq >= mq && s.alinhadas.iter().any(|a| a.0 >= ini && a.0 <= fim))
                .count() as u64;
            if (r["mean"].as_f64().unwrap() - media).abs() > 6e-5 || (r["median"].as_f64().unwrap() - med).abs() > 6e-5 || u(&r["reads"]) != leituras {
                println!("FAIL  {rot}: {texto}: mean {} median {} reads {}; expected {media} {med} {leituras}", r["mean"], r["median"], r["reads"]);
                ok = false;
            }
        }
        t.check(ok, &format!("{rot}: {nbins} bins, means, medians and read counts equal to the per-position sums"));
    }
    t.erro(depth_profile(&d.bam, &["chr7:1-100".to_string()], 0, 0), &["at least 1"], "depth_profile error: bin 0");
    t.erro(depth_profile(&d.bam, &[], 10, 0), &["No region"], "depth_profile error: no region");
    t.erro(depth_profile(&d.bam, &["chr9:1-100".to_string()], 10, 0), &["chr9", "chr7, chr8"], "depth_profile error: unknown sequence");
    t.erro(depth_profile(&d.bam, &["chr7:5000-6000".to_string()], 10, 0), &["starts after the end", "3,000"], "depth_profile error: beyond the end");
}

fn testa_sequencia(t: &mut T, d: &Dados, rng: &mut Rng) {
    // Same genome, four ways: plain, plain + .fai, BGZF, irregular lines.
    let simples = t.p("seq.fa");
    escreve_fasta(&simples, &d.genoma, &NOMES, 61);
    let com_fai = t.p("seq_fai.fa");
    escreve_fasta(&com_fai, &d.genoma, &["7", "8"], 50);
    let mut fai = String::new();
    let mut off = 0usize;
    for (g, n) in d.genoma.iter().zip(["7", "8"]) {
        off += format!(">{n} test sequence\n").len();
        fai.push_str(&format!("{n}\t{}\t{off}\t50\t51\n", g.len()));
        off += g.len() + g.len().div_ceil(50);
    }
    std::fs::write(format!("{com_fai}.fai"), fai).unwrap();
    let gz = t.p("seq.fa.gz");
    let mut linhas = Vec::new();
    for (g, n) in d.genoma.iter().zip(NOMES) {
        linhas.push(format!(">{n}"));
        linhas.extend(g.chunks(80).map(|c| String::from_utf8(c.to_vec()).unwrap()));
    }
    grava(&gz, &linhas);
    let irregular = t.p("irregular.fa");
    let mut s = String::new();
    for (g, n) in d.genoma.iter().zip(NOMES) {
        s.push_str(&format!(">{n}\r\n"));
        let mut i = 0;
        while i < g.len() {
            let k = (1 + rng.below(90) as usize).min(g.len() - i);
            s.push_str(&String::from_utf8(g[i..i + k].to_ascii_lowercase()).unwrap());
            s.push_str("\r\n");
            i += k;
        }
    }
    std::fs::write(&irregular, s).unwrap();
    let mut n = 0;
    let mut ok = true;
    for fa in [&simples, &com_fai, &gz, &irregular] {
        for _ in 0..60 {
            let c = rng.below(2) as usize;
            let a = 1 + rng.below(TAMANHOS[c] as u64);
            let b = (a + rng.below(400)).min(TAMANHOS[c] as u64);
            let nome = if rng.chance(50) { NOMES[c].to_string() } else { NOMES[c][3..].to_string() };
            let esp = std::str::from_utf8(&d.genoma[c][a as usize - 1..b as usize]).unwrap().to_string();
            match fetch_sequence(fa, &format!("{nome}:{a}-{b}")) {
                Ok(x) if x == esp => n += 1,
                Ok(x) => {
                    println!("FAIL  fetch_sequence {fa} {nome}:{a}-{b}: {} bases, expected {}", x.len(), esp.len());
                    ok = false;
                }
                Err(e) => {
                    println!("FAIL  fetch_sequence {fa} {nome}:{a}-{b}: {e}");
                    ok = false;
                }
            }
        }
        let um = fetch_sequence(fa, "chr7:5").unwrap_or_default();
        let fim = fetch_sequence(fa, "chr8:1,990-").unwrap_or_default();
        let todo = fetch_sequence(fa, "8").unwrap_or_default();
        ok &= um.as_bytes() == &d.genoma[0][4..5] && fim.as_bytes() == &d.genoma[1][1989..] && todo.as_bytes() == &d.genoma[1][..];
    }
    t.check(ok && n == 240, &format!("fetch_sequence: {n} random regions, one base, open end and whole sequence equal in plain, .fai, BGZF and irregular FASTA"));
    t.erro(fetch_sequence(&simples, "chr9:1-10"), &["chr9", "seq.fa"], "fetch_sequence error: unknown sequence");
    t.erro(fetch_sequence(&simples, "chr7:2990-3010"), &["ends after the end", "3,000", "chr7:2990-3000"], "fetch_sequence error: past the end");
    t.erro(fetch_sequence(&simples, "chr7:4000-4010"), &["starts after the end"], "fetch_sequence error: start past the end");
    t.erro(fetch_sequence(&t.p("nao_existe.fa"), "chr7:1-10"), &["does not exist"], "fetch_sequence error: missing file");
}

fn testa_erros(t: &mut T, d: &Dados) {
    let o = ReadsOptions::default();
    t.erro(region_reads(&t.p("nao_existe.bam"), "chr7:1-10", &o), &["does not exist", "nao_existe.bam"], "error: missing file");
    let cram = t.p("x.cram");
    std::fs::write(&cram, b"CRAM").unwrap();
    t.erro(region_reads(&cram, "chr7:1-10", &o), &["CRAM", "samtools view -b"], "error: CRAM");
    t.erro(region_reads(&d.bam, "chr7", &o), &["whole sequence", "chr1:1000-1200"], "error: a whole sequence");
    t.erro(region_reads(&d.bam, "chr9:1-100", &o), &["chr9", "chr7, chr8", "leituras.bam"], "error: unknown sequence");
    t.erro(region_reads(&d.bam, "chr7:5000", &o), &["starts after the end", "3,000"], "error: position past the end");
    t.erro(region_reads(&d.bam, "chr7:abc", &o), &["not valid"], "error: unreadable region");
    let grande = t.p("grande.sam");
    grava(&grande, &["@SQ\tSN:chrG\tLN:500000".to_string(), "a\t0\tchrG\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII".to_string()]);
    t.erro(region_reads(&grande, "chrG:1-200000", &o), &["200,000 bases", "at most 100,000", "depth_profile"], "error: window too wide");
    let ruim = t.p("ruim.sam");
    grava(&ruim, &["@SQ\tSN:chr7\tLN:3000".to_string(), "q1\t0\tchr7\t5\t60\t10Q\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII".to_string()]);
    t.erro(region_reads(&ruim, "chr7:1-100", &o), &["q1", "10Q", "ruim.sam"], "error: invalid CIGAR");
    // Hand-made reads: every field at once.
    let mao = t.p("mao.sam");
    let g = &d.genoma[0];
    let refs = |a: usize, b: usize| String::from_utf8(g[a - 1..b].to_vec()).unwrap();
    // r1: 3H2S5M2I5M3D4M100N6M4S at 101: one mismatch at 103 (seq), MD with it.
    let mut seq = String::from("TT");
    let mut m5 = refs(101, 105).into_bytes();
    let orig = m5[2];
    m5[2] = if orig == b'A' { b'C' } else { b'A' };
    seq.push_str(std::str::from_utf8(&m5).unwrap());
    seq.push_str("GG");
    seq.push_str(&refs(106, 110));
    seq.push_str(&refs(114, 117));
    seq.push_str(&refs(218, 223));
    seq.push_str("CCCC");
    let qual: String = (0..seq.len()).map(|i| (b'!' + 20 + (i % 20) as u8) as char).collect();
    let md = format!("2{}7^{}10", orig as char, refs(111, 113));
    grava(&mao, &[
        "@SQ\tSN:chr7\tLN:3000".to_string(),
        format!("r1\t99\tchr7\t101\t42\t3H2S5M2I5M3D4M100N6M4S\t=\t300\t250\t{seq}\t{qual}\tMD:Z:{md}"),
    ]);
    if let Some(v) = t.json(region_reads(&mao, "chr7:90-240", &opcoes(None, 10)), "hand-made read") {
        let r = &v["reads"][0];
        let q103 = u64::from(qual.as_bytes()[4] - 33);
        t.eq(
            &json!([r["start"], r["end"], r["clip_left"], r["clip_right"], r["hard_clip_left"], r["hard_clip_right"], r["clip_left_seq"], r["clip_right_seq"], r["blocks"], r["mismatches"], r["insertions"], r["deletions"], r["skips"], r["mate"]]),
            &json!([101, 223, 2, 4, 3, 0, "TT", "CCCC", [[101, 110], [114, 117], [218, 223]], [[103, (m5[2] as char).to_string(), q103]], [[105, 2, "GG"]], [[111, 3]], [[118, 100]],
                {"chrom": "chr7", "start": 300, "strand": "-", "unmapped": false, "proper_pair": true, "insert_size": 250, "read_number": 1}]),
            "hand-made read: every field",
        );
        t.check(u(&v["coverage"][(111 - 90) as usize][6]) == 1 && u(&v["coverage"][(104 - 90) as usize][7]) == 0 && u(&v["coverage"][(105 - 90) as usize][7]) == 1, "hand-made read: deletion and insertion in the coverage");
        t.check(u(&v["coverage"][(150 - 90) as usize][1]) == 0, "hand-made read: no depth over the intron");
    }
}

// ---------------------------------------------------------------- 2. real data

fn testa_exemplos(t: &mut T, pasta: &str) {
    let p = |n: &str| format!("{pasta}/{n}");
    let fa = p("toy_genome.fa");
    // The pile on chr1:25,000 of the web interface's example.
    let bam = p("pileup_demo.bam");
    let o = ReadsOptions { reference: Some(fa.clone()), ..Default::default() };
    if let Some(v) = t.json(region_reads(&bam, "chr1:25000", &o), "example pileup_demo") {
        confere_pileup(t, &bam, "chr1", &v, &o, "example pileup_demo chr1:24940-25060");
        confere_linhas(t, &v, "example pileup_demo");
        let fora = v["reads"].as_array().unwrap().iter().flat_map(|r| r["mismatches"].as_array().unwrap().clone()).filter(|m| u(&m[0]) != 25000).count();
        let em = v["reads"].as_array().unwrap().iter().flat_map(|r| r["mismatches"].as_array().unwrap().clone()).filter(|m| u(&m[0]) == 25000).count();
        let pl = clinico::pileup(&bam, "chr1", 25000, &clinico::PileupOptions { min_baseq: 0, ..Default::default() }).unwrap();
        let alt = pl.depth - [pl.a, pl.c, pl.g, pl.t][match fetch_sequence(&fa, "chr1:25000").unwrap().as_str() {
            "A" => 0,
            "C" => 1,
            "G" => 2,
            _ => 3,
        }];
        println!("  pileup_demo: {} reads, {} rows, {em} mismatches at 25,000 (pileup without quality filter: {alt} non-reference)", v["reads_shown"], v["rows"]);
        t.check(fora == 0 && em == alt, "example pileup_demo: mismatches only at the SNP, as many as the pileup's non-reference bases");
    }
    // Reads copied from the toy genome: no mismatch anywhere (only the 5S clips differ).
    let bam = p("alignments.sorted.bam");
    let mut leituras = 0;
    let mut ok = true;
    for reg in ["chr1:1-100000", "chr2:1-100000", "chr10:1-100000"] {
        let reg_real = match fetch_sequence(&fa, reg.split(':').next().unwrap()) {
            Ok(s) => format!("{}:1-{}", reg.split(':').next().unwrap(), s.len().min(100_000)),
            Err(_) => continue,
        };
        if let Some(v) = t.json(region_reads(&bam, &reg_real, &ReadsOptions { reference: Some(fa.clone()), max_reads: 100_000, ..Default::default() }), "example alignments") {
            for r in v["reads"].as_array().unwrap() {
                leituras += 1;
                ok &= r["mismatches"].as_array().unwrap().is_empty();
            }
            t.check(avisos(&v).iter().any(|w| w.contains("No index")), "example alignments: no .bai, read whole and said so");
            confere_linhas(t, &v, &format!("example alignments {reg_real}"));
        }
    }
    t.check(ok && leituras > 500, &format!("example alignments: {leituras} reads copied from the genome, no mismatch against the FASTA"));
    // RNA-seq reads (introns): coverage equal to pileup.
    let bam = p("rnaseq_ctrl_1.bam");
    let mut n = 0;
    let mut saltos = 0;
    for reg in ["chr1:1-3000", "chr1:20001-23000", "chr2:5001-8000"] {
        let o = ReadsOptions { max_reads: 100_000, ..Default::default() };
        if let Some(v) = t.json(region_reads(&bam, reg, &o), "example rnaseq") {
            n += confere_pileup(t, &bam, reg.split(':').next().unwrap(), &v, &o, &format!("example rnaseq {reg}"));
            saltos += v["reads"].as_array().unwrap().iter().filter(|r| !r["skips"].as_array().unwrap().is_empty()).count();
            confere_linhas(t, &v, &format!("example rnaseq {reg}"));
        }
    }
    println!("  rnaseq_ctrl_1: {n} positions equal to pileup, {saltos} reads with an intron (N)");
}

/// The aligned (M, =, X) reference positions and query indices of a SAM line.
fn alinhadas_da_linha(linha: &str) -> Option<(String, u16, u8, u64, Vec<(u64, usize)>, Vec<u8>)> {
    let c: Vec<&str> = linha.split('\t').collect();
    let flag: u16 = c[1].parse().ok()?;
    let pos: u64 = c[3].parse().ok()?;
    if c[5] == "*" || flag & 4 != 0 || pos == 0 {
        return None;
    }
    let mut v = Vec::new();
    let (mut r, mut q) = (pos, 0usize);
    let mut n = 0u64;
    for b in c[5].bytes() {
        if b.is_ascii_digit() {
            n = n * 10 + u64::from(b - b'0');
            continue;
        }
        match b {
            b'M' | b'=' | b'X' => {
                for k in 0..n {
                    v.push((r + k, q + k as usize));
                }
                r += n;
                q += n as usize;
            }
            b'D' | b'N' => r += n,
            b'I' | b'S' => q += n as usize,
            _ => {}
        }
        n = 0;
    }
    Some((c[0].to_string(), flag, c[4].parse().unwrap_or(0), pos, v, c[9].as_bytes().to_vec()))
}

fn testa_real(t: &mut T, bam: &str) {
    let inicio = std::time::Instant::now();
    // The stretches the BAM covers.
    let mut r = externo::abre_leitura(bam).unwrap();
    let mut linha = String::new();
    let mut spans: Vec<(String, u64, u64)> = Vec::new();
    let mut com_md = 0u64;
    let mut total = 0u64;
    loop {
        linha.clear();
        if r.read_line(&mut linha).unwrap() == 0 {
            break;
        }
        let l = linha.trim_end();
        if l.starts_with('@') {
            continue;
        }
        let Some((_, _, _, pos, al, _)) = alinhadas_da_linha(l) else { continue };
        total += 1;
        com_md += u64::from(l.contains("\tMD:Z:"));
        let crom = l.split('\t').nth(2).unwrap().to_string();
        let fim = al.last().map(|a| a.0).unwrap_or(pos);
        match spans.last_mut() {
            Some(s) if s.0 == crom && pos <= s.2 + 1000 => s.2 = s.2.max(fim),
            _ => spans.push((crom, pos, fim)),
        }
    }
    let cobertos: u64 = spans.iter().map(|s| s.2 - s.1 + 1).sum();
    println!("  real BAM: {total} mapped reads ({com_md} with MD) in {} stretches, {cobertos} bases ({:.1} s)", spans.len(), inicio.elapsed().as_secs_f64());
    // Windows of 1,000 bases spread over the stretches: ~24,000 positions.
    let alvo = 24usize;
    let passo = (cobertos / alvo as u64).max(1000);
    let mut janelas: Vec<(String, u64, u64)> = Vec::new();
    for s in &spans {
        let mut x = s.1;
        while x + 999 <= s.2 {
            janelas.push((s.0.clone(), x, x + 999));
            x += passo;
        }
    }
    let o = ReadsOptions { max_reads: 1_000_000, ..Default::default() };
    let mut posicoes = 0;
    let mut md_ok = 0u64;
    let mut md_ruim = 0u64;
    let mut leituras = 0u64;
    let mut difs = 0u64;
    let t0 = std::time::Instant::now();
    let mut t_vista = 0f64;
    for (c, a, b) in &janelas {
        let reg = format!("{c}:{a}-{b}");
        let tv = std::time::Instant::now();
        let r = region_reads(bam, &reg, &o);
        t_vista += tv.elapsed().as_secs_f64();
        let Some(v) = t.json(r, &format!("real {reg}")) else { continue };
        posicoes += confere_pileup(t, bam, c, &v, &o, &format!("real {reg}"));
        confere_linhas(t, &v, &format!("real {reg}"));
        leituras += u(&v["reads_shown"]);
        // MD consistency: the window's reference as rebuilt from the MD tags,
        // read against the bases of every read (from the SAM lines).
        if v["mismatch_source"] == json!("MD") {
            let refj = v["reference"].as_str().unwrap().as_bytes().to_vec();
            let mut difs_por: HashMap<(String, u64, u64), HashSet<u64>> = HashMap::new();
            for rd in v["reads"].as_array().unwrap().iter().filter(|rd| rd["mismatch_source"] == json!("MD")) {
                let k = (rd["name"].as_str().unwrap().to_string(), u(&rd["flag"]), u(&rd["start"]));
                let set: HashSet<u64> = rd["mismatches"].as_array().unwrap().iter().map(|m| u(&m[0])).collect();
                difs += set.len() as u64;
                difs_por.insert(k, set);
            }
            indice::consultar(bam, &[reg.clone()], &mut |cab, l| {
                if cab {
                    return Ok(());
                }
                let Some((nome, flag, mapq, pos, al, seq)) = alinhadas_da_linha(l) else { return Ok(()) };
                if flag & 0xF04 != 0 || mapq < o.min_mapq || seq == b"*" {
                    return Ok(());
                }
                let Some(set) = difs_por.get(&(nome, u64::from(flag), pos)) else { return Ok(()) };
                for (p, qi) in al {
                    if p < *a || p > *b {
                        continue;
                    }
                    let rb = refj[(p - a) as usize];
                    let base = seq[qi].to_ascii_uppercase();
                    if rb == b'N' || base == b'N' {
                        continue;
                    }
                    if (base != rb) == set.contains(&p) {
                        md_ok += 1;
                    } else {
                        md_ruim += 1;
                    }
                }
                Ok(())
            })
            .unwrap();
        }
    }
    println!(
        "  real BAM: {} windows, {posicoes} positions equal to pileup, {leituras} reads drawn, {difs} mismatches from MD ({:.1} s; region_reads itself {:.2} s, {:.0} ms per 1,000-base window)",
        janelas.len(),
        t0.elapsed().as_secs_f64(),
        t_vista,
        1000.0 * t_vista / janelas.len().max(1) as f64
    );
    t.check(posicoes >= 20_000, &format!("real BAM: {posicoes} positions compared with pileup"));
    if com_md > 0 {
        t.check(md_ruim == 0 && md_ok > 100_000, &format!("real BAM: MD-rebuilt reference agrees with every aligned base ({md_ok} bases, {md_ruim} disagree)"));
    }
    // The widest window allowed, over the longest stretch.
    if let Some(s) = spans.iter().max_by_key(|s| s.2 - s.1) {
        let fim = (s.1 + navegador::MAX_WINDOW - 1).min(s.2);
        let reg = format!("{}:{}-{}", s.0, s.1, fim);
        let tv = std::time::Instant::now();
        if let Some(j) = t.ok_ou(region_reads(bam, &reg, &ReadsOptions::default()), "real widest window") {
            let v: Value = serde_json::from_str(&j).unwrap();
            println!(
                "  real BAM {reg} ({} bases): {} reads, {} drawn in {} rows, JSON {:.1} MB, {:.2} s",
                fim - s.1 + 1,
                v["reads_total"],
                v["reads_shown"],
                v["rows"],
                j.len() as f64 / 1e6,
                tv.elapsed().as_secs_f64()
            );
            confere_linhas(t, &v, "real widest window");
        }
    }
    // Indexed and streamed: the same view.
    let copia = t.p("real_sem_indice.bam");
    std::fs::copy(bam, &copia).unwrap();
    for (c, a, b) in janelas.iter().step_by((janelas.len() / 3).max(1)).take(3) {
        let reg = format!("{c}:{a}-{b}");
        let x = t.json(region_reads(bam, &reg, &ReadsOptions::default()), "real indexed");
        let y = t.json(region_reads(&copia, &reg, &ReadsOptions::default()), "real streamed");
        if let (Some(x), Some(y)) = (x, y) {
            t.check(sem_modo(&x) == sem_modo(&y), &format!("real BAM {reg}: indexed = streamed ({} reads)", x["reads_shown"]));
        }
    }
    // depth_profile against per-position sums of the aligned bases.
    let regioes: Vec<String> = spans.iter().take(12).map(|s| format!("{}:{}-{}", s.0, s.1, s.2)).collect();
    let t1 = std::time::Instant::now();
    if let Some(v) = t.json(depth_profile(bam, &regioes, 100, 0), "real depth_profile") {
        let mut ok = true;
        let mut nbins = 0;
        for (r, (c, a, b)) in v["regions"].as_array().unwrap().iter().zip(spans.iter()) {
            let mut prof = vec![0u64; (b - a + 1) as usize];
            indice::consultar(bam, &[format!("{c}:{a}-{b}")], &mut |cab, l| {
                if cab {
                    return Ok(());
                }
                let Some((_, flag, _, _, al, _)) = alinhadas_da_linha(l) else { return Ok(()) };
                if flag & 0x704 != 0 {
                    return Ok(());
                }
                for (p, _) in al {
                    if p >= *a && p <= *b {
                        prof[(p - a) as usize] += 1;
                    }
                }
                Ok(())
            })
            .unwrap();
            for x in r["bins"].as_array().unwrap() {
                let (s, e, m) = (u(&x[0]), u(&x[1]), x[2].as_f64().unwrap());
                let soma: u64 = prof[(s - a) as usize..=(e - a) as usize].iter().sum();
                let esp = soma as f64 / (e - s + 1) as f64;
                if (m - esp).abs() > 6e-5 {
                    println!("FAIL  real depth_profile {c}:{s}-{e}: {m}, expected {esp}");
                    ok = false;
                }
                nbins += 1;
            }
        }
        t.check(ok && nbins > 100, &format!("real depth_profile: {nbins} bins over {} regions equal to per-position sums ({:.1} s)", regioes.len(), t1.elapsed().as_secs_f64()));
        let txt = v["text"].as_str().unwrap_or("");
        println!("{}", txt.lines().take(6).collect::<Vec<_>>().join("\n"));
    }
}

// ---------------------------------------------------------------- main

fn main() {
    let dir = std::env::temp_dir().join(format!("mmgen_navegador_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut t = T { ok: 0, falhas: 0, dir: dir.clone() };
    let inicio = std::time::Instant::now();
    let mut rng = Rng(20261116);
    let d = gera_dados(&t, &mut rng);
    println!("synthetic data: {} reads on chr7 (3,000 bp) and chr8 (2,000 bp)", d.reads.len());
    testa_cobertura_e_leituras(&mut t, &d);
    testa_modos(&mut t, &d);
    testa_profundidade(&mut t, &d);
    testa_sequencia(&mut t, &d, &mut rng);
    testa_erros(&mut t, &d);
    if let Ok(p) = std::env::var("NAV_EXAMPLES") {
        println!("examples in {p}");
        testa_exemplos(&mut t, &p);
    }
    if let Ok(b) = std::env::var("NAV_REAL_BAM") {
        println!("real BAM {b}");
        testa_real(&mut t, &b);
    }
    let _ = std::fs::remove_dir_all(&dir);
    println!("navegador tests: {} passed, {} failed ({:.1} s)", t.ok, t.falhas, inicio.elapsed().as_secs_f64());
    if t.falhas > 0 {
        std::process::exit(1);
    }
}
