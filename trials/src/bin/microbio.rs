#![allow(dead_code)]
// Trials of src/microbio.rs: reference sketches, gather, read classification
// and resistance genes.
//
// Every result is compared with an oracle written here in the plainest way,
// on text: the hashes by slicing each k-mer and its reverse complement as
// strings (the hash itself, MurmurHash3, checked against values that
// sourmash 4.9.4 printed); gather by a naive greedy cover over HashSets;
// the read classification by recounting each read's hashes against a HashMap
// of the database; the resistance genes against genes planted at known
// places with a known number of substitutions and indels.
//
// Usage: cargo run --release --bin microbio
//        cargo run --release --bin microbio -- real <research_data> <work dir>
//        cargo run --release --bin microbio -- real-amr-reads <research_data> <work dir>
//
// `real` runs the tools on the data of research_data/MANIFEST.md (the Zymo
// mock communities, the outbreak genomes); `real-amr-reads` the reads mode of
// the resistance genes (reads simulated from real genomes, and the real Zymo
// EVEN reads).
#[path = "../../../src/multimerge_core.rs"]
pub mod multimerge_core;
#[path = "../../../src/bgzf.rs"]
#[allow(dead_code)]
mod bgzf;
#[path = "../../../src/bam.rs"]
#[allow(dead_code)]
mod bam;
#[path = "../../../src/externo.rs"]
mod externo;
#[path = "../../../src/nomes.rs"]
mod nomes;
#[path = "../../../src/microbio.rs"]
mod microbio;

use microbio::*;
use std::collections::{HashMap, HashSet};
use std::io::{Read, Write};
use std::path::PathBuf;

// ---------------------------------------------------------------------
//  The checks
// ---------------------------------------------------------------------

struct T {
    ok: usize,
    total: usize,
    failed: Vec<String>,
    dir: PathBuf,
}

impl T {
    fn check(&mut self, cond: bool, what: impl AsRef<str>) -> bool {
        self.total += 1;
        if cond {
            self.ok += 1;
        } else {
            println!("  FAIL  {}", what.as_ref());
            self.failed.push(what.as_ref().to_string());
        }
        cond
    }

    fn p(&self, name: &str) -> String {
        self.dir.join(name).to_string_lossy().to_string()
    }

    fn write(&self, name: &str, text: &str) -> String {
        let p = self.p(name);
        std::fs::write(&p, text).unwrap();
        p
    }

    /// Compressed by the project's writer (the extension says how).
    fn write_packed(&self, name: &str, text: &str) -> String {
        let p = self.p(name);
        let mut w = externo::Escritor::cria(&p).unwrap();
        w.write_all(text.as_bytes()).unwrap();
        w.termina().unwrap();
        p
    }

    /// An error that must happen, and must say the given things.
    fn fails<R: std::fmt::Debug>(&mut self, what: &str, r: Result<R, String>, pieces: &[&str]) {
        match r {
            Ok(v) => {
                self.check(false, format!("{what}: should have failed, gave {v:?}"));
            }
            Err(e) => {
                let missing: Vec<&&str> = pieces.iter().filter(|p| !e.contains(**p)).collect();
                self.check(missing.is_empty() && e.len() > 40, format!("{what}: the message lacks {missing:?}: {e}"));
            }
        }
    }
}

fn read_text(path: &str) -> String {
    let mut s = String::new();
    externo::abre_leitura(path).unwrap().read_to_string(&mut s).unwrap();
    s
}

fn near(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

// ---------------------------------------------------------------------
//  Random DNA, mutations, reads
// ---------------------------------------------------------------------

/// The high bits of a 64-bit LCG (the low bits repeat too soon).
fn lcg(s: &mut u64) -> u64 {
    *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    *s >> 33
}

fn rnd(s: &mut u64, n: usize) -> usize {
    (lcg(s) % n as u64) as usize
}

fn random_dna(n: usize, s: &mut u64) -> String {
    (0..n).map(|_| char::from(b"ACGT"[rnd(s, 4)])).collect()
}

fn rc_text(t: &str) -> String {
    t.bytes()
        .rev()
        .map(|c| match c {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'a' => 't',
            b'c' => 'g',
            b'g' => 'c',
            b't' => 'a',
            o => char::from(o),
        })
        .collect()
}

fn other_base(b: u8, s: &mut u64) -> u8 {
    loop {
        let c = b"ACGT"[rnd(s, 4)];
        if c != b.to_ascii_uppercase() {
            return c;
        }
    }
}

/// Each base changed with probability `rate`.
fn mutate(seq: &str, rate: f64, s: &mut u64) -> String {
    let p = (rate * 1e6) as u64;
    seq.bytes().map(|b| char::from(if lcg(s) % 1_000_000 < p { other_base(b, s) } else { b })).collect()
}

/// Exactly `n` substitutions, at distinct positions at least `margin` bases
/// from either end.
fn substitute(seq: &str, n: usize, margin: usize, s: &mut u64) -> String {
    let mut b = seq.as_bytes().to_vec();
    let mut used = HashSet::new();
    while used.len() < n {
        let i = margin + rnd(s, b.len() - 2 * margin);
        if used.insert(i) {
            b[i] = other_base(b[i], s);
        }
    }
    String::from_utf8(b).unwrap()
}

fn fasta(records: &[(String, String)], width: usize) -> String {
    let mut o = String::new();
    for (name, seq) in records {
        o.push_str(&format!(">{name}\n"));
        let mut rest = seq.as_str();
        while !rest.is_empty() {
            let n = width.min(rest.len());
            o.push_str(&rest[..n]);
            o.push('\n');
            rest = &rest[n..];
        }
    }
    o
}

fn fastq(records: &[(String, String)]) -> String {
    let mut o = String::new();
    for (name, seq) in records {
        o.push_str(&format!("@{name}\n{seq}\n+\n{}\n", "I".repeat(seq.len())));
    }
    o
}

/// Paired reads from a mixture: each pair comes from source `i` with
/// probability `weights[i]` (so the weights are the shares of the DNA).
/// Names: `p{n}_{label}` (the same in both files).
#[allow(clippy::type_complexity)]
fn simulate_pairs(
    sources: &[(&str, &str)],
    weights: &[f64],
    pairs: usize,
    len: usize,
    insert: usize,
    err: f64,
    s: &mut u64,
) -> (Vec<(String, String)>, Vec<(String, String)>, Vec<usize>) {
    let total: f64 = weights.iter().sum();
    let (mut r1, mut r2, mut origin) = (Vec::new(), Vec::new(), Vec::new());
    for n in 0..pairs {
        let mut x = (lcg(s) % 1_000_000) as f64 / 1e6 * total;
        let mut src = 0;
        while src + 1 < weights.len() && x >= weights[src] {
            x -= weights[src];
            src += 1;
        }
        let g = sources[src].1;
        let start = rnd(s, g.len() - insert);
        let frag = &g[start..start + insert];
        let frag = if lcg(s) % 2 == 0 { frag.to_string() } else { rc_text(frag) };
        let a = mutate(&frag[..len], err, s);
        let b = mutate(&rc_text(&frag[insert - len..]), err, s);
        let name = format!("p{n}_{}", sources[src].0);
        r1.push((format!("{name}/1"), a));
        r2.push((format!("{name}/2"), b));
        origin.push(src);
    }
    (r1, r2, origin)
}

// ---------------------------------------------------------------------
//  The oracles
// ---------------------------------------------------------------------

/// The hashes of the canonical k-mers kept at `max_hash`, with repeats, by
/// slicing text.
fn oracle_hashes(seq: &str, k: usize, max_hash: u64) -> Vec<u64> {
    let up = seq.to_ascii_uppercase();
    let b = up.as_bytes();
    let mut out = Vec::new();
    if b.len() < k {
        return out;
    }
    for i in 0..=b.len() - k {
        let w = &b[i..i + k];
        if !w.iter().all(|c| matches!(c, b'A' | b'C' | b'G' | b'T')) {
            continue;
        }
        let rc = rc_text(std::str::from_utf8(w).unwrap());
        let m: &[u8] = if w <= rc.as_bytes() { w } else { rc.as_bytes() };
        let h = murmur64(m, 42);
        if h <= max_hash {
            out.push(h);
        }
    }
    out
}

fn oracle_set(seqs: &[&str], k: usize, max_hash: u64) -> Vec<u64> {
    let mut v: Vec<u64> = seqs.iter().flat_map(|s| oracle_hashes(s, k, max_hash)).collect();
    v.sort_unstable();
    v.dedup();
    v
}

fn oracle_counts(seqs: &[&str], k: usize, max_hash: u64) -> HashMap<u64, u32> {
    let mut m = HashMap::new();
    for s in seqs {
        for h in oracle_hashes(s, k, max_hash) {
            *m.entry(h).or_insert(0) += 1;
        }
    }
    m
}

/// (name, unique, shared, f_unique_weighted) of each match, by a naive greedy
/// cover with the same tie rule.
fn oracle_gather(sample: &HashMap<u64, u32>, refs: &[(String, HashSet<u64>)], scaled: u64, threshold: u64) -> Vec<(String, u64, u64, f64)> {
    let total: u64 = sample.values().map(|&a| u64::from(a)).sum();
    let mut left: HashSet<u64> = sample.keys().copied().collect();
    let mut used = vec![false; refs.len()];
    let mut out = Vec::new();
    loop {
        let mut best: Option<(usize, u64, u64)> = None;
        for (i, (_, set)) in refs.iter().enumerate() {
            if used[i] {
                continue;
            }
            let unique = set.iter().filter(|h| left.contains(h)).count() as u64;
            let shared = set.iter().filter(|h| sample.contains_key(h)).count() as u64;
            let better = match best {
                None => true,
                Some((b, bu, bs)) => unique > bu || (unique == bu && (shared > bs || (shared == bs && set.len() < refs[b].1.len()))),
            };
            if better {
                best = Some((i, unique, shared));
            }
        }
        let Some((i, unique, shared)) = best else { break };
        if unique == 0 || unique * scaled < threshold {
            break;
        }
        used[i] = true;
        let mut ab = 0u64;
        for h in &refs[i].1 {
            if left.remove(h) {
                ab += u64::from(sample[h]);
            }
        }
        out.push((refs[i].0.clone(), unique, shared, ab as f64 / total as f64));
    }
    out
}

// ---------------------------------------------------------------------
//  The trials
// ---------------------------------------------------------------------

fn trial_hash(t: &mut T) {
    // Printed by sourmash 4.9.4: `hash_murmur(s, 42)`, equal to
    // `mmh3.hash64(s, 42, signed=False)[0]`.
    let vectors: [(&str, u64); 9] = [
        ("", 17305828677633410339),
        ("A", 16750156190880784680),
        ("ACGT", 2597925387403686983),
        ("ACGTACGTACGTACG", 5002160178696126918),
        ("ACGTACGTACGTACGT", 4706917051267373191),
        ("ACGTACGTACGTACGTA", 9520599998506090000),
        ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 2689060699658636553),
        ("ACGTTGCAAGCTTCGAGGATCCATGCATGCC", 3091082158504387503),
        ("GATTACAGATTACAGATTACAGATTACAGATTACAGATTACA", 12203396630237239885),
    ];
    for (s, h) in vectors {
        t.check(murmur64(s.as_bytes(), 42) == h, format!("murmur64({s:?}) = {}, sourmash says {h}", murmur64(s.as_bytes(), 42)));
    }
    // sourmash's MinHash(scaled=...)._max_hash
    for (sc, mh) in [(1u64, u64::MAX), (2, 9223372036854775808), (10, 1844674407370955264), (100, 184467440737095520), (1000, 18446744073709552), (2000, 9223372036854776)] {
        t.check(max_hash_for_scaled(sc) == mh, format!("max_hash for scaled {sc}"));
    }
    // A k=21 sketch at scaled 1 of a sequence with an N, as sourmash made it:
    // 27 hashes, the five smallest printed.
    let seq = "ACGTNACGGTTAGCTAGCTAGGATCGATCGTAGCTAGCTAGCTAGGCATCGA";
    let v = oracle_set(&[seq], 21, u64::MAX);
    t.check(v.len() == 27 && v[..5] == [30486585988286903, 68392065818748454, 614455614217579982, 2644583448264250037, 3348832909232688170], "the text oracle against a sketch made by sourmash");
}

fn trial_sketch(t: &mut T, s: &mut u64) {
    // A genome of three sequences: a long one (two windows of 1 MiB), a short
    // one with N runs and lower case, and one shorter than k.
    let long = random_dna(1_300_000, s);
    let mut dirty: Vec<u8> = random_dna(40_000, s).into_bytes();
    for i in (0..dirty.len()).step_by(997) {
        dirty[i] = b'N';
    }
    for c in dirty[5000..9000].iter_mut() {
        *c = c.to_ascii_lowercase();
    }
    let dirty = String::from_utf8(dirty).unwrap();
    let tiny = "ACGTACGTAC".to_string();
    let recs = vec![("chr long one".to_string(), long.clone()), ("plasmid".to_string(), dirty.clone()), ("bit".to_string(), tiny.clone())];
    let files = [
        t.write("g.fa", &fasta(&recs, 60)),
        t.write_packed("g.fa.gz", &fasta(&recs, 80)),
        t.write_packed("g.fa.lz4", &fasta(&recs, usize::MAX)),
    ];
    for (k, scaled) in [(31u8, 1000u64), (21, 100), (31, 1), (32, 50)] {
        let want = oracle_set(&[&long, &dirty, &tiny], usize::from(k), max_hash_for_scaled(scaled));
        for f in &files {
            let db = t.p("g.sig");
            let sum = sketch_references(&[GenomeInput { files: vec![f.clone()], name: String::new(), host: false }], &db, &SketchOptions { k, scaled });
            let sum = match sum {
                Ok(x) => x,
                Err(e) => {
                    t.check(false, format!("sketch {f}: {e}"));
                    continue;
                }
            };
            let back = read_sketch_db(&db).unwrap();
            let sk = &back.sketches[0];
            t.check(sk.hashes == want, format!("FracMinHash k={k} scaled={scaled} {f}: {} hashes, oracle {}", sk.hashes.len(), want.len()));
            t.check(
                sk.name == "chr long one" && sk.sequences == 3 && sk.bases == (long.len() + dirty.len() + tiny.len()) as u64 && !sk.host,
                format!("name and counts k={k} {f}: {:?} {} {}", sk.name, sk.sequences, sk.bases),
            );
            t.check(sum.genomes[0].n_hashes == want.len() as u64 && back.k == k && back.scaled == scaled, "summary of the sketch");
        }
    }
    // The valid k-mers counted: long + dirty pieces between N's.
    let db = t.p("g.sig");
    let back = read_sketch_db(&db).unwrap();
    let expect_kmers: u64 = [&long, &dirty, &tiny]
        .iter()
        .map(|x| x.as_bytes().windows(32).filter(|w| w.iter().all(|c| b"ACGTacgt".contains(c))).count() as u64)
        .sum();
    t.check(back.sketches[0].kmers == expect_kmers, format!("valid k-mers: {} vs {expect_kmers}", back.sketches[0].kmers));

    // Several genomes, names, host flag, two files making one genome.
    let a = random_dna(30_000, s);
    let b1 = random_dna(20_000, s);
    let b2 = random_dna(15_000, s);
    let fa = t.write("a.fa", &fasta(&[("A genome".into(), a.clone())], 70));
    let fb1 = t.write("b1.fa", &fasta(&[("B chromosome".into(), b1.clone())], 70));
    let fb2 = t.write_packed("b2.fa.gz", &fasta(&[("B plasmid".into(), b2.clone())], 70));
    let db = t.p("two.sig");
    let genomes = vec![
        GenomeInput { files: vec![fa.clone()], name: "Alpha".into(), host: false },
        GenomeInput { files: vec![fb1.clone(), fb2.clone()], name: String::new(), host: true },
    ];
    let sum = sketch_references(&genomes, &db, &SketchOptions { k: 21, scaled: 10 }).unwrap();
    let back = read_sketch_db(&db).unwrap();
    t.check(
        back.sketches.len() == 2
            && back.sketches[0].name == "Alpha"
            && back.sketches[1].name == "B chromosome"
            && back.sketches[1].host
            && back.sketches[1].files == vec![fb1.clone(), fb2.clone()]
            && back.sketches[1].hashes == oracle_set(&[&b1, &b2], 21, max_hash_for_scaled(10))
            && back.sketches[1].sequences == 2,
        "two genomes, one of two files, host flag",
    );
    t.check(describe_sketch_db(&db).unwrap() == SketchSummary { path: db.clone(), ..sum.clone() }, "describe_sketch_db = what was written");
    t.check(sum.json().contains("\"name\":\"B chromosome\"") && sum.json().contains("\"host\":true") && sum.resumo().contains("[host decoy]"), "json and resumo of the summary");

    // Downsampling: scaled 10 -> 100 is the sketch made at 100.
    let db100 = t.p("two100.sig");
    sketch_references(&genomes, &db100, &SketchOptions { k: 21, scaled: 100 }).unwrap();
    let down = load_databases(&[db.clone()], Some(100)).unwrap();
    let direct = read_sketch_db(&db100).unwrap();
    t.check(down.sketches.iter().zip(&direct.sketches).all(|(x, y)| x.hashes == y.hashes) && down.scaled == 100, "downsampling = sketching at the larger scaled");
    t.fails("finer than the database", load_databases(&[db100.clone()], Some(10)), &["scaled = 10", "scaled = 100", "coarser"]);
    let db31 = t.p("k31.sig");
    sketch_references(&genomes[..1], &db31, &SketchOptions { k: 31, scaled: 10 }).unwrap();
    t.fails("k differs", load_databases(&[db.clone(), db31.clone()], None), &["k = 21", "k = 31", &db31]);
    let both = load_databases(&[db.clone(), db100.clone()], None).unwrap();
    t.check(both.scaled == 100 && both.sketches.len() == 4 && both.sketches[0].hashes == direct.sketches[0].hashes, "two databases merged at the coarser scaled");

    // Damaged databases.
    let bytes = std::fs::read(&db).unwrap();
    let mut flipped = bytes.clone();
    flipped[bytes.len() / 2] ^= 0x40;
    std::fs::write(t.p("flipped.sig"), &flipped).unwrap();
    t.fails("a flipped byte", read_sketch_db(&t.p("flipped.sig")), &["damaged", "checksum", "sketch the references again"]);
    std::fs::write(t.p("cut.sig"), &bytes[..bytes.len() - 100]).unwrap();
    t.fails("a cut database", read_sketch_db(&t.p("cut.sig")), &["damaged or incomplete", "cut.sig"]);
    let fq = t.write("notdb.fq", "@r\nACGT\n+\nIIII\n");
    t.fails("not a database", read_sketch_db(&fq), &["not an MMGEN sketch database", "sketch_references"]);
    let mut v2 = bytes.clone();
    v2[7] = b'9';
    std::fs::write(t.p("v9.sig"), &v2).unwrap();
    t.fails("a future version", read_sketch_db(&t.p("v9.sig")), &["version", "MMGSKDB9"]);

    // Bad options and inputs.
    let g = vec![GenomeInput { files: vec![fa.clone()], ..Default::default() }];
    t.fails("k = 33", sketch_references(&g, &t.p("x.sig"), &SketchOptions { k: 33, scaled: 10 }), &["k = 33", "1 to 32"]);
    t.fails("scaled = 0", sketch_references(&g, &t.p("x.sig"), &SketchOptions { k: 31, scaled: 0 }), &["scaled", "1000"]);
    t.fails("no file", sketch_references(&[GenomeInput::default()], &t.p("x.sig"), &SketchOptions::default()), &["no file", "FASTA"]);
    t.fails(
        "missing file",
        sketch_references(&[GenomeInput { files: vec![t.p("nope.fa")], ..Default::default() }], &t.p("x.sig"), &SketchOptions::default()),
        &["nope.fa", "does not exist"],
    );
    t.fails("output over an input", sketch_references(&g, &fa, &SketchOptions::default()), &["also one of the input files"]);
    t.check(read_text(&fa).starts_with(">A genome"), "the input survived");
    let vcf = t.write("v.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\n");
    let out = t.p("v.sig");
    t.fails("a VCF as a genome", sketch_references(&[GenomeInput { files: vec![vcf], ..Default::default() }], &out, &SketchOptions::default()), &["not a FASTA nor a FASTQ", "##fileformat"]);
    t.check(!std::path::Path::new(&out).exists(), "a failed sketch leaves no database behind");
    let hidden = t.p("hidden.fa");
    std::fs::copy(&files[1], &hidden).unwrap();
    t.fails("gzip without .gz", sketch_references(&[GenomeInput { files: vec![hidden], ..Default::default() }], &out, &SketchOptions::default()), &["gzip", ".gz"]);
}

/// The synthetic community of the gather and classification trials.
struct Community {
    names: Vec<String>,
    seqs: Vec<String>,
    host: Vec<bool>,
    files: Vec<String>,
}

fn community(t: &T, s: &mut u64) -> Community {
    let g1 = random_dna(300_000, s);
    let g2 = random_dna(250_000, s);
    let g3 = mutate(&g1, 0.03, s); // a close relative of g1 (ANI ~0.97)
    let g4 = random_dna(400_000, s);
    let g5 = random_dna(200_000, s); // never in the samples
    let g6 = format!("{}{}", &g2[..125_000], random_dna(125_000, s)); // half of g2
    let host = random_dna(500_000, s);
    let names: Vec<String> = ["G1 alpha", "G2 beta", "G3 alpha relative", "G4 gamma", "G5 absent", "G6 beta chimera", "Host"].iter().map(|x| x.to_string()).collect();
    let seqs = vec![g1, g2, g3, g4, g5, g6, host];
    let hostf = vec![false, false, false, false, false, false, true];
    let files: Vec<String> = names
        .iter()
        .zip(&seqs)
        .enumerate()
        .map(|(i, (n, q))| {
            // the host in two files (as chr21 + chrM), the others in one
            if i == 6 {
                t.write("host_a.fa", &fasta(&[(format!("{n} part a"), q[..300_000].to_string())], 60));
                t.write_packed("host_b.fa.gz", &fasta(&[(format!("{n} part b"), q[300_000..].to_string())], 60))
            } else {
                t.write(&format!("ref{i}.fa"), &fasta(&[(n.clone(), q.clone())], 70))
            }
        })
        .collect();
    Community { names, seqs, host: hostf, files }
}

fn genome_inputs(t: &T, c: &Community) -> Vec<GenomeInput> {
    (0..c.names.len())
        .map(|i| GenomeInput {
            files: if i == 6 { vec![t.p("host_a.fa"), c.files[6].clone()] } else { vec![c.files[i].clone()] },
            name: c.names[i].clone(),
            host: c.host[i],
        })
        .collect()
}

fn trial_gather(t: &mut T, c: &Community, s: &mut u64) {
    let db = t.p("community100.sig");
    sketch_references(&genome_inputs(t, c), &db, &SketchOptions { k: 31, scaled: 100 }).unwrap();
    let mh = max_hash_for_scaled(100);
    // The host is two files: no k-mer spans the cut between them.
    let refs: Vec<(String, HashSet<u64>)> = c
        .names
        .iter()
        .zip(&c.seqs)
        .enumerate()
        .map(|(i, (n, q))| {
            let parts: Vec<&str> = if i == 6 { vec![&q[..300_000], &q[300_000..]] } else { vec![q.as_str()] };
            (n.clone(), oracle_set(&parts, 31, mh).into_iter().collect())
        })
        .collect();

    // A sample of known DNA shares: G1 40%, G2 25%, G4 20%, host 15%, with
    // 0.5% sequencing errors.
    let src = [("G1", c.seqs[0].as_str()), ("G2", c.seqs[1].as_str()), ("G4", c.seqs[3].as_str()), ("H", c.seqs[6].as_str())];
    let w = [0.40, 0.25, 0.20, 0.15];
    let (r1, r2, _) = simulate_pairs(&src, &w, 30_000, 150, 400, 0.005, s);
    let f1 = t.write_packed("mix_R1.fq.gz", &fastq(&r1));
    let f2 = t.write("mix_R2.fq", &fastq(&r2));
    let all: Vec<&str> = r1.iter().chain(&r2).map(|x| x.1.as_str()).collect();
    let sample = oracle_counts(&all, 31, mh);

    let out = t.p("gather.tsv");
    let opt = GatherOptions { output: Some(out.clone()), ..Default::default() };
    let g = gather(&[f1.clone(), f2.clone()], &[db.clone()], &opt).unwrap();
    let want = oracle_gather(&sample, &refs, 100, 50_000);
    let got: Vec<(String, u64, u64, f64)> = g.matches.iter().map(|m| (m.name.clone(), m.unique_hashes, m.shared_hashes, m.f_unique_weighted)).collect();
    t.check(
        got.len() == want.len() && got.iter().zip(&want).all(|(a, b)| a.0 == b.0 && a.1 == b.1 && a.2 == b.2 && near(a.3, b.3, 1e-12)),
        format!("gather = the naive greedy cover:\n    got  {got:?}\n    want {want:?}"),
    );
    t.check(g.sample_hashes == sample.len() as u64 && g.sample_abundance == sample.values().map(|&x| u64::from(x)).sum::<u64>(), "the sample sketch = the oracle's");
    // Gather ranks by hashes explained, so the largest genomes come first:
    // the host (500 kb), G4 (400 kb), G1 (300 kb), G2 (250 kb).
    let names: Vec<&str> = g.matches.iter().map(|m| m.name.as_str()).collect();
    t.check(names == ["Host", "G4 gamma", "G1 alpha", "G2 beta"], format!("the organisms of the sample, by hashes explained: {names:?}"));
    // Relative abundance among the microbes: 40 : 25 : 20.
    for (name, want) in [("G1 alpha", 0.40 / 0.85), ("G2 beta", 0.25 / 0.85), ("G4 gamma", 0.20 / 0.85)] {
        let Some(m) = g.matches.iter().find(|m| m.name == name) else {
            t.check(false, format!("{name} not found"));
            continue;
        };
        t.check(m.relative_abundance.is_some_and(|x| near(x, want, 0.02)), format!("relative abundance of {}: {:?}, truth {want:.3}", m.name, m.relative_abundance));
        // G4 is at ~4x: a few percent of its k-mers were never sequenced,
        // which lowers the plain ANI and not the corrected one.
        t.check(m.containment > 0.9 && m.ani > 0.995 && m.ani_adjusted.is_some_and(|a| a > 0.999), format!("{}: containment {} ani {} adj {:?}", m.name, m.containment, m.ani, m.ani_adjusted));
    }
    let host = &g.matches[0];
    t.check(host.host && host.relative_abundance.is_none() && near(g.host_weighted, host.f_unique_weighted, 1e-12), "the host is reported as host");
    // Errors leave about 1 - 0.995^31 of the k-mers unassigned.
    t.check(near(g.unassigned_weighted, 1.0 - 0.995f64.powi(31), 0.03), format!("unassigned weighted {} (errors)", g.unassigned_weighted));
    let sum: f64 = g.matches.iter().map(|m| m.f_unique_weighted).sum::<f64>() + g.unassigned_weighted;
    t.check(near(sum, 1.0, 1e-9), "matches + unassigned = 1");
    let tsv = read_text(&out);
    t.check(tsv.lines().count() == 5 && tsv.starts_with("rank\tname\thost\tintersect_bp\tunique_intersect_bp\t") && tsv.contains("\tG4 gamma\tfalse\t"), "the TSV output");
    // sourmash's meaning: intersect_bp with the whole sample, unique_ with what was left
    t.check(
        g.matches.iter().all(|m| m.intersect_bp == m.shared_hashes * 100 && m.unique_intersect_bp == m.unique_hashes * 100 && m.unique_intersect_bp <= m.intersect_bp),
        "intersect_bp and unique_intersect_bp",
    );
    t.check(g.json().contains("\"name\":\"Host\",\"host\":true") && g.resumo().contains("G2 beta") && g.resumo().contains("host:"), "json and resumo");

    // The two files of a pair together = one file with both.
    let both = t.write("mix_both.fq", &fastq(&r1.iter().chain(&r2).cloned().collect::<Vec<_>>()));
    let g2 = gather(&[both], &[db.clone()], &GatherOptions::default()).unwrap();
    t.check(g2.matches == g.matches, "pairs in two files = one file");

    // min_abundance 2 and a smaller threshold, against the oracle.
    let g3 = gather(&[f1.clone(), f2.clone()], &[db.clone()], &GatherOptions { min_abundance: 2, threshold_bp: 1000, ..Default::default() }).unwrap();
    let filtered: HashMap<u64, u32> = sample.iter().filter(|(_, &a)| a >= 2).map(|(&h, &a)| (h, a)).collect();
    let want3 = oracle_gather(&filtered, &refs, 100, 1000);
    let got3: Vec<(String, u64, u64)> = g3.matches.iter().map(|m| (m.name.clone(), m.unique_hashes, m.shared_hashes)).collect();
    t.check(
        got3 == want3.iter().map(|x| (x.0.clone(), x.1, x.2)).collect::<Vec<_>>() && g3.low_abundance_hashes == (sample.len() - filtered.len()) as u64,
        format!("min_abundance 2: {got3:?} vs {want3:?}"),
    );
    t.check(g3.unassigned_weighted < 0.01 && g3.matches.iter().all(|m| m.ani_adjusted.is_none()), "without singletons almost nothing is unassigned");

    // A divergent strain: reads of G1 with 2% substitutions, as the only
    // organism; the ANI must say ~0.98.
    let g1v = mutate(&c.seqs[0], 0.02, s);
    let (v1, v2, _) = simulate_pairs(&[("V", g1v.as_str())], &[1.0], 8000, 150, 400, 0.0, s);
    let fv = t.write("variant.fq", &fastq(&v1.iter().chain(&v2).cloned().collect::<Vec<_>>()));
    let gv = gather(&[fv], &[db.clone()], &GatherOptions::default()).unwrap();
    let m = &gv.matches[0];
    t.check(m.name == "G1 alpha" && near(m.ani, 0.98, 0.006) && near(m.ani_adjusted.unwrap_or(0.0), 0.98, 0.006), format!("ANI of a 2% divergent strain: {} adj {:?}", m.ani, m.ani_adjusted));
    t.check(gv.matches.len() == 1, format!("the relative G3 is not reported for G1's reads: {:?}", gv.matches.iter().map(|m| &m.name).collect::<Vec<_>>()));

    // Low coverage (0.5x) of G1, no errors: the plain ANI is too low, the
    // corrected one is right.
    let (l1, _, _) = simulate_pairs(&[("L", c.seqs[0].as_str())], &[1.0], 1000, 150, 300, 0.0, s);
    let fl = t.write("low.fq", &fastq(&l1));
    let gl = gather(&[fl], &[db.clone()], &GatherOptions { threshold_bp: 10_000, ..Default::default() }).unwrap();
    let m = &gl.matches[0];
    t.check(
        m.name == "G1 alpha" && m.ani < 0.975 && m.ani_adjusted.is_some_and(|a| a > 0.993) && m.kmer_coverage.is_some_and(|l| near(l, 0.4, 0.1)),
        format!("0.5x coverage: ani {} adjusted {:?} coverage {:?}", m.ani, m.ani_adjusted, m.kmer_coverage),
    );

    // The same at 0.5x for a genome with a 10 kb element in 40 copies (as
    // rDNA or mitochondria): the repeat pulls the mean abundance up, and must
    // not pull the coverage estimate with it.
    let mut rg = random_dna(200_000, s);
    let element = random_dna(10_000, s);
    for _ in 0..40 {
        rg.push_str(&element);
        rg.push_str(&random_dna(200, s));
    }
    let rep_db = t.p("repeat.sig");
    let rep_fa = t.write("repeat.fa", &fasta(&[("Repeat genome".to_string(), rg.clone())], 80));
    sketch_references(&[GenomeInput { files: vec![rep_fa], ..Default::default() }], &rep_db, &SketchOptions { k: 31, scaled: 100 }).unwrap();
    let (q1, _, _) = simulate_pairs(&[("R", rg.as_str())], &[1.0], rg.len() / 300, 150, 300, 0.0, s);
    let fq = t.write("repeat.fq", &fastq(&q1));
    let gr = gather(&[fq], &[rep_db], &GatherOptions { threshold_bp: 10_000, ..Default::default() }).unwrap();
    let m = &gr.matches[0];
    t.check(
        m.ani < 0.975 && m.ani_adjusted.is_some_and(|a| a > 0.993) && m.kmer_coverage.is_some_and(|l| near(l, 0.4, 0.12)) && m.mean_abundance > 1.5,
        format!("0.5x with a 40-copy repeat: ani {} adjusted {:?} coverage {:?} mean abundance {}", m.ani, m.ani_adjusted, m.kmer_coverage, m.mean_abundance),
    );
    // Too shallow to say (a handful of hashes seen twice): no correction.
    let (q2, _, _) = simulate_pairs(&[("L", c.seqs[0].as_str())], &[1.0], 60, 150, 300, 0.0, s);
    let fq2 = t.write("vlow.fq", &fastq(&q2));
    let gq = gather(&[fq2], &[db.clone()], &GatherOptions { threshold_bp: 1000, ..Default::default() }).unwrap();
    t.check(gq.matches.first().is_some_and(|m| m.name == "G1 alpha" && m.ani_adjusted.is_none() && m.kmer_coverage.is_none()), format!("0.05x: no correction: {:?}", gq.matches.first()));

    // An assembly as the sample: G2 cut in five contigs.
    let contigs: Vec<(String, String)> = (0..5).map(|i| (format!("contig{i}"), c.seqs[1][i * 50_000..(i + 1) * 50_000].to_string())).collect();
    let fa = t.write("asm.fa", &fasta(&contigs, 80));
    let ga = gather(&[fa], &[db.clone()], &GatherOptions::default()).unwrap();
    let m = &ga.matches[0];
    t.check(
        !ga.reads && m.name == "G2 beta" && m.containment > 0.99 && m.ani > 0.9997 && m.ani_adjusted.is_none() && ga.matches.len() == 1,
        format!("assembly as sample: {:?} {} {}", m.name, m.containment, m.ani),
    );

    // Nothing above the threshold.
    let gn = gather(&[f2.clone()], &[db.clone()], &GatherOptions { threshold_bp: 100_000_000, ..Default::default() }).unwrap();
    t.check(gn.matches.is_empty() && gn.unassigned == 1.0 && gn.resumo().contains("no reference explains"), "a threshold no reference reaches");
    // Downsampled on the fly to 1000 = a database made at 1000.
    let db1000 = t.p("community1000.sig");
    sketch_references(&genome_inputs(t, c), &db1000, &SketchOptions { k: 31, scaled: 1000 }).unwrap();
    let a = gather(&[f1.clone(), f2.clone()], &[db.clone()], &GatherOptions { scaled: Some(1000), ..Default::default() }).unwrap();
    let b = gather(&[f1.clone(), f2.clone()], &[db1000.clone()], &GatherOptions::default()).unwrap();
    t.check(a.matches == b.matches && a.scaled == 1000, "gather with scaled=1000 on a scaled=100 database = on a scaled=1000 database");
    // The host in its own database file.
    let dbh = t.p("host_only.sig");
    let gi = genome_inputs(t, c);
    sketch_references(&gi[6..], &dbh, &SketchOptions { k: 31, scaled: 100 }).unwrap();
    let dbm = t.p("microbes.sig");
    sketch_references(&gi[..6], &dbm, &SketchOptions { k: 31, scaled: 100 }).unwrap();
    let gs = gather(&[f1.clone(), f2.clone()], &[dbm, dbh], &GatherOptions::default()).unwrap();
    t.check(gs.matches.iter().map(|m| (&m.name, m.unique_hashes)).eq(g.matches.iter().map(|m| (&m.name, m.unique_hashes))), "host decoy in a separate database");

    // Errors.
    t.fails("no sample", gather(&[], &[db.clone()], &GatherOptions::default()), &["No sample"]);
    t.fails("no database", gather(&[f2.clone()], &[], &GatherOptions::default()), &["No sketch database"]);
    t.fails("min_abundance 0", gather(&[f2.clone()], &[db.clone()], &GatherOptions { min_abundance: 0, ..Default::default() }), &["min_abundance"]);
    t.fails("output over the sample", gather(&[f2.clone()], &[db.clone()], &GatherOptions { output: Some(f2.clone()), ..Default::default() }), &["also one of the input files"]);
    let bad = t.write("bad.fq", "@r1\nACGT\n+\nIII\n@r2\nACGT\n+\nIIII\n");
    t.fails("bad FASTQ", gather(&[bad.clone()], &[db.clone()], &GatherOptions::default()), &[&bad, "line 4", "3 characters", "4 bases"]);
    let cutq = t.write("cutq.fq", "@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nII");
    t.fails("qualities cut at the end", gather(&[cutq.clone()], &[db.clone()], &GatherOptions::default()), &[&cutq, "truncated", "line 5", "2 characters for 4 bases"]);
}

/// The oracle of one read (pair): status, references, hits, sampled k-mers.
fn oracle_call(seqs: &[&str], k: usize, mh: u64, index: &HashMap<u64, Vec<usize>>, names: &[String], host: &[bool], min_hits: u32) -> (String, String, u32, u32) {
    let mut counts: HashMap<usize, u32> = HashMap::new();
    let mut kmers = 0u32;
    for s in seqs {
        for h in oracle_hashes(s, k, mh) {
            kmers += 1;
            if let Some(v) = index.get(&h) {
                for &r in v {
                    *counts.entry(r).or_insert(0) += 1;
                }
            }
        }
    }
    let best = counts.values().copied().max().unwrap_or(0);
    if best == 0 || best < min_hits {
        return ("unclassified".into(), "-".into(), best, kmers);
    }
    let mut tied: Vec<usize> = counts.iter().filter(|(_, &c)| c == best).map(|(&r, _)| r).collect();
    tied.sort_unstable();
    if let Some(&h) = tied.iter().find(|&&r| host[r]) {
        return ("host".into(), names[h].clone(), best, kmers);
    }
    if tied.len() == 1 {
        return ("classified".into(), names[tied[0]].clone(), best, kmers);
    }
    ("ambiguous".into(), tied.iter().map(|&r| names[r].as_str()).collect::<Vec<_>>().join(";"), best, kmers)
}

fn trial_classify(t: &mut T, c: &Community, s: &mut u64) {
    let db = t.p("community10.sig");
    sketch_references(&genome_inputs(t, c), &db, &SketchOptions { k: 31, scaled: 10 }).unwrap();
    let back = read_sketch_db(&db).unwrap();
    let mh = back.max_hash;
    let mut index: HashMap<u64, Vec<usize>> = HashMap::new();
    for (r, sk) in back.sketches.iter().enumerate() {
        for &h in &sk.hashes {
            index.entry(h).or_default().push(r);
        }
    }
    // G1 30%, G2 30% (half of it shared with G6), host 20%, and 20% of an
    // organism that is in no database.
    let stranger = random_dna(200_000, s);
    let src = [("G1", c.seqs[0].as_str()), ("G2", c.seqs[1].as_str()), ("H", c.seqs[6].as_str()), ("X", stranger.as_str())];
    let (r1, r2, origin) = simulate_pairs(&src, &[0.3, 0.3, 0.2, 0.2], 12_000, 150, 400, 0.003, s);
    let f1 = t.write("cl_R1.fq", &fastq(&r1));
    let f2 = t.write_packed("cl_R2.fq.gz", &fastq(&r2));
    let (assign, free1, free2, table) = (t.p("assign.tsv"), t.p("free_R1.fq.gz"), t.p("free_R2.fq"), t.p("table.tsv"));
    let opt = ClassifyOptions {
        assignments: Some(assign.clone()),
        host_free: Some(free1.clone()),
        host_free2: Some(free2.clone()),
        output: Some(table.clone()),
        ..Default::default()
    };
    let r = classify_reads(&f1, Some(&f2), &[db.clone()], &opt).unwrap();
    let host: Vec<bool> = back.sketches.iter().map(|x| x.host).collect();
    let names: Vec<String> = back.sketches.iter().map(|x| x.name.clone()).collect();
    let lines: Vec<String> = read_text(&assign).lines().map(str::to_string).collect();
    let mut same = lines.len() == r1.len() + 1 && lines[0] == "read\tstatus\treference\thits\tsampled_kmers";
    let mut first_diff = String::new();
    let mut correct = [0usize; 4];
    let mut total = [0usize; 4];
    for (i, ((a, b), o)) in r1.iter().zip(&r2).zip(&origin).enumerate() {
        let (st, refs, hits, kmers) = oracle_call(&[&a.1, &b.1], 31, mh, &index, &names, &host, 2);
        let id = a.0.trim_end_matches("/1");
        let want = format!("{id}\t{st}\t{refs}\t{hits}\t{kmers}");
        if lines.get(i + 1) != Some(&want) {
            if same {
                first_diff = format!("line {}: got {:?}, oracle {want:?}", i + 2, lines.get(i + 1));
            }
            same = false;
        }
        total[*o] += 1;
        let ok = match o {
            0 => st == "classified" && refs == "G1 alpha",
            1 => (st == "classified" && refs == "G2 beta") || (st == "ambiguous" && refs == "G2 beta;G6 beta chimera"),
            2 => st == "host",
            _ => st == "unclassified",
        };
        correct[*o] += usize::from(ok);
    }
    t.check(same, format!("every read pair = the oracle ({first_diff})"));
    for (o, label) in ["G1", "G2", "host", "stranger"].iter().enumerate() {
        let share = correct[o] as f64 / total[o].max(1) as f64;
        t.check(share > 0.995, format!("reads of {label} called right: {share:.4} of {}", total[o]));
    }
    t.check(r.reads == 12_000 && r.unique + r.ambiguous + r.host + r.unclassified == r.reads && r.classified == r.unique + r.ambiguous, "the counts add up");
    t.check(r.ambiguous > 1000 && r.mean_sampled_kmers > 20.0 && r.warnings.is_empty(), format!("ties between G2 and its chimera: {} (sampled per pair {:.1})", r.ambiguous, r.mean_sampled_kmers));
    // The table: G2's ties go to G2 in proportion to unique reads (G6 has none
    // of its own here, so nearly all of them).
    let row = |n: &str| r.rows.iter().find(|x| x.name == n);
    let (g1, g2) = (row("G1 alpha").unwrap(), row("G2 beta").unwrap());
    t.check(near(g1.fraction_of_classified.unwrap(), 0.5, 0.03) && near(g2.fraction_of_classified.unwrap(), 0.5, 0.03), format!("G1 and G2 half each of the classified: {:?} {:?}", g1.fraction_of_classified, g2.fraction_of_classified));
    let class_sum: f64 = r.rows.iter().filter(|x| !x.host).map(|x| x.reads).sum();
    t.check(near(class_sum, r.classified as f64, 1e-6), "the table adds up to the classified pairs");
    t.check(row("Host").is_some_and(|h| h.host && h.unique_reads == r.host && h.fraction_of_classified.is_none()), "the host row");
    t.check(read_text(&table).lines().count() == r.rows.len() + 1, "the table file");
    // Host-free outputs: every pair that was not host, in order, unchanged.
    let keep: Vec<usize> = (1..lines.len()).filter(|&i| !lines[i].split('\t').nth(1).is_some_and(|x| x == "host")).map(|i| i - 1).collect();
    let want1 = fastq(&keep.iter().map(|&i| r1[i].clone()).collect::<Vec<_>>());
    let want2 = fastq(&keep.iter().map(|&i| r2[i].clone()).collect::<Vec<_>>());
    t.check(read_text(&free1) == want1 && read_text(&free2) == want2 && r.host_free_written == keep.len() as u64, "the host-free reads");
    t.check(std::fs::read(&free1).unwrap().starts_with(&[0x1f, 0x8b]), "the host-free .gz is compressed");
    t.check(r.json().contains("\"paired\":true") && r.resumo().contains("host:") && r.resumo().contains("G1 alpha"), "json and resumo");

    // Single reads: the first file alone, against the oracle.
    let r_single = classify_reads(&f1, None, &[db.clone()], &ClassifyOptions { assignments: Some(assign.clone()), ..Default::default() }).unwrap();
    let lines: Vec<String> = read_text(&assign).lines().skip(1).map(str::to_string).collect();
    let mut ok = lines.len() == r1.len();
    for (a, l) in r1.iter().zip(&lines) {
        let (st, refs, hits, kmers) = oracle_call(&[&a.1], 31, mh, &index, &names, &host, 2);
        ok &= *l == format!("{}\t{st}\t{refs}\t{hits}\t{kmers}", a.0.trim_end_matches("/1"));
    }
    t.check(ok && !r_single.paired, "single reads = the oracle");

    // min_fraction, min_hits and `only`.
    let strict = classify_reads(&f1, Some(&f2), &[db.clone()], &ClassifyOptions { min_fraction: 0.9, ..Default::default() }).unwrap();
    t.check(strict.unclassified > r.unclassified, "min_fraction 0.9 leaves more pairs unclassified");
    let none = classify_reads(&f1, Some(&f2), &[db.clone()], &ClassifyOptions { min_hits: 100_000, ..Default::default() }).unwrap();
    t.check(none.unclassified == none.reads && none.rows.is_empty(), "min_hits above any read");
    let only = classify_reads(&f1, Some(&f2), &[db.clone()], &ClassifyOptions { only: vec!["G1 alpha".into(), "Host".into()], ..Default::default() }).unwrap();
    t.check(only.ambiguous == 0 && only.rows.iter().all(|x| x.name == "G1 alpha" || x.name == "Host") && only.references == 2, "only: two references");
    t.fails("only with an unknown name", classify_reads(&f1, None, &[db.clone()], &ClassifyOptions { only: vec!["Nobody".into()], ..Default::default() }), &["Nobody", "not a reference"]);
    // A coarse database: the warning.
    let coarse = classify_reads(&f1, None, &[db.clone()], &ClassifyOptions { scaled: Some(1000), ..Default::default() }).unwrap();
    t.check(coarse.warnings.len() == 1 && coarse.warnings[0].contains("smaller scaled") && coarse.resumo().contains("warning"), "warning at scaled 1000");

    // Errors.
    let n = 50;
    let short2 = t.write("short_R2.fq", &fastq(&r2[..n]));
    t.fails("fewer reads in R2", classify_reads(&f1, Some(&short2), &[db.clone()], &ClassifyOptions::default()), &["same number of reads", "short_R2.fq", "50"]);
    let mut swapped = r2[..200].to_vec();
    swapped.swap(10, 11);
    let sw = t.write("swapped_R2.fq", &fastq(&swapped));
    let f1s = t.write("first200_R1.fq", &fastq(&r1[..200]));
    t.fails("names out of order", classify_reads(&f1s, Some(&sw), &[db.clone()], &ClassifyOptions::default()), &["not paired", "read 11", "line 41"]);
    let out1 = t.p("o1.fq");
    t.fails("host_free2 missing", classify_reads(&f1, Some(&f2), &[db.clone()], &ClassifyOptions { host_free: Some(out1.clone()), ..Default::default() }), &["host_free2"]);
    t.fails("host_free2 with single reads", classify_reads(&f1, None, &[db.clone()], &ClassifyOptions { host_free: Some(out1.clone()), host_free2: Some(t.p("o2.fq")), ..Default::default() }), &["leave host_free2 out"]);
    t.fails("same output twice", classify_reads(&f1, None, &[db.clone()], &ClassifyOptions { output: Some(out1.clone()), assignments: Some(out1.clone()), ..Default::default() }), &["same file"]);
    t.fails("output over the reads", classify_reads(&f1, None, &[db.clone()], &ClassifyOptions { host_free: Some(f1.clone()), ..Default::default() }), &["also one of the input files"]);
    t.fails("R1 = R2", classify_reads(&f1, Some(&f1), &[db.clone()], &ClassifyOptions::default()), &["same file"]);
    t.fails("min_fraction 2", classify_reads(&f1, None, &[db.clone()], &ClassifyOptions { min_fraction: 2.0, ..Default::default() }), &["min_fraction"]);
    let cut = t.write("cut_R2.fq", &fastq(&r2[..300])[..fastq(&r2[..300]).len() - 40]);
    let f300 = t.write("f300_R1.fq", &fastq(&r1[..300]));
    let free = t.p("free_cut.fq");
    t.fails(
        "R2 truncated",
        classify_reads(&f300, Some(&cut), &[db.clone()], &ClassifyOptions { host_free: Some(free.clone()), host_free2: Some(t.p("free_cut2.fq")), ..Default::default() }),
        &["cut_R2.fq", "truncated"],
    );
    t.check(!std::path::Path::new(&free).exists(), "a failed classification leaves no output behind");
}

/// One planted gene: which allele, where, which strand, and what the
/// alignment must give.
struct Plant {
    gene: &'static str,
    contig: &'static str,
    start: usize,
    end: usize,
    minus: bool,
    identity: f64,
    coverage: f64,
    status: &'static str,
}

fn trial_amr(t: &mut T, s: &mut u64) {
    // The catalog: family blaXYZ with three alleles (-1, -2 two substitutions
    // from -1, -3 thirty from -1), and genes of other families.
    let x1 = random_dna(900, s);
    let x2 = substitute(&x1, 2, 50, s);
    let x3 = substitute(&x1, 30, 50, s);
    // -4 and -5: one substitution each from -1, at two places; a sample
    // carrying both (a novel allele) is matched in part by each.
    let snp = |seq: &str, at: usize, s: &mut u64| {
        let mut b = seq.as_bytes().to_vec();
        b[at] = other_base(b[at], s);
        String::from_utf8(b).unwrap()
    };
    let x4 = snp(&x1, 300, s);
    let x5 = snp(&x1, 600, s);
    let novel = format!("{}{}", &x4[..450], &x5[450..]);
    let gb = random_dna(1200, s);
    let gc = random_dna(600, s);
    let gd = random_dna(1500, s);
    let ge = random_dna(300, s);
    let gp = random_dna(700, s); // a "plus" gene
    let mut gn = random_dna(450, s).into_bytes();
    gn[200] = b'N';
    gn[201] = b'R'; // A or G
    let gn = String::from_utf8(gn).unwrap();
    let gf = random_dna(1000, s); // a bifunctional fusion, listed twice
    let cat_genes: Vec<(&str, &str, &str, &str, &str)> = vec![
        // protein, symbol, family, product, sequence
        ("WP_000001.1", "blaXYZ-1", "blaXYZ", "class_A_beta-lactamase_XYZ-1", &x1),
        ("WP_000002.1", "blaXYZ-2", "blaXYZ", "class_A_beta-lactamase_XYZ-2", &x2),
        ("WP_000003.1", "blaXYZ-3", "blaXYZ", "class_A_beta-lactamase_XYZ-3", &x3),
        ("WP_000011.1", "blaXYZ-4", "blaXYZ", "class_A_beta-lactamase_XYZ-4", &x4),
        ("WP_000012.1", "blaXYZ-5", "blaXYZ", "class_A_beta-lactamase_XYZ-5", &x5),
        ("WP_000004.1", "aacB", "aacB", "aminoglycoside_acetyltransferase", &gb),
        ("WP_000005.1", "tetC", "tetC", "tetracycline_efflux", &gc),
        ("WP_000006.1", "qnrD", "qnrD", "quinolone_protection", &gd),
        ("WP_000007.1", "ermE", "ermE", "rRNA_methyltransferase", &ge),
        ("WP_000008.1", "virP", "virP", "a_virulence_factor", &gp),
        ("AAA00009.1", "fosN", "fosN", "fosfomycin_resistance", &gn),
    ];
    let mut cat = String::new();
    for (i, (p, sym, fam, prod, seq)) in cat_genes.iter().enumerate() {
        cat.push_str(&format!(">{p}|NG_{i:06}.1|1|1|{sym}|{fam}|{prod} NG_{i:06}.1:101-{}\n", 100 + seq.len()));
        for chunk in seq.as_bytes().chunks(80) {
            cat.push_str(&format!("{}\n", String::from_utf8_lossy(chunk).to_lowercase()));
        }
    }
    for part in 1..=2 {
        let sym = if part == 1 { "aacF" } else { "aphF" };
        cat.push_str(&format!(">WP_000010.1|NG_000010.1|{part}|2|{sym}|{sym}|bifunctional_AacF/AphF NG_000010.1:1-1000\n{}\n", gf.to_lowercase()));
    }
    let catalog = t.write("AMR_CDS.fa", &cat);
    let header = "allele\tgene_family\twhitelisted_taxa\tproduct_name\tscope\ttype\tsubtype\tclass\tsubclass\trefseq_protein_accession\trefseq_nucleotide_accession\tcurated_refseq_start\tgenbank_protein_accession\tgenbank_nucleotide_accession\tgenbank_strand\tgenbank_start\tgenbank_stop\trefseq_strand\trefseq_start\trefseq_stop\tpubmed_reference\tblacklisted_taxa\tsynonyms\thierarchy_node\tdb_version";
    let mut table = format!("{header}\n");
    let rows = [
        ("blaXYZ-1", "blaXYZ", "class A beta-lactamase XYZ-1", "core", "AMR", "AMR", "BETA-LACTAM", "CEPHALOSPORIN", "WP_000001.1", "blaXYZ-1"),
        ("blaXYZ-2", "blaXYZ", "class A beta-lactamase XYZ-2", "core", "AMR", "AMR", "BETA-LACTAM", "CARBAPENEM", "WP_000002.1", "blaXYZ-2"),
        ("blaXYZ-3", "blaXYZ", "class A beta-lactamase XYZ-3", "core", "AMR", "AMR", "BETA-LACTAM", "CEPHALOSPORIN", "WP_000003.1", "blaXYZ-3"),
        ("blaXYZ-4", "blaXYZ", "class A beta-lactamase XYZ-4", "core", "AMR", "AMR", "BETA-LACTAM", "CEPHALOSPORIN", "WP_000011.1", "blaXYZ-4"),
        ("blaXYZ-5", "blaXYZ", "class A beta-lactamase XYZ-5", "core", "AMR", "AMR", "BETA-LACTAM", "CEPHALOSPORIN", "WP_000012.1", "blaXYZ-5"),
        ("", "aacB", "aminoglycoside N-acetyltransferase AacB", "core", "AMR", "AMR", "AMINOGLYCOSIDE", "GENTAMICIN", "WP_000004.1", "aacB"),
        ("", "tetC", "tetracycline efflux MFS transporter TetC", "core", "AMR", "AMR", "TETRACYCLINE", "TETRACYCLINE", "WP_000005.1", "tetC"),
        ("", "qnrD", "quinolone resistance pentapeptide repeat protein QnrD", "core", "AMR", "AMR", "QUINOLONE", "QUINOLONE", "WP_000006.1", "qnrD"),
        ("", "ermE", "23S rRNA methyltransferase ErmE", "core", "AMR", "AMR", "MACROLIDE", "ERYTHROMYCIN", "WP_000007.1", "ermE"),
        ("", "virP", "virulence factor VirP", "plus", "VIRULENCE", "VIRULENCE", "", "", "WP_000008.1", "virP"),
        ("", "aacF/aphF", "bifunctional AacF/AphF", "core", "AMR", "AMR", "AMINOGLYCOSIDE", "AMIKACIN/GENTAMICIN", "WP_000010.1", "aacF,aphF"),
        // a point-mutation line of the same protein, which must not be picked
        ("aacB_S83L", "aacB", "aacB mutant", "core", "AMR", "POINT", "QUINOLONE", "QUINOLONE", "WP_000004.1", ""),
    ];
    for (allele, fam, prod, scope, ty, sub, class, subclass, acc, node) in rows {
        table.push_str(&format!("{allele}\t{fam}\t\t{prod}\t{scope}\t{ty}\t{sub}\t{class}\t{subclass}\t{acc}\tNG_X\tNo\tGB_{acc}\tGBN\t+\t1\t2\t+\t1\t2\t\t\t\t{node}\t2026-10-29.1\n"));
    }
    // fosN is only in the table by its GenBank accession
    table.push_str(&format!("\tfosN\t\tFosN fosfomycin resistance\tcore\tAMR\tAMR\tFOSFOMYCIN\tFOSFOMYCIN\t\tNG_X\tNo\tAAA00009.1\tGBN\t+\t1\t2\t+\t1\t2\t\t\t\tfosN\t2026-10-29.1\n"));
    let tablef = t.write("ReferenceGeneCatalog.txt", &table);

    let (genes, unmatched) = load_amr_catalog(&catalog, Some(&tablef), true).unwrap();
    t.check(genes.len() == 12 && unmatched == 0, format!("catalog: {} genes (the fusion as one), {unmatched} without a line", genes.len()));
    let fused = genes.iter().find(|g| g.protein == "WP_000010.1").unwrap();
    t.check(fused.symbol == "aacF/aphF" && fused.class == "AMINOGLYCOSIDE" && fused.subclass == "AMIKACIN/GENTAMICIN", format!("fusion: {fused:?}"));
    let b = genes.iter().find(|g| g.symbol == "aacB").unwrap();
    t.check(b.subtype == "AMR" && b.subclass == "GENTAMICIN" && b.seq == gb.as_bytes(), "the gene line, not the point-mutation line, upper case sequence");
    let f = genes.iter().find(|g| g.symbol == "fosN").unwrap();
    t.check(f.class == "FOSFOMYCIN" && f.product == "FosN fosfomycin resistance", "joined by the GenBank accession");
    let (core_only, _) = load_amr_catalog(&catalog, Some(&tablef), false).unwrap();
    t.check(core_only.len() == 11 && core_only.iter().all(|g| g.scope == "core"), "include_plus = false drops the plus genes");
    let (bare, _) = load_amr_catalog(&catalog, None, true).unwrap();
    t.check(bare.len() == 12 && bare[5].symbol == "aacB" && bare[5].class.is_empty() && bare[5].product == "aminoglycoside acetyltransferase", "without the table: the product of the header");

    // The contigs.
    let mut c1 = random_dna(60_000, s).into_bytes();
    let mut plants: Vec<Plant> = Vec::new();
    let put = |c: &mut Vec<u8>, at: usize, seq: &str, minus: bool| {
        let s2 = if minus { rc_text(seq) } else { seq.to_string() };
        c.splice(at..at + s2.len(), s2.bytes());
        (at + 1, at + s2.len())
    };
    // blaXYZ-2 exact, plus strand
    let (a, z) = put(&mut c1, 5000, &x2, false);
    plants.push(Plant { gene: "blaXYZ-2", contig: "contig1", start: a, end: z, minus: false, identity: 1.0, coverage: 1.0, status: "exact" });
    // aacB with 36 substitutions (97%), minus strand
    let (a, z) = put(&mut c1, 15_000, &substitute(&gb, 36, 30, s), true);
    plants.push(Plant { gene: "aacB", contig: "contig1", start: a, end: z, minus: true, identity: 1.0 - 36.0 / 1200.0, coverage: 1.0, status: "aligned" });
    // tetC with 4 substitutions, a 3-base deletion and a 2-base insertion,
    // apart from each other
    let mut tc = gc.clone().into_bytes();
    for p in [100, 200, 400, 520] {
        tc[p] = other_base(tc[p], s);
    }
    tc.drain(300..303);
    tc.splice(450..450, *b"GA");
    let tc = String::from_utf8(tc).unwrap();
    let (a, z) = put(&mut c1, 25_000, &tc, false);
    plants.push(Plant { gene: "tetC", contig: "contig1", start: a, end: z, minus: false, identity: (600.0 - 3.0 - 4.0) / 602.0, coverage: 1.0, status: "aligned" });
    // qnrD at 85%: under the threshold, never reported
    put(&mut c1, 35_000, &substitute(&gd, 225, 20, s), false);
    // ermE exact, minus strand
    let (a, z) = put(&mut c1, 45_000, &ge, true);
    plants.push(Plant { gene: "ermE", contig: "contig1", start: a, end: z, minus: true, identity: 1.0, coverage: 1.0, status: "exact" });
    // virP (plus) exact
    let (a, z) = put(&mut c1, 50_000, &gp, false);
    plants.push(Plant { gene: "virP", contig: "contig1", start: a, end: z, minus: false, identity: 1.0, coverage: 1.0, status: "exact" });
    // fosN: the catalog has N and R at 200-201; the contig has C and G
    let mut fo = gn.clone().into_bytes();
    fo[200] = b'C';
    fo[201] = b'G';
    let (a, z) = put(&mut c1, 54_000, &String::from_utf8(fo).unwrap(), false);
    plants.push(Plant { gene: "fosN", contig: "contig1", start: a, end: z, minus: false, identity: 1.0, coverage: 1.0, status: "exact" });

    // contig2: the last 60% of qnrD at its very start (cut by the contig
    // end), blaXYZ-3 exact (the second allele of the family, another locus)
    // and the fusion gene on the minus strand.
    let mut c2 = random_dna(20_000, s).into_bytes();
    let tail = &gd[600..];
    c2.splice(0..tail.len(), tail.bytes());
    plants.push(Plant { gene: "qnrD", contig: "contig2", start: 1, end: tail.len(), minus: false, identity: 1.0, coverage: 0.6, status: "partial_contig_end" });
    let (a, z) = put(&mut c2, 8000, &x3, false);
    plants.push(Plant { gene: "blaXYZ-3", contig: "contig2", start: a, end: z, minus: false, identity: 1.0, coverage: 1.0, status: "exact" });
    let (a, z) = put(&mut c2, 12_000, &gf, true);
    plants.push(Plant { gene: "aacF/aphF", contig: "contig2", start: a, end: z, minus: true, identity: 1.0, coverage: 1.0, status: "exact" });
    // contig3: the last 70% of aacB on the minus strand at the contig's end
    // (on the minus strand the gene runs from high to low coordinates, so
    // its first 30% would lie past the end), and 40% of tetC at the start
    // (under min_partial_coverage).
    let mut c3 = random_dna(9_000, s).into_bytes();
    let n3 = c3.len();
    c3.splice(n3 - 840..n3, rc_text(&gb[360..]).bytes());
    plants.push(Plant { gene: "aacB", contig: "contig3", start: n3 - 840 + 1, end: n3, minus: true, identity: 1.0, coverage: 0.7, status: "partial_contig_end" });
    let t40 = &gc[360..];
    c3.splice(0..t40.len(), t40.bytes());
    let contigs = vec![
        ("contig1 the chromosome".to_string(), String::from_utf8(c1.clone()).unwrap()),
        ("contig2".to_string(), String::from_utf8(c2.clone()).unwrap()),
        ("contig3".to_string(), String::from_utf8(c3.clone()).unwrap()),
    ];
    let asm = t.write_packed("isolate.fa.gz", &fasta(&contigs, 60));
    let out = t.p("amr.tsv");
    let r = find_resistance_genes(&[asm.clone()], &catalog, Some(&tablef), &AmrOptions { output: Some(out.clone()), ..Default::default() }).unwrap();
    t.check(r.mode == "assembly" && r.sequences == 3 && r.catalog_genes == 12, "assembly mode");
    for p in &plants {
        let h = r.hits.iter().find(|h| h.gene == p.gene && h.contig == p.contig);
        let ok = h.is_some_and(|h| {
            h.start as usize == p.start
                && h.end as usize == p.end
                && h.strand == if p.minus { '-' } else { '+' }
                && near(h.identity, p.identity, 1e-9)
                && near(h.coverage, p.coverage, 1e-9)
                && h.status == p.status
        });
        t.check(ok, format!("{} on {}: want {}..{} {} id {:.4} cov {:.2} {}; got {:?}", p.gene, p.contig, p.start, p.end, if p.minus { '-' } else { '+' }, p.identity, p.coverage, p.status, h));
    }
    t.check(r.hits.len() == plants.len(), format!("nothing else: {:?}", r.hits.iter().map(|h| (&h.gene, &h.contig, h.identity, h.coverage)).collect::<Vec<_>>()));
    t.check(!r.hits.iter().any(|h| h.gene == "qnrD" && h.contig == "contig1"), "qnrD at 85% is not reported");
    t.check(!r.hits.iter().any(|h| h.gene == "tetC" && h.contig == "contig3"), "40% of a gene at a contig end is not reported");
    let x = r.hits.iter().find(|h| h.gene == "blaXYZ-2").unwrap();
    t.check(x.subclass == "CARBAPENEM" && x.family == "blaXYZ" && x.scope == "core" && x.mismatches == 0 && x.gap_bases == 0, "the best allele at the locus, with its class");
    let tc_hit = r.hits.iter().find(|h| h.gene == "tetC").unwrap();
    t.check(tc_hit.mismatches == 4 && tc_hit.gap_bases == 5 && tc_hit.aligned_length == 602, format!("tetC: mismatches {} gaps {} columns {}", tc_hit.mismatches, tc_hit.gap_bases, tc_hit.aligned_length));
    let tsv = read_text(&out);
    t.check(tsv.lines().count() == plants.len() + 1 && tsv.contains("blaXYZ-2\tblaXYZ\tclass A beta-lactamase XYZ-2\tcore\tAMR\tAMR\tBETA-LACTAM\tCARBAPENEM\t1.0000\t1.0000\tcontig1\t"), "the TSV");
    t.check(r.json().contains("\"gene\":\"aacF/aphF\"") && r.resumo().contains("cut by the contig end") && r.resumo().contains("(plus)"), "json and resumo");
    let core = find_resistance_genes(&[asm.clone()], &catalog, Some(&tablef), &AmrOptions { include_plus: false, ..Default::default() }).unwrap();
    t.check(core.hits.len() == plants.len() - 1 && !core.hits.iter().any(|h| h.gene == "virP"), "include_plus = false");
    let loose = find_resistance_genes(&[asm.clone()], &catalog, Some(&tablef), &AmrOptions { min_identity: 0.8, ..Default::default() }).unwrap();
    let q = loose.hits.iter().find(|h| h.gene == "qnrD" && h.contig == "contig1");
    t.check(q.is_some_and(|h| near(h.identity, 0.85, 0.01)), format!("qnrD at min_identity 0.8: {q:?}"));
    let nopart = find_resistance_genes(&[asm.clone()], &catalog, Some(&tablef), &AmrOptions { min_partial_coverage: 1.0, ..Default::default() }).unwrap();
    t.check(nopart.hits.iter().all(|h| h.status != "partial_contig_end") && nopart.hits.len() == plants.len() - 2, "min_partial_coverage = 1 drops the cut genes");

    // Reads of an isolate (40x, 0.1% errors, both strands): contig1, and
    // contig2 without its piece of qnrD. Reads cannot tell two copies of one
    // gene apart, so each gene here has one copy -- except blaXYZ, whose two
    // alleles must give one line.
    let mut c2r = c2.clone();
    c2r.splice(0..900, random_dna(900, s).bytes());
    let read_sources = [String::from_utf8(c1.clone()).unwrap(), String::from_utf8(c2r).unwrap()];
    let mut reads: Vec<(String, String)> = Vec::new();
    for (ci, cs) in read_sources.iter().enumerate() {
        let (a1, a2, _) = simulate_pairs(&[("c", cs.as_str())], &[1.0], cs.len() * 40 / 300, 150, 350, 0.001, s);
        for (x, y) in a1.into_iter().zip(a2) {
            reads.push((format!("{ci}{}", x.0), x.1));
            reads.push((format!("{ci}{}", y.0), y.1));
        }
    }
    let rf = t.write_packed("isolate_reads.fq.gz", &fastq(&reads));
    let rr = find_resistance_genes(&[rf.clone()], &catalog, Some(&tablef), &AmrOptions::default()).unwrap();
    t.check(rr.mode == "reads", "FASTQ means reads");
    let got: Vec<(&str, f64, f64)> = rr.hits.iter().map(|h| (h.gene.as_str(), h.identity, h.coverage)).collect();
    let find = |g: &str| rr.hits.iter().find(|h| h.gene == g);
    t.check(find("aacB").is_some_and(|h| near(h.identity, 1.0 - 36.0 / 1200.0, 0.002) && h.coverage > 0.99), format!("reads: aacB at 97%: {got:?}"));
    // 600 positions covered (3 of them deleted), 1 insertion: 593 / 601
    t.check(find("tetC").is_some_and(|h| near(h.identity, 593.0 / 601.0, 1e-9) && h.coverage == 1.0 && h.gap_bases == 4), format!("reads: tetC with indels: {:?}", find("tetC")));
    for g in ["ermE", "virP", "fosN", "aacF/aphF"] {
        t.check(find(g).is_some_and(|h| h.identity == 1.0 && h.coverage == 1.0 && h.depth.is_some_and(|d| d > 20.0)), format!("reads: {g} exact: {:?}", find(g)));
    }
    // Two alleles of one family at two loci: the reads disagree where they
    // differ (mixed), and both are reported; blaXYZ-1, -4 and -5 (one or two
    // substitutions from -2) are not.
    t.check(
        rr.hits.iter().filter(|h| h.family == "blaXYZ").count() == 2
            && find("blaXYZ-2").is_some_and(|h| h.status == "exact")
            && find("blaXYZ-3").is_some_and(|h| h.status == "exact")
            && rr.hits.iter().any(|h| h.family == "blaXYZ" && h.mixed),
        format!("reads: both alleles of the family blaXYZ, flagged mixed: {:?}", rr.hits.iter().filter(|h| h.family == "blaXYZ").map(|h| (&h.gene, h.mixed)).collect::<Vec<_>>()),
    );
    t.check(rr.hits.iter().filter(|h| h.family != "blaXYZ").all(|h| !h.mixed), "reads: one copy, not mixed");
    t.check(find("qnrD").is_none(), "reads: qnrD at 85% and its 60% piece are not reported");
    t.check(rr.hits.len() == 8, format!("reads: eight genes: {got:?}"));
    // A novel allele: -1 with the substitution of -4 and the one of -5. Each
    // of -4 and -5 matches some reads best, but the reads agree with one
    // another everywhere: one line, the closest allele at 899/900, not mixed.
    let novel_reads: Vec<(String, String)> = {
        let src = format!("{}{}{}", random_dna(400, s), novel, random_dna(400, s));
        let (a1, a2, _) = simulate_pairs(&[("n", src.as_str())], &[1.0], 400, 150, 350, 0.001, s);
        a1.into_iter().chain(a2).collect()
    };
    let fnv = t.write("novel_reads.fq", &fastq(&novel_reads));
    let rn = find_resistance_genes(&[fnv], &catalog, Some(&tablef), &AmrOptions::default()).unwrap();
    t.check(
        rn.hits.len() == 1 && ["blaXYZ-4", "blaXYZ-5"].contains(&rn.hits[0].gene.as_str()) && near(rn.hits[0].identity, 899.0 / 900.0, 1e-9) && !rn.hits[0].mixed,
        format!("reads: a novel allele gives one line: {:?}", rn.hits.iter().map(|h| (&h.gene, h.identity, h.mixed)).collect::<Vec<_>>()),
    );
    let novel_asm = t.write("novel_asm.fa", &fasta(&[("c".to_string(), format!("{}{}{}", random_dna(500, s), novel, random_dna(500, s)))], 70));
    let an = find_resistance_genes(&[novel_asm], &catalog, Some(&tablef), &AmrOptions::default()).unwrap();
    t.check(
        an.hits.len() == 1 && ["blaXYZ-4", "blaXYZ-5"].contains(&an.hits[0].gene.as_str()) && near(an.hits[0].identity, 899.0 / 900.0, 1e-9),
        format!("assembly: the same novel allele: {:?}", an.hits.iter().map(|h| (&h.gene, h.identity)).collect::<Vec<_>>()),
    );
    // The reads of the gene alone, one allele: the exact allele wins.
    let only_x1: Vec<(String, String)> = {
        let src = format!("{}{}{}", random_dna(400, s), x1, random_dna(400, s));
        let (a1, a2, _) = simulate_pairs(&[("x", src.as_str())], &[1.0], 400, 150, 350, 0.0, s);
        a1.into_iter().chain(a2).collect()
    };
    let fx = t.write("x1_reads.fq", &fastq(&only_x1));
    let rx = find_resistance_genes(&[fx], &catalog, Some(&tablef), &AmrOptions::default()).unwrap();
    t.check(rx.hits.len() == 1 && rx.hits[0].gene == "blaXYZ-1" && rx.hits[0].status == "exact", format!("reads: blaXYZ-1 beats -2 (two substitutions): {:?}", rx.hits.iter().map(|h| (&h.gene, h.identity)).collect::<Vec<_>>()));

    // Errors.
    t.fails("mixed inputs", find_resistance_genes(&[asm.clone(), rf.clone()], &catalog, None, &AmrOptions::default()), &["mix", "FASTA"]);
    let badcat = t.write("bad_cat.fa", ">gene1 some description\nACGTACGT\n");
    t.fails("not an AMR_CDS header", find_resistance_genes(&[asm.clone()], &badcat, None, &AmrOptions::default()), &["bad_cat.fa", "line 1", "protein|nucleotide", "AMR_CDS.fa"]);
    let badtab = t.write("bad_table.txt", "allele\tgene_family\tproduct_name\n");
    t.fails("a table without its columns", find_resistance_genes(&[asm.clone()], &catalog, Some(&badtab), &AmrOptions::default()), &["scope", "ReferenceGeneCatalog.txt"]);
    t.fails("a FASTQ as the catalog", find_resistance_genes(&[asm.clone()], &rf, None, &AmrOptions::default()), &["FASTQ", "AMR_CDS.fa"]);
    t.fails("identity 90", find_resistance_genes(&[asm.clone()], &catalog, None, &AmrOptions { min_identity: 90.0, ..Default::default() }), &["min_identity", "fraction"]);
    t.fails("no input", find_resistance_genes(&[], &catalog, None, &AmrOptions::default()), &["No input"]);
    t.fails("output over the catalog", find_resistance_genes(&[asm.clone()], &catalog, None, &AmrOptions { output: Some(catalog.clone()), ..Default::default() }), &["also one of the input files"]);
    let empty = t.write("empty.fa", "");
    let re = find_resistance_genes(&[empty], &catalog, None, &AmrOptions::default()).unwrap();
    t.check(re.hits.is_empty() && re.sequences == 0 && re.resumo().contains("no gene"), "an empty assembly is not an error");
}

fn trial_reader(t: &mut T) {
    // Malformed sequence files, through the public functions.
    let db = t.p("community100.sig");
    let g = |p: &str| gather(&[p.to_string()], &[db.clone()], &GatherOptions::default());
    let wrapped = t.write("wrapped.fq", "@r1\nACGT\nACGT\n+\nIIIIIIII\n");
    t.fails("wrapped FASTQ", g(&wrapped), &["line 3", "'+'", "seqtk"]);
    let cut = t.write("cut.fq", "@r1\nACGT\n+\nIIII\n@r2\nACG");
    t.fails("truncated FASTQ", g(&cut), &["truncated", "line 5"]);
    let stray = t.write("stray.fq", "@r1\nACGT\n+\nIIII\nACGT\n+\nIIII\n");
    t.fails("stray line", g(&stray), &["line 5", "'@'"]);
    let table = t.write("table.txt", "chrom\tstart\tend\n");
    t.fails("not sequences", g(&table), &["not a FASTA nor a FASTQ", "chrom"]);
    t.fails("missing", g(&t.p("absent.fq")), &["absent.fq", "does not exist"]);
    // Accepted: CRLF, blank lines, a byte-order mark, no final line break,
    // lower case, comments.
    let ok = t.write("ok.fa", "\u{feff}\r\n>s1 x\r\n;comment\r\nacgtACGT\r\n\r\nGGCC\r\n>s2\r\nTTTT");
    let sk = sketch_sample(&[ok], 4, u64::MAX).unwrap();
    let want = oracle_counts(&["acgtACGTGGCC", "TTTT"], 4, u64::MAX);
    let got: HashMap<u64, u32> = sk.hashes.iter().copied().zip(sk.abunds.iter().copied()).collect();
    t.check(got == want && sk.sequences == 2 && sk.bases == 16 && !sk.reads, "a messy FASTA is read as the sequences it holds");
    let okq = t.write("ok.fq", "@r1\r\nACGTTT\r\n+r1\r\n@@@@@@\r\n\r\n@r2\nAC\n+\nII");
    let sk = sketch_sample(&[okq], 3, u64::MAX).unwrap();
    t.check(sk.sequences == 2 && sk.bases == 8 && sk.reads && sk.kmers == 4, "a FASTQ with CRLF, a quality line starting with '@' and no final break");
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(String::as_str) == Some("real") {
        let (data, work) = (args.get(2).expect("research_data folder"), args.get(3).expect("work folder"));
        real::run(data, work);
        return;
    }
    if args.get(1).map(String::as_str) == Some("gather") {
        // ad hoc: microbio gather <database> <file> [file ...]
        let g = gather(&args[3..], &[args[2].clone()], &GatherOptions::default()).unwrap();
        print!("{}", g.resumo());
        for m in &g.matches {
            println!("  {}: coverage {:?}", m.name, m.kmer_coverage);
        }
        return;
    }
    if args.get(1).map(String::as_str) == Some("classify") {
        // ad hoc: microbio classify <database> <R1> [R2] [assignments.tsv]
        let opt = ClassifyOptions { assignments: args.get(5).cloned(), ..Default::default() };
        let c = classify_reads(&args[3], args.get(4).map(String::as_str), &[args[2].clone()], &opt).unwrap();
        print!("{}", c.resumo());
        return;
    }
    if args.get(1).map(String::as_str) == Some("amr") {
        // ad hoc: microbio amr <AMR_CDS.fa> <ReferenceGeneCatalog.txt> <input> [input ...]
        let r = find_resistance_genes(&args[4..], &args[2], Some(&args[3]), &AmrOptions::default()).unwrap();
        print!("{}", r.resumo());
        return;
    }
    if args.get(1).map(String::as_str) == Some("real-amr-reads") {
        let (data, work) = (args.get(2).expect("research_data folder"), args.get(3).expect("work folder"));
        real::amr_reads(data, work);
        return;
    }
    let t0 = std::time::Instant::now();
    let dir = std::env::temp_dir().join(format!("trial_microbio_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut t = T { ok: 0, total: 0, failed: Vec::new(), dir: dir.clone() };
    let mut s = 20261114u64;
    println!("hashing");
    trial_hash(&mut t);
    println!("sketches and the database");
    trial_sketch(&mut t, &mut s);
    let c = community(&t, &mut s);
    println!("gather");
    trial_gather(&mut t, &c, &mut s);
    println!("read classification");
    trial_classify(&mut t, &c, &mut s);
    println!("resistance genes");
    trial_amr(&mut t, &mut s);
    println!("reading");
    trial_reader(&mut t);
    let _ = std::fs::remove_dir_all(&dir);
    println!("{} of {} checks passed ({:.1} s)", t.ok, t.total, t0.elapsed().as_secs_f64());
    if t.ok != t.total {
        std::process::exit(1);
    }
}

/// The real data of the research_data folder (see its MANIFEST.md): the
/// two Zymo mock communities against the reference genomes with the human
/// decoy, and the resistance genes of the outbreak genomes.
mod real {
    use crate::microbio::*;
    use std::collections::{BTreeMap, HashMap};
    use std::time::Instant;

    fn tsv(path: &str) -> Vec<HashMap<String, String>> {
        let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let mut lines = text.lines();
        let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
        lines
            .filter(|l| !l.trim().is_empty())
            .map(|l| head.iter().zip(l.split('\t')).map(|(h, v)| (h.to_string(), v.to_string())).collect())
            .collect()
    }

    /// "Genus species" of an organism name.
    fn species(org: &str) -> String {
        let w: Vec<&str> = org.split_whitespace().collect();
        if w.len() >= 2 && w[0].chars().next().is_some_and(char::is_uppercase) && w[1].chars().all(|c| c.is_ascii_lowercase()) {
            format!("{} {}", w[0], w[1])
        } else {
            org.to_string()
        }
    }

    fn secs(t: Instant) -> String {
        format!("{:.1} s", t.elapsed().as_secs_f64())
    }

    fn pc(x: f64) -> String {
        if x == 0.0 {
            "0".to_string()
        } else if x >= 0.1 {
            format!("{:.2}", 100.0 * x)
        } else {
            format!("{:.5}", 100.0 * x)
        }
    }

    pub fn run(data: &str, work: &str) {
        std::fs::create_dir_all(work).unwrap();
        let p = |rel: &str| format!("{data}/{rel}");
        let w = |name: &str| format!("{work}/{name}");

        // ------------------------------------------------ the references
        let refs = tsv(&p("reference_genomes/genomes.tsv"));
        let mut seen: HashMap<String, usize> = HashMap::new();
        let mut genomes: Vec<GenomeInput> = Vec::new();
        let mut species_of: HashMap<String, String> = HashMap::new();
        for r in &refs {
            let org = r["organism"].clone();
            *seen.entry(org.clone()).or_insert(0) += 1;
            let name = if refs.iter().filter(|x| x["organism"] == org).count() > 1 { format!("{org} {}", r["strain_or_isolate"]) } else { org.clone() };
            species_of.insert(name.clone(), species(&org));
            genomes.push(GenomeInput { files: vec![p(&r["fasta"])], name, host: false });
        }
        let host_name = "Homo sapiens GRCh38 chr21 + chrM (host decoy)".to_string();
        species_of.insert(host_name.clone(), "host".to_string());
        genomes.push(GenomeInput {
            files: vec![p("metagenome/host_decoy/GRCh38.p14_chr21.fna.gz"), p("metagenome/host_decoy/GRCh38.p14_chrMT.fna.gz")],
            name: host_name,
            host: true,
        });
        let db1000 = w("refs_k31_s1000.sig");
        let db10 = w("refs_k31_s10.sig");
        let t = Instant::now();
        let s1 = sketch_references(&genomes, &db1000, &SketchOptions { k: 31, scaled: 1000 }).unwrap();
        println!("== sketch_references, k=31 scaled=1000: {} genomes in {}, {} bytes", s1.genomes.len(), secs(t), std::fs::metadata(&db1000).unwrap().len());
        print!("{}", s1.resumo());
        let t = Instant::now();
        let s10 = sketch_references(&genomes, &db10, &SketchOptions { k: 31, scaled: 10 }).unwrap();
        println!(
            "== sketch_references, k=31 scaled=10: {} hashes in {}, {} bytes",
            s10.genomes.iter().map(|g| g.n_hashes).sum::<u64>(),
            secs(t),
            std::fs::metadata(&db10).unwrap().len()
        );
        let a = load_databases(&[db10.clone()], Some(1000)).unwrap();
        let b = read_sketch_db(&db1000).unwrap();
        println!(
            "   scaled=10 downsampled to 1000 = sketched at 1000: {}",
            a.sketches.iter().zip(&b.sketches).all(|(x, y)| x.hashes == y.hashes)
        );

        // ------------------------------------------------ the two mock communities
        // Zymo's theoretical composition by GENOMIC DNA (D6300 EVEN: 12% each
        // bacterium, 2% each yeast; D6310 LOG: 89.1% down to 0.000089%).
        let even: Vec<(&str, f64)> = vec![
            ("Pseudomonas aeruginosa", 12.0),
            ("Escherichia coli", 12.0),
            ("Salmonella enterica", 12.0),
            ("Limosilactobacillus fermentum", 12.0),
            ("Enterococcus faecalis", 12.0),
            ("Staphylococcus aureus", 12.0),
            ("Listeria monocytogenes", 12.0),
            ("Bacillus subtilis", 12.0),
            ("Saccharomyces cerevisiae", 2.0),
            ("Cryptococcus neoformans", 2.0),
        ];
        let log: Vec<(&str, f64)> = vec![
            ("Listeria monocytogenes", 89.1),
            ("Pseudomonas aeruginosa", 8.9),
            ("Bacillus subtilis", 0.89),
            ("Saccharomyces cerevisiae", 0.89),
            ("Escherichia coli", 0.089),
            ("Salmonella enterica", 0.089),
            ("Limosilactobacillus fermentum", 0.0089),
            ("Enterococcus faecalis", 0.00089),
            ("Cryptococcus neoformans", 0.00089),
            ("Staphylococcus aureus", 0.000089),
        ];
        let mock = "metagenome/mock_community";
        let runs = [
            ("EVEN", p(&format!("{mock}/ERR2984773_zymo_even_D6300_first1M_R1.fastq.gz")), p(&format!("{mock}/ERR2984773_zymo_even_D6300_first1M_R2.fastq.gz")), even),
            ("LOG", p(&format!("{mock}/ERR2935805_zymo_log_D6310_first1M_R1.fastq.gz")), p(&format!("{mock}/ERR2935805_zymo_log_D6310_first1M_R2.fastq.gz")), log),
        ];
        for (label, r1, r2, expected) in &runs {
            println!("\n================ Zymo {label}");
            let t = Instant::now();
            let g = gather(&[r1.clone(), r2.clone()], &[db1000.clone()], &GatherOptions { output: Some(w(&format!("{label}_gather.tsv"))), ..Default::default() }).unwrap();
            println!("== gather ({})", secs(t));
            print!("{}", g.resumo());
            let t = Instant::now();
            let c = classify_reads(
                r1,
                Some(r2),
                &[db10.clone()],
                &ClassifyOptions {
                    output: Some(w(&format!("{label}_classify.tsv"))),
                    host_free: Some(w(&format!("{label}_hostfree_R1.fq.gz"))),
                    host_free2: Some(w(&format!("{label}_hostfree_R2.fq.gz"))),
                    ..Default::default()
                },
            )
            .unwrap();
            println!("== classify_reads, all references ({})", secs(t));
            print!("{}", c.resumo());
            // Restricted to what gather found.
            let only: Vec<String> = g.matches.iter().map(|m| m.name.clone()).collect();
            let t = Instant::now();
            let co = classify_reads(r1, Some(r2), &[db10.clone()], &ClassifyOptions { only: only.clone(), output: Some(w(&format!("{label}_classify_only.tsv"))), ..Default::default() }).unwrap();
            println!("== classify_reads, only the {} references gather found ({})", only.len(), secs(t));
            print!("{}", co.resumo());

            // By species.
            let mut by_g: BTreeMap<String, (f64, Vec<String>)> = BTreeMap::new();
            for m in g.matches.iter().filter(|m| !m.host) {
                let e = by_g.entry(species_of[&m.name].clone()).or_default();
                e.0 += m.relative_abundance.unwrap_or(0.0);
                e.1.push(format!("{} ANI {:.4}/{}", m.name, m.ani, m.ani_adjusted.map_or("-".to_string(), |x| format!("{x:.4}"))));
            }
            let mut by_c: BTreeMap<String, f64> = BTreeMap::new();
            for r in c.rows.iter().filter(|r| !r.host) {
                *by_c.entry(species_of[&r.name].clone()).or_default() += r.fraction_of_classified.unwrap_or(0.0);
            }
            let mut by_co: BTreeMap<String, f64> = BTreeMap::new();
            for r in co.rows.iter().filter(|r| !r.host) {
                *by_co.entry(species_of[&r.name].clone()).or_default() += r.fraction_of_classified.unwrap_or(0.0);
            }
            println!("\n   {label}: % by species (expected = Zymo, genomic DNA)");
            println!("   {:<32} {:>9} {:>9} {:>11} {:>11}", "species", "expected", "gather", "classify", "classify*");
            let mut names: Vec<String> = expected.iter().map(|e| e.0.to_string()).collect();
            for k in by_g.keys().chain(by_c.keys()) {
                if !names.contains(k) {
                    names.push(k.clone());
                }
            }
            for n in &names {
                let e = expected.iter().find(|x| x.0 == n).map(|x| x.1);
                let gg = by_g.get(n).map(|x| x.0).unwrap_or(0.0);
                let cc = by_c.get(n).copied().unwrap_or(0.0);
                let oo = by_co.get(n).copied().unwrap_or(0.0);
                if e.is_none() && cc < 0.0005 && gg == 0.0 {
                    continue;
                }
                println!("   {:<32} {:>9} {:>9} {:>11} {:>11}", n, e.map_or("-".to_string(), |x| format!("{x}")), pc(gg), pc(cc), pc(oo));
            }
            let other: f64 = by_c.iter().filter(|(k, _)| !names.contains(k)).map(|x| x.1).sum();
            println!("   (classify* = only the references gather reported; other species in 'classify': {})", pc(other));
            for (k, v) in &by_g {
                println!("   gather {k}: {}", v.1.join("; "));
            }
        }

        // ------------------------------------------------ resistance genes of the outbreak genomes
        let cat = p("amr/amrfinderplus/2026-10-29.1/AMR_CDS.fa");
        let table = p("amr/amrfinderplus/2026-10-29.1/ReferenceGeneCatalog.txt");
        let (genes, unmatched) = load_amr_catalog(&cat, Some(&table), true).unwrap();
        println!("\n================ resistance genes: catalog of {} genes ({} without a line in the table)", genes.len(), unmatched);
        let family_of: HashMap<String, String> = genes.iter().map(|g| (g.symbol.clone(), g.family.clone())).collect();
        let outbreak = tsv(&p("outbreak/bacteria/genomes.tsv"));
        let exact = tsv(&p("amr/screening_of_downloaded_genomes/exact_allele_hits.tsv"));
        let mut tally = (0usize, 0usize, 0usize, 0usize); // key genes: exact, close allele, missing, total
        let mut exact_tally = (0usize, 0usize, 0usize); // exact screen: found, missed, extra
        for row in &outbreak {
            let acc = &row["assembly_accession"];
            let t = Instant::now();
            let r = find_resistance_genes(&[p(&row["fasta"])], &cat, Some(&table), &AmrOptions { output: Some(w(&format!("{acc}_amr.tsv"))), ..Default::default() }).unwrap();
            println!("\n== {} {} (ST{}), {} hits in {}", acc, row["strain"], row["MLST_ST"], r.hits.len(), secs(t));
            let core: Vec<String> = r
                .hits
                .iter()
                .filter(|h| h.scope == "core")
                .map(|h| if h.status == "exact" { h.gene.clone() } else { format!("{}({:.3}/{:.2})", h.gene, h.identity, h.coverage) })
                .collect();
            println!("   core: {}", core.join(" "));
            let key = &row["key_AMR_genes_best_allele"];
            if !key.starts_with("none") {
                for k in key.split(';') {
                    let (allele, approx) = match k.split_once("(~") {
                        Some((a, _)) => (a.to_string(), true),
                        None => (k.to_string(), false),
                    };
                    tally.3 += 1;
                    let fam = family_of.get(&allele).cloned().unwrap_or_default();
                    let same = r.hits.iter().find(|h| h.gene == allele && h.status == "exact");
                    let close = r.hits.iter().filter(|h| h.family == fam).max_by(|a, b| a.identity.total_cmp(&b.identity));
                    let verdict = match (same, close) {
                        (Some(_), _) => {
                            tally.0 += 1;
                            "same allele, exact".to_string()
                        }
                        (None, Some(h)) => {
                            tally.1 += 1;
                            format!("{} at identity {:.4}, coverage {:.3} ({})", h.gene, h.identity, h.coverage, h.status)
                        }
                        _ => {
                            tally.2 += 1;
                            "MISSING".to_string()
                        }
                    };
                    println!("   key {k:<26} -> {verdict}{}", if approx { "  [screen: no exact allele]" } else { "" });
                }
            } else {
                println!("   key genes: none expected; core hits: {}", core.len());
            }
            // The exact-match screen of the download (a multiset of alleles).
            let mut want: BTreeMap<String, i32> = BTreeMap::new();
            for e in exact.iter().filter(|e| &e["assembly"] == acc) {
                *want.entry(e["allele_symbol"].clone()).or_insert(0) += 1;
            }
            let mut got: BTreeMap<String, i32> = BTreeMap::new();
            for h in r.hits.iter().filter(|h| h.status == "exact") {
                for part in h.gene.split('/') {
                    *got.entry(part.to_string()).or_insert(0) += 1;
                }
            }
            let mut missed = Vec::new();
            let mut extra = Vec::new();
            for (a, &n) in &want {
                let g = got.get(a).copied().unwrap_or(0);
                exact_tally.0 += n.min(g) as usize;
                if g < n {
                    missed.push(format!("{a} x{}", n - g));
                }
            }
            for (a, &n) in &got {
                let e = want.get(a).copied().unwrap_or(0);
                if n > e {
                    extra.push(format!("{a} x{}", n - e));
                }
            }
            exact_tally.1 += missed.len();
            exact_tally.2 += extra.len();
            println!(
                "   exact screen: {} alleles expected; missed: [{}]; extra exact hits: [{}]",
                want.values().sum::<i32>(),
                missed.join(", "),
                extra.join(", ")
            );
            for m in &missed {
                let a = m.split(" x").next().unwrap();
                let fam = family_of.get(a).cloned().unwrap_or_default();
                let near: Vec<String> = r.hits.iter().filter(|h| h.family == fam || h.gene == a).map(|h| format!("{} {:.4}/{:.3} {}:{}", h.gene, h.identity, h.coverage, h.contig, h.start)).collect();
                println!("      {a}: we report {}", if near.is_empty() { "nothing of its family".to_string() } else { near.join("; ") });
            }
        }
        println!(
            "\n== key genes: {} total, {} same allele exact, {} same family other allele, {} missing",
            tally.3, tally.0, tally.1, tally.2
        );
        println!("== exact screen: {} alleles found exact; {} allele lines missed; {} extra lines", exact_tally.0, exact_tally.1, exact_tally.2);
    }

    /// The reads mode of the resistance genes: reads simulated from real
    /// outbreak genomes (compared with the assembly mode on the same genome),
    /// and the real reads of the Zymo EVEN community.
    pub fn amr_reads(data: &str, work: &str) {
        std::fs::create_dir_all(work).unwrap();
        let p = |rel: &str| format!("{data}/{rel}");
        let w = |name: &str| format!("{work}/{name}");
        let cat = p("amr/amrfinderplus/2026-10-29.1/AMR_CDS.fa");
        let table = p("amr/amrfinderplus/2026-10-29.1/ReferenceGeneCatalog.txt");
        let mut s = 7u64;
        for (acc, fasta) in [
            ("KPNIH1", "outbreak/bacteria/klebsiella_pneumoniae/GCF_000281535.2_ASM28153v2_genomic.fna.gz"),
            ("BAA-2146", "outbreak/bacteria/klebsiella_pneumoniae/GCF_000364385.3_ASM36438v3_genomic.fna.gz"),
            ("MRSA252", "outbreak/bacteria/staphylococcus_aureus_mrsa/GCF_000011505.1_ASM1150v1_genomic.fna.gz"),
            // chromosomal blaSHV and plasmid blaSHV-12: two loci of one family
            ("MGH78578", "outbreak/bacteria/klebsiella_pneumoniae/GCF_000016305.1_ASM1630v1_genomic.fna.gz"),
        ] {
            let mut rd = SeqReader::open(&p(fasta)).unwrap();
            let mut rec = SeqRecord::default();
            let mut seqs: Vec<String> = Vec::new();
            while rd.next(&mut rec).unwrap() {
                seqs.push(String::from_utf8(rec.seq.to_ascii_uppercase()).unwrap());
            }
            // 30x of 2 x 150, 0.3% substitutions, every sequence (plasmids too)
            let fq = w(&format!("{acc}_sim_reads.fq.gz"));
            let mut out: Vec<(String, String)> = Vec::new();
            for (i, q) in seqs.iter().enumerate() {
                if q.len() < 500 {
                    continue;
                }
                let n = q.len() * 30 / 300;
                let (a, b, _) = crate::simulate_pairs(&[("s", q.as_str())], &[1.0], n, 150, 400, 0.003, &mut s);
                for (x, y) in a.into_iter().zip(b) {
                    out.push((format!("{i}_{}", x.0), x.1));
                    out.push((format!("{i}_{}", y.0), y.1));
                }
            }
            let mut wtr = crate::externo::Escritor::cria(&fq).unwrap();
            std::io::Write::write_all(&mut wtr, crate::fastq(&out).as_bytes()).unwrap();
            wtr.termina().unwrap();
            let t = Instant::now();
            let asm = find_resistance_genes(&[p(fasta)], &cat, Some(&table), &AmrOptions::default()).unwrap();
            let ta = secs(t);
            let t = Instant::now();
            let rr = find_resistance_genes(&[fq], &cat, Some(&table), &AmrOptions { output: Some(w(&format!("{acc}_reads_amr.tsv"))), ..Default::default() }).unwrap();
            println!("\n== {acc}: assembly mode {} hits ({ta}); reads mode on {} simulated reads, {} hits ({})", asm.hits.len(), out.len(), rr.hits.len(), secs(t));
            // Each gene of the assembly (several copies count once), and what
            // the reads say of it.
            let show = |h: &AmrHit| {
                format!(
                    "{} {:.4}/{:.2}{}{}",
                    h.gene,
                    h.identity,
                    h.coverage,
                    h.depth.map_or(String::new(), |d| format!(" depth {d:.0}")),
                    if h.mixed { " MIXED" } else { "" }
                )
            };
            let mut asm_genes: Vec<&AmrHit> = Vec::new();
            for h in &asm.hits {
                if !asm_genes.iter().any(|x| x.gene == h.gene) {
                    asm_genes.push(h);
                }
            }
            let (mut same, mut family, mut absent) = (0, Vec::new(), Vec::new());
            for h in &asm_genes {
                if rr.hits.iter().any(|r| r.gene == h.gene) {
                    same += 1;
                } else if let Some(r) = rr.hits.iter().find(|r| r.family == h.family) {
                    family.push(format!("{} -> reads {}", show(h), show(r)));
                } else {
                    absent.push(show(h));
                }
            }
            let extra: Vec<String> =
                rr.hits.iter().filter(|r| !asm_genes.iter().any(|h| h.gene == r.gene || h.family == r.family)).map(show).collect();
            println!(
                "   assembly genes: {}; the reads give the same allele for {same}, another allele of the family for {}, nothing for {}; reads-only genes: {}",
                asm_genes.len(),
                family.len(),
                absent.len(),
                extra.len()
            );
            for x in family {
                println!("     other allele  {x}");
            }
            for x in absent {
                println!("     not in reads  {x}");
            }
            for x in extra {
                println!("     reads only    {x}");
            }
            for r in rr.hits.iter().filter(|r| r.mixed) {
                println!("     mixed         {}", show(r));
            }
        }
        // The real reads of the EVEN community.
        let mock = "metagenome/mock_community";
        let (r1, r2) = (
            p(&format!("{mock}/ERR2984773_zymo_even_D6300_first1M_R1.fastq.gz")),
            p(&format!("{mock}/ERR2984773_zymo_even_D6300_first1M_R2.fastq.gz")),
        );
        let t = Instant::now();
        let rr = find_resistance_genes(&[r1, r2], &cat, Some(&table), &AmrOptions { output: Some(w("EVEN_reads_amr.tsv")), ..Default::default() }).unwrap();
        println!("\n== Zymo EVEN, real reads, reads mode ({}, {} alignments)", secs(t), rr.alignments);
        print!("{}", rr.resumo());
    }
}
