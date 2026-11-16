#![allow(dead_code)]
// Trials of src/contaminantes.rs: the adventitious-agent screen.
//
// Every result is compared with an oracle written here in the plainest way:
// the minimizers by hashing every k-mer of every window as text; the index
// by HashMaps of those; the read calls by recounting each read against the
// HashMaps with the rules written out again; SHA-256 against the FIPS 180-4
// test vectors; the calls and the negative controls on simulated samples
// whose truth is known (reads simulated from planted genomes).
//
// Usage: cargo run --release --bin contaminantes
//        cargo run --release --bin contaminantes -- host-index <research_data> <work dir>
//        cargo run --release --bin contaminantes -- spike <research_data> <mycoplasma dir> <lambda.fasta.gz> <work dir>
//        cargo run --release --bin contaminantes -- zymo <research_data> <mycoplasma dir> <lambda.fasta.gz> <work dir>
//
// `host-index`: the human chr21 (+ chrM) host set: memory, time, host
// depletion on simulated reads, multi-pass build. `spike`: the spike-in
// series (limit of detection, false positives, a relative not in the
// database, an integrated viral fragment, a reagent background). `zymo`:
// the ZymoBIOMICS mock communities, and the same reads through
// microbio::classify_reads.
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
#[path = "../../../src/contaminantes.rs"]
mod contaminantes;

use contaminantes::*;
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

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
    fn path(&self, name: &str) -> String {
        self.dir.join(name).to_string_lossy().to_string()
    }
}

// ---------------------------------------------------------------------
//  Random numbers, sequences, files
// ---------------------------------------------------------------------

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: usize) -> usize {
        ((self.next() >> 11) as f64 / (1u64 << 53) as f64 * n as f64) as usize % n.max(1)
    }
    fn f(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn normal(&mut self) -> f64 {
        let u1 = self.f().max(1e-300);
        let u2 = self.f();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
    fn poisson(&mut self, lambda: f64) -> u64 {
        if lambda <= 0.0 {
            return 0;
        }
        if lambda > 60.0 {
            return (lambda + lambda.sqrt() * self.normal()).round().max(0.0) as u64;
        }
        let l = (-lambda).exp();
        let (mut k, mut p) = (0u64, 1.0);
        loop {
            p *= self.f();
            if p <= l {
                return k;
            }
            k += 1;
        }
    }
}

fn random_dna(n: usize, gc: f64, r: &mut Rng) -> Vec<u8> {
    (0..n)
        .map(|_| {
            let x = r.f();
            if x < gc / 2.0 {
                b'G'
            } else if x < gc {
                b'C'
            } else if x < gc + (1.0 - gc) / 2.0 {
                b'A'
            } else {
                b'T'
            }
        })
        .collect()
}

fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        _ => b'N',
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| comp(b)).collect()
}

fn other_base(b: u8, r: &mut Rng) -> u8 {
    let bases = [b'A', b'C', b'G', b'T'];
    loop {
        let x = bases[r.below(4)];
        if x != b.to_ascii_uppercase() {
            return x;
        }
    }
}

fn mutate(s: &[u8], rate: f64, r: &mut Rng) -> Vec<u8> {
    s.iter().map(|&b| if r.f() < rate { other_base(b, r) } else { b }).collect()
}

fn write_fasta(path: &str, recs: &[(String, Vec<u8>)]) {
    let mut o = Vec::new();
    for (n, s) in recs {
        o.extend_from_slice(format!(">{n}\n").as_bytes());
        for c in s.chunks(80) {
            o.extend_from_slice(c);
            o.push(b'\n');
        }
    }
    std::fs::write(path, o).unwrap();
}

fn read_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
    let mut rd = microbio::SeqReader::open(path).unwrap();
    let mut rec = microbio::SeqRecord::default();
    let mut out = Vec::new();
    while rd.next(&mut rec).unwrap() {
        out.push((String::from_utf8_lossy(&rec.head).to_string(), rec.seq.clone()));
    }
    out
}

/// Illumina-like reads: substitutions whose rate grows along the read
/// (0.1% at the first base, 2% at the last: ~0.73% on average), qualities
/// that say so.
#[derive(Clone, Copy)]
struct Sim {
    len: usize,
    insert: f64,
    sd: f64,
    e0: f64,
    e1: f64,
}

const ILLUMINA: Sim = Sim { len: 150, insert: 350.0, sd: 50.0, e0: 0.001, e1: 0.02 };

struct Genome {
    recs: Vec<Vec<u8>>,
    cum: Vec<u64>,
}

impl Genome {
    fn new(recs: Vec<Vec<u8>>) -> Genome {
        let mut cum = vec![0u64];
        for r in &recs {
            cum.push(cum.last().unwrap() + r.len() as u64);
        }
        Genome { recs, cum }
    }
    fn from_fasta(paths: &[String]) -> Genome {
        Genome::new(paths.iter().flat_map(|p| read_fasta(p).into_iter().map(|x| x.1)).collect())
    }
    fn len(&self) -> u64 {
        *self.cum.last().unwrap()
    }
}

fn add_errors(s: &mut [u8], sim: &Sim, r: &mut Rng) -> Vec<u8> {
    let n = s.len();
    let mut q = Vec::with_capacity(n);
    for (i, b) in s.iter_mut().enumerate() {
        let t = i as f64 / (n.max(2) - 1) as f64;
        let p = sim.e0 + (sim.e1 - sim.e0) * t * t;
        if r.f() < p {
            *b = other_base(*b, r);
        }
        let qv = (-10.0 * p.log10()).round().clamp(2.0, 41.0) as u8;
        q.push(33 + qv);
    }
    q
}

/// One read pair from `g` (a region `(record, start, end)` when given). No N
/// in the fragment.
fn sim_pair(g: &Genome, sim: &Sim, region: Option<(usize, usize, usize)>, r: &mut Rng) -> (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>) {
    for _ in 0..1000 {
        let ins = (sim.insert + sim.sd * r.normal()).round().max(sim.len as f64) as usize;
        let (rec, s, e) = match region {
            Some(x) => x,
            None => {
                let x = (r.f() * g.len() as f64) as u64;
                let i = g.cum.partition_point(|&c| c <= x) - 1;
                (i, 0, g.recs[i].len())
            }
        };
        let span = e - s;
        let ins = ins.min(span);
        if ins < sim.len {
            continue;
        }
        let st = s + r.below(span - ins + 1);
        let frag = &g.recs[rec][st..st + ins];
        if frag.iter().any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')) {
            continue;
        }
        let frag: Vec<u8> = if r.below(2) == 0 { frag.to_vec() } else { revcomp(frag) };
        let mut a: Vec<u8> = frag[..sim.len].to_ascii_uppercase();
        let mut b: Vec<u8> = revcomp(&frag)[..sim.len].to_ascii_uppercase();
        let qa = add_errors(&mut a, sim, r);
        let qb = add_errors(&mut b, sim, r);
        return (a, qa, b, qb);
    }
    panic!("no fragment without N could be drawn");
}

/// FASTQ pair writer (plain, or gzip when the name ends in .gz).
struct PairWriter {
    a: Box<dyn Write>,
    b: Box<dyn Write>,
    n: u64,
    prefix: String,
}

fn open_w(path: &str) -> Box<dyn Write> {
    let f = std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(path).unwrap());
    if path.ends_with(".gz") {
        Box::new(flate2::write::GzEncoder::new(f, flate2::Compression::fast()))
    } else {
        Box::new(f)
    }
}

impl PairWriter {
    fn new(p1: &str, p2: &str, prefix: &str) -> PairWriter {
        PairWriter { a: open_w(p1), b: open_w(p2), n: 0, prefix: prefix.to_string() }
    }
    fn put(&mut self, tag: &str, x: (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>)) {
        self.n += 1;
        let name = format!("{}{}:{}", self.prefix, self.n, tag);
        write!(self.a, "@{name} 1:N:0\n").unwrap();
        self.a.write_all(&x.0).unwrap();
        self.a.write_all(b"\n+\n").unwrap();
        self.a.write_all(&x.1).unwrap();
        self.a.write_all(b"\n").unwrap();
        write!(self.b, "@{name} 2:N:0\n").unwrap();
        self.b.write_all(&x.2).unwrap();
        self.b.write_all(b"\n+\n").unwrap();
        self.b.write_all(&x.3).unwrap();
        self.b.write_all(b"\n").unwrap();
    }
    fn put_n(&mut self, tag: &str, g: &Genome, n: u64, region: Option<(usize, usize, usize)>, r: &mut Rng) {
        for _ in 0..n {
            let x = sim_pair(g, &ILLUMINA, region, r);
            self.put(tag, x);
        }
    }
    fn finish(mut self) {
        self.a.flush().unwrap();
        self.b.flush().unwrap();
    }
}

// ---------------------------------------------------------------------
//  The oracle: k-mers, minimizers and the rules, written out plainly
// ---------------------------------------------------------------------

fn naive_code(s: &[u8]) -> Option<u64> {
    let mut c = 0u64;
    for &b in s {
        let v = match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => return None,
        };
        c = (c << 2) | v;
    }
    Some(c)
}

fn naive_hash(kmer: &[u8]) -> Option<u64> {
    let f = naive_code(kmer)?;
    let r = naive_code(&revcomp(kmer))?;
    Some(kmer_hash(f.min(r)))
}

fn naive_low(kmer: &[u8]) -> bool {
    let mut seen = HashSet::new();
    for t in kmer.windows(3) {
        seen.insert(t.to_ascii_uppercase());
    }
    seen.len() < 8
}

/// Every window minimum (hash, position), each window being w consecutive
/// k-mers of a run of A/C/G/T.
fn naive_minimizers(seq: &[u8], k: usize, w: usize) -> Vec<(u64, usize)> {
    let mut out = HashSet::new();
    let mut i = 0;
    while i < seq.len() {
        if naive_code(&seq[i..i + 1]).is_none() {
            i += 1;
            continue;
        }
        let mut j = i;
        while j < seq.len() && naive_code(&seq[j..j + 1]).is_some() {
            j += 1;
        }
        // run [i, j)
        if j - i >= k {
            let kms: Vec<(u64, usize)> = (i..=j - k).map(|p| (naive_hash(&seq[p..p + k]).unwrap(), p)).collect();
            if kms.len() >= w {
                for s in 0..=kms.len() - w {
                    let mut best = kms[s];
                    for &x in &kms[s..s + w] {
                        if x.0 < best.0 {
                            best = x;
                        }
                    }
                    out.insert(best);
                }
            }
        }
        i = j;
    }
    let mut v: Vec<(u64, usize)> = out.into_iter().collect();
    v.sort();
    v
}

fn naive_read_set(seqs: &[&[u8]], k: usize, w: usize) -> Vec<u64> {
    let mut s = HashSet::new();
    for q in seqs {
        for (h, p) in naive_minimizers(q, k, w) {
            if !naive_low(&q[p..p + k]) {
                s.insert(h);
            }
        }
    }
    let mut v: Vec<u64> = s.into_iter().collect();
    v.sort();
    v
}

struct Oracle {
    k: usize,
    w: usize,
    hw: usize,
    /// hash -> references holding it.
    index: HashMap<u64, Vec<u32>>,
    host: HashSet<u64>,
}

impl Oracle {
    fn build(refs: &[(String, Vec<Vec<u8>>, bool)], k: usize, w: usize, hw: usize) -> Oracle {
        let mut index: HashMap<u64, Vec<u32>> = HashMap::new();
        let mut host = HashSet::new();
        for (i, (_, recs, is_host)) in refs.iter().enumerate() {
            for rec in recs {
                if *is_host {
                    for (h, p) in naive_minimizers(rec, k, hw) {
                        if !naive_low(&rec[p..p + k]) {
                            host.insert(h);
                        }
                    }
                } else {
                    for (h, p) in naive_minimizers(rec, k, w) {
                        if !naive_low(&rec[p..p + k]) {
                            let v = index.entry(h).or_default();
                            if !v.contains(&(i as u32)) {
                                v.push(i as u32);
                            }
                        }
                    }
                }
            }
        }
        Oracle { k, w, hw, index, host }
    }

    /// The call of a read (pair) by the documented rules.
    fn call(&self, mates: &[&[u8]], o: &ScreenOptions) -> Outcome {
        let t = naive_read_set(mates, self.k, self.w);
        let h = if self.host.is_empty() { Vec::new() } else { naive_read_set(mates, self.k, self.hw) };
        if t.is_empty() && h.is_empty() {
            return Outcome::NoKmers;
        }
        let host_hits = h.iter().filter(|x| self.host.contains(x)).count() as u32;
        let mut hits: HashMap<u32, u32> = HashMap::new();
        for m in &t {
            if let Some(v) = self.index.get(m) {
                for &r in v {
                    *hits.entry(r).or_insert(0) += 1;
                }
            }
        }
        let best = hits.values().copied().max().unwrap_or(0);
        let best_ref = hits.iter().filter(|x| *x.1 == best).map(|x| *x.0).min().unwrap_or(0);
        let tfrac = if t.is_empty() { 0.0 } else { best as f64 / t.len() as f64 };
        let hfrac = if h.is_empty() { 0.0 } else { host_hits as f64 / h.len() as f64 };
        if host_hits > 0 && host_hits >= o.host_min_hits && hfrac >= tfrac {
            return Outcome::Host;
        }
        if best == 0 || best < o.min_hits {
            return Outcome::Unclassified;
        }
        if tfrac < o.min_fraction {
            return Outcome::Weak(best_ref);
        }
        let mut set: Vec<u32> = hits.iter().filter(|x| *x.1 == best || *x.1 as f64 >= best as f64 * (1.0 - o.unique_margin)).map(|x| *x.0).collect();
        set.sort();
        if set.len() == 1 {
            Outcome::Unique(set[0])
        } else {
            Outcome::Tie(set)
        }
    }
}

// ---------------------------------------------------------------------
//  1. SHA-256
// ---------------------------------------------------------------------

fn trial_sha(t: &mut T) {
    println!("SHA-256");
    t.check(sha256_hex(b"") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", "sha256 of the empty string");
    t.check(sha256_hex(b"abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", "sha256 of abc");
    t.check(
        sha256_hex(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq")
            == "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1",
        "sha256 of the 448-bit vector",
    );
    let million = vec![b'a'; 1_000_000];
    t.check(sha256_hex(&million) == "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0", "sha256 of a million a");
    // In pieces of every size: the same digest.
    let mut r = Rng(7);
    let data: Vec<u8> = (0..100_003).map(|_| r.next() as u8).collect();
    let whole = sha256_hex(&data);
    let mut ok = true;
    for step in [1usize, 3, 55, 56, 63, 64, 65, 127, 1000, 4096] {
        let mut h = Sha256::new();
        for c in data.chunks(step) {
            h.update(c);
        }
        ok &= hex(&h.finalize()) == whole;
    }
    t.check(ok, "sha256 fed in pieces of 1..4096 bytes equals the whole");
    let p = t.path("sha_test.bin");
    std::fs::write(&p, &data).unwrap();
    let (n, h) = sha256_file(&p).unwrap();
    t.check(n == data.len() as u64 && h == whole, "sha256_file equals the in-memory digest");
    // Padding edge cases: lengths 55, 56, 64 against the one-shot digest of
    // farmacogenetica's independent implementation is not available here;
    // check them against each other's prefix property instead (different
    // lengths, different digests).
    let d55 = sha256_hex(&data[..55]);
    let d56 = sha256_hex(&data[..56]);
    let d64 = sha256_hex(&data[..64]);
    t.check(d55 != d56 && d56 != d64 && d55.len() == 64, "sha256 padding lengths give distinct digests");
}

// ---------------------------------------------------------------------
//  2. Minimizers
// ---------------------------------------------------------------------

fn trial_minimizers(t: &mut T) {
    println!("Minimizers");
    let mut r = Rng(11);
    let mut all_ok = true;
    let mut rc_ok = true;
    let mut cases = 0;
    for case in 0..60 {
        let n = 50 + r.below(3000);
        let mut s = random_dna(n, 0.3 + 0.4 * r.f(), &mut r);
        // N runs, lower case, a microsatellite.
        for _ in 0..r.below(4) {
            let p = r.below(n);
            let l = 1 + r.below(40);
            for b in s.iter_mut().skip(p).take(l) {
                *b = b'N';
            }
        }
        for b in s.iter_mut() {
            if r.below(10) == 0 {
                *b = b.to_ascii_lowercase();
            }
        }
        if case % 5 == 0 && n > 200 {
            for j in 0..60 {
                s[100 + j] = if j % 2 == 0 { b'C' } else { b'A' };
            }
        }
        let (k, w) = [(31, 12), (31, 1), (21, 5), (15, 31), (31, 31)][case % 5];
        let mut got = Vec::new();
        let mut dq = std::collections::VecDeque::new();
        minimizers(&s, k, w, &mut dq, |h, p| got.push((h, p)));
        let mut gs = got.clone();
        gs.sort();
        gs.dedup();
        let want = naive_minimizers(&s, k, w);
        all_ok &= gs == want;
        // Consecutive outputs never repeat a position.
        all_ok &= got.windows(2).all(|x| x[0].1 != x[1].1);
        // The reverse complement holds the same minimizer hashes.
        let rc = revcomp(&s);
        let mut a: Vec<u64> = want.iter().map(|x| x.0).collect();
        let mut b: Vec<u64> = naive_minimizers(&rc, k, w).iter().map(|x| x.0).collect();
        a.sort();
        a.dedup();
        b.sort();
        b.dedup();
        rc_ok &= a == b;
        cases += 1;
    }
    t.check(all_ok, format!("minimizers equal the naive window minima in {cases} random sequences (N, lower case, k 15..31, w 1..31)"));
    t.check(rc_ok, "a sequence and its reverse complement have the same minimizers");
    // Low complexity.
    let micro: Vec<u8> = b"CACACACACACACACACACACACACACACAC".to_vec();
    let homo = vec![b'A'; 31];
    let tri: Vec<u8> = b"CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGC".to_vec();
    let rnd = random_dna(31, 0.5, &mut r);
    t.check(low_complexity(&micro, 0, 31) && low_complexity(&homo, 0, 31) && low_complexity(&tri, 0, 31), "microsatellites and homopolymers are low complexity");
    t.check(!low_complexity(&rnd, 0, 31), "a random 31-mer is not low complexity");
    let mut same = true;
    for _ in 0..2000 {
        let s = random_dna(31, 0.2 + 0.6 * r.f(), &mut r);
        let s = if r.below(3) == 0 { [b"ACACACACACACACAC".as_slice(), &s[..15]].concat() } else { s };
        same &= low_complexity(&s, 0, 31) == naive_low(&s);
    }
    t.check(same, "low_complexity equals the naive count of trinucleotides");
    // The hash is a bijection on codes: different k-mers, different hashes.
    let mut seen = HashSet::new();
    let mut distinct = true;
    for c in 0..200_000u64 {
        distinct &= seen.insert(kmer_hash(c));
    }
    t.check(distinct, "kmer_hash gives 200,000 distinct values for 200,000 codes");
}

// ---------------------------------------------------------------------
//  3. A synthetic community: database, read calls, screen
// ---------------------------------------------------------------------

struct Community {
    host: Genome,
    /// A 7%-divergent relative of virus A that is NOT in the database.
    far: Genome,
    host_files: Vec<String>,
    names: Vec<(String, AgentKind)>,
    genomes: Vec<Genome>,
    files: Vec<String>,
}

const REFS: [(&str, AgentKind); 9] = [
    ("virus A", AgentKind::Virus),
    ("virus A2 (relative of A)", AgentKind::Virus),
    ("virus B", AgentKind::Virus),
    ("bacterium", AgentKind::Bacterium),
    ("kitome bacterium", AgentKind::Bacterium),
    ("mycoplasma", AgentKind::Mycoplasma),
    ("fungus", AgentKind::Fungus),
    ("process control", AgentKind::Spike),
    ("bacterium twin (0.1% divergent)", AgentKind::Bacterium),
];

fn community(t: &T) -> Community {
    let mut r = Rng(2026);
    // A host of 3 Mb with a repeated element and a microsatellite region.
    let element = random_dna(3000, 0.42, &mut r);
    let mut host_recs = Vec::new();
    for c in 0..3 {
        let mut s = random_dna(1_000_000, 0.41, &mut r);
        for j in 0..20 {
            let p = 10_000 + j * 45_000;
            s[p..p + 3000].copy_from_slice(&mutate(&element, 0.05, &mut r));
        }
        for j in 0..400 {
            s[500_000 + j] = if j % 2 == 0 { b'T' } else { b'G' };
        }
        if c == 1 {
            for b in s.iter_mut().skip(700_000).take(5000) {
                *b = b'N';
            }
        }
        host_recs.push(s);
    }
    let host_path = t.path("host.fa");
    write_fasta(&host_path, &host_recs.iter().enumerate().map(|(i, s)| (format!("host_chr{} synthetic host", i + 1), s.clone())).collect::<Vec<_>>());
    let host = Genome::new(host_recs);

    let virus_a = random_dna(12_000, 0.45, &mut r);
    let virus_a2 = mutate(&virus_a, 0.04, &mut r);
    let virus_b = random_dna(7_000, 0.5, &mut r);
    let bact = random_dna(600_000, 0.5, &mut r);
    let plasmid = random_dna(40_000, 0.52, &mut r);
    let kitome = random_dna(300_000, 0.62, &mut r);
    let myco = random_dna(400_000, 0.28, &mut r);
    let fun1 = random_dna(500_000, 0.38, &mut r);
    let fun2 = random_dna(300_000, 0.38, &mut r);
    let ctrl = random_dna(6_000, 0.5, &mut r);
    // A near-identical strain of the bacterium, indexed as another reference.
    let twin = Genome::new(vec![mutate(&bact, 0.001, &mut r), mutate(&plasmid, 0.001, &mut r)]);
    let far = Genome::new(vec![mutate(&virus_a, 0.07, &mut r)]);
    let genomes = vec![
        Genome::new(vec![virus_a]),
        Genome::new(vec![virus_a2]),
        Genome::new(vec![virus_b]),
        Genome::new(vec![bact, plasmid]),
        Genome::new(vec![kitome]),
        Genome::new(vec![myco]),
        Genome::new(vec![fun1, fun2]),
        Genome::new(vec![ctrl]),
        twin,
    ];
    let mut files = Vec::new();
    for (i, g) in genomes.iter().enumerate() {
        let p = t.path(&format!("ref{i}.fa"));
        let recs: Vec<(String, Vec<u8>)> = g.recs.iter().enumerate().map(|(j, s)| (format!("seq{i}_{j} {}", REFS[i].0), s.clone())).collect();
        if i == 3 {
            // One reference gzipped, to read .gz.
            let gz = t.path(&format!("ref{i}.fa.gz"));
            let mut w = flate2::write::GzEncoder::new(std::fs::File::create(&gz).unwrap(), flate2::Compression::fast());
            for (n, s) in &recs {
                write!(w, ">{n}\n").unwrap();
                w.write_all(s).unwrap();
                w.write_all(b"\n").unwrap();
            }
            w.finish().unwrap();
            files.push(gz);
        } else {
            write_fasta(&p, &recs);
            files.push(p);
        }
    }
    Community {
        host,
        far,
        host_files: vec![host_path],
        names: REFS.iter().map(|(n, k)| (n.to_string(), *k)).collect(),
        genomes,
        files,
    }
}

fn screen_refs(c: &Community) -> Vec<ScreenReference> {
    let mut v: Vec<ScreenReference> = c
        .names
        .iter()
        .zip(&c.files)
        .map(|((n, k), f)| ScreenReference { files: vec![f.clone()], name: n.clone(), kind: *k })
        .collect();
    v.push(ScreenReference { files: c.host_files.clone(), name: "synthetic host".to_string(), kind: AgentKind::Host });
    v
}

fn trial_database(t: &mut T, c: &Community) -> String {
    println!("Database");
    let db_path = t.path("screen.mmgscreen");
    let opt = ScreenDbOptions::default();
    let s = build_screen_database(&screen_refs(c), &db_path, &opt).unwrap();
    t.check(s.references.len() == 10, "10 references in the database");
    t.check(s.host_minimizers > 0 && s.target_minimizers > 0, "the database has agent and host minimizers");
    let (n, h) = sha256_file(&db_path).unwrap();
    t.check(s.sha256 == h && s.file_bytes == n, "the summary gives the SHA-256 and size of the database file");
    for (i, r) in s.references.iter().enumerate().take(9) {
        let (_, hh) = sha256_file(&c.files[i]).unwrap();
        if !t.check(r.sha256 == vec![hh], format!("SHA-256 of reference {} kept", r.name)) {
            break;
        }
    }
    let db = load_screen_database(&db_path).unwrap();
    t.check(db.summary.sha256 == h, "load_screen_database computes the same file SHA-256");
    t.check(db.summary.references == s.references, "the references read back equal those written");
    t.check(
        db.summary.target_minimizers == s.target_minimizers && db.summary.host_minimizers == s.host_minimizers && db.summary.shared_minimizers == s.shared_minimizers,
        "the counts read back equal those written",
    );
    // Windows: a small genome gets ~20 windows of at least 100 bp.
    let va = &s.references[0];
    t.check(va.window_bp == 600 && va.windows == 20, format!("virus A (12 kb): 20 windows of 600 bp (got {} of {})", va.windows, va.window_bp));
    let ba = &s.references[3];
    t.check(ba.window_bp == 1000 && ba.windows == 640, format!("bacterium (640 kb): 640 windows of 1,000 bp (got {} of {})", ba.windows, ba.window_bp));
    t.check(ba.sequences == 2 && ba.bases == 640_000, "the bacterium has 2 sequences, 640,000 bases (read from .gz)");

    // Every agent minimizer is in the index with its reference and the
    // window of its first occurrence; every host minimizer in the host set.
    let mut ok_idx = true;
    let mut ok_win = true;
    let mut checked = 0u64;
    for (i, g) in c.genomes.iter().enumerate() {
        let mut first: HashMap<u64, u64> = HashMap::new();
        let mut off = 0u64;
        for rec in &g.recs {
            for (h, p) in naive_minimizers(rec, 31, 12) {
                if !naive_low(&rec[p..p + 31]) {
                    let e = first.entry(h).or_insert(off + p as u64);
                    *e = (*e).min(off + p as u64);
                }
            }
            off += rec.len() as u64;
        }
        let info = &s.references[i];
        ok_idx &= first.len() as u64 == info.minimizers;
        for (h, pos) in first {
            let occ = db.references_of(h);
            match occ.iter().find(|x| x.0 == i as u32) {
                Some(&(_, gw)) => ok_win &= u64::from(gw - info.first_window) == pos / u64::from(info.window_bp),
                None => ok_idx = false,
            }
            checked += 1;
        }
    }
    t.check(ok_idx, format!("every agent minimizer ({checked}) is in the index under its reference"));
    t.check(ok_win, "each indexed minimizer points to the window of its first occurrence");
    let mut host_ok = true;
    let mut host_n = 0;
    for rec in &c.host.recs {
        for (h, p) in naive_minimizers(rec, 31, 31) {
            if !naive_low(&rec[p..p + 31]) {
                host_ok &= db.host_contains(h);
                host_n += 1;
            }
        }
    }
    t.check(host_ok, format!("every host minimizer ({host_n}) is in the host set"));
    let mut r = Rng(99);
    let fp = (0..1_000_000).filter(|_| db.host_contains(r.next())).count();
    t.check(fp <= 2, format!("random hashes in the host set: {fp} in 1,000,000 (expected ~0)"));
    // The shared minimizers: A and its relative A2 share some.
    let shared_a = db.summary.shared_minimizers;
    t.check(shared_a > 0, format!("virus A and its 4%-divergent relative share minimizers ({shared_a})"));

    // Multi-pass host set: the same set in 1 pass and with a tiny budget.
    let (one, p1) = host_set_with_budget(&c.host_files, 31, 31, 1 << 24).unwrap();
    let (many, pn) = host_set_with_budget(&c.host_files, 31, 31, 20_000).unwrap();
    let (some, pm) = host_set_with_budget(&c.host_files, 31, 31, 100_000).unwrap();
    // ~190,000 host minimizers: 16 passes of 12,000 at a budget of 20,000
    // (60% of it), 4 at 100,000; no pass overflows and doubles them (the
    // minimizers are spread evenly over the passes by their key).
    t.check(p1 == 1 && pn == 16 && pm == 4, format!("the passes follow the budget: 1, 16 and 4 expected ({p1}, {pn}, {pm})"));
    t.check(one == some, "the host set in 4 passes equals the one in 1 pass");
    t.check(one == many, "the host set built in many passes equals the one built in one pass");
    t.check(one.len() as u64 == db.summary.host_minimizers, "and equals the host set of the database");

    // Damaged databases.
    let bytes = std::fs::read(&db_path).unwrap();
    let cut = t.path("cut.mmgscreen");
    std::fs::write(&cut, &bytes[..bytes.len() / 2]).unwrap();
    let e = load_screen_database(&cut).err().unwrap_or_default();
    t.check(e.contains("damaged or incomplete"), format!("a database cut in half is refused: {e}"));
    let mut flip = bytes.clone();
    let m = flip.len() / 3;
    flip[m] ^= 0x10;
    let fl = t.path("flip.mmgscreen");
    std::fs::write(&fl, &flip).unwrap();
    let e = load_screen_database(&fl).err().unwrap_or_default();
    t.check(e.contains("damaged") || e.contains("SHA-256"), format!("a database with one changed byte is refused: {e}"));
    let e = load_screen_database(&c.files[0]).err().unwrap_or_default();
    t.check(e.contains("not an MMGEN screen database"), "a FASTA given as database is refused");
    // A microbio sketch database given by mistake.
    let sk = t.path("sketch.sig");
    microbio::sketch_references(&[microbio::GenomeInput { files: vec![c.files[0].clone()], name: "a".into(), host: false }], &sk, &microbio::SketchOptions::default()).unwrap();
    let e = load_screen_database(&sk).err().unwrap_or_default();
    t.check(e.contains("sketch database"), format!("a sketch database of sketch_references is recognized: {e}"));

    // Options and inputs.
    let xdb = t.path("x.db");
    let bad = |o: ScreenDbOptions| build_screen_database(&screen_refs(c), &xdb, &o).err().unwrap_or_default();
    let e1 = bad(ScreenDbOptions { k: 30, ..Default::default() });
    let e2 = bad(ScreenDbOptions { minimizer_window: 0, ..Default::default() });
    t.check(e1.contains("odd k"), "an even k is refused");
    t.check(e2.contains("minimizer_window"), "minimizer_window 0 is refused");
    let e = build_screen_database(&screen_refs(c), &c.files[0], &ScreenDbOptions::default()).err().unwrap_or_default();
    t.check(e.contains("also one of the input files"), "the output cannot be an input");
    let mut dup = screen_refs(c);
    dup[1].name = dup[0].name.clone();
    let e = build_screen_database(&dup, &t.path("x.db"), &ScreenDbOptions::default()).err().unwrap_or_default();
    t.check(e.contains("same name"), "two references with the same name are refused");
    let only_host = vec![ScreenReference { files: c.host_files.clone(), name: "h".into(), kind: AgentKind::Host }];
    let e = build_screen_database(&only_host, &t.path("x.db"), &ScreenDbOptions::default()).err().unwrap_or_default();
    t.check(e.contains("Every reference is a host"), "a database of hosts only is refused");
    let missing = vec![ScreenReference { files: vec![t.path("no_such.fa")], name: "x".into(), kind: AgentKind::Virus }];
    let e = build_screen_database(&missing, &t.path("x.db"), &ScreenDbOptions::default()).err().unwrap_or_default();
    t.check(e.contains("does not exist"), "a missing reference file is reported");
    t.check(AgentKind::from_text("Mycoplasma").unwrap() == AgentKind::Mycoplasma && AgentKind::from_text("spike-in").unwrap() == AgentKind::Spike, "kinds are read from text");
    t.check(AgentKind::from_text("prion").err().unwrap_or_default().contains("Use one of"), "an unknown kind says the accepted ones");
    // Default names: the first FASTA header.
    let unnamed = vec![
        ScreenReference { files: vec![c.files[0].clone()], name: String::new(), kind: AgentKind::Virus },
        ScreenReference { files: c.host_files.clone(), name: String::new(), kind: AgentKind::Host },
    ];
    let s2 = build_screen_database(&unnamed, &t.path("unnamed.db"), &ScreenDbOptions::default()).unwrap();
    t.check(s2.references[0].name == "seq0_0 virus A" && s2.references[1].name == "host_chr1 synthetic host", "an unnamed reference takes its first header");
    db_path
}

fn trial_read_calls(t: &mut T, c: &Community, db_path: &str) {
    println!("Read calls");
    let db = load_screen_database(db_path).unwrap();
    let mut refs: Vec<(String, Vec<Vec<u8>>, bool)> = c.genomes.iter().zip(&c.names).map(|(g, n)| (n.0.clone(), g.recs.clone(), false)).collect();
    refs.push(("host".into(), c.host.recs.clone(), true));
    let oracle = Oracle::build(&refs, 31, 12, 31);
    let mut r = Rng(5);
    // The one-pass read minimizers (both schemes, rolling low complexity)
    // equal the naive ones.
    let mut same = true;
    for i in 0..3000 {
        let n = 31 + r.below(300);
        let mut q = random_dna(n, 0.2 + 0.6 * r.f(), &mut r);
        if i % 3 == 0 {
            let p = r.below(n);
            for j in p..(p + 20 + r.below(80)).min(n) {
                q[j] = if (j - p) % (1 + i % 4) == 0 { b'C' } else { b'A' };
            }
        }
        if i % 5 == 0 {
            let p = r.below(n);
            q[p] = b'N';
        }
        if i % 7 == 0 {
            q.make_ascii_lowercase();
        }
        let (mut t1, mut h1) = read_minimizer_sets(&db, &q);
        t1.sort();
        t1.dedup();
        h1.sort();
        h1.dedup();
        same &= t1 == naive_read_set(&[&q], 31, 12) && h1 == naive_read_set(&[&q], 31, 31);
    }
    t.check(same, "the one-pass read minimizers equal the naive ones (3,000 reads with N, lower case, repeats)");
    let opts = [
        ScreenOptions::default(),
        ScreenOptions { unique_margin: 0.0, min_fraction: 0.0, min_hits: 1, ..Default::default() },
        ScreenOptions { unique_margin: 0.5, min_fraction: 0.5, min_hits: 6, host_min_hits: 3, ..Default::default() },
    ];
    let mut agree = 0;
    let mut total = 0;
    let mut first_bad: Option<String> = None;
    for o in &opts {
        for i in 0..1200 {
            let src = i % 10;
            let (a, _, b, _) = match src {
                0..=7 => sim_pair(&c.genomes[src], &ILLUMINA, None, &mut r),
                8 => sim_pair(&c.host, &ILLUMINA, None, &mut r),
                _ => {
                    // A chimera: half host, half bacterium; or a short read.
                    let (h1, q1, _, _) = sim_pair(&c.host, &ILLUMINA, None, &mut r);
                    let (b1, _, b2, q2) = sim_pair(&c.genomes[3], &ILLUMINA, None, &mut r);
                    ([&h1[..75], &b1[75..]].concat(), q1, b2[..40 + r.below(60)].to_vec(), q2)
                }
            };
            let single = i % 7 == 0;
            let got = classify_one(&db, &a, if single { None } else { Some(&b) }, o).outcome;
            let mates: Vec<&[u8]> = if single { vec![&a] } else { vec![&a, &b] };
            let want = oracle.call(&mates, o);
            total += 1;
            if got == want {
                agree += 1;
            } else if first_bad.is_none() {
                first_bad = Some(format!("{got:?} vs oracle {want:?} (source {src})"));
            }
        }
    }
    t.check(agree == total, format!("read calls equal the oracle's in {agree}/{total} reads and pairs (3 option sets){}", first_bad.map_or(String::new(), |x| format!("; first difference: {x}"))));
    // Perfect reads and the expected outcomes.
    let o = ScreenOptions::default();
    let mut counts = vec![0u32; 10];
    for (src, g) in c.genomes.iter().enumerate() {
        let mut n = 0;
        for _ in 0..300 {
            let (a, _, b, _) = sim_pair(g, &ILLUMINA, None, &mut r);
            let got = classify_one(&db, &a, Some(&b), &o).outcome;
            // The bacterium and its twin: a tie between the two is right.
            let twins = src == 3 || src == 8;
            if got == Outcome::Unique(src as u32) || (twins && got == Outcome::Tie(vec![3, 8])) {
                n += 1;
            }
        }
        counts[src] = n;
    }
    t.check(counts.iter().take(9).all(|&n| n >= 290), format!("simulated pairs of each reference go to it, or tie between the twins (of 300: {:?})", &counts[..9]));
    let mut host = 0;
    for _ in 0..1000 {
        let (a, _, b, _) = sim_pair(&c.host, &ILLUMINA, None, &mut r);
        host += u32::from(classify_one(&db, &a, Some(&b), &o).outcome == Outcome::Host);
    }
    t.check(host >= 995, format!("host pairs are host: {host}/1000"));
    let rnd = random_dna(150, 0.5, &mut r);
    let rnd2 = random_dna(150, 0.5, &mut r);
    t.check(classify_one(&db, &rnd, Some(&rnd2), &o).outcome == Outcome::Unclassified, "a random pair is unclassified");
    t.check(classify_one(&db, &vec![b'N'; 150], None, &o).outcome == Outcome::NoKmers, "a read of N has no k-mers");
    let micro: Vec<u8> = (0..150).map(|i| if i % 2 == 0 { b'T' } else { b'G' }).collect();
    t.check(classify_one(&db, &micro, None, &o).outcome == Outcome::NoKmers, "a microsatellite read has no usable k-mers");
}

/// Writes a paired sample made of parts (tag, genome, pairs, region).
fn write_sample(t: &T, name: &str, parts: &[(&str, &Genome, u64, Option<(usize, usize, usize)>)], seed: u64, gz: bool) -> (String, String) {
    let ext = if gz { "fq.gz" } else { "fq" };
    let p1 = t.path(&format!("{name}_R1.{ext}"));
    let p2 = t.path(&format!("{name}_R2.{ext}"));
    let mut w = PairWriter::new(&p1, &p2, &format!("{name}:"));
    let mut r = Rng(seed);
    // Interleave the parts at random, as in a real run.
    let mut left: Vec<u64> = parts.iter().map(|p| p.2).collect();
    let mut total: u64 = left.iter().sum();
    while total > 0 {
        let mut x = r.below(total as usize) as u64;
        let mut i = 0;
        while x >= left[i] {
            x -= left[i];
            i += 1;
        }
        let (tag, g, _, region) = parts[i];
        let pair = sim_pair(g, &ILLUMINA, region, &mut r);
        w.put(tag, pair);
        left[i] -= 1;
        total -= 1;
    }
    w.finish();
    (p1, p2)
}

fn sample(name: &str, f: (String, String)) -> ScreenSample {
    ScreenSample { name: name.to_string(), r1: vec![f.0], r2: vec![f.1] }
}

fn row<'a>(s: &'a SampleResult, name: &str) -> Option<&'a OrganismResult> {
    s.organisms.iter().find(|r| r.name == name)
}

fn trial_screen(t: &mut T, c: &Community, db_path: &str) {
    println!("Screen");
    let g = &c.genomes;
    let host = &c.host;
    // RAW: host + virus A + bacterium + mycoplasma + a vector element of
    // virus B (600 bp, piled) + the process control + the kitome.
    let raw = write_sample(
        t,
        "RAW",
        &[
            ("host", host, 40_000, None),
            ("virusA", &g[0], 40, None),
            ("bact", &g[3], 30, None),
            ("myco", &g[5], 25, None),
            ("vector", &g[2], 300, Some((0, 2000, 2600))),
            ("ctrl", &g[7], 100, None),
            ("kitome", &g[4], 20, None),
        ],
        1,
        true,
    );
    let blank = write_sample(t, "BLANK", &[("host", host, 2000, None), ("kitome", &g[4], 30, None), ("ctrl", &g[7], 100, None)], 2, false);
    let rel = write_sample(t, "REL", &[("host", host, 20_000, None), ("virusA2", &g[1], 60, None)], 3, false);
    let hostonly = write_sample(t, "HOSTONLY", &[("host", host, 20_000, None)], 4, false);
    let relfar = write_sample(t, "RELFAR", &[("host", host, 5_000, None), ("virusAfar", &c.far, 400, None)], 6, false);
    // RAW again, in two lanes (the same reads).
    let lanes = split_lanes(t, &raw, "RAWLANES");
    let spikes = t.path("spikes.tsv");
    std::fs::write(
        &spikes,
        "# the spike-in table\nsample\torganism\tamount\tunit\n*\tprocess control\t100\treads\nRAW\tvirus A\t40\treads\nRAW\tbacterium\t30\tRPM\nREL\tvirus A2 (relative of A)\t60\tcopies\n",
    )
    .unwrap();
    let samples = vec![
        sample("RAW", raw.clone()),
        sample("BLANK", blank),
        sample("REL", rel),
        sample("HOSTONLY", hostonly),
        ScreenSample { name: "RAWLANES".into(), r1: lanes.0, r2: lanes.1 },
        sample("RELFAR", relfar),
    ];
    let o = ScreenOptions {
        negative_controls: vec!["BLANK".into()],
        spike_table: Some(spikes.clone()),
        output_tsv: Some(t.path("report.tsv")),
        output_json: Some(t.path("report.json")),
        output_html: Some(t.path("report.html")),
        ..Default::default()
    };
    let rep = match contaminant_screen(&samples, db_path, &o) {
        Ok(r) => r,
        Err(e) => {
            t.check(false, format!("contaminant_screen ran: {e}"));
            return;
        }
    };
    let s = &rep.samples[0];
    println!("{}", rep.resumo().lines().take(40).collect::<Vec<_>>().join("\n"));
    t.check(s.fragments == 40_515, format!("RAW has 40,515 pairs (got {})", s.fragments));
    t.check(s.host >= 39_900 && s.host <= 40_000, format!("RAW: host pairs set aside {} (40,000 simulated)", s.host));
    let call = |s: &SampleResult, n: &str| row(s, n).map_or(Call::NotDetected, |r| r.call);
    t.check(call(s, "virus A") == Call::Detected, "RAW: virus A detected");
    t.check(call(s, "bacterium") == Call::Detected, "RAW: the bacterium detected");
    t.check(call(s, "mycoplasma") == Call::Detected, "RAW: the mycoplasma detected");
    t.check(call(s, "process control") == Call::Detected, "RAW: the process control detected");
    let vb = row(s, "virus B");
    t.check(
        vb.is_some_and(|r| r.call == Call::BelowThreshold && r.unique_reads >= 250 && r.reasons.iter().any(|x| x.contains("concentrated"))),
        format!("RAW: 300 pairs piled on 600 bp of virus B are below threshold, as concentrated ({:?})", vb.map(|r| (r.call, r.unique_reads, r.windows, r.expected_windows))),
    );
    let kit = row(s, "kitome bacterium");
    t.check(
        kit.is_some_and(|r| r.call == Call::Background && r.blank.as_deref() == Some("BLANK") && r.blank_ratio.is_some_and(|x| x < 10.0)),
        format!("RAW: the kitome bacterium, also in the blank, is background ({:?})", kit.map(|r| (r.call, r.rpm, r.blank_rpm))),
    );
    t.check(call(s, "fungus") == Call::NotDetected && call(s, "virus A2 (relative of A)") != Call::Detected, "RAW: the fungus and the relative A2 are not detected");
    let tw = row(s, "bacterium twin (0.1% divergent)");
    t.check(
        tw.is_some_and(|x| x.call == Call::BelowThreshold && x.reasons.iter().any(|r| r.contains("tie with bacterium") || r.contains("mostly with bacterium"))),
        format!("RAW: the near-identical twin is not called; its tied reads go to the bacterium ({:?})", tw.map(|x| (x.call, x.unique_reads, x.tied_reads, x.won_reads, x.reasons.clone()))),
    );
    let bw = row(s, "bacterium").unwrap();
    t.check(bw.won_reads > 0 && bw.unique_fraction.is_some_and(|f| f > 0.9), format!("RAW: the bacterium wins the ties with its twin ({} unique, {} won)", bw.unique_reads, bw.won_reads));
    t.check(s.detected.len() == 3, format!("RAW: 3 agents detected ({:?})", s.detected));
    let ctrl = row(s, "process control").unwrap();
    t.check(ctrl.recovery.is_some_and(|x| (0.9..=1.02).contains(&x)), format!("RAW: the process control recovery {:?} (100 spiked)", ctrl.recovery));
    t.check(ctrl.kind == AgentKind::Spike && !s.detected.contains(&ctrl.name), "the process control is not listed as a contaminant");
    let va = row(s, "virus A").unwrap();
    t.check(va.recovery.is_some_and(|x| (0.85..=1.02).contains(&x)), format!("RAW: virus A recovery {:?} (40 reads spiked)", va.recovery));
    t.check(va.identity.is_some_and(|x| x > 0.99), format!("RAW: virus A identity {:?}", va.identity));
    let ba = row(s, "bacterium").unwrap();
    let exp = 30.0 * s.fragments as f64 / 1e6;
    t.check(
        ba.recovery.is_some_and(|x| (ba.reads / exp - x).abs() < 1e-9),
        format!("RAW: the bacterium spiked as 30 RPM: recovery = reads / (30 x fragments / 1e6) ({:?})", ba.recovery),
    );
    t.check(ba.windows >= 25 && ba.dispersion.is_some_and(|d| d > 0.9), format!("RAW: 30 bacterial pairs hit {} windows (dispersion {:?})", ba.windows, ba.dispersion));
    // BLANK.
    let b = &rep.samples[1];
    t.check(b.negative_control && call(b, "kitome bacterium") == Call::Detected, "BLANK: the kitome is detected in the blank itself");
    // REL: A2 is there, A is not.
    let rl = &rep.samples[2];
    t.check(call(rl, "virus A2 (relative of A)") == Call::Detected, "REL: the relative A2 is detected");
    t.check(call(rl, "virus A") != Call::Detected, format!("REL: virus A is not detected from A2's reads ({:?})", row(rl, "virus A").map(|r| (r.call, r.unique_reads, r.tied_reads))));
    let a2 = row(rl, "virus A2 (relative of A)").unwrap();
    t.check(a2.spiked_unit.as_deref() == Some("copies") && a2.recovery.is_none(), "a spike in copies has no recovery fraction");
    // HOSTONLY: nothing.
    let h = &rep.samples[3];
    t.check(h.detected.is_empty() && h.organisms.iter().all(|r| r.call != Call::Detected), format!("HOSTONLY: nothing detected ({:?})", h.detected));
    t.check(h.host as f64 >= 0.997 * h.fragments as f64, format!("HOSTONLY: {} of {} host", h.host, h.fragments));
    // The lanes: the same numbers.
    let ln = &rep.samples[4];
    let same = ln.fragments == s.fragments
        && ln.host == s.host
        && ln.organisms.len() == s.organisms.len()
        && ln.organisms.iter().zip(&s.organisms).all(|(a, b)| a.name == b.name && a.unique_reads == b.unique_reads && a.windows == b.windows && a.call == b.call);
    t.check(same, "two lanes give the same result as one file");
    t.check(ln.files.len() == 4, "the lanes sample lists its 4 files");
    // RELFAR: reads of a 7%-divergent relative not in the database.
    let rf = &rep.samples[5];
    let fa = row(rf, "virus A");
    t.check(
        fa.is_some_and(|x| x.call != Call::Detected && x.weak_reads > x.unique_reads),
        format!("RELFAR: a relative not in the database (7% divergent) is not called virus A ({:?})", fa.map(|x| (x.call, x.unique_reads, x.weak_reads))),
    );
    t.check(
        fa.is_some_and(|x| x.unique_reads == 0 || x.reasons.iter().any(|r| r.contains("only weakly"))),
        "RELFAR: the reason says the reads match only weakly (a relative not in the database)",
    );
    // Spike rows: "*" goes to every sample that is not a blank.
    let n_ctrl = rep.spikes.iter().filter(|x| x.organism == "process control").count();
    t.check(n_ctrl == 6, format!("the '*' spike-in of a process control applies to every sample, blanks too ({n_ctrl})"));
    let blank_ctrl = rep.samples[1].organisms.iter().find(|x| x.name == "process control");
    t.check(blank_ctrl.is_some_and(|x| x.call == Call::Detected && x.recovery.is_some_and(|r| r > 0.95)), "the process control is recovered in the blank");
    t.check(rep.lod.iter().all(|l| l.organism != "process control"), "a process control has a recovery, not a limit of detection");
    // Traceability.
    for f in &s.files {
        let (n, h) = sha256_file(&f.path).unwrap();
        if !t.check(f.sha256 == h && f.bytes == n, format!("the report has the SHA-256 of {}", f.path)) {
            break;
        }
    }
    let (_, dbh) = sha256_file(db_path).unwrap();
    t.check(rep.database.sha256 == dbh, "the report has the SHA-256 of the database");
    let json = std::fs::read_to_string(t.path("report.json")).unwrap();
    match serde_json::from_str::<serde_json::Value>(&json) {
        Ok(v) => {
            t.check(v["notice"].as_str().is_some_and(|x| x.contains("GMP") && x.contains("ICH Q5A(R2)")), "the JSON carries the notice (GMP, ICH Q5A(R2))");
            t.check(v["database"]["sha256"].as_str() == Some(dbh.as_str()), "the JSON has the database SHA-256");
            t.check(v["samples"][0]["files"][0]["sha256"].as_str().is_some_and(|x| x.len() == 64), "the JSON has the input SHA-256");
            t.check(v["parameters"]["min_blank_ratio"].as_f64() == Some(10.0) && v["engine"].as_str().is_some_and(|x| x.contains("contaminant")), "the JSON has the parameters and the engine");
            t.check(v["date"].as_str().is_some_and(|x| x.len() == 20 && x.ends_with('Z')), "the JSON has the date (UTC)");
            t.check(v["samples"][0]["organisms"].as_array().is_some_and(|a| a.iter().any(|o| o["call"] == "background")), "the JSON lists the background call");
            t.check(v["lod"].as_array().is_some_and(|a| !a.is_empty()), "the JSON has the limit-of-detection table");
        }
        Err(e) => {
            t.check(false, format!("the JSON report is valid JSON: {e}"));
        }
    }
    let tsv = std::fs::read_to_string(t.path("report.tsv")).unwrap();
    let rows: Vec<&str> = tsv.lines().filter(|l| !l.starts_with('#')).collect();
    let ncol = rows[0].split('\t').count();
    t.check(rows.iter().all(|l| l.split('\t').count() == ncol), format!("every TSV line has the {ncol} columns of the header"));
    t.check(tsv.lines().any(|l| l.starts_with("# database") && l.contains(&dbh)), "the TSV header has the database SHA-256");
    let html = std::fs::read_to_string(t.path("report.html")).unwrap();
    t.check(
        html.contains("For research and teaching") && html.contains("ICH Q5A(R2)") && html.contains(&dbh) && !html.contains("<script"),
        "the HTML report has the notice and the SHA-256, and no script",
    );
    // Determinism: the same report again, but the date and the times.
    let o2 = ScreenOptions { output_tsv: None, output_json: None, output_html: None, ..o.clone() };
    let rep2 = contaminant_screen(&samples, db_path, &o2).unwrap();
    let strip = |r: &ScreenReport| {
        let mut v: serde_json::Value = serde_json::from_str(&r.json()).unwrap();
        v["date"] = serde_json::Value::Null;
        v["seconds"] = serde_json::Value::Null;
        v["text"] = serde_json::Value::Null;
        for s in v["samples"].as_array_mut().unwrap() {
            s["seconds"] = serde_json::Value::Null;
        }
        v
    };
    t.check(strip(&rep) == strip(&rep2), "the same inputs give the same report (but the date and times)");

    // Errors.
    let e = contaminant_screen(&samples, db_path, &ScreenOptions { negative_controls: vec!["NOPE".into()], ..Default::default() }).err().unwrap_or_default();
    t.check(e.contains("not one of the samples"), "an unknown negative control is refused");
    let bad_spike = ScreenOptions { spike_ins: vec![SpikeIn { sample: None, organism: "virus Z".into(), amount: 5.0, unit: "reads".into() }], ..Default::default() };
    let e = contaminant_screen(&samples, db_path, &bad_spike).err().unwrap_or_default();
    t.check(e.contains("not a reference of the screen database"), "an unknown spike-in organism is refused");
    let e = contaminant_screen(&samples, &c.files[0], &ScreenOptions::default()).err().unwrap_or_default();
    t.check(e.contains("not an MMGEN screen database"), "a FASTA as database is refused");
    let out_in = ScreenOptions { output_tsv: Some(raw.0.clone()), ..Default::default() };
    let e = contaminant_screen(&samples, db_path, &out_in).err().unwrap_or_default();
    t.check(e.contains("also one of the input files"), "an output that is an input is refused");
    // Unpaired names and counts.
    let (b1, b2) = (t.path("bad_R1.fq"), t.path("bad_R2.fq"));
    std::fs::write(&b1, "@r1 1\nACGTACGTAC\n+\nIIIIIIIIII\n@r2 1\nACGTACGTAC\n+\nIIIIIIIIII\n").unwrap();
    std::fs::write(&b2, "@r1 2\nACGTACGTAC\n+\nIIIIIIIIII\n@rX 2\nACGTACGTAC\n+\nIIIIIIIIII\n").unwrap();
    let e = contaminant_screen(&[ScreenSample { name: "x".into(), r1: vec![b1.clone()], r2: vec![b2.clone()] }], db_path, &ScreenOptions::default()).err().unwrap_or_default();
    t.check(e.contains("not paired"), format!("reads with different names are not paired: {e}"));
    std::fs::write(&b2, "@r1 2\nACGTACGTAC\n+\nIIIIIIIIII\n").unwrap();
    let e = contaminant_screen(&[ScreenSample { name: "x".into(), r1: vec![b1.clone()], r2: vec![b2.clone()] }], db_path, &ScreenOptions::default()).err().unwrap_or_default();
    t.check(e.contains("same number of reads"), format!("files of a pair with different counts are refused: {e}"));
    let e = contaminant_screen(&[ScreenSample { name: "x".into(), r1: vec![b1.clone()], r2: vec![b1.clone()] }], db_path, &ScreenOptions::default()).err().unwrap_or_default();
    t.check(e.contains("same file"), "the same file as both mates is refused");
    std::fs::write(&spikes, "organism\tamount\nvirus A\tlots\n").unwrap();
    let e = read_spike_table(&spikes).err().unwrap_or_default();
    t.check(e.contains("line 2") && e.contains("not a number"), format!("a spike amount that is not a number is refused with its line: {e}"));
    std::fs::write(&spikes, "name\tvalue\nvirus A\t3\n").unwrap();
    let e = read_spike_table(&spikes).err().unwrap_or_default();
    t.check(e.contains("organism and amount"), "a spike table without the columns is refused");
}

fn split_lanes(t: &T, f: &(String, String), name: &str) -> (Vec<String>, Vec<String>) {
    let read = |p: &str| -> Vec<String> {
        let mut rd = microbio::SeqReader::open(p).unwrap();
        let mut rec = microbio::SeqRecord::default();
        let mut v = Vec::new();
        while rd.next(&mut rec).unwrap() {
            v.push(format!(
                "@{}\n{}\n+\n{}\n",
                String::from_utf8_lossy(&rec.head),
                String::from_utf8_lossy(&rec.seq),
                String::from_utf8_lossy(&rec.qual)
            ));
        }
        v
    };
    let (a, b) = (read(&f.0), read(&f.1));
    let cut = a.len() / 3;
    let mut r1 = Vec::new();
    let mut r2 = Vec::new();
    for (lane, range) in [(1, 0..cut), (2, cut..a.len())] {
        let p1 = t.path(&format!("{name}_L00{lane}_R1.fq"));
        let p2 = t.path(&format!("{name}_L00{lane}_R2.fq"));
        std::fs::write(&p1, a[range.clone()].concat()).unwrap();
        std::fs::write(&p2, b[range].concat()).unwrap();
        r1.push(p1);
        r2.push(p2);
    }
    (r1, r2)
}

/// The limit of detection on a small series: 4 levels x 3 replicates.
fn trial_lod(t: &mut T, c: &Community, db_path: &str) {
    println!("Limit of detection");
    let mut samples = Vec::new();
    let mut spikes = Vec::new();
    let levels = [1u64, 3, 10, 30];
    let mut seed = 100;
    for &l in &levels {
        for rep in 0..3 {
            seed += 1;
            let name = format!("L{l}_r{rep}");
            let f = write_sample(t, &name, &[("host", &c.host, 3000, None), ("virusA", &c.genomes[0], l, None)], seed, false);
            samples.push(sample(&name, f));
            spikes.push(SpikeIn { sample: Some(name), organism: "virus A".into(), amount: l as f64, unit: "reads".into() });
        }
    }
    let rep = contaminant_screen(&samples, db_path, &ScreenOptions { spike_ins: spikes, ..Default::default() }).unwrap();
    t.check(rep.lod.len() == 1 && rep.lod[0].levels.len() == 4, "one LoD table, 4 levels");
    let l = &rep.lod[0];
    let rates: Vec<f64> = l.levels.iter().map(|v| v.rate).collect();
    println!("  virus A: detection rate by level {:?}", l.levels.iter().map(|v| (v.amount, v.detected, v.replicates)).collect::<Vec<_>>());
    t.check(l.levels[0].detected == 0, "1 read is never a detection (min_reads 3)");
    t.check(l.levels[3].detected == 3 && l.levels[2].detected == 3, "10 and 30 reads are detected in every replicate");
    // The LoD follows from the rates.
    let mut want = None;
    for v in l.levels.iter().rev() {
        if v.rate >= 0.95 {
            want = Some(v.amount);
        } else {
            break;
        }
    }
    t.check(l.lod == want && want.is_some(), format!("LoD = lowest level with all higher levels at 95% ({:?}; rates {:?})", l.lod, rates));
    t.check(l.lowest_detected.is_some_and(|x| x <= l.lod.unwrap_or(f64::MAX)), "the lowest level detected is at most the LoD");
    t.check(l.levels.iter().all(|v| v.mean_recovery.is_some_and(|r| r > 0.8 && r <= 1.05) || v.mean_reads == 0.0), "the recovery of the spiked reads is 80-105% at every level");
}

// ---------------------------------------------------------------------
//  Real data
// ---------------------------------------------------------------------

#[cfg(windows)]
fn peak_memory_mb() -> f64 {
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
    extern "system" {
        fn GetCurrentProcess() -> isize;
        fn K32GetProcessMemoryInfo(h: isize, p: *mut Pmc, cb: u32) -> i32;
    }
    let mut p = Pmc {
        cb: std::mem::size_of::<Pmc>() as u32,
        page_fault_count: 0,
        peak_working_set_size: 0,
        working_set_size: 0,
        quota_peak_paged_pool_usage: 0,
        quota_paged_pool_usage: 0,
        quota_peak_non_paged_pool_usage: 0,
        quota_non_paged_pool_usage: 0,
        pagefile_usage: 0,
        peak_pagefile_usage: 0,
    };
    // SAFETY: a plain Win32 call filling a struct of the documented layout.
    let ok = unsafe { K32GetProcessMemoryInfo(GetCurrentProcess(), &mut p, p.cb) };
    if ok == 0 {
        return 0.0;
    }
    p.peak_working_set_size as f64 / 1e6
}

#[cfg(not(windows))]
fn peak_memory_mb() -> f64 {
    let s = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    s.lines().find(|l| l.starts_with("VmHWM")).and_then(|l| l.split_whitespace().nth(1)).and_then(|v| v.parse::<f64>().ok()).map_or(0.0, |kb| kb / 1e3)
}

fn rd(base: &str, rel: &str) -> String {
    Path::new(base).join(rel).to_string_lossy().to_string()
}

/// The research_data references (and the mycoplasmas and lambda), with
/// short unique names.
fn real_refs(data: &str, myco: &str, lambda: &str, skip: &[&str]) -> Vec<ScreenReference> {
    let tsv = std::fs::read_to_string(rd(data, "reference_genomes/genomes.tsv")).unwrap();
    let rows: Vec<Vec<String>> = tsv.lines().skip(1).map(|l| l.split('\t').map(str::to_string).collect()).collect();
    let mut out = Vec::new();
    for r in &rows {
        let (cat, acc, asm, org, strain) = (&r[0], &r[1], &r[2], &r[3], &r[5]);
        let (dir, kind) = if cat.starts_with("bacterium") {
            ("bacteria", AgentKind::Bacterium)
        } else if cat.starts_with("fungus") {
            ("fungi", AgentKind::Fungus)
        } else {
            ("viruses", AgentKind::Virus)
        };
        let dup = rows.iter().filter(|x| &x[3] == org).count() > 1;
        let name = if dup { format!("{org} {strain}") } else { org.clone() };
        if skip.iter().any(|s| name.contains(s)) {
            continue;
        }
        out.push(ScreenReference { files: vec![rd(data, &format!("reference_genomes/{dir}/{acc}_{asm}_genomic.fna.gz"))], name, kind });
    }
    let mut mf: Vec<String> = std::fs::read_dir(myco).unwrap().map(|e| e.unwrap().path().to_string_lossy().to_string()).filter(|p| p.ends_with("_genomic.fna.gz")).collect();
    mf.sort();
    for f in mf {
        let head = read_fasta(&f)[0].0.clone();
        let words: Vec<&str> = head.split_whitespace().collect();
        let name = format!("{} {}", words[1], words[2]);
        out.push(ScreenReference { files: vec![f], name, kind: AgentKind::Mycoplasma });
    }
    out.push(ScreenReference { files: vec![lambda.to_string()], name: "Escherichia phage lambda (process control)".into(), kind: AgentKind::Spike });
    out.push(ScreenReference {
        files: vec![rd(data, "metagenome/host_decoy/GRCh38.p14_chr21.fna.gz"), rd(data, "metagenome/host_decoy/GRCh38.p14_chrMT.fna.gz")],
        name: "Homo sapiens GRCh38 chr21 + chrM".into(),
        kind: AgentKind::Host,
    });
    out
}

fn real_host_index(data: &str, work: &str) {
    std::fs::create_dir_all(work).unwrap();
    let w = |n: &str| rd(work, n);
    let host = vec![rd(data, "metagenome/host_decoy/GRCh38.p14_chr21.fna.gz"), rd(data, "metagenome/host_decoy/GRCh38.p14_chrMT.fna.gz")];
    let refs = vec![
        ScreenReference { files: host.clone(), name: "human chr21 + chrM".into(), kind: AgentKind::Host },
        ScreenReference { files: vec![rd(data, "reference_genomes/viruses/GCF_000845085.1_ViralProj14518_genomic.fna.gz")], name: "Human mastadenovirus C".into(), kind: AgentKind::Virus },
        ScreenReference { files: vec![rd(data, "reference_genomes/bacteria/GCF_000005845.2_ASM584v2_genomic.fna.gz")], name: "Escherichia coli K-12".into(), kind: AgentKind::Bacterium },
    ];
    let m0 = peak_memory_mb();
    let t0 = Instant::now();
    let s = build_screen_database(&refs, &w("host_chr21.mmgscreen"), &ScreenDbOptions::default()).unwrap();
    let secs = t0.elapsed().as_secs_f64();
    let m1 = peak_memory_mb();
    print!("{}", s.resumo());
    let hb = s.references[0].bases;
    let hk = s.references[0].kmers;
    println!("BUILD chr21+chrM: {secs:.1} s; peak working set {m1:.0} MB (before {m0:.0} MB); host minimizers {} for {} bases ({} valid k-mers): {:.2} per 100 valid k-mers; host set {:.1} MB ({:.3} bytes per valid base); bucket bits {}",
        s.host_minimizers, hb, hk, 100.0 * s.host_minimizers as f64 / hk as f64, s.host_index_bytes as f64 / 1e6, s.host_index_bytes as f64 / hk as f64, s.host_bucket_bits);
    // Extrapolation: GRCh38 has ~2.95 billion valid (non-N) bases, ~2.5
    // billion distinct 31-mers.
    let per_kmer = s.host_minimizers as f64 / hk as f64;
    let human_valid = 2.95e9;
    let est = per_kmer * human_valid;
    println!(
        "EXTRAPOLATION full human (2.95 Gb non-N): ~{:.0} M host minimizers (upper bound: chr21 repeats shared with other chromosomes collapse), host set ~{:.0} MB + directory 64 MB; build time ~{:.0} s at this speed (one pass with memory_mb >= {:.0})",
        est / 1e6,
        est * 4.0 / 1e6,
        secs * human_valid / hk as f64,
        est * 8.0 / 0.6 / 1e6
    );
    // Multi-pass equivalence on chr21.
    let t1 = Instant::now();
    let (one, p1) = host_set_with_budget(&host, 31, 31, 1 << 28).unwrap();
    let ta = t1.elapsed().as_secs_f64();
    let t2 = Instant::now();
    let (many, pn) = host_set_with_budget(&host, 31, 31, 400_000).unwrap();
    let tb = t2.elapsed().as_secs_f64();
    println!("MULTIPASS: {p1} pass {ta:.1} s; {pn} passes (budget 400,000 entries = 3.2 MB) {tb:.1} s; identical sets: {}", one == many);
    // Host depletion on simulated reads.
    let db = load_screen_database(&w("host_chr21.mmgscreen")).unwrap();
    let g = Genome::from_fasta(&host);
    let ecoli = Genome::from_fasta(&[rd(data, "reference_genomes/bacteria/GCF_000005845.2_ASM584v2_genomic.fna.gz")]);
    let ad = Genome::from_fasta(&[rd(data, "reference_genomes/viruses/GCF_000845085.1_ViralProj14518_genomic.fna.gz")]);
    let o = ScreenOptions::default();
    let mut r = Rng(21);
    for (len, label) in [(150usize, "2x150"), (100, "2x100"), (75, "2x75")] {
        let sim = Sim { len, insert: (len * 2 + 50) as f64, ..ILLUMINA };
        let n = 200_000;
        let t3 = Instant::now();
        let pairs: Vec<_> = (0..n).map(|_| sim_pair(&g, &sim, None, &mut r)).collect();
        let host_calls: usize = {
            use rayon::prelude::*;
            pairs.par_iter().filter(|p| classify_one(&db, &p.0, Some(&p.2), &o).outcome == Outcome::Host).count()
        };
        let el = t3.elapsed().as_secs_f64();
        let mut other = 0usize;
        let mut ad_host = 0usize;
        for _ in 0..20_000 {
            let p = sim_pair(&ecoli, &sim, None, &mut r);
            other += usize::from(classify_one(&db, &p.0, Some(&p.2), &o).outcome == Outcome::Host);
            let q = sim_pair(&ad, &sim, None, &mut r);
            ad_host += usize::from(classify_one(&db, &q.0, Some(&q.2), &o).outcome == Outcome::Host);
        }
        println!(
            "HOST DEPLETION {label}: {host_calls}/{n} simulated chr21 pairs set aside ({:.3}%); E. coli pairs taken as host {other}/20000, adenovirus pairs taken as host {ad_host}/20000; {:.1} s",
            100.0 * host_calls as f64 / n as f64,
            el
        );
    }
    println!("PEAK working set of the whole run: {:.0} MB", peak_memory_mb());
}

struct SpikeOrg {
    name: &'static str,
    file: String,
    region: Option<(usize, usize, usize)>,
}

fn real_spike(data: &str, myco: &str, lambda: &str, work: &str) {
    std::fs::create_dir_all(work).unwrap();
    let w = |n: &str| rd(work, n);
    // The panel: every reference but Zika MR 766 (spiked as a relative not
    // in the database).
    let refs = real_refs(data, myco, lambda, &["Zika virus MR 766"]);
    let dbp = w("panel.mmgscreen");
    if !Path::new(&dbp).exists() {
        let t0 = Instant::now();
        let s = build_screen_database(&refs, &dbp, &ScreenDbOptions::default()).unwrap();
        println!("panel database: {} references, {:.0} s, peak {:.0} MB", s.references.len(), t0.elapsed().as_secs_f64(), peak_memory_mb());
        print!("{}", s.resumo());
    }
    let host = Genome::from_fasta(&[rd(data, "metagenome/host_decoy/GRCh38.p14_chr21.fna.gz"), rd(data, "metagenome/host_decoy/GRCh38.p14_chrMT.fna.gz")]);
    let v = |f: &str| rd(data, &format!("reference_genomes/viruses/{f}"));
    let b = |f: &str| rd(data, &format!("reference_genomes/bacteria/{f}"));
    let spiked = [
        SpikeOrg { name: "Human mastadenovirus C", file: v("GCF_000845085.1_ViralProj14518_genomic.fna.gz"), region: None },
        SpikeOrg { name: "Hepatitis B virus", file: v("GCF_000861825.2_ViralProj15428_genomic.fna.gz"), region: None },
        SpikeOrg { name: "Norovirus GII", file: v("GCF_003638685.1_ASM363868v1_genomic.fna.gz"), region: None },
        SpikeOrg { name: "Mesomycoplasma hyorhinis", file: rd(myco, "GCF_900476065.1_50465_F02_genomic.fna.gz"), region: None },
        SpikeOrg { name: "Escherichia coli str. K-12 substr. MG1655", file: b("GCF_000005845.2_ASM584v2_genomic.fna.gz"), region: None },
        SpikeOrg { name: "Saccharomyces cerevisiae S288C", file: rd(data, "reference_genomes/fungi/GCF_000146045.2_R64_genomic.fna.gz"), region: None },
        SpikeOrg { name: "Zika virus MR 766 (not in the database)", file: v("GCF_000882815.3_ViralProj36615_genomic.fna.gz"), region: None },
    ];
    let genomes: Vec<Genome> = spiked.iter().map(|s| Genome::from_fasta(&[s.file.clone()])).collect();
    let lam = Genome::from_fasta(&[lambda.to_string()]);
    let pae = Genome::from_fasta(&[b("GCF_000006765.1_ASM676v1_genomic.fna.gz")]);
    let adv = Genome::from_fasta(&[spiked[0].file.clone()]);
    let depth: u64 = 1_000_000;
    let mut samples: Vec<ScreenSample> = Vec::new();
    let mut spikes: Vec<SpikeIn> = Vec::new();
    let mut truth: HashMap<String, Vec<String>> = HashMap::new();
    let bg = (w("bg0_R1.fq.gz"), w("bg0_R2.fq.gz"));
    let t0 = Instant::now();
    if !Path::new(&bg.1).exists() {
        let mut pw = PairWriter::new(&bg.0, &bg.1, "bg0:");
        pw.put_n("host", &host, depth, None, &mut Rng(1000));
        pw.finish();
    }
    for i in 1..3u64 {
        let p = (w(&format!("bg{i}_R1.fq.gz")), w(&format!("bg{i}_R2.fq.gz")));
        if !Path::new(&p.1).exists() {
            let mut pw = PairWriter::new(&p.0, &p.1, &format!("bg{i}:"));
            pw.put_n("host", &host, depth, None, &mut Rng(1000 + i));
            pw.finish();
        }
        samples.push(ScreenSample { name: format!("HOST_ONLY_{i}"), r1: vec![p.0], r2: vec![p.1] });
        truth.insert(format!("HOST_ONLY_{i}"), vec![]);
    }
    samples.push(ScreenSample { name: "HOST_ONLY_0".into(), r1: vec![bg.0.clone()], r2: vec![bg.1.clone()] });
    truth.insert("HOST_ONLY_0".into(), vec![]);
    println!("host backgrounds: {:.0} s", t0.elapsed().as_secs_f64());
    // Blanks: host traces only; host traces + the reagent bacterium.
    let mut r = Rng(77);
    for (name, kit) in [("BLANK_1", 0u64), ("BLANK_2", 50)] {
        let p = (w(&format!("{name}_R1.fq")), w(&format!("{name}_R2.fq")));
        let mut pw = PairWriter::new(&p.0, &p.1, &format!("{name}:"));
        pw.put_n("host", &host, 20_000, None, &mut r);
        pw.put_n("kitome", &pae, kit, None, &mut r);
        pw.finish();
        samples.push(ScreenSample { name: name.into(), r1: vec![p.0], r2: vec![p.1] });
        truth.insert(name.into(), if kit > 0 { vec!["Pseudomonas aeruginosa PAO1".into()] } else { vec![] });
    }
    // The spike series: 6 levels x 5 replicates, every organism at the
    // level (Poisson counts), the process control at 100 reads, and the
    // reagent bacterium at ~50 reads.
    let levels = [1.0f64, 3.0, 10.0, 30.0, 100.0, 1000.0];
    for &l in &levels {
        for rep in 0..5 {
            let name = format!("SPIKE_{l}rpm_r{rep}");
            let p = (w(&format!("{name}_R1.fq")), w(&format!("{name}_R2.fq")));
            let mut pw = PairWriter::new(&p.0, &p.1, &format!("{name}:"));
            for (o, g) in spiked.iter().zip(&genomes) {
                let n = r.poisson(l * depth as f64 / 1e6);
                pw.put_n(o.name, g, n, o.region, &mut r);
                if !o.name.contains("not in the database") {
                    spikes.push(SpikeIn { sample: Some(name.clone()), organism: o.name.to_string(), amount: l, unit: "rpm".into() });
                }
            }
            pw.put_n("lambda", &lam, 100, None, &mut r);
            pw.put_n("kitome", &pae, r.poisson(50.0), None, &mut r);
            pw.finish();
            samples.push(ScreenSample { name: name.clone(), r1: vec![bg.0.clone(), p.0], r2: vec![bg.1.clone(), p.1] });
            spikes.push(SpikeIn { sample: Some(name.clone()), organism: "Escherichia phage lambda (process control)".into(), amount: 100.0, unit: "reads".into() });
            truth.insert(name, spiked.iter().filter(|o| !o.name.contains("not in")).map(|o| o.name.to_string()).collect());
        }
    }
    // The integrated fragment: adenovirus nt 1-4,344 (the left end, as the
    // E1 region integrated in HEK293 cells), without the rest of the genome.
    for &l in &[3.0f64, 10.0, 30.0, 100.0, 1000.0] {
        for rep in 0..3 {
            let name = format!("FRAGMENT_{l}rpm_r{rep}");
            let p = (w(&format!("{name}_R1.fq")), w(&format!("{name}_R2.fq")));
            let mut pw = PairWriter::new(&p.0, &p.1, &format!("{name}:"));
            pw.put_n("adenovirus_E1", &adv, r.poisson(l), Some((0, 0, 4344)), &mut r);
            pw.finish();
            samples.push(ScreenSample { name: name.clone(), r1: vec![bg.0.clone(), p.0], r2: vec![bg.1.clone(), p.1] });
            truth.insert(name, vec![]);
        }
    }
    println!("{} samples written: {:.0} s", samples.len(), t0.elapsed().as_secs_f64());
    let o = ScreenOptions {
        negative_controls: vec!["BLANK_1".into(), "BLANK_2".into()],
        spike_ins: spikes,
        output_tsv: Some(w("spike_report.tsv")),
        output_json: Some(w("spike_report.json")),
        output_html: Some(w("spike_report.html")),
        ..Default::default()
    };
    let t1 = Instant::now();
    let rep = contaminant_screen(&samples, &dbp, &o).unwrap();
    let secs = t1.elapsed().as_secs_f64();
    let frags: u64 = rep.samples.iter().map(|s| s.fragments).sum();
    println!(
        "SCREEN: {} samples, {} pairs in {secs:.0} s ({:.2} M pairs/s), peak working set {:.0} MB",
        rep.samples.len(),
        frags,
        frags as f64 / secs / 1e6,
        peak_memory_mb()
    );
    // Sensitivity by level, LoD.
    for l in &rep.lod {
        let lv: Vec<String> = l.levels.iter().map(|v| format!("{}:{}/{} ({} reads, rec {})", v.amount, v.detected, v.replicates, short_f(v.mean_reads), v.mean_recovery.map_or("-".into(), |x| format!("{:.2}", x)))).collect();
        println!("LOD {} [{}]: lowest detected {:?}, LoD95 {:?} | {}", l.organism, l.unit, l.lowest_detected, l.lod, lv.join("  "));
    }
    // False positives: an agent detected that the sample does not hold.
    let mut fp = Vec::new();
    let mut bgcalls = (0, 0);
    for s in &rep.samples {
        let want = truth.get(&s.name).cloned().unwrap_or_default();
        for x in &s.organisms {
            if x.call == Call::Detected && x.kind.is_agent() && !want.contains(&x.name) {
                fp.push(format!("{}: {} ({} unique, {} windows, rpm {:.1})", s.name, x.name, x.unique_reads, x.windows, x.rpm));
            }
            if x.name == "Pseudomonas aeruginosa PAO1" && s.name.starts_with("SPIKE") {
                bgcalls.1 += 1;
                if x.call == Call::Background {
                    bgcalls.0 += 1;
                }
            }
        }
    }
    println!("FALSE POSITIVES (agent detected, not in the sample): {}", fp.len());
    for f in &fp {
        println!("  {f}");
    }
    println!("REAGENT BACKGROUND: P. aeruginosa called background in {}/{} spiked samples", bgcalls.0, bgcalls.1);
    // The relative not in the database, and the integrated fragment.
    let mut zika = (0u64, 0u64, 0usize);
    for s in rep.samples.iter().filter(|s| s.name.starts_with("SPIKE_1000")) {
        if let Some(x) = s.organisms.iter().find(|x| x.name.starts_with("Zika virus Natal")) {
            zika.0 += x.unique_reads;
            zika.1 += x.weak_reads;
            zika.2 += usize::from(x.call == Call::Detected);
        }
    }
    println!(
        "RELATIVE NOT IN THE DATABASE (Zika MR 766 at 1000 RPM, 5 replicates): Zika Natal unique reads {}, weak reads {}, detected in {} of 5",
        zika.0, zika.1, zika.2
    );
    for s in rep.samples.iter().filter(|s| s.name.starts_with("FRAGMENT")) {
        let x = s.organisms.iter().find(|x| x.name == "Human mastadenovirus C");
        println!(
            "FRAGMENT {}: {}",
            s.name,
            x.map_or("no read".to_string(), |x| format!("{} ({} unique reads, {} of {} windows, expected {:.1}, dispersion {:?})", x.call.name(), x.unique_reads, x.windows, x.genome_windows, x.expected_windows, x.dispersion.map(|d| (d * 100.0).round() / 100.0)))
        );
    }
    // E. coli K-12 vs the O157 Sakai of the database; S. aureus etc.
    let mut sakai = (0usize, 0u64, 0u64);
    for s in rep.samples.iter().filter(|s| s.name.starts_with("SPIKE_1000")) {
        if let Some(x) = s.organisms.iter().find(|x| x.name.contains("O157")) {
            sakai.0 += usize::from(x.call == Call::Detected);
            sakai.1 += x.unique_reads;
            sakai.2 += x.tied_reads;
        }
    }
    println!("RELATIVE IN THE DATABASE (E. coli K-12 spiked; O157:H7 Sakai also indexed), 1000 RPM x5: Sakai detected {} times, unique {} / tied {}", sakai.0, sakai.1, sakai.2);
    let ctrl: Vec<f64> = rep.spikes.iter().filter(|x| x.organism.contains("lambda")).filter_map(|x| x.recovery).collect();
    let mean = ctrl.iter().sum::<f64>() / ctrl.len().max(1) as f64;
    let min = ctrl.iter().cloned().fold(f64::MAX, f64::min);
    println!("PROCESS CONTROL lambda (100 reads): recovery mean {:.3}, min {:.3} over {} samples", mean, min, ctrl.len());
    for s in rep.samples.iter().filter(|s| s.name.starts_with("HOST_ONLY") || s.name.starts_with("BLANK")) {
        println!(
            "{}: {} pairs, host {:.3}%, given to a reference {}, weak {}, detected {:?}",
            s.name,
            s.fragments,
            100.0 * s.host as f64 / s.fragments as f64,
            s.unique + s.tied,
            s.weak,
            s.detected
        );
        for x in s.organisms.iter().take(5) {
            println!("    {} {}: {} unique, {} tied, {} windows", x.call.name(), x.name, x.unique_reads, x.tied_reads, x.windows);
        }
    }
}

fn short_f(x: f64) -> String {
    if x >= 100.0 {
        format!("{x:.0}")
    } else {
        format!("{x:.1}")
    }
}

/// Theoretical composition (genomic DNA %) of ZymoBIOMICS D6300 (EVEN) and
/// D6310 (LOG).
const ZYMO: [(&str, f64, f64); 10] = [
    ("Pseudomonas aeruginosa PAO1", 12.0, 8.9),
    ("Escherichia coli str. K-12 substr. MG1655", 12.0, 0.089),
    ("Salmonella enterica subsp. enterica serovar Typhimurium str. LT2", 12.0, 0.089),
    ("Limosilactobacillus fermentum", 12.0, 0.0089),
    ("Enterococcus faecalis V583", 12.0, 0.00089),
    ("Staphylococcus aureus subsp. aureus", 12.0, 0.000089),
    ("Listeria monocytogenes EGD-e", 12.0, 89.1),
    ("Bacillus subtilis subsp. subtilis str. 168", 12.0, 0.89),
    ("Saccharomyces cerevisiae S288C", 2.0, 0.89),
    ("Cryptococcus neoformans", 2.0, 0.00089),
];

fn zymo_member(name: &str) -> Option<usize> {
    ZYMO.iter().position(|z| name.starts_with(z.0))
}

fn real_zymo(data: &str, myco: &str, lambda: &str, work: &str) {
    std::fs::create_dir_all(work).unwrap();
    let w = |n: &str| rd(work, n);
    let refs = real_refs(data, myco, lambda, &[]);
    let dbp = w("panel_all.mmgscreen");
    if !Path::new(&dbp).exists() {
        let t0 = Instant::now();
        let s = build_screen_database(&refs, &dbp, &ScreenDbOptions::default()).unwrap();
        println!("panel database (all): {} references, {:.0} s, peak {:.0} MB, file {:.0} MB", s.references.len(), t0.elapsed().as_secs_f64(), peak_memory_mb(), s.file_bytes as f64 / 1e6);
    }
    let mc = |n: &str| rd(data, &format!("metagenome/mock_community/{n}"));
    let samples = vec![
        ScreenSample { name: "ZYMO_EVEN".into(), r1: vec![mc("ERR2984773_zymo_even_D6300_first1M_R1.fastq.gz")], r2: vec![mc("ERR2984773_zymo_even_D6300_first1M_R2.fastq.gz")] },
        ScreenSample { name: "ZYMO_LOG".into(), r1: vec![mc("ERR2935805_zymo_log_D6310_first1M_R1.fastq.gz")], r2: vec![mc("ERR2935805_zymo_log_D6310_first1M_R2.fastq.gz")] },
    ];
    let t1 = Instant::now();
    let o = ScreenOptions { output_tsv: Some(w("zymo_report.tsv")), output_json: Some(w("zymo_report.json")), output_html: Some(w("zymo_report.html")), ..Default::default() };
    let rep = contaminant_screen(&samples, &dbp, &o).unwrap();
    println!("SCREEN zymo: {:.0} s, peak {:.0} MB", t1.elapsed().as_secs_f64(), peak_memory_mb());
    for (si, s) in rep.samples.iter().enumerate() {
        let classified: f64 = s.organisms.iter().filter(|x| zymo_member(&x.name).is_some() || x.call != Call::NotDetected).map(|x| x.reads).sum();
        let all_reads: f64 = s.organisms.iter().map(|x| x.reads).sum();
        println!(
            "\n{}: {} pairs; host {}, unique {}, tied {}, weak {}, unclassified {}, no k-mers {}; detected {} agents",
            s.name,
            s.fragments,
            s.host,
            s.unique,
            s.tied,
            s.weak,
            s.unclassified,
            s.no_kmers,
            s.detected.len()
        );
        println!("  {:<62} {:>9} {:>9} {:>9} {:>8} {:>8}  call", "organism", "theory%", "reads%", "unique", "windows", "breadth");
        let _ = classified;
        for (zi, z) in ZYMO.iter().enumerate() {
            let theory = if si == 0 { z.1 } else { z.2 };
            let xs: Vec<&OrganismResult> = s.organisms.iter().filter(|x| x.name.starts_with(z.0)).collect();
            if xs.is_empty() {
                println!("  {:<62} {:>9} {:>9} {:>9} {:>8} {:>8}  {}", z.0, theory, 0, 0, 0, 0, "no read");
            }
            for x in xs {
                println!(
                    "  {:<62} {:>9} {:>9.4} {:>9} {:>8} {:>8.4}  {} (tied {}, won {}, weak {}){}",
                    x.name,
                    theory,
                    100.0 * x.reads / all_reads.max(1.0),
                    x.unique_reads,
                    x.windows,
                    x.breadth,
                    x.call.name(),
                    x.tied_reads,
                    x.won_reads,
                    x.weak_reads,
                    x.reasons.first().map_or(String::new(), |r| format!(" [{}]", r.chars().take(90).collect::<String>()))
                );
            }
            let _ = zi;
        }
        println!("  others with evidence:");
        for x in s.organisms.iter().filter(|x| zymo_member(&x.name).is_none()) {
            println!(
                "  {:<62} {:>9} {:>9.4} {:>9} {:>8} {:>8.4}  {} {}",
                x.name,
                "-",
                100.0 * x.reads / all_reads.max(1.0),
                x.unique_reads,
                x.windows,
                x.breadth,
                x.call.name(),
                x.reasons.first().map_or(String::new(), |r| format!("({})", r.chars().take(110).collect::<String>()))
            );
        }
    }
    // The same reads through microbio::classify_reads (sketch database of
    // the same references, scaled = 10, the host as decoy).
    let skp = w("panel_all_scaled10.sig");
    if !Path::new(&skp).exists() {
        let t0 = Instant::now();
        let gi: Vec<microbio::GenomeInput> = refs.iter().map(|r| microbio::GenomeInput { files: r.files.clone(), name: r.name.clone(), host: r.kind == AgentKind::Host }).collect();
        microbio::sketch_references(&gi, &skp, &microbio::SketchOptions { k: 31, scaled: 10 }).unwrap();
        println!("\nsketch database scaled 10: {:.0} s", t0.elapsed().as_secs_f64());
    }
    for s in &samples {
        let t2 = Instant::now();
        let c = microbio::classify_reads(&s.r1[0], Some(&s.r2[0]), &[skp.clone()], &microbio::ClassifyOptions::default()).unwrap();
        println!(
            "\nclassify_reads {} ({:.0} s): {} pairs, host {}, classified {} (ambiguous {}), unclassified {}; references with >= 1 read: {}, >= 10: {}, >= 100: {}",
            s.name,
            t2.elapsed().as_secs_f64(),
            c.reads,
            c.host,
            c.classified,
            c.ambiguous,
            c.unclassified,
            c.rows.iter().filter(|x| !x.host && x.reads >= 1.0).count(),
            c.rows.iter().filter(|x| !x.host && x.reads >= 10.0).count(),
            c.rows.iter().filter(|x| !x.host && x.reads >= 100.0).count()
        );
        let mine = rep.samples.iter().find(|x| x.name == s.name).unwrap();
        println!("  {:<62} {:>12} {:>12} {:>12}", "reference", "classify", "screen", "screen call");
        for row in c.rows.iter().filter(|x| !x.host) {
            let m = mine.organisms.iter().find(|x| x.name == row.name);
            println!(
                "  {:<62} {:>12.1} {:>12.1} {:>12}",
                row.name,
                row.reads,
                m.map_or(0.0, |x| x.reads),
                m.map_or("-", |x| x.call.name())
            );
        }
    }
    println!("PEAK working set: {:.0} MB", peak_memory_mb());
}

/// Builds the database of a demo folder written by contam_demo.py
/// (references.tsv) and screens its two samples.
fn demo(folder: &str) {
    let f = |n: &str| rd(folder, n);
    let tsv = std::fs::read_to_string(f("references.tsv")).unwrap();
    let refs: Vec<ScreenReference> = tsv
        .lines()
        .skip(1)
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            ScreenReference { files: vec![f(c[2])], name: c[0].to_string(), kind: AgentKind::from_text(c[1]).unwrap() }
        })
        .collect();
    let dbp = f("demo.mmgscreen");
    let s = build_screen_database(&refs, &dbp, &ScreenDbOptions::default()).unwrap();
    print!("{}", s.resumo());
    let samples = vec![
        ScreenSample { name: "raw_material".into(), r1: vec![f("raw_material_R1.fq.gz")], r2: vec![f("raw_material_R2.fq.gz")] },
        ScreenSample { name: "reagent_blank".into(), r1: vec![f("reagent_blank_R1.fq.gz")], r2: vec![f("reagent_blank_R2.fq.gz")] },
    ];
    let o = ScreenOptions {
        negative_controls: vec!["reagent_blank".into()],
        spike_table: Some(f("spike_ins.tsv")),
        output_html: Some(f("screen.html")),
        output_tsv: Some(f("screen.tsv")),
        output_json: Some(f("screen.json")),
        ..Default::default()
    };
    let rep = contaminant_screen(&samples, &dbp, &o).unwrap();
    print!("{}", rep.resumo());
}

/// A synthetic genome of `gb` gigabases (records of 100 Mb, random, GC
/// 41%) indexed as host with `memory_mb`: the memory bound at scale.
fn host_scale(work: &str, gb: f64, memory_mb: u64) {
    std::fs::create_dir_all(work).unwrap();
    let fa = rd(work, &format!("synthetic_{gb}gb.fa"));
    if !Path::new(&fa).exists() {
        let t0 = Instant::now();
        let mut w = std::io::BufWriter::with_capacity(16 << 20, std::fs::File::create(&fa).unwrap());
        let mut r = Rng(31337);
        let n = (gb * 1e9 / 1e8).round() as usize;
        for i in 0..n {
            writeln!(w, ">synthetic_chr{} 100 Mb", i + 1).unwrap();
            let s = random_dna(100_000_000, 0.41, &mut r);
            for c in s.chunks(100) {
                w.write_all(c).unwrap();
                w.write_all(b"\n").unwrap();
            }
        }
        w.flush().unwrap();
        println!("wrote {fa}: {:.0} s", t0.elapsed().as_secs_f64());
    }
    let tiny = rd(work, "tiny_virus.fa");
    let mut r = Rng(5);
    write_fasta(&tiny, &[("tiny".into(), random_dna(5000, 0.5, &mut r))]);
    let refs = vec![
        ScreenReference { files: vec![fa.clone()], name: "synthetic host".into(), kind: AgentKind::Host },
        ScreenReference { files: vec![tiny], name: "tiny virus".into(), kind: AgentKind::Virus },
    ];
    let t0 = Instant::now();
    let s = build_screen_database(&refs, &rd(work, "scale.mmgscreen"), &ScreenDbOptions { memory_mb, ..Default::default() }).unwrap();
    println!(
        "HOST SCALE {gb} Gb, memory_mb {memory_mb}: {} passes, {} host minimizers, host set {:.0} MB, file {:.0} MB, {:.0} s, peak working set {:.0} MB",
        s.host_passes,
        s.host_minimizers,
        s.host_index_bytes as f64 / 1e6,
        s.file_bytes as f64 / 1e6,
        t0.elapsed().as_secs_f64(),
        peak_memory_mb()
    );
    let t1 = Instant::now();
    let db = load_screen_database(&rd(work, "scale.mmgscreen")).unwrap();
    println!("loaded (with its SHA-256) in {:.1} s; peak working set {:.0} MB", t1.elapsed().as_secs_f64(), peak_memory_mb());
    let _ = std::fs::remove_file(rd(work, "scale.mmgscreen"));
    drop(db);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("demo") => return demo(&args[2]),
        Some("hostpasses") => {
            let t0 = Instant::now();
            let (h, p) = host_set_with_budget(&[args[2].clone()], 31, 31, args[3].parse().unwrap()).unwrap();
            println!("{} minimizers, {p} passes, {:.1} s", h.len(), t0.elapsed().as_secs_f64());
            return;
        }
        Some("screen") => {
            // screen <database> <name> <r1> <r2>: one sample, timed.
            let s = ScreenSample { name: args[3].clone(), r1: vec![args[4].clone()], r2: vec![args[5].clone()] };
            let tl = Instant::now();
            let db = load_screen_database(&args[2]).unwrap();
            println!("database loaded (with its SHA-256) in {:.1} s", tl.elapsed().as_secs_f64());
            let tc = Instant::now();
            let mut rd = microbio::SeqReader::open(&args[4]).unwrap();
            let mut rd2 = microbio::SeqReader::open(&args[5]).unwrap();
            let (mut a, mut b) = (microbio::SeqRecord::default(), microbio::SeqRecord::default());
            let mut pairs = Vec::new();
            while rd.next(&mut a).unwrap() && rd2.next(&mut b).unwrap() && pairs.len() < 200_000 {
                pairs.push((a.seq.clone(), b.seq.clone()));
            }
            println!("read 200,000 pairs in {:.1} s", tc.elapsed().as_secs_f64());
            let o = ScreenOptions::default();
            let tk = Instant::now();
            let n: usize = {
                use rayon::prelude::*;
                pairs.par_iter().filter(|p| classify_one(&db, &p.0, Some(&p.1), &o).outcome == Outcome::Host).count()
            };
            println!("classified 200,000 pairs ({n} host) in {:.2} s on {} threads", tk.elapsed().as_secs_f64(), rayon::current_num_threads());
            drop(db);
            let t0 = Instant::now();
            let rep = contaminant_screen(&[s], &args[2], &ScreenOptions::default()).unwrap();
            let el = t0.elapsed().as_secs_f64();
            let f = rep.samples[0].fragments;
            println!("{} pairs in {el:.1} s: {:.0} pairs/s; peak working set {:.0} MB", f, f as f64 / el, peak_memory_mb());
            print!("{}", rep.resumo().lines().skip(3).take(12).collect::<Vec<_>>().join("\n"));
            println!();
            return;
        }
        Some("host-scale") => return host_scale(&args[2], args[3].parse().unwrap(), args[4].parse().unwrap()),
        Some("host-index") => return real_host_index(&args[2], &args[3]),
        Some("spike") => return real_spike(&args[2], &args[3], &args[4], &args[5]),
        Some("zymo") => return real_zymo(&args[2], &args[3], &args[4], &args[5]),
        _ => {}
    }
    let dir = std::env::temp_dir().join(format!("mmgen_contaminantes_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut t = T { ok: 0, total: 0, failed: Vec::new(), dir: dir.clone() };
    let t0 = Instant::now();
    trial_sha(&mut t);
    trial_minimizers(&mut t);
    let c = community(&t);
    let db = trial_database(&mut t, &c);
    trial_read_calls(&mut t, &c, &db);
    trial_screen(&mut t, &c, &db);
    trial_lod(&mut t, &c, &db);
    let _ = std::fs::remove_dir_all(&dir);
    println!("({:.1} s)", t0.elapsed().as_secs_f64());
    println!("contaminantes tests: {} passed, {} failed", t.ok, t.total - t.ok);
    if t.ok != t.total {
        std::process::exit(1);
    }
}
