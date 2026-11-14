#![allow(dead_code)]
// Trials of src/amplicon.rs: 16S amplicon sequence variants (pair merging,
// filtering, UNOISE3, chimeras), the naive Bayesian classifier, diversity,
// and the tumor mutational burden.
//
// Every result is compared with something written here in the plainest
// way: the merging against a scorer that tries every diagonal with the
// published formulas; UNOISE3 against the rule applied with a full edit
// distance; chimeras against planted two-parent chimeras and an ungapped
// brute-force model; the classifier against word probabilities counted
// naively from the training strings, and against held-out sequences of a
// simulated taxonomy; diversity against hand-computed values; PCoA against
// known coordinates; PERMANOVA against the sums of squares written out; the
// TMB against an oracle that builds each coding sequence and translates it.
//
// Usage: cargo run --release --bin amplicon
//        cargo run --release --bin amplicon -- real <research_data> <work dir> [vsearch executable]
//        cargo run --release --bin amplicon -- holdout <research_data> <work dir>
//        cargo run --release --bin amplicon -- tmb <SEQC2 folder> <RefSeq GFF> <work dir>
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
#[path = "../../../src/referencia.rs"]
#[allow(dead_code)]
mod referencia;
#[path = "../../../src/amplicon.rs"]
mod amplicon;

use amplicon::*;
use std::collections::{HashMap, HashSet};
use std::io::{Read, Write};
use std::path::PathBuf;
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

    fn p(&self, name: &str) -> String {
        self.dir.join(name).to_string_lossy().to_string()
    }

    fn write(&self, name: &str, text: &str) -> String {
        let p = self.p(name);
        std::fs::write(&p, text).unwrap();
        p
    }

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
                let s = format!("{v:?}");
                self.check(false, format!("{what}: should have failed, gave {}", &s[..s.len().min(200)]));
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
//  Random DNA, reads
// ---------------------------------------------------------------------

fn lcg(s: &mut u64) -> u64 {
    *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    *s >> 33
}

fn rnd(s: &mut u64, n: usize) -> usize {
    (lcg(s) % n as u64) as usize
}

fn unit(s: &mut u64) -> f64 {
    (lcg(s) as f64) / (1u64 << 31) as f64
}

fn random_dna(n: usize, s: &mut u64) -> Vec<u8> {
    (0..n).map(|_| b"ACGT"[rnd(s, 4)]).collect()
}

fn rc(t: &[u8]) -> Vec<u8> {
    t.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

fn other_base(b: u8, s: &mut u64) -> u8 {
    loop {
        let c = b"ACGT"[rnd(s, 4)];
        if c != b {
            return c;
        }
    }
}

/// `k` substitutions at distinct random positions.
fn mutate(seq: &[u8], k: usize, s: &mut u64) -> Vec<u8> {
    let mut o = seq.to_vec();
    let mut done = HashSet::new();
    while done.len() < k {
        let p = rnd(s, o.len());
        if done.insert(p) {
            o[p] = other_base(o[p], s);
        }
    }
    o
}

fn perr(q: u8) -> f64 {
    10f64.powf(-(q as f64) / 10.0)
}

/// A read with errors drawn from its qualities.
fn sequence_read(template: &[u8], quals: &[u8], s: &mut u64) -> Vec<u8> {
    template.iter().zip(quals).map(|(&b, &q)| if unit(s) < perr(q) { other_base(b, s) } else { b }).collect()
}

/// Qualities of a MiSeq-like read: high, sliding down at the end (more for
/// read 2).
fn read_quals(len: usize, second: bool, s: &mut u64) -> Vec<u8> {
    (0..len)
        .map(|i| {
            let tail = if second { 22.0 } else { 10.0 };
            let base = 37.0 - tail * (i as f64 / len as f64).powi(3);
            (base + (rnd(s, 5) as f64) - 2.0).clamp(2.0, 41.0) as u8
        })
        .collect()
}

/// A simulated pair of `read_len` from `fragment` (read-through into the
/// adapter when the fragment is shorter than the read).
fn simulate_pair(fragment: &[u8], read_len: usize, s: &mut u64) -> (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>) {
    let ad1 = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGATCGGAAGAGCACACGTCTG";
    let ad2 = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCGGAAGAGCGTCGTGTAGG";
    let mut t1 = fragment.to_vec();
    let mut t2 = rc(fragment);
    while t1.len() < read_len {
        t1.extend_from_slice(ad1);
        t2.extend_from_slice(ad2);
    }
    let q1 = read_quals(read_len, false, s);
    let q2 = read_quals(read_len, true, s);
    let r1 = sequence_read(&t1[..read_len], &q1, s);
    let r2 = sequence_read(&t2[..read_len], &q2, s);
    (r1, q1, r2, q2)
}

fn fastq_text(reads: &[(String, Vec<u8>, Vec<u8>)]) -> String {
    let mut o = String::new();
    for (name, seq, q) in reads {
        o.push('@');
        o.push_str(name);
        o.push('\n');
        o.push_str(std::str::from_utf8(seq).unwrap());
        o.push_str("\n+\n");
        o.extend(q.iter().map(|&x| char::from(x + 33)));
        o.push('\n');
    }
    o
}

// ---------------------------------------------------------------------
//  Oracles
// ---------------------------------------------------------------------

fn full_edit(a: &[u8], b: &[u8]) -> usize {
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    for i in 1..=a.len() {
        let mut cur = vec![i; b.len() + 1];
        for j in 1..=b.len() {
            cur[j] = (prev[j - 1] + usize::from(a[i - 1] != b[j - 1])).min(prev[j] + 1).min(cur[j - 1] + 1);
        }
        prev = cur;
    }
    prev[b.len()]
}

/// Every diagonal scored with the published formulas; the best one
/// (highest score, then longer overlap). Returns (score, d, overlap, diffs).
fn oracle_merge(s1: &[u8], q1: &[u8], s2: &[u8], q2: &[u8], min_ov: usize) -> Option<(f64, i64, usize, usize)> {
    // the formulas, evaluated once per pair of qualities
    thread_local! {
        static SCORES: Vec<(f64, f64)> = (0..94 * 94)
            .map(|k| {
                let (p1, p2) = (perr((k / 94) as u8).min(0.75), perr((k % 94) as u8).min(0.75));
                let same = (1.0 - p1) * (1.0 - p2) + p1 * p2 / 3.0;
                ((same / 0.25f64).log2(), ((1.0 - same) / 0.75f64).log2())
            })
            .collect();
    }
    let a = s1.to_vec();
    let b = rc(s2);
    let qb: Vec<u8> = q2.iter().rev().copied().collect();
    let (l1, l2) = (a.len() as i64, b.len() as i64);
    let mut best: Option<(f64, i64, usize, usize)> = None;
    SCORES.with(|sc| {
    for d in -(l2 - 1)..l1 {
        let mut score = 0.0;
        let mut diffs = 0;
        let mut ov = 0;
        for i in 0..l1 {
            let j = i - d;
            if j < 0 || j >= l2 {
                continue;
            }
            ov += 1;
            let (x, y) = (a[i as usize], b[j as usize]);
            if x == b'N' || y == b'N' {
                continue;
            }
            let (agree, differ) = sc[q1[i as usize] as usize * 94 + qb[j as usize] as usize];
            if x == y {
                score += agree;
            } else {
                score += differ;
                diffs += 1;
            }
        }
        if ov < min_ov {
            continue;
        }
        let better = match best {
            None => true,
            Some((bs, _, bo, _)) => score > bs + 1e-9 || ((score - bs).abs() <= 1e-9 && ov > bo),
        };
        if better {
            best = Some((score, d, ov, diffs));
        }
    }
    });
    best
}

/// UNOISE3 written plainly: full edit distance, the skew inequality itself.
fn oracle_unoise(seqs: &[Vec<u8>], abund: &[u64], alpha: f64, min_size: u64) -> Vec<Option<usize>> {
    let n = seqs.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&x, &y| abund[y].cmp(&abund[x]).then_with(|| seqs[x].cmp(&seqs[y])));
    let mut cents: Vec<usize> = Vec::new();
    let mut out = vec![None; n];
    let find = |q: usize, cents: &[usize]| -> Option<usize> {
        let mut best: Option<(usize, usize)> = None;
        for (ci, &c) in cents.iter().enumerate() {
            let d = full_edit(&seqs[q], &seqs[c]);
            let skew = abund[q] as f64 / abund[c] as f64;
            if skew <= 1.0 / 2f64.powf(alpha * d as f64 + 1.0) * (1.0 + 1e-9) && best.map_or(true, |(bd, _)| d < bd) {
                best = Some((d, ci));
            }
        }
        best.map(|(_, ci)| cents[ci])
    };
    for &q in &order {
        if abund[q] < min_size {
            continue;
        }
        match find(q, &cents) {
            Some(c) => out[q] = Some(c),
            None => {
                cents.push(q);
                out[q] = Some(q);
            }
        }
    }
    for &q in &order {
        if abund[q] < min_size {
            out[q] = find(q, &cents);
        }
    }
    out
}

/// Ungapped two-parent model, brute force over (A != B, x): the fewest
/// differences, then the most votes on the weaker side (votes = the other
/// parent's differences on a side, net of the model's). Returns (model,
/// min votes, closest single parent) -- Hamming throughout: with gaps
/// costing 3, sequences that differ by substitutions align without gaps.
fn oracle_model(q: &[u8], parents: &[&[u8]], max_model: usize) -> Option<(usize, i64, usize)> {
    let ham = |a: &[u8], b: &[u8]| a.iter().zip(b).filter(|(x, y)| x != y).count();
    let single = parents.iter().map(|p| ham(q, p)).min().unwrap();
    let mut best: Option<(usize, i64)> = None;
    for (i, a) in parents.iter().enumerate() {
        for (j, b) in parents.iter().enumerate() {
            if i == j {
                continue;
            }
            for x in 0..=q.len() {
                let (la, lb) = (ham(&q[..x], &a[..x]), ham(&q[..x], &b[..x]));
                let (ra, rb) = (ham(&q[x..], &a[x..]), ham(&q[x..], &b[x..]));
                let model = la + rb;
                if model > max_model {
                    continue;
                }
                let votes = (lb as i64 - la as i64).min(ra as i64 - rb as i64);
                if best.map_or(true, |(m, v)| (model, -votes) < (m, -v)) {
                    best = Some((model, votes));
                }
            }
        }
    }
    best.map(|(m, v)| (m, v, single))
}

// ---------------------------------------------------------------------
//  1. Pair merging
// ---------------------------------------------------------------------

fn trials_merge(t: &mut T) {
    println!("pair merging");
    let mut s = 11u64;
    let opt = MergeOptions::default();
    // the posterior qualities against the formulas
    let mut ok = true;
    for q1 in [2u8, 5, 10, 20, 30, 38, 41] {
        for q2 in [2u8, 7, 15, 25, 35, 41] {
            let (p1, p2) = (perr(q1).min(0.75), perr(q2).min(0.75));
            let pa = (p1 * p2 / 3.0) / (1.0 - p1 - p2 + 4.0 * p1 * p2 / 3.0);
            let (lo, hi) = if p1 <= p2 { (p1, p2) } else { (p2, p1) };
            let pd = lo * (1.0 - hi / 3.0) / (lo + hi - 4.0 * lo * hi / 3.0);
            let qa = (-10.0 * pa.log10()).round().clamp(0.0, 93.0) as u8;
            let qd = (-10.0 * pd.log10()).round().clamp(0.0, 93.0) as u8;
            if merged_quality_agree(q1, q2) != qa || merged_quality_differ(q1.max(q2), q1.min(q2)) != qd {
                ok = false;
                println!("    q {q1} {q2}: agree {} vs {qa}, differ {} vs {qd}", merged_quality_agree(q1, q2), merged_quality_differ(q1.max(q2), q1.min(q2)));
            }
        }
    }
    t.check(ok, "posterior qualities of merged bases = Edgar & Flyvbjerg's formulas");
    t.check(merged_quality_agree(40, 40) > 41 && merged_quality_differ(40, 10) < 40 && merged_quality_differ(40, 40) <= 4, "agreeing bases gain quality, disagreeing ones lose it");

    // clean pairs of several fragment lengths
    for &(len, rl) in &[(253usize, 250usize), (300, 250), (400, 250), (480, 250), (253, 150), (390, 200)] {
        let f = random_dna(len, &mut s);
        let r1 = f[..rl].to_vec();
        let r2 = rc(&f)[..rl].to_vec();
        let q = vec![38u8; rl];
        let m = merge_pair(&r1, &q, &r2, &q, &opt);
        let expect_ov = 2 * rl - len;
        match m {
            Ok(m) => {
                t.check(m.seq == f && m.overlap == expect_ov && m.diffs == 0 && !m.staggered, format!("fragment {len} from 2x{rl}: merged back, overlap {expect_ov} (got {})", m.overlap));
                t.check(m.qual.iter().all(|&x| x == 38 || x == 41), format!("fragment {len}: qualities kept outside, 41 (capped) in the overlap"));
            }
            Err(e) => {
                t.check(false, format!("fragment {len} from 2x{rl}: {e:?}"));
            }
        }
    }
    // too short to overlap
    let f = random_dna(520, &mut s);
    let (r1, r2) = (f[..250].to_vec(), rc(&f)[..250].to_vec());
    let q = vec![38u8; 250];
    t.check(merge_pair(&r1, &q, &r2, &q, &opt) == Err(MergeFail::NoAlignment), "fragment 520 from 2x250 (no overlap): no alignment");
    // overlap of 15, under and over min_overlap
    let f = random_dna(485, &mut s);
    let (r1, r2) = (f[..250].to_vec(), rc(&f)[..250].to_vec());
    t.check(merge_pair(&r1, &q, &r2, &q, &opt).map(|m| m.seq) == Ok(f.clone()), "overlap of 15 bases: merged");
    t.check(merge_pair(&r1, &q, &r2, &q, &MergeOptions { min_overlap: 20, ..Default::default() }).is_err(), "overlap of 15 with min_overlap 20: rejected");

    // planted mismatches in the overlap: the better quality wins
    let f = random_dna(253, &mut s);
    let mut r1 = f[..250].to_vec();
    let mut r2 = rc(&f)[..250].to_vec();
    let mut q1 = vec![36u8; 250];
    let mut q2 = vec![36u8; 250];
    // read 1 wrong at 100 with Q8; read 2 wrong at fragment position 150 with Q9
    r1[100] = other_base(r1[100], &mut s);
    q1[100] = 8;
    let j = 252 - 150; // read-2 index of fragment position 150
    r2[j] = other_base(r2[j], &mut s);
    q2[j] = 9;
    match merge_pair(&r1, &q1, &r2, &q2, &opt) {
        Ok(m) => {
            t.check(m.seq == f, "planted mismatches: the higher-quality base is kept on both sides");
            t.check(m.diffs == 2, format!("planted mismatches: 2 differences counted (got {})", m.diffs));
            t.check(m.qual[100] == merged_quality_differ(36, 8) && m.qual[150] == merged_quality_differ(36, 9), "planted mismatches: the disagreement quality");
        }
        Err(e) => {
            t.check(false, format!("planted mismatches: {e:?}"));
        }
    }
    // 11 high-quality mismatches: too many
    let mut r1b = f[..250].to_vec();
    for k in 0..11 {
        let p = 20 + 15 * k;
        r1b[p] = other_base(r1b[p], &mut s);
    }
    let qq = vec![36u8; 250];
    t.check(merge_pair(&r1b, &qq, &rc(&f)[..250], &qq, &opt) == Err(MergeFail::TooManyDiffs), "11 mismatches in the overlap: too many differences");
    t.check(merge_pair(&r1b, &qq, &rc(&f)[..250], &qq, &MergeOptions { max_diffs: 11, ..Default::default() }).is_ok(), "11 mismatches with max_diffs 11: merged");
    t.check(merge_pair(&r1b, &qq, &rc(&f)[..250], &qq, &MergeOptions { max_diffs: 20, max_diff_pct: 4.0, ..Default::default() }) == Err(MergeFail::TooManyDiffs), "11 mismatches in 247 (4.5%) with max_diff_pct 4: too many");
    // N in the overlap: the other read's base
    let mut r1n = f[..250].to_vec();
    r1n[120] = b'N';
    let m = merge_pair(&r1n, &qq, &rc(&f)[..250], &qq, &opt).unwrap();
    t.check(m.seq == f && m.diffs == 0 && m.qual[120] == 36, "an N in read 1: the base of read 2, with its quality, not a difference");

    // staggered: fragment shorter than the reads
    let f = random_dna(200, &mut s);
    let (r1, q1, r2, q2) = {
        let ad1 = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGATCGGAAGAGCACACGTCTG";
        let ad2 = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCGGAAGAGCGTCGTGTAGG";
        let mut a = f.clone();
        a.extend_from_slice(&ad1[..50]);
        let mut b = rc(&f);
        b.extend_from_slice(&ad2[..50]);
        (a, vec![36u8; 250], b, vec![36u8; 250])
    };
    t.check(merge_pair(&r1, &q1, &r2, &q2, &opt) == Err(MergeFail::Staggered), "fragment 200 in 2x250 (read-through): staggered, rejected");
    let m = merge_pair(&r1, &q1, &r2, &q2, &MergeOptions { allow_stagger: true, ..Default::default() });
    t.check(m.as_ref().map(|m| (m.seq.clone(), m.staggered)) == Ok((f.clone(), true)), "with allow_stagger: the fragment, adapters cut");

    // truncation
    let f = random_dna(253, &mut s);
    let r1 = f[..250].to_vec();
    let r2 = rc(&f)[..250].to_vec();
    let mut q2 = vec![36u8; 250];
    for x in q2.iter_mut().skip(160) {
        *x = 2;
    }
    let m = merge_pair(&r1, &vec![36u8; 250], &r2, &q2, &MergeOptions { trunc_qual: Some(2), ..Default::default() }).unwrap();
    t.check(m.overlap == 250 + 160 - 253 && m.seq == f, "trunc_qual 2 cuts read 2 at its first Q2 base: overlap 157");
    let m = merge_pair(&r1, &vec![36u8; 250], &r2, &vec![36u8; 250], &MergeOptions { trunc_len_r1: 240, trunc_len_r2: 160, ..Default::default() }).unwrap();
    t.check(m.overlap == 240 + 160 - 253 && m.seq == f, "trunc_len 240/160 (as the DADA2 tutorial): overlap 147");
    t.check(merge_pair(&r1[..8], &[36u8; 8], &r2, &vec![36u8; 250], &opt) == Err(MergeFail::TooShort), "a read shorter than min_overlap: too short");

    // unrelated reads almost never merge
    let mut merged = 0;
    for _ in 0..2000 {
        let a = random_dna(250, &mut s);
        let b = random_dna(250, &mut s);
        let qa = read_quals(250, false, &mut s);
        let qb = read_quals(250, true, &mut s);
        if merge_pair(&a, &qa, &b, &qb, &opt).is_ok() {
            merged += 1;
        }
    }
    t.check(merged <= 2, format!("2000 unrelated pairs: {merged} merged (at most 2 allowed)"));

    // the brute-force oracle on simulated pairs of every kind
    let mut agree = 0;
    let mut cases = 0;
    let mut first_bad = String::new();
    for it in 0..3000 {
        let len = 160 + rnd(&mut s, 360); // 160..520
        let rl = [150usize, 200, 250][rnd(&mut s, 3)];
        let f = random_dna(len, &mut s);
        let (mut a, mut qa, mut b, mut qb) = simulate_pair(&f, rl, &mut s);
        if it % 7 == 0 {
            // some N
            let p = rnd(&mut s, rl);
            a[p] = b'N';
        }
        if it % 5 == 0 {
            // different read lengths
            let cut = rl - rnd(&mut s, 60);
            b.truncate(cut);
            qb.truncate(cut);
        }
        if it % 11 == 0 {
            // a poor read
            for x in qa.iter_mut().skip(100) {
                *x = 3;
            }
            for x in a.iter_mut().skip(100) {
                if rnd(&mut s, 3) == 0 {
                    *x = other_base(*x, &mut s);
                }
            }
        }
        cases += 1;
        let got = merge_pair(&a, &qa, &b, &qb, &opt);
        let want = oracle_merge(&a, &qa, &b, &qb, opt.min_overlap);
        let l1 = a.len() as i64;
        let l2 = b.len() as i64;
        let fine = match (&got, want) {
            (Err(MergeFail::NoAlignment), None) => true,
            (Err(MergeFail::NoAlignment), Some((sc, ..))) => sc < opt.min_score,
            (_, None) => false,
            (Ok(m), Some((sc, d, ov, df))) => sc >= opt.min_score && m.overlap == ov && m.diffs == df && (m.seq.len() as i64) == d + l2 && (m.score - sc).abs() < 1e-6,
            (Err(MergeFail::TooManyDiffs), Some((sc, _, _, df))) => sc >= opt.min_score && df > opt.max_diffs,
            (Err(MergeFail::Staggered), Some((sc, d, _, df))) => sc >= opt.min_score && df <= opt.max_diffs && (d < 0 || d + l2 < l1),
            (Err(MergeFail::TooShort), _) => false,
        };
        if fine {
            agree += 1;
        } else if first_bad.is_empty() {
            first_bad = format!("case {it}: len {len} reads {}x{}: got {:?}, oracle {:?}", a.len(), b.len(), got.as_ref().map(|m| (m.overlap, m.diffs, m.score)), want);
        }
    }
    t.check(agree == cases, format!("{agree} of {cases} simulated pairs: same overlap, differences and score as the brute-force scorer of every diagonal. {first_bad}"));

    // expected errors
    t.check(near(expected_errors(&[20; 10]), 0.1, 1e-12), "expected errors of 10 x Q20 = 0.1");
    t.check(near(expected_errors(&[10, 30, 40]), 0.1 + 0.001 + 0.0001, 1e-12), "expected errors of Q10 + Q30 + Q40 = 0.1011");
}

// ---------------------------------------------------------------------
//  2. Primers, edit distance
// ---------------------------------------------------------------------

fn trials_primers(t: &mut T) {
    println!("primers and distances");
    let mut s = 23u64;
    let f515 = Primer::new("GTGYCAGCMGCCGCGGTAA", 0.1).unwrap();
    let r806 = Primer::new("GGACTACNVGGGTWTCTAAT", 0.1).unwrap();
    t.check(f515.max_mismatches == 1 && r806.max_mismatches == 2, "10% of 19 and 20 bases: 1 and 2 mismatches allowed");
    let fwd_inst = b"GTGCCAGCAGCCGCGGTAA".to_vec(); // Y=C, M=A
    let rev_inst = b"GGACTACAAGGGTATCTAAT".to_vec(); // N=A V=A W=A
    for spacer in [0usize, 3, 7] {
        let insert = random_dna(253, &mut s);
        let mut read = random_dna(spacer, &mut s);
        read.extend_from_slice(&fwd_inst);
        read.extend_from_slice(&insert);
        read.extend_from_slice(&rc(&rev_inst));
        read.extend_from_slice(&random_dna(spacer, &mut s));
        let e = find_forward_primer(&read, &f515, 10);
        let b = find_reverse_primer(&read, &r806, 10);
        t.check(e == Some(spacer + 19) && b == Some(spacer + 19 + 253), format!("primers after a {spacer}-base spacer: found at their places ({e:?}, {b:?})"));
        let mut one = read.clone();
        one[spacer + 5] = other_base(one[spacer + 5], &mut s);
        t.check(find_forward_primer(&one, &f515, 10) == Some(spacer + 19), "forward primer with 1 mismatch: found");
        one[spacer + 12] = other_base(one[spacer + 12], &mut s);
        t.check(find_forward_primer(&one, &f515, 10).is_none(), "forward primer with 2 mismatches: not found (1 allowed)");
    }
    t.check(find_forward_primer(&random_dna(300, &mut s), &f515, 10).is_none(), "no primer in random DNA");
    t.fails("primer with a wrong letter", Primer::new("GTGYCAGCXGCC", 0.1), &["X", "IUPAC", "position 9"]);

    // bounded edit distance against the full one
    let mut ok = 0;
    let n = 3000;
    for _ in 0..n {
        let a = random_dna(30 + rnd(&mut s, 60), &mut s);
        let mut b = a.clone();
        for _ in 0..rnd(&mut s, 8) {
            match rnd(&mut s, 3) {
                0 => {
                    let p = rnd(&mut s, b.len());
                    b[p] = other_base(b[p], &mut s);
                }
                1 => {
                    let p = rnd(&mut s, b.len() + 1);
                    b.insert(p, b"ACGT"[rnd(&mut s, 4)]);
                }
                _ => {
                    if b.len() > 1 {
                        let p = rnd(&mut s, b.len());
                        b.remove(p);
                    }
                }
            }
        }
        let k = rnd(&mut s, 7);
        let full = full_edit(&a, &b);
        let got = bounded_edit_distance(&a, &b, k);
        if got == if full <= k { Some(full) } else { None } {
            ok += 1;
        }
    }
    t.check(ok == n, format!("bounded edit distance = full edit distance when <= k, None otherwise: {ok} of {n}"));
    // alignment profiles: substitutions exactly where they are; a deletion
    // costs one difference
    let (mut ok_sub, mut ok_del) = (0, 0);
    for _ in 0..500 {
        let a = random_dna(80 + rnd(&mut s, 40), &mut s);
        let b = mutate(&a, rnd(&mut s, 8), &mut s);
        let prof = diff_profile(&a, &b);
        if prof.iter().zip(a.iter().zip(&b)).all(|(&p, (x, y))| p as usize == usize::from(x != y)) {
            ok_sub += 1;
        }
        let k = rnd(&mut s, 5);
        let mut c = mutate(&a, k, &mut s);
        c.remove(10 + rnd(&mut s, c.len() - 20));
        if diff_profile(&a, &c).iter().map(|&x| x as usize).sum::<usize>() <= k + 1 {
            ok_del += 1;
        }
    }
    t.check(ok_sub == 500, format!("alignment profiles of substitutions: a difference at each substituted position and nowhere else ({ok_sub} of 500)"));
    t.check(ok_del == 500, format!("alignment profiles with a deletion and k substitutions: at most k + 1 differences ({ok_del} of 500)"));
}

// ---------------------------------------------------------------------
//  3. UNOISE3 and chimeras
// ---------------------------------------------------------------------

fn trials_unoise(t: &mut T) {
    println!("UNOISE3 and chimeras");
    let mut s = 31u64;
    // small random families against the plain oracle
    let mut agree = 0;
    let rounds = 40;
    for _ in 0..rounds {
        let mut seqs: Vec<Vec<u8>> = Vec::new();
        let mut abund: Vec<u64> = Vec::new();
        let roots = 1 + rnd(&mut s, 4);
        for _ in 0..roots {
            let root = random_dna(60, &mut s);
            let a = 50 + rnd(&mut s, 3000) as u64;
            seqs.push(root.clone());
            abund.push(a);
            for _ in 0..rnd(&mut s, 15) {
                let mut v = mutate(&root, 1 + rnd(&mut s, 4), &mut s);
                if rnd(&mut s, 5) == 0 {
                    v.remove(rnd(&mut s, 60));
                }
                if !seqs.contains(&v) {
                    seqs.push(v);
                    abund.push(1 + rnd(&mut s, (a / 4) as usize) as u64);
                }
            }
        }
        let opt = DenoiseOptions::default();
        let got = unoise(&seqs, &abund, &opt);
        let mapped: Vec<Option<usize>> = got.assigned.iter().map(|a| a.map(|ci| got.centroids[ci])).collect();
        let want = oracle_unoise(&seqs, &abund, opt.alpha, opt.min_size);
        if mapped == want {
            agree += 1;
        } else {
            println!("    unoise differs: got {mapped:?} want {want:?} abund {abund:?}");
        }
    }
    t.check(agree == rounds, format!("UNOISE3 = the rule applied with a full edit distance: {agree} of {rounds} random families"));
    t.check(max_error_distance(80, 10, 2.0) == 1 && max_error_distance(79, 10, 2.0) == 0 && max_error_distance(320, 10, 2.0) == 2, "skew limits: 8x allows 1 difference, 32x allows 2 (boundaries included)");

    // a community with reads of 0.5% errors: the true sequences come back
    let base = random_dna(253, &mut s);
    let truth: Vec<Vec<u8>> = vec![
        random_dna(253, &mut s),
        mutate(&base, 30, &mut s),
        base.clone(),
        mutate(&base, 1, &mut s), // 1 difference, at 1/4 of base's reads: kept
        mutate(&base, 2, &mut s), // 2 differences, at 1/10: kept
        random_dna(253, &mut s),
        random_dna(253, &mut s),
        random_dna(253, &mut s),
    ];
    let reads_of = [6000u64, 3000, 4000, 1000, 400, 800, 300, 200];
    let mut derep: HashMap<Vec<u8>, u64> = HashMap::new();
    for (seq, &n) in truth.iter().zip(&reads_of) {
        for _ in 0..n {
            let r: Vec<u8> = seq.iter().map(|&b| if unit(&mut s) < 0.005 { other_base(b, &mut s) } else { b }).collect();
            *derep.entry(r).or_insert(0) += 1;
        }
    }
    let mut seqs: Vec<Vec<u8>> = derep.keys().cloned().collect();
    seqs.sort();
    let abund: Vec<u64> = seqs.iter().map(|x| derep[x]).collect();
    let den = unoise(&seqs, &abund, &DenoiseOptions::default());
    let got: HashSet<Vec<u8>> = den.centroids.iter().map(|&c| seqs[c].clone()).collect();
    let want: HashSet<Vec<u8>> = truth.iter().cloned().collect();
    t.check(got == want, format!("simulated community (8 sequences, two at 1 and 2 differences, 0.5% errors): {} centroids, {} true, {} spurious", got.len(), got.intersection(&want).count(), got.difference(&want).count()));
    let assigned: u64 = den.assigned.iter().zip(&abund).filter(|(a, _)| a.is_some()).map(|(_, &n)| n).sum();
    let total: u64 = abund.iter().sum();
    t.check(assigned as f64 / total as f64 > 0.97, format!("{:.1}% of the reads assigned to an ASV (map_rare)", 100.0 * assigned as f64 / total as f64));
    // (the reads carry substitutions only: Hamming distance to the truth)
    let ham = |a: &[u8], b: &[u8]| a.iter().zip(b).filter(|(x, y)| x != y).count();
    let wrong: u64 = den
        .assigned
        .iter()
        .enumerate()
        .filter_map(|(u, a)| a.map(|ci| (u, den.centroids[ci])))
        .filter(|&(u, c)| {
            let best = truth.iter().map(|x| ham(&seqs[u], x)).min().unwrap();
            ham(&seqs[u], &seqs[c]) != best
        })
        .map(|(u, _)| abund[u])
        .sum();
    t.check(wrong == 0, format!("every read goes to its closest true sequence ({wrong} reads elsewhere)"));

    // chimeras
    let p1 = random_dna(253, &mut s);
    let p2 = mutate(&p1, 40, &mut s);
    let p3 = random_dna(253, &mut s);
    let mut c1 = p1[..100].to_vec();
    c1.extend_from_slice(&p2[100..]);
    let mut c2 = p2[..170].to_vec();
    c2.extend_from_slice(&p3[170..]);
    let real = mutate(&p1, 3, &mut s); // a real variant, not a chimera
    let mut c_hi = p1[..120].to_vec(); // a chimera too abundant to call
    c_hi.extend_from_slice(&p3[120..]);
    let seqs: Vec<&[u8]> = vec![&p1, &p2, &p3, &c_hi, &real, &c1, &c2];
    let abund = vec![5000u64, 4000, 3000, 500, 120, 100, 90];
    let ch = find_chimeras(&seqs, &abund, &ChimeraOptions::default());
    t.check(ch[0].is_none() && ch[1].is_none() && ch[2].is_none(), "parents are not chimeras");
    match &ch[5] {
        Some(m) => {
            let ok_bp = (0..=100).rev().find(|&x| x < 100 && p1[x] != p2[x]).map_or(0, |x| x + 1) <= m.breakpoint && m.breakpoint <= (100..253).find(|&x| p1[x] != p2[x]).unwrap();
            t.check(m.parent_left == 0 && m.parent_right == 1 && m.model_diffs == 0 && ok_bp, format!("chimera p1|p2 at 100: parents 0 and 1, breakpoint {} between the diagnostic positions", m.breakpoint));
        }
        None => {
            t.check(false, "chimera p1|p2 at 100: not found");
        }
    }
    t.check(ch[6].as_ref().map(|m| (m.parent_left, m.parent_right)) == Some((1, 2)), format!("chimera p2|p3 at 170: parents 1 and 2 ({:?})", ch[6]));
    t.check(ch[4].is_none(), "a real variant with 3 private differences is not a chimera");
    t.check(ch[3].is_none(), "a chimera with 1/6 of its parent's reads is not called (abskew 16)");
    // brute force on substitution-only families
    let mut ok = 0;
    let rounds = 60;
    for _ in 0..rounds {
        let a = random_dna(120, &mut s);
        let b = mutate(&a, 10 + rnd(&mut s, 20), &mut s);
        let c = mutate(&a, 10 + rnd(&mut s, 20), &mut s);
        let x = 20 + rnd(&mut s, 80);
        let mut q = if rnd(&mut s, 2) == 0 { a[..x].to_vec() } else { b[..x].to_vec() };
        q.extend_from_slice(&c[x..]);
        let q = if rnd(&mut s, 3) == 0 { mutate(&q, 1, &mut s) } else { q };
        let seqs: Vec<&[u8]> = vec![&a, &b, &c, &q];
        let abund = vec![1000u64, 900, 800, 10];
        let need = 1 + rnd(&mut s, 3);
        let got = find_chimeras(&seqs, &abund, &ChimeraOptions { max_model_diffs: 1, min_diffs: need, ..Default::default() });
        let want = oracle_model(&q, &[&a, &b, &c], 1);
        let expect = want.map_or(false, |(model, votes, single)| votes >= need as i64 && single > model);
        let fine = match (&got[3], want) {
            (Some(m), Some((model, votes, single))) => expect && m.model_diffs == model && m.closest_parent_diffs == single && m.votes_left.min(m.votes_right) as i64 == votes,
            (None, _) => !expect,
            _ => false,
        };
        if fine {
            ok += 1;
        } else {
            println!("    model: got {:?}, oracle {want:?}", got[3]);
        }
    }
    t.check(ok == rounds, format!("two-parent models = the ungapped brute force: {ok} of {rounds}"));
}

// ---------------------------------------------------------------------
//  4. The whole pipeline on simulated paired FASTQ
// ---------------------------------------------------------------------

fn trials_pipeline(t: &mut T) {
    println!("the ASV pipeline on simulated samples");
    let mut s = 47u64;
    let base = random_dna(253, &mut s);
    let truth: Vec<Vec<u8>> = vec![
        random_dna(253, &mut s),
        base.clone(),
        mutate(&base, 1, &mut s),
        random_dna(253, &mut s),
        random_dna(252, &mut s),
        random_dna(254, &mut s),
        mutate(&base, 25, &mut s),
        random_dna(253, &mut s),
    ];
    // a chimera of truth[3] and truth[7]
    let mut chim = truth[3][..130].to_vec();
    chim.extend_from_slice(&truth[7][130..]);
    // reads per sample per sequence (the last column is the chimera)
    // (the chimera has 60 reads in all: its parents have 2,000 and 2,100, over abskew 16)
    let plan: [[u64; 9]; 3] = [[1500, 1200, 400, 900, 300, 0, 200, 600, 20], [200, 2000, 600, 100, 900, 400, 0, 800, 25], [800, 0, 0, 1000, 600, 300, 500, 700, 15]];
    let mut samples = Vec::new();
    let mut sheet = String::from("sample\tforward\treverse\n");
    for (si, row) in plan.iter().enumerate() {
        let mut r1s = Vec::new();
        let mut r2s = Vec::new();
        let mut k = 0;
        for (ti, &n) in row.iter().enumerate() {
            let f = if ti == 8 { &chim } else { &truth[ti] };
            for _ in 0..n {
                let (a, qa, b, qb) = simulate_pair(f, 200, &mut s);
                r1s.push((format!("r{k} 1:N:0:1"), a, qa));
                r2s.push((format!("r{k} 2:N:0:1"), b, qb));
                k += 1;
            }
        }
        // shuffle the pairs together
        for i in (1..r1s.len()).rev() {
            let j = rnd(&mut s, i + 1);
            r1s.swap(i, j);
            r2s.swap(i, j);
        }
        let n1 = format!("S{}_S1_L001_R1_001.fastq.gz", si + 1);
        let n2 = format!("S{}_S1_L001_R2_001.fastq.gz", si + 1);
        let p1 = t.write_packed(&n1, &fastq_text(&r1s));
        let p2 = t.write_packed(&n2, &fastq_text(&r2s));
        sheet.push_str(&format!("S{}\t{n1}\t{n2}\n", si + 1));
        samples.push(AmpliconSample { name: format!("S{}", si + 1), r1: p1, r2: Some(p2) });
    }
    let sheet_path = t.write("sheet.tsv", &sheet);
    let from_sheet = read_sample_sheet(&sheet_path);
    t.check(from_sheet.as_ref().map(|v| v == &samples).unwrap_or(false), "sample sheet with a header and relative paths: the same samples");
    let files: Vec<String> = samples.iter().flat_map(|x| vec![x.r1.clone(), x.r2.clone().unwrap()]).collect();
    t.check(pair_fastq_files(&files).map(|v| v == samples).unwrap_or(false), "Illumina file names paired into samples S1, S2, S3");

    let fa = t.p("asvs.fasta");
    let tab = t.p("table.tsv");
    let trk = t.p("track.tsv");
    let t0 = Instant::now();
    let r = amplicon_asvs(&samples, &fa, &tab, Some(&trk), &AsvOptions::default());
    let Ok(r) = r else {
        t.check(false, format!("pipeline: {:?}", r.err()));
        return;
    };
    println!("    {:.2} s; {} ASVs; {} chimeras", t0.elapsed().as_secs_f64(), r.asvs.len(), r.chimeras.len());
    let got: HashSet<Vec<u8>> = r.asvs.iter().map(|a| a.sequence.as_bytes().to_vec()).collect();
    let want: HashSet<Vec<u8>> = truth.iter().cloned().collect();
    t.check(got == want, format!("8 true sequences recovered exactly ({} ASVs, {} true)", got.len(), got.intersection(&want).count()));
    t.check(r.chimeras.len() == 1 && r.chimeras[0].sequence.as_bytes() == chim.as_slice(), "the planted chimera is found and removed");
    if let Some(c) = r.chimeras.first() {
        let name_of = |seq: &[u8]| r.asvs.iter().find(|a| a.sequence.as_bytes() == seq).map(|a| a.name.clone()).unwrap_or_default();
        t.check(c.parent_left == name_of(&truth[3]) && c.parent_right == name_of(&truth[7]), format!("its parents are the two planted ones ({} and {})", c.parent_left, c.parent_right));
    }
    // names by reads
    let names_ok = r.asvs.windows(2).all(|w| w[0].reads >= w[1].reads) && r.asvs.iter().enumerate().all(|(i, a)| a.name == format!("ASV_{}", i + 1));
    t.check(names_ok, "ASV_1, ASV_2, ... by decreasing reads");
    // tracking
    let mono = r.samples.iter().all(|x| x.input >= x.merged && x.merged >= x.filtered && x.filtered >= x.denoised && x.denoised >= x.nonchimeric);
    t.check(mono, "tracking: input >= merged >= filtered >= denoised >= non-chimeric");
    let inputs: Vec<u64> = plan.iter().map(|row| row.iter().sum()).collect();
    t.check(r.samples.iter().zip(&inputs).all(|(x, &n)| x.input == n), "tracking: the input of each sample is its number of pairs");
    let kept: f64 = r.samples.iter().map(|x| x.nonchimeric as f64 / x.input as f64).sum::<f64>() / 3.0;
    t.check(kept > 0.6, format!("{:.1}% of the pairs reach the table", 100.0 * kept));
    let track_text = read_text(&trk);
    t.check(track_text.lines().count() == 4 && track_text.starts_with("sample\tinput\tmerged\tfiltered\tdenoised\tnonchimeric"), "tracking table written: header + 3 samples");
    // the table, re-read
    let ft = read_feature_table(&tab).unwrap();
    t.check(ft.samples == vec!["S1", "S2", "S3"] && ft.features.len() == r.asvs.len() && ft.integer, "feature table re-read: 3 samples, one row per ASV, counts");
    let reads_ok = r.asvs.iter().enumerate().all(|(i, a)| ft.counts[i].iter().sum::<f64>() as u64 == a.reads);
    t.check(reads_ok, "table rows add up to each ASV's reads");
    let fasta = amplicon::read_fasta(&fa).unwrap();
    t.check(fasta.len() == r.asvs.len() && fasta.iter().zip(&r.asvs).all(|((id, sq), a)| *id == a.name && sq == a.sequence.as_bytes()), "ASV FASTA = the ASVs, in order");
    // composition: relative abundance in each sample against the plan
    let mut worst = 0.0f64;
    for (si, row) in plan.iter().enumerate() {
        let planted: f64 = row[..8].iter().sum::<u64>() as f64;
        let col: f64 = ft.counts.iter().map(|c| c[si]).sum();
        for (ti, tseq) in truth.iter().enumerate() {
            let want = row[ti] as f64 / planted;
            let gotc = r.asvs.iter().position(|a| a.sequence.as_bytes() == tseq.as_slice()).map_or(0.0, |i| ft.counts[i][si]) / col;
            worst = worst.max((want - gotc).abs());
        }
    }
    t.check(worst < 0.03, format!("relative abundances within 3 points of the plan (worst {:.3})", worst));
    t.check(r.json().starts_with("{\"asv_fasta\":") && r.resumo().contains("ASVs from 3 samples"), "json and resumo");

    // with primers in the reads
    let fwd = b"GTGCCAGCAGCCGCGGTAA";
    let rev = b"GGACTACAAGGGTATCTAAT";
    let mut r1s = Vec::new();
    let mut r2s = Vec::new();
    for k in 0..1500 {
        let ti = [0usize, 1, 3][k % 3];
        let mut f = fwd.to_vec();
        f.extend_from_slice(&truth[ti]);
        f.extend_from_slice(&rc(rev));
        let (a, qa, b, qb) = simulate_pair(&f, 250, &mut s);
        r1s.push((format!("p{k}"), a, qa));
        r2s.push((format!("p{k}"), b, qb));
    }
    let pr = vec![AmpliconSample { name: "P".into(), r1: t.write("p_R1.fq", &fastq_text(&r1s)), r2: Some(t.write("p_R2.fq", &fastq_text(&r2s))) }];
    let opt = AsvOptions { primer_forward: Some("GTGYCAGCMGCCGCGGTAA".into()), primer_reverse: Some("GGACTACNVGGGTWTCTAAT".into()), ..Default::default() };
    match amplicon_asvs(&pr, &t.p("p.fa"), &t.p("p.tsv"), None, &opt) {
        Ok(rp) => {
            let got: HashSet<Vec<u8>> = rp.asvs.iter().map(|a| a.sequence.as_bytes().to_vec()).collect();
            let want: HashSet<Vec<u8>> = [0usize, 1, 3].iter().map(|&i| truth[i].clone()).collect();
            t.check(got == want, "with 515F/806R in the reads: the ASVs are the inserts, primers removed");
            t.check(rp.samples[0].primer_missing < 30, format!("primers found in almost every read ({} missing)", rp.samples[0].primer_missing));
        }
        Err(e) => {
            t.check(false, format!("primers: {e}"));
        }
    }
    // without primers in the reads but asked to remove them: a warning
    let r2 = amplicon_asvs(&samples[..1], &t.p("w.fa"), &t.p("w.tsv"), None, &opt).unwrap();
    t.check(r2.warnings.iter().any(|w| w.contains("primers were missing")), "primers asked but absent: a warning that says so");

    // errors
    let bad_names = vec![AmpliconSample { name: "X".into(), r1: samples[0].r1.clone(), r2: samples[1].r2.clone() }];
    t.fails("mates from different samples", amplicon_asvs(&bad_names, &t.p("x.fa"), &t.p("x.tsv"), None, &AsvOptions::default()), &["named", "same reads in the same order"]);
    let short2 = t.write("short_R2.fq", &fastq_text(&r2s[..10]));
    let sh = vec![AmpliconSample { name: "Y".into(), r1: t.p("p_R1.fq"), r2: Some(short2) }];
    t.fails("read-2 file shorter", amplicon_asvs(&sh, &t.p("y.fa"), &t.p("y.tsv"), None, &AsvOptions::default()), &["more reads than", "same order"]);
    let missing = vec![AmpliconSample { name: "Z".into(), r1: t.p("nothing_R1.fq"), r2: None }];
    t.fails("missing file", amplicon_asvs(&missing, &t.p("z.fa"), &t.p("z.tsv"), None, &AsvOptions::default()), &["nothing_R1.fq", "does not exist"]);
    t.fails("output = input", amplicon_asvs(&samples[..1], &samples[0].r1, &t.p("z.tsv"), None, &AsvOptions::default()), &["also one of the input"]);
    let fa_in = t.write("reads.fa", ">a\nACGTACGTAC\n");
    let fa_s = vec![AmpliconSample { name: "F".into(), r1: fa_in, r2: None }];
    t.fails("FASTA reads", amplicon_asvs(&fa_s, &t.p("f.fa"), &t.p("f.tsv"), None, &AsvOptions::default()), &["FASTA", "qualities"]);
    let bad_q = t.write("badq.fq", "@a\nACGT\n+\nII5I\n");
    let bq = vec![AmpliconSample { name: "Q".into(), r1: bad_q, r2: None }];
    let o64 = AsvOptions { merge: MergeOptions { phred_offset: 64, ..Default::default() }, ..Default::default() };
    t.fails("Phred+33 read as Phred+64", amplicon_asvs(&bq, &t.p("q.fa"), &t.p("q.tsv"), None, &o64), &["below the Phred offset 64", "phred_offset = 33"]);
    let sheet_bad = t.write("sheet_bad.tsv", "S1\tS1_S1_L001_R1_001.fastq.gz\nS1\tS2_S1_L001_R1_001.fastq.gz\n");
    t.fails("sample sheet with a repeated name", read_sample_sheet(&sheet_bad), &["appears twice", "line 2"]);
}

// ---------------------------------------------------------------------
//  5. The classifier
// ---------------------------------------------------------------------

/// A simulated taxonomy: lineage -> sequences (each genus a mutated copy of
/// its family's ancestor, and so on up).
fn simulated_taxonomy(s: &mut u64, len: usize, per_genus: usize) -> Vec<(Vec<String>, Vec<Vec<u8>>)> {
    let root = random_dna(len, s);
    let mut out = Vec::new();
    for p in 0..2 {
        let ph = mutate(&root, len / 7, s);
        for c in 0..2 {
            let cl = mutate(&ph, len / 9, s);
            for o in 0..2 {
                let or = mutate(&cl, len / 12, s);
                for f in 0..2 {
                    let fa = mutate(&or, len / 16, s);
                    for g in 0..3 {
                        let ge = mutate(&fa, len / 20, s);
                        let lin = vec![
                            "Bacteria".to_string(),
                            format!("Phylum{p}"),
                            format!("Class{p}{c}"),
                            format!("Order{p}{c}{o}"),
                            format!("Family{p}{c}{o}{f}"),
                            format!("Genus{p}{c}{o}{f}{g}"),
                        ];
                        let seqs = (0..per_genus).map(|_| mutate(&ge, len / 100, s)).collect();
                        out.push((lin, seqs));
                    }
                }
            }
        }
    }
    out
}

fn trials_classifier(t: &mut T) {
    println!("the naive Bayesian classifier");
    let mut s = 59u64;
    // lineage parsing
    t.check(parse_lineage(b"Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;") == Some(vec!["Bacteria".into(), "Firmicutes".into(), "Bacilli".into(), "Lactobacillales".into(), "Lactobacillaceae".into(), "Lactobacillus".into()]), "DADA2 lineage header");
    t.check(
        parse_lineage(b"RS_GCF_031457235.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae_C;g__Hydrogenophaga;s__Hydrogenophaga laconesensis [locus_tag=x] [location=5..1531]")
            .map(|l| (l.len(), l[0].clone(), l[6].clone()))
            == Some((7, "Bacteria".into(), "Hydrogenophaga laconesensis".into())),
        "GTDB header: rank prefixes removed, tags cut",
    );
    t.check(parse_lineage(b"Bacteria;Bacillota;Clostridia;Eubacteriales;Clostridiaceae;Clostridium sensu stricto;").map(|l| l[5].clone()) == Some("Clostridium sensu stricto".into()), "a name with spaces in a DADA2 header");
    t.check(parse_lineage(b"seq1 some description").is_none(), "a header without lineage");

    // a tiny reference: word probabilities against naive counting
    let tiny = [
        ("Bacteria;P1;C1;O1;F1;G1;", "ACGTACGGTTACCAGTAGCATTTACGACGGATTACAGGACT"),
        ("Bacteria;P1;C1;O1;F1;G1;", "ACGTACGGTTACCAGTAGCATTAACGACGGATTACAGGACT"),
        ("Bacteria;P1;C1;O1;F1;G2;", "TTGACCAGGTACGATGCAGGCATAGGACCATTAGCAGTTAG"),
        ("Bacteria;P1;C1;O1;F2;G3;", "GGGCATTACAGACCAGATTAGCCATAGACGATTACGGACTA"),
        ("Bacteria;P1;C1;O1;F2;G3;", "GGGCATTACAGACCAGATTCGCCATAGACGATTACGGACTA"),
        ("Bacteria;P1;C1;O1;F2;G3;", "GGGCATTACAGACCAGATTAGCCATAGTCGATTACGGACTA"),
        ("Bacteria;P1;C1;O1;F1;", "ACGTACGGTTACCTGTAGCATTTACGACGGATTACCGGACT"),
    ];
    let mut text = String::new();
    for (l, sq) in &tiny {
        text.push_str(&format!(">{l}\n{sq}\n"));
    }
    let tiny_path = t.write("tiny_ref.fa", &text);
    let m = TaxModel::train(&tiny_path, &TrainOptions::default()).unwrap();
    t.check(m.classes() == 4 && m.n_sequences == 7, format!("tiny reference: 4 lineages (one stops at the family), 7 sequences ({}, {})", m.classes(), m.n_sequences));
    let words = |sq: &str| -> HashSet<String> { (0..=sq.len() - 8).map(|i| sq[i..i + 8].to_string()).collect() };
    let classes: Vec<String> = m.lineages.iter().map(|l| l.join(";") + ";").collect();
    let n = tiny.len() as f64;
    let query = "GGGCATTACAGACCAGATTAGCCATAGACGATTACGGTCTAACGTACGG";
    let mut worst = 0.0f64;
    let mut naive_scores = Vec::new();
    for (ci, cl) in classes.iter().enumerate() {
        let members: Vec<&str> = tiny.iter().filter(|(l, _)| l == cl).map(|(_, x)| *x).collect();
        let mut score = 0.0;
        for w in words(query) {
            let nw = tiny.iter().filter(|(_, x)| words(x).contains(&w)).count() as f64;
            let mw = members.iter().filter(|x| words(x).contains(&w)).count() as f64;
            let pw = (nw + 0.5) / (n + 1.0);
            let lp = ((mw + pw) / (members.len() as f64 + 1.0)).ln();
            score += lp;
            let code = w.bytes().fold(0u32, |c, b| (c << 2) | b"ACGT".iter().position(|&x| x == b).unwrap() as u32);
            worst = worst.max((m.log_prob(code, ci) - lp).abs());
        }
        naive_scores.push(score);
    }
    t.check(worst < 1e-12, format!("ln P(w|G) of every query word and lineage = naive counting (worst {worst:.1e})"));
    let full = m.full_scores(query.as_bytes());
    let worst_s = full.iter().zip(&naive_scores).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    t.check(worst_s < 1e-5, format!("whole-query scores = sum of naive log probabilities (worst {worst_s:.1e})"));
    let a = m.classify(query.as_bytes(), &ClassifyOptions::default());
    let best_naive = (0..classes.len()).max_by(|&x, &y| naive_scores[x].partial_cmp(&naive_scores[y]).unwrap()).unwrap();
    t.check(a.best == m.lineages[best_naive], format!("the assignment is the naive argmax ({:?})", a.best));

    // a simulated taxonomy, one sequence of each genus held out
    let tax = simulated_taxonomy(&mut s, 1400, 6);
    let mut train = String::new();
    let mut held: Vec<(Vec<String>, Vec<u8>)> = Vec::new();
    for (lin, seqs) in &tax {
        for (i, q) in seqs.iter().enumerate() {
            if i == 0 {
                held.push((lin.clone(), q.clone()));
            } else {
                train.push_str(&format!(">{};\n{}\n", lin.join(";"), String::from_utf8_lossy(q)));
            }
        }
    }
    let train_path = t.write_packed("sim_train.fa.gz", &train);
    let model_path = t.p("sim.mmgtax");
    let ts = train_classifier(&train_path, Some(&model_path), &TrainOptions::default()).unwrap();
    t.check(ts.classes == 48 && ts.sequences == 240 && ts.taxa_per_rank == vec![1, 2, 4, 8, 16, 48], format!("training summary: 48 genera, 240 sequences, taxa per rank {:?}", ts.taxa_per_rank));
    let loaded = TaxModel::load(&model_path).unwrap();
    let model = TaxModel::train(&train_path, &TrainOptions::default()).unwrap();
    let opt = ClassifyOptions::default();
    let mut right = 0;
    let mut same = true;
    let mut conf_ok = true;
    for (lin, q) in &held {
        let a = model.classify(q, &opt);
        let b = loaded.classify(q, &opt);
        same &= a == b;
        if a.lineage == *lin {
            right += 1;
        }
        conf_ok &= a.confidence.windows(2).all(|w| w[0] >= w[1]);
    }
    t.check(right >= 46, format!("held-out full-length sequences: {right} of 48 assigned to the right genus"));
    t.check(same, "the saved and reloaded model classifies identically");
    t.check(conf_ok, "confidence never rises from one rank to the next");
    // short reads: a 250-base window, confidence 0.5
    let mut right_short = 0;
    let mut family_ok = 0;
    for (lin, q) in &held {
        let a = model.classify(&q[500..750], &ClassifyOptions { min_confidence: 0.5, ..Default::default() });
        if a.lineage == *lin {
            right_short += 1;
        }
        if a.lineage.len() >= 5 && a.lineage[..5] == lin[..5] {
            family_ok += 1;
        }
    }
    t.check(right_short >= 40 && family_ok >= 46, format!("held-out 250-base windows: {right_short} of 48 to the genus, {family_ok} to the family (min_confidence 0.5)"));
    // reverse complement
    let (lin0, q0) = &held[5];
    let rcq = rc(&q0[300..700]);
    let a = model.classify(&rcq, &ClassifyOptions { try_reverse_complement: true, ..Default::default() });
    t.check(a.reverse_complemented && a.lineage == *lin0, "try_reverse_complement: a reversed read is classified on the right strand");
    let a = model.classify(&rcq, &ClassifyOptions::default());
    t.check(!a.reverse_complemented && a.lineage.len() < 6, "without it, the reversed read gets no genus");
    // training on the region between two primers
    let fwd = b"GTGCCAGCAGCCGCGGTAA";
    let rev = b"GGACTACAAGGGTATCTAAT";
    let mut reg_train = String::new();
    let mut reg_held: Vec<(Vec<String>, Vec<u8>)> = Vec::new();
    for (lin, seqs) in &tax {
        for (i, q) in seqs.iter().enumerate() {
            let mut full = q[..300].to_vec();
            full.extend_from_slice(fwd);
            full.extend_from_slice(&q[300..553]);
            full.extend_from_slice(&rc(rev));
            full.extend_from_slice(&q[553..]);
            if i == 0 {
                reg_held.push((lin.clone(), q[300..553].to_vec()));
            } else {
                reg_train.push_str(&format!(">{};\n{}\n", lin.join(";"), String::from_utf8_lossy(&full)));
            }
        }
    }
    reg_train.push_str(&format!(">{};\n{}\n>{};\n{}\n", tax[0].0.join(";"), String::from_utf8_lossy(&tax[0].1[1]), tax[1].0.join(";"), String::from_utf8_lossy(&tax[1].1[2])));
    let reg_path = t.write("region_train.fa", &reg_train);
    let reg_opt = TrainOptions { primer_forward: Some("GTGYCAGCMGCCGCGGTAA".into()), primer_reverse: Some("GGACTACNVGGGTWTCTAAT".into()), ..Default::default() };
    let rs = train_classifier(&reg_path, Some(&t.p("region.mmgtax")), &reg_opt).unwrap();
    let whole = train_classifier(&reg_path, None, &TrainOptions::default()).unwrap();
    t.check(rs.sequences == 240 && rs.without_region == 2 && rs.source.contains("region between"), format!("region training: 240 sequences, the 2 without primer sites left out ({}, {})", rs.sequences, rs.without_region));
    t.check(rs.pairs * 3 < whole.pairs, format!("the region model stores far fewer pairs ({} against {})", rs.pairs, whole.pairs));
    let rm = TaxModel::load(&t.p("region.mmgtax")).unwrap();
    let right_reg = reg_held.iter().filter(|(lin, q)| rm.classify(q, &ClassifyOptions { min_confidence: 0.5, ..Default::default() }).lineage == *lin).count();
    t.check(right_reg >= 44, format!("held-out regions on the region model: {right_reg} of 48 to the genus"));
    t.fails("one primer only", TaxModel::train(&reg_path, &TrainOptions { primer_forward: Some("GTGYCAGCMGCCGCGGTAA".into()), ..Default::default() }), &["both primers"]);
    t.fails("primers found nowhere", TaxModel::train(&reg_path, &TrainOptions { primer_forward: Some("AAAAAAAAAAAAAAAAAAAAAA".into()), primer_reverse: Some("CCCCCCCCCCCCCCCCCCCCCC".into()), ..Default::default() }), &["primers were found in none"]);
    // determinism: same result from the same seed, whatever the order
    let a1 = model.classify(&held[3].1, &opt);
    let _ = model.classify(&held[9].1, &opt);
    let a2 = model.classify(&held[3].1, &opt);
    t.check(a1 == a2, "the bootstrap depends on the sequence and the seed only");

    // classify_asvs, with exact species
    let mut asvs = String::new();
    let mut sp_ref = String::new();
    for (i, (lin, q)) in held.iter().take(6).enumerate() {
        asvs.push_str(&format!(">ASV_{}\n{}\n", i + 1, String::from_utf8_lossy(&q[400..653])));
        // species reference: the full sequence named by species
        sp_ref.push_str(&format!(">ref{i} {} species{i}\n{}\n", lin[5], String::from_utf8_lossy(q)));
    }
    // ASV_2 also in a sequence of another species of the same genus: ambiguous
    sp_ref.push_str(&format!(">refX {} otherspecies\nAAAA{}TTTT\n", held[1].0[5], String::from_utf8_lossy(&held[1].1[400..653])));
    // ASV_3 also in a sequence of another genus: ignored by the genus check
    sp_ref.push_str(&format!(">refY WrongGenus wrongspecies\n{}\n", String::from_utf8_lossy(&held[2].1[400..653])));
    let asv_path = t.write("cl_asvs.fa", &asvs);
    let sp_path = t.write("species.fa", &sp_ref);
    let out = t.p("taxonomy.tsv");
    let res = classify_asvs(&asv_path, &model_path, &out, &ClassifyOptions { min_confidence: 0.5, species_reference: Some(sp_path.clone()), ..Default::default() }).unwrap();
    let sp: Vec<Option<String>> = res.assignments.iter().map(|a| a.species.clone()).collect();
    t.check(sp[0].as_deref() == Some("species0") && sp[1].is_none() && sp[2].as_deref() == Some("species2"), format!("exact species: unique match named, ambiguous left empty, other genus ignored ({sp:?})"));
    let res2 = classify_asvs(&asv_path, &model_path, &out, &ClassifyOptions { min_confidence: 0.5, species_reference: Some(sp_path), max_species: 2, ..Default::default() }).unwrap();
    t.check(res2.assignments[1].species.as_deref() == Some("otherspecies/species1"), "max_species 2: both names, joined by '/'");
    let tt = read_text(&out);
    let head = tt.lines().next().unwrap_or("");
    t.check(head == "ASV\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tKingdom_confidence\tPhylum_confidence\tClass_confidence\tOrder_confidence\tFamily_confidence\tGenus_confidence", format!("taxonomy table header: {head}"));
    t.check(tt.lines().count() == 7, "taxonomy table: one line per ASV");
    t.check(res.resumo().contains("Genus") && res.json().contains("\"assignments\""), "json and resumo");
    // a FASTA reference given directly (trained on the fly)
    let res3 = classify_asvs(&asv_path, &train_path, &t.p("taxonomy2.tsv"), &ClassifyOptions { min_confidence: 0.5, ..Default::default() }).unwrap();
    t.check(res3.assignments.iter().zip(&res.assignments).all(|(a, b)| a.lineage == b.lineage), "the reference FASTA as the model: the same assignments");
    // errors
    let notax = t.write("notax.fa", ">seq1\nACGTACGTACGTAAAAAA\n>seq2\nACGTTTTTTTACGGGAGGA\n");
    t.fails("reference without lineages", TaxModel::train(&notax, &TrainOptions::default()), &["lineage", "seq1"]);
    let bytes = std::fs::read(&model_path).unwrap();
    let cut = t.p("cut.mmgtax");
    std::fs::write(&cut, &bytes[..bytes.len() / 2]).unwrap();
    t.fails("truncated model", TaxModel::load(&cut), &["cut.mmgtax"]);
    t.fails("min_confidence 2", classify_asvs(&asv_path, &model_path, &t.p("x.tsv"), &ClassifyOptions { min_confidence: 2.0, ..Default::default() }), &["min_confidence"]);
}

// ---------------------------------------------------------------------
//  6. Diversity
// ---------------------------------------------------------------------

fn trials_diversity(t: &mut T) {
    println!("diversity");
    let mut s = 71u64;
    let a = alpha_of(&[10, 5, 3, 1, 1, 0]);
    let p = [0.5f64, 0.25, 0.15, 0.05, 0.05];
    let h: f64 = -p.iter().map(|x| x * x.ln()).sum::<f64>();
    t.check(a.depth == 20 && a.observed == 5 && a.singletons == 2 && a.doubletons == 0, "alpha: depth 20, 5 observed, 2 singletons, 0 doubletons");
    t.check(near(a.chao1, 6.0, 1e-12), format!("Chao1 = 5 + 2*1/(2*(0+1)) = 6 (got {})", a.chao1));
    t.check(near(a.shannon, h, 1e-12) && near(a.simpson, 0.66, 1e-12), "Shannon and Simpson (1 - 0.34) by hand");
    t.check(near(a.pielou.unwrap(), h / 5f64.ln(), 1e-12) && near(a.goods_coverage.unwrap(), 0.9, 1e-12), "Pielou = H / ln 5, Good's coverage = 1 - 2/20");
    let one = alpha_of(&[7, 0, 0]);
    t.check(one.pielou.is_none() && one.shannon == 0.0 && one.simpson == 0.0, "a single feature: Shannon 0, Simpson 0, Pielou undefined");
    t.check(near(alpha_of(&[3, 2, 2, 1]).chao1, 4.0 + 0.0, 1e-12) && near(alpha_of(&[1, 1, 1, 2]).chao1, 4.0 + 3.0 * 2.0 / 4.0, 1e-12), "Chao1 with one singleton (no correction) and with three");
    // ln gamma
    let lg_ok = near(ln_gamma(1.0), 0.0, 1e-13) && near(ln_gamma(2.0), 0.0, 1e-13) && near(ln_gamma(0.5), std::f64::consts::PI.sqrt().ln(), 1e-13) && near(ln_gamma(11.0), 3628800f64.ln(), 1e-11) && near(ln_gamma(1001.0), (1..=1000).map(|k| (k as f64).ln()).sum::<f64>(), 1e-8);
    t.check(lg_ok, "ln Gamma at 1, 2, 0.5, 11, 1001");
    // rarefaction
    t.check(near(expected_observed(&[2, 1], 2), 5.0 / 3.0, 1e-12), "expected observed of {a,a,b} at 2 reads = 5/3 (enumerated)");
    let counts = vec![500u64, 200, 100, 50, 20, 10, 5, 3, 2, 1, 1, 1, 4000];
    t.check(near(expected_observed(&counts, counts.iter().sum()), 13.0, 1e-12) && near(expected_observed(&counts, 1), 1.0, 1e-12), "the curve ends at the observed count and starts at 1");
    // the log-gamma path (both c and n large) against the product
    let big = vec![30_000u64, 20_000, 9_000, 7, 1];
    let n = 12_000u64;
    let total: u64 = big.iter().sum();
    let exact: f64 = big.iter().map(|&c| 1.0 - (0..n).map(|j| ((total - c - j) as f64 / (total - j) as f64).ln()).sum::<f64>().exp()).sum();
    t.check(near(expected_observed(&big, n), exact, 1e-9), "expected observed with large counts and depth (log-gamma) = the product");
    let mut acc = 0.0;
    let draws = 3000;
    let mut rng = Rng::new(5);
    let mut sums_ok = true;
    for _ in 0..draws {
        let r = rarefy(&counts, 300, &mut rng);
        sums_ok &= r.iter().sum::<u64>() == 300 && r.iter().zip(&counts).all(|(x, c)| x <= c);
        acc += r.iter().filter(|&&x| x > 0).count() as f64;
    }
    let mc = acc / draws as f64;
    let ex = expected_observed(&counts, 300);
    t.check(sums_ok, "rarefied samples have exactly the depth, never more of a feature than it had");
    t.check(near(mc, ex, 0.08), format!("mean observed over {draws} rarefactions ({mc:.3}) = the exact expectation ({ex:.3})"));
    t.check(rarefy(&counts, 300, &mut Rng::new(9)) == rarefy(&counts, 300, &mut Rng::new(9)), "same seed, same rarefaction");
    // distances
    t.check(near(bray_curtis(&[1.0, 2.0, 3.0], &[3.0, 2.0, 1.0]), 4.0 / 12.0, 1e-15), "Bray-Curtis by hand: 4/12");
    t.check(near(jaccard(&[1.0, 0.0, 1.0, 1.0], &[0.0, 1.0, 1.0, 1.0]), 0.5, 1e-15), "Jaccard by hand: 1 - 2/4");
    // PCoA of points with known coordinates
    let pts: [(f64, f64); 6] = [(0.0, 0.0), (1.0, 0.0), (0.0, 2.0), (3.0, 1.0), (-1.0, -1.0), (2.0, 2.5)];
    let d: Vec<Vec<f64>> = pts.iter().map(|a| pts.iter().map(|b| ((a.0 - b.0) * (a.0 - b.0) + (a.1 - b.1) * (a.1 - b.1)).sqrt()).collect()).collect();
    let pc = pcoa(&d);
    let mut worst = 0.0f64;
    for i in 0..6 {
        for j in 0..6 {
            let e: f64 = pc.coordinates[i].iter().zip(&pc.coordinates[j]).map(|(x, y)| (x - y) * (x - y)).sum::<f64>().sqrt();
            worst = worst.max((e - d[i][j]).abs());
        }
    }
    t.check(pc.coordinates[0].len() == 2 && worst < 1e-9, format!("PCoA of planar points: 2 axes, distances reproduced (worst {worst:.1e})"));
    let tv: f64 = (0..6).flat_map(|i| (i + 1..6).map(move |j| (i, j))).map(|(i, j)| d[i][j] * d[i][j]).sum::<f64>() / 6.0;
    t.check(near(pc.eigenvalues.iter().sum::<f64>(), tv, 1e-9) && near(pc.proportion_explained.iter().sum::<f64>(), 1.0, 1e-12), "eigenvalues add up to the total variance; proportions to 1");
    // Bray-Curtis matrix of random samples: eigenpairs, sum
    let samples: Vec<Vec<f64>> = (0..12).map(|_| (0..30).map(|_| if rnd(&mut s, 3) == 0 { 0.0 } else { rnd(&mut s, 100) as f64 }).collect()).collect();
    let dbc: Vec<Vec<f64>> = samples.iter().map(|a| samples.iter().map(|b| bray_curtis(a, b)).collect()).collect();
    let pc = pcoa(&dbc);
    let tv: f64 = (0..12).flat_map(|i| (i + 1..12).map(move |j| (i, j))).map(|(i, j)| dbc[i][j] * dbc[i][j]).sum::<f64>() / 12.0;
    t.check(near(pc.eigenvalues.iter().sum::<f64>(), tv, 1e-9), "Bray-Curtis PCoA: eigenvalues (negative ones too) add up to the total variance");
    t.check(pc.eigenvalues.windows(2).all(|w| w[0] >= w[1]), "eigenvalues from the largest");
    let pos: f64 = pc.eigenvalues.iter().filter(|&&x| x > 1e-12).sum();
    let prop_ok = pc.eigenvalues.iter().zip(&pc.proportion_explained).all(|(&e, &p)| if e > 1e-12 { near(p, e / pos, 1e-12) } else { p == 0.0 });
    t.check(prop_ok && near(pc.proportion_explained.iter().sum::<f64>(), 1.0, 1e-12), "proportions over the positive eigenvalues, 0 for the negative ones (scikit-bio's convention)");
    let m: Vec<Vec<f64>> = (0..8).map(|i| (0..8).map(|j| ((i * 7 + j * 3) % 11) as f64 + if i == j { 5.0 } else { 0.0 }).collect()).collect();
    let sym: Vec<Vec<f64>> = (0..8).map(|i| (0..8).map(|j| m[i][j] + m[j][i]).collect()).collect();
    let (vals, vecs) = symmetric_eigen(&sym);
    let mut res = 0.0f64;
    for k in 0..8 {
        for i in 0..8 {
            let av: f64 = (0..8).map(|j| sym[i][j] * vecs[j][k]).sum();
            res = res.max((av - vals[k] * vecs[i][k]).abs());
        }
    }
    t.check(res < 1e-9, format!("Jacobi: A v = lambda v for every eigenpair (worst residual {res:.1e})"));
    // PERMANOVA by hand
    let dm: Vec<Vec<f64>> = vec![
        vec![0.0, 0.2, 0.3, 0.25, 0.7, 0.8, 0.75, 0.9],
        vec![0.2, 0.0, 0.25, 0.3, 0.8, 0.7, 0.85, 0.75],
        vec![0.3, 0.25, 0.0, 0.2, 0.75, 0.9, 0.7, 0.8],
        vec![0.25, 0.3, 0.2, 0.0, 0.9, 0.75, 0.8, 0.7],
        vec![0.7, 0.8, 0.75, 0.9, 0.0, 0.3, 0.2, 0.25],
        vec![0.8, 0.7, 0.9, 0.75, 0.3, 0.0, 0.25, 0.2],
        vec![0.75, 0.85, 0.7, 0.8, 0.2, 0.25, 0.0, 0.3],
        vec![0.9, 0.75, 0.8, 0.7, 0.25, 0.2, 0.3, 0.0],
    ];
    let g = vec![0usize, 0, 0, 0, 1, 1, 1, 1];
    let mut sst = 0.0;
    let mut ssw = 0.0;
    for i in 0..8 {
        for j in i + 1..8 {
            sst += dm[i][j] * dm[i][j];
            if g[i] == g[j] {
                ssw += dm[i][j] * dm[i][j];
            }
        }
    }
    let sst = sst / 8.0;
    let ssw = ssw / 4.0; // both groups have 4
    let f_hand = ((sst - ssw) / 1.0) / (ssw / 6.0);
    let (f, r2, pv) = permanova(&dm, &g, 999, 1);
    t.check(near(f, f_hand, 1e-9) && near(r2, (sst - ssw) / sst, 1e-12), format!("PERMANOVA pseudo-F = {f_hand:.4} and R2 by hand (got {f:.4})"));
    t.check(pv < 0.06, format!("separated groups: p = {pv:.3} (35 distinct splits, so p >= 1/35)"));
    let (_, _, pnull) = permanova(&dm, &[0, 1, 0, 1, 0, 1, 0, 1], 999, 1);
    t.check(pnull > 0.2, format!("mixed groups: p = {pnull:.3}"));

    // the file tools
    let table = "#ASV\tA\tB\tC\tD\ttaxonomy\nASV_1\t50\t40\t0\t5\tx\nASV_2\t30\t35\t10\t0\tx\nASV_3\t0\t5\t60\t70\tx\nASV_4\t20\t20\t30\t25\tx\nASV_5\t0\t0\t0\t0\tx\n";
    let tab = t.write("div_table.tsv", table);
    let ft = read_feature_table(&tab).unwrap();
    t.check(ft.samples == vec!["A", "B", "C", "D"] && ft.features.len() == 5 && ft.integer, "BIOM-style table: #ASV header, taxonomy column dropped");
    let ar = alpha_diversity(&tab, Some(&t.p("alpha.tsv")), &AlphaOptions { curves_output: Some(t.p("curves.tsv")), ..Default::default() }).unwrap();
    t.check(ar.rows.len() == 4 && ar.rows[0].observed == 3 && ar.rows[2].depth == 100, "alpha per sample");
    t.check(ar.curves.len() == 4 && ar.curves.iter().all(|c| c.observed.windows(2).all(|w| w[0] <= w[1] + 1e-12)), "rarefaction curves never go down");
    t.check(read_text(&t.p("curves.tsv")).lines().count() > 4, "curves written");
    let ar2 = alpha_diversity(&tab, None, &AlphaOptions { rarefy_depth: Some(100), ..Default::default() }).unwrap();
    t.check(ar2.rows.iter().all(|r| r.depth == 100) && ar2.dropped.is_empty(), "rarefied to 100: every sample has 100 reads");
    t.fails("rarefied above every sample", alpha_diversity(&tab, None, &AlphaOptions { rarefy_depth: Some(101), ..Default::default() }), &["every sample was left out", "smaller rarefaction depth"]);
    let tab2 = t.write("div_table2.tsv", "ASV\tA\tB\tC\nx\t10\t200\t50\ny\t5\t100\t60\n");
    let ar3 = alpha_diversity(&tab2, None, &AlphaOptions { rarefy_depth: Some(100), ..Default::default() }).unwrap();
    t.check(ar3.dropped == vec!["A"] && ar3.rows.len() == 2, format!("rarefied to 100: the sample with 15 reads left out ({:?})", ar3.dropped));
    let meta = t.write("meta.tsv", "#SampleID\tgroup\tweek\nA\tearly\t1\nB\tearly\t2\nC\tlate\t9\nD\tlate\t10\n");
    let br = beta_diversity(&tab, Some(&t.p("bc.tsv")), &BetaOptions { metadata: Some(meta.clone()), pcoa_output: Some(t.p("pcoa.tsv")), ..Default::default() }).unwrap();
    let hand = bray_curtis(&ft.sample_column(0), &ft.sample_column(2));
    t.check(near(br.distances[0][2], hand, 1e-12) && br.samples.len() == 4, "beta_diversity: the Bray-Curtis of the table");
    t.check(br.permanova.as_ref().map(|p| (p.groups.clone(), p.n)) == Some((vec!["early".into(), "late".into()], 4)), "PERMANOVA groups from the metadata");
    t.check(read_text(&t.p("pcoa.tsv")).lines().nth(1).map_or(false, |l| l.starts_with("proportion_explained")), "PCoA table written with the proportions");
    let bj = beta_diversity(&tab, None, &BetaOptions { metric: "jaccard".into(), ..Default::default() }).unwrap();
    t.check(near(bj.distances[0][1], 1.0 - 3.0 / 4.0, 1e-12), "Jaccard of A and B = 1 - 3/4");
    // the k-mer distance sees relatives
    let mut s2 = 5u64;
    let x = random_dna(253, &mut s2);
    let x1 = mutate(&x, 1, &mut s2);
    let y = random_dna(253, &mut s2);
    let fa = t.write("kmer.fa", &format!(">a\n{}\n>b\n{}\n>c\n{}\n", String::from_utf8_lossy(&x), String::from_utf8_lossy(&x1), String::from_utf8_lossy(&y)));
    let kt = t.write("kmer.tsv", "ASV\tS1\tS2\tS3\na\t10\t0\t0\nb\t0\t10\t0\nc\t0\t0\t10\n");
    let bk = beta_diversity(&kt, None, &BetaOptions { metric: "kmer".into(), asv_fasta: Some(fa), ..Default::default() }).unwrap();
    t.check(bk.distances[0][1] < 0.05 && bk.distances[0][2] > 0.5, format!("k-mer distance: close relatives {:.3}, strangers {:.3} (Bray-Curtis says 1 for both)", bk.distances[0][1], bk.distances[0][2]));
    // composition
    let taxo = t.write("tax.tsv", "ASV\tKingdom\tPhylum\tGenus\tKingdom_confidence\nASV_1\tBacteria\tFirmicutes\tLactobacillus\t1\nASV_2\tBacteria\tFirmicutes\t\t1\nASV_3\tBacteria\tBacteroidota\tBacteroides\t1\nASV_4\t\t\t\t0\n");
    let cr = taxa_composition(&tab, &taxo, Some(&t.p("comp.tsv")), &CompositionOptions { rank: "genus".into(), top: 0 }).unwrap();
    t.check(cr.taxa.contains(&"Firmicutes (unclassified)".to_string()) && cr.taxa.contains(&"Unassigned".to_string()), "composition: unclassified at the rank named after the rank above");
    let sums_ok = (0..4).all(|si| near(cr.fractions.iter().map(|r| r[si]).sum::<f64>(), 1.0, 1e-12));
    t.check(sums_ok, "composition: each sample adds up to 1");
    let a_lacto = cr.taxa.iter().position(|x| x == "Lactobacillus").map(|i| cr.fractions[i][0]);
    t.check(a_lacto.map_or(false, |v| near(v, 0.5, 1e-12)), "composition: Lactobacillus is 50/100 of sample A");
    let cr2 = taxa_composition(&tab, &taxo, None, &CompositionOptions { rank: "Phylum".into(), top: 1 }).unwrap();
    t.check(cr2.taxa.len() == 2 && cr2.taxa[1] == "Other", "top 1: the first taxon and Other");
    t.fails("unknown rank", taxa_composition(&tab, &taxo, None, &CompositionOptions { rank: "Order".into(), top: 0 }), &["Order", "Kingdom, Phylum, Genus"]);
    let frac = t.write("frac.tsv", "ASV\tA\tB\nx\t0.5\t0.2\ny\t0.5\t0.8\n");
    t.fails("alpha of proportions", alpha_diversity(&frac, None, &AlphaOptions::default()), &["fractional", "counts"]);
    t.fails("PERMANOVA with one group", beta_diversity(&tab, None, &BetaOptions { metadata: Some(meta), group_column: Some("nothere".into()), ..Default::default() }), &["nothere", "group", "week"]);
}

// ---------------------------------------------------------------------
//  7. Tumor mutational burden
// ---------------------------------------------------------------------

fn translate(c: &[u8]) -> u8 {
    let code = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    let i = |b: u8| b"TCAG".iter().position(|&x| x == b).unwrap();
    code[i(c[0]) * 16 + i(c[1]) * 4 + i(c[2])]
}

fn trials_tmb(t: &mut T) {
    println!("tumor mutational burden");
    let mut s = 83u64;
    let genome = random_dna(3000, &mut s);
    let fa = t.write("tmb_genome.fa", &format!(">chr1\n{}\n>chr2\n{}\n", String::from_utf8_lossy(&genome), String::from_utf8_lossy(&random_dna(500, &mut s))));
    // gene A (+): CDS 101-190 and 301-410 (1-based); gene B (-): 1001-1150, 1201-1260
    let gff = "##gff-version 3\nchr1\tt\tgene\t101\t410\t.\t+\t.\tID=gA;Name=GENEA\nchr1\tt\tmRNA\t101\t410\t.\t+\t.\tID=tA;Parent=gA\nchr1\tt\tCDS\t101\t190\t.\t+\t0\tID=cA1;Parent=tA;gene=GENEA\nchr1\tt\tCDS\t301\t410\t.\t+\t0\tID=cA2;Parent=tA;gene=GENEA\nchr1\tt\tCDS\t1001\t1150\t.\t-\t0\tID=cB1;Parent=tB;gene=GENEB\nchr1\tt\tCDS\t1201\t1260\t.\t-\t0\tID=cB2;Parent=tB;gene=GENEB\n";
    let gff_path = t.write("tmb.gff3", gff);
    // the oracle: genomic positions (0-based) of each coding sequence, in order
    let tx_a: Vec<usize> = (100..190).chain(300..410).collect();
    let tx_b: Vec<usize> = (1200..1260).rev().chain((1000..1150).rev()).collect();
    let comp = |b: u8| match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'A',
    };
    let effect = |g: usize, alt: u8| -> Option<&'static str> {
        for (tx, minus) in [(&tx_a, false), (&tx_b, true)] {
            if let Some(c) = tx.iter().position(|&x| x == g) {
                let cs = c - c % 3;
                if cs + 3 > tx.len() {
                    return None;
                }
                let base = |p: usize, alt_here: bool| {
                    let b = if alt_here { alt } else { genome[tx[p]] };
                    if minus {
                        comp(b)
                    } else {
                        b
                    }
                };
                let r: Vec<u8> = (cs..cs + 3).map(|p| base(p, false)).collect();
                let a: Vec<u8> = (cs..cs + 3).map(|p| base(p, p == c)).collect();
                let (x, y) = (translate(&r), translate(&a));
                return Some(if x == y {
                    "synonymous"
                } else if y == b'*' {
                    "stop_gained"
                } else if x == b'*' {
                    "stop_lost"
                } else if cs == 0 && x == b'M' {
                    "start_lost"
                } else {
                    "missense"
                });
            }
        }
        None
    };
    // pick coding SNVs of each kind, in both genes, one at an exon junction
    let mut variants: Vec<(usize, u8, &'static str)> = Vec::new(); // 0-based pos, alt, expected effect
    for &g in &[100usize, 150, 189, 300, 301, 350, 1259, 1200, 1100, 1001] {
        for &alt in b"ACGT" {
            if alt != genome[g] {
                if let Some(e) = effect(g, alt) {
                    if !variants.iter().any(|v| v.0 == g) {
                        variants.push((g, alt, e));
                    }
                }
            }
        }
    }
    // make sure there is a synonymous one: search third codon positions
    if !variants.iter().any(|v| v.2 == "synonymous") {
        'outer: for c in (2..tx_a.len()).step_by(3) {
            let g = tx_a[c];
            for &alt in b"ACGT" {
                if alt != genome[g] && effect(g, alt) == Some("synonymous") {
                    variants.push((g, alt, "synonymous"));
                    break 'outer;
                }
            }
        }
    }
    let n_syn = variants.iter().filter(|v| v.2 == "synonymous").count();
    let n_coding_snv = variants.len();
    let mut vcf = String::from("##fileformat=VCFv4.2\n##tumor_sample=TUMOR\n##normal_sample=NORMAL\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n");
    let mut rows: Vec<(usize, String)> = Vec::new();
    for &(g, alt, _) in &variants {
        rows.push((g, format!("1\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT:AD:DP:AF\t0/0:30,0:30:0\t0/1:20,15:35:0.43", g + 1, genome[g] as char, alt as char)));
    }
    let snv = |g: usize, filter: &str, fmt: &str| -> (usize, String) {
        let alt = if genome[g] == b'A' { 'C' } else { 'A' };
        (g, format!("1\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:AD:DP:AF\t0/0:30,0:30:0\t{}", g + 1, genome[g] as char, alt, filter, fmt))
    };
    rows.push(snv(250, "PASS", "0/1:20,15:35:0.43")); // intron of A: in region, not coding
    rows.push(snv(700, "PASS", "0/1:20,15:35:0.43")); // intergenic
    rows.push(snv(2500, "PASS", "0/1:20,15:35:0.43")); // outside the BED
    rows.push(snv(800, "weak_evidence", "0/1:20,15:35:0.43")); // not PASS
    rows.push(snv(810, "PASS", "0/1:3,2:5:0.40")); // low depth
    rows.push(snv(820, "PASS", "0/1:98,2:100:0.02")); // low VAF
    rows.push(snv(830, "PASS", "0/1")); // no depth, no VAF: kept
    // indels in CDS: frameshift (A) and in-frame deletion (B); one outside CDS
    rows.push((160, format!("1\t161\t.\t{}\t{}T\t.\tPASS\t.\tGT:AD:DP\t0/0:30,0:30\t0/1:20,10:30", genome[160] as char, genome[160] as char)));
    rows.push((1050, format!("1\t1051\t.\t{}\t{}\t.\tPASS\t.\tGT:AD:DP\t0/0:30,0:30\t0/1:20,10:30", String::from_utf8_lossy(&genome[1050..1054]), genome[1050] as char)));
    rows.push((900, format!("1\t901\t.\t{}\t{}\t.\tPASS\t.\tGT:AD:DP\t0/0:30,0:30\t0/1:20,10:30", String::from_utf8_lossy(&genome[900..902]), genome[900] as char)));
    // multi-allelic: REF = 3 bases, ALT1 = deletion, ALT2 = SNV on the first base (trimmed to an SNV)
    let ma_ref = String::from_utf8_lossy(&genome[600..603]).to_string();
    let first_alt = if genome[600] == b'G' { 'T' } else { 'G' };
    rows.push((600, format!("1\t601\t.\t{ma_ref}\t{},{first_alt}{}\t.\tPASS\t.\tGT:AD:DP:AF\t0/0:30,0,0:30:0,0\t0/1:20,8,9:37:0.22,0.24", &ma_ref[..1], &ma_ref[1..])));
    rows.push((610, format!("1\t611\t.\t{}\t<DEL>\t.\tPASS\t.\tGT:AD:DP\t0/0:30,0:30\t0/1:20,10:30", genome[610] as char)));
    rows.sort_by_key(|r| r.0);
    for (_, l) in &rows {
        vcf.push_str(l);
        vcf.push('\n');
    }
    let vcf_path = t.write("tumor.vcf", &vcf);
    let bed = t.write("region.bed", "chr1\t0\t2000\tcaptured\nchr1\t1500\t1800\toverlap\nchr2\t0\t100\n");
    // default: every PASS variant in the region with depth and VAF
    let out = t.p("tmb_variants.tsv");
    let r = tumor_mutational_burden(&vcf_path, &bed, &TmbOptions { output: Some(out.clone()), ..Default::default() }).unwrap();
    // counted: coding SNVs + intron + intergenic + no-info + 3 indels + multi-allelic (2) = n + 2 + 1 + 3 + 2
    let expect = n_coding_snv + 2 + 1 + 3 + 2;
    t.check(r.sample.as_deref() == Some("TUMOR") && r.bed_bp == 2100, format!("tumor sample from ##tumor_sample, BED merged to 2,100 bases ({})", r.bed_bp));
    t.check(r.counted as usize == expect, format!("default TMB: {} counted, {expect} expected", r.counted));
    t.check(near(r.tmb, expect as f64 / 0.0021, 1e-6), format!("TMB = counted / 0.0021 Mb = {:.1}", r.tmb));
    t.check(r.excluded_filter == 1 && r.excluded_depth == 1 && r.excluded_vaf == 1 && r.outside_region == 1 && r.symbolic == 1, "excluded: 1 not PASS, 1 low depth, 1 low VAF, 1 outside, 1 symbolic");
    t.check(r.without_depth == 1 && r.without_vaf == 1, "the variant without depth or VAF is kept and counted apart");
    t.check(r.counted_indel == 4 && r.counted_snv as usize == n_coding_snv + 3 + 1, format!("indels 4 (3 + the multi-allelic deletion), SNVs {} (the multi-allelic SNV trimmed)", r.counted_snv));
    let lines = read_text(&out).lines().count();
    t.check(lines == 1 + r.alleles as usize, "variant table: one line per ALT allele");
    // coding only
    let rc_ = tumor_mutational_burden(&vcf_path, &bed, &TmbOptions { annotation: Some(gff_path.clone()), coding_only: true, ..Default::default() }).unwrap();
    let cds_bp = 90 + 110 + 150 + 60;
    t.check(rc_.coding_bp == Some(cds_bp) && rc_.region_bp == cds_bp, format!("coding region = CDS inside the BED = {cds_bp} bases ({:?})", rc_.coding_bp));
    t.check(rc_.counted as usize == n_coding_snv + 2, format!("coding only: the {n_coding_snv} coding SNVs + 2 coding indels ({})", rc_.counted));
    // the same annotation as RefSeq writes it: accessions, and a region line naming the chromosome
    let refseq_gff = format!(
        "##gff-version 3\nNC_000001.11\tRefSeq\tregion\t1\t3000\t.\t+\t.\tID=NC_000001.11:1..3000;Name=1;chromosome=1;gbkey=Src;genome=chromosome\n{}",
        gff.lines().skip(1).map(|l| l.replacen("chr1", "NC_000001.11", 1) + "\n").collect::<String>()
    );
    let rs_path = t.write("tmb_refseq.gff3", &refseq_gff);
    let rr = tumor_mutational_burden(&vcf_path, &bed, &TmbOptions { annotation: Some(rs_path), coding_only: true, ..Default::default() }).unwrap();
    t.check(rr.counted == rc_.counted && rr.coding_bp == rc_.coding_bp, "a RefSeq GFF (NC_000001.11 + its region line) gives the same coding TMB");
    let pass_hc = vcf.replacen("\tPASS\t", "\tPASS;HighConf\t", 1);
    let ph = t.write("pass_hc.vcf", &pass_hc);
    let rp = tumor_mutational_burden(&ph, &bed, &TmbOptions::default()).unwrap();
    t.check(rp.counted == r.counted, "FILTER \"PASS;HighConf\" counts as PASS (as bcftools -f PASS)");
    // non-synonymous
    let out2 = t.p("tmb_ns.tsv");
    let rn = tumor_mutational_burden(&vcf_path, &bed, &TmbOptions { annotation: Some(gff_path.clone()), reference: Some(fa.clone()), nonsynonymous_only: true, output: Some(out2.clone()), ..Default::default() }).unwrap();
    t.check(rn.counted as usize == n_coding_snv - n_syn + 2 && rn.synonymous as usize == n_syn, format!("non-synonymous: {} counted, {} synonymous left out (expected {} and {n_syn})", rn.counted, rn.synonymous, n_coding_snv - n_syn + 2));
    let tsv = read_text(&out2);
    let mut cons_ok = 0;
    for &(g, alt, e) in &variants {
        let line = tsv.lines().find(|l| l.split('\t').nth(1) == Some(&(g + 1).to_string()) && l.split('\t').nth(3) == Some(&(alt as char).to_string()));
        if line.map_or(false, |l| l.split('\t').nth(10) == Some(e)) {
            cons_ok += 1;
        }
    }
    t.check(cons_ok == variants.len(), format!("consequence of each coding SNV = the oracle's translation ({cons_ok} of {}; both strands, exon junctions)", variants.len()));
    t.check(tsv.contains("\tframeshift\tGENEA\t") && tsv.contains("\tinframe_indel\tGENEB\t"), "indels: frameshift in GENEA, in-frame deletion in GENEB");
    t.check(rn.ref_mismatch == 0, "no REF mismatch against the genome");
    // Strelka-style counts and chr names without 'chr'
    let strelka = format!(
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n1\t{}\t.\t{}\t{}\t.\tPASS\tSOMATIC\tDP:AU:CU:GU:TU\t30:0,0:0,0:0,0:0,0\t40:{}\n",
        801,
        genome[800] as char,
        if genome[800] == b'A' { 'C' } else { 'A' },
        {
            let r = genome[800];
            let a = if r == b'A' { b'C' } else { b'A' };
            let tier = |b: u8| if b == r { "38,38" } else if b == a { "1,1" } else { "0,0" };
            format!("{}:{}:{}:{}", tier(b'A'), tier(b'C'), tier(b'G'), tier(b'T'))
        }
    );
    let sp = t.write("strelka.vcf", &strelka);
    let rs = tumor_mutational_burden(&sp, &bed, &TmbOptions::default()).unwrap();
    t.check(rs.excluded_vaf == 1 && rs.counted == 0, "Strelka tier counts: VAF 1/39 below 5%, not counted");
    // errors
    let two = vcf.replace("##tumor_sample=TUMOR\n", "").replace("NORMAL\tTUMOR", "N1\tT1");
    let two_path = t.write("two.vcf", &two);
    t.fails("two samples, none marked", tumor_mutational_burden(&two_path, &bed, &TmbOptions::default()), &["N1, T1", "tumor_sample"]);
    t.fails("nonsynonymous without reference", tumor_mutational_burden(&vcf_path, &bed, &TmbOptions { annotation: Some(gff_path.clone()), nonsynonymous_only: true, ..Default::default() }), &["reference"]);
    let bed_other = t.write("other.bed", "chrX\t0\t1000\n");
    t.fails("no sequence in common", tumor_mutational_burden(&vcf_path, &bed_other, &TmbOptions::default()), &["same genome", "chrX"]);
    let bad_bed = t.write("bad.bed", "chr1\t100\n");
    t.fails("bad BED", tumor_mutational_burden(&vcf_path, &bad_bed, &TmbOptions::default()), &["line 1", "BED"]);
    t.check(r.resumo().contains("TMB =") && r.json().contains("\"tmb\":"), "json and resumo");
}

// ---------------------------------------------------------------------
//  Real data
// ---------------------------------------------------------------------

mod real {
    use super::*;

    fn secs(t: Instant) -> String {
        format!("{:.1} s", t.elapsed().as_secs_f64())
    }

    /// Semi-global edit distance: `q` against the best window of `r`.
    fn semiglobal(q: &[u8], r: &[u8]) -> usize {
        let mut prev = vec![0usize; r.len() + 1];
        for i in 1..=q.len() {
            let mut cur = vec![i; r.len() + 1];
            for j in 1..=r.len() {
                cur[j] = (prev[j - 1] + usize::from(q[i - 1] != r[j - 1])).min(prev[j] + 1).min(cur[j - 1] + 1);
            }
            prev = cur;
        }
        *prev.iter().min().unwrap()
    }

    pub fn sop(data: &str, work: &str, vsearch: Option<&str>) {
        std::fs::create_dir_all(work).unwrap();
        let sop = format!("{data}/16s/mothur_miseq_sop/MiSeq_SOP");
        let w = |n: &str| format!("{work}/{n}");
        let samples = read_sample_sheet(&format!("{sop}/stability.files")).unwrap();
        println!("== MiSeq SOP: {} samples", samples.len());
        let t0 = Instant::now();
        let r = amplicon_asvs(&samples, &w("asvs.fasta"), &w("table.tsv"), Some(&w("track.tsv")), &AsvOptions::default()).unwrap();
        println!("amplicon_asvs: {}", secs(t0));
        print!("{}", r.resumo());
        let lens: HashMap<usize, usize> = r.asvs.iter().fold(HashMap::new(), |mut m, a| {
            *m.entry(a.sequence.len()).or_insert(0) += 1;
            m
        });
        let mut lv: Vec<_> = lens.into_iter().collect();
        lv.sort();
        println!("ASV lengths: {lv:?}");
        let mf: [u64; 4] = r.samples.iter().fold([0; 4], |a, x| [a[0] + x.merge_too_short, a[1] + x.merge_no_alignment, a[2] + x.merge_too_many_diffs, a[3] + x.merge_staggered]);
        let ff: [u64; 3] = r.samples.iter().fold([0; 3], |a, x| [a[0] + x.filter_length, a[1] + x.filter_n, a[2] + x.filter_expected_errors]);
        println!("merge failures: too_short {} no_alignment {} too_many_diffs {} staggered {}; filters: length {} N {} expected_errors {}", mf[0], mf[1], mf[2], mf[3], ff[0], ff[1], ff[2]);
        // the Mock community
        let refs = amplicon::read_fasta(&format!("{sop}/HMP_MOCK.v35.fasta")).unwrap();
        let table = read_feature_table(&w("table.tsv")).unwrap();
        let mi = table.samples.iter().position(|x| x == "Mock").unwrap();
        let (mut exact, mut spurious, mut reads_exact, mut reads_sp, mut err_bases, mut bases) = (0, 0, 0u64, 0u64, 0u64, 0u64);
        let mut strains: HashSet<String> = HashSet::new();
        println!("Mock ASVs:");
        for (fi, f) in table.features.iter().enumerate() {
            let n = table.counts[fi][mi] as u64;
            if n == 0 {
                continue;
            }
            let seq = r.asvs.iter().find(|a| &a.name == f).unwrap().sequence.as_bytes().to_vec();
            let hit: Vec<&String> = refs.iter().filter(|(_, x)| x.windows(seq.len()).any(|win| win == seq.as_slice())).map(|(id, _)| id).collect();
            let d = if hit.is_empty() { refs.iter().map(|(_, x)| semiglobal(&seq, x)).min().unwrap() } else { 0 };
            if d == 0 {
                exact += 1;
                reads_exact += n;
                for h in &hit {
                    strains.insert(h.split('.').take(2).collect::<Vec<_>>().join("."));
                }
            } else {
                spurious += 1;
                reads_sp += n;
            }
            err_bases += n * d as u64;
            bases += n * seq.len() as u64;
            println!("  {f:<8} {n:>6} reads  {}", if d == 0 { format!("exact: {}", hit.iter().map(|x| x.as_str()).collect::<Vec<_>>().join(" ")) } else { format!("{d} differences from the closest reference") });
        }
        println!(
            "Mock: {exact} exact ASVs ({} of the 20 strains), {spurious} not exact; reads in exact ASVs {:.2}%; residual error rate {:.2e} per base",
            strains.len(),
            100.0 * reads_exact as f64 / (reads_exact + reads_sp) as f64,
            err_bases as f64 / bases as f64
        );
        // mouse-only table (without the Mock) for diversity
        let mouse: Vec<usize> = (0..table.samples.len()).filter(|&i| i != mi).collect();
        let mut mt = String::from("ASV");
        for &i in &mouse {
            mt.push('\t');
            mt.push_str(&table.samples[i]);
        }
        mt.push('\n');
        for (fi, f) in table.features.iter().enumerate() {
            let row: Vec<String> = mouse.iter().map(|&i| (table.counts[fi][i] as u64).to_string()).collect();
            if row.iter().all(|x| x == "0") {
                continue;
            }
            mt.push_str(f);
            mt.push('\t');
            mt.push_str(&row.join("\t"));
            mt.push('\n');
        }
        std::fs::write(w("table_mouse.tsv"), &mt).unwrap();
        let depths: Vec<u64> = mouse.iter().map(|&i| table.counts.iter().map(|c| c[i] as u64).sum()).collect();
        let min_depth = *depths.iter().min().unwrap();
        println!("\n== diversity of the 19 mouse samples (depths {} to {})", min_depth, depths.iter().max().unwrap());
        let ar = alpha_diversity(&w("table_mouse.tsv"), Some(&w("alpha.tsv")), &AlphaOptions { curves_output: Some(w("curves.tsv")), ..Default::default() }).unwrap();
        print!("{}", ar.resumo());
        let ar2 = alpha_diversity(&w("table_mouse.tsv"), Some(&w("alpha_rarefied.tsv")), &AlphaOptions { rarefy_depth: Some(min_depth), curve_steps: 0, ..Default::default() }).unwrap();
        let mean = |rows: &[AlphaRow], early: bool| {
            let v: Vec<f64> = rows.iter().filter(|r| { let d: u32 = r.sample[3..].parse().unwrap(); (d < 100) == early }).map(|r| r.shannon).collect();
            v.iter().sum::<f64>() / v.len() as f64
        };
        println!("Shannon (rarefied to {min_depth}): early {:.3}, late {:.3}", mean(&ar2.rows, true), mean(&ar2.rows, false));
        let design = format!("{sop}/mouse.time.design");
        for metric in ["braycurtis", "jaccard"] {
            let br = beta_diversity(
                &w("table_mouse.tsv"),
                Some(&w(&format!("{metric}.tsv"))),
                &BetaOptions { metric: metric.into(), metadata: Some(design.clone()), group_column: Some("time".into()), pcoa_output: Some(w(&format!("pcoa_{metric}.tsv"))), ..Default::default() },
            )
            .unwrap();
            print!("{}", br.resumo());
            let pc1: Vec<(String, f64)> = br.samples.iter().zip(&br.pcoa.coordinates).map(|(s, c)| (s.clone(), c[0])).collect();
            let early: Vec<f64> = pc1.iter().filter(|(s, _)| s[3..].parse::<u32>().unwrap() < 100).map(|x| x.1).collect();
            let late: Vec<f64> = pc1.iter().filter(|(s, _)| s[3..].parse::<u32>().unwrap() >= 100).map(|x| x.1).collect();
            println!(
                "  PC1 early {:.3}..{:.3}, late {:.3}..{:.3}",
                early.iter().cloned().fold(f64::INFINITY, f64::min),
                early.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
                late.iter().cloned().fold(f64::INFINITY, f64::min),
                late.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            );
        }
        let bk = beta_diversity(&w("table_mouse.tsv"), None, &BetaOptions { metric: "kmer".into(), asv_fasta: Some(w("asvs.fasta")), metadata: Some(design.clone()), group_column: Some("time".into()), ..Default::default() }).unwrap();
        print!("{}", bk.resumo());
        // taxonomy
        println!("\n== taxonomy");
        let rdp = format!("{data}/16s/reference/rdp19_dada2/rdp_19_toGenus_trainset.fa.gz");
        let v4 = Some(("GTGYCAGCMGCCGCGGTAA".to_string(), "GGACTACNVGGGTWTCTAAT".to_string()));
        for (name, refp, region) in [
            ("rdp19", rdp.clone(), None),
            ("rdp19_v4", rdp.clone(), v4.clone()),
            ("gtdb_r220", format!("{data}/16s/reference/gtdb_r220_dada2/GTDB_bac120_arc53_ssu_r220_genus.fa.gz"), None),
            ("gtdb_r220_v4", format!("{data}/16s/reference/gtdb_r220_dada2/GTDB_bac120_arc53_ssu_r220_genus.fa.gz"), v4.clone()),
        ] {
            let mp = w(&format!("{name}.mmgtax"));
            let (pf, pr) = match region {
                Some((f, r)) => (Some(f), Some(r)),
                None => (None, None),
            };
            println!("-- {name}");
            let ts = train_classifier(&refp, Some(&mp), &TrainOptions { primer_forward: pf, primer_reverse: pr, ..Default::default() }).unwrap();
            print!("{}", ts.resumo());
            let t1 = Instant::now();
            let _ = TaxModel::load(&mp).unwrap();
            println!("  model reload: {}", secs(t1));
            let tx = classify_asvs(&w("asvs.fasta"), &mp, &w(&format!("taxonomy_{name}.tsv")), &ClassifyOptions::default()).unwrap();
            print!("{}", tx.resumo());
            let mock_genera: Vec<String> = table
                .features
                .iter()
                .enumerate()
                .filter(|(fi, _)| table.counts[*fi][mi] > 0.0)
                .map(|(_, f)| {
                    let a = tx.assignments.iter().find(|a| &a.id == f).unwrap();
                    format!("{f}={}", a.lineage.get(5).cloned().unwrap_or_else(|| format!("({})", a.lineage.last().cloned().unwrap_or_default())))
                })
                .collect();
            println!("  Mock ASV genera: {}", mock_genera.join(", "));
            let comp = taxa_composition(&w("table_mouse.tsv"), &w(&format!("taxonomy_{name}.tsv")), Some(&w(&format!("family_{name}.tsv"))), &CompositionOptions { rank: "Family".into(), top: 10 }).unwrap();
            print!("{}", comp.resumo());
        }
        let sp = classify_asvs(
            &w("asvs.fasta"),
            &w("rdp19.mmgtax"),
            &w("taxonomy_rdp19_species.tsv"),
            &ClassifyOptions { species_reference: Some(format!("{data}/16s/reference/rdp19_dada2/rdp_19_toSpecies_trainset.fa.gz")), ..Default::default() },
        )
        .unwrap();
        println!("exact species (RDP 19 toSpecies): {} of {} ASVs", sp.species_assigned, sp.assignments.len());
        if let Some(v) = vsearch {
            compare_vsearch(v, &samples, &r, work);
        }
    }

    fn run(v: &str, args: &[&str]) -> String {
        let o = std::process::Command::new(v).args(args).output().expect("vsearch");
        let e = String::from_utf8_lossy(&o.stderr).to_string();
        if !o.status.success() {
            panic!("vsearch {args:?}: {e}");
        }
        e
    }

    fn fastx(path: &str) -> Vec<(String, Vec<u8>)> {
        let mut out = Vec::new();
        let mut r = SeqIn::open(path).unwrap();
        let mut rec = SeqRec::default();
        while r.next(&mut rec).unwrap() {
            out.push((String::from_utf8_lossy(&rec.head).to_string(), rec.seq.clone()));
        }
        out
    }

    fn compare_vsearch(v: &str, samples: &[AmpliconSample], r: &AsvResult, work: &str) {
        println!("\n== vsearch ({})", run(v, &["--version"]).lines().next().unwrap_or(""));
        let w = |n: &str| format!("{work}/{n}");
        let vd = w("vsearch");
        std::fs::create_dir_all(&vd).unwrap();
        let mut all_filtered = String::new();
        let (mut same_seq, mut both_merged, mut only_us, mut only_them) = (0u64, 0u64, 0u64, 0u64);
        let opt = MergeOptions::default();
        println!("  {:<8} {:>8} {:>8} {:>8} {:>8}", "sample", "merged", "vs.merg", "filtered", "vs.filt");
        for (si, s) in samples.iter().enumerate() {
            let mq = format!("{vd}/{}.merged.fq", s.name);
            let fq = format!("{vd}/{}.filtered.fa", s.name);
            run(v, &["--fastq_mergepairs", &s.r1, "--reverse", s.r2.as_ref().unwrap(), "--fastqout", &mq, "--threads", "4", "--quiet"]);
            run(v, &["--fastq_filter", &mq, "--fastq_maxee", "1", "--fastq_maxns", "0", "--fastaout", &fq, "--relabel", &format!("{}.", s.name), "--quiet"]);
            let merged = fastx(&mq);
            let filtered = fastx(&fq);
            for (_, sq) in &filtered {
                all_filtered.push_str(&format!(">{};sample={}\n{}\n", all_filtered.len(), s.name, String::from_utf8_lossy(sq)));
            }
            println!("  {:<8} {:>8} {:>8} {:>8} {:>8}", s.name, r.samples[si].merged, merged.len(), r.samples[si].filtered, filtered.len());
            // per pair: the same merged sequence?
            let theirs: HashMap<String, Vec<u8>> = merged.into_iter().map(|(h, sq)| (h.split_whitespace().next().unwrap().to_string(), sq)).collect();
            let mut a = SeqIn::open(&s.r1).unwrap();
            let mut b = SeqIn::open(s.r2.as_ref().unwrap()).unwrap();
            let (mut x, mut y) = (SeqRec::default(), SeqRec::default());
            while a.next(&mut x).unwrap() && b.next(&mut y).unwrap() {
                let id = String::from_utf8_lossy(&x.head).split_whitespace().next().unwrap().to_string();
                let q1: Vec<u8> = x.qual.iter().map(|c| c - 33).collect();
                let q2: Vec<u8> = y.qual.iter().map(|c| c - 33).collect();
                let ours = merge_pair(&x.seq, &q1, &y.seq, &q2, &opt).ok().map(|m| m.seq);
                match (ours, theirs.get(&id)) {
                    (Some(o), Some(t)) => {
                        both_merged += 1;
                        if &o == t {
                            same_seq += 1;
                        }
                    }
                    (Some(_), None) => only_us += 1,
                    (None, Some(_)) => only_them += 1,
                    (None, None) => {}
                }
            }
        }
        println!("  pairs merged by both: {both_merged}, identical merged sequence {same_seq} ({:.2}%); merged only by us {only_us}, only by vsearch {only_them}", 100.0 * same_seq as f64 / both_merged as f64);
        let pooled = format!("{vd}/all_filtered.fa");
        std::fs::write(&pooled, all_filtered).unwrap();
        let derep = format!("{vd}/derep.fa");
        run(v, &["--derep_fulllength", &pooled, "--sizeout", "--output", &derep, "--quiet"]);
        let zotus = format!("{vd}/zotus.fa");
        run(v, &["--cluster_unoise", &derep, "--sizein", "--minsize", "8", "--unoise_alpha", "2", "--centroids", &zotus, "--sizeout", "--quiet"]);
        let nonchim = format!("{vd}/zotus_nonchim.fa");
        let chim = format!("{vd}/zotus_chim.fa");
        run(v, &["--uchime3_denovo", &zotus, "--sizein", "--nonchimeras", &nonchim, "--chimeras", &chim, "--quiet"]);
        let z: HashSet<Vec<u8>> = fastx(&zotus).into_iter().map(|x| x.1.to_ascii_uppercase()).collect();
        let zn: HashSet<Vec<u8>> = fastx(&nonchim).into_iter().map(|x| x.1.to_ascii_uppercase()).collect();
        let ours: HashSet<Vec<u8>> = r.asvs.iter().map(|a| a.sequence.as_bytes().to_vec()).collect();
        let ours_all: HashSet<Vec<u8>> = ours.iter().cloned().chain(r.chimeras.iter().map(|c| c.sequence.as_bytes().to_vec())).collect();
        println!(
            "  centroids: ours {}, vsearch --cluster_unoise {}; shared {} (Jaccard {:.3})",
            ours_all.len(),
            z.len(),
            ours_all.intersection(&z).count(),
            ours_all.intersection(&z).count() as f64 / ours_all.union(&z).count() as f64
        );
        println!(
            "  after chimeras: ours {}, vsearch --uchime3_denovo {}; shared {} (Jaccard {:.3}); chimeras ours {}, vsearch {}",
            ours.len(),
            zn.len(),
            ours.intersection(&zn).count(),
            ours.intersection(&zn).count() as f64 / ours.union(&zn).count() as f64,
            r.chimeras.len(),
            z.len() - zn.len()
        );
        // the same denoising from vsearch's own filtered reads, to separate the steps
        let mut seqs: Vec<Vec<u8>> = Vec::new();
        let mut abund: Vec<u64> = Vec::new();
        for (h, sq) in fastx(&derep) {
            let size: u64 = h.split(";size=").nth(1).and_then(|x| x.trim_end_matches(';').parse().ok()).unwrap_or(1);
            seqs.push(sq.to_ascii_uppercase());
            abund.push(size);
        }
        let den = unoise(&seqs, &abund, &DenoiseOptions::default());
        let ourc: HashSet<Vec<u8>> = den.centroids.iter().map(|&c| seqs[c].clone()).collect();
        println!(
            "  UNOISE on vsearch's dereplicated reads: ours {} centroids, vsearch {}; shared {} (Jaccard {:.3})",
            ourc.len(),
            z.len(),
            ourc.intersection(&z).count(),
            ourc.intersection(&z).count() as f64 / ourc.union(&z).count() as f64
        );
        // chimeras on vsearch's centroids
        let zl: Vec<(Vec<u8>, u64)> = fastx(&zotus)
            .into_iter()
            .map(|(h, sq)| (sq.to_ascii_uppercase(), h.split(";size=").nth(1).and_then(|x| x.trim_end_matches(';').parse().ok()).unwrap_or(1)))
            .collect();
        let zs: Vec<&[u8]> = zl.iter().map(|x| x.0.as_slice()).collect();
        let za: Vec<u64> = zl.iter().map(|x| x.1).collect();
        let ch = find_chimeras(&zs, &za, &ChimeraOptions::default());
        let ours_ch: HashSet<Vec<u8>> = ch.iter().enumerate().filter(|(_, c)| c.is_some()).map(|(i, _)| zl[i].0.clone()).collect();
        let theirs_ch: HashSet<Vec<u8>> = z.difference(&zn).cloned().collect();
        println!(
            "  chimeras on vsearch's {} centroids: ours {}, vsearch {}, both {}",
            zl.len(),
            ours_ch.len(),
            theirs_ch.len(),
            ours_ch.intersection(&theirs_ch).count()
        );
        let mut hist: HashMap<String, usize> = HashMap::new();
        for (i, c) in ch.iter().enumerate() {
            if let Some(m) = c {
                let theirs = theirs_ch.contains(&zl[i].0);
                let key = format!(
                    "{} votes>={} closest>={} skew<1/{}",
                    if theirs { "both" } else { "ours only" },
                    m.votes_left.min(m.votes_right).min(6),
                    m.closest_parent_diffs.min(10),
                    if za[m.parent_left].min(za[m.parent_right]) as f64 / za[i] as f64 >= 64.0 { 64 } else { 16 }
                );
                *hist.entry(key).or_insert(0) += 1;
                if std::env::var("CHIM_DEBUG").is_ok() && !theirs {
                    println!(
                        "    ours only: size {} parents {}+{} (sizes {} {}), bp {}, votes {}/{}, closest {} diffs",
                        za[i], m.parent_left, m.parent_right, za[m.parent_left], za[m.parent_right], m.breakpoint, m.votes_left, m.votes_right, m.closest_parent_diffs
                    );
                }
            }
        }
        let mut hv: Vec<_> = hist.into_iter().collect();
        hv.sort();
        for (k, n) in hv {
            println!("    {k}: {n}");
        }
    }

    /// TMB of the SEQC2 HCC1395 somatic truth set (SNVs + indels in the
    /// high-confidence regions), whole and coding-only (RefSeq GFF).
    pub fn tmb(seqc2: &str, gff: &str, work: &str) {
        std::fs::create_dir_all(work).unwrap();
        let snv = format!("{seqc2}/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz");
        let ind = format!("{seqc2}/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz");
        let bed = format!("{seqc2}/High-Confidence_Regions_v1.2.bed");
        // one VCF: the SNV header, the records of both
        let mut text = String::new();
        externo::abre_leitura(&snv).unwrap().read_to_string(&mut text).unwrap();
        let mut ind_text = String::new();
        externo::abre_leitura(&ind).unwrap().read_to_string(&mut ind_text).unwrap();
        for l in ind_text.lines().filter(|l| !l.starts_with('#')) {
            text.push_str(l);
            text.push('\n');
        }
        let both = format!("{work}/seqc2_snv_indel.vcf");
        std::fs::write(&both, text).unwrap();
        let t0 = Instant::now();
        let r = tumor_mutational_burden(&both, &bed, &TmbOptions { output: Some(format!("{work}/seqc2_tmb.tsv")), ..Default::default() }).unwrap();
        println!("== SEQC2 HCC1395, whole HC region ({})", secs(t0));
        print!("{}", r.resumo());
        let t0 = Instant::now();
        let rc = tumor_mutational_burden(&both, &bed, &TmbOptions { annotation: Some(gff.to_string()), coding_only: true, output: Some(format!("{work}/seqc2_tmb_coding.tsv")), ..Default::default() }).unwrap();
        println!("== SEQC2 HCC1395, coding only (CDS of the RefSeq GFF inside the HC region; {})", secs(t0));
        print!("{}", rc.resumo());
        println!("BED {:.1} Mb, CDS inside it {:.2} Mb", rc.bed_bp as f64 / 1e6, rc.coding_bp.unwrap_or(0) as f64 / 1e6);
    }

    /// Hold-out accuracy of the classifier on RDP 19: 1 sequence in 20 left
    /// out of training (from genera with at least 2), classified whole and as
    /// its V4 region (515F-806R).
    pub fn holdout(data: &str, work: &str) {
        std::fs::create_dir_all(work).unwrap();
        let refp = format!("{data}/16s/reference/rdp19_dada2/rdp_19_toGenus_trainset.fa.gz");
        let mut r = SeqIn::open(&refp).unwrap();
        let mut rec = SeqRec::default();
        let mut all: Vec<(Vec<String>, Vec<u8>, String)> = Vec::new();
        while r.next(&mut rec).unwrap() {
            if let Some(l) = parse_lineage(&rec.head) {
                all.push((l, rec.seq.clone(), String::from_utf8_lossy(&rec.head).to_string()));
            }
        }
        let mut per_genus: HashMap<String, usize> = HashMap::new();
        for a in &all {
            *per_genus.entry(a.0.join(";")).or_insert(0) += 1;
        }
        let mut train = String::new();
        let mut test: Vec<(Vec<String>, Vec<u8>)> = Vec::new();
        let mut s = 99u64;
        for a in &all {
            if a.0.len() == 6 && per_genus[&a.0.join(";")] >= 2 && rnd(&mut s, 20) == 0 {
                test.push((a.0.clone(), a.1.clone()));
            } else {
                train.push_str(&format!(">{}\n{}\n", a.2, String::from_utf8_lossy(&a.1)));
            }
        }
        let tp = format!("{work}/holdout_train.fa");
        std::fs::write(&tp, train).unwrap();
        let t0 = Instant::now();
        let m = TaxModel::train(&tp, &TrainOptions::default()).unwrap();
        println!("trained on {} sequences, {} lineages, {}; {} held out", m.n_sequences, m.classes(), secs(t0), test.len());
        let f = Primer::new("GTGYCAGCMGCCGCGGTAA", 0.1).unwrap();
        let rv = Primer::new("GGACTACNVGGGTWTCTAAT", 0.1).unwrap();
        let v4: Vec<Option<Vec<u8>>> = test
            .iter()
            .map(|(_, q)| {
                let up: Vec<u8> = q.iter().map(|c| c.to_ascii_uppercase()).collect();
                let a = find_forward_primer(&up, &f, up.len().saturating_sub(19))?;
                let b = find_reverse_primer(&up[a..], &rv, up.len().saturating_sub(a + 20))?;
                Some(up[a..a + b].to_vec())
            })
            .collect();
        println!("V4 found in {} of {}", v4.iter().filter(|x| x.is_some()).count(), test.len());
        let t1 = Instant::now();
        let mv4 = TaxModel::train(&tp, &TrainOptions { primer_forward: Some("GTGYCAGCMGCCGCGGTAA".into()), primer_reverse: Some("GGACTACNVGGGTWTCTAAT".into()), ..Default::default() }).unwrap();
        println!("V4-region model: {} sequences, {} pairs (whole: {}), {}", mv4.n_sequences, mv4.pairs(), m.pairs(), secs(t1));
        for conf in [0.8, 0.5] {
            for (label, model, region) in [("full length", &m, false), ("V4, whole model", &m, true), ("V4, V4 model", &mv4, true)] {
                let opt = ClassifyOptions { min_confidence: conf, ..Default::default() };
                let res: Vec<(usize, usize, usize)> = test
                    .par_iter()
                    .zip(v4.par_iter())
                    .filter_map(|((lin, q), v)| {
                        let seq = if region { v.as_ref()? } else { q };
                        let a = model.classify(seq, &opt);
                        Some(if a.lineage.len() == 6 {
                            if a.lineage == *lin {
                                (1, 0, 0)
                            } else {
                                (0, 1, 0)
                            }
                        } else if a.lineage.iter().zip(lin).all(|(x, y)| x == y) {
                            (0, 0, 1)
                        } else {
                            (0, 1, 0)
                        })
                    })
                    .collect();
                let n = res.len() as f64;
                let (c, wr, u) = res.iter().fold((0, 0, 0), |a, x| (a.0 + x.0, a.1 + x.1, a.2 + x.2));
                println!("  {label:<12} confidence {conf}: genus right {:.1}%, wrong {:.1}%, stopped above the genus (right so far) {:.1}%  (n = {})", 100.0 * c as f64 / n, 100.0 * wr as f64 / n, 100.0 * u as f64 / n, res.len());
            }
        }
    }

    use rayon::prelude::*;
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("real") => {
            real::sop(args.get(2).expect("research_data folder"), args.get(3).expect("work folder"), args.get(4).map(String::as_str));
            return;
        }
        Some("holdout") => {
            real::holdout(args.get(2).expect("research_data folder"), args.get(3).expect("work folder"));
            return;
        }
        Some("train") => {
            // train <reference> <model out> <max_depth|-> [asvs.fasta taxonomy.tsv]
            let depth = args.get(4).and_then(|x| x.parse().ok());
            let ts = train_classifier(&args[2], Some(&args[3]), &TrainOptions { max_depth: depth, ..Default::default() }).unwrap();
            print!("{}", ts.resumo());
            if let (Some(a), Some(o)) = (args.get(5), args.get(6)) {
                let r = classify_asvs(a, &args[3], o, &ClassifyOptions::default()).unwrap();
                print!("{}", r.resumo());
            }
            return;
        }
        Some("tmb") => {
            real::tmb(args.get(2).expect("SEQC2 folder"), args.get(3).expect("RefSeq GFF"), args.get(4).expect("work folder"));
            return;
        }
        _ => {}
    }
    let dir = std::env::temp_dir().join(format!("mmgen_amplicon_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut t = T { ok: 0, total: 0, failed: Vec::new(), dir: dir.clone() };
    let t0 = Instant::now();
    let sections: [fn(&mut T); 7] = [trials_merge, trials_primers, trials_unoise, trials_pipeline, trials_classifier, trials_diversity, trials_tmb];
    for f in sections {
        let ts = Instant::now();
        let before = t.total;
        f(&mut t);
        println!("    {} checks, {:.1} s", t.total - before, ts.elapsed().as_secs_f64());
    }
    let _ = std::fs::remove_dir_all(&dir);
    println!("{} of {} checks passed ({:.1} s)", t.ok, t.total, t0.elapsed().as_secs_f64());
    if t.ok != t.total {
        std::process::exit(1);
    }
}
