#![allow(dead_code)]
// Trials of the k-mer counter for inputs of any size (src/kmer_count.rs).
//
// Every table is compared, line by line, with a reference written here in the
// most naive way: a `HashMap<String, u32>` filled by slicing the text of each
// record, with the reverse complement done on the text. Nothing in the
// reference uses the code under test. Each case is counted twice -- with
// everything in memory and with the 1 MB budget, which cuts the k-mers in
// runs on disk -- and the two outputs must be the same bytes.
//
// Usage: cargo run --release --bin kmer_count
//        cargo run --release --bin kmer_count -- bench [Mbp] [budget in MB]
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
#[path = "../../../src/kmer.rs"]
mod kmer;
#[path = "../../../src/kmer_count.rs"]
mod kmer_count;

use kmer_count::*;
use std::collections::{BTreeMap, HashMap};
use std::io::{BufRead, Read, Write};
use std::path::{Path, PathBuf};

/// The fixed-seed generator of the other trials, but handing out the HIGH
/// bits of the state. The low bits of a generator of this kind repeat after
/// very few draws -- `(state >> 11) % 4` has a period of 8192 -- and a "random"
/// genome made of them would be the same 8 kb over and over.
fn lcg(s: &mut u64) -> u64 {
    *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    *s >> 33
}

fn random_dna(n: usize, s: &mut u64) -> String {
    (0..n).map(|_| char::from(b"ACGT"[(lcg(s) % 4) as usize])).collect()
}

// ---------------------------------------------------------------------
//  The reference, on text
// ---------------------------------------------------------------------

fn rc_text(t: &str) -> String {
    assert!(t.is_ascii());
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
            other => char::from(other),
        })
        .collect()
}

fn encode_text(t: &str) -> u128 {
    t.chars().fold(0u128, |v, c| {
        v * 4
            + match c {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                'T' => 3,
                _ => panic!("not a base: {c}"),
            }
    })
}

/// Counts by slicing. `records` are the sequences as the file has them (any
/// case, any letter), without white space.
fn reference(records: &[String], k: usize, canonical: bool) -> HashMap<String, u32> {
    let mut m: HashMap<String, u32> = HashMap::with_capacity(records.iter().map(String::len).sum::<usize>().min(1 << 20));
    for r in records {
        let up = r.to_ascii_uppercase();
        if up.len() < k {
            continue;
        }
        for i in 0..=up.len() - k {
            let w = &up[i..i + k];
            if !w.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                continue;
            }
            let key = if canonical {
                let rc = rc_text(w);
                if rc.as_str() < w {
                    rc
                } else {
                    w.to_string()
                }
            } else {
                w.to_string()
            };
            *m.entry(key).or_insert(0) += 1;
        }
    }
    m
}

struct Expected {
    lines: Vec<String>,
    counts: KmerCounts,
}

fn expected(records: &[String], k: u8, canonical: bool, min_count: u32) -> Expected {
    expected_of(&reference(records, usize::from(k), canonical), records, k, canonical, min_count)
}

fn expected_of(m: &HashMap<String, u32>, records: &[String], k: u8, canonical: bool, min_count: u32) -> Expected {
    let mut keys: Vec<&String> = m.keys().collect();
    keys.sort();
    let lines: Vec<String> = keys.iter().filter(|key| m[**key] >= min_count.max(1)).map(|key| format!("{key}\t{}", m[*key])).collect();
    let mut histogram: BTreeMap<u32, u64> = BTreeMap::new();
    for &c in m.values() {
        *histogram.entry(c.min(10_001)).or_insert(0) += 1;
    }
    Expected {
        counts: KmerCounts {
            k,
            canonical,
            sequences: records.len() as u64,
            bases: records.iter().map(|r| r.len() as u64).sum(),
            total_kmers: m.values().map(|&c| u64::from(c)).sum(),
            distinct_kmers: m.len() as u64,
            unique_kmers: m.values().filter(|&&c| c == 1).count() as u64,
            max_count: m.values().copied().max().unwrap_or(0),
            written: lines.len() as u64,
            runs: 0,
            histogram: histogram.into_iter().collect(),
        },
        lines,
    }
}

// ---------------------------------------------------------------------
//  Files
// ---------------------------------------------------------------------

/// 70 bases that would become k-mers of any length if a header were read as
/// sequence.
const BAIT: &str = "ACGTTGCAAGCTTCGAGGATCCATGCATGCCGTAACGTTAGCTAGGCTAACGTACGATCGATTGCACGTAG";

/// Each line takes the next width of the list, in turn.
fn fasta_text(records: &[String], widths: &[usize], eol: &str) -> String {
    let mut o = String::new();
    let mut turn = 0;
    for (i, r) in records.iter().enumerate() {
        o.push_str(&format!(">seq{i} {BAIT} len={}{eol}", r.len()));
        let mut rest = r.as_str();
        while !rest.is_empty() {
            let n = widths[turn % widths.len()].min(rest.len());
            turn += 1;
            o.push_str(&rest[..n]);
            o.push_str(eol);
            rest = &rest[n..];
        }
    }
    o
}

/// Qualities made of the characters most likely to fool a reader: bases, and
/// the marks that start a header, a FASTA record, a comment or a `+` line.
fn fastq_text(records: &[String], eol: &str, s: &mut u64) -> String {
    let marks = b"@+>;ACGTacgtNI!~5";
    let mut o = String::new();
    for (i, r) in records.iter().enumerate() {
        let mut qual: Vec<u8> = (0..r.len()).map(|_| marks[(lcg(s) % marks.len() as u64) as usize]).collect();
        if let Some(first) = qual.first_mut() {
            *first = [b'@', b'+', b'A'][i % 3];
        }
        let bait = if i % 8 == 0 { BAIT } else { "" };
        let plus = if i % 5 == 0 { format!("+read{i} {BAIT}") } else { "+".to_string() };
        o.push_str(&format!("@read{i} {bait}{eol}{r}{eol}{plus}{eol}{}{eol}", String::from_utf8(qual).unwrap()));
    }
    o
}

/// Lower case in stretches, runs of N (shorter than, equal to and longer than
/// the k-mers), and single IUPAC letters, gaps and a U.
fn dirty(clean: &str, s: &mut u64) -> String {
    let mut o = String::with_capacity(clean.len());
    let mut lower = false;
    let mut run_of_n = 0usize;
    for c in clean.chars() {
        let r = lcg(s);
        if r % 97 == 0 {
            lower = !lower;
        }
        if run_of_n > 0 {
            run_of_n -= 1;
            o.push(if lower { 'n' } else { 'N' });
        } else if r % 1201 == 7 {
            run_of_n = [0usize, 1, 29, 30, 31, 499][(lcg(s) % 6) as usize];
            o.push('N');
        } else if r % 311 == 3 {
            o.push(char::from(b"RYKMSWBDHVNryn-*.U"[(lcg(s) % 18) as usize]));
        } else {
            o.push(if lower { c.to_ascii_lowercase() } else { c });
        }
    }
    o
}

struct Trial {
    dir: PathBuf,
    runs_dir: PathBuf,
    passed: usize,
}

impl Trial {
    fn p(&self, name: &str) -> String {
        self.dir.join(name).to_string_lossy().to_string()
    }

    fn write(&self, name: &str, text: &str) -> String {
        let p = self.p(name);
        std::fs::write(&p, text).unwrap();
        p
    }

    /// Compressed by the writer of the project (`.gz` or `.lz4`, by the
    /// name), in small writes so that the LZ4 blocks are small too.
    fn write_packed(&self, name: &str, text: &str) -> String {
        let p = self.p(name);
        let mut w = externo::Escritor::cria(&p).unwrap();
        for piece in text.as_bytes().chunks(4096) {
            w.write_all(piece).unwrap();
        }
        w.termina().unwrap();
        p
    }

    /// A gzip of many tiny members (as BGZF is, in the large). The decoder
    /// hands over one member at a time, so the counter receives the text in
    /// pieces of 1 to 97 bytes: headers, sequences, `+` lines and line breaks
    /// all get cut in the middle.
    fn write_members(&self, name: &str, text: &str, s: &mut u64) -> String {
        let p = self.p(name);
        let mut file = std::io::BufWriter::new(std::fs::File::create(&p).unwrap());
        let bytes = text.as_bytes();
        let mut i = 0;
        while i < bytes.len() {
            let end = (i + 1 + (lcg(s) % 97) as usize).min(bytes.len());
            let mut member = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            member.write_all(&bytes[i..end]).unwrap();
            file.write_all(&member.finish().unwrap()).unwrap();
            i = end;
        }
        file.flush().unwrap();
        drop(file);
        // The premise of this layout: the reader really gets small pieces.
        let mut r = externo::abre_leitura(&p).unwrap();
        let first = r.fill_buf().unwrap().len();
        assert!(first <= 97, "the gzip members are no longer delivered one by one ({first} bytes at once)");
        p
    }

    fn opt(&self, k: u8, canonical: bool, min_count: u32, mem: usize) -> CountOptions {
        CountOptions { k, canonical, min_count, memoria_bytes: mem, dir_temp: Some(self.runs_dir.clone()) }
    }

    fn temp_empty(&self) -> bool {
        std::fs::read_dir(&self.runs_dir).unwrap().next().is_none()
    }
}

fn read_text(path: &str) -> String {
    let mut r = externo::abre_leitura(path).unwrap();
    let mut s = String::new();
    r.read_to_string(&mut s).unwrap();
    s
}

fn to_lines(text: &str, what: &str) -> Vec<String> {
    assert!(text.is_empty() || text.ends_with('\n'), "{what}: the output does not end with a line break");
    assert!(!text.contains('\r'), "{what}: the output has a carriage return");
    text.split_terminator('\n').map(str::to_string).collect()
}

fn first_difference(got: &[String], want: &[String], what: &str) {
    for (i, (g, w)) in got.iter().zip(want.iter()).enumerate() {
        assert!(g == w, "{what}: line {} differs: got {g:?}, expected {w:?}", i + 1);
    }
    assert!(
        got.len() == want.len(),
        "{what}: {} lines, expected {}; the first one beyond the common part: got {:?}, expected {:?}",
        got.len(),
        want.len(),
        got.get(want.len()),
        want.get(got.len())
    );
}

fn same_counts(got: &KmerCounts, want: &KmerCounts, what: &str) {
    assert_eq!(got.k, want.k, "{what}: k");
    assert_eq!(got.canonical, want.canonical, "{what}: canonical");
    assert_eq!(got.sequences, want.sequences, "{what}: sequences");
    assert_eq!(got.bases, want.bases, "{what}: bases");
    assert_eq!(got.total_kmers, want.total_kmers, "{what}: total_kmers");
    assert_eq!(got.distinct_kmers, want.distinct_kmers, "{what}: distinct_kmers");
    assert_eq!(got.unique_kmers, want.unique_kmers, "{what}: unique_kmers");
    assert_eq!(got.max_count, want.max_count, "{what}: max_count");
    assert_eq!(got.written, want.written, "{what}: written");
    assert_eq!(got.runs, want.runs, "{what}: runs");
    for (i, (g, w)) in got.histogram.iter().zip(want.histogram.iter()).enumerate() {
        assert_eq!(g, w, "{what}: histogram entry {i}");
    }
    assert_eq!(got.histogram.len(), want.histogram.len(), "{what}: histogram entries");
    assert!(got == want, "{what}: the summaries differ");
}

/// The lines are `KMER<TAB>count`, with k letters, in strictly increasing
/// order of the value `Kmer128` gives to the k-mer.
fn increasing(lines: &[String], k: u8, what: &str) {
    let mut before: Option<u128> = None;
    for (i, l) in lines.iter().enumerate() {
        let (text, count) = l.split_once('\t').unwrap_or_else(|| panic!("{what}: line {} has no tab: {l:?}", i + 1));
        assert_eq!(text.len(), usize::from(k), "{what}: line {}: k-mer of the wrong length: {l:?}", i + 1);
        assert!(count.parse::<u32>().is_ok_and(|c| c > 0), "{what}: line {}: bad count: {l:?}", i + 1);
        let bits = kmer::Kmer128::encode(text.as_bytes()).bits;
        assert!(before.is_none_or(|b| b < bits), "{what}: line {} is out of order: {l:?}", i + 1);
        before = Some(bits);
    }
}

const BIG: usize = usize::MAX;
const SMALL: usize = 1 << 20;

/// Run files a budget gives: 32 bytes per buffered k-mer, a floor of 1 MB;
/// none when one batch holds everything.
fn runs_for(total: u64, mem: usize) -> usize {
    let per_batch = (mem.max(1 << 20) / 32) as u64;
    if total <= per_batch {
        0
    } else {
        total.div_ceil(per_batch) as usize
    }
}

fn budget_name(mem: usize) -> String {
    if mem == BIG {
        "all in memory".to_string()
    } else {
        format!("{mem} bytes")
    }
}

/// Counts with everything in memory and with the 1 MB budget; both have to be
/// the reference, line by line and field by field.
#[allow(clippy::too_many_arguments)]
fn check(t: &mut Trial, what: &str, files: &[String], records: &[String], k: u8, canonical: bool, min_count: u32, several_runs: bool) -> KmerCounts {
    check_with(t, what, files, &expected(records, k, canonical, min_count), min_count, several_runs)
}

fn check_with(t: &mut Trial, what: &str, files: &[String], want: &Expected, min_count: u32, several_runs: bool) -> KmerCounts {
    let (k, canonical) = (want.counts.k, want.counts.canonical);
    let out = t.p("counts.tsv");
    let mut texts: Vec<String> = Vec::new();
    let mut last = KmerCounts::default();
    for mem in [BIG, SMALL] {
        let ctx = format!("{what} (k={k}, canonical={canonical}, min_count={min_count}, budget: {})", budget_name(mem));
        let got = count_kmers(files, Some(&out), &t.opt(k, canonical, min_count, mem)).unwrap_or_else(|e| panic!("{ctx}: {e}"));
        let text = read_text(&out);
        let lines = to_lines(&text, &ctx);
        first_difference(&lines, &want.lines, &ctx);
        increasing(&lines, k, &ctx);
        let runs = runs_for(want.counts.total_kmers, mem);
        same_counts(&got, &KmerCounts { runs, ..want.counts.clone() }, &ctx);
        if mem == SMALL {
            assert_eq!(got.runs > 1, several_runs, "{ctx}: {} runs for {} k-mers", got.runs, got.total_kmers);
        }
        assert!(t.temp_empty(), "{ctx}: temporary files were left behind");
        texts.push(text);
        last = got;
        t.passed += 1;
    }
    assert!(texts[0] == texts[1], "{what} (k={k}, canonical={canonical}): the output from memory and the output merged from runs differ");
    last
}

/// A counting that has to fail: nothing may be left behind.
fn fails(t: &mut Trial, what: &str, files: &[String], k: u8, mem: usize) -> String {
    let out = t.p("failed.tsv");
    let _ = std::fs::remove_file(&out);
    let e = match count_kmers(files, Some(&out), &t.opt(k, false, 1, mem)) {
        Ok(r) => panic!("{what}: should have failed, but counted {} k-mers", r.total_kmers),
        Err(e) => e,
    };
    assert!(!Path::new(&out).exists(), "{what}: the failed counting left an output behind ({e})");
    assert!(t.temp_empty(), "{what}: the failed counting left temporary files behind ({e})");
    assert!(e.len() > 40, "{what}: the message should explain: {e}");
    t.passed += 1;
    e
}

fn has(e: &str, pieces: &[&str], what: &str) {
    for p in pieces {
        assert!(e.contains(p), "{what}: the message should have {p:?}: {e}");
    }
}

// ---------------------------------------------------------------------
//  Timing (not part of the trials)
// ---------------------------------------------------------------------

fn same_bytes(a: &str, b: &str) -> bool {
    let (mut x, mut y) = (std::fs::File::open(a).unwrap(), std::fs::File::open(b).unwrap());
    let (mut bx, mut by) = (vec![0u8; 1 << 20], vec![0u8; 1 << 20]);
    loop {
        let n = x.read(&mut bx).unwrap();
        if n == 0 {
            return y.read(&mut by).unwrap() == 0;
        }
        if y.read_exact(&mut by[..n]).is_err() || bx[..n] != by[..n] {
            return false;
        }
    }
}

fn bench(t: &Trial, mbp: usize, budget_mb: usize) {
    let mut s = 99u64;
    let fa = t.p("bench.fa");
    let mut w = std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(&fa).unwrap());
    let records = 4usize;
    let per_record = mbp * 1_000_000 / records;
    let mut line = [0u8; 71];
    line[70] = b'\n';
    for r in 0..records {
        writeln!(w, ">chr{}", r + 1).unwrap();
        let mut left = per_record;
        while left > 0 {
            let n = left.min(70);
            for b in line[..n].iter_mut() {
                *b = b"ACGT"[(lcg(&mut s) % 4) as usize];
            }
            w.write_all(&line[..n]).unwrap();
            w.write_all(b"\n").unwrap();
            left -= n;
        }
    }
    w.flush().unwrap();
    drop(w);

    let (k, canonical) = (31u8, true);
    let files = [fa];
    let out = t.p("bench.tsv");
    let start = std::time::Instant::now();
    let r = count_kmers(&files, Some(&out), &t.opt(k, canonical, 1, budget_mb << 20)).unwrap();
    let with_output = start.elapsed().as_secs_f64();
    assert!(t.temp_empty(), "temporary files were left behind");
    assert_eq!(r.bases, (per_record * records) as u64);
    assert_eq!(r.total_kmers, r.bases - (records as u64) * u64::from(k - 1));
    print!("{}", r.resumo());

    let start = std::time::Instant::now();
    let only = count_kmers(&files, None, &t.opt(k, canonical, 1, budget_mb << 20)).unwrap();
    let summary_only = start.elapsed().as_secs_f64();
    assert!(only == KmerCounts { written: 0, ..r.clone() }, "the summary without output differs");

    let out_memory = t.p("bench_memory.tsv");
    let start = std::time::Instant::now();
    let memory = count_kmers(&files, Some(&out_memory), &t.opt(k, canonical, 1, BIG)).unwrap();
    let in_memory = start.elapsed().as_secs_f64();
    assert!(memory == KmerCounts { runs: 0, ..r.clone() }, "the summary from memory differs");
    assert!(same_bytes(&out, &out_memory), "the table merged from runs and the table from memory differ");

    println!(
        "canonical {k}-mers of a random {mbp} Mbp FASTA, budget {budget_mb} MB: {with_output:.1} s writing the table ({} runs, {:.0} MB of output), {summary_only:.1} s for the summary alone; all in memory: {in_memory:.1} s -- same table",
        r.runs,
        std::fs::metadata(&out).unwrap().len() as f64 / 1e6
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let dir = std::env::temp_dir().join(format!("trial_kmer_count_{}", std::process::id()));
    let runs_dir = dir.join("runs");
    std::fs::create_dir_all(&runs_dir).unwrap();
    let mut t = Trial { dir: dir.clone(), runs_dir, passed: 0 };
    if args.get(1).map(String::as_str) == Some("bench") {
        let number = |i: usize, default: usize| args.get(i).and_then(|a| a.parse().ok()).unwrap_or(default);
        bench(&t, number(2, 20), number(3, 64));
        std::fs::remove_dir_all(&dir).unwrap();
        return;
    }
    let mut s = 20261113u64;

    // ------------------------------------------------ the 2-bit functions, against text
    for k in 1..=63u8 {
        for _ in 0..300 {
            let text = random_dna(usize::from(k), &mut s);
            let bits = encode_text(&text);
            assert_eq!(bits, kmer::Kmer128::encode(text.as_bytes()).bits, "the reference encoding is not the one of Kmer128: {text}");
            assert_eq!(decode(bits, k), text, "decode, k={k}");
            assert_eq!(decode(bits, k), kmer::Kmer128 { bits, k }.decode(), "decode is not the one of Kmer128, k={k}");
            let rc = rc_text(&text);
            assert_eq!(decode(reverse_complement(bits, k), k), rc, "reverse_complement of {text}");
            assert_eq!(reverse_complement(bits, k), encode_text(&rc), "reverse_complement of {text}: bits above the k-mer");
            assert_eq!(reverse_complement(reverse_complement(bits, k), k), bits, "reverse_complement twice, {text}");
            let smaller = if rc < text { &rc } else { &text };
            assert_eq!(decode(canonical(bits, k), k), *smaller, "canonical of {text}");
            assert_eq!(canonical(bits, k), canonical(encode_text(&rc), k), "canonical of the two strands of {text}");
            // With an odd k the middle base would have to be its own
            // complement: no k-mer is its own reverse complement.
            if k % 2 == 1 {
                assert_ne!(reverse_complement(bits, k), bits, "odd k, {text}");
            }
        }
        let (a, t_) = ("A".repeat(usize::from(k)), "T".repeat(usize::from(k)));
        assert_eq!(reverse_complement(encode_text(&a), k), encode_text(&t_), "poly-A, k={k}");
        assert_eq!(canonical(encode_text(&t_), k), 0, "poly-T, k={k}");
        t.passed += 1;
    }
    // Palindromes: with an even k, a half followed by its reverse complement.
    for k in (2..=62u8).step_by(2) {
        for _ in 0..100 {
            let half = random_dna(usize::from(k) / 2, &mut s);
            let text = format!("{half}{}", rc_text(&half));
            let bits = encode_text(&text);
            assert_eq!(reverse_complement(bits, k), bits, "palindrome {text}");
            assert_eq!(canonical(bits, k), bits, "palindrome {text}");
        }
        t.passed += 1;
    }
    assert_eq!((decode(0, 0), reverse_complement(77, 0), decode(0b0001_1011, 4)), (String::new(), 0, "ACGT".to_string()));
    t.passed += 1;

    // ------------------------------------------------ small cases with the answer written by hand
    {
        let table = |t: &mut Trial, file: &str, k: u8, canonical: bool, mem: usize| -> (Vec<String>, KmerCounts) {
            let out = t.p("hand.tsv");
            let r = count_kmers(&[file.to_string()], Some(&out), &t.opt(k, canonical, 1, mem)).unwrap();
            t.passed += 1;
            (to_lines(&read_text(&out), "hand"), r)
        };
        // No k-mer goes from one record into the next: together the two
        // records would also give GTAC twice, TACG and ACGT again.
        let two = t.write("two.fa", ">a\nACGTA\n>b\nCGTAC\n");
        for mem in [BIG, SMALL] {
            let (lines, r) = table(&mut t, &two, 4, false, mem);
            assert_eq!(lines, ["ACGT\t1", "CGTA\t2", "GTAC\t1"], "two records");
            assert_eq!((r.sequences, r.bases, r.total_kmers, r.distinct_kmers, r.unique_kmers, r.max_count, r.written, r.runs), (2, 10, 4, 3, 2, 2, 3, 0));
            assert_eq!(r.histogram, [(1, 2), (2, 1)]);
        }
        // The same line broken at every base is the same sequence.
        let broken = t.write("broken.fa", ">a\nA\nC\nG\nT\nA\n>b\nC\r\nG\r\nT\r\nA\r\nC");
        assert_eq!(table(&mut t, &broken, 4, false, BIG).0, ["ACGT\t1", "CGTA\t2", "GTAC\t1"], "one base per line");
        // The two strands.
        let strands = t.write("strands.fa", ">x\nAAAAC\n>y\nGTTTT\n");
        assert_eq!(table(&mut t, &strands, 5, false, BIG).0, ["AAAAC\t1", "GTTTT\t1"]);
        assert_eq!(table(&mut t, &strands, 5, true, BIG).0, ["AAAAC\t2"]);
        // FASTQ: an N, lower case, a quality line that starts with '@' and
        // is made of bases, and a header with bases in it.
        let fq = t.write("hand.fq", "@r1 ACGTACGT\nACGTN\n+\nIIIII\n@r2\nacgta\n+r2\n@+ACG\n");
        let (lines, r) = table(&mut t, &fq, 3, false, BIG);
        assert_eq!(lines, ["ACG\t2", "CGT\t2", "GTA\t1"], "FASTQ by hand");
        assert_eq!(
            r.json(),
            "{\"k\":3,\"canonical\":false,\"sequences\":2,\"bases\":10,\"total_kmers\":5,\"distinct_kmers\":3,\"unique_kmers\":1,\"max_count\":2,\"written\":3,\"runs\":0,\"histogram\":[[1,1],[2,2]]}"
        );
        let text = r.resumo();
        has(&text, &["k = 3", "one strand", "sequences read:    2", "bases read:        10", "k-mers counted:    5", "distinct k-mers:   3", "33.3%", "none"], "resumo");
        // The reference itself, on the same records.
        let by_reference = expected(&["ACGTN".to_string(), "acgta".to_string()], 3, false, 1);
        assert_eq!(by_reference.lines, lines, "the reference disagrees with the answer written by hand");
        same_counts(&r, &by_reference.counts, "FASTQ by hand");
    }

    // ------------------------------------------------ a genome and its reads: every k, both strands modes
    // Records of odd sizes, among them an empty one and some around k = 63.
    let genome: Vec<String> = [40_001usize, 35_000, 25_013, 0, 5, 62, 63, 64, 1].iter().map(|&n| random_dna(n, &mut s)).collect();
    let genome_files = [
        // lines of odd and changing lengths, most of them shorter than k
        t.write("genome.fa", &fasta_text(&genome, &[61, 7, 1, 133, 60], "\n")),
        // CRLF, delivered to the counter a few bytes at a time
        t.write_members("genome_members.fa.gz", &fasta_text(&genome, &[70], "\r\n"), &mut s),
        // each sequence on a single line
        t.write_packed("genome_one_line.fa.lz4", &fasta_text(&genome, &[usize::MAX], "\n")),
    ];
    // Reads from both strands; some shorter than any k tried, some empty.
    let source = random_dna(20_000, &mut s);
    let reads: Vec<String> = (0..1200)
        .map(|i| {
            let len = if i % 97 == 0 { (lcg(&mut s) % 20) as usize } else { 60 + (lcg(&mut s) % 120) as usize };
            let start = (lcg(&mut s) as usize) % (source.len() - len);
            let piece = &source[start..start + len];
            if lcg(&mut s) % 2 == 0 {
                piece.to_string()
            } else {
                rc_text(piece)
            }
        })
        .collect();
    let read_files = [
        t.write("reads.fq", &fastq_text(&reads, "\n", &mut s)),
        t.write_members("reads_members.fq.gz", &fastq_text(&reads, "\r\n", &mut s), &mut s),
        t.write_packed("reads.fq.lz4", &fastq_text(&reads, "\n", &mut s)),
    ];
    let mut turn = 0;
    for k in [1u8, 2, 31, 32, 33, 63] {
        for canonical in [false, true] {
            check(&mut t, "genome, FASTA", &[genome_files[turn % 3].clone()], &genome, k, canonical, 1, true);
            check(&mut t, "reads, FASTQ", &[read_files[(turn + 1) % 3].clone()], &reads, k, canonical, 1, true);
            turn += 1;
        }
    }

    // ------------------------------------------------ FASTA and FASTQ, plain and compressed, in one call
    {
        let files = vec![genome_files[0].clone(), read_files[1].clone(), genome_files[2].clone(), read_files[0].clone(), genome_files[1].clone()];
        let records: Vec<String> = [&genome[..], &reads[..], &genome[..], &reads[..], &genome[..]].concat();
        for (k, canonical) in [(21u8, true), (5, false)] {
            let r = check(&mut t, "mixed files", &files, &records, k, canonical, 1, true);
            assert_eq!(r.sequences, (3 * genome.len() + 2 * reads.len()) as u64);
        }
    }

    // ------------------------------------------------ lower case, runs of N, IUPAC letters, gaps
    {
        let dirty_genome: Vec<String> = [60_000usize, 40_000, 3_000].iter().map(|&n| dirty(&random_dna(n, &mut s), &mut s)).collect();
        let fa = t.write("dirty.fa", &fasta_text(&dirty_genome, &[59, 3, 80], "\n"));
        let fq = t.write_members("dirty.fq.gz", &fastq_text(&dirty_genome, "\n", &mut s), &mut s);
        assert!(dirty_genome[0].contains(&"N".repeat(31)) && dirty_genome[0].contains('r') && dirty_genome[0].contains("a"));
        for (k, canonical) in [(31u8, false), (31, true), (4, true), (63, false)] {
            check(&mut t, "lower case, N and IUPAC, FASTA", &[fa.clone()], &dirty_genome, k, canonical, 1, true);
            check(&mut t, "lower case, N and IUPAC, FASTQ", &[fq.clone()], &dirty_genome, k, canonical, 1, true);
        }
        // All in lower case is the same as all in upper case.
        let lower: Vec<String> = genome.iter().map(|r| r.to_ascii_lowercase()).collect();
        let fa = t.write("lower.fa", &fasta_text(&lower, &[60], "\n"));
        let a = check(&mut t, "all in lower case", &[fa], &lower, 31, true, 1, true);
        let b = check(&mut t, "all in upper case", &[genome_files[0].clone()], &genome, 31, true, 1, true);
        assert!(a == b, "lower and upper case give different summaries");
    }

    // ------------------------------------------------ sequences shorter than k
    {
        let short: Vec<String> = (0..40).map(|n| random_dna(n, &mut s)).collect(); // 0 to 39 bases
        let fa = t.write("short.fa", &fasta_text(&short, &[11], "\n"));
        let fq = t.write("short.fq", &fastq_text(&short, "\n", &mut s));
        for k in [31u8, 39, 40, 63] {
            for file in [&fa, &fq] {
                let r = check(&mut t, "sequences shorter than k", &[file.clone()], &short, k, true, 1, false);
                assert_eq!(r.total_kmers, (0..40u64).map(|n| (n + 1).saturating_sub(u64::from(k))).sum::<u64>());
                assert_eq!((r.sequences, r.runs), (40, 0));
            }
        }
    }

    // ------------------------------------------------ highly repetitive: large counts, the same k-mer in every run
    {
        let rep: Vec<String> =
            vec!["A".repeat(60_000), "ACGT".repeat(10_000), "T".repeat(30_000), "AC".repeat(5_000), "ACGGT".repeat(2_000), random_dna(5_000, &mut s)];
        let fa = t.write("repeats.fa", &fasta_text(&rep, &[usize::MAX, 60, 1000, 7], "\n"));
        let gz = t.write_members("repeats.fa.gz", &fasta_text(&rep, &[73], "\n"), &mut s);
        // The highest count, by hand, where the random record cannot add to
        // it: the run of 60,000 A gives 60,000 - (k - 1) times the k-mer
        // made of A, and with `canonical` the run of 30,000 T joins it.
        for (k, canonical, top) in [(31u8, false, Some(59_970u32)), (31, true, Some(89_940)), (4, false, None), (2, true, None), (63, true, Some(89_876))] {
            for file in [&fa, &gz] {
                let r = check(&mut t, "repeats", &[file.clone()], &rep, k, canonical, 1, true);
                assert!(top.is_none_or(|top| r.max_count == top), "repeats: the highest count is {}, k={k} canonical={canonical}", r.max_count);
                assert!(r.max_count > 50_000, "repeats: the highest count is {}, k={k}", r.max_count);
                assert_eq!(r.histogram.last().map(|h| h.0), Some(10_001), "repeats: everything above 10,000 goes in one last entry");
            }
        }
    }

    // ------------------------------------------------ min_count
    {
        let file = [read_files[0].clone()];
        let mut written = Vec::new();
        let counted = reference(&reads, 15, true);
        for min_count in [0u32, 1, 4, 8, 12, 1_000_000] {
            let want = expected_of(&counted, &reads, 15, true, min_count);
            let r = check_with(&mut t, "min_count", &file, &want, min_count, true);
            written.push(r.written);
        }
        assert_eq!(written[0], written[1], "min_count 0 is the same as 1");
        assert!(written[1] > written[2] && written[2] > written[3] && written[3] > written[4] && written[4] > 0 && written[5] == 0, "min_count: {written:?}");
    }

    // ------------------------------------------------ no output: only the summary
    {
        let want = expected(&reads, 31, true, 1).counts;
        let out = t.p("counts.tsv");
        std::fs::remove_file(&out).unwrap();
        let before = std::fs::read_dir(&t.dir).unwrap().count();
        for mem in [BIG, SMALL] {
            let got = count_kmers(&[read_files[2].clone()], None, &t.opt(31, true, 1, mem)).unwrap();
            same_counts(&got, &KmerCounts { written: 0, runs: runs_for(want.total_kmers, mem), ..want.clone() }, "no output");
            assert!(got.runs > 1 || mem == BIG);
            assert!(t.temp_empty() && std::fs::read_dir(&t.dir).unwrap().count() == before, "no output: something was written");
            t.passed += 1;
        }
    }

    // ------------------------------------------------ compressed output, read back
    {
        let want = expected(&genome, 31, true, 2);
        let whole = expected(&genome, 31, true, 1);
        for (name, magic) in [("counts.tsv.gz", &[0x1Fu8, 0x8B][..]), ("counts.tsv.lz4", &[0x04u8, 0x22, 0x4D, 0x18][..])] {
            for (mem, min_count, want) in [(BIG, 2, &want), (SMALL, 1, &whole)] {
                let out = t.p(name);
                let got = count_kmers(&[genome_files[0].clone()], Some(&out), &t.opt(31, true, min_count, mem)).unwrap();
                assert!(std::fs::read(&out).unwrap().starts_with(magic), "{name} is not compressed");
                first_difference(&to_lines(&read_text(&out), name), &want.lines, name);
                same_counts(&got, &KmerCounts { runs: runs_for(want.counts.total_kmers, mem), ..want.counts.clone() }, name);
                t.passed += 1;
            }
        }
    }

    // ------------------------------------------------ CRLF, blanks, comments, a byte-order mark, no final line break
    {
        let mut fa = String::from("\u{feff}\r\n  \r\n");
        for (i, r) in genome.iter().enumerate() {
            fa.push_str(&format!(">s{i} {BAIT}\r\n"));
            if i % 2 == 0 {
                fa.push_str(&format!(";a comment in the old style {BAIT}\r\n"));
            }
            let mut rest = r.as_str();
            let mut n_line = 0;
            while !rest.is_empty() {
                let n = [59usize, 3, 64][n_line % 3].min(rest.len());
                // blanks before and after the bases, and an empty line now and then
                let (lead, trail) = [("", ""), ("", " "), ("  ", "\t"), ("\t", "  \t ")][n_line % 4];
                fa.push_str(&format!("{lead}{}{trail}\r\n", &rest[..n]));
                if n_line % 5 == 4 {
                    fa.push_str(if n_line % 2 == 0 { "\r\n" } else { " \t\r\n" });
                }
                rest = &rest[n..];
                n_line += 1;
            }
        }
        let fa = fa.trim_end().to_string();
        assert!(!fa.ends_with('\n') && fa.contains("\r\n"));
        let plain = t.write("messy.fa", &fa);
        let members = t.write_members("messy.fa.gz", &fa, &mut s);
        // The file ends on the only base of the last record, with no line break.
        for (file, k) in [(&plain, 31u8), (&members, 31), (&plain, 63), (&members, 2)] {
            check(&mut t, "FASTA with CRLF, blanks, comments and a byte-order mark", &[file.clone()], &genome, k, true, 1, true);
        }

        let mut fq = fastq_text(&reads, "\r\n", &mut s).replace("\r\n@read", "\r\n\r\n \r\n@read");
        fq.truncate(fq.len() - 2);
        assert!(!fq.ends_with('\n') && fq.contains("\r\n\r\n \r\n@read7 "));
        let plain = t.write("messy.fq", &fq);
        let members = t.write_members("messy.fq.gz", &fq, &mut s);
        for (file, k) in [(&plain, 31u8), (&members, 31), (&members, 63), (&plain, 1)] {
            check(&mut t, "FASTQ with CRLF, empty lines between records and no final line break", &[file.clone()], &reads, k, false, 1, true);
        }
    }

    // ------------------------------------------------ a file and the file of its reverse complements
    {
        let forward: Vec<String> = genome.iter().chain(reads.iter().take(300)).map(|r| dirty(r, &mut s)).collect();
        let backward: Vec<String> = forward.iter().map(|r| rc_text(r)).collect();
        let f = t.write("forward.fa", &fasta_text(&forward, &[60], "\n"));
        let b = t.write("backward.fa", &fasta_text(&backward, &[60], "\n"));
        let (out_f, out_b) = (t.p("forward.tsv"), t.p("backward.tsv"));
        for k in [4u8, 31, 32, 63] {
            for mem in [BIG, SMALL] {
                let rf = count_kmers(&[f.clone()], Some(&out_f), &t.opt(k, true, 1, mem)).unwrap();
                let rb = count_kmers(&[b.clone()], Some(&out_b), &t.opt(k, true, 1, mem)).unwrap();
                assert!(rf == rb, "reverse complement, k={k}: the summaries differ");
                assert!(rf.total_kmers > 70_000 && (rf.runs > 1) == (mem == SMALL));
                assert!(read_text(&out_f) == read_text(&out_b), "reverse complement, k={k}: the canonical tables differ");
                t.passed += 1;
            }
        }
        // Without `canonical` the strands are told apart.
        count_kmers(&[f], Some(&out_f), &t.opt(31, false, 1, BIG)).unwrap();
        count_kmers(&[b], Some(&out_b), &t.opt(31, false, 1, BIG)).unwrap();
        assert!(read_text(&out_f) != read_text(&out_b), "the two strands gave the same table without `canonical`");
        t.passed += 1;
    }

    // ------------------------------------------------ nothing to count is not an error
    {
        let zero = |sequences: u64| KmerCounts { k: 31, canonical: true, sequences, ..KmerCounts::default() };
        let cases: Vec<(Vec<String>, u64)> = vec![
            (vec![], 0),
            (vec![t.write("empty.fa", "")], 0),
            (vec![t.write("blank.fa", "\n\r\n   \n\t")], 0),
            (vec![t.write_packed("empty.fq.gz", "")], 0),
            (vec![t.write("headers.fa", ">a\n>b\n>c")], 3),
            (vec![t.write("empty_reads.fq", "@r1\n\n+\n\n@r2\n\n+\n")], 2),
            (vec![t.p("empty.fa"), t.p("headers.fa"), t.p("blank.fa"), t.p("empty_reads.fq")], 5),
        ];
        let out = t.p("nothing.tsv");
        for (files, sequences) in cases {
            for mem in [BIG, SMALL] {
                let got = count_kmers(&files, Some(&out), &t.opt(31, true, 1, mem)).unwrap_or_else(|e| panic!("{files:?}: {e}"));
                same_counts(&got, &zero(sequences), &format!("{files:?}"));
                assert_eq!(std::fs::metadata(&out).unwrap().len(), 0, "{files:?}: the output should be an empty file");
                assert!(got.resumo().contains("no k-mer was counted"));
                t.passed += 1;
            }
        }
        assert_eq!(zero(0).json(), "{\"k\":31,\"canonical\":true,\"sequences\":0,\"bases\":0,\"total_kmers\":0,\"distinct_kmers\":0,\"unique_kmers\":0,\"max_count\":0,\"written\":0,\"runs\":0,\"histogram\":[]}");
    }

    // ------------------------------------------------ the budget: 32 bytes per k-mer, a floor of 1 MB
    {
        let want = expected(&genome, 31, false, 1);
        let total = want.counts.total_kmers;
        assert_eq!(total, 100_023);
        // The trial data itself: a genome that is really random has no
        // 31-mer twice (a generator with a short period would repeat them).
        assert_eq!(want.counts.distinct_kmers, total, "the random genome repeats itself");
        let out = t.p("budget.tsv");
        for (mem, runs) in [(0usize, 4usize), (1, 4), (SMALL, 4), (1_572_864, 3), (2 << 20, 2), (total as usize * 32 - 1, 2), (total as usize * 32, 0), (1 << 30, 0)] {
            let got = count_kmers(&[genome_files[0].clone()], Some(&out), &t.opt(31, false, 1, mem)).unwrap();
            assert_eq!(runs, runs_for(total, mem), "the trial itself: runs expected for {mem} bytes");
            same_counts(&got, &KmerCounts { runs, ..want.counts.clone() }, &format!("budget of {mem} bytes"));
            first_difference(&to_lines(&read_text(&out), "budget"), &want.lines, &format!("budget of {mem} bytes"));
            assert!(t.temp_empty());
            t.passed += 1;
        }
        // Without a folder for the runs, the temporary folder of the system.
        let got = count_kmers(
            &[genome_files[1].clone()],
            Some(&out),
            &CountOptions { k: 31, canonical: false, min_count: 1, memoria_bytes: SMALL, dir_temp: None },
        )
        .unwrap();
        same_counts(&got, &KmerCounts { runs: 4, ..want.counts.clone() }, "runs in the temporary folder of the system");
        first_difference(&to_lines(&read_text(&out), "budget"), &want.lines, "runs in the temporary folder of the system");
        assert!(got.resumo().contains("temporary runs:    4") && got.json().contains("\"runs\":4,"));
        t.passed += 1;
    }

    // ------------------------------------------------ counts saturate at u32::MAX (runs written by hand)
    {
        let runs: Vec<String> = (0..3).map(|i| t.p(&format!("hand_run_{i}.lz4"))).collect();
        let pairs: [&[(u128, u32)]; 3] = [
            &[(encode_text("ACC"), u32::MAX - 3), (encode_text("AGC"), 7)],
            &[(encode_text("ACC"), 10), (encode_text("AGC"), 1), (encode_text("ATA"), u32::MAX)],
            &[(encode_text("ATA"), u32::MAX), (encode_text("ATC"), 1), (encode_text("TTT"), u32::MAX - 1)],
        ];
        for (path, pairs) in runs.iter().zip(pairs) {
            let mut w = RunWriter::create(path).unwrap();
            for &(kmer, count) in pairs {
                w.pair(kmer, count).unwrap();
            }
            w.finish().unwrap();
        }
        let out = t.p("saturated.tsv");
        let mut collector = Collector::new(3, 1, Some(&out)).unwrap();
        merge_runs(&runs, &mut collector).unwrap();
        let r = collector.finish().unwrap();
        assert_eq!(to_lines(&read_text(&out), "saturated"), ["ACC\t4294967295", "AGC\t8", "ATA\t4294967295", "ATC\t1", "TTT\t4294967294"]);
        assert_eq!((r.distinct_kmers, r.unique_kmers, r.max_count, r.written), (5, 1, u32::MAX, 5));
        assert_eq!(r.histogram, [(1, 1), (8, 1), (10_001, 3)]);
        assert!(KmerCounts { total_kmers: 9, ..r }.resumo().contains("4,294,967,295 or more"));
        t.passed += 1;
    }

    // ------------------------------------------------ errors
    {
        let good = genome_files[0].clone();
        for k in [0u8, 64, 255] {
            let e = fails(&mut t, "k out of range", &[good.clone()], k, BIG);
            has(&e, &[&format!("k = {k}"), "from 1 to 63"], "k out of range");
        }

        // Neither FASTA nor FASTQ.
        let vcf = t.write("variants.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tT\t30\tPASS\t.\n");
        let e = fails(&mut t, "a VCF", &[vcf.clone()], 31, BIG);
        has(&e, &[&vcf, "line 1", "##fileformat=VCFv4.2", "FASTA", "FASTQ", "'>'", "'@'"], "a VCF");
        let late = t.write("late.txt", "\n\r\n  \nACGTACGT\n");
        let e = fails(&mut t, "bases without a header", &[late.clone()], 31, BIG);
        has(&e, &[&late, "line 4", "ACGTACGT"], "bases without a header");
        // After good files, and after runs were already on disk.
        let e = fails(&mut t, "a VCF after two good files", &[good.clone(), read_files[0].clone(), vcf.clone()], 31, SMALL);
        has(&e, &[&vcf, "line 1"], "a VCF after two good files");
        // Compressed content under a name that does not say so.
        let hidden = t.p("hidden_gzip.fa");
        std::fs::copy(&genome_files[1], &hidden).unwrap();
        let e = fails(&mut t, "gzip without .gz", &[hidden.clone()], 31, BIG);
        has(&e, &[&hidden, "line 1", "gzip", ".gz"], "gzip without .gz");

        // FASTQ cut in the middle of the second record, which starts at line 5.
        let whole = "@r1\nACGT\n+\nIIII\n@r2\nGGCCA\n+\nIIIII\n";
        for cut in ["@r2", "@r2\n", "@r2\nGG", "@r2\nGGCCA\n", "@r2\nGGCCA\n+", "@r2\nGGCCA\n+\n", "@r2\nGGCCA\n+\nII", "@r2\nGGCCA\n+\nIIII"] {
            let text = format!("@r1\nACGT\n+\nIIII\n{cut}");
            assert!(whole.starts_with(&text) && text != whole);
            let plain = t.write("cut.fq", &text);
            let e = fails(&mut t, "truncated FASTQ", &[plain.clone()], 3, BIG);
            has(&e, &[&plain, "truncated", "line 5"], &format!("truncated FASTQ {cut:?}"));
            let packed = t.write_packed("cut.fq.gz", &text);
            let e = fails(&mut t, "truncated FASTQ, compressed", &[good.clone(), packed.clone()], 3, SMALL);
            has(&e, &[&packed, "truncated", "line 5"], &format!("truncated FASTQ, compressed {cut:?}"));
        }
        // The whole one is fine, with or without the last line break.
        for text in [whole.to_string(), whole.trim_end().to_string()] {
            let r = count_kmers(&[t.write("whole.fq", &text)], None, &t.opt(3, false, 1, BIG)).unwrap();
            assert_eq!((r.sequences, r.bases, r.total_kmers), (2, 9, 5));
            t.passed += 1;
        }

        // FASTQ that is not made of 4-line records.
        let wrapped = t.write("wrapped.fq", "@r1\nACGT\nACGT\n+\nIIIIIIII\n");
        let e = fails(&mut t, "FASTQ with a wrapped sequence", &[wrapped.clone()], 3, BIG);
        has(&e, &[&wrapped, "line 3", "'+'", "ACGT"], "FASTQ with a wrapped sequence");
        let hole = t.write("hole.fq", "@r1\nACGT\n\nIIII\n");
        let e = fails(&mut t, "FASTQ with an empty line for '+'", &[hole.clone()], 3, BIG);
        has(&e, &[&hole, "line 3", "'+'"], "FASTQ with an empty line for '+'");
        let uneven = t.write("uneven.fq", "@r1\nACGT\n+\nIII\n@r2\nAC\n+\nII\n");
        let e = fails(&mut t, "FASTQ with fewer qualities than bases", &[uneven.clone()], 3, BIG);
        has(&e, &[&uneven, "line 4", "3 characters", "4 bases"], "FASTQ with fewer qualities than bases");
        let stray = t.write("stray.fq", "@r1\nACGT\n+\nIIII\nACGT\n");
        let e = fails(&mut t, "FASTQ with a stray line", &[stray.clone()], 3, BIG);
        has(&e, &[&stray, "line 5", "'@'"], "FASTQ with a stray line");

        // A compressed file cut in the middle.
        let bytes = std::fs::read(t.write_packed("to_cut.fa.gz", &fasta_text(&genome, &[60], "\n"))).unwrap();
        let cut = t.p("cut.fa.gz");
        std::fs::write(&cut, &bytes[..bytes.len() / 2]).unwrap();
        let e = fails(&mut t, "gzip cut in the middle", &[cut.clone()], 31, SMALL);
        has(&e, &[&cut, "Error reading"], "gzip cut in the middle");

        // A file that is not there: found before any work, and before the
        // output of an earlier counting is touched.
        let missing = t.p("not_here.fa");
        let e = fails(&mut t, "missing file", &[missing.clone()], 31, BIG);
        has(&e, &[&missing, "Cannot open", "does not exist"], "missing file");
        let kept = t.write("kept.tsv", "the table of an earlier counting\n");
        let e = count_kmers(&[good.clone(), missing.clone()], Some(&kept), &t.opt(31, false, 1, SMALL)).unwrap_err();
        has(&e, &[&missing], "missing second file");
        assert_eq!(read_text(&kept), "the table of an earlier counting\n", "a missing input erased an existing output");
        t.passed += 1;

        // The output over one of the inputs.
        let precious = t.write("precious.fa", ">a\nACGTACGTACGT\n");
        let e = count_kmers(&[good.clone(), precious.clone()], Some(&precious), &t.opt(4, false, 1, BIG)).unwrap_err();
        has(&e, &[&precious, "input"], "output over an input");
        assert_eq!(read_text(&precious), ">a\nACGTACGTACGT\n", "the input was overwritten");
        t.passed += 1;

        // An output that cannot be created.
        let nowhere = t.dir.join("no_such_folder").join("out.tsv").to_string_lossy().to_string();
        let e = count_kmers(&[good.clone()], Some(&nowhere), &t.opt(4, false, 1, BIG)).unwrap_err();
        has(&e, &[&nowhere, "output folder"], "output in a folder that does not exist");
        t.passed += 1;

        // A folder for the runs that cannot be used: an error when runs are
        // needed, and no matter at all when everything fits in memory.
        let not_a_folder = PathBuf::from(t.write("not_a_folder", "x"));
        let out = t.p("failed.tsv");
        let with = |mem: usize| CountOptions { k: 31, canonical: true, min_count: 1, memoria_bytes: mem, dir_temp: Some(not_a_folder.clone()) };
        let e = count_kmers(&[good.clone()], Some(&out), &with(SMALL)).unwrap_err();
        has(&e, &["not_a_folder", "dir_temp"], "unusable folder for the runs");
        assert!(!Path::new(&out).exists(), "the failed counting left an output behind");
        let r = count_kmers(&[good.clone()], Some(&out), &with(BIG)).unwrap();
        assert_eq!((r.runs, r.total_kmers), (0, 100_023));
        t.passed += 2;
    }

    std::fs::remove_dir_all(&dir).unwrap();
    println!("{} k-mer counting tests passed", t.passed);
}
