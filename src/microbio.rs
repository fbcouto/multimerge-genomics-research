// =====================================================================
//  CLINICAL MICROBIOLOGY: WHICH SPECIES, HOW MUCH OF EACH, AND WHICH
//  RESISTANCE GENES
//
//  Four tools; every one streams FASTA/FASTQ of any size (plain, `.gz`,
//  BGZF, `.lz4`):
//
//   1. `sketch_references`   FracMinHash sketch of each reference genome,
//                            saved once in a sketch database file;
//   2. `gather`              which references are in a sample and how much
//                            of it each one explains -- the greedy minimum
//                            set cover of `sourmash gather`;
//   3. `classify_reads`      each read (or pair) to the reference it shares
//                            most sampled k-mers with, as Kraken does but
//                            without a taxonomy tree; host reads apart, and
//                            the host-free reads written out if asked;
//   4. `find_resistance_genes`  the NCBI AMRFinderPlus nucleotide catalog
//                            (AMR_CDS.fa + ReferenceGeneCatalog.txt) against
//                            an assembly (contigs) or reads.
//
//  HASHING -- exactly sourmash's, so that the numbers can be checked against
//  it. The canonical k-mer is the smaller, as text, of the k-mer (upper case)
//  and its reverse complement; its hash is the first 64-bit word of
//  MurmurHash3 x64_128 with seed 42; the k-mer enters a sketch when
//  hash <= max_hash, with max_hash = (2^64 - 1) / scaled computed in f64 and
//  truncated, as sourmash does. A k-mer holding any letter other than A, C,
//  G or T (either case) is skipped. A genome gives here the same hashes as
//  `sourmash sketch dna -p k=31,scaled=1000` gives it.
//
//  FracMinHash is NESTED: the sketch at scaled 1000 is the part of the sketch
//  at scaled 10 whose hashes are <= max_hash(1000). One database built with a
//  small scaled therefore serves any larger one: `gather` and
//  `classify_reads` downsample it when asked (`scaled`).
//
//  READ CLASSIFICATION USES THE SKETCH HASHES, NOT ALL K-MERS (the choice the
//  task left open). A read is classified by the sampled k-mers it holds, so
//  the scaled of the database decides the sensitivity: a 150-base read holds
//  about 120 / scaled of them. scaled = 1000 (the gather default) leaves
//  almost every read without one; scaled = 10 gives ~12 per read and ~24 per
//  pair; scaled = 1 is "all k-mers", for small reference sets (each k-mer
//  costs 12 bytes of memory while classifying). One code path, one file
//  format, and the memory is chosen by the user through scaled.
//
//  THE SKETCH DATABASE FILE (version 1; integers little-endian)
//
//      offset size  field
//      0      8     magic "MMGSKDB1"
//      8      1     k
//      9      1     molecule: 0 = DNA, canonical k-mers
//      10     2     reserved, 0
//      12     4     hash seed (42)
//      16     8     scaled
//      24     8     max_hash
//      32     4     number of sketches N
//      then N sketches, each:
//             1     flags: bit 0 = host decoy
//             4+n   name, UTF-8 (length, then bytes)
//             4+n   source files, joined by '\n'
//             8     sequences    8  bases    8  valid k-mers
//             8     number of hashes H
//             ...   the H hashes in increasing order, as LEB128 varints of
//                   the differences (the first one is the hash itself)
//      then   8     "MMGSKEND"
//             8     FNV-1a 64 of every byte before "MMGSKEND"
//
//  At scaled 1000 a bacterial genome takes ~30 KB.
//
//  RESISTANCE GENES. Each catalog gene is seeded by its shared 21-mers with
//  the contig (or read); the seeds of one gene on one strand are grouped by
//  diagonal, and the gene is aligned to that contig region by a banded local
//  alignment with affine gaps (match +2, mismatch -3, gap 5 + 2 per base, the
//  blastn scores). Identity = identical columns / alignment columns (BLAST's
//  pident); coverage = aligned gene bases / gene length. On an assembly the
//  best allele of each locus is kept (hits that overlap by half the shorter
//  one are one locus). On reads, every read is aligned to each catalog gene
//  it seeds well (at least half the seeds of its best gene) -- end to end,
//  except where it runs past the gene, so that an indel shows instead of
//  being clipped away --, the alignments are piled up per gene, and the
//  consensus gives identity and coverage.
//  Without positions, a locus is told by its reads. A gene half of whose
//  reads also aligned to a better gene already reported is that locus seen
//  through a relative -- or a novel allele that several catalog alleles each
//  match in part -- and gives no line of its own; except as the second
//  allele of the group, when the reported gene's close reads disagree at
//  some position (two alleles are really there, `mixed`) and some reads
//  match the second allele best. At most two alleles per group of
//  near-identical genes: three copies need an assembly. (The catalog's
//  "family" field cannot group: aac(6')-Ib', -AKT, -cr, ... are separate
//  families of near-identical sequences.)
// =====================================================================

use crate::externo::{abre_leitura, Escritor};
use crate::multimerge_core::{multi_merge_sort_stable, multi_merge_sort_stable_by};
use rayon::prelude::*;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
use std::fs::File;
use std::io::{BufRead, ErrorKind, Write};
use std::sync::mpsc::{sync_channel, SyncSender};

// ---------------------------------------------------------------------
//  MURMURHASH3 x64_128 (the first 64-bit word), as sourmash uses it
// ---------------------------------------------------------------------

/// The seed sourmash uses.
pub const SEED: u64 = 42;

#[inline]
fn fmix64(mut k: u64) -> u64 {
    k ^= k >> 33;
    k = k.wrapping_mul(0xff51_afd7_ed55_8ccd);
    k ^= k >> 33;
    k = k.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    k ^= k >> 33;
    k
}

/// MurmurHash3 x64_128 of `data`, first 64-bit word: `mmh3.hash64(data,
/// seed)[0]` in Python, `hash_murmur` in sourmash.
pub fn murmur64(data: &[u8], seed: u64) -> u64 {
    const C1: u64 = 0x87c3_7b91_1142_53d5;
    const C2: u64 = 0x4cf5_ad43_2745_937f;
    let len = data.len();
    let blocks = len / 16;
    let (mut h1, mut h2) = (seed, seed);
    for b in 0..blocks {
        let o = b * 16;
        let mut k1 = u64::from_le_bytes(data[o..o + 8].try_into().unwrap());
        let mut k2 = u64::from_le_bytes(data[o + 8..o + 16].try_into().unwrap());
        k1 = k1.wrapping_mul(C1).rotate_left(31).wrapping_mul(C2);
        h1 ^= k1;
        h1 = h1.rotate_left(27).wrapping_add(h2).wrapping_mul(5).wrapping_add(0x52dc_e729);
        k2 = k2.wrapping_mul(C2).rotate_left(33).wrapping_mul(C1);
        h2 ^= k2;
        h2 = h2.rotate_left(31).wrapping_add(h1).wrapping_mul(5).wrapping_add(0x3849_5ab5);
    }
    let tail = &data[blocks * 16..];
    let rem = tail.len();
    if rem > 8 {
        let mut k2 = 0u64;
        for (i, &b) in tail[8..].iter().enumerate() {
            k2 ^= u64::from(b) << (8 * i);
        }
        h2 ^= k2.wrapping_mul(C2).rotate_left(33).wrapping_mul(C1);
    }
    if rem > 0 {
        let mut k1 = 0u64;
        for (i, &b) in tail[..rem.min(8)].iter().enumerate() {
            k1 ^= u64::from(b) << (8 * i);
        }
        h1 ^= k1.wrapping_mul(C1).rotate_left(31).wrapping_mul(C2);
    }
    h1 ^= len as u64;
    h2 ^= len as u64;
    h1 = h1.wrapping_add(h2);
    h2 = h2.wrapping_add(h1);
    h1 = fmix64(h1);
    h2 = fmix64(h2);
    h1.wrapping_add(h2)
}

/// The largest hash kept at a given scaled, as sourmash computes it.
pub fn max_hash_for_scaled(scaled: u64) -> u64 {
    match scaled {
        0 => 0,
        1 => u64::MAX,
        s => (u64::MAX as f64 / s as f64) as u64,
    }
}

// ---------------------------------------------------------------------
//  BASES
// ---------------------------------------------------------------------

/// 0..3 for A, C, G, T (either case); 4 for anything else.
const CODE: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
};

/// The base in upper case; anything else becomes `N`.
const UPPER: [u8; 256] = {
    let mut t = [b'N'; 256];
    t[b'A' as usize] = b'A';
    t[b'a' as usize] = b'A';
    t[b'C' as usize] = b'C';
    t[b'c' as usize] = b'C';
    t[b'G' as usize] = b'G';
    t[b'g' as usize] = b'G';
    t[b'T' as usize] = b'T';
    t[b't' as usize] = b'T';
    t
};

/// The complement of an upper-case base; `N` stays `N`.
const COMP: [u8; 256] = {
    let mut t = [b'N'; 256];
    t[b'A' as usize] = b'T';
    t[b'C' as usize] = b'G';
    t[b'G' as usize] = b'C';
    t[b'T' as usize] = b'A';
    t
};

/// Bit mask of the bases a letter of a CATALOG gene stands for (IUPAC):
/// A=1, C=2, G=4, T=8; N and the ambiguity codes cover several.
const GENE_MASK: [u8; 256] = {
    let mut t = [0u8; 256];
    let codes: [(u8, u8); 15] = [
        (b'A', 1),
        (b'C', 2),
        (b'G', 4),
        (b'T', 8),
        (b'U', 8),
        (b'R', 5),
        (b'Y', 10),
        (b'S', 6),
        (b'W', 9),
        (b'K', 12),
        (b'M', 3),
        (b'B', 14),
        (b'D', 13),
        (b'H', 11),
        (b'V', 7),
    ];
    let mut i = 0;
    while i < codes.len() {
        t[codes[i].0 as usize] = codes[i].1;
        t[codes[i].0.to_ascii_lowercase() as usize] = codes[i].1;
        i += 1;
    }
    t[b'N' as usize] = 15;
    t[b'n' as usize] = 15;
    t
};

/// Bit mask of a base of a contig or read: only A, C, G, T count; an `N` in
/// the sample never matches anything (a scaffold gap is not evidence).
const SAMPLE_MASK: [u8; 256] = {
    let mut t = [0u8; 256];
    t[b'A' as usize] = 1;
    t[b'a' as usize] = 1;
    t[b'C' as usize] = 2;
    t[b'c' as usize] = 2;
    t[b'G' as usize] = 4;
    t[b'g' as usize] = 4;
    t[b'T' as usize] = 8;
    t[b't' as usize] = 8;
    t
};

// ---------------------------------------------------------------------
//  SMALL HELPERS: messages, JSON, numbers
// ---------------------------------------------------------------------

fn check_input(path: &str) -> Result<(), String> {
    File::open(path).map(|_| ()).map_err(|e| {
        let reason = match e.kind() {
            ErrorKind::NotFound => "the file does not exist".to_string(),
            ErrorKind::PermissionDenied => "permission denied".to_string(),
            _ => e.to_string(),
        };
        format!("Cannot open the input file {path}: {reason}. Check the path and that you are allowed to read the file")
    })
}

fn same_file(a: &str, b: &str) -> bool {
    match (std::fs::canonicalize(a), std::fs::canonicalize(b)) {
        (Ok(x), Ok(y)) => x == y,
        _ => a == b,
    }
}

/// The output must not be one of the inputs: it would be erased before being
/// read.
fn not_an_input(output: &str, inputs: &[&str]) -> Result<(), String> {
    for i in inputs {
        if same_file(i, output) {
            return Err(format!(
                "The output {output} is also one of the input files, and would be erased before being read. Choose another name for the output"
            ));
        }
    }
    Ok(())
}

fn create_output(path: &str) -> Result<Escritor, String> {
    Escritor::cria(path).map_err(|e| format!("{e}. Check that the output folder exists and that you can write to it"))
}

fn write_error(path: &str, e: impl std::fmt::Display) -> String {
    format!("Error writing the output {path}: {e}. Check the free space and the permission to write there")
}

/// Removes the outputs written so far unless disarmed: half a table would
/// pass for a whole one.
struct Cleanup(Vec<String>);

impl Cleanup {
    fn disarm(&mut self) {
        self.0.clear();
    }
}

impl Drop for Cleanup {
    fn drop(&mut self) {
        for p in &self.0 {
            let _ = std::fs::remove_file(p);
        }
    }
}

fn js(s: &str) -> String {
    let mut o = String::with_capacity(s.len() + 2);
    o.push('"');
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            '\r' => o.push_str("\\r"),
            '\t' => o.push_str("\\t"),
            c if (c as u32) < 0x20 => o.push_str(&format!("\\u{:04x}", c as u32)),
            c => o.push(c),
        }
    }
    o.push('"');
    o
}

fn jf(x: f64) -> String {
    if x.is_finite() {
        format!("{x}")
    } else {
        "null".to_string()
    }
}

fn jo(x: Option<f64>) -> String {
    x.map_or_else(|| "null".to_string(), jf)
}

fn jlist(items: &[String]) -> String {
    let v: Vec<String> = items.iter().map(|s| js(s)).collect();
    format!("[{}]", v.join(","))
}

fn thousands(n: u64) -> String {
    let digits = n.to_string();
    let mut o = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            o.push(',');
        }
        o.push(c);
    }
    o
}

fn pct(x: f64) -> String {
    let p = 100.0 * x;
    if p == 0.0 {
        "0%".to_string()
    } else if p.abs() >= 10.0 {
        format!("{p:.1}%")
    } else if p.abs() >= 0.1 {
        format!("{p:.2}%")
    } else {
        format!("{p:.4}%")
    }
}

/// The beginning of a line, safe to show in a message.
fn sample_text(piece: &[u8]) -> String {
    let mut o = String::new();
    for &b in piece.iter().take(40) {
        match b {
            b'"' | b'\\' => {
                o.push('\\');
                o.push(char::from(b));
            }
            0x20..=0x7E => o.push(char::from(b)),
            _ => o.push_str(&format!("\\x{b:02X}")),
        }
    }
    if piece.len() > 40 {
        o.push_str("...");
    }
    o
}

fn trim(b: &[u8]) -> &[u8] {
    let s = b.iter().position(|c| !c.is_ascii_whitespace()).unwrap_or(b.len());
    let e = b.iter().rposition(|c| !c.is_ascii_whitespace()).map_or(s, |p| p + 1);
    &b[s..e.max(s)]
}

/// The read name: the first word of the header, without a trailing `/1` or
/// `/2`.
fn read_id(head: &[u8]) -> &[u8] {
    let first = head.split(|c| c.is_ascii_whitespace()).next().unwrap_or(head);
    if first.len() > 2 && first[first.len() - 2] == b'/' && matches!(first[first.len() - 1], b'1' | b'2') {
        &first[..first.len() - 2]
    } else {
        first
    }
}

// ---------------------------------------------------------------------
//  READING FASTA AND FASTQ, RECORD BY RECORD
// ---------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SeqKind {
    Fasta,
    Fastq,
}

/// One record. A FASTA record of 250 million bases is held whole: the
/// sketching then cuts it in windows, so its cost stays one copy.
#[derive(Clone, Debug, Default)]
pub struct SeqRecord {
    /// The header line without its first character (`>` or `@`).
    pub head: Vec<u8>,
    /// The sequence as written, white space removed.
    pub seq: Vec<u8>,
    /// FASTQ qualities; empty for FASTA.
    pub qual: Vec<u8>,
    /// Line where the record starts, from 1.
    pub line: u64,
}

pub struct SeqReader {
    path: String,
    r: Box<dyn BufRead>,
    kind: Option<SeqKind>,
    line_no: u64,
    buf: Vec<u8>,
    /// The FASTA header that ended the previous record.
    next_head: Option<(Vec<u8>, u64)>,
}

impl SeqReader {
    pub fn open(path: &str) -> Result<SeqReader, String> {
        check_input(path)?;
        let r = abre_leitura(path).map_err(|e| format!("{e}. Check the path and that you are allowed to read the file"))?;
        Ok(SeqReader { path: path.to_string(), r, kind: None, line_no: 0, buf: Vec::new(), next_head: None })
    }

    /// FASTA or FASTQ, known after the first record; `None` for an empty file.
    pub fn kind(&self) -> Option<SeqKind> {
        self.kind
    }

    /// The next line in `buf`, without its line break; `false` at the end.
    fn line(&mut self) -> Result<bool, String> {
        self.buf.clear();
        let n = self.r.read_until(b'\n', &mut self.buf).map_err(|e| {
            format!(
                "Error reading {} near line {}: {e}. If the file is compressed it may be damaged or incomplete: test it (`gzip -t`) or get it again",
                self.path,
                self.line_no + 1
            )
        })?;
        if n == 0 {
            return Ok(false);
        }
        self.line_no += 1;
        while matches!(self.buf.last(), Some(b'\n' | b'\r')) {
            self.buf.pop();
        }
        if self.line_no == 1 && self.buf.starts_with(&[0xEF, 0xBB, 0xBF]) {
            self.buf.drain(..3);
        }
        Ok(true)
    }

    fn not_sequences(&self) -> String {
        let hint = if self.buf.starts_with(&[0x1F, 0x8B]) {
            " The content is gzip-compressed: rename the file so that its name ends in .gz"
        } else if self.buf.starts_with(&[0x04, 0x22, 0x4D, 0x18]) {
            " The content is LZ4-compressed: rename the file so that its name ends in .lz4"
        } else {
            " Give sequences (FASTA or FASTQ, plain, .gz or .lz4), not a table of variants, intervals or alignments"
        };
        format!(
            "{}: this is not a FASTA nor a FASTQ file: line {} starts with \"{}\", and the first line of a FASTA file starts with '>' and of a FASTQ file with '@'.{hint}",
            self.path,
            self.line_no,
            sample_text(&self.buf)
        )
    }

    /// Reads the next record into `rec`; `false` at the end of the file.
    pub fn next(&mut self, rec: &mut SeqRecord) -> Result<bool, String> {
        rec.head.clear();
        rec.seq.clear();
        rec.qual.clear();
        if self.kind.is_none() {
            loop {
                if !self.line()? {
                    return Ok(false);
                }
                let t = trim(&self.buf);
                match t.first() {
                    None => continue,
                    Some(b'>') => self.kind = Some(SeqKind::Fasta),
                    Some(b'@') => self.kind = Some(SeqKind::Fastq),
                    Some(_) => return Err(self.not_sequences()),
                }
                self.next_head = Some((t[1..].to_vec(), self.line_no));
                break;
            }
        }
        match self.kind {
            Some(SeqKind::Fasta) => self.next_fasta(rec),
            _ => self.next_fastq(rec),
        }
    }

    fn next_fasta(&mut self, rec: &mut SeqRecord) -> Result<bool, String> {
        let Some((head, line)) = self.next_head.take() else {
            return Ok(false);
        };
        rec.head = head;
        rec.line = line;
        while self.line()? {
            let t = trim(&self.buf);
            match t.first() {
                None | Some(b';') => continue,
                Some(b'>') => {
                    self.next_head = Some((t[1..].to_vec(), self.line_no));
                    return Ok(true);
                }
                Some(_) => rec.seq.extend(t.iter().filter(|c| !c.is_ascii_whitespace())),
            }
        }
        Ok(true)
    }

    fn truncated(&self, start: u64, what: &str) -> String {
        format!(
            "{}: truncated FASTQ: the file ends in the middle of the record that starts at line {start} ({what}). A FASTQ record has 4 lines (@name, sequence, +, qualities). The file was probably cut during a copy or a download: get it again, or delete the incomplete last record",
            self.path
        )
    }

    fn next_fastq(&mut self, rec: &mut SeqRecord) -> Result<bool, String> {
        let start = match self.next_head.take() {
            Some((h, l)) => {
                rec.head = h;
                l
            }
            None => loop {
                if !self.line()? {
                    return Ok(false);
                }
                let t = trim(&self.buf);
                match t.first() {
                    None => continue,
                    Some(b'@') => {
                        rec.head = t[1..].to_vec();
                        break self.line_no;
                    }
                    Some(_) => {
                        return Err(format!(
                            "{}, line {}: expected a FASTQ record, which starts with '@', found \"{}\". The file is damaged or is not a 4-line FASTQ: check it with the validator, or get it again",
                            self.path,
                            self.line_no,
                            sample_text(t)
                        ))
                    }
                }
            },
        };
        rec.line = start;
        if !self.line()? {
            return Err(self.truncated(start, "its sequence, '+' and quality lines are missing"));
        }
        rec.seq.extend(self.buf.iter().filter(|c| !c.is_ascii_whitespace()));
        if !self.line()? {
            return Err(self.truncated(start, "its '+' line and its quality line are missing"));
        }
        if trim(&self.buf).first() != Some(&b'+') {
            return Err(format!(
                "{}, line {}: expected the '+' line of the FASTQ record that starts at line {start}, found \"{}\". A FASTQ record has exactly 4 lines (@name, sequence, +, qualities); if the sequences of this file are wrapped over several lines, rewrite it with one line per sequence (for instance `seqtk seq -l0`)",
                self.path,
                self.line_no,
                sample_text(&self.buf)
            ));
        }
        if !self.line()? {
            if rec.seq.is_empty() {
                return Ok(true);
            }
            return Err(self.truncated(start, "its quality line is missing"));
        }
        rec.qual.extend(self.buf.iter().filter(|c| !c.is_ascii_whitespace()));
        if rec.qual.len() != rec.seq.len() {
            // Short qualities at the very end of the file: a cut file.
            let at_end = self.r.fill_buf().map(|b| b.is_empty()).unwrap_or(false);
            if at_end && rec.qual.len() < rec.seq.len() {
                return Err(self.truncated(start, &format!("its quality line has {} characters for {} bases", rec.qual.len(), rec.seq.len())));
            }
            return Err(format!(
                "{}, line {}: the quality line has {} characters, but the sequence of the record that starts at line {start} has {} bases, and each base needs one quality character. The file is damaged, cut short or is not a 4-line FASTQ: check it with the validator, or get it again",
                self.path,
                self.line_no,
                rec.qual.len(),
                rec.seq.len()
            ));
        }
        Ok(true)
    }
}

/// FASTA or FASTQ, by the first record; `None` for a file without records.
pub fn sniff_kind(path: &str) -> Result<Option<SeqKind>, String> {
    let mut r = SeqReader::open(path)?;
    let mut rec = SeqRecord::default();
    r.next(&mut rec)?;
    Ok(r.kind())
}

// ---------------------------------------------------------------------
//  BATCHES: a reader thread decompresses and parses while the caller works
// ---------------------------------------------------------------------

type BatchMsg = Result<(Vec<SeqRecord>, SeqKind), String>;

fn reader_loop(paths: &[String], max_recs: usize, max_bases: usize, tx: &SyncSender<BatchMsg>) {
    for p in paths {
        let mut rd = match SeqReader::open(p) {
            Ok(r) => r,
            Err(e) => {
                let _ = tx.send(Err(e));
                return;
            }
        };
        loop {
            let mut recs = Vec::with_capacity(max_recs.min(4096));
            let mut bases = 0usize;
            let mut end = false;
            while recs.len() < max_recs && bases < max_bases {
                let mut rec = SeqRecord::default();
                match rd.next(&mut rec) {
                    Ok(true) => {
                        bases += rec.seq.len();
                        recs.push(rec);
                    }
                    Ok(false) => {
                        end = true;
                        break;
                    }
                    Err(e) => {
                        let _ = tx.send(Err(e));
                        return;
                    }
                }
            }
            // The receiver went away: the caller stopped on an error.
            if !recs.is_empty() && tx.send(Ok((recs, rd.kind().unwrap_or(SeqKind::Fasta)))).is_err() {
                return;
            }
            if end {
                break;
            }
        }
    }
}

/// Calls `f` with batches of records of the files, in order. A batch closes
/// at `max_recs` records or `max_bases` bases.
fn for_each_batch(
    paths: &[String],
    max_recs: usize,
    max_bases: usize,
    mut f: impl FnMut(&[SeqRecord], SeqKind) -> Result<(), String>,
) -> Result<(), String> {
    std::thread::scope(|s| {
        let (tx, rx) = sync_channel::<BatchMsg>(2);
        let owned = paths.to_vec();
        s.spawn(move || reader_loop(&owned, max_recs, max_bases, &tx));
        for msg in rx.iter() {
            let (recs, kind) = msg?;
            f(&recs, kind)?;
        }
        Ok(())
    })
}

/// The two files of paired reads, in batches of `n` pairs. The two files
/// must have the same number of reads, with the same names in the same order.
fn for_each_pair_batch(
    p1: &str,
    p2: &str,
    n: usize,
    mut f: impl FnMut(&[SeqRecord], &[SeqRecord], SeqKind) -> Result<(), String>,
) -> Result<(), String> {
    std::thread::scope(|s| {
        let (tx1, rx1) = sync_channel::<BatchMsg>(2);
        let (tx2, rx2) = sync_channel::<BatchMsg>(2);
        let a = vec![p1.to_string()];
        let b = vec![p2.to_string()];
        s.spawn(move || reader_loop(&a, n, usize::MAX, &tx1));
        s.spawn(move || reader_loop(&b, n, usize::MAX, &tx2));
        let mut done = 0u64;
        loop {
            let (m1, m2) = (rx1.recv(), rx2.recv());
            let (r1, r2, kind) = match (m1, m2) {
                (Err(_), Err(_)) => break,
                (Ok(x), Ok(y)) => {
                    let (r1, kind) = x?;
                    let (r2, _) = y?;
                    (r1, r2, kind)
                }
                (Ok(x), Err(_)) => {
                    let (r1, _) = x?;
                    return Err(unpaired(p1, p2, done, done + r1.len() as u64));
                }
                (Err(_), Ok(y)) => {
                    let (r2, _) = y?;
                    return Err(unpaired(p2, p1, done, done + r2.len() as u64));
                }
            };
            if r1.len() != r2.len() {
                let (longer, shorter, n_short) =
                    if r1.len() > r2.len() { (p1, p2, done + r2.len() as u64) } else { (p2, p1, done + r1.len() as u64) };
                return Err(unpaired(longer, shorter, n_short, n_short + 1));
            }
            for (i, (a, b)) in r1.iter().zip(r2.iter()).enumerate() {
                if read_id(&a.head) != read_id(&b.head) {
                    return Err(format!(
                        "The reads are not paired: read {} is \"{}\" in {p1} (line {}) but \"{}\" in {p2} (line {}). The two files of a pair must list the same reads in the same order: give the two files of the same run, or re-pair them (for instance `repair.sh` of BBTools)",
                        done + i as u64 + 1,
                        sample_text(read_id(&a.head)),
                        a.line,
                        sample_text(read_id(&b.head)),
                        b.line
                    ));
                }
            }
            f(&r1, &r2, kind)?;
            done += r1.len() as u64;
        }
        Ok(())
    })
}

fn unpaired(longer: &str, shorter: &str, n_short: u64, at_least: u64) -> String {
    format!(
        "The two files of the pair do not have the same number of reads: {shorter} ends after {} reads, but {longer} has at least {}. Give the two complete files of the same run (a file cut during the download has fewer reads)",
        thousands(n_short),
        thousands(at_least)
    )
}

// ---------------------------------------------------------------------
//  K-MER HASHES OF A SEQUENCE
// ---------------------------------------------------------------------

#[derive(Default)]
struct HashScratch {
    up: Vec<u8>,
    rc: Vec<u8>,
}

/// Calls `f(hash)` for every canonical k-mer of `seq` whose hash is
/// `<= max_hash`. Returns the number of valid k-mers (all A/C/G/T).
fn select_hashes(seq: &[u8], k: usize, max_hash: u64, sc: &mut HashScratch, mut f: impl FnMut(u64)) -> u64 {
    let n = seq.len();
    if n < k || k == 0 {
        return 0;
    }
    sc.up.clear();
    sc.up.extend(seq.iter().map(|&b| UPPER[usize::from(b)]));
    sc.rc.clear();
    sc.rc.extend(sc.up.iter().rev().map(|&b| COMP[usize::from(b)]));
    let mask = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let top = 2 * (k - 1);
    let (mut fw, mut rv, mut valid) = (0u64, 0u64, 0usize);
    let mut count = 0u64;
    for i in 0..n {
        let c = CODE[usize::from(sc.up[i])];
        if c > 3 {
            valid = 0;
            continue;
        }
        fw = ((fw << 2) | u64::from(c)) & mask;
        rv = (rv >> 2) | (u64::from(3 - c) << top);
        valid += 1;
        if valid >= k {
            count += 1;
            let s = i + 1 - k;
            // The 2-bit order is the text order (A < C < G < T).
            let h = if fw <= rv {
                murmur64(&sc.up[s..s + k], SEED)
            } else {
                let r = n - s - k;
                murmur64(&sc.rc[r..r + k], SEED)
            };
            if h <= max_hash {
                f(h);
            }
        }
    }
    count
}

/// Records longer than this are hashed in windows, in parallel.
const WINDOW: usize = 1 << 20;

/// The selected hashes (with repeats) of a batch of records, and the number
/// of valid k-mers. Long records are cut in windows that overlap by k - 1
/// bases, so that every k-mer is seen exactly once.
fn hash_records(recs: &[SeqRecord], k: usize, max_hash: u64) -> (Vec<u64>, u64) {
    let mut units: Vec<(usize, usize, usize)> = Vec::with_capacity(recs.len());
    for (r, rec) in recs.iter().enumerate() {
        let n = rec.seq.len();
        let mut s = 0;
        while s < n {
            let e = (s + WINDOW + k - 1).min(n);
            units.push((r, s, e));
            s += WINDOW;
        }
    }
    units
        .par_iter()
        .fold(
            || (Vec::new(), 0u64, HashScratch::default()),
            |(mut v, kc, mut sc), &(r, s, e)| {
                let c = select_hashes(&recs[r].seq[s..e], k, max_hash, &mut sc, |h| v.push(h));
                (v, kc + c, sc)
            },
        )
        .map(|(v, kc, _)| (v, kc))
        .reduce(
            || (Vec::new(), 0u64),
            |(mut a, ka), (b, kb)| {
                if a.len() < b.len() {
                    let mut b = b;
                    b.extend_from_slice(&a);
                    return (b, ka + kb);
                }
                a.extend_from_slice(&b);
                (a, ka + kb)
            },
        )
}

// ---------------------------------------------------------------------
//  1. REFERENCE SKETCHES AND THE SKETCH DATABASE
// ---------------------------------------------------------------------

/// One reference genome: one or more FASTA files whose sequences together
/// are the genome (chromosomes, plasmids).
#[derive(Clone, Debug, Default)]
pub struct GenomeInput {
    pub files: Vec<String>,
    /// The name in every report. Empty: the first header of the first file.
    pub name: String,
    /// A host decoy (human chr21 + chrM, ...): reported as host, never as a
    /// microbe, and its reads are the ones left out of the host-free output.
    pub host: bool,
}

#[derive(Clone, Debug)]
pub struct SketchOptions {
    /// k-mer length, 1 to 32 (31 by default; 21 is the other common choice).
    pub k: u8,
    /// One k-mer in `scaled` is kept, on average (1000 by default; 1 keeps
    /// them all).
    pub scaled: u64,
}

impl Default for SketchOptions {
    fn default() -> Self {
        SketchOptions { k: 31, scaled: 1000 }
    }
}

/// A sketch as stored in the database.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RefSketch {
    pub name: String,
    pub files: Vec<String>,
    pub host: bool,
    pub sequences: u64,
    pub bases: u64,
    /// Valid k-mers (all A/C/G/T) seen, with repeats.
    pub kmers: u64,
    /// The hashes kept, distinct, in increasing order.
    pub hashes: Vec<u64>,
}

/// A database read back (several files may be merged into one).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SketchDb {
    pub paths: Vec<String>,
    pub k: u8,
    pub scaled: u64,
    pub max_hash: u64,
    pub sketches: Vec<RefSketch>,
}

/// What `sketch_references` wrote, or what a database holds.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SketchSummary {
    pub path: String,
    pub k: u8,
    pub scaled: u64,
    pub max_hash: u64,
    /// The sketches without their hashes: the field `hashes` is left empty
    /// and `n_hashes` says how many there are.
    pub genomes: Vec<SketchInfo>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SketchInfo {
    pub name: String,
    pub files: Vec<String>,
    pub host: bool,
    pub sequences: u64,
    pub bases: u64,
    pub kmers: u64,
    pub n_hashes: u64,
}

impl SketchSummary {
    pub fn json(&self) -> String {
        let g: Vec<String> = self
            .genomes
            .iter()
            .map(|g| {
                format!(
                    "{{\"name\":{},\"files\":{},\"host\":{},\"sequences\":{},\"bases\":{},\"kmers\":{},\"hashes\":{}}}",
                    js(&g.name),
                    jlist(&g.files),
                    g.host,
                    g.sequences,
                    g.bases,
                    g.kmers,
                    g.n_hashes
                )
            })
            .collect();
        format!(
            "{{\"path\":{},\"k\":{},\"scaled\":{},\"max_hash\":{},\"genomes\":[{}]}}",
            js(&self.path),
            self.k,
            self.scaled,
            self.max_hash,
            g.join(",")
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!(
            "Sketch database {}: k = {}, scaled = {} (one k-mer in {} kept), {} genome(s)\n",
            self.path,
            self.k,
            self.scaled,
            self.scaled,
            self.genomes.len()
        );
        for g in &self.genomes {
            let host = if g.host { "  [host decoy]" } else { "" };
            o.push_str(&format!(
                "  {:>9} hashes  {:>14} bases  {:>5} sequence(s)  {}{host}\n",
                thousands(g.n_hashes),
                thousands(g.bases),
                g.sequences,
                g.name
            ));
            if g.kmers == 0 {
                o.push_str(&format!(
                    "    warning: no k-mer of {} A/C/G/T bases in a row: this genome can never be found\n",
                    self.k
                ));
            }
        }
        o
    }
}

fn info_of(s: &RefSketch) -> SketchInfo {
    SketchInfo {
        name: s.name.clone(),
        files: s.files.clone(),
        host: s.host,
        sequences: s.sequences,
        bases: s.bases,
        kmers: s.kmers,
        n_hashes: s.hashes.len() as u64,
    }
}

const DB_MAGIC: &[u8; 8] = b"MMGSKDB1";
const DB_END: &[u8; 8] = b"MMGSKEND";
const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

fn fnv1a(mut h: u64, data: &[u8]) -> u64 {
    for &b in data {
        h ^= u64::from(b);
        h = h.wrapping_mul(FNV_PRIME);
    }
    h
}

fn put_varint(out: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        out.push((v as u8) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

/// Writes the database as it is produced; the checksum follows the bytes.
struct DbWriter {
    w: std::io::BufWriter<File>,
    fnv: u64,
    path: String,
}

impl DbWriter {
    fn create(path: &str, k: u8, scaled: u64, max_hash: u64, n: u32) -> Result<DbWriter, String> {
        let f = File::create(path).map_err(|e| {
            format!("Error creating the sketch database {path}: {e}. Check that the output folder exists and that you can write to it")
        })?;
        let mut w = DbWriter { w: std::io::BufWriter::with_capacity(1 << 20, f), fnv: FNV_OFFSET, path: path.to_string() };
        let mut head = Vec::with_capacity(36);
        head.extend_from_slice(DB_MAGIC);
        head.push(k);
        head.push(0);
        head.extend_from_slice(&[0, 0]);
        head.extend_from_slice(&(SEED as u32).to_le_bytes());
        head.extend_from_slice(&scaled.to_le_bytes());
        head.extend_from_slice(&max_hash.to_le_bytes());
        head.extend_from_slice(&n.to_le_bytes());
        w.put(&head)?;
        Ok(w)
    }

    fn put(&mut self, bytes: &[u8]) -> Result<(), String> {
        self.fnv = fnv1a(self.fnv, bytes);
        self.w.write_all(bytes).map_err(|e| write_error(&self.path, e))
    }

    fn sketch(&mut self, s: &RefSketch) -> Result<(), String> {
        let mut b = Vec::with_capacity(64 + s.hashes.len() * 7);
        b.push(u8::from(s.host));
        for text in [s.name.as_bytes(), s.files.join("\n").as_bytes()] {
            b.extend_from_slice(&(text.len() as u32).to_le_bytes());
            b.extend_from_slice(text);
        }
        for v in [s.sequences, s.bases, s.kmers, s.hashes.len() as u64] {
            b.extend_from_slice(&v.to_le_bytes());
        }
        let mut prev = 0u64;
        for &h in &s.hashes {
            put_varint(&mut b, h - prev);
            prev = h;
        }
        self.put(&b)
    }

    fn finish(mut self) -> Result<(), String> {
        let fnv = self.fnv;
        self.w.write_all(DB_END).map_err(|e| write_error(&self.path, e))?;
        self.w.write_all(&fnv.to_le_bytes()).map_err(|e| write_error(&self.path, e))?;
        self.w.flush().map_err(|e| write_error(&self.path, e))
    }
}

fn sketch_files(files: &[String], k: usize, max_hash: u64) -> Result<(RefSketch, Option<String>), String> {
    let mut all: Vec<u64> = Vec::new();
    let mut s = RefSketch::default();
    let mut first: Option<String> = None;
    for_each_batch(files, 50_000, 16 << 20, |recs, _| {
        if first.is_none() {
            first = recs.first().map(|r| String::from_utf8_lossy(trim(&r.head)).to_string());
        }
        s.sequences += recs.len() as u64;
        s.bases += recs.iter().map(|r| r.seq.len() as u64).sum::<u64>();
        let (v, kc) = hash_records(recs, k, max_hash);
        s.kmers += kc;
        all.extend_from_slice(&v);
        // Keeps the memory near the number of distinct hashes even for a
        // genome full of repeats.
        if all.len() > (1 << 24) {
            multi_merge_sort_stable(&mut all);
            all.dedup();
        }
        Ok(())
    })?;
    multi_merge_sort_stable(&mut all);
    all.dedup();
    s.hashes = all;
    Ok((s, first))
}

/// Sketches each genome (FracMinHash, canonical k-mers, sourmash's hash) and
/// writes the sketch database to `output`. Each genome is read once, in a
/// stream; only its distinct kept hashes are held in memory.
pub fn sketch_references(genomes: &[GenomeInput], output: &str, opt: &SketchOptions) -> Result<SketchSummary, String> {
    if opt.k == 0 || opt.k > 32 {
        return Err(format!("Invalid k-mer length: k = {}. Use a k from 1 to 32 (31 is the usual one; 21 for more sensitivity to distant relatives)", opt.k));
    }
    if opt.scaled == 0 {
        return Err("Invalid scaled: 0. scaled is how many k-mers there are for each one kept: use 1000 for species identification, 10 for read classification, 1 to keep every k-mer".to_string());
    }
    if genomes.is_empty() {
        return Err("No genome to sketch: give at least one FASTA file".to_string());
    }
    let mut inputs: Vec<&str> = Vec::new();
    for (i, g) in genomes.iter().enumerate() {
        if g.files.is_empty() {
            return Err(format!("Genome {} ({}) has no file: give at least one FASTA file for each genome", i + 1, if g.name.is_empty() { "unnamed" } else { &g.name }));
        }
        for f in &g.files {
            check_input(f)?;
            inputs.push(f);
        }
    }
    not_an_input(output, &inputs)?;
    let max_hash = max_hash_for_scaled(opt.scaled);
    let mut cleanup = Cleanup(vec![output.to_string()]);
    let mut w = DbWriter::create(output, opt.k, opt.scaled, max_hash, genomes.len() as u32)?;
    let mut summary = SketchSummary { path: output.to_string(), k: opt.k, scaled: opt.scaled, max_hash, genomes: Vec::new() };
    for g in genomes {
        let (mut s, first) = sketch_files(&g.files, usize::from(opt.k), max_hash)?;
        s.name = if g.name.trim().is_empty() {
            first.unwrap_or_else(|| g.files[0].clone())
        } else {
            g.name.trim().to_string()
        };
        s.files = g.files.clone();
        s.host = g.host;
        w.sketch(&s)?;
        summary.genomes.push(info_of(&s));
    }
    w.finish()?;
    cleanup.disarm();
    Ok(summary)
}

fn get_varint(d: &[u8], pos: &mut usize) -> Option<u64> {
    let mut v = 0u64;
    let mut shift = 0;
    loop {
        let b = *d.get(*pos)?;
        *pos += 1;
        if shift >= 64 {
            return None;
        }
        v |= u64::from(b & 0x7F) << shift;
        if b < 0x80 {
            return Some(v);
        }
        shift += 7;
    }
}

/// Reads one sketch database file.
pub fn read_sketch_db(path: &str) -> Result<SketchDb, String> {
    check_input(path)?;
    let d = std::fs::read(path).map_err(|e| format!("Error reading the sketch database {path}: {e}"))?;
    let not_db = || {
        format!(
            "{path} is not an MMGEN sketch database (it does not start with \"MMGSKDB1\"). Build one with sketch_references from the FASTA files of the reference genomes"
        )
    };
    if d.len() < 8 || &d[..8] != DB_MAGIC {
        if d.len() >= 7 && &d[..7] == b"MMGSKDB" {
            return Err(format!(
                "{path} is a sketch database of a version this program does not know ({}): update MMGEN, or sketch the references again",
                String::from_utf8_lossy(&d[..8])
            ));
        }
        return Err(not_db());
    }
    let damaged = |what: &str| {
        format!(
            "The sketch database {path} is damaged or incomplete ({what}): it was probably cut during a copy. Copy it again, or sketch the references again"
        )
    };
    if d.len() < 36 + 16 || &d[d.len() - 16..d.len() - 8] != DB_END {
        return Err(damaged("its end mark is missing"));
    }
    let body = &d[..d.len() - 16];
    let stored = u64::from_le_bytes(d[d.len() - 8..].try_into().unwrap());
    if fnv1a(FNV_OFFSET, body) != stored {
        return Err(damaged("the checksum does not match the content"));
    }
    let k = body[8];
    let scaled = u64::from_le_bytes(body[16..24].try_into().unwrap());
    let max_hash = u64::from_le_bytes(body[24..32].try_into().unwrap());
    let n = u32::from_le_bytes(body[32..36].try_into().unwrap()) as usize;
    let seed = u32::from_le_bytes(body[12..16].try_into().unwrap());
    if body[9] != 0 || u64::from(seed) != SEED || k == 0 || k > 32 || scaled == 0 || max_hash != max_hash_for_scaled(scaled) {
        return Err(damaged("its header holds values this program never writes"));
    }
    let mut pos = 36usize;
    let mut sketches = Vec::with_capacity(n);
    let bad = || damaged("a sketch is cut short");
    for _ in 0..n {
        let flags = *body.get(pos).ok_or_else(bad)?;
        pos += 1;
        let mut texts: Vec<String> = Vec::new();
        for _ in 0..2 {
            let len = u32::from_le_bytes(body.get(pos..pos + 4).ok_or_else(bad)?.try_into().unwrap()) as usize;
            pos += 4;
            let t = body.get(pos..pos + len).ok_or_else(bad)?;
            texts.push(String::from_utf8_lossy(t).to_string());
            pos += len;
        }
        let mut nums = [0u64; 4];
        for v in nums.iter_mut() {
            *v = u64::from_le_bytes(body.get(pos..pos + 8).ok_or_else(bad)?.try_into().unwrap());
            pos += 8;
        }
        let count = nums[3] as usize;
        if count > body.len() {
            return Err(bad());
        }
        let mut hashes = Vec::with_capacity(count);
        let mut prev = 0u64;
        for i in 0..count {
            let delta = get_varint(body, &mut pos).ok_or_else(bad)?;
            if i > 0 && delta == 0 {
                return Err(damaged("the hashes are not in increasing order"));
            }
            prev = prev.checked_add(delta).ok_or_else(bad)?;
            hashes.push(prev);
        }
        let files = if texts[1].is_empty() { Vec::new() } else { texts[1].split('\n').map(str::to_string).collect() };
        sketches.push(RefSketch {
            name: std::mem::take(&mut texts[0]),
            files,
            host: flags & 1 == 1,
            sequences: nums[0],
            bases: nums[1],
            kmers: nums[2],
            hashes,
        });
    }
    if pos != body.len() {
        return Err(damaged("there are bytes after the last sketch"));
    }
    Ok(SketchDb { paths: vec![path.to_string()], k, scaled, max_hash, sketches })
}

/// What a sketch database holds (names, sizes), without the hashes.
pub fn describe_sketch_db(path: &str) -> Result<SketchSummary, String> {
    let db = read_sketch_db(path)?;
    Ok(SketchSummary {
        path: path.to_string(),
        k: db.k,
        scaled: db.scaled,
        max_hash: db.max_hash,
        genomes: db.sketches.iter().map(info_of).collect(),
    })
}

/// Reads and merges several databases (same k), downsampled to the larger
/// of their scaled values and `scaled`.
pub fn load_databases(paths: &[String], scaled: Option<u64>) -> Result<SketchDb, String> {
    if paths.is_empty() {
        return Err("No sketch database given: build one with sketch_references and give its path".to_string());
    }
    let mut dbs = Vec::with_capacity(paths.len());
    for p in paths {
        dbs.push(read_sketch_db(p)?);
    }
    let k = dbs[0].k;
    for (db, p) in dbs.iter().zip(paths) {
        if db.k != k {
            return Err(format!(
                "The sketch databases have different k: {} has k = {k} and {p} has k = {}. Hashes of different k never match: sketch all the references with the same k",
                paths[0], db.k
            ));
        }
    }
    let finest = dbs.iter().map(|d| d.scaled).max().unwrap_or(1);
    let target = match scaled {
        Some(0) => return Err("Invalid scaled: 0. Use 1 or more (1000 for species identification, 10 for read classification)".to_string()),
        Some(s) if s < finest => {
            let which = dbs.iter().zip(paths).find(|(d, _)| d.scaled == finest).map(|(_, p)| p.as_str()).unwrap_or("");
            return Err(format!(
                "scaled = {s} was asked, but the sketch database {which} was built with scaled = {finest}: a sketch can be made coarser (larger scaled) but not finer. Use scaled >= {finest}, or sketch the references again with scaled = {s}"
            ));
        }
        Some(s) => s,
        None => finest,
    };
    let max_hash = max_hash_for_scaled(target);
    let mut out = SketchDb { paths: paths.to_vec(), k, scaled: target, max_hash, sketches: Vec::new() };
    for db in dbs {
        for mut s in db.sketches {
            // Nested sketches: downsampling keeps a prefix of the sorted hashes.
            let keep = s.hashes.partition_point(|&h| h <= max_hash);
            s.hashes.truncate(keep);
            out.sketches.push(s);
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------
//  THE SAMPLE: hashes with abundances
// ---------------------------------------------------------------------

/// A sample sketched with abundances: `hashes` distinct and increasing,
/// `abunds[i]` the number of times `hashes[i]` was seen.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SampleSketch {
    pub hashes: Vec<u64>,
    pub abunds: Vec<u32>,
    pub sequences: u64,
    pub bases: u64,
    pub kmers: u64,
    /// Every file was FASTQ: the sample is reads, and the ANI can be corrected
    /// for coverage. An assembly (FASTA) is never corrected.
    pub reads: bool,
}

/// Sketches FASTA/FASTQ files (all together: the two files of a pair go in
/// the same list) keeping the abundance of each hash.
pub fn sketch_sample(paths: &[String], k: u8, max_hash: u64) -> Result<SampleSketch, String> {
    let mut counts: HashMap<u64, u32> = HashMap::new();
    let mut s = SampleSketch { reads: true, ..Default::default() };
    let mut any = false;
    for_each_batch(paths, 50_000, 16 << 20, |recs, kind| {
        any = true;
        if kind != SeqKind::Fastq {
            s.reads = false;
        }
        s.sequences += recs.len() as u64;
        s.bases += recs.iter().map(|r| r.seq.len() as u64).sum::<u64>();
        let (v, kc) = hash_records(recs, usize::from(k), max_hash);
        s.kmers += kc;
        for h in v {
            let c = counts.entry(h).or_insert(0);
            *c = c.saturating_add(1);
        }
        Ok(())
    })?;
    s.reads &= any;
    let mut pairs: Vec<(u64, u32)> = counts.into_iter().collect();
    multi_merge_sort_stable_by(&mut pairs, |a, b| a.0.cmp(&b.0));
    s.hashes = pairs.iter().map(|p| p.0).collect();
    s.abunds = pairs.iter().map(|p| p.1).collect();
    Ok(s)
}

// ---------------------------------------------------------------------
//  2. GATHER: the greedy minimum set cover of sourmash gather
// ---------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct GatherOptions {
    /// Downsample the databases (and the sample) to this scaled; `None`: the
    /// scaled of the databases.
    pub scaled: Option<u64>,
    /// A reference is reported only while it explains at least this many
    /// bases of what is left (unique hashes x scaled); sourmash's default.
    pub threshold_bp: u64,
    /// Hashes seen fewer times in the sample are dropped first (2 removes
    /// most sequencing errors of a deep sample); 1 keeps everything.
    pub min_abundance: u32,
    /// One line per match (TSV; the extension chooses the compression).
    pub output: Option<String>,
}

impl Default for GatherOptions {
    fn default() -> Self {
        GatherOptions { scaled: None, threshold_bp: 50_000, min_abundance: 1, output: None }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct GatherMatch {
    /// 1 for the first reference chosen.
    pub rank: usize,
    pub name: String,
    pub host: bool,
    /// Bases of the reference genome.
    pub genome_bp: u64,
    /// Hashes of the reference (after downsampling).
    pub reference_hashes: u64,
    /// Hashes of the reference found in the sample (in the original sample).
    pub shared_hashes: u64,
    /// Hashes assigned to this reference: shared with what the earlier
    /// matches had not taken.
    pub unique_hashes: u64,
    /// `shared_hashes` x scaled and `unique_hashes` x scaled: bases of the
    /// reference found in the whole sample, and in what was left of it when
    /// this match was chosen (sourmash's `intersect_bp` and
    /// `unique_intersect_bp`).
    pub intersect_bp: u64,
    pub unique_intersect_bp: u64,
    /// Fraction of the reference found in the sample: shared / reference.
    pub containment: f64,
    /// ANI estimated from the containment: containment^(1/k).
    pub ani: f64,
    /// The ANI with the containment corrected for the k-mers a shallow
    /// sequencing did not reach (reads only; `None` when the coverage cannot
    /// be estimated).
    pub ani_adjusted: Option<f64>,
    /// Mean k-mer coverage of this reference in the sample, estimated from the
    /// abundances of its shared hashes (Poisson truncated at zero).
    pub kmer_coverage: Option<f64>,
    /// Shared hashes / sample hashes.
    pub f_orig_query: f64,
    /// Unique hashes / sample hashes: the part of the sample this reference
    /// explains, without the part already explained.
    pub f_unique_to_query: f64,
    /// The same, weighted by abundance: the fraction of the sequenced k-mers
    /// (about the fraction of the DNA) that comes from this reference.
    pub f_unique_weighted: f64,
    /// `f_unique_weighted` among the non-host matches (they add up to 1).
    /// `None` for the host.
    pub relative_abundance: Option<f64>,
    /// Mean and median abundance of the unique hashes (~ k-mer depth).
    pub mean_abundance: f64,
    pub median_abundance: f64,
    /// What is left of the sample, in bases, after this match.
    pub remaining_bp: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct GatherResult {
    pub samples: Vec<String>,
    pub databases: Vec<String>,
    pub k: u8,
    pub scaled: u64,
    pub threshold_bp: u64,
    pub sequences: u64,
    pub bases: u64,
    /// The sample is reads (every file FASTQ).
    pub reads: bool,
    /// Distinct hashes of the sample, and the sum of their abundances, after
    /// `min_abundance`.
    pub sample_hashes: u64,
    pub sample_abundance: u64,
    /// Distinct hashes dropped by `min_abundance`.
    pub low_abundance_hashes: u64,
    pub references: usize,
    pub matches: Vec<GatherMatch>,
    /// Weighted fraction of the sample assigned to host references.
    pub host_weighted: f64,
    /// Fraction of the sample hashes no reference explained, plain and
    /// weighted by abundance (sequencing errors, unknown organisms, and what
    /// fell under the threshold).
    pub unassigned: f64,
    pub unassigned_weighted: f64,
}

/// Indices into `q` of the hashes of `r` that are also in `q` (both sorted).
fn intersect_indices(q: &[u64], r: &[u64]) -> Vec<u32> {
    let mut out = Vec::new();
    if r.len().saturating_mul(16) < q.len() {
        for &h in r {
            if let Ok(i) = q.binary_search(&h) {
                out.push(i as u32);
            }
        }
    } else {
        let (mut i, mut j) = (0, 0);
        while i < q.len() && j < r.len() {
            match q[i].cmp(&r[j]) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    out.push(i as u32);
                    i += 1;
                    j += 1;
                }
            }
        }
    }
    out
}

/// The k-mer coverage lambda of a Poisson truncated at zero whose mean is
/// `m` (m = lambda / (1 - e^-lambda)); `None` when m <= 1.
fn coverage_from_mean(m: f64) -> Option<f64> {
    if !m.is_finite() || m <= 1.0 {
        return None;
    }
    let f = |l: f64| l / -(-l).exp_m1();
    let (mut lo, mut hi) = (1e-12f64, m + 1.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < m {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

/// The k-mer coverage of a reference in the sample, from the abundances of
/// its shared hashes (mean `m`; `n1`, `n2` seen once and twice).
///
/// Deep coverage: from the mean (the k-mers missed are then few, and how
/// many does not matter much). Shallow coverage: the mean is pulled up by the
/// multi-copy k-mers (rDNA, mitochondria, plasmids, transposons), so the
/// coverage comes from the ratio of hashes seen twice to hashes seen once,
/// which a zero-truncated Poisson fixes at lambda / 2 and the few high-copy
/// k-mers barely touch. With too few of them to say: `None`.
fn kmer_coverage_estimate(m: f64, n1: u64, n2: u64) -> Option<f64> {
    match coverage_from_mean(m) {
        Some(l) if l >= 3.0 => Some(l),
        _ if n1 >= 20 && n2 >= 5 => Some(2.0 * n2 as f64 / n1 as f64),
        _ => None,
    }
}

fn median(mut v: Vec<u32>) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    multi_merge_sort_stable(&mut v);
    let n = v.len();
    if n % 2 == 1 {
        f64::from(v[n / 2])
    } else {
        (f64::from(v[n / 2 - 1]) + f64::from(v[n / 2])) / 2.0
    }
}

/// Gather on a sample already sketched (at `db.max_hash`).
///
/// Greedy: while some reference still shares at least `threshold_bp /
/// scaled` hashes with what is left of the sample, the one sharing the most
/// is reported and its hashes are removed from the sample. Ties go to the
/// reference sharing more with the whole sample, then to the smaller
/// reference, then to the earlier one in the database.
pub fn gather_sketch(sample: &SampleSketch, db: &SketchDb, opt: &GatherOptions) -> GatherResult {
    let q = &sample.hashes;
    let qa = &sample.abunds;
    let total_abund: u64 = qa.iter().map(|&a| u64::from(a)).sum();
    let per_ref: Vec<Vec<u32>> = db.sketches.par_iter().map(|s| intersect_indices(q, &s.hashes)).collect();
    let mut assigned = vec![false; q.len()];
    let mut used = vec![false; db.sketches.len()];
    let mut matches: Vec<GatherMatch> = Vec::new();
    let mut remaining = q.len() as u64;
    let k = f64::from(db.k);
    let scaled = db.scaled;
    loop {
        let mut best: Option<(usize, u64)> = None;
        for (r, idx) in per_ref.iter().enumerate() {
            if used[r] {
                continue;
            }
            let unique = idx.iter().filter(|&&i| !assigned[i as usize]).count() as u64;
            let better = match best {
                None => true,
                Some((b, bu)) => {
                    let (s_r, s_b) = (idx.len(), per_ref[b].len());
                    let (n_r, n_b) = (db.sketches[r].hashes.len(), db.sketches[b].hashes.len());
                    unique > bu || (unique == bu && (s_r > s_b || (s_r == s_b && n_r < n_b)))
                }
            };
            if better {
                best = Some((r, unique));
            }
        }
        let Some((r, unique)) = best else { break };
        if unique == 0 || unique.saturating_mul(scaled) < opt.threshold_bp {
            break;
        }
        used[r] = true;
        let s = &db.sketches[r];
        let mut ab: Vec<u32> = Vec::with_capacity(unique as usize);
        for &i in &per_ref[r] {
            if !assigned[i as usize] {
                assigned[i as usize] = true;
                ab.push(qa[i as usize]);
            }
        }
        remaining -= unique;
        let shared = per_ref[r].len() as u64;
        let unique_abund: u64 = ab.iter().map(|&a| u64::from(a)).sum();
        let containment = if s.hashes.is_empty() { 0.0 } else { shared as f64 / s.hashes.len() as f64 };
        let ani = if containment > 0.0 { containment.powf(1.0 / k) } else { 0.0 };
        let (kmer_coverage, ani_adjusted) = if sample.reads && opt.min_abundance <= 1 && shared >= 10 {
            let (mut sum, mut n1, mut n2) = (0u64, 0u64, 0u64);
            for &i in &per_ref[r] {
                let a = qa[i as usize];
                sum += u64::from(a);
                n1 += u64::from(a == 1);
                n2 += u64::from(a == 2);
            }
            match kmer_coverage_estimate(sum as f64 / shared as f64, n1, n2) {
                Some(l) => {
                    let seen = -(-l).exp_m1();
                    let c = (containment / seen).min(1.0);
                    (Some(l), Some(c.powf(1.0 / k)))
                }
                None => (None, None),
            }
        } else {
            (None, None)
        };
        matches.push(GatherMatch {
            rank: matches.len() + 1,
            name: s.name.clone(),
            host: s.host,
            genome_bp: s.bases,
            reference_hashes: s.hashes.len() as u64,
            shared_hashes: shared,
            unique_hashes: unique,
            intersect_bp: shared * scaled,
            unique_intersect_bp: unique * scaled,
            containment,
            ani,
            ani_adjusted,
            kmer_coverage,
            f_orig_query: if q.is_empty() { 0.0 } else { shared as f64 / q.len() as f64 },
            f_unique_to_query: if q.is_empty() { 0.0 } else { unique as f64 / q.len() as f64 },
            f_unique_weighted: if total_abund == 0 { 0.0 } else { unique_abund as f64 / total_abund as f64 },
            relative_abundance: None,
            mean_abundance: if unique == 0 { 0.0 } else { unique_abund as f64 / unique as f64 },
            median_abundance: median(ab),
            remaining_bp: remaining * scaled,
        });
    }
    let microbes: f64 = matches.iter().filter(|m| !m.host).map(|m| m.f_unique_weighted).sum();
    for m in matches.iter_mut() {
        if !m.host {
            m.relative_abundance = Some(if microbes > 0.0 { m.f_unique_weighted / microbes } else { 0.0 });
        }
    }
    let host_weighted: f64 = matches.iter().filter(|m| m.host).map(|m| m.f_unique_weighted).sum();
    let assigned_abund: u64 = assigned.iter().zip(qa).filter(|(a, _)| **a).map(|(_, &x)| u64::from(x)).sum();
    GatherResult {
        samples: Vec::new(),
        databases: db.paths.clone(),
        k: db.k,
        scaled,
        threshold_bp: opt.threshold_bp,
        sequences: sample.sequences,
        bases: sample.bases,
        reads: sample.reads,
        sample_hashes: q.len() as u64,
        sample_abundance: total_abund,
        low_abundance_hashes: 0,
        references: db.sketches.len(),
        matches,
        host_weighted,
        unassigned: if q.is_empty() { 0.0 } else { remaining as f64 / q.len() as f64 },
        unassigned_weighted: if total_abund == 0 { 0.0 } else { (total_abund - assigned_abund) as f64 / total_abund as f64 },
    }
}

/// Which reference genomes are in a sample, and how much of it each one
/// explains: sketches the sample (FASTA/FASTQ; the two files of a pair go
/// together in `samples`) and runs the greedy cover of `gather_sketch`.
pub fn gather(samples: &[String], databases: &[String], opt: &GatherOptions) -> Result<GatherResult, String> {
    if samples.is_empty() {
        return Err("No sample given: give the FASTQ (or FASTA) files of the sample".to_string());
    }
    for s in samples {
        check_input(s)?;
    }
    if opt.min_abundance == 0 {
        return Err("Invalid min_abundance: 0. Use 1 to keep every hash, 2 to drop those seen once (mostly sequencing errors)".to_string());
    }
    if let Some(o) = &opt.output {
        let inputs: Vec<&str> = samples.iter().chain(databases.iter()).map(String::as_str).collect();
        not_an_input(o, &inputs)?;
    }
    let db = load_databases(databases, opt.scaled)?;
    let mut sample = sketch_sample(samples, db.k, db.max_hash)?;
    let mut dropped = 0u64;
    if opt.min_abundance > 1 {
        let keep: Vec<bool> = sample.abunds.iter().map(|&a| a >= opt.min_abundance).collect();
        dropped = keep.iter().filter(|k| !**k).count() as u64;
        let mut it = keep.iter();
        sample.hashes.retain(|_| *it.next().unwrap());
        sample.abunds.retain(|&a| a >= opt.min_abundance);
    }
    let mut r = gather_sketch(&sample, &db, opt);
    r.samples = samples.to_vec();
    r.low_abundance_hashes = dropped;
    if let Some(o) = &opt.output {
        let mut cleanup = Cleanup(vec![o.clone()]);
        let mut w = create_output(o)?;
        w.write_all(r.tsv().as_bytes()).map_err(|e| write_error(o, e))?;
        w.termina().map_err(|e| write_error(o, e))?;
        cleanup.disarm();
    }
    Ok(r)
}

impl GatherMatch {
    fn json(&self) -> String {
        format!(
            "{{\"rank\":{},\"name\":{},\"host\":{},\"genome_bp\":{},\"reference_hashes\":{},\"shared_hashes\":{},\"unique_hashes\":{},\"intersect_bp\":{},\"unique_intersect_bp\":{},\"containment\":{},\"ani\":{},\"ani_adjusted\":{},\"kmer_coverage\":{},\"f_orig_query\":{},\"f_unique_to_query\":{},\"f_unique_weighted\":{},\"relative_abundance\":{},\"mean_abundance\":{},\"median_abundance\":{},\"remaining_bp\":{}}}",
            self.rank,
            js(&self.name),
            self.host,
            self.genome_bp,
            self.reference_hashes,
            self.shared_hashes,
            self.unique_hashes,
            self.intersect_bp,
            self.unique_intersect_bp,
            jf(self.containment),
            jf(self.ani),
            jo(self.ani_adjusted),
            jo(self.kmer_coverage),
            jf(self.f_orig_query),
            jf(self.f_unique_to_query),
            jf(self.f_unique_weighted),
            jo(self.relative_abundance),
            jf(self.mean_abundance),
            jf(self.median_abundance),
            self.remaining_bp
        )
    }
}

fn opt_num(x: Option<f64>, digits: usize) -> String {
    x.map_or_else(String::new, |v| format!("{v:.digits$}"))
}

fn clean_cell(s: &str) -> String {
    s.replace(['\t', '\n', '\r'], " ")
}

impl GatherResult {
    pub fn json(&self) -> String {
        let m: Vec<String> = self.matches.iter().map(GatherMatch::json).collect();
        format!(
            "{{\"samples\":{},\"databases\":{},\"k\":{},\"scaled\":{},\"threshold_bp\":{},\"sequences\":{},\"bases\":{},\"reads\":{},\"sample_hashes\":{},\"sample_abundance\":{},\"low_abundance_hashes\":{},\"references\":{},\"host_weighted\":{},\"unassigned\":{},\"unassigned_weighted\":{},\"matches\":[{}]}}",
            jlist(&self.samples),
            jlist(&self.databases),
            self.k,
            self.scaled,
            self.threshold_bp,
            self.sequences,
            self.bases,
            self.reads,
            self.sample_hashes,
            self.sample_abundance,
            self.low_abundance_hashes,
            self.references,
            jf(self.host_weighted),
            jf(self.unassigned),
            jf(self.unassigned_weighted),
            m.join(",")
        )
    }

    /// One line per match, with a header (what `output` receives).
    pub fn tsv(&self) -> String {
        let mut o = String::from("rank\tname\thost\tintersect_bp\tunique_intersect_bp\tf_unique_weighted\trelative_abundance\tf_unique_to_query\tf_orig_query\tcontainment\tani\tani_adjusted\tkmer_coverage\tmean_abundance\tmedian_abundance\tshared_hashes\tunique_hashes\treference_hashes\tgenome_bp\tremaining_bp\n");
        for m in &self.matches {
            o.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{:.3}\t{:.1}\t{}\t{}\t{}\t{}\t{}\n",
                m.rank,
                clean_cell(&m.name),
                m.host,
                m.intersect_bp,
                m.unique_intersect_bp,
                m.f_unique_weighted,
                opt_num(m.relative_abundance, 6),
                m.f_unique_to_query,
                m.f_orig_query,
                m.containment,
                m.ani,
                opt_num(m.ani_adjusted, 6),
                opt_num(m.kmer_coverage, 3),
                m.mean_abundance,
                m.median_abundance,
                m.shared_hashes,
                m.unique_hashes,
                m.reference_hashes,
                m.genome_bp,
                m.remaining_bp
            ));
        }
        o
    }

    pub fn resumo(&self) -> String {
        let mut o = format!(
            "Species identification (gather): k = {}, scaled = {}, {} reference(s) searched\n",
            self.k, self.scaled, self.references
        );
        o.push_str(&format!(
            "  sample: {} sequences, {} bases, {} distinct hashes, {} k-mer occurrences ({})\n",
            thousands(self.sequences),
            thousands(self.bases),
            thousands(self.sample_hashes),
            thousands(self.sample_abundance),
            if self.reads { "reads" } else { "assembly or FASTA" }
        ));
        if self.low_abundance_hashes > 0 {
            o.push_str(&format!("  {} low-abundance hashes dropped first\n", thousands(self.low_abundance_hashes)));
        }
        if self.matches.is_empty() {
            o.push_str(&format!(
                "  no reference explains {} bases or more of the sample: it holds none of the references, or too little of them (lower threshold_bp, or sequence more)\n",
                thousands(self.threshold_bp)
            ));
            return o;
        }
        o.push_str("  rank  abundance  of sample  containment  ANI     ANI adj.  name\n");
        for m in &self.matches {
            let ab = match m.relative_abundance {
                Some(x) => pct(x),
                None => "host".to_string(),
            };
            o.push_str(&format!(
                "  {:>4}  {:>9}  {:>9}  {:>11}  {:.4}  {:>8}  {}\n",
                m.rank,
                ab,
                pct(m.f_unique_weighted),
                pct(m.containment),
                m.ani,
                m.ani_adjusted.map_or_else(|| "-".to_string(), |a| format!("{a:.4}")),
                m.name
            ));
        }
        if self.host_weighted > 0.0 {
            o.push_str(&format!("  host: {} of the sample\n", pct(self.host_weighted)));
        }
        o.push_str(&format!(
            "  unassigned: {} of the sample k-mers ({} of its distinct hashes): sequencing errors, organisms not in the database, and what fell below the threshold\n",
            pct(self.unassigned_weighted),
            pct(self.unassigned)
        ));
        o.push_str("  abundance = share of the sequenced DNA among the microbes found (compare with the genomic-DNA percentages of a mock community, not with its cell numbers)\n");
        o
    }
}

// ---------------------------------------------------------------------
//  3. READ CLASSIFICATION
// ---------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct ClassifyOptions {
    /// Downsample the databases to this scaled; `None`: their own.
    pub scaled: Option<u64>,
    /// A read (pair) is classified only when its best reference shares at
    /// least this many sampled k-mers with it.
    pub min_hits: u32,
    /// ... and at least this fraction of the read's sampled k-mers (0 to 1;
    /// Kraken's "confidence"). 0 turns it off.
    pub min_fraction: f64,
    /// Use only the references with these names (empty: all). Typically the
    /// references `gather` found, so that reads are not spread over close
    /// relatives that are not in the sample.
    pub only: Vec<String>,
    /// The abundance table (TSV).
    pub output: Option<String>,
    /// One line per read or pair (TSV): name, status, reference(s), hits,
    /// sampled k-mers.
    pub assignments: Option<String>,
    /// The reads that are not host, in the input format (FASTQ or FASTA);
    /// `host_free2` takes the second reads of a pair.
    pub host_free: Option<String>,
    pub host_free2: Option<String>,
}

impl Default for ClassifyOptions {
    fn default() -> Self {
        ClassifyOptions {
            scaled: None,
            min_hits: 2,
            min_fraction: 0.0,
            only: Vec::new(),
            output: None,
            assignments: None,
            host_free: None,
            host_free2: None,
        }
    }
}

/// One line of the abundance table.
#[derive(Clone, Debug, PartialEq)]
pub struct ClassifyRow {
    pub name: String,
    pub host: bool,
    pub genome_bp: u64,
    /// Reads (pairs) whose best reference was this one alone.
    pub unique_reads: u64,
    /// Its share of the reads tied between it and other references, split in
    /// proportion to the unique reads of the tied references (evenly when
    /// none of them has any).
    pub shared_reads: f64,
    /// unique + shared.
    pub reads: f64,
    /// reads / classified non-host reads (`None` for the host).
    pub fraction_of_classified: Option<f64>,
    /// reads / all reads.
    pub fraction_of_all: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ClassifyResult {
    pub input: String,
    pub input2: Option<String>,
    pub databases: Vec<String>,
    pub k: u8,
    pub scaled: u64,
    pub references: usize,
    pub paired: bool,
    /// Reads, or pairs when paired.
    pub reads: u64,
    /// Non-host reads given to a reference (unique + ambiguous).
    pub classified: u64,
    pub unique: u64,
    /// Tied between several references.
    pub ambiguous: u64,
    pub host: u64,
    pub unclassified: u64,
    /// Reads (pairs) written to the host-free output.
    pub host_free_written: u64,
    /// Mean number of sampled k-mers in a read (pair).
    pub mean_sampled_kmers: f64,
    pub rows: Vec<ClassifyRow>,
    pub warnings: Vec<String>,
}

/// The sketch hashes of all references, with the set of references holding
/// each one. Sorted arrays and a directory of buckets: 12 bytes per hash.
struct ClassIndex {
    hashes: Vec<u64>,
    sets: Vec<u32>,
    /// `members[set]`: the references of a set. Sets 0..n are the single
    /// references themselves.
    members: Vec<Vec<u32>>,
    dir: Vec<u32>,
    buckets: u64,
    max_hash: u64,
}

impl ClassIndex {
    /// K-way merge of the sorted sketches (no sort needed).
    fn build(db: &SketchDb, keep: &[bool]) -> ClassIndex {
        let n = db.sketches.len();
        let total: usize = db.sketches.iter().zip(keep).filter(|(_, k)| **k).map(|(s, _)| s.hashes.len()).sum();
        let mut hashes = Vec::with_capacity(total);
        let mut sets = Vec::with_capacity(total);
        let mut members: Vec<Vec<u32>> = (0..n as u32).map(|r| vec![r]).collect();
        let mut interned: HashMap<Vec<u32>, u32> = HashMap::new();
        let mut pos = vec![0usize; n];
        let mut heap: BinaryHeap<Reverse<(u64, u32)>> = BinaryHeap::with_capacity(n);
        for (r, s) in db.sketches.iter().enumerate() {
            if keep[r] && !s.hashes.is_empty() {
                heap.push(Reverse((s.hashes[0], r as u32)));
            }
        }
        let mut cur: Vec<u32> = Vec::new();
        while let Some(Reverse((h, r))) = heap.pop() {
            cur.clear();
            let mut take = |r: u32, heap: &mut BinaryHeap<Reverse<(u64, u32)>>, cur: &mut Vec<u32>| {
                cur.push(r);
                let p = &mut pos[r as usize];
                *p += 1;
                if let Some(&next) = db.sketches[r as usize].hashes.get(*p) {
                    heap.push(Reverse((next, r)));
                }
            };
            take(r, &mut heap, &mut cur);
            while let Some(&Reverse((h2, r2))) = heap.peek() {
                if h2 != h {
                    break;
                }
                heap.pop();
                take(r2, &mut heap, &mut cur);
            }
            let id = if cur.len() == 1 {
                cur[0]
            } else {
                cur.sort_unstable();
                match interned.get(&cur) {
                    Some(&id) => id,
                    None => {
                        let id = members.len() as u32;
                        members.push(cur.clone());
                        interned.insert(cur.clone(), id);
                        id
                    }
                }
            };
            hashes.push(h);
            sets.push(id);
        }
        let buckets = ((hashes.len() / 4).max(1) as u64).next_power_of_two();
        let mut idx = ClassIndex { hashes, sets, members, dir: Vec::new(), buckets, max_hash: db.max_hash };
        let mut dir = Vec::with_capacity(buckets as usize + 1);
        let mut i = 0usize;
        for b in 0..=buckets {
            while i < idx.hashes.len() && idx.bucket(idx.hashes[i]) < b {
                i += 1;
            }
            dir.push(i as u32);
        }
        idx.dir = dir;
        idx
    }

    #[inline]
    fn bucket(&self, h: u64) -> u64 {
        ((u128::from(h) * u128::from(self.buckets)) / (u128::from(self.max_hash) + 1)) as u64
    }

    #[inline]
    fn lookup(&self, h: u64) -> Option<u32> {
        let b = self.bucket(h) as usize;
        let (s, e) = (self.dir[b] as usize, self.dir[b + 1] as usize);
        self.hashes[s..e].binary_search(&h).ok().map(|i| self.sets[s + i])
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum Outcome {
    Unclassified,
    Host(u32),
    Unique(u32),
    Tie(Vec<u32>),
}

struct ReadCall {
    outcome: Outcome,
    hits: u32,
    kmers: u32,
}

#[allow(clippy::too_many_arguments)]
fn classify_one(
    idx: &ClassIndex,
    k: usize,
    host: &[bool],
    mates: &[&[u8]],
    min_hits: u32,
    min_fraction: f64,
    sc: &mut HashScratch,
    counts: &mut Vec<(u32, u32)>,
) -> ReadCall {
    counts.clear();
    let mut kmers = 0u32;
    for m in mates {
        select_hashes(m, k, idx.max_hash, sc, |h| {
            kmers += 1;
            if let Some(set) = idx.lookup(h) {
                for &r in &idx.members[set as usize] {
                    match counts.iter_mut().find(|c| c.0 == r) {
                        Some(c) => c.1 += 1,
                        None => counts.push((r, 1)),
                    }
                }
            }
        });
    }
    let best = counts.iter().map(|c| c.1).max().unwrap_or(0);
    if best == 0 || best < min_hits || f64::from(best) < min_fraction * f64::from(kmers) {
        return ReadCall { outcome: Outcome::Unclassified, hits: best, kmers };
    }
    let mut tied: Vec<u32> = counts.iter().filter(|c| c.1 == best).map(|c| c.0).collect();
    tied.sort_unstable();
    let outcome = if let Some(&h) = tied.iter().find(|&&r| host[r as usize]) {
        // A read that matches the host as well as any microbe is host: it is
        // never released in the host-free output.
        Outcome::Host(h)
    } else if tied.len() == 1 {
        Outcome::Unique(tied[0])
    } else {
        Outcome::Tie(tied)
    };
    ReadCall { outcome, hits: best, kmers }
}

fn push_record(out: &mut Vec<u8>, rec: &SeqRecord, kind: SeqKind) {
    match kind {
        SeqKind::Fastq => {
            out.push(b'@');
            out.extend_from_slice(&rec.head);
            out.push(b'\n');
            out.extend_from_slice(&rec.seq);
            out.extend_from_slice(b"\n+\n");
            out.extend_from_slice(&rec.qual);
            out.push(b'\n');
        }
        SeqKind::Fasta => {
            out.push(b'>');
            out.extend_from_slice(&rec.head);
            out.push(b'\n');
            out.extend_from_slice(&rec.seq);
            out.push(b'\n');
        }
    }
}

/// Classifies each read (or pair) by the sketch hashes it shares with the
/// references: the reference sharing the most wins; reads tied between
/// references are "ambiguous" and their share is split in the table; a read
/// whose best match includes a host decoy is host. Single reads: `input`
/// alone; pairs: `input` and `input2`, judged together.
pub fn classify_reads(input: &str, input2: Option<&str>, databases: &[String], opt: &ClassifyOptions) -> Result<ClassifyResult, String> {
    check_input(input)?;
    if let Some(i2) = input2 {
        check_input(i2)?;
        if same_file(input, i2) {
            return Err(format!("The two files of the pair are the same file ({input}): give the file of the first reads and the file of the second reads"));
        }
    }
    if !(0.0..=1.0).contains(&opt.min_fraction) {
        return Err(format!("Invalid min_fraction: {}. It is a fraction of the read's sampled k-mers, from 0 (off) to 1", opt.min_fraction));
    }
    match (input2.is_some(), opt.host_free.is_some(), opt.host_free2.is_some()) {
        (true, true, false) => {
            return Err("Paired reads with host_free but without host_free2: give host_free2 for the second reads of the pairs".to_string())
        }
        (true, false, true) => return Err("host_free2 was given without host_free: give host_free for the first reads of the pairs".to_string()),
        (false, _, true) => return Err("host_free2 is for the second file of paired reads, and only one read file was given: leave host_free2 out".to_string()),
        _ => {}
    }
    let outputs: Vec<&String> = [&opt.output, &opt.assignments, &opt.host_free, &opt.host_free2].into_iter().flatten().collect();
    let mut inputs: Vec<&str> = vec![input];
    inputs.extend(input2);
    inputs.extend(databases.iter().map(String::as_str));
    for (i, o) in outputs.iter().enumerate() {
        not_an_input(o, &inputs)?;
        if outputs[..i].iter().any(|p| same_file(p, o)) {
            return Err(format!("Two outputs are the same file ({o}): give each output its own name"));
        }
    }

    let mut db = load_databases(databases, opt.scaled)?;
    let mut keep = vec![true; db.sketches.len()];
    if !opt.only.is_empty() {
        for name in &opt.only {
            if !db.sketches.iter().any(|s| &s.name == name) {
                let known: Vec<&str> = db.sketches.iter().map(|s| s.name.as_str()).take(12).collect();
                return Err(format!(
                    "\"{name}\" (in only) is not a reference of the sketch database. Use the names as the database has them; some of them: {}",
                    known.join(" | ")
                ));
            }
        }
        for (i, s) in db.sketches.iter().enumerate() {
            keep[i] = opt.only.iter().any(|n| n == &s.name);
        }
    }
    let idx = ClassIndex::build(&db, &keep);
    // The hashes now live in the index: let the sketches go.
    for s in db.sketches.iter_mut() {
        s.hashes = Vec::new();
    }
    let host: Vec<bool> = db.sketches.iter().map(|s| s.host).collect();
    let n_ref = db.sketches.len();
    let k = usize::from(db.k);

    let mut cleanup = Cleanup(outputs.iter().map(|s| s.to_string()).collect());
    let mut w_assign = match &opt.assignments {
        Some(p) => {
            let mut w = create_output(p)?;
            w.write_all(b"read\tstatus\treference\thits\tsampled_kmers\n").map_err(|e| write_error(p, e))?;
            Some(w)
        }
        None => None,
    };
    let mut w_free = match &opt.host_free {
        Some(p) => Some(create_output(p)?),
        None => None,
    };
    let mut w_free2 = match &opt.host_free2 {
        Some(p) => Some(create_output(p)?),
        None => None,
    };

    let mut unique_counts = vec![0u64; n_ref];
    let mut host_counts = vec![0u64; n_ref];
    let mut ties: HashMap<Vec<u32>, u64> = HashMap::new();
    let mut r = ClassifyResult {
        input: input.to_string(),
        input2: input2.map(str::to_string),
        databases: databases.to_vec(),
        k: db.k,
        scaled: db.scaled,
        references: keep.iter().filter(|x| **x).count(),
        paired: input2.is_some(),
        reads: 0,
        classified: 0,
        unique: 0,
        ambiguous: 0,
        host: 0,
        unclassified: 0,
        host_free_written: 0,
        mean_sampled_kmers: 0.0,
        rows: Vec::new(),
        warnings: Vec::new(),
    };
    let mut sampled_total = 0u64;
    let (min_hits, min_fraction) = (opt.min_hits, opt.min_fraction);
    let names: Vec<String> = db.sketches.iter().map(|s| clean_cell(&s.name)).collect();

    let mut handle = |r1: &[SeqRecord], r2: Option<&[SeqRecord]>, kind: SeqKind| -> Result<(), String> {
        let calls: Vec<ReadCall> = (0..r1.len())
            .into_par_iter()
            .map_init(
                || (HashScratch::default(), Vec::new()),
                |(sc, counts), i| {
                    let mates: Vec<&[u8]> = match r2 {
                        Some(r2) => vec![&r1[i].seq, &r2[i].seq],
                        None => vec![&r1[i].seq],
                    };
                    classify_one(&idx, k, &host, &mates, min_hits, min_fraction, sc, counts)
                },
            )
            .collect();
        let mut lines: Vec<u8> = Vec::new();
        let mut free1: Vec<u8> = Vec::new();
        let mut free2: Vec<u8> = Vec::new();
        for (i, c) in calls.iter().enumerate() {
            r.reads += 1;
            sampled_total += u64::from(c.kmers);
            let (status, refs): (&str, String) = match &c.outcome {
                Outcome::Unclassified => {
                    r.unclassified += 1;
                    ("unclassified", "-".to_string())
                }
                Outcome::Host(h) => {
                    r.host += 1;
                    host_counts[*h as usize] += 1;
                    ("host", names[*h as usize].clone())
                }
                Outcome::Unique(u) => {
                    r.unique += 1;
                    unique_counts[*u as usize] += 1;
                    ("classified", names[*u as usize].clone())
                }
                Outcome::Tie(t) => {
                    r.ambiguous += 1;
                    *ties.entry(t.clone()).or_insert(0) += 1;
                    ("ambiguous", t.iter().map(|&x| names[x as usize].as_str()).collect::<Vec<_>>().join(";"))
                }
            };
            if w_assign.is_some() {
                lines.extend_from_slice(read_id(&r1[i].head));
                lines.extend_from_slice(format!("\t{status}\t{refs}\t{}\t{}\n", c.hits, c.kmers).as_bytes());
            }
            if w_free.is_some() && !matches!(c.outcome, Outcome::Host(_)) {
                r.host_free_written += 1;
                push_record(&mut free1, &r1[i], kind);
                if let Some(r2) = r2 {
                    push_record(&mut free2, &r2[i], kind);
                }
            }
        }
        if let (Some(w), Some(p)) = (w_assign.as_mut(), opt.assignments.as_ref()) {
            w.write_all(&lines).map_err(|e| write_error(p, e))?;
        }
        if let (Some(w), Some(p)) = (w_free.as_mut(), opt.host_free.as_ref()) {
            w.write_all(&free1).map_err(|e| write_error(p, e))?;
        }
        if let (Some(w), Some(p)) = (w_free2.as_mut(), opt.host_free2.as_ref()) {
            w.write_all(&free2).map_err(|e| write_error(p, e))?;
        }
        Ok(())
    };
    match input2 {
        Some(i2) => for_each_pair_batch(input, i2, 16_384, |a, b, kind| handle(a, Some(b), kind))?,
        None => for_each_batch(&[input.to_string()], 16_384, 64 << 20, |a, kind| handle(a, None, kind))?,
    }
    for (w, p) in [(w_assign, &opt.assignments), (w_free, &opt.host_free), (w_free2, &opt.host_free2)] {
        if let (Some(w), Some(p)) = (w, p) {
            w.termina().map_err(|e| write_error(p, e))?;
        }
    }

    // The table: unique reads, plus the share of the tied ones.
    let mut shared = vec![0f64; n_ref];
    let mut tie_list: Vec<(&Vec<u32>, &u64)> = ties.iter().collect();
    tie_list.sort();
    for (set, &n) in tie_list {
        let base: u64 = set.iter().map(|&x| unique_counts[x as usize]).sum();
        for &x in set {
            shared[x as usize] += if base > 0 {
                n as f64 * unique_counts[x as usize] as f64 / base as f64
            } else {
                n as f64 / set.len() as f64
            };
        }
    }
    r.classified = r.unique + r.ambiguous;
    for i in 0..n_ref {
        let reads = if host[i] { host_counts[i] as f64 } else { unique_counts[i] as f64 + shared[i] };
        if reads <= 0.0 {
            continue;
        }
        r.rows.push(ClassifyRow {
            name: db.sketches[i].name.clone(),
            host: host[i],
            genome_bp: db.sketches[i].bases,
            unique_reads: if host[i] { host_counts[i] } else { unique_counts[i] },
            shared_reads: if host[i] { 0.0 } else { shared[i] },
            reads,
            fraction_of_classified: if host[i] {
                None
            } else if r.classified > 0 {
                Some(reads / r.classified as f64)
            } else {
                Some(0.0)
            },
            fraction_of_all: if r.reads > 0 { reads / r.reads as f64 } else { 0.0 },
        });
    }
    r.rows.sort_by(|a, b| a.host.cmp(&b.host).then(b.reads.total_cmp(&a.reads)).then(a.name.cmp(&b.name)));
    r.mean_sampled_kmers = if r.reads > 0 { sampled_total as f64 / r.reads as f64 } else { 0.0 };
    if r.reads > 0 && r.mean_sampled_kmers < 4.0 {
        r.warnings.push(format!(
            "A read{} holds only {:.2} sampled k-mers on average at scaled = {}: most reads cannot be classified. Sketch the references with a smaller scaled for read classification (10 gives ~12 per 150-base read; 1 keeps every k-mer)",
            if r.paired { " pair" } else { "" },
            r.mean_sampled_kmers,
            r.scaled
        ));
    }
    if let Some(o) = &opt.output {
        let mut w = create_output(o)?;
        w.write_all(r.tsv().as_bytes()).map_err(|e| write_error(o, e))?;
        w.termina().map_err(|e| write_error(o, e))?;
    }
    cleanup.disarm();
    Ok(r)
}

impl ClassifyResult {
    pub fn json(&self) -> String {
        let rows: Vec<String> = self
            .rows
            .iter()
            .map(|x| {
                format!(
                    "{{\"name\":{},\"host\":{},\"genome_bp\":{},\"unique_reads\":{},\"shared_reads\":{},\"reads\":{},\"fraction_of_classified\":{},\"fraction_of_all\":{}}}",
                    js(&x.name),
                    x.host,
                    x.genome_bp,
                    x.unique_reads,
                    jf(x.shared_reads),
                    jf(x.reads),
                    jo(x.fraction_of_classified),
                    jf(x.fraction_of_all)
                )
            })
            .collect();
        format!(
            "{{\"input\":{},\"input2\":{},\"databases\":{},\"k\":{},\"scaled\":{},\"references\":{},\"paired\":{},\"reads\":{},\"classified\":{},\"unique\":{},\"ambiguous\":{},\"host\":{},\"unclassified\":{},\"host_free_written\":{},\"mean_sampled_kmers\":{},\"rows\":[{}],\"warnings\":{}}}",
            js(&self.input),
            self.input2.as_deref().map_or_else(|| "null".to_string(), js),
            jlist(&self.databases),
            self.k,
            self.scaled,
            self.references,
            self.paired,
            self.reads,
            self.classified,
            self.unique,
            self.ambiguous,
            self.host,
            self.unclassified,
            self.host_free_written,
            jf(self.mean_sampled_kmers),
            rows.join(","),
            jlist(&self.warnings)
        )
    }

    /// The abundance table (what `output` receives).
    pub fn tsv(&self) -> String {
        let mut o = String::from("name\thost\treads\tunique_reads\tshared_reads\tfraction_of_classified\tfraction_of_all\tgenome_bp\n");
        for x in &self.rows {
            o.push_str(&format!(
                "{}\t{}\t{:.2}\t{}\t{:.2}\t{}\t{:.6}\t{}\n",
                clean_cell(&x.name),
                x.host,
                x.reads,
                x.unique_reads,
                x.shared_reads,
                opt_num(x.fraction_of_classified, 6),
                x.fraction_of_all,
                x.genome_bp
            ));
        }
        o
    }

    pub fn resumo(&self) -> String {
        let unit = if self.paired { "pairs" } else { "reads" };
        let all = self.reads.max(1) as f64;
        let mut o = format!(
            "Read classification: k = {}, scaled = {}, {} reference(s), {:.1} sampled k-mers per {}\n",
            self.k,
            self.scaled,
            self.references,
            self.mean_sampled_kmers,
            if self.paired { "pair" } else { "read" }
        );
        o.push_str(&format!("  {unit}: {}\n", thousands(self.reads)));
        o.push_str(&format!(
            "  classified:   {} ({}), of which {} tied between references\n",
            thousands(self.classified),
            pct(self.classified as f64 / all),
            thousands(self.ambiguous)
        ));
        o.push_str(&format!("  host:         {} ({})\n", thousands(self.host), pct(self.host as f64 / all)));
        o.push_str(&format!("  unclassified: {} ({})\n", thousands(self.unclassified), pct(self.unclassified as f64 / all)));
        if !self.rows.is_empty() {
            o.push_str(&format!("  {:>12}  {:>11}  {:>9}  name\n", unit, "of classif.", "of all"));
            for x in &self.rows {
                o.push_str(&format!(
                    "  {:>12}  {:>11}  {:>9}  {}{}\n",
                    thousands(x.reads.round() as u64),
                    x.fraction_of_classified.map_or_else(|| "host".to_string(), pct),
                    pct(x.fraction_of_all),
                    x.name,
                    if x.host { "  [host]" } else { "" }
                ));
            }
        }
        if self.host_free_written > 0 {
            o.push_str(&format!("  host-free {unit} written: {}\n", thousands(self.host_free_written)));
        }
        for w in &self.warnings {
            o.push_str(&format!("  warning: {w}\n"));
        }
        o
    }
}

// ---------------------------------------------------------------------
//  4. RESISTANCE GENES: the AMRFinderPlus nucleotide catalog
// ---------------------------------------------------------------------

/// Seeds: 21-mers (odd, so no k-mer is its own reverse complement).
const SEED_K: usize = 21;
/// Diagonals searched on each side of the seeds (room for indels).
const BAND: i64 = 24;
/// Seeds of one gene on one strand whose diagonals differ by more than this
/// are two loci.
const DIAG_GAP: i64 = 2 * BAND;
/// blastn's scores.
const MATCH: i32 = 2;
const MISMATCH: i32 = -3;
const GAP_OPEN: i32 = 5;
const GAP_EXTEND: i32 = 2;
/// Reads mode: a read counts for a gene when its alignment has at least this
/// many columns and this identity.
const READ_MIN_COLUMNS: u32 = 30;
const READ_MIN_IDENTITY: f64 = 0.8;

/// One gene of the catalog.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AmrGene {
    /// Protein accession (first field of the AMR_CDS.fa header).
    pub protein: String,
    pub nucleotide: String,
    /// Gene or allele symbol (fifth field): `blaKPC-3`, `mecA`.
    pub symbol: String,
    /// Family (sixth field): `blaKPC`, `mecA`.
    pub family: String,
    pub product: String,
    /// From ReferenceGeneCatalog.txt (empty without it).
    pub scope: String,
    pub element_type: String,
    pub subtype: String,
    pub class: String,
    pub subclass: String,
    /// The coding sequence, upper case (IUPAC codes kept).
    pub seq: Vec<u8>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AmrMode {
    /// FASTQ input means reads; FASTA means an assembly.
    Auto,
    Assembly,
    Reads,
}

#[derive(Clone, Debug)]
pub struct AmrOptions {
    pub min_identity: f64,
    pub min_coverage: f64,
    /// A gene cut by the end of a contig is reported when it covers at least
    /// this much of the gene (status `partial_contig_end`); 1 turns it off.
    pub min_partial_coverage: f64,
    /// Report the "plus" genes of the catalog too (virulence, stress, and the
    /// AMR genes outside the core set); false keeps AMRFinderPlus's default,
    /// core only.
    pub include_plus: bool,
    pub mode: AmrMode,
    /// Reads mode: a gene position counts as covered with this many reads.
    pub min_depth: u32,
    /// The table of hits (TSV).
    pub output: Option<String>,
}

impl Default for AmrOptions {
    fn default() -> Self {
        AmrOptions {
            min_identity: 0.9,
            min_coverage: 0.9,
            min_partial_coverage: 0.5,
            include_plus: true,
            mode: AmrMode::Auto,
            min_depth: 1,
            output: None,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct AmrHit {
    /// The symbol of the best allele (`blaKPC-3`); with identity below 1 the
    /// sample carries a close relative of it.
    pub gene: String,
    pub family: String,
    pub product: String,
    pub scope: String,
    pub element_type: String,
    pub subtype: String,
    pub class: String,
    pub subclass: String,
    pub identity: f64,
    pub coverage: f64,
    /// Assembly: contig name (first word of its header), 1-based inclusive
    /// coordinates on its forward strand, and the strand of the gene. Reads:
    /// empty, 0, 0, '.'.
    pub contig: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub gene_length: u32,
    /// Alignment columns (assembly) or covered positions (reads).
    pub aligned_length: u32,
    pub mismatches: u32,
    pub gap_bases: u32,
    /// "exact" (identity 1, whole gene), "aligned" (passes the thresholds) or
    /// "partial_contig_end" (cut by the end of a contig).
    pub status: String,
    /// Reads mode: mean depth over the covered positions.
    pub depth: Option<f64>,
    /// Reads mode: the reads of this gene disagree at some position (20% or
    /// more of the depth): more than one allele or copy of it is present.
    /// When another catalog allele explains them it is reported too (at most
    /// two alleles per group of near-identical genes); assemble the reads to
    /// see every copy. Always false on an assembly.
    pub mixed: bool,
    pub protein: String,
    pub nucleotide: String,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AmrResult {
    pub inputs: Vec<String>,
    pub catalog: String,
    pub catalog_table: Option<String>,
    /// "assembly" or "reads".
    pub mode: String,
    pub sequences: u64,
    pub bases: u64,
    /// Catalog genes searched (after `include_plus`).
    pub catalog_genes: usize,
    /// Catalog genes without a line in the table (their class is unknown).
    pub catalog_unmatched: usize,
    /// Alignments computed.
    pub alignments: u64,
    pub hits: Vec<AmrHit>,
}

fn catalog_error(path: &str, line: u64, head: &[u8]) -> String {
    format!(
        "{path}, line {line}: this header is not in the AMRFinderPlus AMR_CDS.fa form `protein|nucleotide|part|parts|gene|family|product`: \">{}\". Give the AMR_CDS.fa file of the AMRFinderPlus database (https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/)",
        sample_text(head)
    )
}

struct TableRow {
    allele: String,
    family: String,
    node: String,
    product: String,
    scope: String,
    element_type: String,
    subtype: String,
    class: String,
    subclass: String,
}

/// The lines of ReferenceGeneCatalog.txt, and their indices by protein
/// accession (RefSeq and GenBank).
type CatalogTable = (Vec<TableRow>, HashMap<String, Vec<usize>>);

fn read_catalog_table(path: &str) -> Result<CatalogTable, String> {
    check_input(path)?;
    let mut r = abre_leitura(path)?;
    let mut line = String::new();
    let mut rows = Vec::new();
    let mut by_acc: HashMap<String, Vec<usize>> = HashMap::new();
    let mut cols: Option<Vec<usize>> = None;
    let mut node_col: Option<usize> = None;
    let mut acc_cols = (0usize, 0usize);
    let mut n = 0u64;
    loop {
        line.clear();
        let got = r.read_line(&mut line).map_err(|e| format!("Error reading {path} near line {}: {e}", n + 1))?;
        if got == 0 {
            break;
        }
        n += 1;
        let t = line.trim_end_matches(['\n', '\r']);
        if t.is_empty() {
            continue;
        }
        let f: Vec<&str> = t.split('\t').collect();
        let Some(c) = &cols else {
            let find = |name: &str| f.iter().position(|x| x.trim_start_matches('#') == name);
            let needed = ["allele", "gene_family", "product_name", "scope", "type", "subtype", "class", "subclass"];
            let mut v = Vec::new();
            for name in needed.iter().chain(["refseq_protein_accession", "genbank_protein_accession"].iter()) {
                match find(name) {
                    Some(i) => v.push(i),
                    None => {
                        return Err(format!(
                            "{path}: the column \"{name}\" is not in its first line, so this is not the ReferenceGeneCatalog.txt of AMRFinderPlus. Give the ReferenceGeneCatalog.txt of the same database version as AMR_CDS.fa"
                        ))
                    }
                }
            }
            acc_cols = (v[8], v[9]);
            node_col = find("hierarchy_node");
            cols = Some(v);
            continue;
        };
        let get = |i: usize| f.get(c[i]).map_or(String::new(), |s| s.trim().to_string());
        let row = TableRow {
            allele: get(0),
            family: get(1),
            node: node_col.and_then(|i| f.get(i)).map_or(String::new(), |s| s.trim().to_string()),
            product: get(2),
            scope: get(3),
            element_type: get(4),
            subtype: get(5),
            class: get(6),
            subclass: get(7),
        };
        for a in [acc_cols.0, acc_cols.1] {
            if let Some(acc) = f.get(a).map(|s| s.trim()).filter(|s| !s.is_empty()) {
                by_acc.entry(acc.to_string()).or_default().push(rows.len());
            }
        }
        rows.push(row);
    }
    if cols.is_none() {
        return Err(format!("{path} is empty: give the ReferenceGeneCatalog.txt of the AMRFinderPlus database"));
    }
    Ok((rows, by_acc))
}

/// Reads the nucleotide catalog (AMR_CDS.fa) and, if given, the table of
/// classes (ReferenceGeneCatalog.txt). Returns the genes and how many of
/// them had no line in the table.
pub fn load_amr_catalog(fasta: &str, table: Option<&str>, include_plus: bool) -> Result<(Vec<AmrGene>, usize), String> {
    let mut rd = SeqReader::open(fasta)?;
    let mut rec = SeqRecord::default();
    let mut genes: Vec<AmrGene> = Vec::new();
    // (protein, part, parts) of each gene, to join the two parts of a fusion
    let mut parts: Vec<(u32, u32)> = Vec::new();
    while rd.next(&mut rec)? {
        if rd.kind() == Some(SeqKind::Fastq) {
            return Err(format!("{fasta} is a FASTQ file, not the AMR_CDS.fa catalog (a FASTA file): give the catalog of the AMRFinderPlus database"));
        }
        let first = rec.head.split(|c| c.is_ascii_whitespace()).next().unwrap_or(&[]);
        let f: Vec<&[u8]> = first.split(|&c| c == b'|').collect();
        if f.len() < 7 || f[0].is_empty() || f[4].is_empty() {
            return Err(catalog_error(fasta, rec.line, &rec.head));
        }
        let text = |b: &[u8]| String::from_utf8_lossy(b).to_string();
        let part: u32 = text(f[2]).parse().unwrap_or(1);
        let of: u32 = text(f[3]).parse().unwrap_or(1);
        let seq: Vec<u8> = rec.seq.iter().map(|b| b.to_ascii_uppercase()).collect();
        if seq.is_empty() {
            continue;
        }
        let g = AmrGene {
            protein: text(f[0]),
            nucleotide: text(f[1]),
            symbol: text(f[4]),
            family: text(if f[5].is_empty() { f[4] } else { f[5] }),
            product: text(&f[6..].join(&b'|')).replace('_', " "),
            seq,
            ..AmrGene::default()
        };
        // A bifunctional gene is listed twice, part 1 of 2 and part 2 of 2,
        // with the same sequence: one locus, both names.
        if part > 1 && of > 1 {
            if let (Some(prev), Some(pp)) = (genes.last_mut(), parts.last()) {
                if prev.protein == g.protein && prev.seq == g.seq && pp.1 == of {
                    prev.symbol = format!("{}/{}", prev.symbol, g.symbol);
                    prev.family = format!("{}/{}", prev.family, g.family);
                    continue;
                }
            }
        }
        genes.push(g);
        parts.push((part, of));
    }
    if genes.is_empty() {
        return Err(format!("{fasta} holds no gene: give the AMR_CDS.fa file of the AMRFinderPlus database"));
    }
    let mut unmatched = 0usize;
    if let Some(tp) = table {
        let (rows, by_acc) = read_catalog_table(tp)?;
        for g in genes.iter_mut() {
            let cands = by_acc.get(&g.protein).map(Vec::as_slice).unwrap_or(&[]);
            let pick = cands
                .iter()
                .find(|&&i| {
                    let r = &rows[i];
                    g.symbol == r.allele || g.symbol == r.family || g.symbol == r.node
                })
                .or_else(|| cands.iter().find(|&&i| rows[i].allele.is_empty()))
                .or_else(|| cands.first());
            match pick {
                Some(&i) => {
                    let r = &rows[i];
                    if !r.product.is_empty() {
                        g.product = r.product.clone();
                    }
                    g.scope = r.scope.clone();
                    g.element_type = r.element_type.clone();
                    g.subtype = r.subtype.clone();
                    g.class = r.class.clone();
                    g.subclass = r.subclass.clone();
                }
                None => unmatched += 1,
            }
        }
        if !include_plus {
            genes.retain(|g| g.scope != "plus");
        }
    }
    Ok((genes, unmatched))
}

/// Calls `f(start, forward, reverse)` for each k-mer of A/C/G/T (k <= 32).
fn roll_kmers(seq: &[u8], k: usize, mut f: impl FnMut(usize, u64, u64)) {
    let mask = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let top = 2 * (k - 1);
    let (mut fw, mut rv, mut valid) = (0u64, 0u64, 0usize);
    for (i, &b) in seq.iter().enumerate() {
        let c = CODE[usize::from(b)];
        if c > 3 {
            valid = 0;
            continue;
        }
        fw = ((fw << 2) | u64::from(c)) & mask;
        rv = (rv >> 2) | (u64::from(3 - c) << top);
        valid += 1;
        if valid >= k {
            f(i + 1 - k, fw, rv);
        }
    }
}

/// The canonical 21-mers of the catalog, each with its genes and positions.
struct SeedIndex {
    keys: Vec<u64>,
    /// postings of `keys[i]`: `starts[i]..starts[i + 1]`
    starts: Vec<u32>,
    genes: Vec<u32>,
    /// position << 1 | 1 when the gene's forward k-mer is the canonical one
    pos: Vec<u32>,
    dir: Vec<u32>,
}

const DIR_BITS: u32 = 20;
const DIR_SHIFT: u32 = 2 * SEED_K as u32 - DIR_BITS;

impl SeedIndex {
    fn build(genes: &[AmrGene]) -> SeedIndex {
        let mut all: Vec<(u64, u32, u32)> = genes
            .par_iter()
            .enumerate()
            .flat_map_iter(|(g, gene)| {
                let mut v = Vec::with_capacity(gene.seq.len());
                roll_kmers(&gene.seq, SEED_K, |p, fw, rv| v.push((fw.min(rv), g as u32, ((p as u32) << 1) | u32::from(fw < rv))));
                v
            })
            .collect();
        multi_merge_sort_stable_by(&mut all, |a, b| a.0.cmp(&b.0));
        let mut keys = Vec::new();
        let mut starts = Vec::new();
        let mut gs = Vec::with_capacity(all.len());
        let mut ps = Vec::with_capacity(all.len());
        for (i, &(key, g, p)) in all.iter().enumerate() {
            if i == 0 || all[i - 1].0 != key {
                keys.push(key);
                starts.push(i as u32);
            }
            gs.push(g);
            ps.push(p);
        }
        starts.push(all.len() as u32);
        drop(all);
        let mut dir = Vec::with_capacity((1usize << DIR_BITS) + 1);
        let mut i = 0usize;
        for b in 0..=(1u64 << DIR_BITS) {
            while i < keys.len() && (keys[i] >> DIR_SHIFT) < b {
                i += 1;
            }
            dir.push(i as u32);
        }
        SeedIndex { keys, starts, genes: gs, pos: ps, dir }
    }

    #[inline]
    fn find(&self, key: u64) -> std::ops::Range<usize> {
        let b = (key >> DIR_SHIFT) as usize;
        let (s, e) = (self.dir[b] as usize, self.dir[b + 1] as usize);
        match self.keys[s..e].binary_search(&key) {
            Ok(i) => self.starts[s + i] as usize..self.starts[s + i + 1] as usize,
            Err(_) => 0..0,
        }
    }

    /// Calls `f(gene, minus, diagonal)` for every seed of `seq`. On the plus
    /// strand the diagonal is `seq position - gene position`; on the minus
    /// strand it is the anti-diagonal `seq position + gene position` of the
    /// gene's first base of the k-mer and the sequence's last base of it.
    fn seeds(&self, seq: &[u8], offset: usize, mut f: impl FnMut(u32, bool, i64)) {
        roll_kmers(seq, SEED_K, |c, fw, rv| {
            let range = self.find(fw.min(rv));
            if range.is_empty() {
                return;
            }
            let seq_fw = fw < rv;
            let c = (c + offset) as i64;
            for t in range {
                let p = self.pos[t];
                let gpos = i64::from(p >> 1);
                let minus = seq_fw != (p & 1 == 1);
                let d = if minus { c + (SEED_K as i64 - 1) + gpos } else { c - gpos };
                f(self.genes[t], minus, d);
            }
        });
    }
}

#[derive(Clone, Debug, Default)]
struct Aln {
    score: i32,
    /// gene positions [g0, g1) and window positions [c0, c1), 0-based
    g0: usize,
    g1: usize,
    c0: usize,
    c1: usize,
    matches: u32,
    mismatches: u32,
    ins: u32,
    del: u32,
    /// Reads mode: (gene position, symbol) with 0..3 the base of the read,
    /// 4 a deletion, 5 an insertion after that position.
    ops: Vec<(u32, u8)>,
}

impl Aln {
    fn columns(&self) -> u32 {
        self.matches + self.mismatches + self.ins + self.del
    }
}

const NEG: i32 = i32::MIN / 4;
const STOP: u8 = 0;
const DIAG: u8 = 1;
const FROM_E: u8 = 2;
const FROM_F: u8 = 3;
const E_EXT: u8 = 4;
const F_EXT: u8 = 8;

/// Alignment (affine gaps) of `gene` (rows) against `win` (columns),
/// restricted to the diagonals `j - i` in `[dlo, dhi]`.
///
/// `overlap = None`: local (an assembly: the gene is somewhere in the contig
/// region). `overlap = Some((gene_starts, gene_ends))`: a read, aligned end
/// to end except where it runs past the gene -- before the first row when
/// `gene_starts` (the rows begin at the gene's first base) and after the
/// last row when `gene_ends`. A read is never clipped in the middle of the
/// gene, so an indel shows in the pileup instead of being cut away.
fn align_banded(gene: &[u8], win: &[u8], dlo: i64, dhi: i64, want_ops: bool, overlap: Option<(bool, bool)>) -> Option<Aln> {
    let l = gene.len();
    let w = win.len() as i64;
    if l == 0 || w == 0 || dhi < dlo {
        return None;
    }
    let b = (dhi - dlo + 1) as usize;
    let local = overlap.is_none();
    let (gene_starts, gene_ends) = overlap.unwrap_or((true, true));
    // the floor of every cell: 0 starts a local alignment anywhere
    let floor = if local { 0 } else { NEG };
    // a cell outside the read: its start (j = 0) is free in both modes
    let outside = |j: i64| if local || j == 0 { 0 } else { NEG };
    let mut dirs = vec![STOP; (l + 1) * b];
    let mut hp: Vec<i32> = (0..=b)
        .map(|t| {
            let j = dlo + t as i64;
            if local || j == 0 || (j >= 1 && j <= w && t < b && gene_starts) {
                0
            } else {
                NEG
            }
        })
        .collect();
    let mut fp = vec![NEG; b + 1];
    let (mut hc, mut fc) = (vec![floor; b + 1], vec![NEG; b + 1]);
    let mut best = (if local { 0i32 } else { NEG }, 0usize, 0usize);
    for i in 1..=l {
        let gm = GENE_MASK[usize::from(gene[i - 1])];
        let (mut e, mut h_left) = (NEG, outside(i as i64 + dlo - 1));
        for t in 0..b {
            let j = i as i64 + dlo + t as i64;
            if j < 1 || j > w {
                hc[t] = outside(j);
                fc[t] = NEG;
                e = NEG;
                h_left = outside(j);
                continue;
            }
            let s = if gm & SAMPLE_MASK[usize::from(win[j as usize - 1])] != 0 { MATCH } else { MISMATCH };
            let dg = hp[t] + s;
            let (e_open, e_ext) = (h_left - GAP_OPEN - GAP_EXTEND, e - GAP_EXTEND);
            let (ev, eflag) = if e_ext > e_open { (e_ext, E_EXT) } else { (e_open, 0) };
            let (f_open, f_ext) = (hp[t + 1] - GAP_OPEN - GAP_EXTEND, fp[t + 1] - GAP_EXTEND);
            let (fv, fflag) = if f_ext > f_open { (f_ext, F_EXT) } else { (f_open, 0) };
            let (mut h, mut d) = (floor, STOP);
            if dg > h {
                h = dg;
                d = DIAG;
            }
            if ev > h {
                h = ev;
                d = FROM_E;
            }
            if fv > h {
                h = fv;
                d = FROM_F;
            }
            dirs[i * b + t] = d | eflag | fflag;
            hc[t] = h;
            fc[t] = fv;
            e = ev;
            h_left = h;
            // where an alignment may end: anywhere (local), or at the end of
            // the read, or at the end of the gene
            let may_end = local || j == w || (i == l && gene_ends);
            if may_end && h > best.0 {
                best = (h, i, t);
            }
        }
        hc[b] = floor;
        fc[b] = NEG;
        std::mem::swap(&mut hp, &mut hc);
        std::mem::swap(&mut fp, &mut fc);
    }
    if (local && best.0 <= 0) || (!local && (best.0 <= NEG / 2 || best.1 == 0)) {
        return None;
    }
    let (mut i, mut t) = (best.1, best.2);
    let col = |i: usize, t: usize| (i as i64 + dlo + t as i64).max(0) as usize;
    let mut a = Aln { score: best.0, g1: i, c1: col(i, t), ..Aln::default() };
    let mut state = 0u8; // 0 = H, 1 = E (gap in the gene), 2 = F (gap in the window)
    loop {
        // Row 0, column 0 and the edges of the band hold no alignment.
        if i == 0 || (state == 1 && t == 0) || (state == 2 && t + 1 >= b) {
            break;
        }
        let j = col(i, t);
        if j == 0 {
            break;
        }
        let d = dirs[i * b + t];
        match state {
            0 => match d & 3 {
                DIAG => {
                    let base = win[j - 1];
                    if GENE_MASK[usize::from(gene[i - 1])] & SAMPLE_MASK[usize::from(base)] != 0 {
                        a.matches += 1;
                    } else {
                        a.mismatches += 1;
                    }
                    if want_ops {
                        let c = CODE[usize::from(base)];
                        if c < 4 {
                            a.ops.push(((i - 1) as u32, c));
                        }
                    }
                    i -= 1;
                }
                FROM_E => state = 1,
                FROM_F => state = 2,
                _ => break,
            },
            1 => {
                a.ins += 1;
                if want_ops && a.ops.last() != Some(&((i - 1) as u32, 5)) {
                    a.ops.push(((i - 1) as u32, 5));
                }
                if d & E_EXT == 0 {
                    state = 0;
                }
                t -= 1;
            }
            _ => {
                a.del += 1;
                if want_ops {
                    a.ops.push(((i - 1) as u32, 4));
                }
                if d & F_EXT == 0 {
                    state = 0;
                }
                i -= 1;
                t += 1;
            }
        }
    }
    a.g0 = i;
    a.c0 = col(i, t);
    Some(a)
}

/// A gene found on a contig, before the best allele of each locus is chosen.
#[derive(Clone, Debug)]
struct RawHit {
    gene: u32,
    minus: bool,
    start: usize,
    end: usize,
    aln: Aln,
    identity: f64,
    coverage: f64,
    #[allow(dead_code)]
    partial: bool,
}

impl RawHit {
    fn exact(&self, len: usize) -> bool {
        self.aln.matches == self.aln.columns() && self.aln.g0 == 0 && self.aln.g1 == len
    }
}

/// Best first: exact matches, then identity x coverage, identity, score, and
/// the earlier gene of the catalog.
fn hit_order(a: &RawHit, b: &RawHit, genes: &[AmrGene]) -> std::cmp::Ordering {
    let (la, lb) = (genes[a.gene as usize].seq.len(), genes[b.gene as usize].seq.len());
    b.exact(lb)
        .cmp(&a.exact(la))
        .then((b.identity * b.coverage).total_cmp(&(a.identity * a.coverage)))
        .then(b.identity.total_cmp(&a.identity))
        .then(b.aln.score.cmp(&a.aln.score))
        .then(a.gene.cmp(&b.gene))
}

/// Aligns one gene to the region of the contig its seeds point to.
fn align_on_contig(gene: &[u8], contig: &[u8], minus: bool, dmin: i64, dmax: i64) -> Option<(Aln, usize, usize)> {
    let n = contig.len() as i64;
    let l = gene.len() as i64;
    if !minus {
        let lo = (dmin - BAND).max(0);
        let hi = (dmax + BAND + l).min(n);
        if hi <= lo {
            return None;
        }
        let a = align_banded(gene, &contig[lo as usize..hi as usize], dmin - BAND - lo, dmax + BAND - lo, false, None)?;
        let (s, e) = (lo as usize + a.c0, lo as usize + a.c1);
        Some((a, s, e))
    } else {
        // gene base i sits at contig position A - i (A = anti-diagonal)
        let lo = (dmin - BAND - (l - 1)).max(0);
        let hi = (dmax + BAND + 1).min(n);
        if hi <= lo {
            return None;
        }
        let rc: Vec<u8> = contig[lo as usize..hi as usize].iter().rev().map(|&c| COMP[usize::from(UPPER[usize::from(c)])]).collect();
        let a = align_banded(gene, &rc, hi - 1 - (dmax + BAND), hi - 1 - (dmin - BAND), false, None)?;
        let (s, e) = (hi as usize - a.c1, hi as usize - a.c0);
        Some((a, s, e))
    }
}

/// The seeds of one contig, grouped into clusters (gene, strand, diagonals).
fn contig_clusters(seq: &[u8], idx: &SeedIndex) -> Vec<(u32, bool, i64, i64, u32)> {
    let n = seq.len();
    let mut windows = Vec::new();
    let mut s = 0;
    while s < n {
        windows.push((s, (s + WINDOW + SEED_K - 1).min(n)));
        s += WINDOW;
    }
    let maps: Vec<HashMap<(u32, bool, i64), u32>> = windows
        .par_iter()
        .map(|&(s, e)| {
            let mut m: HashMap<(u32, bool, i64), u32> = HashMap::new();
            idx.seeds(&seq[s..e], s, |g, minus, d| *m.entry((g, minus, d)).or_insert(0) += 1);
            m
        })
        .collect();
    let mut all: Vec<(u32, bool, i64, u32)> = Vec::new();
    if maps.len() == 1 {
        all.extend(maps.into_iter().next().unwrap().into_iter().map(|((g, m, d), c)| (g, m, d, c)));
    } else {
        let mut merged: HashMap<(u32, bool, i64), u32> = HashMap::new();
        for m in maps {
            for (key, c) in m {
                *merged.entry(key).or_insert(0) += c;
            }
        }
        all.extend(merged.into_iter().map(|((g, m, d), c)| (g, m, d, c)));
    }
    multi_merge_sort_stable_by(&mut all, |a, b| (a.0, a.1, a.2).cmp(&(b.0, b.1, b.2)));
    let mut clusters: Vec<(u32, bool, i64, i64, u32)> = Vec::new();
    for (g, m, d, c) in all {
        match clusters.last_mut() {
            Some(last) if last.0 == g && last.1 == m && d - last.3 <= DIAG_GAP => {
                last.3 = d;
                last.4 += c;
            }
            _ => clusters.push((g, m, d, d, c)),
        }
    }
    clusters
}

fn scan_contig(seq: &[u8], idx: &SeedIndex, genes: &[AmrGene], opt: &AmrOptions) -> (Vec<RawHit>, u64) {
    let clusters: Vec<_> = contig_clusters(seq, idx).into_iter().filter(|c| c.4 >= 2).collect();
    let n = seq.len();
    let hits: Vec<RawHit> = clusters
        .par_iter()
        .filter_map(|&(g, minus, dmin, dmax, _)| {
            let gene = &genes[g as usize].seq;
            let l = gene.len();
            let (aln, start, end) = align_on_contig(gene, seq, minus, dmin, dmax)?;
            let columns = aln.columns();
            if columns == 0 {
                return None;
            }
            let identity = f64::from(aln.matches) / f64::from(columns);
            let coverage = (aln.g1 - aln.g0) as f64 / l as f64;
            let (room_start, room_end) = if minus { (n - end, start) } else { (start, n - end) };
            let partial = (aln.g0 > 0 && room_start < aln.g0) || (aln.g1 < l && room_end < l - aln.g1);
            let passes = identity >= opt.min_identity
                && (coverage >= opt.min_coverage || (partial && coverage >= opt.min_partial_coverage));
            passes.then_some(RawHit { gene: g, minus, start, end, aln, identity, coverage, partial })
        })
        .collect();
    let mut hits = hits;
    let aligned = clusters.len() as u64;
    hits.sort_by(|a, b| hit_order(a, b, genes));
    // The best allele of each locus: a hit overlapping a better one by half
    // of the shorter of the two is the same locus.
    let mut kept: Vec<RawHit> = Vec::new();
    for h in hits {
        let same_locus = kept.iter().any(|k| {
            let ov = h.end.min(k.end) as i64 - h.start.max(k.start) as i64;
            ov > 0 && 2 * ov as usize >= (h.end - h.start).min(k.end - k.start)
        });
        if !same_locus {
            kept.push(h);
        }
    }
    kept.sort_by_key(|h| (h.start, h.end));
    (kept, aligned)
}

fn hit_of(g: &AmrGene, identity: f64, coverage: f64) -> AmrHit {
    AmrHit {
        gene: g.symbol.clone(),
        family: g.family.clone(),
        product: g.product.clone(),
        scope: g.scope.clone(),
        element_type: g.element_type.clone(),
        subtype: g.subtype.clone(),
        class: g.class.clone(),
        subclass: g.subclass.clone(),
        identity,
        coverage,
        contig: String::new(),
        start: 0,
        end: 0,
        strand: '.',
        gene_length: g.seq.len() as u32,
        aligned_length: 0,
        mismatches: 0,
        gap_bases: 0,
        status: String::new(),
        depth: None,
        mixed: false,
        protein: g.protein.clone(),
        nucleotide: g.nucleotide.clone(),
    }
}

/// One read aligned to one gene: the gene, the pileup operations, the
/// differences (mismatches + inserted + deleted bases) in the aligned part,
/// and the alignment columns.
type ReadAln = (u32, Vec<(u32, u8)>, u32, u32);

/// Reads mode: the alignments of one read to the catalog genes it seeds well.
fn read_alignments(seq: &[u8], idx: &SeedIndex, genes: &[AmrGene], rc: &mut Vec<u8>) -> (Vec<ReadAln>, u32) {
    let mut seeds: Vec<(u32, bool, i64)> = Vec::new();
    idx.seeds(seq, 0, |g, minus, d| seeds.push((g, minus, d)));
    if seeds.is_empty() {
        return (Vec::new(), 0);
    }
    seeds.sort_unstable();
    // per (gene, strand): seeds, lowest and highest diagonal
    let mut groups: Vec<(u32, bool, i64, i64, u32)> = Vec::new();
    for (g, m, d) in seeds {
        match groups.last_mut() {
            Some(x) if x.0 == g && x.1 == m => {
                x.3 = d;
                x.4 += 1;
            }
            _ => groups.push((g, m, d, d, 1)),
        }
    }
    let best = groups.iter().map(|x| x.4).max().unwrap_or(0);
    let floor = (best / 2).max(2);
    let rl = seq.len() as i64;
    let mut have_rc = false;
    let mut out = Vec::new();
    let mut tried = 0u32;
    for &(g, minus, dmin, dmax, n) in &groups {
        // A read cannot hold a gene region with diagonals more than a read
        // apart: seeds that far apart are repeats inside the gene.
        if n < floor || dmax - dmin > rl {
            continue;
        }
        let gene = &genes[g as usize].seq;
        let l = gene.len() as i64;
        let (win, dlo, dhi): (&[u8], i64, i64) = if minus {
            if !have_rc {
                rc.clear();
                rc.extend(seq.iter().rev().map(|&c| COMP[usize::from(UPPER[usize::from(c)])]));
                have_rc = true;
            }
            (rc.as_slice(), rl - 1 - dmax - BAND, rl - 1 - dmin + BAND)
        } else {
            (seq, dmin - BAND, dmax + BAND)
        };
        // only the gene rows the band can reach inside the read
        let gs = (-dhi).max(0).min(l);
        let ge = (rl - dlo).min(l).max(gs);
        if ge <= gs {
            continue;
        }
        tried += 1;
        let Some(mut a) = align_banded(&gene[gs as usize..ge as usize], win, dlo + gs, dhi + gs, true, Some((gs == 0, ge == l))) else { continue };
        let columns = a.columns();
        if columns < READ_MIN_COLUMNS || f64::from(a.matches) < READ_MIN_IDENTITY * f64::from(columns) {
            continue;
        }
        for op in a.ops.iter_mut() {
            op.0 += gs as u32;
        }
        let diffs = a.mismatches + a.ins + a.del;
        out.push((g, a.ops, diffs, columns));
    }
    (out, tried)
}

/// Per gene position: A, C, G, T, deletion, insertion after it -- from every
/// read aligned to the gene [0..6], from the reads that match this gene best
/// [6..12], and from the reads within 5% of it [12..18].
type Pile = [u32; 18];

/// Reads mode: a position where the reads within 5% of the gene disagree --
/// the second most common base has at least `MIXED_READS` reads and this
/// share of a depth of at least `MIXED_DEPTH` -- says that two alleles (or
/// copies) are present.
const MIXED_SHARE: f64 = 0.25;
const MIXED_DEPTH: u32 = 10;
const MIXED_READS: u32 = 3;

/// What the pile of reads says of one gene, and whether its close reads are
/// mixed. The consensus of a position comes from the reads that match this
/// gene best when they hold at least a quarter of its depth (so the reads of
/// another locus, a close allele, do not vote), and from all its reads
/// otherwise (so a relative seen through another gene's reads shows its
/// differences).
fn evaluate_pileup(g: &AmrGene, pile: &[Pile], min_depth: u32) -> (AmrHit, bool) {
    let mixed = pile.iter().any(|p| {
        let depth: u32 = p[12..17].iter().sum();
        if depth < MIXED_DEPTH {
            return false;
        }
        let mut c: [u32; 5] = p[12..17].try_into().unwrap();
        c.sort_unstable();
        c[3] >= MIXED_READS && f64::from(c[3]) >= MIXED_SHARE * f64::from(depth)
    });
    let mut covered = 0u32;
    let mut matches = 0u32;
    let mut dels = 0u32;
    let mut ins = 0u32;
    let mut depth_sum = 0u64;
    for (i, p) in pile.iter().enumerate() {
        let depth: u32 = p[..5].iter().sum();
        if depth == 0 || depth < min_depth {
            continue;
        }
        covered += 1;
        depth_sum += u64::from(depth);
        let best_depth: u32 = p[6..11].iter().sum();
        let o = if best_depth >= min_depth && 4 * best_depth >= depth { 6 } else { 0 };
        let used: u32 = p[o..o + 5].iter().sum();
        let (best, _) = p[o..o + 5].iter().enumerate().fold((0usize, 0u32), |acc, (s, &c)| if c > acc.1 { (s, c) } else { acc });
        if best == 4 {
            dels += 1;
        } else if GENE_MASK[usize::from(g.seq[i])] & (1u8 << best) != 0 {
            matches += 1;
        }
        if 2 * p[o + 5] > used {
            ins += 1;
        }
    }
    let l = g.seq.len() as f64;
    let columns = covered + ins;
    let mut h = hit_of(g, if columns == 0 { 0.0 } else { f64::from(matches) / f64::from(columns) }, f64::from(covered) / l);
    h.aligned_length = covered;
    h.mismatches = covered - matches - dels;
    h.gap_bases = dels + ins;
    h.depth = Some(if covered == 0 { 0.0 } else { depth_sum as f64 / f64::from(covered) });
    (h, mixed)
}

/// Finds the genes of the AMRFinderPlus nucleotide catalog in an assembly
/// (contigs, FASTA) or in reads (FASTQ; the two files of a pair may both be
/// listed), with identity and coverage of each gene, keeping the best allele
/// of each locus (assembly: by position; reads: by shared reads).
pub fn find_resistance_genes(inputs: &[String], catalog: &str, catalog_table: Option<&str>, opt: &AmrOptions) -> Result<AmrResult, String> {
    if inputs.is_empty() {
        return Err("No input given: give the contigs (FASTA) or the reads (FASTQ) of the isolate".to_string());
    }
    for i in inputs {
        check_input(i)?;
    }
    for (name, v) in [("min_identity", opt.min_identity), ("min_coverage", opt.min_coverage), ("min_partial_coverage", opt.min_partial_coverage)] {
        if !(0.0..=1.0).contains(&v) {
            return Err(format!("Invalid {name}: {v}. It is a fraction, from 0 to 1 (0.9 is the usual value)"));
        }
    }
    if let Some(o) = &opt.output {
        let mut all: Vec<&str> = inputs.iter().map(String::as_str).collect();
        all.push(catalog);
        all.extend(catalog_table);
        not_an_input(o, &all)?;
    }
    let mode = match opt.mode {
        AmrMode::Assembly => AmrMode::Assembly,
        AmrMode::Reads => AmrMode::Reads,
        AmrMode::Auto => {
            let mut kinds = Vec::new();
            for i in inputs {
                if let Some(k) = sniff_kind(i)? {
                    kinds.push((i, k));
                }
            }
            let fq = kinds.iter().filter(|x| x.1 == SeqKind::Fastq).count();
            if fq > 0 && fq < kinds.len() {
                let fa = kinds.iter().find(|x| x.1 == SeqKind::Fasta).map(|x| x.0.as_str()).unwrap_or("");
                return Err(format!(
                    "The inputs mix reads (FASTQ) and an assembly (FASTA, {fa}): give either the contigs of the isolate or its reads, or say which with the mode"
                ));
            }
            if fq > 0 {
                AmrMode::Reads
            } else {
                AmrMode::Assembly
            }
        }
    };
    let (genes, unmatched) = load_amr_catalog(catalog, catalog_table, opt.include_plus)?;
    let idx = SeedIndex::build(&genes);
    let mut res = AmrResult {
        inputs: inputs.to_vec(),
        catalog: catalog.to_string(),
        catalog_table: catalog_table.map(str::to_string),
        mode: if mode == AmrMode::Reads { "reads" } else { "assembly" }.to_string(),
        sequences: 0,
        bases: 0,
        catalog_genes: genes.len(),
        catalog_unmatched: unmatched,
        alignments: 0,
        hits: Vec::new(),
    };
    if mode == AmrMode::Assembly {
        for_each_batch(inputs, 256, 64 << 20, |recs, _| {
            res.sequences += recs.len() as u64;
            res.bases += recs.iter().map(|r| r.seq.len() as u64).sum::<u64>();
            let found: Vec<(Vec<RawHit>, u64)> = recs.par_iter().map(|r| scan_contig(&r.seq, &idx, &genes, opt)).collect();
            for (rec, (hits, aligned)) in recs.iter().zip(found) {
                res.alignments += aligned;
                let name = String::from_utf8_lossy(rec.head.split(|c| c.is_ascii_whitespace()).next().unwrap_or(&[])).to_string();
                for h in hits {
                    let g = &genes[h.gene as usize];
                    let mut x = hit_of(g, h.identity, h.coverage);
                    x.contig = name.clone();
                    x.start = h.start as u64 + 1;
                    x.end = h.end as u64;
                    x.strand = if h.minus { '-' } else { '+' };
                    x.aligned_length = h.aln.columns();
                    x.mismatches = h.aln.mismatches;
                    x.gap_bases = h.aln.ins + h.aln.del;
                    x.status = if h.exact(g.seq.len()) {
                        "exact"
                    } else if h.coverage >= opt.min_coverage {
                        "aligned"
                    } else {
                        "partial_contig_end"
                    }
                    .to_string();
                    res.hits.push(x);
                }
            }
            Ok(())
        })?;
    } else {
        let mut piles: HashMap<u32, Vec<Pile>> = HashMap::new();
        // For each read, the genes it matches best (fewest differences) and
        // all the genes it aligned to, as sets with how many reads had each:
        // they are what tells one locus from another without an assembly.
        let mut best_sets: HashMap<Vec<u32>, u64> = HashMap::new();
        let mut all_sets: HashMap<Vec<u32>, u64> = HashMap::new();
        for_each_batch(inputs, 50_000, 16 << 20, |recs, _| {
            res.sequences += recs.len() as u64;
            res.bases += recs.iter().map(|r| r.seq.len() as u64).sum::<u64>();
            let found: Vec<(Vec<ReadAln>, u32)> =
                recs.par_iter().map_init(Vec::new, |rc, r| read_alignments(&r.seq, &idx, &genes, rc)).collect();
            for (alns, tried) in found {
                res.alignments += u64::from(tried);
                let Some(fewest) = alns.iter().map(|a| a.2).min() else { continue };
                let mut best: Vec<u32> = alns.iter().filter(|a| a.2 == fewest).map(|a| a.0).collect();
                best.sort_unstable();
                best.dedup();
                *best_sets.entry(best).or_insert(0) += 1;
                let mut all: Vec<u32> = alns.iter().map(|a| a.0).collect();
                all.sort_unstable();
                all.dedup();
                *all_sets.entry(all).or_insert(0) += 1;
                for (g, ops, diffs, columns) in alns {
                    let pile = piles.entry(g).or_insert_with(|| vec![[0u32; 18]; genes[g as usize].seq.len()]);
                    let best = diffs == fewest;
                    let close = 20 * diffs <= columns;
                    for (p, s) in ops {
                        let cell = &mut pile[p as usize];
                        let s = usize::from(s);
                        cell[s] += 1;
                        if best {
                            cell[6 + s] += 1;
                        }
                        if close {
                            cell[12 + s] += 1;
                        }
                    }
                }
            }
            Ok(())
        })?;
        let mut cands: Vec<(AmrHit, u32)> = piles
            .par_iter()
            .map(|(&g, pile)| {
                let (mut h, mixed) = evaluate_pileup(&genes[g as usize], pile, opt.min_depth.max(1));
                h.mixed = mixed;
                (h, g)
            })
            .filter(|(h, _)| h.identity >= opt.min_identity && h.coverage >= opt.min_coverage)
            .collect();
        // Best first: exact, then identity x coverage, identity, the longer
        // gene (it explains more of the sample: a shorter allele can look
        // exact when the reads are clipped at its indel), depth, and the
        // earlier gene of the catalog.
        cands.sort_by(|(a, ga), (b, gb)| {
            let ea = a.identity == 1.0 && a.coverage == 1.0;
            let eb = b.identity == 1.0 && b.coverage == 1.0;
            eb.cmp(&ea)
                .then((b.identity * b.coverage).total_cmp(&(a.identity * a.coverage)))
                .then(b.identity.total_cmp(&a.identity))
                .then(b.aligned_length.cmp(&a.aligned_length))
                .then(b.depth.unwrap_or(0.0).total_cmp(&a.depth.unwrap_or(0.0)))
                .then(ga.cmp(gb))
        });
        // One locus, one line. A gene half of whose reads also aligned to a
        // gene already kept belongs to that gene's group: it is the same
        // locus seen through a relative (another allele, a homolog), or a
        // novel allele that different catalog alleles each match in part.
        // It is reported only as the group's second allele, when the reads
        // of the kept gene are mixed (two alleles are really there) and some
        // of its own best-matching reads -- at least 3, and 10% of them -- are
        // not matched as well by a kept gene.
        let sum_sets = |sets: &HashMap<Vec<u32>, u64>, g: u32, other: Option<u32>, kept: Option<&[u32]>| -> u64 {
            sets.iter()
                .filter(|(s, _)| s.binary_search(&g).is_ok())
                .filter(|(s, _)| other.is_none_or(|o| s.binary_search(&o).is_ok()))
                .filter(|(s, _)| kept.is_none_or(|k| !k.iter().any(|x| s.binary_search(x).is_ok())))
                .map(|(_, &n)| n)
                .sum()
        };
        // (hit, gene, this group already has its second allele)
        let mut kept: Vec<(AmrHit, u32, bool)> = Vec::new();
        for (h, g) in cands {
            let all = sum_sets(&all_sets, g, None, None);
            let group: Vec<usize> = (0..kept.len()).filter(|&i| 2 * sum_sets(&all_sets, g, Some(kept[i].1), None) >= all).collect();
            if group.is_empty() {
                kept.push((h, g, false));
                continue;
            }
            let ids: Vec<u32> = kept.iter().map(|k| k.1).collect();
            let own = sum_sets(&best_sets, g, None, None);
            let exclusive = sum_sets(&best_sets, g, None, Some(&ids));
            let open = group.iter().all(|&i| !kept[i].2) && group.iter().any(|&i| kept[i].0.mixed);
            if open && exclusive >= 3 && 10 * exclusive >= own {
                for &i in &group {
                    kept[i].2 = true;
                }
                kept.push((h, g, true));
            }
        }
        let mut hits: Vec<AmrHit> = kept.into_iter().map(|(h, _, _)| h).collect();
        for h in hits.iter_mut() {
            h.status = if h.identity == 1.0 && h.coverage == 1.0 { "exact" } else { "aligned" }.to_string();
        }
        hits.sort_by(|a, b| a.class.cmp(&b.class).then(a.gene.cmp(&b.gene)));
        res.hits = hits;
    }
    if let Some(o) = &opt.output {
        let mut cleanup = Cleanup(vec![o.clone()]);
        let mut w = create_output(o)?;
        w.write_all(res.tsv().as_bytes()).map_err(|e| write_error(o, e))?;
        w.termina().map_err(|e| write_error(o, e))?;
        cleanup.disarm();
    }
    Ok(res)
}

impl AmrHit {
    fn json(&self) -> String {
        format!(
            "{{\"gene\":{},\"family\":{},\"product\":{},\"scope\":{},\"type\":{},\"subtype\":{},\"class\":{},\"subclass\":{},\"identity\":{},\"coverage\":{},\"contig\":{},\"start\":{},\"end\":{},\"strand\":{},\"gene_length\":{},\"aligned_length\":{},\"mismatches\":{},\"gap_bases\":{},\"status\":{},\"depth\":{},\"mixed\":{},\"protein\":{},\"nucleotide\":{}}}",
            js(&self.gene),
            js(&self.family),
            js(&self.product),
            js(&self.scope),
            js(&self.element_type),
            js(&self.subtype),
            js(&self.class),
            js(&self.subclass),
            jf(self.identity),
            jf(self.coverage),
            js(&self.contig),
            self.start,
            self.end,
            js(&self.strand.to_string()),
            self.gene_length,
            self.aligned_length,
            self.mismatches,
            self.gap_bases,
            js(&self.status),
            jo(self.depth),
            self.mixed,
            js(&self.protein),
            js(&self.nucleotide)
        )
    }
}

impl AmrResult {
    pub fn json(&self) -> String {
        let hits: Vec<String> = self.hits.iter().map(AmrHit::json).collect();
        format!(
            "{{\"inputs\":{},\"catalog\":{},\"catalog_table\":{},\"mode\":{},\"sequences\":{},\"bases\":{},\"catalog_genes\":{},\"catalog_unmatched\":{},\"alignments\":{},\"hits\":[{}]}}",
            jlist(&self.inputs),
            js(&self.catalog),
            self.catalog_table.as_deref().map_or_else(|| "null".to_string(), js),
            js(&self.mode),
            self.sequences,
            self.bases,
            self.catalog_genes,
            self.catalog_unmatched,
            self.alignments,
            hits.join(",")
        )
    }

    /// The table of hits (what `output` receives).
    pub fn tsv(&self) -> String {
        let mut o = String::from(
            "gene\tfamily\tproduct\tscope\ttype\tsubtype\tclass\tsubclass\tidentity\tcoverage\tcontig\tstart\tend\tstrand\tgene_length\taligned_length\tmismatches\tgap_bases\tstatus\tdepth\tmixed\tprotein\tnucleotide\n",
        );
        for h in &self.hits {
            o.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                clean_cell(&h.gene),
                clean_cell(&h.family),
                clean_cell(&h.product),
                h.scope,
                h.element_type,
                h.subtype,
                h.class,
                h.subclass,
                h.identity,
                h.coverage,
                clean_cell(&h.contig),
                h.start,
                h.end,
                h.strand,
                h.gene_length,
                h.aligned_length,
                h.mismatches,
                h.gap_bases,
                h.status,
                opt_num(h.depth, 1),
                h.mixed,
                h.protein,
                h.nucleotide
            ));
        }
        o
    }

    pub fn resumo(&self) -> String {
        let mut o = format!(
            "Resistance genes ({}): {} sequences, {} bases searched against {} catalog genes\n",
            self.mode,
            thousands(self.sequences),
            thousands(self.bases),
            thousands(self.catalog_genes as u64)
        );
        if self.catalog_table.is_none() {
            o.push_str("  (no ReferenceGeneCatalog.txt given: class, subclass and scope are unknown)\n");
        } else if self.catalog_unmatched > 0 {
            o.push_str(&format!("  {} catalog genes have no line in the table: their class is unknown\n", self.catalog_unmatched));
        }
        if self.hits.is_empty() {
            o.push_str("  no gene of the catalog passes the thresholds\n");
            return o;
        }
        o.push_str("  gene                 identity  coverage  class / subclass                      where\n");
        for h in &self.hits {
            let place = if h.contig.is_empty() {
                format!("reads, depth {:.1}", h.depth.unwrap_or(0.0))
            } else {
                format!("{}:{}-{} ({})", h.contig, h.start, h.end, h.strand)
            };
            let class = if h.subclass.is_empty() || h.subclass == h.class { h.class.clone() } else { format!("{} / {}", h.class, h.subclass) };
            let flag = match h.status.as_str() {
                "partial_contig_end" => "  [cut by the contig end]",
                _ if h.mixed => "  [the reads hold more than one allele]",
                _ => "",
            };
            o.push_str(&format!(
                "  {:<20} {:>8}  {:>8}  {:<36}  {place}{flag}{}\n",
                h.gene,
                format!("{:.2}%", 100.0 * h.identity),
                format!("{:.2}%", 100.0 * h.coverage),
                class,
                if h.scope == "plus" { "  (plus)" } else { "" }
            ));
        }
        o
    }
}
