// =====================================================================
//  16S rRNA AMPLICONS: FROM PAIRED READS TO SEQUENCE VARIANTS, TAXONOMY
//  AND DIVERSITY; AND THE TUMOR MUTATIONAL BURDEN OF A SOMATIC VCF
//
//  Public tools (every one returns Result<_, String>; the messages name the
//  file and say how to fix it):
//
//   1. `amplicon_asvs`      paired FASTQ of each sample -> amplicon sequence
//                           variants (ASV FASTA), the feature table (ASV x
//                           sample counts) and the per-sample read tracking;
//   2. `train_classifier`   the naive Bayesian classifier of Wang et al. 2007
//                           (the RDP classifier) trained from a reference
//                           FASTA whose headers carry the lineage, saved to a
//                           model file; `classify_asvs` assigns each ASV a
//                           lineage with a bootstrap confidence per rank and,
//                           optionally, the species by exact match;
//   3. `alpha_diversity`, `beta_diversity`, `taxa_composition`: alpha indices
//                           with rarefaction and rarefaction curves; Bray-
//                           Curtis / Jaccard / k-mer distances, principal
//                           coordinates and PERMANOVA; relative abundance per
//                           taxon at a chosen rank;
//   4. `tumor_mutational_burden`  somatic variants per megabase of the
//                           callable region, optionally coding-only or
//                           non-synonymous-only.
//
//  PAIR MERGING (as vsearch --fastq_mergepairs, Edgar & Flyvbjerg 2015).
//  Read 2 is reverse-complemented and every diagonal that leaves at least
//  `min_overlap` bases of overlap is scored with quality-aware log-odds:
//  a column where the two bases agree scores log2(Pa / 0.25), one where they
//  differ log2((1 - Pa) / 0.75), with Pa = (1-p1)(1-p2) + p1 p2 / 3 the
//  probability that two reads of one template show the same base (p = the
//  error probability of each quality, capped at 0.75); an `N` scores 0. The
//  diagonals are tried from the longest overlap down and a diagonal is
//  abandoned as soon as it cannot beat the best one (each column scores at
//  most 2), so a 2 x 250 pair costs about one pass over the overlap. The best
//  diagonal must score at least `min_score` (16, vsearch's default) and have
//  at most `max_diffs` mismatches (10) and `max_diff_pct` percent (100).
//  A diagonal where the fragment is shorter than a read (read-through into
//  the adapter) is "staggered" and rejected unless `allow_stagger`; when
//  allowed, the overhangs are cut. In the overlap the merged base is the
//  agreed base, with the posterior quality -10 log10(p1 p2 / 3 /
//  (1 - p1 - p2 + 4 p1 p2 / 3)); on a disagreement it is the base of the
//  higher quality (read 1 on a tie), with the posterior p1 (1 - p2 / 3) /
//  (p1 + p2 - 4 p1 p2 / 3). Merged qualities are capped at `max_q_out` (41).
//
//  FILTERING. Optional primer removal (IUPAC codes; the forward primer near
//  the start, the reverse complement of the reverse primer near the end, at
//  most `primer_max_error_rate` mismatches, no indels), then expected errors
//  (the sum of 10^(-Q/10) over the merged read, as --fastq_maxee), length
//  range and `N`s.
//
//  DENOISING (UNOISE3, Edgar 2016). The filtered reads of all samples are
//  dereplicated. Unique sequences with at least `min_size` reads (8), from
//  the most abundant down (ties by sequence), become a centroid unless a
//  centroid C already exists with d = edit distance(Q, C) and
//  abundance(Q) / abundance(C) <= 1 / 2^(alpha d + 1) (alpha = 2); then Q is
//  absorbed into the closest such centroid (fewest differences, then the
//  more abundant). With `map_rare`, the uniques below `min_size` are
//  assigned by the same rule (they never become centroids) -- so the feature
//  table counts every read that the model explains as an error of an ASV.
//  The edit distance is computed in a band of width 2 d_max + 1, where d_max
//  is the largest d the abundances allow, and stops as soon as it exceeds
//  d_max: the comparisons that matter are the near-identical ones.
//
//  CHIMERAS (de novo, UCHIME2/UCHIME3 "denoised" mode). The abundance of a
//  centroid here is its UNOISE3 cluster size (the reads of the uniques with
//  at least min_size reads that it absorbed, as the sizes --cluster_unoise
//  writes), not the rare reads mapped later. From the most abundant
//  centroid down, the candidate parents of Q are the centroids already
//  judged non-chimeric with at least `abskew` (16, as --uchime3_denovo)
//  times its abundance. Q is aligned end to end to each
//  candidate (banded edit-distance alignment); for every crossover point x
//  the two-parent model is "the best parent of Q[..x] + the best parent of
//  Q[x..]". The votes of a segment are the differences the other parent has
//  there (net of the model's own): the evidence that the crossover is real.
//  Q is chimeric when the best model (fewest differences, then the most
//  votes on its weaker side) has at most `max_model_diffs` differences (0:
//  a perfect chimera), each segment has at least `min_diffs` votes (3, as
//  --mindiffs) and no single parent explains Q as well. Chimeric centroids
//  and their reads leave the table. Three-parent chimeras are not modelled.
//  The search over parents and crossovers is exhaustive, so it also finds
//  crossovers whose diagnostic differences all lie within a few bases of an
//  end, which UCHIME's chunk-wise search for parents does not reach (on the
//  MiSeq SOP data: every vsearch --uchime3_denovo call, plus 26 such).
//
//  ASV NAMES: ASV_1, ASV_2, ... by decreasing reads in the final table (ties
//  by sequence), so the names are stable for a given input and options.
//
//  TAXONOMY (Wang et al. 2007). Words = the distinct 8-mers (A/C/G/T only;
//  U reads as T) of a sequence. With N training sequences, n(w) of them
//  holding word w, and M_G sequences in genus G of which m(w,G) hold w:
//      P_w = (n(w) + 0.5) / (N + 1),  P(w|G) = (m(w,G) + P_w) / (M_G + 1),
//  and a query scores sum_w log P(w|G). A "genus" is each distinct full
//  lineage (so a lineage that stops at the family is a class of its own).
//  Only the pairs (w, G) with m(w,G) > 0 are stored: log P(w|G) for the
//  others is log P_w - log(M_G + 1), so the score of G is a constant, minus
//  |W| log(M_G + 1), plus the stored bonus log(m + P_w) - log P_w of the
//  words G holds. The assignment uses all the words of the query; then 100
//  bootstrap trials of |W|/8 words drawn with replacement, and the
//  confidence at each rank is the fraction of trials whose winner has the
//  same taxon at that rank (as DADA2's assignTaxonomy). Ties go to one of
//  the tied genera at random, with a seed derived from the sequence, so the
//  result does not depend on the order or the threads. The lineage is kept
//  down to the deepest rank with confidence >= `min_confidence` (0.8 as the
//  RDP classifier; 0.5 is DADA2's minBoot = 50, advised for short reads).
//
//  THE MODEL FILE (version 1; gzip of the following; integers LEB128
//  varints unless said otherwise):
//      "MMGTAX01", k (1 byte), number of ranks + names, N, number of
//      genera, then per genus M_G and its lineage; the source file name;
//      then for each of the 4^k words: n(w), number of genera holding it,
//      and the (genus - previous genus, m(w,G)) pairs; "MMGTAXEN" and the
//      FNV-1a 64 (8 bytes, little-endian) of every byte before it.
//  Strings are a varint length and UTF-8 bytes.
//
//  SPECIES (as DADA2 addSpecies). An ASV gets a species only from the
//  reference sequences that contain it exactly (100% identity over the
//  whole ASV; both strands with `try_reverse_complement`), whose genus is
//  the genus the classifier assigned, and only when they name at most
//  `max_species` species (1: unambiguous; more are joined by '/').
//  Reference headers: ">ID Genus species" (DADA2 assignSpecies files) or a
//  lineage whose last two fields are the genus and the species.
//
//  DIVERSITY. Alpha: observed, Chao1 (bias-corrected, S + F1 (F1 - 1) /
//  (2 (F2 + 1))), Shannon (natural log), Simpson (1 - sum p^2), Pielou
//  (H / ln S), Good's coverage (1 - F1 / N). Rarefaction draws without
//  replacement with a seed derived from the seed and the sample name; the
//  rarefaction curve is the exact expectation of the observed ASVs at each
//  depth (Hurlbert 1971), so it needs no seed. Beta: Bray-Curtis on the
//  counts (or proportions with `relative`), Jaccard on presence, or "kmer":
//  Bray-Curtis between the abundance-weighted 6-mer profiles of the
//  samples' ASVs, a phylogeny-free distance that sees close relatives as
//  close (UniFrac is not implemented). PCoA: eigen-decomposition (Jacobi)
//  of the double-centred matrix -0.5 D^2; the proportion explained is each
//  positive eigenvalue over the sum of the positive ones (scikit-bio's
//  convention; negative eigenvalues are reported and get no axis). PERMANOVA: pseudo-F = (SS_A / (a - 1)) /
//  (SS_W / (N - a)), p = (permutations with F >= observed + 1) /
//  (permutations + 1), labels shuffled with a fixed seed.
//
//  TUMOR MUTATIONAL BURDEN. Each ALT allele of the VCF is a candidate
//  (common prefix and suffix with REF trimmed first, so the type and place
//  are right in multi-allelic records). Counted when: FILTER is PASS (with
//  `pass_only`); the tumor depth and variant allele fraction pass the
//  thresholds (from FORMAT AF / AD / DP, Strelka's tier counts, or INFO
//  AF / DP; a variant without them is kept and counted apart); its REF span
//  touches the BED region; with `coding_only`, it touches a CDS of the
//  annotation; with `nonsynonymous_only`, it changes a protein in at least
//  one transcript (codons translated with the CDS phase and the reference
//  genome; an indel in a CDS always counts; a variant whose codon cannot be
//  completed counts as coding_unknown). TMB = counted / Mb, where the Mb is
//  the BED (merged), or the BED intersected with the CDS when coding_only or
//  nonsynonymous_only.
// =====================================================================

use crate::externo::{abre_leitura, Escritor};
use crate::multimerge_core::multi_merge_sort_stable_by;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, ErrorKind, Read, Write};
use std::sync::OnceLock;

// ---------------------------------------------------------------------
//  SMALL HELPERS: messages, JSON, numbers, random numbers
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

/// Writes a whole text file through the project's writer (the extension
/// chooses plain, .gz or .lz4).
fn write_text(path: &str, text: &str) -> Result<(), String> {
    let mut w = create_output(path)?;
    w.write_all(text.as_bytes()).map_err(|e| write_error(path, e))?;
    w.termina().map_err(|e| write_error(path, e))
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

fn jnums(items: &[f64]) -> String {
    let v: Vec<String> = items.iter().map(|x| jf(*x)).collect();
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

fn pct(part: u64, whole: u64) -> String {
    if whole == 0 {
        return "-".to_string();
    }
    format!("{:.1}%", 100.0 * part as f64 / whole as f64)
}

/// A number for a table: integers without decimals, the rest with up to 6
/// significant decimals.
fn num(x: f64) -> String {
    if !x.is_finite() {
        return "NA".to_string();
    }
    if x == x.trunc() && x.abs() < 1e15 {
        return format!("{}", x as i64);
    }
    let s = format!("{x:.6}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    s.to_string()
}

fn fnv1a(mut h: u64, data: &[u8]) -> u64 {
    for &b in data {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;

/// SplitMix64: small, fast, good enough for subsampling and permutations,
/// and the same on every platform.
#[derive(Clone, Debug)]
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed)
    }
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// Uniform in 0..n (n > 0), by rejection: no modulo bias.
    pub fn below(&mut self, n: u64) -> u64 {
        let zone = u64::MAX - (u64::MAX % n);
        loop {
            let x = self.next_u64();
            if x < zone {
                return x % n;
            }
        }
    }
}

fn seed_for(seed: u64, key: &[u8]) -> u64 {
    let mut r = Rng::new(seed ^ fnv1a(FNV_OFFSET, key));
    r.next_u64()
}

fn sample_text(piece: &[u8]) -> String {
    let mut o = String::new();
    for &b in piece.iter().take(40) {
        match b {
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

/// The first word of a FASTA header, without a usearch-style `;size=N`.
fn fasta_id(head: &[u8]) -> String {
    let first = head.split(|c| c.is_ascii_whitespace()).next().unwrap_or(head);
    let t = String::from_utf8_lossy(first).to_string();
    match t.find(";size=") {
        Some(p) => t[..p].to_string(),
        None => t.trim_end_matches(';').to_string(),
    }
}

/// Stable sort of indices by a comparator (the project's stable sort).
fn sorted_indices<F>(n: usize, cmp: F) -> Vec<u32>
where
    F: Fn(&u32, &u32) -> Ordering + Sync,
{
    let mut idx: Vec<u32> = (0..n as u32).collect();
    multi_merge_sort_stable_by(&mut idx, cmp);
    idx
}

// ---------------------------------------------------------------------
//  BASES
// ---------------------------------------------------------------------

/// The base in upper case (U reads as T); anything else becomes `N`.
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
    t[b'U' as usize] = b'T';
    t[b'u' as usize] = b'T';
    t
};

/// 0..3 for A, C, G, T (either case, U as T); 4 for anything else.
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
    t[b'U' as usize] = 3;
    t[b'u' as usize] = 3;
    t
};

fn complement(b: u8) -> u8 {
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

/// Reverse complement of upper-case bases (anything else becomes `N`).
pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| complement(UPPER[b as usize])).collect()
}

/// Bit mask of the bases an IUPAC letter stands for: A=1, C=2, G=4, T=8.
fn iupac_mask(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => 0,
    }
}

/// The mask of the complementary bases (A<->T, C<->G).
fn comp_mask(m: u8) -> u8 {
    ((m & 1) << 3) | ((m & 8) >> 3) | ((m & 2) << 1) | ((m & 4) >> 1)
}

// ---------------------------------------------------------------------
//  READING FASTA AND FASTQ
// ---------------------------------------------------------------------

#[derive(Clone, Debug, Default)]
pub struct SeqRec {
    /// The header line without its first character (`>` or `@`).
    pub head: Vec<u8>,
    pub seq: Vec<u8>,
    /// FASTQ qualities as written (ASCII); empty for FASTA.
    pub qual: Vec<u8>,
    /// Line where the record starts, from 1.
    pub line: u64,
}

pub struct SeqIn {
    path: String,
    r: Box<dyn BufRead>,
    buf: Vec<u8>,
    line_no: u64,
    /// Some(true) = FASTA, Some(false) = FASTQ; None before the first record.
    fasta: Option<bool>,
    next_head: Option<(Vec<u8>, u64)>,
}

impl SeqIn {
    pub fn open(path: &str) -> Result<SeqIn, String> {
        check_input(path)?;
        let r = abre_leitura(path).map_err(|e| format!("{e}. Check the path and that you are allowed to read the file"))?;
        Ok(SeqIn { path: path.to_string(), r, buf: Vec::new(), line_no: 0, fasta: None, next_head: None })
    }

    pub fn is_fastq(&self) -> bool {
        self.fasta == Some(false)
    }

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
        } else {
            " Give sequences (FASTA or FASTQ, plain, .gz or .lz4)"
        };
        format!(
            "{}: this is not a FASTA nor a FASTQ file: line {} starts with \"{}\", and the first line of a FASTA file starts with '>' and of a FASTQ file with '@'.{hint}",
            self.path,
            self.line_no,
            sample_text(&self.buf)
        )
    }

    /// Reads the next record into `rec`; `false` at the end of the file.
    pub fn next(&mut self, rec: &mut SeqRec) -> Result<bool, String> {
        rec.head.clear();
        rec.seq.clear();
        rec.qual.clear();
        if self.fasta.is_none() {
            loop {
                if !self.line()? {
                    return Ok(false);
                }
                let t = trim(&self.buf);
                match t.first() {
                    None => continue,
                    Some(b'>') => self.fasta = Some(true),
                    Some(b'@') => self.fasta = Some(false),
                    Some(_) => return Err(self.not_sequences()),
                }
                self.next_head = Some((t[1..].to_vec(), self.line_no));
                break;
            }
        }
        if self.fasta == Some(true) {
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
            return Ok(true);
        }
        // FASTQ
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
        let cut = |s: &SeqIn, what: &str| {
            format!(
                "{}: truncated FASTQ: the file ends in the middle of the record that starts at line {start} ({what}). The file was probably cut during a copy or a download: get it again",
                s.path
            )
        };
        if !self.line()? {
            return Err(cut(self, "the sequence, '+' and quality lines are missing"));
        }
        rec.seq.extend(self.buf.iter().filter(|c| !c.is_ascii_whitespace()));
        if !self.line()? {
            return Err(cut(self, "the '+' and quality lines are missing"));
        }
        if trim(&self.buf).first() != Some(&b'+') {
            return Err(format!(
                "{}, line {}: expected the '+' line of the FASTQ record that starts at line {start}, found \"{}\". A FASTQ record has exactly 4 lines (@name, sequence, +, qualities)",
                self.path,
                self.line_no,
                sample_text(&self.buf)
            ));
        }
        if !self.line()? {
            if rec.seq.is_empty() {
                return Ok(true);
            }
            return Err(cut(self, "the quality line is missing"));
        }
        rec.qual.extend(self.buf.iter().filter(|c| !c.is_ascii_whitespace()));
        if rec.qual.len() != rec.seq.len() {
            return Err(format!(
                "{}, line {}: the quality line has {} characters, but the sequence of the record that starts at line {start} has {} bases. The file is damaged or cut short: check it with the validator, or get it again",
                self.path,
                self.line_no,
                rec.qual.len(),
                rec.seq.len()
            ));
        }
        Ok(true)
    }
}

/// Every record of a FASTA (or FASTQ) file: (id, upper-case sequence).
pub fn read_fasta(path: &str) -> Result<Vec<(String, Vec<u8>)>, String> {
    let mut r = SeqIn::open(path)?;
    let mut rec = SeqRec::default();
    let mut out = Vec::new();
    while r.next(&mut rec)? {
        out.push((fasta_id(&rec.head), rec.seq.iter().map(|&b| UPPER[b as usize]).collect()));
    }
    Ok(out)
}

// ---------------------------------------------------------------------
//  PAIR MERGING
// ---------------------------------------------------------------------

/// Highest quality value handled (Phred 93 is '~' at offset 33).
const QN: usize = 94;

struct QualTables {
    /// log2-odds of a column where the bases agree / differ, by (q1, q2).
    agree: Vec<f64>,
    differ: Vec<f64>,
    /// Posterior quality of the merged base, by (q1, q2): agreeing, and
    /// differing (index = (higher q, lower q)).
    q_agree: Vec<u8>,
    q_differ: Vec<u8>,
}

fn err_prob(q: usize) -> f64 {
    10f64.powf(-(q as f64) / 10.0).min(0.75)
}

fn phred_of(p: f64) -> u8 {
    if p <= 0.0 {
        return (QN - 1) as u8;
    }
    let q = (-10.0 * p.log10()).round();
    q.clamp(0.0, (QN - 1) as f64) as u8
}

fn tables() -> &'static QualTables {
    static T: OnceLock<QualTables> = OnceLock::new();
    T.get_or_init(|| {
        let mut t = QualTables { agree: vec![0.0; QN * QN], differ: vec![0.0; QN * QN], q_agree: vec![0; QN * QN], q_differ: vec![0; QN * QN] };
        for a in 0..QN {
            for b in 0..QN {
                let (p1, p2) = (err_prob(a), err_prob(b));
                let same = (1.0 - p1) * (1.0 - p2) + p1 * p2 / 3.0;
                t.agree[a * QN + b] = (same / 0.25).log2();
                t.differ[a * QN + b] = ((1.0 - same).max(1e-300) / 0.75).log2();
                let pa = (p1 * p2 / 3.0) / (1.0 - p1 - p2 + 4.0 * p1 * p2 / 3.0);
                t.q_agree[a * QN + b] = phred_of(pa);
                // a = the higher quality (lower error) base, b = the other
                let (lo, hi) = if p1 <= p2 { (p1, p2) } else { (p2, p1) };
                let pd = lo * (1.0 - hi / 3.0) / (lo + hi - 4.0 * lo * hi / 3.0);
                t.q_differ[a * QN + b] = phred_of(pd);
            }
        }
        t
    })
}

/// Posterior quality of a merged base when both reads agree (Edgar &
/// Flyvbjerg 2015), before the `max_q_out` cap.
pub fn merged_quality_agree(q1: u8, q2: u8) -> u8 {
    tables().q_agree[(q1 as usize).min(QN - 1) * QN + (q2 as usize).min(QN - 1)]
}

/// Posterior quality when they disagree (the base of the higher quality is
/// kept), before the cap.
pub fn merged_quality_differ(q1: u8, q2: u8) -> u8 {
    tables().q_differ[(q1 as usize).min(QN - 1) * QN + (q2 as usize).min(QN - 1)]
}

#[derive(Clone, Debug)]
pub struct MergeOptions {
    /// Shortest overlap accepted (vsearch --fastq_minovlen, 10).
    pub min_overlap: usize,
    /// Most mismatches in the overlap (--fastq_maxdiffs, 10).
    pub max_diffs: usize,
    /// Most mismatches as a percentage of the overlap (--fastq_maxdiffpct, 100).
    pub max_diff_pct: f64,
    /// Lowest alignment score (log2-odds) of the chosen overlap (16).
    pub min_score: f64,
    /// Accept pairs whose fragment is shorter than a read, cutting the
    /// overhangs (--fastq_allowmergestagger).
    pub allow_stagger: bool,
    /// Cut each read at its first base with quality <= this (--fastq_truncqual).
    pub trunc_qual: Option<u8>,
    /// Cut read 1 / read 2 to this many bases first (0 = keep whole); a
    /// shorter read is kept as it is.
    pub trunc_len_r1: usize,
    pub trunc_len_r2: usize,
    /// Quality encoding offset (33; 64 for old Illumina files).
    pub phred_offset: u8,
    /// Highest quality written for a merged base (--fastq_qmaxout, 41).
    pub max_q_out: u8,
}

impl Default for MergeOptions {
    fn default() -> Self {
        MergeOptions {
            min_overlap: 10,
            max_diffs: 10,
            max_diff_pct: 100.0,
            min_score: 16.0,
            allow_stagger: false,
            trunc_qual: None,
            trunc_len_r1: 0,
            trunc_len_r2: 0,
            phred_offset: 33,
            max_q_out: 41,
        }
    }
}

/// A merged pair: the bases and the numeric qualities (no offset).
#[derive(Clone, Debug, PartialEq)]
pub struct MergedRead {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub overlap: usize,
    pub diffs: usize,
    pub score: f64,
    pub staggered: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MergeFail {
    /// A read is shorter than the minimum overlap (after truncation).
    TooShort,
    /// No overlap scores `min_score`: the reads do not overlap (or not
    /// enough), or the overlap is of very low quality.
    NoAlignment,
    /// The best overlap has more mismatches than allowed.
    TooManyDiffs,
    /// The best overlap means the fragment is shorter than a read.
    Staggered,
}

impl MergeFail {
    pub fn name(self) -> &'static str {
        match self {
            MergeFail::TooShort => "too_short",
            MergeFail::NoAlignment => "no_alignment",
            MergeFail::TooManyDiffs => "too_many_diffs",
            MergeFail::Staggered => "staggered",
        }
    }
}

/// The numeric qualities of a FASTQ quality string; `Err(position)` of the
/// first character below the offset.
fn numeric_quals(q: &[u8], offset: u8) -> Result<Vec<u8>, usize> {
    let mut o = Vec::with_capacity(q.len());
    for (i, &c) in q.iter().enumerate() {
        if c < offset {
            return Err(i);
        }
        o.push((c - offset).min((QN - 1) as u8));
    }
    Ok(o)
}

/// Where to cut a read (its new length) by the truncation options.
fn trunc_point(qual: &[u8], trunc_len: usize, trunc_qual: Option<u8>) -> usize {
    let mut n = qual.len();
    if trunc_len > 0 {
        n = n.min(trunc_len);
    }
    if let Some(t) = trunc_qual {
        if let Some(p) = qual[..n].iter().position(|&q| q <= t) {
            n = p;
        }
    }
    n
}

/// Merges one pair. `s1`/`s2` are the reads as sequenced (read 2 is
/// reverse-complemented here), `q1`/`q2` their NUMERIC qualities.
pub fn merge_pair(s1: &[u8], q1: &[u8], s2: &[u8], q2: &[u8], opt: &MergeOptions) -> Result<MergedRead, MergeFail> {
    let t = tables();
    let n1 = trunc_point(q1, opt.trunc_len_r1, opt.trunc_qual).min(s1.len());
    let n2 = trunc_point(q2, opt.trunc_len_r2, opt.trunc_qual).min(s2.len());
    let min_ov = opt.min_overlap.max(1);
    if n1 < min_ov || n2 < min_ov {
        return Err(MergeFail::TooShort);
    }
    let a: Vec<u8> = s1[..n1].iter().map(|&b| UPPER[b as usize]).collect();
    let qa: Vec<u8> = q1[..n1].to_vec();
    let b: Vec<u8> = s2[..n2].iter().rev().map(|&c| complement(UPPER[c as usize])).collect();
    let qb: Vec<u8> = q2[..n2].iter().rev().copied().collect();
    let (l1, l2) = (n1 as i64, n2 as i64);
    // diagonal d: read2rc[j] faces read1[j + d]; overlap = [max(0,d), min(l1, d+l2))
    let overlap = |d: i64| (l1.min(d + l2) - d.max(0)).max(0);
    let d_lo = -(l2 - min_ov as i64);
    let d_hi = l1 - min_ov as i64;
    // Scores diagonal d; abandons it (None) as soon as it cannot beat `bound`
    // (each remaining column scores at most 2).
    let eval = |d: i64, bound: Option<f64>| -> Option<(f64, usize)> {
        let start = d.max(0) as usize;
        let end = l1.min(d + l2) as usize;
        let mut score = 0.0;
        let mut diffs = 0usize;
        for i in start..end {
            let j = (i as i64 - d) as usize;
            let (x, y) = (a[i], b[j]);
            if x != b'N' && y != b'N' {
                let k = qa[i] as usize * QN + qb[j] as usize;
                if x == y {
                    score += t.agree[k];
                } else {
                    score += t.differ[k];
                    diffs += 1;
                }
            }
            if let Some(bs) = bound {
                if score + 2.0 * (end - i - 1) as f64 <= bs {
                    return None;
                }
            }
        }
        Some((score, diffs))
    };
    // Seed: the diagonal with the most shared 6-mers is scored first, so that
    // its score prunes the others (a true overlap of any length is found
    // without scoring every diagonal in full).
    let floor = opt.min_score - 1e-9;
    let mut best: Option<(f64, i64, usize)> = None; // score, d, diffs
    let mut seed_d: Option<i64> = None;
    if n1 >= 6 && n2 >= 6 {
        let code6 = |s: &[u8], i: usize| -> Option<usize> {
            let mut c = 0usize;
            for &x in &s[i..i + 6] {
                let v = CODE[x as usize];
                if v > 3 {
                    return None;
                }
                c = (c << 2) | v as usize;
            }
            Some(c)
        };
        let mut head = vec![-1i32; 4096];
        let mut next = vec![-1i32; n1];
        for i in 0..=n1 - 6 {
            if let Some(c) = code6(&a, i) {
                next[i] = head[c];
                head[c] = i as i32;
            }
        }
        let mut votes = vec![0u16; n1 + n2];
        for j in 0..=n2 - 6 {
            if let Some(c) = code6(&b, j) {
                let mut i = head[c];
                while i >= 0 {
                    votes[(i as i64 - j as i64 + l2) as usize] += 1;
                    i = next[i as usize];
                }
            }
        }
        let mut top = (1u16, 0i64);
        for (k, &v) in votes.iter().enumerate() {
            let d = k as i64 - l2;
            if v > top.0 && d >= d_lo && d <= d_hi {
                top = (v, d);
            }
        }
        if top.0 >= 2 {
            seed_d = Some(top.1);
            let (s, df) = eval(top.1, None).expect("no bound");
            best = Some((s, top.1, df));
        }
    }
    // Then every diagonal: the plateau of the longest overlap, then outwards
    // by decreasing overlap.
    let (p_lo, p_hi) = ((l1 - l2).min(0), (l1 - l2).max(0));
    let mut left = p_lo - 1;
    let mut right = p_hi + 1;
    let mut plateau = p_lo;
    loop {
        let d = if plateau <= p_hi {
            plateau += 1;
            plateau - 1
        } else {
            let ol = if left >= d_lo { overlap(left) } else { -1 };
            let or = if right <= d_hi { overlap(right) } else { -1 };
            if ol < 0 && or < 0 {
                break;
            }
            if ol >= or {
                left -= 1;
                left + 1
            } else {
                right += 1;
                right - 1
            }
        };
        let o = overlap(d);
        if o < min_ov as i64 || seed_d == Some(d) {
            continue;
        }
        // A diagonal must beat the best so far and reach min_score; on equal
        // scores the longer overlap (met first) stays.
        let bound = best.map_or(floor, |x| x.0.max(floor));
        if 2.0 * o as f64 <= bound {
            // every remaining diagonal is at most this long
            break;
        }
        if let Some((score, diffs)) = eval(d, Some(bound)) {
            best = Some((score, d, diffs));
        }
    }
    let Some((score, d, diffs)) = best else {
        return Err(MergeFail::NoAlignment);
    };
    if score < opt.min_score {
        return Err(MergeFail::NoAlignment);
    }
    let o = overlap(d) as usize;
    if diffs > opt.max_diffs || (diffs as f64) * 100.0 > opt.max_diff_pct * o as f64 {
        return Err(MergeFail::TooManyDiffs);
    }
    let staggered = d < 0 || d + l2 < l1;
    if staggered && !opt.allow_stagger {
        return Err(MergeFail::Staggered);
    }
    let len = (d + l2) as usize;
    let mut seq = Vec::with_capacity(len);
    let mut qual = Vec::with_capacity(len);
    for p in 0..len as i64 {
        let in1 = p < l1;
        let in2 = p >= d;
        if in1 && in2 {
            let (i, j) = (p as usize, (p - d) as usize);
            let (x, y, qx, qy) = (a[i], b[j], qa[i], qb[j]);
            if x == b'N' && y == b'N' {
                seq.push(b'N');
                qual.push(qx.min(qy).min(2));
            } else if x == b'N' {
                seq.push(y);
                qual.push(qy);
            } else if y == b'N' {
                seq.push(x);
                qual.push(qx);
            } else if x == y {
                seq.push(x);
                qual.push(t.q_agree[qx as usize * QN + qy as usize]);
            } else if qx >= qy {
                seq.push(x);
                qual.push(t.q_differ[qx as usize * QN + qy as usize]);
            } else {
                seq.push(y);
                qual.push(t.q_differ[qy as usize * QN + qx as usize]);
            }
        } else if in1 {
            seq.push(a[p as usize]);
            qual.push(qa[p as usize]);
        } else {
            let j = (p - d) as usize;
            seq.push(b[j]);
            qual.push(qb[j]);
        }
    }
    for q in qual.iter_mut() {
        *q = (*q).min(opt.max_q_out);
    }
    Ok(MergedRead { seq, qual, overlap: o, diffs, score, staggered })
}

/// Expected errors of a read: the sum of 10^(-Q/10) over its NUMERIC
/// qualities (as vsearch --fastq_maxee).
pub fn expected_errors(qual: &[u8]) -> f64 {
    qual.iter().map(|&q| 10f64.powf(-(q as f64) / 10.0)).sum()
}

// ---------------------------------------------------------------------
//  PRIMERS
// ---------------------------------------------------------------------

/// A primer as bit masks (IUPAC), and its reverse complement.
#[derive(Clone, Debug)]
pub struct Primer {
    pub text: String,
    fwd: Vec<u8>,
    rc: Vec<u8>,
    pub max_mismatches: usize,
}

impl Primer {
    pub fn new(text: &str, max_error_rate: f64) -> Result<Primer, String> {
        let t = text.trim();
        if t.is_empty() {
            return Err("A primer is empty: give its sequence (IUPAC letters, 5' to 3'), or leave the option out".to_string());
        }
        let mut fwd = Vec::with_capacity(t.len());
        for (i, c) in t.bytes().enumerate() {
            let m = iupac_mask(c);
            if m == 0 {
                return Err(format!(
                    "The primer {t} has '{}' at position {}, which is not an IUPAC nucleotide letter (A C G T U R Y S W K M B D H V N). Check the primer sequence",
                    char::from(c),
                    i + 1
                ));
            }
            fwd.push(m);
        }
        if !(0.0..1.0).contains(&max_error_rate) {
            return Err(format!("primer_max_error_rate = {max_error_rate}: give a fraction from 0 to below 1, for instance 0.1"));
        }
        let rc = fwd.iter().rev().map(|&m| comp_mask(m)).collect();
        let max_mismatches = (max_error_rate * fwd.len() as f64 + 1e-9).floor() as usize;
        Ok(Primer { text: t.to_ascii_uppercase(), fwd, rc, max_mismatches })
    }

    pub fn len(&self) -> usize {
        self.fwd.len()
    }
}

fn read_mask(b: u8) -> u8 {
    match b {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        _ => 0,
    }
}

fn mismatches_at(seq: &[u8], masks: &[u8], start: usize, limit: usize) -> usize {
    let mut mm = 0;
    for (k, &m) in masks.iter().enumerate() {
        if read_mask(seq[start + k]) & m == 0 {
            mm += 1;
            if mm > limit {
                break;
            }
        }
    }
    mm
}

/// The forward primer starting within the first `window` + 1 positions:
/// `Some(end of the primer)` for the best match (fewest mismatches, then
/// leftmost) with at most `max_mismatches`.
pub fn find_forward_primer(seq: &[u8], p: &Primer, window: usize) -> Option<usize> {
    let n = p.fwd.len();
    if seq.len() < n {
        return None;
    }
    let mut best: Option<(usize, usize)> = None;
    for s in 0..=window.min(seq.len() - n) {
        let mm = mismatches_at(seq, &p.fwd, s, p.max_mismatches);
        if mm <= p.max_mismatches && best.map_or(true, |(b, _)| mm < b) {
            best = Some((mm, s));
            if mm == 0 {
                break;
            }
        }
    }
    best.map(|(_, s)| s + n)
}

/// The reverse primer (as its reverse complement) ending within the last
/// `window` + 1 positions: `Some(start of it)` for the best match (fewest
/// mismatches, then rightmost).
pub fn find_reverse_primer(seq: &[u8], p: &Primer, window: usize) -> Option<usize> {
    let n = p.rc.len();
    if seq.len() < n {
        return None;
    }
    let mut best: Option<(usize, usize)> = None;
    for back in 0..=window.min(seq.len() - n) {
        let s = seq.len() - n - back;
        let mm = mismatches_at(seq, &p.rc, s, p.max_mismatches);
        if mm <= p.max_mismatches && best.map_or(true, |(b, _)| mm < b) {
            best = Some((mm, s));
            if mm == 0 {
                break;
            }
        }
    }
    best.map(|(_, s)| s)
}

// ---------------------------------------------------------------------
//  DISTANCES BETWEEN SEQUENCES
// ---------------------------------------------------------------------

/// Edit distance (substitutions and single-base indels, unit costs) when it
/// is at most `k`; `None` when it is larger. A band of 2k + 1 diagonals, and
/// the computation stops as soon as every cell of a row exceeds `k`.
pub fn bounded_edit_distance(a: &[u8], b: &[u8], k: usize) -> Option<usize> {
    let (n, m) = (a.len(), b.len());
    if n.abs_diff(m) > k {
        return None;
    }
    if k == 0 {
        return if a == b { Some(0) } else { None };
    }
    let big = (k + 1) as u32;
    let w = 2 * k + 1;
    let mut prev = vec![big; w];
    let mut cur = vec![big; w];
    for j in 0..=k.min(m) {
        prev[j + k] = j as u32;
    }
    for i in 1..=n {
        let mut row_min = big;
        for jj in 0..w {
            let j = i as i64 - k as i64 + jj as i64;
            if j < 0 || j > m as i64 {
                cur[jj] = big;
                continue;
            }
            let j = j as usize;
            let v = if j == 0 {
                (i as u32).min(big)
            } else {
                let mut v = prev[jj] + u32::from(a[i - 1] != b[j - 1]);
                if jj + 1 < w {
                    v = v.min(prev[jj + 1] + 1);
                }
                if jj > 0 {
                    v = v.min(cur[jj - 1] + 1);
                }
                v.min(big)
            };
            cur[jj] = v;
            row_min = row_min.min(v);
        }
        if row_min > k as u32 {
            return None;
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    let d = prev[m + k - n];
    if d <= k as u32 {
        Some(d as usize)
    } else {
        None
    }
}

/// Differences charged to each position of `q` by an end-to-end alignment
/// of `q` with `p` (banded; a mismatch costs 1 and a gap 3, so that
/// homologous regions stay on their diagonal and a divergent region does
/// not pull gaps into an identical one): a mismatch or an extra base of `q`
/// counts 1 at its position; a base of `p` missing from `q` counts 1 at the
/// next position of `q` (the last one at the end).
pub fn diff_profile(q: &[u8], p: &[u8]) -> Vec<u16> {
    const GAP: u32 = 3;
    let (n, m) = (q.len(), p.len());
    let mut prof = vec![0u16; n];
    if n == 0 {
        return prof;
    }
    let k = n.abs_diff(m) + 16;
    let w = 2 * k + 1;
    let big = u32::MAX / 4;
    let inb = |i: usize, j: usize| j + k >= i && j <= i + k;
    let at = |i: usize, j: usize| i * w + (j + k - i);
    let mut dp = vec![big; (n + 1) * w];
    for j in 0..=k.min(m) {
        dp[at(0, j)] = GAP * j as u32;
    }
    for i in 1..=n {
        let lo = i.saturating_sub(k);
        let hi = (i + k).min(m);
        for j in lo..=hi {
            let v = if j == 0 {
                GAP * i as u32
            } else {
                let mut v = big;
                if inb(i - 1, j - 1) {
                    v = dp[at(i - 1, j - 1)] + u32::from(q[i - 1] != p[j - 1]);
                }
                if inb(i - 1, j) {
                    v = v.min(dp[at(i - 1, j)] + GAP);
                }
                if j - 1 >= lo {
                    v = v.min(dp[at(i, j - 1)] + GAP);
                }
                v
            };
            dp[at(i, j)] = v;
        }
    }
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        let here = dp[at(i, j)];
        if i > 0 && j > 0 && inb(i - 1, j - 1) && here == dp[at(i - 1, j - 1)] + u32::from(q[i - 1] != p[j - 1]) {
            if q[i - 1] != p[j - 1] {
                prof[i - 1] += 1;
            }
            i -= 1;
            j -= 1;
        } else if i > 0 && inb(i - 1, j) && here == dp[at(i - 1, j)] + GAP {
            prof[i - 1] += 1;
            i -= 1;
        } else {
            prof[i.min(n - 1)] += 1;
            j -= 1;
        }
    }
    prof
}

// ---------------------------------------------------------------------
//  UNOISE3
// ---------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct DenoiseOptions {
    /// The alpha of the skew rule (2).
    pub alpha: f64,
    /// Uniques with fewer reads never become ASVs (8, as UNOISE3).
    pub min_size: u64,
    /// Assign the uniques below `min_size` to an ASV by the same rule.
    pub map_rare: bool,
}

impl Default for DenoiseOptions {
    fn default() -> Self {
        DenoiseOptions { alpha: 2.0, min_size: 8, map_rare: true }
    }
}

/// The largest number of differences at which a sequence with `a_q` reads
/// is an error of one with `a_c` reads: d <= (log2(a_c / a_q) - 1) / alpha.
/// Negative when not even one difference is allowed.
pub fn max_error_distance(a_c: u64, a_q: u64, alpha: f64) -> i64 {
    if a_q == 0 || a_c == 0 {
        return -1;
    }
    (((a_c as f64 / a_q as f64).log2() - 1.0) / alpha + 1e-9).floor() as i64
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Denoised {
    /// Unique indices of the centroids, from the most abundant.
    pub centroids: Vec<usize>,
    /// For each unique: the centroid (an index into `centroids`) it belongs
    /// to; None = below min_size and not explained by any centroid.
    pub assigned: Vec<Option<usize>>,
    /// Edit distance of each unique to its centroid (0 for a centroid).
    pub distance: Vec<u32>,
}

/// UNOISE3 on dereplicated sequences (`seqs[i]` seen `abund[i]` times).
pub fn unoise(seqs: &[Vec<u8>], abund: &[u64], opt: &DenoiseOptions) -> Denoised {
    let n = seqs.len();
    let alpha = opt.alpha;
    let order = sorted_indices(n, |&x, &y| abund[y as usize].cmp(&abund[x as usize]).then_with(|| seqs[x as usize].cmp(&seqs[y as usize])));
    let mut out = Denoised { centroids: Vec::new(), assigned: vec![None; n], distance: vec![0; n] };
    let min_ratio = 2f64.powf(alpha + 1.0); // a centroid must be this much more abundant for d = 1
    let n_candidates = |cents: &[usize], aq: u64| -> usize {
        let need = aq as f64 * min_ratio * (1.0 - 1e-12);
        cents.partition_point(|&c| abund[c] as f64 >= need)
    };
    let mut rare: Vec<usize> = Vec::new();
    for &q in &order {
        let q = q as usize;
        let aq = abund[q];
        if aq < opt.min_size.max(1) {
            rare.push(q);
            continue;
        }
        let nc = n_candidates(&out.centroids, aq);
        let cents = &out.centroids;
        let hit = (0..nc)
            .into_par_iter()
            .filter_map(|ci| {
                let c = cents[ci];
                let k = max_error_distance(abund[c], aq, alpha);
                if k < 1 {
                    return None;
                }
                bounded_edit_distance(&seqs[q], &seqs[c], k as usize).map(|d| (d, ci))
            })
            .min();
        match hit {
            Some((d, ci)) => {
                out.assigned[q] = Some(ci);
                out.distance[q] = d as u32;
            }
            None => {
                out.assigned[q] = Some(out.centroids.len());
                out.centroids.push(q);
            }
        }
    }
    if opt.map_rare && !out.centroids.is_empty() {
        let cents = &out.centroids;
        let found: Vec<Option<(usize, usize)>> = rare
            .par_iter()
            .map(|&q| {
                let aq = abund[q];
                let nc = n_candidates(cents, aq);
                let mut best: Option<(usize, usize)> = None;
                for (ci, &c) in cents.iter().enumerate().take(nc) {
                    let mut k = max_error_distance(abund[c], aq, alpha);
                    if let Some((bd, _)) = best {
                        k = k.min(bd as i64 - 1);
                    }
                    if k < 1 {
                        continue;
                    }
                    if let Some(d) = bounded_edit_distance(&seqs[q], &seqs[c], k as usize) {
                        best = Some((d, ci));
                        if d <= 1 {
                            break;
                        }
                    }
                }
                best
            })
            .collect();
        for (&q, f) in rare.iter().zip(found) {
            if let Some((d, ci)) = f {
                out.assigned[q] = Some(ci);
                out.distance[q] = d as u32;
            }
        }
    }
    out
}

// ---------------------------------------------------------------------
//  CHIMERAS
// ---------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct ChimeraOptions {
    pub enabled: bool,
    /// Each parent must have at least this many times the reads of the
    /// query (16, as --uchime3_denovo).
    pub abskew: f64,
    /// Most differences between the query and its two-parent model (0: a
    /// perfect chimera, as UCHIME2/3 in denoised mode).
    pub max_model_diffs: usize,
    /// Each segment must carry at least this many differences that the
    /// other parent has and the query does not (3, as --mindiffs).
    pub min_diffs: usize,
}

impl Default for ChimeraOptions {
    fn default() -> Self {
        ChimeraOptions { enabled: true, abskew: 16.0, max_model_diffs: 0, min_diffs: 3 }
    }
}

/// The two-parent model of a chimeric sequence. Indices are positions in
/// the slice given to `find_chimeras`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChimeraModel {
    /// The parent of the left part (query[..breakpoint]) and of the right.
    pub parent_left: usize,
    pub parent_right: usize,
    pub breakpoint: usize,
    pub model_diffs: usize,
    pub closest_parent: usize,
    pub closest_parent_diffs: usize,
    /// Differences of the right parent on the left side, and of the left
    /// parent on the right side, net of the model's own differences there:
    /// the evidence for each segment.
    pub votes_left: usize,
    pub votes_right: usize,
}

/// The best two-parent model of `q` among the candidates, if it makes `q` a
/// chimera.
fn chimera_model(q: &[u8], cands: &[usize], seqs: &[&[u8]], opt: &ChimeraOptions) -> Option<ChimeraModel> {
    let len = q.len();
    let prefix: Vec<Vec<u32>> = cands
        .par_iter()
        .map(|&c| {
            let prof = diff_profile(q, seqs[c]);
            let mut pre = Vec::with_capacity(len + 1);
            let mut s = 0u32;
            pre.push(0);
            for &x in &prof {
                s += u32::from(x);
                pre.push(s);
            }
            pre
        })
        .collect();
    let (mut single, mut single_c) = (u32::MAX, 0usize);
    for (ci, pre) in prefix.iter().enumerate() {
        if pre[len] < single {
            single = pre[len];
            single_c = ci;
        }
    }
    // For each crossover x: the fewest differences on each side; among the
    // parents that reach them, the left parent that differs most from Q on
    // the right and the right parent that differs most on the left (those
    // differences are the votes for the model). Best model: fewest
    // differences, then the most votes on its weaker side.
    let mut best: Option<(u32, i64, usize, usize, usize, i64, i64)> = None; // diffs, -min votes, x, left, right, votes l, votes r
    for x in 0..=len {
        let mut lm = u32::MAX;
        let mut rm = u32::MAX;
        for pre in &prefix {
            lm = lm.min(pre[x]);
            rm = rm.min(pre[len] - pre[x]);
        }
        let model = lm + rm;
        if model as usize > opt.max_model_diffs {
            continue;
        }
        let (mut la, mut a_right) = (usize::MAX, 0u32);
        let (mut ra, mut b_left) = (usize::MAX, 0u32);
        for (ci, pre) in prefix.iter().enumerate() {
            let (l, r) = (pre[x], pre[len] - pre[x]);
            if l == lm && (la == usize::MAX || r > a_right) {
                la = ci;
                a_right = r;
            }
            if r == rm && (ra == usize::MAX || l > b_left) {
                ra = ci;
                b_left = l;
            }
        }
        if la == ra {
            continue;
        }
        let votes_l = b_left as i64 - lm as i64;
        let votes_r = a_right as i64 - rm as i64;
        let key = (model, -votes_l.min(votes_r));
        if best.map_or(true, |b| key < (b.0, b.1)) {
            best = Some((model, -votes_l.min(votes_r), x, la, ra, votes_l, votes_r));
        }
    }
    let (model, _, x, la, ra, vl, vr) = best?;
    let need = opt.min_diffs.max(1) as i64;
    if vl >= need && vr >= need && single > model {
        Some(ChimeraModel {
            parent_left: cands[la],
            parent_right: cands[ra],
            breakpoint: x,
            model_diffs: model as usize,
            closest_parent: cands[single_c],
            closest_parent_diffs: single as usize,
            votes_left: vl as usize,
            votes_right: vr as usize,
        })
    } else {
        None
    }
}

/// De novo chimera check. `seqs` must be in decreasing abundance (`abund`);
/// each one is checked against the more abundant ones already judged
/// non-chimeric.
pub fn find_chimeras(seqs: &[&[u8]], abund: &[u64], opt: &ChimeraOptions) -> Vec<Option<ChimeraModel>> {
    let n = seqs.len();
    let mut out = vec![None; n];
    let mut clean: Vec<usize> = Vec::new();
    for i in 0..n {
        let need = abund[i] as f64 * opt.abskew;
        let cands: Vec<usize> = clean.iter().copied().filter(|&j| abund[j] as f64 >= need).collect();
        if cands.len() >= 2 {
            if let Some(m) = chimera_model(seqs[i], &cands, seqs, opt) {
                out[i] = Some(m);
                continue;
            }
        }
        clean.push(i);
    }
    out
}

// ---------------------------------------------------------------------
//  SAMPLES
// ---------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AmpliconSample {
    pub name: String,
    pub r1: String,
    /// None for single-end (or already merged) reads.
    pub r2: Option<String>,
}

fn join_dir(dir: &std::path::Path, f: &str) -> String {
    let p = std::path::Path::new(f);
    if p.is_absolute() {
        f.to_string()
    } else {
        dir.join(p).to_string_lossy().to_string()
    }
}

/// A sample sheet: one line per sample, tab-separated: name, read-1 file,
/// read-2 file (absent for single-end). File names are relative to the
/// sheet's folder unless absolute. Empty lines and lines starting with '#'
/// are skipped; a first line whose second column has no '.' is a header
/// (mothur's `stability.files` works as it is).
pub fn read_sample_sheet(path: &str) -> Result<Vec<AmpliconSample>, String> {
    check_input(path)?;
    let mut text = String::new();
    abre_leitura(path)?.read_to_string(&mut text).map_err(|e| format!("Error reading the sample sheet {path}: {e}"))?;
    let dir = std::path::Path::new(path).parent().map(|p| p.to_path_buf()).unwrap_or_default();
    let mut out: Vec<AmpliconSample> = Vec::new();
    let mut first = true;
    for (i, line) in text.lines().enumerate() {
        let l = line.trim();
        if l.is_empty() || l.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = if l.contains('\t') { l.split('\t').map(str::trim).filter(|s| !s.is_empty()).collect() } else { l.split_whitespace().collect() };
        let header = first && f.len() >= 2 && !f[1].contains('.');
        first = false;
        if header {
            continue;
        }
        if f.len() < 2 || f.len() > 3 {
            return Err(format!(
                "{path}, line {}: expected 2 or 3 tab-separated columns (sample name, read-1 FASTQ, read-2 FASTQ), found {}: \"{}\". Fix that line",
                i + 1,
                f.len(),
                l
            ));
        }
        let name = f[0].to_string();
        if out.iter().any(|s| s.name == name) {
            return Err(format!("{path}, line {}: the sample name {name} appears twice. Each sample needs its own name", i + 1));
        }
        let r1 = join_dir(&dir, f[1]);
        check_input(&r1).map_err(|e| format!("{path}, line {} (sample {name}): {e}", i + 1))?;
        let r2 = match f.get(2) {
            Some(x) => {
                let r2 = join_dir(&dir, x);
                check_input(&r2).map_err(|e| format!("{path}, line {} (sample {name}): {e}", i + 1))?;
                Some(r2)
            }
            None => None,
        };
        out.push(AmpliconSample { name, r1, r2 });
    }
    if out.is_empty() {
        return Err(format!("The sample sheet {path} lists no sample. Write one line per sample: name, read-1 file, read-2 file (tab-separated)"));
    }
    Ok(out)
}

/// The read-1 / read-2 marker of an Illumina-style file name: the position
/// and length of `R1`/`R2` (or `1`/`2`) between separators.
fn pair_marker(name: &str) -> Option<(usize, usize, u8)> {
    let b = name.as_bytes();
    let sep = |c: u8| matches!(c, b'_' | b'.' | b'-');
    // _R1_ _R1. .R1. -R1_ ...
    for i in 1..b.len().saturating_sub(2) {
        if sep(b[i - 1]) && (b[i] == b'R' || b[i] == b'r') && matches!(b[i + 1], b'1' | b'2') && (i + 2 == b.len() || sep(b[i + 2])) {
            return Some((i, 2, b[i + 1] - b'0'));
        }
    }
    // _1.fastq _2.fq.gz
    for i in 1..b.len().saturating_sub(1) {
        if sep(b[i - 1]) && matches!(b[i], b'1' | b'2') && i + 1 < b.len() && b[i + 1] == b'.' {
            return Some((i, 1, b[i] - b'0'));
        }
    }
    None
}

/// Builds the samples from a list of FASTQ files: `x_R1_001.fastq` pairs
/// with `x_R2_001.fastq` (also `_1`/`_2`); the sample name is the file name
/// up to its first '_' (Illumina's convention), or up to the marker when that
/// would repeat. A file without a mate is a single-end sample.
pub fn pair_fastq_files(paths: &[String]) -> Result<Vec<AmpliconSample>, String> {
    let mut groups: Vec<(String, Option<String>, Option<String>, String)> = Vec::new(); // key, r1, r2, base
    for p in paths {
        check_input(p)?;
        let base = std::path::Path::new(p).file_name().map(|x| x.to_string_lossy().to_string()).unwrap_or_else(|| p.clone());
        let (key, which, stem) = match pair_marker(&base) {
            Some((i, l, w)) => (format!("{}\u{0}{}", &base[..i], &base[i + l..]), w, base[..i].trim_end_matches(['_', '.', '-']).to_string()),
            None => (base.clone(), 1, base.split('.').next().unwrap_or(&base).to_string()),
        };
        let dir = std::path::Path::new(p).parent().map(|d| d.to_string_lossy().to_string()).unwrap_or_default();
        let key = format!("{dir}\u{0}{key}");
        let slot = match groups.iter().position(|g| g.0 == key) {
            Some(s) => s,
            None => {
                groups.push((key, None, None, stem));
                groups.len() - 1
            }
        };
        let g = &mut groups[slot];
        let place = if which == 2 { &mut g.2 } else { &mut g.1 };
        if place.is_some() {
            return Err(format!("{p}: two files claim to be read {which} of the same sample. Give each file once"));
        }
        *place = Some(p.clone());
    }
    let short: Vec<String> = groups.iter().map(|g| g.3.split('_').next().unwrap_or(&g.3).to_string()).collect();
    let mut out = Vec::new();
    for (i, g) in groups.into_iter().enumerate() {
        let unique = short.iter().filter(|s| **s == short[i]).count() == 1;
        let name = if unique { short[i].clone() } else { g.3.clone() };
        match (g.1, g.2) {
            (Some(r1), r2) => out.push(AmpliconSample { name, r1, r2 }),
            (None, Some(r2)) => {
                return Err(format!("{r2} is a read-2 file without its read-1 mate. Give both files of each pair"));
            }
            (None, None) => {}
        }
    }
    if out.iter().map(|s| &s.name).collect::<HashSet<_>>().len() != out.len() {
        return Err("Two samples got the same name from their file names. Use a sample sheet to name them".to_string());
    }
    Ok(out)
}

// ---------------------------------------------------------------------
//  THE ASV PIPELINE
// ---------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct AsvOptions {
    pub merge: MergeOptions,
    /// Primers to remove (IUPAC, 5'->3'); None = the reads have none.
    pub primer_forward: Option<String>,
    pub primer_reverse: Option<String>,
    /// Mismatches allowed in a primer, as a fraction of its length (0.1).
    pub primer_max_error_rate: f64,
    /// How far from the read end a primer may start (10 bases).
    pub primer_window: usize,
    /// Discard a read where a given primer is not found (true).
    pub require_primers: bool,
    /// Most expected errors per read (1.0); None = no limit.
    pub max_ee: Option<f64>,
    /// Length range of the reads kept (after merging and primer removal);
    /// max_length 0 = no limit.
    pub min_length: usize,
    pub max_length: usize,
    /// Most `N`s per read (0).
    pub max_ns: usize,
    pub denoise: DenoiseOptions,
    pub chimera: ChimeraOptions,
}

impl Default for AsvOptions {
    fn default() -> Self {
        AsvOptions {
            merge: MergeOptions::default(),
            primer_forward: None,
            primer_reverse: None,
            primer_max_error_rate: 0.1,
            primer_window: 10,
            require_primers: true,
            max_ee: Some(1.0),
            min_length: 1,
            max_length: 0,
            max_ns: 0,
            denoise: DenoiseOptions::default(),
            chimera: ChimeraOptions::default(),
        }
    }
}

/// What happened to the reads of one sample.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SampleTrack {
    pub sample: String,
    /// Read pairs (or reads) in the files.
    pub input: u64,
    pub merged: u64,
    pub merge_too_short: u64,
    pub merge_no_alignment: u64,
    pub merge_too_many_diffs: u64,
    pub merge_staggered: u64,
    pub primer_missing: u64,
    pub filter_length: u64,
    pub filter_n: u64,
    pub filter_expected_errors: u64,
    /// Reads that passed every filter.
    pub filtered: u64,
    /// Filtered reads explained by an ASV or a chimera.
    pub denoised: u64,
    /// Reads in the final table.
    pub nonchimeric: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Asv {
    pub name: String,
    pub sequence: String,
    /// Reads in the table, over all samples.
    pub reads: u64,
    /// Distinct sequences (the ASV and its error variants) behind it.
    pub uniques: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ChimeraCall {
    pub sequence: String,
    pub reads: u64,
    /// The ASVs whose parts make it, left then right.
    pub parent_left: String,
    pub parent_right: String,
    pub breakpoint: usize,
    pub closest_parent: String,
    pub closest_parent_diffs: usize,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AsvResult {
    pub asv_fasta: String,
    pub feature_table: String,
    pub tracking: Option<String>,
    pub samples: Vec<SampleTrack>,
    /// Distinct filtered sequences over all samples, and how many had at
    /// least min_size reads.
    pub uniques: u64,
    pub uniques_min_size: u64,
    /// Centroids found by UNOISE3, before the chimera check.
    pub centroids: u64,
    pub asvs: Vec<Asv>,
    pub chimeras: Vec<ChimeraCall>,
    /// Filtered reads that no ASV explains (rare and too far from any).
    pub reads_unassigned: u64,
    pub warnings: Vec<String>,
}

enum Outcome {
    Merge(MergeFail),
    NoPrimer,
    Length,
    WithN,
    HighEe,
    Pass(Vec<u8>),
}

struct Prepared {
    fwd: Option<Primer>,
    rev: Option<Primer>,
}

/// Merge (or take the single read), remove primers, filter.
fn process_read(s1: &[u8], q1: &[u8], second: Option<(&[u8], &[u8])>, opt: &AsvOptions, pr: &Prepared) -> Outcome {
    let (mut seq, mut qual) = match second {
        Some((s2, q2)) => match merge_pair(s1, q1, s2, q2, &opt.merge) {
            Ok(m) => (m.seq, m.qual),
            Err(e) => return Outcome::Merge(e),
        },
        None => {
            let n = trunc_point(q1, opt.merge.trunc_len_r1, opt.merge.trunc_qual).min(s1.len());
            (s1[..n].iter().map(|&b| UPPER[b as usize]).collect::<Vec<u8>>(), q1[..n].to_vec())
        }
    };
    if let Some(p) = &pr.fwd {
        match find_forward_primer(&seq, p, opt.primer_window) {
            Some(end) => {
                seq.drain(..end);
                qual.drain(..end);
            }
            None if opt.require_primers => return Outcome::NoPrimer,
            None => {}
        }
    }
    if let Some(p) = &pr.rev {
        match find_reverse_primer(&seq, p, opt.primer_window) {
            Some(start) => {
                seq.truncate(start);
                qual.truncate(start);
            }
            None if opt.require_primers => return Outcome::NoPrimer,
            None => {}
        }
    }
    if seq.len() < opt.min_length.max(1) || (opt.max_length > 0 && seq.len() > opt.max_length) {
        return Outcome::Length;
    }
    if seq.iter().filter(|&&b| b == b'N').count() > opt.max_ns {
        return Outcome::WithN;
    }
    if let Some(m) = opt.max_ee {
        if expected_errors(&qual) > m {
            return Outcome::HighEe;
        }
    }
    Outcome::Pass(seq)
}

struct Derep {
    index: HashMap<Vec<u8>, u32>,
    seqs: Vec<Vec<u8>>,
    total: Vec<u64>,
    /// Per sample: unique -> reads.
    per_sample: Vec<HashMap<u32, u32>>,
}

fn quals_or_error(rec: &SeqRec, path: &str, offset: u8) -> Result<Vec<u8>, String> {
    if rec.qual.is_empty() && !rec.seq.is_empty() {
        return Err(format!("{path}: the record at line {} has no qualities: this is a FASTA file. Give the FASTQ files, which carry the qualities the merging and the expected-error filter need", rec.line));
    }
    numeric_quals(&rec.qual, offset).map_err(|i| {
        let hint = if offset == 64 {
            "the file is Phred+33 (as every recent Illumina file): set phred_offset = 33"
        } else {
            "the file is damaged (Phred+33 qualities start at '!')"
        };
        format!(
            "{path}, record at line {}: the quality character '{}' (position {}) is below the Phred offset {offset}: {hint}",
            rec.line,
            char::from(rec.qual[i]),
            i + 1
        )
    })
}

const BATCH: usize = 16_384;

fn process_sample(si: usize, s: &AmpliconSample, opt: &AsvOptions, pr: &Prepared, dr: &mut Derep) -> Result<SampleTrack, String> {
    let mut tr = SampleTrack { sample: s.name.clone(), ..Default::default() };
    let mut in1 = SeqIn::open(&s.r1)?;
    let mut in2 = match &s.r2 {
        Some(p) => Some(SeqIn::open(p)?),
        None => None,
    };
    let off = opt.merge.phred_offset;
    loop {
        let mut batch: Vec<(Vec<u8>, Vec<u8>, Option<(Vec<u8>, Vec<u8>)>)> = Vec::with_capacity(BATCH);
        let mut end = false;
        while batch.len() < BATCH {
            let mut a = SeqRec::default();
            if !in1.next(&mut a)? {
                end = true;
                if let (Some(r2), Some(p2)) = (in2.as_mut(), s.r2.as_ref()) {
                    let mut b = SeqRec::default();
                    if r2.next(&mut b)? {
                        return Err(format!(
                            "Sample {}: {} has more reads than {} (the read-1 file ended first). The two files of a pair must hold the same reads in the same order: get them again, or re-pair them (for instance with `repair.sh` of BBTools)",
                            s.name, p2, s.r1
                        ));
                    }
                }
                break;
            }
            let qa = quals_or_error(&a, &s.r1, off)?;
            let second = match (in2.as_mut(), s.r2.as_ref()) {
                (Some(r2), Some(p2)) => {
                    let mut b = SeqRec::default();
                    if !r2.next(&mut b)? {
                        return Err(format!(
                            "Sample {}: {} has more reads than {} (the read-2 file ended at read {}). The two files of a pair must hold the same reads in the same order: get them again",
                            s.name,
                            s.r1,
                            p2,
                            tr.input + batch.len() as u64 + 1
                        ));
                    }
                    if read_id(&a.head) != read_id(&b.head) {
                        return Err(format!(
                            "Sample {}: read {} is named \"{}\" in {} (line {}) and \"{}\" in {} (line {}). The two files of a pair must hold the same reads in the same order: check that they are the two halves of one run, or re-pair them",
                            s.name,
                            tr.input + batch.len() as u64 + 1,
                            String::from_utf8_lossy(read_id(&a.head)),
                            s.r1,
                            a.line,
                            String::from_utf8_lossy(read_id(&b.head)),
                            p2,
                            b.line
                        ));
                    }
                    let qb = quals_or_error(&b, p2, off)?;
                    Some((b.seq, qb))
                }
                _ => None,
            };
            batch.push((a.seq, qa, second));
        }
        let outcomes: Vec<Outcome> = batch
            .par_iter()
            .map(|(s1, q1, sec)| process_read(s1, q1, sec.as_ref().map(|(a, b)| (a.as_slice(), b.as_slice())), opt, pr))
            .collect();
        tr.input += batch.len() as u64;
        let counts = &mut dr.per_sample[si];
        for o in outcomes {
            if !matches!(o, Outcome::Merge(_)) {
                tr.merged += 1;
            }
            match o {
                Outcome::Merge(MergeFail::TooShort) => tr.merge_too_short += 1,
                Outcome::Merge(MergeFail::NoAlignment) => tr.merge_no_alignment += 1,
                Outcome::Merge(MergeFail::TooManyDiffs) => tr.merge_too_many_diffs += 1,
                Outcome::Merge(MergeFail::Staggered) => tr.merge_staggered += 1,
                Outcome::NoPrimer => tr.primer_missing += 1,
                Outcome::Length => tr.filter_length += 1,
                Outcome::WithN => tr.filter_n += 1,
                Outcome::HighEe => tr.filter_expected_errors += 1,
                Outcome::Pass(seq) => {
                    tr.filtered += 1;
                    let id = match dr.index.get(&seq) {
                        Some(&id) => id,
                        None => {
                            let id = dr.seqs.len() as u32;
                            dr.index.insert(seq.clone(), id);
                            dr.seqs.push(seq);
                            dr.total.push(0);
                            id
                        }
                    };
                    dr.total[id as usize] += 1;
                    *counts.entry(id).or_insert(0) += 1;
                }
            }
        }
        if end {
            break;
        }
    }
    Ok(tr)
}

fn check_asv_options(opt: &AsvOptions) -> Result<(), String> {
    let m = &opt.merge;
    if m.min_overlap == 0 {
        return Err("min_overlap = 0: give at least 1 (vsearch's default is 10)".to_string());
    }
    if !(m.max_diff_pct >= 0.0) {
        return Err(format!("max_diff_pct = {}: give a percentage from 0 to 100", m.max_diff_pct));
    }
    if m.phred_offset != 33 && m.phred_offset != 64 {
        return Err(format!("phred_offset = {}: FASTQ qualities are offset by 33 (today) or 64 (Illumina 1.3-1.7)", m.phred_offset));
    }
    if !(opt.denoise.alpha > 0.0) {
        return Err(format!("alpha = {}: give a positive number (UNOISE3 uses 2)", opt.denoise.alpha));
    }
    if opt.denoise.min_size == 0 {
        return Err("min_size = 0: give at least 1 (UNOISE3 uses 8)".to_string());
    }
    if let Some(e) = opt.max_ee {
        if !(e >= 0.0) {
            return Err(format!("max_ee = {e}: give a number of expected errors >= 0 (1 is usual), or none for no limit"));
        }
    }
    if opt.max_length > 0 && opt.max_length < opt.min_length {
        return Err(format!("max_length = {} is below min_length = {}: give a range that holds the amplicon", opt.max_length, opt.min_length));
    }
    if !(opt.chimera.abskew >= 1.0) {
        return Err(format!("abskew = {}: give a number >= 1 (UCHIME3 uses 16)", opt.chimera.abskew));
    }
    Ok(())
}

/// The whole pipeline: merge, remove primers, filter, dereplicate, denoise
/// (UNOISE3), remove chimeras, and write the ASV FASTA, the feature table
/// (ASV x sample counts, TSV) and, if asked, the per-sample tracking (TSV).
pub fn amplicon_asvs(samples: &[AmpliconSample], asv_fasta: &str, feature_table: &str, tracking: Option<&str>, opt: &AsvOptions) -> Result<AsvResult, String> {
    if samples.is_empty() {
        return Err("No sample given. Give a sample sheet, or the FASTQ files of each sample".to_string());
    }
    check_asv_options(opt)?;
    let mut names = HashSet::new();
    let mut inputs: Vec<&str> = Vec::new();
    for s in samples {
        if s.name.trim().is_empty() || s.name.contains('\t') {
            return Err(format!("The sample of {} has an empty name or a name with a tab. Give each sample a name", s.r1));
        }
        if !names.insert(s.name.as_str()) {
            return Err(format!("The sample name {} appears twice. Each sample needs its own name", s.name));
        }
        check_input(&s.r1)?;
        inputs.push(&s.r1);
        if let Some(r2) = &s.r2 {
            check_input(r2)?;
            if same_file(r2, &s.r1) {
                return Err(format!("Sample {}: the read-1 and read-2 files are the same file ({r2}). Give the two halves of the pair", s.name));
            }
            inputs.push(r2);
        }
    }
    let mut outs = vec![asv_fasta, feature_table];
    if let Some(t) = tracking {
        outs.push(t);
    }
    for o in &outs {
        not_an_input(o, &inputs)?;
    }
    for (i, a) in outs.iter().enumerate() {
        for b in &outs[i + 1..] {
            if same_file(a, b) {
                return Err(format!("Two outputs have the same name ({a}). Give the ASV FASTA, the feature table and the tracking table different names"));
            }
        }
    }
    let pr = Prepared {
        fwd: opt.primer_forward.as_deref().map(|p| Primer::new(p, opt.primer_max_error_rate)).transpose()?,
        rev: opt.primer_reverse.as_deref().map(|p| Primer::new(p, opt.primer_max_error_rate)).transpose()?,
    };
    let mut dr = Derep { index: HashMap::new(), seqs: Vec::new(), total: Vec::new(), per_sample: vec![HashMap::new(); samples.len()] };
    let mut tracks = Vec::with_capacity(samples.len());
    for (si, s) in samples.iter().enumerate() {
        tracks.push(process_sample(si, s, opt, &pr, &mut dr)?);
    }
    drop(std::mem::take(&mut dr.index));
    let mut warnings = Vec::new();
    for t in &tracks {
        if t.input == 0 {
            warnings.push(format!("Sample {}: the files hold no reads", t.sample));
        } else if t.merged * 2 < t.input && samples.iter().any(|s| s.name == t.sample && s.r2.is_some()) {
            warnings.push(format!(
                "Sample {}: only {} of the pairs merged. If the reads do not reach each other (amplicon longer than read 1 + read 2 - {}), merging cannot work; if their ends are of low quality, try trunc_qual or a larger max_diffs",
                t.sample,
                pct(t.merged, t.input),
                opt.merge.min_overlap
            ));
        }
        if t.primer_missing * 2 > t.merged && t.merged > 0 {
            warnings.push(format!(
                "Sample {}: the primers were missing in {} of the merged reads. If the reads were sequenced from the primers (no primer in the reads, as in the Kozich et al. protocol), leave the primers out",
                t.sample,
                pct(t.primer_missing, t.merged)
            ));
        }
    }
    // Denoise
    let den = unoise(&dr.seqs, &dr.total, &opt.denoise);
    let nc = den.centroids.len();
    // reads of each centroid: all assigned (the table), and those of the
    // uniques with >= min_size reads (UNOISE3's cluster size, which is what
    // the chimera check compares, as --uchime3_denovo after --cluster_unoise)
    let mut cent_reads = vec![0u64; nc];
    let mut cent_size = vec![0u64; nc];
    let mut cent_uniques = vec![0u64; nc];
    for (u, a) in den.assigned.iter().enumerate() {
        if let Some(ci) = a {
            cent_reads[*ci] += dr.total[u];
            cent_uniques[*ci] += 1;
            if dr.total[u] >= opt.denoise.min_size {
                cent_size[*ci] += dr.total[u];
            }
        }
    }
    let cseq = |c: usize| dr.seqs[den.centroids[c]].as_slice();
    let by_size = sorted_indices(nc, |&x, &y| cent_size[y as usize].cmp(&cent_size[x as usize]).then_with(|| cseq(x as usize).cmp(cseq(y as usize))));
    let size_seqs: Vec<&[u8]> = by_size.iter().map(|&c| cseq(c as usize)).collect();
    let size_abund: Vec<u64> = by_size.iter().map(|&c| cent_size[c as usize]).collect();
    let chim_ranked = if opt.chimera.enabled { find_chimeras(&size_seqs, &size_abund, &opt.chimera) } else { vec![None; nc] };
    // per centroid, with the parents as centroid indices
    let mut chim: Vec<Option<ChimeraModel>> = vec![None; nc];
    for (rank, m) in chim_ranked.into_iter().enumerate() {
        if let Some(mut m) = m {
            m.parent_left = by_size[m.parent_left] as usize;
            m.parent_right = by_size[m.parent_right] as usize;
            m.closest_parent = by_size[m.closest_parent] as usize;
            chim[by_size[rank] as usize] = Some(m);
        }
    }
    // names: the non-chimeric, by the reads they will have in the table
    let by_reads = sorted_indices(nc, |&x, &y| cent_reads[y as usize].cmp(&cent_reads[x as usize]).then_with(|| cseq(x as usize).cmp(cseq(y as usize))));
    let mut row_of_centroid: Vec<Option<usize>> = vec![None; nc];
    let mut asvs: Vec<Asv> = Vec::new();
    let mut name_of: Vec<String> = vec![String::new(); nc];
    for &c in &by_reads {
        let c = c as usize;
        if chim[c].is_none() {
            let name = format!("ASV_{}", asvs.len() + 1);
            row_of_centroid[c] = Some(asvs.len());
            name_of[c] = name.clone();
            asvs.push(Asv { name, sequence: String::from_utf8_lossy(cseq(c)).to_string(), reads: 0, uniques: cent_uniques[c] });
        }
    }
    let mut chimeras = Vec::new();
    for &c in &by_reads {
        let c = c as usize;
        if let Some(m) = &chim[c] {
            chimeras.push(ChimeraCall {
                sequence: String::from_utf8_lossy(cseq(c)).to_string(),
                reads: cent_reads[c],
                parent_left: name_of[m.parent_left].clone(),
                parent_right: name_of[m.parent_right].clone(),
                breakpoint: m.breakpoint,
                closest_parent: name_of[m.closest_parent].clone(),
                closest_parent_diffs: m.closest_parent_diffs,
            });
        }
    }
    // the table
    let ns = samples.len();
    let mut table = vec![vec![0u64; ns]; asvs.len()];
    let mut unassigned = 0u64;
    for (si, counts) in dr.per_sample.iter().enumerate() {
        for (&u, &c) in counts {
            match den.assigned[u as usize] {
                Some(ci) => {
                    tracks[si].denoised += u64::from(c);
                    if let Some(row) = row_of_centroid[ci] {
                        table[row][si] += u64::from(c);
                        tracks[si].nonchimeric += u64::from(c);
                    }
                }
                None => unassigned += u64::from(c),
            }
        }
    }
    for (row, a) in asvs.iter_mut().enumerate() {
        a.reads = table[row].iter().sum();
    }
    if asvs.is_empty() {
        warnings.push(format!(
            "No ASV: no sequence reached {} reads after filtering. Check the tracking table: if few reads merged or passed the filters, relax them; with few reads, lower min_size",
            opt.denoise.min_size
        ));
    }
    // write
    let mut cleanup = Cleanup(Vec::new());
    let mut fa = String::with_capacity(asvs.len() * 300);
    for a in &asvs {
        fa.push('>');
        fa.push_str(&a.name);
        fa.push('\n');
        fa.push_str(&a.sequence);
        fa.push('\n');
    }
    cleanup.0.push(asv_fasta.to_string());
    write_text(asv_fasta, &fa)?;
    let mut tt = String::from("ASV");
    for s in samples {
        tt.push('\t');
        tt.push_str(&s.name);
    }
    tt.push('\n');
    for (row, a) in asvs.iter().enumerate() {
        tt.push_str(&a.name);
        for c in &table[row] {
            tt.push('\t');
            tt.push_str(&c.to_string());
        }
        tt.push('\n');
    }
    cleanup.0.push(feature_table.to_string());
    write_text(feature_table, &tt)?;
    if let Some(tp) = tracking {
        let mut t = String::from("sample\tinput\tmerged\tfiltered\tdenoised\tnonchimeric\tkept_pct\tmerge_too_short\tmerge_no_alignment\tmerge_too_many_diffs\tmerge_staggered\tprimer_missing\tfilter_length\tfilter_n\tfilter_expected_errors\n");
        for x in &tracks {
            t.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                x.sample,
                x.input,
                x.merged,
                x.filtered,
                x.denoised,
                x.nonchimeric,
                if x.input > 0 { format!("{:.1}", 100.0 * x.nonchimeric as f64 / x.input as f64) } else { "NA".to_string() },
                x.merge_too_short,
                x.merge_no_alignment,
                x.merge_too_many_diffs,
                x.merge_staggered,
                x.primer_missing,
                x.filter_length,
                x.filter_n,
                x.filter_expected_errors
            ));
        }
        cleanup.0.push(tp.to_string());
        write_text(tp, &t)?;
    }
    cleanup.disarm();
    Ok(AsvResult {
        asv_fasta: asv_fasta.to_string(),
        feature_table: feature_table.to_string(),
        tracking: tracking.map(str::to_string),
        samples: tracks,
        uniques: dr.seqs.len() as u64,
        uniques_min_size: dr.total.iter().filter(|&&t| t >= opt.denoise.min_size).count() as u64,
        centroids: nc as u64,
        asvs,
        chimeras,
        reads_unassigned: unassigned,
        warnings,
    })
}

impl SampleTrack {
    fn json(&self) -> String {
        format!(
            "{{\"sample\":{},\"input\":{},\"merged\":{},\"filtered\":{},\"denoised\":{},\"nonchimeric\":{},\"merge_too_short\":{},\"merge_no_alignment\":{},\"merge_too_many_diffs\":{},\"merge_staggered\":{},\"primer_missing\":{},\"filter_length\":{},\"filter_n\":{},\"filter_expected_errors\":{}}}",
            js(&self.sample),
            self.input,
            self.merged,
            self.filtered,
            self.denoised,
            self.nonchimeric,
            self.merge_too_short,
            self.merge_no_alignment,
            self.merge_too_many_diffs,
            self.merge_staggered,
            self.primer_missing,
            self.filter_length,
            self.filter_n,
            self.filter_expected_errors
        )
    }
}

impl AsvResult {
    pub fn json(&self) -> String {
        let s: Vec<String> = self.samples.iter().map(SampleTrack::json).collect();
        let a: Vec<String> = self
            .asvs
            .iter()
            .map(|a| format!("{{\"name\":{},\"sequence\":{},\"reads\":{},\"uniques\":{}}}", js(&a.name), js(&a.sequence), a.reads, a.uniques))
            .collect();
        let c: Vec<String> = self
            .chimeras
            .iter()
            .map(|c| {
                format!(
                    "{{\"sequence\":{},\"reads\":{},\"parent_left\":{},\"parent_right\":{},\"breakpoint\":{},\"closest_parent\":{},\"closest_parent_diffs\":{}}}",
                    js(&c.sequence),
                    c.reads,
                    js(&c.parent_left),
                    js(&c.parent_right),
                    c.breakpoint,
                    js(&c.closest_parent),
                    c.closest_parent_diffs
                )
            })
            .collect();
        format!(
            "{{\"asv_fasta\":{},\"feature_table\":{},\"tracking\":{},\"samples\":[{}],\"uniques\":{},\"uniques_min_size\":{},\"centroids\":{},\"asvs\":[{}],\"chimeras\":[{}],\"reads_unassigned\":{},\"warnings\":{}}}",
            js(&self.asv_fasta),
            js(&self.feature_table),
            self.tracking.as_deref().map_or("null".to_string(), js),
            s.join(","),
            self.uniques,
            self.uniques_min_size,
            self.centroids,
            a.join(","),
            c.join(","),
            self.reads_unassigned,
            jlist(&self.warnings)
        )
    }

    pub fn resumo(&self) -> String {
        let tot = |f: fn(&SampleTrack) -> u64| self.samples.iter().map(f).sum::<u64>();
        let input = tot(|t| t.input);
        let mut o = format!(
            "{} ASVs from {} samples ({} read pairs). {} distinct sequences after filtering, {} with enough reads; {} centroids, {} chimeras removed.\n",
            self.asvs.len(),
            self.samples.len(),
            thousands(input),
            thousands(self.uniques),
            thousands(self.uniques_min_size),
            self.centroids,
            self.chimeras.len()
        );
        o.push_str(&format!(
            "Reads kept: merged {} ({}), filtered {} ({}), denoised {} ({}), non-chimeric {} ({}).\n",
            thousands(tot(|t| t.merged)),
            pct(tot(|t| t.merged), input),
            thousands(tot(|t| t.filtered)),
            pct(tot(|t| t.filtered), input),
            thousands(tot(|t| t.denoised)),
            pct(tot(|t| t.denoised), input),
            thousands(tot(|t| t.nonchimeric)),
            pct(tot(|t| t.nonchimeric), input)
        ));
        o.push_str(&format!("  {:<16} {:>9} {:>9} {:>9} {:>9} {:>11}\n", "sample", "input", "merged", "filtered", "denoised", "nonchimeric"));
        for t in &self.samples {
            o.push_str(&format!("  {:<16} {:>9} {:>9} {:>9} {:>9} {:>11}\n", t.sample, t.input, t.merged, t.filtered, t.denoised, t.nonchimeric));
        }
        for w in &self.warnings {
            o.push_str(&format!("warning: {w}\n"));
        }
        o.push_str(&format!("Written: {} and {}", self.asv_fasta, self.feature_table));
        if let Some(t) = &self.tracking {
            o.push_str(&format!(" and {t}"));
        }
        o.push('\n');
        o
    }
}

// ---------------------------------------------------------------------
//  TAXONOMY: THE NAIVE BAYESIAN CLASSIFIER (Wang et al. 2007)
// ---------------------------------------------------------------------

/// Rank names given to the levels of a lineage when the reference has none
/// (DADA2's files are Kingdom;Phylum;Class;Order;Family;Genus[;Species]).
pub const DEFAULT_RANKS: [&str; 7] = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"];

#[derive(Clone, Debug)]
pub struct TrainOptions {
    /// Word length (8, as the RDP classifier); 4 to 10.
    pub k: usize,
    /// Keep only this many levels of each lineage (6 = down to the genus);
    /// None = as written.
    pub max_depth: Option<usize>,
    /// Names of the levels; None = Kingdom, Phylum, ... Species.
    pub rank_names: Option<Vec<String>>,
    /// Train only on the region between these primers (as QIIME 2's
    /// feature-classifier extract-reads): a smaller model and a better fit
    /// to reads of that amplicon. Reference sequences where the two primer
    /// sites are not found (15% mismatches allowed) are left out.
    pub primer_forward: Option<String>,
    pub primer_reverse: Option<String>,
}

impl Default for TrainOptions {
    fn default() -> Self {
        TrainOptions { k: 8, max_depth: None, rank_names: None, primer_forward: None, primer_reverse: None }
    }
}

/// The region of a reference sequence between two primers (exclusive).
pub fn amplicon_region<'a>(seq: &'a [u8], fwd: &Primer, rev: &Primer) -> Option<&'a [u8]> {
    let up: Vec<u8> = seq.iter().map(|&b| UPPER[b as usize]).collect();
    let a = find_forward_primer(&up, fwd, up.len().saturating_sub(fwd.len()))?;
    let b = find_reverse_primer(&up[a..], rev, up.len().saturating_sub(a + rev.len()))?;
    if b == 0 {
        return None;
    }
    Some(&seq[a..a + b])
}

/// The lineage of a reference header: DADA2's ">Kingdom;Phylum;...;Genus;"
/// (the whole header), or ">ID lineage" (GTDB, SILVA/QIIME), where the
/// lineage stops at a " [" (GTDB's tags). Rank prefixes such as `g__` are
/// removed; the lineage stops at the first empty level.
pub fn parse_lineage(head: &[u8]) -> Option<Vec<String>> {
    let h = String::from_utf8_lossy(head);
    let h = h.trim();
    let mut parts = h.splitn(2, char::is_whitespace);
    let first = parts.next().unwrap_or("");
    let text: &str = if first.contains(';') {
        h
    } else {
        let rest = parts.next().unwrap_or("").trim();
        match rest.find(" [") {
            Some(p) => &rest[..p],
            None => rest,
        }
    };
    if !text.contains(';') {
        return None;
    }
    let mut out = Vec::new();
    for f in text.split(';') {
        let f = f.trim();
        let b = f.as_bytes();
        let f = if b.len() > 3 && b[0].is_ascii_alphabetic() && &b[1..3] == b"__" { &f[3..] } else { f };
        if f.is_empty() {
            break;
        }
        out.push(f.to_string());
    }
    if out.is_empty() {
        None
    } else {
        Some(out)
    }
}

/// The distinct valid k-mers of `seq` (as 2-bit codes), in order of first
/// appearance. `mark` is a scratch bitset of 4^k bits, left clean.
fn words_of(seq: &[u8], k: usize, mark: &mut [u64], out: &mut Vec<u32>) {
    out.clear();
    let mask = (1u64 << (2 * k)) - 1;
    let (mut code, mut valid) = (0u64, 0usize);
    for &b in seq {
        let c = CODE[b as usize];
        if c > 3 {
            valid = 0;
            code = 0;
            continue;
        }
        code = ((code << 2) | u64::from(c)) & mask;
        valid += 1;
        if valid >= k {
            let w = code as usize;
            if mark[w >> 6] & (1u64 << (w & 63)) == 0 {
                mark[w >> 6] |= 1u64 << (w & 63);
                out.push(w as u32);
            }
        }
    }
    for &w in out.iter() {
        mark[(w as usize) >> 6] &= !(1u64 << (w & 63));
    }
}

/// A trained classifier.
#[derive(Clone, Debug, Default)]
pub struct TaxModel {
    pub k: usize,
    pub ranks: Vec<String>,
    /// The lineage of each class ("genus": each distinct full lineage).
    pub lineages: Vec<Vec<String>>,
    /// M_G: training sequences of each class.
    pub class_size: Vec<u32>,
    /// N: training sequences used.
    pub n_sequences: u32,
    /// n(w): training sequences holding each word.
    pub word_seqs: Vec<u32>,
    /// Records of the reference left out (no lineage, no valid word, or no
    /// primer sites when trained on a region).
    pub skipped: u64,
    /// Of those, the ones without the primer sites (known only right after
    /// training; 0 for a loaded model).
    pub without_region: u64,
    /// The reference file the model was trained from (its name only).
    pub source: String,
    post_start: Vec<u64>,
    post_class: Vec<u32>,
    post_count: Vec<u32>,
    bonus: Vec<f32>,
    pen: Vec<f64>,
    prefix_id: Vec<Vec<u32>>,
}

/// Word counts of one class while training: compacted pairs plus a raw tail.
#[derive(Default)]
struct ClassAcc {
    pairs: Vec<(u32, u32)>,
    raw: Vec<u32>,
}

fn compact(acc: &mut ClassAcc, counter: &mut [u32], touched: &mut Vec<u32>) {
    touched.clear();
    for &(w, c) in &acc.pairs {
        counter[w as usize] = c;
        touched.push(w);
    }
    for &w in &acc.raw {
        if counter[w as usize] == 0 {
            touched.push(w);
        }
        counter[w as usize] += 1;
    }
    acc.pairs = touched.iter().map(|&w| (w, counter[w as usize])).collect();
    for &w in touched.iter() {
        counter[w as usize] = 0;
    }
    acc.raw = Vec::new();
}

fn default_ranks(depth: usize) -> Vec<String> {
    (0..depth).map(|i| DEFAULT_RANKS.get(i).map_or_else(|| format!("Rank{}", i + 1), |s| s.to_string())).collect()
}

impl TaxModel {
    /// Trains from a reference FASTA whose headers carry the lineage.
    pub fn train(reference: &str, opt: &TrainOptions) -> Result<TaxModel, String> {
        if !(4..=10).contains(&opt.k) {
            return Err(format!("k = {}: give a word length from 4 to 10 (the RDP classifier uses 8)", opt.k));
        }
        if opt.max_depth == Some(0) {
            return Err("max_depth = 0: give the number of levels to keep (6 = down to the genus), or leave it out".to_string());
        }
        let k = opt.k;
        let nw = 1usize << (2 * k);
        let region = match (&opt.primer_forward, &opt.primer_reverse) {
            (Some(f), Some(r)) => Some((Primer::new(f, 0.15)?, Primer::new(r, 0.15)?)),
            (None, None) => None,
            _ => return Err("To train on an amplicon region give both primers (primer_forward and primer_reverse), or neither".to_string()),
        };
        let mut no_region = 0u64;
        let mut r = SeqIn::open(reference)?;
        let mut rec = SeqRec::default();
        let mut class_of: HashMap<String, u32> = HashMap::new();
        let mut lineages: Vec<Vec<String>> = Vec::new();
        let mut sizes: Vec<u32> = Vec::new();
        let mut accs: Vec<ClassAcc> = Vec::new();
        let mut word_seqs = vec![0u32; nw];
        let mut counter = vec![0u32; nw];
        let mut touched: Vec<u32> = Vec::new();
        let mut mark = vec![0u64; (nw + 63) / 64];
        let mut words: Vec<u32> = Vec::new();
        let (mut n, mut records, mut skipped) = (0u32, 0u64, 0u64);
        let mut first_bad: Option<(u64, String)> = None;
        while r.next(&mut rec)? {
            records += 1;
            let Some(mut lin) = parse_lineage(&rec.head) else {
                skipped += 1;
                if first_bad.is_none() {
                    first_bad = Some((rec.line, String::from_utf8_lossy(&rec.head).chars().take(60).collect()));
                }
                continue;
            };
            if let Some(d) = opt.max_depth {
                lin.truncate(d);
            }
            let piece: &[u8] = match &region {
                Some((f, rv)) => match amplicon_region(&rec.seq, f, rv) {
                    Some(p) => p,
                    None => {
                        no_region += 1;
                        skipped += 1;
                        continue;
                    }
                },
                None => &rec.seq,
            };
            words_of(piece, k, &mut mark, &mut words);
            if words.is_empty() {
                skipped += 1;
                continue;
            }
            let key = lin.join("\u{1f}");
            let c = match class_of.get(&key) {
                Some(&c) => c,
                None => {
                    let c = lineages.len() as u32;
                    class_of.insert(key, c);
                    lineages.push(lin);
                    sizes.push(0);
                    accs.push(ClassAcc::default());
                    c
                }
            };
            n += 1;
            sizes[c as usize] += 1;
            let acc = &mut accs[c as usize];
            for &w in &words {
                word_seqs[w as usize] += 1;
                acc.raw.push(w);
            }
            if acc.raw.len() >= 1 << 16 {
                compact(acc, &mut counter, &mut touched);
            }
        }
        if n == 0 && no_region > 0 && first_bad.is_none() {
            return Err(format!(
                "{reference}: the primers were found in none of the {records} reference sequences, so nothing is left to train on. Check the primer sequences (5' to 3', the reverse primer as ordered), or train on the whole sequences"
            ));
        }
        if n == 0 {
            let hint = match &first_bad {
                Some((line, h)) => format!(" The header at line {line} is \">{h}\", which carries no lineage."),
                None => String::new(),
            };
            return Err(format!(
                "{reference}: no usable training sequence ({records} records).{hint} The reference must be a FASTA whose headers carry the lineage, as DADA2's training files (\">Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;\") or \">ID d__Bacteria;p__...;g__...\" (GTDB)"
            ));
        }
        if skipped - no_region > records / 2 {
            return Err(format!(
                "{reference}: {skipped} of the {records} records have no lineage in the header (or no valid {k}-mer). Is this the reference training file? The headers must carry the lineage, for instance \">Bacteria;Firmicutes;...;Lactobacillus;\""
            ));
        }
        for acc in accs.iter_mut() {
            compact(acc, &mut counter, &mut touched);
        }
        let mut npost = vec![0u64; nw + 1];
        for acc in &accs {
            for &(w, _) in &acc.pairs {
                npost[w as usize + 1] += 1;
            }
        }
        for i in 1..=nw {
            npost[i] += npost[i - 1];
        }
        let total = npost[nw] as usize;
        let mut post_class = vec![0u32; total];
        let mut post_count = vec![0u32; total];
        let mut cursor: Vec<u64> = npost[..nw].to_vec();
        for (c, acc) in accs.iter().enumerate() {
            for &(w, cnt) in &acc.pairs {
                let p = cursor[w as usize] as usize;
                post_class[p] = c as u32;
                post_count[p] = cnt;
                cursor[w as usize] += 1;
            }
        }
        let depth = lineages.iter().map(Vec::len).max().unwrap_or(0);
        let ranks = match &opt.rank_names {
            Some(rn) => {
                if rn.len() < depth {
                    return Err(format!("rank_names gives {} names, but the lineages have up to {depth} levels. Give a name for every level", rn.len()));
                }
                rn[..depth].to_vec()
            }
            None => default_ranks(depth),
        };
        let mut source = std::path::Path::new(reference).file_name().map(|x| x.to_string_lossy().to_string()).unwrap_or_default();
        if let Some((f, r)) = &region {
            source.push_str(&format!(" (region between {} and {}; {} sequences without it left out)", f.text, r.text, no_region));
        }
        let mut m = TaxModel {
            k,
            ranks,
            lineages,
            class_size: sizes,
            n_sequences: n,
            word_seqs,
            skipped,
            without_region: no_region,
            source,
            post_start: npost,
            post_class,
            post_count,
            ..Default::default()
        };
        m.finish();
        Ok(m)
    }

    /// The derived tables: the bonus of each stored pair, the penalty of
    /// each class, the taxon ids of each lineage prefix.
    fn finish(&mut self) {
        let nn = self.n_sequences as f64;
        let nw = self.word_seqs.len();
        let mut bonus = vec![0f32; self.post_class.len()];
        for w in 0..nw {
            let pw = (self.word_seqs[w] as f64 + 0.5) / (nn + 1.0);
            let lp = pw.ln();
            for p in self.post_start[w] as usize..self.post_start[w + 1] as usize {
                bonus[p] = ((self.post_count[p] as f64 + pw).ln() - lp) as f32;
            }
        }
        self.bonus = bonus;
        self.pen = self.class_size.iter().map(|&m| (m as f64 + 1.0).ln()).collect();
        let mut ids: HashMap<String, u32> = HashMap::new();
        self.prefix_id = self
            .lineages
            .iter()
            .map(|lin| {
                let mut key = String::new();
                lin.iter()
                    .map(|name| {
                        key.push_str(name);
                        key.push('\u{1f}');
                        let next = ids.len() as u32;
                        *ids.entry(key.clone()).or_insert(next)
                    })
                    .collect()
            })
            .collect();
    }

    pub fn classes(&self) -> usize {
        self.lineages.len()
    }

    /// Number of stored (word, class) pairs.
    pub fn pairs(&self) -> usize {
        self.post_class.len()
    }

    /// ln P(w | class) = ln((m(w,G) + P_w) / (M_G + 1)), from the stored
    /// tables (for checks).
    pub fn log_prob(&self, word: u32, class: usize) -> f64 {
        let w = word as usize;
        let pw = (self.word_seqs[w] as f64 + 0.5) / (self.n_sequences as f64 + 1.0);
        let (a, b) = (self.post_start[w] as usize, self.post_start[w + 1] as usize);
        let m = self.post_class[a..b].binary_search(&(class as u32)).map(|i| self.post_count[a + i]).unwrap_or(0);
        ((m as f64 + pw) / (self.class_size[class] as f64 + 1.0)).ln()
    }

    /// The score of every class for a list of words, up to a constant that
    /// is the same for all classes: sum ln P(w|G) - sum ln P_w.
    fn scores(&self, words: &[u32], acc: &mut Vec<f64>) {
        let n = words.len() as f64;
        acc.clear();
        acc.extend(self.pen.iter().map(|&p| -n * p));
        for &w in words {
            let (a, b) = (self.post_start[w as usize] as usize, self.post_start[w as usize + 1] as usize);
            for p in a..b {
                acc[self.post_class[p] as usize] += f64::from(self.bonus[p]);
            }
        }
    }

    /// The scores of every class, exactly as defined (for checks):
    /// sum over the words of ln P(w|G).
    pub fn full_scores(&self, seq: &[u8]) -> Vec<f64> {
        let nw = 1usize << (2 * self.k);
        let mut mark = vec![0u64; (nw + 63) / 64];
        let mut words = Vec::new();
        words_of(seq, self.k, &mut mark, &mut words);
        let mut acc = Vec::new();
        self.scores(&words, &mut acc);
        let base: f64 = words
            .iter()
            .map(|&w| ((self.word_seqs[w as usize] as f64 + 0.5) / (self.n_sequences as f64 + 1.0)).ln())
            .sum();
        acc.iter().map(|x| x + base).collect()
    }

    /// The class with the highest score; ties broken at random.
    fn arg_best(acc: &[f64], rng: &mut Rng) -> (usize, f64) {
        let mut best = f64::NEG_INFINITY;
        let mut chosen = 0usize;
        let mut ties = 0u64;
        for (c, &v) in acc.iter().enumerate() {
            if v > best {
                best = v;
                chosen = c;
                ties = 1;
            } else if v == best {
                ties += 1;
                if rng.below(ties) == 0 {
                    chosen = c;
                }
            }
        }
        (chosen, best)
    }

    /// Classifies one sequence (the `id` of the result is left empty).
    pub fn classify(&self, seq: &[u8], opt: &ClassifyOptions) -> TaxonAssignment {
        let nw = 1usize << (2 * self.k);
        let mut mark = vec![0u64; (nw + 63) / 64];
        let mut words = Vec::new();
        words_of(seq, self.k, &mut mark, &mut words);
        let mut rng = Rng::new(seed_for(opt.seed, seq));
        let mut acc = Vec::new();
        let mut rc = false;
        if words.is_empty() || self.lineages.is_empty() {
            return TaxonAssignment { words: 0, ..Default::default() };
        }
        self.scores(&words, &mut acc);
        let (mut best_c, best_s) = Self::arg_best(&acc, &mut rng);
        if opt.try_reverse_complement {
            let r = revcomp(seq);
            let mut w2 = Vec::new();
            words_of(&r, self.k, &mut mark, &mut w2);
            let mut acc2 = Vec::new();
            self.scores(&w2, &mut acc2);
            // compare on the full log-likelihood: add each side's constant
            let lp = |ws: &[u32]| -> f64 { ws.iter().map(|&w| ((self.word_seqs[w as usize] as f64 + 0.5) / (self.n_sequences as f64 + 1.0)).ln()).sum() };
            let (c2, s2) = Self::arg_best(&acc2, &mut rng);
            if s2 + lp(&w2) > best_s + lp(&words) {
                best_c = c2;
                words = w2;
                rc = true;
            }
        }
        let best = &self.lineages[best_c];
        let ids = &self.prefix_id[best_c];
        let boots = opt.bootstraps.max(1);
        let take = (words.len() / 8).max(1);
        let mut agree = vec![0u32; best.len()];
        let mut sub = Vec::with_capacity(take);
        for _ in 0..boots {
            sub.clear();
            for _ in 0..take {
                sub.push(words[rng.below(words.len() as u64) as usize]);
            }
            self.scores(&sub, &mut acc);
            let (w, _) = Self::arg_best(&acc, &mut rng);
            let wid = &self.prefix_id[w];
            for r in 0..best.len().min(wid.len()) {
                if wid[r] == ids[r] {
                    agree[r] += 1;
                } else {
                    break;
                }
            }
        }
        let confidence: Vec<f64> = agree.iter().map(|&a| a as f64 / boots as f64).collect();
        let depth = confidence.iter().take_while(|&&c| c + 1e-12 >= opt.min_confidence).count();
        TaxonAssignment {
            id: String::new(),
            lineage: best[..depth].to_vec(),
            best: best.clone(),
            confidence,
            species: None,
            reverse_complemented: rc,
            words: words.len(),
        }
    }

    /// Saves the model (gzip; see the format at the top of this file).
    /// Returns the size of the file in bytes.
    pub fn save(&self, path: &str) -> Result<u64, String> {
        let mut b: Vec<u8> = Vec::with_capacity(self.pairs() * 3 + (1 << 20));
        b.extend_from_slice(MODEL_MAGIC);
        b.push(self.k as u8);
        put_varint(&mut b, self.ranks.len() as u64);
        for r in &self.ranks {
            put_str(&mut b, r);
        }
        put_varint(&mut b, u64::from(self.n_sequences));
        put_varint(&mut b, self.skipped);
        put_varint(&mut b, self.lineages.len() as u64);
        for (lin, &m) in self.lineages.iter().zip(&self.class_size) {
            put_varint(&mut b, u64::from(m));
            put_varint(&mut b, lin.len() as u64);
            for name in lin {
                put_str(&mut b, name);
            }
        }
        put_str(&mut b, &self.source);
        for w in 0..self.word_seqs.len() {
            put_varint(&mut b, u64::from(self.word_seqs[w]));
            let (a, e) = (self.post_start[w] as usize, self.post_start[w + 1] as usize);
            put_varint(&mut b, (e - a) as u64);
            let mut prev = 0u32;
            for p in a..e {
                put_varint(&mut b, u64::from(self.post_class[p] - prev));
                prev = self.post_class[p];
                put_varint(&mut b, u64::from(self.post_count[p]));
            }
        }
        b.extend_from_slice(MODEL_END);
        let h = fnv1a(FNV_OFFSET, &b);
        b.extend_from_slice(&h.to_le_bytes());
        let f = File::create(path).map_err(|e| format!("Error creating the model file {path}: {e}. Check that the output folder exists and that you can write to it"))?;
        let mut gz = flate2::write::GzEncoder::new(std::io::BufWriter::with_capacity(1 << 20, f), flate2::Compression::new(6));
        gz.write_all(&b).map_err(|e| write_error(path, e))?;
        let mut w = gz.finish().map_err(|e| write_error(path, e))?;
        w.flush().map_err(|e| write_error(path, e))?;
        drop(w);
        Ok(std::fs::metadata(path).map(|m| m.len()).unwrap_or(0))
    }

    /// Loads a model saved by `save`.
    pub fn load(path: &str) -> Result<TaxModel, String> {
        check_input(path)?;
        let f = File::open(path).map_err(|e| format!("Error opening {path}: {e}"))?;
        let mut b = Vec::new();
        flate2::read::GzDecoder::new(std::io::BufReader::new(f)).read_to_end(&mut b).map_err(|e| {
            format!("{path} is not a classifier model file ({e}). Give a model saved by train_classifier, or the reference FASTA to train from")
        })?;
        let bad = |what: &str| format!("{path}: damaged classifier model ({what}). Train it again from the reference FASTA");
        if b.len() < 24 || &b[..8] != MODEL_MAGIC {
            return Err(format!("{path} is not a classifier model file (it does not start with MMGTAX01). Give a model saved by train_classifier, or the reference FASTA to train from"));
        }
        let end = b.len() - 16;
        if &b[end..end + 8] != MODEL_END {
            return Err(bad("the end marker is missing: the file is incomplete"));
        }
        let h = u64::from_le_bytes(b[end + 8..].try_into().unwrap());
        if fnv1a(FNV_OFFSET, &b[..end + 8]) != h {
            return Err(bad("checksum mismatch"));
        }
        let mut c = Cursor { b: &b[..end], p: 8 };
        let k = c.byte().ok_or_else(|| bad("k"))? as usize;
        if !(4..=10).contains(&k) {
            return Err(bad("k out of range"));
        }
        let nr = c.varint().ok_or_else(|| bad("ranks"))? as usize;
        let mut ranks = Vec::with_capacity(nr);
        for _ in 0..nr {
            ranks.push(c.string().ok_or_else(|| bad("rank name"))?);
        }
        let n_sequences = c.varint().ok_or_else(|| bad("N"))? as u32;
        let skipped = c.varint().ok_or_else(|| bad("skipped"))?;
        let ncl = c.varint().ok_or_else(|| bad("classes"))? as usize;
        let mut lineages = Vec::with_capacity(ncl);
        let mut class_size = Vec::with_capacity(ncl);
        for _ in 0..ncl {
            class_size.push(c.varint().ok_or_else(|| bad("class size"))? as u32);
            let d = c.varint().ok_or_else(|| bad("lineage"))? as usize;
            let mut lin = Vec::with_capacity(d);
            for _ in 0..d {
                lin.push(c.string().ok_or_else(|| bad("lineage name"))?);
            }
            lineages.push(lin);
        }
        let source = c.string().ok_or_else(|| bad("source"))?;
        let nw = 1usize << (2 * k);
        let mut word_seqs = vec![0u32; nw];
        let mut post_start = vec![0u64; nw + 1];
        let mut post_class = Vec::new();
        let mut post_count = Vec::new();
        for w in 0..nw {
            word_seqs[w] = c.varint().ok_or_else(|| bad("word counts"))? as u32;
            let np = c.varint().ok_or_else(|| bad("postings"))? as usize;
            let mut prev = 0u64;
            for _ in 0..np {
                prev += c.varint().ok_or_else(|| bad("posting"))?;
                if prev as usize >= ncl {
                    return Err(bad("class index out of range"));
                }
                post_class.push(prev as u32);
                post_count.push(c.varint().ok_or_else(|| bad("posting count"))? as u32);
            }
            post_start[w + 1] = post_class.len() as u64;
        }
        if c.p != c.b.len() {
            return Err(bad("extra bytes"));
        }
        let mut m = TaxModel { k, ranks, lineages, class_size, n_sequences, word_seqs, skipped, source, post_start, post_class, post_count, ..Default::default() };
        m.finish();
        Ok(m)
    }
}

const MODEL_MAGIC: &[u8; 8] = b"MMGTAX01";
const MODEL_END: &[u8; 8] = b"MMGTAXEN";

fn put_varint(out: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        out.push((v as u8) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

fn put_str(out: &mut Vec<u8>, s: &str) {
    put_varint(out, s.len() as u64);
    out.extend_from_slice(s.as_bytes());
}

struct Cursor<'a> {
    b: &'a [u8],
    p: usize,
}

impl Cursor<'_> {
    fn byte(&mut self) -> Option<u8> {
        let x = *self.b.get(self.p)?;
        self.p += 1;
        Some(x)
    }
    fn varint(&mut self) -> Option<u64> {
        let mut v = 0u64;
        let mut shift = 0;
        loop {
            let x = self.byte()?;
            if shift > 63 {
                return None;
            }
            v |= u64::from(x & 0x7F) << shift;
            if x & 0x80 == 0 {
                return Some(v);
            }
            shift += 7;
        }
    }
    fn string(&mut self) -> Option<String> {
        let n = self.varint()? as usize;
        let s = self.b.get(self.p..self.p + n)?;
        self.p += n;
        String::from_utf8(s.to_vec()).ok()
    }
}

/// True when `path` is a model file saved by `TaxModel::save`.
pub fn is_model_file(path: &str) -> bool {
    let Ok(f) = File::open(path) else {
        return false;
    };
    let mut head = [0u8; 8];
    let mut d = flate2::read::GzDecoder::new(std::io::BufReader::new(f));
    d.read_exact(&mut head).is_ok() && &head == MODEL_MAGIC
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TrainSummary {
    pub reference: String,
    /// What the model says of its origin (file name, and the region when
    /// trained between primers).
    pub source: String,
    pub model: Option<String>,
    pub k: usize,
    pub sequences: u64,
    /// Records left out: no lineage, no valid word, or (with primers) no
    /// primer site; the last ones also counted in `without_region`.
    pub skipped: u64,
    pub without_region: u64,
    pub classes: u64,
    pub ranks: Vec<String>,
    /// Distinct taxa at each rank.
    pub taxa_per_rank: Vec<u64>,
    pub pairs: u64,
    pub model_bytes: Option<u64>,
    pub seconds: f64,
}

impl TrainSummary {
    pub fn json(&self) -> String {
        let t: Vec<String> = self.taxa_per_rank.iter().map(|x| x.to_string()).collect();
        format!(
            "{{\"reference\":{},\"source\":{},\"model\":{},\"k\":{},\"sequences\":{},\"skipped\":{},\"without_region\":{},\"classes\":{},\"ranks\":{},\"taxa_per_rank\":[{}],\"pairs\":{},\"model_bytes\":{},\"seconds\":{}}}",
            js(&self.reference),
            js(&self.source),
            self.model.as_deref().map_or("null".to_string(), js),
            self.k,
            self.sequences,
            self.skipped,
            self.without_region,
            self.classes,
            jlist(&self.ranks),
            t.join(","),
            self.pairs,
            self.model_bytes.map_or("null".to_string(), |x| x.to_string()),
            jf(self.seconds)
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!(
            "Classifier trained from {}: {} sequences, {} lineages ({}-mers), {} stored word-lineage pairs, {:.1} s.\n",
            self.reference,
            thousands(self.sequences),
            thousands(self.classes),
            self.k,
            thousands(self.pairs),
            self.seconds
        );
        let per: Vec<String> = self.ranks.iter().zip(&self.taxa_per_rank).map(|(r, n)| format!("{r} {n}")).collect();
        o.push_str(&format!("Taxa per rank: {}.\n", per.join(", ")));
        if self.without_region > 0 {
            o.push_str(&format!("Trained on the region between the primers: {} reference sequences without both primer sites were left out.\n", thousands(self.without_region)));
        }
        if self.skipped > self.without_region {
            o.push_str(&format!("warning: {} records were left out (no lineage in the header, or no valid word)\n", self.skipped - self.without_region));
        }
        if let (Some(m), Some(b)) = (&self.model, self.model_bytes) {
            o.push_str(&format!("Model written to {m} ({:.1} MB).\n", b as f64 / 1e6));
        }
        o
    }
}

fn summary_of(m: &TaxModel, reference: &str, model: Option<&str>, bytes: Option<u64>, seconds: f64) -> TrainSummary {
    let depth = m.ranks.len();
    let mut taxa_per_rank = vec![0u64; depth];
    for r in 0..depth {
        let set: HashSet<u32> = m.prefix_id.iter().filter_map(|ids| ids.get(r).copied()).collect();
        taxa_per_rank[r] = set.len() as u64;
    }
    TrainSummary {
        reference: reference.to_string(),
        source: m.source.clone(),
        model: model.map(str::to_string),
        k: m.k,
        sequences: u64::from(m.n_sequences),
        skipped: m.skipped,
        without_region: m.without_region,
        classes: m.classes() as u64,
        ranks: m.ranks.clone(),
        taxa_per_rank,
        pairs: m.pairs() as u64,
        model_bytes: bytes,
        seconds,
    }
}

/// Trains the classifier from a reference FASTA and, with `output`, saves
/// the model (train once, then classify fast).
pub fn train_classifier(reference: &str, output: Option<&str>, opt: &TrainOptions) -> Result<TrainSummary, String> {
    if let Some(o) = output {
        not_an_input(o, &[reference])?;
    }
    let t0 = std::time::Instant::now();
    let m = TaxModel::train(reference, opt)?;
    let bytes = match output {
        Some(o) => Some(m.save(o)?),
        None => None,
    };
    Ok(summary_of(&m, reference, output, bytes, t0.elapsed().as_secs_f64()))
}

/// A model file, or a reference FASTA to train from (with the default
/// training options).
pub fn load_classifier(path: &str) -> Result<TaxModel, String> {
    check_input(path)?;
    if is_model_file(path) {
        TaxModel::load(path)
    } else {
        TaxModel::train(path, &TrainOptions::default())
    }
}

#[derive(Clone, Debug)]
pub struct ClassifyOptions {
    /// Bootstrap trials (100).
    pub bootstraps: u32,
    /// A rank is kept when this fraction of the trials agree (0.8 as RDP;
    /// 0.5 as DADA2's minBoot = 50, for short reads).
    pub min_confidence: f64,
    pub seed: u64,
    /// Also try the reverse complement and keep the better strand.
    pub try_reverse_complement: bool,
    /// FASTA of reference sequences named by species, for exact species
    /// assignment (DADA2 addSpecies style).
    pub species_reference: Option<String>,
    /// Most species named for one ASV (1 = only unambiguous ones).
    pub max_species: usize,
}

impl Default for ClassifyOptions {
    fn default() -> Self {
        ClassifyOptions { bootstraps: 100, min_confidence: 0.8, seed: 42, try_reverse_complement: false, species_reference: None, max_species: 1 }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TaxonAssignment {
    pub id: String,
    /// The lineage kept: down to the deepest rank with enough confidence.
    pub lineage: Vec<String>,
    /// The best lineage with all the words, and the bootstrap confidence of
    /// each of its ranks.
    pub best: Vec<String>,
    pub confidence: Vec<f64>,
    /// The species by exact match (epithet; several joined by '/').
    pub species: Option<String>,
    pub reverse_complemented: bool,
    /// Distinct words of the sequence.
    pub words: usize,
}

/// Reference sequences named by species.
struct SpeciesRef {
    seqs: Vec<Vec<u8>>,
    genus: Vec<String>,
    epithet: Vec<String>,
}

fn load_species_reference(path: &str) -> Result<SpeciesRef, String> {
    let mut r = SeqIn::open(path)?;
    let mut rec = SeqRec::default();
    let mut s = SpeciesRef { seqs: Vec::new(), genus: Vec::new(), epithet: Vec::new() };
    let mut records = 0u64;
    while r.next(&mut rec)? {
        records += 1;
        let h = String::from_utf8_lossy(&rec.head).to_string();
        let first = h.split_whitespace().next().unwrap_or("");
        let pair: Option<(String, String)> = if h.contains(';') {
            parse_lineage(&rec.head).and_then(|lin| {
                if lin.len() < 7 {
                    return None;
                }
                let g = lin[lin.len() - 2].clone();
                let sp = lin[lin.len() - 1].clone();
                let ep = match sp.strip_prefix(&format!("{g} ")) {
                    Some(e) => e.to_string(),
                    None => sp.split_whitespace().last().unwrap_or(&sp).to_string(),
                };
                Some((g, ep))
            })
        } else {
            let w: Vec<&str> = h.split_whitespace().collect();
            if w.len() >= 3 && !first.is_empty() {
                Some((w[1].to_string(), w[2].to_string()))
            } else {
                None
            }
        };
        if let Some((g, e)) = pair {
            s.seqs.push(rec.seq.iter().map(|&b| UPPER[b as usize]).collect());
            s.genus.push(g);
            s.epithet.push(e);
        }
    }
    if s.seqs.is_empty() {
        return Err(format!(
            "{path}: none of the {records} records names a species. Give a species reference: headers \">ID Genus species\" (DADA2 assignSpecies files) or a lineage down to the species (\">Bacteria;...;Escherichia;coli;\")"
        ));
    }
    Ok(s)
}

fn genus_matches(assigned: &str, reference: &str) -> bool {
    let a = assigned.trim();
    let r = reference.trim();
    if a.eq_ignore_ascii_case(r) {
        return true;
    }
    a.split('/').any(|x| x.trim().eq_ignore_ascii_case(r)) || r.split('/').any(|x| x.trim().eq_ignore_ascii_case(a))
}

/// Exact species assignment. `genera[i]` is the genus the classifier gave
/// sequence i (None = no genus: no species).
fn assign_species(seqs: &[Vec<u8>], genera: &[Option<String>], sref: &SpeciesRef, try_rc: bool, max_species: usize) -> Vec<Option<String>> {
    const K: usize = 32;
    // queries indexed by their first 32 bases
    let mut index: HashMap<u64, Vec<(usize, bool)>> = HashMap::new();
    let mut short: Vec<(usize, bool)> = Vec::new();
    let mut queries: Vec<Vec<Vec<u8>>> = Vec::with_capacity(seqs.len());
    let code_of = |s: &[u8]| -> Option<u64> {
        let mut c = 0u64;
        for &b in &s[..K] {
            let v = CODE[b as usize];
            if v > 3 {
                return None;
            }
            c = (c << 2) | u64::from(v);
        }
        Some(c)
    };
    for (i, s) in seqs.iter().enumerate() {
        let up: Vec<u8> = s.iter().map(|&b| UPPER[b as usize]).collect();
        let mut forms = vec![up.clone()];
        if try_rc {
            forms.push(revcomp(&up));
        }
        for (fi, f) in forms.iter().enumerate() {
            match if f.len() >= K { code_of(f) } else { None } {
                Some(c) => index.entry(c).or_default().push((i, fi == 1)),
                None => short.push((i, fi == 1)),
            }
        }
        queries.push(forms);
    }
    let hits: Vec<(usize, usize)> = sref
        .seqs
        .par_iter()
        .enumerate()
        .flat_map_iter(|(ri, r)| {
            let mut out = Vec::new();
            if r.len() >= K {
                let mask = u64::MAX; // K = 32 bases fill the 64 bits
                let (mut c, mut valid) = (0u64, 0usize);
                for (p, &b) in r.iter().enumerate() {
                    let v = CODE[b as usize];
                    if v > 3 {
                        valid = 0;
                        continue;
                    }
                    c = ((c << 2) | u64::from(v)) & mask;
                    valid += 1;
                    if valid >= K {
                        let start = p + 1 - K;
                        if let Some(list) = index.get(&c) {
                            for &(qi, rc) in list {
                                let q = &queries[qi][usize::from(rc)];
                                if start + q.len() <= r.len() && &r[start..start + q.len()] == q.as_slice() {
                                    out.push((qi, ri));
                                }
                            }
                        }
                    }
                }
            }
            for &(qi, rc) in &short {
                let q = &queries[qi][usize::from(rc)];
                if !q.is_empty() && r.windows(q.len()).any(|w| w == q.as_slice()) {
                    out.push((qi, ri));
                }
            }
            out
        })
        .collect();
    let mut found: Vec<Vec<usize>> = vec![Vec::new(); seqs.len()];
    for (qi, ri) in hits {
        found[qi].push(ri);
    }
    found
        .iter()
        .enumerate()
        .map(|(qi, refs)| {
            let g = genera[qi].as_ref()?;
            let mut names: Vec<String> = refs.iter().filter(|&&ri| genus_matches(g, &sref.genus[ri])).map(|&ri| sref.epithet[ri].clone()).collect();
            names.sort();
            names.dedup();
            if names.is_empty() || names.len() > max_species.max(1) {
                None
            } else {
                Some(names.join("/"))
            }
        })
        .collect()
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TaxonomyResult {
    pub model: String,
    pub output: String,
    pub ranks: Vec<String>,
    pub assignments: Vec<TaxonAssignment>,
    /// Sequences assigned down to each rank.
    pub assigned_per_rank: Vec<u64>,
    /// The name of the exact-species column ("Species", or "Species_exact"
    /// when the model has a Species rank); None without a species reference.
    pub species_column: Option<String>,
    pub species_assigned: u64,
    pub min_confidence: f64,
    pub seconds: f64,
}

/// Classifies every sequence of a FASTA (the ASVs) with a model file (or a
/// reference FASTA, trained on the fly) and writes the taxonomy table: ASV,
/// one column per rank (empty below the confidence threshold), the exact
/// species if asked, then `<rank>_confidence` for each rank of the best
/// lineage.
pub fn classify_asvs(sequences: &str, model: &str, output: &str, opt: &ClassifyOptions) -> Result<TaxonomyResult, String> {
    let mut ins = vec![sequences, model];
    if let Some(s) = &opt.species_reference {
        ins.push(s);
    }
    not_an_input(output, &ins)?;
    if !(opt.min_confidence >= 0.0 && opt.min_confidence <= 1.0) {
        return Err(format!("min_confidence = {}: give a fraction from 0 to 1 (0.8 as RDP, 0.5 for short reads)", opt.min_confidence));
    }
    if opt.bootstraps == 0 {
        return Err("bootstraps = 0: give the number of bootstrap trials (100 is usual)".to_string());
    }
    let t0 = std::time::Instant::now();
    let queries = read_fasta(sequences)?;
    let m = load_classifier(model)?;
    let mut assignments: Vec<TaxonAssignment> = queries
        .par_iter()
        .map(|(id, s)| {
            let mut a = m.classify(s, opt);
            a.id = id.clone();
            a
        })
        .collect();
    let genus_rank = m.ranks.iter().position(|r| r.eq_ignore_ascii_case("genus"));
    let mut species_column = None;
    if let Some(sp) = &opt.species_reference {
        let sref = load_species_reference(sp)?;
        let genera: Vec<Option<String>> = assignments
            .iter()
            .map(|a| match genus_rank {
                Some(g) => a.lineage.get(g).cloned(),
                None => a.lineage.last().cloned(),
            })
            .collect();
        let seqs: Vec<Vec<u8>> = queries.iter().map(|(_, s)| s.clone()).collect();
        let sp = assign_species(&seqs, &genera, &sref, opt.try_reverse_complement, opt.max_species);
        for (a, s) in assignments.iter_mut().zip(sp) {
            a.species = s;
        }
        species_column = Some(if m.ranks.iter().any(|r| r.eq_ignore_ascii_case("species")) { "Species_exact".to_string() } else { "Species".to_string() });
    }
    let depth = m.ranks.len();
    let mut text = String::from("ASV");
    for r in &m.ranks {
        text.push('\t');
        text.push_str(r);
    }
    if let Some(c) = &species_column {
        text.push('\t');
        text.push_str(c);
    }
    for r in &m.ranks {
        text.push_str(&format!("\t{r}_confidence"));
    }
    text.push('\n');
    let mut assigned_per_rank = vec![0u64; depth];
    for a in &assignments {
        text.push_str(&a.id);
        for r in 0..depth {
            text.push('\t');
            if let Some(n) = a.lineage.get(r) {
                text.push_str(n);
                assigned_per_rank[r] += 1;
            }
        }
        if species_column.is_some() {
            text.push('\t');
            if let Some(s) = &a.species {
                text.push_str(s);
            }
        }
        for r in 0..depth {
            text.push('\t');
            if let Some(c) = a.confidence.get(r) {
                text.push_str(&format!("{c:.2}"));
            }
        }
        text.push('\n');
    }
    write_text(output, &text)?;
    let species_assigned = assignments.iter().filter(|a| a.species.is_some()).count() as u64;
    Ok(TaxonomyResult {
        model: model.to_string(),
        output: output.to_string(),
        ranks: m.ranks.clone(),
        assignments,
        assigned_per_rank,
        species_column,
        species_assigned,
        min_confidence: opt.min_confidence,
        seconds: t0.elapsed().as_secs_f64(),
    })
}

impl TaxonomyResult {
    pub fn json(&self) -> String {
        let a: Vec<String> = self
            .assignments
            .iter()
            .map(|a| {
                format!(
                    "{{\"id\":{},\"lineage\":{},\"best\":{},\"confidence\":{},\"species\":{},\"reverse_complemented\":{},\"words\":{}}}",
                    js(&a.id),
                    jlist(&a.lineage),
                    jlist(&a.best),
                    jnums(&a.confidence),
                    a.species.as_deref().map_or("null".to_string(), js),
                    a.reverse_complemented,
                    a.words
                )
            })
            .collect();
        let per: Vec<String> = self.assigned_per_rank.iter().map(|x| x.to_string()).collect();
        format!(
            "{{\"model\":{},\"output\":{},\"ranks\":{},\"assignments\":[{}],\"assigned_per_rank\":[{}],\"species_column\":{},\"species_assigned\":{},\"min_confidence\":{},\"seconds\":{}}}",
            js(&self.model),
            js(&self.output),
            jlist(&self.ranks),
            a.join(","),
            per.join(","),
            self.species_column.as_deref().map_or("null".to_string(), js),
            self.species_assigned,
            jf(self.min_confidence),
            jf(self.seconds)
        )
    }

    pub fn resumo(&self) -> String {
        let n = self.assignments.len() as u64;
        let mut o = format!("{} sequences classified (confidence >= {}), {:.1} s. Assigned down to each rank:\n", n, self.min_confidence, self.seconds);
        for (r, c) in self.ranks.iter().zip(&self.assigned_per_rank) {
            o.push_str(&format!("  {:<10} {:>6} ({})\n", r, c, pct(*c, n)));
        }
        if let Some(c) = &self.species_column {
            o.push_str(&format!("  {:<10} {:>6} ({}) by exact match\n", c, self.species_assigned, pct(self.species_assigned, n)));
        }
        o.push_str(&format!("Written: {}\n", self.output));
        o
    }
}

// ---------------------------------------------------------------------
//  FEATURE TABLES AND METADATA
// ---------------------------------------------------------------------

/// A feature table: features (ASVs) in rows, samples in columns.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FeatureTable {
    pub features: Vec<String>,
    pub samples: Vec<String>,
    /// counts[feature][sample]
    pub counts: Vec<Vec<f64>>,
    /// Every value is a whole number.
    pub integer: bool,
}

impl FeatureTable {
    pub fn sample_column(&self, s: usize) -> Vec<f64> {
        self.counts.iter().map(|r| r[s]).collect()
    }
    pub fn sample_counts(&self, s: usize) -> Vec<u64> {
        self.counts.iter().map(|r| r[s].round() as u64).collect()
    }
}

/// Reads a feature table (TSV; CSV if the header has commas and no tab):
/// a header line (feature id, then the sample names; a header starting with
/// '#', as BIOM's "#OTU ID", is fine), then one line per feature. Comment
/// lines start with '#' and have no tab. A last column named "taxonomy" is
/// ignored.
pub fn read_feature_table(path: &str) -> Result<FeatureTable, String> {
    check_input(path)?;
    let mut text = String::new();
    abre_leitura(path)?.read_to_string(&mut text).map_err(|e| format!("Error reading the feature table {path}: {e}. Is it a text table?"))?;
    let mut t = FeatureTable { integer: true, ..Default::default() };
    let mut delim = '\t';
    let mut header_done = false;
    let mut ncol = 0usize;
    for (i, raw) in text.lines().enumerate() {
        let line = raw.trim_end_matches('\r');
        if line.trim().is_empty() {
            continue;
        }
        if !header_done {
            if line.starts_with('#') && !line.contains('\t') {
                continue;
            }
            if !line.contains('\t') && line.contains(',') {
                delim = ',';
            }
            let mut f: Vec<String> = line.split(delim).map(|s| s.trim().to_string()).collect();
            if f.last().map_or(false, |x| x.eq_ignore_ascii_case("taxonomy")) {
                f.pop();
            }
            if f.len() < 2 {
                return Err(format!(
                    "{path}, line {}: the header has no sample column. A feature table has a header (feature id, then one column per sample, tab-separated) and one line per feature",
                    i + 1
                ));
            }
            ncol = f.len();
            t.samples = f[1..].to_vec();
            let mut seen = HashSet::new();
            for s in &t.samples {
                if !seen.insert(s.as_str()) {
                    return Err(format!("{path}: the sample {s} appears twice in the header. Each column needs its own sample name"));
                }
            }
            header_done = true;
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split(delim).map(str::trim).collect();
        if f.len() < ncol {
            return Err(format!("{path}, line {}: {} columns, but the header has {ncol}. Every line needs one value per sample", i + 1, f.len()));
        }
        let mut row = Vec::with_capacity(ncol - 1);
        for (j, v) in f[1..ncol].iter().enumerate() {
            let x: f64 = v.parse().map_err(|_| {
                format!("{path}, line {}, sample {}: \"{v}\" is not a number. A feature table holds counts (or proportions)", i + 1, t.samples[j])
            })?;
            if !(x >= 0.0) || !x.is_finite() {
                return Err(format!("{path}, line {}, sample {}: {v} is negative or not finite; counts cannot be", i + 1, t.samples[j]));
            }
            if x != x.trunc() {
                t.integer = false;
            }
            row.push(x);
        }
        t.features.push(f[0].to_string());
        t.counts.push(row);
    }
    if !header_done {
        return Err(format!("{path} is empty: a feature table needs a header and one line per feature"));
    }
    Ok(t)
}

/// A metadata table: the first column holds the sample ids; returns the
/// values of the chosen column (default: the second), by sample.
pub fn read_metadata_column(path: &str, column: Option<&str>) -> Result<(String, HashMap<String, String>), String> {
    check_input(path)?;
    let mut text = String::new();
    abre_leitura(path)?.read_to_string(&mut text).map_err(|e| format!("Error reading the metadata {path}: {e}"))?;
    let mut header: Option<Vec<String>> = None;
    let mut col = 0usize;
    let mut out = HashMap::new();
    for (i, raw) in text.lines().enumerate() {
        let line = raw.trim_end_matches('\r');
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<String> = line.split('\t').map(|s| s.trim().to_string()).collect();
        match &header {
            None => {
                if line.starts_with('#') && !line.contains('\t') {
                    continue;
                }
                if f.len() < 2 {
                    return Err(format!("{path}, line {}: the metadata needs at least two tab-separated columns: the sample id, then the variables", i + 1));
                }
                col = match column {
                    Some(c) => f.iter().position(|h| h.eq_ignore_ascii_case(c)).filter(|&p| p > 0).ok_or_else(|| {
                        format!("{path}: no column named \"{c}\". The columns are: {}", f[1..].join(", "))
                    })?,
                    None => 1,
                };
                header = Some(f);
            }
            Some(_) => {
                if line.starts_with('#') {
                    continue;
                }
                if let Some(v) = f.get(col) {
                    if !v.is_empty() {
                        out.insert(f[0].clone(), v.clone());
                    }
                }
            }
        }
    }
    let h = header.ok_or_else(|| format!("{path} is empty: the metadata needs a header and one line per sample"))?;
    Ok((h[col].clone(), out))
}

// ---------------------------------------------------------------------
//  ALPHA DIVERSITY AND RAREFACTION
// ---------------------------------------------------------------------

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AlphaRow {
    pub sample: String,
    pub depth: u64,
    pub observed: u64,
    pub chao1: f64,
    pub shannon: f64,
    pub simpson: f64,
    /// Shannon / ln(observed); None with fewer than 2 features.
    pub pielou: Option<f64>,
    /// 1 - singletons / depth; None for an empty sample.
    pub goods_coverage: Option<f64>,
    pub singletons: u64,
    pub doubletons: u64,
}

/// The alpha indices of one sample's counts.
pub fn alpha_of(counts: &[u64]) -> AlphaRow {
    let n: u64 = counts.iter().sum();
    let s = counts.iter().filter(|&&c| c > 0).count() as u64;
    let f1 = counts.iter().filter(|&&c| c == 1).count() as u64;
    let f2 = counts.iter().filter(|&&c| c == 2).count() as u64;
    let (mut h, mut d) = (0.0, 0.0);
    if n > 0 {
        for &c in counts.iter().filter(|&&c| c > 0) {
            let p = c as f64 / n as f64;
            h -= p * p.ln();
            d += p * p;
        }
    }
    AlphaRow {
        sample: String::new(),
        depth: n,
        observed: s,
        chao1: s as f64 + (f1 as f64 * (f1 as f64 - 1.0)) / (2.0 * (f2 as f64 + 1.0)),
        shannon: h,
        simpson: if n > 0 { 1.0 - d } else { 0.0 },
        pielou: if s > 1 { Some(h / (s as f64).ln()) } else { None },
        goods_coverage: if n > 0 { Some(1.0 - f1 as f64 / n as f64) } else { None },
        singletons: f1,
        doubletons: f2,
    }
}

/// Subsamples `counts` to `depth` reads without replacement (a partial
/// Fisher-Yates shuffle of the reads). `depth` >= the total returns the
/// counts unchanged.
pub fn rarefy(counts: &[u64], depth: u64, rng: &mut Rng) -> Vec<u64> {
    let total: u64 = counts.iter().sum();
    if depth >= total {
        return counts.to_vec();
    }
    let mut reads: Vec<u32> = Vec::with_capacity(total as usize);
    for (i, &c) in counts.iter().enumerate() {
        reads.extend(std::iter::repeat(i as u32).take(c as usize));
    }
    let mut out = vec![0u64; counts.len()];
    for i in 0..depth as usize {
        let j = i + rng.below((reads.len() - i) as u64) as usize;
        reads.swap(i, j);
        out[reads[i] as usize] += 1;
    }
    out
}

/// ln Gamma(x) for x > 0 (Lanczos, g = 7).
pub fn ln_gamma(x: f64) -> f64 {
    const G: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        let pi = std::f64::consts::PI;
        return (pi / (pi * x).sin().abs()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = G[0];
    let t = x + 7.5;
    for (i, g) in G.iter().enumerate().skip(1) {
        a += g / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// The expected number of features seen when `n` reads are drawn without
/// replacement (Hurlbert 1971): sum over features of 1 - C(N - c, n) / C(N, n).
pub fn expected_observed(counts: &[u64], n: u64) -> f64 {
    let total: u64 = counts.iter().sum();
    if n >= total {
        return counts.iter().filter(|&&c| c > 0).count() as f64;
    }
    let (nt, nn) = (total as f64, n as f64);
    let mut s = 0.0;
    for &c in counts.iter().filter(|&&c| c > 0) {
        if total - c < n {
            s += 1.0;
            continue;
        }
        // C(N-c, n) / C(N, n) = prod_{j<c} (N - n - j) / (N - j)
        //                      = prod_{j<n} (N - c - j) / (N - j)
        let cc = c as f64;
        let ln_ratio = if c <= n && c <= 5000 {
            (0..c).map(|j| ((nt - nn - j as f64) / (nt - j as f64)).ln()).sum::<f64>()
        } else if n <= 5000 {
            (0..n).map(|j| ((nt - cc - j as f64) / (nt - j as f64)).ln()).sum::<f64>()
        } else {
            let cc = c as f64;
            ln_gamma(nt - cc + 1.0) - ln_gamma(nt - cc - nn + 1.0) - ln_gamma(nt + 1.0) + ln_gamma(nt - nn + 1.0)
        };
        s += 1.0 - ln_ratio.exp();
    }
    s
}

#[derive(Clone, Debug)]
pub struct AlphaOptions {
    /// Rarefy every sample to this depth first (samples with fewer reads
    /// are left out); None = use the counts as they are.
    pub rarefy_depth: Option<u64>,
    pub seed: u64,
    /// Points of each rarefaction curve (20); 0 = no curves.
    pub curve_steps: usize,
    /// Deepest point of the curves; None = the deepest sample.
    pub curve_max_depth: Option<u64>,
    /// TSV of the curves (sample, depth, expected observed).
    pub curves_output: Option<String>,
}

impl Default for AlphaOptions {
    fn default() -> Self {
        AlphaOptions { rarefy_depth: None, seed: 42, curve_steps: 20, curve_max_depth: None, curves_output: None }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct RarefactionCurve {
    pub sample: String,
    pub depths: Vec<u64>,
    pub observed: Vec<f64>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AlphaResult {
    pub table: String,
    pub output: Option<String>,
    pub rows: Vec<AlphaRow>,
    pub rarefied_to: Option<u64>,
    /// Samples left out for having fewer reads than the rarefaction depth.
    pub dropped: Vec<String>,
    pub curves: Vec<RarefactionCurve>,
}

fn check_counts(t: &FeatureTable, path: &str, what: &str) -> Result<(), String> {
    if !t.integer {
        return Err(format!(
            "{path} holds fractional values, but {what} needs counts of reads (whole numbers). Give the feature table of counts, not of proportions"
        ));
    }
    Ok(())
}

/// Alpha diversity of every sample of a feature table, with rarefaction if
/// asked, and the rarefaction curves.
pub fn alpha_diversity(table: &str, output: Option<&str>, opt: &AlphaOptions) -> Result<AlphaResult, String> {
    let mut outs: Vec<&str> = output.into_iter().collect();
    if let Some(c) = &opt.curves_output {
        outs.push(c);
    }
    for o in &outs {
        not_an_input(o, &[table])?;
    }
    let t = read_feature_table(table)?;
    check_counts(&t, table, "alpha diversity")?;
    if opt.rarefy_depth == Some(0) {
        return Err("rarefy_depth = 0: give the number of reads to keep per sample (for instance the smallest sample's depth)".to_string());
    }
    let mut rows = Vec::new();
    let mut dropped = Vec::new();
    for (s, name) in t.samples.iter().enumerate() {
        let mut c = t.sample_counts(s);
        if let Some(d) = opt.rarefy_depth {
            let total: u64 = c.iter().sum();
            if total < d {
                dropped.push(name.clone());
                continue;
            }
            c = rarefy(&c, d, &mut Rng::new(seed_for(opt.seed, name.as_bytes())));
        }
        let mut r = alpha_of(&c);
        r.sample = name.clone();
        rows.push(r);
    }
    let mut curves = Vec::new();
    if opt.curve_steps > 0 {
        let deepest = (0..t.samples.len()).map(|s| t.sample_counts(s).iter().sum::<u64>()).max().unwrap_or(0);
        let top = opt.curve_max_depth.unwrap_or(deepest).max(1);
        for (s, name) in t.samples.iter().enumerate() {
            let c = t.sample_counts(s);
            let total: u64 = c.iter().sum();
            let mut depths: Vec<u64> = (1..=opt.curve_steps as u64).map(|i| (top as f64 * i as f64 / opt.curve_steps as f64).round() as u64).filter(|&d| d > 0 && d < total).collect();
            if total > 0 && total <= top {
                depths.push(total);
            }
            depths.dedup();
            let observed = depths.iter().map(|&d| expected_observed(&c, d)).collect();
            curves.push(RarefactionCurve { sample: name.clone(), depths, observed });
        }
    }
    if rows.is_empty() {
        return Err(format!(
            "No sample of {table} has {} reads: every sample was left out. Choose a smaller rarefaction depth",
            opt.rarefy_depth.unwrap_or(0)
        ));
    }
    let mut cleanup = Cleanup(Vec::new());
    if let Some(o) = output {
        let mut s = String::from("sample\tdepth\tobserved\tchao1\tshannon\tsimpson\tpielou\tgoods_coverage\tsingletons\tdoubletons\n");
        for r in &rows {
            s.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                r.sample,
                r.depth,
                r.observed,
                num(r.chao1),
                num(r.shannon),
                num(r.simpson),
                r.pielou.map_or("NA".to_string(), num),
                r.goods_coverage.map_or("NA".to_string(), num),
                r.singletons,
                r.doubletons
            ));
        }
        cleanup.0.push(o.to_string());
        write_text(o, &s)?;
    }
    if let Some(o) = &opt.curves_output {
        let mut s = String::from("sample\tdepth\texpected_observed\n");
        for c in &curves {
            for (d, x) in c.depths.iter().zip(&c.observed) {
                s.push_str(&format!("{}\t{}\t{}\n", c.sample, d, num(*x)));
            }
        }
        cleanup.0.push(o.clone());
        write_text(o, &s)?;
    }
    cleanup.disarm();
    Ok(AlphaResult { table: table.to_string(), output: output.map(str::to_string), rows, rarefied_to: opt.rarefy_depth, dropped, curves })
}

impl AlphaRow {
    fn json(&self) -> String {
        format!(
            "{{\"sample\":{},\"depth\":{},\"observed\":{},\"chao1\":{},\"shannon\":{},\"simpson\":{},\"pielou\":{},\"goods_coverage\":{},\"singletons\":{},\"doubletons\":{}}}",
            js(&self.sample),
            self.depth,
            self.observed,
            jf(self.chao1),
            jf(self.shannon),
            jf(self.simpson),
            jo(self.pielou),
            jo(self.goods_coverage),
            self.singletons,
            self.doubletons
        )
    }
}

impl AlphaResult {
    pub fn json(&self) -> String {
        let r: Vec<String> = self.rows.iter().map(AlphaRow::json).collect();
        let c: Vec<String> = self
            .curves
            .iter()
            .map(|c| {
                let d: Vec<String> = c.depths.iter().map(|x| x.to_string()).collect();
                format!("{{\"sample\":{},\"depths\":[{}],\"observed\":{}}}", js(&c.sample), d.join(","), jnums(&c.observed))
            })
            .collect();
        format!(
            "{{\"table\":{},\"output\":{},\"rows\":[{}],\"rarefied_to\":{},\"dropped\":{},\"curves\":[{}]}}",
            js(&self.table),
            self.output.as_deref().map_or("null".to_string(), js),
            r.join(","),
            self.rarefied_to.map_or("null".to_string(), |x| x.to_string()),
            jlist(&self.dropped),
            c.join(",")
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = match self.rarefied_to {
            Some(d) => format!("Alpha diversity of {} samples, rarefied to {} reads:\n", self.rows.len(), thousands(d)),
            None => format!("Alpha diversity of {} samples:\n", self.rows.len()),
        };
        o.push_str(&format!("  {:<16} {:>8} {:>8} {:>8} {:>8} {:>8} {:>7} {:>7}\n", "sample", "depth", "observed", "chao1", "shannon", "simpson", "pielou", "goods"));
        for r in &self.rows {
            o.push_str(&format!(
                "  {:<16} {:>8} {:>8} {:>8.1} {:>8.3} {:>8.3} {:>7} {:>7}\n",
                r.sample,
                r.depth,
                r.observed,
                r.chao1,
                r.shannon,
                r.simpson,
                r.pielou.map_or("NA".to_string(), |x| format!("{x:.3}")),
                r.goods_coverage.map_or("NA".to_string(), |x| format!("{x:.4}"))
            ));
        }
        if !self.dropped.is_empty() {
            o.push_str(&format!("Left out (fewer reads than the rarefaction depth): {}\n", self.dropped.join(", ")));
        }
        o
    }
}

// ---------------------------------------------------------------------
//  BETA DIVERSITY, PCoA, PERMANOVA
// ---------------------------------------------------------------------

/// Bray-Curtis dissimilarity: sum |a - b| / sum (a + b); 0 for two empty
/// samples.
pub fn bray_curtis(a: &[f64], b: &[f64]) -> f64 {
    let (mut num, mut den) = (0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        num += (x - y).abs();
        den += x + y;
    }
    if den > 0.0 {
        num / den
    } else {
        0.0
    }
}

/// Jaccard distance on presence: 1 - |A and B| / |A or B|; 0 for two empty
/// samples.
pub fn jaccard(a: &[f64], b: &[f64]) -> f64 {
    let (mut both, mut any) = (0u64, 0u64);
    for (x, y) in a.iter().zip(b) {
        let (p, q) = (*x > 0.0, *y > 0.0);
        if p && q {
            both += 1;
        }
        if p || q {
            any += 1;
        }
    }
    if any > 0 {
        1.0 - both as f64 / any as f64
    } else {
        0.0
    }
}

/// Eigenvalues and eigenvectors (columns of the second matrix) of a
/// symmetric matrix, by cyclic Jacobi rotations.
pub fn symmetric_eigen(m: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = m.len();
    let mut a: Vec<Vec<f64>> = m.to_vec();
    let mut v = vec![vec![0.0; n]; n];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    let norm: f64 = a.iter().flatten().map(|x| x * x).sum::<f64>().sqrt();
    for _sweep in 0..100 {
        let off: f64 = (0..n).flat_map(|i| (i + 1..n).map(move |j| (i, j))).map(|(i, j)| a[i][j] * a[i][j]).sum::<f64>().sqrt();
        if off <= 1e-15 * norm.max(1e-300) {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                let apq = a[p][q];
                if apq.abs() <= 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for row in a.iter_mut() {
                    let (akp, akq) = (row[p], row[q]);
                    row[p] = c * akp - s * akq;
                    row[q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for row in v.iter_mut() {
                    let (vkp, vkq) = (row[p], row[q]);
                    row[p] = c * vkp - s * vkq;
                    row[q] = s * vkp + c * vkq;
                }
            }
        }
    }
    ((0..n).map(|i| a[i][i]).collect(), v)
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Pcoa {
    /// Every eigenvalue, from the largest (negative ones included: they add
    /// up to the total variance).
    pub eigenvalues: Vec<f64>,
    /// Each positive eigenvalue over the sum of the positive ones (0 for the
    /// others), as scikit-bio.
    pub proportion_explained: Vec<f64>,
    /// coordinates[sample][axis], one axis per positive eigenvalue.
    pub coordinates: Vec<Vec<f64>>,
    pub negative_eigenvalues: usize,
}

/// Principal coordinates of a distance matrix.
pub fn pcoa(d: &[Vec<f64>]) -> Pcoa {
    let n = d.len();
    if n == 0 {
        return Pcoa::default();
    }
    let a: Vec<Vec<f64>> = d.iter().map(|r| r.iter().map(|x| -0.5 * x * x).collect()).collect();
    let rm: Vec<f64> = a.iter().map(|r| r.iter().sum::<f64>() / n as f64).collect();
    let cm: Vec<f64> = (0..n).map(|j| a.iter().map(|r| r[j]).sum::<f64>() / n as f64).collect();
    let gm: f64 = rm.iter().sum::<f64>() / n as f64;
    let b: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| a[i][j] - rm[i] - cm[j] + gm).collect()).collect();
    let (vals, vecs) = symmetric_eigen(&b);
    let order = sorted_indices(n, |&x, &y| vals[y as usize].partial_cmp(&vals[x as usize]).unwrap_or(Ordering::Equal));
    let eig: Vec<f64> = order.iter().map(|&i| vals[i as usize]).collect();
    let scale: f64 = eig.iter().map(|x| x.abs()).fold(0.0, f64::max);
    let tol = 1e-10 * scale.max(1e-300);
    // proportions over the positive eigenvalues (scikit-bio's convention:
    // the negative ones explain nothing)
    let sum: f64 = eig.iter().filter(|&&x| x > tol).sum();
    let positive: Vec<usize> = (0..n).filter(|&k| eig[k] > tol).collect();
    let mut coords = vec![vec![0.0; positive.len()]; n];
    for (ax, &k) in positive.iter().enumerate() {
        let col = order[k] as usize;
        // sign: the largest component positive, for a reproducible picture
        let big = (0..n).max_by(|&x, &y| vecs[x][col].abs().partial_cmp(&vecs[y][col].abs()).unwrap_or(Ordering::Equal)).unwrap_or(0);
        let sign = if vecs[big][col] < 0.0 { -1.0 } else { 1.0 };
        for (i, row) in coords.iter_mut().enumerate() {
            row[ax] = sign * vecs[i][col] * eig[k].sqrt();
        }
    }
    Pcoa {
        proportion_explained: eig.iter().map(|&x| if sum > 0.0 && x > tol { x / sum } else { 0.0 }).collect(),
        negative_eigenvalues: eig.iter().filter(|&&x| x < -tol).count(),
        eigenvalues: eig,
        coordinates: coords,
    }
}

/// PERMANOVA (Anderson 2001): pseudo-F, R^2 and the permutation p-value.
/// `groups[i]` is the group number of sample i.
pub fn permanova(d: &[Vec<f64>], groups: &[usize], permutations: u32, seed: u64) -> (f64, f64, f64) {
    let n = d.len();
    let a = groups.iter().copied().max().map_or(0, |m| m + 1);
    let mut sizes = vec![0usize; a];
    for &g in groups {
        sizes[g] += 1;
    }
    let ss_t: f64 = (0..n).flat_map(|i| (i + 1..n).map(move |j| (i, j))).map(|(i, j)| d[i][j] * d[i][j]).sum::<f64>() / n as f64;
    let ss_w = |g: &[usize]| -> f64 {
        let mut s = vec![0.0; a];
        for i in 0..n {
            for j in i + 1..n {
                if g[i] == g[j] {
                    s[g[i]] += d[i][j] * d[i][j];
                }
            }
        }
        s.iter().zip(&sizes).map(|(x, &m)| if m > 0 { x / m as f64 } else { 0.0 }).sum()
    };
    let f_of = |w: f64| ((ss_t - w) / (a as f64 - 1.0)) / (w / (n as f64 - a as f64));
    let w0 = ss_w(groups);
    let f = f_of(w0);
    let r2 = if ss_t > 0.0 { (ss_t - w0) / ss_t } else { 0.0 };
    if permutations == 0 {
        return (f, r2, f64::NAN);
    }
    let mut rng = Rng::new(seed);
    let mut g = groups.to_vec();
    let mut hits = 0u32;
    for _ in 0..permutations {
        for i in (1..n).rev() {
            let j = rng.below(i as u64 + 1) as usize;
            g.swap(i, j);
        }
        let fp = f_of(ss_w(&g));
        if fp >= f - 1e-12 * f.abs().max(1.0) {
            hits += 1;
        }
    }
    (f, r2, (f64::from(hits) + 1.0) / (f64::from(permutations) + 1.0))
}

#[derive(Clone, Debug)]
pub struct BetaOptions {
    /// "braycurtis", "jaccard" or "kmer" (Bray-Curtis between the
    /// abundance-weighted 6-mer profiles; needs `asv_fasta`).
    pub metric: String,
    /// Bray-Curtis on proportions instead of counts.
    pub relative: bool,
    pub rarefy_depth: Option<u64>,
    pub seed: u64,
    pub asv_fasta: Option<String>,
    /// Metadata TSV (first column = sample id) and the column of groups,
    /// for PERMANOVA (default column: the second).
    pub metadata: Option<String>,
    pub group_column: Option<String>,
    pub permutations: u32,
    /// TSV of the principal coordinates.
    pub pcoa_output: Option<String>,
}

impl Default for BetaOptions {
    fn default() -> Self {
        BetaOptions {
            metric: "braycurtis".to_string(),
            relative: false,
            rarefy_depth: None,
            seed: 42,
            asv_fasta: None,
            metadata: None,
            group_column: None,
            permutations: 999,
            pcoa_output: None,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Permanova {
    pub column: String,
    pub groups: Vec<String>,
    pub sizes: Vec<usize>,
    pub n: usize,
    pub f: f64,
    pub r2: f64,
    pub p_value: f64,
    pub permutations: u32,
    /// Samples without a value in the metadata column.
    pub excluded: Vec<String>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct BetaResult {
    pub metric: String,
    pub samples: Vec<String>,
    pub distances: Vec<Vec<f64>>,
    pub pcoa: Pcoa,
    pub permanova: Option<Permanova>,
    /// Samples left out: empty, or below the rarefaction depth.
    pub dropped: Vec<String>,
    pub rarefied_to: Option<u64>,
    pub output: Option<String>,
    pub pcoa_output: Option<String>,
}

/// 6-mer profile of a sequence (fractions summing to 1).
fn kmer_profile(s: &[u8]) -> Vec<f64> {
    let mut v = vec![0.0; 4096];
    let (mut code, mut valid, mut n) = (0usize, 0usize, 0.0);
    for &b in s {
        let c = CODE[b as usize];
        if c > 3 {
            valid = 0;
            continue;
        }
        code = ((code << 2) | c as usize) & 4095;
        valid += 1;
        if valid >= 6 {
            v[code] += 1.0;
            n += 1.0;
        }
    }
    if n > 0.0 {
        for x in v.iter_mut() {
            *x /= n;
        }
    }
    v
}

/// Beta diversity between the samples of a feature table, the principal
/// coordinates, and PERMANOVA when a metadata column is given.
pub fn beta_diversity(table: &str, output: Option<&str>, opt: &BetaOptions) -> Result<BetaResult, String> {
    let metric = opt.metric.trim().to_ascii_lowercase().replace(['-', '_', ' '], "");
    let metric = match metric.as_str() {
        "braycurtis" | "bc" | "bray" => "braycurtis",
        "jaccard" => "jaccard",
        "kmer" | "kmerbraycurtis" => "kmer",
        _ => return Err(format!("metric = \"{}\": choose braycurtis, jaccard or kmer", opt.metric)),
    };
    let mut ins = vec![table];
    if let Some(m) = &opt.metadata {
        ins.push(m);
    }
    if let Some(a) = &opt.asv_fasta {
        ins.push(a);
    }
    for o in output.into_iter().chain(opt.pcoa_output.as_deref()) {
        not_an_input(o, &ins)?;
    }
    let t = read_feature_table(table)?;
    if opt.rarefy_depth.is_some() {
        check_counts(&t, table, "rarefaction")?;
    }
    let mut samples = Vec::new();
    let mut cols: Vec<Vec<f64>> = Vec::new();
    let mut dropped = Vec::new();
    for (s, name) in t.samples.iter().enumerate() {
        let mut c = t.sample_column(s);
        let total: f64 = c.iter().sum();
        if total <= 0.0 {
            dropped.push(name.clone());
            continue;
        }
        if let Some(d) = opt.rarefy_depth {
            if (total as u64) < d {
                dropped.push(name.clone());
                continue;
            }
            let r = rarefy(&t.sample_counts(s), d, &mut Rng::new(seed_for(opt.seed, name.as_bytes())));
            c = r.iter().map(|&x| x as f64).collect();
        }
        if opt.relative || metric == "kmer" {
            let tot: f64 = c.iter().sum();
            for x in c.iter_mut() {
                *x /= tot;
            }
        }
        samples.push(name.clone());
        cols.push(c);
    }
    if samples.len() < 2 {
        return Err(format!("{table}: fewer than two samples with reads{}. Beta diversity compares samples", if dropped.is_empty() { String::new() } else { format!(" (left out: {})", dropped.join(", ")) }));
    }
    if metric == "kmer" {
        let fa = opt.asv_fasta.as_deref().ok_or("The kmer distance needs the ASV sequences: give asv_fasta (the FASTA written with the feature table)")?;
        let seqs: HashMap<String, Vec<u8>> = read_fasta(fa)?.into_iter().collect();
        let missing: Vec<&String> = t.features.iter().filter(|f| !seqs.contains_key(*f)).collect();
        if !missing.is_empty() {
            return Err(format!(
                "{} features of {table} are not in {fa} (for instance {}). Give the FASTA written together with this table",
                missing.len(),
                missing[0]
            ));
        }
        let prof: Vec<Vec<f64>> = t.features.par_iter().map(|f| kmer_profile(&seqs[f])).collect();
        cols = cols
            .iter()
            .map(|c| {
                let mut v = vec![0.0; 4096];
                for (f, &x) in c.iter().enumerate() {
                    if x > 0.0 {
                        for (k, p) in prof[f].iter().enumerate() {
                            v[k] += x * p;
                        }
                    }
                }
                v
            })
            .collect();
    }
    let n = samples.len();
    let dist_fn = |a: &[f64], b: &[f64]| if metric == "jaccard" { jaccard(a, b) } else { bray_curtis(a, b) };
    let pairs: Vec<(usize, usize)> = (0..n).flat_map(|i| (i + 1..n).map(move |j| (i, j))).collect();
    let vals: Vec<f64> = pairs.par_iter().map(|&(i, j)| dist_fn(&cols[i], &cols[j])).collect();
    let mut d = vec![vec![0.0; n]; n];
    for (&(i, j), &v) in pairs.iter().zip(&vals) {
        d[i][j] = v;
        d[j][i] = v;
    }
    let pc = pcoa(&d);
    let permanova_res = match &opt.metadata {
        None => None,
        Some(m) => {
            let (column, values) = read_metadata_column(m, opt.group_column.as_deref())?;
            let mut keep = Vec::new();
            let mut excluded = Vec::new();
            let mut labels: Vec<String> = Vec::new();
            for (i, s) in samples.iter().enumerate() {
                match values.get(s) {
                    Some(v) => {
                        keep.push(i);
                        labels.push(v.clone());
                    }
                    None => excluded.push(s.clone()),
                }
            }
            let mut groups: Vec<String> = labels.clone();
            groups.sort();
            groups.dedup();
            if groups.len() < 2 || keep.len() <= groups.len() {
                return Err(format!(
                    "PERMANOVA needs at least two groups and more samples than groups: the column \"{column}\" of {m} gives {} group(s) for {} samples of the table. Check that the sample ids of the metadata match the table's",
                    groups.len(),
                    keep.len()
                ));
            }
            let gi: Vec<usize> = labels.iter().map(|l| groups.iter().position(|g| g == l).unwrap()).collect();
            let sub: Vec<Vec<f64>> = keep.iter().map(|&i| keep.iter().map(|&j| d[i][j]).collect()).collect();
            let (f, r2, p) = permanova(&sub, &gi, opt.permutations, opt.seed);
            let sizes = (0..groups.len()).map(|g| gi.iter().filter(|&&x| x == g).count()).collect();
            Some(Permanova { column, groups, sizes, n: keep.len(), f, r2, p_value: p, permutations: opt.permutations, excluded })
        }
    };
    let mut cleanup = Cleanup(Vec::new());
    if let Some(o) = output {
        let mut s = String::new();
        for x in &samples {
            s.push('\t');
            s.push_str(x);
        }
        s.push('\n');
        for (i, x) in samples.iter().enumerate() {
            s.push_str(x);
            for v in &d[i] {
                s.push('\t');
                s.push_str(&num(*v));
            }
            s.push('\n');
        }
        cleanup.0.push(o.to_string());
        write_text(o, &s)?;
    }
    if let Some(o) = &opt.pcoa_output {
        let axes = pc.coordinates.first().map_or(0, Vec::len);
        let mut s = String::from("sample");
        for k in 0..axes {
            s.push_str(&format!("\tPC{}", k + 1));
        }
        s.push('\n');
        s.push_str("proportion_explained");
        for k in 0..axes {
            s.push_str(&format!("\t{}", num(pc.proportion_explained[k])));
        }
        s.push('\n');
        for (i, x) in samples.iter().enumerate() {
            s.push_str(x);
            for v in &pc.coordinates[i] {
                s.push_str(&format!("\t{}", num(*v)));
            }
            s.push('\n');
        }
        cleanup.0.push(o.clone());
        write_text(o, &s)?;
    }
    cleanup.disarm();
    Ok(BetaResult {
        metric: metric.to_string(),
        samples,
        distances: d,
        pcoa: pc,
        permanova: permanova_res,
        dropped,
        rarefied_to: opt.rarefy_depth,
        output: output.map(str::to_string),
        pcoa_output: opt.pcoa_output.clone(),
    })
}

impl BetaResult {
    pub fn json(&self) -> String {
        let d: Vec<String> = self.distances.iter().map(|r| jnums(r)).collect();
        let c: Vec<String> = self.pcoa.coordinates.iter().map(|r| jnums(r)).collect();
        let p = match &self.permanova {
            None => "null".to_string(),
            Some(p) => format!(
                "{{\"column\":{},\"groups\":{},\"sizes\":[{}],\"n\":{},\"f\":{},\"r2\":{},\"p_value\":{},\"permutations\":{},\"excluded\":{}}}",
                js(&p.column),
                jlist(&p.groups),
                p.sizes.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","),
                p.n,
                jf(p.f),
                jf(p.r2),
                jf(p.p_value),
                p.permutations,
                jlist(&p.excluded)
            ),
        };
        format!(
            "{{\"metric\":{},\"samples\":{},\"distances\":[{}],\"pcoa\":{{\"eigenvalues\":{},\"proportion_explained\":{},\"coordinates\":[{}],\"negative_eigenvalues\":{}}},\"permanova\":{},\"dropped\":{},\"rarefied_to\":{},\"output\":{},\"pcoa_output\":{}}}",
            js(&self.metric),
            jlist(&self.samples),
            d.join(","),
            jnums(&self.pcoa.eigenvalues),
            jnums(&self.pcoa.proportion_explained),
            c.join(","),
            self.pcoa.negative_eigenvalues,
            p,
            jlist(&self.dropped),
            self.rarefied_to.map_or("null".to_string(), |x| x.to_string()),
            self.output.as_deref().map_or("null".to_string(), js),
            self.pcoa_output.as_deref().map_or("null".to_string(), js)
        )
    }

    pub fn resumo(&self) -> String {
        let pe = &self.pcoa.proportion_explained;
        let mut o = format!("Beta diversity ({}) between {} samples.", self.metric, self.samples.len());
        if pe.len() >= 2 {
            o.push_str(&format!(" PCoA: PC1 explains {:.1}%, PC2 {:.1}%", 100.0 * pe[0], 100.0 * pe[1]));
            if pe.len() >= 3 {
                o.push_str(&format!(", PC3 {:.1}%", 100.0 * pe[2]));
            }
            o.push('.');
        }
        if self.pcoa.negative_eigenvalues > 0 {
            o.push_str(&format!(" {} negative eigenvalue(s): the distance is not Euclidean, as usual for Bray-Curtis.", self.pcoa.negative_eigenvalues));
        }
        o.push('\n');
        if let Some(p) = &self.permanova {
            let g: Vec<String> = p.groups.iter().zip(&p.sizes).map(|(g, n)| format!("{g} ({n})")).collect();
            o.push_str(&format!(
                "PERMANOVA by {}: {}; pseudo-F = {:.3}, R2 = {:.3}, p = {} ({} permutations).\n",
                p.column,
                g.join(" vs "),
                p.f,
                p.r2,
                if p.p_value.is_nan() { "NA".to_string() } else { format!("{:.4}", p.p_value) },
                p.permutations
            ));
            if !p.excluded.is_empty() {
                o.push_str(&format!("Not in the metadata (left out of PERMANOVA): {}\n", p.excluded.join(", ")));
            }
        }
        if !self.dropped.is_empty() {
            o.push_str(&format!("Left out: {}\n", self.dropped.join(", ")));
        }
        o
    }
}

// ---------------------------------------------------------------------
//  TAXONOMIC COMPOSITION
// ---------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct CompositionOptions {
    /// The rank to sum at (a column of the taxonomy table: "Phylum",
    /// "Family", "Genus"...).
    pub rank: String,
    /// Keep the `top` most abundant taxa (by mean fraction) and sum the rest
    /// as "Other"; 0 = all.
    pub top: usize,
}

impl Default for CompositionOptions {
    fn default() -> Self {
        CompositionOptions { rank: "Phylum".to_string(), top: 0 }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CompositionResult {
    pub rank: String,
    pub samples: Vec<String>,
    /// From the most abundant (mean fraction); "Other" last when `top` cut.
    pub taxa: Vec<String>,
    /// fractions[taxon][sample]
    pub fractions: Vec<Vec<f64>>,
    /// Features of the table missing from the taxonomy (they count as
    /// "Unassigned").
    pub features_without_taxonomy: Vec<String>,
    pub output: Option<String>,
}

/// Relative abundance of each taxon at a rank, per sample (for stacked
/// bars). A feature unassigned at that rank counts as "<deepest name>
/// (unclassified)", or "Unassigned".
pub fn taxa_composition(table: &str, taxonomy: &str, output: Option<&str>, opt: &CompositionOptions) -> Result<CompositionResult, String> {
    if let Some(o) = output {
        not_an_input(o, &[table, taxonomy])?;
    }
    let t = read_feature_table(table)?;
    check_input(taxonomy)?;
    let mut text = String::new();
    abre_leitura(taxonomy)?.read_to_string(&mut text).map_err(|e| format!("Error reading the taxonomy table {taxonomy}: {e}"))?;
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header: Vec<&str> = lines.next().ok_or_else(|| format!("{taxonomy} is empty: give the taxonomy table written by classify_asvs"))?.split('\t').map(str::trim).collect();
    let rank_cols: Vec<usize> = (1..header.len()).filter(|&j| !header[j].to_ascii_lowercase().ends_with("_confidence")).collect();
    let target = rank_cols.iter().position(|&j| header[j].eq_ignore_ascii_case(opt.rank.trim())).ok_or_else(|| {
        format!(
            "{taxonomy}: no rank named \"{}\". The ranks are: {}",
            opt.rank,
            rank_cols.iter().map(|&j| header[j]).collect::<Vec<_>>().join(", ")
        )
    })?;
    let mut label_of: HashMap<String, String> = HashMap::new();
    for l in lines {
        let f: Vec<&str> = l.split('\t').map(str::trim).collect();
        let get = |j: usize| f.get(j).copied().unwrap_or("");
        let label = {
            let v = get(rank_cols[target]);
            if !v.is_empty() {
                v.to_string()
            } else {
                match rank_cols[..target].iter().rev().map(|&j| get(j)).find(|x| !x.is_empty()) {
                    Some(up) => format!("{up} (unclassified)"),
                    None => "Unassigned".to_string(),
                }
            }
        };
        label_of.insert(f[0].to_string(), label);
    }
    let ns = t.samples.len();
    let mut idx: HashMap<String, usize> = HashMap::new();
    let mut taxa: Vec<String> = Vec::new();
    let mut sums: Vec<Vec<f64>> = Vec::new();
    let mut missing = Vec::new();
    for (fi, f) in t.features.iter().enumerate() {
        let label = match label_of.get(f) {
            Some(l) => l.clone(),
            None => {
                missing.push(f.clone());
                "Unassigned".to_string()
            }
        };
        let k = *idx.entry(label.clone()).or_insert_with(|| {
            taxa.push(label);
            sums.push(vec![0.0; ns]);
            taxa.len() - 1
        });
        for s in 0..ns {
            sums[k][s] += t.counts[fi][s];
        }
    }
    if !t.features.is_empty() && missing.len() == t.features.len() {
        return Err(format!(
            "None of the features of {table} is in {taxonomy}: the names do not match (for instance {}). Give the taxonomy table made from this table's ASVs",
            t.features[0]
        ));
    }
    let totals: Vec<f64> = (0..ns).map(|s| t.counts.iter().map(|r| r[s]).sum()).collect();
    let fr: Vec<Vec<f64>> = sums.iter().map(|r| r.iter().zip(&totals).map(|(x, tot)| if *tot > 0.0 { x / tot } else { 0.0 }).collect()).collect();
    let mean: Vec<f64> = fr.iter().map(|r| r.iter().sum::<f64>() / ns.max(1) as f64).collect();
    let order = sorted_indices(taxa.len(), |&x, &y| {
        mean[y as usize].partial_cmp(&mean[x as usize]).unwrap_or(Ordering::Equal).then_with(|| taxa[x as usize].cmp(&taxa[y as usize]))
    });
    let mut out_taxa = Vec::new();
    let mut out_fr = Vec::new();
    let mut other = vec![0.0; ns];
    let mut has_other = false;
    for (rank, &k) in order.iter().enumerate() {
        if opt.top > 0 && rank >= opt.top {
            has_other = true;
            for s in 0..ns {
                other[s] += fr[k as usize][s];
            }
        } else {
            out_taxa.push(taxa[k as usize].clone());
            out_fr.push(fr[k as usize].clone());
        }
    }
    if has_other {
        out_taxa.push("Other".to_string());
        out_fr.push(other);
    }
    if let Some(o) = output {
        let mut s = String::from("taxon");
        for x in &t.samples {
            s.push('\t');
            s.push_str(x);
        }
        s.push('\n');
        for (name, r) in out_taxa.iter().zip(&out_fr) {
            s.push_str(name);
            for v in r {
                s.push('\t');
                s.push_str(&num(*v));
            }
            s.push('\n');
        }
        write_text(o, &s)?;
    }
    Ok(CompositionResult {
        rank: header[rank_cols[target]].to_string(),
        samples: t.samples,
        taxa: out_taxa,
        fractions: out_fr,
        features_without_taxonomy: missing,
        output: output.map(str::to_string),
    })
}

impl CompositionResult {
    pub fn json(&self) -> String {
        let f: Vec<String> = self.fractions.iter().map(|r| jnums(r)).collect();
        format!(
            "{{\"rank\":{},\"samples\":{},\"taxa\":{},\"fractions\":[{}],\"features_without_taxonomy\":{},\"output\":{}}}",
            js(&self.rank),
            jlist(&self.samples),
            jlist(&self.taxa),
            f.join(","),
            jlist(&self.features_without_taxonomy),
            self.output.as_deref().map_or("null".to_string(), js)
        )
    }

    pub fn resumo(&self) -> String {
        let ns = self.samples.len().max(1) as f64;
        let mut o = format!("Composition at {} rank: {} taxa over {} samples. Mean relative abundance:\n", self.rank, self.taxa.len(), self.samples.len());
        for (name, r) in self.taxa.iter().zip(&self.fractions).take(15) {
            o.push_str(&format!("  {:>6.2}%  {}\n", 100.0 * r.iter().sum::<f64>() / ns, name));
        }
        if self.taxa.len() > 15 {
            o.push_str(&format!("  ... and {} more\n", self.taxa.len() - 15));
        }
        if !self.features_without_taxonomy.is_empty() {
            o.push_str(&format!("warning: {} features have no line in the taxonomy table (counted as Unassigned)\n", self.features_without_taxonomy.len()));
        }
        o
    }
}

// ---------------------------------------------------------------------
//  TUMOR MUTATIONAL BURDEN
// ---------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct TmbOptions {
    /// The tumor column of the VCF; None = `##tumor_sample=` of the header,
    /// a column named TUMOR, or the only sample.
    pub tumor_sample: Option<String>,
    /// Count only FILTER = PASS (true).
    pub pass_only: bool,
    /// Lowest tumor depth (10) and variant allele fraction (0.05).
    pub min_depth: u32,
    pub min_vaf: f64,
    /// GFF3 or GTF with the CDS: needed by coding_only and nonsynonymous_only.
    pub annotation: Option<String>,
    pub coding_only: bool,
    /// Reference genome (FASTA): needed by nonsynonymous_only.
    pub reference: Option<String>,
    pub nonsynonymous_only: bool,
    /// TSV of every ALT allele: counted or why not.
    pub output: Option<String>,
}

impl Default for TmbOptions {
    fn default() -> Self {
        TmbOptions {
            tumor_sample: None,
            pass_only: true,
            min_depth: 10,
            min_vaf: 0.05,
            annotation: None,
            coding_only: false,
            reference: None,
            nonsynonymous_only: false,
            output: None,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TmbResult {
    pub vcf: String,
    pub bed: String,
    pub sample: Option<String>,
    /// Bases of the BED (merged), of the CDS inside it, and the denominator.
    pub bed_bp: u64,
    pub coding_bp: Option<u64>,
    pub region_bp: u64,
    pub megabases: f64,
    pub records: u64,
    pub alleles: u64,
    pub counted: u64,
    pub counted_snv: u64,
    pub counted_mnv: u64,
    pub counted_indel: u64,
    /// Counted per Mb: all, substitutions (SNV + MNV), indels.
    pub tmb: f64,
    pub tmb_snv: f64,
    pub tmb_indel: f64,
    pub excluded_filter: u64,
    pub excluded_depth: u64,
    pub excluded_vaf: u64,
    pub outside_region: u64,
    pub not_coding: u64,
    pub synonymous: u64,
    pub symbolic: u64,
    /// Kept although the depth / VAF could not be read.
    pub without_depth: u64,
    pub without_vaf: u64,
    /// SNV/MNV whose REF differs from the reference genome.
    pub ref_mismatch: u64,
    /// Coding variants whose codon could not be completed (counted).
    pub coding_unknown: u64,
    pub output: Option<String>,
    pub warnings: Vec<String>,
}

/// The name of `name` in a set of known names: itself, or in the other
/// convention (`chr1` / `1`, `chrM` / `MT`).
fn resolve_name(name: &str, known: &HashSet<String>) -> Option<String> {
    use crate::nomes::{converte_nome, Estilo};
    if known.contains(name) {
        return Some(name.to_string());
    }
    for e in [Estilo::ComChr, Estilo::SemChr] {
        if let Some(n) = converte_nome(name, e) {
            if known.contains(&n) {
                return Some(n);
            }
        }
    }
    if name == "chrM" && known.contains("chrMT") {
        return Some("chrMT".to_string());
    }
    None
}

/// Sorts and merges intervals (0-based, half-open).
fn merge_intervals(mut v: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    multi_merge_sort_stable_by(&mut v, |a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut out: Vec<(u64, u64)> = Vec::with_capacity(v.len());
    for (s, e) in v {
        match out.last_mut() {
            Some(last) if s <= last.1 => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

fn overlaps(iv: &[(u64, u64)], s: u64, e: u64) -> bool {
    // first interval ending after s
    let i = iv.partition_point(|x| x.1 <= s);
    i < iv.len() && iv[i].0 < e
}

fn intersect(a: &[(u64, u64)], b: &[(u64, u64)]) -> Vec<(u64, u64)> {
    let (mut i, mut j) = (0, 0);
    let mut out = Vec::new();
    while i < a.len() && j < b.len() {
        let s = a[i].0.max(b[j].0);
        let e = a[i].1.min(b[j].1);
        if s < e {
            out.push((s, e));
        }
        if a[i].1 < b[j].1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    out
}

/// A BED file as merged intervals per sequence.
fn read_bed_region(path: &str) -> Result<HashMap<String, Vec<(u64, u64)>>, String> {
    check_input(path)?;
    let r = abre_leitura(path)?;
    let mut raw: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    for (i, line) in r.lines().enumerate() {
        let line = line.map_err(|e| format!("Error reading {path} near line {}: {e}", i + 1))?;
        let l = line.trim_end_matches('\r');
        if l.trim().is_empty() || l.starts_with('#') || l.starts_with("track") || l.starts_with("browser") {
            continue;
        }
        let f: Vec<&str> = l.split('\t').collect();
        let (s, e) = match (f.get(1).and_then(|x| x.trim().parse::<u64>().ok()), f.get(2).and_then(|x| x.trim().parse::<u64>().ok())) {
            (Some(s), Some(e)) if f.len() >= 3 && e >= s => (s, e),
            _ => {
                return Err(format!(
                    "{path}, line {}: not a BED line (sequence, start, end, tab-separated, start <= end, 0-based): \"{}\". Fix that line",
                    i + 1,
                    sample_text(l.as_bytes())
                ))
            }
        };
        raw.entry(f[0].to_string()).or_default().push((s, e));
    }
    if raw.is_empty() {
        return Err(format!("{path} has no interval: give the BED of the captured (or callable) region"));
    }
    Ok(raw.into_iter().map(|(k, v)| (k, merge_intervals(v))).collect())
}

#[derive(Clone, Debug)]
struct Transcript {
    chrom: String,
    minus: bool,
    /// CDS pieces (0-based, half-open) in transcript order.
    segs: Vec<(u64, u64)>,
    phase0: u64,
    len: u64,
    gene: String,
}

impl Transcript {
    fn tx_pos(&self, g: u64) -> Option<u64> {
        let mut acc = 0;
        for &(s, e) in &self.segs {
            if g >= s && g < e {
                return Some(acc + if self.minus { e - 1 - g } else { g - s });
            }
            acc += e - s;
        }
        None
    }
    fn genome_pos(&self, c: u64) -> Option<u64> {
        let mut acc = 0;
        for &(s, e) in &self.segs {
            let l = e - s;
            if c < acc + l {
                let o = c - acc;
                return Some(if self.minus { e - 1 - o } else { s + o });
            }
            acc += l;
        }
        None
    }
}

struct Annotation {
    txs: Vec<Transcript>,
    /// Per sequence: (start, end, transcript) sorted by start, and the
    /// longest piece.
    index: HashMap<String, (Vec<(u64, u64, u32)>, u64)>,
    /// Merged CDS per sequence.
    cds: HashMap<String, Vec<(u64, u64)>>,
}

fn attr_value(attrs: &str, keys: &[&str]) -> Option<String> {
    for part in attrs.split(';') {
        let p = part.trim();
        if p.is_empty() {
            continue;
        }
        let (k, v) = if let Some(eq) = p.find('=') {
            (&p[..eq], p[eq + 1..].trim())
        } else if let Some(sp) = p.find(' ') {
            (&p[..sp], p[sp + 1..].trim().trim_matches('"'))
        } else {
            continue;
        };
        if keys.iter().any(|x| x.eq_ignore_ascii_case(k.trim())) && !v.is_empty() {
            return Some(v.to_string());
        }
    }
    None
}

/// The CDS of a GFF3 or GTF, on the sequences of the region.
fn read_annotation(path: &str, region: &HashMap<String, Vec<(u64, u64)>>) -> Result<Annotation, String> {
    check_input(path)?;
    let known: HashSet<String> = region.keys().cloned().collect();
    let r = abre_leitura(path)?;
    let mut tx_of: HashMap<(String, String), usize> = HashMap::new();
    let mut txs: Vec<Transcript> = Vec::new();
    let mut phases: Vec<Vec<(u64, u64, u64)>> = Vec::new(); // per tx: (start, end, phase)
    let mut name_cache: HashMap<String, Option<String>> = HashMap::new();
    // RefSeq names its sequences by accession (NC_000001.11); its `region`
    // lines say which chromosome each one is
    let mut alias: HashMap<String, String> = HashMap::new();
    let mut cds_lines = 0u64;
    for (i, line) in r.lines().enumerate() {
        let line = line.map_err(|e| format!("Error reading {path} near line {}: {e}", i + 1))?;
        if line.starts_with("##FASTA") {
            break;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 {
            continue;
        }
        if f[2] == "region" && !known.contains(f[0]) {
            match attr_value(f[8], &["genome"]).as_deref() {
                Some("chromosome") => {
                    if let Some(c) = attr_value(f[8], &["chromosome"]) {
                        alias.insert(f[0].to_string(), c);
                    }
                }
                Some("mitochondrion") => {
                    alias.insert(f[0].to_string(), "MT".to_string());
                }
                _ => {}
            }
            continue;
        }
        if f[2] != "CDS" {
            continue;
        }
        cds_lines += 1;
        let chrom = match name_cache.get(f[0]) {
            Some(c) => c.clone(),
            None => {
                let c = resolve_name(f[0], &known).or_else(|| alias.get(f[0]).and_then(|a| resolve_name(a, &known)));
                name_cache.insert(f[0].to_string(), c.clone());
                c
            }
        };
        let Some(chrom) = chrom else {
            continue;
        };
        let (Ok(s), Ok(e)) = (f[3].trim().parse::<u64>(), f[4].trim().parse::<u64>()) else {
            return Err(format!("{path}, line {}: the start or end of a CDS is not a number. Fix that line", i + 1));
        };
        if s == 0 || e < s {
            return Err(format!("{path}, line {}: CDS from {s} to {e}: GFF coordinates start at 1 and the end is not before the start. Fix that line", i + 1));
        }
        let phase: u64 = f[7].trim().parse().unwrap_or(0) % 3;
        let minus = f[6].trim() == "-";
        let parents: Vec<String> = match attr_value(f[8], &["Parent", "transcript_id"]) {
            Some(p) => p.split(',').map(|x| x.trim().to_string()).collect(),
            None => vec![attr_value(f[8], &["ID"]).unwrap_or_else(|| format!("line{}", i + 1))],
        };
        let gene = attr_value(f[8], &["gene", "gene_name", "Name", "gene_id"]).unwrap_or_else(|| parents[0].clone());
        for p in parents {
            let key = (chrom.clone(), p);
            let t = match tx_of.get(&key) {
                Some(&t) => t,
                None => {
                    txs.push(Transcript { chrom: chrom.clone(), minus, segs: Vec::new(), phase0: 0, len: 0, gene: gene.clone() });
                    phases.push(Vec::new());
                    tx_of.insert(key, txs.len() - 1);
                    txs.len() - 1
                }
            };
            phases[t].push((s - 1, e, phase));
        }
    }
    if cds_lines == 0 {
        return Err(format!("{path}: no CDS line. Give a GFF3 or GTF annotation with the coding sequences (feature type CDS)"));
    }
    let mut index: HashMap<String, (Vec<(u64, u64, u32)>, u64)> = HashMap::new();
    let mut cds_raw: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    for (t, mut pieces) in phases.into_iter().enumerate() {
        multi_merge_sort_stable_by(&mut pieces, |a, b| a.0.cmp(&b.0));
        let tx = &mut txs[t];
        if tx.minus {
            pieces.reverse();
        }
        tx.phase0 = pieces.first().map_or(0, |p| p.2);
        tx.segs = pieces.iter().map(|p| (p.0, p.1)).collect();
        tx.len = tx.segs.iter().map(|(s, e)| e - s).sum();
        let ent = index.entry(tx.chrom.clone()).or_insert_with(|| (Vec::new(), 0));
        for &(s, e) in &tx.segs {
            ent.0.push((s, e, t as u32));
            ent.1 = ent.1.max(e - s);
            cds_raw.entry(tx.chrom.clone()).or_default().push((s, e));
        }
    }
    for v in index.values_mut() {
        multi_merge_sort_stable_by(&mut v.0, |a, b| a.0.cmp(&b.0));
    }
    let cds = cds_raw.into_iter().map(|(k, v)| (k, merge_intervals(v))).collect();
    Ok(Annotation { txs, index, cds })
}

impl Annotation {
    /// Transcripts with a CDS piece over [s, e).
    fn transcripts_at(&self, chrom: &str, s: u64, e: u64) -> Vec<usize> {
        let mut out = Vec::new();
        if let Some((v, maxlen)) = self.index.get(chrom) {
            let mut i = v.partition_point(|x| x.0 < e);
            while i > 0 {
                i -= 1;
                let (a, b, t) = v[i];
                if a + maxlen < s {
                    break;
                }
                if a < e && b > s && !out.contains(&(t as usize)) {
                    out.push(t as usize);
                }
            }
        }
        out.sort_unstable();
        out
    }
}

const GENETIC_CODE: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

/// The amino acid of a codon (standard code); 'X' with a non-ACGT base.
pub fn translate_codon(c: &[u8]) -> u8 {
    let idx = |b: u8| match b {
        b'T' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    };
    match (idx(c[0]), idx(c[1]), idx(c[2])) {
        (Some(a), Some(b), Some(d)) => GENETIC_CODE[a * 16 + b * 4 + d],
        _ => b'X',
    }
}

fn severity(c: &str) -> u8 {
    match c {
        "stop_gained" => 6,
        "frameshift" => 6,
        "stop_lost" => 5,
        "start_lost" => 4,
        "inframe_indel" => 3,
        "missense" => 2,
        "coding_unknown" => 1,
        _ => 0, // synonymous
    }
}

/// The protein effect of base substitutions [(position, ref, alt)] (plus
/// strand) on a transcript; None when no codon can be completed.
fn codon_effect(
    tx: &Transcript,
    subs: &[(u64, u8, u8)],
    genome: &mut crate::referencia::Genoma,
    gname: &str,
) -> Result<(Option<&'static str>, bool), String> {
    let mut codons: Vec<(u64, Vec<(u64, u8, u8)>)> = Vec::new();
    for &(g, r, a) in subs {
        let Some(c) = tx.tx_pos(g) else {
            continue;
        };
        if c < tx.phase0 {
            continue;
        }
        let cs = c - (c - tx.phase0) % 3;
        if cs + 3 > tx.len {
            continue;
        }
        match codons.iter_mut().find(|x| x.0 == cs) {
            Some(x) => x.1.push((c - cs, r, a)),
            None => codons.push((cs, vec![(c - cs, r, a)])),
        }
    }
    if codons.is_empty() {
        return Ok((None, false));
    }
    let mut worst: Option<&'static str> = None;
    let mut mismatch = false;
    for (cs, changes) in &codons {
        let mut rc = [b'N'; 3];
        for (k, slot) in rc.iter_mut().enumerate() {
            let g = tx.genome_pos(cs + k as u64).ok_or_else(|| "internal error: a codon outside its transcript".to_string())?;
            let b = genome.trecho(gname, g, g + 1)?.first().copied().unwrap_or(b'N');
            *slot = if tx.minus { complement(b) } else { b };
        }
        let mut ac = rc;
        for &(off, r, a) in changes {
            let (r, a) = if tx.minus { (complement(r), complement(a)) } else { (r, a) };
            if rc[off as usize] != r {
                mismatch = true;
                rc[off as usize] = r;
            }
            ac[off as usize] = a;
        }
        let (ra, aa) = (translate_codon(&rc), translate_codon(&ac));
        let eff = if ra == aa {
            "synonymous"
        } else if aa == b'*' {
            "stop_gained"
        } else if ra == b'*' {
            "stop_lost"
        } else if *cs == tx.phase0 && ra == b'M' {
            "start_lost"
        } else {
            "missense"
        };
        if worst.map_or(true, |w| severity(eff) > severity(w)) {
            worst = Some(eff);
        }
    }
    Ok((worst, mismatch))
}

/// Tumor depth and VAF of ALT allele `a` (1-based).
fn depth_and_vaf(keys: &[&str], vals: &[&str], info: &str, a: usize, r: &str, alt: &str) -> (Option<u64>, Option<f64>) {
    let get = |k: &str| keys.iter().position(|x| *x == k).and_then(|i| vals.get(i)).filter(|v| **v != "." && !v.is_empty()).copied();
    let nums = |v: &str| -> Vec<Option<f64>> { v.split(',').map(|x| x.trim().parse::<f64>().ok()).collect() };
    let info_get = |k: &str| {
        info.split(';').find_map(|p| {
            let mut it = p.splitn(2, '=');
            if it.next() == Some(k) {
                it.next()
            } else {
                None
            }
        })
    };
    let mut depth = get("DP").and_then(|v| v.parse::<f64>().ok()).map(|x| x as u64);
    let mut vaf = None;
    if let Some(v) = get("AD") {
        let ad = nums(v);
        let tot: f64 = ad.iter().flatten().sum();
        if depth.is_none() && tot > 0.0 {
            depth = Some(tot as u64);
        }
        if let Some(Some(x)) = ad.get(a) {
            if tot > 0.0 {
                vaf = Some(x / tot);
            }
        }
    }
    if let Some(v) = get("AF") {
        if let Some(Some(x)) = nums(v).get(a - 1) {
            vaf = Some(*x);
        }
    }
    if vaf.is_none() && r.len() == 1 && alt.len() == 1 {
        // Strelka SNVs: AU, CU, GU, TU (tier 1 first)
        let tier = |b: &str| get(&format!("{b}U")).and_then(|v| v.split(',').next().and_then(|x| x.parse::<f64>().ok()));
        if let (Some(rc), Some(ac)) = (tier(r), tier(alt)) {
            if rc + ac > 0.0 {
                vaf = Some(ac / (rc + ac));
            }
        }
    }
    if vaf.is_none() {
        // Strelka indels: TAR (ref) and TIR (alt), tier 1 first
        let t1 = |k: &str| get(k).and_then(|v| v.split(',').next().and_then(|x| x.parse::<f64>().ok()));
        if let (Some(rc), Some(ac)) = (t1("TAR"), t1("TIR")) {
            if rc + ac > 0.0 {
                vaf = Some(ac / (rc + ac));
            }
        }
    }
    if vaf.is_none() {
        if let Some(v) = info_get("AF") {
            if let Some(Some(x)) = nums(v).get(a - 1) {
                vaf = Some(*x);
            }
        }
    }
    if depth.is_none() {
        depth = info_get("DP").and_then(|v| v.parse::<f64>().ok()).map(|x| x as u64);
    }
    (depth, vaf)
}

/// Common suffix, then common prefix, removed (keeping one base each).
fn trim_allele(pos: u64, r: &str, a: &str) -> (u64, String, String) {
    let (mut r, mut a) = (r.as_bytes().to_vec(), a.as_bytes().to_vec());
    while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
        r.pop();
        a.pop();
    }
    let mut p = pos;
    while r.len() > 1 && a.len() > 1 && r[0] == a[0] {
        r.remove(0);
        a.remove(0);
        p += 1;
    }
    (p, String::from_utf8_lossy(&r).to_string(), String::from_utf8_lossy(&a).to_string())
}

/// Tumor mutational burden of a somatic VCF over the region of a BED.
pub fn tumor_mutational_burden(vcf: &str, bed: &str, opt: &TmbOptions) -> Result<TmbResult, String> {
    let mut ins = vec![vcf, bed];
    if let Some(a) = &opt.annotation {
        ins.push(a);
    }
    if let Some(r) = &opt.reference {
        ins.push(r);
    }
    if let Some(o) = &opt.output {
        not_an_input(o, &ins)?;
    }
    if (opt.coding_only || opt.nonsynonymous_only) && opt.annotation.is_none() {
        return Err("coding_only and nonsynonymous_only need the annotation (a GFF3 or GTF with the CDS): give `annotation`".to_string());
    }
    if opt.nonsynonymous_only && opt.reference.is_none() {
        return Err("nonsynonymous_only needs the reference genome (FASTA) to translate the codons: give `reference`, or count coding variants only (coding_only)".to_string());
    }
    if !(0.0..=1.0).contains(&opt.min_vaf) {
        return Err(format!("min_vaf = {}: give a fraction from 0 to 1 (0.05 = 5%)", opt.min_vaf));
    }
    check_input(vcf)?;
    let region = read_bed_region(bed)?;
    let known: HashSet<String> = region.keys().cloned().collect();
    let bed_bp: u64 = region.values().flatten().map(|(s, e)| e - s).sum();
    let ann = match &opt.annotation {
        Some(a) => Some(read_annotation(a, &region)?),
        None => None,
    };
    let mut genome = match (&opt.reference, ann.is_some()) {
        (Some(r), true) => Some(crate::referencia::Genoma::abre(r, None)?),
        _ => None,
    };
    let coding_region: Option<HashMap<String, Vec<(u64, u64)>>> = ann.as_ref().map(|a| {
        region
            .iter()
            .map(|(k, v)| (k.clone(), a.cds.get(k).map(|c| intersect(v, c)).unwrap_or_default()))
            .collect()
    });
    let coding_bp = coding_region.as_ref().map(|c| c.values().flatten().map(|(s, e)| e - s).sum::<u64>());
    let use_coding = opt.coding_only || opt.nonsynonymous_only;
    let region_bp = if use_coding { coding_bp.unwrap_or(0) } else { bed_bp };
    let mut res = TmbResult {
        vcf: vcf.to_string(),
        bed: bed.to_string(),
        bed_bp,
        coding_bp,
        region_bp,
        megabases: region_bp as f64 / 1e6,
        output: opt.output.clone(),
        ..Default::default()
    };
    if region_bp == 0 {
        return Err(format!(
            "The region has 0 bases{}: no TMB can be computed. Check that the BED and the annotation use the same genome and sequence names",
            if use_coding { " of CDS inside the BED" } else { "" }
        ));
    }
    // declared before the writer, so that on an error the writer is closed
    // first and the half-written file can then be removed
    let mut cleanup = Cleanup(opt.output.iter().cloned().collect());
    let mut out = match &opt.output {
        Some(o) => {
            let mut w = create_output(o)?;
            w.linha("chrom\tpos\tref\talt\ttype\tfilter\tdepth\tvaf\tin_region\tcoding\tconsequence\tgene\tcounted\treason").map_err(|e| write_error(o, e))?;
            Some(w)
        }
        None => None,
    };
    let r = abre_leitura(vcf)?;
    let mut tumor_col: Option<usize> = None;
    let mut header_tumor: Option<String> = None;
    let mut seen_header = false;
    let mut name_cache: HashMap<String, Option<String>> = HashMap::new();
    let mut vcf_chroms_matched = false;
    let mut vcf_chroms: Vec<String> = Vec::new();
    let mut all_dot_filter = true;
    for (i, line) in r.lines().enumerate() {
        let line = line.map_err(|e| format!("Error reading {vcf} near line {}: {e}. If it is compressed, it may be damaged", i + 1))?;
        let line = line.trim_end_matches('\r');
        if line.starts_with("##") {
            if let Some(v) = line.strip_prefix("##tumor_sample=") {
                header_tumor = Some(v.trim().to_string());
            }
            continue;
        }
        if line.starts_with("#CHROM") {
            let f: Vec<&str> = line.split('\t').collect();
            let names: Vec<String> = f.iter().skip(9).map(|s| s.trim().to_string()).collect();
            tumor_col = if let Some(t) = &opt.tumor_sample {
                Some(names.iter().position(|n| n == t).ok_or_else(|| format!("{vcf}: no sample named {t}. The samples are: {}", names.join(", ")))?)
            } else if let Some(p) = header_tumor.as_ref().and_then(|t| names.iter().position(|n| n == t)) {
                Some(p)
            } else if let Some(p) = names.iter().position(|n| n.eq_ignore_ascii_case("tumor") || n.eq_ignore_ascii_case("tumour")) {
                Some(p)
            } else if names.len() == 1 {
                Some(0)
            } else if names.is_empty() {
                None
            } else {
                return Err(format!(
                    "{vcf} has {} samples ({}) and none is marked as the tumor (no ##tumor_sample= line, no column named TUMOR). Give the tumor sample's name (tumor_sample)",
                    names.len(),
                    names.join(", ")
                ));
            };
            res.sample = tumor_col.map(|c| names[c].clone());
            seen_header = true;
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        if !seen_header {
            return Err(format!("{vcf}: no #CHROM header line before the first record (line {}). Is this a VCF file?", i + 1));
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            return Err(format!("{vcf}, line {}: a VCF record has at least 8 tab-separated columns, this one has {}. Fix that line", i + 1, f.len()));
        }
        let Ok(pos) = f[1].trim().parse::<u64>() else {
            return Err(format!("{vcf}, line {}: POS \"{}\" is not a number. Fix that line", i + 1, f[1]));
        };
        res.records += 1;
        let chrom_raw = f[0];
        let chrom = match name_cache.get(chrom_raw) {
            Some(c) => c.clone(),
            None => {
                let c = resolve_name(chrom_raw, &known);
                if vcf_chroms.len() < 5 && !vcf_chroms.iter().any(|x| x == chrom_raw) {
                    vcf_chroms.push(chrom_raw.to_string());
                }
                name_cache.insert(chrom_raw.to_string(), c.clone());
                c
            }
        };
        if chrom.is_some() {
            vcf_chroms_matched = true;
        }
        let filter = f[6].trim();
        if filter != "." {
            all_dot_filter = false;
        }
        // PASS among the filters (some truth sets write "PASS;HighConf")
        let pass = filter.split(';').any(|x| x.trim() == "PASS");
        let keys: Vec<&str> = f.get(8).map(|x| x.split(':').collect()).unwrap_or_default();
        let vals: Vec<&str> = match tumor_col {
            Some(c) => f.get(9 + c).map(|x| x.split(':').collect()).unwrap_or_default(),
            None => Vec::new(),
        };
        let reference = f[3].trim().to_ascii_uppercase();
        for (ai, alt_raw) in f[4].split(',').enumerate() {
            let alt_raw = alt_raw.trim().to_ascii_uppercase();
            res.alleles += 1;
            let symbolic = alt_raw.is_empty() || alt_raw == "*" || alt_raw == "." || alt_raw.starts_with('<') || alt_raw.contains('[') || alt_raw.contains(']');
            let (p, rf, al) = if symbolic { (pos, reference.clone(), alt_raw.clone()) } else { trim_allele(pos, &reference, &alt_raw) };
            let kind = if symbolic {
                "other"
            } else if rf.len() == 1 && al.len() == 1 {
                "snv"
            } else if rf.len() == al.len() {
                "mnv"
            } else {
                "indel"
            };
            let (depth, vaf) = depth_and_vaf(&keys, &vals, f[7], ai + 1, &rf, &al);
            let p0 = p - 1;
            let span_end = p0 + rf.len().max(1) as u64;
            let in_region = chrom.as_ref().and_then(|c| region.get(c)).map_or(false, |iv| overlaps(iv, p0, span_end));
            let coding = match (&ann, &chrom) {
                (Some(a), Some(c)) => Some(a.cds.get(c).map_or(false, |iv| overlaps(iv, p0, span_end))),
                (Some(_), None) => Some(false),
                _ => None,
            };
            // the consequence, when there is an annotation
            let mut consequence: Option<&'static str> = None;
            let mut gene = String::new();
            if let (Some(a), Some(c), Some(true)) = (&ann, &chrom, coding) {
                let txs = a.transcripts_at(c, p0, span_end);
                if let Some(&t0) = txs.first() {
                    gene = a.txs[t0].gene.clone();
                }
                if kind == "indel" {
                    consequence = Some(if (al.len() as i64 - rf.len() as i64) % 3 != 0 { "frameshift" } else { "inframe_indel" });
                } else if let Some(gm) = genome.as_mut() {
                    let subs: Vec<(u64, u8, u8)> = rf.bytes().zip(al.bytes()).enumerate().filter(|(_, (x, y))| x != y).map(|(k, (x, y))| (p0 + k as u64, x, y)).collect();
                    match gm.nome(c) {
                        Some(gname) => {
                            let mut mism = false;
                            for &t in &txs {
                                let (eff, mm) = codon_effect(&a.txs[t], &subs, gm, &gname)?;
                                mism |= mm;
                                if let Some(e) = eff {
                                    if consequence.map_or(true, |w| severity(e) > severity(w)) {
                                        consequence = Some(e);
                                        gene = a.txs[t].gene.clone();
                                    }
                                }
                            }
                            if mism {
                                res.ref_mismatch += 1;
                            }
                            if consequence.is_none() {
                                consequence = Some("coding_unknown");
                            }
                        }
                        None => consequence = Some("coding_unknown"),
                    }
                } else {
                    consequence = Some("coding");
                }
            }
            let reason: Option<&str> = if symbolic {
                res.symbolic += 1;
                Some("symbolic_allele")
            } else if opt.pass_only && !pass {
                res.excluded_filter += 1;
                Some("not_pass")
            } else if depth.map_or(false, |d| d < u64::from(opt.min_depth)) {
                res.excluded_depth += 1;
                Some("low_depth")
            } else if vaf.map_or(false, |v| v < opt.min_vaf) {
                res.excluded_vaf += 1;
                Some("low_vaf")
            } else if !in_region {
                res.outside_region += 1;
                Some("outside_region")
            } else if use_coding && coding != Some(true) {
                res.not_coding += 1;
                Some("not_coding")
            } else if opt.nonsynonymous_only && consequence == Some("synonymous") {
                res.synonymous += 1;
                Some("synonymous")
            } else {
                None
            };
            if reason.is_none() {
                res.counted += 1;
                match kind {
                    "snv" => res.counted_snv += 1,
                    "mnv" => res.counted_mnv += 1,
                    _ => res.counted_indel += 1,
                }
                if depth.is_none() && opt.min_depth > 0 {
                    res.without_depth += 1;
                }
                if vaf.is_none() && opt.min_vaf > 0.0 {
                    res.without_vaf += 1;
                }
                if consequence == Some("coding_unknown") {
                    res.coding_unknown += 1;
                }
            }
            if let (Some(w), Some(o)) = (out.as_mut(), opt.output.as_deref()) {
                w.linha(&format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    chrom_raw,
                    p,
                    rf,
                    al,
                    kind,
                    filter,
                    depth.map_or("NA".to_string(), |d| d.to_string()),
                    vaf.map_or("NA".to_string(), |v| format!("{v:.4}")),
                    if in_region { "yes" } else { "no" },
                    coding.map_or("NA", |c| if c { "yes" } else { "no" }),
                    consequence.unwrap_or("NA"),
                    if gene.is_empty() { "NA" } else { &gene },
                    if reason.is_none() { "yes" } else { "no" },
                    reason.unwrap_or("")
                ))
                .map_err(|e| write_error(o, e))?;
            }
        }
    }
    if !seen_header {
        return Err(format!("{vcf}: no #CHROM header line. Is this a VCF file?"));
    }
    if res.records > 0 && !vcf_chroms_matched {
        return Err(format!(
            "No sequence of {vcf} ({}) is in the BED {bed} ({}). Check that both use the same genome and sequence names",
            vcf_chroms.join(", "),
            region.keys().take(5).cloned().collect::<Vec<_>>().join(", ")
        ));
    }
    if let (Some(w), Some(o)) = (out.take(), opt.output.as_deref()) {
        w.termina().map_err(|e| write_error(o, e))?;
    }
    cleanup.disarm();
    let mb = res.megabases;
    res.tmb = res.counted as f64 / mb;
    res.tmb_snv = (res.counted_snv + res.counted_mnv) as f64 / mb;
    res.tmb_indel = res.counted_indel as f64 / mb;
    if opt.pass_only && all_dot_filter && res.records > 0 {
        res.warnings.push("Every FILTER is '.', so nothing is PASS and nothing was counted: the VCF was not filtered. Filter it first (for instance FilterMutectCalls), or set pass_only = false".to_string());
    }
    if res.without_depth > 0 || res.without_vaf > 0 {
        res.warnings.push(format!(
            "{} counted variants had no depth and {} no allele fraction in the VCF (no FORMAT DP/AD/AF, Strelka tiers or INFO DP/AF): the depth and VAF filters could not be applied to them",
            res.without_depth, res.without_vaf
        ));
    }
    if res.ref_mismatch > 0 {
        res.warnings.push(format!("{} substitutions have a REF that differs from the reference genome: is it the same genome build as the VCF?", res.ref_mismatch));
    }
    if res.megabases < 0.5 {
        res.warnings.push(format!("The region is only {:.3} Mb: a TMB over so few bases is imprecise (panels for TMB cover about 1 Mb or more)", res.megabases));
    }
    Ok(res)
}

impl TmbResult {
    pub fn json(&self) -> String {
        format!(
            "{{\"vcf\":{},\"bed\":{},\"sample\":{},\"bed_bp\":{},\"coding_bp\":{},\"region_bp\":{},\"megabases\":{},\"records\":{},\"alleles\":{},\"counted\":{},\"counted_snv\":{},\"counted_mnv\":{},\"counted_indel\":{},\"tmb\":{},\"tmb_snv\":{},\"tmb_indel\":{},\"excluded_filter\":{},\"excluded_depth\":{},\"excluded_vaf\":{},\"outside_region\":{},\"not_coding\":{},\"synonymous\":{},\"symbolic\":{},\"without_depth\":{},\"without_vaf\":{},\"ref_mismatch\":{},\"coding_unknown\":{},\"output\":{},\"warnings\":{}}}",
            js(&self.vcf),
            js(&self.bed),
            self.sample.as_deref().map_or("null".to_string(), js),
            self.bed_bp,
            self.coding_bp.map_or("null".to_string(), |x| x.to_string()),
            self.region_bp,
            jf(self.megabases),
            self.records,
            self.alleles,
            self.counted,
            self.counted_snv,
            self.counted_mnv,
            self.counted_indel,
            jf(self.tmb),
            jf(self.tmb_snv),
            jf(self.tmb_indel),
            self.excluded_filter,
            self.excluded_depth,
            self.excluded_vaf,
            self.outside_region,
            self.not_coding,
            self.synonymous,
            self.symbolic,
            self.without_depth,
            self.without_vaf,
            self.ref_mismatch,
            self.coding_unknown,
            self.output.as_deref().map_or("null".to_string(), js),
            jlist(&self.warnings)
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!(
            "TMB = {:.2} mutations/Mb ({} counted over {:.3} Mb{}): substitutions {:.2}/Mb ({} SNV + {} MNV), indels {:.2}/Mb ({}).\n",
            self.tmb,
            self.counted,
            self.megabases,
            if self.coding_bp.is_some() && self.region_bp == self.coding_bp.unwrap_or(0) { " of CDS in the region" } else { "" },
            self.tmb_snv,
            self.counted_snv,
            self.counted_mnv,
            self.tmb_indel,
            self.counted_indel
        );
        o.push_str(&format!(
            "{} records, {} ALT alleles{}. Not counted: {} not PASS, {} low depth, {} low VAF, {} outside the region, {} not coding, {} synonymous, {} symbolic.\n",
            self.records,
            self.alleles,
            self.sample.as_ref().map_or(String::new(), |s| format!(" (tumor sample {s})")),
            self.excluded_filter,
            self.excluded_depth,
            self.excluded_vaf,
            self.outside_region,
            self.not_coding,
            self.synonymous,
            self.symbolic
        ));
        for w in &self.warnings {
            o.push_str(&format!("warning: {w}\n"));
        }
        o.push_str("This is an experimental estimate for teaching and research, not a clinical TMB.\n");
        o
    }
}
