// =====================================================================
//  OUTBREAKS AND SURVEILLANCE: HOW CLOSE GENOMES ARE, WHAT A VIRUS
//  CHANGED, AND WHAT A PCR PANEL WOULD DETECT
//
//  Three tools for clinical and public-health microbiology:
//
//   1. `genome_distances` / `sketch_distances`  the Mash distance and the
//        ANI estimate between every pair of genomes, from FracMinHash
//        sketches (sourmash's hashes, see microbio.rs), a distance matrix
//        (TSV or PHYLIP) and a tree: neighbor joining (Saitou & Nei 1987,
//        with the formulas of Studier & Keppler 1988), midpoint-rooted, or
//        UPGMA. The tree comes as Newick and as a layout a screen can draw
//        as a rectangular tree. `build_tree` and `tree_from_matrix` build
//        the tree from any distance matrix (a SNP matrix, sourmash's).
//   2. `viral_mutations`  each genome of a multi-FASTA aligned to the
//        reference genome, its substitutions, deletions, insertions and
//        missing (N) ranges in reference coordinates, written as Nextclade
//        writes them (C241T, deletion 21765-21770, insertion
//        22204:GAGCCAGAA), and their amino-acid changes per gene from the
//        reference GFF3 (S:D614G, S:H69-, S:214:EPE, ORF1b:P314L).
//   3. `in_silico_pcr`  a panel of primer pairs (and probes) against
//        genomes (the amplicons, with the mismatches of each primer) and
//        against reads (how many reads carry a primer site), and the
//        presence/absence matrix.
//
//  DISTANCES. Jaccard J of two FracMinHash sketches taken at the same
//  scaled = shared hashes / hashes in either. Mash distance
//  D = -ln(2J / (1 + J)) / k (Ondov et al. 2016), 1 when nothing is shared;
//  ANI = (2J / (1 + J))^(1/k), sourmash's point estimate. With scaled = 0
//  the scaled is chosen per run: each genome asks for the largest of 1000,
//  100, 10 and 1 that keeps at least 10,000 of its hashes (a 30 kb virus:
//  1, every k-mer; a 5 Mb bacterium: 100), and the run takes the largest
//  asked for (a short genome among large ones keeps fewer hashes; under
//  100, a warning says so). Sketches are nested, so each genome is hashed once at
//  the scaled its size calls for and then cut down to the common one.
//
//  NEIGHBOR JOINING. At each step, with r nodes left and R_i = sum_j d(i,j),
//  the pair with the smallest Q(i,j) = (r - 2) d(i,j) - R_i - R_j is
//  joined (the first such pair in input order on a tie); the branches are
//  d(i,u) = d(i,j)/2 + (R_i - R_j) / (2(r - 2)) and d(j,u) = d(i,j) -
//  d(i,u); d(u,k) = (d(i,k) + d(j,k) - d(i,j)) / 2. A negative branch is
//  set to 0 and its length moved to the sibling (Kuhner & Felsenstein
//  1994), so that d(i,j) is kept. The unrooted tree is rooted at the
//  middle of its longest leaf-to-leaf path, or, given an outgroup, in the
//  middle of the outgroup's branch (midpoint rooting can put the root on a
//  long branch inside the ingroup: the 20 SARS-CoV-2 lineages are rooted
//  inside Omicron at the midpoint, and Omicron is one clade when rooted on
//  lineage A). UPGMA joins the closest pair,
//  at height d/2, with size-weighted averages. Children are ordered by
//  size (smaller clade first, then by input order): the layout is
//  deterministic. Layout: x = distance from the root, y = leaf order
//  (0, 1, ...); an internal node sits halfway between its first and last
//  child.
//
//  ALIGNMENT OF A VIRAL GENOME TO ITS REFERENCE. Seeds are 17-mers unique
//  in the reference; consecutive seeds on one diagonal make anchors; the
//  best co-linear chain of anchors (length minus diagonal jumps) is kept,
//  each anchor trimmed by 12 bases so that the dynamic programming sees
//  the context; between anchors, and at both ends, a banded global
//  alignment with affine gaps fills in (Gotoh), with Nextclade's scores:
//  match +3, mismatch -1, gap open 6 outside genes, 7 at a codon boundary
//  and 8 inside a codon, gap extension 0. N and the IUPAC codes match any
//  base they stand for. The ends are free: a genome that starts at
//  position 55 is "not covered" from 1 to 54, not deleted; bases past the
//  ends of the reference are ignored. Then every indel whose place is
//  ambiguous (a deletion in a repeat) is put at the equivalent place where
//  opening the gap costs least -- outside genes, then at a codon boundary,
//  then inside a codon -- the leftmost of those: what Nextclade does
//  (checked on 20 SARS-CoV-2 lineages: the A deleted from AAAA before the
//  start codon of N is 28271, not 28274 inside N).A genome that shares more seeds with the reverse complement
//  of the reference is turned around first. A multi-record reference
//  (segmented viruses): each genome record goes to the reference record
//  it shares most seeds with.
//
//  AMINO ACIDS. Each CDS of the GFF3 (its parts joined in coding order, the
//  phase of the first part honoured, minus strand reverse-complemented) is
//  read through the alignment codon by codon. A codon with its 3 bases and
//  no insertion inside is translated (two changes in one codon give one
//  amino-acid change: N:R203K and N:G204R from GGG>AAC). Where indels
//  touch codons, the codons are taken together until the query is back in
//  frame; their translations are aligned (identity scores, a residue
//  placed as far right as an equally good alignment allows: S:G142-,
//  V143-, Y144- and Y145D for the deletion 21987-21995, as Nextclade 3
//  writes it). An indel whose length is not a multiple of 3 is a
//  frameshift: the codons from it to where the frame comes back (or to the
//  end of the gene) are reported as a frameshift, not as changes. A codon
//  with N, or outside the sequenced part, is unknown (X), never a change.
//  An insertion between two codons whose length is a multiple of 3 is an
//  amino-acid insertion after the first codon (S:214:EPE).
//  ORF1ab (a two-part CDS with a -1 ribosomal frameshift, whose name ends
//  in "ab") is reported as two genes, ORF1a up to the slippage and ORF1b
//  from it, each numbered from its own first codon -- Nextclade's names
//  (ORF1b:P314L is ORF1ab:P4715L); the separate pp1a CDS of the same gene
//  is then not reported again. `split_orf1ab = false` keeps ORF1ab whole.
//  Other multi-part CDS (gag-pol, PA-X, spliced genes) keep their name.
//  When two CDS share a gene name (gene=ORF1ab on both the pp1ab and the
//  pp1a CDS of RefSeq), the first word of the product tells them apart.
//
//  IN-SILICO PCR. Primers and probes may hold IUPAC codes; a template base
//  matches when the code includes it (an N in the template never matches).
//  A primer site has at most `max_mismatches` mismatches and none in its
//  last `three_prime_exact` bases (the 3' end, where the polymerase
//  starts); a product needs a forward site on one strand and a reverse
//  site on the other, facing each other, at most `max_amplicon` bases
//  apart (or the assay's maximum size, if larger), with at most
//  `max_total_mismatches` in the two primers together. Every such product
//  is reported (overlapping products of repeated sites too); it counts as
//  a detection when its size is in the assay's range and, if the assay has
//  a probe, the probe binds inside it (either strand, at most
//  `max_probe_mismatches`, no 3' rule). Genomes are linear: a product
//  across the origin of a circular chromosome is not seen. Reads (FASTQ):
//  a read counts for a target when it holds a site of its forward
//  primer, reverse primer or probe (either strand); a sample is positive
//  with at least `min_reads` such reads.
// =====================================================================

use crate::annotation::translate_dna_to_protein;
use crate::externo::{abre_leitura, Escritor};
use crate::microbio::{load_databases, max_hash_for_scaled, murmur64, SeqKind, SeqReader, SeqRecord, SEED};
use crate::multimerge_core::{multi_merge_sort_stable, multi_merge_sort_stable_by};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, ErrorKind, Read, Write};
use std::sync::OnceLock;

// ---------------------------------------------------------------------
//  SMALL HELPERS: messages, JSON, numbers, names
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

/// An output must not be one of the inputs: it would be erased before being
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

/// Writes a whole text file through the project's writer (`.gz` gives BGZF).
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

fn jlist(items: &[String]) -> String {
    let v: Vec<String> = items.iter().map(|s| js(s)).collect();
    format!("[{}]", v.join(","))
}

fn jopt(s: &Option<String>) -> String {
    s.as_deref().map_or_else(|| "null".to_string(), js)
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

fn clean_cell(s: &str) -> String {
    s.replace(['\t', '\n', '\r'], " ")
}

/// A number for a table: up to 8 decimals, trailing zeros removed.
fn fmt_num(x: f64) -> String {
    if !x.is_finite() {
        return "NA".to_string();
    }
    let s = format!("{x:.8}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s.is_empty() || s == "-0" {
        "0".to_string()
    } else {
        s.to_string()
    }
}

fn sample_text(piece: &str) -> String {
    let mut o: String = piece.chars().take(40).collect();
    if piece.chars().count() > 40 {
        o.push_str("...");
    }
    o
}

/// The file name without folders, compression and sequence extensions:
/// `/data/GCF_1.fna.gz` -> `GCF_1`.
fn file_label(path: &str) -> String {
    let base = std::path::Path::new(path).file_name().map_or_else(|| path.to_string(), |s| s.to_string_lossy().to_string());
    let mut b = base.as_str();
    for ext in [".gz", ".bgz", ".lz4"] {
        if b.to_ascii_lowercase().ends_with(ext) {
            b = &b[..b.len() - ext.len()];
            break;
        }
    }
    for ext in [".fasta", ".fna", ".fa", ".fas", ".ffn", ".fastq", ".fq", ".seq"] {
        if b.to_ascii_lowercase().ends_with(ext) {
            b = &b[..b.len() - ext.len()];
            break;
        }
    }
    if b.is_empty() {
        base
    } else {
        b.to_string()
    }
}

/// Names must be unique (in a matrix, a tree, a table): a repeated name gets
/// `_2`, `_3`...
fn make_unique(names: &mut [String], warnings: &mut Vec<String>) {
    let mut seen: HashMap<String, usize> = HashMap::new();
    for n in names.iter() {
        *seen.entry(n.clone()).or_insert(0) += 1;
    }
    let mut used: HashMap<String, usize> = HashMap::new();
    let taken: std::collections::HashSet<String> = names.iter().cloned().collect();
    for n in names.iter_mut() {
        if seen.get(n.as_str()).copied().unwrap_or(0) < 2 {
            continue;
        }
        let c = used.entry(n.clone()).or_insert(0);
        *c += 1;
        if *c == 1 {
            continue;
        }
        let mut k = *c;
        let mut cand = format!("{n}_{k}");
        while taken.contains(&cand) {
            k += 1;
            cand = format!("{n}_{k}");
        }
        warnings.push(format!("the name {n} is used more than once: its copy number {} is called {cand}", *c));
        *n = cand;
    }
}

/// One record of a FASTA or FASTQ file, whole.
struct Record {
    id: String,
    description: String,
    seq: Vec<u8>,
}

fn read_records(path: &str) -> Result<(Vec<Record>, Option<SeqKind>), String> {
    let mut rd = SeqReader::open(path)?;
    let mut rec = SeqRecord::default();
    let mut out = Vec::new();
    while rd.next(&mut rec)? {
        let head = String::from_utf8_lossy(&rec.head).trim().to_string();
        let id = head.split_whitespace().next().unwrap_or("").to_string();
        out.push(Record { id, description: head, seq: std::mem::take(&mut rec.seq) });
    }
    Ok((out, rd.kind()))
}

// ---------------------------------------------------------------------
//  BASES
// ---------------------------------------------------------------------

/// IUPAC bit mask: A=1, C=2, G=4, T=8 (U too); N and the ambiguity codes
/// cover several bases; anything else is 0.
const NUC_MASK: [u8; 256] = {
    let mut t = [0u8; 256];
    let codes: [(u8, u8); 16] = [
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
        (b'N', 15),
    ];
    let mut i = 0;
    while i < codes.len() {
        t[codes[i].0 as usize] = codes[i].1;
        t[codes[i].0.to_ascii_lowercase() as usize] = codes[i].1;
        i += 1;
    }
    t
};

/// Only A, C, G, T (U) count as bases of a template: an N or an ambiguity
/// code in a genome never matches a primer.
const TEMPLATE_MASK: [u8; 256] = {
    let mut t = [0u8; 256];
    t[b'A' as usize] = 1;
    t[b'a' as usize] = 1;
    t[b'C' as usize] = 2;
    t[b'c' as usize] = 2;
    t[b'G' as usize] = 4;
    t[b'g' as usize] = 4;
    t[b'T' as usize] = 8;
    t[b't' as usize] = 8;
    t[b'U' as usize] = 8;
    t[b'u' as usize] = 8;
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

/// A sequence letter in upper case, U as T; anything that is not an IUPAC
/// nucleotide code becomes N.
fn norm_base(b: u8) -> u8 {
    let u = b.to_ascii_uppercase();
    if u == b'U' {
        b'T'
    } else if NUC_MASK[usize::from(u)] != 0 {
        u
    } else {
        b'N'
    }
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        o => o,
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| complement(b)).collect()
}

fn is_acgt(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// The standard genetic code, built once from annotation.rs.
fn codon_table() -> &'static [u8; 64] {
    static T: OnceLock<[u8; 64]> = OnceLock::new();
    T.get_or_init(|| {
        let mut t = [b'X'; 64];
        let bases = [b'A', b'C', b'G', b'T'];
        for (i, slot) in t.iter_mut().enumerate() {
            let c = [bases[i >> 4], bases[(i >> 2) & 3], bases[i & 3]];
            let p = translate_dna_to_protein(std::str::from_utf8(&c).unwrap_or("NNN"));
            *slot = p.bytes().next().unwrap_or(b'X');
        }
        t
    })
}

/// One codon: an ambiguity code (R, Y...) gives the amino acid when every
/// base it stands for gives the same one; an N, a gap or a mixed result
/// gives X.
fn translate_codon(c: &[u8]) -> u8 {
    if c.len() != 3 {
        return b'X';
    }
    let table = codon_table();
    let (a, b, d) = (CODE[usize::from(c[0])], CODE[usize::from(c[1])], CODE[usize::from(c[2])]);
    if a < 4 && b < 4 && d < 4 {
        return table[usize::from(a * 16 + b * 4 + d)];
    }
    if c.iter().any(|&x| matches!(x, b'N' | b'n' | b'-')) {
        return b'X';
    }
    let masks = [NUC_MASK[usize::from(c[0])], NUC_MASK[usize::from(c[1])], NUC_MASK[usize::from(c[2])]];
    if masks.iter().any(|&m| m == 0) {
        return b'X';
    }
    let mut aa: Option<u8> = None;
    for x in 0..4u8 {
        if masks[0] & (1 << x) == 0 {
            continue;
        }
        for y in 0..4u8 {
            if masks[1] & (1 << y) == 0 {
                continue;
            }
            for z in 0..4u8 {
                if masks[2] & (1 << z) == 0 {
                    continue;
                }
                let v = table[usize::from(x * 16 + y * 4 + z)];
                match aa {
                    None => aa = Some(v),
                    Some(p) if p != v => return b'X',
                    _ => {}
                }
            }
        }
    }
    aa.unwrap_or(b'X')
}

fn translate(seq: &[u8]) -> Vec<u8> {
    seq.chunks_exact(3).map(translate_codon).collect()
}

// =====================================================================
//  1. DISTANCES BETWEEN GENOMES AND TREES
// =====================================================================

/// One genome: one or more FASTA files whose sequences together are the
/// genome (chromosome and plasmids, or the contigs of an assembly).
#[derive(Clone, Debug, Default)]
pub struct GenomeFiles {
    pub files: Vec<String>,
    /// The name in the matrix and the tree. Empty: the file name, without
    /// folders and extensions.
    pub name: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TreeMethod {
    NeighborJoining,
    Upgma,
}

impl TreeMethod {
    pub fn from_text(s: &str) -> Result<TreeMethod, String> {
        match s.trim().to_ascii_lowercase().replace(['-', '_', ' '], "").as_str() {
            "nj" | "neighborjoining" | "neighbourjoining" => Ok(TreeMethod::NeighborJoining),
            "upgma" => Ok(TreeMethod::Upgma),
            _ => Err(format!(
                "Unknown tree method \"{s}\": use \"nj\" (neighbor joining, midpoint-rooted: the default, it does not assume a molecular clock) or \"upgma\" (assumes every lineage evolves at the same rate)"
            )),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            TreeMethod::NeighborJoining => "nj",
            TreeMethod::Upgma => "upgma",
        }
    }
}

/// What the matrix file holds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatrixValue {
    /// Mash distance (0 = identical; what the tree is built from).
    Distance,
    /// ANI estimate (1 = identical).
    Ani,
    /// Jaccard of the sketches (1 = identical), what `sourmash compare` gives.
    Jaccard,
}

impl MatrixValue {
    pub fn from_text(s: &str) -> Result<MatrixValue, String> {
        match s.trim().to_ascii_lowercase().as_str() {
            "distance" | "mash" | "mash_distance" => Ok(MatrixValue::Distance),
            "ani" => Ok(MatrixValue::Ani),
            "jaccard" | "similarity" => Ok(MatrixValue::Jaccard),
            _ => Err(format!("Unknown matrix value \"{s}\": use \"distance\" (Mash distance, the default), \"ani\" or \"jaccard\"")),
        }
    }

    fn name(self) -> &'static str {
        match self {
            MatrixValue::Distance => "distance",
            MatrixValue::Ani => "ani",
            MatrixValue::Jaccard => "jaccard",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatrixFormat {
    /// A header line of names, then one line per genome: name and values.
    Tsv,
    /// The number of genomes, then one line per genome: name and values
    /// (relaxed PHYLIP: names of any length, spaces turned into `_`).
    Phylip,
}

impl MatrixFormat {
    pub fn from_text(s: &str) -> Result<MatrixFormat, String> {
        match s.trim().to_ascii_lowercase().as_str() {
            "tsv" | "tab" => Ok(MatrixFormat::Tsv),
            "phylip" | "phy" => Ok(MatrixFormat::Phylip),
            _ => Err(format!("Unknown matrix format \"{s}\": use \"tsv\" or \"phylip\"")),
        }
    }

    /// By the file name: `.phy`/`.phylip` (also before `.gz`) is PHYLIP.
    fn for_path(path: &str) -> MatrixFormat {
        let p = path.to_ascii_lowercase();
        let p = p.strip_suffix(".gz").unwrap_or(&p);
        if p.ends_with(".phy") || p.ends_with(".phylip") || p.ends_with(".dist") {
            MatrixFormat::Phylip
        } else {
            MatrixFormat::Tsv
        }
    }
}

#[derive(Clone, Debug)]
pub struct DistanceOptions {
    /// k-mer length, 1 to 32 (31, sourmash's default; 21 is Mash's). Ignored
    /// for sketch databases, which carry their own k.
    pub k: u8,
    /// One k-mer in `scaled` is kept; 0 = chosen from the genome sizes (see
    /// the top of this file). A sketch database cannot be made finer than it
    /// was built.
    pub scaled: u64,
    /// Each FASTA record is a genome of its own (a multi-FASTA of viral
    /// genomes); otherwise all the records of a genome's files are one
    /// genome.
    pub records_as_genomes: bool,
    pub method: TreeMethod,
    /// Neighbor joining: root on the branch of this genome (an outgroup:
    /// the reference, an older lineage) instead of at the midpoint.
    pub outgroup: Option<String>,
    /// The distance matrix is written here, if given.
    pub matrix_output: Option<String>,
    /// TSV or PHYLIP; `None`: by the name (`.phy`, `.phylip` = PHYLIP).
    pub matrix_format: Option<MatrixFormat>,
    pub matrix_value: MatrixValue,
    /// The tree is written here as Newick, if given.
    pub tree_output: Option<String>,
}

impl Default for DistanceOptions {
    fn default() -> Self {
        DistanceOptions {
            k: 31,
            scaled: 0,
            records_as_genomes: false,
            method: TreeMethod::NeighborJoining,
            outgroup: None,
            matrix_output: None,
            matrix_format: None,
            matrix_value: MatrixValue::Distance,
            tree_output: None,
        }
    }
}

/// One genome as sketched.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct GenomeSketchInfo {
    pub name: String,
    pub files: Vec<String>,
    pub sequences: u64,
    pub bases: u64,
    /// Distinct hashes kept at the scaled of the run.
    pub hashes: u64,
}

/// The tree a screen can draw: each node with its branch to the parent,
/// its x (distance from the root) and y (leaf order).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TreeNode {
    /// The genome name for a leaf; empty for an internal node.
    pub name: String,
    /// The index of the genome (row of the matrix) for a leaf.
    pub taxon: Option<usize>,
    pub parent: Option<usize>,
    /// In drawing order.
    pub children: Vec<usize>,
    /// Length of the branch to the parent (0 for the root).
    pub branch_length: f64,
    pub x: f64,
    pub y: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Tree {
    pub method: TreeMethod,
    pub nodes: Vec<TreeNode>,
    pub root: usize,
    /// The genomes (matrix rows) from top to bottom.
    pub leaf_order: Vec<usize>,
    pub newick: String,
}

#[derive(Clone, Debug)]
pub struct DistanceResult {
    pub k: u8,
    pub scaled: u64,
    /// "fasta" or "sketch database".
    pub source: String,
    pub genomes: Vec<GenomeSketchInfo>,
    /// Hashes shared by each pair (the diagonal: the genome's own hashes).
    pub shared: Vec<Vec<u64>>,
    pub jaccard: Vec<Vec<f64>>,
    pub distance: Vec<Vec<f64>>,
    pub ani: Vec<Vec<f64>>,
    pub tree: Tree,
    pub matrix_output: Option<String>,
    pub matrix_format: MatrixFormat,
    pub matrix_value: MatrixValue,
    pub tree_output: Option<String>,
    pub warnings: Vec<String>,
}

/// Mash distance from the Jaccard of two FracMinHash sketches:
/// D = -ln(2J/(1+J))/k; 1 when nothing is shared (and never above 1).
pub fn mash_distance(j: f64, k: u8) -> f64 {
    if !(j > 0.0) {
        return 1.0;
    }
    if j >= 1.0 {
        return 0.0;
    }
    (-(2.0 * j / (1.0 + j)).ln() / f64::from(k.max(1))).min(1.0)
}

/// ANI estimate from the Jaccard: (2J/(1+J))^(1/k), sourmash's.
pub fn ani_from_jaccard(j: f64, k: u8) -> f64 {
    if !(j > 0.0) {
        return 0.0;
    }
    if j >= 1.0 {
        return 1.0;
    }
    (2.0 * j / (1.0 + j)).powf(1.0 / f64::from(k.max(1)))
}

/// Calls `out.push(hash)` for every canonical k-mer of `seq` whose
/// sourmash hash is `<= max_hash`; returns the number of valid k-mers.
fn hash_kmers(seq: &[u8], k: usize, max_hash: u64, out: &mut Vec<u64>) -> u64 {
    let n = seq.len();
    if n < k || k == 0 {
        return 0;
    }
    let up: Vec<u8> = seq.iter().map(|&b| if CODE[usize::from(b)] < 4 { b"ACGT"[usize::from(CODE[usize::from(b)])] } else { b'N' }).collect();
    let rc: Vec<u8> = up.iter().rev().map(|&b| complement(b)).collect();
    let mask = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let top = 2 * (k - 1);
    let (mut fw, mut rv, mut valid) = (0u64, 0u64, 0usize);
    let mut count = 0u64;
    for i in 0..n {
        let c = CODE[usize::from(up[i])];
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
            let h = if fw <= rv {
                murmur64(&up[s..s + k], SEED)
            } else {
                let r = n - s - k;
                murmur64(&rc[r..r + k], SEED)
            };
            if h <= max_hash {
                out.push(h);
            }
        }
    }
    count
}

const SKETCH_WINDOW: usize = 1 << 20;

/// The distinct kept hashes of a genome (its sequences together), sorted,
/// and its number of valid k-mers.
fn sketch_genome(seqs: &[&[u8]], k: usize, max_hash: u64) -> (Vec<u64>, u64) {
    let mut units: Vec<(usize, usize, usize)> = Vec::new();
    for (r, s) in seqs.iter().enumerate() {
        let n = s.len();
        let mut a = 0;
        while a < n {
            units.push((r, a, (a + SKETCH_WINDOW + k - 1).min(n)));
            a += SKETCH_WINDOW;
        }
    }
    let (mut all, kmers) = units
        .par_iter()
        .map(|&(r, a, b)| {
            let mut v = Vec::new();
            let c = hash_kmers(&seqs[r][a..b], k, max_hash, &mut v);
            (v, c)
        })
        .reduce(
            || (Vec::new(), 0u64),
            |(mut x, cx), (y, cy)| {
                x.extend_from_slice(&y);
                (x, cx + cy)
            },
        );
    multi_merge_sort_stable(&mut all);
    all.dedup();
    (all, kmers)
}

fn shared_count(a: &[u64], b: &[u64]) -> u64 {
    let (mut i, mut j, mut n) = (0, 0, 0u64);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                n += 1;
                i += 1;
                j += 1;
            }
        }
    }
    n
}

/// The scaled a genome of `bases` bases asks for: the largest of 1000, 100,
/// 10 and 1 that keeps at least 10,000 of its hashes.
fn auto_scaled(bases: u64) -> u64 {
    for s in [1000u64, 100, 10] {
        if bases / s >= 10_000 {
            return s;
        }
    }
    1
}

struct Sketched {
    info: GenomeSketchInfo,
    hashes: Vec<u64>,
    /// The scaled the hashes were taken at.
    scaled: u64,
}

/// Distances between genomes given as FASTA files, and the tree.
pub fn genome_distances(genomes: &[GenomeFiles], opt: &DistanceOptions) -> Result<DistanceResult, String> {
    if opt.k == 0 || opt.k > 32 {
        return Err(format!("Invalid k-mer length: k = {}. Use a k from 1 to 32 (31 is sourmash's default, 21 Mash's)", opt.k));
    }
    if genomes.is_empty() {
        return Err("No genome given: give at least one FASTA file (two or more to have distances)".to_string());
    }
    let mut inputs: Vec<&str> = Vec::new();
    for (i, g) in genomes.iter().enumerate() {
        if g.files.is_empty() {
            return Err(format!(
                "Genome {} ({}) has no file: give at least one FASTA file for each genome",
                i + 1,
                if g.name.is_empty() { "unnamed" } else { &g.name }
            ));
        }
        for f in &g.files {
            check_input(f)?;
            inputs.push(f);
        }
    }
    for o in [&opt.matrix_output, &opt.tree_output].into_iter().flatten() {
        not_an_input(o, &inputs)?;
    }
    let k = usize::from(opt.k);
    let mut sk: Vec<Sketched> = Vec::new();
    let mut warnings = Vec::new();
    let sketch_one = |name: String, files: Vec<String>, seqs: &[&[u8]], sk: &mut Vec<Sketched>| {
        let bases: u64 = seqs.iter().map(|s| s.len() as u64).sum();
        let scaled = if opt.scaled == 0 { auto_scaled(bases) } else { opt.scaled };
        let (hashes, _) = sketch_genome(seqs, k, max_hash_for_scaled(scaled));
        sk.push(Sketched {
            info: GenomeSketchInfo { name, files, sequences: seqs.len() as u64, bases, hashes: 0 },
            hashes,
            scaled,
        });
    };
    for g in genomes {
        if opt.records_as_genomes {
            for f in &g.files {
                let (recs, kind) = read_records(f)?;
                if kind == Some(SeqKind::Fastq) {
                    return Err(format!("{f} is a FASTQ file (reads): distances are computed between genomes (FASTA). Assemble the reads first, or give the assembly"));
                }
                if recs.is_empty() {
                    return Err(format!("{f} holds no sequence: give a FASTA file with at least one record"));
                }
                for r in &recs {
                    let name = if r.id.is_empty() { file_label(f) } else { r.id.clone() };
                    sketch_one(name, vec![f.clone()], &[&r.seq], &mut sk);
                }
            }
        } else {
            let mut recs: Vec<Record> = Vec::new();
            for f in &g.files {
                let (r, kind) = read_records(f)?;
                if kind == Some(SeqKind::Fastq) {
                    return Err(format!("{f} is a FASTQ file (reads): distances are computed between genomes (FASTA). Assemble the reads first, or give the assembly"));
                }
                recs.extend(r);
            }
            if recs.is_empty() {
                return Err(format!("{} holds no sequence: give a FASTA file with at least one record", g.files.join(", ")));
            }
            let name = if g.name.trim().is_empty() { file_label(&g.files[0]) } else { g.name.trim().to_string() };
            let seqs: Vec<&[u8]> = recs.iter().map(|r| r.seq.as_slice()).collect();
            sketch_one(name, g.files.clone(), &seqs, &mut sk);
        }
    }
    let scaled = if opt.scaled == 0 { sk.iter().map(|s| s.scaled).max().unwrap_or(1) } else { opt.scaled };
    let mh = max_hash_for_scaled(scaled);
    let mut names: Vec<String> = Vec::new();
    let mut infos = Vec::new();
    let mut hashes = Vec::new();
    for s in sk {
        let mut h = s.hashes;
        let keep = h.partition_point(|&x| x <= mh);
        h.truncate(keep);
        names.push(s.info.name.clone());
        infos.push(s.info);
        hashes.push(h);
    }
    make_unique(&mut names, &mut warnings);
    for (i, n) in infos.iter_mut().zip(&names) {
        i.name = n.clone();
    }
    finish_distances(opt.k, scaled, "fasta", infos, hashes, warnings, opt)
}

/// Distances between the sketches of one or more sketch databases (built by
/// `microbio::sketch_references`); `opt.scaled` (if not 0) makes them
/// coarser.
pub fn sketch_distances(databases: &[String], opt: &DistanceOptions) -> Result<DistanceResult, String> {
    let ins: Vec<&str> = databases.iter().map(String::as_str).collect();
    for o in [&opt.matrix_output, &opt.tree_output].into_iter().flatten() {
        not_an_input(o, &ins)?;
    }
    let db = load_databases(databases, if opt.scaled == 0 { None } else { Some(opt.scaled) })?;
    if db.sketches.is_empty() {
        return Err(format!("The sketch database {} holds no genome: sketch the references again", databases.join(", ")));
    }
    let mut warnings = Vec::new();
    let mut names: Vec<String> = db.sketches.iter().map(|s| s.name.clone()).collect();
    make_unique(&mut names, &mut warnings);
    let mut infos = Vec::new();
    let mut hashes = Vec::new();
    for (s, n) in db.sketches.into_iter().zip(names) {
        infos.push(GenomeSketchInfo { name: n, files: s.files, sequences: s.sequences, bases: s.bases, hashes: 0 });
        hashes.push(s.hashes);
    }
    finish_distances(db.k, db.scaled, "sketch database", infos, hashes, warnings, opt)
}

fn finish_distances(
    k: u8,
    scaled: u64,
    source: &str,
    mut infos: Vec<GenomeSketchInfo>,
    hashes: Vec<Vec<u64>>,
    mut warnings: Vec<String>,
    opt: &DistanceOptions,
) -> Result<DistanceResult, String> {
    let n = infos.len();
    for (i, h) in infos.iter_mut().zip(&hashes) {
        i.hashes = h.len() as u64;
        if h.is_empty() {
            warnings.push(format!(
                "{} keeps no hash at scaled {scaled} (k = {k}): its distance to every genome is 1. Give scaled = 1 (or a smaller scaled) for genomes this short",
                i.name
            ));
        } else if h.len() < 100 {
            warnings.push(format!(
                "{} keeps only {} hashes at scaled {scaled}: its distances are rough. A smaller scaled (1 keeps every k-mer) gives finer distances",
                i.name,
                h.len()
            ));
        }
    }
    let pairs: Vec<(usize, usize)> = (0..n).flat_map(|i| (i + 1..n).map(move |j| (i, j))).collect();
    let counts: Vec<u64> = pairs.par_iter().map(|&(i, j)| shared_count(&hashes[i], &hashes[j])).collect();
    let mut shared = vec![vec![0u64; n]; n];
    let mut jac = vec![vec![0.0f64; n]; n];
    let mut dist = vec![vec![0.0f64; n]; n];
    let mut ani = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        let own = hashes[i].len() as u64;
        shared[i][i] = own;
        jac[i][i] = if own > 0 { 1.0 } else { 0.0 };
        dist[i][i] = 0.0;
        ani[i][i] = if own > 0 { 1.0 } else { 0.0 };
    }
    for (&(i, j), &c) in pairs.iter().zip(&counts) {
        let union = hashes[i].len() as u64 + hashes[j].len() as u64 - c;
        let jj = if union == 0 { 0.0 } else { c as f64 / union as f64 };
        shared[i][j] = c;
        shared[j][i] = c;
        jac[i][j] = jj;
        jac[j][i] = jj;
        let d = mash_distance(jj, k);
        dist[i][j] = d;
        dist[j][i] = d;
        let a = ani_from_jaccard(jj, k);
        ani[i][j] = a;
        ani[j][i] = a;
    }
    let names: Vec<String> = infos.iter().map(|g| g.name.clone()).collect();
    let tree = build_tree(&names, &dist, opt.method, opt.outgroup.as_deref())?;
    let format = opt.matrix_format.unwrap_or_else(|| opt.matrix_output.as_deref().map_or(MatrixFormat::Tsv, MatrixFormat::for_path));
    let mut r = DistanceResult {
        k,
        scaled,
        source: source.to_string(),
        genomes: infos,
        shared,
        jaccard: jac,
        distance: dist,
        ani,
        tree,
        matrix_output: opt.matrix_output.clone(),
        matrix_format: format,
        matrix_value: opt.matrix_value,
        tree_output: opt.tree_output.clone(),
        warnings,
    };
    if format == MatrixFormat::Phylip && opt.matrix_value != MatrixValue::Distance {
        r.warnings.push("a PHYLIP matrix is read by tree programs as distances: this one holds similarities (ANI or Jaccard)".to_string());
    }
    let mut cleanup = Cleanup(Vec::new());
    if let Some(p) = &opt.matrix_output {
        cleanup.0.push(p.clone());
        write_text(p, &r.matrix_text(opt.matrix_value, format))?;
    }
    if let Some(p) = &opt.tree_output {
        cleanup.0.push(p.clone());
        write_text(p, &format!("{}\n", r.tree.newick))?;
    }
    cleanup.disarm();
    Ok(r)
}

/// A PHYLIP name: no white space, no characters that break the format.
fn phylip_name(s: &str) -> String {
    let t: String = s.chars().map(|c| if c.is_whitespace() || matches!(c, '(' | ')' | ':' | ';' | ',' | '[' | ']') { '_' } else { c }).collect();
    if t.is_empty() {
        "_".to_string()
    } else {
        t
    }
}

impl DistanceResult {
    fn values(&self, v: MatrixValue) -> &Vec<Vec<f64>> {
        match v {
            MatrixValue::Distance => &self.distance,
            MatrixValue::Ani => &self.ani,
            MatrixValue::Jaccard => &self.jaccard,
        }
    }

    /// The matrix as written to `matrix_output`.
    pub fn matrix_text(&self, value: MatrixValue, format: MatrixFormat) -> String {
        let m = self.values(value);
        let names: Vec<String> = self.genomes.iter().map(|g| g.name.clone()).collect();
        let mut o = String::new();
        match format {
            MatrixFormat::Tsv => {
                o.push_str("genome");
                for n in &names {
                    o.push('\t');
                    o.push_str(&clean_cell(n));
                }
                o.push('\n');
                for (n, row) in names.iter().zip(m) {
                    o.push_str(&clean_cell(n));
                    for x in row {
                        o.push('\t');
                        o.push_str(&fmt_num(*x));
                    }
                    o.push('\n');
                }
            }
            MatrixFormat::Phylip => {
                o.push_str(&format!("{}\n", names.len()));
                for (n, row) in names.iter().zip(m) {
                    o.push_str(&phylip_name(n));
                    for x in row {
                        o.push(' ');
                        o.push_str(&fmt_num(*x));
                    }
                    o.push('\n');
                }
            }
        }
        o
    }

    pub fn json(&self) -> String {
        let g: Vec<String> = self
            .genomes
            .iter()
            .map(|g| {
                format!(
                    "{{\"name\":{},\"files\":{},\"sequences\":{},\"bases\":{},\"hashes\":{}}}",
                    js(&g.name),
                    jlist(&g.files),
                    g.sequences,
                    g.bases,
                    g.hashes
                )
            })
            .collect();
        let mat = |m: &Vec<Vec<f64>>| {
            let rows: Vec<String> = m.iter().map(|r| format!("[{}]", r.iter().map(|x| jf(*x)).collect::<Vec<_>>().join(","))).collect();
            format!("[{}]", rows.join(","))
        };
        let sh: Vec<String> = self.shared.iter().map(|r| format!("[{}]", r.iter().map(u64::to_string).collect::<Vec<_>>().join(","))).collect();
        format!(
            "{{\"k\":{},\"scaled\":{},\"source\":{},\"names\":{},\"genomes\":[{}],\"shared_hashes\":[{}],\"jaccard\":{},\"distance\":{},\"ani\":{},\"tree\":{},\"matrix_output\":{},\"matrix_format\":{},\"matrix_value\":{},\"tree_output\":{},\"warnings\":{},\"text\":{}}}",
            self.k,
            self.scaled,
            js(&self.source),
            jlist(&self.genomes.iter().map(|g| g.name.clone()).collect::<Vec<_>>()),
            g.join(","),
            sh.join(","),
            mat(&self.jaccard),
            mat(&self.distance),
            mat(&self.ani),
            self.tree.json(),
            jopt(&self.matrix_output),
            js(match self.matrix_format {
                MatrixFormat::Tsv => "tsv",
                MatrixFormat::Phylip => "phylip",
            }),
            js(self.matrix_value.name()),
            jopt(&self.tree_output),
            jlist(&self.warnings),
            js(&self.resumo())
        )
    }

    pub fn resumo(&self) -> String {
        let n = self.genomes.len();
        let mut o = format!(
            "Distances between {n} genome(s) ({}): k = {}, scaled = {} ({}), tree by {}\n",
            self.source,
            self.k,
            self.scaled,
            if self.scaled == 1 { "every k-mer".to_string() } else { format!("one k-mer in {} kept", self.scaled) },
            match self.tree.method {
                TreeMethod::NeighborJoining => "neighbor joining, rooted at the midpoint",
                TreeMethod::Upgma => "UPGMA",
            }
        );
        for g in &self.genomes {
            o.push_str(&format!("  {:>9} hashes  {:>13} bases  {:>5} sequence(s)  {}\n", thousands(g.hashes), thousands(g.bases), g.sequences, g.name));
        }
        if n >= 2 {
            // The closest neighbour of each genome.
            o.push_str("  closest genome of each (Mash distance, ANI):\n");
            for i in 0..n {
                let best = (0..n).filter(|&j| j != i).min_by(|&a, &b| self.distance[i][a].partial_cmp(&self.distance[i][b]).unwrap_or(std::cmp::Ordering::Equal));
                if let Some(j) = best {
                    o.push_str(&format!(
                        "    {}  ->  {}  (D = {}, ANI = {:.4}%)\n",
                        self.genomes[i].name,
                        self.genomes[j].name,
                        fmt_num(self.distance[i][j]),
                        100.0 * self.ani[i][j]
                    ));
                }
            }
        }
        o.push_str("  tree:\n");
        for l in self.tree.resumo().lines() {
            o.push_str("    ");
            o.push_str(l);
            o.push('\n');
        }
        if let Some(p) = &self.matrix_output {
            o.push_str(&format!("  matrix ({}) written to {p}\n", self.matrix_value.name()));
        }
        if let Some(p) = &self.tree_output {
            o.push_str(&format!("  tree (Newick) written to {p}\n"));
        }
        for w in &self.warnings {
            o.push_str(&format!("  warning: {w}\n"));
        }
        o
    }
}

// ---------------------------------------------------------------------
//  Trees
// ---------------------------------------------------------------------

fn check_matrix(names: &[String], d: &[Vec<f64>], whence: &str) -> Result<(), String> {
    let n = names.len();
    if n == 0 {
        return Err(format!("{whence}: no genome, so no tree: give at least one"));
    }
    if d.len() != n || d.iter().any(|r| r.len() != n) {
        return Err(format!("{whence}: the distance matrix must be square, one row and one column per name ({n} names)"));
    }
    for i in 0..n {
        for j in 0..n {
            let x = d[i][j];
            if !x.is_finite() || x < 0.0 {
                return Err(format!(
                    "{whence}: the distance between {} and {} is {x}: distances must be numbers >= 0",
                    names[i], names[j]
                ));
            }
            let y = d[j][i];
            if (x - y).abs() > 1e-9 * x.abs().max(y.abs()).max(1.0) {
                return Err(format!(
                    "{whence}: the matrix is not symmetric: {} -> {} is {x} but {} -> {} is {y}. A distance matrix has the same value on both sides of the diagonal",
                    names[i], names[j], names[j], names[i]
                ));
            }
        }
        if d[i][i] > 1e-12 {
            return Err(format!("{whence}: the distance of {} to itself is {}, not 0: is this a similarity matrix (1 on the diagonal)? Convert it to distances (1 - similarity) first", names[i], d[i][i]));
        }
    }
    Ok(())
}

/// Neighbor joining: the unrooted tree as an adjacency list; nodes 0..n are
/// the taxa, then the internal nodes.
fn nj_unrooted(d: &[Vec<f64>]) -> Vec<Vec<(usize, f64)>> {
    let n = d.len();
    // An n x n matrix of slots: the node made by joining two slots takes the
    // first one's place (and goes to the end of the scanning order).
    let mut m: Vec<Vec<f64>> = d.to_vec();
    let mut node: Vec<usize> = (0..n).collect();
    let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut active: Vec<usize> = (0..n).collect();
    while active.len() > 2 {
        let r = active.len();
        let sums: Vec<f64> = active.iter().map(|&i| active.iter().map(|&j| m[i][j]).sum()).collect();
        let mut best: Option<(f64, usize, usize)> = None;
        for a in 0..r {
            for b in a + 1..r {
                let q = (r as f64 - 2.0) * m[active[a]][active[b]] - sums[a] - sums[b];
                if best.map_or(true, |(bq, _, _)| q < bq) {
                    best = Some((q, a, b));
                }
            }
        }
        let (_, a, b) = best.unwrap_or((0.0, 0, 1));
        let (i, j) = (active[a], active[b]);
        let dij = m[i][j];
        let mut li = dij / 2.0 + (sums[a] - sums[b]) / (2.0 * (r as f64 - 2.0));
        let mut lj = dij - li;
        if li < 0.0 {
            li = 0.0;
            lj = dij;
        } else if lj < 0.0 {
            lj = 0.0;
            li = dij;
        }
        let u = adj.len();
        adj.push(Vec::new());
        let (ni, nj) = (node[i], node[j]);
        adj[u].push((ni, li));
        adj[ni].push((u, li));
        adj[u].push((nj, lj));
        adj[nj].push((u, lj));
        for &k in &active {
            if k != i && k != j {
                let v = (m[i][k] + m[j][k] - dij) / 2.0;
                m[i][k] = v;
                m[k][i] = v;
            }
        }
        node[i] = u;
        active.retain(|&x| x != i && x != j);
        active.push(i);
    }
    if active.len() == 2 {
        let (a, b) = (active[0], active[1]);
        let l = m[a][b].max(0.0);
        let (na, nb) = (node[a], node[b]);
        adj[na].push((nb, l));
        adj[nb].push((na, l));
    }
    adj
}

/// Distances from `from` to every node of an unrooted tree, and the parent
/// of each node on the way.
fn tree_walk(adj: &[Vec<(usize, f64)>], from: usize) -> (Vec<f64>, Vec<Option<usize>>) {
    let mut dist = vec![f64::NAN; adj.len()];
    let mut parent = vec![None; adj.len()];
    dist[from] = 0.0;
    let mut stack = vec![from];
    while let Some(u) = stack.pop() {
        for &(v, l) in &adj[u] {
            if dist[v].is_nan() {
                dist[v] = dist[u] + l;
                parent[v] = Some(u);
                stack.push(v);
            }
        }
    }
    (dist, parent)
}

/// Roots an unrooted tree at the middle of its longest leaf-to-leaf path;
/// returns the root (a new node when the middle falls on a branch).
fn midpoint_root(adj: &mut Vec<Vec<(usize, f64)>>, n: usize) -> usize {
    if n <= 1 {
        return 0;
    }
    let mut best = (0usize, 1usize, -1.0f64);
    for a in 0..n {
        let (dist, _) = tree_walk(adj, a);
        for b in a + 1..n {
            if dist[b] > best.2 * (1.0 + 1e-12) + 1e-300 {
                best = (a, b, dist[b]);
            }
        }
    }
    let (a, b, len) = best;
    let (dist, parent) = tree_walk(adj, a);
    let mut path = vec![b];
    while let Some(p) = parent[*path.last().unwrap_or(&a)] {
        path.push(p);
    }
    path.reverse();
    let half = len.max(0.0) / 2.0;
    let eps = 1e-12 * len.abs().max(1.0);
    let mut i = 0;
    while i + 2 < path.len() && dist[path[i + 1]] < half - eps {
        i += 1;
    }
    let (u, v) = (path[i], path[i + 1]);
    if (half - dist[u]).abs() <= eps && u >= n {
        return u;
    }
    if (dist[v] - half).abs() <= eps && v >= n {
        return v;
    }
    let r = adj.len();
    adj.push(Vec::new());
    adj[u].retain(|&(x, _)| x != v);
    adj[v].retain(|&(x, _)| x != u);
    let lu = (half - dist[u]).max(0.0);
    let lv = (dist[v] - half).max(0.0);
    adj[u].push((r, lu));
    adj[r].push((u, lu));
    adj[v].push((r, lv));
    adj[r].push((v, lv));
    r
}

/// UPGMA: (children, branch length to parent) of every node, and the root.
fn upgma_nodes(d: &[Vec<f64>]) -> (Vec<Vec<usize>>, Vec<f64>, usize) {
    let n = d.len();
    // Slots as in nj_unrooted: n x n.
    let mut m: Vec<Vec<f64>> = d.to_vec();
    let mut node: Vec<usize> = (0..n).collect();
    let mut csize = vec![1.0f64; n];
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut blen = vec![0.0f64; n];
    let mut height = vec![0.0f64; n];
    let mut active: Vec<usize> = (0..n).collect();
    while active.len() > 1 {
        let mut best: Option<(f64, usize, usize)> = None;
        for a in 0..active.len() {
            for b in a + 1..active.len() {
                let x = m[active[a]][active[b]];
                if best.map_or(true, |(bx, _, _)| x < bx) {
                    best = Some((x, a, b));
                }
            }
        }
        let (dij, a, b) = best.unwrap_or((0.0, 0, 1));
        let (i, j) = (active[a], active[b]);
        let (ni, nj) = (node[i], node[j]);
        let u = children.len();
        let h = dij / 2.0;
        children.push(vec![ni, nj]);
        height.push(h);
        blen.push(0.0);
        blen[ni] = (h - height[ni]).max(0.0);
        blen[nj] = (h - height[nj]).max(0.0);
        for &k in &active {
            if k != i && k != j {
                let v = (csize[i] * m[i][k] + csize[j] * m[j][k]) / (csize[i] + csize[j]);
                m[i][k] = v;
                m[k][i] = v;
            }
        }
        csize[i] += csize[j];
        node[i] = u;
        active.retain(|&x| x != i && x != j);
        active.push(i);
    }
    let root = active.first().map_or(0, |&s| node[s]);
    (children, blen, root)
}

/// A Newick name: quoted when it holds characters of the format.
fn newick_name(s: &str) -> String {
    if !s.is_empty() && !s.chars().any(|c| c.is_whitespace() || matches!(c, '(' | ')' | '[' | ']' | '\'' | ':' | ';' | ',')) {
        s.to_string()
    } else {
        format!("'{}'", s.replace('\'', "''"))
    }
}

/// From children lists and branch lengths: the ladderized order, the layout
/// and the Newick text.
fn finish_tree(method: TreeMethod, names: &[String], children: Vec<Vec<usize>>, blen: Vec<f64>, root: usize) -> Tree {
    let total = children.len();
    let n = names.len();
    let mut parent = vec![None; total];
    for (u, ch) in children.iter().enumerate() {
        for &c in ch {
            parent[c] = Some(u);
        }
    }
    // Post-order: leaf count and smallest taxon of each subtree.
    let mut order = Vec::with_capacity(total);
    let mut stack = vec![root];
    while let Some(u) = stack.pop() {
        order.push(u);
        for &c in &children[u] {
            stack.push(c);
        }
    }
    let mut leaves = vec![0usize; total];
    let mut min_taxon = vec![usize::MAX; total];
    for &u in order.iter().rev() {
        if u < n {
            leaves[u] = 1;
            min_taxon[u] = u;
        }
        for &c in &children[u] {
            leaves[u] += leaves[c];
            min_taxon[u] = min_taxon[u].min(min_taxon[c]);
        }
    }
    let mut children = children;
    for ch in children.iter_mut() {
        multi_merge_sort_stable_by(ch, |&a, &b| (leaves[a], min_taxon[a]).cmp(&(leaves[b], min_taxon[b])));
    }
    // Pre-order (children in drawing order): x and the leaf order.
    let mut nodes: Vec<TreeNode> = (0..total)
        .map(|u| TreeNode {
            name: if u < n { names[u].clone() } else { String::new() },
            taxon: if u < n { Some(u) } else { None },
            parent: parent[u],
            children: children[u].clone(),
            branch_length: if u == root { 0.0 } else { blen[u] },
            x: 0.0,
            y: 0.0,
        })
        .collect();
    let mut leaf_order = Vec::with_capacity(n);
    let mut pre = Vec::with_capacity(total);
    let mut stack = vec![root];
    while let Some(u) = stack.pop() {
        pre.push(u);
        if let Some(p) = nodes[u].parent {
            nodes[u].x = nodes[p].x + nodes[u].branch_length;
        }
        if nodes[u].children.is_empty() {
            nodes[u].y = leaf_order.len() as f64;
            leaf_order.push(u);
        }
        for &c in nodes[u].children.iter().rev() {
            stack.push(c);
        }
    }
    for &u in pre.iter().rev() {
        if let (Some(&f), Some(&l)) = (nodes[u].children.first(), nodes[u].children.last()) {
            nodes[u].y = (nodes[f].y + nodes[l].y) / 2.0;
        }
    }
    // Newick, children before parents.
    let mut text: Vec<String> = vec![String::new(); total];
    for &u in pre.iter().rev() {
        let mut s = if nodes[u].children.is_empty() {
            newick_name(&nodes[u].name)
        } else {
            let parts: Vec<String> = nodes[u].children.iter().map(|&c| std::mem::take(&mut text[c])).collect();
            format!("({})", parts.join(","))
        };
        if u != root {
            s.push(':');
            s.push_str(&fmt_num(nodes[u].branch_length));
        }
        text[u] = s;
    }
    let newick = format!("{};", std::mem::take(&mut text[root]));
    let leaf_order = leaf_order.into_iter().map(|u| nodes[u].taxon.unwrap_or(u)).collect();
    Tree { method, nodes, root, leaf_order, newick }
}

/// Roots an unrooted tree in the middle of the branch of leaf `leaf`.
fn outgroup_root(adj: &mut Vec<Vec<(usize, f64)>>, leaf: usize) -> usize {
    let Some(&(v, l)) = adj[leaf].first() else { return leaf };
    let r = adj.len();
    adj.push(Vec::new());
    adj[leaf].retain(|&(x, _)| x != v);
    adj[v].retain(|&(x, _)| x != leaf);
    let h = l.max(0.0) / 2.0;
    adj[leaf].push((r, h));
    adj[r].push((leaf, h));
    adj[v].push((r, h));
    adj[r].push((v, h));
    r
}

/// The tree of a distance matrix: neighbor joining rooted at the midpoint
/// (or on the branch of `outgroup`, one of the names), or UPGMA (rooted by
/// construction: no outgroup).
pub fn build_tree(names: &[String], d: &[Vec<f64>], method: TreeMethod, outgroup: Option<&str>) -> Result<Tree, String> {
    check_matrix(names, d, "Building the tree")?;
    let n = names.len();
    let outgroup = outgroup.map(str::trim).filter(|o| !o.is_empty());
    let out_leaf = match outgroup {
        None => None,
        Some(o) => match names.iter().position(|x| x == o) {
            Some(i) => Some(i),
            None => {
                let shown: Vec<&str> = names.iter().take(10).map(String::as_str).collect();
                return Err(format!(
                    "The outgroup \"{o}\" is not one of the genomes ({}{}): give one of their names exactly, or no outgroup to root at the midpoint",
                    shown.join(", "),
                    if names.len() > 10 { ", ..." } else { "" }
                ));
            }
        },
    };
    match method {
        TreeMethod::NeighborJoining => {
            let mut adj = nj_unrooted(d);
            let root = match out_leaf {
                Some(i) if n >= 2 => outgroup_root(&mut adj, i),
                _ => midpoint_root(&mut adj, n),
            };
            // Orient from the root.
            let total = adj.len();
            let mut children = vec![Vec::new(); total];
            let mut blen = vec![0.0f64; total];
            let mut seen = vec![false; total];
            seen[root] = true;
            let mut stack = vec![root];
            while let Some(u) = stack.pop() {
                for &(v, l) in &adj[u] {
                    if !seen[v] {
                        seen[v] = true;
                        children[u].push(v);
                        blen[v] = l.max(0.0);
                        stack.push(v);
                    }
                }
            }
            Ok(finish_tree(method, names, children, blen, root))
        }
        TreeMethod::Upgma => {
            if out_leaf.is_some() {
                return Err(
                    "A UPGMA tree is rooted by construction (it assumes every lineage evolves at the same rate): an outgroup is not used. Leave the outgroup empty, or build the tree by neighbor joining (\"nj\")"
                        .to_string(),
                );
            }
            let (children, blen, root) = upgma_nodes(d);
            Ok(finish_tree(method, names, children, blen, root))
        }
    }
}

impl Tree {
    /// (parent, child) of every branch.
    pub fn edges(&self) -> Vec<(usize, usize)> {
        let mut e = Vec::new();
        for (u, n) in self.nodes.iter().enumerate() {
            for &c in &n.children {
                e.push((u, c));
            }
        }
        e
    }

    /// The distance along the tree between every two genomes (by matrix
    /// row): for an additive matrix, neighbor joining gives it back.
    pub fn patristic_distances(&self) -> Vec<Vec<f64>> {
        let leaves: Vec<(usize, usize)> = self.nodes.iter().enumerate().filter_map(|(u, n)| n.taxon.map(|t| (t, u))).collect();
        let n = leaves.len();
        let mut out = vec![vec![0.0; n]; n];
        let ancestors = |mut u: usize| {
            let mut v = vec![u];
            while let Some(p) = self.nodes[u].parent {
                v.push(p);
                u = p;
            }
            v
        };
        for a in 0..n {
            let anc_a = ancestors(leaves[a].1);
            for b in a + 1..n {
                let anc_b: std::collections::HashSet<usize> = ancestors(leaves[b].1).into_iter().collect();
                let lca = anc_a.iter().copied().find(|x| anc_b.contains(x)).unwrap_or(self.root);
                let d = self.nodes[leaves[a].1].x + self.nodes[leaves[b].1].x - 2.0 * self.nodes[lca].x;
                let (ta, tb) = (leaves[a].0, leaves[b].0);
                if ta < n && tb < n {
                    out[ta][tb] = d;
                    out[tb][ta] = d;
                }
            }
        }
        out
    }

    pub fn json(&self) -> String {
        let nodes: Vec<String> = self
            .nodes
            .iter()
            .enumerate()
            .map(|(i, n)| {
                format!(
                    "{{\"id\":{i},\"name\":{},\"taxon\":{},\"parent\":{},\"children\":[{}],\"branch_length\":{},\"x\":{},\"y\":{}}}",
                    js(&n.name),
                    n.taxon.map_or_else(|| "null".to_string(), |t| t.to_string()),
                    n.parent.map_or_else(|| "null".to_string(), |t| t.to_string()),
                    n.children.iter().map(usize::to_string).collect::<Vec<_>>().join(","),
                    jf(n.branch_length),
                    jf(n.x),
                    jf(n.y)
                )
            })
            .collect();
        let edges: Vec<String> = self.edges().iter().map(|(a, b)| format!("[{a},{b}]")).collect();
        format!(
            "{{\"method\":{},\"newick\":{},\"root\":{},\"nodes\":[{}],\"edges\":[{}],\"leaf_order\":[{}],\"depth\":{},\"text\":{}}}",
            js(self.method.name()),
            js(&self.newick),
            self.root,
            nodes.join(","),
            edges.join(","),
            self.leaf_order.iter().map(usize::to_string).collect::<Vec<_>>().join(","),
            jf(self.nodes.iter().map(|n| n.x).fold(0.0, f64::max)),
            js(&self.resumo())
        )
    }

    /// The tree as indented text, top to bottom as drawn.
    pub fn resumo(&self) -> String {
        let mut o = String::new();
        let mut stack = vec![(self.root, 0usize)];
        while let Some((u, depth)) = stack.pop() {
            let n = &self.nodes[u];
            let label = if n.children.is_empty() { n.name.clone() } else if u == self.root { "(root)".to_string() } else { "+".to_string() };
            if u == self.root {
                o.push_str(&format!("{label}\n"));
            } else {
                o.push_str(&format!("{}{label}  {}\n", "  ".repeat(depth), fmt_num(n.branch_length)));
            }
            for &c in n.children.iter().rev() {
                stack.push((c, depth + 1));
            }
        }
        o
    }
}

/// Reads a distance (or similarity) matrix: TSV or CSV with a header of
/// names (with or without the names again at the start of each row, as
/// `sourmash compare --csv` writes it), or PHYLIP (square or lower
/// triangle).
pub fn read_distance_matrix(path: &str) -> Result<(Vec<String>, Vec<Vec<f64>>), String> {
    check_input(path)?;
    let mut text = String::new();
    abre_leitura(path)
        .map_err(|e| format!("{e}. Check the path"))?
        .read_to_string(&mut text)
        .map_err(|e| format!("Error reading {path}: {e}. Is it a text table?"))?;
    let lines: Vec<&str> = text.lines().map(str::trim_end).filter(|l| !l.trim().is_empty() && !l.starts_with('#')).collect();
    if lines.is_empty() {
        return Err(format!("{path} is empty: a distance matrix has a line of names and one line per genome"));
    }
    let num = |s: &str, line: usize| -> Result<f64, String> {
        s.trim().parse::<f64>().map_err(|_| format!("{path}, line {line}: \"{}\" is not a number. Every cell after the name must be a distance", sample_text(s.trim())))
    };
    let first: Vec<&str> = lines[0].split_whitespace().collect();
    if first.len() == 1 && first[0].parse::<usize>().is_ok() {
        // PHYLIP
        let n: usize = first[0].parse().unwrap_or(0);
        if lines.len() < n + 1 {
            return Err(format!("{path}: PHYLIP says {n} genomes but has {} lines of distances", lines.len() - 1));
        }
        let mut names = Vec::new();
        let mut d = vec![vec![0.0f64; n]; n];
        for i in 0..n {
            let f: Vec<&str> = lines[i + 1].split_whitespace().collect();
            if f.is_empty() {
                return Err(format!("{path}: PHYLIP row {} is empty", i + 1));
            }
            names.push(f[0].to_string());
            let vals = &f[1..];
            if vals.len() == n {
                for (j, v) in vals.iter().enumerate() {
                    d[i][j] = num(v, i + 2)?;
                }
            } else if vals.len() == i {
                for (j, v) in vals.iter().enumerate() {
                    let x = num(v, i + 2)?;
                    d[i][j] = x;
                    d[j][i] = x;
                }
            } else {
                return Err(format!(
                    "{path}, line {}: {} has {} distances; a PHYLIP row has {n} (square) or {i} (lower triangle)",
                    i + 2,
                    f[0],
                    vals.len()
                ));
            }
        }
        return Ok((names, d));
    }
    let sep = if lines[0].contains('\t') { '\t' } else { ',' };
    let head: Vec<String> = split_fields(lines[0], sep).into_iter().map(|s| s.trim().to_string()).collect();
    let rows: Vec<Vec<String>> = lines[1..].iter().map(|l| split_fields(l, sep)).collect();
    if rows.is_empty() {
        return Err(format!("{path}: only a header line, no distances"));
    }
    let width = rows[0].len();
    let (names, labelled) = if head.len() == width && width == rows.len() {
        (head.clone(), false)
    } else if head.len() == width && width == rows.len() + 1 {
        (head[1..].to_vec(), true)
    } else {
        return Err(format!(
            "{path}: {} names in the header, {} rows of {} cells: a matrix of n genomes has n names in the header and n rows (each starting with the name, or not)",
            head.len() - usize::from(head.len() == width && width == rows.len() + 1),
            rows.len(),
            width
        ));
    };
    let n = names.len();
    let mut d = vec![vec![0.0f64; n]; n];
    for (i, r) in rows.iter().enumerate() {
        if r.len() != width {
            return Err(format!("{path}, line {}: {} cells, the first row has {width}", i + 2, r.len()));
        }
        let vals = if labelled { &r[1..] } else { &r[..] };
        for (j, v) in vals.iter().enumerate() {
            d[i][j] = num(v, i + 2)?;
        }
    }
    Ok((names, d))
}

/// The cells of a line of a table; a cell in double quotes may hold the
/// separator, and `""` stands for a quote (CSV as `sourmash compare` writes
/// it, with names such as "MN908947.3 Severe acute ..., complete genome").
fn split_fields(line: &str, sep: char) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut quoted = false;
    let mut it = line.chars().peekable();
    while let Some(c) = it.next() {
        if quoted {
            if c == '"' {
                if it.peek() == Some(&'"') {
                    cur.push('"');
                    it.next();
                } else {
                    quoted = false;
                }
            } else {
                cur.push(c);
            }
        } else if c == '"' {
            quoted = true;
        } else if c == sep {
            out.push(std::mem::take(&mut cur));
        } else {
            cur.push(c);
        }
    }
    out.push(cur);
    out
}

/// The tree of a matrix file (see `read_distance_matrix`). `similarity`:
/// the cells are similarities (Jaccard, ANI, 1 on the diagonal) and the
/// tree is built from 1 - similarity. `outgroup`: as in `build_tree`. The
/// Newick goes to `output`, if given.
pub fn tree_from_matrix(path: &str, method: TreeMethod, similarity: bool, outgroup: Option<&str>, output: Option<&str>) -> Result<Tree, String> {
    if let Some(o) = output {
        not_an_input(o, &[path])?;
    }
    let (names, mut d) = read_distance_matrix(path)?;
    if similarity {
        for row in d.iter_mut() {
            for x in row.iter_mut() {
                *x = (1.0 - *x).max(0.0);
            }
        }
    }
    check_matrix(&names, &d, path)?;
    let t = build_tree(&names, &d, method, outgroup)?;
    if let Some(o) = output {
        write_text(o, &format!("{}\n", t.newick))?;
    }
    Ok(t)
}

// =====================================================================
//  2. MUTATIONS OF A VIRAL GENOME AGAINST ITS REFERENCE
// =====================================================================

const SEED_K: usize = 17;
const SCORE_MATCH: i32 = 3;
const SCORE_MISMATCH: i32 = -1;
const GAP_OPEN: i32 = 6;
const GAP_OPEN_IN_FRAME: i32 = 7;
const GAP_OPEN_OUT_OF_FRAME: i32 = 8;
const GAP_EXTEND: i32 = 0;
const BAND: i64 = 32;
const ANCHOR_MARGIN: usize = 12;
const CHAIN_LOOKBACK: usize = 2000;
const MAX_DP_CELLS: usize = 150_000_000;

#[derive(Clone, Debug)]
pub struct MutationOptions {
    /// The GFF3 of the reference (its CDS features give the amino-acid
    /// changes). `None`: nucleotide changes only.
    pub annotation: Option<String>,
    /// One TSV per genome is written in this folder (created if needed).
    pub output_dir: Option<String>,
    /// One line per genome, Nextclade-like columns.
    pub summary_output: Option<String>,
    /// ORF1ab reported as ORF1a + ORF1b (Nextclade's names). See the top of
    /// this file.
    pub split_orf1ab: bool,
}

impl Default for MutationOptions {
    fn default() -> Self {
        MutationOptions { annotation: None, output_dir: None, summary_output: None, split_orf1ab: true }
    }
}

/// A gene (CDS) the amino-acid changes are named after.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct GeneInfo {
    pub name: String,
    pub seqid: String,
    pub strand: char,
    /// 1-based, inclusive, in coding order.
    pub parts: Vec<(u64, u64)>,
    /// Codons, the stop codon included when the CDS has one.
    pub codons: u32,
    pub product: String,
    pub protein_id: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NucSubstitution {
    pub pos: u64,
    pub ref_base: char,
    pub alt: char,
    /// The amino-acid changes it takes part in.
    pub aa: Vec<String>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NucDeletion {
    /// 1-based, inclusive.
    pub start: u64,
    pub end: u64,
    pub aa: Vec<String>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NucInsertion {
    /// The reference position the bases follow (0: before the first one).
    pub after: u64,
    pub seq: String,
    pub aa: Vec<String>,
}

/// A substitution (`alt` a residue or `*`) or a deletion (`alt` = `-`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AaChange {
    pub gene: String,
    /// 1-based codon of the gene.
    pub codon: u32,
    pub ref_aa: char,
    pub alt_aa: char,
    /// The nucleotide changes behind it.
    pub nuc: Vec<String>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AaInsertion {
    pub gene: String,
    /// The residues follow this codon (1-based; 0: before the first).
    pub after_codon: u32,
    pub seq: String,
    pub nuc: Vec<String>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FrameShift {
    pub gene: String,
    /// Codons (1-based, inclusive) read out of frame.
    pub codon_start: u32,
    pub codon_end: u32,
    /// The same stretch in reference coordinates.
    pub nuc_start: u64,
    pub nuc_end: u64,
    pub nuc: Vec<String>,
}

/// Codons that could not be read: N, or outside the sequenced part.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct UnknownAa {
    pub gene: String,
    pub codon_start: u32,
    pub codon_end: u32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct GenomeMutations {
    /// The record identifier (first word of the header).
    pub name: String,
    pub description: String,
    pub file: String,
    /// The reference record it was aligned to.
    pub reference: String,
    pub length: u64,
    /// False when it could not be aligned (`problem` says why).
    pub aligned: bool,
    pub problem: Option<String>,
    /// The genome was given as the reverse complement and turned around.
    pub reverse_complement: bool,
    /// First and last reference positions the genome covers (1-based).
    pub alignment_start: u64,
    pub alignment_end: u64,
    /// Share of the reference read as a base (not N, not missing).
    pub coverage: f64,
    pub substitutions: Vec<NucSubstitution>,
    pub deletions: Vec<NucDeletion>,
    pub insertions: Vec<NucInsertion>,
    /// Runs of N, 1-based inclusive.
    pub missing: Vec<(u64, u64)>,
    /// Bases other than A, C, G, T, N (IUPAC codes): (position, code).
    pub ambiguous: Vec<(u64, char)>,
    /// The ends of the reference the genome does not reach.
    pub not_covered: Vec<(u64, u64)>,
    pub aa_substitutions: Vec<AaChange>,
    pub aa_deletions: Vec<AaChange>,
    pub aa_insertions: Vec<AaInsertion>,
    pub frameshifts: Vec<FrameShift>,
    pub unknown_aa: Vec<UnknownAa>,
    /// The TSV written for this genome.
    pub output: Option<String>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct MutationResult {
    pub reference: String,
    /// (identifier, length) of each reference record.
    pub reference_sequences: Vec<(String, u64)>,
    pub annotation: Option<String>,
    pub genes: Vec<GeneInfo>,
    pub genomes: Vec<GenomeMutations>,
    pub output_dir: Option<String>,
    pub summary_output: Option<String>,
    pub warnings: Vec<String>,
}

impl NucSubstitution {
    pub fn label(&self) -> String {
        format!("{}{}{}", self.ref_base, self.pos, self.alt)
    }
}

impl NucDeletion {
    /// `21765-21770`, or `28271` for one base (Nextclade's form).
    pub fn short(&self) -> String {
        if self.start == self.end {
            self.start.to_string()
        } else {
            format!("{}-{}", self.start, self.end)
        }
    }
    pub fn label(&self) -> String {
        format!("del {}", self.short())
    }
}

impl NucInsertion {
    pub fn short(&self) -> String {
        format!("{}:{}", self.after, self.seq)
    }
    pub fn label(&self) -> String {
        format!("ins {}", self.short())
    }
}

impl AaChange {
    pub fn label(&self) -> String {
        format!("{}:{}{}{}", self.gene, self.ref_aa, self.codon, self.alt_aa)
    }
}

impl AaInsertion {
    pub fn label(&self) -> String {
        format!("{}:{}:{}", self.gene, self.after_codon, self.seq)
    }
}

impl FrameShift {
    pub fn label(&self) -> String {
        format!("{}:{}-{}", self.gene, self.codon_start, self.codon_end)
    }
}

impl UnknownAa {
    pub fn label(&self) -> String {
        if self.codon_start == self.codon_end {
            format!("{}:{}", self.gene, self.codon_start)
        } else {
            format!("{}:{}-{}", self.gene, self.codon_start, self.codon_end)
        }
    }
}

fn range_text(r: &(u64, u64)) -> String {
    if r.0 == r.1 {
        r.0.to_string()
    } else {
        format!("{}-{}", r.0, r.1)
    }
}

// ---------------------------------------------------------------------
//  The annotation: CDS features of a GFF3
// ---------------------------------------------------------------------

fn pct_decode(v: &str) -> String {
    if !v.contains('%') {
        return v.to_string();
    }
    let b = v.as_bytes();
    let hex = |c: u8| (c as char).to_digit(16).map(|d| d as u8);
    let mut o: Vec<u8> = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'%' && i + 2 < b.len() {
            if let (Some(x), Some(y)) = (hex(b[i + 1]), hex(b[i + 2])) {
                o.push(x * 16 + y);
                i += 3;
                continue;
            }
        }
        o.push(b[i]);
        i += 1;
    }
    String::from_utf8_lossy(&o).to_string()
}

/// Column 9: GFF3 (`key=value;`) or GTF (`key "value";`).
fn gff_attributes(col9: &str) -> Vec<(String, String)> {
    let mut v = Vec::new();
    for part in col9.split(';') {
        let p = part.trim();
        if p.is_empty() {
            continue;
        }
        match (p.find('='), p.find(' ')) {
            (Some(e), s) if s.map_or(true, |s| e < s) => v.push((p[..e].trim().to_string(), pct_decode(p[e + 1..].trim()))),
            (_, Some(s)) => v.push((p[..s].trim().to_string(), p[s + 1..].trim().trim_matches('"').to_string())),
            _ => {}
        }
    }
    v
}

fn attr<'a>(at: &'a [(String, String)], keys: &[&str]) -> Option<&'a str> {
    keys.iter().find_map(|k| at.iter().find(|(c, _)| c.eq_ignore_ascii_case(k)).map(|(_, v)| v.as_str()).filter(|v| !v.trim().is_empty()))
}

struct RawCds {
    seqid: String,
    strand: char,
    /// (start, end, phase), 1-based inclusive.
    parts: Vec<(u64, u64, u8)>,
    attrs: Vec<(String, String)>,
    line: usize,
}

/// A gene ready to be read through an alignment.
struct Gene {
    info: GeneInfo,
    rec: usize,
    /// 0-based reference positions in coding order.
    positions: Vec<u32>,
    ref_aa: Vec<u8>,
}

fn read_cds(gff: &str) -> Result<Vec<RawCds>, String> {
    check_input(gff)?;
    let rd = abre_leitura(gff).map_err(|e| format!("{e}. Check the path of the annotation"))?;
    let mut order: Vec<String> = Vec::new();
    let mut by_key: HashMap<String, RawCds> = HashMap::new();
    for (i, line) in rd.lines().enumerate() {
        let n = i + 1;
        let line = line.map_err(|e| format!("Error reading the annotation {gff} near line {n}: {e}. If it is compressed it may be damaged: get it again"))?;
        if line.starts_with("##FASTA") {
            break;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 {
            return Err(format!(
                "{gff}, line {n}: a GFF3 line has 9 columns separated by tabs, this one has {} (\"{}\"). Give the GFF3 annotation of the reference (for RefSeq, the *_genomic.gff.gz file)",
                f.len(),
                sample_text(&line)
            ));
        }
        if f[2] != "CDS" {
            continue;
        }
        let (Ok(s), Ok(e)) = (f[3].trim().parse::<u64>(), f[4].trim().parse::<u64>()) else {
            return Err(format!("{gff}, line {n}: start and end of a CDS must be whole numbers, found \"{}\" and \"{}\"", f[3], f[4]));
        };
        if s == 0 || e < s {
            return Err(format!("{gff}, line {n}: a CDS from {s} to {e} (GFF coordinates are 1-based, start <= end)"));
        }
        let strand = match f[6].trim() {
            "+" => '+',
            "-" => '-',
            o => return Err(format!("{gff}, line {n}: a CDS needs a strand, + or -, in column 7; found \"{o}\"")),
        };
        let phase = match f[7].trim() {
            "1" => 1,
            "2" => 2,
            _ => 0,
        };
        let attrs = gff_attributes(f[8]);
        let id = attr(&attrs, &["ID"]).map(str::to_string).unwrap_or_else(|| {
            match (attr(&attrs, &["Parent", "transcript_id"]), attr(&attrs, &["protein_id"])) {
                (Some(p), Some(q)) => format!("{p}|{q}"),
                (Some(p), None) => p.to_string(),
                (None, Some(q)) => q.to_string(),
                (None, None) => format!("line{n}"),
            }
        });
        let key = format!("{}\t{id}", f[0]);
        match by_key.get_mut(&key) {
            Some(r) => r.parts.push((s, e, phase)),
            None => {
                order.push(key.clone());
                by_key.insert(key, RawCds { seqid: f[0].to_string(), strand, parts: vec![(s, e, phase)], attrs, line: n });
            }
        }
    }
    if order.is_empty() {
        return Err(format!(
            "{gff} has no CDS feature (column 3 = CDS): give the GFF3 annotation of the reference genome, the one that describes its genes' coding sequences"
        ));
    }
    Ok(order.into_iter().filter_map(|k| by_key.remove(&k)).collect())
}

/// The genes of the annotation, named and split as described at the top of
/// this file, with their positions and reference proteins.
fn load_genes(gff: &str, ref_ids: &[String], ref_seqs: &[Vec<u8>], split_orf1ab: bool, warnings: &mut Vec<String>) -> Result<Vec<Gene>, String> {
    let mut raw = read_cds(gff)?;
    // Coding order of the parts.
    for r in raw.iter_mut() {
        if r.strand == '+' {
            multi_merge_sort_stable_by(&mut r.parts, |a, b| a.0.cmp(&b.0));
        } else {
            multi_merge_sort_stable_by(&mut r.parts, |a, b| b.1.cmp(&a.1));
        }
    }
    // Names.
    let base: Vec<String> = raw
        .iter()
        .map(|r| attr(&r.attrs, &["gene", "gene_name", "locus_tag", "Name", "product", "ID"]).unwrap_or("CDS").to_string())
        .collect();
    let mut names = base.clone();
    let mut groups: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, b) in base.iter().enumerate() {
        groups.entry(b.as_str()).or_default().push(i);
    }
    for (b, idx) in &groups {
        if idx.len() < 2 {
            continue;
        }
        let alts: Vec<Option<String>> = idx
            .iter()
            .map(|&i| {
                attr(&raw[i].attrs, &["product"])
                    .and_then(|p| p.split_whitespace().next())
                    .map(str::to_string)
                    .filter(|w| w.len() <= 20 && !w.eq_ignore_ascii_case("hypothetical"))
            })
            .collect();
        let distinct = alts.iter().all(Option::is_some) && {
            let mut s: Vec<&String> = alts.iter().flatten().collect();
            s.sort();
            s.dedup();
            s.len() == alts.len()
        };
        for (k, &i) in idx.iter().enumerate() {
            names[i] = if distinct { alts[k].clone().unwrap_or_default() } else { format!("{b}_{}", k + 1) };
        }
    }
    // Sequence of each CDS.
    let mut seq_of: HashMap<&str, usize> = ref_ids.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let mut gff_ids: Vec<&str> = raw.iter().map(|r| r.seqid.as_str()).collect();
    gff_ids.sort();
    gff_ids.dedup();
    if !gff_ids.iter().any(|s| seq_of.contains_key(s)) {
        if ref_ids.len() == 1 && gff_ids.len() == 1 {
            warnings.push(format!(
                "the annotation names its sequence {} and the reference {}: taken as the same sequence",
                gff_ids[0], ref_ids[0]
            ));
            seq_of.insert(gff_ids[0], 0);
        } else {
            return Err(format!(
                "The annotation {gff} describes the sequence(s) {} but the reference holds {}: give the GFF3 of this reference (the names in column 1 must be the FASTA identifiers)",
                gff_ids.join(", "),
                ref_ids.join(", ")
            ));
        }
    }
    struct Piece {
        name: String,
        raw: usize,
        parts: Vec<(u64, u64, u8)>,
        from_split: bool,
    }
    let mut pieces: Vec<Piece> = Vec::new();
    for (i, r) in raw.iter().enumerate() {
        let overlap = r.parts.windows(2).any(|w| if r.strand == '+' { w[1].0 <= w[0].1 } else { w[1].1 >= w[0].0 });
        let slippage = attr(&r.attrs, &["exception"]).is_some_and(|e| e.to_ascii_lowercase().contains("ribosomal slippage"));
        let nm = &names[i];
        if split_orf1ab && r.parts.len() == 2 && (overlap || slippage) && nm.len() > 2 && nm.to_ascii_lowercase().ends_with("ab") {
            let a = nm[..nm.len() - 1].to_string();
            let b = format!("{}{}", &nm[..nm.len() - 2], &nm[nm.len() - 1..]);
            pieces.push(Piece { name: a, raw: i, parts: vec![r.parts[0]], from_split: true });
            pieces.push(Piece { name: b, raw: i, parts: vec![r.parts[1]], from_split: true });
        } else {
            pieces.push(Piece { name: nm.clone(), raw: i, parts: r.parts.clone(), from_split: false });
        }
    }
    // A CDS with the name of a split part, starting where it starts, is the
    // non-frameshifted product of the same gene (pp1a): reported once.
    let mut drop = vec![false; pieces.len()];
    for a in 0..pieces.len() {
        if !pieces[a].from_split {
            continue;
        }
        for b in 0..pieces.len() {
            if a == b || pieces[b].from_split || pieces[b].name != pieces[a].name {
                continue;
            }
            let (ra, rb) = (&raw[pieces[a].raw], &raw[pieces[b].raw]);
            let start = |p: &Piece, r: &RawCds| if r.strand == '+' { p.parts[0].0 } else { p.parts[0].1 };
            if ra.seqid == rb.seqid && ra.strand == rb.strand && start(&pieces[a], ra) == start(&pieces[b], rb) {
                drop[b] = true;
            }
        }
    }
    let mut genes = Vec::new();
    for (pi, p) in pieces.iter().enumerate() {
        if drop[pi] {
            continue;
        }
        let r = &raw[p.raw];
        let Some(&rec) = seq_of.get(r.seqid.as_str()) else {
            warnings.push(format!("the annotation has a CDS ({}) on {}, which is not in the reference: ignored", p.name, r.seqid));
            continue;
        };
        let len = ref_seqs[rec].len() as u64;
        let mut positions: Vec<u32> = Vec::new();
        for &(s, e, _) in &p.parts {
            if e > len {
                return Err(format!(
                    "{gff}, line {}: the CDS {} ends at {e}, past the end of {} ({len} bases): this annotation is not for this reference. Give the GFF3 of the same assembly as the FASTA",
                    r.line, p.name, ref_ids[rec]
                ));
            }
            if r.strand == '+' {
                positions.extend((s - 1) as u32..e as u32);
            } else {
                positions.extend(((s - 1) as u32..e as u32).rev());
            }
        }
        let phase = usize::from(p.parts[0].2).min(positions.len());
        positions.drain(..phase);
        positions.truncate(positions.len() / 3 * 3);
        if positions.is_empty() {
            continue;
        }
        let cds: Vec<u8> = positions
            .iter()
            .map(|&q| {
                let b = ref_seqs[rec][q as usize];
                if r.strand == '+' {
                    b
                } else {
                    complement(b)
                }
            })
            .collect();
        let ref_aa = translate(&cds);
        if let Some(c) = ref_aa[..ref_aa.len().saturating_sub(1)].iter().position(|&a| a == b'*') {
            warnings.push(format!(
                "gene {} has a stop codon at codon {} of the reference itself: check that the annotation is the one of this reference",
                p.name,
                c + 1
            ));
        }
        genes.push(Gene {
            info: GeneInfo {
                name: p.name.clone(),
                seqid: ref_ids[rec].clone(),
                strand: r.strand,
                parts: p.parts.iter().map(|&(s, e, _)| (s, e)).collect(),
                codons: (positions.len() / 3) as u32,
                product: attr(&r.attrs, &["product"]).unwrap_or("").to_string(),
                protein_id: attr(&r.attrs, &["protein_id"]).unwrap_or("").to_string(),
            },
            rec,
            positions,
            ref_aa,
        });
    }
    Ok(genes)
}

// ---------------------------------------------------------------------
//  Alignment: anchors, chain, banded dynamic programming
// ---------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
struct Anchor {
    q: usize,
    r: usize,
    len: usize,
}

/// Calls `f(position, 2-bit code)` for every k-mer of A, C, G, T.
fn kmer_scan(seq: &[u8], k: usize, mut f: impl FnMut(usize, u64)) {
    let mask = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let (mut v, mut valid) = (0u64, 0usize);
    for (i, &b) in seq.iter().enumerate() {
        let c = CODE[usize::from(b)];
        if c > 3 {
            valid = 0;
            v = 0;
            continue;
        }
        v = ((v << 2) | u64::from(c)) & mask;
        valid += 1;
        if valid >= k {
            f(i + 1 - k, v);
        }
    }
}

/// Position of every k-mer of the reference; `u32::MAX` for a repeated one.
fn index_kmers(seq: &[u8], k: usize) -> HashMap<u64, u32> {
    let mut m: HashMap<u64, u32> = HashMap::with_capacity(seq.len());
    kmer_scan(seq, k, |p, v| {
        m.entry(v).and_modify(|x| *x = u32::MAX).or_insert(p as u32);
    });
    m
}

fn find_anchors(q: &[u8], idx: &HashMap<u64, u32>, k: usize) -> Vec<Anchor> {
    let mut anchors: Vec<Anchor> = Vec::new();
    let mut open: HashMap<i64, usize> = HashMap::new();
    kmer_scan(q, k, |p, v| {
        let Some(&r) = idx.get(&v) else { return };
        if r == u32::MAX {
            return;
        }
        let r = r as usize;
        let diag = r as i64 - p as i64;
        if let Some(&ai) = open.get(&diag) {
            let a = &mut anchors[ai];
            if a.q + a.len - k + 1 == p {
                a.len += 1;
                return;
            }
        }
        open.insert(diag, anchors.len());
        anchors.push(Anchor { q: p, r, len: k });
    });
    anchors
}

/// The best co-linear chain (sum of lengths minus diagonal jumps), the
/// overlaps between its anchors cut away; and the bases it covers.
fn best_chain(anchors: &[Anchor]) -> (Vec<Anchor>, usize) {
    let n = anchors.len();
    if n == 0 {
        return (Vec::new(), 0);
    }
    let mut score = vec![0i64; n];
    let mut prev = vec![usize::MAX; n];
    for i in 0..n {
        let a = anchors[i];
        score[i] = a.len as i64;
        for j in i.saturating_sub(CHAIN_LOOKBACK)..i {
            let b = anchors[j];
            if b.q >= a.q || b.r >= a.r || b.q + b.len > a.q + a.len || b.r + b.len > a.r + a.len {
                continue;
            }
            let ov = (b.q + b.len).saturating_sub(a.q).max((b.r + b.len).saturating_sub(a.r));
            if ov >= a.len {
                continue;
            }
            let jump = ((a.r as i64 - a.q as i64) - (b.r as i64 - b.q as i64)).abs();
            let s = score[j] + (a.len - ov) as i64 - jump;
            if s > score[i] {
                score[i] = s;
                prev[i] = j;
            }
        }
    }
    let mut end = 0;
    for i in 1..n {
        if score[i] > score[end] {
            end = i;
        }
    }
    let mut chain = vec![anchors[end]];
    let mut c = end;
    while prev[c] != usize::MAX {
        c = prev[c];
        chain.push(anchors[c]);
    }
    chain.reverse();
    let mut out: Vec<Anchor> = Vec::with_capacity(chain.len());
    for mut a in chain {
        if let Some(p) = out.last() {
            let ov = (p.q + p.len).saturating_sub(a.q).max((p.r + p.len).saturating_sub(a.r));
            if ov >= a.len {
                continue;
            }
            a.q += ov;
            a.r += ov;
            a.len -= ov;
        }
        out.push(a);
    }
    let covered = out.iter().map(|a| a.len).sum();
    for a in out.iter_mut() {
        let m = ANCHOR_MARGIN.min(a.len.saturating_sub(1) / 2);
        a.q += m;
        a.r += m;
        a.len -= 2 * m;
    }
    (out, covered)
}

const TB_STOP: u8 = 0;
const TB_DIAG: u8 = 1;
const TB_FROM_E: u8 = 2;
const TB_FROM_F: u8 = 3;
const TB_E_EXT: u8 = 4;
const TB_F_EXT: u8 = 8;

/// Global alignment of `q` to `r` with affine gaps in a band around the
/// diagonals of the two corners. `gap_open[j]` is the cost of a gap whose
/// first deleted reference base, or the reference base after it, is `r[j]`
/// (length `r.len() + 1`). With `free_start` / `free_end` the gaps before /
/// after the aligned part cost nothing. Returns the operations: `M` (a base
/// against a base), `D` (a reference base against a gap), `I` (a genome base
/// against a gap).
fn align_segment(q: &[u8], r: &[u8], gap_open: &[i32], free_start: bool, free_end: bool) -> Result<Vec<u8>, String> {
    let (lq, lr) = (q.len(), r.len());
    if lq == 0 {
        return Ok(vec![b'D'; lr]);
    }
    if lr == 0 {
        return Ok(vec![b'I'; lq]);
    }
    let diff = lr as i64 - lq as i64;
    let dlo = diff.min(0) - BAND;
    let dhi = diff.max(0) + BAND;
    let w = (dhi - dlo + 1) as usize;
    if (lq + 1).saturating_mul(w) > MAX_DP_CELLS {
        return Err(format!(
            "a stretch of {lq} bases of the genome against {lr} of the reference has no shared 17-mer to anchor it: the genome is too far from this reference (another virus, another genotype?) or holds a very long rearrangement"
        ));
    }
    let neg = i32::MIN / 4;
    let mut tb = vec![TB_STOP; (lq + 1) * w];
    let (mut hp, mut fp) = (vec![neg; w], vec![neg; w]);
    let (mut hc, mut ec, mut fc) = (vec![neg; w], vec![neg; w], vec![neg; w]);
    let mut best = (neg, lq, lr);
    for i in 0..=lq {
        hc.fill(neg);
        ec.fill(neg);
        fc.fill(neg);
        let jlo = (i as i64 + dlo).max(0) as usize;
        let jhi = (i as i64 + dhi).min(lr as i64) as usize;
        for j in jlo..=jhi {
            let o = (j as i64 - i as i64 - dlo) as usize;
            let cell = i * w + o;
            if (i == 0 && j == 0) || (free_start && (i == 0 || j == 0)) {
                hc[o] = 0;
                tb[cell] = TB_STOP;
                if free_end && (i == lq || j == lr) && (0 > best.0 || (0 == best.0 && i + j > best.1 + best.2)) {
                    best = (0, i, j);
                }
                continue;
            }
            let mut t = 0u8;
            let mut e = neg;
            if j > jlo {
                let open = hc[o - 1] - gap_open[j - 1];
                let ext = ec[o - 1] - GAP_EXTEND;
                if ec[o - 1] > neg && ext >= open {
                    e = ext;
                    t |= TB_E_EXT;
                } else {
                    e = open;
                }
            }
            let mut f = neg;
            if i > 0 && o + 1 < w {
                let open = hp[o + 1] - gap_open[j];
                let ext = fp[o + 1] - GAP_EXTEND;
                if fp[o + 1] > neg && ext >= open {
                    f = ext;
                    t |= TB_F_EXT;
                } else {
                    f = open;
                }
            }
            let mut h = neg;
            let mut src = TB_DIAG;
            if i > 0 && j > 0 {
                let s = if NUC_MASK[usize::from(q[i - 1])] & NUC_MASK[usize::from(r[j - 1])] != 0 { SCORE_MATCH } else { SCORE_MISMATCH };
                h = hp[o] + s;
            }
            if e > h {
                h = e;
                src = TB_FROM_E;
            }
            if f > h {
                h = f;
                src = TB_FROM_F;
            }
            hc[o] = h;
            ec[o] = e;
            fc[o] = f;
            tb[cell] = t | src;
            if free_end && (i == lq || j == lr) && (h > best.0 || (h == best.0 && i + j > best.1 + best.2)) {
                best = (h, i, j);
            }
        }
        std::mem::swap(&mut hp, &mut hc);
        std::mem::swap(&mut fp, &mut fc);
    }
    let (bi, bj) = if free_end && best.0 > neg { (best.1, best.2) } else { (lq, lr) };
    let mut ops: Vec<u8> = Vec::with_capacity(lq.max(lr) + 16);
    ops.extend(std::iter::repeat(b'D').take(lr - bj));
    ops.extend(std::iter::repeat(b'I').take(lq - bi));
    let (mut i, mut j) = (bi, bj);
    let mut state = 0u8;
    while i > 0 || j > 0 {
        let o = (j as i64 - i as i64 - dlo) as usize;
        let t = tb[i * w + o];
        match state {
            0 => match t & 3 {
                TB_STOP => break,
                TB_DIAG => {
                    ops.push(b'M');
                    i -= 1;
                    j -= 1;
                }
                TB_FROM_E => state = 1,
                _ => state = 2,
            },
            1 => {
                ops.push(b'D');
                if t & TB_E_EXT == 0 {
                    state = 0;
                }
                j -= 1;
            }
            _ => {
                ops.push(b'I');
                if t & TB_F_EXT == 0 {
                    state = 0;
                }
                i -= 1;
            }
        }
    }
    ops.extend(std::iter::repeat(b'D').take(j));
    ops.extend(std::iter::repeat(b'I').take(i));
    ops.reverse();
    Ok(ops)
}

/// The two rows of the alignment (reference above, genome below; `-` for a
/// gap) from a chain of anchors.
fn align_rows(q: &[u8], r: &[u8], gap_open: &[i32], chain: &[Anchor]) -> Result<(Vec<u8>, Vec<u8>), String> {
    let mut ops: Vec<u8> = Vec::with_capacity(q.len().max(r.len()) + 64);
    let (mut qp, mut rp) = (0usize, 0usize);
    for (i, a) in chain.iter().enumerate() {
        ops.extend(align_segment(&q[qp..a.q], &r[rp..a.r], &gap_open[rp..=a.r], i == 0, false)?);
        ops.extend(std::iter::repeat(b'M').take(a.len));
        qp = a.q + a.len;
        rp = a.r + a.len;
    }
    ops.extend(align_segment(&q[qp..], &r[rp..], &gap_open[rp..=r.len()], chain.is_empty(), true)?);
    let mut rrow = Vec::with_capacity(ops.len());
    let mut qrow = Vec::with_capacity(ops.len());
    let (mut qi, mut ri) = (0usize, 0usize);
    for op in ops {
        match op {
            b'M' => {
                rrow.push(r[ri]);
                qrow.push(q[qi]);
                ri += 1;
                qi += 1;
            }
            b'D' => {
                rrow.push(r[ri]);
                qrow.push(b'-');
                ri += 1;
            }
            _ => {
                rrow.push(b'-');
                qrow.push(q[qi]);
                qi += 1;
            }
        }
    }
    Ok((rrow, qrow))
}

fn is_clean(r: u8, q: u8) -> bool {
    r == q && is_acgt(r)
}

/// Puts every indel whose place is ambiguous at the equivalent place where
/// opening a gap costs least (`gap_open`, by reference position: outside
/// genes, then at a codon boundary, then inside a codon), the leftmost of
/// those. An indel moves only through columns where the genome equals the
/// reference (A, C, G, T).
fn place_indels(rrow: &mut [u8], qrow: &mut [u8], gap_open: &[i32]) {
    let n = rrow.len();
    let aligned = |c: usize, rr: &[u8], qq: &[u8]| rr[c] != b'-' && qq[c] != b'-';
    let Some(lo) = (0..n).find(|&c| aligned(c, rrow, qrow)) else { return };
    let Some(hi) = (0..n).rev().find(|&c| aligned(c, rrow, qrow)) else { return };
    let mut rbefore = rrow[..lo].iter().filter(|&&b| b != b'-').count();
    let mut c = lo;
    while c <= hi {
        let del = qrow[c] == b'-' && rrow[c] != b'-';
        let ins = rrow[c] == b'-' && qrow[c] != b'-';
        if !del && !ins {
            if rrow[c] != b'-' {
                rbefore += 1;
            }
            c += 1;
            continue;
        }
        let c0 = c;
        let mut c1 = c;
        while c1 <= hi && ((del && qrow[c1] == b'-' && rrow[c1] != b'-') || (ins && rrow[c1] == b'-' && qrow[c1] != b'-')) {
            c1 += 1;
        }
        let mut left = 0usize;
        let mut right = 0usize;
        {
            let row: &[u8] = if del { rrow } else { qrow };
            while c0 > lo + left + 1 && is_clean(rrow[c0 - left - 1], qrow[c0 - left - 1]) && row[c0 - left - 1] == row[c1 - left - 1] {
                left += 1;
            }
            while c1 + right < hi && is_clean(rrow[c1 + right], qrow[c1 + right]) && row[c0 + right] == row[c1 + right] {
                right += 1;
            }
        }
        let mut chosen: i64 = -(left as i64);
        let mut cheapest = i32::MAX;
        for o in -(left as i64)..=(right as i64) {
            let p = (rbefore as i64 + o) as usize;
            let cost = gap_open.get(p).copied().unwrap_or(GAP_OPEN);
            if cost < cheapest {
                cheapest = cost;
                chosen = o;
            }
        }
        let row: &mut [u8] = if del { qrow } else { rrow };
        if chosen < 0 {
            for t in 1..=(-chosen) as usize {
                row[c1 - t] = row[c0 - t];
                row[c0 - t] = b'-';
            }
        } else {
            for t in 0..chosen as usize {
                row[c0 + t] = row[c1 + t];
                row[c1 + t] = b'-';
            }
        }
        let new_end = (c1 as i64 + chosen) as usize;
        rbefore = (rbefore as i64 + chosen) as usize + if del { c1 - c0 } else { 0 };
        c = new_end;
    }
}

/// What the alignment says, position by position of the reference.
struct NucCalls {
    start: u64,
    end: u64,
    subs: Vec<NucSubstitution>,
    dels: Vec<NucDeletion>,
    ins: Vec<NucInsertion>,
    missing: Vec<(u64, u64)>,
    ambiguous: Vec<(u64, char)>,
    not_covered: Vec<(u64, u64)>,
    coverage: f64,
    /// The genome base at each reference position: the base, `-` (deleted)
    /// or 0 (not covered).
    base_at: Vec<u8>,
    /// Bases inserted after each reference position (0-based).
    ins_after: HashMap<u32, Vec<u8>>,
}

fn call_nucleotides(rrow: &[u8], qrow: &[u8], ref_len: usize) -> NucCalls {
    let n = rrow.len();
    let aligned = |c: usize| rrow[c] != b'-' && qrow[c] != b'-';
    let lo = (0..n).find(|&c| aligned(c)).unwrap_or(0);
    let hi = (0..n).rev().find(|&c| aligned(c)).unwrap_or(0);
    let mut calls = NucCalls {
        start: 0,
        end: 0,
        subs: Vec::new(),
        dels: Vec::new(),
        ins: Vec::new(),
        missing: Vec::new(),
        ambiguous: Vec::new(),
        not_covered: Vec::new(),
        coverage: 0.0,
        base_at: vec![0u8; ref_len],
        ins_after: HashMap::new(),
    };
    let mut rp = 0usize;
    let mut c = 0usize;
    let mut read = 0usize;
    while c < n {
        let (rb, qb) = (rrow[c], qrow[c]);
        if c < lo || c > hi {
            if rb != b'-' {
                rp += 1;
            }
            c += 1;
            continue;
        }
        if rb != b'-' && qb != b'-' {
            calls.base_at[rp] = qb;
            let pos = rp as u64 + 1;
            if calls.start == 0 {
                calls.start = pos;
            }
            calls.end = pos;
            if qb == b'N' {
                match calls.missing.last_mut() {
                    Some(m) if m.1 + 1 == pos => m.1 = pos,
                    _ => calls.missing.push((pos, pos)),
                }
            } else {
                read += 1;
                if !is_acgt(qb) {
                    calls.ambiguous.push((pos, char::from(qb)));
                } else if qb != rb {
                    calls.subs.push(NucSubstitution { pos, ref_base: char::from(rb), alt: char::from(qb), aa: Vec::new() });
                }
            }
            rp += 1;
            c += 1;
        } else if qb == b'-' {
            let s = rp;
            while c <= hi && qrow[c] == b'-' && rrow[c] != b'-' {
                calls.base_at[rp] = b'-';
                rp += 1;
                c += 1;
            }
            calls.dels.push(NucDeletion { start: s as u64 + 1, end: rp as u64, aa: Vec::new() });
        } else {
            let mut seq = Vec::new();
            while c <= hi && rrow[c] == b'-' && qrow[c] != b'-' {
                seq.push(qrow[c]);
                c += 1;
            }
            if rp > 0 {
                calls.ins_after.insert(rp as u32 - 1, seq.clone());
            }
            calls.ins.push(NucInsertion { after: rp as u64, seq: String::from_utf8_lossy(&seq).to_string(), aa: Vec::new() });
        }
    }
    if calls.start > 1 {
        calls.not_covered.push((1, calls.start - 1));
    }
    if calls.end > 0 && (calls.end as usize) < ref_len {
        calls.not_covered.push((calls.end + 1, ref_len as u64));
    }
    calls.coverage = if ref_len == 0 { 0.0 } else { read as f64 / ref_len as f64 };
    calls
}

// ---------------------------------------------------------------------
//  Amino acids
// ---------------------------------------------------------------------

/// Alignment of two short peptides (identity scores; X matches anything):
/// columns (reference residue, genome residue). Traced back from the end
/// preferring a residue to a gap, which puts residues to the right.
fn align_peptides(r: &[u8], q: &[u8]) -> Vec<(Option<u8>, Option<u8>)> {
    let (m, n) = (r.len(), q.len());
    let w = n + 1;
    let mut h = vec![0i32; (m + 1) * w];
    let gap = -2;
    for i in 0..=m {
        for j in 0..=n {
            if i == 0 && j == 0 {
                continue;
            }
            let mut best = i32::MIN;
            if i > 0 && j > 0 {
                let s = if q[j - 1] == b'X' {
                    0
                } else if q[j - 1] == r[i - 1] {
                    2
                } else {
                    -1
                };
                best = best.max(h[(i - 1) * w + j - 1] + s);
            }
            if i > 0 {
                best = best.max(h[(i - 1) * w + j] + gap);
            }
            if j > 0 {
                best = best.max(h[i * w + j - 1] + gap);
            }
            h[i * w + j] = best;
        }
    }
    let mut cols = Vec::with_capacity(m + n);
    let (mut i, mut j) = (m, n);
    while i > 0 || j > 0 {
        let cur = h[i * w + j];
        if i > 0 && j > 0 {
            let s = if q[j - 1] == b'X' {
                0
            } else if q[j - 1] == r[i - 1] {
                2
            } else {
                -1
            };
            if h[(i - 1) * w + j - 1] + s == cur {
                cols.push((Some(r[i - 1]), Some(q[j - 1])));
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if i > 0 && h[(i - 1) * w + j] + gap == cur {
            cols.push((Some(r[i - 1]), None));
            i -= 1;
            continue;
        }
        cols.push((None, Some(q[j - 1])));
        j -= 1;
    }
    cols.reverse();
    cols
}

struct GeneCalls {
    subs: Vec<(AaChange, (u64, u64))>,
    dels: Vec<(AaChange, (u64, u64))>,
    ins: Vec<(AaInsertion, (u64, u64))>,
    frameshifts: Vec<FrameShift>,
    unknown: Vec<UnknownAa>,
}

fn push_unknown(v: &mut Vec<UnknownAa>, gene: &str, codon: u32) {
    if let Some(u) = v.last_mut() {
        if u.codon_end + 1 == codon {
            u.codon_end = codon;
            return;
        }
    }
    v.push(UnknownAa { gene: gene.to_string(), codon_start: codon, codon_end: codon });
}

/// Reads one gene through the alignment.
fn call_gene(g: &Gene, calls: &NucCalls) -> GeneCalls {
    let pos = &g.positions;
    let minus = g.info.strand == '-';
    let len = pos.len();
    let nc = len / 3;
    let name = g.info.name.as_str();
    let base = |t: usize| -> u8 {
        let b = calls.base_at[pos[t] as usize];
        if minus && b != 0 && b != b'-' {
            complement(b)
        } else {
            b
        }
    };
    let ins = |t: usize| -> Vec<u8> {
        if t + 1 >= len {
            return Vec::new();
        }
        let (a, b) = (pos[t], pos[t + 1]);
        if !minus && b == a + 1 {
            calls.ins_after.get(&a).cloned().unwrap_or_default()
        } else if minus && b + 1 == a {
            calls.ins_after.get(&b).map(|s| revcomp(s)).unwrap_or_default()
        } else {
            Vec::new()
        }
    };
    let span = |c0: usize, c1: usize| -> (u64, u64) {
        let ps = &pos[3 * c0..3 * c1];
        let lo = ps.iter().copied().min().unwrap_or(0) as u64 + 1;
        let hi = ps.iter().copied().max().unwrap_or(0) as u64 + 1;
        (lo, hi)
    };
    let mut out = GeneCalls { subs: Vec::new(), dels: Vec::new(), ins: Vec::new(), frameshifts: Vec::new(), unknown: Vec::new() };
    let mut c = 0usize;
    while c < nc {
        let ts = [3 * c, 3 * c + 1, 3 * c + 2];
        let bs = [base(ts[0]), base(ts[1]), base(ts[2])];
        if bs.iter().any(|&b| b == 0) {
            push_unknown(&mut out.unknown, name, c as u32 + 1);
            c += 1;
            continue;
        }
        let tail = ins(ts[2]);
        let clean = bs.iter().all(|&b| b != b'-') && ins(ts[0]).is_empty() && ins(ts[1]).is_empty() && tail.len() % 3 == 0;
        if clean {
            let aa = translate_codon(&bs);
            if aa == b'X' {
                push_unknown(&mut out.unknown, name, c as u32 + 1);
            } else if aa != g.ref_aa[c] {
                out.subs.push((
                    AaChange { gene: name.to_string(), codon: c as u32 + 1, ref_aa: char::from(g.ref_aa[c]), alt_aa: char::from(aa), nuc: Vec::new() },
                    span(c, c + 1),
                ));
            }
            if !tail.is_empty() {
                let sp = if c + 1 < nc { span(c, c + 2) } else { span(c, c + 1) };
                out.ins.push((
                    AaInsertion { gene: name.to_string(), after_codon: c as u32 + 1, seq: String::from_utf8_lossy(&translate(&tail)).to_string(), nuc: Vec::new() },
                    sp,
                ));
            }
            c += 1;
            continue;
        }
        // Codons touched by indels, taken together until back in frame.
        let start = c;
        let mut qbuf: Vec<u8> = Vec::new();
        let mut unknown = false;
        let mut odd_indel = false;
        let mut run = 0usize;
        loop {
            for t in 3 * c..3 * c + 3 {
                let b = base(t);
                if b == 0 {
                    unknown = true;
                }
                if b == b'-' {
                    run += 1;
                } else {
                    if run % 3 != 0 {
                        odd_indel = true;
                    }
                    run = 0;
                    if b != 0 {
                        qbuf.push(b);
                    }
                }
                let extra = ins(t);
                if extra.len() % 3 != 0 {
                    odd_indel = true;
                }
                qbuf.extend_from_slice(&extra);
            }
            c += 1;
            if qbuf.len() % 3 == 0 || c >= nc {
                break;
            }
        }
        if run % 3 != 0 {
            odd_indel = true;
        }
        let end = c;
        if unknown {
            for k in start..end {
                push_unknown(&mut out.unknown, name, k as u32 + 1);
            }
            continue;
        }
        if odd_indel || qbuf.len() % 3 != 0 {
            let (lo, hi) = span(start, end);
            out.frameshifts.push(FrameShift {
                gene: name.to_string(),
                codon_start: start as u32 + 1,
                codon_end: end as u32,
                nuc_start: lo,
                nuc_end: hi,
                nuc: Vec::new(),
            });
            continue;
        }
        let qaa = translate(&qbuf);
        let raa = &g.ref_aa[start..end];
        let mut ci = start;
        let mut pending: Option<AaInsertion> = None;
        let flush = |p: &mut Option<AaInsertion>, out: &mut GeneCalls, ci: usize| {
            if let Some(a) = p.take() {
                let c0 = ci.saturating_sub(1);
                let c1 = (ci + 1).min(nc);
                out.ins.push((a, span(c0, c1.max(c0 + 1))));
            }
        };
        for (ra, qa) in align_peptides(raa, &qaa) {
            match (ra, qa) {
                (Some(x), Some(y)) => {
                    flush(&mut pending, &mut out, ci);
                    if y == b'X' {
                        push_unknown(&mut out.unknown, name, ci as u32 + 1);
                    } else if x != y {
                        out.subs.push((
                            AaChange { gene: name.to_string(), codon: ci as u32 + 1, ref_aa: char::from(x), alt_aa: char::from(y), nuc: Vec::new() },
                            span(start, end),
                        ));
                    }
                    ci += 1;
                }
                (Some(x), None) => {
                    flush(&mut pending, &mut out, ci);
                    out.dels.push((
                        AaChange { gene: name.to_string(), codon: ci as u32 + 1, ref_aa: char::from(x), alt_aa: '-', nuc: Vec::new() },
                        span(start, end),
                    ));
                    ci += 1;
                }
                (None, Some(y)) => match pending.as_mut() {
                    Some(p) => p.seq.push(char::from(y)),
                    None => {
                        pending = Some(AaInsertion { gene: name.to_string(), after_codon: ci as u32, seq: char::from(y).to_string(), nuc: Vec::new() })
                    }
                },
                (None, None) => {}
            }
        }
        flush(&mut pending, &mut out, ci);
    }
    out
}

// ---------------------------------------------------------------------
//  One genome
// ---------------------------------------------------------------------

struct RefData<'a> {
    ids: &'a [String],
    seqs: &'a [Vec<u8>],
    indexes: &'a [HashMap<u64, u32>],
    gap_open: &'a [Vec<i32>],
    genes: &'a [Gene],
}

fn analyze_genome(rec: &Record, file: &str, rd: &RefData) -> GenomeMutations {
    let q: Vec<u8> = rec.seq.iter().filter(|&&b| !matches!(b, b'-' | b'.' | b'*')).map(|&b| norm_base(b)).collect();
    let mut gm = GenomeMutations {
        name: if rec.id.is_empty() { file_label(file) } else { rec.id.clone() },
        description: rec.description.clone(),
        file: file.to_string(),
        length: q.len() as u64,
        ..Default::default()
    };
    if q.is_empty() {
        gm.problem = Some("the sequence is empty".to_string());
        return gm;
    }
    let qrc = revcomp(&q);
    // The reference record and the strand that share the most.
    let mut best: Option<(usize, bool, Vec<Anchor>, usize)> = None;
    for (ri, idx) in rd.indexes.iter().enumerate() {
        for rc in [false, true] {
            let qq = if rc { &qrc } else { &q };
            let (chain, covered) = best_chain(&find_anchors(qq, idx, SEED_K));
            if covered > 0 && best.as_ref().map_or(true, |b| covered > b.3) {
                best = Some((ri, rc, chain, covered));
            }
        }
    }
    let Some((ri, rc, chain, _)) = best else {
        gm.problem = Some(format!(
            "it shares no {SEED_K}-base stretch with the reference: it is not this virus (or not this segment), or it is mostly N"
        ));
        return gm;
    };
    let qq = if rc { &qrc } else { &q };
    let r = &rd.seqs[ri];
    gm.reference = rd.ids[ri].clone();
    gm.reverse_complement = rc;
    let (mut rrow, mut qrow) = match align_rows(qq, r, &rd.gap_open[ri], &chain) {
        Ok(x) => x,
        Err(e) => {
            gm.problem = Some(e);
            return gm;
        }
    };
    place_indels(&mut rrow, &mut qrow, &rd.gap_open[ri]);
    let calls = call_nucleotides(&rrow, &qrow, r.len());
    gm.aligned = true;
    gm.alignment_start = calls.start;
    gm.alignment_end = calls.end;
    gm.coverage = calls.coverage;
    let mut subs = Vec::new();
    let mut dels = Vec::new();
    let mut ins = Vec::new();
    let mut fs = Vec::new();
    for g in rd.genes.iter().filter(|g| g.rec == ri) {
        let gc = call_gene(g, &calls);
        subs.extend(gc.subs);
        dels.extend(gc.dels);
        ins.extend(gc.ins);
        fs.extend(gc.frameshifts);
        gm.unknown_aa.extend(gc.unknown);
    }
    // Links between nucleotide and amino-acid changes.
    let mut nsubs = calls.subs;
    let mut ndels = calls.dels;
    let mut nins = calls.ins;
    let touches_sub = |p: u64, sp: (u64, u64)| p >= sp.0 && p <= sp.1;
    let touches_del = |d: &NucDeletion, sp: (u64, u64)| d.start <= sp.1 && d.end >= sp.0;
    let touches_ins = |i: &NucInsertion, sp: (u64, u64)| i.after + 1 >= sp.0 && i.after <= sp.1;
    let link = |label: String, sp: (u64, u64), nuc: &mut Vec<String>, nsubs: &mut Vec<NucSubstitution>, ndels: &mut Vec<NucDeletion>, nins: &mut Vec<NucInsertion>| {
        for s in nsubs.iter_mut().filter(|s| touches_sub(s.pos, sp)) {
            nuc.push(s.label());
            s.aa.push(label.clone());
        }
        for d in ndels.iter_mut().filter(|d| touches_del(d, sp)) {
            nuc.push(d.label());
            d.aa.push(label.clone());
        }
        for i in nins.iter_mut().filter(|i| touches_ins(i, sp)) {
            nuc.push(i.label());
            i.aa.push(label.clone());
        }
    };
    for (mut a, sp) in subs {
        let l = a.label();
        link(l, sp, &mut a.nuc, &mut nsubs, &mut ndels, &mut nins);
        gm.aa_substitutions.push(a);
    }
    for (mut a, sp) in dels {
        let l = a.label();
        link(l, sp, &mut a.nuc, &mut nsubs, &mut ndels, &mut nins);
        gm.aa_deletions.push(a);
    }
    for (mut a, sp) in ins {
        let l = a.label();
        link(l, sp, &mut a.nuc, &mut nsubs, &mut ndels, &mut nins);
        gm.aa_insertions.push(a);
    }
    for mut f in fs {
        let l = f.label();
        let sp = (f.nuc_start, f.nuc_end);
        link(l, sp, &mut f.nuc, &mut nsubs, &mut ndels, &mut nins);
        gm.frameshifts.push(f);
    }
    gm.substitutions = nsubs;
    gm.deletions = ndels;
    gm.insertions = nins;
    gm.missing = calls.missing;
    gm.ambiguous = calls.ambiguous;
    gm.not_covered = calls.not_covered;
    gm
}

/// A file name from a genome name.
fn safe_file_name(s: &str) -> String {
    let t: String = s.chars().map(|c| if c.is_ascii_alphanumeric() || matches!(c, '.' | '_' | '-') { c } else { '_' }).collect();
    let t = t.trim_matches('.').to_string();
    if t.is_empty() {
        "genome".to_string()
    } else {
        t
    }
}

/// Every genome of the FASTA files `genomes` (each record a genome) against
/// the reference FASTA: nucleotide changes, and amino-acid changes when the
/// reference GFF3 is given.
pub fn viral_mutations(reference: &str, genomes: &[String], opt: &MutationOptions) -> Result<MutationResult, String> {
    check_input(reference)?;
    if genomes.is_empty() {
        return Err("No genome given: give one or more FASTA files of genomes (each record is a genome)".to_string());
    }
    let mut inputs: Vec<&str> = vec![reference];
    for g in genomes {
        check_input(g)?;
        inputs.push(g);
    }
    if let Some(a) = &opt.annotation {
        check_input(a)?;
        inputs.push(a);
    }
    if let Some(s) = &opt.summary_output {
        not_an_input(s, &inputs)?;
    }
    let mut warnings = Vec::new();
    let (refs, kind) = read_records(reference)?;
    if kind == Some(SeqKind::Fastq) {
        return Err(format!("The reference {reference} is a FASTQ file: give the reference genome as FASTA"));
    }
    if refs.is_empty() {
        return Err(format!("The reference {reference} holds no sequence: give the reference genome as FASTA"));
    }
    let ids: Vec<String> = refs.iter().enumerate().map(|(i, r)| if r.id.is_empty() { format!("sequence{}", i + 1) } else { r.id.clone() }).collect();
    let seqs: Vec<Vec<u8>> = refs.iter().map(|r| r.seq.iter().filter(|&&b| !matches!(b, b'-' | b'.' | b'*')).map(|&b| norm_base(b)).collect()).collect();
    drop(refs);
    let genes = match &opt.annotation {
        Some(a) => load_genes(a, &ids, &seqs, opt.split_orf1ab, &mut warnings)?,
        None => Vec::new(),
    };
    // The cost of opening a gap along each reference record: 6 outside
    // genes, 7 at a codon boundary, 8 inside a codon.
    let mut gap_open: Vec<Vec<i32>> = seqs.iter().map(|s| vec![GAP_OPEN; s.len() + 1]).collect();
    for g in &genes {
        for &p in &g.positions {
            if gap_open[g.rec][p as usize] == GAP_OPEN {
                gap_open[g.rec][p as usize] = GAP_OPEN_OUT_OF_FRAME;
            }
        }
    }
    for g in &genes {
        for codon in g.positions.chunks_exact(3) {
            let lo = codon.iter().copied().min().unwrap_or(0) as usize;
            gap_open[g.rec][lo] = GAP_OPEN_IN_FRAME;
        }
    }
    let indexes: Vec<HashMap<u64, u32>> = seqs.par_iter().map(|s| index_kmers(s, SEED_K)).collect();
    let rd = RefData { ids: &ids, seqs: &seqs, indexes: &indexes, gap_open: &gap_open, genes: &genes };
    // The genomes are read in batches and aligned in parallel: only the
    // batch is in memory, so a file of any number of genomes is fine.
    let mut results: Vec<GenomeMutations> = Vec::new();
    for g in genomes {
        let mut rdr = SeqReader::open(g)?;
        let mut rec = SeqRecord::default();
        let mut batch: Vec<Record> = Vec::new();
        let mut seen = 0usize;
        loop {
            let more = rdr.next(&mut rec)?;
            if more {
                if rdr.kind() == Some(SeqKind::Fastq) {
                    return Err(format!("{g} is a FASTQ file (reads): give assembled genomes (FASTA, one record per genome); for reads, build a consensus first"));
                }
                let head = String::from_utf8_lossy(&rec.head).trim().to_string();
                let id = head.split_whitespace().next().unwrap_or("").to_string();
                batch.push(Record { id, description: head, seq: std::mem::take(&mut rec.seq) });
                seen += 1;
            }
            if batch.len() >= 256 || (!more && !batch.is_empty()) {
                let done: Vec<GenomeMutations> = batch.par_iter().map(|r| analyze_genome(r, g, &rd)).collect();
                results.extend(done);
                batch.clear();
            }
            if !more {
                break;
            }
        }
        if seen == 0 {
            warnings.push(format!("{g} holds no sequence"));
        }
    }
    let mut names: Vec<String> = results.iter().map(|g| g.name.clone()).collect();
    make_unique(&mut names, &mut warnings);
    for (g, n) in results.iter_mut().zip(names) {
        g.name = n;
    }
    for g in &results {
        if let Some(p) = &g.problem {
            warnings.push(format!("{} could not be aligned: {p}", g.name));
        }
    }
    let mut result = MutationResult {
        reference: reference.to_string(),
        reference_sequences: ids.iter().zip(&seqs).map(|(i, s)| (i.clone(), s.len() as u64)).collect(),
        annotation: opt.annotation.clone(),
        genes: genes.iter().map(|g| g.info.clone()).collect(),
        genomes: results,
        output_dir: opt.output_dir.clone(),
        summary_output: opt.summary_output.clone(),
        warnings,
    };
    let mut cleanup = Cleanup(Vec::new());
    if let Some(dir) = &opt.output_dir {
        std::fs::create_dir_all(dir).map_err(|e| format!("Cannot create the output folder {dir}: {e}. Check the path and the permission to write there"))?;
        let mut used: HashMap<String, usize> = HashMap::new();
        for g in result.genomes.iter_mut() {
            let base = safe_file_name(&g.name);
            let k = used.entry(base.to_ascii_lowercase()).or_insert(0);
            *k += 1;
            let fname = if *k == 1 { format!("{base}.tsv") } else { format!("{base}_{k}.tsv") };
            let path = std::path::Path::new(dir).join(fname).to_string_lossy().to_string();
            not_an_input(&path, &inputs)?;
            cleanup.0.push(path.clone());
            write_text(&path, &g.tsv())?;
            g.output = Some(path);
        }
    }
    if let Some(p) = &opt.summary_output {
        cleanup.0.push(p.clone());
        write_text(p, &result.summary_tsv())?;
    }
    cleanup.disarm();
    Ok(result)
}

impl GenomeMutations {
    /// Every change, one per line (what `output_dir` receives).
    pub fn tsv(&self) -> String {
        let mut o = String::from("kind\tlabel\tstart\tend\tref\talt\tgene\tlinked\n");
        let mut line = |kind: &str, label: String, a: u64, b: u64, r: &str, alt: &str, gene: &str, linked: &[String]| {
            o.push_str(&format!("{kind}\t{}\t{a}\t{b}\t{r}\t{alt}\t{gene}\t{}\n", clean_cell(&label), clean_cell(&linked.join(","))));
        };
        for s in &self.substitutions {
            line("substitution", s.label(), s.pos, s.pos, &s.ref_base.to_string(), &s.alt.to_string(), "", &s.aa);
        }
        for d in &self.deletions {
            line("deletion", d.label(), d.start, d.end, "", "-", "", &d.aa);
        }
        for i in &self.insertions {
            line("insertion", i.label(), i.after, i.after, "-", &i.seq, "", &i.aa);
        }
        for m in &self.missing {
            line("missing", format!("missing {}", range_text(m)), m.0, m.1, "", "N", "", &[]);
        }
        for (p, c) in &self.ambiguous {
            line("ambiguous", format!("{c}:{p}"), *p, *p, "", &c.to_string(), "", &[]);
        }
        for m in &self.not_covered {
            line("not_covered", format!("not covered {}", range_text(m)), m.0, m.1, "", "", "", &[]);
        }
        for a in &self.aa_substitutions {
            line("aa_substitution", a.label(), u64::from(a.codon), u64::from(a.codon), &a.ref_aa.to_string(), &a.alt_aa.to_string(), &a.gene, &a.nuc);
        }
        for a in &self.aa_deletions {
            line("aa_deletion", a.label(), u64::from(a.codon), u64::from(a.codon), &a.ref_aa.to_string(), "-", &a.gene, &a.nuc);
        }
        for a in &self.aa_insertions {
            line("aa_insertion", a.label(), u64::from(a.after_codon), u64::from(a.after_codon), "-", &a.seq, &a.gene, &a.nuc);
        }
        for f in &self.frameshifts {
            line("frameshift", f.label(), u64::from(f.codon_start), u64::from(f.codon_end), "", "", &f.gene, &f.nuc);
        }
        for u in &self.unknown_aa {
            line("unknown_aa", u.label(), u64::from(u.codon_start), u64::from(u.codon_end), "", "X", &u.gene, &[]);
        }
        o
    }

    pub fn json(&self) -> String {
        let strs = |v: Vec<String>| jlist(&v);
        let subs: Vec<String> = self
            .substitutions
            .iter()
            .map(|s| format!("{{\"pos\":{},\"ref\":{},\"alt\":{},\"label\":{},\"aa\":{}}}", s.pos, js(&s.ref_base.to_string()), js(&s.alt.to_string()), js(&s.label()), jlist(&s.aa)))
            .collect();
        let dels: Vec<String> = self
            .deletions
            .iter()
            .map(|d| format!("{{\"start\":{},\"end\":{},\"label\":{},\"aa\":{}}}", d.start, d.end, js(&d.short()), jlist(&d.aa)))
            .collect();
        let ins: Vec<String> = self
            .insertions
            .iter()
            .map(|i| format!("{{\"after\":{},\"seq\":{},\"label\":{},\"aa\":{}}}", i.after, js(&i.seq), js(&i.short()), jlist(&i.aa)))
            .collect();
        let aa = |v: &[AaChange]| -> String {
            let x: Vec<String> = v
                .iter()
                .map(|a| {
                    format!(
                        "{{\"gene\":{},\"codon\":{},\"ref\":{},\"alt\":{},\"label\":{},\"nuc\":{}}}",
                        js(&a.gene),
                        a.codon,
                        js(&a.ref_aa.to_string()),
                        js(&a.alt_aa.to_string()),
                        js(&a.label()),
                        jlist(&a.nuc)
                    )
                })
                .collect();
            format!("[{}]", x.join(","))
        };
        let aains: Vec<String> = self
            .aa_insertions
            .iter()
            .map(|a| format!("{{\"gene\":{},\"after_codon\":{},\"seq\":{},\"label\":{},\"nuc\":{}}}", js(&a.gene), a.after_codon, js(&a.seq), js(&a.label()), jlist(&a.nuc)))
            .collect();
        let fsj: Vec<String> = self
            .frameshifts
            .iter()
            .map(|f| {
                format!(
                    "{{\"gene\":{},\"codon_start\":{},\"codon_end\":{},\"nuc_start\":{},\"nuc_end\":{},\"label\":{},\"nuc\":{}}}",
                    js(&f.gene),
                    f.codon_start,
                    f.codon_end,
                    f.nuc_start,
                    f.nuc_end,
                    js(&f.label()),
                    jlist(&f.nuc)
                )
            })
            .collect();
        let ranges = |v: &[(u64, u64)]| format!("[{}]", v.iter().map(|r| format!("[{},{}]", r.0, r.1)).collect::<Vec<_>>().join(","));
        format!(
            "{{\"name\":{},\"description\":{},\"file\":{},\"reference\":{},\"length\":{},\"aligned\":{},\"problem\":{},\"reverse_complement\":{},\"alignment_start\":{},\"alignment_end\":{},\"coverage\":{},\"substitutions\":[{}],\"deletions\":[{}],\"insertions\":[{}],\"missing\":{},\"ambiguous\":{},\"not_covered\":{},\"aa_substitutions\":{},\"aa_deletions\":{},\"aa_insertions\":[{}],\"frameshifts\":[{}],\"unknown_aa\":{},\"output\":{},\"text\":{}}}",
            js(&self.name),
            js(&self.description),
            js(&self.file),
            js(&self.reference),
            self.length,
            self.aligned,
            jopt(&self.problem),
            self.reverse_complement,
            self.alignment_start,
            self.alignment_end,
            jf(self.coverage),
            subs.join(","),
            dels.join(","),
            ins.join(","),
            ranges(&self.missing),
            strs(self.ambiguous.iter().map(|(p, c)| format!("{c}:{p}")).collect()),
            ranges(&self.not_covered),
            aa(&self.aa_substitutions),
            aa(&self.aa_deletions),
            aains.join(","),
            fsj.join(","),
            strs(self.unknown_aa.iter().map(UnknownAa::label).collect()),
            jopt(&self.output),
            js(&self.resumo())
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!("{}  ({} bases)\n", self.name, thousands(self.length));
        if let Some(p) = &self.problem {
            o.push_str(&format!("  not aligned: {p}\n"));
            return o;
        }
        o.push_str(&format!(
            "  covers {} to {} of {}{}, {:.2}% of the reference read\n",
            self.alignment_start,
            self.alignment_end,
            self.reference,
            if self.reverse_complement { " (given reverse-complemented, turned around)" } else { "" },
            100.0 * self.coverage
        ));
        let deleted: u64 = self.deletions.iter().map(|d| d.end - d.start + 1).sum();
        let missing: u64 = self.missing.iter().map(|m| m.1 - m.0 + 1).sum();
        o.push_str(&format!(
            "  {} substitution(s), {} deletion(s) ({} bases), {} insertion(s), {} missing base(s) (N), {} ambiguous base(s)\n",
            self.substitutions.len(),
            self.deletions.len(),
            deleted,
            self.insertions.len(),
            missing,
            self.ambiguous.len()
        ));
        let list = |v: Vec<String>| if v.is_empty() { "-".to_string() } else { v.join(", ") };
        o.push_str(&format!("  substitutions: {}\n", list(self.substitutions.iter().map(NucSubstitution::label).collect())));
        if !self.deletions.is_empty() {
            o.push_str(&format!("  deletions: {}\n", list(self.deletions.iter().map(NucDeletion::short).collect())));
        }
        if !self.insertions.is_empty() {
            o.push_str(&format!("  insertions: {}\n", list(self.insertions.iter().map(NucInsertion::short).collect())));
        }
        if !self.missing.is_empty() {
            o.push_str(&format!("  missing: {}\n", list(self.missing.iter().map(range_text).collect())));
        }
        o.push_str(&format!(
            "  amino acids: {}\n",
            list(
                self.aa_substitutions
                    .iter()
                    .chain(&self.aa_deletions)
                    .map(AaChange::label)
                    .chain(self.aa_insertions.iter().map(AaInsertion::label))
                    .collect()
            )
        ));
        if !self.frameshifts.is_empty() {
            o.push_str(&format!("  frameshifts: {}\n", list(self.frameshifts.iter().map(FrameShift::label).collect())));
        }
        o
    }
}

impl MutationResult {
    /// One line per genome (what `summary_output` receives).
    pub fn summary_tsv(&self) -> String {
        let mut o = String::from(
            "genome\treference\tlength\talignment_start\talignment_end\tcoverage\tsubstitution_count\tdeletion_count\tdeleted_bases\tinsertion_count\tmissing_bases\taa_substitution_count\taa_deletion_count\taa_insertion_count\tframeshift_count\tsubstitutions\tdeletions\tinsertions\tmissing\tnon_acgtn\taa_substitutions\taa_deletions\taa_insertions\tframeshifts\treverse_complement\tproblem\n",
        );
        for g in &self.genomes {
            let deleted: u64 = g.deletions.iter().map(|d| d.end - d.start + 1).sum();
            let missing: u64 = g.missing.iter().map(|m| m.1 - m.0 + 1).sum();
            let j = |v: Vec<String>| clean_cell(&v.join(","));
            o.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                clean_cell(&g.name),
                clean_cell(&g.reference),
                g.length,
                g.alignment_start,
                g.alignment_end,
                g.coverage,
                g.substitutions.len(),
                g.deletions.len(),
                deleted,
                g.insertions.len(),
                missing,
                g.aa_substitutions.len(),
                g.aa_deletions.len(),
                g.aa_insertions.len(),
                g.frameshifts.len(),
                j(g.substitutions.iter().map(NucSubstitution::label).collect()),
                j(g.deletions.iter().map(NucDeletion::short).collect()),
                j(g.insertions.iter().map(NucInsertion::short).collect()),
                j(g.missing.iter().map(range_text).collect()),
                j(g.ambiguous.iter().map(|(p, c)| format!("{c}:{p}")).collect()),
                j(g.aa_substitutions.iter().map(AaChange::label).collect()),
                j(g.aa_deletions.iter().map(AaChange::label).collect()),
                j(g.aa_insertions.iter().map(AaInsertion::label).collect()),
                j(g.frameshifts.iter().map(FrameShift::label).collect()),
                g.reverse_complement,
                clean_cell(g.problem.as_deref().unwrap_or(""))
            ));
        }
        o
    }

    pub fn json(&self) -> String {
        let genes: Vec<String> = self
            .genes
            .iter()
            .map(|g| {
                format!(
                    "{{\"name\":{},\"seqid\":{},\"strand\":{},\"parts\":[{}],\"codons\":{},\"product\":{},\"protein_id\":{}}}",
                    js(&g.name),
                    js(&g.seqid),
                    js(&g.strand.to_string()),
                    g.parts.iter().map(|p| format!("[{},{}]", p.0, p.1)).collect::<Vec<_>>().join(","),
                    g.codons,
                    js(&g.product),
                    js(&g.protein_id)
                )
            })
            .collect();
        let refs: Vec<String> = self.reference_sequences.iter().map(|(i, l)| format!("{{\"id\":{},\"length\":{l}}}", js(i))).collect();
        let gs: Vec<String> = self.genomes.iter().map(GenomeMutations::json).collect();
        format!(
            "{{\"reference\":{},\"reference_sequences\":[{}],\"annotation\":{},\"genes\":[{}],\"genomes\":[{}],\"output_dir\":{},\"summary_output\":{},\"warnings\":{},\"text\":{}}}",
            js(&self.reference),
            refs.join(","),
            jopt(&self.annotation),
            genes.join(","),
            gs.join(","),
            jopt(&self.output_dir),
            jopt(&self.summary_output),
            jlist(&self.warnings),
            js(&self.resumo())
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!(
            "Mutations against {} ({}), {} genome(s){}\n",
            self.reference,
            self.reference_sequences.iter().map(|(i, l)| format!("{i}, {} bases", thousands(*l))).collect::<Vec<_>>().join("; "),
            self.genomes.len(),
            match &self.annotation {
                Some(a) => format!(", {} gene(s) from {a}", self.genes.len()),
                None => ", no annotation (nucleotide changes only)".to_string(),
            }
        );
        for g in &self.genomes {
            for l in g.resumo().lines() {
                o.push_str("  ");
                o.push_str(l);
                o.push('\n');
            }
        }
        if let Some(d) = &self.output_dir {
            o.push_str(&format!("  one table per genome written to {d}\n"));
        }
        if let Some(p) = &self.summary_output {
            o.push_str(&format!("  summary written to {p}\n"));
        }
        for w in &self.warnings {
            o.push_str(&format!("  warning: {w}\n"));
        }
        o
    }
}

// =====================================================================
//  3. IN-SILICO PCR
// =====================================================================

/// One assay of a panel: a primer pair, an optional probe, the expected
/// product size. Sequences 5' -> 3', IUPAC codes allowed.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PcrAssay {
    pub target: String,
    pub forward: String,
    pub reverse: String,
    pub probe: Option<String>,
    pub min_size: Option<u32>,
    pub max_size: Option<u32>,
    /// Where the sequences come from (free text, kept in the reports).
    pub source: String,
}

#[derive(Clone, Debug)]
pub struct PcrOptions {
    /// Mismatches allowed in each primer (3).
    pub max_mismatches: u32,
    /// Bases at the 3' end of a primer that must match (3).
    pub three_prime_exact: u32,
    /// Mismatches allowed in the two primers of a product together (4).
    pub max_total_mismatches: u32,
    /// Mismatches allowed in the probe (3).
    pub max_probe_mismatches: u32,
    /// Longest product searched for when the assay gives no larger maximum
    /// size (2000).
    pub max_amplicon: u32,
    /// Each FASTA record is a genome (a multi-FASTA of viral genomes);
    /// otherwise each FASTA file is one genome.
    pub records_as_genomes: bool,
    /// Reads carrying a site of the assay needed to call a read sample
    /// positive (3).
    pub min_reads: u64,
    pub amplicons_output: Option<String>,
    pub matrix_output: Option<String>,
    pub reads_output: Option<String>,
}

impl Default for PcrOptions {
    fn default() -> Self {
        PcrOptions {
            max_mismatches: 3,
            three_prime_exact: 3,
            max_total_mismatches: 4,
            max_probe_mismatches: 3,
            max_amplicon: 2000,
            records_as_genomes: false,
            min_reads: 3,
            amplicons_output: None,
            matrix_output: None,
            reads_output: None,
        }
    }
}

/// One product.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Amplicon {
    /// The row of the matrix (file or record).
    pub genome: String,
    /// The FASTA record it lies on.
    pub sequence_id: String,
    pub target: String,
    /// `+`: the forward primer reads along the sequence as written; `-`:
    /// along its reverse complement.
    pub strand: char,
    /// On the sequence as written, 1-based, inclusive, primers included.
    pub start: u64,
    pub end: u64,
    pub size: u64,
    pub forward_mismatches: u32,
    pub reverse_mismatches: u32,
    /// `None`: the assay has no probe, or the probe does not bind.
    pub probe_mismatches: Option<u32>,
    /// The size is in the assay's range.
    pub in_range: bool,
    /// In range, and the probe (if any) binds: a positive result.
    pub detected: bool,
    /// The product, from the forward primer to the reverse one.
    pub sequence: String,
}

/// Reads of one sample against one assay.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReadTargetCount {
    pub sample: String,
    pub target: String,
    pub reads: u64,
    pub forward: u64,
    pub reverse: u64,
    pub probe: u64,
    /// Reads holding a whole product (both primers, facing each other).
    pub amplicon: u64,
    /// Reads holding at least one site of the assay.
    pub with_site: u64,
    pub present: bool,
}

/// One row of the presence/absence matrix.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PcrRow {
    pub name: String,
    /// "genome" or "reads".
    pub kind: String,
    pub file: String,
    pub present: Vec<bool>,
    /// Detected products (genome) or reads with a site (reads), per target.
    pub counts: Vec<u64>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PcrResult {
    pub assays: Vec<PcrAssay>,
    pub rows: Vec<PcrRow>,
    pub amplicons: Vec<Amplicon>,
    pub read_counts: Vec<ReadTargetCount>,
    pub max_mismatches: u32,
    pub three_prime_exact: u32,
    pub max_total_mismatches: u32,
    pub max_probe_mismatches: u32,
    pub max_amplicon: u32,
    pub min_reads: u64,
    pub amplicons_output: Option<String>,
    pub matrix_output: Option<String>,
    pub reads_output: Option<String>,
    pub warnings: Vec<String>,
}

/// Products of one assay on one sequence beyond this number are not kept.
const MAX_PRODUCTS: usize = 10_000;

fn clean_oligo(s: &str, what: &str, target: &str, whence: &str) -> Result<String, String> {
    let t: String = s.chars().filter(|c| !c.is_whitespace() && *c != '-').collect::<String>().to_ascii_uppercase().replace('U', "T");
    if t.is_empty() {
        return Err(format!("{whence}: the {what} of {target} is empty: give its sequence 5' to 3'"));
    }
    if let Some(c) = t.chars().find(|&c| !c.is_ascii() || NUC_MASK[c as usize] == 0) {
        return Err(format!(
            "{whence}: the {what} of {target} holds \"{c}\", which is not a nucleotide code. Use A, C, G, T and the IUPAC codes (R Y S W K M B D H V N), 5' to 3'"
        ));
    }
    if t.len() < 8 {
        return Err(format!("{whence}: the {what} of {target} has {} bases; a primer or probe has at least 8 (usually 18 to 30)", t.len()));
    }
    Ok(t)
}

fn checked_assay(a: &PcrAssay, whence: &str) -> Result<PcrAssay, String> {
    let target = a.target.trim();
    if target.is_empty() {
        return Err(format!("{whence}: an assay has no target name: give each primer pair a name (the column of the matrix)"));
    }
    let forward = clean_oligo(&a.forward, "forward primer", target, whence)?;
    let reverse = clean_oligo(&a.reverse, "reverse primer", target, whence)?;
    let probe = match a.probe.as_deref().map(str::trim) {
        None | Some("") | Some("-") | Some("NA") => None,
        Some(p) => Some(clean_oligo(p, "probe", target, whence)?),
    };
    if let (Some(lo), Some(hi)) = (a.min_size, a.max_size) {
        if lo > hi {
            return Err(format!("{whence}: the size range of {target} goes from {lo} to {hi}: the minimum must not exceed the maximum"));
        }
    }
    Ok(PcrAssay { target: target.to_string(), forward, reverse, probe, min_size: a.min_size, max_size: a.max_size, source: a.source.trim().to_string() })
}

/// Reads a panel: a table with a header (tab-separated, or comma-separated
/// if the header has no tab), one assay per line. Columns, by name in any
/// order: `target` (or name, assay), `forward` (fwd, forward_primer),
/// `reverse` (rev, reverse_primer), optional `probe`, `min_size` and
/// `max_size` (or `size` as "60-80" or "72"), `source`. Lines starting with
/// `#` are comments.
pub fn read_pcr_panel(path: &str) -> Result<Vec<PcrAssay>, String> {
    check_input(path)?;
    let rd = abre_leitura(path).map_err(|e| format!("{e}. Check the path of the panel"))?;
    let mut head: Option<(Vec<String>, char)> = None;
    let mut out = Vec::new();
    for (i, line) in rd.lines().enumerate() {
        let n = i + 1;
        let line = line.map_err(|e| format!("Error reading the panel {path} near line {n}: {e}"))?;
        let t = line.trim_end_matches(['\r', '\n']);
        if t.trim().is_empty() || t.trim_start().starts_with('#') {
            continue;
        }
        let Some((cols, sep)) = &head else {
            let sep = if t.contains('\t') { '\t' } else { ',' };
            head = Some((t.split(sep).map(|c| c.trim().trim_matches('"').to_ascii_lowercase()).collect(), sep));
            continue;
        };
        let cells: Vec<&str> = t.split(*sep).map(|c| c.trim().trim_matches('"')).collect();
        let get = |names: &[&str]| -> Option<&str> {
            names.iter().find_map(|nm| cols.iter().position(|c| c == nm)).and_then(|p| cells.get(p).copied()).filter(|s| !s.is_empty())
        };
        let whence = format!("{path}, line {n}");
        let target = get(&["target", "name", "assay", "gene"]);
        let forward = get(&["forward", "fwd", "forward_primer", "f", "primer_f", "left"]);
        let reverse = get(&["reverse", "rev", "reverse_primer", "r", "primer_r", "right"]);
        let (Some(target), Some(forward), Some(reverse)) = (target, forward, reverse) else {
            return Err(format!(
                "{whence}: an assay needs a target, a forward and a reverse primer. The header must name the columns: target, forward, reverse (optional: probe, min_size, max_size or size, source); this header has: {}",
                cols.join(", ")
            ));
        };
        let num = |s: &str| -> Result<u32, String> {
            s.trim().parse::<u32>().map_err(|_| format!("{whence}: \"{}\" is not a size in bases (a whole number)", sample_text(s)))
        };
        let (mut min_size, mut max_size) = (None, None);
        if let Some(s) = get(&["size", "expected_size", "size_range", "amplicon_size"]) {
            match s.split_once('-') {
                Some((a, b)) => {
                    min_size = Some(num(a)?);
                    max_size = Some(num(b)?);
                }
                None => {
                    let v = num(s)?;
                    min_size = Some(v);
                    max_size = Some(v);
                }
            }
        }
        if let Some(s) = get(&["min_size", "min", "minimum_size"]) {
            min_size = Some(num(s)?);
        }
        if let Some(s) = get(&["max_size", "max", "maximum_size"]) {
            max_size = Some(num(s)?);
        }
        let a = PcrAssay {
            target: target.to_string(),
            forward: forward.to_string(),
            reverse: reverse.to_string(),
            probe: get(&["probe", "probe_sequence"]).map(str::to_string),
            min_size,
            max_size,
            source: get(&["source", "reference", "citation"]).unwrap_or("").to_string(),
        };
        out.push(checked_assay(&a, &whence)?);
    }
    if out.is_empty() {
        return Err(format!("{path} holds no assay: give a header line (target, forward, reverse, ...) and one line per primer pair"));
    }
    Ok(out)
}

fn template_masks(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| TEMPLATE_MASK[usize::from(b)]).collect()
}

fn comp_mask(m: u8) -> u8 {
    ((m & 1) << 3) | ((m & 8) >> 3) | ((m & 2) << 1) | ((m & 4) >> 1)
}

fn rc_masks(t: &[u8]) -> Vec<u8> {
    t.iter().rev().map(|&m| comp_mask(m)).collect()
}

fn oligo_masks(s: &str) -> Vec<u8> {
    s.bytes().map(|b| NUC_MASK[usize::from(b)]).collect()
}

/// Every place where the oligo `p` binds `t` as written (both as masks):
/// (start, mismatches); no mismatch in the last `exact3` bases.
fn sites(t: &[u8], p: &[u8], max_mm: u32, exact3: usize) -> Vec<(usize, u32)> {
    let (n, l) = (t.len(), p.len());
    let mut out = Vec::new();
    if l == 0 || n < l {
        return out;
    }
    let e = exact3.min(l);
    'pos: for s in 0..=n - l {
        let w = &t[s..s + l];
        for i in (l - e..l).rev() {
            if w[i] & p[i] == 0 {
                continue 'pos;
            }
        }
        let mut mm = 0u32;
        for i in 0..l - e {
            if w[i] & p[i] == 0 {
                mm += 1;
                if mm > max_mm {
                    continue 'pos;
                }
            }
        }
        out.push((s, mm));
    }
    out
}

/// The oligo binds somewhere in `t`.
fn has_site(t: &[u8], p: &[u8], max_mm: u32, exact3: usize) -> bool {
    let (n, l) = (t.len(), p.len());
    if l == 0 || n < l {
        return false;
    }
    let e = exact3.min(l);
    'pos: for s in 0..=n - l {
        let w = &t[s..s + l];
        for i in (l - e..l).rev() {
            if w[i] & p[i] == 0 {
                continue 'pos;
            }
        }
        let mut mm = 0u32;
        for i in 0..l - e {
            if w[i] & p[i] == 0 {
                mm += 1;
                if mm > max_mm {
                    continue 'pos;
                }
            }
        }
        return true;
    }
    false
}

struct AssayMasks {
    f: Vec<u8>,
    r: Vec<u8>,
    p: Option<Vec<u8>>,
    window: usize,
}

/// (strand, start, end exclusive, forward mm, reverse mm) of every product
/// of one assay on one sequence (masks `t`, reverse complement `trc`).
fn products(t: &[u8], trc: &[u8], am: &AssayMasks, opt: &PcrOptions) -> (Vec<(char, usize, usize, u32, u32)>, bool) {
    let n = t.len();
    let (lf, lr) = (am.f.len(), am.r.len());
    let e3 = opt.three_prime_exact as usize;
    let mm = opt.max_mismatches;
    let fp = sites(t, &am.f, mm, e3);
    let fm = sites(trc, &am.f, mm, e3);
    let rp = sites(t, &am.r, mm, e3);
    let rm = sites(trc, &am.r, mm, e3);
    let mut out = Vec::new();
    let mut capped = false;
    // Forward primer on the sequence as written, reverse on its complement.
    let mut rq: Vec<(usize, u32)> = rm.iter().map(|&(p, m)| (n - p - lr, m)).collect();
    multi_merge_sort_stable_by(&mut rq, |a, b| a.0.cmp(&b.0));
    for &(p, fmm) in &fp {
        let lo = rq.partition_point(|x| x.0 < p);
        for &(q, rmm) in &rq[lo..] {
            let end = q + lr;
            if end - p > am.window {
                break;
            }
            if end < p + lf || fmm + rmm > opt.max_total_mismatches {
                continue;
            }
            if out.len() >= MAX_PRODUCTS {
                capped = true;
                break;
            }
            out.push(('+', p, end, fmm, rmm));
        }
    }
    // Reverse primer on the sequence as written, forward on its complement.
    let mut fq: Vec<(usize, u32)> = fm.iter().map(|&(p, m)| (n - p - lf, m)).collect();
    multi_merge_sort_stable_by(&mut fq, |a, b| a.0.cmp(&b.0));
    for &(p, rmm) in &rp {
        let lo = fq.partition_point(|x| x.0 < p);
        for &(q, fmm) in &fq[lo..] {
            let end = q + lf;
            if end - p > am.window {
                break;
            }
            if end < p + lr || fmm + rmm > opt.max_total_mismatches {
                continue;
            }
            if out.len() >= MAX_PRODUCTS {
                capped = true;
                break;
            }
            out.push(('-', p, end, fmm, rmm));
        }
    }
    multi_merge_sort_stable_by(&mut out, |a, b| (a.1, a.2).cmp(&(b.1, b.2)));
    (out, capped)
}

fn probe_mismatches(product: &[u8], probe: &[u8], max_mm: u32) -> Option<u32> {
    let t = template_masks(product);
    let trc = rc_masks(&t);
    sites(&t, probe, max_mm, 0).into_iter().chain(sites(&trc, probe, max_mm, 0)).map(|(_, m)| m).min()
}

fn in_range(a: &PcrAssay, size: u64) -> bool {
    a.min_size.map_or(true, |m| size >= u64::from(m)) && a.max_size.map_or(true, |m| size <= u64::from(m))
}

/// The panel against genomes (FASTA) and read samples (FASTQ): every
/// product, the reads with sites, and the presence/absence matrix.
pub fn in_silico_pcr(panel: &[PcrAssay], inputs: &[String], opt: &PcrOptions) -> Result<PcrResult, String> {
    if panel.is_empty() {
        return Err("The panel has no assay: give at least one primer pair".to_string());
    }
    if inputs.is_empty() {
        return Err("No genome or read file given: give FASTA files (genomes) or FASTQ files (reads)".to_string());
    }
    let mut assays = Vec::new();
    for (i, a) in panel.iter().enumerate() {
        assays.push(checked_assay(a, &format!("assay {}", i + 1))?);
    }
    let mut seen: HashMap<String, usize> = HashMap::new();
    for a in &assays {
        if seen.insert(a.target.clone(), 0).is_some() {
            return Err(format!("The panel has the target {} twice: give each assay its own name", a.target));
        }
    }
    let mut ins: Vec<&str> = Vec::new();
    for f in inputs {
        check_input(f)?;
        ins.push(f);
    }
    for o in [&opt.amplicons_output, &opt.matrix_output, &opt.reads_output].into_iter().flatten() {
        not_an_input(o, &ins)?;
    }
    let masks: Vec<AssayMasks> = assays
        .iter()
        .map(|a| AssayMasks {
            f: oligo_masks(&a.forward),
            r: oligo_masks(&a.reverse),
            p: a.probe.as_deref().map(oligo_masks),
            window: (opt.max_amplicon as usize).max(a.max_size.unwrap_or(0) as usize),
        })
        .collect();
    let nt = assays.len();
    let mut result = PcrResult {
        assays: assays.clone(),
        max_mismatches: opt.max_mismatches,
        three_prime_exact: opt.three_prime_exact,
        max_total_mismatches: opt.max_total_mismatches,
        max_probe_mismatches: opt.max_probe_mismatches,
        max_amplicon: opt.max_amplicon,
        min_reads: opt.min_reads,
        amplicons_output: opt.amplicons_output.clone(),
        matrix_output: opt.matrix_output.clone(),
        reads_output: opt.reads_output.clone(),
        ..Default::default()
    };
    for f in inputs {
        let mut rd = SeqReader::open(f)?;
        let mut rec = SeqRecord::default();
        if !rd.next(&mut rec)? {
            result.warnings.push(format!("{f} holds no sequence"));
            continue;
        }
        if rd.kind() == Some(SeqKind::Fastq) {
            let counts = pcr_reads(f, &mut rd, rec, &assays, &masks, opt)?;
            let mut row = PcrRow { name: file_label(f), kind: "reads".to_string(), file: f.clone(), present: Vec::new(), counts: Vec::new() };
            for c in counts {
                row.present.push(c.present);
                row.counts.push(c.with_site);
                result.read_counts.push(c);
            }
            result.rows.push(row);
            continue;
        }
        // Genomes: each record on its own, or the file as one genome.
        let label = file_label(f);
        let mut file_row = PcrRow { name: label.clone(), kind: "genome".to_string(), file: f.clone(), present: vec![false; nt], counts: vec![0; nt] };
        loop {
            let id = String::from_utf8_lossy(&rec.head).split_whitespace().next().unwrap_or("").to_string();
            let seq: Vec<u8> = rec.seq.iter().map(|b| b.to_ascii_uppercase()).collect();
            let genome = if opt.records_as_genomes { if id.is_empty() { label.clone() } else { id.clone() } } else { label.clone() };
            let t = template_masks(&seq);
            let trc = rc_masks(&t);
            let found: Vec<(Vec<Amplicon>, bool)> = (0..nt)
                .into_par_iter()
                .map(|ai| {
                    let (prods, capped) = products(&t, &trc, &masks[ai], opt);
                    let a = &assays[ai];
                    let amps = prods
                        .into_iter()
                        .map(|(strand, s, e, fmm, rmm)| {
                            let piece = if strand == '+' { seq[s..e].to_vec() } else { revcomp(&seq[s..e]) };
                            let size = (e - s) as u64;
                            let probe = masks[ai].p.as_ref().and_then(|p| probe_mismatches(&piece, p, opt.max_probe_mismatches));
                            let ok = in_range(a, size);
                            Amplicon {
                                genome: genome.clone(),
                                sequence_id: id.clone(),
                                target: a.target.clone(),
                                strand,
                                start: s as u64 + 1,
                                end: e as u64,
                                size,
                                forward_mismatches: fmm,
                                reverse_mismatches: rmm,
                                probe_mismatches: probe,
                                in_range: ok,
                                detected: ok && (a.probe.is_none() || probe.is_some()),
                                sequence: String::from_utf8_lossy(&piece).to_string(),
                            }
                        })
                        .collect();
                    (amps, capped)
                })
                .collect();
            let mut rec_row = PcrRow { name: genome.clone(), kind: "genome".to_string(), file: f.clone(), present: vec![false; nt], counts: vec![0; nt] };
            for (ai, (amps, capped)) in found.into_iter().enumerate() {
                if capped {
                    result.warnings.push(format!(
                        "{} on {id}: more than {MAX_PRODUCTS} products, the rest not listed: the primers bind almost everywhere (too many mismatches allowed?)",
                        assays[ai].target
                    ));
                }
                for a in amps {
                    if a.detected {
                        rec_row.counts[ai] += 1;
                        rec_row.present[ai] = true;
                    }
                    result.amplicons.push(a);
                }
            }
            if opt.records_as_genomes {
                result.rows.push(rec_row);
            } else {
                for ai in 0..nt {
                    file_row.counts[ai] += rec_row.counts[ai];
                    file_row.present[ai] |= rec_row.present[ai];
                }
            }
            if !rd.next(&mut rec)? {
                break;
            }
        }
        if !opt.records_as_genomes {
            result.rows.push(file_row);
        }
    }
    let mut names: Vec<String> = result.rows.iter().map(|r| r.name.clone()).collect();
    let before = names.clone();
    make_unique(&mut names, &mut result.warnings);
    if names != before {
        for (r, n) in result.rows.iter_mut().zip(&names) {
            r.name = n.clone();
        }
    }
    let mut cleanup = Cleanup(Vec::new());
    if let Some(p) = &opt.amplicons_output {
        cleanup.0.push(p.clone());
        write_text(p, &result.amplicons_tsv())?;
    }
    if let Some(p) = &opt.matrix_output {
        cleanup.0.push(p.clone());
        write_text(p, &result.matrix_tsv())?;
    }
    if let Some(p) = &opt.reads_output {
        cleanup.0.push(p.clone());
        write_text(p, &result.reads_tsv())?;
    }
    cleanup.disarm();
    Ok(result)
}

/// `read_pcr_panel` then `in_silico_pcr`.
pub fn in_silico_pcr_panel(panel: &str, inputs: &[String], opt: &PcrOptions) -> Result<PcrResult, String> {
    let ins: Vec<&str> = inputs.iter().map(String::as_str).collect();
    for o in [&opt.amplicons_output, &opt.matrix_output, &opt.reads_output].into_iter().flatten() {
        not_an_input(o, &[panel])?;
        not_an_input(o, &ins)?;
    }
    let assays = read_pcr_panel(panel)?;
    in_silico_pcr(&assays, inputs, opt)
}

/// Counts, per assay, the reads of a FASTQ file that carry its sites.
fn pcr_reads(
    file: &str,
    rd: &mut SeqReader,
    first: SeqRecord,
    assays: &[PcrAssay],
    masks: &[AssayMasks],
    opt: &PcrOptions,
) -> Result<Vec<ReadTargetCount>, String> {
    let nt = assays.len();
    // Per assay: forward, reverse, probe, amplicon, any.
    let mut tot = vec![[0u64; 5]; nt];
    let mut reads = 0u64;
    let e3 = opt.three_prime_exact as usize;
    let mut batch: Vec<Vec<u8>> = vec![first.seq];
    let mut rec = SeqRecord::default();
    loop {
        let more = loop {
            if batch.len() >= 20_000 {
                break true;
            }
            if !rd.next(&mut rec)? {
                break false;
            }
            batch.push(std::mem::take(&mut rec.seq));
        };
        reads += batch.len() as u64;
        let part: Vec<[u64; 5]> = batch
            .par_iter()
            .fold(
                || vec![[0u64; 5]; nt],
                |mut acc, s| {
                    let t = template_masks(s);
                    let trc = rc_masks(&t);
                    for (ai, am) in masks.iter().enumerate() {
                        let f = has_site(&t, &am.f, opt.max_mismatches, e3) || has_site(&trc, &am.f, opt.max_mismatches, e3);
                        let r = has_site(&t, &am.r, opt.max_mismatches, e3) || has_site(&trc, &am.r, opt.max_mismatches, e3);
                        let p = am.p.as_ref().is_some_and(|p| has_site(&t, p, opt.max_probe_mismatches, 0) || has_site(&trc, p, opt.max_probe_mismatches, 0));
                        let amp = f && r && !products(&t, &trc, am, opt).0.is_empty();
                        let flags = [f, r, p, amp, f || r || p];
                        for (k, &fl) in flags.iter().enumerate() {
                            if fl {
                                acc[ai][k] += 1;
                            }
                        }
                    }
                    acc
                },
            )
            .reduce(
                || vec![[0u64; 5]; nt],
                |mut a, b| {
                    for (x, y) in a.iter_mut().zip(&b) {
                        for k in 0..5 {
                            x[k] += y[k];
                        }
                    }
                    a
                },
            );
        for (x, y) in tot.iter_mut().zip(&part) {
            for k in 0..5 {
                x[k] += y[k];
            }
        }
        batch.clear();
        if !more {
            break;
        }
    }
    let sample = file_label(file);
    Ok(assays
        .iter()
        .zip(&tot)
        .map(|(a, c)| ReadTargetCount {
            sample: sample.clone(),
            target: a.target.clone(),
            reads,
            forward: c[0],
            reverse: c[1],
            probe: c[2],
            amplicon: c[3],
            with_site: c[4],
            present: c[4] >= opt.min_reads.max(1),
        })
        .collect())
}

impl PcrResult {
    pub fn amplicons_tsv(&self) -> String {
        let mut o = String::from(
            "genome\tsequence\ttarget\tstrand\tstart\tend\tsize\tforward_mismatches\treverse_mismatches\ttotal_mismatches\tprobe_mismatches\tin_range\tdetected\tamplicon\n",
        );
        for a in &self.amplicons {
            o.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                clean_cell(&a.genome),
                clean_cell(&a.sequence_id),
                clean_cell(&a.target),
                a.strand,
                a.start,
                a.end,
                a.size,
                a.forward_mismatches,
                a.reverse_mismatches,
                a.forward_mismatches + a.reverse_mismatches,
                a.probe_mismatches.map_or_else(|| "NA".to_string(), |m| m.to_string()),
                a.in_range,
                a.detected,
                a.sequence
            ));
        }
        o
    }

    /// 1 = detected, 0 = not.
    pub fn matrix_tsv(&self) -> String {
        let mut o = String::from("genome\ttype");
        for a in &self.assays {
            o.push('\t');
            o.push_str(&clean_cell(&a.target));
        }
        o.push('\n');
        for r in &self.rows {
            o.push_str(&format!("{}\t{}", clean_cell(&r.name), r.kind));
            for &p in &r.present {
                o.push_str(if p { "\t1" } else { "\t0" });
            }
            o.push('\n');
        }
        o
    }

    pub fn reads_tsv(&self) -> String {
        let mut o = String::from("sample\ttarget\treads\treads_forward\treads_reverse\treads_probe\treads_amplicon\treads_with_site\tpresent\n");
        for c in &self.read_counts {
            o.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                clean_cell(&c.sample),
                clean_cell(&c.target),
                c.reads,
                c.forward,
                c.reverse,
                c.probe,
                c.amplicon,
                c.with_site,
                c.present
            ));
        }
        o
    }

    pub fn json(&self) -> String {
        let assays: Vec<String> = self
            .assays
            .iter()
            .map(|a| {
                format!(
                    "{{\"target\":{},\"forward\":{},\"reverse\":{},\"probe\":{},\"min_size\":{},\"max_size\":{},\"source\":{}}}",
                    js(&a.target),
                    js(&a.forward),
                    js(&a.reverse),
                    jopt(&a.probe),
                    a.min_size.map_or_else(|| "null".to_string(), |v| v.to_string()),
                    a.max_size.map_or_else(|| "null".to_string(), |v| v.to_string()),
                    js(&a.source)
                )
            })
            .collect();
        let rows: Vec<String> = self
            .rows
            .iter()
            .map(|r| {
                format!(
                    "{{\"name\":{},\"type\":{},\"file\":{},\"present\":[{}],\"counts\":[{}]}}",
                    js(&r.name),
                    js(&r.kind),
                    js(&r.file),
                    r.present.iter().map(bool::to_string).collect::<Vec<_>>().join(","),
                    r.counts.iter().map(u64::to_string).collect::<Vec<_>>().join(",")
                )
            })
            .collect();
        let amps: Vec<String> = self
            .amplicons
            .iter()
            .map(|a| {
                format!(
                    "{{\"genome\":{},\"sequence_id\":{},\"target\":{},\"strand\":{},\"start\":{},\"end\":{},\"size\":{},\"forward_mismatches\":{},\"reverse_mismatches\":{},\"probe_mismatches\":{},\"in_range\":{},\"detected\":{},\"sequence\":{}}}",
                    js(&a.genome),
                    js(&a.sequence_id),
                    js(&a.target),
                    js(&a.strand.to_string()),
                    a.start,
                    a.end,
                    a.size,
                    a.forward_mismatches,
                    a.reverse_mismatches,
                    a.probe_mismatches.map_or_else(|| "null".to_string(), |m| m.to_string()),
                    a.in_range,
                    a.detected,
                    js(&a.sequence)
                )
            })
            .collect();
        let reads: Vec<String> = self
            .read_counts
            .iter()
            .map(|c| {
                format!(
                    "{{\"sample\":{},\"target\":{},\"reads\":{},\"forward\":{},\"reverse\":{},\"probe\":{},\"amplicon\":{},\"with_site\":{},\"present\":{}}}",
                    js(&c.sample),
                    js(&c.target),
                    c.reads,
                    c.forward,
                    c.reverse,
                    c.probe,
                    c.amplicon,
                    c.with_site,
                    c.present
                )
            })
            .collect();
        format!(
            "{{\"assays\":[{}],\"targets\":{},\"rows\":[{}],\"amplicons\":[{}],\"read_counts\":[{}],\"max_mismatches\":{},\"three_prime_exact\":{},\"max_total_mismatches\":{},\"max_probe_mismatches\":{},\"max_amplicon\":{},\"min_reads\":{},\"amplicons_output\":{},\"matrix_output\":{},\"reads_output\":{},\"warnings\":{},\"text\":{}}}",
            assays.join(","),
            jlist(&self.assays.iter().map(|a| a.target.clone()).collect::<Vec<_>>()),
            rows.join(","),
            amps.join(","),
            reads.join(","),
            self.max_mismatches,
            self.three_prime_exact,
            self.max_total_mismatches,
            self.max_probe_mismatches,
            self.max_amplicon,
            self.min_reads,
            jopt(&self.amplicons_output),
            jopt(&self.matrix_output),
            jopt(&self.reads_output),
            jlist(&self.warnings),
            js(&self.resumo())
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!(
            "In-silico PCR: {} assay(s) on {} genome(s) / read sample(s); up to {} mismatch(es) per primer ({} in both), none in the last {} 3' base(s); probe up to {}\n",
            self.assays.len(),
            self.rows.len(),
            self.max_mismatches,
            self.max_total_mismatches,
            self.three_prime_exact,
            self.max_probe_mismatches
        );
        let w = self.rows.iter().map(|r| r.name.chars().count()).max().unwrap_or(6).clamp(6, 40);
        o.push_str(&format!("  {:<w$}", "genome"));
        for a in &self.assays {
            o.push_str(&format!("  {}", a.target));
        }
        o.push('\n');
        for r in &self.rows {
            let name: String = r.name.chars().take(w).collect();
            o.push_str(&format!("  {name:<w$}"));
            for (i, a) in self.assays.iter().enumerate() {
                let cell = if r.present[i] {
                    if r.kind == "reads" {
                        format!("+ ({} reads)", r.counts[i])
                    } else if r.counts[i] > 1 {
                        format!("+ (x{})", r.counts[i])
                    } else {
                        "+".to_string()
                    }
                } else if r.kind == "reads" && r.counts[i] > 0 {
                    format!("- ({} reads)", r.counts[i])
                } else {
                    "-".to_string()
                };
                o.push_str(&format!("  {cell:<width$}", width = a.target.chars().count()));
            }
            o.push('\n');
        }
        let shown: Vec<&Amplicon> = self.amplicons.iter().filter(|a| a.detected).take(40).collect();
        if !shown.is_empty() {
            o.push_str("  products detected (first 40):\n");
            for a in shown {
                o.push_str(&format!(
                    "    {}  {}  {}:{}-{} ({}), {} bp, mismatches F{} R{}{}\n",
                    a.genome,
                    a.target,
                    a.sequence_id,
                    a.start,
                    a.end,
                    a.strand,
                    a.size,
                    a.forward_mismatches,
                    a.reverse_mismatches,
                    a.probe_mismatches.map_or_else(String::new, |m| format!(" P{m}"))
                ));
            }
        }
        let off = self.amplicons.iter().filter(|a| !a.detected).count();
        if off > 0 {
            o.push_str(&format!("  {off} other product(s): out of the size range, or the probe does not bind (see the amplicon table)\n"));
        }
        for p in [&self.amplicons_output, &self.matrix_output, &self.reads_output].into_iter().flatten() {
            o.push_str(&format!("  written: {p}\n"));
        }
        for w in &self.warnings {
            o.push_str(&format!("  warning: {w}\n"));
        }
        o
    }
}
