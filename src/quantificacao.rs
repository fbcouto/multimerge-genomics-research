// =====================================================================
//  RNA-SEQ QUANTIFICATION FROM FASTQ: PSEUDOALIGNMENT AND EM
//
//  In the manner of kallisto (Bray, Pimentel, Melsted & Pachter 2016,
//  Nat Biotechnol 34:525-527): no aligner, no genome, no BAM. Two tools:
//
//   - `build_transcript_index`  a transcriptome FASTA (cDNA) -> an index
//                               file: every canonical k-mer (both strands)
//                               -> the set of transcripts that contain it
//                               (its equivalence class), plus the
//                               transcript -> gene map;
//   - `quantify_samples`        FASTQ files (single- or paired-end, plain or
//                               gzip) -> estimated counts and TPM per
//                               transcript and per gene, and a gene x sample
//                               count matrix that
//                               `crate::expressao::differential_expression`
//                               reads as it is.
//
//  PSEUDOALIGNMENT. A read is compatible with the transcripts that contain
//  ALL of its k-mers that are in the index: the intersection of the classes
//  of its k-mers. As in kallisto, k-mers absent from the index (sequencing
//  errors, introns, other species) are ignored; k-mers with a base other than
//  A, C, G, T are skipped; an empty intersection leaves the read unassigned.
//  The two mates of a pair are one fragment: the classes of both are
//  intersected (a mate with no k-mer in the index adds nothing).
//  Unlike kallisto, every k-mer of the read is looked up: kallisto jumps over
//  the k-mers of a de Bruijn contig once it has seen one of them (faster,
//  and nearly always the same answer); here the answer is the exact
//  intersection, and needs no graph.
//
//  EM. Reads are counted per equivalence class; the abundances are the
//  maximum-likelihood estimate found by the EM algorithm of kallisto 0.46
//  (same start, same weights 1 / effective length, same stopping rule: no
//  transcript above 0.01 counts changing by more than 1% after at least 50
//  rounds, then the abundances below 1e-8 are set to zero and one last round
//  runs). Effective length = length - mean fragment length + 1, the mean
//  taken over fragments no longer than the transcript (kallisto's
//  truncation), fragments up to 1,000 bp. Single-end reads need the fragment
//  length (mean and sd: a truncated normal); for pairs it is estimated from
//  the first 10,000 pairs that fall on a single transcript, by placing both
//  mates on it (kallisto only uses pairs whose mates fall on the same graph
//  contig, which favours short fragments in genes with many isoforms; here
//  the whole transcript is used).
//
//  TPM_t = 10^6 * (count_t / efflen_t) / sum_u (count_u / efflen_u). A gene
//  gets the sum of its transcripts' counts and TPMs. The gene count matrix
//  holds the rounded estimated counts (as tximport + DESeqDataSetFromMatrix
//  do).
//
//  Also as kallisto: poly-A tails (10 or more A's at a transcript's end) are
//  left out of the k-mers (the length stays the original one). kallisto's
//  default rule for SINGLE-END reads -- drop the transcripts on which a
//  fragment of the mean length, starting at the read, would run past the
//  end -- is an option (`single_overhang = false`), off by default: on
//  simulated data it improved isoform-level ranking a little but lost reads
//  near transcript ends and worsened gene-level counts (gene MARD 0.06 -> 0.09
//  on human chr21 transcripts); with it on, the results match kallisto's.
//
//  Not done: sequence-specific bias correction, bootstraps, strand-specific
//  filtering, pseudobam.
//
//  THE INDEX IN MEMORY. The k-mers are kept as their mixed values (a
//  bijection of the 2k-bit canonical k-mer, so no information is lost) in one
//  sorted array, with the class of each in a parallel array and a directory
//  of buckets by the top bits: 8 + 4 bytes per k-mer plus about half a byte
//  of directory; with the 2-bit transcript sequences and the classes, the
//  file measured 13.7 MB per million distinct k-mers on human chr21
//  transcripts (24.2 MB; kallisto 0.46's index of the same FASTA: 38.0 MB).
//  A lookup reads one directory entry and scans a bucket of about 8 keys;
//  the k-mers of a read are looked up in phases (all buckets, then all keys,
//  then all classes) so the processor overlaps their cache misses.
//
//  BUILDING. The (k-mer, transcript) pairs are produced in passes, each pass
//  taking the k-mers whose mixed value falls in one slice of the 64-bit
//  range, so the memory needed is set by `memory_mb`, not by the size of the
//  transcriptome. Each pass is sorted (Multimerge_Stable), grouped into
//  classes and streamed to temporary files; the final file is assembled by
//  streaming: the whole index is never in memory while it is built.
//
//  THE INDEX ON DISK (little-endian, format version 1):
//
//     "MMGENTXI"                       8 bytes, the signature
//     u32 format version, u32 k
//     u64 transcripts T, genes G, classes E (the first T are the single
//         transcripts), k-mers N, class data length D, bases B, non-ACGT
//         bases X
//     u32 bucket bits b, u32 flags (bit 0: poly-A tails clipped)
//     text FASTA path, text gene map source, text producer
//                                      (text = u32 byte length + UTF-8)
//     T texts: transcript names
//     T x u32 lengths, T x u32 gene of each transcript
//     G texts gene ids, G texts gene names
//     (E + 1) x u64 class offsets, D x u32 class transcripts (ascending)
//     (2^b + 1) x u32 directory (first key of each bucket)
//     N x u64 keys (mixed canonical k-mers, ascending)
//     N x u32 class of each key
//     (T + 1) x u64 sequence offsets (bases), ceil(B / 4) bytes of 2-bit
//         bases (A=0 C=1 G=2 T=3, first base in the low bits),
//     X x u64 positions of the non-ACGT bases
//     "MMGENEND", u64 file length
//
//  It is written to `<output>.part` and renamed at the end: a failed build
//  never leaves a file that looks complete.
// =====================================================================

use crate::externo::{abre_leitura, Escritor};
use crate::multimerge_core::{multi_merge_sort_stable, multi_merge_sort_stable_by};
use crate::nomes::json_texto;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::time::Instant;

/// kallisto's default k.
pub const DEFAULT_K: usize = 31;
/// Fragments longer than this are not modelled (kallisto's MAX_FRAG_LEN).
pub const MAX_FRAG_LEN: usize = 1000;
/// The version of the index file this code writes and reads.
pub const INDEX_FORMAT_VERSION: u32 = 1;

const ASSINATURA: &[u8; 8] = b"MMGENTXI";
const ASSINATURA_FIM: &[u8; 8] = b"MMGENEND";
const PRODUTOR: &str = "MMGEN quantificacao (pseudoalignment index, format 1)";
/// Reads (or pairs) per batch handed to the threads.
const LOTE: usize = 16_384;
/// The smallest positive double: kallisto's TOLERANCE in the EM.
const TOLERANCIA: f64 = 5e-324;
/// Memory for the pairs of one building pass when `memory_mb` is not given.
const ORCAMENTO_PADRAO_MB: usize = 1024;

// ---------------------------------------------------------------------
//  K-MERS
// ---------------------------------------------------------------------

const fn tabela_bases() -> [u8; 256] {
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
}

/// ASCII -> 2-bit code; 4 for anything that is not A, C, G, T or U.
static BASE: [u8; 256] = tabela_bases();

/// The finalizer of splitmix64: a bijection of u64 (xor-shifts and odd
/// multiplications are invertible), so the mixed value identifies the k-mer.
#[inline]
pub fn mistura(x: u64) -> u64 {
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Calls `f(pos, fwd, rc)` for every k-mer of `seq` made only of A, C, G, T
/// (any case; U reads as T): `pos` is its first base, `fwd` its 2-bit code as
/// read, `rc` the code of its reverse complement. The canonical k-mer is
/// `min(fwd, rc)`; with an odd k a k-mer is never its own reverse complement.
#[inline]
pub fn kmers(seq: &[u8], k: usize, mut f: impl FnMut(usize, u64, u64)) {
    let mascara = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let desloca = 2 * (k - 1);
    let (mut fwd, mut rc, mut validas) = (0u64, 0u64, 0usize);
    for (i, &c) in seq.iter().enumerate() {
        let b = BASE[c as usize];
        if b > 3 {
            validas = 0;
            continue;
        }
        fwd = ((fwd << 2) | b as u64) & mascara;
        rc = (rc >> 2) | (((3 - b) as u64) << desloca);
        validas += 1;
        if validas >= k {
            f(i + 1 - k, fwd, rc);
        }
    }
}

/// The 2-bit code of a k-mer written as text (A, C, G, T only).
pub fn codifica(kmer: &[u8]) -> Option<u64> {
    let mut x = 0u64;
    for &c in kmer {
        let b = BASE[c as usize];
        if b > 3 {
            return None;
        }
        x = (x << 2) | b as u64;
    }
    Some(x)
}

/// The k-mer of a 2-bit code, as text.
pub fn decodifica_kmer(mut x: u64, k: usize) -> String {
    let mut v = vec![b'A'; k];
    for i in (0..k).rev() {
        v[i] = b"ACGT"[(x & 3) as usize];
        x >>= 2;
    }
    String::from_utf8(v).unwrap_or_default()
}

// ---------------------------------------------------------------------
//  A FAST HASH (for hash maps keyed by classes)
// ---------------------------------------------------------------------

const MULT: u64 = 0x9E37_79B9_7F4A_7C15;

#[derive(Default)]
struct Mistura(u64);

impl Hasher for Mistura {
    #[inline]
    fn finish(&self) -> u64 {
        mistura(self.0)
    }
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        for c in bytes.chunks(8) {
            let mut b = [0u8; 8];
            b[..c.len()].copy_from_slice(c);
            self.0 = (self.0.rotate_left(5) ^ u64::from_le_bytes(b)).wrapping_mul(MULT);
        }
    }
    #[inline]
    fn write_u32(&mut self, x: u32) {
        self.write_u64(x as u64)
    }
    #[inline]
    fn write_u64(&mut self, x: u64) {
        self.0 = (self.0.rotate_left(5) ^ x).wrapping_mul(MULT);
    }
    #[inline]
    fn write_usize(&mut self, x: usize) {
        self.write_u64(x as u64)
    }
}

type Rapido = BuildHasherDefault<Mistura>;

fn hash_conjunto(c: &[u32]) -> u64 {
    let mut h = c.len() as u64;
    for &x in c {
        h = (h.rotate_left(5) ^ x as u64).wrapping_mul(MULT);
    }
    mistura(h)
}

/// The transcripts in both sorted lists, into `saida`.
#[inline]
fn intersecao(a: &[u32], b: &[u32], saida: &mut Vec<u32>) {
    saida.clear();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                saida.push(a[i]);
                i += 1;
                j += 1;
            }
        }
    }
}

// ---------------------------------------------------------------------
//  TRANSCRIPT SEQUENCES, 2 BITS PER BASE
// ---------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq)]
struct Sequencias {
    /// 4 bases per byte, base i in bits 2*(i%4).
    bases: Vec<u8>,
    /// Start of each transcript in bases; T + 1 entries when complete.
    inicio: Vec<u64>,
    /// Global positions of the bases that are not A, C, G, T (stored as A).
    ns: Vec<u64>,
    n: u64,
}

impl Sequencias {
    fn novo() -> Sequencias {
        Sequencias { bases: Vec::new(), inicio: vec![0], ns: Vec::new(), n: 0 }
    }

    fn empurra(&mut self, linha: &[u8]) {
        for &c in linha {
            if c == b' ' || c == b'\t' || c == b'\r' {
                continue;
            }
            let mut b = BASE[c as usize];
            if b > 3 {
                self.ns.push(self.n);
                b = 0;
            }
            let i = (self.n % 4) as u32;
            if i == 0 {
                self.bases.push(0);
            }
            if let Some(u) = self.bases.last_mut() {
                *u |= b << (2 * i);
            }
            self.n += 1;
        }
    }

    fn fecha(&mut self) {
        self.inicio.push(self.n);
    }

    fn comprimento(&self, t: usize) -> u64 {
        self.inicio[t + 1] - self.inicio[t]
    }

    #[inline]
    fn base(&self, p: u64) -> u8 {
        (self.bases[(p / 4) as usize] >> (2 * (p % 4))) & 3
    }

    /// The transcript as text (A, C, G, T and N).
    fn decodifica(&self, t: usize, saida: &mut Vec<u8>) {
        self.decodifica_trecho(t, 0, self.comprimento(t), saida);
    }

    /// Bases `ini..fim` of a transcript as text.
    fn decodifica_trecho(&self, t: usize, ini: u64, fim: u64, saida: &mut Vec<u8>) {
        saida.clear();
        let (a, b) = (self.inicio[t] + ini, self.inicio[t] + fim.min(self.comprimento(t)));
        if a >= b {
            return;
        }
        saida.reserve((b - a) as usize);
        for p in a..b {
            saida.push(b"ACGT"[self.base(p) as usize]);
        }
        let i0 = self.ns.partition_point(|&x| x < a);
        for &p in &self.ns[i0..] {
            if p >= b {
                break;
            }
            saida[(p - a) as usize] = b'N';
        }
    }
}

// ---------------------------------------------------------------------
//  READING THE TRANSCRIPTOME
// ---------------------------------------------------------------------

struct Transcritoma {
    nomes: Vec<String>,
    /// (gene id, gene name) read from each header; empty when absent.
    gene_cab: Vec<(String, String)>,
    seqs: Sequencias,
    renomeados: usize,
}

/// The last "(...)" group with no space inside: the gene symbol of a RefSeq
/// header ("Homo sapiens alpha-2-macroglobulin (A2M), transcript variant 1,
/// mRNA").
fn simbolo_entre_parenteses(desc: &str) -> Option<&str> {
    let mut melhor = None;
    let mut resto = desc;
    while let Some(i) = resto.find('(') {
        let depois = &resto[i + 1..];
        let Some(j) = depois.find(')') else { break };
        let dentro = &depois[..j];
        if !dentro.is_empty() && !dentro.contains(char::is_whitespace) && !dentro.contains('(') {
            melhor = Some(dentro);
        }
        resto = &depois[j + 1..];
    }
    melhor
}

/// The transcript name and, when the header says it, the gene id and name:
///
///  - GENCODE: `>ENST...|ENSG...|OTTHUMG|OTTHUMT|NAME-201|NAME|1657|biotype|`
///    -> the first field, gene id the second, gene name the sixth;
///  - Ensembl: `>ENST... cdna chromosome:... gene:ENSG... gene_symbol:NAME`;
///  - RefSeq: `>NM_000014.6 Homo sapiens alpha-2-macroglobulin (A2M), ...`
///    -> gene id and name = the symbol in parentheses;
///  - anything else: the first word, no gene.
pub fn parse_header(cabecalho: &str) -> (String, String, String) {
    let h = cabecalho.trim();
    let primeiro = h.split_whitespace().next().unwrap_or("");
    if primeiro.contains('|') {
        let f: Vec<&str> = primeiro.split('|').collect();
        let pega = |i: usize| f.get(i).map(|s| s.trim()).filter(|s| !s.is_empty() && *s != "-").unwrap_or("").to_string();
        let id = pega(1);
        let nome = pega(5);
        return (f[0].to_string(), id, nome);
    }
    let (mut gid, mut gnome) = (String::new(), String::new());
    for tok in h.split_whitespace().skip(1) {
        if let Some(v) = tok.strip_prefix("gene:") {
            gid = v.to_string();
        } else if let Some(v) = tok.strip_prefix("gene_symbol:") {
            gnome = v.to_string();
        }
    }
    if gid.is_empty() {
        let refseq = ["NM_", "NR_", "XM_", "XR_"].iter().any(|p| primeiro.starts_with(p));
        if refseq {
            let desc = h[primeiro.len()..].trim();
            if let Some(s) = simbolo_entre_parenteses(desc) {
                gid = s.to_string();
                gnome = s.to_string();
            }
        }
    }
    (primeiro.to_string(), gid, gnome)
}

fn le_transcritoma(caminho: &str, unicos: bool) -> Result<Transcritoma, String> {
    let mut r = abre_leitura(caminho)?;
    let mut tr = Transcritoma { nomes: Vec::new(), gene_cab: Vec::new(), seqs: Sequencias::novo(), renomeados: 0 };
    let mut vistos: HashMap<String, usize> = HashMap::new();
    let mut linha: Vec<u8> = Vec::new();
    let mut n = 0usize;
    let mut aberto = false;
    loop {
        linha.clear();
        let lidos = r.read_until(b'\n', &mut linha).map_err(|e| {
            format!("Error reading {caminho} after line {n}: {e}. The file may be truncated or compressed in a format that is not read here (plain text, .gz and .lz4 are).")
        })?;
        if lidos == 0 {
            break;
        }
        n += 1;
        while matches!(linha.last(), Some(b'\n' | b'\r')) {
            linha.pop();
        }
        if linha.first() == Some(&b'>') {
            if aberto {
                tr.seqs.fecha();
            }
            let cab = String::from_utf8_lossy(&linha[1..]).to_string();
            let (mut nome, gid, gnome) = parse_header(&cab);
            if nome.is_empty() {
                return Err(format!("{caminho}, line {n}: a header with no name ('>' and nothing after it); every transcript needs a name"));
            }
            if let Some(&l0) = vistos.get(&nome) {
                if !unicos {
                    return Err(format!(
                        "{caminho}, line {n}: the transcript name \"{nome}\" was already used on line {l0}. Each transcript needs a distinct name (the first word of its header); a file may have been concatenated twice. To keep the repeats as NAME_2, NAME_3..., use make_unique=True."
                    ));
                }
                let mut i = 2;
                loop {
                    let c = format!("{nome}_{i}");
                    if !vistos.contains_key(&c) {
                        nome = c;
                        break;
                    }
                    i += 1;
                }
                tr.renomeados += 1;
            }
            vistos.insert(nome.clone(), n);
            tr.nomes.push(nome);
            tr.gene_cab.push((gid, gnome));
            aberto = true;
        } else if linha.iter().all(|c| c.is_ascii_whitespace()) {
            continue;
        } else if !aberto {
            if linha[0] == b'@' {
                return Err(format!(
                    "{caminho} looks like a FASTQ file of reads (line {n} starts with '@'). The index is built from the transcript sequences: a FASTA file (cDNA) whose records start with '>', such as Ensembl's *.cdna.all.fa.gz or GENCODE's *.transcripts.fa.gz."
                ));
            }
            return Err(format!(
                "{caminho}, line {n}: sequence before the first header. A FASTA file starts with a '>' line naming the transcript (is the file compressed without a .gz name?)"
            ));
        } else {
            tr.seqs.empurra(&linha);
        }
    }
    if aberto {
        tr.seqs.fecha();
    }
    if tr.nomes.is_empty() {
        return Err(format!("{caminho} has no transcript (no line starting with '>'): it is empty or not a FASTA file"));
    }
    for t in 0..tr.nomes.len() {
        if tr.seqs.comprimento(t) > u32::MAX as u64 {
            return Err(format!("{caminho}: the transcript {} has more than 4,294,967,295 bases; this is a genome, not a transcriptome", tr.nomes[t]));
        }
    }
    Ok(tr)
}

// ---------------------------------------------------------------------
//  TRANSCRIPT -> GENE
// ---------------------------------------------------------------------

/// `ENST00000456328.2` -> `ENST00000456328` (only a final `.digits`).
fn sem_versao(s: &str) -> &str {
    if let Some(i) = s.rfind('.') {
        let v = &s[i + 1..];
        if !v.is_empty() && v.bytes().all(|c| c.is_ascii_digit()) && i > 0 {
            return &s[..i];
        }
    }
    s
}

#[derive(Default)]
struct MapaGenes {
    exato: HashMap<String, (String, String)>,
    sem_versao: HashMap<String, (String, String)>,
    tipo: &'static str,
}

impl MapaGenes {
    fn insere(&mut self, tid: &str, gid: &str, gnome: &str) {
        let (tid, gid) = (tid.trim(), gid.trim());
        if tid.is_empty() || gid.is_empty() {
            return;
        }
        let v = (gid.to_string(), gnome.trim().to_string());
        self.sem_versao.entry(sem_versao(tid).to_string()).or_insert_with(|| v.clone());
        self.exato.entry(tid.to_string()).or_insert(v);
    }

    fn procura(&self, nome: &str) -> Option<&(String, String)> {
        self.exato.get(nome).or_else(|| self.sem_versao.get(sem_versao(nome)))
    }
}

fn desescapa(s: &str) -> String {
    if !s.contains('%') {
        return s.to_string();
    }
    let hex = |c: u8| -> Option<u8> { (c as char).to_digit(16).map(|d| d as u8) };
    let b = s.as_bytes();
    let mut o = Vec::with_capacity(b.len());
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

fn atributos_gtf(s: &str) -> Vec<(&str, &str)> {
    s.split(';')
        .filter_map(|p| {
            let p = p.trim();
            if p.is_empty() {
                return None;
            }
            let (k, v) = p.split_once(char::is_whitespace)?;
            Some((k.trim(), v.trim().trim_matches('"')))
        })
        .collect()
}

fn atributos_gff(s: &str) -> Vec<(&str, &str)> {
    s.split(';').filter_map(|p| p.trim().split_once('=')).map(|(k, v)| (k.trim(), v.trim())).collect()
}

/// Sub-transcript features of a GFF3: not transcripts.
const PARTES_DE_TRANSCRITO: [&str; 12] = [
    "exon",
    "CDS",
    "five_prime_UTR",
    "three_prime_UTR",
    "UTR",
    "start_codon",
    "stop_codon",
    "intron",
    "polyA_site",
    "polyA_signal_sequence",
    "stop_codon_redefined_as_selenocysteine",
    "Selenocysteine",
];

struct FeatGff {
    pai: String,
    tid: String,
    gid: String,
    nome: String,
    versao: String,
}

/// A gene map: a GTF (transcript_id -> gene_id, gene_name), a GFF3
/// (a transcript feature -> its top-level ancestor, the gene) or a table
/// (transcript, gene[, gene name]; tab- or comma-separated).
fn le_mapa_genes(caminho: &str) -> Result<MapaGenes, String> {
    let mut r = abre_leitura(caminho)?;
    let mut m = MapaGenes::default();
    let mut linha = String::new();
    let mut n = 0usize;
    let mut tipo: Option<&'static str> = None;
    let mut feats: HashMap<String, FeatGff> = HashMap::new();
    let mut primeira_tabela = true;
    loop {
        linha.clear();
        let lidos = r.read_line(&mut linha).map_err(|e| format!("Error reading {caminho} after line {n}: {e}"))?;
        if lidos == 0 {
            break;
        }
        n += 1;
        let t = linha.trim_end_matches(['\n', '\r']);
        if t.trim().is_empty() || t.starts_with('#') {
            continue;
        }
        if tipo.is_none() {
            let c: Vec<&str> = t.split('\t').collect();
            tipo = Some(if c.len() == 9 && (c[8].contains("transcript_id \"") || c[8].contains("gene_id \"")) {
                "GTF"
            } else if c.len() == 9 && c[8].contains('=') {
                "GFF3"
            } else {
                "table"
            });
        }
        match tipo {
            Some("GTF") => {
                let c: Vec<&str> = t.split('\t').collect();
                if c.len() < 9 {
                    continue;
                }
                let at = atributos_gtf(c[8]);
                let pega = |k: &str| at.iter().find(|(a, _)| *a == k).map(|(_, v)| *v).unwrap_or("");
                let tid = pega("transcript_id");
                if tid.is_empty() {
                    continue;
                }
                let gid = pega("gene_id");
                let mut gnome = pega("gene_name");
                if gnome.is_empty() {
                    gnome = pega("gene_symbol");
                }
                if gnome.is_empty() {
                    gnome = pega("gene");
                }
                m.insere(tid, gid, gnome);
                let tv = pega("transcript_version");
                if !tv.is_empty() && sem_versao(tid) == tid {
                    m.insere(&format!("{tid}.{tv}"), gid, gnome);
                }
            }
            Some("GFF3") => {
                let c: Vec<&str> = t.split('\t').collect();
                if c.len() < 9 {
                    continue;
                }
                let tipo_f = c[2];
                let at = atributos_gff(c[8]);
                let pega = |k: &str| at.iter().find(|(a, _)| *a == k).map(|(_, v)| desescapa(v)).unwrap_or_default();
                let id = pega("ID");
                if id.is_empty() {
                    continue;
                }
                let pai = pega("Parent").split(',').next().unwrap_or("").to_string();
                if pai.is_empty() {
                    if !tipo_f.to_ascii_lowercase().contains("gene") {
                        continue;
                    }
                } else if PARTES_DE_TRANSCRITO.contains(&tipo_f) {
                    continue;
                }
                let mut nome = pega("gene_name");
                if nome.is_empty() {
                    nome = pega("Name");
                }
                if nome.is_empty() {
                    nome = pega("gene");
                }
                let f = FeatGff { pai, tid: pega("transcript_id"), gid: pega("gene_id"), nome, versao: pega("version") };
                feats.entry(id).or_insert(f);
            }
            _ => {
                let s = if t.contains('\t') { '\t' } else { ',' };
                let c: Vec<&str> = t.split(s).map(|x| x.trim().trim_matches('"')).collect();
                if primeira_tabela {
                    primeira_tabela = false;
                    let h = c[0].to_ascii_lowercase();
                    if ["transcript", "transcript_id", "target_id", "txname", "tx", "tx_id", "transcript id"].contains(&h.as_str()) {
                        continue;
                    }
                }
                if c.len() < 2 {
                    return Err(format!(
                        "{caminho}, line {n}: a gene map needs two columns (transcript, gene) and optionally a third (gene name), separated by tabs or commas; it can also be a GTF or GFF3 annotation"
                    ));
                }
                m.insere(c[0], c[1], c.get(2).copied().unwrap_or(""));
            }
        }
    }
    if tipo == Some("GFF3") {
        let tira = |s: &str, prefixos: &[&str]| -> String {
            for p in prefixos {
                if let Some(x) = s.strip_prefix(p) {
                    return x.to_string();
                }
            }
            s.to_string()
        };
        let mut ids: Vec<&String> = feats.keys().collect();
        // a ordem so importa para "o primeiro vence" entre chaves repetidas
        multi_merge_sort_stable_by(&mut ids, |a, b| a.cmp(b));
        for id in ids {
            let f = &feats[id];
            if f.pai.is_empty() {
                continue;
            }
            // sobe ate o ancestral sem pai: o gene
            let mut topo_id = f.pai.clone();
            let mut topo: Option<&FeatGff> = feats.get(&topo_id);
            let mut passos = 0;
            while let Some(p) = topo {
                if p.pai.is_empty() || passos > 16 {
                    break;
                }
                match feats.get(&p.pai) {
                    Some(q) => {
                        topo_id = p.pai.clone();
                        topo = Some(q);
                    }
                    None => {
                        topo_id = p.pai.clone();
                        topo = None;
                    }
                }
                passos += 1;
            }
            let gid = match topo {
                Some(g) if !g.gid.is_empty() => g.gid.clone(),
                _ => tira(&topo_id, &["gene:", "gene-"]),
            };
            let gnome = topo.map(|g| g.nome.clone()).filter(|s| !s.is_empty()).unwrap_or_else(|| gid.clone());
            let chave_id = tira(id, &["transcript:", "rna-"]);
            for chave in [f.tid.as_str(), chave_id.as_str()] {
                if chave.is_empty() {
                    continue;
                }
                m.insere(chave, &gid, &gnome);
                if !f.versao.is_empty() && sem_versao(chave) == chave {
                    m.insere(&format!("{chave}.{}", f.versao), &gid, &gnome);
                }
            }
        }
    }
    m.tipo = tipo.unwrap_or("table");
    if m.exato.is_empty() {
        return Err(format!(
            "{caminho} gave no transcript -> gene pair. Use a GTF (transcript_id and gene_id attributes), a GFF3 (transcripts with a Parent gene) or a table with two columns: transcript, gene."
        ));
    }
    Ok(m)
}

struct Genes {
    de: Vec<u32>,
    ids: Vec<String>,
    nomes: Vec<String>,
    fonte: String,
    sem_gene: usize,
}

fn atribui_genes(nomes_tr: &[String], cab: &[(String, String)], mapa: Option<&str>, avisos: &mut Vec<String>) -> Result<Genes, String> {
    let externo = match mapa {
        Some(p) if !p.trim().is_empty() => Some(le_mapa_genes(p)?),
        _ => None,
    };
    let mut g = Genes { de: Vec::with_capacity(nomes_tr.len()), ids: Vec::new(), nomes: Vec::new(), fonte: String::new(), sem_gene: 0 };
    let mut posicao: HashMap<String, u32> = HashMap::new();
    let mut tem_cab = 0usize;
    for (t, nome) in nomes_tr.iter().enumerate() {
        let (gid, gnome) = match &externo {
            Some(m) => m.procura(nome).cloned().unwrap_or_default(),
            None => cab.get(t).cloned().unwrap_or_default(),
        };
        let (gid, gnome) = if gid.is_empty() {
            g.sem_gene += 1;
            (nome.clone(), String::new())
        } else {
            tem_cab += 1;
            (gid, gnome)
        };
        let i = match posicao.get(&gid) {
            Some(&i) => i,
            None => {
                let i = g.ids.len() as u32;
                posicao.insert(gid.clone(), i);
                g.ids.push(gid);
                g.nomes.push(String::new());
                i
            }
        };
        if g.nomes[i as usize].is_empty() && !gnome.is_empty() {
            g.nomes[i as usize] = gnome;
        }
        g.de.push(i);
    }
    for (i, nm) in g.nomes.iter_mut().enumerate() {
        if nm.is_empty() {
            *nm = g.ids[i].clone();
        }
    }
    match (&externo, mapa) {
        (Some(m), Some(p)) => {
            g.fonte = format!("{p} ({})", m.tipo);
            if g.sem_gene == nomes_tr.len() {
                return Err(format!(
                    "None of the {} transcript names of the FASTA is in the gene map {p} (for example \"{}\"): the map is for another annotation or uses other names. Use the GTF/GFF3 of the same release as the FASTA.",
                    milhar(nomes_tr.len() as u64),
                    nomes_tr[0]
                ));
            }
            if g.sem_gene > 0 {
                avisos.push(format!(
                    "{} transcript(s) are not in the gene map {p}: each of them is counted as its own gene",
                    milhar(g.sem_gene as u64)
                ));
            }
        }
        _ => {
            if tem_cab > 0 {
                g.fonte = "FASTA headers".into();
                if g.sem_gene > 0 {
                    avisos.push(format!(
                        "{} transcript header(s) name no gene: each of those transcripts is counted as its own gene",
                        milhar(g.sem_gene as u64)
                    ));
                }
            } else {
                g.fonte = "none (each transcript is its own gene)".into();
                avisos.push(
                    "The FASTA headers name no gene (Ensembl `gene:`/`gene_symbol:`, GENCODE `|` fields or a RefSeq `(SYMBOL)`): each transcript is its own gene. Give a gene map (GTF, GFF3 or a transcript-gene table) for gene-level results.".into(),
                );
            }
        }
    }
    Ok(g)
}

// ---------------------------------------------------------------------
//  BINARY INPUT AND OUTPUT
// ---------------------------------------------------------------------

struct Bin<W: Write> {
    w: W,
    n: u64,
    buf: Vec<u8>,
}

impl<W: Write> Bin<W> {
    fn novo(w: W) -> Bin<W> {
        Bin { w, n: 0, buf: Vec::with_capacity(1 << 16) }
    }
    fn bytes(&mut self, b: &[u8]) -> std::io::Result<()> {
        self.w.write_all(b)?;
        self.n += b.len() as u64;
        Ok(())
    }
    fn u32(&mut self, x: u32) -> std::io::Result<()> {
        self.bytes(&x.to_le_bytes())
    }
    fn u64(&mut self, x: u64) -> std::io::Result<()> {
        self.bytes(&x.to_le_bytes())
    }
    fn texto(&mut self, s: &str) -> std::io::Result<()> {
        self.u32(s.len() as u32)?;
        self.bytes(s.as_bytes())
    }
    fn u32s(&mut self, v: &[u32]) -> std::io::Result<()> {
        for c in v.chunks(1 << 14) {
            self.buf.clear();
            for x in c {
                self.buf.extend_from_slice(&x.to_le_bytes());
            }
            self.w.write_all(&self.buf)?;
            self.n += self.buf.len() as u64;
        }
        Ok(())
    }
    fn u64s(&mut self, v: &[u64]) -> std::io::Result<()> {
        for c in v.chunks(1 << 13) {
            self.buf.clear();
            for x in c {
                self.buf.extend_from_slice(&x.to_le_bytes());
            }
            self.w.write_all(&self.buf)?;
            self.n += self.buf.len() as u64;
        }
        Ok(())
    }
}

struct LeBin<R: Read> {
    r: R,
    pos: u64,
    total: u64,
    caminho: String,
}

impl<R: Read> LeBin<R> {
    fn truncado(&self, parte: &str) -> String {
        format!(
            "{} is truncated: it ends at byte {} inside the {parte}. The index was not written completely (a full disk, an interrupted copy or download); build it again with build_transcript_index.",
            self.caminho, self.total
        )
    }
    fn exato(&mut self, buf: &mut [u8], parte: &str) -> Result<(), String> {
        if self.pos + buf.len() as u64 > self.total {
            return Err(self.truncado(parte));
        }
        self.r.read_exact(buf).map_err(|e| {
            if e.kind() == std::io::ErrorKind::UnexpectedEof {
                self.truncado(parte)
            } else {
                format!("Error reading {}: {e}", self.caminho)
            }
        })?;
        self.pos += buf.len() as u64;
        Ok(())
    }
    fn u32(&mut self, parte: &str) -> Result<u32, String> {
        let mut b = [0u8; 4];
        self.exato(&mut b, parte)?;
        Ok(u32::from_le_bytes(b))
    }
    fn u64(&mut self, parte: &str) -> Result<u64, String> {
        let mut b = [0u8; 8];
        self.exato(&mut b, parte)?;
        Ok(u64::from_le_bytes(b))
    }
    fn cabe(&self, n: u64, tamanho: u64, parte: &str) -> Result<usize, String> {
        match n.checked_mul(tamanho) {
            Some(b) if self.pos + b <= self.total => Ok(n as usize),
            _ => Err(self.truncado(parte)),
        }
    }
    fn texto(&mut self, parte: &str) -> Result<String, String> {
        let n = self.u32(parte)? as u64;
        let n = self.cabe(n, 1, parte)?;
        let mut b = vec![0u8; n];
        self.exato(&mut b, parte)?;
        String::from_utf8(b).map_err(|_| format!("{}: the {parte} is not valid text; the file is damaged. Build the index again.", self.caminho))
    }
    fn u32s(&mut self, n: u64, parte: &str) -> Result<Vec<u32>, String> {
        let n = self.cabe(n, 4, parte)?;
        let mut v = Vec::with_capacity(n);
        let mut buf = vec![0u8; 1 << 20];
        let mut faltam = n;
        while faltam > 0 {
            let m = faltam.min(buf.len() / 4);
            self.exato(&mut buf[..4 * m], parte)?;
            v.extend(buf[..4 * m].chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])));
            faltam -= m;
        }
        Ok(v)
    }
    fn u64s(&mut self, n: u64, parte: &str) -> Result<Vec<u64>, String> {
        let n = self.cabe(n, 8, parte)?;
        let mut v = Vec::with_capacity(n);
        let mut buf = vec![0u8; 1 << 20];
        let mut faltam = n;
        while faltam > 0 {
            let m = faltam.min(buf.len() / 8);
            self.exato(&mut buf[..8 * m], parte)?;
            v.extend(buf[..8 * m].chunks_exact(8).map(|c| u64::from_le_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]])));
            faltam -= m;
        }
        Ok(v)
    }
    fn bytes(&mut self, n: u64, parte: &str) -> Result<Vec<u8>, String> {
        let n = self.cabe(n, 1, parte)?;
        let mut v = vec![0u8; n];
        self.exato(&mut v, parte)?;
        Ok(v)
    }
}

/// The fixed part of the index file.
struct Cabecalho {
    versao: u32,
    k: usize,
    t: u64,
    g: u64,
    e: u64,
    n: u64,
    d: u64,
    bases: u64,
    x: u64,
    bits: u32,
    fasta: String,
    fonte_genes: String,
    produtor: String,
}

fn le_cabecalho<R: Read>(r: &mut LeBin<R>) -> Result<Cabecalho, String> {
    let mut sig = [0u8; 8];
    if r.total < 8 {
        return Err(format!(
            "{} is not an MMGEN transcript index (it is {} bytes long). Build one from a transcriptome FASTA with build_transcript_index.",
            r.caminho, r.total
        ));
    }
    r.exato(&mut sig, "signature")?;
    if &sig != ASSINATURA {
        return Err(format!(
            "{} is not an MMGEN transcript index (it does not start with the MMGENTXI signature; a kallisto or salmon index is not read here). Build one from a transcriptome FASTA with build_transcript_index.",
            r.caminho
        ));
    }
    let versao = r.u32("header")?;
    if versao != INDEX_FORMAT_VERSION {
        return Err(format!(
            "{} is a transcript index of format version {versao}; this MMGEN reads version {INDEX_FORMAT_VERSION}. Build the index again with this version of MMGEN (build_transcript_index).",
            r.caminho
        ));
    }
    let k = r.u32("header")? as usize;
    let t = r.u64("header")?;
    let g = r.u64("header")?;
    let e = r.u64("header")?;
    let n = r.u64("header")?;
    let d = r.u64("header")?;
    let bases = r.u64("header")?;
    let x = r.u64("header")?;
    let bits = r.u32("header")?;
    let _flags = r.u32("header")?;
    if !(1..=32).contains(&k) || k % 2 == 0 || bits == 0 || bits > 32 || e < t {
        return Err(format!("{}: the header of the index is damaged (k = {k}, {bits} bucket bits). Build the index again.", r.caminho));
    }
    let fasta = r.texto("header")?;
    let fonte_genes = r.texto("header")?;
    let produtor = r.texto("header")?;
    Ok(Cabecalho { versao, k, t, g, e, n, d, bases, x, bits, fasta, fonte_genes, produtor })
}

// ---------------------------------------------------------------------
//  THE INDEX
// ---------------------------------------------------------------------

/// A transcript index loaded in memory.
pub struct TranscriptIndex {
    pub k: usize,
    pub path: String,
    pub fasta: String,
    pub gene_map: String,
    pub bytes: u64,
    pub names: Vec<String>,
    pub lengths: Vec<u32>,
    /// Gene (index into `gene_ids`) of each transcript.
    pub gene_of: Vec<u32>,
    pub gene_ids: Vec<String>,
    pub gene_names: Vec<String>,
    ec_inicio: Vec<u64>,
    ec_dados: Vec<u32>,
    desloca: u32,
    diretorio: Vec<u32>,
    chaves: Vec<u64>,
    classes: Vec<u32>,
    seqs: Sequencias,
}

/// What a read (or pair) pseudoaligns to.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Classe {
    /// No complete k-mer: shorter than k, or N-rich.
    Curta,
    /// k-mers, none of them in the index.
    SemKmer,
    /// k-mers in the index whose classes have no transcript in common.
    Conflito,
    /// A class of the index (an id below the number of transcripts is that
    /// transcript alone).
    Ec(u32),
    /// A set of transcripts that is not a class of the index (an
    /// intersection of classes), ascending.
    Novo(Vec<u32>),
    /// Single-end with `single_overhang` off: there are compatible
    /// transcripts, but on none of them would a fragment of the mean length,
    /// starting at the read, fit.
    Fora,
}

#[derive(Default)]
struct Rascunho {
    ecs: Vec<u32>,
    a: Vec<u32>,
    b: Vec<u32>,
    /// First k-mer of each mate found in the index: (offset in the read,
    /// code as read).
    primeiro: [Option<(usize, u64)>; 2],
    validos: bool,
    mates: u8,
    buf: Vec<u8>,
    /// Per k-mer of the mate being read: mixed canonical value, offset, code
    /// as read, bucket bounds, position of the key (u32::MAX: absent).
    hs: Vec<u64>,
    pos: Vec<u32>,
    fwd: Vec<u64>,
    faixa: Vec<(u32, u32)>,
    onde: Vec<u32>,
}

impl TranscriptIndex {
    pub fn transcripts(&self) -> usize {
        self.names.len()
    }
    pub fn classes(&self) -> usize {
        self.ec_inicio.len() - 1
    }
    pub fn kmers(&self) -> usize {
        self.chaves.len()
    }

    /// The transcripts of a class, ascending.
    #[inline]
    pub fn ec(&self, id: u32) -> &[u32] {
        let i = id as usize;
        &self.ec_dados[self.ec_inicio[i] as usize..self.ec_inicio[i + 1] as usize]
    }

    /// The class of a canonical k-mer (`min(fwd, rc)`), if it is in the index.
    #[inline]
    pub fn lookup(&self, canonico: u64) -> Option<u32> {
        let h = mistura(canonico);
        let b = (h >> self.desloca) as usize;
        let (ini, fim) = (self.diretorio[b] as usize, self.diretorio[b + 1] as usize);
        for (i, &x) in self.chaves[ini..fim].iter().enumerate() {
            if x == h {
                return Some(self.classes[ini + i]);
            }
            if x > h {
                return None;
            }
        }
        None
    }

    /// The sequence of a transcript (A, C, G, T, N).
    pub fn sequence(&self, t: usize) -> String {
        let mut v = Vec::new();
        self.seqs.decodifica(t, &mut v);
        String::from_utf8(v).unwrap_or_default()
    }

    /// Loads an index written by `build_transcript_index`.
    pub fn load(caminho: &str) -> Result<TranscriptIndex, String> {
        let f = File::open(caminho).map_err(|e| format!("Error opening the transcript index {caminho}: {e}"))?;
        let total = f.metadata().map_err(|e| format!("Error reading {caminho}: {e}"))?.len();
        let mut r = LeBin { r: BufReader::with_capacity(1 << 20, f), pos: 0, total, caminho: caminho.to_string() };
        let c = le_cabecalho(&mut r)?;
        let mut names = Vec::with_capacity(r.cabe(c.t, 4, "transcript names")?);
        for _ in 0..c.t {
            names.push(r.texto("transcript names")?);
        }
        let lengths = r.u32s(c.t, "transcript lengths")?;
        let gene_of = r.u32s(c.t, "transcript genes")?;
        let mut gene_ids = Vec::with_capacity(r.cabe(c.g, 4, "gene ids")?);
        for _ in 0..c.g {
            gene_ids.push(r.texto("gene ids")?);
        }
        let mut gene_names = Vec::with_capacity(c.g as usize);
        for _ in 0..c.g {
            gene_names.push(r.texto("gene names")?);
        }
        let ec_inicio = r.u64s(c.e + 1, "class offsets")?;
        let ec_dados = r.u32s(c.d, "classes")?;
        let diretorio = r.u32s((1u64 << c.bits) + 1, "k-mer directory")?;
        let chaves = r.u64s(c.n, "k-mers")?;
        let classes = r.u32s(c.n, "k-mer classes")?;
        let inicio = r.u64s(c.t + 1, "sequence offsets")?;
        let bases = r.bytes(c.bases.div_ceil(4), "sequences")?;
        let ns = r.u64s(c.x, "sequences")?;
        let mut fim = [0u8; 8];
        r.exato(&mut fim, "end signature")?;
        let comprimento = r.u64("end signature")?;
        let danificado = |o: &str| format!("{caminho} is damaged ({o}). Build the index again with build_transcript_index.");
        if &fim != ASSINATURA_FIM || comprimento != total || r.pos != total {
            return Err(danificado("the end of the file is not where the header says"));
        }
        // conferencias baratas: um arquivo corrompido vira erro claro, nao panico
        let (t, e, d, n) = (c.t as usize, c.e as usize, c.d, c.n as usize);
        if ec_inicio[0] != 0 || ec_inicio[e] != d || ec_inicio.windows(2).any(|w| w[0] > w[1]) {
            return Err(danificado("class offsets"));
        }
        if ec_dados.iter().any(|&x| x as usize >= t) {
            return Err(danificado("a class names a transcript that does not exist"));
        }
        if gene_of.iter().any(|&g| g as u64 >= c.g) {
            return Err(danificado("gene of a transcript"));
        }
        if diretorio[0] != 0 || *diretorio.last().unwrap_or(&0) as usize != n || diretorio.windows(2).any(|w| w[0] > w[1]) {
            return Err(danificado("k-mer directory"));
        }
        if classes.iter().any(|&x| x as usize >= e) {
            return Err(danificado("class of a k-mer"));
        }
        if inicio[0] != 0 || inicio[t] != c.bases || inicio.windows(2).any(|w| w[0] > w[1]) || ns.iter().any(|&p| p >= c.bases) {
            return Err(danificado("sequence offsets"));
        }
        for i in 0..t {
            if lengths[i] as u64 != inicio[i + 1] - inicio[i] {
                return Err(danificado("transcript lengths"));
            }
        }
        let _ = c.versao;
        let _ = c.produtor;
        Ok(TranscriptIndex {
            k: c.k,
            path: caminho.to_string(),
            fasta: c.fasta,
            gene_map: c.fonte_genes,
            bytes: total,
            names,
            lengths,
            gene_of,
            gene_ids,
            gene_names,
            ec_inicio,
            ec_dados,
            desloca: 64 - c.bits,
            diretorio,
            chaves,
            classes,
            seqs: Sequencias { bases, n: c.bases, inicio, ns },
        })
    }

    /// The classes of the k-mers of one mate. The lookups go in phases --
    /// every hash, then every bucket, then every key, then every class -- so
    /// that the memory reads of different k-mers do not wait for each other
    /// (the processor overlaps independent cache misses); the result is the
    /// same as looking them up one by one.
    fn coleta(&self, s: &[u8], mate: usize, r: &mut Rascunho) {
        r.hs.clear();
        r.pos.clear();
        r.fwd.clear();
        {
            let (hs, pos, fwd) = (&mut r.hs, &mut r.pos, &mut r.fwd);
            kmers(s, self.k, |p, f, rc| {
                hs.push(mistura(f.min(rc)));
                pos.push(p as u32);
                fwd.push(f);
            });
        }
        if r.hs.is_empty() {
            return;
        }
        r.validos = true;
        r.faixa.clear();
        for &h in &r.hs {
            let b = (h >> self.desloca) as usize;
            r.faixa.push((self.diretorio[b], self.diretorio[b + 1]));
        }
        r.onde.clear();
        for (i, &h) in r.hs.iter().enumerate() {
            let (ini, fim) = r.faixa[i];
            let mut achado = u32::MAX;
            for j in ini..fim {
                let x = self.chaves[j as usize];
                if x >= h {
                    if x == h {
                        achado = j;
                    }
                    break;
                }
            }
            r.onde.push(achado);
        }
        let mut ultimo = u32::MAX;
        let mut achou = false;
        for i in 0..r.onde.len() {
            let j = r.onde[i];
            if j == u32::MAX {
                continue;
            }
            let ec = self.classes[j as usize];
            if !achou {
                achou = true;
                r.primeiro[mate] = Some((r.pos[i] as usize, r.fwd[i]));
            }
            if ec != ultimo {
                ultimo = ec;
                if !r.ecs.contains(&ec) {
                    r.ecs.push(ec);
                }
            }
        }
        if achou {
            r.mates += 1;
        }
    }

    fn pseudoalinha_com(&self, s1: &[u8], s2: Option<&[u8]>, r: &mut Rascunho) -> Classe {
        r.ecs.clear();
        r.primeiro = [None, None];
        r.validos = false;
        r.mates = 0;
        self.coleta(s1, 0, r);
        if let Some(s2) = s2 {
            self.coleta(s2, 1, r);
        }
        if r.ecs.is_empty() {
            return if r.validos { Classe::SemKmer } else { Classe::Curta };
        }
        if r.ecs.len() == 1 {
            return Classe::Ec(r.ecs[0]);
        }
        let Rascunho { ecs, a, b, .. } = r;
        // comeca pela menor classe; empate: a primeira vista (deterministico)
        let mut imin = 0;
        for i in 1..ecs.len() {
            if self.ec(ecs[i]).len() < self.ec(ecs[imin]).len() {
                imin = i;
            }
        }
        let emin = ecs[imin];
        a.clear();
        a.extend_from_slice(self.ec(emin));
        for (i, &e) in ecs.iter().enumerate() {
            if i == imin {
                continue;
            }
            intersecao(a, self.ec(e), b);
            std::mem::swap(a, b);
            if a.is_empty() {
                return Classe::Conflito;
            }
        }
        if a.len() == 1 {
            Classe::Ec(a[0])
        } else if a.len() == self.ec(emin).len() {
            Classe::Ec(emin)
        } else {
            Classe::Novo(a.clone())
        }
    }

    /// Pseudoaligns one read, or a pair (`r2`).
    pub fn pseudoalign(&self, r1: &[u8], r2: Option<&[u8]>) -> Classe {
        let mut r = Rascunho::default();
        self.pseudoalinha_com(r1, r2, &mut r)
    }

    /// The fragment length of a pair, when it pseudoaligns to a single
    /// transcript and both mates can be placed on it (what the fragment-length
    /// estimate uses).
    pub fn pair_fragment_length(&self, r1: &[u8], r2: &[u8]) -> Option<u32> {
        let mut r = Rascunho::default();
        let c = self.pseudoalinha_com(r1, Some(r2), &mut r);
        match (c, r.primeiro) {
            (Classe::Ec(t), [Some(h1), Some(h2)]) if (t as usize) < self.transcripts() => self.fragmento(t, h1, h2, &mut r.buf),
            _ => None,
        }
    }

    /// The transcripts of a pseudoalignment (empty when unassigned).
    pub fn transcripts_of(&self, c: &Classe) -> Vec<u32> {
        match c {
            Classe::Ec(e) => self.ec(*e).to_vec(),
            Classe::Novo(v) => v.clone(),
            _ => Vec::new(),
        }
    }

    /// The fragment length of a pair on transcript `t`, from the first k-mer
    /// of each mate found in the index, as kallisto's `mapPair` does on its
    /// contigs: the mates must lie on opposite strands; |start of the
    /// forward mate - end of the reverse mate|. `None` outside (0, 1000).
    fn fragmento(&self, t: u32, h1: (usize, u64), h2: (usize, u64), buf: &mut Vec<u8>) -> Option<u32> {
        self.seqs.decodifica(t as usize, buf);
        let k = self.k as i64;
        let mut p: [Option<(i64, bool)>; 2] = [None, None];
        let hs = [h1, h2];
        kmers(buf, self.k, |q, f, rc| {
            for m in 0..2 {
                if p[m].is_none() {
                    let (o, x) = hs[m];
                    if f == x {
                        p[m] = Some((q as i64 - o as i64, true));
                    } else if rc == x {
                        p[m] = Some((q as i64 + k + o as i64, false));
                    }
                }
            }
        });
        let ((a, da), (b, db)) = (p[0]?, p[1]?);
        if da == db {
            return None;
        }
        let tl = (a - b).abs();
        if tl > 0 && (tl as usize) < MAX_FRAG_LEN {
            Some(tl as u32)
        } else {
            None
        }
    }

    /// kallisto's single-end rule (on unless --single-overhang): of the
    /// transcripts `ts`, keep those on which a fragment of length `fl`
    /// starting at the read fits. The read is placed by its first k-mer found
    /// in the index, `h` (offset, code as read), as kallisto's findPosition
    /// does: on the forward strand at q, its 1-based start is q - o + 1 and
    /// the transcript is dropped when start + fl > length; on the reverse
    /// strand, the 1-based position of its first base is q + k + o and it is
    /// dropped when that - fl < 0. Only an occurrence near an end can fail,
    /// so only the two end windows of each transcript are read (fl + k bases
    /// each), not the whole transcript. (For a k-mer that occurs more than
    /// once in a transcript, kallisto uses one of the occurrences; here the
    /// transcript is dropped when any occurrence fails.)
    fn cabe_fragmento(&self, ts: &[u32], h: (usize, u64), fl: i64, buf: &mut Vec<u8>, saida: &mut Vec<u32>) {
        saida.clear();
        let k = self.k as i64;
        let o = h.0 as i64;
        for &t in ts {
            let l = self.lengths[t as usize] as i64;
            let mut cai = false;
            // fita direta: q >= l - fl + o
            let a = (l - fl + o).max(0);
            if a + k <= l {
                self.seqs.decodifica_trecho(t as usize, a as u64, l as u64, buf);
                kmers(buf, self.k, |_, f, _| cai |= f == h.1);
            }
            // fita reversa: q < fl - k - o
            let b = fl - k - o;
            if !cai && b > 0 {
                self.seqs.decodifica_trecho(t as usize, 0, (b + k - 1).min(l) as u64, buf);
                kmers(buf, self.k, |_, _, rc| cai |= rc == h.1);
            }
            if !cai {
                saida.push(t);
            }
        }
    }

    /// The multi-transcript classes by their transcript sets.
    fn mapa_ecs(&self) -> HashMap<&[u32], u32, Rapido> {
        let t = self.transcripts();
        let mut m: HashMap<&[u32], u32, Rapido> = HashMap::with_capacity_and_hasher(self.classes() - t, Rapido::default());
        for e in t..self.classes() {
            m.entry(self.ec(e as u32)).or_insert(e as u32);
        }
        m
    }
}

// ---------------------------------------------------------------------
//  BUILDING THE INDEX
// ---------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct IndexOptions {
    /// k-mer length: odd, 7 to 31. Default 31 (kallisto's).
    pub k: usize,
    /// A GTF, a GFF3 or a two/three-column table (transcript, gene[, gene
    /// name]). `None`: the gene of each transcript from the FASTA headers
    /// (Ensembl, GENCODE, RefSeq); a transcript without one is its own gene.
    pub gene_map: Option<String>,
    /// Rename repeated transcript names as NAME_2, NAME_3... instead of
    /// stopping (kallisto's --make-unique).
    pub make_unique: bool,
    /// Memory for the (k-mer, transcript) pairs of one building pass, in MB.
    /// `None`: 1,024. A smaller value means more passes, not a failure.
    pub memory_mb: Option<usize>,
    /// A transcript that ends in 10 or more A's loses that poly-A tail before
    /// its k-mers are taken (as kallisto does): poly-A k-mers would make every
    /// poly-A read compatible with every tailed transcript. The length used
    /// for the effective length stays the original one, as in kallisto.
    /// Default true.
    pub clip_poly_a: bool,
}

impl Default for IndexOptions {
    fn default() -> Self {
        IndexOptions { k: DEFAULT_K, gene_map: None, make_unique: false, memory_mb: None, clip_poly_a: true }
    }
}

#[derive(Clone, Debug, Default)]
pub struct IndexSummary {
    pub fasta: String,
    pub index: String,
    pub format_version: u32,
    pub k: usize,
    pub transcripts: usize,
    pub genes: usize,
    /// "FASTA headers", the gene map file and its kind, or "none".
    pub gene_map: String,
    pub transcripts_without_gene: usize,
    pub bases: u64,
    pub non_acgt_bases: u64,
    /// Transcripts shorter than k: they have no k-mer, and get no read.
    pub short_transcripts: usize,
    /// Distinct canonical k-mers.
    pub kmers: u64,
    /// Equivalence classes, the single transcripts included.
    pub equivalence_classes: u64,
    pub renamed: usize,
    /// Transcripts whose poly-A tail was left out of the k-mers.
    pub poly_a_clipped: usize,
    pub bytes: u64,
    pub passes: usize,
    pub seconds: f64,
    pub warnings: Vec<String>,
}

impl IndexSummary {
    pub fn bytes_per_million_kmers(&self) -> f64 {
        if self.kmers == 0 {
            0.0
        } else {
            self.bytes as f64 / self.kmers as f64 * 1e6
        }
    }

    pub fn json(&self) -> String {
        format!(
            "{{\"fasta\":{},\"index\":{},\"format_version\":{},\"k\":{},\"transcripts\":{},\"genes\":{},\"gene_map\":{},\"transcripts_without_gene\":{},\"bases\":{},\"non_acgt_bases\":{},\"short_transcripts\":{},\"kmers\":{},\"equivalence_classes\":{},\"renamed\":{},\"poly_a_clipped\":{},\"bytes\":{},\"bytes_per_million_kmers\":{},\"passes\":{},\"seconds\":{},\"warnings\":{},\"text\":{}}}",
            json_texto(&self.fasta),
            json_texto(&self.index),
            self.format_version,
            self.k,
            self.transcripts,
            self.genes,
            json_texto(&self.gene_map),
            self.transcripts_without_gene,
            self.bases,
            self.non_acgt_bases,
            self.short_transcripts,
            self.kmers,
            self.equivalence_classes,
            self.renamed,
            self.poly_a_clipped,
            self.bytes,
            jnum((self.bytes_per_million_kmers() * 10.0).round() / 10.0),
            self.passes,
            jnum((self.seconds * 1000.0).round() / 1000.0),
            jtextos(&self.warnings),
            json_texto(&self.resumo())
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!(
            "Transcript index {} (format {}, k = {}): {} transcripts of {} genes (gene map: {}), {} bases; {} distinct k-mers in {} equivalence classes; {} ({:.1} MB per million k-mers)",
            self.index,
            self.format_version,
            self.k,
            milhar(self.transcripts as u64),
            milhar(self.genes as u64),
            self.gene_map,
            milhar(self.bases),
            milhar(self.kmers),
            milhar(self.equivalence_classes),
            mb(self.bytes),
            self.bytes_per_million_kmers() / 1e6
        );
        if self.passes > 0 {
            o.push_str(&format!("; built from {} in {:.1} s, {} pass(es)", self.fasta, self.seconds, self.passes));
        }
        o.push('\n');
        if self.short_transcripts > 0 {
            o.push_str(&format!("  {} transcript(s) shorter than k = {}: no read can be assigned to them\n", self.short_transcripts, self.k));
        }
        for w in &self.warnings {
            o.push_str(&format!("  warning: {w}\n"));
        }
        o
    }
}

/// Removes temporary files (and a half-written index) when building stops.
struct Temporarios(Vec<String>);

impl Drop for Temporarios {
    fn drop(&mut self) {
        for p in &self.0 {
            let _ = std::fs::remove_file(p);
        }
    }
}

/// The multi-transcript classes found so far, by their sets.
struct DicionarioEc {
    mapa: HashMap<u64, u32, Rapido>,
    proximo: Vec<u32>,
    inicio: Vec<u64>,
    dados: Vec<u32>,
    t: u32,
}

impl DicionarioEc {
    fn novo(t: u32) -> DicionarioEc {
        DicionarioEc { mapa: HashMap::default(), proximo: Vec::new(), inicio: vec![0], dados: Vec::new(), t }
    }

    fn fatia(&self, id: u32) -> &[u32] {
        let i = (id - self.t) as usize;
        &self.dados[self.inicio[i] as usize..self.inicio[i + 1] as usize]
    }

    fn id_de(&mut self, conj: &[u32]) -> u32 {
        if conj.len() == 1 {
            return conj[0];
        }
        let h = hash_conjunto(conj);
        let novo = self.t + self.proximo.len() as u32;
        match self.mapa.get(&h) {
            Some(&primeiro) => {
                let mut id = primeiro;
                loop {
                    if self.fatia(id) == conj {
                        return id;
                    }
                    let px = self.proximo[(id - self.t) as usize];
                    if px == u32::MAX {
                        break;
                    }
                    id = px;
                }
                self.proximo[(id - self.t) as usize] = novo;
            }
            None => {
                self.mapa.insert(h, novo);
            }
        }
        self.proximo.push(u32::MAX);
        self.dados.extend_from_slice(conj);
        self.inicio.push(self.dados.len() as u64);
        novo
    }
}

fn bits_para(n: u64) -> u32 {
    let alvo = (n / 8).max(1);
    let mut b = 1;
    while (1u64 << b) < alvo && b < 32 {
        b += 1;
    }
    b
}

/// Which pass a mixed k-mer belongs to: slices of the 64-bit range, so the
/// passes come out in ascending order of the keys.
#[inline]
fn particao(h: u64, passes: u64) -> u64 {
    ((h >> 32) * passes) >> 32
}

/// Where the k-mers of a transcript end: before a final run of A's when the
/// last 10 bases are all A (kallisto's poly-A clipping).
fn fim_sem_poli_a(seq: &[u8], cortar: bool) -> usize {
    let n = seq.len();
    if !cortar || n < 10 || !seq[n - 10..].iter().all(|&c| c == b'A' || c == b'a') {
        return n;
    }
    let mut j = n;
    while j > 0 && (seq[j - 1] == b'A' || seq[j - 1] == b'a') {
        j -= 1;
    }
    j
}

fn erro_escrita(caminho: &str, e: std::io::Error) -> String {
    format!("Error writing {caminho}: {e} (is the disk full?)")
}

/// Builds the pseudoalignment index of a transcriptome FASTA (plain, .gz or
/// .lz4): every canonical k-mer -> the set of transcripts that contain it.
pub fn build_transcript_index(fasta: &str, output: &str, opt: &IndexOptions) -> Result<IndexSummary, String> {
    let inicio = Instant::now();
    let k = opt.k;
    if !(7..=31).contains(&k) || k % 2 == 0 {
        return Err(format!(
            "k = {k}: k must be an odd number from 7 to 31 (kallisto's default, 31, is right for reads of 50 bp or more; use a smaller odd k such as 21 or 23 only for shorter reads). An odd k-mer is never its own reverse complement."
        ));
    }
    recusa_mesmo_arquivo(fasta, output)?;
    if let Some(m) = &opt.gene_map {
        recusa_mesmo_arquivo(m, output)?;
    }
    if output.trim().is_empty() {
        return Err("No output path for the index: give a file name such as transcripts.idx".into());
    }
    let orcamento = (opt.memory_mb.unwrap_or(ORCAMENTO_PADRAO_MB).max(1) as u64) << 20;
    let mut avisos = Vec::new();
    let tr = le_transcritoma(fasta, opt.make_unique)?;
    let genes = atribui_genes(&tr.nomes, &tr.gene_cab, opt.gene_map.as_deref(), &mut avisos)?;
    let nt = tr.nomes.len();
    let seqs = &tr.seqs;

    let temp_chaves = format!("{output}.kmers.tmp");
    let temp_classes = format!("{output}.classes.tmp");
    let parcial = format!("{output}.part");
    let _guarda = Temporarios(vec![temp_chaves.clone(), temp_classes.clone(), parcial.clone()]);

    // --- as passadas ---------------------------------------------------
    let ocorrencias: u64 = (0..nt).map(|t| seqs.comprimento(t).saturating_sub(k as u64 - 1)).sum();
    let por_passe = (orcamento / 40).max(1 << 16);
    let passes = ocorrencias.div_ceil(por_passe).clamp(1, 4096);
    // blocos de ~1 Mb de transcritos para as threads
    let mut blocos: Vec<(usize, usize)> = Vec::new();
    let (mut a, mut acum) = (0usize, 0u64);
    for t in 0..nt {
        acum += seqs.comprimento(t);
        if acum >= 1 << 20 {
            blocos.push((a, t + 1));
            a = t + 1;
            acum = 0;
        }
    }
    if a < nt {
        blocos.push((a, nt));
    }
    let abre = |p: &str| -> Result<Bin<BufWriter<File>>, String> {
        Ok(Bin::novo(BufWriter::with_capacity(1 << 20, File::create(p).map_err(|e| format!("Error creating {p}: {e}"))?)))
    };
    let mut w_chaves = abre(&temp_chaves)?;
    let mut w_classes = abre(&temp_classes)?;
    let mut dic = DicionarioEc::novo(nt as u32);
    let mut n_kmers = 0u64;
    let mut conj: Vec<u32> = Vec::new();
    for p in 0..passes {
        let partes: Vec<Vec<u128>> = blocos
            .par_iter()
            .map(|&(a, b)| {
                let mut v = Vec::new();
                let mut buf = Vec::new();
                for t in a..b {
                    seqs.decodifica(t, &mut buf);
                    let fim = fim_sem_poli_a(&buf, opt.clip_poly_a);
                    kmers(&buf[..fim], k, |_, f, rc| {
                        let h = mistura(f.min(rc));
                        if particao(h, passes) == p {
                            v.push(((h as u128) << 32) | t as u128);
                        }
                    });
                }
                v
            })
            .collect();
        let total: usize = partes.iter().map(Vec::len).sum();
        let mut pares: Vec<u128> = Vec::with_capacity(total);
        for v in partes {
            pares.extend_from_slice(&v);
        }
        multi_merge_sort_stable(&mut pares);
        let mut i = 0;
        while i < pares.len() {
            let h = (pares[i] >> 32) as u64;
            conj.clear();
            let mut j = i;
            while j < pares.len() && (pares[j] >> 32) as u64 == h {
                let t = pares[j] as u32;
                if conj.last() != Some(&t) {
                    conj.push(t);
                }
                j += 1;
            }
            let id = dic.id_de(&conj);
            w_chaves.u64(h).map_err(|e| erro_escrita(&temp_chaves, e))?;
            w_classes.u32(id).map_err(|e| erro_escrita(&temp_classes, e))?;
            n_kmers += 1;
            i = j;
        }
    }
    w_chaves.w.flush().map_err(|e| erro_escrita(&temp_chaves, e))?;
    w_classes.w.flush().map_err(|e| erro_escrita(&temp_classes, e))?;
    drop(w_chaves);
    drop(w_classes);
    if n_kmers >= u32::MAX as u64 {
        return Err(format!(
            "{fasta} has {} distinct k-mers: more than this index format holds (4,294,967,294). Split the transcriptome, or use a protein-coding subset.",
            milhar(n_kmers)
        ));
    }
    if n_kmers == 0 {
        return Err(format!(
            "{fasta} has no k-mer of length {k} made only of A, C, G, T: every transcript is shorter than k or made of N. Check the file, or use a smaller k."
        ));
    }
    if dic.dados.len() as u64 + nt as u64 >= u32::MAX as u64 || dic.proximo.len() as u64 + nt as u64 >= u32::MAX as u64 {
        return Err(format!("{fasta}: too many equivalence classes for this index format"));
    }

    // --- o diretorio, a partir das chaves gravadas -----------------------
    let bits = bits_para(n_kmers);
    let mut diretorio = vec![0u32; (1usize << bits) + 1];
    {
        let f = File::open(&temp_chaves).map_err(|e| format!("Error opening {temp_chaves}: {e}"))?;
        let mut r = LeBin { r: BufReader::with_capacity(1 << 20, f), pos: 0, total: n_kmers * 8, caminho: temp_chaves.clone() };
        let mut faltam = n_kmers;
        while faltam > 0 {
            let m = faltam.min(1 << 17);
            for h in r.u64s(m, "k-mers")? {
                diretorio[(h >> (64 - bits)) as usize + 1] += 1;
            }
            faltam -= m;
        }
        for i in 1..diretorio.len() {
            diretorio[i] += diretorio[i - 1];
        }
    }

    // --- o arquivo final -------------------------------------------------
    let n_ecs = nt as u64 + dic.proximo.len() as u64;
    let d = nt as u64 + dic.dados.len() as u64;
    let curtos = (0..nt).filter(|&t| seqs.comprimento(t) < k as u64).count();
    let cortados = if opt.clip_poly_a {
        let mut buf = Vec::new();
        (0..nt)
            .filter(|&t| {
                seqs.decodifica(t, &mut buf);
                fim_sem_poli_a(&buf, true) < buf.len()
            })
            .count()
    } else {
        0
    };
    {
        let f = File::create(&parcial).map_err(|e| format!("Error creating {parcial}: {e}"))?;
        let mut w = Bin::novo(BufWriter::with_capacity(1 << 20, f));
        let er = |e: std::io::Error| erro_escrita(output, e);
        w.bytes(ASSINATURA).map_err(er)?;
        w.u32(INDEX_FORMAT_VERSION).map_err(er)?;
        w.u32(k as u32).map_err(er)?;
        for x in [nt as u64, genes.ids.len() as u64, n_ecs, n_kmers, d, seqs.n, seqs.ns.len() as u64] {
            w.u64(x).map_err(er)?;
        }
        w.u32(bits).map_err(er)?;
        // flags: bit 0 = poly-A tails clipped
        w.u32(opt.clip_poly_a as u32).map_err(er)?;
        w.texto(fasta).map_err(er)?;
        w.texto(&genes.fonte).map_err(er)?;
        w.texto(PRODUTOR).map_err(er)?;
        for nome in &tr.nomes {
            w.texto(nome).map_err(er)?;
        }
        let comprimentos: Vec<u32> = (0..nt).map(|t| seqs.comprimento(t) as u32).collect();
        w.u32s(&comprimentos).map_err(er)?;
        w.u32s(&genes.de).map_err(er)?;
        for g in &genes.ids {
            w.texto(g).map_err(er)?;
        }
        for g in &genes.nomes {
            w.texto(g).map_err(er)?;
        }
        let mut ofs: Vec<u64> = (0..=nt as u64).collect();
        ofs.extend(dic.inicio[1..].iter().map(|&x| x + nt as u64));
        w.u64s(&ofs).map_err(er)?;
        let singulares: Vec<u32> = (0..nt as u32).collect();
        w.u32s(&singulares).map_err(er)?;
        w.u32s(&dic.dados).map_err(er)?;
        w.u32s(&diretorio).map_err(er)?;
        for (tmp, tamanho) in [(&temp_chaves, n_kmers * 8), (&temp_classes, n_kmers * 4)] {
            let mut f = File::open(tmp).map_err(|e| format!("Error opening {tmp}: {e}"))?;
            let copiados = std::io::copy(&mut f, &mut w.w).map_err(er)?;
            if copiados != tamanho {
                return Err(format!("Error assembling {output}: the temporary file {tmp} has {copiados} bytes, {tamanho} expected"));
            }
            w.n += copiados;
        }
        w.u64s(&seqs.inicio).map_err(er)?;
        w.bytes(&seqs.bases).map_err(er)?;
        w.u64s(&seqs.ns).map_err(er)?;
        w.bytes(ASSINATURA_FIM).map_err(er)?;
        let total = w.n + 8;
        w.u64(total).map_err(er)?;
        w.w.flush().map_err(er)?;
    }
    let _ = std::fs::remove_file(&temp_chaves);
    let _ = std::fs::remove_file(&temp_classes);
    std::fs::rename(&parcial, output).map_err(|e| format!("Error renaming {parcial} to {output}: {e}"))?;
    let bytes = std::fs::metadata(output).map(|m| m.len()).unwrap_or(0);
    if curtos > 0 {
        avisos.push(format!("{curtos} transcript(s) are shorter than k = {k}: they have no k-mer, and no read can be assigned to them"));
    }
    if cortados > 0 {
        avisos.push(format!("{cortados} transcript(s) end in a poly-A tail (10 or more A's): the tail was left out of their k-mers, as kallisto does"));
    }
    if tr.renomeados > 0 {
        avisos.push(format!("{} repeated transcript name(s) were renamed NAME_2, NAME_3... (make_unique)", tr.renomeados));
    }
    if !seqs.ns.is_empty() {
        avisos.push(format!(
            "{} base(s) are not A, C, G or T (N or IUPAC codes): the k-mers that contain them are left out of the index",
            milhar(seqs.ns.len() as u64)
        ));
    }
    Ok(IndexSummary {
        fasta: fasta.to_string(),
        index: output.to_string(),
        format_version: INDEX_FORMAT_VERSION,
        k,
        transcripts: nt,
        genes: genes.ids.len(),
        gene_map: genes.fonte,
        transcripts_without_gene: genes.sem_gene,
        bases: seqs.n,
        non_acgt_bases: seqs.ns.len() as u64,
        short_transcripts: curtos,
        kmers: n_kmers,
        equivalence_classes: n_ecs,
        renamed: tr.renomeados,
        poly_a_clipped: cortados,
        bytes,
        passes: passes as usize,
        seconds: inicio.elapsed().as_secs_f64(),
        warnings: avisos,
    })
}

/// What an index file holds, from its header (fast: the k-mers are not read).
pub fn transcript_index_info(caminho: &str) -> Result<IndexSummary, String> {
    let f = File::open(caminho).map_err(|e| format!("Error opening the transcript index {caminho}: {e}"))?;
    let total = f.metadata().map_err(|e| format!("Error reading {caminho}: {e}"))?.len();
    let mut r = LeBin { r: BufReader::new(f), pos: 0, total, caminho: caminho.to_string() };
    let c = le_cabecalho(&mut r)?;
    Ok(IndexSummary {
        fasta: c.fasta,
        index: caminho.to_string(),
        format_version: c.versao,
        k: c.k,
        transcripts: c.t as usize,
        genes: c.g as usize,
        gene_map: c.fonte_genes,
        bases: c.bases,
        non_acgt_bases: c.x,
        kmers: c.n,
        equivalence_classes: c.e,
        bytes: total,
        ..Default::default()
    })
}

// ---------------------------------------------------------------------
//  FRAGMENT LENGTHS, EFFECTIVE LENGTHS, EM
// ---------------------------------------------------------------------

/// A normal density on 0..1000 (kallisto's truncated Gaussian for
/// single-end reads); unnormalized -- only ratios are used.
pub fn gaussian_fld(media: f64, dp: f64) -> Vec<f64> {
    (0..MAX_FRAG_LEN).map(|i| (-(i as f64 - media).powi(2) / (2.0 * dp * dp)).exp()).collect()
}

/// Effective lengths: length - E[fragment | fragment <= length] + 1, with the
/// fragment-length distribution `fld` over 0..1000 (kallisto's
/// get_frag_len_means + calc_eff_lens). At least 1.
pub fn effective_lengths(comprimentos: &[u32], fld: &[f64]) -> Vec<f64> {
    let mut medias = vec![0.0; fld.len().max(1)];
    let (mut massa, mut soma) = (0.0, 0.0);
    for (i, &f) in fld.iter().enumerate() {
        massa += f * i as f64;
        soma += f;
        medias[i] = if soma > 0.0 { massa / soma } else { 0.0 };
    }
    let ultimo = medias.len() - 1;
    comprimentos
        .iter()
        .map(|&l| {
            let m = medias[(l as usize).min(ultimo)];
            let e = l as f64 - m + 1.0;
            if e >= 1.0 {
                e
            } else {
                1.0
            }
        })
        .collect()
}

#[derive(Clone, Debug)]
pub struct EmOptions {
    pub min_rounds: usize,
    pub max_rounds: usize,
    /// A transcript has converged when its count changes by less than this
    /// fraction (kallisto: 0.01)...
    pub change: f64,
    /// ... or is below this many reads (kallisto: 0.01).
    pub change_limit: f64,
    /// After convergence, counts below a tenth of this are set to zero
    /// (kallisto: 1e-7).
    pub alpha_limit: f64,
}

impl Default for EmOptions {
    fn default() -> Self {
        EmOptions { min_rounds: 50, max_rounds: 10_000, change: 1e-2, change_limit: 1e-2, alpha_limit: 1e-7 }
    }
}

/// The EM of kallisto 0.46 on equivalence-class counts: class `c` holds the
/// transcripts `dados[ofs[c]..ofs[c+1]]` and `cont[c]` reads. Returns the
/// estimated counts per transcript and the rounds run.
pub fn em(n: usize, ofs: &[usize], dados: &[u32], cont: &[f64], eff: &[f64], opt: &EmOptions) -> (Vec<f64>, usize) {
    let mut fixo = vec![0.0; n];
    let mut multi: Vec<usize> = Vec::new();
    for c in 0..cont.len() {
        let ts = &dados[ofs[c]..ofs[c + 1]];
        if cont[c] <= 0.0 || ts.is_empty() {
            continue;
        }
        if ts.len() == 1 {
            fixo[ts[0] as usize] += cont[c];
        } else {
            multi.push(c);
        }
    }
    if n == 0 {
        return (Vec::new(), 0);
    }
    let w: Vec<f64> = eff.iter().map(|&e| 1.0 / e).collect();
    let mut alfa = vec![1.0 / n as f64; n];
    let mut prox = vec![0.0; n];
    let mut ultima = false;
    let mut i = 0;
    while i < opt.max_rounds {
        prox.copy_from_slice(&fixo);
        for &c in &multi {
            let ts = &dados[ofs[c]..ofs[c + 1]];
            let mut den = 0.0;
            for &t in ts {
                den += alfa[t as usize] * w[t as usize];
            }
            if den < TOLERANCIA {
                continue;
            }
            let cn = cont[c] / den;
            for &t in ts {
                let t = t as usize;
                prox[t] += w[t] * alfa[t] * cn;
            }
        }
        let mut mudam = 0usize;
        for t in 0..n {
            if prox[t] > opt.change_limit && ((prox[t] - alfa[t]).abs() / prox[t]) > opt.change {
                mudam += 1;
            }
        }
        std::mem::swap(&mut alfa, &mut prox);
        let para = mudam == 0 && i > opt.min_rounds;
        if ultima {
            break;
        }
        if para {
            ultima = true;
            for a in alfa.iter_mut() {
                if *a < opt.alpha_limit / 10.0 {
                    *a = 0.0;
                }
            }
        }
        i += 1;
    }
    (alfa, i)
}

/// TPM from estimated counts and effective lengths.
pub fn tpm(contagens: &[f64], eff: &[f64]) -> Vec<f64> {
    let rho: Vec<f64> = contagens.iter().zip(eff).map(|(&c, &e)| if e < TOLERANCIA { 0.0 } else { c / e }).collect();
    let s: f64 = rho.iter().sum();
    if s <= 0.0 {
        return vec![0.0; rho.len()];
    }
    rho.iter().map(|&r| r / s * 1e6).collect()
}

// ---------------------------------------------------------------------
//  READING FASTQ IN BATCHES
// ---------------------------------------------------------------------

#[derive(Default)]
struct Lote {
    seq: Vec<u8>,
    fim: Vec<usize>,
    nome: Vec<u8>,
    fim_nome: Vec<usize>,
}

impl Lote {
    fn n(&self) -> usize {
        self.fim.len()
    }
    fn seq(&self, i: usize) -> &[u8] {
        let a = if i == 0 { 0 } else { self.fim[i - 1] };
        &self.seq[a..self.fim[i]]
    }
    fn nome(&self, i: usize) -> &[u8] {
        let a = if i == 0 { 0 } else { self.fim_nome[i - 1] };
        &self.nome[a..self.fim_nome[i]]
    }
}

/// The read name used to match mates: the first word of the header.
fn nome_do_par(cab: &[u8]) -> &[u8] {
    let c = &cab[1.min(cab.len())..];
    let fim = c.iter().position(|&x| x == b' ' || x == b'\t').unwrap_or(c.len());
    &c[..fim]
}

/// Two mates of one pair: the same name, or names that differ only by a
/// final 1 / 2 after '/', '.', '_' or ':' (`read/1` and `read/2`; SRA's
/// `--readids`: `SRR1.7.1` and `SRR1.7.2`).
fn mesmo_par(a: &[u8], b: &[u8]) -> bool {
    if a == b {
        return true;
    }
    let n = a.len();
    n >= 2 && b.len() == n && a[..n - 1] == b[..n - 1] && a[n - 1] == b'1' && b[n - 1] == b'2' && matches!(a[n - 2], b'/' | b'.' | b'_' | b':')
}

struct LeitorFastq {
    caminho: String,
    r: Box<dyn BufRead>,
    linha: u64,
    l: [Vec<u8>; 4],
}

impl LeitorFastq {
    fn le_linha(&mut self, i: usize) -> Result<bool, String> {
        self.l[i].clear();
        let n = self.r.read_until(b'\n', &mut self.l[i]).map_err(|e| {
            format!(
                "Error reading {} after line {}: {e}. The file may be truncated or not a FASTQ file (plain or .gz).",
                self.caminho, self.linha
            )
        })?;
        if n == 0 {
            return Ok(false);
        }
        self.linha += 1;
        while matches!(self.l[i].last(), Some(b'\n' | b'\r')) {
            self.l[i].pop();
        }
        Ok(true)
    }

    /// One record into the batch; false at the end of the file.
    fn proximo(&mut self, lote: &mut Lote, com_nome: bool) -> Result<bool, String> {
        loop {
            if !self.le_linha(0)? {
                return Ok(false);
            }
            if !self.l[0].iter().all(|c| c.is_ascii_whitespace()) {
                break;
            }
        }
        let l0 = self.linha;
        if self.l[0][0] != b'@' {
            if self.l[0][0] == b'>' {
                return Err(format!(
                    "{}, line {l0}: this is a FASTA file; reads are read from FASTQ (records of 4 lines: @name, sequence, +, qualities), plain or .gz",
                    self.caminho
                ));
            }
            let amostra = String::from_utf8_lossy(&self.l[0][..self.l[0].len().min(40)]).to_string();
            return Err(format!(
                "{}, line {l0}: a FASTQ record starts with '@', found \"{amostra}\". The file is not FASTQ, or a record above it has a wrong number of lines.",
                self.caminho
            ));
        }
        for i in 1..4 {
            if !self.le_linha(i)? {
                return Err(format!(
                    "{}: the last record (line {l0}) is incomplete -- a FASTQ record has 4 lines. The file is truncated (an interrupted download or copy?).",
                    self.caminho
                ));
            }
        }
        if self.l[2].first() != Some(&b'+') {
            return Err(format!(
                "{}, line {}: the third line of a FASTQ record starts with '+'; the record at line {l0} is malformed",
                self.caminho,
                l0 + 2
            ));
        }
        if self.l[1].len() != self.l[3].len() {
            return Err(format!(
                "{}, line {l0}: the sequence has {} bases and the quality line {} characters; they must be equal. The file is damaged.",
                self.caminho,
                self.l[1].len(),
                self.l[3].len()
            ));
        }
        lote.seq.extend_from_slice(&self.l[1]);
        lote.fim.push(lote.seq.len());
        if com_nome {
            let nm = nome_do_par(&self.l[0]);
            lote.nome.extend_from_slice(nm);
            lote.fim_nome.push(lote.nome.len());
        }
        Ok(true)
    }
}

fn produtor(caminho: String, tx: SyncSender<Result<Lote, String>>, limite: u64, com_nome: bool) {
    let r = match abre_leitura(&caminho) {
        Ok(r) => r,
        Err(e) => {
            let _ = tx.send(Err(e));
            return;
        }
    };
    let mut leitor = LeitorFastq { caminho, r, linha: 0, l: Default::default() };
    let mut lidos = 0u64;
    loop {
        let mut lote = Lote::default();
        while lote.n() < LOTE && lidos < limite {
            match leitor.proximo(&mut lote, com_nome) {
                Ok(true) => lidos += 1,
                Ok(false) => break,
                Err(e) => {
                    let _ = tx.send(Err(e));
                    return;
                }
            }
        }
        if lote.n() == 0 {
            return;
        }
        let cheio = lote.n() == LOTE;
        if tx.send(Ok(lote)).is_err() || !cheio {
            return;
        }
    }
}

fn recebe(rx: &Receiver<Result<Lote, String>>) -> Result<Option<Lote>, String> {
    match rx.recv() {
        Ok(Ok(l)) => Ok(Some(l)),
        Ok(Err(e)) => Err(e),
        Err(_) => Ok(None),
    }
}

// ---------------------------------------------------------------------
//  QUANTIFICATION
// ---------------------------------------------------------------------

/// One sample: single-end (`r2` = None) or paired-end. `condition`, when every
/// sample has one, goes to a sample table ready for differential expression.
#[derive(Clone, Debug, Default)]
pub struct QuantSample {
    pub name: String,
    pub r1: String,
    pub r2: Option<String>,
    pub condition: Option<String>,
}

#[derive(Clone, Debug)]
pub struct QuantOptions {
    /// Mean fragment length. Required for single-end reads; for pairs, given
    /// = used instead of the estimate.
    pub fragment_length: Option<f64>,
    pub fragment_sd: Option<f64>,
    /// Pairs on a single transcript used to estimate the fragment length
    /// (the first ones in the file). Default 10,000 (kallisto's).
    pub fragment_pairs: usize,
    /// Replaces the gene map stored in the index (GTF, GFF3 or table).
    pub gene_map: Option<String>,
    /// Stop when the two files of a pair are out of step (the read names
    /// differ). Default true.
    pub check_pair_names: bool,
    /// Single-end reads: keep a read even when, on some of its transcripts, a
    /// fragment of the mean length starting at the read would run past the
    /// transcript end. `false` applies kallisto's default rule (drop those
    /// transcripts from the read's class, and the read if none is left);
    /// `true` is kallisto's --single-overhang.
    pub single_overhang: bool,
    /// Only the first reads (pairs) of each sample: for a quick look.
    pub max_reads: Option<u64>,
    /// Keep the reads per equivalence class of every sample in the result
    /// (`SampleQuant::class_counts`): for comparisons with kallisto's
    /// pseudoalignment output; large for a mammalian index.
    pub keep_class_counts: bool,
    pub em: EmOptions,
}

impl Default for QuantOptions {
    fn default() -> Self {
        QuantOptions {
            fragment_length: None,
            fragment_sd: None,
            fragment_pairs: 10_000,
            gene_map: None,
            check_pair_names: true,
            single_overhang: true,
            max_reads: None,
            keep_class_counts: false,
            em: EmOptions::default(),
        }
    }
}

/// The result of one sample.
#[derive(Clone, Debug, Default)]
pub struct SampleQuant {
    pub sample: String,
    pub r1: String,
    pub r2: Option<String>,
    pub condition: Option<String>,
    pub paired: bool,
    /// Reads, or pairs, processed.
    pub reads: u64,
    pub pseudoaligned: u64,
    /// Pseudoaligned to a single transcript (kallisto's n_unique).
    pub unique: u64,
    /// Pseudoaligned to several transcripts.
    pub multi: u64,
    /// With k-mers, none of them in the index.
    pub no_kmer: u64,
    /// With k-mers in the index whose classes share no transcript.
    pub conflicting: u64,
    /// Without a single complete k-mer (shorter than k, or N-rich).
    pub too_short: u64,
    /// Pairs in which only one mate had k-mers in the index.
    pub one_mate_only: u64,
    /// Single-end reads dropped because no compatible transcript could hold
    /// a fragment of the mean length (only with `single_overhang` off).
    pub overhang: u64,
    pub equivalence_classes: usize,
    pub fragment_length_mean: f64,
    pub fragment_length_sd: f64,
    /// "estimated" or "given".
    pub fragment_length_source: String,
    pub fragment_pairs: u64,
    pub em_rounds: usize,
    pub seconds: f64,
    pub folder: String,
    /// Per transcript, in index order.
    pub est_counts: Vec<f64>,
    pub tpm: Vec<f64>,
    pub eff_lengths: Vec<f64>,
    /// Per gene, in the order of `QuantResult::gene_ids`.
    pub gene_counts: Vec<f64>,
    pub gene_tpm: Vec<f64>,
    /// Reads per equivalence class: (transcripts, reads), index classes in id
    /// order, then the new intersections in lexicographic order. Filled only
    /// with `keep_class_counts`.
    pub class_counts: Vec<(Vec<u32>, u64)>,
}

fn pct(x: u64, de: u64) -> f64 {
    if de == 0 {
        0.0
    } else {
        (1000.0 * x as f64 / de as f64).round() / 10.0
    }
}

impl SampleQuant {
    pub fn json(&self) -> String {
        format!(
            "{{\"sample\":{},\"r1\":{},\"r2\":{},\"condition\":{},\"paired\":{},\"reads\":{},\"pseudoaligned\":{},\"pseudoaligned_pct\":{},\"unique\":{},\"unique_pct\":{},\"multi\":{},\"no_kmer\":{},\"conflicting\":{},\"too_short\":{},\"one_mate_only\":{},\"overhang\":{},\"equivalence_classes\":{},\"fragment_length_mean\":{},\"fragment_length_sd\":{},\"fragment_length_source\":{},\"fragment_pairs\":{},\"em_rounds\":{},\"seconds\":{},\"folder\":{}}}",
            json_texto(&self.sample),
            json_texto(&self.r1),
            jopt(self.r2.as_deref()),
            jopt(self.condition.as_deref()),
            self.paired,
            self.reads,
            self.pseudoaligned,
            jnum(pct(self.pseudoaligned, self.reads)),
            self.unique,
            jnum(pct(self.unique, self.reads)),
            self.multi,
            self.no_kmer,
            self.conflicting,
            self.too_short,
            self.one_mate_only,
            self.overhang,
            self.equivalence_classes,
            jnum((self.fragment_length_mean * 100.0).round() / 100.0),
            jnum((self.fragment_length_sd * 100.0).round() / 100.0),
            json_texto(&self.fragment_length_source),
            self.fragment_pairs,
            self.em_rounds,
            jnum((self.seconds * 1000.0).round() / 1000.0),
            json_texto(&self.folder)
        )
    }
}

#[derive(Clone, Debug, Default)]
pub struct QuantResult {
    pub index: String,
    pub k: usize,
    pub output_dir: String,
    pub transcripts: Vec<String>,
    pub gene_ids: Vec<String>,
    pub gene_names: Vec<String>,
    pub gene_map: String,
    pub samples: Vec<SampleQuant>,
    /// (what, path) of every table written.
    pub files: Vec<(String, String)>,
    pub warnings: Vec<String>,
}

impl QuantResult {
    pub fn json(&self) -> String {
        let matriz = |f: &dyn Fn(&SampleQuant) -> &Vec<f64>| -> String {
            format!("[{}]", self.samples.iter().map(|s| jnums(f(s))).collect::<Vec<_>>().join(","))
        };
        let arquivos = self.files.iter().map(|(k, v)| format!("{}:{}", json_texto(k), json_texto(v))).collect::<Vec<_>>().join(",");
        format!(
            "{{\"index\":{},\"k\":{},\"output_dir\":{},\"transcripts\":{},\"genes\":{},\"gene_map\":{},\"samples\":[{}],\"gene_ids\":{},\"gene_names\":{},\"gene_counts\":{},\"gene_tpm\":{},\"files\":{{{}}},\"warnings\":{},\"text\":{}}}",
            json_texto(&self.index),
            self.k,
            json_texto(&self.output_dir),
            self.transcripts.len(),
            self.gene_ids.len(),
            json_texto(&self.gene_map),
            self.samples.iter().map(SampleQuant::json).collect::<Vec<_>>().join(","),
            jtextos(&self.gene_ids),
            jtextos(&self.gene_names),
            matriz(&|s| &s.gene_counts),
            matriz(&|s| &s.gene_tpm),
            arquivos,
            jtextos(&self.warnings),
            json_texto(&self.resumo())
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!(
            "Quantification of {} sample(s) against {} (k = {}, {} transcripts, {} genes; gene map: {}), results in {}\n",
            self.samples.len(),
            self.index,
            self.k,
            milhar(self.transcripts.len() as u64),
            milhar(self.gene_ids.len() as u64),
            self.gene_map,
            self.output_dir
        );
        for s in &self.samples {
            o.push_str(&format!(
                "  {}: {} {} -- pseudoaligned {} ({:.1}%), unique {} ({:.1}%); not assigned: {} without k-mers in the index, {} conflicting, {} too short{}; {} classes; fragment length {:.1} +/- {:.1} ({}{}); EM {} rounds\n",
                s.sample,
                milhar(s.reads),
                if s.paired { "pairs" } else { "reads" },
                milhar(s.pseudoaligned),
                pct(s.pseudoaligned, s.reads),
                milhar(s.unique),
                pct(s.unique, s.reads),
                milhar(s.no_kmer),
                milhar(s.conflicting),
                milhar(s.too_short),
                if s.overhang > 0 { format!(", {} with no transcript long enough for the fragment", milhar(s.overhang)) } else { String::new() },
                milhar(s.equivalence_classes as u64),
                s.fragment_length_mean,
                s.fragment_length_sd,
                s.fragment_length_source,
                if s.fragment_length_source == "estimated" { format!(" from {} pairs", milhar(s.fragment_pairs)) } else { String::new() },
                s.em_rounds
            ));
        }
        for w in &self.warnings {
            o.push_str(&format!("  warning: {w}\n"));
        }
        o
    }
}

struct Leitura {
    reads: u64,
    pseudo: u64,
    unicos: u64,
    sem_kmer: u64,
    conflito: u64,
    curtas: u64,
    um_mate: u64,
    fora: u64,
    ec_contagem: Vec<u64>,
    novos: HashMap<Vec<u32>, u64, Rapido>,
    flens: Vec<u64>,
    pares_frag: u64,
}

struct Res {
    c: Classe,
    frag: u32,
    mates: u8,
}

/// Reads one sample and counts the reads per equivalence class.
fn pseudoalinha_amostra(idx: &TranscriptIndex, s: &QuantSample, opt: &QuantOptions, estima_frag: bool) -> Result<Leitura, String> {
    let paired = s.r2.is_some();
    let limite = opt.max_reads.unwrap_or(u64::MAX);
    let com_nome = paired && opt.check_pair_names;
    let nt = idx.transcripts() as u32;
    let mut lei = Leitura {
        reads: 0,
        pseudo: 0,
        unicos: 0,
        sem_kmer: 0,
        conflito: 0,
        curtas: 0,
        um_mate: 0,
        fora: 0,
        ec_contagem: vec![0; idx.classes()],
        novos: HashMap::default(),
        flens: vec![0; MAX_FRAG_LEN],
        pares_frag: 0,
    };
    std::thread::scope(|escopo| -> Result<(), String> {
        let (tx1, rx1) = sync_channel::<Result<Lote, String>>(2);
        let r1 = s.r1.clone();
        escopo.spawn(move || produtor(r1, tx1, limite, com_nome));
        let rx2 = match &s.r2 {
            Some(r2) => {
                let (tx2, rx2) = sync_channel::<Result<Lote, String>>(2);
                let r2 = r2.clone();
                escopo.spawn(move || produtor(r2, tx2, limite, com_nome));
                Some(rx2)
            }
            None => None,
        };
        loop {
            let l1 = recebe(&rx1)?;
            let l2 = match &rx2 {
                Some(rx) => recebe(rx)?,
                None => None,
            };
            let r2 = s.r2.as_deref().unwrap_or("");
            let descompasso = |curto: &str, longo: &str, n: u64| {
                format!(
                    "{curto} ends after {} reads, but {longo} has more: the two files of a pair must hold the same reads in the same order (is one of them truncated, or from another run?)",
                    milhar(n)
                )
            };
            let (l1, l2) = match (l1, l2) {
                (None, None) => break,
                (Some(a), None) if !paired => (a, None),
                (Some(a), Some(b)) => {
                    if a.n() < b.n() {
                        return Err(descompasso(&s.r1, r2, lei.reads + a.n() as u64));
                    }
                    if b.n() < a.n() {
                        return Err(descompasso(r2, &s.r1, lei.reads + b.n() as u64));
                    }
                    (a, Some(b))
                }
                (Some(_), None) => return Err(descompasso(r2, &s.r1, lei.reads)),
                (None, Some(_)) => return Err(descompasso(&s.r1, r2, lei.reads)),
            };
            let n = l1.n();
            if let Some(l2) = &l2 {
                if com_nome {
                    for i in 0..n {
                        if !mesmo_par(l1.nome(i), l2.nome(i)) {
                            return Err(format!(
                                "Pair {} of sample {}: the read is \"{}\" in {} and \"{}\" in {}. The two files of a pair must list the mates in the same order (as the sequencer writes them); were the files sorted, filtered separately, or mixed up? (check_pair_names=False skips this test)",
                                milhar(lei.reads + i as u64 + 1),
                                s.name,
                                String::from_utf8_lossy(l1.nome(i)),
                                s.r1,
                                String::from_utf8_lossy(l2.nome(i)),
                                s.r2.as_deref().unwrap_or("")
                            ));
                        }
                    }
                }
            }
            let quer = estima_frag && paired && lei.pares_frag < opt.fragment_pairs as u64;
            let filtra = !paired && !opt.single_overhang;
            let fl = opt.fragment_length.unwrap_or(0.0) as i64;
            let l2r = l2.as_ref();
            let res: Vec<Res> = (0..n)
                .into_par_iter()
                .with_min_len(64)
                .map_init(Rascunho::default, |r, i| {
                    let s1 = l1.seq(i);
                    let s2 = l2r.map(|l| l.seq(i));
                    let mut c = idx.pseudoalinha_com(s1, s2, r);
                    if filtra && matches!(c, Classe::Ec(_) | Classe::Novo(_)) {
                        if let Some(h) = r.primeiro[0] {
                            let ts = idx.transcripts_of(&c);
                            let mut buf = std::mem::take(&mut r.buf);
                            let mut fica = std::mem::take(&mut r.b);
                            idx.cabe_fragmento(&ts, h, fl, &mut buf, &mut fica);
                            if fica.len() < ts.len() {
                                c = match fica.len() {
                                    0 => Classe::Fora,
                                    1 => Classe::Ec(fica[0]),
                                    _ => Classe::Novo(fica.clone()),
                                };
                            }
                            r.buf = buf;
                            r.b = fica;
                        }
                    }
                    let mut frag = 0u32;
                    if quer {
                        if let Classe::Ec(t) = c {
                            if t < nt {
                                if let (Some(h1), Some(h2)) = (r.primeiro[0], r.primeiro[1]) {
                                    let mut buf = std::mem::take(&mut r.buf);
                                    frag = idx.fragmento(t, h1, h2, &mut buf).unwrap_or(0);
                                    r.buf = buf;
                                }
                            }
                        }
                    }
                    Res { c, frag, mates: r.mates }
                })
                .collect();
            for x in res {
                lei.reads += 1;
                if paired && x.mates == 1 {
                    lei.um_mate += 1;
                }
                match x.c {
                    Classe::Curta => lei.curtas += 1,
                    Classe::SemKmer => lei.sem_kmer += 1,
                    Classe::Conflito => lei.conflito += 1,
                    Classe::Fora => lei.fora += 1,
                    Classe::Ec(e) => {
                        lei.ec_contagem[e as usize] += 1;
                        lei.pseudo += 1;
                        if e < nt {
                            lei.unicos += 1;
                        }
                    }
                    Classe::Novo(v) => {
                        *lei.novos.entry(v).or_insert(0) += 1;
                        lei.pseudo += 1;
                    }
                }
                if x.frag > 0 && lei.pares_frag < opt.fragment_pairs as u64 {
                    lei.flens[x.frag as usize] += 1;
                    lei.pares_frag += 1;
                }
            }
        }
        Ok(())
    })?;
    Ok(lei)
}

fn valida_amostras(amostras: &[QuantSample], opt: &QuantOptions) -> Result<Vec<QuantSample>, String> {
    if amostras.is_empty() {
        return Err("No sample given: each sample is a name and its FASTQ file(s) -- one file for single-end reads, two (R1, R2) for pairs".into());
    }
    let mut v = Vec::with_capacity(amostras.len());
    let mut vistos: HashMap<String, usize> = HashMap::new();
    for (i, a) in amostras.iter().enumerate() {
        let mut a = a.clone();
        if a.r1.trim().is_empty() {
            return Err(format!("Sample {} has no FASTQ file (r1)", i + 1));
        }
        if a.r2.as_deref().is_some_and(|x| x.trim().is_empty()) {
            a.r2 = None;
        }
        if a.name.trim().is_empty() {
            a.name = nome_de_arquivo(&a.r1);
        }
        a.name = a.name.trim().to_string();
        if a.name.chars().any(|c| matches!(c, '/' | '\\' | ':' | '*' | '?' | '"' | '<' | '>' | '|' | '\t' | '\n' | '\r')) || a.name == "." || a.name == ".." {
            return Err(format!(
                "The sample name \"{}\" cannot be used: it names a folder of results and a column of the count matrix. Use letters, digits, '-', '_' and '.' (no / \\ : * ? \" < > |).",
                a.name
            ));
        }
        if let Some(j) = vistos.insert(a.name.clone(), i) {
            return Err(format!("Samples {} and {} are both named \"{}\": every sample needs a distinct name", j + 1, i + 1, a.name));
        }
        for p in std::iter::once(&a.r1).chain(a.r2.iter()) {
            if std::fs::metadata(p).is_err() {
                return Err(format!("Sample {}: the file {p} does not exist (or cannot be read)", a.name));
            }
        }
        if let Some(r2) = &a.r2 {
            if mesmo_arquivo(&a.r1, r2) {
                return Err(format!(
                    "Sample {}: R1 and R2 are the same file ({r2}). For single-end reads give only r1; for pairs, the two files of the pair.",
                    a.name
                ));
            }
        } else if opt.fragment_length.is_none() {
            return Err(format!(
                "Sample {} is single-end: the fragment length cannot be estimated from single reads, so give fragment_length and fragment_sd (the library's mean insert size and its standard deviation, as kallisto's -l and -s; from the sequencing provider's report or a Bioanalyzer trace -- 200 and 20 are common for Illumina mRNA libraries).",
                a.name
            ));
        }
        v.push(a);
    }
    Ok(v)
}

/// A sample name from a file name: `ctrl1_R1.fastq.gz` -> `ctrl1`.
pub fn nome_de_arquivo(caminho: &str) -> String {
    let base = caminho.rsplit(['/', '\\']).next().unwrap_or(caminho);
    let mut s = base.to_string();
    for ext in [".gz", ".bgz", ".lz4", ".fastq", ".fq", ".txt"] {
        if s.to_ascii_lowercase().ends_with(ext) {
            s.truncate(s.len() - ext.len());
        }
    }
    for suf in ["_R1_001", "_R1", "_1", ".R1", ".1", "_r1"] {
        if s.ends_with(suf) && s.len() > suf.len() {
            s.truncate(s.len() - suf.len());
            break;
        }
    }
    s
}

fn escreve_tabela(caminho: &str, cab: &str, linhas: impl Iterator<Item = String>) -> Result<(), String> {
    let mut w = Escritor::cria(caminho)?;
    let erro = |e: String| format!("Error writing {caminho}: {e}");
    w.linha(cab).map_err(erro)?;
    for l in linhas {
        w.linha(&l).map_err(erro)?;
    }
    w.termina().map_err(erro)
}

fn junta(pasta: &str, nome: &str) -> String {
    std::path::Path::new(pasta).join(nome).to_string_lossy().to_string()
}

/// Quantifies RNA-seq samples from their FASTQ files against an index made
/// by `build_transcript_index`. Writes, in `output_dir`: one folder per
/// sample (`abundance.tsv` as kallisto's -- target_id, length, eff_length,
/// est_counts, tpm --, `abundance_genes.tsv`, `run_info.json`) and the tables
/// of all samples: `gene_counts.tsv` (rounded estimated counts, the matrix
/// `differential_expression` reads), `gene_tpm.tsv`, `transcript_counts.tsv`,
/// `transcript_tpm.tsv`, `gene_names.tsv`, `quant_stats.tsv` and, when every
/// sample has a condition, `samples.tsv`.
pub fn quantify_samples(index: &str, samples: &[QuantSample], output_dir: &str, opt: &QuantOptions) -> Result<QuantResult, String> {
    match (opt.fragment_length, opt.fragment_sd) {
        (Some(m), Some(s)) => {
            if !(m > 0.0 && m < MAX_FRAG_LEN as f64) || !(s > 0.0 && s.is_finite()) {
                return Err(format!(
                    "fragment_length = {m}, fragment_sd = {s}: the mean must be between 0 and {MAX_FRAG_LEN} bp and the standard deviation positive"
                ));
            }
        }
        (Some(_), None) | (None, Some(_)) => {
            return Err("Give both fragment_length and fragment_sd (the mean fragment length and its standard deviation), or neither for paired-end reads (then they are estimated)".into())
        }
        (None, None) => {}
    }
    if opt.fragment_pairs == 0 && opt.fragment_length.is_none() && samples.iter().any(|s| s.r2.is_some()) {
        return Err("fragment_pairs = 0 leaves no pair to estimate the fragment length: give fragment_length and fragment_sd, or a positive fragment_pairs".into());
    }
    if output_dir.trim().is_empty() {
        return Err("No output folder given".into());
    }
    let amostras = valida_amostras(samples, opt)?;
    let mut idx = TranscriptIndex::load(index)?;
    let mut avisos = Vec::new();
    if let Some(m) = opt.gene_map.as_deref().filter(|m| !m.trim().is_empty()) {
        let g = atribui_genes(&idx.names, &[], Some(m), &mut avisos)?;
        idx.gene_of = g.de;
        idx.gene_ids = g.ids;
        idx.gene_names = g.nomes;
        idx.gene_map = g.fonte;
    }
    std::fs::create_dir_all(output_dir).map_err(|e| format!("Error creating the folder {output_dir}: {e}"))?;
    let nt = idx.transcripts();
    let ng = idx.gene_ids.len();
    let mut cache: Option<HashMap<&[u32], u32, Rapido>> = None;
    let mut resultados: Vec<SampleQuant> = Vec::with_capacity(amostras.len());
    let mut arquivos: Vec<(String, String)> = Vec::new();
    for s in &amostras {
        let t0 = Instant::now();
        let paired = s.r2.is_some();
        let estima = paired && opt.fragment_length.is_none();
        let lei = pseudoalinha_amostra(&idx, s, opt, estima)?;
        // --- classes: as do indice, depois as intersecoes novas --------------
        let mut cont_ec = lei.ec_contagem;
        let mut restos: Vec<(Vec<u32>, u64)> = Vec::new();
        if !lei.novos.is_empty() {
            let mapa = cache.get_or_insert_with(|| idx.mapa_ecs());
            for (v, c) in lei.novos {
                match mapa.get(v.as_slice()) {
                    Some(&e) => cont_ec[e as usize] += c,
                    None => restos.push((v, c)),
                }
            }
            let mut ordem: Vec<u32> = (0..restos.len() as u32).collect();
            multi_merge_sort_stable_by(&mut ordem, |a, b| restos[*a as usize].0.cmp(&restos[*b as usize].0));
            let mut r2: Vec<(Vec<u32>, u64)> = Vec::with_capacity(restos.len());
            let mut velhos: Vec<Option<(Vec<u32>, u64)>> = restos.into_iter().map(Some).collect();
            for i in ordem {
                if let Some(x) = velhos[i as usize].take() {
                    r2.push(x);
                }
            }
            restos = r2;
        }
        let mut ofs: Vec<usize> = vec![0];
        let mut dados: Vec<u32> = Vec::new();
        let mut cont: Vec<f64> = Vec::new();
        let mut class_counts: Vec<(Vec<u32>, u64)> = Vec::new();
        for (e, &c) in cont_ec.iter().enumerate() {
            if c > 0 {
                let ts = idx.ec(e as u32);
                dados.extend_from_slice(ts);
                ofs.push(dados.len());
                cont.push(c as f64);
                if opt.keep_class_counts {
                    class_counts.push((ts.to_vec(), c));
                }
            }
        }
        for (v, c) in restos {
            dados.extend_from_slice(&v);
            ofs.push(dados.len());
            cont.push(c as f64);
            if opt.keep_class_counts {
                class_counts.push((v, c));
            }
        }
        let n_classes = cont.len();
        // --- comprimento dos fragmentos ---------------------------------------
        let (fld, media, dp, fonte) = match (opt.fragment_length, opt.fragment_sd) {
            (Some(m), Some(sd)) => (gaussian_fld(m, sd), m, sd, "given"),
            _ => {
                if lei.pares_frag == 0 {
                    return Err(format!(
                        "Sample {}: no read pair could be placed on a single transcript, so the fragment length cannot be estimated ({} pairs read, {} pseudoaligned). Give fragment_length and fragment_sd, or check that the reads come from this transcriptome.",
                        s.name,
                        milhar(lei.reads),
                        milhar(lei.pseudo)
                    ));
                }
                let n: f64 = lei.flens.iter().sum::<u64>() as f64;
                let m: f64 = lei.flens.iter().enumerate().map(|(i, &c)| i as f64 * c as f64).sum::<f64>() / n;
                let v: f64 = lei.flens.iter().enumerate().map(|(i, &c)| (i as f64 - m).powi(2) * c as f64).sum::<f64>() / n;
                (lei.flens.iter().map(|&c| c as f64).collect(), m, v.sqrt(), "estimated")
            }
        };
        let eff = effective_lengths(&idx.lengths, &fld);
        let (alfa, rodadas) = em(nt, &ofs, &dados, &cont, &eff, &opt.em);
        let tp = tpm(&alfa, &eff);
        let mut gc = vec![0.0; ng];
        let mut gt = vec![0.0; ng];
        for t in 0..nt {
            let g = idx.gene_of[t] as usize;
            gc[g] += alfa[t];
            gt[g] += tp[t];
        }
        let pasta = junta(output_dir, &s.name);
        std::fs::create_dir_all(&pasta).map_err(|e| format!("Error creating the folder {pasta}: {e}"))?;
        let mut sq = SampleQuant {
            sample: s.name.clone(),
            r1: s.r1.clone(),
            r2: s.r2.clone(),
            condition: s.condition.clone(),
            paired,
            reads: lei.reads,
            pseudoaligned: lei.pseudo,
            unique: lei.unicos,
            multi: lei.pseudo - lei.unicos,
            no_kmer: lei.sem_kmer,
            conflicting: lei.conflito,
            too_short: lei.curtas,
            one_mate_only: lei.um_mate,
            overhang: lei.fora,
            equivalence_classes: n_classes,
            fragment_length_mean: media,
            fragment_length_sd: dp,
            fragment_length_source: fonte.into(),
            fragment_pairs: if fonte == "estimated" { lei.pares_frag } else { 0 },
            em_rounds: rodadas,
            seconds: 0.0,
            folder: pasta.clone(),
            est_counts: alfa,
            tpm: tp,
            eff_lengths: eff,
            gene_counts: gc,
            gene_tpm: gt,
            class_counts,
        };
        let ab = junta(&pasta, "abundance.tsv");
        escreve_tabela(
            &ab,
            "target_id\tlength\teff_length\test_counts\ttpm",
            (0..nt).map(|t| format!("{}\t{}\t{}\t{}\t{}", idx.names[t], idx.lengths[t], num(sq.eff_lengths[t]), num(sq.est_counts[t]), num(sq.tpm[t]))),
        )?;
        let mut n_tr = vec![0u32; ng];
        for &g in &idx.gene_of {
            n_tr[g as usize] += 1;
        }
        escreve_tabela(
            &junta(&pasta, "abundance_genes.tsv"),
            "gene_id\tgene_name\ttranscripts\test_counts\ttpm",
            (0..ng).map(|g| format!("{}\t{}\t{}\t{}\t{}", idx.gene_ids[g], idx.gene_names[g], n_tr[g], num(sq.gene_counts[g]), num(sq.gene_tpm[g]))),
        )?;
        if lei.reads == 0 {
            avisos.push(format!("Sample {}: {} has no reads", s.name, s.r1));
        } else if lei.pseudo * 10 < lei.reads {
            avisos.push(format!(
                "Sample {}: only {:.1}% of the reads pseudoaligned. The reads may come from another species or another transcriptome, be mostly rRNA or genomic DNA, or the index was built with a different k.",
                s.name,
                pct(lei.pseudo, lei.reads)
            ));
        }
        sq.seconds = t0.elapsed().as_secs_f64();
        std::fs::write(junta(&pasta, "run_info.json"), sq.json()).map_err(|e| format!("Error writing {pasta}/run_info.json: {e}"))?;
        resultados.push(sq);
    }
    // --- as tabelas de todas as amostras ----------------------------------
    let cab_amostras: String = resultados.iter().map(|s| format!("\t{}", s.sample)).collect();
    let p = junta(output_dir, "gene_counts.tsv");
    escreve_tabela(
        &p,
        &format!("gene{cab_amostras}"),
        (0..ng).map(|g| {
            let mut l = idx.gene_ids[g].clone();
            for s in &resultados {
                l.push_str(&format!("\t{}", s.gene_counts[g].round() as u64));
            }
            l
        }),
    )?;
    arquivos.push(("gene_counts".into(), p));
    let p = junta(output_dir, "gene_tpm.tsv");
    escreve_tabela(
        &p,
        &format!("gene{cab_amostras}"),
        (0..ng).map(|g| {
            let mut l = idx.gene_ids[g].clone();
            for s in &resultados {
                l.push('\t');
                l.push_str(&num(s.gene_tpm[g]));
            }
            l
        }),
    )?;
    arquivos.push(("gene_tpm".into(), p));
    let p = junta(output_dir, "transcript_counts.tsv");
    escreve_tabela(
        &p,
        &format!("transcript\tgene{cab_amostras}"),
        (0..nt).map(|t| {
            let mut l = format!("{}\t{}", idx.names[t], idx.gene_ids[idx.gene_of[t] as usize]);
            for s in &resultados {
                l.push('\t');
                l.push_str(&num(s.est_counts[t]));
            }
            l
        }),
    )?;
    arquivos.push(("transcript_counts".into(), p));
    let p = junta(output_dir, "transcript_tpm.tsv");
    escreve_tabela(
        &p,
        &format!("transcript\tgene{cab_amostras}"),
        (0..nt).map(|t| {
            let mut l = format!("{}\t{}", idx.names[t], idx.gene_ids[idx.gene_of[t] as usize]);
            for s in &resultados {
                l.push('\t');
                l.push_str(&num(s.tpm[t]));
            }
            l
        }),
    )?;
    arquivos.push(("transcript_tpm".into(), p));
    let p = junta(output_dir, "gene_names.tsv");
    escreve_tabela(&p, "#gene_id\tgene_name", (0..ng).map(|g| format!("{}\t{}", idx.gene_ids[g], idx.gene_names[g])))?;
    arquivos.push(("gene_names".into(), p));
    let p = junta(output_dir, "quant_stats.tsv");
    escreve_tabela(
        &p,
        "sample\tpaired\treads\tpseudoaligned\tpseudoaligned_pct\tunique\tunique_pct\tmulti\tno_kmer\tconflicting\ttoo_short\tequivalence_classes\tfragment_length_mean\tfragment_length_sd\tfragment_length_source\tem_rounds",
        resultados.iter().map(|s| {
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{}\t{}",
                s.sample,
                s.paired,
                s.reads,
                s.pseudoaligned,
                pct(s.pseudoaligned, s.reads),
                s.unique,
                pct(s.unique, s.reads),
                s.multi,
                s.no_kmer,
                s.conflicting,
                s.too_short,
                s.equivalence_classes,
                s.fragment_length_mean,
                s.fragment_length_sd,
                s.fragment_length_source,
                s.em_rounds
            )
        }),
    )?;
    arquivos.push(("stats".into(), p));
    if resultados.iter().all(|s| s.condition.as_deref().is_some_and(|c| !c.trim().is_empty())) {
        let p = junta(output_dir, "samples.tsv");
        escreve_tabela(
            &p,
            "sample\tcondition",
            resultados.iter().map(|s| format!("{}\t{}", s.sample, s.condition.as_deref().unwrap_or("").trim())),
        )?;
        arquivos.push(("samples".into(), p));
    }
    Ok(QuantResult {
        index: index.to_string(),
        k: idx.k,
        output_dir: output_dir.to_string(),
        transcripts: idx.names.clone(),
        gene_ids: idx.gene_ids.clone(),
        gene_names: idx.gene_names.clone(),
        gene_map: idx.gene_map.clone(),
        samples: resultados,
        files: arquivos,
        warnings: avisos,
    })
}

// ---------------------------------------------------------------------
//  TEACHING: HOW ONE READ IS PSEUDOALIGNED
// ---------------------------------------------------------------------

/// Every k-mer of one read (or pair) with the class it hits, and the
/// intersection: what pseudoalignment does, step by step. Loads the index
/// (meant for small teaching indexes). JSON.
pub fn explain_pseudoalignment(index: &str, read1: &str, read2: Option<&str>) -> Result<String, String> {
    let idx = TranscriptIndex::load(index)?;
    let nomes_de = |ts: &[u32]| -> String {
        let mut v: Vec<String> = ts.iter().take(25).map(|&t| idx.names[t as usize].clone()).collect();
        if ts.len() > 25 {
            v.push(format!("... {} more", ts.len() - 25));
        }
        jtextos(&v)
    };
    let mut mates = Vec::new();
    for (m, s) in std::iter::once(read1).chain(read2).enumerate() {
        let s = s.trim();
        if s.is_empty() {
            return Err(format!("Read {} is empty", m + 1));
        }
        let mut ks = Vec::new();
        kmers(s.as_bytes(), idx.k, |pos, f, rc| {
            let c = idx.lookup(f.min(rc));
            ks.push(format!(
                "{{\"position\":{pos},\"kmer\":{},\"found\":{},\"class\":{},\"transcripts\":{}}}",
                json_texto(&decodifica_kmer(f, idx.k)),
                c.is_some(),
                c.map(|x| x.to_string()).unwrap_or_else(|| "null".into()),
                c.map(|x| nomes_de(idx.ec(x))).unwrap_or_else(|| "[]".into())
            ));
        });
        mates.push(format!("{{\"length\":{},\"kmers\":[{}]}}", s.len(), ks.join(",")));
    }
    let c = idx.pseudoalign(read1.trim().as_bytes(), read2.map(|r| r.trim().as_bytes()));
    let ts = idx.transcripts_of(&c);
    let estado = match &c {
        Classe::Curta => "too_short",
        Classe::SemKmer => "no_kmer",
        Classe::Conflito | Classe::Fora => "conflicting",
        _ if ts.len() == 1 => "unique",
        _ => "multi",
    };
    let mut genes: Vec<u32> = ts.iter().map(|&t| idx.gene_of[t as usize]).collect();
    genes.dedup();
    multi_merge_sort_stable(&mut genes);
    genes.dedup();
    let nomes_genes: Vec<String> = genes.iter().map(|&g| idx.gene_names[g as usize].clone()).collect();
    Ok(format!(
        "{{\"k\":{},\"mates\":[{}],\"status\":{},\"transcripts\":{},\"genes\":{}}}",
        idx.k,
        mates.join(","),
        json_texto(estado),
        nomes_de(&ts),
        jtextos(&nomes_genes)
    ))
}

// ---------------------------------------------------------------------
//  SMALL HELPERS
// ---------------------------------------------------------------------

fn mesmo_arquivo(a: &str, b: &str) -> bool {
    a == b || matches!((std::fs::canonicalize(a), std::fs::canonicalize(b)), (Ok(x), Ok(y)) if x == y)
}

fn recusa_mesmo_arquivo(entrada: &str, saida: &str) -> Result<(), String> {
    if mesmo_arquivo(entrada, saida) {
        return Err(format!(
            "The output `{saida}` is the input file `{entrada}` itself; writing it would destroy the input. Choose a different output path."
        ));
    }
    Ok(())
}

/// A number for a table: shortest exact form, scientific notation for very
/// small or very large values.
fn num(x: f64) -> String {
    if x.is_nan() {
        return "NA".into();
    }
    let a = x.abs();
    if a == 0.0 || (1e-4..1e15).contains(&a) {
        format!("{x}")
    } else {
        format!("{x:e}")
    }
}

fn jnum(x: f64) -> String {
    if x.is_finite() {
        num(x)
    } else {
        "null".into()
    }
}

fn jnums(v: &[f64]) -> String {
    format!("[{}]", v.iter().map(|x| jnum(*x)).collect::<Vec<_>>().join(","))
}

fn jtextos(v: &[String]) -> String {
    format!("[{}]", v.iter().map(|s| json_texto(s)).collect::<Vec<_>>().join(","))
}

fn jopt(s: Option<&str>) -> String {
    s.map(json_texto).unwrap_or_else(|| "null".into())
}

fn milhar(n: u64) -> String {
    let s = n.to_string();
    let mut o = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            o.push(',');
        }
        o.push(c);
    }
    o
}

fn mb(bytes: u64) -> String {
    if bytes >= 1 << 30 {
        format!("{:.2} GB", bytes as f64 / (1u64 << 30) as f64)
    } else {
        format!("{:.1} MB", bytes as f64 / (1u64 << 20) as f64)
    }
}
