// =====================================================================
//  ADVENTITIOUS AGENTS: SCREENING THE READS OF A BIOLOGICAL MATERIAL
//  FOR VIRUSES, BACTERIA, MYCOPLASMA AND FUNGI
//
//  FOR RESEARCH AND TEACHING. Releasing a GMP batch requires a validated
//  method run with validated software: ICH Q5A(R2) admits next-generation
//  sequencing for the detection of adventitious viruses, but the validation
//  (specificity, sensitivity, limit of detection, robustness) is the
//  laboratory's. Every result carries this notice (`NOTICE`).
//
//  Two tools:
//
//   1. `build_screen_database`  indexes reference genomes, each with a kind
//                               (virus, bacterium, mycoplasma, fungus,
//                               other, spike = a spike-in / process control,
//                               host), into one database file;
//   2. `contaminant_screen`     classifies the reads of each sample (single
//                               or paired, several lanes), sets the host
//                               reads apart, and reports for each organism
//                               the evidence and a call (detected / below
//                               threshold / background / not detected),
//                               compared with the negative controls; the
//                               recovery of the spike-ins and the limit of
//                               detection of the run; a traceable report
//                               (TSV, JSON, HTML) with the SHA-256 of every
//                               input and of the database, every parameter,
//                               the date and the engine version.
//
//  WHAT IS REUSED FROM microbio.rs, AND WHY THE INDEX IS NEW. The FASTA/FASTQ
//  reader (`SeqReader`, with its error messages, and through it .gz, BGZF and
//  .lz4) is reused. The sketch database of microbio.rs is not: (a) it keeps
//  hashes without positions, and the breadth of coverage needs to know WHERE
//  in the genome each matched k-mer lies; (b) FracMinHash samples k-mers at
//  random, so a read may hold none of them (Poisson), while minimizers
//  guarantee one sampled k-mer in every window of w k-mers; (c) its hash is
//  MurmurHash over the k-mer as text (a copy per k-mer) -- right for
//  comparing with sourmash, too slow for a human-genome host index; (d) it
//  costs 12 bytes per hash with no compact host-only set. Here a genome of
//  3 Gb fits a host set of ~0.7 GB.
//
//  THE K-MERS. Canonical k-mers (k odd, 31 by default) coded in 2 bits; the
//  hash is a bijection of the code (MurmurHash3's 64-bit finalizer after a
//  XOR), so two different k-mers never share a hash in the target index.
//  Minimizers: in every window of w consecutive k-mers the one with the
//  smallest hash is kept, and stored under its hash mixed once more
//  (`index_key`: a minimum crowds the low end of the hash range, the key
//  spreads the minimizers evenly over the buckets and the passes). A read
//  without errors holds only minimizers of
//  the genome it comes from (its windows are windows of the genome). k-mers
//  of low complexity (fewer than 8 distinct trinucleotides: microsatellites,
//  homopolymers) are dropped on both sides.
//
//  THE TARGET INDEX (every non-host reference). The distinct minimizers of
//  each reference with the genome window where they first occur; the window
//  is `genome_window_bp` (1,000 by default) or, for a small genome, the
//  length / 20 (at least 100 bp), so a virus has ~20 windows. One sorted
//  array of hashes (8 bytes) and one of occurrences (4 bytes: the window
//  when the minimizer is in a single reference, else a pointer to the list
//  of (reference, window)). Density 2/(w+1): with w = 12, ~15% of the
//  k-mers, ~1.9 bytes per reference base in memory.
//
//  THE HOST SET (every host reference). Membership only: the minimizers of
//  a sparser scheme (`host_minimizer_window`, 31 by default: one k-mer in
//  ~16) kept as 32-bit fingerprints in 2^B buckets (B from 8 to 24, by the
//  size): 4 bytes per minimizer, B + 32 bits of the hash compared (a false
//  match per lookup ~ n / 2^(B+32): ~2e-9 for a human genome). It is built
//  in passes over ranges of buckets, as many as `memory_mb` requires, so a
//  human genome indexes on a laptop: ~160 million minimizers, ~0.7 GB.
//
//  A READ (OR PAIR) is: no information (no valid k-mer), host (it shares at
//  least `host_min_hits` host minimizers and a larger fraction of its host
//  minimizers than of its target minimizers -- a read that looks like host
//  as much as like a microbe is host), unclassified (fewer than `min_hits`
//  target minimizers shared with any reference), weak (enough hits but less
//  than `min_fraction` of its minimizers: a distant relative of a
//  reference, or a chimera), unique to one reference (its hits beat every
//  other reference by more than `unique_margin`) or tied between several.
//
//  AN ORGANISM IS DETECTED when, from its own reads -- its UNIQUE reads and
//  the tied reads it wins (a tie set goes to its member with the most
//  unique reads, as the greedy cover of `gather` would give it):
//    - its reads >= min_reads;
//    - distinct genome windows hit >= min_windows;
//    - dispersion >= min_dispersion: windows hit / windows expected if
//      those reads fell at random on the genome, W(1 - (1 - c/W)^n), with c
//      the mean windows per read: a pile of reads on one shared region (a
//      vector element, an integrated viral fragment, a conserved gene) is
//      not a detection;
//    - its reads / (unique + tied) >= min_unique_fraction: a reference whose
//      reads all tie with a better-supported relative is not called;
//    - weak reads / (weak + its reads) <= max_weak_fraction: when most reads
//      that resemble it match it only weakly, the few that pass are the
//      best-matching reads of a relative that is not in the database;
//    - reads per million fragments >= min_rpm.
//  Negative controls (reagent blank, extraction blank): an organism also in
//  a blank is BACKGROUND unless its RPM in the sample is at least
//  `min_blank_ratio` (10 by default, the RPM ratio of Miller et al. 2019,
//  Genome Res 29:831) times its RPM in the blank. Process controls (kind
//  spike) go into the blanks too and are never background.
//
//  THE DATABASE FILE (version 1, little-endian)
//      "MMGSCRN1"  u32 version  u8 k  u8 0  u16 w  u16 host w  u16 0
//      u32 genome window  u32 host bucket bits B  u32 host passes
//      u32 references, each:
//        u8 kind, str name, str files ('\n'), str sha256 ('\n'),
//        str file sizes ('\n'), u64 sequences, u64 bases, u64 k-mers,
//        u64 minimizers, u32 window bp, u32 first window, u32 windows
//        (str = u32 length + UTF-8)
//      u64 n  keys[n] u64 (index_key of the minimizers, increasing)
//      occ[n] u32   u64 m  multi[m] u32
//      u64 h  dir[2^B + 1] u64  fingerprints[h] u32
//      "MMGSCEND"  SHA-256 of every byte before "MMGSCEND" (32 bytes)
// =====================================================================

use crate::microbio::{SeqReader, SeqRecord};
use crate::multimerge_core::{multi_merge_sort, multi_merge_sort_stable};
use rayon::prelude::*;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, VecDeque};
use std::fs::File;
use std::io::{BufReader, BufWriter, ErrorKind, Read, Write};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::{sync_channel, SyncSender};
use std::time::Instant;

/// The notice that goes with every result.
pub const NOTICE: &str = "For research and teaching only. This screen is not a validated method: releasing a GMP batch \
requires a validated method run with validated software. ICH Q5A(R2) admits next-generation sequencing for the detection \
of adventitious viruses, but the validation (specificity, sensitivity, limit of detection, robustness) is the laboratory's \
responsibility, and a detection here must be confirmed by an orthogonal method (PCR, culture). The authors take no \
responsibility for decisions made on the basis of these results.";

/// The name and version of this engine, as reports print it.
pub fn engine_version() -> String {
    format!("MMGEN contaminant screen 1.0 (contaminantes.rs; package {})", env!("CARGO_PKG_VERSION"))
}

const DB_MAGIC: &[u8; 8] = b"MMGSCRN1";
const DB_END: &[u8; 8] = b"MMGSCEND";
const DB_VERSION: u32 = 1;
/// The high bit of an occurrence: a pointer into `multi`, not a window.
const MULTI: u32 = 0x8000_0000;
/// Long sequences are cut in pieces of this size, done in parallel.
const CHUNK: usize = 1 << 20;
/// Pieces handed to the threads at a time (32 MB of sequence).
const CHUNKS_AT_ONCE: usize = 32;
const MIN_WINDOW_BP: u64 = 100;
/// A small genome gets about this many windows.
const SMALL_GENOME_WINDOWS: u64 = 20;
/// A k-mer with fewer distinct trinucleotides is low complexity.
const LOW_COMPLEXITY_TRIMERS: u32 = 8;
const MIX_SEED: u64 = 0x9e37_79b9_7f4a_7c15;
/// A second mixing of a minimizer's hash before it is stored or looked up.
const KEY_SEED: u64 = 0x2545_f491_4f6c_dd1d;

// ---------------------------------------------------------------------
//  KINDS OF REFERENCE
// ---------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum AgentKind {
    Virus,
    #[default]
    Bacterium,
    Mycoplasma,
    Fungus,
    /// Any other agent (a parasite, a protist...).
    Other,
    /// A spike-in or process control: reported, never called a contaminant.
    Spike,
    /// The host (human, bovine, murine, CHO...): its reads are set aside.
    Host,
}

impl AgentKind {
    pub fn from_text(s: &str) -> Result<AgentKind, String> {
        match s.trim().to_ascii_lowercase().replace(['-', ' '], "_").as_str() {
            "virus" | "viruses" | "viral" => Ok(AgentKind::Virus),
            "bacterium" | "bacteria" | "bacterial" => Ok(AgentKind::Bacterium),
            "mycoplasma" | "mycoplasmas" | "mollicute" | "mollicutes" => Ok(AgentKind::Mycoplasma),
            "fungus" | "fungi" | "fungal" | "yeast" | "mold" | "mould" => Ok(AgentKind::Fungus),
            "other" | "parasite" | "protist" | "protozoan" => Ok(AgentKind::Other),
            "spike" | "spike_in" | "spikein" | "control" | "process_control" => Ok(AgentKind::Spike),
            "host" => Ok(AgentKind::Host),
            o => Err(format!(
                "Unknown kind \"{o}\" for a reference. Use one of: virus, bacterium, mycoplasma, fungus, other, spike (a spike-in or process control), host"
            )),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            AgentKind::Virus => "virus",
            AgentKind::Bacterium => "bacterium",
            AgentKind::Mycoplasma => "mycoplasma",
            AgentKind::Fungus => "fungus",
            AgentKind::Other => "other",
            AgentKind::Spike => "spike",
            AgentKind::Host => "host",
        }
    }

    fn code(self) -> u8 {
        match self {
            AgentKind::Virus => 0,
            AgentKind::Bacterium => 1,
            AgentKind::Mycoplasma => 2,
            AgentKind::Fungus => 3,
            AgentKind::Other => 4,
            AgentKind::Spike => 5,
            AgentKind::Host => 6,
        }
    }

    fn from_code(c: u8) -> Option<AgentKind> {
        Some(match c {
            0 => AgentKind::Virus,
            1 => AgentKind::Bacterium,
            2 => AgentKind::Mycoplasma,
            3 => AgentKind::Fungus,
            4 => AgentKind::Other,
            5 => AgentKind::Spike,
            6 => AgentKind::Host,
            _ => return None,
        })
    }

    /// An adventitious agent (not host, not a control).
    pub fn is_agent(self) -> bool {
        !matches!(self, AgentKind::Spike | AgentKind::Host)
    }
}

// ---------------------------------------------------------------------
//  SMALL HELPERS: inputs, outputs, JSON, numbers, dates
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

fn write_error(path: &str, e: impl std::fmt::Display) -> String {
    format!("Error writing the output {path}: {e}. Check the free space and the permission to write there")
}

fn write_file(path: &str, text: &str) -> Result<(), String> {
    let mut f = File::create(path).map_err(|e| {
        format!("Error creating the output {path}: {e}. Check that the output folder exists and that you can write to it")
    })?;
    f.write_all(text.as_bytes()).map_err(|e| write_error(path, e))?;
    f.flush().map_err(|e| write_error(path, e))
}

/// Removes the outputs written so far unless disarmed.
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

/// A number for JSON, rounded to `d` decimals (null when not finite).
fn jr(x: f64, d: i32) -> String {
    if !x.is_finite() {
        return "null".to_string();
    }
    let f = 10f64.powi(d);
    let r = (x * f).round() / f;
    if r == 0.0 {
        "0".to_string()
    } else {
        format!("{r}")
    }
}

fn jo(x: Option<f64>, d: i32) -> String {
    x.map_or_else(|| "null".to_string(), |v| jr(v, d))
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

/// A short number for text: 3 significant digits, no exponent.
fn short(x: f64) -> String {
    if !x.is_finite() {
        return "-".to_string();
    }
    let a = x.abs();
    if a == 0.0 {
        "0".to_string()
    } else if a >= 100.0 {
        format!("{x:.0}")
    } else if a >= 10.0 {
        format!("{x:.1}")
    } else if a >= 1.0 {
        format!("{x:.2}")
    } else if a >= 0.01 {
        format!("{x:.3}")
    } else {
        format!("{x:.2e}")
    }
}

fn clean_cell(s: &str) -> String {
    s.replace(['\t', '\n', '\r'], " ")
}

fn html(s: &str) -> String {
    let mut o = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => o.push_str("&amp;"),
            '<' => o.push_str("&lt;"),
            '>' => o.push_str("&gt;"),
            '"' => o.push_str("&quot;"),
            '\'' => o.push_str("&#39;"),
            c => o.push(c),
        }
    }
    o
}

/// The current date and time, UTC, ISO 8601.
pub fn utc_now() -> String {
    let secs = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0);
    let days = (secs / 86_400) as i64;
    let rem = secs % 86_400;
    // Civil date from days since 1970-01-01 (H. Hinnant's algorithm).
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!("{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}Z", rem / 3600, rem % 3600 / 60, rem % 60)
}

// ---------------------------------------------------------------------
//  SHA-256, in a stream (the input files can be tens of GB)
// ---------------------------------------------------------------------

const K256: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be,
    0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa,
    0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85,
    0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f,
    0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

/// SHA-256 (FIPS 180-4), fed in pieces.
#[derive(Clone)]
pub struct Sha256 {
    state: [u32; 8],
    block: [u8; 64],
    used: usize,
    length: u64,
}

impl Default for Sha256 {
    fn default() -> Self {
        Sha256::new()
    }
}

fn sha_compress(state: &mut [u32; 8], block: &[u8]) {
    let mut w = [0u32; 64];
    for (i, c) in block.chunks_exact(4).enumerate() {
        w[i] = u32::from_be_bytes([c[0], c[1], c[2], c[3]]);
    }
    for i in 16..64 {
        let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
        let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16].wrapping_add(s0).wrapping_add(w[i - 7]).wrapping_add(s1);
    }
    let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = *state;
    for i in 0..64 {
        let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
        let ch = (e & f) ^ (!e & g);
        let t1 = h.wrapping_add(s1).wrapping_add(ch).wrapping_add(K256[i]).wrapping_add(w[i]);
        let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
        let maj = (a & b) ^ (a & c) ^ (b & c);
        let t2 = s0.wrapping_add(maj);
        h = g;
        g = f;
        f = e;
        e = d.wrapping_add(t1);
        d = c;
        c = b;
        b = a;
        a = t1.wrapping_add(t2);
    }
    for (s, v) in state.iter_mut().zip([a, b, c, d, e, f, g, h]) {
        *s = s.wrapping_add(v);
    }
}

impl Sha256 {
    pub fn new() -> Sha256 {
        Sha256 {
            state: [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19],
            block: [0; 64],
            used: 0,
            length: 0,
        }
    }

    pub fn update(&mut self, mut data: &[u8]) {
        self.length = self.length.wrapping_add(data.len() as u64);
        if self.used > 0 {
            let take = (64 - self.used).min(data.len());
            self.block[self.used..self.used + take].copy_from_slice(&data[..take]);
            self.used += take;
            data = &data[take..];
            if self.used == 64 {
                let b = self.block;
                sha_compress(&mut self.state, &b);
                self.used = 0;
            }
        }
        while data.len() >= 64 {
            sha_compress(&mut self.state, &data[..64]);
            data = &data[64..];
        }
        if !data.is_empty() {
            self.block[..data.len()].copy_from_slice(data);
            self.used = data.len();
        }
    }

    pub fn finalize(mut self) -> [u8; 32] {
        let bits = self.length.wrapping_mul(8);
        let pad_len = if self.used < 56 { 56 - self.used } else { 120 - self.used };
        let mut pad = [0u8; 72];
        pad[0] = 0x80;
        self.update(&pad[..pad_len]);
        self.update(&bits.to_be_bytes());
        let mut out = [0u8; 32];
        for (i, s) in self.state.iter().enumerate() {
            out[4 * i..4 * i + 4].copy_from_slice(&s.to_be_bytes());
        }
        out
    }
}

pub fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// SHA-256 of a byte string, in hexadecimal.
pub fn sha256_hex(data: &[u8]) -> String {
    let mut h = Sha256::new();
    h.update(data);
    hex(&h.finalize())
}

/// The size and the SHA-256 (hexadecimal) of a file, read in a stream.
pub fn sha256_file(path: &str) -> Result<(u64, String), String> {
    check_input(path)?;
    let mut f = File::open(path).map_err(|e| format!("Error opening {path}: {e}"))?;
    let mut h = Sha256::new();
    let mut buf = vec![0u8; 1 << 20];
    let mut n = 0u64;
    loop {
        let r = f.read(&mut buf).map_err(|e| format!("Error reading {path} to compute its SHA-256: {e}"))?;
        if r == 0 {
            break;
        }
        n += r as u64;
        h.update(&buf[..r]);
    }
    Ok((n, hex(&h.finalize())))
}

/// Size and SHA-256 of each distinct path, in order.
fn hash_files(paths: &[String]) -> Result<HashMap<String, (u64, String)>, String> {
    let mut out = HashMap::new();
    for p in paths {
        if !out.contains_key(p) {
            let v = sha256_file(p)?;
            out.insert(p.clone(), v);
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------
//  K-MERS AND MINIMIZERS
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

#[inline]
fn fmix64(mut k: u64) -> u64 {
    k ^= k >> 33;
    k = k.wrapping_mul(0xff51_afd7_ed55_8ccd);
    k ^= k >> 33;
    k = k.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    k ^= k >> 33;
    k
}

/// The hash of a canonical k-mer code: a bijection of 64-bit words.
#[inline]
pub fn kmer_hash(code: u64) -> u64 {
    fmix64(code ^ MIX_SEED)
}

/// The key under which a minimizer is stored and looked up: its hash
/// mixed again (a bijection). A minimizer is the SMALLEST hash of its
/// window, so minimizer hashes crowd the low end of the 64-bit range; the
/// key spreads them evenly again, as the bucket directories and the passes
/// of the host set need.
#[inline]
pub fn index_key(minimizer_hash: u64) -> u64 {
    fmix64(minimizer_hash ^ KEY_SEED)
}

/// Calls `f(hash, position)` for the minimizer of every window of `w`
/// consecutive valid k-mers of `seq` (canonical k-mers, k <= 31), once per
/// change of minimizer. A base other than A/C/G/T ends the run. Ties keep
/// the leftmost k-mer.
pub fn minimizers(seq: &[u8], k: usize, w: usize, dq: &mut VecDeque<(u64, usize)>, mut f: impl FnMut(u64, usize)) {
    dq.clear();
    if seq.len() < k || k == 0 || k > 31 || w == 0 {
        return;
    }
    let mask = (1u64 << (2 * k)) - 1;
    let shift = 2 * (k - 1);
    let (mut fw, mut rv) = (0u64, 0u64);
    let mut run = 0usize;
    let mut nk = 0usize;
    let mut last = usize::MAX;
    for (i, &b) in seq.iter().enumerate() {
        let c = CODE[usize::from(b)];
        if c > 3 {
            run = 0;
            nk = 0;
            dq.clear();
            continue;
        }
        let c = u64::from(c);
        fw = ((fw << 2) | c) & mask;
        rv = (rv >> 2) | ((3 - c) << shift);
        run += 1;
        if run < k {
            continue;
        }
        let pos = i + 1 - k;
        let h = kmer_hash(fw.min(rv));
        while let Some(&(bh, _)) = dq.back() {
            if bh > h {
                dq.pop_back();
            } else {
                break;
            }
        }
        dq.push_back((h, pos));
        nk += 1;
        if nk >= w {
            let start = pos + 1 - w;
            while let Some(&(_, p)) = dq.front() {
                if p < start {
                    dq.pop_front();
                } else {
                    break;
                }
            }
            if let Some(&(mh, mp)) = dq.front() {
                if mp != last {
                    last = mp;
                    f(mh, mp);
                }
            }
        }
    }
}

/// Fewer than 8 distinct trinucleotides in the k bases from `pos`: a
/// microsatellite or a homopolymer.
pub fn low_complexity(seq: &[u8], pos: usize, k: usize) -> bool {
    let mut seen = 0u64;
    let mut t = 0usize;
    for (j, &b) in seq[pos..pos + k].iter().enumerate() {
        t = ((t << 2) | usize::from(CODE[usize::from(b)] & 3)) & 63;
        if j >= 2 {
            seen |= 1u64 << t;
        }
    }
    seen.count_ones() < LOW_COMPLEXITY_TRIMERS
}

fn count_valid_kmers(seq: &[u8], k: usize) -> u64 {
    let (mut run, mut n) = (0usize, 0u64);
    for &b in seq {
        if CODE[usize::from(b)] > 3 {
            run = 0;
        } else {
            run += 1;
            if run >= k {
                n += 1;
            }
        }
    }
    n
}

/// The minimizers of one (possibly very long) sequence, low-complexity ones
/// dropped, mapped by `map` (which may drop them too) and handed to `sink`
/// in blocks, in position order. Long sequences are cut in pieces that
/// overlap by k + w - 2 bases (every window lies whole in some piece) and
/// done in parallel; a minimizer of the overlap may come twice.
fn sequence_minimizers<T, M, S, E>(seq: &[u8], k: usize, w: usize, map: M, mut sink: S) -> Result<(), E>
where
    T: Send,
    M: Fn(u64, usize) -> Option<T> + Sync,
    S: FnMut(Vec<T>) -> Result<(), E>,
{
    let span = k + w - 2;
    let piece = |s: usize, e: usize| -> Vec<T> {
        let mut dq = VecDeque::new();
        let mut v = Vec::new();
        let part = &seq[s..e];
        minimizers(part, k, w, &mut dq, |h, p| {
            if !low_complexity(part, p, k) {
                if let Some(x) = map(h, s + p) {
                    v.push(x);
                }
            }
        });
        v
    };
    if seq.len() <= CHUNK + span {
        return sink(piece(0, seq.len()));
    }
    let starts: Vec<usize> = (0..seq.len()).step_by(CHUNK).collect();
    for group in starts.chunks(CHUNKS_AT_ONCE) {
        let parts: Vec<Vec<T>> = group.par_iter().map(|&s| piece(s, (s + CHUNK + span).min(seq.len()))).collect();
        for p in parts {
            sink(p)?;
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------
//  THE DATABASE: references, options, summary
// ---------------------------------------------------------------------

/// One reference genome: one or more FASTA files (chromosomes, plasmids,
/// segments) and its kind.
#[derive(Clone, Debug, Default)]
pub struct ScreenReference {
    pub files: Vec<String>,
    /// The name in every report; empty: the first FASTA header.
    pub name: String,
    pub kind: AgentKind,
}

#[derive(Clone, Debug)]
pub struct ScreenDbOptions {
    /// k-mer length: odd, 15 to 31 (31 by default).
    pub k: u8,
    /// Minimizer window of the agents (k-mers): one k-mer in about
    /// (w + 1) / 2 is kept. 12 by default (~15% of the k-mers).
    pub minimizer_window: u16,
    /// Minimizer window of the host: sparser, for memory. 31 by default
    /// (~6% of the k-mers; a human genome ~0.7 GB).
    pub host_minimizer_window: u16,
    /// Genome window for the breadth of coverage, in bases (1,000 by
    /// default; a small genome gets length / 20, at least 100).
    pub genome_window_bp: u32,
    /// Memory for building the host set; more passes over the host FASTA
    /// when it does not fit (1,024 MB by default).
    pub memory_mb: u64,
}

impl Default for ScreenDbOptions {
    fn default() -> Self {
        ScreenDbOptions { k: 31, minimizer_window: 12, host_minimizer_window: 31, genome_window_bp: 1000, memory_mb: 1024 }
    }
}

/// A reference as the database keeps it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ScreenRefInfo {
    pub name: String,
    pub kind: AgentKind,
    pub files: Vec<String>,
    pub sha256: Vec<String>,
    pub file_bytes: Vec<u64>,
    pub sequences: u64,
    pub bases: u64,
    /// Valid k-mers (all A/C/G/T), with repeats.
    pub kmers: u64,
    /// Distinct minimizers it put in the index.
    pub minimizers: u64,
    /// Breadth windows (none for a host).
    pub window_bp: u32,
    pub first_window: u32,
    pub windows: u32,
}

/// What a database holds.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ScreenDbSummary {
    pub path: String,
    pub sha256: String,
    pub file_bytes: u64,
    pub k: u8,
    pub minimizer_window: u16,
    pub host_minimizer_window: u16,
    pub genome_window_bp: u32,
    pub references: Vec<ScreenRefInfo>,
    /// Distinct minimizers of the agents and controls.
    pub target_minimizers: u64,
    /// ... of which in more than one reference.
    pub shared_minimizers: u64,
    pub host_minimizers: u64,
    pub host_bucket_bits: u32,
    pub host_passes: u32,
    /// Memory the index takes when screening.
    pub target_index_bytes: u64,
    pub host_index_bytes: u64,
    /// Seconds to build (0 when the database was only read).
    pub seconds: f64,
    pub warnings: Vec<String>,
}

impl ScreenDbSummary {
    pub fn json(&self) -> String {
        let refs: Vec<String> = self
            .references
            .iter()
            .map(|r| {
                format!(
                    "{{\"name\":{},\"kind\":{},\"files\":{},\"sha256\":{},\"file_bytes\":[{}],\"sequences\":{},\"bases\":{},\"kmers\":{},\"minimizers\":{},\"window_bp\":{},\"windows\":{}}}",
                    js(&r.name),
                    js(r.kind.name()),
                    jlist(&r.files),
                    jlist(&r.sha256),
                    r.file_bytes.iter().map(|b| b.to_string()).collect::<Vec<_>>().join(","),
                    r.sequences,
                    r.bases,
                    r.kmers,
                    r.minimizers,
                    r.window_bp,
                    r.windows
                )
            })
            .collect();
        format!(
            "{{\"path\":{},\"sha256\":{},\"file_bytes\":{},\"k\":{},\"minimizer_window\":{},\"host_minimizer_window\":{},\"genome_window_bp\":{},\"target_minimizers\":{},\"shared_minimizers\":{},\"host_minimizers\":{},\"host_bucket_bits\":{},\"host_passes\":{},\"target_index_bytes\":{},\"host_index_bytes\":{},\"seconds\":{},\"references\":[{}],\"warnings\":{},\"text\":{},\"notice\":{}}}",
            js(&self.path),
            js(&self.sha256),
            self.file_bytes,
            self.k,
            self.minimizer_window,
            self.host_minimizer_window,
            self.genome_window_bp,
            self.target_minimizers,
            self.shared_minimizers,
            self.host_minimizers,
            self.host_bucket_bits,
            self.host_passes,
            self.target_index_bytes,
            self.host_index_bytes,
            jr(self.seconds, 2),
            refs.join(","),
            jlist(&self.warnings),
            js(&self.resumo()),
            js(NOTICE)
        )
    }

    pub fn resumo(&self) -> String {
        let mut o = format!(
            "Contaminant screen database {}: k = {}, minimizer window {} (agents) / {} (host), genome window {} bp\n  {} reference(s); {} agent minimizers ({} shared by several references), {} host minimizers\n  memory when screening: {:.1} MB (agents) + {:.1} MB (host); file {:.1} MB, SHA-256 {}\n",
            self.path,
            self.k,
            self.minimizer_window,
            self.host_minimizer_window,
            self.genome_window_bp,
            self.references.len(),
            thousands(self.target_minimizers),
            thousands(self.shared_minimizers),
            thousands(self.host_minimizers),
            self.target_index_bytes as f64 / 1e6,
            self.host_index_bytes as f64 / 1e6,
            self.file_bytes as f64 / 1e6,
            self.sha256
        );
        for r in &self.references {
            let win = if r.kind == AgentKind::Host { "host".to_string() } else { format!("{} windows of {} bp", r.windows, r.window_bp) };
            o.push_str(&format!(
                "  {:<10} {:>14} bases  {:>11} minimizers  {:<26} {}\n",
                r.kind.name(),
                thousands(r.bases),
                thousands(r.minimizers),
                win,
                r.name
            ));
        }
        for w in &self.warnings {
            o.push_str(&format!("  warning: {w}\n"));
        }
        o
    }
}

/// The host set: 32-bit fingerprints in buckets of the top `bits` bits.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HostIndex {
    pub bits: u32,
    pub dir: Vec<u64>,
    pub fp: Vec<u32>,
}

impl HostIndex {
    fn empty() -> HostIndex {
        HostIndex { bits: 8, dir: vec![0; 257], fp: Vec::new() }
    }

    #[inline]
    pub fn contains(&self, h: u64) -> bool {
        if self.fp.is_empty() {
            return false;
        }
        let b = (h >> (64 - self.bits)) as usize;
        let f = ((h << self.bits) >> 32) as u32;
        let (s, e) = (self.dir[b] as usize, self.dir[b + 1] as usize);
        self.fp[s..e].binary_search(&f).is_ok()
    }

    pub fn len(&self) -> usize {
        self.fp.len()
    }

    pub fn is_empty(&self) -> bool {
        self.fp.is_empty()
    }

    fn bytes(&self) -> u64 {
        (self.fp.len() * 4 + self.dir.len() * 8) as u64
    }
}

/// A database read back, ready to screen.
pub struct ScreenDb {
    pub summary: ScreenDbSummary,
    k: usize,
    w: usize,
    hw: usize,
    keys: Vec<u64>,
    occ: Vec<u32>,
    multi: Vec<u32>,
    tbits: u32,
    tdir: Vec<u32>,
    /// The reference of each global window.
    win_ref: Vec<u32>,
    total_windows: usize,
    host: HostIndex,
}

impl ScreenDb {
    pub fn references(&self) -> &[ScreenRefInfo] {
        &self.summary.references
    }

    #[inline]
    fn lookup(&self, h: u64) -> Option<usize> {
        if self.keys.is_empty() {
            return None;
        }
        let b = (h >> (64 - self.tbits)) as usize;
        let (s, e) = (self.tdir[b] as usize, self.tdir[b + 1] as usize);
        self.keys[s..e].binary_search(&h).ok().map(|i| s + i)
    }

    /// The references holding a minimizer (by its hash), with the global
    /// window of its first occurrence in each (testing aid).
    pub fn references_of(&self, h: u64) -> Vec<(u32, u32)> {
        let mut v = Vec::new();
        if let Some(i) = self.lookup(index_key(h)) {
            self.occurrences(i, |r, g, _| v.push((r, g)));
        }
        v
    }

    /// Calls `f(reference, global window, slot)` for each occurrence of key
    /// `i`. Slots number the (key, reference) pairs: the key itself for a
    /// key of one reference, n_keys + its place in `multi` otherwise.
    #[inline]
    fn occurrences(&self, i: usize, mut f: impl FnMut(u32, u32, usize)) {
        let o = self.occ[i];
        if o & MULTI == 0 {
            f(self.win_ref[o as usize], o, i);
        } else {
            let off = (o & !MULTI) as usize;
            let n = self.multi[off] as usize;
            for j in 1..=n {
                let g = self.multi[off + j];
                f(self.win_ref[g as usize], g, self.keys.len() + off + j);
            }
        }
    }

    /// Whether a minimizer (by its hash) is in the host set.
    pub fn host_contains(&self, h: u64) -> bool {
        self.host.contains(index_key(h))
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn minimizer_window(&self) -> usize {
        self.w
    }

    pub fn host_minimizer_window(&self) -> usize {
        self.hw
    }
}

fn check_db_options(opt: &ScreenDbOptions) -> Result<(), String> {
    if opt.k < 15 || opt.k > 31 || opt.k % 2 == 0 {
        return Err(format!(
            "Invalid k-mer length: k = {}. Use an odd k from 15 to 31 (31 is the usual one; an odd k means no k-mer is its own reverse complement)",
            opt.k
        ));
    }
    if opt.minimizer_window == 0 || opt.minimizer_window > 256 {
        return Err(format!(
            "Invalid minimizer_window: {}. It is how many consecutive k-mers share one minimizer: use 1 (every k-mer) to 256; 12 is the default",
            opt.minimizer_window
        ));
    }
    if opt.host_minimizer_window == 0 || opt.host_minimizer_window > 1024 {
        return Err(format!(
            "Invalid host_minimizer_window: {}. Use 1 to 1024; 31 is the default (a human genome then takes ~0.7 GB)",
            opt.host_minimizer_window
        ));
    }
    if opt.genome_window_bp < 50 || opt.genome_window_bp > 10_000_000 {
        return Err(format!(
            "Invalid genome_window_bp: {}. It is the size of the genome windows used for the breadth of coverage: use 50 to 10,000,000 bases (1,000 by default)",
            opt.genome_window_bp
        ));
    }
    if opt.memory_mb < 64 {
        return Err(format!("Invalid memory_mb: {}. Give at least 64 MB (1,024 by default) for building the host set", opt.memory_mb));
    }
    Ok(())
}

// ---------------------------------------------------------------------
//  BUILDING: the target index
// ---------------------------------------------------------------------

struct TargetBuilt {
    /// (hash, window in this genome), sorted by hash, one per hash: the
    /// window of its first occurrence.
    list: Vec<(u64, u32)>,
    sequences: u64,
    bases: u64,
    kmers: u64,
    window_bp: u32,
    windows: u32,
    first_header: Option<String>,
}

fn header_text(head: &[u8]) -> String {
    String::from_utf8_lossy(head).trim().to_string()
}

fn index_target(files: &[String], k: usize, w: usize, genome_window: u32) -> Result<TargetBuilt, String> {
    let mut pairs: Vec<(u64, u64)> = Vec::new();
    let mut t = TargetBuilt { list: Vec::new(), sequences: 0, bases: 0, kmers: 0, window_bp: 0, windows: 0, first_header: None };
    let mut offset = 0u64;
    for f in files {
        let mut rd = SeqReader::open(f)?;
        let mut rec = SeqRecord::default();
        while rd.next(&mut rec)? {
            if t.first_header.is_none() {
                t.first_header = Some(header_text(&rec.head));
            }
            t.sequences += 1;
            t.bases += rec.seq.len() as u64;
            t.kmers += count_valid_kmers(&rec.seq, k);
            let off = offset;
            sequence_minimizers(&rec.seq, k, w, |h, p| Some((index_key(h), off + p as u64)), |v| -> Result<(), String> {
                pairs.extend(v);
                Ok(())
            })?;
            offset += rec.seq.len() as u64;
        }
    }
    if t.bases == 0 {
        return Ok(t);
    }
    let wbp = u64::from(genome_window).min(MIN_WINDOW_BP.max(t.bases.div_ceil(SMALL_GENOME_WINDOWS)));
    t.window_bp = wbp as u32;
    t.windows = t.bases.div_ceil(wbp) as u32;
    // (hash, position): the first occurrence of each hash is its lowest
    // position.
    multi_merge_sort_stable(&mut pairs);
    let mut list: Vec<(u64, u32)> = Vec::with_capacity(pairs.len());
    let mut prev: Option<u64> = None;
    for (h, p) in pairs {
        if prev != Some(h) {
            list.push((h, (p / wbp) as u32));
            prev = Some(h);
        }
    }
    list.shrink_to_fit();
    t.list = list;
    Ok(t)
}

/// K-way merge of the sorted lists of the references into one sorted
/// array of keys with their occurrences.
fn merge_targets(lists: &[(usize, Vec<(u64, u32)>)], first_win: &[u32]) -> Result<(Vec<u64>, Vec<u32>, Vec<u32>, u64), String> {
    let total: usize = lists.iter().map(|(_, l)| l.len()).sum();
    let mut keys = Vec::with_capacity(total);
    let mut occ = Vec::with_capacity(total);
    let mut multi: Vec<u32> = Vec::new();
    let mut shared = 0u64;
    let mut pos = vec![0usize; lists.len()];
    let mut heap: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::with_capacity(lists.len());
    for (j, (_, l)) in lists.iter().enumerate() {
        if let Some(&(h, _)) = l.first() {
            heap.push(Reverse((h, j)));
        }
    }
    let mut cur: Vec<u32> = Vec::new();
    while let Some(Reverse((h, j))) = heap.pop() {
        cur.clear();
        let mut take = |j: usize, heap: &mut BinaryHeap<Reverse<(u64, usize)>>, cur: &mut Vec<u32>| {
            let (r, l) = &lists[j];
            cur.push(first_win[*r] + l[pos[j]].1);
            pos[j] += 1;
            if let Some(&(nh, _)) = l.get(pos[j]) {
                heap.push(Reverse((nh, j)));
            }
        };
        take(j, &mut heap, &mut cur);
        while let Some(&Reverse((h2, j2))) = heap.peek() {
            if h2 != h {
                break;
            }
            heap.pop();
            take(j2, &mut heap, &mut cur);
        }
        keys.push(h);
        if cur.len() == 1 {
            occ.push(cur[0]);
        } else {
            shared += 1;
            if multi.len() as u64 + cur.len() as u64 + 1 >= u64::from(MULTI) {
                return Err("The agent references are too large for one screen database (more than 2 billion shared minimizer entries): split them into several databases".to_string());
            }
            occ.push(MULTI | multi.len() as u32);
            multi.push(cur.len() as u32);
            multi.extend_from_slice(&cur);
        }
    }
    Ok((keys, occ, multi, shared))
}

// ---------------------------------------------------------------------
//  BUILDING: the host set, in passes when it does not fit in memory
// ---------------------------------------------------------------------

enum PassErr {
    Overflow,
    Fail(String),
}

impl From<String> for PassErr {
    fn from(e: String) -> Self {
        PassErr::Fail(e)
    }
}

struct HostBuilt {
    index: HostIndex,
    passes: u32,
    stats: Vec<(u64, u64, u64, Option<String>)>,
}

/// Estimated bases in a FASTA file from its size (a compressed genome is
/// about a quarter of its bases).
fn estimated_bases(path: &str) -> u64 {
    let size = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    let c = path.to_ascii_lowercase();
    if c.ends_with(".gz") || c.ends_with(".bgz") || c.ends_with(".lz4") {
        size.saturating_mul(4)
    } else {
        size
    }
}

fn build_host(refs: &[&[String]], k: usize, hw: usize, memory_mb: u64) -> Result<HostBuilt, String> {
    build_host_budget(refs, k, hw, ((memory_mb as usize) << 20) / 8, memory_mb)
}

/// The host set with a budget of `budget` 8-byte entries per pass.
fn build_host_budget(refs: &[&[String]], k: usize, hw: usize, budget: usize, memory_mb: u64) -> Result<HostBuilt, String> {
    if refs.is_empty() {
        return Ok(HostBuilt { index: HostIndex::empty(), passes: 0, stats: Vec::new() });
    }
    let est_bases: u64 = refs.iter().flat_map(|f| f.iter()).map(|p| estimated_bases(p)).sum();
    let est_min = (est_bases as f64 * 2.0 / (hw as f64 + 1.0)).max(1.0);
    let bits = ((est_min / 8.0).log2().ceil() as i64).clamp(8, 24) as u32;
    let budget = budget.max(1024);
    let mut passes = ((est_min / (0.6 * budget as f64)).ceil() as u64).max(1);
    loop {
        // The buffer of a pass is reserved once (a growing vector doubles,
        // and would pass the budget): the expected share of the pass, at
        // most the budget, plus room for one block of minimizers.
        let block = CHUNK * CHUNKS_AT_ONCE * 2 / (hw + 1) / passes as usize + 4096;
        let reserve = ((est_min * 1.2 / passes as f64) as usize).min(budget) + block;
        match host_passes(refs, k, hw, bits, passes, budget, reserve) {
            Ok(b) => return Ok(b),
            Err(PassErr::Overflow) if passes < (1u64 << bits) => passes *= 2,
            Err(PassErr::Overflow) => {
                return Err(format!(
                    "The host set does not fit in memory_mb = {memory_mb} even in {passes} passes: give more memory, or a larger host_minimizer_window"
                ))
            }
            Err(PassErr::Fail(e)) => return Err(e),
        }
    }
}

/// Builds only the host set of `files` with at most `budget` 8-byte
/// entries per pass (a testing aid: a small budget forces several passes).
/// Returns the set and the number of passes.
#[doc(hidden)]
pub fn host_set_with_budget(files: &[String], k: u8, host_window: u16, budget: usize) -> Result<(HostIndex, u32), String> {
    for f in files {
        check_input(f)?;
    }
    let b = build_host_budget(&[files], usize::from(k), usize::from(host_window), budget, (budget as u64 * 8) >> 20)?;
    Ok((b.index, b.passes))
}

fn host_passes(refs: &[&[String]], k: usize, hw: usize, bits: u32, passes: u64, budget: usize, reserve: usize) -> Result<HostBuilt, PassErr> {
    let nb = 1u64 << bits;
    let mut dir = vec![0u64; nb as usize + 1];
    let mut fp: Vec<u32> = Vec::new();
    let mut stats: Vec<(u64, u64, u64, Option<String>)> = vec![(0, 0, 0, None); refs.len()];
    let shift = 64 - bits;
    for p in 0..passes {
        let (lo, hi) = (nb * p / passes, nb * (p + 1) / passes);
        let mut buf: Vec<u64> = Vec::with_capacity(reserve);
        for (ri, files) in refs.iter().enumerate() {
            for f in files.iter() {
                let mut rd = SeqReader::open(f)?;
                let mut rec = SeqRecord::default();
                while rd.next(&mut rec)? {
                    if p == 0 {
                        let s = &mut stats[ri];
                        if s.3.is_none() {
                            s.3 = Some(header_text(&rec.head));
                        }
                        s.0 += 1;
                        s.1 += rec.seq.len() as u64;
                        s.2 += count_valid_kmers(&rec.seq, k);
                    }
                    sequence_minimizers(
                        &rec.seq,
                        k,
                        hw,
                        |h, _| {
                            let key = index_key(h);
                            let b = key >> shift;
                            if b >= lo && b < hi {
                                Some(key >> (32 - bits))
                            } else {
                                None
                            }
                        },
                        |v| -> Result<(), PassErr> {
                            buf.extend(v);
                            if buf.len() >= budget {
                                multi_merge_sort(&mut buf);
                                buf.dedup();
                                if buf.len() as f64 > 0.7 * budget as f64 {
                                    return Err(PassErr::Overflow);
                                }
                            }
                            Ok(())
                        },
                    )?;
                }
            }
        }
        // Integers: the order of equal values cannot be seen, and the
        // unstable form of the engine sorts in place (no second buffer).
        multi_merge_sort(&mut buf);
        buf.dedup();
        fp.reserve_exact(buf.len());
        for e in buf {
            dir[(e >> 32) as usize + 1] += 1;
            fp.push(e as u32);
        }
    }
    for i in 0..nb as usize {
        dir[i + 1] += dir[i];
    }
    Ok(HostBuilt { index: HostIndex { bits, dir, fp }, passes: passes as u32, stats })
}

// ---------------------------------------------------------------------
//  BUILDING: the whole database
// ---------------------------------------------------------------------

struct HashWriter {
    w: BufWriter<File>,
    h: Sha256,
    path: String,
    bytes: u64,
}

impl HashWriter {
    fn put(&mut self, b: &[u8]) -> Result<(), String> {
        self.h.update(b);
        self.bytes += b.len() as u64;
        self.w.write_all(b).map_err(|e| write_error(&self.path, e))
    }
    fn u8(&mut self, v: u8) -> Result<(), String> {
        self.put(&[v])
    }
    fn u16(&mut self, v: u16) -> Result<(), String> {
        self.put(&v.to_le_bytes())
    }
    fn u32(&mut self, v: u32) -> Result<(), String> {
        self.put(&v.to_le_bytes())
    }
    fn u64(&mut self, v: u64) -> Result<(), String> {
        self.put(&v.to_le_bytes())
    }
    fn text(&mut self, s: &str) -> Result<(), String> {
        self.u32(s.len() as u32)?;
        self.put(s.as_bytes())
    }
    fn u64s(&mut self, v: &[u64]) -> Result<(), String> {
        let mut b = Vec::with_capacity(8 * 8192);
        for c in v.chunks(8192) {
            b.clear();
            for x in c {
                b.extend_from_slice(&x.to_le_bytes());
            }
            self.put(&b)?;
        }
        Ok(())
    }
    fn u32s(&mut self, v: &[u32]) -> Result<(), String> {
        let mut b = Vec::with_capacity(4 * 8192);
        for c in v.chunks(8192) {
            b.clear();
            for x in c {
                b.extend_from_slice(&x.to_le_bytes());
            }
            self.put(&b)?;
        }
        Ok(())
    }
}

/// Indexes the references into a screen database file. Agents (virus,
/// bacterium, mycoplasma, fungus, other) and controls (spike) go to the
/// target index, with their genome windows; hosts go to the host set. Each
/// reference file is read once (the host files once per pass), and its
/// SHA-256 is computed at the same time and kept in the database.
pub fn build_screen_database(refs: &[ScreenReference], output: &str, opt: &ScreenDbOptions) -> Result<ScreenDbSummary, String> {
    let t0 = Instant::now();
    check_db_options(opt)?;
    if refs.is_empty() {
        return Err("No reference to index: give at least one FASTA file of an agent to screen for (a virus, a bacterium, a mycoplasma or a fungus), and the host genome".to_string());
    }
    let mut inputs: Vec<&str> = Vec::new();
    for (i, r) in refs.iter().enumerate() {
        if r.files.is_empty() {
            return Err(format!(
                "Reference {} ({}) has no file: give at least one FASTA file for each reference",
                i + 1,
                if r.name.trim().is_empty() { "unnamed" } else { r.name.trim() }
            ));
        }
        for f in &r.files {
            check_input(f)?;
            inputs.push(f);
        }
    }
    if refs.iter().all(|r| r.kind == AgentKind::Host) {
        return Err("Every reference is a host: give also the genomes of the agents to screen for (kind virus, bacterium, mycoplasma, fungus or other)".to_string());
    }
    not_an_input(output, &inputs)?;
    let k = usize::from(opt.k);
    let w = usize::from(opt.minimizer_window);
    let hw = usize::from(opt.host_minimizer_window);
    let mut all_files: Vec<String> = Vec::new();
    for f in &inputs {
        if !all_files.iter().any(|g| g == f) {
            all_files.push(f.to_string());
        }
    }

    let mut infos: Vec<ScreenRefInfo> = refs
        .iter()
        .map(|r| ScreenRefInfo { name: r.name.trim().to_string(), kind: r.kind, files: r.files.clone(), ..Default::default() })
        .collect();
    let mut warnings: Vec<String> = Vec::new();

    let (hashed, built) = std::thread::scope(|s| {
        let hasher = s.spawn(|| hash_files(&all_files));
        let built = (|| -> Result<(Vec<(usize, Vec<(u64, u32)>)>, HostBuilt), String> {
            let mut lists: Vec<(usize, Vec<(u64, u32)>)> = Vec::new();
            for (i, r) in refs.iter().enumerate() {
                if r.kind == AgentKind::Host {
                    continue;
                }
                let t = index_target(&r.files, k, w, opt.genome_window_bp)?;
                let info = &mut infos[i];
                if info.name.is_empty() {
                    info.name = t.first_header.clone().unwrap_or_else(|| r.files[0].clone());
                }
                info.sequences = t.sequences;
                info.bases = t.bases;
                info.kmers = t.kmers;
                info.minimizers = t.list.len() as u64;
                info.window_bp = t.window_bp;
                info.windows = t.windows;
                lists.push((i, t.list));
            }
            let host_refs: Vec<(usize, &[String])> =
                refs.iter().enumerate().filter(|(_, r)| r.kind == AgentKind::Host).map(|(i, r)| (i, r.files.as_slice())).collect();
            let files: Vec<&[String]> = host_refs.iter().map(|(_, f)| *f).collect();
            let hb = build_host(&files, k, hw, opt.memory_mb)?;
            for ((i, r), st) in host_refs.iter().zip(&hb.stats) {
                let info = &mut infos[*i];
                if info.name.is_empty() {
                    info.name = st.3.clone().unwrap_or_else(|| r[0].clone());
                }
                info.sequences = st.0;
                info.bases = st.1;
                info.kmers = st.2;
            }
            Ok((lists, hb))
        })();
        (hasher.join(), built)
    });
    let hashed = hashed.map_err(|_| "Internal error while computing the SHA-256 of the references".to_string())??;
    let (lists, hb) = built?;
    for info in infos.iter_mut() {
        for f in &info.files {
            let (n, h) = hashed.get(f).cloned().unwrap_or_default();
            info.file_bytes.push(n);
            info.sha256.push(h);
        }
    }
    for (i, a) in infos.iter().enumerate() {
        if let Some(b) = infos[..i].iter().find(|b| b.name == a.name) {
            return Err(format!(
                "Two references have the same name, \"{}\" ({} and {}): give each reference its own name",
                a.name,
                b.files.join(", "),
                a.files.join(", ")
            ));
        }
    }
    // Global windows, in the order of the references.
    let mut first_win = vec![0u32; refs.len()];
    let mut total: u64 = 0;
    for (i, info) in infos.iter_mut().enumerate() {
        if info.kind == AgentKind::Host {
            continue;
        }
        first_win[i] = total as u32;
        info.first_window = total as u32;
        total += u64::from(info.windows);
        if total >= u64::from(MULTI) {
            return Err("The agent references have too many genome windows for one database (more than 2 billion): use a larger genome_window_bp, or split the references".to_string());
        }
        if info.kmers == 0 {
            warnings.push(format!("{}: no k-mer of {k} A/C/G/T bases in a row: this reference can never be found", info.name));
        } else if info.minimizers < 50 {
            warnings.push(format!(
                "{}: only {} minimizers ({} bases): a genome this small gives little evidence; a detection needs reads spread over it",
                info.name, info.minimizers, info.bases
            ));
        }
    }
    let (keys, occ, multi, shared) = merge_targets(&lists, &first_win)?;
    drop(lists);
    if hb.index.is_empty() {
        warnings.push("No host reference: host reads cannot be set aside, and every read of the material is screened against the agents. Add the host genome (kind \"host\")".to_string());
    }
    for info in infos.iter_mut().filter(|i| i.kind == AgentKind::Host) {
        info.minimizers = 0;
    }

    let mut cleanup = Cleanup(vec![output.to_string()]);
    let file = File::create(output).map_err(|e| {
        format!("Error creating the screen database {output}: {e}. Check that the output folder exists and that you can write to it")
    })?;
    let mut hwr = HashWriter { w: BufWriter::with_capacity(1 << 20, file), h: Sha256::new(), path: output.to_string(), bytes: 0 };
    hwr.put(DB_MAGIC)?;
    hwr.u32(DB_VERSION)?;
    hwr.u8(opt.k)?;
    hwr.u8(0)?;
    hwr.u16(opt.minimizer_window)?;
    hwr.u16(opt.host_minimizer_window)?;
    hwr.u16(0)?;
    hwr.u32(opt.genome_window_bp)?;
    hwr.u32(hb.index.bits)?;
    hwr.u32(hb.passes)?;
    hwr.u32(infos.len() as u32)?;
    for r in &infos {
        hwr.u8(r.kind.code())?;
        hwr.text(&r.name)?;
        hwr.text(&r.files.join("\n"))?;
        hwr.text(&r.sha256.join("\n"))?;
        hwr.text(&r.file_bytes.iter().map(|b| b.to_string()).collect::<Vec<_>>().join("\n"))?;
        for v in [r.sequences, r.bases, r.kmers, r.minimizers] {
            hwr.u64(v)?;
        }
        hwr.u32(r.window_bp)?;
        hwr.u32(r.first_window)?;
        hwr.u32(r.windows)?;
    }
    hwr.u64(keys.len() as u64)?;
    hwr.u64s(&keys)?;
    hwr.u32s(&occ)?;
    hwr.u64(multi.len() as u64)?;
    hwr.u32s(&multi)?;
    hwr.u64(hb.index.fp.len() as u64)?;
    hwr.u64s(&hb.index.dir)?;
    hwr.u32s(&hb.index.fp)?;
    let body = hwr.h.clone().finalize();
    hwr.put(DB_END)?;
    hwr.put(&body)?;
    hwr.w.flush().map_err(|e| write_error(output, e))?;
    let file_sha = hex(&hwr.h.clone().finalize());
    let file_bytes = hwr.bytes;
    drop(hwr);
    cleanup.disarm();
    let tbits = target_bits(keys.len());
    Ok(ScreenDbSummary {
        path: output.to_string(),
        sha256: file_sha,
        file_bytes,
        k: opt.k,
        minimizer_window: opt.minimizer_window,
        host_minimizer_window: opt.host_minimizer_window,
        genome_window_bp: opt.genome_window_bp,
        references: infos,
        target_minimizers: keys.len() as u64,
        shared_minimizers: shared,
        host_minimizers: hb.index.fp.len() as u64,
        host_bucket_bits: hb.index.bits,
        host_passes: hb.passes,
        target_index_bytes: target_bytes(keys.len(), multi.len(), tbits, total as usize),
        host_index_bytes: hb.index.bytes(),
        seconds: t0.elapsed().as_secs_f64(),
        warnings,
    })
}

fn target_bits(n: usize) -> u32 {
    (((n.max(1) as f64) / 4.0).log2().ceil() as i64).clamp(4, 24) as u32
}

fn target_bytes(n_keys: usize, n_multi: usize, tbits: u32, windows: usize) -> u64 {
    (n_keys * 12 + n_multi * 4 + ((1usize << tbits) + 1) * 4 + windows * 4) as u64
}

// ---------------------------------------------------------------------
//  READING THE DATABASE BACK
// ---------------------------------------------------------------------

struct HashReader {
    r: BufReader<File>,
    h: Sha256,
    path: String,
    pos: u64,
    len: u64,
}

impl HashReader {
    fn damaged(&self, what: &str) -> String {
        format!(
            "The screen database {} is damaged or incomplete ({what}): it was probably cut during a copy. Copy it again, or build it again with build_screen_database",
            self.path
        )
    }
    fn take(&mut self, buf: &mut [u8]) -> Result<(), String> {
        self.r.read_exact(buf).map_err(|e| {
            if e.kind() == ErrorKind::UnexpectedEof {
                self.damaged("it ends too early")
            } else {
                format!("Error reading the screen database {}: {e}", self.path)
            }
        })?;
        self.h.update(buf);
        self.pos += buf.len() as u64;
        Ok(())
    }
    fn u8(&mut self) -> Result<u8, String> {
        let mut b = [0u8; 1];
        self.take(&mut b)?;
        Ok(b[0])
    }
    fn u16(&mut self) -> Result<u16, String> {
        let mut b = [0u8; 2];
        self.take(&mut b)?;
        Ok(u16::from_le_bytes(b))
    }
    fn u32(&mut self) -> Result<u32, String> {
        let mut b = [0u8; 4];
        self.take(&mut b)?;
        Ok(u32::from_le_bytes(b))
    }
    fn u64(&mut self) -> Result<u64, String> {
        let mut b = [0u8; 8];
        self.take(&mut b)?;
        Ok(u64::from_le_bytes(b))
    }
    fn room(&self, n: u64, size: u64) -> Result<usize, String> {
        if n.saturating_mul(size) > self.len.saturating_sub(self.pos) {
            return Err(self.damaged("a table is longer than the file"));
        }
        Ok(n as usize)
    }
    fn text(&mut self) -> Result<String, String> {
        let n = self.u32()?;
        let n = self.room(u64::from(n), 1)?;
        let mut b = vec![0u8; n];
        self.take(&mut b)?;
        Ok(String::from_utf8_lossy(&b).to_string())
    }
    fn u64s(&mut self, n: u64) -> Result<Vec<u64>, String> {
        let n = self.room(n, 8)?;
        let mut v = Vec::with_capacity(n);
        let mut b = vec![0u8; 8 * 65536];
        while v.len() < n {
            let m = (n - v.len()).min(65536);
            self.take(&mut b[..8 * m])?;
            v.extend(b[..8 * m].chunks_exact(8).map(|c| u64::from_le_bytes(c.try_into().unwrap())));
        }
        Ok(v)
    }
    fn u32s(&mut self, n: u64) -> Result<Vec<u32>, String> {
        let n = self.room(n, 4)?;
        let mut v = Vec::with_capacity(n);
        let mut b = vec![0u8; 4 * 65536];
        while v.len() < n {
            let m = (n - v.len()).min(65536);
            self.take(&mut b[..4 * m])?;
            v.extend(b[..4 * m].chunks_exact(4).map(|c| u32::from_le_bytes(c.try_into().unwrap())));
        }
        Ok(v)
    }
}

fn split_lines(s: &str) -> Vec<String> {
    if s.is_empty() {
        Vec::new()
    } else {
        s.split('\n').map(str::to_string).collect()
    }
}

/// Reads a screen database, checking its SHA-256 (the SHA-256 of the whole
/// file is computed in the same pass and goes to the reports).
pub fn load_screen_database(path: &str) -> Result<ScreenDb, String> {
    check_input(path)?;
    let file = File::open(path).map_err(|e| format!("Error opening the screen database {path}: {e}"))?;
    let len = file.metadata().map(|m| m.len()).unwrap_or(0);
    let mut rd = HashReader { r: BufReader::with_capacity(4 << 20, file), h: Sha256::new(), path: path.to_string(), pos: 0, len };
    let mut magic = [0u8; 8];
    if len < 8 || rd.take(&mut magic).is_err() || &magic != DB_MAGIC {
        if &magic[..7] == b"MMGSKDB" {
            return Err(format!(
                "{path} is a sketch database of sketch_references (for gather and classify_reads), not a screen database: build one with build_screen_database"
            ));
        }
        return Err(format!(
            "{path} is not an MMGEN screen database (it does not start with \"MMGSCRN1\"). Build one with build_screen_database from the FASTA files of the agents and of the host"
        ));
    }
    let version = rd.u32()?;
    if version != DB_VERSION {
        return Err(format!("{path} is a screen database of version {version}, which this program does not know: update MMGEN, or build the database again"));
    }
    let k = rd.u8()?;
    let _ = rd.u8()?;
    let w = rd.u16()?;
    let hw = rd.u16()?;
    let _ = rd.u16()?;
    let genome_window_bp = rd.u32()?;
    let bits = rd.u32()?;
    let passes = rd.u32()?;
    let n_refs = rd.u32()?;
    if !(15..=31).contains(&k) || k % 2 == 0 || w == 0 || hw == 0 || !(8..=24).contains(&bits) {
        return Err(rd.damaged("its header holds values this program never writes"));
    }
    let n_refs = rd.room(u64::from(n_refs), 40)?;
    let mut refs = Vec::with_capacity(n_refs);
    for _ in 0..n_refs {
        let code = rd.u8()?;
        let kind = AgentKind::from_code(code).ok_or_else(|| rd.damaged("a reference has an unknown kind"))?;
        let name = rd.text()?;
        let files = split_lines(&rd.text()?);
        let sha256 = split_lines(&rd.text()?);
        let file_bytes: Vec<u64> = split_lines(&rd.text()?).iter().map(|s| s.parse().unwrap_or(0)).collect();
        let sequences = rd.u64()?;
        let bases = rd.u64()?;
        let kmers = rd.u64()?;
        let minimizers = rd.u64()?;
        let window_bp = rd.u32()?;
        let first_window = rd.u32()?;
        let windows = rd.u32()?;
        refs.push(ScreenRefInfo { name, kind, files, sha256, file_bytes, sequences, bases, kmers, minimizers, window_bp, first_window, windows });
    }
    let n_keys = rd.u64()?;
    let keys = rd.u64s(n_keys)?;
    let occ = rd.u32s(n_keys)?;
    let n_multi = rd.u64()?;
    let multi = rd.u32s(n_multi)?;
    let n_host = rd.u64()?;
    let dir = rd.u64s((1u64 << bits) + 1)?;
    let fp = rd.u32s(n_host)?;
    let body = rd.h.clone().finalize();
    let mut end = [0u8; 8];
    rd.take(&mut end)?;
    let mut stored = [0u8; 32];
    rd.take(&mut stored)?;
    if &end != DB_END {
        return Err(rd.damaged("its end mark is missing"));
    }
    if stored != body {
        return Err(rd.damaged("its SHA-256 does not match the content"));
    }
    let mut extra = [0u8; 1];
    if rd.r.read(&mut extra).map_err(|e| format!("Error reading the screen database {path}: {e}"))? != 0 {
        return Err(rd.damaged("there are bytes after its end mark"));
    }
    let file_sha = hex(&rd.h.finalize());

    // Consistency: windows, occurrences, the host directory.
    let mut total = 0u64;
    for r in refs.iter().filter(|r| r.kind != AgentKind::Host) {
        if u64::from(r.first_window) != total {
            return Err(format!(
                "The screen database {path} is damaged or incomplete (the genome windows do not follow each other): build it again"
            ));
        }
        total += u64::from(r.windows);
    }
    let total = total as usize;
    let mut win_ref = vec![0u32; total];
    for (i, r) in refs.iter().enumerate() {
        if r.kind != AgentKind::Host {
            for g in r.first_window as usize..(r.first_window + r.windows) as usize {
                win_ref[g] = i as u32;
            }
        }
    }
    let bad = || format!("The screen database {path} is damaged or incomplete (an occurrence points outside its tables): build it again");
    for &o in &occ {
        if o & MULTI == 0 {
            if o as usize >= total {
                return Err(bad());
            }
        } else {
            let off = (o & !MULTI) as usize;
            let n = *multi.get(off).ok_or_else(bad)? as usize;
            if n < 2 || off + n >= multi.len() || multi[off + 1..=off + n].iter().any(|&g| g as usize >= total) {
                return Err(bad());
            }
        }
    }
    if keys.windows(2).any(|p| p[0] >= p[1]) || dir.windows(2).any(|p| p[0] > p[1]) || dir.last().copied() != Some(fp.len() as u64) {
        return Err(format!("The screen database {path} is damaged or incomplete (a table is out of order): build it again"));
    }
    let tbits = target_bits(keys.len());
    let mut tdir = vec![0u32; (1usize << tbits) + 1];
    for &h in &keys {
        tdir[(h >> (64 - tbits)) as usize + 1] += 1;
    }
    for i in 0..(1usize << tbits) {
        tdir[i + 1] += tdir[i];
    }
    let shared = occ.iter().filter(|&&o| o & MULTI != 0).count() as u64;
    let host = HostIndex { bits, dir, fp };
    let summary = ScreenDbSummary {
        path: path.to_string(),
        sha256: file_sha,
        file_bytes: len,
        k,
        minimizer_window: w,
        host_minimizer_window: hw,
        genome_window_bp,
        references: refs,
        target_minimizers: keys.len() as u64,
        shared_minimizers: shared,
        host_minimizers: host.fp.len() as u64,
        host_bucket_bits: bits,
        host_passes: passes,
        target_index_bytes: target_bytes(keys.len(), multi.len(), tbits, total),
        host_index_bytes: host.bytes(),
        seconds: 0.0,
        warnings: Vec::new(),
    };
    Ok(ScreenDb {
        summary,
        k: usize::from(k),
        w: usize::from(w),
        hw: usize::from(hw),
        keys,
        occ,
        multi,
        tbits,
        tdir,
        win_ref,
        total_windows: total,
        host,
    })
}

/// What a screen database holds (names, kinds, sizes, SHA-256).
pub fn describe_screen_database(path: &str) -> Result<ScreenDbSummary, String> {
    Ok(load_screen_database(path)?.summary)
}

// ---------------------------------------------------------------------
//  CLASSIFYING ONE READ OR PAIR
// ---------------------------------------------------------------------

/// Bits set from several threads.
struct Bits(Vec<AtomicU64>);

impl Bits {
    fn new(n: usize) -> Bits {
        Bits((0..n.div_ceil(64)).map(|_| AtomicU64::new(0)).collect())
    }
    #[inline]
    fn set(&self, i: usize) {
        let (w, b) = (i >> 6, 1u64 << (i & 63));
        // Most bits are set already: a read is cheaper than a write.
        if self.0[w].load(Ordering::Relaxed) & b == 0 {
            self.0[w].fetch_or(b, Ordering::Relaxed);
        }
    }
    #[inline]
    fn get(&self, i: usize) -> bool {
        self.0[i >> 6].load(Ordering::Relaxed) >> (i & 63) & 1 == 1
    }
    fn count(&self, s: usize, e: usize) -> u32 {
        (s..e).filter(|&i| self.get(i)).count() as u32
    }
}

/// The read-level parameters.
#[derive(Clone, Copy, Debug)]
struct ReadParams {
    min_hits: u32,
    min_fraction: f64,
    unique_margin: f64,
    host_min_hits: u32,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Outcome {
    /// No valid k-mer (too short, or N).
    NoKmers,
    Unclassified,
    Host,
    /// Enough hits but a small fraction of the read: a distant relative.
    Weak(u32),
    Unique(u32),
    Tie(Vec<u32>),
}

pub struct FragCall {
    pub outcome: Outcome,
    /// Target minimizers shared with the best reference.
    pub hits: u32,
    /// Distinct target minimizers of the read (pair).
    pub minimizers: u32,
    pub host_hits: u32,
    pub host_minimizers: u32,
    /// Distinct genome windows hit (unique reads).
    pub windows: u32,
}

#[derive(Default)]
struct Scratch {
    dq: VecDeque<(u64, usize, bool)>,
    dqh: VecDeque<(u64, usize, bool)>,
    t: Vec<u64>,
    h: Vec<u64>,
    occ: Vec<(u32, u32, usize)>,
    refs: Vec<(u32, u32)>,
    wins: Vec<u32>,
}

struct Marks {
    slots: Bits,
    win_unique: Bits,
    win_any: Bits,
}

/// The minimizers of a read for the agent scheme (window `w`, into `t`)
/// and the host scheme (window `hw`, into `h`; 0 = none) in one pass, the
/// low-complexity k-mers dropped. The same minimizers as `minimizers` +
/// `low_complexity`, with the k-mer hashed once and the distinct
/// trinucleotides counted as the k-mer slides.
#[allow(clippy::too_many_arguments)]
fn read_minimizers(
    seq: &[u8],
    k: usize,
    w: usize,
    hw: usize,
    dq: &mut VecDeque<(u64, usize, bool)>,
    dqh: &mut VecDeque<(u64, usize, bool)>,
    t: &mut Vec<u64>,
    h: &mut Vec<u64>,
) {
    dq.clear();
    dqh.clear();
    if seq.len() < k {
        return;
    }
    let mask = (1u64 << (2 * k)) - 1;
    let shift = 2 * (k - 1);
    let (mut fw, mut rv) = (0u64, 0u64);
    let (mut run, mut nk) = (0usize, 0usize);
    let (mut last_t, mut last_h) = (usize::MAX, usize::MAX);
    let mut tri = [0u8; 64];
    let mut distinct = 0u32;
    // One scheme: the deque, its window, its last output.
    fn step(dq: &mut VecDeque<(u64, usize, bool)>, w: usize, nk: usize, x: (u64, usize, bool), last: &mut usize, out: &mut Vec<u64>) {
        while let Some(&(bh, _, _)) = dq.back() {
            if bh > x.0 {
                dq.pop_back();
            } else {
                break;
            }
        }
        dq.push_back(x);
        if nk >= w {
            let start = x.1 + 1 - w;
            while let Some(&(_, p, _)) = dq.front() {
                if p < start {
                    dq.pop_front();
                } else {
                    break;
                }
            }
            if let Some(&(mh, mp, lc)) = dq.front() {
                if mp != *last {
                    *last = mp;
                    if !lc {
                        out.push(mh);
                    }
                }
            }
        }
    }
    for (i, &b) in seq.iter().enumerate() {
        let c = CODE[usize::from(b)];
        if c > 3 {
            if run > 0 {
                tri = [0u8; 64];
                distinct = 0;
            }
            run = 0;
            nk = 0;
            dq.clear();
            dqh.clear();
            continue;
        }
        let c = u64::from(c);
        fw = ((fw << 2) | c) & mask;
        rv = (rv >> 2) | ((3 - c) << shift);
        run += 1;
        if run >= 3 {
            let code = (fw & 63) as usize;
            if tri[code] == 0 {
                distinct += 1;
            }
            tri[code] += 1;
            if run > k {
                // The trinucleotide that ends at i + 2 - k left the k-mer.
                let j = i + 2 - k;
                let old = (usize::from(CODE[usize::from(seq[j - 2])]) << 4)
                    | (usize::from(CODE[usize::from(seq[j - 1])]) << 2)
                    | usize::from(CODE[usize::from(seq[j])]);
                tri[old] -= 1;
                if tri[old] == 0 {
                    distinct -= 1;
                }
            }
        }
        if run < k {
            continue;
        }
        let x = (kmer_hash(fw.min(rv)), i + 1 - k, distinct < LOW_COMPLEXITY_TRIMERS);
        nk += 1;
        step(dq, w, nk, x, &mut last_t, t);
        if hw > 0 {
            step(dqh, hw, nk, x, &mut last_h, h);
        }
    }
}

/// Testing aid: the agent and host minimizers of one read, as the
/// classification computes them.
pub fn read_minimizer_sets(db: &ScreenDb, seq: &[u8]) -> (Vec<u64>, Vec<u64>) {
    let (mut dq, mut dqh) = (VecDeque::new(), VecDeque::new());
    let (mut t, mut h) = (Vec::new(), Vec::new());
    read_minimizers(seq, db.k, db.w, if db.host.is_empty() { 0 } else { db.hw }, &mut dq, &mut dqh, &mut t, &mut h);
    (t, h)
}

fn classify_fragment(db: &ScreenDb, a: &[u8], b: Option<&[u8]>, p: &ReadParams, sc: &mut Scratch, marks: Option<&Marks>) -> FragCall {
    sc.t.clear();
    sc.h.clear();
    let hw = if db.host.is_empty() { 0 } else { db.hw };
    for m in std::iter::once(a).chain(b) {
        read_minimizers(m, db.k, db.w, hw, &mut sc.dq, &mut sc.dqh, &mut sc.t, &mut sc.h);
    }
    for x in sc.t.iter_mut().chain(sc.h.iter_mut()) {
        *x = index_key(*x);
    }
    multi_merge_sort(&mut sc.t);
    sc.t.dedup();
    multi_merge_sort(&mut sc.h);
    sc.h.dedup();
    let nt = sc.t.len() as u32;
    let nh = sc.h.len() as u32;
    let mut call = FragCall { outcome: Outcome::NoKmers, hits: 0, minimizers: nt, host_hits: 0, host_minimizers: nh, windows: 0 };
    if nt == 0 && nh == 0 {
        return call;
    }
    call.host_hits = sc.h.iter().filter(|&&h| db.host.contains(h)).count() as u32;
    // A read is host when its host fraction is at least its best target
    // fraction, which is at most (found so far + still to look up) / nt: a
    // read that looks like host stops the lookups as soon as that bound
    // falls to its host fraction (the same call, sooner).
    let hfrac_early = if nh > 0 { f64::from(call.host_hits) / f64::from(nh) } else { 0.0 };
    let host_possible = call.host_hits > 0 && call.host_hits >= p.host_min_hits;
    sc.occ.clear();
    let mut found = 0usize;
    for (j, &h) in sc.t.iter().enumerate() {
        if host_possible && (found + sc.t.len() - j) as f64 / f64::from(nt) <= hfrac_early {
            call.outcome = Outcome::Host;
            return call;
        }
        if let Some(i) = db.lookup(h) {
            found += 1;
            let occ = &mut sc.occ;
            db.occurrences(i, |r, g, s| occ.push((r, g, s)));
        }
    }
    sc.refs.clear();
    for &(r, _, _) in &sc.occ {
        match sc.refs.iter_mut().find(|x| x.0 == r) {
            Some(x) => x.1 += 1,
            None => sc.refs.push((r, 1)),
        }
    }
    let (mut best_ref, mut best) = (u32::MAX, 0u32);
    for &(r, n) in &sc.refs {
        if n > best || (n == best && r < best_ref) {
            best = n;
            best_ref = r;
        }
    }
    call.hits = best;
    let tfrac = if nt > 0 { f64::from(best) / f64::from(nt) } else { 0.0 };
    let hfrac = if nh > 0 { f64::from(call.host_hits) / f64::from(nh) } else { 0.0 };
    if call.host_hits > 0 && call.host_hits >= p.host_min_hits && hfrac >= tfrac {
        call.outcome = Outcome::Host;
        return call;
    }
    if best == 0 || best < p.min_hits {
        call.outcome = Outcome::Unclassified;
        return call;
    }
    if tfrac < p.min_fraction {
        call.outcome = Outcome::Weak(best_ref);
        return call;
    }
    let cut = f64::from(best) * (1.0 - p.unique_margin);
    let mut set: Vec<u32> = sc.refs.iter().filter(|x| x.1 == best || f64::from(x.1) >= cut).map(|x| x.0).collect();
    set.sort_unstable();
    if set.len() == 1 {
        sc.wins.clear();
        for &(r, g, s) in &sc.occ {
            if r == best_ref {
                sc.wins.push(g);
                if let Some(m) = marks {
                    m.slots.set(s);
                    m.win_unique.set(g as usize);
                    m.win_any.set(g as usize);
                }
            }
        }
        multi_merge_sort(&mut sc.wins);
        sc.wins.dedup();
        call.windows = sc.wins.len() as u32;
        call.outcome = Outcome::Unique(best_ref);
    } else {
        if let Some(m) = marks {
            for &(r, g, s) in &sc.occ {
                if set.binary_search(&r).is_ok() {
                    m.slots.set(s);
                    m.win_any.set(g as usize);
                }
            }
        }
        call.outcome = Outcome::Tie(set);
    }
    call
}

/// Classifies one read or pair against a loaded database with the given
/// read-level options (a testing and teaching aid).
pub fn classify_one(db: &ScreenDb, read1: &[u8], read2: Option<&[u8]>, opt: &ScreenOptions) -> FragCall {
    let p = ReadParams { min_hits: opt.min_hits, min_fraction: opt.min_fraction, unique_margin: opt.unique_margin, host_min_hits: opt.host_min_hits };
    let mut sc = Scratch::default();
    classify_fragment(db, read1, read2, &p, &mut sc, None)
}

// ---------------------------------------------------------------------
//  READING THE SAMPLES: single or paired reads, one or several lanes
// ---------------------------------------------------------------------

/// One sample: its name and its read files. `r2` empty: single-end reads;
/// else one second-read file for each first-read file (lanes in order).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ScreenSample {
    pub name: String,
    pub r1: Vec<String>,
    pub r2: Vec<String>,
}

type BatchMsg = Result<(usize, Vec<SeqRecord>), String>;

fn read_files(paths: &[String], n: usize, tx: &SyncSender<BatchMsg>) {
    for (fi, p) in paths.iter().enumerate() {
        let mut rd = match SeqReader::open(p) {
            Ok(r) => r,
            Err(e) => {
                let _ = tx.send(Err(e));
                return;
            }
        };
        loop {
            let mut recs = Vec::with_capacity(n.min(4096));
            let mut end = false;
            while recs.len() < n {
                let mut rec = SeqRecord::default();
                match rd.next(&mut rec) {
                    Ok(true) => recs.push(rec),
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
            if !recs.is_empty() && tx.send(Ok((fi, recs))).is_err() {
                return;
            }
            if end {
                break;
            }
        }
    }
}

/// The read name: the first word of the header, without a trailing /1 or /2.
fn read_id(head: &[u8]) -> &[u8] {
    let first = head.split(|c| c.is_ascii_whitespace()).next().unwrap_or(head);
    if first.len() > 2 && first[first.len() - 2] == b'/' && matches!(first[first.len() - 1], b'1' | b'2') {
        &first[..first.len() - 2]
    } else {
        first
    }
}

fn shown(b: &[u8]) -> String {
    let s = String::from_utf8_lossy(&b[..b.len().min(60)]).to_string();
    if b.len() > 60 {
        s + "..."
    } else {
        s
    }
}

/// Calls `f(first reads, second reads)` with batches of the sample's reads.
fn for_each_fragments(
    s: &ScreenSample,
    n: usize,
    mut f: impl FnMut(&[SeqRecord], Option<&[SeqRecord]>) -> Result<(), String>,
) -> Result<(), String> {
    std::thread::scope(|sc| {
        let (tx1, rx1) = sync_channel::<BatchMsg>(2);
        let a = s.r1.clone();
        sc.spawn(move || read_files(&a, n, &tx1));
        if s.r2.is_empty() {
            for m in rx1.iter() {
                let (_, r1) = m?;
                f(&r1, None)?;
            }
            return Ok(());
        }
        let (tx2, rx2) = sync_channel::<BatchMsg>(2);
        let b = s.r2.clone();
        sc.spawn(move || read_files(&b, n, &tx2));
        let mut done = 0u64;
        loop {
            match (rx1.recv(), rx2.recv()) {
                (Err(_), Err(_)) => break,
                (Ok(x), Ok(y)) => {
                    let (f1, r1) = x?;
                    let (f2, r2) = y?;
                    if f1 != f2 || r1.len() != r2.len() {
                        let fi = f1.min(f2);
                        return Err(format!(
                            "The two files of the pair do not have the same number of reads: {} and {} differ near read {}. Give the two complete files of the same run (a file cut during the download has fewer reads)",
                            s.r1[fi],
                            s.r2[fi],
                            thousands(done + r1.len().min(r2.len()) as u64 + 1)
                        ));
                    }
                    for (i, (x, y)) in r1.iter().zip(r2.iter()).enumerate() {
                        if read_id(&x.head) != read_id(&y.head) {
                            return Err(format!(
                                "The reads are not paired: read {} is \"{}\" in {} (line {}) but \"{}\" in {} (line {}). The two files of a pair must list the same reads in the same order: give the two files of the same run, or re-pair them (for instance `repair.sh` of BBTools)",
                                thousands(done + i as u64 + 1),
                                shown(read_id(&x.head)),
                                s.r1[f1],
                                x.line,
                                shown(read_id(&y.head)),
                                s.r2[f2],
                                y.line
                            ));
                        }
                    }
                    f(&r1, Some(&r2))?;
                    done += r1.len() as u64;
                }
                (Ok(x), Err(_)) => {
                    let (fi, _) = x?;
                    return Err(format!(
                        "The two files of the pair do not have the same number of reads: {} ends after {} reads, but {} has more. Give the two complete files of the same run",
                        s.r2.get(fi).map_or("the second-read file", String::as_str),
                        thousands(done),
                        s.r1[fi]
                    ));
                }
                (Err(_), Ok(y)) => {
                    let (fi, _) = y?;
                    return Err(format!(
                        "The two files of the pair do not have the same number of reads: {} ends after {} reads, but {} has more. Give the two complete files of the same run",
                        s.r1.get(fi).map_or("the first-read file", String::as_str),
                        thousands(done),
                        s.r2[fi]
                    ));
                }
            }
        }
        Ok(())
    })
}

// ---------------------------------------------------------------------
//  THE SCREEN: options and results
// ---------------------------------------------------------------------

/// One spike-in: an organism of the database added to a sample at a known
/// amount. `sample` None: every sample that is not a negative control.
#[derive(Clone, Debug, PartialEq)]
pub struct SpikeIn {
    pub sample: Option<String>,
    pub organism: String,
    pub amount: f64,
    /// "reads" (reads or pairs added), "rpm" (reads per million), or any
    /// other unit ("copies", "copies/mL", "CFU"...): then the recovery is
    /// given as reads per unit.
    pub unit: String,
}

#[derive(Clone, Debug)]
pub struct ScreenOptions {
    /// Names of the samples that are negative controls (blanks).
    pub negative_controls: Vec<String>,
    pub spike_ins: Vec<SpikeIn>,
    /// A spike-in table (TSV: sample, organism, amount, unit), added to
    /// `spike_ins`.
    pub spike_table: Option<String>,
    /// A read (pair) is given to a reference when it shares at least this
    /// many minimizers with it (3)...
    pub min_hits: u32,
    /// ... and at least this fraction of its minimizers (0.2); below it the
    /// read is "weak" (a distant relative, not counted).
    pub min_fraction: f64,
    /// The best reference must beat the second by more than this fraction
    /// of its hits for the read to be unique (0.2); else it is tied.
    pub unique_margin: f64,
    /// A read is host when it shares at least this many host minimizers (1)
    /// and a larger fraction of them than of its target minimizers.
    pub host_min_hits: u32,
    /// Call thresholds, on the reads of an organism (its unique reads and
    /// the tied reads it wins): at least `min_reads` (3) in at least
    /// `min_windows` (2) genome windows, spread with a dispersion of at
    /// least `min_dispersion` (0.5)...
    pub min_reads: u64,
    pub min_windows: u32,
    pub min_dispersion: f64,
    /// ... and at least this fraction of all the reads that match it best
    /// (unique + tied) must be its own (0.25): a reference whose reads
    /// mostly tie with a better-supported relative (a sample strain that is
    /// neither of two indexed strains) is not called.
    pub min_unique_fraction: f64,
    /// When more than this fraction of the reads that resemble an organism
    /// match it only weakly, the few that pass are the best-matching reads
    /// of a relative that is not in the database (0.5).
    pub max_weak_fraction: f64,
    pub min_rpm: f64,
    /// Sample RPM / blank RPM needed for a detection also seen in a blank.
    pub min_blank_ratio: f64,
    /// Detection rate that defines the limit of detection (0.95).
    pub lod_rate: f64,
    pub output_tsv: Option<String>,
    pub output_json: Option<String>,
    pub output_html: Option<String>,
}

impl Default for ScreenOptions {
    fn default() -> Self {
        ScreenOptions {
            negative_controls: Vec::new(),
            spike_ins: Vec::new(),
            spike_table: None,
            min_hits: 3,
            min_fraction: 0.2,
            unique_margin: 0.2,
            host_min_hits: 1,
            min_reads: 3,
            min_windows: 2,
            min_dispersion: 0.5,
            min_unique_fraction: 0.25,
            max_weak_fraction: 0.5,
            min_rpm: 0.0,
            min_blank_ratio: 10.0,
            lod_rate: 0.95,
            output_tsv: None,
            output_json: None,
            output_html: None,
        }
    }
}

impl ScreenOptions {
    fn json(&self) -> String {
        format!(
            "{{\"negative_controls\":{},\"min_hits\":{},\"min_fraction\":{},\"unique_margin\":{},\"host_min_hits\":{},\"min_reads\":{},\"min_windows\":{},\"min_dispersion\":{},\"min_unique_fraction\":{},\"max_weak_fraction\":{},\"min_rpm\":{},\"min_blank_ratio\":{},\"lod_rate\":{},\"spike_table\":{}}}",
            jlist(&self.negative_controls),
            self.min_hits,
            jr(self.min_fraction, 6),
            jr(self.unique_margin, 6),
            self.host_min_hits,
            self.min_reads,
            self.min_windows,
            jr(self.min_dispersion, 6),
            jr(self.min_unique_fraction, 6),
            jr(self.max_weak_fraction, 6),
            jr(self.min_rpm, 6),
            jr(self.min_blank_ratio, 6),
            jr(self.lod_rate, 6),
            self.spike_table.as_deref().map_or_else(|| "null".to_string(), js)
        )
    }

    fn check(&self) -> Result<(), String> {
        let frac = |v: f64, name: &str| -> Result<(), String> {
            if (0.0..=1.0).contains(&v) {
                Ok(())
            } else {
                Err(format!("Invalid {name}: {v}. It is a fraction: give a number from 0 to 1"))
            }
        };
        frac(self.min_fraction, "min_fraction")?;
        frac(self.unique_margin, "unique_margin")?;
        frac(self.min_unique_fraction, "min_unique_fraction")?;
        frac(self.max_weak_fraction, "max_weak_fraction")?;
        if self.lod_rate <= 0.0 || self.lod_rate > 1.0 {
            return Err(format!("Invalid lod_rate: {}. It is the detection rate that defines the limit of detection: give a number above 0 and up to 1 (0.95 is usual)", self.lod_rate));
        }
        if self.min_hits == 0 {
            return Err("Invalid min_hits: 0. A read needs at least 1 minimizer in common with a reference (3 by default)".to_string());
        }
        if !(self.min_dispersion >= 0.0 && self.min_rpm >= 0.0 && self.min_blank_ratio >= 0.0) {
            return Err("min_dispersion, min_rpm and min_blank_ratio cannot be negative".to_string());
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Call {
    Detected,
    Background,
    BelowThreshold,
    NotDetected,
}

impl Call {
    pub fn name(self) -> &'static str {
        match self {
            Call::Detected => "detected",
            Call::Background => "background",
            Call::BelowThreshold => "below_threshold",
            Call::NotDetected => "not_detected",
        }
    }
    fn label(self) -> &'static str {
        match self {
            Call::Detected => "DETECTED",
            Call::Background => "background",
            Call::BelowThreshold => "below threshold",
            Call::NotDetected => "not detected",
        }
    }
}

/// The evidence for one organism in one sample.
#[derive(Clone, Debug, PartialEq)]
pub struct OrganismResult {
    pub name: String,
    pub kind: AgentKind,
    pub genome_bp: u64,
    /// Reads (pairs) given to it alone.
    pub unique_reads: u64,
    /// Reads tied between it and other references.
    pub tied_reads: u64,
    /// Tied reads it wins (it has the most unique reads of the tie).
    pub won_reads: u64,
    /// Its share of the tied reads (in proportion to the unique reads).
    pub shared_reads: f64,
    /// unique + shared.
    pub reads: f64,
    /// Reads whose best match was it, but weak (below min_fraction).
    pub weak_reads: u64,
    /// reads per million fragments of the sample.
    pub rpm: f64,
    pub unique_rpm: f64,
    /// Distinct minimizers of its index matched (by its unique and tied
    /// reads), and of those, the ones only it has.
    pub distinct_kmers: u64,
    pub specific_kmers: u64,
    pub index_kmers: u64,
    /// Genome windows hit by the reads that count for it (its unique
    /// reads, and its tied reads when it wins ties); `windows_any`: by all
    /// its reads.
    pub windows: u32,
    pub windows_any: u32,
    pub genome_windows: u32,
    pub window_bp: u32,
    pub breadth: f64,
    pub breadth_any: f64,
    /// One character per genome window, for the coverage plot: '2' hit by
    /// its unique reads, '1' only by reads it shares (tied), '0' not hit.
    pub window_map: String,
    pub expected_windows: f64,
    pub dispersion: Option<f64>,
    pub unique_fraction: Option<f64>,
    /// Identity estimated from the fraction of the read minimizers matched
    /// (sequencing errors lower it by about the error rate).
    pub identity: Option<f64>,
    pub shared_with: Option<String>,
    pub call: Call,
    pub reasons: Vec<String>,
    pub flags: Vec<String>,
    pub blank: Option<String>,
    pub blank_rpm: Option<f64>,
    pub blank_ratio: Option<f64>,
    pub spiked_amount: Option<f64>,
    pub spiked_unit: Option<String>,
    pub recovery: Option<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct InputFile {
    pub path: String,
    pub bytes: u64,
    pub sha256: String,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SampleResult {
    pub name: String,
    pub negative_control: bool,
    pub paired: bool,
    pub files: Vec<InputFile>,
    /// Reads, or pairs when paired.
    pub fragments: u64,
    pub host: u64,
    pub unique: u64,
    pub tied: u64,
    pub weak: u64,
    pub unclassified: u64,
    pub no_kmers: u64,
    pub mean_minimizers: f64,
    pub organisms: Vec<OrganismResult>,
    /// Agents detected (not controls).
    pub detected: Vec<String>,
    pub seconds: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SpikeResult {
    pub sample: String,
    pub organism: String,
    pub amount: f64,
    pub unit: String,
    pub reads: f64,
    pub unique_reads: u64,
    pub expected_reads: Option<f64>,
    pub recovery: Option<f64>,
    pub reads_per_unit: f64,
    pub call: Call,
}

#[derive(Clone, Debug, PartialEq)]
pub struct LodLevel {
    pub amount: f64,
    pub replicates: u32,
    pub detected: u32,
    pub rate: f64,
    pub mean_reads: f64,
    pub mean_recovery: Option<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct LodResult {
    pub organism: String,
    pub unit: String,
    pub levels: Vec<LodLevel>,
    /// The lowest level detected in at least one replicate.
    pub lowest_detected: Option<f64>,
    /// The lowest level at which this and every higher level reach
    /// `lod_rate`.
    pub lod: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct ScreenReport {
    pub engine: String,
    pub date: String,
    pub database: ScreenDbSummary,
    pub options: ScreenOptions,
    pub samples: Vec<SampleResult>,
    pub spikes: Vec<SpikeResult>,
    pub lod: Vec<LodResult>,
    pub warnings: Vec<String>,
    pub seconds: f64,
}

// ---------------------------------------------------------------------
//  THE SPIKE-IN TABLE
// ---------------------------------------------------------------------

/// Reads a spike-in table: tab-separated (or comma-separated), a header
/// with `organism` and `amount`, and optionally `sample` (empty or `*`:
/// every sample that is not a negative control) and `unit` (default
/// "reads"). Lines starting with '#' are comments.
pub fn read_spike_table(path: &str) -> Result<Vec<SpikeIn>, String> {
    check_input(path)?;
    let text = std::fs::read_to_string(path).map_err(|e| format!("Error reading the spike-in table {path}: {e}. It must be a text file (TSV)"))?;
    let mut header: Option<(Vec<String>, char)> = None;
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        let Some((cols, sep)) = &header else {
            let sep = if line.contains('\t') { '\t' } else { ',' };
            let cols: Vec<String> = line.split(sep).map(|c| c.trim().to_ascii_lowercase()).collect();
            if !cols.iter().any(|c| c == "organism") || !cols.iter().any(|c| c == "amount") {
                return Err(format!(
                    "{path}, line {}: the spike-in table needs a header with the columns organism and amount (and optionally sample and unit), separated by tabs; found \"{}\"",
                    i + 1,
                    line.chars().take(80).collect::<String>()
                ));
            }
            header = Some((cols, sep));
            continue;
        };
        let cells: Vec<&str> = line.split(*sep).map(str::trim).collect();
        let get = |name: &str| cols.iter().position(|c| c == name).and_then(|j| cells.get(j).copied()).unwrap_or("");
        let organism = get("organism");
        if organism.is_empty() {
            return Err(format!("{path}, line {}: the organism is empty. Give the name of the organism as the screen database has it", i + 1));
        }
        let amount: f64 = get("amount").replace(',', ".").parse().map_err(|_| {
            format!("{path}, line {}: the amount \"{}\" is not a number. Write the amount spiked, for instance 30 or 1.5e3", i + 1, get("amount"))
        })?;
        if !(amount.is_finite() && amount > 0.0) {
            return Err(format!("{path}, line {}: the amount must be a number above 0, found {amount}", i + 1));
        }
        let sample = get("sample");
        let unit = get("unit");
        out.push(SpikeIn {
            sample: if sample.is_empty() || sample == "*" { None } else { Some(sample.to_string()) },
            organism: organism.to_string(),
            amount,
            unit: if unit.is_empty() { "reads".to_string() } else { unit.to_string() },
        });
    }
    if header.is_none() {
        return Err(format!("{path}: the spike-in table is empty. It needs a header (sample, organism, amount, unit) and one line per spike-in"));
    }
    Ok(out)
}

// ---------------------------------------------------------------------
//  THE SCREEN: one sample
// ---------------------------------------------------------------------

#[derive(Clone, Default)]
struct RefAcc {
    unique: u64,
    hits: u64,
    mins: u64,
    winreads: u64,
    weak: u64,
}

struct SampleRaw {
    paired: bool,
    fragments: u64,
    host: u64,
    unique: u64,
    tied: u64,
    weak: u64,
    unclassified: u64,
    no_kmers: u64,
    mins: u64,
    acc: Vec<RefAcc>,
    ties: HashMap<Vec<u32>, u64>,
    windows: Vec<u32>,
    windows_any: Vec<u32>,
    window_maps: Vec<String>,
    distinct: Vec<u64>,
    specific: Vec<u64>,
    seconds: f64,
}

fn run_sample(db: &ScreenDb, s: &ScreenSample, p: &ReadParams) -> Result<SampleRaw, String> {
    let t0 = Instant::now();
    let n_refs = db.summary.references.len();
    let marks = Marks {
        slots: Bits::new(db.keys.len() + db.multi.len()),
        win_unique: Bits::new(db.total_windows),
        win_any: Bits::new(db.total_windows),
    };
    let mut raw = SampleRaw {
        paired: !s.r2.is_empty(),
        fragments: 0,
        host: 0,
        unique: 0,
        tied: 0,
        weak: 0,
        unclassified: 0,
        no_kmers: 0,
        mins: 0,
        acc: vec![RefAcc::default(); n_refs],
        ties: HashMap::new(),
        windows: vec![0; n_refs],
        windows_any: vec![0; n_refs],
        window_maps: vec![String::new(); n_refs],
        distinct: vec![0; n_refs],
        specific: vec![0; n_refs],
        seconds: 0.0,
    };
    for_each_fragments(s, 16_384, |r1, r2| {
        let calls: Vec<FragCall> = (0..r1.len())
            .into_par_iter()
            .map_init(Scratch::default, |sc, i| classify_fragment(db, &r1[i].seq, r2.map(|r| r[i].seq.as_slice()), p, sc, Some(&marks)))
            .collect();
        for c in calls {
            raw.fragments += 1;
            raw.mins += u64::from(c.minimizers);
            match c.outcome {
                Outcome::NoKmers => raw.no_kmers += 1,
                Outcome::Unclassified => raw.unclassified += 1,
                Outcome::Host => raw.host += 1,
                Outcome::Weak(r) => {
                    raw.weak += 1;
                    raw.acc[r as usize].weak += 1;
                }
                Outcome::Unique(r) => {
                    raw.unique += 1;
                    let a = &mut raw.acc[r as usize];
                    a.unique += 1;
                    a.hits += u64::from(c.hits);
                    a.mins += u64::from(c.minimizers);
                    a.winreads += u64::from(c.windows);
                }
                Outcome::Tie(set) => {
                    raw.tied += 1;
                    *raw.ties.entry(set).or_insert(0) += 1;
                }
            }
        }
        Ok(())
    })?;
    for (i, r) in db.summary.references.iter().enumerate() {
        if r.kind != AgentKind::Host {
            let (s, e) = (r.first_window as usize, (r.first_window + r.windows) as usize);
            raw.windows[i] = marks.win_unique.count(s, e);
            raw.windows_any[i] = marks.win_any.count(s, e);
            // o mapa das janelas, para o grafico de cobertura do genoma
            raw.window_maps[i] = (s..e)
                .map(|g| if marks.win_unique.get(g) { '2' } else if marks.win_any.get(g) { '1' } else { '0' })
                .collect();
        }
    }
    let nk = db.keys.len();
    for (w, word) in marks.slots.0.iter().enumerate() {
        let mut bits = word.load(Ordering::Relaxed);
        while bits != 0 {
            let slot = w * 64 + bits.trailing_zeros() as usize;
            bits &= bits - 1;
            if slot < nk {
                let r = db.win_ref[db.occ[slot] as usize] as usize;
                raw.distinct[r] += 1;
                raw.specific[r] += 1;
            } else {
                let r = db.win_ref[db.multi[slot - nk] as usize] as usize;
                raw.distinct[r] += 1;
            }
        }
    }
    raw.seconds = t0.elapsed().as_secs_f64();
    Ok(raw)
}

fn expected_windows(windows: u32, reads: u64, per_read: f64) -> f64 {
    if windows == 0 || reads == 0 {
        return 0.0;
    }
    let w = f64::from(windows);
    let p = (per_read.max(1.0) / w).min(1.0);
    w * (1.0 - (1.0 - p).powf(reads as f64))
}

/// The rows of one sample, with the calls before the negative controls.
///
/// Tied reads: for the abundance (`reads`, RPM) a tie set is shared in
/// proportion to the unique reads of its members; for the call it is WON by
/// its member with the most unique reads (the lowest reference on equal
/// counts), as the greedy cover of `gather` would give it: two
/// near-identical references do not hide the organism from each other, and
/// the one that only shares reads with a better-supported relative is not
/// called.
fn organism_rows(db: &ScreenDb, raw: &SampleRaw, o: &ScreenOptions) -> Vec<OrganismResult> {
    let refs = &db.summary.references;
    let n = refs.len();
    let mut shared = vec![0f64; n];
    let mut tied = vec![0u64; n];
    let mut won = vec![0u64; n];
    let mut partner: HashMap<(u32, u32), u64> = HashMap::new();
    let mut sets: Vec<(&Vec<u32>, &u64)> = raw.ties.iter().collect();
    sets.sort();
    for (set, &c) in sets {
        let unique = |x: u32| raw.acc[x as usize].unique;
        let base: u64 = set.iter().map(|&x| unique(x)).sum();
        let winner = set.iter().copied().max_by(|&a, &b| unique(a).cmp(&unique(b)).then(b.cmp(&a))).unwrap_or(set[0]);
        won[winner as usize] += c;
        for &x in set.iter() {
            tied[x as usize] += c;
            shared[x as usize] += if base > 0 {
                c as f64 * unique(x) as f64 / base as f64
            } else if x == winner {
                c as f64
            } else {
                0.0
            };
            for &y in set.iter() {
                if y != x {
                    *partner.entry((x, y)).or_insert(0) += c;
                }
            }
        }
    }
    let frags = raw.fragments.max(1) as f64;
    let mut rows = Vec::new();
    for i in 0..n {
        let r = &refs[i];
        if r.kind == AgentKind::Host {
            continue;
        }
        let a = &raw.acc[i];
        if a.unique == 0 && tied[i] == 0 && a.weak == 0 {
            continue;
        }
        let reads = a.unique as f64 + shared[i];
        let assigned = a.unique + won[i];
        // The windows of the reads that count for it: its unique reads, and
        // its tied reads when it wins ties.
        let windows = if won[i] > 0 { raw.windows_any[i] } else { raw.windows[i] };
        let per_read = if a.unique > 0 { a.winreads as f64 / a.unique as f64 } else { 1.0 };
        let expected = expected_windows(r.windows, assigned, per_read);
        let dispersion = if expected > 0.0 { Some(f64::from(windows) / expected) } else { None };
        let unique_fraction = if a.unique + tied[i] > 0 { Some(assigned as f64 / (a.unique + tied[i]) as f64) } else { None };
        let identity = if a.mins > 0 && a.hits > 0 { Some((a.hits as f64 / a.mins as f64).powf(1.0 / db.k as f64)) } else { None };
        let best_partner = partner
            .iter()
            .filter(|((x, _), _)| *x as usize == i)
            .max_by(|p, q| p.1.cmp(q.1).then(q.0 .1.cmp(&p.0 .1)))
            .map(|((_, y), _)| refs[*y as usize].name.clone());
        let mut row = OrganismResult {
            name: r.name.clone(),
            kind: r.kind,
            genome_bp: r.bases,
            unique_reads: a.unique,
            tied_reads: tied[i],
            won_reads: won[i],
            shared_reads: shared[i],
            reads,
            weak_reads: a.weak,
            rpm: reads / frags * 1e6,
            unique_rpm: a.unique as f64 / frags * 1e6,
            distinct_kmers: raw.distinct[i],
            specific_kmers: raw.specific[i],
            index_kmers: r.minimizers,
            windows,
            windows_any: raw.windows_any[i],
            genome_windows: r.windows,
            window_bp: r.window_bp,
            breadth: if r.windows > 0 { f64::from(windows) / f64::from(r.windows) } else { 0.0 },
            breadth_any: if r.windows > 0 { f64::from(raw.windows_any[i]) / f64::from(r.windows) } else { 0.0 },
            window_map: raw.window_maps[i].clone(),
            expected_windows: expected,
            dispersion,
            unique_fraction,
            identity,
            shared_with: best_partner,
            call: Call::NotDetected,
            reasons: Vec::new(),
            flags: Vec::new(),
            blank: None,
            blank_rpm: None,
            blank_ratio: None,
            spiked_amount: None,
            spiked_unit: None,
            recovery: None,
        };
        call_row(&mut row, o);
        rows.push(row);
    }
    rows
}

fn call_row(row: &mut OrganismResult, o: &ScreenOptions) {
    row.reasons.clear();
    if row.unique_reads == 0 && row.tied_reads == 0 {
        row.call = if row.weak_reads > 0 { Call::BelowThreshold } else { Call::NotDetected };
        if row.weak_reads > 0 {
            row.reasons.push(format!(
                "{} read(s) matched it only weakly (less than {} of their k-mers): possibly a relative that is not in the database, or chimeric reads; not counted",
                thousands(row.weak_reads),
                pct(o.min_fraction)
            ));
        }
        return;
    }
    let partner = row.shared_with.clone().unwrap_or_else(|| "other references".to_string());
    let assigned = row.unique_reads + row.won_reads;
    if assigned == 0 {
        row.reasons.push(format!(
            "no read is its own: its {} read(s) tie with {partner}, which has more specific reads and explains them",
            thousands(row.tied_reads)
        ));
    } else {
        if assigned < o.min_reads {
            row.reasons.push(format!(
                "{assigned} read(s) for it ({} unique, {} tied reads won), fewer than the minimum of {}",
                row.unique_reads, row.won_reads, o.min_reads
            ));
        }
        if row.windows < o.min_windows {
            row.reasons.push(format!("its reads hit {} genome window(s), fewer than the minimum of {}", row.windows, o.min_windows));
        }
        if let Some(d) = row.dispersion {
            if d < o.min_dispersion && row.windows >= o.min_windows.min(1) {
                row.reasons.push(format!(
                    "its reads are concentrated: {} of {} genome windows hit where {} would be expected if {assigned} reads were spread over the genome (dispersion {}, minimum {}); reads piled on one region (a vector element, an integrated fragment, a gene shared with other organisms) are not a detection",
                    row.windows,
                    row.genome_windows,
                    short(row.expected_windows),
                    short(d),
                    short(o.min_dispersion)
                ));
            }
        }
        if let Some(f) = row.unique_fraction {
            if f < o.min_unique_fraction {
                row.reasons.push(format!(
                    "only {} of its reads are its own ({} unique, {} tied, mostly with {partner}, which has more specific reads): what is in the sample is probably {partner} or a strain between the two, not this one",
                    pct(f),
                    row.unique_reads,
                    row.tied_reads
                ));
            }
        }
        let weak_fraction = row.weak_reads as f64 / (row.weak_reads + assigned) as f64;
        if weak_fraction > o.max_weak_fraction {
            row.reasons.push(format!(
                "most reads that resemble it match it only weakly ({} weak, {assigned} of its own): the few that pass are probably the best-matching reads of a relative that is not in the database (more than {} weak)",
                thousands(row.weak_reads),
                pct(o.max_weak_fraction)
            ));
        }
        if row.rpm < o.min_rpm {
            row.reasons.push(format!("{} reads per million, below the minimum of {}", short(row.rpm), short(o.min_rpm)));
        }
    }
    row.call = if row.reasons.is_empty() { Call::Detected } else { Call::BelowThreshold };
    if row.won_reads > 0 && row.call == Call::Detected && row.won_reads > row.unique_reads {
        row.flags.push(format!(
            "most of its reads ({} of {}) tie with {partner}: the two are not well told apart in this sample (the tied reads go to the one with more specific reads)",
            thousands(row.won_reads),
            thousands(row.unique_reads + row.tied_reads)
        ));
    }
    if row.weak_reads > 0 {
        row.flags.push(format!("{} more read(s) matched it only weakly (a relative, or chimeric reads)", thousands(row.weak_reads)));
    }
}

// ---------------------------------------------------------------------
//  THE SCREEN
// ---------------------------------------------------------------------

/// Screens the reads of each sample for the agents of the database: host
/// reads set aside, each read given to the reference it matches best (or
/// tied), and for each organism the evidence (unique and tied reads, RPM,
/// distinct k-mers, breadth over its genome windows) and a call; organisms
/// also in a negative control are background unless far above it;
/// spike-ins get their recovery and the run its limit of detection. The
/// report (TSV, JSON, HTML) records the SHA-256 of every input and of the
/// database, all parameters, the date and the engine version. For research
/// and teaching (see `NOTICE`).
pub fn contaminant_screen(samples: &[ScreenSample], database: &str, opt: &ScreenOptions) -> Result<ScreenReport, String> {
    let t0 = Instant::now();
    opt.check()?;
    if samples.is_empty() {
        return Err("No sample to screen: give at least one sample (a name and its FASTQ files)".to_string());
    }
    let mut inputs: Vec<String> = vec![database.to_string()];
    for (i, s) in samples.iter().enumerate() {
        if s.name.trim().is_empty() {
            return Err(format!("Sample {} has no name: give every sample a name (it is how the report and the negative controls refer to it)", i + 1));
        }
        if samples[..i].iter().any(|t| t.name == s.name) {
            return Err(format!("Two samples are named \"{}\": give each sample its own name", s.name));
        }
        if s.r1.is_empty() {
            return Err(format!("Sample {} has no read file: give its FASTQ (or FASTA) file(s)", s.name));
        }
        if !s.r2.is_empty() && s.r2.len() != s.r1.len() {
            return Err(format!(
                "Sample {} has {} first-read file(s) and {} second-read file(s): give one second-read file for each first-read file (lanes in the same order), or none for single-end reads",
                s.name,
                s.r1.len(),
                s.r2.len()
            ));
        }
        for f in s.r1.iter().chain(s.r2.iter()) {
            check_input(f)?;
            inputs.push(f.clone());
        }
        for (a, b) in s.r1.iter().zip(s.r2.iter()) {
            if same_file(a, b) {
                return Err(format!("Sample {}: the two files of a pair are the same file ({a}): give the file of the first reads and the file of the second reads", s.name));
            }
        }
    }
    for c in &opt.negative_controls {
        if !samples.iter().any(|s| &s.name == c) {
            let names: Vec<&str> = samples.iter().map(|s| s.name.as_str()).collect();
            return Err(format!("The negative control \"{c}\" is not one of the samples ({}): give the name of a sample", names.join(", ")));
        }
    }
    let outputs: Vec<&String> = [&opt.output_tsv, &opt.output_json, &opt.output_html].into_iter().flatten().collect();
    let in_refs: Vec<&str> = inputs.iter().map(String::as_str).collect();
    for (i, out) in outputs.iter().enumerate() {
        not_an_input(out, &in_refs)?;
        if outputs[..i].iter().any(|p| same_file(p, out)) {
            return Err(format!("Two outputs are the same file ({out}): give each output its own name"));
        }
    }
    let mut spikes = opt.spike_ins.clone();
    if let Some(t) = &opt.spike_table {
        spikes.extend(read_spike_table(t)?);
    }

    let db = load_screen_database(database)?;
    let refs = &db.summary.references;
    for sp in &spikes {
        match refs.iter().find(|r| r.name == sp.organism) {
            None => {
                let known: Vec<&str> = refs.iter().filter(|r| r.kind != AgentKind::Host).map(|r| r.name.as_str()).take(15).collect();
                return Err(format!(
                    "The spike-in organism \"{}\" is not a reference of the screen database. Use the names as the database has them; some of them: {}",
                    sp.organism,
                    known.join(" | ")
                ));
            }
            Some(r) if r.kind == AgentKind::Host => {
                return Err(format!("The spike-in organism \"{}\" is a host reference: a spike-in must be an agent or a control of the database", sp.organism))
            }
            _ => {}
        }
        if let Some(sn) = &sp.sample {
            if !samples.iter().any(|s| &s.name == sn) {
                return Err(format!("The spike-in table names the sample \"{sn}\", which is not one of the samples screened: check the name"));
            }
        }
        if !(sp.amount.is_finite() && sp.amount > 0.0) {
            return Err(format!("The spike-in of {} has amount {}: give a number above 0", sp.organism, sp.amount));
        }
    }
    let p = ReadParams { min_hits: opt.min_hits, min_fraction: opt.min_fraction, unique_margin: opt.unique_margin, host_min_hits: opt.host_min_hits };
    let read_files: Vec<String> = inputs[1..].to_vec();

    let (hashed, raws) = std::thread::scope(|s| {
        let hasher = s.spawn(|| hash_files(&read_files));
        let mut raws = Vec::with_capacity(samples.len());
        for smp in samples {
            match run_sample(&db, smp, &p) {
                Ok(r) => raws.push(r),
                Err(e) => return (hasher.join(), Err(e)),
            }
        }
        (hasher.join(), Ok(raws))
    });
    let raws = raws?;
    let hashed = hashed.map_err(|_| "Internal error while computing the SHA-256 of the read files".to_string())??;

    let mut results: Vec<SampleResult> = Vec::with_capacity(samples.len());
    for (smp, raw) in samples.iter().zip(&raws) {
        let files = smp
            .r1
            .iter()
            .chain(smp.r2.iter())
            .map(|f| {
                let (bytes, sha) = hashed.get(f).cloned().unwrap_or_default();
                InputFile { path: f.clone(), bytes, sha256: sha }
            })
            .collect();
        results.push(SampleResult {
            name: smp.name.clone(),
            negative_control: opt.negative_controls.contains(&smp.name),
            paired: raw.paired,
            files,
            fragments: raw.fragments,
            host: raw.host,
            unique: raw.unique,
            tied: raw.tied,
            weak: raw.weak,
            unclassified: raw.unclassified,
            no_kmers: raw.no_kmers,
            mean_minimizers: if raw.fragments > 0 { raw.mins as f64 / raw.fragments as f64 } else { 0.0 },
            organisms: organism_rows(&db, raw, opt),
            detected: Vec::new(),
            seconds: raw.seconds,
        });
    }
    let mut warnings: Vec<String> = Vec::new();

    // Negative controls.
    let blanks: Vec<usize> = (0..results.len()).filter(|&i| results[i].negative_control).collect();
    for i in 0..results.len() {
        if results[i].negative_control {
            continue;
        }
        let mut rows = std::mem::take(&mut results[i].organisms);
        // A process control is added to the blanks too: it is never
        // background.
        for row in rows.iter_mut().filter(|r| r.kind != AgentKind::Spike) {
            let mut worst: Option<(String, f64)> = None;
            for &b in &blanks {
                if let Some(br) = results[b].organisms.iter().find(|x| x.name == row.name && x.reads > 0.0) {
                    if worst.as_ref().is_none_or(|w| br.rpm > w.1) {
                        worst = Some((results[b].name.clone(), br.rpm));
                    }
                }
            }
            if let Some((bname, brpm)) = worst {
                let ratio = if brpm > 0.0 { row.rpm / brpm } else { f64::INFINITY };
                row.blank = Some(bname.clone());
                row.blank_rpm = Some(brpm);
                row.blank_ratio = Some(ratio);
                if row.call == Call::Detected && ratio < opt.min_blank_ratio {
                    row.call = Call::Background;
                    row.reasons.push(format!(
                        "also in the negative control {bname} at {} RPM: the sample has {} RPM, {}x the blank, below the {}x needed to call it above background",
                        short(brpm),
                        short(row.rpm),
                        short(ratio),
                        short(opt.min_blank_ratio)
                    ));
                } else {
                    row.flags.push(format!("also in the negative control {bname} at {} RPM (sample/blank = {}x)", short(brpm), short(ratio)));
                }
            }
        }
        results[i].organisms = rows;
    }

    // Spike-ins.
    let mut spike_rows: Vec<SpikeResult> = Vec::new();
    for sp in &spikes {
        let r = refs.iter().find(|r| r.name == sp.organism).cloned().unwrap_or_default();
        // A spike-in without a sample goes to every sample; an agent spiked
        // for a limit of detection is not put into the blanks, a process
        // control is.
        let targets: Vec<usize> = match &sp.sample {
            Some(n) => results.iter().position(|r| &r.name == n).into_iter().collect(),
            None => (0..results.len()).filter(|&i| r.kind == AgentKind::Spike || !results[i].negative_control).collect(),
        };
        for i in targets {
            let frags = results[i].fragments;
            let sample_name = results[i].name.clone();
            let pos = match results[i].organisms.iter().position(|x| x.name == sp.organism) {
                Some(p) => p,
                None => {
                    results[i].organisms.push(empty_row(&r));
                    results[i].organisms.len() - 1
                }
            };
            let row = &mut results[i].organisms[pos];
            let unit = sp.unit.trim().to_ascii_lowercase();
            let expected = match unit.as_str() {
                "reads" | "read" | "pairs" | "fragments" => Some(sp.amount),
                "rpm" => Some(sp.amount * frags as f64 / 1e6),
                _ => None,
            };
            let recovery = expected.filter(|e| *e > 0.0).map(|e| row.reads / e);
            row.spiked_amount = Some(row.spiked_amount.unwrap_or(0.0) + sp.amount);
            row.spiked_unit = Some(sp.unit.clone());
            row.recovery = recovery;
            spike_rows.push(SpikeResult {
                sample: sample_name,
                organism: sp.organism.clone(),
                amount: sp.amount,
                unit: sp.unit.clone(),
                reads: row.reads,
                unique_reads: row.unique_reads,
                expected_reads: expected,
                recovery,
                reads_per_unit: row.reads / sp.amount,
                call: row.call,
            });
        }
    }
    // The limit of detection is for agents; a process control has its
    // recovery only.
    let agent_rows: Vec<SpikeResult> =
        spike_rows.iter().filter(|x| refs.iter().any(|r| r.name == x.organism && r.kind != AgentKind::Spike)).cloned().collect();
    let lod = limits_of_detection(&agent_rows, opt.lod_rate);
    let fewest = lod.iter().flat_map(|l| l.levels.iter().map(|v| v.replicates)).min();
    if let Some(n) = fewest.filter(|&n| n < 20) {
        warnings.push(format!(
            "Some spike-in level has only {n} replicate(s): with fewer than 20, a detection rate of {} means every replicate, and the limit of detection is only an estimate (a validation study uses 20 or more replicates per level)",
            pct(opt.lod_rate)
        ));
    }

    for r in results.iter_mut() {
        r.organisms.sort_by(|a, b| a.call.cmp(&b.call).then(b.reads.total_cmp(&a.reads)).then(a.name.cmp(&b.name)));
        r.detected = r.organisms.iter().filter(|x| x.call == Call::Detected && x.kind.is_agent()).map(|x| x.name.clone()).collect();
        if r.fragments == 0 {
            warnings.push(format!("Sample {}: no read in its files", r.name));
        } else if r.mean_minimizers < 4.0 {
            warnings.push(format!(
                "Sample {}: its reads hold only {:.1} minimizers on average: reads this short give little evidence",
                r.name, r.mean_minimizers
            ));
        }
        if !r.negative_control && db.host.is_empty() {
            continue;
        }
    }
    if db.host.is_empty() {
        warnings.push("The database has no host reference: host reads were not set aside".to_string());
    }
    if blanks.is_empty() {
        warnings.push("No negative control was given: reagent and laboratory background cannot be told from the sample".to_string());
    }

    let report = ScreenReport {
        engine: engine_version(),
        date: utc_now(),
        database: db.summary.clone(),
        options: ScreenOptions { spike_ins: spikes, ..opt.clone() },
        samples: results,
        spikes: spike_rows,
        lod,
        warnings,
        seconds: t0.elapsed().as_secs_f64(),
    };
    let mut cleanup = Cleanup(outputs.iter().map(|s| s.to_string()).collect());
    if let Some(p) = &opt.output_tsv {
        write_file(p, &report.tsv())?;
    }
    if let Some(p) = &opt.output_json {
        write_file(p, &report.json())?;
    }
    if let Some(p) = &opt.output_html {
        write_file(p, &report.html())?;
    }
    cleanup.disarm();
    Ok(report)
}

fn empty_row(r: &ScreenRefInfo) -> OrganismResult {
    OrganismResult {
        name: r.name.clone(),
        kind: r.kind,
        genome_bp: r.bases,
        unique_reads: 0,
        tied_reads: 0,
        won_reads: 0,
        shared_reads: 0.0,
        reads: 0.0,
        weak_reads: 0,
        rpm: 0.0,
        unique_rpm: 0.0,
        distinct_kmers: 0,
        specific_kmers: 0,
        index_kmers: r.minimizers,
        windows: 0,
        windows_any: 0,
        genome_windows: r.windows,
        window_bp: r.window_bp,
        breadth: 0.0,
        breadth_any: 0.0,
        window_map: String::new(),
        expected_windows: 0.0,
        dispersion: None,
        unique_fraction: None,
        identity: None,
        shared_with: None,
        call: Call::NotDetected,
        reasons: Vec::new(),
        flags: Vec::new(),
        blank: None,
        blank_rpm: None,
        blank_ratio: None,
        spiked_amount: None,
        spiked_unit: None,
        recovery: None,
    }
}

fn limits_of_detection(rows: &[SpikeResult], rate: f64) -> Vec<LodResult> {
    let mut groups: Vec<(String, String, Vec<&SpikeResult>)> = Vec::new();
    for r in rows {
        let unit = r.unit.trim().to_ascii_lowercase();
        match groups.iter_mut().find(|g| g.0 == r.organism && g.1 == unit) {
            Some(g) => g.2.push(r),
            None => groups.push((r.organism.clone(), unit, vec![r])),
        }
    }
    let mut out = Vec::new();
    for (organism, unit, mut v) in groups {
        v.sort_by(|a, b| a.amount.total_cmp(&b.amount));
        let mut levels: Vec<LodLevel> = Vec::new();
        let mut i = 0;
        while i < v.len() {
            let amount = v[i].amount;
            let mut j = i;
            while j < v.len() && (v[j].amount - amount).abs() <= 1e-9 * amount.abs().max(1.0) {
                j += 1;
            }
            let g = &v[i..j];
            let detected = g.iter().filter(|x| x.call == Call::Detected).count() as u32;
            let recs: Vec<f64> = g.iter().filter_map(|x| x.recovery).collect();
            levels.push(LodLevel {
                amount,
                replicates: g.len() as u32,
                detected,
                rate: f64::from(detected) / g.len() as f64,
                mean_reads: g.iter().map(|x| x.reads).sum::<f64>() / g.len() as f64,
                mean_recovery: if recs.is_empty() { None } else { Some(recs.iter().sum::<f64>() / recs.len() as f64) },
            });
            i = j;
        }
        let lowest_detected = levels.iter().find(|l| l.detected > 0).map(|l| l.amount);
        let mut lod = None;
        for l in levels.iter().rev() {
            if l.rate + 1e-12 >= rate {
                lod = Some(l.amount);
            } else {
                break;
            }
        }
        out.push(LodResult { organism, unit, levels, lowest_detected, lod });
    }
    out
}

// ---------------------------------------------------------------------
//  THE REPORTS
// ---------------------------------------------------------------------

impl OrganismResult {
    fn json(&self) -> String {
        format!(
            "{{\"name\":{},\"kind\":{},\"role\":{},\"call\":{},\"genome_bp\":{},\"unique_reads\":{},\"tied_reads\":{},\"won_reads\":{},\"shared_reads\":{},\"reads\":{},\"weak_reads\":{},\"rpm\":{},\"unique_rpm\":{},\"distinct_kmers\":{},\"specific_kmers\":{},\"index_kmers\":{},\"windows\":{},\"windows_any\":{},\"genome_windows\":{},\"window_bp\":{},\"breadth\":{},\"breadth_any\":{},\"window_map\":{},\"expected_windows\":{},\"dispersion\":{},\"unique_fraction\":{},\"identity\":{},\"shared_with\":{},\"reasons\":{},\"flags\":{},\"blank\":{},\"blank_rpm\":{},\"blank_ratio\":{},\"spiked_amount\":{},\"spiked_unit\":{},\"recovery\":{}}}",
            js(&self.name),
            js(self.kind.name()),
            js(if self.kind == AgentKind::Spike { "control" } else { "agent" }),
            js(self.call.name()),
            self.genome_bp,
            self.unique_reads,
            self.tied_reads,
            self.won_reads,
            jr(self.shared_reads, 3),
            jr(self.reads, 3),
            self.weak_reads,
            jr(self.rpm, 4),
            jr(self.unique_rpm, 4),
            self.distinct_kmers,
            self.specific_kmers,
            self.index_kmers,
            self.windows,
            self.windows_any,
            self.genome_windows,
            self.window_bp,
            jr(self.breadth, 6),
            jr(self.breadth_any, 6),
            js(&self.window_map),
            jr(self.expected_windows, 3),
            jo(self.dispersion, 4),
            jo(self.unique_fraction, 4),
            jo(self.identity, 5),
            self.shared_with.as_deref().map_or_else(|| "null".to_string(), js),
            jlist(&self.reasons),
            jlist(&self.flags),
            self.blank.as_deref().map_or_else(|| "null".to_string(), js),
            jo(self.blank_rpm, 4),
            match self.blank_ratio {
                Some(r) if r.is_infinite() => "null".to_string(),
                x => jo(x, 3),
            },
            jo(self.spiked_amount, 6),
            self.spiked_unit.as_deref().map_or_else(|| "null".to_string(), js),
            jo(self.recovery, 4)
        )
    }
}

impl SampleResult {
    fn json(&self) -> String {
        let files: Vec<String> = self
            .files
            .iter()
            .map(|f| format!("{{\"path\":{},\"bytes\":{},\"sha256\":{}}}", js(&f.path), f.bytes, js(&f.sha256)))
            .collect();
        let orgs: Vec<String> = self.organisms.iter().map(OrganismResult::json).collect();
        format!(
            "{{\"name\":{},\"negative_control\":{},\"paired\":{},\"unit\":{},\"files\":[{}],\"fragments\":{},\"host\":{},\"host_fraction\":{},\"unique\":{},\"tied\":{},\"weak\":{},\"unclassified\":{},\"no_kmers\":{},\"mean_minimizers\":{},\"detected\":{},\"seconds\":{},\"organisms\":[{}]}}",
            js(&self.name),
            self.negative_control,
            self.paired,
            js(if self.paired { "pairs" } else { "reads" }),
            files.join(","),
            self.fragments,
            self.host,
            jr(self.host as f64 / self.fragments.max(1) as f64, 6),
            self.unique,
            self.tied,
            self.weak,
            self.unclassified,
            self.no_kmers,
            jr(self.mean_minimizers, 2),
            jlist(&self.detected),
            jr(self.seconds, 2),
            orgs.join(",")
        )
    }
}

impl ScreenReport {
    pub fn json(&self) -> String {
        let samples: Vec<String> = self.samples.iter().map(SampleResult::json).collect();
        let spikes: Vec<String> = self
            .spikes
            .iter()
            .map(|s| {
                format!(
                    "{{\"sample\":{},\"organism\":{},\"amount\":{},\"unit\":{},\"reads\":{},\"unique_reads\":{},\"expected_reads\":{},\"recovery\":{},\"reads_per_unit\":{},\"call\":{}}}",
                    js(&s.sample),
                    js(&s.organism),
                    jr(s.amount, 6),
                    js(&s.unit),
                    jr(s.reads, 3),
                    s.unique_reads,
                    jo(s.expected_reads, 3),
                    jo(s.recovery, 4),
                    jr(s.reads_per_unit, 6),
                    js(s.call.name())
                )
            })
            .collect();
        let lod: Vec<String> = self
            .lod
            .iter()
            .map(|l| {
                let lv: Vec<String> = l
                    .levels
                    .iter()
                    .map(|v| {
                        format!(
                            "{{\"amount\":{},\"replicates\":{},\"detected\":{},\"rate\":{},\"mean_reads\":{},\"mean_recovery\":{}}}",
                            jr(v.amount, 6),
                            v.replicates,
                            v.detected,
                            jr(v.rate, 4),
                            jr(v.mean_reads, 3),
                            jo(v.mean_recovery, 4)
                        )
                    })
                    .collect();
                format!(
                    "{{\"organism\":{},\"unit\":{},\"levels\":[{}],\"lowest_detected\":{},\"lod\":{}}}",
                    js(&l.organism),
                    js(&l.unit),
                    lv.join(","),
                    jo(l.lowest_detected, 6),
                    jo(l.lod, 6)
                )
            })
            .collect();
        let spike_in: Vec<String> = self
            .options
            .spike_ins
            .iter()
            .map(|s| {
                format!(
                    "{{\"sample\":{},\"organism\":{},\"amount\":{},\"unit\":{}}}",
                    s.sample.as_deref().map_or_else(|| "null".to_string(), js),
                    js(&s.organism),
                    jr(s.amount, 6),
                    js(&s.unit)
                )
            })
            .collect();
        format!(
            "{{\"engine\":{},\"date\":{},\"notice\":{},\"database\":{},\"parameters\":{},\"spike_ins\":[{}],\"samples\":[{}],\"spike_results\":[{}],\"lod\":[{}],\"warnings\":{},\"seconds\":{},\"text\":{}}}",
            js(&self.engine),
            js(&self.date),
            js(NOTICE),
            self.database.json(),
            self.options.json(),
            spike_in.join(","),
            samples.join(","),
            spikes.join(","),
            lod.join(","),
            jlist(&self.warnings),
            jr(self.seconds, 2),
            js(&self.resumo())
        )
    }

    /// One line per sample and organism with evidence (or spiked), after a
    /// commented header with the traceability (date, engine, database and
    /// input SHA-256, parameters).
    pub fn tsv(&self) -> String {
        let mut o = String::new();
        o.push_str(&format!("# {}\n# date {}\n", self.engine, self.date));
        o.push_str(&format!("# notice: {}\n", NOTICE));
        o.push_str(&format!("# database {} sha256 {}\n", clean_cell(&self.database.path), self.database.sha256));
        for s in &self.samples {
            for f in &s.files {
                o.push_str(&format!("# input {} {} bytes sha256 {} ({})\n", clean_cell(&f.path), f.bytes, f.sha256, clean_cell(&s.name)));
            }
        }
        o.push_str(&format!("# parameters {}\n", self.options.json()));
        o.push_str("sample\tnegative_control\tfragments\thost_fraction\torganism\tkind\tcall\tunique_reads\ttied_reads\twon_reads\treads\trpm\tdistinct_kmers\tspecific_kmers\twindows\tgenome_windows\tbreadth\texpected_windows\tdispersion\tunique_fraction\tidentity\tweak_reads\tblank\tblank_rpm\tblank_ratio\tspiked_amount\tspiked_unit\trecovery\treasons\n");
        let num = |x: Option<f64>, d: usize| x.filter(|v| v.is_finite()).map_or_else(String::new, |v| format!("{v:.d$}"));
        for s in &self.samples {
            for r in &s.organisms {
                o.push_str(&format!(
                    "{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.3}\t{}\t{}\t{}\t{}\t{:.5}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    clean_cell(&s.name),
                    s.negative_control,
                    s.fragments,
                    s.host as f64 / s.fragments.max(1) as f64,
                    clean_cell(&r.name),
                    r.kind.name(),
                    r.call.name(),
                    r.unique_reads,
                    r.tied_reads,
                    r.won_reads,
                    r.reads,
                    r.rpm,
                    r.distinct_kmers,
                    r.specific_kmers,
                    r.windows,
                    r.genome_windows,
                    r.breadth,
                    r.expected_windows,
                    num(r.dispersion, 3),
                    num(r.unique_fraction, 4),
                    num(r.identity, 4),
                    r.weak_reads,
                    r.blank.as_deref().map_or_else(String::new, clean_cell),
                    num(r.blank_rpm, 3),
                    num(r.blank_ratio, 2),
                    num(r.spiked_amount, 4),
                    r.spiked_unit.as_deref().map_or_else(String::new, clean_cell),
                    num(r.recovery, 4),
                    clean_cell(&r.reasons.join("; "))
                ));
            }
        }
        o
    }

    /// The report as text.
    pub fn resumo(&self) -> String {
        let mut o = format!(
            "Adventitious-agent screen, {} ({})\nDatabase {} (SHA-256 {}), k = {}, {} reference(s)\n{}\n",
            self.date,
            self.engine,
            self.database.path,
            self.database.sha256,
            self.database.k,
            self.database.references.len(),
            NOTICE
        );
        for s in &self.samples {
            let unit = if s.paired { "read pairs" } else { "reads" };
            let all = s.fragments.max(1) as f64;
            o.push_str(&format!(
                "\n{} {}: {} {unit}; host {}, given to a reference {}, unclassified {}{}\n",
                if s.negative_control { "Negative control" } else { "Sample" },
                s.name,
                thousands(s.fragments),
                pct(s.host as f64 / all),
                pct((s.unique + s.tied) as f64 / all),
                pct((s.unclassified + s.weak + s.no_kmers) as f64 / all),
                if s.detected.is_empty() { String::new() } else { format!("; DETECTED: {}", s.detected.join(", ")) }
            ));
            for r in &s.organisms {
                let role = if r.kind == AgentKind::Spike { " [control]" } else { "" };
                o.push_str(&format!(
                    "  {:<16} {}{role} ({}): {} unique + {} tied reads, {} RPM, {} of {} windows (dispersion {}){}\n",
                    r.call.label(),
                    r.name,
                    r.kind.name(),
                    thousands(r.unique_reads),
                    thousands(r.tied_reads),
                    short(r.rpm),
                    r.windows,
                    r.genome_windows,
                    r.dispersion.map_or_else(|| "-".to_string(), short),
                    match (r.spiked_amount, &r.spiked_unit) {
                        (Some(a), Some(u)) => format!("; spiked {} {u}{}", short(a), r.recovery.map_or_else(String::new, |x| format!(", recovery {}", pct(x)))),
                        _ => String::new(),
                    }
                ));
                for why in &r.reasons {
                    o.push_str(&format!("      - {why}\n"));
                }
                for f in &r.flags {
                    o.push_str(&format!("      note: {f}\n"));
                }
            }
        }
        if !self.lod.is_empty() {
            o.push_str("\nSpike-ins and limit of detection\n");
            for l in &self.lod {
                o.push_str(&format!(
                    "  {} ({}): lowest level detected {}, limit of detection ({} of replicates) {}\n",
                    l.organism,
                    l.unit,
                    l.lowest_detected.map_or_else(|| "none".to_string(), short),
                    pct(self.options.lod_rate),
                    l.lod.map_or_else(|| "not reached".to_string(), short)
                ));
                for v in &l.levels {
                    o.push_str(&format!(
                        "    {:>10} {}: {}/{} detected, mean {} reads{}\n",
                        short(v.amount),
                        l.unit,
                        v.detected,
                        v.replicates,
                        short(v.mean_reads),
                        v.mean_recovery.map_or_else(String::new, |r| format!(", recovery {}", pct(r)))
                    ));
                }
            }
        }
        for w in &self.warnings {
            o.push_str(&format!("warning: {w}\n"));
        }
        o
    }

    /// A self-contained HTML page (no script, no external file).
    pub fn html(&self) -> String {
        let mut o = String::with_capacity(64 * 1024);
        o.push_str(HTML_HEAD);
        o.push_str(&format!(
            "<h1>Adventitious-agent screen</h1>\n<p class=\"sub\">{} &middot; {} &middot; database <code>{}</code></p>\n<div class=\"notice\"><strong>For research and teaching.</strong> {}</div>\n",
            html(&self.date),
            html(&self.engine),
            html(&self.database.path),
            html(NOTICE.trim_start_matches("For research and teaching only. "))
        ));
        // The calls: one row per organism, one column per sample.
        let mut names: Vec<(String, AgentKind)> = Vec::new();
        for s in &self.samples {
            for r in &s.organisms {
                if r.call != Call::NotDetected || r.spiked_amount.is_some() {
                    if !names.iter().any(|n| n.0 == r.name) {
                        names.push((r.name.clone(), r.kind));
                    }
                }
            }
        }
        o.push_str("<section class=\"card\"><h2>Calls</h2><p class=\"caption\">Each cell: the call and the reads per million (RPM) of the sample. Negative controls are marked.</p>");
        if names.is_empty() {
            o.push_str("<p>No read was given to any agent of the database.</p>");
        } else {
            o.push_str("<div class=\"scroll\"><table><thead><tr><th>Organism</th><th>Kind</th>");
            for s in &self.samples {
                o.push_str(&format!("<th>{}{}</th>", html(&s.name), if s.negative_control { " <span class=\"tag\">blank</span>" } else { "" }));
            }
            o.push_str("</tr></thead><tbody>");
            for (n, kind) in &names {
                o.push_str(&format!("<tr><td class=\"org\">{}</td><td>{}</td>", html(n), kind.name()));
                for s in &self.samples {
                    match s.organisms.iter().find(|r| &r.name == n) {
                        Some(r) => o.push_str(&format!(
                            "<td class=\"nw\"><span class=\"call {}\">{}</span> <span class=\"num\">{}</span></td>",
                            badge(r).0,
                            badge(r).1,
                            short(r.rpm)
                        )),
                        None => o.push_str("<td class=\"muted\">&ndash;</td>"),
                    }
                }
                o.push_str("</tr>");
            }
            o.push_str("</tbody></table></div>");
        }
        o.push_str("</section>\n");
        // Each sample.
        for s in &self.samples {
            let unit = if s.paired { "read pairs" } else { "reads" };
            let all = s.fragments.max(1) as f64;
            o.push_str(&format!(
                "<section class=\"card\"><h2>{} {}</h2><p class=\"caption\">{} {unit}: host {}, given to a reference {}, weak {}, unclassified {}</p>",
                if s.negative_control { "Negative control" } else { "Sample" },
                html(&s.name),
                thousands(s.fragments),
                pct(s.host as f64 / all),
                pct((s.unique + s.tied) as f64 / all),
                pct(s.weak as f64 / all),
                pct((s.unclassified + s.no_kmers) as f64 / all)
            ));
            if s.organisms.is_empty() {
                o.push_str("<p>No read was given to any agent of the database.</p></section>\n");
                continue;
            }
            o.push_str("<div class=\"scroll\"><table><thead><tr><th>Organism</th><th>Call</th><th class=\"r\">Unique</th><th class=\"r\">Tied</th><th class=\"r\">RPM</th><th class=\"r\">Distinct k-mers</th><th class=\"r\">Windows</th><th class=\"r\">Breadth</th><th class=\"r\">Dispersion</th><th class=\"r\">Identity</th><th>Why</th></tr></thead><tbody>");
            for r in &s.organisms {
                let mut why: Vec<String> = r.reasons.iter().map(|x| html(x)).collect();
                why.extend(r.flags.iter().map(|x| format!("<span class=\"muted\">{}</span>", html(x))));
                if let (Some(a), Some(u)) = (r.spiked_amount, &r.spiked_unit) {
                    why.push(format!(
                        "spiked {} {}{}",
                        short(a),
                        html(u),
                        r.recovery.map_or_else(String::new, |x| format!(", recovery {}", pct(x)))
                    ));
                }
                o.push_str(&format!(
                    "<tr><td class=\"org\">{}<div class=\"muted\">{}</div></td><td><span class=\"call {}\">{}</span></td><td class=\"r\">{}</td><td class=\"r\">{}</td><td class=\"r\">{}</td><td class=\"r\">{}</td><td class=\"r\">{} / {}</td><td class=\"r\">{}</td><td class=\"r\">{}</td><td class=\"r\">{}</td><td class=\"why\">{}</td></tr>",
                    html(&r.name),
                    r.kind.name(),
                    badge(r).0,
                    badge(r).1,
                    thousands(r.unique_reads),
                    thousands(r.tied_reads),
                    short(r.rpm),
                    thousands(r.distinct_kmers),
                    r.windows,
                    r.genome_windows,
                    pct(r.breadth),
                    r.dispersion.map_or_else(|| "&ndash;".to_string(), short),
                    r.identity.map_or_else(|| "&ndash;".to_string(), |x| format!("{:.1}%", 100.0 * x)),
                    why.join("<br>")
                ));
            }
            o.push_str("</tbody></table></div></section>\n");
        }
        if !self.lod.is_empty() {
            o.push_str("<section class=\"card\"><h2>Spike-ins and limit of detection</h2>");
            o.push_str(&format!(
                "<p class=\"caption\">Limit of detection: the lowest level at which this and every higher level are detected in at least {} of the replicates.</p>",
                pct(self.options.lod_rate)
            ));
            o.push_str("<div class=\"scroll\"><table><thead><tr><th>Organism</th><th class=\"r\">Level</th><th class=\"r\">Detected</th><th class=\"r\">Mean reads</th><th class=\"r\">Mean recovery</th></tr></thead><tbody>");
            for l in &self.lod {
                for (j, v) in l.levels.iter().enumerate() {
                    o.push_str(&format!(
                        "<tr><td>{}</td><td class=\"r\">{} {}</td><td class=\"r\">{} / {}</td><td class=\"r\">{}</td><td class=\"r\">{}</td></tr>",
                        if j == 0 {
                            format!(
                                "{}<div class=\"muted\">lowest detected {}; LoD {}</div>",
                                html(&l.organism),
                                l.lowest_detected.map_or_else(|| "none".to_string(), short),
                                l.lod.map_or_else(|| "not reached".to_string(), |x| format!("{} {}", short(x), html(&l.unit)))
                            )
                        } else {
                            String::new()
                        },
                        short(v.amount),
                        html(&l.unit),
                        v.detected,
                        v.replicates,
                        short(v.mean_reads),
                        v.mean_recovery.map_or_else(|| "&ndash;".to_string(), pct)
                    ));
                }
            }
            o.push_str("</tbody></table></div></section>\n");
        }
        // Traceability.
        o.push_str("<section class=\"card\"><h2>Traceability</h2><div class=\"scroll\"><table><tbody>");
        o.push_str(&format!("<tr><th>Engine</th><td>{}</td></tr>", html(&self.engine)));
        o.push_str(&format!("<tr><th>Date (UTC)</th><td>{}</td></tr>", html(&self.date)));
        o.push_str(&format!(
            "<tr><th>Database</th><td><code>{}</code><br>SHA-256 <code>{}</code><br>k = {}, minimizer window {} (agents) / {} (host), genome window {} bp</td></tr>",
            html(&self.database.path),
            self.database.sha256,
            self.database.k,
            self.database.minimizer_window,
            self.database.host_minimizer_window,
            self.database.genome_window_bp
        ));
        o.push_str(&format!("<tr><th>Parameters</th><td><code>{}</code></td></tr>", html(&self.options.json())));
        o.push_str("</tbody></table></div>");
        o.push_str("<h3>Input files</h3><div class=\"scroll\"><table><thead><tr><th>Sample</th><th>File</th><th class=\"r\">Bytes</th><th>SHA-256</th></tr></thead><tbody>");
        for s in &self.samples {
            for f in &s.files {
                o.push_str(&format!(
                    "<tr><td>{}</td><td><code>{}</code></td><td class=\"r\">{}</td><td><code>{}</code></td></tr>",
                    html(&s.name),
                    html(&f.path),
                    thousands(f.bytes),
                    f.sha256
                ));
            }
        }
        o.push_str("</tbody></table></div>");
        o.push_str("<h3>References of the database</h3><div class=\"scroll\"><table><thead><tr><th>Name</th><th>Kind</th><th class=\"r\">Bases</th><th>Files and SHA-256</th></tr></thead><tbody>");
        for r in &self.database.references {
            let files: Vec<String> =
                r.files.iter().zip(&r.sha256).map(|(f, h)| format!("<code>{}</code><br><code class=\"muted\">{}</code>", html(f), h)).collect();
            o.push_str(&format!(
                "<tr><td>{}</td><td>{}</td><td class=\"r\">{}</td><td>{}</td></tr>",
                html(&r.name),
                r.kind.name(),
                thousands(r.bases),
                files.join("<br>")
            ));
        }
        o.push_str("</tbody></table></div></section>\n");
        if !self.warnings.is_empty() {
            o.push_str("<section class=\"card\"><h2>Warnings</h2><ul>");
            for w in &self.warnings {
                o.push_str(&format!("<li>{}</li>", html(w)));
            }
            o.push_str("</ul></section>\n");
        }
        o.push_str("</div></div>\n</body>\n</html>\n");
        o
    }
}

/// The class and label of a call badge: a process control recovered is
/// not an alarm.
fn badge(r: &OrganismResult) -> (&'static str, &'static str) {
    if r.kind == AgentKind::Spike && r.call == Call::Detected {
        ("control", "control recovered")
    } else {
        (r.call.name(), r.call.label())
    }
}

const HTML_HEAD: &str = r#"<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Adventitious agent screen</title>
<style>
:root { color-scheme: light dark;
  --page: #f9f9f7; --surface-1: #fcfcfb; --text-primary: #0b0b0b; --text-secondary: #52514e; --muted: #898781;
  --border: rgba(11,11,11,0.10); --grid: #e1e0d9;
  --detected-bg: #fbe3df; --detected-fg: #8f1d10; --background-bg: #fdf0cf; --background-fg: #6b4a00;
  --below-bg: #ecebe6; --below-fg: #52514e; --notice-bg: #eef4fc; --notice-border: #2a78d6; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --page: #0d0d0d; --surface-1: #1a1a19; --text-primary: #f0efec; --text-secondary: #c3c2b7; --muted: #898781;
  --border: rgba(255,255,255,0.10); --grid: #2c2c2a;
  --detected-bg: #4a1610; --detected-fg: #ffb4a8; --background-bg: #3d2f08; --background-fg: #f5d27a;
  --below-bg: #2c2c2a; --below-fg: #c3c2b7; --notice-bg: #13243a; --notice-border: #3987e5; } }
:root[data-theme="dark"] {
  --page: #0d0d0d; --surface-1: #1a1a19; --text-primary: #f0efec; --text-secondary: #c3c2b7; --muted: #898781;
  --border: rgba(255,255,255,0.10); --grid: #2c2c2a;
  --detected-bg: #4a1610; --detected-fg: #ffb4a8; --background-bg: #3d2f08; --background-fg: #f5d27a;
  --below-bg: #2c2c2a; --below-fg: #c3c2b7; --notice-bg: #13243a; --notice-border: #3987e5; }
html, body { margin: 0; }
body { font-family: system-ui, -apple-system, "Segoe UI", sans-serif; background: var(--page); color: var(--text-primary); }
.page { padding: 24px 16px 48px; box-sizing: border-box; }
.wrap { max-width: 1180px; margin: 0 auto; }
p, h1, .notice, li { overflow-wrap: anywhere; }
h1 { font-size: 22px; font-weight: 600; margin: 0 0 4px; }
h2 { font-size: 16px; font-weight: 600; margin: 0 0 4px; }
h3 { font-size: 14px; font-weight: 600; margin: 16px 0 6px; }
.sub, .caption { color: var(--text-secondary); font-size: 13px; margin: 0 0 12px; }
.notice { background: var(--notice-bg); border-left: 4px solid var(--notice-border); padding: 10px 14px; border-radius: 6px; font-size: 14px; margin: 0 0 16px; }
.card { background: var(--surface-1); border: 1px solid var(--border); border-radius: 10px; padding: 16px; margin-bottom: 16px; }
.scroll { overflow-x: auto; }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th, td { text-align: left; padding: 6px 8px; border-bottom: 1px solid var(--grid); vertical-align: top; overflow-wrap: normal; }
th { color: var(--text-secondary); font-weight: 600; }
tbody th { white-space: nowrap; }
.org { min-width: 150px; }
.nw { white-space: nowrap; }
.r { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
.num { font-variant-numeric: tabular-nums; color: var(--text-secondary); }
.muted { color: var(--muted); font-size: 12px; }
.why { min-width: 220px; font-size: 12px; }
code { font-family: ui-monospace, "Cascadia Code", Consolas, monospace; font-size: 12px; overflow-wrap: anywhere; word-break: break-all; }
.tag { font-size: 11px; color: var(--muted); }
.call { display: inline-block; padding: 1px 7px; border-radius: 9px; font-size: 12px; font-weight: 600; white-space: nowrap; }
.call.detected { background: var(--detected-bg); color: var(--detected-fg); }
.call.background { background: var(--background-bg); color: var(--background-fg); }
.call.below_threshold { background: var(--below-bg); color: var(--below-fg); font-weight: 500; }
.call.not_detected { color: var(--muted); font-weight: 400; }
.call.control { background: var(--notice-bg); color: var(--text-secondary); font-weight: 500; }
</style>
</head>
<body>
<div class="page"><div class="wrap">
"#;
