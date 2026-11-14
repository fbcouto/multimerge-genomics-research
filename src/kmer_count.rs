// =====================================================================
//  K-MER COUNTING FOR INPUTS OF ANY SIZE
//
//  FASTA and FASTQ in; out, one line `KMER<TAB>count` per distinct k-mer, in
//  increasing order of the encoded value (alphabetical: A < C < G < T), plus
//  the abundance spectrum.
//
//      input(s) --> [run 0] [run 1] ... [run n]  --> k-way merge --> output
//                    each batch sorted by Multimerge_Stable and
//                    collapsed into (k-mer, count) pairs
//
//  Memory is bounded by the budget, not by the input. The k-mers are gathered
//  in a batch; a full batch is sorted, collapsed and written to disk as a run
//  (LZ4); at the end a k-way merge joins the runs, ADDING the counts of the
//  k-mers that appear in more than one run. What fits in a single batch never
//  touches the temporary disk.
//
//  The input is read as a stream of bytes, not line by line: a chromosome
//  written on a single line of 250 million letters costs no more memory than
//  one wrapped at 60 columns.
//
//  The encoding is the one of `Kmer128` (src/kmer.rs): 2 bits per base,
//  A=00 C=01 G=10 T=11, first base in the most significant bits, k up to 63.
//  Unlike `Kmer128::encode`, nothing here turns an N into an A: any character
//  other than A, C, G, T (either case) ends the k-mer, and no k-mer contains
//  it.
// =====================================================================

use crate::externo::{abre_leitura, DirTemp, Escritor};
use crate::multimerge_core::multi_merge_sort_stable;
use lz4_flex::frame::{BlockSize, FrameDecoder, FrameEncoder, FrameInfo};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, ErrorKind, Read, Write};
use std::path::{Path, PathBuf};

/// Longest k-mer: 63 bases of 2 bits in a `u128`, the same limit as `Kmer128`.
pub const MAX_K: u8 = 63;

/// The abundance histogram has one exact entry for each count from 1 to this
/// value. Every k-mer seen more often is gathered in a single last entry,
/// whose count is `HISTOGRAM_CAP + 1` and means "more than `HISTOGRAM_CAP`".
pub const HISTOGRAM_CAP: u32 = 10_000;

/// The budget never goes below this.
const MIN_BUDGET: usize = 1 << 20;
/// Memory per buffered k-mer: 16 bytes for the value and 16 for the working
/// buffer of the stable sort.
const BYTES_PER_KMER: usize = 32;
/// One entry of a run: the k-mer (16 bytes) and its count (4), little-endian.
const RUN_ENTRY: usize = 20;
/// Buffer of a run, for writing and for reading. During the merge this (times
/// three: file, compressed block, decompressed block) is all a run costs.
const RUN_BUFFER: usize = 64 << 10;
const OUTPUT_BUFFER: usize = 64 << 10;

// ---------------------------------------------------------------------
//  THE 2-BIT ENCODING
//
//  `reverse_complement`, `canonical` and `decode` are for whoever uses the
//  module (the trials, the Python layer): the counting itself rolls the two
//  strands base by base and never calls them. Hence the `allow(dead_code)`:
//  it keeps the build quiet while nothing else in the crate uses them yet.
// ---------------------------------------------------------------------

/// The low bit of each base.
const LOW_BITS: u128 = 0x5555_5555_5555_5555_5555_5555_5555_5555;

/// Reverse complement of a k-mer in the 2-bit encoding. Only the `k` bases in
/// the low bits are looked at.
#[allow(dead_code)]
pub fn reverse_complement(bits: u128, k: u8) -> u128 {
    if k == 0 {
        return 0;
    }
    let k = u32::from(k.min(64));
    // Complement: A=00 <-> T=11 and C=01 <-> G=10 are a bitwise NOT.
    // Reverse: mirroring the 128 bits puts the bases in reverse order but also
    // swaps the two bits inside each base, so those are swapped back.
    let mirror = (!bits).reverse_bits();
    let bases = ((mirror >> 1) & LOW_BITS) | ((mirror & LOW_BITS) << 1);
    bases >> (128 - 2 * k)
}

/// The smaller of a k-mer and its reverse complement: one name for the two
/// strands of the same piece of DNA.
#[allow(dead_code)]
pub fn canonical(bits: u128, k: u8) -> u128 {
    bits.min(reverse_complement(bits, k))
}

/// Back to text: `k` letters among A, C, G, T.
#[allow(dead_code)]
pub fn decode(bits: u128, k: u8) -> String {
    let mut text = Vec::with_capacity(usize::from(k));
    push_kmer(&mut text, bits, k);
    text.into_iter().map(char::from).collect()
}

fn push_kmer(out: &mut Vec<u8>, bits: u128, k: u8) {
    let k = u32::from(k);
    // 128 bits hold 64 bases: anything asked beyond that can only be A.
    for _ in 64..k {
        out.push(b'A');
    }
    let k = k.min(64);
    if k == 0 {
        return;
    }
    let mut rest = bits << (128 - 2 * k);
    for _ in 0..k {
        out.push(b"ACGT"[(rest >> 126) as usize]);
        rest <<= 2;
    }
}

fn push_number(out: &mut Vec<u8>, mut n: u32) {
    let mut digits = [0u8; 10];
    let mut i = digits.len();
    loop {
        i -= 1;
        digits[i] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    out.extend_from_slice(&digits[i..]);
}

// ---------------------------------------------------------------------
//  OPTIONS AND RESULT
// ---------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct CountOptions {
    /// Length of the k-mers, from 1 to 63.
    pub k: u8,
    /// Count each k-mer together with its reverse complement, under the
    /// smaller of the two. A sequence and its reverse complement then give
    /// the same table.
    pub canonical: bool,
    /// Only k-mers seen at least this many times are written. 0 is the same
    /// as 1. It does not change the summary nor the histogram.
    pub min_count: u32,
    /// Memory budget for the k-mers, in bytes (at least 1 MB is used).
    pub memoria_bytes: usize,
    /// Where the temporary runs go. `None` uses the system temporary folder.
    pub dir_temp: Option<PathBuf>,
}

/// What the counting found.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct KmerCounts {
    pub k: u8,
    pub canonical: bool,
    /// Records read: `>` headers in FASTA, 4-line records in FASTQ.
    pub sequences: u64,
    /// Every character of the sequences: A, C, G, T and anything else (N,
    /// IUPAC letters, gaps). White space is not counted.
    pub bases: u64,
    /// Occurrences counted: the sum of all the counts.
    pub total_kmers: u64,
    pub distinct_kmers: u64,
    /// Distinct k-mers seen exactly once.
    pub unique_kmers: u64,
    /// The highest count. Counts saturate at `u32::MAX`.
    pub max_count: u32,
    /// Lines written to the output; 0 when there is no output.
    pub written: u64,
    /// Run files written to disk; 0 when everything fitted in memory.
    pub runs: usize,
    /// Abundance spectrum: `(count, number of distinct k-mers with that
    /// count)`, sorted by count, over ALL the distinct k-mers whatever
    /// `min_count` is. Exact for counts up to `HISTOGRAM_CAP`; one last entry
    /// with count `HISTOGRAM_CAP + 1` gathers everything above.
    pub histogram: Vec<(u32, u64)>,
}

impl KmerCounts {
    pub fn json(&self) -> String {
        let histogram: Vec<String> = self.histogram.iter().map(|(count, kmers)| format!("[{count},{kmers}]")).collect();
        format!(
            "{{\"k\":{},\"canonical\":{},\"sequences\":{},\"bases\":{},\"total_kmers\":{},\"distinct_kmers\":{},\"unique_kmers\":{},\"max_count\":{},\"written\":{},\"runs\":{},\"histogram\":[{}]}}",
            self.k,
            self.canonical,
            self.sequences,
            self.bases,
            self.total_kmers,
            self.distinct_kmers,
            self.unique_kmers,
            self.max_count,
            self.written,
            self.runs,
            histogram.join(",")
        )
    }

    /// The summary as text, ready to show.
    pub fn resumo(&self) -> String {
        let strands = if self.canonical {
            "canonical (a k-mer and its reverse complement are counted as one)"
        } else {
            "one strand (k-mers as written in the input)"
        };
        let mut o = format!("k-mer counting: k = {}, {strands}\n", self.k);
        o.push_str(&format!("  sequences read:    {}\n", thousands(self.sequences)));
        o.push_str(&format!("  bases read:        {}\n", thousands(self.bases)));
        o.push_str(&format!("  k-mers counted:    {}\n", thousands(self.total_kmers)));
        o.push_str(&format!("  distinct k-mers:   {}\n", thousands(self.distinct_kmers)));
        let share = if self.distinct_kmers > 0 {
            format!(" ({:.1}% of the distinct)", 100.0 * self.unique_kmers as f64 / self.distinct_kmers as f64)
        } else {
            String::new()
        };
        o.push_str(&format!("  seen only once:    {}{share}\n", thousands(self.unique_kmers)));
        let saturated = if self.max_count == u32::MAX { " or more (counts stop there)" } else { "" };
        o.push_str(&format!("  highest count:     {}{saturated}\n", thousands(u64::from(self.max_count))));
        o.push_str(&format!("  lines written:     {}\n", thousands(self.written)));
        if self.runs == 0 {
            o.push_str("  temporary runs:    none (everything fitted in the memory budget)\n");
        } else {
            o.push_str(&format!("  temporary runs:    {} (sorted on disk and merged)\n", self.runs));
        }
        if self.total_kmers == 0 {
            o.push_str(&format!(
                "  no k-mer was counted: no sequence has {} consecutive A, C, G or T; use a smaller k or check the input\n",
                self.k
            ));
        }
        o
    }
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

// ---------------------------------------------------------------------
//  RUNS ON DISK
// ---------------------------------------------------------------------

fn run_error(action: &str, path: &str, e: impl std::fmt::Display) -> String {
    format!(
        "Error {action} the temporary run file {path}: {e}. Check the free space of the temporary folder, or choose another one (dir_temp)"
    )
}

/// Writes one run: (k-mer, count) pairs in increasing order of k-mer, 20
/// bytes each, inside an LZ4 frame.
pub(crate) struct RunWriter {
    w: FrameEncoder<BufWriter<File>>,
    pending: Vec<u8>,
    path: String,
}

impl RunWriter {
    pub(crate) fn create(path: &str) -> Result<RunWriter, String> {
        let file = File::create(path).map_err(|e| run_error("creating", path, e))?;
        // Blocks of 64 KB: it is the block size that sets how much memory the
        // run takes when it is read back in the merge.
        let frame = FrameInfo::new().block_size(BlockSize::Max64KB);
        Ok(RunWriter {
            w: FrameEncoder::with_frame_info(frame, BufWriter::with_capacity(RUN_BUFFER, file)),
            pending: Vec::with_capacity(RUN_BUFFER + RUN_ENTRY),
            path: path.to_string(),
        })
    }

    pub(crate) fn pair(&mut self, kmer: u128, count: u32) -> Result<(), String> {
        self.pending.extend_from_slice(&kmer.to_le_bytes());
        self.pending.extend_from_slice(&count.to_le_bytes());
        if self.pending.len() >= RUN_BUFFER {
            self.write_pending()?;
        }
        Ok(())
    }

    fn write_pending(&mut self) -> Result<(), String> {
        self.w.write_all(&self.pending).map_err(|e| run_error("writing", &self.path, e))?;
        self.pending.clear();
        Ok(())
    }

    /// Closes the LZ4 frame and flushes. Without this the run is truncated.
    pub(crate) fn finish(mut self) -> Result<(), String> {
        self.write_pending()?;
        let RunWriter { w, path, .. } = self;
        let mut file = w.finish().map_err(|e| run_error("writing", &path, e))?;
        file.flush().map_err(|e| run_error("writing", &path, e))
    }
}

struct RunReader {
    r: FrameDecoder<BufReader<File>>,
    path: String,
}

impl RunReader {
    fn open(path: &str) -> Result<RunReader, String> {
        let file = File::open(path).map_err(|e| run_error("opening", path, e))?;
        Ok(RunReader { r: FrameDecoder::new(BufReader::with_capacity(RUN_BUFFER, file)), path: path.to_string() })
    }

    /// The next pair of the run; `None` at its end.
    fn next(&mut self) -> Result<Option<(u128, u32)>, String> {
        let at_end = self.r.fill_buf().map_err(|e| run_error("reading", &self.path, e))?.is_empty();
        if at_end {
            return Ok(None);
        }
        let mut kmer = [0u8; 16];
        let mut count = [0u8; 4];
        self.r
            .read_exact(&mut kmer)
            .and_then(|_| self.r.read_exact(&mut count))
            .map_err(|e| run_error("reading", &self.path, e))?;
        Ok(Some((u128::from_le_bytes(kmer), u32::from_le_bytes(count))))
    }
}

/// Calls `each(k-mer, count)` for every distinct value of a sorted slice.
fn collapse(sorted: &[u128], mut each: impl FnMut(u128, u32) -> Result<(), String>) -> Result<(), String> {
    let mut i = 0;
    while i < sorted.len() {
        let kmer = sorted[i];
        let mut j = i + 1;
        while j < sorted.len() && sorted[j] == kmer {
            j += 1;
        }
        each(kmer, u32::try_from(j - i).unwrap_or(u32::MAX))?;
        i = j;
    }
    Ok(())
}

// ---------------------------------------------------------------------
//  WHERE THE COUNTED K-MERS GO: STATISTICS AND OUTPUT
// ---------------------------------------------------------------------

fn output_error(path: &str, e: impl std::fmt::Display) -> String {
    format!("Error writing the output {path}: {e}. Check the free space and the permission to write there")
}

/// Receives each distinct k-mer once, in increasing order, with its final
/// count.
pub(crate) struct Collector {
    k: u8,
    min_count: u32,
    out: Option<Escritor>,
    out_path: String,
    lines: Vec<u8>,
    distinct: u64,
    unique: u64,
    max_count: u32,
    written: u64,
    /// `spectrum[c]`: distinct k-mers with count `c`; the last position holds
    /// every count above `HISTOGRAM_CAP`.
    spectrum: Vec<u64>,
}

impl Collector {
    /// With `saida`, creates the output file (the extension chooses the
    /// compression: `.gz`, `.lz4` or none).
    pub(crate) fn new(k: u8, min_count: u32, saida: Option<&str>) -> Result<Collector, String> {
        let out = match saida {
            Some(path) => Some(
                Escritor::cria(path)
                    .map_err(|e| format!("{e}. Check that the output folder exists and that you can write to it"))?,
            ),
            None => None,
        };
        Ok(Collector {
            k,
            min_count: min_count.max(1),
            out,
            out_path: saida.unwrap_or_default().to_string(),
            lines: Vec::with_capacity(OUTPUT_BUFFER + 128),
            distinct: 0,
            unique: 0,
            max_count: 0,
            written: 0,
            spectrum: vec![0; HISTOGRAM_CAP as usize + 2],
        })
    }

    pub(crate) fn receive(&mut self, kmer: u128, count: u32) -> Result<(), String> {
        self.distinct += 1;
        self.unique += u64::from(count == 1);
        self.max_count = self.max_count.max(count);
        self.spectrum[count.min(HISTOGRAM_CAP + 1) as usize] += 1;
        if count >= self.min_count && self.out.is_some() {
            push_kmer(&mut self.lines, kmer, self.k);
            self.lines.push(b'\t');
            push_number(&mut self.lines, count);
            self.lines.push(b'\n');
            self.written += 1;
            if self.lines.len() >= OUTPUT_BUFFER {
                self.write_lines()?;
            }
        }
        Ok(())
    }

    fn write_lines(&mut self) -> Result<(), String> {
        if let Some(w) = self.out.as_mut() {
            w.write_all(&self.lines).map_err(|e| output_error(&self.out_path, e))?;
        }
        self.lines.clear();
        Ok(())
    }

    /// Closes the output and hands over what was collected. The fields that
    /// come from the reading (`canonical`, `sequences`, `bases`,
    /// `total_kmers`, `runs`) are left at zero for the caller to fill.
    pub(crate) fn finish(mut self) -> Result<KmerCounts, String> {
        self.write_lines()?;
        if let Some(w) = self.out.take() {
            w.termina().map_err(|e| output_error(&self.out_path, e))?;
        }
        let histogram =
            self.spectrum.iter().enumerate().filter(|(_, kmers)| **kmers > 0).map(|(count, &kmers)| (count as u32, kmers)).collect();
        Ok(KmerCounts {
            k: self.k,
            distinct_kmers: self.distinct,
            unique_kmers: self.unique,
            max_count: self.max_count,
            written: self.written,
            histogram,
            ..KmerCounts::default()
        })
    }
}

/// K-way merge of the runs into the collector: a heap of (k-mer, run), the
/// count of each head kept beside it. The counts of a k-mer present in more
/// than one run are ADDED, saturating at `u32::MAX`. Memory: one small buffer
/// per run.
pub(crate) fn merge_runs(paths: &[String], collector: &mut Collector) -> Result<(), String> {
    let mut readers: Vec<RunReader> = Vec::with_capacity(paths.len());
    let mut counts = vec![0u32; paths.len()];
    let mut heap: BinaryHeap<Reverse<(u128, usize)>> = BinaryHeap::with_capacity(paths.len());
    for (i, path) in paths.iter().enumerate() {
        let mut reader = RunReader::open(path)?;
        if let Some((kmer, count)) = reader.next()? {
            counts[i] = count;
            heap.push(Reverse((kmer, i)));
        }
        readers.push(reader);
    }

    while let Some(Reverse((kmer, first))) = heap.pop() {
        let mut total = 0u32;
        let mut run = first;
        loop {
            total = total.saturating_add(counts[run]);
            if let Some((next, count)) = readers[run].next()? {
                counts[run] = count;
                heap.push(Reverse((next, run)));
            }
            match heap.peek() {
                Some(&Reverse((same, other))) if same == kmer => {
                    heap.pop();
                    run = other;
                }
                _ => break,
            }
        }
        collector.receive(kmer, total)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------
//  FROM BASES TO BATCHES OF K-MERS
// ---------------------------------------------------------------------

/// What a byte of a sequence line is: 0 to 3 for a base (A, C, G, T, either
/// case), `OTHER` for anything else that takes a place in the sequence, and
/// `BLANK` for white space, which takes none.
const OTHER: u8 = 4;
const BLANK: u8 = 5;
const CLASS: [u8; 256] = classes();

const fn classes() -> [u8; 256] {
    let mut t = [OTHER; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t[b' ' as usize] = BLANK;
    t[b'\t' as usize] = BLANK;
    t[b'\r' as usize] = BLANK;
    t[b'\n' as usize] = BLANK;
    t[0x0B] = BLANK;
    t[0x0C] = BLANK;
    t
}

/// The rolling encoder and the batch it fills. A full batch becomes a run on
/// disk.
struct Batches<'a> {
    k: u32,
    canonical: bool,
    /// The 2k low bits.
    mask: u128,
    /// Where the complement of the newest base enters the reverse strand: the
    /// position of the first base, 2 * (k - 1).
    top: u32,
    forward: u128,
    reverse: u128,
    /// Valid bases in a row up to here, counted up to k.
    valid: u32,
    batch: Vec<u128>,
    /// K-mers a batch may hold.
    limit: usize,
    dir_temp: Option<&'a Path>,
    temp: Option<DirTemp>,
    runs: Vec<String>,
    bases: u64,
    total: u64,
}

impl<'a> Batches<'a> {
    fn new(opt: &'a CountOptions) -> Batches<'a> {
        let k = u32::from(opt.k);
        Batches {
            k,
            canonical: opt.canonical,
            mask: (1u128 << (2 * k)) - 1,
            top: 2 * (k - 1),
            forward: 0,
            reverse: 0,
            valid: 0,
            batch: Vec::new(),
            limit: (opt.memoria_bytes.max(MIN_BUDGET) / BYTES_PER_KMER).max(1),
            dir_temp: opt.dir_temp.as_deref(),
            temp: None,
            runs: Vec::new(),
            bases: 0,
            total: 0,
        }
    }

    /// A k-mer never goes from one record into the next.
    fn new_record(&mut self) {
        self.valid = 0;
    }

    /// A piece of sequence, of any length. It continues the piece before it:
    /// a k-mer may start in one and end in another (line breaks, and the
    /// blocks in which the file is read, do not interrupt a sequence).
    fn sequence(&mut self, piece: &[u8]) -> Result<(), String> {
        // Each byte gives at most one k-mer: with this, no `push` below ever
        // has to grow the batch.
        self.reserve(piece.len().min(self.limit - self.batch.len()));
        let (k, mask, top, canonical) = (self.k, self.mask, self.top, self.canonical);
        let (mut forward, mut reverse, mut valid) = (self.forward, self.reverse, self.valid);
        let (mut kmers, mut blanks) = (0u64, 0u64);
        for &byte in piece {
            let class = CLASS[usize::from(byte)];
            if class < OTHER {
                // The new base enters the forward strand at the end and the
                // reverse strand, complemented (3 - code), at the start.
                // After k valid bases nothing older is left in either.
                forward = ((forward << 2) | u128::from(class)) & mask;
                reverse = (reverse >> 2) | (u128::from(3 - class) << top);
                if valid < k {
                    valid += 1;
                }
                if valid == k {
                    // A full batch goes to disk only when one more k-mer
                    // arrives: what fits exactly in a batch is never a run.
                    if self.batch.len() >= self.limit {
                        self.flush()?;
                    }
                    self.batch.push(if canonical { forward.min(reverse) } else { forward });
                    kmers += 1;
                }
            } else if class == OTHER {
                valid = 0;
            } else {
                blanks += 1;
            }
        }
        (self.forward, self.reverse, self.valid) = (forward, reverse, valid);
        self.bases += piece.len() as u64 - blanks;
        self.total += kmers;
        Ok(())
    }

    /// Makes room for `more` k-mers without letting the capacity pass the
    /// budget: a plain `push` doubles it, and could take twice the memory
    /// that was allowed.
    fn reserve(&mut self, more: usize) {
        let needed = self.batch.len() + more; // never above `limit`
        if needed > self.batch.capacity() {
            let target = needed.max(self.batch.capacity().saturating_mul(2)).min(self.limit);
            self.batch.reserve_exact(target - self.batch.len());
        }
    }

    /// Sorts the batch, collapses it into (k-mer, count) pairs and writes it
    /// as a run.
    fn flush(&mut self) -> Result<(), String> {
        multi_merge_sort_stable(&mut self.batch);
        let temp = match &mut self.temp {
            Some(temp) => temp,
            None => {
                let novo = DirTemp::novo(self.dir_temp).map_err(|e| {
                    format!("{e}. The temporary folder must exist and be writable: choose another one (dir_temp)")
                })?;
                self.temp.insert(novo)
            }
        };
        // The name enters the temporary folder before anything is written: if
        // the writing fails halfway, the partial file goes away with the rest.
        let path = temp.arquivo(&format!("kmers_{}.lz4", self.runs.len()));
        self.runs.push(path.clone());
        let mut run = RunWriter::create(&path)?;
        collapse(&self.batch, |kmer, count| run.pair(kmer, count))?;
        run.finish()?;
        self.batch.clear();
        Ok(())
    }
}

// ---------------------------------------------------------------------
//  READING FASTA AND FASTQ
// ---------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Fasta,
    Fastq,
}

/// What the line being read is.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Line {
    /// Nothing but white space so far: its first visible character will tell.
    Start,
    /// A header, a comment or the `+` of a FASTQ: nothing in it is used.
    Skip,
    Sequence,
    Quality,
}

/// The state of the reading of one file. The file arrives in pieces that
/// never contain a line break, and a line may come in several pieces.
struct Reading<'a> {
    path: &'a str,
    kind: Option<Kind>,
    line: Line,
    /// Number of the line being read, from 1.
    line_no: u64,
    /// FASTQ: which of the 4 lines of a record this is, 0 to 3.
    part: u8,
    /// FASTQ: inside a record whose quality line has not ended.
    open: bool,
    /// FASTQ: the line where the current record starts.
    record_line: u64,
    seq_len: u64,
    qual_len: u64,
}

impl Reading<'_> {
    fn piece(&mut self, mut piece: &[u8], batches: &mut Batches, sequences: &mut u64) -> Result<(), String> {
        if self.line == Line::Start {
            // A byte-order mark, which Windows editors put in front of the
            // text, is not part of the first line.
            let unknown = self.kind.is_none();
            let blank = |b: u8| CLASS[usize::from(b)] == BLANK || (unknown && matches!(b, 0xEF | 0xBB | 0xBF));
            let skipped = piece.iter().position(|&b| !blank(b)).unwrap_or(piece.len());
            piece = &piece[skipped..];
            let Some(&first) = piece.first() else {
                return Ok(());
            };
            self.first_character(first, piece, batches, sequences)?;
        }
        match self.line {
            Line::Sequence => {
                let before = batches.bases;
                batches.sequence(piece)?;
                self.seq_len += batches.bases - before;
            }
            Line::Quality => self.qual_len += piece.iter().filter(|&&b| CLASS[usize::from(b)] != BLANK).count() as u64,
            Line::Start | Line::Skip => {}
        }
        Ok(())
    }

    /// The first visible character of a line says what the line is -- and,
    /// on the first line of the file, what the file is.
    fn first_character(&mut self, first: u8, piece: &[u8], batches: &mut Batches, sequences: &mut u64) -> Result<(), String> {
        let kind = match self.kind {
            Some(kind) => kind,
            None => {
                let kind = match first {
                    b'>' => Kind::Fasta,
                    b'@' => Kind::Fastq,
                    _ => return Err(self.not_sequences(piece)),
                };
                self.kind = Some(kind);
                kind
            }
        };
        self.line = match (kind, self.part, first) {
            (Kind::Fasta, _, b'>') => {
                *sequences += 1;
                batches.new_record();
                Line::Skip
            }
            // The old FASTA comment line.
            (Kind::Fasta, _, b';') => Line::Skip,
            (Kind::Fasta, _, _) => Line::Sequence,
            (Kind::Fastq, 2, b'+') => Line::Skip,
            (Kind::Fastq, 2, _) => {
                return Err(format!(
                    "{}, line {}: expected the '+' line of the FASTQ record that starts at line {}, found \"{}\". A FASTQ record has exactly 4 lines (@name, sequence, +, qualities); if the sequences of this file are wrapped over several lines, rewrite it with one line per sequence (for instance `seqtk seq -l0`)",
                    self.path,
                    self.line_no,
                    self.record_line,
                    sample(piece)
                ));
            }
            (Kind::Fastq, _, b'@') => {
                *sequences += 1;
                batches.new_record();
                self.open = true;
                self.record_line = self.line_no;
                self.seq_len = 0;
                self.qual_len = 0;
                Line::Skip
            }
            (Kind::Fastq, _, _) => {
                return Err(format!(
                    "{}, line {}: expected a FASTQ record, which starts with '@', found \"{}\". The file is damaged or is not a 4-line FASTQ: check it with the validator, or get it again",
                    self.path,
                    self.line_no,
                    sample(piece)
                ));
            }
        };
        Ok(())
    }

    fn end_of_line(&mut self) -> Result<(), String> {
        (self.part, self.line) = match self.kind {
            Some(Kind::Fastq) => match (self.part, self.line) {
                // An empty line between two records.
                (0, Line::Start) => (0, Line::Start),
                (0, _) => (1, Line::Sequence),
                (1, _) => (2, Line::Start),
                (2, Line::Start) => {
                    return Err(format!(
                        "{}, line {}: empty line where the '+' line of the FASTQ record that starts at line {} should be. A FASTQ record has exactly 4 lines (@name, sequence, +, qualities): remove the empty line, or get the file again",
                        self.path, self.line_no, self.record_line
                    ));
                }
                (2, _) => (3, Line::Quality),
                _ => {
                    self.check_qualities(false)?;
                    self.open = false;
                    (0, Line::Start)
                }
            },
            _ => (0, Line::Start),
        };
        self.line_no += 1;
        Ok(())
    }

    /// One quality character per base: it is what tells a whole record from
    /// one cut in the middle of its last line.
    fn check_qualities(&self, at_end_of_file: bool) -> Result<(), String> {
        if self.seq_len == self.qual_len {
            return Ok(());
        }
        Err(if at_end_of_file {
            self.truncated(&format!("its quality line has {} characters for {} bases", self.qual_len, self.seq_len))
        } else {
            format!(
                "{}, line {}: the quality line has {} characters, but the sequence of the record (line {}) has {} bases, and each base needs one quality character. The file is damaged or is not a 4-line FASTQ: check it with the validator, or get it again",
                self.path,
                self.line_no,
                self.qual_len,
                self.line_no.saturating_sub(2),
                self.seq_len
            )
        })
    }

    fn end_of_file(&self) -> Result<(), String> {
        if !self.open {
            return Ok(());
        }
        match (self.part, self.line) {
            // The last line of a file does not need a line break.
            (3, _) => self.check_qualities(true),
            (2, Line::Skip) => Err(self.truncated("its quality line is missing")),
            (2, _) => Err(self.truncated("its '+' line and its quality line are missing")),
            _ => Err(self.truncated("its sequence may be incomplete, and its '+' line and its quality line are missing")),
        }
    }

    fn truncated(&self, what: &str) -> String {
        format!(
            "{}: truncated FASTQ: the file ends in the middle of the record that starts at line {} ({what}). A FASTQ record has 4 lines (@name, sequence, +, qualities). The file was probably cut during a copy or a download: get it again, or delete the incomplete last record",
            self.path, self.record_line
        )
    }

    fn not_sequences(&self, piece: &[u8]) -> String {
        let hint = if piece.starts_with(&[0x1F, 0x8B]) {
            " The content is gzip-compressed: rename the file so that its name ends in .gz"
        } else if piece.starts_with(&[0x04, 0x22, 0x4D, 0x18]) {
            " The content is LZ4-compressed: rename the file so that its name ends in .lz4"
        } else {
            " K-mer counting reads FASTA and FASTQ files (plain, .gz or .lz4): give it the sequences, not a table of variants, intervals or alignments"
        };
        format!(
            "{}: this is not a FASTA nor a FASTQ file: line {} starts with \"{}\", and the first line of a FASTA file starts with '>' and of a FASTQ file with '@'.{hint}",
            self.path,
            self.line_no,
            sample(piece)
        )
    }
}

/// The beginning of a line, safe to show in a message.
fn sample(piece: &[u8]) -> String {
    let mut o = String::new();
    for &b in piece.iter().take(30) {
        match b {
            b'"' | b'\\' => {
                o.push('\\');
                o.push(char::from(b));
            }
            0x20..=0x7E => o.push(char::from(b)),
            _ => o.push_str(&format!("\\x{b:02X}")),
        }
    }
    if piece.len() > 30 {
        o.push_str("...");
    }
    o
}

fn read_file(path: &str, batches: &mut Batches, sequences: &mut u64) -> Result<(), String> {
    let mut input = abre_leitura(path)?;
    let mut reading =
        Reading { path, kind: None, line: Line::Start, line_no: 1, part: 0, open: false, record_line: 0, seq_len: 0, qual_len: 0 };
    loop {
        let block = match input.fill_buf() {
            Ok(block) => block,
            Err(e) if e.kind() == ErrorKind::Interrupted => continue,
            Err(e) => {
                return Err(format!(
                    "Error reading {path} near line {}: {e}. If the file is compressed it may be damaged or incomplete: test it (`gzip -t`) or get it again",
                    reading.line_no
                ));
            }
        };
        if block.is_empty() {
            break;
        }
        let size = block.len();
        let mut rest = block;
        while !rest.is_empty() {
            match rest.iter().position(|&b| b == b'\n') {
                Some(end) => {
                    reading.piece(&rest[..end], batches, sequences)?;
                    reading.end_of_line()?;
                    rest = &rest[end + 1..];
                }
                None => {
                    reading.piece(rest, batches, sequences)?;
                    rest = &[];
                }
            }
        }
        input.consume(size);
    }
    reading.end_of_file()
}

// ---------------------------------------------------------------------
//  ENTRY POINT
// ---------------------------------------------------------------------

/// What the reading side of the counting knows.
struct ReadTotals {
    sequences: u64,
    bases: u64,
    total: u64,
    runs: usize,
}

/// Reads every input, in sequence, and delivers the distinct k-mers to the
/// collector in increasing order.
fn count_into(caminhos: &[String], opt: &CountOptions, collector: &mut Collector) -> Result<ReadTotals, String> {
    let mut batches = Batches::new(opt);
    let mut sequences = 0u64;
    for caminho in caminhos {
        read_file(caminho, &mut batches, &mut sequences)?;
    }
    if batches.runs.is_empty() {
        // Everything fitted: sort and deliver, without touching the
        // temporary disk.
        multi_merge_sort_stable(&mut batches.batch);
        collapse(&batches.batch, |kmer, count| collector.receive(kmer, count))?;
    } else {
        if !batches.batch.is_empty() {
            batches.flush()?;
        }
        // Lets go of the last batch before the merge: from here on the
        // memory is one small buffer per run.
        batches.batch = Vec::new();
        merge_runs(&batches.runs, collector)?;
    }
    Ok(ReadTotals { sequences, bases: batches.bases, total: batches.total, runs: batches.runs.len() })
}

fn same_file(a: &str, b: &str) -> bool {
    match (std::fs::canonicalize(a), std::fs::canonicalize(b)) {
        (Ok(x), Ok(y)) => x == y,
        _ => a == b,
    }
}

/// Counts the k-mers of FASTA and FASTQ files (plain, `.gz` or `.lz4`, mixed
/// at will) of any size.
///
/// With `saida`, writes one line `KMER<TAB>count` for each distinct k-mer
/// seen at least `min_count` times, in increasing order of the encoded value
/// (alphabetical); the extension of `saida` chooses the compression. Without
/// it, only the summary is computed.
///
/// - The kind of each file is told by its first non-empty line: `>` FASTA,
///   `@` FASTQ.
/// - A k-mer may span the lines of a FASTA record, never two records.
/// - Upper and lower case are the same. Any character other than A, C, G, T
///   ends the k-mer; white space inside a sequence is ignored.
/// - A file with no sequence at all gives zero counts, not an error.
/// - If the counting fails, the output file is removed: half a table would
///   pass for a whole one.
pub fn count_kmers(caminhos: &[String], saida: Option<&str>, opt: &CountOptions) -> Result<KmerCounts, String> {
    if opt.k == 0 || opt.k > MAX_K {
        return Err(format!(
            "Invalid k-mer length: k = {}. Use a k from 1 to {MAX_K} (a k-mer is kept in 128 bits, 2 per base)",
            opt.k
        ));
    }
    // Before any work, and before the output is touched: an input that is not
    // there, or that the output would overwrite.
    for caminho in caminhos {
        File::open(caminho).map_err(|e| {
            // The two everyday reasons in English, whatever the language of
            // the system.
            let reason = match e.kind() {
                ErrorKind::NotFound => "the file does not exist".to_string(),
                ErrorKind::PermissionDenied => "permission denied".to_string(),
                _ => e.to_string(),
            };
            format!("Cannot open the input file {caminho}: {reason}. Check the path and that you are allowed to read the file")
        })?;
        if let Some(saida) = saida {
            if same_file(caminho, saida) {
                return Err(format!(
                    "The output {saida} is also one of the input files, and would be erased before being read. Choose another name for the output"
                ));
            }
        }
    }

    let mut collector = Collector::new(opt.k, opt.min_count, saida)?;
    let result = match count_into(caminhos, opt, &mut collector) {
        Ok(read) => collector.finish().map(|table| KmerCounts {
            canonical: opt.canonical,
            sequences: read.sequences,
            bases: read.bases,
            total_kmers: read.total,
            runs: read.runs,
            ..table
        }),
        Err(e) => {
            // Closes the output, so that it can be removed.
            drop(collector);
            Err(e)
        }
    };
    if result.is_err() {
        if let Some(saida) = saida {
            let _ = std::fs::remove_file(saida);
        }
    }
    result
}
