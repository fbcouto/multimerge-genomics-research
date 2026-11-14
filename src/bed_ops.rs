// =====================================================================
//  MORE INTERVAL OPERATIONS (BED)
//
//    subtract    "what is left of A once B is removed?"         (bedtools subtract)
//    complement  "which parts of the genome does A not cover?"  (bedtools complement)
//    closest     "which B is nearest to each A, and how far?"   (bedtools closest -d)
//    coverage    "how much of each A does B cover?"             (bedtools coverage)
//    slop        "grow every interval by N bases"               (bedtools slop)
//
//  Same recipe as `intervalos.rs`: each side is sorted by `externo.rs` into a
//  temporary .lz4 file -- so the inputs need not be sorted nor fit in memory --
//  and one linear sweep walks the sorted streams. `slop` is the exception: it
//  needs no order and just streams the inputs.
//
//  BED coordinates: 0-based start, exclusive end. Two intervals overlap when
//  max(starts) < min(ends), so bookended intervals (the end of one equals the
//  start of the other) do NOT overlap, and a zero-length interval overlaps
//  nothing.
//
//  NO SILENT WRONG ANSWERS. An interval whose start is greater than its end is
//  refused (bedtools stops on it too); two files that share no sequence name
//  are refused instead of producing "nothing overlaps"; a run that fails
//  leaves no partial output behind.
// =====================================================================

use crate::externo::{
    abre_leitura, ordenar_arquivos, Criterio, DirTemp, Escritor, Formato, Opcoes, OrdemCromossomos, OrdemNomes,
};
use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BinaryHeap, HashMap};
use std::io::BufRead;

fn fields(line: &str) -> Option<(&str, u64, u64)> {
    let mut c = line.split('\t');
    Some((c.next()?, c.next()?.parse().ok()?, c.next()?.parse().ok()?))
}

/// The same rule the sorter uses for BED.
fn is_header(line: &str) -> bool {
    line.starts_with('#') || line.starts_with("track") || line.starts_with("browser")
}

// ---------------------------------------------------------------------
//  INPUT AND OUTPUT
// ---------------------------------------------------------------------

/// A BED input read one interval at a time, with the current interval parsed.
/// Header lines, empty lines and lines that are not intervals are skipped.
struct Stream {
    reader: Box<dyn BufRead>,
    /// What to call this input in a message: "A", "B", or the path of a file.
    label: String,
    /// True for a file of the user, where the line number means something.
    numbered: bool,
    line_no: usize,
    /// False once the input is exhausted. The fields below are the current
    /// interval while it is true.
    live: bool,
    line: String,
    name_len: usize,
    start: u64,
    end: u64,
}

impl Stream {
    /// Opens the input and loads its first interval.
    fn open(path: &str, label: &str, numbered: bool) -> Result<Stream, String> {
        let mut s = Stream {
            reader: abre_leitura(path)?,
            label: label.to_string(),
            numbered,
            line_no: 0,
            live: false,
            line: String::new(),
            name_len: 0,
            start: 0,
            end: 0,
        };
        s.advance()?;
        Ok(s)
    }

    fn name(&self) -> &str {
        &self.line[..self.name_len]
    }

    /// The current line with new coordinates; columns 4 onward untouched.
    fn with_coordinates(&self, start: u64, end: u64) -> String {
        match self.line.splitn(4, '\t').nth(3) {
            Some(rest) => format!("{}\t{start}\t{end}\t{rest}", self.name()),
            None => format!("{}\t{start}\t{end}", self.name()),
        }
    }

    /// Loads the next interval; `live` turns false at the end of the input.
    fn advance(&mut self) -> Result<(), String> {
        loop {
            self.line.clear();
            let lidos = self.reader.read_line(&mut self.line).map_err(|e| format!("Error reading {}: {e}", self.label))?;
            if lidos == 0 {
                self.live = false;
                return Ok(());
            }
            self.line_no += 1;
            let len = self.line.trim_end_matches(['\n', '\r']).len();
            self.line.truncate(len);
            if is_header(&self.line) {
                continue;
            }
            let Some((name, start, end)) = fields(&self.line) else { continue };
            if start > end {
                let place = if self.numbered {
                    format!("{}, line {}", self.label, self.line_no)
                } else {
                    format!("Input {}", self.label)
                };
                return Err(format!(
                    "{place}: the interval `{}` starts after it ends. In a BED file the start (column 2) cannot be greater than the end (column 3). Fix or remove that line and run again (validate_file points to every line with this problem).",
                    self.line.replace('\t', " ")
                ));
            }
            (self.name_len, self.start, self.end, self.live) = (name.len(), start, end, true);
            return Ok(());
        }
    }
}

/// Sorts `paths` into a temporary file of `temp` and opens the result.
fn sorted(paths: &[String], temp: &mut DirTemp, file: &str, label: &str, op: &Opcoes) -> Result<Stream, String> {
    let path = temp.arquivo(file);
    ordenar_arquivos(paths, &path, Formato::Bed, Criterio::Posicao, op)?;
    Stream::open(&path, label, false)
}

/// The output is created -- and so emptied -- while the inputs are still
/// needed, and it is deleted when the operation fails: an output that is also
/// an input would be destroyed.
fn refuse_overwrite<'a>(saida: &str, inputs: impl Iterator<Item = &'a str>) -> Result<(), String> {
    // An output that does not exist yet cannot be one of the inputs.
    let Ok(out) = std::fs::canonicalize(saida) else { return Ok(()) };
    for input in inputs {
        if std::fs::canonicalize(input).map(|c| c == out).unwrap_or(false) {
            return Err(format!(
                "The output {saida} is also one of the input files, and writing the result would destroy it. Choose a different output path."
            ));
        }
    }
    Ok(())
}

/// Creates the output, runs `body` over it and closes it. When `body` fails
/// the partial output is deleted: a truncated file that looks complete is the
/// worst outcome.
///
/// `table`: the lines are not plain BED beyond the third column (the output
/// of `closest` and `coverage`); a Parquet output then names the extra columns
/// generically instead of calling them `name`, `score`, `strand`.
fn with_output(
    saida: &str,
    table: bool,
    body: impl FnOnce(&mut Escritor) -> Result<usize, String>,
) -> Result<usize, String> {
    let mut w = if table { Escritor::cria_tabela(saida)? } else { Escritor::cria_para(saida, Formato::Bed)? };
    let r = match body(&mut w) {
        Ok(n) => w.termina().map(|_| n),
        Err(e) => {
            drop(w);
            Err(e)
        }
    };
    if r.is_err() {
        let _ = std::fs::remove_file(saida);
    }
    r
}

/// Why two sets of sequence names do not match, and what to do about it.
fn naming_hint(ea: &str, eb: &str, doing: &str) -> String {
    if ea.starts_with("chr") != eb.starts_with("chr") {
        format!("One file uses the `chr` prefix and the other does not: rename the sequences of one of them before {doing} (normalize_chromosomes does it).")
    } else {
        "Check that both files come from the same version of the reference genome (GRCh37/hg19 and GRCh38/hg38, for instance, do not mix).".to_string()
    }
}

// ---------------------------------------------------------------------
//  THE WINDOW OVER B
// ---------------------------------------------------------------------

/// The intervals of B that can still matter to the current interval of A or
/// to a later one. A advances in sorted order and B is pulled in behind it;
/// what can no longer matter leaves. Memory follows the local density of the
/// data -- the B that start before the farthest end A has reached and end
/// after the current start of A -- not the size of the files.
///
/// Two indexes over the same intervals, so that neither entering nor leaving
/// nor listing costs a pass over the whole window: by start (to list them in
/// the order of B) and by end (to find who leaves next).
struct Window {
    b: Stream,
    order: OrdemNomes,
    /// `closest` writes the line of B; the others only need its coordinates.
    keep_lines: bool,
    /// The sequence the window is on, once `entered`.
    chrom: String,
    entered: bool,
    /// (start, arrival number) -> (end, line). B arrives in sorted order, so
    /// this iterates in the order of B.
    by_start: BTreeMap<(u64, u64), (u64, String)>,
    /// (end, start, arrival number), smallest end on top.
    by_end: BinaryHeap<Reverse<(u64, u64, u64)>>,
    arrivals: u64,
    /// Start of the last B taken on this sequence.
    last_start: Option<u64>,
    /// For the check of names at the end: the first name of B, and whether any
    /// sequence of A was found in B.
    first_b: Option<String>,
    common: bool,
}

impl Window {
    fn new(b: Stream, op: &Opcoes, keep_lines: bool) -> Window {
        Window {
            b,
            order: OrdemNomes::nova(&op.ordem),
            keep_lines,
            chrom: String::new(),
            entered: false,
            by_start: BTreeMap::new(),
            by_end: BinaryHeap::new(),
            arrivals: 0,
            last_start: None,
            first_b: None,
            common: false,
        }
    }

    fn is_on(&self, chrom: &str) -> bool {
        self.entered && self.chrom == chrom
    }

    /// A moved to a new sequence: empties the window and skips the B of the
    /// sequences that sort before it. From here on the look-ahead of B is
    /// either on `chrom` or on a later sequence.
    fn enter(&mut self, chrom: &str) -> Result<(), String> {
        self.chrom.clear();
        self.chrom.push_str(chrom);
        self.entered = true;
        self.by_start.clear();
        self.by_end.clear();
        self.last_start = None;
        while self.b.live {
            if self.first_b.is_none() {
                self.first_b = Some(self.b.name().to_string());
            }
            match self.order.compara(self.b.name(), chrom) {
                Ordering::Less => self.b.advance()?,
                Ordering::Equal => {
                    self.common = true;
                    break;
                }
                Ordering::Greater => break,
            }
        }
        Ok(())
    }

    /// True when the next B belongs to the sequence of the window.
    fn b_ahead(&self) -> bool {
        self.b.live && self.b.name() == self.chrom
    }

    /// Moves the look-ahead of B into the window.
    fn take(&mut self) -> Result<(), String> {
        let (start, end) = (self.b.start, self.b.end);
        let line = if self.keep_lines { std::mem::take(&mut self.b.line) } else { String::new() };
        self.by_start.insert((start, self.arrivals), (end, line));
        self.by_end.push(Reverse((end, start, self.arrivals)));
        self.arrivals += 1;
        self.last_start = Some(start);
        self.b.advance()
    }

    /// Removes one interval that ends at or before `s` and returns its end and
    /// line, smallest end first (in the order of B among equal ends). `None`
    /// when every interval left ends after `s`.
    fn expire(&mut self, s: u64) -> Option<(u64, String)> {
        let &Reverse((end, start, arrival)) = self.by_end.peek()?;
        if end > s {
            return None;
        }
        self.by_end.pop();
        self.by_start.remove(&(start, arrival))
    }

    /// For the sweeps that only care about overlaps: brings in every B that
    /// starts before `e` and drops the ones that end at or before `s` -- A is
    /// sorted by start, so they cannot touch this A nor any later one.
    fn reach(&mut self, s: u64, e: u64) -> Result<(), String> {
        while self.b.live && self.b.start < e && self.b_ahead() {
            self.take()?;
        }
        while self.expire(s).is_some() {}
        Ok(())
    }

    /// After `reach`: the B that overlap [s, e), in the order of B.
    fn overlaps(&self, s: u64, e: u64) -> impl Iterator<Item = (u64, u64)> + '_ {
        self.by_start
            .range(..(e, 0))
            .map(|(&(bs, _), &(be, _))| (bs, be))
            .filter(move |&(bs, be)| bs.max(s) < be.min(e))
    }

    /// SILENT NAMING ERROR, as in `bed_intersecao`. When A and B share no
    /// sequence name the sweep finds "nothing overlaps" -- formally right and
    /// almost always wrong: one file says `chr1` and the other `1`, or they
    /// come from different versions of the reference genome. The operation
    /// fails and explains. (`enter` runs for every sequence of A and stops B
    /// on that sequence when B has it, so a shared name cannot go unseen.)
    fn check_names(&self, first_a: Option<&str>, outcome: &str, doing: &str) -> Result<(), String> {
        match (first_a, &self.first_b) {
            (Some(ea), Some(eb)) if !self.common => Err(format!(
                "A and B have no sequence name in common, so {outcome}. A uses names like `{ea}`; B uses names like `{eb}`. {}",
                naming_hint(ea, eb, doing)
            )),
            _ => Ok(()),
        }
    }
}

/// What the three sweeps over A and B have in common: both sides sorted into
/// a temporary directory, and the output guarded against partial results.
/// `step` handles one interval of A and returns how many lines it wrote.
fn sweep(
    a: &[String],
    b: &[String],
    saida: &str,
    keep_lines: bool,
    table: bool,
    op: &Opcoes,
    (outcome, doing): (&str, &str),
    mut step: impl FnMut(&Stream, &mut Window, &mut Escritor) -> Result<usize, String>,
) -> Result<usize, String> {
    refuse_overwrite(saida, a.iter().chain(b).map(String::as_str))?;
    // Declared before the streams: they must close before the directory goes.
    let mut temp = DirTemp::novo(op.dir_temp.as_deref())?;
    let mut sa = sorted(a, &mut temp, "a_sorted.bed.lz4", "A", op)?;
    let sb = sorted(b, &mut temp, "b_sorted.bed.lz4", "B", op)?;
    let mut win = Window::new(sb, op, keep_lines);
    with_output(saida, table, |w| {
        let mut written = 0usize;
        let mut first_a: Option<String> = None;
        while sa.live {
            if !win.is_on(sa.name()) {
                win.enter(sa.name())?;
                if first_a.is_none() {
                    first_a = Some(sa.name().to_string());
                }
            }
            written += step(&sa, &mut win, w)?;
            sa.advance()?;
        }
        win.check_names(first_a.as_deref(), outcome, doing)?;
        Ok(written)
    })
}

// ---------------------------------------------------------------------
//  SUBTRACT
// ---------------------------------------------------------------------

/// Removes from each interval of A every stretch covered by an interval of B
/// and writes what is left: the line of A with the coordinates of each
/// remaining piece, the other columns untouched. An interval of A that B
/// covers entirely writes nothing; one that B does not touch is written as it
/// came. With `remove_whole` (`bedtools subtract -A`) an interval of A that B
/// touches at all is dropped, and the others are written as they came.
///
/// The output follows the sorted order of A. Returns the lines written.
pub fn bed_subtract(a: &[String], b: &[String], saida: &str, remove_whole: bool, op: &Opcoes) -> Result<usize, String> {
    sweep(a, b, saida, false, false, op, ("nothing was subtracted", "subtracting"), |cur, win, w| {
        let (s, e) = (cur.start, cur.end);
        win.reach(s, e)?;
        let mut written = 0usize;
        let mut touched = false;
        // Everything of A before `cursor` is settled: written or removed.
        let mut cursor = s;
        for (bs, be) in win.overlaps(s, e) {
            touched = true;
            if remove_whole {
                break;
            }
            if bs > cursor {
                w.linha(&cur.with_coordinates(cursor, bs))?;
                written += 1;
            }
            cursor = cursor.max(be);
            if cursor >= e {
                break;
            }
        }
        if !touched {
            w.linha(&cur.line)?;
            written += 1;
        } else if !remove_whole && cursor < e {
            w.linha(&cur.with_coordinates(cursor, e))?;
            written += 1;
        }
        Ok(written)
    })
}

// ---------------------------------------------------------------------
//  COVERAGE
// ---------------------------------------------------------------------

/// For each interval of A, writes its line plus four columns: how many
/// intervals of B overlap it, how many of its bases at least one B covers,
/// its length, and the covered fraction with 7 decimals. A zero-length A gets
/// `0 0 0 0.0000000`.
///
/// The output follows the sorted order of A. Returns the lines written.
pub fn bed_coverage(a: &[String], b: &[String], saida: &str, op: &Opcoes) -> Result<usize, String> {
    sweep(a, b, saida, false, true, op, ("no coverage was computed", "computing the coverage"), |cur, win, w| {
        let (s, e) = (cur.start, cur.end);
        win.reach(s, e)?;
        // The overlaps come by start: `upto` is where the covered part of A
        // ends so far, and each B only adds what lies beyond it.
        let (mut count, mut covered, mut upto) = (0u64, 0u64, s);
        for (bs, be) in win.overlaps(s, e) {
            count += 1;
            let (lo, hi) = (bs.max(upto), be.min(e));
            if hi > lo {
                covered += hi - lo;
                upto = hi;
            }
        }
        let length = e - s;
        // Single precision, as bedtools computes it (a `float`): with f64 the
        // seventh decimal differs from bedtools in about one line in a hundred
        // (467/860 is 0.5430233 in f64 and 0.5430232 in bedtools).
        let fraction = if length == 0 { 0.0f32 } else { covered as f32 / length as f32 };
        w.linha(&format!("{}\t{count}\t{covered}\t{length}\t{fraction:.7}", cur.line))?;
        Ok(1)
    })
}

// ---------------------------------------------------------------------
//  CLOSEST
// ---------------------------------------------------------------------

/// For each interval of A, finds the interval(s) of B on the same sequence at
/// the smallest distance and writes one line per closest B: the line of A,
/// the line of B and the distance. Every B tied at that distance is written,
/// in the sorted order of B. When B has nothing on the sequence, the line of
/// A is followed by `.`, `-1`, `-1` and the distance `-1`.
///
/// DISTANCE, as `bedtools closest -d` counts it. For A = [s, e) and
/// B = [bs, be):
///
///     B ends at or before the start of A (be <= s):    s - be + 1
///     B starts at or after the end of A  (bs >= e):    bs - e + 1
///     otherwise (they overlap):                        0
///
/// Bookended intervals are therefore at distance 1, not 0. (The third case
/// also takes a zero-length interval lying strictly inside the other.)
///
/// THE SWEEP. For the current A, every B falls in one of three groups:
///
///   * left behind (be <= s). A is sorted by start, so these stay to the left
///     of every later A as well, and the only ones that can ever be the
///     closest are those with the GREATEST END -- not the last ones read: B
///     can be nested, and a long interval read early may end after a short
///     one read later. They are kept in `left`, all the ties, and the rest
///     is forgotten.
///   * in the window with bs < e: they end after s, so the distance is 0.
///   * in the window with bs >= e: to the right; the closest are all those
///     with the SMALLEST START. B is read until one interval starts at or
///     after e and then until the start changes, so the whole group is in.
///
/// Whatever was not read yet starts after that group and is farther.
///
/// The output follows the sorted order of A. Returns the lines written.
pub fn bed_closest(a: &[String], b: &[String], saida: &str, op: &Opcoes) -> Result<usize, String> {
    // The lines of the B left behind with the greatest end, in the order of B.
    let mut left: Vec<String> = Vec::new();
    let mut left_end = 0u64;
    let names = ("no closest interval was found", "looking for the closest intervals");
    sweep(a, b, saida, true, true, op, names, |cur, win, w| {
        let (s, e) = (cur.start, cur.end);
        // Nothing taken on this sequence yet, so nothing of it was left
        // behind: whatever `left` holds belongs to the previous sequence.
        if win.last_start.is_none() {
            left.clear();
        }
        while win.b_ahead() && win.last_start.map_or(true, |last| last < e || last == win.b.start) {
            win.take()?;
        }
        // Intervals leave by increasing end and, among equal ends, in the
        // order of B; an interval that leaves later never comes before one
        // that left earlier with the same end. So appending keeps `left` in
        // the order of B.
        while let Some((end, line)) = win.expire(s) {
            if left.is_empty() || end > left_end {
                left.clear();
                left_end = end;
            }
            if end == left_end {
                left.push(line);
            }
        }

        let mut written = 0usize;
        for (_, (_, line)) in win.by_start.range(..(e, 0)) {
            w.linha(&format!("{}\t{line}\t0", cur.line))?;
            written += 1;
        }
        if written > 0 {
            return Ok(written);
        }
        let to_left = if left.is_empty() { None } else { Some((s - left_end).saturating_add(1)) };
        let right_start = win.by_start.range((e, 0)..).next().map(|(key, _)| key.0);
        let to_right = right_start.map(|bs| (bs - e).saturating_add(1));
        let Some(distance) = to_left.into_iter().chain(to_right).min() else {
            // Nothing left behind and nothing in the window: B has no
            // interval on this sequence.
            w.linha(&format!("{}\t.\t-1\t-1\t-1", cur.line))?;
            return Ok(1);
        };
        // A tie between the two sides: the left ones come first in B.
        if to_left == Some(distance) {
            for line in &left {
                w.linha(&format!("{}\t{line}\t{distance}", cur.line))?;
                written += 1;
            }
        }
        if to_right == Some(distance) {
            for (_, (_, line)) in win.by_start.range((e, 0)..).take_while(|(key, _)| Some(key.0) == right_start) {
                w.linha(&format!("{}\t{line}\t{distance}", cur.line))?;
                written += 1;
            }
        }
        Ok(written)
    })
}

// ---------------------------------------------------------------------
//  THE GENOME FILE
// ---------------------------------------------------------------------

/// Sequence names and lengths, in the order of the file.
struct Genome {
    names: Vec<String>,
    lengths: Vec<u64>,
    index: HashMap<String, usize>,
}

impl Genome {
    /// Reads a `.fai` index or a chrom.sizes file: the first two columns,
    /// separated by tabs, are the name and the length. Lines that start with
    /// `#` and empty lines are ignored. One line per sequence: it fits in
    /// memory whatever the size of the genome.
    fn read(path: &str) -> Result<Genome, String> {
        let mut reader = abre_leitura(path)?;
        let mut g = Genome { names: Vec::new(), lengths: Vec::new(), index: HashMap::new() };
        let mut line = String::new();
        let mut line_no = 0usize;
        loop {
            line.clear();
            if reader.read_line(&mut line).map_err(|e| format!("Error reading {path}: {e}"))? == 0 {
                break;
            }
            line_no += 1;
            let text = line.trim_end_matches(['\n', '\r']);
            if text.is_empty() || text.starts_with('#') {
                continue;
            }
            let mut cols = text.split('\t');
            let name = cols.next().unwrap_or("");
            let length = cols.next().and_then(|c| c.trim().parse::<u64>().ok());
            let Some(length) = length.filter(|_| !name.is_empty()) else {
                return Err(format!(
                    "{path}, line {line_no}: expected a sequence name and its length separated by a tab, found `{text}`. The genome file is a .fai index (samtools faidx) or a chrom.sizes file, one sequence per line."
                ));
            };
            if g.index.insert(name.to_string(), g.names.len()).is_some() {
                return Err(format!(
                    "{path}, line {line_no}: sequence `{name}` appears twice in the genome file. Keep one line per sequence."
                ));
            }
            g.names.push(name.to_string());
            g.lengths.push(length);
        }
        if g.names.is_empty() {
            return Err(format!(
                "{path} lists no sequence. The genome file needs one line per sequence, with its name and its length separated by a tab (a .fai index or a chrom.sizes file)."
            ));
        }
        Ok(g)
    }

    /// The error for a sequence of the BED that the genome file does not have.
    fn missing(&self, name: &str, path: &str, place: &str) -> String {
        format!(
            "{place} has intervals on sequence `{name}`, which is not in the genome file {path} (its sequences have names like `{}`). Add the sequence and its length to the genome file, or use the one that matches the data. {}",
            self.names[0],
            naming_hint(name, &self.names[0], "running this again")
        )
    }
}

// ---------------------------------------------------------------------
//  COMPLEMENT
// ---------------------------------------------------------------------

/// Writes `name start end` for every stretch of every sequence of the genome
/// file that no interval of A covers, in the order of the genome file. A
/// sequence without intervals in A is written whole; what A has beyond the
/// length of a sequence is ignored. A sequence of A that the genome file does
/// not have is an error. Returns the lines written.
///
/// A is sorted in the order of the genome file (whatever `op.ordem` says), so
/// one pass over the sorted A, sequence by sequence, is enough.
pub fn bed_complement(a: &[String], genome: &str, saida: &str, op: &Opcoes) -> Result<usize, String> {
    refuse_overwrite(saida, a.iter().map(String::as_str).chain([genome]))?;
    let g = Genome::read(genome)?;
    // The sequences that are not in the list sort after all the others: if A
    // has one, the sweep finds it left over at the end.
    let by_genome = Opcoes {
        memoria_bytes: op.memoria_bytes,
        dir_temp: op.dir_temp.clone(),
        ordem: OrdemCromossomos::Lista(g.names.clone()),
    };
    let mut temp = DirTemp::novo(op.dir_temp.as_deref())?;
    let mut sa = sorted(a, &mut temp, "a_sorted.bed.lz4", "A", &by_genome)?;
    with_output(saida, false, |w| {
        let mut written = 0usize;
        for (name, &length) in g.names.iter().zip(&g.lengths) {
            // Everything of the sequence before `cursor` is settled.
            let mut cursor = 0u64;
            while sa.live && sa.name() == name {
                let (s, e) = (sa.start.min(length), sa.end.min(length));
                // A zero-length interval covers no base and must not split a gap.
                if s < e {
                    if s > cursor {
                        w.linha(&format!("{name}\t{cursor}\t{s}"))?;
                        written += 1;
                    }
                    cursor = cursor.max(e);
                }
                sa.advance()?;
            }
            if cursor < length {
                w.linha(&format!("{name}\t{cursor}\t{length}"))?;
                written += 1;
            }
        }
        if sa.live {
            return Err(g.missing(sa.name(), genome, "A"));
        }
        Ok(written)
    })
}

// ---------------------------------------------------------------------
//  SLOP
// ---------------------------------------------------------------------

/// Grows every interval: the start goes down by `left` (not below 0) and the
/// end goes up by `right` (not beyond the length of the sequence in the
/// genome file). No strand handling. The other columns are untouched, and the
/// lines come out in the order they came in: nothing is sorted, the inputs
/// are streamed one after the other. A sequence that the genome file does not
/// have is an error. Returns the lines written.
///
/// An interval that already ended beyond the sequence has its end brought
/// back to the length; one that started beyond it comes out as the
/// zero-length interval at the end of the sequence, so that no line is lost
/// and no output interval is inverted.
pub fn bed_slop(a: &[String], genome: &str, saida: &str, left: u64, right: u64, op: &Opcoes) -> Result<usize, String> {
    // Nothing to sort here: no memory budget, no temporary files, no order.
    let _ = op;
    refuse_overwrite(saida, a.iter().map(String::as_str).chain([genome]))?;
    let g = Genome::read(genome)?;
    // The other operations sort first and so find a missing input before they
    // create the output. Here the same is checked by hand, so that an output
    // that already exists is not emptied for nothing.
    for path in a {
        std::fs::File::open(path).map_err(|e| format!("Error opening {path}: {e}"))?;
    }
    with_output(saida, false, |w| {
        let mut written = 0usize;
        for path in a {
            let mut input = Stream::open(path, path, true)?;
            while input.live {
                let Some(&i) = g.index.get(input.name()) else {
                    return Err(g.missing(input.name(), genome, &format!("{path} (line {})", input.line_no)));
                };
                let end = input.end.saturating_add(right).min(g.lengths[i]);
                let start = input.start.saturating_sub(left).min(end);
                if (start, end) == (input.start, input.end) {
                    w.linha(&input.line)?;
                } else {
                    w.linha(&input.with_coordinates(start, end))?;
                }
                written += 1;
                input.advance()?;
            }
        }
        Ok(written)
    })
}
