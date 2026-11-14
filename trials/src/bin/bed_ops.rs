#![allow(dead_code)]
// Tests of src/bed_ops.rs (subtract, complement, closest, coverage, slop),
// isolated from the rest of the project.
//
// Random inputs come from a fixed-seed generator, so every run is the same.
// Each complete output is compared, line by line, with a brute-force reference
// written here that never calls the module: quadratic loops over vectors in
// memory, and a per-base depth array instead of the sweep. The sorts run with
// an unlimited budget and with the smallest one (1 MB), which puts the runs on
// disk. The tricky semantics also have cases written by hand.
// Usage: cargo run --release --bin bed_ops
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
#[path = "../../../src/bed_ops.rs"]
mod bed_ops;

use bed_ops::*;
use externo::*;
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

fn lcg(s: &mut u64) -> u64 {
    *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    *s >> 11
}

// chr7 exists only in A and chrM only in B.
const A_NAMES: [&str; 5] = ["chr1", "chr2", "chr10", "chrX", "chr7"];
const B_NAMES: [&str; 5] = ["chr1", "chr2", "chr10", "chrX", "chrM"];
// Expected orders, written by hand: they do not use the code under test.
const LEXICAL: [&str; 6] = ["chr1", "chr10", "chr2", "chr7", "chrM", "chrX"];
const NATURAL: [&str; 6] = ["chr1", "chr2", "chr7", "chr10", "chrM", "chrX"];
// The list given to the sorter is `chrX, chr10`; the rest follows in natural order.
const LISTED: [&str; 6] = ["chrX", "chr10", "chr1", "chr2", "chr7", "chrM"];

fn orders() -> Vec<(&'static str, OrdemCromossomos, [&'static str; 6])> {
    vec![
        ("lexical", OrdemCromossomos::Lexica, LEXICAL),
        ("natural", OrdemCromossomos::Natural, NATURAL),
        ("list", OrdemCromossomos::Lista(vec!["chrX".to_string(), "chr10".to_string()]), LISTED),
    ]
}

struct Iv {
    c: &'static str,
    s: u64,
    e: u64,
    /// Columns 4 onward, with the leading tab; empty for a 3-column line.
    extra: String,
    line: String,
}

struct Shape {
    span: u64,
    max_len: u64,
    /// Some zero-length intervals.
    zero: bool,
    /// A few intervals as long as half the sequence.
    giants: bool,
}

/// Random intervals. About a third are built from the previous one, so that
/// the awkward arrangements are frequent whatever the density: bookended,
/// nested, identical, same start and a different end.
fn generate(n: usize, names: &[&'static str], sh: &Shape, tag: &str, seed: &mut u64) -> Vec<Iv> {
    let mut v: Vec<Iv> = Vec::with_capacity(n);
    for i in 0..n {
        let prev = v.last().map(|p| (p.c, p.s, p.e));
        let (c, s, e) = match (lcg(seed) % 16, prev) {
            (0 | 1, Some((c, _, pe))) => (c, pe, pe + 1 + lcg(seed) % sh.max_len),
            (2, Some((c, ps, pe))) if pe > ps => {
                let s = ps + lcg(seed) % (pe - ps);
                (c, s, s + 1 + lcg(seed) % (pe - s))
            }
            (3, Some(p)) => p,
            (4, Some((c, ps, pe))) => (c, ps, pe + 1 + lcg(seed) % sh.max_len),
            _ => {
                let c = names[(lcg(seed) % names.len() as u64) as usize];
                let s = lcg(seed) % sh.span;
                let pick = lcg(seed);
                let len = if sh.giants && pick % 200 == 0 {
                    1 + lcg(seed) % (sh.span / 2)
                } else if sh.zero && pick % 12 == 1 {
                    0
                } else {
                    1 + lcg(seed) % sh.max_len
                };
                (c, s, s + len)
            }
        };
        let extra = if i % 7 == 3 {
            String::new()
        } else {
            format!("\t{tag}{i}\t{}\t{}", lcg(seed) % 1000, if lcg(seed) % 2 == 0 { "+" } else { "-" })
        };
        let line = format!("{c}\t{s}\t{e}{extra}");
        v.push(Iv { c, s, e, extra, line });
    }
    v
}

/// Writes a BED (compression by the extension) with header lines on top and,
/// a third of the way in, lines that are not intervals and must be skipped.
fn write_bed(path: &str, ivs: &[Iv], eol: &str, header: &[&str]) {
    let mut w = Escritor::cria(path).unwrap();
    let mut put = |l: &str| {
        w.write_all(l.as_bytes()).unwrap();
        w.write_all(eol.as_bytes()).unwrap();
    };
    for h in header {
        put(h);
    }
    for (i, x) in ivs.iter().enumerate() {
        if i == ivs.len() / 3 {
            for junk in ["chr1\tnot_a_number\t5\tjunk", "one_column", "", "chr1\t7", "#chr1\t3\t4\tcommented_out"] {
                put(junk);
            }
        }
        put(&x.line);
    }
    w.termina().unwrap();
}

fn read_lines(path: &str) -> Vec<String> {
    let mut r = abre_leitura(path).unwrap_or_else(|e| panic!("{e}"));
    let mut s = String::new();
    std::io::Read::read_to_string(&mut r, &mut s).unwrap();
    // Not `lines()`: it would hide a stray carriage return.
    let mut v: Vec<String> = s.split('\n').map(str::to_string).collect();
    assert_eq!(v.pop().as_deref(), Some(""), "{path}: the last line has no line break");
    v
}

fn same(what: &str, got: &[String], exp: &[String]) {
    for i in 0..got.len().max(exp.len()) {
        if got.get(i) != exp.get(i) {
            panic!(
                "{what}: first difference at line {}\n       got: {:?}\n  expected: {:?}\n  ({} lines got, {} expected)",
                i + 1,
                got.get(i),
                exp.get(i),
                got.len(),
                exp.len()
            );
        }
    }
}

/// Hand-written lines, with spaces standing for the tabs.
fn tabs(lines: &[&str]) -> Vec<String> {
    lines.iter().map(|l| l.replace(' ', "\t")).collect()
}

fn text(lines: &[&str]) -> String {
    tabs(lines).iter().map(|l| format!("{l}\n")).collect()
}

struct Ctx {
    dir: PathBuf,
    /// Where the module is told to put its temporary files: must end empty.
    runs: PathBuf,
    tests: usize,
}

impl Ctx {
    fn p(&self, name: &str) -> String {
        self.dir.join(name).to_string_lossy().to_string()
    }

    fn op(&self, mem: usize, ordem: OrdemCromossomos) -> Opcoes {
        Opcoes { memoria_bytes: mem, dir_temp: Some(self.runs.clone()), ordem }
    }

    fn write(&self, name: &str, content: &str) -> String {
        let p = self.p(name);
        std::fs::write(&p, content).unwrap();
        p
    }

    fn no_temp(&self, what: &str) {
        assert!(std::fs::read_dir(&self.runs).unwrap().next().is_none(), "{what}: a temporary file was left behind");
    }

    /// The operation must succeed, write exactly `exp` and return its length.
    fn check(&mut self, what: &str, r: Result<usize, String>, out: &str, exp: &[String]) {
        let n = r.unwrap_or_else(|e| panic!("{what}: failed with: {e}"));
        same(what, &read_lines(out), exp);
        assert_eq!(n, exp.len(), "{what}: number of lines returned");
        self.no_temp(what);
        self.tests += 1;
    }

    /// The operation must fail, say why, and leave neither output nor temporary.
    fn fails(&mut self, what: &str, r: Result<usize, String>, out: &str, needles: &[&str]) {
        let e = match r {
            Ok(n) => panic!("{what}: should have failed, but wrote {n} lines"),
            Err(e) => e,
        };
        for needle in needles {
            assert!(e.contains(needle), "{what}: the message should mention {needle:?}, and it is: {e}");
        }
        assert!(!Path::new(out).exists(), "{what}: a misleading output was left on disk");
        self.no_temp(what);
        self.tests += 1;
    }
}

// ---------------------------------------------------------------------
//  BRUTE-FORCE REFERENCES
// ---------------------------------------------------------------------

/// What each interval of A must produce, by index of A.
struct Reference {
    subtract: Vec<Vec<String>>,
    subtract_whole: Vec<Vec<String>>,
    coverage: Vec<String>,
    closest: Vec<Vec<String>>,
}

fn reference(a: &[Iv], b: &[Iv]) -> Reference {
    // The B of each sequence in the order the sorter must give: start, end,
    // and the input order among equals (the sort below is stable).
    let mut of: HashMap<&str, Vec<&Iv>> = HashMap::new();
    for y in b {
        of.entry(y.c).or_default().push(y);
    }
    for v in of.values_mut() {
        v.sort_by_key(|y| (y.s, y.e));
    }
    let nothing: Vec<&Iv> = Vec::new();
    let mut r = Reference { subtract: Vec::new(), subtract_whole: Vec::new(), coverage: Vec::new(), closest: Vec::new() };

    for x in a {
        let bs = of.get(x.c).unwrap_or(&nothing);

        // Depth of B over each base of A, from a difference array.
        let n = (x.e - x.s) as usize;
        let mut diff = vec![0i32; n + 1];
        let mut touching = 0u64;
        for y in bs {
            let (lo, hi) = (x.s.max(y.s), x.e.min(y.e));
            if lo < hi {
                touching += 1;
                diff[(lo - x.s) as usize] += 1;
                diff[(hi - x.s) as usize] -= 1;
            }
        }
        let (mut depth, mut covered) = (0i32, 0u64);
        let mut pieces: Vec<String> = Vec::new();
        let mut free_since: Option<usize> = None;
        for p in 0..=n {
            if p < n {
                depth += diff[p];
            }
            let free = p < n && depth == 0;
            if p < n && !free {
                covered += 1;
            }
            match (free, free_since) {
                (true, None) => free_since = Some(p),
                (false, Some(from)) => {
                    pieces.push(format!("{}\t{}\t{}{}", x.c, x.s + from as u64, x.s + p as u64, x.extra));
                    free_since = None;
                }
                _ => {}
            }
        }
        if touching == 0 {
            r.subtract.push(vec![x.line.clone()]);
            r.subtract_whole.push(vec![x.line.clone()]);
        } else {
            r.subtract.push(pieces);
            r.subtract_whole.push(Vec::new());
        }
        // bedtools divides in single precision (a C `float`).
        let fraction = if n == 0 { "0.0000000".to_string() } else { format!("{:.7}", covered as f32 / n as f32) };
        r.coverage.push(format!("{}\t{touching}\t{covered}\t{n}\t{fraction}", x.line));

        // Closest: the distance to every B of the sequence, then the minimum.
        let distance = |y: &Iv| {
            if y.e <= x.s {
                x.s - y.e + 1
            } else if y.s >= x.e {
                y.s - x.e + 1
            } else {
                0
            }
        };
        match bs.iter().map(|y| distance(y)).min() {
            None => r.closest.push(vec![format!("{}\t.\t-1\t-1\t-1", x.line)]),
            Some(best) => r.closest.push(
                bs.iter().filter(|y| distance(y) == best).map(|y| format!("{}\t{}\t{best}", x.line, y.line)).collect(),
            ),
        }
    }
    r
}

/// The indexes of A in the order the sorter must give under a chromosome order.
fn sorted_order(a: &[Iv], table: &[&str; 6]) -> Vec<usize> {
    let rank = |c: &str| table.iter().position(|x| *x == c).unwrap();
    let mut ix: Vec<usize> = (0..a.len()).collect();
    ix.sort_by_key(|&i| (rank(a[i].c), a[i].s, a[i].e));
    ix
}

fn reference_complement(a: &[Iv], genome: &[(&str, u64)]) -> Vec<String> {
    let mut out = Vec::new();
    for &(name, len) in genome {
        let n = len as usize;
        let mut diff = vec![0i32; n + 1];
        for x in a.iter().filter(|x| x.c == name) {
            let (lo, hi) = (x.s.min(len) as usize, x.e.min(len) as usize);
            if lo < hi {
                diff[lo] += 1;
                diff[hi] -= 1;
            }
        }
        let mut depth = 0i32;
        let mut free_since: Option<usize> = None;
        for p in 0..=n {
            if p < n {
                depth += diff[p];
            }
            match (p < n && depth == 0, free_since) {
                (true, None) => free_since = Some(p),
                (false, Some(from)) => {
                    out.push(format!("{name}\t{from}\t{p}"));
                    free_since = None;
                }
                _ => {}
            }
        }
    }
    out
}

fn reference_slop(a: &[Iv], genome: &[(&str, u64)], left: u64, right: u64) -> Vec<String> {
    a.iter()
        .map(|x| {
            let len = genome.iter().find(|g| g.0 == x.c).unwrap().1;
            let mut s = if x.s > left { x.s - left } else { 0 };
            let mut e = x.e + right;
            if e > len {
                e = len;
            }
            if s > e {
                s = e;
            }
            format!("{}\t{s}\t{e}{}", x.c, x.extra)
        })
        .collect()
}

// ---------------------------------------------------------------------
//  RANDOM CASES
// ---------------------------------------------------------------------

fn random_cases(t: &mut Ctx) {
    let mut seed = 20261113u64;
    let out = t.p("out.bed");
    // |A|, |B|, span of the starts, longest interval, zero-length, giants
    let shapes = [
        (300usize, 300usize, 400u64, 30u64, true, false), // small and dense
        (20_000, 20_000, 3_000, 12, true, false),         // dense: deep nesting, many ties
        (25_000, 20_000, 2_000_000, 500, false, false),   // sparse: real distances
        (20_000, 22_000, 200_000, 60, true, true),        // giants that hold the window open
        (22_000, 300, 50_000, 40, false, false),          // B scarce
        (300, 22_000, 50_000, 40, true, false),           // A scarce
        (20_000, 20_000, 1_000, 1, true, false),          // unit and zero lengths: all bookends
    ];
    for (case, &(na, nb, span, max_len, zero, giants)) in shapes.iter().enumerate() {
        let sh = Shape { span, max_len, zero, giants };
        let a = generate(na, &A_NAMES, &sh, "a", &mut seed);
        let b = generate(nb, &B_NAMES, &sh, "b", &mut seed);
        // A comes in two files, the second one gzipped with CRLF line ends and
        // a header of its own; B comes in one .lz4 file.
        let a_files = vec![t.p("a1.bed"), t.p("a2.bed.gz")];
        write_bed(&a_files[0], &a[..na / 2], "\n", &["track name=test", "#chrom\tstart\tend", "browser position chr1:1-1000"]);
        write_bed(&a_files[1], &a[na / 2..], "\r\n", &["track name=second"]);
        let b_files = vec![t.p("b.bed.lz4")];
        write_bed(&b_files[0], &b, "\n", &["#the B side"]);

        // The smallest budget must really put the sort on disk.
        for (files, n) in [(&a_files, na), (&b_files, nb)] {
            let r = ordenar_arquivos(files, &t.p("probe.bed"), Formato::Bed, Criterio::Posicao, &t.op(1, OrdemCromossomos::Lexica))
                .unwrap();
            assert_eq!(r.registros, n, "case {case}: intervals read by the sorter");
            assert!(n < 20_000 || r.lotes >= 2, "case {case}: {n} intervals gave {} runs with the smallest budget", r.lotes);
        }

        let rf = reference(&a, &b);
        for (order_name, ordem, table) in orders() {
            let ix = sorted_order(&a, &table);
            let gather = |per_a: &[Vec<String>]| -> Vec<String> { ix.iter().flat_map(|&i| per_a[i].iter().cloned()).collect() };
            let exp_subtract = gather(&rf.subtract);
            let exp_whole = gather(&rf.subtract_whole);
            let exp_closest = gather(&rf.closest);
            let exp_coverage: Vec<String> = ix.iter().map(|&i| rf.coverage[i].clone()).collect();
            for mem in [usize::MAX, 1] {
                let what = |f: &str| {
                    format!("{f}, random case {case} (|A|={na} |B|={nb} span={span}), {order_name} order, budget {}", if mem == 1 { "1 MB" } else { "unlimited" })
                };
                let op = t.op(mem, ordem.clone());
                t.check(&what("subtract"), bed_subtract(&a_files, &b_files, &out, false, &op), &out, &exp_subtract);
                t.check(&what("subtract -A"), bed_subtract(&a_files, &b_files, &out, true, &op), &out, &exp_whole);
                t.check(&what("coverage"), bed_coverage(&a_files, &b_files, &out, &op), &out, &exp_coverage);
                t.check(&what("closest"), bed_closest(&a_files, &b_files, &out, &op), &out, &exp_closest);
            }
        }

        // The genome file: an order that is neither lexical nor natural, one
        // sequence shorter than the data (clamping), one of a single base, one
        // that A does not have, one of length zero; .fai columns on some lines.
        let genome: Vec<(&str, u64)> =
            vec![("chr2", span / 2), ("chrX", 1), ("chr1", span + 100), ("chrM", 16_569), ("chr10", span), ("chr7", span + 100), ("chrE", 0)];
        let mut g_text = String::from("# name\tlength\n");
        for (i, (name, len)) in genome.iter().enumerate() {
            g_text.push_str(&if i % 2 == 0 { format!("{name}\t{len}\t112\t60\t61\n") } else { format!("{name}\t{len}\n\n") });
        }
        let g_path = t.write("genome.fai", &g_text);

        let exp_complement = reference_complement(&a, &genome);
        for (order_name, ordem, _) in orders().into_iter().take(2) {
            for mem in [usize::MAX, 1] {
                let what = format!("complement, random case {case} (|A|={na} span={span}), {order_name} order, budget {}", if mem == 1 { "1 MB" } else { "unlimited" });
                let op = t.op(mem, ordem.clone());
                t.check(&what, bed_complement(&a_files, &g_path, &out, &op), &out, &exp_complement);
            }
        }
        for (left, right) in [(0u64, 0u64), (7, 0), (0, 13), (25, 250), (span * 4, span * 4)] {
            let what = format!("slop -l {left} -r {right}, random case {case} (|A|={na} span={span})");
            let op = t.op(1, OrdemCromossomos::Lexica);
            t.check(&what, bed_slop(&a_files, &g_path, &out, left, right, &op), &out, &reference_slop(&a, &genome, left, right));
        }
    }
}

// ---------------------------------------------------------------------
//  CASES WRITTEN BY HAND
// ---------------------------------------------------------------------

fn hand_written(t: &mut Ctx) {
    let out = t.p("hand.bed");
    let op = t.op(1 << 30, OrdemCromossomos::Lexica);
    let one = |p: &String| vec![p.clone()];
    let fresh = |p: &str| {
        let _ = std::fs::remove_file(p);
    };

    // ---------------- subtract
    // a1 is cut in the middle; a2 is fully covered; a3 loses two stretches and
    // is bookended on the right; a4 is bookended on both sides and stays; the
    // 3-column line loses its end; a6 comes out exactly as it came (`005`);
    // n loses overlapping, nested and bookended B, and B that starts before
    // it ends exactly where it starts; chr3 is not in B.
    let sub_a = t.write(
        "sub_a.bed",
        &text(&[
            "track",
            "chr3 1 2 a7 2",
            "chr1 100 110",
            "chr2 1000 1100 n 0",
            "chr1 70 80 a4 8",
            "chr1 50 60 a3 7",
            "#chr1 0 5 commented_out",
            "chr1 30 40 a2 6",
            "chr1 x 9",
            "chr2 005 9 a6 1",
            "chr1 10 20 a1 5",
        ]),
    );
    let sub_b = t.write(
        "sub_b.bed",
        &text(&[
            "chr1 105 200",
            "chr2 100 200",
            "chr1 80 90",
            "chr1 60 70",
            "chr1 55 57",
            "chr1 50 52",
            "chr1 25 45",
            "chr1 15 18",
            "chr2 1090 1150",
            "chr2 1030 1040",
            "chr2 1018 1019",
            "chr2 1015 1030",
            "chr2 1010 1020",
            "chr2 900 1000",
        ]),
    );
    let exp_sub = tabs(&[
        "chr1 10 15 a1 5",
        "chr1 18 20 a1 5",
        "chr1 52 55 a3 7",
        "chr1 57 60 a3 7",
        "chr1 70 80 a4 8",
        "chr1 100 105",
        "chr2 005 9 a6 1",
        "chr2 1000 1010 n 0",
        "chr2 1040 1090 n 0",
        "chr3 1 2 a7 2",
    ]);
    t.check("subtract by hand", bed_subtract(&one(&sub_a), &one(&sub_b), &out, false, &op), &out, &exp_sub);
    let exp_whole = tabs(&["chr1 70 80 a4 8", "chr2 005 9 a6 1", "chr3 1 2 a7 2"]);
    t.check("subtract -A by hand", bed_subtract(&one(&sub_a), &one(&sub_b), &out, true, &op), &out, &exp_whole);

    // ---------------- coverage
    // x: two B, one nested, all 10 bases; y: bookended to a B, nothing; z: one
    // base of three; the 3-column line: two bases of three; a zero-length A
    // inside a B; w: the last 10 bases of 100; chr9 is not in B.
    let cov_a = t.write(
        "cov_a.bed",
        &text(&["chr2 0 100 w", "chr1 50 50 zero", "chr1 30 33", "chr1 20 23 z", "chr9 1 5 lonely", "chr1 10 20 y", "chr1 0 10 x"]),
    );
    let cov_b =
        t.write("cov_b.bed", &text(&["chr1 40 60", "chr1 32 40", "chr1 30 31", "chr1 21 22", "chr1 2 4", "chr1 0 10", "chr2 90 200"]));
    let exp_cov = tabs(&[
        "chr1 0 10 x 2 10 10 1.0000000",
        "chr1 10 20 y 0 0 10 0.0000000",
        "chr1 20 23 z 1 1 3 0.3333333",
        "chr1 30 33 2 2 3 0.6666667",
        "chr1 50 50 zero 0 0 0 0.0000000",
        "chr2 0 100 w 1 10 100 0.1000000",
        "chr9 1 5 lonely 0 0 4 0.0000000",
    ]);
    t.check("coverage by hand", bed_coverage(&one(&cov_a), &one(&cov_b), &out, &op), &out, &exp_cov);

    // ---------------- closest
    // c overlaps two B (distance 0, both); b is bookended to `big` (1); for a
    // the closest on the left is `big`, which was read BEFORE `inner` and ends
    // after it; e ties left and right (three B); d has two B with the same
    // start on the right; i has two B that start together on the left, and
    // only the longer one is the closest; g is bookended on both sides (four
    // B at distance 1); f is a zero-length A between them; chr5 is not in B.
    let clo_a = t.write(
        "clo_a.bed",
        &text(&[
            "chr5 1 2 h",
            "chr2 25 25 f",
            "chr1 200 201 d",
            "chr1 120 130 a",
            "chr2 20 30 g",
            "chr1 45 47 c",
            "chr1 500 600 i",
            "chr1 199 201 e",
            "chr1 100 110 b",
        ]),
    );
    let clo_b = t.write(
        "clo_b.bed",
        &text(&[
            "chr2 30 45 r2",
            "chr1 300 400",
            "chr2 10 20 l",
            "chr1 40 50 inner",
            "chr2 30 40 r",
            "chr1 0 100 big",
            "chr1 300 350 far2",
            "chr2 5 20 l2",
        ]),
    );
    let exp_clo = tabs(&[
        "chr1 45 47 c chr1 0 100 big 0",
        "chr1 45 47 c chr1 40 50 inner 0",
        "chr1 100 110 b chr1 0 100 big 1",
        "chr1 120 130 a chr1 0 100 big 21",
        "chr1 199 201 e chr1 0 100 big 100",
        "chr1 199 201 e chr1 300 350 far2 100",
        "chr1 199 201 e chr1 300 400 100",
        "chr1 200 201 d chr1 300 350 far2 100",
        "chr1 200 201 d chr1 300 400 100",
        "chr1 500 600 i chr1 300 400 101",
        "chr2 20 30 g chr2 5 20 l2 1",
        "chr2 20 30 g chr2 10 20 l 1",
        "chr2 20 30 g chr2 30 40 r 1",
        "chr2 20 30 g chr2 30 45 r2 1",
        "chr2 25 25 f chr2 5 20 l2 6",
        "chr2 25 25 f chr2 10 20 l 6",
        "chr2 25 25 f chr2 30 40 r 6",
        "chr2 25 25 f chr2 30 45 r2 6",
        "chr5 1 2 h . -1 -1 -1",
    ]);
    t.check("closest by hand", bed_closest(&one(&clo_a), &one(&clo_b), &out, &op), &out, &exp_clo);

    // ---------------- compressed outputs, read back
    for ext in ["gz", "lz4"] {
        let packed = t.p(&format!("hand.bed.{ext}"));
        t.check(&format!("closest into .{ext}"), bed_closest(&one(&clo_a), &one(&clo_b), &packed, &op), &packed, &exp_clo);
        t.check(&format!("subtract into .{ext}"), bed_subtract(&one(&sub_a), &one(&sub_b), &packed, false, &op), &packed, &exp_sub);
        let magic: &[u8] = if ext == "gz" { &[0x1f, 0x8b] } else { &[0x04, 0x22, 0x4d, 0x18] };
        assert!(std::fs::read(&packed).unwrap().starts_with(magic), "the .{ext} output is not compressed");
    }

    // ---------------- complement
    // The genome file gives the order (chr2 before chr1). On chr1: overlapping
    // and bookended intervals make one covered block, a zero-length interval
    // does not split a gap, and an interval that runs past the end is clamped.
    // On chr2 an interval lies entirely past the end. chrM is not in A and
    // comes out whole; chrE has length zero and has nothing to write.
    let genome = t.write("hand.genome", &text(&["#name length", "chr2 100", "", "chr1 50 6 60 61", "chrM 16", "chrE 0"]));
    let com_a = t.write(
        "com_a.bed",
        &text(&["chr1 45 80", "chr2 200 300", "chr1 15 30", "chr1 40 40", "chr2 90 100", "chr1 30 35", "chr1 10 20", "chr2 0 10"]),
    );
    let exp_com = tabs(&["chr2 10 90", "chr1 0 10", "chr1 35 45", "chrM 0 16"]);
    for (name, ordem) in [("lexical", OrdemCromossomos::Lexica), ("natural", OrdemCromossomos::Natural)] {
        let r = bed_complement(&one(&com_a), &genome, &out, &t.op(1 << 30, ordem));
        t.check(&format!("complement by hand, {name} order"), r, &out, &exp_com);
    }

    // ---------------- slop
    // The input order is kept. Clamping at 0 and at the length; an interval
    // that ends past the sequence is brought back; one that starts past it
    // becomes the zero-length interval at the end; a line that does not
    // change comes out exactly as it came (`005`).
    let slop_genome = t.write("slop.genome", &text(&["chr1 100", "chr2 50"]));
    let slop_a = t.write(
        "slop_a.bed",
        &text(&["chr2 10 20 x", "chr1 30 40 y", "chr1 005 95", "chr1 90 120 z", "chr1 150 160 w", "#c", "chr2 0 50 full"]),
    );
    let slops: [(u64, u64, &[&str]); 4] = [
        (15, 40, &["chr2 0 50 x", "chr1 15 80 y", "chr1 0 100", "chr1 75 100 z", "chr1 100 100 w", "chr2 0 50 full"]),
        (0, 0, &["chr2 10 20 x", "chr1 30 40 y", "chr1 005 95", "chr1 90 100 z", "chr1 100 100 w", "chr2 0 50 full"]),
        (0, 5, &["chr2 10 25 x", "chr1 30 45 y", "chr1 5 100", "chr1 90 100 z", "chr1 100 100 w", "chr2 0 50 full"]),
        (5, 0, &["chr2 5 20 x", "chr1 25 40 y", "chr1 0 95", "chr1 85 100 z", "chr1 100 100 w", "chr2 0 50 full"]),
    ];
    for (left, right, exp) in slops {
        let r = bed_slop(&one(&slop_a), &slop_genome, &out, left, right, &op);
        t.check(&format!("slop -l {left} -r {right} by hand"), r, &out, &tabs(exp));
    }

    // ---------------- empty inputs
    let empty = t.write("empty.bed", "track name=nothing\n#only a header\n");
    let some_a = t.write("some_a.bed", &text(&["chr2 5 9 q", "chr1 10 20 p 1"]));
    let some_b = t.write("some_b.bed", &text(&["chr1 0 5"]));
    let none: Vec<String> = Vec::new();
    for (name, a) in [("empty A", one(&empty)), ("no file in A", Vec::new())] {
        for b in [&some_b, &empty] {
            t.check(&format!("subtract, {name}"), bed_subtract(&a, &one(b), &out, false, &op), &out, &none);
            t.check(&format!("subtract -A, {name}"), bed_subtract(&a, &one(b), &out, true, &op), &out, &none);
            t.check(&format!("coverage, {name}"), bed_coverage(&a, &one(b), &out, &op), &out, &none);
            t.check(&format!("closest, {name}"), bed_closest(&a, &one(b), &out, &op), &out, &none);
        }
        t.check(&format!("slop, {name}"), bed_slop(&a, &slop_genome, &out, 3, 3, &op), &out, &none);
        let whole = tabs(&["chr2 0 100", "chr1 0 50", "chrM 0 16"]);
        t.check(&format!("complement, {name}"), bed_complement(&a, &genome, &out, &op), &out, &whole);
    }
    let as_is = tabs(&["chr1 10 20 p 1", "chr2 5 9 q"]);
    t.check("subtract, empty B", bed_subtract(&one(&some_a), &one(&empty), &out, false, &op), &out, &as_is);
    t.check("subtract -A, empty B", bed_subtract(&one(&some_a), &one(&empty), &out, true, &op), &out, &as_is);
    let zeros = tabs(&["chr1 10 20 p 1 0 0 10 0.0000000", "chr2 5 9 q 0 0 4 0.0000000"]);
    t.check("coverage, empty B", bed_coverage(&one(&some_a), &one(&empty), &out, &op), &out, &zeros);
    let alone = tabs(&["chr1 10 20 p 1 . -1 -1 -1", "chr2 5 9 q . -1 -1 -1"]);
    t.check("closest, empty B", bed_closest(&one(&some_a), &one(&empty), &out, &op), &out, &alone);

    // ---------------- no sequence name in common: `chr1` against `1`
    let with_chr = t.write("with_chr.bed", &text(&["chr1 10 50 a", "chr2 5 9 b"]));
    let no_chr = t.write("no_chr.bed", &text(&["1 20 60 x", "2 1 7 y"]));
    let scaffold = t.write("scaffold.bed", &text(&["scaffold_7 20 60 x"]));
    let (ca, cb) = (one(&with_chr), one(&no_chr));
    let prefix = ["no sequence name in common", "`chr1`", "`1`", "`chr` prefix", "normalize_chromosomes"];
    fresh(&out);
    t.fails("subtract, chr1 x 1", bed_subtract(&ca, &cb, &out, false, &op), &out, &prefix);
    t.fails("subtract -A, chr1 x 1", bed_subtract(&ca, &cb, &out, true, &op), &out, &prefix);
    t.fails("coverage, chr1 x 1", bed_coverage(&ca, &cb, &out, &op), &out, &prefix);
    t.fails("closest, chr1 x 1", bed_closest(&ca, &cb, &out, &op), &out, &prefix);
    // The names differ, but not by the prefix: the hint changes.
    let version = ["no sequence name in common", "`1`", "`scaffold_7`", "version of the reference genome"];
    t.fails("subtract, other names", bed_subtract(&cb, &one(&scaffold), &out, false, &op), &out, &version);
    t.fails("coverage, other names", bed_coverage(&cb, &one(&scaffold), &out, &op), &out, &version);
    t.fails("closest, other names", bed_closest(&cb, &one(&scaffold), &out, &op), &out, &version);
    // One sequence in common is enough: no error, and the others are handled.
    let half = t.write("half.bed", &text(&["chr2 1 7 z", "3 1 2 w"]));
    let r = bed_subtract(&ca, &one(&half), &out, false, &op);
    t.check("subtract, one name in common", r, &out, &tabs(&["chr1 10 50 a", "chr2 7 9 b"]));
    let r = bed_coverage(&ca, &one(&half), &out, &op);
    t.check("coverage, one name in common", r, &out, &tabs(&["chr1 10 50 a 0 0 40 0.0000000", "chr2 5 9 b 1 2 4 0.5000000"]));
    let r = bed_closest(&ca, &one(&half), &out, &op);
    t.check("closest, one name in common", r, &out, &tabs(&["chr1 10 50 a . -1 -1 -1", "chr2 5 9 b chr2 1 7 z 0"]));

    // ---------------- sequences that the genome file does not have
    fresh(&out);
    let off = t.write("off.bed", &text(&["chr1 1 5", "chr2 2 6", "chr9 3 7 lost"]));
    let r = bed_complement(&one(&off), &genome, &out, &op);
    t.fails("complement, sequence not in the genome file", r, &out, &["`chr9`", &genome, "version of the reference genome"]);
    let r = bed_complement(&cb, &genome, &out, &op);
    t.fails("complement, 1 x chr1", r, &out, &["`1`", &genome, "`chr` prefix"]);
    let r = bed_slop(&one(&off), &slop_genome, &out, 1, 1, &op);
    t.fails("slop, sequence not in the genome file", r, &out, &["`chr9`", &slop_genome, &off, "line 3"]);

    // ---------------- a genome file that cannot be used
    let spaces = t.write("spaces.genome", "chr1 100\n");
    let twice = t.write("twice.genome", &text(&["chr1 100", "chr2 50", "chr1 100"]));
    let blank = t.write("blank.genome", "# nothing here\n\n");
    let missing = t.p("no_such.genome");
    let bad_genomes: [(&String, &[&str]); 4] = [
        (&spaces, &["line 1", "separated by a tab", "`chr1 100`"]),
        (&twice, &["line 3", "`chr1`", "twice"]),
        (&blank, &["lists no sequence"]),
        (&missing, &["Error opening", "no_such.genome"]),
    ];
    for (g, needles) in bad_genomes {
        let needles: Vec<&str> = needles.iter().copied().chain([g.as_str()]).collect();
        t.fails("complement, bad genome file", bed_complement(&one(&off), g, &out, &op), &out, &needles);
        t.fails("slop, bad genome file", bed_slop(&one(&off), g, &out, 1, 1, &op), &out, &needles);
    }

    // ---------------- an interval that starts after it ends
    // (B is read only as far as A needs it: the A below reaches past the bad line.)
    let inverted = t.write("inverted.bed", &text(&["chr1 5 9", "chr1 20 10 bad", "chr1 30 40"]));
    let reaching = one(&t.write("reaching.bed", &text(&["chr1 0 5", "chr1 100 200"])));
    let good = one(&some_b);
    let in_a = ["starts after it ends", "`chr1 20 10 bad`", "validate_file", "Input A"];
    let in_b = ["starts after it ends", "`chr1 20 10 bad`", "validate_file", "Input B"];
    let in_file = ["starts after it ends", "`chr1 20 10 bad`", "validate_file", inverted.as_str(), "line 2"];
    t.fails("subtract, inverted A", bed_subtract(&one(&inverted), &good, &out, false, &op), &out, &in_a);
    t.fails("subtract, inverted B", bed_subtract(&reaching, &one(&inverted), &out, false, &op), &out, &in_b);
    t.fails("coverage, inverted B", bed_coverage(&reaching, &one(&inverted), &out, &op), &out, &in_b);
    t.fails("closest, inverted B", bed_closest(&reaching, &one(&inverted), &out, &op), &out, &in_b);
    t.fails("complement, inverted A", bed_complement(&one(&inverted), &genome, &out, &op), &out, &in_a);
    t.fails("slop, inverted A", bed_slop(&one(&inverted), &slop_genome, &out, 1, 1, &op), &out, &in_file);

    // ---------------- an input that does not exist: an output that was
    // already there is not touched, because nothing was computed
    let absent = one(&t.p("no_such.bed"));
    let kept = t.write("kept.bed", "the previous result\n");
    let results = [
        bed_subtract(&absent, &good, &kept, false, &op),
        bed_subtract(&good, &absent, &kept, false, &op),
        bed_coverage(&good, &absent, &kept, &op),
        bed_closest(&absent, &good, &kept, &op),
        bed_complement(&absent, &genome, &kept, &op),
        bed_slop(&absent, &slop_genome, &kept, 1, 1, &op),
    ];
    for r in results {
        let e = r.expect_err("an input that does not exist must be an error");
        assert!(e.contains("no_such.bed"), "the message should name the file: {e}");
        assert_eq!(std::fs::read_to_string(&kept).unwrap(), "the previous result\n", "the previous output was destroyed");
        t.no_temp("missing input");
        t.tests += 1;
    }

    // ---------------- the output must not be one of the inputs
    let before = std::fs::read_to_string(&slop_a).unwrap();
    let results = [
        bed_slop(&one(&slop_a), &slop_genome, &slop_a, 1, 1, &op),
        bed_slop(&one(&slop_a), &slop_genome, &slop_genome, 1, 1, &op),
        bed_complement(&one(&slop_a), &slop_genome, &slop_a, &op),
        bed_subtract(&one(&slop_a), &good, &slop_a, false, &op),
        bed_coverage(&good, &one(&slop_a), &slop_a, &op),
        bed_closest(&good, &one(&slop_a), &slop_a, &op),
    ];
    for r in results {
        let e = r.expect_err("writing over an input must be refused");
        assert!(e.contains("also one of the input files"), "message: {e}");
        assert_eq!(std::fs::read_to_string(&slop_a).unwrap(), before, "the input was destroyed");
        assert_eq!(std::fs::read_to_string(&slop_genome).unwrap(), "chr1\t100\nchr2\t50\n", "the genome file was destroyed");
        t.no_temp("output over an input");
        t.tests += 1;
    }
}

fn main() {
    let dir = std::env::temp_dir().join(format!("bed_ops_test_{}", std::process::id()));
    let runs = dir.join("runs");
    std::fs::create_dir_all(&runs).unwrap();
    let mut t = Ctx { dir, runs, tests: 0 };

    hand_written(&mut t);
    random_cases(&mut t);

    std::fs::remove_dir_all(&t.dir).unwrap();
    println!("{} BED operation tests passed", t.tests);
}
