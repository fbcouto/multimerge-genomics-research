#![allow(dead_code)]
// Trials of src/surto.rs: distances and trees between genomes, the mutations
// of a viral genome against its reference, and in-silico PCR.
//
// Every result is compared with an oracle written here in the plainest way:
// the Jaccard by sets of k-mer strings (and their MurmurHash3, for scaled >
// 1); the trees against hand-worked examples (Wikipedia's neighbor-joining
// and UPGMA examples, a 4-taxon non-additive matrix worked by hand) and
// random additive trees, which neighbor joining must give back exactly;
// the variants against synthetic genomes with planted substitutions,
// indels (their equivalent places found by brute force), N runs, missing
// ends, overhangs and a reverse complement; the amino acids by translating
// both sequences codon by codon (with a genetic code typed here), and by
// hand for the indel cases (S:H69-/V70-, E156-/F157-/R158G, 214:EPE,
// 16:MPLF, frameshifts); the PCR products by a naive search of every
// primer site on both strands.
//
// Usage: cargo run --release --bin surto
//        cargo run --release --bin surto -- real <research_data> <work dir>
//               [nextclade.tsv] [sourmash_kleb.csv] [sourmash_sars.csv]
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
#[path = "../../../src/microbio.rs"]
#[allow(dead_code)]
mod microbio;
#[path = "../../../src/annotation.rs"]
#[allow(dead_code)]
mod annotation;
#[path = "../../../src/surto.rs"]
mod surto;

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::path::PathBuf;
use surto::*;

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

    /// Two sets of labels must be equal; prints the difference.
    fn same_set(&mut self, what: &str, got: &[String], want: &[String]) -> bool {
        let g: BTreeSet<&String> = got.iter().collect();
        let w: BTreeSet<&String> = want.iter().collect();
        let extra: Vec<&&String> = g.difference(&w).collect();
        let missing: Vec<&&String> = w.difference(&g).collect();
        self.check(extra.is_empty() && missing.is_empty() && got.len() == want.len(), format!("{what}: extra {extra:?}, missing {missing:?} (got {}, want {})", got.len(), want.len()))
    }
}

fn near(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * a.abs().max(b.abs()).max(1.0)
}

// ---------------------------------------------------------------------
//  Random DNA
// ---------------------------------------------------------------------

fn lcg(s: &mut u64) -> u64 {
    *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    *s >> 33
}

fn rnd(s: &mut u64, n: usize) -> usize {
    (lcg(s) % n as u64) as usize
}

fn random_dna(n: usize, s: &mut u64) -> Vec<u8> {
    (0..n).map(|_| b"ACGT"[rnd(s, 4)]).collect()
}

fn other_base(b: u8, s: &mut u64) -> u8 {
    loop {
        let c = b"ACGT"[rnd(s, 4)];
        if c != b.to_ascii_uppercase() {
            return c;
        }
    }
}

fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        b'R' => b'Y',
        b'Y' => b'R',
        o => o,
    }
}

fn rc(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| comp(b)).collect()
}

fn fasta(records: &[(String, Vec<u8>)]) -> String {
    let mut o = String::new();
    for (name, seq) in records {
        o.push_str(&format!(">{name} test record\n"));
        for chunk in seq.chunks(70) {
            o.push_str(std::str::from_utf8(chunk).unwrap());
            o.push('\n');
        }
    }
    o
}

fn st(x: &[u8]) -> String {
    String::from_utf8_lossy(x).to_string()
}

// ---------------------------------------------------------------------
//  A genetic code typed here (the oracle does not use annotation.rs)
// ---------------------------------------------------------------------

const AA: &[u8; 64] = b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn iupac(b: u8) -> &'static [u8] {
    match b {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        _ => b"",
    }
}

fn oracle_codon(c: &[u8]) -> u8 {
    if c.iter().any(|&b| b == b'N') {
        return b'X';
    }
    let mut set = BTreeSet::new();
    for &x in iupac(c[0]) {
        for &y in iupac(c[1]) {
            for &z in iupac(c[2]) {
                set.insert(AA[base_index(x).unwrap() * 16 + base_index(y).unwrap() * 4 + base_index(z).unwrap()]);
            }
        }
    }
    if set.len() == 1 {
        *set.iter().next().unwrap()
    } else {
        b'X'
    }
}

fn oracle_translate(seq: &[u8]) -> Vec<u8> {
    seq.chunks_exact(3).map(oracle_codon).collect()
}

/// A random sense codon (no stop).
fn sense_codon(s: &mut u64) -> [u8; 3] {
    loop {
        let c = [b"ACGT"[rnd(s, 4)], b"ACGT"[rnd(s, 4)], b"ACGT"[rnd(s, 4)]];
        if oracle_codon(&c) != b'*' {
            return c;
        }
    }
}

// =====================================================================
//  1. DISTANCES
// =====================================================================

/// The canonical k-mers of the sequences, as text (k-mers with anything but
/// A, C, G, T skipped).
fn kmer_text_set(seqs: &[&[u8]], k: usize) -> HashSet<Vec<u8>> {
    let mut set = HashSet::new();
    for seq in seqs {
        let up: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
        if up.len() < k {
            continue;
        }
        for i in 0..=up.len() - k {
            let w = &up[i..i + k];
            if w.iter().any(|b| !b"ACGT".contains(b)) {
                continue;
            }
            let r = rc(w);
            set.insert(if w <= &r[..] { w.to_vec() } else { r });
        }
    }
    set
}

fn oracle_hash_set(seqs: &[&[u8]], k: usize, scaled: u64) -> HashSet<u64> {
    let mh = microbio::max_hash_for_scaled(scaled);
    kmer_text_set(seqs, k).into_iter().map(|w| microbio::murmur64(&w, 42)).filter(|&h| h <= mh).collect()
}

fn oracle_jaccard(a: &HashSet<u64>, b: &HashSet<u64>) -> f64 {
    let i = a.intersection(b).count() as f64;
    let u = a.union(b).count() as f64;
    if u == 0.0 {
        0.0
    } else {
        i / u
    }
}

fn mutate_rate(seq: &[u8], rate_per_million: u64, s: &mut u64) -> Vec<u8> {
    seq.iter().map(|&b| if lcg(s) % 1_000_000 < rate_per_million { other_base(b, s) } else { b }).collect()
}

fn trial_distances(t: &mut T, s: &mut u64) {
    // Five genomes: a, b (1% from a), c (5% from a), d (a with N runs and
    // lowercase), e (unrelated).
    let a = random_dna(20_000, s);
    let b = mutate_rate(&a, 10_000, s);
    let c = mutate_rate(&a, 50_000, s);
    let mut d = a.clone();
    for i in 1000..1100 {
        d[i] = b'N';
    }
    for x in d[5000..9000].iter_mut() {
        *x = x.to_ascii_lowercase();
    }
    let e = random_dna(15_000, s);
    let genomes: Vec<(String, Vec<u8>)> =
        vec![("ga".into(), a.clone()), ("gb".into(), b.clone()), ("gc".into(), c.clone()), ("gd".into(), d.clone()), ("ge".into(), e.clone())];
    let multi = t.write("five.fasta", &fasta(&genomes));
    for &(k, scaled) in &[(21u8, 1u64), (31, 1), (21, 10), (31, 100), (15, 3)] {
        let opt = DistanceOptions { k, scaled, records_as_genomes: true, ..Default::default() };
        let r = distances_of(&[GenomeFiles { files: vec![multi.clone()], name: String::new() }], &opt);
        let Some(r) = r else {
            t.check(false, format!("distances k={k} scaled={scaled} failed"));
            continue;
        };
        t.check(r.scaled == scaled && r.k == k, format!("k {k} scaled {scaled}: the result says k {} scaled {}", r.k, r.scaled));
        let names: Vec<&str> = r.genomes.iter().map(|g| g.name.as_str()).collect();
        t.check(names == ["ga", "gb", "gc", "gd", "ge"], format!("records as genomes: names {names:?}"));
        let sets: Vec<HashSet<u64>> = genomes.iter().map(|(_, g)| oracle_hash_set(&[g], k as usize, scaled)).collect();
        let mut all_ok = true;
        for i in 0..5 {
            if r.genomes[i].hashes != sets[i].len() as u64 {
                all_ok = false;
                println!("    k {k} scaled {scaled}: genome {i} has {} hashes, the oracle {}", r.genomes[i].hashes, sets[i].len());
            }
            for j in 0..5 {
                let jj = oracle_jaccard(&sets[i], &sets[j]);
                let want_d = if jj == 0.0 { 1.0 } else if jj == 1.0 { 0.0 } else { (-(2.0 * jj / (1.0 + jj)).ln() / k as f64).min(1.0) };
                let want_ani = if jj == 0.0 { 0.0 } else { (2.0 * jj / (1.0 + jj)).powf(1.0 / k as f64) };
                if !(near(r.jaccard[i][j], jj, 1e-12) && near(r.distance[i][j], want_d, 1e-12) && near(r.ani[i][j], want_ani, 1e-12)) {
                    all_ok = false;
                    println!("    k {k} scaled {scaled} ({i},{j}): J {} vs {jj}, D {} vs {want_d}", r.jaccard[i][j], r.distance[i][j]);
                }
                if r.shared[i][j] != sets[i].intersection(&sets[j]).count() as u64 {
                    all_ok = false;
                }
            }
        }
        t.check(all_ok, format!("k {k} scaled {scaled}: every Jaccard, Mash distance, ANI and shared count equals the brute-force sets"));
        // Sanity of the numbers themselves.
        if scaled == 1 {
            let jab = r.jaccard[0][1];
            t.check(r.distance[0][1] > 0.005 && r.distance[0][1] < 0.02, format!("1% divergence gives D = {} (k {k})", r.distance[0][1]));
            t.check(r.ani[0][1] > 0.98 && r.ani[0][1] < 0.995, format!("1% divergence gives ANI = {}", r.ani[0][1]));
            t.check(jab > r.jaccard[0][2], "the 1% genome is closer than the 5% genome");
            t.check(r.distance[0][4] == 1.0 && r.ani[0][4] == 0.0, format!("unrelated genomes: D = {} ANI = {}", r.distance[0][4], r.ani[0][4]));
        }
    }
    // The formulas.
    t.check(mash_distance(1.0, 21) == 0.0 && mash_distance(0.0, 21) == 1.0, "Mash distance at J = 1 and J = 0");
    t.check(near(mash_distance(0.5, 21), -(2.0f64 / 3.0).ln() / 21.0, 1e-15), "Mash distance at J = 0.5");
    t.check(near(ani_from_jaccard(0.5, 31), (2.0f64 / 3.0).powf(1.0 / 31.0), 1e-15), "ANI at J = 0.5");

    // Files as genomes, automatic scaled: 20 kb -> 1; a 150 kb genome -> 10.
    let big = random_dna(150_000, s);
    let big2 = mutate_rate(&big, 5_000, s);
    let fa = t.write("a.fna", &fasta(&[("a1".into(), a[..10_000].to_vec()), ("a2".into(), a[10_000..].to_vec())]));
    let fb = t.write("b.fna", &fasta(&[("b".into(), b.clone())]));
    let fbig = t.write("big.fasta", &fasta(&[("big".into(), big.clone())]));
    let fbig2 = t.write("big2.fasta", &fasta(&[("big2".into(), big2.clone())]));
    let opt = DistanceOptions::default();
    if let Some(r) = distances_of(&[gf(&fa), gf(&fb)], &opt) {
        t.check(r.scaled == 1, format!("20 kb genomes: automatic scaled {} (want 1)", r.scaled));
        t.check(r.genomes[0].name == "a" && r.genomes[0].sequences == 2, format!("a file of 2 records is one genome named after the file: {:?}", r.genomes[0]));
        // Two records of one genome: the k-mers across the junction are not k-mers.
        let sa = oracle_hash_set(&[&a[..10_000], &a[10_000..]], 31, 1);
        let sb = oracle_hash_set(&[&b], 31, 1);
        t.check(near(r.jaccard[0][1], oracle_jaccard(&sa, &sb), 1e-12), "several records of one genome: the Jaccard of their k-mers together");
    }
    if let Some(r) = distances_of(&[gf(&fbig), gf(&fbig2), gf(&fa)], &opt) {
        t.check(r.scaled == 10, format!("150 kb and 20 kb genomes: automatic scaled {} (want 10)", r.scaled));
        let s1 = oracle_hash_set(&[&big], 31, 10);
        let s2 = oracle_hash_set(&[&big2], 31, 10);
        t.check(near(r.jaccard[0][1], oracle_jaccard(&s1, &s2), 1e-12), "automatic scaled: Jaccard at the common scaled equals the oracle");
        t.check(r.warnings.iter().any(|w| w.contains("only")) || r.genomes[2].hashes >= 100, "a small genome at a coarse scaled is warned about");
    }
    // Sketch databases: the same numbers.
    let db = t.p("five.mmgsk");
    let refs: Vec<microbio::GenomeInput> = [&fa, &fb, &fbig, &fbig2]
        .iter()
        .map(|f| microbio::GenomeInput { files: vec![f.to_string()], name: surto_label(f), host: false })
        .collect();
    microbio::sketch_references(&refs, &db, &microbio::SketchOptions { k: 21, scaled: 10 }).unwrap();
    let from_fasta = distances_of(&[gf(&fa), gf(&fb), gf(&fbig), gf(&fbig2)], &DistanceOptions { k: 21, scaled: 10, ..Default::default() });
    match (sketch_distances(&[db.clone()], &DistanceOptions::default()), from_fasta) {
        (Ok(x), Some(y)) => {
            let same = (0..4).all(|i| (0..4).all(|j| near(x.distance[i][j], y.distance[i][j], 1e-12)));
            t.check(same && x.k == 21 && x.scaled == 10, "a sketch database gives the distances of the FASTA files");
            t.check(x.tree.newick == y.tree.newick, "and the same tree");
        }
        (Err(e), _) => {
            t.check(false, format!("sketch_distances: {e}"));
        }
        _ => {}
    }
    match sketch_distances(&[db.clone()], &DistanceOptions { scaled: 100, ..Default::default() }) {
        Ok(x) => {
            let s1 = oracle_hash_set(&[&big], 21, 100);
            let s2 = oracle_hash_set(&[&big2], 21, 100);
            t.check(x.scaled == 100 && near(x.jaccard[2][3], oracle_jaccard(&s1, &s2), 1e-12), "a database made coarser (scaled 100)");
        }
        Err(e) => {
            t.check(false, format!("coarser database: {e}"));
        }
    }
    t.fails("a database made finer", sketch_distances(&[db.clone()], &DistanceOptions { scaled: 2, ..Default::default() }), &["scaled", "10"]);

    // Outputs: TSV and PHYLIP matrices read back, Newick read back.
    let mtsv = t.p("m.tsv");
    let mphy = t.p("m.phy");
    let nwk = t.p("tree.nwk");
    let opt = DistanceOptions {
        k: 21,
        scaled: 1,
        records_as_genomes: true,
        matrix_output: Some(mtsv.clone()),
        tree_output: Some(nwk.clone()),
        ..Default::default()
    };
    if let Some(r) = distances_of(&[gf(&multi)], &opt) {
        match read_distance_matrix(&mtsv) {
            Ok((names, d)) => {
                let ok = names.len() == 5 && (0..5).all(|i| (0..5).all(|j| (d[i][j] - r.distance[i][j]).abs() < 1e-8));
                t.check(ok, "the TSV matrix reads back to the same distances");
            }
            Err(e) => {
                t.check(false, format!("reading the TSV matrix: {e}"));
            }
        }
        let text = std::fs::read_to_string(&nwk).unwrap_or_default();
        match parse_newick(text.trim()) {
            Some(tree) => {
                let pd = r.tree.patristic_distances();
                let names: Vec<String> = r.genomes.iter().map(|g| g.name.clone()).collect();
                let mut ok = true;
                for i in 0..5 {
                    for j in 0..5 {
                        let x = tree.dist(&names[i], &names[j]);
                        if !x.is_some_and(|x| (x - pd[i][j]).abs() < 1e-6) {
                            ok = false;
                        }
                    }
                }
                t.check(ok, "the Newick file reads back to the tree's distances");
            }
            None => {
                t.check(false, format!("the Newick file does not parse: {text}"));
            }
        }
        let j = r.json();
        t.check(j.contains("\"newick\"") && j.contains("\"text\"") && j.contains("\"edges\""), "json has the tree, the edges and the text");
        t.check(r.resumo().contains("closest genome"), "the summary lists the closest genomes");
    }
    let opt = DistanceOptions { k: 21, scaled: 1, records_as_genomes: true, matrix_output: Some(mphy.clone()), ..Default::default() };
    if let Some(r) = distances_of(&[gf(&multi)], &opt) {
        let text = std::fs::read_to_string(&mphy).unwrap_or_default();
        t.check(text.starts_with("5\n"), "PHYLIP chosen by the .phy extension");
        match read_distance_matrix(&mphy) {
            Ok((names, d)) => {
                let ok = names[0] == "ga" && (0..5).all(|i| (0..5).all(|j| (d[i][j] - r.distance[i][j]).abs() < 1e-8));
                t.check(ok, "the PHYLIP matrix reads back");
            }
            Err(e) => {
                t.check(false, format!("reading PHYLIP: {e}"));
            }
        }
    }
    // Errors.
    t.fails("k = 0", genome_distances(&[gf(&fa)], &DistanceOptions { k: 0, ..Default::default() }), &["k = 0", "1 to 32"]);
    t.fails("no genome", genome_distances(&[], &DistanceOptions::default()), &["No genome", "FASTA"]);
    t.fails("missing file", genome_distances(&[gf(&t.p("nope.fa"))], &DistanceOptions::default()), &["nope.fa", "does not exist"]);
    let fq = t.write("reads.fq", "@r1\nACGTACGTAC\n+\nIIIIIIIIII\n");
    t.fails("FASTQ", genome_distances(&[gf(&fq)], &DistanceOptions::default()), &["reads.fq", "FASTQ", "Assemble"]);
    t.fails(
        "output = input",
        genome_distances(&[gf(&fa), gf(&fb)], &DistanceOptions { matrix_output: Some(fa.clone()), ..Default::default() }),
        &["also one of the input files"],
    );
    let empty = t.write("empty.fa", "");
    t.fails("empty FASTA", genome_distances(&[gf(&empty)], &DistanceOptions::default()), &["empty.fa", "no sequence"]);
    // Duplicate names.
    if let Some(r) = distances_of(&[GenomeFiles { files: vec![fa.clone()], name: "x".into() }, GenomeFiles { files: vec![fb.clone()], name: "x".into() }], &DistanceOptions::default()) {
        t.check(r.genomes[0].name == "x" && r.genomes[1].name == "x_2", format!("duplicate names made unique: {} {}", r.genomes[0].name, r.genomes[1].name));
    }
}

/// The name genome_distances gives a file: no folder, no extension.
fn surto_label(p: &str) -> String {
    let b = std::path::Path::new(p).file_name().unwrap().to_string_lossy().to_string();
    b.split('.').next().unwrap().to_string()
}

fn gf(p: &str) -> GenomeFiles {
    GenomeFiles { files: vec![p.to_string()], name: String::new() }
}

fn distances_of(g: &[GenomeFiles], opt: &DistanceOptions) -> Option<DistanceResult> {
    match genome_distances(g, opt) {
        Ok(r) => Some(r),
        Err(e) => {
            println!("    genome_distances failed: {e}");
            None
        }
    }
}

// =====================================================================
//  TREES
// =====================================================================

/// A Newick tree parsed back: leaves with their path to the root.
struct ParsedTree {
    parent: Vec<Option<usize>>,
    len: Vec<f64>,
    leaf: HashMap<String, usize>,
}

impl ParsedTree {
    fn depth(&self, mut u: usize) -> Vec<(usize, f64)> {
        let mut v = vec![(u, 0.0)];
        let mut acc = 0.0;
        while let Some(p) = self.parent[u] {
            acc += self.len[u];
            v.push((p, acc));
            u = p;
        }
        v
    }

    fn dist(&self, a: &str, b: &str) -> Option<f64> {
        let (&x, &y) = (self.leaf.get(a)?, self.leaf.get(b)?);
        let pa = self.depth(x);
        let pb: HashMap<usize, f64> = self.depth(y).into_iter().collect();
        pa.iter().find_map(|(u, da)| pb.get(u).map(|db| da + db))
    }
}

fn parse_newick(s: &str) -> Option<ParsedTree> {
    let b = s.as_bytes();
    let mut t = ParsedTree { parent: Vec::new(), len: Vec::new(), leaf: HashMap::new() };
    let mut stack: Vec<usize> = Vec::new();
    let mut i = 0;
    let mut last: Option<usize> = None;
    let new_node = |t: &mut ParsedTree, parent: Option<usize>| {
        t.parent.push(parent);
        t.len.push(0.0);
        t.parent.len() - 1
    };
    while i < b.len() {
        match b[i] {
            b'(' => {
                let u = new_node(&mut t, stack.last().copied());
                stack.push(u);
                i += 1;
                last = None;
            }
            b',' => {
                i += 1;
                last = None;
            }
            b')' => {
                last = stack.pop();
                i += 1;
            }
            b':' => {
                let st = i + 1;
                i += 1;
                while i < b.len() && !matches!(b[i], b',' | b')' | b';') {
                    i += 1;
                }
                let v: f64 = std::str::from_utf8(&b[st..i]).ok()?.parse().ok()?;
                t.len[last?] = v;
            }
            b';' => break,
            _ => {
                // a name
                let mut name = String::new();
                if b[i] == b'\'' {
                    i += 1;
                    loop {
                        if b[i] == b'\'' {
                            if i + 1 < b.len() && b[i + 1] == b'\'' {
                                name.push('\'');
                                i += 2;
                                continue;
                            }
                            i += 1;
                            break;
                        }
                        name.push(char::from(b[i]));
                        i += 1;
                    }
                } else {
                    while i < b.len() && !matches!(b[i], b':' | b',' | b')' | b';') {
                        name.push(char::from(b[i]));
                        i += 1;
                    }
                }
                if last.is_none() || b[i - 1] != b')' {
                    let u = new_node(&mut t, stack.last().copied());
                    t.leaf.insert(name, u);
                    last = Some(u);
                }
            }
        }
    }
    Some(t)
}

/// Splits of a rooted tree (taxon sets not holding taxon 0) with their
/// lengths; the two branches of a bifurcating root are one split.
fn tree_splits(tr: &Tree, n: usize) -> BTreeMap<Vec<usize>, f64> {
    let mut below: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); tr.nodes.len()];
    // post-order
    let mut order = vec![tr.root];
    let mut i = 0;
    while i < order.len() {
        let u = order[i];
        order.extend(tr.nodes[u].children.iter().copied());
        i += 1;
    }
    for &u in order.iter().rev() {
        if let Some(tx) = tr.nodes[u].taxon {
            below[u].insert(tx);
        }
        let ch = tr.nodes[u].children.clone();
        for c in ch {
            let x = below[c].clone();
            below[u].extend(x);
        }
    }
    let mut out: BTreeMap<Vec<usize>, f64> = BTreeMap::new();
    for (u, node) in tr.nodes.iter().enumerate() {
        if u == tr.root {
            continue;
        }
        let set = &below[u];
        let key: Vec<usize> = if set.contains(&0) { (0..n).filter(|x| !set.contains(x)).collect() } else { set.iter().copied().collect() };
        if key.is_empty() {
            continue;
        }
        *out.entry(key).or_insert(0.0) += node.branch_length;
    }
    out
}

/// A random unrooted binary tree: adjacency (leaves 0..n), and its path
/// distances.
fn random_tree(n: usize, s: &mut u64) -> (Vec<Vec<(usize, f64)>>, Vec<Vec<f64>>) {
    let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let len = |s: &mut u64| 0.05 + (lcg(s) % 1000) as f64 / 1000.0;
    // star of 3
    let c = adj.len();
    adj.push(Vec::new());
    for leaf in 0..3 {
        let l = len(s);
        adj[c].push((leaf, l));
        adj[leaf].push((c, l));
    }
    for leaf in 3..n {
        // split a random edge
        let mut edges = Vec::new();
        for u in 0..adj.len() {
            for &(v, l) in &adj[u] {
                if u < v {
                    edges.push((u, v, l));
                }
            }
        }
        let (u, v, l) = edges[rnd(s, edges.len())];
        adj[u].retain(|x| x.0 != v);
        adj[v].retain(|x| x.0 != u);
        let m = adj.len();
        adj.push(Vec::new());
        let a = l * 0.5;
        adj[u].push((m, a));
        adj[m].push((u, a));
        adj[v].push((m, l - a));
        adj[m].push((v, l - a));
        let ll = len(s);
        adj[m].push((leaf, ll));
        adj[leaf].push((m, ll));
    }
    let mut d = vec![vec![0.0; n]; n];
    for a in 0..n {
        let mut dist = vec![f64::NAN; adj.len()];
        dist[a] = 0.0;
        let mut st = vec![a];
        while let Some(u) = st.pop() {
            for &(v, l) in &adj[u] {
                if dist[v].is_nan() {
                    dist[v] = dist[u] + l;
                    st.push(v);
                }
            }
        }
        for b in 0..n {
            d[a][b] = dist[b];
        }
    }
    (adj, d)
}

fn adjacency_splits(adj: &[Vec<(usize, f64)>], n: usize) -> BTreeMap<Vec<usize>, f64> {
    let mut out = BTreeMap::new();
    for u in 0..adj.len() {
        for &(v, l) in &adj[u] {
            if u >= v {
                continue;
            }
            // leaves on v's side
            let mut seen = vec![false; adj.len()];
            seen[u] = true;
            seen[v] = true;
            let mut st = vec![v];
            let mut set = BTreeSet::new();
            while let Some(x) = st.pop() {
                if x < n {
                    set.insert(x);
                }
                for &(y, _) in &adj[x] {
                    if !seen[y] {
                        seen[y] = true;
                        st.push(y);
                    }
                }
            }
            let key: Vec<usize> = if set.contains(&0) { (0..n).filter(|x| !set.contains(x)).collect() } else { set.into_iter().collect() };
            out.insert(key, l);
        }
    }
    out
}

fn same_splits(a: &BTreeMap<Vec<usize>, f64>, b: &BTreeMap<Vec<usize>, f64>, tol: f64) -> bool {
    a.len() == b.len() && a.iter().all(|(k, v)| b.get(k).is_some_and(|w| (v - w).abs() <= tol))
}

fn names_of(v: &[&str]) -> Vec<String> {
    v.iter().map(|x| x.to_string()).collect()
}

fn layout_ok(tr: &Tree, n: usize) -> bool {
    let nodes = &tr.nodes;
    let leaves: Vec<usize> = (0..nodes.len()).filter(|&u| nodes[u].children.is_empty()).collect();
    if leaves.len() != n || tr.leaf_order.len() != n {
        return false;
    }
    let mut ys: Vec<f64> = leaves.iter().map(|&u| nodes[u].y).collect();
    ys.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if ys.iter().enumerate().any(|(i, &y)| y != i as f64) {
        return false;
    }
    if nodes[tr.root].x != 0.0 || nodes[tr.root].parent.is_some() {
        return false;
    }
    for (u, nd) in nodes.iter().enumerate() {
        if let Some(p) = nd.parent {
            if (nd.x - (nodes[p].x + nd.branch_length)).abs() > 1e-9 || nd.branch_length < 0.0 {
                return false;
            }
            if !nodes[p].children.contains(&u) {
                return false;
            }
        } else if u != tr.root {
            return false;
        }
        if let (Some(&f), Some(&l)) = (nd.children.first(), nd.children.last()) {
            if (nd.y - (nodes[f].y + nodes[l].y) / 2.0).abs() > 1e-12 {
                return false;
            }
        }
    }
    // The leaf order is the y order.
    tr.leaf_order.iter().enumerate().all(|(i, &tx)| nodes.iter().any(|nd| nd.taxon == Some(tx) && nd.y == i as f64)) && tr.edges().len() == nodes.len() - 1
}

fn trial_trees(t: &mut T, s: &mut u64) {
    // Wikipedia's neighbor-joining example (an additive matrix).
    let names = names_of(&["a", "b", "c", "d", "e"]);
    let d = vec![
        vec![0.0, 5.0, 9.0, 9.0, 8.0],
        vec![5.0, 0.0, 10.0, 10.0, 9.0],
        vec![9.0, 10.0, 0.0, 8.0, 7.0],
        vec![9.0, 10.0, 8.0, 0.0, 3.0],
        vec![8.0, 9.0, 7.0, 3.0, 0.0],
    ];
    match build_tree(&names, &d, TreeMethod::NeighborJoining, None) {
        Ok(tr) => {
            let sp = tree_splits(&tr, 5);
            let mut want = BTreeMap::new();
            want.insert(vec![1], 3.0); // b
            want.insert(vec![2], 4.0); // c
            want.insert(vec![3], 2.0); // d
            want.insert(vec![4], 1.0); // e
            want.insert(vec![1, 2, 3, 4], 2.0); // a
            want.insert(vec![2, 3, 4], 3.0); // {a,b} | {c,d,e}
            want.insert(vec![3, 4], 2.0); // {d,e}
            t.check(same_splits(&sp, &want, 1e-9), format!("Wikipedia NJ example: splits {sp:?}"));
            let pd = tr.patristic_distances();
            t.check((0..5).all(|i| (0..5).all(|j| (pd[i][j] - d[i][j]).abs() < 1e-9)), "Wikipedia NJ example: the tree gives the matrix back");
            // Midpoint: the longest path is b-c (10): the root at 5 from both.
            let x = |nm: &str| tr.nodes.iter().find(|n| n.name == nm).map(|n| n.x).unwrap_or(-1.0);
            t.check((x("b") - 5.0).abs() < 1e-9 && (x("c") - 5.0).abs() < 1e-9, format!("midpoint root: b at {} and c at {} (want 5)", x("b"), x("c")));
            t.check(layout_ok(&tr, 5), "Wikipedia NJ example: the layout is consistent");
            t.check(tr.newick.ends_with(';') && tr.newick.matches('(').count() == 4, format!("Newick {}", tr.newick));
        }
        Err(e) => {
            t.check(false, format!("NJ example: {e}"));
        }
    }
    // A 4-taxon NON-additive matrix, worked by hand (see the comments):
    //   R = A 19, B 22, C 22, D 21; Q(AB) = Q(CD) = -33 -> AB (first);
    //   A = 2 + (19 - 22)/4 = 1.25, B = 2.75; d(u,C) = 6.5, d(u,D) = 6;
    //   then C, D, u all tie at Q = -17.5 -> CD: C = 2.75, D = 2.25;
    //   u-v = (6.5 + 6 - 5)/2 = 3.75.
    let names4 = names_of(&["A", "B", "C", "D"]);
    let d4 = vec![vec![0.0, 4.0, 8.0, 7.0], vec![4.0, 0.0, 9.0, 9.0], vec![8.0, 9.0, 0.0, 5.0], vec![7.0, 9.0, 5.0, 0.0]];
    match build_tree(&names4, &d4, TreeMethod::NeighborJoining, None) {
        Ok(tr) => {
            let sp = tree_splits(&tr, 4);
            let mut want = BTreeMap::new();
            want.insert(vec![1, 2, 3], 1.25);
            want.insert(vec![1], 2.75);
            want.insert(vec![2], 2.75);
            want.insert(vec![3], 2.25);
            want.insert(vec![2, 3], 3.75);
            t.check(same_splits(&sp, &want, 1e-9), format!("hand-worked 4-taxon NJ: splits {sp:?}"));
            // longest path B-C = 2.75 + 3.75 + 2.75 = 9.25: root at 4.625
            let x = |nm: &str| tr.nodes.iter().find(|n| n.name == nm).map(|n| n.x).unwrap_or(-1.0);
            t.check((x("B") - 4.625).abs() < 1e-9 && (x("C") - 4.625).abs() < 1e-9, "hand-worked 4-taxon: midpoint at 4.625");
        }
        Err(e) => {
            t.check(false, format!("4-taxon NJ: {e}"));
        }
    }
    // Wikipedia's UPGMA example (5S rRNA of five bacteria).
    let d5 = vec![
        vec![0.0, 17.0, 21.0, 31.0, 23.0],
        vec![17.0, 0.0, 30.0, 34.0, 21.0],
        vec![21.0, 30.0, 0.0, 28.0, 39.0],
        vec![31.0, 34.0, 28.0, 0.0, 43.0],
        vec![23.0, 21.0, 39.0, 43.0, 0.0],
    ];
    match build_tree(&names, &d5, TreeMethod::Upgma, None) {
        Ok(tr) => {
            let pd = tr.patristic_distances();
            t.check((pd[0][1] - 17.0).abs() < 1e-9 && (pd[0][4] - 22.0).abs() < 1e-9 && (pd[2][3] - 28.0).abs() < 1e-9 && (pd[0][2] - 33.0).abs() < 1e-9, format!("UPGMA example: tree distances {pd:?}"));
            let ultra = tr.nodes.iter().filter(|n| n.children.is_empty()).all(|n| (n.x - 16.5).abs() < 1e-9);
            t.check(ultra, "UPGMA example: every leaf at height 16.5");
            let sp = tree_splits(&tr, 5);
            let mut want = BTreeMap::new();
            want.insert(vec![1, 2, 3, 4], 8.5);
            want.insert(vec![1], 8.5);
            want.insert(vec![2], 14.0);
            want.insert(vec![3], 14.0);
            want.insert(vec![4], 11.0);
            want.insert(vec![2, 3, 4], 2.5);
            want.insert(vec![2, 3], 2.5 + 5.5);
            t.check(same_splits(&sp, &want, 1e-9), format!("UPGMA example: splits {sp:?}"));
            t.check(layout_ok(&tr, 5), "UPGMA layout consistent");
        }
        Err(e) => {
            t.check(false, format!("UPGMA example: {e}"));
        }
    }
    // Random additive trees: neighbor joining gives them back exactly.
    let mut ok_topo = 0;
    let mut ok_dist = 0;
    let mut ok_layout = 0;
    let mut ok_mid = 0;
    let mut ok_newick = 0;
    let trials = 300;
    for k in 0..trials {
        let n = 4 + k % 10;
        let (adj, d) = random_tree(n, s);
        let names: Vec<String> = (0..n).map(|i| format!("t{i}")).collect();
        let Ok(tr) = build_tree(&names, &d, TreeMethod::NeighborJoining, None) else { continue };
        if same_splits(&tree_splits(&tr, n), &adjacency_splits(&adj, n), 1e-9) {
            ok_topo += 1;
        }
        let pd = tr.patristic_distances();
        if (0..n).all(|i| (0..n).all(|j| (pd[i][j] - d[i][j]).abs() < 1e-9)) {
            ok_dist += 1;
        }
        if layout_ok(&tr, n) {
            ok_layout += 1;
        }
        // Midpoint: the farthest pair sits at half the diameter.
        let diam = d.iter().flatten().fold(0.0f64, |a, &b| a.max(b));
        let maxx = tr.nodes.iter().filter(|n| n.children.is_empty()).map(|n| n.x).fold(0.0f64, f64::max);
        if (maxx - diam / 2.0).abs() < 1e-9 {
            ok_mid += 1;
        }
        if let Some(pt) = parse_newick(&tr.newick) {
            if (0..n).all(|i| (0..n).all(|j| pt.dist(&names[i], &names[j]).is_some_and(|x| (x - d[i][j]).abs() < 1e-6))) {
                ok_newick += 1;
            }
        }
    }
    t.check(ok_topo == trials, format!("random additive trees: the topology recovered in {ok_topo} of {trials}"));
    t.check(ok_dist == trials, format!("random additive trees: the distances recovered in {ok_dist} of {trials}"));
    t.check(ok_layout == trials, format!("random trees: consistent layout in {ok_layout} of {trials}"));
    t.check(ok_mid == trials, format!("random trees: root at half the diameter in {ok_mid} of {trials}"));
    t.check(ok_newick == trials, format!("random trees: Newick reads back in {ok_newick} of {trials}"));
    // Random additive trees with UPGMA when the tree is ultrametric: build
    // from a random ultrametric (coalescent-like) tree.
    let mut ok_u = 0;
    for _ in 0..100 {
        let n = 3 + rnd(s, 8);
        // random merge order with increasing heights
        let mut clusters: Vec<(Vec<usize>, f64)> = (0..n).map(|i| (vec![i], 0.0)).collect();
        let mut d = vec![vec![0.0; n]; n];
        let mut h = 0.0;
        while clusters.len() > 1 {
            h += 0.1 + (lcg(s) % 100) as f64 / 100.0;
            let a = rnd(s, clusters.len());
            let ca = clusters.remove(a);
            let b = rnd(s, clusters.len());
            let cb = clusters.remove(b);
            for &x in &ca.0 {
                for &y in &cb.0 {
                    d[x][y] = 2.0 * h;
                    d[y][x] = 2.0 * h;
                }
            }
            let mut m = ca.0.clone();
            m.extend(cb.0);
            clusters.push((m, h));
        }
        let names: Vec<String> = (0..n).map(|i| format!("u{i}")).collect();
        if let Ok(tr) = build_tree(&names, &d, TreeMethod::Upgma, None) {
            let pd = tr.patristic_distances();
            if (0..n).all(|i| (0..n).all(|j| (pd[i][j] - d[i][j]).abs() < 1e-9)) && layout_ok(&tr, n) {
                ok_u += 1;
            }
        }
    }
    t.check(ok_u == 100, format!("random ultrametric trees: UPGMA gives them back in {ok_u} of 100"));
    // Random non-additive matrices: a valid tree, no negative branch.
    let mut ok_na = 0;
    for _ in 0..200 {
        let n = 3 + rnd(s, 10);
        let mut d = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in i + 1..n {
                let v = (lcg(s) % 1000) as f64 / 100.0;
                d[i][j] = v;
                d[j][i] = v;
            }
        }
        let names: Vec<String> = (0..n).map(|i| format!("r{i}")).collect();
        for m in [TreeMethod::NeighborJoining, TreeMethod::Upgma] {
            if let Ok(tr) = build_tree(&names, &d, m, None) {
                if layout_ok(&tr, n) && tr.nodes.iter().all(|x| x.branch_length >= 0.0) && parse_newick(&tr.newick).is_some_and(|p| p.leaf.len() == n) {
                    ok_na += 1;
                }
            }
        }
    }
    t.check(ok_na == 400, format!("random non-additive matrices: a valid tree with no negative branch in {ok_na} of 400"));
    // Small cases.
    match build_tree(&names_of(&["solo"]), &[vec![0.0]], TreeMethod::NeighborJoining, None) {
        Ok(tr) => {
            t.check(tr.newick == "solo;" && tr.nodes.len() == 1 && tr.leaf_order == vec![0], format!("one genome: {}", tr.newick));
        }
        Err(e) => {
            t.check(false, format!("one genome: {e}"));
        }
    }
    match build_tree(&names_of(&["x", "y"]), &[vec![0.0, 0.4], vec![0.4, 0.0]], TreeMethod::NeighborJoining, None) {
        Ok(tr) => {
            t.check(tr.newick == "(x:0.2,y:0.2);", format!("two genomes: {}", tr.newick));
        }
        Err(e) => {
            t.check(false, format!("two genomes: {e}"));
        }
    }
    match build_tree(&names_of(&["p", "q", "r", "s"]), &vec![vec![0.0; 4]; 4], TreeMethod::NeighborJoining, None) {
        Ok(tr) => {
            t.check(layout_ok(&tr, 4) && tr.nodes.iter().all(|n| n.branch_length == 0.0), format!("identical genomes: a tree of zero branches {}", tr.newick));
        }
        Err(e) => {
            t.check(false, format!("identical genomes: {e}"));
        }
    }
    // Names that need quoting in Newick.
    match build_tree(&names_of(&["A b", "c:d", "e'f"]), &[vec![0.0, 1.0, 2.0], vec![1.0, 0.0, 2.0], vec![2.0, 2.0, 0.0]], TreeMethod::NeighborJoining, None) {
        Ok(tr) => {
            let p = parse_newick(&tr.newick);
            t.check(p.is_some_and(|p| p.leaf.contains_key("A b") && p.leaf.contains_key("c:d") && p.leaf.contains_key("e'f")), format!("quoted names: {}", tr.newick));
        }
        Err(e) => {
            t.check(false, format!("quoted names: {e}"));
        }
    }
    // Outgroup rooting: the root halves the outgroup's branch.
    match build_tree(&names, &d, TreeMethod::NeighborJoining, Some("c")) {
        Ok(tr) => {
            let x = |nm: &str| tr.nodes.iter().find(|n| n.name == nm).map(|n| n.x).unwrap_or(-1.0);
            let rc_: Vec<usize> = tr.nodes[tr.root].children.clone();
            t.check(rc_.len() == 2 && rc_.iter().any(|&c| tr.nodes[c].name == "c") && (x("c") - 2.0).abs() < 1e-9, format!("outgroup c: the root halves c's branch ({})", tr.newick));
            let pd = tr.patristic_distances();
            t.check((0..5).all(|i| (0..5).all(|j| (pd[i][j] - d[i][j]).abs() < 1e-9)) && layout_ok(&tr, 5), "outgroup rooting keeps the distances and the layout");
        }
        Err(e) => {
            t.check(false, format!("outgroup: {e}"));
        }
    }
    t.fails("an unknown outgroup", build_tree(&names, &d, TreeMethod::NeighborJoining, Some("zz")), &["\"zz\"", "not one of the genomes", "a, b"]);
    t.fails("UPGMA with an outgroup", build_tree(&names, &d5, TreeMethod::Upgma, Some("a")), &["UPGMA", "outgroup", "nj"]);
    let qcsv = t.write("quoted.csv", "\"x1 virus, complete\",\"x2 \"\"strain\"\", complete\"\n1,0.5\n0.5,1\n");
    t.check(
        read_distance_matrix(&qcsv).is_ok_and(|(n, m)| n == names_of(&["x1 virus, complete", "x2 \"strain\", complete"]) && m[0][1] == 0.5),
        "a CSV with quoted names holding commas and quotes",
    );
    // Matrix files.
    let csv = t.write("sm.csv", "a,b,c,d,e\n1,0.9,0.2,0.21,0.2\n0.9,1,0.2,0.2,0.2\n0.2,0.2,1,0.8,0.7\n0.21,0.2,0.8,1,0.6\n0.2,0.2,0.7,0.6,1\n");
    match tree_from_matrix(&csv, TreeMethod::NeighborJoining, true, None, None) {
        Ok(tr) => {
            let pd = tr.patristic_distances();
            t.check(tr.nodes.iter().filter(|n| n.children.is_empty()).count() == 5 && pd[0][1] < pd[0][2], "a sourmash-style CSV of similarities (no row names)");
        }
        Err(e) => {
            t.check(false, format!("similarity CSV: {e}"));
        }
    }
    let lower = t.write("lower.phy", "4\nA\nB 4\nC 8 9\nD 7 9 5\n");
    match tree_from_matrix(&lower, TreeMethod::NeighborJoining, false, None, Some(&t.p("lower.nwk"))) {
        Ok(tr) => {
            t.check(same_splits(&tree_splits(&tr, 4), &{
                let mut w = BTreeMap::new();
                w.insert(vec![1, 2, 3], 1.25);
                w.insert(vec![1], 2.75);
                w.insert(vec![2], 2.75);
                w.insert(vec![3], 2.25);
                w.insert(vec![2, 3], 3.75);
                w
            }, 1e-9), "a lower-triangle PHYLIP gives the hand-worked tree");
            t.check(std::fs::read_to_string(t.p("lower.nwk")).is_ok_and(|x| x.trim() == tr.newick), "tree_from_matrix writes the Newick");
        }
        Err(e) => {
            t.check(false, format!("lower PHYLIP: {e}"));
        }
    }
    let tsv = t.write("labelled.tsv", "genome\tA\tB\tC\nA\t0\t1\t2\nB\t1\t0\t2\nC\t2\t2\t0\n");
    t.check(read_distance_matrix(&tsv).is_ok_and(|(n, d)| n == names_of(&["A", "B", "C"]) && d[0][2] == 2.0), "a TSV with names on both sides");
    let asym = t.write("asym.tsv", "x\tA\tB\nA\t0\t1\nB\t2\t0\n");
    t.fails("asymmetric matrix", tree_from_matrix(&asym, TreeMethod::NeighborJoining, false, None, None), &["not symmetric", "A", "B"]);
    let diag = t.write("diag.tsv", "x\tA\tB\nA\t1\t0.5\nB\t0.5\t1\n");
    t.fails("similarity given as distance", tree_from_matrix(&diag, TreeMethod::NeighborJoining, false, None, None), &["to itself", "1 - similarity"]);
    let neg = t.write("neg.tsv", "x\tA\tB\nA\t0\t-1\nB\t-1\t0\n");
    t.fails("negative distance", tree_from_matrix(&neg, TreeMethod::Upgma, false, None, None), &[">= 0"]);
    let bad = t.write("bad.tsv", "x\tA\tB\nA\t0\tfoo\nB\t1\t0\n");
    t.fails("not a number", read_distance_matrix(&bad), &["bad.tsv", "line 2", "foo"]);
    let ragged = t.write("ragged.tsv", "x\tA\tB\tC\nA\t0\t1\nB\t1\t0\n");
    t.fails("not square", read_distance_matrix(&ragged), &["ragged.tsv", "names"]);
    t.fails("unknown method", TreeMethod::from_text("ml"), &["nj", "upgma"]);
    t.check(TreeMethod::from_text("Neighbor-Joining") == Ok(TreeMethod::NeighborJoining), "method names are forgiving");
}

// =====================================================================
//  2. VIRAL MUTATIONS
// =====================================================================

/// The synthetic reference: genes, designed codons, repeats.
struct Synth {
    seq: Vec<u8>,
    gff: String,
    /// Genes as (name, 0-based positions in coding order, minus strand),
    /// with ORF1ab split (ORF1a, ORF1b), pp1a not reported.
    genes: Vec<(String, Vec<usize>, bool)>,
    /// The same with ORF1ab whole and the pp1a CDS as ORF1a.
    genes_whole: Vec<(String, Vec<usize>, bool)>,
}

const L: usize = 30_000;

fn gene_a_codon(c: usize) -> usize {
    100 + 3 * (c - 1) // 0-based start of codon c of geneA
}

fn build_synth(s: &mut u64) -> Synth {
    let mut seq = random_dna(L, s);
    let put_codons = |seq: &mut Vec<u8>, start: usize, n: usize, s: &mut u64| {
        for c in 0..n {
            let cod = if c == 0 { *b"ATG" } else if c == n - 1 { *b"TAA" } else { sense_codon(s) };
            seq[start + 3 * c..start + 3 * c + 3].copy_from_slice(&cod);
        }
    };
    // geneA: 101..3100 (+), 1000 codons.
    put_codons(&mut seq, 100, 1000, s);
    let mut design = |c: usize, cod: &[u8]| {
        let p = gene_a_codon(c);
        seq[p..p + cod.len()].copy_from_slice(cod);
    };
    design(16, b"GTTGCA");
    design(67, b"GCTATACATGTCTCT"); // A67 I68 H69 V70 S71
    design(155, b"TCTGAGTTCAGAGTT"); // S155 E156 F157 R158 V159
    design(214, b"TTACCC"); // L214 P215
    design(299, b"TTTGGCAAATCT"); // F299 G300 K301 S302
    design(500, b"CAG");
    design(600, b"AGGGGA"); // R600 G601
    design(700, b"AAA");
    design(899, b"TTTCAGCAGCAGCAGTTT"); // F899 Q900..Q903 F904
    // Intergenic repeats.
    seq[3299..3307].copy_from_slice(b"GAAAAAAA");
    seq[3307] = b'C';
    seq[3499..3511].copy_from_slice(b"GACACACACACT");
    // ORF1ab: part 1 4001..7000, part 2 7000..9999 (-1 frameshift), pp1a
    // 4001..7012 (stop at 7010..7012).
    put_codons(&mut seq, 4000, 1000, s);
    seq[6997..7000].copy_from_slice(&[b'C', b'T', b'A']); // part-1 codon 1000 = CTA (L), base 7000 = A
    loop {
        // bases 7001..7014 (0-based 7000..7013) satisfy both frames
        for i in 7000..7014 {
            seq[i] = b"ACGT"[rnd(s, 4)];
        }
        let f2 = |a: usize| oracle_codon(&seq[a..a + 3]) != b'*';
        let pp = (0..3).all(|k| oracle_codon(&seq[7000 + 3 * k..7003 + 3 * k]) != b'*') && matches!(&seq[7009..7012], b"TAA" | b"TAG" | b"TGA");
        if pp && f2(6999) && f2(7002) && f2(7005) && f2(7008) && f2(7011) {
            break;
        }
    }
    for c in 5..999 {
        let cod = sense_codon(s);
        let p = 6999 + 3 * c;
        seq[p..p + 3].copy_from_slice(&cod);
    }
    seq[9996..9999].copy_from_slice(b"TAG");
    // geneM: minus strand 12001..13500 (500 codons).
    let mut m = vec![b'N'; 1500];
    for c in 0..500 {
        let cod = if c == 0 { *b"ATG" } else if c == 499 { *b"TGA" } else { sense_codon(s) };
        m[3 * c..3 * c + 3].copy_from_slice(&cod);
    }
    seq[12000..13500].copy_from_slice(&rc(&m));
    // geneS: spliced, 15001..15600 + 16001..16900 (500 codons).
    let mut sp = vec![b'N'; 1500];
    for c in 0..500 {
        let cod = if c == 0 { *b"ATG" } else if c == 499 { *b"TAA" } else { sense_codon(s) };
        sp[3 * c..3 * c + 3].copy_from_slice(&cod);
    }
    seq[15000..15600].copy_from_slice(&sp[..600]);
    seq[16000..16900].copy_from_slice(&sp[600..]);
    // AAAA running into the start codon of geneS (positions 14998-15001), as
    // the AAAA before N in SARS-CoV-2.
    seq[14996] = b'C';
    seq[14997..15000].copy_from_slice(b"AAA");
    let gff = format!(
        "##gff-version 3\n##sequence-region synth 1 {L}\n\
synth\tTest\tgene\t101\t3100\t.\t+\t.\tID=gene-A;Name=geneA\n\
synth\tTest\tCDS\t101\t3100\t.\t+\t0\tID=cds-A;Parent=gene-A;gene=geneA;product=protein A;protein_id=PA.1\n\
synth\tTest\tCDS\t4001\t7000\t.\t+\t0\tID=cds-1ab;gene=ORF1ab;product=ORF1ab polyprotein;exception=ribosomal slippage;Note=pp1ab%3B -1 frameshift\n\
synth\tTest\tCDS\t7000\t9999\t.\t+\t0\tID=cds-1ab;gene=ORF1ab;product=ORF1ab polyprotein;exception=ribosomal slippage\n\
synth\tTest\tCDS\t4001\t7012\t.\t+\t0\tID=cds-1a;gene=ORF1ab;product=ORF1a polyprotein\n\
synth\tTest\tCDS\t12001\t13500\t.\t-\t0\tID=cds-M;gene=geneM;product=protein M\n\
synth\tTest\tCDS\t16001\t16900\t.\t+\t0\tID=cds-S;gene=geneS;product=protein S\n\
synth\tTest\tCDS\t15001\t15600\t.\t+\t0\tID=cds-S;gene=geneS;product=protein S\n"
    );
    let plus = |a: usize, b: usize| (a - 1..b).collect::<Vec<usize>>();
    let gene_s: Vec<usize> = plus(15001, 15600).into_iter().chain(plus(16001, 16900)).collect();
    let gene_m: Vec<usize> = plus(12001, 13500).into_iter().rev().collect();
    let genes = vec![
        ("geneA".to_string(), plus(101, 3100), false),
        ("ORF1a".to_string(), plus(4001, 7000), false),
        ("ORF1b".to_string(), plus(7000, 9999), false),
        ("geneM".to_string(), gene_m.clone(), true),
        ("geneS".to_string(), gene_s.clone(), false),
    ];
    let whole: Vec<usize> = plus(4001, 7000).into_iter().chain(plus(7000, 9999)).collect();
    let genes_whole = vec![
        ("geneA".to_string(), plus(101, 3100), false),
        ("ORF1ab".to_string(), whole, false),
        ("ORF1a".to_string(), plus(4001, 7012), false),
        ("geneM".to_string(), gene_m, true),
        ("geneS".to_string(), gene_s, false),
    ];
    Synth { seq, gff, genes, genes_whole }
}

/// Codon-by-codon amino-acid changes between the reference and a genome
/// that differs only by substitutions (same length).
fn oracle_aa_subs(refs: &[u8], q: &[u8], genes: &[(String, Vec<usize>, bool)]) -> Vec<String> {
    let mut out = Vec::new();
    for (name, pos, minus) in genes {
        let get = |seq: &[u8], p: usize| if *minus { comp(seq[p]) } else { seq[p] };
        for c in 0..pos.len() / 3 {
            let rc_: Vec<u8> = (0..3).map(|k| get(refs, pos[3 * c + k])).collect();
            let qc: Vec<u8> = (0..3).map(|k| get(q, pos[3 * c + k])).collect();
            let (ra, qa) = (oracle_codon(&rc_), oracle_codon(&qc));
            if qa != b'X' && qa != ra {
                out.push(format!("{name}:{}{}{}", char::from(ra), c + 1, char::from(qa)));
            }
        }
    }
    out
}

/// The cost of opening a gap at a reference position (0-based), as
/// Nextclade sets it: 6 outside genes, 7 at the first base (lowest
/// coordinate) of a codon, 8 inside a codon.
#[derive(Default)]
struct GapCost(HashMap<usize, i32>);

impl GapCost {
    fn get(&self, p: usize) -> i32 {
        self.0.get(&p).copied().unwrap_or(6)
    }
}

fn in_frame_set(genes: &[(String, Vec<usize>, bool)]) -> GapCost {
    let mut m = HashMap::new();
    for (_, pos, _) in genes {
        for &p in pos {
            m.entry(p).or_insert(8);
        }
    }
    for (_, pos, _) in genes {
        for c in pos.chunks_exact(3) {
            m.insert(*c.iter().min().unwrap(), 7);
        }
    }
    GapCost(m)
}

/// Where a deletion of `len` bases at `p0` (0-based) is reported: every
/// place giving the same sequence (brute force), the cheapest to open a gap
/// at, the leftmost of those.
fn oracle_del(refs: &[u8], p0: usize, len: usize, inf: &GapCost) -> usize {
    let lo = p0.saturating_sub(80);
    let hi = (p0 + len + 80).min(refs.len());
    let res = |p: usize| {
        let mut v = refs[lo..p].to_vec();
        v.extend_from_slice(&refs[p + len..hi]);
        v
    };
    let target = res(p0);
    let cands: Vec<usize> = (lo.max(p0.saturating_sub(60))..=(p0 + 60).min(hi - len)).filter(|&p| p >= lo && res(p) == target).collect();
    cands.iter().copied().min_by_key(|&p| (inf.get(p), p)).unwrap()
}

/// Where an insertion of `ins` after `a0` reference bases is reported, and
/// the bases as reported.
fn oracle_ins(refs: &[u8], a0: usize, ins: &[u8], inf: &GapCost) -> (usize, Vec<u8>) {
    let lo = a0.saturating_sub(80);
    let hi = (a0 + 80).min(refs.len());
    let mut q = refs[lo..a0].to_vec();
    q.extend_from_slice(ins);
    q.extend_from_slice(&refs[a0..hi]);
    let n = ins.len();
    let mut cands = Vec::new();
    for a in a0.saturating_sub(60).max(lo)..=(a0 + 60).min(hi) {
        let k = a - lo;
        if q[..k] == refs[lo..a] && q[k + n..] == refs[a..hi] {
            cands.push((a, q[k..k + n].to_vec()));
        }
    }
    cands.iter().min_by_key(|(a, _)| (inf.get(*a), *a)).cloned().unwrap()
}

/// The planted events of one genome, applied from the right so that the
/// positions stay those of the reference.
#[derive(Clone)]
enum Ev {
    Sub(usize, u8),
    Del(usize, usize),
    Ins(usize, Vec<u8>),
}

fn apply(refs: &[u8], evs: &[Ev]) -> Vec<u8> {
    let mut evs = evs.to_vec();
    evs.sort_by_key(|e| std::cmp::Reverse(match e {
        Ev::Sub(p, _) => *p,
        Ev::Del(p, _) => *p,
        Ev::Ins(a, _) => *a,
    }));
    let mut q = refs.to_vec();
    for e in evs {
        match e {
            Ev::Sub(p, b) => q[p] = b,
            Ev::Del(p, l) => {
                q.drain(p..p + l);
            }
            Ev::Ins(a, seq) => {
                q.splice(a..a, seq);
            }
        }
    }
    q
}

fn labels_subs(g: &GenomeMutations) -> Vec<String> {
    g.substitutions.iter().map(|x| x.label()).collect()
}
fn labels_dels(g: &GenomeMutations) -> Vec<String> {
    g.deletions.iter().map(|x| x.short()).collect()
}
fn labels_ins(g: &GenomeMutations) -> Vec<String> {
    g.insertions.iter().map(|x| x.short()).collect()
}
fn labels_aa(g: &GenomeMutations) -> Vec<String> {
    g.aa_substitutions.iter().chain(&g.aa_deletions).map(|x| x.label()).chain(g.aa_insertions.iter().map(|x| x.label())).collect()
}

fn del_label(p: usize, len: usize) -> String {
    if len == 1 {
        format!("{}", p + 1)
    } else {
        format!("{}-{}", p + 1, p + len)
    }
}

fn trial_mutations(t: &mut T, s: &mut u64) {
    let sy = build_synth(s);
    let r = &sy.seq;
    let refp = t.write("synth_ref.fasta", &fasta(&[("synth".into(), r.clone())]));
    let gff = t.write("synth.gff3", &sy.gff);
    let inf = in_frame_set(&sy.genes);
    // ---- genome 1: substitutions (random, plus designed codon changes)
    let mut subs: Vec<(usize, u8)> = Vec::new();
    let mut used = BTreeSet::new();
    let add_sub = |p: usize, b: u8, subs: &mut Vec<(usize, u8)>, used: &mut BTreeSet<usize>| {
        subs.push((p, b));
        used.insert(p);
    };
    let pa = gene_a_codon(500);
    add_sub(pa, b'T', &mut subs, &mut used); // CAG -> TAG: Q500*
    let p600 = gene_a_codon(600);
    add_sub(p600 + 1, b'A', &mut subs, &mut used); // AGG GGA -> AAA CGA: R600K, G601R
    add_sub(p600 + 2, b'A', &mut subs, &mut used);
    add_sub(p600 + 3, b'C', &mut subs, &mut used);
    add_sub(6000, other_base(r[6000], s), &mut subs, &mut used); // ORF1a
    add_sub(8000, other_base(r[8000], s), &mut subs, &mut used); // ORF1b
    add_sub(7004, other_base(r[7004], s), &mut subs, &mut used); // ORF1b codon 2, pp1a tail
    add_sub(12500, other_base(r[12500], s), &mut subs, &mut used); // geneM (minus)
    add_sub(16500, other_base(r[16500], s), &mut subs, &mut used); // geneS exon 2
    add_sub(15100, other_base(r[15100], s), &mut subs, &mut used); // geneS exon 1
    while subs.len() < 70 {
        let p = 200 + rnd(s, L - 400);
        if used.range(p.saturating_sub(40)..p + 40).next().is_none() {
            add_sub(p, other_base(r[p], s), &mut subs, &mut used);
        }
    }
    let g1 = apply(r, &subs.iter().map(|&(p, b)| Ev::Sub(p, b)).collect::<Vec<_>>());
    let want_subs1: Vec<String> = subs.iter().map(|&(p, b)| format!("{}{}{}", char::from(r[p]), p + 1, char::from(b))).collect();
    let want_aa1 = oracle_aa_subs(r, &g1, &sy.genes);
    // ---- genome 2: the indel cases in genes
    let mut ev2: Vec<Ev> = Vec::new();
    let c68 = gene_a_codon(68);
    ev2.push(Ev::Del(c68 + 1, 6)); // H69-, V70- (I68 kept: ATA -> ATC)
    let c156 = gene_a_codon(156);
    ev2.push(Ev::Del(c156 + 1, 6)); // E156-, F157-, R158G
    let c300 = gene_a_codon(300);
    ev2.push(Ev::Del(c300, 6)); // G300-, K301-
    let a214 = gene_a_codon(215); // after codon 214
    ev2.push(Ev::Ins(a214, b"GAGCCAGAA".to_vec()));
    let c16 = gene_a_codon(16);
    ev2.push(Ev::Ins(c16 + 1, b"TCATGCCGCTGT".to_vec())); // 16:MPLF
    let c902 = gene_a_codon(902);
    ev2.push(Ev::Del(c902, 3)); // a Q of the QQQQ repeat
    ev2.push(Ev::Sub(25_000, other_base(r[25_000], s)));
    let g2 = apply(r, &ev2);
    let mut want_dels2 = Vec::new();
    let mut want_ins2 = Vec::new();
    for e in &ev2 {
        match e {
            Ev::Del(p, l) => want_dels2.push(del_label(oracle_del(r, *p, *l, &inf), *l)),
            Ev::Ins(a, sq) => {
                let (aa, sq2) = oracle_ins(r, *a, sq, &inf);
                want_ins2.push(format!("{aa}:{}", st(&sq2)));
            }
            _ => {}
        }
    }
    let q900 = oracle_del(r, c902, 3, &inf);
    let want_aa2: Vec<String> = vec![
        "geneA:H69-".into(),
        "geneA:V70-".into(),
        "geneA:E156-".into(),
        "geneA:F157-".into(),
        "geneA:R158G".into(),
        "geneA:G300-".into(),
        "geneA:K301-".into(),
        "geneA:214:EPE".into(),
        "geneA:16:MPLF".into(),
        format!("geneA:Q{}-", (q900 - 100) / 3 + 1),
    ];
    // ---- genome 3: ends not covered (from 55 to L - 70), a few subs
    let mut g3 = g1[54..L - 70].to_vec();
    g3[100] = g3[100]; // unchanged
    // ---- genome 4: overhangs
    let mut g4 = random_dna(30, s);
    g4.extend_from_slice(r);
    g4.extend(random_dna(25, s));
    // ---- genome 5: N runs and ambiguity codes
    let mut g5 = r.clone();
    for x in g5[2000..2100].iter_mut() {
        *x = b'N';
    }
    g5[10_999] = b'N';
    for x in g5[19_999..20_036].iter_mut() {
        *x = b'N';
    }
    g5[24_999] = if r[24_999] == b'C' || r[24_999] == b'T' { b'Y' } else { b'R' };
    let k700 = gene_a_codon(700);
    g5[k700 + 2] = b'R'; // AAA -> AAR: still K
    // ---- genome 6: genome 1 reverse-complemented
    let g6 = rc(&g1);
    // ---- genome 7: intergenic indels in repeats, long ones
    let ev7 = vec![
        Ev::Del(3303, 2),                       // two A of AAAAAAA
        Ev::Ins(3505, b"AC".to_vec()),          // one more AC in ACACACACAC
        Ev::Del(22_000, 500),                   // a long deletion
        Ev::Ins(26_000, random_dna(300, s)),    // a long insertion
        Ev::Del(28_000, 1),
        Ev::Del(14_999, 1),                     // an A of AAAA(TG): placed before geneS, no frameshift
        Ev::Sub(28_500, other_base(r[28_500], s)),
    ];
    let g7 = apply(r, &ev7);
    let mut want_dels7 = Vec::new();
    let mut want_ins7 = Vec::new();
    for e in &ev7 {
        match e {
            Ev::Del(p, l) => want_dels7.push(del_label(oracle_del(r, *p, *l, &inf), *l)),
            Ev::Ins(a, sq) => {
                let (aa, sq2) = oracle_ins(r, *a, sq, &inf);
                want_ins7.push(format!("{aa}:{}", st(&sq2)));
            }
            _ => {}
        }
    }
    // ---- genome 8: frameshifts. geneM (minus) loses one base; geneS loses
    // one base at codon 250 and gains one in codon 255.
    let pm = 12_750; // inside geneM
    let ps_del = 16_000 + 49 * 3 + 1; // geneS codon 250 (exon-2 codon 50), 2nd base
    let ps_ins = 16_000 + 54 * 3 + 3; // after the 3rd base of geneS codon 255
    let ins_base = (b"ACGT".iter().copied()).find(|&b| b != r[ps_ins - 1] && b != r[ps_ins]).unwrap();
    let ev8 = vec![Ev::Del(pm, 1), Ev::Del(ps_del, 1), Ev::Ins(ps_ins, vec![ins_base]), Ev::Sub(2_000, other_base(r[2_000], s))];
    let g8 = apply(r, &ev8);
    let pm_placed = oracle_del(r, pm, 1, &inf);
    let codon_m = (13_499 - pm_placed) / 3 + 1;
    let ps_placed = oracle_del(r, ps_del, 1, &inf);
    let codon_s = 200 + (ps_placed - 16_000) / 3 + 1;
    // ---- genome 10: identical; 11: unrelated
    let g10 = r.clone();
    let g11 = random_dna(5_000, s);
    // ---- genome 12: lower case, U for T
    let g12: Vec<u8> = g1.iter().map(|&b| if b == b'T' { b'u' } else { b.to_ascii_lowercase() }).collect();
    let records = vec![
        ("g1_subs".to_string(), g1.clone()),
        ("g2_indels".to_string(), g2.clone()),
        ("g3_ends".to_string(), g3.clone()),
        ("g4_overhang".to_string(), g4.clone()),
        ("g5_missing".to_string(), g5.clone()),
        ("g6_revcomp".to_string(), g6.clone()),
        ("g7_repeats".to_string(), g7.clone()),
        ("g8_frameshift".to_string(), g8.clone()),
        ("g10_identical".to_string(), g10.clone()),
        ("g11_unrelated".to_string(), g11.clone()),
        ("g12_lowercase".to_string(), g12.clone()),
    ];
    let gpath = t.write("synth_genomes.fasta", &fasta(&records));
    let outdir = t.p("mut_out");
    let summary = t.p("mut_summary.tsv");
    let opt = MutationOptions { annotation: Some(gff.clone()), output_dir: Some(outdir.clone()), summary_output: Some(summary.clone()), split_orf1ab: true };
    let res = match viral_mutations(&refp, &[gpath.clone()], &opt) {
        Ok(x) => x,
        Err(e) => {
            t.check(false, format!("viral_mutations: {e}"));
            return;
        }
    };
    let gene_names: Vec<&str> = res.genes.iter().map(|g| g.name.as_str()).collect();
    t.check(gene_names == ["geneA", "ORF1a", "ORF1b", "geneM", "geneS"], format!("genes (ORF1ab split, pp1a once): {gene_names:?}"));
    let orf1b = res.genes.iter().find(|g| g.name == "ORF1b");
    t.check(orf1b.is_some_and(|g| g.parts == vec![(7000, 9999)] && g.codons == 1000), "ORF1b = the part after the slippage, 1000 codons");
    let gs = res.genes.iter().find(|g| g.name == "geneS");
    t.check(gs.is_some_and(|g| g.parts == vec![(15001, 15600), (16001, 16900)] && g.codons == 500), "a spliced CDS keeps its two parts in order");
    let by = |n: &str| res.genomes.iter().find(|g| g.name == n).cloned().unwrap_or_default();
    // genome 1
    let r1 = by("g1_subs");
    t.check(r1.aligned && r1.alignment_start == 1 && r1.alignment_end == L as u64 && near(r1.coverage, 1.0, 1e-12), "g1: aligned end to end");
    t.same_set("g1 substitutions", &labels_subs(&r1), &want_subs1);
    t.check(r1.deletions.is_empty() && r1.insertions.is_empty() && r1.missing.is_empty(), "g1: no indel, no missing");
    t.same_set("g1 amino-acid changes (codon by codon)", &labels_aa(&r1), &want_aa1);
    t.check(labels_aa(&r1).contains(&"geneA:Q500*".to_string()), "a stop gained: geneA:Q500*");
    t.check(labels_aa(&r1).contains(&"geneA:R600K".to_string()) && labels_aa(&r1).contains(&"geneA:G601R".to_string()), "three substitutions over two codons: R600K and G601R");
    let r600 = r1.aa_substitutions.iter().find(|a| a.label() == "geneA:R600K");
    t.check(r600.is_some_and(|a| a.nuc.len() == 2), format!("R600K is linked to its two nucleotide changes: {:?}", r600.map(|a| &a.nuc)));
    let sub7004 = r1.substitutions.iter().find(|x| x.pos == 7005);
    t.check(sub7004.is_some_and(|x| x.aa.iter().all(|a| a.starts_with("ORF1b:"))), "a change in the pp1a tail is reported in ORF1b only (split)");
    // genome 2
    let r2 = by("g2_indels");
    t.same_set("g2 deletions (placed as the oracle)", &labels_dels(&r2), &want_dels2);
    t.same_set("g2 insertions (placed as the oracle)", &labels_ins(&r2), &want_ins2);
    t.same_set("g2 amino acids: H69-/V70-, E156-/F157-/R158G, 214:EPE, 16:MPLF...", &labels_aa(&r2), &{
        let mut w = want_aa2.clone();
        w.extend(oracle_aa_subs(r, &apply(r, &[ev2[6].clone()]), &sy.genes));
        w
    });
    t.check(r2.frameshifts.is_empty(), format!("g2: no frameshift ({:?})", r2.frameshifts));
    // genome 3
    let r3 = by("g3_ends");
    t.check(r3.alignment_start == 55 && r3.alignment_end == (L - 70) as u64, format!("g3: covers {}-{}", r3.alignment_start, r3.alignment_end));
    t.check(r3.not_covered == vec![(1, 54), ((L - 69) as u64, L as u64)], format!("g3: not covered {:?}", r3.not_covered));
    t.check(r3.deletions.is_empty() && r3.insertions.is_empty(), "g3: missing ends are not deletions");
    let want3: Vec<String> = want_subs1.iter().filter(|l| {
        let p: usize = l[1..l.len() - 1].parse().unwrap();
        p >= 55 && p <= L - 70
    }).cloned().collect();
    t.same_set("g3 substitutions inside the covered part", &labels_subs(&r3), &want3);
    t.check(r3.unknown_aa.is_empty() || r3.unknown_aa.iter().all(|u| u.gene != "geneA" || u.codon_start > 0), "g3: unknown codons only at the ends");
    // genome 4
    let r4 = by("g4_overhang");
    t.check(r4.alignment_start == 1 && r4.alignment_end == L as u64 && r4.insertions.is_empty() && r4.substitutions.is_empty() && r4.not_covered.is_empty(), format!("g4: overhangs ignored ({:?} {:?})", r4.insertions, r4.substitutions.len()));
    // genome 5
    let r5 = by("g5_missing");
    t.check(r5.missing == vec![(2001, 2100), (11_000, 11_000), (20_000, 20_036)], format!("g5: missing {:?}", r5.missing));
    let amb: Vec<u64> = r5.ambiguous.iter().map(|a| a.0).collect();
    t.check(amb == vec![(k700 + 3) as u64, 25_000], format!("g5: ambiguous {:?}", r5.ambiguous));
    t.check(r5.substitutions.is_empty() && r5.aa_substitutions.is_empty(), format!("g5: N and ambiguity codes are not substitutions ({:?})", labels_aa(&r5)));
    let unk: Vec<String> = r5.unknown_aa.iter().map(|u| u.label()).collect();
    t.check(unk == vec!["geneA:634-667".to_string()], format!("g5: unknown codons {unk:?} (AAR is still K)"));
    t.check(near(r5.coverage, (L - 100 - 1 - 37) as f64 / L as f64, 1e-12), format!("g5: coverage {}", r5.coverage));
    // genome 6
    let r6 = by("g6_revcomp");
    t.check(r6.reverse_complement && labels_subs(&r6) == labels_subs(&r1) && labels_aa(&r6) == labels_aa(&r1), "g6: a reverse complement is turned around, same calls");
    // genome 7
    let r7 = by("g7_repeats");
    t.same_set("g7 deletions (homopolymer, 500 bases, 1 base)", &labels_dels(&r7), &want_dels7);
    t.check(labels_dels(&r7).contains(&"14998".to_string()), format!("g7: the A deleted from AAAA(TG) is placed before the gene (14998), not in its start codon: {:?}", labels_dels(&r7)));
    t.same_set("g7 insertions (tandem repeat, 300 bases)", &labels_ins(&r7), &want_ins7);
    t.check(r7.substitutions.len() == 1, format!("g7: the one substitution, no spurious ones ({:?})", labels_subs(&r7)));
    t.check(r7.aa_substitutions.is_empty() && r7.frameshifts.is_empty(), "g7: intergenic indels change no protein");
    // genome 8
    let r8 = by("g8_frameshift");
    let fs: Vec<String> = r8.frameshifts.iter().map(|f| f.label()).collect();
    t.check(fs.contains(&format!("geneM:{codon_m}-500")), format!("g8: geneM frameshift from codon {codon_m} to the end: {fs:?}"));
    t.check(fs.contains(&format!("geneS:{codon_s}-255")), format!("g8: geneS compensated frameshift {codon_s}-255: {fs:?}"));
    t.check(r8.aa_substitutions.iter().all(|a| a.gene == "geneA"), format!("g8: no amino-acid change called inside a frameshift: {:?}", labels_aa(&r8)));
    // genome 10, 11, 12
    let r10 = by("g10_identical");
    t.check(r10.aligned && r10.substitutions.is_empty() && r10.deletions.is_empty() && r10.insertions.is_empty() && labels_aa(&r10).is_empty() && r10.unknown_aa.is_empty(), "g10: the reference against itself: nothing");
    let r11 = by("g11_unrelated");
    t.check(!r11.aligned && r11.problem.as_deref().is_some_and(|p| p.contains("17-base")), format!("g11: an unrelated sequence is not aligned: {:?}", r11.problem));
    let r12 = by("g12_lowercase");
    t.check(labels_subs(&r12) == labels_subs(&r1) && labels_aa(&r12) == labels_aa(&r1), "g12: lower case and U read as A, C, G, T");
    // Outputs.
    let sum = std::fs::read_to_string(&summary).unwrap_or_default();
    t.check(sum.lines().count() == records.len() + 1 && sum.lines().all(|l| l.split('\t').count() == 26), "summary: one line per genome, 26 columns");
    t.check(sum.contains("del ") == false && sum.contains("geneA:16:MPLF"), "summary: Nextclade-style lists");
    let tsv2 = std::fs::read_to_string(std::path::Path::new(&outdir).join("g2_indels.tsv")).unwrap_or_default();
    t.check(tsv2.contains("aa_insertion\tgeneA:214:EPE") && tsv2.lines().all(|l| l.split('\t').count() == 8), "per-genome TSV with linked changes");
    t.check(r2.output.as_deref().is_some_and(|p| p.ends_with("g2_indels.tsv")), "each genome knows its table");
    let j = res.json();
    t.check(j.contains("\"genomes\"") && j.contains("\"text\"") && j.contains("geneA:H69-"), "json");
    t.check(res.resumo().contains("g2_indels") && res.resumo().contains("could not be aligned"), "summary text");

    // ORF1ab whole.
    let opt2 = MutationOptions { annotation: Some(gff.clone()), split_orf1ab: false, ..Default::default() };
    match viral_mutations(&refp, &[gpath.clone()], &opt2) {
        Ok(x) => {
            let names: Vec<&str> = x.genes.iter().map(|g| g.name.as_str()).collect();
            t.check(names == ["geneA", "ORF1ab", "ORF1a", "geneM", "geneS"], format!("split_orf1ab = false: {names:?}"));
            let g = x.genomes.iter().find(|g| g.name == "g1_subs").unwrap();
            t.same_set("ORF1ab whole: codon-by-codon changes", &labels_aa(g), &oracle_aa_subs(r, &g1, &sy.genes_whole));
        }
        Err(e) => {
            t.check(false, format!("split false: {e}"));
        }
    }
    // No annotation: nucleotides only.
    match viral_mutations(&refp, &[gpath.clone()], &MutationOptions::default()) {
        Ok(x) => {
            let g = x.genomes.iter().find(|g| g.name == "g2_indels").unwrap();
            t.check(x.genes.is_empty() && labels_aa(g).is_empty(), "no annotation: no amino acids");
            // Without genes there is no codon preference: every indel at its leftmost place.
            let no_frame = GapCost::default();
            let want: Vec<String> = ev2
                .iter()
                .filter_map(|e| if let Ev::Del(p, l) = e { Some(del_label(oracle_del(r, *p, *l, &no_frame), *l)) } else { None })
                .collect();
            t.same_set("no annotation: deletions at their leftmost place", &labels_dels(g), &want);
        }
        Err(e) => {
            t.check(false, format!("no annotation: {e}"));
        }
    }
    // A segmented reference: each genome goes to its segment.
    let seg_a = random_dna(4_000, s);
    let seg_b = random_dna(3_000, s);
    let segref = t.write("seg_ref.fasta", &fasta(&[("segA".into(), seg_a.clone()), ("segB".into(), seg_b.clone())]));
    let mut qb = seg_b.clone();
    qb[1500] = other_base(qb[1500], s);
    let mut qa = seg_a.clone();
    qa[100] = other_base(qa[100], s);
    let segq = t.write("seg_q.fasta", &fasta(&[("sampleB".into(), qb), ("sampleA".into(), rc(&qa))]));
    match viral_mutations(&segref, &[segq], &MutationOptions::default()) {
        Ok(x) => {
            let b = &x.genomes[0];
            let a = &x.genomes[1];
            t.check(b.reference == "segB" && labels_subs(b).len() == 1 && b.substitutions[0].pos == 1501, format!("segmented: sampleB on segB {:?}", labels_subs(b)));
            t.check(a.reference == "segA" && a.reverse_complement && a.substitutions.len() == 1 && a.substitutions[0].pos == 101, "segmented: sampleA on segA, turned around");
        }
        Err(e) => {
            t.check(false, format!("segmented: {e}"));
        }
    }
    // 300 genomes in one file (more than one batch of 256): every one
    // aligned, in the order of the file, each with its own substitution.
    let mut many = Vec::new();
    for i in 0..300 {
        let mut q = seg_a.clone();
        let p = 50 + i * 12;
        q[p] = other_base(q[p], s);
        many.push((format!("m{i:03}"), q));
    }
    let manyp = t.write("many.fasta", &fasta(&many));
    match viral_mutations(&segref, &[manyp], &MutationOptions::default()) {
        Ok(x) => {
            let ok = x.genomes.len() == 300
                && x.genomes.iter().enumerate().all(|(i, g)| g.name == format!("m{i:03}") && g.substitutions.len() == 1 && g.substitutions[0].pos == (51 + i * 12) as u64);
            t.check(ok, "300 genomes in batches: all aligned, in order, each with its substitution");
        }
        Err(e) => {
            t.check(false, format!("300 genomes: {e}"));
        }
    }
    // Errors.
    let other_gff = t.write("other.gff3", "##gff-version 3\nchrX\t.\tCDS\t1\t30\t.\t+\t0\tID=a;gene=g\nchrY\t.\tCDS\t1\t30\t.\t+\t0\tID=b;gene=h\n");
    t.fails("a GFF of another genome", viral_mutations(&refp, &[gpath.clone()], &MutationOptions { annotation: Some(other_gff), ..Default::default() }), &["chrX", "synth", "GFF3 of this reference"]);
    let long_gff = t.write("long.gff3", &format!("synth\t.\tCDS\t1\t{}\t.\t+\t0\tID=a;gene=g\n", L + 30));
    t.fails("a CDS past the end", viral_mutations(&refp, &[gpath.clone()], &MutationOptions { annotation: Some(long_gff), ..Default::default() }), &["long.gff3", "line 1", "past the end"]);
    let no_cds = t.write("nocds.gff3", "synth\t.\tgene\t1\t30\t.\t+\t.\tID=a\n");
    t.fails("a GFF without CDS", viral_mutations(&refp, &[gpath.clone()], &MutationOptions { annotation: Some(no_cds), ..Default::default() }), &["nocds.gff3", "no CDS"]);
    let bad_gff = t.write("bad.gff3", "synth CDS 1 30\n");
    t.fails("not a GFF", viral_mutations(&refp, &[gpath.clone()], &MutationOptions { annotation: Some(bad_gff), ..Default::default() }), &["bad.gff3", "line 1", "9 columns"]);
    let fq = t.write("ref.fq", "@r\nACGT\n+\nIIII\n");
    t.fails("a FASTQ reference", viral_mutations(&fq, &[gpath.clone()], &MutationOptions::default()), &["ref.fq", "FASTQ"]);
    t.fails("no genome", viral_mutations(&refp, &[], &MutationOptions::default()), &["No genome"]);
    t.fails("summary = reference", viral_mutations(&refp, &[gpath.clone()], &MutationOptions { summary_output: Some(refp.clone()), ..Default::default() }), &["also one of the input files"]);
}

// =====================================================================
//  3. IN-SILICO PCR
// =====================================================================

/// Naive primer match on text: the template base must be one the code
/// stands for; an N in the template never matches.
fn naive_match(tmpl: &[u8], primer: &[u8], max_mm: u32, exact3: usize) -> Option<u32> {
    let l = primer.len();
    let mut mm = 0;
    for i in 0..l {
        let ok = matches!(tmpl[i], b'A' | b'C' | b'G' | b'T') && iupac(primer[i]).contains(&tmpl[i]);
        if !ok {
            if i >= l - exact3.min(l) {
                return None;
            }
            mm += 1;
        }
    }
    if mm <= max_mm {
        Some(mm)
    } else {
        None
    }
}

fn rc_primer(p: &[u8]) -> Vec<u8> {
    p.iter()
        .rev()
        .map(|&b| match b {
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
            o => o,
        })
        .collect()
}

/// Every product, by brute force: (strand, start0, end0, fmm, rmm).
fn naive_products(seq: &[u8], f: &[u8], r: &[u8], opt: &PcrOptions, window: usize) -> BTreeSet<(char, usize, usize, u32, u32)> {
    let n = seq.len();
    let e3 = opt.three_prime_exact as usize;
    let mm = opt.max_mismatches;
    let rcs = rc(seq);
    let mut out = BTreeSet::new();
    let (lf, lr) = (f.len(), r.len());
    // F on +, R on - (its site on - read as written: R matches rcs at p')
    for p in 0..=n.saturating_sub(lf) {
        let Some(fm) = naive_match(&seq[p..p + lf], f, mm, e3) else { continue };
        for pp in 0..=n.saturating_sub(lr) {
            let Some(rm) = naive_match(&rcs[pp..pp + lr], r, mm, e3) else { continue };
            let q = n - pp - lr;
            if q >= p && q + lr >= p + lf && q + lr - p <= window && fm + rm <= opt.max_total_mismatches {
                out.insert(('+', p, q + lr, fm, rm));
            }
        }
    }
    for p in 0..=n.saturating_sub(lr) {
        let Some(rm) = naive_match(&seq[p..p + lr], r, mm, e3) else { continue };
        for pp in 0..=n.saturating_sub(lf) {
            let Some(fm) = naive_match(&rcs[pp..pp + lf], f, mm, e3) else { continue };
            let q = n - pp - lf;
            if q >= p && q + lf >= p + lr && q + lf - p <= window && fm + rm <= opt.max_total_mismatches {
                out.insert(('-', p, q + lf, fm, rm));
            }
        }
    }
    out
}

fn plant(seq: &mut [u8], at: usize, what: &[u8]) {
    seq[at..at + what.len()].copy_from_slice(what);
}

fn with_mismatches(p: &[u8], positions: &[usize]) -> Vec<u8> {
    let mut v = p.to_vec();
    for &i in positions {
        v[i] = match v[i] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
    }
    v
}

fn trial_pcr(t: &mut T, s: &mut u64) {
    let f = random_dna(20, s);
    let r = random_dna(22, s);
    let probe = random_dna(24, s);
    let rcr = rc(&r);
    let mut g = random_dna(20_000, s);
    // 1. + strand, exact, with the probe inside: 1000..1200
    plant(&mut g, 1000, &f);
    plant(&mut g, 1100, &probe);
    plant(&mut g, 1200 - 22, &rcr);
    // 2. - strand: R at 5000, rc(F) downstream
    plant(&mut g, 5000, &r);
    plant(&mut g, 5150, &rc(&f));
    // 3. two mismatches in F (not at 3'), R exact
    plant(&mut g, 8000, &with_mismatches(&f, &[2, 9]));
    plant(&mut g, 8150, &rcr);
    // 4. a mismatch at the 3' end of F: no product
    plant(&mut g, 11_000, &with_mismatches(&f, &[19]));
    plant(&mut g, 11_160, &rcr);
    // 5. overlapping products: F at 14000 and 14050, R at 14200
    plant(&mut g, 14_000, &f);
    plant(&mut g, 14_050, &f);
    plant(&mut g, 14_200, &rcr);
    // 6. too far for the default window: F at 16000, R at 18500
    plant(&mut g, 16_000, &f);
    plant(&mut g, 18_500 - 22, &rcr);
    // 7. 3 + 2 mismatches: over the total of 4
    plant(&mut g, 3000, &with_mismatches(&f, &[0, 5, 10]));
    plant(&mut g, 3200, &with_mismatches(&rcr, &[18, 20]));
    // 8. an N in the template inside F
    let mut fn_ = f.clone();
    fn_[4] = b'N';
    plant(&mut g, 7000, &fn_);
    plant(&mut g, 7100, &rcr);
    let gpath = t.write("pcr_g.fasta", &fasta(&[("contig1".into(), g.clone())]));
    // IUPAC primer: F with R and N codes.
    let mut f_iupac = f.clone();
    f_iupac[1] = if matches!(f[1], b'A' | b'G') { b'R' } else { b'Y' };
    f_iupac[6] = b'N';
    let assays = vec![
        PcrAssay { target: "T1".into(), forward: st(&f), reverse: st(&r), probe: Some(st(&probe)), min_size: Some(150), max_size: Some(250), source: "test".into() },
        PcrAssay { target: "T2".into(), forward: st(&f_iupac), reverse: st(&r), probe: None, min_size: None, max_size: None, source: String::new() },
    ];
    let opt = PcrOptions { amplicons_output: Some(t.p("amp.tsv")), matrix_output: Some(t.p("pcr_matrix.tsv")), ..Default::default() };
    let res = match in_silico_pcr(&assays, &[gpath.clone()], &opt) {
        Ok(x) => x,
        Err(e) => {
            t.check(false, format!("in_silico_pcr: {e}"));
            return;
        }
    };
    for (ai, a) in assays.iter().enumerate() {
        let want = naive_products(&g, a.forward.as_bytes(), a.reverse.as_bytes(), &opt, 2000);
        let got: BTreeSet<(char, usize, usize, u32, u32)> = res
            .amplicons
            .iter()
            .filter(|x| x.target == a.target)
            .map(|x| (x.strand, (x.start - 1) as usize, x.end as usize, x.forward_mismatches, x.reverse_mismatches))
            .collect();
        t.check(got == want, format!("{}: the products equal the brute-force search ({} vs {})", a.target, got.len(), want.len()));
        let _ = ai;
    }
    let t1: Vec<&Amplicon> = res.amplicons.iter().filter(|x| x.target == "T1").collect();
    let starts: Vec<(char, u64, u64)> = t1.iter().map(|x| (x.strand, x.start, x.end)).collect();
    t.check(starts.contains(&('+', 1001, 1200)), "product 1 on the + strand: 1001-1200");
    t.check(starts.contains(&('-', 5001, 5170)), format!("product 2 on the - strand: {starts:?}"));
    t.check(t1.iter().any(|x| x.start == 8001 && x.forward_mismatches == 2), "two mismatches away from the 3' end are allowed");
    t.check(!t1.iter().any(|x| x.start == 11_001), "a mismatch at the 3' end blocks the product");
    t.check(t1.iter().filter(|x| x.end == 14_222).count() == 2, "overlapping products of two forward sites");
    t.check(!t1.iter().any(|x| x.start == 16_001), "a product longer than the window is not reported");
    t.check(!t1.iter().any(|x| x.start == 3001), "3 + 2 mismatches exceed the total of 4");
    t.check(t1.iter().any(|x| x.start == 7001 && x.forward_mismatches == 1 && !x.in_range), "an N in the template counts as a mismatch (product of 122 bases, out of range)");
    let mut g_n3 = g.clone();
    g_n3[7000 + 18] = b'N';
    let gp_n3 = t.write("pcr_n3.fasta", &fasta(&[("contig1".into(), g_n3)]));
    if let Ok(x) = in_silico_pcr(&assays[..1], &[gp_n3], &PcrOptions::default()) {
        t.check(!x.amplicons.iter().any(|a| a.start == 7001), "an N in the last 3 bases of a primer site blocks it");
    }
    let p1 = t1.iter().find(|x| x.start == 1001).unwrap();
    t.check(p1.probe_mismatches == Some(0) && p1.in_range && p1.detected, "the probe binds product 1, in range: detected");
    t.check(p1.sequence.len() == 200 && p1.sequence.as_bytes()[..20] == f[..] && p1.sequence.as_bytes()[178..] == rcr[..], "the product sequence runs from F to rc(R)");
    let p2 = t1.iter().find(|x| x.start == 5001).unwrap();
    t.check(p2.sequence.as_bytes()[..20] == f[..] && !p2.detected && p2.probe_mismatches.is_none(), "a - strand product is given 5'->3' from F; no probe there: not detected");
    t.check(t1.iter().filter(|x| x.detected).count() == 1, "T1 detected once");
    t.check(res.rows.len() == 1 && res.rows[0].present == vec![true, true] && res.rows[0].counts[0] == 1, format!("matrix row {:?}", res.rows));
    // T2 (IUPAC) finds what T1 finds and more of nothing else.
    let t2: BTreeSet<(u64, u64)> = res.amplicons.iter().filter(|x| x.target == "T2").map(|x| (x.start, x.end)).collect();
    let t1s: BTreeSet<(u64, u64)> = t1.iter().map(|x| (x.start, x.end)).collect();
    t.check(t1s.is_subset(&t2), "an IUPAC primer finds the products of the plain one");
    let m = std::fs::read_to_string(t.p("pcr_matrix.tsv")).unwrap_or_default();
    t.check(m.starts_with("genome\ttype\tT1\tT2\npcr_g\tgenome\t1\t1"), format!("matrix file: {m}"));
    let a = std::fs::read_to_string(t.p("amp.tsv")).unwrap_or_default();
    t.check(a.lines().count() == res.amplicons.len() + 1, "amplicon table: one line per product");
    // Bigger window, larger total.
    let opt2 = PcrOptions { max_amplicon: 3000, max_total_mismatches: 5, ..Default::default() };
    if let Ok(x) = in_silico_pcr(&assays[..1], &[gpath.clone()], &opt2) {
        t.check(x.amplicons.iter().any(|a| a.start == 16_001 && !a.in_range), "with a 3000-base window the long product is found, out of range");
        t.check(x.amplicons.iter().any(|a| a.start == 3001), "with a total of 5 the 3+2 product is found");
        let want = naive_products(&g, &f, &r, &opt2, 3000);
        t.check(x.amplicons.len() == want.len(), "and still equals the brute force");
    }
    // Random primers on random genomes with planted mutated sites: brute force.
    let mut agree = 0;
    for k in 0..30 {
        let f = random_dna(18 + k % 5, s);
        let r = random_dna(18 + k % 4, s);
        let mut g = random_dna(3000, s);
        for _ in 0..6 {
            let at = rnd(s, 2900 - 30);
            let site = with_mismatches(&f, &[rnd(s, f.len()), rnd(s, f.len())]);
            plant(&mut g, at, &if lcg(s) % 2 == 0 { site } else { rc(&site) });
            let at = rnd(s, 2900 - 30);
            let site = with_mismatches(&r, &[rnd(s, r.len())]);
            plant(&mut g, at, &if lcg(s) % 2 == 0 { rc(&site) } else { site });
        }
        let gp = t.write("rand_pcr.fasta", &fasta(&[("c".into(), g.clone())]));
        let a = PcrAssay { target: "X".into(), forward: st(&f), reverse: st(&r), ..Default::default() };
        let o = PcrOptions { max_mismatches: 2, three_prime_exact: 2, max_total_mismatches: 3, ..Default::default() };
        if let Ok(x) = in_silico_pcr(&[a], &[gp], &o) {
            let got: BTreeSet<(char, usize, usize, u32, u32)> = x.amplicons.iter().map(|x| (x.strand, (x.start - 1) as usize, x.end as usize, x.forward_mismatches, x.reverse_mismatches)).collect();
            if got == naive_products(&g, &f, &r, &o, 2000) {
                agree += 1;
            }
        }
    }
    t.check(agree == 30, format!("random primers and planted sites: brute force agrees in {agree} of 30"));
    // Records as genomes.
    let multi = t.write("pcr_multi.fasta", &fasta(&[("rec1".into(), g.clone()), ("rec2".into(), random_dna(5000, s))]));
    if let Ok(x) = in_silico_pcr(&assays, &[multi.clone()], &PcrOptions { records_as_genomes: true, ..Default::default() }) {
        t.check(x.rows.len() == 2 && x.rows[0].name == "rec1" && x.rows[0].present[0] && !x.rows[1].present[0], "records as genomes: one row per record");
    }
    if let Ok(x) = in_silico_pcr(&assays, &[multi], &PcrOptions::default()) {
        t.check(x.rows.len() == 1 && x.rows[0].name == "pcr_multi" && x.rows[0].present[0], "a file as one genome");
    }
    // Reads: simulated 150-base reads, both strands; counts by brute force.
    let mut reads = Vec::new();
    for i in 0..3000 {
        let at = rnd(s, g.len() - 150);
        let mut rd = g[at..at + 150].to_vec();
        if i % 2 == 1 {
            rd = rc(&rd);
        }
        reads.push((format!("r{i}"), rd));
    }
    let fq: String = reads.iter().map(|(n, r)| format!("@{n}\n{}\n+\n{}\n", st(r), "I".repeat(r.len()))).collect();
    let fqp = t.write("reads.fastq", &fq);
    let opt3 = PcrOptions { reads_output: Some(t.p("reads.tsv")), min_reads: 2, ..Default::default() };
    match in_silico_pcr(&assays, &[gpath.clone(), fqp.clone()], &opt3) {
        Ok(x) => {
            let c = x.read_counts.iter().find(|c| c.target == "T1").unwrap();
            let e3 = 3;
            let has = |rd: &[u8], p: &[u8], mm: u32, e: usize| {
                let rr = rc(rd);
                (0..=rd.len() - p.len()).any(|i| naive_match(&rd[i..i + p.len()], p, mm, e).is_some() || naive_match(&rr[i..i + p.len()], p, mm, e).is_some())
            };
            let wf = reads.iter().filter(|(_, rd)| has(rd, &f, 3, e3)).count() as u64;
            let wr = reads.iter().filter(|(_, rd)| has(rd, &r, 3, e3)).count() as u64;
            let wp = reads.iter().filter(|(_, rd)| has(rd, &probe, 3, 0)).count() as u64;
            let wany = reads.iter().filter(|(_, rd)| has(rd, &f, 3, e3) || has(rd, &r, 3, e3) || has(rd, &probe, 3, 0)).count() as u64;
            let wamp = reads.iter().filter(|(_, rd)| !naive_products(rd, &f, &r, &opt3, 2000).is_empty()).count() as u64;
            t.check(c.reads == 3000 && c.forward == wf && c.reverse == wr && c.probe == wp && c.with_site == wany && c.amplicon == wamp,
                format!("reads: F {} vs {wf}, R {} vs {wr}, P {} vs {wp}, any {} vs {wany}, amplicon {} vs {wamp}", c.forward, c.reverse, c.probe, c.with_site, c.amplicon));
            t.check(c.present == (wany >= 2) && x.rows.len() == 2 && x.rows[1].kind == "reads", "a read sample is a row of the matrix");
            let rt = std::fs::read_to_string(t.p("reads.tsv")).unwrap_or_default();
            t.check(rt.lines().count() == 3, "reads table: one line per target");
            t.check(x.json().contains("\"read_counts\"") && x.resumo().contains("reads"), "json and summary");
        }
        Err(e) => {
            t.check(false, format!("reads: {e}"));
        }
    }
    // The panel file.
    let panel = t.write(
        "panel.tsv",
        "# a comment\nname\treverse\tforward\tprobe\tsize\tsource\nA1\tGCGCGACATTCCGAAGAA\tTTACAAACATTGGCCGCAAA\tACAATTTGCCCCCAGCGCTTCAG\t60-80\tLu 2020\nB2\tggactachvgggtwtctaat\tGTGYCAGCMGCCGCGGTAA\t\t\t\n",
    );
    match read_pcr_panel(&panel) {
        Ok(p) => {
            t.check(p.len() == 2 && p[0].target == "A1" && p[0].forward == "TTACAAACATTGGCCGCAAA" && p[0].min_size == Some(60) && p[0].max_size == Some(80) && p[0].source == "Lu 2020", format!("panel row 1: {:?}", p[0]));
            t.check(p[1].reverse == "GGACTACHVGGGTWTCTAAT" && p[1].probe.is_none() && p[1].min_size.is_none(), "panel row 2: upper-cased, no probe, no size");
        }
        Err(e) => {
            t.check(false, format!("panel: {e}"));
        }
    }
    let csv = t.write("panel.csv", "target,forward,reverse,min_size,max_size\nX,ACGTACGTACGTAC,TTTTGGGGCCCCAA,100,200\n");
    t.check(read_pcr_panel(&csv).is_ok_and(|p| p[0].max_size == Some(200)), "a CSV panel");
    let bad = t.write("panel_bad.tsv", "target\tforward\treverse\nX\tACGTJACGTACG\tACGTACGTACGT\n");
    t.fails("a letter that is not IUPAC", read_pcr_panel(&bad), &["panel_bad.tsv", "line 2", "\"J\"", "IUPAC"]);
    let nocol = t.write("panel_nocol.tsv", "name\tprimer1\tprimer2\nX\tACGTACGTAC\tACGTACGTAC\n");
    t.fails("missing columns", read_pcr_panel(&nocol), &["panel_nocol.tsv", "forward", "reverse"]);
    let range = t.write("panel_range.tsv", "target\tforward\treverse\tmin_size\tmax_size\nX\tACGTACGTAC\tACGTACGTAC\t300\t200\n");
    t.fails("min > max", read_pcr_panel(&range), &["300", "200", "minimum"]);
    let short = t.write("panel_short.tsv", "target\tforward\treverse\nX\tACGT\tACGTACGTAC\n");
    t.fails("too short", read_pcr_panel(&short), &["4 bases", "at least 8"]);
    t.fails("duplicate target", in_silico_pcr(&[assays[0].clone(), assays[0].clone()], &[gpath.clone()], &PcrOptions::default()), &["T1", "twice"]);
    t.fails("no input", in_silico_pcr(&assays, &[], &PcrOptions::default()), &["FASTA", "FASTQ"]);
    t.fails("output = input", in_silico_pcr(&assays, &[gpath.clone()], &PcrOptions { matrix_output: Some(gpath.clone()), ..Default::default() }), &["also one of the input files"]);
}

// ---------------------------------------------------------------------

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(String::as_str) == Some("real") {
        let (data, work) = (args.get(2).expect("research_data folder"), args.get(3).expect("work folder"));
        real::run(data, work, args.get(4).map(String::as_str), args.get(5).map(String::as_str), args.get(6).map(String::as_str));
        return;
    }
    let t0 = std::time::Instant::now();
    let dir = std::env::temp_dir().join(format!("trial_surto_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut t = T { ok: 0, total: 0, failed: Vec::new(), dir: dir.clone() };
    let mut s = 20261114u64;
    println!("distances");
    trial_distances(&mut t, &mut s);
    println!("trees");
    trial_trees(&mut t, &mut s);
    println!("viral mutations");
    trial_mutations(&mut t, &mut s);
    println!("in-silico PCR");
    trial_pcr(&mut t, &mut s);
    let _ = std::fs::remove_dir_all(&dir);
    println!("{} of {} checks passed ({:.1} s)", t.ok, t.total, t0.elapsed().as_secs_f64());
    if t.ok != t.total {
        std::process::exit(1);
    }
}

/// The real data of the research_data folder (see its MANIFEST.md).
mod real {
    use crate::surto::*;
    use std::collections::{BTreeSet, HashMap};
    use std::time::Instant;

    fn tsv(path: &str) -> Vec<HashMap<String, String>> {
        let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let mut lines = text.lines().filter(|l| !l.starts_with('#'));
        let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
        lines
            .filter(|l| !l.trim().is_empty())
            .map(|l| head.iter().zip(l.split('\t')).map(|(h, v)| (h.to_string(), v.to_string())).collect())
            .collect()
    }

    fn secs(t: Instant) -> String {
        format!("{:.1} s", t.elapsed().as_secs_f64())
    }

    fn set(s: &str) -> BTreeSet<String> {
        s.split(',').map(str::trim).filter(|x| !x.is_empty()).map(str::to_string).collect()
    }

    /// Markers of each lineage (S gene unless named): (present, absent).
    fn markers(lineage: &str) -> (Vec<&'static str>, Vec<&'static str>) {
        let core_b1 = vec!["S:D614G", "ORF1b:P314L"];
        let mut p: Vec<&'static str> = match lineage {
            "A" => vec!["ORF8:L84S"],
            "B.1" => vec![],
            "B.1.1.7" => vec!["S:N501Y", "S:A570D", "S:P681H", "S:T716I", "S:S982A", "S:D1118H", "N:R203K", "N:G204R"],
            "B.1.351" => vec!["S:D80A", "S:D215G", "S:K417N", "S:E484K", "S:N501Y", "S:A701V", "S:L241-", "S:L242-", "S:A243-"],
            "P.1" => vec!["S:L18F", "S:T20N", "S:P26S", "S:D138Y", "S:R190S", "S:K417T", "S:E484K", "S:N501Y", "S:H655Y", "S:T1027I"],
            "B.1.617.2" | "AY.4" => vec!["S:T19R", "S:L452R", "S:T478K", "S:P681R", "S:D950N", "S:E156-", "S:F157-", "S:R158G"],
            "B.1.427" => vec!["S:S13I", "S:W152C", "S:L452R"],
            "B.1.526" => vec!["S:T95I", "S:D253G"],
            "C.37" => vec!["S:G75V", "S:T76I", "S:L452Q", "S:F490S", "S:T859N", "S:R246-", "S:G252-"],
            "B.1.621" => vec!["S:T95I", "S:R346K", "S:E484K", "S:N501Y", "S:P681H", "S:D950N"],
            "BA.1" => vec!["S:A67V", "S:H69-", "S:V70-", "S:G339D", "S:S371L", "S:S373P", "S:K417N", "S:E484A", "S:N501Y", "S:Y505H", "S:P681H", "S:N969K", "S:214:EPE"],
            "BA.2" => vec!["S:T19I", "S:L24-", "S:A27S", "S:G339D", "S:S371F", "S:K417N", "S:E484A", "S:N501Y", "S:P681H"],
            "BA.5" => vec!["S:H69-", "S:L452R", "S:F486V", "S:N501Y"],
            "BQ.1.1" => vec!["S:R346T", "S:K444T", "S:L452R", "S:N460K", "S:F486V", "S:N501Y"],
            "XBB.1.5" => vec!["S:R346T", "S:N460K", "S:F486P", "S:F490S", "S:N501Y"],
            "EG.5.1" => vec!["S:Q52H", "S:F456L", "S:F486P", "S:N501Y"],
            "BA.2.86" => vec!["S:16:MPLF", "S:I332V", "S:K356T", "S:L452W", "S:N481K", "S:E484K", "S:N501Y", "S:P621S"],
            "JN.1" => vec!["S:L455S", "S:I332V", "S:N481K", "S:N501Y"],
            "KP.3.1.1" => vec!["S:L455S", "S:F456L", "S:Q493E", "S:N501Y", "S:V1104L"],
            _ => vec![],
        };
        if lineage == "A" {
            return (p, core_b1);
        }
        p.extend(core_b1);
        (p, vec![])
    }

    pub fn run(data: &str, work: &str, nextclade: Option<&str>, sm_kleb: Option<&str>, sm_sars: Option<&str>) {
        std::fs::create_dir_all(work).unwrap();
        let p = |rel: &str| format!("{data}/{rel}");
        let w = |name: &str| format!("{work}/{name}");

        // ------------------------------------------- 1. SARS-CoV-2 mutations
        let refp = p("reference_genomes/viruses/GCF_009858895.2_ASM985889v3_genomic.fna.gz");
        let gff = p("reference_genomes/viruses/GCF_009858895.2_ASM985889v3_genomic.gff.gz");
        let sars = p("outbreak/sars_cov_2/sars_cov_2_20_lineages.fasta");
        let meta: HashMap<String, (String, String)> = tsv(&p("outbreak/sars_cov_2/sars_cov_2_20_lineages.metadata.tsv"))
            .into_iter()
            .map(|r| (r["accession"].clone(), (r["pango_lineage"].clone(), r["who_label"].clone())))
            .collect();
        let t0 = Instant::now();
        let opt = MutationOptions {
            annotation: Some(gff.clone()),
            output_dir: Some(w("sars_mutations")),
            summary_output: Some(w("sars_mutations_summary.tsv")),
            split_orf1ab: true,
        };
        let res = viral_mutations(&refp, &[sars.clone()], &opt).unwrap_or_else(|e| panic!("{e}"));
        println!("== SARS-CoV-2: 20 genomes against NC_045512.2 ({}): genes {:?}", secs(t0), res.genes.iter().map(|g| g.name.clone()).collect::<Vec<_>>());
        for wv in &res.warnings {
            println!("   warning: {wv}");
        }
        println!("| genome | lineage | covers | subs | dels | ins | aa subs | aa dels | aa ins | markers found |");
        println!("| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |");
        for g in &res.genomes {
            let (lin, who) = meta.get(&g.name).cloned().unwrap_or_default();
            let labels: BTreeSet<String> = g
                .aa_substitutions
                .iter()
                .chain(&g.aa_deletions)
                .map(|a| a.label())
                .chain(g.aa_insertions.iter().map(|a| a.label()))
                .collect();
            let (present, absent) = markers(&lin);
            let found = present.iter().filter(|m| labels.contains(**m)).count();
            let missing: Vec<&&str> = present.iter().filter(|m| !labels.contains(**m)).collect();
            let wrongly: Vec<&&str> = absent.iter().filter(|m| labels.contains(**m)).collect();
            println!(
                "| {} | {}{} | {}-{} | {} | {} | {} | {} | {} | {} | {}/{}{}{} |",
                g.name,
                lin,
                if who == "-" || who.is_empty() { String::new() } else { format!(" ({who})") },
                g.alignment_start,
                g.alignment_end,
                g.substitutions.len(),
                g.deletions.len(),
                g.insertions.len(),
                g.aa_substitutions.len(),
                g.aa_deletions.len(),
                g.aa_insertions.len(),
                found,
                present.len(),
                if missing.is_empty() { String::new() } else { format!(" missing {missing:?}") },
                if absent.is_empty() { String::new() } else if wrongly.is_empty() { format!(" (and without {absent:?}, as expected)") } else { format!(" WRONGLY has {wrongly:?}") }
            );
        }
        if let Some(nc) = nextclade {
            compare_nextclade(&res, nc, &meta);
        }

        // ------------------------------------------- 2. distances and trees
        let kleb: Vec<HashMap<String, String>> = tsv(&p("outbreak/bacteria/genomes.tsv")).into_iter().filter(|r| r["organism"].contains("Klebsiella")).collect();
        let mut genomes = Vec::new();
        let mut st_of: HashMap<String, String> = HashMap::new();
        for r in &kleb {
            let strain: String = r["strain"].split(';').next().unwrap().trim().chars().map(|c| if c.is_ascii_alphanumeric() { c } else { '_' }).collect();
            let name = format!("{strain}_ST{}", r["MLST_ST"]);
            st_of.insert(name.clone(), r["MLST_ST"].clone());
            genomes.push(GenomeFiles { files: vec![p(&r["fasta"])], name });
        }
        let t0 = Instant::now();
        let opt = DistanceOptions {
            matrix_output: Some(w("kleb_distances.tsv")),
            tree_output: Some(w("kleb_nj.nwk")),
            ..Default::default()
        };
        let d = genome_distances(&genomes, &opt).unwrap_or_else(|e| panic!("{e}"));
        println!("\n== Klebsiella pneumoniae, {} genomes ({}): k = {}, scaled = {}", d.genomes.len(), secs(t0), d.k, d.scaled);
        print_ani(&d);
        println!("{}", d.tree.newick);
        print!("{}", d.tree.resumo());
        let cherries = sister_pairs(&d.tree);
        for (a, b) in [("KPNIH1_ST258", "30684_NJST258_2_ST258"), ("ST307PT01_ST307", "ST307_ST307"), ("PMK1_ST15", "ST15_ST15"), ("ATCC_BAA_2146_ST11", "HS11286_ST11")] {
            let ok = cherries.contains(&(a.to_string(), b.to_string())) || cherries.contains(&(b.to_string(), a.to_string()));
            println!("   {a} and {b} are sister leaves: {ok}  (D = {:.6}, ANI = {:.4}%)", dist_of(&d, a, b).0, 100.0 * dist_of(&d, a, b).1);
        }
        let upg = genome_distances(&genomes, &DistanceOptions { method: TreeMethod::Upgma, ..Default::default() }).unwrap();
        println!("   UPGMA: {}", upg.tree.newick);
        if let Some(csv) = sm_kleb {
            compare_sourmash(&d, csv);
        }
        let _ = st_of;

        let label = |n: &str| meta.get(n).map(|(l, wl)| if wl == "-" { format!("{n} {l}") } else { format!("{n} {l} {wl}") }).unwrap_or(n.to_string());
        let omicron: BTreeSet<String> = meta.iter().filter(|(_, (l, w))| w == "Omicron" && !l.is_empty()).map(|(a, _)| a.clone()).collect();
        let delta: BTreeSet<String> = meta.iter().filter(|(_, (_, w))| w == "Delta").map(|(a, _)| a.clone()).collect();
        for (outgroup, nwk) in [(None, "sars_nj_midpoint.nwk"), (Some("MN985325.1".to_string()), "sars_nj_rooted_on_A.nwk")] {
            let t0 = Instant::now();
            let opt = DistanceOptions {
                records_as_genomes: true,
                outgroup: outgroup.clone(),
                tree_output: Some(w(nwk)),
                matrix_output: Some(w("sars_distances.phy")),
                ..Default::default()
            };
            let d = genome_distances(&[GenomeFiles { files: vec![sars.clone()], name: String::new() }], &opt).unwrap_or_else(|e| panic!("{e}"));
            println!(
                "\n== SARS-CoV-2, 20 genomes ({}): k = {}, scaled = {}, rooted {}",
                secs(t0),
                d.k,
                d.scaled,
                outgroup.as_deref().map_or("at the midpoint".to_string(), |o| format!("on {o} (lineage A)"))
            );
            for l in d.tree.resumo().lines() {
                let t = l.trim_start();
                let ind = &l[..l.len() - t.len()];
                let mut parts = t.splitn(2, "  ");
                let nm = parts.next().unwrap_or("");
                let rest = parts.next().map(|x| format!("  {x}")).unwrap_or_default();
                println!("{ind}{}{rest}", label(nm));
            }
            let clades = clades_of(&d.tree);
            let all: BTreeSet<String> = meta.keys().cloned().collect();
            let rest: BTreeSet<String> = all.difference(&omicron).cloned().collect();
            println!(
                "   the Omicron genomes form one clade: {} (one side of a branch, rooting aside: {})",
                clades.contains(&omicron),
                clades.contains(&omicron) || clades.contains(&rest)
            );
            println!("   the two Delta genomes form one clade: {}", clades.contains(&delta));
            if outgroup.is_none() {
                if let Some(csv) = sm_sars {
                    compare_sourmash(&d, csv);
                }
            }
        }

        // ------------------------------------------- 3. in-silico PCR
        let panel_path = std::env::var("MMG_PANEL").unwrap_or_else(|_| w("demo_panel.tsv"));
        let panel = read_pcr_panel(&panel_path).unwrap_or_else(|e| panic!("{e}"));
        let t0 = Instant::now();
        let cdc: Vec<PcrAssay> = panel.iter().filter(|a| a.target.starts_with("CDC")).cloned().collect();
        let r = in_silico_pcr(&cdc, &[sars.clone()], &PcrOptions { records_as_genomes: true, amplicons_output: Some(w("pcr_sars_amplicons.tsv")), matrix_output: Some(w("pcr_sars_matrix.tsv")), ..Default::default() }).unwrap();
        println!("\n== CDC N1/N2 on the 20 SARS-CoV-2 genomes ({})", secs(t0));
        println!("| genome | lineage | N1 | N2 |");
        println!("| --- | --- | --- | --- |");
        for row in &r.rows {
            let cell = |t: &str| {
                r.amplicons
                    .iter()
                    .filter(|a| a.genome == row.name && a.target == t && a.detected)
                    .map(|a| format!("{}-{} ({} bp; mm F{} R{} P{})", a.start, a.end, a.size, a.forward_mismatches, a.reverse_mismatches, a.probe_mismatches.unwrap_or(99)))
                    .collect::<Vec<_>>()
                    .join("; ")
            };
            println!("| {} | {} | {} | {} |", row.name, meta.get(&row.name).map(|x| x.0.clone()).unwrap_or_default(), cell("CDC_N1"), cell("CDC_N2"));
        }
        let rr = in_silico_pcr(&cdc, &[refp.clone()], &PcrOptions::default()).unwrap();
        for a in &rr.amplicons {
            println!("   on NC_045512.2: {} {}-{} ({} bp), mismatches F{} R{} probe {:?}", a.target, a.start, a.end, a.size, a.forward_mismatches, a.reverse_mismatches, a.probe_mismatches);
        }

        let t0 = Instant::now();
        let refs = tsv(&p("reference_genomes/genomes.tsv"));
        let bact: Vec<(String, String)> = refs.iter().filter(|r| r["category"] == "bacterium").map(|r| (p(&r["fasta"]), r["organism"].clone())).collect();
        let s16: Vec<PcrAssay> = panel.iter().filter(|a| a.target.starts_with("16S")).cloned().collect();
        let files: Vec<String> = bact.iter().map(|b| b.0.clone()).collect();
        let r = in_silico_pcr(&s16, &files, &PcrOptions { amplicons_output: Some(w("pcr_16s_amplicons.tsv")), matrix_output: Some(w("pcr_16s_matrix.tsv")), ..Default::default() }).unwrap();
        println!("\n== 16S primers on {} reference bacteria ({})", files.len(), secs(t0));
        println!("| genome | organism | V4 products | full-length products | sizes |");
        println!("| --- | --- | ---: | ---: | --- |");
        for (row, (_, org)) in r.rows.iter().zip(&bact) {
            let sizes: BTreeSet<u64> = r.amplicons.iter().filter(|a| a.genome == row.name && a.detected).map(|a| a.size).collect();
            println!("| {} | {} | {} | {} | {:?} |", row.name, org.split_whitespace().take(3).collect::<Vec<_>>().join(" "), row.counts[0], row.counts[1], sizes);
        }

        let t0 = Instant::now();
        let mec: Vec<PcrAssay> = panel.iter().filter(|a| a.target == "mecA").cloned().collect();
        let staph = [
            ("outbreak/bacteria/staphylococcus_aureus_mrsa/GCF_000009645.1_ASM964v1_genomic.fna.gz", "N315 (MRSA)"),
            ("outbreak/bacteria/staphylococcus_aureus_mrsa/GCF_000011505.1_ASM1150v1_genomic.fna.gz", "MRSA252 (MRSA)"),
            ("reference_genomes/bacteria/GCF_000013465.1_ASM1346v1_genomic.fna.gz", "USA300_FPR3757 (MRSA)"),
            ("reference_genomes/bacteria/GCF_000013425.1_ASM1342v1_genomic.fna.gz", "NCTC 8325 (MSSA)"),
        ];
        let files: Vec<String> = staph.iter().map(|x| p(x.0)).collect();
        let r = in_silico_pcr(&mec, &files, &PcrOptions { matrix_output: Some(w("pcr_meca_matrix.tsv")), ..Default::default() }).unwrap();
        println!("\n== mecA (Murakami 1991) on S. aureus ({})", secs(t0));
        for (row, (_, lab)) in r.rows.iter().zip(&staph) {
            let a: Vec<String> = r.amplicons.iter().filter(|a| a.genome == row.name).map(|a| format!("{}:{}-{} {} bp mm F{} R{}", a.sequence_id, a.start, a.end, a.size, a.forward_mismatches, a.reverse_mismatches)).collect();
            println!("   {lab}: {}  {:?}", if row.present[0] { "positive" } else { "negative" }, a);
        }
    }

    fn print_ani(d: &DistanceResult) {
        let n = d.genomes.len();
        print!("| ANI % |");
        for i in 0..n {
            print!(" {} |", i + 1);
        }
        println!();
        print!("| --- |");
        for _ in 0..n {
            print!(" ---: |");
        }
        println!();
        for i in 0..n {
            print!("| {} {} |", i + 1, d.genomes[i].name);
            for j in 0..n {
                print!(" {:.2} |", 100.0 * d.ani[i][j]);
            }
            println!();
        }
    }

    fn dist_of(d: &DistanceResult, a: &str, b: &str) -> (f64, f64) {
        let i = d.genomes.iter().position(|g| g.name == a).unwrap_or(0);
        let j = d.genomes.iter().position(|g| g.name == b).unwrap_or(0);
        (d.distance[i][j], d.ani[i][j])
    }

    fn sister_pairs(t: &Tree) -> BTreeSet<(String, String)> {
        let mut out = BTreeSet::new();
        for n in &t.nodes {
            let leaves: Vec<&str> = n.children.iter().filter(|&&c| t.nodes[c].children.is_empty()).map(|&c| t.nodes[c].name.as_str()).collect();
            for a in 0..leaves.len() {
                for b in a + 1..leaves.len() {
                    out.insert((leaves[a].to_string(), leaves[b].to_string()));
                }
            }
        }
        out
    }

    fn clades_of(t: &Tree) -> Vec<BTreeSet<String>> {
        (0..t.nodes.len())
            .map(|u| {
                let mut s = BTreeSet::new();
                let mut st = vec![u];
                while let Some(x) = st.pop() {
                    if t.nodes[x].children.is_empty() {
                        s.insert(t.nodes[x].name.clone());
                    }
                    st.extend(t.nodes[x].children.iter().copied());
                }
                s
            })
            .collect()
    }

    fn compare_sourmash(d: &DistanceResult, csv: &str) {
        let (names, m) = match read_distance_matrix(csv) {
            Ok(x) => x,
            Err(e) => {
                println!("   sourmash matrix: {e}");
                return;
            }
        };
        let idx: HashMap<String, usize> = names.iter().enumerate().map(|(i, n)| (n.split_whitespace().next().unwrap_or("").to_string(), i)).collect();
        let mut worst = 0.0f64;
        let mut n = 0;
        let mut exact = 0;
        for i in 0..d.genomes.len() {
            for j in 0..d.genomes.len() {
                let (Some(&a), Some(&b)) = (idx.get(&d.genomes[i].name), idx.get(&d.genomes[j].name)) else { continue };
                let diff = (d.jaccard[i][j] - m[a][b]).abs();
                worst = worst.max(diff);
                n += 1;
                if diff < 1e-6 {
                    exact += 1;
                }
            }
        }
        println!("   sourmash compare ({csv}): {n} cells matched by name, {exact} equal to 1e-6, largest Jaccard difference {worst:.2e}");
    }

    fn compare_nextclade(res: &MutationResult, path: &str, meta: &HashMap<String, (String, String)>) {
        let rows = super::real::tsv(path);
        let by: HashMap<String, &HashMap<String, String>> = rows.iter().map(|r| (r["seqName"].split_whitespace().next().unwrap_or("").to_string(), r)).collect();
        println!("\n== Against Nextclade 3.29.0 (dataset nextstrain/sars-cov-2/wuhan-hu-1/orfs); ORF9b (not in the RefSeq GFF) and ORF10 (not in Nextclade) left out");
        println!("| genome | lineage | subs | dels | ins | missing | non-ACGTN | aa subs | aa dels | aa ins | start/end |");
        println!("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |");
        let mut all_diff: Vec<String> = Vec::new();
        let mut same_counts = [0usize; 9];
        for g in &res.genomes {
            let Some(r) = by.get(&g.name) else {
                println!("| {} | not in the Nextclade table |", g.name);
                continue;
            };
            let keep = |l: &String| !l.starts_with("ORF9b:") && !l.starts_with("ORF10:");
            let ours: [BTreeSet<String>; 8] = [
                g.substitutions.iter().map(|x| x.label()).collect(),
                g.deletions.iter().map(|x| x.short()).collect(),
                g.insertions.iter().map(|x| x.short()).collect(),
                g.missing.iter().map(|m| if m.0 == m.1 { m.0.to_string() } else { format!("{}-{}", m.0, m.1) }).collect(),
                g.ambiguous.iter().map(|(p, c)| format!("{c}:{p}")).collect(),
                g.aa_substitutions.iter().map(|x| x.label()).filter(keep).collect(),
                g.aa_deletions.iter().map(|x| x.label()).filter(keep).collect(),
                g.aa_insertions.iter().map(|x| x.label()).filter(keep).collect(),
            ];
            let cols = ["substitutions", "deletions", "insertions", "missing", "nonACGTNs", "aaSubstitutions", "aaDeletions", "aaInsertions"];
            let theirs: Vec<BTreeSet<String>> = cols.iter().map(|c| set(r.get(*c).map(String::as_str).unwrap_or("")).into_iter().filter(keep).collect()).collect();
            let mut cells = Vec::new();
            for k in 0..8 {
                let only_ours: Vec<&String> = ours[k].difference(&theirs[k]).collect();
                let only_theirs: Vec<&String> = theirs[k].difference(&ours[k]).collect();
                if only_ours.is_empty() && only_theirs.is_empty() {
                    cells.push(format!("= {}", ours[k].len()));
                    same_counts[k] += 1;
                } else {
                    cells.push(format!("{}/{} (+{} -{})", ours[k].len(), theirs[k].len(), only_ours.len(), only_theirs.len()));
                    all_diff.push(format!("{} {}: MMGEN only {:?}; Nextclade only {:?}", g.name, cols[k], only_ours, only_theirs));
                }
            }
            let se = format!("{}-{}", g.alignment_start, g.alignment_end);
            let se_nc = format!("{}-{}", r["alignmentStart"], r["alignmentEnd"]);
            if se == se_nc {
                same_counts[8] += 1;
            }
            println!(
                "| {} | {} | {} | {} |",
                g.name,
                meta.get(&g.name).map(|x| x.0.clone()).unwrap_or_default(),
                cells.join(" | "),
                if se == se_nc { format!("= {se}") } else { format!("{se} vs {se_nc}") }
            );
        }
        println!("   identical genome by genome: substitutions {}/20, deletions {}/20, insertions {}/20, missing {}/20, non-ACGTN {}/20, aa substitutions {}/20, aa deletions {}/20, aa insertions {}/20, alignment range {}/20",
            same_counts[0], same_counts[1], same_counts[2], same_counts[3], same_counts[4], same_counts[5], same_counts[6], same_counts[7], same_counts[8]);
        for d in all_diff {
            println!("   diff: {d}");
        }
    }
}
