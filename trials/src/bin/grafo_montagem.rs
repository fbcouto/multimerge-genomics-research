#![allow(dead_code)]
// Trials of src/grafo_montagem.rs: the compacted De Bruijn graph (unitigs)
// and the contig overlap graph that screens 3 and 7 draw, with their layout.
//
// The unitigs are compared with what `debruijn::montar_todos_contigs` (the
// assembler of screen 3) extracts; everything else with an oracle written
// here the plainest way, from the unitig sequences as strings: the edges by
// matching every (k-1)-base end against every (k-1)-base start, the degrees,
// the components, the kinds (tip, bubble, repeat, isolated), N50, and that
// every k-mer lies in exactly one unitig. The scaffolds of the overlap graph
// are compared with `scaffolding::fast_hash_scaffold` (screen 7) sequence by
// sequence.
//
// Usage: cargo run --release --bin grafo_montagem
//        cargo run --release --bin grafo_montagem -- <folder>
// where <folder> holds the UI examples (mmgen_examples/ written by
// make_examples) and the k-mers screen 1 makes of them, as raw bytes:
// toy_bacterium_K31.bin, toy_patient_K31.bin, toy_bacterium_reads_K31.bin.
// With GRAFO_JSON_DIR set, the JSON of each case is written there.
#[path = "../../../src/kmer.rs"]
mod kmer;
#[path = "../../../src/debruijn.rs"]
mod debruijn;
#[path = "../../../src/scaffolding.rs"]
mod scaffolding;
#[path = "../../../src/grafo_montagem.rs"]
mod grafo_montagem;

use grafo_montagem::*;
use serde_json::Value;
use std::collections::{BTreeMap, HashMap};
use std::time::Instant;

// ---------------------------------------------------------------------
//  The checks
// ---------------------------------------------------------------------

struct T {
    ok: usize,
    total: usize,
    failed: Vec<String>,
}

impl T {
    fn check(&mut self, cond: bool, what: impl AsRef<str>) -> bool {
        self.total += 1;
        if cond {
            self.ok += 1;
        } else {
            self.failed.push(what.as_ref().to_string());
            println!("  FAIL  {}", what.as_ref());
        }
        cond
    }
    fn eq<A: PartialEq + std::fmt::Debug>(&mut self, a: A, b: A, what: impl AsRef<str>) -> bool {
        let ok = a == b;
        if !ok {
            println!("        got {:?}, expected {:?}", a, b);
        }
        self.check(ok, what)
    }
}

// ---------------------------------------------------------------------
//  Sequences and k-mers
// ---------------------------------------------------------------------

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    fn chance(&mut self, p: f64) -> bool {
        ((self.next() >> 11) as f64 / (1u64 << 53) as f64) < p
    }
}

fn dna(n: usize, rng: &mut Rng) -> Vec<u8> {
    (0..n).map(|_| b"ACGT"[rng.below(4)]).collect()
}

fn other(b: u8, rng: &mut Rng) -> u8 {
    let o: Vec<u8> = b"ACGT".iter().copied().filter(|&c| c != b).collect();
    o[rng.below(3)]
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

fn encode(w: &[u8]) -> u128 {
    w.iter().fold(0u128, |x, &c| {
        (x << 2)
            | match c {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                _ => 3,
            }
    })
}

/// The k-mers of a sequence (the windows), as screen 1 reads them; a
/// circular sequence also gives the k-1 windows across its end.
fn kmers_of(s: &[u8], k: usize, circular: bool) -> Vec<u128> {
    let mut t = s.to_vec();
    if circular {
        t.extend_from_slice(&s[..k - 1]);
    }
    if t.len() < k {
        return Vec::new();
    }
    t.windows(k).map(encode).collect()
}

fn set_of(mut v: Vec<u128>) -> Vec<u128> {
    v.sort_unstable();
    v.dedup();
    v
}

fn json(s: &str) -> Value {
    serde_json::from_str(s).expect("valid JSON")
}

fn u(v: &Value) -> u64 {
    v.as_u64().unwrap_or(u64::MAX)
}

fn save(name: &str, text: &str) {
    if let Ok(dir) = std::env::var("GRAFO_JSON_DIR") {
        let _ = std::fs::create_dir_all(&dir);
        let _ = std::fs::write(format!("{dir}/{name}.json"), text);
    }
}

fn read_fasta(path: &str) -> Vec<String> {
    let text = std::fs::read_to_string(path).unwrap_or_default();
    let mut out: Vec<String> = Vec::new();
    for line in text.lines() {
        if line.starts_with('>') {
            out.push(String::new());
        } else if let Some(last) = out.last_mut() {
            last.push_str(line.trim());
        }
    }
    out
}

fn read_bin(path: &str) -> Option<Vec<u128>> {
    let b = std::fs::read(path).ok()?;
    Some(b.chunks_exact(16).map(|c| u128::from_ne_bytes(c.try_into().unwrap())).collect())
}

/// Reads drawn from a genome, `strands` 1 (forward only) or 2, with a
/// substitution rate of `err` per base.
fn reads(g: &[u8], n: usize, len: usize, err: f64, strands: usize, rng: &mut Rng) -> Vec<Vec<u8>> {
    (0..n)
        .map(|_| {
            let s = rng.below(g.len() - len + 1);
            let mut r = g[s..s + len].to_vec();
            if strands == 2 && rng.chance(0.5) {
                r = revcomp(&r);
            }
            for c in r.iter_mut() {
                if rng.chance(err) {
                    *c = other(*c, rng);
                }
            }
            r
        })
        .collect()
}

fn kmers_of_reads(rs: &[Vec<u8>], k: usize) -> Vec<u128> {
    set_of(rs.iter().flat_map(|r| kmers_of(r, k, false)).collect())
}

// ---------------------------------------------------------------------
//  The oracle: everything again from the unitig sequences, as strings
// ---------------------------------------------------------------------

/// A spelled sequence whose first k-1 bases are its last k-1 bases may be a
/// cycle that the assembler starts anywhere: rotate it to start at its
/// smallest k-mer, as both sides then agree.
fn canonical(s: &str, k: usize) -> String {
    let b = s.as_bytes();
    let k1 = k - 1;
    if b.len() < k || b[..k1] != b[b.len() - k1..] {
        return s.to_string();
    }
    let l = b.len() - k1;
    let cyc = &b[..l];
    // the cycle repeated (a cycle may be shorter than k: AAAAA is a circle of 1)
    let mut dd = Vec::new();
    while dd.len() < 2 * l + k + k1 {
        dd.extend_from_slice(cyc);
    }
    let mut best = 0;
    let mut best_v = u128::MAX;
    for p in 0..l {
        let v = encode(&dd[p..p + k]);
        if v < best_v {
            best_v = v;
            best = p;
        }
    }
    String::from_utf8(dd[best..best + l + k1].to_vec()).unwrap()
}

struct Oracle {
    edges: Vec<(usize, usize)>,
    ins: Vec<usize>,
    outs: Vec<usize>,
    comp: Vec<usize>,
    ncomp: usize,
    comp_size: Vec<usize>,
    circular: Vec<bool>,
    length: Vec<u64>,
    kind: Vec<&'static str>,
}

fn find(p: &mut [usize], mut x: usize) -> usize {
    while p[x] != x {
        p[x] = p[p[x]];
        x = p[x];
    }
    x
}

fn oracle(seqs: &[String], k: usize) -> Oracle {
    let n = seqs.len();
    let k1 = k - 1;
    let pre: Vec<&[u8]> = seqs.iter().map(|s| &s.as_bytes()[..k1]).collect();
    let suf: Vec<&[u8]> = seqs.iter().map(|s| &s.as_bytes()[s.len() - k1..]).collect();
    let mut by_pre: HashMap<&[u8], Vec<usize>> = HashMap::new();
    let mut by_suf: HashMap<&[u8], Vec<usize>> = HashMap::new();
    for i in 0..n {
        by_pre.entry(pre[i]).or_default().push(i);
        by_suf.entry(suf[i]).or_default().push(i);
    }
    let mut edges = Vec::new();
    for i in 0..n {
        if let Some(js) = by_pre.get(suf[i]) {
            for &j in js {
                edges.push((i, j));
            }
        }
    }
    edges.sort_unstable();
    let mut ins = vec![0; n];
    let mut outs = vec![0; n];
    let mut p: Vec<usize> = (0..n).collect();
    for &(a, b) in &edges {
        outs[a] += 1;
        ins[b] += 1;
        let (ra, rb) = (find(&mut p, a), find(&mut p, b));
        if ra != rb {
            p[ra] = rb;
        }
    }
    let mut ids: HashMap<usize, usize> = HashMap::new();
    let mut comp = vec![0; n];
    for i in 0..n {
        let r = find(&mut p, i);
        let next = ids.len();
        comp[i] = *ids.entry(r).or_insert(next);
    }
    let ncomp = ids.len();
    let mut comp_size = vec![0; ncomp];
    for &c in &comp {
        comp_size[c] += 1;
    }
    let self_edge: Vec<bool> = (0..n).map(|i| edges.binary_search(&(i, i)).is_ok()).collect();
    let circular: Vec<bool> = (0..n).map(|i| self_edge[i] && comp_size[comp[i]] == 1).collect();
    let length: Vec<u64> =
        (0..n).map(|i| if circular[i] { (seqs[i].len() - k1) as u64 } else { seqs[i].len() as u64 }).collect();
    let mut kind = vec![""; n];
    let mut groups: BTreeMap<(&[u8], &[u8]), Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        if !circular[i] && comp_size[comp[i]] > 1 && pre[i] != suf[i] && seqs[i].len() < 3 * k {
            groups.entry((pre[i], suf[i])).or_default().push(i);
        }
    }
    for g in groups.values() {
        if g.len() >= 2 {
            for &i in g {
                kind[i] = "bubble";
            }
        }
    }
    for i in 0..n {
        if comp_size[comp[i]] == 1 {
            kind[i] = "isolated";
            continue;
        }
        if kind[i] == "bubble" {
            continue;
        }
        let short = (seqs[i].len() as u64) < 2 * k as u64;
        let tip5 = ins[i] == 0 && outs[i] >= 1 && by_suf[suf[i]].len() >= 2;
        let tip3 = outs[i] == 0 && ins[i] >= 1 && by_pre[pre[i]].len() >= 2;
        if short && (tip5 || tip3) {
            kind[i] = "tip";
        }
    }
    // a repeat is entered by two different paths and left by two: tips do
    // not count, and the branches of a bubble (same start) count as one
    let mut reps = Vec::new();
    for i in 0..n {
        if !kind[i].is_empty() || ins[i] < 2 || outs[i] < 2 {
            continue;
        }
        let mut from: Vec<&[u8]> = edges.iter().filter(|e| e.1 == i && kind[e.0] != "tip").map(|e| pre[e.0]).collect();
        let mut to: Vec<&[u8]> = edges.iter().filter(|e| e.0 == i && kind[e.1] != "tip").map(|e| suf[e.1]).collect();
        from.sort();
        from.dedup();
        to.sort();
        to.dedup();
        if from.len() >= 2 && to.len() >= 2 {
            reps.push(i);
        }
    }
    for i in reps {
        kind[i] = "repeat";
    }
    Oracle { edges, ins, outs, comp, ncomp, comp_size, circular, length, kind }
}

fn n50_of(mut l: Vec<u64>) -> u64 {
    l.sort_unstable_by(|a, b| b.cmp(a));
    let total: u64 = l.iter().sum();
    let mut acc = 0;
    for x in l {
        acc += x;
        if 2 * acc >= total {
            return x;
        }
    }
    0
}

/// The whole graph against the oracle and against montar_todos_contigs.
fn against_oracle(t: &mut T, kmers: &[u128], k: usize, label: &str) -> Value {
    against_oracle_v(t, kmers, k, label, true)
}

fn against_oracle_v(t: &mut T, kmers: &[u128], k: usize, label: &str, verbose: bool) -> Value {
    let input = set_of(kmers.to_vec());
    let seqs = unitig_sequences(kmers, k as u8).unwrap();
    let g = json(&unitig_graph(kmers, k as u8, usize::MAX).unwrap());
    let nodes = g["nodes"].as_array().unwrap();
    let n = seqs.len();
    t.eq(nodes.len(), n, format!("{label}: one node per unitig"));

    // the unitigs are the assembler's contigs
    let contigs = debruijn::montar_todos_contigs(&input, k as u8);
    let mut a: Vec<String> = seqs.iter().map(|s| canonical(s, k)).collect();
    let mut b: Vec<String> = contigs.iter().map(|s| canonical(s, k)).collect();
    a.sort();
    b.sort();
    t.eq(a.len(), b.len(), format!("{label}: as many unitigs as montar_todos_contigs contigs ({})", b.len()));
    t.check(a == b, format!("{label}: the unitigs spell exactly the contigs of montar_todos_contigs"));

    // every k-mer in exactly one unitig
    let mut all: Vec<u128> = seqs.iter().flat_map(|s| kmers_of(s.as_bytes(), k, false)).collect();
    let total = all.len();
    all.sort_unstable();
    all.dedup();
    t.check(total == all.len() && all == input, format!("{label}: every k-mer lies in exactly one unitig ({total})"));

    let o = oracle(&seqs, k);
    let mut e: Vec<(usize, usize)> = g["edges"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| (u(&x["from"]) as usize, u(&x["to"]) as usize))
        .collect();
    e.sort_unstable();
    t.eq(e.len(), o.edges.len(), format!("{label}: number of edges"));
    t.check(e == o.edges, format!("{label}: the edges are exactly the (k-1)-overlaps between unitig ends"));
    t.check(
        g["edges"].as_array().unwrap().iter().all(|x| u(&x["overlap"]) == (k - 1) as u64),
        format!("{label}: every edge overlaps by k-1"),
    );
    let mut ok_deg = true;
    let mut ok_kind = true;
    let mut ok_len = true;
    let mut ok_circ = true;
    let mut ok_kmers = true;
    let mut bad = 0;
    for (i, nd) in nodes.iter().enumerate() {
        let id = u(&nd["id"]) as usize;
        if id != i {
            ok_deg = false;
            continue;
        }
        ok_deg &= u(&nd["in_degree"]) as usize == o.ins[i] && u(&nd["out_degree"]) as usize == o.outs[i];
        if nd["kind"].as_str().unwrap() != o.kind[i] {
            ok_kind = false;
            if bad < 3 {
                println!("        node {i}: kind {:?}, oracle {:?}", nd["kind"], o.kind[i]);
            }
            bad += 1;
        }
        ok_len &= u(&nd["length"]) == o.length[i];
        ok_circ &= nd["circular"].as_bool().unwrap() == o.circular[i];
        let windows = (seqs[i].len() - k + 1) as u64;
        ok_kmers &= u(&nd["kmers"]) == windows;
    }
    t.check(ok_deg, format!("{label}: in- and out-degrees"));
    t.check(ok_kind, format!("{label}: kinds (tip, bubble, repeat, isolated) as the oracle"));
    t.check(ok_len, format!("{label}: lengths"));
    t.check(ok_circ, format!("{label}: circular flags"));
    t.check(ok_kmers, format!("{label}: k-mers per unitig"));
    let s = &g["summary"];
    t.eq(u(&s["components"]) as usize, o.ncomp, format!("{label}: number of components"));
    // same partition: nodes share a component id exactly when the oracle's
    let mut pairs: HashMap<u64, usize> = HashMap::new();
    let mut ok_part = true;
    for (i, nd) in nodes.iter().enumerate() {
        let c = u(&nd["component"]);
        match pairs.get(&c) {
            Some(&oc) => ok_part &= oc == o.comp[i],
            None => {
                pairs.insert(c, o.comp[i]);
            }
        }
    }
    ok_part &= pairs.len() == o.ncomp;
    t.check(ok_part, format!("{label}: components are the oracle's"));
    // components listed from the largest down, with their sizes
    let comps = g["components"].as_array().unwrap();
    let mut ok_cs = true;
    let mut last = u64::MAX;
    for c in comps {
        let len = u(&c["length"]);
        ok_cs &= len <= last;
        last = len;
        let id = u(&c["id"]);
        let members: Vec<usize> = (0..n).filter(|&i| u(&nodes[i]["component"]) == id).collect();
        ok_cs &= members.len() as u64 == u(&c["nodes"]);
        let kms: u64 = members.iter().map(|&i| u(&nodes[i]["kmers"])).sum();
        let sources = members.iter().filter(|&&i| o.ins[i] == 0).count() as u64;
        ok_cs &= len == kms + (k as u64 - 1) * sources;
    }
    t.check(ok_cs, format!("{label}: component lengths, sizes and order"));
    t.eq(u(&s["n50"]), n50_of(o.length.clone()), format!("{label}: N50"));
    t.eq(u(&s["tips"]) as usize, o.kind.iter().filter(|&&x| x == "tip").count(), format!("{label}: tip count"));
    t.eq(
        u(&s["repeats"]) as usize,
        o.kind.iter().filter(|&&x| x == "repeat").count(),
        format!("{label}: repeat count"),
    );
    if verbose { println!(
        "  {label}: {} k-mers, {} unitigs, {} edges, {} components; tips {}, bubbles {} ({} branches), repeats {}, isolated {}, circular {}",
        input.len(),
        n,
        e.len(),
        o.ncomp,
        s["tips"],
        s["bubbles"],
        s["bubble_branches"],
        s["repeats"],
        s["isolated"],
        s["circular_nodes"]
    ); }
    g
}

// ---------------------------------------------------------------------
//  Layout checks
// ---------------------------------------------------------------------

fn all_finite(g: &Value) -> bool {
    let mut ok = true;
    for nd in g["nodes"].as_array().unwrap() {
        ok &= nd["x"].as_f64().is_some_and(f64::is_finite) && nd["y"].as_f64().is_some_and(f64::is_finite);
        for p in nd["path"].as_array().unwrap() {
            ok &= p[0].as_f64().is_some_and(f64::is_finite) && p[1].as_f64().is_some_and(f64::is_finite);
        }
        ok &= !nd["path"].as_array().unwrap().is_empty();
    }
    ok &= g["layout"]["width"].as_f64().is_some_and(f64::is_finite);
    ok
}

/// (the smallest distance between two node centres, between two path points
/// of different nodes), in layout units
fn min_distances(g: &Value) -> (f64, f64) {
    let unit = g["layout"]["unit"].as_f64().unwrap();
    let nodes = g["nodes"].as_array().unwrap();
    let mut c: Vec<(f64, f64)> =
        nodes.iter().map(|n| (n["x"].as_f64().unwrap(), n["y"].as_f64().unwrap())).collect();
    c.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut best_c = f64::INFINITY;
    for i in 0..c.len() {
        for j in (i + 1)..c.len() {
            if c[j].0 - c[i].0 >= best_c {
                break;
            }
            best_c = best_c.min(((c[j].0 - c[i].0).powi(2) + (c[j].1 - c[i].1).powi(2)).sqrt());
        }
    }
    let mut pts: Vec<(f64, f64, usize)> = Vec::new();
    for (ni, n) in nodes.iter().enumerate() {
        for p in n["path"].as_array().unwrap() {
            pts.push((p[0].as_f64().unwrap(), p[1].as_f64().unwrap(), ni));
        }
    }
    pts.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut best_p = f64::INFINITY;
    for i in 0..pts.len() {
        for j in (i + 1)..pts.len() {
            if pts[j].0 - pts[i].0 >= best_p {
                break;
            }
            if pts[i].2 != pts[j].2 {
                best_p = best_p.min(((pts[j].0 - pts[i].0).powi(2) + (pts[j].1 - pts[i].1).powi(2)).sqrt());
            }
        }
    }
    (best_c / unit, best_p / unit)
}

/// The points of a path lie on one circle (relative spread of the radius).
fn on_a_circle(path: &[Value]) -> bool {
    let p: Vec<(f64, f64)> = path.iter().map(|q| (q[0].as_f64().unwrap(), q[1].as_f64().unwrap())).collect();
    let (x0, x1) = p.iter().fold((f64::MAX, f64::MIN), |a, q| (a.0.min(q.0), a.1.max(q.0)));
    let (y0, y1) = p.iter().fold((f64::MAX, f64::MIN), |a, q| (a.0.min(q.1), a.1.max(q.1)));
    let (cx, cy) = ((x0 + x1) / 2.0, (y0 + y1) / 2.0);
    let r: Vec<f64> = p.iter().map(|q| ((q.0 - cx).powi(2) + (q.1 - cy).powi(2)).sqrt()).collect();
    let rm = r.iter().sum::<f64>() / r.len() as f64;
    p.len() >= 8 && r.iter().all(|&x| (x - rm).abs() < 0.03 * rm + 0.2)
}

fn layout_checks(t: &mut T, g: &Value, label: &str) {
    t.check(all_finite(g), format!("{label}: every coordinate is a finite number"));
    if g["nodes"].as_array().unwrap().len() > 1 {
        let (c, p) = min_distances(g);
        t.check(c > 0.2, format!("{label}: node centres apart (closest {c:.2} units)"));
        t.check(p >= 0.37, format!("{label}: no two nodes' points on top of each other (closest {p:.2} units)"));
    }
}

// ---------------------------------------------------------------------
//  Scaffolds of the overlap graph against fast_hash_scaffold
// ---------------------------------------------------------------------

fn scaffolds_match(t: &mut T, seqs: &[String], m: usize, label: &str) -> Value {
    let g = json(&overlap_graph_with_limit(seqs, m, usize::MAX).unwrap());
    let (want, circ) = scaffolding::fast_hash_scaffold(seqs.to_vec(), m);
    let mut got = Vec::new();
    for s in g["scaffolds"].as_array().unwrap() {
        let ids: Vec<usize> = s["contigs"].as_array().unwrap().iter().map(|x| u(x) as usize).collect();
        let mut q = seqs[ids[0]].clone();
        for &i in &ids[1..] {
            q.push_str(&seqs[i][m..]);
        }
        if s["circular"].as_bool().unwrap() {
            q.truncate(q.len() - m);
        }
        got.push(q);
    }
    let ok = got == want;
    let lens: Vec<u64> = g["scaffolds"].as_array().unwrap().iter().map(|s| u(&s["length"])).collect();
    let ok_len = lens == want.iter().map(|s| s.len() as u64).collect::<Vec<_>>();
    let ok_c = u(&g["summary"]["circular_scaffolds"]) as usize == circ;
    t.check(ok && ok_len && ok_c, format!("{label}: scaffolds equal to fast_hash_scaffold's ({} of them, {circ} circular)", want.len()));
    g
}

/// Every overlap of the graph is a real suffix-prefix match, the longest one,
/// and every pair with such a match of at least m bases is an edge.
fn overlaps_brute(t: &mut T, seqs: &[String], m: usize, g: &Value, label: &str) {
    let mut want = Vec::new();
    for (i, a) in seqs.iter().enumerate() {
        for (j, b) in seqs.iter().enumerate() {
            let (a, b) = (a.as_bytes(), b.as_bytes());
            let top = if i == j { a.len() - 1 } else { a.len().min(b.len()) };
            for l in (m..=top).rev() {
                if a[a.len() - l..] == b[..l] {
                    let rule = a.len() >= m && b.len() >= m && a[a.len() - m..] == b[..m];
                    want.push((i, j, l, rule));
                    break;
                }
            }
        }
    }
    let mut got: Vec<(usize, usize, usize, bool)> = g["edges"]
        .as_array()
        .unwrap()
        .iter()
        .map(|e| (u(&e["from"]) as usize, u(&e["to"]) as usize, u(&e["overlap"]) as usize, e["scaffold_rule"].as_bool().unwrap()))
        .collect();
    got.sort_unstable();
    t.check(got == want, format!("{label}: overlaps equal to a brute-force search ({} edges)", want.len()));
}

// ---------------------------------------------------------------------
//  The cases
// ---------------------------------------------------------------------

/// A regression bound for the timings, in seconds (best of 3, wall clock). On
/// an idle i7-7700HQ 5,000 nodes take ~0.3 s; the bound is loose because
/// other jobs on a shared machine can triple the wall time.
const LIMIT: f64 = 2.0;

fn best_of_3(mut f: impl FnMut() -> String) -> (String, f64) {
    let mut best = f64::MAX;
    let mut out = String::new();
    for _ in 0..3 {
        let t0 = Instant::now();
        out = f();
        best = best.min(t0.elapsed().as_secs_f64());
    }
    (out, best)
}

fn node<'a>(g: &'a Value, id: usize) -> &'a Value {
    g["nodes"].as_array().unwrap().iter().find(|n| u(&n["id"]) == id as u64).unwrap()
}

fn kinds(g: &Value) -> BTreeMap<String, usize> {
    let mut m = BTreeMap::new();
    for n in g["nodes"].as_array().unwrap() {
        *m.entry(n["kind"].as_str().unwrap().to_string()).or_insert(0) += 1;
    }
    m
}

fn synthetic(t: &mut T) {
    let mut rng = Rng(2026);
    let k = 31usize;

    println!("== a linear genome");
    let genome = dna(30_000, &mut rng);
    let km = kmers_of(&genome, k, false);
    let s = unitig_graph(&km, k as u8, 2000).unwrap();
    save("linear", &s);
    let g = json(&s);
    t.eq(u(&g["summary"]["nodes"]), 1, "linear: one unitig");
    t.eq(u(&g["summary"]["edges"]), 0, "linear: no edge");
    let n0 = node(&g, 0);
    t.eq(u(&n0["length"]), 30_000, "linear: 30,000 bases");
    t.eq(u(&n0["kmers"]), 29_970, "linear: 29,970 k-mers");
    t.eq(n0["kind"].as_str(), Some("isolated"), "linear: kind isolated");
    t.eq(n0["circular"].as_bool(), Some(false), "linear: not circular");
    t.eq(u(&g["summary"]["n50"]), 30_000, "linear: N50");
    let gc = genome.iter().filter(|&&c| c == b'G' || c == b'C').count() as f64 / 30_000.0;
    t.check((n0["gc"].as_f64().unwrap() - gc).abs() < 1e-4, "linear: GC content");
    t.eq(unitig_sequences(&km, k as u8).unwrap(), vec![String::from_utf8(genome.clone()).unwrap()], "linear: the unitig is the genome");
    layout_checks(t, &g, "linear");

    println!("== a circular plasmid");
    let plasmid = dna(20_000, &mut rng);
    let kp = kmers_of(&plasmid, k, true);
    let s = unitig_graph(&kp, k as u8, 2000).unwrap();
    save("plasmid", &s);
    let g = json(&s);
    let n0 = node(&g, 0);
    t.eq(u(&g["summary"]["nodes"]), 1, "plasmid: one unitig");
    t.eq(n0["circular"].as_bool(), Some(true), "plasmid: circular");
    t.eq(u(&n0["length"]), 20_000, "plasmid: 20,000 bases around");
    t.eq(u(&n0["kmers"]), 20_000, "plasmid: 20,000 k-mers");
    t.eq((u(&n0["in_degree"]), u(&n0["out_degree"])), (1, 1), "plasmid: a self-loop");
    t.eq(n0["kind"].as_str(), Some("isolated"), "plasmid: kind isolated");
    t.eq(g["components"][0]["circular"].as_bool(), Some(true), "plasmid: a circular component");
    t.eq(u(&g["components"][0]["length"]), 20_000, "plasmid: component length 20,000");
    t.check(on_a_circle(n0["path"].as_array().unwrap()), "plasmid: drawn as a ring");
    let sq = unitig_sequences(&kp, k as u8).unwrap();
    t.eq(sq[0].len(), 20_000 + k - 1, "plasmid: spelled with the k-1 bases repeated, as assemble_debruijn");
    t.check(canonical(&sq[0], k) == sq[0], "plasmid: starts at its smallest k-mer");
    layout_checks(t, &g, "plasmid");

    println!("== a repeat longer than k");
    let (a, r, mut b, mut c) = (dna(3000, &mut rng), dna(100, &mut rng), dna(2000, &mut rng), dna(3000, &mut rng));
    // the bases around the two copies differ, so the repeat is exactly R
    if b[1999] == a[2999] {
        b[1999] = other(a[2999], &mut rng);
    }
    if c[0] == b[0] {
        c[0] = other(b[0], &mut rng);
    }
    let rep: Vec<u8> = [&a[..], &r, &b, &r, &c].concat();
    let s = unitig_graph(&kmers_of(&rep, k, false), k as u8, 2000).unwrap();
    save("repeat", &s);
    let g = json(&s);
    t.eq(u(&g["summary"]["nodes"]), 4, "repeat: 4 unitigs");
    t.eq(u(&g["summary"]["edges"]), 4, "repeat: 4 edges");
    let lens: Vec<u64> = g["nodes"].as_array().unwrap().iter().map(|n| u(&n["length"])).collect();
    t.eq(lens, vec![3030, 3030, 2060, 100], "repeat: lengths A+30, C+30, 30+B+30, R");
    let rn = node(&g, 3);
    t.eq((u(&rn["in_degree"]), u(&rn["out_degree"])), (2, 2), "repeat: the repeat unitig has in = out = 2");
    t.eq(rn["kind"].as_str(), Some("repeat"), "repeat: tagged repeat");
    t.eq(u(&g["summary"]["repeats"]), 1, "repeat: one repeat");
    t.check(
        (0..3).all(|i| node(&g, i)["kind"].as_str() == Some("")),
        "repeat: the other three untagged",
    );
    layout_checks(t, &g, "repeat");
    let short: Vec<u8> = [&a[..], &r[..20], &b, &r[..20], &c].concat();
    let g = json(&unitig_graph(&kmers_of(&short, k, false), k as u8, 2000).unwrap());
    t.eq(u(&g["summary"]["nodes"]), 1, "a repeat shorter than k does not branch the graph");

    println!("== SNP bubbles (the healthy genome and a patient together)");
    let mut patient = genome.clone();
    patient[15_000] = other(patient[15_000], &mut rng);
    let union = set_of([kmers_of(&genome, k, false), kmers_of(&patient, k, false)].concat());
    let s = unitig_graph(&union, k as u8, 2000).unwrap();
    save("snp1", &s);
    let g = json(&s);
    t.eq(u(&g["summary"]["nodes"]), 4, "1 SNP: 4 unitigs");
    t.eq(u(&g["summary"]["bubbles"]), 1, "1 SNP: 1 bubble");
    t.eq(u(&g["summary"]["bubble_branches"]), 2, "1 SNP: 2 branches");
    let br: Vec<&Value> = g["nodes"].as_array().unwrap().iter().filter(|n| n["kind"] == "bubble").collect();
    t.check(br.iter().all(|n| u(&n["length"]) == 2 * k as u64 - 1 && u(&n["kmers"]) == k as u64), "1 SNP: each branch has k k-mers, 2k-1 bases");
    t.eq(kinds(&g).get("").copied(), Some(2), "1 SNP: the flanks untagged");
    layout_checks(t, &g, "1 SNP");
    let mut patient3 = genome.clone();
    for p in [5_000, 15_000, 25_000] {
        patient3[p] = other(patient3[p], &mut rng);
    }
    let union3 = set_of([kmers_of(&genome, k, false), kmers_of(&patient3, k, false)].concat());
    let s = unitig_graph(&union3, k as u8, 2000).unwrap();
    save("snp3", &s);
    let g = json(&s);
    t.eq(u(&g["summary"]["nodes"]), 10, "3 SNPs: 10 unitigs");
    t.eq((u(&g["summary"]["bubbles"]), u(&g["summary"]["bubble_branches"])), (3, 6), "3 SNPs: 3 bubbles, 6 branches");
    t.eq(u(&g["summary"]["components"]), 1, "3 SNPs: one component");
    t.eq(u(&g["summary"]["repeats"]), 0, "3 SNPs: the pieces between bubbles are not repeats");
    against_oracle(t, &union3, k, "3 SNPs");
    layout_checks(t, &g, "3 SNPs");

    // two long different stretches between the same two shared pieces:
    // parallel paths, but not a bubble (branches of 3k bases or more)
    let (s0, l1, l2, e0) = (dna(2000, &mut rng), dna(1000, &mut rng), dna(1200, &mut rng), dna(2000, &mut rng));
    let p1: Vec<u8> = [&s0[..], &l1, &e0].concat();
    let p2: Vec<u8> = [&s0[..], &l2, &e0].concat();
    let g = json(&unitig_graph(&set_of([kmers_of(&p1, k, false), kmers_of(&p2, k, false)].concat()), k as u8, 2000).unwrap());
    t.eq(u(&g["summary"]["nodes"]), 4, "long parallel paths: 4 unitigs");
    t.eq(u(&g["summary"]["bubbles"]), 0, "long parallel paths (1,000 and 1,200 bases): not a bubble");

    println!("== tips from reads with an error near an end");
    let g10 = dna(10_000, &mut rng);
    let mut r3 = g10[2000..2040].to_vec();
    r3[35] = other(r3[35], &mut rng);
    let mut r5 = g10[6000..6040].to_vec();
    r5[4] = other(r5[4], &mut rng);
    let km = set_of([kmers_of(&g10, k, false), kmers_of(&r3, k, false), kmers_of(&r5, k, false)].concat());
    let s = unitig_graph(&km, k as u8, 2000).unwrap();
    save("tips", &s);
    let g = json(&s);
    t.eq(u(&g["summary"]["nodes"]), 5, "tips: 3 genome pieces + 2 tips");
    t.eq(u(&g["summary"]["tips"]), 2, "tips: 2 tips");
    let tips: Vec<&Value> = g["nodes"].as_array().unwrap().iter().filter(|n| n["kind"] == "tip").collect();
    t.check(tips.iter().all(|n| u(&n["length"]) == 35 && u(&n["kmers"]) == 5), "tips: each of 5 k-mers, 35 bases");
    let mut ends: Vec<(u64, u64)> = tips.iter().map(|n| (u(&n["in_degree"]), u(&n["out_degree"]))).collect();
    ends.sort_unstable();
    t.eq(ends, vec![(0, 1), (1, 0)], "tips: one 5' tip (in 0) and one 3' tip (out 0)");
    t.eq(kinds(&g).get("").copied(), Some(3), "tips: the genome pieces untagged");
    layout_checks(t, &g, "tips");

    println!("== a chromosome and two plasmids");
    let chrom = dna(40_000, &mut rng);
    let (p1, p2) = (dna(8_000, &mut rng), dna(5_000, &mut rng));
    let km = set_of([kmers_of(&chrom, k, false), kmers_of(&p1, k, true), kmers_of(&p2, k, true)].concat());
    let s = unitig_graph(&km, k as u8, 2000).unwrap();
    save("mixture", &s);
    let g = json(&s);
    t.eq(u(&g["summary"]["components"]), 3, "mixture: 3 components");
    t.eq(u(&g["summary"]["circular_components"]), 2, "mixture: 2 circular");
    let cl: Vec<(u64, bool)> =
        g["components"].as_array().unwrap().iter().map(|c| (u(&c["length"]), c["circular"].as_bool().unwrap())).collect();
    t.eq(cl, vec![(40_000, false), (8_000, true), (5_000, true)], "mixture: lengths and rings, largest first");
    let rings = g["nodes"].as_array().unwrap().iter().filter(|n| n["circular"] == true && on_a_circle(n["path"].as_array().unwrap())).count();
    t.eq(rings, 2, "mixture: the two plasmids drawn as rings");
    layout_checks(t, &g, "mixture");

    println!("== reads with errors, both strands, against the oracle");
    let g20 = dna(20_000, &mut rng);
    let rs = reads(&g20, 1500, 150, 0.005, 2, &mut rng);
    let km_reads = kmers_of_reads(&rs, k);
    let g = against_oracle(t, &km_reads, k, "reads k=31");
    t.check(u(&g["summary"]["tips"]) > 0 && u(&g["summary"]["bubbles"]) > 0, "reads: errors make tips and bubbles");
    save("reads", &serde_json::to_string(&g).unwrap());
    layout_checks(t, &g, "reads k=31");
    for (kk, glen, label) in [(21usize, 8_000usize, "reads k=21"), (5, 1_500, "k=5"), (4, 600, "k=4"), (3, 300, "k=3")] {
        let gg = dna(glen, &mut rng);
        let rs = reads(&gg, if kk > 10 { 400 } else { 30 }, if kk > 10 { 100 } else { 40 }, 0.01, 2, &mut rng);
        let mut km = kmers_of_reads(&rs, kk);
        km.extend(kmers_of(&gg, kk, true));
        let km = set_of(km);
        against_oracle(t, &km, kk, label);
    }
    // many small random cases, tangled ones included
    let mut all_ok = true;
    for case in 0..60 {
        let kk = 3 + case % 7;
        let gg = dna(50 + rng.below(400), &mut rng);
        let mut km = kmers_of(&gg, kk, case % 3 == 0);
        for r in reads(&gg, 5, 20.min(gg.len()), 0.05, 2, &mut rng) {
            km.extend(kmers_of(&r, kk, false));
        }
        let km = set_of(km);
        let mut sub = T { ok: 0, total: 0, failed: Vec::new() };
        against_oracle_v(&mut sub, &km, kk, &format!("k={kk}"), false);
        if !sub.failed.is_empty() {
            all_ok = false;
            println!("        random case {case} (k = {kk}): {:?}", sub.failed);
        }
    }
    t.check(all_ok, "60 random small graphs (k = 3..9, cycles, both strands) match the oracle");

    println!("== edge cases");
    let homo = kmers_of(b"AAAAAAAAAA", 5, false);
    let g = json(&unitig_graph(&homo, 5, 10).unwrap());
    t.check(
        u(&g["summary"]["nodes"]) == 1 && node(&g, 0)["circular"] == true && u(&node(&g, 0)["length"]) == 1,
        "a homopolymer: one k-mer looping on itself, a circle of 1 base",
    );
    t.eq(unitig_sequences(&homo, 5).unwrap(), vec!["AAAAA".to_string()], "the homopolymer spells AAAAA, as assemble_debruijn");
    let mut shuffled = union3.clone();
    shuffled.extend_from_slice(&union3[..100]);
    shuffled.reverse();
    t.eq(
        unitig_graph(&shuffled, k as u8, 2000).unwrap(),
        unitig_graph(&union3, k as u8, 2000).unwrap(),
        "unsorted, repeated k-mers give the same graph",
    );
    let e = unitig_graph(&union3, 21, 2000).unwrap_err();
    t.check(e.contains("another k"), format!("k-mers of k = 31 read as k = 21: {e}"));
    t.check(unitig_graph(&[], 31, 10).unwrap_err().contains("empty"), "no k-mers: an error");
    t.check(unitig_graph(&union3, 2, 10).is_err(), "k = 2: an error");
    t.check(unitig_graph(&union3, 64, 10).is_err(), "k = 64: an error");
    t.check(unitig_graph(&union3, 31, 0).is_err(), "max_nodes = 0: an error");

    println!("== max_nodes: the largest components first");
    let full = json(&unitig_graph(&km_reads, k as u8, usize::MAX).unwrap());
    let total = u(&full["summary"]["nodes"]);
    for max in [50usize, 700, total as usize] {
        let g = json(&unitig_graph(&km_reads, k as u8, max).unwrap());
        let s = &g["summary"];
        let shown = g["nodes"].as_array().unwrap().len() as u64;
        t.eq(shown, (max as u64).min(total), format!("max_nodes {max}: that many nodes shown"));
        t.eq(u(&s["shown_nodes"]) + u(&s["dropped_nodes"]), total, format!("max_nodes {max}: shown + dropped = all"));
        let ids: Vec<u64> = g["components"].as_array().unwrap().iter().map(|c| u(&c["id"])).collect();
        t.check(ids.iter().enumerate().all(|(i, &c)| c == i as u64), format!("max_nodes {max}: the components shown are the largest"));
        let sum: u64 = g["components"].as_array().unwrap().iter().map(|c| u(&c["shown_nodes"])).sum();
        t.eq(sum, shown, format!("max_nodes {max}: component counts add up"));
        let shown_set: std::collections::HashSet<u64> = g["nodes"].as_array().unwrap().iter().map(|n| u(&n["id"])).collect();
        t.check(
            g["edges"].as_array().unwrap().iter().all(|e| shown_set.contains(&u(&e["from"])) && shown_set.contains(&u(&e["to"]))),
            format!("max_nodes {max}: edges only between nodes shown"),
        );
        let full_edges = full["edges"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|e| shown_set.contains(&u(&e["from"])) && shown_set.contains(&u(&e["to"])))
            .count();
        t.eq(g["edges"].as_array().unwrap().len(), full_edges, format!("max_nodes {max}: every edge among them"));
        if let Some(p) = s["partial_component"].as_u64() {
            // the part shown is connected
            let members: Vec<u64> = g["nodes"].as_array().unwrap().iter().filter(|n| u(&n["component"]) == p).map(|n| u(&n["id"])).collect();
            let mut adj: HashMap<u64, Vec<u64>> = HashMap::new();
            for e in g["edges"].as_array().unwrap() {
                adj.entry(u(&e["from"])).or_default().push(u(&e["to"]));
                adj.entry(u(&e["to"])).or_default().push(u(&e["from"]));
            }
            let mut seen = std::collections::HashSet::from([members[0]]);
            let mut stack = vec![members[0]];
            while let Some(v) = stack.pop() {
                for &w in adj.get(&v).map(|x| x.as_slice()).unwrap_or(&[]) {
                    if seen.insert(w) {
                        stack.push(w);
                    }
                }
            }
            t.eq(seen.len(), members.len(), format!("max_nodes {max}: the partial component is connected"));
            t.check(g["text"].as_str().unwrap().contains("only in part"), format!("max_nodes {max}: the text says what was cut"));
        }
        if max as u64 >= total {
            t.check(s["partial_component"].is_null() && u(&s["dropped_nodes"]) == 0, "max_nodes = all: nothing dropped");
        }
        layout_checks(t, &g, &format!("max_nodes {max}"));
    }

    println!("== the overlap graph (screen 7)");
    // as make_examples: a 20 kb plasmid in 4 pieces that overlap by 54
    let pl = dna(20_000, &mut rng);
    let cuts = [0usize, 4_800, 9_900, 15_200, 20_000];
    let pl2: Vec<u8> = [&pl[..], &pl[..]].concat();
    let mut pieces: Vec<String> =
        (0..4).map(|i| String::from_utf8([&pl[cuts[i]..cuts[i + 1]], &pl2[cuts[i + 1]..cuts[i + 1] + 54]].concat()).unwrap()).collect();
    pieces.swap(0, 2);
    pieces.swap(1, 3);
    let g = scaffolds_match(t, &pieces, 54, "plasmid pieces");
    save("overlap_plasmid", &serde_json::to_string(&g).unwrap());
    t.eq(u(&g["summary"]["edges"]), 4, "plasmid pieces: 4 overlaps");
    t.check(g["edges"].as_array().unwrap().iter().all(|e| e["used"] == true && e["scaffold_rule"] == true && u(&e["overlap"]) == 54), "plasmid pieces: all 54 bases, all used");
    t.check(u(&g["summary"]["circular_components"]) == 1 && g["components"][0]["circular"] == true, "plasmid pieces: one circular component");
    t.check(g["scaffolds"].as_array().unwrap().len() == 1 && u(&g["scaffolds"][0]["length"]) == 20_000 && g["scaffolds"][0]["circular"] == true, "plasmid pieces: one circular scaffold of 20,000");
    t.eq(u(&g["components"][0]["length"]), 20_000, "plasmid pieces: component length = its scaffold");
    let all_pts: Vec<Value> = g["nodes"].as_array().unwrap().iter().flat_map(|n| n["path"].as_array().unwrap().clone()).collect();
    t.check(on_a_circle(&all_pts), "plasmid pieces: drawn as one ring");
    layout_checks(t, &g, "plasmid pieces");
    overlaps_brute(t, &pieces, 54, &g, "plasmid pieces");
    let g30 = json(&overlap_graph(&pieces, 30).unwrap());
    t.check(
        u(&g30["summary"]["edges"]) == 4 && u(&g30["summary"]["scaffold_rule_edges"]) == 0 && u(&g30["summary"]["used_edges"]) == 0,
        "plasmid pieces, min_overlap 30: the 54-base overlaps found, none of exactly 30",
    );
    t.check(g30["text"].as_str().unwrap().contains("most common of them is 54 bases"), "plasmid pieces, min_overlap 30: the text suggests 54");
    overlaps_brute(t, &pieces, 30, &g30, "plasmid pieces, 30");
    let lin = dna(24_000, &mut rng);
    let cuts = [0usize, 5_000, 10_400, 15_100, 19_800, 24_000];
    let mut pieces: Vec<String> =
        (0..5).map(|i| String::from_utf8(lin[cuts[i]..(cuts[i + 1] + 30).min(24_000)].to_vec()).unwrap()).collect();
    pieces.swap(0, 4);
    pieces.swap(1, 2);
    let g = scaffolds_match(t, &pieces, 30, "chromosome pieces");
    save("overlap_chromosome", &serde_json::to_string(&g).unwrap());
    t.eq(u(&g["summary"]["edges"]), 4, "chromosome pieces: 4 overlaps");
    t.check(g["scaffolds"].as_array().unwrap().len() == 1 && u(&g["scaffolds"][0]["length"]) == 24_000 && g["scaffolds"][0]["circular"] == false, "chromosome pieces: one linear scaffold of 24,000");
    t.check(u(&g["summary"]["circular_components"]) == 0, "chromosome pieces: not circular");
    layout_checks(t, &g, "chromosome pieces");
    overlaps_brute(t, &pieces, 30, &g, "chromosome pieces");

    // the unitigs as contigs: the overlap graph is the unitig graph
    let useq = unitig_sequences(&km_reads, k as u8).unwrap();
    let og = json(&overlap_graph_with_limit(&useq, k - 1, usize::MAX).unwrap());
    let mut oe: Vec<(u64, u64, u64)> =
        og["edges"].as_array().unwrap().iter().map(|e| (u(&e["from"]), u(&e["to"]), u(&e["overlap"]))).collect();
    oe.sort_unstable();
    let mut ue: Vec<(u64, u64, u64)> =
        full["edges"].as_array().unwrap().iter().map(|e| (u(&e["from"]), u(&e["to"]), u(&e["overlap"]))).collect();
    ue.sort_unstable();
    t.check(oe == ue, format!("the overlap graph of the unitigs (min_overlap k-1) is the unitig graph ({} edges)", ue.len()));
    t.eq(u(&og["summary"]["components"]), u(&full["summary"]["components"]), "... with the same components");
    scaffolds_match(t, &useq, k - 1, "the unitigs as contigs");
    // shortest first: a partial component must still start from its longest contig
    let rev: Vec<String> = useq.iter().rev().cloned().collect();
    let gp = json(&overlap_graph_with_limit(&rev, k - 1, 40).unwrap());
    let p = u(&gp["summary"]["partial_component"]);
    let gfull = json(&overlap_graph_with_limit(&rev, k - 1, usize::MAX).unwrap());
    let longest = gfull["nodes"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|n| u(&n["component"]) == p)
        .max_by(|a, b| u(&a["length"]).cmp(&u(&b["length"])).then(u(&b["id"]).cmp(&u(&a["id"]))))
        .map(|n| u(&n["id"]))
        .unwrap();
    t.check(
        gp["nodes"].as_array().unwrap().iter().any(|n| u(&n["id"]) == longest),
        "overlap graph, max_nodes 40: the partial component is drawn around its longest contig",
    );
    t.eq(gp["nodes"].as_array().unwrap().len(), 40, "overlap graph, max_nodes 40: 40 contigs drawn");

    // the scaffolder's quirks on random contig sets
    let mut all_ok = true;
    for case in 0..80 {
        let m = 3 + rng.below(12);
        let mut seqs: Vec<String> = Vec::new();
        for _ in 0..(1 + rng.below(3)) {
            let len = 60 + rng.below(300);
            let gg = dna(len, &mut rng);
            let circ = rng.chance(0.4);
            let mut pos = 0;
            let total = if circ { len } else { len - m };
            let gg2: Vec<u8> = [&gg[..], &gg[..]].concat();
            while pos < total {
                let end = (pos + 10 + rng.below(80)).min(total);
                seqs.push(String::from_utf8(gg2[pos..end + m].to_vec()).unwrap());
                pos = end;
            }
        }
        for _ in 0..rng.below(4) {
            seqs.push(String::from_utf8(dna(rng.below(2 * m + 5) + 1, &mut rng)).unwrap());
        }
        if rng.chance(0.3) && !seqs.is_empty() {
            let d = seqs[rng.below(seqs.len())].clone();
            seqs.push(d);
        }
        if rng.chance(0.3) {
            seqs.push(String::from_utf8(vec![b'A'; m + rng.below(5)]).unwrap());
        }
        for i in (1..seqs.len()).rev() {
            let j = rng.below(i + 1);
            seqs.swap(i, j);
        }
        let mut sub = T { ok: 0, total: 0, failed: Vec::new() };
        let g = scaffolds_match(&mut sub, &seqs, m, &format!("random contigs {case}"));
        overlaps_brute(&mut sub, &seqs, m, &g, &format!("random contigs {case}"));
        all_ok &= sub.failed.is_empty();
    }
    t.check(all_ok, "80 random contig sets: scaffolds equal to fast_hash_scaffold's, overlaps to brute force");
    t.check(overlap_graph(&[], 30).is_err(), "overlap graph: no contigs is an error");
    t.check(overlap_graph(&pieces, 0).is_err(), "overlap graph: min_overlap 0 is an error");

    println!("== layout: determinism and speed");
    for (label, km) in [("3 SNPs", &union3), ("reads", &km_reads)] {
        let a = unitig_graph(km, k as u8, 3000).unwrap();
        let b = unitig_graph(km, k as u8, 3000).unwrap();
        t.check(a == b, format!("{label}: the same JSON twice"));
    }
    let a = overlap_graph(&pieces, 30).unwrap();
    t.check(a == overlap_graph(&pieces, 30).unwrap(), "overlap graph: the same JSON twice");
    // a graph of more than 5,000 unitigs
    let g50 = dna(40_000, &mut rng);
    let rs = reads(&g50, 2_500, 150, 0.01, 1, &mut rng);
    let big = kmers_of_reads(&rs, k);
    let (s, dt) = best_of_3(|| unitig_graph(&big, k as u8, 5_000).unwrap());
    let g = json(&s);
    save("big5000", &s);
    let nodes_all = u(&g["summary"]["nodes"]);
    t.eq(u(&g["summary"]["shown_nodes"]), 5_000, format!("big graph ({nodes_all} unitigs): 5,000 shown"));
    t.check(dt < LIMIT, format!("big graph: 5,000 nodes laid out in {dt:.2} s (best of 3)"));
    println!("  {} k-mers, {nodes_all} unitigs; 5,000 shown in {dt:.3} s (best of 3)", big.len());
    layout_checks(t, &g, "big graph (5,000 shown)");
    let t0 = Instant::now();
    let s2 = unitig_graph(&big, k as u8, 5_000).unwrap();
    println!("  again: {:.3} s, identical: {}", t0.elapsed().as_secs_f64(), s2 == s);
    t.check(s2 == s, "big graph: the same JSON twice");
    let t0 = Instant::now();
    let gall = json(&unitig_graph(&big, k as u8, usize::MAX).unwrap());
    println!("  all {nodes_all} unitigs: {:.3} s", t0.elapsed().as_secs_f64());
    layout_checks(t, &gall, "big graph (all)");
}

fn examples(t: &mut T, dir: &str) {
    println!("== the UI examples in {dir}");
    let k = 31u8;
    let bac = read_bin(&format!("{dir}/toy_bacterium_K31.bin"));
    let pat = read_bin(&format!("{dir}/toy_patient_K31.bin"));
    let rds = read_bin(&format!("{dir}/toy_bacterium_reads_K31.bin"));
    if let (Some(bac), Some(pat)) = (&bac, &pat) {
        let g = json(&unitig_graph(bac, k, 2000).unwrap());
        t.check(u(&g["summary"]["nodes"]) == 1 && u(&node(&g, 0)["length"]) == 30_000, "toy bacterium: one unitig of 30,000 bases");
        let g = json(&unitig_graph(pat, k, 2000).unwrap());
        t.check(u(&g["summary"]["nodes"]) == 1 && u(&node(&g, 0)["length"]) == 30_000, "toy patient alone: one unitig of 30,000 bases");
        let both = set_of([bac.clone(), pat.clone()].concat());
        let s = unitig_graph(&both, k, 2000).unwrap();
        save("example_bacterium_patient", &s);
        let g = json(&s);
        t.eq(
            (u(&g["summary"]["nodes"]), u(&g["summary"]["bubbles"]), u(&g["summary"]["bubble_branches"])),
            (10, 3, 6),
            "toy bacterium + patient: 10 unitigs, 3 bubbles of 2 branches",
        );
        t.eq(u(&g["summary"]["repeats"]), 0, "toy bacterium + patient: no repeat");
        t.check(
            g["nodes"].as_array().unwrap().iter().filter(|n| n["kind"] == "bubble").all(|n| u(&n["length"]) == 61),
            "toy bacterium + patient: every branch 61 bases (2k-1)",
        );
        against_oracle(t, &both, 31, "toy bacterium + patient");
        layout_checks(t, &g, "toy bacterium + patient");
    } else {
        println!("  (no toy_bacterium_K31.bin / toy_patient_K31.bin: skipped)");
    }
    if let Some(rds) = &rds {
        let g = against_oracle(t, rds, 31, "toy reads");
        for max in [2000usize, 5000] {
            let (s, dt) = best_of_3(|| unitig_graph(rds, k, max).unwrap());
            save(&format!("example_reads_{max}"), &s);
            let gg = json(&s);
            println!("  toy reads, max_nodes {max}: {:.3} s (best of 3); {}", dt, gg["text"].as_str().unwrap());
            t.check(dt < LIMIT, format!("toy reads, max_nodes {max}: {dt:.2} s (best of 3)"));
            layout_checks(t, &gg, &format!("toy reads, max_nodes {max}"));
        }
        let _ = g;
        // the overlap graph of the assembler's contigs, as screen 7 would get them
        let contigs = debruijn::montar_todos_contigs(rds, k);
        let (s, dt) = best_of_3(|| overlap_graph(&contigs, 30).unwrap());
        let g = json(&s);
        println!("  overlap graph of the {} contigs (min_overlap 30, 2,000 drawn): {:.3} s (best of 3); {}", contigs.len(), dt, g["text"].as_str().unwrap());
        t.check(dt < LIMIT, format!("overlap graph of the read contigs: {dt:.2} s (best of 3)"));
        layout_checks(t, &g, "overlap graph of the read contigs");
        scaffolds_match(t, &contigs, 30, "the read contigs");
    } else {
        println!("  (no toy_bacterium_reads_K31.bin: skipped)");
    }
    for (file, m, len, circ) in [("plasmid_contigs.fasta", 54usize, 20_000u64, true), ("chromosome_contigs.fasta", 30, 24_000, false)] {
        let seqs = read_fasta(&format!("{dir}/mmgen_examples/{file}"));
        if seqs.is_empty() {
            println!("  (no {file}: skipped)");
            continue;
        }
        let g = scaffolds_match(t, &seqs, m, file);
        save(&format!("example_{}", file.trim_end_matches(".fasta")), &serde_json::to_string(&g).unwrap());
        t.check(
            g["scaffolds"].as_array().unwrap().len() == 1
                && u(&g["scaffolds"][0]["length"]) == len
                && g["scaffolds"][0]["circular"].as_bool() == Some(circ),
            format!("{file}: one scaffold of {len} bases, circular = {circ}"),
        );
        t.check(u(&g["summary"]["components"]) == 1 && g["components"][0]["circular"].as_bool() == Some(circ), format!("{file}: one component, circular = {circ}"));
        overlaps_brute(t, &seqs, m, &g, file);
        layout_checks(t, &g, file);
        println!("  {file}: {}", g["text"].as_str().unwrap());
    }
}

fn main() {
    let mut t = T { ok: 0, total: 0, failed: Vec::new() };
    synthetic(&mut t);
    if let Some(dir) = std::env::args().nth(1) {
        examples(&mut t, &dir);
    }
    println!("\ngrafo_montagem tests: {} passed, {} failed", t.ok, t.total - t.ok);
    if t.ok != t.total {
        for f in &t.failed {
            println!("  failed: {f}");
        }
        std::process::exit(1);
    }
}
