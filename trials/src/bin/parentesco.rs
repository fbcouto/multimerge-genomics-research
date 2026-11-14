#![allow(dead_code)]
// Tests of src/parentesco.rs: relatedness, trio analysis, runs of
// homozygosity and genetic sex.
//
// Everything is checked against something this file computes on its own,
// from what the simulation WROTE, not by reading the files back:
//  - pedigrees are simulated (founders from an allele-frequency spectrum,
//    children by Mendelian transmission, an MZ twin, half-siblings, a
//    grandchild, first cousins), with genotype errors and missing calls;
//  - a brute-force join (hash maps over the written records, with the absent
//    / BED / different-ALT rules applied by hand) gives every IBS count and
//    every trio genotype table exactly;
//  - the kinship is compared with the pedigree kinship coefficient;
//  - the paternity index is compared with the textbook "obligate paternal
//    allele" formulas;
//  - runs of homozygosity are compared with an independent scan and with the
//    true autozygous segments of a child of first cousins (linked simulation,
//    with crossovers);
//  - genetic sex is compared with the counts of the generator.
//
// Usage:
//   cargo run --release --bin parentesco
//   cargo run --release --bin parentesco -- real relatedness "a.vcf.gz|SAMPLE|callable.bed|label" ...
//   cargo run --release --bin parentesco -- real trio CHILD FATHER [MOTHER|-] [--af file.vcf.gz FIELD] [--denovo out.vcf]
//   cargo run --release --bin parentesco -- real roh PERSON [out.bed]
//   cargo run --release --bin parentesco -- real sex PERSON [GRCh38|GRCh37]
// (a PERSON is "path|sample|bed|label"; empty parts are left out)
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
#[path = "../../../src/nomes.rs"]
#[allow(dead_code)]
mod nomes;
#[path = "../../../src/parentesco.rs"]
mod parentesco;

use parentesco::*;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt::Debug;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

// ---------------------------------------------------------------- random numbers

/// splitmix64: every bit is usable (the LCG of the other trials repeats its
/// low bits every few thousand draws).
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn f(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn chance(&mut self, p: f64) -> bool {
        self.f() < p
    }
    fn poisson(&mut self, l: f64) -> u64 {
        let (lim, mut k, mut p) = ((-l).exp(), 0u64, 1.0);
        loop {
            p *= self.f();
            if p <= lim {
                return k;
            }
            k += 1;
        }
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[self.below(4) as usize]
    }
    fn outra(&mut self, b: u8) -> u8 {
        loop {
            let x = self.base();
            if x != b {
                return x;
            }
        }
    }
    /// An allele frequency: half common, half rare-ish.
    fn freq(&mut self) -> f64 {
        if self.chance(0.5) {
            0.02 + 0.96 * self.f()
        } else {
            0.003 + 0.25 * self.f() * self.f()
        }
    }
}

// ---------------------------------------------------------------- checks

struct T {
    n: usize,
    falhas: Vec<String>,
    dir: PathBuf,
}

impl T {
    fn check(&mut self, ok: bool, what: impl AsRef<str>) {
        self.n += 1;
        if !ok {
            eprintln!("FAILED: {}", what.as_ref());
            self.falhas.push(what.as_ref().to_string());
        }
    }
    fn eq<X: PartialEq + Debug>(&mut self, got: X, want: X, what: &str) {
        let ok = got == want;
        self.check(ok, format!("{what}: got {got:?}, expected {want:?}"));
    }
    fn near(&mut self, got: f64, want: f64, tol: f64, what: &str) {
        self.check((got - want).abs() <= tol, format!("{what}: got {got}, expected {want} +- {tol}"));
    }
    fn err<X: Debug>(&mut self, r: Result<X, String>, needles: &[&str], what: &str) {
        match r {
            Ok(x) => self.check(false, format!("{what}: should have failed, returned {x:?}")),
            Err(e) => {
                let falta: Vec<&&str> = needles.iter().filter(|n| !e.contains(**n)).collect();
                self.check(falta.is_empty(), format!("{what}: the error lacks {falta:?}: {e}"));
            }
        }
    }
    fn p(&self, n: &str) -> String {
        self.dir.join(n).to_string_lossy().to_string()
    }
}

fn grava(path: &str, linhas: &[String]) {
    let mut w = externo::Escritor::cria(path).unwrap_or_else(|e| panic!("{e}"));
    for l in linhas {
        w.write_all(l.as_bytes()).unwrap();
        w.write_all(b"\n").unwrap();
    }
    w.termina().unwrap();
}

fn le_linhas(path: &str) -> Vec<String> {
    let mut r = externo::abre_leitura(path).unwrap();
    let mut s = String::new();
    std::io::Read::read_to_string(&mut r, &mut s).unwrap();
    s.lines().map(str::to_string).collect()
}

fn cabecalho(contigs: &[(&str, u64)], amostras: &[&str]) -> Vec<String> {
    let mut v = vec![
        "##fileformat=VCFv4.2".to_string(),
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
    ];
    for (c, l) in contigs {
        v.push(format!("##contig=<ID={c},length={l}>"));
    }
    let mut h = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string();
    if !amostras.is_empty() {
        h.push_str("\tFORMAT");
        for a in amostras {
            h.push('\t');
            h.push_str(a);
        }
    }
    v.push(h);
    v
}

/// The GT text of an ALT count, in one of its spellings; -1 is missing.
fn gt_txt(c: i8, rng: &mut Rng) -> &'static str {
    match c {
        0 => ["0/0", "0|0"][rng.below(2) as usize],
        1 => ["0/1", "0|1", "1|0", "1/0"][rng.below(4) as usize],
        2 => ["1/1", "1|1"][rng.below(2) as usize],
        _ => ["./.", ".", "0/.", ".|.", "./1"][rng.below(5) as usize],
    }
}

// ---------------------------------------------------------------- the brute-force model

/// A record as WRITTEN in a person's file, with the ALT count the engine
/// must read from it for this person.
#[derive(Clone, Debug)]
struct RegM {
    r: String,
    a: String,
    cod: i8,
    pass: bool,
}

#[derive(Default)]
struct PessoaM {
    /// (canonical chromosome, position) -> the records of the person's FILE there.
    regs: HashMap<(String, u64), Vec<RegM>>,
    /// canonical chromosome -> sorted, disjoint intervals (0-based, end excluded)
    bed: Option<HashMap<String, Vec<(u64, u64)>>>,
}

impl PessoaM {
    fn dentro(&self, c: &str, pos: u64) -> Option<bool> {
        self.bed.as_ref().map(|b| {
            b.get(c).is_some_and(|v| {
                let i = v.partition_point(|iv| iv.1 < pos);
                i < v.len() && v[i].0 < pos
            })
        })
    }
    fn add(&mut self, c: &str, pos: u64, r: &str, a: &str, cod: i8, pass: bool) {
        self.regs.entry((c.to_string(), pos)).or_default().push(RegM { r: r.into(), a: a.into(), cod, pass });
    }
}

fn plain(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}

fn e_snv_txt(r: &str, a: &str) -> bool {
    r.len() == 1 && a.len() == 1 && r != a && r.bytes().chain(a.bytes()).all(|b| b"ACGT".contains(&b))
}

type Chave = (String, u64, String, String);

/// Every biallelic site in the union of the files, with each person's ALT
/// count under the rules, done by hand.
fn junta_bruto(ps: &[&PessoaM], absent_is_ref: bool, pass_only: bool) -> BTreeMap<Chave, Vec<i8>> {
    let mut chaves: BTreeSet<Chave> = BTreeSet::new();
    for p in ps {
        for ((c, pos), v) in &p.regs {
            for r in v {
                if !r.a.contains(',') && plain(&r.r) && plain(&r.a) && !r.r.eq_ignore_ascii_case(&r.a) {
                    chaves.insert((c.clone(), *pos, r.r.to_ascii_uppercase(), r.a.to_ascii_uppercase()));
                }
            }
        }
    }
    let mut out = BTreeMap::new();
    for k in chaves {
        let cods = ps
            .iter()
            .map(|p| {
                let dentro = p.dentro(&k.0, k.1);
                match p.regs.get(&(k.0.clone(), k.1)) {
                    None => match dentro {
                        Some(true) => 0,
                        Some(false) => -1,
                        None => {
                            if absent_is_ref {
                                0
                            } else {
                                -1
                            }
                        }
                    },
                    Some(v) => match v.iter().find(|r| r.r.eq_ignore_ascii_case(&k.2) && r.a.eq_ignore_ascii_case(&k.3)) {
                        Some(r) if dentro != Some(false) && !(pass_only && !r.pass) => r.cod,
                        _ => -1,
                    },
                }
            })
            .collect();
        out.insert(k, cods);
    }
    out
}

/// [sites, hethet, ibs0, het1_hom2, het2_hom1, ibs2] by direct comparison.
fn par_bruto(sitios: &BTreeMap<Chave, Vec<i8>>, i: usize, j: usize) -> [u64; 6] {
    let mut c = [0u64; 6];
    for (k, v) in sitios {
        if !e_snv_txt(&k.2, &k.3) {
            continue;
        }
        let (a, b) = (v[i], v[j]);
        if a < 0 || b < 0 {
            continue;
        }
        c[0] += 1;
        if a == 1 && b == 1 {
            c[1] += 1;
        }
        if (a - b).abs() == 2 {
            c[2] += 1;
        }
        if a == 1 && b != 1 {
            c[3] += 1;
        }
        if b == 1 && a != 1 {
            c[4] += 1;
        }
        if a == b {
            c[5] += 1;
        }
    }
    c
}

fn pares_de(r: &RelatednessReport) -> HashMap<(String, String), &PairResult> {
    r.pairs.iter().map(|p| ((p.a.clone(), p.b.clone()), p)).collect()
}

// =====================================================================
//  SIMULATION A: a pedigree, unlinked sites, errors, several file shapes
// =====================================================================

struct Ped {
    nome: &'static str,
    pai: Option<usize>,
    mae: Option<usize>,
    gemeo: Option<usize>,
}

const PED_A: [Ped; 12] = [
    Ped { nome: "F", pai: None, mae: None, gemeo: None },
    Ped { nome: "M", pai: None, mae: None, gemeo: None },
    Ped { nome: "C1", pai: Some(0), mae: Some(1), gemeo: None },
    Ped { nome: "C2", pai: Some(0), mae: Some(1), gemeo: None },
    Ped { nome: "C1T", pai: None, mae: None, gemeo: Some(2) },
    Ped { nome: "M2", pai: None, mae: None, gemeo: None },
    Ped { nome: "H", pai: Some(0), mae: Some(5), gemeo: None },
    Ped { nome: "S1", pai: None, mae: None, gemeo: None },
    Ped { nome: "G", pai: Some(2), mae: Some(7), gemeo: None },
    Ped { nome: "S2", pai: None, mae: None, gemeo: None },
    Ped { nome: "K", pai: Some(3), mae: Some(9), gemeo: None },
    Ped { nome: "U", pai: None, mae: None, gemeo: None },
];

/// Kinship coefficients of the pedigree (parents come before children).
fn kinship_ped(ped: &[Ped]) -> Vec<Vec<f64>> {
    let n = ped.len();
    let mut k = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..i {
            let v = match (ped[i].gemeo, ped[i].pai, ped[i].mae) {
                (Some(t), _, _) => {
                    if j == t {
                        k[t][t]
                    } else {
                        k[t][j]
                    }
                }
                (None, Some(f), Some(m)) => 0.5 * (k[f][j] + k[m][j]),
                _ => 0.0,
            };
            k[i][j] = v;
            k[j][i] = v;
        }
        k[i][i] = match (ped[i].gemeo, ped[i].pai, ped[i].mae) {
            (Some(t), _, _) => k[t][t],
            (None, Some(f), Some(m)) => 0.5 * (1.0 + k[f][m]),
            _ => 0.5,
        };
    }
    k
}

fn classe_esperada(ped: &[Ped], phi: f64, i: usize, j: usize) -> &'static str {
    let ef = |x: usize| ped[x].gemeo.unwrap_or(x);
    let (a, b) = (ef(i), ef(j));
    let pais = |x: usize| [ped[x].pai, ped[x].mae];
    if phi > 0.354 {
        "duplicate_or_mz_twin"
    } else if phi >= 0.177 {
        if pais(a).contains(&Some(b)) || pais(b).contains(&Some(a)) {
            "parent_offspring"
        } else {
            "full_siblings"
        }
    } else if phi >= 0.0884 {
        "second_degree"
    } else if phi >= 0.0442 {
        "third_degree"
    } else {
        "unrelated"
    }
}

struct SitioA {
    c: usize,
    pos: u64,
    r: String,
    a: String,
    p: f64,
}

struct SimA {
    sitios: Vec<SitioA>,
    /// observed ALT count per person per site (-1 missing)
    obs: Vec<Vec<i8>>,
    modelo: Vec<PessoaM>,
    pessoas: Vec<Person>,
    de_novo: Vec<(String, u64)>,
}

fn simula_a(t: &T, rng: &mut Rng) -> SimA {
    const CROMS: [&str; 3] = ["1", "2", "3"];
    const N: usize = 16_000;
    let mut sitios = Vec::new();
    for c in 0..CROMS.len() {
        let mut pos = 10_000u64;
        for _ in 0..N {
            pos += 20 + rng.below(400);
            let rb = rng.base();
            let (r, a) = if rng.chance(0.85) {
                ((rb as char).to_string(), (rng.outra(rb) as char).to_string())
            } else if rng.chance(0.5) {
                let ins: String = (0..1 + rng.below(3)).map(|_| rng.base() as char).collect();
                ((rb as char).to_string(), format!("{}{ins}", rb as char))
            } else {
                let del: String = (0..1 + rng.below(3)).map(|_| rng.base() as char).collect();
                (format!("{}{del}", rb as char), (rb as char).to_string())
            };
            let p = rng.freq();
            sitios.push(SitioA { c, pos, r, a, p });
        }
    }
    // planted de novo mutations of C1: new SNVs 5 bases after 30 sites
    let mut alvos = BTreeSet::new();
    while alvos.len() < 30 {
        alvos.insert(rng.below(sitios.len() as u64) as usize);
    }
    let mut dn_idx = Vec::new();
    for &i in &alvos {
        let rb = rng.base();
        let s = SitioA {
            c: sitios[i].c,
            pos: sitios[i].pos + 5,
            r: (rb as char).to_string(),
            a: (rng.outra(rb) as char).to_string(),
            p: 0.0,
        };
        sitios.push(s);
        dn_idx.push(sitios.len() - 1);
    }
    let mut ordem: Vec<usize> = (0..sitios.len()).collect();
    ordem.sort_by_key(|&i| (sitios[i].c, sitios[i].pos));
    let de_novo_pos: BTreeSet<(usize, u64)> = dn_idx.iter().map(|&i| (sitios[i].c, sitios[i].pos)).collect();
    let mut s2: Vec<SitioA> = Vec::with_capacity(sitios.len());
    let mut velho: Vec<Option<SitioA>> = sitios.into_iter().map(Some).collect();
    for i in ordem {
        s2.push(velho[i].take().unwrap());
    }
    let sitios = s2;
    let ns = sitios.len();

    // haplotypes, by Mendelian transmission (independent sites)
    let np = PED_A.len();
    let mut hap: Vec<[Vec<u8>; 2]> = Vec::with_capacity(np);
    for (i, pd) in PED_A.iter().enumerate() {
        let h = match (pd.gemeo, pd.pai, pd.mae) {
            (Some(t), _, _) => hap[t].clone(),
            (None, Some(f), Some(m)) => {
                let mut h = [vec![0u8; ns], vec![0u8; ns]];
                for s in 0..ns {
                    h[0][s] = hap[f][rng.below(2) as usize][s];
                    h[1][s] = hap[m][rng.below(2) as usize][s];
                }
                h
            }
            _ => {
                let mut h = [vec![0u8; ns], vec![0u8; ns]];
                for s in 0..ns {
                    h[0][s] = u8::from(rng.chance(sitios[s].p));
                    h[1][s] = u8::from(rng.chance(sitios[s].p));
                }
                h
            }
        };
        hap.push(h);
        if i == 4 {
            // after the twin: the de novo mutations are C1's only (G inherits some)
            for s in 0..ns {
                if de_novo_pos.contains(&(sitios[s].c, sitios[s].pos)) {
                    hap[2][0][s] = 1;
                }
            }
        }
    }
    // observed: errors (0.2%: one allele read wrong) and missing calls (1%)
    let mut obs = vec![vec![0i8; ns]; np];
    for i in 0..np {
        for s in 0..ns {
            let v = (hap[i][0][s] + hap[i][1][s]) as i8;
            obs[i][s] = if rng.chance(0.01) {
                -1
            } else if rng.chance(0.002) {
                match v {
                    1 => 2 * rng.below(2) as i8,
                    _ => 1,
                }
            } else {
                v
            };
        }
    }

    let mut modelo: Vec<PessoaM> = (0..np).map(|_| PessoaM::default()).collect();
    let contigs: Vec<(String, u64)> = CROMS.iter().map(|c| (c.to_string(), 200_000_000)).collect();

    // --- the family file: C1, F, C1T, M; `chr` names; every site; BGZF
    let colunas = [2usize, 0, 4, 1];
    let cont_chr: Vec<(String, u64)> = contigs.iter().map(|(c, l)| (format!("chr{c}"), *l)).collect();
    let cc: Vec<(&str, u64)> = cont_chr.iter().map(|(c, l)| (c.as_str(), *l)).collect();
    let mut linhas = cabecalho(&cc, &colunas.map(|i| PED_A[i].nome));
    for (si, s) in sitios.iter().enumerate() {
        let multi = rng.chance(0.01);
        let alt = if multi { format!("{},{}", s.a, if s.a == "T" { "G" } else { "T" }) } else { s.a.clone() };
        let mut l = format!("chr{}\t{}\t.\t{}\t{}\t50\tPASS\t.\tGT", CROMS[s.c], s.pos, s.r, alt);
        for &i in &colunas {
            l.push('\t');
            l.push_str(gt_txt(obs[i][si], rng));
            modelo[i].add(CROMS[s.c], s.pos, &s.r, &alt, obs[i][si], true);
        }
        linhas.push(l);
    }
    // other sequences after chr3: skipped, whatever their order
    linhas.push("chrUn_KI270742v1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/1\t0/1".into());
    linhas.push("chrM\t50\t.\tA\tG\t50\tPASS\t.\tGT\t1\t1\t1\t1".into());
    linhas.push("chrM\t60\t.\tA\tG\t50\tPASS\t.\tGT\t1\t1\t1\t1".into());
    let fam = t.p("famA.vcf.gz");
    grava(&fam, &linhas);

    // --- one file per person for the others, without `chr`
    let bed_de = |rng: &mut Rng| -> HashMap<String, Vec<(u64, u64)>> {
        let mut b = HashMap::new();
        for c in CROMS {
            let mut v = Vec::new();
            let mut p = rng.below(3000);
            while p < 10_000 + (N as u64) * 420 {
                let ini = p;
                let fim = ini + 5_000 + rng.below(55_000);
                v.push((ini, fim));
                p = fim + 500 + rng.below(8_000);
            }
            b.insert(c.to_string(), v);
        }
        b
    };
    let mut pessoas: Vec<Person> = Vec::new();
    for i in 0..np {
        if colunas.contains(&i) {
            pessoas.push(Person::new(&fam).sample(PED_A[i].nome));
            continue;
        }
        let nome = PED_A[i].nome;
        let ext = match nome {
            "G" => "vcf.gz",
            "K" => "vcf.lz4",
            _ => "vcf",
        };
        let caminho = t.p(&format!("ind_{nome}.{ext}"));
        let bed = matches!(nome, "G" | "K").then(|| bed_de(rng));
        if let Some(b) = &bed {
            let mut bl = Vec::new();
            for c in CROMS {
                for (a, z) in &b[c] {
                    bl.push(format!("chr{c}\t{a}\t{z}"));
                }
            }
            grava(&t.p(&format!("ind_{nome}.bed")), &bl);
        }
        modelo[i].bed = bed.clone();
        let cc: Vec<(&str, u64)> = contigs.iter().map(|(c, l)| (c.as_str(), *l)).collect();
        let mut linhas = cabecalho(&cc, &[nome]);
        for (si, s) in sitios.iter().enumerate() {
            let o = obs[i][si];
            let dentro = bed.as_ref().map(|b| {
                let v = &b[CROMS[s.c]];
                let k = v.partition_point(|iv| iv.1 < s.pos);
                k < v.len() && v[k].0 < s.pos
            });
            let mut escreve = matches!(o, 1 | 2 | -1) || rng.chance(0.02);
            let mut cod = o;
            if dentro == Some(false) {
                // outside the callable regions: junk the engine must ignore
                escreve = rng.chance(0.2);
                cod = 1 + rng.below(2) as i8;
            }
            if !escreve {
                continue;
            }
            let mut alt = s.a.clone();
            let mut gt = gt_txt(cod, rng).to_string();
            let mut filtro = "PASS";
            if nome == "H" && cod > 0 && s.r.len() == 1 && s.a.len() == 1 && rng.chance(0.02) {
                // another ALT allele at the same position
                let b = s.a.as_bytes()[0];
                let mut x = rng.outra(b);
                while x == s.r.as_bytes()[0] {
                    x = rng.outra(b);
                }
                alt = (x as char).to_string();
            }
            if nome == "U" && cod == 1 && rng.chance(0.02) {
                alt = format!("{},{}", s.a, if s.a == "T" { "G" } else { "T" });
            }
            if nome == "S2" && rng.chance(0.03) {
                filtro = "LowQual";
            }
            if nome == "M2" && cod == 2 && rng.chance(0.01) {
                gt = "1".to_string(); // a haploid call on an autosome: missing
                cod = -1;
            }
            linhas.push(format!("{}\t{}\t.\t{}\t{}\t50\t{filtro}\t.\tGT\t{gt}", CROMS[s.c], s.pos, s.r, alt));
            modelo[i].add(CROMS[s.c], s.pos, &s.r, &alt, cod, filtro == "PASS");
            if nome == "M2" && rng.chance(0.003) {
                // a symbolic allele at a position nobody else has
                linhas.push(format!("{}\t{}\t.\tN\t<DEL>\t50\tPASS\t.\tGT\t0/1", CROMS[s.c], s.pos + 1));
                modelo[i].add(CROMS[s.c], s.pos + 1, "N", "<DEL>", 1, true);
            }
        }
        if nome == "H" {
            linhas.push("MT\t10\t.\tA\tG\t50\tPASS\t.\tGT\t1".into());
        }
        grava(&caminho, &linhas);
        let mut p = Person::new(&caminho);
        if bed.is_some() {
            p = p.bed(&t.p(&format!("ind_{nome}.bed")));
        }
        if nome == "U" {
            p = p.label("U");
        }
        pessoas.push(p);
    }
    let de_novo = de_novo_pos.iter().map(|&(c, p)| (CROMS[c].to_string(), p)).collect();
    SimA { sitios, obs, modelo, pessoas, de_novo }
}

fn testa_parentesco(t: &mut T, sim: &SimA) {
    let np = PED_A.len();
    let phi = kinship_ped(&PED_A);
    let refs: Vec<&PessoaM> = sim.modelo.iter().collect();
    for (rotulo, opts) in [
        ("default", JoinOptions::default()),
        ("absent_is_ref=false", JoinOptions { absent_is_ref: false, pass_only: false }),
        ("pass_only", JoinOptions { absent_is_ref: true, pass_only: true }),
    ] {
        let t0 = Instant::now();
        let r = relatedness(&sim.pessoas, &opts).unwrap_or_else(|e| panic!("{rotulo}: {e}"));
        let dt = t0.elapsed();
        let bruto = junta_bruto(&refs, opts.absent_is_ref, opts.pass_only);
        let snv = bruto.keys().filter(|k| e_snv_txt(&k.2, &k.3)).count() as u64;
        t.eq(r.snv_sites, snv, &format!("{rotulo}: SNV sites"));
        t.eq(r.scan.sites, bruto.len() as u64, &format!("{rotulo}: biallelic sites"));
        let nomes: Vec<String> = r.persons.iter().map(|p| p.label.clone()).collect();
        t.eq(nomes, PED_A.iter().map(|p| p.nome.to_string()).collect::<Vec<_>>(), &format!("{rotulo}: labels"));
        // per person
        for (i, p) in r.persons.iter().enumerate() {
            let mut c = [0u64; 4];
            for (k, v) in &bruto {
                if e_snv_txt(&k.2, &k.3) {
                    c[if v[i] < 0 { 3 } else { v[i] as usize }] += 1;
                }
            }
            t.eq([p.hom_ref, p.het, p.hom_alt, p.missing], c, &format!("{rotulo}: person {}", p.label));
        }
        let pares = pares_de(&r);
        t.eq(r.pairs.len(), np * (np - 1) / 2, &format!("{rotulo}: number of pairs"));
        for i in 0..np {
            for j in i + 1..np {
                let (a, b) = (PED_A[i].nome, PED_A[j].nome);
                let p = pares[&(a.to_string(), b.to_string())];
                let c = par_bruto(&bruto, i, j);
                t.eq(
                    [p.sites, p.hethet, p.ibs0, p.het1_hom2, p.het2_hom1, p.ibs2],
                    c,
                    &format!("{rotulo}: counts of {a}-{b}"),
                );
                t.eq(p.ibs1, c[3] + c[4], &format!("{rotulo}: IBS1 of {a}-{b}"));
                if rotulo == "default" {
                    let k = p.kinship.unwrap_or(f64::NAN);
                    let menor = (c[1] + c[3]).min(c[1] + c[4]) as f64;
                    let entre = 0.5 - (4.0 * c[2] as f64 + (c[3] + c[4]) as f64) / (4.0 * menor);
                    t.near(k, entre, 1e-12, &format!("KING between-family (plink2) formula {a}-{b}"));
                    let dentro = (c[1] as f64 - 2.0 * c[2] as f64) / (c[3] + c[4] + 2 * c[1]) as f64;
                    t.near(p.kinship_within.unwrap_or(f64::NAN), dentro, 1e-12, &format!("KING within-family formula {a}-{b}"));
                    t.near(p.kinship_within.unwrap_or(f64::NAN), phi[i][j], 0.02, &format!("within-family kinship {a}-{b} against the pedigree"));
                    t.near(k, phi[i][j], 0.02, &format!("kinship {a}-{b} against the pedigree"));
                    t.eq(p.relationship, classe_esperada(&PED_A, phi[i][j], i, j), &format!("relationship {a}-{b} (kinship {k:.4}, IBS0/HETHET {:?})", p.ibs0_ratio));
                }
            }
        }
        if rotulo == "default" {
            let fs = pares[&("C1".to_string(), "C2".to_string())];
            let po = pares[&("F".to_string(), "C1".to_string())];
            let un = pares[&("F".to_string(), "M".to_string())];
            println!(
                "  simulation A: {} persons, {} SNV sites, {:.2?}; IBS0/HETHET: parent-child {:.4}, full siblings {:.4}, unrelated {:.4}",
                np,
                snv,
                dt,
                po.ibs0_ratio.unwrap(),
                fs.ibs0_ratio.unwrap(),
                un.ibs0_ratio.unwrap()
            );
            // the counts that tell the absent rules apart really changed
            t.check(r.scan.multiallelic > 0 && r.scan.not_plain > 0 && r.scan.other_sequences >= 4, format!("scan counts {:?}", r.scan));
            t.check(r.resumo().contains("parent and child") && r.json().starts_with("{\"persons\":["), "text and JSON");
            t.eq(r.autosomes, 3, "autosomes with data");
            t.check(r.resumo().contains("WARNING: the data cover only 3 autosome(s)"), "few autosomes: the text warns");
        }
    }
    // a person against herself (same file, same sample): duplicate
    let r = relatedness(&[sim.pessoas[0].clone(), sim.pessoas[0].clone(), sim.pessoas[4].clone()], &JoinOptions::default()).unwrap();
    t.eq(r.pairs[0].relationship, "duplicate_or_mz_twin", "a person twice");
    t.eq(r.pairs[0].ibs0 + r.pairs[0].ibs1, 0, "a person twice: every site IBS2");
    t.eq(r.persons[1].label.as_str(), "F_2", "a repeated label gets a suffix");
    t.eq(r.pairs[1].relationship, "parent_offspring", "the twin of the son is genetically a son");
    // the class rule at its edges
    for (k, razao, n, esp) in [
        (Some(0.5), Some(0.0), 10_000, "duplicate_or_mz_twin"),
        (Some(0.27), Some(0.001), 10_000, "parent_offspring"),
        (Some(0.33), Some(0.001), 10_000, "full_siblings"), // IBD2 somewhere: not a parent
        (Some(0.25), Some(0.07), 10_000, "full_siblings"),
        (Some(0.25), Some(0.0299), 10_000, "parent_offspring"),
        (Some(0.12), Some(0.2), 10_000, "second_degree"),
        (Some(0.06), Some(0.35), 10_000, "third_degree"),
        (Some(0.01), Some(0.5), 10_000, "unrelated"),
        (Some(-0.15), Some(0.9), 10_000, "unrelated"),
        (Some(0.25), Some(0.0), 99, "insufficient_data"),
        (None, None, 10_000, "insufficient_data"),
    ] {
        t.eq(classify_relationship(k, razao, n), esp, &format!("class of kinship {k:?}, IBS0/HETHET {razao:?}, {n} sites"));
    }
}

fn mendel_bruto(sitios: &BTreeMap<Chave, Vec<i8>>, c: usize, f: usize, m: Option<usize>) -> (Vec<u64>, Vec<(String, u64)>, u64, u64, u64, u64) {
    let mut tab = vec![0u64; if m.is_some() { 27 } else { 9 }];
    let mut dn = Vec::new();
    let (mut inf_f, mut exc_f, mut inf_m, mut exc_m) = (0, 0, 0, 0);
    for (k, v) in sitios {
        let (cc, ff) = (v[c], v[f]);
        let mm = m.map_or(-1, |m| v[m]);
        match m {
            Some(_) if cc >= 0 && ff >= 0 && mm >= 0 => tab[((ff * 3 + mm) * 3 + cc) as usize] += 1,
            None if cc >= 0 && ff >= 0 => tab[(ff * 3 + cc) as usize] += 1,
            _ => {}
        }
        if m.is_some() && cc == 1 && ff == 0 && mm == 0 {
            dn.push((k.0.clone(), k.1));
        }
        if e_snv_txt(&k.2, &k.3) && (cc == 0 || cc == 2) {
            if ff >= 0 {
                inf_f += 1;
                exc_f += u64::from(ff == 2 - cc);
            }
            if mm >= 0 {
                inf_m += 1;
                exc_m += u64::from(mm == 2 - cc);
            }
        }
    }
    (tab, dn, inf_f, exc_f, inf_m, exc_m)
}

fn testa_trio_a(t: &mut T, sim: &SimA) {
    let pessoa = |n: &str| sim.pessoas[PED_A.iter().position(|p| p.nome == n).unwrap()].clone();
    let idx = |n: &str| PED_A.iter().position(|p| p.nome == n).unwrap();
    // the sites are the union of the files of the persons in the analysis only
    let bruto_de = |nomes: &[&str]| -> BTreeMap<Chave, Vec<i8>> {
        let refs: Vec<&PessoaM> = nomes.iter().map(|n| &sim.modelo[idx(n)]).collect();
        junta_bruto(&refs, true, false)
    };

    // --- trio C1 (child), F, M, with the de novo output
    let saida = t.p("denovo_C1.vcf");
    let op = TrioOptions { de_novo_output: Some(saida.clone()), ..Default::default() };
    let r = trio_analysis(&pessoa("C1"), &pessoa("F"), Some(&pessoa("M")), &op).unwrap();
    let (tab, dn, inf_f, exc_f, inf_m, exc_m) = mendel_bruto(&bruto_de(&["C1", "F", "M"]), 0, 1, Some(2));
    t.eq(&r.genotype_table, &tab, "trio C1/F/M: genotype table");
    let erros: u64 = tab.iter().enumerate().filter(|(i, _)| mendel_code((i / 9) as i8, ((i / 3) % 3) as i8, (i % 3) as i8) > 0).map(|(_, n)| n).sum();
    t.eq(r.errors, erros, "trio C1/F/M: Mendelian errors");
    t.eq(r.sites_checked, tab.iter().sum::<u64>(), "trio C1/F/M: sites checked");
    t.eq(r.errors_by_type.iter().map(|x| x.1).sum::<u64>(), r.errors, "errors by type add up");
    t.eq((r.father_informative, r.father_exclusions, r.mother_informative, r.mother_exclusions), (inf_f, exc_f, inf_m, exc_m), "trio C1/F/M: exclusions");
    t.eq(r.de_novo, dn.len() as u64, "trio C1/F/M: de novo candidates");
    let escritos: Vec<(String, u64)> = le_linhas(&saida)
        .iter()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            assert_eq!(c.len(), 10, "{l}");
            (nomes::canonico(c[0]), c[1].parse().unwrap())
        })
        .collect();
    t.eq(&escritos, &dn, "de novo VCF: exactly the candidates");
    t.eq(r.de_novo_written, dn.len() as u64, "de novo written");
    let plantados = sim.de_novo.iter().filter(|x| dn.contains(x)).count();
    // a planted mutation is missed only when a parent or the child has a missing/erroneous call there
    t.check(plantados >= 26, format!("planted de novo mutations found: {plantados} of 30"));
    let cab = le_linhas(&saida).into_iter().find(|l| l.starts_with("#CHROM")).unwrap();
    t.check(cab.ends_with("\tFORMAT\tC1"), format!("de novo VCF header: {cab}"));
    t.eq(r.verdict, "not_excluded", "trio C1/F/M: verdict");
    t.eq(r.error_level_source, "mother_child", "trio C1/F/M: error level measured on the mother");
    t.eq(r.father_child.relationship, "parent_offspring", "father-child relationship");
    t.eq((r.cpi_sites, r.cpi_sites_without_af, r.log10_cpi), (0, 0, None), "no frequency file: no index");
    t.check(r.resumo().contains("not computed (give a population allele-frequency VCF)"), "no frequency file: the text says how to get the index");

    // --- trio G (with a BED), C1, S1
    let r = trio_analysis(&pessoa("G"), &pessoa("C1"), Some(&pessoa("S1")), &TrioOptions::default()).unwrap();
    let (tab, dn, ..) = mendel_bruto(&bruto_de(&["G", "C1", "S1"]), 0, 1, Some(2));
    t.eq(&r.genotype_table, &tab, "trio G/C1/S1: genotype table");
    t.eq(r.de_novo, dn.len() as u64, "trio G/C1/S1: de novo candidates");

    // --- duo K (a .lz4 with a BED), C2; the right father and a wrong one
    let op = TrioOptions { genotype_error: 0.002, ..Default::default() };
    let r = trio_analysis(&pessoa("K"), &pessoa("C2"), None, &op).unwrap();
    let (tab, _, inf_f, exc_f, ..) = mendel_bruto(&bruto_de(&["K", "C2"]), 0, 1, None);
    t.eq(&r.genotype_table, &tab, "duo K/C2: genotype table");
    t.eq((r.father_informative, r.father_exclusions), (inf_f, exc_f), "duo K/C2: exclusions");
    // duo errors: father 1/1 with child 0/0 (index 2*3+0), father 0/0 with child 1/1 (index 0*3+2)
    t.eq(r.errors_by_type.clone(), vec![(DUO_ERROR_TYPES[0], tab[6]), (DUO_ERROR_TYPES[1], tab[2])], "duo errors by type");
    t.eq(r.errors, tab[6] + tab[2], "duo Mendelian errors");
    t.eq(r.verdict, "not_excluded", &format!("duo K/C2 verdict ({})", r.reason));
    t.eq(r.error_level_source, "genotype_error", "duo: error level from the option");
    let r = trio_analysis(&pessoa("K"), &pessoa("U"), None, &op).unwrap();
    t.eq(r.verdict, "excluded", &format!("duo K/U verdict ({})", r.reason));
    // an uncle as the alleged father (C1 is K's uncle)
    let r = trio_analysis(&pessoa("K"), &pessoa("C1"), None, &op).unwrap();
    t.eq(r.verdict, "excluded", &format!("duo K/uncle verdict ({})", r.reason));
    // wrong mother: the error level is no longer measured on her
    let r = trio_analysis(&pessoa("C1"), &pessoa("F"), Some(&pessoa("U")), &TrioOptions::default()).unwrap();
    t.eq(r.error_level_source, "genotype_error", "trio with a wrong mother: error level");
    t.check(r.reason.contains("mother-child pair does not look like parent and child"), format!("wrong mother noted: {}", r.reason));
}

// =====================================================================
//  SIMULATION C: paternity index against the textbook formulas
// =====================================================================

/// The textbook index, by the obligate paternal allele: with the child
/// homozygous (or heterozygous with a homozygous mother) the father must
/// have given allele x, and PI = P(father transmits x) / freq(x).
fn pi_classico(c: i8, f: i8, m: Option<i8>, p: f64, eps: f64) -> f64 {
    let tf = |x: i8| {
        let t = f64::from(f) / 2.0;
        let t1 = (1.0 - eps) * t + eps * (1.0 - t);
        if x == 1 {
            t1
        } else {
            1.0 - t1
        }
    };
    let freq = |x: i8| if x == 1 { p } else { 1.0 - p };
    match (c, m) {
        (0, _) => tf(0) / freq(0),
        (2, _) => tf(1) / freq(1),
        (1, Some(0)) => tf(1) / freq(1),
        (1, Some(2)) => tf(0) / freq(0),
        (1, Some(_)) => 1.0,
        // duo, child heterozygous: either allele could be the father's
        (1, None) => (tf(1) * freq(0) + tf(0) * freq(1)) / (2.0 * p * (1.0 - p)),
        _ => unreachable!(),
    }
}

fn testa_cpi(t: &mut T, rng: &mut Rng) {
    // GF x GM -> F, B (brothers); F x M -> C; U unrelated
    let nomes = ["GF", "GM", "F", "B", "M", "C", "U"];
    let ns = 20_000usize;
    let mut sitios: Vec<(&str, u64, u8, u8, f64)> = Vec::new();
    for (ci, c) in ["1", "2"].iter().enumerate() {
        let mut pos = 1000u64;
        for _ in 0..ns / 2 {
            pos += 50 + rng.below(1000);
            let rb = rng.base();
            sitios.push((c, pos, rb, rng.outra(rb), rng.freq()));
        }
        let _ = ci;
    }
    let mut h: Vec<[Vec<u8>; 2]> = Vec::new();
    let pais = [None, None, Some((0, 1)), Some((0, 1)), None, Some((2, 4)), None];
    for pp in pais {
        let mut x = [vec![0u8; ns], vec![0u8; ns]];
        for s in 0..ns {
            match pp {
                Some((f, m)) => {
                    x[0][s] = h[f][rng.below(2) as usize][s];
                    x[1][s] = h[m][rng.below(2) as usize][s];
                }
                None => {
                    x[0][s] = u8::from(rng.chance(sitios[s].4));
                    x[1][s] = u8::from(rng.chance(sitios[s].4));
                }
            }
        }
        h.push(x);
    }
    let g = |i: usize, s: usize| (h[i][0][s] + h[i][1][s]) as i8;
    // persons file: C, F, M, B, U, every site, no `chr`
    let cols = [5usize, 2, 4, 3, 6];
    let mut linhas = cabecalho(&[("1", 1 << 30), ("2", 1 << 30)], &cols.map(|i| nomes[i]));
    for (s, x) in sitios.iter().enumerate() {
        let mut l = format!("{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT", x.0, x.1, x.2 as char, x.3 as char);
        for &i in &cols {
            l.push('\t');
            l.push_str(gt_txt(g(i, s), rng));
        }
        linhas.push(l);
    }
    let arq = t.p("trioC.vcf.gz");
    grava(&arq, &linhas);
    // the same, with three planted opposite homozygotes between C and F
    let mut ruim = linhas.clone();
    let mut plantadas = 0;
    let n_cab = ruim.iter().take_while(|l| l.starts_with('#')).count();
    for (s, l) in ruim.iter_mut().skip(n_cab).enumerate() {
        if plantadas < 3 && g(5, s) == 0 && g(2, s) == 1 && g(4, s) == 0 {
            let mut c: Vec<String> = l.split('\t').map(str::to_string).collect();
            c[10] = "1/1".to_string(); // the father becomes 1/1: an exclusion
            *l = c.join("\t");
            plantadas += 1;
        }
    }
    let arq_ruim = t.p("trioC_bad.vcf");
    grava(&arq_ruim, &ruim);

    // frequency file: `chr` names, sites only; AF is wrong on purpose, EUR_AF is right;
    // 5% multiallelic, 3% absent, extra sites nobody has, a chromosome nobody has
    let mut af: HashMap<(String, u64), f64> = HashMap::new();
    let mut la = cabecalho(&[("chr1", 1 << 30), ("chr2", 1 << 30), ("chr3", 1 << 30)], &[]);
    for x in &sitios {
        if rng.chance(0.03) {
            continue;
        }
        let (r, a) = (x.2 as char, x.3 as char);
        if rng.chance(0.05) {
            let outro = (b"ACGT".iter().find(|&&b| b != x.2 && b != x.3).copied().unwrap()) as char;
            if rng.chance(0.5) {
                la.push(format!("chr{}\t{}\t.\t{r}\t{a},{outro}\t.\tPASS\tAF=0.5,0.5;EUR_AF={},0.01", x.0, x.1, x.4));
            } else {
                la.push(format!("chr{}\t{}\t.\t{r}\t{outro},{a}\t.\tPASS\tEUR_AF=0.01,{};AF=0.5,0.5", x.0, x.1, x.4));
            }
        } else {
            la.push(format!("chr{}\t{}\t.\t{r}\t{a}\t.\tPASS\tAF=0.5;EUR_AF={}", x.0, x.1, x.4));
        }
        af.insert((x.0.to_string(), x.1), x.4);
        if rng.chance(0.1) {
            la.push(format!("chr{}\t{}\t.\tA\tC\t.\tPASS\tAF=0.5;EUR_AF=0.3", x.0, x.1 + 1));
        }
    }
    la.push("chr3\t100\t.\tA\tC\t.\tPASS\tEUR_AF=0.2".into());
    let arq_af = t.p("af.vcf.gz");
    grava(&arq_af, &la);

    let pc = Person::new(&arq).sample("C");
    let pf = Person::new(&arq).sample("F");
    let pm = Person::new(&arq).sample("M");
    let pb = Person::new(&arq).sample("B");
    let pu = Person::new(&arq).sample("U");
    t.check(MIN_AF == 0.001, "MIN_AF");
    let esperado = |c: usize, f: usize, m: Option<usize>, eps: f64, arqf: Option<&[String]>| -> (f64, u64, u64) {
        let (mut soma, mut n, mut sem) = (0.0, 0u64, 0u64);
        for (s, x) in sitios.iter().enumerate() {
            let fg = match arqf {
                Some(l) => {
                    let c: Vec<&str> = l[n_cab + s].split('\t').collect();
                    match c[10] {
                        "1/1" | "1|1" => 2,
                        "0/0" | "0|0" => 0,
                        _ => 1,
                    }
                }
                None => g(f, s),
            };
            match af.get(&(x.0.to_string(), x.1)) {
                Some(&p) => {
                    let p = p.clamp(0.001, 0.999);
                    soma += pi_classico(g(c, s), fg, m.map(|m| g(m, s)), p, eps).log10();
                    n += 1;
                }
                None => sem += 1,
            }
        }
        (soma, n, sem)
    };
    let op = |eps: f64| TrioOptions {
        af_path: Some(arq_af.clone()),
        af_field: "EUR_AF".to_string(),
        genotype_error: eps,
        ..Default::default()
    };

    // trio, no error term, clean data: the classic formulas exactly
    let r = trio_analysis(&pc, &pf, Some(&pm), &op(0.0)).unwrap();
    let (l, n, sem) = esperado(5, 2, Some(4), 0.0, None);
    t.eq((r.cpi_sites, r.cpi_sites_without_af), (n, sem), "trio CPI sites");
    t.near(r.log10_cpi.unwrap(), l, 1e-9 * l.abs().max(1.0), "trio log10 CPI, eps = 0");
    t.near(r.probability_of_paternity.unwrap(), 1.0, 1e-12, "trio probability of paternity");
    t.eq(r.verdict, "not_excluded", "trio CPI verdict");
    // the unrelated-man exclusion rate: mean of p^2 (child 0/0) or q^2 (child 1/1)
    let (mut us, mut un) = (0.0, 0u64);
    for (s, x) in sitios.iter().enumerate() {
        if let (Some(&p), c @ (0 | 2)) = (af.get(&(x.0.to_string(), x.1)), g(5, s)) {
            us += if c == 0 { p * p } else { (1.0 - p) * (1.0 - p) };
            un += 1;
        }
    }
    t.near(r.unrelated_expected_rate.unwrap(), us / un as f64, 1e-12, "unrelated exclusion rate");
    println!(
        "  simulation C: trio log10 CPI {:.1} over {} SNVs; unrelated man would be excluded at {:.3} of the informative sites",
        l, n, us / un as f64
    );

    // duo with the error term: the duo formulas, exactly
    for eps in [0.0, 0.001, 0.01] {
        let r = trio_analysis(&pc, &pf, None, &op(eps)).unwrap();
        let (l, n, _) = esperado(5, 2, None, eps, None);
        t.eq(r.cpi_sites, n, "duo CPI sites");
        t.near(r.log10_cpi.unwrap(), l, 1e-9 * l.abs().max(1.0), &format!("duo log10 CPI, eps = {eps}"));
    }
    // three planted exclusions: with eps the index survives, without it is zero
    let pfr = Person::new(&arq_ruim).sample("F");
    let pcr = Person::new(&arq_ruim).sample("C");
    let r = trio_analysis(&pcr, &pfr, None, &op(0.001)).unwrap();
    let (l, ..) = esperado(5, 2, None, 0.001, Some(&ruim));
    t.near(r.log10_cpi.unwrap(), l, 1e-9 * l.abs().max(1.0), "duo with 3 exclusions, eps = 0.001");
    t.check(r.probability_of_paternity.unwrap() > 0.999999, format!("3 errors do not zero the index: {:?}", r.probability_of_paternity));
    t.eq(r.father_exclusions, 3, "three planted exclusions");
    t.eq(r.verdict, "not_excluded", &format!("3 exclusions in ~10000: within the error level ({})", r.reason));
    let r = trio_analysis(&pcr, &pfr, None, &op(0.0)).unwrap();
    t.eq(r.probability_of_paternity, Some(0.0), "without the error term one exclusion zeroes it");
    t.check(r.json().contains("\"log10_cpi\":null"), "-inf is null in JSON");

    // wrong men
    let r = trio_analysis(&pc, &pu, Some(&pm), &op(0.001)).unwrap();
    t.eq(r.verdict, "excluded", &format!("unrelated man ({})", r.reason));
    t.check(r.probability_of_paternity.unwrap() < 1e-6, format!("unrelated man: W = {:?}", r.probability_of_paternity));
    let r = trio_analysis(&pc, &pb, Some(&pm), &op(0.001)).unwrap();
    t.eq(r.verdict, "excluded", &format!("the father's brother ({})", r.reason));
    t.check(r.log10_cpi.unwrap() < -10.0, format!("the father's brother: log10 CPI {:?}", r.log10_cpi));
    t.eq(r.father_child.relationship, "second_degree", "uncle and nephew");
    // the wrong INFO field: AF=0.5 everywhere gives another index
    let r = trio_analysis(&pc, &pf, Some(&pm), &TrioOptions { af_field: "AF".into(), ..op(0.0) }).unwrap();
    t.check((r.log10_cpi.unwrap() - l).abs() > 1.0, "the INFO field is the one asked for");
    let r = trio_analysis(&pc, &pf, Some(&pm), &TrioOptions { af_field: "NOPE".into(), ..op(0.0) }).unwrap();
    t.eq((r.cpi_sites, r.log10_cpi), (0, None), "an INFO field nobody has");
    t.check(r.resumo().contains("no SNV had a population frequency"), "the text says the field matched nothing");
    // formula helpers
    t.near(paternity_index(1, 1, None, 0.3, 0.0), 1.0 / (4.0 * 0.3 * 0.7), 1e-12, "PI het/het duo = 1/4pq");
    t.near(paternity_index(2, 2, Some(1), 0.2, 0.0), 1.0 / 0.2, 1e-12, "PI hom/hom = 1/p");
    t.near(probability_of_paternity(2.0, 0.5), 100.0 / 101.0, 1e-12, "W from log10 CPI");
    t.near(probability_of_paternity(0.0, 0.1), 0.1, 1e-12, "W with CPI 1 equals the prior");
}

// =====================================================================
//  SIMULATION B: runs of homozygosity in a child of first cousins
// =====================================================================

struct SimB {
    croms: Vec<(&'static str, u64)>,
    pos: Vec<Vec<u64>>,
    p: Vec<Vec<f64>>,
}

/// A haplotype: alleles and founder-haplotype labels, per chromosome.
type Hap = (Vec<Vec<u8>>, Vec<Vec<u32>>);

fn fundador_b(sim: &SimB, rotulo: u32, rng: &mut Rng) -> [Hap; 2] {
    let um = |r: u32, rng: &mut Rng| -> Hap {
        let a = sim.p.iter().map(|ps| ps.iter().map(|&p| u8::from(rng.chance(p))).collect()).collect();
        let l = sim.p.iter().map(|ps| vec![r; ps.len()]).collect();
        (a, l)
    };
    [um(rotulo * 2, rng), um(rotulo * 2 + 1, rng)]
}

/// One meiosis with crossovers (on average one per 40 Mb).
fn meiose(sim: &SimB, pai: &[Hap; 2], rng: &mut Rng) -> Hap {
    let mut a = Vec::new();
    let mut l = Vec::new();
    for (c, &(_, len)) in sim.croms.iter().enumerate() {
        let mut cos: Vec<u64> = (0..rng.poisson(len as f64 / 40e6)).map(|_| rng.below(len)).collect();
        cos.sort_unstable();
        let mut h = rng.below(2) as usize;
        let mut k = 0;
        let (mut av, mut lv) = (Vec::new(), Vec::new());
        for (s, &p) in sim.pos[c].iter().enumerate() {
            while k < cos.len() && cos[k] < p {
                h ^= 1;
                k += 1;
            }
            av.push(pai[h].0[c][s]);
            lv.push(pai[h].1[c][s]);
        }
        a.push(av);
        l.push(lv);
    }
    (a, l)
}

/// The runs, by an independent scan: split the calls into chromosome blocks,
/// cut at every het, cut where the gap is too long, keep what is big enough.
fn roh_bruto(
    calls: &[(String, u64, bool)],
    bed: Option<&HashMap<String, Vec<(u64, u64)>>>,
    min_sites: u64,
    min_len: u64,
    max_gap: u64,
) -> Vec<(String, u64, u64, u64)> {
    let dentro = |c: &str, p: u64| bed.map_or(true, |b| b.get(c).is_some_and(|v| v.iter().any(|&(a, z)| a < p && p <= z)));
    let lacuna = |c: &str, a: u64, z: u64| -> u64 {
        match bed {
            None => z - a,
            Some(b) => {
                // walk the bases between a and z, the dumb way: by intervals
                let v = &b[c];
                let mut maior = 0;
                let mut cursor = a; // last covered base (1-based)
                for &(s, e) in v {
                    if e < a || s >= z {
                        continue;
                    }
                    // bases s+1..=e are covered
                    if s + 1 > cursor + 1 {
                        maior = maior.max(s - cursor);
                    }
                    cursor = cursor.max(e);
                }
                maior
            }
        }
    };
    let usados: Vec<&(String, u64, bool)> = calls.iter().filter(|(c, p, _)| dentro(c, *p)).collect();
    let mut blocos: Vec<Vec<&(String, u64, bool)>> = Vec::new();
    for x in usados {
        match blocos.last_mut() {
            Some(b) if b[0].0 == x.0 => b.push(x),
            _ => blocos.push(vec![x]),
        }
    }
    let mut out = Vec::new();
    for b in blocos {
        // hom stretches between hets
        let mut trechos: Vec<Vec<u64>> = vec![Vec::new()];
        for x in &b {
            if x.2 {
                trechos.push(Vec::new());
            } else {
                trechos.last_mut().unwrap().push(x.1);
            }
        }
        for tr in trechos {
            let mut partes: Vec<Vec<u64>> = Vec::new();
            for p in tr {
                match partes.last_mut() {
                    Some(v) if lacuna(&b[0].0, *v.last().unwrap(), p) <= max_gap => v.push(p),
                    _ => partes.push(vec![p]),
                }
            }
            // a gap also cuts between the last hom call and the het that closes it,
            // which does not matter (the het closes the run anyway)
            for v in partes {
                let (ini, fim) = (v[0], *v.last().unwrap());
                if v.len() as u64 >= min_sites && fim - ini + 1 >= min_len {
                    out.push((b[0].0.clone(), ini - 1, fim, v.len() as u64));
                }
            }
        }
    }
    out
}

/// The sliding-window method (plink --homozyg), the slow way: for each call,
/// every window that covers it is counted call by call. Same blocks (by
/// chromosome and by gaps over `max_gap`) as `roh_bruto`.
#[allow(clippy::too_many_arguments)]
fn roh_janela_bruto(
    calls: &[(String, u64, bool)],
    bed: Option<&HashMap<String, Vec<(u64, u64)>>>,
    w: usize,
    h: usize,
    limiar: f64,
    min_sites: u64,
    min_len: u64,
    max_gap: u64,
    densidade: u64,
) -> Vec<(String, u64, u64, u64)> {
    let dentro = |c: &str, p: u64| bed.map_or(true, |b| b.get(c).is_some_and(|v| v.iter().any(|&(a, z)| a < p && p <= z)));
    let lacuna = |c: &str, a: u64, z: u64| -> u64 {
        match bed {
            None => z - a,
            Some(b) => {
                let mut maior = 0;
                let mut cursor = a;
                for &(s, e) in &b[c] {
                    if e < a || s >= z {
                        continue;
                    }
                    if s + 1 > cursor + 1 {
                        maior = maior.max(s - cursor);
                    }
                    cursor = cursor.max(e);
                }
                maior
            }
        }
    };
    let mut blocos: Vec<Vec<&(String, u64, bool)>> = Vec::new();
    for x in calls.iter().filter(|(c, p, _)| dentro(c, *p)) {
        match blocos.last_mut() {
            Some(b) if b[0].0 == x.0 && lacuna(&x.0, b.last().unwrap().1, x.1) <= max_gap => b.push(x),
            _ => blocos.push(vec![x]),
        }
    }
    let mut out = Vec::new();
    for b in blocos {
        let n = b.len();
        if n < w {
            continue;
        }
        let janela_hom = |i: usize| (i..i + w).filter(|&k| b[k].2).count() <= h;
        let cand: Vec<bool> = (0..n)
            .map(|j| {
                let (mut total, mut hom) = (0usize, 0usize);
                for i in j.saturating_sub(w - 1)..=j.min(n - w) {
                    total += 1;
                    if janela_hom(i) {
                        hom += 1;
                    }
                }
                hom > 0 && hom as f64 >= limiar * total as f64
            })
            .collect();
        let mut j = 0;
        while j < n {
            if !cand[j] {
                j += 1;
                continue;
            }
            let ini = j;
            while j < n && cand[j] {
                j += 1;
            }
            let (p0, p1, sites) = (b[ini].1, b[j - 1].1, (j - ini) as u64);
            let comp = p1 - p0 + 1;
            if sites >= min_sites && comp >= min_len && (densidade == 0 || comp <= densidade * sites) {
                out.push((b[0].0.clone(), p0 - 1, p1, sites));
            }
        }
    }
    out
}

fn corridas(r: &RohReport) -> Vec<(String, u64, u64, u64)> {
    r.runs.iter().map(|x| (nomes::canonico(&x.chrom), x.start, x.end, x.sites)).collect()
}

fn testa_roh(t: &mut T, rng: &mut Rng) {
    let croms = vec![("1", 70_000_000u64), ("2", 60_000_000u64)];
    let mut sim = SimB { croms: croms.clone(), pos: Vec::new(), p: Vec::new() };
    for &(_, len) in &croms {
        let mut v = Vec::new();
        let mut q = Vec::new();
        let mut p = 5_000;
        while p < len - 5_000 {
            v.push(p);
            q.push(rng.freq());
            p += 300 + rng.below(2_400);
        }
        sim.pos.push(v);
        sim.p.push(q);
    }
    let refs: Vec<(u8, u8)> = Vec::new();
    let _ = refs;
    let bases: Vec<Vec<(u8, u8)>> = sim.pos.iter().map(|v| v.iter().map(|_| { let r = rng.base(); (r, rng.outra(r)) }).collect()).collect();

    // child of first cousins, with enough autozygosity to test
    let (r_hap, x1, segs) = loop {
        let gf = fundador_b(&sim, 0, rng);
        let gm = fundador_b(&sim, 1, rng);
        let x1 = fundador_b(&sim, 2, rng);
        let y1 = fundador_b(&sim, 3, rng);
        let a = [meiose(&sim, &gf, rng), meiose(&sim, &gm, rng)];
        let b = [meiose(&sim, &gf, rng), meiose(&sim, &gm, rng)];
        let g = [meiose(&sim, &a, rng), meiose(&sim, &x1, rng)];
        let k = [meiose(&sim, &b, rng), meiose(&sim, &y1, rng)];
        let r = [meiose(&sim, &g, rng), meiose(&sim, &k, rng)];
        // true autozygous segments: the two haplotypes carry the same founder haplotype
        let mut segs: Vec<(String, u64, u64)> = Vec::new();
        for (c, &(nome, _)) in sim.croms.iter().enumerate() {
            let mut ini: Option<u64> = None;
            let mut ult = 0;
            for (s, &p) in sim.pos[c].iter().enumerate() {
                let auto = r[0].1[c][s] == r[1].1[c][s];
                match (auto, ini) {
                    (true, None) => ini = Some(p),
                    (false, Some(i)) => {
                        segs.push((nome.to_string(), i, ult));
                        ini = None;
                    }
                    _ => {}
                }
                ult = p;
            }
            if let Some(i) = ini {
                segs.push((nome.to_string(), i, ult));
            }
        }
        let grandes = segs.iter().filter(|s| s.2 - s.1 >= 4_000_000).count();
        let maior = segs.iter().map(|s| s.2 - s.1).max().unwrap_or(0);
        if grandes >= 2 && maior >= 6_000_000 {
            break (r, x1, segs);
        }
    };
    let geno = |h: &[Hap; 2], c: usize, s: usize| (h[0].0[c][s] + h[1].0[c][s]) as i8;

    // R alone (only its variants), R and X1 together (every site, 0/0 included)
    let contigs: Vec<(&str, u64)> = sim.croms.to_vec();
    let mut so_r = cabecalho(&contigs, &["R"]);
    let mut juntos = cabecalho(&contigs, &["X1", "R"]);
    let mut so_x1 = cabecalho(&contigs, &["X1"]);
    let mut calls_r: Vec<(String, u64, bool)> = Vec::new();
    let mut calls_todos: Vec<(String, u64, bool)> = Vec::new();
    for (c, &(nome, _)) in sim.croms.iter().enumerate() {
        for (s, &p) in sim.pos[c].iter().enumerate() {
            let (rb, ab) = (bases[c][s].0 as char, bases[c][s].1 as char);
            let (gr, gx) = (geno(&r_hap, c, s), geno(&x1, c, s));
            if gr > 0 {
                so_r.push(format!("chr{nome}\t{p}\t.\t{rb}\t{ab}\t.\tPASS\t.\tGT\t{}", gt_txt(gr, rng)));
                calls_r.push((nome.to_string(), p, gr == 1));
            }
            if gx > 0 {
                so_x1.push(format!("{nome}\t{p}\t.\t{rb}\t{ab}\t.\tPASS\t.\tGT\t{}", gt_txt(gx, rng)));
            }
            if gr > 0 || gx > 0 {
                juntos.push(format!("{nome}\t{p}\t.\t{rb}\t{ab}\t.\tPASS\t.\tGT\t{}\t{}", gt_txt(gx, rng), gt_txt(gr, rng)));
                calls_todos.push((nome.to_string(), p, gr == 1));
            }
            if s % 97 == 0 && gr > 0 {
                // an indel next to it: not an SNV, ignored
                so_r.push(format!("chr{nome}\t{}\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1", p + 1));
            }
        }
    }
    let arq_r = t.p("roh_R.vcf.gz");
    grava(&arq_r, &so_r);
    let arq_j = t.p("roh_two.vcf");
    grava(&arq_j, &juntos);
    let arq_x = t.p("roh_X1.vcf");
    grava(&arq_x, &so_x1);
    // 1 to 6: the strict method (any het call ends a run); 7 on: the sliding windows
    let op = RohOptions { window_snps: 0, ..Default::default() };

    // 1. against the independent scan, and against the truth
    let saida = t.p("roh_R.bed");
    let rep = runs_of_homozygosity(&Person::new(&arq_r), Some(&saida), &op).unwrap();
    let esp = roh_bruto(&calls_r, None, 50, 1_000_000, 1_000_000);
    t.eq(corridas(&rep), esp.clone(), "ROH of R: the same runs as the independent scan");
    t.eq(rep.non_snv, so_r.iter().filter(|l| l.contains("\tAT\t")).count() as u64, "ROH: indels ignored");
    let lidas: Vec<String> = le_linhas(&saida);
    let esp_bed: Vec<String> = rep.runs.iter().map(|x| format!("{}\t{}\t{}\t{}", x.chrom, x.start, x.end, x.sites)).collect();
    t.eq(lidas, esp_bed, "ROH BED output");
    t.eq(rep.written, rep.count, "ROH lines written");
    let total: u64 = esp.iter().map(|x| x.2 - x.1).sum();
    t.eq(rep.total_length, total, "ROH total length");
    t.near(rep.f_roh, total as f64 / 2_875_001_522.0, 1e-15, "F_ROH over the GRCh38 autosomes");
    t.eq(grch38_autosome_length(), 2_875_001_522, "GRCh38 autosome length");
    t.eq(grch37_autosome_length(), 2_881_033_286, "GRCh37 autosome length");
    let span: u64 = sim.croms.iter().map(|&(n, _)| {
        let v: Vec<u64> = calls_r.iter().filter(|x| x.0 == n).map(|x| x.1).collect();
        v.last().unwrap() - v[0] + 1
    }).sum();
    t.eq(rep.covered_span, span, "ROH covered span");
    t.near(rep.f_roh_covered, total as f64 / span as f64, 1e-15, "F_ROH over the covered span");
    let sobrepoe = |a: (u64, u64), b: (u64, u64)| a.1.min(b.1).saturating_sub(a.0.max(b.0));
    let mut cobertos = 0;
    let grandes: Vec<&(String, u64, u64)> = segs.iter().filter(|s| s.2 - s.1 >= 3_000_000).collect();
    for s in &grandes {
        let cob: u64 = esp.iter().filter(|r| r.0 == s.0).map(|r| sobrepoe((r.1, r.2), (s.1 - 1, s.2))).sum();
        if cob as f64 >= 0.9 * (s.2 - s.1 + 1) as f64 {
            cobertos += 1;
        }
    }
    t.eq(cobertos, grandes.len(), "every true autozygous segment over 3 Mb is >= 90% covered");
    let dentro: u64 = esp
        .iter()
        .map(|r| segs.iter().filter(|s| s.0 == r.0).map(|s| sobrepoe((r.1, r.2), (s.1 - 1, s.2))).sum::<u64>())
        .sum();
    t.check(dentro as f64 >= 0.95 * total as f64, format!("runs fall in true autozygous segments: {dentro} of {total}"));
    let verdade: u64 = segs.iter().map(|s| s.2 - s.1 + 1).sum();
    println!(
        "  simulation B: child of first cousins, {} true autozygous segment(s) ({:.1} Mb); {} run(s) found ({:.1} Mb), F_ROH(covered) {:.4}",
        segs.iter().filter(|s| s.2 - s.1 >= 1_000_000).count(),
        verdade as f64 / 1e6,
        rep.count,
        total as f64 / 1e6,
        rep.f_roh_covered
    );

    // 2. an outbred person: no run
    let rx = runs_of_homozygosity(&Person::new(&arq_x), None, &op).unwrap();
    t.eq(rx.count, 0, "outbred founder: no run of homozygosity");

    // 3. multi-sample: 0/0 calls count as sites
    let rj = runs_of_homozygosity(&Person::new(&arq_j).sample("R"), None, &op).unwrap();
    t.eq(corridas(&rj), roh_bruto(&calls_todos, None, 50, 1_000_000, 1_000_000), "ROH from a multi-sample file");
    t.check(rj.hom_calls > rep.hom_calls, "hom-ref calls of a multi-sample file are used");

    // 4. a callable BED with a 2 Mb hole and a 0.4 Mb hole inside the longest run
    let maior = esp.iter().max_by_key(|x| x.2 - x.1).unwrap().clone();
    let meio = (maior.1 + maior.2) / 2;
    let mut bed: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    for &(n, len) in &sim.croms {
        if n == maior.0 {
            let (b1, b2) = (meio - 1_000_000, meio + 1_000_000);
            let (c1, c2) = (maior.1 + 300_000, maior.1 + 700_000);
            let mut v = vec![(0, c1), (c2, b1), (b2, len)];
            v.retain(|x| x.1 > x.0);
            bed.insert(n.to_string(), v);
        } else {
            bed.insert(n.to_string(), vec![(0, len)]);
        }
    }
    let mut bl = Vec::new();
    for &(n, _) in &sim.croms {
        for &(a, z) in &bed[n] {
            bl.push(format!("{n}\t{a}\t{z}"));
        }
    }
    let arq_bed = t.p("roh_callable.bed");
    grava(&arq_bed, &bl);
    let rb = runs_of_homozygosity(&Person::new(&arq_r).bed(&arq_bed), None, &op).unwrap();
    let esp_b = roh_bruto(&calls_r, Some(&bed), 50, 1_000_000, 1_000_000);
    t.eq(corridas(&rb), esp_b.clone(), "ROH with a callable BED");
    t.check(esp_b.len() == esp.len() + 1, format!("the 2 Mb hole splits a run, the 0.4 Mb one does not: {} -> {}", esp.len(), esp_b.len()));
    t.check(rb.outside_bed > 0, "calls outside the BED are counted");
    // 5. no BED, but no call in 1.5 Mb: a gap
    let (g1, g2) = (meio - 750_000, meio + 750_000);
    let sem: Vec<String> = so_r
        .iter()
        .filter(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            l.starts_with('#') || !(nomes::canonico(c[0]) == maior.0 && (g1..g2).contains(&c[1].parse::<u64>().unwrap()))
        })
        .cloned()
        .collect();
    let arq_gap = t.p("roh_gap.vcf");
    grava(&arq_gap, &sem);
    let calls_gap: Vec<(String, u64, bool)> = calls_r.iter().filter(|x| !(x.0 == maior.0 && (g1..g2).contains(&x.1))).cloned().collect();
    let rg = runs_of_homozygosity(&Person::new(&arq_gap), None, &op).unwrap();
    t.eq(corridas(&rg), roh_bruto(&calls_gap, None, 50, 1_000_000, 1_000_000), "ROH with a 1.5 Mb gap without data");
    let rg2 = runs_of_homozygosity(&Person::new(&arq_gap), None, &RohOptions { max_gap: 2_000_000, ..op.clone() }).unwrap();
    let esp_g2 = roh_bruto(&calls_gap, None, 50, 1_000_000, 2_000_000);
    t.eq(corridas(&rg2), esp_g2.clone(), "ROH with max_gap 2 Mb");
    t.eq(esp_g2.iter().map(|x| (x.1, x.2)).collect::<Vec<_>>(), esp.iter().map(|x| (x.1, x.2)).collect::<Vec<_>>(), "a larger max_gap bridges the 1.5 Mb gap");
    // 6. other thresholds
    for (ms, ml) in [(10, 200_000), (200, 3_000_000), (1, 1)] {
        let o = RohOptions { min_sites: ms, min_length: ml, ..op.clone() };
        let r = runs_of_homozygosity(&Person::new(&arq_r), None, &o).unwrap();
        t.eq(corridas(&r), roh_bruto(&calls_r, None, ms, ml, 1_000_000), &format!("ROH with min_sites {ms}, min_length {ml}"));
    }
    t.check(rep.resumo().contains("F_ROH") && rep.json().contains("\"runs\":[["), "ROH text and JSON");

    // 7. the sliding windows of plink --homozyg (the default), against the slow way
    let jan = RohOptions::default();
    t.eq((jan.window_snps, jan.window_het, jan.window_missing, jan.max_density), (50, 1, 5, 50_000), "window defaults = plink's");
    let rw = runs_of_homozygosity(&Person::new(&arq_r), None, &jan).unwrap();
    let esp_w = roh_janela_bruto(&calls_r, None, 50, 1, 0.05, 50, 1_000_000, 1_000_000, 50_000);
    t.eq(corridas(&rw), esp_w.clone(), "ROH by windows: the same runs as the slow way");
    t.check(rw.method == "window" && rep.method == "strict", "the report says the method");
    t.check(rw.json().contains("\"method\":\"window\"") && rw.resumo().contains("plink --homozyg"), "window method in the JSON and the text");
    let total_w: u64 = esp_w.iter().map(|x| x.2 - x.1).sum();
    // with no genotyping error the two agree, up to the run ends (where few windows are homozygous)
    t.check(
        (total_w as f64 - total as f64).abs() <= 0.01 * total as f64,
        format!("no genotyping error: windows and the strict scan agree within 1%: {total_w} vs {total}"),
    );
    let mut cob_w = 0;
    for s in &grandes {
        let cob: u64 = esp_w.iter().filter(|r| r.0 == s.0).map(|r| sobrepoe((r.1, r.2), (s.1 - 1, s.2))).sum();
        if cob as f64 >= 0.9 * (s.2 - s.1 + 1) as f64 {
            cob_w += 1;
        }
    }
    t.eq(cob_w, grandes.len(), "windows: every true autozygous segment over 3 Mb is >= 90% covered");
    let dentro_w: u64 = esp_w
        .iter()
        .map(|r| segs.iter().filter(|s| s.0 == r.0).map(|s| sobrepoe((r.1, r.2), (s.1 - 1, s.2))).sum::<u64>())
        .sum();
    t.check(dentro_w as f64 >= 0.9 * total_w as f64, format!("window runs fall in true autozygous segments: {dentro_w} of {total_w}"));
    println!("  windows: {} run(s), {:.1} Mb (strict: {} run(s), {:.1} Mb)", rw.count, total_w as f64 / 1e6, rep.count, total as f64 / 1e6);
    let rxw = runs_of_homozygosity(&Person::new(&arq_x), None, &jan).unwrap();
    let calls_x: Vec<(String, u64, bool)> = so_x1
        .iter()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (c[0].to_string(), c[1].parse().unwrap(), c[9].chars().next() != c[9].chars().nth(2))
        })
        .collect();
    t.eq(corridas(&rxw), roh_janela_bruto(&calls_x, None, 50, 1, 0.05, 50, 1_000_000, 1_000_000, 50_000), "windows on the outbred founder = the slow way");
    let rbw = runs_of_homozygosity(&Person::new(&arq_r).bed(&arq_bed), None, &jan).unwrap();
    t.eq(corridas(&rbw), roh_janela_bruto(&calls_r, Some(&bed), 50, 1, 0.05, 50, 1_000_000, 1_000_000, 50_000), "windows with a callable BED");
    for mg in [1_000_000u64, 2_000_000] {
        let rgw = runs_of_homozygosity(&Person::new(&arq_gap), None, &RohOptions { max_gap: mg, ..jan.clone() }).unwrap();
        t.eq(corridas(&rgw), roh_janela_bruto(&calls_gap, None, 50, 1, 0.05, 50, 1_000_000, mg, 50_000), &format!("windows with a 1.5 Mb gap, max_gap {mg}"));
    }
    for (w, h, lim, ms, dens) in [(20u64, 0u64, 0.05, 20u64, 50_000u64), (50, 2, 0.2, 50, 50_000), (100, 1, 0.05, 100, 0), (30, 1, 0.5, 10, 20_000)] {
        let o = RohOptions { window_snps: w, window_het: h, window_threshold: lim, min_sites: ms, max_density: dens, ..jan.clone() };
        let r = runs_of_homozygosity(&Person::new(&arq_r), None, &o).unwrap();
        t.eq(
            corridas(&r),
            roh_janela_bruto(&calls_r, None, w as usize, h as usize, lim, ms, 1_000_000, 1_000_000, dens),
            &format!("windows of {w}, {h} het, threshold {lim}, min_sites {ms}, density {dens}"),
        );
    }

    // 8. one genotyping error in the middle of the longest run: strict splits it, windows do not
    let longa = esp_w.iter().max_by_key(|x| x.2 - x.1).unwrap().clone();
    let meio_l = (longa.1 + longa.2) / 2;
    let alvo = calls_r.iter().filter(|x| x.0 == longa.0 && !x.2 && x.1 >= meio_l).map(|x| x.1).next().unwrap();
    let com_erro: Vec<String> = so_r
        .iter()
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            if !l.starts_with('#') && nomes::canonico(c[0]) == longa.0 && c[1] == alvo.to_string() && c[3].len() == 1 && c[4].len() == 1 {
                let mut c2: Vec<String> = c.iter().map(|s| s.to_string()).collect();
                c2[9] = "0/1".to_string();
                c2.join("\t")
            } else {
                l.clone()
            }
        })
        .collect();
    let arq_err = t.p("roh_error.vcf");
    grava(&arq_err, &com_erro);
    let cobre = |r: &RohReport| r.runs.iter().any(|x| nomes::canonico(&x.chrom) == longa.0 && x.start < meio_l - 1_000_000 && x.end > meio_l + 1_000_000);
    let erro_w = runs_of_homozygosity(&Person::new(&arq_err), None, &jan).unwrap();
    let erro_s = runs_of_homozygosity(&Person::new(&arq_err), None, &op).unwrap();
    t.check(cobre(&erro_w), format!("windows: one wrong het call at {alvo} does not split the run {:?}", longa));
    t.check(!erro_s.runs.iter().any(|x| nomes::canonico(&x.chrom) == longa.0 && x.start < alvo && x.end > alvo), "strict: the wrong het call splits the run");

    // 9. missing calls and the density limit, worked out by hand: 121 hom
    // calls 20 kb apart, the middle one missing
    let mut linhas = cabecalho(&[("1", 10_000_000)], &["S"]);
    for k in 0..121u64 {
        let gt = if k == 60 { "./." } else { "1/1" };
        linhas.push(format!("chr1\t{}\t.\tA\tG\t.\tPASS\t.\tGT\t{gt}", 1_000_000 + 20_000 * k));
    }
    let arq_m = t.p("roh_missing.vcf");
    grava(&arq_m, &linhas);
    let pos = |k: u64| 1_000_000 + 20_000 * k;
    let rm = runs_of_homozygosity(&Person::new(&arq_m), None, &jan).unwrap();
    t.eq(corridas(&rm), vec![("1".to_string(), pos(0) - 1, pos(120), 121)], "one missing call in a window of 50: one run");
    let rm0 = runs_of_homozygosity(&Person::new(&arq_m), None, &RohOptions { window_missing: 0, ..jan.clone() }).unwrap();
    // the calls 58 to 62 have fewer than 3 of their 50 windows free of the missing call
    t.eq(
        corridas(&rm0),
        vec![("1".to_string(), pos(0) - 1, pos(57), 58), ("1".to_string(), pos(63) - 1, pos(120), 58)],
        "no missing call allowed: the run splits around it (threshold 0.05)",
    );
    let rd = runs_of_homozygosity(&Person::new(&arq_m), None, &RohOptions { max_density: 10_000, ..jan.clone() }).unwrap();
    t.eq(rd.count, 0, "a run with one call per 20 kb is dropped at max_density 10 kb");
    t.err(
        runs_of_homozygosity(&Person::new(&arq_m), None, &RohOptions { window_threshold: 1.5, ..jan.clone() }),
        &["window_threshold", "between 0 and 1"],
        "window_threshold over 1",
    );
}

// =====================================================================
//  GENETIC SEX
// =====================================================================

fn testa_sexo(t: &mut T, rng: &mut Rng) {
    #[derive(Clone, Copy, PartialEq)]
    enum X {
        Dois,
        UmHap,
        UmDip,
        Nenhum,
        Poucos,
    }
    #[derive(Clone, Copy, PartialEq)]
    enum Y {
        Sim,
        Nao,
    }
    let auto: Vec<(u64, f64)> = (0..4000).map(|i| (10_000 + i * 5_000, rng.freq())).collect();
    let mut x_sitios: Vec<(u64, f64)> = Vec::new();
    for _ in 0..200 {
        x_sitios.push((10_001 + rng.below(2_689_000), rng.freq())); // PAR1 in both builds
    }
    for _ in 0..150 {
        x_sitios.push((2_699_521 + rng.below(81_959), rng.freq())); // PAR1 of GRCh38 only
    }
    for _ in 0..3000 {
        x_sitios.push((2_800_000 + rng.below(152_000_000), rng.freq()));
    }
    for _ in 0..50 {
        x_sitios.push((155_701_383 + rng.below(329_000), rng.freq())); // PAR2 of GRCh38
    }
    x_sitios.sort_by_key(|x| x.0);
    x_sitios.dedup_by_key(|x| x.0);
    let y_sitios: Vec<(u64, f64)> = (0..400).map(|i| (3_000_000 + i * 40_000, 0.3)).chain((0..30).map(|i| (20_000 + i * 50_000, 0.4))).collect();
    let mut y_sitios = y_sitios;
    y_sitios.sort_by_key(|x| x.0);

    let par38 = |p: u64| (10_001..=2_781_479).contains(&p) || (155_701_383..=156_030_895).contains(&p);
    let par37 = |p: u64| (60_001..=2_699_520).contains(&p) || (154_931_044..=155_260_560).contains(&p);
    let ypar = |p: u64| (10_001..=2_781_479).contains(&p);

    struct Esp {
        x_het: u64,
        x_hom_alt: u64,
        x_hap: u64,
        y_calls: u64,
        y_nonref: u64,
    }
    let dir = t.dir.clone();
    let caso = |nome: &str, x: X, y: Y, x_len: Option<u64>, declara_y: bool, build37: bool, rng: &mut Rng| -> (String, Esp) {
        let mut contigs: Vec<(&str, u64)> = vec![("chr1", 248_956_422)];
        if let Some(l) = x_len {
            contigs.push(("chrX", l));
        }
        if declara_y {
            contigs.push(("chrY", 57_227_415));
        }
        let mut l = cabecalho(&contigs, &[nome]);
        for &(p, f) in &auto {
            let g = u8::from(rng.chance(f)) + u8::from(rng.chance(f));
            if g > 0 {
                l.push(format!("chr1\t{p}\t.\tA\tG\t.\tPASS\t.\tGT\t{}", gt_txt(g as i8, rng)));
            }
        }
        let mut e = Esp { x_het: 0, x_hom_alt: 0, x_hap: 0, y_calls: 0, y_nonref: 0 };
        let par = |p: u64| if build37 { par37(p) } else { par38(p) };
        let mut nx = 0;
        for &(p, f) in &x_sitios {
            if x == X::Nenhum || (x == X::Poucos && nx >= 40) {
                break;
            }
            let diploide = x == X::Dois || x == X::Poucos || par(p) || (p >= 2_699_521 && p <= 2_781_479);
            if diploide {
                let g = u8::from(rng.chance(f)) + u8::from(rng.chance(f));
                if g > 0 {
                    l.push(format!("chrX\t{p}\t.\tC\tT\t.\tPASS\t.\tGT\t{}", gt_txt(g as i8, rng)));
                    nx += 1;
                    if !par(p) {
                        if g == 1 {
                            e.x_het += 1;
                        } else {
                            e.x_hom_alt += 1;
                        }
                    }
                }
            } else if rng.chance(f) {
                nx += 1;
                if x == X::UmHap {
                    l.push(format!("chrX\t{p}\t.\tC\tT\t.\tPASS\t.\tGT\t1"));
                    e.x_hap += 1;
                } else {
                    l.push(format!("chrX\t{p}\t.\tC\tT\t.\tPASS\t.\tGT\t1/1"));
                    e.x_hom_alt += 1;
                }
            }
        }
        if y == Y::Sim {
            for &(p, f) in &y_sitios {
                if ypar(p) {
                    let g = u8::from(rng.chance(f)) + u8::from(rng.chance(f));
                    if g > 0 {
                        l.push(format!("chrY\t{p}\t.\tG\tA\t.\tPASS\t.\tGT\t{}", gt_txt(g as i8, rng)));
                    }
                } else {
                    e.y_calls += 1;
                    let alt = rng.chance(f);
                    e.y_nonref += u64::from(alt);
                    let gt = match (x, alt) {
                        (X::UmHap, true) => "1",
                        (X::UmHap, false) => "0",
                        (_, true) => "1/1",
                        (_, false) => "0/0",
                    };
                    l.push(format!("chrY\t{p}\t.\tG\tA\t.\tPASS\t.\tGT\t{gt}"));
                }
            }
        }
        let caminho = dir.join(format!("sex_{nome}.vcf")).to_string_lossy().to_string();
        grava(&caminho, &l);
        (caminho, e)
    };
    let x38 = Some(156_040_895);
    let casos = [
        ("female", X::Dois, Y::Nao, x38, true, "XX", vec!["no_y_data"]),
        ("male_haploid", X::UmHap, Y::Sim, x38, true, "XY", vec!["haploid_x_calls"]),
        ("male_diploid", X::UmDip, Y::Sim, x38, true, "XY", vec![]),
        ("xxy", X::Dois, Y::Sim, x38, true, "undetermined", vec!["xxy_pattern"]),
        ("x0", X::UmDip, Y::Nao, x38, false, "undetermined", vec!["no_y_data", "x0_pattern"]),
        ("no_x_y", X::Nenhum, Y::Nao, None, false, "undetermined", vec!["no_y_data", "no_x_data"]),
        ("few_x", X::Poucos, Y::Nao, x38, false, "undetermined", vec!["no_y_data", "no_x_data"]),
    ];
    for (nome, x, y, xl, dy, veredito, flags) in casos {
        let (arq, e) = caso(nome, x, y, xl, dy, false, rng);
        let r = genetic_sex(&Person::new(&arq), None, false).unwrap();
        t.eq(r.verdict, veredito, &format!("sex of {nome} ({})", r.reason));
        t.eq(&r.flags, &flags, &format!("flags of {nome}"));
        t.eq((r.x_het, r.x_hom_alt, r.x_haploid, r.y_calls, r.y_non_ref), (e.x_het, e.x_hom_alt, e.x_hap, e.y_calls, e.y_nonref), &format!("sex counts of {nome}"));
        t.eq((r.build, r.build_source), ("GRCh38", if xl.is_some() { "header" } else { "default" }), &format!("build of {nome}"));
        if nome == "no_x_y" {
            t.check(r.reason.contains("no chrX data") && r.reason.contains("no chrY data"), format!("no X/Y said clearly: {}", r.reason));
        }
        if nome == "few_x" {
            t.check(r.reason.contains("only") && r.reason.contains("at least 100"), format!("few X calls: {}", r.reason));
        }
    }
    // the PAR boundary moves with the build: hets in chrX:2,699,521-2,781,479 of a male
    let (arq, e37) = caso("male_grch37", X::UmDip, Y::Sim, Some(155_270_560), true, true, rng);
    let r = genetic_sex(&Person::new(&arq), None, false).unwrap();
    t.eq((r.build, r.build_source), ("GRCh37", "header"), "GRCh37 from the chrX length");
    t.eq(r.x_het, e37.x_het, "GRCh37: hets just past PAR1 are outside the PAR");
    t.check(e37.x_het > 20, "the boundary test has hets");
    let r38 = genetic_sex(&Person::new(&arq), Some("hg38"), false).unwrap();
    t.eq((r38.build, r38.build_source, r38.x_het), ("GRCh38", "argument", 0), "GRCh38 by argument: the same hets are in PAR1");
    t.check(r38.x_par_calls > r.x_par_calls, "more PAR calls with GRCh38");
    t.err(genetic_sex(&Person::new(&arq), Some("hg17"), false), &["hg17", "GRCh38"], "unknown build");
    t.check(r.resumo().contains("Chromosomal sex") && r.json().contains("\"verdict\""), "sex text and JSON");
}

// =====================================================================
//  ORDER, NAMES AND ERRORS
// =====================================================================

fn testa_erros(t: &mut T) {
    let cab = cabecalho(&[], &["S"]);
    let dir = t.dir.clone();
    let um = |linhas: &[&str], nome: &str| -> String {
        let mut l = cab.clone();
        l.extend(linhas.iter().map(|s| s.to_string()));
        let c = dir.join(nome).to_string_lossy().to_string();
        grava(&c, &l);
        c
    };
    let rec = |c: &str, p: u64| format!("{c}\t{p}\t.\tA\tG\t.\tPASS\t.\tGT\t0/1");
    let ok = um(&[&rec("chr1", 100), &rec("chr2", 100)], "ok.vcf");
    let ord = |c: &str| -> String {
        let r = relatedness(&[Person::new(&ok), Person::new(c)], &JoinOptions::default());
        r.map(|_| "ok".to_string()).unwrap_or_else(|e| e)
    };
    let a = um(&[&rec("chr2", 100), &rec("chr1", 100)], "order_chrom.vcf");
    let e = ord(&a);
    t.check(e.contains("sort_genomic") && e.contains("chrom_order=\"natural\"") && e.contains("order_chrom.vcf") && e.contains("line 5"), format!("chromosome order: {e}"));
    let a = um(&[&rec("chr1", 200), &rec("chr1", 100)], "order_pos.vcf");
    let e = ord(&a);
    t.check(e.contains("not sorted by position") && e.contains("sort_genomic"), format!("position order: {e}"));
    let a = um(&[&rec("chr1", 100), &rec("chr10", 100), &rec("chr2", 100)], "order_lexical.vcf");
    t.check(ord(&a).contains("comes after `chr10`"), "lexical order is refused");
    let a = um(&[&rec("chr1", 100), &rec("chr2", 100), &rec("chr1", 300)], "order_twice.vcf");
    t.check(ord(&a).contains("sort_genomic"), "a chromosome in two blocks is refused");
    let a = um(&[&rec("chrM", 5), &rec("1", 100), &rec("2", 50), &rec("chrX", 10), &rec("chrY", 10), &rec("MT", 1), &rec("GL000192.1", 1)], "order_hg19.vcf");
    t.eq(ord(&a), "ok".to_string(), "hg19 order (chrM first), mixed names, scaffolds: accepted");
    let r = relatedness(&[Person::new(&ok), Person::new(&a)], &JoinOptions::default()).unwrap();
    t.eq(r.scan.other_sequences, 3, "chrM, MT and a scaffold are left out and counted");

    // persons and files
    let multi = t.p("multi.vcf");
    let mut l = cabecalho(&[], &["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]);
    l.push("1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1".into());
    grava(&multi, &l);
    t.err(relatedness(&[Person::new(&ok), Person::new(&multi)], &JoinOptions::default()), &["10 samples", "sample=", "A, B, C", "2 more", "multi.vcf"], "multi-sample file without a sample");
    t.err(relatedness(&[Person::new(&ok), Person::new(&multi).sample("Z")], &JoinOptions::default()), &["`Z`", "multi.vcf"], "sample not found");
    t.eq(vcf_samples(&multi).unwrap().len(), 10, "vcf_samples");
    let sitios = t.p("sites_only.vcf");
    let mut l = cabecalho(&[], &[]);
    l.push("1\t100\t.\tA\tG\t.\tPASS\t.".into());
    grava(&sitios, &l);
    t.err(relatedness(&[Person::new(&ok), Person::new(&sitios)], &JoinOptions::default()), &["sites-only"], "sites-only VCF as a person");
    t.err(relatedness(&[Person::new(&ok), Person::new(&t.p("nope.vcf"))], &JoinOptions::default()), &["nope.vcf"], "missing file");
    let bed_txt = t.p("not_a_vcf.vcf");
    grava(&bed_txt, &["chr1\t100\t200".to_string()]);
    t.err(relatedness(&[Person::new(&ok), Person::new(&bed_txt)], &JoinOptions::default()), &["#CHROM", "not_a_vcf.vcf"], "not a VCF");
    let bed_ruim = t.p("bad.bed");
    grava(&bed_ruim, &["chr1\t100\t200".to_string(), "chr1\tabc\t300".to_string()]);
    t.err(relatedness(&[Person::new(&ok), Person::new(&ok).bed(&bed_ruim)], &JoinOptions::default()), &["bad.bed", "line 2"], "bad BED line");
    t.err(relatedness(&[Person::new(&ok)], &JoinOptions::default()), &["at least two"], "one person");
    t.err(
        trio_analysis(&Person::new(&ok), &Person::new(&ok), None, &TrioOptions { de_novo_output: Some(t.p("x.vcf")), ..Default::default() }),
        &["mother"],
        "de novo without a mother",
    );
    t.err(
        trio_analysis(&Person::new(&ok), &Person::new(&ok), Some(&Person::new(&ok)), &TrioOptions { de_novo_output: Some(ok.clone()), ..Default::default() }),
        &["one of the input files"],
        "de novo output over an input",
    );
    t.err(trio_analysis(&Person::new(&ok), &Person::new(&ok), None, &TrioOptions { genotype_error: 0.7, ..Default::default() }), &["genotype_error"], "genotype_error out of range");
    t.err(trio_analysis(&Person::new(&ok), &Person::new(&ok), None, &TrioOptions { prior: 1.0, ..Default::default() }), &["prior"], "prior out of range");
    t.err(runs_of_homozygosity(&Person::new(&ok), None, &RohOptions { min_sites: 0, ..Default::default() }), &["min_sites"], "min_sites 0");
    t.err(runs_of_homozygosity(&Person::new(&ok), Some(&ok), &RohOptions::default()), &["one of the input files"], "ROH output over the input");
    // an error in the middle of the scan leaves no half-written output
    let a = um(&[&rec("chr1", 100), &rec("chr3", 100), &rec("chr2", 100)], "order_late.vcf");
    let saida = t.p("half.vcf");
    let e = trio_analysis(&Person::new(&a), &Person::new(&ok), Some(&Person::new(&ok)), &TrioOptions { de_novo_output: Some(saida.clone()), ..Default::default() });
    t.check(e.is_err() && !std::path::Path::new(&saida).exists(), "a failed trio analysis leaves no output");

    // the GT spellings
    let gts = t.p("gts.vcf");
    let mut l = cabecalho(&[], &["P", "Q"]);
    let casos = [("0/1", 1), ("1|0", 1), ("0|1", 1), ("1/1", 2), ("0/0", 0), ("./.", -1), (".", -1), ("0/.", -1), ("1", -1), ("1/2", -1), ("0/1/1", -1), ("x/1", -1), ("01/1", 2)];
    for (i, (g, _)) in casos.iter().enumerate() {
        l.push(format!("1\t{}\t.\tA\tG\t.\tPASS\t.\tGT:DP\t{g}:5\t0/0:7", 100 + i));
    }
    l.push("1\t500\t.\tA\tG\t.\tPASS\t.\tDP:GT\t5:1/1\t7:0/0".into());
    grava(&gts, &l);
    let r = relatedness(&[Person::new(&gts).sample("P"), Person::new(&gts).sample("Q")], &JoinOptions::default()).unwrap();
    let esp_het = casos.iter().filter(|c| c.1 == 1).count() as u64;
    let esp_alt = casos.iter().filter(|c| c.1 == 2).count() as u64 + 1;
    let esp_falta = casos.iter().filter(|c| c.1 == -1).count() as u64;
    t.eq((r.persons[0].het, r.persons[0].hom_alt, r.persons[0].missing), (esp_het, esp_alt, esp_falta), "GT spellings, GT not first in FORMAT");
}

// =====================================================================
//  REAL DATA
// =====================================================================

fn pessoa_txt(s: &str) -> Person {
    let c: Vec<&str> = s.split('|').collect();
    let mut p = Person::new(c[0]);
    if let Some(x) = c.get(1).filter(|x| !x.is_empty()) {
        p = p.sample(x);
    }
    if let Some(x) = c.get(2).filter(|x| !x.is_empty()) {
        p = p.bed(x);
    }
    if let Some(x) = c.get(3).filter(|x| !x.is_empty()) {
        p = p.label(x);
    }
    p
}

fn real(args: &[String]) {
    let t0 = Instant::now();
    match args.first().map(String::as_str) {
        Some("relatedness") => {
            let ps: Vec<Person> = args[1..].iter().map(|s| pessoa_txt(s)).collect();
            let r = relatedness(&ps, &JoinOptions::default()).unwrap_or_else(|e| panic!("{e}"));
            print!("{}", r.resumo());
        }
        Some("trio") => {
            let filho = pessoa_txt(&args[1]);
            let pai = pessoa_txt(&args[2]);
            let mae = args.get(3).filter(|s| *s != "-" && !s.starts_with("--")).map(|s| pessoa_txt(s));
            let mut op = TrioOptions::default();
            let mut i = 3;
            while i < args.len() {
                match args[i].as_str() {
                    "--af" => {
                        op.af_path = Some(args[i + 1].clone());
                        op.af_field = args[i + 2].clone();
                        i += 3;
                    }
                    "--denovo" => {
                        op.de_novo_output = Some(args[i + 1].clone());
                        i += 2;
                    }
                    "--error" => {
                        op.genotype_error = args[i + 1].parse().unwrap();
                        i += 2;
                    }
                    _ => i += 1,
                }
            }
            let r = trio_analysis(&filho, &pai, mae.as_ref(), &op).unwrap_or_else(|e| panic!("{e}"));
            print!("{}", r.resumo());
            let tipos: Vec<String> = r.errors_by_type.iter().map(|(n, c)| format!("{n}={c}")).collect();
            println!("  [types] {}", tipos.join(" "));
        }
        Some("roh") => {
            let p = pessoa_txt(&args[1]);
            let mut op = RohOptions::default();
            let mut saida = None;
            let mut i = 2;
            while i < args.len() {
                match args[i].as_str() {
                    "--min-length" => op.min_length = args[i + 1].parse().unwrap(),
                    "--min-sites" => op.min_sites = args[i + 1].parse().unwrap(),
                    "--max-gap" => op.max_gap = args[i + 1].parse().unwrap(),
                    s => {
                        saida = Some(s.to_string());
                        i += 1;
                        continue;
                    }
                }
                i += 2;
            }
            let r = runs_of_homozygosity(&p, saida.as_deref(), &op).unwrap_or_else(|e| panic!("{e}"));
            print!("{}", r.resumo());
        }
        Some("sex") => {
            let p = pessoa_txt(&args[1]);
            let pass = args.iter().any(|a| a == "--pass");
            let build = args.get(2).filter(|a| !a.starts_with("--")).map(String::as_str);
            let r = genetic_sex(&p, build, pass).unwrap_or_else(|e| panic!("{e}"));
            print!("{}", r.resumo());
        }
        Some("intersect") => {
            // out.bed in1.bed in2.bed ...: the regions present in every input
            let le = |c: &str| -> BTreeMap<String, Vec<(u64, u64)>> {
                let mut m: BTreeMap<String, Vec<(u64, u64)>> = BTreeMap::new();
                for l in le_linhas(c) {
                    let x: Vec<&str> = l.split('\t').collect();
                    if x.len() >= 3 && !l.starts_with('#') {
                        m.entry(x[0].to_string()).or_default().push((x[1].parse().unwrap(), x[2].parse().unwrap()));
                    }
                }
                for v in m.values_mut() {
                    v.sort_unstable();
                    let mut f: Vec<(u64, u64)> = Vec::new();
                    for &(a, b) in v.iter() {
                        match f.last_mut() {
                            Some(u) if a <= u.1 => u.1 = u.1.max(b),
                            _ => f.push((a, b)),
                        }
                    }
                    *v = f;
                }
                m
            };
            let mut acc = le(&args[2]);
            for c in &args[3..] {
                let b = le(c);
                let mut novo = BTreeMap::new();
                for (k, va) in &acc {
                    let Some(vb) = b.get(k) else { continue };
                    let (mut i, mut j, mut v) = (0, 0, Vec::new());
                    while i < va.len() && j < vb.len() {
                        let (s, e) = (va[i].0.max(vb[j].0), va[i].1.min(vb[j].1));
                        if s < e {
                            v.push((s, e));
                        }
                        if va[i].1 < vb[j].1 {
                            i += 1;
                        } else {
                            j += 1;
                        }
                    }
                    novo.insert(k.clone(), v);
                }
                acc = novo;
            }
            let mut linhas = Vec::new();
            let mut total = 0u64;
            let mut chaves: Vec<&String> = acc.keys().collect();
            chaves.sort_by(|a, b| externo::natural_cmp(a, b));
            for k in chaves {
                for &(a, b) in &acc[k] {
                    linhas.push(format!("{k}\t{a}\t{b}"));
                    total += b - a;
                }
            }
            grava(&args[1], &linhas);
            println!("  {} regions, {:.1} Mb", linhas.len(), total as f64 / 1e6);
        }
        _ => panic!("real: relatedness | trio | roh | sex"),
    }
    println!("  [time] {:.1} s", t0.elapsed().as_secs_f64());
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().map(String::as_str) == Some("real") {
        real(&args[1..]);
        return;
    }
    let dir = std::env::temp_dir().join(format!("test_parentesco_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut t = T { n: 0, falhas: Vec::new(), dir: dir.clone() };
    let mut rng = Rng(20261114);
    let t0 = Instant::now();

    let sim = simula_a(&t, &mut rng);
    testa_parentesco(&mut t, &sim);
    testa_trio_a(&mut t, &sim);
    testa_cpi(&mut t, &mut rng);
    testa_roh(&mut t, &mut rng);
    testa_sexo(&mut t, &mut rng);
    testa_erros(&mut t);

    let _ = std::fs::remove_dir_all(&dir);
    println!("{} of {} checks passed ({:.1} s)", t.n - t.falhas.len(), t.n, t0.elapsed().as_secs_f64());
    if !t.falhas.is_empty() {
        std::process::exit(1);
    }
}
