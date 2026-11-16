//! Ensaios do `src/conjuntos.rs` (teste de conjuntos de genes sobre uma tabela de
//! expressão diferencial).
//!
//!   cargo run --release --bin conjuntos

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
#[path = "../../../src/expressao.rs"]
#[allow(dead_code)]
mod expressao;
#[path = "../../../src/conjuntos.rs"]
mod conjuntos;

use std::path::PathBuf;

struct H {
    ok: usize,
    fail: usize,
}

impl H {
    fn check(&mut self, cond: bool, what: &str) {
        if cond {
            self.ok += 1;
        } else {
            self.fail += 1;
            println!("FAIL  {what}");
        }
    }
}

/// Mann-Whitney z by counting pairs (U = pairs where the member is higher, ties
/// count half), with the same tie-corrected variance: an independent check.
fn z_by_pairs(v: &[f64], members: &[usize]) -> f64 {
    let set: std::collections::HashSet<usize> = members.iter().copied().collect();
    let others: Vec<usize> = (0..v.len()).filter(|i| !set.contains(i)).collect();
    let mut u = 0.0;
    for &a in members {
        for &b in &others {
            if v[a] > v[b] {
                u += 1.0;
            } else if v[a] == v[b] {
                u += 0.5;
            }
        }
    }
    let n = v.len() as f64;
    let (n1, n2) = (members.len() as f64, others.len() as f64);
    let mut counts: std::collections::HashMap<u64, f64> = std::collections::HashMap::new();
    for x in v {
        *counts.entry(x.to_bits()).or_insert(0.0) += 1.0;
    }
    let ties: f64 = counts.values().map(|t| t * t * t - t).sum();
    let var = n1 * n2 / 12.0 * ((n + 1.0) - ties / (n * (n - 1.0)));
    (u - n1 * n2 / 2.0) / var.sqrt()
}

fn lcg(state: &mut u64) -> f64 {
    *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / ((1u64 << 53) as f64)
}

fn main() {
    let mut h = H { ok: 0, fail: 0 };
    let dir: PathBuf = std::env::temp_dir().join(format!("mmgen_conjuntos_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let p = |n: &str| dir.join(n).to_string_lossy().to_string();

    // 1. Ranks with ties, and the rank-sum z against the pair count.
    let (r, t) = conjuntos::ranks(&[3.0, 1.0, 3.0, 2.0]);
    h.check(r == vec![3.5, 1.0, 3.5, 2.0] && t == 6.0, "average ranks and tie term");
    let mut seed = 7u64;
    let mut bad = 0;
    for _ in 0..200 {
        let n = 5 + (lcg(&mut seed) * 60.0) as usize;
        // values on a coarse grid, so that there are ties
        let v: Vec<f64> = (0..n).map(|_| ((lcg(&mut seed) * 8.0) as i64) as f64 - 4.0).collect();
        let k = 1 + (lcg(&mut seed) * (n as f64 - 2.0)) as usize;
        let mut idx: Vec<usize> = (0..n).collect();
        for i in 0..k {
            let j = i + (lcg(&mut seed) * (n - i) as f64) as usize;
            idx.swap(i, j);
        }
        let members = &idx[..k];
        let (r, t) = conjuntos::ranks(&v);
        let z = conjuntos::rank_sum_z(&r, t, members);
        let zp = z_by_pairs(&v, members);
        let ok = match z {
            Some(z) => (z - zp).abs() < 1e-9,
            None => !zp.is_finite(),
        };
        if !ok {
            bad += 1;
        }
    }
    h.check(bad == 0, "rank-sum z equals the pair count in 200 random cases with ties");

    // 2. A DE table: 300 genes; G_UP* strongly up, G_DOWN* down, MIX* half each way.
    let mut tsv = String::from("gene\tname\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n");
    let mut s2 = 11u64;
    for i in 0..300 {
        let (id, name, stat) = if i < 12 {
            (format!("ENSG{i:05}"), format!("G_UP{i}"), 6.0 + lcg(&mut s2))
        } else if i < 20 {
            (format!("ENSG{i:05}"), format!("G_DOWN{i}"), -6.0 - lcg(&mut s2))
        } else if i < 30 {
            (format!("ENSG{i:05}"), format!("MIX{i}"), if i % 2 == 0 { 7.0 } else { -7.0 } + lcg(&mut s2))
        } else {
            (format!("ENSG{i:05}"), format!("BG{i}"), (lcg(&mut s2) - 0.5) * 2.0)
        };
        let padj = if stat.abs() > 4.0 { 0.001 } else { 0.8 };
        tsv.push_str(&format!("{id}\t{name}\t100\t{}\t0.2\t{stat}\t0.01\t{padj}\n", stat / 5.0));
    }
    tsv.push_str("ENSG99999\tALLZERO\t0\tNA\tNA\tNA\tNA\tNA\n");
    let de = p("de.tsv");
    std::fs::write(&de, &tsv).unwrap();
    let gmt = p("sets.gmt");
    std::fs::write(
        &gmt,
        "UP\tgenes up\tg_up0\tG_UP1\tg_up2\tG_UP3\tG_UP4\tG_UP5\n\
         DOWN\tgenes down\tG_DOWN12\tG_DOWN13\tG_DOWN14\tG_DOWN15\n\
         MIXED\tboth ways\tMIX20\tMIX21\tMIX22\tMIX23\tMIX24\tMIX25\n\
         RANDOM\tbackground\tBG40\tBG41\tBG42\tBG43\tBG44\tBG45\tBG46\n\
         TINY\ttoo few\tG_UP0\tNOT_A_GENE\tALLZERO\n",
    )
    .unwrap();
    let opt = conjuntos::SetOptions { statistic: "stat".into(), alpha: 0.1, min_size: 3, require_padj: true };
    let out = p("sets_out.tsv");
    let res = conjuntos::gene_set_test(&de, &gmt, Some(&out), &opt);
    h.check(res.is_ok(), "gene_set_test runs on a GMT");
    let v: serde_json::Value = serde_json::from_str(&res.unwrap_or_else(|e| panic!("{e}"))).unwrap();
    let by: std::collections::HashMap<String, serde_json::Value> =
        v["sets"].as_array().unwrap().iter().map(|s| (s["set"].as_str().unwrap().to_string(), s.clone())).collect();
    h.check(v["genes_ranked"].as_u64() == Some(300), "NA statistics are left out of the ranking (300 of 301)");
    h.check(by["UP"]["found"].as_u64() == Some(6), "symbols matched without case (g_up0 = G_UP0)");
    h.check(by["UP"]["z"].as_f64().is_some_and(|z| z > 3.0) && by["UP"]["fdr"].as_f64().is_some_and(|q| q < 0.01), "UP set: up, FDR < 0.01");
    h.check(by["DOWN"]["z"].as_f64().is_some_and(|z| z < -3.0) && by["DOWN"]["fdr"].as_f64().is_some_and(|q| q < 0.01), "DOWN set: down, FDR < 0.01");
    h.check(
        by["MIXED"]["pvalue"].as_f64().is_some_and(|p| p > 0.05) && by["MIXED"]["fdr_mixed"].as_f64().is_some_and(|q| q < 0.01),
        "MIXED set: not one direction, but changes in either direction",
    );
    h.check(by["RANDOM"]["fdr"].as_f64().is_some_and(|q| q > 0.1), "RANDOM set: not significant");
    h.check(by["UP"]["up"].as_u64() == Some(6) && by["DOWN"]["down"].as_u64() == Some(4), "up / down counts from padj and log2FoldChange");
    h.check(
        by["TINY"]["pvalue"].is_null()
            && by["TINY"]["found"].as_u64() == Some(1)
            && by["TINY"]["missing"].as_array().is_some_and(|m| m.iter().any(|x| x == "NOT_A_GENE") && m.iter().any(|x| x == "ALLZERO (not tested)")),
        "a set with fewer than min_size genes is not tested; the missing and untested genes are listed",
    );
    let table = std::fs::read_to_string(&out).unwrap_or_default();
    h.check(table.starts_with("set\tdescription\tsize\tfound\t") && table.lines().count() == 6, "results table written (header + 5 sets)");

    // 3. A two-column table, a CSV from R (empty first header cell), and errors.
    let two = p("sets.tsv");
    std::fs::write(&two, "set\tgene\nUP\tG_UP0\nUP\tG_UP1\nUP\tG_UP2\nUP\tG_UP3\n").unwrap();
    let v2: serde_json::Value = serde_json::from_str(&conjuntos::gene_set_test(&de, &two, None, &opt).unwrap()).unwrap();
    h.check(v2["sets"][0]["found"].as_u64() == Some(4) && v2["sets"][0]["z"].as_f64().is_some_and(|z| z > 2.0), "two-column (set, gene) table");
    let csv = p("de.csv");
    let mut c = String::from("\"\",\"baseMean\",\"log2FoldChange\",\"lfcSE\",\"stat\",\"pvalue\",\"padj\"\n");
    for line in tsv.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        c.push_str(&format!("\"{}\",{},{},{},{},{},{}\n", f[1], f[2], f[3], f[4], f[5], f[6], f[7]));
    }
    std::fs::write(&csv, c).unwrap();
    let v3: serde_json::Value = serde_json::from_str(&conjuntos::gene_set_test(&csv, &gmt, None, &opt).unwrap()).unwrap();
    h.check(v3["sets"][0]["found"].as_u64() == Some(6), "DESeq2 CSV from R (row names in the first column)");
    let e = conjuntos::gene_set_test(&de, &gmt, None, &conjuntos::SetOptions { statistic: "score".into(), alpha: 0.1, min_size: 3, require_padj: true });
    h.check(e.as_ref().is_err_and(|m| m.contains("no column \"score\"")), "missing statistic column: clear error");
    let e2 = conjuntos::gene_set_test(&de, "no_such_collection", None, &opt);
    h.check(e2.as_ref().is_err_and(|m| m.contains("spinal_cord_injury")), "unknown collection: error names the built-in ones");

    // 3b. Genes without an adjusted p-value (filtered out by DESeq2: barely expressed)
    //     are left out of the ranking by default; require_padj=false keeps them.
    let mut low = tsv.clone();
    for i in 0..50 {
        low.push_str(&format!("ENSL{i:05}\tLOW{i}\t1\t0.01\t2.5\t{}\t0.9\tNA\n", (i as f64 - 25.0) / 1000.0));
    }
    let de_low = p("de_low.tsv");
    std::fs::write(&de_low, &low).unwrap();
    let vy: serde_json::Value = serde_json::from_str(&conjuntos::gene_set_test(&de_low, &gmt, None, &opt).unwrap()).unwrap();
    let vn: serde_json::Value = serde_json::from_str(
        &conjuntos::gene_set_test(&de_low, &gmt, None, &conjuntos::SetOptions { require_padj: false, ..opt.clone() }).unwrap(),
    )
    .unwrap();
    h.check(
        vy["genes_ranked"].as_u64() == Some(300) && vy["genes_without_padj"].as_u64() == Some(50) && vn["genes_ranked"].as_u64() == Some(350),
        "genes without padj: left out by default (300 ranked), kept with require_padj=false (350)",
    );

    // 4. The built-in collection.
    let sci = conjuntos::spinal_cord_injury();
    h.check(sci.len() == 11 && sci.iter().all(|s| s.genes.len() >= 8), "built-in spinal_cord_injury: 11 sets, 8 or more genes each");
    let vb: serde_json::Value = serde_json::from_str(&conjuntos::gene_set_test(&de, "spinal_cord_injury", None, &opt).unwrap()).unwrap();
    h.check(
        vb["gene_sets"].as_str().is_some_and(|s| s.contains("teaching")) && vb["sets"].as_array().map(|a| a.len()) == Some(11),
        "built-in collection runs (no genes found in a synthetic table: nothing tested)",
    );

    let _ = std::fs::remove_dir_all(&dir);
    println!("conjuntos tests: {} passed, {} failed", h.ok, h.fail);
    if h.fail > 0 {
        std::process::exit(1);
    }
}
