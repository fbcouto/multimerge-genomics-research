//! Testes de conjuntos de genes sobre uma tabela de expressão diferencial (a do
//! `differential_expression`, ou a do DESeq2 do R: gene, name, baseMean,
//! log2FoldChange, lfcSE, stat, pvalue, padj).
//!
//! O teste é o de postos de Wilcoxon-Mann-Whitney da estatística dos genes do
//! conjunto contra a dos outros genes testados, como o `wilcoxGST` do limma
//! (teste "competitivo": o conjunto muda mais que o resto?). Dois sentidos:
//! "direcional" (a estatística com sinal: o conjunto sobe ou desce) e "misto"
//! (o valor absoluto: o conjunto muda, para cima ou para baixo, mais que o
//! resto). Como o do limma, supõe genes independentes; genes correlacionados
//! deixam o p menor do que deveria (o CAMERA corrige isso; não é feito aqui).
//!
//! Os conjuntos vêm de um arquivo GMT (nome, descrição, genes...) ou de uma
//! tabela de duas colunas (conjunto, gene), ou de uma coleção embutida
//! (`spinal_cord_injury`): conjuntos de marcadores citados em revisões sobre
//! lesão medular, para ensino e análise exploratória, não assinaturas validadas.
//! Os símbolos são comparados sem caixa (GFAP humano = Gfap de rato e camundongo).

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Write};

use crate::expressao::{benjamini_hochberg, normal_cdf};
use crate::externo;

#[derive(Debug, Clone)]
pub struct GeneSet {
    pub name: String,
    pub description: String,
    pub genes: Vec<String>,
}

/// The built-in collections.
pub const BUILTIN: &[&str] = &["spinal_cord_injury"];

fn set(name: &str, description: &str, genes: &[&str]) -> GeneSet {
    GeneSet { name: name.to_string(), description: description.to_string(), genes: genes.iter().map(|g| g.to_string()).collect() }
}

/// Marker genes often cited in reviews of spinal cord injury and repair, grouped by
/// process. Human symbols; the rodent ones that differ are added (Serpina3n, Chil1,
/// Lyz2). For teaching and exploratory analysis: not validated signatures.
pub fn spinal_cord_injury() -> Vec<GeneSet> {
    vec![
        set("reactive_astrocytes", "Astrocyte reactivity after CNS injury (Escartin et al. 2021, Nat Neurosci)",
            &["GFAP", "VIM", "NES", "CD44", "LCN2", "SERPINA3", "SERPINA3N", "CHI3L1", "CHIL1", "STAT3", "TIMP1", "OSMR"]),
        set("glial_scar_cspgs", "Chondroitin sulfate proteoglycans of the glial scar and their synthesis (Silver & Miller 2004, Nat Rev Neurosci)",
            &["ACAN", "BCAN", "NCAN", "VCAN", "CSPG4", "PTPRZ1", "XYLT1", "CHST11", "CHSY1", "TNC"]),
        set("fibrotic_scar", "Fibrotic scar and extracellular matrix of the lesion core (Soderblom et al. 2013, J Neurosci)",
            &["COL1A1", "COL1A2", "COL3A1", "COL4A1", "FN1", "PDGFRB", "POSTN", "LUM"]),
        set("inflammation_myeloid", "Inflammatory cytokines and microglia/macrophage markers",
            &["IL1B", "TNF", "IL6", "CCL2", "CXCL1", "CD68", "AIF1", "ITGAM", "PTPRC", "CD14", "LYZ", "LYZ2"]),
        set("regeneration_associated", "Regeneration-associated genes of the intrinsic axon growth program (Chandran et al. 2016, Neuron)",
            &["ATF3", "GAP43", "SPRR1A", "SOX11", "JUN", "KLF7", "BASP1", "SMAD1", "STAT3", "TUBA1A"]),
        set("axonal_cytoskeleton", "Neurofilaments and axonal microtubule proteins",
            &["NEFL", "NEFM", "NEFH", "INA", "PRPH", "TUBB3", "TUBA1A", "STMN2", "MAPT"]),
        set("myelin_oligodendrocytes", "Myelin and oligodendrocyte genes",
            &["MBP", "PLP1", "MOG", "MAG", "MOBP", "CNP", "OLIG1", "OLIG2", "SOX10", "MYRF"]),
        set("growth_inhibition_rho", "Myelin-associated growth inhibitors, their receptors and the RhoA/ROCK pathway",
            &["RTN4", "RTN4R", "OMG", "LINGO1", "NGFR", "PTPRS", "RHOA", "ROCK1", "ROCK2"]),
        set("neurotrophins", "Neurotrophins and their receptors",
            &["BDNF", "NTRK2", "NTF3", "NTRK3", "NGF", "NTRK1", "GDNF", "GFRA1"]),
        set("laminin_receptors", "Laminin receptors: integrins a3b1, a6b1, a6b4, a7b1, dystroglycan, Lutheran (Durbeej 2010, Cell Tissue Res)",
            &["ITGA3", "ITGA6", "ITGA7", "ITGB1", "ITGB4", "DAG1", "BCAM", "RPSA"]),
        set("laminin_chains", "Laminin alpha, beta and gamma chains",
            &["LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5", "LAMB1", "LAMB2", "LAMB3", "LAMB4", "LAMC1", "LAMC2", "LAMC3"]),
    ]
}

/// The sets of a built-in collection, or of a GMT file / two-column table.
pub fn load_sets(source: &str) -> Result<(Vec<GeneSet>, String), String> {
    if source == "spinal_cord_injury" {
        return Ok((spinal_cord_injury(), "built-in: spinal_cord_injury (teaching sets)".to_string()));
    }
    if !std::path::Path::new(source).exists() {
        return Err(format!(
            "{source}: no such file, and not a built-in collection ({}). Give a GMT file (name, description, genes, tab-separated) \
             or a two-column table (set, gene).",
            BUILTIN.join(", ")
        ));
    }
    let reader = externo::abre_leitura(source)?;
    let mut sets: Vec<GeneSet> = Vec::new();
    let mut by_name: HashMap<String, usize> = HashMap::new();
    let gmt = source.to_ascii_lowercase().trim_end_matches(".gz").ends_with(".gmt");
    for (i, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| format!("{source}: line {}: {e}", i + 1))?;
        let t = line.trim_end_matches(['\r', '\n']);
        if t.trim().is_empty() || t.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = t.split('\t').map(str::trim).collect();
        if gmt || cols.len() > 2 {
            if cols.len() < 3 {
                return Err(format!("{source}: line {}: a GMT line needs a name, a description and at least one gene.", i + 1));
            }
            let genes: Vec<String> = cols[2..].iter().filter(|g| !g.is_empty()).map(|g| g.to_string()).collect();
            sets.push(GeneSet { name: cols[0].to_string(), description: cols[1].to_string(), genes });
        } else if cols.len() == 2 {
            if i == 0 && (cols[0].eq_ignore_ascii_case("set") || cols[0].eq_ignore_ascii_case("gene_set") || cols[0].eq_ignore_ascii_case("term")) {
                continue; // header
            }
            let k = *by_name.entry(cols[0].to_string()).or_insert_with(|| {
                sets.push(GeneSet { name: cols[0].to_string(), description: String::new(), genes: Vec::new() });
                sets.len() - 1
            });
            sets[k].genes.push(cols[1].to_string());
        } else {
            return Err(format!("{source}: line {}: expected a GMT line or two columns (set, gene).", i + 1));
        }
    }
    if sets.is_empty() {
        return Err(format!("{source}: no gene set read."));
    }
    Ok((sets, source.to_string()))
}

struct DeRow {
    keys: Vec<String>, // upper-case gene id and name
    label: String,
    stat: Option<f64>,
    lfc: Option<f64>,
    padj: Option<f64>,
}

struct DeTable {
    rows: Vec<DeRow>,
    skipped: usize,
    has_padj: bool,
}

fn parse_num(s: &str) -> Option<f64> {
    let t = s.trim().trim_matches('"');
    if t.is_empty() || t.eq_ignore_ascii_case("na") || t.eq_ignore_ascii_case("nan") || t == "." {
        return None;
    }
    t.parse::<f64>().ok().filter(|v| v.is_finite())
}

fn read_results(path: &str, statistic: &str) -> Result<DeTable, String> {
    let reader = externo::abre_leitura(path)?;
    let mut lines = reader.lines();
    let header = loop {
        match lines.next() {
            Some(l) => {
                let l = l.map_err(|e| format!("{path}: {e}"))?;
                if !l.trim().is_empty() && !l.starts_with('#') {
                    break l;
                }
            }
            None => return Err(format!("{path}: empty table.")),
        }
    };
    let sep = if header.contains('\t') { '\t' } else { ',' };
    let cols: Vec<String> = header.trim_end_matches(['\r', '\n']).split(sep).map(|c| c.trim().trim_matches('"').to_string()).collect();
    let find = |names: &[&str]| cols.iter().position(|c| names.iter().any(|n| c.eq_ignore_ascii_case(n)));
    let stat_col = find(&[statistic]).ok_or_else(|| {
        format!("{path}: no column {statistic:?} (columns: {}). A DESeq2-style table has stat; choose statistic=\"log2FoldChange\" otherwise.", cols.join(", "))
    })?;
    // R's write.csv leaves the first header cell empty for the row names (the gene ids).
    let gene_col = find(&["gene", "gene_id", "id", "ensembl", ""]).unwrap_or(0);
    let name_col = find(&["name", "symbol", "gene_name", "gene_symbol"]);
    let lfc_col = find(&["log2FoldChange", "log2_fold_change", "logFC"]);
    let padj_col = find(&["padj", "FDR", "adj.P.Val", "qvalue"]);
    let mut rows = Vec::new();
    let mut skipped = 0usize;
    for (i, line) in lines.enumerate() {
        let line = line.map_err(|e| format!("{path}: line {}: {e}", i + 2))?;
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.trim_end_matches(['\r', '\n']).split(sep).map(|c| c.trim().trim_matches('"')).collect();
        if f.len() <= stat_col.max(gene_col) {
            skipped += 1;
            continue;
        }
        let gene = f[gene_col].to_string();
        let name = name_col.and_then(|c| f.get(c)).map(|s| s.to_string()).filter(|s| !s.is_empty() && s != "NA");
        let mut keys = vec![gene.to_ascii_uppercase()];
        if let Some(n) = &name {
            keys.push(n.to_ascii_uppercase());
        }
        rows.push(DeRow {
            keys,
            label: name.unwrap_or(gene),
            stat: parse_num(f[stat_col]),
            lfc: lfc_col.and_then(|c| f.get(c)).and_then(|s| parse_num(s)),
            padj: padj_col.and_then(|c| f.get(c)).and_then(|s| parse_num(s)),
        });
    }
    Ok(DeTable { rows, skipped, has_padj: padj_col.is_some() })
}

/// Average ranks (1-based) of `v`, and the tie term sum(t³ − t).
pub fn ranks(v: &[f64]) -> (Vec<f64>, f64) {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap_or(std::cmp::Ordering::Equal));
    let mut r = vec![0.0; v.len()];
    let mut ties = 0.0;
    let mut i = 0;
    while i < idx.len() {
        let mut j = i;
        while j + 1 < idx.len() && v[idx[j + 1]] == v[idx[i]] {
            j += 1;
        }
        let avg = (i + j) as f64 / 2.0 + 1.0;
        for k in i..=j {
            r[idx[k]] = avg;
        }
        let t = (j - i + 1) as f64;
        ties += t * t * t - t;
        i = j + 1;
    }
    (r, ties)
}

/// Wilcoxon-Mann-Whitney z of the members against the rest (normal approximation,
/// tie-corrected variance, no continuity correction): positive when the members
/// rank higher.
pub fn rank_sum_z(ranks: &[f64], ties: f64, members: &[usize]) -> Option<f64> {
    let n = ranks.len() as f64;
    let n1 = members.len() as f64;
    let n2 = n - n1;
    if n1 < 1.0 || n2 < 1.0 {
        return None;
    }
    let w: f64 = members.iter().map(|&i| ranks[i]).sum();
    let u = w - n1 * (n1 + 1.0) / 2.0;
    let var = n1 * n2 / 12.0 * ((n + 1.0) - ties / (n * (n - 1.0)));
    (var > 0.0).then(|| (u - n1 * n2 / 2.0) / var.sqrt())
}

#[derive(Debug, Clone, Default)]
pub struct SetOptions {
    pub statistic: String,
    pub alpha: f64,
    pub min_size: usize,
    /// Rank only the genes with an adjusted p-value (when the table has the column):
    /// the ones DESeq2 kept after independent filtering. Without this, thousands of
    /// barely expressed genes with statistics near zero make any set of expressed
    /// genes look as if it changed.
    pub require_padj: bool,
}

struct SetRow {
    name: String,
    description: String,
    size: usize,
    found: Vec<String>,
    missing: Vec<String>,
    mean_lfc: Option<f64>,
    up: usize,
    down: usize,
    z: Option<f64>,
    p: Option<f64>,
    fdr: Option<f64>,
    p_mixed: Option<f64>,
    fdr_mixed: Option<f64>,
}

/// Three significant figures, scientific notation below 0.001.
fn short(x: f64) -> String {
    if x != 0.0 && x.abs() < 0.001 {
        format!("{x:.2e}")
    } else {
        format!("{x:.3}")
    }
}

fn jnum(v: Option<f64>) -> String {
    match v {
        Some(x) if x.is_finite() => format!("{x}"),
        _ => "null".to_string(),
    }
}

fn jstr(s: &str) -> String {
    serde_json::to_string(s).unwrap_or_else(|_| "\"\"".to_string())
}

/// Tests every gene set against a differential-expression table. Writes the
/// results table to `output` when given; returns JSON.
pub fn gene_set_test(results: &str, gene_sets: &str, output: Option<&str>, opt: &SetOptions) -> Result<String, String> {
    let statistic = if opt.statistic.is_empty() { "stat" } else { opt.statistic.as_str() };
    let min_size = opt.min_size.max(1);
    let (sets, source) = load_sets(gene_sets)?;
    let DeTable { rows, skipped, has_padj } = read_results(results, statistic)?;
    // The universe: the genes with a value of the statistic (DESeq2 leaves NA for
    // genes with all counts zero or flagged as outliers) and, by default, with an
    // adjusted p-value (the genes kept by independent filtering).
    let use_padj = opt.require_padj && has_padj;
    let with_stat = rows.iter().filter(|r| r.stat.is_some()).count();
    let tested: Vec<usize> = (0..rows.len()).filter(|&i| rows[i].stat.is_some() && (!use_padj || rows[i].padj.is_some())).collect();
    let filtered_out = with_stat - tested.len();
    if tested.len() < 3 {
        return Err(format!(
            "{results}: only {} gene(s) with a value of {statistic}{}: nothing to rank.",
            tested.len(),
            if use_padj { " and an adjusted p-value (require_padj=False ranks every gene with the statistic)" } else { "" }
        ));
    }
    let pos_of: HashMap<usize, usize> = tested.iter().enumerate().map(|(p, &i)| (i, p)).collect();
    let mut by_key: HashMap<&str, usize> = HashMap::new();
    for (i, r) in rows.iter().enumerate() {
        for k in &r.keys {
            by_key.entry(k.as_str()).or_insert(i);
        }
    }
    let stats: Vec<f64> = tested.iter().map(|&i| rows[i].stat.unwrap()).collect();
    let abs: Vec<f64> = stats.iter().map(|s| s.abs()).collect();
    let (r_dir, t_dir) = ranks(&stats);
    let (r_abs, t_abs) = ranks(&abs);

    let mut out: Vec<SetRow> = Vec::new();
    for s in &sets {
        let mut seen: HashSet<usize> = HashSet::new();
        let mut found = Vec::new();
        let mut missing = Vec::new();
        let mut members: Vec<usize> = Vec::new();
        let mut uniq: HashSet<String> = HashSet::new();
        for g in &s.genes {
            if !uniq.insert(g.to_ascii_uppercase()) {
                continue;
            }
            match by_key.get(g.to_ascii_uppercase().as_str()) {
                Some(&i) if seen.insert(i) => match pos_of.get(&i) {
                    Some(&p) => {
                        members.push(p);
                        found.push(rows[i].label.clone());
                    }
                    None => missing.push(format!("{g} (not tested)")),
                },
                Some(_) => {}
                None => missing.push(g.clone()),
            }
        }
        let lfcs: Vec<f64> = members.iter().filter_map(|&p| rows[tested[p]].lfc).collect();
        let mean_lfc = (!lfcs.is_empty()).then(|| lfcs.iter().sum::<f64>() / lfcs.len() as f64);
        let sig = |p: usize, up: bool| {
            let r = &rows[tested[p]];
            r.padj.is_some_and(|q| q < opt.alpha) && r.lfc.is_some_and(|l| if up { l > 0.0 } else { l < 0.0 })
        };
        let up = members.iter().filter(|&&p| sig(p, true)).count();
        let down = members.iter().filter(|&&p| sig(p, false)).count();
        let (z, p, p_mixed) = if members.len() >= min_size {
            let z = rank_sum_z(&r_dir, t_dir, &members);
            let zm = rank_sum_z(&r_abs, t_abs, &members);
            (z, z.map(|z| 2.0 * (1.0 - normal_cdf(z.abs())).max(0.0)), zm.map(|z| (1.0 - normal_cdf(z)).max(0.0)))
        } else {
            (None, None, None)
        };
        out.push(SetRow {
            name: s.name.clone(),
            description: s.description.clone(),
            size: uniq.len(),
            found,
            missing,
            mean_lfc,
            up,
            down,
            z,
            p,
            fdr: None,
            p_mixed,
            fdr_mixed: None,
        });
    }
    // Benjamini-Hochberg over the sets that were tested, for each kind of test.
    for mixed in [false, true] {
        let idx: Vec<usize> = (0..out.len()).filter(|&i| if mixed { out[i].p_mixed.is_some() } else { out[i].p.is_some() }).collect();
        let ps: Vec<f64> = idx.iter().map(|&i| if mixed { out[i].p_mixed.unwrap() } else { out[i].p.unwrap() }).collect();
        let q = benjamini_hochberg(&ps);
        for (k, &i) in idx.iter().enumerate() {
            if mixed {
                out[i].fdr_mixed = Some(q[k]);
            } else {
                out[i].fdr = Some(q[k]);
            }
        }
    }

    if let Some(path) = output {
        let mut w = externo::Escritor::cria(path)?;
        let mut text = String::from("set\tdescription\tsize\tfound\tmean_log2FoldChange\tup\tdown\tz\tpvalue\tfdr\tdirection\tpvalue_mixed\tfdr_mixed\tgenes\n");
        for r in &out {
            let dir = match (r.z, r.fdr) {
                (Some(z), Some(q)) if q < opt.alpha => if z > 0.0 { "up" } else { "down" },
                (Some(_), _) => "-",
                _ => "too few genes",
            };
            let f = |v: Option<f64>| v.map(|x| format!("{x:.6}")).unwrap_or_else(|| "NA".to_string());
            text.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                r.name, r.description, r.size, r.found.len(), f(r.mean_lfc), r.up, r.down, f(r.z), f(r.p), f(r.fdr), dir,
                f(r.p_mixed), f(r.fdr_mixed), r.found.join(",")
            ));
        }
        w.write_all(text.as_bytes()).map_err(|e| format!("{path}: {e}"))?;
        w.termina()?;
    }

    let n_sig = out.iter().filter(|r| r.fdr.is_some_and(|q| q < opt.alpha)).count();
    let n_mixed = out.iter().filter(|r| r.fdr_mixed.is_some_and(|q| q < opt.alpha)).count();
    let n_tested = out.iter().filter(|r| r.p.is_some()).count();
    let mut lines = vec![format!(
        "{} gene sets ({source}) against {} genes ranked by {statistic} in {results}{}: {n_tested} with at least {min_size} genes found; \
         {n_sig} change in one direction and {n_mixed} change in either direction at FDR < {}.",
        out.len(),
        tested.len(),
        if filtered_out > 0 { format!(" ({filtered_out} without an adjusted p-value left out)") } else { String::new() },
        opt.alpha
    )];
    let mut order: Vec<usize> = (0..out.len()).filter(|&i| out[i].p.is_some()).collect();
    order.sort_by(|&a, &b| out[a].p.partial_cmp(&out[b].p).unwrap_or(std::cmp::Ordering::Equal));
    for &i in order.iter().take(8) {
        let r = &out[i];
        lines.push(format!(
            "  {}: {} of {} genes, z {:.2}, FDR {} ({}; mean log2FC {})",
            r.name,
            r.found.len(),
            r.size,
            r.z.unwrap_or(0.0),
            short(r.fdr.unwrap_or(1.0)),
            if r.z.unwrap_or(0.0) > 0.0 { "up" } else { "down" },
            r.mean_lfc.map(|v| format!("{v:.2}")).unwrap_or_else(|| "NA".to_string())
        ));
    }
    let sets_json: Vec<String> = out
        .iter()
        .map(|r| {
            format!(
                "{{\"set\":{},\"description\":{},\"size\":{},\"found\":{},\"missing\":{},\"genes\":{},\"mean_log2_fold_change\":{},\"up\":{},\"down\":{},\"z\":{},\"pvalue\":{},\"fdr\":{},\"pvalue_mixed\":{},\"fdr_mixed\":{}}}",
                jstr(&r.name),
                jstr(&r.description),
                r.size,
                r.found.len(),
                serde_json::to_string(&r.missing).unwrap_or_else(|_| "[]".into()),
                serde_json::to_string(&r.found).unwrap_or_else(|_| "[]".into()),
                jnum(r.mean_lfc),
                r.up,
                r.down,
                jnum(r.z),
                jnum(r.p),
                jnum(r.fdr),
                jnum(r.p_mixed),
                jnum(r.fdr_mixed)
            )
        })
        .collect();
    Ok(format!(
        "{{\"results\":{},\"gene_sets\":{},\"statistic\":{},\"alpha\":{},\"min_size\":{},\"genes_ranked\":{},\"genes_without_padj\":{},\"rows_skipped\":{},\"output\":{},\"sets\":[{}],\"significant\":{},\"significant_mixed\":{},\"text\":{}}}",
        jstr(results),
        jstr(&source),
        jstr(statistic),
        opt.alpha,
        min_size,
        tested.len(),
        filtered_out,
        skipped,
        output.map(jstr).unwrap_or_else(|| "null".into()),
        sets_json.join(","),
        n_sig,
        n_mixed,
        jstr(&lines.join("\n"))
    ))
}
