













//! BAM+VCF structural/pileup calling for `CYP2D6` and pure-Rust Class I `HLA`
//! (`HLA-A`, `HLA-B`, `HLA-C`) genotyping for `mmgen-pgx` v2.
//!
//! References & licensing:
//! - IPD-IMGT/HLA genomic sequences (`{A,B,C}_gen.fasta`, CC BY-ND 4.0) are
//!   NEVER bundled inside the repository or release tables; they are downloaded
//!   on demand (or resolved from `MMGEN_HLA_DIR`) and cached locally unmodified.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::process::Command;

use serde::Deserialize;

use crate::{bam, bgzf, indice, referencia};

use super::{
    call_combinatorial_gene, format_activity_score, lookup_phenotype_by_counts,
    observe_gene_positions, ref_cpic_value, star_sort_key, Bundle, EvidenceEntry, GeneDef,
    GeneOutcome, GenePhenotypeTable, HapAllele, PosObservation, VcfObservations,
};

// ---------------------------------------------------------------------------
// BAM Genes Manifest Structs (`bam_genes.json`)
// ---------------------------------------------------------------------------

#[derive(Debug, Deserialize)]
pub struct BamGenesBundle {
    pub schema_version: u32,
    pub assembly: String,
    pub data_version: String,
    pub cyp2d6: Cyp2d6BamDef,
    pub hla_genes: Vec<HlaBamDef>,
}

#[derive(Debug, Deserialize)]
pub struct Cyp2d6BamDef {
    #[serde(flatten)]
    pub base: GeneDef,
    #[serde(default)]
    pub sv_alleles: Vec<SvAlleleDef>,
    #[serde(default)]
    pub duplication_alleles: BTreeMap<String, DupAlleleDef>,
    #[serde(default)]
    pub canonical_iupac_defaults: BTreeMap<String, String>,
    pub cnv_regions: CnvRegionsDef,
}

#[derive(Debug, Deserialize)]
pub struct SvAlleleDef {
    pub name: String,
    pub sv_type: String,
    pub clinical_function: Option<String>,
    pub activity_value: Option<f64>,
    pub activity_label: Option<String>,
}

#[derive(Debug, Deserialize)]
pub struct DupAlleleDef {
    pub name: String,
    pub clinical_function: Option<String>,
    pub activity_value: Option<f64>,
    pub activity_label: Option<String>,
}

#[derive(Debug, Deserialize)]
pub struct CnvRegionsDef {
    pub target_whole: GenomicSpanDef,
    pub target_5p: GenomicSpanDef,
    #[serde(rename = "target_ex9")]
    pub target_ex9: GenomicSpanDef,
    pub paralog: GenomicSpanDef,
    pub controls: Vec<NamedSpanDef>,
}

#[derive(Debug, Deserialize, Clone)]
pub struct GenomicSpanDef {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

#[derive(Debug, Deserialize, Clone)]
pub struct NamedSpanDef {
    pub name: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

#[derive(Debug, Deserialize, Clone)]
pub struct HlaBamDef {
    pub gene: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub lookup_method: String,
    pub fasta_filename: String,
    pub fasta_url: String,
    #[serde(default)]
    pub risk_alleles: Vec<String>,
    /// IPD-IMGT/HLA release the tables were validated with (e.g. "3.66.0"). With it, only that
    /// release's file is used: it is looked for in `<cache>/<version>/`, then in `<cache>/`,
    /// and downloaded to `<cache>/<version>/` from `fasta_url` when absent.
    #[serde(default)]
    pub imgt_hla_version: Option<String>,
    /// SHA-256 and size of that release's unmodified `fasta_filename`.
    #[serde(default)]
    pub fasta_sha256: Option<String>,
    #[serde(default)]
    pub fasta_bytes: Option<u64>,
}

// ---------------------------------------------------------------------------
// GRCh38 BAM Header Validation & Contig Mapping
// ---------------------------------------------------------------------------

const GRCH38_LENGTHS: &[(&str, u32)] = &[
    ("1", 248_956_422),
    ("4", 190_214_555),
    ("6", 170_805_979),
    ("7", 159_345_973),
    ("10", 133_797_422),
    ("12", 133_275_309),
    ("14", 107_043_718),
    ("22", 50_818_468),
];

fn canonical_nochr(s: &str) -> &str {
    s.strip_prefix("chr").or_else(|| s.strip_prefix("CHR")).unwrap_or(s)
}

fn resolve_bam_contig(ref_names: &HashMap<String, usize>, refs: &[(String, u32)], wanted: &str) -> Option<String> {
    if let Some(&i) = ref_names.get(wanted) {
        return Some(refs[i].0.clone());
    }
    let alt = if let Some(bare) = wanted.strip_prefix("chr") {
        bare.to_string()
    } else {
        format!("chr{wanted}")
    };
    ref_names.get(&alt).map(|&i| refs[i].0.clone())
}

fn validate_bam_header(bam_path: &str) -> Result<(Vec<(String, u32)>, HashMap<String, usize>), String> {
    let f = File::open(bam_path)
        .map_err(|e| format!("{bam_path}: cannot open BAM file ({e}). Provide a readable coordinate-sorted GRCh38 BAM with .bai index."))?;
    let mut bgzf = bgzf::BgzfReader::novo(BufReader::new(f));
    let cab = bam::le_cabecalho(&mut bgzf)
        .map_err(|e| format!("{bam_path}: invalid or corrupted BAM header ({e})."))?;

    // Reject queryname-sorted BAMs; real-world indexed BAMs (e.g. GIAB Novoalign)
    // can carry SO:unsorted in @HD even after coordinate sorting and .bai indexing.
    for line in cab.texto.lines() {
        if line.starts_with("@HD") {
            for field in line.split('\t') {
                if let Some(so) = field.strip_prefix("SO:") {
                    if so == "queryname" {
                        return Err(format!(
                            "{bam_path}: BAM header reports sort order SO:{so}; mmgen-pgx requires a coordinate-sorted BAM with .bai index."
                        ));
                    }
                }
            }
        }
    }

    let mut ref_map = HashMap::with_capacity(cab.refs.len());
    for (i, (name, len)) in cab.refs.iter().enumerate() {
        ref_map.insert(name.clone(), i);
        let bare = canonical_nochr(name);
        if let Some((_, expected_len)) = GRCH38_LENGTHS.iter().find(|&&(k, _)| k == bare) {
            if len != expected_len {
                return Err(format!(
                    "{bam_path}: BAM reference assembly is not GRCh38 (contig {name} has length {len}, expected GRCh38 length {expected_len})."
                ));
            }
        }
    }

    Ok((cab.refs, ref_map))
}

// ---------------------------------------------------------------------------
// Single-Pass .bai Region Query Collector (Depth, CYP2D6 Pileup, HLA Reads)
// ---------------------------------------------------------------------------

#[derive(Default, Clone)]
struct SitePileup {
    counts: [u32; 4], // A, C, G, T
    clean_depth: u32,
    indel_depth: u32,
}

impl SitePileup {
    fn base_idx(b: u8) -> Option<usize> {
        match b.to_ascii_uppercase() {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        }
    }

    fn total_bases(&self) -> u32 {
        self.counts.iter().sum()
    }
}

struct BamEvidencePool {
    ctrl_names: Vec<String>,
    ctrl_bases: Vec<u64>,
    ctrl_lens: Vec<u64>,
    whole_bases: u64,
    whole_len: u64,
    fivep_bases: u64,
    fivep_len: u64,
    ex9_bases: u64,
    ex9_len: u64,
    paralog_bases: u64,
    paralog_len: u64,
    /// Per CYP2D6 position index (1-based `p.pos` matching SNV base) pileup
    pileup: Vec<SitePileup>,
    /// Fragment observations per qname: vec of (pos_idx, base)
    fragments: HashMap<String, Vec<(usize, u8)>>,
    /// Per HLA gene (0=HLA-A, 1=HLA-B, 2=HLA-C): qname -> Vec<(mate_1_or_2, Vec<u8>)>
    hla_reads: Vec<HashMap<String, Vec<(u8, Vec<u8>)>>>,
}

fn parse_cigar_ops(cigar: &str, out: &mut Vec<(u64, u8)>) -> u64 {
    out.clear();
    let mut num: u64 = 0;
    let mut ref_span: u64 = 0;
    for &b in cigar.as_bytes() {
        if b.is_ascii_digit() {
            num = num * 10 + u64::from(b - b'0');
        } else {
            out.push((num, b));
            if matches!(b, b'M' | b'D' | b'N' | b'=' | b'X') {
                ref_span += num;
            }
            num = 0;
        }
    }
    ref_span
}

fn overlap_len(a0: u64, a1: u64, b0: u64, b1: u64) -> u64 {
    let lo = a0.max(b0);
    let hi = a1.min(b1);
    if hi >= lo { hi - lo + 1 } else { 0 }
}

fn collect_bam_evidence(bam_path: &str, def: &BamGenesBundle) -> Result<BamEvidencePool, String> {
    let (refs, ref_map) = validate_bam_header(bam_path)?;

    let cyp = &def.cyp2d6;
    let cyp_contig = resolve_bam_contig(&ref_map, &refs, &cyp.base.chrom).ok_or_else(|| {
        format!(
            "{bam_path}: BAM header does not contain {} (required for CYP2D6 analysis).",
            cyp.base.chrom
        )
    })?;
    let paralog_contig = resolve_bam_contig(&ref_map, &refs, &cyp.cnv_regions.paralog.chrom)
        .unwrap_or_else(|| cyp_contig.clone());

    let mut regions: Vec<String> = Vec::new();
    let mut ctrl_spans: Vec<(String, String, u64, u64)> = Vec::new();
    for ctrl in &cyp.cnv_regions.controls {
        if let Some(rc) = resolve_bam_contig(&ref_map, &refs, &ctrl.chrom) {
            regions.push(format!("{rc}:{}-{}", ctrl.start, ctrl.end));
            ctrl_spans.push((ctrl.name.clone(), rc, ctrl.start, ctrl.end));
        }
    }
    if ctrl_spans.is_empty() {
        return Err(format!(
            "{bam_path}: BAM contains none of the diploid depth control contigs (e.g. chr4, chr22)."
        ));
    }

    let whole = &cyp.cnv_regions.target_whole;
    let fivep = &cyp.cnv_regions.target_5p;
    let ex9 = &cyp.cnv_regions.target_ex9;
    let paralog = &cyp.cnv_regions.paralog;
    regions.push(format!("{cyp_contig}:{}-{}", whole.start, whole.end));
    regions.push(format!("{paralog_contig}:{}-{}", paralog.start, paralog.end));

    let mut hla_spans: Vec<(String, u64, u64)> = Vec::with_capacity(def.hla_genes.len());
    for hg in &def.hla_genes {
        let hc = resolve_bam_contig(&ref_map, &refs, &hg.chrom).ok_or_else(|| {
            format!(
                "{bam_path}: BAM header does not contain {} (required for {} analysis).",
                hg.chrom, hg.gene
            )
        })?;
        regions.push(format!("{hc}:{}-{}", hg.start, hg.end));
        hla_spans.push((hc, hg.start, hg.end));
    }

    // Build sorted list of CYP2D6 target genomic positions for fast binary-search/slice pileup
    let mut cyp_pos_sorted: Vec<(u64, usize)> = cyp
        .base
        .positions
        .iter()
        .enumerate()
        .map(|(idx, p)| (p.pos, idx))
        .collect();
    cyp_pos_sorted.sort_unstable();

    let mut ctrl_bases: Vec<u64> = vec![0; ctrl_spans.len()];
    let mut whole_bases = 0u64;
    let mut fivep_bases = 0u64;
    let mut ex9_bases = 0u64;
    let mut paralog_bases = 0u64;

    let mut pileup = vec![SitePileup::default(); cyp.base.positions.len()];
    let mut fragments: HashMap<String, Vec<(usize, u8)>> = HashMap::new();
    let mut hla_reads: Vec<HashMap<String, Vec<(u8, Vec<u8>)>>> = vec![HashMap::new(); def.hla_genes.len()];
    let mut cigar_buf: Vec<(u64, u8)> = Vec::with_capacity(16);

    indice::consultar(bam_path, &regions, &mut |is_hdr, line| {
        if is_hdr || line.is_empty() {
            return Ok(());
        }
        let mut fields = line.split('\t');
        let Some(qname) = fields.next() else { return Ok(()); };
        let Some(flag_s) = fields.next() else { return Ok(()); };
        let Some(rname) = fields.next() else { return Ok(()); };
        let Some(pos_s) = fields.next() else { return Ok(()); };
        let Some(mapq_s) = fields.next() else { return Ok(()); };
        let Some(cigar_s) = fields.next() else { return Ok(()); };
        let _ = fields.next(); // RNEXT
        let _ = fields.next(); // PNEXT
        let _ = fields.next(); // TLEN
        let Some(seq_s) = fields.next() else { return Ok(()); };
        let Some(qual_s) = fields.next() else { return Ok(()); };

        let flag: u16 = flag_s.parse().unwrap_or(0);
        // Reject unmapped, secondary, QC-fail, PCR duplicate, supplementary
        if (flag & 0xF04) != 0 {
            return Ok(());
        }
        let Ok(rpos) = pos_s.parse::<u64>() else { return Ok(()); };
        let mapq: u8 = mapq_s.parse().unwrap_or(0);
        let ref_span = parse_cigar_ops(cigar_s, &mut cigar_buf);
        if ref_span == 0 {
            return Ok(());
        }
        let rend = rpos + ref_span - 1;

        // Collect HLA reads regardless of MAPQ (paralogy often lowers BWA MAPQ at HLA exons)
        for (h_idx, (hc, hs, he)) in hla_spans.iter().enumerate() {
            if rname == hc && rend >= *hs && rpos <= *he && seq_s != "*" && seq_s.len() >= 50 {
                let mate = if (flag & 0x80) != 0 { 2 } else { 1 };
                hla_reads[h_idx]
                    .entry(qname.to_string())
                    .or_default()
                    .push((mate, seq_s.as_bytes().to_ascii_uppercase()));
            }
        }

        // Depth controls + CYP2D6 CNV/pileup require MAPQ >= 10
        if mapq < 10 {
            return Ok(());
        }

        // Control depth overlap
        let has_ctrl_hit = ctrl_spans
            .iter()
            .enumerate()
            .any(|(i, (_, rc, cs, ce))| {
                let _ = i;
                rname == rc && rend >= *cs && rpos <= *ce
            });
        let is_cyp_whole = rname == cyp_contig && rend >= whole.start && rpos <= whole.end;
        let is_paralog = rname == paralog_contig && rend >= paralog.start && rpos <= paralog.end;

        if has_ctrl_hit || is_cyp_whole || is_paralog {
            let mut cur_r = rpos;
            for &(len, op) in &cigar_buf {
                match op {
                    b'M' | b'=' | b'X' => {
                        let op_end = cur_r + len - 1;
                        if has_ctrl_hit {
                            for (i, (_, rc, cs, ce)) in ctrl_spans.iter().enumerate() {
                                if rname == rc {
                                    ctrl_bases[i] += overlap_len(cur_r, op_end, *cs, *ce);
                                }
                            }
                        }
                        if is_cyp_whole {
                            whole_bases += overlap_len(cur_r, op_end, whole.start, whole.end);
                            fivep_bases += overlap_len(cur_r, op_end, fivep.start, fivep.end);
                            ex9_bases += overlap_len(cur_r, op_end, ex9.start, ex9.end);
                        }
                        if is_paralog {
                            paralog_bases += overlap_len(cur_r, op_end, paralog.start, paralog.end);
                        }
                        cur_r += len;
                    }
                    b'D' | b'N' => cur_r += len,
                    _ => {}
                }
            }
        }

        // CYP2D6 base pileup + fragment phasing across target positions
        if is_cyp_whole && seq_s != "*" {
            let seq_bytes = seq_s.as_bytes();
            let qual_bytes = if qual_s != "*" { qual_s.as_bytes() } else { &[][..] };
            let mut rp = rpos;
            let mut qp = 0usize;
            let mut frag_pts: Vec<(usize, u8)> = Vec::new();
            let lo_idx = cyp_pos_sorted.partition_point(|&(p, _)| p < rpos);
            let hi_idx = cyp_pos_sorted.partition_point(|&(p, _)| p <= rend);
            if lo_idx < hi_idx {
                let slice = &cyp_pos_sorted[lo_idx..hi_idx];
                for &(len, op) in &cigar_buf {
                    match op {
                        b'M' | b'=' | b'X' => {
                            let op_end = rp + len;
                            let sub_lo = slice.partition_point(|&(p, _)| p < rp);
                            let sub_hi = slice.partition_point(|&(p, _)| p < op_end);
                            for &(p, orig_idx) in &slice[sub_lo..sub_hi] {
                                let qi = qp + (p - rp) as usize;
                                if qi < seq_bytes.len() {
                                    let bq = qual_bytes.get(qi).map(|&q| q.saturating_sub(33)).unwrap_or(30);
                                    if bq >= 15 {
                                        let base = seq_bytes[qi].to_ascii_uppercase();
                                        if let Some(bi) = SitePileup::base_idx(base) {
                                            pileup[orig_idx].counts[bi] += 1;
                                            pileup[orig_idx].clean_depth += 1;
                                            frag_pts.push((orig_idx, base));
                                        }
                                    }
                                }
                            }
                            rp += len;
                            qp += len as usize;
                        }
                        b'I' => {
                            for &(p, orig_idx) in slice {
                                if p + 1 == rp || p == rp {
                                    pileup[orig_idx].indel_depth += 1;
                                }
                            }
                            qp += len as usize;
                        }
                        b'S' => qp += len as usize,
                        b'D' => {
                            let op_end = rp + len;
                            for &(p, orig_idx) in slice {
                                if p + 1 >= rp && p < op_end {
                                    pileup[orig_idx].indel_depth += 1;
                                }
                            }
                            rp += len;
                        }
                        b'N' => rp += len,
                        _ => {}
                    }
                }
            }
            if !frag_pts.is_empty() {
                fragments.entry(qname.to_string()).or_default().extend(frag_pts);
            }
        }
        Ok(())
    })?;

    Ok(BamEvidencePool {
        ctrl_names: ctrl_spans.iter().map(|c| c.0.clone()).collect(),
        ctrl_bases,
        ctrl_lens: ctrl_spans.iter().map(|c| c.3 - c.2 + 1).collect(),
        whole_bases,
        whole_len: whole.end - whole.start + 1,
        fivep_bases,
        fivep_len: fivep.end - fivep.start + 1,
        ex9_bases,
        ex9_len: ex9.end - ex9.start + 1,
        paralog_bases,
        paralog_len: paralog.end - paralog.start + 1,
        pileup,
        fragments,
        hla_reads,
    })
}

// ---------------------------------------------------------------------------
// CYP2D6 BAM + VCF Hybrid/CNV/Diplotype Caller
// ---------------------------------------------------------------------------

fn median_f64(vals: &mut [f64]) -> f64 {
    if vals.is_empty() {
        return 0.0;
    }
    vals.sort_by(|a, b| a.total_cmp(b));
    let n = vals.len();
    if n % 2 == 1 {
        vals[n / 2]
    } else {
        0.5 * (vals[n / 2 - 1] + vals[n / 2])
    }
}

fn merge_bam_pileup_with_vcf(
    gene: &GeneDef,
    obs_list: &mut [PosObservation],
    evidence: &mut [EvidenceEntry],
    pileup: &[SitePileup],
) {
    for (idx, p) in gene.positions.iter().enumerate() {
        if obs_list[idx].origin != "missing" {
            continue;
        }
        let sp = &pileup[idx];
        let is_snv_only = p.grch38_ref.len() == 1
            && p.alt_variants.iter().all(|av| av.vcf_ref.len() == 1 && av.vcf_alt.len() == 1);

        if is_snv_only {
            let dp = sp.total_bases();
            if dp < 6 {
                continue;
            }
            let bases = ["A", "C", "G", "T"];
            let mut active: Vec<(&str, u32)> = bases
                .iter()
                .enumerate()
                .filter_map(|(bi, &b)| {
                    let cnt = sp.counts[bi];
                    if cnt >= 2 && (cnt as f64 / dp as f64) >= 0.18 {
                        Some((b, cnt))
                    } else {
                        None
                    }
                })
                .collect();
            active.sort_by(|a, b| b.1.cmp(&a.1));
            if active.is_empty() {
                continue;
            }
            let rbase = ref_cpic_value(p);
            let (b0, b1, gt) = if active.len() == 1 || (active[0].1 as f64 / dp as f64) >= 0.85 {
                let b = active[0].0;
                let g = if b == rbase { "0/0" } else { "1/1" };
                (b, b, g.to_string())
            } else {
                let a = active[0].0;
                let b = active[1].0;
                let g = if a == rbase || b == rbase { "0/1" } else { "1/2" };
                (a, b, g.to_string())
            };
            let to_hap = |b: &str| {
                let known = p.alt_variants.iter().any(|av| av.cpic_value == b);
                if known || b == rbase {
                    HapAllele::Cpic(b.to_string())
                } else {
                    HapAllele::NonCuratedAlt(b.to_string())
                }
            };
            obs_list[idx].hap_alleles = vec![to_hap(b0), to_hap(b1)];
            obs_list[idx].origin = "bam";
            obs_list[idx].evidence_gt = gt.clone();
            evidence[idx].origin = "bam".to_string();
            evidence[idx].genotype = gt;
        } else if sp.clean_depth >= 6 && sp.indel_depth == 0 {
            let rbase = ref_cpic_value(p).to_string();
            obs_list[idx].hap_alleles = vec![HapAllele::Cpic(rbase.clone()), HapAllele::Cpic(rbase)];
            obs_list[idx].origin = "bam";
            obs_list[idx].evidence_gt = "0/0".to_string();
            evidence[idx].origin = "bam".to_string();
            evidence[idx].genotype = "0/0".to_string();
        }
    }
}

fn cyp_accepted_set(gene: &GeneDef, allele: &super::AlleleDef, pos_id: usize) -> Vec<String> {
    if let Some(av) = allele.variants.iter().find(|v| v.pos_id == pos_id) {
        av.accepted.clone()
    } else {
        vec![ref_cpic_value(&gene.positions[pos_id]).to_string()]
    }
}

fn allele_matches_base(accepted: &[String], ha: &HapAllele) -> bool {
    match ha {
        HapAllele::Cpic(s) => accepted.iter().any(|x| x == s),
        HapAllele::NonCuratedAlt(_) | HapAllele::SpanningDel | HapAllele::Missing => false,
    }
}

/// Extract reliable cis fragment phase links between heterozygous sites `(pos_i, pos_j)`:
/// returns list of `((i, base_i0, base_i1), (j, base_j0, base_j1))` meaning `(base_i0, base_j0)`
/// and `(base_i1, base_j1)` lie on opposite haplotypes in cis.
fn extract_bam_phase_links(
    gene: &GeneDef,
    obs: &[PosObservation],
    fragments: &HashMap<String, Vec<(usize, u8)>>,
) -> Vec<((usize, String, String), (usize, String, String))> {
    let mut het_set: HashMap<usize, (String, String)> = HashMap::new();
    for (idx, o) in obs.iter().enumerate() {
        if o.origin != "missing" && o.hap_alleles.len() == 2 && o.hap_alleles[0] != o.hap_alleles[1] {
            if let (HapAllele::Cpic(a), HapAllele::Cpic(b)) = (&o.hap_alleles[0], &o.hap_alleles[1]) {
                if a.len() == 1 && b.len() == 1 {
                    let _ = gene.positions[idx].pos;
                    het_set.insert(idx, (a.clone(), b.clone()));
                }
            }
        }
    }
    if het_set.len() < 2 {
        return Vec::new();
    }

    // Count (i, j) -> [same_00_11, swap_01_10]
    let mut pair_counts: HashMap<(usize, usize), [u32; 2]> = HashMap::new();
    for pts in fragments.values() {
        // Deduplicate multiple observations of same pos inside one pair
        let mut by_pos: BTreeMap<usize, u8> = BTreeMap::new();
        for &(pidx, base) in pts {
            if let Some((a, b)) = het_set.get(&pidx) {
                let base_arr = [base];
                let b_str = std::str::from_utf8(&base_arr).unwrap_or("");
                if b_str == a || b_str == b {
                    match by_pos.get(&pidx) {
                        None => {
                            by_pos.insert(pidx, base);
                        }
                        Some(&prev) if prev != base => {
                            by_pos.insert(pidx, b'N'); // discordant within fragment
                        }
                        _ => {}
                    }
                }
            }
        }
        let clean: Vec<(usize, u8)> = by_pos.into_iter().filter(|&(_, b)| b != b'N').collect();
        for u in 0..clean.len() {
            for v in (u + 1)..clean.len() {
                let (pi, bi) = clean[u];
                let (pj, bj) = clean[v];
                let (ai0, ai1) = &het_set[&pi];
                let (aj0, aj1) = &het_set[&pj];
                let bi_arr = [bi];
                let bj_arr = [bj];
                let si = std::str::from_utf8(&bi_arr).unwrap_or("");
                let sj = std::str::from_utf8(&bj_arr).unwrap_or("");
                let entry = pair_counts.entry((pi, pj)).or_default();
                if (si == ai0 && sj == aj0) || (si == ai1 && sj == aj1) {
                    entry[0] += 1;
                } else if (si == ai0 && sj == aj1) || (si == ai1 && sj == aj0) {
                    entry[1] += 1;
                }
            }
        }
    }

    let mut out = Vec::new();
    for ((pi, pj), [cis00, cis01]) in pair_counts {
        let tot = cis00 + cis01;
        if tot < 2 {
            continue;
        }
        let (ai0, ai1) = &het_set[&pi];
        let (aj0, aj1) = &het_set[&pj];
        if cis00 >= 2 && (cis00 as f64 / tot as f64) >= 0.85 {
            out.push(((pi, ai0.clone(), ai1.clone()), (pj, aj0.clone(), aj1.clone())));
        } else if cis01 >= 2 && (cis01 as f64 / tot as f64) >= 0.85 {
            out.push(((pi, ai0.clone(), ai1.clone()), (pj, aj1.clone(), aj0.clone())));
        }
    }
    out
}

/// Resolve the core diploid pair `(a1, a2)` for `CYP2D6` using VCF+BAM calls,
/// BAM fragment phase links, VCF `|` phase sets, and canonical IUPAC suballele specificity.
fn resolve_cyp2d6_diploid<'a>(
    cyp: &'a Cyp2d6BamDef,
    obs: &[PosObservation],
    fragments: &HashMap<String, Vec<(usize, u8)>>,
) -> Option<(&'a super::AlleleDef, &'a super::AlleleDef, bool)> {
    let gene = &cyp.base;
    let ref_allele = gene.alleles.iter().find(|a| a.is_reference)?;

    // Candidate single alleles compatible with at least one allele at each observed position
    let mut cands: Vec<&super::AlleleDef> = Vec::new();
    for al in &gene.alleles {
        let mut ok = true;
        for (pid, o) in obs.iter().enumerate() {
            if o.origin == "missing" || o.hap_alleles.is_empty() {
                continue;
            }
            let acc = cyp_accepted_set(gene, al, pid);
            if !o.hap_alleles.iter().any(|ha| allele_matches_base(&acc, ha)) {
                ok = false;
                break;
            }
        }
        if ok {
            cands.push(al);
        }
    }

    let mut valid_pairs: Vec<(&super::AlleleDef, &super::AlleleDef)> = Vec::new();
    for i in 0..cands.len() {
        for j in i..cands.len() {
            let a1 = cands[i];
            let a2 = cands[j];
            let mut ok = true;
            for (pid, o) in obs.iter().enumerate() {
                if o.origin == "missing" || o.hap_alleles.len() != 2 {
                    continue;
                }
                let acc1 = cyp_accepted_set(gene, a1, pid);
                let acc2 = cyp_accepted_set(gene, a2, pid);
                let h0 = &o.hap_alleles[0];
                let h1 = &o.hap_alleles[1];
                let m01 = allele_matches_base(&acc1, h0) && allele_matches_base(&acc2, h1);
                let m10 = allele_matches_base(&acc1, h1) && allele_matches_base(&acc2, h0);
                if !m01 && !m10 {
                    ok = false;
                    break;
                }
            }
            if ok {
                valid_pairs.push((a1, a2));
            }
        }
    }
    if valid_pairs.is_empty() {
        return None;
    }

    // Apply VCF phase_set phasing + BAM fragment phasing
    let bam_links = extract_bam_phase_links(gene, obs, fragments);
    let mut vcf_ps_groups: HashMap<String, Vec<(usize, String, String)>> = HashMap::new();
    for (pid, o) in obs.iter().enumerate() {
        if o.origin == "vcf" && o.phased && o.hap_alleles.len() == 2 && o.hap_alleles[0] != o.hap_alleles[1] {
            if let (Some(ps), HapAllele::Cpic(s0), HapAllele::Cpic(s1)) =
                (&o.phase_set, &o.hap_alleles[0], &o.hap_alleles[1])
            {
                vcf_ps_groups.entry(ps.clone()).or_default().push((pid, s0.clone(), s1.clone()));
            }
        }
    }

    let has_phase = !bam_links.is_empty() || vcf_ps_groups.values().any(|v| v.len() >= 2);
    let phased_pairs: Vec<(&super::AlleleDef, &super::AlleleDef)> = valid_pairs
        .iter()
        .copied()
        .filter(|&(a1, a2)| {
            // Check VCF phase blocks
            for pts in vcf_ps_groups.values() {
                if pts.len() < 2 {
                    continue;
                }
                let fwd = pts.iter().all(|(pid, s0, s1)| {
                    cyp_accepted_set(gene, a1, *pid).iter().any(|x| x == s0)
                        && cyp_accepted_set(gene, a2, *pid).iter().any(|x| x == s1)
                });
                let rev = pts.iter().all(|(pid, s0, s1)| {
                    cyp_accepted_set(gene, a1, *pid).iter().any(|x| x == s1)
                        && cyp_accepted_set(gene, a2, *pid).iter().any(|x| x == s0)
                });
                if !fwd && !rev {
                    return false;
                }
            }
            // Check BAM fragment/mate phase links
            for ((pi, bi0, bi1), (pj, bj0, bj1)) in &bam_links {
                let fwd = cyp_accepted_set(gene, a1, *pi).iter().any(|x| x == bi0)
                    && cyp_accepted_set(gene, a1, *pj).iter().any(|x| x == bj0)
                    && cyp_accepted_set(gene, a2, *pi).iter().any(|x| x == bi1)
                    && cyp_accepted_set(gene, a2, *pj).iter().any(|x| x == bj1);
                let rev = cyp_accepted_set(gene, a1, *pi).iter().any(|x| x == bi1)
                    && cyp_accepted_set(gene, a1, *pj).iter().any(|x| x == bj1)
                    && cyp_accepted_set(gene, a2, *pi).iter().any(|x| x == bi0)
                    && cyp_accepted_set(gene, a2, *pj).iter().any(|x| x == bj0);
                if !fwd && !rev {
                    return false;
                }
            }
            true
        })
        .collect();

    let pool = if !phased_pairs.is_empty() { &phased_pairs } else { &valid_pairs };
    let used_phasing = !phased_pairs.is_empty() && has_phase;

    // Canonical IUPAC specificity penalty + functional certainty tiebreaker
    let score_pair = |a1: &super::AlleleDef, a2: &super::AlleleDef| -> (u32, u8, i64, (u8, u64, String), (u8, u64, String)) {
        let mut penalty = 0u32;
        for (pid, o) in obs.iter().enumerate() {
            if o.origin == "missing" || o.hap_alleles.len() != 2 {
                continue;
            }
            let ref_acc = cyp_accepted_set(gene, ref_allele, pid);
            let acc1 = cyp_accepted_set(gene, a1, pid);
            let acc2 = cyp_accepted_set(gene, a2, pid);
            let pos_key = gene.positions[pid].cpic_pos.to_string();
            let canon = cyp
                .canonical_iupac_defaults
                .get(&pos_key)
                .map(String::as_str)
                .unwrap_or_else(|| ref_cpic_value(&gene.positions[pid]));

            let eval_orient = |x0: &HapAllele, x1: &HapAllele| -> u32 {
                let mut pen = 0u32;
                if acc1.len() > 1 && acc1 != ref_acc {
                    if let HapAllele::Cpic(s) = x0 {
                        if s != canon {
                            pen += 1;
                        }
                    }
                }
                if acc2.len() > 1 && acc2 != ref_acc {
                    if let HapAllele::Cpic(s) = x1 {
                        if s != canon {
                            pen += 1;
                        }
                    }
                }
                pen
            };

            let h0 = &o.hap_alleles[0];
            let h1 = &o.hap_alleles[1];
            let m01 = allele_matches_base(&acc1, h0) && allele_matches_base(&acc2, h1);
            let m10 = allele_matches_base(&acc1, h1) && allele_matches_base(&acc2, h0);
            let p_min = match (m01, m10) {
                (true, true) => eval_orient(h0, h1).min(eval_orient(h1, h0)),
                (true, false) => eval_orient(h0, h1),
                (false, true) => eval_orient(h1, h0),
                (false, false) => 10,
            };
            penalty += p_min;
        }
        let uncertain = u8::from(a1.activity_value.is_none()) + u8::from(a2.activity_value.is_none());
        let freq_sum = -(((a1.max_pop_freq.unwrap_or(0.0) + a2.max_pop_freq.unwrap_or(0.0)) * 1_000_000.0).round() as i64);
        let k1 = star_sort_key(&a1.name, a1.is_reference);
        let k2 = star_sort_key(&a2.name, a2.is_reference);
        let (sk_a, sk_b) = if k1 <= k2 { (k1, k2) } else { (k2, k1) };
        (penalty, uncertain, freq_sum, sk_a, sk_b)
    };

    let mut best = pool[0];
    let mut best_key = score_pair(best.0, best.1);
    for &(a1, a2) in &pool[1..] {
        let k = score_pair(a1, a2);
        if k < best_key {
            best_key = k;
            best = (a1, a2);
        }
    }
    let k0 = star_sort_key(&best.0.name, best.0.is_reference);
    let k1 = star_sort_key(&best.1.name, best.1.is_reference);
    let (h1, h2) = if k0 <= k1 { (best.0, best.1) } else { (best.1, best.0) };
    Some((h1, h2, used_phasing))
}

fn cyp2d6_star_label(a: &super::AlleleDef) -> String {
    a.activity_label
        .as_deref()
        .map(str::trim)
        .filter(|s| !s.is_empty() && *s != "None")
        .map(ToString::to_string)
        .or_else(|| a.activity_value.map(format_activity_score))
        .unwrap_or_else(|| "n/a".to_string())
}

fn build_cyp2d6_outcome_from_labels(
    gene: &GeneDef,
    h1_name: &str,
    h1_label: &str,
    h2_name: &str,
    h2_label: &str,
    copy_number: f64,
    phased: bool,
    method: &str,
    obs: &[PosObservation],
    evidence: Vec<EvidenceEntry>,
    notes: Vec<String>,
    pheno_table: Option<&GenePhenotypeTable>,
) -> GeneOutcome {
    let total = gene.positions.len();
    let read = obs.iter().filter(|o| o.origin != "missing").count();
    let missing = total.saturating_sub(read);
    let (o1, o2) = if star_sort_key(h1_name, h1_name == "*1") <= star_sort_key(h2_name, h2_name == "*1") {
        ((h1_name, h1_label), (h2_name, h2_label))
    } else {
        ((h2_name, h2_label), (h1_name, h1_label))
    };

    let mut counts = BTreeMap::new();
    *counts.entry(o1.1.to_string()).or_insert(0usize) += 1;
    *counts.entry(o2.1.to_string()).or_insert(0usize) += 1;
    let (pheno, act_score, act_label) = lookup_phenotype_by_counts(&gene.lookup_method, &counts, pheno_table);
    let (act_score_final, act_label_final) = if act_score.is_none() {
        let v1: Option<f64> = o1.1.trim_start_matches('≥').parse().ok();
        let v2: Option<f64> = o2.1.trim_start_matches('≥').parse().ok();
        match (v1, v2) {
            (Some(x), Some(y)) => {
                let s = x + y;
                let has_ge = o1.1.starts_with('≥') || o2.1.starts_with('≥');
                let lbl = if has_ge {
                    format!("≥{}", format_activity_score(s))
                } else {
                    format_activity_score(s)
                };
                (Some(s), Some(lbl))
            }
            _ => (None, act_label),
        }
    } else {
        (act_score, act_label)
    };

    GeneOutcome {
        gene: gene.gene.clone(),
        called: true,
        diplotype: Some(format!("{}/{}", o1.0, o2.0)),
        phenotype: pheno,
        activity_score: act_score_final,
        activity_label: act_label_final,
        copy_number: Some(copy_number),
        source: "vcf+bam".to_string(),
        haplotype_1: Some(o1.0.to_string()),
        haplotype_2: Some(o2.0.to_string()),
        not_called_reason: None,
        calling_method: method.to_string(),
        phased,
        positions_total: total,
        positions_read: read,
        positions_missing: missing,
        evidence,
        notes,
        lookup_method: gene.lookup_method.clone(),
        ..Default::default()
    }
}

fn uncalled_cyp(
    gene: &GeneDef,
    method: &str,
    obs: &[PosObservation],
    evidence: Vec<EvidenceEntry>,
    notes: Vec<String>,
    reason: String,
) -> GeneOutcome {
    let mut g = GeneOutcome::uncalled(gene, method, obs, evidence, notes, reason);
    g.source = "vcf+bam".to_string();
    g
}

fn call_cyp2d6_bam_vcf(
    cyp: &Cyp2d6BamDef,
    vcf_obs: &VcfObservations,
    callable: Option<&super::CallableBed>,
    genome: &mut Option<referencia::Genoma>,
    pheno_table: Option<&GenePhenotypeTable>,
    pool: &BamEvidencePool,
) -> Result<GeneOutcome, String> {
    let gene = &cyp.base;
    let (mut obs_list, mut evidence, mut notes) = observe_gene_positions(gene, vcf_obs, callable, genome)?;
    merge_bam_pileup_with_vcf(gene, &mut obs_list, &mut evidence, &pool.pileup);

    let mut ctrl_depths: Vec<f64> = pool
        .ctrl_bases
        .iter()
        .zip(pool.ctrl_lens.iter())
        .map(|(&b, &l)| if l > 0 { b as f64 / l as f64 } else { 0.0 })
        .collect();
    let med_ctrl = median_f64(&mut ctrl_depths);
    if med_ctrl < 5.0 {
        return Ok(uncalled_cyp(
            gene,
            "vcf_bam",
            &obs_list,
            evidence,
            notes,
            format!("Insufficient BAM depth over diploid control regions (median {med_ctrl:.1}x < 5.0x)."),
        ));
    }

    let d_whole = pool.whole_bases as f64 / pool.whole_len.max(1) as f64;
    let d_5p = pool.fivep_bases as f64 / pool.fivep_len.max(1) as f64;
    let d_ex9 = pool.ex9_bases as f64 / pool.ex9_len.max(1) as f64;
    let d_paralog = pool.paralog_bases as f64 / pool.paralog_len.max(1) as f64;

    let cn_whole_raw = 2.0 * d_whole / med_ctrl;
    let cn_5p_raw = 2.0 * d_5p / med_ctrl;
    let cn_ex9_raw = 2.0 * d_ex9 / med_ctrl;
    let _ = &pool.ctrl_names;
    notes.push(format!(
        "BAM CNV depth: whole={cn_whole_raw:.2} ({d_whole:.1}x), 5'={cn_5p_raw:.2} ({d_5p:.1}x), \
         exon9={cn_ex9_raw:.2} ({d_ex9:.1}x), control_median={med_ctrl:.1}x, CYP2D7={d_paralog:.1}x"
    ));

    // Case 1: Homozygous whole-gene deletion (*5/*5)
    if cn_whole_raw < 0.35 {
        return Ok(build_cyp2d6_outcome_from_labels(
            gene, "*5", "0.0", "*5", "0.0", 0.0, true, "bam_cnv_deletion",
            &obs_list, evidence, notes, pheno_table,
        ));
    }

    let read_cnt = obs_list.iter().filter(|o| o.origin != "missing").count();
    let missing_keys: Vec<usize> = gene
        .positions
        .iter()
        .enumerate()
        .filter(|(i, p)| p.is_key && obs_list[*i].origin == "missing")
        .map(|(i, _)| i)
        .collect();
    if read_cnt == 0 || !missing_keys.is_empty() {
        return Ok(uncalled_cyp(
            gene,
            "vcf_bam",
            &obs_list,
            evidence,
            notes,
            format!(
                "CYP2D6 missing {} key position(s) across VCF and BAM.",
                missing_keys.len().max(usize::from(read_cnt == 0))
            ),
        ));
    }

    // Check exon-9 conversion hybrid drop (e.g. *10/*36 or *4/*68)
    let is_exon9_hybrid = (cn_5p_raw >= 1.50 && cn_ex9_raw < 1.35 && (cn_5p_raw - cn_ex9_raw) >= 0.50)
        || (cn_5p_raw >= 2.50 && (cn_5p_raw - cn_ex9_raw) >= 0.65 && (cn_ex9_raw / cn_5p_raw.max(0.01)) <= 0.72);
    if is_exon9_hybrid {
        let ex9_end = cyp.cnv_regions.target_ex9.end;
        // Build masked observation vector where exon 9 sites use homozygous intact allele
        // while 5'/body positions determine the hybrid backbone (*4 -> *68, *10/*1 -> *36).
        let is_star4 = gene
            .positions
            .iter()
            .enumerate()
            .any(|(i, p)| {
                p.pos == 42128945
                    && obs_list[i].hap_alleles.iter().any(|h| matches!(h, HapAllele::Cpic(b) if b == "T"))
            });
        let hybrid_name = if is_star4 { "*68" } else { "*36" };

        // Determine intact haplotype from exon-9 & non-deleted positions
        let mut intact_obs = obs_list.clone();
        for (i, p) in gene.positions.iter().enumerate() {
            if p.pos <= ex9_end && pool.pileup[i].total_bases() >= 4 {
                let sp = &pool.pileup[i];
                let bases = ["A", "C", "G", "T"];
                let best_bi = (0..4).max_by_key(|&bi| sp.counts[bi]).unwrap_or(0);
                let b = bases[best_bi].to_string();
                intact_obs[i].hap_alleles = vec![HapAllele::Cpic(b.clone()), HapAllele::Cpic(b)];
            }
        }
        // Match single best allele using exon-9-intact homozygous sites or non-hybrid positions
        let intact_allele = gene
            .alleles
            .iter()
            .filter(|a| {
                gene.positions.iter().enumerate().all(|(i, _)| {
                    if intact_obs[i].origin == "missing" || intact_obs[i].hap_alleles.is_empty() {
                        return true;
                    }
                    let acc = cyp_accepted_set(gene, a, i);
                    if gene.positions[i].pos <= ex9_end {
                        allele_matches_base(&acc, &intact_obs[i].hap_alleles[0])
                    } else {
                        obs_list[i].hap_alleles.iter().any(|ha| allele_matches_base(&acc, ha))
                    }
                })
            })
            .min_by_key(|a| star_sort_key(&a.name, a.is_reference));

        let (intact_nm, intact_lbl) = if let Some(ia) = intact_allele {
            (ia.name.as_str(), cyp2d6_star_label(ia))
        } else {
            ("*1", "1.0".to_string())
        };
        let cn_rounded = (cn_5p_raw.round()).clamp(2.0, 4.0);
        if cn_rounded >= 3.0 && intact_nm == "*10" && hybrid_name == "*36" {
            if let Some(tandem) = cyp.duplication_alleles.get("*36+*10") {
                return Ok(build_cyp2d6_outcome_from_labels(
                    gene,
                    &tandem.name,
                    tandem.activity_label.as_deref().unwrap_or("0.25"),
                    intact_nm,
                    &intact_lbl,
                    cn_rounded,
                    false,
                    "bam_cnv_exon9_hybrid",
                    &obs_list,
                    evidence,
                    notes,
                    pheno_table,
                ));
            }
        }
        return Ok(build_cyp2d6_outcome_from_labels(
            gene,
            intact_nm,
            &intact_lbl,
            hybrid_name,
            "0.0",
            cn_rounded,
            false,
            "bam_cnv_exon9_hybrid",
            &obs_list,
            evidence,
            notes,
            pheno_table,
        ));
    }

    // Case 2: Heterozygous whole-gene deletion (*5 + one functional/non-functional copy)
    if cn_whole_raw < 1.35 {
        let mut single_obs = obs_list.clone();
        // Enforce majority BAM allele at any slightly noisy sites under CN=1
        for (i, o) in single_obs.iter_mut().enumerate() {
            if o.hap_alleles.len() == 2 && o.hap_alleles[0] != o.hap_alleles[1] {
                let sp = &pool.pileup[i];
                if sp.total_bases() >= 4 {
                    let bases = ["A", "C", "G", "T"];
                    let best_bi = (0..4).max_by_key(|&bi| sp.counts[bi]).unwrap_or(0);
                    let b = HapAllele::Cpic(bases[best_bi].to_string());
                    o.hap_alleles = vec![b.clone(), b];
                }
            }
        }
        if let Some((a1, _, _)) = resolve_cyp2d6_diploid(cyp, &single_obs, &pool.fragments) {
            return Ok(build_cyp2d6_outcome_from_labels(
                gene,
                &a1.name,
                &cyp2d6_star_label(a1),
                "*5",
                "0.0",
                1.0,
                true,
                "bam_cnv_deletion",
                &obs_list,
                evidence,
                notes,
                pheno_table,
            ));
        }
        return Ok(uncalled_cyp(
            gene,
            "bam_cnv_deletion",
            &obs_list,
            evidence,
            notes,
            "CYP2D6 copy number ~1 (*5 deletion) detected, but remaining allele haplotype could not be matched cleanly.".into(),
        ));
    }

    // Case 3: Diploid (CN=2) or Gene Duplication (CN>=3)
    let Some((a1, a2, used_phasing)) = resolve_cyp2d6_diploid(cyp, &obs_list, &pool.fragments) else {
        let vcf_only = call_combinatorial_gene(gene, &obs_list, evidence.clone(), notes.clone(), pheno_table, "vcf_bam", false, None);
        return Ok(GeneOutcome {
            source: "vcf+bam".to_string(),
            copy_number: Some(cn_whole_raw.round().clamp(1.0, 8.0)),
            ..vcf_only
        });
    };

    if cn_whole_raw < 2.55 {
        return Ok(build_cyp2d6_outcome_from_labels(
            gene,
            &a1.name,
            &cyp2d6_star_label(a1),
            &a2.name,
            &cyp2d6_star_label(a2),
            2.0,
            used_phasing,
            "vcf_bam_phased",
            &obs_list,
            evidence,
            notes,
            pheno_table,
        ));
    }

    // Gene duplication (CN >= 3): determine duplicated haplotype from BAM VAF skew at het sites
    let cn_int = (cn_whole_raw.round() as usize).clamp(3, 8);
    let (c1, c2) = if a1.name == a2.name {
        (cn_int - 1, 1usize)
    } else {
        let mut r1 = 0u64;
        let mut r2 = 0u64;
        for (idx, p) in gene.positions.iter().enumerate() {
            let acc1 = cyp_accepted_set(gene, a1, idx);
            let acc2 = cyp_accepted_set(gene, a2, idx);
            if acc1.len() == 1 && acc2.len() == 1 && acc1[0] != acc2[0] && p.grch38_ref.len() == 1 {
                if let (Some(bi1), Some(bi2)) = (
                    SitePileup::base_idx(acc1[0].as_bytes()[0]),
                    SitePileup::base_idx(acc2[0].as_bytes()[0]),
                ) {
                    let sp = &pool.pileup[idx];
                    if sp.counts[bi1] + sp.counts[bi2] >= 8 {
                        r1 += u64::from(sp.counts[bi1]);
                        r2 += u64::from(sp.counts[bi2]);
                    }
                }
            }
        }
        let frac1 = if r1 + r2 > 0 { r1 as f64 / (r1 + r2) as f64 } else { 0.5 };
        notes.push(format!(
            "BAM duplication allele balance for {} vs {}: {:.2} ({r1}/{})",
            a1.name,
            a2.name,
            frac1,
            r1 + r2
        ));
        if cn_int == 3 {
            if frac1 >= 0.5 { (2, 1) } else { (1, 2) }
        } else {
            let best_k1 = (1..cn_int)
                .min_by(|&k, &m| {
                    let dfk = ((k as f64 / cn_int as f64) - frac1).abs();
                    let dfm = ((m as f64 / cn_int as f64) - frac1).abs();
                    dfk.total_cmp(&dfm)
                })
                .unwrap_or(cn_int - 1);
            (best_k1, cn_int - best_k1)
        }
    };

    let format_dup = |al: &super::AlleleDef, copies: usize| -> (String, String) {
        if copies <= 1 {
            return (al.name.clone(), cyp2d6_star_label(al));
        }
        let exact_key = format!("{}x{}", al.name, copies);
        if let Some(d) = cyp.duplication_alleles.get(&exact_key) {
            return (
                d.name.clone(),
                d.activity_label
                    .clone()
                    .or_else(|| d.activity_value.map(format_activity_score))
                    .unwrap_or_else(|| "n/a".to_string()),
            );
        }
        if copies >= 3 {
            let ge_key = format!("{}x≥3", al.name);
            if let Some(d) = cyp.duplication_alleles.get(&ge_key) {
                return (
                    d.name.clone(),
                    d.activity_label
                        .clone()
                        .or_else(|| d.activity_value.map(format_activity_score))
                        .unwrap_or_else(|| "n/a".to_string()),
                );
            }
        }
        let scaled_lbl = al
            .activity_value
            .map(|v| format_activity_score(v * copies as f64))
            .unwrap_or_else(|| "n/a".to_string());
        (exact_key, scaled_lbl)
    };

    let (d1_nm, d1_lbl) = format_dup(a1, c1);
    let (d2_nm, d2_lbl) = format_dup(a2, c2);
    Ok(build_cyp2d6_outcome_from_labels(
        gene,
        &d1_nm,
        &d1_lbl,
        &d2_nm,
        &d2_lbl,
        cn_int as f64,
        used_phasing,
        "bam_cnv_duplication",
        &obs_list,
        evidence,
        notes,
        pheno_table,
    ))
}

// ---------------------------------------------------------------------------
// IPD-IMGT/HLA Cache Resolver & 2-Stage Pure-Rust Class I HLA Caller
// ---------------------------------------------------------------------------

fn resolve_hla_cache_dir(bundle_dir: &Path) -> PathBuf {
    if let Ok(env_dir) = std::env::var("MMGEN_HLA_DIR") {
        let t = env_dir.trim();
        if !t.is_empty() {
            return PathBuf::from(t);
        }
    }
    if let Some(parent) = bundle_dir.parent() {
        return parent.join("hla_cache");
    }
    for home_var in ["HOME", "USERPROFILE"] {
        if let Ok(home) = std::env::var(home_var) {
            return Path::new(&home).join(".cache").join("mmgen").join("hla_cache");
        }
    }
    std::env::temp_dir().join("mmgen_hla_cache")
}

fn validate_imgt_fasta(path: &Path) -> Result<(), String> {
    let f = File::open(path).map_err(|e| format!("{}: cannot open IMGT/HLA FASTA ({e}).", path.display()))?;
    let mut reader = BufReader::new(f);
    let mut first_line = String::new();
    std::io::BufRead::read_line(&mut reader, &mut first_line)
        .map_err(|e| format!("{}: failed reading IMGT/HLA FASTA header ({e}).", path.display()))?;
    if !first_line.starts_with('>') {
        return Err(format!(
            "{}: invalid FASTA file (expected leading '>' header line).",
            path.display()
        ));
    }
    Ok(())
}

fn download_file_curl(url: &str, dest_tmp: &Path) -> Result<(), String> {
    let mut last_err = String::new();
    for prog in ["curl", "curl.exe"] {
        let status = Command::new(prog)
            .args([
                "-fsSL",
                "--connect-timeout",
                "15",
                "--max-time",
                "300",
                "-o",
            ])
            .arg(dest_tmp)
            .arg(url)
            .output();
        match status {
            Ok(out) if out.status.success() => return Ok(()),
            Ok(out) => {
                let err_txt = String::from_utf8_lossy(&out.stderr);
                last_err = format!("{prog} exited with {}: {}", out.status, err_txt.trim());
            }
            Err(e) => {
                last_err = format!("{prog} unavailable ({e})");
            }
        }
    }
    Err(last_err)
}

/// Whether `path` is the pinned release's file (always true when the tables pin nothing).
fn imgt_fasta_matches_pin(path: &Path, hg: &HlaBamDef) -> Result<bool, String> {
    let Some(want) = hg.fasta_sha256.as_deref() else {
        return Ok(true);
    };
    if let Some(bytes) = hg.fasta_bytes {
        let len = std::fs::metadata(path).map_err(|e| format!("{}: cannot read file size ({e}).", path.display()))?.len();
        if len != bytes {
            return Ok(false);
        }
    }
    let data = std::fs::read(path).map_err(|e| format!("{}: cannot read IMGT/HLA FASTA ({e}).", path.display()))?;
    Ok(super::sha256_hex(&data).eq_ignore_ascii_case(want))
}

/// "IPD-IMGT/HLA 3.66.0 (SHA-256 f265427f646d...)", or "" when the tables pin nothing.
fn imgt_pin_label(hg: &HlaBamDef) -> String {
    match (&hg.imgt_hla_version, &hg.fasta_sha256) {
        (Some(v), Some(s)) => format!("IPD-IMGT/HLA {v} (SHA-256 {}...)", &s[..s.len().min(12)]),
        (Some(v), None) => format!("IPD-IMGT/HLA {v}"),
        _ => String::new(),
    }
}

fn resolve_or_download_imgt_fasta(bundle_dir: &Path, hg: &HlaBamDef) -> Result<PathBuf, String> {
    let cache_dir = resolve_hla_cache_dir(bundle_dir);
    let target_dir = match &hg.imgt_hla_version {
        Some(v) => cache_dir.join(v),
        None => cache_dir.clone(),
    };
    let mut candidates = vec![target_dir.join(&hg.fasta_filename)];
    if target_dir != cache_dir {
        candidates.push(cache_dir.join(&hg.fasta_filename));
    }
    // A file of another release is never used: the tables were validated with the pinned one.
    let mut other_release: Vec<String> = Vec::new();
    for cand in &candidates {
        if cand.is_file() {
            validate_imgt_fasta(cand)?;
            if imgt_fasta_matches_pin(cand, hg)? {
                return Ok(cand.clone());
            }
            other_release.push(cand.display().to_string());
        }
    }
    let pin = imgt_pin_label(hg);
    let wrong_note = if other_release.is_empty() {
        String::new()
    } else {
        format!(
            " ('{}' is not that release's {}: its SHA-256 differs, so it was not used.)",
            other_release.join("', '"),
            hg.fasta_filename
        )
    };
    let fasta_path = target_dir.join(&hg.fasta_filename);

    let offline = std::env::var("MMGEN_HLA_OFFLINE")
        .map(|v| matches!(v.trim(), "1" | "true" | "TRUE"))
        .unwrap_or(false);
    if offline {
        return Err(format!(
            "{}: missing cached IPD-IMGT/HLA file '{}'{} (MMGEN_HLA_OFFLINE=1). \
             Download unmodified FASTA from {} and place it in '{}' (or set MMGEN_HLA_DIR=<dir>).{wrong_note}",
            hg.gene,
            hg.fasta_filename,
            if pin.is_empty() { String::new() } else { format!(" of {pin}") },
            hg.fasta_url,
            target_dir.display()
        ));
    }

    if let Err(e) = std::fs::create_dir_all(&target_dir) {
        return Err(format!(
            "{}: cannot create HLA cache directory '{}' ({e}). Set MMGEN_HLA_DIR=<writable_dir>.",
            hg.gene,
            target_dir.display()
        ));
    }
    let tmp_path = target_dir.join(format!("{}.tmp.{}", hg.fasta_filename, std::process::id()));
    let _ = std::fs::remove_file(&tmp_path);
    match download_file_curl(&hg.fasta_url, &tmp_path) {
        Ok(()) => {
            if let Err(e) = validate_imgt_fasta(&tmp_path) {
                let _ = std::fs::remove_file(&tmp_path);
                return Err(e);
            }
            match imgt_fasta_matches_pin(&tmp_path, hg) {
                Ok(true) => {}
                Ok(false) => {
                    let _ = std::fs::remove_file(&tmp_path);
                    return Err(format!(
                        "{}: the {} downloaded from {} is not the file of {pin} that the tables were validated with \
                         (its SHA-256 or size differs); nothing was saved. Try again, or download it by hand into '{}'.",
                        hg.gene,
                        hg.fasta_filename,
                        hg.fasta_url,
                        target_dir.display()
                    ));
                }
                Err(e) => {
                    let _ = std::fs::remove_file(&tmp_path);
                    return Err(e);
                }
            }
            if let Err(e) = std::fs::rename(&tmp_path, &fasta_path) {
                let _ = std::fs::remove_file(&tmp_path);
                return Err(format!("{}: failed saving downloaded IMGT FASTA ({e}).", fasta_path.display()));
            }
            Ok(fasta_path)
        }
        Err(err) => {
            let _ = std::fs::remove_file(&tmp_path);
            Err(format!(
                "{}: missing IPD-IMGT/HLA FASTA '{}'{} and automatic download from {} failed ({err}). \
                 When offline, download unmodified {} manually into '{}' or set MMGEN_HLA_DIR=<dir>.{wrong_note}",
                hg.gene,
                fasta_path.display(),
                if pin.is_empty() { String::new() } else { format!(" of {pin}") },
                hg.fasta_url,
                hg.fasta_filename,
                target_dir.display()
            ))
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct HlaSortKey {
    suffix_rank: u8,
    f1: u32,
    f2: u32,
    f3: u32,
    f4: u32,
    full_name: String,
}

fn parse_hla_sort_key(full_name: &str) -> HlaSortKey {
    let body = full_name.split_once('*').map(|x| x.1).unwrap_or(full_name);
    let mut trimmed = body;
    let mut has_suffix = false;
    while let Some(last) = trimmed.chars().next_back() {
        if last.is_ascii_alphabetic() {
            has_suffix = true;
            trimmed = &trimmed[..trimmed.len() - last.len_utf8()];
        } else {
            break;
        }
    }
    let mut nums = [0u32; 4];
    for (idx, part) in trimmed.split(':').take(4).enumerate() {
        let digits: String = part.chars().take_while(|c| c.is_ascii_digit()).collect();
        nums[idx] = digits.parse::<u32>().unwrap_or(999_999);
    }
    HlaSortKey {
        suffix_rank: u8::from(has_suffix),
        f1: nums[0],
        f2: nums[1],
        f3: nums[2],
        f4: nums[3],
        full_name: full_name.to_string(),
    }
}

fn two_field_name(full_name: &str) -> String {
    let parts: Vec<&str> = full_name.split(':').collect();
    if parts.len() >= 2 {
        format!("{}:{}", parts[0], parts[1])
    } else {
        full_name.to_string()
    }
}

fn strip_hla_gene_prefix(two_field: &str) -> String {
    if let Some((_, rest)) = two_field.split_once('*') {
        format!("*{rest}")
    } else {
        two_field.to_string()
    }
}

struct ImgtEntry {
    sort_key: HlaSortKey,
    two_field: String,
    family_key: String,
    seq: Vec<u8>,
}

fn load_unique_full_imgt(fasta_path: &Path) -> Result<Vec<ImgtEntry>, String> {
    let f = File::open(fasta_path).map_err(|e| format!("{}: cannot open IMGT FASTA ({e}).", fasta_path.display()))?;
    let reader = BufReader::new(f);

    let mut raw_entries: Vec<(HlaSortKey, String, Vec<u8>)> = Vec::new();
    let mut cur_name: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::new();
    let mut max_len = 0usize;

    for line_res in std::io::BufRead::lines(reader) {
        let line = line_res.map_err(|e| format!("{}: read error ({e}).", fasta_path.display()))?;
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        if let Some(hdr) = t.strip_prefix('>') {
            if let Some(full) = cur_name.take() {
                if cur_seq.len() > max_len {
                    max_len = cur_seq.len();
                }
                raw_entries.push((parse_hla_sort_key(&full), full, std::mem::take(&mut cur_seq)));
            }
            if let Some(allele_tok) = hdr.split_whitespace().nth(1) {
                cur_name = Some(allele_tok.to_string());
            }
        } else if cur_name.is_some() {
            cur_seq.extend(t.bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    if let Some(full) = cur_name.take() {
        if cur_seq.len() > max_len {
            max_len = cur_seq.len();
        }
        raw_entries.push((parse_hla_sort_key(&full), full, cur_seq));
    }

    // Real IMGT full-length genomic sequences are ~3.2-3.5kb; keep >=2400bp (unless synthetic tiny FASTA).
    let min_len = if max_len >= 2400 { 2400 } else { 100 };
    raw_entries.retain(|(_, _, s)| s.len() >= min_len);
    raw_entries.sort_by(|a, b| a.0.cmp(&b.0));

    let mut seen_seqs: HashSet<Vec<u8>> = HashSet::new();
    let mut out = Vec::with_capacity(raw_entries.len());
    for (sk, full_name, seq) in raw_entries {
        if seen_seqs.insert(seq.clone()) {
            let tf = two_field_name(&full_name);
            let fam = tf.split(':').next().unwrap_or(&tf).to_string();
            out.push(ImgtEntry { sort_key: sk, two_field: tf, family_key: fam, seq });
        }
    }
    if out.is_empty() {
        return Err(format!("{}: no valid full-length HLA sequences found.", fasta_path.display()));
    }
    Ok(out)
}

#[inline]
fn base_to_2bit(b: u8) -> Option<u64> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

#[inline]
fn revcomp_bytes(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

fn push_2bit_kmers(seq: &[u8], k: usize, mask: u64, counts: &mut HashMap<u64, u32>) {
    let mut val = 0u64;
    let mut valid = 0usize;
    for &b in seq {
        if let Some(code) = base_to_2bit(b) {
            val = ((val << 2) | code) & mask;
            valid += 1;
            if valid >= k {
                *counts.entry(val).or_insert(0) += 1;
            }
        } else {
            val = 0;
            valid = 0;
        }
    }
}

fn eval_hla_allele_status(hg: &HlaBamDef, h1: &str, h2: &str) -> Option<String> {
    if hg.lookup_method == "ALLELE_STATUS" && !hg.risk_alleles.is_empty() {
        let parts: Vec<String> = hg
            .risk_alleles
            .iter()
            .map(|ra| {
                let pos = h1 == ra || h2 == ra;
                format!("{ra} {}", if pos { "positive" } else { "negative" })
            })
            .collect();
        Some(parts.join("; "))
    } else {
        None
    }
}

fn call_hla_gene(
    bundle_dir: &Path,
    hg: &HlaBamDef,
    pairs_map: &HashMap<String, Vec<(u8, Vec<u8>)>>,
) -> Result<GeneOutcome, String> {
    let fasta_path = resolve_or_download_imgt_fasta(bundle_dir, hg)?;
    let entries = load_unique_full_imgt(&fasta_path)?;

    let uncalled_hla = |reason: String, mut notes: Vec<String>| -> GeneOutcome {
        notes.insert(0, reason.clone());
        GeneOutcome {
            gene: hg.gene.clone(),
            called: false,
            diplotype: None,
            phenotype: Some("No Result".to_string()),
            activity_score: None,
            activity_label: None,
            copy_number: None,
            source: "bam".to_string(),
            haplotype_1: None,
            haplotype_2: None,
            not_called_reason: Some(reason),
            calling_method: "imgt_kmer_paired".to_string(),
            phased: false,
            positions_total: 0,
            positions_read: 0,
            positions_missing: 0,
            evidence: Vec::new(),
            notes,
            lookup_method: hg.lookup_method.clone(),
            ..Default::default()
        }
    };

    if pairs_map.is_empty() {
        return Ok(uncalled_hla(
            format!("No primary mapped BAM reads overlapping {} ({}:{}-{}).", hg.gene, hg.chrom, hg.start, hg.end),
            Vec::new(),
        ));
    }

    // Stage 1: K=21 solid k-mers (count >= 3) across reads + revcomps
    const K1: usize = 21;
    let mask1: u64 = (1u64 << (2 * K1)) - 1;
    let mut kmer_counts: HashMap<u64, u32> = HashMap::with_capacity(16_384);
    for rlist in pairs_map.values() {
        for (_, seq) in rlist {
            push_2bit_kmers(seq, K1, mask1, &mut kmer_counts);
            let rc = revcomp_bytes(seq);
            push_2bit_kmers(&rc, K1, mask1, &mut kmer_counts);
        }
    }
    let solid: HashSet<u64> = kmer_counts.into_iter().filter(|&(_, c)| c >= 3).map(|(k, _)| k).collect();
    if solid.is_empty() {
        return Ok(uncalled_hla(
            format!("Insufficient BAM read depth at {} to build solid {K1}-mers.", hg.gene),
            Vec::new(),
        ));
    }

    // Score each allele sequence (best representation per 2-field allele)
    struct Stage1Cand<'a> {
        score_x2: i64, // 2*hits - misses (exact integer equivalent to hits - 0.5*misses)
        neg_max_zero: i64,
        entry: &'a ImgtEntry,
    }

    let mut best_per_2f: HashMap<&str, Stage1Cand<'_>> = HashMap::with_capacity(4096);
    for entry in &entries {
        if entry.seq.len() < K1 {
            continue;
        }
        let nkm = entry.seq.len() - K1 + 1;
        let mut val = 0u64;
        let mut valid = 0usize;
        let mut hits = 0i64;
        let mut misses = 0i64;
        let mut cur_zero = 0i64;
        let mut max_zero = 0i64;
        for (i, &b) in entry.seq.iter().enumerate() {
            if let Some(code) = base_to_2bit(b) {
                val = ((val << 2) | code) & mask1;
                valid += 1;
            } else {
                val = 0;
                valid = 0;
            }
            if i + 1 >= K1 && (i + 1 - K1) % 2 == 0 {
                let hit = valid >= K1 && solid.contains(&val);
                if hit {
                    hits += 1;
                    cur_zero = 0;
                } else {
                    misses += 1;
                    cur_zero += 1;
                    if cur_zero > max_zero {
                        max_zero = cur_zero;
                    }
                }
            }
            let _ = nkm;
        }
        let cand = Stage1Cand { score_x2: 2 * hits - misses, neg_max_zero: -max_zero, entry };
        let replace = match best_per_2f.get(entry.two_field.as_str()) {
            None => true,
            Some(prev) => {
                (cand.score_x2, cand.neg_max_zero) > (prev.score_x2, prev.neg_max_zero)
                    || ((cand.score_x2, cand.neg_max_zero) == (prev.score_x2, prev.neg_max_zero)
                        && cand.entry.sort_key < prev.entry.sort_key)
            }
        };
        if replace {
            best_per_2f.insert(entry.two_field.as_str(), cand);
        }
    }

    let mut ranked: Vec<Stage1Cand<'_>> = best_per_2f.into_values().collect();
    ranked.sort_by(|a, b| {
        b.score_x2
            .cmp(&a.score_x2)
            .then_with(|| b.neg_max_zero.cmp(&a.neg_max_zero))
            .then_with(|| a.entry.sort_key.cmp(&b.entry.sort_key))
    });

    // Diversify: top 40 overall + up to 4 best per 1-field family up to 75 total candidates
    let mut selected: Vec<&ImgtEntry> = Vec::with_capacity(75);
    let mut seen_tf: HashSet<&str> = HashSet::new();
    let mut fam_counts: HashMap<&str, usize> = HashMap::new();
    for item in ranked.iter().take(40) {
        selected.push(item.entry);
        seen_tf.insert(item.entry.two_field.as_str());
        *fam_counts.entry(item.entry.family_key.as_str()).or_insert(0) += 1;
    }
    for item in ranked.iter().skip(40) {
        if selected.len() >= 75 {
            break;
        }
        let fc = fam_counts.get(item.entry.family_key.as_str()).copied().unwrap_or(0);
        if fc < 4 && seen_tf.insert(item.entry.two_field.as_str()) {
            selected.push(item.entry);
            *fam_counts.entry(item.entry.family_key.as_str()).or_insert(0) += 1;
        }
    }

    // Stage 2: S=19 anchor seeds & paired-read ungapped alignment scoring
    const S: usize = 19;
    let mask_s: u64 = (1u64 << (2 * S)) - 1;
    let mut cand_idx: HashMap<u64, Vec<(u16, i32)>> = HashMap::with_capacity(4096);
    for (ci, entry) in selected.iter().enumerate() {
        let mut val = 0u64;
        let mut valid = 0usize;
        for (p, &b) in entry.seq.iter().enumerate() {
            if let Some(code) = base_to_2bit(b) {
                val = ((val << 2) | code) & mask_s;
                valid += 1;
                if valid >= S {
                    cand_idx.entry(val).or_default().push((ci as u16, (p + 1 - S) as i32));
                }
            } else {
                val = 0;
                valid = 0;
            }
        }
    }

    let has_pairs = pairs_map.values().any(|v| v.len() >= 2);
    let nc = selected.len();
    let mut pair_scores: Vec<Vec<f64>> = Vec::with_capacity(pairs_map.len());
    let mut offsets_by_cand: Vec<Vec<i32>> = vec![Vec::new(); nc];

    for rlist in pairs_map.values() {
        if has_pairs && rlist.len() < 2 {
            continue;
        }
        let mut m_best: [Vec<Option<f64>>; 2] = [vec![None; nc], vec![None; nc]];
        let mut any_hit = false;
        let mates_n = if has_pairs { 2 } else { 1 };

        for (mi, (_, seq)) in rlist.iter().take(mates_n).enumerate() {
            let l = seq.len();
            if l < 50 {
                continue;
            }
            let min_ov = 70usize.min((l * 3) / 4);
            let rc = revcomp_bytes(seq);
            let qoffs = [20usize.min(l / 4), (l / 2).saturating_sub(9), l.saturating_sub(40)];

            for orient in [seq.as_slice(), rc.as_slice()] {
                for v in &mut offsets_by_cand {
                    v.clear();
                }
                let mut touched = Vec::new();
                for &qoff in &qoffs {
                    if qoff + S > l {
                        continue;
                    }
                    let mut seed_val = 0u64;
                    let mut ok_seed = true;
                    for &b in &orient[qoff..qoff + S] {
                        if let Some(c) = base_to_2bit(b) {
                            seed_val = (seed_val << 2) | c;
                        } else {
                            ok_seed = false;
                            break;
                        }
                    }
                    if ok_seed {
                        if let Some(hits) = cand_idx.get(&seed_val) {
                            for &(ci, apos) in hits {
                                let off = apos - qoff as i32;
                                let ov_vec = &mut offsets_by_cand[ci as usize];
                                if ov_vec.is_empty() {
                                    touched.push(ci as usize);
                                }
                                if !ov_vec.contains(&off) {
                                    ov_vec.push(off);
                                }
                            }
                        }
                    }
                }

                for ci in touched {
                    let aseq = &selected[ci].seq;
                    let alen = aseq.len() as i32;
                    let mut best_mm = usize::MAX;
                    let mut best_ov = 0usize;
                    for &off in &offsets_by_cand[ci] {
                        let s_i = (-off).max(0) as usize;
                        let e_i = ((l as i32).min(alen - off)).max(0) as usize;
                        if e_i <= s_i {
                            continue;
                        }
                        let ov = e_i - s_i;
                        if ov < min_ov {
                            continue;
                        }
                        let mut mm = 0usize;
                        let a_base = (off + s_i as i32) as usize;
                        for k in 0..ov {
                            if orient[s_i + k] != aseq[a_base + k] {
                                mm += 1;
                            }
                        }
                        if mm < best_mm || (mm == best_mm && ov > best_ov) {
                            best_mm = mm;
                            best_ov = ov;
                        }
                    }
                    if best_ov >= min_ov && best_mm <= (best_ov as f64 * 0.05) as usize {
                        let sc = best_ov as f64 - 8.0 * best_mm as f64;
                        if m_best[mi][ci].is_none_or(|p| sc > p) {
                            m_best[mi][ci] = Some(sc);
                            any_hit = true;
                        }
                    }
                }
            }
        }

        if any_hit {
            let mut combined = Vec::with_capacity(nc);
            for ci in 0..nc {
                let s = if has_pairs {
                    m_best[0][ci].unwrap_or(-20.0) + m_best[1][ci].unwrap_or(-20.0)
                } else {
                    m_best[0][ci].unwrap_or(-20.0)
                };
                combined.push(s);
            }
            pair_scores.push(combined);
        }
    }

    if pair_scores.is_empty() {
        return Ok(uncalled_hla(
            format!("{} reads did not align to any full-length IMGT/HLA allele in {}.", hg.gene, fasta_path.display()),
            Vec::new(),
        ));
    }

    // Joint diploid maximization over all (i <= j) pairs with mild heterozygosity parsimony (-30)
    let mut best_ij = (0usize, 0usize);
    let mut best_tot = f64::NEG_INFINITY;
    for i in 0..nc {
        for j in i..nc {
            let mut tot: f64 = pair_scores.iter().map(|ps| ps[i].max(ps[j])).sum();
            if i != j {
                tot -= 30.0;
            }
            let cmp_better = if tot > best_tot + 1e-6 {
                true
            } else if (tot - best_tot).abs() <= 1e-6 {
                let (bi, bj) = best_ij;
                let (na, nb) = if selected[i].sort_key <= selected[j].sort_key {
                    (&selected[i].sort_key, &selected[j].sort_key)
                } else {
                    (&selected[j].sort_key, &selected[i].sort_key)
                };
                let (oa, ob) = if selected[bi].sort_key <= selected[bj].sort_key {
                    (&selected[bi].sort_key, &selected[bj].sort_key)
                } else {
                    (&selected[bj].sort_key, &selected[bi].sort_key)
                };
                (na, nb) < (oa, ob)
            } else {
                false
            };
            if cmp_better {
                best_tot = tot;
                best_ij = (i, j);
            }
        }
    }

    let (bi, bj) = best_ij;
    let (ea, eb) = if selected[bi].sort_key <= selected[bj].sort_key {
        (selected[bi], selected[bj])
    } else {
        (selected[bj], selected[bi])
    };
    let h1 = strip_hla_gene_prefix(&ea.two_field);
    let h2 = strip_hla_gene_prefix(&eb.two_field);
    let pheno = eval_hla_allele_status(hg, &h1, &h2);
    let pin = imgt_pin_label(hg);
    let notes = vec![format!(
        "IMGT/HLA cache: {}{} ({} unique full alleles, {} candidates, {} aligned pairs, score={best_tot:.1})",
        fasta_path.file_name().and_then(|s| s.to_str()).unwrap_or(&hg.fasta_filename),
        if pin.is_empty() { String::new() } else { format!(", {pin}") },
        entries.len(),
        nc,
        pair_scores.len()
    )];

    Ok(GeneOutcome {
        gene: hg.gene.clone(),
        called: true,
        diplotype: Some(format!("{h1}/{h2}")),
        phenotype: pheno,
        activity_score: None,
        activity_label: None,
        copy_number: Some(2.0),
        source: "bam".to_string(),
        haplotype_1: Some(h1),
        haplotype_2: Some(h2),
        not_called_reason: None,
        calling_method: "imgt_kmer_paired".to_string(),
        phased: false,
        positions_total: 0,
        positions_read: 0,
        positions_missing: 0,
        evidence: Vec::new(),
        notes,
        lookup_method: hg.lookup_method.clone(),
        ..Default::default()
    })
}

/// Call `CYP2D6` (BAM+VCF) and Class I `HLA-A/HLA-B/HLA-C` (BAM+IMGT cache).
pub(crate) fn call_bam_genes(
    bundle: &Bundle,
    bam_path: &str,
    vcf_obs: &VcfObservations,
    callable: Option<&super::CallableBed>,
    genome: &mut Option<referencia::Genoma>,
) -> Result<Vec<GeneOutcome>, String> {
    let bam_def = bundle.bam_genes.as_ref().ok_or_else(|| {
        format!(
            "{}: BAM calling requested, but table bundle lacks bam_genes.json; upgrade or rebuild PGx tables for v2.",
            bundle.dir.display()
        )
    })?;

    let pool = collect_bam_evidence(bam_path, bam_def)?;
    let mut outcomes = Vec::with_capacity(1 + bam_def.hla_genes.len());

    let cyp_outcome = call_cyp2d6_bam_vcf(
        &bam_def.cyp2d6,
        vcf_obs,
        callable,
        genome,
        bundle.phenotypes.genes.get("CYP2D6"),
        &pool,
    )?;
    outcomes.push(cyp_outcome);

    for (idx, hg) in bam_def.hla_genes.iter().enumerate() {
        let out = call_hla_gene(&bundle.dir, hg, &pool.hla_reads[idx])?;
        outcomes.push(out);
    }
    Ok(outcomes)
}
