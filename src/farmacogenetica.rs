//! Native pharmacogenetics engine (`mmgen-pgx`).
//!
//! Calls star alleles, variant activity/status diplotypes, CPIC phenotypes,
//! activity scores, and drug guideline recommendations directly from a GRCh38
//! VCF (plain, `.gz`/BGZF with optional `.tbi`, or `.lz4`), optional callable
//! BED, and optional GRCh38 reference FASTA, using external CC0-1.0 JSON table
//! bundles verified by SHA-256 against `manifest.json`.

#[path = "pgx_bam.rs"]
pub mod pgx_bam;

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::BufRead;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use serde_json::json;

use crate::{externo, indice, nomes, referencia};
pub use pgx_bam::BamGenesBundle;

pub const ENGINE_NAME: &str = "mmgen-pgx";
pub const ENGINE_VERSION: &str = "1.1.0";
pub const SCHEMA_VERSION: &str = "1";
pub const DEFAULT_BUNDLE_VERSION: &str = "cpic-v1.61.0";

/// Execution options for [`run`].
#[derive(Debug, Clone, Default)]
pub struct PgxOptions {
    pub sample: Option<String>,
    pub callable_bed: Option<String>,
    pub bam: Option<String>,
    pub reference: Option<String>,
    pub tables: Option<String>,
}

/// Output produced by [`run`].
#[derive(Debug, Clone)]
pub struct PgxRunResult {
    pub result_json_path: String,
    pub report_html_path: String,
    pub json_text: String,
}

// ---------------------------------------------------------------------------
// Minimal FIPS 180-4 SHA-256 (avoids extra crate dependencies)
// ---------------------------------------------------------------------------

const SHA256_K: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

/// Compute lowercase hex SHA-256 digest of `data`.
pub fn sha256_hex(data: &[u8]) -> String {
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
    ];
    let bit_len = (data.len() as u64).wrapping_mul(8);
    let mut padded = Vec::with_capacity(data.len() + 72);
    padded.extend_from_slice(data);
    padded.push(0x80);
    while (padded.len() % 64) != 56 {
        padded.push(0);
    }
    padded.extend_from_slice(&bit_len.to_be_bytes());

    let mut w = [0u32; 64];
    for chunk in padded.chunks_exact(64) {
        for (i, word) in chunk.chunks_exact(4).enumerate() {
            w[i] = u32::from_be_bytes([word[0], word[1], word[2], word[3]]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16].wrapping_add(s0).wrapping_add(w[i - 7]).wrapping_add(s1);
        }
        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh] = h;
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let temp1 = hh.wrapping_add(s1).wrapping_add(ch).wrapping_add(SHA256_K[i]).wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let temp2 = s0.wrapping_add(maj);
            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(temp1);
            d = c;
            c = b;
            b = a;
            a = temp1.wrapping_add(temp2);
        }
        for (slot, delta) in h.iter_mut().zip([a, b, c, d, e, f, g, hh]) {
            *slot = slot.wrapping_add(delta);
        }
    }
    let mut out = String::with_capacity(64);
    for word in h {
        use std::fmt::Write;
        let _ = write!(out, "{word:08x}");
    }
    out
}

// ---------------------------------------------------------------------------
// External JSON Table Bundle Structures & Verification
// ---------------------------------------------------------------------------

#[derive(Debug, Deserialize)]
pub struct BundleManifest {
    pub bundle_format_version: u32,
    pub engine: String,
    pub engine_version: String,
    pub version: String,
    pub data_version: String,
    pub assembly: String,
    pub files: BTreeMap<String, ManifestFileEntry>,
}

#[derive(Debug, Deserialize)]
pub struct ManifestFileEntry {
    pub sha256: String,
    pub bytes: u64,
    pub version: String,
    pub license: String,
}

#[derive(Debug, Deserialize)]
pub struct GenesBundle {
    pub assembly: String,
    pub data_version: String,
    pub genes: Vec<GeneDef>,
}

#[derive(Debug, Deserialize)]
pub struct GeneDef {
    pub gene: String,
    pub chrom: String,
    pub model: String,
    pub ploidy: String,
    pub lookup_method: String,
    pub reference_allele: String,
    pub positions: Vec<PositionDef>,
    pub alleles: Vec<AlleleDef>,
    /// "frequency": among the diplotypes that fit, the most common in the CPIC population
    /// groups is called (NAT2, as PharmCAT does); absent: the Named Allele Matcher score.
    #[serde(default)]
    pub diplotype_choice: Option<String>,
}

#[derive(Debug, Deserialize)]
pub struct PositionDef {
    pub id: usize,
    pub chrom: String,
    pub pos: u64,
    pub cpic_pos: u64,
    pub rsid: Option<String>,
    pub hgvs: Option<String>,
    pub ref_cpic: String,
    pub grch38_ref: String,
    pub is_key: bool,
    pub alt_variants: Vec<AltVariantDef>,
}

#[derive(Debug, Deserialize)]
pub struct AltVariantDef {
    pub cpic_value: String,
    pub vcf_pos: u64,
    pub vcf_ref: String,
    pub vcf_alt: String,
    pub is_ref: bool,
}

#[derive(Debug, Deserialize)]
pub struct AlleleDef {
    pub name: String,
    pub is_reference: bool,
    pub clinical_function: Option<String>,
    pub activity_value: Option<f64>,
    pub activity_label: Option<String>,
    #[serde(default)]
    pub max_pop_freq: Option<f64>,
    /// Frequency in each CPIC population group (Latino, European, East Asian, ...).
    #[serde(default)]
    pub pop_freqs: BTreeMap<String, f64>,
    pub variants: Vec<AlleleVariantDef>,
}

#[derive(Debug, Deserialize)]
pub struct AlleleVariantDef {
    pub pos_id: usize,
    pub accepted: Vec<String>,
}

#[derive(Debug, Deserialize)]
pub struct PhenotypesBundle {
    pub data_version: String,
    pub genes: BTreeMap<String, GenePhenotypeTable>,
}

#[derive(Debug, Deserialize)]
pub struct GenePhenotypeTable {
    pub lookup_method: String,
    pub rules: Vec<PhenotypeRule>,
}

#[derive(Debug, Deserialize)]
pub struct PhenotypeRule {
    pub counts: BTreeMap<String, usize>,
    pub phenotype: String,
    pub activity_score: Option<f64>,
    #[serde(default)]
    pub activity_label: Option<String>,
}

#[derive(Debug, Deserialize)]
pub struct GuidelinesBundle {
    pub data_version: String,
    pub recommendations: Vec<GuidelineRec>,
}

#[derive(Debug, Deserialize)]
pub struct GuidelineRec {
    pub drug: String,
    pub guideline_id: String,
    pub guideline_name: String,
    pub genes: Vec<String>,
    pub untested_genes: Vec<String>,
    pub lookup: BTreeMap<String, String>,
    pub population: String,
    pub classification: Option<String>,
    pub recommendation: String,
    pub url: Option<String>,
    pub source: String,
}

#[derive(Debug, Deserialize, Default)]
pub struct TraitsBundle {
    #[serde(default)]
    pub data_version: String,
    #[serde(default)]
    pub traits: Vec<TraitDef>,
}

#[derive(Debug, Deserialize)]
pub struct TraitDef {
    pub id: String,
    #[serde(rename = "trait")]
    pub name: String,
    pub model: String,
    pub source: String,
    pub evidence_level: String,
    pub positions: Vec<PositionDef>,
}

/// Verified in-memory pharmacogenetics table bundle.
pub struct Bundle {
    pub dir: PathBuf,
    pub manifest: BundleManifest,
    pub genes: GenesBundle,
    pub bam_genes: Option<BamGenesBundle>,
    pub phenotypes: PhenotypesBundle,
    pub guidelines: GuidelinesBundle,
    pub traits: TraitsBundle,
}

fn push_pgx_candidates(candidates: &mut Vec<PathBuf>, base: &Path) {
    let exact = base.join("mmgen_tools").join("pgx").join(DEFAULT_BUNDLE_VERSION);
    candidates.push(exact);
    let pgx_parent = base.join("mmgen_tools").join("pgx");
    if let Ok(entries) = std::fs::read_dir(&pgx_parent) {
        let mut dirs: Vec<PathBuf> = entries
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| p.is_dir() && p.join("manifest.json").exists())
            .collect();
        dirs.sort();
        dirs.reverse();
        for d in dirs {
            if !candidates.contains(&d) {
                candidates.push(d);
            }
        }
    }
}

/// Locate the unpacked PGx tables folder containing `manifest.json`.
pub fn locate_tables_dir(explicit: Option<&str>) -> Result<PathBuf, String> {
    if let Some(raw) = explicit {
        let p = PathBuf::from(raw);
        if p.is_file() && raw.to_ascii_lowercase().ends_with(".zip") {
            return Err(format!(
                "{raw}: PGx tables must be unpacked into a folder containing manifest.json; \
                 extract this archive into mmgen_tools/pgx/{DEFAULT_BUNDLE_VERSION}/ and pass that directory."
            ));
        }
        if !p.exists() {
            return Err(format!(
                "{raw}: PGx tables directory does not exist. Build or download tables into \
                 mmgen_tools/pgx/{DEFAULT_BUNDLE_VERSION}/ (e.g. via `python trials/python/build_pgx_bundle.py --download`), \
                 or pass tables=<dir>."
            ));
        }
        if !p.join("manifest.json").exists() {
            return Err(format!(
                "{raw}: missing manifest.json in PGx tables directory. Point tables= to the versioned folder \
                 containing manifest.json, genes.json, phenotypes.json, and guidelines.json."
            ));
        }
        return Ok(p);
    }

    if let Ok(env_dir) = std::env::var("MMGEN_PGX_TABLES") {
        let trimmed = env_dir.trim();
        if !trimmed.is_empty() {
            return locate_tables_dir(Some(trimmed));
        }
    }

    let mut candidates = Vec::new();
    if let Ok(cwd) = std::env::current_dir() {
        push_pgx_candidates(&mut candidates, &cwd);
    }
    if let Ok(exe) = std::env::current_exe() {
        for ancestor in exe.ancestors().take(4).skip(1) {
            push_pgx_candidates(&mut candidates, ancestor);
        }
    }
    for home_var in ["HOME", "USERPROFILE"] {
        if let Ok(home) = std::env::var(home_var) {
            candidates.push(Path::new(&home).join(".cache").join("mmgen").join("pgx").join(DEFAULT_BUNDLE_VERSION));
        }
    }

    for cand in &candidates {
        if cand.join("manifest.json").exists() {
            return Ok(cand.clone());
        }
    }

    let searched: Vec<String> = candidates.iter().take(4).map(|p| p.display().to_string()).collect();
    Err(format!(
        "No mmgen-pgx table bundle found (searched: {}). Set MMGEN_PGX_TABLES=<dir>, pass tables=<dir>, \
         or build/install the bundle into mmgen_tools/pgx/{DEFAULT_BUNDLE_VERSION}/.",
        searched.join(", ")
    ))
}

fn read_verified_json<T: for<'de> Deserialize<'de>>(dir: &Path, name: &str, expected: &ManifestFileEntry) -> Result<T, String> {
    let path = dir.join(name);
    let bytes = std::fs::read(&path).map_err(|e| {
        format!(
            "{}: cannot read PGx table file ({e}). Ensure the complete bundle is unpacked in {}.",
            path.display(),
            dir.display()
        )
    })?;
    if bytes.len() as u64 != expected.bytes {
        return Err(format!(
            "{}: size mismatch (expected {} bytes from manifest.json, found {} bytes). Re-download or rebuild the PGx table bundle.",
            path.display(),
            expected.bytes,
            bytes.len()
        ));
    }
    let got_sha = sha256_hex(&bytes);
    if !got_sha.eq_ignore_ascii_case(&expected.sha256) {
        return Err(format!(
            "{}: SHA-256 checksum mismatch (expected {}, got {}). The file may be modified or corrupted; rebuild or re-download the PGx table bundle.",
            path.display(),
            expected.sha256,
            got_sha
        ));
    }
    serde_json::from_slice(&bytes).map_err(|e| format!("{}: invalid JSON format ({e}). Rebuild the PGx table bundle.", path.display()))
}

/// Load and cryptographically verify the PGx JSON bundle.
pub fn load_bundle(explicit: Option<&str>) -> Result<Bundle, String> {
    let dir = locate_tables_dir(explicit)?;
    let manifest_path = dir.join("manifest.json");
    let manifest_bytes = std::fs::read(&manifest_path)
        .map_err(|e| format!("{}: cannot read manifest.json ({e})", manifest_path.display()))?;
    let manifest: BundleManifest = serde_json::from_slice(&manifest_bytes).map_err(|e| {
        format!(
            "{}: invalid manifest.json ({e}). Regenerate tables with build_pgx_bundle.py.",
            manifest_path.display()
        )
    })?;
    if manifest.bundle_format_version != 1 {
        return Err(format!(
            "{}: unsupported bundle_format_version {} (mmgen-pgx {} expects 1).",
            manifest_path.display(),
            manifest.bundle_format_version,
            ENGINE_VERSION
        ));
    }
    if manifest.engine != ENGINE_NAME {
        return Err(format!(
            "{}: manifest engine is {:?}, expected {:?}.",
            manifest_path.display(),
            manifest.engine,
            ENGINE_NAME
        ));
    }
    if manifest.assembly != "GRCh38" {
        return Err(format!(
            "{}: manifest assembly is {:?}; mmgen-pgx only supports GRCh38.",
            manifest_path.display(),
            manifest.assembly
        ));
    }
    let req = |fname: &str| -> Result<&ManifestFileEntry, String> {
        manifest.files.get(fname).ok_or_else(|| {
            format!(
                "{}: manifest.json is missing required file entry {:?}.",
                manifest_path.display(),
                fname
            )
        })
    };
    let genes: GenesBundle = read_verified_json(&dir, "genes.json", req("genes.json")?)?;
    let bam_genes: Option<BamGenesBundle> = match manifest.files.get("bam_genes.json") {
        Some(entry) => Some(read_verified_json(&dir, "bam_genes.json", entry)?),
        None => None,
    };
    let phenotypes: PhenotypesBundle = read_verified_json(&dir, "phenotypes.json", req("phenotypes.json")?)?;
    let guidelines: GuidelinesBundle = read_verified_json(&dir, "guidelines.json", req("guidelines.json")?)?;
    let traits: TraitsBundle = match manifest.files.get("traits.json") {
        Some(entry) => read_verified_json(&dir, "traits.json", entry)?,
        None => TraitsBundle::default(),
    };
    Ok(Bundle { dir, manifest, genes, bam_genes, phenotypes, guidelines, traits })
}

// ---------------------------------------------------------------------------
// Callable BED Index (O(log N) binary search over merged intervals)
// ---------------------------------------------------------------------------

pub(crate) struct CallableBed {
    chroms: HashMap<String, Vec<(u64, u64)>>,
}

impl CallableBed {
    fn load(path: &str, filter_chroms: &HashSet<String>) -> Result<Self, String> {
        if !Path::new(path).exists() {
            return Err(format!("{path}: callable BED file does not exist. Pass a valid BED path or callable_bed=None."));
        }
        let reader = externo::abre_leitura(path)?;
        let mut raw: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        for (line_no, line_res) in reader.lines().enumerate() {
            let line = line_res.map_err(|e| format!("{path}: error reading line {}: {e}", line_no + 1))?;
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') || t.starts_with("track ") || t.starts_with("browser ") {
                continue;
            }
            let mut it = t.split('\t');
            let (Some(c), Some(s), Some(e)) = (it.next(), it.next(), it.next()) else {
                return Err(format!(
                    "{path}: line {} has fewer than 3 tab-separated BED columns (expected chrom, start, end).",
                    line_no + 1
                ));
            };
            let canon = nomes::canonico(c);
            if !filter_chroms.contains(&canon) {
                continue;
            }
            let start = s.parse::<u64>().map_err(|_| format!("{path}: line {} start coordinate {:?} is not an integer.", line_no + 1, s))?;
            let end = e.parse::<u64>().map_err(|_| format!("{path}: line {} end coordinate {:?} is not an integer.", line_no + 1, e))?;
            if end < start {
                return Err(format!("{path}: line {} has end ({end}) < start ({start}).", line_no + 1));
            }
            if end > start {
                raw.entry(canon).or_default().push((start, end));
            }
        }
        for ivs in raw.values_mut() {
            ivs.sort_unstable();
            let mut merged: Vec<(u64, u64)> = Vec::with_capacity(ivs.len());
            for (s, e) in ivs.drain(..) {
                if let Some(last) = merged.last_mut() {
                    if s <= last.1 {
                        last.1 = last.1.max(e);
                        continue;
                    }
                }
                merged.push((s, e));
            }
            *ivs = merged;
        }
        Ok(Self { chroms: raw })
    }

    fn has_chrom(&self, chrom: &str) -> bool {
        self.chroms.contains_key(chrom)
    }

    fn contains(&self, chrom: &str, start0: u64, end0: u64) -> bool {
        let Some(ivs) = self.chroms.get(chrom) else { return false };
        let idx = ivs.partition_point(|&(s, _)| s <= start0);
        idx > 0 && ivs[idx - 1].1 >= end0
    }
}

// ---------------------------------------------------------------------------
// Streaming / .tbi-Indexed VCF Reader + Tan et al. Normalization
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct VcfRecord {
    pos: u64,
    ref_seq: String,
    alts: Vec<String>,
    raw_gt: String,
    phase_set: Option<String>,
}

fn check_assembly_header_line(vcf: &str, line: &str) -> Result<(), String> {
    let lower = line.to_ascii_lowercase();
    if lower.starts_with("##reference=")
        && (lower.contains("grch37") || lower.contains("hg19") || lower.contains("b37") || lower.contains("hs37d5"))
    {
        return Err(format!(
            "{vcf}: header ##reference indicates GRCh37/hg19 ({}). mmgen-pgx v1 requires GRCh38; liftover or realign to GRCh38 first.",
            line.trim()
        ));
    }
    if lower.starts_with("##contig=") {
        if lower.contains("assembly=grch37") || lower.contains("assembly=hg19") || lower.contains("assembly=b37") {
            return Err(format!(
                "{vcf}: contig header indicates GRCh37/hg19 assembly ({}). mmgen-pgx v1 requires GRCh38; liftover or realign to GRCh38 first.",
                line.trim()
            ));
        }
        if let Some(id_start) = line.find("ID=") {
            let rest = &line[id_start + 3..];
            let id = rest.split([',', '>']).next().unwrap_or("").trim();
            let canon = nomes::canonico(id);
            let is_grch37_len = (canon == "1" && (line.contains("length=249250621") || line.contains("length=249,250,621")))
                || (canon == "2" && line.contains("length=243199373"))
                || (canon == "10" && line.contains("length=135534747"))
                || (canon == "19" && line.contains("length=59128983"));
            if is_grch37_len {
                return Err(format!(
                    "{vcf}: contig {id} has GRCh37 length (expected GRCh38; e.g., chr1 is 248,956,422 bp). \
                     mmgen-pgx v1 requires GRCh38; liftover or realign to GRCh38 first."
                ));
            }
        }
    }
    Ok(())
}

fn resolve_sample_column(vcf: &str, chrom_line: &str, requested: Option<&str>) -> Result<(usize, String), String> {
    let cols: Vec<&str> = chrom_line.trim_end_matches(['\r', '\n']).split('\t').collect();
    if cols.len() < 10 {
        return Err(format!(
            "{vcf}: #CHROM header has only {} columns; at least 10 columns (including FORMAT and sample genotype) are required.",
            cols.len()
        ));
    }
    let samples = &cols[9..];
    match requested {
        Some(target) => {
            if let Some(rel) = samples.iter().position(|&s| s == target) {
                Ok((9 + rel, target.to_string()))
            } else {
                let preview: Vec<&str> = samples.iter().copied().take(8).collect();
                let suffix = if samples.len() > 8 { ", ..." } else { "" };
                Err(format!(
                    "{vcf}: sample {:?} was not found in #CHROM header (available: {}{suffix}).",
                    target,
                    preview.join(", ")
                ))
            }
        }
        None => {
            if samples.len() == 1 {
                Ok((9, samples[0].to_string()))
            } else {
                let preview: Vec<&str> = samples.iter().copied().take(8).collect();
                let suffix = if samples.len() > 8 { ", ..." } else { "" };
                Err(format!(
                    "{vcf}: VCF contains {} samples ({}{suffix}); specify `sample=<name>` explicitly.",
                    samples.len(),
                    preview.join(", ")
                ))
            }
        }
    }
}

fn build_query_windows(genes: &[GeneDef], traits: &[TraitDef]) -> HashMap<String, Vec<(u64, u64)>> {
    let mut windows: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    for positions in genes.iter().map(|g| g.positions.as_slice()).chain(traits.iter().map(|t| t.positions.as_slice())) {
        for p in positions {
            let chrom = nomes::canonico(&p.chrom);
            let p_min = p.pos.min(p.cpic_pos).saturating_sub(60).max(1);
            let p_max = p.pos.max(p.cpic_pos) + p.grch38_ref.len().max(1) as u64 + 80;
            windows.entry(chrom).or_default().push((p_min, p_max));
        }
    }
    for list in windows.values_mut() {
        list.sort_unstable();
        let mut merged: Vec<(u64, u64)> = Vec::with_capacity(list.len());
        for (s, e) in list.drain(..) {
            if let Some(last) = merged.last_mut() {
                if s <= last.1 + 2_000 {
                    last.1 = last.1.max(e);
                    continue;
                }
            }
            merged.push((s, e));
        }
        *list = merged;
    }
    windows
}

fn overlaps_window(windows: &HashMap<String, Vec<(u64, u64)>>, chrom: &str, pos: u64, ref_len: usize) -> bool {
    let Some(list) = windows.get(chrom) else { return false };
    let rec_end = pos + ref_len.max(1) as u64 + 60;
    let idx = list.partition_point(|&(s, _)| s <= rec_end);
    for &(s, e) in &list[..idx] {
        if e >= pos.saturating_sub(60) && s <= rec_end {
            return true;
        }
    }
    false
}

fn extract_gt_and_ps(vcf: &str, line: &str, sample_col: usize) -> Result<Option<(String, Option<String>)>, String> {
    let mut fields = line.trim_end_matches(['\r', '\n']).split('\t');
    let chrom = fields.next().unwrap_or("");
    let pos_s = fields.next().unwrap_or("");
    let _id = fields.next();
    let _ref = fields.next();
    let _alt = fields.next();
    let _qual = fields.next();
    let filter = fields.next().unwrap_or(".");
    let _info = fields.next();
    let fmt = fields.next().unwrap_or("");
    let sample_val = line.trim_end_matches(['\r', '\n']).split('\t').nth(sample_col).unwrap_or(".");
    if chrom.is_empty() || pos_s.is_empty() {
        return Ok(None);
    }
    if filter != "." && filter != "PASS" && !filter.is_empty() {
        return Ok(None);
    }
    let mut gt_idx = None;
    let mut ps_idx = None;
    for (idx, tag) in fmt.split(':').enumerate() {
        if tag == "GT" {
            gt_idx = Some(idx);
        } else if tag == "PS" {
            ps_idx = Some(idx);
        }
    }
    let Some(gi) = gt_idx else {
        return Err(format!("{vcf}: record at {chrom}:{pos_s} has FORMAT {:?} without GT field.", fmt));
    };
    let parts: Vec<&str> = sample_val.split(':').collect();
    let gt = parts.get(gi).copied().unwrap_or(".").to_string();
    let ps = ps_idx
        .and_then(|pi| parts.get(pi).copied())
        .filter(|&v| !v.is_empty() && v != ".")
        .map(ToString::to_string);
    Ok(Some((gt, ps)))
}

pub(crate) struct VcfObservations {
    sample: String,
    seen_chroms: HashSet<String>,
    records_by_chrom: HashMap<String, Vec<VcfRecord>>,
}

fn accept_vcf_record(
    vcf: &str,
    line: &str,
    scol: usize,
    windows: &HashMap<String, Vec<(u64, u64)>>,
    seen_chroms: &mut HashSet<String>,
    records_by_chrom: &mut HashMap<String, Vec<VcfRecord>>,
) -> Result<(), String> {
    let trimmed = line.trim_end_matches(['\r', '\n']);
    let mut it = trimmed.split('\t');
    let (Some(raw_chr), Some(pos_s), Some(_id), Some(ref_s), Some(alt_s), Some(_q), Some(filt)) =
        (it.next(), it.next(), it.next(), it.next(), it.next(), it.next(), it.next())
    else {
        return Ok(());
    };
    let canon = nomes::canonico(raw_chr);
    seen_chroms.insert(canon.clone());
    if !windows.contains_key(&canon) {
        return Ok(());
    }
    let pos = pos_s
        .parse::<u64>()
        .map_err(|_| format!("{vcf}: invalid POS {:?} on chromosome {raw_chr}.", pos_s))?;
    if filt != "." && filt != "PASS" && !filt.is_empty() {
        return Ok(());
    }
    if !overlaps_window(windows, &canon, pos, ref_s.len()) {
        return Ok(());
    }
    let Some((gt, ps)) = extract_gt_and_ps(vcf, trimmed, scol)? else {
        return Ok(());
    };
    let alts = if alt_s == "." {
        Vec::new()
    } else {
        alt_s.split(',').map(|s| s.to_ascii_uppercase()).collect()
    };
    records_by_chrom.entry(canon).or_default().push(VcfRecord {
        pos,
        ref_seq: ref_s.to_ascii_uppercase(),
        alts,
        raw_gt: gt,
        phase_set: ps,
    });
    Ok(())
}

fn read_vcf_observations(
    vcf: &str,
    sample_req: Option<&str>,
    windows: &HashMap<String, Vec<(u64, u64)>>,
) -> Result<VcfObservations, String> {
    if !Path::new(vcf).exists() {
        return Err(format!("{vcf}: VCF file does not exist."));
    }
    let mut seen_chroms = HashSet::new();
    let mut records_by_chrom: HashMap<String, Vec<VcfRecord>> = HashMap::new();
    let mut sample_col: Option<(usize, String)> = None;
    let mut header_err: Option<String> = None;

    let tbi_path = format!("{vcf}.tbi");
    if Path::new(&tbi_path).exists() {
        let (_, _, file_chroms) = indice::le_tbi(&tbi_path)?;
        let by_canon: HashMap<String, String> =
            file_chroms.iter().map(|n| (nomes::canonico(n), n.clone())).collect();
        for n in &file_chroms {
            seen_chroms.insert(nomes::canonico(n));
        }
        let mut regions = Vec::new();
        for (canon, list) in windows {
            if let Some(orig) = by_canon.get(canon) {
                for &(s, e) in list {
                    regions.push(format!("{orig}:{s}-{e}"));
                }
            }
        }
        if regions.is_empty() {
            if let Some(first) = file_chroms.first() {
                regions.push(format!("{first}:1-1"));
            } else {
                return Err(format!("{vcf}: .tbi index contains no sequences."));
            }
        }
        indice::consultar(vcf, &regions, &mut |is_hdr, line| {
            if is_hdr {
                if header_err.is_none() {
                    if let Err(e) = check_assembly_header_line(vcf, line) {
                        header_err = Some(e);
                        return Ok(());
                    }
                    if line.starts_with("#CHROM") {
                        match resolve_sample_column(vcf, line, sample_req) {
                            Ok(sc) => sample_col = Some(sc),
                            Err(e) => header_err = Some(e),
                        }
                    }
                }
                return Ok(());
            }
            if let Some(e) = &header_err {
                return Err(e.clone());
            }
            let Some((sc, _)) = sample_col else {
                return Err(format!("{vcf}: VCF records appear before #CHROM header line."));
            };
            accept_vcf_record(vcf, line, sc, windows, &mut seen_chroms, &mut records_by_chrom)
        })?;
        if let Some(e) = header_err {
            return Err(e);
        }
    } else {
        let reader = externo::abre_leitura(vcf)?;
        for (idx, line_res) in reader.lines().enumerate() {
            let line = line_res.map_err(|e| format!("{vcf}: line {}: {e}", idx + 1))?;
            if line.starts_with("##") {
                check_assembly_header_line(vcf, &line)?;
                continue;
            }
            if line.starts_with("#CHROM") {
                sample_col = Some(resolve_sample_column(vcf, &line, sample_req)?);
                continue;
            }
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let Some((sc, _)) = sample_col else {
                return Err(format!("{vcf}: VCF records appear before #CHROM header line."));
            };
            accept_vcf_record(vcf, &line, sc, windows, &mut seen_chroms, &mut records_by_chrom)?;
        }
    }

    let Some((_, sample)) = sample_col else {
        return Err(format!("{vcf}: missing #CHROM header line (not a valid VCF)."));
    };
    for list in records_by_chrom.values_mut() {
        list.sort_by_key(|r| r.pos);
    }
    Ok(VcfObservations { sample, seen_chroms, records_by_chrom })
}

/// Tan et al. (2015) variant normalization, using optional GRCh38 FASTA for left-alignment.
fn normalize_allele(
    chrom: &str,
    mut pos: u64,
    ref_in: &str,
    alt_in: &str,
    mut genome: Option<&mut referencia::Genoma>,
) -> Result<(u64, String, String), String> {
    let mut r = ref_in.as_bytes().to_vec();
    let mut a = alt_in.as_bytes().to_vec();
    if r == a || a.is_empty() || r.is_empty() || a == b"*" || a.starts_with(b"<") {
        return Ok((pos, ref_in.to_string(), alt_in.to_string()));
    }

    let fa_chrom = match genome {
        Some(ref mut g) => g.nome(chrom),
        None => None,
    };

    // Step 1: Trim matching suffix bytes; left-shift into reference genome if FASTA is loaded.
    while r.last() == a.last() {
        if r.len() == 1 || a.len() == 1 {
            let (Some(g), Some(ref seq_name)) = (genome.as_deref_mut(), fa_chrom.as_ref()) else {
                break;
            };
            if pos <= 1 {
                break;
            }
            pos -= 1;
            let base = g.trecho(seq_name, pos - 1, pos)?;
            if base.is_empty() {
                break;
            }
            r.insert(0, base[0]);
            a.insert(0, base[0]);
        }
        r.pop();
        a.pop();
    }
    // Step 2: Trim redundant common prefix bases while keeping at least 1 base in each allele.
    while r.len() >= 2 && a.len() >= 2 && r[0] == a[0] {
        r.remove(0);
        a.remove(0);
        pos += 1;
    }
    Ok((
        pos,
        String::from_utf8_lossy(&r).into_owned(),
        String::from_utf8_lossy(&a).into_owned(),
    ))
}

// ---------------------------------------------------------------------------
// Per-Position Observation Extraction & Evidence Building
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq)]
enum HapAllele {
    Cpic(String),
    NonCuratedAlt(String),
    SpanningDel,
    Missing,
}

#[derive(Debug, Clone)]
struct PosObservation {
    hap_alleles: Vec<HapAllele>,
    phased: bool,
    phase_set: Option<String>,
    evidence_gt: String,
    in_callable: Option<bool>,
    origin: &'static str, // "vcf" | "absent_in_callable" | "missing"
    ref_mismatch_note: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
struct EvidenceEntry {
    chrom: String,
    pos: u64,
    rsid: Option<String>,
    genotype: String,
    in_callable: Option<bool>,
    origin: String,
}

fn find_covering_records<'a>(recs: &'a [VcfRecord], target_pos: u64) -> Vec<&'a VcfRecord> {
    // Records sorted by pos; scan nearby window [target_pos - 60 .. target_pos + 1]
    let start_idx = recs.partition_point(|r| r.pos + 80 < target_pos);
    let mut hits = Vec::new();
    for r in &recs[start_idx..] {
        if r.pos > target_pos + 1 {
            break;
        }
        if r.pos <= target_pos && target_pos < r.pos + r.ref_seq.len().max(1) as u64 {
            hits.push(r);
        }
    }
    hits
}

fn ref_cpic_value(p: &PositionDef) -> &str {
    p.alt_variants
        .iter()
        .find(|av| av.is_ref)
        .map(|av| av.cpic_value.as_str())
        .unwrap_or(p.ref_cpic.as_str())
}

fn vcf_ref_matches_position(rec: &VcfRecord, p: &PositionDef) -> bool {
    if rec.pos > p.pos {
        return false;
    }
    let off = (p.pos - rec.pos) as usize;
    let rec_tail = &rec.ref_seq.as_bytes()[off..];
    let target = p.grch38_ref.as_bytes();
    let cmp_len = rec_tail.len().min(target.len());
    cmp_len > 0 && rec_tail[..cmp_len] == target[..cmp_len]
}

fn decode_record_at_position(
    rec: &VcfRecord,
    p: &PositionDef,
    genome: &mut Option<referencia::Genoma>,
    in_callable: Option<bool>,
) -> Result<PosObservation, String> {
    if rec.pos <= p.pos && !vcf_ref_matches_position(rec, p) {
        return Ok(PosObservation {
            hap_alleles: vec![HapAllele::Missing],
            phased: false,
            phase_set: None,
            evidence_gt: "./.".to_string(),
            in_callable,
            origin: "missing",
            ref_mismatch_note: Some(format!(
                "{}:{} ({}): VCF REF {:?} disagrees with GRCh38 expected {:?}; treated as missing.",
                p.chrom, p.pos, p.rsid.as_deref().unwrap_or("."), rec.ref_seq, p.grch38_ref
            )),
        });
    }

    let phased = rec.raw_gt.contains('|');
    let tokens: Vec<&str> = rec.raw_gt.split(['/', '|']).collect();
    let mut hap_alleles = Vec::with_capacity(tokens.len());
    let mut disp_alleles = Vec::with_capacity(tokens.len());
    let mut any_missing = false;

    for tok in tokens {
        if tok == "." || tok.is_empty() {
            any_missing = true;
            hap_alleles.push(HapAllele::Missing);
            disp_alleles.push(".".to_string());
            continue;
        }
        let Ok(idx) = tok.parse::<usize>() else {
            any_missing = true;
            hap_alleles.push(HapAllele::Missing);
            disp_alleles.push("?".to_string());
            continue;
        };
        if idx == 0 {
            let rv = ref_cpic_value(p).to_string();
            disp_alleles.push(rv.clone());
            hap_alleles.push(HapAllele::Cpic(rv));
            continue;
        }
        let Some(raw_alt) = rec.alts.get(idx - 1) else {
            any_missing = true;
            hap_alleles.push(HapAllele::Missing);
            disp_alleles.push("?".to_string());
            continue;
        };
        if raw_alt == "*" {
            hap_alleles.push(HapAllele::SpanningDel);
            disp_alleles.push("*".to_string());
            continue;
        }
        let (np, nr, na) = normalize_allele(&p.chrom, rec.pos, &rec.ref_seq, raw_alt, genome.as_mut())?;
        if let Some(av) = p
            .alt_variants
            .iter()
            .find(|av| av.vcf_pos == np && av.vcf_ref == nr && av.vcf_alt == na)
        {
            disp_alleles.push(av.cpic_value.clone());
            hap_alleles.push(HapAllele::Cpic(av.cpic_value.clone()));
        } else if rec.pos > p.pos && np > p.pos {
            let rv = ref_cpic_value(p).to_string();
            disp_alleles.push(rv.clone());
            hap_alleles.push(HapAllele::Cpic(rv));
        } else if np == p.pos || (rec.pos == p.pos && p.alt_variants.len() >= 1) {
            let label = format!("{nr}>{na}");
            disp_alleles.push(label.clone());
            hap_alleles.push(HapAllele::NonCuratedAlt(label));
        } else {
            // An upstream deletion/MNP covers this coordinate without matching ref at p.pos.
            hap_alleles.push(HapAllele::SpanningDel);
            disp_alleles.push(format!("{np}:{nr}>{na}"));
        }
    }

    let sep = if phased { "|" } else { "/" };
    let origin = if any_missing || hap_alleles.iter().any(|h| matches!(h, HapAllele::Missing)) {
        "missing"
    } else {
        "vcf"
    };
    Ok(PosObservation {
        hap_alleles,
        phased,
        phase_set: rec.phase_set.clone(),
        evidence_gt: disp_alleles.join(sep),
        in_callable,
        origin,
        ref_mismatch_note: None,
    })
}

fn observe_position_list(
    positions: &[PositionDef],
    ploidy: &str,
    vcf_obs: &VcfObservations,
    callable: Option<&CallableBed>,
    genome: &mut Option<referencia::Genoma>,
) -> Result<(Vec<PosObservation>, Vec<EvidenceEntry>, Vec<String>), String> {
    let mut obs_list = Vec::with_capacity(positions.len());
    let mut evidence = Vec::with_capacity(positions.len());
    let mut notes = Vec::new();

    // Scan matching records across positions; if FASTA is loaded, left-alignment also handles shifted indels.
    for p in positions {
        let canon = nomes::canonico(&p.chrom);
        let recs = vcf_obs.records_by_chrom.get(&canon).map(Vec::as_slice).unwrap_or(&[]);
        let seq_in_vcf = vcf_obs.seen_chroms.contains(&canon);
        let seq_in_bed = callable.is_some_and(|b| b.has_chrom(&canon));

        let end0 = p.pos - 1 + p.grch38_ref.len().max(1) as u64;
        let in_cal = match callable {
            Some(b) => b.contains(&canon, p.pos - 1, end0),
            None => false,
        };
        let in_cal_opt = callable.map(|_| in_cal);

        let mut chosen: Option<PosObservation> = None;
        let has_nonref_obs = |o: &PosObservation, pos_def: &PositionDef| -> bool {
            let rv = ref_cpic_value(pos_def);
            o.hap_alleles.iter().any(|h| match h {
                HapAllele::Cpic(v) => v != rv,
                HapAllele::NonCuratedAlt(_) | HapAllele::SpanningDel => true,
                HapAllele::Missing => false,
            })
        };
        for rec in find_covering_records(recs, p.pos) {
            let decoded = decode_record_at_position(rec, p, genome, in_cal_opt)?;
            if let Some(note) = &decoded.ref_mismatch_note {
                notes.push(note.clone());
            }
            let replace = match &chosen {
                None => true,
                Some(cur) if cur.origin == "missing" && decoded.origin != "missing" => true,
                Some(cur) if !has_nonref_obs(cur, p) && has_nonref_obs(&decoded, p) => true,
                _ => false,
            };
            if replace {
                chosen = Some(decoded);
            }
        }

        if !chosen.as_ref().is_some_and(|c| c.origin != "missing" && has_nonref_obs(c, p)) && genome.is_some() {
            let s_idx = recs.partition_point(|r| r.pos + 80 < p.pos);
            for rec in &recs[s_idx..] {
                if rec.pos > p.pos + 60 {
                    break;
                }
                if rec.pos <= p.pos && p.pos < rec.pos + rec.ref_seq.len().max(1) as u64 {
                    continue;
                }
                for raw_alt in &rec.alts {
                    let (np, _, _) = normalize_allele(&p.chrom, rec.pos, &rec.ref_seq, raw_alt, genome.as_mut())?;
                    if np == p.pos {
                        let decoded = decode_record_at_position(rec, p, genome, in_cal_opt)?;
                        let replace = match &chosen {
                            None => true,
                            Some(cur) if cur.origin == "missing" && decoded.origin != "missing" => true,
                            Some(cur) if !has_nonref_obs(cur, p) && has_nonref_obs(&decoded, p) => true,
                            _ => false,
                        };
                        if replace {
                            chosen = Some(decoded);
                        }
                        break;
                    }
                }
                if chosen.as_ref().is_some_and(|c| c.origin != "missing" && has_nonref_obs(c, p)) {
                    break;
                }
            }
        }

        let obs = match chosen {
            Some(o) => o,
            None if in_cal && (seq_in_vcf || seq_in_bed) => {
                let rv = ref_cpic_value(p).to_string();
                let hap_n = if ploidy == "haploid" { 1 } else { 2 };
                PosObservation {
                    hap_alleles: vec![HapAllele::Cpic(rv.clone()); hap_n],
                    phased: false,
                    phase_set: None,
                    evidence_gt: if hap_n == 1 { rv } else { format!("{rv}/{rv}") },
                    in_callable: Some(true),
                    origin: "absent_in_callable",
                    ref_mismatch_note: None,
                }
            }
            None => PosObservation {
                hap_alleles: vec![HapAllele::Missing],
                phased: false,
                phase_set: None,
                evidence_gt: "./.".to_string(),
                in_callable: in_cal_opt,
                origin: "missing",
                ref_mismatch_note: None,
            },
        };

        evidence.push(EvidenceEntry {
            chrom: p.chrom.clone(),
            pos: p.pos,
            rsid: p.rsid.clone(),
            genotype: obs.evidence_gt.clone(),
            in_callable: obs.in_callable,
            origin: obs.origin.to_string(),
        });
        obs_list.push(obs);
    }

    Ok((obs_list, evidence, notes))
}

fn observe_gene_positions(
    gene: &GeneDef,
    vcf_obs: &VcfObservations,
    callable: Option<&CallableBed>,
    genome: &mut Option<referencia::Genoma>,
) -> Result<(Vec<PosObservation>, Vec<EvidenceEntry>, Vec<String>), String> {
    observe_position_list(&gene.positions, &gene.ploidy, vcf_obs, callable, genome)
}

// ---------------------------------------------------------------------------
// Allele Ordering, Phenotype Lookup & Gene Callers
// ---------------------------------------------------------------------------

fn star_sort_key(name: &str, is_ref: bool) -> (u8, u64, String) {
    if is_ref || name.eq_ignore_ascii_case("Reference") || name.contains("non-responsive") {
        return (0, 0, name.to_string());
    }
    if let Some(after_star) = name.strip_prefix('*') {
        let digits: String = after_star.chars().take_while(|c| c.is_ascii_digit()).collect();
        if let Ok(num) = digits.parse::<u64>() {
            return (1, num, name.to_string());
        }
    }
    (2, u64::MAX, name.to_string())
}

fn sort_pair<'a>(a: &'a AlleleDef, b: &'a AlleleDef) -> (&'a AlleleDef, &'a AlleleDef) {
    if star_sort_key(&a.name, a.is_reference) <= star_sort_key(&b.name, b.is_reference) {
        (a, b)
    } else {
        (b, a)
    }
}

fn lookup_phenotype_and_score(
    gene: &GeneDef,
    hap_alleles: &[&AlleleDef],
    pheno_table: Option<&GenePhenotypeTable>,
) -> (Option<String>, Option<f64>) {
    let Some(table) = pheno_table else {
        return (None, None);
    };
    if table.rules.is_empty() {
        return (None, None);
    }
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for al in hap_alleles {
        let key = if gene.lookup_method == "ACTIVITY_SCORE" {
            al.activity_label.clone().unwrap_or_else(|| "n/a".to_string())
        } else {
            al.clinical_function.clone().unwrap_or_else(|| "Uncertain function".to_string())
        };
        *counts.entry(key).or_insert(0) += 1;
    }
    let (pheno, score, _) = lookup_phenotype_by_counts(&gene.lookup_method, &counts, Some(table));
    (pheno, score)
}

pub(crate) fn format_activity_score(s: f64) -> String {
    let scaled = (s * 4.0).round() as i64;
    if scaled % 2 == 0 {
        format!("{:.1}", (scaled as f64) / 4.0)
    } else {
        format!("{:.2}", (scaled as f64) / 4.0)
    }
}

pub(crate) fn lookup_phenotype_by_counts(
    lookup_method: &str,
    counts: &BTreeMap<String, usize>,
    pheno_table: Option<&GenePhenotypeTable>,
) -> (Option<String>, Option<f64>, Option<String>) {
    let _ = lookup_method;
    let Some(table) = pheno_table else {
        return (None, None, None);
    };
    for rule in &table.rules {
        if &rule.counts == counts {
            let act_label = rule
                .activity_label
                .clone()
                .or_else(|| rule.activity_score.map(format_activity_score));
            return (Some(rule.phenotype.clone()), rule.activity_score, act_label);
        }
    }
    (Some("Indeterminate".to_string()), None, None)
}

/// A diplotype that also fits the positions read but was not chosen (lower score, or the
/// cis/trans configuration not assumed). Written to `x_alternatives`.
#[derive(Debug, Clone, Serialize)]
pub(crate) struct Alternative {
    diplotype: String,
    phenotype: Option<String>,
    activity_score: Option<f64>,
    /// Named Allele Matcher score; `None` for the variant models (no scoring there).
    score: Option<usize>,
    /// Expected frequency (mean over the CPIC population groups), for genes chosen by frequency.
    #[serde(skip_serializing_if = "Option::is_none")]
    frequency: Option<f64>,
    changes_phenotype: bool,
    changes_recommendation: bool,
}

#[derive(Debug, Clone, Default)]
pub(crate) struct GeneOutcome {
    gene: String,
    called: bool,
    diplotype: Option<String>,
    phenotype: Option<String>,
    activity_score: Option<f64>,
    activity_label: Option<String>,
    copy_number: Option<f64>,
    source: String,
    haplotype_1: Option<String>,
    haplotype_2: Option<String>,
    not_called_reason: Option<String>,
    calling_method: String,
    phased: bool,
    positions_total: usize,
    positions_read: usize,
    positions_missing: usize,
    evidence: Vec<EvidenceEntry>,
    notes: Vec<String>,
    lookup_method: String,
    /// Score of the chosen diplotype (Named Allele Matcher), when scoring decided the call.
    score: Option<usize>,
    /// How the call was chosen among several that fit: "score", "frequency" or "trans".
    chosen_by: Option<&'static str>,
    /// Expected frequency of the chosen diplotype, when chosen by frequency.
    frequency: Option<f64>,
    alternatives: Vec<Alternative>,
    /// Gene-level warnings, copied to the result's `warnings` with the gene name in front.
    warnings: Vec<String>,
    missing_key_positions: Vec<String>,
}

impl GeneOutcome {
    fn uncalled(
        gene: &GeneDef,
        method: &str,
        obs: &[PosObservation],
        evidence: Vec<EvidenceEntry>,
        mut notes: Vec<String>,
        reason: String,
    ) -> Self {
        let total = gene.positions.len();
        let read = obs.iter().filter(|o| o.origin != "missing").count();
        let missing = total.saturating_sub(read);
        let phased = obs.iter().any(|o| o.phased);
        notes.retain(|s| !s.is_empty());
        Self {
            gene: gene.gene.clone(),
            called: false,
            diplotype: None,
            phenotype: None,
            activity_score: None,
            activity_label: None,
            copy_number: None,
            source: "vcf".to_string(),
            haplotype_1: None,
            haplotype_2: None,
            not_called_reason: Some(reason),
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
}

fn format_pos_list(gene: &GeneDef, missing_idxs: &[usize]) -> String {
    missing_idxs
        .iter()
        .take(5)
        .map(|&i| {
            let p = &gene.positions[i];
            match &p.rsid {
                Some(rs) => format!("{rs} ({}:{})", p.chrom, p.pos),
                None => format!("{}:{}", p.chrom, p.pos),
            }
        })
        .collect::<Vec<_>>()
        .join(", ")
        + if missing_idxs.len() > 5 { ", ..." } else { "" }
}

/// Evaluate one gene using exact haplotype/diplotype matching (`star_allele` or `variant_genotype`).
fn call_combinatorial_gene(
    gene: &GeneDef,
    obs: &[PosObservation],
    evidence: Vec<EvidenceEntry>,
    notes: Vec<String>,
    pheno_table: Option<&GenePhenotypeTable>,
    method: &str,
    haploid: bool,
) -> GeneOutcome {
    let total = gene.positions.len();
    let read = obs.iter().filter(|o| o.origin != "missing").count();
    let missing = total.saturating_sub(read);

    // Check for any unrecognized non-CPIC allele or spanning deletion at an observed position.
    let mut unrec = Vec::new();
    for (idx, o) in obs.iter().enumerate() {
        for ha in &o.hap_alleles {
            if let HapAllele::NonCuratedAlt(a) = ha {
                let p = &gene.positions[idx];
                unrec.push(format!("{}:{} ({}) ALT {a}", p.chrom, p.pos, p.rsid.as_deref().unwrap_or(".")));
            } else if matches!(ha, HapAllele::SpanningDel) {
                let p = &gene.positions[idx];
                unrec.push(format!("{}:{} (*) spanning deletion", p.chrom, p.pos));
            }
        }
    }
    if !unrec.is_empty() {
        return GeneOutcome::uncalled(
            gene,
            method,
            obs,
            evidence,
            notes,
            format!("Unrecognized allele(s) outside curated CPIC definitions: {}", unrec.join("; ")),
        );
    }

    // Build allele lookup matrix: accepted sets per allele (Vec<&[String]> indexed by pos_id via sparse map).
    let mut var_maps: Vec<HashMap<usize, &[String]>> = Vec::with_capacity(gene.alleles.len());
    for al in &gene.alleles {
        let m: HashMap<usize, &[String]> = al.variants.iter().map(|v| (v.pos_id, v.accepted.as_slice())).collect();
        var_maps.push(m);
    }
    let ref_vals: Vec<&str> = gene.positions.iter().map(ref_cpic_value).collect();

    // Determine whether the sample has any observed non-reference CPIC allele at each position.
    let accepts = |m: &HashMap<usize, &[String]>, pid: usize, val: &str| -> bool {
        match m.get(&pid) {
            Some(list) => list.iter().any(|s| s == val),
            None => val == ref_vals[pid],
        }
    };

    // Single-allele filter + phantom check (an allele with missing defining variants must explain >= 1 observed non-ref allele).
    let mut candidates: Vec<usize> = Vec::new();
    for (a_idx, m) in var_maps.iter().enumerate() {
        let has_missing_def = m.keys().any(|&pid| obs[pid].origin == "missing");
        let mut explains_observed_nonref = false;
        let mut compatible = true;
        for (pid, o) in obs.iter().enumerate() {
            if o.origin == "missing" {
                continue;
            }
            let mut hit_any = false;
            for ha in &o.hap_alleles {
                if let HapAllele::Cpic(val) = ha {
                    if accepts(m, pid, val) {
                        hit_any = true;
                        if val != ref_vals[pid]
                            && m.get(&pid).is_some_and(|list| !list.iter().any(|s| s == ref_vals[pid]))
                        {
                            explains_observed_nonref = true;
                        }
                    }
                }
            }
            if !hit_any {
                compatible = false;
                break;
            }
        }
        if compatible && (!has_missing_def || explains_observed_nonref) {
            candidates.push(a_idx);
        }
    }

    // Enumerate matching haplotypes (if haploid) or diplotypes (if diploid) respecting phase sets.
    let mut matches: Vec<(usize, Option<usize>)> = Vec::new();
    if haploid {
        for &ai in &candidates {
            let ok = obs.iter().enumerate().all(|(pid, o)| match &o.hap_alleles[..] {
                [HapAllele::Missing] => true,
                [HapAllele::Cpic(v), ..] => accepts(&var_maps[ai], pid, v),
                _ => false,
            });
            if ok {
                matches.push((ai, None));
            }
        }
    } else {
        for ii in 0..candidates.len() {
            let ai = candidates[ii];
            let m1 = &var_maps[ai];
            for &aj in &candidates[ii..] {
                let m2 = &var_maps[aj];
                let mut ok = true;
                let mut block_masks: HashMap<&str, u8> = HashMap::new();
                for (pid, o) in obs.iter().enumerate() {
                    if o.origin == "missing" {
                        continue;
                    }
                    let v0 = match &o.hap_alleles[0] {
                        HapAllele::Cpic(v) => v.as_str(),
                        _ => {
                            ok = false;
                            break;
                        }
                    };
                    let v1 = match o.hap_alleles.get(1) {
                        Some(HapAllele::Cpic(v)) => v.as_str(),
                        Some(_) => {
                            ok = false;
                            break;
                        }
                        None => v0,
                    };
                    let o1 = accepts(m1, pid, v0) && accepts(m2, pid, v1);
                    let o2 = accepts(m1, pid, v1) && accepts(m2, pid, v0);
                    if !o1 && !o2 {
                        ok = false;
                        break;
                    }
                    if o.phased && v0 != v1 {
                        let blk = o.phase_set.as_deref().unwrap_or("__GLOBAL__");
                        let cur = block_masks.entry(blk).or_insert(3);
                        *cur &= (u8::from(o1)) | (u8::from(o2) << 1);
                        if *cur == 0 {
                            ok = false;
                            break;
                        }
                    }
                }
                if ok {
                    matches.push((ai, Some(aj)));
                }
            }
        }
    }

    if matches.is_empty() {
        return GeneOutcome::uncalled(
            gene,
            method,
            obs,
            evidence,
            notes,
            format!("Observed {} genotypes do not match any curated CPIC allele/diplotype combination.", gene.gene),
        );
    }

    // Every diplotype that fits the positions read, scored as PharmCAT's Named Allele Matcher
    // does: each allele is worth the positions read that define it, the diplotype the sum.
    struct Cand {
        dip: String,
        pheno: Option<String>,
        act: Option<f64>,
        h1: Option<String>,
        h2: Option<String>,
        sk: (u8, u64, String),
        score: usize,
        freq: Option<f64>,
    }
    let read_mask: Vec<bool> = obs.iter().map(|o| o.origin != "missing").collect();
    let mut cands: Vec<Cand> = Vec::new();
    let mut seen_dips = BTreeSet::new();
    for &(ai, aj_opt) in &matches {
        let (a1, a2_opt) = match aj_opt {
            Some(aj) => {
                let (x, y) = sort_pair(&gene.alleles[ai], &gene.alleles[aj]);
                (x, Some(y))
            }
            None => (&gene.alleles[ai], None),
        };
        let dip = match a2_opt {
            Some(a2) => format!("{}/{}", a1.name, a2.name),
            None => a1.name.clone(),
        };
        if !seen_dips.insert(dip.clone()) {
            continue;
        }
        let haps: Vec<&AlleleDef> = match a2_opt {
            Some(a2) => vec![a1, a2],
            None => vec![a1],
        };
        let (pheno, act) = lookup_phenotype_and_score(gene, &haps, pheno_table);
        let score = haps.iter().map(|a| allele_score(a, &read_mask, read)).sum();
        let freq = diplotype_frequency(a1, a2_opt);
        let sk = star_sort_key(&a1.name, a1.is_reference);
        cands.push(Cand { dip, pheno, act, h1: Some(a1.name.clone()), h2: a2_opt.map(|x| x.name.clone()), sk, score, freq });
    }
    // NAT2 (diplotype_choice "frequency"): the most common diplotype in the CPIC population
    // groups wins, the score breaking ties; without frequencies for any candidate, the score.
    let by_frequency = gene.diplotype_choice.as_deref() == Some("frequency")
        && cands.len() > 1
        && cands.iter().any(|c| c.freq.is_some_and(|f| f > 0.0));
    let fkey = |c: &Cand| c.freq.unwrap_or(-1.0);
    if by_frequency {
        cands.sort_by(|x, y| {
            fkey(y)
                .partial_cmp(&fkey(x))
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| y.score.cmp(&x.score))
                .then_with(|| x.sk.cmp(&y.sk))
                .then_with(|| x.dip.cmp(&y.dip))
        });
    } else {
        cands.sort_by(|x, y| y.score.cmp(&x.score).then_with(|| x.sk.cmp(&y.sk)).then_with(|| x.dip.cmp(&y.dip)));
    }
    let top_score = cands[0].score;
    let top_freq = cands[0].freq;
    let n_top = cands
        .iter()
        .take_while(|c| c.score == top_score && (!by_frequency || fkey(c) == fkey(&cands[0])))
        .count();
    let (top, rest) = cands.split_at(n_top);

    // Ties at the top are all reported (" OR "). When they disagree in phenotype or activity
    // score there is no single phenotype: the gene is still called, without a recommendation.
    let same_clinical = top.iter().all(|c| c.pheno == top[0].pheno && c.act == top[0].act);
    let (phenotype, activity_score) = if same_clinical {
        (top[0].pheno.clone(), top[0].act)
    } else {
        let mut ph: Vec<String> = Vec::new();
        for c in top {
            let p = c.pheno.clone().unwrap_or_else(|| "Indeterminate".to_string());
            if !ph.contains(&p) {
                ph.push(p);
            }
        }
        (Some(ph.join(" OR ")), None)
    };
    let dip_joined = top.iter().map(|c| c.dip.as_str()).collect::<Vec<_>>().join(" OR ");
    let (h1, h2) = if top.len() == 1 { (top[0].h1.clone(), top[0].h2.clone()) } else { (None, None) };
    let phased_flag = obs.iter().filter(|o| o.origin == "vcf").all(|o| o.phased || o.hap_alleles.windows(2).all(|w| w[0] == w[1]));

    let describe = |c: &Cand| describe_call(&c.dip, c.pheno.as_deref(), c.act, &gene.lookup_method);
    let mut notes = notes;
    let mut warnings = Vec::new();
    if !same_clinical {
        let list = top.iter().map(describe).collect::<Vec<_>>().join("; ");
        notes.push(format!(
            "The highest-scoring diplotypes tie and differ in phenotype or activity score: {list}. No recommendation is given for {}.",
            gene.gene
        ));
        warnings.push(format!(
            "tied diplotypes give different phenotypes or activity scores ({list}); no recommendation was given. Phased genotypes would settle it."
        ));
    }
    if !rest.is_empty() && by_frequency {
        let best_by_score = cands.iter().max_by(|x, y| x.score.cmp(&y.score).then_with(|| y.sk.cmp(&x.sk))).map(|c| (c.dip.clone(), c.score));
        let by_score_note = match best_by_score {
            Some((d, s)) if s > top_score => format!(" By score alone, {d} (score {s}) would be chosen."),
            _ => String::new(),
        };
        notes.push(format!(
            "Chosen by population frequency, as PharmCAT does for {}: {dip_joined} is the most common of the diplotypes that fit \
             (expected frequency {}, the mean over the CPIC population groups; the next candidate {}).{by_score_note}",
            gene.gene,
            format_freq(top_freq),
            format_freq(rest[0].freq)
        ));
    } else if !rest.is_empty() {
        notes.push(format!(
            "Chosen by score, as PharmCAT's Named Allele Matcher does (each allele counts the positions read that define it; \
             the reference allele and *1 count every position read): {dip_joined} scores {top_score}, the next candidate {}.",
            rest[0].score
        ));
    }
    if !rest.is_empty() {
        let shown: Vec<String> = rest.iter().take(MAX_ALTERNATIVES_IN_NOTE).map(describe).collect();
        let more = if rest.len() > MAX_ALTERNATIVES_IN_NOTE {
            format!("; and {} more (see x_alternatives)", rest.len() - MAX_ALTERNATIVES_IN_NOTE)
        } else {
            String::new()
        };
        notes.push(format!("Also possible with the positions read: {}{more}", shown.join("; ")));
    }
    let alternatives: Vec<Alternative> = rest
        .iter()
        .take(MAX_ALTERNATIVES)
        .map(|c| Alternative {
            diplotype: c.dip.clone(),
            phenotype: c.pheno.clone(),
            activity_score: c.act,
            score: Some(c.score),
            frequency: if by_frequency { c.freq } else { None },
            changes_phenotype: !top.iter().any(|t| t.pheno == c.pheno && t.act == c.act),
            changes_recommendation: false, // filled in by flag_alternatives(), which sees the guidelines
        })
        .collect();

    GeneOutcome {
        gene: gene.gene.clone(),
        called: true,
        diplotype: Some(dip_joined),
        phenotype,
        activity_score,
        haplotype_1: h1,
        haplotype_2: h2,
        not_called_reason: None,
        calling_method: method.to_string(),
        phased: phased_flag,
        positions_total: total,
        positions_read: read,
        positions_missing: missing,
        evidence,
        notes,
        lookup_method: gene.lookup_method.clone(),
        source: "vcf".to_string(),
        copy_number: None,
        activity_label: None,
        score: Some(top_score),
        chosen_by: match (rest.is_empty(), by_frequency) {
            (true, _) => None,
            (false, true) => Some("frequency"),
            (false, false) => Some("score"),
        },
        frequency: if by_frequency { top_freq } else { None },
        alternatives,
        warnings,
        missing_key_positions: Vec::new(),
    }
}

/// Expected frequency of a diplotype: in each CPIC population group with a frequency for
/// both alleles, 2·f(a)·f(b) (f(a)² when a = b; f(a) for a haploid call), then the mean over
/// those groups. `None` when no group has both.
fn diplotype_frequency(a: &AlleleDef, b: Option<&AlleleDef>) -> Option<f64> {
    let mut sum = 0.0;
    let mut n = 0usize;
    for (pop, fa) in &a.pop_freqs {
        let v = match b {
            None => *fa,
            Some(b) => match b.pop_freqs.get(pop) {
                Some(fb) if a.name == b.name => fa * fb,
                Some(fb) => 2.0 * fa * fb,
                None => continue,
            },
        };
        sum += v;
        n += 1;
    }
    (n > 0).then(|| sum / n as f64)
}

fn format_freq(f: Option<f64>) -> String {
    match f {
        Some(v) if v >= 0.001 => format!("{v:.3}"),
        Some(v) => format!("{v:.1e}"),
        None => "unknown".to_string(),
    }
}

/// How many alternatives go into `x_alternatives`, and how many are spelled out in the note.
const MAX_ALTERNATIVES: usize = 50;
const MAX_ALTERNATIVES_IN_NOTE: usize = 8;

/// Named Allele Matcher score of one allele: the number of positions read that define it.
///
/// The reference allele is defined at every position, so it is worth every position read
/// (PharmCAT: "the reference allele will always have the maximum score"). So is `*1`, the
/// PharmVar reference for function: ClinPGx, the source of PharmCAT's definitions, writes
/// it out at every position. This matters only in CYP2C19, where the GRCh38 sequence is
/// `*38` and the CPIC table lists `*1` by its one difference (rs3758581); scored by that
/// one position, `*2/*35` would outscore `*1/*2` in a sample such as NA12878.
fn allele_score(al: &AlleleDef, read_mask: &[bool], n_read: usize) -> usize {
    if al.is_reference || al.name == "*1" {
        n_read
    } else {
        al.variants.iter().filter(|v| read_mask.get(v.pos_id).copied().unwrap_or(false)).count()
    }
}

/// "*2/*35 (Poor Metabolizer)", with the activity score when the gene is looked up by it.
fn describe_call(dip: &str, pheno: Option<&str>, act: Option<f64>, lookup_method: &str) -> String {
    let p = pheno.unwrap_or("no phenotype in the tables");
    match act {
        Some(a) if lookup_method == "ACTIVITY_SCORE" => format!("{dip} ({p}, activity score {})", format_activity_score(a)),
        _ => format!("{dip} ({p})"),
    }
}

/// Evaluate `variant_activity` (`DPYD`) or `variant_status` (`CFTR`, `RYR1`, `CACNA1S`, `MT-RNR1`).
fn call_variant_model_gene(
    gene: &GeneDef,
    obs: &[PosObservation],
    evidence: Vec<EvidenceEntry>,
    mut notes: Vec<String>,
    pheno_table: Option<&GenePhenotypeTable>,
    method: &str,
) -> GeneOutcome {
    let total = gene.positions.len();
    let read = obs.iter().filter(|o| o.origin != "missing").count();
    let missing = total.saturating_sub(read);
    let haploid = gene.ploidy == "haploid";

    // Extract active non-reference CPIC variants per haplotype slot or unphased copy.
    #[derive(Clone)]
    struct VarCopy<'a> {
        allele: &'a AlleleDef,
        phased: bool,
        phase_set: Option<String>,
        strand: usize, // 0 or 1 when phased
    }

    // Map each observed non-reference position -> set of non-ref values seen at that position (with strand info).
    let mut nonref_obs_at_pos: HashMap<usize, Vec<(String, usize, bool, Option<String>)>> = HashMap::new();
    for (pid, o) in obs.iter().enumerate() {
        if o.origin == "missing" {
            continue;
        }
        let rv = ref_cpic_value(&gene.positions[pid]);
        for (slot, ha) in o.hap_alleles.iter().enumerate() {
            match ha {
                HapAllele::Cpic(val) if val != rv => {
                    nonref_obs_at_pos
                        .entry(pid)
                        .or_default()
                        .push((val.clone(), slot, o.phased, o.phase_set.clone()));
                }
                HapAllele::NonCuratedAlt(alt) => {
                    let p = &gene.positions[pid];
                    notes.push(format!(
                        "{}:{} ({}): non-CPIC alternate allele {alt} does not match any curated {} variant.",
                        p.chrom,
                        p.pos,
                        p.rsid.as_deref().unwrap_or("."),
                        gene.gene
                    ));
                }
                _ => {}
            }
        }
    }

    // For haploid MT-RNR1, reject heteroplasmy (het calls at a variant site).
    if haploid {
        for (&pid, hits) in &nonref_obs_at_pos {
            if obs[pid].hap_alleles.len() > 1 && obs[pid].hap_alleles.iter().any(|h| matches!(h, HapAllele::Cpic(v) if v == ref_cpic_value(&gene.positions[pid]))) {
                return GeneOutcome::uncalled(
                    gene,
                    method,
                    obs,
                    evidence,
                    notes,
                    format!(
                        "{}:{} ({}): mitochondrial heteroplasmy observed; haploid MT-RNR1 calling requires homoplasmy.",
                        gene.positions[pid].chrom,
                        gene.positions[pid].pos,
                        gene.positions[pid].rsid.as_deref().unwrap_or(".")
                    ),
                );
            }
            let _ = hits;
        }
    }

    // Greedy multi-variant matching (e.g. DPYD HapB3, CFTR F508C,S1251N) then single-variant alleles.
    let mut multi_then_single: Vec<&AlleleDef> = gene.alleles.iter().filter(|a| !a.is_reference).collect();
    multi_then_single.sort_by(|a, b| b.variants.len().cmp(&a.variants.len()).then_with(|| a.name.cmp(&b.name)));

    let mut remaining = nonref_obs_at_pos;
    let mut var_copies: Vec<VarCopy<'_>> = Vec::new();
    for al in multi_then_single {
        if al.variants.is_empty() {
            continue;
        }
        loop {
            let can_take = al.variants.iter().all(|v| {
                remaining
                    .get(&v.pos_id)
                    .is_some_and(|hits| hits.iter().any(|h| v.accepted.iter().any(|ac| ac == &h.0)))
            });
            if !can_take {
                break;
            }
            let mut phased_all = true;
            let mut blk: Option<String> = None;
            let mut strand = 0usize;
            for (vi, v) in al.variants.iter().enumerate() {
                let hits = remaining.get_mut(&v.pos_id).unwrap();
                let pos_idx = hits.iter().position(|h| v.accepted.iter().any(|ac| ac == &h.0)).unwrap();
                let (_, s, ph, ps) = hits.remove(pos_idx);
                if vi == 0 {
                    phased_all = ph;
                    blk = ps;
                    strand = s;
                } else if !ph || s != strand {
                    phased_all = false;
                }
            }
            var_copies.push(VarCopy { allele: al, phased: phased_all, phase_set: blk, strand });
            if haploid {
                break;
            }
        }
    }

    let ref_allele = gene.alleles.iter().find(|a| a.is_reference).unwrap_or(&gene.alleles[0]);

    if var_copies.is_empty() {
        let (pheno, act_score) = if haploid {
            lookup_phenotype_and_score(gene, &[ref_allele], pheno_table)
        } else {
            lookup_phenotype_and_score(gene, &[ref_allele, ref_allele], pheno_table)
        };
        let dip = if haploid {
            ref_allele.name.clone()
        } else {
            format!("{}/{}", ref_allele.name, ref_allele.name)
        };
        return GeneOutcome {
            gene: gene.gene.clone(),
            called: true,
            diplotype: Some(dip),
            phenotype: pheno,
            activity_score: act_score,
            activity_label: act_score.map(format_activity_score),
            copy_number: None,
            source: "vcf".to_string(),
            haplotype_1: Some(ref_allele.name.clone()),
            haplotype_2: if haploid { None } else { Some(ref_allele.name.clone()) },
            not_called_reason: None,
            calling_method: method.to_string(),
            phased: true,
            positions_total: total,
            positions_read: read,
            positions_missing: missing,
            evidence,
            notes,
            lookup_method: gene.lookup_method.clone(),
            ..Default::default()
        };
    }

    if haploid {
        if var_copies.len() > 1 {
            return GeneOutcome::uncalled(
                gene,
                method,
                obs,
                evidence,
                notes,
                format!("Multiple non-reference variants found in haploid {}: cannot resolve single allele.", gene.gene),
            );
        }
        let v = var_copies[0].allele;
        let (pheno, act_score) = lookup_phenotype_and_score(gene, &[v], pheno_table);
        return GeneOutcome {
            gene: gene.gene.clone(),
            called: true,
            diplotype: Some(v.name.clone()),
            phenotype: pheno,
            activity_score: act_score,
            activity_label: act_score.map(format_activity_score),
            copy_number: None,
            source: "vcf".to_string(),
            haplotype_1: Some(v.name.clone()),
            haplotype_2: None,
            not_called_reason: None,
            calling_method: method.to_string(),
            phased: true,
            positions_total: total,
            positions_read: read,
            positions_missing: missing,
            evidence,
            notes,
            lookup_method: gene.lookup_method.clone(),
            ..Default::default()
        };
    }

    // Diploid variant_activity / variant_status with 1 or >=2 variant copies:
    if var_copies.len() == 1 {
        let v = var_copies[0].allele;
        let (a1, a2) = sort_pair(ref_allele, v);
        let (pheno, act_score) = lookup_phenotype_and_score(gene, &[a1, a2], pheno_table);
        return GeneOutcome {
            gene: gene.gene.clone(),
            called: true,
            diplotype: Some(format!("{}/{}", a1.name, a2.name)),
            phenotype: pheno,
            activity_score: act_score,
            activity_label: act_score.map(format_activity_score),
            copy_number: None,
            source: "vcf".to_string(),
            haplotype_1: Some(a1.name.clone()),
            haplotype_2: Some(a2.name.clone()),
            not_called_reason: None,
            calling_method: method.to_string(),
            phased: true,
            positions_total: total,
            positions_read: read,
            positions_missing: missing,
            evidence,
            notes,
            lookup_method: gene.lookup_method.clone(),
            ..Default::default()
        };
    }

    // >= 2 variant copies: check phasing or evaluate cis vs trans equivalence.
    let single_block = var_copies.iter().all(|c| c.phased)
        && var_copies.windows(2).all(|w| w[0].phase_set == w[1].phase_set);

    let hap_label = |vs: &[&AlleleDef]| -> String {
        if gene.lookup_method == "ACTIVITY_SCORE" {
            let mut min_act = f64::INFINITY;
            for a in vs {
                if let Some(val) = a.activity_value {
                    if val < min_act {
                        min_act = val;
                    }
                } else {
                    return "n/a".to_string();
                }
            }
            if min_act == 0.0 {
                "0.0".to_string()
            } else if (min_act - 0.5).abs() < 1e-6 {
                "0.5".to_string()
            } else {
                "1.0".to_string()
            }
        } else {
            // For status genes: MH associated > Normal function; ivacaftor responsive > non-responsive
            if vs.iter().any(|a| {
                a.clinical_function.as_deref().is_some_and(|f| {
                    f.contains("Malignant Hyperthermia") || f == "ivacaftor responsive"
                })
            }) {
                if vs.iter().any(|a| a.clinical_function.as_deref() == Some("ivacaftor responsive")) {
                    "ivacaftor responsive".to_string()
                } else {
                    "Malignant Hyperthermia associated".to_string()
                }
            } else if vs.iter().any(|a| a.clinical_function.as_deref() == Some("Uncertain function")) {
                "Uncertain function".to_string()
            } else {
                ref_allele.clinical_function.as_deref().unwrap_or("Normal function").to_string()
            }
        }
    };

    let eval_counts = |strand0: &[&AlleleDef], strand1: &[&AlleleDef]| -> (Option<String>, Option<f64>) {
        let l0 = if strand0.is_empty() {
            if gene.lookup_method == "ACTIVITY_SCORE" {
                ref_allele.activity_label.as_deref().unwrap_or("1.0").to_string()
            } else {
                ref_allele.clinical_function.as_deref().unwrap_or("Normal function").to_string()
            }
        } else {
            hap_label(strand0)
        };
        let l1 = if strand1.is_empty() {
            if gene.lookup_method == "ACTIVITY_SCORE" {
                ref_allele.activity_label.as_deref().unwrap_or("1.0").to_string()
            } else {
                ref_allele.clinical_function.as_deref().unwrap_or("Normal function").to_string()
            }
        } else {
            hap_label(strand1)
        };
        let mut counts = BTreeMap::new();
        *counts.entry(l0).or_insert(0usize) += 1;
        *counts.entry(l1).or_insert(0usize) += 1;
        if let Some(tab) = pheno_table {
            for r in &tab.rules {
                if r.counts == counts {
                    return (Some(r.phenotype.clone()), r.activity_score);
                }
            }
        }
        (Some("Indeterminate".to_string()), None)
    };

    let fmt_s = |vs: &[&AlleleDef]| -> String {
        if vs.is_empty() {
            ref_allele.name.clone()
        } else if vs.len() == 1 {
            vs[0].name.clone()
        } else {
            let mut names: Vec<&str> = vs.iter().map(|a| a.name.as_str()).collect();
            names.sort_unstable();
            names.join("+")
        }
    };

    if single_block {
        let s0: Vec<&AlleleDef> = var_copies.iter().filter(|c| c.strand == 0).map(|c| c.allele).collect();
        let s1: Vec<&AlleleDef> = var_copies.iter().filter(|c| c.strand == 1).map(|c| c.allele).collect();
        let n0 = fmt_s(&s0);
        let n1 = fmt_s(&s1);
        let (h1, h2) = if n0 <= n1 { (n0, n1) } else { (n1, n0) };
        let (pheno, act_score) = eval_counts(&s0, &s1);
        return GeneOutcome {
            gene: gene.gene.clone(),
            called: true,
            diplotype: Some(format!("{h1}/{h2}")),
            phenotype: pheno,
            activity_score: act_score,
            activity_label: act_score.map(format_activity_score),
            copy_number: None,
            source: "vcf".to_string(),
            haplotype_1: Some(h1),
            haplotype_2: Some(h2),
            not_called_reason: None,
            calling_method: method.to_string(),
            phased: true,
            positions_total: total,
            positions_read: read,
            positions_missing: missing,
            evidence,
            notes,
            lookup_method: gene.lookup_method.clone(),
            ..Default::default()
        };
    }

    // Unphased >= 2 variants: check all assignments into two haplotypes (feasible for <= 8 variants; typical 2).
    if var_copies.len() <= 8 {
        let split = |mask: u32| -> (Vec<&AlleleDef>, Vec<&AlleleDef>) {
            let mut s0 = Vec::new();
            let mut s1 = Vec::new();
            for (i, c) in var_copies.iter().enumerate() {
                if (mask >> i) & 1 == 0 {
                    s0.push(c.allele);
                } else {
                    s1.push(c.allele);
                }
            }
            (s0, s1)
        };
        let mut results = Vec::new();
        for mask in 0u32..(1u32 << var_copies.len()) {
            let (s0, s1) = split(mask);
            results.push(eval_counts(&s0, &s1));
        }
        let first = &results[0];
        let all_agree = results.iter().all(|r| r == first);
        if !all_agree {
            // cis and trans disagree. Activity-score genes (DPYD) follow CPIC and PharmCAT: the
            // variants are taken to be on different copies, so that the two lowest-function
            // ones both count (the lowest activity score of all assignments). Status genes
            // take the plain trans assignment. The other outcome goes to the alternatives.
            let chosen_mask = if gene.lookup_method == "ACTIVITY_SCORE" {
                let key = |m: u32| results[m as usize].1.unwrap_or(f64::INFINITY);
                (0u32..(1u32 << var_copies.len()))
                    .min_by(|&a, &b| key(a).partial_cmp(&key(b)).unwrap_or(std::cmp::Ordering::Equal).then(a.cmp(&b)))
                    .unwrap_or(0)
            } else {
                (0..var_copies.len()).filter(|i| i % 2 == 1).fold(0u32, |m, i| m | (1 << i))
            };
            let label_of = |mask: u32| -> (String, String) {
                let (s0, s1) = split(mask);
                let (n0, n1) = (fmt_s(&s0), fmt_s(&s1));
                if n0 <= n1 { (n0, n1) } else { (n1, n0) }
            };
            let (h1, h2) = label_of(chosen_mask);
            let (pheno, act_score) = results[chosen_mask as usize].clone();
            let mut alternatives: Vec<Alternative> = Vec::new();
            for mask in 0u32..(1u32 << var_copies.len()) {
                let r = &results[mask as usize];
                if r.0 == pheno && r.1 == act_score {
                    continue;
                }
                if alternatives.iter().any(|a| a.phenotype == r.0 && a.activity_score == r.1) {
                    continue;
                }
                let (a1, a2) = label_of(mask);
                alternatives.push(Alternative {
                    diplotype: format!("{a1}/{a2}"),
                    phenotype: r.0.clone(),
                    activity_score: r.1,
                    score: None,
                    frequency: None,
                    changes_phenotype: true,
                    changes_recommendation: false,
                });
            }
            let vnames: Vec<&str> = var_copies.iter().map(|c| c.allele.name.as_str()).collect();
            let mut notes = notes;
            notes.push(if gene.lookup_method == "ACTIVITY_SCORE" {
                format!(
                    "Unphased: {} were taken to be on different copies (in trans), so that the lowest-function variants all count, \
                     as CPIC and PharmCAT do; on the same copy (in cis) the activity score would differ.",
                    vnames.join(" and ")
                )
            } else {
                format!("Unphased: {} were taken to be on different copies (in trans).", vnames.join(" and "))
            });
            let shown: Vec<String> = alternatives
                .iter()
                .map(|a| describe_call(&a.diplotype, a.phenotype.as_deref(), a.activity_score, &gene.lookup_method))
                .collect();
            notes.push(format!("Also possible with the positions read: {}", shown.join("; ")));
            return GeneOutcome {
                gene: gene.gene.clone(),
                called: true,
                diplotype: Some(format!("{h1}/{h2}")),
                phenotype: pheno,
                activity_score: act_score,
                activity_label: act_score.map(format_activity_score),
                copy_number: None,
                source: "vcf".to_string(),
                haplotype_1: Some(h1),
                haplotype_2: Some(h2),
                not_called_reason: None,
                calling_method: method.to_string(),
                phased: false,
                positions_total: total,
                positions_read: read,
                positions_missing: missing,
                evidence,
                notes,
                lookup_method: gene.lookup_method.clone(),
                alternatives,
                chosen_by: Some("trans"),
                ..Default::default()
            };
        }
        let mut names: Vec<&str> = var_copies.iter().map(|c| c.allele.name.as_str()).collect();
        names.sort_unstable();
        let dip = if var_copies.len() == 2 && var_copies[0].allele.name == var_copies[1].allele.name {
            format!("{}/{}", names[0], names[1])
        } else {
            names.join(" AND ")
        };
        let (h1, h2) = if dip.contains('/') {
            (Some(names[0].to_string()), Some(names[1].to_string()))
        } else {
            (None, None)
        };
        return GeneOutcome {
            gene: gene.gene.clone(),
            called: true,
            diplotype: Some(dip),
            phenotype: first.0.clone(),
            activity_score: first.1,
            activity_label: first.1.map(format_activity_score),
            copy_number: None,
            source: "vcf".to_string(),
            haplotype_1: h1,
            haplotype_2: h2,
            not_called_reason: None,
            calling_method: method.to_string(),
            phased: false,
            positions_total: total,
            positions_read: read,
            positions_missing: missing,
            evidence,
            notes,
            lookup_method: gene.lookup_method.clone(),
            ..Default::default()
        };
    }

    GeneOutcome::uncalled(
        gene,
        method,
        obs,
        evidence,
        notes,
        format!("Too many unphased {} variants ({}) to resolve safely.", gene.gene, var_copies.len()),
    )
}

fn evaluate_gene(
    gene: &GeneDef,
    vcf_obs: &VcfObservations,
    callable: Option<&CallableBed>,
    pheno_table: Option<&GenePhenotypeTable>,
    genome: &mut Option<referencia::Genoma>,
) -> Result<GeneOutcome, String> {
    let (obs, evidence, notes) = observe_gene_positions(gene, vcf_obs, callable, genome)?;
    let method = format!("vcf_{}", gene.model);
    let total = gene.positions.len();
    let read = obs.iter().filter(|o| o.origin != "missing").count();

    if read == 0 {
        let canon = nomes::canonico(&gene.chrom);
        let reason = if !vcf_obs.seen_chroms.contains(&canon) && !callable.is_some_and(|b| b.has_chrom(&canon)) {
            format!("Sequence {} is absent from the VCF and callable BED.", gene.chrom)
        } else if callable.is_none() {
            format!(
                "No {} positions found in VCF and no `callable_bed` was provided (sparse VCFs require callable_bed to distinguish homozygous reference from missing).",
                gene.gene
            )
        } else {
            format!("All {} positions of {} fall outside callable regions or have missing genotypes.", total, gene.gene)
        };
        return Ok(GeneOutcome::uncalled(gene, &method, &obs, evidence, notes, reason));
    }

    let missing_keys: Vec<usize> = gene
        .positions
        .iter()
        .enumerate()
        .filter(|(i, p)| p.is_key && obs[*i].origin == "missing")
        .map(|(i, _)| i)
        .collect();
    // A key position that was not read does not stop the call (PharmCAT drops it from the
    // match, as any missing position); the alleles defined only there are not tested, which
    // the note and the warning say.
    let mut notes = notes;
    let mut key_warning = None;
    if !missing_keys.is_empty() {
        let hint = if callable.is_none() {
            " Pass `callable_bed=` for variant-only VCFs, or supply explicit genotypes."
        } else {
            ""
        };
        let list = format_pos_list(gene, &missing_keys);
        notes.push(format!(
            "{} key position(s) not read for {}: {list}. The call uses only the positions read; alleles defined there could not be tested.",
            missing_keys.len(),
            gene.gene,
        ));
        key_warning = Some(format!(
            "key position(s) not read ({list}); the alleles defined there were not tested, so the call may miss them.{hint}"
        ));
    }

    // Detect haploid vs diploid on sex_x (G6PD) and haploid MT-RNR1:
    let obs_has_haploid = obs.iter().any(|o| o.origin == "vcf" && o.hap_alleles.len() == 1);
    let obs_has_het_diploid = obs
        .iter()
        .any(|o| o.origin == "vcf" && o.hap_alleles.len() == 2 && o.hap_alleles[0] != o.hap_alleles[1]);
    let haploid = match gene.ploidy.as_str() {
        "haploid" => true,
        "sex_x" => obs_has_haploid && !obs_has_het_diploid,
        _ => false,
    };

    let mut outcome = match gene.model.as_str() {
        "variant_activity" | "variant_status" => call_variant_model_gene(gene, &obs, evidence, notes, pheno_table, &method),
        _ => call_combinatorial_gene(gene, &obs, evidence, notes, pheno_table, &method, haploid),
    };
    outcome.missing_key_positions = missing_keys
        .iter()
        .map(|&i| {
            let p = &gene.positions[i];
            p.rsid.clone().unwrap_or_else(|| format!("{}:{}", p.chrom, p.pos))
        })
        .collect();
    if let (true, Some(w)) = (outcome.called, key_warning) {
        outcome.warnings.push(w);
    }
    Ok(outcome)
}

/// Marks the alternatives that would change a recommendation and writes the gene warning
/// when an alternative changes the phenotype or a recommendation. Needs every gene's call,
/// because a recommendation can depend on more than one gene.
fn flag_alternatives(outcomes: &mut [GeneOutcome], recs: &[GuidelineRec]) {
    let light: Vec<GeneOutcome> = outcomes
        .iter()
        .map(|o| GeneOutcome {
            evidence: Vec::new(),
            notes: Vec::new(),
            alternatives: Vec::new(),
            warnings: Vec::new(),
            ..o.clone()
        })
        .collect();
    let signature = |os: &[GeneOutcome], mine: &[&GuidelineRec]| -> BTreeSet<String> {
        match_guidelines(os, mine.iter().copied())
            .into_iter()
            .map(|d| format!("{}\t{}\t{}\t{:?}\t{:?}", d.drug, d.guideline, d.recommendation, d.classification, d.x_population))
            .collect()
    };
    for gi in 0..outcomes.len() {
        if !outcomes[gi].called || outcomes[gi].alternatives.is_empty() {
            continue;
        }
        let gname = outcomes[gi].gene.clone();
        let mine: Vec<&GuidelineRec> = recs
            .iter()
            .filter(|r| r.genes.contains(&gname) || r.untested_genes.contains(&gname) || r.lookup.contains_key(&gname))
            .collect();
        let base = signature(&light, &mine);
        let mut work = light.clone();
        for alt in outcomes[gi].alternatives.iter_mut() {
            work[gi].diplotype = Some(alt.diplotype.clone());
            work[gi].phenotype = alt.phenotype.clone();
            work[gi].activity_score = alt.activity_score;
            work[gi].activity_label = alt.activity_score.map(format_activity_score);
            alt.changes_recommendation = signature(&work, &mine) != base;
        }

        let g = &outcomes[gi];
        let changing: Vec<&Alternative> = g.alternatives.iter().filter(|a| a.changes_phenotype || a.changes_recommendation).collect();
        if changing.is_empty() {
            continue;
        }
        let what = match (
            changing.iter().any(|a| a.changes_phenotype),
            changing.iter().any(|a| a.changes_recommendation),
        ) {
            (true, true) => "the phenotype or a recommendation",
            (true, false) => "the phenotype",
            _ => "a recommendation",
        };
        let mut list: Vec<String> = changing
            .iter()
            .take(5)
            .map(|a| describe_call(&a.diplotype, a.phenotype.as_deref(), a.activity_score, &g.lookup_method))
            .collect();
        if changing.len() > 5 {
            list.push(format!("and {} more", changing.len() - 5));
        }
        let how = match g.chosen_by {
            Some("frequency") => "as the most common diplotype in the CPIC population groups",
            Some("trans") => "assuming the variants are in trans",
            _ => "as the highest-scoring diplotype",
        };
        let called = describe_call(g.diplotype.as_deref().unwrap_or("-"), g.phenotype.as_deref(), g.activity_score, &g.lookup_method);
        let settle = if g.missing_key_positions.is_empty() && g.positions_missing == 0 {
            "Phased genotypes would settle it."
        } else {
            "Phased genotypes, or reading the missing positions, would settle it."
        };
        let verb = if changing.len() == 1 { "is" } else { "are" };
        let warning = format!(
            "{called} was chosen {how}, but {} {verb} also possible with the positions read and would change {what}. {settle}",
            list.join("; ")
        );
        outcomes[gi].warnings.push(warning);
    }
}

// ---------------------------------------------------------------------------
// Guideline Drug Matching & HTML Report Generation
// ---------------------------------------------------------------------------

#[derive(Debug, Serialize)]
struct DrugEntry {
    drug: String,
    #[serde(rename = "x_guideline")]
    guideline: String,
    genes: Vec<String>,
    #[serde(rename = "x_phenotype")]
    phenotype: String,
    recommendation: String,
    #[serde(rename = "strength")]
    classification: Option<String>,
    source: String,
    url: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    x_population: Option<String>,
}

fn match_guidelines<'r>(outcomes: &[GeneOutcome], recs: impl IntoIterator<Item = &'r GuidelineRec>) -> Vec<DrugEntry> {
    let by_gene: HashMap<&str, &GeneOutcome> = outcomes.iter().map(|o| (o.gene.as_str(), o)).collect();
    let mut out = Vec::new();
    for rec in recs {
        if rec.lookup.get("_algorithm").map(String::as_str) == Some("warfarin") {
            let all_called = rec.genes.iter().all(|g| by_gene.get(g.as_str()).is_some_and(|o| o.called));
            if all_called {
                let summary_parts: Vec<String> = rec
                    .genes
                    .iter()
                    .filter_map(|g| {
                        by_gene.get(g.as_str()).map(|o| {
                            format!(
                                "{}: {}",
                                g,
                                o.phenotype.as_deref().or(o.diplotype.as_deref()).unwrap_or("Called")
                            )
                        })
                    })
                    .collect();
                out.push(DrugEntry {
                    drug: rec.drug.clone(),
                    guideline: rec.guideline_name.clone(),
                    genes: rec.genes.clone(),
                    phenotype: summary_parts.join("; "),
                    recommendation: rec.recommendation.clone(),
                    classification: rec.classification.clone(),
                    source: rec.source.clone(),
                    url: rec.url.clone(),
                    x_population: if rec.population == "general" { None } else { Some(rec.population.clone()) },
                });
            }
            continue;
        }

        if rec
            .untested_genes
            .iter()
            .any(|g| by_gene.get(g.as_str()).is_some_and(|o| o.called))
        {
            continue;
        }
        let mut any_called = false;
        let mut all_match = true;
        let mut pheno_parts = Vec::new();
        for (gene_sym, expected) in &rec.lookup {
            let Some(go) = by_gene.get(gene_sym.as_str()) else {
                all_match = false;
                break;
            };
            let matched_val = if !go.called {
                if expected == "No Result" {
                    Some("No Result".to_string())
                } else {
                    None
                }
            } else if go.lookup_method == "ACTIVITY_SCORE" {
                let s = go
                    .activity_label
                    .clone()
                    .or_else(|| go.activity_score.map(format_activity_score))
                    .unwrap_or_else(|| "n/a".to_string());
                if s == *expected {
                    Some(s)
                } else {
                    None
                }
            } else if go.lookup_method == "ALLELE_STATUS" {
                match &go.phenotype {
                    Some(p) if !p.is_empty() => {
                        if p.split(';').map(str::trim).any(|tok| tok == expected.as_str()) {
                            Some(expected.clone())
                        } else {
                            None
                        }
                    }
                    _ => None,
                }
            } else {
                match &go.phenotype {
                    Some(p) if !p.is_empty() && p == expected => Some(p.clone()),
                    _ => None,
                }
            };
            let Some(got) = matched_val else {
                all_match = false;
                break;
            };
            if go.called {
                any_called = true;
                let label = if go.lookup_method == "ACTIVITY_SCORE" {
                    format!(
                        "{}: {} (AS {})",
                        gene_sym,
                        go.phenotype.as_deref().unwrap_or(""),
                        got
                    )
                } else {
                    format!("{gene_sym}: {got}")
                };
                pheno_parts.push(label);
            } else {
                pheno_parts.push(format!("{gene_sym}: No Result"));
            }
        }
        if all_match && any_called && !rec.recommendation.trim().is_empty() {
            out.push(DrugEntry {
                drug: rec.drug.clone(),
                guideline: rec.guideline_name.clone(),
                genes: rec.genes.clone(),
                phenotype: pheno_parts.join("; "),
                recommendation: rec.recommendation.clone(),
                classification: rec.classification.clone(),
                source: rec.source.clone(),
                url: rec.url.clone(),
                x_population: if rec.population == "general" { None } else { Some(rec.population.clone()) },
            });
        }
    }
    out
}

fn html_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

// ---------------------------------------------------------------------------
// Genetic Traits Evaluator (VCF-based, Open Loci)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize)]
struct TraitOutcome {
    #[serde(rename = "trait")]
    trait_name: String,
    genotype: String,
    interpretation: String,
    evidence: String,
    source: String,
    x_id: String,
    x_called: bool,
    x_positions_read: usize,
    x_positions_total: usize,
    x_evidence: Vec<EvidenceEntry>,
    x_notes: Vec<String>,
}

fn hap_allele_str(h: &HapAllele) -> Option<&str> {
    match h {
        HapAllele::Cpic(s) | HapAllele::NonCuratedAlt(s) => Some(s.as_str()),
        HapAllele::SpanningDel | HapAllele::Missing => None,
    }
}

fn diploid_bases(obs: &PosObservation) -> Option<(&str, &str)> {
    if obs.origin == "missing" || obs.hap_alleles.len() != 2 {
        return None;
    }
    Some((hap_allele_str(&obs.hap_alleles[0])?, hap_allele_str(&obs.hap_alleles[1])?))
}

fn format_trait_genotype(positions: &[PositionDef], obs_list: &[PosObservation]) -> String {
    positions
        .iter()
        .zip(obs_list.iter())
        .map(|(p, o)| {
            let label = p.rsid.as_deref().unwrap_or("site");
            format!("{label} {}", o.evidence_gt)
        })
        .collect::<Vec<_>>()
        .join("; ")
}

fn count_base(pair: (&str, &str), target: &str) -> usize {
    usize::from(pair.0 == target) + usize::from(pair.1 == target)
}

fn evaluate_lactase_trait(
    def: &TraitDef,
    obs_list: &[PosObservation],
    evidence: Vec<EvidenceEntry>,
    mut notes: Vec<String>,
) -> TraitOutcome {
    let total = def.positions.len();
    let read = obs_list.iter().filter(|o| o.origin != "missing").count();
    let genotype = format_trait_genotype(&def.positions, obs_list);
    let primary_rs = def.positions.first().and_then(|p| p.rsid.as_deref()).unwrap_or("rs4988235");

    let (x_called, interpretation, ev_summary) = match obs_list.first().and_then(diploid_bases) {
        None => (
            false,
            format!("Not called: primary variant {primary_rs} (MCM6) is missing or outside callable regions"),
            format!("Insufficient coverage ({read}/{total} loci called)"),
        ),
        Some((a, b)) => {
            let a_n = count_base((a, b), "A");
            let g_n = count_base((a, b), "G");
            let interp = if a_n == 2 {
                "Lactase persistence likely (homozygous MCM6 -13910*T, rs4988235 A/A; dominant trait)".to_string()
            } else if a_n == 1 {
                "Lactase persistence likely (heterozygous MCM6 -13910*T carrier, rs4988235; dominant trait)".to_string()
            } else if g_n == 2 {
                "Lactase non-persistence predicted at European MCM6 -13910C>T (rs4988235 G/G); non-European LCT persistence variants not excluded".to_string()
            } else {
                format!(
                    "Uncertain: non-standard allele at MCM6 {primary_rs} ({a}/{b}); not classified by European -13910*T rule"
                )
            };
            if let Some((t1, t2)) = obs_list.get(1).and_then(diploid_bases) {
                let tag_t = count_base((t1, t2), "T");
                if (a_n > 0) != (tag_t > 0) {
                    notes.push(format!(
                        "Linkage note: primary rs4988235 ({a}/{b}) and tag rs182549 ({t1}/{t2}) alleles are discordant."
                    ));
                }
            }
            (a_n > 0 || g_n == 2, interp, format!("{} ({read}/{total} loci called)", def.evidence_level))
        }
    };

    TraitOutcome {
        trait_name: def.name.clone(),
        genotype,
        interpretation,
        evidence: ev_summary,
        source: def.source.clone(),
        x_id: def.id.clone(),
        x_called,
        x_positions_read: read,
        x_positions_total: total,
        x_evidence: evidence,
        x_notes: notes,
    }
}

fn evaluate_eye_color_trait(
    def: &TraitDef,
    obs_list: &[PosObservation],
    evidence: Vec<EvidenceEntry>,
    notes: Vec<String>,
) -> TraitOutcome {
    const MODIFIERS: &[(&str, &str)] = &[
        ("rs1800407", "T"),
        ("rs12896399", "T"),
        ("rs16891982", "G"),
        ("rs1393350", "A"),
        ("rs12203592", "T"),
    ];
    let total = def.positions.len();
    let read = obs_list.iter().filter(|o| o.origin != "missing").count();
    let genotype = format_trait_genotype(&def.positions, obs_list);

    let (x_called, interpretation, ev_summary) = match obs_list.first().and_then(diploid_bases) {
        None => (
            false,
            "Not called: primary locus HERC2 rs12913832 is missing or outside callable regions".to_string(),
            format!("Insufficient coverage ({read}/{total} SNPs called)"),
        ),
        Some((h1, h2)) => {
            let g_n = count_base((h1, h2), "G");
            let a_n = count_base((h1, h2), "A");
            let mut mod_called = 0usize;
            let mut light_mods = 0usize;
            let mut oca2_t = 0usize;
            for (idx, (rs, light_base)) in MODIFIERS.iter().enumerate() {
                let obs_opt = def
                    .positions
                    .iter()
                    .zip(obs_list.iter())
                    .find(|(p, _)| p.rsid.as_deref() == Some(*rs))
                    .map(|(_, o)| o)
                    .or_else(|| obs_list.get(idx + 1));
                if let Some(pair) = obs_opt.and_then(diploid_bases) {
                    mod_called += 1;
                    let cnt = count_base(pair, light_base);
                    light_mods += cnt;
                    if *rs == "rs1800407" {
                        oca2_t = cnt;
                    }
                }
            }
            let missing_mods = MODIFIERS.len().saturating_sub(mod_called);
            let cov_note = if missing_mods > 0 {
                format!("; {read}/{total} SNPs called ({missing_mods} modifier SNP(s) missing/outside callable BED)")
            } else {
                format!("; {total}/{total} SNPs called")
            };
            if g_n == 2 {
                let interp = if light_mods >= 3 || oca2_t >= 1 {
                    format!(
                        "Predicted iris color: blue / light (HERC2 rs12913832 G/G with {light_mods} light-modifier allele(s){cov_note})"
                    )
                } else {
                    format!(
                        "Predicted iris color: blue / light most likely, green or hazel possible (HERC2 rs12913832 G/G with {light_mods} light-modifier allele(s){cov_note})"
                    )
                };
                (true, interp, format!("{} ({read}/{total} SNPs called)", def.evidence_level))
            } else if a_n == 2 {
                (
                    true,
                    format!(
                        "Predicted iris color: brown / dark most likely (HERC2 rs12913832 A/A with {light_mods} light-modifier allele(s){cov_note})"
                    ),
                    format!("{} ({read}/{total} SNPs called)", def.evidence_level),
                )
            } else if g_n == 1 && a_n == 1 {
                let interp = if oca2_t >= 1 && light_mods >= 3 {
                    format!(
                        "Predicted iris color: intermediate (green / hazel) or light brown (HERC2 rs12913832 A/G with OCA2 rs1800407*T and {light_mods} light-modifier allele(s){cov_note})"
                    )
                } else if light_mods >= 2 {
                    format!(
                        "Predicted iris color: brown or intermediate (hazel / green) (HERC2 rs12913832 A/G with {light_mods} light-modifier allele(s){cov_note})"
                    )
                } else {
                    format!(
                        "Predicted iris color: brown most likely, intermediate possible (HERC2 rs12913832 A/G with {light_mods} light-modifier allele(s){cov_note})"
                    )
                };
                (true, interp, format!("{} ({read}/{total} SNPs called)", def.evidence_level))
            } else {
                (
                    false,
                    format!("Uncertain: non-standard alleles at HERC2 rs12913832 ({h1}/{h2})"),
                    format!("Insufficient or non-standard alleles ({read}/{total} SNPs called)"),
                )
            }
        }
    };

    TraitOutcome {
        trait_name: def.name.clone(),
        genotype,
        interpretation,
        evidence: ev_summary,
        source: def.source.clone(),
        x_id: def.id.clone(),
        x_called,
        x_positions_read: read,
        x_positions_total: total,
        x_evidence: evidence,
        x_notes: notes,
    }
}

fn evaluate_traits(
    traits: &[TraitDef],
    vcf_obs: &VcfObservations,
    callable: Option<&CallableBed>,
    genome: &mut Option<referencia::Genoma>,
) -> Result<Vec<TraitOutcome>, String> {
    let mut out = Vec::with_capacity(traits.len());
    for def in traits {
        let (obs_list, evidence, mut notes) = observe_position_list(&def.positions, "diploid", vcf_obs, callable, genome)?;
        for (p, o) in def.positions.iter().zip(obs_list.iter()) {
            let rs = p.rsid.as_deref().unwrap_or("site");
            if o.origin == "vcf" && o.in_callable == Some(false) {
                notes.push(format!("{rs} observed in VCF outside callable BED regions."));
            } else if o.origin == "missing" {
                notes.push(format!("{rs} missing or outside callable BED regions."));
            }
        }
        let outcome = match def.model.as_str() {
            "lactase_mcm6" => evaluate_lactase_trait(def, &obs_list, evidence, notes),
            "eye_color_6snp" => evaluate_eye_color_trait(def, &obs_list, evidence, notes),
            other => {
                return Err(format!("Unsupported trait model {other:?} in traits.json for trait {:?}.", def.id));
            }
        };
        out.push(outcome);
    }
    Ok(out)
}

fn render_html_report(
    sample: &str,
    data_version: &str,
    outcomes: &[GeneOutcome],
    drugs: &[DrugEntry],
    traits: &[TraitOutcome],
    bam_note: Option<&str>,
    warnings: &[String],
) -> String {
    let called_n = outcomes.iter().filter(|g| g.called).count();
    let not_called_n = outcomes.len() - called_n;
    let mut gene_rows = String::new();
    for g in outcomes {
        let status = if g.called { "<span class=\"ok\">Called</span>" } else { "<span class=\"warn\">Not Called</span>" };
        let dip = g.diplotype.as_deref().unwrap_or("-");
        let pheno = g.phenotype.as_deref().unwrap_or("-");
        let as_str = g
            .activity_label
            .clone()
            .or_else(|| g.activity_score.map(format_activity_score))
            .unwrap_or_else(|| "-".to_string());
        let cn_str = g
            .copy_number
            .map(|c| {
                if (c - c.round()).abs() < 0.05 {
                    format!("{:.0}", c.round())
                } else {
                    format!("{c:.2}")
                }
            })
            .unwrap_or_else(|| "-".to_string());
        let joined_notes = g.notes.join(" ");
        let note_str = if g.called { joined_notes.as_str() } else { g.not_called_reason.as_deref().unwrap_or("") };
        gene_rows.push_str(&format!(
            "<tr><td><strong>{}</strong></td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}/{}</td><td>{}</td></tr>\n",
            html_escape(&g.gene),
            status,
            html_escape(dip),
            html_escape(pheno),
            html_escape(&as_str),
            html_escape(&cn_str),
            html_escape(&g.source),
            g.positions_read,
            g.positions_total,
            html_escape(note_str)
        ));
    }

    let mut drug_rows = String::new();
    for d in drugs {
        let pop = d.x_population.as_deref().unwrap_or("general");
        drug_rows.push_str(&format!(
            "<tr><td><strong>{}</strong></td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>\n",
            html_escape(&d.drug),
            html_escape(&d.genes.join(", ")),
            html_escape(&d.phenotype),
            html_escape(pop),
            html_escape(d.classification.as_deref().unwrap_or("-")),
            html_escape(&d.recommendation),
        ));
    }

    let mut trait_rows = String::new();
    for tr in traits {
        trait_rows.push_str(&format!(
            "<tr><td><strong>{}</strong></td><td><code>{}</code></td><td>{}</td><td>{}</td><td>{}</td></tr>\n",
            html_escape(&tr.trait_name),
            html_escape(&tr.genotype),
            html_escape(&tr.interpretation),
            html_escape(&tr.evidence),
            html_escape(&tr.source),
        ));
    }

    let bam_banner = bam_note
        .map(|n| format!("<p class=\"note\">{}</p>", html_escape(n)))
        .unwrap_or_default();
    let warnings_html = if warnings.is_empty() {
        String::new()
    } else {
        let items: String = warnings.iter().map(|w| format!("<li>{}</li>", html_escape(w))).collect();
        format!("<h2>Warnings</h2>\n<ul class=\"warnlist\">{items}</ul>\n")
    };

    format!(
        "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n\
         <title>MMGEN Pharmacogenetics Report - {sample}</title>\n\
         <style>body{{font-family:system-ui,-apple-system,sans-serif;max-width:1180px;margin:2rem auto;padding:0 1rem;color:#1f2937}}\
         h1,h2{{color:#111827}} .meta,.note{{background:#f3f4f6;padding:.75rem 1rem;border-radius:6px;margin-bottom:1rem}}\
         table{{border-collapse:collapse;width:100%;margin-bottom:2rem;font-size:.9rem}} th,td{{border:1px solid #d1d5db;padding:.5rem .6rem;text-align:left;vertical-align:top}}\
         th{{background:#f9fafb}} .ok{{color:#065f46;font-weight:600}} .warn{{color:#92400e;font-weight:600}}</style>\n</head>\n<body>\n\
         <h1>MMGEN Pharmacogenetics Report (<code>mmgen-pgx</code>)</h1>\n\
         <div class=\"meta\"><strong>Sample:</strong> {sample} | <strong>Assembly:</strong> GRCh38 | <strong>Data:</strong> {data_version} | \
         <strong>Genes Called:</strong> {called_n}/{total} ({not_called_n} not called) | <strong>Recommendations:</strong> {n_drugs} | <strong>Traits:</strong> {n_traits}</div>\n\
         {bam_banner}\n\
         <p class=\"note\">Experimental software for teaching and research only: not a medical diagnosis or treatment guidance.</p>\n\
         {warnings_html}\
         <h2>Gene Calls</h2>\n\
         <table><thead><tr><th>Gene</th><th>Status</th><th>Diplotype</th><th>Phenotype</th><th>Activity Score</th><th>Copy Number</th><th>Source</th><th>Positions</th><th>Reason / Notes</th></tr></thead>\n\
         <tbody>\n{gene_rows}</tbody></table>\n\
         <h2>CPIC Guideline Recommendations</h2>\n\
         <table><thead><tr><th>Drug</th><th>Genes</th><th>Phenotype / Lookup</th><th>Population</th><th>Classification</th><th>Recommendation</th></tr></thead>\n\
         <tbody>\n{drug_rows}</tbody></table>\n\
         <h2>Genetic Traits</h2>\n\
         <table><thead><tr><th>Trait</th><th>Genotype</th><th>Interpretation</th><th>Evidence</th><th>Source</th></tr></thead>\n\
         <tbody>\n{trait_rows}</tbody></table>\n\
         </body>\n</html>\n",
        sample = html_escape(sample),
        data_version = html_escape(data_version),
        total = outcomes.len(),
        n_drugs = drugs.len(),
        n_traits = traits.len(),
    )
}

// ---------------------------------------------------------------------------
// Public Entry Point
// ---------------------------------------------------------------------------

/// Execute the native `mmgen-pgx` v1 pharmacogenetic pipeline.
pub fn run(vcf: &str, output_dir: &str, opts: &PgxOptions) -> Result<PgxRunResult, String> {
    if let Some(bam_path) = &opts.bam {
        if !Path::new(bam_path).exists() {
            return Err(format!("{bam_path}: BAM file does not exist."));
        }
    }
    let bundle = load_bundle(opts.tables.as_deref())?;
    std::fs::create_dir_all(output_dir).map_err(|e| format!("{output_dir}: cannot create output directory ({e})."))?;

    let mut target_chroms = HashSet::new();
    for g in &bundle.genes.genes {
        target_chroms.insert(nomes::canonico(&g.chrom));
    }
    for t in &bundle.traits.traits {
        for p in &t.positions {
            target_chroms.insert(nomes::canonico(&p.chrom));
        }
    }
    let callable = match &opts.callable_bed {
        Some(p) => Some(CallableBed::load(p, &target_chroms)?),
        None => None,
    };
    let mut genome = match &opts.reference {
        Some(r) => Some(referencia::Genoma::abre(r, Some(Path::new(output_dir)))?),
        None => None,
    };

    let windows = build_query_windows(&bundle.genes.genes, &bundle.traits.traits);
    let vcf_obs = read_vcf_observations(vcf, opts.sample.as_deref(), &windows)?;
    let trait_outcomes = evaluate_traits(&bundle.traits.traits, &vcf_obs, callable.as_ref(), &mut genome)?;

    let mut outcomes = Vec::with_capacity(bundle.genes.genes.len());
    for g in &bundle.genes.genes {
        outcomes.push(evaluate_gene(
            g,
            &vcf_obs,
            callable.as_ref(),
            bundle.phenotypes.genes.get(&g.gene),
            &mut genome,
        )?);
    }
    outcomes.sort_by(|a, b| a.gene.cmp(&b.gene));

    let bam_note = if let Some(bam_path) = opts.bam.as_deref() {
        let bam_outcomes = pgx_bam::call_bam_genes(&bundle, bam_path, &vcf_obs, callable.as_ref(), &mut genome)?;
        outcomes.extend(bam_outcomes);
        None
    } else if bundle.bam_genes.is_some() {
        Some(
            "No BAM file was provided; BAM-dependent genes (CYP2D6, HLA-A, HLA-B, HLA-C) were skipped (pass bam=<sorted_grch38.bam> to enable them)."
                .to_string(),
        )
    } else {
        None
    };

    flag_alternatives(&mut outcomes, &bundle.guidelines.recommendations);
    let drugs = match_guidelines(&outcomes, &bundle.guidelines.recommendations);

    let result_json_path = Path::new(output_dir).join("result.json").to_string_lossy().to_string();
    let report_html_path = Path::new(output_dir).join("report.html").to_string_lossy().to_string();

    let called_n = outcomes.iter().filter(|g| g.called).count();
    let not_called_n = outcomes.len() - called_n;
    let unique_drugs: BTreeSet<&str> = drugs.iter().map(|d| d.drug.as_str()).collect();

    let gene_alleles = |g: &GeneOutcome| -> Vec<String> {
        if !g.called {
            return Vec::new();
        }
        match (&g.haplotype_1, &g.haplotype_2) {
            (Some(h1), Some(h2)) => vec![h1.clone(), h2.clone()],
            (Some(h1), None) => vec![h1.clone()],
            _ => match &g.diplotype {
                Some(d) if !d.contains(" OR ") && d.contains(" AND ") => {
                    d.split(" AND ").map(|s| s.trim().to_string()).collect()
                }
                _ => Vec::new(),
            },
        }
    };

    let genes_json: Vec<serde_json::Value> = outcomes
        .iter()
        .map(|g| {
            let mut obj = serde_json::Map::new();
            obj.insert("gene".into(), json!(g.gene));
            obj.insert("called".into(), json!(g.called));
            obj.insert("diplotype".into(), json!(g.diplotype));
            obj.insert("alleles".into(), json!(gene_alleles(g)));
            obj.insert("phenotype".into(), json!(g.phenotype));
            obj.insert("activity_score".into(), json!(g.activity_score));
            obj.insert("source".into(), json!(g.source));
            if let Some(cn) = g.copy_number {
                let cn_val = if (cn - cn.round()).abs() < 0.05 {
                    json!(cn.round() as u64)
                } else {
                    json!((cn * 100.0).round() / 100.0)
                };
                obj.insert("copy_number".into(), cn_val);
            } else {
                obj.insert("copy_number".into(), serde_json::Value::Null);
            }
            obj.insert("positions_read".into(), json!(g.positions_read));
            obj.insert("positions_missing".into(), json!(g.positions_missing));
            obj.insert("evidence".into(), serde_json::to_value(&g.evidence).unwrap_or(json!([])));
            obj.insert("not_called_reason".into(), json!(g.not_called_reason));
            obj.insert("notes".into(), json!(g.notes));
            obj.insert("x_calling_method".into(), json!(g.calling_method));
            obj.insert("x_phased".into(), json!(g.phased));
            obj.insert("x_positions_total".into(), json!(g.positions_total));
            if g.haplotype_1.is_some() || g.haplotype_2.is_some() {
                obj.insert("x_haplotype_1".into(), json!(g.haplotype_1));
                obj.insert("x_haplotype_2".into(), json!(g.haplotype_2));
            }
            if let Some(s) = g.score {
                obj.insert("x_score".into(), json!(s));
            }
            if let Some(c) = g.chosen_by {
                obj.insert("x_chosen_by".into(), json!(c));
            }
            if let Some(fq) = g.frequency {
                obj.insert("x_frequency".into(), json!(fq));
            }
            if !g.alternatives.is_empty() {
                obj.insert("x_alternatives".into(), serde_json::to_value(&g.alternatives).unwrap_or(json!([])));
            }
            if !g.missing_key_positions.is_empty() {
                obj.insert("x_missing_key_positions".into(), json!(g.missing_key_positions));
            }
            serde_json::Value::Object(obj)
        })
        .collect();

    let mut warnings: Vec<String> = Vec::new();
    if opts.callable_bed.is_none() {
        warnings.push(
            "No callable_bed was provided; absent positions were treated as missing rather than homozygous reference.".to_string(),
        );
    }
    if let Some(ref bn) = bam_note {
        warnings.push(bn.clone());
    }
    for g in &outcomes {
        for n in &g.notes {
            if n.contains("disagrees with GRCh38") {
                warnings.push(format!("{}: {n}", g.gene));
            }
        }
        for w in &g.warnings {
            warnings.push(format!("{}: {w}", g.gene));
        }
    }

    // The HLA calls depend on the IPD-IMGT/HLA release too: name it next to the CPIC version.
    let mut data_version = bundle.manifest.data_version.clone();
    if opts.bam.is_some() {
        let hla_version = bundle.bam_genes.as_ref().and_then(|b| b.hla_genes.iter().find_map(|h| h.imgt_hla_version.clone()));
        if let Some(v) = hla_version {
            data_version = format!("{data_version}; IPD-IMGT/HLA {v}");
        }
    }

    let text_summary = format!(
        "{} v{} ({}) analyzed {} pharmacogenes on GRCh38 for sample {}: {}/{} genes called ({} not called); {} CPIC recommendation(s) across {} drug(s).",
        ENGINE_NAME,
        ENGINE_VERSION,
        data_version,
        outcomes.len(),
        vcf_obs.sample,
        called_n,
        outcomes.len(),
        not_called_n,
        drugs.len(),
        unique_drugs.len(),
    );

    let mut root = serde_json::Map::new();
    root.insert("schema_version".into(), json!(SCHEMA_VERSION));
    root.insert("engine".into(), json!(ENGINE_NAME));
    root.insert("engine_version".into(), json!(ENGINE_VERSION));
    root.insert("data_version".into(), json!(data_version));
    root.insert("sample".into(), json!(vcf_obs.sample));
    root.insert("genes".into(), serde_json::Value::Array(genes_json));
    root.insert("drugs".into(), serde_json::to_value(&drugs).unwrap_or(json!([])));
    root.insert("traits".into(), serde_json::to_value(&trait_outcomes).unwrap_or(json!([])));
    root.insert("warnings".into(), json!(warnings));
    let mut x_files = serde_json::Map::new();
    x_files.insert("report_html".into(), json!(report_html_path));
    x_files.insert("result_json".into(), json!(result_json_path));
    if opts.bam.is_some() {
        let mut outside_lines = Vec::new();
        for g in &outcomes {
            if matches!(g.gene.as_str(), "CYP2D6" | "HLA-A" | "HLA-B") && g.called {
                if let Some(ref dip) = g.diplotype {
                    outside_lines.push(format!("{}\t{}", g.gene, dip));
                }
            }
        }
        let outside_path = Path::new(&result_json_path)
            .with_file_name("outside_calls.tsv")
            .to_string_lossy()
            .to_string();
        let outside_text = if outside_lines.is_empty() {
            String::new()
        } else {
            outside_lines.join("\n") + "\n"
        };
        std::fs::write(&outside_path, outside_text.as_bytes())
            .map_err(|e| format!("{outside_path}: error writing outside_calls.tsv ({e})."))?;
        x_files.insert("pharmcat_outside_calls".into(), json!(outside_path));
    }

    root.insert("report_html".into(), json!(report_html_path));
    root.insert("text".into(), json!(text_summary));
    root.insert("x_assembly".into(), json!("GRCh38"));
    root.insert("x_files".into(), serde_json::Value::Object(x_files));
    root.insert(
        "x_summary".into(),
        json!({
            "genes_total": outcomes.len(),
            "genes_called": called_n,
            "genes_not_called": not_called_n,
            "drugs_with_recommendation": unique_drugs.len(),
        }),
    );
    if let Some(ref bn) = bam_note {
        root.insert("x_bam_note".into(), json!(bn));
    }

    let json_text = serde_json::to_string_pretty(&serde_json::Value::Object(root))
        .map_err(|e| format!("Failed to serialize PGx result JSON: {e}"))?
        + "\n";
    let html_text = render_html_report(
        &vcf_obs.sample,
        &data_version,
        &outcomes,
        &drugs,
        &trait_outcomes,
        bam_note.as_deref(),
        &warnings,
    );

    std::fs::write(&result_json_path, json_text.as_bytes())
        .map_err(|e| format!("{result_json_path}: error writing result.json ({e})."))?;
    std::fs::write(&report_html_path, html_text.as_bytes())
        .map_err(|e| format!("{report_html_path}: error writing report.html ({e})."))?;

    Ok(PgxRunResult { result_json_path, report_html_path, json_text })
}
