#![allow(dead_code)]
// Tests of src/pharmcat.rs: the PharmCAT input prepared from a variants-only
// VCF, checked position by position against hand-worked expectations, through
// a plain file and through a BGZF file with its .tbi (the two must agree).
// With PHARMCAT_JAR set (and Java installed), PharmCAT itself is run on it.
// Usage: cargo run --release --bin pharmcat
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
#[path = "../../../src/nomes.rs"]
mod nomes;
#[path = "../../../src/indice.rs"]
mod indice;
#[path = "../../../src/pharmcat.rs"]
mod pharmcat;

use std::io::Read;

static mut DONE: usize = 0;
static mut FAILED: usize = 0;

fn check(name: &str, cond: bool, detail: &str) {
    unsafe {
        DONE += 1;
        if !cond {
            FAILED += 1;
        }
    }
    if cond {
        println!("  ok    {name}");
    } else {
        println!("  FAIL  {name}  {detail}");
    }
}

fn write(path: &str, text: &str) {
    std::fs::write(path, text).unwrap();
}

fn read_text(path: &str) -> String {
    let mut s = String::new();
    externo::abre_leitura(path).unwrap().read_to_string(&mut s).unwrap();
    s
}

/// The GT written for each position, by "chrom:pos".
fn gts(path: &str) -> Vec<(String, String, String)> {
    read_text(path)
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (format!("{}:{}", f[0], f[1]), f[4].to_string(), f[9].to_string())
        })
        .collect()
}

fn main() {
    let dir = std::env::temp_dir().join(format!("mmg_pharmcat_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = |n: &str| dir.join(n).to_string_lossy().to_string();

    // The positions file: PharmCAT's layout (chr names, multiallelic lists).
    let positions = "##fileformat=VCFv4.2\n##source=test positions\n##INFO=<ID=PX,Number=.,Type=String,Description=\"Gene\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n\
chr1\t100\trsA\tG\tA\t.\tPASS\tPX=G1\tGT\t0/0\n\
chr1\t200\trsB\tT\tC,G\t.\tPASS\tPX=G1\tGT\t0/0\n\
chr1\t300\trsC\tCAT\tC,CATAT\t.\tPASS\tPX=G1\tGT\t0/0\n\
chr1\t400\trsD\tA\tG\t.\tPASS\tPX=G1\tGT\t0/0\n\
chr1\t500\trsE\tC\tT\t.\tPASS\tPX=G1\tGT\t0/0\n\
chr1\t600\trsF\tG\tT\t.\tPASS\tPX=G1\tGT\t0/0\n\
chr1\t650\trsM\tC\tT\t.\tPASS\tPX=G1\tGT\t0/0\n\
chr1\t700\trsG\tTA\tT\t.\tPASS\tPX=G1\tGT\t0/0\n\
chr1\t800\trsH\tA\tC\t.\tPASS\tPX=G1\tGT\t0/0\n\
chr1\t900\trsI\tG\tA\t.\tPASS\tPX=G1\tGT\t0/0\n\
chr2\t100\trsJ\tC\tT\t.\tPASS\tPX=G2\tGT\t0/0\n\
chr2\t5000\trsK\tA\tG\t.\tPASS\tPX=G2\tGT\t0/0\n\
chrX\t100\trsL\tG\tA\t.\tPASS\tPX=G3\tGT\t0/0\n";
    write(&p("positions.vcf"), positions);

    // The person: `1`/`2` names (no chr), three samples; no chrX at all.
    let person = "##fileformat=VCFv4.2\n##contig=<ID=1,length=10000>\n##contig=<ID=2,length=10000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tother\tme\tthird\n\
1\t100\t.\tG\tA\t50\tPASS\t.\tGT:DP\t0/0:10\t0/1:30\t1/1:20\n\
1\t200\t.\tT\tG\t50\tPASS\t.\tGT\t0/1\t1|0\t0/0\n\
1\t300\t.\tCATAT\tCAT\t50\tPASS\t.\tGT\t0/1\t1/1\t0/0\n\
1\t400\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\n\
1\t500\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\t./.\t0/0\n\
1\t650\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\n\
1\t700\t.\tT\tTC\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\n\
1\t900\t.\tG\tA\t50\tPASS\t.\tGT\t0/1\t1\t0/0\n\
2\t100\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\t0/0\t0/0\n";
    write(&p("person.vcf"), person);
    // Callable: chr1 1-850 and chr2 1-1000 (so chr1:900 has a record anyway,
    // chr1:800 is absent inside, chr2:5000 absent outside).
    write(&p("callable.bed"), "chr1\t0\t850\nchr2\t0\t1000\n");

    let op = pharmcat::Opcoes { amostra: Some("me".into()), chamavel: Some(p("callable.bed")), ausente_e_ref: true };

    println!("== positions file");
    let pos = pharmcat::le_posicoes(&p("positions.vcf")).unwrap();
    check("13 positions, multiallelic lists kept", pos.lista.len() == 13 && pos.lista[1].alts == vec!["C", "G"], "");

    println!("== a plain VCF, sample `me`, with a callable BED");
    let r = pharmcat::preparar(&p("person.vcf"), &p("positions.vcf"), &p("out_plain.vcf"), &op).unwrap();
    let got = gts(&p("out_plain.vcf"));
    let expect = [
        ("chr1:100", "A", "0/1"),         // plain het
        ("chr1:200", "C,G", "2|0"),       // second ALT of the list, phased kept
        ("chr1:300", "C,CATAT", "1/1"),   // longer REF reconciled: CATAT>CAT == CAT>C
        ("chr1:400", "G,T", "0/2"),       // allele not in the list: added as the 2nd ALT
        ("chr1:500", "T", "./."),         // the person's record has no genotype
        ("chr1:600", "T", "0/0"),         // absent, inside the callable region
        ("chr1:650", "T", "./."),         // the person's REF (A) is not the position's (C)
        ("chr1:700", "T,TCA", "0/2"),     // T>TC written over REF TA is TA>TCA: added
        ("chr1:800", "C", "0/0"),         // absent, inside
        ("chr1:900", "A", "1"),           // haploid call kept
        ("chr2:100", "T", "0/0"),         // a reference record
        ("chr2:5000", "G", "./."),        // absent, outside the callable region
        ("chrX:100", "A", "./."),         // a sequence the file does not have
    ];
    for (k, alt, gt) in expect {
        let row = got.iter().find(|(kk, _, _)| kk == k);
        check(&format!("{k}: ALT {alt}, GT {gt}"), row.map_or(false, |(_, a, g)| a == alt && g == gt), &format!("{row:?}"));
    }
    check("the counts", r.com_variante == 6 && r.registro_referencia == 1 && r.referencia_preenchida == 2
        && r.faltando_fora_da_regiao == 1 && r.faltando_sem_sequencia == 1 && r.faltando_no_registro == 1
        && r.ref_divergente == 1 && r.alelos_acrescentados == 2 && r.amostra == "me" && !r.indexado
        && r.posicoes == 13, &r.json());
    check("chrX reported as absent", r.sequencias_ausentes == vec!["chrX".to_string()], &format!("{:?}", r.sequencias_ausentes));
    check("the JSON parses as an object", r.json().starts_with('{') && r.json().ends_with('}') && r.json().contains("\"ref_mismatch\":1"), &r.json());
    let header = read_text(&p("out_plain.vcf"));
    check("the header: positions file lines, the MMGEN line and the sample column `me`",
          header.contains("##source=test positions") && header.contains("##MMGEN=") && header.contains("FORMAT\tme\n"), "");

    println!("== the same through BGZF + .tbi");
    let mut w = externo::Escritor::cria(&p("person.vcf.gz")).unwrap();
    for l in person.lines() {
        w.linha(l).unwrap();
    }
    w.termina().unwrap();
    indice::indexar(&p("person.vcf.gz"), None).unwrap();
    let r2 = pharmcat::preparar(&p("person.vcf.gz"), &p("positions.vcf"), &p("out_indexed.vcf"), &op).unwrap();
    check("indexed: the same genotypes as the plain file", gts(&p("out_indexed.vcf")) == got, "");
    check("indexed: reported as indexed, same counts", r2.indexado && r2.lidas() == r.lidas() && r2.faltando_sem_sequencia == 1, &r2.json());

    println!("== without a BED");
    let op2 = pharmcat::Opcoes { amostra: Some("me".into()), chamavel: None, ausente_e_ref: true };
    pharmcat::preparar(&p("person.vcf"), &p("positions.vcf"), &p("out_nobed.vcf"), &op2).unwrap();
    let g = gts(&p("out_nobed.vcf"));
    let gt_of = |v: &Vec<(String, String, String)>, k: &str| v.iter().find(|(kk, _, _)| kk == k).map(|x| x.2.clone()).unwrap_or_default();
    check("absent_is_ref: chr2:5000 becomes 0/0", gt_of(&g, "chr2:5000") == "0/0", "");
    check("but a sequence absent from the file stays missing", gt_of(&g, "chrX:100") == "./.", "");
    let op3 = pharmcat::Opcoes { amostra: Some("me".into()), chamavel: None, ausente_e_ref: false };
    pharmcat::preparar(&p("person.vcf"), &p("positions.vcf"), &p("out_noref.vcf"), &op3).unwrap();
    let g = gts(&p("out_noref.vcf"));
    check("absent_is_ref=false: absent positions are missing", gt_of(&g, "chr1:600") == "./." && gt_of(&g, "chr2:5000") == "./.", "");

    println!("== samples and errors");
    let op4 = pharmcat::Opcoes { amostra: None, chamavel: Some(p("callable.bed")), ausente_e_ref: true };
    let r4 = pharmcat::preparar(&p("person.vcf"), &p("positions.vcf"), &p("out_first.vcf"), &op4).unwrap();
    check("no sample given: the first one (`other`)", r4.amostra == "other" && gt_of(&gts(&p("out_first.vcf")), "chr1:100") == "0/0", "");
    let e = pharmcat::preparar(&p("person.vcf"), &p("positions.vcf"), &p("x.vcf"),
        &pharmcat::Opcoes { amostra: Some("nobody".into()), chamavel: None, ausente_e_ref: true }).err().unwrap_or_default();
    check("a wrong sample: the error lists the samples", e.contains("nobody") && e.contains("other, me, third"), &e);
    write(&p("sites.vcf"), "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tG\tA\t50\tPASS\t.\n");
    let e = pharmcat::preparar(&p("sites.vcf"), &p("positions.vcf"), &p("x.vcf"), &op2).err().unwrap_or_default();
    check("a sites-only VCF is refused", e.contains("no sample column"), &e);
    let e = pharmcat::preparar(&p("person.vcf"), &p("person.vcf"), &p("x.vcf"), &op2).err().unwrap_or_default();
    check("positions from a file that is not one: still read (it is a VCF)", e.is_empty(), &e);
    write(&p("bad_positions.txt"), "chr1 100\n");
    let e = pharmcat::preparar(&p("person.vcf"), &p("bad_positions.txt"), &p("x.vcf"), &op2).err().unwrap_or_default();
    check("a positions file that is not a VCF: refused with the fix", e.contains("pharmcat_positions"), &e);
    let e = pharmcat::preparar(&p("missing.vcf"), &p("positions.vcf"), &p("x.vcf"), &op2).err().unwrap_or_default();
    check("a missing input file", e.contains("does not exist"), &e);

    println!("== the calls table and the positions next to the jar");
    write(&p("x.report.tsv"), "PharmCAT 3.4.0\nGene\tSource Diplotype\tPhenotype\tActivity Score\n\
CYP2C19\t*1/*2\tIntermediate Metabolizer\t\nDPYD\tReference/Reference\tNormal Metabolizer\t2.0\n");
    let (v, c) = pharmcat::le_chamadas(&p("x.report.tsv")).unwrap();
    check("version and two genes, snake_case columns",
          v == "PharmCAT 3.4.0" && c.len() == 2 && c[0].campos[1] == ("source_diplotype".to_string(), "*1/*2".to_string())
          && c[1].campos[3] == ("activity_score".to_string(), "2.0".to_string()), &format!("{v}"));
    std::fs::create_dir_all(dir.join("pc")).unwrap();
    for n in ["pharmcat-3.4.0-all.jar", "pharmcat_positions_3.3.0.vcf", "pharmcat_positions_3.4.0.vcf", "pharmcat_positions_3.4.0.uniallelic.vcf"] {
        write(&dir.join("pc").join(n).to_string_lossy(), "");
    }
    let found = pharmcat::posicoes_junto_do_jar(&dir.join("pc").join("pharmcat-3.4.0-all.jar").to_string_lossy()).unwrap_or_default();
    check("the positions of the same version are found next to the jar", found.ends_with("pharmcat_positions_3.4.0.vcf"), &found);

    println!("== PharmCAT's report in the common schema (pgx_result.schema.json, version 1)");
    write(&p("x.report.json"), r#"{"pharmcatVersion":"3.4.0","dataVersion":"2026-07-13",
 "genes":{
  "CYP2C19":{"chr":"chr10","callSource":"MATCHER","messages":[],
   "sourceDiplotypes":[{"allele1":{"name":"*1"},"allele2":{"name":"*2"},"label":"*1/*2","phenotypes":["Intermediate Metabolizer"],"activityScore":null}],
   "variants":[{"chromosome":"chr10","position":94781859,"dbSnpId":"rs4244285","call":"G/A"},
               {"chromosome":"chr10","position":94761900,"dbSnpId":"rs12248560","call":null}]},
  "CYP2C9":{"chr":"chr10","callSource":"MATCHER","messages":[],
   "sourceDiplotypes":[{"allele1":{"name":"*1"},"allele2":{"name":"*1"},"phenotypes":["Normal Metabolizer"],"activityScore":"2.0"}],
   "variants":[{"chromosome":"chr10","position":94942290,"dbSnpId":"rs1799853","call":"C/C"}]},
  "SLCO1B1":{"chr":"chr12","callSource":"MATCHER","messages":[],
   "sourceDiplotypes":[{"allele1":{"name":"*1"},"allele2":{"name":"*15"},"phenotypes":["Decreased Function"]},
                       {"allele1":{"name":"*1"},"allele2":{"name":"*45"},"phenotypes":["Decreased Function"]}],
   "variants":[{"chromosome":"chr12","position":21178615,"dbSnpId":"rs2306283","call":"A/G"}]},
  "G6PD":{"chr":"chrX","callSource":"MATCHER","messages":[],
   "sourceDiplotypes":[{"allele1":{"name":"Unknown"},"allele2":null,"phenotypes":["No Result"]}],
   "variants":[{"chromosome":"chrX","position":154532046,"dbSnpId":null,"call":null}]},
  "CYP2D6":{"chr":null,"callSource":"NONE","messages":[],
   "sourceDiplotypes":[{"allele1":{"name":"Unknown"},"allele2":{"name":"Unknown"},"phenotypes":["No Result"]}],"variants":[]},
  "TPMT":{"chr":"chr6","callSource":"OUTSIDE","messages":[{"message":"Called outside PharmCAT."}],
   "sourceDiplotypes":[{"allele1":{"name":"*1"},"allele2":{"name":"*3A"},"phenotypes":["Intermediate Metabolizer"]}],"variants":[]}},
 "unannotatedGeneCalls":[{"geneSymbol":"F5","chr":"chr1","callSource":"MATCHER","messages":[],
   "sourceDiplotypes":[{"allele1":{"name":"rs6025 reference (C)"},"allele2":{"name":"rs6025 reference (C)"},"phenotypes":[]}],
   "variants":[{"chromosome":"chr1","position":169549811,"dbSnpId":"rs6025","call":"C/C"}]}],
 "drugs":{
  "CPIC Guideline Annotation":{"clopidogrel":{"name":"clopidogrel","guidelines":[{"url":"https://example.org/g1","annotations":[
     {"drugRecommendation":"Use an alternative","classification":"Strong","lookupKey":[{"CYP2C19":"Intermediate Metabolizer"}]}]}]},
   "abacavir":{"name":"abacavir","guidelines":[{"url":"https://example.org/g2","annotations":[]}]}},
  "DPWG Guideline Annotation":{"clopidogrel":{"name":"clopidogrel","guidelines":[{"url":null,"annotations":[
     {"drugRecommendation":"Choose another drug","classification":null,"lookupKey":[],
      "genotypes":[{"diplotypes":[{"gene":"CYP2C19"}]}]}]}]}}}}"#);
    let extras = serde_json::json!({"prepared": {"positions": 6}});
    let common = pharmcat::esquema_comum(&p("x.report.json"), "S1", &["chrX".to_string()], extras, Some("x.report.html"), 2)
        .unwrap_or_else(|e| e);
    let v: serde_json::Value = serde_json::from_str(&common).unwrap_or(serde_json::Value::Null);
    let gene = |g: &str| v["genes"].as_array().and_then(|a| a.iter().find(|x| x["gene"] == g)).cloned().unwrap_or_default();
    check("schema version 1, engine pharmcat, its versions and the sample",
          v["schema_version"] == "1" && v["engine"] == "pharmcat" && v["engine_version"] == "3.4.0"
          && v["data_version"] == "PharmCAT data 2026-07-13" && v["sample"] == "S1", &common);
    let c19 = gene("CYP2C19");
    check("CYP2C19 called: diplotype, alleles, phenotype, one position read and one not",
          c19["called"] == true && c19["diplotype"] == "*1/*2" && c19["alleles"] == serde_json::json!(["*1", "*2"])
          && c19["phenotype"] == "Intermediate Metabolizer" && c19["positions_read"] == 1 && c19["positions_missing"] == 1
          && c19["not_called_reason"].is_null() && c19["source"] == "vcf", &c19.to_string());
    check("the evidence keeps each position, its rsID and the genotype ('./.' when not read)",
          c19["evidence"][0]["pos"] == 94781859 && c19["evidence"][0]["rsid"] == "rs4244285" && c19["evidence"][0]["genotype"] == "G/A"
          && c19["evidence"][1]["genotype"] == "./." && c19["evidence"][0]["in_callable"].is_null(), &c19.to_string());
    check("the label is composed when PharmCAT gives none; the activity score is read from text",
          gene("CYP2C9")["diplotype"] == "*1/*1" && gene("CYP2C9")["activity_score"] == 2.0, &gene("CYP2C9").to_string());
    check("two possible diplotypes: both, joined by OR, and the phenotype said once",
          gene("SLCO1B1")["diplotype"] == "*1/*15 OR *1/*45" && gene("SLCO1B1")["phenotype"] == "Decreased Function"
          && gene("SLCO1B1")["activity_score"].is_null(), &gene("SLCO1B1").to_string());
    check("a gene on a sequence missing from the file: not called, and the reason names the sequence",
          gene("G6PD")["called"] == false && gene("G6PD")["diplotype"].is_null()
          && gene("G6PD")["not_called_reason"].as_str().unwrap_or("").contains("chrX"), &gene("G6PD").to_string());
    check("CYP2D6: not called, and the reason says PharmCAT needs it as an outside call",
          gene("CYP2D6")["called"] == false && gene("CYP2D6")["not_called_reason"].as_str().unwrap_or("").contains("outside call"),
          &gene("CYP2D6").to_string());
    check("an outside call keeps its source and PharmCAT's message",
          gene("TPMT")["source"] == "outside" && gene("TPMT")["notes"] == serde_json::json!(["Called outside PharmCAT."]), &gene("TPMT").to_string());
    check("a gene PharmCAT calls without a guideline (F5, kept apart in its report) is there too, without a phenotype",
          gene("F5")["called"] == true && gene("F5")["diplotype"] == "rs6025 reference (C)/rs6025 reference (C)"
          && gene("F5")["phenotype"].is_null(), &gene("F5").to_string());
    let drugs = v["drugs"].as_array().cloned().unwrap_or_default();
    check("drug recommendations: CPIC with its strength and link, DPWG with the gene from the genotypes; none without a text",
          drugs.len() == 2 && drugs[0]["drug"] == "clopidogrel" && drugs[0]["source"] == "CPIC" && drugs[0]["strength"] == "Strong"
          && drugs[0]["url"] == "https://example.org/g1" && drugs[0]["genes"] == serde_json::json!(["CYP2C19"])
          && drugs[1]["source"] == "DPWG" && drugs[1]["genes"] == serde_json::json!(["CYP2C19"]) && drugs[1]["strength"].is_null(),
          &v["drugs"].to_string());
    check("warnings for the missing sequence and the matcher; report path, summary, and the run under x_pharmcat",
          v["warnings"].as_array().map_or(0, |a| a.len()) == 2 && v["report_html"] == "x.report.html"
          && v["text"].as_str().unwrap_or("").contains("5 of 7 genes called") && v["x_pharmcat"]["prepared"]["positions"] == 6
          && v["traits"] == serde_json::json!([]), &common);
    let e = pharmcat::esquema_comum(&p("x.report.tsv"), "S1", &[], serde_json::json!({}), None, 0).err().unwrap_or_default();
    check("a file that is not a PharmCAT JSON report: a clear error", e.contains("is not a PharmCAT JSON report"), &e);

    // PharmCAT itself, when available: on the synthetic person (whose positions
    // are not the real ones, so nothing is called) and, with PHARMCAT_PERSON
    // (+ PHARMCAT_BED), on a real genome.
    if let Ok(jar) = std::env::var("PHARMCAT_JAR") {
        println!("== PharmCAT itself ({jar})");
        let rr = pharmcat::rodar(&p("person.vcf"), &p("pc_out"), &jar, None, &op, "java", None, Some("me"));
        match rr {
            Ok(rr) => {
                check("PharmCAT ran: its version, and the calls table and HTML report written",
                      rr.versao.starts_with("PharmCAT")
                      && rr.arquivos.iter().any(|(k, f)| k == "calls_tsv" && std::path::Path::new(f).exists())
                      && rr.arquivos.iter().any(|(k, f)| k == "report_html" && std::path::Path::new(f).exists()),
                      &rr.versao);
                check("the run reports its command and timings", rr.comando[1] == "-jar" && rr.segundos_pharmcat > 0.0, "");
            }
            Err(e) => check("PharmCAT ran", false, &e),
        }
        let e = pharmcat::rodar(&p("person.vcf"), &p("pc_out2"), &jar, None, &op, "no_such_java_program", None, None)
            .err().unwrap_or_default();
        check("Java not found: the error says how to install it", e.contains("Java was not found") && e.contains("adoptium"), &e);
        if let Ok(person) = std::env::var("PHARMCAT_PERSON") {
            let op5 = pharmcat::Opcoes { amostra: None, chamavel: std::env::var("PHARMCAT_BED").ok(), ausente_e_ref: true };
            match pharmcat::rodar(&person, &p("pc_real"), &jar, None, &op5, "java", None, None) {
                Ok(rr) => {
                    let gene = |g: &str| rr.chamadas.iter().find(|c| c.campos.first().map_or(false, |(_, v)| v == g))
                        .and_then(|c| c.campos.iter().find(|(k, _)| k == "source_diplotype").map(|(_, v)| v.clone()));
                    println!("        {} genes; CYP2C19 {:?}, CYP2C9 {:?}, DPYD {:?}, SLCO1B1 {:?}; prepared: {} read of {}",
                             rr.chamadas.len(), gene("CYP2C19"), gene("CYP2C9"), gene("DPYD"), gene("SLCO1B1"),
                             rr.preparo.lidas(), rr.preparo.posicoes);
                    check("a real genome: PharmCAT calls the genes (CYP2C19 among them)", rr.chamadas.len() >= 15 && gene("CYP2C19").is_some(), "");
                    let common: serde_json::Value = rr.esquema_comum().ok().and_then(|s| serde_json::from_str(&s).ok()).unwrap_or_default();
                    let called = common["genes"].as_array().map_or(0, |a| a.iter().filter(|g| g["called"] == true).count());
                    let c19 = common["genes"].as_array().and_then(|a| a.iter().find(|g| g["gene"] == "CYP2C19")).cloned().unwrap_or_default();
                    println!("        common schema: {called} genes called, {} drug recommendations",
                             common["drugs"].as_array().map_or(0, |a| a.len()));
                    check("a real genome in the common schema: the called genes and CYP2C19 agree with the calls table",
                          called == rr.chamadas.len() && c19["diplotype"].as_str() == gene("CYP2C19").as_deref()
                          && common["drugs"].as_array().map_or(0, |a| a.len()) > 10, "");
                }
                Err(e) => check("PharmCAT on the real genome", false, &e),
            }
        }
    } else {
        println!("== PharmCAT itself: skipped (set PHARMCAT_JAR to the pharmcat-<version>-all.jar)");
    }

    let _ = std::fs::remove_dir_all(&dir);
    let (done, failed) = unsafe { (DONE, FAILED) };
    println!("\n{} of {} checks passed", done - failed, done);
    if failed > 0 {
        std::process::exit(1);
    }
}
