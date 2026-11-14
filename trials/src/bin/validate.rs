#![allow(dead_code)]
// Ensaios do validador (src/validador.rs) e dos nomes de sequencia (src/nomes.rs).
// Cada defeito e plantado numa linha conhecida; o ensaio confere o codigo do
// problema, a gravidade e a LINHA apontada.
// Uso: cargo run --release --bin validate [arquivo real]...
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
#[path = "../../../src/validador.rs"]
mod validador;

use externo::Formato;
use nomes::*;
use validador::*;

fn acha<'a>(r: &'a Relatorio, codigo: &str) -> &'a Problema {
    r.problemas.iter().find(|p| p.codigo == codigo).unwrap_or_else(|| panic!("missing problem `{codigo}` in:\n{}", r.resumo()))
}

fn sem(r: &Relatorio, codigo: &str) {
    assert!(r.problemas.iter().all(|p| p.codigo != codigo), "there should be no `{codigo}`:\n{}", r.resumo());
}

fn main() {
    let dir = std::env::temp_dir().join(format!("teste_valida_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = |n: &str| dir.join(n).to_string_lossy().to_string();
    let grava = |n: &str, t: &str| {
        std::fs::write(p(n), t).unwrap();
        p(n)
    };
    let valida = |n: &str, t: &str| validar_arquivo(&grava(n, t), None, None).unwrap();
    let mut ensaios = 0;
    let mut jsons: Vec<String> = Vec::new();

    // ------------------------------------------------ VCF limpo
    let cab = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n##contig=<ID=chr2,length=242193529>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let r = valida("bom.vcf", &format!("{cab}chr1\t10\ta\tA\tT\t30\tPASS\t.\nchr1\t10\tb\tG\tC,T\t.\tPASS\t.\nchr2\t5\tc\tAT\t<DEL>\t9.5\tPASS\t.\n"));
    assert!(r.valido() && r.problemas.is_empty(), "{}", r.resumo());
    assert_eq!((r.linhas, r.registros, r.sequencias, r.nomenclatura, r.ordenado), (7, 3, 2, "chr", Some(true)));
    jsons.push(r.json());
    ensaios += 1;

    // ------------------------------------------------ VCF com um defeito por linha
    let ruim = [
        "##source=sem_fileformat",                    // 1  aviso: primeira linha
        "##contig=<ID=chr1,length=10>",               // 2
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", // 3
        "chr1\t50\ta\tA\tT\t30\tPASS\t.",             // 4
        "chr1 60 b A T 30 PASS .",                    // 5  erro: espacos no lugar de tab
        "chr1\tabc\tc\tA\tT\t30\tPASS\t.",            // 6  erro: POS
        "chr1\t70\td\tAXG\tT\t30\tPASS\t.",           // 7  erro: REF
        "chr1\t80\te\tA\tT\talto\tPASS\t.",           // 8  erro: QUAL
        "",                                           // 9  aviso: linha vazia
        "chr1\t20\tf\tA\tT\t1\tPASS\t.",              // 10 info: desordenado (20 depois de 80)
        "1\t90\tg\tA\tT\t1\tPASS\t.",                 // 11 erro: mistura chr1 e 1; aviso: contig nao declarado
        "chr1\t95\th\tA\tT\t1\tPASS",                 // 12 erro: 7 colunas
        "chr1\t99\ti\tA\tT?\t1\tPASS\t.",             // 13 erro: ALT
    ];
    let r = valida("ruim.vcf", &format!("{}\r\n", ruim.join("\r\n")));
    for (codigo, gravidade, linha) in [
        ("vcf_missing_fileformat", Gravidade::Aviso, 1),
        ("vcf_columns", Gravidade::Erro, 5),
        ("vcf_pos", Gravidade::Erro, 6),
        ("vcf_ref", Gravidade::Erro, 7),
        ("vcf_qual", Gravidade::Erro, 8),
        ("blank_line", Gravidade::Aviso, 9),
        ("unsorted", Gravidade::Info, 10),
        ("vcf_contig", Gravidade::Aviso, 11),
        ("mixed_names", Gravidade::Erro, 4),
        ("vcf_alt", Gravidade::Erro, 13),
    ] {
        let q = acha(&r, codigo);
        assert_eq!((q.gravidade, q.linha), (gravidade, linha), "{codigo}:\n{}", r.resumo());
    }
    assert!(acha(&r, "vcf_columns").mensagem.contains("TABS"), "a linha 5 e de espacos");
    assert_eq!(acha(&r, "vcf_columns").ocorrencias, 2, "linhas 5 e 12");
    assert_eq!(acha(&r, "crlf").ocorrencias, 13);
    assert!(acha(&r, "mixed_names").mensagem.contains("`chr1` (line 4) and `1` (line 11)"));
    assert_eq!((r.valido(), r.ordenado, r.nomenclatura), (false, Some(false), "mixed"));
    // do mais grave para o menos grave
    assert!(r.problemas.windows(2).all(|w| w[0].gravidade <= w[1].gravidade));
    jsons.push(r.json());
    ensaios += 1;

    // sem a linha #CHROM; e validacao parcial
    let r = valida("sem_chrom.vcf", "##fileformat=VCFv4.2\nchr1\t5\ta\tA\tT\t1\tPASS\t.\n");
    assert_eq!(acha(&r, "vcf_missing_chrom").linha, 2);
    let r = validar_arquivo(&p("ruim.vcf"), None, Some(6)).unwrap();
    assert!(!r.completo && r.linhas == 6);
    sem(&r, "vcf_ref");
    ensaios += 2;

    // ------------------------------------------------ BED
    let r = valida("x.bed", "track name=t\nchr1\t10\t20\tn\nchr1\t30\t25\tn\nchr1\t-5\t9\nchr2 1 2\nchr2\t7\n");
    assert_eq!(acha(&r, "bed_inverted").linha, 3);
    assert_eq!(acha(&r, "bed_coordinate").linha, 4);
    let c = acha(&r, "bed_columns");
    assert_eq!((c.linha, c.ocorrencias), (5, 2));
    assert!(c.mensagem.contains("TABS"));
    let r = valida("ok.bed", "chr1\t0\t10\nchr1\t5\t7\nchrX\t1\t2\n");
    assert!(r.valido() && r.ordenado == Some(true));
    ensaios += 2;

    // ------------------------------------------------ GFF
    let r = valida(
        "x.gff",
        "##gff-version 3\nchr1\tf\tgene\t10\t20\t.\t+\t.\tID=a\nchr1\tf\tgene\t0\t5\t.\t+\t.\tID=b\nchr1\tf\tCDS\t30\t40\t.\tx\t.\tID=c\nchr1\tf\tgene\t9\t3\t.\t-\t.\tID=d\nchr1\tf\tgene\t50\t60\t.\t+\t.\n##FASTA\n>chr1\nACGT 123 lixo\n",
    );
    assert_eq!(acha(&r, "gff_zero_based").linha, 3);
    assert_eq!(acha(&r, "gff_strand").linha, 4);
    assert_eq!(acha(&r, "gff_phase").linha, 4);
    assert_eq!(acha(&r, "gff_inverted").linha, 5);
    assert_eq!(acha(&r, "gff_columns").linha, 6);
    assert_eq!((acha(&r, "gff_fasta").linha, acha(&r, "gff_fasta").gravidade), (7, Gravidade::Info));
    assert_eq!(r.registros, 5, "a secao FASTA nao conta como registro");
    ensaios += 1;

    // ------------------------------------------------ SAM
    let sam = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t50\t60\t4M\t*\t0\t0\tACGT\tIIII\nr2\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIII\nr3\t0\tchr9\t30\t60\t4M\t*\t0\t0\tACGT\tIIII\nr4\t0\tchr1\t40\nr5\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";
    let r = valida("x.sam", sam);
    assert_eq!(acha(&r, "sam_quality").linha, 4);
    assert_eq!(acha(&r, "sam_rname").linha, 5);
    assert_eq!(acha(&r, "sam_columns").linha, 6);
    assert_eq!(acha(&r, "sam_so_false").gravidade, Gravidade::Erro);
    let r = valida("sem_sq.sam", "r1\t0\tchr1\t50\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
    assert_eq!(acha(&r, "sam_missing_sq").gravidade, Gravidade::Aviso);
    assert!(r.valido());
    ensaios += 2;

    // ------------------------------------------------ FASTA: "erro na linha 45: caractere invalido"
    let mut fa = String::from(">seq1 descricao\n");
    for _ in 0..43 {
        fa.push_str("ACGTACGTNNACGT\n");
    }
    fa.push_str("AC7TACGT\n"); // linha 45
    fa.push_str(">seq2\n>seq1\nACGT\n"); // seq2 vazia (linha 46); seq1 repetido (linha 47)
    let r = valida("x.fasta", &fa);
    let c = acha(&r, "fasta_character");
    assert_eq!(c.linha, 45);
    assert!(c.mensagem.contains("`7`") && c.mensagem.contains("column 3"), "{}", c.mensagem);
    assert_eq!(acha(&r, "fasta_empty").linha, 46);
    assert_eq!(acha(&r, "fasta_duplicate_id").linha, 47);
    assert_eq!((r.registros, r.ordenado, r.nomenclatura), (3, None, "-"));
    let r = valida("y.fa", "ACGT\n>a\nAC GT\n");
    assert_eq!(acha(&r, "fasta_missing_header").linha, 1);
    assert!(acha(&r, "fasta_character").mensagem.contains("a space"));
    let r = valida("z.fa.gz.fasta", ">p1\nMKV*LLA-\n>p2\nacgtn\n");
    assert!(r.valido() && r.problemas.is_empty(), "proteina, parada, lacuna e minusculas valem: {}", r.resumo());
    jsons.push(valida("x2.fasta", &fa).json());
    ensaios += 3;

    // ------------------------------------------------ FASTQ
    let r = valida("x.fastq", "@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIII\n@r3\nAC1T\n-\nIIII\n@r4\nAC\n");
    assert_eq!(acha(&r, "fastq_length").linha, 8);
    assert_eq!(acha(&r, "fastq_base").linha, 10);
    assert_eq!(acha(&r, "fastq_separator").linha, 11);
    let t = acha(&r, "fastq_truncated");
    assert!(t.linha == 14 && t.mensagem.contains("2 line(s) missing"), "{}", t.mensagem);
    let r = valida("y.fq", "@r1\nACGT\n+\nhhhh\n@r2\nAC\n+\nef\n");
    assert_eq!(acha(&r, "fastq_phred64").gravidade, Gravidade::Aviso);
    let r = valida("z.fq", "@r1\nACGT\n+\nII#5\n\n");
    assert!(r.valido() && r.problemas.is_empty(), "{}", r.resumo());
    let r = valida("w.fq", "r1\nACGT\n+\nIIII\n");
    assert_eq!(acha(&r, "fastq_header").linha, 1);
    ensaios += 4;

    // ------------------------------------------------ arquivo vazio, so cabecalho, binario, formato desconhecido
    assert_eq!(acha(&valida("v.vcf", ""), "empty").gravidade, Gravidade::Erro);
    assert_eq!(acha(&valida("c.vcf", cab), "no_records").gravidade, Gravidade::Aviso);
    std::fs::write(p("b.bed"), [b'c', b'h', b'r', 0xff, 0xfe, b'\t', b'1', b'\t', b'2', b'\n']).unwrap();
    assert_eq!(acha(&validar_arquivo(&p("b.bed"), None, None).unwrap(), "not_text").linha, 1);
    assert!(validar_arquivo(&grava("dados.txt", "x"), None, None).is_err());
    assert!(validar_arquivo(&p("dados.txt"), Some(Tipo::Bed), None).is_ok());
    assert!(validar_arquivo(&p("nao_existe.vcf"), None, None).is_err());
    assert_eq!(Tipo::pelo_nome("a/b.FASTQ.GZ"), Some(Tipo::Fastq));
    ensaios += 4;

    // ------------------------------------------------ nomes: reconhecer e converter
    for (n, e) in [("1", true), ("22", true), ("X", true), ("MT", true), ("2L", true), ("XVI", true), ("GL000192.1", false), ("Un_KI270742v1", false), ("NC_000001.11", false), ("", false), ("1_random", false)] {
        assert_eq!(parece_cromossomo(n), e, "{n}");
    }
    for (de, estilo, para) in [
        ("1", Estilo::ComChr, Some("chr1")),
        ("MT", Estilo::ComChr, Some("chrM")),
        ("chrMT", Estilo::ComChr, Some("chrM")),
        ("chr1", Estilo::ComChr, None),
        ("chr1", Estilo::SemChr, Some("1")),
        ("chrM", Estilo::SemChr, Some("MT")),
        ("M", Estilo::SemChr, Some("MT")),
        ("X", Estilo::SemChr, None),
        ("chrUn_KI270742v1", Estilo::SemChr, None),
        ("GL000192.1", Estilo::ComChr, None),
    ] {
        assert_eq!(converte_nome(de, estilo).as_deref(), para, "{de}");
    }
    assert_eq!((canonico("chr7"), canonico("7"), canonico("chrM"), canonico("M"), canonico("scaf_1")), ("7".into(), "7".into(), "MT".into(), "MT".into(), "scaf_1".into()));
    ensaios += 1;

    // ------------------------------------------------ normalizar: os quatro formatos, cabecalho junto
    let vcf = format!("{cab}chr1\t10\ta\tA\tT\t30\tPASS\t.\nchrM\t3\tb\tA\tT\t1\tPASS\t.\nchrUn_x\t3\tc\tA\tT\t1\tPASS\t.\n");
    grava("n.vcf", &vcf);
    let r = normalizar_cromossomos(&p("n.vcf"), &p("n_sem.vcf"), Formato::Vcf, Estilo::SemChr, None).unwrap();
    let sem_chr = std::fs::read_to_string(p("n_sem.vcf")).unwrap();
    assert!(sem_chr.contains("##contig=<ID=1,length=248956422>") && sem_chr.contains("\n1\t10\ta") && sem_chr.contains("\nMT\t3\tb") && sem_chr.contains("\nchrUn_x\t3\tc"));
    assert_eq!((r.registros, r.alterados, r.mantidos.clone()), (3, 2, vec!["chrUn_x".to_string()]));
    assert!(r.convertidos.contains(&("chrM".to_string(), "MT".to_string())));
    // ida e volta devolve o original
    normalizar_cromossomos(&p("n_sem.vcf"), &p("n_volta.vcf"), Formato::Vcf, Estilo::ComChr, None).unwrap();
    assert_eq!(std::fs::read_to_string(p("n_volta.vcf")).unwrap(), vcf);
    jsons.push(r.json());
    // o arquivo que misturava convencoes deixa de misturar
    normalizar_cromossomos(&p("ruim.vcf"), &p("ruim_chr.vcf"), Formato::Vcf, Estilo::ComChr, None).unwrap();
    let r = validar_arquivo(&p("ruim_chr.vcf"), None, None).unwrap();
    sem(&r, "mixed_names");
    sem(&r, "crlf");
    assert_eq!(r.nomenclatura, "chr");

    grava("n.sam", "@HD\tVN:1.6\n@SQ\tSN:1\tLN:9\n@SQ\tSN:MT\tLN:5\nr1\t0\t1\t5\t60\t4M\t=\t9\t0\tACGT\tIIII\nr2\t0\tMT\t2\t60\t4M\t1\t7\t0\tACGT\tIIII\nr3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n");
    normalizar_cromossomos(&p("n.sam"), &p("n_chr.sam"), Formato::Sam, Estilo::ComChr, None).unwrap();
    assert_eq!(
        std::fs::read_to_string(p("n_chr.sam")).unwrap(),
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:9\n@SQ\tSN:chrM\tLN:5\nr1\t0\tchr1\t5\t60\t4M\t=\t9\t0\tACGT\tIIII\nr2\t0\tchrM\t2\t60\t4M\tchr1\t7\t0\tACGT\tIIII\nr3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
    );
    grava("n.gff", "##gff-version 3\n##sequence-region 1 1 1000\n1\tf\tgene\t1\t9\t.\t+\t.\tID=a\n");
    normalizar_cromossomos(&p("n.gff"), &p("n_chr.gff"), Formato::Gff, Estilo::ComChr, None).unwrap();
    assert_eq!(std::fs::read_to_string(p("n_chr.gff")).unwrap(), "##gff-version 3\n##sequence-region chr1 1 1000\nchr1\tf\tgene\t1\t9\t.\t+\t.\tID=a\n");
    // mapa explicito: converte o que a regra nao sabe (acesso RefSeq) e tem precedencia
    grava("mapa.txt", "# de para\nNC_000001.11\tchr1\nGL1\n2 dois\n");
    grava("n.bed", "NC_000001.11\t1\t2\n2\t3\t4\nGL1\t5\t6\n");
    let r = normalizar_cromossomos(&p("n.bed"), &p("n_chr.bed"), Formato::Bed, Estilo::ComChr, Some(&p("mapa.txt"))).unwrap();
    assert_eq!(std::fs::read_to_string(p("n_chr.bed")).unwrap(), "chr1\t1\t2\ndois\t3\t4\nGL1\t5\t6\n");
    assert_eq!(r.mantidos, vec!["GL1".to_string()]);
    ensaios += 5;

    // ------------------------------------------------ referencia: GRCh37 x GRCh38 pelos tamanhos
    let v37 = grava("g37.vcf", "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n##contig=<ID=2,length=243199373>\n##contig=<ID=MT,length=16569>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t5\ta\tA\tT\t1\tPASS\t.\n");
    let v38 = grava("g38.vcf", cab);
    let s38 = grava("g38.sam", "@HD\tVN:1.6\n@SQ\tSN:1\tLN:248956422\tAS:GRCh38\n@SQ\tSN:2\tLN:242193529\n@SQ\tSN:MT\tLN:16569\n");
    let fai37 = grava("hg19.fai", "chr1\t249250621\t52\t60\t61\nchr2\t243199373\t253404903\t60\t61\nchrM\t16569\t1\t60\t61\n");
    assert_eq!(referencia_do_arquivo(&v37).unwrap().montagem, Some("GRCh37 (hg19)"));
    assert_eq!(referencia_do_arquivo(&v38).unwrap().montagem, Some("GRCh38 (hg38)"));
    let r = referencia_do_arquivo(&s38).unwrap();
    assert_eq!((r.montagem, r.declarada.as_deref(), r.sequencias.len()), (Some("GRCh38 (hg38)"), Some("GRCh38"), 3));
    assert_eq!(referencia_do_arquivo(&fai37).unwrap().sequencias.len(), 3);

    let c = comparar_referencias(&v37, &v38).unwrap();
    assert_eq!((c.compativeis, c.diferentes.len(), c.iguais), (Some(false), 2, 0));
    assert!(c.veredito.contains("DIFFERENT") && c.veredito.contains("249250621") && c.veredito.contains("GRCh37") && c.veredito.contains("GRCh38"), "{}", c.veredito);
    jsons.push(c.json());
    // mesma referencia, convencoes de nome diferentes
    let c = comparar_referencias(&v38, &s38).unwrap();
    assert_eq!((c.compativeis, c.iguais, c.convencoes_diferentes, c.so_em_b), (Some(true), 2, true, 1));
    assert!(c.veredito.contains("normalize_chromosomes"), "{}", c.veredito);
    // VCF x .fai da mesma montagem
    let c = comparar_referencias(&v37, &fai37).unwrap();
    assert_eq!((c.compativeis, c.iguais, c.convencoes_diferentes), (Some(true), 3, true));
    // quem nao declara tamanho: sem veredito
    let c = comparar_referencias(&v37, &grava("mudo.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")).unwrap();
    assert_eq!(c.compativeis, None);
    jsons.push(referencia_do_arquivo(&s38).unwrap().json());
    ensaios += 4;

    // os JSON gerados, um por linha, para conferir com um analisador de verdade
    std::fs::write(p("saidas.jsonl"), jsons.join("\n")).unwrap();
    println!("JSONL: {}", p("saidas.jsonl"));
    println!("{ensaios} validator and naming tests passed");

    // ------------------------------------------------ arquivos reais, se informados
    for arq in std::env::args().skip(1) {
        let t = std::time::Instant::now();
        let limite = if arq.contains(".vcf") { Some(3_000_000) } else { None };
        let r = validar_arquivo(&arq, None, limite).unwrap();
        print!("\n[{:.1} s] {}", t.elapsed().as_secs_f64(), r.resumo());
        let t = std::time::Instant::now();
        let f = referencia_do_arquivo(&arq).unwrap();
        println!("[{:.1} s] reference: assembly={:?} declared={:?} sequences with a length={}", t.elapsed().as_secs_f64(), f.montagem, f.declarada, f.sequencias.len());
    }
    let reais: Vec<String> = std::env::args().skip(1).filter(|a| a.contains(".vcf")).collect();
    if reais.len() == 2 {
        println!("\ncomparison: {}", comparar_referencias(&reais[0], &reais[1]).unwrap().veredito);
    }
}
