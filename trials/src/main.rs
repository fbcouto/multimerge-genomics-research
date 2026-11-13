#![allow(dead_code)]
// Validacao de src/externo.rs e src/intervalos.rs, isolada do resto do projeto.
// Cada resultado e comparado, linha a linha, com uma referencia independente:
// o sort estavel da std (ordenacao) ou forca bruta (fusao e intersecao).
// Uso: cargo run --release
#[path = "../../src/multimerge_core.rs"]
pub mod multimerge_core;
#[path = "../../src/externo.rs"]
mod externo;
#[path = "../../src/intervalos.rs"]
mod intervalos;

use externo::*;
use intervalos::*;
use std::path::PathBuf;

fn lcg(s: &mut u64) -> u64 {
    *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    *s >> 11
}

const CROMS: [&str; 8] = ["chr1", "chr2", "chr10", "chr21", "chrM", "chrX", "chrY", "chr3"];
// Ordens esperadas, escritas a mao: nao usam o codigo sob teste.
const NATURAL: [&str; 8] = ["chr1", "chr2", "chr3", "chr10", "chr21", "chrM", "chrX", "chrY"];
const LISTA: [&str; 3] = ["chrX", "chr10", "chr1"]; // o resto vai para o fim, em ordem natural

fn posto(ordem: &OrdemCromossomos, c: &str) -> (u32, usize, String) {
    let nat = NATURAL.iter().position(|x| *x == c).unwrap_or(99);
    match ordem {
        OrdemCromossomos::Lexica => (0, 0, c.to_string()),
        OrdemCromossomos::Natural | OrdemCromossomos::Cabecalho => (0, nat, String::new()),
        OrdemCromossomos::Lista(l) => match l.iter().position(|x| x == c) {
            Some(p) => (p as u32, 0, String::new()),
            None => (u32::MAX, nat, String::new()),
        },
    }
}

fn op(mem: usize, dir: &PathBuf, ordem: OrdemCromossomos) -> Opcoes {
    Opcoes { memoria_bytes: mem, dir_temp: Some(dir.clone()), ordem }
}

fn le(p: &str) -> Vec<String> {
    let mut r = abre_leitura(p).unwrap();
    let mut s = String::new();
    std::io::Read::read_to_string(&mut r, &mut s).unwrap();
    s.lines().map(str::to_string).collect()
}

fn col(l: &str, i: usize) -> &str {
    l.split("\t").nth(i).unwrap()
}

fn num(l: &str, i: usize) -> u64 {
    col(l, i).parse().unwrap()
}

type Iv = (String, u64, u64, String);

fn main() {
    let dir = std::env::temp_dir().join(format!("teste_v2_{}", std::process::id()));
    let lotes_dir = dir.join("lotes");
    std::fs::create_dir_all(&lotes_dir).unwrap();
    let p = |n: &str| dir.join(n).to_string_lossy().to_string();
    let mut s = 7u64;
    let mut ensaios = 0;
    let ordens = || {
        vec![
            OrdemCromossomos::Lexica,
            OrdemCromossomos::Natural,
            OrdemCromossomos::Lista(LISTA.iter().map(|x| x.to_string()).collect()),
        ]
    };
    let vazio = |d: &PathBuf| std::fs::read_dir(d).unwrap().next().is_none();

    // ---------------- BED: ordenacao
    for n in [0usize, 1, 3000, 60_000] {
        let linhas: Vec<String> = (0..n)
            .map(|i| {
                let c = CROMS[(lcg(&mut s) % 8) as usize];
                let ini = lcg(&mut s) % 500;
                format!("{c}\t{ini}\t{}\tr{i}\t{}", ini + lcg(&mut s) % 40, lcg(&mut s) % 9)
            })
            .collect();
        let mut txt = String::from("track name=teste\n#comentario\n");
        for l in &linhas {
            txt.push_str(l);
            txt.push_str("\n");
        }
        txt.push_str("chr1\tnao_numero\t5\n\nlinha_curta\n"); // 2 invalidas, 1 vazia
        std::fs::write(p("in.bed"), &txt).unwrap();
        for ordem in ordens() {
            let mut esp = linhas.clone();
            esp.sort_by(|a, b| {
                posto(&ordem, col(a, 0))
                    .cmp(&posto(&ordem, col(b, 0)))
                    .then(num(a, 1).cmp(&num(b, 1)))
                    .then(num(a, 2).cmp(&num(b, 2)))
            });
            for (mem, saida) in [(usize::MAX, "out.bed"), (1, "out.bed"), (1, "out.bed.gz")] {
                let r = ordenar_arquivos(
                    &[p("in.bed")],
                    &p(saida),
                    Formato::Bed,
                    Criterio::Posicao,
                    &op(mem, &lotes_dir, ordem.clone()),
                )
                .unwrap();
                let got = le(&p(saida));
                assert_eq!(&got[..2], ["track name=teste", "#comentario"], "cabecalho BED");
                assert!(got[2..] == esp[..], "BED difere: n={n} {ordem:?} mem={mem}");
                assert_eq!((r.registros, r.ignoradas), (n, 2));
                assert!(vazio(&lotes_dir), "sobrou temporario");
                ensaios += 1;
            }
        }
    }

    // ---------------- GFF
    {
        let n = 40_000;
        let linhas: Vec<String> = (0..n)
            .map(|i| {
                let c = CROMS[(lcg(&mut s) % 8) as usize];
                let ini = 1 + lcg(&mut s) % 300;
                format!("{c}\tfonte\tgene\t{ini}\t{}\t.\t+\t.\tID=g{i}", ini + lcg(&mut s) % 30)
            })
            .collect();
        // `###` espalhados (como num GFF do Ensembl) e uma diretiva no meio
        let mut txt = String::from("##gff-version 3\n");
        for (i, l) in linhas.iter().enumerate() {
            txt.push_str(l);
            txt.push_str("\n");
            if i % 7 == 0 {
                txt.push_str("###\n");
            }
            if i == 1000 {
                txt.push_str("##sequence-region chr2 1 999\n");
            }
        }
        std::fs::write(p("in.gff"), txt).unwrap();
        let nat = OrdemCromossomos::Natural;
        let mut esp = linhas.clone();
        // mesmo inicio: o mais longo (fim maior) primeiro -- pai antes do filho
        esp.sort_by(|a, b| {
            posto(&nat, col(a, 0))
                .cmp(&posto(&nat, col(b, 0)))
                .then(num(a, 3).cmp(&num(b, 3)))
                .then(num(b, 4).cmp(&num(a, 4)))
        });
        for mem in [usize::MAX, 1] {
            ordenar_arquivos(&[p("in.gff")], &p("out.gff"), Formato::Gff, Criterio::Posicao, &op(mem, &lotes_dir, nat.clone()))
                .unwrap();
            let got = le(&p("out.gff"));
            assert_eq!(&got[..2], ["##gff-version 3", "##sequence-region chr2 1 999"], "diretivas no topo, sem ###");
            assert!(got[2..] == esp[..], "GFF difere mem={mem}");
            ensaios += 1;
        }
        // gene, mRNA e exon com o mesmo inicio: pai antes do filho
        std::fs::write(
            p("familia.gff"),
            "chr1\tf\texon\t100\t200\t.\t+\t.\tID=e1\nchr1\tf\tgene\t100\t900\t.\t+\t.\tID=g\nchr1\tf\tmRNA\t100\t900\t.\t+\t.\tID=m\nchr1\tf\texon\t50\t60\t.\t+\t.\tID=e0\n",
        )
        .unwrap();
        ordenar_arquivos(&[p("familia.gff")], &p("out.gff"), Formato::Gff, Criterio::Posicao, &op(1 << 30, &lotes_dir, nat.clone())).unwrap();
        let ids: Vec<String> = le(&p("out.gff")).iter().map(|l| col(l, 8).to_string()).collect();
        assert_eq!(ids, ["ID=e0", "ID=g", "ID=m", "ID=e1"], "pai antes do filho");
        ensaios += 1;
        std::fs::write(p("fasta.gff"), "##gff-version 3\nchr1\tf\tgene\t1\t2\t.\t+\t.\tID=a\n##FASTA\n>chr1\nACGT\n").unwrap();
        assert!(ordenar_arquivos(&[p("fasta.gff")], &p("out.gff"), Formato::Gff, Criterio::Posicao, &op(1 << 30, &lotes_dir, nat))
            .is_err());
    }

    // ---------------- SAM: coordenada (ordem do @SQ) e nome
    {
        // ordem do cabecalho, de proposito nem lexica nem natural
        let sq = ["chr2", "chr10", "chr1", "chrM"];
        let n = 50_000;
        let linhas: Vec<String> = (0..n)
            .map(|i| {
                let k = lcg(&mut s) % 5;
                let (rn, pos) = if k == 4 { ("*", 0) } else { (sq[k as usize], 1 + lcg(&mut s) % 400) };
                format!("read{}\t0\t{rn}\t{pos}\t60\t4M\t*\t0\t0\tACGT\tIIII\tXI:i:{i}", lcg(&mut s) % 3000)
            })
            .collect();
        let cab = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr2\tLN:1000\n@SQ\tSN:chr10\tLN:1000\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chrM\tLN:1000\n@PG\tID:x\n";
        std::fs::write(p("in.sam"), format!("{cab}{}\n", linhas.join("\n"))).unwrap();
        let posto_sq = |c: &str| if c == "*" { 99 } else { sq.iter().position(|x| *x == c).unwrap() };
        let mut esp = linhas.clone();
        esp.sort_by(|a, b| posto_sq(col(a, 2)).cmp(&posto_sq(col(b, 2))).then(num(a, 3).cmp(&num(b, 3))));
        let mut esp_nome = linhas.clone();
        let n_read = |l: &str| col(l, 0)[4..].parse::<u64>().unwrap();
        esp_nome.sort_by(|a, b| n_read(a).cmp(&n_read(b)));
        let cabec = OrdemCromossomos::Cabecalho;
        for mem in [usize::MAX, 1] {
            ordenar_arquivos(&[p("in.sam")], &p("out.sam"), Formato::Sam, Criterio::Posicao, &op(mem, &lotes_dir, cabec.clone()))
                .unwrap();
            let got = le(&p("out.sam"));
            assert_eq!(got[0], "@HD\tVN:1.6\tSO:coordinate");
            assert_eq!(got[5], "@PG\tID:x");
            assert!(got[6..] == esp[..], "SAM coordenada difere mem={mem}");
            ordenar_arquivos(&[p("in.sam")], &p("out.sam"), Formato::Sam, Criterio::Nome, &op(mem, &lotes_dir, cabec.clone()))
                .unwrap();
            let got = le(&p("out.sam"));
            assert_eq!(got[0], "@HD\tVN:1.6\tSO:queryname");
            assert!(got[6..] == esp_nome[..], "SAM nome difere mem={mem}");
            ensaios += 2;
        }
        // SAM sem @HD ganha um
        std::fs::write(p("sem_hd.sam"), "@SQ\tSN:chr1\tLN:9\nr1\t0\tchr1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\n").unwrap();
        ordenar_arquivos(&[p("sem_hd.sam")], &p("out.sam"), Formato::Sam, Criterio::Posicao, &op(1 << 30, &lotes_dir, cabec))
            .unwrap();
        assert_eq!(le(&p("out.sam"))[0], "@HD\tVN:1.6\tSO:coordinate");
    }

    // ---------------- VCF: varios arquivos, ordem do ##contig, QUAL, criterio invalido
    {
        let cab = "##fileformat=VCFv4.2\n##contig=<ID=chr10,length=5>\n##contig=<ID=chr2>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let quals = ["30", "30.5", "99", ".", "0", "12.25", "NaN", "60"];
        let gera = |n: usize, s: &mut u64, tag: &str| -> Vec<String> {
            (0..n)
                .map(|i| {
                    let c = ["chr2", "chr10", "chr1"][(lcg(s) % 3) as usize];
                    let q = quals[(lcg(s) % 8) as usize];
                    format!("{c}\t{}\t{tag}{i}\tA\tT\t{q}\tPASS\t.", lcg(s) % 200)
                })
                .collect()
        };
        let (l1, l2) = (gera(30_000, &mut s, "x"), gera(20_000, &mut s, "y"));
        std::fs::write(p("in1.vcf"), format!("{cab}{}\n", l1.join("\n"))).unwrap();
        // segundo arquivo: CRLF, cabecalho proprio que nao deve sair
        std::fs::write(p("in2.vcf"), format!("##outro\r\n{}\r\n", l2.join("\r\n"))).unwrap();
        let todas: Vec<String> = l1.iter().chain(l2.iter()).cloned().collect();
        let posto_c = |c: &str| ["chr10", "chr2", "chr1"].iter().position(|x| *x == c).unwrap();
        let mut esp = todas.clone();
        esp.sort_by(|a, b| posto_c(col(a, 0)).cmp(&posto_c(col(b, 0))).then(num(a, 1).cmp(&num(b, 1))));
        let q = |l: &str| col(l, 5).parse::<f32>().unwrap_or(0.0);
        let mut esp_q = todas.clone();
        esp_q.sort_by(|a, b| q(b).total_cmp(&q(a)));
        let entradas = [p("in1.vcf"), p("in2.vcf")];
        for mem in [usize::MAX, 1] {
            let r = ordenar_arquivos(&entradas, &p("out.vcf"), Formato::Vcf, Criterio::Posicao, &op(mem, &lotes_dir, OrdemCromossomos::Cabecalho))
                .unwrap();
            let got = le(&p("out.vcf"));
            assert_eq!(got[3], "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
            assert!(got[4..] == esp[..], "VCF por cabecalho difere mem={mem}");
            assert_eq!(r.registros, todas.len());
            ordenar_arquivos(&entradas, &p("out.vcf"), Formato::Vcf, Criterio::Qualidade, &op(mem, &lotes_dir, OrdemCromossomos::Lexica))
                .unwrap();
            assert!(le(&p("out.vcf"))[4..] == esp_q[..], "VCF por QUAL difere mem={mem}");
            ensaios += 2;
        }
        assert!(ordenar_arquivos(&entradas, &p("out.vcf"), Formato::Bed, Criterio::Nome, &op(1, &lotes_dir, OrdemCromossomos::Lexica)).is_err());
        // sem cabecalho: sai o minimo
        std::fs::write(p("nu.vcf"), "chr1\t5\tx\tA\tT\t1\tPASS\t.\n").unwrap();
        ordenar_arquivos(&[p("nu.vcf")], &p("out.vcf"), Formato::Vcf, Criterio::Posicao, &op(1 << 30, &lotes_dir, OrdemCromossomos::Lexica))
            .unwrap();
        assert_eq!(le(&p("out.vcf"))[0], "##fileformat=VCFv4.2");
        assert!(ordenar_arquivos(&["nao_existe.vcf".to_string()], &p("out.vcf"), Formato::Vcf, Criterio::Posicao, &op(1, &lotes_dir, OrdemCromossomos::Lexica)).is_err());
    }

    // ---------------- BED: fundir e intersecao, contra forca bruta
    for (na, nb, larg) in [(0usize, 50usize, 30u64), (2000, 0, 30), (2500, 2500, 12), (4000, 3000, 60), (3000, 4000, 1)] {
        let gera = |n: usize, s: &mut u64, tag: &str| -> Vec<Iv> {
            (0..n)
                .map(|i| {
                    let c = CROMS[(lcg(s) % 4) as usize].to_string();
                    let ini = lcg(s) % 3000;
                    (c, ini, ini + 1 + lcg(s) % larg, format!("{tag}{i}"))
                })
                .collect()
        };
        let (a, b) = (gera(na, &mut s, "a"), gera(nb, &mut s, "b"));
        let grava = |nome: &str, v: &[Iv]| {
            let t: String = v.iter().map(|(c, i, f, t)| format!("{c}\t{i}\t{f}\t{t}\n")).collect();
            std::fs::write(p(nome), t).unwrap();
        };
        grava("a.bed", &a);
        grava("b.bed", &b);
        let ordem = OrdemCromossomos::Natural;
        let chave = |x: &Iv| (posto(&ordem, &x.0), x.1, x.2);

        // fundir
        for dist in [0u64, 7] {
            let mut ord = a.clone();
            ord.sort_by(|x, y| chave(x).cmp(&chave(y)));
            let mut esp: Vec<(String, u64, u64)> = Vec::new();
            for (c, i, f, _) in &ord {
                match esp.last_mut() {
                    Some((ca, _, fa)) if ca == c && *i <= *fa + dist => *fa = (*fa).max(*f),
                    _ => esp.push((c.clone(), *i, *f)),
                }
            }
            let esp_txt: Vec<String> = esp.iter().map(|(c, i, f)| format!("{c}\t{i}\t{f}")).collect();
            for mem in [usize::MAX, 1] {
                let (lidos, gravados) = bed_fundir(&[p("a.bed")], &p("m.bed"), dist, &op(mem, &lotes_dir, ordem.clone())).unwrap();
                assert!(le(&p("m.bed")) == esp_txt, "fundir difere na={na} dist={dist} mem={mem}");
                assert_eq!((lidos, gravados), (na, esp.len()));
                ensaios += 1;
            }
        }

        // intersecao: forca bruta O(|A|*|B|), com A e B em ordem
        let mut ao = a.clone();
        ao.sort_by(|x, y| chave(x).cmp(&chave(y)));
        let mut bo = b.clone();
        bo.sort_by(|x, y| chave(x).cmp(&chave(y)));
        for modo in [ModoIntersecao::Regioes, ModoIntersecao::Unicos, ModoIntersecao::Ausentes] {
            let mut esp: Vec<String> = Vec::new();
            for (c, i, f, t) in &ao {
                let mut tocou = false;
                for (cb, ib, fb, _) in &bo {
                    if cb == c && i.max(ib) < f.min(fb) {
                        tocou = true;
                        if modo == ModoIntersecao::Regioes {
                            esp.push(format!("{c}\t{}\t{}\t{t}", i.max(ib), f.min(fb)));
                        }
                    }
                }
                if (modo == ModoIntersecao::Unicos && tocou) || (modo == ModoIntersecao::Ausentes && !tocou) {
                    esp.push(format!("{c}\t{i}\t{f}\t{t}"));
                }
            }
            for mem in [usize::MAX, 1] {
                let n = bed_intersecao(&[p("a.bed")], &[p("b.bed")], &p("x.bed"), modo, &op(mem, &lotes_dir, ordem.clone())).unwrap();
                let got = le(&p("x.bed"));
                assert!(got == esp, "intersecao difere na={na} nb={nb} {modo:?} mem={mem}: {} x {}", got.len(), esp.len());
                assert_eq!(n, esp.len());
                assert!(vazio(&lotes_dir), "sobrou temporario");
                ensaios += 1;
            }
        }
    }

    // ---------------- erro silencioso de nomenclatura: chr1 contra 1
    {
        std::fs::write(p("com_chr.bed"), "chr1\t10\t50\ta\nchr2\t5\t9\tb\n").unwrap();
        std::fs::write(p("sem_chr.bed"), "1\t20\t60\tx\n2\t1\t7\ty\n").unwrap();
        for modo in [ModoIntersecao::Regioes, ModoIntersecao::Unicos, ModoIntersecao::Ausentes] {
            let e = bed_intersecao(&[p("com_chr.bed")], &[p("sem_chr.bed")], &p("x.bed"), modo, &op(1 << 30, &lotes_dir, OrdemCromossomos::Natural))
                .unwrap_err();
            assert!(e.contains("`chr1`") && e.contains("`1`") && e.contains("`chr` prefix"), "mensagem: {e}");
            assert!(!std::path::Path::new(&p("x.bed")).exists(), "saida enganosa ficou no disco");
            assert!(vazio(&lotes_dir));
        }
        // nomes diferentes sem ser questao de prefixo: a dica muda
        std::fs::write(p("outro.bed"), "scaffold_7\t20\t60\tx\n").unwrap();
        let e = bed_intersecao(&[p("sem_chr.bed")], &[p("outro.bed")], &p("x.bed"), ModoIntersecao::Unicos, &op(1 << 30, &lotes_dir, OrdemCromossomos::Natural))
            .unwrap_err();
        assert!(e.contains("version of the reference genome"), "mensagem: {e}");
        // um cromossomo em comum basta: sem erro
        std::fs::write(p("meio.bed"), "chr2\t1\t7\tz\n3\t1\t2\tw\n").unwrap();
        assert_eq!(bed_intersecao(&[p("com_chr.bed")], &[p("meio.bed")], &p("x.bed"), ModoIntersecao::Unicos, &op(1 << 30, &lotes_dir, OrdemCromossomos::Natural)).unwrap(), 1);
        ensaios += 5;
    }

    // ---------------- natural_cmp: ordem total sobre nomes dificeis
    let nomes = ["chr1", "chr01", "chr2", "chr10", "chr10_alt", "chrX", "chr", "1", "10", "2", "MT", "scaffold_9", "scaffold_10", ""];
    for a in nomes {
        for b in nomes {
            assert_eq!(natural_cmp(a, b), natural_cmp(b, a).reverse(), "antissimetria {a} {b}");
            for c in nomes {
                if natural_cmp(a, b).is_le() && natural_cmp(b, c).is_le() {
                    assert!(natural_cmp(a, c).is_le(), "transitividade {a} {b} {c}");
                }
            }
        }
    }
    assert!(natural_cmp("chr2", "chr10").is_lt() && natural_cmp("scaffold_9", "scaffold_10").is_lt());
    assert_eq!(Formato::pelo_nome("x/y.vcf.gz"), Some(Formato::Vcf));
    assert_eq!(Formato::pelo_nome("a.GTF"), Some(Formato::Gff));

    std::fs::remove_dir_all(&dir).unwrap();
    println!("{ensaios} tests identical to the reference (BED, GFF, SAM, VCF, merge, intersect)");
}
