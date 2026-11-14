// Validation of the container formats: BGZF (src/bgzf.rs), the tabix and BAI
// indexes (src/indice.rs) and BAM (src/bam.rs).
//
//   cargo run --release --bin formats
//
// Everything is checked against something independent of the code under test:
// BGZF files are read back with flate2's gzip decoder, region queries are
// compared with a linear scan of the same records, and BAM round trips are
// compared with the SAM text they came from.

#[path = "../../../src/multimerge_core.rs"]
pub mod multimerge_core;
#[path = "../../../src/bgzf.rs"]
mod bgzf;
#[path = "../../../src/bam.rs"]
mod bam;
#[path = "../../../src/externo.rs"]
#[allow(dead_code)]
mod externo;
#[path = "../../../src/indice.rs"]
mod indice;
#[path = "../../../src/intervalos.rs"]
#[allow(dead_code)]
mod intervalos;

use bgzf::{BgzfParLeitor, BgzfReader, BgzfWriter};
use std::io::{BufRead, Read, Write};

// ---------------------------------------------------------------------
//  Harness
// ---------------------------------------------------------------------

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 11
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
    fn chance(&mut self, pct: u64) -> bool {
        self.below(100) < pct
    }
}

struct Tests {
    n: usize,
    dir: std::path::PathBuf,
}

impl Tests {
    fn ok(&mut self, name: &str) {
        self.n += 1;
        println!("ok {:>3}  {name}", self.n);
    }
    fn path(&self, name: &str) -> String {
        self.dir.join(name).to_string_lossy().to_string()
    }
}

fn write_file(path: &str, data: &[u8]) {
    std::fs::write(path, data).unwrap_or_else(|e| panic!("cannot write {path}: {e}"));
}

fn gunzip(path: &str) -> Vec<u8> {
    let mut out = Vec::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_end(&mut out)
        .unwrap_or_else(|e| panic!("gzip decoder failed on {path}: {e}"));
    out
}

fn first_difference(a: &[u8], b: &[u8]) -> String {
    let i = a.iter().zip(b).position(|(x, y)| x != y).unwrap_or(a.len().min(b.len()));
    format!("lengths {} and {}, first difference at byte {i}", a.len(), b.len())
}

fn same_lines(what: &str, got: &[String], want: &[String]) {
    for (i, (g, w)) in got.iter().zip(want).enumerate() {
        assert!(g == w, "{what}: line {} differs\n   got: {g}\n  want: {w}", i + 1);
    }
    assert!(got.len() == want.len(), "{what}: {} lines, expected {}", got.len(), want.len());
}

// ---------------------------------------------------------------------
//  BGZF
// ---------------------------------------------------------------------

/// Text-like data that compresses, with some incompressible stretches.
fn sample_data(rng: &mut Lcg, len: usize, random_pct: u64) -> Vec<u8> {
    let mut v = Vec::with_capacity(len);
    while v.len() < len {
        if rng.chance(random_pct) {
            for _ in 0..rng.below(5000) + 1 {
                v.push(rng.next() as u8);
            }
        } else {
            let line = format!("chr{}\t{}\trs{}\tA\tG\t{}\tPASS\tDP={}\n", rng.below(22) + 1, rng.below(1 << 28), rng.next(), rng.below(99), rng.below(500));
            v.extend_from_slice(line.as_bytes());
        }
    }
    v.truncate(len);
    v
}

fn bgzf_write(path: &str, data: &[u8], pieces: &mut Lcg) {
    let mut w = BgzfWriter::novo(std::fs::File::create(path).unwrap(), 6);
    let mut i = 0;
    while i < data.len() {
        let n = (pieces.below(200_000) as usize + 1).min(data.len() - i);
        w.write_all(&data[i..i + n]).unwrap();
        i += n;
    }
    w.termina().unwrap();
}

/// Walks the blocks of a BGZF file by the size each one declares.
fn block_sizes(raw: &[u8]) -> Vec<usize> {
    let mut sizes = Vec::new();
    let mut i = 0;
    while i < raw.len() {
        assert!(raw[i] == 0x1f && raw[i + 1] == 0x8b && raw[i + 3] == 4, "block at byte {i} is not BGZF");
        assert!(&raw[i + 12..i + 14] == b"BC", "block at byte {i} has no BC field");
        let size = u16::from_le_bytes([raw[i + 16], raw[i + 17]]) as usize + 1;
        sizes.push(size);
        i += size;
    }
    assert!(i == raw.len(), "the block sizes do not add up to the file size");
    sizes
}

fn test_bgzf(t: &mut Tests) {
    let mut rng = Lcg(2026);
    let sizes = [0usize, 1, 100, 0xff00 - 1, 0xff00, 0xff00 + 1, 3 * 0xff00, 1_000_003, 9_500_000];
    for (case, &len) in sizes.iter().enumerate() {
        for random_pct in [0u64, 30, 100] {
            let data = sample_data(&mut rng, len, random_pct);
            let path = t.path(&format!("bgzf_{case}_{random_pct}.gz"));
            bgzf_write(&path, &data, &mut rng);

            // 1. An ordinary gzip reader gets the same bytes back.
            let back = gunzip(&path);
            assert!(back == data, "gzip read-back differs ({len} bytes, {random_pct}% random): {}", first_difference(&back, &data));

            // 2. Well-formed blocks, each within 64 KiB, and the end marker.
            let raw = std::fs::read(&path).unwrap();
            let blocks = block_sizes(&raw);
            assert!(blocks.iter().all(|&s| s <= 65536), "a block is larger than 64 KiB");
            assert!(raw.ends_with(&bgzf::MARCA_DE_FIM), "the end-of-file marker is missing");
            assert!(bgzf::e_bgzf(&path), "e_bgzf does not recognise its own file");

            // 3. Both readers of the module.
            let mut one = Vec::new();
            BgzfReader::novo(std::fs::File::open(&path).unwrap()).read_to_end(&mut one).unwrap();
            assert!(one == data, "BgzfReader differs: {}", first_difference(&one, &data));
            let mut par = Vec::new();
            BgzfParLeitor::novo(std::fs::File::open(&path).unwrap()).read_to_end(&mut par).unwrap();
            assert!(par == data, "BgzfParLeitor differs: {}", first_difference(&par, &data));
        }
        t.ok(&format!("BGZF round trip, {len} bytes (compressible, mixed, random)"));
    }

    // Virtual offsets: note where random stretches start, then seek to them
    // in another order and read the same bytes.
    let data = sample_data(&mut rng, 6_000_000, 20);
    let path = t.path("seek.gz");
    bgzf_write(&path, &data, &mut rng);
    let mut r = BgzfReader::novo(std::fs::File::open(&path).unwrap());
    let mut marks: Vec<(u64, usize)> = Vec::new();
    let mut at = 0usize;
    let mut buf = vec![0u8; 70_000];
    while at < data.len() {
        marks.push((r.pos_virtual(), at));
        let n = (rng.below(70_000) as usize + 1).min(data.len() - at);
        r.read_exact(&mut buf[..n]).unwrap();
        assert!(buf[..n] == data[at..at + n], "sequential read differs at byte {at}");
        at += n;
    }
    let end = r.pos_virtual();
    assert!(marks.len() > 100);
    for _ in 0..2000 {
        let (voff, at) = marks[rng.below(marks.len() as u64) as usize];
        r.busca(voff).unwrap();
        assert!(r.pos_virtual() == voff || voff & 0xffff == 0, "pos_virtual after seek: {} instead of {voff}", r.pos_virtual());
        let n = (rng.below(70_000) as usize + 1).min(data.len() - at);
        r.read_exact(&mut buf[..n]).unwrap();
        assert!(buf[..n] == data[at..at + n], "read after seek to {voff} differs (byte {at} of the data)");
    }
    r.busca(end).unwrap();
    assert!(r.read(&mut buf).unwrap() == 0, "reading at the final offset should give end of file");
    t.ok("BGZF virtual offsets: 2000 random seeks read the right bytes");

    // Lines through BufRead, across block borders.
    let mut lines = 0usize;
    let mut r = BgzfReader::novo(std::fs::File::open(&path).unwrap());
    let mut line = Vec::new();
    let mut total = 0usize;
    loop {
        line.clear();
        let n = r.read_until(b'\n', &mut line).unwrap();
        if n == 0 {
            break;
        }
        assert!(line[..] == data[total..total + n], "line {lines} differs");
        total += n;
        lines += 1;
    }
    assert!(total == data.len() && lines > 1000);
    t.ok("BGZF line reading across block borders");

    // Damage and truncation are reported, not swallowed.
    let mut raw = std::fs::read(&path).unwrap();
    let middle = raw.len() / 2;
    raw[middle] ^= 0x55;
    let damaged = t.path("damaged.gz");
    write_file(&damaged, &raw);
    let mut sink = Vec::new();
    let e = BgzfReader::novo(std::fs::File::open(&damaged).unwrap()).read_to_end(&mut sink);
    assert!(e.is_err(), "a damaged block was read without error");
    let e = BgzfParLeitor::novo(std::fs::File::open(&damaged).unwrap()).read_to_end(&mut sink);
    assert!(e.is_err(), "a damaged block was read without error by the parallel reader");
    raw[middle] ^= 0x55;
    let cut = t.path("cut.gz");
    write_file(&cut, &raw[..raw.len() - 5000]);
    let e = BgzfReader::novo(std::fs::File::open(&cut).unwrap()).read_to_end(&mut sink);
    assert!(e.is_err(), "a truncated file was read without error");
    t.ok("BGZF damaged and truncated files give an error");

    // Plain gzip is refused with an explanation.
    let plain = t.path("plain.gz");
    let mut enc = flate2::write::GzEncoder::new(std::fs::File::create(&plain).unwrap(), flate2::Compression::fast());
    enc.write_all(&data[..100_000]).unwrap();
    enc.finish().unwrap();
    assert!(!bgzf::e_bgzf(&plain));
    let e = BgzfReader::novo(std::fs::File::open(&plain).unwrap()).read_to_end(&mut sink).unwrap_err().to_string();
    assert!(e.contains("plain gzip"), "unexpected message for a plain gzip file: {e}");
    t.ok("plain gzip is told apart from BGZF");
}

fn main() {
    let dir = std::env::temp_dir().join(format!("mmgen_formats_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut t = Tests { n: 0, dir: dir.clone() };

    test_bgzf(&mut t);
    indice_tests::run(&mut t);
    bam_tests::run(&mut t);

    let _ = std::fs::remove_dir_all(&dir);
    println!("{} format tests passed", t.n);
}

mod indice_tests {
    use super::*;
    use externo::Formato;

    /// One record of a test file: the line, and the 0-based half-open span the
    /// reference implementation uses to decide whether a region touches it.
    pub struct Rec {
        pub line: String,
        pub chrom: String,
        pub beg: i64,
        pub end: i64,
    }

    pub fn touches(r: &Rec, chrom: &str, beg: i64, end: i64) -> bool {
        r.chrom == chrom && r.beg < end && r.end.max(r.beg + 1) > beg
    }

    fn make_vcf(rng: &mut Lcg, chroms: &[&str], per_chrom: usize, step: u64, noisy: bool) -> (Vec<String>, Vec<Rec>) {
        let mut header = vec!["##fileformat=VCFv4.2".to_string()];
        for c in chroms {
            header.push(format!("##contig=<ID={c}>"));
        }
        header.push("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string());
        let bases = [b'A', b'C', b'G', b'T'];
        let mut recs = Vec::new();
        for c in chroms {
            let mut pos = 1 + rng.below(step) as i64;
            for _ in 0..per_chrom {
                pos += rng.below(step) as i64; // repeated positions happen
                let ref_len = if rng.chance(5) { rng.below(300) as usize + 2 } else { 1 };
                let reference: String = (0..ref_len).map(|_| bases[rng.below(4) as usize] as char).collect();
                let plain_end = pos - 1 + ref_len as i64;
                let mut end = plain_end;
                let mut info = format!("DP={}", rng.below(900));
                if rng.chance(2) {
                    // A symbolic deletion: INFO/END stretches the record.
                    let e = pos + rng.below(200_000) as i64 + 1;
                    info = format!("SVTYPE=DEL;END={e};{info}");
                    end = e;
                }
                if noisy {
                    info.push_str(&format!(";X={:016x}{:016x}", rng.next(), rng.next()));
                }
                let alt = if end > plain_end { "<DEL>" } else { "T" };
                let line = format!("{c}\t{pos}\trs{}\t{reference}\t{alt}\t{}\tPASS\t{info}", rng.below(1 << 30), rng.below(90));
                recs.push(Rec { line, chrom: c.to_string(), beg: pos - 1, end });
            }
        }
        (header, recs)
    }

    fn make_bed(rng: &mut Lcg, chroms: &[&str], per_chrom: usize) -> (Vec<String>, Vec<Rec>) {
        let header = vec!["track name=test".to_string(), "# a comment".to_string()];
        let mut recs = Vec::new();
        for c in chroms {
            let mut start = 0i64;
            for i in 0..per_chrom {
                start += rng.below(3000) as i64;
                let len = if rng.chance(3) { rng.below(2_000_000) } else { rng.below(1500) } as i64;
                let line = format!("{c}\t{start}\t{}\tfeat{i}\t{}\t+", start + len, rng.below(1000));
                recs.push(Rec { line, chrom: c.to_string(), beg: start, end: start + len });
            }
        }
        (header, recs)
    }

    fn make_gff(rng: &mut Lcg, chroms: &[&str], per_chrom: usize) -> (Vec<String>, Vec<Rec>) {
        let header = vec!["##gff-version 3".to_string()];
        let mut recs = Vec::new();
        for c in chroms {
            let mut start = 1i64;
            for i in 0..per_chrom {
                start += rng.below(5000) as i64;
                let end = start + if rng.chance(10) { rng.below(400_000) } else { rng.below(3000) } as i64;
                let line = format!("{c}\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID=g{i}");
                recs.push(Rec { line, chrom: c.to_string(), beg: start - 1, end });
            }
        }
        (header, recs)
    }

    fn write_bgzf_lines(path: &str, header: &[String], recs: &[Rec]) {
        let mut w = BgzfWriter::novo(std::fs::File::create(path).unwrap(), 6);
        for l in header {
            writeln!(w, "{l}").unwrap();
        }
        for r in recs {
            writeln!(w, "{}", r.line).unwrap();
        }
        w.termina().unwrap();
    }

    pub fn query(path: &str, regions: &[String]) -> Result<(Vec<String>, Vec<String>), String> {
        let (mut header, mut lines) = (Vec::new(), Vec::new());
        let n = indice::consultar(path, regions, &mut |is_header, l| {
            if is_header {
                header.push(l.to_string())
            } else {
                lines.push(l.to_string())
            }
            Ok(())
        })?;
        assert!(n as usize == lines.len(), "consultar returned {n} but delivered {} records", lines.len());
        Ok((header, lines))
    }

    /// Compares `rounds` random region queries with a linear scan.
    pub fn check_regions(what: &str, path: &str, recs: &[Rec], chroms: &[&str], span: i64, rounds: usize, rng: &mut Lcg) {
        for round in 0..rounds {
            let c = chroms[rng.below(chroms.len() as u64) as usize];
            let (text, beg, end) = match round % 4 {
                0 => (c.to_string(), 0, i64::MAX),
                1 => {
                    let p = rng.below(span as u64) as i64 + 1;
                    (format!("{c}:{p}"), p - 1, p)
                }
                _ => {
                    let a = rng.below(span as u64) as i64 + 1;
                    let len = if rng.chance(50) { rng.below(2000) } else { rng.below(span as u64 / 3) } as i64;
                    (format!("{c}:{a}-{}", a + len), a - 1, a + len)
                }
            };
            let want: Vec<String> = recs.iter().filter(|r| touches(r, c, beg, end)).map(|r| r.line.clone()).collect();
            let (_, got) = query(path, &[text.clone()]).unwrap_or_else(|e| panic!("{what}: query {text} failed: {e}"));
            same_lines(&format!("{what}, region {text}"), &got, &want);
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn test_format(
        t: &mut Tests,
        rng: &mut Lcg,
        name: &str,
        formato: Formato,
        header: Vec<String>,
        recs: Vec<Rec>,
        chroms: &[&str],
        rounds: usize,
    ) {
        let ext = match formato {
            Formato::Vcf => "vcf",
            Formato::Bed => "bed",
            _ => "gff",
        };
        let path = t.path(&format!("{}.{ext}.gz", name.replace(' ', "_")));
        write_bgzf_lines(&path, &header, &recs);
        let resumo = indice::indexar(&path, None).unwrap_or_else(|e| panic!("{name}: indexing failed: {e}"));
        assert!(resumo.tipo == "tbi" && resumo.indice == format!("{path}.tbi"));
        assert!(resumo.registros as usize == recs.len(), "{name}: indexed {} of {} records", resumo.registros, recs.len());
        assert!(resumo.sequencias == chroms.len() && resumo.ignoradas == 0);
        let span = recs.iter().map(|r| r.end).max().unwrap_or(1) + 10_000;
        check_regions(name, &path, &recs, chroms, span, rounds, rng);

        // The header comes back, and the index survives a write/read cycle.
        let (got_header, _) = query(&path, &[chroms[0].to_string()]).unwrap();
        same_lines(&format!("{name}: header"), &got_header, &header);
        let (idx, conf, nomes) = indice::le_tbi(&resumo.indice).unwrap();
        assert!(nomes == chroms.iter().map(|c| c.to_string()).collect::<Vec<_>>(), "{name}: sequence names in the index");
        let again = t.path("again.tbi");
        indice::grava_tbi(&idx, &conf, &nomes, &again).unwrap();
        assert!(indice::le_tbi(&again).unwrap() == (idx, conf, nomes), "{name}: index changed after a write/read cycle");
        assert!(bgzf::e_bgzf(&resumo.indice), "{name}: the .tbi is not BGZF");
        t.ok(&format!("tabix index of {name}: {} records, {rounds} region queries equal a linear scan", recs.len()));
    }

    pub fn run(t: &mut Tests) {
        let mut rng = Lcg(777);

        // The line parser, against values worked out by hand.
        let vcf = indice::Tabix::de(Formato::Vcf).unwrap();
        assert!(vcf.intervalo("chr1\t100\t.\tACGT\tA\t.\t.\t.") == Some(("chr1", 99, 103)));
        assert!(vcf.intervalo("chr1\t100\t.\tA\t<DEL>\t.\t.\tSVTYPE=DEL;END=5000") == Some(("chr1", 99, 5000)));
        assert!(vcf.intervalo("chr1\t100\t.\tA\tT\t.\t.\tEND=5000;DP=3") == Some(("chr1", 99, 5000)));
        assert!(vcf.intervalo("chr1\t100\t.\tA\tT\t.\t.\tBEND=5000;DP=3") == Some(("chr1", 99, 100)));
        assert!(vcf.intervalo("chr1\tx\t.\tA\tT").is_none());
        let bed = indice::Tabix::de(Formato::Bed).unwrap();
        assert!(bed.intervalo("chr2\t10\t20\tname") == Some(("chr2", 10, 20)));
        assert!(bed.intervalo("chr2\t10") == Some(("chr2", 10, 11)));
        let gff = indice::Tabix::de(Formato::Gff).unwrap();
        assert!(gff.intervalo("chr3\tsrc\tgene\t11\t20\t.\t+\t.\tID=x") == Some(("chr3", 10, 20)));
        assert!(indice::Tabix::de(Formato::Sam).is_err());
        t.ok("tabix line parser: VCF (REF length, INFO/END), BED, GFF");

        // Bin numbers, against the definition: the smallest bin whose window
        // holds the whole interval.
        for _ in 0..200_000 {
            let beg = rng.below(1 << 29) as i64;
            let bits = rng.below(28) + 1;
            let end = (beg + 1 + rng.below(1 << bits) as i64).min(1 << 29);
            let mut want = 0u32;
            for (level, first) in [(5u32, 4681u32), (4, 585), (3, 73), (2, 9), (1, 1)] {
                let size = 1i64 << (14 + 3 * (5 - level));
                if beg / size == (end - 1) / size {
                    want = first + (beg / size) as u32;
                    break;
                }
            }
            assert!(bam::reg2bin(beg, end) == want, "reg2bin({beg}, {end}) = {} instead of {want}", bam::reg2bin(beg, end));
        }
        assert!(bam::reg2bin(-1, 0) == 4680);
        t.ok("reg2bin equals the definition on 200000 random intervals");

        let chroms = ["chr1", "chr2", "chr10", "chrX"];
        let (h, r) = make_vcf(&mut rng, &chroms, 60_000, 90, false);
        test_format(t, &mut rng, "dense VCF", Formato::Vcf, h, r, &chroms, 240);
        let (h, r) = make_vcf(&mut rng, &chroms[..2], 40_000, 4000, true);
        test_format(t, &mut rng, "sparse VCF that compresses badly", Formato::Vcf, h, r, &chroms[..2], 240);
        let (h, r) = make_vcf(&mut rng, &["1"], 300, 50, false);
        test_format(t, &mut rng, "tiny VCF", Formato::Vcf, h, r, &["1"], 120);
        let (h, r) = make_bed(&mut rng, &chroms, 30_000);
        test_format(t, &mut rng, "BED", Formato::Bed, h, r, &chroms, 240);
        let (h, r) = make_gff(&mut rng, &chroms[..3], 30_000);
        test_format(t, &mut rng, "GFF", Formato::Gff, h, r, &chroms[..3], 240);

        // Several regions at once: each record once, in file order.
        let (h, r) = make_vcf(&mut rng, &chroms, 5000, 100, false);
        let path = t.path("multi.vcf.gz");
        write_bgzf_lines(&path, &h, &r);
        indice::indexar(&path, None).unwrap();
        let regions = ["chrX:1000-9000", "chr2:5,000-20,000", "chrX:8000-12000", "chr2:19000-", "chr1:777"].map(String::from);
        let spans = [("chr1", 776, 777), ("chr2", 4999, i64::MAX), ("chrX", 999, 12000)];
        let want: Vec<String> = r.iter().filter(|x| spans.iter().any(|s| touches(x, s.0, s.1, s.2))).map(|x| x.line.clone()).collect();
        let (_, got) = query(&path, &regions).unwrap();
        same_lines("several regions", &got, &want);
        t.ok("several overlapping regions: each record once, in file order");

        // Errors that must explain themselves.
        let e = query(&path, &["chr7:1-100".to_string()]).unwrap_err();
        assert!(e.contains("chr7") && e.contains("chr1"), "unknown sequence: {e}");
        let e = query(&path, &["chr1:abc".to_string()]).unwrap_err();
        assert!(e.contains("not valid"), "bad region: {e}");
        let e = query(&path, &["chr1:500-100".to_string()]).unwrap_err();
        assert!(e.contains("not valid"), "inverted region: {e}");
        std::fs::remove_file(format!("{path}.tbi")).unwrap();
        let e = query(&path, &["chr1".to_string()]).unwrap_err();
        assert!(e.contains("index_file"), "missing index: {e}");
        t.ok("region errors: unknown sequence, malformed region, missing index");

        let mut shuffled: Vec<Rec> = r;
        shuffled.swap(10, 4000);
        let path = t.path("unsorted.vcf.gz");
        write_bgzf_lines(&path, &h, &shuffled);
        let e = indice::indexar(&path, None).unwrap_err();
        assert!(e.contains("not sorted") && e.contains("sort_genomic") && e.contains("line "), "unsorted: {e}");
        let (h2, r2) = make_vcf(&mut rng, &["chr1", "chr2", "chr1"], 50, 100, false);
        write_bgzf_lines(&path, &h2, &r2);
        let e = indice::indexar(&path, None).unwrap_err();
        assert!(e.contains("not sorted"), "split sequence: {e}");
        let plain = t.path("plain.vcf.gz");
        let mut enc = flate2::write::GzEncoder::new(std::fs::File::create(&plain).unwrap(), flate2::Compression::fast());
        enc.write_all(b"#CHROM\tPOS\nchr1\t1\t.\tA\tT\t.\t.\t.\n").unwrap();
        enc.finish().unwrap();
        let e = indice::indexar(&plain, None).unwrap_err();
        assert!(e.contains("BGZF"), "plain gzip: {e}");
        let big = t.path("big.vcf.gz");
        let far = Rec { line: "chr1\t600000000\t.\tA\tT\t.\t.\t.".into(), chrom: "chr1".into(), beg: 0, end: 0 };
        write_bgzf_lines(&big, &[], &[far]);
        let e = indice::indexar(&big, None).unwrap_err();
        assert!(e.contains("512 Mbp"), "position beyond the index range: {e}");
        t.ok("indexing errors: unsorted, split sequence, plain gzip, beyond 512 Mbp");

        // A header-only file still gets a valid index.
        let empty = t.path("empty.vcf.gz");
        write_bgzf_lines(&empty, &h, &[]);
        let resumo = indice::indexar(&empty, None).unwrap();
        assert!(resumo.registros == 0 && resumo.sequencias == 0);
        let e = query(&empty, &["chr1".to_string()]).unwrap_err();
        assert!(e.contains("not in"), "query on an empty file: {e}");
        t.ok("header-only file: valid empty index");
    }
}

mod bam_tests {
    use super::indice_tests::{check_regions, query, Rec};
    use super::*;
    use bam::{BamEscritor, BamLeitor};

    const REFS: [(&str, u32); 4] = [("chr1", 3_000_000), ("chr2", 2_000_000), ("chrM", 16_569), ("HLA-A*01:01", 3_503)];

    fn header() -> Vec<String> {
        let mut h = vec!["@HD\tVN:1.6\tSO:coordinate".to_string()];
        for (n, l) in REFS {
            h.push(format!("@SQ\tSN:{n}\tLN:{l}"));
        }
        h.push("@RG\tID:g1\tSM:sample one\tPL:ILLUMINA".to_string());
        h.push("@PG\tID:test\tPN:formats\tCL:made up".to_string());
        h.push("@CO\ta free comment, with\ttabs".to_string());
        h
    }

    fn sam_to_bam(path: &str, lines: &[String]) -> std::io::Result<()> {
        let mut w = BamEscritor::novo(std::fs::File::create(path)?);
        for l in lines {
            w.write_all(l.as_bytes())?;
            w.write_all(b"\n")?;
        }
        w.termina().map(|_| ())
    }

    fn bam_to_sam(path: &str) -> std::io::Result<Vec<String>> {
        let mut text = String::new();
        BamLeitor::novo(std::fs::File::open(path)?)?.read_to_string(&mut text)?;
        Ok(text.lines().map(str::to_string).collect())
    }

    fn pick<'a>(rng: &mut Lcg, options: &[&'a str]) -> &'a str {
        options[rng.below(options.len() as u64) as usize]
    }

    /// Optional fields in the exact text form samtools prints, so that a
    /// round trip must return the line unchanged.
    fn random_tags(rng: &mut Lcg) -> String {
        let mut s = String::new();
        let ints = ["0", "1", "-1", "127", "128", "-128", "-129", "255", "256", "32767", "32768", "-32768", "-32769", "65535", "65536", "2147483647", "-2147483648", "4294967295"];
        let floats = ["0", "1", "-2.5", "3.14159", "1e+10", "1e-05", "0.0001", "123456", "1.23457e+08", "-0.333333", "inf"];
        if rng.chance(70) {
            s.push_str(&format!("\tNM:i:{}", pick(rng, &ints)));
        }
        if rng.chance(40) {
            s.push_str(&format!("\tXA:A:{}", pick(rng, &["a", "Z", "!", "~", "7"])));
        }
        if rng.chance(40) {
            s.push_str(&format!("\tXF:f:{}", pick(rng, &floats)));
        }
        if rng.chance(50) {
            s.push_str(&format!("\tRG:Z:{}", pick(rng, &["g1", "", "with space", "a:b:c", "x=y;z"])));
        }
        if rng.chance(20) {
            s.push_str(&format!("\tXH:H:{}", pick(rng, &["1AE301", "", "00FF"])));
        }
        if rng.chance(30) {
            let (kind, values): (&str, &[&str]) = match rng.below(7) {
                0 => ("c", &["-128", "0", "127"]),
                1 => ("C", &["0", "255", "7"]),
                2 => ("s", &["-32768", "32767", "5"]),
                3 => ("S", &["0", "65535"]),
                4 => ("i", &["-2147483648", "2147483647", "0"]),
                5 => ("I", &["4294967295", "0", "12"]),
                _ => ("f", &["1.5", "-0.25", "1e+20"]),
            };
            s.push_str(&format!("\tXB:B:{kind}"));
            for _ in 0..rng.below(6) {
                s.push_str(&format!(",{}", pick(rng, values)));
            }
        }
        s
    }

    fn random_cigar(rng: &mut Lcg, read_len: usize) -> (String, i64) {
        // Splits the read among M, I, S and =/X; D and N only advance the reference.
        let mut left = read_len;
        let (mut s, mut ref_len) = (String::new(), 0i64);
        if rng.chance(20) && left > 4 {
            let n = rng.below(4) as usize + 1;
            s.push_str(&format!("{n}S"));
            left -= n;
        }
        while left > 0 {
            let n = (rng.below(40) as usize + 1).min(left);
            let op = pick(rng, &["M", "M", "M", "=", "X", "I"]);
            s.push_str(&format!("{n}{op}"));
            left -= n;
            if op != "I" {
                ref_len += n as i64;
                if left > 0 && rng.chance(25) {
                    let d = if rng.chance(10) { rng.below(30_000) } else { rng.below(20) } as i64 + 1;
                    s.push_str(&format!("{d}{}", if d > 20 { "N" } else { "D" }));
                    ref_len += d;
                }
            }
        }
        (s, ref_len)
    }

    /// A coordinate-sorted set of alignments: mapped reads on each reference
    /// (some flagged unmapped but placed), then reads with no position.
    fn random_alignments(rng: &mut Lcg, per_ref: usize, unplaced: usize) -> Vec<Rec> {
        let iupac = b"ACGTNACGTACGTRYKMSWBDHV=";
        let mut recs = Vec::new();
        let mut n = 0usize;
        let read = |rng: &mut Lcg, len: usize| -> (String, String) {
            let seq: String = (0..len).map(|_| iupac[rng.below(iupac.len() as u64) as usize] as char).collect();
            let qual: String = if rng.chance(10) { "*".into() } else { (0..len).map(|_| (33 + rng.below(60) as u8) as char).collect() };
            (seq, qual)
        };
        for (ri, (name, len)) in REFS.iter().enumerate() {
            let mut pos = 1i64;
            for _ in 0..per_ref {
                pos += rng.below(u64::from(*len) / per_ref as u64 + 1) as i64;
                if pos >= i64::from(*len) {
                    break;
                }
                n += 1;
                let read_len = rng.below(150) as usize + 1;
                let (seq, qual) = read(rng, read_len);
                let unmapped_placed = rng.chance(4);
                let (cigar, ref_len) = if unmapped_placed { ("*".to_string(), 0) } else { random_cigar(rng, read_len) };
                let flag = if unmapped_placed { 4 + 1 + 64 } else { [0, 16, 99, 147, 83, 163, 2048, 256][rng.below(8) as usize] };
                let (rnext, pnext) = match rng.below(4) {
                    0 => ("*".to_string(), 0),
                    1 => ("=".to_string(), pos + rng.below(500) as i64),
                    _ => (REFS[(ri + 1) % REFS.len()].0.to_string(), rng.below(10_000) as i64 + 1),
                };
                let line = format!(
                    "read{n}\t{flag}\t{name}\t{pos}\t{}\t{cigar}\t{rnext}\t{pnext}\t{}\t{seq}\t{qual}{}",
                    rng.below(61),
                    rng.below(2000) as i64 - 1000,
                    random_tags(rng)
                );
                recs.push(Rec { line, chrom: name.to_string(), beg: pos - 1, end: pos - 1 + ref_len.max(1) });
            }
        }
        for _ in 0..unplaced {
            n += 1;
            let read_len = rng.below(100) as usize + 1;
            let (seq, qual) = read(rng, read_len);
            let line = format!("read{n}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}{}", random_tags(rng));
            recs.push(Rec { line, chrom: "*".to_string(), beg: -1, end: 0 });
        }
        recs
    }

    pub fn run(t: &mut Tests) {
        let mut rng = Lcg(4242);

        // 1. One record, byte by byte, against the layout in the SAM/BAM
        //    specification worked out by hand.
        let path = t.path("one.bam");
        let lines = vec![
            "@SQ\tSN:chr1\tLN:1000".to_string(),
            "r1\t99\tchr1\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGATACTG\t*\tNM:i:1".to_string(),
        ];
        sam_to_bam(&path, &lines).unwrap();
        let raw = gunzip(&path);
        let mut want: Vec<u8> = Vec::new();
        want.extend_from_slice(b"BAM\x01");
        want.extend_from_slice(&20i32.to_le_bytes());
        want.extend_from_slice(b"@SQ\tSN:chr1\tLN:1000\n");
        want.extend_from_slice(&1i32.to_le_bytes());
        want.extend_from_slice(&5i32.to_le_bytes());
        want.extend_from_slice(b"chr1\0");
        want.extend_from_slice(&1000i32.to_le_bytes());
        let mut rec: Vec<u8> = Vec::new();
        rec.extend_from_slice(&0i32.to_le_bytes()); // refID
        rec.extend_from_slice(&6i32.to_le_bytes()); // pos, 0-based
        rec.push(3); // l_read_name
        rec.push(30); // mapq
        rec.extend_from_slice(&4681u16.to_le_bytes()); // bin
        rec.extend_from_slice(&5u16.to_le_bytes()); // n_cigar_op
        rec.extend_from_slice(&99u16.to_le_bytes()); // flag
        rec.extend_from_slice(&17i32.to_le_bytes()); // l_seq
        rec.extend_from_slice(&0i32.to_le_bytes()); // next refID
        rec.extend_from_slice(&36i32.to_le_bytes()); // next pos
        rec.extend_from_slice(&39i32.to_le_bytes()); // tlen
        rec.extend_from_slice(b"r1\0");
        for (n, op) in [(8u32, 0u32), (2, 1), (4, 0), (1, 2), (3, 0)] {
            rec.extend_from_slice(&(n << 4 | op).to_le_bytes());
        }
        rec.extend_from_slice(&[0x88, 0x14, 0x18, 0x11, 0x14, 0x41, 0x81, 0x28, 0x40]); // TTAGATAAAGGATACTG
        rec.extend_from_slice(&[0xff; 17]);
        rec.extend_from_slice(b"NMC\x01");
        want.extend_from_slice(&(rec.len() as i32).to_le_bytes());
        want.extend_from_slice(&rec);
        assert!(raw == want, "BAM bytes differ from the hand-built record: {}", first_difference(&raw, &want));
        assert!(std::fs::read(&path).unwrap().ends_with(&bgzf::MARCA_DE_FIM), "BAM without the BGZF end marker");
        same_lines("one record back to SAM", &bam_to_sam(&path).unwrap(), &lines);
        t.ok("BAM encoding of one record equals the bytes worked out by hand");

        // 2. Round trips of random alignments with every kind of field.
        for (case, (per_ref, unplaced)) in [(2000usize, 300usize), (40_000, 5000), (1, 0), (0, 50)].into_iter().enumerate() {
            let recs = random_alignments(&mut rng, per_ref, unplaced);
            let mut lines = header();
            lines.extend(recs.iter().map(|r| r.line.clone()));
            let path = t.path(&format!("rt{case}.bam"));
            sam_to_bam(&path, &lines).unwrap_or_else(|e| panic!("writing BAM failed: {e}"));
            let back = bam_to_sam(&path).unwrap_or_else(|e| panic!("reading BAM failed: {e}"));
            same_lines("SAM -> BAM -> SAM", &back, &lines);
            t.ok(&format!("SAM -> BAM -> SAM returns the same {} lines", lines.len()));

            // 3. The BAI index: region queries against a linear scan.
            if per_ref >= 2000 {
                let resumo = indice::indexar(&path, None).unwrap_or_else(|e| panic!("indexing the BAM failed: {e}"));
                assert!(resumo.tipo == "bai" && resumo.registros as usize == recs.len());
                let names: Vec<&str> = REFS.iter().map(|r| r.0).collect();
                check_regions("BAM", &path, &recs, &names, 3_000_000, 300, &mut rng);
                let idx = indice::le_bai(&resumo.indice).unwrap();
                assert!(idx.sem_coordenada as usize == unplaced, "{} unplaced reads in the index, expected {unplaced}", idx.sem_coordenada);
                for (i, (name, _)) in REFS.iter().enumerate() {
                    let of_ref: Vec<&Rec> = recs.iter().filter(|r| r.chrom == *name).collect();
                    let unmapped = of_ref.iter().filter(|r| r.line.split('\t').nth(1) == Some("69")).count() as u64;
                    let meta = idx.refs[i].meta.unwrap_or([0; 4]);
                    assert!(meta[2] == of_ref.len() as u64 - unmapped && meta[3] == unmapped, "{name}: mapped/unmapped counts in the index: {meta:?}");
                }
                let again = t.path("again.bai");
                indice::grava_bai(&idx, &again).unwrap();
                assert!(indice::le_bai(&again).unwrap() == idx, "BAI changed after a write/read cycle");
                let (h, _) = query(&path, &["chrM".to_string()]).unwrap();
                same_lines("BAM header from a query", &h, &header());
                // A sequence name with ':' in it is still found.
                let want: Vec<String> = recs.iter().filter(|r| r.chrom == "HLA-A*01:01").map(|r| r.line.clone()).collect();
                let (_, got) = query(&path, &["HLA-A*01:01".to_string()]).unwrap();
                same_lines("region named HLA-A*01:01", &got, &want);
                t.ok(&format!("BAI index: {} alignments, 300 region queries equal a linear scan", recs.len()));
            }
        }

        // 4. Text that is not in samtools' canonical form comes back canonical.
        let path = t.path("canon.bam");
        let given = vec![
            "@SQ\tSN:chr1\tLN:1000".to_string(),
            "q\t0\tchr1\t5\t9\t4M\tchr1\t9\t0\tacgt\tIIII\tXF:f:1e10\tXG:f:0.00001\tXI:i:+7\tXJ:f:123456789".to_string(),
        ];
        let canon = vec![
            "@SQ\tSN:chr1\tLN:1000".to_string(),
            "q\t0\tchr1\t5\t9\t4M\t=\t9\t0\tACGT\tIIII\tXF:f:1e+10\tXG:f:1e-05\tXI:i:7\tXJ:f:1.23457e+08".to_string(),
        ];
        sam_to_bam(&path, &given).unwrap();
        same_lines("canonical form", &bam_to_sam(&path).unwrap(), &canon);
        t.ok("non-canonical text (lower case, RNEXT by name, float forms) comes back as samtools prints it");

        // 5. A CIGAR with more than 65535 operations goes through the CG tag.
        let ops = 70_000usize;
        let cigar: String = (0..ops).map(|i| if i % 2 == 0 { "1M" } else { "1I" }).collect();
        let seq = "A".repeat(ops);
        let long = vec![
            "@SQ\tSN:chr1\tLN:3000000".to_string(),
            format!("long\t0\tchr1\t100\t60\t{cigar}\t*\t0\t0\t{seq}\t*\tNM:i:3"),
        ];
        let path = t.path("long.bam");
        sam_to_bam(&path, &long).unwrap();
        let raw = gunzip(&path);
        let rec = &raw[4 + 4 + 23 + 4 + 4 + 5 + 4 + 4..];
        assert!(u16::from_le_bytes([rec[12], rec[13]]) == 2, "a long CIGAR must leave 2 placeholder operations in the fixed field");
        assert!(raw.windows(4).any(|w| w == b"CGBI"), "the CG tag is missing");
        same_lines("long CIGAR", &bam_to_sam(&path).unwrap(), &long);
        let info = bam::info(rec).unwrap();
        assert!(info.pos == 99 && info.fim == 99 + 35_000, "reference span of the long CIGAR: {}..{}", info.pos, info.fim);
        t.ok("CIGAR with 70000 operations: CG tag written, text unchanged, reference span right");

        // 6. A header without @SQ in the text gets them from the binary list.
        let path = t.path("nosq.bam");
        let mut w = BgzfWriter::novo(std::fs::File::create(&path).unwrap(), 6);
        let cab = bam::Cabecalho { texto: "@HD\tVN:1.6\n@CO\tno sequences here\n".into(), refs: vec![("chrZ".into(), 77), ("chrW".into(), 5)] };
        bam::grava_cabecalho(&mut w, &cab).unwrap();
        w.termina().unwrap();
        let want = ["@HD\tVN:1.6", "@SQ\tSN:chrZ\tLN:77", "@SQ\tSN:chrW\tLN:5", "@CO\tno sequences here"].map(String::from);
        same_lines("header without @SQ", &bam_to_sam(&path).unwrap(), &want);
        // A header-only SAM becomes a header-only BAM.
        let path = t.path("header_only.bam");
        sam_to_bam(&path, &header()).unwrap();
        same_lines("header-only BAM", &bam_to_sam(&path).unwrap(), &header());
        t.ok("headers: @SQ rebuilt from the binary list; header-only file");

        // 7. Errors name the line and the reason.
        let path = t.path("bad.bam");
        let bad = |line: &str| -> String {
            let lines = vec!["@SQ\tSN:chr1\tLN:1000".to_string(), "ok\t0\tchr1\t1\t0\t1M\t*\t0\t0\tA\tI".to_string(), line.to_string()];
            sam_to_bam(&path, &lines).expect_err("an invalid line was accepted").to_string()
        };
        let e = bad("x\t0\tchr9\t1\t0\t1M\t*\t0\t0\tA\tI");
        assert!(e.contains("line 3") && e.contains("chr9") && e.contains("@SQ"), "undeclared sequence: {e}");
        let e = bad("x\t0\tchr1\t1\t0\t1Q\t*\t0\t0\tA\tI");
        assert!(e.contains("line 3") && e.contains("CIGAR"), "bad CIGAR: {e}");
        let e = bad("x\t0\tchr1\t1\t0\t1M\t*\t0\t0\tAC\tI");
        assert!(e.contains("line 3") && e.contains("same length"), "SEQ and QUAL lengths: {e}");
        let e = bad("x\t0\tchr1\t1\t0\t1M\t*\t0\t0\tA");
        assert!(e.contains("line 3") && e.contains("QUAL"), "missing column: {e}");
        let e = bad("x\t0\tchr1\t1\t0\t1M\t*\t0\t0\tA\tI\tNM:i:abc");
        assert!(e.contains("line 3") && e.contains("NM:i:abc"), "bad tag: {e}");
        let e = bad("x\t0\tchr1\t1\t0\t1M\t*\t0\t0\tA\tI\tNM:i:99999999999");
        assert!(e.contains("32 bits"), "integer too large: {e}");
        let e = bad("@CO\ta late header line");
        assert!(e.contains("line 3") && e.contains("header"), "late header line: {e}");
        let no_ln = vec!["@SQ\tSN:chr1".to_string(), "x\t0\tchr1\t1\t0\t1M\t*\t0\t0\tA\tI".to_string()];
        let e = sam_to_bam(&path, &no_ln).unwrap_err().to_string();
        assert!(e.contains("LN"), "@SQ without LN: {e}");
        t.ok("writing BAM: errors name the line and the reason");

        let good = t.path("rt0.bam");
        let raw = std::fs::read(&good).unwrap();
        let cut = t.path("cut.bam");
        write_file(&cut, &raw[..raw.len() * 2 / 3]);
        assert!(bam_to_sam(&cut).is_err(), "a truncated BAM was read without error");
        let not_bam = t.path("not.bam");
        let mut w = BgzfWriter::novo(std::fs::File::create(&not_bam).unwrap(), 6);
        w.write_all(b"@HD\tVN:1.6\nthis is SAM text inside BGZF\n").unwrap();
        w.termina().unwrap();
        let e = bam_to_sam(&not_bam).unwrap_err().to_string();
        assert!(e.contains("not a BAM"), "not a BAM: {e}");
        t.ok("reading BAM: truncated file and wrong signature are errors");

        // 8. An unsorted BAM is refused by the indexer, with the way out.
        let mut recs = random_alignments(&mut rng, 500, 0);
        recs.swap(3, 400);
        let mut lines = header();
        lines.extend(recs.iter().map(|r| r.line.clone()));
        let path = t.path("unsorted.bam");
        sam_to_bam(&path, &lines).unwrap();
        let e = indice::indexar(&path, None).unwrap_err();
        assert!(e.contains("not sorted") && e.contains("sort_genomic") && e.contains("alignment "), "unsorted BAM: {e}");
        t.ok("indexing an unsorted BAM is refused with the way out");
    }
}
