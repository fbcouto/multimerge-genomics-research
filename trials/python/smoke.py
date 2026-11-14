# Smoke test of the installed extension: no external data, no numpy, any
# system, a few seconds. The CI runs it on every wheel, and it is the quick way
# to check an installation:
#     python trials/python/smoke.py
import gzip
import os
import random
import struct
import sys
import tempfile

import multimerge_genomics as mmg

failures = []
done = 0


def check(name, cond, detail=""):
    global done
    done += 1
    if not cond:
        failures.append(name)
    print(("  ok    " if cond else "  FAIL  ") + name + ("" if cond else f"  {detail}"))


info = mmg.engine_info()
print("engine_info:", info)
check("engine_info has the expected keys", all(k in info for k in ("version", "gpu_compiled", "gpu", "gpu_memory_mb", "threads")))

random.seed(42)
T = tempfile.mkdtemp(prefix="mmg_smoke_")
p = lambda name: os.path.join(T, name)


def write(name, lines, gz=False):
    opener = gzip.open if gz else open
    with opener(p(name), "wt", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return p(name)


def read(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        return [l for l in f.read().split("\n") if l]


# ---------------------------------------------------------------- VCF
header = ["##fileformat=VCFv4.2", "##contig=<ID=1,length=248956422>", "##contig=<ID=2,length=242193529>",
          "##contig=<ID=10,length=133797422>", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
records = [f"{c}\t{pos}\trs{i}\tA\tG\t{random.randrange(1, 90)}\tPASS\t." for i, (c, pos) in
           enumerate((random.choice(["1", "2", "10"]), random.randrange(1, 5000)) for _ in range(20000))]
a = write("a.vcf", header + records[:12000])
b = write("b.vcf.gz", header + records[12000:], gz=True)

# Sort two files (text + gz) with a minimal budget: forces runs on disk.
n, runs, skipped = mmg.sort_genomic([a, b], p("sorted.vcf"), chrom_order="natural", memory_mb=1)
ref = sorted(records, key=lambda l: (int(l.split("\t")[0]), int(l.split("\t")[1])))
out = [l for l in read(p("sorted.vcf")) if not l.startswith("#")]
check("sort_genomic: identical to Python's stable sort", out == ref and n == len(records), f"n={n}")
check("sort_genomic: used runs on disk", runs > 1, f"runs={runs}")
check("sort_genomic: header preserved", [l for l in read(p("sorted.vcf")) if l.startswith("#")] == header)

# The function behind screen 8 of the interface: merges several files into one output.
n2 = mmg.multimerge_vcf([a, b], "POS", p("merged.vcf"))
ref_lex = sorted(records, key=lambda l: (l.split("\t")[0], int(l.split("\t")[1])))
check("multimerge_vcf: merges BOTH files", n2 == len(records) and [l for l in read(p("merged.vcf")) if not l.startswith("#")] == ref_lex, f"n={n2}")

# Validator: defects planted on known lines.
bad = header + records[:50]
bad[10] = bad[10].replace("\t", " ", 1)
fields = bad[20].split("\t"); fields[1] = "12k"; bad[20] = "\t".join(fields)
r = mmg.validate_file(write("bad.vcf", bad))
codes = {q["code"]: q["line"] for q in r["problems"]}
check("validate_file: finds the space-separated column on line 11", codes.get("vcf_columns") == 11, str(codes))
check("validate_file: finds the invalid POS on line 21", codes.get("vcf_pos") == 21, str(codes))
check("validate_file: the bad file is not valid; the good one is", r["valid"] is False and mmg.validate_file(p("sorted.vcf"))["valid"] is True)
check("validate_file: sorted / unsorted", mmg.validate_file(p("sorted.vcf"))["sorted"] is True and mmg.validate_file(a)["sorted"] is False)
check("validate_file: the report is in English", "is not an integer" in r["text"] and "fix:" in r["text"], r["text"][:300])

# Chromosome names: there and back.
rn = mmg.normalize_chromosomes(a, p("a_chr.vcf"), style="chr")
mmg.normalize_chromosomes(p("a_chr.vcf"), p("a_back.vcf"), style="no_chr")
check("normalize_chromosomes: round trip", read(p("a_back.vcf")) == read(a) and rn["changed"] == 12000 and read(p("a_chr.vcf"))[-1].startswith("chr"))

# Reference genome.
header37 = [l.replace("248956422", "249250621") for l in header]
c37 = write("c37.vcf", header37 + records[:10])
check("file_reference: recognizes GRCh38", "GRCh38" in (mmg.file_reference(a)["assembly"] or ""))
check("compare_references: GRCh38 x GRCh37 are incompatible", mmg.compare_references(a, c37)["compatible"] is False)

# ---------------------------------------------------------------- BED
bed_a = write("a.bed", ["chr1\t10\t20\tx", "chr1\t15\t30\ty", "chr1\t30\t40\tz", "chr2\t5\t9\tw", "chr1\t100\t110\tk"])
bed_b = write("b.bed", ["chr1\t18\t35\tq", "chr2\t100\t200\tr"])
n_read, n_written = mmg.bed_merge([bed_a], p("merged.bed"))
check("bed_merge", (n_read, n_written) == (5, 3) and read(p("merged.bed")) == ["chr1\t10\t40", "chr1\t100\t110", "chr2\t5\t9"], str(read(p("merged.bed"))))
nu = mmg.bed_intersect([bed_a], [bed_b], p("u.bed"), mode="unique")
check("bed_intersect (unique)", nu == 3 and sorted(read(p("u.bed"))) == ["chr1\t10\t20\tx", "chr1\t15\t30\ty", "chr1\t30\t40\tz"], str(read(p("u.bed"))))
try:
    mmg.bed_intersect([bed_a], [write("nochr.bed", ["1\t1\t2\tx"])], p("x.bed"))
    check("bed_intersect: refuses chr1 against 1", False, "no error was raised")
except RuntimeError as e:
    check("bed_intersect: refuses chr1 against 1", "in common" in str(e), str(e)[:80])
mmg.bed_subtract([bed_a], [bed_b], p("sub.bed"))
check("bed_subtract", read(p("sub.bed")) == ["chr1\t10\t18\tx", "chr1\t15\t18\ty", "chr1\t35\t40\tz", "chr1\t100\t110\tk", "chr2\t5\t9\tw"], str(read(p("sub.bed"))))
mmg.bed_closest([write("one.bed", ["chr1\t50\t60\tm"])], [bed_b], p("closest.bed"))
check("bed_closest: distance of separated intervals", read(p("closest.bed")) == ["chr1\t50\t60\tm\tchr1\t18\t35\tq\t16"], str(read(p("closest.bed"))))
mmg.bed_coverage([write("one2.bed", ["chr1\t10\t30\tm"])], [bed_b], p("cov.bed"))
check("bed_coverage", read(p("cov.bed")) == ["chr1\t10\t30\tm\t1\t12\t20\t0.6000000"], str(read(p("cov.bed"))))
genome = write("genome.txt", ["chr1\t120", "chr2\t50"])
mmg.bed_complement([bed_a], genome, p("comp.bed"))
check("bed_complement", read(p("comp.bed")) == ["chr1\t0\t10", "chr1\t40\t100", "chr1\t110\t120", "chr2\t0\t5", "chr2\t9\t50"], str(read(p("comp.bed"))))
mmg.bed_slop([bed_a], genome, p("slop.bed"), both=15)
check("bed_slop: clamped at 0 and at the sequence end", read(p("slop.bed"))[0] == "chr1\t0\t35\tx" and read(p("slop.bed"))[4] == "chr1\t85\t120\tk", str(read(p("slop.bed"))))

# ---------------------------------------------------------------- VCF tools
st = mmg.vcf_stats(p("sorted.vcf"))
check("vcf_stats: every A>G record is a SNP and a transition", st["records"] == st["snps"] == st["transitions"] == len(records) and st["transversions"] == 0, str(st)[:120])
fr = mmg.vcf_filter(p("sorted.vcf"), p("filtered.vcf"), expression="QUAL >= 50 && CHROM == \"2\"")
want = [l for l in ref if int(l.split("\t")[5]) >= 50 and l.startswith("2\t")]
check("vcf_filter: expression", fr["kept"] == len(want) and [l for l in read(p("filtered.vcf")) if not l.startswith("#")] == want, str(fr))
try:
    mmg.vcf_filter(p("sorted.vcf"), p("x.vcf"), expression="QUAL >> 3")
    check("vcf_filter: a mistake in the expression is reported", False, "no error")
except RuntimeError as e:
    check("vcf_filter: a mistake in the expression is reported with its position", "position" in str(e), str(e)[:120])
multi = write("multi.vcf", ["##fileformat=VCFv4.2", '##INFO=<ID=AF,Number=A,Type=Float,Description="f">', "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                            "1\t100\t.\tCTT\tCT,CTTT\t50\tPASS\tAF=0.1,0.2"])
nr = mmg.vcf_normalize(multi, p("norm.vcf"))
check("vcf_normalize: splits and trims", [l for l in read(p("norm.vcf")) if not l.startswith("#")] ==
      ["1\t100\t.\tCT\tC\t50\tPASS\tAF=0.1", "1\t100\t.\tC\tCT\t50\tPASS\tAF=0.2"] and nr["split"] == 1, str(read(p("norm.vcf"))[-2:]))
# With the reference genome: a T deleted or inserted in TTTTT (96-100), written
# at its right end, moves to the left end of the run, after the A at 95.
genome = write("genome.fa", [">1", "ACGTACGTAC" * 9 + "GGGCATTTTTGCA" + "ACGT" * 10])
run_t = write("run_t.vcf", read(multi)[:-1] + ["1\t99\t.\tTT\tT,TTT\t50\tPASS\tAF=0.1,0.2"])
nr = mmg.vcf_normalize(run_t, p("norm_ref.vcf"), reference=genome)
check("vcf_normalize(reference=...): left-aligns, as bcftools norm -f", [l for l in read(p("norm_ref.vcf")) if not l.startswith("#")] ==
      ["1\t95\t.\tAT\tA\t50\tPASS\tAF=0.1", "1\t95\t.\tA\tAT\t50\tPASS\tAF=0.2"] and nr["left_aligned"] == 2 and nr["ref_mismatch"] == 0,
      str(read(p("norm_ref.vcf"))[-2:]) + str(nr))
# TYPE as in bcftools: == every ALT allele, ~ at least one.
mixed = write("mixed.vcf", header[:1] + header[-1:] + ["1\t10\ta\tA\tG\t50\tPASS\t.", "1\t20\tb\tA\tG,AT\t50\tPASS\t.", "1\t30\tc\tA\tAT\t50\tPASS\t."])
ids = lambda expr: [l.split("\t")[2] for l in read(p("t.vcf")) if not l.startswith("#")] if mmg.vcf_filter(mixed, p("t.vcf"), expression=expr) else None
every, some, none = ids('TYPE=="snp"'), ids('TYPE~"snp"'), ids('TYPE!~"snp"')
check("vcf_filter: TYPE==\"snp\" is every allele, TYPE~\"snp\" at least one", every == ["a"] and some == ["a", "b"] and none == ["c"],
      f"{every} {some} {none}")

# ---------------------------------------------------------------- containers: BGZF, index, BAM, Parquet
gz = p("sorted.vcf.gz")
mmg.sort_genomic([a, b], gz, chrom_order="natural", index=True)
check("sort_genomic(index=True): the .gz is ordinary gzip and has a .tbi", read(gz)[len(header):] == ref and os.path.exists(gz + ".tbi"))
got = mmg.query_region(gz, "2:1000-2000")
check("query_region", got == [l for l in ref if l.startswith("2\t") and 1000 <= int(l.split("\t")[1]) <= 2000], f"{len(got)} records")
sam = write("x.sam", ["@HD\tVN:1.6", "@SQ\tSN:chr1\tLN:5000", "r1\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tNM:i:0",
                      "r2\t16\tchr1\t50\t60\t2M1I2M\t*\t0\t0\tACGTT\t*\tXS:f:1.5"])
mmg.sort_genomic([sam], p("x.bam"), index=True)
mmg.convert_file(p("x.bam"), p("back.sam"))
back = read(p("back.sam"))
check("SAM -> sorted BAM -> SAM", back[2:] == [read(sam)[3], read(sam)[2]] and "SO:coordinate" in back[0], str(back))
check("query_region on a BAM", [l.split("\t")[0] for l in mmg.query_region(p("x.bam"), "chr1:100-100")] == ["r1"])
mmg.convert_file(p("sorted.vcf"), p("sorted.vcf.parquet"))
mmg.convert_file(p("sorted.vcf.parquet"), p("from_parquet.vcf"))
check("VCF -> Parquet -> VCF", read(p("from_parquet.vcf")) == read(p("sorted.vcf")))
# QUAL as a number (the default) or as the text it was: 12.50 comes back as 12.5, or as 12.50.
quals = write("quals.vcf", header[:1] + header[-1:] + ["1\t10\ta\tA\tG\t12.50\tPASS\t.", "1\t20\tb\tA\tG\t.\tPASS\t.", "1\t30\tc\tA\tG\t30.0\tPASS\t."])
mmg.convert_file(quals, p("q.vcf.parquet"))
mmg.convert_file(p("q.vcf.parquet"), p("q_number.vcf"))
before = mmg.parquet_options()
mmg.parquet_options(numbers="text")
mmg.convert_file(quals, p("q_text.vcf.parquet"))
mmg.convert_file(p("q_text.vcf.parquet"), p("q_text.vcf"))
after = mmg.parquet_options(numbers="number")
try:
    mmg.parquet_options(numbers="float")
    refused = False
except ValueError:
    refused = True
check("parquet_options: QUAL as a number by default, or exactly as written",
      before == {"numbers": "number"} and after == {"numbers": "number"} and refused
      and [l.split("\t")[5] for l in read(p("q_number.vcf")) if not l.startswith("#")] == ["12.5", ".", "30"]
      and read(p("q_text.vcf")) == read(quals), str(read(p("q_text.vcf"))[-3:]))

# ---------------------------------------------------------------- FASTQ and k-mers
ADAPTER = "AGATCGGAAGAGC"
fq_reads = []
for i in range(500):
    s = "".join(random.choice("ACGT") for _ in range(60)) + ADAPTER + "ACGTACGTAC"
    fq_reads += [f"@r{i}", s, "+", "I" * len(s)]
fq = write("x.fastq", fq_reads)
fs = mmg.fastq_stats(fq)
check("fastq_stats: Phred+33, reads and adapter", fs["phred_offset"] == 33 and fs["reads"] == 500 and fs["likely_adapter"]["sequence"] == ADAPTER, str(fs["likely_adapter"]))
tr = mmg.fastq_trim(fq, p("trimmed.fastq.gz"), adapter="auto", min_length=20)
check("fastq_trim: removes the adapter", tr["reads_with_adapter"] == 500 and all(len(l) == 60 for l in read(p("trimmed.fastq.gz"))[1::4]), str(tr))
mmg.fastq_report(fq, p("report.html"))
page = open(p("report.html"), encoding="utf-8").read()
check("fastq_report: writes a self-contained page with the numbers", page.startswith("<!doctype html>") and '"reads":500' in page and "<script src" not in page)
kc = mmg.count_kmers([fq], p("kmers.tsv"), k=21, memory_mb=1)
check("count_kmers: canonical counts add up", kc["total_kmers"] == 500 * (83 - 21 + 1) and sum(int(l.split("\t")[1]) for l in read(p("kmers.tsv"))) == kc["total_kmers"], str(kc)[:120])

# ---------------------------------------------------------------- k-mers
check("128-bit k-mer: encodes and decodes", "ACGTACGTACGTTTGACCA | Bits" in mmg.test_giant_kmer("ACGTACGTACGTTTGACCA"))
inside = [random.getrandbits(62) for _ in range(5000)]
outside = [random.getrandbits(62) | (1 << 62) for _ in range(5000)]
u128 = lambda vs: b"".join(struct.pack("<QQ", v, 0) for v in vs)
bloom = mmg.build_bloom_filter(u128(inside), 1, 5)
check("Bloom filter: not all zero", any(bloom))
check("Bloom filter: CPU and auto give the same filter", bloom == mmg.build_bloom_filter(u128(inside), 1, 5, device="cpu"))
check("query: what went in is not a mutation", mmg.extract_mutations(u128(inside), bloom, 5) == b"")
mut = mmg.extract_mutations(u128(inside[:100] + outside), bloom, 5)
check("query: what did not go in is a mutation", len(mut) // 16 >= 4990 and mut[:16] == struct.pack("<QQ", outside[0], 0), f"{len(mut) // 16}")
x = sorted(random.getrandbits(63) for _ in range(3000)); y = sorted(random.getrandbits(63) for _ in range(2000))
u64 = lambda vs: struct.pack(f"<{len(vs)}Q", *vs)
check("sort_merge_bytes: correct merge", mmg.sort_merge_bytes(u64(x), u64(y)) == u64(sorted(x + y)))
if info["gpu"]:
    check("GPU: filter identical to the CPU's", bloom == mmg.build_bloom_filter(u128(inside), 1, 5, device="gpu"))
    check("GPU: identical merge", mmg.sort_merge_bytes(u64(x), u64(y), device="gpu") == u64(sorted(x + y)))
else:
    try:
        mmg.build_bloom_filter(u128(inside), 1, 5, device="gpu")
        check("no GPU: asking for the GPU raises an error", False, "no error was raised")
    except RuntimeError:
        check("no GPU: asking for the GPU raises an error", True)

print()
print(f"{done - len(failures)} of {done} checks passed")
if failures:
    print("FAILED:", failures)
    sys.exit(1)
