



























































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
for _call, _want in ((lambda: mmg.read_fasta_batch(["no/such/file.fasta"], 31), OSError),
                     (lambda: mmg.read_fastq_batch(["no/such/file.fastq"], 20, 31), OSError),
                     (lambda: mmg.read_fasta_batch([__file__], 0), ValueError)):
    try:
        _call()
        _got = None
    except BaseException as _err:  # a Rust panic would come as PanicException (a BaseException)
        _got = _err
    check(f"k-mer readers: a bad input is a {_want.__name__}, not a crash", isinstance(_got, _want), repr(_got))
from multimerge_genomics import ui as _ui
from multimerge_genomics import exercises as _exercises
check("the class exercises come with the package (mmgen-exercises)",
      set(_exercises.EXERCISES) == {"snps", "broken_vcf", "adapter", "paternity", "expression"})
check("the web interface comes with the package (mmgen-ui)", os.path.exists(_ui.app_path())
      and "@app.cell" in open(_ui.app_path(), encoding="utf-8").read() and "marimo" in _ui.command())

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
# The genome browser's data: the reads of a window, their depth, and the reference bases.
_ref_fa = write("browser_ref.fa", [">chr1", "A" * 99 + "ACGTA" + "C" * 896])
_view = mmg.region_reads(p("x.bam"), "chr1:95-110", reference=_ref_fa)
_r1 = next(r for r in _view["reads"] if r["name"] == "r1")
check("region_reads: the read, its place and its row, the coverage, the reference",
      _r1["start"] == 100 and _r1["end"] == 104 and _r1["mismatches"] == [] and _view["reference"].startswith("AAAAA")
      and any(c[0] == 100 and c[1] == 1 for c in _view["coverage"]), str(_view)[:300])
_dp = mmg.depth_profile(p("x.bam"), "chr1:100-104", bin=5)
check("depth_profile: one read over five bases is depth 1", abs(_dp["regions"][0]["bins"][0][2] - 1.0) < 1e-9, str(_dp)[:200])
check("fetch_sequence", mmg.fetch_sequence(_ref_fa, "chr1:100-104") == "ACGTA")
# The assembly graph: one random sequence is one isolated unitig; two contigs that overlap make one scaffold.
_rnd = random.Random(11)
_seq = "".join(_rnd.choice("ACGT") for _ in range(400))
_g = mmg.unitig_graph(mmg.read_fasta_batch([write("one.fa", [">s", _seq])], 31), 31)
check("unitig_graph: one sequence, one isolated unitig of its length", _g["summary"]["nodes"] == 1 and _g["nodes"][0]["length"] == 400
      and _g["nodes"][0]["kind"] == "isolated", str(_g["summary"]))
_og = mmg.overlap_graph([_seq[:250], _seq[230:]], 20)
check("overlap_graph: two contigs that overlap by 20 bases, joined", len(_og["edges"]) == 1 and _og["edges"][0]["used"], str(_og)[:300])
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
# The engine sorts by radix from 32,768 elements on, and every check above is
# smaller: these make each system of the CI go through that road (64-bit
# values; and the 128-bit k-mers, by their exact key with k = 31 and by a
# prefix key with k = 45).
_big_x = [random.getrandbits(64) for _ in range(100_000)]
_big_y = [random.getrandbits(64) for _ in range(80_000)]
check("sort_merge_bytes, 180,000 unsorted values: the radix road", mmg.sort_merge_bytes(u64(_big_x), u64(_big_y)) == u64(sorted(_big_x + _big_y)))
_code = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}


def _kmer_windows(seq, k):
    out, bits, run, mask = [], 0, 0, (1 << (2 * k)) - 1
    for ch in seq:
        if ch in "Nn":
            bits, run = 0, 0
            continue
        bits = ((bits << 2) | _code[ch]) & mask
        run += 1
        if run >= k:
            out.append(bits)
    return out


_recs = []
for _r in range(2):
    _s = "".join(random.choice("ACGTACGTACGTacgt") for _ in range(70_000))
    _recs.append(_s[:30_000] + "N" * 25 + _s[30_000:])
_bigfa = write("radix.fa", sum(([f">r{_i}"] + [_s[j:j + 60] for j in range(0, len(_s), 60)] for _i, _s in enumerate(_recs)), []))
for _k in (31, 45):
    _want = sorted(set(w for _s in _recs for w in _kmer_windows(_s, _k)))
    _got = mmg.read_fasta_batch([_bigfa], _k)
    check(f"read_fasta_batch, k={_k}, {len(_want):,} k-mers: the radix road gives the same k-mers",
          _got == b"".join(struct.pack("<QQ", w & (2**64 - 1), w >> 64) for w in _want), f"{len(_got) // 16} vs {len(_want)}")
if info["gpu"]:
    check("GPU: filter identical to the CPU's", bloom == mmg.build_bloom_filter(u128(inside), 1, 5, device="gpu"))
    check("GPU: identical merge", mmg.sort_merge_bytes(u64(x), u64(y), device="gpu") == u64(sorted(x + y)))
else:
    try:
        mmg.build_bloom_filter(u128(inside), 1, 5, device="gpu")
        check("no GPU: asking for the GPU raises an error", False, "no error was raised")
    except RuntimeError:
        check("no GPU: asking for the GPU raises an error", True)

# ---------------------------------------------------------------- human genetics, research, microbiology (experimental)
NEW = ["relatedness", "trio_analysis", "runs_of_homozygosity", "genetic_sex", "build_panel", "panel_from_frequencies",
       "estimate_ancestry", "polygenic_score", "clinvar_findings", "gene_regions", "variant_lookup", "pileup",
       "pharmcat_prepare", "pharmcat_run", "count_reads", "differential_expression", "crispr_count", "crispr_test",
       "sketch_references", "gather", "classify_reads", "find_resistance_genes", "genome_distances", "sketch_distances",
       "tree_from_matrix", "viral_mutations", "in_silico_pcr", "amplicon_asvs", "train_classifier", "classify_asvs",
       "alpha_diversity", "beta_diversity", "taxa_composition", "tumor_mutational_burden", "gene_set_test", "pgx_cohort",
       "genotype_association", "association_power", "build_screen_database", "describe_screen_database",
       "contaminant_screen", "build_transcript_index", "transcript_index_info", "quantify_samples", "explain_pseudoalignment",
       "ancestry_pca", "region_reads", "depth_profile", "fetch_sequence", "unitig_graph", "overlap_graph"]
check("every function of the newer modules is there", all(hasattr(mmg, f) for f in NEW), [f for f in NEW if not hasattr(mmg, f)])
rng = random.Random(7)
DNA = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
COMP = str.maketrans("ACGT", "TGCA")
revcomp = lambda s: s.translate(COMP)[::-1]


def mutate(s, rate):
    return "".join(rng.choice("ACGT".replace(c, "")) if rng.random() < rate else c for c in s)


# A trio on 3,000 SNPs: the child gets one allele from each parent.
sites = sorted(rng.sample(range(10_000, 40_000_000), 3000))
genotype = lambda: rng.choice(["0/0", "0/1", "0/1", "1/1"])
father, mother, stranger = ([genotype() for _ in sites] for _ in range(3))
child = ["/".join(sorted([rng.choice(f.split("/")), rng.choice(m.split("/"))])) for f, m in zip(father, mother)]


def person_vcf(name, gts):
    return write(f"{name}.vcf", ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=248956422>",
                                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                                 f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{name}"]
                 + [f"chr1\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT\t{g}" for pos, g in zip(sites, gts)])


vc, vf, vm, vs = (person_vcf(n, g) for n, g in (("child", child), ("father", father), ("mother", mother), ("stranger", stranger)))
tr = mmg.trio_analysis(vc, vf, vm)
check("trio_analysis: the true father is not excluded, no Mendelian error", tr["verdict"] == "not_excluded" and tr["errors"] == 0, str(tr)[:200])
check("trio_analysis: a stranger is excluded", mmg.trio_analysis(vc, vs, vm)["verdict"] == "excluded")
rel = sorted(x["relationship"] for x in mmg.relatedness([vc, vf, vm])["pairs"])
check("relatedness: two parent-child pairs and the unrelated parents", rel == ["parent_offspring", "parent_offspring", "unrelated"], str(rel))
check("the human-genome results carry the notice", "no responsibility" in tr.get("notice", ""))

# RNA-seq: 300 genes, 3 + 3 samples, one gene 8x up in the treated.
names = [f"ctrl_{i}" for i in (1, 2, 3)] + [f"trt_{i}" for i in (1, 2, 3)]
rows = []
for g in range(300):
    mean = rng.uniform(50, 2000)
    rows.append(f"gene{g}\t" + "\t".join(str(max(0, int(rng.gauss(mean * (8 if g == 7 and s.startswith("trt") else 1), mean ** 0.5 * 2))))
                                          for s in names))
counts = write("rna_counts.tsv", ["gene\t" + "\t".join(names)] + rows)
design = write("rna_samples.tsv", ["sample\tcondition"] + [f"{s}\t{s.split('_')[0]}" for s in names])
de = mmg.differential_expression(counts, design, "condition", "ctrl")
pl_ = de["plot"]
up = [g for g, l, q in zip(pl_["gene"], pl_["log2_fold_change"], pl_["padj"]) if l is not None and q is not None and l > 2 and q < 0.01]
check("differential_expression: the planted gene is up, ~8x", up == ["gene7"], str(up))
de_table = p("rna_de.tsv")
mmg.differential_expression(counts, design, "condition", "ctrl", output=de_table)
top_genes = [g for g, _ in sorted(zip(pl_["gene"], pl_["stat"]), key=lambda t: -(t[1] if t[1] is not None else -1e9))[:6]]
sets = write("rna_sets.gmt", ["TOP\tthe six highest statistics\t" + "\t".join(top_genes), "NOTHING\tnone of these genes\tNOPE1\tNOPE2\tNOPE3"])
gs = mmg.gene_set_test(de_table, sets)
by_set = {s["set"]: s for s in gs["sets"]}
check("gene_set_test: a set of the highest statistics goes up; a set of absent genes is not tested",
      by_set["TOP"]["z"] is not None and by_set["TOP"]["z"] > 3 and by_set["TOP"]["fdr"] < 0.05 and by_set["NOTHING"]["pvalue"] is None,
      str(gs["sets"])[:300])
check("gene_set_test: the built-in spinal cord injury collection loads (11 sets)", len(mmg.gene_set_test(de_table)["sets"]) == 11)

# Genotype x outcome: the simulated teaching cohort of the package (a planted
# effect in rs12190259), and the power of a 5-person study.
from multimerge_genomics.demos import association as _assoc_demo
_assoc_demo.write_demo(p("assoc"))
_ga = mmg.genotype_association(p("assoc/cohort.vcf.gz"), p("assoc/phenotypes.tsv"), "recovery_score",
                               covariates=["age", "sex", "center", "anc_afr", "anc_eas"], output=p("assoc/out.tsv"))
_ga_hits = [r["id"] for r in _ga["rows"] if r["fdr"] is not None and r["fdr"] < 0.05]
check("genotype_association: 115 people, 40 alleles, the planted variant the only FDR < 0.05, and the notice",
      _ga["people"]["analysed"] == 115 and _ga["variants"]["tested"] == 40 and _ga_hits == ["rs12190259"]
      and "notice" in _ga and os.path.exists(p("assoc/out.tsv")), f"{_ga['people']} {_ga_hits}")
_ap = mmg.association_power(0.25, effect=2.0, n=5)
check("association_power: 5 people cannot detect an odds ratio of 2", _ap["power"] < 0.1 and _ap["n_needed"] > 100, str(_ap.get("text"))[:200])

# Contaminant screen: the simulated raw material of the package and its blank.
import csv as _csv
from multimerge_genomics.demos import contaminants as _contam_demo
_cd = _contam_demo.write_demo(p("contam"))
_crefs = list(_csv.DictReader(open(p("contam/references.tsv"), encoding="utf-8"), delimiter="\t"))
_cdb = mmg.build_screen_database([p("contam/" + r["file"]) for r in _crefs], p("contam/demo.mmgscreen"),
                                 kinds=[r["kind"] for r in _crefs], names=[r["name"] for r in _crefs])
_cs = mmg.contaminant_screen([("raw_material", p("contam/raw_material_R1.fq.gz"), p("contam/raw_material_R2.fq.gz")),
                              {"name": "reagent_blank", "r1": p("contam/reagent_blank_R1.fq.gz"), "r2": p("contam/reagent_blank_R2.fq.gz")}],
                             p("contam/demo.mmgscreen"), negative_controls=["reagent_blank"], spike_ins=p("contam/spike_ins.tsv"),
                             output_html=p("contam/screen.html"))
_calls = {s["name"]: {o["name"]: o["call"] for o in s["organisms"]} for s in _cs["samples"]}
check("contaminant_screen: the calls the demo plants (detected, background, below threshold), a report and the notice",
      all(_calls[s].get(o) == c for s, exp in _cd["expected"].items() for o, c in exp.items())
      and os.path.exists(p("contam/screen.html")) and "notice" in _cs and len(_cdb["sha256"]) == 64, str(_calls)[:300])

# RNA-seq from FASTQ: the simulated transcriptome, quantified and compared.
from multimerge_genomics.demos import rnaseq_fastq as _quant_demo
_quant_demo.write_demo(p("quant"))
_qi = mmg.build_transcript_index(p("quant/transcripts.fa"), p("quant/transcripts.idx"))
_qs = [{"name": n, "r1": p(f"quant/{n}_R1.fastq.gz"), "r2": p(f"quant/{n}_R2.fastq.gz"), "condition": n.split("_")[0]}
       for n in ("control_1", "control_2", "control_3", "treated_1", "treated_2", "treated_3")]
_qr = mmg.quantify_samples(p("quant/transcripts.idx"), _qs, p("quant/out"))
_qde = mmg.differential_expression(_qr["files"]["gene_counts"], _qr["files"]["samples"], "condition", "control")
check("quantify_samples: 269 transcripts, every pair pseudoaligned, and differential expression finds the 3 up and 3 down",
      _qi["transcripts"] == 269 and min(s["pseudoaligned_pct"] for s in _qr["samples"]) > 99.5 and _qde["up"] == 3 and _qde["down"] == 3,
      f"{_qi['transcripts']} {[s['pseudoaligned_pct'] for s in _qr['samples']]} {_qde['up']} {_qde['down']}")
_qe = mmg.explain_pseudoalignment(p("quant/transcripts.idx"), "GGTCTCCCAACCGATCACGTTGTGTTAAGTGTTTTAACGCGATACGGGTCATTGAGATAGATTTACACGAGACTA")
check("explain_pseudoalignment: a read of Gene80 fits its two isoforms", _qe["status"] == "multi" and _qe["genes"] == ["Gene80"], str(_qe)[:200])

# CRISPR screen: 40 genes x 4 guides; the guides of ESS1 drop out.
crows = []
for g in ["ESS1"] + [f"G{i}" for i in range(39)]:
    for k in range(4):
        before = [rng.randint(300, 700) for _ in range(2)]
        after = [max(1, b // 20) if g == "ESS1" else rng.randint(300, 700) for b in before]
        crows.append(f"{g}_{k}\t{g}\t" + "\t".join(map(str, before + after)))
ct = write("crispr_counts.txt", ["sgRNA\tGene\tc1\tc2\tt1\tt2"] + crows)
cr = mmg.crispr_test(ct, ["t1", "t2"], control=["c1", "c2"])
check("crispr_test: the depleted gene ranks first in negative selection", cr["rank_plot"]["neg"]["gene"][0] == "ESS1", str(cr["rank_plot"]["neg"]["gene"][:3]))

# Microbiology: two genomes and a close relative of the first.
gA, gB = DNA(60_000), DNA(60_000)
gA2 = mutate(gA, 0.01)
fa = {n: write(f"{n}.fa", [f">{n}", s]) for n, s in (("A", gA), ("A2", gA2), ("B", gB))}
db = p("refs.sig")
mmg.sketch_references([[fa["A"]], [fa["B"]]], db, names=["A", "B"], scaled=10)
reads = []
for i in range(2000):
    j = rng.randrange(0, len(gA) - 150)
    reads += [f"@r{i}", gA[j:j + 150], "+", "I" * 150]
ga = mmg.gather([write("sample.fastq", reads)], [db], threshold_bp=5000)
check("gather: the reads come from genome A", [m["name"] for m in ga["matches"]][:1] == ["A"], str([m["name"] for m in ga["matches"]]))
gd = mmg.genome_distances([fa["A"], fa["A2"], fa["B"]], names=["A", "A2", "B"], tree_output=p("tree.nwk"))
d = gd["distance"]
check("genome_distances: A is nearest A2, and the tree is written", d[0][1] < d[0][2] and gd["ani"][0][1] > 0.98 and os.path.exists(p("tree.nwk")), str(d))

# A virus: one substitution in the second codon of a gene (GCT -> TCT, A2S).
ref = list(DNA(3000))
ref[100:106] = list("ATGGCT")
ref = "".join(ref)
virus = ref[:103] + "T" + ref[104:]
vr = write("virus_ref.fa", [">ref1", ref])
vg = write("virus_sample.fa", [">sample1", virus])
vgff = write("virus_ref.gff3", ["##gff-version 3", "ref1\t.\tCDS\t101\t1000\t.\t+\t0\tID=cds1;gene=G1"])
vm_ = mmg.viral_mutations(vr, [vg], annotation=vgff)
g0 = vm_["genomes"][0]
check("viral_mutations: the substitution and its amino-acid change", [s["label"] for s in g0["substitutions"]] == ["G104T"] and
      [a["label"] for a in g0["aa_substitutions"]] == ["G1:A2S"], str(g0)[:300])

# In-silico PCR: a primer pair that gives a 300-base product.
panel = write("panel.tsv", ["target\tforward\treverse\tsize", f"T1\t{ref[200:220]}\t{revcomp(ref[480:500])}\t250-350"])
pc = mmg.in_silico_pcr(panel, [vr])
check("in_silico_pcr: one product of 300 bases", [a["size"] for a in pc["amplicons"]] == [300] and pc["rows"][0]["present"] == [True], str(pc["amplicons"])[:200])

# 16S: a classifier trained on three genera; diversity of a small table.
lineage = {"Lactobacillus": "Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;",
           "Bacteroides": "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;",
           "Escherichia": "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;"}
template = {g: DNA(450) for g in lineage}
ref16 = write("ref16s.fa", [x for g in lineage for i in range(6) for x in (f">{lineage[g]}", mutate(template[g], 0.02))])
model = p("c.model")
mmg.train_classifier(ref16, model)
asv = write("asv.fa", [">ASV_1", mutate(template["Bacteroides"], 0.02)])
cl = mmg.classify_asvs(asv, model, p("taxonomy.tsv"))
check("classify_asvs: the right genus", cl["assignments"][0]["lineage"][-1:] == ["Bacteroides"], str(cl["assignments"][0])[:200])
tab = write("table.tsv", ["ASV\tS1\tS2\tS3\tS4", "ASV_1\t100\t90\t5\t0", "ASV_2\t0\t10\t80\t100", "ASV_3\t50\t50\t50\t50"])
meta = write("meta.tsv", ["sample\tgroup", "S1\tx", "S2\tx", "S3\ty", "S4\ty"])
al = mmg.alpha_diversity(tab)
check("alpha_diversity: observed ASVs per sample", [r["observed"] for r in al["rows"]] == [2, 3, 3, 2], str([r["observed"] for r in al["rows"]]))
be = mmg.beta_diversity(tab, metadata=meta, permutations=99)
check("beta_diversity: Bray-Curtis S1-S4 and PERMANOVA groups", abs(be["distances"][0][3] - (1 - 2 * 50 / 300)) < 1e-9 and be["permanova"]["groups"] == ["x", "y"],
      str(be["distances"][0]))

# Tumor mutational burden: 40 somatic SNVs in a 1-Mb region.
tvcf = write("tumor.vcf", ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=248956422>",
                           '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                           '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
                           '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
                           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR"]
             + [f"chr1\t{1_000_001 + 20_000 * i}\t.\tC\tT\t60\tPASS\t.\tGT:AD:DP\t0/1:20,10:30" for i in range(40)])
tbed = write("panel.bed", ["chr1\t1000000\t2000000"])
tm = mmg.tumor_mutational_burden(tvcf, tbed)
check("tumor_mutational_burden: 40 per Mb", tm["counted"] == 40 and abs(tm["tmb"] - 40.0) < 1e-9, str(tm)[:200])

# Pharmacogenetics (mmgen-pgx v1): native Rust engine using external JSON bundle fixture.
import json as _pgx_json

_fx_tables = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fixtures", "pgx_tables")
_pgx_vcf = write("pgx_phased.vcf", [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=248956422>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA_SYN",
    "chr1\t100\trs1000\tA\tG\t50\tPASS\t.\tGT:PS\t1|0:10",
    "chr1\t120\trs1001\tC\tA\t50\tPASS\t.\tGT:PS\t0|1:10",
    "chr2\t200\trs2000\tC\tT\t50\tPASS\t.\tGT\t0/0",
    "chr3\t300\trs3000\tG\tA\t50\tPASS\t.\tGT:PS\t1|0:30",
    "chr3\t320\trs3001\tT\tC\t50\tPASS\t.\tGT:PS\t0|1:30",
    "chrX\t400\trs4000\tG\tA\t50\tPASS\t.\tGT\t1",
    "chrM\t500\trs5000\tA\t.\t50\tPASS\t.\tGT\t0",
    "chrM\t520\t.\tC\tCC\t50\tPASS\t.\tGT\t1",
])
_pgx_res = mmg.mmgen_pgx_run(_pgx_vcf, p("pgx_out_1"), tables=_fx_tables)
_pgx_res2 = mmg.pgx_run(_pgx_vcf, p("pgx_out_2"), engine="mmgen-pgx", tables=_fx_tables)
_by_g = {g["gene"]: g for g in _pgx_res["genes"]}
check("mmgen_pgx_run: calls phased star allele, variant_activity, haploid chrX and chrM indel",
      _pgx_res["engine"] == "mmgen-pgx"
      and _pgx_res["schema_version"] == "1"
      and _by_g["SYNG_STAR"]["diplotype"] == "*2/*3"
      and _by_g["SYNG_ACT"]["diplotype"] == "varA/varB"
      and _by_g["SYNG_X"]["diplotype"] == "A-"
      and _by_g["SYNG_MITO"]["diplotype"] == "520insC"
      and _by_g["SYNG_MITO"]["called"] is True,
      str(_by_g)[:300])
_pgx_disk = _pgx_json.load(open(_pgx_res["x_files"]["result_json"], encoding="utf-8"))
check("mmgen_pgx_run: returns notice on Python dict without leaking notice into result.json",
      "no responsibility" in _pgx_res.get("notice", "")
      and "notice" not in _pgx_disk
      and all(k in _pgx_disk for k in ("schema_version", "engine", "engine_version", "data_version", "sample", "genes", "drugs", "traits", "warnings", "report_html", "text"))
      and len(_pgx_disk["traits"]) == 2
      and os.path.exists(_pgx_res["report_html"]))
check("pgx_run: engine='mmgen-pgx' agrees with mmgen_pgx_run and omits file-specific differences",
      _pgx_res2["genes"] == _pgx_res["genes"] and _pgx_res2["drugs"] == _pgx_res["drugs"])
# The contract of both engines ships with the package (validated against it in
# test_ui and validate_pgx; here, without jsonschema, its required keys).
with open(mmg.PGX_SCHEMA_PATH, encoding="utf-8") as _fh:
    _pgx_schema = _pgx_json.load(_fh)
check("pgx_result.schema.json ships with the package, and the result has every key it requires",
      _pgx_schema.get("properties", {}).get("schema_version", {}).get("const") == "1"
      and all(k in _pgx_disk for k in _pgx_schema["required"])
      and all(k in g for g in _pgx_disk["genes"] for k in _pgx_schema["$defs"]["gene"]["required"]))
try:
    _b37 = write("pgx_grch37.vcf", ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=249250621>",
                                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"])
    mmg.mmgen_pgx_run(_b37, p("pgx_out_err"), tables=_fx_tables)
    check("mmgen_pgx_run: GRCh37 VCF raises error", False, "no error raised")
except RuntimeError as _e:
    check("mmgen_pgx_run: GRCh37 VCF raises clear error", "GRCh37" in str(_e) and "GRCh38" in str(_e), str(_e))
_coh = mmg.pgx_cohort([_pgx_vcf], p("pgx_cohort"), tables=_fx_tables, drugs="syndrug_a, gabapentin")
check("pgx_cohort: one person, a medicine with a recommendation and one without a guideline, and the notice",
      any(m["medicine"] == "syndrug_a" and m["status"] == "recommendation" for m in _coh["medicines"])
      and _coh["medicines_without_guideline"] == ["gabapentin"] and os.path.exists(_coh["files"]["medicines"])
      and "no responsibility" in _coh.get("notice", ""), str(_coh)[:300])

# Pharmacogenetics v2 (BAM+VCF CYP2D6 + Class I HLA + outside_calls.tsv)
def _synth_hla_dna(seed, length=240):
    x = seed ^ 0x9E3779B97F4A7C15
    bases = "ACGT"
    out = []
    for _ in range(length):
        x = (x * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)
        out.append(bases[(x >> 32) & 3])
    return "".join(out)

_hla_py_dir = p("pgx_v2_hla")
os.makedirs(_hla_py_dir, exist_ok=True)
_hla_specs = {
    "A_gen.fasta": [("A*01:01:01:01", 101), ("A*31:01:02:01", 202), ("A*02:01:01:01", 303)],
    "B_gen.fasta": [("B*07:02:01:01", 404), ("B*57:01:01:01", 505), ("B*35:01:01:01", 606)],
    "C_gen.fasta": [("C*04:01:01:01", 707), ("C*12:03:01:01", 808), ("C*07:01:01:01", 909)],
}
_hla_seqs = {}
for _fn, _als in _hla_specs.items():
    _lines = []
    _hla_seqs[_fn] = []
    for _i, (_al, _sd) in enumerate(_als):
        _s = _synth_hla_dna(_sd, 240)
        _hla_seqs[_fn].append(_s)
        _lines.append(f">HLA:HLA_PY_{_fn[0]}{_i} {_al} 240 bp\n{_s}\n")
    with open(os.path.join(_hla_py_dir, _fn), "w", encoding="utf-8") as _fh:
        _fh.writelines(_lines)

_v2_sam_lines = [
    "@HD\tVN:1.6\tSO:unsorted",
    "@SQ\tSN:chr4\tLN:190214555",
    "@SQ\tSN:chr6\tLN:170805979",
    "@SQ\tSN:chr22\tLN:50818468",
]
for _cpos in range(73404287, 73420484, 1000):
    _clen = min(1000, 73421482 - _cpos + 1)
    _cseq = "A" * _clen
    _cqual = "I" * _clen
    for _cp in range(20):
        _v2_sam_lines.append(f"ctrl_{_cpos}_{_cp}\t0\tchr4\t{_cpos}\t60\t{_clen}M\t*\t0\t0\t{_cseq}\t{_cqual}")

for _tag, _bpos, _s1, _s2 in [
    ("hlaA", 29936100, _hla_seqs["A_gen.fasta"][0], _hla_seqs["A_gen.fasta"][1]),
    ("hlaC", 31265100, _hla_seqs["C_gen.fasta"][0], _hla_seqs["C_gen.fasta"][1]),
    ("hlaB", 31350100, _hla_seqs["B_gen.fasta"][0], _hla_seqs["B_gen.fasta"][1]),
]:
    for _hi, _seq in enumerate((_s1, _s2), start=1):
        for _p in range(4):
            _qn = f"{_tag}_{_hi}_{_p}"
            _v2_sam_lines.append(f"{_qn}\t65\tchr6\t{_bpos + _p}\t60\t100M\t*\t0\t0\t{_seq[:100]}\t{'I'*100}")
            _v2_sam_lines.append(f"{_qn}\t129\tchr6\t{_bpos + 140 + _p}\t60\t100M\t*\t0\t0\t{_seq[140:240]}\t{'I'*100}")

_s2 = [b"G", b"A", b"C", b"C", b"T", b"G", b"G"]
_s4 = [b"C", b"G", b"T", b"C", b"T", b"G", b"A"]
for _htag, _hb in (("h1", _s2), ("h2", _s4)):
    for _r in range(10):
        _ex9 = bytearray(b"A" * 301)
        _ex9[42126611 - 42126500] = _hb[0][0]
        _mid = bytearray(b"A" * 1199)
        _mid[42127941 - 42126801] = _hb[1][0]
        _fivep = bytearray(b"A" * 2811)
        _fivep[42128945 - 42128000] = _hb[2][0]
        _fivep[42129130 - 42128000] = _hb[3][0]
        _fivep[42129809 - 42128000] = _hb[4][0]
        _fivep[42129819 - 42128000] = _hb[5][0]
        _fivep[42130692 - 42128000] = _hb[6][0]
        _qn = f"{_htag}_{_r}"
        _v2_sam_lines.append(f"{_qn}_e\t65\tchr22\t42126500\t60\t301M\t*\t0\t0\t{_ex9.decode()}\t{'I'*301}")
        _v2_sam_lines.append(f"{_qn}_b\t65\tchr22\t42126801\t60\t1199M\t*\t0\t0\t{_mid.decode()}\t{'I'*1199}")
        _v2_sam_lines.append(f"{_qn}_b\t129\tchr22\t42128000\t60\t2811M\t*\t0\t0\t{_fivep.decode()}\t{'I'*2811}")

_v2_sam = write("pgx_v2.sam", _v2_sam_lines)
_v2_bam = p("pgx_v2.bam")
mmg.sort_genomic([_v2_sam], _v2_bam, index=True)
_v2_vcf = write("pgx_v2.vcf", [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr22,length=50818468>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA_V2",
    "chr22\t42126611\t.\tC\tG\t50\tPASS\t.\tGT\t0/1",
    "chr22\t42127941\t.\tG\tA\t50\tPASS\t.\tGT\t0/1",
    "chr22\t42128945\t.\tC\tT\t50\tPASS\t.\tGT\t0/1",
    "chr22\t42130692\t.\tG\tA\t50\tPASS\t.\tGT\t0/1",
])

_old_hla_dir = os.environ.get("MMGEN_HLA_DIR")
_old_hla_off = os.environ.get("MMGEN_HLA_OFFLINE")
os.environ["MMGEN_HLA_DIR"] = _hla_py_dir
os.environ["MMGEN_HLA_OFFLINE"] = "1"
try:
    _v2_res = mmg.mmgen_pgx_run(_v2_vcf, p("pgx_v2_out"), tables=_fx_tables, bam=_v2_bam)
    _v2_wrap = mmg.pgx_run(_v2_vcf, p("pgx_v2_wrap"), engine="mmgen-pgx", tables=_fx_tables, bam=_v2_bam)
finally:
    if _old_hla_dir is None:
        os.environ.pop("MMGEN_HLA_DIR", None)
    else:
        os.environ["MMGEN_HLA_DIR"] = _old_hla_dir
    if _old_hla_off is None:
        os.environ.pop("MMGEN_HLA_OFFLINE", None)
    else:
        os.environ["MMGEN_HLA_OFFLINE"] = _old_hla_off

_v2_by_g = {g["gene"]: g for g in _v2_res["genes"]}
_v2_oc = _v2_res["x_files"].get("pharmcat_outside_calls", "")
_v2_oc_txt = open(_v2_oc, encoding="utf-8").read() if _v2_oc and os.path.exists(_v2_oc) else ""
_v2_drugs = {d["drug"] for d in _v2_res["drugs"]}
check(
    "mmgen_pgx_run v2 BAM: calls CYP2D6 *2/*4 CN=2, HLA-A/B/C, drug guidelines & outside_calls.tsv",
    _v2_by_g.get("CYP2D6", {}).get("diplotype") == "*2/*4"
    and _v2_by_g.get("CYP2D6", {}).get("copy_number") == 2
    and _v2_by_g.get("CYP2D6", {}).get("source") == "vcf+bam"
    and _v2_by_g.get("CYP2D6", {}).get("x_phased") is True
    and _v2_by_g.get("HLA-A", {}).get("diplotype") == "*01:01/*31:01"
    and _v2_by_g.get("HLA-B", {}).get("diplotype") == "*07:02/*57:01"
    and _v2_by_g.get("HLA-C", {}).get("diplotype") == "*04:01/*12:03"
    and {"abacavir", "carbamazepine"}.issubset(_v2_drugs)
    and "CYP2D6\t*2/*4\n" in _v2_oc_txt
    and "HLA-A\t*01:01/*31:01\n" in _v2_oc_txt
    and "HLA-B\t*07:02/*57:01\n" in _v2_oc_txt
    and "HLA-C" not in _v2_oc_txt
    and _v2_wrap["genes"] == _v2_res["genes"]
    and "pharmcat_outside_calls" not in _pgx_res["x_files"],
    str(_v2_by_g.get("CYP2D6"))[:200] + " | " + repr(_v2_oc_txt),
)

print()
print(f"{done - len(failures)} of {done} checks passed")
if failures:
    print("FAILED:", failures)
    sys.exit(1)
