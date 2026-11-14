# Cross-check of the file formats and tools against independent software:
# what this library writes is read back by other readers, and what other
# writers produce is read by this library. No network, no external data.
#
#     pip install pyarrow polars biopython bamnostic cutadapt
#     python trials/python/crosscheck.py
#
# Each group needs its reference library; a missing one is reported as
# "skipped" and the rest still runs. The Linux-only references (samtools,
# bcftools, bedtools, tabix, jellyfish) are used when they are on the PATH.
import collections
import gzip
import importlib.util
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile

import multimerge_genomics as mmg

failures, skipped = [], []
done = 0


def check(name, cond, detail=""):
    global done
    done += 1
    if not cond:
        failures.append(name)
    print(("  ok    " if cond else "  FAIL  ") + name + ("" if cond else f"  {detail}"))


def have(module):
    if importlib.util.find_spec(module) is None:
        skipped.append(module)
        print(f"  skip  {module} is not installed")
        return False
    return True


def tool(name):
    path = shutil.which(name)
    if not path:
        skipped.append(name)
        print(f"  skip  {name} is not on the PATH")
    return path


random.seed(2026)
T = tempfile.mkdtemp(prefix="mmg_cross_")
p = lambda name: os.path.join(T, name)


def write(name, lines):
    with open(p(name), "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return p(name)


def lines_of(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", newline="") as f:
        return [l.rstrip("\r\n") for l in f if l.strip()]


def base(n):
    return "".join(random.choice("ACGT") for _ in range(n))


# ---------------------------------------------------------------- test data
CHROMS = [("chr1", 2_000_000), ("chr2", 1_500_000), ("chr10", 900_000), ("chrX", 700_000)]
vcf_header = ["##fileformat=VCFv4.2"] + [f"##contig=<ID={c},length={n}>" for c, n in CHROMS] + [
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
    '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
]
vcf_records = []
for i in range(30000):
    c, n = random.choice(CHROMS)
    ref = base(random.choice([1, 1, 1, 2, 4]))
    alts = ",".join(base(random.choice([1, 1, 3])) for _ in range(random.choice([1, 1, 2])))
    qual = random.choice([".", str(random.randrange(1, 99)), f"{random.uniform(0, 99):.2f}"])
    info = f"DP={random.randrange(1, 300)};AF={','.join(str(round(random.random(), 3)) for _ in alts.split(','))}" + (";DB" if random.random() < 0.3 else "")
    gts = "\t".join(random.choice(["0/0", "0/1", "1/1", "./.", "0|1"]) for _ in range(2))
    vcf_records.append(f"{c}\t{random.randrange(1, n)}\trs{i}\t{ref}\t{alts}\t{qual}\t{random.choice(['PASS', 'q10', '.'])}\t{info}\tGT\t{gts}")
vcf_in = write("in.vcf", vcf_header + vcf_records)
natural = lambda c: (0, int(c[3:])) if c[3:].isdigit() else (1, c)
order = {c: i for i, (c, _) in enumerate(CHROMS)}
vcf_sorted = sorted(vcf_records, key=lambda l: (order[l.split("\t")[0]], int(l.split("\t")[1])))

# ---------------------------------------------------------------- 1. BGZF
print("\nBGZF (.gz written by this library)")
gz = p("sorted.vcf.gz")
mmg.sort_genomic([vcf_in], gz, memory_mb=1)
check("the .gz reads back with Python's gzip", lines_of(gz) == vcf_header + vcf_sorted)
if have("Bio"):
    from Bio import bgzf

    with open(gz, "rb") as f:
        blocks = list(bgzf.BgzfBlocks(f))
    check("Biopython walks every BGZF block", len(blocks) > 3 and all(b[1] <= 65536 for b in blocks), f"{len(blocks)} blocks")
    check("the last block is the BGZF end marker", blocks[-1][3] == 0)
    with bgzf.BgzfReader(gz, "rb") as r:
        parts = []
        while chunk := r.read(1 << 20):
            parts.append(chunk)
        text = b"".join(parts).decode()
    check("Biopython's BGZF reader gives the same text", text.splitlines() == vcf_header + vcf_sorted)
    # A BGZF file written by Biopython is read by this library.
    bio_gz = p("bio.vcf.gz")
    with bgzf.BgzfWriter(bio_gz, "wb") as w:
        w.write(("\n".join(vcf_header + vcf_records) + "\n").encode())
    n, _, _ = mmg.sort_genomic([bio_gz], p("from_bio.vcf"))
    check("a BGZF file written by Biopython sorts like the plain text", lines_of(p("from_bio.vcf")) == vcf_header + vcf_sorted and n == len(vcf_records))
plain_gz = p("plain.vcf.gz")
with gzip.open(plain_gz, "wt", newline="\n") as f:
    f.write("\n".join(vcf_header + vcf_records) + "\n")
r = mmg.validate_file(plain_gz)
check("a plain (non-BGZF) gzip is still read", r["records"] == len(vcf_records), str(r["records"]))
try:
    mmg.index_file(plain_gz)
    check("indexing a plain gzip is refused", False, "no error")
except RuntimeError as e:
    check("indexing a plain gzip is refused, with the way out", "BGZF" in str(e) and "bgzip" in str(e), str(e)[:100])

# ---------------------------------------------------------------- 2. tabix index and region queries
print("\ntabix index and region queries")
info = mmg.index_file(gz)
check("index_file writes x.vcf.gz.tbi", info["kind"] == "tbi" and os.path.exists(gz + ".tbi") and info["records"] == len(vcf_records), str(info))


def span(line):
    f = line.split("\t")
    beg = int(f[1]) - 1
    return f[0], beg, beg + len(f[3])


def scan(records, chrom, beg, end):
    return [l for l in records if (lambda c, b, e: c == chrom and b < end and max(e, b + 1) > beg)(*span(l))]


ok = True
for _ in range(300):
    c, n = random.choice(CHROMS)
    a = random.randrange(1, n)
    e = a + random.choice([0, 10, 1000, 50_000])
    got = mmg.query_region(gz, f"{c}:{a}-{e}")
    if got != scan(vcf_sorted, c, a - 1, e):
        ok = False
        print("    differs at", f"{c}:{a}-{e}")
        break
check("300 random region queries equal a linear scan", ok)
got = mmg.query_region(gz, ["chr2:1-200000", "chrX", "chr2:150000-300000"], output=p("regions.vcf"))
want = scan(vcf_sorted, "chr2", 0, 300000) + scan(vcf_sorted, "chrX", 0, 10**12)
check("several regions to a file: header, then each record once in file order", lines_of(p("regions.vcf")) == vcf_header + want and got == len(want))
tabix = tool("tabix")
if tabix:
    out = subprocess.run([tabix, gz, "chr1:100000-300000"], capture_output=True, text=True).stdout.splitlines()
    check("tabix reads our .tbi and gives the same records", out == mmg.query_region(gz, "chr1:100000-300000"))

# ---------------------------------------------------------------- 3. BAM
print("\nBAM")
sam_header = ["@HD\tVN:1.6\tSO:coordinate"] + [f"@SQ\tSN:{c}\tLN:{n}" for c, n in CHROMS] + ["@RG\tID:g1\tSM:s1"]
sam_records = []
for i in range(20000):
    c, n = random.choice(CHROMS)
    L = random.randrange(20, 151)
    pos = random.randrange(1, n - 500)
    cigar = random.choice([f"{L}M", f"5S{L - 5}M", f"{L // 2}M3D{L - L // 2}M", f"{L // 2}M2I{L - L // 2 - 2}M", f"{L // 2}M1000N{L - L // 2}M"])
    unmapped = random.random() < 0.03
    flag = 4 if unmapped else random.choice([0, 16, 99, 147, 83, 163])
    qual = "".join(chr(33 + random.randrange(2, 41)) for _ in range(L))
    tags = f"\tNM:i:{random.randrange(0, 9)}\tAS:i:{random.randrange(0, 150)}\tRG:Z:g1" + (f"\tXS:f:{random.choice(['1.5', '-0.25', '3e+10'])}" if random.random() < 0.2 else "")
    sam_records.append(f"r{i}\t{flag}\t{c}\t{pos}\t{0 if unmapped else random.randrange(0, 61)}\t{'*' if unmapped else cigar}\t=\t{pos + 200}\t{0 if unmapped else 350}\t{base(L)}\t{qual}{tags}")
sam_in = write("in.sam", sam_header + sam_records)
bam = p("sorted.bam")
mmg.sort_genomic([sam_in], bam, memory_mb=1, index=True)
check("sort_genomic(..., index=True) writes the BAM and its .bai", os.path.exists(bam) and os.path.exists(bam + ".bai"))
mmg.convert_file(bam, p("back.sam"))
sam_sorted = sorted(sam_records, key=lambda l: (order[l.split("\t")[2]], int(l.split("\t")[3])))
back = lines_of(p("back.sam"))
check("BAM back to SAM: the same alignments, sorted", [l for l in back if not l.startswith("@")] == sam_sorted)
check("BAM back to SAM: @HD says SO:coordinate", back[0].startswith("@HD") and "SO:coordinate" in back[0], back[0])
if have("bamnostic"):
    import bamnostic as bs

    f = bs.AlignmentFile(bam, "rb")
    reads = list(f)
    ok, first_bad = len(reads) == len(sam_sorted), None
    for r, l in zip(reads, sam_sorted):
        x = l.split("\t")
        same = (r.query_name == x[0] and r.flag == int(x[1]) and r.reference_start == int(x[3]) - 1 and r.mapping_quality == int(x[4])
                and (r.cigarstring or "*") == x[5] and r.query_sequence == x[9])
        if not same:
            ok, first_bad = False, (x[0], r.query_name, r.flag, r.reference_start, r.cigarstring)
            break
    check("bamnostic (an independent BAM reader) reads every field back", ok, str(first_bad))
    f.close()
    # bamnostic does not stretch a spliced read over its skipped region (the
    # N of the CIGAR), so it misses such reads in a region that falls inside
    # the intron; samtools, and this library, include them (the SAM spec
    # counts N on the reference). Those reads are left out of this comparison
    # and checked against samtools below, when it is installed.
    spliced = {l.split("\t")[0] for l in sam_records if "N" in l.split("\t")[5]}
    f = bs.AlignmentFile(bam, "rb", index_filename=bam + ".bai")
    ok = True
    for _ in range(100):
        c, n = random.choice(CHROMS)
        a = random.randrange(0, n)
        e = a + random.choice([1, 500, 20000])
        theirs = sorted(r.query_name for r in f.fetch(c, a, e) if r.query_name not in spliced)
        ours = sorted(x for x in (l.split("\t")[0] for l in mmg.query_region(bam, f"{c}:{a + 1}-{e}")) if x not in spliced)
        if theirs != ours:
            ok = False
            print("    differs at", c, a, e, len(theirs), len(ours))
            break
    check("bamnostic uses our .bai: 100 region queries agree with ours (spliced reads aside)", ok)
    f.close()
samtools = tool("samtools")
if samtools:
    r = subprocess.run([samtools, "quickcheck", "-v", bam], capture_output=True, text=True)
    check("samtools quickcheck accepts the BAM", r.returncode == 0, r.stderr[:200])
    out = subprocess.run([samtools, "view", bam], capture_output=True, text=True).stdout.splitlines()
    check("samtools view prints the same SAM lines", out == sam_sorted, f"{len(out)} lines")
    idx = subprocess.run([samtools, "idxstats", bam], capture_output=True, text=True).stdout
    check("samtools idxstats reads our .bai", "chr1" in idx, idx[:200])
    ok = True
    for _ in range(100):
        c, n = random.choice(CHROMS)
        a = random.randrange(1, n)
        reg = f"{c}:{a}-{a + random.choice([0, 500, 20000])}"
        theirs = subprocess.run([samtools, "view", bam, reg], capture_output=True, text=True).stdout.splitlines()
        if theirs != mmg.query_region(bam, reg):
            ok = False
            print("    differs at", reg)
            break
    check("samtools view with our .bai: 100 region queries give the same alignments", ok)

# ---------------------------------------------------------------- 4. Parquet
print("\nParquet")
pq = p("variants.vcf.parquet")
info = mmg.convert_file(gz, pq)
check("convert_file to Parquet counts every record", info["records"] == len(vcf_records), str(info))
mmg.convert_file(pq, p("back.vcf"))


def number_as_stored(text):
    # QUAL is stored as a number, so it comes back in its shortest form:
    # 12.50 as 12.5, 30.0 as 30.
    if text == ".":
        return text
    r = repr(float(text))
    return r[:-2] if r.endswith(".0") else r


expected = ["\t".join(f[:5] + [number_as_stored(f[5])] + f[6:]) for f in (l.split("\t") for l in vcf_sorted)]
check("Parquet back to VCF: the same lines, header included (QUAL as a number)", lines_of(p("back.vcf")) == vcf_header + expected)
if have("pyarrow"):
    import pyarrow.parquet as papq

    t = papq.read_table(pq)
    check("pyarrow reads it: one column per VCF field and per sample",
          t.column_names == ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1", "S2"], str(t.column_names))
    check("pyarrow: POS is int64 and QUAL float64 with nulls for '.'",
          str(t.schema.field("POS").type) == "int64" and str(t.schema.field("QUAL").type) == "double"
          and t.column("QUAL").null_count == sum(1 for l in vcf_records if l.split("\t")[5] == "."))
    meta = papq.read_metadata(pq).metadata
    check("pyarrow: the header is in the file metadata", meta[b"mmgen.header"].decode().split("\n") == vcf_header and meta[b"mmgen.format"] == b"vcf")
    check("pyarrow: positions in file order", t.column("POS").to_pylist() == [int(l.split("\t")[1]) for l in vcf_sorted])
    # A table written by pyarrow, with its own types and names, is read by this library.
    import pyarrow as pa

    bed = [(c, s, s + random.randrange(1, 900), f"f{i}") for i, (c, s) in enumerate((random.choice(CHROMS)[0], random.randrange(0, 600000)) for _ in range(5000))]
    tbl = pa.table({"chrom": pa.array([b[0] for b in bed]).dictionary_encode(), "start": pa.array([b[1] for b in bed], pa.int32()),
                    "end": pa.array([b[2] for b in bed], pa.uint32()), "name": [b[3] for b in bed]})
    papq.write_table(tbl, p("pyarrow.bed.parquet"))
    n, _, _ = mmg.sort_genomic([p("pyarrow.bed.parquet")], p("from_pyarrow.bed"))
    want = [f"{c}\t{s}\t{e}\t{nm}" for c, s, e, nm in sorted(bed, key=lambda b: (b[0], b[1], b[2]))]
    check("a BED table written by pyarrow (int32, dictionary text) sorts correctly", lines_of(p("from_pyarrow.bed")) == want and n == len(bed))
if have("polars"):
    import polars as pl

    df = pl.read_parquet(pq)
    check("polars reads it: same rows and values", df.height == len(vcf_records) and df["ID"].to_list() == [l.split("\t")[2] for l in vcf_sorted])
    gff = pl.DataFrame({"seqid": ["chr2", "chr1", "chr1"], "source": ["x"] * 3, "type": ["gene", "gene", "exon"],
                        "start": [500, 100, 100], "end": [900, 400, 200], "score": [None, 1.5, None],
                        "strand": ["+", "-", "-"], "phase": [".", ".", "0"], "attributes": ["ID=a", "ID=b", "Parent=b"]})
    gff.write_parquet(p("polars.gff.parquet"))
    mmg.sort_genomic([p("polars.gff.parquet")], p("from_polars.gff"))
    got = [l for l in lines_of(p("from_polars.gff")) if not l.startswith("#")]
    check("a GFF table written by polars sorts correctly (parent before child)",
          got == ["chr1\tx\tgene\t100\t400\t1.5\t-\t.\tID=b", "chr1\tx\texon\t100\t200\t.\t-\t0\tParent=b", "chr2\tx\tgene\t500\t900\t.\t+\t.\tID=a"], str(got))
r = mmg.validate_file(pq)
check("validate_file reads a Parquet VCF", r["format"] == "vcf" and r["records"] == len(vcf_records), str({k: r[k] for k in ("format", "records")}))
# The option: QUAL kept as the text it was, for a byte-for-byte round trip.
mmg.parquet_options(numbers="text")
try:
    mmg.convert_file(gz, p("text.vcf.parquet"))
finally:
    mmg.parquet_options(numbers="number")
mmg.convert_file(p("text.vcf.parquet"), p("back_text.vcf"))
check("parquet_options(numbers='text'): Parquet back to VCF gives the original lines", lines_of(p("back_text.vcf")) == vcf_header + vcf_sorted)
if have("pyarrow"):
    t = papq.read_table(p("text.vcf.parquet"))
    check("pyarrow: with numbers='text', QUAL is a text column with the original values",
          str(t.schema.field("QUAL").type) in ("large_string", "string") and t.column("QUAL").to_pylist() == [l.split("\t")[5] for l in vcf_sorted])

# ---------------------------------------------------------------- 5. FASTQ against cutadapt
print("\nFASTQ trimming against cutadapt")
ADAPTER = "AGATCGGAAGAGC"
reads = []
for i in range(4000):
    L = random.randrange(40, 151)
    s = list(base(L))
    if random.random() < 0.4:
        cut = random.randrange(0, L)
        tail = (ADAPTER * 20)[: L - cut]
        s[cut:] = tail
    q = [random.randrange(2, 41) if j < L - 15 else random.randrange(2, 20) for j in range(L)]
    reads.append((f"read{i}", "".join(s), "".join(chr(33 + v) for v in q)))
fq = write("in.fastq", [x for r in reads for x in (f"@{r[0]}", r[1], "+", r[2])])
stats = mmg.fastq_stats(fq)
check("fastq_stats finds the planted adapter", stats["likely_adapter"] and stats["likely_adapter"]["sequence"] == ADAPTER, str(stats["likely_adapter"]))
check("fastq_stats counts reads and bases", stats["reads"] == len(reads) and stats["bases"] == sum(len(r[1]) for r in reads))
r = mmg.fastq_trim(fq, p("trimmed.fastq"), adapter="auto", quality_cutoff=20, min_length=20)
check("fastq_trim(adapter='auto') uses the detected adapter", r["adapter"] == ADAPTER, str(r["adapter"]))
cutadapt = shutil.which("cutadapt") or (os.path.join(os.path.dirname(sys.executable), "cutadapt.exe") if os.path.exists(os.path.join(os.path.dirname(sys.executable), "cutadapt.exe")) else None)
if cutadapt and importlib.util.find_spec("cutadapt"):
    cmd = [sys.executable, "-m", "cutadapt", "--no-indels", "-j", "1", "-a", ADAPTER, "-q", "20", "-m", "20", "-o", p("cutadapt.fastq"), fq]
    subprocess.run(cmd, capture_output=True, check=True)
    check("the result is identical to cutadapt --no-indels", lines_of(p("trimmed.fastq")) == lines_of(p("cutadapt.fastq")))
else:
    skipped.append("cutadapt")
    print("  skip  cutadapt is not installed")
html = p("report.html")
mmg.fastq_report(fq, html)
page = open(html, encoding="utf-8").read()
check("fastq_report writes a self-contained page with the numbers", page.startswith("<!doctype html>") and '"reads": ' not in page[:100]
      and f'"reads":{len(reads)}' in page.replace(" ", "") and "http" not in page.split("<script")[1].split(">")[0])

# ---------------------------------------------------------------- 6. k-mer counting
print("\nk-mer counting")
fa = write("genome.fa", [">a", base(5000), base(3000), ">b", base(4000)])
seqs = []
for l in lines_of(fa):
    if l.startswith(">"):
        seqs.append("")
    else:
        seqs[-1] += l
comp = str.maketrans("ACGT", "TGCA")
for k in (5, 21, 31):
    want = collections.Counter(min(s[i:i + k], s[i:i + k].translate(comp)[::-1]) for s in seqs for i in range(len(s) - k + 1))
    r = mmg.count_kmers([fa], p(f"k{k}.tsv"), k=k, memory_mb=1)
    got = {l.split("\t")[0]: int(l.split("\t")[1]) for l in lines_of(p(f"k{k}.tsv"))}
    check(f"count_kmers k={k}: identical to a Python count of canonical k-mers", got == want and r["distinct_kmers"] == len(want))
jellyfish = tool("jellyfish")
if jellyfish:
    subprocess.run([jellyfish, "count", "-m", "21", "-s", "1M", "-C", "-o", p("jf.jf"), fa], check=True)
    out = subprocess.run([jellyfish, "dump", "-c", p("jf.jf")], capture_output=True, text=True).stdout.split("\n")
    jf = {l.split(" ")[0]: int(l.split(" ")[1]) for l in out if l}
    ours = {l.split("\t")[0]: int(l.split("\t")[1]) for l in lines_of(p("k21.tsv"))}
    check("jellyfish -C gives the same canonical counts", jf == ours)

# ---------------------------------------------------------------- 7. BED operations and VCF tools against bedtools/bcftools
bedtools = tool("bedtools")
if bedtools:
    print("\nBED operations against bedtools")
    beds = []
    for name in ("a", "b"):
        rows = sorted((c, s, s + random.randrange(1, 3000)) for c, s in ((random.choice(CHROMS)[0], random.randrange(0, 600000)) for _ in range(4000)))
        beds.append(write(f"{name}.bed", [f"{c}\t{s}\t{e}" for c, s, e in rows]))
    genome = write("genome.txt", [f"{c}\t{n}" for c, n in sorted(CHROMS)])
    run = lambda *args: subprocess.run([bedtools, *args], capture_output=True, text=True, check=True).stdout.splitlines()
    mmg.bed_subtract([beds[0]], [beds[1]], p("sub.bed"))
    check("bed_subtract equals bedtools subtract", lines_of(p("sub.bed")) == run("subtract", "-a", beds[0], "-b", beds[1]))
    mmg.bed_closest([beds[0]], [beds[1]], p("closest.bed"))
    check("bed_closest equals bedtools closest -d", lines_of(p("closest.bed")) == run("closest", "-d", "-a", beds[0], "-b", beds[1]))
    mmg.bed_coverage([beds[0]], [beds[1]], p("cov.bed"))
    check("bed_coverage equals bedtools coverage", lines_of(p("cov.bed")) == run("coverage", "-a", beds[0], "-b", beds[1]))
    mmg.bed_complement([beds[0]], genome, p("comp.bed"))
    check("bed_complement equals bedtools complement", lines_of(p("comp.bed")) == run("complement", "-i", beds[0], "-g", genome))
    mmg.bed_slop([beds[0]], genome, p("slop.bed"), both=250)
    check("bed_slop equals bedtools slop -b", lines_of(p("slop.bed")) == run("slop", "-b", "250", "-i", beds[0], "-g", genome))
bcftools = tool("bcftools")
if bcftools:
    print("\nVCF tools against bcftools")
    run = lambda *args: [l for l in subprocess.run([bcftools, *args], capture_output=True, text=True, check=True).stdout.splitlines() if not l.startswith("#")]
    # bcftools rewrites numbers on output (QUAL 76.40 as 76.4, AF=1.0 as AF=1);
    # vcf_filter copies the lines as they were. The records are compared by ID.
    ids = lambda lines: [l.split("\t")[2] for l in lines]
    mmg.vcf_filter(gz, p("f.vcf"), expression="QUAL>=30 && INFO/DP>100")
    mine = [l for l in lines_of(p("f.vcf")) if not l.startswith("#")]
    theirs = run("view", "-H", "-i", "QUAL>=30 && INFO/DP>100", gz)
    check("vcf_filter keeps the records bcftools view -i keeps, in the same order", ids(mine) == ids(theirs) and len(mine) > 1000,
          f"{len(mine)} against {len(theirs)}")
    mmg.vcf_filter(gz, p("r.vcf"), regions="chr2:100000-400000")
    mine = [l for l in lines_of(p("r.vcf")) if not l.startswith("#")]
    theirs = run("view", "-H", "-r", "chr2:100000-400000", gz)
    check("vcf_filter(regions=...) keeps the records bcftools view -r keeps, in the same order", ids(mine) == ids(theirs) and len(mine) > 100,
          f"{len(mine)} against {len(theirs)}")

    # TYPE and FILTER, with == and ~. The alleles are ones both tools classify
    # the same way (bcftools trims before classifying, so AC>AT is a SNP to it
    # and an MNP here): SNPs on a one-base REF, pure insertions and deletions,
    # and MNPs whose every base differs.
    def allele(ref):
        kind = random.choice(["snp", "ins"] if len(ref) == 1 else ["del", "ins", "mnp"])
        if kind == "snp":
            return random.choice([b for b in "ACGT" if b != ref])
        if kind == "ins":
            return ref + base(random.randrange(1, 4))
        if kind == "del":
            return ref[:random.randrange(1, len(ref))]
        return "".join(random.choice([b for b in "ACGT" if b != x]) for x in ref)

    typed = []
    for i in range(6000):
        ref = base(random.choice([1, 1, 2, 3]))
        alts = "." if random.random() < 0.05 else ",".join(dict.fromkeys(allele(ref) for _ in range(random.choice([1, 1, 2, 3]))))
        filt = random.choice(["PASS", ".", "q10", "s50", "q10;s50"])
        typed.append(f"chr1\t{1000 + 10 * i}\tt{i}\t{ref}\t{alts}\t30\t{filt}\t.")
    tvcf = write("typed.vcf", ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=2000000>",
                               '##FILTER=<ID=q10,Description="q10">', '##FILTER=<ID=s50,Description="s50">',
                               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"] + typed)
    for expr in ['TYPE="snp"', 'TYPE~"snp"', 'TYPE!="snp"', 'TYPE!~"snp"', 'TYPE="indel"', 'TYPE~"indel"', 'TYPE!~"indel"',
                 'TYPE="mnp"', 'TYPE~"mnp"', 'TYPE="ref"', 'TYPE="snp" || TYPE="indel"', 'TYPE~"snp" && TYPE~"indel"',
                 'FILTER="q10"', 'FILTER~"q10"', 'FILTER!~"q10"', 'FILTER~"q10;s50"', 'FILTER="PASS"']:
        mmg.vcf_filter(tvcf, p("typed_out.vcf"), expression=expr)
        mine = [l for l in lines_of(p("typed_out.vcf")) if not l.startswith("#")]
        theirs = run("view", "-H", "-i", expr, tvcf)
        check(f"vcf_filter {expr} equals bcftools view -i", mine == theirs, f"{len(mine)} against {len(theirs)}")

    # Left-alignment: vcf_normalize(reference=...) against bcftools norm -f -m-,
    # on indels written anywhere inside the repeats of a random genome.
    def repeats(n):
        out = []
        while len(out) < n:
            unit = base(random.choice([1, 1, 2, 3]))
            out.extend(unit * random.randrange(2, 9) if random.random() < 0.6 else base(random.randrange(1, 8)))
        return "".join(out[:n])

    genome = {"chr1": repeats(300_000), "chr2": repeats(100_000)}
    with open(p("ref.fa"), "w", newline="\n") as f:
        for name, seq in genome.items():
            f.write(f">{name}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    indels = []
    for name, seq in genome.items():
        for pos0 in sorted(random.sample(range(500, len(seq) - 500), 3000)):
            a = random.randrange(0, 3)          # bases of context on the left
            ref = seq[pos0 - a:pos0 + 1]
            alts = []
            for _ in range(random.choice([1, 1, 2])):
                l = random.randrange(1, 6)
                if random.random() < 0.5:        # deletion of the next l bases
                    alts.append((seq[pos0 - a:pos0 + 1], l))
                else:                            # insertion: a copy of what follows, or random
                    alts.append((seq[pos0 - a:pos0 + 1] + (seq[pos0 + 1:pos0 + 1 + l] if random.random() < 0.7 else base(l)), 0))
            longest = max(l for _, l in alts)
            ref = seq[pos0 - a:pos0 + 1 + longest]
            alt = ",".join(dict.fromkeys(s + seq[pos0 + 1 + l:pos0 + 1 + longest] for s, l in alts))
            if alt != ref and ref not in alt.split(","):
                indels.append(f"{name}\t{pos0 - a + 1}\tv{len(indels)}\t{ref}\t{alt}\t30\tPASS\t.")
    ivcf = write("indels.vcf", ["##fileformat=VCFv4.2"] + [f"##contig=<ID={n},length={len(s)}>" for n, s in genome.items()]
                 + ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"] + indels)
    r = mmg.vcf_normalize(ivcf, p("ours_norm.vcf"), reference=p("ref.fa"))
    key = lambda l: tuple(l.split("\t")[:5])
    mine = sorted(key(l) for l in lines_of(p("ours_norm.vcf")) if not l.startswith("#"))
    theirs = sorted(key(l) for l in run("norm", "-f", p("ref.fa"), "-m-", ivcf))
    check("vcf_normalize(reference=...) equals bcftools norm -f -m- (CHROM, POS, ID, REF, ALT)",
          mine == theirs and r["ref_mismatch"] == 0 and r["left_aligned"] > 1000,
          f"{sum(1 for a, b in zip(mine, theirs) if a != b)} differ; {r}")

shutil.rmtree(T, ignore_errors=True)
print(f"\n{done - len(failures)} of {done} cross-checks passed" + (f"; skipped: {', '.join(sorted(set(skipped)))}" if skipped else ""))
if failures:
    print("FAILED:", "; ".join(failures))
    sys.exit(1)
