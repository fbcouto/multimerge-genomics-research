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
