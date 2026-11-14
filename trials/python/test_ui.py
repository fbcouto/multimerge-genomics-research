# Test of the marimo interface without a browser: runs the cells of each screen
# as functions, with stand-in forms -- screen 8 (VCF) in detail, the k-mer flow
# of screens 1, 3, 4, 5, 6 and 7, and every task of the file tools (screens
# 9-14), whose outputs are compared with the direct call to the API. Screen 2
# (clustering) and the navigation are covered by test_browser.py, which drives a
# real browser.
# Usage: python trials/python/test_ui.py ui_main.py
# Needs marimo, polars, numpy and the extension installed.
import asyncio
import glob
import html as html_lib
import importlib.util
import inspect
import os
import random
import re
import sys
import tempfile

import marimo as mo
import numpy as np
import polars as pl

import multimerge_genomics as mmg

failures = []
done = 0


def check(name, cond, detail=""):
    global done
    done += 1
    if not cond:
        failures.append(name)
    print(("  ok    " if cond else "  FAIL  ") + name + ("" if cond else f"  {detail}"))


ui_path = os.path.abspath(sys.argv[1])
spec = importlib.util.spec_from_file_location("ui_main", ui_path)
ui = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ui)

# Work in an empty folder: the screens list the files of the current directory.
T = tempfile.mkdtemp(prefix="mmg_interface_")
os.chdir(T)
random.seed(8)
header = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]


def make_records(n, tag):
    return [f"{random.choice(['chr1', 'chr2', 'chr10'])}\t{random.randrange(1, 100000)}\t{tag}{i}\tA\tT\t{random.randrange(1, 99)}\tPASS\t." for i in range(n)]


def write_vcf(name, recs):
    with open(name, "w", newline="\n") as f:
        f.write("\n".join(header + recs) + "\n")


def read_vcf(name):
    with open(name, newline="\n") as f:
        return [l for l in f.read().split("\n") if l and not l.startswith("#")]


ra, rb, rc = make_records(3000, "a"), make_records(2000, "b"), make_records(1000, "c")
write_vcf("sample_a.vcf", ra)
write_vcf("sample_b.vcf", rb)
write_vcf("sample_c.vcf", rc)
open("table.parquet", "wb").close()


def run(cell, **refs):
    """Runs a cell; returns (output, definitions) or the exception raised."""
    try:
        r = cell.run(**refs)
        if inspect.isawaitable(r):
            async def wait():
                return await r
            r = asyncio.run(wait())
        return r
    except BaseException as e:  # mo.stop raises marimo's own exception
        if e.__cause__ is not None:  # marimo wraps the cell's error: show the inner one
            e.args = (repr(e.__cause__),)
        return e


def text_of(output):
    """The HTML of a marimo output, with the escapes undone."""
    t = html_lib.unescape(getattr(output, "text", None) or str(output))
    # The components keep their content as JSON inside an attribute.
    t = re.sub(r"\\u([0-9a-fA-F]{4})", lambda m: chr(int(m.group(1), 16)), t)
    return t.replace(r"\n", "\n")


class Form:
    def __init__(self, **value):
        self.value = value


print("== the whole notebook runs (home page)")
r = run(ui.app)
check("app.run() without error", not isinstance(r, BaseException), repr(r)[:300])

print("== screen 8: form")
r = run(ui.vcf_page_ui, get_page=lambda: "vcf", glob=glob, mo=mo, set_page=lambda _: None)
check("the form cell runs", not isinstance(r, BaseException), repr(r)[:300])
if not isinstance(r, BaseException):
    html = text_of(r[0])
    check("lists the .vcf files of the folder", all(n in html for n in ("sample_a.vcf", "sample_b.vcf", "sample_c.vcf")), html[:200])
    check("no longer offers .parquet", "table.parquet" not in html)


def backend(**value):
    return run(ui.vcf_page_backend, asyncio=asyncio, form_vcf=Form(**value), get_page=lambda: "vcf", mmg=mmg, mo=mo, os=os)


print("== screen 8: three files at once (the old defect kept only the last one)")
r = backend(files=["sample_a.vcf", "sample_b.vcf", "sample_c.vcf"], sort_by="POS", output_name="Merged.vcf")
check("the cell runs", not isinstance(r, BaseException), repr(r)[:300])
ref = sorted(ra + rb + rc, key=lambda l: (l.split("\t")[0], int(l.split("\t")[1])))
out = read_vcf("Merged.vcf") if os.path.exists("Merged.vcf") else []
check("the output has the 6000 variants of the THREE files", len(out) == 6000, f"{len(out)}")
check("identical to Python's stable sort", out == ref)
if not isinstance(r, BaseException):
    html = text_of(r[0])
    check("the panel shows the result", "6,000" in html and "Merged.vcf" in html, html[:300])

print("== screen 8: by QUAL")
r = backend(files=["sample_a.vcf", "sample_b.vcf"], sort_by="QUAL", output_name="ByQual.vcf")
quals = [int(l.split("\t")[5]) for l in read_vcf("ByQual.vcf")] if os.path.exists("ByQual.vcf") else []
check("5000 variants, QUAL descending", len(quals) == 5000 and quals == sorted(quals, reverse=True), f"{len(quals)}")

print("== screen 8: guards")
before = read_vcf("sample_a.vcf")
r = backend(files=["sample_a.vcf", "sample_b.vcf"], sort_by="POS", output_name="sample_a.vcf")
check("output equal to an input: not overwritten", read_vcf("sample_a.vcf") == before)
check("output equal to an input: the cell stops with a warning",
      isinstance(r, BaseException) or "cannot be one of the input files" in text_of(r[0]), repr(r)[:200])
r = backend(files=[], sort_by="POS", output_name="x.vcf")
check("no file: nothing is written", not os.path.exists("x.vcf"))
r = backend(files=["does_not_exist.vcf"], sort_by="POS", output_name="y.vcf")
check("missing file: error shown in the panel, nothing crashes",
      not isinstance(r, BaseException) and "ERROR" in text_of(r[0]), repr(r)[:300])

# ---------------------------------------------------------------------------
# The other screens: the k-mer flow from end to end, from FASTA to scaffold.
# ---------------------------------------------------------------------------
random.seed(1)
K = 55
genome = "".join(random.choice("ACGT") for _ in range(30000))
patient = list(genome)
for pos in (5000, 15000, 25000):
    patient[pos] = {"A": "C", "C": "G", "G": "T", "T": "A"}[patient[pos]]
patient = "".join(patient)


def write_fasta(name, ident, seq):
    with open(name, "w", newline="\n") as f:
        f.write(f">{ident}\n")
        for i in range(0, len(seq), 70):
            f.write(seq[i:i + 70] + "\n")


write_fasta("reference.fasta", "reference", genome)
write_fasta("patient.fasta", "patient", patient)


def number(html, label):
    # The first number after the label, in the text without the HTML tags.
    m = re.search(re.escape(label) + r"\D*?([\d,]+)", re.sub(r"<[^>]+>", "", html))
    return int(m.group(1).replace(",", "")) if m else None


print("== screen 1: QC / ingestion (FASTA -> k-mers -> .arrow)")
for name in ("reference", "patient"):
    r = run(ui.qc_page_backend, asyncio=asyncio, get_page=lambda: "qc", mmg=mmg, mo=mo, pl=pl,
            form_qc=Form(files=[f"{name}.fasta"], kmer_size=K, min_quality=20))
    ok = not isinstance(r, BaseException) and "successfully" in text_of(r[0]) and os.path.exists(f"{name}_K{K}.arrow")
    check(f"{name}.fasta -> {name}_K{K}.arrow", ok, repr(r)[:300] if isinstance(r, BaseException) else text_of(r[0])[:300])
    if ok:
        n = number(text_of(r[0]), "extracted:")
        check(f"  {n} k-mers (expected {30000 - K + 1})", n == 30000 - K + 1, str(n))

print("== screen 3: assembly (De Bruijn)")
r = run(ui.assembly_page_backend, asyncio=asyncio, get_page=lambda: "assembly", mmg=mmg, mo=mo, os=os, pl=pl,
        form_assembly=Form(kmer_file=f"reference_K{K}.arrow"))
ok = not isinstance(r, BaseException) and "successfully" in text_of(r[0]) and os.path.exists(f"MMGEN_Contigs_K{K}.fasta")
check("assembles and writes the contigs", ok, repr(r)[:300] if isinstance(r, BaseException) else text_of(r[0])[:300])
if ok:
    contigs = [l.strip() for l in open(f"MMGEN_Contigs_K{K}.fasta") if not l.startswith(">")]
    longest = max(contigs, key=len)
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    revcomp = "".join(comp[b] for b in reversed(longest))
    check(f"  the longest contig ({len(longest)} bp) is a stretch of the original genome", longest in genome or revcomp in genome, longest[:60])
    check("  and it rebuilds the whole genome (a sequence with no repeats)", len(longest) == len(genome), f"{len(longest)} of {len(genome)}; {len(contigs)} contigs")

print("== screen 4: variant calling (Bloom)")
r = run(ui.variant_page_backend, asyncio=asyncio, get_page=lambda: "variant", mmg=mmg, mo=mo, pl=pl,
        form_variant=Form(reference=f"reference_K{K}.arrow", patient=f"patient_K{K}.arrow", filter_mb=16, num_hashes=5))
ok = not isinstance(r, BaseException) and "scan complete" in text_of(r[0])
check("the scan runs", ok, repr(r)[:300] if isinstance(r, BaseException) else text_of(r[0])[:300])
if ok:
    n = number(text_of(r[0]), "patient):")
    check(f"  {n} mutated k-mers (3 SNPs x {K} = {3 * K})", n == 3 * K, str(n))
r = run(ui.variant_page_backend, asyncio=asyncio, get_page=lambda: "variant", mmg=mmg, mo=mo, pl=pl,
        form_variant=Form(reference=f"reference_K{K}.arrow", patient=f"reference_K{K}.arrow", filter_mb=16, num_hashes=5))
check("  reference == patient: stops with a warning (this is where kind=\"warning\" used to break)",
      not isinstance(r, BaseException) or type(r).__name__ == "MarimoStopError", repr(r)[:200])

print("== screen 5: annotation (translation)")
r = run(ui.annotation_page_backend, asyncio=asyncio, get_page=lambda: "annotation", mmg=mmg, mo=mo,
        form_annotation=Form(dna="ATGGCCAAATAG"))
check("ATGGCCAAATAG -> MAK*", not isinstance(r, BaseException) and "MAK" in text_of(r[0]) and "4 amino acids" in text_of(r[0]),
      repr(r)[:300] if isinstance(r, BaseException) else text_of(r[0])[:400])

print("== screen 6: Bloom filter")
r = run(ui.bloom_page_backend, asyncio=asyncio, get_page=lambda: "bloom", mmg=mmg, mo=mo, np=np, pl=pl,
        form_bloom=Form(kmer_file=f"reference_K{K}.arrow", filter_mb=16, num_hashes=5))
ok = not isinstance(r, BaseException) and "successfully" in text_of(r[0])
check("builds the filter", ok, repr(r)[:300] if isinstance(r, BaseException) else text_of(r[0])[:300])
if ok:
    bits = number(text_of(r[0]), "Bits set:")
    expected = 5 * (30000 - K + 1)
    check(f"  {bits} bits set (5 per k-mer: close to {expected})", bits is not None and 0.99 * expected <= bits <= expected, str(bits))

print("== screen 7: scaffolding")
r = run(ui.scaffold_page_backend, asyncio=asyncio, get_page=lambda: "scaffold", mmg=mmg, mo=mo,
        form_scaffold=Form(contigs_fasta=f"MMGEN_Contigs_K{K}.fasta", overlap=K - 1))
check("stitches the contigs", not isinstance(r, BaseException) and "Scaffolding complete" in text_of(r[0]),
      repr(r)[:300] if isinstance(r, BaseException) else text_of(r[0])[:300])

# ---------------------------------------------------------------------------
# Screens 9-14, the file tools. Each screen writes its result with the same
# call a script would make: the test makes that call itself and compares the
# two files, so a wrong argument between the form and the engine shows up.
# ---------------------------------------------------------------------------
import builtins
import gzip

helpers = run(ui.file_helpers, mmg=mmg, mo=mo, os=os)
check("the helpers cell runs", not isinstance(helpers, BaseException), repr(helpers)[:300])
POOL = dict(helpers[1], asyncio=asyncio, glob=glob, mmg=mmg, mo=mo, np=np, os=os, pl=pl)


class Choice:
    """Stands in for a radio, a dropdown or a button: only `.value` is read."""
    def __init__(self, value=None):
        self.value = value


def cell_run(cell, **given):
    pool = dict(POOL, **given)
    return run(cell, **{k: pool[k] for k in cell.refs if not hasattr(builtins, k)})


def raw(r):
    """The HTML the cell put on the page (a mo.stop comes back as the output too)."""
    if isinstance(r, BaseException):
        return repr(r)
    return text_of(r[0])


def shown(r):
    """The text of the page, without the tags (the code panel is highlighted
    token by token, so a call only reads as one string once they are gone)."""
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", raw(r)))


def stopped(r, message):
    """The cell stopped with this warning (and so wrote nothing)."""
    return not isinstance(r, BaseException) and message in shown(r) and "ERROR" not in shown(r)


def body(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", newline="\n") as f:
        return [l for l in f.read().split("\n") if l and not l.startswith("#")]


def same_output(name, screen_file, api_file):
    ok = os.path.exists(screen_file) and os.path.exists(api_file) and body(screen_file) == body(api_file)
    check(name, ok, f"{screen_file}: {os.path.exists(screen_file)}")


random.seed(14)
CHROMS = {"chr1": 300_000, "chr2": 200_000, "chr10": 100_000}
with open("genome.fai", "w", newline="\n") as f:
    for c, n in CHROMS.items():
        f.write(f"{c}\t{n}\t0\t60\t61\n")


def write_bed(name, n):
    with open(name, "w", newline="\n") as f:
        for _ in range(n):
            c = random.choice(list(CHROMS))
            s = random.randrange(0, CHROMS[c] - 2000)
            f.write(f"{c}\t{s}\t{s + random.randrange(50, 1500)}\tname{_}\n")


write_bed("peaks.bed", 400)
write_bed("genes.bed", 250)

print("== screens 9-14: every task builds its form")
for cell, page, task_name, tasks in (
    (ui.sort_page_form, "sort", "sort_task", ["Sort & merge", "Index & query a region"]),
    (ui.validate_page_form, "validate", "validate_task", ["Validate a file", "Rename chromosomes (chr1 ↔ 1)", "Compare reference genomes"]),
    (ui.vcftools_page_form, "vcftools", "vcftools_task", ["Statistics", "Filter", "Normalize"]),
    (ui.fastq_page_form, "fastq", "fastq_task", ["Quality report", "Trim reads"]),
):
    for t in tasks:
        r = cell_run(cell, get_page=lambda p=page: p, **{task_name: Choice(t)}, interval_op=Choice("merge"),
                     **{k: None for k in cell.refs if k.startswith("btn_refresh")})
        check(f"{page}: '{t}'", not isinstance(r, BaseException), repr(r)[:300])
for op in ("merge", "intersect", "subtract", "complement", "closest", "coverage", "slop"):
    r = cell_run(ui.sort_page_form, get_page=lambda: "sort", sort_task=Choice("Interval operations (BED)"),
                 interval_op=Choice(op), btn_refresh_sort=None)
    check(f"sort: interval operation '{op}'", not isinstance(r, BaseException) and "peaks.bed" in shown(r), repr(r)[:300])
for cell, page in ((ui.kmercount_page_form, "kmercount"), (ui.convert_page_form, "convert")):
    r = cell_run(cell, get_page=lambda p=page: p, **{k: None for k in cell.refs if k.startswith("btn_refresh")})
    check(f"{page}: the form", not isinstance(r, BaseException), repr(r)[:300])
r = cell_run(ui.sort_page_form, get_page=lambda: "sort", sort_task=Choice("Sort & merge"), interval_op=Choice("merge"), btn_refresh_sort=None)
check("sort: offers the VCF and BED files of the folder and the .fai as a chromosome order",
      all(n in shown(r) for n in ("sample_a.vcf", "peaks.bed", "the order of genome.fai")), shown(r)[:300])


def sort_backend(task, op="merge", **value):
    return cell_run(ui.sort_page_backend, get_page=lambda: "sort", sort_task=Choice(task), interval_op=Choice(op), form_sort=Form(**value))


print("== screen 9: sort & merge, with the index")
r = sort_backend("Sort & merge", files=["sample_a.vcf", "sample_b.vcf", "sample_c.vcf"], by="position", order="natural",
                 memory=1, index=True, output="")
check("writes sample_a.sorted.vcf.gz and its .tbi", os.path.exists("sample_a.sorted.vcf.gz") and os.path.exists("sample_a.sorted.vcf.gz.tbi"), shown(r)[-300:])
mmg.sort_genomic(["sample_a.vcf", "sample_b.vcf", "sample_c.vcf"], "api_sorted.vcf.gz", chrom_order="natural", memory_mb=1)
same_output("  the same records as sort_genomic(chrom_order='natural')", "sample_a.sorted.vcf.gz", "api_sorted.vcf.gz")
check("  the panel shows 6,000 records and the Python call", "6,000" in shown(r) and "mmg.sort_genomic(" in shown(r) and "index=True" in shown(r), shown(r)[-400:])
r = sort_backend("Sort & merge", files=["sample_a.vcf"], by="position", order="default", memory=0, index=False, output="sample_a.vcf")
check("  output equal to an input: stops with a warning", stopped(r, "cannot be one of the input files"), shown(r)[:200])

print("== screen 9: index & query a region")
expected = [l for l in body("api_sorted.vcf.gz") if l.split("\t")[0] == "chr2" and 1000 <= int(l.split("\t")[1]) <= 30000]
r = sort_backend("Index & query a region", file="sample_a.sorted.vcf.gz", regions="chr2:1000-30000", output="")
check(f"  {len(expected)} records of chr2:1000-30000, shown on the page", f"{len(expected):,} record(s) in chr2:1000-30000" in shown(r), shown(r)[-300:])
check("  the existing index is reused", "Using the index" in shown(r), shown(r)[:300])
r = sort_backend("Index & query a region", file="sample_a.sorted.vcf.gz", regions="chr2:1000-30000", output="region.vcf")
check("  saved to region.vcf: the same records", os.path.exists("region.vcf") and body("region.vcf") == expected, shown(r)[-300:])
os.remove("sample_a.sorted.vcf.gz.tbi")
r = sort_backend("Index & query a region", file="sample_a.sorted.vcf.gz", regions="", output="")
check("  without the index: builds it", "Index written" in shown(r) and os.path.exists("sample_a.sorted.vcf.gz.tbi"), shown(r)[:300])

print("== screen 9: the seven interval operations, against the direct call")
api = {
    "merge": lambda o: mmg.bed_merge(["peaks.bed"], o, distance=100),
    "intersect": lambda o: mmg.bed_intersect(["peaks.bed"], ["genes.bed"], o, mode="unique"),
    "subtract": lambda o: mmg.bed_subtract(["peaks.bed"], ["genes.bed"], o, remove_whole=True),
    "complement": lambda o: mmg.bed_complement(["peaks.bed"], "genome.fai", o),
    "closest": lambda o: mmg.bed_closest(["peaks.bed"], ["genes.bed"], o),
    "coverage": lambda o: mmg.bed_coverage(["peaks.bed"], ["genes.bed"], o),
    "slop": lambda o: mmg.bed_slop(["peaks.bed"], "genome.fai", o, left=50, right=300),
}
form_values = {
    "merge": dict(distance=100), "intersect": dict(b=["genes.bed"], mode="unique"),
    "subtract": dict(b=["genes.bed"], whole=True), "complement": dict(genome="genome.fai"),
    "closest": dict(b=["genes.bed"]), "coverage": dict(b=["genes.bed"]), "slop": dict(genome="genome.fai", left=50, right=300),
}
for op, call in api.items():
    r = sort_backend("Interval operations (BED)", op, a=["peaks.bed"], output="", **form_values[op])
    call(f"api_{op}.bed")
    same_output(f"  {op}: peaks.{op}.bed is what the direct call writes", f"peaks.{op}.bed", f"api_{op}.bed")
r = sort_backend("Interval operations (BED)", "intersect", a=["peaks.bed"], b=[], mode="regions", output="")
check("  intersect without B: stops with a warning", stopped(r, "Select at least one B file"), shown(r)[:200])


def validate_backend(task, **value):
    return cell_run(ui.validate_page_backend, get_page=lambda: "validate", validate_task=Choice(task), form_validate=Form(**value))


print("== screen 10: validation, names, references")
contigs = "".join(f"##contig=<ID={c},length={n}>\n" for c, n in CHROMS.items())
with open("clean.vcf", "w", newline="\n") as f:
    f.write("##fileformat=VCFv4.2\n" + contigs + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for i, p in enumerate(range(100, 5000, 100)):
        f.write(f"chr1\t{p}\tv{i}\tA\tG\t50\tPASS\t.\n")
with open("broken.vcf", "w", newline="\n") as f:
    f.write("##fileformat=VCFv4.2\n" + contigs + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    f.write("chr1\t500\ta\tA\tG\t50\tPASS\t.\nchr1\t100\tb\tA\tG\t50\tPASS\t.\n1\t200\tc\tA\tG\t50\tPASS\t.\n")
r = validate_backend("Validate a file", file="clean.vcf", max_lines=0)
check("  a clean VCF: no problem found", "No problem found" in shown(r), shown(r)[:300])
r = validate_backend("Validate a file", file="broken.vcf", max_lines=0)
check("  a broken VCF: errors, with the line", "error(s)" in shown(r) and re.search(r"line \d", shown(r)) is not None, shown(r)[:600])
check("  and the next steps name the screens", "Rename chromosomes" in shown(r) and "9. Sorting" in shown(r), shown(r)[-500:])
r = validate_backend("Rename chromosomes (chr1 ↔ 1)", file="clean.vcf", style="no_chr", mapping="none", output="")
check("  rename to 1, X, MT: clean.nochr.vcf", os.path.exists("clean.nochr.vcf") and all(l.startswith("1\t") for l in body("clean.nochr.vcf")), shown(r)[:300])
with open("hg38.fai", "w", newline="\n") as f:
    f.write("chr1\t248956422\t0\t0\t0\nchr2\t242193529\t0\t0\t0\n")
with open("hg19.fai", "w", newline="\n") as f:
    f.write("1\t249250621\t0\t0\t0\n2\t243199373\t0\t0\t0\n")
r = validate_backend("Compare reference genomes", a="hg38.fai", b="hg19.fai")
check("  GRCh38 against GRCh37: DIFFERENT, with both assemblies named",
      "DIFFERENT reference genomes" in shown(r) and "GRCh38" in shown(r) and "GRCh37" in shown(r), shown(r)[:500])


def vcftools_backend(task, **value):
    return cell_run(ui.vcftools_page_backend, get_page=lambda: "vcftools", vcftools_task=Choice(task), form_vcftools=Form(**value))


print("== screen 11: VCF statistics, filter, normalization")
r = vcftools_backend("Statistics", file="sample_a.vcf")
check("  statistics: 3,000 records, 3,000 SNPs", "3,000" in shown(r) and "SNP 3,000" in shown(r), shown(r)[:400])
r = vcftools_backend("Filter", file="sample_a.vcf", expression="QUAL >= 50", regions="chr1, chr10", samples="", output="")
mmg.vcf_filter("sample_a.vcf", "api_filtered.vcf", expression="QUAL >= 50", regions=["chr1", "chr10"])
same_output("  filter: sample_a.filtered.vcf is what vcf_filter writes", "sample_a.filtered.vcf", "api_filtered.vcf")
r = vcftools_backend("Filter", file="sample_a.vcf", expression="QUAL >= ", regions="", samples="", output="bad.vcf")
check("  a broken expression: the error is shown on the page", not isinstance(r, BaseException) and "ERROR" in shown(r), shown(r)[:300])
# A file in the wrong format (5 columns: no REF, ALT, FILTER, INFO) used to give
# a green panel full of zeros and "Lines skipped: 5", with no reason.
with open("five_columns.vcf", "w", newline="\n") as f:
    f.write("##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tQUAL\tGT\n1\t100\trs1\t10.5\t0/0\n1\t200\trs2\t50.2\t0/1\n")
for task, extra in (("Statistics", {}), ("Filter", dict(expression="QUAL > 1", regions="", samples="", output=""))):
    r = vcftools_backend(task, file="five_columns.vcf", **extra)
    check(f"  {task.lower()} of a 5-column file: says nothing was read, and the validator's reason",
          "Nothing was read" in shown(r) and "needs at least 8" in shown(r), shown(r)[:500])
r = sort_backend("Sort & merge", files=["five_columns.vcf"], by="position", order="default", memory=0, index=False, output="")
check("  sort of a 5-column file: the same explanation", "Nothing to sort" in shown(r) and "needs at least 8" in shown(r), shown(r)[:500])
r = vcftools_backend("Statistics", file="sample_a.vcf")
check("  a valid file: no such panel", "skipped" not in shown(r).replace("Lines skipped (not valid records): 0", ""), shown(r)[:500])

random.seed(5)
left = "".join(random.choice("ACG") for _ in range(100))
genome_seq = left + "C" + "TTTTT" + "G" + "".join(random.choice("ACGT") for _ in range(400))
with open("small_genome.fa", "w", newline="\n") as f:
    f.write(">chr1\n" + "\n".join(genome_seq[i:i + 60] for i in range(0, len(genome_seq), 60)) + "\n")
with open("indel.vcf", "w", newline="\n") as f:
    f.write(f"##fileformat=VCFv4.2\n##contig=<ID=chr1,length={len(genome_seq)}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    f.write("chr1\t105\tdel\tTT\tT\t50\tPASS\t.\n")
r = vcftools_backend("Normalize", file="indel.vcf", split=True, reference="small_genome.fa", output="")
check("  normalize with the genome: the deletion moves to the left of the T run (POS 101, CT>C)",
      os.path.exists("indel.norm.vcf") and body("indel.norm.vcf") == ["chr1\t101\tdel\tCT\tC\t50\tPASS\t."], body("indel.norm.vcf") if os.path.exists("indel.norm.vcf") else shown(r)[:300])
check("  and the panel counts it", "Indels left-aligned: 1" in shown(r), shown(r)[:600])


def fastq_backend(task, **value):
    return cell_run(ui.fastq_page_backend, get_page=lambda: "fastq", fastq_task=Choice(task), form_fastq=Form(**value))


print("== screen 12: FASTQ report and trimming")
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
for mate in (1, 2):
    with gzip.open(f"pair_R{mate}.fastq.gz", "wt", newline="\n") as f:
        for i in range(2000):
            insert = "".join(random.choice("ACGT") for _ in range(random.randrange(40, 100)))
            seq = (insert + ADAPTER)[:100]
            f.write(f"@p{i}\n{seq}\n+\n{''.join(chr(33 + random.randrange(5, 41)) for _ in seq)}\n")
r = fastq_backend("Quality report", file="pair_R1.fastq.gz", max_reads=0, output="")
check("  report: pair_R1.report.html, the adapter found, the page embedded",
      os.path.exists("pair_R1.report.html") and "Illumina Universal" in shown(r) and "iframe" in raw(r).lower(), shown(r)[:400])
r = fastq_backend("Trim reads", file="pair_R1.fastq.gz", file2="pair_R2.fastq.gz", adapter="auto", custom="", quality=20,
                  min_length=20, trim_n=True, output="", output2="")
mmg.fastq_trim("pair_R1.fastq.gz", "api_R1.fastq.gz", input2="pair_R2.fastq.gz", output2="api_R2.fastq.gz", adapter="auto",
               quality_cutoff=20, trim_n=True, min_length=20)
same_output("  trim, paired: R1 is what fastq_trim writes", "pair_R1.trimmed.fastq.gz", "api_R1.fastq.gz")
same_output("  trim, paired: R2 is what fastq_trim writes", "pair_R2.trimmed.fastq.gz", "api_R2.fastq.gz")
check("  the panel says pairs and the adapter used", "Pairs in / kept" in shown(r) and "AGATCGGAAGAGC" in shown(r), shown(r)[:500])
r = fastq_backend("Trim reads", file="pair_R1.fastq.gz", file2="none", adapter="custom", custom="ACGTX", quality=0,
                  min_length=0, trim_n=False, output="x.fastq", output2="")
check("  a custom adapter with a wrong letter: stops with a warning", stopped(r, "A, C, G, T and N only"), shown(r)[:200])

print("== screen 13: k-mer counting")
r = cell_run(ui.kmercount_page_backend, get_page=lambda: "kmercount",
             form_kmercount=Form(files=["pair_R1.fastq.gz", "pair_R2.fastq.gz"], k=21, canonical=True, min_count=2, output="kmers.tsv.gz"))
api_counts = mmg.count_kmers(["pair_R1.fastq.gz", "pair_R2.fastq.gz"], "api_kmers.tsv.gz", k=21, min_count=2)
same_output("  the table is what count_kmers writes", "kmers.tsv.gz", "api_kmers.tsv.gz")
check("  the panel shows the totals, the spectrum chart and its table",
      f"{api_counts['distinct_kmers']:,}" in shown(r) and "<svg" in raw(r) and "Abundance spectrum" in shown(r) and "times seen" in shown(r), shown(r)[:400])

print("== screen 14: formats & conversion")


def convert_backend(**value):
    return cell_run(ui.convert_page_backend, get_page=lambda: "convert", form_convert=Form(**value))


with open("quals.vcf", "w", newline="\n") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    f.write("chr1\t10\ta\tA\tG\t12.50\tPASS\t.\nchr1\t20\tb\tC\tT\t30.0\tPASS\t.\n")
r = convert_backend(file="quals.vcf", target=".parquet", numbers="text", output="")
check("  VCF -> quals.vcf.parquet, with QUAL as text", os.path.exists("quals.vcf.parquet") and "parquet_options(numbers=\"text\")" in shown(r), shown(r)[:400])
check("  the process setting is put back afterwards", mmg.parquet_options()["numbers"] == "number", mmg.parquet_options())
r = convert_backend(file="quals.vcf.parquet", target="text", numbers="number", output="quals_back.vcf")
check("  and back to text: the same bytes (12.50, 30.0)", open("quals_back.vcf", "rb").read() == open("quals.vcf", "rb").read(), body("quals_back.vcf"))
with open("reads.sam", "w", newline="\n") as f:
    f.write("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n")
    for i in range(20):
        f.write(f"q{i}\t0\tchr1\t{10 + i * 30}\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n")
r = convert_backend(file="reads.sam", target=".bam", numbers="number", output="")
check("  SAM -> reads.bam, 20 records", os.path.exists("reads.bam") and "22 (20 records)" in shown(r), shown(r)[:300])
r = convert_backend(file="quals.vcf", target="text", numbers="number", output="")
check("  already in that container: stops with a warning", stopped(r, "already in that container"), shown(r)[:200])

print()
print(f"{done - len(failures)} of {done} checks passed")
if failures:
    print("FAILED:", failures)
    sys.exit(1)
