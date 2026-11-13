# Test of the marimo interface without a browser: runs the cells of each screen
# as functions, with stand-in forms -- screen 8 (VCF) in detail and the k-mer
# flow of screens 1, 3, 4, 5, 6 and 7. Screen 2 (clustering) and the navigation
# are covered by test_browser.py, which drives a real browser.
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

print()
print(f"{done - len(failures)} of {done} checks passed")
if failures:
    print("FAILED:", failures)
    sys.exit(1)
