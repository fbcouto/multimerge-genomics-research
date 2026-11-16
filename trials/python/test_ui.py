# Test of the marimo interface without a browser: runs the cells of each screen
# as functions, with stand-in forms -- screen 8 (VCF) in detail, the k-mer flow
# of screens 1, 3, 4, 5, 6 and 7, every task of the file tools (screens 9-14),
# whose outputs are compared with the direct call to the API, and every example
# of every screen (picked, the form filled in, run). Screen 2 (clustering) and
# the navigation are covered by test_browser.py, which drives a real browser.
# Usage: python trials/python/test_ui.py [python/multimerge_genomics/ui_main.py]
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


# The interface comes with the package: python/multimerge_genomics/ui_main.py (the default).
ui_path = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "python", "multimerge_genomics", "ui_main.py"))
spec = importlib.util.spec_from_file_location("ui_main", ui_path)
ui = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ui)
# The figures from the same tree as the interface (the installed package may be older).
_fspec = importlib.util.spec_from_file_location("mmgen_figures", os.path.join(os.path.dirname(ui_path), "figures.py"))
figs = importlib.util.module_from_spec(_fspec)
_fspec.loader.exec_module(figs)

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


class Choice:
    """Stands in for a radio, a dropdown or a button: only `.value` is read."""
    def __init__(self, value=None):
        self.value = value


print("== the whole notebook runs (home page)")
r = run(ui.app)
check("app.run() without error", not isinstance(r, BaseException), repr(r)[:300])

import builtins
import gzip

helpers = run(ui.file_helpers, mmg=mmg, mo=mo, os=os, re=re)
check("the helpers cell runs", not isinstance(helpers, BaseException), repr(helpers)[:300])
examples = run(ui.example_files, EXAMPLES_DIR=helpers[1]["EXAMPLES_DIR"], mmg=mmg, mo=mo, os=os, pl=pl)
check("the examples cell runs", not isinstance(examples, BaseException), repr(examples)[:300])
POOL = dict(helpers[1], **examples[1], asyncio=asyncio, figs=figs, glob=glob, mmg=mmg, mo=mo, np=np, os=os, pl=pl, re=re)

# The AI (screen 29 and the "Ask the AI" panel) with a stand-in for Claude:
# no key, nothing leaves the machine. It answers by the kind of request.
from types import SimpleNamespace as NS

from multimerge_genomics import ai as ai_lib
from multimerge_genomics.ai import assistant as ai_assistant_mod

AI_REQUESTS = []
PLAN_ANSWER = "\n".join([
    "1. **The problem**: the two files name chromosomes differently.", "",
    "2. **The path**: screen 10, then screen 9.", "",
    "3. **Watch out**: the reference genome.", "",
    "4. **The script**:",
    "```python",
    "import multimerge_genomics as mmg",
    "mmg.normalize_chromosomes('ensembl_names.vcf', 'renamed.vcf', style='chr')",
    "mmg.sort_genomic(['renamed.vcf'], 'renamed.sorted.vcf.gz', index=True)",
    "```", ""])


def ai_message(text="", tools=(), stop="end_turn"):
    content = [NS(type="text", text=text)] if text else []
    content += [NS(type="tool_use", id=f"tu_{i}", name=n, input=inp) for i, (n, inp) in enumerate(tools)]
    return NS(content=content, stop_reason=stop, stop_details=None, model="claude-opus-5-5",
              usage=NS(input_tokens=1200, output_tokens=150, cache_read_input_tokens=20000, cache_creation_input_tokens=0))


def ai_reply(kw):
    last = kw["messages"][-1]["content"]
    if "tools" in kw:  # the agent: one round of calls, then the report
        if isinstance(last, str):
            return ai_message("First I validate the file.", [("validate_file", {"path": "with_problems.vcf"})], stop="tool_use")
        return ai_message("Report: the validator found the planted problems.")
    if "My problem" in last:
        return ai_message(PLAN_ANSWER)
    return ai_message("The tool read the file; the numbers look as expected for a small simulated file.")


class AIStream:
    def __init__(self, msg):
        self.msg = msg

    async def __aenter__(self):
        return self

    async def __aexit__(self, *a):
        return False

    @property
    def text_stream(self):
        async def gen():
            text = "".join(b.text for b in self.msg.content if b.type == "text")
            for i in range(0, len(text), 9):
                yield text[i:i + 9]
        return gen()

    async def get_final_message(self):
        return self.msg


class AIClient:
    def __init__(self):
        class Messages:
            def create(self, **kw):
                AI_REQUESTS.append(kw)
                return ai_reply(kw)

            def stream(self, **kw):
                AI_REQUESTS.append(kw)
                return AIStream(ai_reply(kw))
        self.beta = NS(messages=Messages())


ai_assistant_mod.client_factory = lambda key, is_async: AIClient()


class Holder:
    """Stands in for a mo.state: get() and set(value or function)."""
    def __init__(self, value=None):
        self.value = value

    def get(self):
        return self.value

    def set(self, v):
        self.value = v(self.value) if callable(v) else v


AI_SETTINGS = {"key": "", "model": "claude-opus-5-5", "effort": "medium", "level": "student", "consent": True, "human": False}
AGENT, AI_REV = Holder(), Holder(0)
ai_help = run(ui.ai_helpers, error_panel=POOL["error_panel"], mo=mo, os=os, warn=POOL["warn"])
check("the AI helpers cell runs", not isinstance(ai_help, BaseException) and ai_help[1]["ai_lib"] is not None, repr(ai_help)[:300])
POOL.update(ai_help[1], ai_settings=AI_SETTINGS, get_agent=AGENT.get, set_agent=AGENT.set, get_ai_rev=AI_REV.get, set_ai_rev=AI_REV.set)


def cell_run(cell, **given):
    """Runs a cell with what it refers to: the helpers, the modules and `given`.
    The examples lists default to "none picked", the refresh buttons to None,
    the data folders to "."."""
    pool = dict(POOL, **given)
    for k in cell.refs:
        if k.endswith("_example") and k not in pool:
            pool[k] = Choice("")
        if k.startswith("btn_refresh") and k not in pool:
            pool[k] = None
        if k.endswith("_folder") and k not in pool:
            pool[k] = Choice(".")
    pool.setdefault("PHARMCAT_DIR", "mmgen_tools/pharmcat")
    pool.setdefault("pgx_downloaded", 0)
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


print("== screen 8: form")
r = cell_run(ui.vcf_page_form, get_page=lambda: "vcf")
check("the form cell runs", not isinstance(r, BaseException), repr(r)[:300])
if not isinstance(r, BaseException):
    html = text_of(r[0])
    check("lists the .vcf files of the folder", all(n in html for n in ("sample_a.vcf", "sample_b.vcf", "sample_c.vcf")), html[:200])
    check("no longer offers .parquet", "table.parquet" not in html)
    # mo.ui.dictionary draws its fields inside a JSON tree (marimo-json-output).
    check("a clean form, not a JSON tree", "marimo-json-output" not in raw(r) and "marimo-multiselect" in raw(r), html[:300])


def backend(**value):
    return cell_run(ui.vcf_page_backend, get_page=lambda: "vcf", form_vcf=Form(**value))


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
    r = cell_run(ui.qc_page_backend, get_page=lambda: "qc", form_qc=Form(files=[f"{name}.fasta"], kmer_size=K, min_quality=20))
    ok = not isinstance(r, BaseException) and "successfully" in text_of(r[0]) and os.path.exists(f"{name}_K{K}.arrow")
    check(f"{name}.fasta -> {name}_K{K}.arrow", ok, repr(r)[:300] if isinstance(r, BaseException) else text_of(r[0])[:300])
    if ok:
        n = number(text_of(r[0]), "extracted:")
        check(f"  {n} k-mers (expected {30000 - K + 1})", n == 30000 - K + 1, str(n))

print("== screen 3: assembly (De Bruijn)")
r = cell_run(ui.assembly_page_backend, get_page=lambda: "assembly", form_assembly=Form(kmer_file=f"reference_K{K}.arrow"))
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
r = cell_run(ui.variant_page_backend, get_page=lambda: "variant",
             form_variant=Form(reference=f"reference_K{K}.arrow", patient=f"patient_K{K}.arrow", filter_mb=16, num_hashes=5))
ok = not isinstance(r, BaseException) and "scan complete" in text_of(r[0])
check("the scan runs", ok, repr(r)[:300] if isinstance(r, BaseException) else text_of(r[0])[:300])
if ok:
    n = number(text_of(r[0]), "patient):")
    check(f"  {n} mutated k-mers (3 SNPs x {K} = {3 * K})", n == 3 * K, str(n))
r = cell_run(ui.variant_page_backend, get_page=lambda: "variant",
             form_variant=Form(reference=f"reference_K{K}.arrow", patient=f"reference_K{K}.arrow", filter_mb=16, num_hashes=5))
check("  reference == patient: stops with a warning (this is where kind=\"warning\" used to break)",
      not isinstance(r, BaseException) or type(r).__name__ == "MarimoStopError", repr(r)[:200])

print("== screen 5: annotation (translation)")
r = cell_run(ui.annotation_page_backend, get_page=lambda: "annotation", form_annotation=Form(dna="ATGGCCAAATAG"))
check("ATGGCCAAATAG -> MAK*", not isinstance(r, BaseException) and "MAK" in text_of(r[0]) and "4 amino acids" in text_of(r[0]),
      repr(r)[:300] if isinstance(r, BaseException) else text_of(r[0])[:400])

print("== screen 6: Bloom filter")
r = cell_run(ui.bloom_page_backend, get_page=lambda: "bloom", form_bloom=Form(kmer_file=f"reference_K{K}.arrow", filter_mb=16, num_hashes=5))
ok = not isinstance(r, BaseException) and "successfully" in text_of(r[0])
check("builds the filter", ok, repr(r)[:300] if isinstance(r, BaseException) else text_of(r[0])[:300])
if ok:
    bits = number(text_of(r[0]), "Bits set:")
    expected = 5 * (30000 - K + 1)
    check(f"  {bits} bits set (5 per k-mer: close to {expected})", bits is not None and 0.99 * expected <= bits <= expected, str(bits))

print("== screen 7: scaffolding")
r = cell_run(ui.scaffold_page_backend, get_page=lambda: "scaffold", form_scaffold=Form(contigs_fasta=f"MMGEN_Contigs_K{K}.fasta", overlap=K - 1))
check("stitches the contigs", not isinstance(r, BaseException) and "Scaffolding complete" in text_of(r[0]),
      repr(r)[:300] if isinstance(r, BaseException) else text_of(r[0])[:300])

# ---------------------------------------------------------------------------
# Screens 9-14, the file tools. Each screen writes its result with the same
# call a script would make: the test makes that call itself and compares the
# two files, so a wrong argument between the form and the engine shows up.
# ---------------------------------------------------------------------------
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
      f"{api_counts['distinct_kmers']:,}" in shown(r) and "<svg" in raw(r) and "Spectrum (log scale)" in shown(r) and "Spectrum (linear)" in shown(r) and "times seen" in shown(r), shown(r)[:400])

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

print("== the rail of tools and the header strip")
r = cell_run(ui.workbench_rail, get_page=lambda: "sort", set_page=lambda _: None)
check("  the rail lists the 31 screens as links, the open one marked",
      not isinstance(r, BaseException) and raw(r).count('class="mmg-item') == 32 and 'class="mmg-item active" href="#sort"' in raw(r)
      and raw(r).count("mmg-rail-bridge") == 1, raw(r)[:300])
_rail = r[1]
_rail_pages = {p for _, _, items in _rail["RAIL_GROUPS"] for p, _, _ in items}
_menu_pages = {p for p, _ in examples[1]["EXAMPLES"]} | {"guide"}   # every screen has examples but the guide
check("  every screen is in the rail, and only those", _rail_pages == _menu_pages, sorted(_rail_pages ^ _menu_pages))
r = cell_run(ui.main_menu, RAIL_GROUPS=_rail["RAIL_GROUPS"], get_page=lambda: "home")
check("  the start page: the families with their screens as links, the folder, the engine",
      not isinstance(r, BaseException) and raw(r).count('class="mmg-card"') == 5 and raw(r).count("data-page=") == 31
      and "running on:" in shown(r) and os.path.basename(os.getcwd()) in shown(r), raw(r)[:300])
r = cell_run(ui.workbench_header, RAIL_GROUPS=_rail["RAIL_GROUPS"], RAIL_HUMAN=_rail["RAIL_HUMAN"], get_page=lambda: "pgx")
check("  the header strip: group, screen, folder and the human-genome chip",
      not isinstance(r, BaseException) and "Human genetics" in shown(r) and "19. Pharmacogenomics" in shown(r)
      and os.path.basename(os.getcwd()) in shown(r) and "human genome" in shown(r), shown(r)[:300])
for _side in ("ui_theme.css", "ui_head.html"):
    check(f"  {_side} is next to the interface", os.path.exists(os.path.join(os.path.dirname(ui.__file__), _side)))

print("== screen 15: the guide to the tests")
r = cell_run(ui.guide_page_ui, get_page=lambda: "guide", set_page=lambda _: None)
check("  the guide shows the table of tests, their blind spots and the notice",
      not isinstance(r, BaseException) and "Does not see" in shown(r) and "Optical genome mapping" in shown(r)
      and "for teaching and research only" in shown(r), shown(r)[:300])

# ---------------------------------------------------------------------------
# The examples: for every screen and task, each example fills in the form,
# and the form, submitted as it is, runs without a warning or an error.
# ---------------------------------------------------------------------------
print("== the examples of every screen: picked, filled in, run")
SCREENS = {
    "qc": (ui.qc_page_form, ui.qc_page_backend, "form_qc", None),
    "cluster": (ui.cluster_page_form, ui.cluster_page_backend, "form_cluster", None),
    "assembly": (ui.assembly_page_form, ui.assembly_page_backend, "form_assembly", None),
    "variant": (ui.variant_page_form, ui.variant_page_backend, "form_variant", None),
    "annotation": (ui.annotation_page_form, ui.annotation_page_backend, "form_annotation", None),
    "bloom": (ui.bloom_page_form, ui.bloom_page_backend, "form_bloom", None),
    "scaffold": (ui.scaffold_page_form, ui.scaffold_page_backend, "form_scaffold", None),
    "vcf": (ui.vcf_page_form, ui.vcf_page_backend, "form_vcf", None),
    "sort": (ui.sort_page_form, ui.sort_page_backend, "form_sort", "sort_task"),
    "validate": (ui.validate_page_form, ui.validate_page_backend, "form_validate", "validate_task"),
    "vcftools": (ui.vcftools_page_form, ui.vcftools_page_backend, "form_vcftools", "vcftools_task"),
    "fastq": (ui.fastq_page_form, ui.fastq_page_backend, "form_fastq", "fastq_task"),
    "kmercount": (ui.kmercount_page_form, ui.kmercount_page_backend, "form_kmercount", None),
    "convert": (ui.convert_page_form, ui.convert_page_backend, "form_convert", None),
    "findings": (ui.findings_page_form, ui.findings_page_backend, "form_findings", "findings_task"),
    "pgx": (ui.pgx_page_form, ui.pgx_page_backend, "form_pgx", "pgx_engine"),
    "lookup": (ui.lookup_page_form, ui.lookup_page_backend, "form_lookup", "lookup_task"),
    "assoc": (ui.assoc_page_form, ui.assoc_page_backend, "form_assoc", "assoc_task"),
    "contam": (ui.contam_page_form, ui.contam_page_backend, "form_contam", "contam_task"),
    "kin": (ui.kin_page_form, ui.kin_page_backend, "form_kin", "kin_task"),
    "chrom": (ui.chrom_page_form, ui.chrom_page_backend, "form_chrom", "chrom_task"),
    "anc": (ui.anc_page_form, ui.anc_page_backend, "form_anc", "anc_task"),
    "rnaseq": (ui.rnaseq_page_form, ui.rnaseq_page_backend, "form_rnaseq", "rnaseq_task"),
    "crispr": (ui.crispr_page_form, ui.crispr_page_backend, "form_crispr", "crispr_task"),
    "meta": (ui.meta_page_form, ui.meta_page_backend, "form_meta", "meta_task"),
    "amr": (ui.amr_page_form, ui.amr_page_backend, "form_amr", None),
    "outbreak": (ui.outbreak_page_form, ui.outbreak_page_backend, "form_outbreak", "outbreak_task"),
    "amplicon": (ui.amplicon_page_form, ui.amplicon_page_backend, "form_amplicon", "amplicon_task"),
    "tmb": (ui.tmb_page_form, ui.tmb_page_backend, "form_tmb", None),
    "ai": (ui.ai_page_form, ui.ai_page_backend, "form_ai", "ai_task"),
}
import shutil
RESEARCH_DEMO = os.environ.get("MMG_RESEARCH_DEMO")
if RESEARCH_DEMO and os.path.exists(RESEARCH_DEMO):
    shutil.copy(RESEARCH_DEMO, os.path.basename(RESEARCH_DEMO))
# The 16S examples (screen 27) need their own pack (MMG_16S_DEMO; published with the release since 2.3.0).
DEMO_16S = os.environ.get("MMG_16S_DEMO")
if DEMO_16S and os.path.exists(DEMO_16S):
    shutil.copy(DEMO_16S, os.path.basename(DEMO_16S))
# The human-genome examples need the demo data (MMG_HUMAN_DEMO = the path of
# mmgen_demo_human_v1.zip) and, for screen 19, PharmCAT (PHARMCAT_JAR = the
# pharmcat-<version>-all.jar, with its positions file next to it).
import shutil
HUMAN_DEMO = os.environ.get("MMG_HUMAN_DEMO")
if HUMAN_DEMO and os.path.exists(HUMAN_DEMO):
    shutil.copy(HUMAN_DEMO, os.path.basename(HUMAN_DEMO))
# MMGEN's own engine (screen 19) needs its tables (MMGEN_PGX_TABLES = their
# folder), and its BAM example the pgx demo pack (MMG_PGX_DEMO) and the
# IPD-IMGT/HLA sequences (MMGEN_HLA_DIR = the folder with A_gen.fasta...).
PGX_TABLES = os.environ.get("MMGEN_PGX_TABLES")
PGX_DEMO = os.environ.get("MMG_PGX_DEMO")
if PGX_DEMO and os.path.exists(PGX_DEMO):
    shutil.copy(PGX_DEMO, os.path.basename(PGX_DEMO))
HLA_DIR = os.environ.get("MMGEN_HLA_DIR")
_dl = cell_run(ui.pgx_page_download, get_page=lambda: "pgx", btn_get_hla=Choice(0), btn_get_pgx_tables=Choice(0),
               btn_get_pharmcat=Choice(0))
POOL.update(pgx_tables_dir=_dl[1]["pgx_tables_dir"], pgx_hla_dir=_dl[1]["pgx_hla_dir"])
if PGX_TABLES and os.path.exists(os.path.join(PGX_TABLES, "bam_genes.json")):
    # The screen downloads the IPD-IMGT/HLA release and the table zip by their SHA-256: the
    # constants in ui_main.py must be those of the tables (bam_genes.json pins the release).
    import hashlib, json, re
    _src = open(ui.__file__, encoding="utf-8").read()
    _ver = re.search(r'PGX_HLA_VERSION = "([^"]+)"', _src).group(1)
    _shas = dict(re.findall(r'"([ABC]_gen\.fasta)": "([0-9a-f]{64})"', _src))
    _hg = json.load(open(os.path.join(PGX_TABLES, "bam_genes.json"), encoding="utf-8"))["hla_genes"]
    _pins = {g["fasta_filename"]: (g.get("imgt_hla_version"), g.get("fasta_sha256")) for g in _hg}
    check("screen 19: the HLA release and SHA-256 sums in the screen are the ones the tables pin",
          _pins == {f: (_ver, s) for f, s in _shas.items()}, f"screen {_ver} {_shas} / tables {_pins}")
    _zip = os.path.join(os.path.dirname(os.path.abspath(PGX_TABLES)), "mmgen_pgx_tables_cpic-v1.61.0.zip")
    if os.path.exists(_zip):
        _want = re.search(r'PGX_TABLES_SHA256 = "([0-9a-f]{64})"', _src).group(1)
        _got = hashlib.sha256(open(_zip, "rb").read()).hexdigest()
        check("screen 19: PGX_TABLES_SHA256 is the SHA-256 of the table zip next to the tables", _got == _want, f"{_got} != {_want}")
JAR = os.environ.get("PHARMCAT_JAR")
if JAR and os.path.exists(JAR):
    os.makedirs("mmgen_tools/pharmcat", exist_ok=True)
    for f in os.listdir(os.path.dirname(JAR)):
        if f.endswith("-all.jar") or (f.startswith("pharmcat_positions_") and f.endswith(".vcf")):
            shutil.copy(os.path.join(os.path.dirname(JAR), f), os.path.join("mmgen_tools/pharmcat", f))
OPS = ("merge", "intersect", "subtract", "complement", "closest", "coverage", "slop")
# What the result of some examples must say (in the page text, without tags).
EXPECT = {
    ("qc", None, 0): ["Distinct k-mers extracted: 29,970"],
    ("qc", None, 1): ["K-mer spectrum of the reads", "Genome size, from the spectrum"],
    ("cluster", None, 0): ["Groups found: 3", "species_A (20)", "species_B (20)", "species_C (20)"],
    ("cluster", None, 1): ["kpneumoniae_KPNIH1_KPC-3", "saureus_USA300_MRSA"],
    ("assembly", None, 0): ["Longest contig: 30,000 bp", "The assembly graph", "Assembly graph (k = 31)", "isolated (a whole sequence)"],
    ("variant", None, 0): ["only in the patient): 93"],
    ("annotation", None, 0): ["MSTHEGENEISREAD"],
    ("annotation", None, 1): ["MSTHEGENEIS"],
    ("annotation", None, 2): ["MSTMKVKMKFLEKLI"],
    ("scaffold", None, 0): ["Longest continuous sequence: 20,000 bases", "Circular genomes detected: 1", "Contigs joined by their ends"],
    ("scaffold", None, 1): ["Longest continuous sequence: 24,000 bases", "Circular genomes detected: 0"],
    ("vcf", None, 0): ["Variants written: 600"],
    ("vcf", None, 1): ["Variants written: 600"],
    ("sort", "Sort & merge", 0): ["Records written: 450", ".tbi"],
    ("sort", "Sort & merge", 3): ["Records written: 600", ".bai"],
    ("sort", "Index & query a region", 0): ["record(s) in chr1:10000-20000"],
    ("sort", "merge", 0): ["Intervals read"],
    ("validate", "Validate a file", 0): ["4 error(s)", "line"],
    ("validate", "Validate a file", 1): ["needs at least 8"],
    ("validate", "Validate a file", 2): ["No problem found"],
    ("validate", "Rename chromosomes (chr1 ↔ 1)", 0): ["MT → chrM"],
    ("validate", "Compare reference genomes", 0): ["DIFFERENT reference genomes", "GRCh37"],
    ("validate", "Compare reference genomes", 1): ["Same reference"],
    ("vcftools", "Statistics", 0): ["500 in 3 sequence(s), with 3 sample(s)"],
    ("vcftools", "Statistics", 1): ["Nothing was read", "needs at least 8"],
    ("vcftools", "Filter", 2): ["Samples: 3 → 1"],
    ("vcftools", "Normalize", 1): ["Indels left-aligned: 3"],
    ("fastq", "Quality report", 0): ["Illumina Universal"],
    ("fastq", "Trim reads", 1): ["Pairs in / kept"],
    ("kmercount", None, 0): ["peak:", "What the spectrum says", "Genome size (from the spectrum)", "Spectrum (linear)"],
    ("convert", None, 1): ['parquet_options(numbers="text")'],
    ("findings", "ClinVar findings", 0): ["F11", "Pathogenic", "1 in one copy"],
    ("findings", "Gene panel from an annotation", 0): ["BRCA1", "Genes found: 13 of 13"],
    ("pgx", "PharmCAT", 0): ["CYP2C19", "Positions read", "PharmCAT data", "Genes called"],
    ("pgx", "MMGEN's engine (mmgen-pgx)", 0): ["CYP2C19", "mmgen-pgx", "Genes called", "G6PD", "absent from the VCF"],
    ("pgx", "MMGEN's engine (mmgen-pgx)", 1): ["CYP2D6", "*2/*4", "copy number", "*35:08/*38:01", "outside-call file",
                                                "CYP2D6 copies along the gene",
                                                "IPD-IMGT/HLA 3.66.0"],
    ("pgx", "MMGEN's engine (mmgen-pgx)", 2): ["CYP2C19", "*1/*2", "also possible", "*2/*35 (Poor Metabolizer)", "the evidence matrix",
                                                "rs4244285",
                                                "was chosen as the highest-scoring diplotype"],
    ("pgx", "A group of people and a protocol's medicines (mmgen-pgx)", 0): [
        "Pharmacogenomics of a group", "HG001", "HG007", "No CPIC guideline in the tables for", "gabapentin", "enoxaparin",
        "gene not called", "Each medicine for each person", "clopidogrel", "Each gene for each person", "Each person's ancestry"],
    ("lookup", "Look up variants in relatives (cascade testing)", 0): ["carries it", "does not carry it"],
    ("lookup", "Allele counts at a position (BAM)", 0): ["Depth", "Is it real? The usual checks", "Strand bias", "Where in the reads",
                                                         "chr1:24,940-25,060"],
    ("assoc", "Test variants against an outcome", 0): ["People analysed", "115", "2 genotyped but not in the table", "1 without age",
                                                       "rs12190259", "QQ plot of the p-values", "λGC", "Effect of the alleles with the smallest p", "What this study can detect"],
    ("assoc", "Test variants against an outcome", 1): ["rs12190259", "rs2030324"],
    ("assoc", "Test variants against an outcome", 2): ["yes is the event; 74 with it, 43 without", "firth 1", "odds ratio [95% CI]",
                                                       "with the rarer outcome"],
    ("assoc", "Test variants against an outcome", 3): ["p (family-wise)", "Freedman-Lane"],
    ("assoc", "Study power: how many people?", 0): ["People needed for 80% power", "164", "none in the range searched", "simulated"],
    ("assoc", "Study power: how many people?", 1): ["Power against the size of the study", "237"],
    ("assoc", "Study power: how many people?", 2): ["11,025", "5.00e-08"],
    ("rnaseq", "Build a transcript index", 0): ["Transcript index", "269", "140", "equivalence classes"],
    ("rnaseq", "Quantify from FASTQ (pseudoalignment)", 0): ["Quantified from the reads", "100.0%", "control_1", "treated_3",
                                                             "gene_counts.tsv", "The 15 most expressed genes", "(estimated)"],
    ("rnaseq", "Quantify from FASTQ (pseudoalignment)", 1): ["Quantified from the reads", "(given)", "single"],
    ("rnaseq", "See how one read is pseudoaligned", 0): ["unique", "DEMO0080-202", "Gene80", "each k-mer"],
    ("rnaseq", "Differential expression", 2): ["Gene7", "Gene19", "Gene33", "Gene12", "Gene48", "Gene71"],
    ("contam", "Screen samples for adventitious agents", 0): ["🔴 detected", "⚪ background", "🟡 below threshold", "✅ recovered",
                                                              "where the reads fall on each genome", "4/20 windows",
                                                              "Synthetic mycoplasma", "concentrated", "reagent_blank (negative control)"],
    ("contam", "Screen samples for adventitious agents", 1): ["Synthetic reagent bacterium", "screen_without_blank.html",
                                                              "No negative control was given"],
    ("contam", "Build a screening database", 0): ["Screening database", "Synthetic mycoplasma", "minimizers", "host"],
    ("kin", "Relatedness of several people", 0): ["parent and child", "HG002", "Kinship against IBS0"],
    ("kin", "Trio and paternity test", 0): ["is not excluded", "Combined paternity index"],
    ("kin", "Trio and paternity test", 1): ["is excluded"],
    ("kin", "Trio and paternity test", 2): ["is not excluded"],
    ("chrom", "Genetic sex", 0): ["XX"],
    ("chrom", "Genetic sex", 1): ["XY"],
    ("chrom", "Genetic sex", 2): ["Undetermined"],
    ("chrom", "Runs of homozygosity", 0): ["Runs: 2", "as plink --homozyg", "Runs of homozygosity along the chromosomes"],
    ("chrom", "Runs of homozygosity", 2): ["Runs: 2", "strict: any heterozygous call ends a run"],
    ("anc", "Ancestry", 0): ["HG005", "EAS", "best single group", "STRUCTURE-style"],
    ("anc", "Polygenic score", 0): ["percentile in EUR", "PGS000004"],
    ("anc", "Build a reference panel", 0): ["Panel written", "AFR, AMR, EAS, EUR, SAS"],
    ("rnaseq", "Count reads per gene", 0): ["Reads counted", "ctrl_1", "treat_2"],
    ("rnaseq", "Differential expression", 0): ["gene03", "gene10", "Volcano plot", "MA plot", "Samples: PCA of the 500 most variable genes",
                                                "Distances between samples", "p-values of the genes tested"],
    ("rnaseq", "Differential expression", 1): ["3,013 up, 2,539 down", "CACNB2", "Volcano plot"],
    ("rnaseq", "Test gene sets", 0): ["Gene sets", "GLUCOCORTICOID_RESPONSE", "⬆️ up", "HOUSEKEEPING", "competitive"],
    ("rnaseq", "Test gene sets", 1): ["built-in: spinal_cord_injury", "laminin_receptors", "glial_scar_cspgs", "regeneration_associated"],
    ("crispr", "Count guides in reads", 0): ["Guides counted", "58.1%"],
    ("crispr", "Rank genes (MAGeCK RRA)", 0): ["ACTR8", "Negative selection"],
    ("crispr", "Rank genes (MAGeCK RRA)", 1): ["RPL3", "7,115 genes"],
    ("meta", "Which species, and how much", 0): ["Pseudomonas aeruginosa", "Staphylococcus aureus", "Saccharomyces cerevisiae"],
    ("meta", "Which species, and how much", 1): ["Listeria monocytogenes"],
    ("meta", "Classify the reads", 0): ["Reads classified", "Escherichia coli"],
    ("meta", "Sketch reference genomes", 0): ["Sketches written", "USA300"],
    ("amr", None, 0): ["blaKPC-3", "CARBAPENEM"],
    ("amr", None, 1): ["blaNDM-1"],
    ("amr", None, 2): ["mecA"],
    ("outbreak", "Distances and a tree", 0): ["Tree built", "rooted on", "B.1.617.2_Delta", "Newick"],
    ("outbreak", "Distances and a tree", 1): ["KPNIH1 (ST258)", "30684/NJST258_2 (ST258)", "Closest pair"],
    ("outbreak", "Distances and a tree", 2): ["saureus_USA300_MRSA", "ANI (%) of every pair"],
    ("outbreak", "Mutations of viral genomes", 0): ["Mutations found", "20 of 20 aligned", "S:D614G", "ORF1b:P314L", "S:N501Y",
                                                     "Mutations along the genome"],
    ("outbreak", "In-silico PCR", 0): ["PCR done", "CDC_N1", "CDC_N2", "40 of 100", "A_MN985325.1"],
    ("outbreak", "In-silico PCR", 1): ["mecA", "16S_V4", "533"],
    ("amplicon", "Amplicon variants (ASVs)", 0): ["ASVs found", "Mock", "F3D150", "ASV_1"],
    ("amplicon", "Taxonomy", 0): ["ASVs classified", "Muribaculaceae", "Composition of each sample", "Taxonomy (share of the reads)"],
    ("amplicon", "Diversity", 0): ["Diversity measured", "PERMANOVA", "Early", "Late", "Principal coordinates", "Rarefaction curves",
                                   "95% ellipse"],
    ("amplicon", "Diversity", 1): ["PERMANOVA", "Jaccard"],
    ("amplicon", "Train a classifier", 0): ["Classifier trained", "my_v4_classifier.model"],
    ("tmb", None, 0): ["TMB measured", "mutations per megabase", "HCC1395"],
    ("tmb", None, 1): ["coding mutations per megabase", "not in a coding sequence"],
    ("ai", "Describe your problem", 0): ["The AI's answer", "name chromosomes differently", "tokens in"],
    ("ai", "Let the AI run a task", 0): ["The task went to the AI"],
}
# The warnings a form stops with; an example must not hit any of them.
STOPS = ("Select ", "cannot be", "must be", "Type the", "Type an", "is needed")
# Words of a stop that are, in these examples, the engine's own warning (what the example shows).
STOPS_ALLOWED = {("contam", "Screen samples for adventitious agents", 1): {"cannot be"}}
before = set(os.listdir("."))
EXAMPLES = examples[1]["EXAMPLES"]
check("every screen has examples", {p for p, _ in EXAMPLES} == set(SCREENS), sorted({p for p, _ in EXAMPLES} ^ set(SCREENS)))
for (page, task), items in EXAMPLES.items():
    form_cell, backend_cell, form_name, task_var = SCREENS[page]
    if page == "pgx" and task.startswith("A group"):
        backend_cell = ui.pgx_group_backend  # the group has a cell of its own
    if page == "rnaseq" and task in ("Build a transcript index", "Quantify from FASTQ (pseudoalignment)", "See how one read is pseudoaligned"):
        backend_cell = ui.rnaseq_fastq_backend  # the tasks that start from the reads
    picker = next(k for k in form_cell.refs if k.endswith("_example"))
    for i, (label, _about, _values) in enumerate(items):
        flat = [x for v in _values.values() for x in (v if isinstance(v, list) else [v])]
        human = any(isinstance(x, str) and x.startswith("mmgen_examples/human") for x in flat)
        if human and not (HUMAN_DEMO and os.path.exists(HUMAN_DEMO)):
            print(f"  skip  {page}: {label} (set MMG_HUMAN_DEMO to the demo data zip)")
            continue
        if any(isinstance(x, str) and x.startswith("mmgen_examples/research_16s/") for x in flat) \
                and not (DEMO_16S and os.path.exists(DEMO_16S)):
            print(f"  skip  {page}: {label} (set MMG_16S_DEMO to the 16S demo data zip)")
            continue
        research = any(isinstance(x, str) and x.startswith("mmgen_examples/research/") for x in flat)
        if research and not (RESEARCH_DEMO and os.path.exists(RESEARCH_DEMO)):
            print(f"  skip  {page}: {label} (set MMG_RESEARCH_DEMO to the research demo zip)")
            continue
        if page == "pgx" and task == "PharmCAT" and not (JAR and os.path.exists(JAR)):
            print(f"  skip  {page}: {label} (set PHARMCAT_JAR)")
            continue
        if page == "pgx" and task != "PharmCAT" and not (PGX_TABLES and os.path.isdir(PGX_TABLES)):
            print(f"  skip  {page}: {label} (set MMGEN_PGX_TABLES to the mmgen-pgx tables)")
            continue
        if any(isinstance(x, str) and x.startswith("mmgen_examples/pgx/") for x in flat)                 and not (PGX_DEMO and os.path.exists(PGX_DEMO) and HLA_DIR and os.path.isdir(HLA_DIR)):
            print(f"  skip  {page}: {label} (set MMG_PGX_DEMO to the pgx demo zip and MMGEN_HLA_DIR to the HLA sequences)")
            continue
        given = {"get_page": lambda p=page: p, picker: Choice(str(i)), "interval_op": Choice(task if task in OPS else "merge")}
        if task_var:
            given[task_var] = Choice("Interval operations (BED)" if task in OPS else task)
        name = f"{page}{' / ' + task if task else ''}: {label}"
        r = cell_run(form_cell, **given)
        if isinstance(r, BaseException) or form_name not in r[1]:
            check(name + " (the form)", False, repr(r)[:300])
            continue
        about_shown = label in shown(r)
        given[form_name] = Form(**r[1][form_name].element.value)
        rb = cell_run(backend_cell, **{k: v for k, v in given.items() if k in backend_cell.refs})
        page_text = shown(rb)
        missing = [e for e in EXPECT.get((page, task, i), []) if e not in page_text]
        # The engine's own HTML report is embedded at the end: its texts are not the screen's warnings.
        hit = [s for s in STOPS if s in page_text.split("The engine's full report")[0]
               and s not in STOPS_ALLOWED.get((page, task, i), ())]
        check(name, about_shown and not isinstance(rb, BaseException) and "ERROR" not in page_text and not hit and not missing,
              f"about shown: {about_shown}; stops: {hit}; missing: {missing}; {page_text[:400]}")
check("the examples write only inside mmgen_examples/", set(os.listdir(".")) - before <= {"mmgen_examples"},
      sorted(set(os.listdir(".")) - before))
check("the example files are written once, not on every pick", examples[1]["make_examples"]() is False)

# ---------------------------------------------------------------------------
# The AI: screen 29 and the panel at the bottom of the other screens, with
# the stand-in for Claude above (nothing is sent anywhere).
# ---------------------------------------------------------------------------
print("== screen 29 and the AI panel")
ai_page = lambda: "ai"  # noqa: E731
r = cell_run(ui.ai_page_ui, get_page=ai_page, set_page=lambda _: None)
check("  screen 29 opens with its three tasks", not isinstance(r, BaseException) and "Describe your problem" in shown(r)
      and "Let the AI run a task" in shown(r) and "The AI can be wrong" in shown(r), shown(r)[:300])
AI_SETTINGS.update(consent=False)
r = cell_run(ui.ai_page_status, get_page=ai_page)
check("  the status line says the AI is off before the agreement", "AI: off" in shown(r), shown(r)[:200])
# Settings: saved into the dict, the redraw counter moves, nothing is sent.
n_req, rev = len(AI_REQUESTS), AI_REV.value
settings = Choice("Settings and the session record")
r = cell_run(ui.ai_page_form, get_page=ai_page, ai_task=settings)
check("  the settings form has the key, the model and the two agreements",
      not isinstance(r, BaseException) and "form_ai" in r[1] and "Anthropic API key" in shown(r) and "human-genome results" in shown(r), shown(r)[:300])
r = cell_run(ui.ai_page_backend, get_page=ai_page, ai_task=settings,
             form_ai=Form(key=" sk-test ", model="claude-sonnet-5-5", effort="high", level="researcher", consent=True, human=True))
check("  saving the settings: kept in memory, redrawn, nothing sent",
      "Saved for this session" in shown(r) and AI_SETTINGS["key"] == "sk-test" and AI_SETTINGS["model"] == "claude-sonnet-5-5"
      and AI_SETTINGS["human"] and AI_REV.value == rev + 1 and len(AI_REQUESTS) == n_req, (AI_SETTINGS, shown(r)[:200]))
r = cell_run(ui.ai_page_backend, get_page=ai_page, ai_task=settings,
             form_ai=Form(key="", model="claude-opus-5-5", effort="medium", level="student", consent=False, human=True))
check("  the human-genome agreement needs the first one", not AI_SETTINGS["consent"] and not AI_SETTINGS["human"], AI_SETTINGS)
AI_SETTINGS.update(consent=True, human=False, key="")

# The panel: a session record as the engine, a validation run on screen 10.
rec = ai_lib.Recorder(mmg)
val = Choice(None)
r = cell_run(ui.ai_panel_ui, get_page=lambda: "validate", mmg=rec, set_page=lambda _: None)
check("  the panel is at the bottom of a screen, with the form", not isinstance(r, BaseException) and "form_ai_ask" in r[1]
      and "Ask the AI about this result" in raw(r) and "What is sent to Anthropic" in raw(r), raw(r)[:300])
form_ask = r[1]["form_ai_ask"]
for page in ("home", "guide", "ai"):
    rr = cell_run(ui.ai_panel_ui, get_page=lambda p=page: p, mmg=rec, set_page=lambda _: None)
    check(f"  no panel on '{page}'", not isinstance(rr, BaseException) and rr[0] is None and "form_ai_ask" not in rr[1], repr(rr)[:200])
rr = cell_run(ui.ai_panel_ui, get_page=lambda: "validate", mmg=mmg, set_page=lambda _: None)
check("  no panel with an engine that has no session record", not isinstance(rr, BaseException) and rr[0] is None, repr(rr)[:200])
ask = Form(question="", level="student")
r = cell_run(ui.ai_panel_answer, get_page=lambda: "validate", mmg=rec, form_ai_ask=ask)
check("  asking before any result: says to run a task first", stopped(r, "Run a task on this screen first"), shown(r)[:200])
rec.context = "validate"
vrep = rec.validate_file(os.path.abspath("mmgen_examples/with_problems.vcf"))
n_req = len(AI_REQUESTS)
r = cell_run(ui.ai_panel_answer, get_page=lambda: "validate", mmg=rec, form_ai_ask=ask)
req = AI_REQUESTS[-1] if len(AI_REQUESTS) > n_req else {}
user = req.get("messages", [{}])[0].get("content", "")
check("  the answer is shown, with the model and the tokens", "numbers look as expected" in shown(r) and "validate_file" in shown(r)
      and "tokens in" in shown(r), shown(r)[:300])
check("  what went out: the file name, never the folder", "with_problems.vcf" in user and os.path.abspath(".") not in user
      and "mmgen_examples" not in user and str(vrep["errors"]) in user, user[:400])
check("  the request: model, effort and fallbacks from the settings", req.get("model") == "claude-opus-5-5"
      and req.get("output_config") == {"effort": "medium"} and req.get("fallbacks") == "default", {k: req.get(k) for k in ("model", "fallbacks")})
AI_SETTINGS.update(consent=False)
r = cell_run(ui.ai_panel_answer, get_page=lambda: "validate", mmg=rec, form_ai_ask=ask)
check("  without the agreement nothing is sent", stopped(r, "The AI is off") and len(AI_REQUESTS) == n_req + 1, shown(r)[:200])
r = cell_run(ui.ai_panel_ui, get_page=lambda: "validate", mmg=rec, set_page=lambda _: None)
check("  and the panel offers the settings instead of the form", "The AI is off" in raw(r) and "Open the AI settings" in raw(r), raw(r)[:300])
AI_SETTINGS.update(consent=True)
# A human-genome result needs the second agreement, and the answer carries the notice.
rec.calls.append(ai_lib.Call("relatedness", (["family.vcf"],), {}, result={"pairs": [{"a": "child", "b": "father", "kinship": 0.25}]},
                             context="kin"))
r = cell_run(ui.ai_panel_answer, get_page=lambda: "kin", mmg=rec, form_ai_ask=ask)
check("  a human-genome result is not sent without the second agreement", stopped(r, "human-genome result"), shown(r)[:200])
AI_SETTINGS.update(human=True)
r = cell_run(ui.ai_panel_answer, get_page=lambda: "kin", mmg=rec, form_ai_ask=Form(question="Is he the father?", level="researcher"))
check("  with it: the answer, and the notice that the authors take no responsibility",
      "take no responsibility" in shown(r) and "Is he the father?" in shown(r) and "human-genome result" in AI_REQUESTS[-1]["messages"][0]["content"],
      shown(r)[:300])
AI_SETTINGS.update(human=False)

# Describe your problem: the plan, the script checked against the catalog, the downloads.
plan = Choice("Describe your problem")
r = cell_run(ui.ai_page_backend, get_page=ai_page, ai_task=plan, mmg=rec,
             form_ai=Form(problem="My two VCFs share nothing.", folder="mmgen_examples", level="student"))
check("  the plan is shown, with the script and the answer to download", "name chromosomes differently" in shown(r)
      and "mmgen_plan.py" in raw(r) and "mmgen_plan.md" in raw(r) and "Check the script" not in shown(r), shown(r)[:300])
check("  the AI saw the names of the example files", "cohort.vcf" in AI_REQUESTS[-1]["messages"][0]["content"])
r = cell_run(ui.ai_page_backend, get_page=ai_page, ai_task=plan, mmg=rec, form_ai=Form(problem="  ", folder=".", level="student"))
check("  an empty problem: a warning, nothing sent", stopped(r, "Type your problem first"), shown(r)[:200])

# The agent: a round proposed, nothing run until approved, then the report.
agent_task = Choice("Let the AI run a task")
AGENT.value = None
r = cell_run(ui.ai_page_backend, get_page=ai_page, ai_task=agent_task, mmg=rec,
             form_ai=Form(task="Validate with_problems.vcf.", folder="mmgen_examples", human=False, rounds=4))
held = AGENT.value
check("  the agent starts and proposes a round", "The task went to the AI" in shown(r) and held is not None and held[0].pending
      and held[0].pending[0].name == "validate_file", shown(r)[:300])
n_calls = len(rec.calls)
rec.context = "ai"  # the interface sets it when screen 29 opens
r = cell_run(ui.ai_page_agent_view, get_page=ai_page, ai_task=agent_task, get_agent=lambda: held)
check("  the proposal is shown with the buttons; nothing has run", "Proposed calls" in shown(r) and "mmg.validate_file" in shown(r)
      and "btn_agent_run" in r[1] and len(rec.calls) == n_calls, shown(r)[:300])
r = cell_run(ui.ai_page_agent_act, btn_agent_run=Choice(True), btn_agent_decline=Choice(False), get_agent=lambda: held)
agent = AGENT.value[0]
check("  approved: the call runs (and is recorded) and the AI writes its report",
      agent.done and "validator found" in agent.final_text and rec.calls[-1].name == "validate_file" and rec.calls[-1].page == "ai"
      and AGENT.value[1] == 1, (agent.done, agent.final_text, rec.calls[-1].name))
r = cell_run(ui.ai_page_agent_view, get_page=ai_page, ai_task=agent_task, get_agent=lambda: AGENT.value)
check("  the view shows the result and the report", "The agent's report" in shown(r) and "mmgen_agent_transcript.md" in raw(r), shown(r)[:300])

# The session record: every call, as a table and as a script.
r = cell_run(ui.ai_page_session, get_page=ai_page, ai_task=settings, mmg=rec)
check("  the session record lists the calls and offers the script", "The session record" in shown(r) and "validate_file" in raw(r)
      and "mmgen_session.py" in raw(r), shown(r)[:300])
check("  the script repeats the calls", f"r = mmg.validate_file({os.path.abspath('mmgen_examples/with_problems.vcf')!r})" in rec.script(),
      rec.script()[:500])

print()
print(f"{done - len(failures)} of {done} checks passed")
if failures:
    print("FAILED:", failures)
    sys.exit(1)
