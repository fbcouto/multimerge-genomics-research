# Test of the marimo interface IN A REAL BROWSER (the system's Edge or Chrome,
# driven by Playwright): starts `marimo run`, enters each screen through the
# menu, fills in and submits the forms, and checks what shows on the page and
# what was written to disk. It covers what test_ui.py cannot see: the
# navigation, the rendered components, the websocket and screen 2 (clustering).
#
# Usage: python trials/python/test_browser.py ui_main.py [folder_for_screenshots]
#
# Needs marimo, polars, numpy, scikit-learn, playwright and the extension
# installed. `playwright install` is not needed: the test uses the browser the
# machine already has (MMG_BROWSER=msedge|chrome picks one).
# It works in a temporary folder, with data it creates itself.
import contextlib
import gzip
import os
import random
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import urllib.request

try:
    from playwright.sync_api import sync_playwright
except ImportError:
    sys.exit("Playwright is missing: pip install playwright")

import polars as pl

import multimerge_genomics as mmg

for _stream in (sys.stdout, sys.stderr):
    _stream.reconfigure(errors="replace")  # Windows consoles that are not UTF-8

failures = []
done = 0


def check(name, cond, detail=""):
    global done
    done += 1
    if not cond:
        failures.append(name)
    print(("  ok    " if cond else "  FAIL  ") + name + ("" if cond else f"  {str(detail)[:300]}"))


if len(sys.argv) < 2:
    sys.exit("Usage: python trials/python/test_browser.py ui_main.py [folder_for_screenshots]")
UI_PATH = os.path.abspath(sys.argv[1])
SHOTS = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else None
if SHOTS:
    os.makedirs(SHOTS, exist_ok=True)

# ---------------------------------------------------------------------------
# The data, in an empty folder: the screens list the files of the current directory.
# ---------------------------------------------------------------------------
T = tempfile.mkdtemp(prefix="mmg_browser_")
shutil.copy(UI_PATH, os.path.join(T, "ui_main.py"))
os.chdir(T)
random.seed(1)
K = 55
GENOME_LEN = 40000  # the "Generate simulated patient" button only plants its mutations above 30,000 bp
genome = "".join(random.choice("ACGT") for _ in range(GENOME_LEN))
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}

# Screen 4 makes the patient from a file with exactly this name.
with open("ecoli.fasta", "w", newline="\n") as f:
    f.write(">synthetic_reference\n")
    for i in range(0, GENOME_LEN, 70):
        f.write(genome[i:i + 70] + "\n")

with open("reads.fastq", "w", newline="\n") as f:
    for i in range(60):
        start = random.randrange(0, GENOME_LEN - 120)
        f.write(f"@read_{i}\n{genome[start:start + 120]}\n+\n{'I' * 120}\n")

vcf_header = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]


def make_records(n, tag):
    return [f"{random.choice(['chr1', 'chr2', 'chr10'])}\t{random.randrange(1, 100000)}\t{tag}{i}\tA\tT\t{random.randrange(1, 99)}\tPASS\t." for i in range(n)]


def read_vcf(name):
    with open(name, newline="\n") as f:
        return [l for l in f.read().split("\n") if l and not l.startswith("#")]


BATCHES = {"sample_a.vcf": make_records(3000, "a"), "sample_b.vcf": make_records(2000, "b"), "sample_c.vcf": make_records(1000, "c")}
for _name, _recs in BATCHES.items():
    with open(_name, "w", newline="\n") as f:
        f.write("\n".join(vcf_header + _recs) + "\n")

# Screen 2: two text files of 31-mers, and two .arrow files in the old format
# (one 64-bit k-mer per row), which is what that screen can read.
for _name in ("species_A.txt", "species_B.txt"):
    _seq = "".join(random.choice("ACGT") for _ in range(2030))
    with open(_name, "w", newline="\n") as f:
        f.write("\n".join(_seq[i:i + 31] for i in range(2000)) + "\n")
pl.DataFrame({"K-mers": [random.getrandbits(62) for _ in range(1500)]}).write_ipc("legacy_A.arrow")
pl.DataFrame({"K-mers": [random.getrandbits(62) for _ in range(1700)]}).write_ipc("legacy_B.arrow")

# Screens 9-14: intervals and a genome file, a VCF with defects, two .fai of
# different assemblies, an indel inside a repeat with its genome, a pair of
# FASTQ files with an adapter, and a VCF whose QUAL is not in its shortest form.
CHROMS = {"chr1": 300_000, "chr2": 200_000, "chr10": 100_000}
with open("genome.fai", "w", newline="\n") as f:
    f.writelines(f"{c}\t{n}\t0\t60\t61\n" for c, n in CHROMS.items())
for _name, _n in (("peaks.bed", 300), ("genes.bed", 200)):
    with open(_name, "w", newline="\n") as f:
        for _i in range(_n):
            _c = random.choice(list(CHROMS))
            _s = random.randrange(0, CHROMS[_c] - 2000)
            f.write(f"{_c}\t{_s}\t{_s + random.randrange(50, 1500)}\tname{_i}\n")
_contigs = "".join(f"##contig=<ID={c},length={n}>\n" for c, n in CHROMS.items())
with open("broken.vcf", "w", newline="\n") as f:
    f.write("##fileformat=VCFv4.2\n" + _contigs + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    f.write("chr1\t500\ta\tA\tG\t50\tPASS\t.\nchr1\t100\tb\tA\tG\t50\tPASS\t.\n1\t200\tc\tA\tG\t50\tPASS\t.\n")
with open("hg38.fai", "w", newline="\n") as f:
    f.write("chr1\t248956422\t0\t0\t0\nchr2\t242193529\t0\t0\t0\n")
with open("hg19.fai", "w", newline="\n") as f:
    f.write("1\t249250621\t0\t0\t0\n2\t243199373\t0\t0\t0\n")
_left = "".join(random.choice("ACG") for _ in range(100))
SMALL_GENOME = _left + "C" + "TTTTT" + "G" + "".join(random.choice("ACGT") for _ in range(400))
with open("small_genome.fa", "w", newline="\n") as f:
    f.write(">chr1\n" + "\n".join(SMALL_GENOME[i:i + 60] for i in range(0, len(SMALL_GENOME), 60)) + "\n")
with open("indel.vcf", "w", newline="\n") as f:
    f.write(f"##fileformat=VCFv4.2\n##contig=<ID=chr1,length={len(SMALL_GENOME)}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    f.write("chr1\t105\tdel\tTT\tT\t50\tPASS\t.\n")
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
for _mate in (1, 2):
    with gzip.open(f"pair_R{_mate}.fastq.gz", "wt", newline="\n") as f:
        for _i in range(2000):
            _insert = "".join(random.choice("ACGT") for _ in range(random.randrange(40, 100)))
            _seq = (_insert + ADAPTER)[:100]
            f.write(f"@p{_i}\n{_seq}\n+\n{''.join(chr(33 + random.randrange(5, 41)) for _ in _seq)}\n")
with open("five_columns.vcf", "w", newline="\n") as f:  # not a VCF: no REF, ALT, FILTER, INFO
    f.write("##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tQUAL\tGT\n1\t100\trs1\t10.5\t0/0\n1\t200\trs2\t50.2\t0/1\n")
with open("quals.vcf", "w", newline="\n") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    f.write("chr1\t10\ta\tA\tG\t12.50\tPASS\t.\nchr1\t20\tb\tC\tT\t30.0\tPASS\t.\n")


def records(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", newline="\n") as f:
        return [l for l in f.read().split("\n") if l and not l.startswith("#")]

# ---------------------------------------------------------------------------
# The server
# ---------------------------------------------------------------------------
with socket.socket() as _s:
    _s.bind(("127.0.0.1", 0))
    PORT = _s.getsockname()[1]
URL = f"http://127.0.0.1:{PORT}/"
log = open(os.path.join(T, "server.log"), "w", encoding="utf-8", errors="replace")
server = subprocess.Popen(
    [sys.executable, "-m", "marimo", "run", "ui_main.py", "--headless", "--port", str(PORT), "--no-token"],
    cwd=T, stdout=log, stderr=subprocess.STDOUT,
)


def shut_down():
    server.terminate()
    try:
        server.wait(timeout=15)
    except subprocess.TimeoutExpired:
        server.kill()
    log.close()


def log_tail():
    log.flush()
    with open(os.path.join(T, "server.log"), encoding="utf-8", errors="replace") as f:
        return f.read()[-1500:]


_up = False
for _ in range(90):
    if server.poll() is not None:
        break
    try:
        _up = urllib.request.urlopen(URL, timeout=3).status == 200
        break
    except Exception:
        time.sleep(1)
if not _up:
    print(log_tail())
    shut_down()
    sys.exit("`marimo run` did not answer.")

# ---------------------------------------------------------------------------
# Browser helpers
# ---------------------------------------------------------------------------
MENU = {
    1: ("1. QC & Ingestion", "Data Ingestion & QC"),
    2: ("2. Taxonomic Clustering", "Taxonomic Clustering"),
    3: ("3. Genome Assembly", "Genome Assembly"),
    4: ("4. Variant Calling", "Variant Calling"),
    5: ("5. Functional Annotation", "Functional Annotation"),
    6: ("6. Bloom Filters", "Probabilistic Bloom Filters"),
    7: ("7. Genome Scaffolding", "Genome Scaffolding"),
    8: ("8. VCF Multimerge", "VCF Multimerge"),
    9: ("9. Sort, Index & Intervals", "9. Sorting, Indexing & Intervals"),
    10: ("10. Validation & Chromosome Names", "10. Validation & Chromosome Names"),
    11: ("11. VCF Stats, Filter & Normalize", "11. VCF Statistics, Filtering & Normalization"),
    12: ("12. FASTQ Quality & Trimming", "12. FASTQ Quality Report & Trimming"),
    13: ("13. K-mer Counting", "13. K-mer Counting"),
    14: ("14. Formats & Conversion", "14. Formats & Conversion"),
}

# Words that would mean a piece of the interface is still in Portuguese.
PORTUGUESE = re.compile(r"\b(ficheiros?|arquivos?|voltar|aguardando|selecione|m[oó]dulo|mem[oó]ria|iniciar|n[aã]o|s[aã]o|voc[eê]|sucesso|conclu[ií]d[oa])\b", re.I)

# The visible text of the page, going through the shadow roots of the components.
JS_TEXT = """
() => {
  const parts = [];
  const walk = (node) => {
    if (node.nodeType === Node.TEXT_NODE) {
      const t = node.textContent.trim();
      if (t) parts.push(t);
      return;
    }
    if (node.nodeType === Node.ELEMENT_NODE) {
      const tag = node.tagName.toLowerCase();
      if (tag === 'script' || tag === 'style' || tag === 'option') return;
      const st = getComputedStyle(node);
      if (st.display === 'none' || st.visibility === 'hidden') return;
      if (node.shadowRoot) walk(node.shadowRoot);
    } else if (node.nodeType !== Node.DOCUMENT_FRAGMENT_NODE) {
      return;
    }
    for (const child of node.childNodes) walk(child);
  };
  walk(document.body);
  return parts.join(' | ');
}
"""

js_errors = []
n_shot = 0


@contextlib.contextmanager
def page(browser):
    """A new tab: in `marimo run`, every tab is a session, which starts at the menu."""
    pg = browser.new_page(viewport={"width": 1280, "height": 1500})
    pg.on("pageerror", lambda e: js_errors.append(str(e)))
    pg.on("console", lambda m: js_errors.append(m.text) if m.type == "error" else None)
    pg.goto(URL, wait_until="domcontentloaded")
    menu_button(pg, 1).wait_for(timeout=90000)
    try:
        yield pg
    finally:
        pg.close()


def menu_button(pg, n):
    return pg.get_by_role("button", name=MENU[n][0])


def back(pg):
    return pg.get_by_role("button", name="Back to Menu")


def enter(pg, n):
    menu_button(pg, n).click()
    back(pg).wait_for(timeout=30000)
    pg.wait_for_timeout(600)


def shot(pg, name):
    global n_shot
    if SHOTS:
        n_shot += 1
        pg.wait_for_timeout(300)
        pg.screenshot(path=os.path.join(SHOTS, f"{n_shot:02d}_{name}.png"), full_page=True)


def text_of(pg):
    return pg.evaluate(JS_TEXT)


def wait_for(pg, pattern, seconds=120):
    """Waits for a text (regex) to show up; returns the page text, or None."""
    try:
        pg.get_by_text(re.compile(pattern)).first.wait_for(timeout=seconds * 1000)
        pg.wait_for_timeout(500)
        return text_of(pg)
    except Exception:
        return None


def number(t, label):
    m = re.search(re.escape(label) + r"\D*?(\d[\d,]*)", t or "")
    return int(m.group(1).replace(",", "")) if m else None


def multiselect(pg, options):
    pg.locator("button:has-text('Select...')").first.click()
    pg.wait_for_timeout(400)
    for op in options:
        pg.get_by_role("option", name=op, exact=True).click()
        pg.wait_for_timeout(150)
    pg.keyboard.press("Escape")
    pg.wait_for_timeout(300)


def choose(pg, index, label):
    pg.locator("select:visible").nth(index).select_option(label=label)
    pg.wait_for_timeout(200)


def submit(pg, label):
    pg.get_by_role("button", name=label).click()


def radio(pg, name):
    """Picks a task on screens 9-14; the form below is rebuilt for it."""
    pg.get_by_role("radio", name=name, exact=True).click()
    pg.wait_for_timeout(1200)


def placeholder(pg, text):
    return pg.get_by_placeholder(text, exact=True)


def code_text(pg):
    """The last code block of the page: the Python call under a result."""
    blocks = pg.locator("pre:visible")
    return blocks.last.inner_text() if blocks.count() else ""


# The icon-only buttons that are visible: for the diagnosis, when the table's
# "next page" button is not found.
JS_BUTTONS = """
() => [...document.querySelectorAll('button')].filter(b => b.getBoundingClientRect().width > 0 && !b.innerText.trim())
        .map(b => [b.getAttribute('aria-label'), b.getAttribute('title'), b.querySelector('svg')?.getAttribute('class')].join('/')).join(' ; ')
"""


def next_page(pg):
    """The table's "next page" button, whose name changes with the marimo version."""
    for loc in (pg.get_by_role("button", name=re.compile("next page", re.I)),
                pg.locator("button:has(svg.lucide-chevron-right):visible")):
        if loc.count():
            return loc.first
    return None


def ingest(browser, files):
    """Screen 1 with the given files; returns the page text with the result."""
    with page(browser) as pg:
        enter(pg, 1)
        multiselect(pg, files)
        submit(pg, "Start batch ingestion")
        t = wait_for(pg, "Batch processed successfully|ERROR", 300)
        shot(pg, "screen1_" + files[0].split(".")[0])
        return t


def launch(p):
    preferred = os.environ.get("MMG_BROWSER")
    for channel in ([preferred] if preferred else ["msedge", "chrome"]):
        try:
            return p.chromium.launch(channel=channel, headless=True), channel
        except Exception:
            continue
    return p.chromium.launch(headless=True), "Playwright's chromium"


# ---------------------------------------------------------------------------
# The script
# ---------------------------------------------------------------------------
try:
    # The version comes from the metadata: importing marimo here swaps asyncio's
    # event loop on Windows for one that cannot launch subprocesses, and
    # Playwright needs to launch its own.
    from importlib.metadata import version
    engine = mmg.engine_info()
    WHERE = "GPU" if engine["gpu"] else "CPU"
    with sync_playwright() as p:
        browser, channel = launch(p)
        print(f"marimo {version('marimo')}, browser: {channel} {browser.version}, engine on the {WHERE}")

        print("== navigation")
        with page(browser) as pg:
            check("tab title: MMGEN Suite", pg.title() == "MMGEN Suite", pg.title())
            check(f"the menu shows the {len(MENU)} modules", all(menu_button(pg, n).is_visible() for n in MENU))
            t = text_of(pg)
            check("the menu says where the engine runs", "running on:" in t and WHERE in t, t[-200:])
            check("menu: no Portuguese left on the page", not PORTUGUESE.search(t), PORTUGUESE.findall(t))
            shot(pg, "menu")
            for n, (_, title) in MENU.items():
                enter(pg, n)
                t = text_of(pg)
                check(f"screen {n}: opens, with its title drawn as a heading", pg.locator("h1", has_text=title).count() >= 1, t[:200])
                check(f"screen {n}: the menu leaves the page", not any(menu_button(pg, m).is_visible() for m in MENU))
                check(f"screen {n}: no Portuguese left on the page", not PORTUGUESE.search(t), PORTUGUESE.findall(t))
                shot(pg, f"screen{n}_empty")
                back(pg).click()
                try:
                    menu_button(pg, 1).wait_for(timeout=30000)
                except Exception:
                    pass
                pg.wait_for_timeout(400)
                check(f"screen {n}: 'Back to Menu' brings the menu back and removes the screen",
                      menu_button(pg, 1).is_visible() and not back(pg).is_visible() and pg.locator("h1", has_text=title).count() == 0)

        print("== screen 1: ingestion (FASTA and FASTQ -> k-mers -> .arrow)")
        t = ingest(browser, ["ecoli.fasta"])
        n = number(t, "extracted:")
        check(f"ecoli.fasta -> ecoli_K{K}.arrow with {GENOME_LEN - K + 1} k-mers", n == GENOME_LEN - K + 1 and os.path.exists(f"ecoli_K{K}.arrow"), n)
        t = ingest(browser, ["reads.fastq"])
        n = number(t, "extracted:")
        check("reads.fastq -> k-mers", bool(n) and os.path.exists(f"reads_K{K}.arrow"), n)
        with page(browser) as pg:
            enter(pg, 1)
            submit(pg, "Start batch ingestion")
            check("submitting with no file: a warning, no error", wait_for(pg, "Select at least one file", 30) is not None, text_of(pg)[-300:])

        print("== screen 4: generate the simulated patient")
        with page(browser) as pg:
            enter(pg, 4)
            pg.get_by_role("button", name="Generate simulated patient").click()
            check("the button confirms on the screen", wait_for(pg, "Done!", 30) is not None, text_of(pg)[-400:])
            check("and writes mutated_patient.fasta", os.path.exists("mutated_patient.fasta"))
            shot(pg, "screen4_patient")
        t = ingest(browser, ["mutated_patient.fasta"])
        check(f"mutated_patient.fasta -> mutated_patient_K{K}.arrow", number(t, "extracted:") == GENOME_LEN - K + 1, number(t, "extracted:"))

        print("== screen 3: assembly (De Bruijn)")
        with page(browser) as pg:
            enter(pg, 3)
            choose(pg, 0, f"ecoli_K{K}.arrow")
            submit(pg, "Build the De Bruijn graph")
            t = wait_for(pg, "Graph resolved successfully|ERROR", 300)
            shot(pg, "screen3_assembly")
            n = number(t, "Longest contig:")
            check(f"rebuilds the {GENOME_LEN} bp in a single contig", n == GENOME_LEN, n)
            contigs = [l.strip() for l in open(f"MMGEN_Contigs_K{K}.fasta") if not l.startswith(">")] if os.path.exists(f"MMGEN_Contigs_K{K}.fasta") else [""]
            longest = max(contigs, key=len)
            check("the contig written to disk is the original genome", longest == genome or "".join(COMP[b] for b in reversed(longest)) == genome, f"{len(longest)} bp")
            check("the contig table shows up", "NODE_1" in (t or ""))

        print("== screen 4: variant calling (Bloom)")
        with page(browser) as pg:
            enter(pg, 4)
            choose(pg, 0, f"ecoli_K{K}.arrow")
            choose(pg, 1, f"mutated_patient_K{K}.arrow")
            submit(pg, "Run the search")
            t = wait_for(pg, "Genome scan complete|ERROR", 300)
            shot(pg, "screen4_variants")
            n = number(t, "patient):")
            check(f"3 SNPs -> {3 * K} mutated k-mers", n == 3 * K, n)
            check(f"the panel says it ran on the {WHERE}, which is true", f"on the {WHERE}" in (t or ""), (t or "")[-300:])
        with page(browser) as pg:
            enter(pg, 4)
            choose(pg, 0, f"ecoli_K{K}.arrow")
            choose(pg, 1, f"ecoli_K{K}.arrow")
            submit(pg, "Run the search")
            check("the same file in both fields: a warning", wait_for(pg, "cannot be the same file", 30) is not None, text_of(pg)[-300:])

        print("== screen 5: annotation (translation)")
        with page(browser) as pg:
            enter(pg, 5)
            pg.locator("textarea:visible").first.fill("ATGGCCAAATAG")
            submit(pg, "Translate to protein")
            t = wait_for(pg, "Translation complete|ERROR", 60)
            shot(pg, "screen5_annotation")
            check("ATGGCCAAATAG -> MAK*, 4 amino acids", t is not None and "MAK" in t and "4 amino acids" in t, (t or "")[-400:])
            colour = pg.locator("b", has_text="*").first.evaluate("e => getComputedStyle(e).color") if pg.locator("b", has_text="*").count() else None
            check("the stop codon is drawn as a red asterisk", colour == "rgb(255, 0, 0)" and "MAK**" not in (t or ""), colour)

        print("== screen 6: Bloom filter")
        with page(browser) as pg:
            enter(pg, 6)
            choose(pg, 0, f"ecoli_K{K}.arrow")
            submit(pg, "Build the Bloom filter")
            t = wait_for(pg, "built successfully", 120)
            shot(pg, "screen6_bloom")
            n = number(t, "Bits set:")
            expected = 5 * (GENOME_LEN - K + 1)
            check(f"{n} bits set (5 per k-mer: close to {expected})", n is not None and 0.99 * expected <= n <= expected, n)

        print("== screen 7: scaffolding")
        with page(browser) as pg:
            enter(pg, 7)
            choose(pg, 0, f"MMGEN_Contigs_K{K}.fasta")
            submit(pg, "Close the genome")
            t = wait_for(pg, "Scaffolding complete|ERROR", 120)
            shot(pg, "screen7_scaffold")
            check("stitches the contigs", number(t, "continuous sequence:") == GENOME_LEN, number(t, "continuous sequence:"))

        print("== screen 8: VCF")
        with page(browser) as pg:
            enter(pg, 8)
            multiselect(pg, list(BATCHES))
            pg.locator("input[type='text']:visible").fill("Browser_POS.vcf")
            submit(pg, "Start the Multimerge engine")
            t = wait_for(pg, "VCF Multimerge complete|MULTIMERGE ENGINE ERROR", 120)
            shot(pg, "screen8_vcf")
            everything = [r for recs in BATCHES.values() for r in recs]
            ref = sorted(everything, key=lambda l: (l.split("\t")[0], int(l.split("\t")[1])))
            out = read_vcf("Browser_POS.vcf") if os.path.exists("Browser_POS.vcf") else []
            check("three files: the panel shows 6,000 variants", number(t, "Variants written:") == 6000, number(t, "Variants written:"))
            check("the output is identical to Python's stable sort", out == ref, len(out))

            choose(pg, 0, "QUAL")
            pg.locator("input[type='text']:visible").fill("Browser_QUAL.vcf")
            submit(pg, "Start the Multimerge engine")
            t = wait_for(pg, "Browser_QUAL.vcf", 120)
            quals = [int(l.split("\t")[5]) for l in read_vcf("Browser_QUAL.vcf")] if os.path.exists("Browser_QUAL.vcf") else []
            check("by QUAL: 6,000 variants in descending order", len(quals) == 6000 and quals == sorted(quals, reverse=True), len(quals))

            before = read_vcf("sample_a.vcf")
            pg.locator("input[type='text']:visible").fill("sample_a.vcf")
            submit(pg, "Start the Multimerge engine")
            check("output equal to an input: a warning", wait_for(pg, "cannot be one of the input files", 30) is not None, text_of(pg)[-300:])
            check("output equal to an input: not overwritten", read_vcf("sample_a.vcf") == before)

        print("== screen 9: sort & merge, index & query, intervals")
        with page(browser) as pg:
            enter(pg, 9)
            multiselect(pg, list(BATCHES))
            choose(pg, 1, "natural: chr1, chr2, chr10")
            pg.get_by_role("checkbox", name=re.compile("Also write the index")).click()
            submit(pg, "🚀 Sort")
            t = wait_for(pg, "Sorted|ERROR", 120)
            shot(pg, "screen9_sort")
            mmg.sort_genomic(list(BATCHES), "api_sorted.vcf.gz", chrom_order="natural")
            check("sort: 6,000 records, and the .tbi", number(t, "Records written:") == 6000 and os.path.exists("sample_a.sorted.vcf.gz.tbi"), (t or "")[-300:])
            check("sort: the same records as the direct call",
                  os.path.exists("sample_a.sorted.vcf.gz") and records("sample_a.sorted.vcf.gz") == records("api_sorted.vcf.gz"))
            check("sort: the Python call under the result", "mmg.sort_genomic(" in code_text(pg) and "index=True" in code_text(pg), code_text(pg))

            radio(pg, "Index & query a region")
            choose(pg, 0, "sample_a.sorted.vcf.gz")
            placeholder(pg, "chr1:1000000-2000000, chr2").fill("chr2:1000-30000")
            submit(pg, "Index and query")
            t = wait_for(pg, r"record\(s\) in chr2|ERROR", 60)
            shot(pg, "screen9_query")
            expected = [l for l in records("api_sorted.vcf.gz") if l.split("\t")[0] == "chr2" and 1000 <= int(l.split("\t")[1]) <= 30000]
            check(f"query: {len(expected)} records in chr2:1000-30000", f"{len(expected):,} record(s) in chr2:1000-30000" in (t or ""), (t or "")[-300:])

            radio(pg, "Interval operations (BED)")
            choose(pg, 0, "intersect: A against B (bedtools intersect)")
            pg.get_by_text("B: BED file(s):").first.wait_for(timeout=30000)
            multiselect(pg, ["peaks.bed"])
            multiselect(pg, ["genes.bed"])
            choose(pg, 1, "unique: the lines of A that touch B (-u)")
            submit(pg, "Run intersect")
            t = wait_for(pg, "Intersect done|ERROR", 60)
            shot(pg, "screen9_intersect")
            mmg.bed_intersect(["peaks.bed"], ["genes.bed"], "api_intersect.bed", mode="unique")
            check("intersect -u: the same lines as the direct call",
                  os.path.exists("peaks.intersect.bed") and records("peaks.intersect.bed") == records("api_intersect.bed"), (t or "")[-300:])

        print("== screen 10: validation, reference genomes")
        with page(browser) as pg:
            enter(pg, 10)
            choose(pg, 0, "broken.vcf")
            submit(pg, "Validate")
            t = wait_for(pg, r"error\(s\)|ERROR", 60)
            shot(pg, "screen10_validate")
            check("validate: the errors with their lines, and the screens that fix them",
                  t is not None and "line" in t and "Rename chromosomes" in t and "9. Sorting" in t, (t or "")[-400:])
            radio(pg, "Compare reference genomes")
            choose(pg, 0, "hg38.fai")
            choose(pg, 1, "hg19.fai")
            submit(pg, "Compare")
            t = wait_for(pg, "DIFFERENT reference genomes|ERROR", 60)
            check("compare: GRCh38 against GRCh37 is DIFFERENT", t is not None and "GRCh38" in t and "GRCh37" in t, (t or "")[-300:])

        print("== screen 11: VCF statistics, filter and normalization")
        with page(browser) as pg:
            enter(pg, 11)
            choose(pg, 0, "five_columns.vcf")
            submit(pg, "Count")
            t = wait_for(pg, "Nothing was read|ERROR", 60)
            shot(pg, "screen11_wrong_format")
            check("statistics of a file in the wrong format: says why nothing was read",
                  t is not None and "Nothing was read" in t and "needs at least 8" in t, (t or text_of(pg))[-400:])
            radio(pg, "Filter")
            choose(pg, 0, "sample_a.vcf")
            placeholder(pg, 'QUAL >= 30 && FILTER == "PASS"').fill("QUAL >= 50")
            submit(pg, "🔍 Filter")
            t = wait_for(pg, "Filtered|ERROR", 60)
            shot(pg, "screen11_filter")
            mmg.vcf_filter("sample_a.vcf", "api_filtered.vcf", expression="QUAL >= 50")
            check("filter: the same records as vcf_filter",
                  os.path.exists("sample_a.filtered.vcf") and records("sample_a.filtered.vcf") == records("api_filtered.vcf"), (t or "")[-300:])
            radio(pg, "Normalize")
            choose(pg, 0, "indel.vcf")
            choose(pg, 1, "small_genome.fa")
            submit(pg, "🧹 Normalize")
            t = wait_for(pg, "Normalized|ERROR", 60)
            shot(pg, "screen11_normalize")
            check("normalize with the genome: the deletion moves left (POS 101, CT>C)",
                  os.path.exists("indel.norm.vcf") and records("indel.norm.vcf") == ["chr1\t101\tdel\tCT\tC\t50\tPASS\t."], (t or "")[-300:])

        print("== screen 12: FASTQ report and trimming")
        with page(browser) as pg:
            enter(pg, 12)
            choose(pg, 0, "pair_R1.fastq.gz")
            submit(pg, "Measure")
            t = wait_for(pg, "Report page|ERROR", 120)
            frame = pg.frame_locator("iframe").first
            try:
                frame.locator("svg").first.wait_for(timeout=30000)
                charts = frame.locator("svg").count()
            except Exception:
                charts = 0
            shot(pg, "screen12_report")
            check(f"report: the adapter found, and the report's charts drawn inside the screen ({charts} svg)",
                  t is not None and "Illumina Universal" in t and charts >= 3, (t or "")[-300:])
            radio(pg, "Trim reads")
            choose(pg, 0, "pair_R1.fastq.gz")
            choose(pg, 1, "pair_R2.fastq.gz")
            submit(pg, "✂️ Trim")
            t = wait_for(pg, "Trimmed|ERROR", 120)
            shot(pg, "screen12_trim")
            mmg.fastq_trim("pair_R1.fastq.gz", "api_R1.fastq.gz", input2="pair_R2.fastq.gz", output2="api_R2.fastq.gz",
                           adapter="auto", quality_cutoff=20, min_length=20)
            check("trim, paired: R1 and R2 as fastq_trim writes them",
                  all(os.path.exists(n) for n in ("pair_R1.trimmed.fastq.gz", "pair_R2.trimmed.fastq.gz"))
                  and records("pair_R1.trimmed.fastq.gz") == records("api_R1.fastq.gz")
                  and records("pair_R2.trimmed.fastq.gz") == records("api_R2.fastq.gz"), (t or "")[-300:])

        print("== screen 13: k-mer counting")
        with page(browser) as pg:
            enter(pg, 13)
            multiselect(pg, ["pair_R1.fastq.gz", "pair_R2.fastq.gz"])
            placeholder(pg, "blank: the summary only, no table (for example kmers.tsv.gz)").fill("kmers.tsv.gz")
            submit(pg, "Count")
            t = wait_for(pg, "Counted|ERROR", 120)
            shot(pg, "screen13_kmers")
            api = mmg.count_kmers(["pair_R1.fastq.gz", "pair_R2.fastq.gz"], "api_kmers.tsv.gz", k=31)
            check("count: the totals, and the table as count_kmers writes it",
                  number(t, "Distinct k-mers:") == api["distinct_kmers"] and records("kmers.tsv.gz") == records("api_kmers.tsv.gz"), (t or "")[-300:])
            check("count: the spectrum chart is drawn", pg.locator("svg[aria-label^='Abundance spectrum']").count() >= 1)

        print("== screen 14: formats & conversion")
        with page(browser) as pg:
            enter(pg, 14)
            choose(pg, 0, "quals.vcf")
            choose(pg, 1, "Parquet (.parquet)")
            pg.get_by_role("radio", name="as text, exactly as written (a byte-for-byte round trip)").click()
            submit(pg, "Convert")
            t = wait_for(pg, "Converted|ERROR", 60)
            shot(pg, "screen14_convert")
            back_ok = False
            if os.path.exists("quals.vcf.parquet"):
                mmg.convert_file("quals.vcf.parquet", "quals_back.vcf")
                back_ok = open("quals_back.vcf", "rb").read() == open("quals.vcf", "rb").read()
            check("convert to Parquet with QUAL as text: back to text, the same bytes", back_ok, (t or "")[-300:])
            pg.get_by_role("button", name="Refresh the file list").click()
            pg.wait_for_timeout(1500)
            options = pg.locator("select:visible").first.locator("option").all_inner_texts()
            check("'Refresh the file list' offers the new file", "quals.vcf.parquet" in options, options)

        print("== screen 2: clustering")
        with page(browser) as pg:
            enter(pg, 2)
            choose(pg, 0, f"ecoli_K{K}.arrow")
            choose(pg, 1, f"mutated_patient_K{K}.arrow")
            submit(pg, "Load and merge the k-mers")
            t = wait_for(pg, "not accepted here yet|Engine finished", 60)
            shot(pg, "screen2_128bit_arrow")
            check("Module 1 .arrow (128 bits): refused with an explanation, instead of '2 k-mers' of garbage",
                  t is not None and "not accepted here yet" in t, (t or text_of(pg))[-300:])
        with page(browser) as pg:
            enter(pg, 2)
            choose(pg, 0, "legacy_A.arrow")
            choose(pg, 1, "legacy_B.arrow")
            submit(pg, "Load and merge the k-mers")
            t = wait_for(pg, "Engine finished|not accepted here yet", 120)
            check("old .arrow (one 64-bit k-mer per row): 3,200 k-mers", number(t, "Engine finished.") == 3200, number(t, "Engine finished."))
        with page(browser) as pg:
            enter(pg, 2)
            choose(pg, 0, "species_A.txt")
            choose(pg, 1, "species_B.txt")
            submit(pg, "Load and merge the k-mers")
            t = wait_for(pg, "Engine finished", 120)
            check("two .txt files of 31-mers: 4,000 k-mers", number(t, "Engine finished.") == 4000, number(t, "Engine finished."))
            submit(pg, "Cluster the k-mers")
            t = wait_for(pg, "Result:", 300)
            shot(pg, "screen2_clustering")
            check("HDBSCAN runs and the table shows up", t is not None and "Cluster_ID" in t, (t or text_of(pg))[-300:])
            # The table fetches each page from the server: it only works if the
            # component is still alive on the Python side.
            first = pg.locator("tbody tr").first.inner_text()
            nxt = next_page(pg)
            if nxt is not None:
                nxt.click()
                pg.wait_for_timeout(1500)
            check("the table changes page", nxt is not None and pg.locator("tbody tr").first.inner_text() != first,
                  first if nxt is not None else pg.evaluate(JS_BUTTONS))

        check("no JavaScript error on the page", not js_errors, js_errors[:3])
        browser.close()
except Exception as e:
    check("the script reached the end", False, repr(e))
    print(log_tail())
finally:
    shut_down()
    os.chdir(os.path.dirname(UI_PATH))
    shutil.rmtree(T, ignore_errors=True)

print()
print(f"{done - len(failures)} of {done} checks passed")
if SHOTS:
    print(f"screenshots in {SHOTS}")
if failures:
    print("FAILED:", failures)
    sys.exit(1)
