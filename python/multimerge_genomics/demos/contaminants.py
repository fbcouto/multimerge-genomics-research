"""Teaching demo for the MMGEN contaminant screen (contaminantes.rs).

write_demo(folder, research_data=None) writes a small, deterministic set of
files (about 4 MB; gzip with mtime=0, the same bytes on every run with the
same zlib):

    references/          the genomes to build the screen database from
    references.tsv       name, kind and file of each reference
    raw_material_R1.fq.gz, raw_material_R2.fq.gz
                         a "raw material" (a cell-culture harvest): host
                         reads with two agents spiked at a low level, a
                         reagent contaminant, a viral fragment piled on one
                         region (as a vector element would be), and the
                         process control
    reagent_blank_R1.fq.gz, reagent_blank_R2.fq.gz
                         the reagent blank: host traces, the reagent
                         contaminant and the process control
    spike_ins.tsv        what was spiked, and how much
    README.md            what to run and what to expect

Genomes. Without `research_data` every genome is synthetic (random
sequence with the GC content of the kind it stands for), named
"synthetic ..." so that nobody takes it for a real organism. With the path
of MMGEN's research_data folder, real pieces are used where that folder has
them: 2 Mb of human chr21 (GRCh38, from 30,000,000) and the whole chrM as
the host, the human adenovirus C and hepatitis B virus genomes, the first
300 kb of Pseudomonas aeruginosa PAO1 (the reagent contaminant) and of
Escherichia coli K-12 MG1655, and the first 230 kb of Saccharomyces
cerevisiae chromosome I. The mycoplasma and the process control are always
synthetic (research_data has no mycoplasma genome).

Python standard library only. The random numbers come from SplitMix64,
written here, so the files do not depend on the version of `random`.

For research and teaching: this is simulated data, and the screen is not a
validated method (see README.md).
"""

import gzip
import hashlib
import io
import os

DEMO_VERSION = 1

MASK = (1 << 64) - 1

NOTICE = (
    "For research and teaching only. The screen is not a validated method: releasing a GMP batch "
    "requires a validated method run with validated software (ICH Q5A(R2) admits next-generation "
    "sequencing for adventitious viruses, but the validation is the laboratory's)."
)


class SplitMix64:
    def __init__(self, seed):
        self.s = seed & MASK

    def next(self):
        self.s = (self.s + 0x9E3779B97F4A7C15) & MASK
        z = self.s
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK
        return z ^ (z >> 31)

    def uniform(self):
        return (self.next() >> 11) / float(1 << 53)

    def below(self, n):
        return int(self.uniform() * n) % max(n, 1)


COMPLEMENT = bytes.maketrans(b"ACGTacgtN", b"TGCAtgcaN")


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def random_dna(n, gc, rng):
    """n bases with GC fraction `gc`: SHAKE-256 of a seed drawn from `rng`
    gives n bytes, and each byte picks a base through a 256-entry table."""
    g = int(round(gc * 128))
    table = bytes([ord("G")] * g + [ord("C")] * g + [ord("A")] * (128 - g) + [ord("T")] * (128 - g))
    seed = rng.next().to_bytes(8, "little")
    return hashlib.shake_256(b"mmgen-demo" + seed).digest(n).translate(table)


def mutate(seq, rate, rng):
    s = bytearray(seq)
    i = _gap(rate, rng)
    while i < len(s):
        s[i] = _other(s[i], rng)
        i += 1 + _gap(rate, rng)
    return bytes(s)


def _other(b, rng):
    while True:
        x = b"ACGT"[rng.below(4)]
        if x != b:
            return x


def _gap(rate, rng):
    # Geometric number of error-free bases before the next error.
    import math

    u = max(rng.uniform(), 1e-300)
    return int(math.log(u) / math.log(1.0 - rate))


def add_errors(read, rng):
    """Illumina-like substitutions: 0.3% over the first two thirds of the
    read, 1.5% over the last third; the qualities say so."""
    n = len(read)
    cut = 2 * n // 3
    s = bytearray(read)
    for lo, hi, rate in ((0, cut, 0.003), (cut, n, 0.015)):
        i = lo + _gap(rate, rng)
        while i < hi:
            s[i] = _other(s[i], rng)
            i += 1 + _gap(rate, rng)
    qual = b"I" * cut + b"5" * (n - cut)  # Q40 then Q20
    return bytes(s), qual


def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    recs, name, seq = [], None, []
    with opener(path, "rb") as f:
        for line in f:
            line = line.strip()
            if line.startswith(b">"):
                if name is not None:
                    recs.append((name, b"".join(seq)))
                name, seq = line[1:].decode(), []
            elif line:
                seq.append(line)
    if name is not None:
        recs.append((name, b"".join(seq)))
    return recs


def fasta_bytes(records):
    out = io.BytesIO()
    for name, seq in records:
        out.write(b">" + name.encode() + b"\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i : i + 80] + b"\n")
    return out.getvalue()


def write_gz(path, data):
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=1) as z:
            z.write(data)


class Genome:
    def __init__(self, records):
        self.recs = [s.upper() for _, s in records]
        self.cum = [0]
        for s in self.recs:
            self.cum.append(self.cum[-1] + len(s))

    def pair(self, rng, read_len=150, insert=330, region=None):
        for _ in range(1000):
            ins = insert + rng.below(81) - 40
            if region is None:
                x = rng.below(self.cum[-1])
                r = 0
                while self.cum[r + 1] <= x:
                    r += 1
                s, e = 0, len(self.recs[r])
            else:
                r, s, e = region
            ins = min(ins, e - s)
            if ins < read_len:
                continue
            st = s + rng.below(e - s - ins + 1)
            frag = self.recs[r][st : st + ins]
            if frag.strip(b"ACGT"):
                continue
            if rng.below(2):
                frag = revcomp(frag)
            a, qa = add_errors(frag[:read_len], rng)
            b, qb = add_errors(revcomp(frag)[:read_len], rng)
            return a, qa, b, qb
        raise ValueError("no fragment without N")


def write_pairs(folder, name, parts, seed):
    """parts: (tag, genome, pairs, region). The pairs are interleaved at
    random, as in a run."""
    rng = SplitMix64(seed)
    left = [p[2] for p in parts]
    r1, r2 = io.BytesIO(), io.BytesIO()
    n = 0
    total = sum(left)
    while total:
        x = rng.below(total)
        i = 0
        while x >= left[i]:
            x -= left[i]
            i += 1
        tag, genome, _, region = parts[i]
        a, qa, b, qb = genome.pair(rng, region=region)
        n += 1
        head = "%s:%d:%s" % (name, n, tag)
        r1.write(b"@" + head.encode() + b" 1:N:0\n" + a + b"\n+\n" + qa + b"\n")
        r2.write(b"@" + head.encode() + b" 2:N:0\n" + b + b"\n+\n" + qb + b"\n")
        left[i] -= 1
        total -= 1
    write_gz(os.path.join(folder, name + "_R1.fq.gz"), r1.getvalue())
    write_gz(os.path.join(folder, name + "_R2.fq.gz"), r2.getvalue())
    return n


def write_demo(folder, research_data=None):
    """Writes the demo into `folder` (created if needed). Returns a dict with
    the paths written and what the screen should find."""
    os.makedirs(os.path.join(folder, "references"), exist_ok=True)
    rng = SplitMix64(20261116)
    real = research_data is not None

    def rd(*parts):
        return os.path.join(research_data, *parts)

    refs = []  # (name, kind, file name, records)
    if real:
        chr21 = read_fasta(rd("metagenome", "host_decoy", "GRCh38.p14_chr21.fna.gz"))[0][1]
        chrm = read_fasta(rd("metagenome", "host_decoy", "GRCh38.p14_chrMT.fna.gz"))[0][1]
        host = [("chr21_30000001_32000000 Homo sapiens GRCh38 chr21, 2 Mb piece", chr21[30_000_000:32_000_000].upper()), ("chrM Homo sapiens mitochondrion", chrm.upper())]
        host_name = "Human host (2 Mb of chr21 + chrM)"
        vdir = rd("reference_genomes", "viruses")
        adv = read_fasta(os.path.join(vdir, "GCF_000845085.1_ViralProj14518_genomic.fna.gz"))
        hbv = read_fasta(os.path.join(vdir, "GCF_000861825.2_ViralProj15428_genomic.fna.gz"))
        pae = read_fasta(rd("reference_genomes", "bacteria", "GCF_000006765.1_ASM676v1_genomic.fna.gz"))[0]
        eco = read_fasta(rd("reference_genomes", "bacteria", "GCF_000005845.2_ASM584v2_genomic.fna.gz"))[0]
        yea = read_fasta(rd("reference_genomes", "fungi", "GCF_000146045.2_R64_genomic.fna.gz"))[0]
        virus1 = ("Human mastadenovirus C", adv)
        virus2 = ("Hepatitis B virus", hbv)
        kitome = ("Pseudomonas aeruginosa PAO1 (first 300 kb)", [(pae[0] + " (first 300 kb)", pae[1][:300_000])])
        bact2 = ("Escherichia coli K-12 MG1655 (first 300 kb)", [(eco[0] + " (first 300 kb)", eco[1][:300_000])])
        fungus = ("Saccharomyces cerevisiae chrI (first 230 kb)", [(yea[0], yea[1][:230_000])])
    else:
        hs = random_dna(2_000_000, 0.41, rng)
        element = random_dna(3000, 0.42, rng)
        hs = bytearray(hs)
        for j in range(25):
            p = 20_000 + j * 75_000
            hs[p : p + 3000] = mutate(element, 0.05, rng)
        host = [("synthetic_host_1 synthetic host genome (GC 41%, with a repeated element)", bytes(hs)), ("synthetic_host_mito synthetic mitochondrion", random_dna(16_500, 0.44, rng))]
        host_name = "Synthetic host"
        virus1 = ("Synthetic DNA virus (36 kb)", [("synthetic_dna_virus synthetic DNA virus genome", random_dna(36_000, 0.55, rng))])
        virus2 = ("Synthetic small virus (3.2 kb)", [("synthetic_small_virus synthetic small virus genome", random_dna(3_200, 0.48, rng))])
        kitome = ("Synthetic reagent bacterium (300 kb)", [("synthetic_reagent_bacterium synthetic GC-rich bacterium", random_dna(300_000, 0.66, rng))])
        bact2 = ("Synthetic enteric bacterium (300 kb)", [("synthetic_enteric_bacterium synthetic bacterium", random_dna(300_000, 0.51, rng))])
        fungus = ("Synthetic yeast (230 kb)", [("synthetic_yeast synthetic yeast chromosome", random_dna(230_000, 0.38, rng))])
    myco = ("Synthetic mycoplasma (250 kb, GC 28%)", [("synthetic_mycoplasma synthetic mycoplasma-like genome (GC 28%)", random_dna(250_000, 0.28, rng))])
    control = ("Synthetic process control (5 kb)", [("synthetic_process_control synthetic process-control sequence", random_dna(5_000, 0.50, rng))])

    refs = [
        (host_name, "host", "host.fa.gz", host),
        (virus1[0], "virus", "virus_1.fa.gz", virus1[1]),
        (virus2[0], "virus", "virus_2.fa.gz", virus2[1]),
        (kitome[0], "bacterium", "bacterium_reagent.fa.gz", kitome[1]),
        (bact2[0], "bacterium", "bacterium_2.fa.gz", bact2[1]),
        (myco[0], "mycoplasma", "mycoplasma.fa.gz", myco[1]),
        (fungus[0], "fungus", "fungus.fa.gz", fungus[1]),
        (control[0], "spike", "process_control.fa.gz", control[1]),
    ]
    lines = ["name\tkind\tfile"]
    for name, kind, fname, recs in refs:
        write_gz(os.path.join(folder, "references", fname), fasta_bytes(recs))
        lines.append("%s\t%s\treferences/%s" % (name, kind, fname))
    with open(os.path.join(folder, "references.tsv"), "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")

    g = {name: Genome(recs) for name, _, _, recs in refs}
    hostg, v1, v2, kit, my, ctl = (g[host_name], g[virus1[0]], g[virus2[0]], g[kitome[0]], g[myco[0]], g[control[0]])
    # The virus-2 piece that stands for a vector element: 500 bp.
    piece = (0, 1000, 1500)
    raw = [
        ("host", hostg, 20_000, None),
        ("virus1", v1, 15, None),
        ("mycoplasma", my, 12, None),
        ("reagent", kit, 4, None),
        ("vector_piece", v2, 40, piece),
        ("control", ctl, 50, None),
    ]
    blank = [("host", hostg, 2_900, None), ("reagent", kit, 50, None), ("control", ctl, 50, None)]
    n_raw = write_pairs(folder, "raw_material", raw, 1)
    n_blank = write_pairs(folder, "reagent_blank", blank, 2)
    with open(os.path.join(folder, "spike_ins.tsv"), "w", newline="\n") as f:
        f.write("sample\torganism\tamount\tunit\n")
        f.write("raw_material\t%s\t15\treads\n" % virus1[0])
        f.write("raw_material\t%s\t12\treads\n" % myco[0])
        f.write("*\t%s\t50\treads\n" % control[0])
    readme = README % {
        "notice": NOTICE,
        "source": "real pieces of research_data where available (see references.tsv)" if real else "synthetic genomes only (random sequence)",
        "virus1": virus1[0],
        "virus2": virus2[0],
        "myco": myco[0],
        "kitome": kitome[0],
        "control": control[0],
        "n_raw": n_raw,
        "n_blank": n_blank,
    }
    with open(os.path.join(folder, "README.md"), "w", newline="\n") as f:
        f.write(readme)
    return {
        "references": os.path.join(folder, "references.tsv"),
        "samples": {
            "raw_material": [os.path.join(folder, "raw_material_R1.fq.gz"), os.path.join(folder, "raw_material_R2.fq.gz")],
            "reagent_blank": [os.path.join(folder, "reagent_blank_R1.fq.gz"), os.path.join(folder, "reagent_blank_R2.fq.gz")],
        },
        "spike_ins": os.path.join(folder, "spike_ins.tsv"),
        "negative_controls": ["reagent_blank"],
        "expected": {
            "raw_material": {
                virus1[0]: "detected",
                myco[0]: "detected",
                kitome[0]: "background",
                virus2[0]: "below_threshold",
                control[0]: "detected",
            },
            "reagent_blank": {kitome[0]: "detected", control[0]: "detected"},
        },
    }


README = """# Contaminant screen: teaching demo

%(notice)s

Genomes: %(source)s. All the reads are simulated (2 x 150, Illumina-like
errors).

## Files

- `references/` and `references.tsv`: the genomes of the screen database,
  with their kind (host, virus, bacterium, mycoplasma, fungus, spike = the
  process control).
- `raw_material_R1.fq.gz`, `raw_material_R2.fq.gz`: %(n_raw)s read pairs of
  a "raw material" (say, a cell-culture harvest): mostly host, plus
  - %(virus1)s, 15 pairs (spiked);
  - %(myco)s, 12 pairs (spiked);
  - %(kitome)s, 4 pairs (a reagent contaminant);
  - 40 pairs from a 500-bp piece of %(virus2)s only (as a vector element or
    an integrated fragment would give);
  - %(control)s, 50 pairs (the process control).
- `reagent_blank_R1.fq.gz`, `reagent_blank_R2.fq.gz`: %(n_blank)s pairs of the
  reagent blank: host traces, the reagent contaminant (50 pairs) and the
  process control (50 pairs).
- `spike_ins.tsv`: what was spiked and how much.

## Run it

```python
import csv, multimerge_genomics as mmg

refs = list(csv.DictReader(open("references.tsv"), delimiter="\\t"))
mmg.build_screen_database([r["file"] for r in refs], "demo.mmgscreen",
                          kinds=[r["kind"] for r in refs], names=[r["name"] for r in refs])
report = mmg.contaminant_screen(
    [("raw_material", "raw_material_R1.fq.gz", "raw_material_R2.fq.gz"),
     ("reagent_blank", "reagent_blank_R1.fq.gz", "reagent_blank_R2.fq.gz")],
    "demo.mmgscreen", negative_controls=["reagent_blank"], spike_ins="spike_ins.tsv",
    output_html="screen.html", output_tsv="screen.tsv", output_json="screen.json")
print(report["text"])
```

## What to expect

In `raw_material`:

- **%(virus1)s: detected**, and **%(myco)s: detected**: few reads, but
  specific to them and spread over their genomes;
- **%(kitome)s: background**: it is also in the reagent blank, at a much
  higher rate (reads per million): the sample is not 10 times above the
  blank;
- **%(virus2)s: below threshold**: 40 pairs, all on the same 500 bp: the
  reads are concentrated in a few genome windows where reads of a real
  infection would be spread over the whole genome;
- %(control)s: the process control, recovered at about 100%%.

In `reagent_blank`: the reagent contaminant and the process control.

Open `screen.html` for the report: the calls, the evidence of each
organism (unique and tied reads, reads per million, distinct k-mers,
genome windows, dispersion), and the traceability (the SHA-256 of every
file and of the database, the parameters, the date, the engine version).
"""


if __name__ == "__main__":
    import sys

    out = write_demo(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
    print(out["expected"])
