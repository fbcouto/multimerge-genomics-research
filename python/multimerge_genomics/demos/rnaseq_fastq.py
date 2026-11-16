"""A small, deterministic teaching demo for RNA-seq quantification from FASTQ
(build_transcript_index -> quantify_samples -> differential_expression).

Standard library only. `write_demo(folder)` writes, in `folder`:

  transcripts.fa          269 synthetic transcripts of 140 genes (Ensembl-style
                          headers: `gene:` and `gene_symbol:`, so no gene map
                          is needed); many genes have several isoforms that
                          share exons (exon skipping), two genes are paralogs
                          (2% divergent copies of another gene)
  gene_map.tsv            the same transcript -> gene -> name map, as a table
                          (for tools that need one)
  control_1_R1.fastq.gz   paired-end reads, 75 bp, 3 control + 3 treated
  ... treated_3_R2.fastq.gz   samples, ~20,000 pairs each, fragments of
                          200 +/- 25 bp, 0.2% substitution errors
  samples.tsv             sample -> condition (control / treated)
  truth.tsv               the abundance drawn for every transcript in every
                          sample, and the fragments simulated from it
  README.md               what is in the folder and what to expect

Every run writes the same bytes (fixed seed; gzip with no file name and
mtime=0). About 7 MB in all.

What to expect (the planted changes):
  - Gene7, Gene19, Gene33 go UP 4x in treated (log2 fold change ~ +2);
  - Gene12, Gene48, Gene71 go DOWN 4x (log2 fold change ~ -2);
  - Gene90 SWITCHES ISOFORMS: its two isoforms use mutually exclusive
    middle exons of the same length; 80% of its molecules are isoform -201
    in control and 20% in treated (the reverse for -202); the gene's total
    does not change: differential expression at the gene level finds nothing
    there, the transcript table shows the switch;
  - everything else changes only by biological noise (about 10%).
With the default settings, differential_expression (condition "condition",
reference "control") on the gene_counts.tsv written by quantify_samples
should call the six planted genes (padj < 0.1, the right sign) and no
other; nearly every pair (over 99.9%) pseudoaligns; the estimated fragment
length is close to 200.
"""

import gzip
import math
import os
import random

DEMO_VERSION = 1

SEED = 2026
GENES = 140
READ_LEN = 75
PAIRS = 20_000
FRAG_MEAN = 200.0
FRAG_SD = 25.0
ERROR = 0.002
UP = (7, 19, 33)
DOWN = (12, 48, 71)
SWITCH = 90
PARALOGS = {41: 40, 101: 100}  # gene -> the gene it is a 2%-divergent copy of
SAMPLES = ["control_1", "control_2", "control_3", "treated_1", "treated_2", "treated_3"]

_COMP = str.maketrans("ACGT", "TGCA")


def _rc(s):
    return s.translate(_COMP)[::-1]


def _random_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _mutate(rng, s, rate):
    out = list(s)
    for i, c in enumerate(out):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != c])
    return "".join(out)


def _transcriptome(rng):
    """[(transcript, gene id, gene name, sequence)] and the transcripts of each gene."""
    exons_of = {}
    txs = []
    for g in range(1, GENES + 1):
        if g in PARALOGS:
            exons = [_mutate(rng, e, 0.02) for e in exons_of[PARALOGS[g]]]
        elif g in UP + DOWN + (SWITCH,):
            # planted genes: long enough for 200-bp fragments to cover them well
            n_ex = rng.randint(3, 6)
            exons = [_random_seq(rng, rng.randint(150, 400)) for _ in range(n_ex)]
        else:
            n_ex = rng.randint(1, 6)
            exons = [_random_seq(rng, rng.randint(80, 400)) for _ in range(n_ex)]
        exons_of[g] = exons
        n_ex = len(exons)
        isoforms = [tuple(range(n_ex))]
        if g == SWITCH:
            # two isoforms with mutually exclusive middle exons of the same
            # length: the switch changes which exon is used, not the length
            n_mid = rng.randint(150, 300)
            exons = [
                _random_seq(rng, rng.randint(200, 400)),
                _random_seq(rng, n_mid),
                _random_seq(rng, n_mid),
                _random_seq(rng, rng.randint(200, 400)),
            ]
            exons_of[g] = exons
            n_ex = 4
            isoforms = [(0, 1, 3), (0, 2, 3)]
        elif n_ex >= 3:
            for _ in range(rng.randint(1, 4)):
                iso = tuple(i for i in range(n_ex) if i in (0, n_ex - 1) or rng.random() > 0.4)
                if iso not in isoforms:
                    isoforms.append(iso)
        for k, iso in enumerate(isoforms):
            seq = "".join(exons[i] for i in iso)
            txs.append((f"DEMO{g:04d}-{201 + k}", f"DEMOG{g:04d}", f"Gene{g}", seq))
    return txs


def _eff_len(length):
    num = den = 0.0
    for f in range(1, min(length, 999) + 1):
        p = math.exp(-((f - FRAG_MEAN) ** 2) / (2 * FRAG_SD * FRAG_SD))
        num += p * (length - f + 1)
        den += p
    return num / den if den > 0 else 1.0


def _gzip_writer(path):
    raw = open(path, "wb")
    return raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=6)


def _errors(rng, s):
    """Substitutions at rate ERROR, jumping from one error to the next."""
    if ERROR <= 0:
        return s
    out = None
    i = int(math.log(1.0 - rng.random()) / math.log(1.0 - ERROR))
    while i < len(s):
        if out is None:
            out = list(s)
        out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
        i += 1 + int(math.log(1.0 - rng.random()) / math.log(1.0 - ERROR))
    return "".join(out) if out is not None else s


def write_demo(folder):
    """Writes the demo into `folder` (created if needed). Returns the list of
    files written."""
    os.makedirs(folder, exist_ok=True)
    rng = random.Random(SEED)
    txs = _transcriptome(rng)
    written = []

    path = os.path.join(folder, "transcripts.fa")
    with open(path, "w", newline="\n") as f:
        for name, gid, gname, seq in txs:
            f.write(f">{name} cdna synthetic gene:{gid} gene_symbol:{gname}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")
    written.append(path)
    path = os.path.join(folder, "gene_map.tsv")
    with open(path, "w", newline="\n") as f:
        f.write("transcript\tgene\tgene_name\n")
        for name, gid, gname, _ in txs:
            f.write(f"{name}\t{gid}\t{gname}\n")
    written.append(path)

    # abundances: log-normal per transcript; planted changes per gene
    base = [math.exp(rng.gauss(3.0, 1.2)) for _ in txs]
    gene_num = [int(gid[5:]) for _, gid, _, _ in txs]
    for i, g in enumerate(gene_num):
        if g in UP + DOWN + (SWITCH,):
            base[i] = max(base[i], 60.0)  # planted genes well expressed
    eff = [_eff_len(len(t[3])) for t in txs]
    switch_tx = [i for i, g in enumerate(gene_num) if g == SWITCH]

    truth_cols = []
    for s, sample in enumerate(SAMPLES):
        treated = sample.startswith("treated")
        noise = {g: math.exp(rng.gauss(0.0, 0.1)) for g in set(gene_num)}
        theta = []
        for i, g in enumerate(gene_num):
            x = base[i] * noise[g]
            if treated and g in UP:
                x *= 4.0
            elif treated and g in DOWN:
                x *= 0.25
            theta.append(x)
        if switch_tx:
            total = sum(base[i] for i in switch_tx) * noise[SWITCH]
            a, b = (0.2, 0.8) if treated else (0.8, 0.2)
            theta[switch_tx[0]] = total * a
            theta[switch_tx[1]] = total * b
        weights = [th * e if len(t[3]) >= 31 else 0.0 for th, e, t in zip(theta, eff, txs)]
        cum = []
        acc = 0.0
        for w in weights:
            acc += w
            cum.append(acc)
        n_pairs = int(PAIRS * (0.9 + 0.2 * rng.random()))
        counts = [0] * len(txs)
        p1 = os.path.join(folder, f"{sample}_R1.fastq.gz")
        p2 = os.path.join(folder, f"{sample}_R2.fastq.gz")
        raw1, w1 = _gzip_writer(p1)
        raw2, w2 = _gzip_writer(p2)
        qual = "I" * READ_LEN
        for r in range(n_pairs):
            u = rng.random() * acc
            lo, hi = 0, len(cum) - 1
            while lo < hi:
                mid = (lo + hi) // 2
                if cum[mid] <= u:
                    lo = mid + 1
                else:
                    hi = mid
            t = lo
            seq = txs[t][3]
            length = len(seq)
            f = 0
            for _ in range(100):
                x = int(round(rng.gauss(FRAG_MEAN, FRAG_SD)))
                if 1 <= x <= length and x < 1000:
                    f = x
                    break
            if f == 0:
                f = min(length, 999)
            start = rng.randint(0, length - f)
            frag = seq[start:start + f]
            if rng.random() < 0.5:
                frag = _rc(frag)
            n = min(READ_LEN, f)
            a = _errors(rng, frag[:n])
            b = _errors(rng, _rc(frag)[:n])
            w1.write(f"@{sample}.{r + 1}/1\n{a}\n+\n{qual[:n]}\n".encode())
            w2.write(f"@{sample}.{r + 1}/2\n{b}\n+\n{qual[:n]}\n".encode())
            counts[t] += 1
        for w, raw in ((w1, raw1), (w2, raw2)):
            w.close()
            raw.close()
        written += [p1, p2]
        tot = sum(th for th, t in zip(theta, txs) if len(t[3]) >= 31)
        truth_cols.append(([th / tot * 1e6 if len(t[3]) >= 31 else 0.0 for th, t in zip(theta, txs)], counts))

    path = os.path.join(folder, "samples.tsv")
    with open(path, "w", newline="\n") as f:
        f.write("sample\tcondition\n")
        for sample in SAMPLES:
            f.write(f"{sample}\t{'treated' if sample.startswith('treated') else 'control'}\n")
    written.append(path)

    path = os.path.join(folder, "truth.tsv")
    with open(path, "w", newline="\n") as f:
        head = ["transcript", "gene", "gene_name", "length", "planted"]
        for sample in SAMPLES:
            head += [f"{sample}_true_tpm", f"{sample}_fragments"]
        f.write("\t".join(head) + "\n")
        for i, (name, gid, gname, seq) in enumerate(txs):
            g = gene_num[i]
            planted = "up_4x" if g in UP else "down_4x" if g in DOWN else "isoform_switch" if g == SWITCH else ""
            row = [name, gid, gname, str(len(seq)), planted]
            for tpm, counts in truth_cols:
                row += [f"{tpm[i]:.3f}", str(counts[i])]
            f.write("\t".join(row) + "\n")
    written.append(path)

    n_genes_multi = sum(1 for g in set(gene_num) if gene_num.count(g) > 1)
    path = os.path.join(folder, "README.md")
    with open(path, "w", newline="\n") as f:
        f.write(
            "# RNA-seq quantification demo (synthetic data)\n\n"
            f"- `transcripts.fa`: {len(txs)} transcripts of {GENES} genes ({n_genes_multi} genes with several isoforms "
            "sharing exons; Gene41 and Gene101 are 2%-divergent copies of Gene40 and Gene100). The headers name the gene "
            "(`gene:` and `gene_symbol:`), so no gene map is needed; `gene_map.tsv` holds the same map as a table.\n"
            f"- Six paired-end samples (`*_R1.fastq.gz`, `*_R2.fastq.gz`): {READ_LEN}-bp reads, about {PAIRS:,} pairs each, "
            f"fragments of {FRAG_MEAN:.0f} +/- {FRAG_SD:.0f} bp, {ERROR * 100:.1f}% sequencing errors.\n"
            "- `samples.tsv`: the condition of each sample. `truth.tsv`: the abundance drawn for each transcript and the "
            "fragments simulated from it, per sample.\n\n"
            "## Steps\n\n"
            "1. `build_transcript_index(\"transcripts.fa\", \"demo.idx\")`\n"
            "2. `quantify_samples(\"demo.idx\", samples, \"quant\")` with the six samples (name, R1, R2, condition)\n"
            "3. `differential_expression(\"quant/gene_counts.tsv\", \"quant/samples.tsv\", \"condition\", \"control\")`\n\n"
            "## What to expect\n\n"
            "- Nearly every pair (over 99.9%) pseudoaligns: the reads come only from these transcripts. The fragment "
            "length is estimated close to 200 bp.\n"
            "- Gene7, Gene19 and Gene33 go up 4x in treated (log2 fold change near +2); Gene12, Gene48 and Gene71 go down 4x "
            "(near -2). These six should be the significant genes (padj < 0.1).\n"
            "- Gene90 switches isoforms: its two isoforms use mutually exclusive middle exons of the same length; "
            "DEMO0090-201 is 80% of the gene in control and 20% in treated, DEMO0090-202 the reverse. The gene's total "
            "does not change, so the gene-level test does not find it; compare the two isoforms in "
            "`quant/transcript_tpm.tsv`. (Had the isoforms differed in length, the gene's read count would have changed "
            "with the switch even at a constant number of molecules: reads are proportional to molecules x length.)\n"
            "- Estimated counts can be compared with the truth in `truth.tsv`: close at the gene level; at the isoform "
            "level, less so for isoforms that differ by little sequence.\n\n"
            "All sequences and reads are simulated; they do not come from any organism.\n"
        )
    written.append(path)
    return written


if __name__ == "__main__":
    import sys

    out = sys.argv[1] if len(sys.argv) > 1 else "quant_demo"
    files = write_demo(out)
    total = sum(os.path.getsize(p) for p in files)
    print(f"{len(files)} files, {total / 1e6:.1f} MB in {out}")
