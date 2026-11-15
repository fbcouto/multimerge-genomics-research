"""Exercises for a class: every student gets data of their own, made from a
seed (the course and the student's name), and the teacher gets the answer key,
computed by MMGEN itself on each student's data.

    mmgen-exercises make --students students.txt --out class_folder [--course "Genomics 2026"]
                         [--exercises snps,broken_vcf,adapter,paternity,expression] [--zip]
    mmgen-exercises grade --key class_folder/_teacher/answer_key.csv --answers class_folder [--out grades.csv]

`students.txt` has one name per line. Each student gets a folder with the data,
a README.md with the questions and an answers.csv to fill in; the teacher gets
`_teacher/answer_key.csv`. The same course and name always give the same data,
so a lost folder can be made again. Every file is simulated: no real person or
organism. From Python: `make_class(...)` and `grade(...)`.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import random
import re
import sys
import unicodedata
import zipfile

from . import _native as mmg

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _rng(course: str, student: str, exercise: str) -> random.Random:
    digest = hashlib.sha256(f"{course}\x1f{student}\x1f{exercise}".encode("utf-8")).hexdigest()
    return random.Random(int(digest[:16], 16))


def _slug(name: str) -> str:
    s = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode().lower()
    s = re.sub(r"[^a-z0-9]+", "_", s).strip("_")
    return s or "student"


def _write(path: str, text: str) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)


def _write_gz(path: str, text: str) -> None:
    # mtime=0: the same data give the same bytes
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
        gz.write(text.encode("utf-8"))


def _fasta(name: str, seq: str) -> str:
    return f">{name}\n" + "".join(seq[i:i + 80] + "\n" for i in range(0, len(seq), 80))


def _dna(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def _q(qid, text, answer, kind="int", tolerance=0.0, accepted=()):
    return {"id": qid, "text": text, "answer": answer, "kind": kind, "tolerance": tolerance, "accepted": list(accepted)}


# ---------------------------------------------------------------------------
# The exercises: each writes its files into `folder` and returns
# (title, screens, intro, files, questions). The answers come from MMGEN.
# ---------------------------------------------------------------------------


def _ex_snps(rng, folder):
    n = rng.randint(2, 6)
    length = rng.randint(30, 50) * 1000
    genome = _dna(rng, length)
    patient = list(genome)
    for p in sorted(rng.sample(range(100, length - 100, 100), n)):
        patient[p] = rng.choice([b for b in "ACGT" if b != genome[p]])
    ref, pat = os.path.join(folder, "reference.fasta"), os.path.join(folder, "patient.fasta")
    _write(ref, _fasta("reference", genome))
    _write(pat, _fasta("patient", "".join(patient)))
    # as screen 4 counts them: the patient's k-mers (k = 31) that are not in the reference
    kr, kp = mmg.read_fasta_batch([ref], 31), mmg.read_fasta_batch([pat], 31)
    known = {kr[i:i + 16] for i in range(0, len(kr), 16)}
    novel = sum(1 for i in range(0, len(kp), 16) if kp[i:i + 16] not in known)
    return ("Find the planted SNPs", "1 and 4",
            "A healthy genome and a patient's copy of it, which carries a few single-base substitutions (SNPs).",
            ["reference.fasta", "patient.fasta"],
            [_q("snps_kmers", "Turn both files into k-mers on screen 1 (k = 31), then compare them on screen 4 "
                              "(reference: the healthy genome). How many of the patient's k-mers are not in the reference?", novel),
             _q("snps_count", "How many SNPs does the patient carry? (Each SNP changes k of the patient's k-mers.)", n)])


def _ex_broken_vcf(rng, folder):
    lens = {"chr1": 248956422, "chr2": 242193529}
    recs = []
    for chrom in ("chr1", "chr2"):
        for pos in sorted(rng.sample(range(10_000, 5_000_000), rng.randint(14, 20))):
            ref = rng.choice("ACGT")
            alt = rng.choice([b for b in "ACGT" if b != ref])
            recs.append([chrom, str(pos), f"rs{rng.randint(1000, 9_999_999)}", ref, alt, str(rng.randint(20, 99)), "PASS", "."])
    # three kinds of defect, each in a random place
    lines = ["\t".join(r) for r in recs]
    i = rng.randrange(3, len(lines) - 3)
    lines[i], lines[i + 1] = lines[i + 1], lines[i]                       # two records out of order
    j = rng.choice([k for k in range(len(lines)) if k not in (i, i + 1)])
    lines[j] = lines[j].replace("\t", " ")                                 # spaces instead of tabs
    k = rng.choice([k for k in range(len(lines)) if k not in (i, i + 1, j) and lines[k].startswith("chr2")])
    lines[k] = lines[k].replace("chr2", "2", 1)                            # another chromosome naming
    header = ["##fileformat=VCFv4.2"] + [f"##contig=<ID={c},length={n}>" for c, n in lens.items()] + \
             ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    path = os.path.join(folder, "lab_results.vcf")
    _write(path, "\n".join(header + lines) + "\n")
    rep = mmg.validate_file(path)
    problems = [p for p in rep["problems"] if p.get("line")]
    return ("What is wrong with this VCF?", "10",
            "A VCF that another program refuses. Find out why before you try to fix it.",
            ["lab_results.vcf"],
            [_q("vcf_kinds", "Validate the file on screen 10. How many kinds of problem does the validator report "
                             "(errors, warnings and notes together)?", len(rep["problems"])),
             _q("vcf_first_line", "On which line of the file is the first problem?", min(p["line"] for p in problems)),
             _q("vcf_errors", "How many of those kinds are errors (not warnings or notes)?",
                sum(1 for p in rep["problems"] if p.get("severity") == "error"))])


_ADAPTERS = {
    "Illumina TruSeq": ("AGATCGGAAGAGC", ("truseq", "illumina universal", "agatcggaagagc")),
    "Nextera": ("CTGTCTCTTATACACATCT", ("nextera", "ctgtctcttatacacatct")),
    "Illumina small RNA": ("TGGAATTCTCGGGTGCCAAGG", ("small rna", "smallrna", "tggaattctcgggtgccaagg")),
}


def _ex_adapter(rng, folder):
    name = rng.choice(sorted(_ADAPTERS))
    seq, aliases = _ADAPTERS[name]
    genome = _dna(rng, 20_000)
    out = []
    for r in range(2000):
        size = rng.randint(30, 140)  # the DNA fragment; shorter than the read, the read runs into the adapter
        a = rng.randrange(0, len(genome) - size)
        read = (genome[a:a + size] + seq + _dna(rng, 100))[:100]
        qual = ["I"] * 100
        if rng.random() < 0.25:  # a low-quality end
            for x in range(100 - rng.randint(5, 25), 100):
                qual[x] = "#"
        out.append(f"@read{r + 1}\n{read}\n+\n{''.join(qual)}\n")
    path = os.path.join(folder, "reads.fastq.gz")
    _write_gz(path, "".join(out))
    trimmed = os.path.join(folder, ".key_trimmed.fastq")
    t = mmg.fastq_trim(path, trimmed, adapter=seq, quality_cutoff=20, min_length=25)
    os.remove(trimmed)
    return ("Which adapter, and what is left after trimming?", "12",
            "2,000 reads of 100 bases. Many DNA fragments were shorter than the read, so the sequencer read on into "
            "the adapter; some reads also end in low-quality bases.",
            ["reads.fastq.gz"],
            [_q("adapter_name", "Make the quality report on screen 12: which adapter was used? (Its name or its "
                                "sequence.)", name, kind="text", accepted=aliases + (name.lower(),)),
             _q("adapter_reads_left", "Trim the reads on screen 12 with that adapter, a quality cutoff of 20 and a "
                                      "minimum length of 25. How many reads are left?", t["reads_out"], tolerance=0.02)])


def _ex_paternity(rng, folder):
    n = 4000
    freqs = [rng.uniform(0.05, 0.95) for _ in range(n)]

    def founder():
        return [tuple(int(rng.random() < p) for _ in range(2)) for p in freqs]

    def child(a, b):
        return [(x[rng.randrange(2)], y[rng.randrange(2)]) for x, y in zip(a, b)]

    grandfather, grandmother, mother, stranger = founder(), founder(), founder(), founder()
    father, uncle = child(grandfather, grandmother), child(grandfather, grandmother)
    kid = child(father, mother)
    roles = {"mother": mother, "child": kid, "father": father, "uncle": uncle, "stranger": stranger}
    labels = [f"person_{i}" for i in range(1, 6)]
    rng.shuffle(labels)
    name = dict(zip(roles, labels))
    order = sorted(roles, key=lambda r: name[r])
    path = os.path.join(folder, "family.vcf")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n" + "".join(f"##contig=<ID=chr{c},length=200000000>\n" for c in (1, 2, 3)) +
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(name[r] for r in order) + "\n")
        per = -(-n // 3)  # sites per chromosome, positions increasing within each
        for i in range(n):
            chrom, pos = 1 + i // per, 10_000 + (i % per) * 40_000
            ref, alt = rng.sample("ACGT", 2)
            fh.write(f"chr{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t" +
                     "\t".join(f"{roles[r][i][0]}|{roles[r][i][1]}" for r in order) + "\n")
    rel = mmg.relatedness([path])
    kin = {frozenset((p["a"], p["b"])): p["kinship"] for p in rel["pairs"]}
    k_father = kin[frozenset((name["child"], name["father"]))]
    candidates = sorted(name[r] for r in ("father", "uncle", "stranger"))
    return ("Who is the father?", "16",
            f"Five simulated people (no real person): {name['mother']} is the mother and {name['child']} the child; "
            f"{', '.join(candidates)} are the candidates. One of them is the father, another is his brother.",
            ["family.vcf"],
            [_q("father", "Measure the relatedness of everyone on screen 16. Which candidate is the child's father?",
                name["father"], kind="text"),
             _q("father_kinship", "What is the kinship coefficient between the child and the father? (Two decimals.)",
                round(k_father, 2), kind="float", tolerance=0.02),
             _q("uncle", "Which other candidate is related to the child? (A second-degree relative: the father's brother.)",
                name["uncle"], kind="text")])


def _ex_expression(rng, folder):
    genes = [f"GENE{i:04d}" for i in range(1, 2001)]
    changed = rng.sample(genes, rng.randint(4, 8))
    up = set(changed[:len(changed) // 2])
    samples = ["control_1", "control_2", "control_3", "treated_1", "treated_2", "treated_3"]
    depth = [rng.uniform(0.8, 1.25) for _ in samples]
    rows = []
    for g in genes:
        mu = rng.lognormvariate(5.5, 1.0)
        fold = rng.uniform(4, 8) if g in up else (1 / rng.uniform(4, 8) if g in changed else 1.0)
        counts = []
        for j, _s in enumerate(samples):
            m = mu * depth[j] * (fold if j >= 3 else 1.0)
            lam = rng.gammavariate(1 / 0.05, m * 0.05)  # negative binomial: gamma-Poisson, dispersion 0.05
            counts.append(max(0, int(round(rng.gauss(lam, lam ** 0.5)))) if lam > 30 else _poisson(rng, lam))
        rows.append(g + "\t" + "\t".join(map(str, counts)))
    cpath, spath = os.path.join(folder, "counts.tsv"), os.path.join(folder, "samples.tsv")
    _write(cpath, "gene\t" + "\t".join(samples) + "\n" + "\n".join(rows) + "\n")
    _write(spath, "sample\tcondition\n" + "".join(f"{s}\t{s.split('_')[0]}\n" for s in samples))
    out = os.path.join(folder, ".key_de.tsv")
    d = mmg.differential_expression(cpath, spath, "condition", "control", output=out, alpha=0.1)
    with open(out, encoding="utf-8") as fh:
        table = list(csv.DictReader(fh, delimiter="\t"))
    os.remove(out)
    sig = [r for r in table if r["padj"] not in ("", "NA", "nan") and float(r["padj"]) < 0.1]
    top = max(sig, key=lambda r: abs(float(r["log2FoldChange"])))["gene"] if sig else ""
    return ("Which genes changed?", "22",
            "Read counts of 2,000 genes in three control and three treated samples (simulated). A few genes respond "
            "to the treatment.",
            ["counts.tsv", "samples.tsv"],
            [_q("de_up", "Compare treated against control on screen 22 (task Differential expression; false-discovery "
                         "rate 0.1). How many genes go up?", d["up"]),
             _q("de_down", "How many genes go down?", d["down"]),
             _q("de_top", "Which gene changed the most (the largest absolute log2 fold change among the significant ones)?",
                top, kind="text")])


def _poisson(rng, lam):
    if lam <= 0:
        return 0
    import math
    limit, k, p = math.exp(-lam), 0, 1.0
    while True:
        p *= rng.random()
        if p <= limit:
            return k
        k += 1


EXERCISES = {
    "snps": _ex_snps,
    "broken_vcf": _ex_broken_vcf,
    "adapter": _ex_adapter,
    "paternity": _ex_paternity,
    "expression": _ex_expression,
}

# ---------------------------------------------------------------------------
# Making the class
# ---------------------------------------------------------------------------


def make_class(students, out: str, course: str = "MMGEN class", exercises=None, zip_each: bool = False) -> dict:
    """Writes one folder per student (data, README.md with the questions,
    answers.csv to fill in) and the teacher's answer key in `_teacher/`.
    `students`: a list of names, or the path of a text file with one per line."""
    if isinstance(students, str):
        with open(students, encoding="utf-8") as fh:
            students = [line.strip() for line in fh if line.strip() and not line.startswith("#")]
    exercises = list(exercises or EXERCISES)
    unknown = [e for e in exercises if e not in EXERCISES]
    if unknown:
        raise ValueError(f"unknown exercise(s): {', '.join(unknown)}; the exercises are: {', '.join(EXERCISES)}")
    if not students:
        raise ValueError("no student given")
    slugs, seen = {}, {}
    for s in students:
        base = _slug(s)
        seen[base] = seen.get(base, 0) + 1
        slugs[s] = base if seen[base] == 1 else f"{base}_{seen[base]}"
    os.makedirs(os.path.join(out, "_teacher"), exist_ok=True)
    key_rows = []
    for student in students:
        sdir = os.path.join(out, slugs[student])
        os.makedirs(sdir, exist_ok=True)
        parts = [f"# {course}: exercises for {student}\n",
                 "Every file here was made for you: your data are not the same as your colleagues'. Answer in "
                 "`answers.csv` (one answer per question, in the `answer` column) and hand it in. The data are "
                 "simulated: no real person or organism.\n",
                 "Open the files with MMGEN (`mmgen-ui` in this folder, or the Colab notebook); the screen to use is "
                 "given with each exercise.\n"]
        answer_rows = []
        for n, ex in enumerate(exercises, 1):
            edir = os.path.join(sdir, ex)
            os.makedirs(edir, exist_ok=True)
            title, screens, intro, files, questions = EXERCISES[ex](_rng(course, student, ex), edir)
            parts.append(f"\n## {n}. {title}\n\n{intro}\n\nFiles: " + ", ".join(f"`{ex}/{f}`" for f in files) +
                         f". Screen(s): {screens}.\n")
            for q in questions:
                parts.append(f"* **{q['id']}**: {q['text']}")
                answer_rows.append({"question": q["id"], "answer": ""})
                key_rows.append({"student": student, "folder": slugs[student], "exercise": ex, "question": q["id"],
                                 "answer": q["answer"], "kind": q["kind"], "tolerance": q["tolerance"],
                                 "accepted": "|".join(q["accepted"])})
        _write(os.path.join(sdir, "README.md"), "\n".join(parts) + "\n")
        with open(os.path.join(sdir, "answers.csv"), "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["question", "answer"])
            w.writeheader()
            w.writerows(answer_rows)
        if zip_each:
            with zipfile.ZipFile(os.path.join(out, slugs[student] + ".zip"), "w", zipfile.ZIP_DEFLATED) as z:
                for base, _dirs, names in os.walk(sdir):
                    for f in sorted(names):
                        full = os.path.join(base, f)
                        z.write(full, os.path.relpath(full, out))
    key_path = os.path.join(out, "_teacher", "answer_key.csv")
    with open(key_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(key_rows[0]))
        w.writeheader()
        w.writerows(key_rows)
    _write(os.path.join(out, "_teacher", "README.md"),
           f"# {course}: the teacher's folder\n\n"
           f"{len(students)} student(s), {len(exercises)} exercise(s): {', '.join(exercises)}.\n\n"
           "* Give each student their own folder (or the `.zip` made with `--zip`). Do not hand out this `_teacher/` folder.\n"
           "* `answer_key.csv` holds every student's answers, computed by MMGEN on that student's data.\n"
           "* Grade the `answers.csv` the students hand back, put back into their folders:\n\n"
           f"      mmgen-exercises grade --key {key_path} --answers {out}\n\n"
           "  or a single CSV with the columns `student,question,answer`.\n"
           "* The same course name and student name always give the same data: a lost folder can be made again.\n")
    return {"students": len(students), "exercises": exercises, "questions": len(key_rows), "folder": out, "key": key_path}


# ---------------------------------------------------------------------------
# Grading
# ---------------------------------------------------------------------------


def _correct(row: dict, given: str) -> bool:
    given = str(given if given is not None else "").strip()
    if not given:
        return False
    kind, expected = row["kind"], str(row["answer"])
    if kind in ("int", "float"):
        try:
            g, e = float(given.replace(",", ".")), float(expected)
        except ValueError:
            return False
        tol = float(row.get("tolerance") or 0)
        if kind == "float":
            return abs(g - e) <= tol
        return abs(g - e) <= tol * max(abs(e), 1) if tol else g == e
    norm = lambda x: re.sub(r"\s+", " ", x.strip().lower())
    accepted = [norm(expected)] + [norm(a) for a in (row.get("accepted") or "").split("|") if a]
    g = norm(given)
    return any(g == a or (len(a) >= 5 and a in g) for a in accepted)


def grade(key: str, answers: str, out: str | None = None) -> dict:
    """Grades the answers against the key. `answers`: the class folder (each
    student's answers.csv in their folder) or one CSV with student, question,
    answer. Writes `out` (one row per answer) when given; returns the scores."""
    with open(key, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    given = {}
    if os.path.isdir(answers):
        for folder in {r["folder"] for r in rows}:
            path = os.path.join(answers, folder, "answers.csv")
            if os.path.exists(path):
                with open(path, encoding="utf-8") as fh:
                    for a in csv.DictReader(fh):
                        given[(folder, a.get("question", "").strip())] = a.get("answer", "")
    else:
        names = {r["student"]: r["folder"] for r in rows}
        with open(answers, encoding="utf-8") as fh:
            for a in csv.DictReader(fh):
                folder = names.get(a.get("student", "").strip(), _slug(a.get("student", "")))
                given[(folder, a.get("question", "").strip())] = a.get("answer", "")
    results, scores = [], {}
    for r in rows:
        g = given.get((r["folder"], r["question"]), "")
        ok = _correct(r, g)
        results.append({"student": r["student"], "exercise": r["exercise"], "question": r["question"],
                        "expected": r["answer"], "given": g, "correct": ok})
        s = scores.setdefault(r["student"], [0, 0])
        s[0] += ok
        s[1] += 1
    if out:
        with open(out, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(results[0]))
            w.writeheader()
            w.writerows(results)
    return {"students": {k: {"correct": v[0], "questions": v[1], "percent": round(100 * v[0] / v[1], 1)}
                         for k, v in scores.items()},
            "answers": results}


# ---------------------------------------------------------------------------
# The command line
# ---------------------------------------------------------------------------


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="mmgen-exercises", description="Exercises with data of their own for every student, and the answer key.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("make", help="make the students' folders and the answer key")
    m.add_argument("--students", required=True, help="a text file with one student name per line")
    m.add_argument("--out", required=True, help="the class folder to write")
    m.add_argument("--course", default="MMGEN class", help="the course name (part of every seed: change it for a new class)")
    m.add_argument("--exercises", default=",".join(EXERCISES), help=f"comma-separated, among: {', '.join(EXERCISES)}")
    m.add_argument("--zip", action="store_true", help="also write one .zip per student, to hand out")
    g = sub.add_parser("grade", help="grade the answers against the key")
    g.add_argument("--key", required=True, help="_teacher/answer_key.csv")
    g.add_argument("--answers", required=True, help="the class folder (answers.csv in each student's folder) or one CSV with student,question,answer")
    g.add_argument("--out", help="write every answer, right or wrong, to this CSV")
    args = ap.parse_args(argv)
    if args.cmd == "make":
        r = make_class(args.students, args.out, args.course, [e.strip() for e in args.exercises.split(",") if e.strip()], args.zip)
        print(f"{r['students']} student folder(s), {r['questions']} question(s) in all, in {r['folder']}; answer key: {r['key']}")
    else:
        r = grade(args.key, args.answers, args.out)
        for s, v in r["students"].items():
            print(f"{s}: {v['correct']}/{v['questions']} ({v['percent']}%)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
