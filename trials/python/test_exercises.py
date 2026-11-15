# Test of the class exercises (mmgen-exercises): a class is made, every
# exercise is solved the way a student would solve it on the screens (screen
# 4's Bloom filter, screen 10's validator, screen 12's report and trimming,
# screen 16's relatedness, screen 22's differential expression), the answers
# are graded against the key, and the data are checked to be the same for the
# same student and different between students.
# Usage: python trials/python/test_exercises.py
import csv
import filecmp
import os
import sys
import tempfile

import multimerge_genomics as mmg
from multimerge_genomics import exercises as ex

failures, done = [], 0


def check(name, cond, detail=""):
    global done
    done += 1
    if not cond:
        failures.append(name)
    print(("  ok    " if cond else "  FAIL  ") + name + ("" if cond else f"  {detail}"))


T = tempfile.mkdtemp(prefix="mmg_exercises_")
os.chdir(T)
students = ["Ana Souza", "Bruno Lima", "João Pereira"]
with open("students.txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(students) + "\n")

print("== making the class")
r = ex.make_class("students.txt", "class", course="Genomics 2026", zip_each=True)
check("one folder and one zip per student, and the teacher's key",
      all(os.path.isdir(f"class/{s}") and os.path.exists(f"class/{s}.zip") for s in ("ana_souza", "bruno_lima", "joao_pereira"))
      and os.path.exists("class/_teacher/answer_key.csv"), os.listdir("class"))
check("5 exercises x 3 students; 13 questions each", r["questions"] == 39 and len(r["exercises"]) == 5, r)
key = list(csv.DictReader(open("class/_teacher/answer_key.csv", encoding="utf-8")))
check("the README of a student lists every question", all(q["question"] in open("class/ana_souza/README.md", encoding="utf-8").read()
                                                         for q in key if q["folder"] == "ana_souza"))
ex.make_class("students.txt", "again", course="Genomics 2026")
same = all(filecmp.cmp(f"class/ana_souza/{p}", f"again/ana_souza/{p}", shallow=False)
           for p in ("snps/patient.fasta", "broken_vcf/lab_results.vcf", "adapter/reads.fastq.gz", "paternity/family.vcf", "expression/counts.tsv"))
check("the same course and student give the same bytes", same)
check("two students get different data",
      not filecmp.cmp("class/ana_souza/snps/patient.fasta", "class/bruno_lima/snps/patient.fasta", shallow=False)
      and not filecmp.cmp("class/ana_souza/paternity/family.vcf", "class/bruno_lima/paternity/family.vcf", shallow=False))
ex.make_class(["Ana Souza"], "other_course", course="Genomics 2027")
check("another course gives the same student other data",
      not filecmp.cmp("class/ana_souza/snps/patient.fasta", "other_course/ana_souza/snps/patient.fasta", shallow=False))
try:
    ex.make_class(["x"], "bad", exercises=["snps", "nope"])
    check("an unknown exercise is refused", False)
except ValueError as err:
    check("an unknown exercise is refused, naming the exercises", "nope" in str(err) and "paternity" in str(err), err)


def expected(folder, q):
    return next(r["answer"] for r in key if r["folder"] == folder and r["question"] == q)


print("== solving every exercise as a student would, and grading")
for folder in ("ana_souza", "bruno_lima", "joao_pereira"):
    d = f"class/{folder}"
    ans = {}
    # screens 1 and 4: k-mers, a Bloom filter of the reference, the patient's k-mers that are not in it
    ref = mmg.read_fasta_batch([f"{d}/snps/reference.fasta"], 31)
    pat = mmg.read_fasta_batch([f"{d}/snps/patient.fasta"], 31)
    bloom = mmg.build_bloom_filter(ref, 64, 5)
    novel = len(mmg.extract_mutations(pat, bloom, 5)) // 16
    ans["snps_kmers"], ans["snps_count"] = novel, novel // 31
    # screen 10
    rep = mmg.validate_file(f"{d}/broken_vcf/lab_results.vcf")
    ans["vcf_kinds"] = len(rep["problems"])
    ans["vcf_first_line"] = min(p["line"] for p in rep["problems"] if p.get("line"))
    ans["vcf_errors"] = sum(p["severity"] == "error" for p in rep["problems"])
    # screen 12: the report names the adapter; trimming with it
    st = mmg.fastq_stats(f"{d}/adapter/reads.fastq.gz")
    likely = st.get("likely_adapter") or {}
    ans["adapter_name"], seq = likely.get("name", ""), likely.get("sequence")
    t = mmg.fastq_trim(f"{d}/adapter/reads.fastq.gz", f"{d}/adapter/trimmed.fastq", adapter=seq or "auto",
                       quality_cutoff=20, min_length=25)
    ans["adapter_reads_left"] = t["reads_out"]
    # screen 16
    rel = mmg.relatedness([f"{d}/paternity/family.vcf"])
    readme = open(f"{d}/README.md", encoding="utf-8").read()
    child = readme.split(" the child;")[0].split(" and ")[-1].strip()
    mine = [p for p in rel["pairs"] if child in (p["a"], p["b"])]
    other = lambda p: p["b"] if p["a"] == child else p["a"]
    father = max(mine, key=lambda p: p["kinship"] if "mother" not in p["relationship"] else -1)
    cands = [p for p in mine if "is the mother" not in readme or other(p) not in readme.split(" is the mother")[0]]
    by_kin = sorted(cands, key=lambda p: -p["kinship"])
    # the mother and the father are both parent/child with the child; the README names the mother
    mother = readme.split("(no real person): ")[1].split(" is the mother")[0]
    rest = [p for p in by_kin if other(p) != mother]
    ans["father"], ans["father_kinship"], ans["uncle"] = other(rest[0]), round(rest[0]["kinship"], 2), other(rest[1])
    # screen 22
    de = mmg.differential_expression(f"{d}/expression/counts.tsv", f"{d}/expression/samples.tsv", "condition", "control",
                                     output=f"{d}/expression/de.tsv", alpha=0.1)
    tab = list(csv.DictReader(open(f"{d}/expression/de.tsv", encoding="utf-8"), delimiter="\t"))
    sig = [r for r in tab if r["padj"] not in ("", "NA", "nan") and float(r["padj"]) < 0.1]
    ans["de_up"], ans["de_down"] = de["up"], de["down"]
    ans["de_top"] = max(sig, key=lambda r: abs(float(r["log2FoldChange"])))["gene"]
    check(f"{folder}: the Bloom-filter count of screen 4 is the key's (each SNP: 31 k-mers)",
          str(ans["snps_kmers"]) == expected(folder, "snps_kmers") and ans["snps_kmers"] == 31 * ans["snps_count"],
          (ans["snps_kmers"], expected(folder, "snps_kmers")))
    check(f"{folder}: the father is the one with ~0.25, the uncle ~0.125",
          0.2 < ans["father_kinship"] < 0.3 and 0.08 < rest[1]["kinship"] < 0.17, (ans["father_kinship"], rest[1]["kinship"]))
    check(f"{folder}: the report's adapter is accepted for the key's", ex._correct(
        next(r for r in key if r["folder"] == folder and r["question"] == "adapter_name"), ans["adapter_name"]),
        (ans["adapter_name"], expected(folder, "adapter_name")))
    with open(f"{d}/answers.csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["question", "answer"])
        for q, a in ans.items():
            w.writerow([q, a])

# one student gets two answers wrong
rows = list(csv.reader(open("class/bruno_lima/answers.csv", encoding="utf-8")))
for r in rows:
    if r[0] in ("snps_count", "father"):
        r[1] = "999" if r[0] == "snps_count" else "nobody"
with open("class/bruno_lima/answers.csv", "w", encoding="utf-8", newline="") as fh:
    csv.writer(fh).writerows(rows)
g = ex.grade("class/_teacher/answer_key.csv", "class", out="grades.csv")
sc = g["students"]
check("solved on the screens' tools: 13/13", sc["Ana Souza"]["correct"] == 13 and sc["João Pereira"]["correct"] == 13, sc)
check("two wrong answers: 11/13", sc["Bruno Lima"]["correct"] == 11, sc)
check("grades.csv has one row per answer", len(list(csv.DictReader(open("grades.csv", encoding="utf-8")))) == 39)
# one CSV for the whole class, and the command line
with open("all.csv", "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["student", "question", "answer"])
    for r in key:
        if r["student"] == "Ana Souza":
            w.writerow([r["student"], r["question"], r["answer"]])
g2 = ex.grade("class/_teacher/answer_key.csv", "all.csv")
check("one CSV for the class: the key's own answers score 100%, missing students 0", g2["students"]["Ana Souza"]["percent"] == 100.0
      and g2["students"]["Bruno Lima"]["correct"] == 0, g2["students"])
check("numbers with a comma, extra spaces and case are accepted",
      ex._correct({"kind": "float", "answer": "0.25", "tolerance": "0.02"}, " 0,26 ")
      and ex._correct({"kind": "text", "answer": "Nextera", "accepted": "nextera|ctgtctcttatacacatct"}, "  NEXTERA ")
      and not ex._correct({"kind": "int", "answer": "93", "tolerance": "0"}, "94"))
rc = ex.main(["grade", "--key", "class/_teacher/answer_key.csv", "--answers", "class"])
check("the command line grades too", rc == 0)

print()
print(f"{done - len(failures)} of {done} checks passed")
if failures:
    print("FAILED:", failures)
    sys.exit(1)
