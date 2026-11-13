# End-to-end test of the Python API with real data.
# Usage: python e2e.py <data dir> <work dir>
# Runs the same on Linux and Windows. In the data directory it looks for the two
# dbSNP b151 common files (GRCh38 and GRCh37) and the RefSeq GFF, by the short
# names used on Colab or by their original NCBI names.
import gzip
import os
import random
import sys
import threading
import time

import multimerge_genomics as mmg

DADOS, T = sys.argv[1], sys.argv[2]
os.makedirs(T, exist_ok=True)


def acha(*nomes):
    for nome in nomes:
        if os.path.exists(f"{DADOS}/{nome}"):
            return f"{DADOS}/{nome}"
    sys.exit(f"Missing in {DADOS}: one of {nomes}")


V38 = acha("dbsnp38.vcf.gz", "dbsnp151_common_all.vcf.gz")
V37 = acha("dbsnp37.vcf.gz", "dbsnp151_common_all_GRCh37.vcf.gz")
GFF = acha("refseq.gff.gz", "refseq_GRCh38p14_genomic.gff.gz")
falhas = []
feitos = 0


def confere(nome, cond, detalhe=""):
    global feitos
    feitos += 1
    if cond:
        print(f"  ok   {nome}")
    else:
        falhas.append(nome)
        print(f"  FAIL {nome} {detalhe}")


def relogio(f, *a, **k):
    t = time.perf_counter()
    r = f(*a, **k)
    return r, time.perf_counter() - t


def abre(caminho, modo):
    # Sempre LF e UTF-8: no Windows o modo texto trocaria os fins de linha.
    ab = gzip.open if caminho.endswith(".gz") else open
    return ab(caminho, modo + "t", newline="\n", encoding="utf-8")


def linhas(caminho):
    with abre(caminho, "r") as f:
        return f.read().split("\n")


def registros(caminho):
    return [l for l in linhas(caminho) if l and not l.startswith("#")]


def chave_natural(c):
    return (0, int(c), "") if c.isdigit() else (1, 0, c)


print("== 1. validate_file (real VCF, 2 million lines)")
r, dt = relogio(mmg.validate_file, V38, max_lines=2_000_000)
print(r["text"])
print(f"  [{dt:.2f} s]")
confere("returns a dictionary with the promised keys",
        all(k in r for k in ("valid", "errors", "warnings", "lines", "records", "complete", "sorted", "naming", "problems", "text")),
        str(sorted(r)))
confere("dbSNP is valid", r["valid"] is True, str(r["errors"]))
confere("naming without chr", r["naming"] == "no_chr", str(r["naming"]))
confere("line limit respected", r["complete"] is False and r["lines"] <= 2_000_000, str(r["lines"]))

print("== 2. file_reference / compare_references")
r38 = mmg.file_reference(V38)
r37 = mmg.file_reference(V37)
print("  38:", {k: v for k, v in r38.items() if k != "sequences"})
print("  37:", {k: v for k, v in r37.items() if k != "sequences"})
c = mmg.compare_references(V38, V37)
print("  verdict:", c["verdict"])
confere("38 x 37 incompatible", c["compatible"] is False, str(c["compatible"]))
c2 = mmg.compare_references(V38, V38)
confere("38 x 38 compatible", c2["compatible"] is True, str(c2))

print("== 3. sample: 1 in every 40 records of the whole VCF")
AM = f"{T}/amostra.vcf"


def amostra(origem, destino, passo):
    # O cabecalho inteiro e uma linha a cada `passo`, em bytes crus.
    with gzip.open(origem, "rb") as f, open(destino, "wb") as g:
        for nr, linha in enumerate(f, 1):
            if nr % passo == 0 or linha.startswith(b"#"):
                g.write(linha)


_, dt = relogio(amostra, V38, AM, 40)
orig = linhas(AM)
cab = [l for l in orig if l.startswith("#")]
regs = [l for l in orig if l and not l.startswith("#")]
croms = []
for l in regs:
    c0 = l.split("\t", 1)[0]
    if not croms or croms[-1] != c0:
        croms.append(c0)
print(f"  {len(regs)} records, {len(croms)} chromosomes, {dt:.1f} s to extract")

print("== 4. normalize_chromosomes: round trip")
CHR, VOLTA = f"{T}/amostra_chr.vcf", f"{T}/amostra_volta.vcf"
n, dt = relogio(mmg.normalize_chromosomes, AM, CHR, style="chr")
print("  ", {k: (v if not isinstance(v, list) else v[:6]) for k, v in n.items()}, f"[{dt:.2f} s]")
confere("every record changed", n["changed"] == len(regs), f"{n['changed']} != {len(regs)}")
v = mmg.validate_file(CHR)
confere("the chr result is valid", v["naming"] == "chr" and v["valid"], str(v["naming"]))
confere("MT became chrM", any(l.startswith("chrM\t") for l in registros(CHR)) == ("MT" in croms))
mmg.normalize_chromosomes(CHR, VOLTA, style="no_chr")
confere("the round trip returns the original file", linhas(VOLTA) == orig)

print("== 5. sort_genomic: shuffled -> sorted, with runs on disk")
random.seed(7)
emb = regs[:]
random.shuffle(emb)
EMB, ORD = f"{T}/embaralhado.vcf", f"{T}/ordenado.vcf"
with abre(EMB, "w") as f:
    f.write("\n".join(cab + emb) + "\n")
ref = sorted(emb, key=lambda l: (chave_natural(l.split("\t", 2)[0]), int(l.split("\t", 2)[1])))
(nreg, lotes, ign), dt = relogio(mmg.sort_genomic, [EMB], ORD, chrom_order="natural", memory_mb=16)
print(f"  records={nreg} runs={lotes} skipped={ign} [{dt:.2f} s]")
confere("used runs on disk", lotes > 1, str(lotes))
confere("record count", nreg == len(regs))
confere("identical to Python's stable sort", registros(ORD) == ref)
confere("header preserved", [l for l in linhas(ORD) if l.startswith("#")] == cab)
v = mmg.validate_file(ORD)
confere("the validator confirms the order", v["sorted"] is True, str(v["sorted"]))
v = mmg.validate_file(EMB)
confere("the validator flags the shuffled file", v["sorted"] is False and any(p["code"] == "unsorted" for p in v["problems"]))

print("== 5b. compressed output and two inputs")
ORDGZ = f"{T}/ordenado2.vcf.gz"
meio = len(emb) // 2
P1, P2 = f"{T}/p1.vcf", f"{T}/p2.vcf.gz"
with abre(P1, "w") as f:
    f.write("\n".join(cab + emb[:meio]) + "\n")
with abre(P2, "w") as f:
    f.write("\n".join(cab + emb[meio:]) + "\n")
(nreg, lotes, ign), dt = relogio(mmg.sort_genomic, [P1, P2], ORDGZ, chrom_order="natural")
confere("two inputs (text + gz) -> identical .gz", registros(ORDGZ) == ref, f"{nreg} {lotes}")
n_v, dt = relogio(mmg.multimerge_vcf, [P1, P2], "POS", f"{T}/legado.vcf", 16)
confere("multimerge_vcf returns the count", n_v == len(regs), str(n_v))

print("== 6. BED: merge and intersect against a Python reference")
random.seed(11)


def gera_bed(n, croms_bed, maxpos, maxlen):
    out = []
    for i in range(n):
        c0 = random.choice(croms_bed)
        a = random.randrange(maxpos)
        out.append((c0, a, a + 1 + random.randrange(maxlen), f"r{i}"))
    return out


def grava_bed(caminho, ivs):
    with abre(caminho, "w") as f:
        for iv in ivs:
            f.write("\t".join(map(str, iv)) + "\n")


def le_bed(caminho):
    return [tuple(l.split("\t")) for l in linhas(caminho) if l]


CB = ["chr1", "chr2", "chr10", "chrX"]
A = gera_bed(40_000, CB, 2_000_000, 300)
B = gera_bed(30_000, CB[:3], 2_000_000, 500)
BA, BB = f"{T}/a.bed", f"{T}/b.bed"
grava_bed(BA, A)
grava_bed(BB, B)


def ref_fundir(ivs, d):
    out = []
    for c0 in sorted(set(i[0] for i in ivs)):
        atual = None
        for _, a, b, _ in sorted(i for i in ivs if i[0] == c0):
            if atual and a <= atual[1] + d:
                atual[1] = max(atual[1], b)
            else:
                if atual:
                    out.append((c0, str(atual[0]), str(atual[1])))
                atual = [a, b]
        if atual:
            out.append((c0, str(atual[0]), str(atual[1])))
    return out


for d in (0, 100):
    (lidos, grav), dt = relogio(mmg.bed_merge, [BA], f"{T}/fund{d}.bed", distance=d, memory_mb=1)
    confere(f"bed_merge d={d} ({lidos} -> {grav})", le_bed(f"{T}/fund{d}.bed") == ref_fundir(A, d) and lidos == len(A))

bpc = {}
for c0, a, b, _ in B:
    bpc.setdefault(c0, []).append((a, b))


def toca(iv):
    return [(a, b) for a, b in bpc.get(iv[0], []) if a < iv[2] and iv[1] < b]


# So uma fatia de A no brute force: 4000 x 30000.
A2 = A[:4000]
BA2 = f"{T}/a2.bed"
grava_bed(BA2, A2)
n_u = mmg.bed_intersect([BA2], [BB], f"{T}/u.bed", mode="unique")
n_v2 = mmg.bed_intersect([BA2], [BB], f"{T}/v.bed", mode="absent")
n_r = mmg.bed_intersect([BA2], [BB], f"{T}/r.bed", mode="regions")
ref_u = sorted((i[0], str(i[1]), str(i[2]), i[3]) for i in A2 if toca(i))
ref_v = sorted((i[0], str(i[1]), str(i[2]), i[3]) for i in A2 if not toca(i))
ref_r = sorted((i[0], str(max(i[1], a)), str(min(i[2], b)), i[3]) for i in A2 for a, b in toca(i))
confere(f"bed_intersect unique ({n_u})", sorted(le_bed(f"{T}/u.bed")) == ref_u and n_u == len(ref_u))
confere(f"bed_intersect absent ({n_v2})", sorted(le_bed(f"{T}/v.bed")) == ref_v and n_v2 == len(ref_v))
confere(f"bed_intersect regions ({n_r})", sorted(le_bed(f"{T}/r.bed")) == ref_r and n_r == len(ref_r))

print("== 7. errors say what happened")
for nome, f in [
    ("missing file", lambda: mmg.sort_genomic([f"{T}/does_not_exist.vcf"], f"{T}/x.vcf")),
    ("unknown format", lambda: mmg.sort_genomic([EMB], f"{T}/x.vcf", format="xls")),
    ("unknown sort key", lambda: mmg.sort_genomic([EMB], f"{T}/x.vcf", by="colour")),
    ("unknown style", lambda: mmg.normalize_chromosomes(AM, f"{T}/x.vcf", style="purple")),
    ("validating a missing file", lambda: mmg.validate_file(f"{T}/does_not_exist.vcf")),
    ("intersect with no name in common", lambda: mmg.bed_intersect([BA], [AM.replace('.vcf', '_semchr.bed')], f"{T}/x.bed")),
]:
    if nome.startswith("intersect"):
        grava_bed(AM.replace(".vcf", "_semchr.bed"), [("1", 10, 20, "x"), ("2", 5, 9, "y")])
    try:
        f()
        confere(nome, False, "no error was raised")
    except RuntimeError as e:
        print(f"       -> {str(e)[:150]}")
        confere(nome, len(str(e)) > 10)

print("== 8. the GIL is released while sorting")
tiques = [0]
fim = threading.Event()


def conta():
    while not fim.is_set():
        tiques[0] += 1
        time.sleep(0.005)


th = threading.Thread(target=conta)
th.start()
_, dt = relogio(mmg.sort_genomic, [EMB], f"{T}/gil.vcf", chrom_order="natural", memory_mb=16)
fim.set()
th.join()
print(f"  {tiques[0]} ticks in {dt:.2f} s")
confere("the other thread kept running while sorting", tiques[0] >= dt / 0.005 * 0.3, str(tiques[0]))

print("== 9. real GFF: validate and sort")
v, dt = relogio(mmg.validate_file, GFF)
print(f"  {v['records']} records, valid={v['valid']}, errors={v['errors']}, warnings={v['warnings']} [{dt:.1f} s]")
(nreg, lotes, ign), dt = relogio(mmg.sort_genomic, [GFF], f"{T}/refseq_ord.gff", memory_mb=256)
print(f"  sorted: records={nreg} runs={lotes} skipped={ign} [{dt:.1f} s]")
confere("GFF: every record", nreg == v["records"], f"{nreg} != {v['records']}")
v2 = mmg.validate_file(f"{T}/refseq_ord.gff")
confere("the sorted GFF passes the validator as sorted", v2["sorted"] is True and v2["records"] == nreg, str(v2["sorted"]))

print()
print(f"{feitos - len(falhas)} of {feitos} checks passed")
if falhas:
    print("FAILED:", falhas)
    sys.exit(1)
