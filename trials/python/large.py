# Real-size sorting through Python: the two dbSNP files (18.7 GB of text,
# 75 million variants), with a small budget and with the automatic one.
# Usage: python large.py <data dir> <work dir>
import os
import sys
import time

import multimerge_genomics as mmg

DADOS, T = sys.argv[1], sys.argv[2]
os.makedirs(T, exist_ok=True)
entradas = [f"{DADOS}/dbsnp38.vcf.gz", f"{DADOS}/dbsnp37.vcf.gz"]
print(mmg.engine_info(), flush=True)

for rotulo, memoria in (("1 GB", 1024), ("automatic", None)):
    saida = f"{T}/tudo_{memoria}.vcf"
    t = time.perf_counter()
    reg, lotes, ign = mmg.sort_genomic(entradas, saida, chrom_order="natural", memory_mb=memoria)
    dt = time.perf_counter() - t
    gb = os.path.getsize(saida) / 2**30
    print(f"budget {rotulo}: {reg:,} records, {lotes} runs, {gb:.2f} GB written in {dt:.1f} s ({gb * 1024 / dt:.0f} MB/s)", flush=True)
    t = time.perf_counter()
    v = mmg.validate_file(saida)
    print(f"  validator: {v['records']:,} records, sorted={v['sorted']}, errors={v['errors']}, warnings={v['warnings']} [{time.perf_counter() - t:.1f} s]", flush=True)
    for p in v["problems"][:5]:
        print("   ", p["severity"], p["code"], "line", p["line"], "x", p["occurrences"], "-", p["message"][:140], flush=True)
    os.remove(saida)
