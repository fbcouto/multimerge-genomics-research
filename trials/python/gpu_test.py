# Test of the GPU modules: results checked against a numpy reference, GPU
# against CPU bit for bit, and timings (best of 3, after a warm-up run).
# Usage: python gpu_test.py [n_millions ...]
import sys
import time

import numpy as np

import multimerge_genomics as mmg

falhas = []
feitos = 0


def confere(nome, cond, detalhe=""):
    global feitos
    feitos += 1
    print(("  ok   " if cond else "  FAIL ") + nome + ("" if cond else " " + str(detalhe)))
    if not cond:
        falhas.append(nome)


def melhor(f, vezes=5):
    f()  # aquecimento descartado
    ts = []
    for _ in range(vezes):
        t = time.perf_counter()
        r = f()
        ts.append(time.perf_counter() - t)
    return r, min(ts) * 1000, sorted(ts)[len(ts) // 2] * 1000


M64 = np.uint64(0xFFFFFFFFFFFFFFFF)


def hash_np(k, semente):
    k = k ^ np.uint64(semente)
    k ^= k >> np.uint64(33)
    k *= np.uint64(0xFF51AFD7ED558CCD)
    k ^= k >> np.uint64(33)
    k *= np.uint64(0xC4CEB9FE1A85EC53)
    k ^= k >> np.uint64(33)
    return k


def bloom_np(kmers64, nbytes, nh):
    bits = np.uint64((nbytes // 4) * 32)
    pal = np.zeros(nbytes // 4, dtype=np.uint32)
    for i in range(nh):
        b = hash_np(kmers64.copy(), (i * 0x9E3779B9) & 0xFFFFFFFF) % bits
        np.bitwise_or.at(pal, (b // np.uint64(32)).astype(np.int64), (np.uint32(1) << (b % np.uint64(32)).astype(np.uint32)))
    return pal.tobytes()


def consulta_np(kmers64, filtro, nh):
    pal = np.frombuffer(filtro, dtype=np.uint32)
    bits = np.uint64(len(pal) * 32)
    ok = np.ones(len(kmers64), dtype=bool)
    for i in range(nh):
        b = hash_np(kmers64.copy(), (i * 0x9E3779B9) & 0xFFFFFFFF) % bits
        ok &= (pal[(b // np.uint64(32)).astype(np.int64)] >> (b % np.uint64(32)).astype(np.uint32)) & np.uint32(1) == 1
    return ok


def dobra(k128):
    v = k128.view(np.uint64).reshape(-1, 2)  # little-endian: [baixo, alto]
    return v[:, 0] ^ v[:, 1]


info = mmg.engine_info()
print("engine_info:", info)
MODOS = ["cpu"] + (["gpu"] if info["gpu"] else [])

rng = np.random.default_rng(20261113)
tamanhos = [int(float(a) * 1e6) for a in sys.argv[1:]] or [1_000_000, 10_000_000]

np.seterr(over="ignore")
for n in tamanhos:
    print(f"== n = {n:,}")
    # ---- fusao de dois blocos u64
    a = rng.integers(0, 2**63, n, dtype=np.uint64)
    b = rng.integers(0, 2**63, n, dtype=np.uint64)
    a[: n // 10] = b[: n // 10]  # chaves repetidas entre os dois lados
    ba, bb = a.tobytes(), b.tobytes()
    t = time.perf_counter()
    ref = np.sort(np.concatenate([a, b]), kind="stable")
    t_np = (time.perf_counter() - t) * 1000
    saidas = {}
    for modo in MODOS:
        kw = {"device": modo}
        r, tmin, tmed = melhor(lambda: mmg.sort_merge_bytes(ba, bb, **kw), 3)
        saidas[modo] = r
        confere(f"sort_merge_bytes [{modo}] == numpy   min {tmin:8.1f} ms  median {tmed:8.1f} ms  (numpy {t_np:.0f} ms)",
                np.array_equal(np.frombuffer(r, dtype=np.uint64), ref))

    # ---- filtro de Bloom
    k128 = rng.integers(0, 2**63, 2 * n, dtype=np.uint64)  # n k-mers de 128 bits
    ref_k = k128.tobytes()
    amostra = k128.copy()
    mut = rng.choice(n, size=max(1, n // 1000), replace=False)
    amostra.reshape(-1, 2)[mut, 0] ^= np.uint64(0x5DEECE66D)  # 0,1% de k-mers mutados
    am_k = amostra.tobytes()
    k64 = dobra(k128)
    am64 = dobra(amostra)
    MB, NH = 64, 5
    filtros = {}
    for modo in MODOS:
        kw = {"device": modo}
        f, tmin, tmed = melhor(lambda: mmg.build_bloom_filter(ref_k, MB, NH, **kw), 3)
        filtros[modo] = f
        print(f"  bloom build [{modo}]  min {tmin:8.1f} ms  median {tmed:8.1f} ms")
    if n <= 2_000_000:
        confere("filter == numpy reference (bit for bit)", all(f == bloom_np(k64, MB << 20, NH) for f in filtros.values()))
    if len(filtros) > 1:
        confere("GPU filter == CPU filter (bit for bit)", filtros["gpu"] == filtros["cpu"])
    filtro = next(iter(filtros.values()))
    confere("the filter is not all zero", any(filtro[:4096]) or np.frombuffer(filtro, dtype=np.uint8).any())
    esperado = am64[~consulta_np(am64, filtro, NH)] if n <= 2_000_000 else None
    muts = {}
    for modo in MODOS:
        kw = {"device": modo}
        m, tmin, tmed = melhor(lambda: mmg.extract_mutations(am_k, filtro, NH, **kw), 3)
        muts[modo] = m
        nm = len(m) // 16
        print(f"  bloom query [{modo}]  min {tmin:8.1f} ms  median {tmed:8.1f} ms   {nm} mutants of {len(mut)} planted")
        confere(f"mutants [{modo}]: found nearly all the planted ones and nothing else", 0.95 * len(mut) <= nm <= len(mut), nm)
        if esperado is not None:
            confere(f"mutants [{modo}] == numpy reference", np.array_equal(dobra(np.frombuffer(m, dtype=np.uint64)), esperado))
    if len(muts) > 1:
        confere("GPU mutants == CPU mutants", muts["gpu"] == muts["cpu"])

print()
print(f"{feitos - len(falhas)} of {feitos} checks passed")
if falhas:
    print("FAILED:", falhas)
    sys.exit(1)
