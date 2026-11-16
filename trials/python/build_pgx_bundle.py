#!/usr/bin/env python3
"""Build versioned, CC0-derived JSON lookup tables for the mmgen-pgx engine.

This maintainer script parses an official CPIC PostgreSQL dump (CC0 1.0 Public
Domain) and GRCh38 reference sequence windows to build a deterministic,
checksummed JSON table bundle:

  mmgen_tools/pgx/<version>/
    manifest.json
    genes.json
    phenotypes.json
    guidelines.json

And packages it as a reproducible zip archive:
  mmgen_pgx_tables_<version>.zip

Usage:
  python trials/python/build_pgx_bundle.py \
      --dump path/to/cpic_db_dump-v1.61.0.sql.gz \
      --out mmgen_tools/pgx/cpic-v1.61.0

Optional flags:
  --download              Download official CPIC dump from files.cpicpgx.org
  --cpic-version v1.61.0  CPIC dump version (default: v1.61.0)
  --fasta GRCh38.fa       Optional local indexed FASTA (uses .fai)
  --ref-cache cache.json  Cache file for GRCh38 sequence windows
  --zip out.zip           Output zip path (default: <parent>/mmgen_pgx_tables_<ver>.zip)
  --verify-vcf pos.vcf    Optional reference VCF to cross-check variant coordinates
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import sys
import tempfile
import urllib.request
import zipfile
from collections import defaultdict
from datetime import date

TARGET_GENES = [
    "ABCG2",
    "CACNA1S",
    "CFTR",
    "CYP2B6",
    "CYP2C19",
    "CYP2C9",
    "CYP3A5",
    "CYP4F2",
    "DPYD",
    "G6PD",
    "IFNL3",
    "MT-RNR1",
    "NAT2",
    "NUDT15",
    "RYR1",
    "SLCO1B1",
    "TPMT",
    "UGT1A1",
    "VKORC1",
]

BAM_VCF_GENES = ["CYP2D6"]
BAM_HLA_GENES = ["HLA-A", "HLA-B", "HLA-C"]
ALL_PGX_GENES = TARGET_GENES + BAM_VCF_GENES + BAM_HLA_GENES

GENE_MODELS = {
    "DPYD": "variant_activity",
    "CFTR": "variant_status",
    "RYR1": "variant_status",
    "CACNA1S": "variant_status",
    "MT-RNR1": "variant_status",
    "ABCG2": "variant_genotype",
    "VKORC1": "variant_genotype",
    "IFNL3": "variant_genotype",
}

# Genes whose diplotype is chosen by population frequency rather than by the Named Allele
# Matcher score. PharmCAT does this for NAT2 only (its "unphased-priority" exemptions): the
# NAT2 reference allele, *1, is the GRCh38 haplotype, rare in every population, and a score
# that counts every position read for it prefers it over the common *4, *5, *6 alleles.
# The engine computes the choice from the CPIC frequencies below (pop_freqs), not from
# PharmCAT's list.
FREQUENCY_CHOICE_GENES = {"NAT2"}

IUPAC_MAP = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"],
}

TABLES_TO_EXTRACT = {
    "gene",
    "allele",
    "allele_definition",
    "sequence_location",
    "allele_location_value",
    "gene_result",
    "gene_result_lookup",
    "gene_result_diplotype",
    "guideline",
    "drug",
    "recommendation",
    "pair",
}


# ---------------------------------------------------------------------------
# Streaming CPIC PostgreSQL dump parser
# ---------------------------------------------------------------------------

def iter_sql_statements(path: str):
    """Yield complete SQL statements from a .sql or .sql.gz file."""
    opener = gzip.open if path.endswith(".gz") else open
    buf = []
    in_str = False
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not buf:
                stripped = line.lstrip()
                if not stripped or stripped.startswith("--") or stripped.startswith("SET ") or stripped.startswith("SELECT "):
                    continue
            buf.append(line)
            i = 0
            n = len(line)
            while i < n:
                ch = line[i]
                if in_str:
                    if ch == "\\":
                        i += 2
                        continue
                    if ch == "'":
                        if i + 1 < n and line[i + 1] == "'":
                            i += 2
                            continue
                        in_str = False
                else:
                    if ch == "'":
                        in_str = True
                    elif ch == ";":
                        yield "".join(buf)
                        buf = []
                        break
                i += 1
    if buf:
        yield "".join(buf)


def parse_sql_value(tok: str):
    tok = tok.strip()
    if tok.upper() == "NULL":
        return None
    if tok.lower() == "true":
        return True
    if tok.lower() == "false":
        return False
    if tok.startswith("'") and tok.endswith("'"):
        inner = tok[1:-1]
        return inner.replace("''", "'")
    if re.fullmatch(r"-?\d+", tok):
        return int(tok)
    if re.fullmatch(r"-?\d+\.\d+", tok):
        return float(tok)
    return tok


def parse_insert(stmt: str):
    """Parse an INSERT INTO cpic.<table> (...) VALUES (...), (...); or COPY."""
    m = re.match(
        r"^\s*INSERT\s+INTO\s+(?:cpic\.)?([a-z_]+)\s*\(([^)]+)\)\s*VALUES\s*(.+);\s*$",
        stmt,
        re.DOTALL | re.IGNORECASE,
    )
    if not m:
        return None, None, []
    tbl = m.group(1).lower()
    if tbl not in TABLES_TO_EXTRACT:
        return None, None, []
    cols = [c.strip().strip('"').lower() for c in m.group(2).split(",")]
    body = m.group(3)
    rows = []
    i = 0
    n = len(body)
    while i < n:
        if body[i] != "(":
            i += 1
            continue
        i += 1
        vals = []
        tok_start = i
        in_str = False
        while i < n:
            ch = body[i]
            if in_str:
                if ch == "\\":
                    i += 2
                    continue
                if ch == "'":
                    if i + 1 < n and body[i + 1] == "'":
                        i += 2
                        continue
                    in_str = False
            else:
                if ch == "'":
                    in_str = True
                elif ch == ",":
                    vals.append(parse_sql_value(body[tok_start:i]))
                    tok_start = i + 1
                elif ch == ")":
                    vals.append(parse_sql_value(body[tok_start:i]))
                    i += 1
                    break
            i += 1
        if len(vals) == len(cols):
            rows.append(dict(zip(cols, vals)))
    return tbl, cols, rows


def load_cpic_dump(path: str) -> dict[str, list[dict]]:
    """Load required CPIC tables from pg_dump SQL (COPY or INSERT format)."""
    tables: dict[str, list[dict]] = defaultdict(list)
    opener = gzip.open if path.endswith(".gz") else open

    # Pass 1: handle PostgreSQL COPY blocks in a single fast stream
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        cur_tbl = None
        cur_cols = None
        for line in fh:
            if cur_tbl is None:
                if line.startswith("COPY "):
                    m = re.match(r'^COPY\s+(?:"?[a-z_]+"?\.)?"?([a-z_]+)"?\s*\(([^)]+)\)\s+FROM\s+stdin;', line, re.I)
                    if m and m.group(1).lower() in TABLES_TO_EXTRACT:
                        cur_tbl = m.group(1).lower()
                        cur_cols = [c.strip().strip('"').lower() for c in m.group(2).split(",")]
            else:
                if line.rstrip("\r\n") == r"\.":
                    cur_tbl = None
                    cur_cols = None
                    continue
                raw_fields = line.rstrip("\r\n").split("\t")
                row = {}
                for col, rf in zip(cur_cols, raw_fields):
                    if rf == r"\N":
                        row[col] = None
                    elif rf == "t":
                        row[col] = True
                    elif rf == "f":
                        row[col] = False
                    elif re.fullmatch(r"-?\d+", rf):
                        row[col] = int(rf)
                    else:
                        rf = (
                            rf.replace(r"\\", "\\")
                            .replace(r"\t", "\t")
                            .replace(r"\n", "\n")
                            .replace(r"\r", "\r")
                        )
                        row[col] = rf
                tables[cur_tbl].append(row)
    if tables:
        return tables

    # Fallback: SQL INSERT format (e.g. --inserts dumps)
    for stmt in iter_sql_statements(path):
        tbl, _, rows = parse_insert(stmt)
        if tbl:
            tables[tbl].extend(rows)
    return tables


# ---------------------------------------------------------------------------
# GRCh38 Reference Window Fetching & Tan et al. (2015) VCF Normalization
# ---------------------------------------------------------------------------

def canonical_chrom(raw: str) -> str:
    c = raw.strip()
    if not c.startswith("chr"):
        c = "chr" + c
    if c == "chrMT":
        c = "chrM"
    return c


class FastaReader:
    """Minimal .fai-indexed local GRCh38 FASTA reader."""

    def __init__(self, fasta_path: str):
        fai_path = fasta_path + ".fai"
        if not os.path.exists(fai_path):
            raise FileNotFoundError(f"{fai_path}: FASTA index (.fai) not found")
        self.fasta_path = fasta_path
        self.index = {}
        with open(fai_path, encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) >= 5:
                    name = canonical_chrom(parts[0])
                    self.index[name] = (int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]))

    def fetch(self, chrom: str, start_1based: int, end_1based: int) -> str:
        c = canonical_chrom(chrom)
        if c not in self.index:
            raise KeyError(f"{self.fasta_path}: chromosome {c} not in .fai")
        length, offset, line_bases, line_bytes = self.index[c]
        s0 = max(0, start_1based - 1)
        e0 = min(length, end_1based)
        out = []
        with open(self.fasta_path, "rb") as fh:
            pos = s0
            while pos < e0:
                line_idx = pos // line_bases
                in_line = pos % line_bases
                take = min(e0 - pos, line_bases - in_line)
                byte_off = offset + line_idx * line_bytes + in_line
                fh.seek(byte_off)
                out.append(fh.read(take).decode("ascii").upper())
                pos += take
        return "".join(out)


class RefWindowProvider:
    """Fetches and caches GRCh38 flanking reference windows."""

    def __init__(self, fasta: str | None, cache_path: str):
        self.fasta = FastaReader(fasta) if fasta else None
        self.cache_path = cache_path
        self.cache: dict[str, str] = {}
        self.spans: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
        if cache_path and os.path.exists(cache_path):
            with open(cache_path, encoding="utf-8") as fh:
                self.cache = json.load(fh)
            for k, v in self.cache.items():
                if ":" in k and "-" in k:
                    c, rng = k.split(":", 1)
                    s, e = rng.split("-", 1)
                    self.spans[c].append((int(s), int(e), v))

    def save(self):
        if self.cache_path:
            os.makedirs(os.path.dirname(os.path.abspath(self.cache_path)), exist_ok=True)
            with open(self.cache_path, "w", encoding="utf-8") as fh:
                json.dump(self.cache, fh, indent=2, sort_keys=True)

    def prefetch_span(self, chrom: str, start_1based: int, end_1based: int):
        """Ensure [start_1based, end_1based] is cached so every SNV can verify GRCh38 offline."""
        chrom = canonical_chrom(chrom)
        if self.fasta is not None:
            return
        for s0, e0, _ in self.spans.get(chrom, []):
            if s0 <= start_1based and end_1based <= e0:
                return
        key = f"{chrom}:{start_1based}-{end_1based}"
        chunks = []
        cur = start_1based
        while cur <= end_1based:
            ce = min(end_1based, cur + 400000 - 1)
            chunks.append(self._fetch_http(chrom, cur, ce))
            cur = ce + 1
        seq = "".join(chunks)
        self.cache[key] = seq
        self.spans[chrom].append((start_1based, end_1based, seq))

    def _fetch_http(self, chrom: str, start_1based: int, end_1based: int) -> str:
        # Fetch from UCSC hg38 API, falling back to Ensembl REST API
        ucsc_url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start_1based - 1};end={end_1based}"
        try:
            req = urllib.request.Request(ucsc_url, headers={"User-Agent": "mmgen-pgx-builder/1.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                seq = json.loads(resp.read().decode())["dna"].upper()
                if len(seq) == end_1based - start_1based + 1:
                    return seq
        except Exception:
            pass
        ens_chr = chrom.replace("chr", "")
        if ens_chr == "M":
            ens_chr = "MT"
        ens_url = f"https://rest.ensembl.org/sequence/region/human/{ens_chr}:{start_1based}..{end_1based}:1?content-type=text/plain"
        req = urllib.request.Request(ens_url, headers={"User-Agent": "mmgen-pgx-builder/1.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            return resp.read().decode("ascii").strip().upper()

    def get(self, chrom: str, start_1based: int, end_1based: int) -> str:
        chrom = canonical_chrom(chrom)
        if self.fasta is not None:
            return self.fasta.fetch(chrom, start_1based, end_1based)
        for s0, e0, span_seq in self.spans.get(chrom, []):
            if s0 <= start_1based and end_1based <= e0:
                return span_seq[start_1based - s0 : end_1based - s0 + 1]
        key = f"{chrom}:{start_1based}-{end_1based}"
        if key in self.cache:
            return self.cache[key]
        seq = self._fetch_http(chrom, start_1based, end_1based)
        self.cache[key] = seq
        self.spans[chrom].append((start_1based, end_1based, seq))
        return seq


def normalize_variant(seq: str, win_start: int, pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Tan et al. (2015) VCF parsimonious left-alignment."""
    ref = ref.upper()
    alt = alt.upper()
    if ref == alt:
        return pos, ref, alt
    # Step 1: right-trim & left-shift
    while True:
        if ref[-1] == alt[-1]:
            if len(ref) == 1 or len(alt) == 1:
                if pos <= win_start:
                    break
                pos -= 1
                prev = seq[pos - win_start]
                ref = prev + ref
                alt = prev + alt
            ref = ref[:-1]
            alt = alt[:-1]
        else:
            break
    # Step 2: left-trim while keeping at least 1 base
    while len(ref) >= 2 and len(alt) >= 2 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1
    return pos, ref, alt


def expand_cpic_code(ref_code: str, raw_code: str) -> list[str]:
    """Expand a single CPIC allele_location_value string into atomic CPIC codes."""
    if len(ref_code) == 1 and ref_code in IUPAC_MAP and raw_code in IUPAC_MAP:
        return list(IUPAC_MAP[raw_code])
    if " or " in raw_code:
        return [x.strip() for x in raw_code.split(" or ") if x.strip()]
    if ";" in raw_code:
        return [x.strip() for x in raw_code.split(";") if x.strip()]
    return [raw_code]


def cpic_to_vcf_alleles(
    chrom: str,
    cpic_pos: int,
    ref_code: str,
    code: str,
    ref_provider: RefWindowProvider,
) -> list[tuple[int, str, str]]:
    """Convert one atomic CPIC allele_location_value string into normalized GRCh38 (pos, ref, alt) tuples."""
    # Verify SNVs directly against GRCh38 (allowing IUPAC ambiguity in ref_code such as CYP2D6 *1 R/S)
    if len(ref_code) == 1 and ref_code in IUPAC_MAP and code in IUPAC_MAP:
        actual = ref_provider.get(chrom, cpic_pos, cpic_pos)
        if actual not in IUPAC_MAP[ref_code]:
            raise ValueError(f"{chrom}:{cpic_pos}: CPIC ref_code={ref_code} differs from GRCh38={actual}")
        return [(cpic_pos, actual, base) for base in IUPAC_MAP[code]]

    # Fetch generous 120bp window around cpic_pos
    win_start = max(1, cpic_pos - 60)
    win_end = cpic_pos + 80
    seq = ref_provider.get(chrom, win_start, win_end)

    def base_at(p: int) -> str:
        return seq[p - win_start]

    # Case 1: semicolon-separated or ' or ' alternatives
    if ";" in code or " or " in code:
        out = []
        for sub in expand_cpic_code(ref_code, code):
            out.extend(cpic_to_vcf_alleles(chrom, cpic_pos, ref_code, sub, ref_provider))
        return out

    # If code equals reference code
    if code == ref_code:
        if len(ref_code) == 1 and ref_code in "ACGT":
            return [(cpic_pos, ref_code, ref_code)]
        if ref_code == "del":
            b = base_at(cpic_pos)
            return [(cpic_pos, b, b)]
        # Verify repeat or multi-base ref_code against GRCh38
        m_rep = re.fullmatch(r"([ACGT]+)\((\d+)\)", ref_code)
        if m_rep:
            unit, k = m_rep.group(1), int(m_rep.group(2))
            expected = unit * k
            actual = seq[cpic_pos - win_start : cpic_pos - win_start + len(expected)]
            if actual != expected:
                raise ValueError(f"{chrom}:{cpic_pos}: expected repeat {expected}, found {actual} in GRCh38")
            b = base_at(cpic_pos)
            return [(cpic_pos, b, b)]
        if ref_code.startswith("del") and all(c in "ACGT" for c in ref_code[3:]):
            # MNV where ref is denoted delTC -> actual GRCh38 is TC
            expected = ref_code[3:]
            actual = seq[cpic_pos - win_start : cpic_pos - win_start + len(expected)]
            if actual != expected:
                raise ValueError(f"{chrom}:{cpic_pos}: expected MNV ref {expected}, found {actual}")
            return [(cpic_pos, expected, expected)]
        if all(c in "ACGT" for c in ref_code):
            actual = seq[cpic_pos - win_start : cpic_pos - win_start + len(ref_code)]
            if actual != ref_code:
                raise ValueError(f"{chrom}:{cpic_pos}: expected {ref_code}, found {actual}")
            return [(cpic_pos, ref_code[0], ref_code[0])]
        raise ValueError(f"{chrom}:{cpic_pos}: unhandled ref_code {ref_code!r}")

    # Non-reference alleles:
    # A) Deletion: code = del<SEQ> and ref_code = <SEQ>
    if code.startswith("del") and not code.startswith("delins") and all(c in "ACGT" for c in code[3:]):
        deleted = code[3:]
        actual = seq[cpic_pos - win_start : cpic_pos - win_start + len(deleted)]
        if actual != deleted:
            raise ValueError(f"{chrom}:{cpic_pos}: expected deleted seq {deleted}, found {actual} in GRCh38")
        anchor_pos = cpic_pos - 1
        anchor_base = base_at(anchor_pos)
        return [normalize_variant(seq, win_start, anchor_pos, anchor_base + deleted, anchor_base)]

    # B) Insertion after cpic_pos: ref_code == "del" and code = ins<SEQ>
    if ref_code == "del" and code.startswith("ins") and all(c in "ACGT" for c in code[3:]):
        inserted = code[3:]
        anchor_base = base_at(cpic_pos)
        return [normalize_variant(seq, win_start, cpic_pos, anchor_base, anchor_base + inserted)]

    # C) Single- or multi-base duplication: code == ref_code + ref_code (e.g. A->AA or AC->ACAC)
    if len(ref_code) >= 1 and all(c in "ACGT" for c in ref_code) and code == ref_code * 2:
        actual = seq[cpic_pos - win_start : cpic_pos - win_start + len(ref_code)]
        if actual != ref_code:
            raise ValueError(f"{chrom}:{cpic_pos}: dup expected ref {ref_code}, found {actual}")
        anchor_pos = cpic_pos - 1
        anchor_base = base_at(anchor_pos)
        return [normalize_variant(seq, win_start, anchor_pos, anchor_base + ref_code, anchor_base + code)]

    # D) Repeat expansion/contraction: ref_code = UNIT(k) or implicit k=1 UNIT, code = UNIT(m)
    m_alt = re.fullmatch(r"([ACGT]+)\((\d+)\)", code)
    if m_alt:
        m_ref = re.fullmatch(r"([ACGT]+)\((\d+)\)", ref_code)
        if m_ref and m_ref.group(1) == m_alt.group(1):
            unit, k_ref = m_ref.group(1), int(m_ref.group(2))
        elif all(c in "ACGT" for c in ref_code) and m_alt.group(1) == ref_code:
            unit, k_ref = ref_code, 1
        else:
            unit, k_ref = None, 0
        if unit is not None:
            k_alt = int(m_alt.group(2))
            ref_rep = unit * k_ref
            alt_rep = unit * k_alt
            actual = seq[cpic_pos - win_start : cpic_pos - win_start + len(ref_rep)]
            if actual != ref_rep:
                raise ValueError(f"{chrom}:{cpic_pos}: repeat expected {ref_rep}, found {actual}")
            anchor_pos = cpic_pos - 1
            anchor_base = base_at(anchor_pos)
            return [normalize_variant(seq, win_start, anchor_pos, anchor_base + ref_rep, anchor_base + alt_rep)]

    # E) MNV: ref_code = del<SEQ1>, code = ins<SEQ2>
    if ref_code.startswith("del") and code.startswith("ins"):
        old_seq = ref_code[3:]
        new_seq = code[3:]
        actual = seq[cpic_pos - win_start : cpic_pos - win_start + len(old_seq)]
        if actual != old_seq:
            raise ValueError(f"{chrom}:{cpic_pos}: MNV expected {old_seq}, found {actual}")
        return [normalize_variant(seq, win_start, cpic_pos, old_seq, new_seq)]

    # F) delins: single-base ref replaced by multi-base alt (e.g. MT-RNR1 m.961delTinsCC)
    if len(ref_code) == 1 and code.startswith("delins") and all(c in "ACGT" for c in code[6:]):
        actual = base_at(cpic_pos)
        if actual != ref_code:
            raise ValueError(f"{chrom}:{cpic_pos}: delins expected {ref_code}, found {actual}")
        return [normalize_variant(seq, win_start, cpic_pos, ref_code, code[6:])]

    raise ValueError(f"{chrom}:{cpic_pos}: cannot normalize ref={ref_code!r} -> alt={code!r}")


# ---------------------------------------------------------------------------
# Bundle Builders: genes.json, phenotypes.json, guidelines.json, manifest.json
# ---------------------------------------------------------------------------

def max_pop_freq(freq_json_str: str | None) -> float | None:
    if not freq_json_str:
        return None
    try:
        obj = json.loads(freq_json_str)
    except Exception:
        return None
    vals = [float(v) for v in obj.values() if v is not None]
    return max(vals) if vals else None


def pop_freqs(freq_json_str: str | None) -> dict[str, float]:
    """The allele's frequency in each CPIC population group, e.g. {"European": 0.419, ...}."""
    if not freq_json_str:
        return {}
    try:
        obj = json.loads(freq_json_str)
    except Exception:
        return {}
    return {k: round(float(v), 6) for k, v in sorted(obj.items()) if v is not None}


def is_functional_allele(clinical_func: str | None) -> bool:
    if not clinical_func:
        return True  # e.g. VKORC1, IFNL3, CYP4F2 (single/small site or genotype)
    norm = clinical_func.strip().lower()
    return norm not in ("normal function", "iv/normal", "normal risk of aminoglycoside-induced hearing loss", "uncertain function", "unknown function", "uncertain", "none")


def build_traits_entries() -> list[dict]:
    """Open, peer-reviewed GRCh38 trait loci (no proprietary/patented regression matrices)."""
    def snv_pos(pid: int, chrom: str, pos: int, rsid: str, ref: str, alts: list[str], is_key: bool) -> dict:
        ordered = [ref] + [b for b in alts if b != ref] + [b for b in ("A", "C", "G", "T") if b != ref and b not in alts]
        return {
            "id": pid,
            "chrom": chrom,
            "pos": pos,
            "cpic_pos": pos,
            "rsid": rsid,
            "hgvs": f"{chrom}:g.{pos}{ref}>{alts[0]}",
            "ref_cpic": ref,
            "grch38_ref": ref,
            "is_key": is_key,
            "alt_variants": [
                {
                    "cpic_value": b,
                    "vcf_pos": pos,
                    "vcf_ref": ref,
                    "vcf_alt": b,
                    "is_ref": (b == ref),
                }
                for b in ordered
            ],
        }

    return [
        {
            "id": "lactase_persistence",
            "trait": "Lactase persistence",
            "model": "lactase_mcm6",
            "source": "Enattah et al. 2002 (Nat Genet 30:233); OMIM 223100",
            "evidence_level": "Strong (OMIM 223100)",
            "positions": [
                snv_pos(0, "chr2", 135851076, "rs4988235", "G", ["A", "C", "T"], True),
                snv_pos(1, "chr2", 135859184, "rs182549", "C", ["T", "A", "G"], False),
            ],
        },
        {
            "id": "eye_color",
            "trait": "Eye color (HERC2 + 6-SNP panel)",
            "model": "eye_color_6snp",
            "source": "Eiberg et al. 2008; Sturm et al. 2008 (Am J Hum Genet 82:424); Liu et al. 2009",
            "evidence_level": "Strong (HERC2 rs12913832) + Moderate (5 modifiers)",
            "positions": [
                snv_pos(0, "chr15", 28120472, "rs12913832", "A", ["G"], True),
                snv_pos(1, "chr15", 27985172, "rs1800407", "C", ["T"], False),
                snv_pos(2, "chr14", 92307319, "rs12896399", "G", ["T"], False),
                snv_pos(3, "chr5", 33951588, "rs16891982", "C", ["G"], False),
                snv_pos(4, "chr11", 89277878, "rs1393350", "G", ["A"], False),
                snv_pos(5, "chr6", 396321, "rs12203592", "C", ["T"], False),
            ],
        },
    ]


def _parse_numeric_activity(raw) -> float | None:
    if raw is None:
        return None
    s = str(raw).strip().lstrip("≥")
    return float(s) if re.fullmatch(r"-?\d+(\.\d+)?", s) else None


def build_gene_entry(
    sym: str,
    g: dict,
    locs_by_id: dict,
    defs_by_id: dict,
    al_by_gene_name: dict,
    al_by_def_id: dict,
    alv_by_def: dict,
    locs_by_gene: dict,
    ref_provider: RefWindowProvider,
) -> dict:
    """Build one normalized gene entry (positions + alleles) from CPIC metadata."""
    chrom = canonical_chrom(g["chr"])
    ploidy = "haploid" if sym == "MT-RNR1" else ("sex_x" if sym == "G6PD" else "diploid")
    model = GENE_MODELS.get(sym, "star_allele")
    lookup_method = g["lookupmethod"] or "NONE"

    ref_defs = [
        d for d in defs_by_id.values()
        if d["genesymbol"] == sym and d.get("matchesreferencesequence") and not d.get("structuralvariation")
    ]
    if len(ref_defs) != 1:
        raise RuntimeError(f"{sym}: expected 1 reference allele_definition, got {len(ref_defs)}")
    ref_def = ref_defs[0]
    ref_vals = alv_by_def[ref_def["id"]]

    gene_defs = [
        d for d in defs_by_id.values()
        if d["genesymbol"] == sym and not d.get("structuralvariation") and len(alv_by_def[d["id"]]) > 0
    ]
    loc_ids = sorted(locs_by_gene[sym], key=lambda lid: (locs_by_id[lid]["position"], lid))
    if loc_ids:
        p_min = min(locs_by_id[lid]["position"] for lid in loc_ids)
        p_max = max(locs_by_id[lid]["position"] for lid in loc_ids)
        ref_provider.prefetch_span(chrom, max(1, p_min - 100), p_max + 120)

    def _allele_row(d: dict) -> dict | None:
        return al_by_gene_name.get((sym, d["name"])) or al_by_def_id.get(d["id"])

    pos_codes = defaultdict(set)
    pos_key_score = defaultdict(float)
    for d in gene_defs:
        al = _allele_row(d)
        fq = max_pop_freq(al.get("frequency") if al else None)
        clin_fn = al.get("clinicalfunctionalstatus") if al else None
        act_val = al.get("activityvalue") if al else None
        d_map = alv_by_def[d["id"]]
        diffs = []
        for lid in loc_ids:
            code = d_map.get(lid)
            rv = ref_vals.get(lid)
            if code is not None and rv is not None:
                for atomic in expand_cpic_code(rv, code):
                    pos_codes[lid].add(atomic)
                if code != rv:
                    diffs.append(lid)
        if fq is not None and (is_functional_allele(clin_fn) or act_val in ("0.0", "0.25", "0.5")):
            for lid in diffs:
                if fq > pos_key_score[lid]:
                    pos_key_score[lid] = fq

    norm_ref_by_lid = {}
    positions_out = []
    for idx, lid in enumerate(loc_ids):
        loc = locs_by_id[lid]
        cpic_pos = int(loc["position"])
        rv = ref_vals[lid]
        ref_tuples = cpic_to_vcf_alleles(chrom, cpic_pos, rv, rv, ref_provider)
        anchor_pos, anchor_ref, _ = ref_tuples[0]
        if len(rv) == 1 and rv in IUPAC_MAP:
            norm_ref = ref_provider.get(chrom, cpic_pos, cpic_pos)
            pos_codes[lid].add(norm_ref)
        else:
            norm_ref = rv
        norm_ref_by_lid[lid] = norm_ref

        all_alts = []
        seen_triples = set()
        for code in sorted(pos_codes[lid]):
            tuples = cpic_to_vcf_alleles(chrom, cpic_pos, rv, code, ref_provider)
            for vp, vr, va in tuples:
                triple = (code, vp, vr, va)
                if triple not in seen_triples:
                    seen_triples.add(triple)
                    all_alts.append({
                        "cpic_value": code,
                        "vcf_pos": vp,
                        "vcf_ref": vr,
                        "vcf_alt": va,
                        "is_ref": (code == norm_ref),
                    })

        non_ref_positions = sorted({x["vcf_pos"] for x in all_alts if not x["is_ref"]})
        display_pos = non_ref_positions[0] if len(non_ref_positions) == 1 else cpic_pos
        if len(non_ref_positions) == 1 and display_pos != anchor_pos:
            anchor_ref = ref_provider.get(chrom, display_pos, display_pos)
            for av in all_alts:
                if av["is_ref"]:
                    av["vcf_pos"] = display_pos
                    av["vcf_ref"] = anchor_ref
                    av["vcf_alt"] = anchor_ref

        grch38_ref_exact = ref_provider.get(chrom, display_pos, display_pos + len(anchor_ref) - 1)
        if grch38_ref_exact != anchor_ref:
            raise ValueError(f"{chrom}:{display_pos}: anchor_ref={anchor_ref} differs from GRCh38={grch38_ref_exact}")

        is_key = (pos_key_score[lid] >= 0.005) or (len(loc_ids) <= 6)
        positions_out.append({
            "id": idx,
            "chrom": chrom,
            "pos": display_pos,
            "cpic_pos": cpic_pos,
            "rsid": loc.get("dbsnpid"),
            "hgvs": loc.get("chromosomelocation"),
            "ref_cpic": norm_ref,
            "grch38_ref": anchor_ref,
            "vcf_ref_base": anchor_ref,
            "is_key": is_key,
            "max_functional_maf": round(pos_key_score[lid], 6),
            "alt_variants": all_alts,
        })

    alleles_out = []
    sig_to_names = defaultdict(list)
    for d in sorted(gene_defs, key=lambda x: (not bool(x.get("matchesreferencesequence")), x["name"])):
        al = _allele_row(d)
        clin_fn = al.get("clinicalfunctionalstatus") if al else None
        act_str = al.get("activityvalue") if al else None
        act_num = float(act_str) if act_str and re.fullmatch(r"-?\d+(\.\d+)?", str(act_str)) else None
        fq = max_pop_freq(al.get("frequency") if al else None)
        d_map = alv_by_def[d["id"]]

        sparse_vars = []
        full_sig = []
        for idx, lid in enumerate(loc_ids):
            rv = ref_vals[lid]
            norm_ref = norm_ref_by_lid[lid]
            raw_code = d_map.get(lid, rv)
            accepted = expand_cpic_code(rv, raw_code)
            full_sig.append(tuple(sorted(accepted)))
            if accepted != [norm_ref]:
                sparse_vars.append({"pos_id": idx, "accepted": accepted})

        sig_to_names[tuple(full_sig)].append(d["name"])
        alleles_out.append({
            "name": d["name"],
            "is_reference": bool(d.get("matchesreferencesequence")),
            "clinical_function": clin_fn,
            "activity_value": act_num,
            "activity_label": str(act_str) if act_str is not None else None,
            "max_pop_freq": round(fq, 6) if fq is not None else None,
            "pop_freqs": pop_freqs(al.get("frequency") if al else None),
            "variants": sparse_vars,
        })

    dups = [names for names in sig_to_names.values() if len(names) > 1]
    if dups:
        raise RuntimeError(f"{sym}: duplicate allele signatures found: {dups}")

    entry_extra = {"diplotype_choice": "frequency"} if sym in FREQUENCY_CHOICE_GENES else {}
    return {
        **entry_extra,
        "gene": sym,
        "chrom": chrom,
        "model": model,
        "ploidy": ploidy,
        "lookup_method": lookup_method,
        "reference_allele": ref_def["name"],
        "positions": positions_out,
        "alleles": alleles_out,
    }


def build_bam_genes_section(
    t: dict,
    genes_meta: dict,
    locs_by_id: dict,
    defs_by_id: dict,
    al_by_gene_name: dict,
    al_by_def_id: dict,
    alv_by_def: dict,
    locs_by_gene: dict,
    ref_provider: RefWindowProvider,
    hla_pins: dict[str, dict],
) -> dict:
    """Build v2 BAM+VCF (CYP2D6) and Class I HLA metadata (without bundling IPD-IMGT sequences)."""
    cyp2d6 = build_gene_entry(
        "CYP2D6",
        genes_meta["CYP2D6"],
        locs_by_id,
        defs_by_id,
        al_by_gene_name,
        al_by_def_id,
        alv_by_def,
        locs_by_gene,
        ref_provider,
    )

    sv_type_map = {"*5": "deletion", "*36": "exon9_hybrid", "*68": "exon9_hybrid", "*61": "exon9_hybrid", "*63": "exon9_hybrid", "*13": "exon_rearrangement"}
    sv_alleles = []
    seen_sv = set()
    for d in t["allele_definition"]:
        if d["genesymbol"] == "CYP2D6" and (d.get("structuralvariation") or d["name"] in sv_type_map):
            nm = d["name"]
            if nm in seen_sv:
                continue
            seen_sv.add(nm)
            al = al_by_gene_name.get(("CYP2D6", nm))
            act_str = al.get("activityvalue") if al else None
            sv_alleles.append({
                "name": nm,
                "sv_type": sv_type_map.get(nm, "structural"),
                "clinical_function": al.get("clinicalfunctionalstatus") if al else "No function",
                "activity_value": _parse_numeric_activity(act_str),
                "activity_label": str(act_str) if act_str is not None else "0.0",
            })
    if "*36" not in seen_sv and ("CYP2D6", "*36") in al_by_gene_name:
        al = al_by_gene_name[("CYP2D6", "*36")]
        act_str = al.get("activityvalue")
        sv_alleles.append({
            "name": "*36",
            "sv_type": "exon9_hybrid",
            "clinical_function": al.get("clinicalfunctionalstatus"),
            "activity_value": _parse_numeric_activity(act_str),
            "activity_label": str(act_str) if act_str is not None else None,
        })
    sv_alleles.sort(key=lambda x: x["name"])

    dup_alleles = {}
    for a in sorted(t["allele"], key=lambda x: x["name"]):
        if a["genesymbol"] != "CYP2D6":
            continue
        nm = a["name"]
        if re.fullmatch(r"\*\d+(x2|x3|x≥3)", nm) or nm in ("*36+*10", "*36+*10x2", "*68+*4", "*68+*2", "*13+*1", "*13+*2"):
            act_str = a.get("activityvalue")
            dup_alleles[nm] = {
                "name": nm,
                "clinical_function": a.get("clinicalfunctionalstatus"),
                "activity_value": _parse_numeric_activity(act_str),
                "activity_label": str(act_str) if act_str is not None else None,
            }

    cyp2d6["sv_alleles"] = sv_alleles
    cyp2d6["duplication_alleles"] = dup_alleles
    # Canonical PharmVar *4.001 IUPAC resolutions (prefer canonical haplotype over rare suballeles when unphased)
    cyp2d6["canonical_iupac_defaults"] = {
        "42126611": "G",
        "42127941": "G",
        "42129809": "C",
        "42129819": "T",
        "42130692": "A",
    }
    cyp2d6["cnv_regions"] = {
        "target_whole": {"chrom": "chr22", "start": 42126499, "end": 42130810},
        "target_5p": {"chrom": "chr22", "start": 42128000, "end": 42130810},
        "target_ex9": {"chrom": "chr22", "start": 42126500, "end": 42126800},
        "paralog": {"chrom": "chr22", "start": 42139576, "end": 42144483},
        "controls": [
            {"name": "ALB", "chrom": "chr4", "start": 73404287, "end": 73421482},
            {"name": "EGFR_3p", "chrom": "chr7", "start": 55140000, "end": 55211628},
            {"name": "RPP30", "chrom": "chr10", "start": 90871974, "end": 90908556},
            {"name": "GAPDH", "chrom": "chr12", "start": 6534517, "end": 6538371},
            {"name": "VDR", "chrom": "chr12", "start": 47841537, "end": 47904994},
            {"name": "RPPH1", "chrom": "chr14", "start": 20343071, "end": 20343411},
            {"name": "CYP2D8P", "chrom": "chr22", "start": 42149883, "end": 42155016},
        ],
    }

    hla_bounds = {
        "HLA-A": (29936000, 29950000, "A_gen.fasta"),
        "HLA-B": (31350000, 31362000, "B_gen.fasta"),
        "HLA-C": (31265000, 31276000, "C_gen.fasta"),
    }
    hla_genes = []
    for sym in BAM_HLA_GENES:
        s_pos, e_pos, fname = hla_bounds[sym]
        pin = hla_pins[fname]
        risk_alleles = sorted({
            r["result"].split()[0]
            for r in t["gene_result"]
            if r["genesymbol"] == sym and "*" in (r.get("result") or "")
        })
        meta_g = genes_meta.get(sym, {})
        hla_genes.append({
            "gene": sym,
            "chrom": "chr6",
            "start": s_pos,
            "end": e_pos,
            "lookup_method": meta_g.get("lookupmethod") or "NONE",
            "fasta_filename": fname,
            "fasta_url": pin["url"],
            "imgt_hla_version": pin["version"],
            "fasta_sha256": pin["sha256"],
            "fasta_bytes": pin["bytes"],
            "risk_alleles": risk_alleles,
        })

    first_pin = next(iter(hla_pins.values()))
    imgt_hla = {
        "name": "IPD-IMGT/HLA Database",
        "version": first_pin["version"],
        "repository": "https://github.com/ANHIG/IMGTHLA",
        "branch": first_pin["version"].replace(".", ""),
        "license": "CC BY-ND 4.0",
        "note": "The FASTAs are downloaded unmodified and checked against fasta_sha256; they are never part of this bundle.",
    }
    return {"cyp2d6": cyp2d6, "hla_genes": hla_genes, "imgt_hla": imgt_hla}


def imgt_hla_pins(version: str, hla_dir: str | None) -> dict[str, dict]:
    """SHA-256 and size of the three IPD-IMGT/HLA genomic FASTAs of one release.

    The official mirror on GitHub (ANHIG/IMGTHLA) keeps one branch per release, named by the
    version without dots (3.66.0 -> 3660), and serves each file in seconds (the EBI server took
    minutes). The files are hashed, never copied into the bundle (their licence is CC BY-ND).
    With --hla-dir, the files already there are hashed instead of downloaded.
    """
    import hashlib
    branch = version.replace(".", "")
    pins = {}
    for fname in ("A_gen.fasta", "B_gen.fasta", "C_gen.fasta"):
        url = f"https://raw.githubusercontent.com/ANHIG/IMGTHLA/{branch}/fasta/{fname}"
        if hla_dir:
            with open(os.path.join(hla_dir, fname), "rb") as fh:
                data = fh.read()
        else:
            print(f"Downloading {url} to hash it")
            req = urllib.request.Request(url, headers={"User-Agent": "mmgen-pgx-builder/1.0"})
            with urllib.request.urlopen(req, timeout=120) as resp:
                data = resp.read()
        if not data.startswith(b">HLA:"):
            raise RuntimeError(f"{fname}: not an IPD-IMGT/HLA FASTA (expected '>HLA:' at the start)")
        pins[fname] = {"version": version, "url": url, "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}
    return pins


def build_bundle(args):
    dump_path = args.dump
    if args.download and not os.path.exists(dump_path):
        url = f"https://files.cpicpgx.org/data/database/cpic_db_dump-{args.cpic_version}.sql.gz"
        print(f"Downloading {url} -> {dump_path}")
        os.makedirs(os.path.dirname(os.path.abspath(dump_path)), exist_ok=True)
        urllib.request.urlretrieve(url, dump_path)

    if not os.path.exists(dump_path):
        raise FileNotFoundError(f"{dump_path}: CPIC dump not found (pass --dump <file> or --download)")

    print(f"Loading CPIC dump: {dump_path}")
    t = load_cpic_dump(dump_path)
    target_set = set(TARGET_GENES)

    # Exclude synonymous G6PD rs2230037 (chrX:154532439, p.Y437Y) where GRCh38
    # carries 'A' while CPIC B(ref) lists ancestral 'G', and 'Mediterranean
    # Haplotype' which differs from 'Mediterranean, Dallas, Panama, Sassari,
    # Cagliari, Birmingham' solely by rs2230037 (clinical calling uses rs5030868).
    excluded_rsids = {("G6PD", "rs2230037")}
    excluded_alleles = {("G6PD", "Mediterranean Haplotype")}

    vcf_target_set = set(TARGET_GENES)
    all_target_set = set(ALL_PGX_GENES)
    genes_meta = {g["symbol"]: g for g in t["gene"] if g["symbol"] in all_target_set}
    locs_by_id = {l["id"]: l for l in t["sequence_location"]}
    defs_by_id = {
        d["id"]: d for d in t["allele_definition"]
        if d["genesymbol"] in all_target_set and (d["genesymbol"], d["name"]) not in excluded_alleles
    }
    al_by_gene_name = {(a["genesymbol"], a["name"]): a for a in t["allele"] if a["genesymbol"] in all_target_set}
    al_by_def_id = {}
    for a in t["allele"]:
        if a["genesymbol"] in all_target_set and a["definitionid"] is not None:
            # Prefer non-duplicated single-copy allele names when multiple rows share definitionid
            if a["definitionid"] not in al_by_def_id or ("x" not in a["name"] and "+" not in a["name"]):
                al_by_def_id[a["definitionid"]] = a

    default_cache = os.path.join(tempfile.gettempdir(), "mmgen_pgx_grch38_windows_cache.json")
    cache_path = args.ref_cache or default_cache
    ref_provider = RefWindowProvider(args.fasta, cache_path)

    # Group allele_location_value per gene
    alv_by_def = defaultdict(dict)
    locs_by_gene = defaultdict(set)
    for alv in t["allele_location_value"]:
        adef = defs_by_id.get(alv["alleledefinitionid"])
        if not adef or adef.get("structuralvariation"):
            continue
        loc_id = alv["locationid"]
        loc = locs_by_id[loc_id]
        if (adef["genesymbol"], loc.get("dbsnpid")) in excluded_rsids:
            continue
        val = alv["variantallele"]
        if val is None:
            continue
        alv_by_def[adef["id"]][loc_id] = val
        locs_by_gene[adef["genesymbol"]].add(loc_id)

    gene_entries = [
        build_gene_entry(
            sym,
            genes_meta[sym],
            locs_by_id,
            defs_by_id,
            al_by_gene_name,
            al_by_def_id,
            alv_by_def,
            locs_by_gene,
            ref_provider,
        )
        for sym in TARGET_GENES
    ]

    bam_genes_payload = build_bam_genes_section(
        t,
        genes_meta,
        locs_by_id,
        defs_by_id,
        al_by_gene_name,
        al_by_def_id,
        alv_by_def,
        locs_by_gene,
        ref_provider,
        imgt_hla_pins(args.imgt_hla_version, args.hla_dir),
    )
    ref_provider.save()

    # Optional cross-check against a reference VCF (maintainer verification only)
    if args.verify_vcf and os.path.exists(args.verify_vcf):
        verify_against_vcf(gene_entries, args.verify_vcf)

    # -----------------------------------------------------------------------
    # Phenotypes table (phenotypes.json)
    # -----------------------------------------------------------------------
    gr_by_id = {r["id"]: r for r in t["gene_result"]}
    pheno_genes = {}
    for sym in ALL_PGX_GENES:
        g = genes_meta.get(sym)
        if not g:
            continue
        rules = []
        seen_rule_keys = set()
        for lk in t["gene_result_lookup"]:
            gr = gr_by_id.get(lk["phenotypeid"])
            if not gr or gr["genesymbol"] != sym:
                continue
            counts = json.loads(lk["lookupkey"])
            canon_key = json.dumps(counts, sort_keys=True)
            raw_as = str(gr["activityscore"]) if gr.get("activityscore") is not None else None
            as_val = _parse_numeric_activity(raw_as)
            entry = (canon_key, gr["result"], raw_as)
            if entry not in seen_rule_keys:
                seen_rule_keys.add(entry)
                rules.append({
                    "counts": counts,
                    "phenotype": gr["result"],
                    "activity_score": as_val,
                    "activity_label": raw_as,
                })
        pheno_genes[sym] = {
            "lookup_method": g["lookupmethod"] or "NONE",
            "rules": sorted(rules, key=lambda x: json.dumps(x["counts"], sort_keys=True)),
        }

    # -----------------------------------------------------------------------
    # Guidelines table (guidelines.json)
    # -----------------------------------------------------------------------
    drugs_by_id = {d["drugid"]: d for d in t["drug"]}
    gds_by_id = {g["id"]: g for g in t["guideline"]}
    # Pick guideline rows covering ALL_PGX_GENES, separating active conditions from 'No Result' fallbacks
    rec_entries = []
    seen_recs = set()
    for r in t["recommendation"]:
        if not r.get("lookupkey"):
            continue
        lk = json.loads(r["lookupkey"])
        if not lk:
            continue
        active = [g for g, val in lk.items() if g in all_target_set and val != "No Result"]
        if not active:
            continue
        untested = [g for g, val in lk.items() if g not in all_target_set or val == "No Result"]
        if any(g not in all_target_set and lk[g] != "No Result" for g in untested):
            continue
        text = (r.get("drugrecommendation") or "").strip()
        if not text or text.lower() == "no recommendation":
            continue
        gd = gds_by_id[r["guidelineid"]]
        dname = drugs_by_id[r["drugid"]]["name"]
        rec_key = (
            dname,
            gd["id"],
            r.get("population") or "",
            json.dumps(lk, sort_keys=True),
        )
        if rec_key in seen_recs:
            continue
        seen_recs.add(rec_key)
        active_lookup = {g: lk[g] for g in sorted(active)}
        rec_entries.append({
            "drug": dname,
            "guideline_id": str(gd["id"]),
            "guideline_name": gd["name"],
            "genes": sorted(active),
            "untested_genes": sorted(untested),
            "lookup": active_lookup,
            "population": (r.get("population") or "general").strip(),
            "classification": r.get("classification"),
            "recommendation": text,
            "url": gd.get("url"),
            "source": "CPIC",
        })

    # Add guideline metadata for Warfarin (GD 100425), which uses CPIC clinical dosing algorithm instead of static rows
    warfarin_gd = gds_by_id.get(100425)
    if warfarin_gd:
        rec_entries.append({
            "drug": "warfarin",
            "guideline_id": "100425",
            "guideline_name": warfarin_gd["name"],
            "genes": ["CYP2C9", "CYP4F2", "VKORC1"],
            "untested_genes": [],
            "lookup": {"_algorithm": "warfarin"},
            "population": "general",
            "classification": "Strong",
            "recommendation": "Use a pharmacogenetic warfarin dosing algorithm integrating CYP2C9 star alleles, VKORC1 rs9923231 (-1639G>A), and CYP4F2*3 with clinical factors (see CPIC Warfarin Guideline).",
            "url": warfarin_gd.get("url"),
            "source": "CPIC",
        })

    rec_entries.sort(key=lambda x: (x["drug"], x["population"], json.dumps(x["lookup"], sort_keys=True)))

    # -----------------------------------------------------------------------
    # Build deterministic open Traits catalog (Lactase + Eye Color 6-SNP panel)
    # -----------------------------------------------------------------------
    traits_entries = build_traits_entries()

    # -----------------------------------------------------------------------
    # Write deterministic JSON bundle + manifest + zip
    # -----------------------------------------------------------------------
    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)
    bundle_ver = f"cpic-{args.cpic_version}"
    traits_ver = "v1.0.0"
    data_ver = f"CPIC {args.cpic_version}; Open Traits {traits_ver}"
    build_date = args.date or date.today().isoformat()
    source_url = f"https://files.cpicpgx.org/data/database/cpic_db_dump-{args.cpic_version}.sql.gz"

    payloads = {
        "genes.json": {"schema_version": 1, "assembly": "GRCh38", "data_version": f"CPIC {args.cpic_version}", "genes": gene_entries},
        "bam_genes.json": {
            "schema_version": 1,
            "assembly": "GRCh38",
            "data_version": f"CPIC {args.cpic_version}",
            **bam_genes_payload,
        },
        "phenotypes.json": {"schema_version": 1, "data_version": f"CPIC {args.cpic_version}", "genes": pheno_genes},
        "guidelines.json": {"schema_version": 1, "data_version": f"CPIC {args.cpic_version}", "recommendations": rec_entries},
        "traits.json": {"schema_version": 1, "assembly": "GRCh38", "data_version": f"Open Traits {traits_ver}", "traits": traits_entries},
    }

    file_meta = {}
    for fname, obj in payloads.items():
        fpath = os.path.join(out_dir, fname)
        raw = (json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=True) + "\n").encode("utf-8")
        with open(fpath, "wb") as fh:
            fh.write(raw)
        sha = hashlib.sha256(raw).hexdigest()
        is_traits = fname == "traits.json"
        file_meta[fname] = {
            "source": "Peer-reviewed open GWAS/OMIM loci (Enattah 2002; Eiberg 2008; Sturm 2008; Liu 2009)" if is_traits else "CPIC Database Dump + GRCh38 Reference",
            "version": traits_ver if is_traits else args.cpic_version,
            "date": build_date,
            "license": "Apache-2.0" if is_traits else "CC0-1.0",
            "url": "https://omim.org/entry/223100" if is_traits else source_url,
            "sha256": sha,
            "bytes": len(raw),
        }

    manifest = {
        "bundle_format_version": 1,
        "engine": "mmgen-pgx",
        "engine_version": "1.0.0",
        "version": bundle_ver,
        "data_version": data_ver,
        "date": build_date,
        "assembly": "GRCh38",
        "source": "CPIC Database Dump + Open Traits Catalog",
        "license": "CC0-1.0 / Apache-2.0",
        "url": source_url,
        "files": file_meta,
    }
    manifest_bytes = (json.dumps(manifest, indent=2, sort_keys=True, ensure_ascii=True) + "\n").encode("utf-8")
    manifest_path = os.path.join(out_dir, "manifest.json")
    with open(manifest_path, "wb") as fh:
        fh.write(manifest_bytes)
    manifest_sha = hashlib.sha256(manifest_bytes).hexdigest()

    zip_path = args.zip or os.path.join(os.path.dirname(os.path.abspath(out_dir)), f"mmgen_pgx_tables_{bundle_ver}.zip")
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for fname in ["manifest.json", "genes.json", "bam_genes.json", "phenotypes.json", "guidelines.json", "traits.json"]:
            fpath = os.path.join(out_dir, fname)
            zi = zipfile.ZipInfo(filename=fname, date_time=(2026, 11, 15, 0, 0, 0))
            zi.compress_type = zipfile.ZIP_DEFLATED
            with open(fpath, "rb") as fh:
                zf.writestr(zi, fh.read())

    with open(zip_path, "rb") as fh:
        zip_sha = hashlib.sha256(fh.read()).hexdigest()

    total_pos = sum(len(g["positions"]) for g in gene_entries)
    total_key = sum(sum(1 for p in g["positions"] if p["is_key"]) for g in gene_entries)
    total_alleles = sum(len(g["alleles"]) for g in gene_entries)
    print(f"\nBuilt PGx bundle in: {out_dir}")
    print(f"  data_version    : {data_ver}")
    print(f"  genes           : {len(gene_entries)} (positions={total_pos}, key={total_key}, alleles={total_alleles})")
    print(f"  recommendations : {len(rec_entries)}")
    print(f"  manifest.json   : sha256={manifest_sha}")
    print(f"  zip archive     : {zip_path} (sha256={zip_sha})")


def verify_against_vcf(gene_entries: list[dict], vcf_path: str):
    """Optional maintainer cross-check against a reference positions VCF."""
    sites = defaultdict(list)
    with open(vcf_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            chrom = canonical_chrom(parts[0])
            pos = int(parts[1])
            ref = parts[3].upper()
            alts = parts[4].upper().split(",")
            sites[(chrom, pos)].append((ref, alts))
    missing = []
    for g in gene_entries:
        for p in g["positions"]:
            for av in p["alt_variants"]:
                if av["is_ref"]:
                    continue
                recs = sites.get((p["chrom"], av["vcf_pos"]), [])
                matched = any(av["vcf_ref"] == r or av["vcf_ref"].startswith(r) or r.startswith(av["vcf_ref"]) for r, _ in recs)
                if not matched:
                    missing.append((g["gene"], p["chrom"], av["vcf_pos"], av["vcf_ref"], av["vcf_alt"], av["cpic_value"], p["rsid"]))
    print(f"Reference VCF cross-check ({vcf_path}): {len(missing)} unmatched variant coordinates")
    for m in missing[:25]:
        print("  DIFF:", m)


def main():
    p = argparse.ArgumentParser(description="Build mmgen-pgx JSON table bundle from CPIC dump")
    p.add_argument("--dump", default=r"cpic_db_dump-v1.61.0.sql.gz", help="Path to CPIC SQL dump (.sql or .sql.gz)")
    p.add_argument("--download", action="store_true", help="Download CPIC dump from files.cpicpgx.org if missing")
    p.add_argument("--cpic-version", default="v1.61.0", help="CPIC version tag (e.g. v1.61.0)")
    p.add_argument("--date", default="2026-11-15", help="Build date YYYY-MM-DD for reproducibility")
    p.add_argument("--out", default=r"mmgen_tools\pgx\cpic-v1.61.0", help="Output bundle directory")
    p.add_argument("--zip", default=None, help="Output zip archive path")
    p.add_argument("--fasta", default=None, help="Optional GRCh38 FASTA with .fai")
    p.add_argument("--ref-cache", default=None, help="JSON cache for GRCh38 reference windows")
    p.add_argument("--verify-vcf", default=None, help="Optional reference VCF to verify variant coordinates")
    p.add_argument("--imgt-hla-version", default="3.66.0",
                   help="IPD-IMGT/HLA release pinned in bam_genes.json (its GitHub branch is the version without dots)")
    p.add_argument("--hla-dir", default=None,
                   help="Folder with that release's A_gen.fasta, B_gen.fasta, C_gen.fasta (hashed, not copied); "
                        "without it they are downloaded from github.com/ANHIG/IMGTHLA to be hashed")
    build_bundle(p.parse_args())


if __name__ == "__main__":
    main()
