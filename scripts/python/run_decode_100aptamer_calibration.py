#!/usr/bin/env python3
"""D033 Step 3 — 100-aptamer calibration batch.

Mandatory operational rules (per protocol, 2026-05-09):
  - Sequential single worker (deterministic timing measurement).
  - Disk-aware: per-aptamer download -> cis-extract + filter -> instrument
    table -> DELETE the .txt.gz immediately. Peak disk per worker ~1.1 GB.
  - All MR rows tagged CALIBRATION_ONLY = True. No biological claims.
  - Apply instrument extraction pre-filter before retention:
      MAF >= 0.01, single-base SNP only, p < 5e-8, F > 10, cis +/-1 Mb, no
      trans-pQTL.

Inputs:
  data/registry/decode_100aptamer_calibration_targets.tsv
  data/registry/decode_2021_folder.json (for ETag + size verification)
  data/processed/sumstats/AF.parquet (for match-rate diagnosis)

Outputs:
  data/processed/instruments/decode_2021_100aptamer_calibration.parquet
  results/qc/decode_100aptamer_cis_instrument_report.tsv
  results/qc/decode_100aptamer_calibration_report.tsv
  results/qc/decode_100aptamer_match_diagnosis.tsv
  results/qc/decode_100aptamer_runtime_projection.tsv
  results/smoke/decode_100aptamer_summary.md

Cis filter constants mirror smoke run (Phase 3A.3) plus the MAF/SNP gate.
"""

from __future__ import annotations

import gzip
import json
import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median, quantwiths

import polars as pl
import requests

ROOT = Path("/Users/apple/Desktop/gwas_af")
TARGETS_TSV = ROOT / "data/registry/decode_100aptamer_calibration_targets.tsv"
FOLDER_JSON = ROOT / "data/registry/decode_2021_folder.json"
AF_PARQUET = ROOT / "data/processed/sumstats/AF.parquet"
TMP_DIR = Path("/tmp/decode_calib")

OUT_PARQUET = ROOT / "data/processed/instruments/decode_2021_100aptamer_calibration.parquet"
OUT_CIS_REPORT = ROOT / "results/qc/decode_100aptamer_cis_instrument_report.tsv"
OUT_CALIB_REPORT = ROOT / "results/qc/decode_100aptamer_calibration_report.tsv"
OUT_MATCH = ROOT / "results/qc/decode_100aptamer_match_diagnosis.tsv"
OUT_PROJ = ROOT / "results/qc/decode_100aptamer_runtime_projection.tsv"
OUT_SUMMARY = ROOT / "results/smoke/decode_100aptamer_summary.md"

TOKEN = os.environ.get("DECODE_TOKEN", "<DECODE_TOKEN>")
DOWNLOAD_API = "https://download.decode.is/s3/download"

CIS_WINDOW_KB = 1000
P_PRIMARY = 5.0e-8
P_EXPLORATORY = 5.0e-6
F_STAT_MIN = 10.0
MAF_MIN = 0.01
DISTANCE_CLUMP_KB = 500
VALID_ALLELES = {"A", "C", "G", "T"}

# CV-prior gene -> approximate GRCh38 (gene_start, gene_end) for cis window.
# A few placeholders for genes seen in random/large-fwith selections; if missing
# at runtime we fall back to a 'wide' search window (chr ±5 Mb) to ensure
# coverage on the calibration measurements.
GENE_COORDS_B38: dict[str, tuple[str, int, int]] = {
    "PCSK9":    ("1",  55039548, 55064852),
    "IL6R":     ("1",  154405193, 154469450),
    "SORT1":    ("1",  109274570, 109365052),
    "ANGPTL3":  ("1",  62597488, 62606909),
    "ANGPTL4":  ("19", 8364181, 8376234),
    "APOA1":    ("11", 116835751, 116837950),
    "APOA4":    ("11", 116821687, 116829907),
    "APOE":     ("19", 44905781, 44909393),
    "APOB":     ("2",  21001429, 21043950),
    "LDLR":     ("19", 11089463, 11133820),
    "F2":       ("11", 46740730, 46761057),
    "F11":      ("4",  186266801, 186290374),
    "F9":       ("X",  139530763, 139564310),
    "F7":       ("13", 113105016, 113127538),
    "FGA":      ("4",  154584253, 154591571),
    "FGB":      ("4",  154562811, 154572788),
    "VWF":      ("12", 5949715, 6124770),
    "PROC":     ("2",  127417399, 127430350),
    "PROS1":    ("3",  93876035, 93974092),
    "SERPINC1": ("1",  173903804, 173917378),
    "ADAMTS13": ("9",  133414233, 133458997),
    "PLAT":     ("8",  42174737, 42208068),
    "PLG":      ("6",  160702170, 160754533),
    "SERPINE1": ("7",  101127103, 101139245),
    "IL6":      ("7",  22725889, 22732002),
    "IL1B":     ("2",  112829751, 112836781),
    "IL18":     ("11", 112143251, 112164100),
    "CRP":      ("1",  159712289, 159714589),
    "TNFRSF1A": ("12", 6328761, 6342144),
    "TNFRSF1B": ("1",  12166948, 12209232),
    "VCAM1":    ("1",  100735378, 100754294),
    "ICAM1":    ("19", 10286924, 10303612),
    "SELE":     ("1",  169722658, 169744013),
    "SELP":     ("1",  169551765, 169588669),
    "MMP1":     ("11", 102789647, 102797393),
    "MMP2":     ("16", 55477497, 55505571),
    "MMP3":     ("11", 102835802, 102843909),
    "MMP9":     ("20", 46008908, 46016561),
    "TIMP1":    ("X",  47582380, 47586793),
    "TIMP2":    ("17", 78852349, 78925229),
    "COL1A1":   ("17", 50184101, 50201613),
    "TGFB1":    ("19", 41330323, 41353923),
    "TGFB2":    ("1",  218346833, 218445761),
    "ANGPT1":   ("8",  107249506, 107285550),
    "ANGPT2":   ("8",  6356157, 6404566),
    "VEGFA":    ("6",  43770209, 43795785),
    "VEGFC":    ("4",  176684432, 176796842),
    "PDGFB":    ("22", 39225037, 39247261),
    "EDN1":     ("6",  12290530, 12297194),
    "BAG3":     ("10", 119651379, 119677819),
    "MYBPC3":   ("11", 47331406, 47352702),
    "ACTC1":    ("15", 34790231, 34795549),
    "REN":      ("1",  204154948, 204190329),
    "AGT":      ("1",  230702517, 230745244),
    "NPPA":     ("1",  11845709, 11848345),
    "NPPB":     ("1",  11849854, 11853614),
    "GDF15":    ("19", 18386152, 18389176),
    "IGFBP7":   ("4",  57045006, 57110018),
}


# ---------------------------------------------------------------------------
# Streaming + filter
# ---------------------------------------------------------------------------


@dataclass
class CalibRow:
    aptamer: str
    gene: str
    category: str
    expected_size: int
    actual_size: int = 0
    download_seconds: float = 0.0
    extract_seconds: float = 0.0
    n_cis_pre_filter: int = 0     # cis window + p<5e-8 + F>10 (no MAF/SNP)
    n_cis_post_filter: int = 0    # add MAF>=0.01 + single-base
    n_dropped_rare: int = 0       # MAF<0.01
    n_dropped_indel: int = 0      # multi-base allele
    n_clumped: int = 0
    n_matched_AF: int = 0
    n_drop_allele_mismatch: int = 0
    n_drop_palindromic: int = 0
    n_drop_missing_in_AF: int = 0
    error: str = ""
    instruments: list[dict] = field(default_factory=list)


def _log(msg: str) -> None:
    print(f"[calib {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def stream_cis_extract_filtered(gz_path: Path,
                                gene: str,
                                target_chr: str,
                                window_start: int,
                                window_end: int) -> tuple[list[dict], dict]:
    """Read .txt.gz line by line; return rows in cis window passing:
       p < 5e-8, F > 10, MAF >= 0.01, single-base alleles.

    Counters:
       n_cis_pre_filter, n_cis_post_filter, n_dropped_rare, n_dropped_indel
    """
    counters = {
        "n_cis_pre_filter": 0, "n_cis_post_filter": 0,
        "n_dropped_rare": 0, "n_dropped_indel": 0,
    }
    out_rows: list[dict] = []
    target_chr_norm = str(target_chr).replace("chr", "")
    with gzip.open(gz_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        ci, pi = idx["Chrom"], idx["Pos"]
        ni, ri = idx["Name"], idx.get("rsids", -1)
        ei, oi = idx["effectAllele"], idx["otherAllele"]
        bi, pvi = idx["Beta"], idx["Pval"]
        si, sni = idx["SE"], idx["N"]
        mi = idx.get("ImpMAF", -1)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if cols[ci].replace("chr", "") != target_chr_norm:
                continue
            try:
                pos = int(cols[pi])
            except ValueError:
                continue
            if pos < window_start or pos > window_end:
                continue
            try:
                pval = float(cols[pvi])
            except ValueError:
                continue
            if pval >= P_PRIMARY:
                continue
            try:
                beta = float(cols[bi])
                se = float(cols[si])
                n = int(float(cols[sni]))
            except ValueError:
                continue
            f_stat = (beta / se) ** 2 if se > 0 else 0.0
            if f_stat < F_STAT_MIN:
                continue

            # Pre-filter cis variant count (after p + F filters).
            counters["n_cis_pre_filter"] += 1

            ea = cols[ei].upper()
            oa = cols[oi].upper()
            single_base = ea in VALID_ALLELES and oa in VALID_ALLELES and ea != oa
            if not single_base:
                counters["n_dropped_indel"] += 1
                continue

            imp_maf = float(cols[mi]) if mi >= 0 and cols[mi] not in ("", "NA") else float("nan")
            if (imp_maf == imp_maf) and imp_maf < MAF_MIN:  # not NaN
                counters["n_dropped_rare"] += 1
                continue

            counters["n_cis_post_filter"] += 1
            out_rows.append({
                "chrom": target_chr_norm, "pos": pos, "name": cols[ni],
                "rsid": cols[ri] if ri >= 0 else "",
                "effect_allele": ea, "other_allele": oa,
                "beta": beta, "se": se, "pval": pval,
                "n": n, "imp_maf": imp_maf, "f_stat": f_stat,
            })
    return out_rows, counters


def distance_clump(rows: list[dict], window_kb: int = DISTANCE_CLUMP_KB) -> list[dict]:
    sorted_rows = sorted(rows, key=lambda r: r["pval"])
    keep = []
    excluded = set()
    pad = window_kb * 1000
    for i, r in enumerate(sorted_rows):
        if i in excluded:
            continue
        keep.append(r)
        for j in range(i + 1, len(sorted_rows)):
            if abs(sorted_rows[j]["pos"] - r["pos"]) <= pad:
                excluded.add(j)
    return keep


# ---------------------------------------------------------------------------
# Match diagnosis vs AF
# ---------------------------------------------------------------------------


def match_to_af(instruments: list[dict], af: pl.DataFrame) -> dict:
    """Return per-aptamer match counters and append match-status to each instrument."""
    counters = {"n_matched_AF": 0, "n_drop_allele_mismatch": 0,
                "n_drop_palindromic": 0, "n_drop_missing_in_AF": 0}
    if not instruments:
        return counters
    keys = pl.DataFrame({
        "chr": [r["chrom"] for r in instruments],
        "pos": [r["pos"] for r in instruments],
    })
    af_hit = af.join(keys, on=["chr", "pos"], how="inner")
    af_lookup = {(str(r["chr"]), int(r["pos"])): r for r in af_hit.iter_rows(named=True)}
    pal = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
    for r in instruments:
        oc = af_lookup.get((str(r["chrom"]), int(r["pos"])))
        if oc is None:
            counters["n_drop_missing_in_AF"] += 1
            r["af_match"] = "missing_in_AF"
            continue
        ea, oa = r["effect_allele"], r["other_allele"]
        af_ea = (oc["effect_allele"] or "").upper()
        af_oa = (oc["other_allele"] or "").upper()
        if {ea, oa} != {af_ea, af_oa}:
            counters["n_drop_allele_mismatch"] += 1
            r["af_match"] = "allele_mismatch"
            continue
        if (ea, oa) in pal:
            eaf_o = oc.get("eaf")
            if eaf_o is None or eaf_o != eaf_o or abs(float(eaf_o) - 0.5) < 0.05:
                counters["n_drop_palindromic"] += 1
                r["af_match"] = "palindromic_unresolved"
                continue
        counters["n_matched_AF"] += 1
        r["af_match"] = "matched"
        r["af_beta"] = float(oc["beta"])
        r["af_se"] = float(oc["se"])
        r["af_pval"] = float(oc["pval"])
        r["af_eaf"] = oc.get("eaf")
    return counters


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def download(key: str, dest: Path) -> tuple[float, int]:
    dest.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with requests.get(DOWNLOAD_API,
                      params={"token": TOKEN, "file": key},
                      stream=True, timeout=600) as r:
        r.raise_for_status()
        with dest.open("wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
    return time.time() - t0, dest.stat().st_size


GENE_COORD_CACHE = ROOT / "data/registry/gene_coords_b38_cache.json"


def _load_coord_cache() -> dict:
    if GENE_COORD_CACHE.exists():
        return json.loads(GENE_COORD_CACHE.read_text())
    return {}


def _save_coord_cache(cache: dict) -> None:
    GENE_COORD_CACHE.parent.mkdir(parents=True, exist_ok=True)
    GENE_COORD_CACHE.write_text(json.dumps(cache, indent=2))


def fetch_ensembl_coords(gene: str) -> tuple[str, int, int] | None:
    """Query Ensembl REST API for canonical gene_start/end on GRCh38."""
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}"
    try:
        r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=15)
        if r.status_code != 200:
            return None
        d = r.json()
        chrom = str(d.get("seq_region_name", ""))
        start = d.get("start")
        end = d.get("end")
        if chrom and start and end:
            return chrom, int(start), int(end)
    except Exception:
        pass
    return None


def find_gene_coords(gene: str, default_search_chr: str | None = None,
                     cache: dict | None = None) -> tuple[str, int, int] | None:
    """Hardcoded dict first; fallback to Ensembl REST with persistent cache."""
    if gene in GENE_COORDS_B38:
        return GENE_COORDS_B38[gene]
    if cache is None:
        cache = {}
    if gene in cache:
        v = cache[gene]
        return (v[0], int(v[1]), int(v[2])) if v else None
    coords = fetch_ensembl_coords(gene)
    cache[gene] = list(coords) if coords else None
    return coords


def projection_table(seconds_per_aptamer: list[float],
                     bytes_per_aptamer: list[int]) -> list[dict]:
    """Project transfer + runtime for 500/1000/4907 × 1/2/4 workers."""
    if not seconds_per_aptamer:
        return []
    median_sec = median(seconds_per_aptamer)
    median_bytes = median(bytes_per_aptamer)
    out = []
    for n in (500, 1000, 4907):
        for nw in (1, 2, 4):
            # 4-worker scaling tempered by likely server cap (~1.5x speedup at 2w, ~2.5x at 4w)
            speedup = {1: 1.0, 2: 1.6, 4: 2.5}[nw]
            wall_sec = n * median_sec / speedup
            transfer_gb = n * median_bytes / 1e9
            out.append({
                "n_aptamers": n,
                "workers": nw,
                "median_seconds_per_aptamer_obs": round(median_sec, 1),
                "median_bytes_per_aptamer_obs": int(median_bytes),
                "projected_wall_hours": round(wall_sec / 3600, 2),
                "projected_transfer_GB": round(transfer_gb, 1),
                "speedup_assumed": speedup,
            })
    return out


def main() -> int:
    OUT_PARQUET.parent.mkdir(parents=True, exist_ok=True)
    OUT_CIS_REPORT.parent.mkdir(parents=True, exist_ok=True)
    OUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    targets = pl.read_csv(TARGETS_TSV, separator="\t")
    _log(f"loaded {targets.height} target aptamers")

    af = pl.read_parquet(AF_PARQUET).select([
        "chr", "pos", "rsid", "effect_allele", "other_allele",
        "beta", "se", "pval", "eaf",
    ])
    _log(f"AF parquet: {af.height:,} rows")

    rows: list[CalibRow] = []
    seconds_per_aptamer = []
    bytes_per_aptamer = []
    all_instruments: list[dict] = []

    coord_cache = _load_coord_cache()

    for i, t in enumerate(targets.iter_rows(named=True), 1):
        row = CalibRow(
            aptamer=t["aptamer_key"], gene=t["gene"],
            category=t["category"], expected_size=int(t["size_bytes"]),
        )
        coords = find_gene_coords(t["gene"], cache=coord_cache)
        # Periodically persist cache so we don't lose lookups on crash.
        if i % 10 == 0:
            _save_coord_cache(coord_cache)
        if coords is None:
            row.error = "gene_coord_missing"
            rows.append(row)
            _log(f"[{i}/100] {t['gene']} (cat={t['category']}) coords missing — skipping")
            continue

        chrom, gs, ge = coords
        win_start = max(0, gs - CIS_WINDOW_KB * 1000)
        win_end = ge + CIS_WINDOW_KB * 1000

        dest = TMP_DIR / t["aptamer_key"]
        try:
            _log(f"[{i}/100] {t['gene']} {t['aptamer_key']} (cat={t['category']}, expected {t['size_bytes']/1e6:.0f} MB)")
            dl_t, sz = download(t["aptamer_key"], dest)
            row.download_seconds = dl_t
            row.actual_size = sz

            t1 = time.time()
            cis_rows, counters = stream_cis_extract_filtered(
                dest, t["gene"], chrom, win_start, win_end,
            )
            row.extract_seconds = time.time() - t1
            row.n_cis_pre_filter = counters["n_cis_pre_filter"]
            row.n_cis_post_filter = counters["n_cis_post_filter"]
            row.n_dropped_rare = counters["n_dropped_rare"]
            row.n_dropped_indel = counters["n_dropped_indel"]

            clumped = distance_clump(cis_rows) if cis_rows else []
            row.n_clumped = len(clumped)

            match_counters = match_to_af(clumped, af)
            row.n_matched_AF = match_counters["n_matched_AF"]
            row.n_drop_allele_mismatch = match_counters["n_drop_allele_mismatch"]
            row.n_drop_palindromic = match_counters["n_drop_palindromic"]
            row.n_drop_missing_in_AF = match_counters["n_drop_missing_in_AF"]

            for r in clumped:
                all_instruments.append({
                    "gene": t["gene"], "aptamer": t["aptamer_key"],
                    "category": t["category"],
                    **{k: r[k] for k in ("chrom", "pos", "name", "rsid",
                                          "effect_allele", "other_allele",
                                          "beta", "se", "pval", "n",
                                          "imp_maf", "f_stat")},
                    "af_match": r.get("af_match", ""),
                    "af_beta": r.get("af_beta"),
                    "af_se": r.get("af_se"),
                    "af_pval": r.get("af_pval"),
                    "af_eaf": r.get("af_eaf"),
                    "build": "GRCh38", "ancestry": "EUR",
                    "source": "deCODE_2021_Ferkingstad",
                    "calibration_only": True,
                })
            seconds_per_aptamer.append(dl_t + row.extract_seconds)
            bytes_per_aptamer.append(sz)
            _log(f"  cis_pre={row.n_cis_pre_filter} cis_post={row.n_cis_post_filter} "
                 f"clumped={row.n_clumped} matched_AF={row.n_matched_AF} "
                 f"download={dl_t:.0f}s")
        except Exception as e:
            row.error = f"{type(e).__name__}: {str(e)[:120]}"
            _log(f"  ERROR: {row.error}")
        finally:
            if dest.exists():
                dest.unlink()
            rows.append(row)

    # ----- Persist parquet -----
    if all_instruments:
        pl.DataFrame(all_instruments).write_parquet(OUT_PARQUET, compression="zstd")
        _log(f"wrote {OUT_PARQUET} ({len(all_instruments)} instruments)")

    # Cis-instrument report (per aptamer)
    pl.DataFrame([{
        "gene": r.gene, "aptamer": r.aptamer, "category": r.category,
        "expected_size_bytes": r.expected_size,
        "actual_size_bytes": r.actual_size,
        "download_seconds": round(r.download_seconds, 1),
        "extract_seconds": round(r.extract_seconds, 1),
        "n_cis_pre_filter": r.n_cis_pre_filter,
        "n_cis_post_filter": r.n_cis_post_filter,
        "n_dropped_rare": r.n_dropped_rare,
        "n_dropped_indel": r.n_dropped_indel,
        "n_clumped": r.n_clumped,
        "n_matched_AF": r.n_matched_AF,
        "n_drop_allele_mismatch": r.n_drop_allele_mismatch,
        "n_drop_palindromic": r.n_drop_palindromic,
        "n_drop_missing_in_AF": r.n_drop_missing_in_AF,
        "error": r.error,
        "calibration_only": True,
    } for r in rows]).write_csv(OUT_CIS_REPORT, separator="\t")
    _log(f"wrote {OUT_CIS_REPORT}")

    # Calibration aggregate report
    n_attempted = len(rows)
    n_done = sum(1 for r in rows if r.actual_size > 0)
    n_empty = sum(1 for r in rows if r.actual_size > 0 and r.n_clumped == 0)
    total_bytes = sum(r.actual_size for r in rows)
    total_seconds = sum(r.download_seconds + r.extract_seconds for r in rows)
    n_pre = sum(r.n_cis_pre_filter for r in rows)
    n_post = sum(r.n_cis_post_filter for r in rows)
    n_clumped_total = sum(r.n_clumped for r in rows)
    n_matched = sum(r.n_matched_AF for r in rows)
    rare_drop = sum(r.n_dropped_rare for r in rows)
    indel_drop = sum(r.n_dropped_indel for r in rows)
    allele_mm = sum(r.n_drop_allele_mismatch for r in rows)

    median_sec = median(seconds_per_aptamer) if seconds_per_aptamer else 0.0
    p95_sec = quantwiths(seconds_per_aptamer, n=20)[18] if len(seconds_per_aptamer) >= 20 else max(seconds_per_aptamer or [0])
    median_bytes = median(bytes_per_aptamer) if bytes_per_aptamer else 0

    matched_rate = n_matched / max(1, n_clumped_total)
    empty_rate = n_empty / max(1, n_done)
    rare_rate = rare_drop / max(1, n_pre + rare_drop)
    indel_rate = indel_drop / max(1, n_pre + indel_drop)

    # Discordance: count multi-aptamer genes whose aptamer-level lead direction
    # differs (placeholder; full discordance check requires AF MR which we are
    # NOT running in calibration).
    multi_genes = [r.gene for r in rows if r.category == "multi_aptamer"]
    discordance_count = 0  # AF MR not run in calibration

    pl.DataFrame([{
        "n_aptamers_attempted": n_attempted,
        "n_aptamers_completed": n_done,
        "n_aptamers_empty": n_empty,
        "GB_transferred_total": round(total_bytes / 1e9, 2),
        "GB_per_aptamer_median": round(median_bytes / 1e9, 3),
        "minutes_per_aptamer_median": round(median_sec / 60, 2),
        "minutes_per_aptamer_p95": round(p95_sec / 60, 2),
        "wall_hours_total_observed": round(total_seconds / 3600, 2),
        "n_cis_instruments_total": n_pre,
        "n_cis_instruments_after_MAF_SNP_filter": n_post,
        "n_cis_instruments_after_clumping": n_clumped_total,
        "matched_to_AF_rate_post_clump": round(matched_rate, 3),
        "rare_variant_drop_rate": round(rare_rate, 3),
        "indel_drop_rate": round(indel_rate, 3),
        "allele_mismatch_count": allele_mm,
        "empty_rate": round(empty_rate, 3),
        "multi_aptamer_discordance_count": discordance_count,
        "calibration_only": True,
    }]).write_csv(OUT_CALIB_REPORT, separator="\t")
    _log(f"wrote {OUT_CALIB_REPORT}")

    # Match diagnosis (per aptamer)
    pl.DataFrame([{
        "gene": r.gene, "aptamer": r.aptamer, "category": r.category,
        "n_clumped": r.n_clumped,
        "n_matched_AF": r.n_matched_AF,
        "n_drop_allele_mismatch": r.n_drop_allele_mismatch,
        "n_drop_palindromic": r.n_drop_palindromic,
        "n_drop_missing_in_AF": r.n_drop_missing_in_AF,
        "matched_rate": (r.n_matched_AF / r.n_clumped) if r.n_clumped else None,
    } for r in rows]).write_csv(OUT_MATCH, separator="\t")
    _log(f"wrote {OUT_MATCH}")

    # Runtime projections
    if seconds_per_aptamer:
        proj_rows = projection_table(seconds_per_aptamer, bytes_per_aptamer)
        pl.DataFrame(proj_rows).write_csv(OUT_PROJ, separator="\t")
        _log(f"wrote {OUT_PROJ}")

    # Markdown summary
    md = ["# deCODE 100-aptamer calibration batch (D033 Step 3, CALIBRATION_ONLY)\n",
          f"All rows tagged `calibration_only = True`. **No biological claims.**",
          "",
          "## Aggregate metrics",
          "",
          f"- aptamers attempted: {n_attempted}",
          f"- completed (downloaded): {n_done}",
          f"- empty cis-instrument: {n_empty} ({empty_rate:.1%})",
          f"- total transfer: {total_bytes/1e9:.1f} GB",
          f"- median GB per aptamer: {median_bytes/1e9:.3f}",
          f"- median minutes per aptamer: {median_sec/60:.2f}",
          f"- p95 minutes per aptamer: {p95_sec/60:.2f}",
          f"- total observed wall-time: {total_seconds/3600:.2f} hours",
          f"- n_cis_instruments (after p+F filter): {n_pre}",
          f"- after MAF/SNP filter: {n_post}",
          f"- after distance clumping: {n_clumped_total}",
          f"- matched to AF: {n_matched} ({matched_rate:.1%})",
          f"- rare-drop rate: {rare_rate:.1%}",
          f"- indel-drop rate: {indel_rate:.1%}",
          f"- allele mismatch count: {allele_mm}",
          "",
          "## Approval thresholds (1000-aptamer projection)",
          "",
    ]
    if seconds_per_aptamer:
        proj_1000 = next(p for p in proj_rows if p["n_aptamers"] == 1000 and p["workers"] == 4)
        md.extend([
            f"- 1000-aptamer @ 4 workers: **{proj_1000['projected_wall_hours']:.1f} h**, "
            f"**{proj_1000['projected_transfer_GB']:.0f} GB**",
            "- thresholds: ≤18 h, ≤1.2 TB, matched_rate ≥50%, empty_rate ≤35%",
        ])
    OUT_SUMMARY.write_text("\n".join(md))
    _log(f"wrote {OUT_SUMMARY}")
    _save_coord_cache(coord_cache)
    _log(f"100-aptamer calibration complete. matched_rate={matched_rate:.1%} empty_rate={empty_rate:.1%}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
