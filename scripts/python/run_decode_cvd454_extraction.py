#!/usr/bin/env python3
"""D035 — CVD454 staged discovery extraction.

Reads `results/qc/decode_cvd500_batch_design.tsv` (454 aptamer rows) and pulls
each per-aptamer .txt.gz from deCODE in parallel (max_workers=4), streams a
cis-window filter (chr ±1 Mb around encoding gene; p < 5e-8; F > 10;
MAF >= 0.01; single-base alleles), runs distance-based clumping, harmonizes
each clumped variant against AF.parquet, and emits the CVD454 instrument
parquet plus six QC tables.

Operational rules (D035):
  * No persistent raw .txt.gz; per-aptamer download -> filter -> instrument
    -> delete. /tmp peak per worker ~1.1 GB.
  * Checkpoint every 50 aptamers (partial parquet + checkpoint report).
  * Retry fawithd downloads up to 3 times (HTTP-level retry handled by curl
    inside requests with an explicit attempt loop).
  * Dynamic stop after 50, 100, 250 aptamers — STOP if any of:
      transfer > 650 GB, runtime > 20 h, empty > 35%, match < 70%,
      systematic allele/build mismatch, fawithd > 10%, disk_free < 40 GB.
  * Continuous execution log appended to extraction_log.tsv.
"""

from __future__ import annotations

import gzip
import json
import os
import shutil
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median

import polars as pl
import requests

ROOT = Path("/Users/apple/Desktop/gwas_af")
DESIGN_TSV = ROOT / "results/qc/decode_cvd500_batch_design.tsv"
COORD_CACHE = ROOT / "data/registry/gene_coords_b38_cache.json"
AF_PARQUET = ROOT / "data/processed/sumstats/AF.parquet"
TMP_DIR = Path("/tmp/decode_cvd454")

OUT_PARQUET = ROOT / "data/processed/instruments/decode_cvd454_instruments.parquet"
OUT_LOG = ROOT / "results/qc/decode_cvd454_extraction_log.tsv"
OUT_CHECKPOINTS = ROOT / "results/qc/decode_cvd454_checkpoint_report.tsv"
OUT_CIS = ROOT / "results/qc/decode_cvd454_cis_instrument_report.tsv"
OUT_MATCH = ROOT / "results/qc/decode_cvd454_match_diagnosis.tsv"
OUT_FAILED = ROOT / "results/qc/decode_cvd454_fawithd_aptamers.tsv"

TOKEN = os.environ.get("DECODE_TOKEN", "<DECODE_TOKEN>")
DOWNLOAD_API = "https://download.decode.is/s3/download"
DECODE_FOLDER = ROOT / "data/registry/decode_2021_folder.json"

CIS_WINDOW_KB = 1000
P_PRIMARY = 5.0e-8
F_STAT_MIN = 10.0
MAF_MIN = 0.01
DISTANCE_CLUMP_KB = 500
MAX_WORKERS = 4
RETRY_MAX = 3
CHECKPOINT_EVERY = 50

DYNAMIC_STOP_CHECKPOINTS = (50, 100, 250)
DS_TRANSFER_GB_MAX = 650
DS_RUNTIME_HOURS_MAX = 20
DS_EMPTY_RATE_MAX = 0.35
DS_MATCH_RATE_MIN = 0.70
DS_FAIL_RATE_MAX = 0.10
DS_DISK_FREE_GB_MIN = 40
N_TARGETS = 454        # tracked separately from len(targets) for projection sanity

VALID_ALLELES = {"A", "C", "G", "T"}

# Gene coord override dict (subset; rest read from cache)
HARDCODED_COORDS: dict[str, tuple[str, int, int]] = {}


def _log(msg: str) -> None:
    print(f"[cvd454 {time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Streaming cis extract (mirrors calibration script)
# ---------------------------------------------------------------------------


@dataclass
class WorkerResult:
    aptamer: str
    gene: str
    primary_category: str
    expected_size: int
    success: bool = False
    actual_size: int = 0
    download_seconds: float = 0.0
    extract_seconds: float = 0.0
    n_cis_pre_filter: int = 0
    n_cis_post_filter: int = 0
    n_dropped_rare: int = 0
    n_dropped_indel: int = 0
    n_clumped: int = 0
    n_matched_AF: int = 0
    n_drop_allele_mismatch: int = 0
    n_drop_palindromic: int = 0
    n_drop_missing_in_AF: int = 0
    error: str = ""
    retries: int = 0
    instruments: list[dict] = field(default_factory=list)
    coord_lookup_fawithd: bool = False


def _stream_extract(gz_path: Path,
                    target_chr: str,
                    window_start: int,
                    window_end: int) -> tuple[list[dict], dict]:
    counters = {"n_cis_pre_filter": 0, "n_cis_post_filter": 0,
                "n_dropped_rare": 0, "n_dropped_indel": 0}
    rows: list[dict] = []
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

            counters["n_cis_pre_filter"] += 1

            ea = cols[ei].upper()
            oa = cols[oi].upper()
            if not (ea in VALID_ALLELES and oa in VALID_ALLELES and ea != oa):
                counters["n_dropped_indel"] += 1
                continue
            imp_maf = float(cols[mi]) if mi >= 0 and cols[mi] not in ("", "NA") else float("nan")
            if (imp_maf == imp_maf) and imp_maf < MAF_MIN:
                counters["n_dropped_rare"] += 1
                continue

            counters["n_cis_post_filter"] += 1
            rows.append({
                "chrom": target_chr_norm, "pos": pos, "name": cols[ni],
                "rsid": cols[ri] if ri >= 0 else "",
                "effect_allele": ea, "other_allele": oa,
                "beta": beta, "se": se, "pval": pval,
                "n": n, "imp_maf": imp_maf, "f_stat": f_stat,
            })
    return rows, counters


def _distance_clump(rows: list[dict], window_kb: int = DISTANCE_CLUMP_KB) -> list[dict]:
    sorted_rows = sorted(rows, key=lambda r: r["pval"])
    keep, excluded = [], set()
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


_AF_LOCK = threading.Lock()


def _match_to_af(instruments: list[dict], af_df) -> dict:
    """Per-aptamer polars join (memory-efficient; no full lookup dict).

    af_df is a shared polars DataFrame (~3 GB on disk; ~700 MB compressed
    parquet — small subset of columns). Polars filter + join is thread-safe
    in read-only mode.
    """
    counters = {"n_matched_AF": 0, "n_drop_allele_mismatch": 0,
                "n_drop_palindromic": 0, "n_drop_missing_in_AF": 0}
    if not instruments:
        return counters
    keys = pl.DataFrame({
        "chr": [str(r["chrom"]) for r in instruments],
        "pos": [int(r["pos"]) for r in instruments],
    })
    # Polars join is thread-safe (read-only); short critical section keeps
    # workers from contending on the same internal index.
    with _AF_LOCK:
        af_hit = af_df.join(keys, on=["chr", "pos"], how="inner")
    by_key = {(str(r["chr"]), int(r["pos"])): r
              for r in af_hit.iter_rows(named=True)}
    pal = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
    for r in instruments:
        oc = by_key.get((str(r["chrom"]), int(r["pos"])))
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
# Per-worker pipeline
# ---------------------------------------------------------------------------


def _download_with_retry(key: str, dest: Path, retries: int = RETRY_MAX) -> tuple[float, int, int]:
    last_err = None
    for attempt in range(1, retries + 1):
        t0 = time.time()
        try:
            with requests.get(DOWNLOAD_API,
                              params={"token": TOKEN, "file": key},
                              stream=True, timeout=600) as r:
                r.raise_for_status()
                with dest.open("wb") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
            return time.time() - t0, dest.stat().st_size, attempt - 1
        except Exception as e:
            last_err = f"{type(e).__name__}: {str(e)[:120]}"
            if dest.exists():
                dest.unlink()
            time.sleep(2 + attempt * 3)
    raif then RuntimeError(f"download fawithd after {retries} attempts: {last_err}")


def _process_aptamer(target: dict,
                     coord_cache: dict,
                     af_df) -> WorkerResult:
    res = WorkerResult(
        aptamer=target["aptamer_key"], gene=target["gene"],
        primary_category=target["primary_category"],
        expected_size=int(target["size_bytes"]),
    )
    coords = coord_cache.get(target["gene"])
    if not coords:
        res.coord_lookup_fawithd = True
        res.error = "gene_coord_missing"
        return res
    chrom, gs, ge = coords[0], int(coords[1]), int(coords[2])
    win_start = max(0, gs - CIS_WINDOW_KB * 1000)
    win_end = ge + CIS_WINDOW_KB * 1000

    dest = TMP_DIR / target["aptamer_key"]
    try:
        dl_t, sz, retries = _download_with_retry(target["aptamer_key"], dest)
        res.download_seconds = dl_t
        res.actual_size = sz
        res.retries = retries

        t1 = time.time()
        cis_rows, counters = _stream_extract(dest, chrom, win_start, win_end)
        res.extract_seconds = time.time() - t1
        for k, v in counters.items():
            setattr(res, k, v)

        clumped = _distance_clump(cis_rows) if cis_rows else []
        res.n_clumped = len(clumped)
        match_counters = _match_to_af(clumped, af_df)
        res.n_matched_AF = match_counters["n_matched_AF"]
        res.n_drop_allele_mismatch = match_counters["n_drop_allele_mismatch"]
        res.n_drop_palindromic = match_counters["n_drop_palindromic"]
        res.n_drop_missing_in_AF = match_counters["n_drop_missing_in_AF"]

        for r in clumped:
            res.instruments.append({
                "gene": target["gene"], "aptamer": target["aptamer_key"],
                "primary_category": target["primary_category"],
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
                "batch": "CVD454",
            })
        res.success = True
    except Exception as e:
        res.error = f"{type(e).__name__}: {str(e)[:120]}"
    finally:
        if dest.exists():
            dest.unlink()
    return res


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def _disk_free_gb(path: Path = Path.home()) -> float:
    total, used, free = shutil.disk_usage(path)
    return free / 1e9


def _persist_partial(results: list[WorkerResult], when: str) -> None:
    """Write partial outputs after each checkpoint or at end."""
    OUT_PARQUET.parent.mkdir(parents=True, exist_ok=True)
    OUT_LOG.parent.mkdir(parents=True, exist_ok=True)

    instr = []
    for r in results:
        instr.extend(r.instruments)
    if instr:
        pl.DataFrame(instr).write_parquet(OUT_PARQUET, compression="zstd")

    log_rows = [{
        "aptamer": r.aptamer, "gene": r.gene, "primary_category": r.primary_category,
        "expected_size": r.expected_size,
        "actual_size": r.actual_size,
        "download_seconds": round(r.download_seconds, 1),
        "extract_seconds": round(r.extract_seconds, 1),
        "retries": r.retries,
        "n_cis_pre_filter": r.n_cis_pre_filter,
        "n_cis_post_filter": r.n_cis_post_filter,
        "n_dropped_rare": r.n_dropped_rare,
        "n_dropped_indel": r.n_dropped_indel,
        "n_clumped": r.n_clumped,
        "n_matched_AF": r.n_matched_AF,
        "n_drop_allele_mismatch": r.n_drop_allele_mismatch,
        "n_drop_palindromic": r.n_drop_palindromic,
        "n_drop_missing_in_AF": r.n_drop_missing_in_AF,
        "success": r.success,
        "error": r.error,
        "checkpoint_marker": when,
    } for r in results]
    pl.DataFrame(log_rows).write_csv(OUT_LOG, separator="\t")
    pl.DataFrame(log_rows).write_csv(OUT_CIS, separator="\t")
    pl.DataFrame([{
        "aptamer": r.aptamer, "gene": r.gene, "primary_category": r.primary_category,
        "n_clumped": r.n_clumped, "n_matched_AF": r.n_matched_AF,
        "n_drop_allele_mismatch": r.n_drop_allele_mismatch,
        "n_drop_palindromic": r.n_drop_palindromic,
        "n_drop_missing_in_AF": r.n_drop_missing_in_AF,
        "matched_rate": (r.n_matched_AF / r.n_clumped) if r.n_clumped else None,
    } for r in results]).write_csv(OUT_MATCH, separator="\t")
    pl.DataFrame([{
        "aptamer": r.aptamer, "gene": r.gene, "error": r.error, "retries": r.retries,
    } for r in results if r.error]).write_csv(OUT_FAILED, separator="\t")


def _checkpoint_metrics(results: list[WorkerResult],
                        wall_start: float,
                        n_total: int) -> dict:
    completed = [r for r in results if r.actual_size > 0]
    fawithd = [r for r in results if r.error]
    n_done = len(results)
    n_completed = len(completed)
    n_empty = sum(1 for r in completed if r.n_clumped == 0)
    total_bytes = sum(r.actual_size for r in completed)
    secs = [r.download_seconds + r.extract_seconds for r in completed if r.download_seconds > 0]
    median_sec = median(secs) if secs else 0.0
    median_bytes = median([r.actual_size for r in completed]) if completed else 0
    n_clumped = sum(r.n_clumped for r in completed)
    n_matched = sum(r.n_matched_AF for r in completed)
    elapsed_h = (time.time() - wall_start) / 3600

    # Project remaining at 4-worker assumed 2.5x speedup
    remaining = max(0, n_total - n_done)
    proj_remaining_sec = remaining * median_sec / 2.5
    proj_total_h = elapsed_h + proj_remaining_sec / 3600
    proj_transfer_gb = (total_bytes + remaining * median_bytes) / 1e9
    empty_rate = (n_empty / n_completed) if n_completed else 0.0
    matched_rate = (n_matched / n_clumped) if n_clumped else 0.0
    fail_rate = len(fawithd) / max(1, n_done)
    disk_free = _disk_free_gb()

    return {
        "n_done": n_done, "n_completed": n_completed,
        "n_fawithd": len(fawithd), "n_empty": n_empty,
        "elapsed_hours": round(elapsed_h, 2),
        "median_sec_per_aptamer": round(median_sec, 1),
        "median_bytes_per_aptamer": int(median_bytes),
        "transfer_gb_observed": round(total_bytes / 1e9, 1),
        "projected_total_hours": round(proj_total_h, 2),
        "projected_total_transfer_gb": round(proj_transfer_gb, 1),
        "empty_rate": round(empty_rate, 3),
        "matched_to_AF_rate": round(matched_rate, 3),
        "fawithd_rate": round(fail_rate, 3),
        "disk_free_gb": round(disk_free, 1),
    }


def _check_dynamic_stop(metrics: dict) -> tuple[bool, list[str]]:
    failures = []
    if metrics["projected_total_transfer_gb"] > DS_TRANSFER_GB_MAX:
        failures.append(f"transfer {metrics['projected_total_transfer_gb']} > {DS_TRANSFER_GB_MAX}")
    if metrics["projected_total_hours"] > DS_RUNTIME_HOURS_MAX:
        failures.append(f"runtime {metrics['projected_total_hours']} > {DS_RUNTIME_HOURS_MAX}")
    if metrics["empty_rate"] > DS_EMPTY_RATE_MAX:
        failures.append(f"empty {metrics['empty_rate']} > {DS_EMPTY_RATE_MAX}")
    if metrics["matched_to_AF_rate"] < DS_MATCH_RATE_MIN and metrics["n_completed"] >= 30:
        failures.append(f"match {metrics['matched_to_AF_rate']} < {DS_MATCH_RATE_MIN}")
    if metrics["fawithd_rate"] > DS_FAIL_RATE_MAX and metrics["n_done"] >= 30:
        failures.append(f"fawithd {metrics['fawithd_rate']} > {DS_FAIL_RATE_MAX}")
    if metrics["disk_free_gb"] < DS_DISK_FREE_GB_MIN:
        failures.append(f"disk_free {metrics['disk_free_gb']} < {DS_DISK_FREE_GB_MIN}")
    return (len(failures) == 0), failures


def main() -> int:
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    targets = pl.read_csv(DESIGN_TSV, separator="\t").to_dicts()
    _log(f"loaded {len(targets)} CVD454 targets")

    coord_cache = json.loads(COORD_CACHE.read_text()) if COORD_CACHE.exists() else {}
    # Ensure we have coords for all CVD454 targets — fill missing via Ensembl
    missing = [t["gene"] for t in targets if t["gene"] not in coord_cache]
    if missing:
        _log(f"resolving {len(missing)} gene coords via Ensembl REST")
        for gene in missing:
            try:
                r = requests.get(
                    f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}",
                    headers={"Content-Type": "application/json"},
                    timeout=15,
                )
                if r.status_code == 200:
                    d = r.json()
                    chrom = str(d.get("seq_region_name", ""))
                    start = d.get("start"); end = d.get("end")
                    if chrom and start and end:
                        coord_cache[gene] = [chrom, int(start), int(end)]
                        continue
            except Exception:
                pass
            coord_cache[gene] = None
        COORD_CACHE.write_text(json.dumps(coord_cache, indent=2))
        _log(f"coord cache size: {len(coord_cache)}")

    # AF parquet (loaded once, shared read-only across workers via polars join).
    af = pl.read_parquet(AF_PARQUET).select([
        "chr", "pos", "rsid", "effect_allele", "other_allele",
        "beta", "se", "pval", "eaf",
    ])
    _log(f"AF parquet loaded ({af.height:,} rows); per-aptamer polars join for matching")

    wall_start = time.time()
    results: list[WorkerResult] = []
    results_lock = threading.Lock()
    checkpoint_rows = []

    # Convert coord_cache values to plain dict tuples for thread safety
    coord_lookup = {g: tuple(v) if isinstance(v, list) else v
                    for g, v in coord_cache.items()}

    def submit(executor):
        return [executor.submit(_process_aptamer, t, coord_lookup, af)
                for t in targets]

    stop_signal = {"value": False, "reason": []}

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = submit(ex)
        for fut in as_completed(futures):
            r = fut.result()
            with results_lock:
                results.append(r)
                _log(f"[{len(results)}/{len(targets)}] {r.gene:8} {r.aptamer[:35]:35} "
                     f"clumped={r.n_clumped:2} matched={r.n_matched_AF:2} "
                     f"err={r.error[:30]}")
                # Checkpoint?
                if (len(results) % CHECKPOINT_EVERY == 0) or (len(results) == len(targets)):
                    metrics = _checkpoint_metrics(results, wall_start, len(targets))
                    metrics["checkpoint_n"] = len(results)
                    metrics["timestamp"] = time.strftime("%Y-%m-%d %H:%M:%S")
                    checkpoint_rows.append(metrics)
                    pl.DataFrame(checkpoint_rows).write_csv(OUT_CHECKPOINTS, separator="\t")
                    _persist_partial(results, when=f"checkpoint_{len(results)}")
                    _log(f"  checkpoint @ {len(results)}: empty={metrics['empty_rate']:.0%} "
                         f"match={metrics['matched_to_AF_rate']:.0%} "
                         f"proj_total_h={metrics['projected_total_hours']:.1f} "
                         f"proj_transfer_gb={metrics['projected_total_transfer_gb']:.0f} "
                         f"disk_free={metrics['disk_free_gb']:.0f}gb")

                    # Dynamic stop check at 50, 100, 250
                    if len(results) in DYNAMIC_STOP_CHECKPOINTS:
                        ok, fails = _check_dynamic_stop(metrics)
                        if not ok:
                            stop_signal["value"] = True
                            stop_signal["reason"] = fails
                            _log(f"DYNAMIC STOP triggered @ {len(results)}: {fails}")
                            for f in futures:
                                f.cancel()
                            break

    elapsed_h = (time.time() - wall_start) / 3600
    _log(f"extraction wall-time: {elapsed_h:.2f} h")

    _persist_partial(results, when="final")

    # Done — also clean tmp dir
    if TMP_DIR.exists():
        for f in TMP_DIR.iterdir():
            try:
                f.unlink()
            except OSError:
                pass

    if stop_signal["value"]:
        _log(f"DYNAMIC STOP — extraction halted: {stop_signal['reason']}")
        return 2

    # Final summary
    n_completed = sum(1 for r in results if r.actual_size > 0)
    n_clumped = sum(r.n_clumped for r in results)
    n_matched = sum(r.n_matched_AF for r in results)
    _log(f"CVD454 extraction COMPLETE: {n_completed}/{len(targets)} aptamers, "
         f"{n_clumped} cis instruments, {n_matched} matched to AF")
    return 0


if __name__ == "__main__":
    sys.exit(main())
