#!/usr/bin/env python3
"""D038 + D039 — CVD150 Wave 1 selective coloc with sanity validation.

D039 modes:
  --preflight   : run only IL6R × AF and NPPA × CES, emit preflight_validation.tsv
                  Numerical sanity:
                    PP sums to 1.0 ± 1e-6
                    PP.H3 not negative / NaN
                    n_variants reasonable (>= 500 in ±500 kb region)
                    sign-flip allele logic visible in counters
                    no underflow in log-ABF
                    no p-filter inside region
                  Optional: R coloc.abf cross-check (skip with note if R missing)
  --run-wave1   : run all 15 D038 pairs

Priors (D039):
  primary  prior_model=default_quant_cc : p1=p2=1e-4 p12=1e-5 W_pqtl=0.0225 W_outcome=0.04
  sensitivity prior_model=flat_W04      : p1=p2=1e-4 p12=5e-6 W_pqtl=0.04   W_outcome=0.04
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import polars as pl
import requests

ROOT = Path("/Users/apple/Desktop/gwas_af")
INSTR_PARQUET = ROOT / "data/processed/instruments/decode_cvd150_instruments.parquet"
COORD_CACHE = ROOT / "data/registry/gene_coords_b38_cache.json"
CHAIN_PATH = ROOT / "data/raw/annotations/hg38ToHg19.over.chain.gz"
TMP_DIR = Path("/tmp/decode_wave1_coloc")

OUT_DIR = ROOT / "results/coloc"
OUT_RES = OUT_DIR / "decode_cvd150_wave1_coloc_results.tsv"
OUT_REG = OUT_DIR / "decode_cvd150_wave1_region_qc.tsv"
OUT_ART = OUT_DIR / "decode_cvd150_wave1_artifact_flags.tsv"
OUT_MD = OUT_DIR / "decode_cvd150_wave1_coloc_summary.md"
OUT_PREFLIGHT = OUT_DIR / "decode_cvd150_wave1_preflight_validation.tsv"

TOKEN = os.environ.get("DECODE_TOKEN", "<DECODE_TOKEN>")
DOWNLOAD_API = "https://download.decode.is/s3/download"
REGION_KB_PRIMARY = 500
PALINDROMIC = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
PALINDROMIC_AMBIG_WINDOW = 0.05

PRIOR_PRIMARY = {
    "name": "default_quant_cc",
    "p1": 1.0e-4, "p2": 1.0e-4, "p12": 1.0e-5,
    "W_pqtl": 0.0225, "W_outcome": 0.04,
}
PRIOR_SENSITIVITY = {
    "name": "flat_W04",
    "p1": 1.0e-4, "p2": 1.0e-4, "p12": 5.0e-6,
    "W_pqtl": 0.04, "W_outcome": 0.04,
}

OUTCOMES = {
    "AF":                   {"path": ROOT / "data/processed/sumstats/AF.parquet",                       "build": "GRCh38"},
    "HF_overall":           {"path": ROOT / "data/processed/sumstats/HF_2024_EUR_HF_overall.parquet",   "build": "GRCh37"},
    "HF_ni_HFpEF":          {"path": ROOT / "data/processed/sumstats/HF_2024_EUR_HF_ni_HFpEF.parquet", "build": "GRCh37"},
    "ischemic_stroke":      {"path": ROOT / "data/processed/sumstats/STROKE_AIS_EUR.parquet",          "build": "GRCh38"},
    "cardioembolic_stroke": {"path": ROOT / "data/processed/sumstats/STROKE_CES_EUR.parquet",          "build": "GRCh38"},
    "large_artery_stroke":  {"path": ROOT / "data/processed/sumstats/STROKE_LAS_EUR.parquet",          "build": "GRCh38"},
    "small_vessel_stroke":  {"path": ROOT / "data/processed/sumstats/STROKE_SVS_EUR.parquet",          "build": "GRCh38"},
}

WAVE1_PAIRS: list[tuple[str, str, str]] = [
    ("15602_43_IL6R_IL_6_sRa.txt.gz",  "IL6R",  "AF"),
    ("15602_43_IL6R_IL_6_sRa.txt.gz",  "IL6R",  "HF_overall"),
    ("15602_43_IL6R_IL_6_sRa.txt.gz",  "IL6R",  "ischemic_stroke"),
    ("15602_43_IL6R_IL_6_sRa.txt.gz",  "IL6R",  "cardioembolic_stroke"),
    ("13126_52_DSC2_DSC2.txt.gz",      "DSC2",  "AF"),
    ("13126_52_DSC2_DSC2.txt.gz",      "DSC2",  "small_vessel_stroke"),
    ("13126_52_DSC2_DSC2.txt.gz",      "DSC2",  "HF_ni_HFpEF"),
    ("17224_12_OGN_MIME.txt.gz",       "OGN",   "AF"),
    ("17224_12_OGN_MIME.txt.gz",       "OGN",   "ischemic_stroke"),
    ("17224_12_OGN_MIME.txt.gz",       "OGN",   "cardioembolic_stroke"),
    ("17224_12_OGN_MIME.txt.gz",       "OGN",   "large_artery_stroke"),
    ("5443_62_NPPA_ANP.txt.gz",        "NPPA",  "AF"),
    ("5443_62_NPPA_ANP.txt.gz",        "NPPA",  "cardioembolic_stroke"),
    ("2480_58_TIMP3_TIMP_3.txt.gz",    "TIMP3", "AF"),
    ("2480_58_TIMP3_TIMP_3.txt.gz",    "TIMP3", "cardioembolic_stroke"),
]

PREFLIGHT_PAIRS = [
    ("15602_43_IL6R_IL_6_sRa.txt.gz", "IL6R", "AF"),
    ("5443_62_NPPA_ANP.txt.gz",       "NPPA", "cardioembolic_stroke"),
]


def _log(msg: str) -> None:
    print(f"[wave1 {time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Region pQTL extraction
# ---------------------------------------------------------------------------


def download_aptamer(key: str, dest: Path) -> int:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1_000_000:
        return dest.stat().st_size
    _log(f"  downloading {key}")
    t0 = time.time()
    with requests.get(DOWNLOAD_API,
                      params={"token": TOKEN, "file": key},
                      stream=True, timeout=600) as r:
        r.raise_for_status()
        with dest.open("wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
    sz = dest.stat().st_size
    _log(f"  done {sz/1e6:.0f} MB in {time.time()-t0:.0f}s")
    return sz


def extract_region(gz_path: Path, target_chr: str,
                   start: int, end: int) -> list[dict]:
    """Extract ALL variants in chr:start-end region. NO p-value filter (D025)."""
    rows: list[dict] = []
    target_chr_norm = str(target_chr).replace("chr", "")
    with gzip.open(gz_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        ci, pi = idx["Chrom"], idx["Pos"]
        ni, ri = idx["Name"], idx.get("rsids", -1)
        ei, oi = idx["effectAllele"], idx["otherAllele"]
        bi, pvi = idx["Beta"], idx["Pval"]
        si = idx["SE"]
        mi = idx.get("ImpMAF", -1)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if cols[ci].replace("chr", "") != target_chr_norm:
                continue
            try:
                pos = int(cols[pi])
            except ValueError:
                continue
            if pos < start or pos > end:
                continue
            try:
                beta = float(cols[bi])
                se = float(cols[si])
                pval = float(cols[pvi])
            except ValueError:
                continue
            if se <= 0:
                continue
            ea = cols[ei].upper()
            oa = cols[oi].upper()
            imp_maf = float(cols[mi]) if mi >= 0 and cols[mi] not in ("", "NA") else float("nan")
            rows.append({
                "chrom": target_chr_norm, "pos": pos,
                "name": cols[ni], "rsid": cols[ri] if ri >= 0 else "",
                "ea": ea, "oa": oa, "beta": beta, "se": se,
                "pval": pval, "imp_maf": imp_maf,
            })
    return rows


def liftover_b37_for_region(pqtl_rows: list[dict]) -> list[dict]:
    from pyliftover import LiftOver
    lo = LiftOver(str(CHAIN_PATH))
    out = []
    for r in pqtl_rows:
        cv = lo.convert_coordinate("chr" + str(r["chrom"]), int(r["pos"]))
        if cv:
            r2 = dict(r)
            r2["pos_b37"] = int(cv[0][1])
            out.append(r2)
    return out


def harmonize_region(pqtl_rows: list[dict],
                     outcome_df: pl.DataFrame,
                     build: str) -> tuple[list[dict], dict]:
    pos_key = "pos_b37" if build == "GRCh37" else "pos"
    counters = {
        "n_pqtl_in_region": len(pqtl_rows),
        "n_outcome_matched": 0,
        "n_dropped_allele_mismatch": 0,
        "n_dropped_palindromic": 0,
        "n_dropped_missing_outcome": 0,
        "n_sign_flipped": 0,
    }
    if not pqtl_rows:
        return [], counters

    keys = pl.DataFrame({
        "chr": [str(r["chrom"]) for r in pqtl_rows],
        "pos": [int(r[pos_key]) for r in pqtl_rows],
    })
    outcome_subset = outcome_df.join(keys, on=["chr", "pos"], how="inner")
    by_key = {(str(r["chr"]), int(r["pos"])): r
              for r in outcome_subset.iter_rows(named=True)}

    aligned = []
    for r in pqtl_rows:
        oc = by_key.get((str(r["chrom"]), int(r[pos_key])))
        if oc is None:
            counters["n_dropped_missing_outcome"] += 1
            continue
        ea_e, oa_e = r["ea"], r["oa"]
        ea_o = (oc["effect_allele"] or "").upper()
        oa_o = (oc["other_allele"] or "").upper()
        if {ea_e, oa_e} != {ea_o, oa_o}:
            counters["n_dropped_allele_mismatch"] += 1
            continue
        sign_flip = (ea_e == oa_o and oa_e == ea_o)
        if (ea_e, oa_e) in PALINDROMIC:
            eaf_o = oc.get("eaf")
            if eaf_o is None or eaf_o != eaf_o or abs(float(eaf_o) - 0.5) < PALINDROMIC_AMBIG_WINDOW:
                counters["n_dropped_palindromic"] += 1
                continue
        beta_o = float(oc["beta"]) * (-1.0 if sign_flip else 1.0)
        if sign_flip:
            counters["n_sign_flipped"] += 1
        aligned.append({
            **r,
            "beta_outcome": beta_o,
            "se_outcome": float(oc["se"]),
            "pval_outcome": float(oc["pval"]),
            "eaf_outcome": oc.get("eaf"),
        })
    counters["n_outcome_matched"] = len(aligned)
    return aligned, counters


# ---------------------------------------------------------------------------
# coloc.abf (Giambartolomei 2014, Wakefield 2009 ABF)
# ---------------------------------------------------------------------------


def _logsumexp(xs: list[float]) -> float:
    finite = [x for x in xs if x != float("-inf") and not math.isnan(x)]
    if not finite:
        return float("-inf")
    m = max(finite)
    return m + math.log(sum(math.exp(x - m) for x in finite))


def log_abf(z: float, V: float, W: float) -> float:
    """Wakefield 2009 ABF, log scale.  log ABF = 0.5 * log(1-r) + 0.5 * z² * r ; r = W/(V+W)."""
    if V <= 0:
        return float("-inf")
    r = W / (V + W)
    if r <= 0 or r >= 1:
        return float("-inf")
    return 0.5 * math.log(1.0 - r) + 0.5 * (z * z) * r


def coloc_abf(rows: list[dict], prior: dict) -> dict:
    n = len(rows)
    if n == 0:
        return {"PP.H0": float("nan"), "PP.H1": float("nan"),
                "PP.H2": float("nan"), "PP.H3": float("nan"),
                "PP.H4": float("nan"), "n_variants": 0,
                "lead_pqtl_pos": None, "lead_pqtl_rsid": "",
                "lead_pqtl_pval": None, "lead_outcome_pos": None,
                "lead_outcome_rsid": "", "lead_outcome_pval": None,
                "log_L_xy": None, "log_L_distinct": None,
                "h3_negative_flag": False, "underflow_flag": False,
                "underflow_count": 0}

    W_x = prior["W_pqtl"]; W_y = prior["W_outcome"]
    p1 = prior["p1"]; p2 = prior["p2"]; p12 = prior["p12"]

    log_abf_x, log_abf_y = [], []
    underflow_count = 0
    for r in rows:
        beta_x = r["beta"]; se_x = r["se"]
        beta_y = r["beta_outcome"]; se_y = r["se_outcome"]
        z_x = beta_x / se_x if se_x > 0 else 0.0
        z_y = beta_y / se_y if se_y > 0 else 0.0
        V_x = se_x * se_x; V_y = se_y * se_y
        ax = log_abf(z_x, V_x, W_x)
        ay = log_abf(z_y, V_y, W_y)
        if ax == float("-inf") or ay == float("-inf"):
            underflow_count += 1
        log_abf_x.append(ax)
        log_abf_y.append(ay)

    L_x = _logsumexp(log_abf_x)
    L_y = _logsumexp(log_abf_y)
    L_xy = _logsumexp([log_abf_x[i] + log_abf_y[i] for i in range(n)])

    log_sum_x_sum_y = L_x + L_y
    h3_negative = False
    if log_sum_x_sum_y > L_xy:
        diff = math.exp(L_xy - log_sum_x_sum_y)
        if diff < 1.0:
            L_distinct = log_sum_x_sum_y + math.log(1.0 - diff)
        else:
            L_distinct = float("-inf")
            h3_negative = True
    else:
        L_distinct = float("-inf")
        h3_negative = True

    log_p1 = math.log(p1); log_p2 = math.log(p2); log_p12 = math.log(p12)
    log_terms = [
        0.0,                                  # H0 (prior_H0 normalized in 5-hyp posterior)
        log_p1 + L_x,                        # H1
        log_p2 + L_y,                        # H2
        log_p1 + log_p2 + L_distinct,        # H3
        log_p12 + L_xy,                      # H4
    ]
    log_total = _logsumexp(log_terms)
    pp = [math.exp(x - log_total) for x in log_terms]

    lead_pqtl_idx = min(range(n), key=lambda i: rows[i]["pval"])
    lead_outcome_idx = min(range(n), key=lambda i: rows[i]["pval_outcome"])

    return {
        "PP.H0": pp[0], "PP.H1": pp[1], "PP.H2": pp[2],
        "PP.H3": pp[3], "PP.H4": pp[4],
        "n_variants": n,
        "lead_pqtl_pos": rows[lead_pqtl_idx]["pos"],
        "lead_pqtl_rsid": rows[lead_pqtl_idx]["rsid"],
        "lead_pqtl_pval": rows[lead_pqtl_idx]["pval"],
        "lead_outcome_pos": rows[lead_outcome_idx]["pos"],
        "lead_outcome_rsid": rows[lead_outcome_idx]["rsid"],
        "lead_outcome_pval": rows[lead_outcome_idx]["pval_outcome"],
        "log_L_xy": L_xy, "log_L_distinct": L_distinct,
        "h3_negative_flag": h3_negative,
        "underflow_flag": underflow_count > 0,
        "underflow_count": underflow_count,
    }


# ---------------------------------------------------------------------------
# Optional R coloc cross-check
# ---------------------------------------------------------------------------


def try_r_coloc(rows: list[dict], prior: dict) -> dict | None:
    rscript = shutil.which("Rscript")
    if not rscript:
        return None
    if not rows:
        return None
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        in_csv = td / "in.csv"
        with in_csv.open("w") as f:
            f.write("snp,beta_x,se_x,beta_y,se_y\n")
            for i, r in enumerate(rows):
                f.write(f"v{i},{r['beta']},{r['se']},{r['beta_outcome']},{r['se_outcome']}\n")
        rcode = f"""
suppressMessages({{
  ok <- requireNamespace("coloc", quietly = TRUE)
  if (!ok) {{
    cat('{{"r_coloc_available": false}}')
    quit(status = 0)
  }}
  library(coloc)
  d <- read.csv("{in_csv}")
  ds1 <- list(beta = d$beta_x, varbeta = d$se_x^2, type = "quant",
              sdY = sqrt({prior['W_pqtl']}), snp = d$snp,
              N = rep(35559, nrow(d)))
  ds2 <- list(beta = d$beta_y, varbeta = d$se_y^2, type = "cc",
              s = 0.1, snp = d$snp, N = rep(100000, nrow(d)))
  res <- coloc.abf(ds1, ds2, p1 = {prior['p1']}, p2 = {prior['p2']}, p12 = {prior['p12']})
  pp <- res$summary
  cat(sprintf('{{"r_coloc_available": true, "PP.H0": %.6f, "PP.H1": %.6f, "PP.H2": %.6f, "PP.H3": %.6f, "PP.H4": %.6f}}',
              pp["PP.H0.abf"], pp["PP.H1.abf"], pp["PP.H2.abf"], pp["PP.H3.abf"], pp["PP.H4.abf"]))
}})
"""
        try:
            result = subprocess.run(
                [rscript, "-e", rcode],
                capture_output=True, text=True, timeout=180,
            )
            log_path = Path("/tmp/decode_wave1_coloc/r_coloc_debug.log")
            with log_path.open("a") as logf:
                logf.write(f"\n=== R run rc={result.returncode}\n")
                logf.write(f"STDOUT: {result.stdout!r}\n")
                logf.write(f"STDERR: {result.stderr[:1500]!r}\n")
            if result.returncode != 0:
                return {"r_coloc_available": False, "error": result.stderr[:200]}
            try:
                return json.loads(result.stdout)
            except json.JSONDecodeError as je:
                return {"r_coloc_available": False, "error": f"JSON parse: {je} | stdout={result.stdout[:200]!r}"}
        except subprocess.TimeoutExpired as e:
            return {"r_coloc_available": False, "error": str(e)[:200]}


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def init_environment(target_aptamers: set[str]):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    instr = pl.read_parquet(INSTR_PARQUET)
    coord_cache = json.loads(COORD_CACHE.read_text()) if COORD_CACHE.exists() else {}
    _log(f"loaded {instr.height} CVD150 instruments")

    leads, apt_paths = {}, {}
    for apt in sorted(target_aptamers):
        sub = instr.filter(pl.col("aptamer") == apt).sort("pval")
        if sub.is_empty():
            _log(f"WARN: no instruments for {apt}")
            continue
        row = sub.row(0, named=True)
        leads[apt] = {
            "chrom": str(row["chrom"]), "pos": int(row["pos"]),
            "rsid": row.get("rsid", ""), "gene": row["gene"],
            "ea": row["effect_allele"], "oa": row["other_allele"],
            "beta": float(row["beta"]), "se": float(row["se"]), "pval": float(row["pval"]),
        }
        _log(f"  lead {apt} {row['gene']} chr{row['chrom']}:{row['pos']} p={row['pval']:.2e}")
        dest = TMP_DIR / apt
        try:
            download_aptamer(apt, dest)
            apt_paths[apt] = dest
        except Exception as e:
            _log(f"FAIL download {apt}: {e}")
    return leads, apt_paths, coord_cache, instr, {"region": {}, "region_b37": {}, "outcome": {}}


def build_pair_rows(apt, gene, outcome_label, *, pqtl_region_cache,
                    pqtl_region_b37_cache, apt_paths, leads, outcome_cache):
    if apt not in leads or apt not in apt_paths:
        return None, None
    meta = OUTCOMES[outcome_label]
    if outcome_label not in outcome_cache:
        outcome_cache[outcome_label] = pl.read_parquet(meta["path"]).select([
            "chr", "pos", "rsid", "effect_allele", "other_allele",
            "beta", "se", "pval", "eaf",
        ])
    if apt not in pqtl_region_cache:
        lead = leads[apt]
        start = max(0, lead["pos"] - REGION_KB_PRIMARY * 1000)
        end = lead["pos"] + REGION_KB_PRIMARY * 1000
        pqtl_region_cache[apt] = extract_region(apt_paths[apt], lead["chrom"], start, end)
        _log(f"  region {apt} {gene}: {len(pqtl_region_cache[apt])} variants in ±{REGION_KB_PRIMARY}kb")
    pqtl_rows = pqtl_region_cache[apt]
    if meta["build"] == "GRCh37":
        if apt not in pqtl_region_b37_cache:
            pqtl_region_b37_cache[apt] = liftover_b37_for_region(pqtl_rows)
            _log(f"  liftover b38→b37 {apt}: {len(pqtl_region_b37_cache[apt])}/{len(pqtl_rows)}")
        pqtl_for_outcome = pqtl_region_b37_cache[apt]
    else:
        pqtl_for_outcome = pqtl_rows
    aligned, counters = harmonize_region(pqtl_for_outcome, outcome_cache[outcome_label], meta["build"])
    return aligned, counters


def annotate_artifact(gene, lead_pos, coord_cache, instr) -> dict:
    intragenic = False
    if lead_pos is not None:
        gc = coord_cache.get(gene)
        if gc:
            try:
                gs = int(gc[1]); ge = int(gc[2])
                intragenic = gs <= int(lead_pos) <= ge
            except (TypeError, ValueError):
                intragenic = False
    n_aptamers = instr.filter(pl.col("gene") == gene)["aptamer"].n_unique()
    multi = "single" if n_aptamers == 1 else f"multi_n={n_aptamers}"
    return {"intragenic": intragenic, "multi_aptamer_status": multi}


def run_pair(apt, gene, outcome_label, prior, *, pqtl_region_cache,
             pqtl_region_b37_cache, apt_paths, leads, outcome_cache,
             coord_cache, instr):
    aligned, counters = build_pair_rows(
        apt, gene, outcome_label,
        pqtl_region_cache=pqtl_region_cache,
        pqtl_region_b37_cache=pqtl_region_b37_cache,
        apt_paths=apt_paths, leads=leads, outcome_cache=outcome_cache,
    )
    if aligned is None:
        empty = coloc_abf([], prior)
        return ({"gene": gene, "aptamer": apt, "outcome": outcome_label,
                 "outcome_build": OUTCOMES[outcome_label]["build"], **empty,
                 "prior_model": prior["name"], "p1": prior["p1"], "p2": prior["p2"],
                 "p12": prior["p12"], "W_pqtl": prior["W_pqtl"],
                 "W_outcome": prior["W_outcome"],
                 "lead_pqtl_intragenic": False, "multi_aptamer_status": "missing",
                 "scope": "CVD150_Wave1_NOT_proteome_wide"},
                {"gene": gene, "aptamer": apt, "outcome": outcome_label,
                 "n_pqtl_in_region": 0, "n_outcome_matched": 0,
                 "n_dropped_allele_mismatch": 0, "n_dropped_palindromic": 0,
                 "n_dropped_missing_outcome": 0, "n_sign_flipped": 0},
                {"intragenic": False, "multi_aptamer_status": "missing"},
                aligned)
    coloc = coloc_abf(aligned, prior)
    art = annotate_artifact(gene, coloc.get("lead_pqtl_pos"), coord_cache, instr)
    res_row = {
        "gene": gene, "aptamer": apt, "outcome": outcome_label,
        "outcome_build": OUTCOMES[outcome_label]["build"],
        **coloc,
        "prior_model": prior["name"], "p1": prior["p1"], "p2": prior["p2"],
        "p12": prior["p12"], "W_pqtl": prior["W_pqtl"],
        "W_outcome": prior["W_outcome"],
        "lead_pqtl_intragenic": art["intragenic"],
        "multi_aptamer_status": art["multi_aptamer_status"],
        "scope": "CVD150_Wave1_NOT_proteome_wide",
    }
    qc_row = {"gene": gene, "aptamer": apt, "outcome": outcome_label, **counters}
    return res_row, qc_row, art, aligned


def cleanup_tmp():
    for f in TMP_DIR.iterdir():
        try:
            f.unlink()
        except OSError:
            pass


def write_results(coloc_rows, qc_rows, artifact_rows):
    pl.DataFrame(coloc_rows).write_csv(OUT_RES, separator="\t")
    pl.DataFrame(qc_rows).write_csv(OUT_REG, separator="\t")
    pl.DataFrame(artifact_rows).write_csv(OUT_ART, separator="\t")
    _log(f"wrote {OUT_RES} ({len(coloc_rows)} rows)")


def write_summary_md(coloc_rows):
    primary = [r for r in coloc_rows if r["prior_model"] == PRIOR_PRIMARY["name"]]
    sens = {(r["gene"], r["outcome"]): r for r in coloc_rows
            if r["prior_model"] == PRIOR_SENSITIVITY["name"]}
    strong = [r for r in primary if r["PP.H4"] is not None and not math.isnan(r["PP.H4"]) and r["PP.H4"] >= 0.80]
    moderate = [r for r in primary if r["PP.H4"] is not None and not math.isnan(r["PP.H4"]) and 0.50 <= r["PP.H4"] < 0.80]
    weak = [r for r in primary if r["PP.H4"] is not None and not math.isnan(r["PP.H4"]) and r["PP.H4"] < 0.50]
    ld_warn = [r for r in primary if r["PP.H3"] is not None and r["PP.H4"] is not None
               and not math.isnan(r["PP.H3"]) and not math.isnan(r["PP.H4"])
               and r["PP.H3"] >= 0.50 and r["PP.H3"] > r["PP.H4"]]

    md = [
        "# CVD150 Wave 1 colocalization (D038/D039, NOT proteome-wide)",
        "",
        f"**Scope:** {len(primary)} target-outcome pairs, primary + sensitivity prior models.",
        f"**Region:** lead cis-pQTL ±{REGION_KB_PRIMARY} kb.",
        f"**Primary prior** ({PRIOR_PRIMARY['name']}): p1=p2={PRIOR_PRIMARY['p1']}, p12={PRIOR_PRIMARY['p12']}, "
        f"W_pqtl={PRIOR_PRIMARY['W_pqtl']}, W_outcome={PRIOR_PRIMARY['W_outcome']}",
        f"**sensitivity** ({PRIOR_SENSITIVITY['name']}): p12={PRIOR_SENSITIVITY['p12']}, "
        f"W=0.04 (both)",
        "",
        "## Per-pair results (primary prior; sensitivity PP.H4 in parentheses)",
        "",
        "| Gene | Outcome | n_present | PP.H0 | PP.H1 | PP.H2 | PP.H3 | PP.H4 (sens) | Lead pQTL | Lead outcome p | Intragenic | Multi-apt |",
        "|---|---|---:|---:|---:|---:|---:|---|---|---|:-:|:-:|",
    ]
    for r in primary:
        s = sens.get((r["gene"], r["outcome"]))
        sens_h4 = f" ({s['PP.H4']:.3f})" if s and s["PP.H4"] is not None and not math.isnan(s["PP.H4"]) else ""
        h4_str = f"**{r['PP.H4']:.3f}**{sens_h4}" if r["PP.H4"] is not None and not math.isnan(r["PP.H4"]) else "—"
        md.append(
            f"| {r['gene']} | {r['outcome']} | {r['n_variants']} | "
            f"{r['PP.H0']:.3f} | {r['PP.H1']:.3f} | {r['PP.H2']:.3f} | "
            f"{r['PP.H3']:.3f} | {h4_str} | "
            f"{r['lead_pqtl_rsid']} (p={r['lead_pqtl_pval']:.1e}) | "
            f"{r['lead_outcome_pval']:.1e} | "
            f"{'✓' if r['lead_pqtl_intragenic'] else '✗'} | "
            f"{r['multi_aptamer_status']} |"
        )
    md.extend([
        "",
        "## Summary by PP.H4 strength (primary prior)",
        "",
        f"- **Strong (PP.H4 ≥ 0.80):** {len(strong)} pairs",
        f"- **Moderate (0.50 ≤ PP.H4 < 0.80):** {len(moderate)} pairs",
        f"- **Weak (PP.H4 < 0.50):** {len(weak)} pairs",
        f"- **LD-confounding warning (PP.H3 ≥ 0.50 and PP.H3 > PP.H4):** {len(ld_warn)} pairs",
        "",
        "## Operasyonel kural",
        "",
        "Final classification (GREEN/RED/AMBER/BLUE/PURPLE) **STILL LOCKED**.",
        "Wave 2 recommendations: KNG1 × AIS/CES, F11 × AIS/CES, PCSK9 × HFrEF, F13B + Wave 1 subtype expansions.",
    ])
    OUT_MD.write_text("\n".join(md))
    _log(f"wrote {OUT_MD}")


def cmd_preflight():
    target_apts = {p[0] for p in PREFLIGHT_PAIRS}
    leads, apt_paths, coord_cache, instr, caches = init_environment(target_apts)
    audit_dir = OUT_DIR / "preflight_audit"
    audit_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for apt, gene, outcome in PREFLIGHT_PAIRS:
        res, qc, _art, aligned = run_pair(
            apt, gene, outcome, PRIOR_PRIMARY,
            pqtl_region_cache=caches["region"],
            pqtl_region_b37_cache=caches["region_b37"],
            apt_paths=apt_paths, leads=leads,
            outcome_cache=caches["outcome"],
            coord_cache=coord_cache, instr=instr,
        )
        if aligned:
            audit_path = audit_dir / f"aligned_{gene}_{outcome}.tsv"
            pl.DataFrame(aligned).write_csv(audit_path, separator="\t")
            _log(f"  audit dump: {audit_path} ({len(aligned)} rows)")
        r_check = try_r_coloc(aligned or [], PRIOR_PRIMARY) if aligned else None
        if r_check and not r_check.get("r_coloc_available"):
            _log(f"  R coloc error: {r_check.get('error', '?')}")
        diff_max = None
        r_status = "unavailable"
        if r_check and r_check.get("r_coloc_available"):
            r_status = "available"
            diff_max = max(
                abs(r_check["PP.H0"] - res["PP.H0"]),
                abs(r_check["PP.H1"] - res["PP.H1"]),
                abs(r_check["PP.H2"] - res["PP.H2"]),
                abs(r_check["PP.H3"] - res["PP.H3"]),
                abs(r_check["PP.H4"] - res["PP.H4"]),
            )
        rows.append({
            "gene": gene, "outcome": outcome,
            "n_variants": res["n_variants"],
            "PP.H0": res["PP.H0"], "PP.H1": res["PP.H1"],
            "PP.H2": res["PP.H2"], "PP.H3": res["PP.H3"], "PP.H4": res["PP.H4"],
            "h3_negative_flag": res.get("h3_negative_flag"),
            "underflow_flag": res.get("underflow_flag"),
            "underflow_count": res.get("underflow_count"),
            "n_sign_flipped": qc.get("n_sign_flipped"),
            "n_dropped_allele_mismatch": qc.get("n_dropped_allele_mismatch"),
            "n_dropped_palindromic": qc.get("n_dropped_palindromic"),
            "n_dropped_missing_outcome": qc.get("n_dropped_missing_outcome"),
            "r_coloc_status": r_status,
            "r_python_max_abs_diff": diff_max,
            "r_diff_passes_1e3": (diff_max is not None and diff_max < 1e-3) if r_status == "available" else None,
            "p1": PRIOR_PRIMARY["p1"], "p2": PRIOR_PRIMARY["p2"],
            "p12": PRIOR_PRIMARY["p12"],
            "W_pqtl": PRIOR_PRIMARY["W_pqtl"], "W_outcome": PRIOR_PRIMARY["W_outcome"],
        })

    # Validation
    failures = []
    for r in rows:
        tag = f"{r['gene']}-{r['outcome']}"
        pp_sum = (r["PP.H0"] or 0) + (r["PP.H1"] or 0) + (r["PP.H2"] or 0) + (r["PP.H3"] or 0) + (r["PP.H4"] or 0)
        if abs(pp_sum - 1.0) > 1e-6:
            failures.append(f"{tag}: PP sum = {pp_sum:.6f} (expected 1.0±1e-6)")
        for h in ("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"):
            v = r[h]
            if v is None or (isinstance(v, float) and (math.isnan(v) or v < 0)):
                failures.append(f"{tag}: {h} invalid ({v})")
        if r["n_variants"] is not None and r["n_variants"] < 500:
            failures.append(f"{tag}: only {r['n_variants']} variants in region (<500 minimum)")
        if r.get("h3_negative_flag"):
            failures.append(f"{tag}: H3 distinct likelihood went negative (numerical instability)")
        if r.get("underflow_flag"):
            failures.append(f"{tag}: ABF underflow detected on {r.get('underflow_count', '?')} variants")
        if r["r_coloc_status"] == "available" and r["r_diff_passes_1e3"] is False:
            failures.append(f"{tag}: R vs Python coloc max|ΔPP| = {r['r_python_max_abs_diff']:.4f} (>1e-3)")

    pl.DataFrame(rows).write_csv(OUT_PREFLIGHT, separator="\t")
    # NOTE: TMP_DIR not cleaned in --preflight so re-run debug doesn't re-download.

    print("---PREFLIGHT VALIDATION SUMMARY---", flush=True)
    for r in rows:
        ppstr = (f"H0={r['PP.H0']:.3f} H1={r['PP.H1']:.3f} H2={r['PP.H2']:.3f} "
                 f"H3={r['PP.H3']:.3f} H4={r['PP.H4']:.3f}")
        print(f"  {r['gene']:5s} × {r['outcome']:25s}  n={r['n_variants']:5d}  {ppstr}", flush=True)
        print(f"    sign_flipped={r['n_sign_flipped']}  palindromic_dropped={r['n_dropped_palindromic']}  "
              f"underflow={r['underflow_count']}  R={r['r_coloc_status']}  Δmax={r['r_python_max_abs_diff']}", flush=True)
    if failures:
        print("\nVALIDATION FAILURES:", flush=True)
        for f in failures:
            print(f"  - {f}", flush=True)
        print("\nFAIL — do not run --run-wave1", flush=True)
        sys.exit(2)
    print("\nPASS — proceed with --run-wave1", flush=True)
    return 0


def cmd_run_wave1():
    target_apts = {p[0] for p in WAVE1_PAIRS}
    leads, apt_paths, coord_cache, instr, caches = init_environment(target_apts)

    coloc_rows, qc_rows, artifact_rows = [], [], []
    for apt, gene, outcome in WAVE1_PAIRS:
        for prior in (PRIOR_PRIMARY, PRIOR_SENSITIVITY):
            res, qc, art, _aligned = run_pair(
                apt, gene, outcome, prior,
                pqtl_region_cache=caches["region"],
                pqtl_region_b37_cache=caches["region_b37"],
                apt_paths=apt_paths, leads=leads,
                outcome_cache=caches["outcome"],
                coord_cache=coord_cache, instr=instr,
            )
            coloc_rows.append(res)
            if prior["name"] == PRIOR_PRIMARY["name"]:
                qc_rows.append(qc)
                artifact_rows.append({
                    "gene": gene, "aptamer": apt, "outcome": outcome,
                    "lead_pqtl_intragenic": art["intragenic"],
                    "multi_aptamer_status": art["multi_aptamer_status"],
                    "lead_pqtl_protein_altering": "TBD_VEP_required",
                    "aptamer_artifact_risk": ("low" if art["intragenic"]
                                              and art["multi_aptamer_status"] == "single"
                                              else "review"),
                })
        _log(f"  pair {gene}-{outcome} done (primary + sensitivity)")

    write_results(coloc_rows, qc_rows, artifact_rows)
    write_summary_md(coloc_rows)
    cleanup_tmp()
    return 0


def main() -> int:
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--preflight", action="store_true")
    g.add_argument("--run-wave1", action="store_true")
    args = p.parse_args()
    if args.preflight:
        return cmd_preflight()
    return cmd_run_wave1()


if __name__ == "__main__":
    sys.exit(main())
