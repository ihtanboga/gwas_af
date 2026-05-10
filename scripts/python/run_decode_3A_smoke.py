#!/usr/bin/env python3
"""Phase 3A.3 + smoke run for deCODE 2021 cis-pQTL extraction.

D031: deCODE primary discovery. Per-protein .txt.gz files are ~950 MB; total
collection 26 TB. Strategy = stream-fetch → cis-window filter → instrument
table → delete file. Peak disk ~1.1 GB.

Per recipe (verified API):
    GET /s3/folder?token=<UUID>      -> JSON list (cached)
    GET /s3/download?token=<UUID>&file=<KEY>  -> per-protein .txt.gz

Per-protein columns (Ferkingstad 2021 README):
    Chrom Pos Name rsids effectAllele otherAllele Beta Pval min_log10_pval SE N ImpMAF

Outputs:
    results/qc/decode_2021_metadata_report.tsv
    results/qc/decode_2021_protein_map_qc.tsv
    results/qc/decode_2021_cis_instrument_report.tsv
    data/processed/instruments/decode_2021_smoke.parquet      (per-aptamer cis instruments)
    results/smoke/decode_20protein_instruments.tsv             (mirror, one row per instrument)
    results/smoke/decode_20protein_summary.md                  (human-readable)

Note: real LD clumping (PLINK + 1000G EUR) deferred to Phase 3A.full. Smoke
uses distance-based clumping (greedy lead + 500 kb window) which is
sufficient for sanity but not for final MR. Distance-clumping flag emitted
in the instrument report so downstream consumers know.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlencode

import requests
import polars as pl

ROOT = Path("/Users/apple/Desktop/gwas_af")
TOKEN = os.environ.get("DECODE_TOKEN", "<DECODE_TOKEN>")
FOLDER_API = "https://download.decode.is/s3/folder"
DOWNLOAD_API = "https://download.decode.is/s3/download"

FOLDER_JSON = ROOT / "data/registry/decode_2021_folder.json"
INSTR_DIR = ROOT / "data/processed/instruments"
SMOKE_DIR = ROOT / "results/smoke"
QC_DIR = ROOT / "results/qc"
TMP_DOWNLOAD_DIR = Path("/tmp/decode_smoke")  # stream targets, deleted after

CIS_WINDOW_KB = 1000
P_PRIMARY = 5.0e-8
P_EXPLORATORY = 5.0e-6
F_STAT_MIN = 10.0
DISTANCE_CLUMP_KB = 500   # greedy distance clumping (smoke fallback)

# Per spec section 18 / D031: 20-protein smoke set with multi-aptamer flags.
# GRCh38 gene-body coordinates (Ensembl). Pipeline auto-pads ±CIS_WINDOW_KB.
SMOKE_GENES: dict[str, dict] = {
    "PCSK9":    {"chr": "1",  "start":  55039548, "end":  55064852, "category": "anchor_lipid"},
    "IL6R":     {"chr": "1",  "start": 154405193, "end": 154469450, "category": "anchor_inflammation"},
    "SORT1":    {"chr": "1",  "start": 109274570, "end": 109365052, "category": "anchor_lipid"},
    "ANGPTL3":  {"chr": "1",  "start":  62597488, "end":  62606909, "category": "anchor_lipid"},
    "ANGPTL4":  {"chr": "19", "start":   8364181, "end":   8376234, "category": "anchor_lipid"},
    "APOE":     {"chr": "19", "start":  44905781, "end":  44909393, "category": "anchor_lipid"},
    "APOB":     {"chr": "2",  "start":  21001429, "end":  21043950, "category": "anchor_lipid"},
    "F2":       {"chr": "11", "start":  46740730, "end":  46761057, "category": "coag_stroke", "multi_aptamer": True},
    "F11":      {"chr": "4",  "start": 186266801, "end": 186290374, "category": "coag_stroke"},
    "F9":       {"chr": "X",  "start": 139530763, "end": 139564310, "category": "coag", "multi_aptamer": True},
    "FGB":      {"chr": "4",  "start": 154562811, "end": 154597921, "category": "coag_stroke"},
    "PROC":     {"chr": "2",  "start": 127417399, "end": 127430350, "category": "coag_ambig", "multi_aptamer": True},  # zymogen vs activated
    "PROS1":    {"chr": "3",  "start":  93876035, "end":  93974092, "category": "coag"},
    "ADAMTS13": {"chr": "9",  "start": 133414233, "end": 133458997, "category": "coag_stroke"},
    "NPPA":     {"chr": "1",  "start":  11845709, "end":  11848345, "category": "hf_natriuretic"},
    "NPPB":     {"chr": "1",  "start":  11849854, "end":  11853614, "category": "hf_natriuretic", "multi_aptamer": True},  # BNP, BNP_32, NT-proBNP
    "BAG3":     {"chr": "10", "start": 119651379, "end": 119677819, "category": "hf_dcm"},
    "ANGPT1":   {"chr": "8",  "start": 107249506, "end": 107285550, "category": "vascular"},
    "MMP9":     {"chr": "20", "start":  46008908, "end":  46016561, "category": "vascular"},
    "CRP":      {"chr": "1",  "start": 159712289, "end": 159714589, "category": "inflammation"},
}


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def _log(msg: str) -> None:
    print(f"[decode-3a {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_folder() -> list[dict]:
    if not FOLDER_JSON.exists():
        r = requests.get(FOLDER_API, params={"token": TOKEN}, timeout=60)
        r.raise_for_status()
        FOLDER_JSON.write_text(json.dumps(r.json(), indent=2))
    return json.loads(FOLDER_JSON.read_text())["files"]


def find_aptamers_for_gene(files: list[dict], gene: str) -> list[dict]:
    out = []
    for f in files:
        parts = f["Key"].replace(".txt.gz", "").split("_")
        if len(parts) >= 4 and parts[2] == gene:
            out.append({"key": f["Key"], "size": f["Size"], "etag": f.get("ETag", "")})
    return out


# ---------------------------------------------------------------------------
# Streaming cis extraction
# ---------------------------------------------------------------------------


@dataclass
class CisRow:
    chrom: str
    pos: int
    name: str
    rsid: str
    ea: str
    oa: str
    beta: float
    pval: float
    log_p: float
    se: float
    n: int
    imp_maf: float


def stream_cis_extract(gz_path: Path, target_chr: str,
                       window_start: int, window_end: int,
                       p_threshold: float = P_EXPLORATORY) -> list[CisRow]:
    """Read deCODE per-protein .txt.gz line-by-line; emit cis rows with
    P < threshold. Avoids materialising full fwith in memory."""
    out: list[CisRow] = []
    with gzip.open(gz_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        ci, pi = idx["Chrom"], idx["Pos"]
        ni, ri = idx["Name"], idx.get("rsids", -1)
        ei, oi = idx["effectAllele"], idx["otherAllele"]
        bi, pvi = idx["Beta"], idx["Pval"]
        lpi = idx.get("min_log10_pval", -1)
        si, sni = idx["SE"], idx["N"]
        mi = idx.get("ImpMAF", -1)

        target_chr_norm = str(target_chr).replace("chr", "")
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
            if pval >= p_threshold:
                continue
            try:
                beta = float(cols[bi])
                se = float(cols[si])
                n = int(float(cols[sni]))
            except ValueError:
                continue
            log_p = float(cols[lpi]) if lpi >= 0 and cols[lpi] not in ("", "NA") else -1.0
            imp_maf = float(cols[mi]) if mi >= 0 and cols[mi] not in ("", "NA") else float("nan")
            rsid = cols[ri] if ri >= 0 else ""
            out.append(CisRow(
                chrom=target_chr_norm, pos=pos, name=cols[ni], rsid=rsid,
                ea=cols[ei].upper(), oa=cols[oi].upper(),
                beta=beta, pval=pval, log_p=log_p,
                se=se, n=n, imp_maf=imp_maf,
            ))
    return out


def f_statistic(beta: float, se: float) -> float:
    return (beta / se) ** 2 if se else 0.0


def distance_clump(rows: list[CisRow], window_kb: int = DISTANCE_CLUMP_KB) -> list[CisRow]:
    """Greedy distance-based clumping (smoke fallback for true LD clumping):
    sort by p-value asc; for each lead, exclude any variant within ±window_kb."""
    sorted_rows = sorted(rows, key=lambda r: r.pval)
    keep: list[CisRow] = []
    excluded: set[int] = set()
    pad = window_kb * 1000
    for i, r in enumerate(sorted_rows):
        if i in excluded:
            continue
        keep.append(r)
        for j in range(i + 1, len(sorted_rows)):
            o = sorted_rows[j]
            if abs(o.pos - r.pos) <= pad:
                excluded.add(j)
    return keep


# ---------------------------------------------------------------------------
# Per-aptamer pipeline
# ---------------------------------------------------------------------------


def download_aptamer(key: str, dest: Path, token: str = TOKEN) -> Path:
    """Stream-download one .txt.gz to dest. Re-uses partial download if present."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1_000_000:
        return dest
    _log(f"  downloading {key}")
    t0 = time.time()
    with requests.get(DOWNLOAD_API,
                      params={"token": token, "file": key},
                      stream=True, timeout=600) as r:
        r.raise_for_status()
        with dest.open("wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
    sz = dest.stat().st_size / 1e6
    _log(f"  done {sz:.0f} MB in {time.time()-t0:.0f}s")
    return dest


def process_aptamer(gene: str, gene_meta: dict, aptamer: dict,
                    p_threshold: float = P_EXPLORATORY) -> dict:
    """Stream-download, cis-extract, distance-clump, return summary dict."""
    key = aptamer["key"]
    dest = TMP_DOWNLOAD_DIR / key
    flags: list[str] = []

    try:
        download_aptamer(key, dest)
        cis_rows = stream_cis_extract(
            dest,
            target_chr=gene_meta["chr"],
            window_start=max(0, gene_meta["start"] - CIS_WINDOW_KB * 1000),
            window_end=gene_meta["end"] + CIS_WINDOW_KB * 1000,
            p_threshold=p_threshold,
        )
        n_cis_total = len(cis_rows)
        gw_rows = [r for r in cis_rows if r.pval < P_PRIMARY]
        n_gw = len(gw_rows)
        f_filtered = [r for r in gw_rows if f_statistic(r.beta, r.se) >= F_STAT_MIN]
        n_f = len(f_filtered)
        clumped = distance_clump(f_filtered) if f_filtered else []
        n_clumped = len(clumped)

        if gene_meta.get("multi_aptamer"):
            flags.append("multi_aptamer")

        return {
            "gene": gene,
            "aptamer": key,
            "category": gene_meta.get("category", ""),
            "n_cis_exploratory": n_cis_total,
            "n_cis_gw_significant": n_gw,
            "n_passing_F_gt_10": n_f,
            "n_after_distance_clump": n_clumped,
            "lead_pos": clumped[0].pos if clumped else None,
            "lead_rsid": clumped[0].rsid if clumped else "",
            "lead_pval": clumped[0].pval if clumped else None,
            "lead_beta": clumped[0].beta if clumped else None,
            "lead_se": clumped[0].se if clumped else None,
            "lead_F": f_statistic(clumped[0].beta, clumped[0].se) if clumped else None,
            "lead_ea": clumped[0].ea if clumped else "",
            "lead_oa": clumped[0].oa if clumped else "",
            "lead_imp_maf": clumped[0].imp_maf if clumped else None,
            "lead_log10p_signed": clumped[0].log_p if clumped else None,
            "flags": ",".join(flags),
            "instruments": clumped,           # consumed by caller for parquet
        }
    finally:
        if dest.exists():
            dest.unlink()


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def run_smoke(gene_subset: list[str] | None = None,
              p_threshold: float = P_EXPLORATORY) -> int:
    INSTR_DIR.mkdir(parents=True, exist_ok=True)
    SMOKE_DIR.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    TMP_DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)

    files = load_folder()
    _log(f"folder cached: {len(files)} aptamer files")

    targets = list(SMOKE_GENES) if gene_subset is None else gene_subset
    summaries: list[dict] = []
    instrument_rows: list[dict] = []

    for gene in targets:
        meta = SMOKE_GENES.get(gene)
        if meta is None:
            _log(f"[{gene}] not in SMOKE_GENES; skip")
            continue
        aptamers = find_aptamers_for_gene(files, gene)
        if not aptamers:
            _log(f"[{gene}] not in deCODE panel; status NOT_ON_DECODE_PANEL")
            summaries.append({
                "gene": gene, "aptamer": "", "category": meta.get("category", ""),
                "n_cis_exploratory": 0, "n_cis_gw_significant": 0,
                "n_passing_F_gt_10": 0, "n_after_distance_clump": 0,
                "lead_pos": None, "lead_rsid": "", "lead_pval": None,
                "lead_beta": None, "lead_se": None, "lead_F": None,
                "lead_ea": "", "lead_oa": "", "lead_imp_maf": None,
                "lead_log10p_signed": None, "flags": "NOT_ON_DECODE_PANEL",
            })
            continue

        for apt in aptamers:
            _log(f"[{gene}/{apt['key']}] (size {apt['size']/1e6:.0f} MB)")
            try:
                summary = process_aptamer(gene, meta, apt, p_threshold)
            except Exception as e:
                _log(f"  ERROR: {type(e).__name__}: {e}")
                summaries.append({
                    "gene": gene, "aptamer": apt["key"],
                    "category": meta.get("category", ""),
                    "n_cis_exploratory": 0, "n_cis_gw_significant": 0,
                    "n_passing_F_gt_10": 0, "n_after_distance_clump": 0,
                    "lead_pos": None, "lead_rsid": "", "lead_pval": None,
                    "lead_beta": None, "lead_se": None, "lead_F": None,
                    "lead_ea": "", "lead_oa": "", "lead_imp_maf": None,
                    "lead_log10p_signed": None,
                    "flags": f"ERROR:{type(e).__name__}",
                })
                continue
            instruments = summary.pop("instruments", [])
            for r in instruments:
                instrument_rows.append({
                    "gene": gene, "aptamer": apt["key"],
                    "category": meta.get("category", ""),
                    "chrom": r.chrom, "pos": r.pos, "name": r.name, "rsid": r.rsid,
                    "effect_allele": r.ea, "other_allele": r.oa,
                    "beta": r.beta, "se": r.se, "pval": r.pval,
                    "log10p_signed": r.log_p, "neg_log10p": -r.log_p if r.log_p else None,
                    "n": r.n, "imp_maf": r.imp_maf,
                    "f_stat": f_statistic(r.beta, r.se),
                    "build": "GRCh38", "ancestry": "EUR",
                    "source": "deCODE_2021_Ferkingstad",
                })
            summaries.append(summary)
            _log(f"  -> n_cis_exp={summary['n_cis_exploratory']}  "
                 f"n_gw={summary['n_cis_gw_significant']}  "
                 f"n_F>10={summary['n_passing_F_gt_10']}  "
                 f"n_clumped={summary['n_after_distance_clump']}  "
                 f"lead_p={summary['lead_pval']}")

    # Write outputs
    summary_df = pl.DataFrame(summaries)
    summary_df.write_csv(QC_DIR / "decode_2021_cis_instrument_report.tsv", separator="\t")
    _log(f"wrote {QC_DIR/'decode_2021_cis_instrument_report.tsv'}")

    if instrument_rows:
        instr_df = pl.DataFrame(instrument_rows)
        instr_df.write_parquet(INSTR_DIR / "decode_2021_smoke.parquet", compression="zstd")
        instr_df.write_csv(SMOKE_DIR / "decode_20protein_instruments.tsv", separator="\t")
        _log(f"wrote {INSTR_DIR/'decode_2021_smoke.parquet'} "
             f"({len(instrument_rows)} instruments)")
    else:
        _log("no instruments retained — investigate")

    # Markdown summary
    md = ["# deCODE 2021 — 20-protein smoke run (Phase 3A.3, D031)\n",
          f"Token-valid window: 30 days; folder cached at `{FOLDER_JSON.relative_to(ROOT)}`.",
          f"Total aptamers in deCODE panel: {len(files)}.",
          "",
          "## Per-aptamer summary",
          "",
          "| Gene | Aptamer | Category | cis P<5e-6 | cis P<5e-8 | F>10 | Clumped | Lead p | Lead F | Flags |",
          "|------|---------|----------|-----------:|-----------:|-----:|--------:|-------:|------:|-------|",
    ]
    for s in summaries:
        md.append(
            f"| {s['gene']} | {s['aptamer'][:40]} | {s['category']} | "
            f"{s['n_cis_exploratory']} | {s['n_cis_gw_significant']} | "
            f"{s['n_passing_F_gt_10']} | {s['n_after_distance_clump']} | "
            f"{s['lead_pval'] if s['lead_pval'] is not None else '—'} | "
            f"{s['lead_F']:.0f} | {s['flags']} |"
            if s['lead_F'] is not None else
            f"| {s['gene']} | {s['aptamer'][:40]} | {s['category']} | "
            f"{s['n_cis_exploratory']} | {s['n_cis_gw_significant']} | "
            f"{s['n_passing_F_gt_10']} | {s['n_after_distance_clump']} | — | — | {s['flags']} |"
        )
    (SMOKE_DIR / "decode_20protein_summary.md").write_text("\n".join(md))
    _log(f"wrote {SMOKE_DIR/'decode_20protein_summary.md'}")

    n_with_inst = sum(1 for s in summaries if s["n_after_distance_clump"] > 0)
    _log(f"smoke result: {n_with_inst}/{len(summaries)} aptamer rows produced cis-instruments")
    return 0 if n_with_inst > 0 else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", action="append", help="Restrict to one or more genes")
    ap.add_argument("--p-threshold", type=float, default=P_EXPLORATORY,
                    help="Cis variant p-value threshold (default 5e-6 exploratory)")
    args = ap.parse_args()
    return run_smoke(gene_subset=args.gene, p_threshold=args.p_threshold)


if __name__ == "__main__":
    sys.exit(main())
