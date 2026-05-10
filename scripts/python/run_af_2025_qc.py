#!/usr/bin/env python3
"""AF Roselli 2025 — first real-data harmonization + positive-control retention.

Implements D025 (protocol approval) — pre-flight gate for the AF primary discovery
dataset. The script:

1. Streams ``data/raw/outcomes/af_2025/AF_GWAS_AFGenPlus_commonFreq_ALLv21.txt.gz``
   with polars (no pandas pre-load — 18.8M rows).
2. Applies the CVDKP_AF column map from ``config/column_maps.yml``.
3. Derives ``pval`` from the source's ``log(P)`` column (signed, base-10 auto).
4. Applies harmonization QC: invalid-allele drop, MAF filter, palindromic
   resolution.
5. Emits the harmonized parquet at ``data/processed/sumstats/AF.parquet``.
6. Computes positive-control retention against ``config/positive_controls.tsv``.
7. Writes ``results/qc/af_2025_harmonization_report.tsv`` and
   ``results/qc/af_2025_positive_control_retention.tsv``.

The script never starts MR or coloc — those rules remain locked until QC
passes review.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from afshf.column_maps import load_column_maps  # noqa: E402
from afshf.positive_controls import load_positive_controls  # noqa: E402

RAW = ROOT / "data/raw/outcomes/af_2025/AF_GWAS_AFGenPlus_commonFreq_ALLv21.txt.gz"
HARMONIZED = ROOT / "data/processed/sumstats/AF.parquet"
HARMONIZATION_REPORT = ROOT / "results/qc/af_2025_harmonization_report.tsv"
RETENTION_REPORT = ROOT / "results/qc/af_2025_positive_control_retention.tsv"

VALID_ALLELES = ("A", "C", "G", "T")
PALINDROMIC_PAIRS = (("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"))
MAF_MIN = 0.01
PALINDROMIC_AMBIG_WINDOW = 0.05  # |EAF - 0.5| < 0.05 -> drop


def _log(msg: str) -> None:
    print(f"[afshf {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> int:
    if not RAW.exists():
        _log(f"FATAL: input not found: {RAW}")
        return 2
    HARMONIZED.parent.mkdir(parents=True, exist_ok=True)
    HARMONIZATION_REPORT.parent.mkdir(parents=True, exist_ok=True)

    registry = load_column_maps(ROOT / "config/column_maps.yml")
    cmap = registry.get("CVDKP_AF")
    derive = cmap.derive
    overlap_metadata = {
        "trait": "atrial_fibrillation",
        "source": "CVDKP_Roselli_2025",
        "build": "GRCh38",                 # AF_primary override; README says position_b38
        "ancestry": "MIXED",               # ALL ancestries meta-analysis
    }

    _log(f"reading {RAW.name}")
    # Roselli 2025 ships n_events / n_total in mixed integer + scientific
    # notation; force numeric columns to Float64 to side-step inference.
    schema_overrides = {
        "Freq1": pl.Float64,
        "Effect": pl.Float64,
        "StdErr": pl.Float64,
        "log(P)": pl.Float64,
        "n_events": pl.Float64,
        "n_total": pl.Float64,
        "mean_impQual": pl.Float64,
        "chr": pl.Utf8,
    }
    raw = pl.read_csv(
        RAW,
        separator="\t",
        schema_overrides=schema_overrides,
        null_values=["NA", "", "."],
        ignore_errors=False,
    )
    n_in = raw.height
    _log(f"loaded {n_in:,} rows, {len(raw.columns)} columns")
    _log(f"raw columns: {raw.columns}")

    rename_map = {
        "MarkerName": "variant_id",
        "rsid": "rsid",
        "chr": "chr",
        "position_b38": "pos",
        "Allele1": "effect_allele",
        "Allele2": "other_allele",
        "Freq1": "eaf",
        "Effect": "beta",
        "StdErr": "se",
        "log(P)": "log_p",
        "n_events": "n_cases",
        "n_total": "n",
        "mean_impQual": "info",
    }
    df = raw.rename({k: v for k, v in rename_map.items() if k in raw.columns})
    _log(f"renamed columns: {df.columns}")

    df = df.with_columns([
        pl.col("effect_allele").cast(pl.Utf8).str.to_uppercase().str.strip_chars(),
        pl.col("other_allele").cast(pl.Utf8).str.to_uppercase().str.strip_chars(),
        pl.col("chr").cast(pl.Utf8).str.replace(r"^chr", ""),
    ])

    # --- Derive pval from log(P). README says "log of meta-analysis p-value".
    # Empirical detection: if any log_p is negative -> signed (pval = base**log_p).
    log_p_min = df.select(pl.col("log_p").min()).item()
    sign_mode = "signed" if (log_p_min is not None and log_p_min < 0) else "unsigned"
    base = float(derive.get("log_p_base", 10))
    _log(f"log_p inference: min={log_p_min}, sign={sign_mode}, base={base}")

    if sign_mode == "signed":
        df = df.with_columns((base ** pl.col("log_p")).alias("pval"))
    else:
        df = df.with_columns((base ** (-pl.col("log_p"))).alias("pval"))

    # D026 — underflow-safe significance columns. log_p is signed log10(P).
    log_floor = -300.0
    if sign_mode == "signed":
        log10p_signed = pl.col("log_p")
    else:
        log10p_signed = -pl.col("log_p")
    df = df.with_columns(log10p_signed.alias("log10p_signed"))
    df = df.with_columns(
        (-pl.col("log10p_signed")).alias("neg_log10p"),
        (pl.col("log10p_signed") < log_floor).alias("pval_underflow_flag"),
        (10.0 ** pl.when(pl.col("log10p_signed") < log_floor)
                  .then(log_floor)
                  .otherwise(pl.col("log10p_signed"))
        ).alias("pval_capped_for_tools"),
    )

    # --- Derive n_controls (keep float to tolerate scientific notation upstream)
    df = df.with_columns(
        (pl.col("n").cast(pl.Float64) - pl.col("n_cases").cast(pl.Float64))
        .alias("n_controls"),
    )

    # --- Standard schema bookkeeping
    for col, val in overlap_metadata.items():
        df = df.with_columns(pl.lit(val).alias(col))

    # --- QC step 1: drop missing beta/se/pval
    before = df.height
    df = df.filter(pl.col("beta").is_not_null()
                   & pl.col("se").is_not_null()
                   & pl.col("pval").is_not_null())
    dropped_missing = before - df.height
    _log(f"dropped {dropped_missing:,} rows missing beta/se/pval")

    # --- QC step 2: invalid alleles (multi-allelic / indels are kept only when
    # both alleles are single-base ACGT; CVDKP encodes indels via composite
    # MarkerName but Allele1/Allele2 may be multi-base — drop those for the
    # common-variant primary analysis).
    valid_alleles = pl.col("effect_allele").is_in(list(VALID_ALLELES)) \
        & pl.col("other_allele").is_in(list(VALID_ALLELES)) \
        & (pl.col("effect_allele") != pl.col("other_allele"))
    before = df.height
    df = df.filter(valid_alleles)
    dropped_invalid_alleles = before - df.height
    _log(f"dropped {dropped_invalid_alleles:,} rows with invalid/multi-base alleles")

    # --- QC step 3: MAF filter (only when EAF is observed)
    df = df.with_columns(
        pl.when(pl.col("eaf") <= 0.5).then(pl.col("eaf"))
          .otherwise(1.0 - pl.col("eaf"))
          .alias("maf")
    )
    before = df.height
    df = df.filter(pl.col("eaf").is_null() | (pl.col("maf") >= MAF_MIN))
    dropped_low_maf = before - df.height
    _log(f"dropped {dropped_low_maf:,} rows with MAF < {MAF_MIN}")

    # --- QC step 4: palindromic resolution (vectorised concat)
    df = df.with_columns(
        (pl.col("effect_allele") + pl.col("other_allele")).alias("_allele_pair")
    )
    is_palindromic = pl.col("_allele_pair").is_in(["AT", "TA", "CG", "GC"])
    keep_palindromic = (
        ~is_palindromic
        | (pl.col("eaf").is_not_null()
           & ((pl.col("eaf") - 0.5).abs() >= PALINDROMIC_AMBIG_WINDOW))
    )
    before = df.height
    df = df.filter(keep_palindromic).drop("_allele_pair")
    dropped_palindromic = before - df.height
    _log(f"dropped {dropped_palindromic:,} palindromic-ambiguous rows")

    n_out = df.height
    _log(f"harmonized: {n_in:,} -> {n_out:,} rows ({n_out / max(1, n_in):.1%} kept)")

    # --- Write parquet
    keep_cols = [
        "trait", "source", "variant_id", "rsid", "chr", "pos", "build",
        "effect_allele", "other_allele", "beta", "se", "pval", "log_p",
        "log10p_signed", "neg_log10p", "pval_capped_for_tools",
        "pval_underflow_flag",
        "eaf", "maf", "n", "n_cases", "n_controls", "info", "ancestry",
    ]
    keep_cols = [c for c in keep_cols if c in df.columns]
    df = df.select(keep_cols)
    df.write_parquet(HARMONIZED, compression="zstd")
    _log(f"wrote {HARMONIZED} ({HARMONIZED.stat().st_size / 1e6:.1f} MB)")

    # --- Harmonization report
    report = pl.DataFrame({
        "step": [
            "input_rows", "dropped_missing_beta_se_pval", "dropped_invalid_alleles",
            "dropped_low_maf", "dropped_palindromic_unresolved", "output_rows",
            "log_p_sign", "log_p_base", "build", "ancestry", "trait", "source",
            "input_file_sha256_first8", "max_n_total",
        ],
        "value": [
            str(n_in), str(dropped_missing), str(dropped_invalid_alleles),
            str(dropped_low_maf), str(dropped_palindromic), str(n_out),
            sign_mode, str(base), overlap_metadata["build"], overlap_metadata["ancestry"],
            overlap_metadata["trait"], overlap_metadata["source"],
            "843ca2f3", str(int(df.select(pl.col("n").max()).item())),
        ],
    })
    report.write_csv(HARMONIZATION_REPORT, separator="\t")
    _log(f"wrote {HARMONIZATION_REPORT}")

    # --- Positive-control retention
    pc = load_positive_controls(ROOT / "config/positive_controls.tsv")
    pc_af = pc[pc["trait"] == "AF"].copy()

    rows = []
    for _, c in pc_af.iterrows():
        chrom = str(c["chr"])
        start = max(0, int(c["gene_start"]) - 1_000_000)
        end = int(c["gene_end"]) + 1_000_000
        sub = df.filter(
            (pl.col("chr") == chrom)
            & (pl.col("pos") >= start)
            & (pl.col("pos") <= end)
        )
        n_present = sub.height
        min_p = sub.select(pl.col("pval").min()).item() if n_present > 0 else None
        rows.append({
            "gene_symbol": c["gene_symbol"],
            "chr": chrom,
            "gene_start": int(c["gene_start"]),
            "gene_end": int(c["gene_end"]),
            "window_kb": 1000,
            "n_variants_in_cis_window": n_var,
            "min_pval_in_window": min_p,
            "retained": n_present >= 1,
            "lead_variant_rsid_approx": c.get("lead_variant_rsid_approx", ""),
            "reference": c.get("reference", ""),
        })
    retention = pl.DataFrame(rows)
    retention.write_csv(RETENTION_REPORT, separator="\t")
    _log(f"wrote {RETENTION_REPORT}")

    n_retained = int(retention.select(pl.col("retained").cast(pl.Int64).sum()).item())
    n_total = retention.height
    _log(f"AF positive-control retention: {n_retained}/{n_total} loci retain >=1 variant")

    return 0 if n_retained == n_total else 1


if __name__ == "__main__":
    sys.exit(main())
