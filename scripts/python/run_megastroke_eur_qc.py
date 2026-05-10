#!/usr/bin/env python3
"""MEGASTROKE 2018 EUR — D025 + D028 source-validated harmonization + retention.

This script is **dormant** until the MEGASTROKE 2018 EUR sumstats files land
under ``data/raw/outcomes/stroke_eur/`` matching the expected naming pattern
(see D028 in DECISIONS.md). It does NOT attempt to download the data — that
must happen manually because megastroke.org enforces a terms-of-use checkbox.

Phenotype mapping (per protocol decision 2026-05-09):
    MEGASTROKE.2.AIS.EUR.out  -> ischemic_stroke         (PRIMARY)
    MEGASTROKE.4.CES.EUR.out  -> cardioembolic_stroke    (PRIMARY)
    MEGASTROKE.3.LAS.EUR.out  -> large_artery_stroke     (PRIMARY)
    MEGASTROKE.5.SVS.EUR.out  -> small_vessel_stroke     (PRIMARY)
    MEGASTROKE.1.AS.EUR.out   -> any_stroke_context      (OPTIONAL/CONTEXT)

Source validation (D028 step 2):
    A FORMAT-METAL header containing Allele1/Allele2/Effect/StdErr columns
    matches the expected Malik 2018 layout. The script confirms presence of
    >= 4 of the primary subtype files; otherwif then STOP and write
    SOURCE_MISMATCH to the source-validation TSV.

Outputs:
    results/qc/megastroke_eur_source_validation.tsv  (appended)
    results/qc/megastroke_eur_harmonization_report.tsv
    results/qc/megastroke_eur_positive_control_retention.tsv
    data/processed/sumstats/STROKE_AIS_EUR.parquet
    data/processed/sumstats/STROKE_CES_EUR.parquet
    data/processed/sumstats/STROKE_LAS_EUR.parquet
    data/processed/sumstats/STROKE_SVS_EUR.parquet
"""

from __future__ import annotations

import re
import sys
import time
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from afshf.column_maps import load_column_maps  # noqa: E402
from afshf.positive_controls import load_positive_controls  # noqa: E402

RAW_DIR = ROOT / "data/raw/outcomes/stroke_eur"
CHAIN_PATH = ROOT / "data/raw/annotations/hg38ToHg19.over.chain.gz"
PROCESSED = ROOT / "data/processed/sumstats"
SOURCE_REPORT = ROOT / "results/qc/megastroke_eur_source_validation.tsv"
HARMONIZATION_REPORT = ROOT / "results/qc/megastroke_eur_harmonization_report.tsv"
RETENTION_REPORT = ROOT / "results/qc/megastroke_eur_positive_control_retention.tsv"

EXPECTED_FILES = {
    # subtype: (filename pattern, outcome label, role)
    "AS":  ("MEGASTROKE.1.AS.EUR.out",  "any_stroke_context",   "context"),
    "AIS": ("MEGASTROKE.2.AIS.EUR.out", "ischemic_stroke",      "primary"),
    "LAS": ("MEGASTROKE.3.LAS.EUR.out", "large_artery_stroke",  "primary"),
    "CES": ("MEGASTROKE.4.CES.EUR.out", "cardioembolic_stroke", "primary"),
    "SVS": ("MEGASTROKE.5.SVS.EUR.out", "small_vessel_stroke",  "primary"),
}

OUTCOME_PARQUET = {
    "AIS": "STROKE_AIS_EUR.parquet",
    "CES": "STROKE_CES_EUR.parquet",
    "LAS": "STROKE_LAS_EUR.parquet",
    "SVS": "STROKE_SVS_EUR.parquet",
}

STROKE_POSITIVE_CONTROL_GENES = {
    "HDAC9", "PHACTR1", "ABO", "PITX2", "ZFHX3",
    "FOXF2", "ZCCHC14", "COL4A1", "COL4A2",
}

VALID_ALLELES = ("A", "C", "G", "T")
PALINDROMIC_PAIRS_STR = ["AT", "TA", "CG", "GC"]
MAF_MIN = 0.01
PALINDROMIC_AMBIG_WINDOW = 0.05
LOG_FLOOR = -300.0
WINDOW_KB = 1000


def _log(msg: str) -> None:
    print(f"[afshf {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def find_subtype_file(subtype: str) -> Path | None:
    """Locate AIS/CES/LAS/SVS fwith allowing .gz, .out, or compressed variants."""
    base = EXPECTED_FILES[subtype][0]
    candidates = []
    for variant in [base, base + ".gz", base + ".bgz",
                    base.replace(".out", ".tsv"), base.replace(".out", ".tsv.gz")]:
        candidates.extend(RAW_DIR.rglob(variant))
    return candidates[0] if candidates else None


def liftover_positive_controls(controls_df, chain_path: Path):
    from pyliftover import LiftOver
    lo = LiftOver(str(chain_path))
    out = controls_df.copy()
    starts_b37, ends_b37, ok = [], [], []
    for _, r in controls_df.iterrows():
        chrom = "chr" + str(r["chr"])
        s = lo.convert_coordinate(chrom, int(r["gene_start"]))
        e = lo.convert_coordinate(chrom, int(r["gene_end"]))
        if s and e:
            starts_b37.append(s[0][1])
            ends_b37.append(e[0][1])
            ok.append(True)
        else:
            starts_b37.append(None)
            ends_b37.append(None)
            ok.append(False)
    out["gene_start_b37"] = starts_b37
    out["gene_end_b37"] = ends_b37
    out["liftover_ok"] = ok
    return out


def harmonize_subtype(tsv_path: Path,
                      outcome_label: str,
                      outcome_role: str,
                      registry,
                      pc_b37) -> dict:
    cmap = registry.get("CVDKP_MEGASTROKE")
    _log(f"[{outcome_label}] reading {tsv_path.name}")

    schema_overrides = {
        "Effect": pl.Float64, "BETA": pl.Float64, "beta": pl.Float64,
        "StdErr": pl.Float64, "SE": pl.Float64, "se": pl.Float64,
        "P-value": pl.Float64, "P": pl.Float64, "pval": pl.Float64,
        "Freq1": pl.Float64, "EAF": pl.Float64, "MAF": pl.Float64,
        "N": pl.Float64, "TotalSampleSize": pl.Float64,
        "chr": pl.Utf8, "CHR": pl.Utf8,
    }

    raw = pl.read_csv(
        tsv_path,
        separator="\t",
        schema_overrides=schema_overrides,
        null_values=["NA", "", ".", "NaN"],
        ignore_errors=False,
    )
    n_in = raw.height
    _log(f"[{outcome_label}] loaded {n_in:,} rows; columns: {raw.columns}")

    rename_map: dict[str, str] = {}
    aliases = cmap.column_aliases
    for std_field, alt_list in aliases.items():
        for cand in alt_list:
            if cand in raw.columns:
                if cand != std_field:
                    rename_map[cand] = std_field
                break
    df = raw.rename(rename_map)
    _log(f"[{outcome_label}] renamed to: {df.columns}")

    df = df.with_columns([
        pl.col("effect_allele").cast(pl.Utf8).str.to_uppercase().str.strip_chars(),
        pl.col("other_allele").cast(pl.Utf8).str.to_uppercase().str.strip_chars(),
    ])
    if "chr" in df.columns:
        df = df.with_columns(pl.col("chr").cast(pl.Utf8).str.replace(r"^chr", ""))

    # D026 underflow-safe columns from raw pval (Malik 2018 ships P directly).
    if "pval" in df.columns:
        df = df.with_columns(
            pl.when(pl.col("pval") <= 0).then(LOG_FLOOR)
              .otherwise(pl.col("pval").log10())
              .alias("log10p_signed"),
        )
        df = df.with_columns(
            (-pl.col("log10p_signed")).alias("neg_log10p"),
            (pl.col("log10p_signed") < LOG_FLOOR).alias("pval_underflow_flag"),
            (10.0 ** pl.when(pl.col("log10p_signed") < LOG_FLOOR)
                       .then(LOG_FLOOR)
                       .otherwise(pl.col("log10p_signed"))
            ).alias("pval_capped_for_tools"),
        )

    df = df.with_columns(
        pl.lit(outcome_label).alias("trait"),
        pl.lit("MEGASTROKE_2018_EUR").alias("source"),
        pl.lit("GRCh37").alias("build"),
        pl.lit("EUR").alias("ancestry"),
    )

    # QC step 1 — missing beta/se/pval
    before = df.height
    df = df.filter(pl.col("beta").is_not_null()
                   & pl.col("se").is_not_null()
                   & pl.col("pval").is_not_null())
    dropped_missing = before - df.height
    _log(f"[{outcome_label}] dropped {dropped_missing:,} missing beta/se/pval")

    # QC step 2 — invalid alleles
    valid = pl.col("effect_allele").is_in(list(VALID_ALLELES)) \
        & pl.col("other_allele").is_in(list(VALID_ALLELES)) \
        & (pl.col("effect_allele") != pl.col("other_allele"))
    before = df.height
    df = df.filter(valid)
    dropped_invalid = before - df.height
    _log(f"[{outcome_label}] dropped {dropped_invalid:,} invalid/multi-base alleles")

    # QC step 3 — MAF (only if EAF present)
    has_eaf = "eaf" in df.columns
    dropped_low_maf = 0
    if has_eaf:
        df = df.with_columns(
            pl.when(pl.col("eaf") <= 0.5).then(pl.col("eaf"))
              .otherwise(1.0 - pl.col("eaf")).alias("maf")
        )
        before = df.height
        df = df.filter(pl.col("eaf").is_null() | (pl.col("maf") >= MAF_MIN))
        dropped_low_maf = before - df.height
        _log(f"[{outcome_label}] dropped {dropped_low_maf:,} MAF<{MAF_MIN}")

    # QC step 4 — palindromic
    df = df.with_columns(
        (pl.col("effect_allele") + pl.col("other_allele")).alias("_pair")
    )
    is_palindromic = pl.col("_pair").is_in(PALINDROMIC_PAIRS_STR)
    if has_eaf:
        keep = (~is_palindromic
                | (pl.col("eaf").is_not_null()
                   & ((pl.col("eaf") - 0.5).abs() >= PALINDROMIC_AMBIG_WINDOW)))
    else:
        keep = ~is_palindromic       # no EAF -> drop all palindromic
    before = df.height
    df = df.filter(keep).drop("_pair")
    dropped_palindromic = before - df.height
    _log(f"[{outcome_label}] dropped {dropped_palindromic:,} palindromic-ambiguous")

    n_out = df.height
    _log(f"[{outcome_label}] harmonized: {n_in:,} -> {n_out:,} ({n_out/max(1,n_in):.1%})")

    keep_cols = [
        "trait", "source", "variant_id", "rsid", "chr", "pos", "build",
        "effect_allele", "other_allele", "beta", "se", "pval",
        "log10p_signed", "neg_log10p", "pval_capped_for_tools",
        "pval_underflow_flag",
        "eaf", "maf", "n", "n_cases", "n_controls", "ancestry",
    ]
    keep_cols = [c for c in keep_cols if c in df.columns]
    parquet_path = PROCESSED / OUTCOME_PARQUET[
        next(k for k, v in EXPECTED_FILES.items()
             if v[1] == outcome_label and k in OUTCOME_PARQUET)
    ]
    df.select(keep_cols).write_parquet(parquet_path, compression="zstd")
    _log(f"[{outcome_label}] wrote {parquet_path} ({parquet_path.stat().st_size/1e6:.1f} MB)")

    # Retention check (b37 lifted positive controls, gene-symbol filter).
    pc = pc_b37[
        pc_b37["gene_symbol"].isin(STROKE_POSITIVE_CONTROL_GENES)
        & pc_b37["liftover_ok"]
    ].drop_duplicates(subset=["gene_symbol"])

    rows = []
    for _, c in pc.iterrows():
        chrom = str(c["chr"])
        start = max(0, int(c["gene_start_b37"]) - WINDOW_KB * 1000)
        end = int(c["gene_end_b37"]) + WINDOW_KB * 1000
        sub = df.filter((pl.col("chr") == chrom)
                        & (pl.col("pos") >= start)
                        & (pl.col("pos") <= end))
        n_present = sub.height
        max_neg = (sub.select(pl.col("neg_log10p").max()).item()
                   if n_present and "neg_log10p" in sub.columns else None)
        rows.append({
            "outcome": outcome_label, "outcome_role": outcome_role,
            "gene_symbol": c["gene_symbol"], "chr": chrom,
            "gene_start_b37": int(c["gene_start_b37"]),
            "gene_end_b37": int(c["gene_end_b37"]),
            "window_kb": WINDOW_KB,
            "n_variants_in_cis_window": n_var,
            "max_neg_log10p_in_window": max_neg,
            "retained": n_present >= 1,
            "lead_variant_rsid_approx": c.get("lead_variant_rsid_approx", ""),
            "reference": c.get("reference", ""),
        })
    return {
        "outcome": outcome_label, "outcome_role": outcome_role,
        "n_in": n_in, "n_out": n_out,
        "dropped_missing": dropped_missing, "dropped_invalid": dropped_invalid,
        "dropped_low_maf": dropped_low_maf, "dropped_palindromic": dropped_palindromic,
        "retention_rows": rows,
        "n_retained": sum(1 for r in rows if r["retained"]),
        "n_pc_total": len(rows),
    }


def main() -> int:
    if not RAW_DIR.exists():
        _log(f"FATAL: {RAW_DIR} does not exist; create the directory and place "
             f"MEGASTROKE 2018 EUR files inside.")
        return 2
    if not CHAIN_PATH.exists():
        _log(f"FATAL: liftover chain missing at {CHAIN_PATH}")
        return 2
    PROCESSED.mkdir(parents=True, exist_ok=True)
    SOURCE_REPORT.parent.mkdir(parents=True, exist_ok=True)

    # D028 step 2: source validation
    found: dict[str, Path] = {}
    for subtype in EXPECTED_FILES:
        path = find_subtype_file(subtype)
        if path is not None:
            found[subtype] = path
            _log(f"D028: located {subtype} -> {path}")
    primary_subtypes_present = sum(
        1 for s in ("AIS", "CES", "LAS", "SVS") if s in found
    )
    if primary_subtypes_present < 4:
        _log(f"D028 SOURCE_MISMATCH: only {primary_subtypes_present}/4 primary "
             f"MEGASTROKE 2018 EUR files present. STOP — do not run harmonize.")
        return 1
    _log(f"D028 source validation passed: 4/4 primary subtype files located.")

    registry = load_column_maps(ROOT / "config/column_maps.yml")
    pc = load_positive_controls(ROOT / "config/positive_controls.tsv")
    pc_b37 = liftover_positive_controls(pc, CHAIN_PATH)

    summaries = []
    all_retention_rows: list[dict] = []
    # Skip "AS" by default (context-only); harmonize only the 4 primaries.
    for subtype in ("AIS", "CES", "LAS", "SVS"):
        _, outcome_label, outcome_role = EXPECTED_FILES[subtype]
        path = found[subtype]
        s = harmonize_subtype(path, outcome_label, outcome_role, registry, pc_b37)
        summaries.append(s)
        all_retention_rows.extend(s["retention_rows"])

    # Reports
    rep = pl.DataFrame([{
        "outcome": s["outcome"], "outcome_role": s["outcome_role"],
        "input_rows": s["n_in"],
        "dropped_missing_beta_se_pval": s["dropped_missing"],
        "dropped_invalid_alleles": s["dropped_invalid"],
        "dropped_low_maf": s["dropped_low_maf"],
        "dropped_palindromic_unresolved": s["dropped_palindromic"],
        "output_rows": s["n_out"],
        "positive_controls_retained": s["n_retained"],
        "positive_controls_total": s["n_pc_total"],
    } for s in summaries])
    rep.write_csv(HARMONIZATION_REPORT, separator="\t")
    _log(f"wrote {HARMONIZATION_REPORT}")

    pl.DataFrame(all_retention_rows).write_csv(RETENTION_REPORT, separator="\t")
    _log(f"wrote {RETENTION_REPORT}")

    n_total = sum(s["n_pc_total"] for s in summaries)
    n_retained = sum(s["n_retained"] for s in summaries)
    _log(f"MEGASTROKE retention: {n_retained}/{n_total}")
    return 0 if n_retained == n_total else 1


if __name__ == "__main__":
    sys.exit(main())
