#!/usr/bin/env python3
"""HERMES 2.0 / 2024 EUR — first real-data harmonization + retention QC.

Implements D025 for the HF discovery panel (per protocol decision 2026-05-09):

  Pheno1_EUR -> HF_overall      (primary)
  Pheno2_EUR -> HF_nonischemic  (primary)
  Pheno3_EUR -> HF_ni_HFrEF     (primary)
  Pheno4_EUR -> HF_ni_HFpEF     (exploratory, D001 REVISED)

Bundle layout (per README):
  HERMES2_GWAS_HF_EUR.zip
    Pheno1_EUR/
      FORMAT-METAL_Pheno1_EUR.tsv.gz
      FORMAT-METAL_Pheno1_EUR.tsv.gz.tbi
      FORMAT-METAL_Pheno1_EUR.tsv.gz.md5
    Pheno2_EUR/ ...
    Pheno3_EUR/ ...
    Pheno4_EUR/ ...

Build = GRCh37 (pos_b37). For the QC retention check we lift the b38 positive
control coordinates over to b37 with pyliftover. A full b37->b38 liftover of
the variant table is deferred to Phase 1 (when AF/HF/stroke must be merged).

Effect column scale: A1_beta is log OR. logP is unsigned -log10(P).
The pipeline reuses ``apply_column_map(CVDKP_HERMES)`` -> standard schema +
D026 underflow-safe columns (log10p_signed, neg_log10p, pval_capped_for_tools,
pval_underflow_flag).
"""

from __future__ import annotations

import sys
import time
import zipfile
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from afshf.column_maps import load_column_maps  # noqa: E402
from afshf.positive_controls import load_positive_controls  # noqa: E402

ZIP_PATH = ROOT / "data/raw/outcomes/hf_2024_eur/HERMES2_GWAS_HF_EUR.zip"
EXTRACT_DIR = ROOT / "data/raw/outcomes/hf_2024_eur/extracted"
CHAIN_PATH = ROOT / "data/raw/annotations/hg38ToHg19.over.chain.gz"

PROCESSED = ROOT / "data/processed/sumstats"
HARMONIZATION_REPORT = ROOT / "results/qc/hf_2024_eur_harmonization_report.tsv"
RETENTION_REPORT = ROOT / "results/qc/hf_2024_eur_positive_control_retention.tsv"

PHENOTYPE_MAP = {
    "Pheno1": ("HF_overall",      "primary"),
    "Pheno2": ("HF_nonischemic",  "primary"),
    "Pheno3": ("HF_ni_HFrEF",     "primary"),
    "Pheno4": ("HF_ni_HFpEF",     "exploratory"),
}

# HF positive controls — gene-symbol filter; all entries with trait == HF_overall
# in config/positive_controls.tsv apply to every HF subtype QC.
HF_POSITIVE_CONTROL_GENES = {
    "BAG3", "TTN", "NPPA", "NPPB", "FLNC", "TBX5",
    "MYH7", "LMNA", "PLN", "TNNT2",
}

VALID_ALLELES = ("A", "C", "G", "T")
PALINDROMIC_PAIRS = (("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"))
MAF_MIN = 0.01
PALINDROMIC_AMBIG_WINDOW = 0.05
LOG_FLOOR = -300.0
WINDOW_KB = 1000


def _log(msg: str) -> None:
    print(f"[afshf {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def liftover_positive_controls(controls_df, chain_path: Path):
    """Convert positive control gene_start / gene_end b38 -> b37 in-place.

    Returns a DataFrame with new ``gene_start_b37`` / ``gene_end_b37`` columns.
    Drops rows for which liftover fails on either boundary.
    """
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


def harmonize_phenotype(tsv_path: Path,
                         outcome_label: str,
                         outcome_role: str,
                         registry,
                         positive_controls_b37) -> dict:
    """Stream a single FORMAT-METAL phenotype TSV.gz and emit:
       - harmonized parquet at data/processed/sumstats/{outcome_label}.parquet
       - per-phenotype harmonization counters
       - per-phenotype positive-control retention rows
    """
    cmap = registry.get("CVDKP_HERMES")

    _log(f"[{outcome_label}] reading {tsv_path.name}")
    schema_overrides = {
        "chr": pl.Utf8,
        "pos_b37": pl.Int64,
        "A1": pl.Utf8,
        "A2": pl.Utf8,
        "A1_beta": pl.Float64,
        "A1_freq": pl.Float64,
        "se": pl.Float64,
        "pval": pl.Float64,
        "logP": pl.Float64,
        "N_case": pl.Float64,
        "N_total": pl.Float64,
        "isq_het": pl.Float64,
        "p_het": pl.Float64,
    }
    raw = pl.read_csv(
        tsv_path,
        separator="\t",
        schema_overrides=schema_overrides,
        null_values=["NA", "", "."],
        ignore_errors=False,
    )
    n_in = raw.height
    _log(f"[{outcome_label}] loaded {n_in:,} rows; columns: {raw.columns}")

    rename_map = {
        "#key": "variant_id",
        "rsID": "rsid",
        "chr": "chr",
        "pos_b37": "pos",
        "A1": "effect_allele",
        "A2": "other_allele",
        "A1_beta": "beta",
        "A1_freq": "eaf",
        "se": "se",
        "pval": "pval",
        "logP": "log10p",       # unsigned -log10(P) per README
        "N_case": "n_cases",
        "N_total": "n",
    }
    df = raw.rename({k: v for k, v in rename_map.items() if k in raw.columns})

    df = df.with_columns([
        pl.col("effect_allele").cast(pl.Utf8).str.to_uppercase().str.strip_chars(),
        pl.col("other_allele").cast(pl.Utf8).str.to_uppercase().str.strip_chars(),
        pl.col("chr").cast(pl.Utf8).str.replace(r"^chr", ""),
    ])

    # D026 columns from log10p (unsigned -log10P).
    df = df.with_columns(
        (-pl.col("log10p")).alias("log10p_signed"),
        pl.col("log10p").alias("neg_log10p"),
        (pl.col("log10p") > -LOG_FLOOR).alias("pval_underflow_flag"),
        (10.0 ** pl.when(pl.col("log10p") > -LOG_FLOOR)
                  .then(LOG_FLOOR)
                  .otherwise(-pl.col("log10p"))
        ).alias("pval_capped_for_tools"),
    )

    # n_controls = N_total - N_case
    df = df.with_columns(
        (pl.col("n").cast(pl.Float64) - pl.col("n_cases").cast(pl.Float64)).alias("n_controls"),
    )

    # Bookkeeping metadata.
    df = df.with_columns(
        pl.lit(outcome_label).alias("trait"),
        pl.lit("CVDKP_HERMES_2024_EUR").alias("source"),
        pl.lit("GRCh37").alias("build"),
        pl.lit("EUR").alias("ancestry"),
    )

    # QC step 1: drop missing beta/se/pval
    before = df.height
    df = df.filter(pl.col("beta").is_not_null()
                   & pl.col("se").is_not_null()
                   & pl.col("pval").is_not_null())
    dropped_missing = before - df.height
    _log(f"[{outcome_label}] dropped {dropped_missing:,} rows missing beta/se/pval")

    # QC step 2: invalid alleles
    valid_alleles = pl.col("effect_allele").is_in(list(VALID_ALLELES)) \
        & pl.col("other_allele").is_in(list(VALID_ALLELES)) \
        & (pl.col("effect_allele") != pl.col("other_allele"))
    before = df.height
    df = df.filter(valid_alleles)
    dropped_invalid_alleles = before - df.height
    _log(f"[{outcome_label}] dropped {dropped_invalid_alleles:,} invalid/multi-base allele rows")

    # QC step 3: MAF filter
    df = df.with_columns(
        pl.when(pl.col("eaf") <= 0.5).then(pl.col("eaf"))
          .otherwise(1.0 - pl.col("eaf"))
          .alias("maf")
    )
    before = df.height
    df = df.filter(pl.col("eaf").is_null() | (pl.col("maf") >= MAF_MIN))
    dropped_low_maf = before - df.height
    _log(f"[{outcome_label}] dropped {dropped_low_maf:,} rows with MAF < {MAF_MIN}")

    # QC step 4: palindromic — vectorised string concat replaces map_elements UDF
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
    _log(f"[{outcome_label}] dropped {dropped_palindromic:,} palindromic-ambiguous rows")

    n_out = df.height
    _log(f"[{outcome_label}] harmonized: {n_in:,} -> {n_out:,} ({n_out/max(1,n_in):.1%} kept)")

    # Write parquet
    keep_cols = [
        "trait", "source", "variant_id", "rsid", "chr", "pos", "build",
        "effect_allele", "other_allele", "beta", "se", "pval", "log10p",
        "log10p_signed", "neg_log10p", "pval_capped_for_tools",
        "pval_underflow_flag",
        "eaf", "maf", "n", "n_cases", "n_controls", "ancestry",
    ]
    keep_cols = [c for c in keep_cols if c in df.columns]
    df = df.select(keep_cols)
    parquet_path = PROCESSED / f"HF_2024_EUR_{outcome_label}.parquet"
    df.write_parquet(parquet_path, compression="zstd")
    _log(f"[{outcome_label}] wrote {parquet_path} ({parquet_path.stat().st_size/1e6:.1f} MB)")

    # Positive-control retention using b37-lifted coordinates.
    pc = positive_controls_b37[
        positive_controls_b37["gene_symbol"].isin(HF_POSITIVE_CONTROL_GENES)
        & positive_controls_b37["liftover_ok"]
    ].drop_duplicates(subset=["gene_symbol"]).copy()

    retention_rows = []
    for _, c in pc.iterrows():
        chrom = str(c["chr"])
        start = max(0, int(c["gene_start_b37"]) - WINDOW_KB * 1000)
        end = int(c["gene_end_b37"]) + WINDOW_KB * 1000
        sub = df.filter(
            (pl.col("chr") == chrom)
            & (pl.col("pos") >= start)
            & (pl.col("pos") <= end)
        )
        n_present = sub.height
        if n_var:
            min_neg_log10p = sub.select(pl.col("neg_log10p").max()).item()
        else:
            min_neg_log10p = None
        retention_rows.append({
            "outcome": outcome_label,
            "outcome_role": outcome_role,
            "gene_symbol": c["gene_symbol"],
            "chr": chrom,
            "gene_start_b37": int(c["gene_start_b37"]),
            "gene_end_b37": int(c["gene_end_b37"]),
            "window_kb": WINDOW_KB,
            "n_variants_in_cis_window": n_var,
            "max_neg_log10p_in_window": min_neg_log10p,
            "retained": n_present >= 1,
            "lead_variant_rsid_approx": c.get("lead_variant_rsid_approx", ""),
            "reference": c.get("reference", ""),
        })
    n_retained = sum(1 for r in retention_rows if r["retained"])
    n_total = len(retention_rows)
    _log(f"[{outcome_label}] retention: {n_retained}/{n_total} HF positive controls retained")

    return {
        "outcome": outcome_label,
        "outcome_role": outcome_role,
        "n_in": n_in,
        "n_out": n_out,
        "dropped_missing": dropped_missing,
        "dropped_invalid_alleles": dropped_invalid_alleles,
        "dropped_low_maf": dropped_low_maf,
        "dropped_palindromic": dropped_palindromic,
        "n_retained": n_retained,
        "n_pc_total": n_total,
        "retention_rows": retention_rows,
        "max_n_total": int(df.select(pl.col("n").max()).item()) if df.height else 0,
        "max_n_cases": int(df.select(pl.col("n_cases").max()).item()) if df.height else 0,
    }


def main() -> int:
    if not ZIP_PATH.exists():
        _log(f"FATAL: input not found: {ZIP_PATH}")
        return 2
    if not CHAIN_PATH.exists():
        _log(f"FATAL: liftover chain not found: {CHAIN_PATH}")
        return 2

    PROCESSED.mkdir(parents=True, exist_ok=True)
    HARMONIZATION_REPORT.parent.mkdir(parents=True, exist_ok=True)
    EXTRACT_DIR.mkdir(parents=True, exist_ok=True)

    # Extract zip — the CVDKP bundle nests an inner zip; unwrap both layers.
    _log(f"opening {ZIP_PATH.name}")
    with zipfile.ZipFwith(ZIP_PATH) as zf:
        members = zf.namelist()
        _log(f"outer zip contains {len(members)} entries")
        zf.extractall(EXTRACT_DIR)

    inner_zips = list(EXTRACT_DIR.rglob("*.zip"))
    for inner in inner_zips:
        _log(f"unwrapping inner zip {inner.relative_to(EXTRACT_DIR)}")
        with zipfile.ZipFwith(inner) as zf:
            zf.extractall(inner.parent)
        # Keep inner zip for reproducibility but mark it processed via mtime
        # (no rename needed; downstream globbing uses .tsv.gz extension only).

    registry = load_column_maps(ROOT / "config/column_maps.yml")

    pc = load_positive_controls(ROOT / "config/positive_controls.tsv")
    pc_b37 = liftover_positive_controls(pc, CHAIN_PATH)
    n_lift_fawithd = int((~pc_b37["liftover_ok"]).sum())
    if n_lift_fawithd:
        _log(f"WARN: {n_lift_fawithd} positive controls fawithd liftover; check chain coverage")

    summaries = []
    all_retention_rows: list[dict] = []
    for pheno_code, (outcome_label, outcome_role) in PHENOTYPE_MAP.items():
        # The inner zip may nest phenotype dirs at any depth; locate the file.
        candidates = list(EXTRACT_DIR.rglob(f"FORMAT-METAL_{pheno_code}_EUR.tsv.gz"))
        if not candidates:
            _log(f"WARN: missing FORMAT-METAL_{pheno_code}_EUR.tsv.gz under {EXTRACT_DIR}; skipping")
            continue
        tsv_path = candidates[0]
        _log(f"[{outcome_label}] resolved {tsv_path}")
        s = harmonize_phenotype(tsv_path, outcome_label, outcome_role, registry, pc_b37)
        summaries.append(s)
        all_retention_rows.extend(s["retention_rows"])

    # Reports
    if summaries:
        rep = pl.DataFrame([{
            "outcome": s["outcome"],
            "outcome_role": s["outcome_role"],
            "input_rows": s["n_in"],
            "dropped_missing_beta_se_pval": s["dropped_missing"],
            "dropped_invalid_alleles": s["dropped_invalid_alleles"],
            "dropped_low_maf": s["dropped_low_maf"],
            "dropped_palindromic_unresolved": s["dropped_palindromic"],
            "output_rows": s["n_out"],
            "max_n_total": s["max_n_total"],
            "max_n_cases": s["max_n_cases"],
            "positive_controls_retained": s["n_retained"],
            "positive_controls_total": s["n_pc_total"],
        } for s in summaries])
        rep.write_csv(HARMONIZATION_REPORT, separator="\t")
        _log(f"wrote {HARMONIZATION_REPORT}")

    if all_retention_rows:
        ret = pl.DataFrame(all_retention_rows)
        ret.write_csv(RETENTION_REPORT, separator="\t")
        _log(f"wrote {RETENTION_REPORT}")

    # Pass / fail summary
    all_ok = True
    for s in summaries:
        if s["outcome_role"] == "primary" and s["n_retained"] < s["n_pc_total"]:
            _log(f"WARN: primary outcome {s['outcome']} retention {s['n_retained']}/{s['n_pc_total']}")
            # Allow partial retention but flag for review:
            if s["n_retained"] == 0:
                all_ok = False
        if s["outcome_role"] == "exploratory":
            _log(f"INFO: exploratory outcome {s['outcome']} retention {s['n_retained']}/{s['n_pc_total']}")

    # Disk hygiene: drop the unzipped tree once parquets are emitted.
    # Original ZIP_PATH is kept for reproducibility (re-extract on demand).
    import shutil
    if all_ok and EXTRACT_DIR.exists():
        try:
            shutil.rmtree(EXTRACT_DIR)
            _log(f"cleaned up {EXTRACT_DIR}")
        except OSError as e:
            _log(f"WARN: could not clean {EXTRACT_DIR}: {e}")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
