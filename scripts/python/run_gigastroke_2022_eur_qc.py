#!/usr/bin/env python3
"""GIGASTROKE 2022 EUR (Mishra et al., Nat Genet 2022, PMID 36180795).

D029 primary stroke source — harmonised summary statistics from GWAS Catalog
(GCST90104540=AIS, GCST90104541=CES, GCST90104542=LAS, GCST90104543=SVS).
All four phenotypes are GRCh38 — no liftover needed for retention.

D028 source validation: meta.yaml verifies trait + ancestry + build before
the harmonization step runs.

Outputs:
    results/qc/gigastroke_2022_eur_source_validation.tsv
    results/qc/gigastroke_2022_eur_harmonization_report.tsv
    results/qc/gigastroke_2022_eur_positive_control_retention.tsv
    data/processed/sumstats/STROKE_AIS_EUR.parquet
    data/processed/sumstats/STROKE_CES_EUR.parquet
    data/processed/sumstats/STROKE_LAS_EUR.parquet
    data/processed/sumstats/STROKE_SVS_EUR.parquet
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import polars as pl
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from afshf.column_maps import load_column_maps  # noqa: E402
from afshf.positive_controls import load_positive_controls  # noqa: E402

RAW_DIR = ROOT / "data/raw/outcomes/stroke_eur/gigastroke"
PROCESSED = ROOT / "data/processed/sumstats"
SOURCE_REPORT = ROOT / "results/qc/gigastroke_2022_eur_source_validation.tsv"
HARMONIZATION_REPORT = ROOT / "results/qc/gigastroke_2022_eur_harmonization_report.tsv"
RETENTION_REPORT = ROOT / "results/qc/gigastroke_2022_eur_positive_control_retention.tsv"

# GCST -> (outcome label, role, parquet name, expected case count)
PHENOTYPES = {
    "GCST90104540": ("ischemic_stroke",        "primary",     "STROKE_AIS_EUR.parquet",  62100),
    "GCST90104541": ("cardioembolic_stroke",   "primary",     "STROKE_CES_EUR.parquet",  10804),
    "GCST90104542": ("large_artery_stroke",    "primary",     "STROKE_LAS_EUR.parquet",   6399),
    "GCST90104543": ("small_vessel_stroke",    "primary",     "STROKE_SVS_EUR.parquet",   6811),
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


def validate_source(meta_path: Path,
                    expected_trait: str,
                    expected_n_case: int) -> dict:
    """D028 source validation against meta.yaml."""
    with open(meta_path) as f:
        m = yaml.safe_load(f)
    traits = m.get("trait_description") or []
    ancestries = []
    n_total = 0
    for s in m.get("samples", []):
        ancestries.extend(s.get("sample_ancestry_category", []))
        n_total += int(s.get("sample_size", 0))
    build = m.get("genome_assembly", "")
    is_harmonised = m.get("is_harmonised", False)
    pmid = m.get("gwas_catalog_api", "").rsplit("/", 1)[-1]
    md5 = m.get("data_file_md5sum", "")
    trait_match = expected_trait.lower() in [t.lower() for t in traits]
    ancestry_match = "European" in ancestries
    build_match = build == "GRCh38"

    return {
        "trait_listed": traits,
        "expected_trait": expected_trait,
        "trait_match": trait_match,
        "ancestry": ancestries,
        "ancestry_match": ancestry_match,
        "build": build,
        "build_match": build_match,
        "n_total_meta": n_total,
        "is_harmonised": bool(is_harmonised),
        "md5_meta": md5,
        "validation_pass": trait_match and ancestry_match and build_match,
    }


def harmonize_gigastroke(tsv_path: Path,
                         outcome_label: str,
                         outcome_role: str,
                         parquet_name: str,
                         registry,
                         positive_controls) -> dict:
    cmap = registry.get("GIGASTROKE_GWAS_Catalog")
    _log(f"[{outcome_label}] reading {tsv_path.name}")

    schema_overrides = {
        "chromosome": pl.Utf8,
        "base_pair_location": pl.Int64,
        "effect_allele": pl.Utf8,
        "other_allele": pl.Utf8,
        "beta": pl.Float64,
        "standard_error": pl.Float64,
        "effect_allele_frequency": pl.Float64,
        "p_value": pl.Float64,
        "odds_ratio": pl.Float64,
        "rsid": pl.Utf8,
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
    for std_field, alt_list in cmap.column_aliases.items():
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
        pl.col("chr").cast(pl.Utf8).str.replace(r"^chr", ""),
    ])

    # D026 underflow-safe columns from p_value (Mishra ships P directly).
    df = df.with_columns(
        pl.when(pl.col("pval") <= 0).then(LOG_FLOOR)
          .otherwise(pl.col("pval").log10()).alias("log10p_signed"),
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
        pl.lit("GIGASTROKE_2022_EUR").alias("source"),
        pl.lit("GRCh38").alias("build"),
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
    is_pal = pl.col("_pair").is_in(PALINDROMIC_PAIRS_STR)
    keep_palindromic = (~is_pal | (pl.col("eaf").is_not_null()
                                   & ((pl.col("eaf") - 0.5).abs() >= PALINDROMIC_AMBIG_WINDOW)))
    before = df.height
    df = df.filter(keep_palindromic).drop("_pair")
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
    parquet_path = PROCESSED / parquet_name
    df.select(keep_cols).write_parquet(parquet_path, compression="zstd")
    _log(f"[{outcome_label}] wrote {parquet_path} ({parquet_path.stat().st_size/1e6:.1f} MB)")

    # Retention check — positive_controls.tsv is GRCh38, GIGASTROKE is GRCh38, no liftover.
    pc = positive_controls[
        positive_controls["gene_symbol"].isin(STROKE_POSITIVE_CONTROL_GENES)
    ].drop_duplicates(subset=["gene_symbol"])

    rows = []
    for _, c in pc.iterrows():
        chrom = str(c["chr"])
        start = max(0, int(c["gene_start"]) - WINDOW_KB * 1000)
        end = int(c["gene_end"]) + WINDOW_KB * 1000
        sub = df.filter((pl.col("chr") == chrom)
                        & (pl.col("pos") >= start)
                        & (pl.col("pos") <= end))
        n_present = sub.height
        max_neg = sub.select(pl.col("neg_log10p").max()).item() if n_present else None
        rows.append({
            "outcome": outcome_label, "outcome_role": outcome_role,
            "gene_symbol": c["gene_symbol"], "chr": chrom,
            "gene_start": int(c["gene_start"]), "gene_end": int(c["gene_end"]),
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
        _log(f"FATAL: {RAW_DIR} missing")
        return 2
    PROCESSED.mkdir(parents=True, exist_ok=True)
    SOURCE_REPORT.parent.mkdir(parents=True, exist_ok=True)

    # D028 source validation per meta.yaml
    val_rows = []
    for accession, (outcome, role, parquet, expected_n) in PHENOTYPES.items():
        meta_path = RAW_DIR / f"{accession}.h.tsv.gz-meta.yaml"
        if not meta_path.exists():
            _log(f"WARN: missing {meta_path}")
            continue
        v = validate_source(meta_path, expected_trait=outcome.replace("_", " "),
                            expected_n_case=expected_n)
        v["accession"] = accession
        v["outcome"] = outcome
        v["outcome_role"] = role
        val_rows.append(v)
        if not v["validation_pass"]:
            _log(f"D028 FAIL: {accession} {outcome} validation_pass=False -> {v}")

    pl.DataFrame([{
        "accession": v["accession"], "outcome": v["outcome"], "outcome_role": v["outcome_role"],
        "trait_listed": "|".join(v["trait_listed"]), "expected_trait": v["expected_trait"],
        "trait_match": v["trait_match"], "ancestry": "|".join(v["ancestry"]),
        "ancestry_match": v["ancestry_match"], "build": v["build"],
        "build_match": v["build_match"], "n_total_meta": v["n_total_meta"],
        "is_harmonised": v["is_harmonised"], "md5_meta": v["md5_meta"],
        "validation_pass": v["validation_pass"],
    } for v in val_rows]).write_csv(SOURCE_REPORT, separator="\t")
    _log(f"wrote {SOURCE_REPORT}")

    if not all(v["validation_pass"] for v in val_rows):
        _log("FATAL: D028 source validation fawithd for at least one phenotype.")
        return 1

    registry = load_column_maps(ROOT / "config/column_maps.yml")
    pc = load_positive_controls(ROOT / "config/positive_controls.tsv")

    summaries = []
    all_retention_rows: list[dict] = []
    for accession, (outcome, role, parquet, _) in PHENOTYPES.items():
        tsv_path = RAW_DIR / f"{accession}.h.tsv.gz"
        if not tsv_path.exists():
            _log(f"WARN: missing {tsv_path}; skipping")
            continue
        s = harmonize_gigastroke(tsv_path, outcome, role, parquet, registry, pc)
        summaries.append(s)
        all_retention_rows.extend(s["retention_rows"])

    if summaries:
        pl.DataFrame([{
            "outcome": s["outcome"], "outcome_role": s["outcome_role"],
            "input_rows": s["n_in"],
            "dropped_missing_beta_se_pval": s["dropped_missing"],
            "dropped_invalid_alleles": s["dropped_invalid"],
            "dropped_low_maf": s["dropped_low_maf"],
            "dropped_palindromic_unresolved": s["dropped_palindromic"],
            "output_rows": s["n_out"],
            "positive_controls_retained": s["n_retained"],
            "positive_controls_total": s["n_pc_total"],
        } for s in summaries]).write_csv(HARMONIZATION_REPORT, separator="\t")
        _log(f"wrote {HARMONIZATION_REPORT}")

    if all_retention_rows:
        pl.DataFrame(all_retention_rows).write_csv(RETENTION_REPORT, separator="\t")
        _log(f"wrote {RETENTION_REPORT}")

    n_total = sum(s["n_pc_total"] for s in summaries)
    n_retained = sum(s["n_retained"] for s in summaries)
    _log(f"GIGASTROKE retention: {n_retained}/{n_total}")
    return 0 if n_retained == n_total else 1


if __name__ == "__main__":
    sys.exit(main())
