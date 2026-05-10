"""End-to-end smoke test: column_map -> harmonize -> positive control retention.

Mirrors the user-requested checklist:
  1. Mini sumstats with realistic source headers can be ingested.
  2. Palindromic SNP filter behaves correctly.
  3. OR -> beta conversion is applied where required.
  4. AF positive-control loci (PITX2, ZFHX3, KCNN3, ...) survive harmonize().
"""

from math import isclose, log
from pathlib import Path

import pandas as pd
import pytest

from afshf.column_maps import apply_column_map, load_column_maps
from afshf.harmonize import STANDARD_COLUMNS, harmonize_sumstats
from afshf.positive_controls import (
    check_positive_control_retention,
    load_positive_controls,
    summarise_retention,
)


CONFIG_DIR = Path(__file__).resolve().parents[1] / "config"
FIXTURES = Path(__file__).resolve().parent / "synthetic_data"


@pytest.fixture(scope="module")
def registry():
    return load_column_maps(CONFIG_DIR / "column_maps.yml")


@pytest.fixture(scope="module")
def positive_controls():
    return load_positive_controls(CONFIG_DIR / "positive_controls.tsv")


# ---------------------------------------------------------------------------
# Positive control registry sanity
# ---------------------------------------------------------------------------


def test_positive_controls_has_expected_genes(positive_controls):
    af_genes = set(positive_controls.loc[positive_controls["trait"] == "AF",
                                         "gene_symbol"])
    assert {"PITX2", "ZFHX3", "KCNN3", "PRRX1", "CAV1", "TBX5"}.issubset(af_genes)

    hf_genes = set(positive_controls.loc[positive_controls["trait"] == "HF_overall",
                                         "gene_symbol"])
    assert {"BAG3", "TTN", "NPPA", "NPPB", "FLNC"}.issubset(hf_genes)

    stroke_genes = set(positive_controls.loc[
        positive_controls["trait"].isin(
            ["large_artery_stroke", "cardioembolic_stroke"]),
        "gene_symbol"])
    assert {"HDAC9", "PHACTR1", "ABO", "PITX2"}.issubset(stroke_genes)


# ---------------------------------------------------------------------------
# CVDKP_AF end-to-end smoke
# ---------------------------------------------------------------------------


def test_cvdkp_af_smoke_pipeline(registry, positive_controls):
    raw = pd.read_csv(FIXTURES / "cvdkp_af_mini.tsv", sep="\t")
    standardised = apply_column_map(
        raw, registry.get("CVDKP_AF"),
        trait="atrial_fibrillation", source="CVDKP_AF",
    )

    # Build all 18 standard columns the harmonizer expects.
    for col in STANDARD_COLUMNS:
        if col not in standardised.columns:
            standardised[col] = pd.NA

    cleaned, report = harmonize_sumstats(standardised)

    # The truly missing-beta row should be dropped.
    assert "rs99999992" not in set(cleaned["rsid"].dropna())
    assert report.dropped_missing_betas >= 1

    # Positive-control retention: every AF positive locus has at least one variant.
    result = check_positive_control_retention(
        cleaned, positive_controls, trait="AF", window_kb=1000,
    )
    assert result.total_loci > 0, "AF positive controls registry is empty"
    assert result.all_retained, f"missing AF loci: {result.missing}"

    # Cardioembolic stroke positives (PITX2, ZFHX3, ABO) live on the same file
    # because the AF lead variants colocate with those loci.
    ces_result = check_positive_control_retention(
        cleaned, positive_controls, trait="cardioembolic_stroke", window_kb=1000,
    )
    assert ces_result.retained >= 2  # PITX2 + ZFHX3 at minimum


def test_summarise_retention_table_shape(registry, positive_controls):
    raw = pd.read_csv(FIXTURES / "cvdkp_af_mini.tsv", sep="\t")
    standardised = apply_column_map(
        raw, registry.get("CVDKP_AF"),
        trait="atrial_fibrillation", source="CVDKP_AF",
    )
    for col in STANDARD_COLUMNS:
        if col not in standardised.columns:
            standardised[col] = pd.NA
    cleaned, _ = harmonize_sumstats(standardised)

    results = [
        check_positive_control_retention(cleaned, positive_controls, trait=t)
        for t in ["AF", "HF_overall", "cardioembolic_stroke"]
    ]
    table = summarise_retention(results)
    assert {"trait", "total_loci", "retained", "all_retained"}.issubset(table.columns)
    assert len(table) == 3


# ---------------------------------------------------------------------------
# GWAS Catalog OR-only smoke (verifies OR -> beta wiring)
# ---------------------------------------------------------------------------


def test_gwas_catalog_or_to_beta_smoke(registry, positive_controls):
    raw = pd.read_csv(FIXTURES / "gwas_catalog_or_mini.tsv", sep="\t")
    standardised = apply_column_map(
        raw, registry.get("GWAS_Catalog"),
        trait="atrial_fibrillation", source="GWAS_Catalog",
    )
    for col in STANDARD_COLUMNS:
        if col not in standardised.columns:
            standardised[col] = pd.NA

    # Sanity: rs2200733 originally OR=1.154 -> beta = log(1.154)
    rs = standardised.set_index("rsid")
    assert isclose(rs.loc["rs2200733", "beta"], log(1.154), abs_tol=1e-9)

    cleaned, _ = harmonize_sumstats(standardised)
    result = check_positive_control_retention(
        cleaned, positive_controls, trait="AF", window_kb=1000,
    )
    # Mini fixture only contains 3 AF lead variants (PITX2, ZFHX3, KCNN3) so
    # only those three loci can possibly retain.
    assert "PITX2" not in result.missing
    assert "ZFHX3" not in result.missing
    assert "KCNN3" not in result.missing


# ---------------------------------------------------------------------------
# Palindromic filter behaves under the column map (uses harmonize.py rules)
# ---------------------------------------------------------------------------


def test_palindromic_filter_in_smoke():
    """Inject a palindromic ambiguous variant and confirm it is dropped."""
    df = pd.DataFrame([
        {
            "trait": "AF", "source": "synthetic", "variant_id": "x", "rsid": "rs_pal_amb",
            "chr": "1", "pos": 1234, "build": "GRCh38",
            "effect_allele": "A", "other_allele": "T",
            "beta": 0.05, "se": 0.01, "or": None, "pval": 1e-3, "eaf": 0.49,
            "n": 1000, "n_cases": 100, "n_controls": 900, "ancestry": "EUR",
        },
        {
            "trait": "AF", "source": "synthetic", "variant_id": "y", "rsid": "rs_pal_keep",
            "chr": "1", "pos": 1235, "build": "GRCh38",
            "effect_allele": "A", "other_allele": "T",
            "beta": 0.05, "se": 0.01, "or": None, "pval": 1e-3, "eaf": 0.12,
            "n": 1000, "n_cases": 100, "n_controls": 900, "ancestry": "EUR",
        },
    ])
    cleaned, report = harmonize_sumstats(df)
    rsid_kept = set(cleaned["rsid"])
    assert "rs_pal_keep" in rsid_kept
    assert "rs_pal_amb" not in rsid_kept
    assert report.dropped_palindromic_unresolved == 1
