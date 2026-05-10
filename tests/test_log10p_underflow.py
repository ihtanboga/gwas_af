"""D026 — underflow-safe significance columns.

Protocol rule (2026-05-09):
    log_p = -1140.76 must NOT become missing.
    pval_capped_for_tools = 1e-300
    neg_log10p = 1140.76
    threshold decisions use neg_log10p, not pval=0.
"""

from math import isclose, log10

import numpy as np
import pandas as pd

from afshf.column_maps import (
    PVAL_CAPPED_FLOOR,
    derive_log10p_columns,
)


# ---------------------------------------------------------------------------
# Protocol regression (verbatim)
# ---------------------------------------------------------------------------


def test_signed_log_p_extreme_underflow_PITX2_case():
    """PITX2 in Roselli 2025: log_p ≈ -1140.76. Underflow must be flagged
    and the row must NOT be missing."""
    df = pd.DataFrame({"log_p": [-1140.76, -10.0, -0.30]})
    out = derive_log10p_columns(df, source_signed_col="log_p",
                                source_unsigned_col=None)

    # Row 0 — extreme underflow case.
    assert pd.notna(out.loc[0, "log10p_signed"])
    assert isclose(out.loc[0, "log10p_signed"], -1140.76)
    assert isclose(out.loc[0, "neg_log10p"], 1140.76)
    assert isclose(out.loc[0, "pval_capped_for_tools"], 1e-300)
    assert bool(out.loc[0, "pval_underflow_flag"]) is True

    # Row 1 — strong but not underflow.
    assert isclose(out.loc[1, "neg_log10p"], 10.0)
    assert isclose(out.loc[1, "pval_capped_for_tools"], 1e-10, rel_tol=1e-9)
    assert bool(out.loc[1, "pval_underflow_flag"]) is False

    # Row 2 — non-significant.
    assert isclose(out.loc[2, "neg_log10p"], 0.30)
    assert isclose(out.loc[2, "pval_capped_for_tools"], 10 ** -0.30,
                   rel_tol=1e-9)
    assert bool(out.loc[2, "pval_underflow_flag"]) is False


def test_unsigned_log10p_input_ukbppp_style():
    """UKB-PPP style: log10p column is unsigned -log10(P)."""
    df = pd.DataFrame({"log10p": [301.0, 50.0, 0.30]})
    out = derive_log10p_columns(df, source_signed_col=None,
                                source_unsigned_col="log10p")

    # 301 -> p ≈ 1e-301 -> below 1e-300 floor -> capped to 1e-300, flag True.
    assert isclose(out.loc[0, "log10p_signed"], -301.0)
    assert isclose(out.loc[0, "neg_log10p"], 301.0)
    assert isclose(out.loc[0, "pval_capped_for_tools"], 1e-300)
    assert bool(out.loc[0, "pval_underflow_flag"]) is True

    # 50 -> p = 1e-50, no underflow.
    assert isclose(out.loc[1, "pval_capped_for_tools"], 1e-50, rel_tol=1e-9)
    assert bool(out.loc[1, "pval_underflow_flag"]) is False


def test_floor_constant_is_1e_minus_300():
    assert PVAL_CAPPED_FLOOR == 1e-300


def test_natural_log_signed_input():
    """If a source ships ln(P), base_is_natural=True converts to log10."""
    # ln(1e-10) ≈ -23.026
    df = pd.DataFrame({"log_p": [-23.025850929940457]})
    out = derive_log10p_columns(df, source_signed_col="log_p",
                                source_unsigned_col=None,
                                base_is_natural=True)
    # log10p_signed should be ≈ -10.
    assert isclose(out.loc[0, "log10p_signed"], -10.0, abs_tol=1e-9)
    assert isclose(out.loc[0, "neg_log10p"], 10.0, abs_tol=1e-9)
    assert isclose(out.loc[0, "pval_capped_for_tools"], 1e-10, rel_tol=1e-9)


def test_no_columns_added_when_no_source_present():
    df = pd.DataFrame({"foo": [1, 2, 3]})
    out = derive_log10p_columns(df, source_signed_col="log_p",
                                source_unsigned_col="log10p")
    assert "log10p_signed" not in out.columns
    assert "neg_log10p" not in out.columns


# ---------------------------------------------------------------------------
# Protocol rule 4: significance threshold decisions use neg_log10p
# ---------------------------------------------------------------------------


def test_threshold_decisions_use_neg_log10p_not_pval_zero():
    """When pval underflows to 0, neg_log10p still gives a correct decision."""
    df = pd.DataFrame({
        "log_p": [-1140.76, -10.0, -0.30],
    })
    out = derive_log10p_columns(df, source_signed_col="log_p")

    threshold_p = 0.05
    threshold_neg_log10p = -log10(threshold_p)         # ≈ 1.301

    # Use neg_log10p path.
    significant_neg = out["neg_log10p"] > threshold_neg_log10p

    assert bool(significant_neg.iloc[0]) is True       # PITX2: 1140 > 1.301
    assert bool(significant_neg.iloc[1]) is True       # 10 > 1.301
    assert bool(significant_neg.iloc[2]) is False      # 0.3 < 1.301

    # Pval-zero path would still be technically OK here (0 < 0.05 is True for
    # significant), but downstream R/coloc tools cannot accept pval=0. The
    # capped column always works in their place.
    assert (out["pval_capped_for_tools"] >= PVAL_CAPPED_FLOOR).all()


# ---------------------------------------------------------------------------
# Pipeline-level: apply_column_map populates D026 columns
# ---------------------------------------------------------------------------


def test_apply_column_map_for_cvdkp_af_yields_d026_columns():
    """End-to-end: a CVDKP_AF-like row gets the four new columns."""
    from pathlib import Path

    from afshf.column_maps import apply_column_map, load_column_maps

    registry = load_column_maps(
        Path(__file__).resolve().parents[1] / "config" / "column_maps.yml"
    )
    cmap = registry.get("CVDKP_AF")
    raw = pd.DataFrame([{
        "MarkerName": "4:111720761:T:C", "rsid": "rs2200733",
        "chr": "4", "position_b38": 111720761,
        "Allele1": "T", "Allele2": "C",
        "Freq1": 0.12, "Effect": 0.143, "StdErr": 0.012,
        "log(P)": -1140.76, "n_events": 123790, "n_total": 1319851,
        "mean_impQual": 0.99,
    }])
    out = apply_column_map(raw, cmap, trait="atrial_fibrillation",
                           source="CVDKP_AF")
    assert {"log10p_signed", "neg_log10p", "pval_capped_for_tools",
            "pval_underflow_flag"}.issubset(out.columns)
    row = out.iloc[0]
    assert isclose(row["neg_log10p"], 1140.76)
    assert isclose(row["pval_capped_for_tools"], 1e-300)
    assert bool(row["pval_underflow_flag"]) is True


def test_classify_target_uses_neg_log10p_when_pval_underflows():
    """Protocol rule 6: threshold decisions use neg_log10p when pval=0."""
    from afshf.intervention_direction import TargetOutcomeEffect
    from afshf.scoring import classify_target

    # PITX2-like: AF MR has p underflowed to 0, but neg_log10p ≈ 1140.
    af = TargetOutcomeEffect(
        target_id="PITX2", outcome="AF",
        beta=-0.20, se=0.01, p=0.0, fdr=0.0, pph4=0.95,
        action_direction=1, theta=-0.20, or_intervention=0.819,
        benefit_flag=False, harm_flag=False,
        neg_log10p=1140.0,
    )
    # Add an HF harm signal (normal, no underflow).
    hf = TargetOutcomeEffect(
        target_id="PITX2", outcome="HF_overall",
        beta=0.10, se=0.02, p=1e-6, fdr=1e-6, pph4=0.85,
        action_direction=1, theta=0.10, or_intervention=1.105,
        benefit_flag=False, harm_flag=False,
        neg_log10p=6.0,
    )
    label = classify_target(af, [hf])
    # AF tier strong (neg_log10p>>1.301 + PPH4=0.95 -> strong) and HF harm
    # significant -> RED-1.
    assert label == "RED-1_AF_BENEFIT_HF_HARM"


def test_classify_target_falls_back_to_p_when_neg_log10p_missing():
    """Backward compatibility: existing TargetOutcomeEffect with no
    neg_log10p still uses raw p/fdr."""
    from afshf.intervention_direction import TargetOutcomeEffect
    from afshf.scoring import classify_target

    af = TargetOutcomeEffect(
        target_id="X", outcome="AF",
        beta=-0.20, se=0.01, p=1e-9, fdr=1e-9, pph4=0.95,
        action_direction=1, theta=-0.20, or_intervention=0.819,
        benefit_flag=False, harm_flag=False,
        # neg_log10p left at default None
    )
    hf = TargetOutcomeEffect(
        target_id="X", outcome="HF_overall",
        beta=0.10, se=0.02, p=1e-6, fdr=1e-6, pph4=0.85,
        action_direction=1, theta=0.10, or_intervention=1.105,
        benefit_flag=False, harm_flag=False,
    )
    assert classify_target(af, [hf]) == "RED-1_AF_BENEFIT_HF_HARM"


def test_apply_column_map_for_ukbppp_yields_d026_columns():
    """UKB-PPP style: log10p (unsigned) -> D026 columns."""
    from pathlib import Path

    from afshf.column_maps import apply_column_map, load_column_maps

    registry = load_column_maps(
        Path(__file__).resolve().parents[1] / "config" / "column_maps.yml"
    )
    cmap = registry.get("UKBPPP")
    raw = pd.DataFrame([{
        "CHROM": "1", "GENPOS": 55039974, "ID": "rs11591147",
        "ALLELE0": "G", "ALLELE1": "T", "A1FREQ": 0.015,
        "BETA": -1.234, "SE": 0.045, "LOG10P": 100.5, "N": 54219,
        "INFO": 0.99,
    }])
    out = apply_column_map(raw, cmap, source="UKBPPP", trait="PCSK9_pQTL")
    assert "neg_log10p" in out.columns
    assert isclose(out.iloc[0]["neg_log10p"], 100.5)
    # 10**(-100.5) is well above the 1e-300 floor.
    assert bool(out.iloc[0]["pval_underflow_flag"]) is False
