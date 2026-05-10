"""D005 REVISED — AF action-direction tier system.

Tiers:
    strong   : FDR_AF < 0.05 AND PPH4_AF >= 0.80
    moderate : p_AF   < 0.05 AND PPH4_AF >= 0.50  (FDR not required)
    weak     : p_AF   < 0.05 AND PPH4_AF <  0.50
    none     : p_AF   >= 0.05 OR beta missing
"""

from afshf.intervention_direction import (
    PRIMARY_CLASS_ELIGIBLE_TIERS,
    af_signal_strength,
)


def test_strong_requires_fdr_and_strong_coloc():
    assert af_signal_strength(p_target_to_AF=1e-9, fdr_target_to_AF=1e-7,
                              pph4_AF=0.85) == "strong"


def test_moderate_when_fdr_misses_but_coloc_moderate():
    # FDR >= 0.05 -> not strong; p < 0.05 + PPH4 in [0.50, 0.80) -> moderate
    assert af_signal_strength(p_target_to_AF=1e-3, fdr_target_to_AF=0.10,
                              pph4_AF=0.65) == "moderate"


def test_moderate_when_pph4_strong_but_fdr_missing():
    # FDR not provided; p < 0.05 + PPH4 >= 0.50 -> moderate
    assert af_signal_strength(p_target_to_AF=1e-4, fdr_target_to_AF=None,
                              pph4_AF=0.85) == "moderate"


def test_weak_when_coloc_below_moderate():
    assert af_signal_strength(p_target_to_AF=1e-3, fdr_target_to_AF=0.001,
                              pph4_AF=0.20) == "weak"


def test_none_when_p_above_threshold():
    assert af_signal_strength(p_target_to_AF=0.20, fdr_target_to_AF=0.20,
                              pph4_AF=0.95) == "none"


def test_none_when_p_missing():
    assert af_signal_strength(p_target_to_AF=None, fdr_target_to_AF=None,
                              pph4_AF=0.95) == "none"


def test_primary_class_eligibility():
    assert "strong" in PRIMARY_CLASS_ELIGIBLE_TIERS
    assert "moderate" in PRIMARY_CLASS_ELIGIBLE_TIERS
    assert "weak" not in PRIMARY_CLASS_ELIGIBLE_TIERS
    assert "none" not in PRIMARY_CLASS_ELIGIBLE_TIERS
