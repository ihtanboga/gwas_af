"""Spec section 18 — scoring formulas (Section 10.7)."""

from math import isclose

from afshf.intervention_direction import TargetOutcomeEffect
from afshf.scoring import (
    aggregate_target_score,
    af_benefit_score,
    coloc_weight,
    harm_score,
    replication_weight,
    shared_benefit_score,
    therapeutic_balance,
)


SEVERITY = {
    "AF": 1.0,
    "HF_overall": 1.2,
    "HF_ni_HFrEF": 1.3,
    "HF_ni_HFpEF": 1.2,
    "ischemic_stroke": 1.5,
    "cardioembolic_stroke": 1.4,
    "large_artery_stroke": 1.3,
    "small_vessel_stroke": 1.3,
}


def _eff(outcome, theta, p, pph4):
    # Beta = theta when action_direction = +1, else -theta. We do not need to
    # round-trip beta in the score tests.
    return TargetOutcomeEffect(
        target_id="X", outcome=outcome,
        beta=theta, se=0.05, p=p, fdr=p, pph4=pph4,
        action_direction=1, theta=theta, or_intervention=1.0,
        benefit_flag=theta < 0 and pph4 >= 0.5 and p < 0.05,
        harm_flag=theta > 0 and pph4 >= 0.5 and p < 0.05,
    )


def test_coloc_weight_buckets():
    assert coloc_weight(0.95) == 1.0
    assert coloc_weight(0.65) == 0.6
    assert coloc_weight(0.10) == 0.2
    assert coloc_weight(None) == 0.2


def test_replication_weight_buckets():
    assert replication_weight("replicated") == 1.0
    assert replication_weight("partial") == 0.6
    assert replication_weight("none") == 0.3


def test_af_benefit_score_zero_when_no_protection():
    eff = _eff("AF", theta=0.10, p=1e-3, pph4=0.9)  # harm direction
    assert af_benefit_score(eff) == 0.0


def test_af_benefit_score_positive_when_protective():
    eff = _eff("AF", theta=-0.10, p=1e-5, pph4=0.9)
    expected = 5.0 * 1.0 * 0.3  # -log10(1e-5) * coloc_strong * replication_none
    assert isclose(af_benefit_score(eff), expected)


def test_harm_score_uses_severity_weight():
    eff = _eff("HF_overall", theta=0.20, p=1e-4, pph4=0.85)
    expected = 4.0 * 1.0 * 0.3 * SEVERITY["HF_overall"]
    assert isclose(harm_score(eff, SEVERITY["HF_overall"]), expected)


def test_shared_benefit_score_zero_when_outcome_harm():
    eff = _eff("HF_overall", theta=0.20, p=1e-4, pph4=0.9)
    assert shared_benefit_score(eff) == 0.0


def test_therapeutic_balance_penalises_harm():
    af = _eff("AF", theta=-0.10, p=1e-5, pph4=0.9)
    hf = _eff("HF_overall", theta=0.10, p=1e-3, pph4=0.85)

    af_b = af_benefit_score(af)
    hf_h = harm_score(hf, SEVERITY["HF_overall"])
    bal = therapeutic_balance(af_b, 0.0, hf_h, lambda_harm=1.5)
    assert isclose(bal, af_b - 1.5 * hf_h)


def test_aggregate_target_score_components_present():
    af = _eff("AF", theta=-0.20, p=1e-8, pph4=0.95)
    outcomes = [
        _eff("HF_overall", theta=0.10, p=1e-3, pph4=0.85),
        _eff("cardioembolic_stroke", theta=-0.15, p=1e-4, pph4=0.85),
    ]
    score = aggregate_target_score(
        target_id="X",
        af_effect=af,
        outcome_effects=outcomes,
        severity_weights=SEVERITY,
    )
    assert score.af_benefit > 0.0
    assert score.shared_benefit_total > 0.0
    assert score.harm_total > 0.0
    assert "AF_benefit" in score.components
    assert "HF_overall_harm" in score.components
    assert "cardioembolic_stroke_benefit" in score.components
