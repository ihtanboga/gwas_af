"""Spec section 18 — Tests 4 and 5 (target classification)."""

from afshf.intervention_direction import TargetOutcomeEffect
from afshf.scoring import classify_target


def _eff(outcome, theta, p, pph4, fdr=None):
    return TargetOutcomeEffect(
        target_id="X", outcome=outcome,
        beta=theta, se=0.05, p=p,
        fdr=fdr if fdr is not None else p,
        pph4=pph4,
        action_direction=1,
        theta=theta,
        or_intervention=1.0,
        benefit_flag=False, harm_flag=False,
    )


def test_red1_af_benefit_hf_harm():
    """Spec Test 4: theta_AF<0 + theta_HF>0 with strong stats and coloc."""
    af = _eff("AF", theta=-0.20, p=1e-9, pph4=0.85, fdr=1e-9)
    hf = _eff("HF_overall", theta=0.15, p=1e-6, pph4=0.85, fdr=1e-6)
    label = classify_target(af, [hf])
    assert label == "RED-1_AF_BENEFIT_HF_HARM"


def test_green2_shared_benefit():
    """Spec Test 5: AF benefit + HF + stroke benefit, all FDR<0.05 PPH4>0.8."""
    af = _eff("AF", theta=-0.20, p=1e-9, pph4=0.95, fdr=1e-9)
    hf = _eff("HF_overall", theta=-0.10, p=1e-4, pph4=0.90, fdr=1e-4)
    stroke = _eff("cardioembolic_stroke", theta=-0.12, p=1e-4, pph4=0.85, fdr=1e-4)
    label = classify_target(af, [hf, stroke])
    assert label == "GREEN-2_SHARED_BENEFIT"


def test_red2_af_benefit_stroke_harm():
    af = _eff("AF", theta=-0.10, p=1e-7, pph4=0.85, fdr=1e-7)
    stroke = _eff("ischemic_stroke", theta=0.20, p=1e-5, pph4=0.85, fdr=1e-5)
    label = classify_target(af, [stroke])
    assert label == "RED-2_AF_BENEFIT_STROKE_HARM"


def test_red3_multi_harm():
    af = _eff("AF", theta=-0.10, p=1e-7, pph4=0.85, fdr=1e-7)
    hf = _eff("HF_overall", theta=0.10, p=1e-4, pph4=0.85, fdr=1e-4)
    stroke = _eff("ischemic_stroke", theta=0.15, p=1e-4, pph4=0.85, fdr=1e-4)
    label = classify_target(af, [hf, stroke])
    assert label == "RED-3_AF_BENEFIT_MULTI_HARM"


def test_green1_neutral_safety():
    af = _eff("AF", theta=-0.10, p=1e-7, pph4=0.85, fdr=1e-7)
    hf = _eff("HF_overall", theta=0.01, p=0.50, pph4=0.10, fdr=0.50)
    stroke = _eff("ischemic_stroke", theta=-0.01, p=0.50, pph4=0.10, fdr=0.50)
    label = classify_target(af, [hf, stroke])
    assert label == "GREEN-1_AF_BENEFIT_NEUTRAL_SAFETY"


def test_amber1_unresolved_safety():
    af = _eff("AF", theta=-0.10, p=1e-7, pph4=0.85, fdr=1e-7)
    # Significant on p but coloc weak -> unresolved
    hf = _eff("HF_overall", theta=0.10, p=1e-3, pph4=0.10, fdr=1e-3)
    label = classify_target(af, [hf])
    assert label == "AMBER-1_AF_BENEFIT_UNRESOLVED"


def test_purple_hfstroke_benefit_af_harm():
    af = _eff("AF", theta=0.20, p=1e-8, pph4=0.85, fdr=1e-8)
    hf = _eff("HF_overall", theta=-0.10, p=1e-4, pph4=0.90, fdr=1e-4)
    label = classify_target(af, [hf])
    assert label == "PURPLE_HF_STROKE_BENEFIT_AF_HARM"


def test_blue_non_af_shared_mediator():
    af = _eff("AF", theta=-0.05, p=0.50, pph4=0.10, fdr=0.50)  # AF weak
    hf = _eff("HF_overall", theta=-0.10, p=1e-4, pph4=0.90, fdr=1e-4)
    label = classify_target(af, [hf])
    assert label == "BLUE_NON_AF_SHARED_MEDIATOR"


def test_unclassified_when_nothing_significant():
    af = _eff("AF", theta=-0.05, p=0.50, pph4=0.10, fdr=0.50)
    hf = _eff("HF_overall", theta=0.01, p=0.50, pph4=0.10, fdr=0.50)
    label = classify_target(af, [hf])
    assert label == "UNCLASSIFIED"
