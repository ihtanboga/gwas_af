"""Protocol decisions (D001 / D005 / D008 REVISED + D024) — regression tests."""

from afshf.intervention_direction import TargetOutcomeEffect
from afshf.scoring import classify_target, exploratory_concerns


def _eff(outcome, theta, p, pph4, fdr=None):
    fdr = fdr if fdr is not None else p
    return TargetOutcomeEffect(
        target_id="X", outcome=outcome,
        beta=theta, se=0.05, p=p, fdr=fdr, pph4=pph4,
        action_direction=1, theta=theta, or_intervention=1.0,
        benefit_flag=False, harm_flag=False,
    )


# ---------------------------------------------------------------------------
# D001 REVISED — HFpEF cannot trigger primary RED/GREEN classes
# ---------------------------------------------------------------------------


def test_hfpef_alone_does_not_trigger_red1():
    af = _eff("AF", theta=-0.20, p=1e-9, pph4=0.95, fdr=1e-9)
    hfpef_harm = _eff("HF_ni_HFpEF", theta=0.15, p=1e-6, pph4=0.85, fdr=1e-6)
    label = classify_target(af, [hfpef_harm])
    assert label != "RED-1_AF_BENEFIT_HF_HARM"
    # Without harm in primary HF set and no shared benefit, falls through to
    # GREEN-1 (HFpEF treated as exploratory only).
    assert label == "GREEN-1_AF_BENEFIT_NEUTRAL_SAFETY"


def test_hfref_still_triggers_red1():
    af = _eff("AF", theta=-0.20, p=1e-9, pph4=0.95, fdr=1e-9)
    hfref_harm = _eff("HF_ni_HFrEF", theta=0.15, p=1e-6, pph4=0.85, fdr=1e-6)
    label = classify_target(af, [hfref_harm])
    assert label == "RED-1_AF_BENEFIT_HF_HARM"


def test_hf_overall_still_triggers_red1():
    af = _eff("AF", theta=-0.20, p=1e-9, pph4=0.95, fdr=1e-9)
    hf_harm = _eff("HF_overall", theta=0.10, p=1e-4, pph4=0.85, fdr=1e-4)
    label = classify_target(af, [hf_harm])
    assert label == "RED-1_AF_BENEFIT_HF_HARM"


def test_hfpef_signal_is_reported_as_exploratory_concern():
    af = _eff("AF", theta=-0.20, p=1e-9, pph4=0.95, fdr=1e-9)
    hfpef_harm = _eff("HF_ni_HFpEF", theta=0.15, p=1e-3, pph4=0.65, fdr=1e-3)
    flags = exploratory_concerns(af, [hfpef_harm])
    assert "hfpef_signal" in flags
    assert "harm" in flags["hfpef_signal"]


def test_hfpef_benefit_also_appears_as_exploratory():
    af = _eff("AF", theta=-0.20, p=1e-9, pph4=0.95, fdr=1e-9)
    hfpef_benefit = _eff("HF_ni_HFpEF", theta=-0.15, p=1e-3, pph4=0.85, fdr=1e-3)
    flags = exploratory_concerns(af, [hfpef_benefit])
    assert "benefit" in flags["hfpef_signal"]


# ---------------------------------------------------------------------------
# D005 REVISED — weak AF tier cannot reach primary GREEN/RED
# ---------------------------------------------------------------------------


def test_weak_af_tier_does_not_trigger_green2():
    """AF p<0.05 but PPH4<0.50 -> weak tier; should not promote to GREEN-2."""
    af = _eff("AF", theta=-0.10, p=1e-3, pph4=0.20, fdr=1e-3)  # weak
    hf_benefit = _eff("HF_overall", theta=-0.10, p=1e-4, pph4=0.85, fdr=1e-4)
    label = classify_target(af, [hf_benefit])
    # AF tier is weak -> not eligible for GREEN-2; HF benefit alone routes BLUE.
    assert label == "BLUE_NON_AF_SHARED_MEDIATOR"


def test_weak_af_tier_emits_exploratory_flag():
    af = _eff("AF", theta=-0.10, p=1e-3, pph4=0.20, fdr=1e-3)
    flags = exploratory_concerns(af, [])
    assert flags.get("weak_af_direction") == ["only"]


def test_moderate_af_tier_reaches_green2():
    # p<0.05, PPH4=0.55 -> moderate; HF benefit at PPH4>=0.80 -> GREEN-2
    af = _eff("AF", theta=-0.10, p=1e-3, pph4=0.55, fdr=0.10)
    hf_benefit = _eff("HF_overall", theta=-0.10, p=1e-4, pph4=0.85, fdr=1e-4)
    label = classify_target(af, [hf_benefit])
    assert label == "GREEN-2_SHARED_BENEFIT"


# ---------------------------------------------------------------------------
# D008 REVISED — target_class is severity-objective
# ---------------------------------------------------------------------------


def test_classify_target_does_not_use_severity_weights():
    """classify_target() takes no severity_weights argument; class invariant
    under any severity scheme. Regression: ensure signature is severity-free."""
    import inspect
    sig = inspect.signature(classify_target)
    assert "severity_weights" not in sig.parameters


# ---------------------------------------------------------------------------
# D024 — disease-prior tag derivation from positive controls
# ---------------------------------------------------------------------------


def test_disease_prior_tags_from_positive_controls():
    from pathlib import Path

    from afshf.annotate_druggability import (
        DARK_BUT_DRUGGABLE,
        DruggabilityRecord,
        KNOWN_AF_TARGET,
        KNOWN_CVD_TARGET,
        KNOWN_HF_TARGET,
        KNOWN_STROKE_TARGET,
        annotate_disease_priors,
    )
    from afshf.positive_controls import load_positive_controls

    pc = load_positive_controls(
        Path(__file__).resolve().parents[1] / "config" / "positive_controls.tsv"
    )

    # PITX2 is registered for AF + cardioembolic_stroke -> two trait tags + CVD.
    pitx2 = annotate_disease_priors(
        DruggabilityRecord(target_id="PITX2", gene_symbol="PITX2"),
        positive_controls=pc,
    )
    assert KNOWN_AF_TARGET in pitx2.disease_prior_tags
    assert KNOWN_STROKE_TARGET in pitx2.disease_prior_tags
    assert KNOWN_CVD_TARGET in pitx2.disease_prior_tags

    # BAG3 is registered for HF only.
    bag3 = annotate_disease_priors(
        DruggabilityRecord(target_id="BAG3", gene_symbol="BAG3"),
        positive_controls=pc,
    )
    assert KNOWN_HF_TARGET in bag3.disease_prior_tags
    assert KNOWN_AF_TARGET not in bag3.disease_prior_tags
    assert KNOWN_CVD_TARGET in bag3.disease_prior_tags

    # Random non-registered target marked as Tdark -> dark_but_druggable.
    novel = annotate_disease_priors(
        DruggabilityRecord(target_id="UNKN1", gene_symbol="UNKN1",
                           pharos_tdl="Tdark"),
        positive_controls=pc,
    )
    assert DARK_BUT_DRUGGABLE in novel.disease_prior_tags


def test_known_drug_target_flag_from_clinical_phase():
    from afshf.annotate_druggability import (
        DruggabilityRecord,
        annotate_disease_priors,
    )

    record = annotate_disease_priors(
        DruggabilityRecord(target_id="PCSK9", gene_symbol="PCSK9",
                           max_clinical_phase=4, known_drugs_count=3),
    )
    assert record.known_drug_target is True

    record2 = annotate_disease_priors(
        DruggabilityRecord(target_id="UNKN2", gene_symbol="UNKN2",
                           max_clinical_phase=0, known_drugs_count=0),
    )
    assert record2.known_drug_target is False
