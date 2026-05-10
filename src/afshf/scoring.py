"""Scoring and target classification (spec sections 5, 10, 13).

The pipeline produces, per target:

    AF_benefit_score      : evidence weight for predicted AF benefit
    shared_benefit_score  : evidence weight per non-AF outcome with predicted benefit
    harm_score            : evidence weight per outcome with predicted harm
    therapeutic_balance   : AF_benefit + sum(shared_benefit) - lambda * sum(harm)

`therapeutic_balance` is for *ranking only*. Class labels (GREEN/AMBER/RED/...)
are decided by the threshold-based ``classify_target`` function.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import isfinite, log10
from typing import Iterable, Optional

from .intervention_direction import (
    PRIMARY_CLASS_ELIGIBLE_TIERS,
    TargetOutcomeEffect,
    af_signal_strength,
)


# D001 REVISED — HFpEF is exploratory. Primary class triggers must NOT depend
# on HFpEF alone. The full HF set still appears in plots/tables.
HF_OUTCOMES_PRIMARY = {"HF_overall", "HF_nonischemic", "HF_ni_HFrEF"}
HF_OUTCOMES_EXPLORATORY = {"HF_ni_HFpEF"}
HF_OUTCOMES = HF_OUTCOMES_PRIMARY | HF_OUTCOMES_EXPLORATORY

STROKE_OUTCOMES = {
    "ischemic_stroke",
    "cardioembolic_stroke",
    "large_artery_stroke",
    "small_vessel_stroke",
}


# ---------------------------------------------------------------------------
# Component scores
# ---------------------------------------------------------------------------


def _bounded_logp(p: float | None, cap: float) -> float:
    if p is None or not isfinite(p) or p <= 0.0:
        return 0.0
    return min(cap, -log10(p))


def coloc_weight(pph4: Optional[float],
                 strong: float = 1.0,
                 strong_threshold: float = 0.80,
                 moderate: float = 0.6,
                 moderate_threshold: float = 0.50,
                 weak: float = 0.2) -> float:
    if pph4 is None or not isfinite(pph4):
        return weak
    if pph4 >= strong_threshold:
        return strong
    if pph4 >= moderate_threshold:
        return moderate
    return weak


def replication_weight(status: str,
                       replicated: float = 1.0,
                       partial: float = 0.6,
                       none: float = 0.3) -> float:
    return {"replicated": replicated, "partial": partial, "none": none}.get(status, none)


def af_benefit_score(effect_af: TargetOutcomeEffect,
                     replication_status: str = "none",
                     max_logp: float = 10.0) -> float:
    """Score the AF-benefit signal. Zero when AF effect is not protective."""
    if not (isfinite(effect_af.theta) and effect_af.theta < 0):
        return 0.0
    return (
        _bounded_logp(effect_af.p, max_logp)
        * coloc_weight(effect_af.pph4)
        * replication_weight(replication_status)
    )


def harm_score(effect: TargetOutcomeEffect,
               severity_weight: float,
               replication_status: str = "none",
               max_logp: float = 10.0) -> float:
    if not (isfinite(effect.theta) and effect.theta > 0):
        return 0.0
    return (
        _bounded_logp(effect.p, max_logp)
        * coloc_weight(effect.pph4)
        * replication_weight(replication_status)
        * severity_weight
    )


def shared_benefit_score(effect: TargetOutcomeEffect,
                         replication_status: str = "none",
                         max_logp: float = 10.0) -> float:
    if not (isfinite(effect.theta) and effect.theta < 0):
        return 0.0
    return (
        _bounded_logp(effect.p, max_logp)
        * coloc_weight(effect.pph4)
        * replication_weight(replication_status)
    )


# ---------------------------------------------------------------------------
# Target-level aggregate
# ---------------------------------------------------------------------------


@dataclass
class TargetScore:
    target_id: str
    af_benefit: float = 0.0
    shared_benefit_total: float = 0.0
    harm_total: float = 0.0
    therapeutic_balance: float = 0.0
    components: dict[str, float] = field(default_factory=dict)


def therapeutic_balance(af_benefit: float,
                        shared_benefit_total: float,
                        harm_total: float,
                        lambda_harm: float = 1.5) -> float:
    return af_benefit + shared_benefit_total - lambda_harm * harm_total


def aggregate_target_score(target_id: str,
                           af_effect: TargetOutcomeEffect,
                           outcome_effects: Iterable[TargetOutcomeEffect],
                           severity_weights: dict[str, float],
                           replication_status: dict[str, str] | None = None,
                           lambda_harm: float = 1.5,
                           max_logp: float = 10.0) -> TargetScore:
    replication_status = replication_status or {}
    af_score = af_benefit_score(
        af_effect, replication_status.get("AF", "none"), max_logp=max_logp
    )

    shared_total = 0.0
    harm_total = 0.0
    components: dict[str, float] = {"AF_benefit": af_score}

    for effect in outcome_effects:
        if effect.outcome == "AF":
            continue
        rep = replication_status.get(effect.outcome, "none")
        sev = severity_weights.get(effect.outcome, 1.0)
        sb = shared_benefit_score(effect, rep, max_logp=max_logp)
        hs = harm_score(effect, sev, rep, max_logp=max_logp)
        shared_total += sb
        harm_total += hs
        components[f"{effect.outcome}_benefit"] = sb
        components[f"{effect.outcome}_harm"] = hs

    balance = therapeutic_balance(af_score, shared_total, harm_total, lambda_harm)
    return TargetScore(
        target_id=target_id,
        af_benefit=af_score,
        shared_benefit_total=shared_total,
        harm_total=harm_total,
        therapeutic_balance=balance,
        components=components,
    )


# ---------------------------------------------------------------------------
# Classification (spec section 5)
# ---------------------------------------------------------------------------


def _significant_at(effect: TargetOutcomeEffect,
                    fdr_threshold: float,
                    pph4_threshold: float) -> bool:
    """D026 — prefer ``neg_log10p`` for the significance gate when available.

    A pval=0 underflow is not unusual at canonical loci (e.g. PITX2 in AF).
    Using -log10(p) keeps these signals visible even when raw p has rounded
    to zero in storage.
    """
    if not isfinite(effect.theta):
        return False

    nlp = effect.neg_log10p
    if nlp is not None and isfinite(nlp):
        # neg_log10p significance: nlp > -log10(threshold).
        if nlp <= -log10_safe(fdr_threshold):
            return False
    else:
        fdr = effect.fdr if effect.fdr is not None else effect.p
        if fdr is None or not isfinite(fdr) or fdr >= fdr_threshold:
            return False

    if effect.pph4 is None or not isfinite(effect.pph4):
        return False
    if effect.pph4 < pph4_threshold:
        return False
    return True


def log10_safe(x: float) -> float:
    """log10 that returns -inf for non-positive values without raising."""
    from math import log10
    if x is None or not isfinite(x) or x <= 0.0:
        return float("-inf")
    return log10(x)


def classify_target(af_effect: TargetOutcomeEffect,
                    outcome_effects: Iterable[TargetOutcomeEffect],
                    fdr_threshold: float = 0.05,
                    pph4_threshold: float = 0.50,
                    pph4_strong: float = 0.80) -> str:
    """Return the spec class label for a target.

    D005 REVISED — primary GREEN/RED labels require AF tier in
    PRIMARY_CLASS_ELIGIBLE_TIERS (`strong` or `moderate`).
    D001 REVISED — HFpEF is exploratory; it does NOT trigger RED/GREEN.
    D008 REVISED — severity weights MUST NOT influence target_class.

    Order of precedence:
        PURPLE  > BLUE  > UNCLASSIFIED   (when AF tier is none/weak)
        RED-3 > RED-1 / RED-2 > GREEN-2 > GREEN-1 > AMBER-1
                                          (when AF tier is strong/moderate)
    """

    outcome_effects = list(outcome_effects)
    af_tier = af_signal_strength(
        p_target_to_AF=af_effect.p,
        fdr_target_to_AF=af_effect.fdr,
        pph4_AF=af_effect.pph4,
        p_threshold=0.05,
        fdr_threshold=fdr_threshold,
        pph4_strong=pph4_strong,
        pph4_moderate=pph4_threshold,
    )
    af_protective = isfinite(af_effect.theta) and af_effect.theta < 0
    af_eligible_protective = af_protective and af_tier in PRIMARY_CLASS_ELIGIBLE_TIERS
    af_eligible_harm = (
        isfinite(af_effect.theta)
        and af_effect.theta > 0
        and af_tier in PRIMARY_CLASS_ELIGIBLE_TIERS
    )

    hf_harm = [
        e for e in outcome_effects
        if e.outcome in HF_OUTCOMES_PRIMARY      # D001 — HFpEF excluded
        and isfinite(e.theta) and e.theta > 0
        and _significant_at(e, fdr_threshold, pph4_threshold)
    ]
    stroke_harm = [
        e for e in outcome_effects
        if e.outcome in STROKE_OUTCOMES
        and isfinite(e.theta) and e.theta > 0
        and _significant_at(e, fdr_threshold, pph4_threshold)
    ]
    hf_benefit = [
        e for e in outcome_effects
        if e.outcome in HF_OUTCOMES_PRIMARY      # D001 — HFpEF excluded
        and isfinite(e.theta) and e.theta < 0
        and _significant_at(e, fdr_threshold, pph4_strong)
    ]
    stroke_benefit = [
        e for e in outcome_effects
        if e.outcome in STROKE_OUTCOMES
        and isfinite(e.theta) and e.theta < 0
        and _significant_at(e, fdr_threshold, pph4_strong)
    ]

    # PURPLE: AF harm but HF/stroke benefit
    if af_eligible_harm and (hf_benefit or stroke_benefit):
        return "PURPLE_HF_STROKE_BENEFIT_AF_HARM"

    # BLUE: AF tier not eligible but shared HF/stroke benefit elsewhere
    if not af_eligible_protective and (hf_benefit or stroke_benefit):
        return "BLUE_NON_AF_SHARED_MEDIATOR"

    if not af_eligible_protective:
        return "UNCLASSIFIED"

    # Below this point the target has an eligible AF-protective signal.
    if hf_harm and stroke_harm:
        return "RED-3_AF_BENEFIT_MULTI_HARM"
    if hf_harm:
        return "RED-1_AF_BENEFIT_HF_HARM"
    if stroke_harm:
        return "RED-2_AF_BENEFIT_STROKE_HARM"
    if hf_benefit or stroke_benefit:
        return "GREEN-2_SHARED_BENEFIT"

    # No harm and no shared-benefit -> evaluate neutrality vs unresolved.
    # AMBER-1 inspects PRIMARY HF + stroke outcomes only (HFpEF goes through
    # the exploratory_concerns channel instead).
    primary_safety_set = HF_OUTCOMES_PRIMARY | STROKE_OUTCOMES
    weakly_significant_unresolved = any(
        e.outcome in primary_safety_set
        and isfinite(e.p) and e.p < 0.05
        and (e.pph4 is None or e.pph4 < pph4_threshold)
        for e in outcome_effects
    )
    if weakly_significant_unresolved:
        return "AMBER-1_AF_BENEFIT_UNRESOLVED"
    return "GREEN-1_AF_BENEFIT_NEUTRAL_SAFETY"


# ---------------------------------------------------------------------------
# Exploratory concerns — never affect class, always reported in master table
# ---------------------------------------------------------------------------


def exploratory_concerns(af_effect: TargetOutcomeEffect,
                         outcome_effects: Iterable[TargetOutcomeEffect],
                         *,
                         p_threshold: float = 0.05,
                         pph4_moderate: float = 0.50,
                         pph4_strong: float = 0.80) -> dict[str, list[str]]:
    """Return per-target exploratory flags surfaced alongside target_class.

    Shape:
        {
          "hfpef_signal":           ["harm" | "benefit"]   if HFpEF significant,
          "weak_af_direction":      [True]                 if AF tier == weak,
          "ld_confounding_warning": ["AF" | outcome name]  if PPH3 high vs PPH4 low,
        }
    """
    outcome_effects = list(outcome_effects)
    out: dict[str, list[str]] = {}

    for e in outcome_effects:
        if e.outcome in HF_OUTCOMES_EXPLORATORY:
            if (isfinite(e.p) and e.p < p_threshold and
                e.pph4 is not None and e.pph4 >= pph4_moderate):
                direction = "harm" if e.theta > 0 else "benefit"
                out.setdefault("hfpef_signal", []).append(direction)

    af_tier = af_signal_strength(
        p_target_to_AF=af_effect.p,
        fdr_target_to_AF=af_effect.fdr,
        pph4_AF=af_effect.pph4,
        p_threshold=p_threshold,
        fdr_threshold=0.05,
        pph4_strong=pph4_strong,
        pph4_moderate=pph4_moderate,
    )
    if af_tier == "weak":
        out["weak_af_direction"] = ["only"]

    return out
