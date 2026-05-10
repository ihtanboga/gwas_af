"""Therapeutic intervention direction logic.

Implements the core spec rules (Section 7):

    beta_T_AF = MR beta for genetically increased target -> AF log(OR)
    d_T       = -sign(beta_T_AF)
        d_T = -1 -> therapeutic strategy is to inhibit / decrease the target
        d_T = +1 -> therapeutic strategy is to activate / increase the target
    theta_T_O = d_T * beta_T_O   for every outcome O
    OR_T_O    = exp(theta_T_O)

Interpretation:
    theta_T_AF       < 0 -> predicted AF benefit
    theta_T_HF       > 0 -> predicted HF harm
    theta_T_stroke   > 0 -> predicted stroke harm
    theta_T_O        < 0 -> predicted benefit for outcome O

Off-axis decomposition (Section 7 secondary analysis):

    expected_AF_mediated_effect_on_outcome = theta_T_AF * beta_AF_to_outcome
    excess_off_axis_effect = theta_T_outcome - expected_AF_mediated_effect_on_outcome
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, isfinite
from typing import Literal, Optional


ActionDirection = Literal[-1, 1, 0]
AfSignalStrength = Literal["strong", "moderate", "weak", "none"]


@dataclass(frozen=True)
class TargetOutcomeEffect:
    """Predicted effect of the chosen therapeutic intervention on one outcome."""

    target_id: str
    outcome: str
    beta: float                       # MR beta on outcome (raised target -> outcome)
    se: float
    p: float
    fdr: Optional[float]
    pph4: Optional[float]
    action_direction: ActionDirection  # d_T derived from AF
    theta: float                       # d_T * beta
    or_intervention: float             # exp(theta)
    benefit_flag: bool
    harm_flag: bool
    # D026 — underflow-safe significance signal. When provided, scoring routines
    # prefer ``neg_log10p`` over raw ``p`` so a pval=0 underflow does not
    # swithntly drop a credible signal.
    neg_log10p: Optional[float] = None


def derive_action_direction(beta_target_to_AF: float,
                            p_target_to_AF: float | None = None,
                            af_min_p_for_direction: float = 0.05) -> ActionDirection:
    """Return d_T for a target given its MR effect on AF.

    Rule:
        d_T = -sign(beta_T_AF)
        If beta is exactly zero or AF effect is not significant (when a p-value
        is provided), return 0 to indicate "no actionable direction yet".
    """

    if not isfinite(beta_target_to_AF):
        return 0
    if beta_target_to_AF == 0.0:
        return 0
    if p_target_to_AF is not None and isfinite(p_target_to_AF):
        if p_target_to_AF >= af_min_p_for_direction:
            return 0
    return -1 if beta_target_to_AF > 0 else 1


def af_signal_strength(p_target_to_AF: float | None,
                       fdr_target_to_AF: float | None,
                       pph4_AF: float | None,
                       *,
                       p_threshold: float = 0.05,
                       fdr_threshold: float = 0.05,
                       pph4_strong: float = 0.80,
                       pph4_moderate: float = 0.50) -> AfSignalStrength:
    """Tier the AF signal per spec D005 REVISED.

        strong   : FDR_AF < fdr_threshold  AND PPH4_AF >= pph4_strong
        moderate : p_AF   < p_threshold    AND PPH4_AF >= pph4_moderate
                   (FDR not required; e.g. exploratory FDR or genome-wide FDR fail)
        weak     : p_AF   < p_threshold    AND PPH4_AF <  pph4_moderate
        none     : p_AF   >= p_threshold   OR  beta missing
    """
    p_ok = p_target_to_AF is not None and isfinite(p_target_to_AF) and p_target_to_AF < p_threshold
    if not p_ok:
        return "none"

    fdr_ok = (
        fdr_target_to_AF is not None
        and isfinite(fdr_target_to_AF)
        and fdr_target_to_AF < fdr_threshold
    )
    pph4 = pph4_AF if (pph4_AF is not None and isfinite(pph4_AF)) else 0.0

    if fdr_ok and pph4 >= pph4_strong:
        return "strong"
    if pph4 >= pph4_moderate:
        return "moderate"
    return "weak"


PRIMARY_CLASS_ELIGIBLE_TIERS: frozenset[AfSignalStrength] = frozenset({"strong", "moderate"})


def project_effect(target_id: str,
                   outcome: str,
                   beta_outcome: float,
                   se_outcome: float,
                   p_outcome: float,
                   action_direction: ActionDirection,
                   fdr_outcome: Optional[float] = None,
                   pph4: Optional[float] = None,
                   benefit_p_threshold: float = 0.05,
                   benefit_pph4_threshold: float = 0.50) -> TargetOutcomeEffect:
    """Project the MR beta onto the chosen intervention direction.

    A "benefit" or "harm" flag fires only when both the statistical signal
    *and* a colocalization gate are satisfied. Callers may override the
    gates per-outcome.
    """

    theta = action_direction * beta_outcome if action_direction != 0 else float("nan")
    or_intervention = exp(theta) if isfinite(theta) else float("nan")

    significant = (p_outcome < benefit_p_threshold) if isfinite(p_outcome) else False
    coloc_ok = (pph4 is not None and pph4 >= benefit_pph4_threshold)
    has_signal = significant and coloc_ok and isfinite(theta)

    benefit_flag = bool(has_signal and theta < 0)
    harm_flag = bool(has_signal and theta > 0)

    return TargetOutcomeEffect(
        target_id=target_id,
        outcome=outcome,
        beta=beta_outcome,
        se=se_outcome,
        p=p_outcome,
        fdr=fdr_outcome,
        pph4=pph4,
        action_direction=action_direction,
        theta=theta,
        or_intervention=or_intervention,
        benefit_flag=benefit_flag,
        harm_flag=harm_flag,
    )


def expected_af_mediated_effect(theta_target_AF: float,
                                beta_AF_to_outcome: float) -> float:
    """Expected effect on outcome that flows purely through AF.

    Per Section 7:
        expected = theta_T_AF * beta_AF_to_outcome
    """
    return theta_target_AF * beta_AF_to_outcome


def excess_off_axis_effect(theta_target_outcome: float,
                           theta_target_AF: float,
                           beta_AF_to_outcome: float) -> float:
    """Residual outcome effect not explained by AF mediation."""
    expected = expected_af_mediated_effect(theta_target_AF, beta_AF_to_outcome)
    return theta_target_outcome - expected
