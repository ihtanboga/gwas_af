"""Spec section 18 — Tests 1 and 2 (therapeutic intervention direction)."""

from math import exp, isclose

from afshf.intervention_direction import (
    derive_action_direction,
    project_effect,
    expected_af_mediated_effect,
    excess_off_axis_effect,
)


# ---------------------------------------------------------------------------
# Test 1 — increase-of-target raises AF
# ---------------------------------------------------------------------------


def test_action_direction_inhibit_when_target_raises_AF():
    """Spec Test 1: beta_AF=+0.20, beta_HF=-0.10 -> d=-1, AF benefit, HF harm."""
    beta_AF = 0.20
    beta_HF = -0.10

    d = derive_action_direction(beta_AF, p_target_to_AF=1e-9)
    assert d == -1

    af_eff = project_effect(
        target_id="X", outcome="AF",
        beta_outcome=beta_AF, se_outcome=0.02, p_outcome=1e-9,
        action_direction=d, pph4=0.9,
    )
    assert isclose(af_eff.theta, -0.20)
    assert af_eff.benefit_flag is True
    assert af_eff.harm_flag is False

    hf_eff = project_effect(
        target_id="X", outcome="HF_overall",
        beta_outcome=beta_HF, se_outcome=0.03, p_outcome=1e-3,
        action_direction=d, pph4=0.9,
    )
    assert isclose(hf_eff.theta, 0.10)
    assert hf_eff.benefit_flag is False
    assert hf_eff.harm_flag is True
    assert isclose(hf_eff.or_intervention, exp(0.10))


# ---------------------------------------------------------------------------
# Test 2 — increase-of-target lowers AF
# ---------------------------------------------------------------------------


def test_action_direction_activate_when_target_lowers_AF():
    """Spec Test 2: beta_AF=-0.30, beta_stroke=-0.20 -> d=+1, AF + stroke benefit."""
    beta_AF = -0.30
    beta_stroke = -0.20

    d = derive_action_direction(beta_AF, p_target_to_AF=1e-12)
    assert d == 1

    af_eff = project_effect(
        target_id="Y", outcome="AF",
        beta_outcome=beta_AF, se_outcome=0.04, p_outcome=1e-12,
        action_direction=d, pph4=0.95,
    )
    assert isclose(af_eff.theta, -0.30)
    assert af_eff.benefit_flag is True

    stroke_eff = project_effect(
        target_id="Y", outcome="cardioembolic_stroke",
        beta_outcome=beta_stroke, se_outcome=0.05, p_outcome=1e-4,
        action_direction=d, pph4=0.85,
    )
    assert isclose(stroke_eff.theta, -0.20)
    assert stroke_eff.benefit_flag is True
    assert stroke_eff.harm_flag is False


def test_direction_zero_when_AF_signal_is_not_significant():
    d = derive_action_direction(beta_target_to_AF=0.10, p_target_to_AF=0.30)
    assert d == 0


def test_off_axis_decomposition():
    # Target reduces AF (theta_AF<0). AF in turn raises stroke risk (beta>0).
    # Expected AF-mediated stroke effect therefore is negative.
    theta_AF = -0.40
    beta_AF_to_stroke = 0.50
    expected = expected_af_mediated_effect(theta_AF, beta_AF_to_stroke)
    assert isclose(expected, -0.20)

    # Observed theta on stroke is even more protective than expected -> negative excess.
    excess = excess_off_axis_effect(
        theta_target_outcome=-0.30,
        theta_target_AF=theta_AF,
        beta_AF_to_outcome=beta_AF_to_stroke,
    )
    assert isclose(excess, -0.10)
