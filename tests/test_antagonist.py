"""Phase 6c checks: the pair reduces to the single-muscle model when the
triceps is off; co-contraction equilibria are stationary; all three
switches satisfy the Camlibel et al. continuity condition with disjoint
support; the closed-form force-velocity ratio matches finite differences
and equals 2 when the triceps is off."""

import numpy as np
import pytest

from muscle_activation_control import antagonist, joint

THETA = np.deg2rad(60.0)


def test_pair_reduces_to_single_muscle_when_triceps_off():
    # at THETA_REF the triceps fiber is at l_0 (no passive force) and a_t=0
    x = np.array([THETA, 0.3, 0.2, 0.0])
    d_pair = antagonist.state_derivative(x, [0.25, 0.0])
    d_single = joint.state_derivative(x[:3], 0.25)
    assert np.allclose(d_pair[:3], d_single, rtol=1e-12, atol=1e-12)
    assert d_pair[3] == 0.0


def test_cocontraction_equilibrium_is_stationary():
    a_t = 0.2
    a_b = antagonist.find_equilibrium(THETA, a_t)
    assert a_b > antagonist.find_equilibrium(THETA, 0.0)
    f0 = antagonist.state_derivative([THETA, 0.0, a_b, a_t], [a_b, a_t])
    assert np.allclose(f0, 0.0, atol=1e-9)


@pytest.mark.parametrize("a_t", [0.0, 0.05, 0.3])
def test_three_switches_satisfy_continuity_with_disjoint_support(a_t):
    a_b = antagonist.find_equilibrium(THETA, a_t)
    v = antagonist.verify_switches(THETA, a_b, a_t)
    for k in ("bra_activation", "tri_activation", "velocity"):
        assert v[k]["residual"] < 1e-6
    assert v["disjoint"]
    assert abs(v["bra_activation"]["e"][2]) > 1.0
    assert abs(v["velocity"]["e"][1]) > 0.1
    if a_t == 0.0:
        # triceps at zero activation: its activation switch still exists
        # (tau_a depends on the branch even at a=0) but contributes no torque
        assert abs(v["tri_activation"]["e"][3]) > 1.0


def test_force_velocity_ratio_closed_form_matches_fd_and_reduces_to_two():
    for a_t in (0.0, 0.1, 0.3):
        a_b = antagonist.find_equilibrium(THETA, a_t)
        _, _, ratio = antagonist.damping_contributions(THETA, a_b, a_t)
        x = np.array([THETA, 0, a_b, a_t]); u = np.array([a_b, a_t])
        a_fl, _, _ = antagonist.linearize(x, u, "activating", "activating", "flexing")
        a_ex, _, _ = antagonist.linearize(x, u, "activating", "activating", "extending")
        assert a_ex[1, 1] / a_fl[1, 1] == pytest.approx(ratio, rel=1e-4)
        if a_t == 0.0:
            assert ratio == pytest.approx(2.0, rel=1e-6)
        else:
            assert 0.5 < ratio < 2.0


def test_balanced_cocontraction_cancels_the_kink():
    a_bal = antagonist.balanced_cocontraction(THETA)
    assert a_bal is not None and 0 < a_bal < 1
    a_b = antagonist.find_equilibrium(THETA, a_bal)
    d_b, d_t, ratio = antagonist.damping_contributions(THETA, a_b, a_bal)
    assert ratio == pytest.approx(1.0, abs=1e-6)
    v = antagonist.verify_switches(THETA, a_b, a_bal)
    assert abs(v["velocity"]["e"][1]) < 1e-4 * abs(d_b / joint.FOREARM_HAND_I_ELBOW)
