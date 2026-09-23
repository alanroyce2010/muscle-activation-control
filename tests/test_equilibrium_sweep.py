"""Phase 6a checks: the closed-form activation-branch ratio matches the
linearization code, the constant-torque held-load equilibrium is a true
equilibrium of the loaded model, and the Phase 4 result reproduces from
the generalized helpers."""

import numpy as np
import pytest

from muscle_activation_control import joint, linearize, muscle, prediction


def _rho_a(a_star):
    return muscle.TAU_DEACT / (muscle.TAU_ACT * (0.5 + 1.5 * a_star) ** 2)


@pytest.mark.parametrize("a_star", [0.0, 0.133, 0.5, 0.884, 1.0])
def test_activation_ratio_closed_form_matches_linearize(a_star):
    theta = np.deg2rad(65.0)
    a_act, _ = linearize.linearize(theta, a_star, "activating", "shortening")
    a_de, _ = linearize.linearize(theta, a_star, "deactivating", "shortening")
    assert a_act[2, 2] / a_de[2, 2] == pytest.approx(_rho_a(a_star), rel=1e-12)


def test_activation_ratio_unity_point():
    a_unity = (np.sqrt(muscle.TAU_DEACT / muscle.TAU_ACT) - 0.5) / 1.5
    assert _rho_a(a_unity) == pytest.approx(1.0, rel=1e-12)
    assert 0.88 < a_unity < 0.89


def test_velocity_ratio_is_exactly_two_everywhere():
    for theta_deg, a_star in [(20, 0.05), (65, 0.133), (110, 0.7), (65, 0.95)]:
        theta = np.deg2rad(theta_deg)
        a_len, _ = linearize.linearize(theta, a_star, "activating", "lengthening")
        a_sh, _ = linearize.linearize(theta, a_star, "activating", "shortening")
        assert a_len[1, 1] / a_sh[1, 1] == pytest.approx(2.0, rel=1e-12)


def test_held_load_equilibrium_is_stationary():
    theta, a_star = np.deg2rad(65.0), 0.6
    tau_ext = prediction.external_torque_for_equilibrium(theta, a_star)
    traj = prediction.true_trajectory(theta, a_star, 0.0, 20, 0.01, tau_ext)
    assert np.allclose(traj[-1], [theta, 0.0, a_star], atol=1e-9)


def test_reproduces_phase4_headline_numbers():
    # Phase 4 headline, delta_u=+0.01, 50 ms horizon, theta* = 60 deg with
    # the pixel-level PCHIP moment arm (2026-09-11; was 0.70 / 0.13 / 0.73 at
    # 65 deg with the by-eye table): wrong-activation 71%, wrong-velocity
    # 14%, both wrong 75% of the true swing.
    theta = np.deg2rad(60.0)
    a_star = linearize.find_equilibrium(theta)
    ef = prediction.error_fractions(theta, a_star, 0.01, 5, 0.01)
    assert ef["wrong_activation"] == pytest.approx(0.711, abs=0.01)
    assert ef["wrong_velocity"] == pytest.approx(0.138, abs=0.01)
    assert ef["both_wrong"] == pytest.approx(0.749, abs=0.01)
    assert ef["matched"] < 0.02
