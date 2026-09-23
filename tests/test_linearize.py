"""Phase 2 checks: equilibrium validity, branch-aware Jacobian structure,
ZOH discretization sanity. The branch-derivative formulas themselves are
verified against sympy in scripts/_verify_branch_derivatives.py (a
one-off script, not part of this test suite, since it doesn't import the
package -- it's an independent symbolic cross-check of the hand algebra
that these formulas were then transcribed from).
"""

import numpy as np
import pytest

from muscle_activation_control import joint, linearize, muscle


THETA_STAR = np.deg2rad(60.0)


def test_equilibrium_is_isometric_torque_balance():
    a_star = linearize.find_equilibrium(THETA_STAR)
    assert 0.0 < a_star < 1.0
    tau_m = joint.isometric_torque(THETA_STAR, a_star)
    tau_g = joint.gravity_torque(THETA_STAR)
    assert tau_m == pytest.approx(tau_g, rel=1e-8)


def test_equilibrium_activation_derivative_is_zero_at_u_equals_a():
    a_star = linearize.find_equilibrium(THETA_STAR)
    assert muscle.activation_derivative(a_star, a_star) == pytest.approx(0.0, abs=1e-12)


def test_a_matrix_kinematic_row_exact():
    a_star = linearize.find_equilibrium(THETA_STAR)
    a_cont, _ = linearize.linearize(THETA_STAR, a_star, "activating", "shortening")
    assert a_cont[0, 0] == 0.0
    assert a_cont[0, 1] == 1.0
    assert a_cont[0, 2] == 0.0


def test_b_matrix_only_activation_row_nonzero():
    a_star = linearize.find_equilibrium(THETA_STAR)
    _, b_cont = linearize.linearize(THETA_STAR, a_star, "activating", "shortening")
    assert b_cont[0, 0] == 0.0
    assert b_cont[1, 0] == 0.0
    assert b_cont[2, 0] != 0.0


def test_activation_branch_asymmetry_matches_time_constant_ratio():
    # docs/theory/phase2_linearization.md section 3: the two activation
    # branches should differ by roughly tau_deact/tau_act scaled by a
    # (0.5+1.5a*)^2 factor -- not equal, and not wildly inconsistent with
    # that closed form.
    a_star = linearize.find_equilibrium(THETA_STAR)
    a_cont_act, b_cont_act = linearize.linearize(THETA_STAR, a_star, "activating", "shortening")
    a_cont_deact, b_cont_deact = linearize.linearize(
        THETA_STAR, a_star, "deactivating", "shortening"
    )
    eig_activating = a_cont_act[2, 2]
    eig_deactivating = a_cont_deact[2, 2]
    assert eig_activating != pytest.approx(eig_deactivating)
    # eig_activating = -1/(tau_act*base), eig_deactivating = -base/tau_deact
    # ratio = eig_activating/eig_deactivating = tau_deact / (tau_act * base**2)
    base = 0.5 + 1.5 * a_star
    expected_ratio = muscle.TAU_DEACT / (muscle.TAU_ACT * base**2)
    actual_ratio = eig_activating / eig_deactivating
    assert actual_ratio == pytest.approx(expected_ratio, rel=1e-8)


def test_velocity_branch_katz_asymmetry_is_exactly_2x():
    a_star = linearize.find_equilibrium(THETA_STAR)
    a_cont_short, _ = linearize.linearize(THETA_STAR, a_star, "activating", "shortening")
    a_cont_long, _ = linearize.linearize(THETA_STAR, a_star, "activating", "lengthening")
    d_short = a_cont_short[1, 1]
    d_long = a_cont_long[1, 1]
    assert d_short != 0.0
    assert d_long / d_short == pytest.approx(2.0, rel=1e-6)


def test_zoh_discretization_approaches_identity_as_ts_shrinks():
    a_star = linearize.find_equilibrium(THETA_STAR)
    a_cont, b_cont = linearize.linearize(THETA_STAR, a_star, "activating", "shortening")
    a_d, _ = linearize.discretize_zoh(a_cont, b_cont, 1e-8)
    assert np.allclose(a_d, np.eye(3), atol=1e-4)


def test_conditioning_sweep_runs_all_branch_combinations():
    a_star, results = linearize.conditioning_sweep(THETA_STAR, [0.001, 0.01])
    assert len(results) == 4 * 2  # 4 branch combos x 2 sample times
    for r in results:
        assert np.all(np.isfinite(r["a_d"]))
