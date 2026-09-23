"""Phase 6b mechanics checks (not the scientific outcome, which is read
off scripts/phase6b_closed_loop.py): the branch-forced general
linearizer agrees with the analytic one at the equilibrium; the
'blended' central-difference linearization is the average of the two
one-sided slopes at each kink; the affine ZOH reduces to the plain one
at an equilibrium; and each MPC variant actually regulates."""

import numpy as np
import pytest

from muscle_activation_control import joint, linearize

THETA = np.deg2rad(60.0)


@pytest.mark.parametrize("act,vel", linearize.BRANCH_COMBINATIONS)
def test_general_linearizer_matches_analytic_at_equilibrium(act, vel):
    a_star = linearize.find_equilibrium(THETA)
    a_ref, b_ref = linearize.linearize(THETA, a_star, act, vel)
    a_fd, b_fd, f0 = linearize.linearize_general([THETA, 0.0, a_star], a_star, act, vel)
    assert np.allclose(a_fd, a_ref, rtol=1e-5, atol=1e-6)
    assert np.allclose(b_fd, b_ref, rtol=1e-5, atol=1e-6)
    assert np.allclose(f0, 0.0, atol=1e-8)


def test_blended_linearization_averages_the_kinks():
    a_star = linearize.find_equilibrium(THETA)
    x_eq = [THETA, 0.0, a_star]
    a_bl, b_bl, _ = linearize.linearize_general(x_eq, a_star, None, None)
    a_as, _ = linearize.linearize(THETA, a_star, "activating", "shortening")
    a_al, _ = linearize.linearize(THETA, a_star, "activating", "lengthening")
    a_ds, b_ds = linearize.linearize(THETA, a_star, "deactivating", "shortening")
    b_as = linearize.linearize(THETA, a_star, "activating", "shortening")[1]
    assert a_bl[1, 1] == pytest.approx(0.5 * (a_as[1, 1] + a_al[1, 1]), rel=1e-5)
    assert a_bl[2, 2] == pytest.approx(0.5 * (a_as[2, 2] + a_ds[2, 2]), rel=1e-5)
    assert b_bl[2, 0] == pytest.approx(0.5 * (b_as[2, 0] + b_ds[2, 0]), rel=1e-5)
    # kink-free entries unaffected
    assert a_bl[1, 0] == pytest.approx(a_as[1, 0], rel=1e-5)
    assert a_bl[1, 2] == pytest.approx(a_as[1, 2], rel=1e-5)


def test_forced_branches_do_not_change_default_model():
    x = np.array([np.deg2rad(70.0), 0.4, 0.3])
    d0 = joint.state_derivative(x, 0.5)
    assert np.allclose(d0, joint.state_derivative(x, 0.5, "activating", "shortening"))
    x = np.array([np.deg2rad(70.0), -0.4, 0.3])
    d0 = joint.state_derivative(x, 0.1)
    assert np.allclose(d0, joint.state_derivative(x, 0.1, "deactivating", "lengthening"))


def test_affine_zoh_reduces_to_plain_at_equilibrium():
    a_star = linearize.find_equilibrium(THETA)
    a_c, b_c = linearize.linearize(THETA, a_star, "activating", "shortening")
    a_d, b_d = linearize.discretize_zoh(a_c, b_c, 0.01)
    a_d2, b_d2, c_d = linearize.discretize_affine(a_c, b_c, np.zeros(3), 0.01)
    assert np.allclose(a_d, a_d2) and np.allclose(b_d, b_d2) and np.allclose(c_d, 0.0)


cvxpy = pytest.importorskip("cvxpy")
from muscle_activation_control import mpc  # noqa: E402


@pytest.mark.parametrize("mode", mpc.MODES)
def test_each_variant_regulates_a_small_step(mode):
    a0 = linearize.find_equilibrium(THETA)
    ref = lambda t: np.where(np.asarray(t) >= 0.1, np.deg2rad(62.0), THETA)
    ctrl = mpc.BranchAwareMPC(mode)
    log = mpc.simulate(ctrl, ref, 1.5, [THETA, 0.0, a0])
    final_err_deg = np.rad2deg(abs(log["x"][-1, 0] - np.deg2rad(62.0)))
    assert final_err_deg < 0.3
    assert np.all(np.isfinite(log["u"]))


def test_trajectory_mode_assigns_per_step_branches():
    a0 = linearize.find_equilibrium(THETA)
    ctrl = mpc.BranchAwareMPC("trajectory")
    preview = np.full(ctrl.N + 1, THETA + np.deg2rad(3.0))
    ctrl.compute_control([THETA, 0.0, a0], a0, preview)          # seeds the plan
    _, info = ctrl.compute_control([THETA, 0.0, a0], a0, preview)
    assert len(info["branches"]) == ctrl.N
    # a flexion step planned from rest: the first steps must be activating/shortening
    assert info["branches"][0] == ("activating", "shortening")


def test_tracking_rule_picks_direction_at_rest():
    ctrl = mpc.BranchAwareMPC("tracking")
    a0 = linearize.find_equilibrium(THETA)
    assert ctrl.select_branch([THETA, 0.0, a0], a0, THETA + 0.01) == ("activating", "shortening")
    assert ctrl.select_branch([THETA, 0.0, a0], a0, THETA - 0.01) == ("deactivating", "lengthening")
    assert ctrl.select_branch([THETA, -0.5, a0], a0 + 0.1, THETA + 0.01) == ("activating", "lengthening")


def test_hybrid_mode_regulates_a_small_step():
    # coarse switch grid (never / immediately) keeps this fast; the full
    # grid is exercised by scripts/phase7_realism.py
    a0 = linearize.find_equilibrium(THETA)
    ref = lambda t: np.where(np.asarray(t) >= 0.05, np.deg2rad(62.0), THETA)
    ctrl = mpc.BranchAwareMPC("hybrid", switch_grid=[0, 10, 20])
    log = mpc.simulate(ctrl, ref, 0.8, [THETA, 0.0, a0])
    assert np.rad2deg(abs(log["x"][-1, 0] - np.deg2rad(62.0))) < 0.5
    assert len(ctrl.hyb_violations) == len(log["u"])


def test_realistic_simulator_reduces_to_ideal_when_all_factors_off():
    a0 = linearize.find_equilibrium(THETA)
    ref = lambda t: np.where(np.asarray(t) >= 0.1, np.deg2rad(62.0), THETA)
    ideal = mpc.simulate(mpc.BranchAwareMPC("fixed_naive"), ref, 0.5, [THETA, 0.0, a0])
    real = mpc.simulate_realistic(mpc.BranchAwareMPC("fixed_naive"), ref, 0.5, [THETA, 0.0, a0])
    assert np.allclose(real["x"][:, 0], ideal["x"][:, 0], atol=1e-6)


def test_delay_predictor_recovers_undelayed_response_on_nominal_plant():
    a0 = linearize.find_equilibrium(THETA)
    ref = lambda t: np.where(np.asarray(t) >= 0.3, np.deg2rad(62.0), THETA)
    x0 = [THETA, 0.0, a0]
    base = mpc.simulate_realistic(mpc.BranchAwareMPC("fixed_naive"), ref, 0.8, x0)
    comp = mpc.simulate_realistic(mpc.BranchAwareMPC("fixed_naive"), ref, 0.8, x0,
                                  delay_steps=3, compensate_delay=True)
    raw = mpc.simulate_realistic(mpc.BranchAwareMPC("fixed_naive"), ref, 0.8, x0, delay_steps=3)
    err_comp = np.max(np.abs(comp["x"][:, 0] - base["x"][:, 0]))
    err_raw = np.max(np.abs(raw["x"][:, 0] - base["x"][:, 0]))
    assert err_comp < np.deg2rad(0.01)
    assert err_raw > 10 * err_comp


def test_blended_activation_mode_averages_only_the_activation_row():
    a0 = linearize.find_equilibrium(THETA)
    A_m, B_m, _ = linearize.linearize_general([THETA, 0.0, a0], a0, None, "shortening")
    A_act, B_act = linearize.linearize(THETA, a0, "activating", "shortening")
    A_dea, B_dea = linearize.linearize(THETA, a0, "deactivating", "shortening")
    assert np.allclose(A_m[2], 0.5 * (A_act[2] + A_dea[2]), rtol=1e-4, atol=1e-6)
    assert np.allclose(B_m[2], 0.5 * (B_act[2] + B_dea[2]), rtol=1e-4, atol=1e-6)
    assert np.allclose(A_m[:2], A_act[:2], rtol=1e-4, atol=1e-6)
