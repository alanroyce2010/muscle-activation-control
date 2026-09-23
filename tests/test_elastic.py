"""Phase 7 section 4 checks: the elastic-tendon brachialis keeps the
two-switch structure (the velocity switch moves into the fibre row),
satisfies the continuity condition with disjoint support, and approaches
the rigid-tendon model as the tendon stiffens."""

import numpy as np
import pytest

from muscle_activation_control import elastic, linearize

TH = np.deg2rad(60.0)


def test_tendon_curve_is_c1_at_toe_transition():
    et = 0.609 * elastic.EPS0_T
    lo, hi = elastic.tendon_force_norm(et - 1e-7), elastic.tendon_force_norm(et + 1e-7)
    assert lo == pytest.approx(hi, abs=1e-5)
    s_lo = (elastic.tendon_force_norm(et) - elastic.tendon_force_norm(et - 1e-6)) / 1e-6
    s_hi = (elastic.tendon_force_norm(et + 1e-6) - elastic.tendon_force_norm(et)) / 1e-6
    assert s_lo == pytest.approx(s_hi, rel=1e-3)
    assert elastic.tendon_strain(float(elastic.tendon_force_norm(0.01))) == pytest.approx(0.01, rel=1e-9)


def test_equilibrium_is_stationary_and_on_the_fibre_switch():
    l_bar, a = elastic.find_equilibrium(TH)
    x = np.array([TH, 0.0, l_bar, a])
    assert np.allclose(elastic.state_derivative(x, a), 0.0, atol=1e-9)
    assert elastic.switching_function(x) == pytest.approx(0.0, abs=1e-12)


def test_switches_satisfy_continuity_with_disjoint_support():
    l_bar, a = elastic.find_equilibrium(TH)
    x = np.array([TH, 0.0, l_bar, a])
    A_as, b_as = elastic.linearize(x, a, "activating", "shortening")
    A_ds, b_ds = elastic.linearize(x, a, "deactivating", "shortening")
    A_al, b_al = elastic.linearize(x, a, "activating", "lengthening")
    e_a, r_a = elastic.continuity(A_as - A_ds, b_as - b_ds, np.array([0, 0, 0, -1.0]), 1.0)
    c_v = elastic.switching_gradient(x)
    e_v, r_v = elastic.continuity(A_al - A_as, b_al - b_as, c_v, 0.0)
    assert r_a < 1e-5 and r_v < 1e-5
    assert np.flatnonzero(np.abs(e_a) > 1e-6).tolist() == [3]     # activation row only
    assert np.flatnonzero(np.abs(e_v) > 1e-6).tolist() == [2]     # fibre-length row only
    # the thetaddot row carries no branch dependence at all in the elastic model
    assert np.allclose(A_al[1], A_as[1], atol=1e-6) and np.allclose(A_ds[1], A_as[1], atol=1e-6)


def test_stiff_tendon_approaches_rigid_model():
    prev = elastic.EPS0_T
    try:
        elastic.EPS0_T = 0.0004
        _, a = elastic.find_equilibrium(TH)
    finally:
        elastic.EPS0_T = prev
    assert a == pytest.approx(linearize.find_equilibrium(TH), rel=2e-2)
