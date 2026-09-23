"""Primary Phase 1 validation target: CLAUDE.md section 6.3, first block.

Fully self-contained checks (no literature dependency beyond the
equations already transcribed in CLAUDE.md 6.1) -- these are the ones
Phase 1 should actually gate on.
"""

import numpy as np
import pytest

from muscle_activation_control import muscle


def test_active_force_length_peak():
    assert muscle.active_force_length(1.0) == pytest.approx(1.0)


def test_active_force_length_symmetric_gaussian():
    f_low = muscle.active_force_length(0.5)
    f_high = muscle.active_force_length(1.5)
    expected = np.exp(-0.25 / 0.45)
    assert f_low == pytest.approx(expected, abs=1e-6)
    assert f_high == pytest.approx(expected, abs=1e-6)
    assert expected == pytest.approx(0.57375, abs=1e-4)


def test_passive_force_length_boundary():
    # l_bar = 1 + eps0_M = 1.6: F_PE = (e^kPE - 1)/(e^kPE - 1) = 1 exactly
    assert muscle.passive_force_length(1.6) == pytest.approx(1.0, abs=1e-9)


def test_passive_force_length_clipped_below_optimal():
    # Unclipped Eq. 3 goes negative for l_bar < 1; must be clipped to 0.
    assert muscle.passive_force_length(0.8) == 0.0
    assert muscle.passive_force_length(1.0) == pytest.approx(0.0, abs=1e-9)


def test_isometric_force_equals_activation_times_force_length():
    # v_norm=0 must give F_bar = a*f_l exactly, both branches agree there.
    for a in (0.3, 0.7, 1.0):
        for l_bar in (0.9, 1.0, 1.1):
            f_l = muscle.active_force_length(l_bar)
            f_ce = muscle.contractile_force(a, l_bar, 0.0)
            assert f_ce == pytest.approx(a * f_l, abs=1e-8)


def test_max_shortening_velocity_gives_zero_force():
    # At a=1, l_bar=1 (f_l=1): F->0 as v_norm -> -V_MAX_M.
    f_ce = muscle.contractile_force(1.0, 1.0, -muscle.V_MAX_M)
    assert f_ce == pytest.approx(0.0, abs=1e-6)


def test_beyond_max_shortening_velocity_clips_to_zero():
    # v_norm more negative than -V_MAX_M is unphysical; must clip, not go negative.
    f_ce = muscle.contractile_force(1.0, 1.0, -2.0 * muscle.V_MAX_M)
    assert f_ce == 0.0


def test_eccentric_plateau():
    # As v_norm -> +infinity, F_bar -> a*f_l*F_LEN (saturates, per the
    # closed-form derivation in muscle.py -- see its docstring).
    f_ce_large_v = muscle.contractile_force(1.0, 1.0, 1e6)
    assert f_ce_large_v == pytest.approx(muscle.F_LEN, abs=1e-3)


def test_activation_time_constants_at_extremes():
    # u > a (activating): tau = tau_act*(0.5+1.5a)
    assert muscle.activation_tau(a=0.0, u=1.0) == pytest.approx(
        muscle.TAU_ACT * 0.5
    )
    # u <= a (deactivating): tau = tau_deact/(0.5+1.5a)
    assert muscle.activation_tau(a=1.0, u=0.0) == pytest.approx(
        muscle.TAU_DEACT / 2.0
    )


def test_activation_derivative_signs():
    assert muscle.activation_derivative(a=0.0, u=1.0) > 0  # activating
    assert muscle.activation_derivative(a=1.0, u=0.0) < 0  # deactivating
    assert muscle.activation_derivative(a=0.5, u=0.5) == pytest.approx(0.0)
