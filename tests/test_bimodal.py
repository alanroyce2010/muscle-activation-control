"""Phase 3 check: do the two equilibrium kinks (docs/theory/
phase2_linearization.md section 3) actually satisfy Camlibel, Heemels,
Schumacher (2008)'s continuity condition for bimodal piecewise-linear
systems, individually? Verified algebraically here, not just claimed in
prose -- see src/muscle_activation_control/bimodal.py docstring.
"""

import numpy as np
import pytest

from muscle_activation_control import bimodal, linearize

THETA_STAR = np.deg2rad(60.0)


def test_activation_switch_satisfies_continuity_condition():
    a_star = linearize.find_equilibrium(THETA_STAR)
    result = bimodal.verify_activation_switch(THETA_STAR, a_star)
    assert result["residual"] < 1e-9


def test_velocity_switch_satisfies_continuity_condition():
    a_star = linearize.find_equilibrium(THETA_STAR)
    result = bimodal.verify_velocity_switch(THETA_STAR, a_star)
    assert result["residual"] < 1e-9


def test_activation_switch_e_only_touches_activation_row():
    a_star = linearize.find_equilibrium(THETA_STAR)
    result = bimodal.verify_activation_switch(THETA_STAR, a_star)
    e = result["e"]
    assert e[0] == pytest.approx(0.0, abs=1e-9)
    assert e[1] == pytest.approx(0.0, abs=1e-9)
    assert e[2] != 0.0


def test_velocity_switch_e_only_touches_thetadot_row():
    a_star = linearize.find_equilibrium(THETA_STAR)
    result = bimodal.verify_velocity_switch(THETA_STAR, a_star)
    e = result["e"]
    assert e[0] == pytest.approx(0.0, abs=1e-9)
    assert e[1] != 0.0
    assert e[2] == pytest.approx(0.0, abs=1e-9)


def test_the_two_switches_are_structurally_independent():
    a_star = linearize.find_equilibrium(THETA_STAR)
    result = bimodal.verify_switch_independence(THETA_STAR, a_star)
    assert result["disjoint"]


def test_velocity_switch_e_matches_katz_asymmetry():
    # Katz 1939 (docs/theory/phase2_linearization.md section 3): the
    # lengthening slope is 2x the shortening slope, so e_v (=A_len-A_short
    # at the thetadot/thetadot entry) should equal the shortening-branch
    # slope itself (since 2x - 1x = 1x).
    a_star = linearize.find_equilibrium(THETA_STAR)
    a_short, _ = linearize.linearize(THETA_STAR, a_star, "activating", "shortening")
    result = bimodal.verify_velocity_switch(THETA_STAR, a_star)
    assert result["e"][1] == pytest.approx(a_short[1, 1], rel=1e-6)
