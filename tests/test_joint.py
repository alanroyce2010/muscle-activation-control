"""Secondary Phase 1 validation target: CLAUDE.md section 6.3, second
block / docs/theory/single_joint_dynamics.md section 6.

These checks are largely "did I transcribe the digitized moment-arm
table and the fiber-length reference-angle convention correctly," not an
independent test of the muscle force model -- see the honesty note in
docs/theory/single_joint_dynamics.md section 6 about the residual
circularity here. The force-model test lives in test_muscle.py.
"""

import numpy as np
import pytest

from muscle_activation_control import joint


def test_moment_arm_endpoints_match_digitized_table():
    assert joint.moment_arm(0.0) == pytest.approx(0.00891, abs=1e-9)
    assert joint.moment_arm(joint.THETA_MAX) == pytest.approx(0.02881, abs=1e-9)


def test_moment_arm_monotonically_increasing():
    # Murray (1995): brachialis moment arm increases roughly monotonically
    # with flexion, unlike biceps/brachioradialis which peak mid-ROM.
    theta = np.linspace(joint.THETA_MIN, joint.THETA_MAX, 200)
    r = joint.moment_arm(theta)
    assert np.all(np.diff(r) >= -1e-9)


def test_moment_arm_peak_in_target_range():
    # CLAUDE.md 6.3: peak 2.0-3.5 cm, occurring above ~100 deg.
    theta = np.linspace(joint.THETA_MIN, joint.THETA_MAX, 200)
    r = joint.moment_arm(theta)
    peak_r = r.max()
    peak_theta = theta[np.argmax(r)]
    assert 0.020 <= peak_r <= 0.035
    assert np.rad2deg(peak_theta) > 100.0


def test_moment_arm_average_near_holzbaur_target():
    # Holzbaur (2005) ma_avg = 1.8 cm (a different model, with different
    # wrapping). The pixel-level Murray model-panel curve averages 2.18 cm
    # (docs/theory/phase7_limitations.md section 1); the band brackets both.
    theta = np.linspace(joint.THETA_MIN, joint.THETA_MAX, 2000)
    r = joint.moment_arm(theta)
    avg_r = np.trapezoid(r, theta) / (joint.THETA_MAX - joint.THETA_MIN)
    assert 0.017 <= avg_r <= 0.024


def test_fiber_length_at_reference_angle_is_optimal():
    l_bar = joint.fiber_length_norm(joint.THETA_REF)
    assert l_bar == pytest.approx(1.0, abs=1e-9)


def test_fiber_length_stays_in_plausible_range():
    # Fiber shouldn't stretch/shrink wildly over a 130 deg ROM for a
    # muscle whose moment arm (a few cm) is much smaller than its
    # excursion-integrated length changes -- sanity bound, not a precise
    # target.
    theta = np.linspace(joint.THETA_MIN, joint.THETA_MAX, 200)
    l_bar = joint.fiber_length_norm(theta)
    assert np.all(l_bar > 0.3)
    assert np.all(l_bar < 2.0)


def test_isometric_torque_is_nonnegative_and_bounded():
    theta = np.linspace(joint.THETA_MIN, joint.THETA_MAX, 200)
    tau = joint.isometric_torque(theta, a=1.0)
    assert np.all(tau >= 0.0)
    # Loose upper bound: F0_M * max possible passive+active scale (say 2x)
    # times the largest plausible moment arm (3.5 cm) -- just guards
    # against a units/sign bug blowing this up by orders of magnitude.
    assert np.all(tau < joint.F0_M * 2.0 * 0.035)


def test_isometric_torque_zero_activation_is_passive_only():
    theta = np.array([joint.THETA_REF])
    tau_active = joint.isometric_torque(theta, a=1.0)
    tau_passive = joint.isometric_torque(theta, a=0.0)
    # At the reference angle l_bar=1 (optimal length), passive force is
    # ~0 (F_PE(1)=0), so a=0 torque should be near zero there.
    assert tau_passive[0] == pytest.approx(0.0, abs=1e-6)
    assert tau_active[0] > tau_passive[0]


def test_gravity_torque_zero_at_full_extension():
    assert joint.gravity_torque(0.0) == pytest.approx(0.0, abs=1e-9)


def test_gravity_torque_positive_and_increasing_near_extension():
    # sin(theta) increasing for theta in [0, 90deg]
    t1 = joint.gravity_torque(np.deg2rad(30.0))
    t2 = joint.gravity_torque(np.deg2rad(60.0))
    assert 0.0 < t1 < t2
