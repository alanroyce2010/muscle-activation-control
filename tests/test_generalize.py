"""Phase 7 generalization checks (docs/theory/phase7_limitations.md
sections 2-3): the parameterized brachialis reproduces joint.py exactly;
biceps and brachioradialis plants have the same two-switch structure;
every Thelen constant set keeps rho_v = 2 exactly and rho_a in closed
form; rescaling the triceps moment arm keeps r = dl_MT/dtheta."""

import numpy as np
import pytest

from muscle_activation_control import antagonist, bimodal, generalize, joint, linearize, muscle

TH = np.deg2rad(60.0)


def test_bra_spec_reproduces_joint_model():
    p = generalize.Plant(generalize.PLANTS["BRA"]())
    a = p.equilibrium(TH)
    assert a == pytest.approx(linearize.find_equilibrium(TH), rel=1e-9)
    for act, vel in linearize.BRANCH_COMBINATIONS:
        jv = "flexing" if vel == "shortening" else "extending"
        A, B = p.linearize([TH, 0.0, a], a, act, jv)
        A0, B0 = linearize.linearize(TH, a, act, vel)
        assert np.allclose(A, A0, rtol=1e-5, atol=1e-6)
        assert np.allclose(B, B0, rtol=1e-5, atol=1e-6)


@pytest.mark.parametrize("name", ["BIC", "BRD"])
def test_other_flexors_have_the_same_switch_structure(name):
    p = generalize.Plant(generalize.PLANTS[name]())
    a = p.equilibrium(TH)
    x = [TH, 0.0, a]
    A_af, b_af = p.linearize(x, a, "activating", "flexing")
    A_df, b_df = p.linearize(x, a, "deactivating", "flexing")
    A_ae, b_ae = p.linearize(x, a, "activating", "extending")
    e_a, res_a = bimodal._solve_continuity(A_af - A_df, b_af - b_df, np.array([0, 0, -1.0]), 1.0)
    e_v, res_v = bimodal._solve_continuity(A_ae - A_af, b_ae - b_af, np.array([0, 1.0, 0]), 0.0)
    assert res_a < 1e-6 and res_v < 1e-6
    assert abs(e_a[0]) < 1e-9 and abs(e_a[1]) < 1e-9 and abs(e_a[2]) > 1.0     # activation row only
    assert abs(e_v[0]) < 1e-9 and abs(e_v[2]) < 1e-9 and abs(e_v[1]) > 0.1     # velocity row only
    assert A_ae[1, 1] / A_af[1, 1] == pytest.approx(2.0, rel=1e-4)


@pytest.mark.parametrize("name", ["young", "old", "opensim"])
def test_constant_sets_keep_rho_v_two_and_rho_a_closed_form(name):
    prev = muscle.use_constant_set(name)
    try:
        a = linearize.find_equilibrium(TH)
        A_as, _ = linearize.linearize(TH, a, "activating", "shortening")
        A_al, _ = linearize.linearize(TH, a, "activating", "lengthening")
        A_ds, _ = linearize.linearize(TH, a, "deactivating", "shortening")
        assert A_al[1, 1] / A_as[1, 1] == pytest.approx(2.0, rel=1e-9)
        rho = A_as[2, 2] / A_ds[2, 2]
        assert rho == pytest.approx(muscle.TAU_DEACT / (muscle.TAU_ACT * (0.5 + 1.5 * a) ** 2), rel=1e-12)
    finally:
        muscle.restore_constants(prev)
    assert muscle.TAU_ACT == 0.015 and muscle.F_LEN == 1.4


def test_triceps_rescaling_keeps_length_table_consistent():
    try:
        antagonist.set_tri_moment_arm_scale(0.81)
        th = np.deg2rad(np.linspace(20.0, 110.0, 10))
        h = 1e-5
        dl = (antagonist.tri_musculotendon_length(th + h) - antagonist.tri_musculotendon_length(th - h)) / (2 * h)
        assert np.allclose(dl, antagonist.tri_moment_arm(th), rtol=5e-3)
        assert antagonist.tri_musculotendon_length(joint.THETA_REF) == pytest.approx(antagonist.L_MT_TRI_REF)
    finally:
        antagonist.set_tri_moment_arm_scale(1.0)
    assert antagonist.TRI_R_SCALE == 1.0
