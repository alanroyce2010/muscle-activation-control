"""Phase 2: equilibrium finding, branch-aware linearization, ZOH discretization.

See docs/theory/phase2_linearization.md for the full derivation --
read that first, especially section 3. Short version: this model has two
kinks (activation-ODE branch, force-velocity branch) that sit exactly on
every physically valid equilibrium, unlike the smooth pneumatic pressure
ODE in origami-arm-control. Naive central finite differences would
silently blend across a kink and return a meaningless slope, so the two
branch-sensitive Jacobian entries are computed analytically (verified
against sympy in scripts/_verify_branch_derivatives.py) instead of by
finite-differencing straight through x*.
"""

from itertools import product

import numpy as np
from scipy.optimize import brentq
from scipy.signal import cont2discrete

from muscle_activation_control import joint, muscle

FD_STEP = 1e-6  # finite-difference step for the two kink-free partials


def find_equilibrium(theta_star):
    """a* (and u*=a*) such that thetaddot=0 and adot=0 at (theta_star, 0, a*).

    Isometric equilibrium: a*_f_L(l_bar(theta_star))*F0_M*r(theta_star) =
    gravity_torque(theta_star). Solved generally via root-finding (not
    assuming l_bar(theta_star)=1 the way the closed-form in the theory
    doc's illustrative calc does), so this works for any theta_star.
    """

    def residual(a):
        return joint.isometric_torque(theta_star, a) - joint.gravity_torque(theta_star)

    a_star = brentq(residual, 0.0, 1.0)
    return a_star


def _mechanical_row(theta_star, a_star, velocity_branch):
    """[d(thetaddot)/d(theta), d(thetaddot)/d(thetadot), d(thetaddot)/d(a)]
    at (theta_star, thetadot=0, a_star). theta and a partials are safe
    with central finite differences (see docs/theory/phase2_linearization.md
    section 3); the thetadot partial is branch-dependent (Katz 1939
    asymmetry) and computed analytically.
    """
    h = FD_STEP

    def thetaddot_of(theta, thetadot, a):
        tau_m = joint.muscle_torque(theta, thetadot, a)
        tau_g = joint.gravity_torque(theta)
        return (tau_m - tau_g) / joint.FOREARM_HAND_I_ELBOW

    d_dtheta = (
        thetaddot_of(theta_star + h, 0.0, a_star)
        - thetaddot_of(theta_star - h, 0.0, a_star)
    ) / (2 * h)
    d_da = (
        thetaddot_of(theta_star, 0.0, a_star + h)
        - thetaddot_of(theta_star, 0.0, a_star - h)
    ) / (2 * h)

    r_star = joint.moment_arm(theta_star)
    l_bar_star = joint.fiber_length_norm(theta_star)
    f_l_star = muscle.active_force_length(l_bar_star)
    k_star = (0.25 + 0.75 * a_star) * muscle.V_MAX_M
    a_f_l = a_star * f_l_star

    if velocity_branch == "shortening":
        df_ce_dv = a_f_l * (1.0 + 1.0 / muscle.A_F) / k_star
    elif velocity_branch == "lengthening":
        df_ce_dv = a_f_l * (2.0 + 2.0 / muscle.A_F) / k_star
    else:
        raise ValueError(velocity_branch)

    dv_dthetadot = -r_star / joint.L0_M
    d_dthetadot = (r_star * joint.F0_M * df_ce_dv * dv_dthetadot) / joint.FOREARM_HAND_I_ELBOW

    return d_dtheta, d_dthetadot, d_da


def _activation_row(a_star, activation_branch):
    """[d(adot)/d(a), d(adot)/d(u)] at (a_star, u=a_star). Both branches
    agree adot=0 there, but the slope depends on which branch -- see
    docs/theory/phase2_linearization.md section 3, verified against
    sympy in scripts/_verify_branch_derivatives.py.
    """
    if activation_branch == "activating":
        base = 0.5 + 1.5 * a_star
        d_da = -1.0 / (muscle.TAU_ACT * base)
        d_du = 1.0 / (muscle.TAU_ACT * base)
    elif activation_branch == "deactivating":
        base = 0.5 + 1.5 * a_star
        d_da = -base / muscle.TAU_DEACT
        d_du = base / muscle.TAU_DEACT
    else:
        raise ValueError(activation_branch)
    return d_da, d_du


def linearize(theta_star, a_star, activation_branch, velocity_branch):
    """Continuous-time A (3x3), B (3x1) at (theta_star, 0, a_star, u=a_star)."""
    d_dtheta, d_dthetadot, d_da_mech = _mechanical_row(theta_star, a_star, velocity_branch)
    d_da_act, d_du_act = _activation_row(a_star, activation_branch)

    a_cont = np.array(
        [
            [0.0, 1.0, 0.0],
            [d_dtheta, d_dthetadot, d_da_mech],
            [0.0, 0.0, d_da_act],
        ]
    )
    b_cont = np.array([[0.0], [0.0], [d_du_act]])
    return a_cont, b_cont


def linearize_general(x, u, activation_branch=None, velocity_branch=None, h=FD_STEP):
    """Central finite-difference Jacobians (A 3x3, B 3x1) and the drift
    f0 = f(x, u) of joint.state_derivative at ANY (x, u), not just an
    equilibrium.

    With both branches forced, the branch-forced model is smooth in a
    neighbourhood of any point (including points *on* a switching
    manifold), so central differences are exact-to-O(h^2) for that
    branch; at an equilibrium this reproduces linearize() (tested in
    tests/test_mpc.py). With a branch left as None, the model's own
    comparison is used inside the difference stencil -- at a point on
    that switch the stencil straddles the kink and the corresponding
    Jacobian entry comes out as the AVERAGE of the two one-sided slopes.
    That is exactly the silent "blended" linearization this project
    warns against, made available deliberately so Phase 6b can measure
    what it costs in closed loop.
    """
    x = np.asarray(x, dtype=float)

    def f(xx, uu):
        return joint.state_derivative(xx, uu, activation_branch, velocity_branch)

    f0 = f(x, u)
    a_cont = np.zeros((3, 3))
    for i in range(3):
        dx = np.zeros(3)
        dx[i] = h
        a_cont[:, i] = (f(x + dx, u) - f(x - dx, u)) / (2 * h)
    b_cont = ((f(x, u + h) - f(x, u - h)) / (2 * h)).reshape(3, 1)
    return a_cont, b_cont, f0


def discretize_affine(a_cont, b_cont, f0, t_s):
    """ZOH discretization of the affine model xdot = A dx + B du + f0
    (f0 = drift at the linearization point, zero at an equilibrium):
    dx[k+1] = A_d dx[k] + B_d du[k] + c_d. Done by treating f0 as a
    second, constant input column so scipy's exact ZOH handles it."""
    b_aug = np.hstack([b_cont, np.asarray(f0, dtype=float).reshape(3, 1)])
    a_d, b_aug_d, _, _, _ = cont2discrete(
        (a_cont, b_aug, np.eye(3), np.zeros((3, 2))), t_s, method="zoh"
    )
    return a_d, b_aug_d[:, :1], b_aug_d[:, 1]


def discretize_zoh(a_cont, b_cont, t_s):
    """ZOH discretization, matching origami-arm-control's mpc_control.py
    pattern (scipy.signal.cont2discrete, method='zoh')."""
    c = np.eye(3)
    d = np.zeros((3, 1))
    a_d, b_d, _, _, _ = cont2discrete((a_cont, b_cont, c, d), t_s, method="zoh")
    return a_d, b_d


BRANCH_COMBINATIONS = list(
    product(["activating", "deactivating"], ["shortening", "lengthening"])
)


def conditioning_sweep(theta_star, t_s_values):
    """max|A_d| for every (activation_branch, velocity_branch, t_s)
    combination -- the diagnostic origami-arm-control's mpc_control.md
    section 2 used to catch the pneumatic pressure blowup.
    """
    a_star = find_equilibrium(theta_star)
    results = []
    for activation_branch, velocity_branch in BRANCH_COMBINATIONS:
        a_cont, b_cont = linearize(theta_star, a_star, activation_branch, velocity_branch)
        for t_s in t_s_values:
            a_d, b_d = discretize_zoh(a_cont, b_cont, t_s)
            results.append(
                {
                    "activation_branch": activation_branch,
                    "velocity_branch": velocity_branch,
                    "t_s": t_s,
                    "max_abs_a_d": float(np.max(np.abs(a_d))),
                    "a_cont": a_cont,
                    "a_d": a_d,
                }
            )
    return a_star, results
