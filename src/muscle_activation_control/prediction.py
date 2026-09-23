"""One-shot branch-linearized prediction vs. true nonlinear evolution,
around an equilibrium (theta*, 0, a*) held by gravity plus an optional
constant external torque.

Generalizes the helpers in scripts/phase4_quantify_branch_consequence.py
(which are left untouched as the Phase 4 record) to arbitrary theta*,
a*, and a constant load torque tau_ext -- see
docs/theory/phase6a_equilibrium_sweep.md section 1.3 for why a constant
torque leaves linearize.linearize() valid unchanged.
"""

import numpy as np
from scipy.integrate import solve_ivp

from muscle_activation_control import joint, linearize


def external_torque_for_equilibrium(theta_star, a_star):
    """Constant tau_ext (N*m, positive = flexion) that makes
    (theta_star, 0, a_star) an equilibrium of the loaded joint."""
    # thetaddot = (tau_m - tau_g - tau_ext)/I = 0  =>  tau_ext = tau_m - tau_g
    return joint.isometric_torque(theta_star, a_star) - joint.gravity_torque(theta_star)


def true_trajectory(theta_star, a_star, delta_u, n_steps, t_s, tau_ext=0.0):
    """Integrate the true nonlinear model from the equilibrium under a
    held excitation step u = a* + delta_u, with constant load tau_ext
    opposing flexion. Returns (n_steps+1, 3) states sampled every t_s."""
    u = a_star + delta_u
    x0 = np.array([theta_star, 0.0, a_star])
    t_eval = np.arange(n_steps + 1) * t_s

    def rhs(t, x):
        dx = joint.state_derivative(x, u)
        dx[1] -= tau_ext / joint.FOREARM_HAND_I_ELBOW
        return dx

    sol = solve_ivp(rhs, (0, n_steps * t_s), x0, t_eval=t_eval, rtol=1e-10, atol=1e-12)
    return sol.y.T


def linear_prediction(theta_star, a_star, delta_u, n_steps, t_s,
                      activation_branch, velocity_branch):
    """Propagate the fixed ZOH branch-linearization forward n_steps under
    the same held excitation step. Returns (n_steps+1, 3) states."""
    a_cont, b_cont = linearize.linearize(theta_star, a_star, activation_branch, velocity_branch)
    a_d, b_d = linearize.discretize_zoh(a_cont, b_cont, t_s)
    b_d = b_d.flatten()
    dx = np.zeros((n_steps + 1, 3))
    delta_x = np.zeros(3)
    for k in range(n_steps):
        dx[k] = delta_x
        delta_x = a_d @ delta_x + b_d * delta_u
    dx[n_steps] = delta_x
    return np.array([theta_star, 0.0, a_star]) + dx


def matched_branch(delta_u):
    """Branch pair the true trajectory sits on for a held excitation step
    (activating + shortening for delta_u > 0, the mirror for < 0), as
    long as it does not cross a switch during the horizon -- callers
    should check that with branch_occupancy()."""
    if delta_u > 0:
        return ("activating", "shortening")
    return ("deactivating", "lengthening")


def other_branches(matched):
    act, vel = matched
    wrong_act = ("deactivating" if act == "activating" else "activating", vel)
    wrong_vel = (act, "lengthening" if vel == "shortening" else "shortening")
    return {"wrong_activation": wrong_act, "wrong_velocity": wrong_vel,
            "both_wrong": (wrong_act[0], wrong_vel[1])}


def branch_occupancy(true_traj, a_star, delta_u, tol=1e-9):
    """Fraction of samples (excluding t=0, which sits exactly on both
    switches) on which the true trajectory is on the *matched* side of
    each switch: activating (u > a) and shortening (thetadot > 0) for
    delta_u > 0, the mirror for delta_u < 0.

    Tolerance: under a held excitation step, a(t) converges onto u
    asymptotically, so u - a reaches integrator precision (~1e-12) within
    a few time constants; a strict comparison would then report a
    spurious "crossing" that is really the trajectory settling onto the
    switching manifold. A sample counts as a genuine crossing only if it
    is on the wrong side by more than tol."""
    u = a_star + delta_u
    a = true_traj[1:, 2]
    thetadot = true_traj[1:, 1]
    sgn = 1.0 if delta_u > 0 else -1.0
    on_act = np.mean(sgn * (u - a) >= -tol)
    on_short = np.mean(sgn * thetadot >= -tol)
    return float(on_act), float(on_short)


def error_fractions(theta_star, a_star, delta_u, n_steps, t_s, tau_ext=0.0):
    """Final-time theta prediction error of each branch pair, as a fraction
    of the true swing |theta_true(T) - theta*|. Also returns the branch
    occupancy so the caller can tell whether 'matched' really was."""
    true_traj = true_trajectory(theta_star, a_star, delta_u, n_steps, t_s, tau_ext)
    swing = abs(true_traj[-1, 0] - theta_star)
    matched = matched_branch(delta_u)
    branches = {"matched": matched, **other_branches(matched)}
    out = {}
    for name, br in branches.items():
        pred = linear_prediction(theta_star, a_star, delta_u, n_steps, t_s, *br)
        out[name] = abs(pred[-1, 0] - true_traj[-1, 0]) / swing
    occ_act, occ_short = branch_occupancy(true_traj, a_star, delta_u)
    out["swing_deg"] = float(np.rad2deg(swing))
    out["on_matched_activation_frac"] = occ_act
    out["on_matched_velocity_frac"] = occ_short
    return out
