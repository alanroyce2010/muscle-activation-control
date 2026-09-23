"""Phase 6c: brachialis + lumped (lateral+medial) triceps antagonist pair.

Derivation and parameter provenance: docs/theory/phase6c_antagonist_pair.md.
State x = [theta, thetadot, a_b, a_t], input u = [u_b, u_t]. Brachialis
comes unchanged from joint.py; the triceps is built here with the same
muscle.py Thelen (2003) relations, a Murray (1995) digitized moment arm,
and a constant-width pennation model (alpha_0 = 9 deg).

Branch arguments (all default to the model's own switching logic):
  act_b, act_t: "activating"/"deactivating" per muscle
  vel: "flexing" (thetadot>0: brachialis shortening, triceps lengthening)
       or "extending" (the reverse) -- ONE switch for both muscles.
"""

from itertools import product

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.signal import cont2discrete

from muscle_activation_control import joint, muscle

# --- lumped triceps (TRIlat + TRImed), Holzbaur (2005) Table 1 ---
TRI_F0_M = 1248.6      # N  (624.3 + 624.3)
TRI_L0_M = 0.114       # m
TRI_LS_T = 0.0945      # m  (mean of 9.8 and 9.1 cm; only fixes the reference length)
TRI_ALPHA_0 = np.deg2rad(9.0)
TRI_W = TRI_L0_M * np.sin(TRI_ALPHA_0)   # constant fiber "width"

# --- triceps extension moment arm magnitude, Murray (1995) Fig. 4 Model panel ---
# pixel-level digitization (2026-09-11), 5-degree grid, magnitude; see
# data/murray1995_fig4_digitized.csv. Agrees with the original by-eye table
# [2.8, 3.0, 2.9, 2.75, 2.6, 2.4, 2.25, 2.2] cm within 0.13 cm everywhere.
_TRI_R_DEG = np.arange(0, 131, 5, dtype=float)
_TRI_R_M = np.array([2.886, 2.915, 2.964, 2.966, 2.925, 2.886, 2.809, 2.818, 2.818, 2.789, 2.775, 2.775, 2.714, 2.644, 2.615, 2.538, 2.470, 2.407, 2.305, 2.257, 2.310, 2.310, 2.257, 2.259, 2.210, 2.160, 2.160]) * 1e-2
# Holzbaur-consistent scaling (docs/theory/phase7_limitations.md section 2):
# Holzbaur (2005) ma_avg 2.1 cm / this curve's ROM mean. Applied only when
# a caller sets TRI_R_SCALE (sensitivity analysis); default 1 = Murray model.
TRI_R_SCALE = 1.0
_TRI_R_RAD = np.deg2rad(_TRI_R_DEG)


_TRI_PCHIP = PchipInterpolator(_TRI_R_RAD, _TRI_R_M)


def tri_moment_arm(theta):
    """|r_tri(theta)|, m (extension; the sign is applied in tri_torque)."""
    theta = np.clip(np.asarray(theta, dtype=float), _TRI_R_RAD[0], _TRI_R_RAD[-1])
    return TRI_R_SCALE * _TRI_PCHIP(theta)


# l_MT_tri lengthens with flexion: dl_MT/dtheta = +r_tri
_GRID = np.linspace(joint.THETA_MIN, joint.THETA_MAX, 2001)


def _rebuild_length_table():
    global _RG, _CUM, _DELTA_FROM_REF
    _RG = tri_moment_arm(_GRID)
    _CUM = np.concatenate(([0.0], np.cumsum(0.5 * (_RG[1:] + _RG[:-1]) * np.diff(_GRID))))
    _DELTA_FROM_REF = _CUM - np.interp(joint.THETA_REF, _GRID, _CUM)


_rebuild_length_table()


def set_tri_moment_arm_scale(scale):
    """Scale the triceps moment arm and rebuild l_MT(theta) consistently
    (r = dl_MT/dtheta must stay the derivative of the length table)."""
    global TRI_R_SCALE
    TRI_R_SCALE = float(scale)
    _rebuild_length_table()
L_MT_TRI_REF = TRI_LS_T + TRI_L0_M * np.cos(TRI_ALPHA_0)   # fiber at l_0 (pennated) at THETA_REF


def tri_musculotendon_length(theta):
    return L_MT_TRI_REF + np.interp(np.asarray(theta, dtype=float), _GRID, _DELTA_FROM_REF)


def tri_fiber_state(theta):
    """(l_bar, cos_alpha) from the constant-width rigid-tendon geometry."""
    proj = tri_musculotendon_length(theta) - TRI_LS_T
    l_m = np.sqrt(proj ** 2 + TRI_W ** 2)
    return l_m / TRI_L0_M, proj / l_m


def tri_fiber_velocity_norm(theta, thetadot):
    _, cos_a = tri_fiber_state(theta)
    return tri_moment_arm(theta) * np.asarray(thetadot, dtype=float) * cos_a / TRI_L0_M


def _tri_branch(vel):
    return None if vel is None else ("lengthening" if vel == "flexing" else "shortening")


def _bra_branch(vel):
    return None if vel is None else ("shortening" if vel == "flexing" else "lengthening")


def tri_torque(theta, thetadot, a, vel=None):
    """Extension torque magnitude of the lumped triceps, N*m (>= 0)."""
    l_bar, cos_a = tri_fiber_state(theta)
    v_norm = tri_fiber_velocity_norm(theta, thetadot)
    f_mt = muscle.muscle_tendon_force(a, l_bar, v_norm, TRI_F0_M, alpha=np.arccos(cos_a),
                                      velocity_branch=_tri_branch(vel))
    return f_mt * tri_moment_arm(theta)


def state_derivative(x, u, act_b=None, act_t=None, vel=None):
    theta, thetadot, a_b, a_t = x
    u_b, u_t = u
    tau = (joint.muscle_torque(theta, thetadot, a_b, velocity_branch=_bra_branch(vel))
           - tri_torque(theta, thetadot, a_t, vel) - joint.gravity_torque(theta))
    return np.array([thetadot, tau / joint.FOREARM_HAND_I_ELBOW,
                     muscle.activation_derivative(a_b, u_b, branch=act_b),
                     muscle.activation_derivative(a_t, u_t, branch=act_t)])


def find_equilibrium(theta_star, a_t_star):
    """a_b* balancing gravity plus the held triceps co-contraction a_t*."""
    def residual(a_b):
        return (joint.isometric_torque(theta_star, a_b) - tri_torque(theta_star, 0.0, a_t_star)
                - joint.gravity_torque(theta_star))
    return brentq(residual, 0.0, 1.0)


BRANCHES = list(product(["activating", "deactivating"], ["activating", "deactivating"],
                        ["flexing", "extending"]))


def linearize(x, u, act_b=None, act_t=None, vel=None, h=1e-6):
    """Central-FD (A 4x4, B 4x2, f0) on the branch-forced pair model, the
    same construction as linearize.linearize_general."""
    x = np.asarray(x, dtype=float)
    u = np.asarray(u, dtype=float)

    def f(xx, uu):
        return state_derivative(xx, uu, act_b, act_t, vel)

    f0 = f(x, u)
    a = np.zeros((4, 4))
    b = np.zeros((4, 2))
    for i in range(4):
        dx = np.zeros(4); dx[i] = h
        a[:, i] = (f(x + dx, u) - f(x - dx, u)) / (2 * h)
    for i in range(2):
        du = np.zeros(2); du[i] = h
        b[:, i] = (f(x, u + du) - f(x, u - du)) / (2 * h)
    return a, b, f0


def discretize_zoh(a, b, t_s):
    a_d, b_d, _, _, _ = cont2discrete((a, b, np.eye(4), np.zeros((4, 2))), t_s, method="zoh")
    return a_d, b_d


# --- Camlibel/Heemels/Schumacher continuity check, multi-input version ---
def _solve_continuity(delta_a, delta_b, c, d):
    """delta_A = e c^T, delta_b = e d^T with c having exactly one nonzero
    entry (as in bimodal._solve_continuity, extended to a d vector)."""
    nz = np.flatnonzero(c)
    if len(nz) != 1:
        raise ValueError("c must have exactly one nonzero entry")
    j = nz[0]
    e = delta_a[:, j] / c[j]
    res = max(np.max(np.abs(delta_a - np.outer(e, c))), np.max(np.abs(delta_b - np.outer(e, d))))
    return e, res


def verify_switches(theta_star, a_b_star, a_t_star):
    x = np.array([theta_star, 0.0, a_b_star, a_t_star])
    u = np.array([a_b_star, a_t_star])
    base = ("activating", "activating", "flexing")
    out = {}
    specs = {
        "bra_activation": (0, ("deactivating", "activating", "flexing"), np.array([0, 0, -1.0, 0]), np.array([1.0, 0])),
        "tri_activation": (1, ("activating", "deactivating", "flexing"), np.array([0, 0, 0, -1.0]), np.array([0, 1.0])),
        "velocity": (2, ("activating", "activating", "extending"), np.array([0, 1.0, 0, 0]), np.array([0.0, 0])),
    }
    a1, b1, _ = linearize(x, u, *base)
    for name, (_, other, c, d) in specs.items():
        a2, b2, _ = linearize(x, u, *other)
        e, res = _solve_continuity(a1 - a2, b1 - b2, c, d)
        out[name] = {"c": c, "d": d, "e": e, "residual": res}
    rows = {"bra_activation": 2, "tri_activation": 3, "velocity": 1}
    disjoint = all(
        np.all(np.abs(np.delete(out[n]["e"], rows[n])) < 1e-7 * max(1.0, np.abs(out[n]["e"][rows[n]])))
        for n in out)
    out["disjoint"] = disjoint
    return out


def damping_contributions(theta_star, a_b_star, a_t_star):
    """Closed-form shortening-side damping torque slopes D_b, D_t (N*m*s)
    and the predicted force-velocity ratio, extending-side over
    flexing-side damping, (2 D_b + D_t)/(D_b + 2 D_t) -- the same
    convention as Phase 2's lengthening/shortening = 2 for one muscle.
    docs/theory/phase6c_antagonist_pair.md section 2."""
    slope = 1.0 + 1.0 / muscle.A_F
    r_b = joint.moment_arm(theta_star)
    f_l_b = muscle.active_force_length(joint.fiber_length_norm(theta_star))
    k_b = (0.25 + 0.75 * a_b_star) * muscle.V_MAX_M
    d_b = r_b ** 2 * joint.F0_M * a_b_star * f_l_b * slope / (k_b * joint.L0_M)
    r_t = tri_moment_arm(theta_star)
    l_bar_t, cos_t = tri_fiber_state(theta_star)
    f_l_t = muscle.active_force_length(l_bar_t)
    k_t = (0.25 + 0.75 * a_t_star) * muscle.V_MAX_M
    d_t = r_t ** 2 * cos_t ** 2 * TRI_F0_M * a_t_star * f_l_t * slope / (k_t * TRI_L0_M)
    ratio = (2 * d_b + d_t) / (d_b + 2 * d_t)
    return d_b, d_t, ratio


def balanced_cocontraction(theta_star):
    """a_t* at which D_b = D_t (force-velocity kink cancels), if it exists in (0, 1)."""
    def g(a_t):
        a_b = find_equilibrium(theta_star, a_t)
        d_b, d_t, _ = damping_contributions(theta_star, a_b, a_t)
        return d_b - d_t
    # upper bound: the largest a_t* that still admits an equilibrium (a_b* < 1)
    for hi in (0.9, 0.7, 0.5, 0.3, 0.2):
        try:
            find_equilibrium(theta_star, hi)
            break
        except ValueError:
            continue
    else:
        return None
    try:
        return brentq(g, 1e-4, hi)
    except ValueError:
        return None
