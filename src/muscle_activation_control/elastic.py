"""Phase 7 section 4: brachialis with an elastic tendon.

docs/theory/phase7_limitations.md section 4. Thelen (2003) Eq. 5 tendon
(eps0_T = 0.04, k_toe = 3, F_toe = 0.33, eps_toe = 0.609 eps0_T,
k_lin = 1.712/eps0_T), brachialis geometry and moment arm from joint.py,
zero pennation. State x = [theta, thetadot, l_bar_M, a], input u.

With an elastic tendon Thelen's force-velocity law runs in its published
direction: (theta, l_bar_M) fix the tendon force F_T, the contractile force
is F_ce = F_T - F_PE(l_bar_M), and Eq. 6-7 give the fibre velocity.
Every equilibrium needs dl_M/dt = 0, i.e. F_ce = a f_L: exactly the
boundary between Eq. 7's shortening and lengthening branches, now a fibre
switch with switching function y = F_ce(theta, l_M) - a f_L(l_M).

Branch arguments: act ("activating"/"deactivating") and fib
("shortening": y <= 0, "lengthening": y > 0); None = the model's own logic.
"""

import numpy as np
from scipy.optimize import brentq
from scipy.signal import cont2discrete

from muscle_activation_control import joint, muscle

EPS0_T = 0.04
K_TOE = 3.0
F_TOE = 0.33


def _eps_toe():
    return 0.609 * EPS0_T


def _k_lin():
    return 1.712 / EPS0_T


def tendon_force_norm(eps):
    """Thelen 2003 Eq. 5, normalized by F0; zero when slack (eps <= 0)."""
    eps = np.asarray(eps, dtype=float)
    et = _eps_toe()
    toe = F_TOE * (np.exp(K_TOE * eps / et) - 1.0) / (np.exp(K_TOE) - 1.0)
    lin = _k_lin() * (eps - et) + F_TOE
    return np.where(eps <= 0, 0.0, np.where(eps <= et, toe, lin))


def tendon_strain(f):
    """Inverse of tendon_force_norm for f > 0."""
    et = _eps_toe()
    if f <= F_TOE:
        return et / K_TOE * np.log(1.0 + f * (np.exp(K_TOE) - 1.0) / F_TOE)
    return et + (f - F_TOE) / _k_lin()


def forces(theta, l_bar):
    """(F_T, F_ce) normalized by F0, rigid-geometry l_MT(theta), alpha = 0."""
    l_t = joint.musculotendon_length(theta) - l_bar * joint.L0_M
    f_t = tendon_force_norm((l_t - joint.LS_T) / joint.LS_T)
    return f_t, f_t - muscle.passive_force_length(l_bar)


def switching_function(x):
    theta, _, l_bar, a = x
    _, f_ce = forces(theta, l_bar)
    return f_ce - a * muscle.active_force_length(l_bar)


def fibre_velocity(theta, l_bar, a, fib=None):
    """dl_bar/dt (optimal lengths per second), Thelen 2003 Eq. 6-7, with
    F_ce kept inside [0, 0.95 a f_L F_len] for the velocity formula (the
    formula is singular at the eccentric plateau; small-signal runs stay far
    from both limits)."""
    _, f_ce = forces(theta, l_bar)
    f_l = muscle.active_force_length(l_bar)
    afl = a * f_l
    f_ce = np.clip(f_ce, 0.0, 0.95 * afl * muscle.F_LEN)
    k = (0.25 + 0.75 * a) * muscle.V_MAX_M
    b_s = afl + f_ce / muscle.A_F
    b_l = (2.0 + 2.0 / muscle.A_F) * (afl * muscle.F_LEN - f_ce) / (muscle.F_LEN - 1.0)
    if fib is None:
        fib = "shortening" if f_ce <= afl else "lengthening"
    b = b_s if fib == "shortening" else b_l
    return k * (f_ce - afl) / b


def state_derivative(x, u, act=None, fib=None, tau_ext=0.0):
    theta, thetadot, l_bar, a = x
    f_t, _ = forces(theta, l_bar)
    tau = f_t * joint.F0_M * joint.moment_arm(theta) - joint.gravity_torque(theta) - tau_ext
    return np.array([thetadot, tau / joint.FOREARM_HAND_I_ELBOW,
                     fibre_velocity(theta, l_bar, a, fib),
                     muscle.activation_derivative(a, u, branch=act)])


def find_equilibrium(theta_star):
    """(l_bar*, a*) with zero velocity, zero fibre velocity, torque balance."""
    f_t = float(joint.gravity_torque(theta_star) / (joint.F0_M * joint.moment_arm(theta_star)))
    l_t = joint.LS_T * (1.0 + tendon_strain(f_t))
    l_bar = float((joint.musculotendon_length(theta_star) - l_t) / joint.L0_M)
    a = (f_t - float(muscle.passive_force_length(l_bar))) / float(muscle.active_force_length(l_bar))
    return l_bar, a


def linearize(x, u, act, fib, h=1e-7):
    """Branch-forced central differences (4x4 A, 4x1 B)."""
    x = np.asarray(x, dtype=float)
    f = lambda xx, uu: state_derivative(xx, uu, act, fib)
    A = np.zeros((4, 4))
    for i in range(4):
        dx = np.zeros(4); dx[i] = h
        A[:, i] = (f(x + dx, u) - f(x - dx, u)) / (2 * h)
    B = ((f(x, u + h) - f(x, u - h)) / (2 * h)).reshape(4, 1)
    return A, B


def switching_gradient(x, h=1e-7):
    x = np.asarray(x, dtype=float)
    g = np.zeros(4)
    for i in range(4):
        dx = np.zeros(4); dx[i] = h
        g[i] = (switching_function(x + dx) - switching_function(x - dx)) / (2 * h)
    return g


def continuity(delta_a, delta_b, c, d):
    """General Camlibel et al. check: e = (A1 - A2) v / (c^T v) for a v with
    c^T v != 0 (here v = c), then the residual of A1 - A2 = e c^T, b1 - b2 = e d."""
    v = c / np.dot(c, c)
    e = delta_a @ v
    res = max(np.max(np.abs(delta_a - np.outer(e, c))), np.max(np.abs(delta_b.ravel() - e * d)))
    return e, res


def discretize_zoh(A, B, t_s):
    a_d, b_d, _, _, _ = cont2discrete((A, B, np.eye(4), np.zeros((4, 1))), t_s, method="zoh")
    return a_d, b_d


BRANCHES = [(a, f) for a in ("activating", "deactivating") for f in ("shortening", "lengthening")]
