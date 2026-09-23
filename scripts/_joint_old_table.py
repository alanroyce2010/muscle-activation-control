"""Single-joint elbow + brachialis dynamics.

Matches docs/theory/single_joint_dynamics.md exactly -- read that first.
Brachialis (BRA) parameters from Holzbaur et al. (2005) / arm26.osim
(CLAUDE.md 6.2). Moment arm r(theta) digitized from Murray et al. (1995)
Fig. 4 Model panel (docs/theory/single_joint_dynamics.md section 3).
Forearm+hand pendulum parameters from de Leva (1996) Table 4, combined
via the parallel axis theorem (docs/theory/single_joint_dynamics.md
section 2).
"""

import numpy as np

from muscle_activation_control import muscle

# --- Brachialis (BRA) parameters, Holzbaur (2005) / arm26.osim, CLAUDE.md 6.2 ---
F0_M = 987.26  # N, max isometric force
L0_M = 0.0858  # m, optimal fiber length
LS_T = 0.0535  # m, tendon slack length (unused once tendon is rigid)
ALPHA_0 = 0.0  # rad, pennation angle (BRA is unpennated)

THETA_MIN = 0.0  # rad, full extension
THETA_MAX = np.deg2rad(130.0)  # rad, full flexion
THETA_REF = np.deg2rad(65.0)  # rad, reference angle where l_M = L0_M (mid-ROM
# convention, docs/theory/single_joint_dynamics.md section 3)

# --- Moment arm r(theta), digitized from Murray (1995) Fig. 4 Model panel ---
# docs/theory/single_joint_dynamics.md section 3. Degrees -> meters.
_R_THETA_DEG = np.array([0, 20, 40, 60, 80, 100, 120, 130], dtype=float)
_R_THETA_M = np.array([0.6, 1.0, 1.4, 1.9, 2.3, 2.7, 3.0, 3.05], dtype=float) * 1e-2
_R_THETA_RAD = np.deg2rad(_R_THETA_DEG)


def moment_arm(theta):
    """r(theta), m. Piecewise-linear interpolation through the digitized
    Murray (1995) table -- deliberately not a fitted closed-form curve,
    see docs/theory/single_joint_dynamics.md section 3."""
    theta = np.asarray(theta, dtype=float)
    return np.interp(theta, _R_THETA_RAD, _R_THETA_M)


# Precompute l_MT(theta) - l_MT(THETA_REF) on a fine grid via cumulative
# trapezoidal integration of -r(theta), then interpolate. r(theta) = -dl_MT/dtheta,
# so l_MT(theta) = l_MT(THETA_REF) - integral_{THETA_REF}^{theta} r(phi) dphi.
_THETA_GRID = np.linspace(THETA_MIN, THETA_MAX, 2001)
_R_GRID = moment_arm(_THETA_GRID)
_NEG_R_CUMINT = np.concatenate(([0.0], np.cumsum(
    -0.5 * (_R_GRID[1:] + _R_GRID[:-1]) * np.diff(_THETA_GRID)
)))
# _NEG_R_CUMINT[i] = integral_{THETA_MIN}^{_THETA_GRID[i]} (-r(phi)) dphi
_DELTA_L_MT_FROM_REF = _NEG_R_CUMINT - np.interp(THETA_REF, _THETA_GRID, _NEG_R_CUMINT)

L_MT_REF = L0_M + LS_T  # l_MT(THETA_REF), from l_M(THETA_REF)=L0_M, rigid tendon


def musculotendon_length(theta):
    """l_MT(theta), m."""
    theta = np.asarray(theta, dtype=float)
    delta = np.interp(theta, _THETA_GRID, _DELTA_L_MT_FROM_REF)
    return L_MT_REF + delta


def fiber_length_norm(theta):
    """l_bar(theta) = l_M(theta)/L0_M, rigid tendon, alpha=0 (cos(alpha)=1)."""
    l_m = musculotendon_length(theta) - LS_T
    return l_m / L0_M


def fiber_velocity_norm(theta, thetadot):
    """v_norm = V^M/L0_M (optimal fiber lengths / s), from joint kinematics.

    docs/theory/single_joint_dynamics.md section 4:
    V^M = dl_M/dtheta * thetadot = -r(theta)*thetadot.
    """
    theta = np.asarray(theta, dtype=float)
    thetadot = np.asarray(thetadot, dtype=float)
    v_raw = -moment_arm(theta) * thetadot  # m/s
    return v_raw / L0_M


def muscle_torque(theta, thetadot, a, velocity_branch=None):
    """tau_muscle(theta,thetadot,a) = F^MT(theta,thetadot,a) * r(theta), N*m.
    velocity_branch: None (model's own comparison) or a forced branch,
    see muscle.contractile_force."""
    l_bar = fiber_length_norm(theta)
    v_norm = fiber_velocity_norm(theta, thetadot)
    f_mt = muscle.muscle_tendon_force(a, l_bar, v_norm, F0_M, alpha=ALPHA_0,
                                      velocity_branch=velocity_branch)
    return f_mt * moment_arm(theta)


# --- Forearm+hand pendulum, de Leva (1996) Table 4, combined at the elbow ---
# docs/theory/single_joint_dynamics.md section 2.
FOREARM_HAND_MASS = 1.628  # kg
FOREARM_HAND_LC = 0.1815  # m, elbow to combined COM
FOREARM_HAND_I_ELBOW = 0.0763  # kg*m^2, about elbow axis
G = 9.81  # m/s^2


def gravity_torque(theta):
    """tau_gravity(theta) = m*g*l_c*sin(theta), opposing flexion.

    Gravity configuration: upper arm fixed, hanging vertically at the
    side; theta=0 is forearm hanging straight down (stable, zero-torque
    equilibrium). docs/theory/single_joint_dynamics.md section 1.
    """
    theta = np.asarray(theta, dtype=float)
    return FOREARM_HAND_MASS * G * FOREARM_HAND_LC * np.sin(theta)


def state_derivative(state, u, activation_branch=None, velocity_branch=None):
    """dx/dt for x=[theta, thetadot, a], input u = neural excitation.

    docs/theory/single_joint_dynamics.md section 5. The two branch
    arguments default to the model's own switching logic; forcing them
    is only for linearizing on a chosen side of a switch (Phase 6).
    """
    theta, thetadot, a = state
    tau_m = muscle_torque(theta, thetadot, a, velocity_branch=velocity_branch)
    tau_g = gravity_torque(theta)
    thetaddot = (tau_m - tau_g) / FOREARM_HAND_I_ELBOW
    adot = muscle.activation_derivative(a, u, branch=activation_branch)
    return np.array([thetadot, thetaddot, adot])


def isometric_torque(theta, a=1.0):
    """Isometric joint torque (thetadot=0) at activation a, across theta.
    The Phase 1 validation quantity, CLAUDE.md 6.3 / docs/theory section 6."""
    return muscle_torque(theta, 0.0, a)
