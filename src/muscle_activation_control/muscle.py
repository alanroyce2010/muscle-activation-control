"""Thelen (2003) Hill-type muscle model, rigid-tendon variant.

Equations and constants match CLAUDE.md section 6.1 exactly (which in
turn transcribes Thelen, D.G. (2003), "Adjustment of muscle mechanics
model parameters to simulate dynamic contractions in older adults,"
J Biomech Eng 125(1):70-77, Appendix -- see
docs/literature/thelen2003_muscle_model.pdf). Young-adult constants
throughout; do not reuse for the paper's old-adult variant without
re-deriving.

All "l_bar" / "v_norm" quantities are normalized: length by optimal
fiber length (L0_M), velocity by optimal fiber lengths per second. Force
quantities returned by this module are normalized by F0_M (max isometric
force) unless a muscle-specific F0 is passed in, in which case the
output is in newtons.
"""

import numpy as np

TAU_ACT = 0.015  # s, young adult (Thelen 2003 Table 1)
TAU_DEACT = 0.050  # s, young adult

GAMMA = 0.45  # active force-length Gaussian shape factor
K_PE = 5.0  # passive force-length exponential shape factor
EPS0_M = 0.6  # passive muscle strain at F0_M

A_F = 0.25  # force-velocity shape factor, young adult
F_LEN = 1.4  # eccentric plateau (max normalized lengthening force), young adult
V_MAX_M = 10.0  # optimal fiber lengths / s, young adult

# Alternative constant sets (docs/theory/phase7_limitations.md section 3.3).
# Function defaults read the module constants at call time, so
# use_constant_set() switches every downstream module (joint, linearize,
# antagonist) in one place; "young" is the project default.
CONSTANT_SETS = {
    "young":   dict(TAU_ACT=0.015, TAU_DEACT=0.050, GAMMA=0.45, K_PE=5.0, EPS0_M=0.6,
                    A_F=0.25, F_LEN=1.4, V_MAX_M=10.0),     # Thelen 2003, young adult
    "old":     dict(TAU_ACT=0.015, TAU_DEACT=0.060, GAMMA=0.45, K_PE=5.0, EPS0_M=0.5,
                    A_F=0.25, F_LEN=1.8, V_MAX_M=8.0),      # Thelen 2003 Table 1, old adult
    "opensim": dict(TAU_ACT=0.010, TAU_DEACT=0.040, GAMMA=0.5, K_PE=4.0, EPS0_M=0.6,
                    A_F=0.3, F_LEN=1.8, V_MAX_M=10.0),      # OpenSim Thelen2003Muscle (John 2011)
}


def use_constant_set(name):
    """Install a named constant set module-wide; returns the previous values
    so callers can restore them (use in try/finally)."""
    g = globals()
    previous = {k: g[k] for k in CONSTANT_SETS["young"]}
    g.update(CONSTANT_SETS[name])
    return previous


def restore_constants(previous):
    globals().update(previous)


def activation_tau(a, u, tau_act=None, tau_deact=None, branch=None):
    """Activation/deactivation time constant tau_a(a,u), Thelen 2003 Eq. 2.

    branch=None uses the model's own comparison (u > a). branch=
    "activating"/"deactivating" forces that branch regardless of (u, a)
    -- used by Phase 6 to linearize on a chosen side of the u=a switch
    (linearize.linearize_general) and to build the antagonist-pair
    model; never changes the default behaviour."""
    tau_act = TAU_ACT if tau_act is None else tau_act
    tau_deact = TAU_DEACT if tau_deact is None else tau_deact
    a = np.asarray(a, dtype=float)
    u = np.asarray(u, dtype=float)
    tau_activating = tau_act * (0.5 + 1.5 * a)
    tau_deactivating = tau_deact / (0.5 + 1.5 * a)
    if branch is None:
        return np.where(u > a, tau_activating, tau_deactivating)
    if branch == "activating":
        return tau_activating
    if branch == "deactivating":
        return tau_deactivating
    raise ValueError(branch)


def activation_derivative(a, u, tau_act=None, tau_deact=None, branch=None):
    """da/dt, Thelen 2003 Eq. 1. See activation_tau for `branch`."""
    return (np.asarray(u, dtype=float) - np.asarray(a, dtype=float)) / activation_tau(
        a, u, tau_act, tau_deact, branch
    )


def active_force_length(l_bar, gamma=None):
    """f_l(l_bar), Thelen 2003 Eq. 4 -- normalized active force-length curve."""
    gamma = GAMMA if gamma is None else gamma
    l_bar = np.asarray(l_bar, dtype=float)
    return np.exp(-((l_bar - 1.0) ** 2) / gamma)


def passive_force_length(l_bar, k_pe=None, eps0_m=None):
    """F_PE(l_bar), Thelen 2003 Eq. 3 -- normalized passive force-length curve.

    Eq. 3 as published is a single exponential with no explicit l_bar<1
    case, but it goes negative there (unphysical -- passive tissue can't
    push). Clipped at zero for l_bar<1; this clip is a standard,
    physically-motivated implementation detail, not something Thelen
    (2003)'s text states explicitly.
    """
    k_pe = K_PE if k_pe is None else k_pe
    eps0_m = EPS0_M if eps0_m is None else eps0_m
    l_bar = np.asarray(l_bar, dtype=float)
    raw = (np.exp(k_pe * (l_bar - 1.0) / eps0_m) - 1.0) / (np.exp(k_pe) - 1.0)
    return np.maximum(raw, 0.0)


def contractile_force(a, l_bar, v_norm, a_f=None, f_len=None, v_max_m=None,
                      branch=None):
    """Normalized active (contractile-element) force F_bar_M, given
    activation, normalized fiber length, and normalized fiber velocity.

    Closed-form inversion of Thelen 2003 Eq. 6-7 (force-velocity, given
    as V^M(F_bar_M) there) for F_bar_M(V^M) -- derived by hand for this
    project (CLAUDE.md 6.1 / docs/theory/single_joint_dynamics.md flagged
    this inversion as possibly needing a numerical solve; it has a closed
    form on each branch, algebraically straightforward, so no root-finding
    is needed). v_norm <= 0 is shortening (concentric), v_norm > 0 is
    lengthening (eccentric) -- sign convention matches
    docs/theory/single_joint_dynamics.md section 4 (V^M = -r(theta)*thetadot).

    Clipped at zero: for v_norm more negative than the theoretical max
    shortening velocity (-V_max_M at a=1), the unclipped algebra gives a
    negative force, which is unphysical (concentric contraction can't
    produce negative -- i.e. resistive -- active force in this model).

    branch=None uses the model's own comparison (v_norm <= 0);
    "shortening"/"lengthening" forces that branch (see activation_tau).
    """
    a_f = A_F if a_f is None else a_f
    f_len = F_LEN if f_len is None else f_len
    v_max_m = V_MAX_M if v_max_m is None else v_max_m
    a = np.asarray(a, dtype=float)
    l_bar = np.asarray(l_bar, dtype=float)
    v_norm = np.asarray(v_norm, dtype=float)

    f_l = active_force_length(l_bar)
    a_f_l = a * f_l
    k = (0.25 + 0.75 * a) * v_max_m

    f_shortening = a_f_l * (v_norm + k) / (k - v_norm / a_f)

    c = v_norm * (2.0 + 2.0 / a_f) / (f_len - 1.0)
    f_lengthening = a_f_l * (c * f_len + k) / (k + c)

    if branch is None:
        f = np.where(v_norm <= 0, f_shortening, f_lengthening)
    elif branch == "shortening":
        f = f_shortening
    elif branch == "lengthening":
        f = f_lengthening
    else:
        raise ValueError(branch)
    return np.maximum(f, 0.0)


def muscle_tendon_force(a, l_bar, v_norm, f0_m, alpha=0.0, velocity_branch=None):
    """Rigid-tendon muscle-tendon force F^MT, in the same units as f0_m
    (pass f0_m in newtons for a force in newtons, or 1.0 for normalized).

    CLAUDE.md 6.1: F^MT = (F_PE(l_bar) + a*f_L*f_V) * F0_M * cos(alpha) --
    here `contractile_force` returns the combined a*f_L*f_V term directly
    (see its docstring), so this is F^MT = (F_PE + F_bar_ce) * F0_M * cos(alpha).
    """
    f_pe = passive_force_length(l_bar)
    f_ce = contractile_force(a, l_bar, v_norm, branch=velocity_branch)
    return (f_pe + f_ce) * f0_m * np.cos(alpha)
