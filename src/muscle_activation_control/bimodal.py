"""Phase 3: formalize the two equilibrium kinks (docs/theory/phase2_linearization.md
section 3) in Camlibel, Heemels, Schumacher (2008)'s bimodal piecewise-linear
system language.

Their system class (Automatica 44(5):1261-1267, Eqs. 1-3):
    xdot = A1*x + b1*u  if y<=0
    xdot = A2*x + b2*u  if y>0
    y = c^T*x + d*u
required continuous across {y=0}, equivalent to existence of e with
    A1 - A2 = e*c^T
    b1 - b2 = e*d

Our system has TWO independent switches, not their one -- this module
checks each individually against their continuity condition (each *is*
literally an instance of their bimodal system class on its own, verified
here algebraically, not just asserted), and separately checks whether the
two switches' effects are structurally independent (act on disjoint
matrix entries), which is what would let the 4-mode combination decompose
cleanly rather than requiring the full multimodal/"conewise" extension
their own paper's future-work section points to but which this project
has not pulled or verified (see docs/theory/phase3_bimodal_formalization.md).
"""

import numpy as np

from muscle_activation_control import linearize


def _solve_continuity(delta_a, delta_b, c, d, tol=1e-9):
    """Given delta_A = A1-A2, delta_b = b1-b2, and the switch's (c,d),
    solve for e such that delta_A = e*c^T and delta_b = e*d (Camlibel et
    al. Eq. 3), and return (e, max residual) -- residual should be ~0 if
    the switch genuinely satisfies their continuity condition.

    Assumes c has exactly one nonzero entry (true for both our switches),
    so e is solved directly rather than via a general least-squares fit --
    a general fit would silently "succeed" (small residual) even for a
    switch that isn't really of this form, by finding the best rank-one
    approximation; solving directly and checking the residual is the
    honest version of this check.
    """
    nonzero = np.flatnonzero(c)
    if len(nonzero) != 1:
        raise ValueError("this helper assumes c has exactly one nonzero entry")
    j = nonzero[0]
    e = delta_a[:, j] / c[j]
    residual_a = delta_a - np.outer(e, c)
    residual_b = delta_b.flatten() - e * d
    max_residual = max(np.max(np.abs(residual_a)), np.max(np.abs(residual_b)))
    return e, max_residual


def verify_activation_switch(theta_star, a_star):
    """Activation-ODE switch: y = u - a (c=[0,0,-1], d=1)."""
    a_act, b_act = linearize.linearize(theta_star, a_star, "activating", "shortening")
    a_deact, b_deact = linearize.linearize(theta_star, a_star, "deactivating", "shortening")
    c = np.array([0.0, 0.0, -1.0])
    d = 1.0
    e, residual = _solve_continuity(a_act - a_deact, b_act - b_deact, c, d)
    return {"c": c, "d": d, "e": e, "residual": residual}


def verify_velocity_switch(theta_star, a_star):
    """Force-velocity switch: y = thetadot (c=[0,1,0], d=0)."""
    a_len, b_len = linearize.linearize(theta_star, a_star, "activating", "lengthening")
    a_short, b_short = linearize.linearize(theta_star, a_star, "activating", "shortening")
    c = np.array([0.0, 1.0, 0.0])
    d = 0.0
    e, residual = _solve_continuity(a_len - a_short, b_len - b_short, c, d)
    return {"c": c, "d": d, "e": e, "residual": residual}


def verify_switch_independence(theta_star, a_star):
    """Check whether the two switches' effects on the Jacobian occupy
    disjoint entries -- if so, the 4-mode combination is exactly the
    direct sum of two independent bimodal switches (simpler, fully
    formalizable with the single paper this project has), rather than a
    genuinely coupled multimodal/conewise system (which would need the
    unverified follow-up paper, see the module docstring).
    """
    act_result = verify_activation_switch(theta_star, a_star)
    vel_result = verify_velocity_switch(theta_star, a_star)
    e_a = act_result["e"]
    e_v = vel_result["e"]
    # e_a should only touch the 'a' row (index 2), e_v only the 'thetadot' row (index 1)
    overlap = np.abs(e_a[1]) > 1e-9 or np.abs(e_v[2]) > 1e-9
    return {
        "e_activation": e_a,
        "e_velocity": e_v,
        "disjoint": not overlap,
    }
