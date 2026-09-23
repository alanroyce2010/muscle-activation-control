"""Phase 6c analogue of prediction.py for the antagonist pair: held
brachialis excitation step, triceps excitation held at its equilibrium."""

import numpy as np
from scipy.integrate import solve_ivp

from muscle_activation_control import antagonist


def true_trajectory(theta_star, a_b, a_t, delta_u_b, n_steps, t_s):
    u = np.array([a_b + delta_u_b, a_t])
    x0 = np.array([theta_star, 0.0, a_b, a_t])
    t_eval = np.arange(n_steps + 1) * t_s
    sol = solve_ivp(lambda t, x: antagonist.state_derivative(x, u), (0, n_steps * t_s), x0,
                    t_eval=t_eval, rtol=1e-10, atol=1e-12)
    return sol.y.T


def linear_prediction(theta_star, a_b, a_t, delta_u_b, n_steps, t_s, branch):
    x0 = np.array([theta_star, 0.0, a_b, a_t])
    u0 = np.array([a_b, a_t])
    a_c, b_c, _ = antagonist.linearize(x0, u0, *branch)
    a_d, b_d = antagonist.discretize_zoh(a_c, b_c, t_s)
    du = np.array([delta_u_b, 0.0])
    dx = np.zeros((n_steps + 1, 4))
    d = np.zeros(4)
    for k in range(n_steps):
        dx[k] = d
        d = a_d @ d + b_d @ du
    dx[n_steps] = d
    return x0 + dx


def matched_branch(delta_u_b):
    # triceps excitation is held exactly at a_t*, so it sits ON its switch;
    # "deactivating" is what the model's own comparison returns at equality
    if delta_u_b > 0:
        return ("activating", "deactivating", "flexing")
    return ("deactivating", "deactivating", "extending")


def error_fractions(theta_star, a_b, a_t, delta_u_b, n_steps, t_s):
    true = true_trajectory(theta_star, a_b, a_t, delta_u_b, n_steps, t_s)
    swing = abs(true[-1, 0] - theta_star)
    m = matched_branch(delta_u_b)
    flip = {"activating": "deactivating", "deactivating": "activating",
            "flexing": "extending", "extending": "flexing"}
    variants = {
        "matched": m,
        "wrong_bra_activation": (flip[m[0]], m[1], m[2]),
        "wrong_tri_activation": (m[0], flip[m[1]], m[2]),
        "wrong_velocity": (m[0], m[1], flip[m[2]]),
        "all_wrong": (flip[m[0]], flip[m[1]], flip[m[2]]),
    }
    out = {}
    for name, br in variants.items():
        pred = linear_prediction(theta_star, a_b, a_t, delta_u_b, n_steps, t_s, br)
        out[name] = abs(pred[-1, 0] - true[-1, 0]) / swing
    out["swing_deg"] = float(np.rad2deg(swing))
    sgn = 1.0 if delta_u_b > 0 else -1.0
    out["on_matched_bra_act"] = float(np.mean(sgn * (a_b + delta_u_b - true[1:, 2]) >= -1e-9))
    out["on_matched_vel"] = float(np.mean(sgn * true[1:, 1] >= -1e-9))
    out["a_t_drift"] = float(np.max(np.abs(true[:, 3] - a_t)))
    return out
