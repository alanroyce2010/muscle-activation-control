"""Phase 7 section 7 (docs/theory/phase7_limitations.md): why feedback
washes the branch ambiguity out. For the four branch models and the
blended (central-difference) model at the representative equilibrium,
ZOH at T_s = 10 ms, full-state output:

  - Vinnicombe nu-gap delta_nu(P_i, P_j) (pointwise chordal distance,
    maximized over the unit circle, winding condition checked);
  - generalized stability margin b(P_i, K_i) of the LQR gain designed on
    each model with the MPC's weights (Q = diag(1000, 1, 0), R = 200; the
    unconstrained MPC with a DARE terminal cost is this LQR);
  - the direct check: spectral radius of A_j - B_j K_i for every pair.

delta_nu(P_i, P_j) < b(P_i, K_i) certifies that K_i stabilizes P_j (a
statement about the linear branch models, not a proof for the
piecewise-affine plant). Writes results/phase7_nugap.json."""

import json
import sys
from pathlib import Path

import cvxpy as cp
import numpy as np
from scipy.linalg import solve_discrete_are

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from muscle_activation_control import linearize  # noqa: E402

TH = np.deg2rad(60.0)
T_S = 0.010
Q = np.diag([1000.0, 1.0, 0.0])
R = np.array([[200.0]])
OMEGA = np.concatenate([np.logspace(-7, np.log10(np.pi), 20000)])


def models():
    a = linearize.find_equilibrium(TH)
    out = {}
    for act, vel in linearize.BRANCH_COMBINATIONS:
        A, B = linearize.linearize(TH, a, act, vel)
        out[f"{act[:5]}/{vel[:5]}"] = linearize.discretize_zoh(A, B, T_S)
    A, B, _ = linearize.linearize_general([TH, 0.0, a], a, None, None)
    out["blended"] = linearize.discretize_zoh(A, B, T_S)
    return out


def freq(Ad, Bd, w):
    z = np.exp(1j * w)
    return np.array([np.linalg.solve(zi * np.eye(3) - Ad, Bd[:, 0]) for zi in z])   # (len(w), 3)


def chordal(p1, p2):
    """Pointwise nu-gap distance for single-input, full-state-output plants:
    sine of the angle between the graph lines span([p; 1]) in C^4."""
    g1 = np.concatenate([p1, [1.0]]); g2 = np.concatenate([p2, [1.0]])
    c = abs(np.vdot(g1, g2)) / (np.linalg.norm(g1) * np.linalg.norm(g2))
    return float(np.sqrt(max(0.0, 1.0 - c * c)))


def winding_ok(P1, P2):
    """wno(1 + P2^H P1) = 0 over the full circle (stable plants), and no zero on it."""
    w = np.concatenate([-OMEGA[::-1], OMEGA])
    f1 = np.concatenate([np.conj(P1[::-1]), P1]); f2 = np.concatenate([np.conj(P2[::-1]), P2])
    d = 1.0 + np.einsum("ij,ij->i", np.conj(f2), f1)
    ph = np.unwrap(np.angle(d))
    return bool(abs(ph[-1] - ph[0]) < np.pi) and bool(np.min(np.abs(d)) > 1e-9)


def lqr(Ad, Bd):
    P = solve_discrete_are(Ad, Bd, Q, R)
    return np.linalg.solve(R + Bd.T @ P @ Bd, Bd.T @ P @ Ad)     # u = -K x


def margin(Pw, K):
    """b(P, Kc) with positive-feedback controller Kc = -K (u = Kc y):
    min over w of |1 - Kc P| / (sqrt(1 + ||P||^2) sqrt(1 + ||Kc||^2))."""
    Kc = -K[0]
    num = np.abs(1.0 - Pw @ Kc)
    den = np.sqrt(1.0 + np.sum(np.abs(Pw) ** 2, axis=1)) * np.sqrt(1.0 + Kc @ Kc)
    return float(np.min(num / den))


def cqlf_rate(mats, iters=30):
    """Smallest gamma (bisection) for which a common quadratic Lyapunov
    function P > 0 with M^T P M <= gamma^2 P exists for every M in mats:
    a certificate of stability under ARBITRARY switching among the modes
    when gamma < 1 (upper bound on the joint spectral radius). Verified a
    posteriori with eigenvalues."""
    lo, hi, best = 0.0, 1.0, None
    for _ in range(iters):
        g = 0.5 * (lo + hi)
        P = cp.Variable((3, 3), symmetric=True)
        cons = [P >> np.eye(3)] + [M.T @ P @ M - g ** 2 * P << -1e-9 * np.eye(3) for M in mats]
        prob = cp.Problem(cp.Minimize(cp.trace(P)), cons)
        try:
            prob.solve(solver=cp.CLARABEL)
        except cp.error.SolverError:
            prob.solve(solver=cp.SCS)
        ok = prob.status in ("optimal", "optimal_inaccurate") and P.value is not None
        if ok:
            Pv = 0.5 * (P.value + P.value.T)
            ok = np.min(np.linalg.eigvalsh(Pv)) > 0 and all(
                np.max(np.linalg.eigvalsh(M.T @ Pv @ M - g ** 2 * Pv)) <= 1e-7 for M in mats)
        if ok:
            hi, best = g, g
        else:
            lo = g
    return best


def main():
    ms = models()
    names = list(ms)
    Pw = {n: freq(*ms[n], OMEGA) for n in names}
    K = {n: lqr(*ms[n]) for n in names}
    out = dict(nugap={}, margin={}, spectral_radius={}, certified={})
    print("generalized stability margin b(P_i, K_i):")
    for n in names:
        out["margin"][n] = margin(Pw[n], K[n])
        print(f"  {n:12s} b = {out['margin'][n]:.3f}")
    print("nu-gap delta(P_i, P_j) [winding ok] / spectral radius of A_j - B_j K_i:")
    for i in names:
        for j in names:
            if i == j:
                continue
            d = max(chordal(p1, p2) for p1, p2 in zip(Pw[i], Pw[j]))
            ok = winding_ok(Pw[i], Pw[j])
            Ad, Bd = ms[j]
            rho = float(max(abs(np.linalg.eigvals(Ad - Bd @ K[i]))))
            out["nugap"][f"{i}|{j}"] = dict(delta=d, winding_ok=ok)
            out["spectral_radius"][f"K_{i} on {j}"] = rho
            out["certified"][f"K_{i} on {j}"] = bool(ok and d < out["margin"][i])
            print(f"  {i:12s} -> {j:12s} delta = {d:.3f} [{'ok' if ok else 'FAIL'}]  "
                  f"{'<' if d < out['margin'][i] else '>='} b_i;  rho(A_j - B_j K_i) = {rho:.4f}")
    n_cert = sum(out["certified"].values()); n_stab = sum(v < 1 for v in out["spectral_radius"].values())
    print(f"certified by nu-gap: {n_cert}/{len(out['certified'])}; directly stable: {n_stab}/{len(out['spectral_radius'])}")
    branches = [n for n in names if n != "blended"]
    print("common quadratic Lyapunov function (arbitrary switching among the 4 branch plants):")
    out["cqlf"] = {}
    for i in names:
        g = cqlf_rate([ms[j][0] - ms[j][1] @ K[i] for j in branches])
        out["cqlf"][f"K_{i}"] = g
        print(f"  fixed design K_{i:12s}: gamma = {g if g is None else round(g, 4)}  -> {'CERTIFIED' if g is not None and g < 1 else 'not certified'}")
    g4 = cqlf_rate([ms[j][0] - ms[j][1] @ K[j] for j in branches])
    out["cqlf"]["sign test, correct classification (4 matched modes)"] = g4
    print(f"  sign-test rule, correct classification (4 matched modes): gamma = {g4 if g4 is None else round(g4, 4)}"
          f"  -> {'CERTIFIED' if g4 is not None and g4 < 1 else 'not certified'}")
    g = cqlf_rate([ms[j][0] - ms[j][1] @ K[i] for i in branches for j in branches])
    out["cqlf"]["sign test, any misclassification (16 modes)"] = g
    print(f"  sign-test rule, any design on any plant (16 modes): gamma = {g if g is None else round(g, 4)}"
          f"  -> {'CERTIFIED' if g is not None and g < 1 else 'not certified'}")
    (ROOT / "results").mkdir(exist_ok=True)
    json.dump(out, open(ROOT / "results" / "phase7_nugap.json", "w"), indent=1)


if __name__ == "__main__":
    main()
