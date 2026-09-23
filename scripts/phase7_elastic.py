"""Phase 7 section 4 (docs/theory/phase7_limitations.md): elastic-tendon
brachialis at the representative posture. Structure (continuity,
disjoint support), eigenvalues per branch, ZOH conditioning, and the
Phase 4 linearize-and-hold protocol against the true nonlinear model.
Writes results/phase7_elastic.json."""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from muscle_activation_control import elastic, linearize  # noqa: E402

TH = np.deg2rad(60.0)
T_S = 0.010
HOR = {50: 5, 100: 10, 200: 20, 500: 50}


def main():
    l_bar, a = elastic.find_equilibrium(TH)
    x0 = np.array([TH, 0.0, l_bar, a])
    print(f"elastic-tendon equilibrium at 60 deg: l_bar* = {l_bar:.4f}, a* = {a:.4f} "
          f"(rigid: a* = {linearize.find_equilibrium(TH):.4f}); switching function y = {elastic.switching_function(x0):.1e}")
    mats = {br: elastic.linearize(x0, a, *br) for br in elastic.BRANCHES}
    A_as, b_as = mats[("activating", "shortening")]
    A_ds, b_ds = mats[("deactivating", "shortening")]
    A_al, b_al = mats[("activating", "lengthening")]
    e_a, r_a = elastic.continuity(A_as - A_ds, b_as - b_ds, np.array([0, 0, 0, -1.0]), 1.0)
    c_v = elastic.switching_gradient(x0)
    e_v, r_v = elastic.continuity(A_al - A_as, b_al - b_as, c_v, 0.0)
    ratio_fib = A_as[2, 2] / A_al[2, 2]
    print(f"activation switch: e = {np.round(e_a, 3)}, residual {r_a:.1e}")
    print(f"fibre switch: c = grad y = {np.round(c_v, 3)}, e = {np.round(e_v, 3)}, residual {r_v:.1e}")
    print(f"fibre-row slope ratio shortening/lengthening d(l_dot)/d(l_bar): {ratio_fib:.4f} (predicted 2)")
    eig = {f"{a_}/{f_}": np.sort(np.real(np.linalg.eigvals(m[0]))).tolist() for (a_, f_), m in mats.items()}
    for k, v in eig.items():
        print(f"  eigenvalues {k:24s}: {np.round(v, 2)}")
    cond = max(np.max(np.abs(elastic.discretize_zoh(*m, ts)[0]))
               for m in mats.values() for ts in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1))
    print(f"max|A_d| over all branches, T_s 1-100 ms: {cond:.2f}")
    tab = {}
    for du in (0.01, -0.01):
        u = a + du
        for ms, n in HOR.items():
            sol = solve_ivp(lambda t, x: elastic.state_derivative(x, u), (0, n * T_S), x0,
                            t_eval=np.arange(1, n + 1) * T_S, rtol=1e-10, atol=1e-12, method="LSODA")
            th_true = sol.y[0, -1]
            swing = abs(th_true - TH)
            y = np.array([elastic.switching_function(sol.y[:, i]) for i in range(sol.y.shape[1])])
            sgn = 1.0 if du > 0 else -1.0
            occ = (float(np.mean(sgn * (u - sol.y[3]) >= -1e-9)), float(np.mean(-sgn * y >= -1e-9)))
            act_m = "activating" if du > 0 else "deactivating"
            fib_m = "shortening" if du > 0 else "lengthening"
            flip = {"activating": "deactivating", "deactivating": "activating",
                    "shortening": "lengthening", "lengthening": "shortening"}
            cases = {"matched": (act_m, fib_m), "wrong_activation": (flip[act_m], fib_m),
                     "wrong_fibre": (act_m, flip[fib_m]), "both_wrong": (flip[act_m], flip[fib_m])}
            row = {"swing_deg": float(np.rad2deg(swing)), "occupancy": occ}
            for k, br in cases.items():
                Ad, Bd = elastic.discretize_zoh(*mats[br], T_S)
                d = np.zeros(4)
                for _ in range(n):
                    d = Ad @ d + Bd[:, 0] * du
                row[k] = 100 * abs(TH + d[0] - th_true) / swing
            tab[f"{du:+.2f}@{ms}"] = row
            print(f"  {du:+.2f} @ {ms:3d} ms: swing {row['swing_deg']:.3f} deg; matched {row['matched']:.1f}%, "
                  f"wrong act {row['wrong_activation']:.1f}%, wrong fibre {row['wrong_fibre']:.1f}%, "
                  f"both {row['both_wrong']:.1f}%; occupancy {occ}")
    wrong = [tab[f"+0.01@{ms}"][k] for ms in HOR for k in ("wrong_activation", "wrong_fibre", "both_wrong")]
    out = dict(l_bar=l_bar, a_star=a, e_act=e_a.tolist(), e_fib=e_v.tolist(), c_fib=c_v.tolist(),
               residuals=[r_a, r_v], fibre_ratio=ratio_fib, eigenvalues=eig, max_abs_Ad=cond,
               phase4=tab, range_activating_pct=[min(wrong), max(wrong)])
    print(f"activating-step wrong-branch range: {min(wrong):.1f}-{max(wrong):.1f}%")
    (ROOT / "results").mkdir(exist_ok=True)
    json.dump(out, open(ROOT / "results" / "phase7_elastic.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
