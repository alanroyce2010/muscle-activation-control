"""Every number the manuscripts quote about the single-muscle model, from
one place, so papers/poster can be re-synced after any model change
(written 2026-09-11 after the moment-arm correction). Prints plain text;
also writes results/headline_numbers.json."""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from muscle_activation_control import bimodal, joint, linearize, muscle, prediction  # noqa: E402

TH = np.deg2rad(60.0)
T_S = 0.010


def main():
    out = {}
    th = np.linspace(joint.THETA_MIN, joint.THETA_MAX, 2001)
    r = joint.moment_arm(th)
    tau = joint.isometric_torque(th, 1.0)
    out["moment_arm"] = dict(peak_cm=100 * r.max(), peak_deg=float(np.rad2deg(th[np.argmax(r)])),
                             avg_cm=100 * np.trapezoid(r, th) / (th[-1] - th[0]),
                             at_star_cm=100 * float(joint.moment_arm(TH)))
    out["torque"] = dict(peak_Nm=float(tau.max()), peak_deg=float(np.rad2deg(th[np.argmax(tau)])),
                         gravity_peak_Nm=float(joint.gravity_torque(np.pi / 2)))
    a_star = linearize.find_equilibrium(TH)
    rho_a = muscle.TAU_DEACT / (muscle.TAU_ACT * (0.5 + 1.5 * a_star) ** 2)
    act = bimodal.verify_activation_switch(TH, a_star)
    vel = bimodal.verify_velocity_switch(TH, a_star)
    eigs = []
    for br in linearize.BRANCH_COMBINATIONS:
        A, _ = linearize.linearize(TH, a_star, *br)
        eigs += list(np.real(np.linalg.eigvals(A)))
    cond = max(np.max(np.abs(linearize.discretize_zoh(*linearize.linearize(TH, a_star, *br), ts)[0]))
               for br in linearize.BRANCH_COMBINATIONS for ts in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1))
    grav = [linearize.find_equilibrium(np.deg2rad(d)) for d in range(10, 126, 5)]
    out["equilibrium"] = dict(a_star=a_star, rho_a=rho_a, e_act=list(act["e"]), e_vel=list(vel["e"]),
                              residual_act=act["residual"], residual_vel=vel["residual"],
                              eig_min=min(eigs), eig_max=max(eigs), max_abs_Ad=cond,
                              gravity_a_min=min(grav), gravity_a_max=max(grav),
                              a_unity=(np.sqrt(muscle.TAU_DEACT / muscle.TAU_ACT) - 0.5) / 1.5)
    tab = {}
    for du in (0.01, -0.01):
        for ms, n in ((50, 5), (100, 10), (200, 20), (500, 50)):
            ef = prediction.error_fractions(TH, a_star, du, n, T_S)
            tab[f"{du:+.2f}@{ms}"] = {k: (100 * v if k in ("matched", "wrong_activation", "wrong_velocity", "both_wrong") else v)
                                      for k, v in ef.items()}
    out["phase4"] = tab
    wrong = [tab[f"+0.01@{ms}"][k] for ms in (50, 100, 200, 500)
             for k in ("wrong_activation", "wrong_velocity", "both_wrong")]
    out["phase4_range_pct_activating"] = [min(wrong), max(wrong)]
    (ROOT / "results").mkdir(exist_ok=True)
    json.dump(out, open(ROOT / "results" / "headline_numbers.json", "w"), indent=1, default=float)
    m, t, e = out["moment_arm"], out["torque"], out["equilibrium"]
    print(f"moment arm: peak {m['peak_cm']:.2f} cm at {m['peak_deg']:.0f} deg, avg {m['avg_cm']:.2f} cm, r(theta*) {m['at_star_cm']:.2f} cm")
    print(f"torque: peak {t['peak_Nm']:.1f} N m at {t['peak_deg']:.0f} deg; gravity peak {t['gravity_peak_Nm']:.2f} N m")
    print(f"a*(theta*=60) = {e['a_star']:.4f}; rho_a = {e['rho_a']:.2f}; a_unity = {e['a_unity']:.3f}")
    print(f"e_act = {np.round(e['e_act'], 2)}  e_vel = {np.round(e['e_vel'], 2)}  residuals {e['residual_act']:.1e} {e['residual_vel']:.1e}")
    print(f"branch eigenvalues {e['eig_min']:.1f} .. {e['eig_max']:.1f} /s; max|A_d| over T_s 1-100 ms {e['max_abs_Ad']:.2f}")
    print(f"gravity-only a* over 10-125 deg: {e['gravity_a_min']:.3f} .. {e['gravity_a_max']:.3f}")
    print("phase4 (% of true response):  horizon  swing  matched  wrongact  wrongvel  both")
    for key, v in tab.items():
        print(f"   {key:>9s}  {v['swing_deg']:.3f} deg  {v['matched']:5.1f}  {v['wrong_activation']:6.1f}  {v['wrong_velocity']:6.1f}  {v['both_wrong']:6.1f}")
    print(f"activating-step wrong-branch range: {out['phase4_range_pct_activating'][0]:.1f} .. {out['phase4_range_pct_activating'][1]:.1f} %")


if __name__ == "__main__":
    main()
