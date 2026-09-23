"""Phase 7 sections 2-3 (docs/theory/phase7_limitations.md): triceps
moment-arm sensitivity, other elbow flexors, group torque validation
against measured maximum isometric moments (Holzbaur 2005 Fig. 5A), and
alternative Thelen constant sets. Writes results/phase7_generalization.json
and results/phase7_torque_validation.{png,pdf}."""

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from muscle_activation_control import (antagonist, bimodal, generalize, joint, linearize,  # noqa: E402
                                       muscle, prediction, prediction_pair)

TH = np.deg2rad(60.0)
T_S = 0.010
HOR = {50: 5, 100: 10, 200: 20, 500: 50}
RES = ROOT / "results"


def rho_a(a):
    return muscle.TAU_DEACT / (muscle.TAU_ACT * (0.5 + 1.5 * a) ** 2)


def section2_triceps():
    print("### section 2: triceps moment arm ###")
    g = np.linspace(joint.THETA_MIN, joint.THETA_MAX, 2001)
    antagonist.set_tri_moment_arm_scale(1.0)
    mean_cm = 100 * np.trapezoid(antagonist.tri_moment_arm(g), g) / (g[-1] - g[0])
    s_h = 2.1 / mean_cm
    r90 = 100 * float(antagonist.tri_moment_arm(np.pi / 2))
    print(f"pixel Murray TRI ROM mean {mean_cm:.2f} cm, r(90) {r90:.2f} cm; Holzbaur-consistent scale s = 2.1/{mean_cm:.2f} = {s_h:.3f}"
          f" -> r(90) {s_h * r90:.2f} cm (Holzbaur Fig. 3A model point ~2.0 cm)")
    out = dict(rom_mean_cm=mean_cm, r90_cm=r90, s_holzbaur=s_h, balance={})
    for s in (s_h, 0.9, 1.0):
        antagonist.set_tri_moment_arm_scale(s)
        # scan the feasible co-contraction range (a_b* <= 1) for the ratio curve
        a_ts, ratios = [], []
        for a_t in np.linspace(0.0, 0.95, 96):
            try:
                a_b = antagonist.find_equilibrium(TH, a_t)
            except ValueError:
                break
            a_ts.append(a_t); ratios.append(antagonist.damping_contributions(TH, a_b, a_t)[2])
        a_ts, ratios = np.array(a_ts), np.array(ratios)
        rec = dict(a_t_max_feasible=float(a_ts[-1]), ratio_min=float(ratios.min()),
                   ratio_min_at=float(a_ts[np.argmin(ratios)]), balance=None)
        a_bal = antagonist.balanced_cocontraction(TH)
        if a_bal is not None:
            a_b = antagonist.find_equilibrium(TH, a_bal)
            errs = {ms: prediction_pair.error_fractions(TH, a_b, a_bal, 0.01, n, T_S) for ms, n in ((50, 5), (500, 50))}
            rec["balance"] = dict(a_t=a_bal, a_b=a_b, wrong_vel_50=100 * errs[50]["wrong_velocity"],
                                  wrong_vel_500=100 * errs[500]["wrong_velocity"],
                                  all_wrong_50=100 * errs[50]["all_wrong"])
            print(f"  s = {s:.3f}: balance a_t* = {a_bal:.4f}, a_b* = {a_b:.4f}; wrong-velocity error at balance "
                  f"{100 * errs[50]['wrong_velocity']:.1f}% (50 ms), {100 * errs[500]['wrong_velocity']:.1f}% (500 ms)")
        else:
            print(f"  s = {s:.3f}: NO balance within feasible co-contraction (a_t up to {a_ts[-1]:.2f}); "
                  f"ratio falls from 2 to a minimum of {ratios.min():.3f} at a_t = {a_ts[np.argmin(ratios)]:.2f}")
        out["balance"][f"{s:.3f}"] = rec
    antagonist.set_tri_moment_arm_scale(1.0)
    return out


def section31_flexors():
    print("### section 3.1: other elbow flexors, single-muscle plants at 60 deg ###")
    out = {}
    for name, make in generalize.PLANTS.items():
        p = generalize.Plant(make())
        a = p.equilibrium(TH)
        x = [TH, 0.0, a]
        A_af, b_af = p.linearize(x, a, "activating", "flexing")
        A_df, b_df = p.linearize(x, a, "deactivating", "flexing")
        A_ae, b_ae = p.linearize(x, a, "activating", "extending")
        e_a, res_a = bimodal._solve_continuity(A_af - A_df, b_af - b_df, np.array([0, 0, -1.0]), 1.0)
        e_v, res_v = bimodal._solve_continuity(A_ae - A_af, b_ae - b_af, np.array([0, 1.0, 0]), 0.0)
        eig = max(np.max(np.real(np.linalg.eigvals(p.linearize(x, a, ac, jv)[0])))
                  for ac in ("activating", "deactivating") for jv in ("flexing", "extending"))
        errs = {}
        for du in (0.01, -0.01):
            for ms, n in HOR.items():
                errs[f"{du:+.2f}@{ms}"] = p.error_fractions(TH, du, n, T_S)
        wrong = [100 * errs[f"+0.01@{ms}"][k] for ms in HOR for k in ("wrong_activation", "wrong_velocity", "both_wrong")]
        occ_ok = all(min(v["occupancy"]) == 1.0 for v in errs.values())
        out[name] = dict(a_star=a, rho_a=rho_a(a), rho_v=A_ae[1, 1] / A_af[1, 1],
                         e_act=list(e_a), e_vel=list(e_v), residual=max(res_a, res_v),
                         slowest_eig=eig, capacity_60_Nm=float(p.muscle_torque(TH, 0.0, 1.0)),
                         range_activating_pct=[min(wrong), max(wrong)], occupancy_ok=occ_ok,
                         errors={k: {kk: (100 * vv if kk in ("matched", "wrong_activation", "wrong_velocity", "both_wrong") else vv)
                                     for kk, vv in v.items()} for k, v in errs.items()})
        e50 = errs["+0.01@50"]
        print(f"  {name}: a* = {a:.4f}, rho_a = {rho_a(a):.2f}, rho_v = {A_ae[1, 1] / A_af[1, 1]:.4f}, residual {max(res_a, res_v):.1e}, "
              f"slowest eig {eig:+.2f}/s; 50 ms: wrong act {100 * e50['wrong_activation']:.1f}%, "
              f"wrong vel {100 * e50['wrong_velocity']:.1f}%, both {100 * e50['both_wrong']:.1f}%; "
              f"activating range {min(wrong):.1f}-{max(wrong):.1f}%; occupancy ok: {occ_ok}")
    return out


def load_fig5a():
    rows = [r for r in csv.reader(open(ROOT / "data" / "holzbaur2005_fig5a_digitized.csv"))
            if r and not r[0].startswith("#") and r[0] != "source"]
    data = {}
    for src, direction, deg, nm in rows:
        data.setdefault((src, direction), []).append((float(deg), float(nm)))
    return {k: np.array(v) for k, v in data.items()}


def section32_validation():
    print("### section 3.2: group torque validation (Holzbaur 2005 Fig. 5A measured data) ###")
    meas = load_fig5a()
    th = np.deg2rad(np.linspace(0, 130, 131))
    deg = np.rad2deg(th)
    parts = {"BRA": [generalize.bra_spec()], "BIC": generalize.bic_specs(), "BRD": [generalize.brd_spec()]}
    flex = {k: generalize.group_isometric_torque(v, th) for k, v in parts.items()}
    flex_total = sum(flex.values())
    s_h = 2.1 / (100 * np.trapezoid(antagonist.tri_moment_arm(th), th) / (th[-1] - th[0]))
    ext = {"s=1": -generalize.group_isometric_torque(generalize.tri_specs(1.0), th),
           f"s={s_h:.2f}": -generalize.group_isometric_torque(generalize.tri_specs(s_h), th)}
    out = dict(flexion=dict(peak_Nm=float(flex_total.max()), peak_deg=float(deg[np.argmax(flex_total)]),
                            parts_peak={k: float(v.max()) for k, v in flex.items()}),
               extension={k: dict(peak_Nm=float(v.max()), peak_deg=float(deg[np.argmax(v)])) for k, v in ext.items()},
               comparison={})
    print(f"  flexor group (BRA+BIC+BRD, a=1) peak {flex_total.max():.1f} N m at {deg[np.argmax(flex_total)]:.0f} deg "
          f"(parts: " + ", ".join(f"{k} {v.max():.1f}" for k, v in flex.items()) + ")")
    for k, v in ext.items():
        print(f"  triceps (3 heads, a=1, {k}) peak {v.max():.1f} N m at {deg[np.argmax(v)]:.0f} deg")
    for (src, direction), arr in meas.items():
        model = np.interp(arr[:, 0], deg, flex_total) if direction == "flexion" else None
        if direction == "flexion":
            rms = float(np.sqrt(np.mean((model - arr[:, 1]) ** 2)))
            out["comparison"][f"{src}/{direction}"] = dict(rms_Nm=rms, measured_peak=float(arr[:, 1].max()))
            print(f"  vs {src} {direction}: RMS difference {rms:.1f} N m (measured peak {arr[:, 1].max():.0f} N m)")
        else:
            for k, v in ext.items():
                mv = np.interp(arr[:, 0], deg, v)
                rms = float(np.sqrt(np.mean((mv - np.abs(arr[:, 1])) ** 2)))
                out["comparison"][f"{src}/{direction}/{k}"] = dict(rms_Nm=rms, measured_peak=float(np.abs(arr[:, 1]).max()))
                print(f"  vs {src} {direction} ({k}): RMS difference {rms:.1f} N m (measured peak {np.abs(arr[:, 1]).max():.0f} N m)")
    # figure, plot theme
    plt.style.use(str(ROOT / "scripts" / "icra.mplstyle"))
    fig, axes = plt.subplots(1, 2, figsize=(7.16, 2.6))
    ax = axes[0]
    ax.plot(deg, flex_total, color="#0072B2", label="model group (BRA+BIC+BRD)")
    for (k, v), ls in zip(flex.items(), ("--", "-.", ":")):
        ax.plot(deg, v, color="#999999", ls=ls, lw=0.9, label=f"model {k}")
    for src, mk in (("Buchanan", "s"), ("Amis", "^")):
        arr = meas[(src, "flexion")]
        ax.plot(arr[:, 0], arr[:, 1], mk, color="#000000", ms=3.2, mfc="white", ls="none", label=f"{src} (measured)")
    ax.set_xlabel("elbow flexion (deg)"); ax.set_ylabel("max isometric flexion moment (N m)")
    ax.set_title("Flexors", loc="left"); ax.legend(frameon=False, fontsize=6)
    ax = axes[1]
    for (k, v), c, ls in zip(ext.items(), ("#0072B2", "#E69F00"), ("-", "--")):
        ax.plot(deg, v, color=c, ls=ls, label=f"model triceps, {k}")
    for src, mk in (("Buchanan", "s"), ("Amis", "^")):
        arr = meas[(src, "extension")]
        ax.plot(arr[:, 0], np.abs(arr[:, 1]), mk, color="#000000", ms=3.2, mfc="white", ls="none", label=f"{src} (measured)")
    ax.set_xlabel("elbow flexion (deg)"); ax.set_ylabel("max isometric extension moment (N m)")
    ax.set_title("Extensors", loc="left"); ax.legend(frameon=False, fontsize=6)
    fig.tight_layout()
    fig.savefig(RES / "phase7_torque_validation.png", dpi=300)
    fig.savefig(RES / "phase7_torque_validation.pdf")
    return out


def section33_constant_sets():
    print("### section 3.3: Thelen constant sets (brachialis plant, 60 deg) ###")
    out = {}
    for name in ("young", "old", "opensim"):
        prev = muscle.use_constant_set(name)
        try:
            a = linearize.find_equilibrium(TH)
            A_as, _ = linearize.linearize(TH, a, "activating", "shortening")
            A_al, _ = linearize.linearize(TH, a, "activating", "lengthening")
            A_ds, _ = linearize.linearize(TH, a, "deactivating", "shortening")
            errs = {f"{du:+.2f}@{ms}": prediction.error_fractions(TH, a, du, n, T_S)
                    for du in (0.01, -0.01) for ms, n in HOR.items()}
            wrong = [100 * errs[f"+0.01@{ms}"][k] for ms in HOR for k in ("wrong_activation", "wrong_velocity", "both_wrong")]
            e50 = errs["+0.01@50"]
            out[name] = dict(a_star=a, rho_a=A_as[2, 2] / A_ds[2, 2], rho_v=A_al[1, 1] / A_as[1, 1],
                             range_activating_pct=[min(wrong), max(wrong)],
                             wrong_act_50=100 * e50["wrong_activation"], wrong_vel_50=100 * e50["wrong_velocity"],
                             both_50=100 * e50["both_wrong"],
                             deact_wrong_act_50=100 * errs["-0.01@50"]["wrong_activation"])
            print(f"  {name:8s}: a* = {a:.4f}, rho_a = {A_as[2, 2] / A_ds[2, 2]:.2f}, rho_v = {A_al[1, 1] / A_as[1, 1]:.6f}; "
                  f"50 ms wrong act {100 * e50['wrong_activation']:.1f}%, wrong vel {100 * e50['wrong_velocity']:.1f}%, "
                  f"both {100 * e50['both_wrong']:.1f}%; activating range {min(wrong):.1f}-{max(wrong):.1f}%; "
                  f"deactivating wrong act 50 ms {100 * errs['-0.01@50']['wrong_activation']:.0f}%")
        finally:
            muscle.restore_constants(prev)
    return out


def main():
    RES.mkdir(exist_ok=True)
    out = dict(triceps=section2_triceps(), flexors=section31_flexors(),
               validation=section32_validation(), constant_sets=section33_constant_sets())
    json.dump(out, open(RES / "phase7_generalization.json", "w"), indent=1, default=float)
    print(f"wrote {RES / 'phase7_generalization.json'}")


if __name__ == "__main__":
    main()
