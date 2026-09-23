"""Phase 6c (docs/theory/phase6c_antagonist_pair.md): antagonist pair.
1. triceps validation (moment arm average vs Holzbaur, isometric torque),
2. switch structure at co-contraction equilibria (continuity, disjointness,
   closed-form force-velocity ratio vs finite differences),
3. wrong-branch one-step prediction error vs co-contraction level."""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from muscle_activation_control import antagonist, joint, muscle, prediction_pair  # noqa: E402

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
THETA_STAR = np.deg2rad(60.0)
T_S = 0.010
HORIZONS = {50: 5, 100: 10, 200: 20, 500: 50}
A_T_SWEEP = [0.0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5]
KEYS = ["matched", "wrong_bra_activation", "wrong_tri_activation", "wrong_velocity", "all_wrong"]


def rho_a(a):
    return muscle.TAU_DEACT / (muscle.TAU_ACT * (0.5 + 1.5 * a) ** 2)


def validation():
    print("### 1. triceps validation ###")
    th = np.linspace(joint.THETA_MIN, joint.THETA_MAX, 1301)
    r = antagonist.tri_moment_arm(th)
    r_avg = np.trapezoid(r, th) / (th[-1] - th[0])
    print(f"digitized r_tri: peak {r.max()*100:.2f} cm at {np.rad2deg(th[np.argmax(r)]):.0f} deg, "
          f"ROM average {r_avg*100:.2f} cm (Holzbaur 2005 ma_avg 2.1 cm: {100*(r_avg/0.021-1):+.0f}%)")
    l_bar, cos_a = antagonist.tri_fiber_state(th)
    print(f"l_bar_tri over ROM: {l_bar.min():.3f} .. {l_bar.max():.3f}; cos(alpha): {cos_a.min():.4f} .. {cos_a.max():.4f}")
    tau_t = antagonist.tri_torque(th, 0.0, 1.0)
    tau_b = joint.isometric_torque(th, 1.0)
    print(f"isometric extension torque at a=1: peak {tau_t.max():.1f} N*m at {np.rad2deg(th[np.argmax(tau_t)]):.0f} deg; "
          f"brachialis flexion peak {tau_b.max():.1f} N*m at {np.rad2deg(th[np.argmax(tau_b)]):.0f} deg")
    tau_t0 = antagonist.tri_torque(th, 0.0, 0.0)
    print(f"passive triceps torque at a=0: max {tau_t0.max():.2f} N*m at {np.rad2deg(th[np.argmax(tau_t0)]):.0f} deg\n")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(np.rad2deg(th), r * 100, label="triceps (extension, magnitude)")
    ax[0].plot(np.rad2deg(th), joint.moment_arm(th) * 100, label="brachialis (flexion)")
    ax[0].set_xlabel("elbow flexion (deg)"); ax[0].set_ylabel("moment arm (cm)"); ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[1].plot(np.rad2deg(th), tau_b, label="brachialis, a=1")
    ax[1].plot(np.rad2deg(th), -tau_t, label="triceps (lat+med), a=1")
    ax[1].plot(np.rad2deg(th), -tau_t0, "--", label="triceps passive, a=0")
    ax[1].plot(np.rad2deg(th), -joint.gravity_torque(th), ":", label="gravity (forearm+hand)")
    ax[1].set_xlabel("elbow flexion (deg)"); ax[1].set_ylabel("isometric torque (N*m, +flexion)"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(RESULTS_DIR / "phase6c_triceps_validation.png", dpi=150)


def structure():
    print("### 2. switch structure at theta*=60 deg ###")
    print(f"{'a_t*':>6s} {'a_b*':>7s} {'rho_a(b)':>8s} {'rho_a(t)':>8s} {'D_b':>8s} {'D_t':>8s} {'FV ratio pred':>13s} {'FV ratio FD':>11s} {'max resid':>10s} disjoint")
    for a_t in A_T_SWEEP:
        a_b = antagonist.find_equilibrium(THETA_STAR, a_t)
        v = antagonist.verify_switches(THETA_STAR, a_b, a_t)
        d_b, d_t, ratio = antagonist.damping_contributions(THETA_STAR, a_b, a_t)
        x = np.array([THETA_STAR, 0, a_b, a_t]); u = np.array([a_b, a_t])
        a_flex, _, _ = antagonist.linearize(x, u, "activating", "activating", "flexing")
        a_ext, _, _ = antagonist.linearize(x, u, "activating", "activating", "extending")
        ratio_fd = a_ext[1, 1] / a_flex[1, 1]
        res = max(v[k]["residual"] for k in ("bra_activation", "tri_activation", "velocity"))
        print(f"{a_t:6.2f} {a_b:7.4f} {rho_a(a_b):8.2f} {rho_a(a_t):8.2f} {d_b:8.4f} {d_t:8.4f} {ratio:13.4f} {ratio_fd:11.4f} {res:10.2e} {v['disjoint']}")
    a_bal = antagonist.balanced_cocontraction(THETA_STAR)
    if a_bal is not None:
        a_b = antagonist.find_equilibrium(THETA_STAR, a_bal)
        print(f"force-velocity kink cancels (D_b = D_t) at a_t* = {a_bal:.4f} (a_b* = {a_b:.4f})")
    else:
        print("no balanced co-contraction level in (0,1)")
    v = antagonist.verify_switches(THETA_STAR, antagonist.find_equilibrium(THETA_STAR, 0.1), 0.1)
    for k in ("bra_activation", "tri_activation", "velocity"):
        print(f"  e_{k} (a_t*=0.1) = {np.array2string(v[k]['e'], precision=3)}")
    print()
    return a_bal


def quantify(a_bal):
    print("### 3. wrong-branch prediction error vs co-contraction (theta*=60 deg, delta_u_b=+/-0.01) ###")
    sweep = sorted(set(A_T_SWEEP + ([round(a_bal, 4)] if a_bal else [])))
    res = {du: {h: {k: [] for k in KEYS} for h in HORIZONS} for du in (0.01, -0.01)}
    for du in (0.01, -0.01):
        print(f"delta_u_b = {du:+.2f}")
        print(f"{'a_t*':>7s} {'a_b*':>7s} | " + "  ".join(f"{h}ms m/wb/wt/wv/all" for h in HORIZONS) + " | occ  a_t drift")
        for a_t in sweep:
            a_b = antagonist.find_equilibrium(THETA_STAR, a_t)
            row, occ, drift = [], (1.0, 1.0), 0.0
            for h, n in HORIZONS.items():
                ef = prediction_pair.error_fractions(THETA_STAR, a_b, a_t, du, n, T_S)
                for k in KEYS:
                    res[du][h][k].append(ef[k])
                row.append(f"{ef['matched']:.3f}/{ef['wrong_bra_activation']:.2f}/{ef['wrong_tri_activation']:.1e}/"
                           f"{ef['wrong_velocity']:.2f}/{ef['all_wrong']:.2f}")
                occ = (min(occ[0], ef["on_matched_bra_act"]), min(occ[1], ef["on_matched_vel"]))
                drift = max(drift, ef["a_t_drift"])
            flag = "" if min(occ) == 1.0 else "  <-- crosses a switch"
            print(f"{a_t:7.4f} {a_b:7.4f} | " + "  ".join(row) + f" | ({occ[0]:.2f},{occ[1]:.2f}) {drift:.1e}{flag}")
        print()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
    for ax, key in zip(axes, ["wrong_bra_activation", "wrong_velocity", "all_wrong"]):
        for h in HORIZONS:
            ax.plot(sweep, 100 * np.array(res[0.01][h][key]), "-o", ms=3, label=f"{h} ms")
        ax.plot(sweep, 100 * np.array(res[0.01][500]["matched"]), "k--", lw=1, label="matched (500 ms)")
        if a_bal:
            ax.axvline(a_bal, color="gray", ls=":", lw=1)
        ax.set_title(key.replace("_", " ")); ax.set_xlabel("triceps co-contraction a_t*"); ax.grid(alpha=0.3)
    axes[0].set_ylabel("final theta error, % of true swing"); axes[0].legend(fontsize=7)
    fig.suptitle("Phase 6c: antagonist pair, brachialis excitation step +0.01 at theta*=60 deg (dotted: D_b = D_t)")
    fig.tight_layout(); fig.savefig(RESULTS_DIR / "phase6c_cocontraction_sweep.png", dpi=150)
    print(f"saved plots to {RESULTS_DIR}")


def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    validation()
    a_bal = structure()
    quantify(a_bal)


if __name__ == "__main__":
    main()
