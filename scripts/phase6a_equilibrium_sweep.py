"""Phase 6a (docs/theory/phase6a_equilibrium_sweep.md): repeat the Phase 4
wrong-branch prediction-error measurement across (A) every gravity-only
equilibrium over the ROM and (B) a held-load sweep of a* at theta*=60 deg,
to test whether the "it matters" finding holds at every equilibrium and
how it scales with a* (closed-form prediction: activation-branch ratio
rho_a(a*) = tau_deact/(tau_act (0.5+1.5a*)^2), unity at a*=0.884).
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from muscle_activation_control import linearize, muscle, prediction  # noqa: E402

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
T_S = 0.010
HORIZONS = {50: 5, 100: 10, 200: 20, 500: 50}
DELTA_U = 0.01
BRANCH_KEYS = ["matched", "wrong_activation", "wrong_velocity", "both_wrong"]


def rho_a(a_star):
    return muscle.TAU_DEACT / (muscle.TAU_ACT * (0.5 + 1.5 * a_star) ** 2)


def a_star_unity():
    return (np.sqrt(muscle.TAU_DEACT / muscle.TAU_ACT) - 0.5) / 1.5


def sweep(points, label, delta_u):
    """points: list of (theta_star, a_star, tau_ext). Returns dict of arrays
    indexed [horizon][branch] plus occupancy flags."""
    res = {h: {k: [] for k in BRANCH_KEYS} for h in HORIZONS}
    occ = []
    print(f"=== {label}, delta_u={delta_u:+.3f} ===")
    print(f"{'theta*':>7s} {'a*':>6s} {'rho_a':>6s} {'tau_ext':>8s} | "
          + " ".join(f"{h}ms: m/wa/wv/bw" for h in HORIZONS) + " | on matched (act,vel)")
    for theta_star, a_star, tau_ext in points:
        row = []
        occ_worst = (1.0, 1.0)
        for h_ms, n in HORIZONS.items():
            ef = prediction.error_fractions(theta_star, a_star, delta_u, n, T_S, tau_ext)
            for k in BRANCH_KEYS:
                res[h_ms][k].append(ef[k])
            row.append(f"{ef['matched']:.3f}/{ef['wrong_activation']:.2f}/"
                       f"{ef['wrong_velocity']:.2f}/{ef['both_wrong']:.2f}")
            oa, os_ = ef["on_matched_activation_frac"], ef["on_matched_velocity_frac"]
            occ_worst = (min(occ_worst[0], oa), min(occ_worst[1], os_))
        occ.append(occ_worst)
        flag = "" if min(occ_worst) == 1.0 else "  <-- crosses a switch within horizon"
        print(f"{np.rad2deg(theta_star):7.1f} {a_star:6.3f} {rho_a(a_star):6.2f} {tau_ext:8.3f} | "
              + " ".join(row) + f" | ({occ_worst[0]:.2f},{occ_worst[1]:.2f}){flag}")
    print()
    for h in res:
        for k in res[h]:
            res[h][k] = np.array(res[h][k])
    return res, occ


def plot(xvals, xlabel, res, title, out_name, vline=None):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
    for ax, key in zip(axes, ["wrong_activation", "wrong_velocity", "both_wrong"]):
        for h in HORIZONS:
            ax.plot(xvals, 100 * res[h][key], "-o", ms=3, label=f"{h} ms")
        ax.plot(xvals, 100 * res[max(HORIZONS)]["matched"], "k--", lw=1,
                label="matched (500 ms)")
        if vline is not None:
            ax.axvline(vline, color="gray", ls=":", lw=1)
        ax.set_title(key.replace("_", " "))
        ax.set_xlabel(xlabel)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("final theta error, % of true swing")
    axes[0].legend(fontsize=7)
    fig.suptitle(title)
    fig.tight_layout()
    out = RESULTS_DIR / out_name
    fig.savefig(out, dpi=150)
    print(f"saved: {out}\n")


def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    print(f"rho_a = 1 at a* = {a_star_unity():.4f}\n")

    # Sweep A: gravity-only equilibria over the ROM
    thetas = np.deg2rad(np.arange(12.5, 123, 5))   # off the 65-deg passive kink (l_bar = 1)
    pts_a = []
    for th in thetas:
        a_s = linearize.find_equilibrium(th)
        pts_a.append((th, a_s, 0.0))
    print(f"gravity-only a* range over ROM: {min(p[1] for p in pts_a):.4f} .. "
          f"{max(p[1] for p in pts_a):.4f}\n")
    for du in (DELTA_U, -DELTA_U):
        res, _ = sweep(pts_a, "Sweep A: gravity-only equilibria", du)
        plot(np.rad2deg(thetas), "theta* (deg)", res,
             f"Sweep A: gravity-only equilibria, delta_u={du:+.2f}, T_s=10 ms",
             f"phase6a_gravity_sweep_{'pos' if du > 0 else 'neg'}.png")

    # Sweep B: held constant-torque load at theta*=60 deg, a* swept
    theta65 = np.deg2rad(60.0)
    a_vals = np.array([0.02, 0.05, 0.1, 0.133, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
                       0.884, 0.95])
    pts_b = [(theta65, a_s, prediction.external_torque_for_equilibrium(theta65, a_s))
             for a_s in a_vals]
    for du in (DELTA_U, -DELTA_U):
        res, _ = sweep(pts_b, "Sweep B: held load, theta*=60 deg", du)
        plot(a_vals, "a*", res,
             f"Sweep B: held-load equilibria at theta*=60 deg, delta_u={du:+.2f}",
             f"phase6a_load_sweep_{'pos' if du > 0 else 'neg'}.png",
             vline=a_star_unity())


if __name__ == "__main__":
    main()
