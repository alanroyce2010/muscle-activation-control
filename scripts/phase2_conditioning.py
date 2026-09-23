"""Phase 2 conditioning check: does this model's ZOH-discretized A matrix
blow up the way origami-arm-control's pneumatic pressure state did
(max|A_d| ~ 900,000 at a practical sample time, mpc_control.md section 2)?

See docs/theory/phase2_linearization.md for the full derivation.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from muscle_activation_control import linearize  # noqa: E402

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"

THETA_STAR = np.deg2rad(60.0)
T_S_VALUES_MS = [1, 2, 5, 10, 20, 50, 100]
T_S_VALUES = [t / 1000.0 for t in T_S_VALUES_MS]


def main():
    RESULTS_DIR.mkdir(exist_ok=True)

    a_star, results = linearize.conditioning_sweep(THETA_STAR, T_S_VALUES)

    print("=== Phase 2: equilibrium ===")
    print(f"theta* = 60 deg, thetadot* = 0, a* = u* = {a_star:.4f}")
    print()

    print("=== continuous-time A matrix eigenvalues, all 4 branch combinations ===")
    seen = set()
    for r in results:
        key = (r["activation_branch"], r["velocity_branch"])
        if key in seen:
            continue
        seen.add(key)
        eigs = np.linalg.eigvals(r["a_cont"])
        print(f"{r['activation_branch']:>12s} / {r['velocity_branch']:>11s}: "
              f"eigenvalues = {np.round(eigs, 2)}  "
              f"(1/s; time constants: {np.round(1.0/np.abs(eigs[np.abs(eigs)>1e-9]), 4)} s)")
    print()

    print("=== max|A_d| vs. sample time, all 4 branch combinations ===")
    print(f"{'T_s (ms)':>10s}  {'act/short':>12s}  {'act/long':>12s}  "
          f"{'deact/short':>12s}  {'deact/long':>12s}")
    by_key = {}
    for r in results:
        by_key.setdefault(r["t_s"], {})[
            (r["activation_branch"], r["velocity_branch"])
        ] = r["max_abs_a_d"]

    for t_s in T_S_VALUES:
        row = by_key[t_s]
        print(
            f"{t_s*1000:>10.1f}  "
            f"{row[('activating','shortening')]:>12.3g}  "
            f"{row[('activating','lengthening')]:>12.3g}  "
            f"{row[('deactivating','shortening')]:>12.3g}  "
            f"{row[('deactivating','lengthening')]:>12.3g}"
        )

    # For reference: origami-arm-control's pneumatic pressure case hit
    # max|A_d| ~ 900,000 at its practical sample time (mpc_control.md sec 2).
    pneumatic_reference = 900_000
    print(f"\nreference: origami-arm-control pneumatic pressure case, "
          f"max|A_d| ~ {pneumatic_reference:,} at its practical T_s")

    worst_max = max(r["max_abs_a_d"] for r in results)
    print(f"worst max|A_d| found here, across all branches and T_s tested: "
          f"{worst_max:.3g}")
    if worst_max > pneumatic_reference / 10:
        print("=> comparable order of magnitude to the pneumatic blowup: "
              "SUGGESTS TRANSFER")
    elif worst_max > 100:
        print("=> elevated but well below the pneumatic blowup: "
              "SUGGESTS PARTIAL TRANSFER")
    else:
        print("=> stays well-conditioned at all sample times tested: "
              "SUGGESTS NON-TRANSFER")

    fig, ax = plt.subplots(figsize=(7, 5))
    for (activation_branch, velocity_branch) in linearize.BRANCH_COMBINATIONS:
        vals = [by_key[t_s][(activation_branch, velocity_branch)] for t_s in T_S_VALUES]
        ax.plot(
            [t * 1000 for t in T_S_VALUES],
            vals,
            marker="o",
            label=f"{activation_branch}/{velocity_branch}",
        )
    ax.axhline(pneumatic_reference, color="red", linestyle="--",
               label="origami-arm-control pneumatic reference (~900,000)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("sample time T_s (ms)")
    ax.set_ylabel("max|A_d| entry")
    ax.set_title("ZOH-discretized A matrix conditioning vs. sample time")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out_path = RESULTS_DIR / "phase2_conditioning.png"
    fig.savefig(out_path, dpi=150)
    print(f"\nsaved plot: {out_path}")


if __name__ == "__main__":
    main()
