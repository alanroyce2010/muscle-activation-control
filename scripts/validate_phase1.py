"""Phase 1 validation run: CLAUDE.md section 4 step 1 / section 6.3.

Computes the isometric torque-angle curve tau(theta) = F^MT(theta;V=0,a=1)
* r(theta), checks it against the combined target in CLAUDE.md 6.3, and
saves a plot + printed report to results/.

Run with pytest (tests/test_muscle.py, tests/test_joint.py) for the
pass/fail numeric checks -- this script is for the human-readable
plot + summary, not a second copy of those checks.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from muscle_activation_control import joint, muscle  # noqa: E402

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"


def main():
    RESULTS_DIR.mkdir(exist_ok=True)

    theta_deg = np.linspace(0, 130, 300)
    theta_rad = np.deg2rad(theta_deg)

    r = joint.moment_arm(theta_rad)
    l_bar = joint.fiber_length_norm(theta_rad)
    f_l = muscle.active_force_length(l_bar)
    f_pe = muscle.passive_force_length(l_bar)
    tau_isometric = joint.isometric_torque(theta_rad, a=1.0)

    avg_r = np.trapezoid(r, theta_rad) / (theta_rad[-1] - theta_rad[0])
    peak_r = r.max()
    peak_theta_deg = theta_deg[np.argmax(r)]

    print("=== Phase 1 validation: isometric torque-angle, brachialis ===")
    print(f"moment arm: peak {peak_r*100:.2f} cm at {peak_theta_deg:.0f} deg "
          f"(target: 2.0-3.5 cm, >100 deg)")
    print(f"moment arm: average {avg_r*100:.2f} cm over 0-130 deg "
          f"(target: ~1.8 cm, Holzbaur 2005 ma_avg)")
    print(f"fiber length l_bar: min {l_bar.min():.3f}, max {l_bar.max():.3f} "
          f"(=1.0 at theta=65 deg by construction)")
    print(f"active force-length f_L: min {f_l.min():.3f}, max {f_l.max():.3f} "
          f"(peaks at 1.0 near theta=65 deg)")
    print(f"isometric torque: peak {tau_isometric.max():.2f} N*m "
          f"at {theta_deg[np.argmax(tau_isometric)]:.0f} deg")

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))

    ax = axes[0, 0]
    ax.plot(theta_deg, r * 100, label="this model (digitized+interp)")
    ax.axhline(1.8, color="gray", linestyle="--", label="Holzbaur ma_avg=1.8cm")
    ax.set_xlabel("elbow flexion angle (deg)")
    ax.set_ylabel("moment arm (cm)")
    ax.set_title("Brachialis moment arm r(theta)")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.plot(theta_deg, l_bar)
    ax.axvline(65, color="gray", linestyle="--", label="reference angle (l_bar=1)")
    ax.set_xlabel("elbow flexion angle (deg)")
    ax.set_ylabel("normalized fiber length l_bar")
    ax.set_title("Fiber length vs joint angle")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.plot(theta_deg, f_l, label="active f_L")
    ax.plot(theta_deg, f_pe, label="passive F_PE")
    ax.set_xlabel("elbow flexion angle (deg)")
    ax.set_ylabel("normalized force")
    ax.set_title("Force-length components (at l_bar(theta))")
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    ax.plot(theta_deg, tau_isometric)
    ax.set_xlabel("elbow flexion angle (deg)")
    ax.set_ylabel("isometric torque (N*m)")
    ax.set_title("Isometric elbow flexion torque, a=1")

    fig.tight_layout()
    out_path = RESULTS_DIR / "phase1_torque_angle_validation.png"
    fig.savefig(out_path, dpi=150)
    print(f"\nsaved plot: {out_path}")

    peak_ok = 0.020 <= peak_r <= 0.035 and peak_theta_deg > 100
    avg_ok = 0.016 <= avg_r <= 0.022
    print(f"\nmoment arm peak check: {'PASS' if peak_ok else 'FAIL'}")
    print(f"moment arm average check: {'PASS' if avg_ok else 'FAIL'}")


if __name__ == "__main__":
    main()
