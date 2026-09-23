"""Phase 6b (docs/theory/phase6b_closed_loop.md): closed-loop MPC on the
true nonlinear single-joint Hill plant, four prediction-model variants
that differ ONLY in how the branch is chosen at the reference
equilibrium. Same weights, horizon, task, plant, initial condition."""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from muscle_activation_control import linearize, mpc  # noqa: E402

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
THETA0 = np.deg2rad(60.0)
T_END = 4.3
SEGMENTS = {"step_up": (0.3, 1.3), "step_down": (1.3, 2.3), "sinusoid": (2.3, 4.3)}
PREVIEW_S = 0.2  # = horizon * T_s


def reference(t):
    t = np.asarray(t, dtype=float)
    ref = np.full(t.shape, 60.0)
    ref = np.where(t >= 0.3, 63.0, ref)
    ref = np.where(t >= 1.3, 57.0, ref)
    sin_mask = t >= 2.3
    ref = np.where(sin_mask, 60.0 + 3.0 * np.sin(2 * np.pi * 0.5 * (t - 2.3)), ref)
    return np.deg2rad(ref)


def segment_metrics(log):
    out = {}
    t = log["t"]
    err = np.rad2deg(log["x"][:, 0] - log["theta_ref"])
    res = np.rad2deg(log["residual_theta"])
    for name, (t0, t1) in SEGMENTS.items():
        m = (t >= t0) & (t < t1)
        mr = m[:-1]
        seg = {"rms_track_deg": float(np.sqrt(np.mean(err[m] ** 2))),
               "peak_track_deg": float(np.max(np.abs(err[m]))),
               "rms_residual_deg": float(np.sqrt(np.mean(res[mr] ** 2)))}
        if name.startswith("step"):
            # settle/offset judged up to PREVIEW_S before the next reference
            # change, since the controller previews the reference and starts
            # moving early (by design, not a tracking failure)
            ms = m & (t < t1 - PREVIEW_S)
            target = np.rad2deg(log["theta_ref"][ms][-1])
            th = np.rad2deg(log["x"][ms, 0])
            within = np.abs(th - target) <= 0.2
            settle_idx = next((i for i in range(len(within)) if within[i:].all()), None)
            seg["settle_s"] = None if settle_idx is None else float(t[ms][settle_idx] - t0)
            seg["final_offset_deg"] = float(th[-1] - target)
            direction = np.sign(target - 60.0)
            seg["overshoot_deg"] = float(np.maximum(0.0, direction * (th - target)).max())
        out[name] = seg
    return out


def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    a0 = linearize.find_equilibrium(THETA0)
    x0 = np.array([THETA0, 0.0, a0])
    logs, table = {}, {}
    for mode in mpc.MODES:
        ctrl = mpc.BranchAwareMPC(mode)
        log = mpc.simulate(ctrl, reference, T_END, x0)
        logs[mode] = log
        table[mode] = {"overall": mpc.metrics(log), **segment_metrics(log)}
        if mode in ("tracking", "current_state", "trajectory"):
            br = log["branch"]
            table[mode]["branch_counts"] = {
                f"{a}/{v}": sum(1 for b in br if b == (a, v))
                for a in ("activating", "deactivating") for v in ("shortening", "lengthening")}
        if mode == "trajectory":
            table[mode]["traj_iterations_mean"] = float(np.mean(ctrl.traj_iters))
            table[mode]["traj_iterations_max"] = int(np.max(ctrl.traj_iters))
            table[mode]["traj_fallbacks"] = ctrl.traj_fallbacks
        print(f"--- {mode} ---")
        print(json.dumps(table[mode], indent=1))

    print("\n=== summary (deg) ===")
    print(f"{'variant':>14s} {'rms track':>10s} {'peak track':>10s} {'rms resid':>10s} "
          f"{'peak resid':>10s} {'up settle':>10s} {'up ovs':>8s} {'up off':>8s} {'dn settle':>10s} {'dn ovs':>8s} {'dn off':>8s} {'sin rms':>8s}")
    for mode, m in table.items():
        o, up, dn, sn = m["overall"], m["step_up"], m["step_down"], m["sinusoid"]
        fmt = lambda v: "  n/a" if v is None else f"{v:.3f}"
        print(f"{mode:>14s} {o['rms_track_deg']:10.3f} {o['peak_track_deg']:10.3f} "
              f"{o['rms_residual_deg']:10.4f} {o['peak_residual_deg']:10.4f} "
              f"{fmt(up['settle_s']):>10s} {up['overshoot_deg']:8.3f} {up['final_offset_deg']:+8.3f} "
              f"{fmt(dn['settle_s']):>10s} {dn['overshoot_deg']:8.3f} {dn['final_offset_deg']:+8.3f} {sn['rms_track_deg']:8.3f}")

    with open(RESULTS_DIR / "phase6b_metrics.json", "w") as f:
        json.dump(table, f, indent=1)

    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    t = logs["tracking"]["t"]
    axes[0].plot(t, np.rad2deg(logs["tracking"]["theta_ref"]), "k:", lw=1.5, label="reference")
    for mode, log in logs.items():
        axes[0].plot(t, np.rad2deg(log["x"][:, 0]), lw=1.2, label=mode)
        axes[1].plot(t[:-1], log["u"], lw=1.0, label=mode)
        axes[2].plot(t[:-1], np.rad2deg(log["residual_theta"]), lw=1.0, label=mode)
    axes[0].set_ylabel("theta (deg)")
    axes[1].set_ylabel("excitation u")
    axes[2].set_ylabel("one-step residual |dtheta| (deg)")
    axes[2].set_yscale("log")
    axes[2].set_xlabel("time (s)")
    for ax in axes:
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8, ncol=6)
    fig.suptitle("Phase 6b: closed-loop MPC, four branch-selection rules, identical weights/horizon/plant")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "phase6b_closed_loop.png", dpi=150)

    # zoom on the step-up transient
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    m = (t >= 0.25) & (t <= 1.3)
    axes[0].plot(t[m], np.rad2deg(logs["tracking"]["theta_ref"][m]), "k:", lw=1.5, label="reference")
    for mode, log in logs.items():
        axes[0].plot(t[m], np.rad2deg(log["x"][m, 0]), lw=1.3, label=mode)
        axes[1].plot(t[:-1][m[:-1]], np.rad2deg(log["residual_theta"][m[:-1]]), lw=1.0, label=mode)
    axes[0].set_ylabel("theta (deg)")
    axes[1].set_ylabel("one-step residual (deg)")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("time (s)")
    axes[0].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.suptitle("Phase 6b: +3 deg step, zoom")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "phase6b_step_zoom.png", dpi=150)
    print(f"saved plots to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
