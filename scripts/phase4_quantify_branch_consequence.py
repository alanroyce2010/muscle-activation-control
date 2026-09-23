"""Phase 4 (CLAUDE.md section 4.2 step 2): does the branch ambiguity from
Phase 2/3 matter for a realistic control design, or is it a curiosity?

Scenario: a successive-linearization MPC's internal prediction model
typically linearizes once per control step and holds that linear model
fixed over the prediction horizon (the same simplification
origami-arm-control's MPCController uses, CLAUDE.md section 2). Simulate
what happens if that one-shot linearization, taken at this project's
equilibrium, picks the "wrong" branch: propagate a step excitation change
through (a) the true nonlinear model and (b) each of the 4 fixed
discrete-time linearizations from Phase 2, and compare, in physical units
(degrees of joint angle, N*m of torque), over a plausible MPC horizon.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from muscle_activation_control import joint, linearize  # noqa: E402

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"

THETA_STAR = np.deg2rad(60.0)
T_S = 0.010  # 10 ms, a representative control rate from the Phase 2 sweep
HORIZON_STEPS = 20  # 200 ms, a plausible MPC prediction horizon
DELTA_U_MAG = 0.01  # step change in excitation from equilibrium -- kept
# small deliberately (this is the small-signal regime linearization is
# actually supposed to be valid in; DELTA_U_MAG=0.05 was tried first and
# produced excursions large enough (>25 deg over 500ms) that generic
# linearization breakdown swamped the branch-choice effect being
# isolated here -- see docs/theory/phase4_quantification.md section 2 for
# that finding and why this smaller value was kept as the primary result
# instead.


def true_trajectory(a_star, delta_u, n_steps, t_s):
    u = a_star + delta_u
    x0 = np.array([THETA_STAR, 0.0, a_star])
    t_eval = np.arange(n_steps + 1) * t_s

    def rhs(t, x):
        return joint.state_derivative(x, u)

    sol = solve_ivp(rhs, (0, n_steps * t_s), x0, t_eval=t_eval, rtol=1e-10, atol=1e-12)
    return sol.y.T  # shape (n_steps+1, 3)


def linear_prediction(theta_star, a_star, delta_u, n_steps, t_s, activation_branch, velocity_branch):
    a_cont, b_cont = linearize.linearize(theta_star, a_star, activation_branch, velocity_branch)
    a_d, b_d = linearize.discretize_zoh(a_cont, b_cont, t_s)
    b_d = b_d.flatten()
    dx = np.zeros((n_steps + 1, 3))
    delta_x = np.zeros(3)
    for k in range(n_steps):
        dx[k] = delta_x
        delta_x = a_d @ delta_x + b_d * delta_u
    dx[n_steps] = delta_x
    x0 = np.array([theta_star, 0.0, a_star])
    return x0 + dx


def run_scenario(a_star, delta_u, label):
    print(f"=== {label}: delta_u = {delta_u:+.3f} (u* = {a_star:.3f}) ===")
    true_traj = true_trajectory(a_star, delta_u, HORIZON_STEPS, T_S)
    t_ms = np.arange(HORIZON_STEPS + 1) * T_S * 1000

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    ax.plot(t_ms, np.rad2deg(true_traj[:, 0] - THETA_STAR), "k-", linewidth=2.5,
             label="true nonlinear")

    errors_deg = {}
    final_torque_true = joint.isometric_torque(true_traj[-1, 0], true_traj[-1, 2])
    print(f"  true final state: theta={np.rad2deg(true_traj[-1,0]):.3f} deg, "
          f"thetadot={true_traj[-1,1]:.4f} rad/s, a={true_traj[-1,2]:.4f}")

    for act_branch, vel_branch in linearize.BRANCH_COMBINATIONS:
        pred_traj = linear_prediction(THETA_STAR, a_star, delta_u, HORIZON_STEPS, T_S,
                                        act_branch, vel_branch)
        ax.plot(t_ms, np.rad2deg(pred_traj[:, 0] - THETA_STAR), "--",
                 label=f"{act_branch}/{vel_branch}")
        err_deg = np.rad2deg(np.abs(pred_traj[:, 0] - true_traj[:, 0]))
        errors_deg[(act_branch, vel_branch)] = err_deg

        pred_torque = joint.isometric_torque(pred_traj[-1, 0], np.clip(pred_traj[-1, 2], 0, 1))
        torque_err = abs(pred_torque - final_torque_true)
        print(f"  {act_branch:>12s}/{vel_branch:>11s}: final theta error = "
              f"{err_deg[-1]:.4f} deg, max error = {err_deg.max():.4f} deg, "
              f"final torque error ~= {torque_err:.3f} N*m")

    ax.set_xlabel("time (ms)")
    ax.set_ylabel("theta - theta* (deg)")
    ax.set_title(f"Predicted vs true trajectory ({label})")
    ax.legend(fontsize=7)

    ax = axes[1]
    for key, err in errors_deg.items():
        ax.plot(t_ms, err, "--", label=f"{key[0]}/{key[1]}")
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("|theta_predicted - theta_true| (deg)")
    ax.set_title("Prediction error vs. true nonlinear trajectory")
    ax.legend(fontsize=7)

    fig.tight_layout()
    out_path = RESULTS_DIR / f"phase4_{label}.png"
    fig.savefig(out_path, dpi=150)
    print(f"  saved: {out_path}\n")
    return errors_deg


def horizon_sweep(a_star, delta_u, horizons_steps):
    print(f"=== horizon sweep, delta_u = {delta_u:+.3f} ===")
    print(f"{'horizon (ms)':>14s}  {'true swing (deg)':>18s}  {'matched err':>12s}  "
          f"{'wrong-act err':>14s}  {'wrong-vel err':>14s}  {'both-wrong err':>15s}")
    for n_steps in horizons_steps:
        true_traj = true_trajectory(a_star, delta_u, n_steps, T_S)
        true_swing = np.rad2deg(true_traj[-1, 0] - THETA_STAR)
        matched_branch = ("activating", "shortening") if delta_u > 0 else ("deactivating", "lengthening")
        wrong_act = ("deactivating" if matched_branch[0] == "activating" else "activating", matched_branch[1])
        wrong_vel = (matched_branch[0], "lengthening" if matched_branch[1] == "shortening" else "shortening")
        both_wrong = (wrong_act[0], wrong_vel[1])

        def final_err(branch):
            pred = linear_prediction(THETA_STAR, a_star, delta_u, n_steps, T_S, *branch)
            return np.rad2deg(abs(pred[-1, 0] - true_traj[-1, 0]))

        print(f"{n_steps*T_S*1000:>14.0f}  {true_swing:>18.4f}  "
              f"{final_err(matched_branch):>12.5f}  {final_err(wrong_act):>14.5f}  "
              f"{final_err(wrong_vel):>14.5f}  {final_err(both_wrong):>15.5f}")
    print()


def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    a_star = linearize.find_equilibrium(THETA_STAR)

    print("### small-signal horizon sweep (the primary Phase 4 result) ###\n")
    horizon_sweep(a_star, DELTA_U_MAG, [5, 10, 20, 50])
    horizon_sweep(a_star, -DELTA_U_MAG, [5, 10, 20, 50])

    print("### full trajectory plots, small-signal magnitude ###\n")
    run_scenario(a_star, DELTA_U_MAG, "activating_step_small")
    run_scenario(a_star, -DELTA_U_MAG, "deactivating_step_small")

    print("### full trajectory plots, larger-magnitude excitation (for contrast --")
    print("### shows generic linearization breakdown swamping the branch effect) ###\n")
    run_scenario(a_star, 0.05, "activating_step_large")
    run_scenario(a_star, -0.05, "deactivating_step_large")


if __name__ == "__main__":
    main()
