"""Poster figures for the IROS 2026 Neuromuscular Robotics workshop
(poster #5). Regenerates the key results at poster scale (large type,
fixed categorical palette, text outlined as paths) as SVG + PNG into
docs/poster/figures/. Re-computes everything from the model; nothing is
read from the results/*.png files."""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from muscle_activation_control import (antagonist, joint, linearize, mpc, muscle,  # noqa: E402
                                       prediction, prediction_pair)

JOURNAL = "--journal" in sys.argv   # journal-scale type + PDF, into docs/manuscript/figures
ICRA = "--icra" in sys.argv         # IEEE two-column sizes (exact column widths), into docs/icra2027/figures
OUT = ROOT / "docs" / ("icra2027" if ICRA else ("manuscript" if JOURNAL else "poster")) / "figures"
STYLE = ROOT / "scripts" / "icra.mplstyle"   # copy of ~/.claude/plot-theme/icra.mplstyle

# exact final sizes (inches) for the ICRA set: 3.45 = column, 7.16 = page width
ICRA_SIZES = {
    "fig1_kinks": (7.16, 2.5), "fig2_prediction": (3.45, 2.6),
    "fig3_equilibria": (7.16, 2.5), "fig4_closedloop": (3.45, 2.6),
    "fig5_pair": (7.16, 2.5), "fig6_validation": (3.45, 1.9),
}
# Okabe-Ito categorical cycle per ~/.claude/plot-color-theme.md: proposed/ours
# always blue, baselines orange/green/purple, vermilion for accents/annotations,
# black for reference/ground truth; linestyles paired per series for grayscale.
BLUE, ORANGE, GREEN, PURPLE, VERM = "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00"
C1, C2, C3, C4 = BLUE, ORANGE, GREEN, PURPLE     # fixed categorical order
ACCENT = VERM
TRUTH = "#000000"                                # reference / ground truth
INK, MUTED, GRID = "#222222", "#555555", "#B0B0B0"
THETA = np.deg2rad(60.0)
T_S = 0.010

def fs(n):
    """font size: poster value n, journal-scale (about half), or ICRA (8 pt-ish)."""
    if ICRA:
        return max(5.5, round(n * 0.36, 1))
    return n if not JOURNAL else max(7, round(n * 0.5))


plt.style.use(str(STYLE))
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": fs(22),
    "axes.titlesize": fs(24), "axes.labelsize": fs(22), "legend.fontsize": fs(19),
    "xtick.labelsize": fs(20), "ytick.labelsize": fs(20),
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6 if ICRA else (0.8 if JOURNAL else 1.4), "lines.linewidth": 1.1 if ICRA else (1.8 if JOURNAL else 3.2),
    "grid.color": GRID, "grid.linewidth": 0.4 if ICRA else (0.6 if JOURNAL else 1.0),
    "svg.fonttype": "path", "figure.dpi": 100, "savefig.dpi": 220,
})


def save(fig, name):
    if ICRA and name in ICRA_SIZES:
        fig.set_size_inches(*ICRA_SIZES[name])
    fig.tight_layout(pad=0.3 if ICRA else 1.08)
    fig.savefig(OUT / f"{name}.svg")
    fig.savefig(OUT / f"{name}.png")
    fig.savefig(OUT / f"{name}.pdf")   # vector, for the LaTeX (gemini) poster and the papers
    plt.close(fig)
    print("saved", name)


def label_end(ax, x, y, text, color, dx=8, dy=0, va="center"):
    ax.annotate(text, (x[-1], y[-1]), xytext=(dx, dy), textcoords="offset points",
                color=color, fontsize=fs(19), fontweight="bold", va=va, ha="left")


# ---------- F1: the two kinks, at the Phase 2 equilibrium ----------
def fig_kinks():
    a_star = linearize.find_equilibrium(THETA)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.0))
    ax = axes[0]
    y = np.linspace(-0.06, 0.06, 401)          # y = u - a
    adot = muscle.activation_derivative(a_star, a_star + y)
    ax.plot(y[y <= 0], adot[y <= 0], color=ORANGE, ls="--", label="deactivating: u < a")
    ax.plot(y[y >= 0], adot[y >= 0], color=BLUE, label="activating: u > a")
    ax.plot([0], [0], "o", ms=fs(13), color=ACCENT, zorder=5)
    ax.annotate("equilibrium\n$u^\\ast = a^\\ast$", (0, 0), xytext=(-0.055, 2.2), fontsize=fs(19), color=ACCENT,
                arrowprops=dict(arrowstyle="-", color=ACCENT, lw=1.5))
    rho = muscle.TAU_DEACT / (muscle.TAU_ACT * (0.5 + 1.5 * a_star) ** 2)
    ax.text(0.012, -0.75, f"slope ratio {rho:.1f}x", fontsize=fs(20), color=INK)
    ax.axhline(0, color=GRID, lw=1); ax.axvline(0, color=GRID, lw=1)
    ax.set_xlabel(r"excitation minus activation, $u-a$")
    ax.set_ylabel(r"$\dot a$ (1/s)")
    ax.set_title("Activation switch (Winters 1995)", loc="left")
    ax.legend(frameon=False, loc="upper left")

    ax = axes[1]
    v = np.linspace(-1.0, 1.0, 401)             # normalized fiber velocity
    f = muscle.contractile_force(a_star, 1.0, v) / a_star
    ax.plot(v[v <= 0], f[v <= 0], color=BLUE, label="shortening: v < 0")
    ax.plot(v[v >= 0], f[v >= 0], color=ORANGE, ls="--", label="lengthening: v > 0")
    ax.plot([0], [1], "o", ms=fs(13), color=ACCENT, zorder=5)
    ax.annotate("equilibrium\n$v^\\ast = 0$", (0, 1), xytext=(-0.95, 1.22), fontsize=fs(19), color=ACCENT,
                arrowprops=dict(arrowstyle="-", color=ACCENT, lw=1.5))
    ax.text(0.3, 0.88, "slope ratio 2x\n(exactly, for any $a^\\ast$)", fontsize=fs(20), color=INK)
    ax.axvline(0, color=GRID, lw=1)
    ax.set_xlabel(r"normalized fibre velocity $v/v_{\max}$")
    ax.set_ylabel(r"force / $(a f_L)$")
    ax.set_title("Force-velocity switch (Katz 1939)", loc="left")
    ax.legend(frameon=False, loc="lower right")
    save(fig, "fig1_kinks")


# ---------- F2: one-shot prediction, Phase 4 ----------
def fig_prediction():
    a_star = linearize.find_equilibrium(THETA)
    n = 50
    t = np.arange(n + 1) * T_S * 1000
    true = prediction.true_trajectory(THETA, a_star, 0.01, n, T_S)
    fig, ax = plt.subplots(figsize=(11, 5.4))
    tt = np.rad2deg(true[:, 0] - THETA)
    ax.plot(t, tt, color=TRUTH, lw=1.8 if ICRA else (2.8 if JOURNAL else 5), label="true nonlinear model")
    label_end(ax, t, tt, "true model", TRUTH, dy=5)
    cases = [(("activating", "shortening"), "matched branch", BLUE, "-"),
             (("deactivating", "shortening"), "wrong activation branch", ORANGE, "--"),
             (("activating", "lengthening"), "wrong velocity branch", GREEN, "-."),
             (("deactivating", "lengthening"), "both wrong", PURPLE, ":")]
    for br, name, c, ls in cases:
        p = prediction.linear_prediction(THETA, a_star, 0.01, n, T_S, *br)
        y = np.rad2deg(p[:, 0] - THETA)
        ax.plot(t, y, color=c, ls=ls, label=name)
        label_end(ax, t, y, name, c, dy=-11 if name == "matched branch" else 0)
    ax.set_xlabel("time after excitation step (ms)")
    ax.set_ylabel("joint angle change (deg)")
    ax.set_title("One-shot linear predictions vs the true model", loc="left")
    ax.grid(axis="y")
    ax.set_xlim(0, 700)
    ax.set_ylim(-0.15, 6.2)
    save(fig, "fig2_prediction")


# ---------- F3: every equilibrium, Phase 6a held-load sweep ----------
def fig_equilibria():
    a_vals = np.array([0.02, 0.05, 0.1, 0.133, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.884, 0.95])
    a_unity = (np.sqrt(muscle.TAU_DEACT / muscle.TAU_ACT) - 0.5) / 1.5
    res = {h: {k: [] for k in ("wrong_activation", "wrong_velocity", "both_wrong")} for h in (5, 50)}
    for a in a_vals:
        tau = prediction.external_torque_for_equilibrium(THETA, a)
        for h in (5, 50):
            ef = prediction.error_fractions(THETA, a, 0.01, h, T_S, tau)
            for k in res[h]:
                res[h][k].append(100 * ef[k])
    fig, axes = plt.subplots(1, 3, figsize=(16, 7.0), sharey=True)
    titles = {"wrong_activation": "wrong activation", "wrong_velocity": "wrong velocity",
              "both_wrong": "both wrong"}
    for ax, k in zip(axes, titles):
        y50, y500 = np.array(res[5][k]), np.array(res[50][k])
        ax.plot(a_vals, y50, "-o", color=C1, ms=fs(7) if (JOURNAL or ICRA) else 7, label="50 ms horizon")
        ax.plot(a_vals, y500, "--s", color=C2, ms=fs(7) if (JOURNAL or ICRA) else 7, label="500 ms horizon")
        # no end labels: they collided near a* = 0.95; the legend names the series
        ax.axvline(a_unity, color=ACCENT, ls=":", lw=2)
        ax.set_title(titles[k], loc="left")
        ax.set_xlabel(r"equilibrium activation $a^\ast$")
        ax.grid(axis="y")
        ax.set_xlim(0, 1.0)
        ax.set_xticks([0, 0.3, 0.6, 0.9])   # end ticks of adjacent panels collided ("1.20.0")
        ax.set_ylim(0, 96)
    axes[0].set_ylabel("prediction error, % of true swing")
    axes[0].text(a_unity - 0.03, 91, r"$a^\ast = 0.884$: branches agree", color=ACCENT, fontsize=fs(17), ha="right")
    axes[1].legend(frameon=False, loc="upper left")
    save(fig, "fig3_equilibria")


# ---------- F4: closed loop, Phase 6b step ----------
def fig_closedloop():
    """Both steps of the Phase 6b task (60 -> 63 -> 57 deg), four branch
    rules, identical weights; each panel zoomed near its target so the
    overshoots are visible. Relinearization (the rule the papers
    recommend) is the blue series."""
    from phase6b_closed_loop import reference as task_ref
    a0 = linearize.find_equilibrium(THETA)
    rules = {"fixed_naive": ("naive fixed branch", GREEN, "-."),
             "blended": ("blended (central diff.)", ORANGE, "--"),
             "tracking": ("sign test", PURPLE, (0, (1.2, 1.2))),
             "current_state": ("relinearize at state", BLUE, "-")}
    logs = {m: mpc.simulate(mpc.BranchAwareMPC(m), task_ref, 2.1, [THETA, 0.0, a0]) for m in rules}
    t = logs["tracking"]["t"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.9))
    windows = [(0.25, 0.85, 62.0, 63.8, "+3 deg step"), (1.25, 1.85, 56.4, 58.2, "-6 deg step")]
    for ax, (t0, t1, y0, y1, title) in zip(axes, windows):
        ax.plot(t, np.rad2deg(logs["tracking"]["theta_ref"]), color=TRUTH, ls=":", lw=plt.rcParams["lines.linewidth"] * 0.8,
                label="reference")
        for m, (name, c, ls) in rules.items():
            ax.plot(t, np.rad2deg(logs[m]["x"][:, 0]), color=c, ls=ls, label=name)
        ax.set_xlim(t0, t1); ax.set_ylim(y0, y1)
        ax.set_title(title, loc="left")
        ax.set_xlabel("time (s)")
        ax.grid(axis="y")
    axes[0].set_ylabel("elbow angle (deg)")
    axes[1].legend(frameon=False, loc="upper right", fontsize=fs(15), handlelength=2.2)
    save(fig, "fig4_closedloop")
    return {m: (np.rad2deg(logs[m]["x"][(t >= 0.3) & (t < 1.3), 0]).max() - 63.0,
                57.0 - np.rad2deg(logs[m]["x"][(t >= 1.3) & (t < 2.1), 0]).min()) for m in rules}


# ---------- F5: antagonist pair, Phase 6c ----------
def fig_pair():
    a_bal = antagonist.balanced_cocontraction(THETA)
    sweep = sorted({0.0, 0.02, 0.05, 0.1, round(a_bal, 4), 0.2, 0.3, 0.5})
    res = {h: {k: [] for k in ("wrong_bra_activation", "wrong_velocity", "all_wrong")} for h in (5, 50)}
    ratios = []
    for a_t in sweep:
        a_b = antagonist.find_equilibrium(THETA, a_t)
        ratios.append(antagonist.damping_contributions(THETA, a_b, a_t)[2])
        for h in (5, 50):
            ef = prediction_pair.error_fractions(THETA, a_b, a_t, 0.01, h, T_S)
            for k in res[h]:
                res[h][k].append(100 * ef[k])
    fig, axes = plt.subplots(1, 3, figsize=(16, 7.4))
    # same ratio with the Holzbaur-consistent triceps moment arm (x0.81),
    # docs/theory/phase7_limitations.md section 9.2: no balance point
    s_h = 0.810
    antagonist.set_tri_moment_arm_scale(s_h)
    sweep_h = [a for a in np.linspace(0.0, 0.7, 15)]
    ratios_h = [antagonist.damping_contributions(THETA, antagonist.find_equilibrium(THETA, a_t), a_t)[2]
                for a_t in sweep_h]
    antagonist.set_tri_moment_arm_scale(1.0)
    ax = axes[0]
    # green, not orange: orange means "500 ms" in the other two panels
    ms_pt = fs(7) if (JOURNAL or ICRA) else 7
    ax.plot(sweep, ratios, "-o", color=BLUE, ms=ms_pt, label="Murray triceps")
    ax.plot(sweep_h, ratios_h, "--s", color=GREEN, ms=ms_pt * 0.8, label="triceps x0.81")
    ax.legend(frameon=False, fontsize=fs(15), loc="upper right")
    ax.axhline(1.0, color=GRID, lw=1.5); ax.axhline(2.0, color=GRID, lw=1.5)
    ax.axvline(a_bal, color=ACCENT, ls=":", lw=2)
    ax.text(a_bal + 0.03, 1.22, f"kink cancels at\n$a_t^\\ast$ = {a_bal:.3f}", color=ACCENT, fontsize=fs(18))
    ax.set_title("slope ratio", loc="left")
    ax.set_ylabel("extending / flexing damping")
    ax.set_xlabel(r"triceps co-contraction $a_t^\ast$")
    ax.set_ylim(0.8, 2.1)
    ax.set_xlim(0, 0.72)
    titles = {"wrong_velocity": "wrong velocity", "all_wrong": "all three wrong"}
    for ax, k in zip(axes[1:], titles):
        y50, y500 = np.array(res[5][k]), np.array(res[50][k])
        ax.plot(sweep, y50, "-o", color=C1, ms=fs(7) if (JOURNAL or ICRA) else 7, label="50 ms")
        ax.plot(sweep, y500, "--s", color=C2, ms=fs(7) if (JOURNAL or ICRA) else 7, label="500 ms")
        # no end labels: both curves end near 1-3 % at a_t* = 0.5 and the
        # labels collided; the legend in the centre panel names them
        ax.axvline(a_bal, color=ACCENT, ls=":", lw=2)
        ax.set_title(titles[k], loc="left")
        ax.set_xlabel(r"triceps co-contraction $a_t^\ast$")
        ax.set_xlim(0, 0.6)
        ax.grid(axis="y")
    axes[1].set_ylabel("prediction error, % of true swing")
    axes[1].legend(frameon=False, loc="upper right")
    save(fig, "fig5_pair")
    return a_bal


# ---------- F6: model validation, Phase 1 ----------
# Male-specimen (anatomical) brachialis moment arm, pixel-level digitization
# (2026-09-11) of Murray et al. (1995) Fig. 4 "Male Specimen" panel, valid
# measured range 25-115 deg (below 25 deg the curve meets BRD); reproduces
# Murray's Table 2 male spread (58.9% vs 58%). data/murray1995_fig4_digitized.csv.
_MALE_DEG = np.array([25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115], dtype=float)
_MALE_CM = np.array([1.354, 1.495, 1.566, 1.671, 1.817, 1.910, 2.013, 2.144, 2.296, 2.403, 2.510, 2.639, 2.755, 2.864, 2.984, 3.067, 3.145, 3.223, 3.286])


def fig_validation():
    th = np.linspace(joint.THETA_MIN, joint.THETA_MAX, 600)
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.9))
    ax = axes[0]
    ax.plot(np.rad2deg(th), joint.moment_arm(th) * 100, color=C1, zorder=3,
            label="this model (interp.)")
    ax.plot(joint._R_THETA_DEG[::2], joint._R_THETA_M[::2] * 100, "o", color=C1,
            ms=fs(7), mfc="white", zorder=4, label="digitized pts (Murray, model)")
    ax.plot(_MALE_DEG, _MALE_CM, "s", ls=(0, (2.6, 1.8)), color=MUTED,
            ms=max(2.2, fs(6) * 0.6), lw=0.9 if ICRA else 1.4,
            zorder=2, label="Murray male specimen")
    ax.set_xlabel("elbow flexion (deg)"); ax.set_ylabel("moment arm (cm)")
    ax.set_ylim(0.4, 3.6)
    ax.set_title(r"Moment arm $r(\theta)$", loc="left")
    ax.legend(frameon=False, fontsize=fs(15), loc="lower right", handlelength=1.6,
              borderaxespad=0.2, labelspacing=0.25)
    ax = axes[1]
    ax.plot(np.rad2deg(th), joint.isometric_torque(th, 1.0), color=C1, label="brachialis, a = 1")
    ax.plot(np.rad2deg(th), joint.gravity_torque(th), color=TRUTH, ls=":", lw=2.5, label="forearm + hand gravity")
    ax.set_xlabel("elbow flexion (deg)"); ax.set_ylabel("torque (N m)")
    ax.set_title(r"Torque$-$angle, $a=1$", loc="left")
    ax.legend(frameon=False, fontsize=fs(15))
    save(fig, "fig6_validation")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    fig_kinks(); fig_prediction(); fig_equilibria(); fig_validation()
    a_bal = fig_pair()
    ov = fig_closedloop()
    print("balance", a_bal, "overshoot", ov)
