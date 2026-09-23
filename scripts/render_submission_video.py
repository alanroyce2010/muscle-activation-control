"""Render one ICRA simulation video per scene from the repository's actual model.

Usage: .venv/bin/python scripts/render_submission_video.py
       .venv/bin/python scripts/render_submission_video.py --preview
Scene MP4s are saved in results/submission_video/scenes and replaced on reruns.
Cached trajectories are reused when verified unless --recompute is supplied.
FFmpeg is required for video export.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", "/tmp/muscle-video-matplotlib")
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "scripts")]
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.animation import FFMpegWriter
from matplotlib.patches import Circle
import numpy as np
from muscle_activation_control import antagonist, linearize, mpc, muscle, prediction
from phase6b_closed_loop import reference, segment_metrics

OUT = ROOT / "results" / "submission_video"
BG, PANEL, INK, MUTED, GRID = "#fdf0d5", "#fff8e9", "#003049", "#405967", "#d2d6cc"
NAVY, RED, STEEL, MAROON = "#003049", "#c1121f", "#669bbc", "#780000"
MOTION_GAIN = 12.0
MODES = ["fixed_naive", "blended", "tracking", "current_state"]
LABELS = ["Fixed branch", "Blended", "Sign-test branch", "Relinearize at state"]
COLORS = [MAROON, RED, STEEL, NAVY]
STYLES = [":", "--", "-.", "-"]
THETA = np.deg2rad(60.)
SCENES = [("intro", 12), ("switches", 22), ("prediction", 22),
          ("raise", 24), ("lower", 20), ("track", 24), ("pair", 18)]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13,
    "text.color": INK, "axes.labelcolor": MUTED, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.edgecolor": GRID, "axes.facecolor": PANEL,
    "figure.facecolor": BG, "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": GRID, "grid.alpha": .5, "axes.titleweight": "bold",
    "savefig.facecolor": BG})


def source_hash():
    files = sorted((ROOT / "src/muscle_activation_control").glob("*.py"))
    files += [ROOT / "scripts/phase6b_closed_loop.py", Path(__file__)]
    return hashlib.sha256(b"".join(p.read_bytes() for p in files)).hexdigest()


def data(recompute=False):
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / "trajectories.npz"
    manifest = OUT / "metrics.json"
    digest = source_hash()
    if cache.exists() and manifest.exists() and not recompute:
        metrics = json.loads(manifest.read_text())
        if metrics.get("source_sha256") == digest:
            print("Using verified trajectory cache", flush=True)
            with np.load(cache) as stored:
                return dict(stored), metrics
    a0 = linearize.find_equilibrium(THETA)
    d = {"a0": np.array(a0)}
    metrics = {"source_sha256": digest, "controllers": {},
        "configuration": {"sample_s": .01, "horizon_steps": 20,
            "q_theta": 1000, "q_velocity": 1, "r": 200,
            "plant": "rigid-tendon brachialis, nominal parameters, no noise or delay"}}
    for mode in MODES:
        print(f"Simulating {mode} ...", flush=True)
        log = mpc.simulate(mpc.BranchAwareMPC(mode), reference, 4.3, [THETA, 0., a0])
        for key in ["t", "x", "u", "theta_ref"]:
            d[f"{mode}_{key}"] = log[key]
        metrics["controllers"][mode] = segment_metrics(log)
    true = prediction.true_trajectory(THETA, a0, .01, 50, .01)
    occupancy = prediction.branch_occupancy(true, a0, .01)
    assert min(occupancy) == 1., "Prediction demo must stay on the matched branches"
    d["pred_true"] = true
    d["pred_matched"] = prediction.linear_prediction(THETA, a0, .01, 50, .01,
                                                    "activating", "shortening")
    d["pred_wrong"] = prediction.linear_prediction(THETA, a0, .01, 50, .01,
                                                  "deactivating", "lengthening")
    metrics["prediction"] = {"excitation_step": .01, "horizon_s": .5,
        "wrong_error_deg": float(np.rad2deg(abs(d["pred_wrong"][-1, 0] - true[-1, 0]))),
        "true_swing_deg": float(np.rad2deg(true[-1, 0] - THETA))}
    d["pair_at"] = np.linspace(0, .5, 101)
    d["pair_ab"] = np.array([antagonist.find_equilibrium(THETA, a) for a in d["pair_at"]])
    d["pair_ratio"] = np.array([antagonist.damping_contributions(THETA, b, a)[2]
                               for b, a in zip(d["pair_ab"], d["pair_at"])])
    metrics["pair"] = {"triceps_moment_arm_scale": 1.,
        "balanced_activation": antagonist.balanced_cocontraction(THETA)}
    for v in d.values():
        assert np.all(np.isfinite(v)), "Non-finite simulation data"
    np.savez_compressed(cache, **d)
    manifest.write_text(json.dumps(metrics, indent=2) + "\n")
    return d, metrics


def text(fig, x, y, s, size=16, color=INK, **kw):
    return fig.text(x, y, s, fontsize=size, color=color, **kw)


def base(title, subtitle, chapter):
    fig = plt.figure(figsize=(12.8, 7.2), dpi=100)
    # text(fig, .045, .945, "MUSCLE ACTIVATION CONTROL", 11, NAVY, weight="bold")
    # text(fig, .955, .945, "COMPUTATIONAL STUDY  /  SIMULATION", 10, MUTED, ha="right")
    text(fig, .045, .877, title, 27, weight="bold")
    text(fig, .045, .827, subtitle, 13, MUTED)
    text(fig, .045, .036, chapter, 10, MUTED)
    # text(fig, .955, .036, "Hill-type elbow model  |  ICRA 2027", 10, MUTED, ha="right")
    caption = text(fig, .5, .106, "", 16, ha="center")
    return fig, caption


def axes(fig, rect, xlabel, ylabel):
    ax = fig.add_axes(rect)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.tick_params(labelsize=11)
    ax.grid(True)
    return ax


class Arm:
    """Schematic with optional angular magnification and a true-scale inset.

    Main angle = 60 degrees + gain * (physical angle - 60 degrees).
    The target uses the same transform. All numerical readouts remain physical.
    """
    def __init__(self, fig, rect, color=NAVY, pair=False, gain=1.):
        self.gain = gain
        self.ax = fig.add_axes(rect)
        self.ax.set(xlim=(-.5, 1.15), ylim=(-1.05, 1.0), aspect="equal")
        self.ax.axis("off")
        self.ax.plot([0, 0], [0, .83], color=MUTED, lw=15, solid_capstyle="round")
        self.ghost, = self.ax.plot([], [], color=INK, lw=2, ls="--", alpha=.8, zorder=6)
        self.ghost.set_path_effects([pe.Stroke(linewidth=4, foreground=BG), pe.Normal()])
        self.forearm, = self.ax.plot([], [], color=color, lw=14, solid_capstyle="round")
        self.sweep, = self.ax.plot([], [], color=color, lw=2, alpha=.6, ls=":")
        self.tip, = self.ax.plot([], [], "o", color=color, ms=6, zorder=7)
        self.flexor, = self.ax.plot([], [], color=RED, lw=7, solid_capstyle="round")
        self.extensor, = self.ax.plot([], [], color=MAROON, lw=7, solid_capstyle="round")
        self.extensor.set_visible(pair)
        self.ax.add_patch(Circle((0, 0), .075, facecolor=BG, edgecolor=INK, lw=2, zorder=8))
        self.readout = self.ax.text(.5, .01, "", transform=self.ax.transAxes,
                                    ha="center", fontsize=12, color=INK)
        # if gain != 1:
        #     self.ax.plot([.78, .78], [.48, .72], color=MUTED, lw=4, solid_capstyle="round")
        #     self.true_forearm, = self.ax.plot([], [], color=color, lw=4, solid_capstyle="round")
        #     self.true_target, = self.ax.plot([], [], color=INK, lw=1, ls="--", zorder=6)
            # self.ax.text(.78, .84, "Actual 1×", ha="center", fontsize=8, color=MUTED)

    def update(self, angle, activation, target=None, antagonist_a=0):
        display_angle = THETA + self.gain * (angle - THETA)
        x, y = .95 * np.sin(display_angle), -.95 * np.cos(display_angle)
        self.forearm.set_data([0, x], [0, y])
        self.tip.set_data([x], [y])
        sweep = np.linspace(THETA, display_angle, 40)
        self.sweep.set_data(.95*np.sin(sweep), -.95*np.cos(sweep))
        self.flexor.set_data([.10, .32*x], [.60, .32*y])
        self.flexor.set_alpha(.2 + .8 * np.clip(activation, 0, 1))
        self.extensor.set_data([-.12, -.15, .15*x], [.60, -.09, .15*y])
        self.extensor.set_alpha(.2 + .8 * np.clip(antagonist_a, 0, 1))
        if target is not None:
            displayed_target = THETA + self.gain * (target - THETA)
            self.ghost.set_data([0, .95*np.sin(displayed_target)], [0, -.95*np.cos(displayed_target)])
        # if self.gain != 1:
        #     # self.true_forearm.set_data([.78, .78+.24*np.sin(angle)], [.48, .48-.24*np.cos(angle)])
        #     if target is not None:
        #         self.true_target.set_data([.78, .78+.24*np.sin(target)], [.48, .48-.24*np.cos(target)])
        self.readout.set_text(f"\n{np.rad2deg(angle):.2f}°   |   a = {activation:.3f}")


def intro(d, m):
    fig, cap = base("",
        "", "")
    arm = Arm(fig, [.06, .24, .37, .54], gain=MOTION_GAIN)
    # text(fig, .46, .69, "A simulated muscle-driven elbow", 22, weight="bold")
    # for y, s in zip([.60, .52, .44], ["",
    #         "",
    #         ""]):
    #     text(fig, .46, y, s, 17, MUTED)
    # text(fig, .46, .31, "", 13, MUTED)
    # text(fig, .46, .23, "", 13, MUTED)
    # cap.set_text("")
    def update(p):
        t = 4.3*p
        x = np.array([np.interp(t, d["current_state_t"], d["current_state_x"][:, j]) for j in range(3)])
        arm.update(x[0], x[2], float(reference(t)))
    return fig, update


def switches(d, m):
    fig, cap = base("", "", "")
    a0 = float(d["a0"])
    ax1 = axes(fig, [.085, .32, .365, .40], "Excitation − activation, u − a", "Activation rate (1/s)")
    ax2 = axes(fig, [.59, .32, .365, .40], "Normalized fibre velocity (lengths/s)", "Active force / activation")
    xx = np.linspace(-.035, .035, 301)
    vv = np.linspace(-.5, .5, 301)
    for ax, x, y in [(ax1, xx, muscle.activation_derivative(a0, a0+xx)),
                     (ax2, vv, muscle.contractile_force(a0, 1., vv)/a0)]:
        for mask, col, ls in [(x <= 0, RED, "--"), (x >= 0, NAVY, "-")]:
            ax.plot(x[mask], y[mask], color=col, ls=ls, lw=3)
        ax.axvline(0, color=MUTED, ls=":")
    ax1.set_title("Activation switch", loc="left", fontsize=17)
    ax2.set_title("Force–velocity switch", loc="left", fontsize=17)
    dot1, = ax1.plot([], [], "o", ms=10, color=INK)
    dot2, = ax2.plot([], [], "o", ms=10, color=INK)
    text(fig, .27, .215, "u = a at equilibrium", 18, NAVY, ha="center")
    text(fig, .77, .215, "v = 0 at equilibrium", 18, NAVY, ha="center")
    # cap.set_text("Two switching surfaces → four candidate branch linearizations.")
    def update(p):
        q = -np.cos(4*np.pi*min(p/.75, 1)) if p < .75 else 0.
        x, v = .03*q, .45*q
        dot1.set_data([x], [muscle.activation_derivative(a0, a0+x)])
        dot2.set_data([v], [muscle.contractile_force(a0, 1., v)/a0])
    return fig, update


def prediction_scene(d, m):
    fig, cap = base("",
        "",
        "")
    ax = axes(fig, [.095, .29, .61, .46], "Simulation time (ms)", "Angle change (deg)")
    tt = np.arange(51)*10
    ys = [np.rad2deg(d[k][:, 0]-THETA) for k in ["pred_true", "pred_matched", "pred_wrong"]]
    lines = []
    for y, col, ls, lab in zip(ys, [MAROON, NAVY, RED], ["-", "--", ":"],
                              ["Nonlinear simulation", "Matched branches", "Both branches wrong"]):
        line, = ax.plot([], [], color=col, ls=ls, lw=3, label=lab)
        lines.append(line)
    ax.set(xlim=(0, 500), ylim=(-.04, max(max(y) for y in ys)*1.13))
    ax.legend(loc="upper left", frameon=False, fontsize=12, labelcolor=INK)
    value = text(fig, .82, .58, "", 30, RED, ha="center", weight="bold")
    # text(fig, .82, .49, "wrong-branch error\n(degrees)", 14, MUTED, ha="center")
    # text(fig, .82, .34, "Same initial state\nSame excitation", 14, MUTED, ha="center")
    # cap.set_text("The matched local model stays close; the wrong branches underpredict motion.")
    # text(fig, .5, .185, "0.5 s of simulated motion, slowed for inspection", 12, MUTED, ha="center")
    def update(p):
        n = min(51, int(min(p/.84, 1)*50)+1)
        for ln, y in zip(lines, ys):
            ln.set_data(tt[:n], y[:n])
        value.set_text(f"{abs(ys[2][n-1]-ys[0][n-1]):.3f}°")
    return fig, update


def control(d, m, name):
    spec = {"raise": (.10, .75, "", "step_up"), # Raising: the blended model overshoots., 
            "lower": (1.15, 1.80, "", "step_down"), # Lowering: the linearization point matters., step_down
            "track": (2.3, 4.3, "", "sinusoid")} # Tracking: compare the full motion., sinusoid
    t0, t1, title, metric = spec[name]
    fig, cap = base(title,
        "",
        "")
    arms, numbers = [], []
    for i, (col, lab, mode) in enumerate(zip(COLORS, LABELS, MODES)):
        x = .04 + i*.24
        label_color = INK if col == STEEL else col
        text(fig, x+.105, .737, lab, 15, label_color, ha="center", weight="bold")
        arms.append(Arm(fig, [x, .44, .215, .285], col, gain=MOTION_GAIN))
        sm = m["controllers"][mode][metric]
        val = sm["rms_track_deg"] if name == "track" else sm["overshoot_deg"]
        label = "RMS error" if name == "track" else "Overshoot"
        numbers.append(text(fig, x+.105, .407, f"{label}: {val:.3f}°", 14, label_color, ha="center"))
    ax = axes(fig, [.075, .205, .405, .145], "Simulation time (s)", "Angle (deg)")
    err = axes(fig, [.59, .205, .365, .145], "Simulation time (s)", "Error (deg)")
    tt = d["fixed_naive_t"]
    mask = (tt >= t0) & (tt <= t1)
    ax.plot(tt[mask], np.rad2deg(d["fixed_naive_theta_ref"][mask]), color=INK, ls="--", lw=1.5)
    err.axhline(0, color=INK, ls="--", lw=1)
    ax.set(xlim=(t0, t1), ylim=(56.4, 64.1) if name != "raise" else (59.5, 64.1))
    err.set(xlim=(t0, t1), ylim=(-.85, .85) if name == "track" else (-3.3, 3.3))
    lines = [(ax.plot([], [], color=c, ls=s, lw=2)[0], err.plot([], [], color=c, ls=s, lw=2)[0])
             for c, s in zip(COLORS, STYLES)]
    cursors = [a.axvline(t0, color=MUTED, alpha=.6) for a in [ax, err]]
    # For step scenes, a separate inset magnifies the settled transient error.
    if name != "track":
        err.set_ylim(-.8, .8)
        err.set_title("Magnified error; initial step error outside scale", fontsize=10, loc="left")
    clock = text(fig, .5, .795, "", 11, MUTED, ha="center")
    # text(fig, .5, .765, "ARM MOTION ×12 ABOUT 60°  •  INSETS: ACTUAL SCALE  •  ALL NUMBERS: ACTUAL DEGREES",
    #      10, MAROON, ha="center", weight="bold")
    captions = {"raise": "", # Feedback maintains stability; the transient still depends on the model.
        "lower": "", # Reading the branch at the reference is not the same as relinearizing at the state.
        "track": ""} # At these baseline weights, state relinearization has the lowest tracking RMS.
    cap.set_text(captions[name])
    cap.set_fontsize(14)
    def update(p):
        q = min(p/.85, 1.)
        now = t0 + (t1-t0)*q
        ref = float(reference(now))
        # clock.set_text(f"t = {now:.2f} s  |  slowed playback  |  dashed: reference  |  schematic geometry")
        for arm, mode, (line, eline) in zip(arms, MODES, lines):
            x = d[f"{mode}_x"]
            arm.update(np.interp(now, tt, x[:, 0]), np.interp(now, tt, x[:, 2]), ref)
            show = mask & (tt <= now)
            line.set_data(tt[show], np.rad2deg(x[show, 0]))
            eline.set_data(tt[show], np.rad2deg(x[show, 0]-d[f"{mode}_theta_ref"][show]))
        for c in cursors:
            c.set_xdata([now, now])
    return fig, update


def pair(d, m):
    fig, cap = base("",
        "",
        "")
    arm = Arm(fig, [.035, .28, .28, .47], pair=True)
    arm.readout.set_visible(False)
    ax = axes(fig, [.43, .31, .505, .42], "Triceps equilibrium activation", "Velocity slope ratio")
    ax.set(xlim=(0, .5), ylim=(.85, 2.1))
    ax.axhline(1, color=INK, ls="--", lw=1)
    ax.text(.49, 1.045, "equal slopes", color=MUTED, ha="right", fontsize=12)
    ax.plot(d["pair_at"], d["pair_ratio"], color=GRID, lw=2)
    line, = ax.plot([], [], color=NAVY, lw=3)
    dot, = ax.plot([], [], "o", color=NAVY, ms=9)
    val = text(fig, .17, .25, "", 14, MUTED, ha="center")
    text(fig, .68, .21, "Extending-side damping / flexing-side damping", 12, MUTED, ha="center")
    # cap.set_text("Cancellation depends on muscle geometry; activation switches remain.")
    cap.set_fontsize(15)
    def update(p):
        n = min(100, int(min(p/.86, 1)*100))
        ab, at, ratio = d["pair_ab"][n], d["pair_at"][n], d["pair_ratio"][n]
        arm.update(THETA, ab, antagonist_a=at)
        val.set_text(f"Flexor a = {ab:.3f}\nExtensor a = {at:.3f}")
        line.set_data(d["pair_at"][:n+1], d["pair_ratio"][:n+1])
        dot.set_data([at], [ratio])
    return fig, update


def end(d, m):
    # fig, cap = base("Read the branch. Update the operating point.",
    #     "A computational result for linearization-based muscle control", "06 / TAKEAWAY")
    # for y, n, title, body in [(.67, "01", "Equilibrium sits on switching surfaces.", "A single blended Jacobian hides the branch structure."),
    #         (.49, "02", "Prediction mismatch survives as transient error.", "The controller comparisons use the same plant, weights and horizon."),
    #         (.31, "03", "Relinearize at the measured state.", "Use explicit branches; include transport delay in the prediction model.")]:
    #     text(fig, .065, y, n, 28, NAVY, weight="bold")
    #     text(fig, .14, y+.01, title, 21, weight="bold")
    #     text(fig, .14, y-.045, body, 15, MUTED)
    # cap.set_text("Simulation study • baseline weighting shown • hardware validation remains future work")
    # cap.set_fontsize(13)
    return 
# fig, lambda p: None


def make_scene(name, d, m):
    if name in ("raise", "lower", "track"):
        return control(d, m, name)
    return {"intro": intro, "switches": switches, "prediction": prediction_scene,
            "pair": pair}[name](d, m)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview", action="store_true", help="Only generate stills and contact sheet")
    parser.add_argument("--recompute", action="store_true")
    parser.add_argument("--fps", type=int, default=24)
    args = parser.parse_args()
    if args.fps < 20:
        parser.error("Submission video requires at least 20 fps")
    d, m = data(args.recompute)
    previews = []
    for name, duration in SCENES:
        fig, update = make_scene(name, d, m)
        update(.78)
        path = OUT / f"preview_{name}.png"
        fig.savefig(path, dpi=100)
        previews.append(path)
        plt.close(fig)
    sheet, axs = plt.subplots(4, 2, figsize=(16, 18), dpi=100)
    for ax, path in zip(axs.flat, previews):
        ax.imshow(plt.imread(path))
        ax.axis("off")
    sheet.subplots_adjust(0, 0, 1, 1, .015, .015)
    sheet.savefig(OUT / "storyboard.png")
    plt.close(sheet)
    if args.preview:
        print(f"Previews saved: {OUT}", flush=True)
        return
    scene_dir = OUT / "scenes"
    scene_dir.mkdir(parents=True, exist_ok=True)
    for idx, (name, duration) in enumerate(SCENES, start=1):
        video = scene_dir / f"{idx:02d}_{name}.mp4"
        writer = FFMpegWriter(fps=args.fps, codec="libx264", bitrate=850,
            metadata={"title": f"Hill-type muscle control: {name}",
                      "comment": "Simulation trajectories from the accompanying computational model"},
            extra_args=["-pix_fmt", "yuv420p", "-preset", "medium", "-maxrate", "950k",
                        "-bufsize", "1900k", "-movflags", "+faststart"])
        fig, update = make_scene(name, d, m)
        try:
            print(f"Rendering {name}: {duration} s", flush=True)
            with writer.saving(fig, str(video), dpi=100):
                for frame in range(duration*args.fps):
                    update(frame / max(1, duration*args.fps-1))
                    writer.grab_frame()
        finally:
            plt.close(fig)
        assert video.stat().st_size < 20_000_000, f"{video.name} exceeds the 20 MB submission limit"
        print(f"Saved {video} ({video.stat().st_size/1e6:.2f} MB)", flush=True)


if __name__ == "__main__":
    main()
