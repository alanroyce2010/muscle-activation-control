"""Export separate, synchronized arm-only MP4 assets for video editing."""
import json
from pathlib import Path
import zipfile

import render_submission_video as video
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import numpy as np

OUT = video.OUT / "arm_clips"
FPS = 24


def scene(color, pair=False):
    fig = plt.figure(figsize=(7.2, 7.2), dpi=100, facecolor=video.BG)
    arm = video.Arm(fig, [.07, .08, .86, .88], color=color,
                    pair=pair, gain=1 if pair else video.MOTION_GAIN)
    # Keep only the main arm, muscle lines and target: no insets/readouts.
    visible = {arm.ax.lines[0], arm.forearm, arm.flexor, arm.extensor, arm.ghost}
    for line in arm.ax.lines:
        if line not in visible:
            line.set_visible(False)
    for label in arm.ax.texts:
        label.set_visible(False)
    footer = ("Simulation · fixed 60° · co-contraction equilibrium sweep" if pair else
              "Simulation · motion ×12 about 60° · playback 0.1×")
    video.text(fig, .5, .025, footer, 10, video.MUTED, ha="center")
    return fig, arm


def export(name, fig, update, seconds):
    path = OUT / f"{name}.mp4"
    writer = FFMpegWriter(fps=FPS, codec="libx264", bitrate=700,
                         extra_args=["-pix_fmt", "yuv420p", "-movflags", "+faststart"])
    print(f"Rendering {name}: {seconds} s", flush=True)
    with writer.saving(fig, str(path), dpi=100):
        for frame in range(round(seconds*FPS)):
            update(frame / FPS)
            writer.grab_frame()
    plt.close(fig)
    return path


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data, metrics = video.data()
    paths = []
    for mode, color in zip(video.MODES, video.COLORS):
        fig, arm = scene(color)
        t, x = data[f"{mode}_t"], data[f"{mode}_x"]
        def update(playback, arm=arm, t=t, x=x):
            now = playback / 10
            arm.update(np.interp(now, t, x[:, 0]), np.interp(now, t, x[:, 2]),
                       float(video.reference(now)))
        paths.append(export(mode, fig, update, 43))
    fig, arm = scene(video.NAVY, pair=True)
    def pair_update(playback):
        at = .5*min(playback/17, 1.)
        ab = np.interp(at, data["pair_at"], data["pair_ab"])
        arm.update(video.THETA, ab, antagonist_a=at)
    paths.append(export("antagonist_cocontraction", fig, pair_update, 18))
    notes = OUT / "README.txt"
    notes.write_text(
        "ARM-ONLY VIDEO ASSETS\n\n"
        "720 x 720, 24 fps, H.264 MP4, cream background, no audio.\n"
        "Four controller clips share exactly the same timing (43 seconds).\n"
        "Place them at the same start time to compare controllers.\n\n"
        "fixed_naive.mp4: fixed branch\nblended.mp4: blended Jacobian\n"
        "tracking.mp4: sign-test branch at reference equilibrium\n"
        "current_state.mp4: measured-state relinearization\n\n"
        "Playback time = 10 x simulation time. Target changes to 63 degrees\n"
        "at video 00:03, to 57 degrees at 00:13, then to sinusoidal tracking\n"
        "at 00:23. MPC preview can move the arm before a target change.\n"
        "The visible angle is 60 + 12 x (simulated angle - 60) degrees.\n"
        "The dashed target uses the same magnification. This is visual\n"
        "magnification only; model inputs and simulation results are unchanged.\n"
        "Red muscle opacity encodes activation. Anatomy is schematic.\n\n"
        "antagonist_cocontraction.mp4 (18 seconds): fixed true 60-degree\n"
        "posture while flexor and extensor activation rise. This is a sequence\n"
        "of equilibria, not a dynamic movement; opacity changes are intentional.\n\n"
        "Regenerate: .venv/bin/python scripts/render_arm_clips.py\n")
    archive = video.OUT / "arm_animations.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for path in [*paths, notes]:
            z.write(path, path.name)
    print(f"Saved five clips and {archive}", flush=True)


if __name__ == "__main__":
    main()
