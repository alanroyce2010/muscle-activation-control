"""Pixel-level digitization of Murray, Delp & Buchanan (1995) Fig. 4
(elbow flexion/extension moment arms, forearm neutral), replacing the
original by-eye table (docs/theory/phase7_limitations.md section 1).

Renders the figure page at 600 dpi from the local (git-ignored) PDF,
calibrates each panel's axes from its tick marks by least squares, and
tracks each curve at 1-degree steps by nearest-neighbour continuity on
dark-pixel run centres (dash gaps carried; label text excluded by angle
window). Checks every track against Murray's own Table 2 percent
difference (max - min)/max over 25-120 deg, writes
data/murray1995_fig4_digitized.csv on a 5-degree grid, and saves
overlays to results/digitization/ for visual verification.

Requires pdftoppm (poppler) and the PDF at
docs/literature/murray1995_elbow_moment_arms.pdf.
"""

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parent.parent
PDF = ROOT / "docs" / "literature" / "murray1995_elbow_moment_arms.pdf"
OUT_CSV = ROOT / "data" / "murray1995_fig4_digitized.csv"
OUT_DIR = ROOT / "results" / "digitization"

DARK = 110                      # grey level below which a pixel counts as ink
GRID = np.arange(0, 131, 5)     # output grid, degrees

# panel crops as page fractions (y0, y1, x0, x1) and tick calibration:
# x ticks (column px at 20..120 deg), y ticks (row px, cm values)
PANELS = {
    "model": dict(crop=(0.49, 0.72, 0.10, 0.90),
                  xt=([20, 40, 60, 80, 100, 120], [1178, 1451, 1723, 1999, 2272, 2542]),
                  yt=([6, 4, 2, 0, -2], [251, 462, 667, 871, 1079])),
    "male": dict(crop=(0.08, 0.32, 0.10, 0.90),
                 xt=([20, 40, 60, 80, 100, 120], [1173, 1446, 1718, 1997, 2269, 2543]),
                 yt=([6, 4, 2, 0, -2], [362, 569, 773, 978, 1182])),
}

# (panel, name): seed (deg, cm), valid angle range, excluded angle window
# (label text), window half-width in cm, Murray Table 2 value (%)
TRACKS = {
    ("model", "BRA"): dict(seed=(70, 2.41), rng=(0, 130), excl=(116, 131), tol=0.16, table2=48),
    ("model", "BIC"): dict(seed=(70, 3.03), rng=(0, 130), excl=(116, 131), tol=0.16, table2=43),
    ("model", "BRD"): dict(seed=(90, 5.20), rng=(0, 130), excl=None, tol=0.16, table2=55),
    ("model", "TRI"): dict(seed=(70, -2.62), rng=(0, 130), excl=(116, 131), tol=0.16, table2=25),
    # male specimen: below ~25 deg the BRA track meets the steeply rising
    # BRD curve, so the valid measured range used here is 25-116 deg
    ("male", "BRA"): dict(seed=(60, 2.14), rng=(25, 116), excl=None, tol=0.13, table2=58),
}

# End segments inside label windows, read from cluster listings at
# 120/125/130 deg (see phase7 doc): the curves continue as drawn there.
END_OVERRIDES = {
    ("model", "BRA"): {120: 2.88, 125: 2.88, 130: 2.88},
    ("model", "BIC"): {120: 3.49, 125: 3.33, 130: 3.05},
    ("model", "TRI"): {120: -2.21, 125: -2.16, 130: -2.16},
}


def render_page():
    tmp = Path(tempfile.mkdtemp())
    subprocess.run(["pdftoppm", "-f", "6", "-l", "6", "-r", "600", "-gray", "-png",
                    str(PDF), str(tmp / "p")], check=True)
    f = sorted(tmp.glob("p*.png"))[0]
    return np.array(Image.open(f).convert("L"))


def panel_mask(page, crop):
    H, W = page.shape
    y0, y1, x0, x1 = crop
    return page[int(H * y0):int(H * y1), int(W * x0):int(W * x1)] < DARK


def calibration(spec):
    kx, bx = np.polyfit(*spec["xt"], 1)
    ky, by = np.polyfit(*spec["yt"], 1)
    return kx, bx, ky, by


def run_centres(P, deg, cal, half=2, rmin=40, rmax=1400):
    kx, bx, ky, by = cal
    ax_row = int(round(by))
    c = int(round(kx * deg + bx))
    out = []
    for cc in range(c - half, c + half + 1):
        if not 0 <= cc < P.shape[1]:
            continue
        col = P[rmin:rmax, cc].copy()
        col[ax_row - rmin - 5:ax_row - rmin + 6] = False      # x-axis line
        idx = np.where(col)[0] + rmin
        if len(idx) == 0:
            continue
        runs = np.split(idx, np.where(np.diff(idx) > 2)[0] + 1)
        out += [(r.mean() - by) / ky for r in runs if 2 <= len(r) <= 30]
    return np.array(out)


def track(P, cal, seed, rng, excl, tol):
    pts = {}
    for direction in (+1, -1):
        d, v = seed
        while rng[0] <= d <= rng[1]:
            if excl is None or not (excl[0] <= d <= excl[1]):
                c = run_centres(P, d, cal)
                near = c[np.abs(c - v) < tol] if len(c) else c
                if len(near):
                    v = float(np.median(near))
                    pts[d] = v
            d += direction
    return pts


def main():
    if not PDF.exists():
        sys.exit(f"missing {PDF} (git-ignored literature; see docs/literature/SOURCES.md)")
    page = render_page()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    masks = {p: panel_mask(page, s["crop"]) for p, s in PANELS.items()}
    cals = {p: calibration(s) for p, s in PANELS.items()}
    tables, checks = {}, {}
    for (panel, name), spec in TRACKS.items():
        pts = track(masks[panel], cals[panel], spec["seed"], spec["rng"], spec["excl"], spec["tol"])
        pts.update(END_OVERRIDES.get((panel, name), {}))
        ds = np.array(sorted(pts)); vs = np.array([pts[d] for d in ds])
        m = (ds >= 25) & (ds <= 120)
        a = np.abs(vs[m])
        checks[(panel, name)] = (100 * (a.max() - a.min()) / a.max(), spec["table2"])
        lo, hi = spec["rng"]
        grid = GRID[(GRID >= lo) & (GRID <= hi)]
        tables[(panel, name)] = (grid, np.interp(grid, ds, vs), ds, vs)
        print(f"{panel:5s} {name}: {len(ds)} tracked degrees; Table 2 check "
              f"{checks[(panel, name)][0]:.1f}% vs Murray {spec['table2']}%; "
              f"ROM mean |r| {np.mean(np.abs(np.interp(np.arange(lo, hi + 1), ds, vs))):.2f} cm")
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["panel", "muscle", "elbow_flexion_deg", "moment_arm_cm"])
        for (panel, name), (grid, vals, _, _) in tables.items():
            for g, v in zip(grid, vals):
                w.writerow([panel, name, int(g), f"{v:.3f}"])
    print(f"wrote {OUT_CSV}")
    colors = {"BRA": (0, 170, 0), "BIC": (0, 120, 255), "BRD": (230, 0, 0), "TRI": (200, 0, 200)}
    for panel in PANELS:
        P = masks[panel]; kx, bx, ky, by = cals[panel]
        img = Image.fromarray((~P * 255).astype(np.uint8)).convert("RGB")
        dr = ImageDraw.Draw(img)
        for (pn, name), (_, _, ds, vs) in tables.items():
            if pn != panel:
                continue
            for d, v in zip(ds, vs):
                cx, cy = kx * d + bx, ky * v + by
                dr.ellipse([cx - 5, cy - 5, cx + 5, cy + 5], fill=colors[name])
        img.resize((img.width // 3, img.height // 3)).save(OUT_DIR / f"murray_fig4_{panel}_tracks.png")
    for key, (got, want) in checks.items():
        assert abs(got - want) <= 3.5, f"{key}: Table 2 check {got:.1f}% vs {want}%"
    print("all Table 2 checks within 3.5 percentage points")


if __name__ == "__main__":
    main()
