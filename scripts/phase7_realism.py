"""Phase 7 sections 5-6 (docs/theory/phase7_limitations.md): closed-loop
realism sweep and hybrid mode-sequence MPC on the Phase 6b task (steps
60 -> 63 -> 57 deg, then a 0.5 Hz +/-3 deg sinusoid about 60 deg).

One factor at a time from the Phase 6b baseline: control rate, angle noise
with estimated activation, plant-model activation mismatch, transport
delay. Rules: naive fixed, blended, sign test, and (baseline and one noisy
case) the hybrid enumeration MPC. Runs in parallel; writes
results/phase7_realism.json."""

import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from muscle_activation_control import linearize, mpc  # noqa: E402
from phase6b_closed_loop import THETA0, T_END, reference, segment_metrics  # noqa: E402

RULES = ("fixed_naive", "blended", "tracking")
# `--only RULE` runs just that rule over the same conditions (no hybrid jobs)
# and writes results/phase7_realism_<RULE>.json; used to add relinearization
# ("current_state") after the red-team review asked for it.
ONLY = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None
if ONLY:
    RULES = (ONLY,)


def jobs():
    js = []
    for ts in (0.005, 0.010, 0.020, 0.040):
        for r in RULES:
            js.append(dict(group="rate", rule=r, t_s=ts))
    for sigma in (0.0, 0.05, 0.2):
        for seed in (range(1) if sigma == 0.0 else range(5)):
            for r in RULES:
                js.append(dict(group="noise", rule=r, noise=sigma, est=True, seed=seed, vel_diff=True))
    for name, pt in (("tau x0.75", (0.01125, 0.0375)), ("tau x1.25", (0.01875, 0.0625)), ("old-adult plant", (0.015, 0.060))):
        for r in RULES:
            js.append(dict(group="mismatch", rule=r, plant=name, plant_tau=pt, est=True))
    for d in (2, 5):
        for r in RULES:
            js.append(dict(group="delay", rule=r, delay=d))
    if ONLY:
        return js
    js.append(dict(group="hybrid", rule="hybrid"))
    js.append(dict(group="hybrid", rule="tracking"))
    js.append(dict(group="hybrid-noise", rule="hybrid", noise=0.05, est=True, seed=0, vel_diff=True))
    js.append(dict(group="hybrid-noise", rule="tracking", noise=0.05, est=True, seed=0, vel_diff=True))
    return js


def run(job):
    t0 = time.time()
    ts = job.get("t_s", 0.010)
    ctrl = mpc.BranchAwareMPC(job["rule"], t_s=ts, horizon=int(round(0.2 / ts)))
    a0 = linearize.find_equilibrium(THETA0)
    log = mpc.simulate_realistic(ctrl, reference, T_END, [THETA0, 0.0, a0],
                                 noise_deg=job.get("noise", 0.0), estimate_activation=job.get("est", False),
                                 plant_tau=job.get("plant_tau"), delay_steps=job.get("delay", 0),
                                 seed=job.get("seed", 0), velocity_from_angle=job.get("vel_diff", False))
    out = dict(job)
    out.pop("plant_tau", None)
    out["diverged"] = log["diverged"]
    out["seconds"] = time.time() - t0
    if not log["diverged"]:
        seg = segment_metrics(log)
        out.update(up_overshoot=seg["step_up"]["overshoot_deg"], up_settle=seg["step_up"]["settle_s"],
                   dn_overshoot=seg["step_down"]["overshoot_deg"], dn_settle=seg["step_down"]["settle_s"],
                   sin_rms=seg["sinusoid"]["rms_track_deg"], rms=mpc.metrics(log)["rms_track_deg"],
                   branch_changes=int(sum(b1 != b2 for b1, b2 in zip(log["branch"][:-1], log["branch"][1:]))))
        if job["rule"] == "hybrid":
            out["hyb_violation_steps"] = int(np.sum(np.array(ctrl.hyb_violations) > 0))
    return out


def main():
    js = jobs()
    print(f"{len(js)} closed-loop runs")
    with Pool(min(20, len(js))) as pool:
        res = pool.map(run, js, chunksize=1)
    out_name = f"phase7_realism_{ONLY}.json" if ONLY else "phase7_realism.json"
    json.dump(res, open(ROOT / "results" / out_name, "w"), indent=1, default=float)

    def agg(rows, key):
        v = [r[key] for r in rows if not r["diverged"] and r.get(key) is not None]
        return (np.mean(v), np.std(v)) if v else (np.nan, np.nan)

    def show(label, rows):
        n_div = sum(r["diverged"] for r in rows)
        ov, _ = agg(rows, "up_overshoot"); dv, _ = agg(rows, "dn_overshoot")
        st, _ = agg(rows, "up_settle"); sr, ss = agg(rows, "sin_rms"); bc, _ = agg(rows, "branch_changes")
        print(f"  {label:34s} up ovs {ov:6.3f}  settle {st:5.2f}  dn ovs {dv:6.3f}  sin RMS {sr:6.3f} +/- {ss:5.3f}"
              f"  branch changes {bc:6.1f}  diverged {n_div}/{len(rows)}")

    for group in (("rate", "noise", "mismatch", "delay") if ONLY else ("rate", "noise", "mismatch", "delay", "hybrid", "hybrid-noise")):
        print(f"=== {group} ===")
        rows = [r for r in res if r["group"] == group]
        keys = sorted({(r["rule"], r.get("t_s"), r.get("noise"), r.get("plant"), r.get("delay")) for r in rows},
                      key=lambda k: tuple(str(x) for x in k))
        for k in keys:
            sub = [r for r in rows if (r["rule"], r.get("t_s"), r.get("noise"), r.get("plant"), r.get("delay")) == k]
            lab = " ".join(str(x) for x in k if x is not None)
            show(lab, sub)
        for r in rows:
            if r["rule"] == "hybrid" and not r["diverged"]:
                print(f"  hybrid: steps whose applied plan had consistency violations: {r['hyb_violation_steps']}; "
                      f"wall time {r['seconds']:.0f} s")


if __name__ == "__main__":
    main()
