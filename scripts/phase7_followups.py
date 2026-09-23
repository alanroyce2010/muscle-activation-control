"""Phase 7 section 10 (docs/theory/phase7_limitations.md): follow-up runs
after the review of the ICRA draft, on the Phase 6b task (steps 60 -> 63
-> 57 deg, then a 0.5 Hz +/-3 deg sinusoid about 60 deg).

  - delay_comp: transport delay of 20 or 50 ms WITH a predictor (the
    measured state rolled forward through the queued inputs by the
    controller's nominal model, reference preview shifted to match), on
    the nominal plant and on a plant whose activation time constants are
    x1.25 (the predictor's model is then wrong);
  - weights: the undelayed baseline under three other cost weightings.

Rules: naive, blended, sign test, relinearize, and blended_act (the sign
test with only the activation row averaged). blended_act under the full
realism sweep runs through `phase7_realism.py --only blended_act`.
Writes results/phase7_followups.json."""

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

RULES = ("fixed_naive", "blended", "tracking", "current_state", "blended_act")
# (q_theta, q_thetadot, r); "base" is the weighting of every earlier run
WEIGHTS = {"base": (1000.0, 1.0, 200.0), "R x0.25": (1000.0, 1.0, 50.0),
           "R x4": (1000.0, 1.0, 800.0), "q_thetadot x20": (1000.0, 20.0, 200.0)}


def jobs():
    js = []
    for d in (2, 5):
        for plant, pt in (("nominal", None), ("tau x1.25", (0.01875, 0.0625))):
            for r in RULES:
                js.append(dict(group="delay_comp", rule=r, delay=d, plant=plant, plant_tau=pt))
    for w in WEIGHTS:
        if w != "base":
            for r in RULES:
                js.append(dict(group="weights", rule=r, weights=w))
    return js


def run(job):
    t0 = time.time()
    qt, qd, r = WEIGHTS[job.get("weights", "base")]
    ctrl = mpc.BranchAwareMPC(job["rule"], q_theta=qt, q_thetadot=qd, r=r)
    a0 = linearize.find_equilibrium(THETA0)
    log = mpc.simulate_realistic(ctrl, reference, T_END, [THETA0, 0.0, a0], plant_tau=job.get("plant_tau"),
                                 delay_steps=job.get("delay", 0),
                                 compensate_delay=job["group"] == "delay_comp")
    out = {k: v for k, v in job.items() if k != "plant_tau"}
    out["diverged"] = log["diverged"]
    out["seconds"] = time.time() - t0
    if not log["diverged"]:
        seg = segment_metrics(log)
        out.update(up_overshoot=seg["step_up"]["overshoot_deg"], up_settle=seg["step_up"]["settle_s"],
                   dn_overshoot=seg["step_down"]["overshoot_deg"], dn_settle=seg["step_down"]["settle_s"],
                   sin_rms=seg["sinusoid"]["rms_track_deg"],
                   branch_changes=int(sum(b1 != b2 for b1, b2 in zip(log["branch"][:-1], log["branch"][1:]))))
    return out


def main():
    js = jobs()
    print(f"{len(js)} closed-loop runs")
    with Pool(min(20, len(js))) as pool:
        res = pool.map(run, js, chunksize=1)
    json.dump(res, open(ROOT / "results" / "phase7_followups.json", "w"), indent=1, default=float)
    for group, key in (("delay_comp", lambda r: f"delay {10 * r['delay']} ms, {r['plant']}"),
                       ("weights", lambda r: r["weights"])):
        print(f"=== {group} ===")
        rows = [r for r in res if r["group"] == group]
        for cond in dict.fromkeys(key(r) for r in rows):
            print(f"  {cond}")
            for r in rows:
                if key(r) != cond:
                    continue
                if r["diverged"]:
                    print(f"    {r['rule']:14s} DIVERGED")
                else:
                    print(f"    {r['rule']:14s} up ovs {r['up_overshoot']:6.3f}  dn ovs {r['dn_overshoot']:6.3f}"
                          f"  sin RMS {r['sin_rms']:6.3f}  branch changes {r['branch_changes']:4d}")


if __name__ == "__main__":
    main()
