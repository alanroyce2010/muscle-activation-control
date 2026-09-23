"""Diagnostic (2026-09-11): eigenvalues of the branch Jacobians at the
gravity-only equilibria, old (by-eye) vs corrected moment-arm table, with
one-sided theta derivatives to expose the passive force-length kink at
l_bar = 1 (65 deg by construction) and the fold of the equilibrium family."""
import importlib.util, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from muscle_activation_control import joint as joint_new, linearize, muscle

spec = importlib.util.spec_from_file_location("joint_old", ROOT / "scripts" / "_joint_old_table.py")
joint_old = importlib.util.module_from_spec(spec); spec.loader.exec_module(joint_old)

def jac(jm, th, a, act, vel, side):
    """3x3 Jacobian; theta column one-sided (side=+1 forward, -1 backward, 0 central)."""
    linearize.joint = jm
    A, _ = linearize.linearize(th, a, act, vel)
    h = 1e-6
    f = lambda t: (jm.muscle_torque(t, 0.0, a) - jm.gravity_torque(t)) / jm.FOREARM_HAND_I_ELBOW
    if side == +1: A[1, 0] = (f(th + h) - f(th)) / h
    elif side == -1: A[1, 0] = (f(th) - f(th - h)) / h
    return A

for name, jm in (("OLD by-eye table", joint_old), ("NEW pixel table", joint_new)):
    linearize.joint = jm
    print(f"===== {name} =====")
    for deg in (45, 55, 60, 62, 64, 65, 66, 68, 70, 75, 85):
        th = np.deg2rad(deg); a = linearize.find_equilibrium(th)
        lb = jm.fiber_length_norm(th)
        out = []
        for side in (-1, +1):
            eg = []
            for act, vel in linearize.BRANCH_COMBINATIONS:
                ev = np.linalg.eigvals(jac(jm, th, a, act, vel, side))
                eg.append(np.max(np.real(ev)))
            out.append((min(eg), max(eg)))
        print(f"{deg:3d} deg  a*={a:.4f}  l_bar={lb:.4f}  slowest eig (backward theta) {out[0][1]:+.3f}  (forward theta) {out[1][1]:+.3f} /s")
    th = np.deg2rad(65.0); a = linearize.find_equilibrium(th)
    for side, lab in ((-1, "backward (extension side)"), (+1, "forward (flexion side)"), (0, "central (what linearize.py uses)")):
        A = jac(jm, th, a, "activating", "shortening", side)
        print(f"  65 deg  d(thetaddot)/dtheta {lab}: {A[1,0]:+.3f} /s^2; eigs act/short: {np.round(np.sort(np.real(np.linalg.eigvals(A))),3)}")
    grid = np.deg2rad(np.linspace(20, 120, 1001))
    astar = np.array([linearize.find_equilibrium(t) for t in grid])
    print(f"  gravity-only a*(theta) peaks at {np.rad2deg(grid[np.argmax(astar)]):.1f} deg (a*={astar.max():.4f}) -> fold of the equilibrium family")
linearize.joint = joint_new
