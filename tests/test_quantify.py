"""Phase 4 directional sanity check: for a step excitation change that
keeps the true trajectory on one side of the activation switch the whole
time (see scripts/phase4_quantify_branch_consequence.py), the
branch-matched one-shot linear prediction should track the true nonlinear
trajectory better than the mismatched branch's prediction. Not pinning
exact error magnitudes here (those are read off scripts/phase4_*.py's
printed output / plots for the write-up) -- just the qualitative relation
a "does the branch choice matter" claim depends on.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from phase4_quantify_branch_consequence import (  # noqa: E402
    HORIZON_STEPS,
    T_S,
    linear_prediction,
    true_trajectory,
)
from muscle_activation_control import linearize  # noqa: E402

THETA_STAR = np.deg2rad(60.0)


def test_matched_activation_branch_outperforms_mismatched():
    a_star = linearize.find_equilibrium(THETA_STAR)
    delta_u = 0.05  # activating step: true trajectory stays on u>a throughout
    true_traj = true_trajectory(a_star, delta_u, HORIZON_STEPS, T_S)

    matched = linear_prediction(THETA_STAR, a_star, delta_u, HORIZON_STEPS, T_S,
                                  "activating", "shortening")
    mismatched = linear_prediction(THETA_STAR, a_star, delta_u, HORIZON_STEPS, T_S,
                                     "deactivating", "shortening")

    err_matched = abs(matched[-1, 0] - true_traj[-1, 0])
    err_mismatched = abs(mismatched[-1, 0] - true_traj[-1, 0])
    assert err_matched < err_mismatched
