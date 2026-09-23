# Phase 4: Does the Branch Ambiguity Matter for Control?

> **Revision note (2026-09-11, Phase 7):** the table below was computed
> at the old 65 deg equilibrium with the by-eye moment-arm table (and the
> "13-88%" summary overstated its own table, whose maximum was 73%). At
> the corrected 60 deg posture (`phase7_limitations.md` section 9.1) the
> activating-step table reads, as % of the true response (matched / wrong
> activation / wrong velocity / both): 50 ms 1.1 / 71.1 / 13.8 / 74.9,
> 100 ms 1.3 / 57.3 / 24.1 / 66.9, 200 ms 1.2 / 37.6 / 35.1 / 58.4, 500 ms
> 0.9 / 13.7 / 42.9 / 50.7: **wrong-branch range 13.7-74.9%**.
> Deactivating steps: wrong activation 241% at 50 ms. The conclusion
> ("it matters") is unchanged.

CLAUDE.md section 4.2 step 2. Builds on Phase 2 (`phase2_linearization.md`)
and Phase 3 (`phase3_bimodal_formalization.md`) -- the branch ambiguity is
established and formalized; this asks whether it's big enough to matter
for a realistic control design, or a curiosity that washes out. Code:
`scripts/phase4_quantify_branch_consequence.py`,
`tests/test_quantify.py` (1/1 pass, part of 35/35 total).

## 1. Setup

Scenario: a successive-linearization MPC's internal prediction model
typically linearizes once per control step and holds that linear model
fixed over the prediction horizon (the same simplification
`origami-arm-control`'s `MPCController` uses, CLAUDE.md section 2). Start
exactly at the Phase 2 equilibrium (`theta*=65deg, thetadot*=0,
a*=u*=0.133`), apply a small step change in excitation `delta_u`, and
compare:

- the **true** nonlinear trajectory (`scipy.integrate.solve_ivp` on the
  actual model from `joint.py`)
- the **discrete-time linear prediction** using each of the 4 branch
  combinations' `A_d, B_d` (ZOH at `T_s=10ms`, from Phase 2), held fixed
  and propagated forward `k` steps: `dx[k+1] = A_d@dx[k] + B_d*delta_u`

`delta_u=+0.01` keeps the true trajectory on the `activating/shortening`
branch throughout (excitation stays above activation, joint keeps
flexing); `delta_u=-0.01` keeps it on `deactivating/lengthening`
throughout. This lets "matched branch" vs. "mismatched branch" be
compared unambiguously, since the true trajectory doesn't cross either
switching manifold during the horizon.

## 2. Primary result: small-signal regime (delta_u=0.01)

This is the regime linearization is *supposed* to be valid in -- small
perturbations, moderate horizons. Final-time `theta` prediction error
(degrees), `delta_u=+0.01`:

| horizon | true swing | matched branch | wrong-activation | wrong-velocity | both wrong |
|---|---|---|---|---|---|
| 50 ms  | 0.105 deg | 0.00087 | 0.073 (70% of swing) | 0.013 (13%) | 0.077 (73%) |
| 100 ms | 0.434 deg | 0.00389 | 0.243 (56%) | 0.099 (23%) | 0.283 (65%) |
| 200 ms | 1.451 deg | 0.00800 | 0.538 (37%) | 0.500 (34%) | 0.833 (57%) |
| 500 ms | 5.397 deg | 0.06499 | 0.875 (16%) | 2.469 (46%) | 2.873 (53%) |

(`delta_u=-0.01` gives the mirror-image pattern -- see script output,
same qualitative conclusion.)

**The branch ambiguity matters.** Across every horizon tested, using the
wrong activation branch, wrong velocity branch, or both introduces a
prediction error that is a **substantial fraction of the true signal
itself** -- ranging from 13% to 88% of the actual swing being predicted,
never negligible. The matched branch's own error (from linearizing at
all, not from branch choice) stays 1-2 orders of magnitude smaller than
any mismatched branch's error at every horizon tested. See
`results/phase4_activating_step_small.png` for the trajectory plot -- the
matched branch (blue) tracks the true nonlinear curve (black) almost
exactly; all three mismatched branches visibly diverge within the first
50-100ms.

**A secondary, unplanned finding**: the *relative* importance of the two
kinks shifts with horizon. Wrong-activation error is proportionally
largest at short horizons (70% at 50ms) and shrinks in relative terms by
500ms (16%); wrong-velocity error does the opposite (13% at 50ms, growing
to 46% at 500ms). This means the larger raw eigenvalue spread found in
Phase 2 (~6.8x for activation vs. exactly 2x for force-velocity) does
**not** straightforwardly predict which kink matters more for a given
control-relevant output (joint angle) at a given horizon -- the
force-velocity kink enters the mechanical equation (and hence `theta`)
directly, while the activation kink's effect is mediated through `a` and
then through muscle force, a more indirect path whose relative
contribution to `theta` error changes over the horizon. Not fully
explained here (would need tracing the closed-form step response of each
branch's transfer function from `u` to `theta`) -- flagged as a specific,
bounded follow-up if there's room for it, not blocking the headline
finding.

## 3. Contrast: larger-signal regime (delta_u=0.05)

Run for comparison (`results/phase4_*_large.png`), not the headline
result: at this magnitude the true trajectory swings substantially
(up to ~26 deg over 500ms in an earlier, longer-horizon run at this same
magnitude), which is large enough that **generic one-shot-linearization
breakdown** (the linear model no longer tracking the true nonlinear
trajectory at all, for any branch) becomes entangled with the
branch-choice question this section is trying to isolate. At
`delta_u=-0.05`, the "matched" branch (`deactivating/lengthening`) was
no longer even the best predictor at 200ms (`deactivating/shortening`
did better, 0.13 deg vs. 1.60 deg) -- a sign that at this magnitude the
true trajectory's local branch assignment may not hold for the whole
held-fixed horizon the way it reliably does at `delta_u=0.01`. **This is
why section 2's small-signal result, not this one, is the primary
finding** -- it isolates the effect actually being asked about.

## 4. Answer to the section 1 question

**The non-smoothness matters.** In the regime where equilibrium
linearization is normally trusted (small perturbations, the same
`T_s=10ms` found well-conditioned in Phase 2), a successive-linearization
MPC's internal prediction model -- if it silently picks one branch at
this equilibrium, e.g. via whatever a naive `u > a` floating-point
comparison happens to return, or via central finite differences blending
across the kink the way `linearize.py`'s docstring warned against --
would carry a prediction error that is a large fraction of the actual
signal, at every horizon tested. This is not a mathematical curiosity;
it is a concrete, quantified reason a controller designer working with
Hill-type muscle models needs to track which side of each switch the
current operating point is on, something the smooth pneumatic pressure
model in `origami-arm-control` never required.
