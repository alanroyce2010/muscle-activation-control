# Phase 6b: Does the Branch Ambiguity Survive Feedback?

> **Revision note (2026-09-11, Phase 7):** rerun with the corrected
> moment arm and the task recentred on the 60 deg posture (60 -> 63 -> 57
> deg, sinusoid 60 +/- 3 deg; all reference equilibria now open-loop
> stable). Blended overshoot 0.654 deg vs 0.072-0.081 for the
> branch-aware rules (8.1-9.1x), settling 0.23 vs 0.10 s; naive fixed
> sinusoid RMS 0.200 vs sign test 0.153 deg (31% higher; was 38%); the
> sign test still overshoots the down step (0.374 deg); per-step
> iteration fell back on 157 of 430 steps. No variant unstable or offset.
> Phase 7 adds a hybrid enumeration MPC, a realism sweep and a Lyapunov
> certificate (`phase7_limitations.md` sections 6, 5, 9.5).

Tier 2 item T2.b of `docs/completion_plan.md`, and the "closed-loop MPC
test of branch-tracking against a naive fixed linearization" the
submitted abstract promised as future work. Phases 4 and 6a measured
**open-loop** one-shot prediction error. Feedback re-measures the state
every `T_s` and re-solves, so a prediction model that is wrong by 70% at
50 ms might still produce an acceptable closed-loop controller. This
phase builds the controller and checks. Code:
`src/muscle_activation_control/mpc.py`, `scripts/phase6b_closed_loop.py`,
`tests/test_mpc.py`.

## 1. Controller

Successive-linearization MPC over the full state `x=[theta, thetadot,
a]` with scalar input `u` (excitation), mirroring the
`origami-arm-control` `MPCController` QP structure (CLAUDE.md section
2) but with no inner loop and no reduced state -- Phase 2 showed there
is no stiffness to reduce away, so the honest baseline is the plain
full-state version.

At each control step `k`, with measured `x_k`, previous input `u_{k-1}`,
and a previewed reference `theta_ref(t_k + j T_s)`, `j=0..N`:

1. Linearization point `(x_lin, u_lin)`: the **reference equilibrium**
   `x_lin = [theta_ref(t_k), 0, a_ref]`, `u_lin = a_ref`, with
   `a_ref = find_equilibrium(theta_ref)`. This is the standard
   equilibrium-linearization design (the origami controller does the
   same, `_equilibrium_slow(x_d)`), and it is exactly the point at which
   the model is non-smooth: `u_lin = a_ref` sits on the activation switch
   and `thetadot = 0` on the force-velocity switch. **The controller
   cannot avoid choosing a branch here.**
2. `(A, B)` on the chosen branch, ZOH-discretized at `T_s` to `(A_d, B_d)`.
3. QP: `min sum_j (x_j - x_ref_j)^T Q (x_j - x_ref_j) + R (u_j - u_lin)^2
   + terminal DARE cost`, subject to
   `x_{j+1} = x_lin + A_d (x_j - x_lin) + B_d (u_j - u_lin)`,
   `x_0 = x_k`, `0 <= u_j <= 1`. Apply `u_0`.

`T_s = 10 ms`, `N = 20` (200 ms horizon), `Q = diag(q_theta, q_thetadot,
0)`, `R` tuned once (section 3) and then held identical across variants,
so the only thing that differs between variants is step 2.

## 2. Variants (the experiment)

| variant | branch used at the reference equilibrium |
|---|---|
| `fixed_naive` | whatever the model's own comparisons return *at exact equality*: `u > a` is False so `deactivating`; `v_norm <= 0` is True so `shortening`. This is what a designer gets by calling a finite-difference linearizer on the nonlinear model with one-sided differences, or by evaluating the model's branch logic at the equilibrium. |
| `blended` | central finite differences of the nonlinear model straight through both kinks -- the silent failure mode `linearize.py`'s docstring warns about. Each kinked Jacobian entry becomes the *average* of its two one-sided values (verified in `tests/test_mpc.py`). |
| `tracking` | branch chosen each step from the current operating point: `activating` if `u_{k-1} > a_k`, `shortening` if `thetadot_k > 0`; ties (within `1e-6`) broken by the direction the reference is asking for (`theta_ref > theta_k` means flex, so activate/shorten). This is the "track which side of each switch the operating point is on" the abstract recommends. Same linearization *point* as the other two; only the branch differs. |
| `current_state` | full successive linearization at the measured `(x_k, u_{k-1})` instead of the reference equilibrium, with the same branch rule as `tracking` for the at-rest ties. Included to separate "wrong branch" from "wrong linearization point": off the manifolds this variant's branch is unambiguous, so if it matches `tracking` the branch, not the point, is what matters. |
| `trajectory` (added after the first run) | per-horizon-step branch, read off the previous QP's planned `(x, u)` sequence (shifted one step), then re-solved until the assignment read off the new plan agrees with the one it was built with (successive piecewise linearization, at most 4 iterations; if it cycles, fall back to the `tracking` branch for all steps). The principled version of "track the switch", over the whole horizon rather than at the current instant. |

All four variants run the same task on the same true nonlinear plant
(`joint.state_derivative` integrated with `solve_ivp` between control
updates), from the same initial condition, with the same weights and
horizon.

## 3. Task and metrics

Reference (degrees, `theta` about the Phase 2/4 equilibrium):
hold 65 for 0.3 s; step to 68 at 0.3 s; step to 62 at 1.3 s; hold to
2.3 s; then a 0.5 Hz sinusoid 65 +/- 3 to 4.3 s. Steps of 3 deg are
small enough to stay near the small-signal regime Phase 4 isolated
(a 0.01 excitation step moved the joint ~1.5 deg in 200 ms), large
enough to be a realistic posture adjustment.

Metrics per variant:

- **closed-loop tracking**: RMS and peak `|theta - theta_ref|` (deg),
  per segment (step up, step down, sinusoid) and overall; settling time
  to within 0.2 deg after each step; overshoot.
- **internal model error**: one-step-ahead residual
  `|theta_pred(k+1|k) - theta_true(k+1)|`, where `theta_pred` is the
  controller's own prediction from `x_k` under the `u_k` it actually
  applied. This is the closed-loop analogue of the Phase 4/6a quantity,
  and the direct measure of "how wrong is the model the QP is optimizing
  over."
- **effort**: RMS `u`, number of saturated steps, QP failures.

## 4. What the outcomes would mean

- If `fixed_naive`/`blended` track about as well as `tracking`: feedback
  washes the ambiguity out at this rate and these gains; the finding is
  real for prediction/estimation but not a closed-loop design hazard at
  `T_s=10 ms`. Report as such; it narrows the abstract's practical claim.
- If `tracking` tracks measurably better (lower RMS, less overshoot, or
  smaller residuals) with identical weights: the ambiguity survives
  feedback and the abstract's recommendation is justified in closed loop.
- Either way, `current_state` vs `tracking` says whether the fix is
  "choose the branch" (cheap: a sign test) or "relinearize every step"
  (more expensive, and still needs the sign test at rest).

## 5. Results

Run 2026-08-24 on the remote (`python3 scripts/phase6b_closed_loop.py`;
`tests/test_mpc.py` 14/14). Plots: `results/phase6b_closed_loop.png`
(full task), `results/phase6b_step_zoom.png` (+3 deg step); numbers in
`results/phase6b_metrics.json`. Weights `Q=diag(1000, 1, 0)`, `R=200`,
`N=20`, `T_s=10 ms`, identical for every variant.

### 5.1 Two things the first runs taught, before any comparison was valid

1. **Chattering of the sign-test rule.** With a strict sign test the
   `tracking` branch flips every step at rest (`u-a` and `thetadot`
   hover at solver precision around zero on the manifold), producing a
   +/-0.01 sawtooth in `u`. Fixed with deadbands (`|u-a| > 2e-3`,
   `|thetadot| > 0.02 rad/s`, else the reference direction decides,
   else keep the previous branch): sawtooth amplitude 1e-15 afterwards.
   Worth recording because it is the first practical cost of "track the
   switch": a sign test alone is not implementable.
2. **A spurious steady-state offset, my bug, not the plant's.** The
   first runs showed `fixed_naive` settling 0.26 deg above a 68 deg
   reference (fully at rest, `thetadot ~ 1e-14`), the other variants
   0.004-0.016 deg. This looked like a branch effect and was not: the
   terminal DARE cost was centred on `[theta_ref, 0, 0]` instead of the
   full equilibrium `[theta_ref, 0, a_ref]`, and since `P` couples `a`
   with `theta` (differently per branch) that biased the fixed point.
   With the reference corrected, every one of the 4 fixed branches held
   at 68 deg settles to an offset below 1e-13 deg. Recorded because a
   less careful pass would have reported "wrong branch causes 0.26 deg
   steady-state error" as a finding. (Also noted along the way: the 68
   deg reference equilibrium is open-loop *unstable*, eigenvalue +0.2 to
   +0.39 /s depending on branch, which is why it was so sensitive.)

### 5.2 The comparison (deg)

| variant | RMS track, all | sinusoid RMS | step-up settle (s) / overshoot | step-down settle (s) / overshoot | RMS 1-step residual | peak residual |
|---|---|---|---|---|---|---|
| `fixed_naive`   | 0.496 | 0.214 | 0.10 / 0.07 | 0.24 / 0.00 | 0.0012 | 0.0033 |
| `blended`       | 0.421 | 0.194 | **0.24 / 0.62** | 0.27 / 0.08 | 0.0034 | 0.0240 |
| `tracking`      | 0.460 | **0.155** | 0.10 / 0.08 | 0.14 / **0.31** | 0.0046 | 0.0351 |
| `current_state` | 0.483 | **0.156** | 0.10 / 0.08 | 0.16 / 0.01 | **0.0007** | 0.0075 |
| `trajectory`    | 0.459 | 0.157 | 0.10 / 0.14 | 0.16 / 0.37 | 0.0046 | 0.0351 |

(Overall RMS is dominated by the unavoidable error during the two 3 deg
step transitions and is nearly the same for all; the sinusoid column and
the transient columns are where the variants differ. Steady-state offset
after both steps is below 0.003 deg for every variant. No QP failures,
no saturation.)

- **Feedback washes most of it out.** The 70-88% open-loop one-step
  prediction errors of Phase 4/6a become closed-loop differences of at
  most 0.06 deg RMS on a 3 deg sinusoid and 0.55 deg of overshoot on a
  3 deg step. No variant is unstable, none has steady-state error. On
  the abstract's practical claim this is a narrowing: the ambiguity is
  a prediction/estimation problem first, a control-performance problem
  second, and not a stability problem at `T_s=10 ms` with these gains.
- **Where it does not wash out.** The `blended` (central-difference)
  model, the failure mode the abstract warns about most, is the worst
  transient by a wide margin: 8x the step-up overshoot of every other
  variant (0.62 vs 0.07-0.08 deg) and 2.4x slower settling, because its
  activation gain is the average of a fast and a slow branch and
  matches neither. The naive fixed branch has 38% higher sinusoid RMS
  than the state-matched rules (0.214 vs 0.155 deg).
- **Matched-to-current-state is not uniformly best.** `tracking`
  overshoots the step-*down* by 0.31 deg where `fixed_naive` overshoots
  0.00: on the way down the plan has to *re-activate* to brake, i.e.
  cross the activation switch mid-horizon, and a single branch held over
  200 ms is wrong for the second half of the plan whichever branch it
  is. The naive branch happens to be the braking-phase branch. This is
  the mechanism behind Phase 4's horizon-dependence finding, seen from
  the controller's side.
- **Per-step branch assignment did not help.** `trajectory` ended within
  0.002 deg RMS of `tracking` on every metric, with a worse step-down
  overshoot (0.37 deg); the fixed-point iteration cycled on 135 of 430
  control steps (mean 2.1 solves/step, max 4) and fell back. The QP
  solution is discontinuous in the branch assignment, so plain
  iteration does not converge; the right tool for the over-the-horizon
  version of this problem is a hybrid/PWA MPC (mixed-integer, or an
  explicit enumeration of the 4^N assignments with pruning), not a
  successive-linearization loop. Not built here.
- **Branch vs. point.** `current_state` (relinearize at the measured
  state, branch by sign test) has the smallest internal model error
  (RMS residual 0.0007 deg, 7x below `tracking`) but the same closed-loop
  tracking as `tracking`. So the branch choice, not the linearization
  point, is what moves the closed-loop numbers; and the sign test at
  rest is needed either way (`current_state` sits on both manifolds
  whenever the joint is at rest at its setpoint, i.e. most of the time).

### 5.3 Answer

The branch ambiguity survives feedback as a **transient-quality**
effect, not a stability or steady-state one, at `T_s=10 ms`: choosing
the branch from the measured operating point (a sign test with a
deadband, essentially free) buys ~40% lower periodic tracking error over
a naive fixed choice, and avoiding the blended central-difference model
avoids an 8x overshoot penalty. Holding *any* single branch over a 200
ms horizon is the remaining source of error, and fixing that properly
needs hybrid MPC. The sentence for the abstract's future-work line is
therefore: "in closed loop at 10 ms the ambiguity degrades transient
tracking (up to 8x overshoot for a blended model, ~40% RMS for a naive
fixed branch) but not stability; a sign-test branch rule recovers most
of it."
