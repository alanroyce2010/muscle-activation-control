# Phase 7: Closing the Project's Stated Limitations

Started 2026-09-11 at the user's request ("go through the current limits
of the project and implement them"). Written before the code, per the
project's derive-before-code rule; each item states what would count as
the limitation being closed and what the closed forms predict, so the
numerical results can confirm or refute them. Results are appended per
section as they come in. No em dashes by convention.

## 0. Inventory and triage

The limitations as the manuscripts and `docs/PROJECT_WALKTHROUGH.md`
section 12 state them, sorted by whether simulation can close them:

| limitation | closable in simulation? | section |
|---|---|---|
| moment arms digitized by eye from published figures | yes: pixel-level re-digitization, checked against the source's own tables | 1 |
| triceps moment arm 27% above Holzbaur (the "one soft number") | yes: second anchor, measured-torque check, sensitivity | 2 |
| one muscle; Thelen young-adult constants only | yes: other elbow flexors, other constant sets | 3 |
| no torque validation (no brachialis-only measurement exists) | partly: a flexor-group model can be compared with measured group moments | 3 |
| rigid tendon | yes: elastic-tendon variant (Thelen Eq. 5 tendon) | 4 |
| exact state measurement, no noise, no delay, one control rate, no parameter mismatch | yes: closed-loop realism sweep | 5 |
| horizon-level branch assignment left to "hybrid MPC, not built" | yes: mode-sequence enumeration MPC | 6 |
| no explanation of *why* feedback washes the ambiguity out | yes: nu-gap vs generalized stability margin | 7 |
| hardware validation; EMG test of the co-contraction prediction; reflexes, short-range stiffness, history dependence; systematic literature review | no: experimental or a change of model class; stays future work | 8 |

## 1. Digitization audit (done 2026-09-11, before any other item)

**Why first.** Every number in the project passes through the brachialis
moment arm `r(theta)`. The original table
(`single_joint_dynamics.md` section 3) was read by eye at low resolution
with a stated error of +/-0.15-0.2 cm per point.

**Method.** Murray et al. (1995) Fig. 4 rendered at 600 dpi; axes
calibrated from tick marks by least squares (13.656 px per degree,
103.25 px per cm, model panel); every curve tracked at 1-degree steps by
nearest-neighbour continuity on dark-pixel run centres, with dash gaps
carried and label text boxes excluded; tracks verified by overlay.

**Objective check against the source.** Murray's Table 2 gives, for the
model panel, the percent difference (max - min)/max over 25-120 degrees.
The pixel tracks reproduce all four:

| muscle | Murray Table 2 (model) | pixel digitization | original eyeball table |
|---|---|---|---|
| BRA | 48% | 48.0% | about 63% |
| BIC | 43% | 41.9% | (not digitized) |
| BRD | 55% | 55.9% | (not digitized) |
| TRI | 25% | 22.0% | about 25% |

**Finding.** The original brachialis table reads low below 90 degrees by
more than its stated error: 0.89 vs 0.6 cm at 0 deg, 1.38 vs 1.0 at 20,
1.85 vs 1.4 at 40 (0.45 cm, twice the stated tolerance), 2.23 vs 1.9 at
60; it agrees above 90 degrees. The original triceps table agrees with
the pixel track within 0.13 cm everywhere. **Decision: replace both
tables with the pixel tracks** (5-degree grid) and rerun every phase.
This is a correction to the data, recorded as a revision note in
`single_joint_dynamics.md` section 3.

**Predicted consequences.** At 65 deg the moment arm rises from 2.00 to
2.32 cm (+16%), so the gravity-only equilibrium activation falls from
0.133 to about 0.115 and `rho_a` rises from 6.8 to about 7.4 (it depends
on `a*` only). The range-of-motion average rises from 1.94 to 2.18 cm,
now 21% above Holzbaur's 1.8 cm (which is a different model with
different wrapping; Murray's own model, not Holzbaur's, is what the
tables transcribe). Structural results (switch placement, exact
continuity, disjoint support, `rho_v = 2`) cannot change. Quantitative
results should move by a few percentage points; the headline claims
(13-88%, never below 20%, 8x overshoot) are predicted to survive in
shape, with updated numbers.

## 2. Triceps moment arm: the one soft number

Two independent anchors disagree with the Murray model-panel triceps
curve (ROM mean 2.60 cm): Holzbaur (2005) Table 1 `ma_avg = 2.1 cm` and
Holzbaur Fig. 3A, whose model point at 90 deg reads about 2.0 cm (the
Murray curve gives 2.31 cm there). Two readings of "the truth":

- **Holzbaur-consistent**: scale the Murray shape by `s = 2.1/2.60 =
  0.81` (at 90 deg: 1.87 cm, within 0.15 cm of Holzbaur's Fig. 3A point).
- **Murray model as digitized**: `s = 1`.

**A measured-data tiebreaker.** Holzbaur Fig. 5A reproduces measured
maximum isometric elbow *extension* moments (Amis et al.; Buchanan et
al.): about 35-52 N m, peaking near 50-52 N m at 80-100 deg, with
Holzbaur's own model peaking at about 43 N m (17% low, stated in their
text). A full triceps (long + lateral + medial, F_0 = 2047 N, same elbow
moment arm for all heads per Holzbaur) at full activation predicts a
peak of roughly `2047 N x r x f_L`: about 53 N m at `s = 1` and 43 N m at
`s = 0.81`. **Prediction:** the unscaled curve matches measured strength,
the scaled curve matches Holzbaur's model; the moment-arm uncertainty is
therefore bracketed by `s` in [0.81, 1] and the co-contraction balance
point should be reported as a range across it.

## 3. Generalization

### 3.1 Other elbow flexors

Same single-muscle plant with brachialis replaced by biceps (long + short
heads lumped, shoulder fixed neutral; Holzbaur Table 1: F_0 = 1059.9 N,
l_0 = 11.6/13.2 cm, both heads share one elbow moment arm) or
brachioradialis (F_0 = 261.3 N, l_0 = 17.3 cm), moment arms from the
pixel tracks of section 1. Predictions: the structure is identical
(model form); `a*` at 65 deg scales with `m g l_c sin(theta) / (F_0 r)`:
biceps about 0.08, brachioradialis about 0.26, brachialis about 0.115,
giving `rho_a` about 8.9, 4.5 and 7.4. Wrong-activation error should rank
biceps > brachialis > brachioradialis; `rho_v = 2` for all.

### 3.2 A torque validation at the group level

No brachialis-only torque measurement exists, but maximum isometric
*group* flexion moments do (Holzbaur Fig. 5A, Amis and Buchanan: about
40-79 N m, peaking at 76-79 N m near 60-100 deg). The three flexors with
published parameters and Murray moment arms (BRA + BIC + BRD) summed at
full activation give a model group curve to compare. This is the
torque-level validation the project lacked; it tests force capacity and
moment-arm magnitude together. The fibre-length reference convention
(optimal length at mid-range) is a modelling choice that shapes the
curve; the comparison is on magnitude and peak location.

### 3.3 Other Thelen constant sets

| set | tau_act / tau_deact | A_f | F_len | V_max | gamma | k_PE | eps0_M |
|---|---|---|---|---|---|---|---|
| paper, young (current) | 15 / 50 ms | 0.25 | 1.4 | 10 | 0.45 | 5 | 0.6 |
| paper, old adult | 15 / 60 ms | 0.25 | 1.8 | 8 | 0.45 | 5 | 0.5 |
| OpenSim Thelen2003Muscle | 10 / 40 ms | 0.3 | 1.8 | 10 | 0.5 | 4 | 0.6 |

Predictions: `rho_v = (2 + 2/A_f)/(1 + 1/A_f) = 2` exactly for every set
(`A_f` and `F_len` both cancel); `rho_a` base ratio 3.33 (young) vs 4.0
(old and OpenSim), so the activation kink is *larger* in both
alternatives. Requires `muscle.py` defaults to bind late (module
constants read at call time) so a constant set can be swapped in one
place; the young-adult default must reproduce every existing test.

## 4. Elastic tendon

Thelen (2003) Eq. 5 tendon (`eps0_T = 0.04`, `k_toe = 3`,
`F_toe = 0.33`, `eps_toe = 0.609 eps0_T`, `k_lin = 1.712/eps0_T`),
brachialis tendon slack length 5.35 cm, pennation zero. State
`[theta, thetadot, l_M, a]`. With an elastic tendon the natural direction
of Thelen's force-velocity law is the published one: given
`(theta, l_M)` the tendon force `F_T` is known, the contractile force is
`F_ce = F_T - F_PE(l_M)`, and Eq. 6-7 give the fibre velocity directly:

```
dl_M/dt = (0.25 + 0.75a) V_max (F_ce - a f_L) / b,
b = a f_L + F_ce/A_f                      if F_ce <= a f_L  (shortening)
b = (2 + 2/A_f)(a f_L F_len - F_ce)/(F_len - 1)   otherwise (lengthening)
```

**Predictions.** (i) Any equilibrium needs `dl_M/dt = 0`, i.e.
`F_ce = a f_L`: exactly the boundary of Eq. 7's two branches, so the
velocity switch persists at every equilibrium, now as a *fibre*
switch with switching function `y = F_ce(theta, l_M) - a f_L(l_M)`,
whose `c` vector has several non-zero entries. (ii) Both branches give
`dl_M/dt = 0` on the switch and the slope `d(dl_M/dt)/dy` differs by
exactly 2 (b doubles on the lengthening side), so the continuity
condition holds with `e` supported only on the `l_M` row. (iii) The
activation switch's `e` stays on the `a` row, so supports remain
disjoint and the direct sum persists. (iv) The tendon adds a fast
fibre-length mode whose rate is set by tendon stiffness over the
force-velocity slope; since that slope doubles across the switch, the
fast mode's rate changes by a factor of about 2 across the velocity
branch. (v) The continuity check needs the general solver
`e = (A1 - A2) v / (c^T v)` for any `v` with `c^T v != 0`, then the
residual `||A1 - A2 - e c^T||` (the existing helper assumes one non-zero
`c` entry).

## 5. Closed-loop realism

Rerun the Phase 6b task with the three rules that matter (naive fixed,
blended, sign test), varying one factor at a time from the Phase 6b
baseline:

- **Control rate**: `T_s` in {5, 10, 20, 40} ms, horizon held at 200 ms.
- **Measurement noise and estimated activation**: angle noise sigma in
  {0, 0.05, 0.2} deg, velocity from a filtered difference of the noisy
  angle, activation *not measured* but estimated by integrating the
  model's activation ODE from the commanded input; 5 noise seeds.
- **Plant-model mismatch**: plant activation time constants off by
  +/-25% (and the old-adult set) while the controller keeps the young
  set; this makes the activation estimate wrong in a structured way.
- **Transport delay**: 0, 20, 50 ms on the input, not modelled by the
  controller.

Predictions: the blended penalty is a model bias, not a noise effect, so
it persists under every factor; the sign test's advantage over the naive
branch should shrink with noise (branch misclassification near rest) but
its deadband should keep it from chattering; delay degrades all rules
and may destabilize at 50 ms; no factor should reverse the ranking at
moderate levels.

## 6. Hybrid MPC by mode-sequence enumeration

The residual Phase 6b error is holding one branch over the horizon. A
proper fix assigns a branch per horizon step. With no MIQP solver on the
remote, use structured enumeration: each switching signal gets at most
one switch per horizon, at a time `k` from a grid
{0, 1, 2, 3, 5, 7, 10, 14, N} (N meaning "never"), giving 81 candidate
mode sequences (current-side branch before `k`, opposite after). Solve
the QP for each; keep candidates whose own plan is self-consistent (the
sign sequence read off the plan matches the assumed one outside the
deadband); apply the lowest-cost consistent plan (fall back to the sign
test if none is consistent). Prediction: it anticipates the braking
re-activation and removes the sign test's step-down overshoot
(0.31 deg in Phase 6b) while keeping its step-up and sinusoid
performance.

> **Revision note (2026-09-11):** the implementation (`mpc.py`,
> `_compute_hybrid`) does not filter to consistent candidates with a
> sign-test fallback as designed above; it ranks all 81 candidates by the
> number of self-consistency violations, then by cost, and applies the
> best. That is why some applied plans still carry violations (46 of 430
> steps, section 9.7). The papers state the implemented rule.

## 7. Why feedback washes it out: nu-gap vs stability margin

For discrete-time plants `P_i(z)` (branch models, input `u`, full state
output) and the unconstrained LQR gain `K` that the MPC reduces to near
the equilibrium, Vinnicombe's result says `K` stabilizes every plant
within nu-gap `delta_nu(P_1, P_2) < b(P_1, K)` of the design plant,
where `b(P, K) = || [P; I] (I - K P)^{-1} [-K, I] ||_inf^{-1}`. Compute
`delta_nu` between each pair of branch models (and blended vs each) and
`b` for each design branch. Prediction: every `delta_nu` is below the
corresponding `b`, which explains the Phase 6b observation that no rule
destabilizes; the blended model sits between the branches in `delta_nu`
but still carries the largest transient error, because the nu-gap bounds
stability, not performance. Caveat to state: this is a statement about
the linear branch models, not a proof for the piecewise-affine plant.

## 8. What stays future work

Hardware validation of the sign-test rule; an EMG test of the
co-contraction prediction; reflexes, short-range stiffness and history
dependence (a change of model class, not a parameter); a systematic
database literature review; a certificate of closed-loop stability for
the sampled-data piecewise-affine loop.

## 9. Results

Run 2026-09-11 on the remote; suite 81 tests. Numbers come from
`scripts/headline_numbers.py`, `scripts/phase7_*.py` and the rerun
Phase 1-6 scripts (JSON in `results/`).

### 9.1 Digitization audit, and what it uncovered

Confirmed as predicted: all five pixel tracks reproduce Murray's Table 2
(BRA 48.0 vs 48, BIC 41.9 vs 43, BRD 55.9 vs 55, TRI 23.4 vs 25, male
BRA 58.9 vs 58%). Brachialis peak 2.88 cm at 110 deg, range-of-motion
average 2.18 cm.

**Two further problems surfaced while rerunning, both fixed.**

1. *The old 65 deg equilibrium was open-loop unstable.* Rechecking the
   eigenvalues of the ORIGINAL model (by-eye table) at 65 deg with
   one-sided theta derivatives gives a slowest eigenvalue of +0.17/s
   (central difference; +0.02 and +0.31/s one-sided). The earlier record
   ("all real, negative, between -11.6 and -95.3/s", CLAUDE.md 4.1 and
   both manuscripts) omitted this slow mode: -11.6 was the mechanical
   fast mode. The old fold of the gravity-only equilibrium family (peak
   of `a*(theta)`) sat at 53.4 deg, so 65 deg was on the unstable side;
   consistent with Phase 6b's observation that 68 deg was unstable. None
   of the branch results depended on stability (the prediction protocol
   compares linear models of the same equilibrium), but the stated
   property was wrong. `scripts/_diag_eigs.py` reproduces both models.
2. *65 deg is a triple coincidence.* By the reference-length convention
   the fibre is exactly at optimal length there, where the clipped
   passive curve `F_PE` has its own kink (one-sided stiffness entries
   -5.04 vs -1.96 /s^2 in the corrected model, so the central difference
   used for the theta column blended them: the very failure mode this
   project warns about); it was a knot of the 5-degree table
   (piecewise-linear interpolation puts an artificial stiffness kink at
   every knot); and with the corrected table it lies 5 deg from the new
   fold (69.9 deg).

**Fixes.** Moment arms now use a monotone C1 (PCHIP) interpolant through
the pixel tables (no knot kinks; data points reproduced exactly). The
representative posture moves to **60 deg**: fibre at l_bar = 1.023
(f_L = 0.9988, essentially the peak; passive curve smooth), off every
knot, 10 deg on the stable side of the fold. The closed-loop task is
recentred (60 -> 63 -> 57 deg, sinusoid 60 +/- 3 deg, all on the stable
side); the range-of-motion sweep grid moves to 12.5, 17.5, ..., 122.5 deg
(off the 65 deg passive kink).

**Corrected headline numbers (60 deg).** `a* = 0.1126`, `rho_a = 7.45`,
`e_act = [0, 0, 86.29]`, `e_vel = [0, -12.63, 0]`, residuals exactly 0;
all four branch Jacobians stable, eigenvalues -99.7 to -0.3/s (the slow
mode is the near-fold mechanical mode); `max|A_d| = 7.72` over T_s 1-100
ms (non-transfer unchanged). Activating-step prediction table (% of true
response): 50 ms 1.1 / 71.1 / 13.8 / 74.9 (matched / wrong act. / wrong
vel. / both), 100 ms 1.3 / 57.3 / 24.1 / 66.9, 200 ms 1.2 / 37.6 / 35.1 /
58.4, 500 ms 0.9 / 13.7 / 42.9 / 50.7: **wrong-branch range 13.7-74.9%**
(was quoted as 13-88%). Deactivating steps: wrong activation 241% at
50 ms. Peak isometric torque 25.8 N m at 94 deg; gravity-only `a*` over
the range of motion 0.028-0.115.

### 9.2 Triceps moment arm

Pixel TRI mean 2.59 cm; Holzbaur-consistent scale `s = 2.1/2.59 = 0.810`
gives 1.87 cm at 90 deg against Holzbaur Fig. 3A's model point of about
2.0 cm, so the two anchors agree with each other. **The measured-data
tiebreaker does not discriminate** (prediction refuted): full triceps
(three heads) at full activation peaks at 55.6 N m (s = 1) or 45.4 N m
(s = 0.81) against measured peaks of 51-52 N m, but the RMS difference
to Amis is 8.9 vs 9.3 N m (a tie) and to Buchanan 19.1 vs 14.6 N m
(favouring s = 0.81); the peak angle is off in both (model 16-55 deg,
measured 80-100 deg), which points at the uniform optimal-length
reference convention rather than the moment arm.

**Consequence for co-contraction (prediction refuted, claim weakened).**
With s = 1 the force-velocity kink cancels at `a_t* = 0.264`
(`a_b* = 0.514`; wrong-velocity error 0.8% at 50 ms there). With
s = 0.90 or s = 0.81 **no balance exists within feasible co-contraction**:
the ratio falls from 2 to a minimum of 1.008 (s = 0.9) or 1.068
(s = 0.81) at the highest feasible triceps activation. The robust
statement is that co-contraction shrinks the velocity kink from 2 to
within 7% of 1; exact cancellation depends on the triceps moment arm,
which the literature does not pin down to the required accuracy.

### 9.3 Generalization

*Other flexors (predictions confirmed).* At 60 deg:

| plant | a* | rho_a | rho_v | residual | 50 ms wrong act. / vel. / both | activating range |
|---|---|---|---|---|---|---|
| brachialis | 0.113 | 7.45 | 2.0000 | 5e-10 | 71.1 / 13.8 / 74.9% | 13.7-74.9% |
| biceps (two heads) | 0.082 | 8.57 | 2.0000 | 6e-10 | 73.3 / 13.0 / 76.7% | 12.7-76.7% |
| brachioradialis | 0.261 | 4.20 | 2.0000 | 2e-10 | 59.6 / 9.2 / 63.2% | 9.2-63.2% |

Same structure in all three (disjoint e supports, exact rho_v = 2); the
wrong-activation ranking biceps > brachialis > brachioradialis follows
`rho_a(a*)` as predicted.

*Group torque validation (new).* Brachialis + biceps + brachioradialis at
full activation peak at **73.6 N m at 94 deg**; Buchanan et al.'s measured
flexion moments peak at 79 N m at 100 deg and the **RMS difference over
their range is 3.6 N m**. Amis et al.'s curve (peaking at 76 N m near 60
deg) differs in shape (RMS 17.7 N m); Holzbaur's own model differs from
Amis by a similar amount. This is the torque-level validation the
project previously lacked.

*Constant sets (predictions confirmed).* rho_v = 2.000000 for all three;
rho_a 7.45 (young), 8.95 (old adult), 9.01 (OpenSim): the activation
kink is larger in both alternatives. Activating-step range 13.7-74.9%
(young), 16.4-79.1% (old), 10.3-73.0% (OpenSim); deactivating wrong
activation at 50 ms 241%, 296%, 222%.

### 9.4 Elastic tendon (predictions confirmed)

Equilibrium at 60 deg: `l_bar* = 1.013`, `a* = 0.1132` (rigid 0.1126),
switching function 4e-15. Activation switch e on the activation row
only (residual 5e-12); fibre switch with `c = grad y = [-6.74, 0, -25.97,
-1]` (several non-zero entries, as predicted) and e on the fibre-length
row only (residual 8e-8, finite-difference level); **the thetaddot row
carries no branch dependence at all** in the elastic model; fibre-row
slope ratio exactly 2.0000. The tendon adds a fast fibre mode (-139.9/s
on the shortening branch) that merges with the mechanical mode into an
oscillatory pair (real part -38.3/s) on the lengthening branch: the
branch changes the mode's character, not just its rate (prediction iv
only partly right). Conditioning `max|A_d| = 67` (rigid 7.7; pneumatic
900,000). Linearize-and-hold: activating range **14.3-80.0%** (rigid
13.7-74.9%), matched branch at most 1.2%; deactivating wrong activation
269% at 50 ms. The rigid-tendon assumption does not drive the result.

### 9.5 Why feedback washes it out

*Nu-gap (prediction refuted).* LQR stability margins `b(P_i, K_i)` are
0.060-0.081; nu-gaps between activation branches are 0.71-0.73, between
velocity branches 0.06-0.12, blended vs branches 0.24-0.59. Only 1 of 20
design/plant pairs is certified by `delta_nu < b`, yet all 20 are
directly stable (spectral radius 0.82-0.98). The nu-gap treats the branch
difference as unstructured and is far too conservative here.

*Common quadratic Lyapunov function (partial certificate).* Under
*arbitrary* switching of the plant among its four branch modes, a
quadratic Lyapunov function certifies the unconstrained local loop for 4
of the 5 fixed designs (contraction factors 0.928-0.999; blended 0.987);
the activating/lengthening design and the sign-test rule under arbitrary
misclassification (16 modes) are not certified. **The sign-test rule
with correct classification (its four matched modes, switching
arbitrarily) is certified, contraction factor 0.9215.** So the local
loop is provably stable under arbitrary mode switching for the sign-test
rule as long as it reads the branch correctly, and for 4 of 5 fixed
designs; what stays empirical is robustness to misclassification inside
the deadband, which the realism sweep (section 9.6) probes directly.
Caveat: these are statements about the unconstrained local loop (the
MPC's input bounds inactive near the equilibrium) and the discrete
branch models; switching within a sample is not modelled.

### 9.6 Closed-loop realism (section 5 predictions: mostly confirmed)

`scripts/phase7_realism.py`, 64 runs, `results/phase7_realism.json`.
Step-up overshoot and sinusoid RMS (deg), naive (N) / blended (B) / sign
test (S):

| condition | overshoot N / B / S | sinusoid RMS N / B / S |
|---|---|---|
| baseline, 10 ms | 0.07 / 0.65 / 0.08 | 0.20 / 0.20 / 0.15 |
| T_s = 5 ms | 0.07 / 0.61 / 0.08 | 0.20 / 0.19 / 0.16 |
| T_s = 40 ms | 0.08 / 0.99 / 0.14 | 0.22 / 0.27 / 0.13 |
| noise 0.05 deg (5 seeds) | 0.18 / 0.83 / 0.14 | 0.27 / 0.26 / 0.18 |
| noise 0.2 deg (5 seeds) | 0.75 / 1.12 / 0.55 | 0.67 / 0.53 / 0.47 |
| plant tau x0.75 | 0.01 / 0.45 / 0.00 | 0.20 / 0.16 / 0.17 |
| plant tau x1.25 | 0.17 / 0.84 / 0.19 | 0.21 / 0.23 / 0.16 |
| old-adult plant | 0.15 / 0.80 / 0.18 | 0.20 / 0.22 / 0.15 |
| delay 20 ms | 0.11 / 1.19 / 0.27 | 0.18 / 0.30 / 0.14 |
| delay 50 ms | 2.98 / diverged / diverged | 11.6 / diverged / diverged |

(Noise rows use velocity from a first-order filtered difference of the
noisy angle, 20 ms time constant, and activation estimated by integrating
the nominal activation ODE from the commands; mismatch rows estimate
activation the same way; all other rows measure the state.)

- **Confirmed:** the blended model overshoots most in every condition
  (9-12x the naive branch across rates); the sign test tracks the
  sinusoid best in all conditions but one (plant tau x0.75, where the
  blended model beats it by 0.004 deg). Under noise the sign test's
  branch choice chatters (147 branch changes in 4.3 s at 0.05 deg, 263 at
  0.2 deg; 16 without noise) yet it keeps the lowest sinusoid error: the
  deadband keeps the chatter harmless.
- **Refuted in part:** an unmodelled transport delay overturns the
  ranking. At 20 ms the sign test overshoots the downward step by
  1.43 deg (naive 0.00), because it reads the branch from a state that no
  longer describes the plant; at 50 ms the sign-test and blended loops
  diverge and the naive loop degrades to 11.6 deg RMS. Delay must be in
  the prediction model; that is ordinary MPC practice, not specific to the
  branch problem, but the sign test is the rule most exposed to it.

> **Revision note (2026-09-11, after the red-team review of the ICRA
> draft):** the table above omits the downward step, and the first bullet
> overstates the sign test. The review asked for step-down overshoot and
> for relinearization at the measured state ("current_state", the fourth
> rule) under the same conditions: `scripts/phase7_realism.py --only
> current_state`, 20 runs, `results/phase7_realism_current_state.json`.
> Full table, overshoot up / down and sinusoid RMS (deg), naive (N),
> blended (B), sign test (S), relinearize (R):
>
> | condition | up N / B / S / R | down N / B / S / R | sinusoid N / B / S / R |
> |---|---|---|---|
> | baseline, 10 ms | 0.07 / 0.65 / 0.08 / 0.07 | 0.00 / 0.15 / 0.37 / 0.02 | 0.200 / 0.196 / 0.153 / 0.149 |
> | T_s = 5 ms | 0.07 / 0.61 / 0.08 / 0.07 | 0.00 / 0.11 / 0.33 / 0.02 | 0.198 / 0.189 / 0.161 / 0.152 |
> | T_s = 20 ms | 0.07 / 0.76 / 0.09 / 0.07 | 0.00 / 0.20 / 0.56 / 0.00 | 0.204 / 0.213 / 0.143 / 0.138 |
> | T_s = 40 ms | 0.08 / 0.98 / 0.14 / 0.14 | 0.00 / 0.41 / 0.58 / 0.00 | 0.215 / 0.274 / 0.133 / 0.119 |
> | estimated a, no noise | 0.00 / 0.74 / 0.00 / 0.00 | 0.00 / 0.24 / 0.17 / 0.00 | 0.188 / 0.208 / 0.139 / 0.132 |
> | noise 0.05 deg | 0.18 / 0.83 / 0.14 / 0.15 | 0.00 / 0.25 / 0.20 / 0.00 | 0.273 / 0.256 / 0.175 / 0.195 |
> | noise 0.2 deg | 0.75 / 1.11 / 0.55 / 0.63 | 0.00 / 0.27 / 0.15 / 0.00 | 0.667 / 0.525 / 0.470 / 0.586 |
> | plant tau x0.75 | 0.01 / 0.45 / 0.00 / 0.00 | 0.00 / 0.11 / 0.33 / 0.01 | 0.200 / 0.163 / 0.167 / 0.151 |
> | plant tau x1.25 | 0.17 / 0.84 / 0.19 / 0.19 | 0.00 / 0.19 / 0.32 / 0.03 | 0.205 / 0.234 / 0.155 / 0.150 |
> | old-adult plant | 0.15 / 0.80 / 0.18 / 0.18 | 0.00 / 0.13 / 0.31 / 0.01 | 0.203 / 0.224 / 0.152 / 0.147 |
> | delay 20 ms | 0.11 / 1.19 / 0.27 / 0.26 | 0.00 / 1.02 / 1.42 / 0.02 | 0.184 / 0.298 / 0.142 / 0.132 |
> | delay 50 ms | 2.97 / div / div / 1.88 | 3.83 / div / div / 2.01 | 11.6 / div / div / 2.32 |
>
> What this changes:
> - The naive branch never overshoots downward (0.00 everywhere before
>   the 50 ms delay); the sign test always does (0.15-0.58 deg).
> - Relinearization uses the same branch rule as the sign test and
>   differs only in the linearization point, yet overshoots downward by
>   at most 0.03 deg and tracks the sinusoid best in every noise-free
>   condition. So the old Phase 6b line "the branch, not the linearization
>   point, moves the closed-loop numbers" is wrong for the downward step:
>   the point matters as much as the branch there.
> - Under angle noise the sign test tracks better than relinearization
>   (0.175 vs 0.195 at 0.05 deg; 0.470 vs 0.586 at 0.2 deg), plausibly
>   because relinearization moves its linearization point with the noise.
>   Branch changes reach 269 (sign test) and 276 (relinearize) in 4.3 s.
> - Relinearization is the only branch-reading rule that survives the
>   50 ms delay (about 2 deg), and at 20 ms it does not overshoot
>   downward (0.02 vs the sign test's 1.42).
> - Practical recommendation, as now stated in the papers: relinearize
>   at the measured state each step with the branch read from it and a
>   deadbanded sign test at rest; model transport delay.

### 9.7 Hybrid mode-sequence MPC (section 6 prediction: partly confirmed)

81 candidate mode sequences per step (switch grid {0, 1, 2, 3, 5, 7, 10,
14, 20}). Baseline: downward overshoot **0.37 -> 0.24 deg** vs the sign
test, step-up overshoot 0.089 vs 0.081, sinusoid RMS 0.159 vs 0.153; 46 of
430 applied plans still carried self-consistency violations. With
0.05 deg noise (one seed): downward overshoot **0.35 -> 0.14 deg**, sinusoid
RMS 0.134 vs 0.174. It anticipates the braking re-activation as
predicted, but does not reach the naive branch's zero downward overshoot,
and it costs 81 QPs per step (158 s for a 4.3 s run in cvxpy/Python,
about 37x slower than real time on this machine, but the single-QP sign test itself runs 4.7x slower than real time in this Python code, so the fair figure is 7.7x the sign test's wall time). A mixed-integer
formulation or a compiled solver would be needed to deploy it; the
enumeration shows the payoff is real but moderate.

> **Revision note (2026-09-11):** against the fourth rule the hybrid MPC
> looks worse than stated above: relinearizing at the measured state
> (one QP per step) gives downward overshoot 0.016 deg at baseline and
> 0.000 with 0.05 deg noise (five seeds), against the hybrid's 0.236 and
> 0.135 (one seed), with sinusoid RMS 0.149 vs 0.159 at baseline. The
> enumeration repairs part of the held-branch plan at the reference
> point; moving the linearization point repairs more at 1/81 of the cost.

### 9.8 Status of the section 0 inventory

| limitation | status |
|---|---|
| by-eye moment arms | **closed**: pixel-level digitization reproducing Murray's Table 2; exposed and fixed two further problems (section 9.1) |
| triceps moment arm | **characterized, not resolved**: literature does not pin it down; co-contraction cancellation holds for one of two values |
| one muscle / one constant set | **closed**: 3 flexors, 3 constant sets, same structure |
| no torque validation | **closed at group level**: 3.6 N m RMS vs Buchanan |
| rigid tendon | **closed**: elastic tendon keeps the structure and the error |
| exact state, no noise/delay, one rate, no mismatch | **closed**, with one negative finding (unmodelled delay >= 20 ms) |
| horizon-level branch assignment | **addressed**: enumeration MPC helps (0.37 -> 0.24 deg), not real-time; relinearization at the measured state does better (0.02 deg) at one QP per step |
| why feedback washes it out | **partly**: nu-gap fails; quadratic Lyapunov certificate holds for correct classification and 4 of 5 fixed designs |
| hardware, EMG, reflexes/history, systematic review, sampled-data PWA certificate | future work (section 8) |

## 10. Follow-up runs after the review (2026-09-11; predictions written before the runs)

Hardware validation is not available to this project (user, 2026-09-11).
Three limitations the ICRA draft states can still be closed in
simulation. Code: `mpc.simulate_realistic(..., compensate_delay=True)`,
mode `"blended_act"` in `mpc.py`, `scripts/phase7_followups.py`
(`results/phase7_followups.json`) and `scripts/phase7_realism.py --only
blended_act` (`results/phase7_realism_blended_act.json`).

### 10.1 Delay-aware MPC

Delay was the one negative closed-loop finding (section 9.6 revision
note). Standard remedy, a predictor: at step `k` the controller rolls the
measured state forward through the `d` inputs already sent but not yet
applied, with its nominal model, and plans against the reference preview
starting at `t_{k+d}`, when its input lands. On the nominal plant the
predictor is exact, so every rule must reproduce its undelayed baseline;
that is a correctness check (tested in `tests/test_mpc.py`), not a
finding. The informative case is a predictor with the wrong model: plant
activation time constants x1.25. **Prediction:** every rule stays stable
at 20 and 50 ms, the undelayed pattern returns (blended worst upward,
naive zero downward, sign test overshooting downward, relinearization
best on the sinusoid), and the numbers move toward the undelayed
tau x1.25 row. If a rule still diverges at 50 ms with a 25% model error,
compensation is not enough and the papers say so.

### 10.2 Weighting sensitivity

The comparison used one weighting, Q = diag(1000, 1, 0), R = 200. Rerun
the undelayed baseline with R x0.25 (aggressive), R x4 (conservative) and
q_thetadot x20 (velocity-damped). **Prediction:** the pattern is
structural, not a tuning artefact: under every weighting the blended
model overshoots the upward step most, the naive branch never overshoots
downward, and relinearization overshoots downward less than the sign
test. Magnitudes should scale: a smaller R leans harder on the model, so
wrong-model penalties grow; q_thetadot x20 damps every overshoot.

### 10.3 Smoothing only the activation switch

At any equilibrium, De Groote et al.'s tanh-smoothed activation has the
averaged activation row for any smoothing constant (ICRA Discussion;
checked against OpenSim's `DeGrooteFregly2016Muscle` source). New rule
`blended_act`: the sign test with its activation row replaced by that
average, so it differs from the sign test in the activation row alone and
from `blended` in the velocity row alone. **Prediction** from the papers'
own mechanism sentence ("its activation gain averages a fast and a slow
branch and matches neither"): most of the blended model's 0.65 deg upward
overshoot comes from the activation row, so `blended_act` overshoots
upward far more than the sign test (0.08) and close to `blended`. If it
stays near the sign test instead, that sentence is wrong and the velocity
blending carries the cost.

### 10.4 Results (35 runs, `results/phase7_followups.json`)

Overshoot up / down and sinusoid RMS (deg); N naive, B blended, S sign
test, R relinearize, BA blended activation only.

| condition | N | B | S | R | BA |
|---|---|---|---|---|---|
| delay 20 or 50 ms, predictor, nominal plant | 0.07 / 0.00 / 0.200 | 0.65 / 0.15 / 0.196 | 0.08 / 0.37 / 0.153 | 0.07 / 0.02 / 0.149 | 0.42 / 0.48 / 0.169 |
| delay 20 ms, predictor, plant tau x1.25 | 0.14 / 0.00 / 0.202 | 0.84 / 0.22 / 0.233 | 0.16 / 0.34 / 0.153 | 0.15 / 0.03 / 0.149 | 0.60 / 0.58 / 0.190 |
| delay 50 ms, predictor, plant tau x1.25 | 0.18 / 0.00 / 0.205 | 0.87 / 0.28 / 0.242 | 0.21 / 0.44 / 0.156 | 0.20 / 0.04 / 0.152 | 0.63 / 0.61 / 0.197 |
| R x0.25 | 0.07 / 0.00 / 0.163 | 0.68 / 0.00 / 0.174 | 0.08 / 0.23 / 0.137 | 0.06 / 0.01 / 0.133 | 0.46 / 0.28 / 0.152 |
| R x4 | 0.06 / 0.00 / 0.305 | 0.53 / 0.10 / 0.208 | 0.06 / 0.42 / 0.182 | 0.05 / 0.00 / 0.198 | 0.32 / 0.41 / 0.177 |
| q_thetadot x20 | 0.00 / 0.00 / 0.537 | 0.00 / 0.00 / 0.337 | 0.00 / 0.00 / 0.374 | 0.00 / 0.00 / 0.420 | 0.00 / 0.00 / 0.370 |

- **10.1 confirmed.** With the nominal plant the predictor reproduces the
  undelayed baseline to the printed digit at both delays (the
  correctness check). With a 25% model error in the predictor no rule
  diverges at 20 or 50 ms (uncompensated, blended and sign test diverge
  at 50 ms) and the undelayed pattern returns: blended worst upward,
  naive zero downward, sign test 0.34-0.44 downward, relinearization
  at most 0.04 downward and best on the sinusoid.
- **10.2 partly confirmed.** The overshoot pattern survives a fourfold
  change of R either way (blended upward 0.53-0.68; naive downward 0;
  relinearization downward at most 0.01 against the sign test's
  0.23-0.42). The sinusoid ranking does not: at R x4 the sign test
  (0.182) and blended-activation (0.177) track better than relinearization
  (0.198), and with q_thetadot x20 no rule overshoots, every rule tracks
  2-3.6x worse than at baseline, and the blended model tracks best
  (0.337). So "relinearization tracks best" holds at the baseline and
  aggressive weightings only; the papers say so.
- **10.3 confirmed.** Activation-only blending overshoots upward 0.42 deg,
  five times the sign test and (0.42-0.08)/(0.65-0.08) = 60% of the
  blended model's excess, and downward 0.48 deg, more than any other
  rule. The activation row carries most of the blended model's cost, and a
  controller built on a tanh-smoothed activation model inherits it.
- **blended_act under the realism sweep** (`phase7_realism.py --only
  blended_act`, 20 runs, `results/phase7_realism_blended_act.json`):
  upward overshoot 0.23-0.90 deg, above the sign test and below the
  blended model in every undelayed condition; downward 0.40-1.10 deg,
  above both in every undelayed condition (T_s 40 ms: 1.10); unmodelled
  delay 20 ms: 0.87 up / 1.69 down; 50 ms: diverges. Smoothing only the
  activation switch is not a cheap substitute for reading the branch.
