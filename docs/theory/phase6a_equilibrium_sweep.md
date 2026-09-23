# Phase 6a: Does the Branch Ambiguity Matter at *Every* Equilibrium?

> **Revision note (2026-09-11, Phase 7):** rerun with the pixel-level
> moment arm (PCHIP), the range-of-motion grid moved to 12.5-122.5 deg in
> 5-degree steps (off the 65 deg passive-curve kink) and the held-load
> sweep at 60 deg (`phase7_limitations.md` section 9.1). Gravity-only a*
> spans 0.028-0.115. Activating steps: both-wrong error never below 75% at
> 50 ms and 23% at 500 ms across the range of motion; on the held-load
> sweep never below 24% at 50 ms and 34% at 500 ms (wrong activation falls
> 78% -> 0% at a* = 0.884 -> 7% at 0.95; wrong velocity rises 2% -> 28% at
> 50 ms). **The "never below 20%" claim holds for activating steps only**:
> for deactivating steps at 500 ms it reaches 14% at the lightest load
> (a* = 0.02) and 22% across the range of motion, while at 50 ms
> deactivating errors stay at or above 30% (298% on the range of motion).

Tier 2 item T2.a of `docs/completion_plan.md`. Phase 4
(`phase4_quantification.md`) quantified the wrong-branch prediction error
at **one** equilibrium (`theta*=65 deg`, `a*=0.133`, gravity load only).
The abstract's structural claim is about every equilibrium; this phase
extends the quantitative claim to match. Code:
`scripts/phase6a_equilibrium_sweep.py`, `tests/test_equilibrium_sweep.py`.

## 1. What is known in closed form before running anything

### 1.1 Activation switch: the ratio is a function of `a*` only

From `phase2_linearization.md` section 3 (and `linearize.py`,
`_activation_row`), the two one-sided partials of `adot` w.r.t. `a` at
`u*=a*` are

```
activating:    d(adot)/da = -1 / (tau_act * (0.5 + 1.5 a*))
deactivating:  d(adot)/da = -(0.5 + 1.5 a*) / tau_deact
```

so their ratio is

```
rho_a(a*) = tau_deact / (tau_act * (0.5 + 1.5 a*)^2)
```

with `tau_act = 15 ms`, `tau_deact = 50 ms` (Thelen 2003, young adult):

| `a*` | `rho_a` | note |
|---|---|---|
| 0.00 | 13.33 | maximum (fully relaxed) |
| 0.133 | 6.81 | Phase 2/4 equilibrium |
| 0.30 | 3.70 | |
| 0.50 | 2.13 | |
| 0.884 | 1.00 | the two branches momentarily agree: `a* = (sqrt(tau_deact/tau_act) - 0.5)/1.5` |
| 1.00 | 0.83 | branches swap order (deactivation faster than activation) |

The ratio depends on nothing but `a*` and the two Thelen time constants:
not on the joint, the load, the moment arm, or `theta*`. **Two
predictions to test numerically below:** (i) wrong-activation-branch
prediction error should be largest at low `a*` and shrink toward
`a* ~ 0.88`; (ii) it should be *nonzero again* above 0.88 (the branches
cross, they do not merge), though small.

### 1.2 Force-velocity switch: the ratio is exactly 2, everywhere

From `linearize.py`, `_mechanical_row`, the one-sided slopes of
contractile force w.r.t. normalized fiber velocity at `v=0` are

```
shortening:   dF/dv = a* f_L * (1 + 1/A_f) / k
lengthening:  dF/dv = a* f_L * (2 + 2/A_f) / k,     k = (0.25 + 0.75 a*) V_max
```

`F_len` cancels out of the lengthening slope at `v=0` (it multiplies
both numerator and denominator of the Thelen Eq. 7 inverse there), so the
ratio is exactly 2 for every `a*`, `theta*`, and parameter choice within
Thelen's family. Unlike the activation ratio it never approaches 1. The
*absolute* size of both slopes scales with `a* f_L(theta*)`, so the
force-velocity kink's effect on the joint-angle prediction is expected
to grow with activation and with proximity to optimal fiber length.

### 1.3 What sets `a*`

Gravity-only equilibria: `a*(theta*)` solves
`tau_iso(theta*, a*) = m g l_c sin(theta*)`. The forearm+hand gravity
torque peaks at `m g l_c = 2.90 N*m` (`theta*=90 deg`), while brachialis
alone can produce ~25 N*m at `a=1` (Phase 1), so **gravity-only
equilibria all have low activation** (`a*` well below 0.2 across the ROM).
This is the postural regime a prosthesis or exosuit controller spends
most of its time in, and it is exactly where `rho_a` is largest.

To reach higher `a*` the joint must hold a load. The sweep uses a
**constant external torque** `tau_ext` chosen so that `(theta*, a*)` is
an equilibrium: `tau_ext = tau_iso(theta*, a*) - tau_g(theta*)`. A
constant torque adds nothing to the Jacobian, so `linearize.linearize(
theta_star, a_star, ...)` is already the correct branch-linearization at
such a point without modification; only the true-trajectory integrator
needs the extra term. (A hand-held mass would instead add a
`theta`-dependent term to `d(thetaddot)/d(theta)`; a constant torque was
chosen to isolate the `a*` dependence, and is physically a cable/pulley
weight or a torque motor. Stated so that nobody later reads this sweep
as "holding a dumbbell.")

## 2. Sweep design

Same protocol as Phase 4 (`T_s=10 ms`, step `delta_u=+/-0.01` from
equilibrium, one-shot linearization held fixed over the horizon, error
= `|theta_pred(T) - theta_true(T)|` reported as a fraction of the true
swing `|theta_true(T) - theta*|`), repeated over:

- **Sweep A (gravity only):** `theta*` from 10 to 125 deg, `a*` from the
  gravity balance.
- **Sweep B (held load):** `theta* = 65 deg` fixed, `a*` from 0.02 to
  0.95.

For each equilibrium and each horizon `T` in {50, 100, 200, 500} ms:
error fractions for matched / wrong-activation / wrong-velocity /
both-wrong branches, plus a check that the true trajectory actually stays
on the matched branch for the whole horizon (fraction of samples with
`u > a` and with `thetadot > 0`). Equilibria where it does not are
flagged in the output rather than silently included.

## 3. Results

(Filled in after the run; see below.)

Run 2026-08-24 on the remote (`python3 scripts/phase6a_equilibrium_sweep.py`,
9/9 tests in `tests/test_equilibrium_sweep.py`, including a check that the
generalized helpers reproduce Phase 4's headline numbers at `theta*=65 deg`
to within 2 percentage points). Every equilibrium in both sweeps stayed on
its matched branch for the full 500 ms horizon (occupancy check 1.00/1.00
throughout, with the tolerance described in `prediction.branch_occupancy`),
so "matched" vs. "wrong" is unambiguous everywhere below. Plots:
`results/phase6a_gravity_sweep_{pos,neg}.png`,
`results/phase6a_load_sweep_{pos,neg}.png`.

### 3.1 Sweep A: gravity-only equilibria (`theta*` = 10-125 deg)

`a*` ranges only from 0.049 (10 deg) to 0.133 (50-65 deg) and back to
0.099 (125 deg); `rho_a` from 6.8 to 10.2. Final-time error as a fraction
of the true swing, `delta_u=+0.01`:

| `theta*` | `a*` | horizon | matched | wrong activation | wrong velocity | both wrong |
|---|---|---|---|---|---|---|
| 10 deg  | 0.049 | 50 ms  | 0.5% | 77% | 1%  | 77% |
| 10 deg  | 0.049 | 500 ms | 0.9% | 16% | 10% | 23% |
| 65 deg  | 0.133 | 50 ms  | 0.8% | 70% | 13% | 73% |
| 65 deg  | 0.133 | 500 ms | 1.2% | 16% | 46% | 53% |
| 120 deg | 0.101 | 50 ms  | 1.4% | 72% | 18% | 77% |
| 120 deg | 0.101 | 500 ms | 2.7% | 14% | 48% | 55% |

- **Wrong-activation error is flat across the ROM** (70-77% at 50 ms,
  56-65% at 100 ms, 37-49% at 200 ms, 12-19% at 500 ms) because `a*`
  barely varies, exactly as section 1.1 predicts: it depends on `a*`
  only, not on `theta*`.
- **Wrong-velocity error grows with `theta*`** (1% at 10 deg to 18% at
  120 deg at 50 ms; 10% to 49% at 500 ms), tracking the growth of the
  absolute force-velocity slope with `a* f_L(theta*) r(theta*)`
  (section 1.2). Near full extension the joint is so lightly loaded that
  the force-velocity kink is nearly invisible at short horizons.
- **Both-wrong error is never below 23%** at any `theta*` or horizon.

`delta_u=-0.01` (deactivating step) is *not* a quantitative mirror image,
correcting the "mirror-image pattern" shorthand in
`phase4_quantification.md` section 2: the wrong-activation error is
**225-330% of the swing at 50 ms** and 112-184% at 100 ms across the ROM,
because the wrong branch there is the *faster* one (`rho_a` times faster
deactivation than reality), so it overpredicts the response instead of
underpredicting it. The qualitative conclusion is the same; the
magnitude on the deactivating side is 3-4x worse. The matched-branch
error also grows to 15-32% by 500 ms on this side (generic
linearization error, larger because the true trajectory travels further
when unloading), so at 500 ms the activation-branch choice no longer
dominates on the deactivating side either -- the force-velocity choice
does (32-49%).

### 3.2 Sweep B: held constant-torque load, `theta*=65 deg`, `a*` = 0.02-0.95

`tau_ext` ranges from -2.2 N*m (an assisting load, `a*=0.02`) to +16.1
N*m (`a*=0.95`). `delta_u=+0.01`, 50 ms and 500 ms horizons:

| `a*` | `rho_a` | wrong act (50 ms) | wrong vel (50 ms) | both (50 ms) | wrong act (500 ms) | wrong vel (500 ms) | both (500 ms) |
|---|---|---|---|---|---|---|---|
| 0.02  | 11.87 | 78% | 2%  | 79% | 20% | 21% | 37% |
| 0.133 | 6.81  | 70% | 13% | 73% | 16% | 46% | 53% |
| 0.30  | 3.69  | 56% | 19% | 64% | 12% | 53% | 57% |
| 0.50  | 2.13  | 38% | 22% | 51% | 7%  | 57% | 59% |
| 0.70  | 1.39  | 19% | 24% | 38% | 3%  | 60% | 61% |
| 0.884 | 1.00  | 0%  | 25% | 25% | 1%  | 63% | 63% |
| 0.95  | 0.90  | 7%  | 25% | 20% | 2%  | 64% | 63% |

(matched-branch error at most 1.7% at 50 ms and 2.0% at 500 ms across the
whole sweep.)

- **Prediction (i) confirmed:** wrong-activation error decreases
  monotonically with `a*`, from 78% to exactly 0% at `a*=0.884`, the
  closed-form unity point of `rho_a`.
- **Prediction (ii) confirmed:** above the unity point the branches have
  swapped order and the error reappears (7% at `a*=0.95`), small because
  `rho_a` only reaches 0.83 at `a*=1`.
- **The force-velocity error does the opposite** -- 2% to 25% at 50 ms,
  21% to 64% at 500 ms -- so the two kinks trade off: at low activation
  the activation switch dominates, at high activation the force-velocity
  switch does. **The both-wrong error is never below 20% at 50 ms and
  never below 37% at 500 ms, at any `a*` tested.**
- Deactivating side (`delta_u=-0.01`): same structure with larger
  magnitudes (wrong-activation 364% of swing at `a*=0.02`, 50 ms;
  wrong-velocity up to 159% at `a*=0.95`, 500 ms).

## 4. Answer

**The Phase 4 finding holds at every equilibrium tested, and its `a*`
dependence is exactly the closed-form one.** There is no equilibrium at
which silently picking one branch is safe: the activation ambiguity is
worst precisely in the low-activation postural regime that gravity-only
equilibria all occupy (`a*` at most 0.133 for this forearm), and where
it does vanish (`a* ~ 0.88`, i.e. near-maximal effort against a 15 N*m
load) the force-velocity ambiguity has grown to take its place. This is
stronger than what the abstract currently claims (one equilibrium,
13-88%); the honest headline is now "13-88% at the representative
gravity-only equilibrium; never below 20% of the true swing at any
equilibrium from `a*=0.02` to `0.95`, and 2-4x worse for deactivating
steps than activating ones."

Caveat: everything here is still open-loop one-shot prediction error,
the same quantity Phase 4 measured. Whether feedback washes it out in
closed loop is Phase 6b's question, not answered here.
