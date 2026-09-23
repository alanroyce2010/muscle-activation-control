# Phase 6c: Antagonist Pair (Brachialis + Triceps)

> **Revision note (2026-09-11, Phase 7):** rerun at 60 deg with the
> pixel-level moment arms. With the Murray-model triceps curve the
> force-velocity kink cancels at a_t* = 0.264 (a_b* = 0.514). **This
> cancellation is not robust**: with the Holzbaur-consistent triceps
> moment arm (scaled by 0.81) no balance exists within feasible
> co-contraction and the ratio only falls to 1.068 (`phase7_limitations.md`
> section 9.2). The robust claim is that co-contraction shrinks the kink
> from 2 to within about 7% of 1. At very high co-contraction (a_t* = 0.5,
> a_b* = 0.87) every branch error falls to about 2%, because brachialis
> activation approaches the 0.884 point where rho_a = 1.

Tier 2 item T2.c of `docs/completion_plan.md`; the "antagonist pair"
future work the submitted abstract promised. Question: does the
equilibrium non-smoothness structure (Phase 2/3) survive adding a second
muscle, does the direct-sum decomposition survive, and how does
co-contraction change the size of each kink? Code:
`src/muscle_activation_control/antagonist.py`,
`scripts/phase6c_antagonist_pair.py`, `tests/test_antagonist.py`.

## 1. Second muscle: lumped lateral + medial triceps

**Choice**: the two single-joint heads of triceps, lateral (`TRIlat`)
and medial (`TRImed`), lumped into one extensor. The long head
(`TRIlong`) crosses the shoulder and is excluded for exactly the reason
biceps was excluded in favour of brachialis (CLAUDE.md 6.2): a
single-joint muscle's torque is `F * r(theta)` with no shoulder angle to
hold fixed. Murray et al. (1995) p. 520 report "no substantial
differences" between the moment-arm curves of the three heads, which is
what makes lumping the two single-joint heads reasonable.

Parameters, Holzbaur et al. (2005) Table 1 (read from the PDF, section
6.4 of CLAUDE.md):

| | `TRIlat` | `TRImed` | lumped extensor used here |
|---|---|---|---|
| PCSA (cm^2) | 4.5 | 4.5 | 9.0 |
| peak force `F_0^M` (N) | 624.3 | 624.3 | **1248.6** (sum) |
| optimal fiber length `l_0^M` (cm) | 11.4 | 11.4 | **11.4** |
| tendon slack length (cm) | 9.8 | 9.1 | 9.45 (mean; unused with a rigid tendon except to fix the reference length, section 1.2) |
| pennation `alpha_0` (deg) | 9 | 9 | **9** |
| `ma_avg` (cm) | -2.1 | -2.1 | -2.1 (negative = extension) |

Same Thelen (2003) young-adult activation and force-velocity constants
as brachialis (`tau_act=15 ms`, `tau_deact=50 ms`, `A_f=0.25`,
`F_len=1.4`, `V_max=10 l_0/s`, `gamma=0.45`, `k_PE=5`, `eps_0^M=0.6`).

### 1.1 Triceps moment arm `r_tri(theta)`

Same method as brachialis (`single_joint_dynamics.md` section 3):
digitized from Murray et al. (1995) Fig. 4, **Model** panel, `TRI`
curve (p. 519; the text on p. 520 pins the shape: the model curve "peaks
(~3.0 cm) in the first 25 deg of flexion, and decreases by 25% between
25 and 120 deg", and Table 2 gives the model's max-min percent
difference as 25%). Visual reads, +/-0.15-0.2 cm each, extension
moment arm reported as a positive magnitude:

| theta (deg) | 0 | 20 | 40 | 60 | 80 | 100 | 120 | 130 |
|---|---|---|---|---|---|---|---|---|
| `r_tri` (cm) | 2.8 | 3.0 | 2.9 | 2.75 | 2.6 | 2.4 | 2.25 | 2.2 |

Trapezoidal average over 0-130 deg: **~2.6 cm**, against Holzbaur
(2005)'s `ma_avg = 2.1 cm` for both heads -- a **24% gap**, three times
the 8% gap the brachialis digitization had. The two sources are
different models (Murray 1995's triceps path is the medial head over a
polynomial-fit tendon excursion; Holzbaur 2005's is a wrap-cylinder
path in a different musculoskeletal model), so some disagreement is
expected, but this one is not within digitization error. **Decision:
use the Murray digitization**, for consistency with how brachialis was
done and because it is the curve with a published shape; and record the
consequence honestly: triceps torque and its force-velocity damping
contribution both scale linearly with `r_tri` (damping with `r_tri^2`),
so the numbers in section 3 that depend on triceps strength are
uncertain by roughly that factor. **None of the structural results
(section 2) depend on `r_tri` at all.**

### 1.2 Rigid-tendon kinematics with pennation

Triceps lengthens with flexion: `dl^MT_tri/dtheta = +r_tri(theta)` (the
opposite sign to brachialis). Constant-width pennation model (Zajac
1989; Millard 2013 rigid-tendon form), `w = l_0^M sin(alpha_0)`:

```
l^M     = sqrt( (l^MT - l_s^T)^2 + w^2 )           (closed form, no solve)
cos(alpha) = (l^MT - l_s^T) / l^M
v^M     = v^MT * cos(alpha)                        (from d/dt of l^M sin(alpha) = w)
F^MT    = F^M * cos(alpha)
```

With `alpha_0 = 9 deg`, `cos(alpha)` stays within about 0.985-0.99
across the ROM: a ~1% effect, kept exact because it costs three lines.
Reference length: same mid-ROM convention as brachialis, `l^M_tri(65
deg) = l_0^M`, so `l^MT_tri(theta) = l_s^T + l_0^M cos(alpha_0) +
integral_65deg^theta r_tri(phi) dphi`.

### 1.3 Pair dynamics

State `x = [theta, thetadot, a_b, a_t]`, input `u = [u_b, u_t]`:

```
I thetaddot = tau_bra(theta, thetadot, a_b) - tau_tri(theta, thetadot, a_t) - m g l_c sin(theta)
adot_b = (u_b - a_b)/tau_a(a_b, u_b),   adot_t = (u_t - a_t)/tau_a(a_t, u_t)
tau_tri = F^MT_tri * r_tri,   v^M_tri = +r_tri thetadot cos(alpha)      (flexion lengthens triceps)
tau_bra as before,            v^M_bra = -r_bra thetadot                  (flexion shortens brachialis)
```

Equilibria: `thetadot* = 0`, `u_b* = a_b*`, `u_t* = a_t*`, and torque
balance. With two muscles the equilibrium set at a given `theta*` is a
one-parameter family: fix the co-contraction level `a_t*` and solve
`tau_bra(theta*, 0, a_b*) = tau_g(theta*) + tau_tri(theta*, 0, a_t*)`
for `a_b*` (brentq, as before).

## 2. Switch structure: what changes and what does not

Three scalar switching functions, not four:

| switch | `y` | `c^T` | `d` |
|---|---|---|---|
| brachialis activation | `u_b - a_b` | `[0, 0, -1, 0]` | `[1, 0]` |
| triceps activation | `u_t - a_t` | `[0, 0, 0, -1]` | `[0, 1]` |
| force-velocity (both muscles) | `thetadot` | `[0, 1, 0, 0]` | `[0, 0]` |

**The two muscles' force-velocity kinks share one switching function.**
`thetadot > 0` means brachialis shortening *and* triceps lengthening;
`thetadot < 0` the reverse. They cannot switch independently, so this
is one bimodal switch, and the pair system has `2^3 = 8` modes, not 16.
(`docs/completion_plan.md` worried the two force-velocity `e` vectors
would overlap on the `thetadot` row and break the direct sum; they do
overlap, but as a single switch, so nothing breaks.)

Predictions to verify numerically:

1. Each of the three switches satisfies Camlibel et al.'s continuity
   condition exactly (`A1 - A2 = e c^T`, `b1 - b2 = e d`), as in Phase 3.
2. `e_b` touches only the `a_b` row, `e_t` only the `a_t` row, `e_v` only
   the `thetadot` row: disjoint support, so the 8-mode system is again
   the direct sum of independent bimodal switches, and the single-switch
   2008 Automatica paper still suffices.
3. Activation ratios are unchanged and per-muscle: `rho_a(a_b*)` and
   `rho_a(a_t*)` from Phase 6a section 1.1.
4. **The force-velocity ratio is no longer 2.** Let `D_i > 0` be muscle
   `i`'s damping contribution on its own *shortening* side,

   ```
   D_i = r_i^2 cos^2(alpha_i) F_0,i a_i* f_L,i (1 + 1/A_f) / (k_i l_0,i),   k_i = (0.25 + 0.75 a_i*) V_max
   ```

   (lengthening side is `2 D_i`, Katz). Then

   ```
   thetadot > 0 (flexing):    d(thetaddot)/d(thetadot) = -(D_b + 2 D_t)/I
   thetadot < 0 (extending):  d(thetaddot)/d(thetadot) = -(2 D_b + D_t)/I
   ratio (extending / flexing, the Phase 2 lengthening/shortening convention)
         = (2 D_b + D_t)/(D_b + 2 D_t)   in (1/2, 2)
   ```

   so `e_v[1] = (D_b - D_t)/I`: the kink is the *difference* of the two
   muscles' damping contributions. It is the single-muscle value 2 when
   `D_t = 0`, **exactly 1 (no kink) when `D_b = D_t`**, and below 1 when
   the triceps dominates. (First draft of this section had the ratio
   written upside down; caught by the test against finite differences,
   corrected in place.) Co-contraction does not remove the activation
   kinks (it adds one) but it can tune the force-velocity kink away.
   Physically: what one muscle loses in lengthening-side damping the
   other gains, and at balanced damping the joint sees the same total
   damping in both directions.

## 3. Quantification plan

Same protocol as Phase 4/6a (`T_s = 10 ms`, held excitation step, one-shot
branch-linearization vs. true nonlinear evolution, `theta` error as a
fraction of the true swing), at `theta* = 65 deg`, with a brachialis
excitation step `delta_u_b = +/-0.01` and the triceps excitation held at
`u_t = a_t*`. Co-contraction sweep `a_t*` in {0, 0.02, 0.05, 0.1, 0.2,
0.3}. For each: the matched 8-mode branch, the three single-wrong-switch
branches, and all-wrong. Note the triceps activation switch is a
different animal here: with `u_t` held at `a_t*`, the true trajectory
sits exactly *on* `y_t = 0` for the whole horizon, so neither triceps
activation branch is "wrong" in the sense of Phase 4 -- both predict
`adot_t = 0` -- and the wrong-triceps-branch error should be zero. That
is a prediction, not an assumption; the script checks it.

## 4. Results

Run 2026-08-24 on the remote (`python3 scripts/phase6c_antagonist_pair.py`;
`tests/test_antagonist.py` 7/7, full suite 65/65). Plots:
`results/phase6c_triceps_validation.png`,
`results/phase6c_cocontraction_sweep.png`.

### 4.1 Triceps validation

- Digitized `r_tri`: peak 3.00 cm at 20 deg, ROM average **2.66 cm**,
  i.e. **+27%** against Holzbaur (2005)'s `ma_avg = 2.1 cm` (section 1.1
  flagged this before the run; it is the one soft number in this phase).
- Fiber length over the ROM `l_bar` 0.72-1.24, `cos(alpha)` 0.976-0.992.
- Isometric extension torque at `a=1`: peak **34.8 N*m at 40 deg**
  (brachialis flexion peak 25.2 N*m at 108 deg, Phase 1). Two heads of
  triceps at ~35 N*m is in the range of whole-triceps maximal elbow
  extension in the literature (roughly 40-70 N*m for young adults with
  all three heads), so the magnitude is plausible but not a precise
  match to any one study; a proper validation target for the triceps
  was not pulled from the literature in this phase.
- Passive triceps torque at `a=0`: at most 1.19 N*m, at full flexion.

### 4.2 Switch structure (`theta* = 65 deg`)

| `a_t*` | `a_b*` | `rho_a(a_b*)` | `rho_a(a_t*)` | `D_b` | `D_t` | FV ratio, closed form | FV ratio, finite diff. | max continuity residual | disjoint |
|---|---|---|---|---|---|---|---|---|---|
| 0.00 | 0.133 | 6.81 | 13.33 | 0.875 | 0.000 | 2.0000 | 2.0000 | 1.2e-9 | yes |
| 0.02 | 0.167 | 5.92 | 11.87 | 1.024 | 0.297 | 1.4497 | 1.4497 | 1.0e-9 | yes |
| 0.05 | 0.218 | 4.88 | 10.08 | 1.212 | 0.684 | 1.2050 | 1.2050 | 7.9e-10 | yes |
| 0.10 | 0.303 | 3.66 | 7.89 | 1.460 | 1.209 | 1.0645 | 1.0645 | 5.5e-10 | yes |
| **0.1555** | **0.396** | | | **equal** | **equal** | **1.0000** | | | yes |
| 0.20 | 0.472 | 2.28 | 5.21 | 1.798 | 1.965 | 0.9708 | 0.9708 | 2.9e-10 | yes |
| 0.30 | 0.641 | 1.56 | 3.69 | 2.019 | 2.483 | 0.9336 | 0.9336 | 1.8e-10 | yes |
| 0.50 | 0.980 | 0.86 | 2.13 | 2.290 | 3.145 | 0.9004 | 0.9004 | 7.7e-11 | yes |

All four section 2 predictions hold:

1. Continuity residuals are at finite-difference precision (1e-9 to
   1e-11; the Phase 3 check was exact because it used analytic
   Jacobians, this one is central differences on the branch-forced
   model, see `antagonist.linearize`).
2. Disjoint support at every co-contraction level, e.g. at `a_t*=0.1`:
   `e_bra = [0, 0, 50.8, 0]`, `e_tri = [0, 0, 0, 89.6]`, `e_v = [0, 3.28,
   0, 0]`. **The 8-mode pair system is the direct sum of three
   independent bimodal switches**, so the 2008 Automatica single-switch
   framework still covers it with no appeal to the conewise follow-up.
3. Activation ratios are per-muscle and unchanged in form.
4. The force-velocity ratio is the closed form to four decimals, 2 with
   the triceps off, and **crosses 1 at `a_t* = 0.1555` (`a_b* = 0.396`)**,
   where the kink cancels. Above that, extension is *less* damped than
   flexion (ratio 0.90 at `a_t*=0.5`).

### 4.3 Wrong-branch prediction error vs co-contraction (`delta_u_b = +0.01`)

Final `theta` error as a fraction of the true swing; `wt` = wrong
triceps-activation branch, which equals the matched error exactly at
every point (its excitation is held on its switch, both branches predict
`adot_t = 0`; `a_t` drift in the true trajectory 0.0), as section 3
predicted.

| `a_t*` | `a_b*` | horizon | matched | wrong bra. activation | wrong velocity | all wrong |
|---|---|---|---|---|---|---|
| 0.00   | 0.133 | 50 ms  | 0.9% | 70% | 13% | 73% |
| 0.00   | 0.133 | 500 ms | 2.7% | 12% | 43% | 51% |
| 0.10   | 0.303 | 50 ms  | 0.9% | 55% | 4%  | 56% |
| 0.1555 | 0.396 | 50 ms  | 0.9% | 47% | 1%  | 47% |
| 0.1555 | 0.396 | 500 ms | 5.0% | 12% | 5%  | 12% |
| 0.20   | 0.472 | 50 ms  | 0.8% | 40% | 1%  | 39% |
| 0.20   | 0.472 | 500 ms | 4.4% | 10% | 1%  | 7%  |
| 0.30   | 0.641 | 500 ms | 3.6% | 7%  | 5%  | 1%  |
| 0.50   | 0.980 | 50 ms  | 0.4% | 9%  | 7%  | 17% |
| 0.50   | 0.980 | 500 ms | 2.7% | 2%  | 12% | 13% |

- **Co-contraction shrinks the activation kink's consequence**, because
  it raises `a_b*` and `rho_a(a_b*)` falls: wrong-brachialis-activation
  error 70% -> 9% (50 ms) across the sweep, tracking Phase 6a's `a*`
  dependence, now reached by co-contracting instead of by loading.
- **Co-contraction can cancel the force-velocity kink's consequence.**
  Wrong-velocity error falls from 13% to about 1% at 50 ms near the
  balance point and stays below 5% at 500 ms there (vs 43% with the
  triceps off), then grows again slowly. The minimum in the sweep sits
  at `a_t*=0.2` rather than exactly 0.1555 because the true trajectory
  leaves the equilibrium where the balance was computed.
- **Nothing makes all-wrong negligible at short horizons**: the minimum
  over the sweep at 50 ms is 17% (`a_t*=0.5`); the all-wrong error can
  get to ~1% only at one horizon (500 ms, `a_t*=0.3`) where the two
  kinks' errors happen to cancel each other.
- Deactivating side (`delta_u_b=-0.01`): same structure, larger
  magnitudes (wrong-brachialis-activation 225% at `a_t*=0`, 85% at the
  balance point, 9% at `a_t*=0.5`, all at 50 ms), as in Phase 6a.

## 5. Answer

Adding the antagonist changes the size of the kinks, not their
existence or their structure. Every equilibrium of the pair still sits
on three switching manifolds; the local model is still a direct sum of
independent bimodal switches, so the formal result of Phase 3 carries
over without modification. What is new: (i) a third switch (the
triceps' own activation kink) appears, with its own `rho_a(a_t*)` that
is *largest* at low co-contraction, and (ii) the force-velocity kink is
no longer a fixed factor of 2 but `(2 D_b + D_t)/(D_b + 2 D_t)`, which
co-contraction can drive to 1. Since real posture control uses
co-contraction, this gives a concrete, testable mechanism by which the
nervous system (or a controller) could make a muscle-driven joint
*locally smoother* in its velocity dynamics -- at the cost of the
activation-switch ambiguity of a second muscle. Whether that is a
biological strategy or a coincidence of this model is well outside this
project's scope, but it is the kind of prediction the equilibrium
non-smoothness view produces that the smooth-model view cannot.

The soft spot, stated again: the triceps moment arm used here averages
27% above Holzbaur (2005)'s value, so the balance point `a_t*=0.1555`
is uncertain by roughly that factor (a smaller `r_tri` pushes it to
higher co-contraction). The structural results do not depend on it.
