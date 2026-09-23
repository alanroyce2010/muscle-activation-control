# Phase 2: Linearization and ZOH Conditioning Check

> **Revision note (2026-09-11, Phase 7):** the numbers below use the
> original by-eye moment-arm table at theta* = 65 deg. Phase 7
> (`phase7_limitations.md` section 9.1) replaced that table with a
> pixel-level digitization (PCHIP interpolation) and found that the old
> 65 deg equilibrium was **open-loop unstable** (slowest eigenvalue
> +0.17/s, a slow mode the eigenvalue statement here omitted) and sat on
> the clipped passive curve's kink (l_bar = 1). The representative posture
> is now 60 deg: a* = 0.1126, rho_a = 7.45, branch eigenvalues -99.7 to
> -0.3/s (all stable), max|A_d| = 7.72 over T_s 1-100 ms. The conditioning
> conclusion (informative non-transfer) and the two-switch finding are
> unchanged.

CLAUDE.md section 4, step 2. This is the actual core question of the
project: does the discretized full-state model blow up numerically the
way `origami-arm-control`'s pneumatic pressure state did
(`mpc_control.md` section 2, `max|A_d| ~ 900,000` at a practical sample
time)? Builds on Phase 1 (`src/muscle_activation_control/{muscle,joint}.py`).

## 1. Nonlinear model recap

`x = [theta, thetadot, a]`, input `u` = excitation:

```
thetadot_dot = thetaddot = (tau_muscle(theta,thetadot,a) - tau_gravity(theta)) / I_elbow
adot         = (u-a) / tau_a(a,u)
```

(the trivial first row, `d(theta)/dt = thetadot`, is state 1's own
derivative, included for completeness of the 3-state system).

## 2. Equilibrium point

Pick `theta* = 65 deg` (the reference angle from Phase 1 where `l_bar=1`
by construction, so `f_L=1` and `F_PE=0` there -- keeps the equilibrium
calculation simple and avoids conflating "where on the force-length curve
am I" with "what does linearizing at equilibrium look like"), `thetadot*
= 0`. Isometric equilibrium (`F_ce = a*f_L` at `V=0`, regardless of
branch -- see section 3) requires:

```
tau_muscle(theta*, 0, a*) = tau_gravity(theta*)
a* * f_L(l_bar(theta*)) * F0_M * r(theta*) = m*g*l_c*sin(theta*)
a* = tau_gravity(theta*) / (F0_M * r(theta*))          [since f_L(l_bar=1)=1]
```

`adot=0` additionally requires `u* = a*` exactly -- this is not a free
choice, it's forced by the activation ODE's structure (both branches of
`tau_a(a,u)` give `adot=0` only when `u=a`). This is where the
non-smoothness in section 3 becomes unavoidable: **every physically valid
equilibrium of this system sits exactly on the activation ODE's
switching manifold.**

Exact numeric value of `a*` computed in code
(`src/muscle_activation_control/linearize.py`, not hand-derived here, to
avoid the kind of rounding slip already caught once this session --
CLAUDE.md 6.3's `f_L(0.5)` value).

## 3. Non-smoothness at equilibrium (the real finding here)

This is a structural difference from `origami-arm-control`'s pneumatic
case, and it is worth stating plainly: **the pneumatic pressure ODE,
`Pdot = -lambda_P(P-P_cmd)`, is linear and smooth everywhere -- its
equilibrium linearization was completely unambiguous. This project's
model is not smooth at any of its own equilibria, in two independent
ways, both traceable to well-established classical muscle-mechanics
asymmetries cited by Thelen (2003) itself, not to modeling artifacts or
bugs:**

**(a) Activation dynamics** (Winters 1995, via Thelen 2003 Eq. 2):
`tau_a(a,u)` switches between `tau_act*(0.5+1.5a)` (`u>a`) and
`tau_deact/(0.5+1.5a)` (`u<=a`). At the equilibrium (`u*=a*`), the
function value (`adot=0`) is continuous across the switch, but its
*slope* is not, since `tau_act != tau_deact`. Verified symbolically
(`scripts/_verify_branch_derivatives.py`, sympy) against these hand
formulas:

```
activating branch  (u>a):  d(adot)/da|eq = -1/(tau_act*(0.5+1.5a*))
                            d(adot)/du|eq =  1/(tau_act*(0.5+1.5a*))
deactivating branch (u<=a): d(adot)/da|eq = -(0.5+1.5a*)/tau_deact
                            d(adot)/du|eq =  (0.5+1.5a*)/tau_deact
```

At `tau_act=15ms, tau_deact=50ms` these differ by roughly a factor of
`tau_deact/tau_act * (0.5+1.5a*)^2` -- several-fold, not a rounding-level
difference. **There is no single "the" fast eigenvalue at this
equilibrium; there are two, depending on which direction a perturbation
pushes the excitation relative to activation.**

**(b) Force-velocity relation** (Katz 1939, via Thelen 2003 Eq. 6-7):
the contractile force `F_ce(v_norm)` switches between a "shortening"
branch (`v_norm<=0`) and "lengthening" branch (`v_norm>0`). Both agree
`F_ce=a*f_L` at `v_norm=0` (continuous), but Thelen (2003) explicitly
states the lengthening-side slope is *designed* to be twice the
shortening-side slope at zero velocity, citing Katz (1939) -- confirmed
symbolically here, exactly `2x`, not approximately. Since any equilibrium
requires `thetadot*=0` (hence `v_norm*=0` too, via
`v_norm=-r(theta)*thetadot/L0_M`), **this kink is also unavoidable at
equilibrium.**

Fortunately, tracing through the chain rule shows this second kink only
affects one entry of the mechanical row (`d(thetaddot)/d(thetadot)`) --
`d(thetaddot)/d(theta)` and `d(thetaddot)/d(a)` turn out to be safe to
compute with ordinary central finite differences, because at fixed
`thetadot*=0`, `v_norm` stays *exactly* zero for any perturbation of
`theta` or `a` alone (the kink is only crossed by perturbing `thetadot`
itself). See `linearize.py` docstring for the full argument.

**Practical consequence:** "the" linearization at this equilibrium is not
a single matrix -- it's up to 4 nearby linearizations (2 activation
branches x 2 velocity branches). All 4 are computed and reported below,
not just one picked arbitrarily. **Naive central finite differences
straddling either kink would silently return some physically meaningless
blend of the two branches -- flagging this explicitly because it is
exactly the kind of quiet, plausible-looking error this project's working
notes (CLAUDE.md section 2) already warn about repeating.**

## 4. Linearized matrices

```
A = [[0,                    1,                          0                ],
     [d(thetaddot)/d(theta), d(thetaddot)/d(thetadot)*, d(thetaddot)/d(a)],
     [0,                    0,                          d(adot)/d(a)*    ]]

B = [[0], [0], [d(adot)/d(u)*]]
```

(`*` = branch-dependent, section 3). Off-branch entries via central
finite differences (`src/muscle_activation_control/linearize.py`).

## 5. ZOH discretization and conditioning check

`scipy.signal.cont2discrete((A,B,C,D), T_s, method="zoh")`, matching the
pattern in `origami-arm-control/src/origami_arm_control/control/mpc_control.py`
(CLAUDE.md section 2). Sweep across candidate control-rate sample times
(e.g. `T_s in {1, 2, 5, 10, 20, 50} ms`, spanning what's actually
plausible for a real-time controller) and report `max|A_d|` for each of
the 4 branch combinations at each `T_s` -- the same diagnostic
`mpc_control.md` section 2 used to catch the pneumatic pressure blowup.

Numeric results: `scripts/phase2_conditioning.py` output, not
hand-computed here.

## 6. What would count as "transfer," "partial transfer," or "non-transfer"

Per CLAUDE.md section 4 step 4:

- **Clean transfer**: `max|A_d|` blows up at practical `T_s` (mirroring
  the pneumatic `~900,000` finding) *and* the reduced/algebraic-
  substitution model (Phase 3) fixes it the same way.
- **Partial transfer**: it blows up, but the reduced model's fix looks
  different in some material way (e.g. because of the branch
  non-smoothness in section 3 -- there may not be a single well-defined
  "fast state" to substitute out the way pressure was).
- **Informative non-transfer**: it does *not* blow up at practical
  `T_s` -- i.e. activation dynamics is not numerically stiff the way
  pressure was, despite structural (algebraic-entry) similarity. Still a
  reportable, useful result per CLAUDE.md section 2's explicit framing.

Not pre-judging which outcome this is -- that's what running section 5's
sweep is for.
