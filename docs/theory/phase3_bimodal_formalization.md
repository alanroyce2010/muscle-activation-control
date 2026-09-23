# Phase 3: Formalizing the Equilibrium Kinks as Bimodal Piecewise-Linear Switches

> **Revision note (2026-09-11, Phase 7):** the e vectors below were
> computed at the old 65 deg equilibrium with the by-eye moment-arm table.
> At the corrected 60 deg representative posture (`phase7_limitations.md`
> section 9.1): e_act = [0, 0, 86.29], e_vel = [0, -12.63, 0], residuals
> exactly 0, supports disjoint. The same structure holds for biceps and
> brachioradialis plants, all three Thelen constant sets, and the
> elastic-tendon model, where the velocity switch moves to the fibre row
> (section 9.3-9.4 there).

CLAUDE.md section 4.2 step 1. Builds on `docs/theory/phase2_linearization.md`
section 3 (the two kinks) -- read that first. Code:
`src/muscle_activation_control/bimodal.py`, verified by
`tests/test_bimodal.py` (6/6 pass, part of the project's 34/34 total).

## 1. The framework

**Camlibel, M.K., Heemels, W.P.M.H., Schumacher, J.M. (2008)**, "A full
characterization of stabilizability of bimodal piecewise linear systems
with scalar inputs," *Automatica* 44(5):1261-1267
(`docs/literature/camlibel2008_bimodal_piecewise_linear_stabilizability.pdf`,
read in full). Their system class:

```
xdot = A1*x + b1*u   if y <= 0
xdot = A2*x + b2*u   if y > 0
y = c^T*x + d*u
```

required *continuous* across the switching hyperplane `{y=0}` --
equivalent to existence of a vector `e` such that

```
A1 - A2 = e*c^T
b1 - b2 = e*d
```

(their Eq. 3). **What this project uses from that paper**: this
system-class definition and continuity characterization only -- *not*
their main result (Theorem 2.3, a stabilizability characterization via
geometric control theory, an existence question about whether *some*
input drives the state to the origin). That's a different question than
this project asks (CLAUDE.md section 1: does the branch ambiguity matter
*quantitatively* for a specific, practical control design, not whether
the plant is abstractly stabilizable).

**Scope honesty, restated from CLAUDE.md 4.2**: their formal system class
is strictly bimodal -- one scalar switching function `y`, two modes. This
project's linearized system has *two independent* switching functions
(section 2 below), hence up to 4 modes -- a genuine generalization their
own concluding remarks flag as future work ("multi-modal systems"),
pointing to a follow-up "conewise linear systems" paper whose publication
details this project has **not** independently verified or read. Section
3 below shows our specific case doesn't need that follow-up paper anyway,
because the two switches turn out to be structurally independent -- but
that's a property of *this* system, not a general result being claimed
from the unread paper.

## 2. Each switch, verified individually

Both kinks from `phase2_linearization.md` section 3 are checked against
the continuity condition above, **algebraically, not asserted** --
`bimodal.py`'s `_solve_continuity` solves for `e` directly (not a
least-squares fit that would silently "succeed" with a small but nonzero
residual even for a switch that isn't really this form) and reports the
residual, which must be exactly (to floating point) zero for a genuine
instance of this structure.

**Activation switch** (`y = u - a`, so `c = [0,0,-1]`, `d = 1`):

```
e_activation = [0, 0, 81.30]   (residual: 0.0, exact)
```

Only the third entry is nonzero -- meaning the activation switch's effect
on the linearization is confined entirely to the `(a, a)` and `(a, u)`
entries of `A` and `B` (the activation row), touching neither the `theta`
row nor the `thetadot` row. This matches the mechanical-row argument in
`phase2_linearization.md` section 3 directly (at fixed `thetadot*=0`, the
mechanical row's `theta` and `a` partials don't depend on activation
branch at all).

**Force-velocity switch** (`y = thetadot`, so `c = [0,1,0]`, `d = 0`):

```
e_velocity = [0, -11.47, 0]   (residual: 0.0, exact)
```

Only the second entry (the `thetadot` row's own `thetadot` column) is
nonzero. `d=0` here (the switching variable `thetadot` doesn't depend on
`u`), so the continuity condition's `b1-b2=e*d` reduces to `b1=b2`
exactly -- confirmed: `B` never depends on the velocity branch at all
(`linearize.py`'s `_activation_row` is the only source of `B`'s nonzero
entry). `e_velocity[1] = -11.47` exactly equals the shortening-branch
value of `A[1,1]` itself (`tests/test_bimodal.py::
test_velocity_switch_e_matches_katz_asymmetry`) -- the direct algebraic
consequence of the lengthening slope being exactly 2x the shortening
slope (Katz 1939, `phase2_linearization.md` section 3): `2x - 1x = 1x`.

## 3. Structural independence -- why the 4-mode case doesn't need the unread follow-up paper

`e_activation` and `e_velocity` have **disjoint support**: the activation
switch's `e` is nonzero only at index 2 (the `a` row), the velocity
switch's `e` is nonzero only at index 1 (the `thetadot` row).
(`tests/test_bimodal.py::test_the_two_switches_are_structurally_independent`,
verified, not assumed.) This means the full 4-mode linearized system
decomposes exactly as the **direct sum of two independent bimodal
switches** acting on non-overlapping matrix entries -- not a generically
coupled multimodal system where the switches interact. Each of the 4
`(A,b)` pairs computed in `linearize.py` can be written

```
A(act_branch, vel_branch) = A_base
    + [activating? e_activation*[0,0,1] : 0]
    + [lengthening? e_velocity*[0,1,0] : 0]
```

i.e. the two kinks' effects simply add, with no cross term. This is a
genuine finding about *this* system (brachialis + single pendulum joint
+ Thelen 2003 dynamics), not a general claim about all Hill-muscle
systems -- a system where, e.g., the force-velocity relation's shape
parameters themselves depended on activation in a way that coupled the
two switches' Jacobian contributions would not have this clean
decomposition, and would need the genuinely multimodal/conewise theory
this project has deliberately not reached for.

## 4. What this buys the project

A precise, literature-grounded name and formal justification for what
Phase 2 found empirically: "the" linearization at this equilibrium is
provably (not just observed to be) the direct sum of two independent
bimodal-piecewise-linear switches, each individually a verified instance
of a well-studied system class. This sets up Phase 4: quantifying whether
the resulting `~81.3` and `~11.5` (1/s, roughly) Jacobian-entry spreads
translate into a control-relevant prediction discrepancy at a realistic
sample time, or wash out.
