# Single-Joint Dynamics: Elbow + Brachialis

Phase 1 derivation (CLAUDE.md section 4, step 1). Builds directly on
CLAUDE.md section 6 (formulation, parameters, validation targets) --
read that first. No simulation code written yet; this is the
derive-before-code checkpoint per section 7.

## 1. System and coordinates

Single revolute joint (elbow), one muscle (brachialis, `BRA`), one rigid
body (forearm + hand) rotating under gravity. Matches CLAUDE.md section 3
("one muscle against a known load, pendulum-style").

- `theta`: elbow flexion angle. `theta = 0` at full extension, `theta =
  130 deg` (2.269 rad) at full flexion -- same convention as the
  `r_elbow_flex` coordinate in Holzbaur (2005)/`arm26.osim` (CLAUDE.md
  6.2), so the ROM bound is literature-consistent, not invented.
- Gravity configuration: **upper arm fixed, hanging vertically at the
  side** (shoulder neutral). At `theta = 0`, the forearm hangs straight
  down, aligned with the upper arm -- gravity torque is zero there
  (stable equilibrium with no muscle activation). Flexion swings the
  forearm up through horizontal (`theta = 90 deg`, maximum gravity
  torque) toward the upper arm (`theta = 130 deg`). This is a modeling
  choice, not something pinned by the moment-arm literature in section 6
  (which reports isometric capacity, not a specific task) -- **flagging
  for review**, since "arm hanging at the side, curling up" is the most
  standard single-joint elbow setup in the Hill-muscle-control
  literature, but it is a choice, not a derived fact.
- `a`: brachialis activation, `[0, 1]`, dynamics per CLAUDE.md 6.1.

State `x = [theta, thetadot, a]`. Input `u` = neural excitation, `[0,
1]`.

## 2. Gravity torque

Forearm + hand modeled as a rigid pendulum about the elbow axis:

```
tau_gravity(theta) = m * g * l_c * sin(theta)
```

opposing flexion for `0 < theta < 180 deg` (always true here, ROM caps at
130 deg).

**Pinned to de Leva, P. (1996)**, "Adjustments to Zatsiorsky-Seluyanov's
segment inertia parameters," *J Biomech* 29(9):1223-1230 (freely
available, not paywalled; saved at
`docs/literature/deleva1996_segment_inertia_parameters.pdf`) -- primary
source, superseding the earlier Winter (2009) textbook placeholder in
this section's first draft. de Leva's Table 4 reports forearm and hand as
**separate** segments (not a combined "forearm+hand" segment the way
Winter's simplified table does), from gamma-ray-scan measurements of 100
male and 15 female Caucasian subjects (Zatsiorsky et al. 1990a),
adjusted by de Leva to reference segment endpoints to joint centers.
Values for males (`m_body = 73.0 kg`, stature `1.741 m`, de Leva Table 4
caption -- these are the actual measured subject-sample values, not a
round assumed number the way the first draft's 75 kg was):

```
Forearm: length (EJC to WJC)     = 268.9 mm
         mass fraction            = 1.62% of body mass
         CM position from EJC     = 45.74% of forearm length
         sagittal radius of gyration (about own CM) = 27.6% of forearm length

Hand:    length (WJC to 3rd MET) = 86.2 mm
         mass fraction            = 0.61% of body mass
         CM position from WJC     = 79.00% of hand length
         sagittal radius of gyration (about own CM) = 62.8% of hand length
```

(EJC = elbow joint center, WJC = wrist joint center, both cross-checked
against de Leva's own Fig. 1 diagram, which states the forearm CM
position unambiguously as "45.74%" from the elbow for males -- exact
match to the table value, so the table row alignment used here is
verified, not just OCR-trusted. The two radius-of-gyration figures above
were read from the same table row/column position but were **not**
independently cross-checked against a second, unambiguous source the way
the mass and CM figures were via Fig. 1 -- flagging that distinction
honestly, since this document is trying not to repeat the
overconfident-secondhand-number mistake from CLAUDE.md 6.1's revision
note.)

Combining forearm + hand into one rigid body about the elbow axis via the
parallel axis theorem (`I_segment_cm = m*(r_sagittal)^2`, then `I_about_elbow
= I_cm + m*d^2` where `d` is each segment's own CM distance from the
elbow):

```
m_forearm = 0.0162 * 73.0 = 1.1826 kg
m_hand    = 0.0061 * 73.0 = 0.4453 kg
m         = m_forearm + m_hand = 1.628 kg

d_forearm (elbow to forearm CM) = 0.4574 * 0.2689 m = 0.1230 m
d_hand    (elbow to hand CM)    = 0.2689 + 0.79*0.0862 = 0.3370 m
                                   (forearm length + hand CM offset from WJC)

l_c = (m_forearm*d_forearm + m_hand*d_hand) / m
    = (1.1826*0.1230 + 0.4453*0.3370) / 1.628 = 0.1815 m

I_forearm_cm = m_forearm*(0.276*0.2689)^2 = 0.00651 kg*m^2
I_hand_cm    = m_hand*(0.628*0.0862)^2    = 0.00131 kg*m^2
I_elbow = [I_forearm_cm + m_forearm*d_forearm^2] + [I_hand_cm + m_hand*d_hand^2]
        = 0.02440 + 0.05189 = 0.0763 kg*m^2
```

Result: `m = 1.628 kg`, `l_c = 0.1815 m`, `I_elbow = 0.0763 kg*m^2` -- for
reference, the discarded Winter-table placeholder from the first draft of
this section gave `m=1.65 kg, l_c=0.177 m, I_elbow=0.0805 kg*m^2`, all
within about 3-6% of these de Leva-derived values, a reassuring
cross-check between two independent anthropometric sources even though
only one of them is now the pinned, citable number.

These are parameters for a "representative" subject, not a validation
target -- section 6's torque-angle target is normalized by the muscle's
own force capacity and doesn't depend on these at all (see section 6
below). Free to change (e.g. to female parameters, also in de Leva Table
4) once real simulation runs make a specific choice matter.

## 3. Muscle-tendon geometry: `r(theta)`

This is the one piece CLAUDE.md section 6.2 flagged as unresolved: exact
`r(theta)` requires the `arm26.osim` wrap-cylinder geometry (radius,
pose) that Phase 0 didn't extract, because brachialis wraps over the
elbow in that model and a straight origin-insertion line would cut
through bone at high flexion.

**Resolution: don't reconstruct the wrap geometry. Digitize the actual
published `BRA` moment-arm-vs-angle curve** from Murray (1995) Fig. 4
(`docs/literature/murray1995_elbow_moment_arms.pdf`, p. 519) and use it
directly, rather than fitting a 2-parameter algebraic shape to just the
peak and average summary numbers (the original v1 draft of this section
did the latter; superseded here after re-reading the actual figure).

Fig. 4 has three panels (Male Specimen, Female Specimen, Model), each
plotting `BRA` moment arm (cm) vs. elbow flexion angle (deg), gridlines
every 2 cm / 20 deg. Points below are read directly off the **Model**
panel -- the one computed over the full 0-130 deg range (the two
anatomical specimens' plotted ranges are narrower, limited by how far
each cadaver's tendon could be excursion-tested) and the one whose shape
is most directly analogous to what this project needs (a smooth
model-computed `r(theta)`, not a single specimen's noisy measurement).
**These are visual reads off a printed figure, not extracted from a
data table -- treat each point as accurate to roughly ±0.15-0.2 cm, not
exact:**

| theta (deg) | r (cm), digitized from Model panel |
|---|---|
| 0   | 0.6 |
| 20  | 1.0 |
| 40  | 1.4 |
| 60  | 1.9 |
| 80  | 2.3 |
| 100 | 2.7 |
| 120 | 3.0 |
| 130 | 3.05 (extrapolated plateau, consistent with the paper's text: "increases throughout the range of motion to a maximum of approximately 3 cm") |

Cross-check against the Male Specimen panel (anatomical, not
model-computed): visually consistent with this same shape and reaches
the same ~3 cm by 120 deg -- expected, since the paper reports ICC > 0.99
between the male specimen and model `BRA` curves specifically.

> **Revision note (2026-09-11): the table above is superseded.** A
> pixel-level re-digitization of the same Model panel at 600 dpi
> (`scripts/digitize_murray_fig4.py`, `data/murray1995_fig4_digitized.csv`,
> method in `docs/theory/phase7_limitations.md` section 1) shows this
> by-eye table reads low below 90 deg by more than its stated
> +/-0.15-0.2 cm: 0.89 vs 0.6 cm at 0 deg, 1.38 vs 1.0 at 20, 1.85 vs 1.4
> at 40, 2.23 vs 1.9 at 60; it agrees above 90 deg. The objective check
> is Murray's own Table 2 (model BRA spread 48% over 25-120 deg): the
> pixel curve gives 48.0%, this table about 63%. The code now uses the
> pixel curve (5-degree grid, peak 2.88 cm above 110 deg, range-of-motion
> average 2.18 cm, i.e. 21% above Holzbaur's 1.8 cm rather than the 8%
> reported below). Every downstream number was rerun.

**v1 implementation: piecewise-linear interpolation through this table**
(`numpy.interp` or equivalent), not a fitted closed-form function --
honest about the data being a small set of read-off points rather than a
curve with a known functional form, and trivial to replace with a better
digitization later if it turns out to matter.

Trapezoidal-rule average of the table over `[0, 130 deg]`: **~1.94 cm**,
within ~8% of Holzbaur (2005)'s independently-reported `ma_avg = 1.8 cm`
(CLAUDE.md 6.3) -- a real, non-circular cross-check, since this table was
read off Murray (1995)'s figure without reference to Holzbaur's number,
and the two sources agree to within plausible digitization error. This
replaces the previous linear-ramp draft's 11% average-value gap with an
~8% one, using actual curve shape instead of an assumed one.

**Remaining judgment call, now smaller:** the digitization above is a
visual read of a printed figure, not a precise data extraction -- if the
torque-angle validation check (section 6 below) turns up a discrepancy
that doesn't wash out within the ±0.2 cm/point reading uncertainty
already flagged, the next-best fix is extracting the actual
`arm26.osim` wrap-cylinder geometry (CLAUDE.md 6.2), not re-reading the
figure more carefully.

> **Revision note (2026-09-07):** for the manuscripts' validation
> figure, the **male-specimen** (anatomical, tendon-excursion) `BRA`
> curve was also digitized from the same Fig. 4 (dash-dot curve, Male
> Specimen panel), as an independent comparison against the model-panel
> digitization above: (15, 1.3), (25, 1.45), (40, 1.65), (60, 2.0),
> (80, 2.5), (100, 2.95), (120, 3.25) (deg, cm), +/-0.15-0.2 cm per
> point, measured range ~15-125 deg. Two anchors from Murray's own text
> and tables: the male max is "approximately 3 cm" (our read: 3.25 at
> the last measured point), and Table 2's (max-min)/max for this
> specimen is 58% (our reads give 60% -- consistent within the stated
> reading error). Murray reports ICC > 0.99 between this specimen curve
> and their model curve, which is what the figure now shows visually.
> Lives in `scripts/poster_figures.py` (`_MALE_DEG`/`_MALE_CM`).

`l^MT(theta)` follows from `r(theta) = -dl^MT/dtheta`, integrated
numerically (piecewise-linear `r` integrates to a piecewise-quadratic
`l^MT`, straightforward in code, not worth hand-deriving symbolically
here):

```
l^MT(theta) = l^MT(0) - integral_0^theta r(phi) dphi
```

`l^MT(0)` (muscle-tendon length at full extension) is not independently
known -- fixed by requiring the fiber to sit at its optimal length
`l_0^M` at some reference angle. Following the usual convention that
optimal fiber length corresponds to a mid-ROM posture (Holzbaur-style
models are typically built this way), set `l^M(65 deg) = l_0^M =
0.0858 m` (CLAUDE.md 6.2) and solve for `l^MT(0)` via the rigid-tendon
relation `l^M = (l^MT - l_s^T)/cos(alpha)`, `alpha=0`:

```
l^MT(65deg) = l_0^M + l_s^T = 0.0858 + 0.0535 = 0.1393 m
l^MT(0) = l^MT(65deg) + integral_0^65deg r(phi) dphi
```

(Numeric evaluation deferred to implementation -- this is an equation to
code, not a number to hand-derive here.)

## 4. Fiber velocity

Rigid tendon means muscle-tendon velocity is entirely fiber velocity
(pennation 0, so `cos(alpha)=1` and there's no pennation-angle-rate term
either):

```
V^M = dl^M/dt = dl^M/dtheta * thetadot = -r(theta) * thetadot
```

(negative because increasing `theta`, i.e. flexing, shortens the muscle
belly, i.e. `dl^MT/dtheta = -r(theta) < 0`.) This is the direct link
between joint angular velocity and the force-velocity relation in
CLAUDE.md 6.1 -- `V^M` here plugs into the inverted force-velocity
solve needed to get `F^MT` from `(theta, thetadot, a)`.

## 5. Equation of motion

```
I_elbow * thetaddot = tau_muscle(theta, thetadot, a) - m*g*l_c*sin(theta)
tau_muscle = F^MT(theta, thetadot, a) * r(theta)
```

No joint damping term added -- nothing in section 6 or here motivates a
specific value, and section 3's minimal-version discipline says not to
add unmotivated terms. Add one later only if simulation reveals a
concrete numerical need (e.g. chattering at rest).

Together with the activation ODE (CLAUDE.md 6.1) and `V^M` above, this
is the full state-space model:

```
d/dt [theta]    = [thetadot]
d/dt [thetadot] = [tau_muscle(theta,thetadot,a)/I_elbow - (m g l_c/I_elbow) sin(theta)]
d/dt [a]        = [(u-a)/tau_a(a,u)]
```

## 6. Validation plan (next step, not yet run)

Per CLAUDE.md section 4 step 1 and section 6.3: compute the **isometric**
torque-angle curve `tau(theta) = F^MT(theta; V^M=0, a=1) * r(theta)`
across `theta in [0, 130deg]` and check against the combined target in
6.3 (peak `2.0-3.5 cm` moment arm above ~100 deg, `ma_avg ~= 1.8 cm`,
monotonic non-mid-peaking shape). `r(theta)` is now a digitization of the
actual Murray (1995) curve (section 3 above), not a shape fit to summary
statistics, so the moment-arm part of this check is largely "did I
transcribe the table correctly," not a real independent test -- a
genuinely independent geometric validation would need a third source
Holzbaur (2005)'s own model wasn't built from (not pursued here, out of
scope for a 1-page abstract). The part of the check that **is**
independent, and the actual test of the muscle *force* model (section
6.1's equations, which the moment-arm digitization has no bearing on), is
whether `F^MT(theta;V^M=0,a=1)` (active + passive fiber force, no moment
arm involved yet) stays sensibly bounded and peaks near `l̄^M=1`, i.e.
near `theta=65 deg` by construction in section 3 -- that's the real
force-model check, separate from the moment-arm-shape check.

Not implemented in code yet -- this document is the checkpoint before
writing that code.
