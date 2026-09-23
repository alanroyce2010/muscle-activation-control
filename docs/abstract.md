# Non-Smooth Equilibria in Hill-Type Muscle Models: A Quantified Consequence for Linearization-Based Control

*Workshop abstract draft, reframed 2026-07-31 (second pass) around the
project's actual finding rather than its original motivating hypothesis.
See the revision note in CLAUDE.md section 1 for why. Content drawn
directly from CLAUDE.md sections 4.1.1, 4.2.1, 4.2.2 and the underlying
theory docs in `docs/theory/`; all numbers below are taken from those
documents, not re-derived here. Convert to the target venue's template as
needed; this is the content draft, not a formatted camera-ready file.
The template-fitted, submission-ready version is
`docs/abstract_IROS2026_neuromuscular_robotics.docx`.*

## Abstract

Controllers for muscle-actuated systems (prosthetics, exosuits,
neuromuscular robots) routinely rely on equilibrium-linearization-based
design (LQR, successive-linearization MPC), which implicitly assumes a
single well-defined local linear model at each operating point. Whether
Hill-type muscle models actually admit one has not been examined. This is
distinct from Yeo, Verheul, Herzog, and Sueda's (2023) recent finding of
Hill-type numerical instability, which is a different phenomenon
(negative-stiffness eigenvalue instability arising from the descending
limb of the force-length curve, plus history-dependence effects the
model's state-space assumption misses) and does not address activation
dynamics or its consequence for equilibrium-linearization-based control
at all. Using a validated single-joint model (brachialis, rigid-tendon
Thelen (2003) Hill dynamics, Holzbaur et al. (2005) muscle parameters,
Murray et al. (1995) moment-arm data, de Leva (1996) forearm+hand
inertial parameters), we show Hill-type models do not admit a single
well-defined linearization at equilibrium.

**Every equilibrium of a Hill-type single-muscle system is structurally
non-smooth**, in two independent, literature-grounded ways. First,
activation dynamics switches its time constant depending on whether
excitation exceeds activation (Winters, 1995); any equilibrium requires
excitation to exactly equal activation (`u*=a*`), sitting precisely on
that switch. Second, the force-velocity relation switches between
shortening and lengthening branches with the lengthening-side slope built
to be exactly twice the shortening-side slope at zero velocity (Katz,
1939, as implemented in Thelen 2003); any equilibrium requires zero joint
velocity, sitting precisely on that switch too. Both are general facts
about any finite-time-constant Hill-type activation model at any of its
equilibria, not artifacts of this project's specific joint or muscle
choice.

We formalize this using Camlibel, Heemels, and Schumacher's (2008)
framework for bimodal piecewise-linear systems with continuous vector
fields. Each switch, checked individually, satisfies their continuity
condition (`A₁-A₂=ec^T`, `b₁-b₂=ed`) *exactly*: verified algebraically,
not approximately (residual `0.0`). The two switches' effects further
occupy disjoint entries of the linearized system matrix (one touches only
the activation state's own row, the other only the joint-velocity state's
row), so the full four-mode local system decomposes exactly as the direct
sum of two independent bimodal switches, a property of this system rather
than a general claim requiring the fully multimodal ("conewise")
extension of their theory.

Does this non-smoothness matter, or is it a curiosity? We quantify it
directly: for a representative successive-linearization MPC-style
one-step prediction (`T_s=10 ms`, small-signal excitation perturbation,
deliberately the regime in which linearization is normally trusted), a
prediction model that silently uses the *wrong* branch accrues a
joint-angle error equal to **13-88% of the true predicted signal**,
across horizons from 50-500 ms, and 1-2 orders of magnitude larger than
the correctly-matched branch's own error at every horizon tested.
**The answer is yes, it matters.**

## Key takeaways

1. **Hill-type equilibria are provably non-smooth, not just observed to
   be.** Both switches individually satisfy Camlibel et al.'s (2008)
   continuity condition exactly (residual `0.0`), and their effects are
   structurally independent (disjoint matrix entries): a verified
   property of this system, not an assumption or a numerical artifact.
2. **The non-smoothness is practically significant, not a curiosity.** A
   controller built on standard equilibrium-linearization machinery that
   does not track which side of each switch the operating point is on
   will carry a quantifiably large model-plant mismatch (13-88% of the
   true predicted signal), even in the small-signal regime such
   linearizations are normally trusted in.

*Figure: `results/phase4_activating_step_small.png`. True nonlinear
joint-angle trajectory vs. the four branch-linearizations' one-step-ahead
predictions, small-signal excitation step at the equilibrium. The
matched branch tracks the truth almost exactly; the other three visibly
diverge within 50-100 ms.*

## Post-submission addendum (2026-08-24, not part of the submitted text)

The two items listed as future work in the submitted abstract have been
done (CLAUDE.md section 4.3; `docs/theory/phase6{a,b,c}_*.md`), plus a
sweep over all equilibria. For the poster or a follow-on paper:

- **Every equilibrium**: the activation-branch mismatch is closed-form in
  `a*` alone, `tau_deact/(tau_act(0.5+1.5a*)^2)`, from 13.3x at rest to
  1x at `a*=0.884`; the force-velocity mismatch is exactly 2x everywhere
  for a single muscle. The combined wrong-branch prediction error is
  never below 20% of the true swing at any equilibrium tested
  (`a*=0.02-0.95`), and 3-4x larger for deactivating than activating
  steps.
- **Closed loop**: at 10 ms with successive-linearization MPC, feedback
  removes the stability and steady-state consequences but not the
  transient ones: a central-difference (blended) linearization gives 8x
  the step overshoot of a branch-matched one, and a naive fixed branch
  38% higher periodic-tracking RMS; a sign-test branch rule with a
  deadband recovers most of it. Holding one branch over the horizon is
  the remaining error, which calls for hybrid MPC.
- **Antagonist pair**: the structure survives (3 switches, 8 modes, still
  a direct sum of bimodal switches), and the force-velocity kink becomes
  `(2D_b+D_t)/(D_b+2D_t)`, which co-contraction can drive to 1: at
  `a_t*=0.16` the joint's velocity dynamics are locally smooth. A
  possible functional role for co-contraction that the smooth-model view
  cannot express.

## References

- Thelen, D.G. (2003). Adjustment of muscle mechanics model parameters to
  simulate dynamic contractions in older adults. *J Biomech Eng*
  125(1):70-77.
- Winters, J.M. (1995). An improved muscle-reflex actuator for use in
  large-scale neuromusculoskeletal models. *Ann Biomed Eng* 23:359-374.
- Katz, B. (1939). The relation between force and speed in muscular
  contraction. *J Physiol* 96:45-64.
- Holzbaur, K.R.S., Murray, W.M., Delp, S.L. (2005). A model of the upper
  extremity for simulating musculoskeletal surgery and analyzing
  neuromuscular control. *Ann Biomed Eng* 33(6):829-840.
- Murray, W.M., Delp, S.L., Buchanan, T.S. (1995). Variation of muscle
  moment arms with elbow and forearm position. *J Biomech* 28(5):513-525.
- de Leva, P. (1996). Adjustments to Zatsiorsky-Seluyanov's segment
  inertia parameters. *J Biomech* 29(9):1223-1230.
- Camlibel, M.K., Heemels, W.P.M.H., Schumacher, J.M. (2008). A full
  characterization of stabilizability of bimodal piecewise linear systems
  with scalar inputs. *Automatica* 44(5):1261-1267.
- Yeo, S-H., Verheul, J., Herzog, W., Sueda, S. (2023). Numerical
  instability of Hill-type muscle models. *J R Soc Interface*
  20(199):20220430.
- Millard, M., Uchida, T., Seth, A., Delp, S.L. (2013). Flexing
  computational muscle: modeling and simulation of musculotendon
  dynamics. *J Biomech Eng* 135(2):021005.
- Zajac, F.E. (1989). Muscle and tendon: properties, models, scaling, and
  application to biomechanics and motor control. *Crit Rev Biomed Eng*
  17(4):359-411.
