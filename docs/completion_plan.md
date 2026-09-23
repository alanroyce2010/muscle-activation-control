# Completion plan (2026-08-24)

Assessment of where the project stands against CLAUDE.md section 9, and
the ordered work needed to call it complete. Two tiers: Tier 1 closes out
the abstract deliverable as it stands; Tier 2 is the work the abstract
itself promised as future work, needed if the poster is accepted.

## 0. State of the project

Done and verified (2026-08-24 check):

- Phases 0-5 complete. 35/35 tests pass on the remote
  (`abhiyaan-cu@cu:/home/abhiyaan-cu/Alan/muscle-activation-control`).
- Abstract exists in three forms: `docs/abstract.md` (content draft),
  `docs/abstract_IROS2026_neuromuscular_robotics.docx` (template-filled,
  author/affiliation/email now filled in), and
  `docs/POSTER_AlanRoyceGabriel_Samuel.pdf` (rendered 1-page check copy,
  dated 2026-08-02).
- Related-work check (CLAUDE.md 4.2.3) done; Yeo et al. 2023 cited.

Loose ends found:

1. **No git history.** `main` has zero commits; every file is untracked.
   Nothing is versioned, including the submitted abstract.
2. **Local environment not runnable.** `python3` here has no scipy /
   matplotlib and the package is not installed, so `pytest` fails at
   collection locally (5 import errors). Only the remote works.
3. **`pyproject.toml` is stale**: description still states the pre-pivot
   hypothesis; dependencies omit `sympy`
   (`scripts/_verify_branch_derivatives.py`), `python-docx`
   (`scripts/fill_abstract_template.py`), and `pytest`.
4. **Doc drift in CLAUDE.md**: section 6.4 does not list
   `deleva1996_segment_inertia_parameters.pdf` although
   `single_joint_dynamics.md` and the abstract cite it; section 9 refers
   to `docs/abstract_IROS2026_neuromuscular_robotics.pdf` (does not
   exist; the rendered copy is `POSTER_AlanRoyceGabriel_Samuel.pdf`) and
   says author fields were left as placeholders (they are now filled);
   section 9 does not record the 4.2.3 related-work check as a milestone.
5. **`docs/theory/single_joint_dynamics.md` section 6** still reads
   "Validation plan (next step, not yet run)" although Phase 1 validation
   ran and passed.
6. **The abstract commits to future work** the repo has not started:
   "antagonist pair, closed-loop MPC test of branch-tracking against a
   naive fixed linearization" (Conclusions section of the poster).
7. Phase 4's quantification is at **one** equilibrium (`theta*=65 deg`,
   `a*=0.133`), while the abstract's claim is about **every** equilibrium.
   The structural argument covers that, but the quantitative one does not
   yet.

## Tier 3 (added 2026-08-26): poster accepted, poster + flash pitch built

Poster #5, 1st Workshop on Neuromuscular Robotics, IROS 2026, 1 October
2026, Pittsburgh. `docs/poster/`: `poster_A0.pdf` (portrait A0),
`05_AlanRoyceGabriel_Samuel.pptx` (flash pitch, **due 21 September**),
`flash_pitch.md`. Details in CLAUDE.md section 9, Phase 5.

## Tier 4 (2026-08-26): journal manuscript

Target: J R Soc Interface (`docs/journal_options.md`). Draft compiled:
`docs/manuscript/main.pdf` (11 pages, rsproca_new template, tectonic).
Open: placeholders (DOI, funding, acknowledgements) and the section 3
pre-submission list in `journal_options.md`. Tier 1 items T1.1 (commit
and tag) and T1.2 (`pyproject.toml`) are now done; T1.3 doc drift and
T1.5 conewise citation remain.

## Tier 1: close out the deliverable (housekeeping, ~half a day)

Order matters: version first, then fix, then re-verify.

- [ ] T1.1 `git init` is already done; make the first commit of the
      current tree (add `.gitignore` for `__pycache__`, `.pytest_cache`,
      `.DS_Store`, `*.egg-info`). Tag it `abstract-submitted-2026-08-02`
      so the submitted state is recoverable. (Needs the user's go-ahead
      to commit.)
- [ ] T1.2 Make the local checkout runnable: `pip install -e .` plus the
      missing deps, so `pytest` passes locally as well as on the remote.
      Update `pyproject.toml`: add `sympy`, `pytest`; add `python-docx`
      under an optional `[project.optional-dependencies] abstract` group;
      rewrite `description` to the pivoted question.
- [ ] T1.3 Sync docs to reality (all in-place with `> **Revision note**`
      callouts, per the section 2 convention, `--` not em-dash):
      - CLAUDE.md 6.4: add de Leva (1996) entry.
      - CLAUDE.md 9 / Phase 5: point at `POSTER_AlanRoyceGabriel_Samuel.pdf`,
        note author fields are filled, add 4.2.3 as a checked milestone.
      - `single_joint_dynamics.md` section 6: mark as run, link
        `scripts/validate_phase1.py` and `results/phase1_torque_angle_validation.png`.
- [ ] T1.4 Re-run `pytest` locally and on the remote after T1.2; rsync
      (no `--delete`) so the remote mirrors the committed tree.
- [ ] T1.5 Verify the one unverified citation flagged in CLAUDE.md 4.2:
      Camlibel, Heemels, Schumacher, "Algebraic necessary and sufficient
      conditions for the controllability of conewise linear systems,"
      IEEE TAC (believed 53(3):762-774, 2008 -- confirm from the primary
      source before it is ever cited). Not needed for the abstract; needed
      before Tier 2b cites it.

## Tier 2: Phase 6 -- the work the abstract promised

> **Status (2026-08-24): T2.a, T2.b, T2.c done.** Results in CLAUDE.md
> section 4.3 and `docs/theory/phase6{a,b,c}_*.md`; code in
> `src/muscle_activation_control/{prediction,mpc,antagonist,prediction_pair}.py`;
> 30 new tests (65/65 total, run on the remote). Deviations from the plan
> below, all recorded in the theory docs: T2.b grew a fifth variant
> (per-step branch assignment) and two bug-fix lessons; T2.c's worry
> about overlapping force-velocity `e` vectors resolved itself (they
> share one switching function, so it is one switch, not two). T2.d
> remains deprioritized. Tier 1 items T1.1 (commit/tag), T1.3 (doc
> drift in section 9's Phase 5 entry, `single_joint_dynamics.md`
> section 6) and T1.5 (conewise citation) are still open; T1.2's
> `pyproject.toml` fix was done as a side effect (cvxpy needed).

Follow section 7's discipline for each item: derive (theory doc), then
implement, then validate, with a review stop between. Ordered by
value-per-effort and by how directly each item defends the submitted
claim.

### T2.a Quantify across all equilibria (cheap, defends the "every equilibrium" claim)

- Derive: the activation-branch ratio is closed-form,
  `tau_deact / (tau_act * (0.5 + 1.5 a*)^2)`, so it ranges from 13.3x at
  `a*=0` to 0.83x at `a*=1` (crossing 1x near `a*=0.88`, where the two
  branches momentarily agree). The FV ratio is a constant 2x. Write this
  up as `docs/theory/phase6a_equilibrium_sweep.md`.
- Implement: sweep the isometric equilibrium over the ROM
  (`theta*` in 10-125 deg, which sets `a*` through gravity balance), repeat
  the Phase 4 one-step prediction comparison at each, and plot the
  wrong-branch error fraction vs `theta*` / `a*`.
- Validate: the `theta*=65 deg` point must reproduce the 13-88% figure.
- Deliverable: one figure and one paragraph; strengthens the poster.
- Risk: at high `a*` the activation branches nearly agree, so the
  wrong-branch error may shrink there -- report it either way; it
  sharpens the claim ("matters most at low activation, i.e. exactly the
  postural regime prosthetic/exosuit controllers spend most time in").

### T2.b Closed-loop test: branch-tracking MPC vs naive fixed-branch MPC

This is the open scientific question. Phase 4 is an open-loop prediction
result; feedback may or may not wash the mismatch out. The outcome could
weaken or strengthen the abstract's practical claim, and must be reported
faithfully either way.

- Derive: `docs/theory/phase6b_closed_loop.md`. Successive-linearization
  MPC over the full state `[theta, thetadot, a]` with input `u` (no
  reduced model -- Phase 2 showed there is nothing to reduce). Two
  prediction-model variants: (i) *branch-tracking*, select the
  linearization each step from `sign(u_k - a_k)` and `sign(thetadot_k)`;
  (ii) *naive fixed*, one branch chosen once at the equilibrium. Also
  (iii) *central-difference blended*, the "silent" failure mode
  `linearize.py`'s docstring warns about. Task: small setpoint steps
  (+/- 2-5 deg) and a slow sinusoid around `theta*=65 deg`, `T_s=10 ms`,
  horizon 200-500 ms, constraints `0<=u<=1`.
- Implement: `src/muscle_activation_control/mpc.py`, reusing the
  `MPCController` shape from `origami-arm-control` (cvxpy QP,
  `cont2discrete(..., "zoh")`, equilibrium feedforward). Add `cvxpy` to
  dependencies. Run on the remote.
- Validate: closed-loop RMS tracking error, settling time, and
  per-step one-step prediction residual for each variant. Sanity check:
  the branch-tracking variant's residual should match Phase 4's matched
  branch error.
- Deliverable: `docs/theory/phase6b_closed_loop.md`,
  `results/phase6b_*.png`, tests in `tests/test_mpc.py`.

### T2.c Antagonist pair (brachialis + triceps)

- Derive: add triceps (`TRIlong` or lumped triceps from `arm26.osim`,
  Holzbaur 2005 Table 1) with its own activation state. State becomes
  `[theta, thetadot, a_bra, a_tri]`, inputs `[u_bra, u_tri]`. Each
  equilibrium now sits on 4 switches (2 activation, 2 force-velocity), up
  to 16 modes. The key structural check: are the four `e` vectors still
  of disjoint support (so the system stays a direct sum of bimodal
  switches, and the 2008 Automatica paper still suffices)? Expected yes
  (each activation `e` touches only its own `a_i` row; both FV `e`s touch
  only the `thetadot` row -- **those two overlap**, so the FV pair may
  *not* be a direct sum, which would be the first genuinely conewise
  case and would need the T1.5 citation).
- Implement: generalize `joint.py`, `linearize.py`, `bimodal.py` to N
  muscles. Co-contraction level becomes a free equilibrium parameter;
  sweep it.
- Validate: extend `test_bimodal.py` for the disjoint-support check;
  Phase 1 torque-angle validation for triceps against Murray 1995.
- Deliverable: `docs/theory/phase6c_antagonist_pair.md`, plus tests.
- Also re-tests the `origami-arm-control` finding that model-free RL
  fails at multi-actuator coordination, if the RL comparison (T2.d) is
  ever pursued.

### T2.d Deprioritized (only if the venue asks or time allows)

- PID / LQR / MPC / RL controller comparison (original section 5).
- Original Phase 3 algebraic-substitution reduced model (no stiffness
  problem to fix; Phase 2 answered this).
- Compliant tendon, tendon hysteresis (Bouc-Wen).

## Suggested sequence

1. Tier 1 in full (T1.1-T1.5), one commit per item.
2. T2.a (a day; directly defends the submitted claim).
3. T2.b (2-3 days; the actual open question).
4. T2.c (2-3 days) if the poster is accepted and there is room for a
   second figure or a follow-on paper.
5. Update CLAUDE.md section 9 after each phase, same as before.
