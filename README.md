# Muscle Activation Control

Non-smooth equilibria in Hill-type muscle models and their quantified
consequence for linearization-based control, studied on a single elbow
joint (brachialis, with an optional triceps antagonist).

## Motivation

Controllers for muscle-actuated systems (prosthetics, exosuits,
neuromuscular robots) commonly use equilibrium-linearization-based design
(LQR, successive-linearization MPC), which implicitly assumes a single
well-defined local linear model at each operating point. This project
shows that Hill-type muscle models do not admit one.

Every equilibrium of a Hill-type single-muscle system is structurally
non-smooth, in two independent, literature-grounded ways:

1. **Activation dynamics** switches its time constant depending on
   whether excitation exceeds activation (Winters, 1995); any
   equilibrium requires excitation to exactly equal activation
   (`u* = a*`), sitting precisely on that switch.
2. **The force-velocity relation** switches between shortening and
   lengthening branches, with the lengthening-side slope built to be
   exactly twice the shortening-side slope at zero velocity (Katz,
   1939, as implemented in Thelen 2003); any equilibrium requires zero
   joint velocity, sitting precisely on that switch too.

We formalize this with Camlibel, Heemels, and Schumacher's (2008)
framework for bimodal piecewise-linear systems, verify both switches
satisfy their continuity condition exactly (residual `0.0`), and
quantify the practical consequence: a successive-linearization MPC-style
one-step prediction that uses the wrong branch accrues a joint-angle
error equal to **13-88% of the true predicted signal** across horizons
from 50-500 ms.

See [`docs/abstract.md`](docs/abstract.md) for the full write-up and
[`docs/theory/`](docs/theory) for the phase-by-phase derivations.

## Model

- Rigid-tendon Thelen (2003) Hill-type activation/contraction dynamics
- Holzbaur et al. (2005) brachialis muscle parameters
- Murray et al. (1995) digitized moment-arm data
- de Leva (1996) forearm+hand inertial parameters
- Optional lumped triceps antagonist for closed-loop / co-contraction
  studies

## Repository layout

```
src/muscle_activation_control/  Core library (muscle, joint, linearization,
                                 bimodal formalization, MPC, antagonist pair)
docs/theory/                    Phase-by-phase derivations and theory notes
docs/abstract.md                Workshop abstract draft with headline results
scripts/                        Analysis, validation, and figure-generation scripts
results/                        Generated figures, metrics, and headline numbers
tests/                          pytest unit tests
data/                           Digitized source data (Holzbaur 2005, Murray 1995)
```

## Setup

This project uses [uv](https://docs.astral.sh/uv/) for dependency
management.

```bash
uv sync
```

Or with plain pip:

```bash
pip install -e ".[dev]"
```

## Running

```bash
# Run the test suite
pytest

# Reproduce a specific phase's figures/results, e.g.
python scripts/validate_phase1.py
python scripts/phase6b_closed_loop.py
```

Each `scripts/phase*.py` script corresponds to a `docs/theory/phase*.md`
derivation and writes its output into `results/`.
