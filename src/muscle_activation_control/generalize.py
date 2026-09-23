"""Phase 7 section 3: other elbow muscles as parameterized Hill muscles.

docs/theory/phase7_limitations.md section 3. A MuscleSpec is one
musculotendon compartment (Holzbaur et al. 2005 Table 1 parameters, rigid
tendon, constant-width pennation) with a moment arm from the pixel-level
Murray et al. (1995) Fig. 4 digitization (data/murray1995_fig4_digitized.csv,
PCHIP interpolation, magnitude). A plant is a list of compartments that
share ONE activation state (e.g. the two biceps heads) acting on the
forearm+hand pendulum from joint.py, so single-muscle plants for biceps
and brachioradialis reuse exactly the brachialis plant's structure.

Conventions shared with joint.py: fibre at optimal length at
joint.THETA_REF (65 deg) for every compartment; flexion positive; a
flexor shortens with flexion, an extensor lengthens. Velocity branch names
are per compartment ("shortening"/"lengthening") but a plant is forced
with the joint-level label "flexing" / "extending" (one switch, thetadot).
"""

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from muscle_activation_control import joint, linearize, muscle

DATA = Path(__file__).resolve().parents[2] / "data" / "murray1995_fig4_digitized.csv"


def murray_moment_arm(name, panel="model"):
    """(theta_rad, |r| m) from the pixel-level digitization CSV."""
    rows = [r for r in csv.DictReader(open(DATA)) if r["panel"] == panel and r["muscle"] == name]
    deg = np.array([float(r["elbow_flexion_deg"]) for r in rows])
    cm = np.abs(np.array([float(r["moment_arm_cm"]) for r in rows]))
    return np.deg2rad(deg), cm * 1e-2


@dataclass
class MuscleSpec:
    name: str
    f0: float                 # N
    l0: float                 # m, optimal fibre length
    ls: float                 # m, tendon slack length
    alpha0_deg: float         # pennation at optimal length
    r_theta: np.ndarray       # rad
    r_m: np.ndarray           # m, moment-arm magnitude
    extensor: bool = False
    r_scale: float = 1.0
    _pchip: object = field(init=False, repr=False)

    def __post_init__(self):
        self._pchip = PchipInterpolator(self.r_theta, self.r_m)
        self.w = self.l0 * np.sin(np.deg2rad(self.alpha0_deg))
        self.l_mt_ref = self.ls + self.l0 * np.cos(np.deg2rad(self.alpha0_deg))
        g = np.linspace(joint.THETA_MIN, joint.THETA_MAX, 2001)
        rg = self.moment_arm(g)
        cum = np.concatenate(([0.0], np.cumsum(0.5 * (rg[1:] + rg[:-1]) * np.diff(g))))
        sgn = 1.0 if self.extensor else -1.0          # dl_MT/dtheta = sgn * r
        self._grid, self._dl = g, sgn * (cum - np.interp(joint.THETA_REF, g, cum))

    def moment_arm(self, theta):
        th = np.clip(np.asarray(theta, dtype=float), self.r_theta[0], self.r_theta[-1])
        return self.r_scale * self._pchip(th)

    def fibre(self, theta):
        """(l_bar, cos_alpha), constant-width rigid-tendon geometry."""
        proj = self.l_mt_ref + np.interp(np.asarray(theta, dtype=float), self._grid, self._dl) - self.ls
        l_m = np.sqrt(proj ** 2 + self.w ** 2)
        return l_m / self.l0, proj / l_m

    def torque(self, theta, thetadot, a, joint_vel=None):
        """Joint torque, flexion positive. joint_vel: None, "flexing", "extending"."""
        l_bar, cos_a = self.fibre(theta)
        r = self.moment_arm(theta)
        sgn = 1.0 if self.extensor else -1.0
        v_norm = sgn * r * np.asarray(thetadot, dtype=float) * cos_a / self.l0
        branch = None
        if joint_vel is not None:
            shortening = (joint_vel == "flexing") != self.extensor
            branch = "shortening" if shortening else "lengthening"
        f = muscle.muscle_tendon_force(a, l_bar, v_norm, self.f0, alpha=np.arccos(cos_a),
                                       velocity_branch=branch)
        return (-1.0 if self.extensor else 1.0) * f * r


def _table(name):
    return murray_moment_arm(name)


def bra_spec():
    """Identical to joint.py's brachialis (same table incl. cumulative max)."""
    return MuscleSpec("BRA", joint.F0_M, joint.L0_M, joint.LS_T, 0.0, joint._R_THETA_RAD, joint._R_THETA_M)


def bic_specs():
    th, r = _table("BIC")
    return [MuscleSpec("BIClong", 624.3, 0.116, 0.272, 0.0, th, r),
            MuscleSpec("BICshort", 435.6, 0.132, 0.192, 0.0, th, r)]


def brd_spec():
    th, r = _table("BRD")
    return MuscleSpec("BRD", 261.3, 0.173, 0.133, 0.0, th, r)


def tri_specs(scale=1.0):
    th, r = _table("TRI")
    return [MuscleSpec("TRIlong", 798.5, 0.134, 0.143, 12.0, th, r, extensor=True, r_scale=scale),
            MuscleSpec("TRIlat", 624.3, 0.114, 0.098, 9.0, th, r, extensor=True, r_scale=scale),
            MuscleSpec("TRImed", 624.3, 0.114, 0.091, 9.0, th, r, extensor=True, r_scale=scale)]


PLANTS = {"BRA": lambda: [bra_spec()], "BIC": bic_specs, "BRD": lambda: [brd_spec()]}


class Plant:
    """Compartments sharing one activation, driving the forearm+hand pendulum.
    State [theta, thetadot, a], input u; optional constant load tau_ext."""

    def __init__(self, specs, tau_ext=0.0):
        self.specs = specs
        self.tau_ext = tau_ext

    def muscle_torque(self, theta, thetadot, a, joint_vel=None):
        return sum(s.torque(theta, thetadot, a, joint_vel) for s in self.specs)

    def f(self, x, u, act=None, joint_vel=None):
        th, thd, a = x
        tau = self.muscle_torque(th, thd, a, joint_vel) - joint.gravity_torque(th) - self.tau_ext
        return np.array([thd, tau / joint.FOREARM_HAND_I_ELBOW,
                         muscle.activation_derivative(a, u, branch=act)])

    def equilibrium(self, theta):
        res = lambda a: self.muscle_torque(theta, 0.0, a) - joint.gravity_torque(theta) - self.tau_ext
        return brentq(res, 0.0, 1.0)

    def linearize(self, x, u, act, joint_vel, h=1e-6):
        x = np.asarray(x, dtype=float)
        f = lambda xx, uu: self.f(xx, uu, act, joint_vel)
        A = np.zeros((3, 3))
        for i in range(3):
            dx = np.zeros(3); dx[i] = h
            A[:, i] = (f(x + dx, u) - f(x - dx, u)) / (2 * h)
        B = ((f(x, u + h) - f(x, u - h)) / (2 * h)).reshape(3, 1)
        return A, B

    def error_fractions(self, theta, delta_u, n_steps, t_s=0.010):
        """Phase 4 protocol: final-angle error of each branch pair as a
        fraction of the true swing, for a held excitation step."""
        a0 = self.equilibrium(theta)
        x0 = np.array([theta, 0.0, a0])
        u = a0 + delta_u
        sol = solve_ivp(lambda t, x: self.f(x, u), (0, n_steps * t_s), x0,
                        t_eval=np.arange(1, n_steps + 1) * t_s, rtol=1e-10, atol=1e-12)
        th_true = sol.y[0, -1]
        swing = abs(th_true - theta)
        sgn = 1.0 if delta_u > 0 else -1.0
        occupancy = (float(np.mean(sgn * (u - sol.y[2]) >= -1e-9)),
                     float(np.mean(sgn * sol.y[1] >= -1e-9)))
        act_m = "activating" if delta_u > 0 else "deactivating"
        vel_m = "flexing" if delta_u > 0 else "extending"
        flip = {"activating": "deactivating", "deactivating": "activating",
                "flexing": "extending", "extending": "flexing"}
        cases = {"matched": (act_m, vel_m), "wrong_activation": (flip[act_m], vel_m),
                 "wrong_velocity": (act_m, flip[vel_m]), "both_wrong": (flip[act_m], flip[vel_m])}
        out = {"a_star": a0, "swing_deg": float(np.rad2deg(swing)), "occupancy": occupancy}
        for k, (act, jv) in cases.items():
            A, B = self.linearize(x0, a0, act, jv)
            Ad, Bd = linearize.discretize_zoh(A, B, t_s)
            d = np.zeros(3)
            for _ in range(n_steps):
                d = Ad @ d + Bd[:, 0] * delta_u
            out[k] = abs(theta + d[0] - th_true) / swing
        return out


def group_isometric_torque(specs, theta, a=1.0):
    return sum(s.torque(theta, 0.0, a) for s in specs)
