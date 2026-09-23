"""Phase 6b: successive-linearization MPC on the single-joint Hill model,
with the branch-selection rule at the (non-smooth) reference equilibrium
as the experimental variable.

See docs/theory/phase6b_closed_loop.md. Mirrors the QP structure of
origami-arm-control's MPCController (cvxpy, ZOH via cont2discrete,
equilibrium feedforward, DARE terminal cost) with no inner loop and the
full state [theta, thetadot, a] predicted directly.
"""

import numpy as np
import cvxpy as cp
from scipy.integrate import solve_ivp
from scipy.linalg import solve_discrete_are, sqrtm

from muscle_activation_control import joint, linearize, muscle

MODES = ("fixed_naive", "blended", "tracking", "current_state", "trajectory", "hybrid", "blended_act")
NAIVE_BRANCH = ("deactivating", "shortening")  # what `u > a` / `v <= 0` return at equality


class BranchAwareMPC:
    def __init__(self, mode, t_s=0.010, horizon=20, q_theta=1000.0, q_thetadot=1.0,
                 r=200.0, u_min=0.0, u_max=1.0, act_deadband=2e-3, vel_deadband=0.02,
                 ref_deadband=np.deg2rad(0.05), fixed_branch=None, traj_max_iter=4,
                 switch_grid=None):
        """mode: one of MODES, or "fixed" with an explicit fixed_branch=
        (activation_branch, velocity_branch) held for every step
        (diagnostic: lets each of the 4 branches be run as a fixed
        linearization, not just the naive one)."""
        if mode == "fixed":
            if fixed_branch is None:
                raise ValueError("mode='fixed' needs fixed_branch")
            self.fixed_branch = tuple(fixed_branch)
        elif mode == "fixed_naive":
            self.fixed_branch = NAIVE_BRANCH
        elif mode not in MODES:
            raise ValueError(mode)
        self.mode = mode
        self.t_s = t_s
        self.N = horizon
        self.Q = np.diag([q_theta, q_thetadot, 0.0])
        self.R = float(r)
        self.u_min, self.u_max = u_min, u_max
        # Deadbands for the sign-test branch rule. With a strict sign test
        # (deadbands -> 0) the rule chatters at rest: at the equilibrium
        # u-a and thetadot hover around zero at integrator/QP precision and
        # the chosen model flips every step, which shows up as a +/-0.01
        # sawtooth in u (first Phase 6b run, recorded in
        # docs/theory/phase6b_closed_loop.md section 5). Inside the
        # deadband the direction the reference asks for decides; if that
        # is also negligible the previous choice is kept.
        self.act_deadband = act_deadband
        self.vel_deadband = vel_deadband
        self.ref_deadband = ref_deadband
        self._last_branch = ("activating", "shortening")
        self._eq_cache = {}
        # 'trajectory' mode: re-solve until the per-step branch assignment
        # read off the plan agrees with the one the plan was built with
        # (successive piecewise linearization). Without this loop the
        # assignment flips between consecutive solves and u chatters
        # (second Phase 6b run, docs/theory/phase6b_closed_loop.md s.5).
        self.traj_max_iter = traj_max_iter
        self.traj_iters = []          # iterations used per control step
        self.traj_fallbacks = 0       # steps where it cycled and fell back
        # 'hybrid' mode (docs/theory/phase7_limitations.md section 6): switch
        # times per switching signal, N meaning "never switch in the horizon"
        self.switch_grid = (sorted(set([0, 1, 2, 3, 5, 7, 10, 14, horizon]))
                            if switch_grid is None else list(switch_grid))
        self.hyb_violations = []      # consistency violations of the applied plan
        self._build_problem()

    # ---- QP (built once, parameters updated per step: DPP-compliant) ----
    def _build_problem(self):
        """Time-varying affine prediction model over the horizon
        (per-step parameters), so a mode may assign a different branch to
        each horizon step; time-invariant modes just repeat one model."""
        N = self.N
        self.p_A = [cp.Parameter((3, 3)) for _ in range(N)]
        self.p_B = [cp.Parameter(3) for _ in range(N)]
        self.p_d = [cp.Parameter(3) for _ in range(N)]     # affine term folded in
        self.p_x0 = cp.Parameter(3)
        self.p_xref = cp.Parameter((N + 1, 3))
        self.p_ulin = cp.Parameter()
        self.p_Pf_sqrt = cp.Parameter((3, 3))
        self.p_Pf_xref = cp.Parameter(3)   # = Pf_sqrt @ xref[N], pre-multiplied (DPP: no param*param)
        self.x = cp.Variable((N + 1, 3))
        self.u = cp.Variable(N)
        q_sqrt = np.sqrt(np.diag(self.Q))
        cost = 0
        cons = [self.x[0] == self.p_x0]
        for j in range(N):
            cost += cp.sum_squares(cp.multiply(q_sqrt, self.x[j] - self.p_xref[j]))
            cost += self.R * cp.square(self.u[j] - self.p_ulin)
            cons.append(self.x[j + 1] == self.p_A[j] @ self.x[j] + self.p_B[j] * self.u[j] + self.p_d[j])
            cons += [self.u[j] >= self.u_min, self.u[j] <= self.u_max]
        cost += cp.sum_squares(self.p_Pf_sqrt @ self.x[N] - self.p_Pf_xref)
        self.problem = cp.Problem(cp.Minimize(cost), cons)
        assert self.problem.is_dcp(dpp=True)
        self._plan = None  # last (x_plan (N+1,3), u_plan (N,)) for 'trajectory' mode

    # ---- model selection (the experiment) ----
    def _a_ref(self, theta_ref):
        key = round(float(theta_ref), 12)
        if key not in self._eq_cache:
            self._eq_cache[key] = linearize.find_equilibrium(theta_ref)
        return self._eq_cache[key]

    def select_branch(self, x, u_prev, theta_ref, last=None):
        """Branch rule for 'tracking', 'current_state' and (per horizon
        step) 'trajectory': the side of each switch the operating point
        (x, u_prev) is on, with ties broken by the direction the reference
        asks for, and a full tie keeping the previous choice `last`."""
        last = self._last_branch if last is None else last
        theta, thetadot, a = x
        want_flex = theta_ref - theta
        y_act = u_prev - a
        if abs(y_act) > self.act_deadband:
            act = "activating" if y_act > 0 else "deactivating"
        elif abs(want_flex) > self.ref_deadband:
            act = "activating" if want_flex > 0 else "deactivating"
        else:
            act = last[0]
        if abs(thetadot) > self.vel_deadband:
            vel = "shortening" if thetadot > 0 else "lengthening"
        elif abs(want_flex) > self.ref_deadband:
            vel = "shortening" if want_flex > 0 else "lengthening"
        else:
            vel = last[1]
        self._last_branch = (act, vel)
        return act, vel

    def prediction_model(self, x, u_prev, theta_ref_preview):
        """Return (x_lin, u_lin, models, branches) where models is a list of
        N (A_d, B_d, c_d) tuples (one per horizon step) and branches the
        branch used at each step."""
        theta_ref = float(theta_ref_preview[0])
        a_ref = self._a_ref(theta_ref)
        x_eq = np.array([theta_ref, 0.0, a_ref])
        f0_eq = joint.state_derivative(x_eq, a_ref)

        def eq_model(branch):
            if None in branch:   # an unforced switch comes out averaged
                a_c, b_c, f0 = linearize.linearize_general(x_eq, a_ref, *branch)
            else:
                a_c, b_c = linearize.linearize(theta_ref, a_ref, *branch)
                f0 = f0_eq
            return linearize.discretize_affine(a_c, b_c, f0, self.t_s)

        if self.mode in ("fixed_naive", "fixed"):
            branches = [self.fixed_branch] * self.N
            x_lin, u_lin = x_eq, a_ref
            models = [eq_model(self.fixed_branch)] * self.N
        elif self.mode == "blended":
            branches = [(None, None)] * self.N
            x_lin, u_lin = x_eq, a_ref
            models = [eq_model((None, None))] * self.N
        elif self.mode == "tracking":
            br = self.select_branch(x, u_prev, theta_ref)
            branches = [br] * self.N
            x_lin, u_lin = x_eq, a_ref
            models = [eq_model(br)] * self.N
        elif self.mode == "current_state":
            br = self.select_branch(x, u_prev, theta_ref)
            branches = [br] * self.N
            a_c, b_c, f0 = linearize.linearize_general(x, u_prev, *br)
            x_lin, u_lin = np.asarray(x, dtype=float), float(u_prev)
            models = [linearize.discretize_affine(a_c, b_c, f0, self.t_s)] * self.N
        elif self.mode == "blended_act":
            # the sign test with only the activation row averaged: the
            # Jacobian a tanh-smoothed activation model (De Groote et al.
            # 2016) has at any equilibrium, for any smoothing constant
            # (docs/theory/phase7_limitations.md section 10.3)
            br = (None, self.select_branch(x, u_prev, theta_ref)[1])
            branches = [br] * self.N
            x_lin, u_lin = x_eq, a_ref
            models = [eq_model(br)] * self.N
        else:  # trajectory: branch per horizon step from the previous plan, shifted
            x_lin, u_lin = x_eq, a_ref
            branches = self._branches_from_plan(x, u_prev, theta_ref_preview, self._plan, shift=True)
            models = self._models_for(branches, eq_model)
        self._eq_model = eq_model
        return x_lin, u_lin, models, branches

    def _solve(self, x_lin, u_lin, models, xref):
        for j, (a_d, b_d, c_d) in enumerate(models):
            self.p_A[j].value = a_d
            self.p_B[j].value = b_d[:, 0]
            self.p_d[j].value = x_lin - a_d @ x_lin - b_d[:, 0] * u_lin + c_d
        a_N, b_N, _ = models[-1]
        pf = solve_discrete_are(a_N, b_N, self.Q, np.array([[self.R]]))
        pf_sqrt = np.real(sqrtm(pf))
        self.p_Pf_sqrt.value = pf_sqrt
        self.p_Pf_xref.value = pf_sqrt @ xref[self.N]
        try:
            self.problem.solve(solver=cp.CLARABEL)
        except cp.error.SolverError:
            self.problem.solve(solver=cp.SCS)
        if self.u.value is None:
            raise RuntimeError(f"MPC QP did not solve (status={self.problem.status})")

    def _models_for(self, branches, eq_model):
        cache, models = {}, []
        for br in branches:
            if br not in cache:
                cache[br] = eq_model(br)
            models.append(cache[br])
        return models

    def _branches_from_plan(self, x, u_prev, theta_ref_preview, plan, shift):
        """Per-horizon-step branch assignment read off a planned (x, u)
        sequence. shift=True treats `plan` as last step's plan (advance it
        one step, measured x replaces its first entry); shift=False reads a
        plan just computed for this step."""
        if plan is None:
            br = self.select_branch(x, u_prev, float(theta_ref_preview[0]))
            return [br] * self.N
        x_plan, u_plan = plan
        if shift:
            x_plan = np.vstack([x_plan[1:], x_plan[-1:]])
            u_plan = np.concatenate([u_plan[1:], u_plan[-1:]])
        x_plan = x_plan.copy()
        x_plan[0] = x
        u_cmp = np.concatenate([[u_prev], u_plan[:-1]])  # input in force entering step j
        branches, last = [], self._last_branch
        for j in range(self.N):
            last = self.select_branch(x_plan[j], u_cmp[j], float(theta_ref_preview[j]), last)
            branches.append(last)
        self._last_branch = branches[0]
        return branches

    # ---- control ----
    def compute_control(self, x, u_prev, theta_ref_preview):
        """theta_ref_preview: array of N+1 reference angles from t_k on.
        Returns (u0, info) with info['x_pred_next'] the controller's own
        one-step prediction under u0 (for the residual metric)."""
        theta_ref_preview = np.asarray(theta_ref_preview, dtype=float)
        x = np.asarray(x, dtype=float)
        if self.mode == "hybrid":
            return self._compute_hybrid(x, u_prev, theta_ref_preview)
        x_lin, u_lin, models, branches = self.prediction_model(x, u_prev, theta_ref_preview)
        xref = np.zeros((self.N + 1, 3))
        xref[:, 0] = theta_ref_preview[: self.N + 1]
        # Reference is the FULL equilibrium state [theta_ref, 0, a_ref] at
        # every preview point. Q has zero weight on a, but the terminal
        # DARE cost couples a with theta, so centring it on a=0 instead of
        # a_ref biases the closed-loop fixed point (first Phase 6b run:
        # a spurious, branch-dependent steady-state offset of up to 0.26
        # deg -- docs/theory/phase6b_closed_loop.md section 5).
        xref[:, 2] = [self._a_ref(th) for th in theta_ref_preview[: self.N + 1]]
        self.p_x0.value = x
        self.p_xref.value = xref
        self.p_ulin.value = u_lin

        n_iter = 0
        while True:
            n_iter += 1
            self._solve(x_lin, u_lin, models, xref)
            if self.mode != "trajectory":
                break
            plan = (self.x.value.copy(), self.u.value.copy())
            new_branches = self._branches_from_plan(x, u_prev, theta_ref_preview, plan, shift=False)
            if new_branches == branches:
                break
            if n_iter >= self.traj_max_iter:
                # cycling: fall back to the measured-state branch for all steps
                self.traj_fallbacks += 1
                br = self.select_branch(x, u_prev, float(theta_ref_preview[0]))
                branches = [br] * self.N
                models = self._models_for(branches, self._eq_model)
                self._solve(x_lin, u_lin, models, xref)
                break
            branches = new_branches
            models = self._models_for(branches, self._eq_model)
        if self.mode == "trajectory":
            self.traj_iters.append(n_iter)
        self._plan = (self.x.value.copy(), self.u.value.copy())
        u0 = float(np.clip(self.u.value[0], self.u_min, self.u_max))
        a_d, b_d, _ = models[0]
        x_pred_next = a_d @ x + b_d[:, 0] * u0 + self.p_d[0].value
        return u0, {"branch": branches[0], "branches": branches,
                    "x_pred_next": x_pred_next, "u_lin": u_lin}


    # ---- hybrid: mode-sequence enumeration ----
    def _violations(self, xp, up, branches):
        """Horizon steps whose planned signs contradict the assumed branch
        (outside the deadbands, where either branch is acceptable)."""
        n = 0
        for j, (act, vel) in enumerate(branches):
            y = up[j] - xp[j, 2]
            if abs(y) > self.act_deadband and (y > 0) != (act == "activating"):
                n += 1
            v = 0.5 * (xp[j, 1] + xp[j + 1, 1])
            if abs(v) > self.vel_deadband and (v > 0) != (vel == "shortening"):
                n += 1
        return n

    def _compute_hybrid(self, x, u_prev, theta_ref_preview):
        """Each switching signal may switch once per horizon, at a step from
        switch_grid (N = never); every candidate mode sequence is solved as
        a QP and the plan with the fewest self-consistency violations, then
        the lowest cost, is applied."""
        theta_ref = float(theta_ref_preview[0])
        a_ref = self._a_ref(theta_ref)
        x_eq = np.array([theta_ref, 0.0, a_ref])
        key = round(theta_ref, 12)
        if getattr(self, "_hyb_key", None) != key:
            f0_eq = joint.state_derivative(x_eq, a_ref)
            self._hyb_models, self._hyb_pf = {}, {}
            for br in linearize.BRANCH_COMBINATIONS:
                a_c, b_c = linearize.linearize(theta_ref, a_ref, *br)
                m = linearize.discretize_affine(a_c, b_c, f0_eq, self.t_s)
                self._hyb_models[br] = m
                pf = solve_discrete_are(m[0], m[1], self.Q, np.array([[self.R]]))
                self._hyb_pf[br] = np.real(sqrtm(pf))
            self._hyb_key = key
        xref = np.zeros((self.N + 1, 3))
        xref[:, 0] = theta_ref_preview[: self.N + 1]
        xref[:, 2] = [self._a_ref(th) for th in theta_ref_preview[: self.N + 1]]
        self.p_x0.value = x
        self.p_xref.value = xref
        self.p_ulin.value = a_ref
        act0, vel0 = self.select_branch(x, u_prev, theta_ref)
        flip = {"activating": "deactivating", "deactivating": "activating",
                "shortening": "lengthening", "lengthening": "shortening"}
        best = None
        for ka in self.switch_grid:
            for kv in self.switch_grid:
                branches = [(act0 if j < ka else flip[act0], vel0 if j < kv else flip[vel0])
                            for j in range(self.N)]
                for j, br in enumerate(branches):
                    a_d, b_d, c_d = self._hyb_models[br]
                    self.p_A[j].value = a_d
                    self.p_B[j].value = b_d[:, 0]
                    self.p_d[j].value = x_eq - a_d @ x_eq - b_d[:, 0] * a_ref + c_d
                pf_sqrt = self._hyb_pf[branches[-1]]
                self.p_Pf_sqrt.value = pf_sqrt
                self.p_Pf_xref.value = pf_sqrt @ xref[self.N]
                try:
                    self.problem.solve(solver=cp.CLARABEL)
                except cp.error.SolverError:
                    self.problem.solve(solver=cp.SCS)
                if self.u.value is None:
                    continue
                cand = (self._violations(self.x.value, self.u.value, branches), self.problem.value)
                if best is None or cand < best[0]:
                    best = (cand, branches, self.x.value.copy(), self.u.value.copy())
        if best is None:
            raise RuntimeError("hybrid MPC: no candidate QP solved")
        (viol, _), branches, xp, up = best
        self.hyb_violations.append(viol)
        self._last_branch = branches[0]
        self._plan = (xp, up)
        u0 = float(np.clip(up[0], self.u_min, self.u_max))
        a_d, b_d, c_d = self._hyb_models[branches[0]]
        x_pred_next = a_d @ x + b_d[:, 0] * u0 + (x_eq - a_d @ x_eq - b_d[:, 0] * a_ref + c_d)
        return u0, {"branch": branches[0], "branches": branches,
                    "x_pred_next": x_pred_next, "u_lin": a_ref}


def simulate_realistic(controller, reference_fn, t_end, x0, noise_deg=0.0, vel_filter_tau=0.02,
                       estimate_activation=False, plant_tau=None, delay_steps=0, seed=0,
                       diverge_deg=30.0, velocity_from_angle=None, compensate_delay=False):
    """Closed loop with the realism factors of docs/theory/phase7_limitations.md
    section 5: Gaussian angle noise (deg), velocity from a first-order
    filtered difference of the noisy angle, activation either measured or
    estimated by integrating the controller's nominal activation ODE from the
    commands it sent, plant activation time constants plant_tau =
    (tau_act, tau_deact) differing from the controller's, and a transport
    delay of delay_steps samples the controller does not model. Velocity is
    measured exactly unless velocity_from_angle (default: whenever there is
    angle noise), in which case it is the filtered difference. Tracking
    metrics use the true state.

    compensate_delay (section 10.1): the controller rolls the measured
    state forward through the inputs already sent but not yet applied,
    with its nominal model, and plans against the reference from the time
    its own input lands (a predictor; exact when the plant is nominal)."""
    if velocity_from_angle is None:
        velocity_from_angle = noise_deg > 0
    rng = np.random.default_rng(seed)
    t_s, N = controller.t_s, controller.N
    n_steps = int(round(t_end / t_s))
    x = np.asarray(x0, dtype=float).copy()
    u_prev = float(x[2])
    queue = [u_prev] * delay_steps
    a_hat, v_hat, th_prev = float(x[2]), float(x[1]), float(x[0])
    alpha = t_s / (vel_filter_tau + t_s)
    t = np.arange(n_steps + 1) * t_s
    log = {"t": t, "x": np.full((n_steps + 1, 3), np.nan), "u": np.full(n_steps, np.nan),
           "theta_ref": reference_fn(t), "residual_theta": np.zeros(n_steps), "branch": [],
           "diverged": False}
    log["x"][0] = x

    def plant_rhs(_, xx, u):
        dx = joint.state_derivative(xx, u)
        if plant_tau is not None:
            dx[2] = muscle.activation_derivative(xx[2], u, tau_act=plant_tau[0], tau_deact=plant_tau[1])
        return dx

    for k in range(n_steps):
        th_m = x[0] + np.deg2rad(noise_deg) * rng.standard_normal()
        if k > 0:
            v_hat = (1 - alpha) * v_hat + alpha * (th_m - th_prev) / t_s
        th_prev = th_m
        x_meas = np.array([th_m, v_hat if velocity_from_angle else x[1],
                           a_hat if estimate_activation else x[2]])
        x_ctrl, t_ctrl = x_meas, k * t_s
        if compensate_delay and delay_steps:
            for u_q in queue:   # applied at steps k .. k+d-1, in order
                x_ctrl = solve_ivp(lambda _, xx, uu: joint.state_derivative(xx, uu), (0, t_s), x_ctrl,
                                   args=(u_q,), rtol=1e-9, atol=1e-11).y[:, -1]
            t_ctrl = (k + delay_steps) * t_s
        preview = reference_fn(t_ctrl + np.arange(N + 1) * t_s)
        u_k, info = controller.compute_control(x_ctrl, u_prev, preview)
        queue.append(u_k)
        u_apply = queue.pop(0)
        x = solve_ivp(plant_rhs, (0, t_s), x, args=(u_apply,), rtol=1e-9, atol=1e-11).y[:, -1]
        if estimate_activation:
            # integrate the input the plant actually received (equal to u_k
            # without delay; every earlier run with estimation had no delay)
            a_hat = float(solve_ivp(lambda _, aa: [muscle.activation_derivative(aa[0], u_apply)],
                                    (0, t_s), [a_hat], rtol=1e-9, atol=1e-11).y[0, -1])
        log["x"][k + 1] = x
        log["u"][k] = u_k
        log["residual_theta"][k] = abs(info["x_pred_next"][0] - x[0])
        log["branch"].append(info["branch"])
        u_prev = u_k
        if (not np.all(np.isfinite(x))) or abs(np.rad2deg(x[0] - log["theta_ref"][k + 1])) > diverge_deg:
            log["diverged"] = True
            break
    return log


def simulate(controller, reference_fn, t_end, x0, u0=None):
    """Closed loop on the true nonlinear plant. reference_fn(t) -> theta_ref
    (rad), vectorized over t. Returns a dict of arrays."""
    t_s, N = controller.t_s, controller.N
    n_steps = int(round(t_end / t_s))
    x = np.asarray(x0, dtype=float).copy()
    u_prev = float(x[2]) if u0 is None else float(u0)
    log = {"t": np.arange(n_steps + 1) * t_s, "x": np.zeros((n_steps + 1, 3)),
           "u": np.zeros(n_steps), "theta_ref": np.zeros(n_steps + 1),
           "residual_theta": np.zeros(n_steps), "residual_x": np.zeros((n_steps, 3)),
           "branch": []}
    log["x"][0] = x
    log["theta_ref"] = reference_fn(log["t"])
    for k in range(n_steps):
        t_k = k * t_s
        preview = reference_fn(t_k + np.arange(N + 1) * t_s)
        u_k, info = controller.compute_control(x, u_prev, preview)
        sol = solve_ivp(lambda t, xx: joint.state_derivative(xx, u_k), (0, t_s), x,
                        rtol=1e-9, atol=1e-11)
        x = sol.y[:, -1]
        log["x"][k + 1] = x
        log["u"][k] = u_k
        log["residual_x"][k] = info["x_pred_next"] - x
        log["residual_theta"][k] = abs(info["x_pred_next"][0] - x[0])
        log["branch"].append(info["branch"])
        u_prev = u_k
    return log


def metrics(log, settle_tol_deg=0.2):
    """Tracking / model-error / effort summary in degrees."""
    err = np.rad2deg(log["x"][:, 0] - log["theta_ref"])
    res = np.rad2deg(log["residual_theta"])
    return {
        "rms_track_deg": float(np.sqrt(np.mean(err ** 2))),
        "peak_track_deg": float(np.max(np.abs(err))),
        "rms_residual_deg": float(np.sqrt(np.mean(res ** 2))),
        "peak_residual_deg": float(np.max(res)),
        "rms_u": float(np.sqrt(np.mean(log["u"] ** 2))),
        "saturated_steps": int(np.sum((log["u"] <= 1e-9) | (log["u"] >= 1 - 1e-9))),
    }
