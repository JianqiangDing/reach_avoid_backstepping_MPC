"""Phase 1.3 shared logic (self-contained; imports nothing from the main repo).

Pure functions used by dubins_synth.ipynb / manipulator_synth.ipynb to:
  * define the frozen example dynamics / sets symbolically (mirrors the .m
    example definition and python/systems.py, but lives inside revision/);
  * lambdify the closed loop  xdot = f(x) + g(x) u_opt(x), the reach-avoid
    certificate V(x), and the safe/target sets from an imported controller;
  * sample the synthesis box and simulate a closed-form rollout.

The reach-avoid set used everywhere is the certificate zero-superlevel set
{ x : V(x) >= 0 }; the controller and certificate come from the MATLAB export
(revision/controllers/<example>_clean.py).
"""

from __future__ import annotations

import numpy as np
import sympy as sp


# ----------------------------------------------------------------------------
# Frozen example definitions (state symbols, dynamics, output map, sets).
# These mirror solve_*_clean.m / python systems exactly; the controller and the
# synthesis box (X_S_eff) are supplied separately by the notebook.
# ----------------------------------------------------------------------------
def dubins_system():
    x1, x2, th, v = sp.symbols("x1 x2 th v")
    y1, y2 = sp.symbols("y1 y2")
    f = sp.Matrix([v * sp.cos(th), v * sp.sin(th), 0, 0])
    g = sp.Matrix([[0, 0], [0, 0], [1, 0], [0, 1]])
    h = sp.Matrix([x1, x2])
    h_raw = -(y1**4 + y2**4 - 16) * (y1**4 + y2**4 - 4)
    target_y = (y2) ** 2 + ((y1 + 1.7) / 0.5) ** 2 - 0.4
    safe_y = 1e-3 * (-target_y + 300) * h_raw
    return dict(
        name="dubins", state=[x1, x2, th, v], y=[y1, y2], f=f, g=g, h=h,
        safe_y=safe_y, target_y=target_y,
        channels=["omega", "a"], proj=[(0, 1, "x1", "x2"), (2, 3, "theta", "v")],
        integ="ivp", t_max=10.0, dt=0.02, sing_tol=None,
    )


def manipulator_system(consts=None):
    c = dict(m1=1.0, m2=1.0, l1=4.0, l2=4.0, lc1=2.0, lc2=2.0, I1=0.02, I2=0.02, g=9.81)
    if consts:
        c.update({k: float(v) for k, v in consts.items() if k in c})
    x1, x2, x3, x4 = sp.symbols("x1 x2 x3 x4")
    y1, y2 = sp.symbols("y1 y2")
    m1, m2, l1, l2, lc1, lc2, I1, I2, gg = (
        c["m1"], c["m2"], c["l1"], c["l2"], c["lc1"], c["lc2"], c["I1"], c["I2"], c["g"]
    )
    M11 = I1 + I2 + m1 * lc1**2 + m2 * (l1**2 + lc2**2 + 2 * l1 * lc2 * sp.cos(x2))
    M12 = m2 * (lc2**2 + l1 * lc2 * sp.cos(x2)) + I2
    M22 = m2 * lc2**2 + I2
    Minv = sp.Matrix([[M11, M12], [M12, M22]]).inv()
    C = sp.Matrix([[-m2 * l1 * lc2 * sp.sin(x2) * x4, -m2 * l1 * lc2 * sp.sin(x2) * (x3 + x4)],
                   [m2 * l1 * lc2 * sp.sin(x2) * x3, sp.Integer(0)]])
    G = sp.Matrix([(m1 * gg * lc1 + m2 * gg * l1) * sp.cos(x1) + m2 * gg * lc2 * sp.cos(x1 + x2),
                   m2 * gg * lc2 * sp.cos(x1 + x2)])
    f = sp.Matrix([x3, x4]).col_join(Minv @ (-C @ sp.Matrix([x3, x4]) - G))
    g = sp.Matrix([[0, 0], [0, 0]]).col_join(Minv)
    h = sp.Matrix([l1 * sp.cos(x1) + l2 * sp.cos(x1 + x2),
                   l1 * sp.sin(x1) + l2 * sp.sin(x1 + x2)])
    safe_y = -((4 * (y1 - 2) - 2 * y2**3) ** 2) + 0.8 * y2**3 + 10
    target_y = ((y1 - 2 - 3.5) ** 2 / 1.2**2) + ((y2 - 1.8) ** 2 / 0.4**2) - 2
    return dict(
        name="manipulator", state=[x1, x2, x3, x4], y=[y1, y2], f=f, g=g, h=h,
        safe_y=safe_y, target_y=target_y,
        channels=["tau1", "tau2"], proj=[(0, 1, "q1", "q2"), (2, 3, "dq1", "dq2")],
        integ="rk4", t_max=20.0, dt=0.01, sing_tol=0.05,  # stop if |sin(x2)| < sing_tol
    )


# ----------------------------------------------------------------------------
# Lambdify closed loop + sets + certificate for an imported controller.
# ----------------------------------------------------------------------------
def make_callables(sys, u_opt, certificate_opt):
    """u_opt: list/Matrix of sympy exprs (state coords); certificate_opt: sympy expr.

    Returns numpy callables f(*state) for the closed loop, each input channel,
    the safe/target sets, the output map h, and the reach-avoid certificate V.
    """
    st = sys["state"]
    y = sys["y"]
    sub_y2x = {y[0]: sys["h"][0], y[1]: sys["h"][1]}
    safe_x = sp.sympify(sys["safe_y"]).subs(sub_y2x)
    target_x = sp.sympify(sys["target_y"]).subs(sub_y2x)
    u_mat = sp.Matrix(list(u_opt))
    cl = sys["f"] + sys["g"] @ u_mat
    return dict(
        cl=sp.lambdify(st, cl, "numpy"),
        u=[sp.lambdify(st, ui, "numpy") for ui in u_mat],
        safe=sp.lambdify(st, safe_x, "numpy"),
        target=sp.lambdify(st, target_x, "numpy"),
        h=sp.lambdify(st, sys["h"], "numpy"),
        V=sp.lambdify(st, sp.sympify(certificate_opt), "numpy"),
    )


# ----------------------------------------------------------------------------
# Sampling + evaluation helpers.
# ----------------------------------------------------------------------------
def sample_box(lo, hi, n, seed=0):
    """Uniform samples in the box [lo, hi]; returns (n, dim)."""
    rng = np.random.default_rng(seed)
    lo = np.asarray(lo, float)
    hi = np.asarray(hi, float)
    return lo + (hi - lo) * rng.random((n, len(lo)))


def evaluate(cb, X):
    """Vectorized eval over X (n, dim). Returns V, safe, target, |u| (m, n)."""
    cols = [X[:, k] for k in range(X.shape[1])]
    nX = X.shape[0]
    V = np.broadcast_to(np.atleast_1d(np.squeeze(cb["V"](*cols))).astype(float), (nX,)).copy()
    safe = np.broadcast_to(np.atleast_1d(np.squeeze(cb["safe"](*cols))).astype(float), (nX,)).copy()
    target = np.broadcast_to(np.atleast_1d(np.squeeze(cb["target"](*cols))).astype(float), (nX,)).copy()
    U = np.array([np.broadcast_to(np.atleast_1d(uf(*cols)).flatten().astype(float), (nX,))
                  for uf in cb["u"]])  # (m, n)
    return V, safe, target, U


def _rk4_step(cl, x, dt):
    def d(xx):
        return np.array(cl(*xx)).flatten()
    k1 = d(x)
    k2 = d(x + 0.5 * dt * k1)
    k3 = d(x + 0.5 * dt * k2)
    k4 = d(x + dt * k3)
    return x + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def simulate(sys, cb, x0, t_max=None, dt=None, rtol=1e-8, atol=1e-10):
    """Closed-loop rollout from x0. Returns dict with t, X (n,dim), reached,
    stayed_safe, max_abs_u (m,), U (m,n), stop_reason."""
    n = len(sys["state"])
    x0 = np.asarray(x0, float)
    t_max = sys["t_max"] if t_max is None else float(t_max)
    dt = sys["dt"] if dt is None else float(dt)
    stop = "t_max"
    if sys["integ"] == "ivp":
        from scipy.integrate import solve_ivp

        def rhs(t, x):
            dd = np.array(cb["cl"](*x)).flatten()
            return dd if np.all(np.isfinite(dd)) else np.zeros(n)

        def ev(t, x):
            return float(np.squeeze(cb["target"](*x)))
        ev.terminal = True
        ev.direction = -1
        sol = solve_ivp(rhs, (0.0, t_max), x0, method="DOP853", events=ev,
                        rtol=rtol, atol=atol, max_step=dt)
        X = sol.y.T
        t = sol.t
        reached = sol.status == 1
        if reached:
            stop = "reached_target"
    else:  # rk4 with kinematic-singularity stop (manipulator)
        steps = int(t_max / dt)
        xs = [x0.copy()]
        t = [0.0]
        x = x0.copy()
        reached = False
        for k in range(steps):
            if sys["sing_tol"] is not None and abs(np.sin(x[1])) < sys["sing_tol"]:
                stop = "near_singularity"
                break
            if float(np.squeeze(cb["target"](*x))) <= 0:
                reached = True
                stop = "reached_target"
                break
            x = _rk4_step(cb["cl"], x, dt)
            if not np.all(np.isfinite(x)):
                stop = "diverged"
                break
            xs.append(x.copy())
            t.append((k + 1) * dt)
        X = np.array(xs)
        t = np.array(t)
    if len(X) == 0:
        return dict(t=t, X=X, reached=False, stayed_safe=False,
                    max_abs_u=np.array([np.nan] * len(cb["u"])), U=np.empty((len(cb["u"]), 0)),
                    stop_reason="empty")
    nX = X.shape[0]
    safe_vals = np.broadcast_to(np.atleast_1d(np.squeeze(cb["safe"](*X.T))), (nX,))
    stayed_safe = bool(np.all(safe_vals >= 0))
    U = np.array([np.broadcast_to(np.atleast_1d(uf(*X.T)).flatten(), (nX,)) for uf in cb["u"]])  # (m, n)
    max_abs_u = np.nanmax(np.abs(U), axis=1) if U.size else np.array([np.nan] * len(cb["u"]))
    return dict(t=t, X=X, reached=bool(reached), stayed_safe=stayed_safe,
                max_abs_u=max_abs_u, U=U, stop_reason=stop)
