"""Shared symbolic system definitions + closed-loop simulation for the examples.

Used by the experiment notebooks in ../scripts/ (e.g. simulate_compare.ipynb) so the
dynamics/sets are defined once. Each `*_system()` returns a dict with the symbolic
dynamics f(x), g(x), output h(x), safe/target sets, the FL sampling box, control-channel
labels, the reference input bound, and the integration settings.

Controllers (u_opt) are loaded separately from the MATLAB exports in ../controllers/ and
combined here into the closed loop  xdot = f(x) + g(x) u(x).
"""

from __future__ import annotations

import numpy as np
import sympy as sp


def dubins_system():
    x1, x2, th, v = sp.symbols("x1 x2 th v")
    y1, y2 = sp.symbols("y1 y2")
    state = [x1, x2, th, v]
    f = sp.Matrix([v * sp.cos(th), v * sp.sin(th), 0, 0])
    g = sp.Matrix([[0, 0], [0, 0], [1, 0], [0, 1]])
    h = sp.Matrix([x1, x2])
    h_raw = -(y1**4 + y2**4 - 16) * (y1**4 + y2**4 - 4)
    target_y = (y2) ** 2 + ((y1 + 1.7) / 0.5) ** 2 - 0.4
    safe_y = 1e-3 * (-target_y + 300) * h_raw
    return dict(
        name="dubins_car", state=state, y=[y1, y2], f=f, g=g, h=h,
        safe_y=safe_y, target_y=target_y,
        lo=np.array([-2.0, -2.0, 2 * np.pi / 3, 0.1]),
        hi=np.array([2.0, 2.0, 4 * np.pi / 3, 1.0]),
        channels=["u1 (omega)", "u2 (a)"], u_bound=[5.0, 5.0],
        integ="ivp", t_max=10.0, dt=0.02, sing_tol=None,
    )


def manipulator_system():
    x1, x2, x3, x4 = sp.symbols("x1 x2 x3 x4")
    y1, y2 = sp.symbols("y1 y2")
    state = [x1, x2, x3, x4]
    m1, m2 = 1.0, 1.0
    l1, l2 = 4.0, 4.0
    lc1, lc2 = 2.0, 2.0
    I1, I2 = 0.02, 0.02
    g_grav = 9.81
    M11 = I1 + I2 + m1 * lc1**2 + m2 * (l1**2 + lc2**2 + 2 * l1 * lc2 * sp.cos(x2))
    M12 = m2 * (lc2**2 + l1 * lc2 * sp.cos(x2)) + I2
    M22 = m2 * lc2**2 + I2
    Minv = sp.Matrix([[M11, M12], [M12, M22]]).inv()
    C = sp.Matrix([[-m2 * l1 * lc2 * sp.sin(x2) * x4, -m2 * l1 * lc2 * sp.sin(x2) * (x3 + x4)],
                   [m2 * l1 * lc2 * sp.sin(x2) * x3, sp.Integer(0)]])
    G = sp.Matrix([(m1 * g_grav * lc1 + m2 * g_grav * l1) * sp.cos(x1) + m2 * g_grav * lc2 * sp.cos(x1 + x2),
                   m2 * g_grav * lc2 * sp.cos(x1 + x2)])
    f = sp.Matrix([x3, x4]).col_join(Minv @ (-C @ sp.Matrix([x3, x4]) - G))
    g = sp.Matrix([[0, 0], [0, 0]]).col_join(Minv)
    h = sp.Matrix([l1 * sp.cos(x1) + l2 * sp.cos(x1 + x2),
                   l1 * sp.sin(x1) + l2 * sp.sin(x1 + x2)])
    safe_y = -((4 * (y1 - 2) - 2 * y2**3) ** 2) + 0.8 * y2**3 + 10
    target_y = ((y1 - 2 - 3.5) ** 2 / 1.2**2) + ((y2 - 1.8) ** 2 / 0.4**2) - 2
    return dict(
        name="manipulator", state=state, y=[y1, y2], f=f, g=g, h=h,
        safe_y=safe_y, target_y=target_y,
        lo=np.array([-2.0, 0.8, -0.5, -0.5]),
        hi=np.array([0.5, np.pi - 0.8, 0.5, 0.5]),
        channels=["tau1", "tau2"], u_bound=[500.0, 500.0],
        integ="rk4", t_max=20.0, dt=0.01, sing_tol=0.05,  # stop if |sin(x2)| < sing_tol
    )


SYSTEMS = {"dubins_car": dubins_system, "manipulator": manipulator_system}


def make_callables(sys, u_opt):
    """Lambdify closed-loop dynamics + sets for a given controller u_opt (sympy Matrix)."""
    st = sys["state"]; y = sys["y"]
    sub_y2x = {y[0]: sys["h"][0], y[1]: sys["h"][1]}
    safe_x = sp.sympify(sys["safe_y"]).subs(sub_y2x)
    target_x = sp.sympify(sys["target_y"]).subs(sub_y2x)
    cl = sys["f"] + sys["g"] @ sp.Matrix(u_opt)
    return dict(
        cl=sp.lambdify(st, cl, "numpy"),
        u=[sp.lambdify(st, ui, "numpy") for ui in sp.Matrix(u_opt)],
        safe=sp.lambdify(st, safe_x, "numpy"),
        target=sp.lambdify(st, target_x, "numpy"),
        h=sp.lambdify(st, sys["h"], "numpy"),
    )


def simulate(sys, cb, x0, max_step=None, rtol=1e-8, atol=1e-10):
    """Closed-loop simulate from x0. Returns dict: t, X (n_steps x dim), reached, stayed_safe, max_abs_u.

    Parameters
    ----------
    sys : dict
        System description (from `dubins_system()` / `manipulator_system()`).
    cb : dict
        Lambdified closed-loop callables (from `make_callables`).
    x0 : array-like
        Initial state.
    max_step : float, optional
        Integrator max step size. If None, uses `sys["dt"]` (the default 0.02s
        for Dubins, 0.0025s for the manipulator). Passing a smaller value
        forces finer time resolution — useful when verifying that observed
        behaviour (e.g. transient safe-set violations) isn't a step-size
        artifact.
    rtol, atol : float
        Relative / absolute tolerance for the adaptive integrator (DOP853 for
        Dubins, otherwise fixed-step RK4).
    """
    n = len(sys["state"]); x0 = np.asarray(x0, float)
    dt = sys["dt"] if max_step is None else float(max_step)
    if sys["integ"] == "ivp":
        from scipy.integrate import solve_ivp

        def rhs(t, x):
            d = np.array(cb["cl"](*x)).flatten()
            return d if np.all(np.isfinite(d)) else np.zeros(n)

        def ev(t, x):
            return float(np.squeeze(cb["target"](*x)))
        ev.terminal = True; ev.direction = -1
        sol = solve_ivp(rhs, (0.0, sys["t_max"]), x0, method="DOP853", events=ev,
                        rtol=rtol, atol=atol, max_step=dt, dense_output=False)
        X = sol.y.T; t = sol.t
        reached = sol.status == 1  # terminal event fired
    else:  # rk4 with kinematic-singularity stop
        from functional import simulate as rk4_step
        max_steps = int(sys["t_max"] / dt)
        xs = [x0.copy()]; t = [0.0]; x = x0.reshape(1, -1)
        reached = False
        for k in range(max_steps):
            if sys["sing_tol"] is not None and abs(np.sin(x[0, 1])) < sys["sing_tol"]:
                break
            if float(np.squeeze(cb["target"](*x.T))) <= 0:
                reached = True; break
            x = rk4_step(cb["cl"], x, dt, 3)
            if not np.all(np.isfinite(x)):
                break
            xs.append(x.flatten()); t.append((k + 1) * dt)
        X = np.array(xs); t = np.array(t)
    if len(X) == 0:
        return dict(t=t, X=X, reached=False, stayed_safe=False, max_abs_u=np.nan)
    safe_vals = np.atleast_1d(np.squeeze(cb["safe"](*X.T)))
    stayed_safe = bool(np.all(safe_vals >= 0))
    U = np.array([np.atleast_1d(uf(*X.T)).flatten() for uf in cb["u"]])  # (m, n_steps)
    max_abs_u = np.nanmax(np.abs(U), axis=1) if U.size else np.array([np.nan] * len(cb["u"]))
    return dict(t=t, X=X, reached=bool(reached), stayed_safe=stayed_safe,
                max_abs_u=max_abs_u, U=U)
