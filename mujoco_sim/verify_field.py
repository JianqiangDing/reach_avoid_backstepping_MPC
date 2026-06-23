"""Milestone 0, part 2: confirm the single-integrator k1 field is a reach-avoid flow.

Integrates the output-space single-integrator dynamics y_dot = k1(y) from safe
starts and checks that streamlines flow into the target without leaving the safe
set. For the holonomic double integrator the k1 field should itself be reach-avoid.
Run: python verify_field.py  (triggers the MATLAB synthesis via synth.synthesize).
"""

import os
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scene_def as S

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures", "milestone0")


def k1_callable(bundle):
    """Lambdify k1_opt (list[2] in y1,y2) into a callable y -> ndarray(2)."""
    f0 = sp.lambdify((S.Y1, S.Y2), bundle["k1_opt"][0], "numpy")
    f1 = sp.lambdify((S.Y1, S.Y2), bundle["k1_opt"][1], "numpy")
    return lambda y: np.array([float(f0(y[0], y[1])), float(f1(y[0], y[1]))])


def _safe_starts(scene, n, max_tries=100000):
    # psi is a product polynomial with arbitrary scale, so test its SIGN (> 0),
    # not a magnitude threshold.
    safe = S.safe_func(scene); phi = S.target_func(scene)
    cx, cy = scene["workspace"]["center"]; ax, ay = scene["workspace"]["semi"]
    rng = np.random.default_rng(0); pts = []
    for _ in range(max_tries):
        if len(pts) >= n:
            break
        p = np.array([rng.uniform(cx - ax, cx + ax), rng.uniform(cy - ay, cy + ay)])
        if safe(p[0], p[1]) > 0 and phi(p[0], p[1]) > 0:   # strictly safe and not already in target
            pts.append(p)
    if len(pts) < n:
        raise RuntimeError(f"only sampled {len(pts)}/{n} safe starts")
    return np.array(pts)


def field_reach_avoid(scene, k1, n_starts=40, t_max=20.0, dt=0.02):
    """Integrate y_dot=k1(y) from safe starts; report reach / leave-safe fractions."""
    safe = S.safe_func(scene); phi = S.target_func(scene)
    reached = left = 0
    for p0 in _safe_starts(scene, n_starts):
        p = p0.copy()
        for _ in range(int(t_max / dt)):
            if safe(p[0], p[1]) < 0:
                left += 1; break
            if phi(p[0], p[1]) <= 0:
                reached += 1; break
            p = p + dt * k1(p)
        # ran out of time without reaching -> neither reached nor left
    return dict(success_frac=reached / n_starts, left_safe_frac=left / n_starts)


def plot_field(scene, k1, out=None):
    os.makedirs(FIGDIR, exist_ok=True)
    out = out or os.path.join(FIGDIR, "k1_field.png")
    cx, cy = scene["workspace"]["center"]; ax, ay = scene["workspace"]["semi"]
    gx, gy = np.meshgrid(np.linspace(cx - ax, cx + ax, 26), np.linspace(cy - ay, cy + ay, 26))
    U = np.zeros_like(gx); Vv = np.zeros_like(gy)
    for i in range(gx.shape[0]):
        for j in range(gx.shape[1]):
            d = k1(np.array([gx[i, j], gy[i, j]])); U[i, j], Vv[i, j] = d
    psi = S.safe_func(scene)(gx, gy); phi = S.target_func(scene)(gx, gy)
    fig, axx = plt.subplots(figsize=(6, 6))
    axx.streamplot(gx, gy, U, Vv, density=1.2, color="0.4")
    axx.contour(gx, gy, psi, levels=[0], colors="green")
    axx.contour(gx, gy, phi, levels=[0], colors="red")
    axx.set_aspect("equal"); axx.set_title("single-integrator k1 field + sets")
    axx.set_xlabel("y1"); axx.set_ylabel("y2")
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
    return out


def main():
    import synth
    sc = S.SCENE
    bundle = synth.synthesize(sc)
    k1 = k1_callable(bundle)
    r = field_reach_avoid(sc, k1)
    out = plot_field(sc, k1)
    print(f"k1 field: success_frac={r['success_frac']:.2f} left_safe_frac={r['left_safe_frac']:.2f}")
    print(f"figure -> {out}")
    print("PASS" if (r["success_frac"] > 0.8 and r["left_safe_frac"] < 0.05) else "FAIL")


if __name__ == "__main__":
    main()
