"""Milestone 0, part 1: confirm the composed single-polynomial sets represent the scene.

Plots the zero-level sets of the composed psi (safe) and phi (target), and checks
on a grid that sign(psi) matches the boolean "inside workspace AND outside every
obstacle". Run: python verify_sets.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scene_def as S

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures", "milestone0")


def _grid(scene, n):
    cx, cy = scene["workspace"]["center"]; ax, ay = scene["workspace"]["semi"]
    gx, gy = np.meshgrid(np.linspace(cx - ax - 0.05, cx + ax + 0.05, n),
                         np.linspace(cy - ay - 0.05, cy + ay + 0.05, n))
    return gx, gy


def set_consistency(scene=S.SCENE, n=120):
    """Fraction of grid points where sign(psi) disagrees with the boolean safe set."""
    gx, gy = _grid(scene, n)
    poly_safe = S.safe_func(scene)(gx, gy) >= 0
    cx, cy = scene["workspace"]["center"]; ax, ay = scene["workspace"]["semi"]
    truth = ((gx - cx) / ax) ** 2 + ((gy - cy) / ay) ** 2 <= 1
    for ob in scene["obstacles"]:
        ox, oy = ob["center"]
        truth &= (gx - ox) ** 2 + (gy - oy) ** 2 >= ob["radius"] ** 2
    mism = float(np.mean(poly_safe != truth))
    return dict(agree=(mism == 0.0), mismatch_frac=mism)


def plot_sets(scene=S.SCENE, out=None):
    os.makedirs(FIGDIR, exist_ok=True)
    out = out or os.path.join(FIGDIR, "sets.png")
    gx, gy = _grid(scene, 300)
    psi = S.safe_func(scene)(gx, gy)
    phi = S.target_func(scene)(gx, gy)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.contourf(gx, gy, (psi >= 0).astype(float), levels=[0.5, 1.5], colors=["#bfe3c0"])
    ax.contour(gx, gy, psi, levels=[0], colors="green")
    ax.contour(gx, gy, phi, levels=[0], colors="red")
    sx, sy = scene["start"]
    ax.plot(sx, sy, "b*", ms=12, label="start")
    ax.set_aspect("equal"); ax.legend()
    ax.set_title("safe {psi>=0} (green) + target {phi<=0} (red)")
    ax.set_xlabel("y1"); ax.set_ylabel("y2")
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
    return out


def main():
    sc = S.SCENE
    S.validate_scene(sc)
    r = set_consistency(sc)
    out = plot_sets(sc)
    print(f"validate: OK   consistency: agree={r['agree']} mismatch_frac={r['mismatch_frac']}")
    print(f"figure -> {out}")
    print("PASS" if r["agree"] else "FAIL")


if __name__ == "__main__":
    main()
