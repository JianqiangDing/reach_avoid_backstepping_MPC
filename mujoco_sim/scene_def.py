"""Planar reach-avoid scene definition for the Franka task-space testbed.

Defines (parametrically) the horizontal task plane at z0 and the reach-avoid
geometry -- elliptical workspace bound, target ellipse, circular obstacles --
placed inside the Franka's DEXTEROUS region (away from the +x arm-extension
boundary where the double-integrator abstraction saturates; see the P0.5 probe).

Both sets are represented as a SINGLE polynomial in the output y=(y1,y2):
  - target  phi(y) <= 0 : the goal ellipse.
  - safe    psi(y) >= 0 : w(y) * prod_i o_i(y), where w is one elliptical
    workspace factor and o_i = (y-c_i)^T(y-c_i) - r_i^2 is "outside obstacle i".
The product is exact (no spurious lobes) iff obstacles are pairwise disjoint and
strictly interior to the workspace ellipse -- enforced by validate_scene().

Coordinates are world (x, y) in the plane at height z0.
"""

from __future__ import annotations

import os
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))

# output-space symbols (shared by the composed polynomials)
Y1, Y2 = sp.symbols("y1 y2")


# default scene: elliptical workspace in the dexterous region (x<=0.62 to avoid
# the +x extension boundary; home EE is at (0.554, 0, 0.625)). Target on the left,
# two obstacles straddling a central corridor between start (right) and target.
SCENE = dict(
    z0=0.625,
    workspace=dict(center=(0.46, 0.0), semi=(0.16, 0.26)),   # ellipse: x in [0.30,0.62], y in [-0.26,0.26]
    start=(0.56, 0.0),
    target=dict(center=(0.36, 0.0), radii=(0.045, 0.045)),
    obstacles=[dict(center=(0.46, 0.09), radius=0.05),
               dict(center=(0.46, -0.09), radius=0.05)],
    a_max=1.0,
)


# ---- single-polynomial set composition --------------------------------------
def _ws_factor(scene, y1, y2):
    cx, cy = scene["workspace"]["center"]; ax, ay = scene["workspace"]["semi"]
    return 1 - ((y1 - cx) / ax) ** 2 - ((y2 - cy) / ay) ** 2


def compose_target_poly(scene=SCENE):
    """phi(y) <= 0 inside the target ellipse (single quadratic)."""
    cx, cy = scene["target"]["center"]; rx, ry = scene["target"]["radii"]
    return sp.expand(((Y1 - cx) / rx) ** 2 + ((Y2 - cy) / ry) ** 2 - 1)


def compose_safe_poly(scene=SCENE):
    """psi(y) >= 0 inside the safe set = workspace ellipse minus the obstacles."""
    psi = _ws_factor(scene, Y1, Y2)
    for ob in scene["obstacles"]:
        ox, oy = ob["center"]; r = ob["radius"]
        psi = psi * ((Y1 - ox) ** 2 + (Y2 - oy) ** 2 - r ** 2)
    return sp.expand(psi)


def target_func(scene=SCENE):
    """Vectorized callable f(y1, y2) -> phi value(s)."""
    return sp.lambdify((Y1, Y2), compose_target_poly(scene), "numpy")


def safe_func(scene=SCENE):
    """Vectorized callable f(y1, y2) -> psi value(s)."""
    return sp.lambdify((Y1, Y2), compose_safe_poly(scene), "numpy")


def target_phi(y, scene=SCENE):
    return float(target_func(scene)(y[0], y[1]))


def safe_psi(y, scene=SCENE):
    return float(safe_func(scene)(y[0], y[1]))


def validate_scene(scene=SCENE, n_boundary=72):
    """Ensure obstacles are pairwise disjoint and obstacles+target are strictly
    interior to the workspace ellipse (so the product safe set is exact)."""
    cx, cy = scene["workspace"]["center"]; ax, ay = scene["workspace"]["semi"]

    def w(px, py):
        return 1 - ((px - cx) / ax) ** 2 - ((py - cy) / ay) ** 2

    th = np.linspace(0, 2 * np.pi, n_boundary, endpoint=False)
    obs = scene["obstacles"]
    for i, ob in enumerate(obs):
        ox, oy = ob["center"]; r = ob["radius"]
        if not np.all(w(ox + r * np.cos(th), oy + r * np.sin(th)) > 0):
            raise ValueError(f"obstacle {i} not strictly inside workspace ellipse")
    for i in range(len(obs)):
        for j in range(i + 1, len(obs)):
            ci = np.array(obs[i]["center"]); cj = np.array(obs[j]["center"])
            if np.linalg.norm(ci - cj) <= obs[i]["radius"] + obs[j]["radius"]:
                raise ValueError(f"obstacles {i} and {j} overlap")
    tcx, tcy = scene["target"]["center"]; trx, try_ = scene["target"]["radii"]
    if not np.all(w(tcx + trx * np.cos(th), tcy + try_ * np.sin(th)) > 0):
        raise ValueError("target not strictly inside workspace ellipse")
    return True


def _geom(s):
    return "    " + s + "\n"


def generate_mjcf(scene=SCENE, path=None):
    z0 = scene["z0"]
    ws = scene["workspace"]; cx, cy = ws["center"]; hx, hy = ws["semi"]
    tg = scene["target"]; tcx, tcy = tg["center"]; trx, try_ = tg["radii"]
    sx, sy = scene["start"]

    body = ""
    # translucent task plane (workspace region)
    body += _geom(f'<geom name="taskplane" type="box" pos="{cx} {cy} {z0}" '
                  f'size="{hx} {hy} 0.002" contype="0" conaffinity="0" '
                  f'rgba="0.5 0.55 0.6 0.12"/>')
    # workspace outline (4 thin bars)
    for nm, px, py, sxh, syh in [
        ("ws_n", cx, cy + hy, hx, 0.004), ("ws_s", cx, cy - hy, hx, 0.004),
        ("ws_e", cx + hx, cy, 0.004, hy), ("ws_w", cx - hx, cy, 0.004, hy)]:
        body += _geom(f'<geom name="{nm}" type="box" pos="{px} {py} {z0}" '
                      f'size="{sxh} {syh} 0.004" contype="0" conaffinity="0" '
                      f'rgba="0.35 0.4 0.5 0.8"/>')
    # target ellipse (flat ellipsoid)
    body += _geom(f'<geom name="target" type="ellipsoid" pos="{tcx} {tcy} {z0}" '
                  f'size="{trx} {try_} 0.003" contype="0" conaffinity="0" '
                  f'rgba="0.2 0.8 0.3 0.55"/>')
    # obstacles (red posts spanning z0 +/- 0.1)
    for i, ob in enumerate(scene["obstacles"]):
        ox, oy = ob["center"]
        body += _geom(f'<geom name="obstacle{i}" type="cylinder" pos="{ox} {oy} {z0}" '
                      f'size="{ob["radius"]} 0.10" contype="0" conaffinity="0" '
                      f'rgba="0.85 0.2 0.2 0.7"/>')
    # start marker (blue sphere)
    body += _geom(f'<geom name="start" type="sphere" pos="{sx} {sy} {z0}" '
                  f'size="0.02" contype="0" conaffinity="0" rgba="0.2 0.4 0.9 0.9"/>')

    xml = f"""<mujoco model="franka reach-avoid scene">
  <include file="panda.xml"/>

  <statistic center="0.3 0 0.5" extent="1.1"/>

  <visual>
    <headlight diffuse="0.6 0.6 0.6" ambient="0.3 0.3 0.3" specular="0 0 0"/>
    <rgba haze="0.15 0.25 0.35 1"/>
    <global azimuth="120" elevation="-20" offwidth="1280" offheight="960"/>
  </visual>

  <asset>
    <texture type="skybox" builtin="gradient" rgb1="0.3 0.5 0.7" rgb2="0 0 0" width="512" height="3072"/>
    <texture type="2d" name="groundplane" builtin="checker" mark="edge" rgb1="0.2 0.3 0.4" rgb2="0.1 0.2 0.3"
      markrgb="0.8 0.8 0.8" width="300" height="300"/>
    <material name="groundplane" texture="groundplane" texuniform="true" texrepeat="5 5" reflectance="0.2"/>
  </asset>

  <worldbody>
    <light pos="0 0 1.5" dir="0 0 -1" directional="true"/>
    <geom name="floor" size="0 0 0.05" type="plane" material="groundplane"/>
{body}  </worldbody>
</mujoco>
"""
    path = path or os.path.join(HERE, "models", "franka_emika_panda", "scene_reachavoid.xml")
    with open(path, "w") as f:
        f.write(xml)
    return path


if __name__ == "__main__":
    validate_scene(SCENE)
    p = generate_mjcf()
    print(f"validated + wrote {p}")
    s = np.array(SCENE["start"])
    print(f"start safe_psi={safe_psi(s):.4f} (>0 ok)  target_phi={target_phi(s):.4f} (>0 = not reached)")
    tc = np.array(SCENE["target"]["center"])
    print(f"target-center safe_psi={safe_psi(tc):.4f}  target_phi={target_phi(tc):.4f} (<0 = inside target)")
