"""Planar reach-avoid scene definition for the Franka task-space testbed.

Defines (parametrically) the horizontal task plane at z0 and the reach-avoid
geometry -- workspace bound, target ellipse, circular obstacles -- placed inside
the Franka's DEXTEROUS region (away from the +x arm-extension boundary where the
double-integrator abstraction saturates; see the P0.5 probe).

Provides:
  - SCENE: the default scene config (numbers chosen in the dexterous box)
  - numpy set evaluators: target_phi(y)<=0 inside target, safe_psi(y)>=0 inside safe
  - generate_mjcf(): writes scene_reachavoid.xml (Franka + floor + decorative
    geoms for the sets). Geoms are contype/conaffinity=0 (visual only): obstacles
    are enforced by the controller's safe-set constraint, not physical contact.

Coordinates are world (x, y) in the plane at height z0.
"""

from __future__ import annotations

import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


# default scene: workspace box centered in the dexterous region, x<=0.62 to avoid
# the +x extension boundary (home EE is at (0.554, 0, 0.625)).
SCENE = dict(
    z0=0.625,
    workspace=dict(center=(0.46, 0.0), half=(0.16, 0.28)),   # x in [0.30,0.62], y in [-0.28,0.28]
    start=(0.554, 0.0),
    target=dict(center=(0.36, 0.19), radii=(0.07, 0.06)),
    obstacles=[dict(center=(0.47, 0.10), radius=0.05),
               dict(center=(0.45, -0.12), radius=0.05)],
    a_max=1.0,
)


def target_phi(y, scene=SCENE):
    """phi(y) <= 0 inside the target ellipse."""
    cx, cy = scene["target"]["center"]; rx, ry = scene["target"]["radii"]
    return ((y[0] - cx) / rx) ** 2 + ((y[1] - cy) / ry) ** 2 - 1.0


def safe_psi(y, scene=SCENE):
    """psi(y) >= 0 inside the safe set (inside workspace AND outside every obstacle)."""
    cx, cy = scene["workspace"]["center"]; hx, hy = scene["workspace"]["half"]
    inside = min(hx ** 2 - (y[0] - cx) ** 2, hy ** 2 - (y[1] - cy) ** 2)
    val = inside
    for ob in scene["obstacles"]:
        ox, oy = ob["center"]
        val = min(val, (y[0] - ox) ** 2 + (y[1] - oy) ** 2 - ob["radius"] ** 2)
    return val


def _geom(s):
    return "    " + s + "\n"


def generate_mjcf(scene=SCENE, path=None):
    z0 = scene["z0"]
    ws = scene["workspace"]; cx, cy = ws["center"]; hx, hy = ws["half"]
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
    p = generate_mjcf()
    print(f"wrote {p}")
    # quick sanity: start should be safe and not yet at target
    s = np.array(SCENE["start"])
    print(f"start safe_psi={safe_psi(s):.4f} (>0 ok)  target_phi={target_phi(s):.4f} (>0 = not reached)")
    tc = np.array(SCENE["target"]["center"])
    print(f"target-center safe_psi={safe_psi(tc):.4f}  target_phi={target_phi(tc):.4f} (<0 = inside target)")
