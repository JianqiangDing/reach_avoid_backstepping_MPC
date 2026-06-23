# Franka Reach-Avoid P1 — Milestone 0 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and verify, in 2D, that a Franka planar reach-avoid scene yields a sound single-polynomial safe/target representation and a single-integrator reach-avoid vector field — the gate that must pass before any double-integrator / MPC / MuJoCo work.

**Architecture:** Extend `mujoco_sim/scene_def.py` to compose single-polynomial safe/target sets (elliptical workspace × obstacle-exterior products) with validation. Add a parameterized MATLAB SOS synthesis (`matlab/example_franka_planar.m`) driven from Python (`mujoco_sim/synth.py`). Two verification scripts (`verify_sets.py`, `verify_field.py`) produce the 2D figures and numeric pass/fail checks.

**Tech Stack:** Python 3.11 (conda env `rab_mpc`: numpy, scipy, sympy, matplotlib), MATLAB + SOSTOOLS + MOSEK (existing repo synthesis toolchain), pytest.

## Global Constraints

- Run Python with the `rab_mpc` env: `/home/jianqiang/miniconda3/envs/rab_mpc/bin/python3` (or `conda activate rab_mpc`). All code/comments in English.
- All new code is self-contained under `mujoco_sim/` and `matlab/`; no runtime dependency on the paper repo.
- Sets are SINGLE polynomials in output `y=(y1,y2)`: target `φ(y)≤0` (ellipse); safe `ψ(y)=w(y)·∏ᵢoᵢ(y)≥0` where `w` is one **elliptical** workspace factor `1−((y1−cx)/ax)²−((y2−cy)/ay)²` and `oᵢ=(y1−oxᵢ)²+(y2−oyᵢ)²−rᵢ²`. Workspace is an ellipse, never a box (box-as-product-of-slabs has spurious corner lobes).
- Scene validity (enforced): obstacles pairwise disjoint AND strictly interior to the workspace ellipse; target strictly interior to the workspace ellipse.
- Fixed constants: `z0=0.625`, `a_max=1.0`, dexterous box `x∈[0.30,0.62], y∈[-0.28,0.28]` (all geometry inside it).
- MATLAB tasks require MATLAB + SOSTOOLS + MOSEK on the path. Model `example_franka_planar.m` on the existing `matlab/example_dubins_car.m` (same call sequence: define `fx_sym, gx_sym, hx_sym, x_vars_sym, y_vars_sym, safe_set_sym, target_set_sym, lb, ub`, degrees, `mu_val`, `samples_num`, then synthesise + export).
- Tests live in `mujoco_sim/tests/`, run with pytest. If pytest is missing: `/home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 -m pip install pytest`.
- **GATE:** do not begin controllers / MPC / MuJoCo closed loop until `verify_sets.py` and `verify_field.py` both pass (Task 8).

---

### Task 1: Single-polynomial set composition + scene validation in `scene_def.py`

**Files:**
- Modify: `mujoco_sim/scene_def.py` (workspace box→ellipse; add `Y1,Y2`, `validate_scene`, `compose_target_poly`, `compose_safe_poly`, `safe_func`, `target_func`; rewrite `safe_psi`/`target_phi` to use them)
- Test: `mujoco_sim/tests/test_scene_def.py`

**Interfaces:**
- Produces:
  - `Y1, Y2` — module-level `sympy.Symbol`s for `y1, y2`.
  - `compose_target_poly(scene=SCENE) -> sympy.Expr` (φ in `Y1,Y2`).
  - `compose_safe_poly(scene=SCENE) -> sympy.Expr` (ψ in `Y1,Y2`).
  - `target_func(scene=SCENE) -> callable(y1,y2)->ndarray`, `safe_func(scene=SCENE) -> callable(y1,y2)->ndarray` (vectorized, via `sympy.lambdify`).
  - `validate_scene(scene=SCENE) -> True` or raises `ValueError`.
  - `safe_psi(y, scene=SCENE) -> float`, `target_phi(y, scene=SCENE) -> float` (scalar wrappers).
  - `SCENE["workspace"]` now uses key `semi=(ax,ay)` (was `half`).

- [ ] **Step 1: Write the failing test**

```python
# mujoco_sim/tests/test_scene_def.py
import numpy as np
import pytest
import scene_def as S


def test_safe_sign_inside_outside():
    sc = S.SCENE
    # a point in the workspace, clear of obstacles -> safe > 0
    assert S.safe_psi((0.54, 0.0), sc) > 0
    # inside obstacle 0 -> unsafe < 0
    ox, oy = sc["obstacles"][0]["center"]
    assert S.safe_psi((ox, oy), sc) < 0
    # outside the workspace ellipse -> unsafe < 0
    assert S.safe_psi((1.2, 0.0), sc) < 0


def test_target_sign():
    sc = S.SCENE
    tc = sc["target"]["center"]
    assert S.target_phi(tc, sc) < 0          # target center is inside target
    assert S.target_phi((0.54, 0.0), sc) > 0  # start is not in target


def test_product_matches_truth_on_grid():
    sc = S.SCENE
    safe = S.safe_func(sc)
    cx, cy = sc["workspace"]["center"]; ax, ay = sc["workspace"]["semi"]
    gx, gy = np.meshgrid(np.linspace(cx - ax, cx + ax, 60),
                         np.linspace(cy - ay, cy + ay, 60))
    poly_safe = safe(gx, gy) >= 0
    in_ws = ((gx - cx) / ax) ** 2 + ((gy - cy) / ay) ** 2 <= 1
    out_obs = np.ones_like(in_ws, bool)
    for ob in sc["obstacles"]:
        ox, oy = ob["center"]
        out_obs &= (gx - ox) ** 2 + (gy - oy) ** 2 >= ob["radius"] ** 2
    truth = in_ws & out_obs
    # product sign agrees with the boolean safe region everywhere on the grid
    assert np.array_equal(poly_safe, truth)


def test_validate_rejects_overlap_and_exterior():
    assert S.validate_scene(S.SCENE) is True
    bad_overlap = {**S.SCENE, "obstacles": [
        {"center": (0.46, 0.0), "radius": 0.06},
        {"center": (0.49, 0.0), "radius": 0.06}]}
    with pytest.raises(ValueError):
        S.validate_scene(bad_overlap)
    bad_exterior = {**S.SCENE, "obstacles": [{"center": (0.46, 0.30), "radius": 0.05}]}
    with pytest.raises(ValueError):
        S.validate_scene(bad_exterior)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd mujoco_sim && /home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 -m pytest tests/test_scene_def.py -v`
Expected: FAIL (`AttributeError`/`KeyError`: `compose_safe_poly`, `safe_func`, `validate_scene`, or `semi` not present).

- [ ] **Step 3: Write minimal implementation**

In `mujoco_sim/scene_def.py`, set the workspace to an ellipse and add the composition/validation API. Replace the `SCENE` workspace entry and the numpy `target_phi`/`safe_psi` with:

```python
import sympy as sp

Y1, Y2 = sp.symbols("y1 y2")

# in SCENE: workspace=dict(center=(0.46, 0.0), semi=(0.16, 0.26)),
# (replace the old `half=(...)` key with `semi=(...)`)

def _ws_factor(scene, y1, y2):
    cx, cy = scene["workspace"]["center"]; ax, ay = scene["workspace"]["semi"]
    return 1 - ((y1 - cx) / ax) ** 2 - ((y2 - cy) / ay) ** 2


def compose_target_poly(scene=SCENE):
    cx, cy = scene["target"]["center"]; rx, ry = scene["target"]["radii"]
    return sp.expand(((Y1 - cx) / rx) ** 2 + ((Y2 - cy) / ry) ** 2 - 1)


def compose_safe_poly(scene=SCENE):
    psi = _ws_factor(scene, Y1, Y2)
    for ob in scene["obstacles"]:
        ox, oy = ob["center"]; r = ob["radius"]
        psi = psi * ((Y1 - ox) ** 2 + (Y2 - oy) ** 2 - r ** 2)
    return sp.expand(psi)


def target_func(scene=SCENE):
    return sp.lambdify((Y1, Y2), compose_target_poly(scene), "numpy")


def safe_func(scene=SCENE):
    return sp.lambdify((Y1, Y2), compose_safe_poly(scene), "numpy")


def target_phi(y, scene=SCENE):
    return float(target_func(scene)(y[0], y[1]))


def safe_psi(y, scene=SCENE):
    return float(safe_func(scene)(y[0], y[1]))


def validate_scene(scene=SCENE, n_boundary=72):
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
```

Also update `target` / `obstacles` in `SCENE` if `validate_scene` rejects the defaults; the validated defaults are: `target=dict(center=(0.34, 0.13), radii=(0.06, 0.05))`, `obstacles=[dict(center=(0.47, 0.10), radius=0.05), dict(center=(0.45, -0.12), radius=0.05)]`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd mujoco_sim && /home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 -m pytest tests/test_scene_def.py -v`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add mujoco_sim/scene_def.py mujoco_sim/tests/test_scene_def.py
git commit -m "feat(scene): single-polynomial safe/target composition + scene validation"
```

---

### Task 2: Render the elliptical workspace in the MuJoCo scene

**Files:**
- Modify: `mujoco_sim/scene_def.py` (`generate_mjcf`: draw the workspace as an ellipse outline instead of a box; use `semi`)
- Test: `mujoco_sim/tests/test_scene_mjcf.py`

**Interfaces:**
- Consumes: `SCENE["workspace"]["semi"]` (Task 1).
- Produces: `generate_mjcf(scene=SCENE, path=None) -> path` writing a loadable `scene_reachavoid.xml`.

- [ ] **Step 1: Write the failing test**

```python
# mujoco_sim/tests/test_scene_mjcf.py
import mujoco
import scene_def as S


def test_mjcf_loads_and_has_geoms():
    path = S.generate_mjcf(S.SCENE)
    m = mujoco.MjModel.from_xml_path(path)
    names = {mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_GEOM, i) for i in range(m.ngeom)}
    assert "target" in names
    assert "obstacle0" in names
    assert "ws_ring" in names  # elliptical workspace outline
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd mujoco_sim && /home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 -m pytest tests/test_scene_mjcf.py -v`
Expected: FAIL (`ws_ring` geom absent; current code emits `ws_n/ws_s/ws_e/ws_w` box bars).

- [ ] **Step 3: Write minimal implementation**

In `generate_mjcf`, replace the 4 box-bar workspace geoms and the translucent box with an elliptical workspace: a flat translucent ellipsoid region + a thin ellipsoid "ring" marker. Replace the workspace block with:

```python
    ax, ay = ws["semi"]
    # translucent workspace region (flat ellipsoid)
    body += _geom(f'<geom name="taskplane" type="ellipsoid" pos="{cx} {cy} {z0}" '
                  f'size="{ax} {ay} 0.002" contype="0" conaffinity="0" '
                  f'rgba="0.5 0.55 0.6 0.12"/>')
    # workspace outline (slightly larger flat ellipsoid, low alpha, reads as a ring)
    body += _geom(f'<geom name="ws_ring" type="ellipsoid" pos="{cx} {cy} {z0-0.001}" '
                  f'size="{ax+0.004} {ay+0.004} 0.001" contype="0" conaffinity="0" '
                  f'rgba="0.35 0.4 0.5 0.5"/>')
```

(Change `hx, hy = ws["half"]` to `ax, ay = ws["semi"]`; delete the 4-bar loop.)

- [ ] **Step 4: Run test to verify it passes**

Run: `cd mujoco_sim && /home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 -m pytest tests/test_scene_mjcf.py -v`
Expected: PASS.

- [ ] **Step 5: Regenerate the scene file and commit**

```bash
cd mujoco_sim && /home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 scene_def.py
git add mujoco_sim/scene_def.py mujoco_sim/tests/test_scene_mjcf.py mujoco_sim/models/franka_emika_panda/scene_reachavoid.xml
git commit -m "feat(scene): elliptical workspace rendering"
```

---

### Task 3: `verify_sets.py` — 2D level-set figure + grid consistency check (Milestone 0, part 1)

**Files:**
- Create: `mujoco_sim/verify_sets.py`
- Test: `mujoco_sim/tests/test_verify_sets.py`

**Interfaces:**
- Consumes: `scene_def.{SCENE, validate_scene, safe_func, target_func}` (Task 1).
- Produces:
  - `set_consistency(scene=SCENE, n=120) -> dict` with keys `agree` (bool), `mismatch_frac` (float).
  - `plot_sets(scene=SCENE, out=...) -> out` saving `figures/milestone0/sets.png`.
  - CLI `python verify_sets.py` runs validate + consistency + plot, prints PASS/FAIL.

- [ ] **Step 1: Write the failing test**

```python
# mujoco_sim/tests/test_verify_sets.py
import scene_def as S
import verify_sets as V


def test_set_consistency_passes_for_default_scene():
    r = V.set_consistency(S.SCENE, n=120)
    assert r["agree"] is True
    assert r["mismatch_frac"] == 0.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd mujoco_sim && /home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 -m pytest tests/test_verify_sets.py -v`
Expected: FAIL (`ModuleNotFoundError: verify_sets`).

- [ ] **Step 3: Write minimal implementation**

```python
# mujoco_sim/verify_sets.py
"""Milestone 0, part 1: confirm the composed single-polynomial sets represent the scene."""
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
    ax.set_aspect("equal"); ax.legend(); ax.set_title("safe {ψ≥0} (green) + target {φ≤0} (red)")
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
```

- [ ] **Step 4: Run test + the CLI**

Run: `cd mujoco_sim && /home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 -m pytest tests/test_verify_sets.py -v && /home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 verify_sets.py`
Expected: test PASS; CLI prints `validate: OK ... agree=True ... PASS` and writes `figures/milestone0/sets.png`.

- [ ] **Step 5: Commit**

```bash
git add mujoco_sim/verify_sets.py mujoco_sim/tests/test_verify_sets.py
git commit -m "feat(milestone0): 2D set-composition verification"
```

---

### Task 4: MATLAB toolchain smoke test

**Files:**
- No code change. Documentation/verification step that the existing synthesis toolchain runs.

**Interfaces:**
- Produces: confidence that SOSTOOLS + MOSEK are on the MATLAB path before writing the new example.

- [ ] **Step 1: Run the existing Dubins synthesis in MATLAB**

From the repo's `matlab/` directory in MATLAB:
```matlab
example_dubins_car
```
Expected: completes without error and writes a fresh `controllers/sop_bounded_control_dubins_car_*.py`. If it errors on `sosprogram`/`mosek`, fix the MATLAB path (add SOSTOOLS, configure MOSEK) before proceeding — Tasks 5-7 depend on it.

- [ ] **Step 2: Confirm the export loads in Python**

Run: `cd /home/jianqiang/Downloads/reach_avoid_backstepping_MPC/controllers && /home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 -c "import sop_bounded_control_dubins_car_result as c; print(type(c.u_opt), type(c.certificate_opt), type(c.k1_opt))"`
Expected: prints `<class 'list'> <class 'sympy...'> <class 'list'>`.

- [ ] **Step 3: Commit (note only)**

No file change; record the gate in the commit log of the next task. (Skip if nothing to commit.)

---

### Task 5: `matlab/example_franka_planar.m` — parameterized double-integrator synthesis

**Files:**
- Create: `matlab/example_franka_planar.m`
- Create (generated at run time by Task 6): `matlab/franka_planar_scene.m` (scene-specific symbolic sets + a_max + degrees)

**Interfaces:**
- Consumes: the synthesis + export functions used by `example_dubins_car.m` (`reach_avoid_controller`, `solvesop_bounded_control`, `export_to_python`).
- Produces: on run, `controllers/sop_bounded_control_franka_planar_*.py` with `u_opt` (list[2]), `certificate_opt` (Expr in x1..x4), `k1_opt` (list[2] in y1,y2).

- [ ] **Step 1: Write `example_franka_planar.m` (copy + adapt `example_dubins_car.m`)**

Open `matlab/example_dubins_car.m` and copy it to `matlab/example_franka_planar.m`. Replace ONLY the system-definition block with the planar double integrator and load the scene-specific sets from `franka_planar_scene.m`:

```matlab
% --- planar task-space double integrator: x=[px,py,vx,vy], u=[ax,ay] ---
syms x1 x2 x3 x4 real        % px, py, vx, vy
syms y1 y2 real
x_vars_sym = [x1; x2; x3; x4];
y_vars_sym = [y1; y2];
fx_sym = [x3; x4; 0; 0];
gx_sym = [0 0; 0 0; 1 0; 0 1];
hx_sym = [x1; x2];

franka_planar_scene;          % defines safe_set_sym, target_set_sym, a_max, ds, dv, mu_val, samples_num
lb = [-a_max; -a_max];
ub = [ a_max;  a_max];
```

Keep the rest of `example_dubins_car.m` (the calls to `reach_avoid_controller`, `solvesop_bounded_control`, `export_to_python`) unchanged, but set the export name to `franka_planar`. Match the exact variable names and call signatures found in `example_dubins_car.m` (read it first — it is the source of truth for the API).

- [ ] **Step 2: Create a default `franka_planar_scene.m` to test the example standalone**

```matlab
% default scene (overwritten by synth.py for other scenes)
syms y1 y2 real
safe_set_sym = (1 - ((y1-0.46)/0.16)^2 - ((y2-0.0)/0.26)^2) ...
             * ((y1-0.47)^2 + (y2-0.10)^2 - 0.05^2) ...
             * ((y1-0.45)^2 + (y2+0.12)^2 - 0.05^2);
target_set_sym = ((y1-0.34)/0.06)^2 + ((y2-0.13)/0.05)^2 - 1;
a_max = 1.0;
ds = 4; dv = 4; mu_val = 0.1; samples_num = 1000;
```

- [ ] **Step 3: Run the synthesis in MATLAB**

From `matlab/`:
```matlab
example_franka_planar
```
Expected: completes and writes `controllers/sop_bounded_control_franka_planar_*.py`. If SOS is infeasible, reduce target difficulty / increase `samples_num` or degrees `ds,dv`, and record what was needed (this is the per-scene feasibility risk).

- [ ] **Step 4: Confirm the export loads and has the right symbols**

Run: `cd controllers && /home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 -c "import importlib,glob,os; f=sorted(glob.glob('sop_bounded_control_franka_planar*'))[-1][:-3]; c=importlib.import_module(f); print('syms', c.certificate_opt.free_symbols); print('k1', len(c.k1_opt))"`
Expected: certificate free symbols ⊆ {x1,x2,x3,x4}; `k1` has length 2.

- [ ] **Step 5: Commit**

```bash
git add matlab/example_franka_planar.m matlab/franka_planar_scene.m
git commit -m "feat(matlab): planar double-integrator reach-avoid synthesis example"
```

---

### Task 6: `synth.py` — drive per-scene synthesis from Python + cache

**Files:**
- Create: `mujoco_sim/synth.py`
- Modify: `python/matlab_runner.py` only if a thin "run a named script" helper is missing (otherwise import as-is)
- Test: `mujoco_sim/tests/test_synth.py`

**Interfaces:**
- Consumes: `scene_def.{compose_safe_poly, compose_target_poly}` (Task 1); `matlab/example_franka_planar.m` (Task 5).
- Produces:
  - `scene_to_matlab(scene) -> str` (the `franka_planar_scene.m` contents: MATLAB symbolic sets + a_max + degrees, via `sympy` → MATLAB strings).
  - `scene_hash(scene) -> str`.
  - `load_bundle(py_module_path) -> dict(u_opt, certificate_opt, k1_opt)` (SymPy).
  - `synthesize(scene, force=False) -> dict` (writes the scene file, runs MATLAB, returns the loaded bundle; caches on `scene_hash`).

- [ ] **Step 1: Write the failing test (pure parts: matlab-string gen, hashing, bundle loading via a fixture)**

```python
# mujoco_sim/tests/test_synth.py
import os
import textwrap
import scene_def as S
import synth


def test_scene_to_matlab_contains_sets_and_bound():
    src = synth.scene_to_matlab(S.SCENE)
    assert "safe_set_sym" in src and "target_set_sym" in src
    assert "a_max = 1.0" in src
    # the elliptical workspace and an obstacle factor appear
    assert "y1" in src and "y2" in src


def test_scene_hash_stable_and_sensitive():
    h1 = synth.scene_hash(S.SCENE)
    h2 = synth.scene_hash({**S.SCENE, "a_max": 2.0})
    assert h1 == synth.scene_hash(S.SCENE)
    assert h1 != h2


def test_load_bundle_from_fixture(tmp_path):
    mod = tmp_path / "fixture_ctrl.py"
    mod.write_text(textwrap.dedent('''
        from sympy import symbols
        x1, x2, x3, x4 = symbols("x1 x2 x3 x4")
        y1, y2 = symbols("y1 y2")
        u_opt = [x3, x4]
        certificate_opt = 1 - x1**2 - x2**2 - x3**2 - x4**2
        k1_opt = [-y1, -y2]
    '''))
    b = synth.load_bundle(str(mod))
    assert len(b["u_opt"]) == 2 and len(b["k1_opt"]) == 2
    assert b["certificate_opt"].free_symbols
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd mujoco_sim && /home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 -m pytest tests/test_synth.py -v`
Expected: FAIL (`ModuleNotFoundError: synth`).

- [ ] **Step 3: Write minimal implementation**

```python
# mujoco_sim/synth.py
"""Drive the per-scene MATLAB SOS synthesis for the planar double integrator."""
import hashlib
import importlib.util
import json
import os
import sympy as sp
import scene_def as S

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
MATLAB_DIR = os.path.join(REPO, "matlab")
CONTROLLERS = os.path.join(REPO, "controllers")


def scene_to_matlab(scene=S.SCENE, ds=4, dv=4, mu_val=0.1, samples_num=1000):
    safe = sp.octave_code(S.compose_safe_poly(scene))
    targ = sp.octave_code(S.compose_target_poly(scene))
    return (
        "syms y1 y2 real\n"
        f"safe_set_sym = {safe};\n"
        f"target_set_sym = {targ};\n"
        f"a_max = {float(scene['a_max'])};\n"
        f"ds = {ds}; dv = {dv}; mu_val = {mu_val}; samples_num = {samples_num};\n"
    )


def scene_hash(scene=S.SCENE):
    key = json.dumps(scene, sort_keys=True, default=list)
    return hashlib.sha1(key.encode()).hexdigest()[:12]


def load_bundle(py_module_path):
    spec = importlib.util.spec_from_file_location("ctrl_bundle", py_module_path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return dict(u_opt=mod.u_opt, certificate_opt=mod.certificate_opt, k1_opt=mod.k1_opt)


def synthesize(scene=S.SCENE, force=False, ds=4, dv=4):
    S.validate_scene(scene)
    cache = os.path.join(HERE, "synth_cache", scene_hash(scene) + ".py")
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    if os.path.exists(cache) and not force:
        return load_bundle(cache)
    with open(os.path.join(MATLAB_DIR, "franka_planar_scene.m"), "w") as f:
        f.write(scene_to_matlab(scene, ds=ds, dv=dv))
    # run MATLAB (headless) on example_franka_planar; matlab must be on PATH
    import subprocess
    subprocess.run(["matlab", "-batch", "example_franka_planar"], cwd=MATLAB_DIR, check=True)
    import glob
    latest = sorted(glob.glob(os.path.join(CONTROLLERS, "sop_bounded_control_franka_planar*.py")))[-1]
    import shutil; shutil.copy(latest, cache)
    return load_bundle(cache)
```

- [ ] **Step 4: Run the pure-part tests**

Run: `cd mujoco_sim && /home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 -m pytest tests/test_synth.py -v`
Expected: PASS (3 tests). (`synthesize()` itself is exercised in Task 7, which needs MATLAB.)

- [ ] **Step 5: Commit**

```bash
git add mujoco_sim/synth.py mujoco_sim/tests/test_synth.py
git commit -m "feat(synth): python driver for per-scene double-integrator synthesis + cache"
```

---

### Task 7: `verify_field.py` — single-integrator k1 field check (Milestone 0, part 2)

**Files:**
- Create: `mujoco_sim/verify_field.py`
- Test: `mujoco_sim/tests/test_verify_field.py`

**Interfaces:**
- Consumes: `synth.synthesize` (Task 6) → bundle with `k1_opt` (list[2] in `y1,y2`); `scene_def.{safe_func, target_func, SCENE}`.
- Produces:
  - `k1_callable(bundle) -> callable(y)->ndarray(2)` (lambdified `k1_opt`).
  - `field_reach_avoid(scene, k1, n_starts=40, t_max=20.0, dt=0.02) -> dict(success_frac, left_safe_frac)` — integrate `ẏ=k1(y)` from safe starts; success = enters target, failure = leaves safe set.
  - `plot_field(scene, k1, out=...)` saving `figures/milestone0/k1_field.png`.
  - CLI `python verify_field.py` → synthesize default scene, run check, plot, print PASS/FAIL.

- [ ] **Step 1: Write the failing test (use an analytic stand-in field so the test is MATLAB-free)**

```python
# mujoco_sim/tests/test_verify_field.py
import numpy as np
import scene_def as S
import verify_field as F


def test_field_reach_avoid_on_gradient_field():
    # analytic reach-avoid-ish field: descend target gradient, repelled by obstacles
    sc = S.SCENE
    tc = np.array(sc["target"]["center"])
    def k1(y):
        v = tc - y                      # pull toward target
        for ob in sc["obstacles"]:
            d = y - np.array(ob["center"]); n = np.linalg.norm(d) + 1e-6
            v = v + 0.02 * d / n**3      # push from obstacles
        s = np.linalg.norm(v) + 1e-9
        return v / s * 0.3
    r = F.field_reach_avoid(sc, k1, n_starts=24, t_max=30.0, dt=0.02)
    assert r["success_frac"] > 0.8
    assert r["left_safe_frac"] < 0.2
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd mujoco_sim && /home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 -m pytest tests/test_verify_field.py -v`
Expected: FAIL (`ModuleNotFoundError: verify_field`).

- [ ] **Step 3: Write minimal implementation**

```python
# mujoco_sim/verify_field.py
"""Milestone 0, part 2: confirm the single-integrator k1 field is a reach-avoid flow."""
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
    f0 = sp.lambdify((S.Y1, S.Y2), bundle["k1_opt"][0], "numpy")
    f1 = sp.lambdify((S.Y1, S.Y2), bundle["k1_opt"][1], "numpy")
    return lambda y: np.array([float(f0(y[0], y[1])), float(f1(y[0], y[1]))])


def _safe_starts(scene, n):
    safe = S.safe_func(scene); phi = S.target_func(scene)
    cx, cy = scene["workspace"]["center"]; ax, ay = scene["workspace"]["semi"]
    rng = np.random.default_rng(0); pts = []
    while len(pts) < n:
        p = np.array([rng.uniform(cx - ax, cx + ax), rng.uniform(cy - ay, cy + ay)])
        if safe(p[0], p[1]) > 1e-3 and phi(p[0], p[1]) > 0:   # safe and not already in target
            pts.append(p)
    return np.array(pts)


def field_reach_avoid(scene, k1, n_starts=40, t_max=20.0, dt=0.02):
    safe = S.safe_func(scene); phi = S.target_func(scene)
    reached = left = 0
    for p0 in _safe_starts(scene, n_starts):
        p = p0.copy(); ok = True
        for _ in range(int(t_max / dt)):
            if safe(p[0], p[1]) < 0:
                left += 1; ok = False; break
            if phi(p[0], p[1]) <= 0:
                reached += 1; ok = False; break
            p = p + dt * k1(p)
        # ran out of time without reaching counts as neither reached nor left
    n = n_starts
    return dict(success_frac=reached / n, left_safe_frac=left / n)


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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd mujoco_sim && /home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 -m pytest tests/test_verify_field.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add mujoco_sim/verify_field.py mujoco_sim/tests/test_verify_field.py
git commit -m "feat(milestone0): single-integrator k1 reach-avoid field check"
```

---

### Task 8: Milestone 0 gate — run the real synthesis end-to-end and decide

**Files:**
- No new code. Runs the real MATLAB synthesis and the two verification scripts on the default scene; records the result.

**Interfaces:**
- Consumes: Tasks 1-7.

- [ ] **Step 1: Run the full Milestone-0 verification (needs MATLAB)**

```bash
cd /home/jianqiang/Downloads/reach_avoid_backstepping_MPC/mujoco_sim
/home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 verify_sets.py
/home/jianqiang/miniconda3/envs/rab_mpc/bin/python3 verify_field.py   # triggers synth.synthesize (MATLAB)
```
Expected: `verify_sets.py` prints `PASS`; `verify_field.py` synthesises (MATLAB runs once), prints `success_frac` high and `left_safe_frac` ~0 → `PASS`, and writes `figures/milestone0/{sets.png,k1_field.png}`.

- [ ] **Step 2: Inspect the two figures**

Open `mujoco_sim/figures/milestone0/sets.png` and `k1_field.png`. Confirm by eye: safe region = workspace ellipse minus the two obstacle disks; target ellipse inside; k1 streamlines flow into the target and bend around obstacles without crossing the green boundary.

- [ ] **Step 3: Record the gate outcome**

- If PASS: write a one-line note in the plan/commit and proceed to the Phase-2 plan (controllers + MuJoCo closed loop).
- If FAIL (SOS infeasible, or k1 field leaves the safe set): STOP. Capture which scene/parameters failed; options are (a) adjust degrees `ds,dv` / `samples_num`, (b) simplify the scene (fewer/positioned obstacles), (c) revisit whether the DI k1 field needs a different reach-avoid term. Do NOT start Phase 2 until this passes.

- [ ] **Step 4: Commit the milestone artifacts**

```bash
git add mujoco_sim/figures/milestone0 mujoco_sim/synth_cache
git commit -m "chore(milestone0): record passing set + k1-field verification"
```

---

## Self-review

- **Spec coverage:** Milestone 0 (§5.1) → Tasks 1,3,5,6,7,8. Single-polynomial composition (§5) → Task 1. Elliptical workspace correctness (§5) → Tasks 1,2. Parameterized per-scene synthesis (§6) → Tasks 5,6. Validation strategy "Milestone 0 gate" (§10) → Task 8. Plant/controllers/sim_loop/metrics/sweeps (§4,§7,§8 remainder) → **deferred to the Phase-2 plan** (post-gate), by design.
- **Placeholder scan:** no TBD/TODO; every code/test step has concrete code. Scene default numbers are concrete and validated by Task 1's tests.
- **Type consistency:** `compose_safe_poly`/`compose_target_poly`/`safe_func`/`target_func`/`validate_scene`/`scene_hash`/`load_bundle`/`synthesize`/`k1_callable`/`field_reach_avoid` names are consistent across Tasks 1,3,6,7. `SCENE["workspace"]["semi"]` used consistently after Task 1.

## Out of scope (Phase-2 plan, after the gate passes)

`controllers.py` (RA / RA-MPC / Vanilla) · `sim_loop.py` (MuJoCo closed loop using the built `FrankaPlanarPlant`) · `metrics.py` · `viz.py` (closed-loop trajectory + video via `view_sim.py`) · `sweep.py` (initial-condition / mismatch sweeps) · headline comparison (decision §9-E: ~20 initial conditions, RA-MPC reach-avoid success vs Vanilla failures on obstacle scenes).
