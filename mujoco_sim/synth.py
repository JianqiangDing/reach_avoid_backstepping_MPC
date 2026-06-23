"""Drive the per-scene MATLAB SOS synthesis for the planar double integrator.

Generates a scene-specific MATLAB file (symbolic safe/target sets + a_max + degrees),
runs `example_franka_planar.m`, and loads the exported SymPy controller bundle
(u_opt / certificate_opt / k1_opt). Results are cached on a scene hash.
"""

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
    """MATLAB source defining the scene-specific symbolic sets + parameters."""
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
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return dict(u_opt=mod.u_opt, certificate_opt=mod.certificate_opt, k1_opt=mod.k1_opt)


def synthesize(scene=S.SCENE, force=False, ds=4, dv=4):
    """Run the per-scene MATLAB synthesis (cached) and return the controller bundle."""
    S.validate_scene(scene)
    cache = os.path.join(HERE, "synth_cache", scene_hash(scene) + ".py")
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    if os.path.exists(cache) and not force:
        return load_bundle(cache)
    with open(os.path.join(MATLAB_DIR, "franka_planar_scene.m"), "w") as f:
        f.write(scene_to_matlab(scene, ds=ds, dv=dv))
    import subprocess
    subprocess.run(["matlab", "-batch", "example_franka_planar"], cwd=MATLAB_DIR, check=True)
    import glob
    exported = sorted(glob.glob(os.path.join(CONTROLLERS, "sop_bounded_control_franka_planar*.py")))
    if not exported:
        raise RuntimeError("synthesis produced no controller export (SOS infeasible?)")
    import shutil
    shutil.copy(exported[-1], cache)
    return load_bundle(cache)
