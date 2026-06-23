"""Visualize the Franka planar plant in MuJoCo.

Two modes:
  live   : interactive window (orbit/zoom, watch in real time). Needs a display.
  record : offscreen render to an mp4 (works headless via EGL + ffmpeg).

Demo controller = a task-space PD circle tracker that outputs an EE acceleration u
(the same interface the reach-avoid law will use). The circle is slow/small so the
commanded acceleration stays well inside the dexterous envelope in all directions.

Examples:
  python view_sim.py --mode live
  python view_sim.py --mode record --out figures/franka_circle.mp4 --seconds 8
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCENES = {
    "plain": os.path.join(HERE, "models", "franka_emika_panda", "scene.xml"),
    "reachavoid": os.path.join(HERE, "models", "franka_emika_panda", "scene_reachavoid.xml"),
}


def circle_ref(t, c, r, T):
    w = 2 * np.pi / T
    p = np.array([c[0] + r * np.cos(w * t), c[1] + r * np.sin(w * t)])
    v = np.array([-r * w * np.sin(w * t), r * w * np.cos(w * t)])
    a = np.array([-r * w * w * np.cos(w * t), -r * w * w * np.sin(w * t)])
    return p, v, a


class CircleTracker:
    """Outputs task-space acceleration u to trace a horizontal circle (demo only)."""

    def __init__(self, plant, center=None, r=0.10, T=5.0, kp=60.0, kd=16.0, a_max=0.6):
        p0, _ = plant.ee_full()
        self.c = np.asarray(center, float) if center is not None else p0[:2].copy()
        self.r, self.T, self.kp, self.kd, self.a_max = r, T, kp, kd, a_max

    def __call__(self, t, x):
        p, v = x[:2], x[2:]
        pr, vr, ar = circle_ref(t, self.c, self.r, self.T)
        u = ar + self.kp * (pr - p) + self.kd * (vr - v)
        return np.clip(u, -self.a_max, self.a_max)


def make_camera(mujoco):
    cam = mujoco.MjvCamera()
    cam.azimuth = 135.0
    cam.elevation = -20.0
    cam.distance = 2.2
    cam.lookat[:] = [0.3, 0.0, 0.5]
    return cam


def run_live(scene_path, center, seconds):
    # force the on-screen GL backend (avoid egl/osmesa being inherited)
    os.environ["MUJOCO_GL"] = "glfw"
    import mujoco
    import mujoco.viewer
    from plant import FrankaPlanarPlant

    plant = FrankaPlanarPlant(scene_path, ctrl_hz=100.0)
    ctrl = CircleTracker(plant, center=center)
    plant.reset_home()
    # seconds <= 0  -> run until the user closes the window
    deadline = None if seconds is None or seconds <= 0 else (time.time() + seconds)
    with mujoco.viewer.launch_passive(plant.m, plant.d) as viewer:
        cam = make_camera(mujoco)
        viewer.cam.azimuth, viewer.cam.elevation = cam.azimuth, cam.elevation
        viewer.cam.distance, viewer.cam.lookat[:] = cam.distance, cam.lookat
        while viewer.is_running() and (deadline is None or time.time() < deadline):
            tic = time.time()
            t = plant.d.time
            plant.step(ctrl(t, plant.state()))
            viewer.sync()
            dt = plant.ctrl_dt - (time.time() - tic)
            if dt > 0:
                time.sleep(dt)
    print("live viewer closed.")
    sys.stdout.flush()
    # MuJoCo's passive viewer can segfault during GL/GLFW teardown at interpreter
    # shutdown on some Linux/conda setups; bypass destructors with a hard exit.
    os._exit(0)


def run_record(scene_path, center, out, seconds, fps, w, h):
    seconds = 8.0 if not seconds or seconds <= 0 else seconds
    os.environ.setdefault("MUJOCO_GL", "egl")
    import subprocess
    import mujoco
    from plant import FrankaPlanarPlant

    plant = FrankaPlanarPlant(scene_path, ctrl_hz=100.0)
    ctrl = CircleTracker(plant, center=center)
    plant.reset_home()
    cam = make_camera(mujoco)
    renderer = mujoco.Renderer(plant.m, h, w)

    out = out if os.path.isabs(out) else os.path.join(HERE, out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cmd = ["ffmpeg", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24",
           "-s", f"{w}x{h}", "-r", str(fps), "-i", "-", "-an",
           "-vcodec", "libx264", "-pix_fmt", "yuv420p", out]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.DEVNULL)

    n_steps = int(seconds / plant.ctrl_dt)
    every = max(1, int(round(1.0 / (fps * plant.ctrl_dt))))   # render every k control steps
    nframes = 0
    for k in range(n_steps):
        t = plant.d.time
        plant.step(ctrl(t, plant.state()))
        if k % every == 0:
            renderer.update_scene(plant.d, camera=cam)
            proc.stdin.write(renderer.render().tobytes())
            nframes += 1
    proc.stdin.close()
    proc.wait()
    renderer.close()
    print(f"wrote {nframes} frames -> {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["live", "record"], default="live")
    ap.add_argument("--scene", choices=list(SCENES), default="reachavoid")
    ap.add_argument("--out", default="figures/franka_circle.mp4")
    ap.add_argument("--seconds", type=float, default=0.0,
                    help="live: 0 = run until window closed; record: 0 = 8s")
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--width", type=int, default=640)
    ap.add_argument("--height", type=int, default=480)
    args = ap.parse_args()

    scene_path = SCENES[args.scene]
    center = None
    if args.scene == "reachavoid":
        import scene_def
        center = scene_def.SCENE["workspace"]["center"]   # trace the demo circle inside the workspace

    if args.mode == "live":
        run_live(scene_path, center, args.seconds)
    else:
        run_record(scene_path, center, args.out, args.seconds, args.fps, args.width, args.height)


if __name__ == "__main__":
    main()
