"""Abstraction-fidelity probe (P0.5).

Question: when we command a planar EE acceleration u to the Franka-under-OSC, how
faithfully does the realized planar motion match the ideal double integrator
    p_dot = v,  v_dot = u ?

We run open-loop command profiles, compare against the ideal DI rollout under the
SAME u(t), and report gain/isotropy, effective lag/bandwidth, DI tracking RMSE,
height/orientation drift, and torque-saturation margin -> verdict.

Run:  python mujoco_sim/probe_abstraction.py
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plant import FrankaPlanarPlant

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL = os.path.join(HERE, "models", "franka_emika_panda", "panda.xml")
FIGDIR = os.path.join(HERE, "figures", "probe")
os.makedirs(FIGDIR, exist_ok=True)


# ---- ideal double integrator (exact ZOH under piecewise-constant u) ----------
def ideal_di(u_seq, dt, p0, v0):
    p, v = np.array(p0, float), np.array(v0, float)
    P, V = [p.copy()], [v.copy()]
    for u in u_seq:
        p = p + v * dt + 0.5 * u * dt * dt
        v = v + u * dt
        P.append(p.copy()); V.append(v.copy())
    return np.array(P), np.array(V)


def run_signal(plant: FrankaPlanarPlant, u_seq):
    """Drive plant with u_seq (ZOH per step). Return actual + ideal planar rollouts."""
    plant.reset_home()
    x0 = plant.state()
    p0, v0 = x0[:2], x0[2:]
    P, V, Z, ORI, SAT = [p0.copy()], [v0.copy()], [0.0], [0.0], [0.0]
    for u in u_seq:
        r = plant.step(u)
        P.append(r["p"][:2].copy()); V.append(r["v"][:2].copy())
        Z.append(r["z_drift"]); ORI.append(r["ori_drift"]); SAT.append(r["sat_frac"])
    Pi, Vi = ideal_di(u_seq, plant.ctrl_dt, p0, v0)
    return dict(P=np.array(P), V=np.array(V), Z=np.array(Z), ORI=np.array(ORI),
                SAT=np.array(SAT), Pi=Pi, Vi=Vi, dt=plant.ctrl_dt)


# ---- individual tests --------------------------------------------------------
def test_gain_isotropy(plant, a0=1.0, dur=0.6):
    n = int(dur / plant.ctrl_dt)
    dirs = np.linspace(0, 2 * np.pi, 8, endpoint=False)
    rows = []
    for th in dirs:
        uhat = np.array([np.cos(th), np.sin(th)])
        u_seq = [a0 * uhat] * n
        out = run_signal(plant, u_seq)
        V = out["V"]
        a_est = np.diff(V, axis=0) / plant.ctrl_dt          # realized accel per step
        a_mean = a_est[2:].mean(axis=0)                      # skip initial transient
        gain = float(a_mean @ uhat / a0)
        cross = float(np.linalg.norm(a_mean - (a_mean @ uhat) * uhat) / a0)
        rows.append((float(th), gain, cross, float(out["Z"].max()),
                     float(out["ORI"].max()), float(out["SAT"].max())))
    return rows


def test_chirp(plant, a0=1.0, f0=0.2, f1=3.0, dur=5.0, axis=1):
    """Chirp along a dexterous axis (default +y) to measure bandwidth/lag in-envelope."""
    n = int(dur / plant.ctrl_dt)
    t = np.arange(n) * plant.ctrl_dt
    f = f0 + (f1 - f0) * t / dur
    ua = a0 * np.sin(2 * np.pi * f * t)
    u_seq = [(np.array([0.0, uai]) if axis == 1 else np.array([uai, 0.0])) for uai in ua]
    out = run_signal(plant, u_seq)
    # effective velocity-gain vs frequency in a few bands (actual vs ideal amplitude)
    Vx, Vix = out["V"][1:, axis], out["Vi"][1:, axis]
    bands = [(0.2, 0.6), (0.6, 1.2), (1.2, 2.0), (2.0, 3.0)]
    bandinfo = []
    for lo, hi in bands:
        mask = (f >= lo) & (f < hi)
        if mask.sum() < 5:
            continue
        amp_act = np.sqrt(2) * Vx[mask].std()
        amp_idl = np.sqrt(2) * Vix[mask].std()
        g = amp_act / amp_idl if amp_idl > 1e-9 else np.nan
        bandinfo.append((0.5 * (lo + hi), float(g)))
    return out, f, bandinfo


def test_maneuver(plant, a0=0.6, dur=1.5):
    """A curved maneuver kept inside the dexterous region (gentle -x bias + y oscillation)."""
    n = int(dur / plant.ctrl_dt)
    t = np.arange(n) * plant.ctrl_dt
    ux = -0.3 * a0 * (1 - np.cos(2 * np.pi * 0.5 * t))   # always <= 0 (toward -x)
    uy = a0 * np.sin(2 * np.pi * 0.5 * t)
    u_seq = [np.array([ux[k], uy[k]]) for k in range(n)]
    out = run_signal(plant, u_seq)
    pos_rmse = float(np.sqrt(((out["P"] - out["Pi"]) ** 2).sum(axis=1).mean()))
    vel_rmse = float(np.sqrt(((out["V"] - out["Vi"]) ** 2).sum(axis=1).mean()))
    return out, pos_rmse, vel_rmse


# ---- driver ------------------------------------------------------------------
def main():
    print(f"model: {MODEL}")
    deg = 180.0 / np.pi

    # --- self-check: command zero accel, expect EE to hold ---
    plant = FrankaPlanarPlant(MODEL, ctrl_hz=100.0)
    print(f"home TCP z0={plant.z0:.4f}  ctrl_hz={plant.ctrl_hz}  n_sub={plant.n_sub}")
    out0 = run_signal(plant, [np.zeros(2)] * int(1.0 / plant.ctrl_dt))
    drift_xy = float(np.linalg.norm(out0["P"][-1] - out0["P"][0]))
    print("\n[SELF-CHECK u=0, 1s]")
    print(f"  planar drift  = {drift_xy*1000:.3f} mm   (want < 2 mm)")
    print(f"  z drift max   = {out0['Z'].max()*1000:.3f} mm")
    print(f"  ori drift max = {out0['ORI'].max()*deg:.3f} deg")

    # manipulability at home (planar position Jacobian)
    jacp, _ = plant._tcp_jac()
    Jpl = jacp[:2, :]
    manip = float(np.sqrt(max(np.linalg.det(Jpl @ Jpl.T), 0.0)))
    print(f"planar manipulability sqrt(det(J J^T)) at home = {manip:.3f}")

    # --- gain / isotropy ---
    rows = test_gain_isotropy(plant, a0=1.0, dur=0.6)
    gains = np.array([r[1] for r in rows]); cross = np.array([r[2] for r in rows])
    sat = np.array([r[5] for r in rows])
    dex = sat < 0.9                                  # dexterous (non-saturating) directions
    print("\n[GAIN / ISOTROPY  a0=1.0 m/s^2, 8 dirs]   (* = kinematic-boundary dir)")
    print("   dir(deg)  gain   cross-leak  zdrift(mm)  ori(deg)  sat")
    for r in rows:
        flag = " " if r[5] < 0.9 else "*"
        print(f" {flag} {r[0]*deg:6.0f}   {r[1]:.3f}   {r[2]:.3f}      "
              f"{r[3]*1000:6.2f}     {r[4]*deg:5.2f}   {r[5]:.2f}")
    print(f"  DEXTEROUS dirs ({dex.sum()}/8): gain mean={gains[dex].mean():.3f} "
          f"min={gains[dex].min():.3f} max={gains[dex].max():.3f}  cross-leak max={cross[dex].max():.3f}")
    bnd = np.where(~dex)[0]
    if len(bnd):
        print(f"  BOUNDARY dirs (saturate at a0=1): "
              + ", ".join(f"{rows[i][0]*deg:.0f}deg(gain {rows[i][1]:.2f})" for i in bnd))

    # --- chirp / bandwidth ---
    out_c, f, bandinfo = test_chirp(plant, a0=1.0)
    print("\n[CHIRP 0.2->3 Hz]  velocity-gain (actual/ideal) per band:")
    for fc, g in bandinfo:
        print(f"   ~{fc:.1f} Hz : {g:.3f}")

    # --- maneuver / DI tracking ---
    out_m, pos_rmse, vel_rmse = test_maneuver(plant, a0=1.0, dur=2.0)
    print("\n[MANEUVER 2s curved]")
    print(f"  position RMSE vs ideal DI = {pos_rmse*1000:.3f} mm")
    print(f"  velocity RMSE vs ideal DI = {vel_rmse*1000:.3f} mm/s")
    print(f"  z drift max = {out_m['Z'].max()*1000:.3f} mm   "
          f"ori max = {out_m['ORI'].max()*deg:.3f} deg   sat max = {out_m['SAT'].max():.2f}")

    # --- control-rate sweep (in-envelope: +y, a0=1.0) ---
    print("\n[CONTROL-RATE SWEEP  const accel a0=1.0, dir=+y (dexterous)]")
    for hz in (50, 100, 200, 500):
        p = FrankaPlanarPlant(MODEL, ctrl_hz=float(hz))
        n = int(0.6 / p.ctrl_dt)
        out = run_signal(p, [np.array([0.0, 1.0])] * n)
        a_est = np.diff(out["V"], axis=0)[2:, 1] / p.ctrl_dt
        print(f"   {hz:4d} Hz : realized ay={a_est.mean():.3f} (cmd 1.0)  "
              f"sat={out['SAT'].max():.2f}  zdrift={out['Z'].max()*1000:.2f}mm")

    # --- plots ---
    _plots(rows, out_c, f, out_m)

    # --- verdict (evaluated WITHIN the dexterous envelope) ---
    print("\n========== VERDICT (within dexterous envelope) ==========")
    ok_gain = (gains[dex].min() >= 0.85) and (gains[dex].max() <= 1.15)
    ok_iso = cross[dex].max() < 0.15
    ok_track = pos_rmse < 0.005
    ok_z = out_m["Z"].max() < 0.002
    ok_ori = out_m["ORI"].max() < 0.0175
    ok_sat = out_m["SAT"].max() < 1.0
    for name, ok in [("dexterous-dir gain in [0.85,1.15]", ok_gain),
                     ("dexterous-dir isotropy (leak<0.15)", ok_iso),
                     ("in-envelope DI pos RMSE<5mm", ok_track), ("z drift<2mm", ok_z),
                     ("ori drift<1deg", ok_ori), ("no saturation on in-envelope maneuver", ok_sat)]:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    overall = all([ok_gain, ok_iso, ok_track, ok_z, ok_ori, ok_sat])
    print(f"\n  OVERALL: {'PASS - double integrator is faithful inside the dexterous workspace' if overall else 'NEEDS ATTENTION'}")
    print(f"  NOTE: {int((~dex).sum())}/8 directions hit the kinematic boundary (arm extension);")
    print(f"        site the task workspace in the dexterous region and bound a_max accordingly.")
    print(f"\n  figures -> {FIGDIR}")


def _plots(rows, out_c, f, out_m):
    deg = 180.0 / np.pi
    # gain polar
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111, projection="polar")
    th = np.array([r[0] for r in rows] + [rows[0][0]])
    g = np.array([r[1] for r in rows] + [rows[0][1]])
    ax.plot(th, g, "o-"); ax.plot(th, np.ones_like(th), "k--", alpha=0.4)
    ax.set_title("realized accel gain vs direction"); ax.set_ylim(0, 1.2)
    fig.savefig(os.path.join(FIGDIR, "gain_polar.png"), dpi=130, bbox_inches="tight"); plt.close(fig)

    # chirp velocity actual vs ideal
    fig, ax = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    t = np.arange(len(out_c["V"]) - 1) * out_c["dt"]
    ax[0].plot(t, out_c["Vi"][1:, 0], label="ideal vx", lw=1)
    ax[0].plot(t, out_c["V"][1:, 0], label="actual vx", lw=1, alpha=0.8)
    ax[0].legend(); ax[0].set_ylabel("vx (m/s)"); ax[0].set_title("chirp: ideal vs actual planar velocity")
    ax[1].plot(t, f); ax[1].set_ylabel("freq (Hz)"); ax[1].set_xlabel("t (s)")
    fig.savefig(os.path.join(FIGDIR, "chirp.png"), dpi=130, bbox_inches="tight"); plt.close(fig)

    # maneuver path
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].plot(out_m["Pi"][:, 0], out_m["Pi"][:, 1], label="ideal DI", lw=2)
    ax[0].plot(out_m["P"][:, 0], out_m["P"][:, 1], "--", label="actual Franka", lw=2)
    ax[0].set_aspect("equal"); ax[0].legend(); ax[0].set_title("maneuver path (x-y)")
    ax[0].set_xlabel("x (m)"); ax[0].set_ylabel("y (m)")
    tm = np.arange(len(out_m["Z"])) * out_m["dt"]
    ax[1].plot(tm, out_m["Z"] * 1000, label="z drift (mm)")
    ax[1].plot(tm, out_m["ORI"] * deg, label="ori drift (deg)")
    ax[1].plot(tm, out_m["SAT"] * 100, label="torque sat (%)")
    ax[1].legend(); ax[1].set_xlabel("t (s)"); ax[1].set_title("regulation drift & saturation")
    fig.savefig(os.path.join(FIGDIR, "maneuver.png"), dpi=130, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    main()
