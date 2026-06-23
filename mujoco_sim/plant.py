"""Franka Panda planar task-space plant for the abstraction-fidelity probe.

The reach-avoid controller is designed for an idealized *planar double integrator*
in task space:  p_dot = v,  v_dot = u,  with u = (u_x, u_y) a commanded EE
acceleration in a horizontal plane at fixed height z0.

This module wraps the full 7-DoF Franka (MuJoCo Menagerie) behind that interface.
The robot's "built-in controller" is a standard Operational Space Controller (OSC):
  - the outer reach-avoid command u_xy is held Zero-Order-Hold over one control
    period (ctrl_hz),
  - the inner OSC recomputes joint torques at every physics step (fast inner loop),
    realizing the planar acceleration while regulating EE height z0, EE orientation,
    and a nullspace posture (the redundant 7th DoF).

The probe then measures how faithfully the realized planar EE motion matches the
ideal double integrator under the same u(t).

Coordinates: planar task = world (x, y); regulated height = world z.
"""

from __future__ import annotations

import numpy as np
import mujoco


ARM_DOF = 7  # first 7 joints are the arm; 7,8 are the gripper fingers


class FrankaPlanarPlant:
    def __init__(
        self,
        model_path: str,
        ctrl_hz: float = 100.0,
        tcp_offset: float = 0.0,           # TCP point below the hand origin (0 = hand origin, no moment arm)
        # OSC task gains (inner loop)
        kp_xy: float = 0.0, kd_xy: float = 0.0,   # x,y are acceleration-commanded (no PD)
        kp_z: float = 900.0, kd_z: float = 60.0,
        kp_rot: float = 600.0, kd_rot: float = 50.0,
        kp_null: float = 25.0, kd_null: float = 10.0,
        meas_noise_pos: float = 0.0,
        meas_noise_vel: float = 0.0,
    ):
        self.m = mujoco.MjModel.from_xml_path(model_path)
        self.d = mujoco.MjData(self.m)
        self.ctrl_hz = float(ctrl_hz)
        self.ctrl_dt = 1.0 / self.ctrl_hz
        self.n_sub = max(1, int(round(self.ctrl_dt / self.m.opt.timestep)))
        self.tcp_offset = float(tcp_offset)

        self.kp_z, self.kd_z = kp_z, kd_z
        self.kp_rot, self.kd_rot = kp_rot, kd_rot
        self.kp_null, self.kd_null = kp_null, kd_null
        self.meas_noise_pos = meas_noise_pos
        self.meas_noise_vel = meas_noise_vel
        self._rng = np.random.default_rng(0)

        self.hand_bid = mujoco.mj_name2id(self.m, mujoco.mjtObj.mjOBJ_BODY, "hand")
        assert self.hand_bid >= 0, "hand body not found"

        # Make the 7 arm position-servo actuators inert; we inject torque via qfrc_applied.
        for i in range(ARM_DOF):
            self.m.actuator_gainprm[i, :] = 0.0
            self.m.actuator_biasprm[i, :] = 0.0
        # Arm torque limits (from the model's forcerange) for saturation reporting/clipping.
        self.tau_limit = np.abs(self.m.actuator_forcerange[:ARM_DOF, 1]).copy()

        self._J_prev = None
        self.reset_home()
        # Regulation targets captured at the home pose.
        self.z0 = float(self._tcp_pos()[2])
        self.R_des = self.d.xmat[self.hand_bid].reshape(3, 3).copy()
        self.q_home = self.d.qpos[:ARM_DOF].copy()

    # ---- state access -------------------------------------------------------
    def reset_home(self):
        mujoco.mj_resetDataKeyframe(self.m, self.d, 0)  # "home" keyframe
        self.d.qvel[:] = 0.0
        self.d.qfrc_applied[:] = 0.0
        self._J_prev = None
        mujoco.mj_forward(self.m, self.d)

    def _tcp_pos(self) -> np.ndarray:
        # TCP = hand origin shifted by tcp_offset along the hand's local -z (toward fingers).
        R = self.d.xmat[self.hand_bid].reshape(3, 3)
        return self.d.xpos[self.hand_bid] + R @ np.array([0.0, 0.0, self.tcp_offset])

    def _tcp_jac(self):
        jacp = np.zeros((3, self.m.nv))
        jacr = np.zeros((3, self.m.nv))
        mujoco.mj_jac(self.m, self.d, jacp, jacr, self._tcp_pos(), self.hand_bid)
        return jacp[:, :ARM_DOF], jacr[:, :ARM_DOF]

    def ee_full(self):
        """Return (p3, v3) of the TCP in world frame."""
        jacp, _ = self._tcp_jac()
        p = self._tcp_pos()
        v = jacp @ self.d.qvel[:ARM_DOF]
        return p, v

    def state(self) -> np.ndarray:
        """Idealized planar state x = [px, py, vx, vy] (with optional measurement noise)."""
        p, v = self.ee_full()
        x = np.array([p[0], p[1], v[0], v[1]])
        if self.meas_noise_pos:
            x[:2] += self._rng.normal(0, self.meas_noise_pos, 2)
        if self.meas_noise_vel:
            x[2:] += self._rng.normal(0, self.meas_noise_vel, 2)
        return x

    # ---- OSC inner loop -----------------------------------------------------
    def _osc_torque(self, u_xy: np.ndarray) -> np.ndarray:
        m, d = self.m, self.d
        jacp, jacr = self._tcp_jac()
        J = np.vstack([jacp, jacr])          # (6, 7)
        qvel = d.qvel[:ARM_DOF]

        # arm mass matrix (7x7)
        full = np.zeros((m.nv, m.nv))
        mujoco.mj_fullM(m, d, full)
        M = full[:ARM_DOF, :ARM_DOF]
        Minv = np.linalg.inv(M)

        # task-space inertia Lambda = (J Minv J^T)^-1  (regularized)
        JMinvJt = J @ Minv @ J.T
        Lambda = np.linalg.inv(JMinvJt + 1e-6 * np.eye(6))

        # desired task acceleration: [ax, ay, az_reg, arot_reg(3)]
        p, v = self.ee_full()
        omega = jacr @ qvel
        a_pos = np.array([u_xy[0], u_xy[1],
                          self.kp_z * (self.z0 - p[2]) - self.kd_z * v[2]])
        a_rot = self.kp_rot * self._ori_err() - self.kd_rot * omega
        a_des = np.concatenate([a_pos, a_rot])

        # J_dot * qvel via finite difference of the task Jacobian (fast inner loop)
        Jdq = np.zeros(6)
        if self._J_prev is not None:
            Jdq = ((J - self._J_prev) / m.opt.timestep) @ qvel
        self._J_prev = J.copy()

        F = Lambda @ (a_des - Jdq)            # operational-space force
        tau = J.T @ F

        # dynamically-consistent nullspace posture (keep redundant DoF near home)
        Jbar = Minv @ J.T @ Lambda            # (7,6)
        N = np.eye(ARM_DOF) - J.T @ Jbar.T     # (7,7)
        tau_null = self.kp_null * (self.q_home - d.qpos[:ARM_DOF]) - self.kd_null * qvel
        tau = tau + N @ tau_null

        # gravity + Coriolis/centrifugal compensation (joint space)
        tau = tau + d.qfrc_bias[:ARM_DOF]
        # passive joint-damping compensation (qfrc_bias excludes passive damping);
        # a real Cartesian-impedance controller compensates joint friction/damping.
        tau = tau + self.m.dof_damping[:ARM_DOF] * qvel
        return np.clip(tau, -self.tau_limit, self.tau_limit)

    def step(self, u_xy) -> dict:
        """Hold u_xy (ZOH) over one control period; recompute OSC torque each physics step."""
        u_xy = np.asarray(u_xy, float)
        sat_frac = 0.0
        for _ in range(self.n_sub):
            tau = self._osc_torque(u_xy)
            self.d.qfrc_applied[:ARM_DOF] = tau
            sat_frac = max(sat_frac, float(np.max(np.abs(tau) / self.tau_limit)))
            mujoco.mj_step(self.m, self.d)
        self.d.qfrc_applied[:ARM_DOF] = 0.0
        p, v = self.ee_full()
        return dict(
            t=float(self.d.time), p=p.copy(), v=v.copy(),
            z_drift=float(p[2] - self.z0),
            ori_drift=float(np.linalg.norm(self._ori_err())),
            sat_frac=sat_frac,
        )

    def _ori_err(self) -> np.ndarray:
        """World-frame orientation error driving R_cur -> R_des (consistent with world jacr)."""
        R = self.d.xmat[self.hand_bid].reshape(3, 3)
        Rd = self.R_des
        return 0.5 * (np.cross(R[:, 0], Rd[:, 0])
                      + np.cross(R[:, 1], Rd[:, 1])
                      + np.cross(R[:, 2], Rd[:, 2]))
