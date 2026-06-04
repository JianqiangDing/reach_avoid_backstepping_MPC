# Design: Update the 3 manipulator example notebooks (unified feasibility)

**Date:** 2026-06-04
**Topic:** Sub-project 2 of the manipulator parallel — apply the same unified feasibility
criterion (and fresh controllers) to the three `example_manipulator_*` notebooks that was
applied to the dubins ones.

## Unified feasibility (all three, reach window `T_max = 10 s`)

A start is **feasible** iff the method safely drives the acrobot end-effector into the
target region within the window:

- **Vanilla MPC**: `solve_mpc` with terminal `target_ca(x_N) <= 0` is solvable (reach the
  target within its `N=20`, `dt=0.01` → 0.2 s horizon; equivalent to the 10 s receding loop
  since the target terminal forces each step to reach in-horizon). **No code change beyond a
  re-run** — its feasibility is already target-terminal solvability.
- **RA-MPC**: the receding closed loop (cert terminal `cert_ca_vec(x_N) >= 0`) reaches the
  target safely within 10 s, using **continuous-ZOH** stepping (plant integrated continuously
  with the ZOH torque; the MPC plans with the coarse `dt` model). Use the existing
  `solve_mpc_fast` (pre-compiled NLP) for tractability.
- **Unconstrained**: the closed loop under the synthesized controller `u_opt` reaches the
  target safely within 10 s.

## Changes per notebook

### `example_manipulator_vanilla_mpc.ipynb`
- No controller import (self-contained). Feasibility already = `solve_mpc(target)` solvable.
- **Re-run only** to refresh outputs (+ ensure `data/` exists for any savez).

### `example_manipulator_reach_avoid_mpc.ipynb`
- **Repoint import**: `sop_bounded_control_acrobot_result_20260317_222858` →
  `sop_bounded_control_acrobot_result`.
- **Feasibility (cell 3)**: replace the single-solve `_check_one` with a continuous-ZOH
  receding check (cert terminal, `solve_mpc_fast`, `T_max=10 s`, plant stepped by fine RK4
  substeps of the true dynamics under ZOH torque; safe along the way; reaches target).
- Re-run.

### `example_manipulator_unconstrained_reach_avoid.ipynb`
- **Switch controller**: drop `from k1_acrobot_cdc2026 import (...)` + the in-notebook
  k1→backstepping construction; import `u_opt, certificate_opt, k1_opt` from the freshly
  synthesized `sop_bounded_control_acrobot_unconstrained` (parallel to dubins). Build
  `cl_dyn = f_sym + g_mat @ u_opt`.
- **Feasibility (cells 3–4)**: replace cert-sign (`cert_vals >= 0`) with closed-loop reaches
  target safely within `T_max=10 s` (continuous ODE under `u_opt`).
- Re-run.

## Reused manipulator settings (from the notebooks, unchanged)
- 4-state `[x1,x2,x3,x4] = [q1,q2,dq1,dq2]`; output = end-effector `[4cos(q1+q2)+4cos(q1),
  4sin(...)+...]` (l1=l2=4); full M-C-G dynamics.
- MPC: `dt_mpc=0.01`, `N_hor=20`, `u_max=500`; target/safe sets in y-space; `TARGET_POS≈[5.8,1.9]`.
- Sampling: y-space sample → `inverse_kinematics` → joint space; `N_candidates_max=100`, seed 42.

## Error handling
- Missing `data/` dir for savez → create it (as for dubins).
- Controller file missing → clear assert.

## Testing / verification
Run each in place via nbclient (rab_mpc kernel); verify 0 cell errors, the feasibility
counts print, figures embed, and (RA-MPC/uncon) the closed loops reach the target. The
manipulator RA-MPC continuous-ZOH receding feasibility may be slow (dt=0.01); use
`solve_mpc_fast` and the 100-candidate cap.

## Out of scope
- Sub-project 3 (`make_paper_figs_manipulator.ipynb`) — separate.
- No MATLAB / controller-synthesis changes (done in sub-project 1).
