# Design: Manipulator (acrobot) controller-synthesis notebook

**Date:** 2026-06-04
**Topic:** Sub-project 1 of the manipulator parallel — a `notebooks/` notebook that drives
the MATLAB synthesis to generate the acrobot *unconstrained* and *constrained* control laws.
Mirrors `notebooks/synthesize_dubins_controllers.ipynb`
(spec: `2026-06-04-dubins-controller-synthesis-notebook-design.md`).

## Goal

One notebook that, when run, synthesizes both acrobot control laws by calling the existing
MATLAB pipeline via `matlab -batch`, then loads and displays the two resulting control laws.
Prerequisite for sub-projects 2 (update the manipulator example notebooks) and 3
(`make_paper_figs_manipulator.ipynb`).

## Approach

Drive the **existing, unmodified** `matlab/example_manipulator.m` via `matlab -batch`,
reusing `python/matlab_runner.py` (`find_matlab`, `repo_root`, `run_one`). No new `.m` files;
no edits to existing files. One MATLAB run produces both control laws:

- **unconstrained** (vanilla baseline) → `controllers/sop_bounded_control_acrobot_unconstrained.py`
- **constrained** (bounded-control SOP result) → `controllers/sop_bounded_control_acrobot_result.py`

Each exported module defines `u_opt` (list of 2 exprs), `k1_opt`, `certificate_opt` over the
4-state `[x1,x2,x3,x4] = [q1,q2,dq1,dq2]` (the MATLAB 5-state `x5=x1+x2` aux is substituted back
to 4-state before export). The certificates carry the `-delta/lambda` shift (the delta fix is
already in the shared `reach_avoid_controller.m` / `solvesop_bounded_control.m`).

## Cell layout (same as the dubins synthesis notebook)

1. **md** — title, overview, prerequisites (MATLAB on PATH + SOSTOOLS + Mosek; rab_mpc kernel;
   the solve is slower than dubins — `samples_num=10000`, full manipulator dynamics).
2. **bootstrap** — 4-level repo-root search; put `python/` + `controllers/` on `sys.path`;
   resolve `MATLAB_DIR`.
3. **locate MATLAB** — `find_matlab()` (PATH → fallback `/usr/local/bin/matlab`).
4. **run synthesis** — snapshot the two target files' mtimes; `run_one(matlab, MATLAB_DIR,
   "example_manipulator")`; print `status / design_s / sop_solve_s / wall_s` + stdout tail;
   raise on `status != "ok"`.
5. **confirm exports** — assert both files exist and mtime advanced; print path / size / mtime.
6. **display unconstrained** — importlib-load `sop_bounded_control_acrobot_unconstrained`;
   print `u_opt` (tau1, tau2), `k1_opt`, `certificate_opt`, free symbols, parameter header.
7. **display constrained** — same for `sop_bounded_control_acrobot_result`.
8. **md** — summary; note unconstrained = vanilla baseline, constrained = bounded-control SOP
   result (torque bounds `lb=[-5500;-700]`, `ub=[5500;700]` in the MATLAB synthesis).

## Error handling

- MATLAB not found → explicit message.
- Synthesis failure → surface `run_one` parsed error + stdout tail, raise.
- Missing / stale exports → assert against the expected path / pre-run mtime.

## Testing / verification

Run in place via nbclient (rab_mpc kernel). Verify: 0 cell errors; both controller files
written fresh in `controllers/`; the displayed control laws have the expected free symbols
(`{x1,x2,x3,x4}` for `u_opt`/`certificate_opt`, `{y1,y2}` for `k1_opt`); no free `delta` symbol.
First smoke-test the load/display logic against the archived acrobot controllers.

## Out of scope

- Sub-projects 2 (example-notebook updates) and 3 (comparison figures) — separate specs.
- No MATLAB source changes; controllers loaded/synthesized as-is.
