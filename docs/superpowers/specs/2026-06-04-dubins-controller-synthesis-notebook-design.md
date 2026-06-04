# Design: Dubins controller-synthesis notebook

**Date:** 2026-06-04
**Topic:** A `notebooks/` notebook that drives the MATLAB synthesis to generate the Dubins-car *unconstrained* and *constrained* control laws.

## Goal

Provide a single notebook in `notebooks/` that, when run, synthesizes both Dubins
control laws by calling the existing MATLAB pipeline via `matlab -batch`, then
loads and displays the two resulting control laws. The previous
`revision/controller_synthesis/` approach has been archived; this is a clean
orchestration notebook over the canonical `matlab/example_dubins_car.m`.

## Approach (chosen)

- Drive the **existing, unmodified** `matlab/example_dubins_car.m` via
  `matlab -batch`, reusing `python/matlab_runner.py` (`find_matlab`, `repo_root`,
  `run_one`).
- **No new `.m` files. No edits to existing files.** Parameters (`lb/ub = ±5`,
  `ds = dv = 4`, `mu_val = 0.1`, `samples_num = 1000`) stay hard-coded in the
  `.m` as the single source of truth.
- One MATLAB run produces both control laws at once:
  - **unconstrained** (vanilla baseline) →
    `controllers/sop_bounded_control_dubins_car_unconstrained.py`
  - **constrained** (bounded-control SOP result) →
    `controllers/sop_bounded_control_dubins_car_result.py`
- Each exported module defines `u_opt` (list of 2 exprs), `k1_opt`,
  `certificate_opt`, plus a `# Parameters` header comment.

Rejected alternatives:
- *Parametrized notebook + thin solve script* — reintroduces the complexity the
  user just archived; not needed.
- *`matlab.engine` in-process* — engine not installed in `rab_mpc`; repo
  convention is `matlab -batch` subprocess.

## Components & data flow

```
notebooks/synthesize_dubins_controllers.ipynb
   │  (uses)
   ├── python/matlab_runner.py  → find_matlab / repo_root / run_one
   │        │  matlab -batch "addpath(matlab/); example_dubins_car; ..."
   │        ▼
   ├── matlab/example_dubins_car.m   (unchanged)
   │        │  reach_avoid_controller + solvesop_bounded_control + export_to_python
   │        ▼
   └── controllers/
         ├── sop_bounded_control_dubins_car_unconstrained.py   (unconstrained)
         └── sop_bounded_control_dubins_car_result.py          (constrained)
              │  (loaded back via importlib)
              ▼
        notebook prints u_opt / k1_opt / certificate_opt for both
```

## Notebook cell layout

1. **Markdown** — title + overview + prerequisites (MATLAB on PATH with
   SOSTOOLS + Mosek; run with the `rab_mpc` kernel).
2. **Repo-root bootstrap** — standard 4-level upward search for a dir containing
   `controllers/`; insert `python/` and `controllers/` onto `sys.path`; resolve
   `MATLAB_DIR = ROOT/matlab`. Matches the existing notebooks' bootstrap cell.
3. **Locate MATLAB** — `find_matlab()`; if `None`, raise with a clear
   install/PATH message.
4. **Run synthesis** — record both target files' mtime (if present);
   `run_one(matlab, MATLAB_DIR, "example_dubins_car")` (default 3600 s timeout); print
   `status / design_s / sop_solve_s / wall_s` and the stdout tail; raise on
   `status != "ok"`, surfacing `result["error"]`.
5. **Confirm exports** — assert both files exist and their mtime advanced past
   the pre-run timestamp; print path / size / mtime for each.
6. **Display unconstrained law** — `importlib` load
   `sop_bounded_control_dubins_car_unconstrained`; sympy-ify and print `u_opt[0]`,
   `u_opt[1]`, `k1_opt`, `certificate_opt`, free symbols, and the Parameters header.
7. **Display constrained law** — same for
   `sop_bounded_control_dubins_car_result`.
8. **Markdown summary** — one paragraph: unconstrained = vanilla baseline (no
   input-bound enforcement); constrained = bounded-control SOP result over the
   sampled state set.

## Error handling

- MATLAB not found → explicit message (how to put `matlab` on PATH).
- Synthesis failure → print `run_one` parsed `error` + stdout tail, raise.
- Missing exports → assert with the expected absolute path.
- Stale exports → mtime comparison against the pre-run snapshot.

## Testing / verification

`rab_mpc` lacks nbconvert, so the notebook is run interactively (its kernel).
Before declaring done, smoke-test with the `rab_mpc` interpreter
(`/home/jianqiang/miniconda3/envs/rab_mpc/bin/python`):

- the bootstrap + `importlib` load/display logic (cells 2, 6, 7) against the
  already-committed controller files, and
- the `find_matlab` / `run_one` wiring (dry, without necessarily running a full
  multi-minute solve).

## Out of scope

- No vector-field plots or closed-loop rollout (scope = generate + display).
- No parameter sweeps, no slack-weight selection, no MATLAB source changes.
