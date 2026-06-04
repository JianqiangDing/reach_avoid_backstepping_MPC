# Experiment notebooks

Standalone experiment notebooks for the L-CSS revision live here, kept separate
from the per-controller demonstration notebooks in `../notebooks/`. They are
**Jupyter notebooks** so each step's output can be inspected directly.

Each notebook starts with a small bootstrap cell that puts `../python/` on
`sys.path`; reusable logic shared across these experiments lives in `../python/`
(e.g. `matlab_runner.py`) rather than being copied into each notebook. Notebooks
import the exported controllers from `../controllers/`, reuse the same
initial-condition sampling and the unified success criterion, and write outputs
to `../data/`.

**Unified success criterion** (shared by all experiments): an initial state
`x0` in the safe set but outside the target set is *successful* if its
closed-loop trajectory can be solved all the way until it enters the target set
while remaining inside the safe set throughout — which, for the reach-avoid MPC,
is the empirical signature of the recursive-feasibility guarantee (feasible at
every step until the target). The success rate is the fraction of sampled
initial conditions that are successful.

## Notebooks

- `run_matlab_timing.ipynb` — batch-runs the MATLAB synthesis examples and
  records solve times (offline `design` + `sop_solve`, plus subprocess
  wall-clock). Uses `../python/matlab_runner.py`. _(implemented)_

## Planned notebooks (not yet implemented)

- `sampling_period_sweep.ipynb` — sweep the sampling period
  `T_s in {0.01, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8}` s on a fixed sampled
  initial-condition set, report closed-loop success rate vs `T_s` (CSV + curve).
- `disturbance_sweep.ipynb` — additive disturbance on the closed-loop dynamics
  (`xdot = f + g u_c + d`, `||d||_inf <= dbar`), sweep `dbar` and report success
  rate vs disturbance level (for the response letter).
