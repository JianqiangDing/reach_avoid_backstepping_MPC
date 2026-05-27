# Experiment scripts

Standalone, runnable experiment scripts for the L-CSS revision live here, kept
separate from the per-controller demonstration notebooks in `../notebooks/`.

Each script imports the exported controllers from `../controllers/` and the
shared helpers from `../python/`, reuses the same initial-condition sampling and
the unified success criterion, and writes its outputs to `../data/`.

**Unified success criterion** (shared by all experiments): an initial state
`x0` in the safe set but outside the target set is *successful* if its
closed-loop trajectory can be solved all the way until it enters the target set
while remaining inside the safe set throughout. The success rate is the fraction
of sampled initial conditions that are successful.

## Planned scripts (not yet implemented)

- `sampling_period_sweep.py` — sweep the sampling period
  `T_s in {0.01, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8}` s on a fixed sampled
  initial-condition set, report closed-loop success rate vs `T_s` (CSV + curve).
- `timing_benchmark.py` — offline reach-avoid-set synthesis time and online
  per-step MPC solve time, across the three controllers.
- `disturbance_sweep.py` — additive disturbance on the closed-loop dynamics
  (`xdot = f + g u_c + d`, `||d||_inf <= dbar`), sweep `dbar` and report success
  rate vs disturbance level (for the response letter).
