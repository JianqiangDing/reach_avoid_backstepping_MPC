# Dubins demo — paper-figure notebook

Self-contained bundle for the **Dubins car constrained-vs-unconstrained**
demonstration described in [`../../DUBINS_DEMO_PLAN.md`](../../DUBINS_DEMO_PLAN.md).

## Layout

```
dubins_demo/
├── README.md                          # this file
├── dubins_demo.ipynb                   # main notebook (§1–§9 of the plan)
├── matlab/                             # new MATLAB scripts for this demo
│   ├── build_dubins_demo.m             # §4: synthesize con + unc controllers
│   ├── solvesop_bounded_control_hard_xi.m  # §8: hard-bound SOP variant
│   └── test_hardbound_dubins.m         # §8: driver
└── data/                               # generated outputs (cached controllers, hard-bound result)
```

## Existing repo-wide helpers used (not duplicated here)

- `matlab/reach_avoid_controller.m`
- `matlab/solve_vanilla_k1_controller_xi.m`
- `matlab/solvesop_bounded_control_slack_xi.m`
- `matlab/compute_poly_bounds_sampling.m`
- `python/systems.py` (Dubins closed-loop integrator)
- `python/matlab_runner.py` (`run_matlab` helper)

The notebook prepends the repo's `matlab/` directory to MATLAB's path automatically.

## Running

Open `dubins_demo.ipynb` with the `rab_mpc` kernel. §1 sets up paths and parameters; §2–§3
use existing CSVs in the repo `data/` so they run instantly. §4 calls MATLAB to synthesize the
two controllers (cached, so re-runs are fast); §5–§7 visualize and simulate. §8 runs the
hard-bound SOP for the strict-subset demonstration; §9 ties everything together.

## Parameters (defaults)

| param | value | from |
|---|---|---|
| `μ` | 0.1 | paper baseline (`REVISION_NOTES.md` §4) |
| `xi0` | 10 | leverage sweet spot (`REVISION_NOTES.md` §6.4) |
| `ub` | 20 | between con sample-claim (15.8) and unc dense max (~55) |

All overridable in §1 of the notebook.
