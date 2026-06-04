# Design: Dubins three-method comparison / paper-figure notebook

**Date:** 2026-06-04
**Topic:** A self-contained notebook that regenerates the Dubins paper figures in one
shot — one per-method figure for each of the three methods (its feasible region + its
single shared-x0 trajectory) and one combined figure (RA-MPC sampling-time sweep +
shared-x0 control-input comparison).

## Goal

One script, run once, produces all Dubins comparison figures for the paper, reusing
the **exact settings** of the three individual notebooks. It can use a single shared
initial point to compare the unconstrained controller's and RA-MPC's control inputs.

## Approach

Self-contained notebook `notebooks/make_paper_figs_dubins.ipynb` (mirrors the archived
`archive/revision/experiments/make_paper_figs_dubins.ipynb`, rebuilt for the current
controllers, MPC parameters, and feasibility criteria). It duplicates the small amount
of solver/rollout logic the individual notebooks use (the archived precedent was also
self-contained); no shared module / refactor of the existing notebooks.

### Reused settings (identical to the 3 individual notebooks)

- Controllers: `controllers/sop_bounded_control_dubins_car_result.py` (RA-MPC / constrained,
  provides `certificate_opt`, `k1_opt`) and `..._unconstrained.py` (uncon, provides `u_opt`,
  `k1_opt`).
- Dynamics `f,g`, output `h=[x1,x2]`, sets `safe_set_x`, `target_set_x` (y-space defs).
- MPC: `dt=0.05`, `N_hor=25`, `u1_max=u2_max=5`, costs `Q_y=5, Qf_y=80, R_u=0.05`,
  `TARGET_POS=[-1.7, 0]`, RK4 discretisation `F_disc`.
- Sampling: `seed=42`, bounds `[-2,-2,2π/3,-1]..[2,2,4π/3,1]`, first 300 candidates
  (safe & outside target).
- Feasibility criteria (final, per the individual notebooks):
  - **Vanilla MPC**: `solve_mpc` with terminal `target(x_N)<=0` is solvable.
  - **RA-MPC**: receding-horizon closed loop with terminal `cert(x_N)>=0` reaches the
    target SAFELY within `T_max=10 s` (F_disc step + per-step warm-started solve).
  - **Unconstrained**: closed-loop ODE under `u_opt` reaches the target SAFELY within 10 s.

### Parameters (top-of-notebook)

- `X0_SHARED = None` — shared initial point; `None` → auto-pick.
- `OUTDIR = <repo>/figures` — figures saved here as `.pdf` + `.png` (dpi=200), also inline.
- `RECOMPUTE_FEAS = False` — reuse cached classification when params match.
- `N_SWEEP_CAND = 150` — candidate subset used for the dt sweep (smaller for speed).
- `DT_SWEEP = [0.025, 0.05, 0.125, 0.25]` with `N = round(1.25/dt) = [50,25,10,5]`
  (fixed horizon time `dt*N = 1.25 s`, matching the archived sweep config).

## Components / cell layout

1. **md** — title, overview, parameter block description.
2. **bootstrap** — `sys.path` (controllers/, python/), create `OUTDIR`.
3. **setup** — load both controllers; symbols/dynamics/sets; numpy + CasADi lambdify;
   `sympy2casadi`; `safe_ca/cert_ca/target_ca`; `F_disc`; `solve_mpc(x0, terminal, dt, N, U_warm)`
   (terminal ∈ {"target","cert"}); `rollout_mpc(x0, terminal, dt, N)` (receding ZOH, returns
   `t,u,x`); `rollout_uncon(x0)` (continuous ODE, returns `t,u,x`); the three feasibility
   predicates.
4. **feasibility** — sample 300 candidates; classify under all 3 methods (joblib parallel);
   cache to `figures/feasibility_cache.npz` keyed by a params signature (reuse unless
   `RECOMPUTE_FEAS` or signature mismatch). Print the three success rates.
5. **shared x0** — use `X0_SHARED` if set; else auto-pick a candidate feasible under ALL
   three methods (iterate farthest→nearest from target, accept first all-3-success). Print
   the chosen x0 (+ candidate id).
6. **rollouts** — from the shared x0: `rollout_mpc(x0,"target")`, `rollout_mpc(x0,"cert")`,
   `rollout_uncon(x0)`; keep `(t,u,x)` for each.
7. **Per-method figures** — ONE figure per method (`fig_dubins_vanilla`, `fig_dubins_rampc`,
   `fig_dubins_uncon`, each `.{pdf,png}`), aligned with the archived c5 layout. Each shows,
   in y-space: the safe/target sets; that method's feasible (blue ○) / infeasible (red ×)
   initial states; and that method's single closed-loop trajectory from the shared x0
   (start/end markers). For the two reach-avoid methods (RA-MPC, uncon) also overlay the
   `k1(y)` vector field as a faint background (vanilla MPC has no `k1`, so no field).
   Per-figure title/legend shows the method's success %.
8. **dt sweep** — for each `(dt,N)` in `DT_SWEEP`, recompute the RA-MPC closed-loop success
   rate (receding, cert terminal, `T_max=10 s`, reach target safely) over the **first
   `N_SWEEP_CAND`** of the same 300-candidate set (joblib parallel); cache to
   `figures/dt_sweep_cache.npz`. (At `dt=0.05` this ≈ the main RA-MPC rate, modulo the
   150-vs-300 subset.)
9. **Fig compare** (`fig_dubins_compare.{pdf,png}`): 3 stacked subplots (archived c6 layout) —
   - top: RA-MPC success rate ρ[%] vs sampling time T=dt (decreasing), vertical line at base dt=0.05;
   - middle: ω(=u1) — uncon continuous line vs RA-MPC ZOH step, ±u1_max band;
   - bottom: a(=u2) — same, ±u2_max band. (control comparison: uncon + RA-MPC only.)
10. **md** — summary listing the saved files + the three success rates.

## Output

`figures/` (created if absent): `fig_dubins_vanilla`, `fig_dubins_rampc`, `fig_dubins_uncon`,
`fig_dubins_compare`, each as `.pdf` (paper) + `.png` (preview), dpi=200. Caches:
`feasibility_cache.npz`, `dt_sweep_cache.npz`.

## Error handling

- MATLAB not needed (pure Python/CasADi). MATLAB-free.
- Missing controller file → clear assert with the expected path.
- Auto-pick finds no all-3-success candidate → assert with guidance (lower expectations /
  set X0_SHARED manually).
- Cache signature mismatch → recompute (never silently use stale numbers).

## Testing / verification

Run in place via nbclient (rab_mpc kernel). Verify: 0 cell errors; the three success
rates match the individual notebooks (Vanilla 81.7%, RA-MPC 95.7%, Uncon 31.3% at base
dt); the dt-sweep ρ is monotone-ish decreasing in dt; shared x0 is feasible under all
three; all eight figure files (4×{pdf,png}) exist and are non-empty.

## Out of scope

- No refactor of the three individual notebooks into a shared module.
- No new controller synthesis (MATLAB); controllers are loaded as-is.
- The archived top-subplot "sweep loaded from CSV" is replaced by an in-script sweep
  (the archived sweep CSVs are from a different/old setup and are not reused).
