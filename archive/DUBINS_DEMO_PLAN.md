# Dubins Demo Notebook — Detailed Plan

Status: **draft (2026-05-30)** — awaiting parameter confirmation before implementation.
Target deliverable: `scripts/dubins_demo.ipynb` + supporting MATLAB drivers.
Source data: REVISION_NOTES.md §§4–7; `data/codesign_sop_dubins_car.csv`, `data/lambda_codesign_dubins_car.csv`.

## 1. Goal

Produce the paper's Dubins constrained-vs-unconstrained demonstration cleanly, under the
**methodologically correct framing**:

> *Same `μ`, same `xi0`. The only difference between the two controllers is whether the SOP
> applies the input-bound constraints during `k1` synthesis.*

The notebook is **modular by parameter choice** — each "knob" (`xi0`, `ub`) gets its own
justification section before being used downstream. This lets the user re-pick parameters
and immediately see consequences without re-running the whole pipeline.

**Two paper claims that the notebook must demonstrate visually and quantitatively:**

1. **Bound-respect claim** — the constrained controller satisfies `|u| ≤ ub` while the
   unconstrained one does not. Visualized in §3 (range plot), §5 (field heatmaps with
   `ub` overlay), §7 (closed-loop `|u(t)|`).
2. **Subset claim** — `{V_con ≥ 0} ⊊ {V_unc ≥ 0}` (the constrained reach-avoid set is a
   strict subset of the unconstrained one, because the additional input constraint can only
   remove states, not add them). Visualized in §6 (4-color overlay on the slack-SOP
   controllers — first pass) and §8 (the same 4-color overlay on the hard-bound controller
   — the clean demonstration).

## 2. Why this framing (vs the alternatives we explored)

Three framings were considered (see REVISION_NOTES §§6–7):

| framing | how it works | verdict |
|---|---|---|
| (a) different-`μ` for con vs unc | unc at `μ=0.1`, con at `μ=0.0068` (Dubins §4) | rejected — conflates *parameter tuning* with *constraint enforcement*; not the paper's algorithmic comparison |
| (b) same-`μ`, baseline `xi0` | `μ=0.1`, `xi0=1e-8` for both | rejected — SOP has zero leverage on `u` at baseline `λ` (§6.1), so con ≡ unc; nothing to show |
| **(c) same-`μ`, raised `xi0`** | `μ=0.1`, `xi0=10` (=100·μ) for both | **chosen** — at `xi0=10` the SOP has real leverage (k1's effect ≈84% — §6.4), so con and unc are genuinely different controllers from the *same* algorithm with vs without the bound constraint |

Framing (c) is what the paper claims and what this notebook will demonstrate.

## 3. Parameter choices (defaults)

| parameter | value | justification (§ ref) |
|---|---|---|
| `μ` | **0.1** | Paper's committed Dubins choice; `ra/safe ≈ 0.89` (§4) — large unc reach-avoid set |
| `xi0` | **10** | `=100·μ`; SOP leverage 84% (§6.4); con sample-claim drops to 15.8 vs unc dense 55 |
| `ub` | **20** | Between con sample-claim (15.8) and unc dense max (55); clean separation |
| `dv` | 4 | Paper's choice (degree of `k1` polynomial in `y`) |
| `ds` | 4 | Paper's choice (SOS auxiliary degree) |
| sample count | 1000 | Paper's choice for SOP sample set |
| dense eval `N` | 20000 | Box samples for dense statistics |

All defaults are **overridable in §1 of the notebook**, so the user can rerun any section
under a different parameter without editing code.

## 4. Notebook sections (detailed)

### §1 — Bootstrap & parameters

**Purpose.** Make every downstream choice flow from a single central dict.

**Steps.**
1. `sys.path` bootstrap (`controllers/`, `python/`).
2. Imports: `numpy`, `pandas`, `matplotlib`, `sympy`, `os`, `subprocess`, `re`.
3. Define `PARAMS = {mu: 0.1, xi0: 10.0, ub: 20.0, dv: 4, ds: 4, samples_num: 1000, N_dense: 20000}`.
4. Define paths: `ROOT`, `DATA`, `MATLAB_DIR`.
5. `MATLAB = shutil.which("matlab") or "/usr/local/bin/matlab"`.
6. Helper: `run_matlab(script_name, regenerate=False)` — calls `matlab -batch` if the
   corresponding CSV is missing or `regenerate=True`.

**Verification.** Print `PARAMS`, confirm DATA dir exists, list available CSVs.

---

### §2 — Justify `xi0` (leverage selection)

**Purpose.** Show why `xi0=10` is the sweet spot.

**Inputs.** `data/codesign_sop_dubins_car.csv` (4 rows: xi0 ∈ {1e-8, 0.1, 1, 10}).

**Steps.**
1. Load CSV.
2. Plot 2-panel: left = distinguishability `%` vs `xi0` (log x); right = con sample-claim
   `(ub+slack)` and unc dense `max|u₁|` vs `xi0`.
3. Mark `xi0 = PARAMS["xi0"]` with a vertical line on both panels.
4. Print table: xi0 | distinguishability | sample-claim | unc dense | raset.

**Verification.** At `xi0=10`, distinguishability ≥ 80% AND sample-claim < 0.5 · unc dense.

**Decision point.** User can change `PARAMS["xi0"]` and rerun; later sections pick up the new value.

---

### §3 — Justify `ub` (bound selection)

**Purpose.** Show why `ub=20` is between con's achievable and unc's range.

**Inputs.** At chosen `xi0`, the con sample-claim and unc dense range (from §2 data).

**Steps.**
1. Horizontal bar chart: con sample-claim `[-15.8, 15.8]` and unc dense `[-55, 52]`.
2. Overlay candidate `ub` lines: `{15, 20, 30, 40}` at `±ub`.
3. Annotate each candidate: "con respects (sample): yes/no", "unc respects (dense): yes/no".
4. Print: chosen `ub`, expected demo outcome ("con respects at samples, unc violates densely").

**Verification.** Chosen `ub` lies strictly between con sample-claim and unc dense max.

---

### §4 — Synthesize the two controllers at `(μ, xi0, ub)`

**Purpose.** Get `u_con`, `u_unc`, `V_con`, `V_unc` as evaluable Python objects.

**New MATLAB driver:** `scripts/dubins_demo/matlab/build_dubins_demo.m`
- Inputs: `mu`, `xi0`, `ub` (passed via env vars or command-line args).
- Calls `solvesop_bounded_control_slack_xi(...)` once with these params.
- Returns: `ux_con`, `ux_uncon`, `cert_con`, `cert_uncon`, `slack_opt`, `k1_lambda`.
- Exports to `scripts/dubins_demo/data/controllers_<mu>_<xi0>_<ub>.mat`.

**Notebook steps.**
1. Check if cached `.mat`/`.py` for this `(mu, xi0, ub)` exists; if not, call `build_dubins_demo.m`.
2. Load: get `u_con(x)`, `u_unc(x)`, `V_con(x)`, `V_unc(x)` as `sympy` expressions or as
   `lambdify`-ready callables.
3. Convert to fast NumPy-callable form (per dim).

**Verification.**
- Eval `u_con` at 10000 box samples, confirm sample-claim matches CSV value.
- Eval `V_con` and `V_unc` over box, confirm raset numbers match CSV.

---

### §5 — Controller field comparison

**Purpose.** Visualize where each controller violates `ub`.

**Steps.**
1. Build 2-D grids over `(x1, x2)` at fixed `(θ, v)` slices. Default slices:
   `(θ=π, v=0.5)`, `(θ=π/2, v=0.7)`, `(θ=4π/3, v=0.3)` (3 slices).
2. For each slice, compute `|u_con|`, `|u_unc|` over the grid.
3. Heatmap (2-panel per slice): left = `|u_unc|`, right = `|u_con|`.
4. Overlay: contour at `|u| = ub`, red-shade regions where `|u| > ub`.
5. Overlay: boundary of `{V_unc ≥ 0}` (the reach-avoid set).

**Verification.** unc has large red regions; con has thin red strips between samples (the
overfit, shown honestly).

---

### §6 — Reach-avoid set comparison (subset visualization + slack-SOP)

**Purpose.** Visualize and quantify whether `{V_con ≥ 0} ⊆ {V_unc ≥ 0}` — i.e. whether
the constrained certified reach-avoid set is a (strict) subset of the unconstrained one.
**This is one of the two central paper claims** ("constraint shrinks the certified set"),
so it gets explicit visual treatment here.

**Why the subset is expected.** Logically, adding input constraints to the reach-avoid
problem can only shrink the set of states that admit a feasible solution. So we expect
the constrained certified set to lie inside the unconstrained one — every `x0` certified
under the input bound is also certified without it. Visualizing this directly is the
strongest evidence for the paper's second claim.

**Steps.**
1. Same 2-D slices as §5.
2. Compute `V_con` and `V_unc` over the grid of each slice.
3. **4-color overlay (the subset visualization).**
   - **Green**: `V_con ≥ 0 AND V_unc ≥ 0` — both certify (common reach-avoid region).
   - **Blue**: `V_unc ≥ 0 AND V_con < 0` — only unc certifies; **these are the states the
     input bound "removes" from the reach-avoid set**. If con ⊊ unc cleanly, this region
     is non-empty. *This is the visual proof of the "smaller set" claim.*
   - **Red**: `V_con ≥ 0 AND V_unc < 0` — only con certifies; **unexpected if con ⊆ unc**.
     Non-zero red region means the SOS shapes don't strictly nest (different k1 ⇒ different
     cert geometry). Flag and report.
   - **White / gray**: neither certifies (outside both reach-avoid sets).
4. Compute subset metrics over `N_dense` box samples:
   - `subset_rate = |{V_con ≥ 0} ∩ {V_unc ≥ 0}| / |{V_con ≥ 0}|`
     → fraction of con-cert states that are also unc-certified. Should be ≈ 1 if con ⊆ unc.
   - `removed_frac = |{V_unc ≥ 0} ∩ {V_con < 0}| / |{V_unc ≥ 0}|`
     → fraction of unc-cert states that the input bound removes. Should be > 0 if the
     constraint actually shrinks the set.
   - `unexpected_frac = |{V_con ≥ 0} ∩ {V_unc < 0}| / |{V_con ≥ 0}|`
     → fraction of con-cert states NOT certified by unc. Should be ≈ 0; non-zero values
     indicate SOS-shape mismatch (the slack mechanism can produce con-certified states
     outside unc's cert by paying slack).
   - existing `raset_con`, `raset_unc`, `ra/safe_con`, `ra/safe_unc`.
5. Report and interpret:
   - **`subset_rate ≈ 1` + `removed_frac > 0`**: con ⊊ unc cleanly. **Paper's claim shown.** ✓
   - **`subset_rate ≈ 1` + `removed_frac ≈ 0`**: con ≈ unc essentially; the constraint
     didn't shrink the set (slack absorbs the bound).
   - **`unexpected_frac > 0`**: the slack SOP produces a con-cert region partly outside
     unc's region — flag this as a "shape mismatch" and emphasize §8 (hard-bound) as the
     cleaner test.

**Verification.** All three metrics (`subset_rate`, `removed_frac`, `unexpected_frac`) computed
and printed; 4-color visualization clear in at least one slice.

**Note.** With slack SOP we expect `subset_rate` near 1 but `removed_frac` possibly small.
The clean, strict subset demonstration is the goal of §8 (hard-bound SOP).

---

### §7 — Closed-loop simulation

**Purpose.** Show controller behavior along trajectories, not just at static field samples.

**Steps.**
1. Sample 20 initial conditions `x0` uniformly from `{V_unc ≥ 0}` (the reach-avoid set).
2. For each `x0`, simulate both con and unc closed-loop dynamics using
   `python/systems.py` (Dubins integrator, `solve_ivp` with terminal target event).
3. Per trajectory: collect `(t, x(t), u(t))`; record reach (yes/no), max `|u(t)|`,
   fraction of timesteps with `|u(t)| ≤ ub`.
4. Plot:
   - `(x1, x2)` trajectories for both controllers (20 trajectories overlaid, color by controller).
   - `|u(t)|` vs `t` (separate axes for u1, u2), with `ub` line.
5. Summary table:
   - reach rate (con, unc)
   - mean fraction of timesteps respecting `ub`
   - mean max `|u(t)|`

**Verification.** con's bound-respect rate at the sample-evaluated points is higher than unc's
(in line with the PAC guarantee).

---

### §8 — Hard-bound SOP variant (the strict-subset demonstration)

**Purpose.** Show the strict subset relationship `{V_con ≥ 0} ⊊ {V_unc ≥ 0}` cleanly —
the unambiguous form of the paper's "smaller set" claim. This is §6's subset analysis
done with a controller that does NOT absorb constraint violations into slack.

**Hypothesis.** The slack mechanism keeps the certificate feasible everywhere by absorbing
constraint violation into `s_i`, which can produce a `V_con` that is *not* strictly nested in
`V_unc` (the `unexpected_frac` region in §6). If we drop the slack and force `|u_i| ≤ ub` as hard
constraints, infeasibility in violating regions should naturally:
1. Shrink `{V_con ≥ 0}` (smaller set), AND
2. Force `{V_con ≥ 0}` to lie *inside* `{V_unc ≥ 0}` (the strict subset, no `unexpected_frac`).

**New MATLAB file:** `scripts/dubins_demo/matlab/solvesop_bounded_control_hard_xi.m`
- Copy of `matlab/solvesop_bounded_control_slack_xi.m` (repo-wide helper) with:
  - Drop slack variables `s_i` from the SOS program.
  - Replace soft `sosineq(prog, num/den - lb(i) + s_i)` with hard `sosineq(prog, num/den - lb(i))`.
  - Objective: still `delta` only (no `Σ s_i`).
  - Add try/catch around `sossolve` to detect infeasibility.

**New driver:** `scripts/dubins_demo/matlab/test_hardbound_dubins.m`
- Run `solvesop_bounded_control_hard_xi` at `(mu=0.1, xi0=10, ub=20)`.
- Report: solve status (`feasible/infeasible/unknown`), `raset`, dense `max|u|`,
  sample-level violation count.
- Write `scripts/dubins_demo/data/hardbound_<mu>_<xi0>_<ub>.csv`.

**Notebook steps.**
1. Trigger MATLAB run if CSV missing.
2. Load `V_con_hard` and `u_con_hard`.
3. **Repeat the §6 4-color subset visualization with `V_con_hard` in place of `V_con_slack`.**
   Same slices, same color scheme; this is the headline figure for the set-size claim.
4. Compute the same three subset metrics with the hard-bound controller:
   - `subset_rate_hard`, `removed_frac_hard`, `unexpected_frac_hard`.
5. Tabulate side-by-side with §6's slack-SOP metrics:

   | metric | slack SOP (§6) | hard SOP (§8) | expectation |
   |---|---|---|---|
   | `subset_rate` | ? | ? | hard closer to 1 |
   | `removed_frac` | ? | ? | hard larger |
   | `unexpected_frac` | ? | ? | hard ≈ 0 |
   | `raset` | 0.42 | ? | hard < slack |

6. Three branches:

   - **`raset_hard < raset_slack` AND `unexpected_frac_hard ≈ 0`:** paper's strict-subset
     claim holds cleanly. Use the §8 4-color figure as the headline; finalize §9.
   - **`raset_hard ≈ raset_slack`:** the cert shape doesn't depend on the slack mechanism;
     the "smaller set" needs a different argument. Reframe the paper's second claim around
     something we can show (e.g., a quantitative `|u|`-vs-set-size Pareto curve).
   - **Hard SOP infeasible:** the chosen `ub` is too tight for hard enforcement at this
     `(μ, xi0)`. Try `ub ∈ {25, 30}` or drop back to the slack version and reframe.

**Verification.** SOS solver status printed clearly; subset metrics computed and tabulated;
hard-SOP outcome documented; the 4-color figure produced.

---

### §9 — Final comparison table + paper-ready figures

**Purpose.** Single-page summary suitable for the paper.

**Contents.**
1. Two-row comparison table:

   | | unconstrained | constrained (slack) | constrained (hard, §8) |
   |---|---|---|---|
   | `μ`, `xi0` | 0.1, 10 | 0.1, 10 | 0.1, 10 |
   | SOP bound | — | `ub = 20` | `ub = 20` |
   | sample max `|u₁|` | ~55 | ≤20 ✓ | ≤20 ✓ |
   | dense max `|u₁|` | ~55 | 92 | ? |
   | `raset` (frac of box) | 0.42 | 0.42 | ? |
   | `ra/safe` | 0.89 | 0.89 | ? |
   | `subset_rate` (con ⊆ unc) | — | ? | ? (expect ≈ 1) |
   | `removed_frac` (states constraint removes) | — | ? | ? (expect > 0) |
   | `unexpected_frac` | — | ? | ? (expect ≈ 0) |
   | reach rate (closed loop) | (§7 number) | (§7 number) | — |

2. Selected figures (3–4 total): one field comparison (§5), one set comparison (§6 or §8),
   one trajectory plot (§7).

3. Caveats list (honest with reviewers):
   - Bound respect is at sample level (PAC); dense overshoot present.
   - Set-size separation: shown by hard SOP / not shown (depending on §8 result).

---

## 5. Files

All code and data for **this test only** live under a single dedicated folder
`scripts/dubins_demo/`. Existing repo-wide helpers (in `matlab/`, `python/`) are referenced
via `addpath`/`sys.path` — not duplicated.

### Folder layout

```
scripts/dubins_demo/
├── README.md                              # short overview + pointer to ../../DUBINS_DEMO_PLAN.md
├── dubins_demo.ipynb                      # main notebook (§1–§9)
├── matlab/
│   ├── build_dubins_demo.m                # §4: synth con + unc at (μ, xi0, ub), export controllers
│   ├── solvesop_bounded_control_hard_xi.m # §8: hard-bound SOP variant (no slack)
│   └── test_hardbound_dubins.m            # §8: driver that calls the hard-bound solver
└── data/
    ├── controllers_<μ>_<xi0>_<ub>.mat     # cached con + unc controllers + certificates (§4 output)
    └── hardbound_<μ>_<xi0>_<ub>.csv       # §8 result (raset, subset metrics, dense |u|)
```

(Filenames inside this folder drop the `dubins_demo_` prefix — the folder itself provides the
namespace. `<μ>_<xi0>_<ub>` are the active parameter values, so changing them creates a new
cache file without overwriting prior runs.)

### Existing repo files referenced (read only — not duplicated)

| path | used by |
|---|---|
| `data/codesign_sop_dubins_car.csv` | §2 (xi0 leverage table) |
| `data/lambda_codesign_dubins_car.csv` | §2 context (feasibility across xi0) |
| `data/separation_dubins_xi*.csv` | §7 reference (existing simulation data) |
| `matlab/reach_avoid_controller.m` | §4 (symbolic synthesis) |
| `matlab/solve_vanilla_k1_controller_xi.m` | §4 (vanilla k1) and §8 |
| `matlab/solvesop_bounded_control_slack_xi.m` | §4 (slack-SOP controller) |
| `matlab/compute_poly_bounds_sampling.m` | §4 (dense `|u|` bounds) |
| `python/systems.py` | §7 (Dubins closed-loop integrator) |
| `python/matlab_runner.py` | `run_matlab` helper in §1 |

### Path bootstrap (notebook §1)

```python
import os, sys
HERE = os.path.dirname(os.path.abspath("__file__"))   # scripts/dubins_demo/
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
for sub in ("python",):
    p = os.path.join(ROOT, sub)
    if p not in sys.path: sys.path.insert(0, p)
DATA_REPO  = os.path.join(ROOT, "data")        # existing CSVs (read-only here)
DATA_LOCAL = os.path.join(HERE, "data")        # this demo's own outputs
MATLAB_REPO  = os.path.join(ROOT, "matlab")    # existing matlab/
MATLAB_LOCAL = os.path.join(HERE, "matlab")    # this demo's matlab/
os.makedirs(DATA_LOCAL, exist_ok=True)
```

MATLAB runs done from the notebook prepend both `matlab/` directories to the path:
```python
def run_matlab(script_name):
    cmd = [MATLAB_BIN, "-batch",
        f"addpath('{MATLAB_REPO}'); addpath('{MATLAB_LOCAL}'); {script_name}; exit(0);"]
    subprocess.run(cmd, check=True)
```

## 6. Execution order

1. **§1–§3** first (existing CSVs, fast iterate on `xi0`/`ub` choice).
2. **§4** (synthesize at chosen params, **cache** the controllers — re-running notebook should
   not re-synthesize unless params change).
3. **§5–§7** (visualization + closed loop) — these consume §4 cache.
4. **§8** (the new hard-bound experiment) — independent; can run in parallel with §5–§7.
5. **§9** (summary) — consumes everything above.

## 7. Success criteria

| section | success means |
|---|---|
| §2 | distinguishability ≥ 80% at chosen `xi0` |
| §3 | con sample-claim `< ub < unc dense max` |
| §4 | sample-claim from §4 controllers matches §2 CSV value within 1% |
| §5 | unc has visible regions with `|u| > ub`; con's same regions are smaller (or zero at samples) |
| §6 | `raset` matches CSV; 4-color subset visualization rendered in ≥1 slice; `subset_rate`, `removed_frac`, `unexpected_frac` all computed and printed |
| §7 | both reach rates reported; con's `ub`-respect rate > unc's |
| §8 | hard SOP outcome documented; 4-color subset figure produced; subset metrics tabulated side-by-side with slack-SOP (§6); strict-subset verdict reached (holds / doesn't / infeasible) |
| §9 | one paper-ready table, 3–4 selected figures |

## 8. Open decisions (before implementation)

1. **`ub` default**: confirm `20`. Alternatives in `(16, 55)` also work; smaller `ub` makes con's
   demo more impressive but stresses the hard SOP feasibility in §8.
2. **Field slices in §5/§6**: confirm default 3 slices `(θ, v) ∈ {(π, 0.5), (π/2, 0.7), (4π/3, 0.3)}`.
3. **Closed-loop sim in §7**: 20 ICs uniformly from the RA set, or a fixed grid? Default: random
   sampling with `rng(seed)` for reproducibility.
4. ~~**§8 priority**: high or optional?~~ — **resolved: high.** The subset claim is one of the
   two paper claims (§1), so §8 is the centerpiece of the set-side demonstration and must be
   implemented before §9. The §6 subset analysis on the slack-SOP controllers is the first pass;
   §8 provides the clean, strict-subset version on the hard-bound controllers.
5. **Cache format for §4 controllers**: `.mat` (faster, MATLAB-native) or `.py` (sympy-loadable,
   editable)? `.py` matches the existing `controllers/` convention; `.mat` is faster.

## 9. Open risks

- **§4 caching pitfall**: if the user changes a parameter not in the filename (e.g., `dv`),
  the cache won't invalidate. Mitigation: include all SOP-affecting params in the cache filename
  hash, or print a warning if `PARAMS` differs from cached file's metadata.
- **§8 hard-SOP infeasibility**: realistic at tight `ub`. Mitigation: §8 must handle infeasibility
  gracefully (don't crash the notebook) and report the closest feasible `ub`.
- **Dubins closed-loop fragility (§7)**: existing `data/separation_dubins_xi*.csv` showed reach
  rates of only 3–4/50 (§7). The notebook should report these honestly and not over-claim closed-
  loop success.

---

*End of plan. Awaiting confirmation on §8 (open decisions) before implementation begins.*
