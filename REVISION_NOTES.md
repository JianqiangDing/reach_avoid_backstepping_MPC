# LCSS Revision — Experiment Progress Notes

Branch: `lcss-revision`. Goal: make the reach-avoid backstepping controllers respect
the input bounds (Dubins `|u| <= 5`, manipulator `|tau| <= 500`) and characterize when
this is feasible. Verification is by **sampling** throughout (chosen for simplicity and
consistency with the scenario-optimization nature of the synthesis).

## 1. Feedback-linearizability / decoupling-matrix region

The feedback law `u = A(x)^{-1} b(x)` blows up where the decoupling matrix `A(x)` is
(near-)singular. Singular loci (from `det A`):

| system            | det A           | singular locus              |
|-------------------|-----------------|-----------------------------|
| Dubins            | `-v`            | `v = 0`                     |
| manipulator       | `∝ sin(x2)`     | elbow angle `x2 ∈ {0, pi}`  |
| double integrator | `1`             | none                        |

Tooling: `matlab/check_decoupling_invertibility.m` + `matlab/run_decoupling_check.m`
write `data/decoupling_<sys>.csv`; `scripts/check_decoupling_matrix.ipynb` visualizes.

**Checker fix:** the manipulator must be sampled on the physical manifold `x5 = x1 + x2`.
Sampling `x5` independently fabricates off-manifold configurations (`x5 - x1 ≈ 0, pi`)
that the real system never visits and that spuriously flag as singular. Added via the
optional `constraint_fn` argument to the checker.

**Restricted sampling regions** (in the example scripts) to stay nonsingular:
- Dubins:      `v ∈ [0.1, 1.0]`     (was `[-1, 1]`)
- manipulator: `x2 ∈ [0.8, pi-0.8]` (was `[0.8, pi]`); `x5 = x1 + x2` on-manifold

After restriction all three systems are well-conditioned (min singular value ~0.1–1,
condition number ~1–12).

## 2. Control-bound re-synthesis on the restricted regions

| system       | before (old region) | after (FL region) | limit |
|--------------|---------------------|-------------------|-------|
| Dubins `u1`  | `[-30477, 18337]`   | `[-71, 92]`       | ±5    |
| manip `tau1` | `[-5258, 6164]`     | `[-5757, 6042]`   | ±500  |

Excluding the singularity helped Dubins ~330x (removed the `1/v` amplification) but
neither example is within bound. The residual blowup is **not** from `A^{-1}` (the
controller denominators equal `det A`, which is bounded over the restricted regions) —
it is from the numerator `b`.

## 3. k1 regularization does NOT help

Tested an L2 (Tikhonov) penalty on the `k1` coefficients in `solve_k1_controller_sop.m`.
With a strong weight the solver drove `||k1||^2 -> ~2e-8` (`k1 ≈ 0`), yet the control
bounds were unchanged. The control magnitude is **k1-independent**: `k1` only sets the
top-level virtual output velocity, while the bulk of `u` cancels the drift and stabilizes
the relative-degree chain per the certificate (high-degree safe-set gradients, `mu`/`lambda`).
The regularization was reverted. More precisely, `u` is *affine* in `k1` but the leverage is
far too weak; §5–§6 confirm this (slack incentive, and the full `λ`/co-design study).

## 4. mu_val is the effective lever (with a reach-avoid set tradeoff)

`matlab/sweep_mu_dubins.m` sweeps `mu` and reports both the control bound and the
reach-avoid set size (fraction of the safe set certified by the solved certificate).
Dubins:

| mu      | u1            | within ±5 | reach-avoid set (frac of safe) |
|---------|---------------|-----------|--------------------------------|
| 0.006   | `[-4.40, 3.99]` | yes     | 0.20                           |
| ~0.0068 | ~`±5`         | boundary  | ~0.21                          |
| 0.01    | `[-7.33, 6.65]` | no      | 0.28                           |
| 0.05    | `[-36.7, 39.3]` | no      | 0.76                           |
| 0.1     | `[-73.3, 78.5]` | no      | 0.89                           |

`u` scales ~linearly with `mu` (`u1 ≈ ±733·mu`). Strict `|u| <= 5` is achievable at
`mu ≈ 0.0068`, but the certified reach-avoid set then covers only ~21% of the safe set
(vs ~89% at `mu = 0.1`). The tradeoff is fundamental and partly **physical**: with
`|u| <= 5` the car cannot steer every safe state to the target.

The committed paper controller (`mu = 0.1`) does **not** strictly satisfy `|u| <= 5` over
its region either — its reported "validity" is the trajectory/sampling-based claim
(the tested initial conditions succeed and stay within bounds).

## 5. Soft-constraint (slack) SOP — definitive test of whether the SOP can bound `u`

Built a slack variant (**new files; originals untouched**): `matlab/solve_k1_controller_sop_slack.m`,
`matlab/solvesop_bounded_control_slack.m`, `matlab/test_slack_dubins.m`, `matlab/test_slack_manipulator.m`.

Two changes vs the original SOP:
1. **Sample the whole reach-avoid region** (safe, outside target, `cert ≥ 0`); drop the original's
   "already-feasible" pre-filter. *Root cause this exposed:* `solvesop_bounded_control.m:92‑94`
   keeps only samples where the **unconstrained** controller already satisfies the bound
   (`valid_indices`), so the original SOP **never constrains the violating region** — the main
   reason constrained ≈ unconstrained.
2. **Relax** `|u_i| ≤ ub_i` to `|u_i| ≤ ub_i + s_i`, `s_i ≥ 0`, and **minimize `delta + Σ s_i`**.
   Always feasible; `s_i` reports the achievable-bound margin (`achievable = ub + s`).

Constraints actually added: Dubins 398 samples → 1592 control scalar ineqs + 398 (3×3) certificate
LMIs; manipulator 446 samples → 1784 + 446.

Results (requested ub: Dubins ±5, manip ±500):

| system          | slack `s*`     | achievable `ub+s` | slack-controller **dense** range        | vs unconstrained |
|-----------------|----------------|-------------------|------------------------------------------|------------------|
| Dubins (mu=0.1) | `[75.8, 5.0]`  | `[±80.8, ±10]`    | `u1∈[-70.9, 92.1]`, `u2∈[-11.3, 8.6]`    | **identical**    |
| manip (mu=15)   | `[8420, 1809]` | `[±8920, ±2309]`  | `tau1∈[-16200, 10480]`                   | **worse**        |

- **Dubins — no improvement.** Even sampling the whole region *and* explicitly minimizing the
  control slack, the controller is pointwise identical to the unconstrained one → the SOP-on-`k1`
  cannot reduce `|u|` (strongest-possible confirmation of §3).
- **Manipulator — worse.** The degree-2 `k1` over 446 sparse 5-D samples **overfits**: slack at the
  samples is ±8920 but the dense range is ±16200 (blows up between samples, amplified by `A^{-1}`).
  The whole-region scenario approach can *degrade* the controller when samples are sparse.
- **Slack underestimates the true range** (Dubins 80.8 vs dense 92; manip 8920 vs dense 16200) — the
  sample-generalization gap. Use the **dense** `compute_poly_bounds_sampling` range as the achievable
  bound, not `s*`.

Net: the slack reformulation is well-posed and a good diagnostic but **confirms, not fixes** the
limitation. Four independent checks now agree (control-field range, pointwise field diff, closed-loop
simulation, and this incentivized slack solve): the input-constraint SOP does **not** yield a distinct,
more-bound-respecting controller. **`mu` is the only lever** (§4); do not adopt the slack controllers.

## 6. Can `λ` (hence `k1`) bound `u`? — full `λ`/co-design study

Follow-up to §3/§5. **New files (originals untouched):**
`matlab/solve_vanilla_k1_controller_xi.m` (the SOS floor `xi0` exposed as an argument),
`matlab/solvesop_bounded_control_slack_xi.m`, and tests `matlab/test_lambda_effect.m`,
`test_lambda_sweep.m`, `test_lambda_codesign.m`, `test_codesign_sop.m`, `test_manip_collapse.m`.
Visualized in `scripts/compare_lambda_effect.ipynb` (6 sections; data in `data/lambda_*.csv`,
`codesign_sop_*.csv`, `manip_collapse.csv`).

**6.1 `λ` enters `u` only via `0.5·λ·(Lf h − k1)`, and `λ` is pinned at the `1e-8` floor.**
For relative degree 2, `b_i = μ·Dψ_i + Σ_j J_k1(i,j)·Lf h_j + 0.5·λ·(Lf h_i − k1_i) − Lf² h_i`,
`u = A⁻¹ b`. `solve_vanilla_k1_controller` only imposes `λ ≥ xi0 = 1e-8` and minimizes `δ`, so the
solved `λ` sits at the floor (Dubins `1.1e-8`, manip `1.7e-7`). `test_lambda_effect`: `u(λ=computed)`
vs `u(λ=0)` differ by ≤1e-8 % of `|u|` (bit-identical); even `k1=0` (value **and** Jacobian) changes
`u` by ≤4e-6 % (Dubins) / ≤0.002 % (manip). At the solved point `k1`/`λ` are immaterial.

**6.2 Raising `λ` to `μ`-scale does NOT activate the (fixed) `k1`.** `test_lambda_sweep` holds the
obtained `k1` and sweeps `λ`. `μ` multiplies the *large* safe-set gradient `Dψ`; `λ` multiplies the
*small* `(Lf h − k1)`, so equal gains ≠ equal magnitudes. At `λ=μ`, k1's effect on `u` is still
~1e-6 % (Dubins) / ~0.04 % (manip); it reaches 1 % only at `λ≈30·μ` (manip), never within `λ≤1e3` (Dubins).

**6.3 Co-design (raise `xi0`, re-solve `k1`) activates `k1` but never *reduces* `|u|`.**
`test_lambda_codesign` raises the `xi0` floor and re-solves `k1` (λ stays a decision variable, not
fixed). Feasible up to `100·μ`; the re-solved `k1` grows with `λ` so k1's effect climbs (manip 8 % at
`λ=μ`, ~100 % at `λ=100·μ`). **But `max|u|` never shrinks** — flat (set by `μ·Dψ`), then blows up
(Dubins 356→8109 at `λ=100`; manip →2.7e13). Making `k1` matter only *adds* control effort.

**6.4 Co-design + slack SOP — samples still can't bound `u`; the manipulator set collapses.**
`test_codesign_sop` runs the whole-region slack SOP at the raised `λ`. **Dubins**: raising `xi0`
makes constrained ≠ unconstrained (84 % distinct at `λ=10`) and drops the *sample-claimed* bound 5×
(`ub+slack` 80→16), **but the dense `|u1|` stays ±92** — the SOP **overfits** the 242 samples (the §5
gap, now amplified by the leverage). **Manipulator**: raising `λ` even to `0.1` (≪ μ=15) **collapses
the reach-avoid set to empty** (`raset` 0.051→0; the sampler returns "0 reach-avoid samples").

**Mechanism of the manipulator collapse (measured — `test_manip_collapse.m` → `data/manip_collapse.csv`).**
For relative degree 2 the certificate is `V = safe(h) − Σ_i (1/(2μ_i))·(Lf h_i − k1_i)²`, and `λ`
enters **only** through the solved vanilla `k1` (`ks(1)=k1`). Crucially `k1` is a polynomial in the
**output** `y=h(x)` (end-effector position) only, while `Lf h_i = −l1·sin(x1)·x3 − l2·sin(x5)·(x3+x4)`
depends on the joint **velocities** `x3,x4` — which `k1` structurally cannot represent. Sweeping the
floor `xi0` (over the FL box, restricted to `{safe≥0}`, ~1296 samples):

| xi0 (≈λ) | max\|k1\| | median penalty | max penalty | min V | reach-avoid frac |
|----------|-----------|----------------|-------------|-------|------------------|
| 1e-8 (baseline) | 34   | 1.2   | 39    | −38    | 0.046 |
| 1e-4            | 5512 | 2.8e4 | 1.0e6 | −1.0e6 | 0     |
| 1e-2            | 4522 | 1.8e4 | 6.9e5 | −6.9e5 | 0     |
| 0.1             | 7841 | 5.5e4 | 2.1e6 | −2.1e6 | 0     |

Reading the table — the collapse is sharp, not gradual:
1. **`k1` jumps ~160× the instant `λ` leaves the `1e-8` floor** (34 → ~5500). At `λ≈0` the vanilla
   solve only minimizes `δ` and leaves `k1` tiny; once `λ` is a *genuinely* positive number the
   reaching condition `∇S·k1 ≳ λ·S` forces `k1` large enough to dominate `λ·S` over the safe set.
   The exact value of `λ∈[1e-4, 0.5]` barely matters — `k1` sits in the "large" regime (~5–8k) throughout.
2. **The penalty scales like `k1²`** (when `|k1| ≫ |Lf h|`, `(Lf h − k1)² ≈ k1²`): max penalty goes
   `39 → 1.0e6`, i.e. `≈ (5512/34)² ≈ 2.6e4×`, matching the squared `k1` jump.
3. **The penalty (~1e6) dwarfs the safe-set value** (`safe_m = −(4(y1−2)−2y2³)² + 0.8·y2³ + 10`, so
   `safe ≲ 10`), so `V = safe − penalty ≈ −1e6` **everywhere** → `{V≥0}` is empty.

Why the manipulator and not Dubins: the manipulator's set is **already razor-thin at baseline**
(`min V = −38`, only **4.6 %** of the box has `V≥0`), so it has essentially no margin to absorb *any*
`k1` growth — even `λ=1e-4` (≪ μ=15) empties it. Dubins keeps a fat margin (41 %) and its `k1` stays
moderate, so it tolerates `λ` up to ~10. (μ does not fight this — μ is *fixed* here; only `k1`,
through `λ`, changes.)

**Net.** Five independent checks agree (control-field range §2, pointwise diff, closed-loop sim, slack
SOP §5, and this `λ`/co-design study §6): the input-constraint SOP cannot produce a distinct, more-
bound-respecting controller, and `λ`/certificate-shaping cannot bound `u` while keeping a valid
reach-avoid set. **`μ` is the only lever** (§4).

## 7. Feedback-level constrained-vs-unconstrained SEPARATION test (2026-05-30)

Question: at `mu = 0.1`, can we exhibit an `x0` in `X_RA` where the CONSTRAINED feedback
controller reaches the target respecting a (freely chosen) bound while the UNCONSTRAINED
one violates it — i.e. a feedback-level demonstration of C1? Raise `xi0` to 50x-100x `mu`
so `k1` has leverage (§6), synthesize both controllers at the same `lambda`, simulate the
closed loop `dx/dt = f + g·u`, and search initial conditions.

New files (originals untouched, no export to `controllers/`):
`matlab/test_separation_dubins.m`, `matlab/test_xi_sweep_dubins.m`;
data in `data/separation_dubins_xi*.csv`, `data/xi_sweep_dubins.csv`, `data/sep*_xi*.csv`.

**7.1 Direct `X_S\X_T` sampling is INFEASIBLE.** Sampling the constrained slack SOP directly
over `X_S\X_T` (skipping the vanilla pre-solve) makes Mosek infeasible (`feasratio ≈ -0.9997`,
residual `3.5e4`): the certificate LMI `M ⪰ 0` (i.e. `V = psi − Σ(1/2μ)(η−k1)² ≥ 0`) is a
HARD per-sample constraint and cannot hold over the whole `X_S\X_T` (`X_RA` is a strict
subset). The two controllers then come out as near-identical garbage. **=> the sampling region
must be a certificate-feasible subset; the "first vanilla solve" that defines that region is
necessary, not redundant.**

**7.2 `xi0` sweep `[5,6,7,8,9,10]` (=50x..100x `mu`): NO clean separation at any `lambda`.**
(vanilla-cert-region slack SOP; 50 random `x0` in `X_RA` per `xi0`; demo bound `ub = 5`.)

| xi0=λ | con achiev \|u1\| (=ub+slack) | dist% | **con dense \|u1\|** | unc dense \|u1\| | reach con | reach unc | sep1 | sep2 |
|-------|------------------------------|-------|----------------------|------------------|-----------|-----------|------|------|
| 5  | 39.6 | 23% | **76.2** | 76.6 | 3/50 | 2/50 | 0 | 0 |
| 6  | 35.1 | 33% | **71.3** | 78.9 | 4/50 | 1/50 | 0 | 0 |
| 7  | 28.9 | 57% | **70.4** | 57.6 | 3/50 | 1/50 | 0 | 0 |
| 8  | 25.6 | 68% | **83.0** | 66.0 | 3/50 | 1/50 | 0 | 0 |
| 9  | 16.1 | 90% | **85.3** | 58.9 | 3/50 | 1/50 | 0 | 0 |
| 10 | 15.8 | 86% | **84.2** | 67.8 | 3/50 | 0/50 | 0 | 0 |

- **con dense `|u1|` ≈ 70-85 at EVERY `lambda`** — never shrinks. The slack (sample-claimed
  bound) drops `39.6 → 15.8`, but the dense reality stays `~80` (the §6.4 overfit). **=> the
  constrained controller respects NO set-wide bound `< 80`; "con respects `ū` over `X_RA`" is
  impossible for any `ū < 80`.**
- reach&safe rates are very low (con 3-4/50, unc 0-2/50): the raw feedback `u = b/(−v)` is
  singular at `v=0` and fragile in continuous closed loop. Trajectories that DO reach are
  benign (small `|u|` for BOTH controllers, so no separation along them); trajectories with a
  large-`|u|` gap almost never reach. `sep1` (two-trajectory) = `sep2` (uncon-along-con-traj
  counterfactual) = **0** for all six `lambda`.

**7.3 The single marginal hit is not credible.** An earlier 120-`x0` run at `xi0=10` found
exactly ONE separation: `x0=[-1.572, 0.581, 3.637, 0.590]`, con peak `|u|=[4.62, 4.92] ≤ 5`,
unc peak `|u1|=5.29 > 5` — ~1/120, margin ~0.3, cherry-picked; it does not survive a denser
search.

**Net (7th independent check).** The feedback-level "constrained respects / unconstrained
violates" separation is **not achievable at `mu = 0.1` across the whole 50x-100x `mu` `lambda`
range**, consistent with §3-§6: the `k1`-SOP cannot produce a bound-respecting feedback
controller. **=> the input-constraint demonstration must live at the MPC level** (the MPC
enforces `|u| ≤ ub` as hard constraints; the input-constrained reach-avoid set is the
**terminal set** enabling recursive feasibility), and **C1 should be framed as "a terminal set
consistent with the input constraints", NOT "a bound-respecting feedback law".**

## 8. Manipulator μ-sweep (parallel to §4, 2026-05-30)

**Why this test.** §7 pivoted the input-constraint demo to the MPC level. Before committing to
that, check whether a different `mu` can give the manipulator a "sufficiently fat" reach-avoid
set, then revisit `k1`/`λ` adjustment on top of that fatter baseline. New script
`matlab/sweep_mu_manipulator.m` solves vanilla `k1` once (μ-independent), binds it into `u(x;mu)`
and `V(x;mu)`, and scans μ — reporting **`ra/safe`** = fraction of the *safe set* certified
(the meaningful denominator, matching §4) and the **dense** `|τ|` range over the certified RA set.

| μ | ra/safe | τ1 range (RA set) | max\|τ\| | within ±500 |
|---|---|---|---|---|
| 5             | 0.54 | `[-2370, 1660]`   | `2.4e3` | no |
| 10            | 0.66 | `[-5460, 4990]`   | `5.5e3` | no |
| **15 (paper)** | **0.72** | **`[-10160, 8120]`** | **`1.0e4`** | no |
| 30            | 0.84 | `[-47000, 33700]` | `4.7e4` | no |
| 50            | 0.89 | `[-3.5e5, 7.5e4]` | `3.5e5` | no |
| 100           | 0.94 | `[-6.9e5, 3.3e5]` | `6.9e5` | no |
| 200           | 0.96 | `[-1.4e6, 1.2e7]` | `1.2e7` | no |
| 500           | 0.99 | `[-3.4e6, 3.1e7]` | `3.1e7` | no |
| 1000          | 1.00 | `[-6.8e6, 6.2e7]` | `6.2e7` | no |

**Three findings:**

1. **The §6.4 "razor-thin margin (5%)" framing was on the wrong denominator.** That 5% is
   *fraction of the box*; the **safe set is only ~6% of the box** to begin with (1908/30000 samples
   satisfy `safe≥0`). On the *safe-set* fraction — the meaningful one, matching §4 — the manipulator's
   baseline RA set is **72%**, comparable to Dubins' 89% at μ=0.1. What makes the §6.4 collapse
   catastrophic is `min V ≈ -38` (close to zero relative to penalty magnitudes ~1e6 once `k1`
   jumps), **not** a small set-fraction. The "redundancy" is in `min V`, not in set-fraction.

2. **μ↔ra/safe monotonic for the manipulator** (analog of §4 Dubins confirmed): 54%→100% as μ
   goes 5→1000. *Larger μ → larger reach-avoid set* holds structurally for both systems.

3. **μ↔\|τ\| tradeoff is far steeper than Dubins.** Where Dubins paid `u1 ≈ ±733·μ` (linear),
   manipulator max\|τ\| grows super-linearly: `2.4e3 → 1.0e4 → 3.5e5 → 1.2e7` for `μ = 5 → 15 → 50 → 200`.
   **No tested μ keeps `|τ| ≤ 500`** — even μ=5 is 5× over while only 54% ra/safe; the
   strict-bound regime (analog of Dubins' `μ≈0.0068`) would be at μ<5, giving an unhelpfully small set.

**Implication for the "fatten-then-adjust" plan.** Step 1 (fatten the set) is already largely
achieved at the paper's `μ=15` (72% ra/safe ≈ Dubins' fat-margin regime); pushing to 89% costs
max\|τ\|≈3.5e5. Step 2 (use raised `λ` + co-design `k1` to drive `|τ|→500`) would require the SOP
to cancel **20× excess at μ=15 (10160→500)** or **700× excess at μ=50**. But §7 Dubins — with only
~16× excess at fat margin (41%) and full leverage — **could not shrink dense `|u1|`** (overfit at
70-85). The manipulator's structurally worse conditions (deg-2 `k1` vs deg-5 `Dψ`; output-only `k1`
vs velocity-dependent `Lf h`; §6.4) plus the larger excess factor make a feedback-level
bound-respecting result strictly less achievable than Dubins, not more.

**Net.** §8 corrects the §6.4 framing (manipulator baseline margin is decent, ~72% of safe set —
its fragility is in `min V`, not in set fraction) and reproduces the §4 tradeoff for the manipulator,
but **does not open a path to a bound-respecting feedback controller.** Recommendation: hold §7's
conclusion — MPC-level demonstration with **`μ = 15`** as the unconstrained baseline; the
input-constraint terminal-set framing of C1 covers what feedback-level bounding cannot deliver.

## 9. Dubins per-sample slack demo + output-space comparison (2026-05-30)

Followup to §6 / §7. A self-contained demo bundle at `scripts/dubins_demo/` exercises the
slack-SOP and hard-SOP machinery end-to-end, with the per-sample slack distribution as the
selection tool for `ub_demo` and a closed-loop output-space comparison.

**New artifacts (all under `scripts/dubins_demo/`):**
- `dubins_demo.ipynb` — 32 cells, §0 math formulation, §1–§4 SOP synthesis + slack distribution,
  §5 `ub_demo` from p80, §6 hard SOP on filtered subset, §7 output-space difference, §7.2
  4-color RA-set overlay, §8 closed-loop trajectories from a green-region $x_0$.
- `matlab/sample_n_valid.m` — iteratively samples the FL box until exactly `n_valid` lie inside
  the RA region (drops the original's "initial random count" semantics).
- `matlab/solve_k1_controller_sop_slack_persample.m` — per-sample slack solver, takes
  `slack_weight` parameter; objective is `δ + slack_weight · Σ s_{j,i}` with one scalar slack
  per (sample, channel) shared between the upper and lower bound (KKT: $s^* = \max(0, |u| - \text{ub})$).
- `matlab/build_dubins_demo.m`, `build_dubins_demo_subset.m`, `build_dubins_demo_subset_hard.m`
  — env-driven drivers (MU, XI0, UB1/UB2, N_VALID, SLACK_WEIGHT, SAMPLES_CSV, TAG).
- `/tmp/gen_dubins_demo_nb.py` + `/tmp/exec_nb_with_outputs.py` — gen script + headless
  executor that preserves outputs of unchanged code cells (source matching), runs only
  empty-output cells, and supports `--force / --force-all` for explicit re-execution.
- Notebook checks cached `meta_*.csv` against current PARAMS before deciding to call MATLAB,
  so changing `slack_weight` / `n_valid_samples` (which don't appear in the filename suffix)
  still triggers regeneration.

**Settled PARAMS for the demo:** μ=0.1, xi0=10, ub=5, n_valid_samples=250, dv=4, ds=4,
slack_weight=**100**, coverage_pct=80, sim_max_step=1e-4.

### 9.1 `slack_weight` matters: MOSEK is suboptimal at the default `w=1`

At `slack_weight=1` the per-sample slack SOP returns a numerically-suboptimal solution that
trades slack for marginally smaller `δ`. KKT diagnostic (§5.4 in the notebook):

- After the subset SOP at `(μ=0.1, xi0=10, ub_demo≈5.02)` on 157 kept samples, MOSEK
  reported `ch1 max slack = 1.05`, `δ = 108.5`.
- But `max|u_orig(x_j)| = 5.007` over the kept samples; re-using the original `k1` in the
  subset SOP would give `sum_slack = 0`, `δ = 109.08`, total objective ≈ 109.08.
- MOSEK's returned solution: `δ + Σ s = 108.5 + 1.24 = 109.74`, **0.7 unit worse**.

Sweeping `w ∈ {1, 10, 100, 1000}` at the same `(μ, xi0, ub_demo)`:

| w | `δ` | `ch1 max slack` | `dense max|u₁|` |
|---|---|---|---|
| 1 | 108.5 | 1.04 | 92 |
| 10 | 109.8 | 0.20 | 93 |
| **100** | **111.1** | **0.0005** | **74** |
| 1000 | 203.3 (blow-up) | 0.0022 | 76 |

`w=100` is the sweet spot. Beyond ~100 the SDP conditioning fails (`δ` doubles, slack stops
improving). Empirical heuristic: keep `w · sum_slack` within ~3 orders of magnitude of `δ`
(here `δ` ~10²), i.e. `w ∈ [10, 10⁴]` is the workable range.

**Side benefit observed:** at `w=100` the *dense* max|u| drops from ~92 to ~74 versus
`w=1`. Heuristic: pinning slack at KKT tightens the constraint set the polynomial `k_1`
must satisfy *at the SOP samples*, reducing the polynomial's "swing room" between samples
— a soft regularizer.

### 9.2 Honest `p80` of the slack distribution at `w=100`

With the weighted solver the §4 slack distribution becomes:

| channel | max | p99 | p90 | **p80** | p50 |
|---|---|---|---|---|---|
| u1 (ω) | 11.1 | 5.97 | 0.080 | **0.0083** | 0.008 |
| u2 (a) | 3.84 | 0.97 | 0.008 | **0.0081** | 0.008 |

(With `w=1` p80 was 0.021 — inflated by ~3.5× because MOSEK overstated slacks.)

So `ub_demo_i = 5 + p80(s_{·,i}^*) ≈ 5.008` for both channels — essentially `ub_request`
itself once the outliers are dropped.

### 9.3 Hard SOP on the kept subset is feasible

`build_dubins_demo_subset_hard.m` solves:

```
min δ  s.t.  CBF SOS,  M_j(k_1)⪰0 ∀j ∈ kept,  -ub_demo_i ≤ u_i(x_j) ≤ +ub_demo_i,  δ≥0
```

At `(μ=0.1, xi0=10, ub_demo≈5.008, 156 kept samples)`: **feasible**, `δ = 112.5`,
`feasratio = -0.96`, `pinf = dinf = 0`. So strict bounds at the 156 kept samples are jointly
achievable.

### 9.4 Output-space RA-set overlay reaffirms §6.4 / §7 negative finding

On 4000 dense samples in $X_S \setminus X_T$ (§7.2), 4-color membership for the **hard-SOP**
controller (vs vanilla `k_1`):

| | count | % |
|---|---|---|
| both certify (green) | 2238 | 55.95% |
| only unc certifies (blue) | 55 | **1.38%** |
| only con certifies (red) | 1321 | **33.02%** |
| neither | 386 | 9.65% |

Subset metrics: `subset_rate = 0.629`, `removed_frac = 0.024`, `unexpected_frac = **0.371**`.
Even with **hard SOP + filtered samples + tight `ub_demo ≈ ub`**, the constrained certificate
region is **NOT** a subset of the unconstrained one — 37% of con-cert area lies outside
unc-cert area. This **confirms** the structural negative finding from §6 / §7: the SOS
`minimize δ` objective inherently favors a `k_1` with larger `{V ≥ 0}` than vanilla, and
input constraints don't shrink that. Removal is < 3%.

### 9.5 Closed-loop simulations: controllers drive the trajectory OUT of FL region

§8 picks a green $x_0$ (both cert ≥ 0, far from $X_T$) and simulates both controllers with
`scipy.solve_ivp(DOP853, rtol=1e-8, atol=1e-10, max_step=1e-4)`. Result for one
representative `x_0 = (1.85, 1.00, 3.07, 0.39)`:

| diagnostic | con (hard) | unc (vanilla) |
|---|---|---|
| reached target | ✓ | ✓ |
| stayed_safe | ✗ | ✗ |
| **min v(t)** | 0.185 | 0.391 |
| **max v(t)** | **23.02** | **45.08** |
| max\|u₁(t)\| | 293.7 | 41.2 |
| max\|u₂(t)\| | 102.5 | 330 |
| v at peak\|u₁\| | **0.185** (FL lower bound!) | 11.14 (way outside FL) |
| % time \|u₁\| ≤ ub_demo₁ | 96.8% | 65.8% |
| % time \|u₂\| ≤ ub_demo₂ | 44.6% | 19.6% |

**Critical observation:** the FL region is $v \in [0.1, 1.0]$. Both closed loops drive `v`
to **23–45** — *two orders of magnitude past the FL upper bound*. The controllers were
synthesized only on samples inside the FL region; outside it, $u$ is unconstrained
extrapolation of the polynomial expression, plus the `1/v` factor in $A^{-1}$ spikes
whenever `v` returns near 0.

`max|u|` is **converged** across `max_step ∈ {1e-3, 1e-4}` (291 → 293, < 1% change), so the
spikes are **real** polynomial-extrapolation values, not step-size artifacts.

**Cascade:** u₂ (acceleration) is large → `v̇ = u₂` accelerates `v` to ~45 in a few seconds
→ trajectory leaves FL region → `u₁` enters extrapolation regime, including the `1/v`
spike when `v` later passes back through the lower bound.

The closed-loop simulation therefore **does not** stay inside the region the controller was
designed for. `V(x) ≥ 0` was *supposed* to keep the trajectory in the RA set (which is
inside the FL region), but `V` is a polynomial-approximate sufficient condition and the
dense closed loop can escape.

### 9.6 What this means for the paper

- The §6 / §7 / §8 cumulative finding **does not soften**: even with the cleanest possible
  per-sample slack distribution (`w=100`), the honest `ub_demo` (5.008), the strictest hard
  SOP on the filtered "easy 80%" subset, the constrained certified region still extends
  beyond the unconstrained one (`unexpected_frac ≈ 37%`).
- The closed-loop simulation is *not a fair test* of "controller respects bound" because
  the simulated trajectory leaves the controller's domain of definition (the FL region).
  Any claim like "constrained controller stays bounded along trajectories" requires either
  saturating `u` at execution time, restricting the simulation, or rejecting states whose
  closed loop exits FL.
- These are framework-level problems with the polynomial-SOS reach-avoid approach when
  combined with a sampled input-bound SOP — they are *not* fixed by any parameter tuning
  in this revision cycle. The §7 conclusion (pivot to MPC-level demonstration) stands.

## Open decisions

- **Strict bound** (`mu ≈ 0.0068`, small reach-avoid set ~21%) vs **large set**
  (`mu = 0.1`, ~89%, empirical/trajectory-based bound claim) — for the paper.
- ~~Whether `lambda` / certificate-shaping can enlarge the at-±5 region~~ — **answered (§6): no.**
  Raising `λ` does not reduce `|u|`; it overfits (Dubins) or collapses the reach-avoid set (manip).
- ~~Whether a feedback-level constrained-vs-unconstrained separation is demonstrable (Dubins,
  `lambda` in 50x-100x `mu`)~~ — **answered (§7): no.** Pivot the input-constraint demonstration
  to the MPC level; reframe C1 as a terminal-set contribution. **(decision pending: adopt the
  MPC-level demo / reframe C1 in the manuscript.)**
- ~~Whether a larger μ for the manipulator can fatten the set enough for k1/λ adjustment to bound
  `|τ|`~~ — **answered (§8): no.** Larger μ does fatten the set (monotonic 54%→100%) but the
  μ↔|τ| tradeoff is super-linear (max|τ| reaches `1e7` by μ=200), and §7's overfit wall is more
  severe at the manipulator's `|τ|` scales than for Dubins. Keep `μ=15` as the baseline.
- Controller files in `controllers/` are kept at the committed (paper) versions; the
  FL-region edits in the example scripts are **not yet re-exported** to controllers.
- ~~Whether the per-sample slack SOP can deliver a controller that respects `ub` along
  closed-loop trajectories~~ — **answered (§9): no.** Even with `slack_weight = 100`
  (KKT-honest slack), `coverage_pct = 80` (`ub_demo ≈ ub`), filtered subset + hard SOP,
  and `max_step = 1e-4` (converged), closed-loop drives `v` to 23–45 (FL upper bound is
  1.0), pushing the controller into extrapolation where `1/v` spikes |u| to ~300. The SOS
  certificate `V ≥ 0` is a polynomial-approximate sufficient condition; it does *not* in
  practice prevent the closed loop from leaving the FL region.
