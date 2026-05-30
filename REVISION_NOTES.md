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

## Open decisions

- **Strict bound** (`mu ≈ 0.0068`, small reach-avoid set ~21%) vs **large set**
  (`mu = 0.1`, ~89%, empirical/trajectory-based bound claim) — for the paper.
- ~~Whether `lambda` / certificate-shaping can enlarge the at-±5 region~~ — **answered (§6): no.**
  Raising `λ` does not reduce `|u|`; it overfits (Dubins) or collapses the reach-avoid set (manip).
- ~~Whether a feedback-level constrained-vs-unconstrained separation is demonstrable (Dubins,
  `lambda` in 50x-100x `mu`)~~ — **answered (§7): no.** Pivot the input-constraint demonstration
  to the MPC level; reframe C1 as a terminal-set contribution. **(decision pending: adopt the
  MPC-level demo / reframe C1 in the manuscript.)**
- Controller files in `controllers/` are kept at the committed (paper) versions; the
  FL-region edits in the example scripts are **not yet re-exported** to controllers.
