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
far too weak; §5 confirms this with the strongest incentive (minimize the control slack).

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

## Open decisions

- **Strict bound** (`mu ≈ 0.0068`, small reach-avoid set ~21%) vs **large set**
  (`mu = 0.1`, ~89%, empirical/trajectory-based bound claim) — for the paper.
- Whether `lambda` / certificate-shaping can enlarge the at-±5 region beyond ~21%.
- Controller files in `controllers/` are kept at the committed (paper) versions; the
  FL-region edits in the example scripts are **not yet re-exported** to controllers.
