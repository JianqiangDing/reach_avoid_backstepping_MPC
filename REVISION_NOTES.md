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
The regularization was reverted.

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

## Open decisions

- **Strict bound** (`mu ≈ 0.0068`, small reach-avoid set ~21%) vs **large set**
  (`mu = 0.1`, ~89%, empirical/trajectory-based bound claim) — for the paper.
- Whether `lambda` / certificate-shaping can enlarge the at-±5 region beyond ~21%.
- Controller files in `controllers/` are kept at the committed (paper) versions; the
  FL-region edits in the example scripts are **not yet re-exported** to controllers.
