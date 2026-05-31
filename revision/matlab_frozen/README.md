# `revision/matlab_frozen/` — frozen MATLAB toolchain (minimal)

**Frozen copies** of only the MATLAB library functions the revision experiments
actually call, vendored here so revision is self-contained: revision MATLAB drivers
`addpath` **only this directory** (plus their own `revision/<phase>/`), never the
main repo's `matlab/`. A change in the main repo cannot alter revision results.

Kept **minimal**: only the closed dependency set of the current phase is vendored,
not the whole toolchain. Files are added when a later phase first needs them.

Snapshot from `matlab/` and `scripts/dubins_demo/matlab/`, 2026-05-31.

## Phase 1.2 (per-sample slack diagnostic) — closed set

| File | Role |
|---|---|
| `reach_avoid_controller.m` | backstepping → `u`, `k1`, certificate, `cert_term_dict`, `p`, `r_deg` (self-contained: all helpers are local subfunctions) |
| `solve_vanilla_k1_controller_xi.m` | vanilla `k1` + `lambda` at leverage `xi0` (calls `poly2sym`, `sym2pvar`) |
| `sample_n_valid.m` | sample N valid reach-avoid states (safe, outside target, certificate ≥ 0) |
| `solve_k1_controller_sop_slack_persample.m` | **per-sample slack SOP**: each (sample j, channel i) slack `s_{j,i}`, min `delta + w·Σ s_{j,i}`; reach-avoid stays hard (calls `detect_trigonometric_terms`, `sym2pvar`, `poly2sym`) |
| `detect_trigonometric_terms.m` | find trig terms to dummy-substitute for SOS |
| `sym2pvar.m` | sym → SOSTOOLS `pvar` |
| `poly2sym.m` | SOSTOOLS poly → sym |

Not vendored (not needed by Phase 1.2): controller export, hard / per-channel
solvers, bound estimators, plotting, decoupling check. Phase 1.3 will add the
hard-SOP solver and `export_to_python.m` when it lands.

External deps (SOSTOOLS, Mosek, MATLAB Symbolic Toolbox) are environment packages.
