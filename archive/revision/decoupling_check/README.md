# Phase 1.1 — Decoupling matrix A(x) invertibility check

## Purpose

Confirm that the backstepping FL coordinates are uniformly valid over the synthesis region. Locations where $A(\bm x)$ is near-singular make the SOP-synthesized $u^*(\bm x)$ blow up whenever the closed loop touches them (Dubins has $1/v$, manipulator has $1/\sin x_2$).

## Self-contained

These notebooks are **fully self-contained**: they do not read main-repo `data/` or import main-repo `python/` / `matlab/` at runtime. Each notebook:

- defines the system ($f, g, h$, safe/target sets) in-notebook;
- derives $A(\bm x)$ and $\det A(\bm x)$ symbolically (SymPy);
- provides a vectorized Python evaluator `A_metrics` and **validates it against the in-notebook symbolic $A$** (NumPy SVD), so correctness does not depend on any external file;
- generates its own samples (NumPy, fixed seed) and writes all outputs to `revision/data/`.

(The earlier "consume the main-repo `data/decoupling_*.csv`" approach was dropped: revision experiments must own their code and data so that changes in the main repo cannot silently affect revision results.)

## The two notebooks

| Notebook | Output |
|---|---|
| `dubins_detA.ipynb` | `revision/data/phase1_1_outputs_dubins.json`, `detA_dubins_output_space.png` |
| `manipulator_detA.ipynb` | `revision/data/phase1_1_outputs_manipulator.json`, `detA_manipulator_output_space.png` |

Each notebook is structured as: §0 parameter table, §1 system + $A$ derivation, §1-SymPy cell, §2 setup + `A_metrics` + validation against symbolic $A$, §3 invertibility over $\mathcal X_S\setminus\mathcal X_T$, §4 output-space coverage colored by $|\det A|$, §5 verdict, §6 export JSON.

### Manipulator region relaxation

The manipulator notebook reports the **relaxed synthesis region** $\mathcal X_S^{\rm eff}$ (elbow margin $0.3$, full shoulder range $[-\pi,\pi]$), which covers ~99% of the reachable safe set vs ~17% for the original paper box $x_2\in[0.8,\pi-0.8]$, while keeping $\sigma_{\min}(A)\approx 7\epsilon$. Dubins needs no relaxation ($\det A = -v$ is position-independent and $v\in[0.1,1]$ spans the full band).

## Verdict

Threshold $\epsilon = 10^{-2}$:

- $\min_{\mathcal X_S\setminus\mathcal X_T}\sigma_{\min}(A) \ge \epsilon$ → **PASS**;
- otherwise → **NEEDS_SHRINK** (tighten the region toward the singular-direction margin).

## Phase 1.1 exit contract

One JSON per example with `verdict`, `epsilon`, `self_contained`, `validated_against`, `stats`, and `X_S_eff_def`. Phase 1.2 reads only the `X_S_eff_def` field.
