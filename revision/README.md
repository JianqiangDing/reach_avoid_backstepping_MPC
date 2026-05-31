# `revision/` — L-CSS revision experiment workspace

This directory holds all new experiment code and outputs for the L-CSS revision of paper 26-0694. It is **fully self-contained**: every experiment reads and writes only files under `revision/`, and does **not** read the main repo's `data/` or import/call its `python/` / `matlab/` / `controllers/` at runtime. The main repo is left untouched.

Top-level plan: `../EXPERIMENT_AUDIT.md`.

## Subdirectories

| Path | Phase | Contents |
|---|---|---|
| `decoupling_check/` | Phase 1.1 | A(x) invertibility check (consumes main-repo `data/decoupling_*.csv`, writes `revision/data/phase1_1_outputs_*.json`) |
| `slack_diagnostic/` | Phase 1.2 | per-sample slack diagnostic SOP (selects `u_max_eff` + filters the sample set; **does not produce the final controller**) |
| `controller_synthesis/` | Phase 1.3 | hard-constrained SOP (no slack) synthesizing the final $u^*, V, k_1$, exported to `controllers/` |
| `controllers/` | Phase 1.3 output | revision-resynthesized `.py` controllers (the only controller input for Phase 2/3 experiments) |
| `success_module/` | Phase 2 | unified closed-loop success module (implemented after Phase 1 completes) |
| `experiments/` | Phase 3 | the §4.2–4.6 experiment notebooks (implemented after Phase 2 completes) |
| `data/` | all-phase outputs | npz / csv / mat / json / png |

## Self-containment principle

Revision experiments must own their code and data so that changes in the main repo cannot silently affect revision results (which would undermine the reliability of every downstream experiment). Therefore:

- **No main-repo data reads** at runtime (no `data/decoupling_*.csv`, etc.). Inputs are generated inside `revision/` (e.g. samples drawn in-notebook).
- **No main-repo code imports/calls** at runtime. Needed logic is re-derived or re-implemented inside `revision/`:
  - A(x) is derived symbolically in-notebook and a Python evaluator is validated against it (Phase 1.1);
  - for Phase 1.2/1.3, the SOP solvers and any paper controllers needed for reproduction are kept as **frozen copies** under `revision/`, not called from main-repo `matlab/` / `controllers/`.
- **Correctness is validated against an in-notebook source of truth** (the symbolic derivation), not against a main-repo CSV.

## Run conventions

- Each notebook resolves the repo root only to locate `revision/data/` for outputs; it does not read other top-level directories.
- All revision outputs go to `revision/data/`.
- The main repo is never modified.
