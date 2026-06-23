# P1 — Franka Reach-Avoid in MuJoCo (high-fidelity testbed)

Status: **Draft for review** · Date: 2026-06-23 · Repo: `reach_avoid_backstepping_MPC`

## 1. Goal & scope

Build the first closed-loop milestone (P1) of a MuJoCo high-fidelity testbed that
demonstrates the reach-avoid backstepping + SOS / RA-MPC algorithm executing on a
realistic Franka Panda, as a stepping stone toward real lab hardware.

**P1 in scope:**
- Franka end-effector modeled as a **planar task-space double integrator** (validated
  abstraction, §4).
- **Oracle** (ground-truth) safe/target/obstacle geometry — no vision yet.
- **Per-scene** certificate synthesis (parameterized MATLAB SOS) for whatever
  oracle sets are chosen.
- Closed loop in MuJoCo comparing **RA-MPC vs Vanilla MPC vs closed-form RA**.
- Metrics (reach-avoid success, trajectories, control effort) + rendered video.

**Explicitly out of P1** (later phases): vision-built sets (P3), full 3D motion +
grasping (P2), TurtleBot (P4), real hardware (P5). See
`memory/project-mujoco-reachavoid-testbed-plan.md`.

## 2. Background

The repo synthesizes, for given dynamics `ẋ = f(x)+g(x)u`, output `y=h(x)`, safe set
`{ψ(y)≥0}`, target `{φ(y)≤0}`, bounded `u`:
- a reach-avoid **certificate** `V(x)` and a closed-form bounded controller `u*(x)`
  (backstepping + SOS, MATLAB; exported to `controllers/*.py` as SymPy);
- an **RA-MPC** that uses `V(x_N)≥0` as terminal constraint (CasADi/IPOPT);
- a **Vanilla MPC** baseline with target terminal constraint `φ(h(x_N))≤0`.

The current Python "simulation" steps the *idealized* model (`F_disc`, RK4). P1's core
move is to **replace idealized stepping with MuJoCo physics** while the controller's
internal model stays idealized — i.e., a model-mismatch + ZOH-discrete-control test.

**Hardware abstraction (decided):** operate at the device's exposed interface, not its
low-level dynamics. Franka → task-space Cartesian control; the 7-DoF arm dynamics,
redundancy, and joint limits are handled by the robot's own controller (an OSC layer
in sim). Our reach-avoid law only reads EE state and commands task-space acceleration.

## 3. Validated foundation (already built)

The P0.5 abstraction-fidelity spike (see `memory/project-franka-di-abstraction-validated.md`)
established and shipped:
- `mujoco_sim/plant.py` — `FrankaPlanarPlant`: loads Menagerie Franka, OSC inner loop
  (task-space acceleration interface) with gravity/Coriolis **and joint-damping**
  compensation, regulates EE height `z0` + orientation + nullspace posture; outer
  command `u` held ZOH at `ctrl_hz`, inner OSC recomputed each physics step.
- `mujoco_sim/scene_def.py` — parameterized planar scene (workspace, target ellipse,
  obstacles) + numpy `safe_psi(y)≥0`, `target_phi(y)≤0`; emits `scene_reachavoid.xml`.
- `mujoco_sim/view_sim.py` — live viewer + offscreen mp4 recording.
- `mujoco_sim/probe_abstraction.py` — the fidelity probe.

**Validated result:** the planar double integrator is a faithful abstraction **inside
the dexterous workspace** (gain ≈1.00, decoupled, sub-mm tracking vs ideal DI), and
**degrades toward the +x arm-extension boundary** (wrist ±12 N·m saturates). Home EE is
at `(0.554, 0, 0.625)`; the dexterous box is roughly `x∈[0.30,0.62], y∈[-0.28,0.28]` at
`z0=0.625`.

**Constraints this imposes on P1 (hard):**
1. Site the task workspace + all sets inside the dexterous region (`x ≲ 0.62`).
2. Bound the control acceleration `a_max` to what the wrist sustains in-envelope.
3. Keep joint-damping compensation in the plant.

## 4. System model (P1)

Planar task-space double integrator at fixed height `z0`:
- State `x = [p_x, p_y, v_x, v_y]` (EE planar position + velocity).
- Dynamics `ṗ = v`, `v̇ = u`; control `u = (u_x,u_y)` = task-space acceleration,
  bounded `‖u‖∞ ≤ a_max`.
- Control-affine: `f = [v_x,v_y,0,0]ᵀ`, `g = [0;0;I₂]`; output `y = h(x) = (p_x,p_y)`.
- No trigonometric terms → SOS synthesis is simpler than the existing Dubins/manipulator
  examples; relative degree 2 (position ← velocity ← acceleration) suits backstepping.

`z0`, EE orientation, and the redundant DoF are regulated by the plant's OSC (not part
of the reach-avoid state).

## 5. Set construction (a single polynomial each)

Both sets are represented as a SINGLE polynomial in output space `y`, matching the
existing MATLAB pipeline (which takes one safe polynomial and one target polynomial).
A preprocessing step in `scene_def.py` COMPOSES them from the scene:
- **Target** `φ(y) ≤ 0`: the goal ellipse (a single quadratic) —
  `compose_target_poly(scene) -> φ(y)`.
- **Safe set** `ψ(y) ≥ 0`: composed as a PRODUCT
  `ψ(y) = w(y) · ∏_i o_i(y)`, where `w(y) = 1 − ((y₁−c_x)/a_x)² − ((y₂−c_y)/a_y)² ≥ 0`
  is a single **elliptical** workspace factor and
  `o_i(y) = (y−c_i)ᵀ(y−c_i) − r_i² ≥ 0` outside obstacle i —
  `compose_safe_poly(scene) -> ψ(y)`. The workspace is an ellipse, NOT a box: a box as
  a product of two slabs `(a_x²−(y₁−c_x)²)(a_y²−(y₂−c_y)²)` has spurious positive lobes
  in the diagonal exterior corners, breaking exactness.

**Why the product is an exact single polynomial here:** if every obstacle is pairwise
disjoint and lies strictly inside the workspace, then at any point at most one factor is
negative, so `ψ ≥ 0` holds *exactly* on the safe region (no spurious positive lobes).
`scene_def` MUST validate this (disjoint + interior) when composing. Degree grows as
`2(k+1)` for k obstacles (k=2 → degree 6), within SOS range — the Dubins example already
uses a degree-8 product safe set.

Geometry is decorative in MuJoCo (`contype/conaffinity=0`): obstacles are enforced by
the controller's safe-set constraint, not physical contact.

**CRITICAL PRECONDITION (learned from the Milestone-0 run, 2026-06-23).** The
backstepping reach-avoid construction requires the safe-set / seed function to be
**navigation-function-like**: it must attain its global maximum at a point **inside the
target set**, and have **no other critical points** (∇=0) inside the safe set. The plain
product `ψ = w·∏ᵢoᵢ` does NOT satisfy this — its maximum sits mid-workspace and it has
interior saddles/extrema between obstacles — so the synthesised `k1` flow converges to
the wrong place. The set-construction step must therefore RESHAPE the function (see the
Milestone-0 findings in the plan doc) so its unique interior maximum lies in the target.

### 5.1 Milestone 0 — set & vector-field verification (GATES everything else)

Per the agreed sequencing, before any double-integrator / MPC / MuJoCo work we verify in
2D that a scene yields a sound reach-avoid design. Only when this passes do we build the
rest.
1. **Set check (2D):** plot the zero-level sets of composed `ψ(y)` and `φ(y)` over the
   workspace; confirm `{ψ≥0}` = workspace minus obstacles and `{φ≤0}` = the goal ellipse,
   and that the disjoint/interior validation holds.
2. **Single-integrator field check (2D):** synthesize for the scene, then plot the
   output-space single-integrator field `ẏ = k1(y)` (the backstepping virtual control) on
   a grid, overlaid with the sets; confirm streamlines flow into the target while staying
   in `{ψ≥0}`. For the **holonomic** double integrator the output is fully actuated through
   the cascade, so the `k1` field should itself be a valid reach-avoid flow.
   - Caveat: for the *nonholonomic* Dubins example the exported k1 field can leave the safe
     set — there the guarantee lives in `V`/full `u`, not k1 (see
     `memory/project-dubins-k1-flow-not-reachavoid`). The DI is the well-behaved case, and
     this check is exactly what confirms it.

This milestone also subsumes "verify the MATLAB pipeline generalizes to the double
integrator" (step 2 requires a successful synthesis). Deliverable: a `verify_sets.py` /
`verify_field.py` producing the two 2D figures.

## 6. Per-scene synthesis (parameterized)

New `matlab/example_franka_planar.m`: takes the DI dynamics + chosen sets + `a_max`,
runs `reach_avoid_controller.m` (backstepping) + `solvesop_bounded_control.m` (SOS
bounding), exports `u_opt / certificate_opt / k1_opt` to `controllers/`.

Driven from Python by extending `python/matlab_runner.py`:
- `mujoco_sim/synth.py`: `synthesize(scene) -> ControllerBundle`, hashing the scene
  (sets + a_max + model degrees) to a cache key; re-runs MATLAB only on cache miss.
- Output bundle: SymPy `u_opt`, `certificate_opt`, `k1_opt` (+ metadata).

**Risks:** SOS may be infeasible for some set choices; synthesis takes ~minutes; sets
must be polynomial. The synthesis driver must surface infeasibility clearly (no silent
fallback). Verifying the pipeline generalizes to the planar double integrator (it has
only run on Dubins/manipulator) is part of **Milestone 0** (§5.1, step 2).

## 7. Controllers (decoupled interface)

Unified interface, independent of the simulator:

```
class Controller:
    def reset(self): ...
    def compute(self, x, t) -> u   # x = [px,py,vx,vy] measured; returns task-space accel u
```

Three implementations, all consuming a synthesized bundle for the active scene:
- **RA (closed form):** evaluate `u_opt(x)`. No online optimization.
- **RA-MPC:** CasADi/IPOPT NLP; internal model = idealized DI (RK4); terminal
  `V(x_N) ≥ 0`; safe-set inequalities at each node; `‖u‖∞ ≤ a_max`; output-tracking cost.
- **Vanilla MPC:** same NLP with target terminal `φ(h(x_N)) ≤ 0`, no certificate.

Closed loop: read EE state from `FrankaPlanarPlant` → `Controller.compute` → task-space
`u` → `plant.step(u)` (OSC realizes it) → repeat at `ctrl_hz`.

## 8. Module layout, interfaces, data flow

```
mujoco_sim/
  plant.py          [built]  FrankaPlanarPlant: OSC task-space DI interface + mismatch knobs
  scene_def.py      [built]  Scene config + safe/target evaluators + MJCF generation
                             [extend] compose_safe_poly / compose_target_poly (single-poly
                             composition) + disjoint/interior validation
  view_sim.py       [built]  live viewer + mp4 recording (will replay closed loop)
  probe_abstraction.py [built] abstraction-fidelity probe
  verify_sets.py    [new]    Milestone 0: 2D level-set plot of composed ψ/φ
  verify_field.py   [new]    Milestone 0: 2D single-integrator k1 vector-field plot
  synth.py          [new]    drive per-scene MATLAB SOS synthesis (extends matlab_runner) + cache
  controllers.py    [new]    RA / RA-MPC / Vanilla, unified compute(x,t)->u (CasADi for MPC)
  sim_loop.py       [new]    closed-loop runner: plant x controller -> trajectory + success flags
  sweep.py          [new]    batch over initial conditions / scenes / mismatch knobs -> metrics
  metrics.py        [new]    reach-avoid success, safety margin, control effort, saturation
  viz.py            [new]    output-space trajectory plots (path + sets) + video hooks
matlab/
  example_franka_planar.m [new]  parameterized DI + sets -> synthesis -> export
```

**Unit boundaries (each independently testable):**
- `Scene` (data + set evals) — no sim/MATLAB deps.
- `synth` (Scene → controller bundle) — MATLAB dep, cached, side-effect = exported file.
- `Controller` (bundle + measured x → u) — no MuJoCo dep; testable on the ideal DI.
- `FrankaPlanarPlant` (u → physics + measured x) — no controller/SOS dep.
- `sim_loop` composes the above; `metrics`/`viz` consume its trajectory.

**Data flow:** `scene_def → synth (MATLAB) → bundle → Controller`; `Controller ↔ FrankaPlanarPlant`
closed loop in `sim_loop`; trajectory → `metrics`/`viz`/`view_sim`.

## 9. Decisions still needed (for review)

- **A. Set representation (RESOLVED per review).** Safe & target are each a SINGLE
  polynomial. Safe = `workspace · ∏(outside-obstacle)`, exact because obstacles are
  pairwise disjoint & interior to the workspace (`scene_def` validates this); no
  multi-inequality SOS extension needed. Composition + the 2D set/field checks are
  **Milestone 0** (§5.1) and gate everything else.
- **B. `a_max` value.** Pick from the probe envelope (in-envelope a0=1 had ~45% torque
  peak; propose `a_max = 1.0 m/s²` to start, tighten if boundary scenes saturate).
- **C. Task-plane height.** Keep `z0 = 0.625` (validated EE-home height, geometry floats)
  for P1; move to a real tabletop in P2 (needs re-validating dexterity at that height).
  (Recommended: keep 0.625 for P1.)
- **D. MPC tuning for the DI.** Horizon `N`, `δt`, output weights `Q_y/Q_{f,y}`, `R_u`,
  IPOPT tolerances — start from the repo's existing values and adapt.
- **E. Scene difficulty / success criteria.** How many initial conditions, what counts
  as the headline result (e.g., RA-MPC success rate ≥ X% while Vanilla fails on a
  fraction of obstacle scenes).

## 10. Validation strategy

- **Milestone 0 (gate, do first):** composed `ψ`/`φ` 2D level-set figure is correct;
  synthesis succeeds for the DI; single-integrator `k1` field flows into the target and
  stays in `{ψ≥0}` (§5.1). Nothing else starts until this passes.
- **Unit:** set evaluators (point in/out) + disjoint/interior validation, OSC self-check
  (already passing: u=0 → EE holds), controller on ideal DI (reaches target, respects
  bounds).
- **Integration:** closed-loop rollout in MuJoCo from several initial EE positions in a
  dexterous scene; assert reach (enter target) ∧ safe (ψ≥0 throughout); compare the
  three controllers; record video.
- **Mismatch sweeps:** vary control rate, sensor noise, mass/damping perturbation;
  report success-rate degradation.

## 11. Risks

| Risk | Mitigation |
|---|---|
| SOS infeasible / slow per scene | Milestone 0 catches it early per scene; cache; surface infeasibility |
| Pipeline never run on a double integrator | Milestone 0 step 2 verifies it on the DI |
| Product safe set has spurious lobes (obstacles overlap / exit workspace) | `scene_def` validates disjoint + interior; Milestone 0 visual check |
| `k1` field not reach-avoid (as in nonholonomic Dubins) | Milestone 0 field check; DI is holonomic, expected clean |
| Kinematic-boundary saturation | keep sets in dexterous box; bound `a_max` |
| RA-MPC not real-time | report solve times; tune N/δt; warm-start |
| Floating task plane looks unnatural | acceptable for P1; tabletop in P2 |

## 12. Out of scope (future)

P2 3D + simplified grasp · P3 vision-built sets · P4 TurtleBot (unicycle via `cmd_vel`)
· P5 real hardware (interface already aligned).
