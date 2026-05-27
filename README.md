# Reach-Avoid Backstepping MPC

<img src="reach_avoid_mpc.png" width="800" alt="Reach-Avoid MPC"/>

## Overview

This repository implements a two-phase framework for synthesizing controllers that drive a nonlinear system into a _target set_ while staying inside a _safe set_, subject to bounded control inputs.

Given system dynamics $\dot{x} = f(x) + g(x)u$, an output map $y = h(x)$, a safe set $\{\psi(y) \geq 0\}$, and a target set $\{\phi(y) \leq 0\}$, the MATLAB scripts use backstepping and Sum-of-Squares (SOS) polynomial optimization to compute:

- A reach-avoid **certificate** $V(x) \geq 0$ certifying that any trajectory starting with $V(x_0) \geq 0$ will reach the target without leaving the safe set.
- A **bounded controller** $u^*(x) \in [u_\mathrm{lb},\, u_\mathrm{ub}]$ that makes $V$ a valid Lyapunov-like certificate.

Control input bounds are enforced via Scenario Optimization Programming (SOP).

Three methods are compared across all example systems:

| Method                        | Description                                                                                                                                                  |
| ----------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Unconstrained Reach-Avoid** | Closed-form certificate controller $u^*(x)$ synthesised by MATLAB; operates entirely inside the certificate region $\{V(x)\geq 0\}$; no online optimisation. |
| **Reach-Avoid MPC**           | MPC with certificate terminal constraint $V(x_N)\geq 0$; extends guaranteed reach-avoid coverage to states outside the certificate region ($V(x_0)<0$).      |
| **Vanilla MPC**               | Standard tracking MPC with hard target-set terminal constraint $\phi(h(x_N))\leq 0$; no safety certificate; baseline for comparison.                         |

## Systems

| System                   | State dim | Inputs                                  | Notebook prefix               |
| ------------------------ | --------- | --------------------------------------- | ----------------------------- |
| Double integrator        | 2         | scalar $u \in [-20, 20]$                | `example_double_integrator_*` |
| Dubins car               | 4         | angular rate $\omega$, acceleration $a$ | `example_dubins_car_*`        |
| 2-DoF planar manipulator | 6         | joint torques $\tau_1, \tau_2$          | `example_manipulator_*`       |

## Repository Structure

```
.
├── README.md  reach_avoid_mpc.png  .gitignore
│
├── matlab/                           # MATLAB synthesis pipeline
│   ├── reach_avoid_controller.m      #   Main backstepping certificate synthesis
│   ├── solvesop_bounded_control.m    #   Bounded control via SOP + SOS
│   ├── solve_k1_controller_sop.m     #   k1 controller SOP solve
│   ├── solve_vanilla_k1_controller.m #   Unconstrained (vanilla) k1 controller
│   ├── export_to_python.m            #   Export symbolic results → controllers/*.py
│   ├── example_double_integrator.m   #   Double integrator
│   ├── example_dubins_car.m          #   Dubins car
│   ├── example_manipulator.m         #   2-DoF planar manipulator
│   └── …                             #   Utilities (poly2sym, sym2pvar, plotPolyL, …)
│
├── python/                           # Shared Python helpers
│   ├── functional.py                 #   Simulation utilities and color helpers
│   └── acrobot_inverse_kinematics.py #   IK solver for the 2-link manipulator
│
├── controllers/                      # Exported controllers (active, imported by notebooks)
│   ├── sop_bounded_control_ex1_debug_20260314_213424.py        # double integrator MPC
│   ├── sop_bounded_control_dubins_car_result_20260316_211504.py        # Dubins RA-MPC
│   ├── sop_bounded_control_unconstrained_controller_20260316_211252.py # Dubins unconstrained
│   ├── sop_bounded_control_acrobot_result_20260317_222858.py           # manipulator RA-MPC
│   └── k1_acrobot_cdc2026.py                                           # manipulator k1 controller
│
├── notebooks/                        # Experiment notebooks (run from this dir)
│   ├── example_double_integrator_no_mpc.ipynb  # Unconstrained Reach-Avoid
│   ├── example_double_integrator_mpc.ipynb     # Reach-Avoid MPC
│   ├── example_dubins_car_unconstrained_reach_avoid.ipynb   # Unconstrained Reach-Avoid
│   ├── example_dubins_car_reach_avoid_mpc.ipynb             # Reach-Avoid MPC
│   ├── example_dubins_car_vanilla_mpc.ipynb                 # Vanilla MPC
│   ├── example_manipulator_unconstrained_reach_avoid.ipynb  # Unconstrained Reach-Avoid
│   ├── example_manipulator_reach_avoid_mpc.ipynb            # Reach-Avoid MPC
│   └── example_manipulator_vanilla_mpc.ipynb                # Vanilla MPC
│
├── data/                             # Cached trajectories (loaded by notebooks)
│   ├── traj_controls_reach_avoid_mpc.npz
│   ├── traj_controls_vanilla_mpc.npz
│   └── traj_controls_ra_mpc_acrobot.npz
│
├── scripts/                          # Standalone experiment scripts (see scripts/README.md)
└── archive/                          # Superseded material
    └── example_manipulator_cdc2025.ipynb
```

Each notebook starts with a small bootstrap cell that puts `controllers/` and
`python/` on `sys.path` and points `DATA` at `data/`, so notebooks run correctly
from the `notebooks/` directory.

## Workflow

```
MATLAB
  example_XX.m
      │  defines f, g, h, safe set, target set, control bounds
      ▼
  reach_avoid_controller.m            ← backstepping design
      │  produces symbolic k1(x)
      ▼
  solvesop_bounded_control.m          ← SOP + SOS bounding
      │  produces u*(x), V(x)
      ▼
  export_to_python.m
      │  writes  controllers/sop_bounded_control_XXX_<timestamp>.py
      ▼
Python
  example_XX_unconstrained_reach_avoid.ipynb  ← Unconstrained Reach-Avoid
  example_XX_reach_avoid_mpc.ipynb            ← Reach-Avoid MPC
  example_XX_vanilla_mpc.ipynb               ← Vanilla MPC
```

Each exported Python file contains three SymPy expressions:

```python
u_opt         # list[2]  — bounded controller [u1, u2]
certificate_opt  # Expr   — reach-avoid certificate V(x)
k1_opt        # list[2]  — backstepping intermediate signal
```

## Dependencies

### MATLAB

- MATLAB R2022a+
- Symbolic Math Toolbox
- [SOSTOOLS](https://github.com/oxfordcontrol/SOSTOOLS) (SOS programming)

### Python

```
python >= 3.10
sympy
numpy
scipy
casadi      # MPC NLP solver (IPOPT backend)
matplotlib
```

Install Python dependencies:

```bash
pip install sympy numpy scipy casadi matplotlib
```

Or with conda:

```bash
conda env create -f environment.yml   # if provided
conda activate rab_mpc
```

## Running the Notebooks

1. Run the desired MATLAB example script (from `matlab/`) to generate the
   controller Python file; the export is written to `controllers/`:

   ```matlab
   % in MATLAB, from the matlab/ directory
   example_double_integrator   % synthesises controller and exports to controllers/
   example_dubins_car
   example_manipulator
   ```

2. Open the corresponding Jupyter notebook (from `notebooks/`) and run all cells:
   ```bash
   jupyter notebook notebooks/example_double_integrator_mpc.ipynb
   ```

The notebooks are self-contained after the controller Python file exists; the
bootstrap cell at the top resolves `controllers/`, `python/`, and `data/`.

## Experiment Settings

### System & Set Parameters

| System             | State $x$                                     | Output $y = h(x)$                                                                   | Safe set                                            | Target set                                                                                 | Control bounds                             |
| ------------------ | --------------------------------------------- | ----------------------------------------------------------------------------------- | --------------------------------------------------- | ------------------------------------------------------------------------------------------ | ------------------------------------------ |
| Double integrator  | $[x_1, x_2]$                                  | $y = x_1$                                                                           | $1 - x_1^2 \geq 0$                                  | $x_1^2 - 0.01 \leq 0$                                                                      | $u\in[-20,20]$                             |
| Dubins car         | $[x_1, x_2, \theta, v]$                       | $[x_1, x_2]$                                                                        | Annular ring: $4 \leq y_1^4+y_2^4 \leq 16$          | Ellipse centred at $(-1.7, 0)$: $\frac{(y_1+1.7)^2}{0.1}+\frac{y_2^2}{0.4}\leq 1$          | $\omega\in[-5,5]$ rad/s, $a\in[-5,5]$ m/s² |
| Planar manipulator | $[q_1, q_2, \dot{q}_1, \dot{q}_2, q_1{+}q_2]$ | EE position $(l_1\cos q_1 + l_2\cos(q_1{+}q_2),\ l_1\sin q_1 + l_2\sin(q_1{+}q_2))$ | Workspace polytope ($\vert\sin q_2\vert \geq 0.15$) | Ellipse centred at $(5.5, 1.8)$: $\frac{(y_1-5.5)^2}{2.88}+\frac{(y_2-1.8)^2}{0.32}\leq 1$ | $\tau_{1,2}\in[-500,500]$ N·m              |

Physical constants — manipulator: $m_i=1$ kg, $l_i=4$ m, $I_i=0.02$ kg·m² ($i=1,2$).

### Controller Settings Per Notebook

| Notebook                                        | Method                    | $\delta t$ (s) | $N$ | $Q_y$                 | $Q_{f,y}$              | $R_u$            | Terminal constraint   |
| ----------------------------------------------- | ------------------------- | :------------: | :-: | --------------------- | ---------------------- | ---------------- | --------------------- |
| `example_double_integrator_no_mpc`              | Unconstrained Reach-Avoid |      0.01      |  —  | —                     | —                      | —                | $V(x)\geq 0$ (always) |
| `example_double_integrator_mpc`                 | Reach-Avoid MPC           |      0.01      | 20  | $Q=\text{diag}(10,1)$ | $Q_f=\text{diag}(1,1)$ | $R=0.01$         | $V(x_N)\geq 0$        |
| `example_dubins_car_unconstrained_reach_avoid`  | Unconstrained Reach-Avoid |      0.02      |  —  | —                     | —                      | —                | $V(x)\geq 0$ (always) |
| `example_dubins_car_reach_avoid_mpc`            | Reach-Avoid MPC           |      0.05      | 25  | diag(5, 5)            | diag(80, 80)           | diag(0.05, 0.05) | $V(x_N)\geq 0$        |
| `example_dubins_car_vanilla_mpc`                | Vanilla MPC               |      0.05      | 25  | diag(5, 5)            | diag(80, 80)           | diag(0.05, 0.05) | $\phi(h(x_N))\leq 0$  |
| `example_manipulator_unconstrained_reach_avoid` | Unconstrained Reach-Avoid |      0.01      |  —  | —                     | —                      | —                | $V(x)\geq 0$ (always) |
| `example_manipulator_reach_avoid_mpc`           | Reach-Avoid MPC           |      0.01      | 20  | diag(5, 5)            | diag(80, 80)           | diag(1e-3, 1e-3) | $V(x_N)\geq 0$        |
| `example_manipulator_vanilla_mpc`               | Vanilla MPC               |      0.01      | 20  | diag(5, 5)            | diag(50, 50)           | diag(1e-3, 1e-3) | $\phi(h(x_N))\leq 0$  |

All MPC problems solved with IPOPT (tolerance $10^{-4}$, max 3000 iterations) via CasADi.

## MPC Formulation

The MPC solves at each step $t$ (Dubins car and manipulator use output-tracking cost; double integrator uses state cost):

$$
\min_{x_{0:N},\, u_{0:N-1}} \sum_{k=0}^{N-1} {\lVert h(x_k) - y^{\ast}\rVert}_{Q_y}^2 + {\lVert u_k\rVert}_{R_u}^2 + {\lVert h(x_N) - y^{\ast}\rVert}_{Q_{f,y}}^2
\quad \text{(output cost, Dubins / manipulator)}
$$

$$
\min_{x_{0:N},\, u_{0:N-1}} \sum_{k=0}^{N-1} x_k^{\top} Q x_k + u_k^{\top} R u_k + x_N^{\top} Q_f x_N
\quad \text{(state cost, double integrator)}
$$

subject to:

- $x_{k+1} = F_\mathrm{RK4}(x_k, u_k)$ (RK4 discretisation, step $\delta t$)
- $u_\mathrm{lb} \leq u_k \leq u_\mathrm{ub}$ (control bounds)
- $\psi(h(x_k)) \geq 0,\; k = 0,\ldots,N$ (safe set at every node)
- $V(x_N) \geq 0$ (Reach-Avoid MPC) **or** $\phi(h(x_N)) \leq 0$ (Vanilla MPC) — terminal constraint, see table above
