# Summary

Code and data of the article *Reach-Avoid Model Predictive Control with
Guaranteed Recursive Feasibility via Input Constrained Backstepping* by
Jianqiang Ding, Nishant Jayesh Bhave and Shankar A. Deka.

## Intent

The code implements the reach-avoid model predictive controller (RA-MPC) of
the article, whose terminal set is an input-constrained reach-avoid set built
by backstepping from an exponential control guidance-barrier function, and
reproduces every certificate, simulation, table, figure and video of the
article. The committed data allow every figure and number to be regenerated
without rerunning the experiments.

## Components

| Component | Location | Content |
| --- | --- | --- |
| Certificate synthesis (MATLAB) | `experiments/strict_xra/solve_vanilla_k1_*.m`, `experiments/strict_xra/design/`, `matlab/` | sum-of-squares programs (12) and (14) of the article for the virtual controller `k1` and the reach-avoid set, solved with SOSTOOLS and MOSEK; helpers between symbolic and polynomial objects |
| Controller construction (Python) | `experiments/strict_xra/design/cert.py`, `ra_mpc_cert.py`, `python/ra_mpc_dubins.py`, `controllers/` | ECGBF `Psi`, backstepping controller `u^c`, time scaling, Algorithm 1, vanilla MPC and the shared CasADi/IPOPT MPC pieces |
| Example 1, Dubins car | `experiments/strict_xra/final/` | three controllers on 300 common initial states, sampled-data checks of Lemma 2, success rate over the sampling period, timings, figures, `summary.txt`; notes in `PAPER_DUBINS.md` |
| Example 2, SCARA manipulator | `experiments/strict_xra/arm/`, `arm_top/` | workcell layout, eight kitting tasks under common torque bounds, the same experiments as Example 1 for each task; notes in `PAPER_ARM_TOP.md` |
| MuJoCo scene and videos | `experiments/strict_xra/scara/`, `docs/media/` | SCARA scene, Fig. 6 of the article, kitting-cycle and red-task videos; notes in `PAPER_SCARA.md` |
| Additional material not in the article | `experiments/strict_xra/franka/` | a Panda-like arm in MuJoCo with the same certificate construction, paused after the first task; notes in `PAPER_FRANKA.md` |
| Documentation | `README.md`, `experiments/strict_xra/README.md`, `REPORT.md`, `STUDY_14_DIRECT.md` | requirements, reproduction steps, folder index and the design decisions behind the experiments |

## Reproduction

See `README.md` for the software versions, the conda environments
(`environment.yml`, `environment-mujoco.yml`) and the commands that reproduce
each example in order, with run times on an Apple M5.

## Licence

The code and data are released under the MIT License (`LICENSE`). The
dependencies keep their own licences, listed in `README.md`.
