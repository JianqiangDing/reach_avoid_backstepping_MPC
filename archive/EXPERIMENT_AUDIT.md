# L-CSS 修订实验计划 — 论文 26-0694

> 对应 paper 仓库 `LCSS_Revsion/REVISION_PLAN.md` §4 / §7。
> 仓库 `reach_avoid_backstepping_MPC`，分支 `lcss-revision`。
> 论文例：Ex1 = Dubins car；Ex2 = two-link manipulator。
> 修订：2026-05-31。

---

## 0. 工作约定

- 所有修订实验代码与产出落顶层新目录 `revision/`，且**完全自包含**：运行时**不读主仓 `data/`、不导入/调用主仓 `python/` / `matlab/` / `controllers/`**。理由——若 revision 依赖主仓活文件，主仓一变就会悄悄改变 revision 结果，破坏后续所有实验的可靠性。需要的逻辑在 `revision/` 内重新实现：A(x) 在 notebook 内符号推导 + 自采样 + 对符号 A 校验；Phase 1.2/1.3 的 SOP 求解器与复现所需的论文控制器以**冻结副本**形式放在 `revision/`。所有产出写 `revision/data/`，不修改主仓。
- **测试与可视化一律用 Jupyter notebook**（`.ipynb`），不写裸 `.py` 实验脚本；共享逻辑（无 UI 的纯函数）写成 `.py` 模块供 notebook `import`。
- **按 phase 串行推进**：当前 phase 不全部跑通并核对，不进入下一 phase。本文档只详写当前 phase 1，后续 phase 给 stub，phase 1 收尾时再回填细节。

```
revision/
├── README.md
├── decoupling_check/               # Phase 1.1
│   ├── lib.py                      #   det A(x) 与采样的纯函数
│   ├── dubins_detA.ipynb
│   └── manipulator_detA.ipynb
├── slack_diagnostic/               # Phase 1.2  —— 仅诊断，不产出最终控制器
│   ├── solvesop_bounded_control_slack.m   #  带 per-sample slack 的 SOP
│   ├── solve_dubins_slack.m
│   ├── solve_manipulator_slack.m
│   ├── dubins_slack_diag.ipynb     #   读 .m 输出，画 slack 分布，选 ub_eff + 过滤采样
│   └── manipulator_slack_diag.ipynb
├── controller_synthesis/           # Phase 1.3  —— 最终无 slack 控制器（实验真正用）
│   ├── solvesop_bounded_control_clean.m   #  与主仓等价的硬约束 SOP
│   ├── solve_dubins_clean.m
│   ├── solve_manipulator_clean.m
│   ├── dubins_synth.ipynb          #   读 .m 输出 + Phase 1.2 选定的 ub_eff & 过滤集
│   └── manipulator_synth.ipynb
├── controllers/                    # Phase 1.3 导出的 revision 控制器 .py（实验输入）
├── data/                           # 所有 revision 产出
├── success_module/                 # Phase 2  ——细节待 Phase 1 收尾后补
└── experiments/                    # Phase 3  ——细节待 Phase 2 收尾后补
```

---

## 1. Phase 1.1 — Decoupling matrix A(x) 可逆性核验

### 1.1 目的

确认 backstepping FL 坐标在合成区间 $\mathcal X_S$ 上一致有效。$A(\bm x)$ 出现奇异/接近奇异的位置将导致后续 SOP 输出的 $u^*(\bm x)$ 在闭环触及该位置时爆炸（Dubins 含 $1/v$，manipulator 含 $1/\sin q_2$）。Phase 1.1 的结论决定 Phase 1.2 SOP 用哪个状态盒。

### 1.2 自包含实现（in-notebook 解析 A + 自采样 + 符号校验）

revision/decoupling_check/ 的两个 notebook **完全自包含**，不读主仓 CSV、不导入主仓代码：

- 在 notebook 内用 SymPy 推导 $A(\bm x)$ 与 $\det A(\bm x)$；对 manipulator 必须代入流形 $x_5 = x_1 + x_2$ 后 $\det A$ 才化简为 $\propto \sin x_2$；
- 提供向量化 Python 评估器 `A_metrics`（2×2 闭式 SVD 出 $|\det A|$ / $\sigma_{\min}$ / cond），并**对 in-notebook 符号 A 校验**（NumPy SVD `A_sym_fun`，误差 ~机器精度）——正确性不依赖任何外部文件；
- 自己用 NumPy 采样（固定 seed），safe/target 隶属由 notebook 内的多项式直接评估；
- 产出全部写 `revision/data/`。

> 早期"消费主仓 `data/decoupling_*.csv`"的做法已废弃（见 §0 自包含约定）。

### 1.3 Notebook 测试 `dubins_detA.ipynb` / `manipulator_detA.ipynb`

每个 notebook 的 cell 顺序：

1. §0 参数表（系统、输出、safe/target、采样盒、$\epsilon$、$N$、seed）；
2. §1 系统 + $A$ 解析推导（markdown）+ SymPy cell（符号 A，并 `lambdify` 出 `A_sym_fun`）；
3. §2 setup：定位 repo root（仅为定位 `revision/data/` 输出）、定义 `A_metrics`、**对符号 A 校验**、自采样；
4. §3 $\mathcal X_S\setminus\mathcal X_T$ 上 `|det A|` / `sigma_min` / `cond` 统计 + 全点 $\sigma_{\min}\ge\epsilon$ 确认；
5. §4 输出空间覆盖图：$\mathcal X_S\setminus\mathcal X_T$ 采样点按 $|\det A|$ 着色 + safe/target 基准线（存 `revision/data/detA_<example>_output_space.png`）；
6. §5 判据 cell：打印 `VERDICT: PASS | NEEDS_SHRINK`；
7. §6 出口 cell：写 `revision/data/phase1_1_outputs_<example>.json`，含 `verdict`, `self_contained`, `validated_against`, `stats`, `X_S_eff_def`。

manipulator 额外报告放松盒 $\mathcal X_S^{\rm eff}$（$x_2$ 边距 0.3、$x_1\in[-\pi,\pi]$，覆盖 ~99% 可达 safe 集）。

### 1.4 判据

阈值 $\epsilon = 10^{-2}$（依例可在 notebook 里调整）。

- $\min_{\bm x\in\mathcal X_S} |\det A| \ge \epsilon$ 且 $\min \sigma_{\min} \ge \epsilon$ → **PASS**，Phase 1.2 用当前 $\mathcal X_S$。
- 否则 → 在 notebook 里定位 arg min 区域，给出 $\mathcal X_S^{\text{eff}} \subset \mathcal X_S$ 的显式约束（如 Dubins 加 $v \ge v_{\min}$），并把 $v_{\min}$ 选成使 $\min|\det A| = \epsilon$ 的阈值。这个 $\mathcal X_S^{\text{eff}}$ 写回 notebook 的最后一个 cell 作为后续输入。

### 1.5 Phase 1.1 出口

每例 notebook 最后输出一个明确的"$\mathcal X_S^{\text{eff}}$ 定义"（如果 PASS，则 $\mathcal X_S^{\text{eff}} = \mathcal X_S$）。这是 Phase 1.2 的唯一输入。

---

## 2. Phase 1.2 — Slack 诊断（不产出最终控制器）

### 2.1 目的

slack 不是控制器设计变量，**仅用作诊断工具**，回答两件事：

1. **采样点过滤**：在 $\mathcal X_S^{\text{eff}}$ 上的候选采样集 $\mathcal S$ 里，哪些点对 SOP 求解"难"（强行满足名义 bound 时局部需要的松弛大）→ 这些点保留还是剔除；
2. **bound 选择**：在 $\mathcal S$ 上 SOP 能保证的最紧 effective bound $u_{\max}^{\text{eff}}$ 是多少。

诊断的输出（过滤后采样集 $\mathcal S^{\text{filtered}}$ + 选定 $u_{\max}^{\text{eff}}$）作为 Phase 1.3 输入；**Phase 1.2 不产出任何被实验调用的控制器**。

### 2.2 诊断 SOP 形式（MATLAB）

`revision/slack_diagnostic/solvesop_bounded_control_slack.m`：主仓 `matlab/solvesop_bounded_control.m` 的**冻结副本**（拷进 `revision/`，运行时不引用主仓原件），控制约束改为**每采样点一个 slack $s_j \ge 0$**：

$$
-u_{\max} - s_j \le u^*(\bm x_j) \le u_{\max} + s_j, \quad s_j \ge 0, \quad j=1,\dots,|\mathcal S|.
$$

目标：$\min \sum_j s_j$（或可选 $\min \max_j s_j$，作为 sweep 选项）。返回值除原 `u_opt, V_opt, k1_opt` 外，追加：

- `s_vec`：长度 $|\mathcal S|$ 的 slack 数组；
- `solver_status`, `wallclock_s`。

`solve_dubins_slack.m` / `solve_manipulator_slack.m`：调用上面，输入 $\mathcal X_S^{\text{eff}}$（来自 Phase 1.1）+ 名义 $u_{\max}$ + 初始采样集 $\mathcal S_0$，输出 `revision/data/slack_diag_<example>.mat`（含 `S0_points`, `s_vec`, `u_opt_struct`, `V_opt_struct`, `solver_status`）。

> **注**：`u_opt`/`V_opt` 在这里只用于诊断对比（看看带 slack 后形态合理性），**不导出 .py、不进入实验**。

### 2.3 Diag Notebook `dubins_slack_diag.ipynb` / `manipulator_slack_diag.ipynb`

cell 顺序：

1. bootstrap；
2. 读 `revision/data/slack_diag_<example>.mat`；
3. **slack 分布分析**：
   - histogram of $s_j$；
   - $s_j$ 按从大到小排序的 elbow 图；
   - $s_j$ vs $\|\bm x_j - \bm x_{\text{center}}\|$ 散点（看是否难点都聚在边界）；
   - $s_j$ vs $|\det A(\bm x_j)|$（用 Phase 1.1 数据）散点（看是否难点都在接近奇异）。
4. **采样点过滤决策**：
   - 设阈值 $s_{\text{drop}}$（默认 `percentile(s_vec, 95)`，notebook 里可调），剔除 $s_j > s_{\text{drop}}$ 的点 → 得 $\mathcal S^{\text{filtered}}$；
   - 在 notebook 里画"被剔除点"在状态空间的散布，确认其位置合理（应聚在 $A$ 接近奇异 / 边界附近）。
5. **bound 选择**：
   - 以 $\mathcal S^{\text{filtered}}$ 上的 $\max_j s_j$（设为 $s^*$）定 $u_{\max}^{\text{eff}} = u_{\max} - s^*$（即：把名义 bound 收紧 $s^*$，使得在过滤后的集合上 SOP 不再需要 slack）；
   - 检查 $s^* / u_{\max}$：≤5% / 5–30% / >30% 三档分别提示"几乎兑现名义 / 弱化但可接受 / 需重新过滤或回 Phase 1.1"。
6. **导出 Phase 1.3 输入**：写 `revision/data/phase1_2_outputs_<example>.json`：
   ```json
   {
     "X_S_eff_def": "...",
     "S_filtered_npz": "revision/data/S_filtered_<example>.npz",
     "u_max_nom": ...,
     "s_drop_threshold": ...,
     "s_star": ...,
     "u_max_eff": ...
   }
   ```
   过滤后采样点写 `revision/data/S_filtered_<example>.npz`。

### 2.4 判据 / 出口

- 出口物：$\mathcal S^{\text{filtered}}$ + $u_{\max}^{\text{eff}}$；
- 通过条件：过滤后被剔除点 ≤ 10% 且 $u_{\max}^{\text{eff}} / u_{\max} \ge 0.7$；否则 notebook 里给"回 Phase 1.1 收缩 $\mathcal X_S^{\text{eff}}$"的明确提示，不放行 Phase 1.3。

### 2.5 Phase 1.2 绝不做的事

- 不导出 slack 版控制器 `.py` 进 `revision/controllers/`；
- 不让任何 Phase 3 实验直接调用 slack 版控制器；
- 不把 slack 当成"反正小所以可以忍"的常驻松弛。

---

## 3. Phase 1.3 — 最终 k1 控制器合成（无 slack 硬约束 SOP）

### 3.1 目的

用 Phase 1.2 给出的 ($\mathcal S^{\text{filtered}}$, $u_{\max}^{\text{eff}}$)，跑**硬约束** SOP——与主仓 `matlab/solvesop_bounded_control.m` 完全等价的形式，没有 slack——合成实验真正用的 $u^*, V, k_1$，导出到 `revision/controllers/`。

### 3.2 MATLAB 入口

`revision/controller_synthesis/solvesop_bounded_control_clean.m`：主仓 `matlab/solvesop_bounded_control.m` 的**冻结副本**（函数级等价，运行时不引用主仓原件），仅签名上接受 $\mathcal S^{\text{filtered}}$ 与 $u_{\max}^{\text{eff}}$ 而不是默认 $\mathcal S$ 与 $u_{\max}$。控制约束：

$$
-u_{\max}^{\text{eff}} \le u^*(\bm x_j) \le u_{\max}^{\text{eff}}, \quad \forall \bm x_j \in \mathcal S^{\text{filtered}}.
$$

`solve_dubins_clean.m` / `solve_manipulator_clean.m`：调用上面 + 读 Phase 1.2 的 `phase1_2_outputs_<example>.json` 拿参数。

### 3.3 导出

`solve_*_clean.m` 在求解后调用 `revision/controller_synthesis/export_to_python.m`（主仓 `matlab/export_to_python.m` 的冻结副本，导出目标改为 `revision/controllers/`），把 `u_opt, certificate_opt, k1_opt` 写到 `revision/controllers/<example>_clean_<timestamp>.py`，同名 `.json` 登记：

```json
{
  "u_max_nom": ...,
  "u_max_eff": ...,
  "S_filtered_size": ...,
  "X_S_eff_def": "...",
  "solver_status": "...",
  "wallclock_s": ...
}
```

### 3.4 Notebook 验证 `dubins_synth.ipynb` / `manipulator_synth.ipynb`

cell 顺序：

1. bootstrap；
2. `import` 新导出的 `revision/controllers/<example>_clean_<timestamp>` 拿 `u_opt`；
3. **bound 兑现核验**：在 $\mathcal X_S^{\text{eff}}$ 上 dense 采样 $N=10^5$（独立 seed，不与 $\mathcal S^{\text{filtered}}$ 重合），评估 $|u^*(\bm x)|$，输出 `max`、`p99`、`p99.9`；预期 `max ≤ u_max_eff`（允许 1e-6 数值容差），否则报警；
4. **certificate $V$ 形状对比**：与主仓对应 controller 的 $V$ 等高线叠加，确认 RA 集面积可比；
5. **闭式律 sanity rollout**：从 $\mathcal X_S^{\text{eff}}\setminus\mathcal X_T$ 内选 5 个代表性 $\bm x_0$，用 $u_{\text{opt}}$ 作 closed-form 控制律跑 DOP853 闭环，记录 $\max_t |u(t)|$、是否进入 $\mathcal X_T$、是否离开 $\mathcal X_S^{\text{eff}}$；
6. 通过/失败标记写 `revision/data/phase1_3_outputs_<example>.json`，含最终选用的控制器路径，作为 Phase 2 的唯一控制器输入。

### 3.5 判据 / 出口

- bound 兑现：`max|u*| ≤ u_max_eff + 1e-6` → PASS；
- 5 条 sanity rollout：≥4 条成功（短 horizon 内进入目标且不离开 $\mathcal X_S^{\text{eff}}$）→ PASS；
- 否则诊断后回 Phase 1.2 调过滤阈值或回 Phase 1.1。

### 3.6 Phase 1 整体出口

`revision/data/phase1_outputs_<example>.json`（由 Phase 1.3 notebook 写出）含三件物：

1. `controller_py_path`：最终 `revision/controllers/<example>_clean_<timestamp>.py`
2. `u_max_eff`：最终采用的输入 bound
3. `X_S_eff_def`：最终采用的合成区间

**这三件物是 Phase 2 / 3 的唯一接口**；后续实验只读它们，不回看 Phase 1.2 的诊断数据。

---

## 4. Phase 2 — 统一 success 模块（stub，待 Phase 1 完成后回填）

### 4.1 范围预告（不在此 phase 实现）

- `revision/success_module/success.py` 提供闭环 rollout + recursive-feasibility 逐步核验，返回 `(feas_0, feas_all, safe_all, reached)` 四元；
- 三分类：`outside_RoA` / `success` / `theorem_violation`；
- 采样集生成：`sample_ex1`（≥300 candidate）、`sample_ex2`（≥500 valid candidate）；
- 控制器适配类：`ClosedFormController` / `MPCController` / `VanillaMPCController`，输入参数 = Phase 1 出口 JSON。

### 4.2 何时回填本节

Phase 1.3 两个例子的 `phase1_outputs_<example>.json` 都生成且 PASS 之后，开 phase 2 的细节会写，包括接口签名、采样集 fixture、控制器适配类、单测（notebook 形式）。

---

## 5. Phase 3 — 实验 notebooks（stub，待 Phase 2 完成后回填）

### 5.1 范围预告

- §4.2 复现 Ex1（`revision/experiments/reproduce_ex1.ipynb`）
- §4.3 加密 Ex2（`revision/experiments/reproduce_ex2.ipynb`）
- §4.4 计时（`revision/experiments/timing_benchmark.ipynb` + MATLAB 离线 tic/toc 写 csv）
- §4.5 采样周期扫描（`revision/experiments/sampling_period_sweep.ipynb`）
- §4.6 扰动扫描（`revision/experiments/disturbance_sweep.ipynb`）

### 5.2 何时回填本节

Phase 2 的 `success_module` notebook 单测全绿后回填。

---

## 6. 当前活跃任务（Phase 1.1）

唯一活跃工作：

1. 建 `revision/` 骨架 + `README.md`；
2. 实现 `revision/decoupling_check/lib.py`；
3. 跑 `dubins_detA.ipynb` 与 `manipulator_detA.ipynb` 到判据 PASS 或给出 $\mathcal X_S^{\text{eff}}$ 收缩方案；
4. 两个 notebook 输出 cell 里固化 $\mathcal X_S^{\text{eff}}$ 定义，作为 Phase 1.2 的唯一输入。

Phase 1.2 的 MATLAB / notebook 在 Phase 1.1 两个 notebook 都 PASS 之后再开始。

---

## 7. 不在本轮的工作

- §4.7 参数扰动 sweep（仅记录方法，不实现）；
- §4.8 仓库公开 / requirements.txt / run_all / LICENSE / CITATION；
- 回填主仓 notebooks 调用 revision 模块（保留主仓为论文版可复现基线）；
- double-integrator no-MPC notebook 旧导入修复（不在论文）。

---

## 附录 A — 主仓代码现状速查（仅作参考，revision 运行时不引用）

| 例  | 控制器           | `notebooks/`                                    | 主仓 `controllers/*.py`                                        | `data/*.npz`                    |
| --- | ---------------- | ----------------------------------------------- | -------------------------------------------------------------- | ------------------------------- |
| Ex1 | Unconstrained RA | `example_dubins_car_unconstrained_reach_avoid`  | `sop_bounded_control_unconstrained_controller_20260316_211252` | `traj_controls_reach_avoid_mpc` |
| Ex1 | RA-MPC           | `example_dubins_car_reach_avoid_mpc`            | `sop_bounded_control_dubins_car_result_20260316_211504`        | `traj_controls_reach_avoid_mpc` |
| Ex1 | Vanilla MPC      | `example_dubins_car_vanilla_mpc`                | （无，目标集终端约束内置）                                     | `traj_controls_vanilla_mpc`     |
| Ex2 | Unconstrained RA | `example_manipulator_unconstrained_reach_avoid` | `k1_acrobot_cdc2026`                                           | `traj_controls_ra_mpc_acrobot`  |
| Ex2 | RA-MPC           | `example_manipulator_reach_avoid_mpc`           | `sop_bounded_control_acrobot_result_20260317_222858`           | `traj_controls_ra_mpc_acrobot`  |
| Ex2 | Vanilla MPC      | `example_manipulator_vanilla_mpc`               | （无）                                                         | —                               |

本表仅供查阅论文版的对应关系。需要复现这些控制器时，把对应 `.py` 以**冻结副本**拷入 `revision/controllers/` 再用，**不在运行时引用主仓原件**（见 §0 自包含约定）。
