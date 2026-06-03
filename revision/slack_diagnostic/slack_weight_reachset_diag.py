"""Phase 1.2 diagnostic: certified reach-avoid funnel size vs slack_weight (dubins).

The per-sample slack SOP returns (delta, lambda). For the reduced single-integrator
system ydot = k1, the certified reach-avoid funnel is

    F(delta/lambda) = { y in X_S\X_T : psi(y) := safe(y) >= delta/lambda },

because along the flow  d/dt safe = grad(safe).k1 >= lambda*safe - delta, which is
positive (barrier grows -> trajectory pushed into the target while staying safe)
ONLY where safe > delta/lambda. A large delta/lambda shrinks -- and once it exceeds
max_{X_S\X_T} safe, EMPTIES -- that funnel.

This script reads slack_weight_sweep_dubins.csv (written by sweep_slack_weight_dubins.m),
densely samples X_S\X_T in output space, and reports, for each slack_weight (and the
vanilla / no-per-sample-constraints reference), the fraction of dense points lying in
F(delta/lambda). Run sweep_slack_weight_dubins.m (MATLAB) first to produce the CSV.

Self-contained: defines the dubins safe/target sets locally (mirrors synth_lib); reads
only revision/data/; writes a plot to revision/data/.
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sympy as sp

REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # revision/
DATA = os.path.join(REV, "data")
CSV = os.path.join(DATA, "slack_weight_sweep_dubins.csv")

# ---- dubins safe/target sets in output coords (mirror synth_lib; self-contained) ----
y1, y2 = sp.symbols("y1 y2")
h_raw = -(y1**4 + y2**4 - 16) * (y1**4 + y2**4 - 4)
target_y = (y2) ** 2 + ((y1 + 1.7) / 0.5) ** 2 - 0.4
safe_y = 1e-3 * (-target_y + 300) * h_raw
safe_f = sp.lambdify((y1, y2), safe_y, "numpy")
targ_f = sp.lambdify((y1, y2), target_y, "numpy")

# ---- dense X_S\X_T sample (safe >= 0 AND outside target) ----------------------------
N = 4_000_000
LO, HI = -2.2, 2.2  # output-space window covering the safe annulus (|y|^4 <= 16 -> |y|<=2)
rng = np.random.default_rng(0)
Y = rng.uniform(LO, HI, size=(2, N))
sv = np.asarray(safe_f(Y[0], Y[1]), float)
tv = np.asarray(targ_f(Y[0], Y[1]), float)
in_set = (sv >= 0) & (tv > 0)  # X_S \ X_T
safe_in = sv[in_set]
safe_max = float(safe_in.max())
print(
    f"dense X_S\\X_T: {safe_in.size:,} of {N:,} samples in set;  "
    f"safe(y) range = [{safe_in.min():.4g}, {safe_max:.4g}]"
)
print(f"=> funnel is EMPTY whenever delta/lambda > max safe = {safe_max:.4g}\n")

# ---- read the sweep + report the funnel fraction per slack_weight -------------------
df = pd.read_csv(CSV)
print(
    f"{'slack_weight':>14} {'delta':>12} {'lambda':>12} {'delta/lambda':>14} "
    f"{'frac in funnel':>16} {'verdict':>10}"
)
print("-" * 84)
fracs, labels, taus = [], [], []
for _, r in df.iterrows():
    tau = float(r["delta_over_lambda"])
    frac = float(np.mean(safe_in >= tau))  # fraction of X_S\X_T inside {safe >= tau}
    wlabel = "vanilla" if np.isnan(r["slack_weight"]) else f"{r['slack_weight']:g}"
    verdict = "EMPTY" if frac == 0.0 else ("FULL" if frac > 0.99 else "partial")
    print(
        f"{wlabel:>14} {r['delta']:>12.4g} {r['lambda']:>12.4g} {tau:>14.4g} "
        f"{100 * frac:>15.2f}% {verdict:>10}"
    )
    fracs.append(100 * frac)
    labels.append(wlabel)
    taus.append(tau)

# ---- plot: funnel fraction vs slack_weight (vanilla shown as a reference line) -------
wmask = ~df["slack_weight"].isna().values
ws = df["slack_weight"].values[wmask]
fr = np.array(fracs)[wmask]
van_frac = np.array(fracs)[~wmask]
fig, ax = plt.subplots(figsize=(7.2, 5.0))
ax.semilogx(ws, fr, "o-", color="#2166ac", lw=1.8, label="per-sample slack SOP")
if van_frac.size:
    ax.axhline(
        van_frac[0], color="#1a9850", ls="--", lw=1.5,
        label=f"vanilla (no per-sample constraints) = {van_frac[0]:.1f}%",
    )
ax.set_xlabel("slack_weight")
ax.set_ylabel(r"% of $X_S\setminus X_T$ inside funnel $\{$safe $\geq \delta/\lambda\}$")
ax.set_ylim(-3, 103)
ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=9)
ax.set_title("dubins: certified reach-avoid funnel size vs slack_weight")
fig.tight_layout()
out_png = os.path.join(DATA, "slack_weight_reachset_dubins.png")
fig.savefig(out_png, dpi=130, bbox_inches="tight")
print(f"\nwrote {out_png}")
