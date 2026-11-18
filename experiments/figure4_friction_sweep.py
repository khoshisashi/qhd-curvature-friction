"""Fig. 4: P(F <= 0.01) against the friction at fixed objective-phase exposure, ten points on the
seed-11 landscape and the coarser grid on all five draws. Reads the rows that the sections table,
dense and draws of experiments.table6_friction_sweep store under results/, and writes
figures/fig2_friction.pdf, the file name the paper source uses."""
import json

import matplotlib.pyplot as plt
import numpy as np

from core.figures import COL, RESULTS, frame, save
from core.report import header
from experiments.table6_friction_sweep import COARSE


def stored(section):
    return json.loads((RESULTS / f"friction_sweep_{section}.json").read_text())


seed_11 = sorted((r for r in stored("table")["rows"] + stored("dense")["rows"] if r["rho"] <= 1),
                 key=lambda r: r["rho"])
coarse = {11: [r for r in seed_11 if r["rho"] in COARSE]}
coarse.update({int(seed): draw["rows"] for seed, draw in stored("draws")["draws"].items()})

fig, ax = plt.subplots(figsize=(COL, 1.54))
fig.subplots_adjust(left=0.16, right=0.975, bottom=0.215, top=0.955)
ax.semilogx([r["rho"] for r in seed_11], [r["p01"] for r in seed_11], "k-o", ms=1.8, lw=1.1,
            label="seed-11", zorder=4)
for k, seed in enumerate(sorted(coarse)):
    ax.semilogx([r["rho"] for r in coarse[seed]], [r["p01"] for r in coarse[seed]], "-", color="0.66", lw=1.1,
                marker="o", ms=1.8, mfc="0.66", mec="none", label="samples" if k == 1 else None, zorder=2)
ax.axvline(1.0, color="0.35", lw=0.6, ls=":", zorder=0)
ax.axvspan(0.032, 0.25, color="0.90", zorder=0)
ax.text(0.16, 0.10, "plateau", fontsize=8, color="0.10", ha="center")
ax.text(1.15, 0.10, "critical\ndamping", fontsize=8, color="0.10")
ax.set_xticks([0.002, 0.01, 0.1, 1, 4])
ax.set_xticklabels(["0.002", "0.01", "0.1", "1", "4"])
ax.set_xlim(0.0017, 4.8)
ax.set_ylim(0, 1.06)
ax.set_xlabel(r"$\rho$", labelpad=1.5)
ax.set_ylabel(r"$P(\tilde F\leq0.01)$", labelpad=1.5)
ax.legend(loc="lower left", frameon=False, handlelength=1.2, fontsize=6, borderpad=0.1, labelspacing=0.25,
          bbox_to_anchor=(0.02, 0.01))
frame(ax)

header("Fig. 4. Friction sweep at fixed objective-phase exposure")
save(fig, "fig2_friction")
best = max(seed_11, key=lambda r: r["p01"])
print(f"  seed 11, {len(seed_11)} frictions from {seed_11[0]['rho']:g} to {seed_11[-1]['rho']:g}: "
      f"P(F<=0.01) largest at rho = {best['rho']:g}")
print("  largest on the coarser grid, seeds " + ", ".join(str(seed) for seed in sorted(coarse)) + ": "
      + ", ".join(f"{max(coarse[seed], key=lambda r: r['p01'])['rho']:g}" for seed in sorted(coarse)))
plateau = [r["p01"] for r in seed_11 if 0.032 <= r["rho"] <= 0.25]
print(f"  on the plateau 0.032 <= rho <= 0.25, P(F<=0.01) varies by {np.ptp(plateau):.4f}")
