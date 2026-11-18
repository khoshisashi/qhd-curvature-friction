"""Fig. 3: raw residual of the four families on the thirty targets, sorted by the as-drawn value.
Reads results/representation_gap.json, which experiments.table2_representation_gap writes, and
writes figures/fig4_gap.pdf, the file name the paper source uses."""
import json

import matplotlib.pyplot as plt
import numpy as np

from core.figures import COL, RESULTS, frame, save
from core.report import header

families = json.loads((RESULTS / "representation_gap.json").read_text())["families"]
drawn, orth, five, vardiag = (np.array(families[name]["raw_gap"])
                              for name in ("M3_drawn", "M3_orth", "M5_drawn", "M3_vardiag"))
order = np.argsort(drawn)[::-1]
index = np.arange(1, len(drawn) + 1)

fig, ax = plt.subplots(figsize=(COL, 2.0))
fig.subplots_adjust(left=0.145, right=0.975, bottom=0.24, top=0.96)
ax.plot(index, drawn[order], "o", ms=2.9, color="k", label="$M=3$, as drawn", zorder=4)
ax.plot(index, orth[order], "o", ms=2.9, mfc="none", mec="k", mew=0.7, label="$M=3$, orthogonal", zorder=3)
ax.plot(index, vardiag[order], "^", ms=3.0, mfc="none", mec="0.35", mew=0.7,
        label="$M=3$, variable diag.", zorder=2)
ax.plot(index, five[order], "s", ms=2.6, mfc="none", mec="0.6", mew=0.7, label="$M=5$, as drawn", zorder=1)
ax.set_xlim(0.2, len(drawn) + 0.8)
ax.set_ylim(0, 0.066)
ax.set_xticks([1, 10, 20, 30])
ax.set_xlabel("target, sorted by the as-drawn residual", labelpad=1.5)
ax.set_ylabel(r"raw residual $\Delta F$", labelpad=1.5)
ax.legend(loc="upper right", frameon=False, handlelength=1.0, fontsize=5.7, borderpad=0.1,
          labelspacing=0.22, handletextpad=0.4)
frame(ax)

header("Fig. 3. Raw residual on the thirty targets")
save(fig, "fig4_gap")
print(f"  as drawn above orthogonal in {int((drawn > orth).sum())} of {len(drawn)}, "
      f"above variable diagonal in {int((drawn > vardiag).sum())}, "
      f"above five factors in {int((drawn > five).sum())}")
smallest = min(x.min() for x in (drawn, orth, five, vardiag))
print(f"  largest raw residual {drawn.max():.4f}, smallest {smallest:.4f}")
