"""Fig. 1: the positive semidefinite lens, the design box of the entry parameterisation and the
image of a factor family, in the slice Sigma_22 = 0.05. Writes figures/fig3_feasible.pdf, the
file name the paper source uses, at the IEEEtran conference column width."""
import pathlib

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon, Rectangle
from scipy.spatial import ConvexHull

from core.report import header

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 7.0,
    "axes.labelsize": 7.0,
    "axes.linewidth": 0.6,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.0,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "pdf.fonttype": 42,
})

COL = 252.0 / 72.27   # IEEEtran conference column width, 252 pt, in inches
OUT = pathlib.Path(__file__).resolve().parents[1] / "figures" / "fig3_feasible.pdf"

rng = np.random.default_rng(7)
S22 = 0.05
BOX11, BOX12 = (0.02, 0.20), (-0.06, 0.10)

a = rng.uniform(*BOX11, 600000)
b = rng.uniform(*BOX12, 600000)
infeasible = 1 - (np.abs(b) <= np.sqrt(a * S22)).mean()

delta = np.array([0.018, 0.018])
V = np.array([[0.34, 0.10, 0.22],
              [0.09, 0.30, -0.19]])
X = rng.uniform(0, 1, size=(1500000, 3))
Sm = np.einsum("nm,im,jm->nij", X, V, V) + np.diag(delta)
feasible_factor = (np.abs(Sm[:, 0, 1]) <= np.sqrt(Sm[:, 0, 0] * Sm[:, 1, 1])).mean()
sel = np.abs(Sm[:, 1, 1] - S22) < 0.0035
pts = np.stack([Sm[sel, 0, 0], Sm[sel, 0, 1]], 1)
poly = pts[ConvexHull(pts).vertices]

fig, ax = plt.subplots(figsize=(COL, 2.30))
fig.subplots_adjust(left=0.155, right=0.985, bottom=0.185, top=0.965)

xg = np.linspace(0, 0.225, 500)
env = np.sqrt(xg * S22)

ax.add_patch(Rectangle((BOX11[0], BOX12[0]), BOX11[1] - BOX11[0], BOX12[1] - BOX12[0],
                       facecolor="none", edgecolor="0.45", hatch="////", lw=0.0, zorder=1))
ax.fill_between(xg, -env, env, color="0.90", lw=0, zorder=2)
ax.plot(xg, env, "k-", lw=0.9, zorder=4)
ax.plot(xg, -env, "k-", lw=0.9, zorder=4)
ax.add_patch(Rectangle((BOX11[0], BOX12[0]), BOX11[1] - BOX11[0], BOX12[1] - BOX12[0],
                       facecolor="none", edgecolor="k", ls="--", lw=1.0, zorder=5))
ax.add_patch(Polygon(poly, closed=True, facecolor="0.40", edgecolor="k", lw=0.6, zorder=6))

ax.annotate("hypercube register", (0.055, 0.093), xytext=(0.083, 0.108), fontsize=6.3,
            arrowprops=dict(arrowstyle="-", lw=0.5, color="k"))
ax.annotate("infeasible", (0.033, 0.074), xytext=(0.023, 0.118), fontsize=6.3,
            arrowprops=dict(arrowstyle="->", lw=0.5, color="k"))
ax.annotate("infeasible", (0.038, -0.049), xytext=(0.052, -0.079), fontsize=6.3,
            arrowprops=dict(arrowstyle="->", lw=0.5, color="k"))
ax.annotate("factor register", (0.115, 0.010), xytext=(0.140, -0.048), fontsize=6.3, color="k",
            arrowprops=dict(arrowstyle="->", lw=0.5, color="k"))
ax.annotate(r"$|\Sigma_{12}|=\sqrt{\Sigma_{11}\Sigma_{22}}$", (0.209, 0.1023), xytext=(0.150, 0.124),
            fontsize=6.3, ha="center", arrowprops=dict(arrowstyle="-", lw=0.5, color="k"))

ax.set_xlim(0, 0.225)
ax.set_ylim(-0.105, 0.140)
ax.set_xlabel(r"$\Sigma_{11}$", labelpad=1.0)
ax.set_ylabel(r"$\Sigma_{12}$", labelpad=1.0)
ax.set_xticks([0.05, 0.10, 0.15, 0.20])
ax.set_yticks([-0.10, -0.05, 0.0, 0.05, 0.10])
for side in ("top", "right"):
    ax.spines[side].set_visible(False)

OUT.parent.mkdir(exist_ok=True)
fig.savefig(OUT)
header("Fig. 1. Feasible slice")
print(f"wrote figures/{OUT.name}")
print(f"  infeasible share of the box in this slice: {infeasible:.3f}")
print(f"  feasible fraction of the factor image: {feasible_factor:.5f}")
print(f"  hull vertices of the factor polygon: {len(poly)}")
