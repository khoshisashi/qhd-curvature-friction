"""Sec. V-A, text: what the infeasible points of the d = 2 entry register do to the objective.

This illustration keeps the set-up of the original script it comes from, which differs from the
benchmark of Sec. IV: a grid of width 4 with 32 points per direction, a tent-shaped claim, and an
objective restricted to the central window |x_i| <= 0.5. The covariance entries are the three
design variables, on a register of m = 4 bits each.
"""
import itertools
import sys

import numpy as np

from core.report import Published, header

D, NX, WIDTH, RATE, HORIZON = 2, 32, 4.0, 0.05, 1.0
BOX = np.array([(0.02, 0.20), (0.02, 0.20), (-0.06, 0.10)])
SIGMA_TARGET = np.array([[0.10, 0.03], [0.03, 0.06]])

x1 = (np.arange(NX) - NX // 2) * WIDTH / NX
k1 = 2 * np.pi * np.fft.fftfreq(NX, d=WIDTH / NX)
K = np.stack(np.meshgrid(k1, k1, indexing="ij"))
X = np.stack(np.meshgrid(x1, x1, indexing="ij"))


def symbol(sigma):
    quad = 0.5 * sum(sigma[i, j] * K[i] * K[j] for i in range(D) for j in range(D))
    drift = sum((RATE - sigma[i, i] / 2) * K[i] for i in range(D))
    return RATE + quad - 1j * drift


def propagate(u, sigma):
    f = np.exp(-symbol(sigma) * HORIZON)
    mag = np.abs(f)
    f = np.where(mag > 1.0, f / mag, f)
    return np.fft.ifftn(f * np.fft.fftn(u))


def decode(z):
    s = BOX[:, 0] + z * (BOX[:, 1] - BOX[:, 0])
    return np.array([[s[0], s[2]], [s[2], s[1]]])


u0 = np.ones_like(X[0])
window = np.ones_like(X[0], dtype=bool)
for i in range(D):
    u0 = u0 * np.clip(1.0 - np.abs(X[i]) / 0.8, 0.0, None)
    window &= np.abs(X[i]) <= WIDTH / 8
u0 /= np.linalg.norm(u0)
u_target = propagate(u0, SIGMA_TARGET)

header("Sec. V-A. Infeasible points of the d = 2 entry register")
n = 2 ** 4
grid = np.linspace(0.0, 1.0, n)
points = [np.array([grid[a], grid[b], grid[c]]) for a, b, c in itertools.product(range(n), repeat=3)]
obj = np.array([np.linalg.norm((propagate(u0, decode(z)) - u_target)[window]) ** 2 for z in points])
feasible = np.array([np.linalg.eigvalsh(decode(z)).min() >= 0 for z in points])
growth = np.array([np.abs(np.exp(-symbol(decode(z)) * HORIZON)).max() for z in points])
rescaled = (obj - obj.min()) / (obj.max() - obj.min())

n_bad = int((~feasible).sum())
print(f"register points                         {len(points)}")
print(f"infeasible points                       {n_bad}")
print(f"  with max_k |exp(-a(k) T)| > 1         {int((growth[~feasible] > 1).sum())}")
print(f"feasible points with that maximum > 1   {int((growth[feasible] > 1).sum())}")
print(f"rescaled objective, feasible minimum    {rescaled[feasible].min():.4f}")
print(f"rescaled objective on infeasible points min {rescaled[~feasible].min():.4f}  "
      f"median {np.median(rescaled[~feasible]):.4f}  max {rescaled[~feasible].max():.4f}")
pub = Published()
pub.add("smallest rescaled objective on the infeasible points", rescaled[~feasible].min(), 0.091, 3)
sys.exit(1 if pub.report() else 0)
