"""Proposition 2, numerical check: the measurement distributions for F and F + b agree pointwise.
Seed-11 landscape, m = 5, critical damping. Exits with a non-zero status if they differ by more
than 1e-9."""
import sys

import numpy as np

from core.landscape import curvature, objective
from core.qhd import EPSILON, evolve, register_objective, rescale
from core.report import header

header("Gauge invariance of the measurement distribution")
F, arrays = objective(6, 2, 11)
mu = curvature(F)[0][0]
nu = 2 * np.sqrt(mu)
horizon = np.log(1 / EPSILON) / nu
f_tilde = rescale(register_objective(arrays, 5)[0])
p0 = evolve(f_tilde, nu, horizon, 3000)
print("largest pointwise difference between the distributions for F and F + b")
worst = 0.0
for b in (1.0, 1e3):
    gap = float(np.abs(evolve(f_tilde + b, nu, horizon, 3000) - p0).max())
    worst = max(worst, gap)
    print(f"  b = {b:<7g} {gap:.1e}")
print("printed in the paper: 1.3e-14 at b = 1 and 1.2e-11 at b = 1e3")
print("the differences are rounding errors of the phase, so the last digits depend on the platform")
sys.exit(1 if worst > 1e-9 else 0)
