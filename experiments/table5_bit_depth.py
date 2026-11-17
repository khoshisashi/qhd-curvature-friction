"""Table V: QHD on the seed-11 landscape against the design bit depth, at critical damping under the
terminal rule lambda(S) = 1 / epsilon. The run is repeated with 10^4 steps to show convergence."""
import sys

import numpy as np

from core.landscape import XI_STAR, curvature, objective
from core.qhd import EPSILON, evolve, register_objective, rescale, summary
from core.report import Published, header

PRINTED = {4: (0.9624, 0.0555, 0.0667), 5: (0.9876, 0.0149, 0.0323), 6: (0.9936, 0.0207, 0.0159)}

header("Table V. QHD against design bit depth, seed-11 landscape")
F, arrays = objective(6, 2, 11)
mu = curvature(F)[0][0]
nu = 2 * np.sqrt(mu)
horizon = np.log(1 / EPSILON) / nu
print(f"mu = {mu:.4f}, nu = {nu:.4f}, S = {horizon:.3f}")
print(f"{'m':>2} {'points':>8} {'P(F<=0.05)':>11} {'|mode-xi*|':>11} {'spacing':>8} {'volume of F<=0.05':>18} {'max change at 1e4 steps':>24}")
pub = Published()
for m in (4, 5, 6):
    f_grid, grid = register_objective(arrays, m)
    f_tilde = rescale(f_grid)
    a = summary(evolve(f_tilde, nu, horizon, 3000), f_tilde, grid, XI_STAR)
    b = summary(evolve(f_tilde, nu, horizon, 10000), f_tilde, grid, XI_STAR)
    change = max(abs(a[k] - b[k]) for k in ("p05", "dist"))
    spacing = 1.0 / (2 ** m - 1)
    print(f"{m:>2} {f_tilde.size:>8} {a['p05']:>11.4f} {a['dist']:>11.4f} {spacing:>8.4f} "
          f"{(f_tilde <= 0.05).mean():>18.3f} {change:>24.1e}")
    p05, dist, sp = PRINTED[m]
    pub.add(f"m = {m} P(F <= 0.05)", a["p05"], p05, 4)
    pub.add(f"m = {m} distance", a["dist"], dist, 4)
    pub.add(f"m = {m} spacing", spacing, sp, 4)
sys.exit(1 if pub.report() else 0)
