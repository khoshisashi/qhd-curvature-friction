"""Secs. IV and VI, text: numerical checks behind the coarse grids and the curvature measurement.

Spatial refinement at d = 4 and d = 6, the dependence of the curvature on d, the truncation to the
heaviest modes, and the sensitivity of the curvature to the scan and to the finite-difference step.
All landscapes use seed 11. The refined grids with more than 20000 modes keep the 20000 heaviest;
the comparison across dimensions uses all modes.
"""
import sys

import numpy as np

from core.landscape import N_KEEP, curvature, factor_model, objective
from core.report import Published, header

header("Refinement, truncation and dimension checks")
pub = Published()

print("curvature against dimension, N_x = 4, all modes")
print(f"{'d':>3} {'modes':>10} {'mu':>8} {'kappa':>8}")
for d, printed in ((4, 0.015), (6, 0.156), (8, 0.095), (10, 0.076)):
    ev = curvature(objective(d, 2, 11, n_keep=4 ** d)[0])[0]
    print(f"{d:>3} {4 ** d:>10,} {ev[0]:>8.4f} {ev[-1] / ev[0]:>8.2f}", flush=True)
    pub.add(f"mu at d = {d}", ev[0], printed, 3)

print()
print("spatial refinement at d = 4: Hessian spectrum")
for n_sp in (2, 3, 4, 5):
    ev = curvature(objective(4, n_sp, 11)[0])[0]
    print(f"  N_x = {2 ** n_sp:>2}: " + "  ".join(f"{v:.4f}" for v in ev), flush=True)
    pub.add(f"d = 4, N_x = {2 ** n_sp}, smallest eigenvalue", ev[0], 0.0150, 4)
    if n_sp >= 3:
        pub.add(f"d = 4, N_x = {2 ** n_sp}, second eigenvalue", ev[1], 0.494, 3)
        pub.add(f"d = 4, N_x = {2 ** n_sp}, largest eigenvalue", ev[2], 1.486, 3)
print("  the two larger eigenvalues are compared with (0.494, 1.486) from N_x = 8 on;")
print("  at N_x = 4 they are 2 to 3% lower")

print()
print("two resolutions at d = 6")
for n_sp, printed in ((2, 0.1560), (3, 0.1557)):
    ev = curvature(objective(6, n_sp, 11)[0])[0]
    print(f"  N_x = {2 ** n_sp}: mu = {ev[0]:.4f}", flush=True)
    pub.add(f"mu at d = 6, N_x = {2 ** n_sp}", ev[0], printed, 4)

print()
print(f"truncation to the {N_KEEP} heaviest modes at d = 6, N_x = 8")
weight = factor_model(6, 3, 11)[0]
share = np.sort(weight)[::-1][:N_KEEP].sum() / weight.sum()
print(f"  share of the energy kept: {100 * share:.2f}%")
pub.add("energy kept, per cent", 100 * share, 99.88, 2)
F_cut = objective(6, 3, 11)[0]
F_all = objective(6, 3, 11, n_keep=weight.size)[0]
g = np.linspace(0.0, 1.0, 13)
pts = [np.array([a, b, c]) for a in g for b in g for c in g]
cut, full = np.array([F_cut(z) for z in pts]), np.array([F_all(z) for z in pts])
shift = np.abs(cut / (cut.max() - cut.min()) - full / (full.max() - full.min())).max()
print(f"  largest change of the rescaled landscape on the 13-point scan: {shift:.1e}")
print(f"  minimiser with all modes: {np.round(curvature(F_all)[1], 4)}, truncated: {np.round(curvature(F_cut)[1], 4)}")

print()
print("sensitivity of the curvature at d = 6, N_x = 4")
F = objective(6, 2, 11)[0]
ev, z, f_lo, f_hi = curvature(F)
for step in (2e-2, 5e-3, 1.25e-3):
    print(f"  finite-difference step {step:g}: " + "  ".join(f"{v:.5f}" for v in curvature(F, step=step)[0]))
for n_scan in (13, 17, 25):
    print(f"  scan with {n_scan} points per axis: F_max = {curvature(F, n_scan=n_scan)[3]:.9e}")
print(f"  mu with the scan minimum {ev[0]:.5f}, with F_min = 0 {ev[0] * (f_hi - f_lo) / f_hi:.5f}, "
      f"relative change {100 * f_lo / f_hi:.3f}%")
sys.exit(1 if pub.report() else 0)
