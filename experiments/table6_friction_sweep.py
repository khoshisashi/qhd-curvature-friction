"""Table VI and the statements of Sec. V-E: friction sweep at fixed objective-phase exposure.

  python -m experiments.table6_friction_sweep            Table VI, seed-11 landscape
  python -m experiments.table6_friction_sweep dense      four further frictions on the seed-11 landscape
  python -m experiments.table6_friction_sweep draws      four further draws (seeds 12 to 15), coarser grid
  python -m experiments.table6_friction_sweep steps      sensitivity to the integration step, all five draws
  python -m experiments.table6_friction_sweep periodic   the periodic spectral Laplacian, three draws

All runs use m = 5 bits per variable. rho is the friction in units of critical damping. The sections
also store their rows under results/, where the script of Fig. 4 reads those of the first three.
"""
import json
import pathlib
import sys

import numpy as np

from core.landscape import XI_STAR, curvature, objective
from core.qhd import ALPHA_OBJ, evolve, n_steps, register_objective, rescale, schedule, summary
from core.report import Published, header

PRINTED = {  # rho: P(F <= 0.01), P(F <= 0.05), E[F]
    0.004: (0.552, 0.947, 0.0152), 0.016: (0.925, 0.998, 0.0038), 0.032: (0.961, 0.998, 0.0023),
    0.125: (0.959, 0.997, 0.0018), 0.25: (0.956, 0.996, 0.0019), 1: (0.832, 0.988, 0.0061),
    4: (0.151, 0.747, 0.0387),
}
COARSE = (0.004, 0.016, 0.0625, 0.25, 1)
DENSE = (0.002, 0.008, 0.0625, 0.5)
UPPER = (0.25, 0.5, 0.75, 1, 1.5, 2, 4)     # frictions of the runs with a fixed number of steps
RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"
REFINED = 0.02 / 3     # integration step of the refinement check
REFINE = (0.016, 0.032, 0.0625, 0.125, 0.25)     # its frictions: the peak, its neighbours, the plateau
COLUMNS = (f"{'rho':>7} {'S':>8} {'lambda(S)':>10} {'steps':>7} {'ds':>7} {'P(F<=0.01)':>11} "
           f"{'P(F<=0.05)':>11} {'E[F]':>8} {'TV':>6}")


def landscape(seed):
    F, arrays = objective(6, 2, seed)
    mu = curvature(F)[0][0]
    f_grid, grid = register_objective(arrays, 5)
    return mu, rescale(f_grid), grid


def run(mu, f_tilde, grid, rho, floor=3000, ds_max=0.02, steps=None, periodic=False):
    nu, lam_end, horizon = schedule(rho, mu)
    steps = steps or n_steps(horizon, floor, ds_max)
    out = summary(evolve(f_tilde, nu, horizon, steps, periodic=periodic), f_tilde, grid, XI_STAR)
    out.update(rho=rho, horizon=horizon, lam_end=lam_end, steps=steps, ds=horizon / steps)
    return out


def line(r):
    return (f"{r['rho']:>7g} {r['horizon']:>8.2f} {r['lam_end']:>10.0f} {r['steps']:>7d} {r['ds']:>7.4f} "
            f"{r['p01']:>11.4f} {r['p05']:>11.4f} {r['mean']:>8.5f} {r['tv']:>6.3f}")


def store(section, content):
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / f"friction_sweep_{section}.json").write_text(json.dumps(content, indent=1) + "\n")


def table():
    header("Table VI. Friction sweep at fixed exposure, seed-11 landscape")
    mu, f_tilde, grid = landscape(11)
    nu, lam_end, horizon = schedule(1, mu)
    print(f"mu = {mu:.4f}, 2 sqrt(mu) = {2 * np.sqrt(mu):.4f}, "
          f"exposure E = {6 * ALPHA_OBJ * (lam_end - 1) / nu:.4e} with alpha_obj = {ALPHA_OBJ}")
    print(COLUMNS)
    pub = Published()
    rows = []
    for rho in PRINTED:
        r = run(mu, f_tilde, grid, rho)
        rows.append(r)
        print(line(r), flush=True)
        for key, printed, dec in zip(("p01", "p05", "mean"), PRINTED[rho], (3, 3, 4)):
            pub.add(f"rho = {rho:g} {key}", r[key], printed, dec)
    store("table", dict(seed=11, mu=mu, rows=rows))
    return pub.report()


def dense():
    header("Sec. V-E. Further frictions on the seed-11 landscape")
    mu, f_tilde, grid = landscape(11)
    print(COLUMNS)
    rows = []
    for rho in DENSE:
        rows.append(run(mu, f_tilde, grid, rho))
        print(line(rows[-1]), flush=True)
    store("dense", dict(seed=11, mu=mu, rows=rows))
    return 0


def draws():
    header("Sec. V-E. Four further draws on the coarser grid")
    stored = {}
    for seed in (12, 13, 14, 15):
        mu, f_tilde, grid = landscape(seed)
        print(f"seed {seed}: mu = {mu:.4f}")
        print(COLUMNS)
        rows = [run(mu, f_tilde, grid, rho) for rho in COARSE]
        stored[seed] = dict(mu=mu, rows=rows)
        for r in rows:
            print(line(r), flush=True)
        best_p = max(rows, key=lambda r: r["p01"])
        best_e = min(rows, key=lambda r: r["mean"])
        crit = rows[-1]
        print(f"  P(F<=0.01) peaks at rho = {best_p['rho']:g}, E[F] is smallest at rho = {best_e['rho']:g}; "
              f"against rho = 1: {100 * (best_p['p01'] - crit['p01']):+.1f} points, "
              f"E[F] smaller by a factor of {crit['mean'] / best_e['mean']:.2f}")
        print()
    store("draws", dict(draws=stored))
    return 0


def steps():
    header("Sec. V-E. Sensitivity to the number of integration steps")
    print("3 x 10^3 against 10^4 steps: five landscapes, seven frictions, three statistics")
    worst, stored = {}, {}
    for seed in (11, 12, 13, 14, 15):
        mu, f_tilde, grid = landscape(seed)
        if seed == 11:
            seed_11 = mu, f_tilde, grid
        print(f"seed {seed}: mu = {mu:.4f}")
        print(COLUMNS + f" {'change':>8}")
        pairs = []
        for rho in UPPER:
            a = run(mu, f_tilde, grid, rho, steps=3000)
            b = run(mu, f_tilde, grid, rho, steps=10000)
            pairs.append((a, b))
            change = max(abs(a[k] - b[k]) for k in ("p01", "p05", "mean"))
            print(f"{line(b)} {change:>8.1e}", flush=True)
            worst[seed] = max(worst.get(seed, 0.0), change)
        stored[seed] = dict(mu=mu, rows_3000=[a for a, _ in pairs], rows_10000=[b for _, b in pairs])
    print("largest change of the three statistics: "
          + ", ".join(f"seed {seed} {change:.2e}" for seed, change in worst.items()))
    print(f"over the five landscapes: {max(worst.values()):.2e}")
    print()
    print("step refined from 0.02 to 0.02 / 3 around the peak and on the plateau, seed-11 landscape")
    print(COLUMNS + f" {'change P':>9} {'E[F]':>8}")
    worst_p = worst_e = 0.0
    refinement = []
    for rho in REFINE:
        a = run(*seed_11, rho)
        b = run(*seed_11, rho, ds_max=REFINED)
        refinement.append(dict(rho=rho, base=a, refined=b))
        dp = max(abs(a["p01"] - b["p01"]), abs(a["p05"] - b["p05"]))
        de = abs(a["mean"] - b["mean"])
        print(f"{line(b)} {dp:>9.1e} {de:>8.1e}", flush=True)
        worst_p, worst_e = max(worst_p, dp), max(worst_e, de)
    print(f"largest change of P: {worst_p:.2e}, of E[F]: {worst_e:.2e}; the refined runs peak at "
          f"rho = {max(refinement, key=lambda r: r['refined']['p01'])['rho']:g}")
    store("steps", dict(draws=stored, refinement=refinement))
    pub = Published()
    pub.add("3e3 against 1e4 steps, largest change, in 1e-4", max(worst.values()) * 1e4, 1.9, 1)
    pub.add("refined step, largest change of P, in 1e-4", worst_p * 1e4, 2.1, 1)
    pub.add("refined step, largest change of E[F], in 1e-6", worst_e * 1e6, 8.1, 1)
    return pub.report()


def periodic():
    header("Sec. VI. Operator dependence, three landscapes")
    print("periodic spectral Laplacian at 10^4 steps, and the interval Laplacian where it peaks")
    stored = {}
    for seed in (11, 12, 13):
        mu, f_tilde, grid = landscape(seed)
        print(f"seed {seed}: mu = {mu:.4f}")
        print(COLUMNS)
        rows = []
        for rho in UPPER:
            rows.append(run(mu, f_tilde, grid, rho, steps=10000, periodic=True))
            print(line(rows[-1]), flush=True)
        best = max(rows, key=lambda r: r["p01"])
        interval = run(mu, f_tilde, grid, best["rho"], steps=10000)
        stored[seed] = dict(mu=mu, rows=rows, interval_at_peak=interval)
        print(f"  P(F<=0.01) peaks at rho = {best['rho']:g} with {best['p01']:.5f}; "
              f"the interval Laplacian gives {interval['p01']:.4f} there")
        print()
    store("periodic", dict(draws=stored))
    return 0


if __name__ == "__main__":
    section = sys.argv[1] if len(sys.argv) > 1 else "table"
    sections = {"table": table, "dense": dense, "draws": draws, "steps": steps, "periodic": periodic}
    sys.exit(1 if sections[section]() else 0)
