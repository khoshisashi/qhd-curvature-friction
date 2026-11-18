"""Stored rows of this repository against the result files of the original runs.

original/ holds the four result files of the runs behind the paper. This script compares their
numbers with the rows that the experiments store under results/ and prints, for each group of
values, the largest difference, absolute or relative to the original value. It exits with a
non-zero status if a group goes beyond its tolerance.

Not compared, because the experiments here do not store them: run times, condition numbers and the
exposure integral.
"""
import json
import pathlib
import sys

import numpy as np

from core.report import header
from experiments.table6_friction_sweep import UPPER

ROOT = pathlib.Path(__file__).resolve().parents[1]
STATISTICS = (("P001", "p01"), ("P005", "p05"), ("EF", "mean"))     # names in the original files and here


def load(folder, name):
    return json.loads((ROOT / folder / name).read_text())


class Comparison:
    """Largest difference of each group of values, against a tolerance."""

    def __init__(self):
        self.rows = []

    def add(self, label, ours, theirs, tolerance, relative=False):
        ours = np.asarray(ours, dtype=float).ravel()
        theirs = np.asarray(theirs, dtype=float).ravel()
        if ours.shape != theirs.shape:
            raise ValueError(f"{label}: {ours.size} values here, {theirs.size} in the original file")
        diff = np.abs(ours - theirs) / (np.abs(theirs) if relative else 1.0)
        kind = "relative" if relative else "absolute"
        self.rows.append((label, ours.size, float(diff.max()), tolerance, kind))

    def report(self):
        width = max(len(r[0]) for r in self.rows)
        print(f"{'':<{width}}  {'values':>6}  {'largest difference':>18}  {'tolerance':>9}")
        beyond = 0
        for label, n, worst, tolerance, kind in self.rows:
            flag = "" if worst <= tolerance else "  BEYOND"
            beyond += bool(flag)
            print(f"{label:<{width}}  {n:>6}  {worst:>9.1e} {kind}  {tolerance:>9.0e}{flag}")
        print()
        print(f"{sum(r[1] for r in self.rows)} values in {len(self.rows)} groups, "
              f"{beyond} groups beyond their tolerance")
        return beyond


def representation_gap(cmp):
    """Rows 1 to 3 of Table II, from the run on a GPU."""
    theirs = load("original", "qhd_gap_results.json")
    ours = load("results", "representation_gap.json")
    names = ("M3_drawn", "M3_orth", "M5_drawn")
    for label, here, there, tolerance in (("normalised gap", "gap", "gaps", 1e-6),
                                          ("raw residual", "raw_gap", "raw_gap", 1e-6),
                                          ("range", "range", "raw_range", 1e-9)):
        cmp.add(f"Table II, rows 1 to 3: {label}", [ours["families"][n][here] for n in names],
                [theirs[there][n] for n in names], tolerance, relative=True)
    cmp.add("Table II: cone minimum", ours["cone_residual"], theirs["cone_residual"], 1e-10)


def variable_diagonal(cmp):
    """Row 4 of Table II. The original script does not subtract the cone minimum, which is added
    back to the stored residual here."""
    theirs = load("original", "vardiag_gap.json")
    ours = load("results", "representation_gap.json")
    row = ours["families"]["M3_vardiag"]
    best = np.array(row["raw_gap"]) + np.array(ours["cone_residual"])
    cmp.add("Table II, row 4: box minimum", best, theirs["vardiag_raw"], 1e-6, relative=True)
    cmp.add("Table II, row 4: minimum over range", best / np.array(row["range"]), theirs["vardiag_norm"],
            1e-6, relative=True)


def friction_bracket(cmp):
    """Table VI, Fig. 4 and the refinement of the step, from the bracket of the friction on a GPU."""
    theirs = load("original", "qhd_bracket_results.json")
    table = load("results", "friction_sweep_table.json")
    draws = load("results", "friction_sweep_draws.json")["draws"]
    seed_11 = {r["rho"]: r for r in table["rows"] + load("results", "friction_sweep_dense.json")["rows"]}

    t = theirs["seed11"]
    rows = [seed_11[rho] for rho in t["r"]]
    cmp.add("seed 11, ten frictions: statistics",
            [[r[k] for r in rows] for k in ("p01", "p05", "mean", "tv")],
            [t[k] for k in ("P001", "P005", "EF", "TV_to_uniform")], 1e-9)
    cmp.add("seed 11, ten frictions: S and lambda(S)", [[r[k] for r in rows] for k in ("horizon", "lam_end")],
            [t["S"], t["lam_end"]], 1e-9, relative=True)
    cmp.add("seed 11, ten frictions: steps", [r["steps"] for r in rows], t["steps"], 0)

    ours, there = dict(stat=[], sched=[], steps=[], mu=[]), dict(stat=[], sched=[], steps=[], mu=[])
    for seed, d in sorted(theirs["coarse_all"].items()):
        rows = [seed_11[rho] for rho in d["r"]] if seed == "11" else draws[seed]["rows"]
        if [r["rho"] for r in rows] != d["r"]:
            raise ValueError(f"seed {seed}: the frictions of the coarser grid differ")
        ours["stat"] += [[r[k] for r in rows] for k in ("p01", "mean")]
        there["stat"] += [d["P001"], d["EF"]]
        ours["sched"] += [[r[k] for r in rows] for k in ("horizon", "lam_end")]
        there["sched"] += [d["S"], d["lam_end"]]
        ours["steps"] += [r["steps"] for r in rows]
        there["steps"] += d["steps"]
        ours["mu"].append(table["mu"] if seed == "11" else draws[seed]["mu"])
        there["mu"].append(d["mu"])
    cmp.add("five draws, coarser grid: statistics", ours["stat"], there["stat"], 1e-9)
    cmp.add("five draws, coarser grid: S and lambda(S)", ours["sched"], there["sched"], 1e-9, relative=True)
    cmp.add("five draws, coarser grid: steps", ours["steps"], there["steps"], 0)
    cmp.add("five draws: curvature", ours["mu"], there["mu"], 1e-9, relative=True)

    refinement = {r["rho"]: r for r in load("results", "friction_sweep_steps.json")["refinement"]}
    ours, there = dict(stat=[], ds=[], steps=[]), dict(stat=[], ds=[], steps=[])
    for rho, pair in theirs["step_refinement"].items():
        for mine, original in zip((refinement[float(rho)]["base"], refinement[float(rho)]["refined"]), pair):
            ours["stat"] += [mine[k] for _, k in STATISTICS]
            there["stat"] += [original[k] for k, _ in STATISTICS]
            ours["ds"].append(mine["ds"])
            there["ds"].append(original["dt"])
            ours["steps"].append(mine["steps"])
            there["steps"].append(original["steps"])
    cmp.add("refined step, three frictions: statistics", ours["stat"], there["stat"], 1e-9)
    cmp.add("refined step, three frictions: step", ours["ds"], there["ds"], 1e-9, relative=True)
    cmp.add("refined step, three frictions: steps", ours["steps"], there["steps"], 0)


def first_sweep(cmp):
    """The comparison of step counts, the periodic Laplacian and Table V, from the first sweep on a GPU."""
    theirs = load("original", "qhd_dirichlet_results.json")
    steps = load("results", "friction_sweep_steps.json")["draws"]
    periodic = load("results", "friction_sweep_periodic.json")["draws"]
    depths = load("results", "bit_depth.json")["depths"]
    if list(theirs["friction"]) != list(UPPER):
        raise ValueError("the frictions of the original sweep differ from those of the section steps")

    sweep = theirs["sweep_dirichlet"]
    seeds = sorted(sweep)
    for n in (3000, 10000):
        cmp.add(f"five draws, seven frictions, {n} steps: statistics",
                [[r[k] for r in steps[s][f"rows_{n}"]] for s in seeds for _, k in STATISTICS],
                [sweep[s][f"{k}_{n}"] for s in seeds for k, _ in STATISTICS], 1e-9)
    cmp.add("five draws, seven frictions: S and lambda(S)",
            [[r[k] for r in steps[s]["rows_10000"]] for s in seeds for k in ("horizon", "lam_end")],
            [sweep[s][k] for s in seeds for k in ("S", "lam_end")], 1e-9, relative=True)
    cmp.add("five draws, seven frictions: curvature", [steps[s]["mu"] for s in seeds],
            [sweep[s]["mu"] for s in seeds], 1e-9, relative=True)

    seeds = sorted(theirs["sweep_periodic"])
    cmp.add("periodic Laplacian, three draws: P(F<=0.01)",
            [[r["p01"] for r in periodic[s]["rows"]] for s in seeds],
            [theirs["sweep_periodic"][s] for s in seeds], 1e-9)

    keys = [(m, n) for m in (4, 5, 6) for n in (3000, 10000)]
    cmp.add("Table V, both step counts: P(F<=0.05)", [depths[str(m)][f"steps_{n}"]["p05"] for m, n in keys],
            [theirs["bit_depth"][f"{m}_{n}"]["P005"] for m, n in keys], 1e-9)
    # the original measures the distance from the minimiser of the fit, this repository from xi*
    cmp.add("Table V, both step counts: distance", [depths[str(m)][f"steps_{n}"]["dist"] for m, n in keys],
            [theirs["bit_depth"][f"{m}_{n}"]["err"] for m, n in keys], 1e-6)


def main():
    header("Stored rows against the result files of the original runs")
    cmp = Comparison()
    for part in (representation_gap, variable_diagonal, friction_bracket, first_sweep):
        part(cmp)
    return 1 if cmp.report() else 0


if __name__ == "__main__":
    sys.exit(main())
