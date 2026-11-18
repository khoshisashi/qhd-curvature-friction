# qhd-curvature-friction

Code and run logs for

> K. Hoshisashi, "Curvature-Calibrated Friction and an Exact Gauge Symmetry for Quantum Hamiltonian
> Descent," in *Proc. IEEE Symposium Series on Computational Intelligence (SSCI)*, Gold Coast,
> Australia, Feb. 2027.

Everything runs classically. The optimiser is integrated as a Schrödinger equation on the design
register alone, whose dimension `2^(mM)` does not grow with the physical state. No quantum hardware
and no circuit simulation are involved.

## Quick start

```
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
python -m experiments.table1_feasible_fraction     # a few seconds
bash run_all.sh                                    # everything, about two hours and a quarter on one core
```

The experiments are run as modules from the root of the repository. Those for the tables and for the
statements in the text print each computed value next to the value printed in the paper, and exit
with a non-zero status if one of them differs by more than one unit in the last printed digit. The
sections `dense`, `draws` and `periodic` of the sweep print the quantities behind statements of
Secs. V-E and VI without such a comparison. The output of every experiment is stored under `logs/`.
The experiment of Table II and the sections `table`, `dense` and `draws` of the sweep also store
their rows under `results/`, from which Figs. 3 and 4 are drawn.

## Layout

```
core/            shared modules
  landscape.py     benchmark of Sec. IV and curvature measurement
  qhd.py           evolution on the design register
  report.py        run header and comparison with the printed values
  figures.py       style shared by the figure scripts
experiments/     one script per table, figure or group of statements
logs/            stored output of every experiment
results/         rows from which Figs. 3 and 4 are drawn
run_all.sh       runs every experiment and rewrites logs/ and results/
```

The figure scripts write to `figures/`, which is created on the first run and not tracked.

## What produces what

| Paper | `python -m experiments.` | Log in `logs/` | Time |
|---|---|---|---|
| Table I, feasible fraction of the entry parameterisation | `table1_feasible_fraction` | `table1_feasible_fraction.txt` | 2 s |
| Sec. V-A, infeasible points of the `d = 2` entry register | `entry_register_d2` | `entry_register_d2.txt` | 1 s |
| Table II and the paired comparisons of Secs. V-B and V-C, thirty targets outside the family | `table2_representation_gap` | `table2_representation_gap.txt` | 35 min |
| Tables III and IV, thirty factor structures, as drawn and orthogonalised | `table3_table4_conditioning` | `table3_table4_conditioning.txt` | 20 s |
| Table V, design bit depth | `table5_bit_depth` | `table5_bit_depth.txt` | 10 min |
| Table VI, friction sweep on the seed-11 landscape | `table6_friction_sweep` | `table6_friction_sweep.txt` | 8 min |
| Sec. V-E and Fig. 4, the frictions 0.002, 0.008, 0.0625 and 0.5 | `table6_friction_sweep dense` | `table6_friction_sweep_dense.txt` | 12 min |
| Sec. V-E, four further draws | `table6_friction_sweep draws` | `table6_friction_sweep_draws.txt` | 40 min |
| Sec. V-E, sensitivity to the integration step, all five draws | `table6_friction_sweep steps` | `table6_friction_sweep_steps.txt` | 20 min |
| Sec. V-E, gauge invariance (Proposition 2) | `gauge_invariance` | `gauge_invariance.txt` | 1 min |
| Sec. VI, periodic Laplacian | `table6_friction_sweep periodic` | `table6_friction_sweep_periodic.txt` | 1 min |
| Secs. IV and VI, refinement, truncation and dimension | `refinement_and_dimension` | `refinement_and_dimension.txt` | 4 min |
| Fig. 1, feasible slice | `figure1_feasible_slice` | `figure1_feasible_slice.txt` | 2 s |
| Fig. 3, raw residual on the thirty targets | `figure3_representation_gap` | `figure3_representation_gap.txt` | 1 s |
| Fig. 4, friction sweep on the five draws | `figure4_friction_sweep` | `figure4_friction_sweep.txt` | 1 s |

The second column is what follows `python -m experiments.` on the command line, for example
`python -m experiments.table6_friction_sweep dense`. The times are rounded wall-clock times on one
core of a shared x86-64 machine, with Python 3.13. Fig. 3 needs the result file of Table II, and
Fig. 4 those of the sections `table`, `dense` and `draws` of the sweep; all four are stored.

## Not in this repository

| Paper | Missing |
|---|---|
| Sec. V-B, the two Frobenius errors, `9.2e-2` with fitted directions and `0.93` for a fixed frame | script and run log |
| Fig. 2 | a schematic, without data |

## Conventions the numbers depend on

- **Random draws.** All draws use `numpy.random.default_rng(seed)`. The thirty factor structures of
  Tables III and IV are the seeds 11 to 40, the diagonal drawn before the loadings. The seed-11
  landscape is the first of them, and the four further draws of the sweep are the seeds 12 to 15.
  Table I uses seed 0, one stream across the dimensions.
- **Representation gap.** The thirty targets come from one stream, seed 2, the matrix `B` of a
  target before its diagonal. Target `t = 1, ..., 30` is paired with the factor structure of seed
  `10 + t`, so the three-factor frames are those of Tables III and IV; the scan of a box uses seed
  `99 + t` and the start of the cone minimisation seed `t - 1`. A box is scanned at its corners and
  at 4000 random points, and the minimum is refined by Nelder-Mead from the three best of them. The
  range of a fixed-diagonal family is that of the scan. The box of the variable diagonal has nine
  variables and a minimiser on its boundary: there the minimum is refined from the twelve best
  points by bounded quasi-Newton descent followed by Nelder-Mead, and the range runs from the
  largest value of the scan to that minimum.
- **Curvature.** The objective range is the difference between the largest and the smallest value
  of a 13-point-per-axis scan of the box. With the smallest value replaced by zero, its exact value
  for an in-family target, the curvature changes by under 0.1%.
- **Truncation.** A landscape of the refinement checks with more than 20000 modes keeps the 20000
  heaviest: `d = 4` at `N_x = 16` and `32`, where they carry over 99.99% of the energy, and `d = 6`
  at `N_x = 8`, where they carry 99.88%. Everything else uses all modes.
- **Register objective.** At each bit depth the objective is rescaled to the unit interval over
  the register points, and the tolerance sets `F <= 0.01` and `F <= 0.05` refer to that scale.
- **Schedule.** With `rho` the friction in units of critical damping, `nu = 2 rho sqrt(mu)`,
  `lambda(S) = 1 + rho (1/epsilon - 1)` and `S = log(lambda(S)) / nu`, which holds the exposure at
  its value under the terminal rule at critical damping.
- **Integration.** Strang splitting with kinetic half steps around the potential step and `lambda`
  taken at the step midpoint. The kinetic factor is exact, through the type-I discrete sine
  transform that diagonalises the interval Laplacian. The step is 0.02 or below, with a floor of
  3000 steps, and the refinement check of Sec. V-E divides the largest step by three.
- **`d = 2` illustration.** `experiments/entry_register_d2.py` keeps the set-up of the working
  script it comes from, a grid of width 4 with 32 points per direction, a tent-shaped claim and an
  objective restricted to a central window, which differs from the benchmark of Sec. IV.

## Provenance

The repository was reorganised for the public release in November 2026. `core/landscape.py`, the
experiments `table1_feasible_fraction`, `entry_register_d2` and `table2_representation_gap` and the
three figure scripts carry over the numerics of the working scripts and notebooks of the project.
`core/qhd.py` and the experiments built on it were rewritten for the release from the protocol
stated in the paper. The representation gap and the two largest sweeps of the paper were first run
on a GPU; the versions here run on a CPU, and the logs show that they reproduce the printed entries
of Tables II, V and VI. The variable-diagonal family, the fourth row of Table II, does not come
from a working script: it is implemented here from its description in the paper, with the search
stated under the conventions above, and it reproduces the printed row and the ratios quoted with
it. The reorganisation and the rewriting were done with the assistance of generative AI (Claude,
Anthropic).

## Licence and citation

The code is released under the MIT licence. To cite it, please cite the paper.

```bibtex
@inproceedings{hoshisashi2027curvature,
  author    = {Hoshisashi, Kentaro},
  title     = {Curvature-Calibrated Friction and an Exact Gauge Symmetry for Quantum {H}amiltonian Descent},
  booktitle = {Proc. IEEE Symposium Series on Computational Intelligence (SSCI)},
  address   = {Gold Coast, Australia},
  year      = {2027},
  note      = {To appear}
}
```
