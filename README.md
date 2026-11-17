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
bash run_all.sh                                    # everything, about an hour and a half on one core
```

The experiments are run as modules from the root of the repository. Those for the tables and for the
statements in the text print each computed value next to the value printed in the paper, and exit
with a non-zero status if one of them differs by more than one unit in the last printed digit. The
sections `dense`, `draws` and `periodic` of the sweep print the quantities behind statements of
Secs. V-E and VI without such a comparison. The output of every experiment is stored under `logs/`.

## Layout

```
core/            shared modules
  landscape.py     benchmark of Sec. IV and curvature measurement
  qhd.py           evolution on the design register
  report.py        run header and comparison with the printed values
experiments/     one script per table or group of statements
logs/            stored output of every experiment
run_all.sh       runs every experiment and rewrites logs/
```

The figure script writes to `figures/`, which is created on the first run and not tracked.

## What produces what

| Paper | `python -m experiments.` | Log in `logs/` | Time |
|---|---|---|---|
| Table I, feasible fraction of the entry parameterisation | `table1_feasible_fraction` | `table1_feasible_fraction.txt` | 2 s |
| Sec. V-A, infeasible points of the `d = 2` entry register | `entry_register_d2` | `entry_register_d2.txt` | 1 s |
| Tables III and IV, thirty factor structures, as drawn and orthogonalised | `table3_table4_conditioning` | `table3_table4_conditioning.txt` | 20 s |
| Table V, design bit depth | `table5_bit_depth` | `table5_bit_depth.txt` | 10 min |
| Table VI, friction sweep on the seed-11 landscape | `table6_friction_sweep` | `table6_friction_sweep.txt` | 8 min |
| Sec. V-E, the frictions 0.002 and 0.0625 | `table6_friction_sweep dense` | `table6_friction_sweep_dense.txt` | 8 min |
| Sec. V-E, four further draws | `table6_friction_sweep draws` | `table6_friction_sweep_draws.txt` | 40 min |
| Sec. V-E, sensitivity to the integration step, all five draws | `table6_friction_sweep steps` | `table6_friction_sweep_steps.txt` | 20 min |
| Sec. V-E, gauge invariance (Proposition 2) | `gauge_invariance` | `gauge_invariance.txt` | 1 min |
| Sec. VI, periodic Laplacian | `table6_friction_sweep periodic` | `table6_friction_sweep_periodic.txt` | 1 min |
| Secs. IV and VI, refinement, truncation and dimension | `refinement_and_dimension` | `refinement_and_dimension.txt` | 4 min |
| Fig. 1, feasible slice | `figure1_feasible_slice` | `figure1_feasible_slice.txt` | 2 s |

The second column is what follows `python -m experiments.` on the command line, for example
`python -m experiments.table6_friction_sweep dense`. The times are rounded wall-clock times on one
core of a shared x86-64 machine, with Python 3.13.

## Not yet in this repository

| Paper | Missing |
|---|---|
| Table II and Fig. 3, the representation gap over thirty out-of-family targets, with the paired counts and ratios quoted in Secs. V-B and V-C | script and run log |
| Fig. 4 | plotting script; the numbers are those of `experiments/table6_friction_sweep.py` |
| Fig. 2 | a schematic, without data |

## Conventions the numbers depend on

- **Random draws.** All draws use `numpy.random.default_rng(seed)`. The thirty factor structures of
  Tables III and IV are the seeds 11 to 40, the diagonal drawn before the loadings. The seed-11
  landscape is the first of them, and the four further draws of the sweep are the seeds 12 to 15.
  Table I uses seed 0, one stream across the dimensions.
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

The repository was reorganised for the public release in November 2026. `core/landscape.py` and the
experiments `table1_feasible_fraction`, `entry_register_d2` and `figure1_feasible_slice` carry over
the numerics of the working scripts of the project unchanged. `core/qhd.py` and the experiments
built on it were rewritten for the release from the protocol stated in the paper. The two largest
sweeps of the paper were first run on a GPU; the versions here run on a CPU, and the logs show that
they reproduce the printed entries of Tables V and VI. The reorganisation and the rewriting were
done with the assistance of generative AI (Claude, Anthropic).

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
