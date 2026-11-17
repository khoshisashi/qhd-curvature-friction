"""Tables III and IV: conditioning of the factor landscape across thirty draws of the loadings at
d = 6, M = 3 (seeds 11 to 40), as drawn and after a QR orthogonalisation at unchanged column norms."""
import sys

import numpy as np

from core.landscape import curvature, objective
from core.report import Published, header

SEEDS = range(11, 41)
PRINTED_III = {  # min, 25%, median, 75%, max
    "mu": ([0.029, 0.084, 0.161, 0.237, 0.415], 3),
    "kappa": ([3.75, 7.07, 10.12, 20.21, 79.18], 2),
    "nu": ([0.343, 0.578, 0.803, 0.974, 1.289], 3),
}
PRINTED_IV = {  # minimum, median, maximum
    "sqrt(mu) ratio": ([1.00, 1.21, 2.47], 2),
    "kappa as drawn": ([3.75, 10.12, 79.18], 2),
    "kappa orthogonal": ([2.87, 8.54, 63.54], 2),
}


def five(x):
    return [x.min(), np.percentile(x, 25), np.median(x), np.percentile(x, 75), x.max()]


def three(x):
    return [x.min(), np.median(x), x.max()]


header("Tables III and IV. Thirty factor structures at d = 6, M = 3")
rows = []
print(f"{'seed':>4} {'mu':>8} {'kappa':>7} {'nu':>7} | {'mu orth':>8} {'kappa':>7} {'sqrt ratio':>10}")
for seed in SEEDS:
    ev = curvature(objective(6, 2, seed)[0])[0]
    ev_o = curvature(objective(6, 2, seed, orth=True)[0])[0]
    rows.append((ev[0], ev[-1], ev_o[0], ev_o[-1]))
    print(f"{seed:>4} {ev[0]:>8.4f} {ev[-1] / ev[0]:>7.2f} {2 * np.sqrt(ev[0]):>7.4f} | "
          f"{ev_o[0]:>8.4f} {ev_o[-1] / ev_o[0]:>7.2f} {np.sqrt(ev_o[0] / ev[0]):>10.4f}")
mu, big, mu_o, big_o = np.array(rows).T
kappa, kappa_o, nu = big / mu, big_o / mu_o, 2 * np.sqrt(mu)
ratio = np.sqrt(mu_o / mu)

print()
print("Table III            min      25%   median      75%      max")
for name, x, dec in (("mu", mu, 4), ("kappa", kappa, 2), ("nu = 2 sqrt(mu)", nu, 4)):
    print(f"{name:<16}" + "".join(f"{v:>9.{dec}f}" for v in five(x)))
print(f"kappa spans a factor of {kappa.max() / kappa.min():.2f}, nu a factor of {nu.max() / nu.min():.3f}")
print()
print("Table IV             min   median      max")
for name, x in (("sqrt(mu) ratio", ratio), ("kappa as drawn", kappa), ("kappa orthogonal", kappa_o)):
    print(f"{name:<16}" + "".join(f"{v:>9.4f}" for v in three(x)))
worst = int(np.argmin(mu_o / mu))
print(f"orthogonalisation raises mu in {int((mu_o > mu).sum())} of {len(mu)} pairs; "
      f"smallest mu ratio {mu_o[worst] / mu[worst]:.5f} at seed {SEEDS[worst]}")

pub = Published()
for key, x in (("mu", mu), ("kappa", kappa), ("nu", nu)):
    printed, dec = PRINTED_III[key]
    for tag, c, p in zip(("min", "25%", "median", "75%", "max"), five(x), printed):
        pub.add(f"Table III {key} {tag}", c, p, dec)
for key, x in (("sqrt(mu) ratio", ratio), ("kappa as drawn", kappa), ("kappa orthogonal", kappa_o)):
    printed, dec = PRINTED_IV[key]
    for tag, c, p in zip(("min", "median", "max"), three(x), printed):
        pub.add(f"Table IV {key} {tag}", c, p, dec)
pub.add("text: kappa spans a factor of", kappa.max() / kappa.min(), 21, 0)
pub.add("text: nu spans a factor of", nu.max() / nu.min(), 3.75, 2)
sys.exit(1 if pub.report() else 0)
