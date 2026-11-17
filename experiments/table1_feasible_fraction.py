"""Table I: fraction of the entry-parameterised design hypercube that decodes to a positive
semidefinite covariance. Diagonal entries are uniform in [0.02, 0.20] and off-diagonal entries in
[-0.06, 0.10], with 2e5 samples per dimension."""
import sys

import numpy as np

from core.report import Published, header

PRINTED = {2: (0.896, 3), 3: (0.635, 3), 4: (0.335, 3), 5: (0.127, 3), 6: (0.033, 3), 8: (0.0008, 4)}
N_SAMPLES = 200000

header("Table I. Feasible fraction of the entry parameterisation")
rng = np.random.default_rng(0)
pub = Published()
print(f"{'d':>3} {'M':>4} {'feasible':>10} {'hits':>8} {'std. error':>11}")
for d in (2, 3, 4, 5, 6, 8):
    hits = 0
    for _ in range(N_SAMPLES // 1000):
        diag = rng.uniform(0.02, 0.20, size=(1000, d))
        off = rng.uniform(-0.06, 0.10, size=(1000, d * (d - 1) // 2))
        sigma = np.zeros((1000, d, d))
        iu = np.triu_indices(d, 1)
        sigma[:, iu[0], iu[1]] = off
        sigma = sigma + np.transpose(sigma, (0, 2, 1))
        sigma[:, np.arange(d), np.arange(d)] = diag
        hits += int((np.linalg.eigvalsh(sigma)[:, 0] >= 0).sum())
    frac = hits / N_SAMPLES
    err = np.sqrt(frac * (1 - frac) / N_SAMPLES)
    print(f"{d:>3} {d * (d + 1) // 2:>4} {frac:>10.5f} {hits:>8} {err:>11.1e}")
    pub.add(f"d = {d}", frac, *PRINTED[d])
sys.exit(1 if pub.report() else 0)
