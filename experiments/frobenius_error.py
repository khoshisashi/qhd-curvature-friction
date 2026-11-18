"""Sec. V-B: Frobenius error of the factor family on the thirty targets of Table II.

The family is Sigma = diag(delta) + sum_mu xi_mu v_mu v_mu^T, the diagonal as drawn. For a fixed
random frame the loadings are those of Table II and the coefficients are fitted on the box
[0, 1]^M. With fitted directions the loadings are free as well, which leaves the best positive
semidefinite approximation of rank M to the target minus the diagonal. Both errors are Frobenius
distances from the target, relative to its Frobenius norm.

  python -m experiments.frobenius_error
"""
import sys

import numpy as np
from scipy.optimize import lsq_linear

from core.report import Published, header
from experiments.table2_representation_gap import N_TARGETS, frame, targets

FACTORS = (3, 5)
KINDS = {"fitted": "fitted directions", "frame": "fixed random frame"}


def fixed_frame(residual, V):
    """Distance of the residual from sum_mu xi_mu v_mu v_mu^T, minimised over xi in [0, 1]^M."""
    atoms = np.stack([np.outer(v, v).ravel() for v in V.T], axis=1)
    xi = lsq_linear(atoms, residual.ravel(), bounds=(0.0, 1.0), method="bvls").x
    return np.linalg.norm(residual - (V * xi) @ V.T)


def fitted_directions(residual, n_factors):
    """The same distance with free loadings: the largest eigenvalues are kept where positive."""
    w, U = np.linalg.eigh(residual)
    top = U[:, -n_factors:]
    return np.linalg.norm(residual - (top * np.clip(w[-n_factors:], 0.0, None)) @ top.T)


def main():
    header("Sec. V-B. Frobenius error of the factor family on the thirty targets of Table II")
    error = {(M, kind): [] for M in FACTORS for kind in KINDS}
    for t, (S_target, _) in enumerate(targets(N_TARGETS)):
        norm = np.linalg.norm(S_target)
        for M in FACTORS:
            V, delta = frame(M, 11 + t)
            residual = S_target - np.diag(delta)
            error[M, "fitted"].append(fitted_directions(residual, M) / norm)
            error[M, "frame"].append(fixed_frame(residual, V) / norm)

    print(f"{'relative Frobenius error':<28}{'min':>8}{'median':>8}{'max':>8}")
    for (M, kind), values in error.items():
        print(f"{f'M={M}, {KINDS[kind]}':<28}{min(values):>8.3f}{np.median(values):>8.3f}{max(values):>8.3f}")

    pub = Published()
    pub.add("Sec. V-B: median error at M=3, fitted directions", np.median(error[3, "fitted"]), 0.12, 2)
    pub.add("Sec. V-B: median error at M=3, fixed random frame", np.median(error[3, "frame"]), 0.79, 2)
    return 1 if pub.report() else 0


if __name__ == "__main__":
    sys.exit(main())
