"""Benchmark of Sec. IV of the paper.

A constant-coefficient d-asset diffusion in log price on a periodic spectral grid. The covariance is
in factor form, Sigma(xi) = diag(delta) + sum_mu xi_mu v_mu v_mu^T, so the Fourier symbol of the
generator is affine in the design variables, a(k, xi) = a_0(k) + sum_mu xi_mu c_mu(k). The objective
is the full-domain squared residual against an in-family target, evaluated mode by mode.
"""
import numpy as np
from scipy.optimize import minimize

RATE = 0.05                               # short rate r
HORIZON = 1.0                             # T
WIDTH = 6.0                               # width of the periodic grid
XI_STAR = np.array([0.62, 0.28, 0.75])    # design coordinates of the in-family target
N_KEEP = 20000                            # heaviest modes kept on grids larger than this


def factor_model(d, n_sp, seed, orth=False, n_factors=3):
    """Mode weights, a_0(k) and c_mu(k) for one draw of the factor structure.

    The grid has 2**n_sp points per direction. With orth=True the loadings are replaced by their
    QR factor at unchanged column norms.
    """
    rng = np.random.default_rng(seed)
    nx = 2 ** n_sp
    x1 = (np.arange(nx) - nx // 2) * WIDTH / nx
    k1 = 2 * np.pi * np.fft.fftfreq(nx, d=WIDTH / nx)
    kg = np.stack(np.meshgrid(*[k1] * d, indexing="ij")).reshape(d, -1)
    xg = np.stack(np.meshgrid(*[x1] * d, indexing="ij")).reshape(d, -1)
    u0 = (np.exp(-(xg ** 2).sum(0) / (2 * 0.8 ** 2))
          - 0.6 * np.exp(-((xg - 0.9) ** 2).sum(0) / (2 * 0.8 ** 2)))
    u0 /= np.linalg.norm(u0)
    weight = np.abs(np.fft.fftn(u0.reshape([nx] * d)).ravel() / np.sqrt(u0.size)) ** 2
    delta = rng.uniform(0.01, 0.04, d)
    v = rng.normal(size=(d, n_factors)) / np.sqrt(d) * 0.55
    if orth:
        q, _ = np.linalg.qr(v)
        v = q * np.linalg.norm(v, axis=0)
    s0 = (delta[:, None] * kg ** 2).sum(0)
    p_mu = np.stack([(v[:, m] @ kg) ** 2 for m in range(n_factors)])
    t0 = kg.sum(0)
    d0 = (delta[:, None] * kg).sum(0)
    e_mu = np.stack([((v[:, m] ** 2)[:, None] * kg).sum(0) for m in range(n_factors)])
    a0 = RATE + 0.5 * s0 - 1j * (RATE * t0 - 0.5 * d0)
    c_mu = 0.5 * p_mu + 0.5j * e_mu
    return weight, a0, c_mu


def objective(d, n_sp, seed, orth=False, n_keep=N_KEEP, xi_target=XI_STAR):
    """F(xi) for the in-family target, with the arrays it is built from.

    Returns (F, (w, a0, c, target)), restricted to the n_keep heaviest modes when the grid is larger.
    """
    weight, a0, c_mu = factor_model(d, n_sp, seed, orth=orth, n_factors=len(xi_target))
    keep = np.argsort(weight)[::-1][:n_keep] if weight.size > n_keep else np.arange(weight.size)
    w, a0, c = weight[keep], a0[keep], c_mu[:, keep]

    def amplitude(xi):
        v = np.exp(-(a0 + np.asarray(xi) @ c) * HORIZON)
        mag = np.abs(v)
        return np.where(mag > 1.0, v / mag, v)

    target = amplitude(xi_target)

    def F(xi):
        return float((w * np.abs(amplitude(xi) - target) ** 2).sum())

    return F, (w, a0, c, target)


def curvature(F, n_var=3, n_scan=13, step=5e-3):
    """Hessian spectrum of the rescaled objective at its minimiser.

    The objective range is taken from an n_scan-point-per-axis scan of the box, the minimiser from
    Nelder-Mead, and the Hessian from central differences.
    Returns (eigenvalues, minimiser, scan minimum, scan maximum).
    """
    g = np.linspace(0.0, 1.0, n_scan)
    scan = [F(np.array([a, b, c])) for a in g for b in g for c in g]
    f_lo, f_hi = min(scan), max(scan)
    scale = f_hi - f_lo
    opt = minimize(F, np.full(n_var, 0.5), method="Nelder-Mead",
                   options=dict(xatol=1e-8, fatol=1e-17, maxiter=4000))
    z, f0 = opt.x, opt.fun
    hess = np.zeros((n_var, n_var))
    for i in range(n_var):
        ei = np.eye(n_var)[i] * step
        hess[i, i] = (F(z + ei) - 2 * f0 + F(z - ei)) / step ** 2 / scale
        for j in range(i + 1, n_var):
            ej = np.eye(n_var)[j] * step
            hess[i, j] = hess[j, i] = (F(z + ei + ej) - F(z + ei - ej)
                                       - F(z - ei + ej) + F(z - ei - ej)) / (4 * step ** 2) / scale
    return np.linalg.eigvalsh(hess), z, f_lo, f_hi
