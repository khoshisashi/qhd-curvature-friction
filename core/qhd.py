"""Quantum Hamiltonian Descent on the design register, integrated classically (Secs. II and IV).

H(s) = -Lambda / (2 lambda(s)) + lambda(s) F, lambda(s) = exp(nu s), on M = 3 design variables of
m bits each. Lambda is built from h^-2 tridiag(1, -2, 1) without wrap-around, which the type-I
discrete sine transform diagonalises, so the kinetic factor of the Strang splitting is exact.
"""
import numpy as np
from scipy.fft import dstn

from .landscape import HORIZON

EPSILON = 4.5e-5        # terminal-lambda parameter
ALPHA_OBJ = 2.51        # subnormalisation of the objective block encoding, enters the exposure only


def register_objective(arrays, m_bits, chunk=2048):
    """Objective on the endpoint-inclusive register grid. Returns (F of shape (G, G, G), grid)."""
    w, a0, c, target = arrays
    n = 2 ** m_bits
    grid = np.linspace(0.0, 1.0, n)
    pts = np.stack(np.meshgrid(grid, grid, grid, indexing="ij"), -1).reshape(-1, 3)
    out = np.empty(len(pts))
    for i in range(0, len(pts), chunk):
        v = np.exp(-(a0[None, :] + pts[i:i + chunk] @ c) * HORIZON)
        mag = np.abs(v)
        v = np.where(mag > 1.0, v / mag, v)
        out[i:i + chunk] = (w[None, :] * np.abs(v - target[None, :]) ** 2).sum(1)
    return out.reshape(n, n, n), grid


def rescale(f_grid):
    """Objective rescaled to the unit interval over the register points."""
    return (f_grid - f_grid.min()) / (f_grid.max() - f_grid.min())


def interval_spectrum(n):
    """Eigenvalues of h^-2 tridiag(1, -2, 1) on n endpoint-inclusive points, h = 1 / (n - 1)."""
    h = 1.0 / (n - 1)
    j = np.arange(1, n + 1)
    return -(4.0 / h ** 2) * np.sin(j * np.pi / (2 * (n + 1))) ** 2


def schedule(rho, mu, eps=EPSILON):
    """Friction, terminal lambda and horizon at fixed objective-phase exposure.

    The exposure is held at its value under the terminal rule lambda(S) = 1 / eps at critical
    damping, which gives lambda(S) = 1 + rho (1 / eps - 1) whatever alpha_obj.
    """
    nu = rho * 2.0 * np.sqrt(mu)
    lam_end = 1.0 + rho * (1.0 / eps - 1.0)
    return nu, lam_end, np.log(lam_end) / nu


def n_steps(horizon, floor=3000, ds_max=0.02):
    """Integration steps: step at ds_max or below, with a floor on the number of steps."""
    return max(floor, int(np.ceil(horizon / ds_max)))


def evolve(f_tilde, nu, horizon, steps, periodic=False):
    """Final measurement distribution from the uniform superposition.

    Strang splitting with kinetic half steps around the potential step and lambda at the step
    midpoint. Adjacent kinetic half steps are merged in the spectral domain. With periodic=True the
    interval Laplacian is replaced by the periodic spectral one.
    """
    n = f_tilde.shape[0]
    ds = horizon / steps
    if periodic:
        e1 = -(2 * np.pi * np.fft.fftfreq(n, d=1.0 / n)) ** 2
        forward, inverse = np.fft.fftn, np.fft.ifftn
    else:
        e1 = interval_spectrum(n)

        def forward(a):
            return dstn(a, type=1, norm="ortho")

        inverse = forward
    e3 = e1[:, None, None] + e1[None, :, None] + e1[None, None, :]
    lam = np.exp(nu * (np.arange(steps) + 0.5) * ds)
    psi = np.full((n, n, n), 1.0 / np.sqrt(n ** 3), dtype=complex)
    phi = forward(psi) * np.exp(1j * e3 * ds / (4.0 * lam[0]))
    for k in range(steps):
        psi = inverse(phi)
        psi *= np.exp(-1j * lam[k] * f_tilde * ds)
        phi = forward(psi)
        inv_lam = 1.0 / lam[k] + (1.0 / lam[k + 1] if k + 1 < steps else 0.0)
        phi *= np.exp(1j * e3 * ds * inv_lam / 4.0)
    p = np.abs(inverse(phi)) ** 2
    return p / p.sum()


def summary(p, f_tilde, grid, xi_star):
    """P(F <= 0.01), P(F <= 0.05), E[F], distance of the mode from xi_star, total variation to uniform."""
    i = np.unravel_index(int(np.argmax(p)), p.shape)
    mode = np.array([grid[i[0]], grid[i[1]], grid[i[2]]])
    return dict(p01=float(p[f_tilde <= 0.01].sum()), p05=float(p[f_tilde <= 0.05].sum()),
                mean=float((p * f_tilde).sum()), dist=float(np.linalg.norm(mode - xi_star)),
                tv=float(0.5 * np.abs(p - 1.0 / p.size).sum()))
